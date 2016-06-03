
#include <fstream>

#include <limits>

#include <boost/lexical_cast.hpp>

#include "scsi/constants.h"
#include "scsi/moment2.h"

#include "scsi/moment2.h"
#include "scsi/moment2_sup.h"
#include "scsi/rf_cavity.h"
#include "scsi/chg_stripper.h"

#include "scsi/h5loader.h"


namespace {

} // namespace

Moment2State::Moment2State(const Config& c)
    :StateBase(c)
    ,ref()
    ,real()
    ,moment0(maxsize, 0e0)
    ,state(boost::numeric::ublas::identity_matrix<double>(maxsize))
{
    /* check to see if "cstate" is defined.
     * If "cstate" is defined then we are simulating one of many charge states.
     * If not take IonZ, moment0, and state directly from config variables.
     * If so, then
     *   for IonZ expect the config vector "IonChargeStates" and index with "cstate" (0 index).
     *   append the value of "cstate" to the vector and matrix variable names.
     *     eg. cstate=1, vector_variable=S -> looks for variable "S1".
     */
    double icstate_f = 0.0;
    bool multistate = c.tryGet<double>("cstate", icstate_f);
    size_t icstate = (size_t)icstate_f;

    std::string vectorname(c.get<std::string>("vector_variable", "moment0"));
    std::string matrixname(c.get<std::string>("matrix_variable", "initial"));
    std::vector<double> ics;

    if(multistate) {
        ics = c.get<std::vector<double> >("IonChargeStates");
        if(ics.empty())
            throw std::invalid_argument("IonChargeStates w/ length 0");
        if(icstate>=ics.size())
            throw std::invalid_argument("IonChargeStates[cstate] is out of bounds");
        ref.IonZ = ics[icstate];

        std::string icstate_s(boost::lexical_cast<std::string>(icstate));
        vectorname  += icstate_s;
        matrixname  += icstate_s;
    }

    try{
        const std::vector<double>& I = c.get<std::vector<double> >(vectorname);
        if(I.size()!=moment0.size())
            throw std::invalid_argument("Initial moment0 size mis-match");
        std::copy(I.begin(), I.end(), moment0.begin());
    }catch(key_error&){
        if(multistate)
            throw std::invalid_argument(vectorname+" not defined");
        // default to zeros
    }catch(boost::bad_any_cast&){
        throw std::invalid_argument("'initial' has wrong type (must be vector)");
    }

    try{
        const std::vector<double>& I = c.get<std::vector<double> >(matrixname);
        if(I.size()!=state.size1()*state.size2())
            throw std::invalid_argument("Initial state size mis-match");
        std::copy(I.begin(), I.end(), state.data().begin());
    }catch(key_error&){
        if(multistate)
            throw std::invalid_argument(matrixname+" not defined");
        // default to identity
    }catch(boost::bad_any_cast&){
        throw std::invalid_argument("'initial' has wrong type (must be vector)");
    }

    ref.IonEs      = c.get<double>("IonEs", 0e0),
    ref.IonEk      = c.get<double>("IonEk", 0e0);

    ref.SampleIonK = (ref.IonEs != 0e0)? 2e0*M_PI/(ref.beta*SampleLambda) : 2e0*M_PI/SampleLambda;

    real           = ref;

    real.phis      = moment0[PS_S];
    real.IonEk    += moment0[PS_PS]*MeVtoeV;

    ref.recalc();
    real.recalc();

    if(!multistate) {
        real.IonZ = ref.IonZ       = c.get<double>("IonZ", 0e0);
    } else {
        ref.IonZ  = ics[0];
        real.IonZ = ics[icstate];
    }
}

Moment2State::~Moment2State() {}

Moment2State::Moment2State(const Moment2State& o, clone_tag t)
    :StateBase(o, t)
    ,ref(o.ref)
    ,real(o.real)
    ,moment0(o.moment0)
    ,state(o.state)
{}

void Moment2State::assign(const StateBase& other)
{
    const Moment2State *O = dynamic_cast<const Moment2State*>(&other);
    if(!O)
        throw std::invalid_argument("Can't assign State: incompatible types");
    ref     = O->ref;
    real    = O->real;
    moment0 = O->moment0;
    state   = O->state;
    StateBase::assign(other);
}

void Moment2State::show(std::ostream& strm) const
{
    int j, k;

    strm << std::scientific << std::setprecision(8)
         << "\nState:\n  energy [eV] =\n" << std::setw(20) << real.IonEk << "\n  moment0 =\n    ";
    for (k = 0; k < Moment2State::maxsize; k++)
        strm << std::scientific << std::setprecision(8) << std::setw(16) << moment0(k);
    strm << "\n  state =\n";
    for (j = 0; j < Moment2State::maxsize; j++) {
        strm << "    ";
        for (k = 0; k < Moment2State::maxsize; k++) {
            strm << std::scientific << std::setprecision(8) << std::setw(16) << state(j, k);
        }
        if (j < Moment2State::maxsize-1) strm << "\n";
    }
}

bool Moment2State::getArray(unsigned idx, ArrayInfo& Info) {
    unsigned I=0;
    if(idx==I++) {
        Info.name = "state";
        Info.ptr = &state(0,0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 2;
        Info.dim[0] = state.size1();
        Info.dim[1] = state.size2();
        return true;
    } else if(idx==I++) {
        Info.name = "moment0";
        Info.ptr = &moment0(0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim[0] = moment0.size();
        return true;
    } else if(idx==I++) {
        Info.name = "ref_IonZ";
        Info.ptr = &ref.IonZ;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_IonEs";
        Info.ptr = &ref.IonEs;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_IonW";
        Info.ptr = &ref.IonW;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_gamma";
        Info.ptr = &ref.gamma;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_beta";
        Info.ptr = &ref.beta;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_bg";
        Info.ptr = &ref.bg;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_SampleIonK";
        Info.ptr = &ref.SampleIonK;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_phis";
        Info.ptr = &ref.phis;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_IonEk";
        Info.ptr = &ref.IonEk;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_IonZ";
        Info.ptr = &real.IonZ;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_IonEs";
        Info.ptr = &real.IonEs;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_IonW";
        Info.ptr = &real.IonW;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_gamma";
        Info.ptr = &real.gamma;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_beta";
        Info.ptr = &real.beta;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_bg";
        Info.ptr = &real.bg;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_SampleIonK";
        Info.ptr = &real.SampleIonK;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_phis";
        Info.ptr = &real.phis;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_IonEk";
        Info.ptr = &real.IonEk;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    }
    return StateBase::getArray(idx-I, Info);
}

Moment2ElementBase::Moment2ElementBase(const Config& c)
    :ElementVoid(c)
    ,transfer(state_t::maxsize, state_t::maxsize)
    ,transfer_raw(boost::numeric::ublas::identity_matrix<double>(state_t::maxsize))
    ,misalign(boost::numeric::ublas::identity_matrix<double>(state_t::maxsize))
    ,misalign_inv(state_t::maxsize, state_t::maxsize)
    ,scratch(state_t::maxsize, state_t::maxsize)
{

    inverse(misalign_inv, misalign);

    // spoil to force recalculation of energy dependent terms
    last_Kenergy_in = last_Kenergy_out = std::numeric_limits<double>::quiet_NaN();
}

Moment2ElementBase::~Moment2ElementBase() {}

void Moment2ElementBase::assign(const ElementVoid *other)
{
    const Moment2ElementBase *O = static_cast<const Moment2ElementBase*>(other);
    length = O->length;
    transfer = O->transfer;
    transfer_raw = O->transfer_raw;
    misalign = O->misalign;
    ElementVoid::assign(other);

    // spoil to force recalculation of energy dependent terms
    last_Kenergy_in = last_Kenergy_out = std::numeric_limits<double>::quiet_NaN();
}

void Moment2ElementBase::show(std::ostream& strm) const
{
    using namespace boost::numeric::ublas;
    ElementVoid::show(strm);
    strm<<"Length "<<length<<"\n"
          "Transfer: "<<transfer<<"\n"
          "Transfer Raw: "<<transfer_raw<<"\n"
          "Mis-align: "<<misalign<<"\n";
}

void Moment2ElementBase::get_misalign(state_t& ST)
{
    value_mat R,
              R_inv,
              scl     = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize),
              scl_inv = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize),
              T       = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize),
              T_inv   = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

    double dx    = conf().get<double>("dx",    0e0)*MtoMM,
           dy    = conf().get<double>("dy",    0e0)*MtoMM,
           pitch = conf().get<double>("pitch", 0e0),
           yaw   = conf().get<double>("yaw",   0e0),
           tilt  = conf().get<double>("tilt",  0e0);

    scl(state_t::PS_S, state_t::PS_S)   /= -ST.real.SampleIonK;
    scl(state_t::PS_PS, state_t::PS_PS) /= sqr(ST.real.beta)*ST.real.gamma*ST.ref.IonEs/MeVtoeV;

    inverse(scl_inv, scl);

    // Translate to center of element.
    T(state_t::PS_S,  6) = -length/2e0*MtoMM;
    T(state_t::PS_PS, 6) = 1e0;
    inverse(T_inv, T);

    RotMat(dx, dy, pitch, yaw, tilt, R);

    misalign = prod(T, scl);
    misalign = prod(R, misalign);
    misalign = prod(T_inv, misalign);
    misalign = prod(scl_inv, misalign);

    // Can not use inverse or transpose of R.
    RotMat(-dx, -dy, -pitch, -yaw, -tilt, R_inv);

    // Translate to center of element.
    T(state_t::PS_S,  6) = length/2e0*MtoMM;
    T(state_t::PS_PS, 6) = 1e0;
    inverse(T_inv, T);

    misalign_inv = prod(T, scl);
    misalign_inv = prod(R_inv, misalign_inv);
    misalign_inv = prod(T_inv, misalign_inv);
    misalign_inv = prod(scl_inv, misalign_inv);
}

void Moment2ElementBase::advance(StateBase& s)
{
    state_t&  ST = static_cast<state_t&>(s);
    using namespace boost::numeric::ublas;

    // IonEk is Es + E_state; the latter is set by user.
    ST.real.recalc();

    if(ST.real.IonEk!=last_Kenergy_in) {
        // need to re-calculate energy dependent terms

        recompute_matrix(ST); // updates transfer and last_Kenergy_out

        get_misalign(ST);

        scratch  = prod(transfer, misalign);
        transfer = prod(misalign_inv, scratch);

        ST.real.recalc();
    }

    // recompute_matrix only called when ST.IonEk != last_Kenergy_in.
    // Matrix elements are scaled with particle energy.

    ST.pos += length;

    ST.ref.phis   += ST.ref.SampleIonK*length*MtoMM;
    ST.real.phis  += ST.real.SampleIonK*length*MtoMM;
    ST.real.IonEk  = last_Kenergy_out;

    ST.moment0 = prod(transfer, ST.moment0);

    noalias(scratch)  = prod(transfer, ST.state);
    noalias(ST.state) = prod(scratch, trans(transfer));
}

void Moment2ElementBase::recompute_matrix(state_t& ST)
{
    // Default, for passive elements.

    transfer = transfer_raw;

    std::string t_name = type_name(); // C string -> C++ string.
//    if (t_name != "drift") {
//        // Scale matrix elements.
//        for (unsigned k = 0; k < 2; k++) {
//            transfer(2*k, 2*k+1) *= ST.bg_ref/ST.bg1;
//            transfer(2*k+1, 2*k) *= ST.bg1/ST.bg_ref;
//        }
//    }

    last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
}

namespace {

struct ElementSource : public Moment2ElementBase
{
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementSource(const Config& c): base_t(c), istate(c) {}

    void misalign1(state_t& ST);

    virtual void advance(StateBase& s)
    {
        state_t& ST = static_cast<state_t&>(s);
        // Replace state with our initial values
        ST.assign(istate);
    }

    virtual void show(std::ostream& strm) const
    {
        ElementVoid::show(strm);
        strm<<"Initial: "<<istate.state<<"\n";
    }

    state_t istate;
    // note that 'transfer' is not used by this element type

    virtual ~ElementSource() {}

    virtual const char* type_name() const {return "source";}
};

struct ElementMark : public Moment2ElementBase
{
    // Transport (identity) matrix for Marker.
    typedef Moment2ElementBase     base_t;
    typedef typename base_t::state_t state_t;

    ElementMark(const Config& c): base_t(c) {length = 0e0;}
    virtual ~ElementMark() {}
    virtual const char* type_name() const {return "marker";}
};

struct ElementBPM : public Moment2ElementBase
{
    // Transport (identity) matrix for BPM.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    Particle state;

    ElementBPM(const Config& c): base_t(c) {length = 0e0;}
    virtual ~ElementBPM() {}
    virtual const char* type_name() const {return "bpm";}
};

struct ElementDrift : public Moment2ElementBase
{
    // Transport matrix for Drift.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementDrift(const Config& c) : base_t(c) {}
    virtual ~ElementDrift() {}
    virtual const char* type_name() const {return "drift";}

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.
        transfer_raw = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

        double L = length*MtoMM; // Convert from [m] to [mm].

        transfer_raw(state_t::PS_X, state_t::PS_PX) = L;
        transfer_raw(state_t::PS_Y, state_t::PS_PY) = L;
        transfer_raw(state_t::PS_S, state_t::PS_PS) =
                -2e0*M_PI/(SampleLambda*ST.real.IonEs/MeVtoeV*cube(ST.real.bg))*L;

        transfer = transfer_raw;

        last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
    }
};

struct ElementOrbTrim : public Moment2ElementBase
{
    // Transport matrix for Orbit Trim.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementOrbTrim(const Config& c) : base_t(c) {length = 0e0;}
    virtual ~ElementOrbTrim() {}
    virtual const char* type_name() const {return "orbtrim";}

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.
        transfer = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

        double theta_x = conf().get<double>("theta_x", 0e0),
               theta_y = conf().get<double>("theta_y", 0e0);

        transfer = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
        transfer(state_t::PS_PX, 6) = theta_x*ST.ref.IonZ/ST.real.IonZ;
        transfer(state_t::PS_PY, 6) = theta_y*ST.ref.IonZ/ST.real.IonZ;

        last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
    }
};

struct ElementSBend : public Moment2ElementBase
{
    // Transport matrix for Gradient Sector Bend; with edge focusing (cylindrical coordinates).
    // Note, TLM only includes energy offset for the orbit; not the transport matrix.

    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;
    ElementSBend(const Config& c) : base_t(c) {}
    virtual ~ElementSBend() {}
    virtual const char* type_name() const {return "sbend";}

    virtual void advance(StateBase& s)
    {
        double    phis_temp, di_bg, Ek00, beta00, gamma00, IonK_Bend, dphis_temp;
        state_t&  ST = static_cast<state_t&>(s);
        using namespace boost::numeric::ublas;

        // IonEk is Es + E_state; the latter is set by user.
        ST.real.recalc();

        if(ST.real.IonEk!=last_Kenergy_in) {
            // need to re-calculate energy dependent terms

            recompute_matrix(ST); // updates transfer and last_Kenergy_out

            get_misalign(ST);

            scratch  = prod(transfer, misalign);
            transfer = prod(misalign_inv, scratch);

            ST.real.recalc();
        }

        // recompute_matrix only called when ST.IonEk != last_Kenergy_in.
        // Matrix elements are scaled with particle energy.

        ST.pos += length;

        phis_temp = ST.moment0[state_t::PS_S];

        ST.moment0 = prod(transfer, ST.moment0);

        ST.ref.phis   += ST.ref.SampleIonK*length*MtoMM;

        dphis_temp = ST.moment0[state_t::PS_S] - phis_temp;

        std::string HdipoleFitMode = conf().get<std::string>("HdipoleFitMode", "1");
        if (HdipoleFitMode != "1") {
            di_bg     = conf().get<double>("bg");
            // Dipole reference energy.
            Ek00      = (sqrt(sqr(di_bg)+1e0)-1e0)*ST.ref.IonEs;
            gamma00   = (Ek00+ST.ref.IonEs)/ST.ref.IonEs;
            beta00    = sqrt(1e0-1e0/sqr(gamma00));
            IonK_Bend = 2e0*M_PI/(beta00*SampleLambda);

            // J.B.: this is odd.
//            ST.real.phis  += IonK_Bend*length*MtoMM + dphis_temp;
            ST.real.phis  += ST.real.SampleIonK*length*MtoMM + dphis_temp;
        } else
            ST.real.phis  += ST.real.SampleIonK*length*MtoMM + dphis_temp;

        ST.real.IonEk  = last_Kenergy_out;

        noalias(scratch)  = prod(transfer, ST.state);
        noalias(ST.state) = prod(scratch, trans(transfer));
    }

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.

        transfer_raw = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

        double L     = conf().get<double>("L")*MtoMM,
               phi   = conf().get<double>("phi")*M_PI/180e0,
               phi1  = conf().get<double>("phi1")*M_PI/180e0,
               phi2  = conf().get<double>("phi2")*M_PI/180e0,
               K     = conf().get<double>("K", 0e0)/sqr(MtoMM),
               qmrel = (ST.real.IonZ-ST.ref.IonZ)/ST.ref.IonZ;

        std::string HdipoleFitMode = conf().get<std::string>("HdipoleFitMode", "1");
        if (HdipoleFitMode != "1") {
            double dip_bg    = conf().get<double>("bg"),
                   // Dipole reference energy.
                   dip_Ek    = (sqrt(sqr(dip_bg)+1e0)-1e0)*ST.ref.IonEs,
                   dip_gamma = (dip_Ek+ST.ref.IonEs)/ST.ref.IonEs,
                   dip_beta  = sqrt(1e0-1e0/sqr(dip_gamma)),
                   d         = (ST.ref.gamma-dip_gamma)/(sqr(dip_beta)*dip_gamma) - qmrel,
                   dip_IonK  = 2e0*M_PI/(dip_beta*SampleLambda);

            GetSBendMatrix(L, phi, phi1, phi2, K, ST.ref.IonEs, ST.ref.gamma, qmrel,
                           dip_beta, dip_gamma, d, dip_IonK, transfer_raw);
        } else
            GetSBendMatrix(L, phi, phi1, phi2, K, ST.ref.IonEs, ST.ref.gamma, qmrel,
                           ST.ref.beta, ST.ref.gamma, - qmrel, ST.ref.SampleIonK, transfer_raw);

        transfer = transfer_raw;

        last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
    }
};

struct ElementQuad : public Moment2ElementBase
{
    // Transport matrix for Quadrupole; K = B2/Brho.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementQuad(const Config& c) : base_t(c) {}
    virtual ~ElementQuad() {}
    virtual const char* type_name() const {return "quadrupole";}

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.
        transfer_raw = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

        double Brho  = ST.real.beta*(ST.real.IonEk+ST.real.IonEs)/(C0*ST.real.IonZ),
               K     = conf().get<double>("B2")/Brho/sqr(MtoMM),
               L     = conf().get<double>("L")*MtoMM,
               dx    = conf().get<double>("dx", 0e0)*MtoMM,
               dy    = conf().get<double>("dy", 0e0)*MtoMM,
               pitch = conf().get<double>("pitch", 0e0),
               yaw   = conf().get<double>("yaw", 0e0),
               tilt  = conf().get<double>("tilt", 0e0);

        // Horizontal plane.
        GetQuadMatrix(L,  K, (unsigned)state_t::PS_X, transfer_raw);
        // Vertical plane.
        GetQuadMatrix(L, -K, (unsigned)state_t::PS_Y, transfer_raw);
        // Longitudinal plane.
//        transfer_raw(state_t::PS_S, state_t::PS_S) = L;

        transfer_raw(state_t::PS_S, state_t::PS_PS) =
                -2e0*M_PI/(SampleLambda*ST.real.IonEs/MeVtoeV*cube(ST.real.bg))*L;

        transfer = transfer_raw;

        RotMat(dx, dy, pitch, yaw, tilt, misalign);
        inverse(misalign_inv, misalign);

        last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
    }
};

struct ElementSolenoid : public Moment2ElementBase
{
    // Transport (identity) matrix for a Solenoid; K = B/(2 Brho).
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementSolenoid(const Config& c) : base_t(c) {}
    virtual ~ElementSolenoid() {}
    virtual const char* type_name() const {return "solenoid";}

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.

        transfer = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

        double Brho  = ST.real.beta*(ST.real.IonEk+ST.real.IonEs)/(C0*ST.real.IonZ),
               K     = conf().get<double>("B")/(2e0*Brho)/MtoMM,
               L     = conf().get<double>("L")*MtoMM;

        GetSolMatrix(L, K, transfer);

        transfer(state_t::PS_S, state_t::PS_PS) =
                -2e0*M_PI/(SampleLambda*ST.real.IonEs/MeVtoeV*cube(ST.real.bg))*L;

        get_misalign(ST);

        scratch  = prod(transfer, misalign);
        transfer = prod(misalign_inv, scratch);

        last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
    }
};

struct ElementEDipole : public Moment2ElementBase
{
    // Transport matrix for Electrostatic Dipole with edge focusing.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementEDipole(const Config& c) : base_t(c) {}
    virtual ~ElementEDipole() {}
    virtual const char* type_name() const {return "edipole";}

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.

        transfer_raw = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

        value_mat R;

        bool   ver         = conf().get<double>("ver") == 1.0;
        double L           = conf().get<double>("L")*MtoMM,
               phi         = conf().get<double>("phi")*M_PI/180e0,
               // fit to TLM unit.            
               fringe_x    = conf().get<double>("fringe_x", 0e0)/MtoMM,
               fringe_y    = conf().get<double>("fringe_y", 0e0)/MtoMM,
               kappa       = conf().get<double>("asym_fac", 0e0),
               // spher: cylindrical - 0, spherical - 1.
               spher       = conf().get<double>("spher"),
               rho         = L/phi,
               eta0        = (ST.real.gamma-1e0)/2e0,
               // magnetic - 0, electrostatic - 1.
               h           = 1e0,
               Kx          = (1e0-spher+sqr(1e0+2e0*eta0))/sqr(rho),
               Ky          = spher/sqr(rho),
               dip_beta    = conf().get<double>("beta"),
               dip_gamma   = 1e0/sqrt(1e0-sqr(dip_beta)),
               delta_KZ    = ST.ref.IonZ/ST.real.IonZ - 1e0,
               SampleIonK  = 2e0*M_PI/(ST.real.beta*SampleLambda);

        GetEBendMatrix(L, phi, fringe_x, fringe_y, kappa, Kx, Ky, ST.ref.IonEs, ST.ref.beta, ST.real.gamma,
                       eta0, h, dip_beta, dip_gamma, delta_KZ, SampleIonK, transfer_raw);

        if (ver) {
            // Rotate transport matrix by 90 degrees.
            R = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
            R(state_t::PS_X,  state_t::PS_X)   =  0e0;
            R(state_t::PS_PX, state_t::PS_PX)  =  0e0;
            R(state_t::PS_Y,  state_t::PS_Y)   =  0e0;
            R(state_t::PS_PY, state_t::PS_PY)  =  0e0;
            R(state_t::PS_X,  state_t::PS_Y)   = -1e0;
            R(state_t::PS_PX, state_t::PS_PY)  = -1e0;
            R(state_t::PS_Y,  state_t::PS_X)   =  1e0;
            R(state_t::PS_PY,  state_t::PS_PX) =  1e0;

            scratch = prod(R, transfer_raw);
            transfer_raw = prod(scratch, trans(R));
        }

        transfer = transfer_raw;

        last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
    }
};

struct ElementEQuad : public Moment2ElementBase
{
    // Transport matrix for Electrostatic Quadrupole.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementEQuad(const Config& c) : base_t(c) {}
    virtual ~ElementEQuad() {}
    virtual const char* type_name() const {return "equad";}

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.
        // V0 [V] electrode voltage and R [m] electrode half-distance.
        transfer_raw = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

        double Brho = ST.real.beta*(ST.real.IonEk+ST.real.IonEs)/(C0*ST.real.IonZ),
               V0   = conf().get<double>("V"),
               R    = conf().get<double>("radius"),
               K    = 2e0*V0/(C0*ST.real.beta*sqr(R))/Brho/sqr(MtoMM),
               L    = conf().get<double>("L")*MtoMM;

        // Horizontal plane.
        GetQuadMatrix(L,  K, (unsigned)state_t::PS_X, transfer_raw);
        // Vertical plane.
        GetQuadMatrix(L, -K, (unsigned)state_t::PS_Y, transfer_raw);
        // Longitudinal plane.
//        transfer_raw(state_t::PS_S, state_t::PS_S) = L;

        transfer_raw(state_t::PS_S, state_t::PS_PS) =
                -2e0*M_PI/(SampleLambda*ST.real.IonEs/MeVtoeV*cube(ST.real.bg))*L;

        transfer = transfer_raw;

        last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
    }
};

struct ElementGeneric : public Moment2ElementBase
{
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementGeneric(const Config& c)
        :base_t(c)
    {
        std::vector<double> I = c.get<std::vector<double> >("transfer");
        if(I.size()>this->transfer_raw.data().size())
            throw std::invalid_argument("Initial transfer size too big");
        std::copy(I.begin(), I.end(), this->transfer_raw.data().begin());
    }
    virtual ~ElementGeneric() {}
    virtual const char* type_name() const {return "generic";}
};

} // namespace

void registerMoment2()
{
    Machine::registerState<Moment2State>("MomentMatrix2");

    Machine::registerElement<ElementSource                 >("MomentMatrix2", "source");

    Machine::registerElement<ElementMark                   >("MomentMatrix2", "marker");

    Machine::registerElement<ElementBPM                    >("MomentMatrix2", "bpm");

    Machine::registerElement<ElementDrift                  >("MomentMatrix2", "drift");

    Machine::registerElement<ElementOrbTrim                >("MomentMatrix2", "orbtrim");

    Machine::registerElement<ElementSBend                  >("MomentMatrix2", "sbend");

    Machine::registerElement<ElementQuad                   >("MomentMatrix2", "quadrupole");

    Machine::registerElement<ElementSolenoid               >("MomentMatrix2", "solenoid");

    Machine::registerElement<ElementRFCavity               >("MomentMatrix2", "rfcavity");

    Machine::registerElement<ElementStripper               >("MomentMatrix2", "stripper");

    Machine::registerElement<ElementEDipole                >("MomentMatrix2", "edipole");

    Machine::registerElement<ElementEQuad                  >("MomentMatrix2", "equad");

    Machine::registerElement<ElementGeneric                >("MomentMatrix2", "generic");
}
