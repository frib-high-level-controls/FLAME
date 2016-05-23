
#include <fstream>

#include <limits>

#include <boost/numeric/ublas/lu.hpp>
#include <boost/lexical_cast.hpp>

#include "scsi/constants.h"
#include "scsi/moment2.h"

#include "scsi/moment2.h"
#include "scsi/moment2_sup.h"
#include "scsi/rf_cavity.h"
#include "scsi/chg_stripper.h"

#include "scsi/h5loader.h"


namespace {
// http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
// by LU-decomposition.
void inverse(Moment2ElementBase::value_t& out, const Moment2ElementBase::value_t& in)
{
    using boost::numeric::ublas::permutation_matrix;
    using boost::numeric::ublas::lu_factorize;
    using boost::numeric::ublas::lu_substitute;
    using boost::numeric::ublas::identity_matrix;

    Moment2ElementBase::value_t scratch(in); // copy
    permutation_matrix<size_t> pm(scratch.size1());
    if(lu_factorize(scratch, pm)!=0)
        throw std::runtime_error("Failed to invert matrix");
    out.assign(identity_matrix<double>(scratch.size1()));
    //out = identity_matrix<double>(scratch.size1());
    lu_substitute(scratch, pm, out);
}

} // namespace

std::ostream& operator<<(std::ostream& strm, const Particle& P)
{
    strm <<std::setprecision(8)<<std::setw(14)
      <<"IonZ="<<P.IonZ
      <<" IonQ="<<P.IonQ
      <<" IonEs="<<P.IonEs
      <<" IonEk="<<P.IonEk
      <<" SampleIonK="<<P.SampleIonK
      <<" phis="<<P.phis
      <<" IonW="<<P.IonW
      <<" gamma="<<P.gamma
      <<" beta="<<P.beta
      <<" bg="<<P.bg
      ;
    return strm;
}

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

        std::vector<double> nchg = c.get<std::vector<double> >("NCharge");
        if(nchg.size()!=ics.size())
            throw std::invalid_argument("NCharge[] and IonChargeStates[] must have equal length");
        ref.IonQ = nchg[icstate];

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
    ref.recalc();

    real           = ref;

    real.phis      = moment0[PS_S];
    real.IonEk    += moment0[PS_PS]*MeVtoeV;

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
        Info.name = "ref_IonQ";
        Info.ptr = &ref.IonQ;
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
    ,transfer(boost::numeric::ublas::identity_matrix<double>(state_t::maxsize))
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
    misalign_inv = O->misalign_inv;
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

void Moment2ElementBase::advance(StateBase& s)
{
    double   phis_temp;
    state_t& ST = static_cast<state_t&>(s);
    using namespace boost::numeric::ublas;

    // IonEk is Es + E_state; the latter is set by user.
    ST.real.recalc();

    if(ST.real.IonEk!=last_Kenergy_in) {
        // need to re-calculate energy dependent terms

        recompute_matrix(ST); // updates transfer and last_Kenergy_out

        noalias(scratch)  = prod(misalign, transfer);
        noalias(transfer) = prod(scratch, misalign_inv);

        ST.real.recalc();
    }

    // recompute_matrix only called when ST.IonEk != last_Kenergy_in.
    // Matrix elements are scaled with particle energy.

    ST.pos += length;

    std::string t_name = type_name(); // C string -> C++ string.
    if ((t_name != "rfcavity") && (t_name != "sbend")) {
        ST.ref.phis   += ST.ref.SampleIonK*length*MtoMM;
        ST.real.phis  += ST.real.SampleIonK*length*MtoMM;
        ST.real.IonEk  = last_Kenergy_out;
    } else if (t_name == "sbend")
        phis_temp = ST.moment0[state_t::PS_S];

    ST.moment0 = prod(transfer, ST.moment0);

    if (t_name == "rfcavity") {
        ST.moment0[state_t::PS_S]  = ST.real.phis - ST.ref.phis;
        ST.moment0[state_t::PS_PS] = (ST.real.IonEk-ST.ref.IonEk)/MeVtoeV;
    } else if (t_name == "sbend") {
        ST.ref.phis   += ST.ref.SampleIonK*length*MtoMM;

        double dphis_temp = ST.moment0[state_t::PS_S] - phis_temp;
         ST.real.phis  += ST.real.SampleIonK*length*MtoMM + dphis_temp;

        ST.real.IonEk  = last_Kenergy_out;
    }

    noalias(scratch) = prod(transfer, ST.state);
    noalias(ST.state) = prod(scratch, trans(transfer));
}

void Moment2ElementBase::recompute_matrix(state_t& ST)
{
    // Default no-op
    transfer = transfer_raw;

    last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
}

namespace {

struct ElementSource : public Moment2ElementBase
{
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementSource(const Config& c): base_t(c), istate(c) {}

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

    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;
    ElementSBend(const Config& c) : base_t(c) {}
    virtual ~ElementSBend() {}
    virtual const char* type_name() const {return "sbend";}

    virtual void recompute_matrix(state_t& ST)
    {
        double L    = conf().get<double>("L")*MtoMM,
               phi  = conf().get<double>("phi")*M_PI/180e0,
               phi1 = conf().get<double>("phi1")*M_PI/180e0,
               phi2 = conf().get<double>("phi2")*M_PI/180e0,
               rho  = L/phi,
               K    = conf().get<double>("K", 0e0)/sqr(MtoMM),
               Kx   = K + 1e0/sqr(rho),
               Ky   = -K,
               dx   = 0e0,
               sx   = 0e0;

        typename Moment2ElementBase::value_t edge1, edge2;

        // Edge focusing.
        GetEdgeMatrix(rho, phi1, edge1);
        // Horizontal plane.
        GetQuadMatrix(L, Kx, (unsigned)state_t::PS_X, transfer_raw);
        // Vertical plane.
        GetQuadMatrix(L, Ky, (unsigned)state_t::PS_Y, transfer_raw);

        // Include dispersion.
        if (Kx == 0e0) {
            dx = sqr(L)/2e0;
            sx = L;
        } else if (Kx > 0e0) {
            dx = (1e0-cos(sqrt(Kx)*L))/Kx;
            sx = sin(sqrt(Kx)*L)/sqrt(Kx);
        } else {
            dx = (1e0-cosh(sqrt(-Kx)*L))/Kx;
            sx = sin(sqrt(Kx)*L)/sqrt(Kx);
        }

        transfer_raw(state_t::PS_X,  state_t::PS_PS) = dx/(rho*sqr(ST.ref.beta)*ST.ref.gamma*ST.ref.IonEs/MeVtoeV);
        transfer_raw(state_t::PS_PX, state_t::PS_PS) = sx/(rho*sqr(ST.ref.beta)*ST.ref.gamma*ST.ref.IonEs/MeVtoeV);
        transfer_raw(state_t::PS_S,  state_t::PS_X)  = sx/rho*ST.ref.SampleIonK;
        transfer_raw(state_t::PS_S,  state_t::PS_PX) = dx/rho*ST.ref.SampleIonK;
        // Low beta approximation.
        transfer_raw(state_t::PS_S,  state_t::PS_PS) =
                ((L-sx)/(Kx*sqr(rho))-L/sqr(ST.ref.gamma))*ST.ref.SampleIonK
                /(sqr(ST.ref.beta)*ST.ref.gamma*ST.ref.IonEs/MeVtoeV);

        double qmrel = (ST.real.IonZ-ST.ref.IonZ)/ST.ref.IonZ;

        // Add dipole terms.
        transfer_raw(state_t::PS_X,  6) = -dx/rho*qmrel;
        transfer_raw(state_t::PS_PX, 6) = -sx/rho*qmrel;
        transfer_raw(state_t::PS_S,  6) = -(L-sx)/(Kx*sqr(rho))*ST.ref.SampleIonK*qmrel;

        // Edge focusing.
        GetEdgeMatrix(rho, phi2, edge2);

        transfer_raw = prod(transfer_raw, edge1);
        transfer_raw = prod(edge2, transfer_raw);

        // Longitudinal plane.
        // For total path length.
//        transfer_raw(state_t::PS_S,  state_t::PS_S) = L;

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

        double Brho = ST.real.beta*(ST.real.IonEk+ST.real.IonEs)/(C0*ST.real.IonZ),
               K    = conf().get<double>("B2")/Brho/sqr(MtoMM),
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
        transfer_raw = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

        double Brho = ST.real.beta*(ST.real.IonEk+ST.real.IonEs)/(C0*ST.real.IonZ),
               K    = conf().get<double>("B")/(2e0*Brho)/MtoMM,
               L    = conf().get<double>("L")*MtoMM;      // Convert from [m] to [mm].

        GetSolMatrix(L, K, transfer_raw);

        transfer_raw(state_t::PS_S, state_t::PS_PS) =
                -2e0*M_PI/(SampleLambda*ST.real.IonEs/MeVtoeV*cube(ST.real.bg))*L;

        transfer = transfer_raw;

        last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
    }
};

struct ElementEDipole : public Moment2ElementBase
{
    // Transport matrix for Electrostatic Dipole; with edge focusing.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementEDipole(const Config& c) : base_t(c) {}
    virtual ~ElementEDipole() {}
    virtual const char* type_name() const {return "edipole";}

    virtual void recompute_matrix(state_t& ST)
    {
        //double L = c.get<double>("L")*MtoMM;

//        mscpTracker.position=fribnode.position+fribnode.length;
//        double gamma=(Ek[ii_state]+FRIBPara.ionEs)/FRIBPara.ionEs;
//        double beta=Math.sqrt(1.0-1.0/(gamma*gamma));
//        double SampleionK=2*Math.PI/(beta*SampleLamda*1e3); //rad/mm, marking the currut using ionK
//        double length=fribnode.length; // m
//        double bangle=fribnode.attribute[1]; // degree
//        double beta00=fribnode.attribute[2]; // The voltage is set that the particle with beta_EB would be bended ideally
//        double n1=fribnode.attribute[3]; // 0 is for cylindrical, and 1 is for spherical
//        int vertical=(int) fribnode.attribute[4]; // 0 for horizontal, 1 is for vertical ebend matrix
//        double fringeX=fribnode.attribute[5]; // X fringe
//        double fringeY=fribnode.attribute[6]; // Y fringe
//        double kappa=(int) fribnode.attribute[7]; // voltage asymetry factor
//        double gamma0=tlmPara.Gamma_tab.get(lattcnt+1)[1];
//        double beta0=tlmPara.Beta_tab.get(lattcnt+1)[1];
//        double h=1.0; // 0 for magnetic and 1 for electrostatic
//        double rh0=length/(bangle/180*Math.PI);

//        double Fy_temp=TransVector[ii_state].getElem(4);
//        double dFy_temp=TransVector[ii_state].getElem(4)-Fy_temp;
//        Fy_abs[ii_state]=Fy_abs[ii_state]+SampleionK*1000*fribnode.length+dFy_temp;

    }
   virtual void assign(const ElementVoid *other)
   {
        throw std::logic_error("assign/reconfigure() not implemented for rf cavity");
        // TODO: can't copy 'inf'.  Need to parse once in ctor.
//        const ElementRFCavity *O = static_cast<const ElementRFCavity*>(other);
//        Moment2ElementBase::assign(O);
//        CavData = O->CavData;
//        CavTLMLineTab = O->CavTLMLineTab;
//        phi_ref = O->phi_ref;
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

    Machine::registerElement<ElementSource                 >("MomentMatrix2",   "source");

    Machine::registerElement<ElementMark                   >("MomentMatrix2",   "marker");

    Machine::registerElement<ElementBPM                    >("MomentMatrix2",   "bpm");

    Machine::registerElement<ElementDrift                  >("MomentMatrix2",   "drift");

    Machine::registerElement<ElementOrbTrim                >("MomentMatrix2",   "orbtrim");

    Machine::registerElement<ElementSBend                  >("MomentMatrix2",   "sbend");

    Machine::registerElement<ElementQuad                   >("MomentMatrix2",   "quadrupole");

    Machine::registerElement<ElementSolenoid               >("MomentMatrix2",   "solenoid");

    Machine::registerElement<ElementRFCavity               >("MomentMatrix2",   "rfcavity");

    Machine::registerElement<ElementStripper               >("MomentMatrix2",   "stripper");

    Machine::registerElement<ElementEDipole                >("MomentMatrix2",   "edipole");

    Machine::registerElement<ElementEDipole                >("MomentMatrix2",   "equad");

    Machine::registerElement<ElementGeneric                >("MomentMatrix2",   "generic");
}
