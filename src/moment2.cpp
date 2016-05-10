
#include <fstream>

#include <limits>

#include <boost/numeric/ublas/lu.hpp>
#include <boost/lexical_cast.hpp>

#include "scsi/constants.h"
#include "scsi/moment2.h"

#include "scsi/rf_cavity.h"

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

void Moment2ElementBase::advance(StateBase& s)
{
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
    if (t_name != "rfcavity") {
        ST.ref.phis   += ST.ref.SampleIonK*length*MtoMM;
        ST.real.phis  += ST.real.SampleIonK*length*MtoMM;
        ST.real.IonEk  = last_Kenergy_out;
    }

    ST.moment0 = prod(transfer, ST.moment0);

    if (t_name == "rfcavity") {
        ST.moment0[state_t::PS_S]  = ST.real.phis - ST.ref.phis;
        ST.moment0[state_t::PS_PS] = (ST.real.IonEk-ST.ref.IonEk)/MeVtoeV;
    }

    noalias(scratch) = prod(transfer, ST.state);
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

void GetEdgeMatrix(const double rho, const double phi, typename Moment2ElementBase::value_t &M)
{
    typedef typename Moment2ElementBase::state_t state_t;

    M = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

    M(state_t::PS_PX, state_t::PS_X) =  tan(phi)/rho;
    M(state_t::PS_PY, state_t::PS_Y) = -tan(phi)/rho;
}

void GetQuadMatrix(const double L, const double K, const unsigned ind, typename Moment2ElementBase::value_t &M)
{
    // 2D quadrupole transport matrix.
    double sqrtK,
           psi,
           cs,
           sn;

    if (K > 0e0) {
        // Focusing.
        sqrtK = sqrt(K);
        psi = sqrtK*L;
        cs = ::cos(psi);
        sn = ::sin(psi);

        M(ind, ind) = M(ind+1, ind+1) = cs;
        if (sqrtK != 0e0)
            M(ind, ind+1) = sn/sqrtK;
        else
            M(ind, ind+1) = L;
        if (sqrtK != 0e0)
            M(ind+1, ind) = -sqrtK*sn;
        else
            M(ind+1, ind) = 0e0;
    } else {
        // Defocusing.
        sqrtK = sqrt(-K);
        psi = sqrtK*L;
        cs = ::cosh(psi);
        sn = ::sinh(psi);

        M(ind, ind) = M(ind+1, ind+1) = cs;
        if (sqrtK != 0e0)
            M(ind, ind+1) = sn/sqrtK;
        else
            M(ind, ind+1) = L;
        if (sqrtK != 0e0)
            M(ind+1, ind) = sqrtK*sn;
        else
            M(ind+1, ind) = 0e0;
    }
}

void GetSolMatrix(const double L, const double K, typename Moment2ElementBase::value_t &M)
{
    typedef typename Moment2ElementBase::state_t state_t;

    double C = ::cos(K*L),
           S = ::sin(K*L);

    M(state_t::PS_X, state_t::PS_X)
            = M(state_t::PS_PX, state_t::PS_PX)
            = M(state_t::PS_Y, state_t::PS_Y)
            = M(state_t::PS_PY, state_t::PS_PY)
            = sqr(C);

    if (K != 0e0)
        M(state_t::PS_X, state_t::PS_PX) = S*C/K;
    else
        M(state_t::PS_X, state_t::PS_PX) = L;
    M(state_t::PS_X, state_t::PS_Y) = S*C;
    if (K != 0e0)
        M(state_t::PS_X, state_t::PS_PY) = sqr(S)/K;
    else
        M(state_t::PS_X, state_t::PS_PY) = 0e0;

    M(state_t::PS_PX, state_t::PS_X) = -K*S*C;
    M(state_t::PS_PX, state_t::PS_Y) = -K*sqr(S);
    M(state_t::PS_PX, state_t::PS_PY) = S*C;

    M(state_t::PS_Y, state_t::PS_X) = -S*C;
    if (K != 0e0)
        M(state_t::PS_Y, state_t::PS_PX) = -sqr(S)/K;
    else
        M(state_t::PS_Y, state_t::PS_PX) = 0e0;
    if (K != 0e0)
        M(state_t::PS_Y, state_t::PS_PY) = S*C/K;
    else
        M(state_t::PS_Y, state_t::PS_PY) = L;

    M(state_t::PS_PY, state_t::PS_X) = K*sqr(S);
    M(state_t::PS_PY, state_t::PS_PX) = -S*C;
    M(state_t::PS_PY, state_t::PS_Y) = -K*S*C;

    // Longitudinal plane.
    // For total path length.
//        M(state_t::PS_S, state_t::PS_S) = L;
}

struct ElementSource : public Moment2ElementBase
{
    typedef Moment2ElementBase base_t;
    typedef typename base_t::state_t state_t;
    ElementSource(const Config& c)
        :base_t(c), istate(c)
    {}

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
    // Transport (identity) matrix for a Marker.
    typedef Moment2ElementBase base_t;
    typedef typename base_t::state_t state_t;
    ElementMark(const Config& c) :base_t(c){
        length = 0e0;
    }
    virtual ~ElementMark() {}
    virtual const char* type_name() const {return "marker";}
};

struct ElementDrift : public Moment2ElementBase
{
    // Transport matrix for a Drift.
    typedef Moment2ElementBase base_t;
    typedef typename base_t::state_t state_t;
    ElementDrift(const Config& c)
        :base_t(c)
    {
    }
    virtual ~ElementDrift() {}

    virtual void recompute_matrix(state_t& ST)
    {
        double L = length*MtoMM; // Convert from [m] to [mm].

        this->transfer_raw(state_t::PS_X, state_t::PS_PX) = L;
        this->transfer_raw(state_t::PS_Y, state_t::PS_PY) = L;
        transfer_raw(state_t::PS_S, state_t::PS_PS) =
                -2e0*M_PI/(SampleLambda*ST.real.IonEs/MeVtoeV*cube(ST.real.bg))*L;

        transfer = transfer_raw;

        last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
    }

    virtual const char* type_name() const {return "drift";}
};

struct ElementSBend : public Moment2ElementBase
{
    // Transport matrix for a Gradient Sector Bend (cylindrical coordinates).

    // *** Add entrance and exit angles.

    typedef Moment2ElementBase base_t;
    typedef typename base_t::state_t state_t;
    ElementSBend(const Config& c)
        :base_t(c)
    {
    }
    virtual ~ElementSBend() {}

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
        GetQuadMatrix(L, Kx, (unsigned)state_t::PS_X, this->transfer_raw);
        // Vertical plane.
        GetQuadMatrix(L, Ky, (unsigned)state_t::PS_Y, this->transfer_raw);

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

        this->transfer_raw(state_t::PS_X,  state_t::PS_PS) = dx/(rho*sqr(ST.ref.beta)*ST.ref.gamma*ST.ref.IonEs/MeVtoeV);
        this->transfer_raw(state_t::PS_PX, state_t::PS_PS) = sx/(rho*sqr(ST.ref.beta)*ST.ref.gamma*ST.ref.IonEs/MeVtoeV);
        this->transfer_raw(state_t::PS_S,  state_t::PS_X)  = sx/rho*ST.ref.SampleIonK;
        this->transfer_raw(state_t::PS_S,  state_t::PS_PX) = dx/rho*ST.ref.SampleIonK;
        // Low beta approximation.
        this->transfer_raw(state_t::PS_S,  state_t::PS_PS) =
                ((L-sx)/(Kx*sqr(rho))-L/sqr(ST.ref.gamma))*ST.ref.SampleIonK
                /(sqr(ST.ref.beta)*ST.ref.gamma*ST.ref.IonEs/MeVtoeV);

        double qmrel = (ST.real.IonZ-ST.ref.IonZ)/ST.ref.IonZ;

        // Add dipole terms.
        this->transfer_raw(state_t::PS_X,  6) = -dx/rho*qmrel;
        this->transfer_raw(state_t::PS_PX, 6) = -sx/rho*qmrel;
        // Check expression.
        this->transfer_raw(state_t::PS_S,  6) =
                -((L-sx)/(Kx*sqr(rho))-L/sqr(ST.ref.gamma)+L/sqr(ST.ref.gamma))*ST.ref.SampleIonK*qmrel;

        // Edge focusing.
        GetEdgeMatrix(rho, phi2, edge2);

        transfer_raw = prod(transfer_raw, edge1);
        transfer_raw = prod(edge2, transfer_raw);

        // Longitudinal plane.
        // For total path length.
//        this->transfer_raw(state_t::PS_S,  state_t::PS_S) = L;

        transfer = transfer_raw;

        last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
    }

    virtual const char* type_name() const {return "sbend";}
};

struct ElementQuad : public Moment2ElementBase
{
    // Transport matrix for a Quadrupole; K = B2/Brho.
    typedef Moment2ElementBase base_t;
    typedef typename base_t::state_t state_t;
    ElementQuad(const Config& c)
        :base_t(c)
    {
    }
    virtual ~ElementQuad() {}

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.
        this->transfer_raw = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

        double Brho = ST.real.beta*(ST.real.IonEk+ST.real.IonEs)/(C0*ST.real.IonZ),
               K    = conf().get<double>("B2")/Brho/sqr(MtoMM),
               L    = conf().get<double>("L")*MtoMM;

        // Horizontal plane.
        GetQuadMatrix(L,  K, (unsigned)state_t::PS_X, this->transfer_raw);
        // Vertical plane.
        GetQuadMatrix(L, -K, (unsigned)state_t::PS_Y, this->transfer_raw);
        // Longitudinal plane.
//        this->transfer_raw(state_t::PS_S, state_t::PS_S) = L;

        transfer_raw(state_t::PS_S, state_t::PS_PS) =
                -2e0*M_PI/(SampleLambda*ST.real.IonEs/MeVtoeV*cube(ST.real.bg))*L;

        transfer = transfer_raw;

        last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
    }

    virtual const char* type_name() const {return "quadrupole";}
};

struct ElementSolenoid : public Moment2ElementBase
{
    // Transport (identity) matrix for a Solenoid; K = B0/(2 Brho).
    typedef Moment2ElementBase base_t;
    typedef typename base_t::state_t state_t;
    ElementSolenoid(const Config& c)
        :base_t(c)
    {
    }
    virtual ~ElementSolenoid() {}
    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.
        this->transfer_raw = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

        double Brho = ST.real.beta*(ST.real.IonEk+ST.real.IonEs)/(C0*ST.real.IonZ),
               K    = conf().get<double>("B")/(2e0*Brho)/MtoMM,
               L    = conf().get<double>("L")*MtoMM;      // Convert from [m] to [mm].

        GetSolMatrix(L, K, this->transfer_raw);

        transfer_raw(state_t::PS_S, state_t::PS_PS) =
                -2e0*M_PI/(SampleLambda*ST.real.IonEs/MeVtoeV*cube(ST.real.bg))*L;

        transfer = transfer_raw;

        last_Kenergy_in = last_Kenergy_out = ST.real.IonEk; // no energy gain
    }

    virtual const char* type_name() const {return "solenoid";}
};


struct ElementStripper : public Moment2ElementBase
{
    // Transport (identity) matrix for a Charge Stripper.
    typedef Moment2ElementBase base_t;
    typedef typename base_t::state_t state_t;
    ElementStripper(const Config& c)
        :base_t(c)
    {
        // Identity matrix.
    }
    virtual ~ElementStripper() {}

    virtual const char* type_name() const {return "stripper";}
};

struct ElementEDipole : public Moment2ElementBase
{
    // Transport matrix for an Electric Dipole.
    typedef Moment2ElementBase base_t;
    typedef typename base_t::state_t state_t;
    ElementEDipole(const Config& c)
        :base_t(c)
    {
        //double L = c.get<double>("L")*MtoMM;

    }
    virtual ~ElementEDipole() {}

    virtual const char* type_name() const {return "edipole";}
};

struct ElementGeneric : public Moment2ElementBase
{
    typedef Moment2ElementBase base_t;
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

    Machine::registerElement<ElementDrift                  >("MomentMatrix2",   "drift");

    Machine::registerElement<ElementSBend                  >("MomentMatrix2",   "sbend");

    Machine::registerElement<ElementQuad                   >("MomentMatrix2",   "quadrupole");

    Machine::registerElement<ElementSolenoid               >("MomentMatrix2",   "solenoid");

    Machine::registerElement<ElementRFCavity               >("MomentMatrix2",   "rfcavity");

    Machine::registerElement<ElementStripper               >("MomentMatrix2",   "stripper");

    Machine::registerElement<ElementEDipole                >("MomentMatrix2",   "edipole");

    Machine::registerElement<ElementGeneric                >("MomentMatrix2",   "generic");
}
