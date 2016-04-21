
#include <fstream>

#include <limits>

#include <boost/numeric/ublas/lu.hpp>

#include "scsi/constants.h"
#include "scsi/moment2.h"

#include "scsi/h5loader.h"

void CavDataType::RdData(const std::string &FileName)
{
    std::string       line;
    double            s, Elong;
    std::stringstream str;
    std::fstream      inf;

    inf.open(FileName.c_str(), std::ifstream::in);
    if (!inf.is_open()) {
        std::cerr << "*** RdData: failed to open " << FileName << "\n";
        exit(1);
    }
    while (getline(inf, line) && !inf.fail()) {
        str.str(line);
        str >> s >> Elong;
        this->s.push_back(s), this->Elong.push_back(Elong);
        // Convert from [mm] to [m].
//        this->s[this->s.size()-1] /= MtoMM;
    }
    inf.close();

    if (false) {
        std::cout << "\n";
        for (size_t k = 0; k < this->s.size(); k++)
            this->show(std::cout, k);
    }
}


void CavDataType::show(std::ostream& strm, const int k) const
{
    strm << std::scientific << std::setprecision(5)
         << std::setw(13) << this->s[k] << std::setw(13) << this->Elong[k] << "\n";
}


void CavDataType::show(std::ostream& strm) const
{
    for (unsigned int k = 0; k < this->s.size(); k++)
        this->show(strm, k);
}


class CavTLMLineType {
public:
    std::vector<double> s;         // Longitudinal position [m].
    std::vector<std::string> Elem;
    std::vector<double> E0,
                        T,
                        S,
                        Accel;

    void clear(void);
    void set(const double, const std::string &, const double,
             const double, const double, const double);
    void show(std::ostream& strm, const int) const;
    void show(std::ostream& strm) const;
};


void CavTLMLineType::clear(void)
{
    this->s.clear(); this->Elem.clear(); this->E0.clear();
    this->T.clear(); this->S.clear(); this->Accel.clear();
}


void CavTLMLineType::set(const double s, const std::string &Elem, const double E0,
                         const double T, const double S, const double Accel)
{
    this->s.push_back(s); this->Elem.push_back(Elem); this->E0.push_back(E0);
    this->T.push_back(T); this->S.push_back(S); this->Accel.push_back(Accel);
}


void CavTLMLineType::show(std::ostream& strm, const int k) const
{
    strm << std::fixed << std::setprecision(5)
         << std::setw(9) << this->s[k] << std::setw(10) << this->Elem[k]
         << std::setw(9) << this->T[k] << std::setw(9) << this->S[k]
         << std::setw(9) << this->Accel[k] << "\n";
}


void CavTLMLineType::show(std::ostream& strm) const
{
    for (unsigned int k = 0; k < this->s.size(); k++)
        this->show(strm, k);
}


// Long. sampling frequency [Hz]; must be set to RF Cavity frequency.
# define SampleFreq   80.5e6
// Sampling distance [m].
# define SampleLambda (C0/SampleFreq*MtoMM)

const int MaxNCav    = 5;

CavDataType         CavData2[MaxNCav];
std::stringstream   CavTLMstream2[MaxNCav];
CavTLMLineType      CavTLMLineTab2[MaxNCav];


void PrtVec(const std::vector<double> &a)
{
    for (size_t k = 0; k < a.size(); k++)
        std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << a[k];
    std::cout << "\n";
}


typedef boost::numeric::ublas::vector<double> value_vec;
typedef boost::numeric::ublas::matrix<double> value_mat;

enum {maxsize = 7};

typedef boost::numeric::ublas::vector<double,
                boost::numeric::ublas::bounded_array<double, maxsize>
> vector_t;

typedef boost::numeric::ublas::matrix<double,
                boost::numeric::ublas::row_major,
                boost::numeric::ublas::bounded_array<double, maxsize*maxsize>
> matrix_t;


void PrtVec(const vector_t &a)
{
    for (size_t k = 0; k < a.size(); k++)
        std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << a[k];
    std::cout << "\n";
}


void PrtMat(const value_mat &M)
{
    for (size_t j = 0; j < M.size1(); j++) {
        for (size_t k = 0; k < M.size2(); k++)
            std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << M(j, k);
        std::cout << "\n";
    }
}


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
//    ,pos(c.get<double>("L", 0e0))
    ,pos(0e0)
    ,FyAbs(0e0)
    ,Fy_absState(0e0)
    ,EkState(0e0)
    ,moment0(maxsize, 0e0)
    ,state(boost::numeric::ublas::identity_matrix<double>(maxsize))
{
    try{
        const std::vector<double>& I = c.get<std::vector<double> >("moment0");
        if(I.size()>moment0.size())
            throw std::invalid_argument("Initial moment0 size too big");
        std::copy(I.begin(), I.end(), moment0.begin());
    }catch(key_error&){
        // default to zeros
    }catch(boost::bad_any_cast&){
        throw std::invalid_argument("'initial' has wrong type (must be vector)");
    }

    try{
        const std::vector<double>& I = c.get<std::vector<double> >("initial");
        if(I.size()>state.data().size())
            throw std::invalid_argument("Initial state size too big");
        std::copy(I.begin(), I.end(), state.data().begin());
    }catch(key_error&){
        // default to identity
    }catch(boost::bad_any_cast&){
        throw std::invalid_argument("'initial' has wrong type (must be vector)");
    }

    double Erest     = c.get<double>("IonEs", 0e0), // Rest energy.
           Ekinetic0 = c.get<double>("IonEk", 0e0);
           IonZ      = c.get<double>("IonZ", 0e0);

    gamma_ref      = (Erest+Ekinetic0)/Erest;
    beta_ref       = sqrt(1e0-1e0/sqr(gamma_ref));
    bg_ref         = beta_ref*gamma_ref;

    SampleIonK_ref = 2e0*M_PI/(beta_ref*SampleLambda);

    Ekinetic  = Ekinetic0;
}

Moment2State::~Moment2State() {}

Moment2State::Moment2State(const Moment2State& o, clone_tag t)
    :StateBase(o, t)
    ,pos(o.pos)
    ,Ekinetic(o.Ekinetic)
    ,SampleIonK_ref(o.SampleIonK_ref)
    ,SampleIonK(o.SampleIonK)
    ,moment0(o.moment0)
    ,state(o.state)
{}

void Moment2State::assign(const StateBase& other)
{
    const Moment2State *O = dynamic_cast<const Moment2State*>(&other);
    if(!O)
        throw std::invalid_argument("Can't assign State: incompatible types");
    pos             = O->pos;
    Ekinetic        = O->Ekinetic;
    SampleIonK_ref  = O->SampleIonK_ref;
    SampleIonK      = O->SampleIonK;
    FyAbs           = O->FyAbs;
    Fy_absState     = O->Fy_absState;
    EkState         = O->EkState;
    gamma_ref       = O->gamma_ref;
    beta_ref        = O->beta_ref;
    gamma           = O->gamma;
    beta            = O->beta;
    bg_ref          = O->bg_ref;
    bg1             = O->bg1;
    moment0         = O->moment0;
    state           = O->state;
    StateBase::assign(other);
}

void Moment2State::show(std::ostream& strm) const
{
    strm<<"State: energy="<<Ekinetic<<" moment0="<<moment0<<" state="<<state<<"\n";
}

bool Moment2State::getArray(unsigned idx, ArrayInfo& Info) {
    if(idx==0) {
        Info.name = "state";
        Info.ptr = &state(0,0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 2;
        Info.dim[0] = state.size1();
        Info.dim[1] = state.size2();
        return true;
    } else if(idx==1) {
        Info.name = "moment0";
        Info.ptr = &moment0(0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim[0] = moment0.size();
        return true;
    } else if(idx==2) {
        Info.name = "pos";
        Info.ptr = &pos;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==3) {
        Info.name = "FyAbs";
        Info.ptr = &FyAbs;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==4) {
        Info.name = "Fy_absState";
        Info.ptr = &Fy_absState;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==5) {
        Info.name = "EkState";
        Info.ptr = &EkState;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==6) {
        Info.name = "Ekinetic";
        Info.ptr = &Ekinetic;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==7) {
        Info.name = "SampleIonK_ref";
        Info.ptr = &SampleIonK_ref;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==8) {
        Info.name = "SampleIonK";
        Info.ptr = &SampleIonK;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==9) {
        Info.name = "gamma_ref";
        Info.ptr = &gamma_ref;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==10) {
        Info.name = "beta_ref";
        Info.ptr = &beta_ref;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==11) {
        Info.name = "gamma";
        Info.ptr = &gamma;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==12) {
        Info.name = "beta";
        Info.ptr = &beta;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==13) {
        Info.name = "bg_ref";
        Info.ptr = &bg_ref;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==14) {
        Info.name = "bg1";
        Info.ptr = &bg1;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    }
    return StateBase::getArray(idx-15, Info);
}

Moment2ElementBase::Moment2ElementBase(const Config& c)
    :ElementVoid(c)
    ,transfer(state_t::maxsize, state_t::maxsize)
    ,transfer_raw(boost::numeric::ublas::identity_matrix<double>(state_t::maxsize))
    ,misalign(boost::numeric::ublas::identity_matrix<double>(state_t::maxsize))
    ,misalign_inv(state_t::maxsize, state_t::maxsize)
    ,scratch(state_t::maxsize, state_t::maxsize)
{
    length = c.get<double>("L", 0e0);
    Erest = c.get<double>("IonEs");

    inverse(misalign_inv, misalign);

    // spoil to force recalculation of energy dependent terms
    last_Kenergy_in = last_Kenergy_out = std::numeric_limits<double>::quiet_NaN();
}

Moment2ElementBase::~Moment2ElementBase() {}

void Moment2ElementBase::assign(const ElementVoid *other)
{
    const Moment2ElementBase *O = static_cast<const Moment2ElementBase*>(other);
    length = O->length;
    FSampLength = O->FSampLength;
    Erest = O->Erest;
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
          "FSampLength "<<FSampLength<<"\n"
          "Erest "<<Erest<<"\n"
          "Transfer: "<<transfer<<"\n"
          "Transfer Raw: "<<transfer_raw<<"\n"
          "Mis-align: "<<misalign<<"\n";
}

void Moment2ElementBase::advance(StateBase& s)
{
    state_t& ST = static_cast<state_t&>(s);
    using namespace boost::numeric::ublas;

    if(ST.Ekinetic!=last_Kenergy_in) {

        // need to re-calculate energy dependent terms

        recompute_matrix(ST); // updates transfer and last_Kenergy_out

        noalias(scratch) = prod(misalign, transfer);
        noalias(transfer) = prod(scratch, misalign_inv);
    }

    ST.pos += length;

    // recompute_matrix only called when ST.Ekinetic != last_Kenergy_in.
    // Matrix elements are scaled with particle energy.
    ST.gamma          = (Erest+ST.EkState)/Erest;   // Approximate (E_k = m0*v^2/2 vs. p*c0).
    ST.beta           = sqrt(1e0-1e0/sqr(ST.gamma));
    ST.bg1            = ST.beta*ST.gamma;
    ST.SampleIonK     = 2e0*M_PI/(ST.beta*SampleLambda);

    // Evaluate momentum compaction.
    const double R56 = -2e0*M_PI/(SampleLambda*Erest/MeVtoeV*cube(ST.beta*ST.gamma))*length*MtoMM;

    std::string t_name = type_name(); // C string -> C++ string.
    if (t_name != "rfcavity") {
        ST.FyAbs       += ST.SampleIonK_ref*length*MtoMM;
        ST.Fy_absState += ST.SampleIonK*length*MtoMM;
        ST.Ekinetic     = last_Kenergy_out;

        transfer(state_t::PS_S, state_t::PS_PS) = R56;
    }

    ST.moment0 = prod(transfer, ST.moment0);

    if (t_name == "rfcavity") {
        ST.moment0[state_t::PS_S]  = ST.Fy_absState - ST.FyAbs;
        ST.moment0[state_t::PS_PS] = (ST.EkState-ST.Ekinetic)/MeVtoeV;
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

    transfer(state_t::PS_S, state_t::PS_PS) *= cube(ST.bg_ref/ST.bg1);

    last_Kenergy_in = last_Kenergy_out = ST.Ekinetic; // no energy gain
}

namespace {

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
        length = phase_factor = 0e0;
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
        double L = length*MtoMM; // Convert from [m] to [mm].

        this->transfer_raw(state_t::PS_X, state_t::PS_PX) = L;
        this->transfer_raw(state_t::PS_Y, state_t::PS_PY) = L;
    }
    virtual ~ElementDrift() {}

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
        double L   = c.get<double>("L")*MtoMM,
               phi = c.get<double>("phi"),               // [rad].
               rho = L/phi,
               K   = c.get<double>("K", 0e0)/sqr(MtoMM), // [1/m^2].
               Kx  = K + 1e0/sqr(rho),
               Ky  = -K;

        // Horizontal plane.
        GetQuadMatrix(L, Kx, (unsigned)state_t::PS_X, this->transfer_raw);
        // Vertical plane.
        GetQuadMatrix(L, Ky, (unsigned)state_t::PS_Y, this->transfer_raw);
        // Longitudinal plane.
        // For total path length.
//        this->transfer_raw(state_t::PS_S,  state_t::PS_S) = L;
    }
    virtual ~ElementSBend() {}

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

        double Brho = ST.beta*(ST.EkState+ST.IonEs)/(C0*ST.IonZ),
               K    = conf().get<double>("B2")/Brho/sqr(MtoMM),
               L    = conf().get<double>("L")*MtoMM;

        // Horizontal plane.
        GetQuadMatrix(L,  K, (unsigned)state_t::PS_X, this->transfer_raw);
        // Vertical plane.
        GetQuadMatrix(L, -K, (unsigned)state_t::PS_Y, this->transfer_raw);
        // Longitudinal plane.
//        this->transfer_raw(state_t::PS_S, state_t::PS_S) = L;

        transfer = transfer_raw;

        last_Kenergy_in = last_Kenergy_out = ST.Ekinetic; // no energy gain
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

        double Brho = ST.beta*(ST.EkState+ST.IonEs)/(C0*ST.IonZ),
               K    = conf().get<double>("B")/(2e0*Brho)/MtoMM,
               L    = conf().get<double>("L")*MtoMM;      // Convert from [m] to [mm].

        GetSolMatrix(L, K, this->transfer_raw);

        transfer = transfer_raw;

        last_Kenergy_in = last_Kenergy_out = ST.Ekinetic; // no energy gain
    }

    virtual const char* type_name() const {return "solenoid";}
};


// RF Cavity beam dynamics functions.

typedef boost::numeric::ublas::matrix<double> value_mat;


// Global constants.

// Phase space dimension; including vector for orbit/1st moment.
# define PS_Dim Moment2State::maxsize // Set to 7; to include orbit.

// Mpultipole level: 0 only include focusing and defocusing effects,
//                   1 include dipole terms,
//                   2 include quadrupole terms.
const int MpoleLevel = 2;

typedef typename Moment2ElementBase::state_t state_t;


double PwrSeries(const double beta,
                 const double a0, const double a1, const double a2, const double a3,
                 const double a4, const double a5, const double a6, const double a7,
                 const double a8, const double a9)
{
    int    k;
    double f;

    const int    n   = 10;
    const double a[] = {a0, a1, a2, a3, a4, a5, a6, a7, a8, a9};

    f = a[0];
    for (k = 1; k < n; k++)
        f += a[k]*pow(beta, k);

    return f;
}


void TransFacts(const int cavilabel, double beta, const int gaplabel, const double EfieldScl,
                double &Ecen, double &T, double &Tp, double &S, double &Sp, double &V0)
{
    // Evaluate Electric field center, transit factors [T, T', S, S'] and cavity field.
    std::vector<double> vec;

    switch (cavilabel) {
    case 41:
        if (beta < 0.025 || beta > 0.08) {
            std::cerr << "*** GetTransitFac: beta out of Range " << beta << "\n";
            exit(1);
        }
        switch (gaplabel) {
        case 0:
            // One gap evaluation.
            Ecen = 120.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = PwrSeries(beta, -4.109, 399.9, -1.269e4, 1.991e5, -1.569e6, 4.957e6, 0.0, 0.0, 0.0, 0.0);
            Sp   = PwrSeries(beta, 61.98, -1.073e4, 4.841e5, 9.284e6, 8.379e7, -2.926e8, 0.0, 0.0, 0.0, 0.0);
            V0   = 0.98477*EfieldScl;
            break;
        case 1:
            // Two gap calculation, first gap.
            Ecen = 0.0006384*pow(beta, -1.884) + 86.69;
            T    = PwrSeries(beta, 0.9232, -123.2, 3570, -5.476e4, 4.316e5, -1.377e6, 0.0, 0.0, 0.0, 0.0);
            Tp   = PwrSeries(beta, 1.699, 924.7, -4.062e4, 7.528e5, -6.631e6, 2.277e7, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp   = PwrSeries(beta, -1.571, 25.59, 806.6, -2.98e4, 3.385e5, -1.335e6, 0.0, 0.0, 0.0, 0.0);
            V0   = 0.492385*EfieldScl;
            break;
        case 2:
            // Two gap calculation, second gap.
            Ecen = -0.0006384*pow(beta, -1.884) + 33.31;
            T    = PwrSeries(beta, -0.9232, 123.2, -3570, 5.476e4, -4.316e5, 1.377e6, 0.0, 0.0, 0.0, 0.0);
            Tp   = PwrSeries(beta, -1.699, -924.7, 4.062e4, -7.528e5, 6.631e6, -2.277e7, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp    = PwrSeries(beta, -1.571, 25.59, 806.6, -2.98e4, 3.385e5, -1.335e6, 0.0, 0.0, 0.0, 0.0);
            V0   = 0.492385*EfieldScl;
            break;
        default:
            std::cerr << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            exit(1);
        }
        break;
    case 85:
        if (beta < 0.05 || beta > 0.25) {
            std::cerr << "*** GetTransitFac: beta out of range " << beta << "\n";
            exit(1);
        }
        switch (gaplabel) {
          case 0:
            Ecen = 150.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = PwrSeries(beta, -6.811, 343.9, -6385, 6.477e4, -3.914e5, 1.407e6, -2.781e6, 2.326e6, 0.0, 0.0);
            Sp   = PwrSeries(beta, 162.7, -1.631e4, 4.315e5, -5.344e6, 3.691e7, -1.462e8, 3.109e8, -2.755e8, 0.0, 0.0);
            V0   = 1.967715*EfieldScl;
            break;
        case 1:
            Ecen = 0.0002838*pow(beta, -2.13) + 76.5;
            T    = 0.0009467*pow(beta, -1.855) - 1.002;
            Tp   = PwrSeries(beta, 24.44, -334, 2468, -1.017e4, 2.195e4, -1.928e4, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp   = -0.0009751*pow(beta, -1.898) + 0.001568;
            V0   = 0.9838574*EfieldScl;
            break;
        case 2:
            Ecen = -0.0002838*pow(beta, -2.13) + 73.5;
            T    = -0.0009467*pow(beta, -1.855) + 1.002;
            Tp   = PwrSeries(beta,  24.44, 334,  2468, 1.017e4, -2.195e4, 1.928e4, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp   = -0.0009751*pow(beta, -1.898) + 0.001568;
            V0   = 0.9838574*EfieldScl;
            break;
        default:
            std::cerr << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            exit(1);
        }
        break;
    default:
        std::cerr << "*** GetTransitFac: undef. cavity type" << "\n";
        exit(1);
    }

    // Convert from [mm] to [m].
//    Ecen /= MtoMM;
}


void TransitFacMultipole(const int cavi, const std::string &flabel, const double IonK,
                         double &T, double &S)
{

    if ((cavi == 1) && (IonK < 0.025 || IonK > 0.055)) {
        std::cerr << "*** TransitFacMultipole: IonK out of Range" << "\n";
        exit(1);
    } else if ((cavi == 2) && (IonK < 0.006 || IonK > 0.035)) {
        std::cerr << "*** TransitFacMultipole: IonK out of Range" << "\n";
    }

    if (flabel == "CaviMlp_EFocus1") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 1.256386e+02, -3.108322e+04, 3.354464e+06, -2.089452e+08, 8.280687e+09, -2.165867e+11,
                          3.739846e+12, -4.112154e+13, 2.613462e14, -7.316972e14);
            S = PwrSeries(IonK, 1.394183e+02, -3.299673e+04, 3.438044e+06, -2.070369e+08, 7.942886e+09, -2.013750e+11,
                         3.374738e+12, -3.605780e+13, 2.229446e+14, -6.079177e+14);
            break;
        case 2:
            T = PwrSeries(IonK, -9.450041e-01, -3.641390e+01, 9.926186e+03, -1.449193e+06, 1.281752e+08, -7.150297e+09,
                          2.534164e+11, -5.535252e+12, 6.794778e+13, -3.586197e+14);
            S = PwrSeries(IonK, 9.928055e-02, -5.545119e+01, 1.280168e+04, -1.636888e+06, 1.279801e+08, -6.379800e+09,
                          2.036575e+11, -4.029152e+12, 4.496323e+13, -2.161712e+14);
            break;
        }
    } else if (flabel == "CaviMlp_EFocus2") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 1.038803e+00, -9.121320e+00, 8.943931e+02, -5.619149e+04, 2.132552e+06, -5.330725e+07,
                          8.799404e+08, -9.246033e+09, 5.612073e+10, -1.499544e+11);
            S = PwrSeries(IonK, 1.305154e-02, -2.585211e+00, 2.696971e+02, -1.488249e+04, 5.095765e+05, -1.154148e+07,
                          1.714580e+08, -1.604935e+09, 8.570757e+09, -1.983302e+10);
            break;
        case 2:
            T = PwrSeries(IonK, 9.989307e-01, 7.299233e-01, -2.932580e+02, 3.052166e+04, -2.753614e+06, 1.570331e+08,
                          -5.677804e+09, 1.265012e+11, -1.584238e+12, 8.533351e+12);
            S = PwrSeries(IonK, -3.040839e-03, 2.016667e+00, -4.313590e+02, 5.855139e+04, -4.873584e+06, 2.605444e+08,
                          -8.968899e+09, 1.923697e+11, -2.339920e+12, 1.233014e+13);
            break;
        }
    } else if (flabel == "CaviMlp_EDipole") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, -1.005885e+00, 1.526489e+00, -1.047651e+02, 1.125013e+04, -4.669147e+05, 1.255841e+07,
                          -2.237287e+08, 2.535541e+09, -1.656906e+10, 4.758398e+10);
            S = PwrSeries(IonK, -2.586200e-02, 5.884367e+00, -6.407538e+02, 3.888964e+04, -1.488484e+06, 3.782592e+07,
                          -6.361033e+08, 6.817810e+09, -4.227114e+10, 1.155597e+11);
            break;
        case 2:
            T = PwrSeries(IonK, -9.999028e-01, -6.783669e-02, 1.415756e+02, -2.950990e+03, 2.640980e+05, -1.570742e+07,
                          5.770450e+08, -1.303686e+10, 1.654958e+11, -9.030017e+11);
            S = PwrSeries(IonK, 2.108581e-04, -3.700608e-01, 2.851611e+01, -3.502994e+03, 2.983061e+05, -1.522679e+07,
                          4.958029e+08, -1.002040e+10, 1.142835e+11, -5.617061e+11);
            break;
        }
    } else if (flabel == "CaviMlp_EQuad") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 1.038941e+00, -9.238897e+00, 9.127945e+02, -5.779110e+04, 2.206120e+06, -5.544764e+07,
                          9.192347e+08, -9.691159e+09, 5.896915e+10, -1.578312e+11);
            S = PwrSeries(IonK, 1.248096e-01, -2.923507e+01, 3.069331e+03, -1.848380e+05, 7.094882e+06, -1.801113e+08,
                          3.024208e+09, -3.239241e+10, 2.008767e+11, -5.496217e+11);
            break;
        case 2:
            T = PwrSeries(IonK, 1.000003e+00, -1.015639e-03, -1.215634e+02, 1.720764e+01, 3.921401e+03, 2.674841e+05,
                          -1.236263e+07, 3.128128e+08, -4.385795e+09, 2.594631e+10);
            S = PwrSeries(IonK, -1.756250e-05, 2.603597e-01, -2.551122e+00, -4.840638e+01, -2.870201e+04, 1.552398e+06,
                          -5.135200e+07, 1.075958e+09, -1.277425e+10, 6.540748e+10);
            break;
        }
    } else if (flabel == "CaviMlp_HMono") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 1.703336e+00, -1.671357e+02, 1.697657e+04, -9.843253e+05, 3.518178e+07, -8.043084e+08,
                          1.165760e+10, -1.014721e+11, 4.632851e+11, -7.604796e+11);
            S = PwrSeries(IonK, 1.452657e+01, -3.409550e+03, 3.524921e+05, -2.106663e+07, 8.022856e+08,
                          -2.019481e+10, 3.360597e+11, -3.565836e+12, 2.189668e+13, -5.930241e+13);
            break;
        case 2:
            T = PwrSeries(IonK, 1.003228e+00, -1.783406e+00, 1.765330e+02, -5.326467e+04, 4.242623e+06, -2.139672e+08,
                          6.970488e+09, -1.411958e+11, 1.617248e+12, -8.000662e+12);
            S = PwrSeries(IonK, -1.581533e-03, 1.277444e+00, -2.742508e+02, 3.966879e+04, -3.513478e+06, 1.962939e+08,
                          -6.991916e+09, 1.539708e+11, -1.910236e+12, 1.021016e+13);
            break;
        }
    } else if (flabel == "CaviMlp_HDipole") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 6.853803e-01, 7.075414e+01, -7.117391e+03, 3.985674e+05, -1.442888e+07, 3.446369e+08,
                          -5.420826e+09, 5.414689e+10, -3.116216e+11, 7.869717e+11);
            S = PwrSeries(IonK, 1.021102e+00, -2.441117e+02, 2.575274e+04, -1.569273e+06, 6.090118e+07, -1.562284e+09,
                          2.649289e+10, -2.864139e+11, 1.791634e+12, -4.941947e+12);
            break;
        case 2:
            T = PwrSeries(IonK, 1.014129e+00, -8.016304e+00, 1.631339e+03, -2.561826e+05, 2.115355e+07, -1.118723e+09,
                          3.821029e+10, -8.140248e+11, 9.839613e+12, -5.154137e+13);
            S = PwrSeries(IonK, -4.688714e-03, 3.299051e+00, -8.101936e+02, 1.163814e+05, -1.017331e+07, 5.607330e+08,
                          -1.967300e+10, 4.261388e+11, -5.194592e+12, 2.725370e+13);
            break;
        }
    } else if (flabel == "CaviMlp_HQuad") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, -1.997432e+00, 2.439177e+02, -2.613724e+04, 1.627837e+06, -6.429625e+07, 1.676173e+09,
                          -2.885455e+10, 3.163675e+11, -2.005326e+12, 5.600545e+12);
            S = PwrSeries(IonK, -2.470704e+00, 5.862902e+02, -6.135071e+04, 3.711527e+06, -1.431267e+08, 3.649414e+09,
                          -6.153570e+10, 6.617859e+11, -4.119861e+12, 1.131390e+13);
            break;
        case 2:
            T = PwrSeries(IonK, -1.000925e+00, 5.170302e-01, 9.311761e+01, 1.591517e+04, -1.302247e+06, 6.647808e+07,
                          -2.215417e+09, 4.603390e+10, -5.420873e+11, 2.764042e+12);
            S = PwrSeries(IonK, 3.119419e-04, -4.540868e-01, 5.433028e+01, -7.571946e+03, 6.792565e+05, -3.728390e+07,
                          1.299263e+09, -2.793705e+10, 3.377097e+11, -1.755126e+12);
            break;
        }
    } else {
        std::cerr << "*** TransitFacMultipole: undef. multipole type " << flabel << "\n";
        exit(1);
    }
}


double GetCavPhase(const int cavi, state_t &ST, const double IonFys, const double multip)
{
    /* If the cavity is not at full power, the method gives synchrotron
     * phase slightly different from the nominal value.                 */

    double IonEk, Fyc;

    IonEk = (ST.IonW_ref-ST.IonEs)/MeVtoeV;

    switch (cavi) {
    case 1:
        Fyc = 4.394*pow(IonEk, -0.4965) - 4.731;
        break;
    case 2:
        Fyc = 5.428*pow(IonEk, -0.5008) + 1.6;
        break;
    case 3:
        Fyc = 22.35*pow(IonEk, -0.5348) + 2.026;
        break;
    case 4:
        Fyc = 41.43*pow(IonEk, -0.5775) + 2.59839;
        break;
    case 5:
        Fyc = 5.428*pow(IonEk, -0.5008) + 1.6;
        break;
    default:
        std::cerr << "*** GetCavPhase: undef. cavity type" << "\n";
        exit(1);
    }

    return IonFys - Fyc - ST.FyAbs*multip;
}


void GetCavMatParams(const int cavi, std::stringstream &inf,
                     const double beta_tab[], const double gamma_tab[], const double IonK[])
{
    // Evaluate time transit factors and acceleration.

    std::string       line, Elem, Name;
    double            s, Length, Aper, E0, T, S, Accel;
    std::stringstream str;

    inf.clear();
    inf.seekg(0, inf.beg);

    CavTLMLineTab2[cavi-1].clear();

    s = CavData2[cavi-1].s[0];
    while (getline(inf, line) && !inf.fail()) {
        T = 0e0, S = 0e0, Accel = 0e0;
        if (line[0] == '%') {
            // Comment.
        } else {
            str.str(line);
            str >> Elem >> Name >> Length >> Aper;

            s += Length;

            if ((Elem != "drift") && (Elem != "AccGap"))
                str >> E0;
            else
                E0 = 0e0;

            if (Elem == "drift") {
            } else if (Elem == "EFocus1") {
                if (s < 0e0) {
                    // First gap. By reflection 1st Gap EFocus1 is 2nd gap EFocus2.
                    TransitFacMultipole(cavi, "CaviMlp_EFocus2", IonK[0], T, S);
                    // First gap *1, transverse E field the same.
                    S = -S;
                } else {
                    // Second gap.
                    TransitFacMultipole(cavi, "CaviMlp_EFocus1", IonK[1], T, S);
                }
            } else if (Elem == "EFocus2") {
                if (s < 0e0) {
                    // First gap.
                    TransitFacMultipole(cavi, "CaviMlp_EFocus1", IonK[0], T, S);
                    S = -S;
                } else {
                    // Second gap.
                    TransitFacMultipole(cavi, "CaviMlp_EFocus2", IonK[1], T, S);
                }
            } else if (Elem == "EDipole") {
                if (MpoleLevel >= 1) {
                    if (s < 0e0) {
                        TransitFacMultipole(cavi, "CaviMlp_EDipole", IonK[0], T, S);
                        // First gap *1, transverse E field the same.
                        S = -S;
                    } else {
                        // Second gap.
                        TransitFacMultipole(cavi, "CaviMlp_EDipole", IonK[1], T, S);
                    }
                }
            } else if (Elem == "EQuad") {
                if (MpoleLevel >= 2) {
                    if (s < 0e0) {
                        // First gap.
                        TransitFacMultipole(cavi, "CaviMlp_EQuad", IonK[0], T, S);
                        S = -S;
                    } else {
                        // Second Gap
                        TransitFacMultipole(cavi, "CaviMlp_EQuad", IonK[1], T, S);
                    }
                }
            } else if (Elem == "HMono") {
                if (MpoleLevel >= 2) {
                    if (s < 0e0) {
                        // First gap.
                        TransitFacMultipole(cavi, "CaviMlp_HMono", IonK[0], T, S);
                        T = -T;
                    } else {
                        // Second Gap
                        TransitFacMultipole(cavi, "CaviMlp_HMono", IonK[1], T, S);
                    }
                }
            } else if (Elem == "HDipole") {
                if (MpoleLevel >= 1) {
                    if (s < 0e0) {
                        // First gap.
                        TransitFacMultipole(cavi, "CaviMlp_HDipole", IonK[0], T, S);
                        T = -T;
                    }  else {
                        // Second gap.
                        TransitFacMultipole(cavi, "CaviMlp_HDipole", IonK[1], T, S);
                    }
                }
            } else if (Elem == "HQuad") {
                if (MpoleLevel >= 2) {
                    if (s < 0e0) {
                        // First gap.
                        TransitFacMultipole(cavi, "CaviMlp_HQuad", IonK[0], T, S);
                        T = -T;
                    } else {
                        // Second gap.
                        TransitFacMultipole(cavi, "CaviMlp_HQuad", IonK[1], T, S);
                    }
                }
            } else if (Elem == "AccGap") {
                if (s < 0e0) {
                    // First gap.
                    Accel = (beta_tab[0]*gamma_tab[0])/((beta_tab[1]*gamma_tab[1]));
                } else {
                    // Second gap.
                    Accel = (beta_tab[1]*gamma_tab[1])/((beta_tab[2]*gamma_tab[2]));
                }
            } else {
                std::cerr << "*** GetCavMatParams: undef. multipole element " << Elem << "\n";
                exit(1);
            }

            CavTLMLineTab2[cavi-1].set(s, Elem, E0, T, S, Accel);
        }
    }

    if (false) {
        std::cout << "\n";
        CavTLMLineTab2[cavi-1].show(std::cout);
    }
}


void EvalGapModel(const double dis, const double IonW0, state_t &ST, const double IonFy0,
                  const double k, const double Lambda, const double Ecen,
                  const double T, const double S, const double Tp, const double Sp, const double V0,
                  double &IonW_f, double &IonFy_f)
{
    double Iongamma_f, IonBeta_f, k_f;

    IonW_f     = IonW0 + ST.IonZ*V0*T*cos(IonFy0+k*Ecen)*MeVtoeV - ST.IonZ*V0*S*sin(IonFy0+k*Ecen)*MeVtoeV;
    Iongamma_f = IonW_f/ST.IonEs;
    IonBeta_f  = sqrt(1e0-1e0/sqr(Iongamma_f));
    k_f        = 2e0*M_PI/(IonBeta_f*Lambda);

    IonFy_f = IonFy0 + k*Ecen + k_f*(dis-Ecen)
              + ST.IonZ*V0*k*(Tp*sin(IonFy0+k*Ecen)+Sp*cos(IonFy0+k*Ecen))/(2e0*(IonW0-ST.IonEs)/MeVtoeV);
}


void GenCavMat(const int cavi, const double dis, const double EfieldScl, const double TTF_tab[],
               const double beta_tab[], const double gamma_tab[], const double Lambda,
               state_t &ST, const double IonFys[], std::stringstream &inf,
               const double Rm, value_mat &M)
{
    /* RF cavity model, transverse only defocusing.
     * 2-gap matrix model.                                            */

    std::string       line, Elem, Name;
    int               seg, n;
    double            Length, Aper, Efield, s, k_s[3];
    double            Ecens[2], Ts[2], Ss[2], V0s[2], ks[2], L1, L2, L3;
    double            beta, gamma, kfac, V0, T, S, kfdx, kfdy, dpy, Accel, IonFy;
    value_mat         Idmat, Mlon_L1, Mlon_K1, Mlon_L2;
    value_mat         Mlon_K2, Mlon_L3, Mlon, Mtrans, Mprob;
    std::stringstream str;

    const double IonA = 1e0;

    using boost::numeric::ublas::prod;

    Idmat = boost::numeric::ublas::identity_matrix<double>(PS_Dim);

    inf.clear();
    inf.seekg(0, inf.beg);

    k_s[0] = 2e0*M_PI/(beta_tab[0]*Lambda);
    k_s[1] = 2e0*M_PI/(beta_tab[1]*Lambda);
    k_s[2] = 2e0*M_PI/(beta_tab[2]*Lambda);

    // Longitudinal model:  Drift-Kick-Drift, dis: total lenghth centered at 0,
    // Ecens[0] & Ecens[1]: Electric Center position where accel kick applies, Ecens[0] < 0
    // TTFtab:              2*6 vector, Ecens, T Tp S Sp, V0;

    Ecens[0] = TTF_tab[0];
    Ts[0]    = TTF_tab[1];
    Ss[0]    = TTF_tab[3];
    V0s[0]   = TTF_tab[5];
    ks[0]    = 0.5*(k_s[0]+k_s[1]);
    L1       = dis + Ecens[0];       //try change dis/2 to dis 14/12/12

    Mlon_L1 = Idmat;
    Mlon_K1 = Idmat;
    // Pay attention, original is -
    Mlon_L1(4, 5) = -2e0*M_PI/Lambda*(1e0/cube(beta_tab[0]*gamma_tab[0])*MeVtoeV/ST.IonEs*L1);
    // Pay attention, original is -k1-k2
    Mlon_K1(5, 4) = -ST.IonZ*V0s[0]*Ts[0]*sin(IonFys[0]+ks[0]*L1)-ST.IonZ*V0s[0]*Ss[0]*cos(IonFys[0]+ks[0]*L1);

    Ecens[1] = TTF_tab[6];
    Ts[1]    = TTF_tab[7];
    Ss[1]    = TTF_tab[9];
    V0s[1]   = TTF_tab[11];
    ks[1]    = 0.5*(k_s[1]+k_s[2]);
    L2       = Ecens[1] - Ecens[0];

    Mlon_L2 = Idmat;
    Mlon_K2 = Idmat;

    Mlon_L2(4, 5) = -2e0*M_PI/Lambda*(1e0/cube(beta_tab[1]*gamma_tab[1])*MeVtoeV/ST.IonEs*L2); //Problem is Here!!
    Mlon_K2(5, 4) = -ST.IonZ*V0s[1]*Ts[1]*sin(IonFys[1]+ks[1]*Ecens[1])-ST.IonZ*V0s[1]*Ss[1]*cos(IonFys[1]+ks[1]*Ecens[1]);

    L3 = dis - Ecens[1]; //try change dis/2 to dis 14/12/12

    Mlon_L3       = Idmat;
    Mlon_L3(4, 5) = -2e0*M_PI/Lambda*(1e0/cube(beta_tab[2]*gamma_tab[2])*MeVtoeV/ST.IonEs*L3);

    Mlon = Idmat;
    Mlon = prod(Mlon_K1, Mlon_L1);
    Mlon = prod(Mlon_L2, Mlon);
    Mlon = prod(Mlon_K2, Mlon);
    Mlon = prod(Mlon_L3, Mlon);

    // Transverse model
    // Drift-FD-Drift-LongiKick-Drift-FD-Drift-0-Drift-FD-Drift-LongiKick-Drift-FD-Drift

    seg    = 0;

    Mtrans = Idmat;
    Mprob  = Idmat;

    beta   = beta_tab[0];
    gamma  = gamma_tab[0];
    IonFy  = IonFys[0];
    kfac   = k_s[0];

    V0 = 0e0, T = 0e0, S = 0e0, kfdx = 0e0, kfdy = 0e0, dpy = 0e0;

    s = CavData2[cavi-1].s[0];
    n = 0;
    while (getline(inf, line) && !inf.fail()) {
        if (line[0] == '%') {
            // Comment.
        } else {
            n++;
            str.str(line);
            str >> Elem >> Name >> Length >> Aper;

            s += Length;

            if ((Elem != "drift") && (Elem != "AccGap"))
                str >> Efield;
            else
                Efield = 0e0;

            if (false)
                printf("%9.5f %8s %8s %9.5f %9.5f %9.5f\n",
                       s, Elem.c_str(), Name.c_str(), Length, Aper, Efield);

            Mprob = Idmat;
            if (Elem == "drift") {
                IonFy = IonFy + kfac*Length;

                Mprob(0, 1) = Length;
                Mprob(2, 3) = Length;
                Mtrans      = prod(Mprob, Mtrans);
            } else if (Elem == "EFocus1") {
                V0   = CavTLMLineTab2[cavi-1].E0[n-1]*EfieldScl;
                T    = CavTLMLineTab2[cavi-1].T[n-1];
                S    = CavTLMLineTab2[cavi-1].S[n-1];
                kfdx = ST.IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy))/Rm;
                kfdy = kfdx;

                Mprob(1, 0) = kfdx;
                Mprob(3, 2) = kfdy;
                Mtrans      = prod(Mprob, Mtrans);
            } else if (Elem == "EFocus2") {
                V0   = CavTLMLineTab2[cavi-1].E0[n-1]*EfieldScl;
                T    = CavTLMLineTab2[cavi-1].T[n-1];
                S    = CavTLMLineTab2[cavi-1].S[n-1];
                kfdx = ST.IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy))/Rm;
                kfdy = kfdx;

                Mprob(1, 0) = kfdx;
                Mprob(3, 2) = kfdy;
                Mtrans      = prod(Mprob, Mtrans);
            } else if (Elem == "EDipole") {
                if (MpoleLevel >= 1) {
                    V0  = CavTLMLineTab2[cavi-1].E0[n-1]*EfieldScl;
                    T   = CavTLMLineTab2[cavi-1].T[n-1];
                    S   = CavTLMLineTab2[cavi-1].S[n-1];
                    dpy = ST.IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy));

                    Mprob(3, 6) = dpy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "EQuad") {
                if (MpoleLevel >= 2) {
                    V0   = CavTLMLineTab2[cavi-1].E0[n-1]*EfieldScl;
                    T    = CavTLMLineTab2[cavi-1].T[n-1];
                    S    = CavTLMLineTab2[cavi-1].S[n-1];
                    kfdx =  ST.IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy))/Rm;
                    kfdy = -kfdx;

                    Mprob(1, 0) = kfdx;
                    Mprob(3, 2) = kfdy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "HMono") {
                if (MpoleLevel >= 2) {
                    V0   = CavTLMLineTab2[cavi-1].E0[n-1]*EfieldScl;
                    T    = CavTLMLineTab2[cavi-1].T[n-1];
                    S    = CavTLMLineTab2[cavi-1].S[n-1];
                    kfdx = -MU0*C0*ST.IonZ*V0/beta/gamma/IonA/AU*(T*cos(IonFy+M_PI/2e0)-S*sin(IonFy+M_PI/2e0))/Rm;
                    kfdy = kfdx;

                    Mprob(1, 0) = kfdx;
                    Mprob(3, 2) = kfdy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "HDipole") {
                if (MpoleLevel >= 1) {
                    V0  = CavTLMLineTab2[cavi-1].E0[n-1]*EfieldScl;
                    T   = CavTLMLineTab2[cavi-1].T[n-1];
                    S   = CavTLMLineTab2[cavi-1].S[n-1];
                    dpy = -MU0*C0*ST.IonZ*V0/beta/gamma/IonA/AU*(T*cos(IonFy+M_PI/2e0)-S*sin(IonFy+M_PI/2e0));

                    Mprob(3, 6) = dpy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "HQuad") {
                if (MpoleLevel >= 2) {
                    if (s < 0e0) {
                        // First gap.
                        beta  = (beta_tab[0]+beta_tab[1])/2e0;
                        gamma = (gamma_tab[0]+gamma_tab[1])/2e0;
                    } else {
                        beta  = (beta_tab[1]+beta_tab[2])/2e0;
                        gamma = (gamma_tab[1]+gamma_tab[2])/2e0;
                    }
                    V0   = CavTLMLineTab2[cavi-1].E0[n-1]*EfieldScl;
                    T    = CavTLMLineTab2[cavi-1].T[n-1];
                    S    = CavTLMLineTab2[cavi-1].S[n-1];
                    kfdx = -MU0*C0*ST.IonZ*V0/beta/gamma/IonA/AU*(T*cos(IonFy+M_PI/2e0)-S*sin(IonFy+M_PI/2e0))/Rm;
                    kfdy = -kfdx;

                    Mprob(1, 0) = kfdx;
                    Mprob(3, 2) = kfdy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "AccGap") {
                //IonFy = IonFy + ST.IonZ*V0s[0]*kfac*(TTF_tab[2]*sin(IonFy)
                //        + TTF_tab[4]*cos(IonFy))/2/((gamma-1)*ST.IonEs/MeVtoeV); //TTF_tab[2]~Tp
                seg    = seg + 1;
                beta   = beta_tab[seg];
                gamma  = gamma_tab[seg];
                kfac   = 2e0*M_PI/(beta*Lambda);
                Accel  = CavTLMLineTab2[cavi-1].Accel[n-1];

                Mprob(1, 1) = Accel;
                Mprob(3, 3) = Accel;
                Mtrans      = prod(Mprob, Mtrans);
            } else {
                std::cerr << "*** GenCavMat: undef. multipole type " << Elem << "\n";
                exit(1);
            }
//            std::cout << Elem << "\n";
//            PrtMat(Mprob);
        }
    }

//    inf.close();

    M = Mtrans;

    M(4, 4) = Mlon(4, 4);
    M(4, 5) = Mlon(4, 5);
    M(5, 4) = Mlon(5, 4);
    M(5, 5) = Mlon(5, 5);
}


void GetCavMat(const int cavi, const int cavilabel, const double Rm, state_t &ST,
               const double EfieldScl, const double IonFyi_s,
               const double IonEk_s, const double fRF, value_mat &M)
{
    int    n;
    double IonLambda, Ecen[2], T[2], Tp[2], S[2], Sp[2], V0[2];
    double dis, IonW_s[3], IonFy_s[3], gamma_s[3], beta_s[3], IonK_s[3];
    double IonK[2];

    IonLambda  = C0/fRF*MtoMM;

    IonW_s[0]  = IonEk_s + ST.IonEs;
    IonFy_s[0] = IonFyi_s;
    gamma_s[0] = IonW_s[0]/ST.IonEs;
    beta_s[0]  = sqrt(1e0-1e0/sqr(gamma_s[0]));
    IonK_s[0]  = 2e0*M_PI/(beta_s[0]*IonLambda);

    n   = CavData2[cavi-1].s.size();
    dis = (CavData2[cavi-1].s[n-1]-CavData2[cavi-1].s[0])/2e0;

    TransFacts(cavilabel, beta_s[0], 1, EfieldScl, Ecen[0], T[0], Tp[0], S[0], Sp[0], V0[0]);
    EvalGapModel(dis, IonW_s[0], ST, IonFy_s[0], IonK_s[0], IonLambda,
                Ecen[0], T[0], S[0], Tp[0], Sp[0], V0[0], IonW_s[1], IonFy_s[1]);
    gamma_s[1] = IonW_s[1]/ST.IonEs;
    beta_s[1]  = sqrt(1e0-1e0/sqr(gamma_s[1]));
    IonK_s[1]  = 2e0*M_PI/(beta_s[1]*IonLambda);

    TransFacts(cavilabel, beta_s[1], 2, EfieldScl, Ecen[1], T[1], Tp[1], S[1], Sp[1], V0[1]);
    EvalGapModel(dis, IonW_s[1], ST, IonFy_s[1], IonK_s[1], IonLambda,
                Ecen[1], T[1], S[1], Tp[1], Sp[1], V0[1], IonW_s[2], IonFy_s[2]);
    gamma_s[2] = IonW_s[2]/ST.IonEs;
    beta_s[2]  = sqrt(1e0-1e0/sqr(gamma_s[2]));
    IonK_s[2]  = 2e0*M_PI/(beta_s[2]*IonLambda);

    Ecen[0] = Ecen[0] - dis;

    double TTF_tab[] = {Ecen[0], T[0], Tp[0], S[0], Sp[0], V0[0], Ecen[1], T[1], Tp[1], S[1], Sp[1], V0[1]};
    IonK[0] = (IonK_s[0]+IonK_s[1])/2e0;
    IonK[1] = (IonK_s[1]+IonK_s[2])/2e0;

    if (false) {
        printf("\n GetCavMat:\n");
        printf("IonK  : %15.10f %15.10f %15.10f\n", IonK_s[0], IonK_s[1], IonK_s[2]);
        printf("IonK  : %15.10f %15.10f\n", IonK[0], IonK[1]);
        printf("beta  : %15.10f %15.10f %15.10f\n", beta_s[0], beta_s[1], beta_s[2]);
        printf("gamma : %15.10f %15.10f %15.10f\n", gamma_s[0], gamma_s[1], gamma_s[2]);
        printf("Ecen  : %15.10f %15.10f\n", Ecen[0], Ecen[1]);
        printf("T     : %15.10f %15.10f\n", T[0], T[1]);
        printf("Tp    : %15.10f %15.10f\n", Tp[0], Tp[1]);
        printf("S     : %15.10f %15.10f\n", S[0], S[1]);
        printf("Sp    : %15.10f %15.10f\n", Sp[0], Sp[1]);
        printf("V0    : %15.10f %15.10f\n", V0[0], V0[1]);
    }

    GetCavMatParams(cavi, CavTLMstream2[cavi-1], beta_s, gamma_s, IonK);
    GenCavMat(cavi, dis, EfieldScl, TTF_tab, beta_s, gamma_s, IonLambda, ST, IonFy_s, CavTLMstream2[cavi-1], Rm, M);
}


void GetCavBoost(const CavDataType &CavData, const double IonW0,
                 const double IonFy0, const double IonK0, const double IonZ,
                 const double IonEs, const double fRF,
                 const double EfieldScl, double &IonW, double &IonFy, double &accIonW)
{
    int    n = CavData.s.size(),
           k;

    double dis = CavData.s[n-1] - CavData.s[0],
           dz  = dis/(n-1),
           IonLambda, IonK, IonFylast, IonGamma, IonBeta;

    IonLambda = C0/fRF*MtoMM;


    IonFy = IonFy0;
    IonK  = IonK0;
    IonW  = IonW0;
    for (k = 0; k < n-1; k++) {
        IonFylast = IonFy;
        IonFy += IonK*dz;
        IonW  += IonZ*EfieldScl*(CavData.Elong[k]+CavData.Elong[k+1])/2e0
                 *cos((IonFylast+IonFy)/2e0)*dz/MtoMM;
        IonGamma = IonW/IonEs;
        IonBeta = sqrt(1e0-1e0/sqr(IonGamma));
        if ((IonW-IonEs) < 0e0) {
            IonW = IonEs;
            IonBeta = 0e0;
        }
        IonK = 2e0*M_PI/(IonBeta*IonLambda);
    }

    accIonW = IonW - IonW0;
}


void PropagateLongRFCav(Config &conf, state_t &ST)
{
    std::string CavType;
    int         cavi;
    double      fRF, multip, caviIonK, IonFys, EfieldScl, caviFy, IonFy_i, IonFy_o;
    double      IonW_o, accIonW;

    CavType = conf.get<std::string>("cavtype");
    if (CavType == "0.041QWR") {
        cavi = 1;
    } else if (conf.get<std::string>("cavtype") == "0.085QWR") {
        cavi = 2;
    } else {
        std::cerr << "*** PropagateLongRFCav: undef. cavity type: " << CavType << "\n";
        exit(1);
    }

    fRF       = conf.get<double>("f");
    multip    = fRF/SampleFreq;
    caviIonK  = 2e0*M_PI*fRF/(ST.beta_ref*C0)/MtoMM;
    IonFys    = conf.get<double>("phi")*M_PI/180e0;  // Synchrotron phase [rad].
    EfieldScl = conf.get<double>("scl_fac");         // Electric field scale factor.

    caviFy = GetCavPhase(cavi, ST, IonFys, multip);

    IonFy_i = multip*ST.FyAbs + caviFy;
    conf.set<double>("phi_ref", caviFy);

    // For the reference particle, evaluate the change of:
    // kinetic energy, absolute phase, beta, and gamma.
    GetCavBoost(CavData2[cavi-1], ST.IonW_ref, IonFy_i, caviIonK, ST.IonZ,
                ST.IonEs, fRF, EfieldScl, IonW_o, IonFy_o, accIonW);

    ST.IonW_ref        = IonW_o;
    ST.gamma_ref       = ST.IonW_ref/ST.IonEs;
    ST.beta_ref        = sqrt(1e0-1e0/sqr(ST.gamma_ref));
    ST.SampleIonK_ref  = 2e0*M_PI/(ST.beta_ref*SampleLambda);
    ST.FyAbs          += (IonFy_o-IonFy_i)/multip;
    ST.Ekinetic        = ST.IonW_ref - ST.IonEs;
    ST.gamma           = (ST.EkState+ST.IonEs)/ST.IonEs;
    ST.beta            = sqrt(1e0-1e0/sqr(ST.gamma));
}


void InitRFCav(const Config &conf, state_t &ST, double &accIonW,
               double &avebeta, double &avegamma, value_mat &M)
{
    std::string CavType;
    int         cavi, cavilabel, multip;
    double      Rm, IonFy_i, Ek_i, IonW, fRF, CaviIonK, EfieldScl;
    double      IonW_o, IonFy_o;

    CavType = conf.get<std::string>("cavtype");
    if (CavType == "0.041QWR") {
        cavi       = 1;
        cavilabel  = 41;
        multip     = 1;
        Rm         = 17e0;
    } else if (conf.get<std::string>("cavtype") == "0.085QWR") {
        cavi       = 2;
        cavilabel  = 85;
        multip     = 1;
        Rm         = 17e0;
    } else {
        std::cerr << "*** InitRFCav: undef. cavity type: " << CavType << "\n";
        exit(1);
    }

    IonFy_i = multip*ST.Fy_absState + conf.get<double>("phi_ref");
    Ek_i    = ST.EkState;
    IonW    = ST.EkState + ST.IonEs;

    avebeta   = ST.beta;
    avegamma  = ST.gamma;
    fRF       = conf.get<double>("f");
    CaviIonK  = 2e0*M_PI*fRF/(ST.beta*C0*MtoMM);
    //double SampleIonK = 2e0*M_PI/(ST.beta*C0/SampleFreq*MtoMM);
    EfieldScl  = conf.get<double>("scl_fac");         // Electric field scale factor.

    GetCavBoost(CavData2[cavi-1], IonW, IonFy_i, CaviIonK, ST.IonZ, ST.IonEs,
                fRF, EfieldScl, IonW_o, IonFy_o, accIonW);

    ST.EkState     = IonW_o - ST.IonEs;
    IonW           = ST.EkState + ST.IonEs;
    ST.gamma       = IonW/ST.IonEs;
    ST.beta        = sqrt(1e0-1e0/sqr(ST.gamma));
    ST.Fy_absState += (IonFy_o-IonFy_i)/multip;
    avebeta        += ST.beta;
    avebeta        /= 2e0;
    avegamma       += ST.gamma;
    avegamma       /= 2e0;

    GetCavMat(cavi, cavilabel, Rm, ST, EfieldScl, IonFy_i, Ek_i, fRF, M);
}


struct ElementRFCavity : public Moment2ElementBase
{
    // Transport matrix for an RF Cavity.
    typedef Moment2ElementBase base_t;
    typedef typename base_t::state_t state_t;
    ElementRFCavity(const Config& c)
        :base_t(c)
    {
        std::string cav_type = c.get<std::string>("cavtype");
        double L             = c.get<double>("L")*MtoMM,         // Convert from [m] to [mm].
               phi_ref       = c.get<double>("phi_ref", 0e0);

        this->transfer_raw(state_t::PS_X, state_t::PS_PX) = L;
        this->transfer_raw(state_t::PS_Y, state_t::PS_PY) = L;
        // For total path length.
//        this->transfer(state_t::PS_S, state_t::PS_S)  = L;
    }
    virtual ~ElementRFCavity() {}

    virtual void recompute_matrix(state_t& ST)
    {
        transfer = transfer_raw;

        last_Kenergy_in = ST.Ekinetic;

        PropagateLongRFCav(conf(), ST);

        last_Kenergy_out = ST.Ekinetic;

        // Define initial conditions.
        double accIonW, avebeta, avegamma;

        InitRFCav(conf(), ST, accIonW, avebeta, avegamma, transfer);
   }

    virtual const char* type_name() const {return "rfcavity";}
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
