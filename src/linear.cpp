
#include <string>
#include <algorithm>

#include "scsi/linear.h"
#include "scsi/moment.h"
#include "scsi/state/vector.h"
#include "scsi/state/matrix.h"

MatrixState::MatrixState(const Config& c)
    :StateBase(c)
    ,state(boost::numeric::ublas::identity_matrix<double>(6))
{
    try{
        const std::vector<double>& I = c.get<std::vector<double> >("initial");
        if(I.size()<state.data().size())
            throw std::invalid_argument("Initial state size too big");
        std::copy(I.begin(), I.end(), state.data().begin());
    }catch(key_error&){
        // default to identity
    }catch(boost::bad_any_cast&){
        throw std::invalid_argument("'initial' has wrong type (must be vector)");
    }
}

MatrixState::~MatrixState() {}

void MatrixState::show(std::ostream& strm) const
{
    strm<<"State: "<<state<<"\n";
}

bool MatrixState::getArray(unsigned idx, ArrayInfo& Info) {
    if(idx==0) {
        Info.name = "state";
        Info.ptr = &state(0,0);
        Info.ndim = 2;
        Info.dim[0] = state.size1();
        Info.dim[1] = state.size2();
        return true;
    }
    return false;
}

VectorState::VectorState(const Config& c)
    :StateBase(c)
    ,state(6, 0.0)
{
    try{
        const std::vector<double>& I = c.get<std::vector<double> >("initial");
        if(I.size()<state.size())
            throw std::invalid_argument("Initial state size too big");
        std::copy(I.begin(), I.end(), state.begin());
    }catch(key_error&){
    }catch(boost::bad_any_cast&){
    }
}

VectorState::~VectorState() {}

void VectorState::show(std::ostream& strm) const
{
    strm<<"State: "<<state<<"\n";
}

bool VectorState::getArray(unsigned idx, ArrayInfo& Info) {
    if(idx==0) {
        Info.name = "state";
        Info.ptr = &state(0);
        Info.ndim = 1;
        Info.dim[0] = state.size();
        return true;
    }
    return false;
}

namespace {

template<typename Base>
void Get2by2Matrix(const double L, const double K, const unsigned ind, typename Base::value_t &M)
{
    // Transport matrix for one plane for a Quadrupole.
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

template<typename Base>
struct ElementSource : public Base
{
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementSource(const Config& c)
        :base_t(c), istate(c)
    {}

    virtual void advance(StateBase& s) const
    {
        state_t& ST = static_cast<state_t&>(s);
        // Replace state with our initial values
        ST.state = this->istate.state;
        ST.ionZ  = this->istate.ionZ;
        ST.ionEs = this->istate.ionEs;
        ST.ionEk = this->istate.ionEk;
        ST.ionW  = this->istate.ionW;
        ST.Brho  = this->istate.Brho;
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

template<typename Base>
struct ElementMark : public Base
{
    // Transport (identity) matrix for a Marker.
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementMark(const Config& c)
        :base_t(c)
    {
        // Identity matrix.
    }
    virtual ~ElementMark() {}

    virtual const char* type_name() const {return "marker";}
};

template<typename Base>
struct ElementDrift : public Base
{
    // Transport matrix for a Drift.
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementDrift(const Config& c)
        :base_t(c)
    {
        double L = c.get<double>("L");

        this->transfer(state_t::PS_X, state_t::PS_PX) = L;
        this->transfer(state_t::PS_Y, state_t::PS_PY) = L;
        this->transfer(state_t::PS_S, state_t::PS_S)  = L;
    }
    virtual ~ElementDrift() {}

    virtual const char* type_name() const {return "drift";}
};

template<typename Base>
struct ElementSBend : public Base
{
    // Transport matrix for a Gradient Sector Bend (cylindrical coordinates).
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementSBend(const Config& c)
        :base_t(c)
    {
        double L   = c.get<double>("L",   0e0),
               phi = c.get<double>("phi", 0e0), // [rad].
               rho = L/phi,
               K   = c.get<double>("K",   0e0), // [1/m^2].
               Kx  = K + 1e0/sqr(rho),
               Ky  = -K;

        // Horizontal plane.
        Get2by2Matrix<Base>(L, Kx, (unsigned)state_t::PS_X, this->transfer);
        // Vertical plane.
        Get2by2Matrix<Base>(L, Ky, (unsigned)state_t::PS_Y, this->transfer);
        // Longitudinal plane.
        this->transfer(state_t::PS_S,  state_t::PS_S) = L;
    }
    virtual ~ElementSBend() {}

    virtual const char* type_name() const {return "sbend";}
};

template<typename Base>
struct ElementQuad : public Base
{
    // Transport matrix for a Quadrupole (Cartesian coordinates).
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementQuad(const Config& c)
        :base_t(c)
    {
        double L    = c.get<double>("L"),
               K    = c.get<double>("K", 0e0);

        // Horizontal plane.
        Get2by2Matrix<Base>(L,  K, (unsigned)state_t::PS_X, this->transfer);
        // Vertical plane.
        Get2by2Matrix<Base>(L, -K, (unsigned)state_t::PS_Y, this->transfer);
        // Longitudinal plane.
        this->transfer(state_t::PS_S, state_t::PS_S) = L;
    }
    virtual ~ElementQuad() {}

    virtual const char* type_name() const {return "quadrupole";}
};

template<typename Base>
struct ElementSolenoid : public Base
{
    // Transport (identity) matrix for a Solenoid; K = B0/(B*rho).
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementSolenoid(const Config& c)
        :base_t(c)
    {
        double L = c.get<double>("L"),
               B = c.get<double>("B"),
               K = B/c.get<double>("Brho"),
               C = ::cos(K*L),
               S = ::sin(K*L);

        this->transfer(state_t::PS_X, state_t::PS_X)
                = this->transfer(state_t::PS_PX, state_t::PS_PX)
                = this->transfer(state_t::PS_Y, state_t::PS_Y)
                = this->transfer(state_t::PS_PY, state_t::PS_PY)
                = sqr(C);

        if (K != 0e0)
            this->transfer(state_t::PS_X, state_t::PS_PX) = S*C/K;
        else
            this->transfer(state_t::PS_X, state_t::PS_PX) = L;
        this->transfer(state_t::PS_X, state_t::PS_Y) = S*C;
        if (K != 0e0)
            this->transfer(state_t::PS_X, state_t::PS_PY) = sqr(S)/K;
        else
            this->transfer(state_t::PS_X, state_t::PS_PY) = 0e0;

        this->transfer(state_t::PS_PX, state_t::PS_X) = -K*S*C;
        this->transfer(state_t::PS_PX, state_t::PS_Y) = -K*sqr(S);
        this->transfer(state_t::PS_PX, state_t::PS_PY) = S*C;

        this->transfer(state_t::PS_Y, state_t::PS_X) = -S*C;
        if (K != 0e0)
            this->transfer(state_t::PS_Y, state_t::PS_PX) = -sqr(S)/K;
        else
            this->transfer(state_t::PS_Y, state_t::PS_PX) = 0e0;
        if (K != 0e0)
            this->transfer(state_t::PS_Y, state_t::PS_PY) = S*C/K;
        else
            this->transfer(state_t::PS_Y, state_t::PS_PY) = L;

        this->transfer(state_t::PS_PY, state_t::PS_X) = K*sqr(S);
        this->transfer(state_t::PS_PY, state_t::PS_PX) = -S*C;
        this->transfer(state_t::PS_PY, state_t::PS_Y) = -K*S*C;

        // Longitudinal plane.
        this->transfer(state_t::PS_S, state_t::PS_S) = L;
    }
    virtual ~ElementSolenoid() {}

    virtual const char* type_name() const {return "solenoid";}
};

template<typename Base>
struct ElementRFCavity : public Base
{
    // Transport matrix for an RF Cavity.
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementRFCavity(const Config& c)
        :base_t(c)
    {
        std::string cav_type = c.get<std::string>("cavtype");
        double L             = c.get<double>("L"),
        phi                  = c.get<double>("phi");

        this->transfer(state_t::PS_X, state_t::PS_PX) = L;
        this->transfer(state_t::PS_Y, state_t::PS_PY) = L;
//        this->transfer(state_t::PS_S, state_t::PS_S)  = L;
    }
    virtual ~ElementRFCavity() {}

    virtual const char* type_name() const {return "rfcavity";}
};

template<typename Base>
struct ElementEDipole : public Base
{
    // Transport matrix for an Electric Dipole.
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementEDipole(const Config& c)
        :base_t(c)
    {
        double L = c.get<double>("L");

    }
    virtual ~ElementEDipole() {}

    virtual const char* type_name() const {return "edipole";}
};

template<typename Base>
struct ElementGeneric : public Base
{
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementGeneric(const Config& c)
        :base_t(c)
    {
        std::vector<double> I = c.get<std::vector<double> >("transfer");
        if(I.size()<this->transfer.data().size())
            throw std::invalid_argument("Initial transfer size too big");
        std::copy(I.begin(), I.end(), this->transfer.data().begin());
    }
    virtual ~ElementGeneric() {}

    virtual const char* type_name() const {return "generic";}
};

} // namespace

void registerLinear()
{
    Machine::registerState<VectorState>("Vector");
    Machine::registerState<MatrixState>("TransferMatrix");
    Machine::registerState<MatrixState>("MomentMatrix");

    Machine::registerElement<ElementSource<LinearElementBase<VectorState> > >("Vector",         "source");
    Machine::registerElement<ElementSource<LinearElementBase<MatrixState> > >("TransferMatrix", "source");
    Machine::registerElement<ElementSource<MomentElementBase>               >("MomentMatrix",   "source");

    Machine::registerElement<ElementMark<LinearElementBase<VectorState> > >("Vector",         "marker");
    Machine::registerElement<ElementMark<LinearElementBase<MatrixState> > >("TransferMatrix", "marker");
    Machine::registerElement<ElementMark<MomentElementBase>               >("MomentMatrix",   "marker");

    Machine::registerElement<ElementDrift<LinearElementBase<VectorState> > >("Vector",         "drift");
    Machine::registerElement<ElementDrift<LinearElementBase<MatrixState> > >("TransferMatrix", "drift");
    Machine::registerElement<ElementDrift<MomentElementBase>               >("MomentMatrix",   "drift");

    Machine::registerElement<ElementSBend<LinearElementBase<VectorState> > >("Vector",         "sbend");
    Machine::registerElement<ElementSBend<LinearElementBase<MatrixState> > >("TransferMatrix", "sbend");
    Machine::registerElement<ElementSBend<MomentElementBase>               >("MomentMatrix",   "sbend");

    Machine::registerElement<ElementQuad<LinearElementBase<VectorState> > >("Vector",         "quadrupole");
    Machine::registerElement<ElementQuad<LinearElementBase<MatrixState> > >("TransferMatrix", "quadrupole");
    Machine::registerElement<ElementQuad<MomentElementBase>               >("MomentMatrix",   "quadrupole");

    Machine::registerElement<ElementSolenoid<LinearElementBase<VectorState> > >("Vector",         "solenoid");
    Machine::registerElement<ElementSolenoid<LinearElementBase<MatrixState> > >("TransferMatrix", "solenoid");
    Machine::registerElement<ElementSolenoid<MomentElementBase>               >("MomentMatrix",   "solenoid");

    Machine::registerElement<ElementRFCavity<LinearElementBase<VectorState> > >("Vector",         "rfcavity");
    Machine::registerElement<ElementRFCavity<LinearElementBase<MatrixState> > >("TransferMatrix", "rfcavity");
    Machine::registerElement<ElementRFCavity<MomentElementBase>               >("MomentMatrix",   "rfcavity");

    Machine::registerElement<ElementEDipole<LinearElementBase<VectorState> > >("Vector",         "edipole");
    Machine::registerElement<ElementEDipole<LinearElementBase<MatrixState> > >("TransferMatrix", "edipole");
    Machine::registerElement<ElementEDipole<MomentElementBase>               >("MomentMatrix",   "edipole");

    Machine::registerElement<ElementGeneric<LinearElementBase<VectorState> > >("Vector",         "generic");
    Machine::registerElement<ElementGeneric<LinearElementBase<MatrixState> > >("TransferMatrix", "generic");
    Machine::registerElement<ElementGeneric<MomentElementBase>               >("MomentMatrix",   "generic");
}
