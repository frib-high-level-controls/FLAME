
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
void get(typename Base::value_t &M)
{
    M(0, 0) = 1e0;
}

template<typename Base>
struct ElementSource : public Base
{
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementSource(const Config& c)
        :base_t(c)
        ,ivect(c.get<std::vector<double> >("initial",
                                           std::vector<double>()))
    {}

    virtual void advance(StateBase& s) const
    {
        state_t& ST = static_cast<state_t&>(s);
        if(ivect.size()==0)
            return; // use defaults
        // Replace state with our initial values
        if(ST.state.data().size()!=ivect.size())
            throw std::invalid_argument("Initial state size incorrect");
        std::copy(ivect.begin(), ivect.end(), ST.state.data().begin());
    }

    std::vector<double> ivect;

    virtual ~ElementSource() {}

    virtual const char* type_name() const {return "source";}
};

template<typename Base>
struct ElementDrift : public Base
{
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementDrift(const Config& c)
        :base_t(c)
    {
        double len = c.get<double>("L");
        this->transfer(state_t::PS_X, state_t::PS_PX) = len;
        this->transfer(state_t::PS_Y, state_t::PS_PY) = len;
        this->transfer(state_t::PS_S, state_t::PS_S)  = len;
    }
    virtual ~ElementDrift() {}

    virtual const char* type_name() const {return "drift";}
};

template<typename Base>
struct ElementSBend : public Base
{
    // Sector Bend.
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementSBend(const Config& c)
        :base_t(c)
    {
        double L      = c.get<double>("L", 0e0),
               phi    = c.get<double>("phi", 0e0), // [rad].
               K      = c.get<double>("K", 0e0),   // [1/m^2].
               rho    = L/phi,
               Kx     = K + 1e0/sqr(rho),
               Ky     = -K,
               sqrtK,
               psi,
               cs,
               sn;

        // Horizontal plane.
        if (Kx > 0e0) {
            sqrtK = sqrt(Kx);
            psi = sqrtK*L;
            cs = ::cos(fabs(psi));
            sn = ::sin(fabs(psi));
            this->transfer(state_t::PS_X,  state_t::PS_X)  = cs;
            this->transfer(state_t::PS_X,  state_t::PS_PX) = sn/sqrtK;
            this->transfer(state_t::PS_PX, state_t::PS_X)  = -sqrtK*sn;
            this->transfer(state_t::PS_PX, state_t::PS_PX) = cs;
        } else {
            sqrtK = sqrt(-Kx);
            psi = sqrtK*L;
            cs = ::cosh(fabs(psi));
            sn = ::sinh(fabs(psi));
            this->transfer(state_t::PS_X,  state_t::PS_X)  = cs;
            this->transfer(state_t::PS_X,  state_t::PS_PX) = sn/sqrtK;
            this->transfer(state_t::PS_PX, state_t::PS_X)  = sqrtK*sn;
            this->transfer(state_t::PS_PX, state_t::PS_PX) = cs;
        }

        // Vertical plane.
        if (Ky > 0e0) {
            sqrtK = sqrt(Ky);
            psi = sqrtK*L;
            cs = ::cos(fabs(psi));
            sn = ::sin(fabs(psi));
            this->transfer(state_t::PS_Y,  state_t::PS_Y)  = cs;
            this->transfer(state_t::PS_Y,  state_t::PS_PY) = sn/sqrtK;
            this->transfer(state_t::PS_PY, state_t::PS_Y)  = -sqrtK*sn;
            this->transfer(state_t::PS_PY, state_t::PS_PY) = cs;
        } else {
            sqrtK = sqrt(-Ky);
            psi = sqrtK*L;
            cs = ::cosh(fabs(psi));
            sn = ::sinh(fabs(psi));
            this->transfer(state_t::PS_Y,  state_t::PS_Y)  = cs;
            this->transfer(state_t::PS_Y,  state_t::PS_PY) = sn/sqrtK;
            this->transfer(state_t::PS_PY, state_t::PS_Y)  = sqrtK*sn;
            this->transfer(state_t::PS_PY, state_t::PS_PY) = cs;
        }
        // Longitudinal plane.
        this->transfer(state_t::PS_S,  state_t::PS_S) = L;
    }
    virtual ~ElementSBend() {}

    virtual const char* type_name() const {return "sbend";}
};

template<typename Base>
struct ElementQuad : public Base
{
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementQuad(const Config& c)
        :base_t(c)
    {
        double L    = c.get<double>("L"),
               K    = c.get<double>("K", 0e0),
               aK   = fabs(K),
               sK   = sqrt(aK),
               sKL  = sK*L,
               cos  = ::cos(sKL),
               sin  = ::sin(sKL),
               cosh = ::cosh(sKL),
               sinh = ::sinh(sKL);
        unsigned Find, Dind;

        if(K < 0e0) {
            // defocus in X, focus in Y
            Find = state_t::PS_Y;
            Dind = state_t::PS_X;
        } else {
            // focus in X, defocus in Y
            Find = state_t::PS_X;
            Dind = state_t::PS_Y;
        }

        this->transfer(Find, Find) = this->transfer(Find+1, Find+1) = cos;
        if (sK != 0e0)
            this->transfer(Find, Find+1) = sin/sK;
        else
            this->transfer(Find, Find+1) = L;
        if (sK != 0e0)
            this->transfer(Find+1, Find) = -sK*sin;
        else
            this->transfer(Find+1, Find) = 0e0;
        this->transfer(Dind, Dind) = this->transfer(Dind+1, Dind+1) = cosh;
        if (sK != 0e0)
            this->transfer(Dind, Dind+1) = sinh/sK;
        else
            this->transfer(Dind, Dind+1) = L;
        if (sK != 0e0)
            this->transfer(Dind+1, Dind) = sK*sinh;
        else
            this->transfer(Dind+1, Dind) = 0e0;

        this->transfer(state_t::PS_S, state_t::PS_S) = L;
    }
    virtual ~ElementQuad() {}

    virtual const char* type_name() const {return "quad";}
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

    Machine::registerElement<ElementDrift<LinearElementBase<VectorState> > >("Vector",         "drift");
    Machine::registerElement<ElementDrift<LinearElementBase<MatrixState> > >("TransferMatrix", "drift");
    Machine::registerElement<ElementDrift<MomentElementBase>               >("MomentMatrix",   "drift");

    Machine::registerElement<ElementSBend<LinearElementBase<VectorState> > >("Vector",         "sbend");
    Machine::registerElement<ElementSBend<LinearElementBase<MatrixState> > >("TransferMatrix", "sbend");
    Machine::registerElement<ElementSBend<MomentElementBase>               >("MomentMatrix",   "sbend");

    Machine::registerElement<ElementQuad<LinearElementBase<VectorState> > >("Vector",         "quad");
    Machine::registerElement<ElementQuad<LinearElementBase<MatrixState> > >("TransferMatrix", "quad");
    Machine::registerElement<ElementQuad<MomentElementBase>               >("MomentMatrix",   "quad");

    Machine::registerElement<ElementGeneric<LinearElementBase<VectorState> > >("Vector",         "generic");
    Machine::registerElement<ElementGeneric<LinearElementBase<MatrixState> > >("TransferMatrix", "generic");
    Machine::registerElement<ElementGeneric<MomentElementBase>               >("MomentMatrix",   "generic");
}
