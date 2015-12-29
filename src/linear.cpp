
#include <algorithm>

#include "scsi/linear.h"
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

template<typename State>
struct LinearSource : LinearElementBase<State>
{
    typedef LinearElementBase<State> base_t;
    LinearSource(const Config& c)
        :base_t(c)
        ,ivect(c.get<std::vector<double> >("initial",
                                           std::vector<double>()))
    {}

    virtual void advance(StateBase& s) const
    {
        State& ST = static_cast<State&>(s);
        if(ivect.size()==0)
            return; // use defaults
        // Replace state with our initial values
        if(ST.state.data().size()!=ivect.size())
            throw std::invalid_argument("Initial state size incorrect");
        std::copy(ivect.begin(), ivect.end(), ST.state.data().begin());
    }

    std::vector<double> ivect;

    virtual ~LinearSource() {}

    virtual const char* type_name() const {return "source";}
};

template<typename State>
struct LinearDrift : LinearElementBase<State>
{
    typedef LinearElementBase<State> base_t;
    LinearDrift(const Config& c)
        :base_t(c)
    {
        this->transfer(State::L_X, State::P_X) = c.get<double>("length");
        this->transfer(State::L_Y, State::P_Y) = c.get<double>("length");
        this->transfer(State::L_Z, State::L_Z) = c.get<double>("length");
    }
    virtual ~LinearDrift() {}

    virtual const char* type_name() const {return "drift";}
};

template<typename State>
struct LinearThinDipole : LinearElementBase<State>
{
    typedef LinearElementBase<State> base_t;
    LinearThinDipole(const Config& c)
        :base_t(c)
    {
        double angle = c.get<double>("angle"), // in rad.
               P = c.get<double>("radius", 1.0),
               off = c.get<double>("vertical", 0.0)!=0.0 ? State::L_Y : State::L_X ,
               cos = ::cos(angle),
               sin = ::sin(angle);

        this->transfer(off,off) = this->transfer(off+1,off+1) = cos;
        this->transfer(off,off+1) = P*sin;
        this->transfer(off+1,off) = -sin/P;
    }
    virtual ~LinearThinDipole() {}

    virtual const char* type_name() const {return "dipole";}
};

template<typename State>
struct LinearThinQuad : LinearElementBase<State>
{
    typedef LinearElementBase<State> base_t;
    LinearThinQuad(const Config& c)
        :base_t(c)
    {
        double L  = c.get<double>("length"),
               K  = c.get<double>("strength", 1.0),
               aK = fabs(K),
               sK = sqrt(aK),
               sKL=sK*L,
               cos = ::cos(sKL),
               sin = ::sin(sKL),
               cosh = ::cosh(sKL),
               sinh = ::sinh(sKL);
        unsigned Fdir, Ddir;

        if(K<0.0) {
            // defocus in X, focus in Y
            Fdir = State::L_Y;
            Ddir = State::L_X;
        } else {
            // focus in X, defocus in Y
            Fdir = State::L_X;
            Ddir = State::L_Y;
        }

        this->transfer(Fdir,Fdir) = this->transfer(Fdir+1,Fdir+1) = cos;
        this->transfer(Fdir,Fdir+1) = sin/sK;
        this->transfer(Fdir+1,Fdir) = sK*sin;

        this->transfer(Ddir,Ddir) = this->transfer(Ddir+1,Ddir+1) = cosh;
        this->transfer(Ddir,Ddir+1) = sinh/sK;
        this->transfer(Ddir+1,Ddir) = sK*sinh;
    }
    virtual ~LinearThinQuad() {}

    virtual const char* type_name() const {return "quad";}
};
template<typename State>
struct LinearGeneric : LinearElementBase<State>
{
    typedef LinearElementBase<State> base_t;
    LinearGeneric(const Config& c)
        :base_t(c)
    {
        std::vector<double> I = c.get<std::vector<double> >("transfer");
        if(I.size()<this->transfer.data().size())
            throw std::invalid_argument("Initial transfer size too big");
        std::copy(I.begin(), I.end(), this->transfer.data().begin());
    }
    virtual ~LinearGeneric() {}

    virtual const char* type_name() const {return "generic";}
};

} // namespace

void registerLinear()
{
    Machine::registerState<VectorState>("Vector");
    Machine::registerState<MatrixState>("TransferMatrix");

    Machine::registerElement<LinearSource<VectorState> >("Vector", "source");
    Machine::registerElement<LinearSource<MatrixState> >("TransferMatrix", "source");

    Machine::registerElement<LinearDrift<VectorState> >("Vector", "drift");
    Machine::registerElement<LinearDrift<MatrixState> >("TransferMatrix", "drift");

    Machine::registerElement<LinearThinDipole<VectorState> >("Vector", "dipole");
    Machine::registerElement<LinearThinDipole<MatrixState> >("TransferMatrix", "dipole");

    Machine::registerElement<LinearThinQuad<VectorState> >("Vector", "quad");
    Machine::registerElement<LinearThinQuad<MatrixState> >("TransferMatrix", "quad");

    Machine::registerElement<LinearGeneric<VectorState> >("Vector", "generic");
    Machine::registerElement<LinearGeneric<MatrixState> >("TransferMatrix", "generic");
}
