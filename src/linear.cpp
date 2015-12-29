
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
        double len = c.get<double>("length");
        this->transfer(state_t::L_X, state_t::P_X) = len;
        this->transfer(state_t::L_Y, state_t::P_Y) = len;
        this->transfer(state_t::L_Z, state_t::L_Z) = len;
    }
    virtual ~ElementDrift() {}

    virtual const char* type_name() const {return "drift";}
};

template<typename Base>
struct ElementThinDipole : public Base
{
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementThinDipole(const Config& c)
        :base_t(c)
    {
        double angle = c.get<double>("angle"), // in rad.
               P = c.get<double>("radius", 1.0),
               off = c.get<double>("vertical", 0.0)!=0.0 ? state_t::L_Y : state_t::L_X ,
               cos = ::cos(angle),
               sin = ::sin(angle);

        this->transfer(off,off) = this->transfer(off+1,off+1) = cos;
        this->transfer(off,off+1) = P*sin;
        this->transfer(off+1,off) = -sin/P;
    }
    virtual ~ElementThinDipole() {}

    virtual const char* type_name() const {return "dipole";}
};

template<typename Base>
struct ElementThinQuad : public Base
{
    typedef Base base_t;
    typedef typename base_t::state_t state_t;
    ElementThinQuad(const Config& c)
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
            Fdir = state_t::L_Y;
            Ddir = state_t::L_X;
        } else {
            // focus in X, defocus in Y
            Fdir = state_t::L_X;
            Ddir = state_t::L_Y;
        }

        this->transfer(Fdir,Fdir) = this->transfer(Fdir+1,Fdir+1) = cos;
        this->transfer(Fdir,Fdir+1) = sin/sK;
        this->transfer(Fdir+1,Fdir) = sK*sin;

        this->transfer(Ddir,Ddir) = this->transfer(Ddir+1,Ddir+1) = cosh;
        this->transfer(Ddir,Ddir+1) = sinh/sK;
        this->transfer(Ddir+1,Ddir) = sK*sinh;
    }
    virtual ~ElementThinQuad() {}

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

    Machine::registerElement<ElementSource<LinearElementBase<VectorState> > >("Vector", "source");
    Machine::registerElement<ElementSource<LinearElementBase<MatrixState> > >("TransferMatrix", "source");
    Machine::registerElement<ElementSource<MomentElementBase> >("MomentMatrix", "source");

    Machine::registerElement<ElementDrift<LinearElementBase<VectorState> > >("Vector", "drift");
    Machine::registerElement<ElementDrift<LinearElementBase<MatrixState> > >("TransferMatrix", "drift");
    Machine::registerElement<ElementDrift<MomentElementBase> >("MomentMatrix", "drift");

    Machine::registerElement<ElementThinDipole<LinearElementBase<VectorState> > >("Vector", "dipole");
    Machine::registerElement<ElementThinDipole<LinearElementBase<MatrixState> > >("TransferMatrix", "dipole");
    Machine::registerElement<ElementThinDipole<MomentElementBase> >("MomentMatrix", "dipole");

    Machine::registerElement<ElementThinQuad<LinearElementBase<VectorState> > >("Vector", "quad");
    Machine::registerElement<ElementThinQuad<LinearElementBase<MatrixState> > >("TransferMatrix", "quad");
    Machine::registerElement<ElementThinQuad<MomentElementBase> >("MomentMatrix", "quad");

    Machine::registerElement<ElementGeneric<LinearElementBase<VectorState> > >("Vector", "generic");
    Machine::registerElement<ElementGeneric<LinearElementBase<MatrixState> > >("TransferMatrix", "generic");
    Machine::registerElement<ElementGeneric<MomentElementBase> >("MomentMatrix", "generic");
}
