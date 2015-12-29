
#include "scsi/moment.h"

MomentElementBase::MomentElementBase(const Config& c)
    :ElementVoid(c)
    ,transfer(boost::numeric::ublas::identity_matrix<double>(6))
{}

MomentElementBase::~MomentElementBase() {}

void MomentElementBase::show(std::ostream& strm) const
{
    using namespace boost::numeric::ublas;
    ElementVoid::show(strm);
    strm<<"Transfer: "<<transfer<<"\n";
    strm<<"TransferT: "<<trans(transfer)<<"\n";
}

void MomentElementBase::advance(StateBase& s) const
{
    state_t& ST = static_cast<state_t&>(s);
    using namespace boost::numeric::ublas;

    state_t::value_t temp(state_t::maxsize,state_t::maxsize); //TODO: attach to some obj to avoid realloc?
    noalias(temp) = prod(transfer, ST.state);
    noalias(ST.state) = prod(temp, trans(transfer));
}

namespace {
struct MomentSource : MomentElementBase
{
    MomentSource(const Config& c)
        :MomentElementBase(c)
        ,ivect(c.get<std::vector<double> >("initial",
                                           std::vector<double>()))
    {

    }

    virtual ~MomentSource() {}

    virtual void advance(StateBase& s) const
    {
        MomentElementBase::state_t& ST = static_cast<MomentElementBase::state_t&>(s);
        if(ivect.size()==0)
            return; // use defaults
        // Replace state with our initial values
        if(ST.state.data().size()!=ivect.size())
            throw std::invalid_argument("Initial state size incorrect");
        std::copy(ivect.begin(), ivect.end(), ST.state.data().begin());
    }

    std::vector<double> ivect;

    virtual const char* type_name() const {return "source";}
};

struct MomentDipole : MomentElementBase
{
    MomentDipole(const Config& c)
        :MomentElementBase(c)
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
    virtual ~MomentDipole() {}

    virtual const char* type_name() const {return "drift";}
};

struct MomentDrift : MomentElementBase
{
    MomentDrift(const Config& c)
        :MomentElementBase(c)
    {
        this->transfer(state_t::L_X, state_t::P_X) = c.get<double>("length");
        this->transfer(state_t::L_Y, state_t::P_Y) = c.get<double>("length");
        this->transfer(state_t::L_Z, state_t::L_Z) = c.get<double>("length");
    }
    virtual ~MomentDrift() {}

    virtual const char* type_name() const {return "drift";}
};

struct GenericMomentElement : MomentElementBase
{
    GenericMomentElement(const Config& c)
        :MomentElementBase(c)
    {
        const std::vector<double>& I = c.get<std::vector<double> >("transfer");
        if(I.size()<transfer.data().size())
            throw std::invalid_argument("Initial transfer size too big");
        std::copy(I.begin(), I.end(), transfer.data().begin());
    }
    virtual ~GenericMomentElement() {}

    virtual const char* type_name() const {return "generic";}
};
}

void registerMoment()
{
    Machine::registerState<MatrixState>("MomentMatrix");

    Machine::registerElement<GenericMomentElement>("MomentMatrix", "generic");
    Machine::registerElement<MomentSource>("MomentMatrix", "source");
    Machine::registerElement<MomentDrift>("MomentMatrix", "drift");
    Machine::registerElement<MomentDipole>("MomentMatrix", "dipole");
}
