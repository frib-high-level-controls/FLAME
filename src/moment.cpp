
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

    state_t::value_t temp;
    noalias(temp) = prod(transfer, ST.state);
    ST.state = prod(temp, trans(transfer));
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

    Machine::registerElement<GenericMomentElement >("MomentMatrix", "generic");
}
