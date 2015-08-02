
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
struct GenericMomentElement : MomentElementBase
{
    GenericMomentElement(const Config& c)
        :MomentElementBase(c)
    {
        try{
            const std::vector<double>& I = c.get<const std::vector<double>&>("initial");
            if(I.size()<transfer.data().size())
                throw std::invalid_argument("Initial transfer size too big");
            std::copy(I.begin(), I.end(), transfer.data().begin());
        }catch(key_error&){
        }catch(boost::bad_any_cast&){
        }
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
