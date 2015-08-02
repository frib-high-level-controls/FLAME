
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
struct NOOPMomentElement : MomentElementBase
{
    NOOPMomentElement(const Config& c)
        :MomentElementBase(c)
    {
        // leave the default (identity)
    }

    virtual const char* type_name() const {return "noop";}
};
}

void registerMoment()
{
    Machine::registerState<MatrixState>("MomentMatrix");

    Machine::registerElement<NOOPMomentElement >("MomentMatrix", "noop");
}
