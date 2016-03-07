
#include "scsi/moment.h"

MomentState::MomentState(const Config& c)
    :StateBase(c)
    ,moment0(maxsize, 0.0)
    ,state(boost::numeric::ublas::identity_matrix<double>(maxsize))
{
    try{
        const std::vector<double>& I = c.get<std::vector<double> >("moment0");
        if(I.size()>moment0.size())
            throw std::invalid_argument("Initial state size too big");
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
}

MomentState::~MomentState() {}

void MomentState::assign(const StateBase& other)
{
    const MomentState *O = dynamic_cast<const MomentState*>(&other);
    if(!O)
        throw std::invalid_argument("Can't assign State: incompatible types");
    moment0 = O->moment0;
    state = O->state;
    StateBase::assign(other);
}

void MomentState::show(std::ostream& strm) const
{
    strm<<"State: moment0="<<moment0<<" state="<<state<<"\n";
}

bool MomentState::getArray(unsigned idx, ArrayInfo& Info) {
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
    }
    return StateBase::getArray(idx-2, Info);
}

MomentElementBase::MomentElementBase(const Config& c)
    :ElementVoid(c)
    ,transfer(boost::numeric::ublas::identity_matrix<double>(state_t::maxsize))
    ,scratch(state_t::maxsize, state_t::maxsize)
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

    ST.moment0 = prod(transfer, ST.moment0);

    noalias(scratch) = prod(transfer, ST.state);
    noalias(ST.state) = prod(scratch, trans(transfer));
}

void registerMoment()
{
}
