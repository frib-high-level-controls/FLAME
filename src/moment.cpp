
#include "scsi/moment.h"

MomentState::MomentState(const Config& c)
    :StateBase(c)
    ,IonEs(c.get<double>("IonEs", 0))
    ,IonEk(c.get<double>("IonEk", 0))
    ,IonW(c.get<double >("IonW", 0))
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

    IonEs          = c.get<double>("IonEs", 0e0), // Rest energy.
    IonEk          = c.get<double>("IonEk", 0e0);
    IonZ           = c.get<double>("IonZ", 0e0);

    IonW           = IonEs + IonEk;
}

MomentState::~MomentState() {}

MomentState::MomentState(const MomentState& o, clone_tag t)
    :StateBase(o, t)
    ,IonZ(o.IonZ)
    ,IonEs(o.IonEs)
    ,IonEk(o.IonEk)
    ,IonW(o.IonW)
    ,moment0(o.moment0)
    ,state(o.state)
{}

void MomentState::assign(const StateBase& other)
{
    const MomentState *O = dynamic_cast<const MomentState*>(&other);
    if(!O)
        throw std::invalid_argument("Can't assign State: incompatible types");
    IonZ     = O->IonZ;
    IonEs    = O->IonEs;
    IonEk    = O->IonEk;
    IonW     = O->IonW;
    moment0  = O->moment0;
    state    = O->state;
    StateBase::assign(other);
}

void MomentState::show(std::ostream& strm, int level) const
{
    strm<<"State: moment0="<<moment0<<" state="<<state<<"\n";
}

bool MomentState::getArray(unsigned idx, ArrayInfo& Info) {
    unsigned I=0;
    if(idx==I++) {
        Info.name = "state";
        Info.ptr = &state(0,0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 2;
        Info.dim[0] = state.size1();
        Info.dim[1] = state.size2();
        return true;
    } else if(idx==I++) {
        Info.name = "moment0";
        Info.ptr = &moment0(0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim[0] = moment0.size();
        return true;
    } else if(idx==I++) {
        Info.name = "IonZ";
        Info.ptr = &IonZ;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "IonEs";
        Info.ptr = &IonEs;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "IonEk";
        Info.ptr = &IonEk;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "IonW";
        Info.ptr = &IonW;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    }
    return StateBase::getArray(idx-I, Info);
}

MomentElementBase::MomentElementBase(const Config& c)
    :ElementVoid(c)
    ,transfer(boost::numeric::ublas::identity_matrix<double>(state_t::maxsize))
    ,scratch(state_t::maxsize, state_t::maxsize)
{}

MomentElementBase::~MomentElementBase() {}

void MomentElementBase::show(std::ostream& strm, int level) const
{
    using namespace boost::numeric::ublas;
    ElementVoid::show(strm, level);
    strm<<"Transfer: "<<transfer<<"\n";
    strm<<"TransferT: "<<trans(transfer)<<"\n";
}

void MomentElementBase::advance(StateBase& s)
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
