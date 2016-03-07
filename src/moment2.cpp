
#include "scsi/moment2.h"

#include "scsi/h5loader.h"

Moment2State::Moment2State(const Config& c)
    :StateBase(c)
    ,energy(c.get<double>("energy", 0.0))
    ,do_recalc_energy(true)
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

Moment2State::~Moment2State() {}

void Moment2State::assign(const StateBase& other)
{
    const Moment2State *O = dynamic_cast<const Moment2State*>(&other);
    if(!O)
        throw std::invalid_argument("Can't assign State: incompatible types");
    energy = O->energy;
    do_recalc_energy = O->do_recalc_energy;
    moment0 = O->moment0;
    state = O->state;
    StateBase::assign(other);
}

void Moment2State::show(std::ostream& strm) const
{
    strm<<"State: energy="<<energy<<" moment0="<<moment0<<" state="<<state<<"\n";
}

bool Moment2State::getArray(unsigned idx, ArrayInfo& Info) {
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

Moment2ElementBase::Moment2ElementBase(const Config& c)
    :ElementVoid(c)
    ,do_recalc_energy(true)
    ,energy_gain(0.0)
    ,do_recalc_transfer(true)
    ,transfer(boost::numeric::ublas::identity_matrix<double>(state_t::maxsize))
    ,scratch(state_t::maxsize, state_t::maxsize)
{}

Moment2ElementBase::~Moment2ElementBase() {}

void Moment2ElementBase::show(std::ostream& strm) const
{
    using namespace boost::numeric::ublas;
    ElementVoid::show(strm);
    strm<<"Energy gain "<<energy_gain<<"\n";
    strm<<"Transfer: "<<transfer<<"\n";
    strm<<"TransferT: "<<trans(transfer)<<"\n";
}

void Moment2ElementBase::advance(StateBase& s)
{
    state_t& ST = static_cast<state_t&>(s);
    using namespace boost::numeric::ublas;

    if(do_recalc_energy || ST.do_recalc_energy) {
        do_recalc_energy = false;
        do_recalc_transfer = true;
        ST.do_recalc_energy = true;
        recalc_energy_gain(ST);
    }
    if(do_recalc_transfer) {
        do_recalc_transfer = false;
        recalc_transfer(ST);
    }

    ST.energy += energy_gain;
    ST.moment0 = prod(transfer, ST.moment0);

    noalias(scratch) = prod(transfer, ST.state);
    noalias(ST.state) = prod(scratch, trans(transfer));
}

void Moment2ElementBase::assign(const ElementVoid *other)
{
    const Moment2ElementBase *O = static_cast<const Moment2ElementBase*>(other);
    do_recalc_energy = O->do_recalc_energy;
    energy_gain = O->energy_gain;
    do_recalc_transfer = O->do_recalc_transfer;
    transfer = O->transfer;
    ElementVoid::assign(other);
}

namespace {

struct Moment2Passive : public Moment2ElementBase
{
    double L;

    Moment2Passive(const Config& conf)
        :Moment2ElementBase(conf)
        ,L(conf.get<double>("L"))
    {}

    virtual void recalc_energy_gain(state_t& s)
    {
        energy_gain = 0.0;
    }

    virtual void recalc_transfer(state_t& s)
    {
        transfer(0,0) = L*42.0 + energy_gain;
    }

    virtual const char* type_name() const {return "passive";}
};

struct Moment2RF : public Moment2ElementBase
{
    double L, foo;
    boost::numeric::ublas::matrix<double> fmap;

    Moment2RF(const Config& conf)
        :Moment2ElementBase(conf)
        ,L(conf.get<double>("L"))
        ,foo(conf.get<double>("foo"))
    {
        H5Loader loader(conf.get<std::string>("cavityfile"));
        fmap = loader.load("fieldmap");
    }

    virtual void recalc_energy_gain(state_t& s)
    { // recalc_energy()
        energy_gain = foo/2.0;
    }

    virtual void recalc_transfer(state_t& s)
    {
        transfer(0,0) = L*42.0 + energy_gain;
    }

    virtual const char* type_name() const {return "rfcavity";}
};

} // namespace

void registerMoment2()
{
    Machine::registerState<Moment2State>("MomentMatrix2");
    Machine::registerElement<Moment2Passive>("MomentMatrix2", "passive");
    Machine::registerElement<Moment2RF>("MomentMatrix2", "rfcaviy");
}
