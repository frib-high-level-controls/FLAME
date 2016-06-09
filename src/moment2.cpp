
#include <fstream>

#include <limits>

#include <boost/lexical_cast.hpp>

#include "scsi/constants.h"
#include "scsi/moment2.h"

#include "scsi/moment2.h"
#include "scsi/moment2_sup.h"
#include "scsi/rf_cavity.h"
#include "scsi/chg_stripper.h"

#include "scsi/h5loader.h"


namespace {

template<typename ARR>
bool load_storage(ARR& to, const Config& conf, const std::string& name, bool T=true)
{
    try{
        const std::vector<double>& val(conf.get<std::vector<double> >(name));
        if(to.size()!=val.size()) {
            std::ostringstream strm;
            strm<<"Array "<<name<<" must have "<<to.size()<<" elements, not "<<val.size();
            throw std::invalid_argument(strm.str());
        }
        std::copy(val.begin(), val.end(), to.begin());

        return true;
    }catch(key_error&){
        if(T)
            throw std::invalid_argument(name+" not defined");
        else
            return false;
        // default to identity
    }catch(boost::bad_any_cast&){
        if(T)
            throw std::invalid_argument(name+" has wrong type (must be vector)");
        else
            return false;
    }
}

} // namespace

std::ostream& operator<<(std::ostream& strm, const Particle& P)
{
    strm <<std::setprecision(8)<<std::setw(14)
      <<"IonZ="<<P.IonZ
      <<" IonQ="<<P.IonQ
      <<" IonEs="<<P.IonEs
      <<" IonEk="<<P.IonEk
      <<" SampleIonK="<<P.SampleIonK
      <<" phis="<<P.phis
      <<" IonW="<<P.IonW
      <<" gamma="<<P.gamma
      <<" beta="<<P.beta
      <<" bg="<<P.bg
      ;
    return strm;
}

Moment2State::Moment2State(const Config& c)
    :StateBase(c)
    ,ref()
    ,real()
    ,moment0_env(maxsize, 0e0)
    ,moment0_rms(maxsize, 0e0)
    ,moment1_env(boost::numeric::ublas::identity_matrix<double>(maxsize))
{
    // hack.  getArray() promises that returned pointers will remain valid for our lifetime.
    // This may not be true if std::vectors are resized.
    // Reserve for up to 10 states and hope for the best...
    // either need to provide real limit to max. states, or change getArray() iface
    real.reserve(10);

    double icstate_f = 0.0;
    bool have_cstate = c.tryGet<double>("cstate", icstate_f);
    size_t icstate = (size_t)icstate_f;

    std::string vectorname(c.get<std::string>("vector_variable", "moment0"));
    std::string matrixname(c.get<std::string>("matrix_variable", "initial"));

    std::vector<double> ics, nchg;
    bool have_ics = c.tryGet<std::vector<double> >("IonChargeStates", ics);

    ref.IonEs      = c.get<double>("IonEs", 0e0),
    ref.IonEk      = c.get<double>("IonEk", 0e0);
    ref.recalc();

    if(!have_ics) {
        ref.IonZ = c.get<double>("IonZ", 0e0);
        ref.IonQ = c.get<double>("IonQ", 1e0);

        ics.push_back(ref.IonZ);
        nchg.push_back(ref.IonQ);

    } else {
        if(ics.empty())
            throw std::invalid_argument("IonChargeStates w/ length 0");
        if(icstate>=ics.size())
            throw std::invalid_argument("IonChargeStates[cstate] is out of bounds");

        nchg = c.get<std::vector<double> >("NCharge");
        if(nchg.size()!=ics.size())
            throw std::invalid_argument("NCharge[] and IonChargeStates[] must have equal length");

        ref.IonZ = ics[0];
        ref.IonQ = nchg[0];
    }

    /* Possible configurations
     * 1. Neither 'cstate' nor 'IonChargeStates' defined (empty Config).
     *    No charge states, must go through source element to be useful
     * 2. 'IonChargeStates' defined, but not 'cstate'.
     *    Load all charge states
     * 3. 'cstate' and 'IonChargeStates' defined.
     *    Load a single charge state
     */
    if(!have_cstate && !have_ics) {
        // no-op
    } else if(!have_cstate && have_ics) {
        // many charge states
        icstate = 0;

    } else if(have_cstate && have_ics) {
        // single charge state

        // drop other than selected state
        ics[0]  = ics[icstate];
        nchg[0] = nchg[icstate];
        ics.resize(1);
        nchg.resize(1);

    } else {
        throw std::invalid_argument("Moment2State: must define IonChargeStates and NCharge when cstate is set");
    }

    if(have_ics) {
        real.resize(ics.size());
        moment0.resize(ics.size());
        moment1.resize(ics.size());

        double totalQ = 0.0;
        for(size_t i=0; i<ics.size(); i++) {
            std::string num(boost::lexical_cast<std::string>(icstate+i));

            moment0[i].resize(maxsize);
            moment1[i].resize(maxsize, maxsize);
            moment1[i] = boost::numeric::ublas::identity_matrix<double>(maxsize);

            load_storage(moment0[i].data(), c, vectorname+num);
            load_storage(moment1[i].data(), c, matrixname+num);

            real[i] = ref;

            real[i].IonZ = ics[i];
            real[i].IonQ = nchg[i];

            real[i].phis      = moment0[i][PS_S];
            real[i].IonEk    += moment0[i][PS_PS]*MeVtoeV;

            real[i].recalc();

            moment0_env += moment0[i]*real[i].IonQ;
            totalQ += real[i].IonQ;
        }

        moment0_env /= totalQ;
        moment1_env = moment1[0];
    } else {
        real.resize(1); // hack, ensure at least one element so getArray() can return some pointer
        real[0] = ref;

        moment0.resize(1);
        moment1.resize(1);
        moment0[0].resize(maxsize);
        moment1[0].resize(maxsize, maxsize);

        load_storage(moment0[0].data(), c, vectorname, false);
        load_storage(moment1[0].data(), c, matrixname, false);

        moment0_env = moment0[0];
        moment1_env = moment1[0];
    }

    calc_rms();
}

Moment2State::~Moment2State() {}

void Moment2State::calc_rms()
{
    // moment0_env already updated

    //TODO: avoid recalc of total charge
    double totQ = 0.0;
    for(size_t n=0; n<real.size(); n++) {
        totQ += real[n].IonQ;
    }

    for(size_t j=0; j<maxsize; j++) {
        double variance = 0.0;
        for(size_t n=0; n<moment0.size(); n++) {
            const double Q = real[n].IonQ;
            const double diff = moment0[n][j]-moment0_env[j];

            variance += Q*(moment1[n](j,j) + diff*diff); // Q * (sum squares + square of the sum)
        }

        moment0_rms[j] = sqrt(variance/totQ);
    }
}

Moment2State::Moment2State(const Moment2State& o, clone_tag t)
    :StateBase(o, t)
    ,ref(o.ref)
    ,real(o.real)
    ,moment0(o.moment0)
    ,moment1(o.moment1)
    ,moment0_env(o.moment0_env)
    ,moment0_rms(o.moment0_rms)
    ,moment1_env(o.moment1_env)
{}

void Moment2State::assign(const StateBase& other)
{
    const Moment2State *O = dynamic_cast<const Moment2State*>(&other);
    if(!O)
        throw std::invalid_argument("Can't assign State: incompatible types");
    ref     = O->ref;
    real    = O->real;
    moment0 = O->moment0;
    moment1 = O->moment1;
    moment0_env = O->moment0_env;
    moment0_rms = O->moment0_rms;
    moment1_env = O->moment1_env;
    StateBase::assign(other);
}

void Moment2State::show(std::ostream& strm, int level) const
{
    int j, k;

    if(real.empty()) {
        strm<<"\nState: empty\n";
        return;
    }

    strm << std::scientific << std::setprecision(8)
         << "\nState:\n  energy [eV] =\n" << std::setw(20) << real[0].IonEk << "\n  moment0 mean =\n    ";
    for (k = 0; k < Moment2State::maxsize; k++)
        strm << std::scientific << std::setprecision(8) << std::setw(16) << moment0_env(k);
    strm << std::scientific << std::setprecision(8)
         << "\nmoment0 rms =\n    ";
    for (k = 0; k < Moment2State::maxsize; k++)
        strm << std::scientific << std::setprecision(8) << std::setw(16) << moment0_rms(k);
    strm << "\n  state =\n";
    for (j = 0; j < Moment2State::maxsize; j++) {
        strm << "    ";
        for (k = 0; k < Moment2State::maxsize; k++) {
            strm << std::scientific << std::setprecision(8) << std::setw(16) << moment1[0](j, k);
        }
        if (j < Moment2State::maxsize-1) strm << "\n";
    }
}

bool Moment2State::getArray(unsigned idx, ArrayInfo& Info) {
    unsigned I=0;
    if(idx==I++) {
        Info.name = "state";
        Info.ptr = &moment1_env(0,0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 2;
        Info.dim[0] = moment1_env.size1();
        Info.dim[1] = moment1_env.size2();
        return true;
    } else if(idx==I++) {
        Info.name = "moment0";
        Info.ptr = &moment0_env(0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim[0] = moment0_env.size();
        return true;
    } else if(idx==I++) {
        Info.name = "moment0_rms";
        Info.ptr = &moment0_rms(0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim[0] = moment0_rms.size();
        return true;
    } else if(idx==I++) {
        Info.name = "ref_IonZ";
        Info.ptr = &ref.IonZ;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_IonQ";
        Info.ptr = &ref.IonQ;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_IonEs";
        Info.ptr = &ref.IonEs;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_IonW";
        Info.ptr = &ref.IonW;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_gamma";
        Info.ptr = &ref.gamma;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_beta";
        Info.ptr = &ref.beta;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_bg";
        Info.ptr = &ref.bg;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_SampleIonK";
        Info.ptr = &ref.SampleIonK;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_phis";
        Info.ptr = &ref.phis;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "ref_IonEk";
        Info.ptr = &ref.IonEk;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_IonZ";
        Info.ptr = &real[0].IonZ;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_IonEs";
        Info.ptr = &real[0].IonEs;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_IonW";
        Info.ptr = &real[0].IonW;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_gamma";
        Info.ptr = &real[0].gamma;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_beta";
        Info.ptr = &real[0].beta;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_bg";
        Info.ptr = &real[0].bg;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_SampleIonK";
        Info.ptr = &real[0].SampleIonK;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_phis";
        Info.ptr = &real[0].phis;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    } else if(idx==I++) {
        Info.name = "real_IonEk";
        Info.ptr = &real[0].IonEk;
        Info.type = ArrayInfo::Double;
        Info.ndim = 0;
        return true;
    }
    return StateBase::getArray(idx-I, Info);
}

Moment2ElementBase::Moment2ElementBase(const Config& c)
    :ElementVoid(c)
    ,misalign(boost::numeric::ublas::identity_matrix<double>(state_t::maxsize))
    ,misalign_inv(boost::numeric::ublas::identity_matrix<double>(state_t::maxsize))
    ,scratch(state_t::maxsize, state_t::maxsize)
{
}

Moment2ElementBase::~Moment2ElementBase() {}

void Moment2ElementBase::assign(const ElementVoid *other)
{
    const Moment2ElementBase *O = static_cast<const Moment2ElementBase*>(other);
    last_Kenergy_in = O->last_Kenergy_in;
    last_Kenergy_out = O->last_Kenergy_out;
    transfer = O->transfer;
    transfer_raw = O->transfer_raw;
    misalign = O->misalign;
    misalign_inv = O->misalign_inv;
    ElementVoid::assign(other);
}

void Moment2ElementBase::show(std::ostream& strm, int level) const
{
    using namespace boost::numeric::ublas;
    ElementVoid::show(strm, level);
    /*
    strm<<"Length "<<length<<"\n"
          "Transfer: "<<transfer<<"\n"
          "Transfer Raw: "<<transfer_raw<<"\n"
          "Mis-align: "<<misalign<<"\n";
          */
}

void Moment2ElementBase::get_misalign(const state_t &ST, const Particle &real, value_t &M, value_t &IM) const
{
    state_t::matrix_t R, R_inv,
              scl     = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize),
              scl_inv = scl,
              T       = scl,
              T_inv   = scl;

    double dx    = conf().get<double>("dx",    0e0)*MtoMM,
           dy    = conf().get<double>("dy",    0e0)*MtoMM,
           pitch = conf().get<double>("pitch", 0e0),
           yaw   = conf().get<double>("yaw",   0e0),
           tilt  = conf().get<double>("tilt",  0e0);

    scl(state_t::PS_S, state_t::PS_S)   /= -real.SampleIonK;
    scl(state_t::PS_PS, state_t::PS_PS) /= sqr(real.beta)*real.gamma*ST.ref.IonEs/MeVtoeV;

    inverse(scl_inv, scl);

    // Translate to center of element.
    T(state_t::PS_S,  6) = -length/2e0*MtoMM;
    T(state_t::PS_PS, 6) = 1e0;
    inverse(T_inv, T);

    RotMat(dx, dy, pitch, yaw, tilt, R);

    M = prod(T, scl);
    M = prod(R, M);
    M = prod(T_inv, M);
    M = prod(scl_inv, M);

    // J.B. Bug in TLM: should be inverse.
    RotMat(-dx, -dy, -pitch, -yaw, -tilt, R_inv);

    // Translate to center of element.
    T(state_t::PS_S,  6) = length/2e0*MtoMM;
    T(state_t::PS_PS, 6) = 1e0;
    inverse(T_inv, T);

    IM = prod(T, scl);
    IM = prod(R_inv, IM);
    IM = prod(T_inv, IM);
    IM = prod(scl_inv, IM);
}

void Moment2ElementBase::advance(StateBase& s)
{
    state_t&  ST = static_cast<state_t&>(s);
    using namespace boost::numeric::ublas;

    // IonEk is Es + E_state; the latter is set by user.
    ST.recalc();

    if(!check_cache(ST))
    {
        // need to re-calculate energy dependent terms

        resize_cache(ST);

        recompute_matrix(ST); // updates transfer and last_Kenergy_out

        ST.recalc();
    }

    // recompute_matrix only called when ST.IonEk != last_Kenergy_in.
    // Matrix elements are scaled with particle energy.

    ST.pos += length;

    ST.ref.phis   += ST.ref.SampleIonK*length*MtoMM;

    for(size_t k=0; k<last_Kenergy_in.size(); k++) {
        ST.real[k].phis  += ST.real[k].SampleIonK*length*MtoMM;
        ST.real[k].IonEk  = last_Kenergy_out[k];

        ST.moment0[k] = prod(transfer[k], ST.moment0[k]);

        scratch  = prod(transfer[k], ST.moment1[k]);
        ST.moment1[k] = prod(scratch, trans(transfer[k]));
    }
}

bool Moment2ElementBase::check_cache(const state_t& ST) const
{
    if(last_Kenergy_in.size()!=ST.size()) return false; // different # of charge states

    for(size_t k=0; k<last_Kenergy_in.size(); k++) {
        if(last_Kenergy_in[k]!=ST.real[k].IonEk) return false;
    }
    return true;
}

bool Moment2ElementBase::resize_cache(const state_t& ST)
{
    last_Kenergy_in.resize(ST.real.size());
    last_Kenergy_out.resize(ST.real.size());
    transfer.resize(ST.real.size(), boost::numeric::ublas::identity_matrix<double>(state_t::maxsize));
    transfer_raw.resize(ST.real.size(), boost::numeric::ublas::identity_matrix<double>(state_t::maxsize));
    misalign.resize(ST.real.size(), boost::numeric::ublas::identity_matrix<double>(state_t::maxsize));
    misalign_inv.resize(ST.real.size(), boost::numeric::ublas::identity_matrix<double>(state_t::maxsize));
}

void Moment2ElementBase::recompute_matrix(state_t& ST)
{
    // Default, for passive elements.

    for(size_t k=0; k<last_Kenergy_in.size(); k++) {
        transfer[k] = transfer_raw[k];

        last_Kenergy_in[k] = last_Kenergy_out[k] = ST.real[k].IonEk; // no energy gain
    }
}

namespace {

struct ElementSource : public Moment2ElementBase
{
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementSource(const Config& c): base_t(c), istate(c) {}

    virtual void advance(StateBase& s)
    {
        state_t& ST = static_cast<state_t&>(s);
        // Replace state with our initial values
        ST.assign(istate);
    }

    virtual void show(std::ostream& strm, int level) const
    {
        ElementVoid::show(strm, level);
        strm<<"Initial: "<<istate.moment0_env<<"\n";
    }

    state_t istate;
    // note that 'transfer' is not used by this element type

    virtual ~ElementSource() {}

    virtual const char* type_name() const {return "source";}
};

struct ElementMark : public Moment2ElementBase
{
    // Transport (identity) matrix for Marker.
    typedef Moment2ElementBase     base_t;
    typedef typename base_t::state_t state_t;

    ElementMark(const Config& c): base_t(c) {length = 0e0;}
    virtual ~ElementMark() {}
    virtual const char* type_name() const {return "marker";}
};

struct ElementBPM : public Moment2ElementBase
{
    // Transport (identity) matrix for BPM.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    Particle state;

    ElementBPM(const Config& c): base_t(c) {length = 0e0;}
    virtual ~ElementBPM() {}
    virtual const char* type_name() const {return "bpm";}
};

struct ElementDrift : public Moment2ElementBase
{
    // Transport matrix for Drift.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementDrift(const Config& c) : base_t(c) {}
    virtual ~ElementDrift() {}
    virtual const char* type_name() const {return "drift";}

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.

        const double L = length*MtoMM; // Convert from [m] to [mm].

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
            transfer[i](state_t::PS_X, state_t::PS_PX) = L;
            transfer[i](state_t::PS_Y, state_t::PS_PY) = L;
            transfer[i](state_t::PS_S, state_t::PS_PS) =
                -2e0*M_PI/(SampleLambda*ST.real[i].IonEs/MeVtoeV*cube(ST.real[i].bg))*L;

            last_Kenergy_in[i] = last_Kenergy_out[i] = ST.real[i].IonEk; // no energy gain
        }
    }
};

struct ElementOrbTrim : public Moment2ElementBase
{
    // Transport matrix for Orbit Trim.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementOrbTrim(const Config& c) : base_t(c) {length = 0e0;}
    virtual ~ElementOrbTrim() {}
    virtual const char* type_name() const {return "orbtrim";}

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.
        double theta_x = conf().get<double>("theta_x", 0e0),
               theta_y = conf().get<double>("theta_y", 0e0);

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
            transfer[i](state_t::PS_PX, 6) = theta_x*ST.ref.IonZ/ST.real[i].IonZ;
            transfer[i](state_t::PS_PY, 6) = theta_y*ST.ref.IonZ/ST.real[i].IonZ;

            last_Kenergy_in[i] = last_Kenergy_out[i] = ST.real[i].IonEk; // no energy gain

            get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            noalias(scratch)  = prod(transfer[i], misalign[i]);
            noalias(transfer) = prod(misalign_inv[i], scratch);
        }
    }
};

struct ElementSBend : public Moment2ElementBase
{
    // Transport matrix for Gradient Sector Bend; with edge focusing (cylindrical coordinates).
    // Note, TLM only includes energy offset for the orbit; not the transport matrix.

    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    unsigned HdipoleFitMode;

    ElementSBend(const Config& c) : base_t(c), HdipoleFitMode(0) {

        std::istringstream strm(c.get<std::string>("HdipoleFitMode", "1"));
        strm>>HdipoleFitMode;
        if(!strm.eof() && strm.fail())
            throw std::runtime_error("HdipoleFitMode must be an integer");
    }
    virtual ~ElementSBend() {}
    virtual const char* type_name() const {return "sbend";}

    virtual void advance(StateBase& s)
    {
        double    di_bg, Ek00, beta00, gamma00, IonK_Bend;
        state_t&  ST = static_cast<state_t&>(s);
        using namespace boost::numeric::ublas;

        // IonEk is Es + E_state; the latter is set by user.
        ST.recalc();

        if(!check_cache(ST)) {
            // need to re-calculate energy dependent terms
            resize_cache(ST);

            recompute_matrix(ST); // updates transfer and last_Kenergy_out

            ST.recalc();
        }

        // recompute_matrix only called when ST.IonEk != last_Kenergy_in.
        // Matrix elements are scaled with particle energy.

        ST.pos += length;

        ST.ref.phis += ST.ref.SampleIonK*length*MtoMM;

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            double phis_temp = ST.moment0[i][state_t::PS_S];

            ST.moment0[i] = prod(transfer[i], ST.moment0[i]);

            double dphis_temp = ST.moment0[i][state_t::PS_S] - phis_temp;

            if (HdipoleFitMode != 1) {
                di_bg     = conf().get<double>("bg");
                // Dipole reference energy.
                Ek00      = (sqrt(sqr(di_bg)+1e0)-1e0)*ST.ref.IonEs;
                gamma00   = (Ek00+ST.ref.IonEs)/ST.ref.IonEs;
                beta00    = sqrt(1e0-1e0/sqr(gamma00));
                IonK_Bend = 2e0*M_PI/(beta00*SampleLambda);

                // J.B.: this is odd.
    //            ST.real.phis  += IonK_Bend*length*MtoMM + dphis_temp;
                ST.real[i].phis  += ST.real[i].SampleIonK*length*MtoMM + dphis_temp;
            } else
                ST.real[i].phis  += ST.real[i].SampleIonK*length*MtoMM + dphis_temp;

            ST.real[i].IonEk  = last_Kenergy_out[i];

            noalias(scratch)  = prod(transfer[i], ST.moment1[i]);
            noalias(ST.moment1[i]) = prod(scratch, trans(transfer[i]));
        }
    }

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.

        double L     = conf().get<double>("L")*MtoMM,
               phi   = conf().get<double>("phi")*M_PI/180e0,
               phi1  = conf().get<double>("phi1")*M_PI/180e0,
               phi2  = conf().get<double>("phi2")*M_PI/180e0,
               K     = conf().get<double>("K", 0e0)/sqr(MtoMM);

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            double qmrel = (ST.real[i].IonZ-ST.ref.IonZ)/ST.ref.IonZ;

            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

            if (HdipoleFitMode != 1) {
                double dip_bg    = conf().get<double>("bg"),
                       // Dipole reference energy.
                       dip_Ek    = (sqrt(sqr(dip_bg)+1e0)-1e0)*ST.ref.IonEs,
                       dip_gamma = (dip_Ek+ST.ref.IonEs)/ST.ref.IonEs,
                       dip_beta  = sqrt(1e0-1e0/sqr(dip_gamma)),
                       d         = (ST.ref.gamma-dip_gamma)/(sqr(dip_beta)*dip_gamma) - qmrel,
                       dip_IonK  = 2e0*M_PI/(dip_beta*SampleLambda);

                GetSBendMatrix(L, phi, phi1, phi2, K, ST.ref.IonEs, ST.ref.gamma, qmrel,
                               dip_beta, dip_gamma, d, dip_IonK, transfer[i]);
            } else
                GetSBendMatrix(L, phi, phi1, phi2, K, ST.ref.IonEs, ST.ref.gamma, qmrel,
                               ST.ref.beta, ST.ref.gamma, - qmrel, ST.ref.SampleIonK, transfer[i]);

            get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            noalias(scratch)     = prod(transfer[i], misalign[i]);
            noalias(transfer[i]) = prod(misalign_inv[i], scratch);

            last_Kenergy_in[i] = last_Kenergy_out[i] = ST.real[i].IonEk; // no energy gain
        }
    }
};

struct ElementQuad : public Moment2ElementBase
{
    // Transport matrix for Quadrupole; K = B2/Brho.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementQuad(const Config& c) : base_t(c) {}
    virtual ~ElementQuad() {}
    virtual const char* type_name() const {return "quadrupole";}

    virtual void recompute_matrix(state_t& ST)
    {
        const double B2= conf().get<double>("B2"),
                     L = conf().get<double>("L")*MtoMM;

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            // Re-initialize transport matrix.
            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

            double Brho = ST.real[i].beta*(ST.real[i].IonEk+ST.real[i].IonEs)/(C0*ST.real[i].IonZ),
                   K = B2/Brho/sqr(MtoMM);

            // Horizontal plane.
            GetQuadMatrix(L,  K, (unsigned)state_t::PS_X, transfer[i]);
            // Vertical plane.
            GetQuadMatrix(L, -K, (unsigned)state_t::PS_Y, transfer[i]);
            // Longitudinal plane.
    //        transfer_raw(state_t::PS_S, state_t::PS_S) = L;

            transfer[i](state_t::PS_S, state_t::PS_PS) =
                    -2e0*M_PI/(SampleLambda*ST.real[i].IonEs/MeVtoeV*cube(ST.real[i].bg))*L;

            get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            noalias(scratch)     = prod(transfer[i], misalign[i]);
            noalias(transfer[i]) = prod(misalign_inv[i], scratch);

            last_Kenergy_in[i] = last_Kenergy_out[i] = ST.real[i].IonEk; // no energy gain
        }
    }
};

struct ElementSolenoid : public Moment2ElementBase
{
    // Transport (identity) matrix for a Solenoid; K = B/(2 Brho).
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementSolenoid(const Config& c) : base_t(c) {}
    virtual ~ElementSolenoid() {}
    virtual const char* type_name() const {return "solenoid";}

    virtual void recompute_matrix(state_t& ST)
    {
        const double B = conf().get<double>("B"),
                     L = conf().get<double>("L")*MtoMM;      // Convert from [m] to [mm].

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            // Re-initialize transport matrix.
            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

            double Brho = ST.real[i].beta*(ST.real[i].IonEk+ST.real[i].IonEs)/(C0*ST.real[i].IonZ),
                   K    = B/(2e0*Brho)/MtoMM;

            GetSolMatrix(L, K, transfer[i]);

            transfer[i](state_t::PS_S, state_t::PS_PS) =
                    -2e0*M_PI/(SampleLambda*ST.real[i].IonEs/MeVtoeV*cube(ST.real[i].bg))*L;

            get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            noalias(scratch)     = prod(transfer[i], misalign[i]);
            noalias(transfer[i]) = prod(misalign_inv[i], scratch);

            last_Kenergy_in[i] = last_Kenergy_out[i] = ST.real[i].IonEk; // no energy gain
        }
    }
};

struct ElementEDipole : public Moment2ElementBase
{
    // Transport matrix for Electrostatic Dipole with edge focusing.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementEDipole(const Config& c) : base_t(c) {}
    virtual ~ElementEDipole() {}
    virtual const char* type_name() const {return "edipole";}

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.

        //value_mat R;

        bool   ver         = conf().get<double>("ver") == 1.0;
        double L           = conf().get<double>("L")*MtoMM,
               phi         = conf().get<double>("phi")*M_PI/180e0,
               // fit to TLM unit.            
               fringe_x    = conf().get<double>("fringe_x", 0e0)/MtoMM,
               fringe_y    = conf().get<double>("fringe_y", 0e0)/MtoMM,
               kappa       = conf().get<double>("asym_fac", 0e0),
               // spher: cylindrical - 0, spherical - 1.
               spher       = conf().get<double>("spher"),
               rho         = L/phi,
               // magnetic - 0, electrostatic - 1.
               h           = 1e0,
               Ky          = spher/sqr(rho),
               dip_beta    = conf().get<double>("beta"),
               dip_gamma   = 1e0/sqrt(1e0-sqr(dip_beta));

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            double eta0        = (ST.real[i].gamma-1e0)/2e0,
                   Kx          = (1e0-spher+sqr(1e0+2e0*eta0))/sqr(rho),
                   delta_KZ    = ST.ref.IonZ/ST.real[i].IonZ - 1e0,
                   SampleIonK  = 2e0*M_PI/(ST.real[i].beta*SampleLambda);

            transfer_raw[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);


            GetEBendMatrix(L, phi, fringe_x, fringe_y, kappa, Kx, Ky, ST.ref.IonEs, ST.ref.beta, ST.real[i].gamma,
                           eta0, h, dip_beta, dip_gamma, delta_KZ, SampleIonK, transfer_raw[i]);

            if (ver) {
                // Rotate transport matrix by 90 degrees.
                value_t
                R = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
                R(state_t::PS_X,  state_t::PS_X)   =  0e0;
                R(state_t::PS_PX, state_t::PS_PX)  =  0e0;
                R(state_t::PS_Y,  state_t::PS_Y)   =  0e0;
                R(state_t::PS_PY, state_t::PS_PY)  =  0e0;
                R(state_t::PS_X,  state_t::PS_Y)   = -1e0;
                R(state_t::PS_PX, state_t::PS_PY)  = -1e0;
                R(state_t::PS_Y,  state_t::PS_X)   =  1e0;
                R(state_t::PS_PY,  state_t::PS_PX) =  1e0;

                noalias(scratch)  = prod(transfer[i], misalign[i]);
                noalias(transfer[i]) = prod(misalign_inv[i], scratch);
                //TODO: no-op code?  results are unconditionally overwritten
            }

            transfer[i] = transfer_raw[i];

            get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            noalias(scratch)     = prod(transfer[i], misalign[i]);
            noalias(transfer[i]) = prod(misalign_inv[i], scratch);

            last_Kenergy_in[i] = last_Kenergy_out[i] = ST.real[i].IonEk; // no energy gain
        }
    }
   virtual void assign(const ElementVoid *other)
   {
        throw std::logic_error("assign/reconfigure() not implemented for rf cavity");
        // TODO: can't copy 'inf'.  Need to parse once in ctor.
//        const ElementRFCavity *O = static_cast<const ElementRFCavity*>(other);
//        Moment2ElementBase::assign(O);
//        CavData = O->CavData;
//        CavTLMLineTab = O->CavTLMLineTab;
//        phi_ref = O->phi_ref;
   }

};

struct ElementEQuad : public Moment2ElementBase
{
    // Transport matrix for Electrostatic Quadrupole.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementEQuad(const Config& c) : base_t(c) {}
    virtual ~ElementEQuad() {}
    virtual const char* type_name() const {return "equad";}

    virtual void recompute_matrix(state_t& ST)
    {
        const double V0   = conf().get<double>("V"),
                     R    = conf().get<double>("radius"),
                     L    = conf().get<double>("L")*MtoMM;

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            // Re-initialize transport matrix.
            // V0 [V] electrode voltage and R [m] electrode half-distance.
            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

            double Brho = ST.real[i].beta*(ST.real[i].IonEk+ST.real[i].IonEs)/(C0*ST.real[i].IonZ),
                   K    = 2e0*V0/(C0*ST.real[i].beta*sqr(R))/Brho/sqr(MtoMM);

            // Horizontal plane.
            GetQuadMatrix(L,  K, (unsigned)state_t::PS_X, transfer[i]);
            // Vertical plane.
            GetQuadMatrix(L, -K, (unsigned)state_t::PS_Y, transfer[i]);
            // Longitudinal plane.
            //        transfer_raw(state_t::PS_S, state_t::PS_S) = L;

            transfer[i](state_t::PS_S, state_t::PS_PS) =
                    -2e0*M_PI/(SampleLambda*ST.real[i].IonEs/MeVtoeV*cube(ST.real[i].bg))*L;

            get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            noalias(scratch)     = prod(transfer[i], misalign[i]);
            noalias(transfer[i]) = prod(misalign_inv[i], scratch);

            last_Kenergy_in[i] = last_Kenergy_out[i] = ST.real[i].IonEk; // no energy gain
        }
    }
};

struct ElementGeneric : public Moment2ElementBase
{
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    value_t proto;

    ElementGeneric(const Config& c)
        :base_t(c)
        ,proto(state_t::maxsize, state_t::maxsize)
    {
        load_storage(proto.data(), c, "transfer");
    }
    virtual ~ElementGeneric() {}

    virtual void recompute_matrix(state_t& ST)
    {
        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            transfer[i] = proto;

            last_Kenergy_in[i] = last_Kenergy_out[i] = ST.real[i].IonEk; // no energy gain
        }
    }

    virtual const char* type_name() const {return "generic";}
};

} // namespace

void registerMoment2()
{
    Machine::registerState<Moment2State>("MomentMatrix2");

    Machine::registerElement<ElementSource                 >("MomentMatrix2", "source");

    Machine::registerElement<ElementMark                   >("MomentMatrix2", "marker");

    Machine::registerElement<ElementBPM                    >("MomentMatrix2", "bpm");

    Machine::registerElement<ElementDrift                  >("MomentMatrix2", "drift");

    Machine::registerElement<ElementOrbTrim                >("MomentMatrix2", "orbtrim");

    Machine::registerElement<ElementSBend                  >("MomentMatrix2", "sbend");

    Machine::registerElement<ElementQuad                   >("MomentMatrix2", "quadrupole");

    Machine::registerElement<ElementSolenoid               >("MomentMatrix2", "solenoid");

    Machine::registerElement<ElementRFCavity               >("MomentMatrix2", "rfcavity");

    Machine::registerElement<ElementStripper               >("MomentMatrix2", "stripper");

    Machine::registerElement<ElementEDipole                >("MomentMatrix2", "edipole");

    Machine::registerElement<ElementEQuad                  >("MomentMatrix2", "equad");

    Machine::registerElement<ElementGeneric                >("MomentMatrix2", "generic");
}
