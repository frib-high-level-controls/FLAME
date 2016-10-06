
#include <fstream>

#include <limits>

#include <boost/lexical_cast.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "flame/constants.h"

#include "flame/moment.h"
#include "flame/moment_sup.h"
#include "flame/rf_cavity.h"
#include "flame/chg_stripper.h"

#include "flame/h5loader.h"

#define sqr(x)  ((x)*(x))
#define cube(x) ((x)*(x)*(x))

namespace ub = boost::numeric::ublas;

namespace {

// ARR should be an array-like object (std::vector or or ublas vector or matrix storage)
// fill 'to' using config.get<>(name)
// 'T' selects throw error (true) or return boolean (false)
template<typename ARR>
bool load_storage(ARR& to, const Config& conf, const std::string& name, bool T=true)
{
    try{
        const std::vector<double>& val(conf.get<std::vector<double> >(name));
        if(to.size()!=val.size()) {
            throw std::invalid_argument(SB()<<"Array "<<name<<" must have "<<to.size()<<" elements, not "<<val.size());
        }
        std::copy(val.begin(), val.end(), to.begin());

        return true;
    }catch(std::exception&){
        if(T)
            throw std::invalid_argument(name+" not defined or has wrong type (must be vector)");
        else
            return false;
        // default to identity
    }
}

} // namespace

std::ostream& operator<<(std::ostream& strm, const Particle& P)
{
    strm
      <<"IonZ="<<std::scientific << std::setprecision(10)<<P.IonZ
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

MomentState::MomentState(const Config& c)
    :StateBase(c)
    ,ref()
    ,real()
    ,moment0_env(maxsize, 0e0)
    ,moment0_rms(maxsize, 0e0)
    ,moment1_env(boost::numeric::ublas::identity_matrix<double>(maxsize))
{
    sim_mode = (int)c.get<double>("sim_mode", 0e0);

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
        throw std::invalid_argument("MomentState: must define IonChargeStates and NCharge when cstate is set");
    }

    if(have_ics) {
        real.resize(ics.size());
        moment0.resize(ics.size());
        moment1.resize(ics.size());

        for(size_t i=0; i<ics.size(); i++) {
            std::string num(boost::lexical_cast<std::string>(icstate+i));

            moment0[i].resize(maxsize);
            moment1[i].resize(maxsize, maxsize);
            moment1[i] = boost::numeric::ublas::identity_matrix<double>(maxsize);

            load_storage(moment0[i].data(), c, vectorname+num);

            if ((int)sim_mode != 2){
                load_storage(moment1[i].data(), c, matrixname+num);
            }

            real[i] = ref;

            real[i].IonZ = ics[i];
            real[i].IonQ = nchg[i];

            real[i].phis      = moment0[i][PS_S];
            real[i].IonEk    += moment0[i][PS_PS]*MeVtoeV;

            real[i].recalc();
        }
    } else {
        real.resize(1); // hack, ensure at least one element so getArray() can return some pointer
        real[0] = ref;

        moment0.resize(1);
        moment1.resize(1);
        moment0[0].resize(maxsize);
        moment1[0].resize(maxsize, maxsize);

        load_storage(moment0[0].data(), c, vectorname, false);
        load_storage(moment1[0].data(), c, matrixname, false);
    }

    calc_rms();
}

MomentState::~MomentState() {}

void MomentState::calc_rms()
{
    assert(real.size()>0);
    assert(moment0_env.size()==maxsize);
    assert(moment0_rms.size()==maxsize);
    assert(moment1_env.size1()==maxsize);
    assert(moment1_env.size2()==maxsize);

    double totQ = 0.0;
    for(size_t n=0; n<real.size(); n++) {
        totQ += real[n].IonQ;
        if(n==0)
            noalias(moment0_env)  = moment0[n]*real[n].IonQ;
        else
            noalias(moment0_env) += moment0[n]*real[n].IonQ;
    }
    moment0_env /= totQ;

    // Zero orbit terms.
    moment1_env = ub::zero_matrix<double>(state_t::maxsize);
    boost::numeric::ublas::slice S(0, 1, 6);
    for(size_t n=0; n<real.size(); n++) {
        const double Q = real[n].IonQ;
        state_t::vector_t m0diff(moment0[n]-moment0_env);

        ub::project(moment1_env, S, S) += ub::project(Q*(moment1[n]+ub::outer_prod(m0diff, m0diff)), S, S);
    }
    moment1_env /= totQ;

    for(size_t j=0; j<maxsize; j++) {
        moment0_rms[j] = sqrt(moment1_env(j,j));
    }
}

MomentState::MomentState(const MomentState& o, clone_tag t)
    :StateBase(o, t)
    ,ref(o.ref)
    ,real(o.real)
    ,moment0(o.moment0)
    ,moment1(o.moment1)
    ,moment0_env(o.moment0_env)
    ,moment0_rms(o.moment0_rms)
    ,moment1_env(o.moment1_env)
    ,sim_mode(o.sim_mode)
{}

void MomentState::assign(const StateBase& other)
{
    const MomentState *O = dynamic_cast<const MomentState*>(&other);
    if(!O)
        throw std::invalid_argument("Can't assign State: incompatible types");
    ref     = O->ref;
    real    = O->real;
    moment0 = O->moment0;
    moment1 = O->moment1;
    moment0_env = O->moment0_env;
    moment0_rms = O->moment0_rms;
    moment1_env = O->moment1_env;
    sim_mode = O->sim_mode;
    StateBase::assign(other);
}

void MomentState::show(std::ostream& strm, int level) const
{
    if(real.empty()) {
        strm<<"State: empty";
        return;
    }

    if(level<=0) {
        strm<<"State: moment0 mean="<<moment0_env;
    }
    if(level>=1) {
        strm << std::scientific << std::setprecision(8)
             << "\nState:\n  energy [eV] =\n" << std::setw(20) << real[0].IonEk << "\n  moment0 mean =\n    ";
        for (size_t k = 0; k < MomentState::maxsize; k++)
            strm << std::scientific << std::setprecision(10) << std::setw(18) << moment0_env(k) << ",";
        strm << std::scientific << std::setprecision(10)
             << "\n  moment0 rms =\n    ";
        for (size_t k = 0; k < MomentState::maxsize; k++)
            strm << std::scientific << std::setprecision(10) << std::setw(18) << moment0_rms(k) << ",";
        strm << "\n  moment1 mean =\n";
        for (size_t j = 0; j < MomentState::maxsize; j++) {
            strm << "    ";
            for (size_t k = 0; k < MomentState::maxsize; k++) {
                strm << std::scientific << std::setprecision(10) << std::setw(18) << moment1_env(j, k) << ",";
            }
            if (j < MomentState::maxsize-1) strm << "\n";
        }
    }
    if(level>=2) {
        strm<< "\n  Reference state:\n    "<<ref<<"\n";
        for(size_t k=0; k<size(); k++) {
            strm<<"  Charge state "<<k<<"\n"
                  "    "<<real[k]<<"\n"
                  "    moment0 "<<moment0[k]<<"\n";
            if(level>=3)
                strm<<
                  "    moment1 "<<moment1[k]<<"\n";
        }
    }
}

bool MomentState::getArray(unsigned idx, ArrayInfo& Info) {
    unsigned I=0;
    if(idx==I++) {
        Info.name = "moment1_env";
        Info.ptr = &moment1_env(0,0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 2;
        Info.dim[0] = moment1_env.size1();
        Info.dim[1] = moment1_env.size2();
        Info.stride[0] = sizeof(double)*moment1_env.size2();
        Info.stride[1] = sizeof(double);
        return true;
    } else if(idx==I++) {
        /* Slight evilness here
         * moment0 is vector of ublas::matrix
         * We assume ublax::matrix uses storage bounded_array<>, and that this storage
         * is really a C array, which means that everything is part of one big allocation.
         * Further we assume that all entries in the vector have the same shape.
         * If this isn't the case, then SIGSEGV here we come...
         */
        static_assert(sizeof(moment1[0])>=sizeof(double)*maxsize*maxsize,
                      "storage assumption violated");
        Info.name = "moment1";
        Info.ptr = &moment1[0](0,0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 3;
        Info.dim[0] = moment1[0].size1();
        Info.dim[1] = moment1[0].size2();
        Info.dim[2] = moment1.size();
        Info.stride[0] = sizeof(double)*moment1_env.size2();
        Info.stride[1] = sizeof(double);
        Info.stride[2] = sizeof(moment1[0]);
        return true;
    } else if(idx==I++) {
        Info.name = "moment0_env";
        Info.ptr = &moment0_env(0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim[0] = moment0_env.size();
        Info.stride[0] = sizeof(double);
        return true;
    } else if(idx==I++) {
        Info.name = "moment0_rms";
        Info.ptr = &moment0_rms(0);
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim[0] = moment0_rms.size();
        Info.stride[0] = sizeof(double);
        return true;
    } else if(idx==I++) {
        // more evilness here, see above
        static_assert(sizeof(moment0[0])>=sizeof(double)*maxsize,
                "storage assumption violated");
        Info.name = "moment0";
        Info.ptr = &moment0[0][0];
        Info.type = ArrayInfo::Double;
        Info.ndim = 2;
        Info.dim[0] = moment0[0].size();
        Info.dim[1] = moment0.size();
        Info.stride[0] = sizeof(double);
        Info.stride[1] = sizeof(moment0[0]);
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
        Info.name = "IonZ";
        Info.ptr  = &real[0].IonZ;
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim   [0] = real.size();
        Info.stride[0] = sizeof(real[0]);
        // Note: this array is discontigious as we reference a single member from a Particle[]
        return true;
    } else if(idx==I++) {
        Info.name = "IonEs";
        Info.ptr  = &real[0].IonEs;
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim   [0] = real.size();
        Info.stride[0] = sizeof(real[0]);
        return true;
    } else if(idx==I++) {
        Info.name = "IonW";
        Info.ptr  = &real[0].IonW;
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim   [0] = real.size();
        Info.stride[0] = sizeof(real[0]);
        return true;
    } else if(idx==I++) {
        Info.name = "gamma";
        Info.ptr  = &real[0].gamma;
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim   [0] = real.size();
        Info.stride[0] = sizeof(real[0]);
        return true;
    } else if(idx==I++) {
        Info.name = "beta";
        Info.ptr  = &real[0].beta;
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim   [0] = real.size();
        Info.stride[0] = sizeof(real[0]);
        return true;
    } else if(idx==I++) {
        Info.name = "bg";
        Info.ptr  = &real[0].bg;
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim   [0] = real.size();
        Info.stride[0] = sizeof(real[0]);
        return true;
    } else if(idx==I++) {
        Info.name = "SampleIonK";
        Info.ptr  = &real[0].SampleIonK;
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim   [0] = real.size();
        Info.stride[0] = sizeof(real[0]);
        return true;
    } else if(idx==I++) {
        Info.name = "phis";
        Info.ptr  = &real[0].phis;
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim   [0] = real.size();
        Info.stride[0] = sizeof(real[0]);
        return true;
    } else if(idx==I++) {
        Info.name = "IonEk";
        Info.ptr  = &real[0].IonEk;
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim   [0] = real.size();
        Info.stride[0] = sizeof(real[0]);
        return true;
    } else if(idx==I++) {
        Info.name = "IonQ";
        Info.ptr  = &real[0].IonQ;
        Info.type = ArrayInfo::Double;
        Info.ndim = 1;
        Info.dim   [0] = real.size();
        Info.stride[0] = sizeof(real[0]);
        // Note: this array is discontigious as we reference a single member from a Particle[]
        return true;
    } else if(idx==I++) {
        Info.name = "sim_mode";
        Info.ptr = &sim_mode;
        Info.type = ArrayInfo::Sizet;
        Info.ndim = 0;
        return true;
    }
    return StateBase::getArray(idx-I, Info);
}

MomentElementBase::MomentElementBase(const Config& c)
    :ElementVoid(c)
    ,dx   (c.get<double>("dx",    0e0)*MtoMM)
    ,dy   (c.get<double>("dy",    0e0)*MtoMM)
    ,pitch(c.get<double>("pitch", 0e0))
    ,yaw  (c.get<double>("yaw",   0e0))
    ,roll (c.get<double>("roll",  0e0))
    ,skipcache(c.get<double>("skipcache", 0.0)!=0.0)
    ,strict3Dvelocity(c.get<double>("strict3Dvelocity", 0.0)!=0.0)
    ,scratch(state_t::maxsize, state_t::maxsize)
{
}

MomentElementBase::~MomentElementBase() {}

void MomentElementBase::assign(const ElementVoid *other)
{
    const MomentElementBase *O = static_cast<const MomentElementBase*>(other);
    last_real_in = O->last_real_in;
    last_real_out= O->last_real_out;
    last_ref_in  = O->last_ref_in;
    last_ref_out = O->last_ref_out;
    transfer = O->transfer;
    misalign = O->misalign;
    misalign_inv = O->misalign_inv;
    dx = O->dx;
    dy = O->dy;
    pitch = O->pitch;
    yaw   = O->yaw;
    roll  = O->roll;
    skipcache = O->skipcache;
    strict3Dvelocity = O->strict3Dvelocity;
    ElementVoid::assign(other);
}

void MomentElementBase::show(std::ostream& strm, int level) const
{
    using namespace boost::numeric::ublas;
    ElementVoid::show(strm, level);
    /*
    strm<<"Length "<<length<<"\n"
          "Transfer: "<<transfer<<"\n"
          "Mis-align: "<<misalign<<"\n";
          */
}

void MomentElementBase::get_misalign(const state_t &ST, const Particle &real, value_t &M, value_t &IM) const
{
    state_t::matrix_t R,
              scl     = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize),
              scl_inv = scl,
              T       = scl,
              T_inv   = scl,
              R_inv   = scl;

    scl(state_t::PS_S, state_t::PS_S)   /= -real.SampleIonK;
    scl(state_t::PS_PS, state_t::PS_PS) /= sqr(real.beta)*real.gamma*ST.ref.IonEs/MeVtoeV;

    inverse(scl_inv, scl);

    // Translate to center of element.
    T(state_t::PS_S,  6) = -length/2e0*MtoMM;
    T(state_t::PS_PS, 6) = 1e0;
    inverse(T_inv, T);

    RotMat(dx, dy, pitch, yaw, roll, R);

    M = prod(T, scl);
    M = prod(R, M);
    M = prod(T_inv, M);
    M = prod(scl_inv, M);

    inverse(R_inv, R);

    // Translate to center of element.
    T(state_t::PS_S,  6) = length/2e0*MtoMM;
    T(state_t::PS_PS, 6) = 1e0;
    inverse(T_inv, T);

    IM = prod(T, scl);
    IM = prod(R_inv, IM);
    IM = prod(T_inv, IM);
    IM = prod(scl_inv, IM);
}

void MomentElementBase::advance(StateBase& s)
{
    state_t&  ST = static_cast<state_t&>(s);
    using namespace boost::numeric::ublas;

    // IonEk is Es + E_state; the latter is set by user.
    ST.recalc();

    if ((int)ST.sim_mode == 1)
    {
        // limit to longitudinal run
    
        ST.pos += length;
        ST.ref.phis   += ST.ref.SampleIonK*length*MtoMM;

    } else {

        if(!check_cache(ST)) 
        {
            // need to re-calculate energy dependent terms

            last_ref_in = ST.ref;
            last_real_in = ST.real;
            resize_cache(ST);

            recompute_matrix(ST); // updates transfer and last_Kenergy_out

            ST.recalc();

            for(size_t k=0; k<last_real_in.size(); k++)
                ST.real[k].phis  += ST.real[k].SampleIonK*length*MtoMM;
            ST.ref.phis   += ST.ref.SampleIonK*length*MtoMM;

            last_ref_out = ST.ref;
            last_real_out = ST.real;
        } else {
            ST.ref = last_ref_out;
            assert(last_real_out.size()==ST.real.size()); // should be true if check_cache() -> true
            std::copy(last_real_out.begin(),
                      last_real_out.end(),
                      ST.real.begin());
        }

        ST.pos += length;

        for(size_t k=0; k<last_real_in.size(); k++) {

            if (strict3Dvelocity) {
                double sclz3d = sqr(tan(ST.moment0[k][state_t::PS_PX])) + sqr(tan(ST.moment0[k][state_t::PS_PY])) + 1e0;
                transfer[k](state_t::PS_S, 6) =
                        -2e0*M_PI/(SampleLambda*ST.real[k].IonEs/MeVtoeV*cube(ST.real[k].bg))*(1e0/sclz3d-1e0)*ST.real[k].IonEk/MeVtoeV*length*MtoMM;
            }

            ST.moment0[k] = prod(transfer[k], ST.moment0[k]);

            if (ST.sim_mode == 2){
                ST.moment1[k] = prod(transfer[k], ST.moment1[k]); 
            } else {
                scratch  = prod(transfer[k], ST.moment1[k]);
                ST.moment1[k] = prod(scratch, trans(transfer[k]));
            }
        }

        ST.calc_rms();
    }
}

bool MomentElementBase::check_cache(const state_t& ST) const
{
    return !skipcache && !strict3Dvelocity
            && last_real_in.size()==ST.size()
            && last_ref_in==ST.ref
            && std::equal(last_real_in.begin(),
                          last_real_in.end(),
                          ST.real.begin());
}

void MomentElementBase::resize_cache(const state_t& ST)
{
    transfer.resize(ST.real.size(), boost::numeric::ublas::identity_matrix<double>(state_t::maxsize));
    misalign.resize(ST.real.size(), boost::numeric::ublas::identity_matrix<double>(state_t::maxsize));
    misalign_inv.resize(ST.real.size(), boost::numeric::ublas::identity_matrix<double>(state_t::maxsize));
}

void MomentElementBase::recompute_matrix(state_t& ST)
{
    // Default, initialize as no-op

    for(size_t k=0; k<last_real_in.size(); k++) {
        transfer[k] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
    }
}

namespace {

struct ElementSource : public MomentElementBase
{
    typedef ElementSource            self_t;
    typedef MomentElementBase       base_t;
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

    virtual void assign(const ElementVoid *other) {
        base_t::assign(other);
        const self_t* O=static_cast<const self_t*>(other);
        istate.assign(O->istate);
    }
};

struct ElementMark : public MomentElementBase
{
    // Transport (identity) matrix for Marker.
    typedef ElementMark            self_t;
    typedef MomentElementBase     base_t;
    typedef typename base_t::state_t state_t;

    ElementMark(const Config& c): base_t(c) {length = 0e0;}
    virtual ~ElementMark() {}
    virtual const char* type_name() const {return "marker";}

    virtual void assign(const ElementVoid *other) { base_t::assign(other); }
};

struct ElementBPM : public MomentElementBase
{
    // Transport (identity) matrix for BPM.
    typedef ElementBPM               self_t;
    typedef MomentElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementBPM(const Config& c): base_t(c) {length = 0e0;}
    virtual ~ElementBPM() {}
    virtual const char* type_name() const {return "bpm";}

    virtual void assign(const ElementVoid *other) { base_t::assign(other); }
};

struct ElementDrift : public MomentElementBase
{
    // Transport matrix for Drift.
    typedef ElementDrift             self_t;
    typedef MomentElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementDrift(const Config& c) : base_t(c) {}
    virtual ~ElementDrift() {}
    virtual const char* type_name() const {return "drift";}

    virtual void assign(const ElementVoid *other) { base_t::assign(other); }

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.

        const double L = length*MtoMM; // Convert from [m] to [mm].

        for(size_t i=0; i<last_real_in.size(); i++) {

            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
            transfer[i](state_t::PS_X, state_t::PS_PX) = L;
            transfer[i](state_t::PS_Y, state_t::PS_PY) = L;
            transfer[i](state_t::PS_S, state_t::PS_PS) = 
                    -2e0*M_PI/(SampleLambda*ST.real[i].IonEs/MeVtoeV*cube(ST.real[i].bg))*L;
        }
    }
};

struct ElementOrbTrim : public MomentElementBase
{
    // Transport matrix for Orbit Trim.
    typedef ElementOrbTrim           self_t;
    typedef MomentElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementOrbTrim(const Config& c) : base_t(c) {length = 0e0;}
    virtual ~ElementOrbTrim() {}
    virtual const char* type_name() const {return "orbtrim";}

    virtual void assign(const ElementVoid *other) { base_t::assign(other); }

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.
        double theta_x = conf().get<double>("theta_x", 0e0),
               theta_y = conf().get<double>("theta_y", 0e0);

        for(size_t i=0; i<last_real_in.size(); i++) {
            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
            transfer[i](state_t::PS_PX, 6) = theta_x*ST.real[i].IonZ/ST.ref.IonZ;
            transfer[i](state_t::PS_PY, 6) = theta_y*ST.real[i].IonZ/ST.ref.IonZ;

            get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            noalias(scratch)  = prod(transfer[i], misalign[i]);
            noalias(transfer[i]) = prod(misalign_inv[i], scratch);
        }
    }
};

struct ElementSBend : public MomentElementBase
{
    // Transport matrix for Gradient Sector Bend; with edge focusing (cylindrical coordinates).
    // Note, TLM only includes energy offset for the orbit; not the transport matrix.

    typedef ElementSBend             self_t;
    typedef MomentElementBase       base_t;
    typedef typename base_t::state_t state_t;

    unsigned HdipoleFitMode;
    double L, phi, phi1, phi2, K, Ftype;


    ElementSBend(const Config& c) : base_t(c), HdipoleFitMode(0) {

        std::istringstream strm(c.get<std::string>("HdipoleFitMode", "1"));
        strm>>HdipoleFitMode;
        if(!strm.eof() && strm.fail())
            throw std::runtime_error("HdipoleFitMode must be an integer");

        L     = conf().get<double>("L")*MtoMM;
        phi   = conf().get<double>("phi")*M_PI/180e0;
        phi1  = conf().get<double>("phi1")*M_PI/180e0;
        phi2  = conf().get<double>("phi2")*M_PI/180e0;
        K     = conf().get<double>("K", 0e0)/sqr(MtoMM);
        Ftype  = conf().get<double>("FringeType", 1e0);
        if ( ((float)(int)Ftype != Ftype) or (Ftype > 2e0 or Ftype < 0e0) )
            throw std::runtime_error("sbend: FringeType must be 0, 1, or 2");

    }
    virtual ~ElementSBend() {}
    virtual const char* type_name() const {return "sbend";}

    virtual void assign(const ElementVoid *other) {
        base_t::assign(other);
        const self_t* O=static_cast<const self_t*>(other);
        HdipoleFitMode = O->HdipoleFitMode;
        L = O->L;
        phi = O->phi;
        phi1 = O->phi1;
        phi2 = O->phi2;
        K   = O->K;
        Ftype = O->Ftype;
    }

    virtual void advance(StateBase& s)
    {
        state_t&  ST = static_cast<state_t&>(s);
        using namespace boost::numeric::ublas;

        // IonEk is Es + E_state; the latter is set by user.
        ST.recalc();

        if ((int)ST.sim_mode == 1)
        {
            // limit to longitudinal run
        
            ST.pos += length;
            ST.ref.phis   += ST.ref.SampleIonK*length*MtoMM;

        } else {


            if(!check_cache(ST)) {
                // need to re-calculate energy dependent terms
                last_ref_in = ST.ref;
                last_real_in = ST.real;
                resize_cache(ST);

                recompute_matrix(ST); // updates transfer and last_Kenergy_out

                ST.recalc();
                last_ref_out = ST.ref;
                last_real_out = ST.real;
            } else {
                ST.ref = last_ref_out;
                assert(last_real_out.size()==ST.real.size()); // should be true if check_cache() -> true
                std::copy(last_real_out.begin(),
                          last_real_out.end(),
                          ST.real.begin());
            }

            ST.pos += length;

            ST.ref.phis += ST.ref.SampleIonK*length*MtoMM;

            for(size_t i=0; i<last_real_in.size(); i++) {
                double phis_temp = ST.moment0[i][state_t::PS_S];

                /*
                if ((int)Ftype != 2) {
                    ST.moment0[i]          = prod(transfer[i], ST.moment0[i]);
                    noalias(scratch)       = prod(transfer[i], ST.moment1[i]);
                    noalias(ST.moment1[i]) = prod(scratch, trans(transfer[i]));

                } else {
                */

                // 2nd order Edge focusing.
                if ((int)Ftype == 2 and phi1 != 0.0){
                    fringe2nd(true, i, phi1, ST);
                }

                ST.moment0[i]          = prod(transfer[i], ST.moment0[i]);
    
                if (ST.sim_mode == 2){
                    ST.moment1[i] = prod(transfer[i], ST.moment1[i]); 
                } else {
                    noalias(scratch)       = prod(transfer[i], ST.moment1[i]);
                    noalias(ST.moment1[i]) = prod(scratch, trans(transfer[i]));
                }

                if ((int)Ftype == 2 and phi2 != 0.0){
                    fringe2nd(false, i, phi2, ST);
                }

                //}

                double dphis_temp = ST.moment0[i][state_t::PS_S] - phis_temp;

                if (HdipoleFitMode != 1) {

                    //double    di_bg, Ek00, beta00, gamma00, IonK_Bend;
                    //di_bg     = conf().get<double>("bg");
                    // Dipole reference energy.
                    //Ek00      = (sqrt(sqr(di_bg)+1e0)-1e0)*ST.ref.IonEs;
                    //gamma00   = (Ek00+ST.ref.IonEs)/ST.ref.IonEs;
                    //beta00    = sqrt(1e0-1e0/sqr(gamma00));
                    //IonK_Bend = 2e0*M_PI/(beta00*SampleLambda);

                    // J.B.: this is odd.
                    //ST.real[i].phis  += IonK_Bend*length*MtoMM + dphis_temp;
                    
                    ST.real[i].phis  += ST.real[i].SampleIonK*length*MtoMM + dphis_temp;
                } else
                    ST.real[i].phis  += ST.real[i].SampleIonK*length*MtoMM + dphis_temp;
            }

            ST.calc_rms();
        }
    }

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.

        for(size_t i=0; i<last_real_in.size(); i++) {
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


            if ((int)Ftype == 1) {
                // 1st order Edge focusing.
                double rho = L/phi;

                MomentState::matrix_t edge1, edge2;

                GetEdgeMatrix(rho, phi1, edge1);
                GetEdgeMatrix(rho, phi2, edge2);

                noalias(scratch)     = prod(transfer[i], edge1);
                noalias(transfer[i]) = prod(edge2, scratch);
            }

            get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            noalias(scratch)     = prod(transfer[i], misalign[i]);
            noalias(transfer[i]) = prod(misalign_inv[i], scratch);
        }
    }

    void fringe2nd(int ent, int i, double phif, state_t& ST)
    {
        // 2nd order Edge focusing.

        double h0 = phi/L,
               tan1 = tan(phif),
               tan2 = sqr(tan1),
               sec2 = tan2 + 1,
               aper = conf().get<double>("aper")*MtoMM,
               theta = 2e0*aper*K*h0*(1+sqr(sin(phif)))/cos(phif),
               tan1y = tan(phif-theta),
               qmrel = (ST.real[i].IonZ-ST.ref.IonZ)/ST.ref.IonZ,
               d = 0e0,
               dpdw = 1e0/(sqr(ST.ref.beta)*ST.ref.gamma*ST.ref.IonEs/MeVtoeV);

        double xx  = sqr(ST.moment0[i](state_t::PS_X)),
               xy  = ST.moment0[i](state_t::PS_X)*ST.moment0[i](state_t::PS_Y),
               yy  = sqr(ST.moment0[i](state_t::PS_Y)),
               xpx = ST.moment0[i](state_t::PS_X)*ST.moment0[i](state_t::PS_PX),
               xpy = ST.moment0[i](state_t::PS_X)*ST.moment0[i](state_t::PS_PY),
               ypx = ST.moment0[i](state_t::PS_Y)*ST.moment0[i](state_t::PS_PX),
               ypy = ST.moment0[i](state_t::PS_Y)*ST.moment0[i](state_t::PS_PY),
               xdp = ST.moment0[i](state_t::PS_X)*(ST.moment0[i](state_t::PS_PS)*dpdw +d-qmrel),
               ydp = ST.moment0[i](state_t::PS_Y)*(ST.moment0[i](state_t::PS_PS)*dpdw +d-qmrel);

        xx  += ST.moment1[i](state_t::PS_X, state_t::PS_X),
        xy  += ST.moment1[i](state_t::PS_X, state_t::PS_Y),
        yy  += ST.moment1[i](state_t::PS_Y, state_t::PS_Y),
        xpx += ST.moment1[i](state_t::PS_X, state_t::PS_PX),
        xpy += ST.moment1[i](state_t::PS_X, state_t::PS_PY),
        ypx += ST.moment1[i](state_t::PS_Y, state_t::PS_PX),
        ypy += ST.moment1[i](state_t::PS_Y, state_t::PS_PY),
        xdp += ST.moment1[i](state_t::PS_X, state_t::PS_PS)*dpdw,
        ydp += ST.moment1[i](state_t::PS_Y, state_t::PS_PS)*dpdw;

        MomentState::matrix_t edge1 = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);;

        if (ent) {
            edge1(state_t::PS_PX, state_t::PS_X) =  h0*tan1;
            edge1(state_t::PS_PY, state_t::PS_Y) = -h0*tan1y;

            edge1(state_t::PS_X,  6) = h0*(-tan2*xx + sec2*yy)/2e0;
            edge1(state_t::PS_PX, 6) = h0*( tan2*(xpx-ypy) + h0*(tan1/2e0+tan1*tan2)*yy - tan1*xdp);

            edge1(state_t::PS_Y,  6) = h0*( tan2*xy);
            edge1(state_t::PS_PY, 6) = h0*(-tan2*xpy - sec2*ypx + (tan1-theta*(1+sqr(tan1y)))*ydp);

        } else {
            edge1(state_t::PS_PX, state_t::PS_X) =  h0*tan1;
            edge1(state_t::PS_PY, state_t::PS_Y) = -h0*tan1y;

            edge1(state_t::PS_X,  6) = h0*( tan2*xx-sec2*yy)/2e0;
            edge1(state_t::PS_PX, 6) = h0*( h0*tan1*tan2*(xx-yy)/2e0 - tan2*(xpx-ypy) - tan1*xdp);

            edge1(state_t::PS_Y,  6) = h0*(-tan2*xy);
            edge1(state_t::PS_PY, 6) = h0*( h0*sec2*tan1*xy + tan2*xpy + sec2*ypx +(tan1-theta*(1+sqr(tan1y)))*ydp);
        }

        ST.moment0[i]    = prod(edge1, ST.moment0[i]);

        if (ST.sim_mode == 2){
            ST.moment1[i] = prod(edge1, ST.moment1[i]); 
        } else {
            noalias(scratch) = prod(edge1, ST.moment1[i]);
            noalias(ST.moment1[i])      = prod(scratch, trans(edge1));
        }
    }
};

struct ElementQuad : public MomentElementBase
{
    // Transport matrix for Quadrupole; K = B2/Brho.
    typedef ElementQuad              self_t;
    typedef MomentElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementQuad(const Config& c) : base_t(c) {}
    virtual ~ElementQuad() {}
    virtual const char* type_name() const {return "quadrupole";}

    virtual void assign(const ElementVoid *other) { base_t::assign(other); }

    virtual void recompute_matrix(state_t& ST)
    {
        const double B2= conf().get<double>("B2"),
                     L = conf().get<double>("L")*MtoMM;

        for(size_t i=0; i<last_real_in.size(); i++) {
            // Re-initialize transport matrix.
            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

            double Brho = ST.real[i].beta*(ST.real[i].IonEk+ST.real[i].IonEs)/(C0*ST.real[i].IonZ),
                   K = B2/Brho/sqr(MtoMM);


            // Horizontal plane.
            GetQuadMatrix(L,  K, (unsigned)state_t::PS_X, transfer[i]);
            // Vertical plane.
            GetQuadMatrix(L, -K, (unsigned)state_t::PS_Y, transfer[i]);
            // Longitudinal plane.
            transfer[i](state_t::PS_S, state_t::PS_PS) =
                    -2e0*M_PI/(SampleLambda*ST.real[i].IonEs/MeVtoeV*cube(ST.real[i].bg))*L;

            get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            noalias(scratch)     = prod(transfer[i], misalign[i]);
            noalias(transfer[i]) = prod(misalign_inv[i], scratch);
        }
    }
};

struct ElementSext : public MomentElementBase
{
    // Transport matrix for Sextupole; K = B3/Brho.
    typedef ElementSext              self_t;
    typedef MomentElementBase       base_t;
    typedef typename base_t::state_t state_t;

    double B3, L;
    int step;
    bool thinlens, dstkick;

    ElementSext(const Config& c) : base_t(c) {
        B3= conf().get<double>("B3");
        L = conf().get<double>("L")*MtoMM;
        step = conf().get<double>("step", 1.0);
        thinlens = conf().get<double>("thinlens", 0.0) == 1.0;
        dstkick = conf().get<double>("dstkick", 1.0) == 1.0;
    }
    virtual ~ElementSext() {}
    virtual const char* type_name() const {return "sextupole";}

    virtual void assign(const ElementVoid *other) {
        base_t::assign(other);
        const self_t* O=static_cast<const self_t*>(other);
        B3 = O->B3;
        L = O->L;
        step = O->step;
        thinlens = O->thinlens;
        dstkick = O->dstkick;
    }

    virtual void advance(StateBase& s)
    {
        state_t&  ST = static_cast<state_t&>(s);
        using namespace boost::numeric::ublas;

        ST.recalc();

        last_ref_in = ST.ref;
        last_real_in = ST.real;
        resize_cache(ST);

        const double dL = L/step;

        for(size_t k=0; k<last_real_in.size(); k++) {

            transfer[k] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

            double Brho = ST.real[k].Brho(),
                   K = B3/Brho/cube(MtoMM);

            for(int i=0; i<step; i++){
                double Dx = ST.moment0[k][state_t::PS_X],
                       Dy = ST.moment0[k][state_t::PS_Y],
                       D2x = ST.moment1[k](state_t::PS_X, state_t::PS_X),
                       D2y = ST.moment1[k](state_t::PS_Y, state_t::PS_Y),
                       D2xy = ST.moment1[k](state_t::PS_X, state_t::PS_Y);


                GetSextMatrix(dL,  K, Dx, Dy, D2x, D2y, D2xy, thinlens, dstkick, transfer[k]);

                transfer[k](state_t::PS_S, state_t::PS_PS) =
                        -2e0*M_PI/(SampleLambda*ST.real[k].IonEs/MeVtoeV*cube(ST.real[k].bg))*dL;

                if (strict3Dvelocity) {
                    //@- introduce 3d velocity integrator
                    double sclz3d = sqr(tan(ST.moment0[k][state_t::PS_PX])) + sqr(tan(ST.moment0[k][state_t::PS_PY])) + 1e0;
                    transfer[k](state_t::PS_S, 6) =
                        transfer[k](state_t::PS_S, state_t::PS_PS)*(1e0/sclz3d-1e0)*ST.real[k].IonEk/MeVtoeV;
                }

                get_misalign(ST, ST.real[k], misalign[k], misalign_inv[k]);
                noalias(scratch)     = prod(transfer[k], misalign[k]);
                noalias(transfer[k]) = prod(misalign_inv[k], scratch);

                ST.moment0[k] = prod(transfer[k], ST.moment0[k]);

                scratch  = prod(transfer[k], ST.moment1[k]);
                ST.moment1[k] = prod(scratch, trans(transfer[k]));
            }
        }

        ST.recalc();

        for(size_t k=0; k<last_real_in.size(); k++)
            ST.real[k].phis  += ST.real[k].SampleIonK*length*MtoMM;
        ST.ref.phis   += ST.ref.SampleIonK*length*MtoMM;

        last_ref_out = ST.ref;
        last_real_out = ST.real;

        ST.pos += length;

        ST.calc_rms();
    }
};

struct ElementSolenoid : public MomentElementBase
{
    // Transport (identity) matrix for a Solenoid; K = B/(2 Brho).
    typedef ElementSolenoid          self_t;
    typedef MomentElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementSolenoid(const Config& c) : base_t(c) {}
    virtual ~ElementSolenoid() {}
    virtual const char* type_name() const {return "solenoid";}

    virtual void assign(const ElementVoid *other) { base_t::assign(other); }

    virtual void recompute_matrix(state_t& ST)
    {
        const double B = conf().get<double>("B"),
                     L = conf().get<double>("L")*MtoMM;      // Convert from [m] to [mm].

        for(size_t i=0; i<last_real_in.size(); i++) {
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
        }
    }
};

struct ElementEDipole : public MomentElementBase
{
    // Transport matrix for Electrostatic Dipole with edge focusing.
    typedef ElementEDipole           self_t;
    typedef MomentElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementEDipole(const Config& c) : base_t(c) {}
    virtual ~ElementEDipole() {}
    virtual const char* type_name() const {return "edipole";}

    virtual void assign(const ElementVoid *other) { base_t::assign(other); }

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

        for(size_t i=0; i<last_real_in.size(); i++) {
            double eta0        = (ST.real[i].gamma-1e0)/2e0,
                   Kx          = (1e0-spher+sqr(1e0+2e0*eta0))/sqr(rho),
                   delta_KZ    = ST.real[i].IonZ/ST.ref.IonZ - 1e0,
                   SampleIonK  = 2e0*M_PI/(ST.real[i].beta*SampleLambda);

            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

            GetEBendMatrix(L, phi, fringe_x, fringe_y, kappa, Kx, Ky, ST.ref.IonEs, ST.ref.beta, ST.real[i].gamma,
                           eta0, h, dip_beta, dip_gamma, delta_KZ, SampleIonK, transfer[i]);

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

                noalias(scratch)     = prod(R, transfer[i]);
                noalias(transfer[i]) = prod(scratch, trans(R));
                //TODO: no-op code?  results are unconditionally overwritten
            }

            get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            noalias(scratch)     = prod(transfer[i], misalign[i]);
            noalias(transfer[i]) = prod(misalign_inv[i], scratch);
        }
    }
};

struct ElementEQuad : public MomentElementBase
{
    // Transport matrix for Electrostatic Quadrupole.
    typedef ElementEQuad             self_t;
    typedef MomentElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementEQuad(const Config& c) : base_t(c) {}
    virtual ~ElementEQuad() {}
    virtual const char* type_name() const {return "equad";}

    virtual void assign(const ElementVoid *other) { base_t::assign(other); }

    virtual void recompute_matrix(state_t& ST)
    {
        const double V0   = conf().get<double>("V"),
                     R    = conf().get<double>("radius"),
                     L    = conf().get<double>("L")*MtoMM;

        for(size_t i=0; i<last_real_in.size(); i++) {
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
            //        transfer(state_t::PS_S, state_t::PS_S) = L;

            transfer[i](state_t::PS_S, state_t::PS_PS) =
                    -2e0*M_PI/(SampleLambda*ST.real[i].IonEs/MeVtoeV*cube(ST.real[i].bg))*L;

            get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            noalias(scratch)     = prod(transfer[i], misalign[i]);
            noalias(transfer[i]) = prod(misalign_inv[i], scratch);
        }
    }
};

/*
struct ElementRFQcell : public MomentElementBase
{
    // Transport matrix for Electrostatic Quadrupole.
    typedef ElementRFQcell           self_t;
    typedef MomentElementBase        base_t;
    typedef typename base_t::state_t state_t;

    ElementRFQcell(const Config& c) : base_t(c) {}
    virtual ~ElementRFQcell() {}
    virtual const char* type_name() const {return "rfqcell";}

    virtual void assign(const ElementVoid *other) { base_t::assign(other); }

    virtual void advance(StateBase& s)
    {
        state_t&  ST = static_cast<state_t&>(s);
        using namespace boost::numeric::ublas;

        // IonEk is Es + E_state; the latter is set by user.
        ST.recalc();

        last_ref_in = ST.ref;
        last_real_in = ST.real;
        resize_cache(ST);

        //recompute_matrix(ST); // updates transfer and last_Kenergy_out

        const double V0   = conf().get<double>("V"),
                     R0   = conf().get<double>("r0"),
                     phi0   = conf().get<double>("phi0"),
                     A10   = conf().get<double>("A10"),
                     L    = conf().get<double>("L")*MtoMM;
        const int step = conf().get<double>("step", 1.0);

        const double dL = L/step;


        for(size_t k=0; k<last_real_in.size(); k++) {

            double phi = phi0 + M_PI/(ST.real[k].beta*SampleLambda)*dL,
                   loc = dL/2e0; // Reference orbit position in RFQcell

            //Brho = ST.real[k].beta*(ST.real[k].IonEk+ST.real[k].IonEs)/(C0*ST.real[k].IonZ),
            //K    = 2e0*V0/(C0*ST.real[k].beta*sqr(R))/Brho/sqr(MtoMM);


            //--------------

            transfer[k] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
           
            for(int i=0; i<step; i++){

                double Erho = sqr(ST.real[k].beta)*(ST.real[k].IonEk+ST.real[k].IonEs)/ST.real[k].IonZ,
                       Accl = sqr(M_PI/L)*A10/4e0*cos(M_PI*loc/L), // modulation term
                       I0 = 1e0;
                double Kx = 2e0*V0*dL/Erho*(-1e0/sqr(R0)-Accl)/sqr(MtoMM)*sin(phi),
                       Ky = 2e0*V0*dL/Erho*( 1e0/sqr(R0)-Accl)/sqr(MtoMM)*sin(phi),
                       Kz = V0*dL*A10*I0*(M_PI/L)*sin(M_PI*loc/L)*cos(phi)/Erho;

                // Horizontal plane.
                GetQuadMatrix(L,  Kx, (unsigned)state_t::PS_X, transfer[k]);
                // Vertical plane.
                GetQuadMatrix(L,  Ky, (unsigned)state_t::PS_Y, transfer[k]);



                transfer[k](state_t::PS_S, state_t::PS_PS) =
                        -2e0*M_PI/(SampleLambda*ST.real[k].IonEs/MeVtoeV*cube(ST.real[k].bg))*dL;

                get_misalign(ST, ST.real[k], misalign[k], misalign_inv[k]);

                noalias(scratch)     = prod(transfer[k], misalign[k]);
                noalias(transfer[k]) = prod(misalign_inv[k], scratch);

                //--------------

                ST.moment0[k] = prod(transfer[k], ST.moment0[k]);

                scratch  = prod(transfer[k], ST.moment1[k]);
                ST.moment1[k] = prod(scratch, trans(transfer[k]));

                phi += 2e0*M_PI/(ST.real[k].beta*SampleLambda)*dL;
                loc += dL;

                ST.recalc();

            }
        }

        for(size_t k=0; k<last_real_in.size(); k++)
            ST.real[k].phis  += ST.real[k].SampleIonK*length*MtoMM;
        ST.ref.phis   += ST.ref.SampleIonK*length*MtoMM;

        last_ref_out = ST.ref;
        last_real_out = ST.real;

        ST.pos += length;

        ST.calc_rms();
    }
};
*/


} // namespace

void registerMoment()
{
    Machine::registerState<MomentState>("MomentMatrix");

    Machine::registerElement<ElementSource                 >("MomentMatrix", "source");

    Machine::registerElement<ElementMark                   >("MomentMatrix", "marker");

    Machine::registerElement<ElementBPM                    >("MomentMatrix", "bpm");

    Machine::registerElement<ElementDrift                  >("MomentMatrix", "drift");

    Machine::registerElement<ElementOrbTrim                >("MomentMatrix", "orbtrim");

    Machine::registerElement<ElementSBend                  >("MomentMatrix", "sbend");

    Machine::registerElement<ElementQuad                   >("MomentMatrix", "quadrupole");

    Machine::registerElement<ElementSext                   >("MomentMatrix", "sextupole");

    Machine::registerElement<ElementSolenoid               >("MomentMatrix", "solenoid");

    Machine::registerElement<ElementRFCavity               >("MomentMatrix", "rfcavity");

    Machine::registerElement<ElementStripper               >("MomentMatrix", "stripper");

    Machine::registerElement<ElementEDipole                >("MomentMatrix", "edipole");

    Machine::registerElement<ElementEQuad                  >("MomentMatrix", "equad");

//    Machine::registerElement<ElementRFQcell                 >("MomentMatrix", "rfqcell");
}
