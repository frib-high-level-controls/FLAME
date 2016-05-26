
#include <fstream>

#include <limits>

#include <boost/numeric/ublas/lu.hpp>
#include <boost/lexical_cast.hpp>

#include "scsi/constants.h"
#include "scsi/moment2.h"

#include "scsi/moment2.h"
#include "scsi/moment2_sup.h"
#include "scsi/rf_cavity.h"
#include "scsi/chg_stripper.h"

#include "scsi/h5loader.h"


namespace {
// http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
// by LU-decomposition.
void inverse(Moment2ElementBase::value_t& out, const Moment2ElementBase::value_t& in)
{
    using boost::numeric::ublas::permutation_matrix;
    using boost::numeric::ublas::lu_factorize;
    using boost::numeric::ublas::lu_substitute;
    using boost::numeric::ublas::identity_matrix;

    Moment2ElementBase::value_t scratch(in); // copy
    permutation_matrix<size_t> pm(scratch.size1());
    if(lu_factorize(scratch, pm)!=0)
        throw std::runtime_error("Failed to invert matrix");
    out.assign(identity_matrix<double>(scratch.size1()));
    //out = identity_matrix<double>(scratch.size1());
    lu_substitute(scratch, pm, out);
}

template<typename ARR>
void load_storage(ARR& to, const Config& conf, const std::string& name)
{
    try{
        const std::vector<double>& val(conf.get<std::vector<double> >(name));
        if(to.size()!=val.size()) {
            std::ostringstream strm;
            strm<<"Array "<<name<<" must have "<<to.size()<<" elements, not "<<val.size();
            throw std::invalid_argument(strm.str());
        }
        std::copy(val.begin(), val.end(), to.begin());

    }catch(key_error&){
        throw std::invalid_argument(name+" not defined");
        // default to identity
    }catch(boost::bad_any_cast&){
        throw std::invalid_argument(name+" has wrong type (must be vector)");
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

        for(size_t i=0; i<ics.size(); i++) {
            std::string num(boost::lexical_cast<std::string>(i));

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
        }
    } else {
        real.resize(1); // hack, ensure at least one element so getArray() can return some pointer
    }
}

Moment2State::~Moment2State() {}

Moment2State::Moment2State(const Moment2State& o, clone_tag t)
    :StateBase(o, t)
    ,ref(o.ref)
    ,real(o.real)
    ,moment0(o.moment0)
    ,moment1(o.moment1)
    ,moment0_env(o.moment0_env)
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
    moment1_env = O->moment1_env;
    StateBase::assign(other);
}

void Moment2State::show(std::ostream& strm) const
{
    int j, k;

    if(real.empty()) {
        strm<<"\nState: empty\n";
        return;
    }

    strm << std::scientific << std::setprecision(8)
         << "\nState:\n  energy [eV] =\n" << std::setw(20) << real[0].IonEk << "\n  moment0 =\n    ";
    for (k = 0; k < Moment2State::maxsize; k++)
        strm << std::scientific << std::setprecision(8) << std::setw(16) << moment0[0](k);
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
    ,scratch(state_t::maxsize, state_t::maxsize)
{
    inverse(misalign_inv, misalign);
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

void Moment2ElementBase::show(std::ostream& strm) const
{
    using namespace boost::numeric::ublas;
    ElementVoid::show(strm);
    /*
    strm<<"Length "<<length<<"\n"
          "Transfer: "<<transfer<<"\n"
          "Transfer Raw: "<<transfer_raw<<"\n"
          "Mis-align: "<<misalign<<"\n";
          */
}

void Moment2ElementBase::advance(StateBase& s)
{
    state_t& ST = static_cast<state_t&>(s);
    using namespace boost::numeric::ublas;

    // IonEk is Es + E_state; the latter is set by user.
    ST.recalc();

    std::cout<<"Advance Element "<<index<<" '"<<name<<"'\n";

    bool cache_hit = ST.real.size()==last_Kenergy_in.size();
    if(cache_hit) {
        for(size_t n=0; n<ST.real.size(); n++) {
            if(ST.real[n].IonEk!=last_Kenergy_in[n]) {
                cache_hit = false;
                break;
            }
        }
    }

    if(!cache_hit)
    {
        // need to re-calculate energy dependent terms

        last_Kenergy_in.resize(ST.real.size());
        last_Kenergy_out.resize(ST.real.size());
        transfer.resize(ST.real.size(), boost::numeric::ublas::identity_matrix<double>(state_t::maxsize));
        transfer_raw.resize(ST.real.size(), boost::numeric::ublas::identity_matrix<double>(state_t::maxsize));

        recompute_matrix(ST); // updates transfer and last_Kenergy_out

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            // TODO: broken somehow?
            //noalias(scratch)  = prod(misalign, transfer[i]);
            //noalias(transfer[i]) = prod(scratch, misalign_inv);

            ST.real[i].recalc();
        }
    }

    // recompute_matrix only called when ST.IonEk != last_Kenergy_in.
    // Matrix elements are scaled with particle energy.

    ST.pos += length;

    const char * const t_name = type_name();
    bool isrf  = strcmp(t_name, "rfcavity")==0,
         isbend= strcmp(t_name, "sbend")==0;

    if(!isrf)
        ST.ref.phis   += ST.ref.SampleIonK*length*MtoMM;

    std::fill(ST.moment0_env.begin(), ST.moment0_env.end(), 0.0);
    std::fill(ST.moment1_env.data().begin(), ST.moment1_env.data().end(), 0.0);
    double totalQ = 0.0;

    for(size_t i=0; i<last_Kenergy_in.size(); i++) {
        if(!isrf && !isbend) {
            ST.real[i].phis  += ST.real[i].SampleIonK*length*MtoMM;
            ST.real[i].IonEk  = last_Kenergy_out[i];
        }

        double   phis_temp;
        if(isbend)
            phis_temp = ST.moment0[i][state_t::PS_S];

        std::cout<<"moment0 in  "<<ST.moment0[i]
               <<"\ntransfer    "<<transfer[i]<<"\n";

        ST.moment0[i] = prod(transfer[i], ST.moment0[i]);

        std::cout<<"moment0 out "<<ST.moment0[i]<<"\n";

        if(isrf) {
            ST.moment0[i][state_t::PS_S]  = ST.real[i].phis - ST.ref.phis;
            ST.moment0[i][state_t::PS_PS] = (ST.real[i].IonEk-ST.ref.IonEk)/MeVtoeV;

        } else if(isbend) {
            ST.real[i].phis  += ST.real[i].SampleIonK*length*MtoMM + ST.moment0[i][state_t::PS_S] - phis_temp;

            ST.real[i].IonEk  = last_Kenergy_out[i];

        }

        noalias(scratch) = prod(transfer[i], ST.moment1[i]);
        noalias(ST.moment1[i]) = prod(scratch, trans(transfer[i]));

        ST.moment0_env += ST.moment0[i]*ST.real[i].IonQ;
        ST.moment1_env += ST.moment1[i]*ST.real[i].IonQ;
        totalQ += ST.real[i].IonQ;
    }

    ST.moment0_env /= totalQ;
    ST.moment1_env /= totalQ;
}

void Moment2ElementBase::recompute_matrix(state_t& ST)
{
    // Default no-op

    for(size_t i=0; i<last_Kenergy_in.size(); i++) {
        last_Kenergy_out[i] = ST.real[i].IonEk;
    }

    last_Kenergy_in = last_Kenergy_out; // no energy gain
    transfer = transfer_raw;
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

    virtual void show(std::ostream& strm) const
    {
        ElementVoid::show(strm);
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
        double L = length*MtoMM; // Convert from [m] to [mm].

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            transfer_raw[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
            transfer_raw[i](state_t::PS_X, state_t::PS_PX) = L;
            transfer_raw[i](state_t::PS_Y, state_t::PS_PY) = L;
            transfer_raw[i](state_t::PS_S, state_t::PS_PS) =
                    -2e0*M_PI/(SampleLambda*ST.real[i].IonEs/MeVtoeV*cube(ST.real[i].bg))*L;

            transfer[i] = transfer_raw[i];

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
        double theta_x = conf().get<double>("theta_x", 0e0),
               theta_y = conf().get<double>("theta_y", 0e0);

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
            transfer[i](state_t::PS_PX, 6) = theta_x*ST.ref.IonZ/ST.real[i].IonZ;
            transfer[i](state_t::PS_PY, 6) = theta_y*ST.ref.IonZ/ST.real[i].IonZ;

            last_Kenergy_in[i] = last_Kenergy_out[i] = ST.real[i].IonEk; // no energy gain
        }
    }
};

struct ElementSBend : public Moment2ElementBase
{
    // Transport matrix for Gradient Sector Bend; with edge focusing (cylindrical coordinates).

    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;
    ElementSBend(const Config& c) : base_t(c) {}
    virtual ~ElementSBend() {}
    virtual const char* type_name() const {return "sbend";}

    virtual void recompute_matrix(state_t& ST)
    {
        const
        double L    = conf().get<double>("L")*MtoMM,
               phi  = conf().get<double>("phi")*M_PI/180e0,
               phi1 = conf().get<double>("phi1")*M_PI/180e0,
               phi2 = conf().get<double>("phi2")*M_PI/180e0,
               rho  = L/phi,
               K    = conf().get<double>("K", 0e0)/sqr(MtoMM),
               Kx   = K + 1e0/sqr(rho),
               Ky   = -K;

        value_t edge1(boost::numeric::ublas::identity_matrix<double>(state_t::maxsize)),
                edge2(edge1), base(edge2);

        // Edge focusing.
        GetEdgeMatrix(rho, phi1, edge1);
        // Horizontal plane.
        GetQuadMatrix(L, Kx, (unsigned)state_t::PS_X, base);
        // Vertical plane.
        GetQuadMatrix(L, Ky, (unsigned)state_t::PS_Y, base);

        // Include dispersion.
        double dx, sx;
        if (Kx == 0e0) {
            dx = sqr(L)/2e0;
            sx = L;
        } else if (Kx > 0e0) {
            dx = (1e0-cos(sqrt(Kx)*L))/Kx;
            sx = sin(sqrt(Kx)*L)/sqrt(Kx);
        } else {
            dx = (1e0-cosh(sqrt(-Kx)*L))/Kx;
            sx = sin(sqrt(Kx)*L)/sqrt(Kx);
        }

        base(state_t::PS_X,  state_t::PS_PS) = dx/(rho*sqr(ST.ref.beta)*ST.ref.gamma*ST.ref.IonEs/MeVtoeV);
        base(state_t::PS_PX, state_t::PS_PS) = sx/(rho*sqr(ST.ref.beta)*ST.ref.gamma*ST.ref.IonEs/MeVtoeV);
        base(state_t::PS_S,  state_t::PS_X)  = sx/rho*ST.ref.SampleIonK;
        base(state_t::PS_S,  state_t::PS_PX) = dx/rho*ST.ref.SampleIonK;
        // Low beta approximation.
        base(state_t::PS_S,  state_t::PS_PS) =
                ((L-sx)/(Kx*sqr(rho))-L/sqr(ST.ref.gamma))*ST.ref.SampleIonK
                /(sqr(ST.ref.beta)*ST.ref.gamma*ST.ref.IonEs/MeVtoeV);

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            double qmrel = (ST.real[i].IonZ-ST.ref.IonZ)/ST.ref.IonZ;

            transfer_raw[i] = base;

            // Add dipole terms.
            transfer_raw[i](state_t::PS_X,  6) = -dx/rho*qmrel;
            transfer_raw[i](state_t::PS_PX, 6) = -sx/rho*qmrel;
            transfer_raw[i](state_t::PS_S,  6) = -(L-sx)/(Kx*sqr(rho))*ST.ref.SampleIonK*qmrel;

            // Edge focusing.
            GetEdgeMatrix(rho, phi2, edge2);

            transfer_raw[i] = prod(transfer_raw[i], edge1);
            transfer_raw[i] = prod(edge2, transfer_raw[i]);

            // Longitudinal plane.
            // For total path length.
    //        transfer_raw(state_t::PS_S,  state_t::PS_S) = L;

            transfer[i] = transfer_raw[i];

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
            transfer_raw[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

            double Brho = ST.real[i].beta*(ST.real[i].IonEk+ST.real[i].IonEs)/(C0*ST.real[i].IonZ),
                   K = B2/Brho/sqr(MtoMM);

            // Horizontal plane.
            GetQuadMatrix(L,  K, (unsigned)state_t::PS_X, transfer_raw[i]);
            // Vertical plane.
            GetQuadMatrix(L, -K, (unsigned)state_t::PS_Y, transfer_raw[i]);
            // Longitudinal plane.
    //        transfer_raw(state_t::PS_S, state_t::PS_S) = L;

            transfer_raw[i](state_t::PS_S, state_t::PS_PS) =
                    -2e0*M_PI/(SampleLambda*ST.real[i].IonEs/MeVtoeV*cube(ST.real[i].bg))*L;

            transfer[i] = transfer_raw[i];

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
            transfer_raw[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

            double Brho = ST.real[i].beta*(ST.real[i].IonEk+ST.real[i].IonEs)/(C0*ST.real[i].IonZ),
                   K    = B/(2e0*Brho)/MtoMM;

            GetSolMatrix(L, K, transfer_raw[i]);

            transfer_raw[i](state_t::PS_S, state_t::PS_PS) =
                    -2e0*M_PI/(SampleLambda*ST.real[i].IonEs/MeVtoeV*cube(ST.real[i].bg))*L;

            transfer[i] = transfer_raw[i];

            last_Kenergy_in[i] = last_Kenergy_out[i] = ST.real[i].IonEk; // no energy gain
        }
    }
};

struct ElementEDipole : public Moment2ElementBase
{
    // Transport matrix for Electrostatic Dipole; with edge focusing.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    ElementEDipole(const Config& c) : base_t(c) {}
    virtual ~ElementEDipole() {}
    virtual const char* type_name() const {return "edipole";}

    virtual void recompute_matrix(state_t& ST)
    {
        //double L = c.get<double>("L")*MtoMM;

//        mscpTracker.position=fribnode.position+fribnode.length;
//        double gamma=(Ek[ii_state]+FRIBPara.ionEs)/FRIBPara.ionEs;
//        double beta=Math.sqrt(1.0-1.0/(gamma*gamma));
//        double SampleionK=2*Math.PI/(beta*SampleLamda*1e3); //rad/mm, marking the currut using ionK
//        double length=fribnode.length; // m
//        double bangle=fribnode.attribute[1]; // degree
//        double beta00=fribnode.attribute[2]; // The voltage is set that the particle with beta_EB would be bended ideally
//        double n1=fribnode.attribute[3]; // 0 is for cylindrical, and 1 is for spherical
//        int vertical=(int) fribnode.attribute[4]; // 0 for horizontal, 1 is for vertical ebend matrix
//        double fringeX=fribnode.attribute[5]; // X fringe
//        double fringeY=fribnode.attribute[6]; // Y fringe
//        double kappa=(int) fribnode.attribute[7]; // voltage asymetry factor
//        double gamma0=tlmPara.Gamma_tab.get(lattcnt+1)[1];
//        double beta0=tlmPara.Beta_tab.get(lattcnt+1)[1];
//        double h=1.0; // 0 for magnetic and 1 for electrostatic
//        double rh0=length/(bangle/180*Math.PI);

//        double Fy_temp=TransVector[ii_state].getElem(4);
//        double dFy_temp=TransVector[ii_state].getElem(4)-Fy_temp;
//        Fy_abs[ii_state]=Fy_abs[ii_state]+SampleionK*1000*fribnode.length+dFy_temp;

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
            transfer_raw[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

            double Brho = ST.real[i].beta*(ST.real[i].IonEk+ST.real[i].IonEs)/(C0*ST.real[i].IonZ),
                   K    = 2e0*V0/(C0*ST.real[i].beta*sqr(R))/Brho/sqr(MtoMM);

            // Horizontal plane.
            GetQuadMatrix(L,  K, (unsigned)state_t::PS_X, transfer_raw[i]);
            // Vertical plane.
            GetQuadMatrix(L, -K, (unsigned)state_t::PS_Y, transfer_raw[i]);
            // Longitudinal plane.
            //        transfer_raw(state_t::PS_S, state_t::PS_S) = L;

            transfer_raw[i](state_t::PS_S, state_t::PS_PS) =
                    -2e0*M_PI/(SampleLambda*ST.real[i].IonEs/MeVtoeV*cube(ST.real[i].bg))*L;

            transfer[i] = transfer_raw[i];

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

    Machine::registerElement<ElementSource                 >("MomentMatrix2",   "source");

    Machine::registerElement<ElementMark                   >("MomentMatrix2",   "marker");

    Machine::registerElement<ElementBPM                    >("MomentMatrix2",   "bpm");

    Machine::registerElement<ElementDrift                  >("MomentMatrix2",   "drift");

    Machine::registerElement<ElementOrbTrim                >("MomentMatrix2",   "orbtrim");

    Machine::registerElement<ElementSBend                  >("MomentMatrix2",   "sbend");

    Machine::registerElement<ElementQuad                   >("MomentMatrix2",   "quadrupole");

    Machine::registerElement<ElementSolenoid               >("MomentMatrix2",   "solenoid");

    Machine::registerElement<ElementRFCavity               >("MomentMatrix2",   "rfcavity");

    Machine::registerElement<ElementStripper               >("MomentMatrix2",   "stripper");

    Machine::registerElement<ElementEDipole                >("MomentMatrix2",   "edipole");

    Machine::registerElement<ElementEDipole                >("MomentMatrix2",   "equad");

    Machine::registerElement<ElementGeneric                >("MomentMatrix2",   "generic");
}
