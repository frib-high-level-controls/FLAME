#ifndef RF_CAVITY_H
#define RF_CAVITY_H

#endif // RF_CAVITY_H

#include <boost/numeric/ublas/matrix.hpp>

#include "moment.h"
#include "util.h"

// Phase space dimension; including vector for orbit/1st moment.
#define PS_Dim MomentState::maxsize // Set to 7; to include orbit.

#ifdef DEFPATH
    #define defpath DEFPATH
#else
    #define defpath "."
#endif


class CavDataType {
// Cavity on-axis longitudinal electric field vs. s.
public:
    std::vector<double> s,     // s coordinate [m]
                        Elong; // Longitudinal Electric field [V/m].

    void RdData(std::istream &inf);
    void show(std::ostream&, const int) const;
    void show(std::ostream&) const;
};


class CavTLMLineType {
public:
    std::vector<double> s;         // Longitudinal position [m].
    std::vector<std::string> Elem;
    std::vector<double> E0,
                        T,
                        S,
                        Accel;

    void clear(void);
    void set(const double, const std::string &, const double,
             const double, const double, const double);
    void show(const int) const;
    void show() const;
};

struct ElementRFCavity : public MomentElementBase
{
    // Transport matrix for an RF Cavity.
    typedef ElementRFCavity          self_t;
    typedef MomentElementBase        base_t;
    typedef typename base_t::state_t state_t;

    struct RawParams {
        std::string name, type;
        double length, aperature, E0;
        // vector Tfit and Sfit always have ten elements
        std::vector<double> Tfit, Sfit;
    };
    std::vector<RawParams> lattice; // from thinlenlon_*.txt

    numeric_table mlptable, // from CaviMlp_*.txt
                  CavData; // from axisData_*.txt

    std::string CavType,
                DataPath,
                DataFile;

    std::vector<double> SynAccTab;

    bool have_RefNrm,
         have_SynComplex,
         have_EkLim,
         have_NrLim;

    double RefNrm; // Reference scale factor q0*1.0/m0

    std::vector<double> SynComplex, // Fitting model coefficients
                        EkLim,      // Limits for incident energy
                        NrLim;      // Limits for normalization factor q*scl/m

    double calFitPow(double kfac, const std::vector<double>& Tfit) const;
    static std::map<std::string,boost::shared_ptr<Config> > CavConfMap;

    std::vector<CavTLMLineType> CavTLMLineTab; // from lattice, for each charge state
    double fRF,    // RF frequency [Hz]
           IonFys, // Synchrotron phase [rad].
           phi_ref,
           cRm;
    int cavi;
    bool forcettfcalc;

    unsigned MpoleLevel,
             EmitGrowth;

    ElementRFCavity(const Config& c);

    void LoadCavityFile(const Config& c);

    void GetCavMatParams(const int cavi,
                         const double beta_tab[], const double gamma_tab[], const double IonK[],
                         CavTLMLineType& lineref) const;

    void GetCavMat(const int cavi, const int cavilabel, const double Rm, Particle &real,
                   const double EfieldScl, const double IonFyi_s,
                   const double IonEk_s, state_t::matrix_t &M,
                   CavTLMLineType &linetab) const;

    void GetCavMatGeneric(Particle &real, const double EfieldScl, const double IonFyi_s,
                                const double IonEk_s, state_t::matrix_t &M, CavTLMLineType &linetab) const;

    void GenCavMat2(const int cavi, const double dis, const double EfieldScl, const double TTF_tab[],
                   const double beta_tab[], const double gamma_tab[], const double Lambda,
                   Particle &real, const double IonFys[], const double Rm, state_t::matrix_t &M,
                   const CavTLMLineType& linetab) const;

    void PropagateLongRFCav(Particle &ref, double &phi_ref) const;

    void calRFcaviEmitGrowth(const state_t::matrix_t &matIn, Particle &state, const int n,
                             const double betaf, const double gamaf,
                             const double aveX2i, const double cenX, const double aveY2i, const double cenY,
                             state_t::matrix_t &matOut);

    void InitRFCav(Particle &real, state_t::matrix_t &M, CavTLMLineType &linetab);

    void GetCavBoost(const numeric_table &CavData, Particle &state, const double IonFy0,
                     const double EfieldScl, double &IonFy) const;

    void TransFacts(const int cavilabel, double beta, const double CaviIonK, const int gaplabel, const double EfieldScl,
                    double &Ecen, double &T, double &Tp, double &S, double &Sp, double &V0) const;

    void TransitFacMultipole(const int cavi, const std::string &flabel, const double CaviIonK,
                             double &T, double &S) const;

    virtual ~ElementRFCavity() {}

    virtual void assign(const ElementVoid *other) {
        base_t::assign(other);
        const self_t* O=static_cast<const self_t*>(other);
        // *all* member variables must be assigned here or reconfigure() will result in inconsistancy
        lattice       = O->lattice;
        mlptable      = O->mlptable;
        CavData       = O->CavData;
        CavTLMLineTab = O->CavTLMLineTab;
        fRF           = O->fRF;
        IonFys        = O->IonFys;
        phi_ref       = O->phi_ref;
        MpoleLevel    = O->MpoleLevel;
        EmitGrowth    = O->EmitGrowth;
        cRm           = O->cRm;
        cavi          = O->cavi;
        forcettfcalc  = O->forcettfcalc;
    }

    virtual void advance(StateBase& s)
    {
        state_t&  ST = static_cast<state_t&>(s);
        using namespace boost::numeric::ublas;

        double x0[2], x2[2], s0[2];

        // IonEk is Es + E_state; the latter is set by user.
        ST.recalc();

        if(!check_cache(ST) && !ST.retreat) {
            last_ref_in = ST.ref;
            last_real_in = ST.real;
            resize_cache(ST);
            // need to re-calculate energy dependent terms

            std::string newtype = conf().get<std::string>("cavtype");
            if (CavType != newtype){
                lattice.clear();
                SynAccTab.clear();
                LoadCavityFile(conf());
            } else if (CavType == "Generic") {
                std::string newfile = conf().get<std::string>("Eng_Data_Dir", defpath);
                newfile += "/" + conf().get<std::string>("datafile");
                if (DataFile != newfile) {
                    lattice.clear();
                    SynAccTab.clear();
                    LoadCavityFile(conf());
                }
            }

            recompute_matrix(ST); // updates transfer and last_Kenergy_out

            for(size_t i=0; i<last_real_in.size(); i++)
                get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            ST.recalc();

            last_ref_out = ST.ref;
            last_real_out = ST.real;
        } else if(ST.retreat){
            if (!check_backward(ST))
                throw std::runtime_error(SB()<<
                    "Backward propagation error at " << ST.next_elem << ": beam state does not match to the previous propagation.");

            ST.ref.phis -= (last_ref_out.phis - last_ref_in.phis);
            ST.ref.IonEk = last_ref_in.IonEk;
            for(size_t k=0; k<last_real_in.size(); k++) {
                ST.real[k].phis -= (last_real_out[k].phis - last_real_in[k].phis);
                ST.real[k].IonEk = last_real_in[k].IonEk;
                get_misalign(ST, ST.real[k], misalign[k], misalign_inv[k]);
            }

            ST.recalc();

        } else {
            ST.ref = last_ref_out;
            assert(last_real_out.size()==ST.real.size()); // should be true if check_cache() -> true
            std::copy(last_real_out.begin(),
                      last_real_out.end(),
                      ST.real.begin());
        }
        // note that calRFcaviEmitGrowth() assumes real[] isn't changed after this point

        if(!ST.retreat){
            // Forward propagation
            ST.pos += length;
            for(size_t i=0; i<last_real_in.size(); i++) {
                ST.moment0[i] = prod(misalign[i], ST.moment0[i]);

                // Inconsistency in TLM; orbit at entrace should be used to evaluate emittance growth.
                x0[0]  = ST.moment0[i][state_t::PS_X];
                x0[1]  = ST.moment0[i][state_t::PS_Y];
                x2[0]  = ST.moment1[i](0, 0);
                x2[1]  = ST.moment1[i](2, 2);

                // reset extra parameters in transfer matrix
                transfer[i](state_t::PS_S, 6) = 0.0;
                transfer[i](state_t::PS_PS, 6) = 0.0;

                ST.moment0[i] = prod(transfer[i], ST.moment0[i]);

                // combine new z and zp centroid to transfer matrix for backward propagation
                s0[0] = (ST.real[i].phis - ST.ref.phis);
                s0[1] = (ST.real[i].IonEk - ST.ref.IonEk)/MeVtoeV;
                transfer[i](state_t::PS_S, 6) = - ST.moment0[i][state_t::PS_S] + s0[0];
                transfer[i](state_t::PS_PS, 6) = - ST.moment0[i][state_t::PS_PS] + s0[1];

                // insert new z and zp centroid
                ST.moment0[i][state_t::PS_S]  = s0[0];
                ST.moment0[i][state_t::PS_PS] = s0[1];

                ST.moment0[i] = prod(misalign_inv[i], ST.moment0[i]);

                scratch = prod(misalign[i], ST.moment1[i]);
                ST.moment1[i] = prod(scratch, trans(misalign[i]));

                scratch = prod(transfer[i], ST.moment1[i]);
                ST.moment1[i] = prod(scratch, trans(transfer[i]));

                if (EmitGrowth) {
                    calRFcaviEmitGrowth(ST.moment1[i], ST.ref, i, ST.real[i].beta, ST.real[i].gamma, x2[0], x0[0], x2[1], x0[1], scratch);
                    ST.moment1[i] = scratch;
                }

                scratch = prod(misalign_inv[i], ST.moment1[i]);
                ST.moment1[i] = prod(scratch, trans(misalign_inv[i]));

                scratch = prod(transfer[i], misalign[i]);
                scratch = prod(misalign_inv[i], scratch);
                ST.transmat[i] = scratch;
            }
        } else {
            // Backward propagation
            ST.pos -= length;
            value_t invmat = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
            for(size_t i=0; i<last_real_in.size(); i++) {
                scratch = prod(transfer[i], misalign[i]);
                scratch = prod(misalign_inv[i], scratch);

                inverse(invmat, scratch);

                ST.moment0[i] = prod(invmat, ST.moment0[i]);

                scratch  = prod(invmat, ST.moment1[i]);
                ST.moment1[i] = prod(scratch, trans(invmat));
                ST.transmat[i] = invmat;
            }

        }

        ST.last_caviphi0 = fmod(phi_ref*180e0/M_PI, 360e0); // driven phase [degree]
        ST.calc_rms();
    }

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix. and update ST.ref and ST.real[]

        CavTLMLineTab.resize(last_real_in.size());

        PropagateLongRFCav(ST.ref, phi_ref);

        for(size_t i=0; i<last_real_in.size(); i++) {
            // TODO: 'transfer' is overwritten in InitRFCav()?
            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
            transfer[i](state_t::PS_X, state_t::PS_PX) = length;
            transfer[i](state_t::PS_Y, state_t::PS_PY) = length;

            // J.B. Bug in TLM.
            double SampleIonK = ST.real[i].SampleIonK;

            InitRFCav(ST.real[i], transfer[i], CavTLMLineTab[i]);

            // J.B. Bug in TLM.
            ST.real[i].SampleIonK = SampleIonK;
        }
   }

    virtual const char* type_name() const {return "rfcavity";}
};

