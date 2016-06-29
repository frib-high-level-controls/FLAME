#ifndef RF_CAVITY_H
#define RF_CAVITY_H

#endif // RF_CAVITY_H

#include <boost/numeric/ublas/matrix.hpp>

#include "moment.h"
#include "util.h"

// Phase space dimension; including vector for orbit/1st moment.
# define PS_Dim MomentState::maxsize // Set to 7; to include orbit.


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
    void show(std::ostream& strm, const int) const;
    void show(std::ostream& strm) const;
};

struct ElementRFCavity : public MomentElementBase
{
    // Transport matrix for an RF Cavity.
    typedef ElementRFCavity          self_t;
    typedef MomentElementBase       base_t;
    typedef typename base_t::state_t state_t;

    struct RawParams {
        std::string name, type;
        double length, aperature, E0;
    };
    std::vector<RawParams> lattice; // from axisData_*.txt

    numeric_table mlptable, // from CaviMlp_*.txt
                  CavData; // from thinlenlon_*.txt

    std::vector<CavTLMLineType> CavTLMLineTab; // from lattice, for each charge state
    double phi_ref;
    int MpoleLevel;
    int cavi;
    bool forcettfcalc;

    ElementRFCavity(const Config& c);

    void GetCavMatParams(const int cavi,
                         const double beta_tab[], const double gamma_tab[], const double IonK[],
                         CavTLMLineType& lineref) const;

    void GetCavMat(const int cavi, const int cavilabel, const double Rm, Particle &real,
                   const double EfieldScl, const double IonFyi_s,
                   const double IonEk_s, const double fRF, state_t::matrix_t &M,
                   CavTLMLineType &linetab) const;

    void GenCavMat2(const int cavi, const double dis, const double EfieldScl, const double TTF_tab[],
                   const double beta_tab[], const double gamma_tab[], const double Lambda,
                   Particle &real, const double IonFys[], const double Rm, state_t::matrix_t &M,
                   const CavTLMLineType& linetab) const;

    void PropagateLongRFCav(Particle &ref);

    void InitRFCav(Particle &real, state_t::matrix_t &M, CavTLMLineType &linetab);

    void GetCavBoost(const numeric_table &CavData, Particle &state, const double IonFy0, const double fRF,
                     const double EfieldScl, double &IonFy) const;

    void TransFacts(const int cavilabel, double beta, const double CaviIonK, const int gaplabel, const double EfieldScl,
                    double &Ecen, double &T, double &Tp, double &S, double &Sp, double &V0) const;

    void TransitFacMultipole(const int cavi, const std::string &flabel, const double CaviIonK,
                             double &T, double &S) const;

    virtual ~ElementRFCavity() {}

    virtual void assign(const ElementVoid *other) {
        base_t::assign(other);
        const self_t* O=static_cast<const self_t*>(other);
        lattice       = O->lattice;
        mlptable      = O->mlptable;
        CavData       = O->CavData;
        CavTLMLineTab = O->CavTLMLineTab;
        phi_ref       = O->phi_ref;
        MpoleLevel    = O->MpoleLevel;
        cavi          = O->cavi;
        forcettfcalc  = O->forcettfcalc;
    }

    virtual void advance(StateBase& s)
    {
        state_t&  ST = static_cast<state_t&>(s);
        using namespace boost::numeric::ublas;

        // IonEk is Es + E_state; the latter is set by user.
        ST.recalc();

        if(!check_cache(ST)) {
            last_ref_in = ST.ref;
            last_real_in = ST.real;
            resize_cache(ST);
            state_t::matrix_t dummy_mat;

            // need to re-calculate energy dependent terms

            for(size_t i=0; i<last_real_in.size(); i++)
                get_misalign(ST, ST.real[i], misalign[i], dummy_mat);
            // only use entrance misalignment matrix, exit is dummy.

            recompute_matrix(ST); // updates transfer and last_Kenergy_out

            for(size_t i=0; i<last_real_in.size(); i++)
                get_misalign(ST, ST.real[i], dummy_mat, misalign_inv[i]);
            // only use exit misalignment matrix, entrance is dummy.

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

        for(size_t i=0; i<last_real_in.size(); i++) {
            ST.moment0[i] = prod(misalign[i], ST.moment0[i]);
            ST.moment0[i] = prod(transfer[i], ST.moment0[i]);

            ST.moment0[i][state_t::PS_S]  = ST.real[i].phis - ST.ref.phis;
            ST.moment0[i][state_t::PS_PS] = (ST.real[i].IonEk-ST.ref.IonEk)/MeVtoeV;

            ST.moment0[i] = prod(misalign_inv[i], ST.moment0[i]);

            scratch  = prod(misalign[i], ST.moment1[i]);
            ST.moment1[i] = prod(scratch, trans(misalign[i]));

            scratch  = prod(transfer[i], ST.moment1[i]);
            ST.moment1[i] = prod(scratch, trans(transfer[i]));

            scratch  = prod(misalign_inv[i], ST.moment1[i]);
            ST.moment1[i] = prod(scratch, trans(misalign_inv[i]));
        }

        ST.calc_rms();
    }

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix. and update ST.ref and ST.real[]

        CavTLMLineTab.resize(last_real_in.size());

        PropagateLongRFCav(ST.ref);

        for(size_t i=0; i<last_real_in.size(); i++) {
            // TODO: 'transfer' is overwritten in InitRFCav()?
            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
            transfer[i](state_t::PS_X, state_t::PS_PX) = length;
            transfer[i](state_t::PS_Y, state_t::PS_PY) = length;

            // J.B. Bug in TLM.
            //double SampleIonK = ST.real[i].SampleIonK;


            InitRFCav(ST.real[i], transfer[i], CavTLMLineTab[i]);

            // J.B. Bug in TLM.
            //ST.real[i].SampleIonK = SampleIonK;
        }
   }

    virtual const char* type_name() const {return "rfcavity";}
};

