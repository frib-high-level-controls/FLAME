#ifndef RF_CAVITY_H
#define RF_CAVITY_H

#endif // RF_CAVITY_H

#include <boost/numeric/ublas/matrix.hpp>

#include "moment2.h"

// Phase space dimension; including vector for orbit/1st moment.
# define PS_Dim Moment2State::maxsize // Set to 7; to include orbit.


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

struct MLPtable {
    typedef Moment2State::matrix_t value_t;

    typedef std::map<std::string, size_t> colnames_t;
    colnames_t colnames;

    value_t table;

    void read(const std::string&);
};

struct ElementRFCavity : public Moment2ElementBase
{
    // Transport matrix for an RF Cavity.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    struct RawParams {
        std::string name, type;
        double length, aperature, E0;
    };
    std::vector<RawParams> lattice; // from axisData_*.txt

    CavDataType    CavData; // from thinlenlon_*.txt
    std::vector<CavTLMLineType> CavTLMLineTab; // from lattice, for each charge state
    double phi_ref;

    MLPtable mlptable;

    ElementRFCavity(const Config& c);

    void GetCavMatParams(const int cavi,
                         const double beta_tab[], const double gamma_tab[], const double IonK[],
                         CavTLMLineType& lineref) const;

    void GetCavMat(const int cavi, const int cavilabel, const double Rm, Particle &real,
                   const double EfieldScl, const double IonFyi_s,
                   const double IonEk_s, const double fRF, state_t::matrix_t &M,
                   CavTLMLineType &linetab);

    void GenCavMat(const int cavi, const double dis, const double EfieldScl, const double TTF_tab[],
                   const double beta_tab[], const double gamma_tab[], const double Lambda,
                   Particle &real, const double IonFys[], const double Rm, state_t::matrix_t &M,
                   const CavTLMLineType& linetab) const;

    void PropagateLongRFCav(Particle &ref);

    void InitRFCav(Particle &real, state_t::matrix_t &M, CavTLMLineType &linetab);

    void GetCavBoost(const CavDataType &CavData, Particle &state, const double IonFy0, const double fRF,
                     const double EfieldScl, double &IonFy);

    void TransFacts(const int cavilabel, double beta, const double CaviIonK, const int gaplabel, const double EfieldScl,
                    double &Ecen, double &T, double &Tp, double &S, double &Sp, double &V0);

    void TransitFacMultipole(const int cavi, const std::string &flabel, const double CaviIonK,
                             double &T, double &S);

    virtual ~ElementRFCavity() {}

    virtual void advance(StateBase& s)
    {
        state_t&  ST = static_cast<state_t&>(s);
        using namespace boost::numeric::ublas;

        // IonEk is Es + E_state; the latter is set by user.
        ST.recalc();

        if(!check_cache(ST)) {
            // need to re-calculate energy dependent terms

            recompute_matrix(ST); // updates transfer and last_Kenergy_out

            for(size_t i=0; i<last_Kenergy_in.size(); i++)
                get_misalign(ST, ST.real[i], misalign[i], misalign_inv[i]);

            ST.recalc();
        }

        // recompute_matrix only called when ST.IonEk != last_Kenergy_in.
        // Matrix elements are scaled with particle energy.

        ST.pos += length;

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
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
    }

    virtual void recompute_matrix(state_t& ST)
    {
        // Re-initialize transport matrix.


        CavTLMLineTab.resize(last_Kenergy_in.size());

        ElementRFCavity::PropagateLongRFCav(ST.ref);

        for(size_t i=0; i<last_Kenergy_in.size(); i++) {
            // TODO: 'transfer' is overwritten in InitRFCav()?
            transfer[i] = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
            transfer[i](state_t::PS_X, state_t::PS_PX) = length;
            transfer[i](state_t::PS_Y, state_t::PS_PY) = length;

            last_Kenergy_in[i] = ST.real[i].IonEk;
            // J.B. Bug in TLM.
            double SampleIonK = ST.real[i].SampleIonK;


            ElementRFCavity::InitRFCav(ST.real[i], transfer[i], CavTLMLineTab[i]);

            // J.B. Bug in TLM.
            ST.real[i].SampleIonK = SampleIonK;

            last_Kenergy_out[i] = ST.real[i].IonEk;
        }
   }

    virtual const char* type_name() const {return "rfcavity";}
};

