#ifndef RF_CAVITY_H
#define RF_CAVITY_H

#endif // RF_CAVITY_H

#include <boost/numeric/ublas/matrix.hpp>


// Phase space dimension; including vector for orbit/1st moment.
# define PS_Dim Moment2State::maxsize // Set to 7; to include orbit.

// Mpultipole level: 0 only include focusing and defocusing effects,
//                   1 include dipole terms,
//                   2 include quadrupole terms.
const int MpoleLevel = 2;


class CavDataType {
// Cavity on-axis longitudinal electric field vs. s.
public:
    std::vector<double> s,     // s coordinate [m]
                        Elong; // Longitudinal Electric field [V/m].

    void RdData(const std::string&);
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


struct ElementRFCavity : public Moment2ElementBase
{
    // Transport matrix for an RF Cavity.
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;

    struct RawParams {
        std::string name, type;
        double length, aperature, E0;
    };
    std::vector<RawParams> lattice;

    CavDataType    CavData;
    CavTLMLineType CavTLMLineTab;
    double phi_ref;

    ElementRFCavity(const Config& c);

    void GetCavMatParams(const int cavi,
                         const double beta_tab[], const double gamma_tab[], const double IonK[],
                         CavTLMLineType& lineref) const;

    void GetCavMat(const int cavi, const int cavilabel, const double Rm, Particle &real,
                   const double EfieldScl, const double IonFyi_s,
                   const double IonEk_s, const double fRF, state_t::matrix_t &M);

    void GenCavMat(const int cavi, const double dis, const double EfieldScl, const double TTF_tab[],
                   const double beta_tab[], const double gamma_tab[], const double Lambda,
                   Particle &real, const double IonFys[], const double Rm, state_t::matrix_t &M);

    void PropagateLongRFCav(const Config &conf, Particle &ref);

    void InitRFCav(const Config &conf, Particle &real, state_t::matrix_t &M);

    void GetCavBoost(const CavDataType &CavData, Particle &state, const double IonFy0, const double fRF,
                     const double EfieldScl, double &IonFy);

    virtual ~ElementRFCavity() {}

    virtual void recompute_matrix(state_t& ST)
    {
        std::cout<<"Recompute Element "<<index<<" '"<<name<<"' ref: "<<ST.ref<<"\n";

        transfer = transfer_raw;

        last_Kenergy_in = ST.real.IonEk;

        this->ElementRFCavity::PropagateLongRFCav(conf(), ST.ref);

        last_Kenergy_out = ST.real.IonEk;

        this->ElementRFCavity::InitRFCav(conf(), ST.real, transfer);
   }

    virtual const char* type_name() const {return "rfcavity";}
};

