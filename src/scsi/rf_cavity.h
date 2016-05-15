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


void TransitFacMultipole(const int cavi, const std::string &flabel, const double IonK,
                         double &T, double &S);

void TransFacts(const int cavilabel, double beta, const int gaplabel, const double EfieldScl,
                double &Ecen, double &T, double &Tp, double &S, double &Sp, double &V0);

void EvalGapModel(const double dis, const double IonW0, Particle &real, const double IonFy0,
                  const double k, const double Lambda, const double Ecen,
                  const double T, const double S, const double Tp, const double Sp, const double V0,
                  double &IonW_f, double &IonFy_f);

double GetCavPhase(const int cavi, Particle ref, const double IonFys, const double multip);

struct ElementRFCavity : public Moment2ElementBase
{
    // Transport matrix for an RF Cavity.
    typedef Moment2ElementBase base_t;
    typedef typename base_t::state_t state_t;

    CavDataType    CavData;
    std::fstream   inf;
    CavTLMLineType CavTLMLineTab;

    ElementRFCavity(const Config& c)
        :base_t(c)
    {
        std::string cav_type = c.get<std::string>("cavtype");
        double L             = c.get<double>("L")*MtoMM;         // Convert from [m] to [mm].

        std::string CavType      = conf().get<std::string>("cavtype");
        std::string Eng_Data_Dir = conf().get<std::string>("Eng_Data_Dir", "");

        if (CavType == "0.041QWR") {
            CavData.RdData(Eng_Data_Dir+"/axisData_41.txt");
            inf.open((Eng_Data_Dir+"/Multipole41/thinlenlon_41.txt").c_str(), std::ifstream::in);
        } else if (conf().get<std::string>("cavtype") == "0.085QWR") {
            CavData.RdData(Eng_Data_Dir+"/axisData_85.txt");
            inf.open((Eng_Data_Dir+"/Multipole85/thinlenlon_85.txt").c_str(), std::ifstream::in);
        } else {
            std::ostringstream strm;
            strm << "*** InitRFCav: undef. cavity type: " << CavType << "\n";
            throw std::runtime_error(strm.str());
        }

        this->transfer_raw(state_t::PS_X, state_t::PS_PX) = L;
        this->transfer_raw(state_t::PS_Y, state_t::PS_PY) = L;
        // For total path length.
//        this->transfer(state_t::PS_S, state_t::PS_S)  = L;
    }

    void GetCavMatParams(const int cavi,
                         const double beta_tab[], const double gamma_tab[], const double IonK[]);

    void GetCavMat(const int cavi, const int cavilabel, const double Rm, Particle &real,
                   const double EfieldScl, const double IonFyi_s,
                   const double IonEk_s, const double fRF, value_mat &M);

    void GenCavMat(const int cavi, const double dis, const double EfieldScl, const double TTF_tab[],
                   const double beta_tab[], const double gamma_tab[], const double Lambda,
                   Particle &real, const double IonFys[], const double Rm, value_mat &M);

    void PropagateLongRFCav(Config &conf, Particle &ref);

    void InitRFCav(const Config &conf, Particle &real, double &accIonW,
                   double &avebeta, double &avegamma, value_mat &M);

    void GetCavBoost(const CavDataType &CavData, Particle &state, const double IonFy0, const double fRF,
                     const double EfieldScl, double &IonFy, double &accIonW);

    virtual ~ElementRFCavity() {}

    virtual void recompute_matrix(state_t& ST)
    {
        transfer = transfer_raw;

        last_Kenergy_in = ST.real.IonEk;

        this->ElementRFCavity::PropagateLongRFCav(conf(), ST.ref);

        last_Kenergy_out = ST.real.IonEk;

        // Define initial conditions.
        double accIonW, avebeta, avegamma;

        this->ElementRFCavity::InitRFCav(conf(), ST.real, accIonW, avebeta, avegamma, transfer);
   }

    virtual const char* type_name() const {return "rfcavity";}
};

