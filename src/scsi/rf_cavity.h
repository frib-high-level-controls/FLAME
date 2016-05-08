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

