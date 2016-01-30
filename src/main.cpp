
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <math.h>

#include <vector>

#include "scsi/config.h"

#include <scsi/base.h>
//#include <scsi/linear.h>
#include <scsi/moment.h>
#include <scsi/state/vector.h>
#include <scsi/state/matrix.h>


extern int glps_debug;

// Phase-space units.
#if false
    // Use [m, rad, m, rad, rad, eV/u].
    #define MtoMM   1e0
    #define MeVtoeV 1e0
#else
    // Use [mm, rad, mm, rad, rad, MeV/u].
    #define MtoMM   1e3
    #define MeVtoeV 1e6
#endif

// Global constansts.

// Speed of light [m/s].
# define C0           2.99792458e8
// Atomic mass unit [eV/c^2].
# define AU           (931.49432e6/MeVtoeV)
// Vacuum permeability.
# define MU0          4e0*M_PI*1e-7
// Long. sampling frequency [Hz]; must be set to RF Cavity frequency.
# define SampleFreq   80.5e6
// Sampling distance [m].
# define SampleLambda C0/SampleFreq*MtoMM

// Phase space dimension.
# define PS_Dim       6

// Global constansts and parameters.

// Mpultipole level: 0 only include focusing and defocusing effects,
//                   1 include dipole terms,
//                   2 include quadrupole terms.
const int  MpoleLevel = 0;

// Charge stripper parameters.
const int    Stripper_n                 = 5;            // Number of charge states.
const double Stripper_IonZ              = 78e0/238e0,
             Stripper_IonMass           = 238e0,
             Stripper_IonProton         = 92e0,
             Stripper_IonChargeStates[] = {76e0/238e0, 77e0/238e0, 78e0/238e0,
                                           79e0/238e0, 80e0/238e0},
             // Energy dependance. Unit for last parmeter ???
             Stripper_E0Para[]          = {16.348e6/MeVtoeV, 1.00547, -0.10681},
             // Thickness [microns], thickness variation [%], reference energy [MeV/u].
             StripperPara[]             = {3e0, 20e0, 16.623e6/MeVtoeV};


typedef boost::numeric::ublas::vector<double> value_vec;
typedef boost::numeric::ublas::matrix<double> value_mat;


class LongTabType {
// Table for longitudinal initialization for reference particle.
public:
    std::vector<double> s,     // Longitudinal position [m].
                        Ek,    // Kinetic energy [eV/u].
                        FyAbs, // Synchrotron phase [rad].
                        Beta,  // Relativistic factor beta.
                        Gamma; // Relativistic factor gamma.

    void set(const double, const double, const double, const double, const double);
    void show(std::ostream& strm, const int) const;
};


class CavDataType {
// Cavity on-axis longitudinal electric field vs. s.
public:
    std::vector<double> s,     // s coordinate [m]
                        Elong; // Longitudinal Electric field [V/m].

    void RdData(const std::string);
    void show(std::ostream& strm, const int) const;
    void show(std::ostream& strm) const;
};


class CavTLMLineType {
public:
    std::vector<double> s;         // Longitudinal position [m].
    std::vector<std::string> Elem;
    std::vector<double> E0,
                        T,
                        S,
                        Accel;

    void set(const double, const std::string &, const double,
             const double, const double, const double);
    void show(std::ostream& strm, const int) const;
    void show(std::ostream& strm) const;
};


class RFCavType {
public:
    std::string CavType;

};


// Global variables.

const int MaxNCav = 5;

CavDataType         CavData[MaxNCav];
LongTabType         LongTab;
std::vector<double> CavPhases;
CavTLMLineType      CavTLMLineTab[MaxNCav];
RFCavType           RFCav[MaxNCav];

const std::string HomeDir = "/home/johan/tlm_workspace/TLM_JB";


void LongTabType::set(const double s, const double Ek, const double FyAbs,
                      const double Beta, const double Gamma)
{
    this->s.push_back(s); this->Ek.push_back(Ek); this->FyAbs.push_back(FyAbs);
    this->Beta.push_back(Beta); this->Gamma.push_back(Gamma);
}


void LongTabType::show(std::ostream& strm, const int k) const
{
    strm << std::scientific << std::setprecision(10)
         << std::setw(18) << this->s[k]
         << std::setw(18) << this->Ek[k]
         << std::setw(18) << this->FyAbs[k]
         << std::setw(18) << this->Beta[k]
         << std::setw(18) << this->Gamma[k] << "\n";
}


void CavDataType::RdData(const std::string FileName)
{
    std::string       line;
    int               k;
    double            s, Elong;
    std::stringstream str;
    std::fstream      inf;

    inf.open(FileName.c_str(), std::ifstream::in);
    if (!inf.is_open()) {
        std::cerr << "*** RdData: failed to open " << FileName << "\n";
        exit(1);
    }
    while (getline(inf, line)) {
        str.str(line);
        str >> s >> Elong;
        this->s.push_back(s), this->Elong.push_back(Elong);
        // Convert from [mm] to [m].
//        this->s[this->s.size()-1] /= MtoMM;
    }
    inf.close();

    if (false) {
        std::cout << "\n";
        for (k = 0; k < this->s.size(); k++)
            this->show(std::cout, k);
    }
}


void CavDataType::show(std::ostream& strm, const int k) const
{
    strm << std::scientific << std::setprecision(5)
         << std::setw(13) << this->s[k] << std::setw(13) << this->Elong[k] << "\n";
}


void CavDataType::show(std::ostream& strm) const
{
    for (unsigned int k = 0; k < this->s.size(); k++)
        this->show(strm, k);
}


void CavTLMLineType::set(const double s, const std::string &Elem, const double E0,
                         const double T, const double S, const double Accel)
{
    this->s.push_back(s); this->Elem.push_back(Elem); this->E0.push_back(E0);
    this->T.push_back(T); this->S.push_back(S); this->Accel.push_back(Accel);
}


void CavTLMLineType::show(std::ostream& strm, const int k) const
{
    strm << std::fixed << std::setprecision(5)
         << std::setw(9) << this->s[k] << std::setw(10) << this->Elem[k]
         << std::setw(9) << this->T[k] << std::setw(9) << this->S[k]
         << std::setw(9) << this->Accel[k] << "\n";
}


void CavTLMLineType::show(std::ostream& strm) const
{
    for (unsigned int k = 0; k < this->s.size(); k++)
        this->show(strm, k);
}


void PrtVec(const std::vector<double> &a)
{
    int k;

    for (k = 0; k < PS_Dim; k++)
        std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << a[k];
    std::cout << "\n";
}


void PrtMat(const value_mat &M)
{
    int j, k;

    for (j = 0; j < PS_Dim; j++) {
        for (k = 0; k < PS_Dim; k++)
            std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << M(j, k);
        std::cout << "\n";
    }
}


void PrtLat(const Machine &sim)
{
    Machine::p_elements_t::const_iterator it;
    MomentElementBase                     *ElemPtr;
    Config                                D;
    std::auto_ptr<StateBase>              state(sim.allocState(D));
    ElementVoid                           *elem;

    std::cout << "\nPrtLat:\n";

//    std::cout << "  Machine configuration:\n" << sim << "\n";

    std::cout << "  Simulation type:    " << sim.p_info.name
              << "\n  Number of elements: " << sim.p_elements.size() << "\n";

    for (it = sim.p_elements.begin(); it != sim.p_elements.end(); ++it) {
        elem = *it;

        std::cout << "\n" << elem->type_name() << " " << elem->name << "\n";

        ElemPtr = dynamic_cast<MomentElementBase *>(elem);
        assert(ElemPtr != NULL);
        PrtMat(ElemPtr->transfer);
    }
}


double calGauss(double in, const double Q_ave, const double d)
{
    // Gaussian distribution.
    return 1e0/sqrt(2e0*M_PI)/d*exp(-0.5e0*sqr(in-Q_ave)/sqr(d));
}


void calStripperCharge(const double IonProton, const double beta,
                       double &Q_ave, double &d)
{
    // Use Baron's formula for carbon foil.
    double Q_ave1, Y;

    Q_ave1 = IonProton*(1e0-exp(-83.275*(beta/pow(IonProton, 0.447))));
    Q_ave  = Q_ave1*(1e0-exp(-12.905+0.2124*IonProton-0.00122*sqr(IonProton)));
    Y      = Q_ave1/IonProton;
    d      = sqrt(Q_ave1*(0.07535+0.19*Y-0.2654*sqr(Y)));
}


void ChargeStripper(const double IonMass, const double IonProton, const double beta,
                    const int nChargeStates, const double IonChargeStates[],
                    double chargeAmount_Baron[])
{
    int    k;
    double Q_ave, d;

    calStripperCharge(IonProton, beta, Q_ave, d);
    for (k = 0; k < nChargeStates; k++)
        chargeAmount_Baron[k] = calGauss(IonChargeStates[k]*IonMass, Q_ave, d);
}


void PropagateLongStripper(const Config &conf, const int n, double &IonZ, const double IonEs,
                           double &IonW, double &SampleIonK, double &IonBeta)
{
    double IonEk, IonGamma;
    double chargeAmount_Baron[Stripper_n];

    IonZ = Stripper_IonZ;
    ChargeStripper(Stripper_IonMass, Stripper_IonProton, IonBeta,
                   Stripper_n, Stripper_IonChargeStates,
                   chargeAmount_Baron);
    // Evaluate change in reference particle energy due to stripper model energy straggling.
    IonEk      = (LongTab.Ek[n-2]-StripperPara[2])*Stripper_E0Para[1] + Stripper_E0Para[0];
    IonW       = IonEk + IonEs;
    IonGamma   = IonW/IonEs;
    IonBeta    = sqrt(1e0-1e0/sqr(IonGamma));
    SampleIonK = 2e0*M_PI/(IonBeta*SampleLambda);

    LongTab.set(LongTab.s[n-2], IonEk, LongTab.FyAbs[n-2], IonBeta, IonGamma);

    // chargeAmount = fribstripper.chargeAmount_Baron;
}


void EvalGapModel(const double dis, const double IonW0, const double IonEs, const double IonFy0,
                  const double k, const double IonZ, const double Lambda, const double Ecen,
                  const double T, const double S, const double Tp, const double Sp, const double V0,
                  double &IonW_f, double &IonFy_f)
{
    double Iongamma_f, IonBeta_f, k_f;

    IonW_f     = IonW0 + IonZ*V0*T*cos(IonFy0+k*Ecen) - IonZ*V0*S*sin(IonFy0+k*Ecen);
    Iongamma_f = IonW_f/IonEs;
    IonBeta_f  = sqrt(1e0-1e0/sqr(Iongamma_f));
    k_f        = 2e0*M_PI/(IonBeta_f*Lambda);

    IonFy_f = IonFy0 + k*Ecen + k_f*(dis-Ecen)
              + IonZ*V0*k*(Tp*sin(IonFy0+k*Ecen)+Sp*cos(IonFy0+k*Ecen))/(2e0*(IonW0-IonEs));
}


double GetCavPhase(const int cavi, const double IonEk, const double IonFys,
                   const double FyAbs, const double multip)
{
    /* If the cavity is not at full power, the method gives synchrotron
     * phase slightly different from the nominal value.                 */

    double Fyc;

    switch (cavi) {
    case 1:
        Fyc = 4.394*pow(IonEk, -0.4965) - 4.731;
        break;
    case 2:
        Fyc = 5.428*pow(IonEk, -0.5008) + 1.6;
        break;
    case 3:
        Fyc = 22.35*pow(IonEk, -0.5348) + 2.026;
        break;
    case 4:
        Fyc = 41.43*pow(IonEk, -0.5775) + 2.59839;
        break;
    case 5:
        Fyc = 5.428*pow(IonEk, -0.5008) + 1.6;
        break;
    default:
        std::cerr << "*** GetCavPhase: undef. cavity type" << "\n";
        exit(1);
    }

    return IonFys - Fyc - FyAbs*multip;
}

void GetCavBoost(const CavDataType &CavData, const double IonW0,
                 const double IonFy0, const double IonK0, const double IonZ,
                 const double IonEs, const double fRF,
                 const double EfieldScl, double &IonW, double &IonFy)
{
    int    n = CavData.s.size(),
           k;

    double dis = CavData.s[n-1] - CavData.s[0],
           dz  = dis/(n-1),
           IonLambda, IonK, IonFylast, IonGamma, IonBeta;

    IonLambda = C0/fRF*MtoMM;

    IonFy = IonFy0;
    IonK  = IonK0;
    IonW  = IonW0;
    for (k = 0; k < n-1; k++) {
        IonFylast = IonFy;
        IonFy += IonK*dz;
        IonW  += IonZ*EfieldScl*(CavData.Elong[k]+CavData.Elong[k+1])/(2e0*MeVtoeV)
                 *cos((IonFylast+IonFy)/2e0)*dz/MtoMM;
        IonGamma = IonW/IonEs;
        IonBeta = sqrt(1e0-1e0/sqr(IonGamma));
        if ((IonW-IonEs) < 0e0) {
            IonW = IonEs;
            IonBeta = 0e0;
        }
        IonK = 2e0*M_PI/(IonBeta*IonLambda);
    }
}


void PropagateLongRFCav(const Config &conf, const int n, const double IonZ, const double IonEs, double &IonW,
                        double &SampleIonK, double &IonBeta)
{
    std::string CavType;
    int         cavi;
    double      fRF, multip, caviIonK, IonFys, EfieldScl, caviFy, IonFy_i, IonFy_o;
    double      IonW_o, IonGamma;

    CavType = conf.get<std::string>("cavtype");
    if (CavType == "0.041QWR") {
        cavi = 1;
    } else if (conf.get<std::string>("cavtype") == "0.085QWR") {
        cavi = 2;
    } else {
        std::cerr << "*** PropagateLongRFCav: undef. cavity type: " << CavType << "\n";
        exit(1);
    }

    fRF       = conf.get<double>("f");
    multip    = fRF/SampleFreq;
    caviIonK  = 2e0*M_PI*fRF/(IonBeta*C0)/MtoMM;
    IonFys    = conf.get<double>("phi")*M_PI/180e0; // Synchrotron phase [rad].
    EfieldScl = conf.get<double>("scl_fac");       // Electric field scale factor.

    caviFy = GetCavPhase(cavi, IonW-IonEs, IonFys, LongTab.FyAbs[n-2], multip);

    IonFy_i = multip*LongTab.FyAbs[n-2] + caviFy;
    CavPhases.push_back(caviFy);

    if (false)
        std::cout << std::scientific << std::setprecision(10)
                  << "CavPhase: " << std::setw(3) << CavPhases.size()
                  << std::setw(18) << CavPhases[CavPhases.size()-1] << "\n";

    // For the reference particle, evaluate the change of:
    // kinetic energy, absolute phase, beta, and gamma.
    GetCavBoost(CavData[cavi-1], IonW, IonFy_i, caviIonK, IonZ,
                IonEs, fRF, EfieldScl, IonW_o, IonFy_o);
    IonW       = IonW_o;
    IonGamma   = IonW/IonEs;
    IonBeta    = sqrt(1e0-1e0/sqr(IonGamma));
    SampleIonK = 2e0*M_PI/(IonBeta*SampleLambda);

    LongTab.set(LongTab.s[n-2]+conf.get<double>("L"), IonW-IonEs,
            LongTab.FyAbs[n-2]+(IonFy_o-IonFy_i)/multip, IonBeta, IonGamma);
}


void InitLong(const Machine &sim)
{
    /* Longitudinal initialization for reference particle.
     * Evaluate beam energy and cavity loaded phase along the lattice. */
    int                                   n;
    double                                IonGamma, IonBeta, SampleIonK;
    double                                IonW, IonZ, IonEs;
    Machine::p_elements_t::const_iterator it;

    Config                   D;
    std::auto_ptr<StateBase> state(sim.allocState(D));
    // Propagate through first element.
    sim.propagate(state.get(), 0, 1);

    IonZ  = state->IonZ;
    IonEs = state->IonEs/MeVtoeV;
    IonW  = state->IonW/MeVtoeV;

    IonGamma   = IonW/IonEs;
    IonBeta    = sqrt(1e0-1e0/sqr(IonGamma));
    SampleIonK = 2e0*M_PI/(IonBeta*SampleLambda);

    std::cout << "\n" << "InitLong:" << "\n";
    std::cout << std::fixed << std::setprecision(5)
              << "  IonEs [Mev/u] = " << IonEs << ", IonEk [Mev/u] = " << state->IonEk/MeVtoeV
              << ", IonW [Mev/u] = " << IonW << "\n";

    n = 1;
    LongTab.set(0e0, state->IonEk/MeVtoeV, 0e0, IonBeta, IonGamma);

    it = sim.p_elements.begin();
    // Skip over state.
    it++;
    do {
        ElementVoid*  elem   = *it;
        const Config& conf   = elem->conf();
        std::string   t_name = elem->type_name(); // C string -> C++ string.

        if (t_name == "marker") {
        } else if (t_name == "drift") {
            n++;
            LongTab.set(LongTab.s[n-2]+conf.get<double>("L"), LongTab.Ek[n-2],
                    LongTab.FyAbs[n-2]+SampleIonK*conf.get<double>("L")*MtoMM,
                    LongTab.Beta[n-2], LongTab.Gamma[n-2]);
        } else if (t_name == "sbend") {
            n++;
            LongTab.set(LongTab.s[n-2]+conf.get<double>("L"), LongTab.Ek[n-2],
                    LongTab.FyAbs[n-2]+SampleIonK*conf.get<double>("L")*MtoMM,
                    LongTab.Beta[n-2], LongTab.Gamma[n-2]);
        } else if (t_name == "quadrupole") {
            n++;
            LongTab.set(LongTab.s[n-2]+conf.get<double>("L"), LongTab.Ek[n-2],
                    LongTab.FyAbs[n-2]+SampleIonK*conf.get<double>("L")*MtoMM,
                    LongTab.Beta[n-2], LongTab.Gamma[n-2]);
        } else if (t_name == "solenoid") {
            n++;
            LongTab.set(LongTab.s[n-2]+conf.get<double>("L"), LongTab.Ek[n-2],
                    LongTab.FyAbs[n-2]+SampleIonK*conf.get<double>("L")*MtoMM,
                    LongTab.Beta[n-2], LongTab.Gamma[n-2]);
        } else if (t_name == "rfcavity") {
            n++;
            PropagateLongRFCav(conf, n, IonZ, IonEs, IonW, SampleIonK, IonBeta);
        } else if (t_name == "stripper") {
            // Evaluate change in reference particle energy and multi-charge states, and charge.
            n++;
            PropagateLongStripper(conf, n, IonZ, IonEs, IonW, SampleIonK, IonBeta);
        }

        if (false) {
            std::cout << std::setw(10) << std::left << t_name << std::setw(25)
                      << elem->name << std::internal;
            LongTab.show(std::cout, n-1);
        }
        it++;
    } while (it != sim.p_elements.end());
}


double PwrSeries(const double beta,
                 const double a0, const double a1, const double a2, const double a3,
                 const double a4, const double a5, const double a6, const double a7)
{
    int    k;
    double f;

    const int    n   = 8;
    const double a[] = {a0, a1, a2, a3, a4, a5, a6, a7};

    f = a[0];
    for (k = 1; k < n; k++)
        f += a[k]*pow(beta, k);

    return f;
}


void TransFacts(const int cavilabel, double beta, const int gaplabel, const double EfieldScl,
                double &Ecen, double &T, double &Tp, double &S, double &Sp, double &V0)
{
    // Evaluate Electric field center, transit factors [T, T', S, S'] and cavity field.
    std::vector<double> vec;

    switch (cavilabel) {
    case 41:
        if (beta < 0.025 || beta > 0.08) {
            std::cerr << "*** GetTransitFac: beta out of Range " << beta << "\n";
            exit(1);
        }
        switch (gaplabel) {
        case 0:
            // One gap evaluation.
            Ecen = 120.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = PwrSeries(beta, -4.109, 399.9, -1.269e4, 1.991e5, -1.569e6, 4.957e6, 0.0, 0.0);
            Sp   = PwrSeries(beta, 61.98, -1.073e4, 4.841e5, 9.284e6, 8.379e7, -2.926e8, 0.0, 0.0);
            V0   = 0.98477*EfieldScl;
            break;
        case 1:
            // Two gap calculation, first gap.
            Ecen = 0.0006384*pow(beta, -1.884) + 86.69;
            T    = PwrSeries(beta, 0.9232, -123.2, 3570, -5.476e4, 4.316e5, -1.377e6, 0.0, 0.0);
            Tp   = PwrSeries(beta, 1.699, 924.7, -4.062e4, 7.528e5, -6.631e6, 2.277e7, 0.0, 0.0);
            S    = 0.0;
            Sp   = PwrSeries(beta, -1.571, 25.59, 806.6, -2.98e4, 3.385e5, -1.335e6, 0.0, 0.0);
            V0   = 0.492385*EfieldScl;
            break;
        case 2:
            // Two gap calculation, second gap.
            Ecen = -0.0006384*pow(beta, -1.884) + 33.31;
            T    = PwrSeries(beta, -0.9232, 123.2, -3570, 5.476e4, -4.316e5, 1.377e6, 0.0, 0.0);
            Tp   = PwrSeries(beta, -1.699, -924.7, 4.062e4, -7.528e5, 6.631e6, -2.277e7, 0.0, 0.0);
            S    = 0.0;
            Sp    = PwrSeries(beta, -1.571, 25.59, 806.6, -2.98e4, 3.385e5, -1.335e6, 0.0, 0.0);
            V0   = 0.492385*EfieldScl;
            break;
        default:
            std::cerr << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            exit(1);
        }
        break;
    case 85:
        if (beta < 0.05 || beta > 0.25) {
            std::cerr << "*** GetTransitFac: beta out of range " << beta << "\n";
            exit(1);
        }
        switch (gaplabel) {
          case 0:
            Ecen = 150.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = PwrSeries(beta, -6.811, 343.9, -6385, 6.477e4, -3.914e5, 1.407e6, -2.781e6, 2.326e6);
            Sp   = PwrSeries(beta, 162.7, -1.631e4, 4.315e5, -5.344e6, 3.691e7, -1.462e8, 3.109e8, -2.755e8);
            V0   = 1.967715*EfieldScl;
            break;
        case 1:
            Ecen = 0.0002838*pow(beta, -2.13) + 76.5;
            T    = 0.0009467*pow(beta, -1.855) - 1.002;
            Tp   = PwrSeries(beta, 24.44, -334, 2468, -1.017e4, 2.195e4, -1.928e4, 0.0, 0.0);
            S    = 0.0;
            Sp   = -0.0009751*pow(beta, -1.898) + 0.001568;
            V0   = 0.9838574*EfieldScl;
            break;
        case 2:
            Ecen = -0.0002838*pow(beta, -2.13) + 73.5;
            T    = -0.0009467*pow(beta, -1.855) + 1.002;
            Tp   = PwrSeries(beta,  24.44, 334,  2468, 1.017e4, -2.195e4, 1.928e4, 0.0, 0.0);
            S    = 0.0;
            Sp   = -0.0009751*pow(beta, -1.898) + 0.001568;
            V0   = 0.9838574*EfieldScl;
            break;
        default:
            std::cerr << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            exit(1);
        }
        break;
    default:
        std::cerr << "*** GetTransitFac: undef. cavity type" << "\n";
        exit(1);
    }

    // Convert from [mm] to [m].
//    Ecen /= MtoMM;
}


double PwrSeries(const double beta,
                 const double a0, const double a1, const double a2, const double a3,
                 const double a4, const double a5, const double a6, const double a7,
                 const double a8, const double a9)
{
    int    k;
    double f;

    const int    n   = 10;
    const double a[] = {a0, a1, a2, a3, a4, a5, a6, a7, a8, a9};

    f = a[0];
    for (k = 1; k < n; k++)
        f += a[k]*pow(beta, k);

    return f;
}


void TransitFacMultipole(const int cavi, const std::string &flabel, const double IonK,
                         double &T, double &S)
{

    if ((cavi == 1) && (IonK < 0.025 || IonK > 0.055)) {
        std::cerr << "*** TransitFacMultipole: IonK out of Range" << "\n";
        exit(1);
    } else if ((cavi == 2) && (IonK < 0.006 || IonK > 0.035)) {
        std::cerr << "*** TransitFacMultipole: IonK out of Range" << "\n";
    }

    if (flabel == "CaviMlp_EFocus1") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 1.256386e+02, -3.108322e+04, 3.354464e+06, -2.089452e+08, 8.280687e+09, -2.165867e+11,
                          3.739846e+12, -4.112154e+13, 2.613462e14, -7.316972e14);
            S = PwrSeries(IonK, 1.394183e+02, -3.299673e+04, 3.438044e+06, -2.070369e+08, 7.942886e+09, -2.013750e+11,
                         3.374738e+12, -3.605780e+13, 2.229446e+14, -6.079177e+14);
            break;
        case 2:
            T = PwrSeries(IonK, -9.450041e-01, -3.641390e+01, 9.926186e+03, -1.449193e+06, 1.281752e+08, -7.150297e+09,
                          2.534164e+11, -5.535252e+12, 6.794778e+13, -3.586197e+14);
            S = PwrSeries(IonK, 9.928055e-02, -5.545119e+01, 1.280168e+04, -1.636888e+06, 1.279801e+08, -6.379800e+09,
                          2.036575e+11, -4.029152e+12, 4.496323e+13, -2.161712e+14);
            break;
        }
    } else if (flabel == "CaviMlp_EFocus2") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 1.038803e+00, -9.121320e+00, 8.943931e+02, -5.619149e+04, 2.132552e+06, -5.330725e+07,
                          8.799404e+08, -9.246033e+09, 5.612073e+10, -1.499544e+11);
            S = PwrSeries(IonK, 1.305154e-02, -2.585211e+00, 2.696971e+02, -1.488249e+04, 5.095765e+05, -1.154148e+07,
                          1.714580e+08, -1.604935e+09, 8.570757e+09, -1.983302e+10);
            break;
        case 2:
            T = PwrSeries(IonK, 9.989307e-01, 7.299233e-01, -2.932580e+02, 3.052166e+04, -2.753614e+06, 1.570331e+08,
                          -5.677804e+09, 1.265012e+11, -1.584238e+12, 8.533351e+12);
            S = PwrSeries(IonK, -3.040839e-03, 2.016667e+00, -4.313590e+02, 5.855139e+04, -4.873584e+06, 2.605444e+08,
                          -8.968899e+09, 1.923697e+11, -2.339920e+12, 1.233014e+13);
            break;
        }
    } else if (flabel == "CaviMlp_EDipole") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, -1.005885e+00, 1.526489e+00, -1.047651e+02, 1.125013e+04, -4.669147e+05, 1.255841e+07,
                          -2.237287e+08, 2.535541e+09, -1.656906e+10, 4.758398e+10);
            S = PwrSeries(IonK, -2.586200e-02, 5.884367e+00, -6.407538e+02, 3.888964e+04, -1.488484e+06, 3.782592e+07,
                          -6.361033e+08, 6.817810e+09, -4.227114e+10, 1.155597e+11);
            break;
        case 2:
            T = PwrSeries(IonK, -9.999028e-01, -6.783669e-02, 1.415756e+02, -2.950990e+03, 2.640980e+05, -1.570742e+07,
                          5.770450e+08, -1.303686e+10, 1.654958e+11, -9.030017e+11);
            S = PwrSeries(IonK, 2.108581e-04, -3.700608e-01, 2.851611e+01, -3.502994e+03, 2.983061e+05, -1.522679e+07,
                          4.958029e+08, -1.002040e+10, 1.142835e+11, -5.617061e+11);
            break;
        }
    } else if (flabel == "CaviMlp_EQuad") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 1.038941e+00, -9.238897e+00, 9.127945e+02, -5.779110e+04, 2.206120e+06, -5.544764e+07,
                          9.192347e+08, -9.691159e+09, 5.896915e+10, -1.578312e+11);
            S = PwrSeries(IonK, 1.248096e-01, -2.923507e+01, 3.069331e+03, -1.848380e+05, 7.094882e+06, -1.801113e+08,
                          3.024208e+09, -3.239241e+10, 2.008767e+11, -5.496217e+11);
            break;
        case 2:
            T = PwrSeries(IonK, 1.000003e+00, -1.015639e-03, -1.215634e+02, 1.720764e+01, 3.921401e+03, 2.674841e+05,
                          -1.236263e+07, 3.128128e+08, -4.385795e+09, 2.594631e+10);
            S = PwrSeries(IonK, -1.756250e-05, 2.603597e-01, -2.551122e+00, -4.840638e+01, -2.870201e+04, 1.552398e+06,
                          -5.135200e+07, 1.075958e+09, -1.277425e+10, 6.540748e+10);
            break;
        }
    } else if (flabel == "CaviMlp_HMono") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 1.703336e+00, -1.671357e+02, 1.697657e+04, -9.843253e+05, 3.518178e+07, -8.043084e+08,
                          1.165760e+10, -1.014721e+11, 4.632851e+11, -7.604796e+11);
            S = PwrSeries(IonK, 1.452657e+01, -3.409550e+03, 3.524921e+05, -2.106663e+07, 8.022856e+08, 3.360597e+11,
                          -2.019481e+10, -3.565836e+12, 2.189668e+13, -5.930241e+13);
            break;
        case 2:
            T = PwrSeries(IonK, 1.003228e+00, -1.783406e+00, 1.765330e+02, -5.326467e+04, 4.242623e+06, -2.139672e+08,
                          6.970488e+09, -1.411958e+11, 1.617248e+12, -8.000662e+12);
            S = PwrSeries(IonK, -1.581533e-03, 1.277444e+00, -2.742508e+02, 3.966879e+04, -3.513478e+06, 1.962939e+08,
                          -6.991916e+09, 1.539708e+11, -1.910236e+12, 1.021016e+13);
            break;
        }
    } else if (flabel == "CaviMlp_HDipole") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 6.853803e-01, 7.075414e+01, -7.117391e+03, 3.985674e+05, -1.442888e+07, 3.446369e+08,
                          -5.420826e+09, 5.414689e+10, -3.116216e+11, 7.869717e+11);
            S = PwrSeries(IonK, 1.021102e+00, -2.441117e+02, 2.575274e+04, -1.569273e+06, 6.090118e+07, -1.562284e+09,
                          2.649289e+10, -2.864139e+11, 1.791634e+12, -4.941947e+12);
            break;
        case 2:
            T = PwrSeries(IonK, 1.014129e+00, -8.016304e+00, 1.631339e+03, -2.561826e+05, 2.115355e+07, -1.118723e+09,
                          3.821029e+10, -8.140248e+11, 9.839613e+12, -5.154137e+13);
            S = PwrSeries(IonK, -4.688714e-03, 3.299051e+00, -8.101936e+02, 1.163814e+05, -1.017331e+07, 5.607330e+08,
                          -1.967300e+10, 4.261388e+11, -5.194592e+12, 2.725370e+13);
            break;
        }
    } else if (flabel == "CaviMlp_HQuad") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, -1.997432e+00, 2.439177e+02, -2.613724e+04, 1.627837e+06, -6.429625e+07, 1.676173e+09,
                          -2.885455e+10, 3.163675e+11, -2.005326e+12, 5.600545e+12);
            S = PwrSeries(IonK, -2.470704e+00, 5.862902e+02, -6.135071e+04, 3.711527e+06, -1.431267e+08, 3.649414e+09,
                          -6.153570e+10, 6.617859e+11, -4.119861e+12, 1.131390e+13);
            break;
        case 2:
            T = PwrSeries(IonK, -1.000925e+00, 5.170302e-01, 9.311761e+01, 1.591517e+04, -1.302247e+06, 6.647808e+07,
                          -2.215417e+09, 4.603390e+10, -5.420873e+11, 2.764042e+12);
            S = PwrSeries(IonK, 3.119419e-04, -4.540868e-01, 5.433028e+01, -7.571946e+03, 6.792565e+05, -3.728390e+07,
                          1.299263e+09, -2.793705e+10, 3.377097e+11, -1.755126e+12);
            break;
        }
    } else {
        std::cerr << "*** TransitFacMultipole: undef. multipole type " << flabel << "\n";
        exit(1);
    }
}


void GetCavMatParams(const int cavi, const std::string &CavTLMLine,
                     const double beta_tab[], const double gamma_tab[], const double IonK[])
{
    // Evaluate time transit factors and acceleration.
    std::string       line, Elem, Name;
    double            s, Length, Aper, E0, T, S, Accel;
    std::stringstream str;
    std::fstream      inf;

    inf.open((HomeDir+CavTLMLine).c_str(), std::ifstream::in);
    if (!inf.is_open()) {
        std::cerr << "*** GetCavMatParams: failed to open " << CavTLMLine << "\n";
        exit(1);
    }

    s = CavData[cavi-1].s[0];
    while (getline(inf, line)) {
        T = 0e0, S = 0e0, Accel = 0e0;
        if (line[0] == '%') {
            // Comment.
        } else {
            str.str(line);
            str >> Elem >> Name >> Length >> Aper;

            s += Length;

            if ((Elem != "drift") && (Elem != "AccGap"))
                str >> E0;
            else
                E0 = 0e0;

            if (Elem == "drift") {
            } else if (Elem == "EFocus1") {
                if (s < 0e0) {
                    // First gap. By reflection 1st Gap EFocus1 is 2nd gap EFocus2.
                    TransitFacMultipole(cavi, "CaviMlp_EFocus2", IonK[0], T, S);
                    // First gap *1, transverse E field the same.
                    S = -S;
                } else {
                    // Second gap.
                    TransitFacMultipole(cavi, "CaviMlp_EFocus1", IonK[1], T, S);
                }
            } else if (Elem == "EFocus2") {
                if (s < 0e0) {
                    // First gap.
                    TransitFacMultipole(cavi, "CaviMlp_EFocus1", IonK[0], T, S);
                    S = -S;
                } else {
                    // Second gap.
                    TransitFacMultipole(cavi, "CaviMlp_EFocus2", IonK[1], T, S);
                }
            } else if (Elem == "EDipole") {
                if (MpoleLevel >= 1) {
                    if (s < 0e0) {
                        TransitFacMultipole(cavi, "CaviMlp_EDipole", IonK[0], T, S);
                        // First gap *1, transverse E field the same.
                        S = -S;
                    } else {
                        // Second gap.
                        TransitFacMultipole(cavi, "CaviMlp_EDipole", IonK[1], T, S);
                    }
                }
            } else if (Elem == "EQuad") {
                if (MpoleLevel >= 2) {
                    if (s < 0e0) {
                        // First gap.
                        TransitFacMultipole(cavi, "CaviMlp_EQuad", IonK[0], T, S);
                        S = -S;
                    } else {
                        // Second Gap
                        TransitFacMultipole(cavi, "CaviMlp_EQuad", IonK[1], T, S);
                    }
                }
            } else if (Elem == "HMono") {
                if (MpoleLevel >= 2) {
                    if (s < 0e0) {
                        // First gap.
                        TransitFacMultipole(cavi, "CaviMlp_HMono", IonK[0], T, S);
                        T = -T;
                    } else {
                        // Second Gap
                        TransitFacMultipole(cavi, "CaviMlp_HMono", IonK[1], T, S);
                    }
                }
            } else if (Elem == "HDipole") {
                if (MpoleLevel >= 1) {
                    if (s < 0e0) {
                        // First gap.
                        TransitFacMultipole(cavi, "CaviMlp_HDipole", IonK[0], T, S);
                        T = -T;
                    }  else {
                        // Second gap.
                        TransitFacMultipole(cavi, "CaviMlp_HDipole", IonK[1], T, S);
                    }
                }
            } else if (Elem == "HQuad") {
                if (MpoleLevel >= 2) {
                    if (s < 0e0) {
                        // First gap.
                        TransitFacMultipole(cavi, "CaviMlp_HQuad", IonK[0], T, S);
                        T = -T;
                    } else {
                        // Second gap.
                        TransitFacMultipole(cavi, "CaviMlp_HQuad", IonK[1], T, S);
                    }
                }
            } else if (Elem == "AccGap") {
                if (s < 0e0) {
                    // First gap.
                    Accel = (beta_tab[0]*gamma_tab[0])/((beta_tab[1]*gamma_tab[1]));
                } else {
                    // Second gap.
                    Accel = (beta_tab[1]*gamma_tab[1])/((beta_tab[2]*gamma_tab[2]));
                }
            } else {
                std::cerr << "*** GetCavMatParams: undef. multipole element " << Elem << "\n";
                exit(1);
            }

            CavTLMLineTab[cavi-1].set(s, Elem, E0, T, S, Accel);
        }
    }

    if (false) {
        std::cout << "\n";
        CavTLMLineTab[cavi-1].show(std::cout);
    }

    inf.close();
}


void GenCavMat(const int cavi, const double dis, const double EfieldScl, const double TTF_tab[],
               const double beta_tab[], const double gamma_tab[], const double Lambda,
               const double IonZ, const double IonEs, const double IonFys[], const std::string &CavTLMLine,
               const double Rm, value_mat &M)
{
    /* RF cavity model, transverse only defocusing.
     * 2-gap matrix model.                                            */

    std::string       line, Elem, Name;
    int               seg, n;
    double            Length, Aper, Efield, s, k_s[3];
    double            Ecens[2], Ts[2], Ss[2], V0s[2], ks[2], L1, L2, L3;
    double            beta, gamma, kfac, V0, T, S, kfdx, kfdy, dpy, Accel, IonFy;
    value_mat         Idmat, Mlon_L1, Mlon_K1, Mlon_L2;
    value_mat         Mlon_K2, Mlon_L3, Mlon, Mtrans, Mprob;
    std::stringstream str;
    std::fstream      inf;

    const double IonA = 1e0;

    using boost::numeric::ublas::prod;

    Idmat = boost::numeric::ublas::identity_matrix<double>(6);

    inf.open((HomeDir+CavTLMLine).c_str(), std::ifstream::in);
    if (!inf.is_open()) {
        std::cerr << "*** GenCavMat: failed to open " << CavTLMLine << "\n";
        exit(1);
    }

    k_s[0] = 2e0*M_PI/(beta_tab[0]*Lambda);
    k_s[1] = 2e0*M_PI/(beta_tab[1]*Lambda);
    k_s[2] = 2e0*M_PI/(beta_tab[2]*Lambda);

    // Longitudinal model: Drift-Kick-Drift, dis: total lenghth centered at 0,
    // Ecens[0] & Ecens[1]: Electric Center position where accel kick applies, Ecens[0] < 0
    // TTFtab: 2*6 vector, Ecens, T Tp S Sp, V0;

    Ecens[0] = TTF_tab[0];
    Ts[0]    = TTF_tab[1];
    Ss[0]    = TTF_tab[3];
    V0s[0]   = TTF_tab[5];
    ks[0]    = 0.5*(k_s[0]+k_s[1]);
    L1       = dis + Ecens[0];       //try change dis/2 to dis 14/12/12

    Mlon_L1 = Idmat;
    Mlon_K1 = Idmat;
    // Pay attention, original is -
    Mlon_L1(4, 5) = -2e0*M_PI/Lambda*(1e0/cube(beta_tab[0]*gamma_tab[0])/IonEs*L1);
    // Pay attention, original is -k1-k2
    Mlon_K1(5, 4) = -IonZ*V0s[0]*Ts[0]*sin(IonFys[0]+ks[0]*L1)-IonZ*V0s[0]*Ss[0]*cos(IonFys[0]+ks[0]*L1);

    Ecens[1] = TTF_tab[6];
    Ts[1]    = TTF_tab[7];
    Ss[1]    = TTF_tab[9];
    V0s[1]   = TTF_tab[11];
    ks[1]    = 0.5*(k_s[1]+k_s[2]);
    L2      = Ecens[1] - Ecens[0];

    Mlon_L2 = Idmat;
    Mlon_K2 = Idmat;

    Mlon_L2(4, 5) = -2e0*M_PI/Lambda*(1e0/cube(beta_tab[1]*gamma_tab[1])/IonEs*L2); //Problem is Here!!
    Mlon_K2(5, 4) = -IonZ*V0s[1]*Ts[1]*sin(IonFys[1]+ks[1]*Ecens[1])-IonZ*V0s[1]*Ss[1]*cos(IonFys[1]+ks[1]*Ecens[1]);

    L3 = dis - Ecens[1]; //try change dis/2 to dis 14/12/12

    Mlon_L3       = Idmat;
    Mlon_L3(4, 5) = -2e0*M_PI/Lambda*(1e0/cube(beta_tab[2]*gamma_tab[2])/IonEs*L3);

    Mlon = Idmat;
    Mlon = prod(Mlon_K1, Mlon_L1);
    Mlon = prod(Mlon_L2, Mlon);
    Mlon = prod(Mlon_K2, Mlon);
    Mlon = prod(Mlon_L3, Mlon);

    // Transverse model
    // Drift-FD-Drift-LongiKick-Drift-FD-Drift-0-Drift-FD-Drift-LongiKick-Drift-FD-Drift

    seg    = 0;

    Mtrans = Idmat;
    Mprob  = Idmat;

    beta   = beta_tab[0];
    gamma  = gamma_tab[0];
    IonFy  = IonFys[0];
    kfac   = k_s[0];

    V0 = 0e0, T = 0e0, S = 0e0, kfdx = 0e0, kfdy = 0e0, dpy = 0e0;

    s = CavData[cavi-1].s[0];
    n = 0;
    while (getline(inf, line)) {
        if (line[0] == '%') {
            // Comment.
        } else {
            n++;
            str.str(line);
            str >> Elem >> Name >> Length >> Aper;

            s += Length;

            if ((Elem != "drift") && (Elem != "AccGap"))
                str >> Efield;
            else
                Efield = 0e0;

            if (false)
                printf("%9.5f %8s %8s %9.5f %9.5f %9.5f\n",
                       s, Elem.c_str(), Name.c_str(), Length, Aper, Efield);

            Mprob = Idmat;
            if (Elem == "drift") {
                IonFy = IonFy + kfac*Length;

                Mprob(0, 1) = Length;
                Mprob(2, 3) = Length;
                Mtrans      = prod(Mprob, Mtrans);
            } else if (Elem == "EFocus1") {
                V0   = CavTLMLineTab[cavi-1].E0[n-1]*EfieldScl;
                T    = CavTLMLineTab[cavi-1].T[n-1];
                S    = CavTLMLineTab[cavi-1].S[n-1];
                kfdx = IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy))/Rm;
                kfdy = kfdx;

                Mprob(1, 0) = kfdx;
                Mprob(3, 2) = kfdy;
                Mtrans      = prod(Mprob, Mtrans);
            } else if (Elem == "EFocus2") {
                V0   = CavTLMLineTab[cavi-1].E0[n-1]*EfieldScl;
                T    = CavTLMLineTab[cavi-1].T[n-1];
                S    = CavTLMLineTab[cavi-1].S[n-1];
                kfdx = IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy))/Rm;
                kfdy = kfdx;

                Mprob(1, 0) = kfdx;
                Mprob(3, 2) = kfdy;
                Mtrans      = prod(Mprob, Mtrans);
            } else if (Elem == "EDipole") {
                if (MpoleLevel >= 1) {
                    V0  = CavTLMLineTab[cavi-1].E0[n-1]*EfieldScl;
                    T   = CavTLMLineTab[cavi-1].T[n-1];
                    S   = CavTLMLineTab[cavi-1].S[n-1];
                    dpy = IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy));

                    Mprob(3, 6) = dpy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "EQuad") {
                if (MpoleLevel >= 2) {
                    V0   = CavTLMLineTab[cavi-1].E0[n-1]*EfieldScl;
                    T    = CavTLMLineTab[cavi-1].T[n-1];
                    S    = CavTLMLineTab[cavi-1].S[n-1];
                    kfdx =  IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy))/Rm;
                    kfdy = -kfdx;

                    Mprob(1, 0) = kfdx;
                    Mprob(3, 2) = kfdy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "HMono") {
                if (MpoleLevel >= 2) {
                    V0   = CavTLMLineTab[cavi-1].E0[n-1]*EfieldScl;
                    T    = CavTLMLineTab[cavi-1].T[n-1];
                    S    = CavTLMLineTab[cavi-1].S[n-1];
                    kfdx = -MU0*C0*IonZ*V0/beta/gamma/IonA/AU*(T*cos(IonFy+M_PI/2e0)-S*sin(IonFy+M_PI/2e0))/Rm;
                    kfdy = kfdx;

                    Mprob(1, 0) = kfdx;
                    Mprob(3, 2) = kfdy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "HDipole") {
                if (MpoleLevel >= 1) {
                    V0  = CavTLMLineTab[cavi-1].E0[n-1]*EfieldScl;
                    T   = CavTLMLineTab[cavi-1].T[n-1];
                    S   = CavTLMLineTab[cavi-1].S[n-1];
                    dpy = -MU0*C0*IonZ*V0/beta/gamma/IonA/AU*(T*cos(IonFy+M_PI/2e0)-S*sin(IonFy+M_PI/2e0));

                    Mprob(3, 6) = dpy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "HQuad") {
                if (MpoleLevel >= 2) {
                    if (s < 0e0) {
                        // First gap.
                        beta  = (beta_tab[0]+beta_tab[1])/2e0;
                        gamma = (gamma_tab[0]+gamma_tab[1])/2e0;
                    } else {
                        beta  = (beta_tab[1]+beta_tab[2])/2e0;
                        gamma = (gamma_tab[1]+gamma_tab[2])/2e0;
                    }
                    V0   = CavTLMLineTab[cavi-1].E0[n-1]*EfieldScl;
                    T    = CavTLMLineTab[cavi-1].T[n-1];
                    S    = CavTLMLineTab[cavi-1].S[n-1];
                    kfdx = -MU0*C0*IonZ*V0/beta/gamma/IonA/AU*(T*cos(IonFy+M_PI/2e0)-S*sin(IonFy+M_PI/2e0))/Rm;
                    kfdy =  kfdx;

                    Mprob(1, 0) = kfdx;
                    Mprob(3, 2) = kfdy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "AccGap") {
                //IonFy = IonFy + IonZ*V0s[0]*kfac*(TTF_tab[2]*sin(IonFy)
                //        + TTF_tab[4]*cos(IonFy))/2/((gamma-1)*IonEs); //TTF_tab[2]~Tp
                seg    = seg + 1;
                beta   = beta_tab[seg];
                gamma  = gamma_tab[seg];
                kfac   = 2e0*M_PI/(beta*Lambda);
                Accel  = CavTLMLineTab[cavi-1].Accel[n-1];

                Mprob(1, 1) = Accel;
                Mprob(3, 3) = Accel;
                Mtrans      = prod(Mprob, Mtrans);
            } else {
                std::cerr << "*** GenCavMat: undef. multipole type " << Elem << "\n";
                exit(1);
            }
        }
    }

    inf.close();

    M = Mtrans;

    M(4, 4) = Mlon(4, 4);
    M(4, 5) = Mlon(4, 5);
    M(5, 4) = Mlon(5, 4);
    M(5, 5) = Mlon(5, 5);
}


void GetCavMat(const int cavi, const int cavilabel, const double Rm,
               const std::string &CavTLMLine,
               const double IonZ, const double IonEs, const double EfieldScl, const double IonFyi_s,
               const double IonEk_s, const double fRF, value_mat &M)
{
    int    n;
    double IonLambda, Ecen[2], T[2], Tp[2], S[2], Sp[2], V0[2];
    double dis, IonW_s[3], IonFy_s[3], gamma_s[3], beta_s[3], IonK_s[3];
    double IonK[2];

    IonLambda  = C0/fRF*MtoMM;

    IonW_s[0]  = IonEk_s + IonEs;
    IonFy_s[0] = IonFyi_s;
    gamma_s[0] = IonW_s[0]/IonEs;
    beta_s[0]  = sqrt(1e0-1e0/sqr(gamma_s[0]));
    IonK_s[0]  = 2e0*M_PI/(beta_s[0]*IonLambda);

    n   = CavData[cavi-1].s.size();
    dis = (CavData[cavi-1].s[n-1]-CavData[cavi-1].s[0])/2e0;

    TransFacts(cavilabel, beta_s[0], 1, EfieldScl, Ecen[0], T[0], Tp[0], S[0], Sp[0], V0[0]);
    EvalGapModel(dis, IonW_s[0], IonEs, IonFy_s[0], IonK_s[0], IonZ, IonLambda,
                Ecen[0], T[0], S[0], Tp[0], Sp[0], V0[0], IonW_s[1], IonFy_s[1]);
    gamma_s[1] = IonW_s[1]/IonEs;
    beta_s[1]  = sqrt(1e0-1e0/sqr(gamma_s[1]));
    IonK_s[1]  = 2e0*M_PI/(beta_s[1]*IonLambda);

    TransFacts(cavilabel, beta_s[1], 2, EfieldScl, Ecen[1], T[1], Tp[1], S[1], Sp[1], V0[1]);
    EvalGapModel(dis, IonW_s[1], IonEs, IonFy_s[1], IonK_s[1], IonZ, IonLambda,
                Ecen[1], T[1], S[1], Tp[1], Sp[1], V0[1], IonW_s[2], IonFy_s[2]);
    gamma_s[2] = IonW_s[2]/IonEs;
    beta_s[2]  = sqrt(1e0-1e0/sqr(gamma_s[2]));
    IonK_s[2]  = 2e0*M_PI/(beta_s[2]*IonLambda);

    Ecen[0] = Ecen[0] - dis;

    double TTF_tab[] = {Ecen[0], T[0], Tp[0], S[0], Sp[0], V0[0], Ecen[1], T[1], Tp[1], S[1], Sp[1], V0[1]};
    IonK[0] = (IonK_s[0]+IonK_s[1])/2e0;
    IonK[1] = (IonK_s[1]+IonK_s[2])/2e0;

    if (false) {
        printf("\n GetCavMat:\n");
        printf("IonK  : %15.10f %15.10f %15.10f\n", IonK_s[0], IonK_s[1], IonK_s[2]);
        printf("IonK  : %15.10f %15.10f\n", IonK[0], IonK[1]);
        printf("beta  : %15.10f %15.10f %15.10f\n", beta_s[0], beta_s[1], beta_s[2]);
        printf("gamma : %15.10f %15.10f %15.10f\n", gamma_s[0], gamma_s[1], gamma_s[2]);
        printf("Ecen  : %15.10f %15.10f\n", Ecen[0], Ecen[1]);
        printf("T     : %15.10f %15.10f\n", T[0], T[1]);
        printf("Tp    : %15.10f %15.10f\n", Tp[0], Tp[1]);
        printf("S     : %15.10f %15.10f\n", S[0], S[1]);
        printf("Sp    : %15.10f %15.10f\n", Sp[0], Sp[1]);
        printf("V0    : %15.10f %15.10f\n", V0[0], V0[1]);
    }

    GetCavMatParams(cavi, CavTLMLine, beta_s, gamma_s, IonK);
    GenCavMat(cavi, dis, EfieldScl, TTF_tab, beta_s, gamma_s, IonLambda, IonZ, IonEs, IonFy_s, CavTLMLine, Rm, M);

    std::cout << "\n";
    PrtMat(M);
}


void InitRFCav(const Config &conf, const int CavCnt,
               const double IonZ, const double IonEs, double &IonW, double &EkState,
               double &Fy_absState, double &beta, double &gamma)
{
    std::string CavType, CavTLMLine;
    int         cavi, cavilabel, multip;
    double      Rm, IonFy_i, Ek_i, avebeta, avegamma, fRF, CaviIonK, SampleIonK, EfieldScl;
    double      IonW_o, IonFy_o, accIonW;
    value_mat   M;

    CavType = conf.get<std::string>("cavtype");
    if (CavType == "0.041QWR") {
        cavi       = 1;
        cavilabel  = 41;
        multip     = 1;
        Rm         = 17e0;
        CavTLMLine = "/data/Multipole41/thinlenlon_41.txt";
    } else if (conf.get<std::string>("cavtype") == "0.085QWR") {
        cavi       = 2;
        cavilabel  = 85;
        multip     = 1;
        Rm         = 17e0;
        CavTLMLine = "/data/Multipole85/thinlenlon_85.txt";
    } else {
        std::cerr << "*** InitLong: undef. cavity type: " << CavType << "\n";
        exit(1);
    }

    IonFy_i = multip*Fy_absState + CavPhases[CavCnt-1];
    Ek_i    = EkState;
    IonW    = EkState + IonEs;

    avebeta    = beta;
    avegamma   = gamma;
    fRF        = conf.get<double>("f");
    CaviIonK   = 2e0*M_PI*fRF/(beta*C0*MtoMM);
    SampleIonK = 2e0*M_PI/(beta*C0/SampleFreq*MtoMM);
    EfieldScl  = conf.get<double>("scl_fac");         // Electric field scale factor.

    GetCavBoost(CavData[cavi-1], IonW, IonFy_i, CaviIonK, IonZ, IonEs,
                fRF, EfieldScl, IonW_o, IonFy_o);

    accIonW      = IonW_o - IonW;
    IonW         = IonW_o;
    EkState      = IonW - IonEs;
    IonW         = EkState + IonEs;
    gamma        = IonW/IonEs;
    beta         = sqrt(1e0-1e0/sqr(gamma));
    avebeta      = (avebeta+beta)/2e0;
    avegamma     = (avegamma+gamma)/2e0;
    Fy_absState += (IonFy_o-IonFy_i)/multip;

    GetCavMat(cavi, cavilabel, Rm, CavTLMLine, IonZ, IonEs, EfieldScl, IonFy_i, Ek_i, fRF, M);

//    TransVector[ii_state] = CaviMatrix.times(TransVector[ii_state]);
//    TransVector[ii_state](4, Fy_abs[ii_state]-tlmPara.Fy_abs_tab.get(lattcnt+1)[1]);
//    TransVector[ii_state](5, Ek[ii_state]-tlmPara.Ek_tab.get(lattcnt+1)[1]);

//    double aveX2i = ThetaMatrix[state].getElem(0, 0);
//    double aveY2i = ThetaMatrix[state].getElem(2, 2);
//    double IonFys = CavTLMLineTab[cavi-1].S[n-1]/180e0*M_PI;
//    double E0TL = accIonW/cos(IonFys)/tlmPara.IonChargeStates[state];
//    ThetaMatrix[state] = ThetaMatrix[state].conjugateTrans(CaviMatrix);
//    ThetaMatrix[state] = calRFcaviEmitGrowth(
//                ThetaMatrix[state], tlmPara.IonChargeStates[state], E0TL,
//                avebeta, avegamma, beta, gamma, caviLambda, IonFys, aveX2i,
//                TransVector[state].getElem(0), aveY2i, TransVector[state].getElem(2));
}


void InitLattice(Machine &sim, const double IonZ, const std::vector<double> BaryCenter)
{
    // Evaluate transport matrices for given beam initial conditions.
    typedef MatrixState state_t;

    std::stringstream               strm;
    int                             CavCnt;
    double                          IonW, IonEs, IonEk, s, L, beta, gamma;
    double                          SampleionK, R56, Brho, K, EkState, Fy_absState;
    Machine::p_elements_t::iterator it;
    MomentElementBase               *ElemPtr;
//    LinearElementBase<MatrixState>  *ElemPtr;

    Config                   D;
    std::auto_ptr<StateBase> state(sim.allocState(D));
    // Propagate through first element (beam initial conditions).
    sim.propagate(state.get(), 0, 1);

    IonEk = state->IonEk/MeVtoeV;
    IonEs = state->IonEs/MeVtoeV;
    IonW  = state->IonW/MeVtoeV;

    // Define initial conditions.
    Fy_absState = BaryCenter[state_t::PS_S];
    EkState     = IonEk + BaryCenter[state_t::PS_PS];

//    TransVector[state] = get(k).centerVector;
//    ThetaMatrix[state] = get(k).thetaMatrix;
//    Fy_absState      = get(k).centerVector.getElem(4);

    std::cout << "\n" << "InitLattice:" << "\n";
    std::cout << std::fixed << std::setprecision(5)
              << "  IonZ = " << IonZ
              << "  IonEs [Mev/u] = " << IonEs << ", IonEk [Mev/u] = " << IonEk
              << ", IonW [Mev/u] = " << IonW << "\n";

    s = 0e0;
    it = sim.p_elements.begin();
    // Skip over state.
    it++;
    CavCnt = 0;
    do {
        ElementVoid*  elem   = *it;
        const Config& conf   = elem->conf();
        Config        newconf(conf);
        std::string   t_name = elem->type_name(); // C string -> C++ string.

//        ElemPtr = dynamic_cast<LinearElementBase<MatrixState> *>(elem);
        ElemPtr = dynamic_cast<MomentElementBase *>(elem);
        assert(ElemPtr != NULL);

        if (t_name != "marker") {
            L = conf.get<double>("L");
            s += L;
        }

        gamma      = (EkState+IonEs)/IonEs;
        beta       = sqrt(1e0-1e0/sqr(gamma));
        SampleionK = 2e0*M_PI/(beta*SampleLambda);

        // Evaluate momentum compaction.
        R56 = -2e0*M_PI/(SampleLambda*IonEs*cube(beta*gamma))*L*MeVtoeV;

        if (t_name == "marker") {
        } else if (t_name == "drift") {
            Fy_absState += SampleionK*L;
            ElemPtr->transfer(state_t::PS_S, state_t::PS_PS) = R56;
        } else if (t_name == "sbend") {
            Fy_absState += SampleionK*L;
            ElemPtr->transfer(state_t::PS_S, state_t::PS_PS) = R56;
        } else if (t_name == "quadrupole") {
        } else if (t_name == "solenoid") {
            Fy_absState += SampleionK*L;
            Brho = beta*(EkState+IonEs)/(C0*IonZ);
            // Scale B field.
            K = conf.get<double>("B")/(2e0*Brho);
            size_t elem_index = ElemPtr->index;
            Config newconf(sim[elem_index]->conf());
            newconf.set<double>("K", K);
            sim.reconfigure(elem_index, newconf);
            // Re-initialize after re-allocation.
            elem = *it;
            ElemPtr = dynamic_cast<MomentElementBase *>(elem);

            ElemPtr->transfer(state_t::PS_S, state_t::PS_PS) = R56;
        } else if (t_name == "rfcavity") {
            CavCnt++;
            InitRFCav(conf, CavCnt, IonZ, IonEs, IonW, EkState, Fy_absState, beta, gamma);
//            PrtMat(ElemPtr->transfer);
        }
        it++;
    } while (it != sim.p_elements.end());
}


//void PropagateState(const Machine &sim, value_mat &S)
void PropagateState(const Machine &sim)
{
    Machine::p_elements_t::const_iterator it;
    MomentElementBase                     *ElemPtr;
    Config                                D;
    std::auto_ptr<StateBase>              state(sim.allocState(D));
    ElementVoid                           *elem;
    std::fstream                          outf;

    outf.open("test_jb.out", std::ofstream::out);

    std::cout << "\n" << "PropagateState:" << "\n";

    for (it = sim.p_elements.begin(); it != sim.p_elements.end(); ++it) {
        elem = *it;

        std::cout << "\n" << elem->type_name() << " " << elem->name << "\n";

//        ElemPtr = dynamic_cast<MomentElementBase *>(elem);
//        assert(ElemPtr != NULL);
//        PrtMat(ElemPtr->transfer);

//        sim.propagate(state.get(), elem->index, 1);
        elem->advance(*state);

        MatrixState* StatePtr = dynamic_cast<MatrixState*>(state.get());
//        std::cout << "\n" << *state;
        PrtMat(StatePtr->state);
    }

    outf.close();
}


void BaryCenterUnitFix(std::vector<double> &BaryCenter)
{
    // Change units from [mm, rad, mm, rad, rad, MeV/u] to [m, rad, m, rad, rad, eV/u].
    typedef MatrixState state_t;

    BaryCenter[state_t::PS_X]  /= MtoMM;
    BaryCenter[state_t::PS_Y]  /= MtoMM;
    BaryCenter[state_t::PS_PS] *= MeVtoeV;
}


void StateUnitFix(Machine &sim)
{
    /* Change units from:
     *   [mm, rad, mm, rad, rad, MeV/u]
     * to:
     *   [m, rad, m, rad, rad, eV/u].
     * Notes: State is only changed locally.
     *        Etot could be scaled wit p0. */

    int                      j, k;
    Config                   D;
    std::auto_ptr<StateBase> state(sim.allocState(D));

//    std::cout << "Machine configuration\n" << sim << "\n";

    // Propagate through first element (beam initial conditions).
    sim.propagate(state.get(), 0, 1);

    MatrixState* StatePtr = dynamic_cast<MatrixState*>(state.get());

    std::cout << "\n";
    PrtMat(StatePtr->state);

    for (j = 0; j < PS_Dim-2; j++)
        for (k = 0; k < PS_Dim-2; k++) {
            if (j % 2 == 0) StatePtr->state(j, k) /= MtoMM;
            if (k % 2 == 0) StatePtr->state(j, k) /= MtoMM;
        }

    for (j = 0; j < PS_Dim; j++)
        StatePtr->state(j, PS_Dim-1) *= MeVtoeV;

    std::cout << "\n";
    PrtMat(StatePtr->state);
}


//Machine::propagate_jb(StateBase* S, size_t start, size_t max) const
//{
//    const size_t nelem = p_elements.size();

//    for(size_t i=start; S->next_elem<nelem && i<max; i++)
//    {
//        ElementVoid* E = p_elements[S->next_elem];
//        S->next_elem++;
//        E->advance(*S);
//        if(p_trace)
//            (*p_trace) << "After "<< i << " " << *S;
//    }
//}


int main(int argc, char *argv[])
{
    int                 k;
    std::vector<double> ChgState;
    std::vector<double> BaryCenter[2];

    FILE *in = stdin;
    if(argc>1) {
        in = fopen(argv[1], "r");
        if (!in) {
            fprintf(stderr, "Failed to open %s\n", argv[1]);
            return 2;
        }
    }

    glps_debug = 0; // 0 or 1.

    std::auto_ptr<Config> conf;

    try {
        GLPSParser P;
        conf.reset(P.parse(in));
        fprintf(stderr, "Parsing succeeds\n");
    } catch(std::exception& e) {
        fprintf(stderr, "Parse error: %s\n", e.what());
        fclose(in);
        return 1;
    }

    // Register state and element types.
    registerLinear();
    registerMoment();

    Machine sim(*conf);
    if (false)
        sim.set_trace(&std::cout);
    else
        sim.set_trace(NULL);

    CavData[0].RdData(HomeDir+"/data/axisData_41.txt");
    CavData[1].RdData(HomeDir+"/data/axisData_85.txt");

    // Turn trace on/off.
    if (false) sim.set_trace(&std::cout); else sim.set_trace(NULL);

//    StateUnitFix(sim);

    // Charge states.
    ChgState      = conf->get<std::vector<double> >("IonChargeStates");
    BaryCenter[0] = conf->get<std::vector<double> >("BaryCenter1");
    BaryCenter[1] = conf->get<std::vector<double> >("BaryCenter2");

//    value_mat S1 = conf->get<value_mat>("S1");
    const std::vector<double>& S1vec = conf->get<std::vector<double> >("S1");

    value_mat S1(PS_Dim, PS_Dim);
    if (S1vec.size() > S1.data().size())
        throw std::invalid_argument("Initial state size too big");
    std::copy(S1vec.begin(), S1vec.end(), S1.data().begin());


    std::cout << "\n" << "Ion charge states:\n";
    for (k = 0; k < 2; k++)
        std::cout << std::fixed << std::setprecision(5) << std::setw(9) << ChgState[k];
    std::cout << "\n";
    std::cout << "\nBarycenter:\n";
    PrtVec(BaryCenter[0]);
    PrtVec(BaryCenter[1]);
    std::cout << "\nBeam envelope:\n";
    PrtMat(S1);

    InitLong(sim);

    InitLattice(sim, ChgState[0], BaryCenter[0]);
    InitLattice(sim, ChgState[1], BaryCenter[1]);

//    PrtLat(sim);

//    PropagateState(sim, S);
//    PropagateState(sim);

//     for n states
//       get Chg_n, S_n
//       Initialize Lattice(Chg_n, S_n)
//       PropagateState(sim, S);


//        std::cout<<"# Reduced lattice\n";
//        GLPSPrint(std::cout, *conf);
//        std::cout<<"\n";

//        std::cerr<<"Generic AST:\n";
//        std::cerr << *conf;
//        std::cerr<<"GLPS:\n";
//        GLPSPrint(std::cout, *conf);

    if (false) {
        try {
            Machine sim(*conf);
            sim.set_trace(&std::cout);

            std::cout << "# Machine configuration\n" << sim << "\n\n";

            Config D;
            std::auto_ptr<StateBase> state(sim.allocState(D));
            sim.propagate(state.get());

            std::cout << "\n# Final " << *state << "\n";
        } catch(std::exception& e) {
            std::cerr << "Simulation error: " << e.what() << "\n";
            fclose(in);
            return 1;
        }
    }

    fprintf(stderr, "Done\n");
    fclose(in);
    return 0;
}
