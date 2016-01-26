
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


// Global constansts.

// Speed of light [m/s].
# define c0           2.99792458e8
// Atomic mass unit [eV/c^2].
# define u            931.49432e6
// Vacuum permeability.
# define mu0          4e0*M_PI*1e-7
// Long. sampling frequency [Hz]; must be set to RF freq.
# define SampleFreq   80.5e6
# define SampleLambda c0/SampleFreq

// Global constansts and parameters.

// Charge stripper parameters.
const int    Stripper_n                 = 5;   // Number of charge states.
const double Stripper_IonZ              = 78e0/238e0,
             Stripper_IonMass           = 238e0,
             Stripper_IonProton         = 92e0,
             Stripper_IonChargeStates[] = {76e0/238e0, 77e0/238e0, 78e0/238e0,
                                           79e0/238e0, 80e0/238e0},
             // Energy dependance. Unit for last parmeter ???
             Stripper_E0Para[]          = {16.348e6, 1.00547, -0.10681},
             // Thickness [microns], thickness variation [%], reference energy [eV/u].
             StripperPara[]             = {3e0, 20e0, 16.623e6};


typedef boost::numeric::ublas::matrix<double> value_t;


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
};


// Global variables.

CavDataType         CavData[5];
LongTabType         LongTab;
std::vector<double> CavPhases;


void LongTabType::set(const double s, const double Ek, const double FyAbs,
                      const double Beta, const double Gamma)
{
    this->s.push_back(s);    this->Ek.push_back(Ek); this->FyAbs.push_back(FyAbs);
    this->Beta.push_back(Beta); this->Gamma.push_back(Gamma);
}


void LongTabType::show(std::ostream& strm, const int k) const
{
    strm << std::scientific << std::setprecision(10)
         << std::setw(18) << this->s[k]
         << std::setw(18) << this->Ek[k]*1e-6
         << std::setw(18) << this->FyAbs[k]
         << std::setw(18) << this->Beta[k]
         << std::setw(18) << this->Gamma[k] << "\n";
}


void CavDataType::RdData(const std::string FileName)
{
    std::string  line;
    int          k;
    double       s, Elong;
    std::fstream inf;

    inf.open(FileName.c_str(), std::ifstream::in);
    if (!inf.is_open()) {
        std::cerr << "Failed to open " << FileName << "\n";
        exit(1);
    }
    while (getline(inf, line)) {
        sscanf(line.c_str(), "%lf %lf", &s, &Elong);
        this->s.push_back(s), this->Elong.push_back(Elong);
        // Convert from [mm] to [m].
        this->s[this->s.size()-1] *= 1e-3;
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


void PrtMat(const value_t &M)
{
    int j, k;

    const int n = 6;

    std::cout << "\n";
    for (j = 0; j < n; j++) {
        for (k = 0; k < n; k++)
            std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << M(j, k);
        std::cout << "\n";
    }
}


void prt_lat(Machine &sim)
{
  std::stringstream strm;

  std::cout << "# Machine configuration\n" << sim << "\n\n";

  strm << "sim_type: " << sim.p_info.name << "\n#Elements: "
       << sim.p_elements.size() << "\n";
  for (Machine::p_elements_t::const_iterator it = sim.p_elements.begin(),
       end = sim.p_elements.end(); it != end; ++it) {
      (*it)->show(strm);
      std::cout << strm;
  }
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


double GetCavPhase(const int cavi, const double IonEk, const double IonFys,
                   const double FyAbs, const double multip)
{
    /* If the cavity is not at full power, the method gives synchrotron
     * phase slightly different from the nominal value.                 */

    double Fyc;

    switch (cavi) {
    case 1:
        Fyc = 4.394*pow(IonEk*1e-6, -0.4965) - 4.731;
        break;
    case 2:
        Fyc = 5.428*pow(IonEk*1e-6, -0.5008) + 1.6;
        break;
    case 3:
        Fyc = 22.35*pow(IonEk*1e-6, -0.5348) + 2.026;
        break;
    case 4:
        Fyc = 41.43*pow(IonEk*1e-6, -0.5775) + 2.59839;
        break;
    case 5:
        Fyc = 5.428*pow(IonEk*1e-6, -0.5008) + 1.6;
        break;
    default:
        std::cerr << "*** GetCavPhase: undef. cavity type" << "\n";
        exit(1);
    }

    return IonFys - Fyc - FyAbs*multip;
}

void GetCavBoost(const CavDataType &CavData, const double IonW0,
                 const double IonFy0, const double IonK0, const double IonZ,
                 const double IonEs, const double IonLamda,
                 const double EfieldScl, double &IonW, double &IonFy)
{
    int    n = CavData.s.size(),
           k;

    double dis = CavData.s[n-1] - CavData.s[0],
           dz  = dis/(n-1),
           IonK, IonFylast, IonGamma, IonBeta;

    IonFy = IonFy0;
    IonK  = IonK0;
    for (k = 0; k < n-1; k++) {
        IonFylast = IonFy;
        IonFy += IonK*dz;
        IonW  += IonZ*EfieldScl*(CavData.Elong[k]+CavData.Elong[k+1])/2e0
                 *cos((IonFylast+IonFy)/2e0)*dz;
        IonGamma = IonW/IonEs;
        IonBeta = sqrt(1e0-1e0/sqr(IonGamma));
        if ((IonW-IonEs) < 0e0) {
            IonW = IonEs;
            IonBeta = 0e0;
        }
        IonK = 2e0*M_PI/(IonBeta*IonLamda);
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


void PropagateLongRFCav(const Config &conf, const int n, const double IonZ, const double IonEs, double &IonW,
                        double &SampleIonK, double &IonBeta)
{
    std::string CavType;
    int         cavi;
    double      fRF, multip, caviIonK, IonFys, EfieldScl, caviFy, IonFy_i, IonFy_o, IonGamma;

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
    caviIonK  = 2e0*M_PI*fRF/(IonBeta*c0);
    IonFys    = conf.get<double>("phi")*M_PI/180e0; // Synchrotron phase [rad].
    EfieldScl = conf.get<double>("scl_fac");       // Electric field scale factor.

    caviFy = GetCavPhase(cavi, IonW-IonEs, IonFys, LongTab.FyAbs[n-2], multip);

    IonFy_i = multip*LongTab.FyAbs[n-2] + caviFy;
    CavPhases.push_back(caviFy);
    if (false)
        std::cout << std::scientific << std::setprecision(10)
                  << "CavPhase: " << std::setw(3) << CavPhases.size()
                  << std::setw(18) << CavPhases[CavPhases.size()-1] << "\n";
    // Evaluate change of reference particle kinetic energy, absolute phase, beta, and gamma.
    GetCavBoost(CavData[cavi-1], IonW, IonFy_i, caviIonK, IonZ,
                IonEs, c0/fRF, EfieldScl, IonW, IonFy_o);
    IonGamma   = IonW/IonEs;
    IonBeta    = sqrt(1e0-1e0/sqr(IonGamma));
    SampleIonK = 2e0*M_PI/(IonBeta*SampleLambda);

    LongTab.set(LongTab.s[n-2]+conf.get<double>("L"), IonW-IonEs,
            LongTab.FyAbs[n-2]+(IonFy_o-IonFy_i)/multip, IonBeta, IonGamma);
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
    IonEk = (LongTab.Ek[n-2]-StripperPara[2])*Stripper_E0Para[1] + Stripper_E0Para[0];
    IonW        = IonEk + IonEs;
    IonGamma    = IonW/IonEs;
    IonBeta     = sqrt(1e0-1e0/sqr(IonGamma));
    SampleIonK  = 2e0*M_PI/(IonBeta*SampleLambda);

    LongTab.set(LongTab.s[n-2], IonEk, LongTab.FyAbs[n-2], IonBeta, IonGamma);

    //            chargeAmount = fribstripper.chargeAmount_Baron;
}


void InitLong(Machine &sim)
{
    /* Longitudinal initialization for reference particle.
     * Evaluate beam energy and cavity loaded phase along the lattice. */
    int                             n;
    double                          IonGamma, IonBeta, SampleIonK;
    double                          IonW, IonZ, IonEs;
    Machine::p_elements_t::iterator it;

    Config                   D;
    std::auto_ptr<StateBase> state(sim.allocState(D));
    // Propagate through first element.
    sim.propagate(state.get(), 0, 1);

    IonZ  = state->IonZ;
    IonEs = state->IonEs;
    IonW  = state->IonW;

    IonGamma   = IonW/IonEs;
    IonBeta    = sqrt(1e0-1e0/sqr(IonGamma));
    SampleIonK = 2e0*M_PI/(IonBeta*SampleLambda);

    std::cout << "\n" << "InitLong:" << "\n\n";
    std::cout << std::scientific << std::setprecision(5)
              << "IonEs [Mev/u] = " << IonEs*1e-6 << ", IonEk [Mev/u] = " << state->IonEk*1e-6
              << ", IonW [Mev/u] = " << IonW*1e-6 << "\n";

    n = 1;
    LongTab.set(0e0, state->IonEk, 0e0, IonBeta, IonGamma);

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
                    LongTab.FyAbs[n-2]+SampleIonK*conf.get<double>("L"),
                    LongTab.Beta[n-2], LongTab.Gamma[n-2]);
        } else if (t_name == "sbend") {
            n++;
            LongTab.set(LongTab.s[n-2]+conf.get<double>("L"), LongTab.Ek[n-2],
                    LongTab.FyAbs[n-2]+SampleIonK*conf.get<double>("L"),
                    LongTab.Beta[n-2], LongTab.Gamma[n-2]);
        } else if (t_name == "quadrupole") {
            n++;
            LongTab.set(LongTab.s[n-2]+conf.get<double>("L"), LongTab.Ek[n-2],
                    LongTab.FyAbs[n-2]+SampleIonK*conf.get<double>("L"),
                    LongTab.Beta[n-2], LongTab.Gamma[n-2]);
        } else if (t_name == "solenoid") {
            n++;
            LongTab.set(LongTab.s[n-2]+conf.get<double>("L"), LongTab.Ek[n-2],
                    LongTab.FyAbs[n-2]+SampleIonK*conf.get<double>("L"),
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
            std::cout << std::setw(10) << std::left << t_name << std::internal;
            LongTab.show(std::cout, n-1);
        }
        it++;
    } while (it != sim.p_elements.end());
}


void GetTransicFac(const int cavilabel, double beta, const int gaplabel, const double E_fac_ad,
                 double &Ecen, double &T, double &Tp, double &S, double &Sp, double &V0)
{
    // Evaluate Electric field center, transit factors [T, T', S, S'] and cavity field.

    switch (cavilabel) {
    case 41:
        if (beta < 0.025 || beta > 0.08) {
            std::cerr << "*** GetTransicFac: beta out of Range" << "\n";
            exit(1);
        }
        switch (gaplabel) {
        case 0:
            // One gap evaluation.
            Ecen = 120.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = 4.957e6*pow(beta, 5.0) - 1.569e6*pow(beta, 4.0) + 1.991e5*pow(beta, 3.0)
                      - 1.269e4*pow(beta, 2.0) + 399.9*beta - 4.109;
            Sp   = -2.926e8*pow(beta, 5.0) + 8.379e7*pow(beta, 4.0) + 9.284e6*pow(beta, 3.0)
                      + 4.841e5*pow(beta, 2.0) - 1.073e4*beta + 61.98;
            V0   = 0.98477*E_fac_ad;
            break;
        case 1:
            // Two gap calculation, first gap.
            Ecen = 0.0006384*pow(beta, -1.884) + 86.69;
            T    = -1.377e6*pow(beta, 5.0) + 4.316e5*pow(beta, 4.0) - 5.476e4*pow(beta, 3.0)
                      + 3570*pow(beta, 2.0) - 123.2*beta + 0.9232;
            Tp   = 2.277e7*pow(beta, 5.0) - 6.631e6 *pow(beta, 4.0) + 7.528e5*pow(beta, 3.0)
                      - 4.062e4*pow(beta, 2.0) + 924.7*beta + 1.699;
            S    = 0.0;
            Sp   =-1.335e6*pow(beta, 5.0) + 3.385e5*pow(beta, 4.0) - 2.98e4*pow(beta, 3.0)
                      + 806.6*pow(beta, 2.0) + 25.59*beta - 1.571;
            V0   = 0.492385*E_fac_ad;
            break;
        case 2:
            // Two gap calculation, second gap.
            Ecen = -0.0006384*pow(beta, -1.884) + 33.31;
            T    = 1.377e6*pow(beta, 5.0) - 4.316e5*pow(beta, 4.0) + 5.476e4*pow(beta, 3.0)
                      - 3570*pow(beta, 2.0) + 123.2*beta - 0.9232;
            Tp   = -2.277e7*pow(beta, 5.0) + 6.631e6*pow(beta, 4.0) - 7.528e5*pow(beta, 3.0)
                      + 4.062e4*pow(beta, 2.0) - 924.7*beta - 1.699;
            S    = 0.0;
            Sp   =-1.335e6*pow(beta, 5.0) +  3.385e5*pow(beta, 4.0) - 2.98e4*pow(beta, 3.0)
                      + 806.6*pow(beta, 2.0) + 25.59*beta - 1.571;
            V0   = 0.492385*E_fac_ad;
            break;
        default:
            std::cerr << "*** GetTransicFac: undef. number of gaps " << gaplabel << "\n";
            exit(1);
        }
        break;
    case 85:
        if (beta < 0.05 || beta > 0.25) {
            std::cerr << "*** GetTransicFac: beta out of range" << "\n";
            exit(1);
        }
        switch (gaplabel) {
          case 0:
            Ecen = 150.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = 2.326e6*pow(beta, 7.0) - 2.781e6*pow(beta, 6.0) + 1.407e6*pow(beta, 5.0)
                      - 3.914e5*pow(beta, 4.0) + 6.477e4*pow(beta, 3.0) - 6385*pow(beta, 2.0)
                      + 343.9*beta - 6.811;
            Sp   =-2.755e8*pow(beta,7.0) + 3.109e8*pow(beta,6.0) - 1.462e8*pow(beta, 5.0)
                      + 3.691e7*pow(beta, 4.0) - 5.344e6*pow(beta, 3.0) + 4.315e5*pow(beta, 2.0)
                      - 1.631e4*beta + 162.7;
            V0   = 1.967715*E_fac_ad;
            break;
        case 1:
            Ecen = 0.0002838*pow(beta, -2.13) + 76.5;
            T    = 0.0009467*pow(beta, -1.855) - 1.002;
            Tp   = -1.928e4*pow(beta, 5.0) + 2.195e4*pow(beta, 4.0) - 1.017e4*pow(beta, 3.0)
                      + 2468*pow(beta, 2.0) - 334*beta + 24.44;
            S    = 0.0;
            Sp   =-0.0009751*pow(beta, -1.898) + 0.001568;
            V0   = 0.9838574*E_fac_ad;
            break;
        case 2:
            Ecen = -0.0002838*pow(beta, -2.13) + 73.5;
            T    = -0.0009467*pow(beta, -1.855) + 1.002;
            Tp   = 1.928e4*pow(beta, 5.0) - 2.195e4*pow(beta, 4.0) + 1.017e4*pow(beta, 3.0)
                      - 2468*pow(beta, 2.0) + 334*beta - 24.44;
            S    = 0.0;
            Sp   =-0.0009751*pow(beta, -1.898) + 0.001568;
            V0   = 0.9838574*E_fac_ad;
            break;
        default:
            std::cerr << "*** GetTransicFac: undef. number of gaps " << gaplabel << "\n";
            exit(1);
        }
        break;
    case 29:
        if (beta < 0.15 || beta > 0.4) {
            std::cerr << "*** GetTransicFac: beta out of range" << "\n";
            exit(1);
        }
        switch (gaplabel) {
        case 0:
            Ecen = 150.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = 76.54*pow(beta, 5.0) - 405.6*pow(beta, 4.0) + 486*pow(beta, 3.0)
                      - 248*pow(beta, 2.0) + 58.08*beta - 4.285;
            Sp   =-2.025e6*pow(beta,7.0) + 4.751e6*pow(beta,6.0) - 4.791e6*pow(beta, 5.0)
                      + 2.695e6*pow(beta, 4.0) - 9.127e5*pow(beta, 3.0) + 1.854e5*pow(beta, 2.0)
                      - 2.043e4*beta + 888;
            V0   = 2.485036*E_fac_ad;
            break;
        case 1:
            Ecen = 0.01163*pow(beta, -2.001) + 91.77;
            T    = 0.02166*pow(beta, -1.618) - 1.022;
            Tp   = 1.389e4*pow(beta, 5.0) - 2.147e4*pow(beta, 4.0) + 1.313e4*pow(beta, 3.0)
                      - 3917*pow(beta, 2.0) + 534.7*beta - 11.25;
            S    = 0.0;
            Sp   =-454.4*pow(beta, 5.0) + 645.1*pow(beta, 4.0) - 343.9*pow(beta, 3.0)
                      + 78.77*pow(beta, 2.0) - 4.409*beta - 0.8283;
            V0   = 1.242518*E_fac_ad;
        case 2:
            Ecen = -0.01163*pow(beta, -2.001) + 58.23;
            T    = -0.02166*pow(beta, -1.618) + 1.022;
            Tp   = -1.389e4*pow(beta, 5.0) + 2.147e4*pow(beta, 4.0) - 1.313e4*pow(beta, 3.0)
                      + 3917*pow(beta, 2.0) - 534.7*beta + 11.25;
            S    = 0.0;
            Sp   =-454.4*pow(beta, 5.0) + 645.1*pow(beta, 4.0) - 343.9*pow(beta, 3.0)
                      + 78.77*pow(beta, 2.0) - 4.409*beta - 0.8283;
            V0   = 1.242518*E_fac_ad;
            break;
        default:
            std::cerr << "*** GetTransicFac: undef. number of gaps " << gaplabel << "\n";
            exit(1);
        }
        break;
    case 53:
        if (beta < 0.3 || beta > 0.6) {
            std::cerr << "*** GetTransicFac: beta out of range" << "\n";
            exit(1);
        }
        switch (gaplabel) {
        case 0:
            Ecen = 250.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = -52.93*pow(beta, 5.0) + 84.12*pow(beta, 4.0) - 17.73*pow(beta, 3.0)
                      - 38.49*pow(beta, 2.0) + 26.64*beta - 4.222;
            Sp   =-4.167e4*pow(beta, 5.0) + 1.075e5*pow(beta, 4.0) - 1.111e5*pow(beta, 3.0)
                      + 5.702e4*pow(beta, 2.0) - 1.413e4*beta - 1261;
            V0   = 4.25756986*E_fac_ad;
            break;
        case 1:
            Ecen = 0.01219*pow(beta, -2.348) + 137.8;
            T    = 0.04856*pow(beta, -1.68) - 1.018;
            Tp   = 1641*pow(beta, 5.0) - 4109*pow(beta, 4.0) + 4081*pow(beta, 3.0)
                      - 1973*pow(beta, 2.0) + 422.8*beta - 3.612;
            S    = 0.0;
            Sp   =-0.03969*pow(beta, -1.775) + 0.009034;
            V0   = 2.12878493*E_fac_ad;
            break;
        case 2:
            Ecen = -0.01219*pow(beta, -2.348) + 112.2;
            T    = -0.04856*pow(beta, -1.68) + 1.018;
            Tp   = -1641*pow(beta, 5.0) + 4109*pow(beta, 4.0) - 4081*pow(beta, 3.0)
                      + 1973*pow(beta, 2.0) - 422.8*beta + 3.612;
            S    = 0.0;
            Sp   = -0.03969*pow(beta, -1.775) + 0.009034;
            V0   = 2.12878493*E_fac_ad;
            break;
        default:
            std::cerr << "*** GetTransicFac: undef. number of gaps " << gaplabel << "\n";
            exit(1);
        }
        break;
    default:
        std::cerr << "*** GetTransicFac: undef. cavity type" << "\n";
        exit(1);
    }

    // Convert from [mm] to [m].
    Ecen *= 1e-3;
}


void CavityMatrix(double dis, double E_fac_ad, double TTF_tab[], double beta_tab[], double gamma_tab[],
                  double lamda, double ionZ, double ionFy0, double ionFyc,
                  std::string *thinlenLine, double Rm)
{
    /* RF cavity model, transverse only defocusing
     * 2-gap matrix model                                            */

    int    seg;
    double ki_s, kc_s, kf_s;
    double Ecen1, T_1, S_1, V0_1, k_1, L1;
    double Ecen2, T_2, S_2, V0_2, k_2, L2;
    double L3;
    double beta, gamma, ionFy, k, V0, T, S, kfdx, kfdy, dpy, acc;


//    PhaseMatrix matrix = PhaseMatrix.identity();

//    ki_s = 2*M_PI/(beta_tab[0]*lamda*1e3); // rad/mm
//    kc_s = 2*M_PI/(beta_tab[1]*lamda*1e3); // rad/mm
//    kf_s = 2*M_PI/(beta_tab[2]*lamda*1e3); // rad/mm

//    // Longitudinal model: Drift-Kick-Drift, dis: total lenghth centered at 0,
//    // Ecen1&Ecen2: Electric Center position where acc kick applies, Ecen1<0
//    // TTFtab: 2*6 vector, Ecen, T Tp S Sp, V0;

//    Ecen1 = TTF_tab[0];
//    T_1   = TTF_tab[1];
//    S_1   = TTF_tab[3];
//    V0_1  = TTF_tab[5];
//    k_1   = 0.5*(ki_s+kc_s);
//    L1    = dis + Ecen1; //try change dis/2 to dis 14/12/12
//    PhaseMatrix Mlon_L1 = PhaseMatrix.identity();
//    PhaseMatrix Mlon_K1 = PhaseMatrix.identity();
//    // Pay attention, original is -
//    Mlon_L1.setElem(4, 5,
//                    -2*M_PI/lamda/1e3*1e0/cube(beta_tab[0]*gamma_tab[0])/FRIBPara.ionEs*L1));
//     Pay attention, original is -k1-k2
//    Mlon_K1.setElem(5, 4, -ionZ*V0_1*T_1*Math.sin(ionFy0+k_1*L1)-ionZ*V0_1*S_1*Math.cos(ionFy0+k_1*L1));

//    Ecen2 = TTF_tab[6];
//    T_2   = TTF_tab[7];
//    S_2   = TTF_tab[9];
//    V0_2  = TTF_tab[11];
//    k_2   = 0.5*(kc_s+kf_s);
//    L2    = Ecen2-Ecen1;
//    PhaseMatrix Mlon_L2 = PhaseMatrix.identity();
//    PhaseMatrix Mlon_K2 = PhaseMatrix.identity();
//    Mlon_L2.setElem(4, 5,
//                    -2*M_PI/lamda/1e3*1e0/cube(beta_tab[1]*gamma_tab[1])/FRIBPara.ionEs*L2)); //Problem is Here!!
//    Mlon_K2.setElem(5, 4, -ionZ*V0_2*T_2*Math.sin(ionFyc+k_2*Ecen2)-ionZ*V0_2*S_2*Math.cos(ionFyc+k_2*Ecen2));


//    L3 = dis-Ecen2; //try change dis/2 to dis 14/12/12
//    PhaseMatrix Mlon_L3 = PhaseMatrix.identity();
//    Mlon_L3.setElem(4, 5,
//                    -2*M_PI/lamda/1e3*1e0/cube(beta_tab[2]*gamma_tab[2])/FRIBPara.ionEs*L3));

//    PhaseMatrix Mlon = PhaseMatrix.identity();
//    Mlon = Mlon_L3.times(Mlon_K2.times(Mlon_L2.times(Mlon_K1.times(Mlon_L1))));

//    // Transverse model
//    // Drift-FD-Drift-LongiKick-Drift-FD-Drift-0-Drift-FD-Drift-LongiKick-Drift-FD-Drift

//    PhaseMatrix Mtrans = PhaseMatrix.identity();
//    PhaseMatrix Mprob = PhaseMatrix.identity();
//    beta  = beta_tab[0];
//    gamma = gamma_tab[0];
//    seg   = 0;
//    ionFy = ionFy0;
//    k     = ki_s;
//    V0    = 0.0;
//    T     = 0.0;
//    S     = 0.0;
//    kfdx  = 0.0;
//    kfdy  = 0.0;
//    dpy   = 0.0;

//    for(TlmNode fribnode:thinlenLine) {
//        if (fribnode.label! = null) {
//            if (fribnode.label == "drift") {
//                ionFy  = ionFy + k*fribnode.length; // lenghth already in mm
//                Mprob  = PhaseMatrix.identity();
//                Mprob.setElem(0, 1, fribnode.length);
//                Mprob.setElem(2, 3, fribnode.length);
//                Mtrans = Mprob.times(Mtrans);
//            } else if (fribnode.label == "EFocus1") {
//                V0     = fribnode.attribute[1]*E_fac_ad;
//                T      = fribnode.attribute[2];
//                S      = fribnode.attribute[3];
//                kfdx   = ionZ*V0/(beta*beta)/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy)-S*Math.sin(ionFy))/Rm;
//                kfdy   = ionZ*V0/(beta*beta)/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy)-S*Math.sin(ionFy))/Rm;
//                Mprob = PhaseMatrix.identity();
//                Mprob.setElem(1, 0, kfdx);
//                Mprob.setElem(3, 2, kfdy);
//                Mtrans = Mprob.times(Mtrans);
//            } else if (fribnode.label == "EFocus2") {
//                V0     = fribnode.attribute[1]*E_fac_ad;
//                T      = fribnode.attribute[2];
//                S      = fribnode.attribute[3];
//                kfdx   = ionZ*V0/(beta*beta)/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy)-S*Math.sin(ionFy))/Rm;
//                kfdy   = ionZ*V0/(beta*beta)/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy)-S*Math.sin(ionFy))/Rm;
//                Mprob  = PhaseMatrix.identity();
//                Mprob.setElem(1, 0, kfdx);
//                Mprob.setElem(3, 2, kfdy);
//                Mtrans = Mprob.times(Mtrans);
//            } else if (fribnode.label == "EDipole") {
//                if (multipoleLevel<1) break;
//                V0     = fribnode.attribute[1]*E_fac_ad;
//                T      = fribnode.attribute[2];
//                S      = fribnode.attribute[3];
//                dpy    = ionZ*V0/(beta*beta)/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy)-S*Math.sin(ionFy));
//                Mprob  = PhaseMatrix.identity();
//                Mprob.setElem(3, 6, dpy);
//                Mtrans = Mprob.times(Mtrans);
//            } else if (fribnode.label == "EQuad") {
//                if (multipoleLevel<2) break;
//                V0     = fribnode.attribute[1]*E_fac_ad;
//                T      = fribnode.attribute[2];
//                S      = fribnode.attribute[3];
//                kfdx   = ionZ*V0/(beta*beta)/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy)-S*Math.sin(ionFy))/Rm;
//                kfdy   = -ionZ*V0/(beta*beta)/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy)-S*Math.sin(ionFy))/Rm;
//                Mprob  = PhaseMatrix.identity();
//                Mprob.setElem(1, 0, kfdx);
//                Mprob.setElem(3, 2, kfdy);
//                Mtrans = Mprob.times(Mtrans);
//            } else if (fribnode.label == "HMono") {
//                if (multipoleLevel<2) break;
//                V0     = fribnode.attribute[1]*E_fac_ad;
//                T      = fribnode.attribute[2];
//                S      = fribnode.attribute[3];
//                kfdx   = -FRIBPara.MU0*FRIBPara.C*ionZ*V0/beta/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy+M_PI/2)-S*Math.sin(ionFy+M_PI/2))/Rm;
//                kfdy   = -FRIBPara.MU0*FRIBPara.C*ionZ*V0/beta/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy+M_PI/2)-S*Math.sin(ionFy+M_PI/2))/Rm;
//                Mprob  = PhaseMatrix.identity();
//                Mprob.setElem(1, 0, kfdx);
//                Mprob.setElem(3, 2, kfdy);
//                Mtrans = Mprob.times(Mtrans);
//            } else if (fribnode.label == "HDipole") {
//                if (multipoleLevel<1) break;
//                V0     = fribnode.attribute[1]*E_fac_ad;
//                T      = fribnode.attribute[2];
//                S      = fribnode.attribute[3];
//                dpy    = -FRIBPara.MU0*FRIBPara.C*ionZ*V0/beta/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy+M_PI/2)-S*Math.sin(ionFy+M_PI/2));
//                Mprob  = PhaseMatrix.identity();
//                Mprob.setElem(3, 6, dpy);
//                Mtrans = Mprob.times(Mtrans);
//            } else if (fribnode.label == "HQuad") {
//                if (multipoleLevel < 2) break;
//                if (fribnode.position < 0) {
//                    // First gap.
//                    beta  = (beta_tab[0]+beta_tab[1])/2.0;
//                    gamma = (gamma_tab[0]+gamma_tab[1])/2.0;
//                } else {
//                    beta = (beta_tab[1]+beta_tab[2])/2.0;
//                    gamma = (gamma_tab[1]+gamma_tab[2])/2.0;
//                }
//                V0     = fribnode.attribute[1]*E_fac_ad;
//                T      = fribnode.attribute[2];
//                S      = fribnode.attribute[3];
//                kfdx   = -FRIBPara.MU0*FRIBPara.C*ionZ*V0/beta/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy+M_PI/2)-S*Math.sin(ionFy+M_PI/2))/Rm;
//                kfdy   = FRIBPara.MU0*FRIBPara.C*ionZ*V0/beta/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy+M_PI/2)-S*Math.sin(ionFy+M_PI/2))/Rm;
//                Mprob  = PhaseMatrix.identity();
//                Mprob.setElem(1, 0, kfdx);
//                Mprob.setElem(3, 2, kfdy);
//                Mtrans = Mprob.times(Mtrans);
//            } else if (fribnode.label == "AccGap") {
//                //ionFy = ionFy + ionZ*V0_1*k*(TTF_tab[2]*Math.sin(ionFy)
//                //        + TTF_tab[4]*Math.cos(ionFy))/2/((gamma-1)*FRIBPara.ionEs); //TTF_tab[2]~Tp
//                seg    = seg + 1;
//                beta   = beta_tab[seg];
//                gamma  = gamma_tab[seg];
//                k      = 2*M_PI/(beta*lamda*1e3); // rad/mm
//                acc    = fribnode.attribute[1];
//                Mprob  = PhaseMatrix.identity();
//                Mprob.setElem(1, 1, acc);
//                Mprob.setElem(3, 3, acc);
//                Mtrans = Mprob.times(Mtrans);
//            else {
//                std::cerr << "*** CavityMatrix: undef. multipole type " << fribnode.label << "\n";
//                exit(1);
//            }
//        }
//    }

//    matrix = Mtrans;
//    matrix.setElem(4, 4, Mlon.getElem(4, 4));
//    matrix.setElem(4, 5, Mlon.getElem(4, 5));
//    matrix.setElem(5, 4, Mlon.getElem(5, 4));
//    matrix.setElem(5, 5, Mlon.getElem(5, 5));

//    return matrix;
}


void CavityTLMMatrix(const double Rm, const double ionZ, const double IonEs, const double E_fac_ad,
                     const double ionFyi_s, const  double ionEk_s,
                     const double ionLamda, int cavi)
{
    int    k;
    double ionWi_s, gammai_s, betai_s, ionKi_s;

    ionWi_s  = ionEk_s + IonEs;
    gammai_s = ionWi_s/IonEs;
    betai_s  = sqrt(1e0-1e0/sqr(gammai_s));
    ionKi_s  = 2e0*M_PI/(betai_s*ionLamda);

//	    axisData_41_ad(:,cExRe:9)=FRIBPara.QWR1.E_fac_ad*axisData_41_ad(:,cExRe:9);
//	    axisData_41_1_ad=axisData_41_ad(1:FRIBPara.planeSize_z_41,:);
//	    axisData_41_2_ad=axisData_41_ad(FRIBPara.planeSize_z_41:2*FRIBPara.planeSize_z_41-1,:);
//

//    GetTransicFac(cavilabel, betai_s, 1 , E_fac_ad);

//    util.ArrayList<double[]> axisData_1_ad = new util.ArrayList<double[]>();
//    for (k = 0; k < round((axisData.size()-1)/2.0); k++) {
//        double[] entry = {axisData.get(k)[0],axisData.get(k)[1]*E_fac_ad};
//        axisData_1_ad.add(entry);
//    }
//    output = calTransfac(axisData_1_ad,ionKi_s);

//    double Ecen_1 = output[0]; //output = {Ecen,T,Tp,S,Sp,V0};
//    double T_1 = output[1];
//    double Tp_1 = output[2];
//    double S_1 = output[3];
//    double Sp_1 = output[4];
//    double V0_1 = output[5];
//    double dis = (axisData.get(axisData.size()-1)[0]-axisData.get(0)[0])/2; // should divided by 2!!!
//    output = calGapModel(dis,ionWi_s,ionFyi_s,ionKi_s,ionZ,ionLamda,Ecen_1,T_1,S_1,Tp_1,Sp_1,V0_1);
//    double ionWc_s = output[0];
//    double ionFyc_s = output[1];
//    double gammac_s = ionWc_s/FRIBPara.ionEs;
//    double betac_s = sqrt(1.0-1.0/(gammac_s*gammac_s));
//    double ionKc_s = 2*M_PI/(betac_s*ionLamda*1e3); //rad/mm

//    output = GetTransicFac(cavilabel,betac_s,2,E_fac_ad);
//    util.ArrayList<double[]> axisData_2_ad = new util.ArrayList<double[]>();
//    for (int k = 0; k < round((axisData.size()-1)/2.0); k++) {
//        double[] entry = {axisData.get(k)[0],axisData.get(k)[1]*E_fac_ad};
//        axisData_2_ad.add(entry);
//    }
//    output = calTransfac(axisData_2_ad,ionKc_s);
//    double Ecen_2 = output[0]; //output = {Ecen,T,Tp,S,Sp,V0};
//    double T_2 = output[1];
//    double Tp_2 = output[2];
//    double S_2 = output[3];
//    double Sp_2 = output[4];
//    double V0_2 = output[5];
//    output = calGapModel(dis,ionWc_s,ionFyc_s,ionKc_s,ionZ,ionLamda,Ecen_2,T_2,S_2,Tp_2,Sp_2,V0_2);
//    double ionWf_s = output[0];
//    double gammaf_s = ionWf_s/FRIBPara.ionEs;
//    double betaf_s = sqrt(1.0-1.0/(gammaf_s*gammaf_s));
//    double ionKf_s = 2*M_PI/(betaf_s*ionLamda*1e3); //rad/mm

//    Ecen_1 = Ecen_1-dis; // Switch Ecen_1 axis centered at cavity center// dis/2 to dis 14/12/12

//    double[] TTF_tab = {Ecen_1,T_1,Tp_1,S_1,Sp_1,V0_1,Ecen_2,T_2,Tp_2,S_2,Sp_2,V0_2};
//    double[] beta_tab = {betai_s,betac_s,betaf_s};
//    double[] gamma_tab = {gammai_s,gammac_s,gammaf_s};
//    double ionK_1 = (ionKi_s+ionKc_s)/2;
//    double ionK_2 = (ionKc_s+ionKf_s)/2;		    // Bug found 2015-01-07

//    PhaseMatrix matrix = PhaseMatrix.identity();
//    util.ArrayList<TlmNode> thinlenLine  =  new util.ArrayList<TlmNode>();
//    thinlenLine = calCavity_thinlenLine(beta_tab,gamma_tab,cavi,ionK_1,ionK_2);
//    matrix = calCavity_Matrix(dis,E_fac_ad,TTF_tab,beta_tab,gamma_tab,ionLamda,ionZ,ionFyi_s,ionFyc_s,thinlenLine,Rm);

//    return matrix;
}


void InitRFCav(const Config &conf, const int &CavCnt, const double IonZ, const double IonEs, double &IonW,
               const double IonChargeState, double &EkState, const double Fy_absState,
               double &beta, double &gamma)
{
    std::string CavType;
    int         cavi, multip;
    double      Rm, IonFy_i, Ek_i, avebeta, avegamma, fRF, CaviIonK, SampleIonK, EfieldScl;
    double      IonW_o, IonFy_o, accIonW;

    CavType = conf.get<std::string>("cavtype");
    if (CavType == "0.041QWR") {
        cavi   = 1;
        multip = 1;
        Rm     = 17e0;
    } else if (conf.get<std::string>("cavtype") == "0.085QWR") {
        cavi   = 2;
        multip = 1;
        Rm     = 17e0;
    } else {
        std::cerr << "*** InitLong: undef. cavity type: "
                  << CavType << "\n";
        exit(1);
    }

    IonFy_i = multip*Fy_absState + CavPhases[CavCnt-1];
    Ek_i    = EkState;
    IonW    = EkState + IonEs;

    avebeta    = beta;
    avegamma   = gamma;
    fRF        = conf.get<double>("f");
    CaviIonK   = 2e0*M_PI*fRF/(beta*c0);
    SampleIonK = 2e0*M_PI/(beta*c0/SampleFreq);
    EfieldScl  = conf.get<double>("scl_fac");   // Electric field scale factor.

    GetCavBoost(CavData[cavi-1], IonW, IonFy_i, CaviIonK, IonChargeState, IonEs,
                c0/fRF, EfieldScl, IonW_o, IonFy_o);

    accIonW = IonW_o - IonW;
    IonW    = IonW_o;
    EkState = IonW - IonEs;
    IonW    = EkState + IonEs;
    gamma   = IonW/IonEs;
    beta    = sqrt(1e0-1e0/sqr(gamma));
    avebeta = (avebeta+beta)/2e0;
    avegamma = (avegamma+gamma)/2e0;

//    CaviMatrix = calCavityTLM_Matrix(
//                Rm, IonChargeState, EfieldScl, IonFy_i, Ek_i, caviLamda, cavi);

//    TransVector[ii_state] = CaviMatrix.times(TransVector[ii_state]);
//    TransVector[ii_state].setElem(4, Fy_abs[ii_state]-tlmPara.Fy_abs_tab.get(lattcnt+1)[1]);
//    TransVector[ii_state].setElem(5, Ek[ii_state]-tlmPara.Ek_tab.get(lattcnt+1)[1]);

//    double aveX2i = ThetaMatrix[state].getElem(0, 0);
//    double aveY2i = ThetaMatrix[state].getElem(2, 2);
//    double IonFys = fribnode.attribute[3]/180e0*M_PI;
//    double E0TL = accIonW/cos(IonFys)/tlmPara.IonChargeStates[state];
//    ThetaMatrix[state] = ThetaMatrix[state].conjugateTrans(CaviMatrix);
//    ThetaMatrix[state] = calRFcaviEmitGrowth(
//                ThetaMatrix[state], tlmPara.IonChargeStates[state], E0TL,
//                avebeta, avegamma, beta, gamma, caviLamda, IonFys, aveX2i,
//                TransVector[state].getElem(0), aveY2i, TransVector[state].getElem(2));
}


void InitLattice(Machine &sim, const double P)
{
    // Evaluate transport matrices for given beam initial conditions.
    typedef MatrixState state_t;

    std::stringstream               strm;
    int                             CavCnt;
    double                          IonW, IonEs, IonEk, IonZ, IonLambda, s, L, beta, gamma;
    double                          R56, Brho, K, Ek_ini, EkState;
    Machine::p_elements_t::iterator it;
    MomentElementBase               *ElemPtr;
//    LinearElementBase<MatrixState> *ElemPtr;

    const double QWR1f = 80.5e6, QWR1Lambda = c0/QWR1f;

    Config                   D;
    std::auto_ptr<StateBase> state(sim.allocState(D));
    // Propagate through first element (beam initial conditions).
    sim.propagate(state.get(), 0, 1);

    IonZ  = state->IonZ;
    IonEk = state->IonEk;
    IonEs = state->IonEs;
    IonW  = state->IonW;

    IonLambda = QWR1Lambda;

    // Define initial conditions.
    EkState          = IonEk + P;
//    TransVector[state] = get(k).centerVector;
//    ThetaMatrix[state] = get(k).thetaMatrix;
//    Fy_absState      = get(k).centerVector.getElem(4);

    std::cout << "\n" << "InitLattice:" << "\n\n";
//    std::cout << *state << "\n";
    std::cout << std::scientific << std::setprecision(5)
              << "IonEs [Mev/u] = " << IonEs*1e-6 << ", IonEk [Mev/u] = " << state->IonEk*1e-6
              << ", IonW [Mev/u] = " << IonW*1e-6 << ", IonLambda [m^-1] = " << IonLambda << "\n";

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

        std::cout << t_name << " " << elem->name << "\n";

        if (t_name != "marker") {
            L = conf.get<double>("L");
            s += L;
        }

        gamma = (IonEk+IonEs)/IonEs;
        beta = sqrt(1e0-1e0/sqr(gamma));
        // Evaluate momentum compaction.
        R56 = -2e0*M_PI/(SampleLambda*IonEs*cube(beta*gamma))*L;
        // Convert from [m] and [eV/u] to [mm] and [MeV/u].
        R56 = -2e0*M_PI/(SampleLambda*IonEs*cube(beta*gamma))*L*1e-3;

        if (t_name == "marker") {
        } else if (t_name == "drift") {
            ElemPtr->transfer(state_t::PS_S, state_t::PS_PS) = R56;
            PrtMat(ElemPtr->transfer);
        } else if (t_name == "sbend") {
            ElemPtr->transfer(state_t::PS_S, state_t::PS_PS) = R56;
            PrtMat(ElemPtr->transfer);
        } else if (t_name == "quadrupole") {
        } else if (t_name == "solenoid") {
            Brho = c0*IonZ/(beta*IonW);
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
            PrtMat(ElemPtr->transfer);
        } else if (t_name == "rfcavity") {
            CavCnt++;
//            InitRFCav(conf, CavCnt, IonZ, IonEs, IonW, IonChargeState, EkState,
//                      Fy_absState, beta, gamma);
        }
        it++;
    } while (it != sim.p_elements.end());
}


void StateUnitFix(Machine &sim)
{
    Config                   D;
    std::auto_ptr<StateBase> state(sim.allocState(D));

//    std::cout << "# Machine configuration\n" << sim << "\n\n";

    // Propagate through first element (beam initial conditions).
    sim.propagate(state.get(), 0, 1);

    MatrixState* StatePtr = dynamic_cast<MatrixState*>(state.get());
//    std::cout << "\n" << *state;
    PrtMat(StatePtr->state);
    StatePtr->state(0, 0) = 1e0;
    PrtMat(StatePtr->state);
}


void PropagateState(Machine &sim)
{
    Machine::p_elements_t::iterator it;
    MomentElementBase               *ElemPtr;
    Config                          D;
    std::auto_ptr<StateBase>        state(sim.allocState(D));
    ElementVoid                     *elem;
    std::fstream                    outf;

    outf.open("test_jb.out", std::ofstream::out);

    //    std::cout << "# Machine configuration\n" << sim << "\n\n";

//    StateUnitFix(sim);
//    exit(1);

    for (it = sim.p_elements.begin(); it != sim.p_elements.end(); ++it) {
        elem = *it;

        ElemPtr = dynamic_cast<MomentElementBase *>(elem);
        assert(ElemPtr != NULL);

        PrtMat(ElemPtr->transfer);
//        sim.propagate(state.get(), elem->index, 1);
        elem->advance(*state);
        MatrixState* StatePtr = dynamic_cast<MatrixState*>(state.get());
    //    std::cout << "\n" << *state;
        PrtMat(StatePtr->state);
    }

    outf.close();
}



int main(int argc, char *argv[])
{
  Machine::p_elements_t::const_iterator it;

  const std::string HomeDir = "/home/johan/tlm_workspace/TLM_JB";

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

  std::vector<double> chgstate = conf->get<std::vector<double> >("ionChargeStates");
  std::cout << chgstate[0] << " " << chgstate[1] << "\n";

  // register state and element types
  registerLinear();
  registerMoment();

  Machine sim(*conf);
  sim.set_trace(&std::cout);

  CavData[0].RdData(HomeDir+"/data/axisData_41.txt");
  CavData[1].RdData(HomeDir+"/data/axisData_85.txt");

  // Turn trace on/off.
  if (false) sim.set_trace(&std::cout); else sim.set_trace(NULL);

  InitLong(sim);

//  InitLattice(sim);

  PropagateState(sim);

// for n states
//   get Chg_n, S_n
//   Initialize Lattice(Chg_n, S_n)
//   PropagateState(sim, S);
//

  //    it = sim.p_elements.begin();

  std::cout << "\nLattice parameters:\n";
  std::cout << "sim_type: " << sim.p_info.name
	    << "\n#Elements: " << sim.p_elements.size()
	    << "\n";

//      prt_lat(sim);

  //    std::cout<<"# Reduced lattice\n";
  //    GLPSPrint(std::cout, *conf);
  //    std::cout<<"\n";

  //    std::cerr<<"Generic AST:\n";
  //    std::cerr << *conf;
  //    std::cerr<<"GLPS:\n";
  //    GLPSPrint(std::cout, *conf);

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
