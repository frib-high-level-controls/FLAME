
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

#if true
    // Use [m, rad, m, rad, rad, eV/u].
    #define MeVtoeV 1e0
#else
    // Use [mm, rad, mm, rad, rad, MeV/u].
    #define MeVtoeV 1e6
#endif

// Global constansts.

// Speed of light [m/s].
# define c0           2.99792458e8
// Atomic mass unit [eV/c^2].
//# define AU           931.49432e6
// Vacuum permeability.
# define mu0          4e0*M_PI*1e-7
// Long. sampling frequency [Hz]; must be set to RF freq.
# define SampleFreq   80.5e6
# define SampleLambda c0/SampleFreq

// Phase space dimension.
# define PS_Dim       6

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
    IonEk = (LongTab.Ek[n-2]-StripperPara[2])*Stripper_E0Para[1] + Stripper_E0Para[0];
    IonW        = IonEk + IonEs;
    IonGamma    = IonW/IonEs;
    IonBeta     = sqrt(1e0-1e0/sqr(IonGamma));
    SampleIonK  = 2e0*M_PI/(IonBeta*SampleLambda);

    LongTab.set(LongTab.s[n-2], IonEk, LongTab.FyAbs[n-2], IonBeta, IonGamma);

    //            chargeAmount = fribstripper.chargeAmount_Baron;
}


void calGapModel(const double dis, const double IonW0, const double IonEs, const double IonFy0,
                 const double k, const double IonZ, const double lamda, const double Ecen,
                 const double T, const double S, const double Tp, const double Sp, const double V0,
                 double &IonW_f, double &IonFy_f)
{
    double Iongamma_f, IonBeta_f, k_f;

    IonW_f     = IonW0 + IonZ*V0*T*cos(IonFy0+k*Ecen) - IonZ*V0*S*sin(IonFy0+k*Ecen);
    Iongamma_f = IonW_f/IonEs;
    IonBeta_f  = sqrt(1e0-1e0/sqr(Iongamma_f));
//    k_f        = 2e0*M_PI/(IonBeta_f*lamda*1e3); // rad/mm
    k_f        = 2e0*M_PI/(IonBeta_f*lamda*1e3);

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
    IonW  = IonW0;
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
    IonEs = state->IonEs;
    IonW  = state->IonW;

    IonGamma   = IonW/IonEs;
    IonBeta    = sqrt(1e0-1e0/sqr(IonGamma));
    SampleIonK = 2e0*M_PI/(IonBeta*SampleLambda);

    std::cout << "\n" << "InitLong:" << "\n";
    std::cout << std::fixed << std::setprecision(5)
              << "  IonEs [Mev/u] = " << IonEs*1e-6 << ", IonEk [Mev/u] = " << state->IonEk*1e-6
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


void GetTransitFac(const int cavilabel, double beta, const int gaplabel, const double EfieldScl,
                   double &Ecen, double &T, double &Tp, double &S, double &Sp, double &V0)
{
    // Evaluate Electric field center, transit factors [T, T', S, S'] and cavity field.

    switch (cavilabel) {
    case 41:
        if (beta < 0.025 || beta > 0.08) {
            std::cerr << "*** GetTransitFac: beta out of Range" << "\n";
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
            V0   = 0.98477*EfieldScl;
            break;
        case 1:
            // Two gap calculation, first gap.
            Ecen = 0.0006384*pow(beta, -1.884) + 86.69;
            T    = -1.377e6*pow(beta, 5.0) + 4.316e5*pow(beta, 4.0) - 5.476e4*pow(beta, 3.0)
                   + 3570*pow(beta, 2.0) - 123.2*beta + 0.9232;
            Tp   = 2.277e7*pow(beta, 5.0) - 6.631e6 *pow(beta, 4.0) + 7.528e5*pow(beta, 3.0)
                   - 4.062e4*pow(beta, 2.0) + 924.7*beta + 1.699;
            S    = 0.0;
            Sp   = -1.335e6*pow(beta, 5.0) + 3.385e5*pow(beta, 4.0) - 2.98e4*pow(beta, 3.0)
                   + 806.6*pow(beta, 2.0) + 25.59*beta - 1.571;
            V0   = 0.492385*EfieldScl;
            break;
        case 2:
            // Two gap calculation, second gap.
            Ecen = -0.0006384*pow(beta, -1.884) + 33.31;
            T    = 1.377e6*pow(beta, 5.0) - 4.316e5*pow(beta, 4.0) + 5.476e4*pow(beta, 3.0)
                   - 3570*pow(beta, 2.0) + 123.2*beta - 0.9232;
            Tp   = -2.277e7*pow(beta, 5.0) + 6.631e6*pow(beta, 4.0) - 7.528e5*pow(beta, 3.0)
                   + 4.062e4*pow(beta, 2.0) - 924.7*beta - 1.699;
            S    = 0.0;
            Sp   = -1.335e6*pow(beta, 5.0) +  3.385e5*pow(beta, 4.0) - 2.98e4*pow(beta, 3.0)
                   + 806.6*pow(beta, 2.0) + 25.59*beta - 1.571;
            V0   = 0.492385*EfieldScl;
            break;
        default:
            std::cerr << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            exit(1);
        }
        break;
    case 85:
        if (beta < 0.05 || beta > 0.25) {
            std::cerr << "*** GetTransitFac: beta out of range" << "\n";
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
            Sp   = -2.755e8*pow(beta,7.0) + 3.109e8*pow(beta,6.0) - 1.462e8*pow(beta, 5.0)
                   + 3.691e7*pow(beta, 4.0) - 5.344e6*pow(beta, 3.0) + 4.315e5*pow(beta, 2.0)
                   - 1.631e4*beta + 162.7;
            V0   = 1.967715*EfieldScl;
            break;
        case 1:
            Ecen = 0.0002838*pow(beta, -2.13) + 76.5;
            T    = 0.0009467*pow(beta, -1.855) - 1.002;
            Tp   = -1.928e4*pow(beta, 5.0) + 2.195e4*pow(beta, 4.0) - 1.017e4*pow(beta, 3.0)
                   + 2468*pow(beta, 2.0) - 334*beta + 24.44;
            S    = 0.0;
            Sp   = -0.0009751*pow(beta, -1.898) + 0.001568;
            V0   = 0.9838574*EfieldScl;
            break;
        case 2:
            Ecen = -0.0002838*pow(beta, -2.13) + 73.5;
            T    = -0.0009467*pow(beta, -1.855) + 1.002;
            Tp   = 1.928e4*pow(beta, 5.0) - 2.195e4*pow(beta, 4.0) + 1.017e4*pow(beta, 3.0)
                   - 2468*pow(beta, 2.0) + 334*beta - 24.44;
            S    = 0.0;
            Sp   = -0.0009751*pow(beta, -1.898) + 0.001568;
            V0   = 0.9838574*EfieldScl;
            break;
        default:
            std::cerr << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            exit(1);
        }
        break;
    case 29:
        if (beta < 0.15 || beta > 0.4) {
            std::cerr << "*** GetTransitFac: beta out of range" << "\n";
            exit(1);
        }
        switch (gaplabel) {
        case 0:
            Ecen = 150.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = 76.54*pow(beta, 5.0) - 405.6*pow(beta, 4.0) + 486*pow(beta, 3.0)
                   - 248*pow(beta, 2.0) + 58.08*beta - 4.285;
            Sp   = -2.025e6*pow(beta,7.0) + 4.751e6*pow(beta,6.0) - 4.791e6*pow(beta, 5.0)
                   + 2.695e6*pow(beta, 4.0) - 9.127e5*pow(beta, 3.0) + 1.854e5*pow(beta, 2.0)
                   - 2.043e4*beta + 888;
            V0   = 2.485036*EfieldScl;
            break;
        case 1:
            Ecen = 0.01163*pow(beta, -2.001) + 91.77;
            T    = 0.02166*pow(beta, -1.618) - 1.022;
            Tp   = 1.389e4*pow(beta, 5.0) - 2.147e4*pow(beta, 4.0) + 1.313e4*pow(beta, 3.0)
                   - 3917*pow(beta, 2.0) + 534.7*beta - 11.25;
            S    = 0.0;
            Sp   = -454.4*pow(beta, 5.0) + 645.1*pow(beta, 4.0) - 343.9*pow(beta, 3.0)
                   + 78.77*pow(beta, 2.0) - 4.409*beta - 0.8283;
            V0   = 1.242518*EfieldScl;
        case 2:
            Ecen = -0.01163*pow(beta, -2.001) + 58.23;
            T    = -0.02166*pow(beta, -1.618) + 1.022;
            Tp   = -1.389e4*pow(beta, 5.0) + 2.147e4*pow(beta, 4.0) - 1.313e4*pow(beta, 3.0)
                   + 3917*pow(beta, 2.0) - 534.7*beta + 11.25;
            S    = 0.0;
            Sp   = -454.4*pow(beta, 5.0) + 645.1*pow(beta, 4.0) - 343.9*pow(beta, 3.0)
                   + 78.77*pow(beta, 2.0) - 4.409*beta - 0.8283;
            V0   = 1.242518*EfieldScl;
            break;
        default:
            std::cerr << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            exit(1);
        }
        break;
    case 53:
        if (beta < 0.3 || beta > 0.6) {
            std::cerr << "*** GetTransitFac: beta out of range" << "\n";
            exit(1);
        }
        switch (gaplabel) {
        case 0:
            Ecen = 250.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = -52.93*pow(beta, 5.0) + 84.12*pow(beta, 4.0) - 17.73*pow(beta, 3.0)
                   - 38.49*pow(beta, 2.0) + 26.64*beta - 4.222;
            Sp   = -4.167e4*pow(beta, 5.0) + 1.075e5*pow(beta, 4.0) - 1.111e5*pow(beta, 3.0)
                   + 5.702e4*pow(beta, 2.0) - 1.413e4*beta - 1261;
            V0   = 4.25756986*EfieldScl;
            break;
        case 1:
            Ecen = 0.01219*pow(beta, -2.348) + 137.8;
            T    = 0.04856*pow(beta, -1.68) - 1.018;
            Tp   = 1641*pow(beta, 5.0) - 4109*pow(beta, 4.0) + 4081*pow(beta, 3.0)
                   - 1973*pow(beta, 2.0) + 422.8*beta - 3.612;
            S    = 0.0;
            Sp   = -0.03969*pow(beta, -1.775) + 0.009034;
            V0   = 2.12878493*EfieldScl;
            break;
        case 2:
            Ecen = -0.01219*pow(beta, -2.348) + 112.2;
            T    = -0.04856*pow(beta, -1.68) + 1.018;
            Tp   = -1641*pow(beta, 5.0) + 4109*pow(beta, 4.0) - 4081*pow(beta, 3.0)
                   + 1973*pow(beta, 2.0) - 422.8*beta + 3.612;
            S    = 0.0;
            Sp   = -0.03969*pow(beta, -1.775) + 0.009034;
            V0   = 2.12878493*EfieldScl;
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
    Ecen *= 1e-3;
}


void CavityMatrix(double dis, double EfieldScl, double TTF_tab[], double beta_tab[], double gamma_tab[],
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
//                V0     = fribnode.attribute[1]*EfieldScl;
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
//                V0     = fribnode.attribute[1]*EfieldScl;
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
//                V0     = fribnode.attribute[1]*EfieldScl;
//                T      = fribnode.attribute[2];
//                S      = fribnode.attribute[3];
//                dpy    = ionZ*V0/(beta*beta)/gamma/FRIBPara.ionA/FRIBPara.AU
//                         *(T*Math.cos(ionFy)-S*Math.sin(ionFy));
//                Mprob  = PhaseMatrix.identity();
//                Mprob.setElem(3, 6, dpy);
//                Mtrans = Mprob.times(Mtrans);
//            } else if (fribnode.label == "EQuad") {
//                if (multipoleLevel<2) break;
//                V0     = fribnode.attribute[1]*EfieldScl;
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
//                V0     = fribnode.attribute[1]*EfieldScl;
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
//                V0     = fribnode.attribute[1]*EfieldScl;
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
//                V0     = fribnode.attribute[1]*EfieldScl;
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


void CavTLMMat(const int cavilabel, const double Rm, const double IonZ, const double IonEs,
                  const double EfieldScl, const double IonFyi_s, const double IonEk_s,
                  const double IonLamda)
{
    int    k;
    double IonWi_s, gammai_s, betai_s, IonKi_s;
    double Ecen[2], T[2], Tp[2], S[2], Sp[2], V0[2];
    double dis, IonW_s[2], IonFy_s[2], gamma_s[2], beta_s[2], IonK_s[2];

    const int NGaps = 1;

    IonWi_s  = IonEk_s + IonEs;
    gammai_s = IonWi_s/IonEs;
    betai_s  = sqrt(1e0-1e0/sqr(gammai_s));
    IonKi_s  = 2e0*M_PI/(betai_s*IonLamda);

    GetTransitFac(cavilabel, betai_s, NGaps, EfieldScl, Ecen[0], T[0], Tp[0], S[0], Sp[0], V0[0]);
//    double dis = (axisData.get(axisData.size()-1)[0]-axisData.get(0)[0])/2; // should divided by 2!!!
    calGapModel(dis, IonWi_s, IonEs, IonFyi_s, IonKi_s, IonZ, IonLamda,
                Ecen[0], T[0], S[0], Tp[0], Sp[0], V0[0], IonW_s[0], IonFy_s[0]);
    gamma_s[0] = IonW_s[0]/IonEs;
    beta_s[0]  = sqrt(1e0-1e0/sqr(gamma_s[0]));
    IonK_s[0]  = 2*M_PI/(beta_s[0]*IonLamda);

    GetTransitFac(cavilabel, betai_s, NGaps, EfieldScl, Ecen[1], T[1], Tp[1], S[1], Sp[1], V0[1]);
    calGapModel(dis, IonWi_s, IonEs, IonFy_s[1], IonK_s[1], IonZ, IonLamda,
                Ecen[1], T[1], S[1], Tp[1], Sp[1], V0[1], IonW_s[1], IonFy_s[1]);
    gamma_s[1] = IonW_s[1]/IonEs;
    beta_s[1]  = sqrt(1e0-1e0/sqr(gamma_s[1]));
    IonK_s[1]  = 2*M_PI/(beta_s[1]*IonLamda);

//    Ecen_1 = Ecen_1-dis; // Switch Ecen_1 axis centered at cavity center// dis/2 to dis 14/12/12

//    double[] TTF_tab = {Ecen_1, T_1, Tp_1, S_1, Sp_1, V0_1, Ecen_2, T_2, Tp_2, S_2, Sp_2, V0_2};
//    double[] beta_tab = {betai_s, betac_s, betaf_s};
//    double[] gamma_tab = {gammai_s, gammac_s, gammaf_s};
//    double IonK_1 = (IonKi_s+IonKc_s)/2;
//    double IonK_2 = (IonKc_s+IonKf_s)/2;		    // Bug found 2015-01-07

//    PhaseMatrix matrix = PhaseMatrix.identity();
//    util.ArrayList<TlmNode> thinlenLine  =  new util.ArrayList<TlmNode>();
//    thinlenLine = calCavity_thinlenLine(beta_tab, gamma_tab, cavi, IonK_1, IonK_2);
//    matrix = calCavity_Matrix(dis, EfieldScl, TTF_tab, beta_tab, gamma_tab, IonLamda, IonZ, IonFyi_s,
//    IonFyc_s, thinlenLine, Rm);

//    return matrix;
}


void InitRFCav(const Config &conf, const int CavCnt, const double IonZ, const double IonEs, double &IonW,
               const double ChgState,
               double &EkState, double &Fy_absState, double &beta, double &gamma)
{
    std::string CavType;
    int         cavi, cavilabel, multip;
    double      Rm, IonFy_i, Ek_i, avebeta, avegamma, fRF, CaviIonK, SampleIonK, EfieldScl;
    double      IonW_o, IonFy_o, accIonW;

    CavType = conf.get<std::string>("cavtype");
    if (CavType == "0.041QWR") {
        cavi      = 1;
        cavilabel = 41;
        multip    = 1;
        Rm        = 17e0;
    } else if (conf.get<std::string>("cavtype") == "0.085QWR") {
        cavi      = 2;
        cavilabel = 85;
        multip    = 1;
        Rm        = 17e0;
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

    GetCavBoost(CavData[cavi-1], IonW, IonFy_i, CaviIonK, ChgState, IonEs,
                c0/fRF, EfieldScl, IonW_o, IonFy_o);

    accIonW      = IonW_o - IonW;
    IonW         = IonW_o;
    EkState      = IonW - IonEs;
    IonW         = EkState + IonEs;
    gamma        = IonW/IonEs;
    beta         = sqrt(1e0-1e0/sqr(gamma));
    avebeta      = (avebeta+beta)/2e0;
    avegamma     = (avegamma+gamma)/2e0;
    Fy_absState += (IonFy_o-IonFy_i)/multip;

//    CaviMatrix = calCavityTLM_Matrix(Rm, ChgState, EfieldScl, IonFy_i, Ek_i, caviLamda, cavi);
    CavTLMMat(cavilabel, Rm, ChgState, IonEs, EfieldScl, IonFy_i, Ek_i, c0/fRF);

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


void InitLattice(Machine &sim, const double ChgState, const std::vector<double> BaryCenter)
{
    // Evaluate transport matrices for given beam initial conditions.
    typedef MatrixState state_t;

    std::stringstream               strm;
    int                             CavCnt;
    double                          IonW, IonEs, IonEk, IonZ, IonLambda, s, L, beta, gamma;
    double                          SampleionK, R56, Brho, K, Ek_ini, EkState, Fy_absState;
    Machine::p_elements_t::iterator it;
    MomentElementBase               *ElemPtr;
//    LinearElementBase<MatrixState>  *ElemPtr;

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
    Fy_absState = BaryCenter[state_t::PS_S];
    EkState     = IonEk + BaryCenter[state_t::PS_PS];

//    TransVector[state] = get(k).centerVector;
//    ThetaMatrix[state] = get(k).thetaMatrix;
//    Fy_absState      = get(k).centerVector.getElem(4);

    std::cout << "\n" << "InitLattice:" << "\n";
    std::cout << std::fixed << std::setprecision(5)
              << "  IonEs [Mev/u] = " << IonEs*1e-6 << ", IonEk [Mev/u] = " << state->IonEk*1e-6
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

        if (t_name != "marker") {
            L = conf.get<double>("L");
            s += L;
        }

        gamma      = (EkState+IonEs)/IonEs;
        beta       = sqrt(1e0-1e0/sqr(gamma));
        SampleionK = 2e0*M_PI/(beta*SampleLambda);

        // Evaluate momentum compaction.
//        R56 = -2e0*M_PI/(SampleLambda*IonEs*cube(beta*gamma))*L;
        // Convert from [m] and [eV/u] to [mm] and [MeV/u].
        R56 = -2e0*M_PI/(SampleLambda*IonEs*cube(beta*gamma))*L*1e6; // Convert from [eV/u] to [MeV/u].

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
            Brho = beta*(EkState+IonEs)/(c0*IonZ);
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
            InitRFCav(conf, CavCnt, IonZ, IonEs, IonW, ChgState, EkState,
                      Fy_absState, beta, gamma);
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
            if (j % 2 == 0) StatePtr->state(j, k) *= 1e-3;
            if (k % 2 == 0) StatePtr->state(j, k) *= 1e-3;
        }

    for (j = 0; j < PS_Dim; j++)
        StatePtr->state(j, PS_Dim-1) *= 1e6;

    std::cout << "\n";
    PrtMat(StatePtr->state);
}


int main(int argc, char *argv[])
{
    typedef MatrixState state_t;

    int                 k;
    std::vector<double> ChgState;
    std::vector<double> BaryCenter[2];

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
    // Change units from [mm, rad, mm, rad, rad, MeV/u] to [m, rad, m, rad, rad, eV/u].
    BaryCenter[0][state_t::PS_X]  *= 1e-3;
    BaryCenter[0][state_t::PS_Y]  *= 1e-3;
    BaryCenter[0][state_t::PS_PS] *= 1e6;
    BaryCenter[1][state_t::PS_X]  *= 1e-3;
    BaryCenter[1][state_t::PS_Y]  *= 1e-3;
    BaryCenter[1][state_t::PS_PS] *= 1e6;


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
