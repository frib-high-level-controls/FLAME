
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <math.h>

#include "scsi/config.h"

#include <scsi/base.h>
//#include <scsi/linear.h>
#include <scsi/moment.h>
#include <scsi/state/vector.h>
#include <scsi/state/matrix.h>


extern int glps_debug;


typedef boost::numeric::ublas::matrix<double> value_t;


// Global constants and parameters.
const double c0           = 2.99792458e8,   // Speed of light [m/s].
             u            = 931.49432e6,    // Atomic mass unit [eV/c^2].
             mu0          = 4e0*M_PI*1e-7,  // Vacuum permeability.
             SampleFreq   = 80.5e6,         // Long. sampling frequency [Hz];
                                           // must be set to RF freq.
             SampleLambda = c0/SampleFreq;

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


const int MaxTab = 2000; // Max number of entries.

class TableType {
// Table for longitudinal initialization for reference particle.
public:
    int    n;         // Length of table.
    double s[MaxTab], // Longitudinal position [m].
    Ek[MaxTab],       // Kinetic energy [eV/u].
    FyAbs[MaxTab],    // Synchrotron phase [rad].
    Beta[MaxTab],     // Relativistic factor beta.
    Gamma[MaxTab];    // Relativistic factor gamma.

    void zero(const int);
    void set(const int, const double, const double, const double,
             const double, const double);
    void show(std::ostream& strm, const int) const;
};


const int MaxCavData = 1000; // Max number of data points.

class CavDataType {
// Cavity on-axis longitudinal electric field vs. s.
public:
    int    n;                 // Lenght of data.
    double s[MaxCavData],     // s coordinate [m]
    Elong[MaxCavData]; // Longitudinal Electric field [V/m].

    void RdData(const std::string);
    void show(std::ostream& strm, const int) const;
};


// Global variables.
CavDataType CavData[2];
TableType   LongTab;


void TableType::zero(const int k)
{
    this->s[k] = 0e0; this->Ek[k] = 0e0; this->FyAbs[k] = 0e0;
    this->Beta[k] = 0e0; this->Gamma[k] = 0e0;
}


void TableType::set(const int k, const double s, const double Ek,
                    const double FyAbs, const double Beta, const double Gamma)
{
    this->s[k] = s; this->Ek[k] = Ek; this->FyAbs[k] = FyAbs; this->Beta[k]  = Beta;
    this->Gamma[k] = Gamma;
}


void TableType::show(std::ostream& strm, const int k) const
{
    strm << std::scientific << std::setprecision(10)
         << std::setw(18) << this->s[k]
         << std::setw(18) << this->Ek[k]*1e-6
         << std::setw(18) << this->FyAbs[k]
         << std::setw(18) << this->Beta[k]
         << std::setw(18) << this->Gamma[k];
}


void CavDataType::RdData(const std::string FileName)
{
    int          k;
    std::string  line;
    std::fstream inf;

    inf.open(FileName.c_str(), std::ifstream::in);
    if (!inf.is_open()) {
        std::cerr << "Failed to open " << FileName << "\n";
        exit(1);
    }
    this->n = 0;
    while (getline(inf, line)) {
        this->n++;
        if (n > MaxCavData) {
            std::cout << "*** RdData: max data exceeded";
            exit(1);
        }
        sscanf(line.c_str(), "%lf %lf",
               &this->s[this->n-1], &this->Elong[this->n-1]);
        // Convert from [mm] to [m].
        this->s[this->n-1] *= 1e-3;
    }
    inf.close();

    if (false) {
        std::cout << "\n";
        for (k = 0; k < this->n; k++) {
            this->show(std::cout, k);
            std::cout << "\n";
        }
    }
}


void CavDataType::show(std::ostream& strm, const int k) const
{
    strm << std::scientific << std::setprecision(5)
         << std::setw(13) << this->s[k] << std::setw(13) << this->Elong[k];
}


void prt_lat(Machine &sim)
{
  std::stringstream strm;

  std::cout << "# Machine configuration\n" << sim << "\n\n";

  strm << "sim_type: " << sim.p_info.name << "\n#Elements: "
       << sim.p_elements.size() << "\n";
  for(Machine::p_elements_t::const_iterator it = sim.p_elements.begin(),
        end = sim.p_elements.end(); it != end; ++it)
      (*it)->show(strm);
  std::cout << strm;
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


double calFindPhaseTable_simplify(const int cavi, const double IonEk,
                                  const double IonFys, const double FyAbs,
                                  const double multip)
{
    // If the cavity is not at full power, the method gives synchrotron
    // phase slightly different from the nominal value.

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
        std::cout << "*** calFindPhaseTable_simplify: undef. cavity type"
                  << "\n";
    }

    return IonFys - Fyc - FyAbs*multip;
}

void calGapTrace(const CavDataType &CavData, const double IonW0,
                 const double IonFy0, const double IonK0, const double IonZ,
                 const double IonEs, const double IonLamda,
                 const double E_fac_ad, double &IonW, double &IonFy)
{
    int    n = CavData.n,
           k;

    double dis = CavData.s[n-1] - CavData.s[0],
           dz  = dis/(n-1),
           IonK, IonFylast, Iongamma, IonBeta;

    IonFy = IonFy0;
    IonK  = IonK0;
    for (k = 0; k < n-1; k++) {
        IonFylast = IonFy;
        IonFy += IonK*dz;
        IonW  += IonZ*E_fac_ad*(CavData.Elong[k]+CavData.Elong[k+1])/2e0
                *cos((IonFylast+IonFy)/2e0)*dz; // MeV
        Iongamma = IonW/IonEs;
        IonBeta = sqrt(1e0-1e0/sqr(Iongamma));
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
    // Use Baron's formulafor carbon foil.
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
    double      fRF, multip, caviIonK, IonFys, E_fac_ad, caviFy, IonFy_i, IonFy_o, Iongamma;

    CavType = conf.get<std::string>("cavtype");
    if (CavType == "0.041QWR") {
        cavi = 1;
    } else if (conf.get<std::string>("cavtype") == "0.085QWR") {
        cavi = 2;
    } else {
        std::cout << "*** InitLong: undef. cavity type: "
                  << CavType << "\n";
        exit(1);
    }

    fRF      = conf.get<double>("f");
    multip   = fRF/SampleFreq;
    caviIonK = 2e0*M_PI/(IonBeta*c0/fRF);
    IonFys   = conf.get<double>("phi")*M_PI/180e0; // Synchrotron phase.
    E_fac_ad = conf.get<double>("scl_fac");        // Electric field scale factor.

    caviFy = calFindPhaseTable_simplify(
                cavi, IonW-IonEs, IonFys, LongTab.FyAbs[n-1], multip);

    IonFy_i = multip*LongTab.FyAbs[n-1] + caviFy;
    // Evaluate change of reference particle kinetic energy,
    // absolute phase, beta, and gamma.
    calGapTrace(CavData[cavi-1], IonW, IonFy_i, caviIonK, IonZ,
                IonEs, c0/fRF, E_fac_ad, IonW, IonFy_o);
    Iongamma   = IonW/IonEs;
    IonBeta    = sqrt(1e0-1e0/sqr(Iongamma));
    SampleIonK = 2e0*M_PI/(IonBeta*SampleLambda);

    LongTab.set(n-1, LongTab.s[n-2]+conf.get<double>("L"), IonW-IonEs,
            LongTab.FyAbs[n-2]+(IonFy_o-IonFy_i)/multip, IonBeta, Iongamma);
}


void PropagateLongStripper(const Config &conf, const int n, double &IonZ, const double IonEs,
                           double &IonW, double &SampleIonK, double &IonBeta)
{
    double IonEk, Iongamma;
    double chargeAmount_Baron[Stripper_n];

    IonZ = Stripper_IonZ;
    ChargeStripper(Stripper_IonMass, Stripper_IonProton, IonBeta,
                   Stripper_n, Stripper_IonChargeStates,
                   chargeAmount_Baron);
    // Evaluate change in reference particle energy due to stripper
    // model energy straggling.
    IonEk = (LongTab.Ek[n-2]-StripperPara[2])*Stripper_E0Para[1] + Stripper_E0Para[0];
    IonW        = IonEk + IonEs;
    Iongamma    = IonW/IonEs;
    IonBeta     = sqrt(1e0-1e0/sqr(Iongamma));
    SampleIonK  = 2e0*M_PI/(IonBeta*SampleLambda);

    LongTab.set(n-1, LongTab.s[n-2], IonEk, LongTab.FyAbs[n-2], IonBeta, Iongamma);

    //            chargeAmount = fribstripper.chargeAmount_Baron;
}


void InitLong(Machine &sim)
{
    // Longitudinal initialization for reference particle.
    // Evaluate beam energy and cavity loaded phase along the lattice.
    int                             n;
    double                          Iongamma, IonBeta, SampleIonK;
    double                          IonW, IonZ, IonEs;
    Machine::p_elements_t::iterator it;

    Config                   D;
    std::auto_ptr<StateBase> state(sim.allocState(D));
    // Propagate through first element.
    sim.propagate(state.get(), 0, 1);

    IonZ  = state->IonZ;
    IonEs = state->IonEs;
    IonW  = state->IonW;

    Iongamma   = IonW/IonEs;
    IonBeta    = sqrt(1e0-1e0/sqr(Iongamma));
    SampleIonK = 2e0*M_PI/(IonBeta*SampleLambda);

    std::cout << "\n" << "InitLong:" << "\n\n";
    std::cout << std::scientific << std::setprecision(5)
              << "IonEs [Mev/u] = " << IonEs*1e-6 << ", IonEk [Mev/u] = " << state->IonEk*1e-6
              << ", IonW [Mev/u] = " << IonW*1e-6 << "\n";

    n = 1;
    LongTab.set(n-1, 0e0, state->IonEk, 0e0, IonBeta, Iongamma);

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
            LongTab.set(n-1, LongTab.s[n-2]+conf.get<double>("L"), LongTab.Ek[n-2],
                    LongTab.FyAbs[n-2]+SampleIonK*conf.get<double>("L"),
                    LongTab.Beta[n-2], LongTab.Gamma[n-2]);
        } else if (t_name == "sbend") {
            n++;
            LongTab.set(n-1, LongTab.s[n-2]+conf.get<double>("L"), LongTab.Ek[n-2],
                    LongTab.FyAbs[n-2]+SampleIonK*conf.get<double>("L"),
                    LongTab.Beta[n-2], LongTab.Gamma[n-2]);
        } else if (t_name == "quadrupole") {
            n++;
            LongTab.set(n-1, LongTab.s[n-2]+conf.get<double>("L"), LongTab.Ek[n-2],
                    LongTab.FyAbs[n-2]+SampleIonK*conf.get<double>("L"),
                    LongTab.Beta[n-2], LongTab.Gamma[n-2]);
        } else if (t_name == "solenoid") {
            n++;
            LongTab.set(n-1, LongTab.s[n-2]+conf.get<double>("L"), LongTab.Ek[n-2],
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
            std::cout << "\n";
        }
        it++;
    } while (it != sim.p_elements.end());
}


void PrtMat(const value_t &M)
{
    int j, k;

    const int n = 6;

    std::cout << "\n";
    for (j = 0; j < n; j++) {
        for (k = 0; k < n; k++)
            std::cout << std::scientific << std::setprecision(5)
                      << std::setw(13) << M(j, k);
        std::cout << "\n";
    }
}


void InitLattice(Machine &sim)
{
    // Evaluate transport matrices for given beam intial conditions.
    typedef MatrixState state_t;

    std::stringstream               strm;
    double                          IonW, IonEs, IonEk, IonZ, IonLambda, s, L, beta, gamma;
    double                          R56, Brho, K;
    Machine::p_elements_t::iterator it;
    MomentElementBase *ElemPtr;
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

    std::cout << "\n" << "InitLattice:" << "\n\n";
//    std::cout << *state << "\n";
    std::cout << std::scientific << std::setprecision(5)
              << "IonEs [Mev/u] = " << IonEs*1e-6 << ", IonEk [Mev/u] = " << state->IonEk*1e-6
              << ", IonW [Mev/u] = " << IonW*1e-6 << ", IonLambda [m^-1] = " << IonLambda << "\n";

    s = 0e0;
    it = sim.p_elements.begin();
    // Skip over state.
    it++;
    do {
        ElementVoid*  elem   = *it;
        const Config& conf   = elem->conf();
        Config        newconf(conf);
        std::string   t_name = elem->type_name(); // C string -> C++ string.

//        ElemPtr = dynamic_cast<LinearElementBase<MatrixState> *>(elem);
        ElemPtr = dynamic_cast<MomentElementBase *>(elem);
        assert(ElemPtr != NULL);

        std::cout << t_name << "\n";

        if (t_name != "marker") {
            L = conf.get<double>("L");
            s += L;
        }

        gamma = (IonEk+IonEs)/IonEs;
        beta = sqrt(1e0-1e0/sqr(gamma));
        // Evaluate momentum compaction.
        R56 = -2e0*M_PI/(SampleLambda*IonEs*cube(beta*gamma))*L;

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
        }
        it++;
    } while (it != sim.p_elements.end());
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

  InitLattice(sim);

  //    it = sim.p_elements.begin();

  std::cout << "\nLattice parameters:\n";
  std::cout << "sim_type: " << sim.p_info.name
	    << "\n#Elements: " << sim.p_elements.size()
	    << "\n";

  //    prt_lat(sim);

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
