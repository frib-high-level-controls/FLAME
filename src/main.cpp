
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <math.h>

#include "scsi/config.h"

#include <scsi/base.h>
#include <scsi/linear.h>
#include <scsi/state/vector.h>
#include <scsi/state/matrix.h>


extern int glps_debug;


const int MaxTab = 2000; // Max number of entries.

class TableType {
// Table for longitudinal initialization for reference particle.
public:
    int    n;              // Length of table.
    double s[MaxTab],      // Longitudinal position [m].
           Ek[MaxTab],     // Kinetic energy [eV/u].
           Fy_abs[MaxTab], // Synchrotron phase [rad].
           Beta[MaxTab],   // Relativistic factor beta.
           Gamma[MaxTab];  // Relativistic factor gamma.

    void zero(const int);
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


void TableType::zero(const int k)
{
    this->s[k] = 0e0; this->Fy_abs[k] = 0e0; this->Ek[k] = 0e0; this->Beta[k] = 0e0;
    this->Gamma[k] = 0e0;
}


void TableType::show(std::ostream& strm, const int k) const
{
    strm << std::scientific << std::setprecision(10)
         << std::setw(18) << this->s[k]
         << std::setw(18) << this->Ek[k]*1e-6
         << std::setw(18) << this->Fy_abs[k]
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
        std::cerr << "Failed to open " << FileName << std::endl;
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
        std::cout << std::endl;
        for (k = 0; k < this->n; k++) {
            this->show(std::cout, k);
            std::cout << std::endl;
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

    strm << "sim_type: " << sim.p_info.name << "\n#Elements: " << sim.p_elements.size() << "\n";
    for(Machine::p_elements_t::const_iterator it = sim.p_elements.begin(),
        end = sim.p_elements.end(); it != end; ++it)
    {
        (*it)->show(strm);
    }
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


double calFindPhaseTable_simplify(const int cavi, const double ionEk, const double ionFys,
                                  const double Fy_abs, const double multip)
{
    // Note, if the cavity is not at full power, the method gives synchrotron phase slightly
    // different from the value at set point.

    double Fy_c;

    // ionEk = ionW - ionEs

    switch (cavi) {
    case 1:
        Fy_c = 4.394*pow(ionEk*1e-6, -0.4965) - 4.731;
        break;
    case 2:
        Fy_c = 5.428*pow(ionEk*1e-6, -0.5008) + 1.6;
        break;
    case 3:
        Fy_c = 22.35*pow(ionEk*1e-6, -0.5348) + 2.026;
        break;
    case 4:
        Fy_c = 41.43*pow(ionEk*1e-6, -0.5775) + 2.59839;
        break;
    case 5:
        Fy_c = 5.428*pow(ionEk*1e-6, -0.5008) + 1.6;
        break;
    default:
        std::cout << "*** calFindPhaseTable_simplify: undef. cavity type" << std::endl;
    }

    return ionFys - Fy_c - Fy_abs*multip;
}

void calGapTrace(const CavDataType &CavData, const double ionW0, const double ionFy0,
                 const double ionK0, const double ionZ, const double ionEs, const double ionLamda,
                 const double E_fac_ad, double &ionW, double &ionFy)
{
    int    n = CavData.n,
           k;

    double dis   = CavData.s[n-1] - CavData.s[0],
           dz    = dis/(n-1),
           ionK, ionFylast, iongamma, ionBeta;

    ionFy = ionFy0;
    ionK  = ionK0;
    for (k = 0; k < n-1; k++) {
        ionFylast = ionFy;
        ionFy += ionK*dz;
        ionW  += ionZ*E_fac_ad*(CavData.Elong[k]+CavData.Elong[k+1])/2e0
                 *cos((ionFylast+ionFy)/2e0)*dz; // MeV
        iongamma = ionW/ionEs;
        ionBeta = sqrt(1e0-1e0/sqr(iongamma));
        if ((ionW-ionEs) < 0e0) {
            ionW = ionEs;
            ionBeta = 0e0;
        }
        ionK = 2e0*M_PI/(ionBeta*ionLamda);
    }
}


double calGauss(double in, const double Q_ave, const double d)
{
    // Gaussian distribution.
    return 1e0/sqrt(2e0*M_PI)/d*exp(-0.5e0*sqr(in-Q_ave)/sqr(d));
}


void calStripperCharge(const double ionProton, const double beta,
                       double &Q_ave, double &d)
{
    // Use Baron's formulafor carbon foil.
    double Q_ave1, Y;

    Q_ave1 = ionProton*(1e0-exp(-83.275*(beta/pow(ionProton, 0.447))));
    Q_ave  = Q_ave1*(1e0-exp(-12.905+0.2124*ionProton-0.00122*sqr(ionProton)));
    Y      = Q_ave1/ionProton;
    d      = sqrt(Q_ave1*(0.07535+0.19*Y-0.2654*sqr(Y)));
}


void ChargeStripper(const double ionMass, const double ionProton, const double beta,
                    const int nChargeStates, const double ionChargeStates[],
                    double chargeAmount_Baron[])
{
    int    k;
    double Q_ave, d;

    calStripperCharge(ionProton, beta, Q_ave, d);
    for (k = 0; k < nChargeStates; k++)
        chargeAmount_Baron[k] = calGauss(ionChargeStates[k]*ionMass, Q_ave, d);
}


void init_long(Machine &sim, const CavDataType CavData[])
{
    // Longitudinal initialization for reference particle.
    std::string    CavType;
    int       n, cavi;
    double    iongamma, ionBeta, SampleLamda, SampleionK, ionFy_i;
    double    multip, caviionK, ionFys, E_fac_ad, fRF, caviFy, ionW, ionZ, ionEs, ionFy_o, ionEk1;
    TableType tab;

    const double c0         = 2.99792458e8,   // Speed of light [m/s].
                 u          = 931.49432e6,    // Atomic mass unit [eV/c^2].
                 mu0        = 4e0*M_PI*1e-7,  // Vacuum permeability.
                 SampleFref = 80.5e6;         // Long. sampling frequency [Hz]; must be set to RF freq.

    // Charge stripper parameters.
    const int    Stripper_n                 = 5;
    const double Stripper_ionZ              = 78e0/238e0,
                 Stripper_ionMass           = 238e0,
                 Stripper_ionProton         = 92e0,
                 Stripper_ionChargeStates[] = {76e0/238e0, 77e0/238e0, 78e0/238e0, 79e0/238e0, 80e0/238e0},
                 Stripper_E0Para[]          = {16.348e6, 1.00547, -0.10681}, // Unit for last ???
                 StripperPara[]             = {3e0, 20e0, 16.623e6}; // thickness [microns],
                                                                      // thickness variation [%],
                                                                      // reference energy [eV/u].

    double chargeAmount_Baron[Stripper_n];

    Config D;
    std::auto_ptr<StateBase> state(sim.allocState(D));
    // Propagate through first element.
    sim.propagate(state.get(), 0, 1);

    ionW  = state->ionW;
    ionEs = state->ionEs;
    ionZ  = state->ionZ;

    iongamma    = ionW/ionEs;
    ionBeta     = sqrt(1e0-1e0/sqr(iongamma));
    SampleLamda = c0/SampleFref;
    SampleionK  = 2e0*M_PI/(ionBeta*SampleLamda); // [m^-1].

    std::cout << *state << std::endl;
    std::cout << std::scientific << std::setprecision(5)
              << "ionEs = " << ionEs  << ", ionEk = " << state->ionEk
              << ", ionW = " << ionW << std::endl << std::endl;

    n               = 1;
    tab.zero(n-1);
    tab.Ek[n-1]     = state->ionEk;
    tab.Beta[n-1]   = ionBeta;
    tab.Gamma[n-1]  = iongamma;

    for(Machine::p_elements_t::iterator it = sim.p_elements.begin(),
        end = sim.p_elements.end(); it != end; ++it)
    {
        ElementVoid* elem = *it;
        const Config & conf = elem->conf();
        std::string t_name = elem->type_name(); // C string -> C++ string.

//        std::cout << t_name << std::endl;

        if (t_name == "marker") {
        } else if (t_name == "drift") {
            n++;
            tab.s[n-1]      = tab.s[n-2] + conf.get<double>("L");
            tab.Fy_abs[n-1] = tab.Fy_abs[n-2] + SampleionK*conf.get<double>("L");
            tab.Ek[n-1]     = tab.Ek[n-2];
            tab.Beta[n-1]   = tab.Beta[n-2];
            tab.Gamma[n-1]  = tab.Gamma[n-2];
        } else if (t_name == "sbend") {
            n++;
            tab.s[n-1]      = tab.s[n-2] + conf.get<double>("L");
            tab.Fy_abs[n-1] = tab.Fy_abs[n-2] + SampleionK*conf.get<double>("L");
            tab.Ek[n-1]     = tab.Ek[n-2];
            tab.Beta[n-1]   = tab.Beta[n-2];
            tab.Gamma[n-1]  = tab.Gamma[n-2];
        } else if (t_name == "quadrupole") {
            n++;
            tab.s[n-1]      = tab.s[n-2] + conf.get<double>("L");
            tab.Fy_abs[n-1] = tab.Fy_abs[n-2] + SampleionK*conf.get<double>("L");
            tab.Ek[n-1]     = tab.Ek[n-2];
            tab.Beta[n-1]   = tab.Beta[n-2];
            tab.Gamma[n-1]  = tab.Gamma[n-2];
        } else if (t_name == "solenoid") {
            n++;
            tab.s[n-1]      = tab.s[n-2] + conf.get<double>("L");
            tab.Fy_abs[n-1] = tab.Fy_abs[n-2] + SampleionK*conf.get<double>("L");
            tab.Ek[n-1]     = tab.Ek[n-2];
            tab.Beta[n-1]   = tab.Beta[n-2];
            tab.Gamma[n-1]  = tab.Gamma[n-2];
        } else if (t_name == "rfcavity") {
            n++;
            CavType = conf.get<std::string>("cavtype");
            if (CavType == "0.041QWR") {
                cavi = 1;
            } else if (conf.get<std::string>("cavtype") == "0.085QWR") {
                cavi = 2;
            } else {
                std::cout << "*** init_long: undef. cavity type: " << CavType << std::endl;
                exit(1);
            }
            fRF      = conf.get<double>("f");
            multip   = fRF/SampleFref;
            caviionK = 2e0*M_PI/(ionBeta*c0/fRF);
             // Synchrotron phase.
            ionFys   = conf.get<double>("phi")*M_PI/180e0;
            // Electric field scale factor.
            E_fac_ad = conf.get<double>("scl_fac");

            caviFy = calFindPhaseTable_simplify(
                        cavi, ionW-ionEs, ionFys, tab.Fy_abs[n-1], multip);

            ionFy_i = multip*tab.Fy_abs[n-1] + caviFy;
            // Evaluate change of reference particle kinetic energy, absolute phase, beta,
            // and gamma.
            calGapTrace(CavData[cavi-1], ionW, ionFy_i, caviionK, ionZ, ionEs, c0/fRF, E_fac_ad,
                        ionW, ionFy_o);
            iongamma = ionW/ionEs;
            ionBeta = sqrt(1e0-1e0/sqr(iongamma));
            SampleionK = 2e0*M_PI/(ionBeta*SampleLamda);
            tab.s[n-1] = tab.s[n-2] + conf.get<double>("L");
            tab.Ek[n-1] = ionW - ionEs;
            tab.Fy_abs[n-1] = tab.Fy_abs[n-2] + (ionFy_o-ionFy_i)/multip;
            tab.Beta[n-1] = ionBeta;
            tab.Gamma[n-1] = iongamma;
        } else if (t_name == "stripper") {
            // For the longitudinal plane evaluate change in reference and multi-charge states,
            // and charge amounts.
            n++;
            ionZ = Stripper_ionZ;
            ChargeStripper(Stripper_ionMass, Stripper_ionProton, ionBeta,
                           Stripper_n, Stripper_ionChargeStates, chargeAmount_Baron);
            // Evaluate change in reference particle energy due to stripper model energy straggling.
            ionEk1 = (tab.Ek[n-2]-StripperPara[2])*Stripper_E0Para[1] + Stripper_E0Para[0];
            ionW = ionEk1 + ionEs;
            iongamma = ionW/ionEs;
            ionBeta = sqrt(1e0-1e0/sqr(iongamma));
            SampleionK = 2e0*M_PI/(ionBeta*SampleLamda);
            tab.s[n-1] = tab.s[n-2];    // Length is zero.
            tab.Ek[n-1] = ionW - ionEs;
            tab.Fy_abs[n-1] = tab.Fy_abs[n-2];
            tab.Beta[n-1] = ionBeta;
            tab.Gamma[n-1] = iongamma;
//            chargeAmount = fribstripper.chargeAmount_Baron;
        }
        std::cout << std::setw(10) << std::left << t_name << std::internal;
        tab.show(std::cout, n-1);
        std::cout << std::endl;
//        if (t_name == "stripper") break;
    }
}


int main(int argc, char *argv[])
{
    CavDataType Data[2];
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

    Data[0].RdData(HomeDir+"/data/axisData_41.txt");
    Data[1].RdData(HomeDir+"/data/axisData_85.txt");

    // Turn trac on/off.
    if (false) sim.set_trace(&std::cout); else sim.set_trace(NULL);

    init_long(sim, Data);

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

