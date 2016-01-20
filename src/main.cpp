
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
    otherwise:
        std::cout << "*** calFindPhaseTable_simplify: undef. cavity type" << std::endl;
    }

    return ionFys - Fy_c - Fy_abs*multip;
}


void RdData(const std::string FileName, int &n, double Data[][2])
{
    int          k;
    std::string  line;
    std::fstream inf;

    std::cout << "RdData:" << std::endl;
    inf.open(FileName.c_str(), std::ifstream::in);
    n = 0;
    while (getline(inf, line).eof()) {
        n++;
        sscanf(line.c_str(), "%lf %lf", &Data[n-1][0], &Data[n-1][1]);
        // Convert from [mm] to [m].
        Data[n-1][0] *= 1e-3;
    }
    inf.close();

    std::cout << std::endl;
    for (k = 0; k < n; k++)
        std::cout << std::scientific << std::setprecision(5)
                  << std::setw(13) << Data[k][0] << std::setw(13) << Data[k][1] << std::endl;
}


void calGapTrace(const double axisData[], const double ionW0,
                 const double ionFy0, double &ionK, const double ionZ,
                 const double ionEs, const double ionLamda, const double E_fac_ad)
{
    // axisData position: mm ; longitudinal electric field: V/m.
    int    axisStep, ii;
    double ionFy = ionFy0,
           ionW  = ionW0,
           dis,
           dz;

//    dis = axisData.get(axisData.size()-1)[0] - axisData.get(0)[0];
//    dz = dis/(axisData.size()-1); //mm

//    axisStep = axisData.size();
    for (ii = 0; ii < axisStep-1; ii++) {
        double ionFylast = ionFy;
        ionFy += ionK*dz;
//        ionW += ionZ*E_fac_ad*(axisData.get(ii)[1]+axisData.get(ii+1)[1])/2/1e6
//                *Math.cos((ionFylast+ionFy)/2)*dz/1000; //MeV
        double iongamma = ionW/ionEs;
        double ionBeta = sqrt(1e0-1e0/sqr(iongamma));
        if ((ionW-ionEs)<0)
        {
            ionW = ionEs;
            ionBeta = 0;
        }
        ionK = 2*M_PI/(ionBeta*ionLamda*1e3); // rad/mm
    }

//    double[] output = {ionW, ionFy};
//    return output;
}


void init_long(Machine &sim)
{
    /* Ek_tab:      reference particle kinetic energy [MeV/u] vs. s [m].
     * Fy_abs_tab:  reference particle phase [rad] vs. s [m].
     * Beta_tab:    reference particle beta vs. s [m].
     * Gamma_tab:   reference particle beta vs. s [m].                               */

    const int max_tab = 2000;

    int    n, cavi;
    double iongamma, ionBeta, SampleLamda, SampleionK, ionFy_i;
    double multip, caviionK, ionFys, E_fac_ad;
    double s_tab[max_tab], Ek_tab[max_tab], Fy_abs_tab[max_tab], Beta_tab[max_tab];
    double Gamma_tab[max_tab], DrivenPhase_tab[max_tab];

    const double c0         = 2.99792458e8,   // Speed of light [m/s].
                 u          = 931.49432e6,    // Atomic mass unit [eV/c^2].
                 mu0        = 4e0*M_PI*1e-7,  // Vacuum permeability.
                 SampleFref = 80.5e6;         // Long. sampling frequency [Hz]; must be set to RF freq.

    // Turn trac on/off.
    if (false)
        sim.set_trace(&std::cout);
    else
        sim.set_trace(NULL);

    Config D;
    std::auto_ptr<StateBase> state(sim.allocState(D));
    // Propagate through first element.
    sim.propagate(state.get(), 0, 1);

    n = 0;
    iongamma    = state->ionW/state->ionEs;
    ionBeta     = sqrt(1e0-1e0/sqr(iongamma));
    SampleLamda = c0/SampleFref/1e6;
    SampleionK  = 2e0*M_PI/(ionBeta*SampleLamda); // [m^-1].

    std::cout << *state << std::endl;
    std::cout << std::scientific << std::setprecision(5)
              << "ionEs = " << state->ionEs  << ", ionEk = " << state->ionEk
              << ", ionW = " << state->ionW << std::endl << std::endl;

    s_tab[n]      = 0e0;
    Ek_tab[n]     = state->ionEk;
    Fy_abs_tab[n] = 0e0;
    Beta_tab[n]   = ionBeta;
    Gamma_tab[n]  = iongamma;

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
            s_tab[n]      = s_tab[n-1] + conf.get<double>("L");
            Fy_abs_tab[n] = Fy_abs_tab[n-1] + SampleionK*conf.get<double>("L");
            Ek_tab[n]     = Ek_tab[n-1];
            Beta_tab[n]   = Beta_tab[n-1];
            Gamma_tab[n]  = Gamma_tab[n-1];
        } else if (t_name == "sbend") {
        } else if (t_name == "quadrupole") {
        } else if (t_name == "solenoid") {
            n++;
            s_tab[n]      = s_tab[n-1] + conf.get<double>("L");
            Fy_abs_tab[n] = Fy_abs_tab[n-1] + SampleionK*conf.get<double>("L");
            Ek_tab[n]     = Ek_tab[n-1];
            Beta_tab[n]   = Beta_tab[n-1];
            Gamma_tab[n]  = Gamma_tab[n-1];
        } else if (t_name == "rfcavity") {
            n++;
            if ((conf.get<double>("L") == 0.24) &&
                (conf.get<double>("f") == 80.5e6)) {
                // 0.041QWR.
                cavi = 1;
//                axisData = "axisData_41.txt";
            } else if ((conf.get<double>("L") == 0.30) &&
                       (conf.get<double>("f") == 80.5e6)) {
                // 0.085QWR.
                cavi = 2;
//                axisData = "axisData_85.txt";
            } else {
                std::cout << "*** init_long: undef. cavity type" << std::endl;
                exit(1);
            }
            multip   = conf.get<double>("f")/SampleFref;
            caviionK = 2e0*M_PI/(ionBeta*c0/conf.get<double>("f"));
             // Synchrotron phase.
            ionFys   = conf.get<double>("f")*M_PI/180e0;
            // Electric field scale factor.
            E_fac_ad = conf.get<double>("scl_fac");

            DrivenPhase_tab[n] = calFindPhaseTable_simplify(
                        cavi, state->ionEk, ionFys, Fy_abs_tab[n-1], multip);

//            ionFy_i = multip*Fy_abs + caviFy_tab;
//            // Calculate influence of the cavity of reference particle kinetic energy, absolute phase, beta and gamma
//            // The data is then recorded for later TLM tracking use
//            double[] output = calGapTrace(axisData, ionW, ionFy_i, caviionK, fribPara.ionZ, FRIBPara.ionEs,
//                                        caviLamda, E_fac_ad);
//            // calGapTrace(axisData_41_ad,FRIBPara.Wion,FRIBPara.Fy,FRIBPara.ionK,FRIBPara.Zion,FRIBPara.Esion,FRIBPara.QWR1.lamda);
//            ionW = output[0];
//            double ionFy_o = output[1];
//            ionEk = ionW-FRIBPara.ionEs;
//            entry = new double[]{tlmnode.position+tlmnode.length,ionEk};
//            fribPara.Ek_tab.add(entry);
//            iongamma = ionW/FRIBPara.ionEs;
//            ionBeta = Math.sqrt(1-1/(iongamma*iongamma));
//            SampleionK = 2*Math.PI/(ionBeta*SampleLamda*1e3); // rad/mm
//            Fy_abs = Fy_abs+(ionFy_o-ionFy_i)/multip;
//            entry = new double[]{tlmnode.position+tlmnode.length,Fy_abs};
//            fribPara.Fy_abs_tab.add(entry);
//            entry = new double[]{tlmnode.position+tlmnode.length,ionBeta};
//            fribPara.Beta_tab.add(entry);
//            entry = new double[]{tlmnode.position+tlmnode.length,iongamma};
//            fribPara.Gamma_tab.add(entry);


        }
        std::cout << std::scientific << std::setprecision(5)
                  << std::setw(10) << std::left << t_name << std::internal
                  << std::setw(13) << s_tab[n]
                  << std::setw(13) << Fy_abs_tab[n]
                  << std::setw(13) << Ek_tab[n]
                  << std::setw(13) << Beta_tab[n]
                  << std::setw(13) << Gamma_tab[n]
                  << std::setw(13) << DrivenPhase_tab[n] << std::endl;
    }
}


int main(int argc, char *argv[])
{
    const int MaxNCav = 10;

    int    n;
    double Data[MaxNCav][2];
    Machine::p_elements_t::const_iterator it;

    const std::string HomeDir = "/home/johan/tlm_workspace/TLM_JB";

    FILE *in = stdin;
    if(argc>1) {
        in = fopen(argv[1], "r");
        if(!in) {
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

    std::cout << "RdData:" << std::endl;
    RdData(HomeDir+"axisData_41.txt", n, Data);

    init_long(sim);

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

