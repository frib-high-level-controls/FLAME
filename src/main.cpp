
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <iomanip>

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


void init_long(Machine &sim)
{
    /* Ek_tab:      reference particle kinetic energy [MeV/u] vs. s [m].
     * Fy_abs_tab:  reference particle phase [rad] vs. s [m].
     * Beta_tab:    reference particle beta vs. s [m].
     * Gamma_tab:   reference particle beta vs. s [m].                               */

    const int max_tab = 2000;

    int    n;
    double iongamma, ionBeta, SampleLamda, SampleionK;
    double s_tab[max_tab], Ek_tab[max_tab], Fy_abs_tab[max_tab], Beta_tab[max_tab];
    double Gamma_tab[max_tab];

    const double c0         = 2.99792458e8,   // Speed of light [m/s].
                 u          = 931.49432e6,    // Atomic mass unit [eV/c^2].
                 mu0        = 4e0*M_PI*1e-7,  // Vacuum permeability.
                 SampleFref = 80.5e6;         // Long. sampling frequency [Hz].

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

    std::cout << "ionEs = " << state->ionEs  << ", ionEk = " << state->ionEk
              << ", ionW = " << state->ionW << std::endl;
    std::cout << *state << std::endl;

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
            s_tab[n] = s_tab[n-1] + conf.get<double>("L");
            Fy_abs_tab[n] = Fy_abs_tab[n-1] + SampleionK*conf.get<double>("L");
            Ek_tab[n] = Ek_tab[n-1];
            Beta_tab[n] = Beta_tab[n-1];
            Gamma_tab[n] = Gamma_tab[n-1];
        } else if (t_name == "sbend") {
        } else if (t_name == "quadrupole") {
        } else if (t_name == "solenoid") {
            n++;
            s_tab[n] = s_tab[n-1] + conf.get<double>("L");
            Fy_abs_tab[n] = Fy_abs_tab[n-1] + SampleionK*conf.get<double>("L");
            Ek_tab[n] = Ek_tab[n-1];
            Beta_tab[n] = Beta_tab[n-1];
            Gamma_tab[n] = Gamma_tab[n-1];
        } else if (t_name == "rfcavity") {
//            std::cout << "cavity: L =" << conf.get<double>("L")
//                      << ", phi = " << conf.get<double>("phi") << std::endl;
        }
        std::cout << std::scientific << std::setprecision(5)
                  << std::setw(8) << std::left << t_name << std::internal
                  << std::setw(13) << s_tab[n]
                  << std::setw(13) << Fy_abs_tab[n]
                  << std::setw(13) << Ek_tab[n]
                  << std::setw(13) << Beta_tab[n]
                  << std::setw(13) << Gamma_tab[n] << std::endl;
    }
}


int main(int argc, char *argv[])
{
    Machine::p_elements_t::const_iterator it;

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

