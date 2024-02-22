#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <ctime>

#include <boost/numeric/ublas/io.hpp>

#include <flame/constants.h>
#include <flame/base.h>
#include <flame/moment.h>
#include <flame/chg_stripper.h>
#include <flame/state/vector.h>
#include <flame/state/matrix.h>

#include <flame/moment_sup.h>


typedef MomentState state_t;
typedef state_t::vector_t value_vec;
typedef state_t::matrix_t value_mat;

extern int glps_debug;


static
void PrtVec(const MomentState::vector_t &a)
{
    for (size_t k = 0; k < a.size(); k++)
        std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << a[k];
    std::cout << "\n";
}

static
void PrtMat(const value_mat &M)
{
    for (size_t j = 0; j < M.size1(); j++) {
        for (size_t k = 0; k < M.size2(); k++)
            std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << M(j, k);
        std::cout << "\n";
    }
}

void PrtState(const state_t& ST)
{
    MomentState::vector_t CenofChg, BeamRMS;

    if (true)
        for (size_t k = 0; k < ST.size(); k++) {
            std::cout << "\nState: "<<k<<" s = "<<std::fixed << std::setprecision(5) << ST.pos << "\n"
                      <<"\n Ref:  "<<ST.ref
                     <<"\n Real: "<<ST.real[k]
                       <<"\n moment0\n";
            PrtVec(ST.moment0[k]);
            std::cout << " moment1\n";
            PrtMat(ST.moment1[k]);
        }

    std::cout<<"\nCenter:        "<< std::scientific << std::setprecision(10)
             << std::setw(18)<<ST.moment0_env<<"\n";
    std::cout<<"RMS Beam Size: "<< std::scientific << std::setprecision(10)
             << std::setw(18)<<ST.moment0_rms<<"\n";
    PrtMat(ST.moment1_env);
}

static
void prt_initial_cond(Machine &sim,
                      state_t &ST)
{
    sim.propagate(&ST, 0, 1);
    PrtState(ST);
}

// static
// void PrtOut(std::ofstream &outf1, std::ofstream &outf2, std::ofstream &outf3, const state_t& ST)
// {
//     outf1 << std::scientific << std::setprecision(14) << std::setw(22) << ST.pos;
//     for (size_t j = 0; j < ST.size(); j++)
//         for (int k = 0; k < 6; k++)
//             outf1 << std::scientific << std::setprecision(14) << std::setw(22) << ST.moment0[j][k];
//     for (int k = 0; k < 6; k++)
//         outf1 << std::scientific << std::setprecision(14) << std::setw(22) << ST.moment0_env[k];
//     outf1 << "\n";

//     outf2 << std::scientific << std::setprecision(14) << std::setw(22) << ST.pos;
//     for (size_t j = 0; j < ST.size(); j++)
//         for (int k = 0; k < 6; k++)
//             outf2 << std::scientific << std::setprecision(14) << std::setw(22) << sqrt(ST.moment1[j](k, k));
//     for (int k = 0; k < 6; k++)
//         outf2 << std::scientific << std::setprecision(14) << std::setw(22) << sqrt(ST.moment1_env(k, k));
//     outf2 << "\n";

//     outf3 << std::scientific << std::setprecision(14)
//           << std::setw(22) << ST.pos << std::setw(22) << ST.ref.phis << std::setw(22) << ST.ref.IonEk << "\n";
// }

static
void propagate(const Config &conf)
{
    // Propagate element-by-element for each charge state.
    Machine                  sim(conf);
    std::unique_ptr<StateBase> state(sim.allocState());
    state_t                  *StatePtr = dynamic_cast<state_t*>(state.get());
    std::ofstream            outf1, outf2, outf3;

    if(!StatePtr) throw std::runtime_error("Only sim_type MomentMatrix is supported");

    outf1.open("moment0.txt",   std::ios::out);
    outf2.open("moment1.txt",   std::ios::out);
    outf3.open("ref_orbit.txt", std::ios::out);

    prt_initial_cond(sim, *StatePtr);

    clock_t tStamp[2];

    tStamp[0] = clock();

    prt_initial_cond(sim, *StatePtr);

    Machine::iterator it = sim.begin()+1;
    while (it != sim.end()) {
        ElementVoid* elem = *it;

        elem->advance(*state);
        ++it;

//        PrtOut(outf1, outf2, outf3, *StatePtr);

//        PrtState(*StatePtr);
    }

    outf1.close();
    outf2.close();
    outf3.close();

    tStamp[1] = clock();

    PrtState(*StatePtr);

    std::cout << std::fixed << std::setprecision(5)
              << "\npropagate: " << double(tStamp[1]-tStamp[0])/CLOCKS_PER_SEC << " sec" << "\n";
}


int main(int argc, char *argv[])
{
    try {
        std::unique_ptr<Config> conf;

        if(argc>2)

        glps_debug = 0; // 0 or 1.

        try {
            GLPSParser P;
            conf.reset(P.parse_file(argc>1 ? argv[1] : NULL));
            fprintf(stderr, "Parsing succeeds\n");
        } catch(std::exception& e) {
            fprintf(stderr, "Parse error: %s\n", e.what());
            return 1;
        }

//        std::cout<<"# Reduced lattice\n";
//        GLPSPrint(std::cout, *conf);
//        std::cout<<"\n";

        registerMoment();

        propagate(*conf);

        return 0;
    } catch(std::exception& e) {
        std::cerr << "Main exception: " << e.what() << "\n";
        Machine::registeryCleanup();
        return 1;
    }
}
