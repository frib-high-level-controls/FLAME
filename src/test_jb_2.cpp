

#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <ctime>

#include <boost/numeric/ublas/io.hpp>

#include <scsi/constants.h>
#include <scsi/base.h>
#include <scsi/moment2.h>
#include <scsi/chg_stripper.h>
#include <scsi/state/vector.h>
#include <scsi/state/matrix.h>

#include <scsi/moment2_sup.h>


typedef Moment2State state_t;
typedef state_t::vector_t value_vec;
typedef state_t::matrix_t value_mat;

extern int glps_debug;

static
void PrtVec(const Moment2State::vector_t &a)
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

void PrtState(const state_t& StatePtr)
{
    size_t k;

    for (k = 0; k < StatePtr.size(); k++) {
        std::cout << "\nState: "<<k<<" s = "<<std::fixed << std::setprecision(3) << StatePtr.pos << "\n"
                          <<"\n Ref:  "<<StatePtr.ref
                          <<"\n Real: "<<StatePtr.real[k]
                          <<"\n moment0\n";
        PrtVec(StatePtr.moment0[k]);
        std::cout << " moment1\n";
        PrtMat(StatePtr.moment1[k]);
        std::cout << "\n";
    }

    std::cout<<"Center: "<< std::scientific << std::setprecision(10)
             << std::setw(18)<<StatePtr.moment0_env<<"\n\n";
}

std::vector<double> GetChgState(Config &conf)
{
    return conf.get<std::vector<double> >("IonChargeStates");
}

static
void prt_initial_cond(Machine &sim,
                      state_t &ST)
{
    sim.propagate(&ST, 0, 1);
    PrtState(ST);
}

static
void propagate1(const Config &conf)
{
    // Propagate element-by-element for each charge state.
    size_t                                     elem_no = (size_t)-1;
    std::vector<double>                        ChgState;

    Machine sim(conf);
    std::auto_ptr<StateBase> state(sim.allocState());
    state_t *StatePtr = dynamic_cast<state_t*>(state.get());
    if(!StatePtr) throw std::runtime_error("Only sim_type MomentMatrix2 is supported");

    prt_initial_cond(sim, *StatePtr);

    clock_t tStamp[2];

    tStamp[0] = clock();

    prt_initial_cond(sim, *StatePtr);

    Machine::iterator it = sim.begin()+1;
    while (it != sim.end()) {
        ElementVoid* elem   = *it;
        Moment2ElementBase* ELEM = static_cast<Moment2ElementBase*>(elem);
        unsigned idx = elem->index;
        std::string  t_name = elem->type_name(); // C string -> C++ string.
        std::cout<<"At element "<<idx<<" "<<elem->name<<"\n";

        if (t_name == "stripper") {
            elem_no = elem->index;
            std::cout << "\nElement no: " << elem_no << "\n";
            std::cout << elem->conf() << "\n";
            //Stripper_GetMat(conf, *StatePtr);

            while (it != sim.end()) {
                elem   = *it;
                ELEM = static_cast<Moment2ElementBase*>(elem);
                unsigned idx = elem->index;
                std::cout<<"At element "<<idx<<" "<<elem->name<<"\n";

                elem->advance(*state);
                for (size_t k = 0; k < ELEM->transfer.size(); k++) {
                    std::cout<<"Transfer "<<k<<"\n";
                    PrtMat(ELEM->transfer[k]);
                }
                ++it;
                std::cout<<"After element "<<idx<<"\n";
                PrtState(*StatePtr);
            }

            break;
        }

        elem->advance(*state);
        for (size_t k = 0; k < ELEM->transfer.size(); k++) {
            std::cout<<"Transfer "<<k<<"\n";
            PrtMat(ELEM->transfer[k]);
        }
        ++it;
        std::cout<<"After element "<<idx<<"\n\n";
        PrtState(*StatePtr);
    }

    tStamp[1] = clock();

    PrtState(*StatePtr);

    std::cout << std::fixed << std::setprecision(5)
              << "\npropagate: " << double(tStamp[1]-tStamp[0])/CLOCKS_PER_SEC << " sec" << "\n";
}


int main(int argc, char *argv[])
{
    try {
        std::auto_ptr<Config> conf;

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

        registerMoment2();

//        propagate(conf);
        propagate1(*conf);

        return 0;
    } catch(std::exception& e) {
        std::cerr << "Main exception: " << e.what() << "\n";
        Machine::registeryCleanup();
        return 1;
    }
}
