

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


typedef Moment2State state_t;

typedef boost::numeric::ublas::vector<double> value_vec;
typedef boost::numeric::ublas::matrix<double> value_mat;

extern int glps_debug;


void PrtVec(const std::vector<double> &a)
{
    for (size_t k = 0; k < a.size(); k++)
        std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << a[k];
    std::cout << "\n";
}


void PrtVec(const Moment2State::vector_t &a)
{
    for (size_t k = 0; k < a.size(); k++)
        std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << a[k];
    std::cout << "\n";
}


void PrtMat(const value_mat &M)
{
    for (size_t j = 0; j < M.size1(); j++) {
        for (size_t k = 0; k < M.size2(); k++)
            std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << M(j, k);
        std::cout << "\n";
    }
}


void PrtState(std::vector<state_t*> StatePtr)
{
    int k;

    for (k = 0; k < StatePtr.size(); k++) {
        std::cout << "\n";
        PrtVec(StatePtr[k]->moment0);
        std::cout << "\n";
        PrtMat(StatePtr[k]->state);
        std::cout << std::fixed << std::setprecision(3) << "\ns [m] = " << StatePtr[k]->pos << "\n";
    }
}


std::vector<double> GetChgState(Config &conf)
{
    return conf.get<std::vector<double> >("IonChargeStates");
}


std::vector<double> GetNChg(Config &conf)
{
    return conf.get<std::vector<double> >("NCharge");
}


void prt_initial_cond(std::vector<boost::shared_ptr<Machine> > &sim, std::vector<double> ChgState,
                      std::vector<boost::shared_ptr<StateBase> > &ST)
{
    int     k;
    state_t *StatePtr;

    for (k = 0; k < ChgState.size(); k++) {
        std::cout << "\nIon charge state:\n"
                  << "\nIonZ = " << std::fixed << std::setprecision(5) << std::setw(9) << ChgState[k] << "\n";
        // Propagate through first element (beam initial conditions).
        sim[k]->propagate(ST[k].get(), 0, 1);

        StatePtr = dynamic_cast<state_t*>(ST[k].get());

        std::cout << "\nBarycenter:\n";
        PrtVec(StatePtr->moment0);
        std::cout << "\nBeam envelope:\n";
        PrtMat(StatePtr->state);

        std::cout << std::fixed << std::setprecision(3) << "\ns [m] = " << StatePtr->pos << "\n";
    }
}


void propagate(std::auto_ptr<Config> &conf)
{
    int                        k, nChgStates;
    boost::shared_ptr<Machine> sim;
    std::vector<double>        ChgState;
    std::vector<state_t*>      StatePtr;

    nChgStates = GetChgState(*conf).size();
    ChgState.resize(nChgStates);
    ChgState = GetChgState(*conf);

    for (k = 0; k < nChgStates; k++) {
        conf->set<double>("cstate", k);

        sim = boost::shared_ptr<Machine> (new Machine(*conf));
        sim->set_trace(NULL);

        std::auto_ptr<StateBase> state(sim->allocState());

        StatePtr.push_back(dynamic_cast<state_t*>(state.get()));
        if(!StatePtr[k]) throw std::runtime_error("Only sim_type MomentMatrix2 is supported");

        std::cout << std::fixed << std::setprecision(3) << "\ns [m] = " << StatePtr[k]->pos << "\n";

        // Propagate through first element (beam initial conditions).
        sim->propagate(state.get(), 0, 1);

        clock_t tStamp[2];
        tStamp[0] = clock();
        sim->propagate(state.get(), 1);
        tStamp[1] = clock();

        std::cout << std::fixed << std::setprecision(5)
                  << "\npropagate: " << double(tStamp[1]-tStamp[0])/CLOCKS_PER_SEC << " sec" << "\n";
    }

    PrtState(StatePtr);
}


void propagate1(std::auto_ptr<Config> &conf)
{
    // Propagate element-by-element for each charge state.
    int                                        k, nChgStates, elem_no;
    std::vector<double>                        ChgState;
    std::vector<boost::shared_ptr<Machine> >   sim;
    std::vector<boost::shared_ptr<StateBase> > state;
    std::vector<state_t*>                      StatePtr;
    std::vector<Machine::iterator>             it;

    nChgStates = GetChgState(*conf).size();
    ChgState.resize(nChgStates);
    ChgState = GetChgState(*conf);

    for (k = 0; k < nChgStates; k++) {
        conf->set<double>("cstate", k);

        sim.push_back(boost::shared_ptr<Machine> (new Machine(*conf)));
        sim[k]->set_trace(NULL);

        state.push_back(boost::shared_ptr<StateBase> (sim[k]->allocState()));

        StatePtr.push_back(dynamic_cast<state_t*>(state[k].get()));
        if(!StatePtr[k]) throw std::runtime_error("Only sim_type MomentMatrix2 is supported");

        it.push_back(sim[k]->begin());
    }

    prt_initial_cond(sim, ChgState, state);

    clock_t tStamp[2];

    tStamp[0] = clock();

    while (it[0] != sim[0]->end()) {
        ElementVoid* elem   = *it[0];
        std::string  t_name = elem->type_name(); // C string -> C++ string.

        if (t_name == "stripper") {
            elem_no = it[0] - sim[0]->begin();
            std::cout << "\nElement no: " << elem_no << "\n";
            std::cout << conf->get<Config::vector_t>("elements")[elem_no] << "\n";
            Stripper_GetMat(*conf, sim, state, ChgState);

            StatePtr.clear();
            it.clear();
            for (k = 0; k < state.size(); k++) {
                StatePtr.push_back(dynamic_cast<state_t*>(state[k].get()));
                if(!StatePtr[k]) throw std::runtime_error("Only sim_type MomentMatrix2 is supported");

                it.push_back(sim[k]->begin()+elem_no+1);
            }

            while (it[0] != sim[0]->end()) {
                for (k = 0; k < state.size(); k++) {
                    (*it[k])->advance(*state[k]);
                    ++it[k];
                }
            }

            break;
        }

        for (k = 0; k < nChgStates; k++) {
            (*it[k])->advance(*state[k]);
            ++it[k];
        }
    }

    tStamp[1] = clock();

    PrtState(StatePtr);

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
        propagate1(conf);

        return 0;
    } catch(std::exception& e) {
        std::cerr << "Main exception: " << e.what() << "\n";
        Machine::registeryCleanup();
        return 1;
    }
}
