

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
void PrtVec(const std::vector<double> &a)
{
    for (size_t k = 0; k < a.size(); k++)
        std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << a[k];
    std::cout << "\n";
}


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

static
void PrtElement(Moment2ElementBase *elem)
{
    for(unsigned k = 0; k<elem->transfer.size(); k++) {
        std::cout<<"Element "<<elem->index<<" \""<<elem->name<<"\" State "<<k<<" Transfer\n";
        PrtMat(elem->transfer[k]);
    }
}

static
void PrtState(state_t* StatePtr)
{
    unsigned k;

    for (k = 0; k < StatePtr->size(); k++) {
        std::cout << "\nState: "<<k<<" s = "<<std::fixed << std::setprecision(3) << StatePtr->pos << "\n"
                  <<"\n Ref:  "<<StatePtr->ref
                  <<"\n Real: "<<StatePtr->real[k]
                  <<"\n moment0\n";
        PrtVec(StatePtr->moment0[k]);
        std::cout << " moment1\n";
        PrtMat(StatePtr->moment1[k]);
        std::cout << "\n";
    }
}

static
std::vector<double> GetChgState(Config &conf)
{
    return conf.get<std::vector<double> >("IonChargeStates");
}


static
std::vector<double> GetNChg(Config &conf)
{
    return conf.get<std::vector<double> >("NCharge");
}


static
void prt_initial_cond(Machine &sim,
                      state_t &ST)
{
    unsigned     k;

    // Propagate through first element (beam initial conditions).
    sim.propagate(&ST, 0, 1);

    for (k = 0; k < ST.size(); k++) {
        std::cout << "\nIon charge state:\n"
                  << "\nIonZ = " << std::fixed << std::setprecision(5) << std::setw(9) << ST.real[k].IonZ << "\n";

        std::cout << "\nBarycenter:\n";
        PrtVec(ST.moment0[k]);
        std::cout << "\nBeam envelope:\n";
        PrtMat(ST.moment1[k]);

        std::cout << std::fixed << std::setprecision(3) << "\ns [m] = " << ST.pos << "\n";
    }
}


static
void propagate(std::auto_ptr<Config> &conf)
{
    int                        k, nChgStates;
    boost::shared_ptr<Machine> sim;
    std::vector<double>        ChgState;
    std::auto_ptr<state_t>     StatePtr;

    nChgStates = GetChgState(*conf).size();
    ChgState.resize(nChgStates);
    ChgState = GetChgState(*conf);

    for (k = 0; k < nChgStates; k++) {
        conf->set<double>("cstate", k);

        sim = boost::shared_ptr<Machine> (new Machine(*conf));
        sim->set_trace(NULL);

        std::auto_ptr<StateBase> state(sim->allocState());

        state_t* StatePtr(dynamic_cast<state_t*>(state.get()));
        if(!StatePtr) throw std::runtime_error("Only sim_type MomentMatrix2 is supported");

        std::cout << std::fixed << std::setprecision(3) << "\ns [m] = " << StatePtr->pos << "\n";

        // Propagate through first element (beam initial conditions).
        sim->propagate(state.get(), 0, 1);

        clock_t tStamp[2];
        tStamp[0] = clock();
        sim->propagate(state.get(), 1);
        tStamp[1] = clock();

        std::cout << std::fixed << std::setprecision(5)
                  << "\npropagate: " << double(tStamp[1]-tStamp[0])/CLOCKS_PER_SEC << " sec" << "\n";

        PrtState(StatePtr);
    }
}


static
void propagate1(const Config &conf)
{
    // Propagate element-by-element for each charge state.
    unsigned                                   k, elem_no;
    std::vector<double>                        ChgState;

    Machine sim(conf);
    std::auto_ptr<StateBase> state(sim.allocState());
    state_t *StatePtr = dynamic_cast<state_t*>(state.get());
    if(!StatePtr) throw std::runtime_error("Only sim_type MomentMatrix2 is supported");

    prt_initial_cond(sim, *StatePtr);

    clock_t tStamp[2];

    tStamp[0] = clock();

    Machine::iterator it;

    for(it = sim.begin(); it!=sim.end(); ++it)
    {
        ElementVoid* elem   = *it;

        if(strcmp(elem->type_name(), "stripper")==0) {
            elem_no = elem->index;
            std::cout << "\nElement no: " << elem_no << "\n";
            std::cout << elem->conf() << "\n";
            break;
        }

        elem->advance(*state);
        PrtElement(static_cast<Moment2ElementBase*>(elem));

        std::cout<<"After element "<<elem->index<<"\n";
        PrtState(StatePtr);
    }

    {
        Moment2State::vector_t cent, rms;
        GetCenofChg(conf, *StatePtr, cent, rms);

        std::cout << "\nCenter: "<<k<<"\n";
        PrtVec(cent);
        std::cout << "\nRMS: "<<k<<"\n";
        PrtVec(rms);
    }

    Stripper_GetMat(conf, *StatePtr, ChgState);


    for(; it!=sim.end(); ++it) {
        ElementVoid* elem = *it;
        elem->advance(*state);
        std::cout<<"Element "<<elem->index<<" \""<<elem->name<<"\" State "<<k<<" Transfer\n";

        std::cout<<"After element "<<elem->index<<"\n";
        PrtState(StatePtr);
    }

    tStamp[1] = clock();

    PrtState(StatePtr);

    {
        Moment2State::vector_t cent, rms;
        GetCenofChg(conf, *StatePtr, cent, rms);

        std::cout << "\nCenter: "<<k<<"\n";
        PrtVec(cent);
        std::cout << "\nRMS: "<<k<<"\n";
        PrtVec(rms);
    }

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
