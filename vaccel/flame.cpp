#include <stdio.h>
#include <iostream>

#include <epicsExport.h>
#include <epicsExit.h>

#include "iocshelper.h"
#include "flame.h"

SimGlobal_t SimGlobal;

Sim::Sim(const std::string& n)
    :name(n)
    ,stop(false)
    ,doSim(false)
    ,valid(false)
    ,level(0)
    ,worker(*this, "flame",
            epicsThreadGetStackSize(epicsThreadStackSmall),
            epicsThreadPriorityMedium)
{
    scanIoInit(&aftersim);

    memset(&last_run, 0, sizeof(last_run));
    last_duration = -1;

    worker.start();
}

Sim::~Sim() {
    {
        Guard G(lock);
        stop = true;
    }
    event.signal();
    worker.exitWait();

    for(size_t i=0; i<machines.size(); i++) {
        delete machines[i];
    }
    machines.clear();
}

void Sim::run()
{
    Guard G(lock);
    while(!stop) {
        if(!doSim) {
            UnGuard U(G);
            event.wait();
            continue;
        }

        doSim = false;

        epicsTimeStamp start;
        epicsTimeGetCurrent(&start);

        try {
            Config empty;

            for(size_t i=0; i<machines.size(); i++) {
                std::auto_ptr<StateBase> state(machines[i]->allocState(empty));

                machines[i]->propagate(state.get());
            }

            for(measures_t::const_iterator it = measures.begin(), end=measures.end(); it!=end; ++it)
            {
                it->second->reduce();
            }

            valid = true;
        }catch(std::exception& e){
            last_msg = e.what();
            valid = false;
        }

        epicsTimeGetCurrent(&last_run);
        last_duration = epicsTimeDiffInSeconds(&last_run, &start);

        scanIoRequest(aftersim);
    }
}

void Sim::queueSim()
{
    if(!doSim) {
        if(test_debug(1))
            errlogPrintf("%s: do sim\n", name.c_str());
        doSim = true;
        event.signal();

    } else {
        if(test_debug(1))
            errlogPrintf("%s: sim already pending\n", name.c_str());
    }
}

bool Sim::test_debug(unsigned lvl)
{
    return lvl<=level;
}

/**
 * @brief Add a new simulation instance
 * @param name unique internal name for this instance
 * @param lattice file name
 */
void flamePrepare(const char *name, const char *lattice)
{
    if(!name) return;
    try {
        if(SimGlobal.sims.find(name)!=SimGlobal.sims.end())
            throw std::runtime_error("Sim name already in use");

        std::auto_ptr<Config> conf;
        {
            GLPSParser parser;
            conf.reset(parser.parse_file(lattice));
        }

        std::auto_ptr<Sim> sim(new Sim(name));

        //const std::string& ST(conf->get<std::string>("sim_type"));

        std::vector<double> cstates;
        if(conf->tryGet<std::vector<double> >("IonChargeStates", cstates)) {
            if(cstates.size()==0)
                throw std::runtime_error("Found empty IonChargeStates[]");
            printf("Found %u change states\n", (unsigned)cstates.size());

            sim->ncharge = conf->get<std::vector<double> >("NCharge");

            if(sim->ncharge.size() != cstates.size())
                throw std::runtime_error("length of IonChargeStates and NCharge must match");

            sim->total_charge = 0.0;
            for(size_t i=0; i<sim->ncharge.size(); i++)
                sim->total_charge += sim->ncharge[i];

            sim->machines.resize(cstates.size(), NULL);

            for(size_t i=0; i<cstates.size(); i++) {
                Config sconf(conf->new_scope());
                sconf.set<double>("ChangeState", double(i));

                sim->machines[i] = new Machine(sconf);
            }

        } else {
            sim->machines.resize(1, NULL);
            sim->ncharge.resize(1, 1.0); // default weight  1
            sim->total_charge = 1.0;

            sim->machines[0] = new Machine(*conf);
        }

        Guard G(SimGlobal.lock);
        SimGlobal.sims[name] = sim.get();
        sim.release();

    }catch(std::exception& e){
        fprintf(stderr, "Error: %s\n", e.what());
    }
}

void flameDebug(const char *name, int level)
{
    if(!name) return;
    try {
        Sim *sim;
        {
            Guard G(SimGlobal.lock);
            if(!find(SimGlobal.sims, name, sim))
                throw std::runtime_error("Unknown sim instance name");
        }
        Guard G(sim->lock);
        sim->level = std::max(0, level);

    }catch(std::exception& e){
        fprintf(stderr, "Error: %s\n", e.what());
    }
}

void flameShowConfig(const char *name)
{
    if(!name) return;
    try {
        Sim *sim;
        {
            Guard G(SimGlobal.lock);
            if(!find(SimGlobal.sims, name, sim))
                throw std::runtime_error("Unknown sim instance name");
        }
        Guard G(sim->lock);

        GLPSPrint(std::cout, sim->machines[0]->conf());

    }catch(std::exception& e){
        fprintf(stderr, "Error: %s\n", e.what());
    }
}

void flameShow(const char *name)
{
    if(!name) return;
    try {
        Sim *sim;
        {
            Guard G(SimGlobal.lock);
            if(!find(SimGlobal.sims, name, sim))
                throw std::runtime_error("Unknown sim instance name");
        }
        Guard G(sim->lock);

        std::cout<<*sim->machines[0]<<"\n";

    }catch(std::exception& e){
        fprintf(stderr, "Error: %s\n", e.what());
    }
}

namespace {

//! Stop worker threads
void flameCleanup(void*)
{
    //Guard G(SimGlobal.lock);

    for(SimGlobal_t::sims_t::iterator it=SimGlobal.sims.begin(), end=SimGlobal.sims.end();
        it!=end; ++it)
    {
        Sim* sim = it->second;

        {
            Guard G(sim->lock);
            sim->stop = true;
        }
        sim->event.signal();
        sim->worker.exitWait();
    }
}

void flameRegistrar()
{
    registerLinear();
    registerMoment();

    iocshRegister<const char*, const char*, &flamePrepare>("flamePrepare", "name", "lattice");
    iocshRegister<const char*, int, &flameDebug>("flameDebug", "name", "level");
    iocshRegister<const char*, &flameShowConfig>("flameShowConfig", "name");
    iocshRegister<const char*, &flameShow>("flameShow", "name");

    epicsAtExit(&flameCleanup, NULL);
}

} // namespace

epicsExportRegistrar(flameRegistrar);
