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
            std::unique_ptr<StateBase> state(machine->allocState());

            machine->propagate(state.get());

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
        Guard G(SimGlobal.lock);

        if(SimGlobal.sims.find(name)!=SimGlobal.sims.end())
            throw std::runtime_error("Sim name already in use");

        std::unique_ptr<Config> conf;
        {
            GLPSParser parser;
            conf.reset(parser.parse_file(lattice));
        }

        std::unique_ptr<Sim> sim(new Sim(name));

        sim->machine.reset(new Machine(*conf));

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

        GLPSPrint(std::cout, sim->machine->conf());

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

        std::cout<<*sim->machine<<"\n";

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
