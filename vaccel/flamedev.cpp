#include <typeinfo>

#include <waveformRecord.h>
#include <boRecord.h>
#include <aiRecord.h>
#include <longinRecord.h>

#include <menuFtype.h>

#include "flame.h"
#include "flame/state/vector.h"

long Sim::io_aftersim(int, dbCommon* prec, IOSCANPVT* io)
{
    if(!prec->dpvt) return 0;
    SimDev *priv = (SimDev*)prec->dpvt;
    try {
        *io = priv->sim->aftersim;
        return 0;
    }catch(std::exception& e){
        errlogPrintf("%s: error in %s: %s\n", prec->name, __FUNCTION__, e.what());
        return -1;
    }
}

VIOCObserver* Sim::get_measure(size_t index)
{
    VIOCObserver *ret;
    if(!find(measures, index, ret))
    {
        ElementVoid *elem = machine->at(index);
        assert(elem->observer()==NULL);
        ret = new VIOCObserver;
        elem->set_observer(ret);
        measures[index] = ret;
    }
    return ret;
}

bool SimDev::test_debug(unsigned level)
{
    return level<=sim->level || level<prec->tpro;
}

static
void init_global(dbCommon *prec, const char *link)
{
    try {
        std::unique_ptr<SimDev> priv(new SimDev);
        priv->prec = prec;

        if(!find(SimGlobal.sims, link, priv->sim))
            throw std::runtime_error("No such sim");

        prec->dpvt = priv.release();
    }catch(std::exception& e){
        fprintf(stderr, "%s: init error: %s\n", prec->name, e.what());
    }
}

static
long init_global_wf(waveformRecord *prec)
{
    assert(prec->inp.type==INST_IO);
    if(prec->ftvl!=menuFtypeCHAR) return 0;
    init_global((dbCommon*)prec, prec->inp.value.instio.string);
    return 0;
}

static
long read_last_msg(waveformRecord *prec)
{
    TRY(SimDev) {
        Guard G(priv->sim->lock);
        size_t tocopy = std::min(size_t(prec->nelm), priv->sim->last_msg.size());
        memcpy(prec->bptr, priv->sim->last_msg.c_str(), tocopy);
        prec->nord = tocopy;
        return 0;
    }CATCH_ALARM()
}

DSET6(waveform, LastMsg, init_global_wf, Sim::io_aftersim, read_last_msg);

static
long init_global_bo(boRecord *prec)
{
    assert(prec->out.type==INST_IO);
    init_global((dbCommon*)prec, prec->out.value.instio.string);
    return 0;
}

static
long write_do_sim(boRecord *prec)
{
    TRY(SimDev) {
        Guard G(priv->sim->lock);
        priv->sim->queueSim();

        return 0;
    }CATCH_ALARM()
}

DSET6(bo, DoSim, init_global_bo, NULL, write_do_sim);

static
long init_global_ai(aiRecord *prec)
{
    assert(prec->inp.type==INST_IO);
    init_global((dbCommon*)prec, prec->inp.value.instio.string);
    return 0;
}

static
long read_global_runtime(aiRecord *prec)
{
    TRY(SimDev) {
        Guard G(priv->sim->lock);

        prec->val = priv->sim->last_duration;
        prec->udf = 0;

        if(prec->aslo!=0.0) prec->val *= prec->aslo;
        prec->val += prec->aoff;

        if(prec->linr && prec->eslo!=0.0) prec->val = prec->val*prec->eslo + prec->eoff;

        if(prec->tse==epicsTimeEventDeviceTime)
            prec->time = priv->sim->last_run;

        return 2;
    }CATCH_ALARM()
}

DSET6(ai, Runtime, init_global_ai, Sim::io_aftersim, read_global_runtime);

static
long init_global_longin(longinRecord *prec)
{
    assert(prec->inp.type==INST_IO);
    init_global((dbCommon*)prec, prec->inp.value.instio.string);
    return 0;
}

static
long read_count(longinRecord *prec)
{
    prec->val += 1;
    return 0;
}

DSET6(longin, Runcount, init_global_longin, Sim::io_aftersim, read_count);
