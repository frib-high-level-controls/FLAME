#include <typeinfo>

#include <waveformRecord.h>
#include <boRecord.h>
#include <aiRecord.h>

#include <menuFtype.h>

#include "flame.h"
#include "scsi/state/vector.h"
#include "scsi/moment.h"

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

VMeasure* Sim::get_measure(size_t index)
{
    VMeasure *ret;
    if(!find(measures, index, ret))
    {
        ret = new VMeasure;
        ret->sim = this;
        ret->element_index = index;
        ret->observers.resize(machines.size());
        for(size_t i=0; i<machines.size(); i++)
        {
            ElementVoid *elem = machines[i]->at(index);
            assert(elem->observer()==NULL);
            ret->observers[i] = new VIOCObserver;
            elem->set_observer(ret->observers[i]);
        }
        measures[index] = ret;
    }
    return ret;
}

void VMeasure::reduce()
{
    // we break the StateBase abstraction here for reasons of performance, and
    // my not being able to think of any abstract interface for this operation...

    if(!reduced.get())
    {
        StateBase *S = observers[0]->last.get();
        // first run
        if(!S)
            return; // Nothing observed yet???

        if(dynamic_cast<VectorState*>(S)) {
            stype = SVector;
        } else if(dynamic_cast<MomentState*>(S)) {
            stype = SMoment;
        } else {
            stype = SUnknown;
            if(observers.size()!=1)
                errlogPrintf("%s: VMeassure Can't reduce unsupported state type: %s\n",
                             sim->name.c_str(), typeid(S).name());
        }

        reduced.reset(observers[0]->last->clone());

    } else {
        // not first run
        reduced->assign(*observers[0]->last);
    }

    // check for trival case, and unknown StateBase sub-classes
    if(observers.size()==1 || stype==SUnknown)
        return;

    if(stype==SVector) {
        VectorState *ST = static_cast<VectorState*>(reduced.get());

        ST->state *= sim->ncharge[0];

        for(size_t i=1; i<observers.size(); i++)
        {
            VectorState *O = static_cast<VectorState*>(observers[i]->last.get());
            ST->state += sim->ncharge[i] * O->state;
        }

        ST->state /= sim->total_charge;

    } else if(stype==SMoment) {
        MomentState *ST = static_cast<MomentState*>(reduced.get());

        ST->moment0 *= sim->ncharge[0];

        for(size_t i=1; i<observers.size(); i++)
        {
            MomentState *O = static_cast<MomentState*>(observers[i]->last.get());
            ST->moment0 += sim->ncharge[i] * O->moment0;
        }

        ST->moment0 /= sim->total_charge;

        //TODO: reduce ST->state?

    } else
        throw std::logic_error("oops in VMeasure::reduce");
}

bool SimDev::test_debug(unsigned level)
{
    return level<=sim->level || level<prec->tpro;
}

static
void init_global(dbCommon *prec, const char *link)
{
    try {
        std::auto_ptr<SimDev> priv(new SimDev);
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
