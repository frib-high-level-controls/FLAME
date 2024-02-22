#include <sstream>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <aoRecord.h>
#include <aiRecord.h>

#include "flame.h"

namespace {

template<typename REC> // REC is aoRecord or aiRecord
long setting_init_a(REC *prec, DBLINK *plink)
{
    try {
        // "simname elementname param"  implies inst==0
        // "simname elementname[inst] param"
        static boost::regex linkpat("(\\S+) ([^\\s\\[]+)(?:\\[(\\d+)\\])? (\\S+)");

        boost::cmatch M;
        if(!boost::regex_match(plink->value.instio.string, M, linkpat))
            throw std::runtime_error("Bad link string");

        std::unique_ptr<SimDevSetting> priv(new SimDevSetting);
        priv->prec = (dbCommon*)prec;

        if(!find(SimGlobal.sims, M.str(1), priv->sim))
            throw std::runtime_error("No such simulation instance");

        size_t inst = 0;
        if(M[3].matched)
            inst = boost::lexical_cast<size_t>(M.str(3));

        ElementVoid* elem = priv->sim->machine->find(M.str(2), inst);
        if(!elem)
            throw std::runtime_error("No such element");
        priv->element_index = elem->index;

        priv->param = M.str(4);

        // fetch current value to initialize VAL
        // autosave may overwrite
        double initval = elem->conf().get<double>(priv->param);

        if (prec->aslo) initval *= prec->aslo;
        initval += prec->aoff;

        if(prec->linr && prec->eslo!=0) initval = initval*prec->eslo + prec->eoff;

        prec->val = initval;
        prec->udf = 0;
        // Some minor mischief
        // we don't want the initial UDF alarm, as VAL is indeed defined
        // however, processing all settings immediately (PINI=YES)
        // is redundant, and slows down initialization unnecessarily.
        // so just clear it now.
        prec->sevr = prec->stat = 0;

        prec->dpvt = priv.release();
        return 2;
    }catch(std::exception& e){
        fprintf(stderr, "%s: init error: %s\n", prec->name, e.what());
        return -1;
    }
}

long setting_init_ao(aoRecord *prec)
{
    return setting_init_a(prec, &prec->out);
}

long setting_init_ai(aiRecord *prec)
{
    long ret = setting_init_a(prec, &prec->inp);
    if(ret==2) ret=0;
    return ret;
}

static
long setting_change_ao(aoRecord *prec)
{
    TRY(SimDevSetting) {
        Guard G(priv->sim->lock);

        double value = prec->val;
        if(prec->linr) {
            if(prec->eslo!=0) value = (value - prec->eoff) / prec->eslo;
        }
        value -= prec->aoff;
        if (prec->aslo != 0) value /= prec->aslo;

        if(priv->test_debug(2))
            errlogPrintf("%s: set %u %s.%s = %f\n",
                        prec->name,
                        (unsigned)priv->element_index,
                        priv->sim->machine->at(priv->element_index)->name.c_str(),
                        priv->param.c_str(),
                        value);

        {
            ElementVoid* elem = priv->sim->machine->at(priv->element_index);
            Config conf(elem->conf());
            conf.set<double>(priv->param, value);

            priv->sim->machine->reconfigure(priv->element_index, conf);
        }

        return 0;
    }CATCH_ALARM()
}

static
long setting_readback_ai(aiRecord *prec)
{
    TRY(SimDevSetting) {
        Guard G(priv->sim->lock);

        ElementVoid* elem = priv->sim->machine->at(priv->element_index);

        double curval = elem->conf().get<double>(priv->param);

        if (prec->aslo) curval *= prec->aslo;
        curval += prec->aoff;

        if(prec->linr && prec->eslo!=0) curval = curval*prec->eslo + prec->eoff;

        prec->val = curval;

        return 2;
    }CATCH_ALARM()
}

} // namespace

DSET6(ao, Setting, setting_init_ao, NULL, setting_change_ao);
DSET6(ai, Setting, setting_init_ai, Sim::io_aftersim, setting_readback_ai);
