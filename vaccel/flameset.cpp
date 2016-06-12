#include <sstream>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <aoRecord.h>

#include "flame.h"

static
long setting_init_ao(aoRecord *prec)
{
    try {
        // "simname elementname param"  implies inst==0
        // "simname elementname[inst] param"
        static boost::regex linkpat("(\\S+) ([^\\s\\[]+)(?:\\[(\\d+)\\])? (\\S+)");

        boost::cmatch M;
        if(!boost::regex_match(prec->out.value.instio.string, M, linkpat))
            throw std::runtime_error("Bad link string");

        std::auto_ptr<SimDevSetting> priv(new SimDevSetting);
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
//        if(!elem->conf().has<double>(priv->param))
//            throw std::runtime_error("No such parameter for this element");

        prec->dpvt = priv.release();
        return 2;
    }catch(std::exception& e){
        fprintf(stderr, "%s: init error: %s\n", prec->name, e.what());
        return -1;
    }
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

DSET6(ao, Setting, setting_init_ao, NULL, setting_change_ao);
