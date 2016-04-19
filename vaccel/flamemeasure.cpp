#include <sstream>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <aiRecord.h>
#include <waveformRecord.h>
#include <menuFtype.h>

#include "flame.h"

void VIOCObserver::view(const ElementVoid* elem, const StateBase* state)
{
    if(!last.get())
        last.reset(state->clone());
    else
        last->assign(*state);
}

struct SimDevMeasScalar : public SimDev
{
    std::string param;
    unsigned param_index; //!< index of named parameter in StateBase::getArray
    size_t param_offset; //!< value offset in parmeter array

    VMeasure *measure;
    double *pvalue;
};

static
long measure_init_common(dbCommon *prec, const char *link)
{
    try {
        // "simname elementname param"  implies inst==0 and off==0
        // "simname elementname[inst] param[off]"
        static boost::regex linkpat("(\\S+) ([^\\s\\[]+)(?:\\[(\\d+)\\])? ([^\\s\\[]+)(?:\\[(\\d+)\\])?");

        boost::cmatch M;
        if(!boost::regex_match(link, M, linkpat))
            throw std::runtime_error("Bad link string");

        std::auto_ptr<SimDevMeasScalar> priv(new SimDevMeasScalar);
        priv->prec = (dbCommon*)prec;

        if(!find(SimGlobal.sims, M.str(1), priv->sim))
            throw std::runtime_error("No such simulation instance");

        size_t inst = 0;
        if(M[3].matched)
            inst = boost::lexical_cast<size_t>(M.str(3));

        priv->param_offset = 0;
        if(M[5].matched)
            priv->param_offset = boost::lexical_cast<size_t>(M.str(5));

        ElementVoid* elem = priv->sim->machines[0]->find(M.str(2), inst);
        if(!elem)
            throw std::runtime_error("No such element");
        priv->element_index = elem->index;

        priv->measure = priv->sim->get_measure(elem->index);

        priv->param = M.str(4);
        priv->pvalue = NULL;

        unsigned idx;
        {
            Config empty;
            std::auto_ptr<StateBase> state(priv->sim->machines[0]->allocState(empty));

            for(idx=0; true; idx++) {
                StateBase::ArrayInfo info;
                if(!state->getArray(idx, info))
                    throw std::runtime_error("state has no parameter");

                if(info.name==priv->param) {
                    if(info.type!=StateBase::ArrayInfo::Double)
                        throw std::runtime_error("Requested parameter must be type Double");

                    size_t max=0;
                    if(info.ndim)
                        for(size_t dim=0; dim<info.ndim; dim++)
                            max += info.dim[dim];
                    else
                        max = 1;

                    if(priv->param_offset>=max)
                        throw std::runtime_error("Requested parameter offset is out of range");

                    priv->param_index = idx;
                    break;
                }
            }
        }

        prec->dpvt = priv.release();
        return 0;
    }catch(std::exception& e){
        fprintf(stderr, "%s: init error: %s\n", prec->name, e.what());
        return -1;
    }
}

static
long measure_init_ai(aiRecord *prec)
{
    return measure_init_common((dbCommon*)prec, prec->inp.value.instio.string);
}

static
long measure_read_ai(aiRecord *prec)
{
    TRY(SimDevMeasScalar) {
        Guard G(priv->sim->lock);
        bool valid = priv->sim->valid;

        double value = 0.0;

        {
            StateBase::ArrayInfo info;

            if(priv->pvalue) {
                value = *priv->pvalue;

            } else if(priv->measure->reduced.get() &&
                      priv->measure->reduced->getArray(priv->param_index, info))
            {
                double * const arr = (double*)info.ptr;
                priv->pvalue = &arr[priv->param_offset];
                value = *priv->pvalue;

            } else {
                valid = false;
            }

        }

        prec->val = value/priv->sim->total_charge;
        prec->udf = 0;

        if(prec->aslo!=0.0) prec->val *= prec->aslo;
        prec->val += prec->aoff;

        if(prec->linr && prec->eslo!=0.0) prec->val = prec->val*prec->eslo + prec->eoff;

        if(!valid)
            (void)recGblSetSevr(prec, READ_ALARM, INVALID_ALARM);

        return 2;
    } CATCH_ALARM()
}

DSET6(ai, State, measure_init_ai, Sim::io_aftersim, measure_read_ai);
