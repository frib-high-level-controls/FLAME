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

    unsigned ndim;
    size_t dim[2];

    VIOCObserver *measure;
};

static
long measure_init_common(dbCommon *prec, const char *link)
{
    try {
        // "simname elementname param"  implies inst==0 and ndim==0
        // "simname elementname param[i]"
        // "simname elementname param[i,j]"
        // "simname elementname[inst] param"
        // "simname elementname[inst] param[i]"
        // "simname elementname[inst] param[i,j]"
        // matches between 3 and 6 groups
        static boost::regex linkpat("(\\S+) ([^\\s\\[]+)(?:\\[(\\d+)\\])? ([^\\s\\[]+)(?:\\[(\\d+)(?:,(\\d+))?\\])?");

        boost::cmatch M;
        if(!boost::regex_match(link, M, linkpat))
            throw std::runtime_error("Bad link string");

        std::unique_ptr<SimDevMeasScalar> priv(new SimDevMeasScalar);
        priv->prec = (dbCommon*)prec;

        if(!find(SimGlobal.sims, M.str(1), priv->sim))
            throw std::runtime_error("No such simulation instance");

        size_t inst = 0;
        if(M[3].matched)
            inst = boost::lexical_cast<size_t>(M.str(3));

        priv->ndim = 0;
        if(M[5].matched) {
            priv->ndim++;
            priv->dim[0] = boost::lexical_cast<size_t>(M.str(5));
        }
        if(M[6].matched) {
            priv->ndim++;
            priv->dim[1] = boost::lexical_cast<size_t>(M.str(6));
        }

        ElementVoid* elem = priv->sim->machine->find(M.str(2), inst);
        if(!elem)
            throw std::runtime_error("No such element");
        priv->element_index = elem->index;

        priv->measure = priv->sim->get_measure(elem->index);

        priv->param = M.str(4);

        unsigned idx;
        {
            Config empty;
            std::unique_ptr<StateBase> state(priv->sim->machine->allocState(empty));

            for(idx=0; true; idx++) {
                StateBase::ArrayInfo info;
                if(!state->getArray(idx, info))
                    throw std::runtime_error("state has no parameter");

                if(info.name==priv->param) {
                    if(info.type!=StateBase::ArrayInfo::Double)
                        throw std::runtime_error("Requested parameter must be type Double");
                    else if(info.ndim!=priv->ndim)
                        throw std::runtime_error("Requested parameter has different cardinality");

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

            if(priv->measure->last.get() &&
               priv->measure->last->getArray(priv->param_index, info) &&
               priv->ndim==info.ndim &&
               info.inbounds(priv->dim))
            {
                double * const arr = info.get<double>(priv->dim);
                value = *arr;

            } else {
                valid = false;
            }

        }

        prec->val = value;
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
