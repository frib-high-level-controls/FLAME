#include <sstream>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <aiRecord.h>
#include <waveformRecord.h>
#include <menuFtype.h>

#include "flame.h"

struct SimDevOrbit : public SimDev
{
    std::string param;
    unsigned param_index; //!< index of named parameter in StateBase::getArray
    size_t param_offset; //!< value offset in parmeter array

    Machine::lookup_iterator start,end;
    size_t nelems;

    std::vector<double*> pvalues;
};

static
long orbit_init_common(dbCommon *prec, const char *link)
{
    try {
        // "simname elementname param"  implies off==0
        // "simname elementname param[off]"
        static boost::regex linkpat("(\\S+) ([^\\s\\[]+) ([^\\s\\[]+)(?:\\[(\\d+)\\])?");

        boost::cmatch M;
        if(!boost::regex_match(link, M, linkpat))
            throw std::runtime_error("Bad link string");

        std::auto_ptr<SimDevOrbit> priv(new SimDevOrbit);
        priv->prec = (dbCommon*)prec;

        if(!find(SimGlobal.sims, M.str(1), priv->sim))
            throw std::runtime_error("No such simulation instance");

        priv->param_offset = 0;
        if(M[4].matched)
            priv->param_offset = boost::lexical_cast<size_t>(M.str(4));

        std::pair<Machine::lookup_iterator, Machine::lookup_iterator> range = priv->sim->machine->equal_range_type(M.str(2));

        if(range.first==range.second)
            throw std::runtime_error("No such elements");

        priv->start = range.first;
        priv->end = range.second;

        size_t N=0;
        for(Machine::lookup_iterator it = range.first; it!=range.second; ++it, N++) {
            ElementVoid *elem = *it;
            priv->sim->get_measure(elem->index);
        }
        priv->nelems = N;
        priv->pvalues.resize(N, NULL);
        assert(N>0);

        priv->param = M.str(3);

        unsigned idx;
        {
            Config empty;
            std::auto_ptr<StateBase> state(priv->sim->machine->allocState(empty));

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
long orbit_init_wf(waveformRecord *prec)
{
    assert(prec->ftvl==menuFtypeDOUBLE);
    return orbit_init_common((dbCommon*)prec, prec->inp.value.instio.string);
}

static
long orbit_read_wf(waveformRecord *prec)
{
    TRY(SimDevOrbit) {
        Guard G(priv->sim->lock);
        bool valid = priv->sim->valid;

        double *pbuf = (double*)prec->bptr;

        size_t i = 0;

        for(Machine::lookup_iterator it = priv->start, end = priv->end;
            it!=end; ++it, i++)
        {
            if(i==prec->nelm) {
                (void)recGblSetSevr(prec, READ_ALARM, INVALID_ALARM);
                break;
            }
            if(priv->pvalues[i]) {
                // use cached value of pvalues[]
                pbuf[i] = *priv->pvalues[i];
                continue;
            }

            ElementVoid *elem = *it;
            VIOCObserver *meas = priv->sim->measures[elem->index];
            StateBase::ArrayInfo info;

            if(meas->last.get() &&
               meas->last->getArray(priv->param_index, info))
            {
                // fill cached pvalues[]
                double * const arr = (double*)info.ptr;
                priv->pvalues[i] = &arr[priv->param_offset];
                pbuf[i] = *priv->pvalues[i];

            } else {
                valid = false;
            }
        }

        prec->nord = i;

        if(!valid)
            (void)recGblSetSevr(prec, READ_ALARM, INVALID_ALARM);

        return 0;
    }CATCH_ALARM()
}

DSET6(waveform, SOrbit, orbit_init_wf, Sim::io_aftersim, orbit_read_wf);

static
long orbit_index_wf(waveformRecord *prec)
{
    TRY(SimDevOrbit) {
        Guard G(priv->sim->lock);

        double *pbuf = (double*)prec->bptr;
        size_t i = 0;

        for(Machine::lookup_iterator it = priv->start, end = priv->end;
            it!=end; ++it, i++)
        {
            if(i==prec->nelm) {
                (void)recGblSetSevr(prec, READ_ALARM, INVALID_ALARM);
                break;
            }
            ElementVoid *elem = *it;
            pbuf[i] = double(elem->index);
        }

        prec->nord = i;

        return 0;
    }CATCH_ALARM()
}

DSET6(waveform, IOrbit, orbit_init_wf, Sim::io_aftersim, orbit_index_wf);
