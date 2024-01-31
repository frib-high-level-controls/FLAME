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
    unsigned ndim;
    size_t dim[2];

    Machine::lookup_iterator start,end;
    size_t nelems;
};

static
long orbit_init_common(dbCommon *prec, const char *link)
{
    try {
        // "simname elementname param"      ndim==0
        // "simname elementname param[i]"   ndim==1
        // "simname elementname param[i,j]" ndim==2
        // matches between 3 and 5 groups
        static boost::regex linkpat("(\\S+) ([^\\s\\[]+) ([^\\s\\[]+)(?:\\[(\\d+)(?:,(\\d+))?\\])?");

        boost::cmatch M;
        if(!boost::regex_match(link, M, linkpat))
            throw std::runtime_error("Bad link string");

        std::unique_ptr<SimDevOrbit> priv(new SimDevOrbit);
        priv->prec = (dbCommon*)prec;

        if(!find(SimGlobal.sims, M.str(1), priv->sim))
            throw std::runtime_error("No such simulation instance");

        priv->ndim = 0;
        if(M[4].matched) {
            priv->ndim++;
            priv->dim[0] = boost::lexical_cast<size_t>(M.str(4));
        }
        if(M[5].matched) {
            priv->ndim++;
            priv->dim[1] = boost::lexical_cast<size_t>(M.str(5));
        }

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
        assert(N>0);

        priv->param = M.str(3);

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
            ElementVoid * const elem = *it;

            if(i==prec->nelm) {
                valid = false;
                break;
            }

            VIOCObserver *meas = priv->sim->measures[elem->index];
            StateBase::ArrayInfo info;

            if(meas->last.get() &&
               meas->last->getArray(priv->param_index, info) &&
               priv->ndim==info.ndim &&
               info.inbounds(priv->dim))
            {
                double * const arr = info.get<double>(priv->dim);
                pbuf[i] = *arr;

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
