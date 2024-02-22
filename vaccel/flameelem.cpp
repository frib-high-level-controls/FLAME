#include <sstream>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <aiRecord.h>
#include <longinRecord.h>

#include "flame.h"

static
void elem_init_common(dbCommon *prec, const char *link)
{
    try {
        // "simname elementname"  implies inst==0
        // "simname elementname[inst]"
        static boost::regex linkpat("(\\S+) ([^\\s\\[]+)(?:\\[(\\d+)\\])?");

        boost::cmatch M;
        if(!boost::regex_match(link, M, linkpat))
            throw std::runtime_error("Bad link string");

        std::unique_ptr<SimDev> priv(new SimDev);
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

        prec->dpvt = priv.release();
    }catch(std::exception& e){
        fprintf(stderr, "%s: init error: %s\n", prec->name, e.what());
    }
}

static
long elem_init_li(longinRecord *prec)
{
    elem_init_common((dbCommon*)prec, prec->inp.value.instio.string);
    return 0;
}

static
long elem_read_index(longinRecord *prec)
{
    TRY(SimDev) {
        Guard G(priv->sim->lock);

        prec->val = priv->element_index;

        return 0;
    }CATCH_ALARM()
}

DSET6(longin, Index, elem_init_li, NULL, elem_read_index);
