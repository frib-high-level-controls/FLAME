#ifndef FLAME_H
#define FLAME_H

#include <stdarg.h>

#include <epicsMutex.h>
#include <epicsEvent.h>
#include <epicsThread.h>
#include <epicsGuard.h>
#include <epicsTime.h>
#include <epicsExport.h>
#include <errlog.h>

#include <dbScan.h>
#include <devSup.h>
#include <dbCommon.h>

#include <devSup.h>
#include <recGbl.h>
#include <alarm.h>

#include "scsi/base.h"

typedef epicsGuard<epicsMutex> Guard;
typedef epicsGuardRelease<epicsMutex> UnGuard;

struct Sim;

struct SimGlobal_t {
    epicsMutex lock;
    typedef std::map<std::string, Sim*> sims_t;
    sims_t sims;
};
extern SimGlobal_t SimGlobal;

struct VIOCObserver : public Observer
{
    VIOCObserver() {}
    virtual ~VIOCObserver() {}

    std::auto_ptr<StateBase> last;

    virtual void view(const ElementVoid* elem, const StateBase* state);
};

struct VMeasure
{
    Sim *sim;
    size_t element_index;
    std::vector<VIOCObserver*> observers;

    std::auto_ptr<StateBase> reduced;

    void reduce();

    enum {
        SUnknown, SVector, SMoment,
    } stype;
};

struct Sim : public epicsThreadRunable {
    Sim(const std::string& n);
    virtual ~Sim();

    const std::string name;

    bool stop;
    bool doSim;
    bool valid;
    unsigned level;

    epicsMutex lock;
    epicsEvent event;
    epicsThread worker;

    IOSCANPVT aftersim;
    epicsTimeStamp last_run;
    double last_duration;
    std::string last_msg;

    std::vector<Machine*> machines;
    std::vector<double> ncharge; //!< charge of each charge state, used in weighted average of measurements
    double total_charge;

    typedef std::map<size_t, VMeasure*> measures_t;
    measures_t measures;

    VMeasure *get_measure(size_t index);

    virtual void run();

    void queueSim();

    bool test_debug(unsigned level);

    static long io_aftersim(int, dbCommon* prec, IOSCANPVT* io);
};

struct SimDev {
    Sim *sim;
    dbCommon *prec;

    size_t element_index;

    SimDev() :sim(NULL), prec(NULL), element_index(0) {}

    bool test_debug(unsigned level);
};

struct SimDevSetting : public SimDev
{
    std::string param;

    SimDevSetting() {}
};

#define TRY(TYPE) if(!prec->dpvt) return 0; TYPE *priv = (TYPE*)prec->dpvt; try

#define CATCH_ALARM() catch(std::exception& e) { errlogPrintf("%s: error %s\n", prec->name, e.what()); (void)recGblSetSevr(prec, READ_ALARM, INVALID_ALARM); return -1; }

template<typename R>
struct dset6 {
    long num;
    long (*report)(int);
    long (*init)(int);
    long (*init_record)(R*);
    long (*get_iointr_info)(int, dbCommon* prec, IOSCANPVT* io);
    long (*readwrite)(R*);
    long (*junk)(R*);
};
// defines dset devFLAMERECNAME
// eg. DSET6(ao,Foo, ...) -> devFLAMEaoFoo
#define DSET6(REC, NAME, IR, IO, RW) static dset6<REC##Record> devFLAME##REC##NAME = {6, NULL, NULL, IR, IO, RW, NULL}; \
    epicsExportAddress(dset, devFLAME##REC##NAME)

#endif // FLAME_H
