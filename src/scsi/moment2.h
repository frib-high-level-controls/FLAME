#ifndef SCSI_MOMENT_H
#define SCSI_MOMENT_H

#include <ostream>
#include <limits>
#include <iomanip>
#include <math.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "base.h"


// Long. sampling frequency [Hz]; must be set to RF Cavity frequency.
# define SampleFreq   80.5e6
// Sampling distance [m].
# define SampleLambda (C0/SampleFreq*MtoMM)


/** @brief Simulation state which include only a matrix
 */

struct Particle {
    double IonZ,        // Charge state.
           IonQ,        // Ion Charge
           IonEs,       // Rest energy.
           IonW,        // Total energy.
           gamma,       // Gamma for ion.
           beta,        // Beta for ion.
           bg,          // Beta*gamma;
           SampleIonK,
           phis,        // Synchrotron phase.
           IonEk;       // Kinetic energy.

    Particle() {
        phis = 0.0;
        // initially spoil
        IonZ = IonQ = IonEs = IonW
        = gamma = beta = bg
        = SampleIonK = IonEk
        = std::numeric_limits<double>::quiet_NaN();
    }

    // call after changing IonEs or IonEk
    void recalc() {
        IonW       = IonEs + IonEk;
        gamma      = (IonEs != 0e0)? IonW/IonEs : 1e0;
        beta       = sqrt(1e0-1e0/sqr(gamma));
        bg         = (beta != 0e0)? beta*gamma : 1e0;
        SampleIonK = 2e0*M_PI/(beta*SampleLambda);
    }
};

std::ostream& operator<<(std::ostream&, const Particle&);

struct Moment2State : public StateBase
{
    enum {maxsize=7};
    enum param_t {
        PS_X, PS_PX, PS_Y, PS_PY, PS_S, PS_PS,
        PS_QQ // ???
    };

    Moment2State(const Config& c);
    virtual ~Moment2State();

    typedef boost::numeric::ublas::vector<double,
                    boost::numeric::ublas::bounded_array<double, maxsize>
    > vector_t;

    typedef boost::numeric::ublas::matrix<double,
                    boost::numeric::ublas::row_major,
                    boost::numeric::ublas::bounded_array<double, maxsize*maxsize>
    > matrix_t;

    virtual void assign(const StateBase& other);

    virtual void show(std::ostream& strm, int level) const;

    Particle ref;

    // all three must have the same length, which is the # of charge states
    std::vector<Particle> real;
    std::vector<vector_t> moment0;
    std::vector<matrix_t> moment1;

    vector_t moment0_env, moment0_rms;
    matrix_t moment1_env;

    virtual bool getArray(unsigned idx, ArrayInfo& Info);

    virtual Moment2State* clone() const {
        return new Moment2State(*this, clone_tag());
    }

    void recalc() {
        for(size_t i=0; i<real.size(); i++) real[i].recalc();
    }

    void calc_rms();

    inline size_t size() const { return real.size(); } // # of charge states

protected:
    Moment2State(const Moment2State& o, clone_tag);
};

/** @brief An Element which propagates the statistical moments of a bunch
 *transfer_raw
 */
struct Moment2ElementBase : public ElementVoid
{
    typedef Moment2State state_t;

    Moment2ElementBase(const Config& c);
    virtual ~Moment2ElementBase();

    virtual void advance(StateBase& s);

    virtual void recompute_matrix(state_t&);

    virtual void show(std::ostream& strm, int level) const;

    typedef boost::numeric::ublas::matrix<double,
                    boost::numeric::ublas::row_major,
                    boost::numeric::ublas::bounded_array<double, state_t::maxsize*state_t::maxsize>
    > value_t;

    std::vector<double> last_Kenergy_in, last_Kenergy_out;
    //! final transfer matrix
    std::vector<value_t> transfer, transfer_raw;
    value_t misalign, misalign_inv;

    virtual void assign(const ElementVoid *other);

protected:
    // scratch space to avoid temp. allocation in advance()
    // An Element can't be shared between multiple threads
    state_t::matrix_t scratch;
};

#endif // SCSI_MOMENT_H
