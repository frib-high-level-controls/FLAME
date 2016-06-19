#ifndef FLAME_MOMENT_H
#define FLAME_MOMENT_H

#include <ostream>
#include <limits>
#include <iomanip>
#include <math.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "base.h"

#include "constants.h"

// Long. sampling frequency [Hz]; must be set to RF Cavity frequency.
# define SampleFreq   80.5e6
// Sampling distance [m].
# define SampleLambda (C0/SampleFreq*MtoMM)

//! Extra information about a bunch not encoded the vector or matrix of MomentState
//!
struct Particle {
    double IonZ,        //!< Charge state.
           IonQ,        //!< Ion Charge
           IonEs,       //!< Rest energy.
           IonW,        //!< Total energy. (dependent)
           gamma,       //!< Gamma for ion. (dependent)
           beta,        //!< Beta for ion. (dependent)
           bg,          //!< Beta*gamma. (dependent)
           SampleIonK,  //!< Sample rate; different RF Cavity due to RF frequenies. (dependent)
           phis,        //!< Absolute synchrotron phase [rad].
           IonEk;       //!< Kinetic energy.

    Particle() {
        phis = 0.0;
        // initially spoil
        IonZ = IonQ = IonEs = IonW
        = gamma = beta = bg
        = SampleIonK = IonEk
        = std::numeric_limits<double>::quiet_NaN();
    }

    //! Recalculate dependent (cached) values.
    //! Call after changing IonEs or IonEk
    void recalc() {
        IonW       = IonEs + IonEk;
        gamma      = (IonEs != 0e0)? IonW/IonEs : 1e0;
        beta       = sqrt(1e0-1e0/(gamma*gamma));
        bg         = (beta != 0e0)? beta*gamma : 1e0;
        SampleIonK = 2e0*M_PI/(beta*SampleLambda);
    }
};

std::ostream& operator<<(std::ostream&, const Particle&);

inline bool operator==(const Particle& lhs, const Particle& rhs)
{
    // compare only independent variables.
    return lhs.IonEk==rhs.IonEk && lhs.IonEs==rhs.IonEs
            && lhs.IonZ==rhs.IonZ && lhs.IonQ==rhs.IonQ
            && lhs.phis==rhs.phis;
}

inline bool operator!=(const Particle& lhs, const Particle& rhs)
{
    return !(lhs==rhs);
}

/** State for sim_type=MomentMatrix
 *
 * Represents a set of charge states
 *
 * @see @ref simmoment
 */
struct MomentState : public StateBase
{
    enum {maxsize=7};
    enum param_t {
        PS_X, PS_PX, PS_Y, PS_PY, PS_S, PS_PS,
        PS_QQ // ???
    };

    MomentState(const Config& c);
    virtual ~MomentState();

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

    virtual MomentState* clone() const {
        return new MomentState(*this, clone_tag());
    }

    void recalc() {
        ref.recalc();
        for(size_t i=0; i<real.size(); i++) real[i].recalc();
    }

    void calc_rms();

    inline size_t size() const { return real.size(); } //!< # of charge states

protected:
    MomentState(const MomentState& o, clone_tag);
};

/** @brief An Element which propagates the statistical moments of a bunch
 */
struct MomentElementBase : public ElementVoid
{
    typedef MomentState state_t;

    typedef state_t::matrix_t value_t;

    MomentElementBase(const Config& c);
    virtual ~MomentElementBase();

    void get_misalign(const state_t& ST, const Particle& real, value_t& M, value_t& IM) const;

    virtual void advance(StateBase& s);

    //! Return true if previously calculated 'transfer' matricies may be reused
    //! Should compare new input state against values used when 'transfer' was
    //! last computed
    virtual bool check_cache(const state_t& S) const;
    //! Helper to resize our std::vector s to match the # of charge states
    //! in the provided new input state.
    void resize_cache(const state_t& ST);

    //! recalculate 'transfer' taking into consideration the provided input state
    virtual void recompute_matrix(state_t& ST);

    virtual void show(std::ostream& strm, int level) const;

    Particle last_ref_in, last_ref_out;
    std::vector<Particle> last_real_in, last_real_out;
    //! final transfer matricies
    std::vector<value_t> transfer;
    std::vector<value_t> misalign, misalign_inv;

    //! constituents of misalign
    double dx, dy, pitch, yaw, tilt;

    //! If set, check_cache() will always return false
    bool skipcache;

    virtual void assign(const ElementVoid *other) =0;

protected:
    // scratch space to avoid temp. allocation in advance()
    // An Element can't be shared between multiple threads
    state_t::matrix_t scratch;
};

#endif // FLAME_MOMENT_H
