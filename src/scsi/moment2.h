#ifndef SCSI_MOMENT_H
#define SCSI_MOMENT_H

#include <ostream>
#include <math.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "base.h"

/** @brief Simulation state which include only a matrix
 */
struct Moment2State : public StateBase
{
    enum {maxsize=7};
    enum param_t {
        PS_X, PS_PX, PS_Y, PS_PY, PS_S, PS_PS
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

    void assign(const StateBase& other);

    virtual void show(std::ostream& strm) const;

    double energy;

    bool do_recalc_energy;
    vector_t moment0;
    matrix_t state; // TODO: better name

    virtual bool getArray(unsigned idx, ArrayInfo& Info);

    virtual Moment2State* clone() const {
        return new Moment2State(*this, clone_tag());
    }

protected:
    Moment2State(const Moment2State& o, clone_tag);
};

/** @brief An Element which propagates the statistical moments of a bunch
 *
 */
struct Moment2ElementBase : public ElementVoid
{
    typedef Moment2State state_t;

    Moment2ElementBase(const Config& c);
    virtual ~Moment2ElementBase();

    virtual void advance(StateBase& s);
    virtual void recalc_energy_gain(state_t& s)=0;
    virtual void recalc_transfer(state_t& s)=0;

    virtual void show(std::ostream& strm) const;

    typedef boost::numeric::ublas::matrix<double> value_t;

    bool do_recalc_energy;
    double energy_gain;
    bool do_recalc_transfer;
    value_t transfer;
    //value_t transferT;

    virtual void assign(const ElementVoid *other);

private:
    // scratch space to avoid temp. allocation in advance()
    // An Element can't be shared between multiple threads
    //TODO: non-const advance()
    mutable state_t::matrix_t scratch;
};

#endif // SCSI_MOMENT_H
