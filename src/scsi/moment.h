#ifndef SCSI_MOMENT_H
#define SCSI_MOMENT_H

#include <ostream>
#include <math.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "base.h"

/** @brief Simulation state which include only a matrix
 */
struct MomentState : public StateBase
{
    enum {maxsize=7};
    enum param_t {
        PS_X, PS_PX, PS_Y, PS_PY, PS_S, PS_PS
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

    void assign(const StateBase& other);

    virtual void show(std::ostream& strm) const;

    double IonEs,
           IonEk,
           IonW;

    vector_t moment0; //!< the state vector (0th moment)
    matrix_t state; //!< The state matrix (1st moment)

    virtual bool getArray(unsigned idx, ArrayInfo& Info);

    virtual MomentState* clone() const {
        return new MomentState(*this, clone_tag());
    }

protected:
    MomentState(const MomentState& o, clone_tag);
};

/** @brief An Element which propagates the statistical moments of a bunch
 *
 */
struct MomentElementBase : public ElementVoid
{
    typedef MomentState state_t;

    MomentElementBase(const Config& c);
    virtual ~MomentElementBase();

    virtual void advance(StateBase& s) const;

    virtual void show(std::ostream& strm) const;

    typedef boost::numeric::ublas::matrix<double> value_t;

    value_t transfer; //!< The transfer matrix

    virtual void assign(const ElementVoid *other)
    {
        const MomentElementBase *O = static_cast<const MomentElementBase*>(other);
        transfer = O->transfer;
        ElementVoid::assign(other);
    }

private:
    // scratch space to avoid temp. allocation in advance()
    // An Element can't be shared between multiple threads
    //TODO: non-const advance()
    mutable state_t::matrix_t scratch;
};

#endif // SCSI_MOMENT_H
