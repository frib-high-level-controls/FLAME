#ifndef SCSI_LINEAR_H
#define SCSI_LINEAR_H

#include <ostream>
#include <math.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "base.h"

/** @brief An Element based on a simple Transfer matrix
 *
 * @f$ state = transfer * state @f$
 * @tparam one of VectorState or MatrixState
 */
template<typename State>
struct LinearElementBase : public ElementVoid
{
    typedef State state_t;

    LinearElementBase(const Config& c)
        :ElementVoid(c)
        ,transfer(boost::numeric::ublas::identity_matrix<double>(state_t::maxsize))
    {}
    virtual ~LinearElementBase() {}

    virtual void advance(StateBase& s)
    {
        State& ST = static_cast<State&>(s);
        advanceT(ST);
    }

    virtual void show(std::ostream& strm) const
    {
        ElementVoid::show(strm);
        strm<<"Transfer: "<<transfer<<"\n";
    }

    typedef boost::numeric::ublas::matrix<double> value_t;

    value_t transfer; //!< The transfer matrix

    virtual void assign(const ElementVoid *other)
    {
        const LinearElementBase *O = static_cast<const LinearElementBase*>(other);
        transfer = O->transfer;
        ElementVoid::assign(other);
    }

private:
    void advanceT(State& s)
    {
        using boost::numeric::ublas::prod;
        s.pos += length;
        s.state = prod(transfer, s.state);
    }
};

#endif // SCSI_LINEAR_H
