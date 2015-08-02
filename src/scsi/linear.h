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
 * x' = M * x
 */
template<typename State>
struct LinearElementBase : public ElementVoid
{
    typedef State state_t;

    LinearElementBase(const Config& c)
        :ElementVoid(c)
        ,transfer(boost::numeric::ublas::identity_matrix<double>(6))
    {}
    virtual ~LinearElementBase() {}

    virtual void advance(StateBase& s) const
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

    value_t transfer;
private:
    void advanceT(State& s) const
    {
        using boost::numeric::ublas::prod;
        s.state = prod(transfer, s.state);
    }
};

#endif // SCSI_LINEAR_H
