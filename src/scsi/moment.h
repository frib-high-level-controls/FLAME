#ifndef SCSI_MOMENT_H
#define SCSI_MOMENT_H

#include <ostream>
#include <math.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "base.h"
#include "state/matrix.h"

/** @brief An Element which propogates the statistical moments of a bunch
 *
 */
struct MomentElementBase : public ElementVoid
{
    typedef MatrixState state_t;

    MomentElementBase(const Config& c);
    virtual ~MomentElementBase();

    virtual void advance(StateBase& s) const;

    virtual void show(std::ostream& strm) const;

    typedef boost::numeric::ublas::matrix<double> value_t;

    value_t transfer;
    //value_t transferT;
};

#endif // SCSI_MOMENT_H
