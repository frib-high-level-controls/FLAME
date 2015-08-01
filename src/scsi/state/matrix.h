#ifndef SCSI_STATE_MATRIX_HPP
#define SCSI_STATE_MATRIX_HPP

#include <ostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include "../base.h"

/** @brief Simulation state which include only a matrix
 */
struct MatrixState : public StateBase
{
    enum {maxsize=6};

    MatrixState(const Config& c);
    virtual ~MatrixState();

    typedef boost::numeric::ublas::matrix<double,
                    boost::numeric::ublas::row_major,
                    boost::numeric::ublas::bounded_array<double, maxsize*maxsize>
    > value_t;

    virtual void show(std::ostream& strm) const;

    value_t state;

    static const char* type_name();

    virtual bool getArray(unsigned idx, ArrayInfo& Info);
};

#endif // SCSI_STATE_MATRIX_HPP
