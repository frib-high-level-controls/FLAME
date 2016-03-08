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
    enum {maxsize=7};
    enum param_t {
        PS_X, PS_PX, PS_Y, PS_PY, PS_S, PS_PS
    };

    MatrixState(const Config& c);
    virtual ~MatrixState();

    void assign(const StateBase& other);

    typedef boost::numeric::ublas::matrix<double,
                    boost::numeric::ublas::row_major,
                    boost::numeric::ublas::bounded_array<double, maxsize*maxsize>
    > value_t;

    virtual void show(std::ostream& strm) const;

    value_t state;

    virtual bool getArray(unsigned idx, ArrayInfo& Info);

    virtual MatrixState* clone() const {
        return new MatrixState(*this, clone_tag());
    }

protected:
    MatrixState(const MatrixState& o, clone_tag);
};

#endif // SCSI_STATE_MATRIX_HPP
