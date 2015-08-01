
#include "scsi/linear.h"
#include "scsi/state/vector.h"
#include "scsi/state/matrix.h"

MatrixState::MatrixState(const Config& c)
    :StateBase(c)
    ,state(boost::numeric::ublas::identity_matrix<double>(6))
{}

MatrixState::~MatrixState() {}

void MatrixState::show(std::ostream& strm) const
{
    strm<<"State: "<<state<<"\n";
}

bool MatrixState::getArray(unsigned idx, ArrayInfo& Info) {
    if(idx==0) {
        Info.name = "state";
        Info.ptr = &state(0,0);
        Info.ndim = 2;
        Info.dim[0] = state.size1();
        Info.dim[1] = state.size2();
        return true;
    }
    return false;
}

VectorState::VectorState(const Config& c)
    :StateBase(c)
    ,state(6, 0.0)
{}

VectorState::~VectorState() {}

void VectorState::show(std::ostream& strm) const
{
    strm<<"State: "<<state<<"\n";
}

bool VectorState::getArray(unsigned idx, ArrayInfo& Info) {
    if(idx==0) {
        Info.name = "state";
        Info.ptr = &state(0);
        Info.ndim = 1;
        Info.dim[0] = state.size();
        return true;
    }
    return false;
}

void registerLinear()
{
    Machine::registerState<VectorState>("Vector");
    Machine::registerState<MatrixState>("TransferMatrix");

    Machine::registerElement<LinearDrift<VectorState> >("Vector", "drift");
    Machine::registerElement<LinearDrift<MatrixState> >("TransferMatrix", "drift");

    Machine::registerElement<LinearThinDipole<VectorState> >("Vector", "dipole");
    Machine::registerElement<LinearThinDipole<MatrixState> >("TransferMatrix", "dipole");

    Machine::registerElement<LinearThinQuad<VectorState> >("Vector", "quad");
    Machine::registerElement<LinearThinQuad<MatrixState> >("TransferMatrix", "quad");
}
