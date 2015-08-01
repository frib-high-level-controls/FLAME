
#include "scsi/linear.h"
#include "scsi/state/vector.h"
#include "scsi/state/matrix.h"

MatrixState::MatrixState(const Config& c)
    :StateBase(c)
    ,state(boost::numeric::ublas::identity_matrix<double>(2))
{}

MatrixState::~MatrixState() {}

void MatrixState::show(std::ostream& strm) const
{
    strm<<"State: "<<state<<"\n";
}

const char* MatrixState::type_name() {return "TransferMatrix";}

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
    ,state(2, 0.0)
{}

VectorState::~VectorState() {}

void VectorState::show(std::ostream& strm) const
{
    strm<<"State: "<<state<<"\n";
}

const char* VectorState::type_name() {return "Vector";}

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
    Machine::registerState<VectorState>();
    Machine::registerState<MatrixState>();

    Machine::registerElement<LinearDrift<VectorState> >("drift");
    Machine::registerElement<LinearDrift<MatrixState> >("drift");

    Machine::registerElement<LinearThinDipole<VectorState> >("dipole");
    Machine::registerElement<LinearThinDipole<MatrixState> >("dipole");

    Machine::registerElement<LinearThinQuad<VectorState> >("quad");
    Machine::registerElement<LinearThinQuad<MatrixState> >("quad");
}
