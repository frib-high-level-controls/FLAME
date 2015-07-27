
#include "scsi/linear.h"

LinearMatrixState::LinearMatrixState(const Config& c)
    :StateBase(c)
    ,state(boost::numeric::ublas::identity_matrix<double>(2))
{}

LinearMatrixState::~LinearMatrixState() {}

void LinearMatrixState::show(std::ostream& strm) const
{
    strm<<"State: "<<state<<"\n";
}

const char* LinearMatrixState::type_name() {return "Linear1DTransfer";}

bool LinearMatrixState::getArray(unsigned idx, ArrayInfo& Info) {
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

LinearVectorState::LinearVectorState(const Config& c)
    :StateBase(c)
    ,state(2, 0.0)
{}

LinearVectorState::~LinearVectorState() {}

void LinearVectorState::show(std::ostream& strm) const
{
    strm<<"State: "<<state<<"\n";
}

const char* LinearVectorState::type_name() {return "Linear1D";}

bool LinearVectorState::getArray(unsigned idx, ArrayInfo& Info) {
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
    Machine::registerState<LinearVectorState>();
    Machine::registerState<LinearMatrixState>();

    Machine::registerElement<LinearDrift<LinearVectorState> >("drift");
    Machine::registerElement<LinearDrift<LinearMatrixState> >("drift");

    Machine::registerElement<LinearThinDipole<LinearVectorState> >("dipole");
    Machine::registerElement<LinearThinDipole<LinearMatrixState> >("dipole");

    Machine::registerElement<LinearThinQuad<LinearVectorState> >("quad");
    Machine::registerElement<LinearThinQuad<LinearMatrixState> >("quad");
}
