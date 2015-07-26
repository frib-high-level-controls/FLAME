
#include "scsi/linear.h"

LinearMatrixState::~LinearMatrixState() {}

LinearVectorState::~LinearVectorState() {}

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
