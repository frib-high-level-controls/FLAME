#ifndef CAV_TTF_H
#define CAV_TTF_H

#endif // CAV_TTF_H

#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <string>

// Phase space dimension; including vector for orbit/1st moment.
# define PS_Dim Moment2State::maxsize // Set to 7; to include orbit.

void calTransfac(const int cavi, const double IonK,
                 double &Ecenter, double &Trans, double &Transp, double &Strans, double &Stransp, double &V0);
