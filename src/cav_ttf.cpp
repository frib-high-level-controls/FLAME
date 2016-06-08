
#include <boost/numeric/ublas/lu.hpp>

#include "scsi/constants.h"
#include "scsi/moment2.h"
#include "scsi/cav_ttf.h"


void calTransfac(const int cavi, const double IonK,
                 double &Ecenter, double &T, double &Tp, double &S, double &Sp, double &V0)
{
    int    n, k;
    double eml, em_mom;

    std::vector<double> z, EM;

    n = z.size();

    // Start at zero.
    for (k = 0; k < n; k++)
        z[k] -= z[0];

    eml = 0e0, em_mom = 0e0;
    for (k = 0; k < n-1; k++) {
        eml    += (fabs(EM[k])+fabs(EM[k+1]))/2e0*(z[k+1]-z[k]);
        em_mom += (z[k]+z[k+1])/2e0*(fabs(EM[k])+fabs(EM[k+1]))/2e0*(z[k+1]-z[k]);
    }
    Ecenter = em_mom/eml;

    T = 0;
    for (k = 0; k < n-1; k++)
        T += (EM[k]+EM[k+1])/2e0*cos(IonK*(z[k]+z[k+1])/2e0)*(z[k+1]-z[k]);
    T /= eml;

    Tp = 0;
    for (k = 0; k < n-1; k++)
        Tp -= ((z[k]+z[k+1])/2e0)*(EM[k]+EM[k+1])/2e0*sin(IonK*(z[k]+z[k+1])/2e0)*(z[k+1]-z[k]);
    Tp /= eml;

    S = 0e0;
    for (k = 0; k < n-1; k++)
        S += (EM[k]+EM[k+1])/2e0*sin(k*(z[k]+z[k+1])/2e0)*(z[k+1]-z[k]);
    S /= eml;

    Sp = 0e0;
    for (k = 0; k < n-1; k++) {
        Sp += (z[k]+z[k+1])/2e0*(EM[k]+EM[k+1])/2e0*cos(IonK*(z[k]+z[k+1])/2e0)*(z[k+1]-z[k]);
    }
    Sp /= eml;

    V0 = eml/MeVtoeV/MtoMM;
}

