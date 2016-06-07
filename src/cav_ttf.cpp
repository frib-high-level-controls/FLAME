
#include <boost/numeric/ublas/lu.hpp>

#include "scsi/constants.h"
#include "scsi/moment2.h"
#include "scsi/cav_ttf.h"

void calTransfac(const int cavi, const double IonK,
                 double &Ecenter, double &T, double &Tp, double &S, double &Sp, double &V0)
{
//    int    k;
//    double ezData_shift, ecenterup, ecenterdown, transup, transpup;
//    double dis[n], dz[n], entry[n];

//    dis = axisData.get(axisData.size()-1)[0] - axisData.get(0)[0];
//    dz  = dis/(axisData.size()-1); //mm

//    // rewrite to let axisData starts at 0
//    ezData_shift = axisData.get(0)[0];
//    k = 0;
//    for (k = 0; k < axisData.size(); k++) {
//        entry = {axisData.get(k)[0]-ezData_shift, axisData.get(k)[1]};
//        axisData.set(k, entry);
//    }

//    // calculate for electric geometry center from start point
//    ecenterup   = 0e0;
//    ecenterdown = 0e0;
//    for (k = 0; k < axisData.size()-1; k++) {
//        ecenterup +=
//                (fabs(axisData.get(k)[1])+fabs(axisData.get(k+1)[1]))/2e0
//                *(axisData.get(k)[0]+axisData.get(k+1)[0])/2*dz;
//        ecenterdown += (fabs(axisData.get(k)[1])+fabs(axisData.get(k+1)[1]))/2e0*dz;
//    }
//    Ecenter = ecenterup/ecenterdown;

//    // previous was axisData.size()-1, I think that's a bug
//    for (k = 0; k < axisData.size(); k++) {
//        entry = {axisData.get(k)[0]-Ecenter, axisData.get(k)[1]};
//        axisData.set(k, entry);
//    }

//    // calculate for time transit factor
//    transup = 0;
//    for (k = 0; k < axisData.size()-1; k++) {
//        transup +=
//                (axisData.get(k)[1]+axisData.get(k+1)[1])/2e0
//                *cos(IonK*(axisData.get(k)[0]+axisData.get(k+1)[0])/2e0)*dz;
//    }
//    T = transup/ecenterdown;

//    // calculate for time transit factor prime dTrans/dk
//    transpup = 0;
//    for (k = 0; k < axisData.size()-1; k++) {
//        transpup -=
//                ((axisData.get(k)[0]+axisData.get(k+1)[0])/2e0)*(axisData.get(k)[1]+axisData.get(k+1)[1])/2e0
//                *sin(IonK*(axisData.get(k)[0]+axisData.get(k+1)[0])/2)*dz;
//    }
//    Tp = transpup/ecenterdown;

//    // calculate for time S transit factor
//    stransup = 0e0;
//    for (k = 0; k < axisData.size()-1; k++) {
//        stransup +=
//                (axisData.get(k)[1]+axisData.get(k+1)[1])/2e0
//                *sin(k*(axisData.get(k)[0]+axisData.get(k+1)[0])/2e0)*dz;
//    }
//    S = stransup/ecenterdown;

//    // calculate for time S transit factor prime dStrans/dk
//    stranspup = 0e0;
//    for (k = 0; k < axisData.size()-1; k++) {
//        stranspup +=
//                (axisData.get(k)[0]+axisData.get(k+1)[0])/2e0
//                *(axisData.get(k)[1]+axisData.get(k+1)[1])/2e0
//                *cos(IonK*(axisData.get(k)[0]+axisData.get(k+1)[0])/2)*dz;
//    }
//    Sp = stranspup/ecenterdown;

//    V0 = ecenterdown/1e6/1000; //MV
}

