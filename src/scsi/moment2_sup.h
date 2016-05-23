#ifndef MOMENT2_SUP_H
#define MOMENT2_SUP_H

#endif // MOMENT2_SUP_H


extern bool HdipoleFitMode;


void GetQuadMatrix(const double L, const double K, const unsigned ind, typename Moment2ElementBase::value_t &M);

void GetEdgeMatrix(const double rho, const double phi, typename Moment2ElementBase::value_t &M);

void GetEEdgeMatrix(const double rho, const double phi, typename Moment2ElementBase::value_t &M);

void GetSBendMatrix(const double L, const double phi, const double phi1, const double phi2, const double Kx, const double Ky,
                    const double IonEs, const double ref_gamma, const double qmrel, const double dip_beta,
                    const double dip_gamma, const double d, const double dip_IonK, typename Moment2ElementBase::value_t &M);

void GetSolMatrix(const double L, const double K, typename Moment2ElementBase::value_t &M);
