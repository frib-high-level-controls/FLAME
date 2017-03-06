#ifndef MOMENT2_SUP_H
#define MOMENT2_SUP_H

#include <vector>
#include "base.h"
#include "moment.h"

void inverse(MomentElementBase::value_t& out, const MomentElementBase::value_t& in);

void RotMat(const double dx, const double dy,
            const double theta_x, const double theta_y, const double theta_z,
            typename MomentElementBase::value_t &R);

void GetQuadMatrix(const double L, const double K, const unsigned ind, typename MomentElementBase::value_t &M);

void GetSextMatrix(const double L, const double K, double Dx, double Dy,
                   const double D2x, const double D2y, const double D2xy, const bool thinlens, const bool dstkick, typename MomentElementBase::value_t &M);

void GetEdgeMatrix(const double rho, const double phi, typename MomentElementBase::value_t &M);

void GetEEdgeMatrix(const double fringe_x, const double fringe_y, const double kappa, typename MomentElementBase::value_t &M);

void GetSBendMatrix(const double L, const double phi, const double phi1, const double phi2, const double K,
                    const double IonEs, const double ref_gamma, const double qmrel,
                    const double dip_beta, const double dip_gamma, const double d, const double dip_IonK, typename MomentElementBase::value_t &M);

void GetSolMatrix(const double L, const double K, typename MomentElementBase::value_t &M);


void GetEBendMatrix(const double L, const double phi, const double fringe_x, const double fringe_y, const double kappa, const double Kx, const double Ky,
                    const double IonEs, const double ref_beta, const double real_gamma, const double eta0, const double h, const double dip_beta,
                    const double dip_gamma, const double delta_KZ, const double SampleIonK, typename MomentElementBase::value_t &M);

#endif // MOMENT2_SUP_H
