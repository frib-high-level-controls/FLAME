#ifndef MOMENT2_SUP_H
#define MOMENT2_SUP_H

#endif // MOMENT2_SUP_H


void inverse(Moment2ElementBase::value_t& out, const Moment2ElementBase::value_t& in);

void PrtVec1(const Moment2State::vector_t &a);

void PrtMat1(const value_mat &M);

void RotMat(const double dx, const double dy,
            const double theta_x, const double theta_y, const double theta_z,
            typename Moment2ElementBase::value_t &R);

void GetQuadMatrix(const double L, const double K, const unsigned ind, typename Moment2ElementBase::value_t &M);

void GetEdgeMatrix(const double rho, const double phi, typename Moment2ElementBase::value_t &M);

void GetEEdgeMatrix(const double fringe_x, const double fringe_y, const double kappa, typename Moment2ElementBase::value_t &M);

void GetSBendMatrix(const double L, const double phi, const double phi1, const double phi2, const double K,
                    const double IonEs, const double ref_gamma, const double qmrel,
                    const double dip_beta, const double dip_gamma, const double d, const double dip_IonK, typename Moment2ElementBase::value_t &M);

void GetSolMatrix(const double L, const double K, typename Moment2ElementBase::value_t &M);


void GetEBendMatrix(const double L, const double phi, const double fringe_x, const double fringe_y, const double kappa, const double Kx, const double Ky,
                    const double IonEs, const double ref_beta, const double real_gamma, const double eta0, const double h, const double dip_beta,
                    const double dip_gamma, const double delta_KZ, const double SampleIonK, typename Moment2ElementBase::value_t &M);
