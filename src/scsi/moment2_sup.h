#ifndef MOMENT2_SUP_H
#define MOMENT2_SUP_H

#endif // MOMENT2_SUP_H


void GetEdgeMatrix(const double rho, const double phi, typename Moment2ElementBase::value_t &M);

void GetQuadMatrix(const double L, const double K, const unsigned ind, typename Moment2ElementBase::value_t &M);

void GetSolMatrix(const double L, const double K, typename Moment2ElementBase::value_t &M);
