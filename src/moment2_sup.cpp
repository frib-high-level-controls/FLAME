
#include "scsi/constants.h"
#include "scsi/moment2.h"


// HdipoleFitMode = true    dipole strength adjusted to beam energy.
bool HdipoleFitMode = true;


void GetQuadMatrix(const double L, const double K, const unsigned ind, typename Moment2ElementBase::value_t &M)
{
    // 2D quadrupole transport matrix.
    double sqrtK,
           psi,
           cs,
           sn;

    if (K > 0e0) {
        // Focusing.
        sqrtK = sqrt(K);
        psi = sqrtK*L;
        cs = ::cos(psi);
        sn = ::sin(psi);

        M(ind, ind) = M(ind+1, ind+1) = cs;
        if (sqrtK != 0e0)
            M(ind, ind+1) = sn/sqrtK;
        else
            M(ind, ind+1) = L;
        if (sqrtK != 0e0)
            M(ind+1, ind) = -sqrtK*sn;
        else
            M(ind+1, ind) = 0e0;
    } else {
        // Defocusing.
        sqrtK = sqrt(-K);
        psi = sqrtK*L;
        cs = ::cosh(psi);
        sn = ::sinh(psi);

        M(ind, ind) = M(ind+1, ind+1) = cs;
        if (sqrtK != 0e0)
            M(ind, ind+1) = sn/sqrtK;
        else
            M(ind, ind+1) = L;
        if (sqrtK != 0e0)
            M(ind+1, ind) = sqrtK*sn;
        else
            M(ind+1, ind) = 0e0;
    }
}


void GetEdgeMatrix(const double rho, const double phi, typename Moment2ElementBase::value_t &M)
{
    typedef typename Moment2ElementBase::state_t state_t;

    M = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

    M(state_t::PS_PX, state_t::PS_X) =  tan(phi)/rho;
    M(state_t::PS_PY, state_t::PS_Y) = -tan(phi)/rho;
}


void GetSBendMatrix(const double L, const double phi, const double phi1, const double phi2, const double K,
                    const double IonEs, const double ref_gamma, const double qmrel, const double dip_beta,
                    const double dip_gamma, const double d, const double dip_IonK, typename Moment2ElementBase::value_t &M)
{
    typedef typename Moment2ElementBase::state_t state_t;

    double  rho    = L/phi,
            Kx     = K + 1e0/sqr(rho),
            Ky     = -K,
            dx     = 0e0,
            sx     = 0e0,
            dip_bg = dip_beta*dip_gamma;

    typename Moment2ElementBase::value_t edge1, edge2;

    // Edge focusing.
    GetEdgeMatrix(rho, phi1, edge1);
    // Horizontal plane.
    GetQuadMatrix(L, Kx, (unsigned)state_t::PS_X, M);
    // Vertical plane.
    GetQuadMatrix(L, Ky, (unsigned)state_t::PS_Y, M);

    // Include dispersion.
    if (Kx == 0e0) {
        dx = sqr(L)/2e0;
        sx = L;
    } else if (Kx > 0e0) {
        dx = (1e0-cos(sqrt(Kx)*L))/Kx;
        sx = sin(sqrt(Kx)*L)/sqrt(Kx);
    } else {
        dx = (1e0-cosh(sqrt(-Kx)*L))/Kx;
        sx = sin(sqrt(Kx)*L)/sqrt(Kx);
    }

    M(state_t::PS_X,  state_t::PS_PS) = dx/(rho*sqr(dip_beta)*dip_gamma*IonEs/MeVtoeV);
    M(state_t::PS_PX, state_t::PS_PS) = sx/(rho*sqr(dip_beta)*dip_gamma*IonEs/MeVtoeV);

    M(state_t::PS_S,  state_t::PS_X)  = sx/rho*dip_IonK;
    M(state_t::PS_S,  state_t::PS_PX) = dx/rho*dip_IonK;

    // Low beta approximation.
    M(state_t::PS_S,  state_t::PS_PS) =
            ((L-sx)/(Kx*sqr(rho))-L/sqr(ref_gamma))*dip_IonK
            /(sqr(dip_beta)*dip_gamma*IonEs/MeVtoeV);

    // Add dipole terms.
    M(state_t::PS_S,  6) = ((L-sx)/(Kx*sqr(rho))*d-L/sqr(ref_gamma)*(d+qmrel))*dip_IonK;

    // Add dipole terms.
    M(state_t::PS_X,  6) = dx/rho*d;
    M(state_t::PS_PX, 6) = sx/rho*d;

    // Edge focusing.
    GetEdgeMatrix(rho, phi2, edge2);

    M = prod(M, edge1);
    M = prod(edge2, M);

    // Longitudinal plane.
    // For total path length.
    //        M(state_t::PS_S,  state_t::PS_S) = L;
}


void GetSolMatrix(const double L, const double K, typename Moment2ElementBase::value_t &M)
{
    typedef typename Moment2ElementBase::state_t state_t;

    double C = ::cos(K*L),
           S = ::sin(K*L);

    M(state_t::PS_X, state_t::PS_X)
            = M(state_t::PS_PX, state_t::PS_PX)
            = M(state_t::PS_Y, state_t::PS_Y)
            = M(state_t::PS_PY, state_t::PS_PY)
            = sqr(C);

    if (K != 0e0)
        M(state_t::PS_X, state_t::PS_PX) = S*C/K;
    else
        M(state_t::PS_X, state_t::PS_PX) = L;
    M(state_t::PS_X, state_t::PS_Y) = S*C;
    if (K != 0e0)
        M(state_t::PS_X, state_t::PS_PY) = sqr(S)/K;
    else
        M(state_t::PS_X, state_t::PS_PY) = 0e0;

    M(state_t::PS_PX, state_t::PS_X) = -K*S*C;
    M(state_t::PS_PX, state_t::PS_Y) = -K*sqr(S);
    M(state_t::PS_PX, state_t::PS_PY) = S*C;

    M(state_t::PS_Y, state_t::PS_X) = -S*C;
    if (K != 0e0)
        M(state_t::PS_Y, state_t::PS_PX) = -sqr(S)/K;
    else
        M(state_t::PS_Y, state_t::PS_PX) = 0e0;
    if (K != 0e0)
        M(state_t::PS_Y, state_t::PS_PY) = S*C/K;
    else
        M(state_t::PS_Y, state_t::PS_PY) = L;

    M(state_t::PS_PY, state_t::PS_X) = K*sqr(S);
    M(state_t::PS_PY, state_t::PS_PX) = -S*C;
    M(state_t::PS_PY, state_t::PS_Y) = -K*S*C;

    // Longitudinal plane.
    // For total path length.
//        M(state_t::PS_S, state_t::PS_S) = L;
}
