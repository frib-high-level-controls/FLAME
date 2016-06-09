
#include <boost/numeric/ublas/lu.hpp>

#include "scsi/constants.h"
#include "scsi/moment2.h"


// http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
// by LU-decomposition.
void inverse(Moment2ElementBase::value_t& out, const Moment2ElementBase::value_t& in)
{
    using boost::numeric::ublas::permutation_matrix;
    using boost::numeric::ublas::lu_factorize;
    using boost::numeric::ublas::lu_substitute;
    using boost::numeric::ublas::identity_matrix;

    Moment2ElementBase::value_t scratch(in); // copy
    permutation_matrix<size_t> pm(scratch.size1());
    if(lu_factorize(scratch, pm)!=0)
        throw std::runtime_error("Failed to invert matrix");
    out.assign(identity_matrix<double>(scratch.size1()));
    //out = identity_matrix<double>(scratch.size1());
    lu_substitute(scratch, pm, out);
}


void PrtVec1(const Moment2State::vector_t &a)
{
    for (size_t k = 0; k < a.size(); k++)
        std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << a[k];
    std::cout << "\n";
}


void PrtMat1(const value_mat &M)
{
    for (size_t j = 0; j < M.size1(); j++) {
        for (size_t k = 0; k < M.size2(); k++)
            std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << M(j, k);
        std::cout << "\n";
    }
}


void RotMat(const double dx, const double dy,
            const double theta_x, const double theta_y, const double theta_z,
            typename Moment2ElementBase::value_t &R)
{
    typedef typename Moment2ElementBase::state_t state_t;

    value_mat T = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

    R = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

#if 1

    // Left-handed coordinate system => theta_y -> -theta_y.

    double m11 =  cos(-theta_y)*cos(theta_z),
           m12 =  sin(theta_x)*sin(-theta_y)*cos(theta_z) + cos(theta_x)*sin(theta_z),
           m13 = -cos(theta_x)*sin(-theta_y)*cos(theta_z) + sin(theta_x)*sin(theta_z),

           m21 = -cos(-theta_y)*sin(theta_z),
           m22 = -sin(theta_x)*sin(-theta_y)*sin(theta_z) + cos(theta_x)*cos(theta_z),
           m23 =  cos(theta_x)*sin(-theta_y)*sin(theta_z) + sin(theta_x)*cos(theta_z),

           m31 =  sin(-theta_y),
           m32 = -sin(theta_x)*cos(-theta_y),
           m33 =  cos(theta_x)*cos(-theta_y);

    R(0, 0) = m11, R(0, 2) = m12, R(0, 4) = m13;
    R(2, 0) = m21, R(2, 2) = m22, R(2, 4) = m23;
    R(4, 0) = m31, R(4, 2) = m32, R(4, 4) = m33;

    R(1, 1) = m11, R(1, 3) = m12, R(1, 5) = m13;
    R(3, 1) = m21, R(3, 3) = m22, R(3, 5) = m23;
    R(5, 1) = m31, R(5, 3) = m32, R(5, 5) = m33;

#else

    value_mat Rx = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize),
              Ry = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize),
              Rz = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

    Rx(2, 2) =  cos(theta_x), Rx(2, 4) =  sin(theta_x);
    Rx(4, 2) = -sin(theta_x), Rx(4, 4) =  cos(theta_x);

    Ry(0, 0) =  cos(theta_y), Ry(0, 4) =  sin(theta_y);
    Ry(4, 0) = -sin(theta_y), Ry(4, 4) =  cos(theta_y);

    Rz(0, 0) =  cos(theta_z), Rz(0, 2) =  sin(theta_z);
    Rz(2, 0) = -sin(theta_z), Rz(2, 2) =  cos(theta_z),

    Rx(3, 3) =  cos(theta_x), Rx(3, 5) =  sin(theta_x);
    Rx(5, 3) = -sin(theta_x), Rx(5, 5) =  cos(theta_x),

    Ry(1, 1) =  cos(theta_y), Ry(1, 5) =  sin(theta_y);
    Ry(5, 1) = -sin(theta_y), Ry(5, 5) =  cos(theta_y);

    Rz(1, 1) =  cos(theta_z); Rz(1, 3) =  sin(theta_z);
    // J.B. Bug in TLM: should be Rz(3, 1) vs. Rz(2, 1).
    Rz(2, 1) = -sin(theta_z), Rz(3, 3) =  cos(theta_z);

    R = prod(Ry, Rx);
    R = prod(Rz, R);

#endif

    T(0, 6) = -dx, T(2, 6) = -dy;

    R = prod(R, T);
}


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


void GetEEdgeMatrix(const double fringe_x, const double fringe_y, const double kappa, typename Moment2ElementBase::value_t &M)
{
    // Edge focusing for electrostatic dipole.
    typedef typename Moment2ElementBase::state_t state_t;

    M = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

    M(state_t::PS_PX, state_t::PS_X)  = fringe_x;
    M(state_t::PS_PX, state_t::PS_PX) = sqrt(1e0+kappa);
    M(state_t::PS_PY, state_t::PS_Y)  = fringe_y;
    M(state_t::PS_PS, state_t::PS_PS) = 1e0+kappa;
}


void GetSBendMatrix(const double L, const double phi, const double phi1, const double phi2, const double K,
                    const double IonEs, const double ref_gamma, const double qmrel,
                    const double dip_beta, const double dip_gamma, const double d, const double dip_IonK, typename Moment2ElementBase::value_t &M)
{
    typedef typename Moment2ElementBase::state_t state_t;

    value_mat edge1, edge2, R;

    double  rho = L/phi,
            Kx  = K + 1e0/sqr(rho),
            Ky  = -K,
            dx  = 0e0,
            sx  = 0e0;

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
    M(state_t::PS_X,  6) = dx/rho*d;
    M(state_t::PS_PX, 6) = sx/rho*d;

    // Edge focusing.
    GetEdgeMatrix(rho, phi1, edge1);
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

    M(state_t::PS_PX, state_t::PS_X)  = -K*S*C;
    M(state_t::PS_PX, state_t::PS_Y)  = -K*sqr(S);
    M(state_t::PS_PX, state_t::PS_PY) = S*C;

    M(state_t::PS_Y, state_t::PS_X)   = -S*C;
    if (K != 0e0)
        M(state_t::PS_Y, state_t::PS_PX) = -sqr(S)/K;
    else
        M(state_t::PS_Y, state_t::PS_PX) = 0e0;
    if (K != 0e0)
        M(state_t::PS_Y, state_t::PS_PY) = S*C/K;
    else
        M(state_t::PS_Y, state_t::PS_PY) = L;

    M(state_t::PS_PY, state_t::PS_X)  = K*sqr(S);
    M(state_t::PS_PY, state_t::PS_PX) = -S*C;
    M(state_t::PS_PY, state_t::PS_Y)  = -K*S*C;

    // Longitudinal plane.
    // For total path length.
//        M(state_t::PS_S, state_t::PS_S) = L;
}


void GetEBendMatrix(const double L, const double phi, const double fringe_x, const double fringe_y, const double kappa, const double Kx, const double Ky,
                    const double IonEs, const double ref_beta, const double real_gamma, const double eta0, const double h,
                    const double dip_beta, const double dip_gamma, const double delta_KZ, const double SampleIonK, typename Moment2ElementBase::value_t &M)
{
    typedef typename Moment2ElementBase::state_t state_t;

    value_mat edge;

    double  rho = L/phi,
            scl = (real_gamma - 1e0)*IonEs/MeVtoeV,
            dx  = 0e0,
            sx  = 0e0;

    // Horizontal plane.
    GetQuadMatrix(L, Kx, (unsigned)state_t::PS_X, M);
    // Vertical plane.
    GetQuadMatrix(L, Ky, (unsigned)state_t::PS_Y, M);

    // Include dispersion.
    if (Kx == 0e0) {
        dx = 0e0;
        sx = 0e0;
    } else if (Kx > 0e0) {
        dx = (1e0-cos(sqrt(Kx)*L))/(rho*Kx);
        sx = sin(sqrt(Kx)*L)/sqrt(Kx);
    } else {
        dx = (1e0-cosh(sqrt(-Kx)*L))/(rho*Kx);
        sx = sin(sqrt(Kx)*L)/sqrt(Kx);
    }

    double Nk   = (sqr(1e0+2e0*eta0)+h)/(2e0*(1e0+eta0)*(1e0+2e0*eta0)),
           Nt   = 1e0 + h/sqr(1e0+2e0*eta0),
           CorT = -real_gamma/(1e0+real_gamma),
           tx   = sx/rho*Nt*CorT,
           txp  = dx*Nt*CorT,
           tzp  = (-L/(2e0*(1e0+eta0)*(1e0+2e0*eta0))+((L-sx)/Kx/sqr(rho))*Nk*Nt)*CorT;

    M(state_t::PS_X,  state_t::PS_PS) = dx*Nk/scl;
    M(state_t::PS_PX, state_t::PS_PS) = sx/rho*Nk/scl;

    M(state_t::PS_S,  state_t::PS_X)  = -tx*SampleIonK;
    M(state_t::PS_S,  state_t::PS_PX) = -txp*SampleIonK;
    // Low beta approximation.
    M(state_t::PS_S,  state_t::PS_PS) = -tzp*SampleIonK/scl;

    // Add dipole terms.
    double delta_K = sqr(ref_beta/dip_beta) - 1e0;

    M(state_t::PS_X,  6) = dx*Nk*(delta_K+delta_KZ);
    M(state_t::PS_PX, 6) = sx/rho*Nk*(delta_K+delta_KZ);

    // Edge focusing.
    GetEEdgeMatrix(fringe_x, fringe_y ,kappa, edge);

    M = prod(M, edge);
    M = prod(edge, M);

    // Longitudinal plane.
    // For total path length.
    //        M(state_t::PS_S,  state_t::PS_S) = L;
}
