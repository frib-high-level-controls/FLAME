
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "flame/constants.h"
#include "flame/moment.h"

#define sqr(x)  ((x)*(x))
#define cube(x) ((x)*(x)*(x))

#ifdef DEFPATH
    #define defpath DEFPATH
#else
    #define defpath "."
#endif

std::map<std::string,boost::shared_ptr<Config> > CurveMap;

// http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
// by LU-decomposition.
void inverse(MomentElementBase::value_t& out, const MomentElementBase::value_t& in)
{
    using boost::numeric::ublas::permutation_matrix;
    using boost::numeric::ublas::lu_factorize;
    using boost::numeric::ublas::lu_substitute;
    using boost::numeric::ublas::identity_matrix;

    MomentElementBase::value_t scratch(in); // copy
    permutation_matrix<size_t> pm(scratch.size1());
    if(lu_factorize(scratch, pm)!=0)
        throw std::runtime_error("Failed to invert matrix");
    out.assign(identity_matrix<double>(scratch.size1()));
    //out = identity_matrix<double>(scratch.size1());
    lu_substitute(scratch, pm, out);
}

void RotMat(const double dx, const double dy,
            const double theta_x, const double theta_y, const double theta_z,
            typename MomentElementBase::value_t &R)
{
    typedef typename MomentElementBase::state_t state_t;

    MomentState::matrix_t T = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

    R = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

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

    T(0, 6) = -dx, T(2, 6) = -dy;

    R = prod(R, T);
}


void GetQuadMatrix(const double L, const double K, const unsigned ind, typename MomentElementBase::value_t &M)
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

void GetSextMatrix(const double L, const double K3, double Dx, double Dy,
                   const double D2x, const double D2y, const double D2xy, const bool thinlens, const bool dstkick, typename MomentElementBase::value_t &M)
{
    typedef typename MomentElementBase::state_t state_t;
    // 2D sextupole transport matrix.
    double sqrtK, psi, cs, sn, ch, sh,
           dr = sqrt(sqr(Dx)+sqr(Dx));

    if (thinlens) {

        //thin-lens (drift-kick-drift) model

        MomentState::matrix_t T = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
        MomentState::matrix_t P = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);
        MomentState::matrix_t scratch;

        T(state_t::PS_X, state_t::PS_PX) = L/2e0;
        T(state_t::PS_Y, state_t::PS_PY) = L/2e0;

        P(state_t::PS_PX, state_t::PS_X) = -K3*L*Dx;
        P(state_t::PS_PX, state_t::PS_Y) =  K3*L*Dy;

        P(state_t::PS_PY, state_t::PS_X) = K3*L*Dy;
        P(state_t::PS_PY, state_t::PS_Y) = K3*L*Dx;

        scratch = prod(P,T);
        M = prod(T,scratch);

    } else {

        //thick-lens model

        sqrtK = sqrt(fabs(K3)*dr);
        psi = sqrtK*L;

        // Focusing or Defocusing switch
        Dx *= copysign(1e0,K3);
        Dy *= copysign(1e0,K3);

        cs = ::cos(psi);
        sn = ::sin(psi);
        ch = ::cosh(psi);
        sh = ::sinh(psi);


        if (sqrtK != 0e0) {
            M(state_t::PS_X, state_t::PS_X) = M(state_t::PS_PX, state_t::PS_PX) = ((dr+Dx)*cs+(dr-Dx)*ch)/(2e0*dr);
            M(state_t::PS_X, state_t::PS_PX) = ((dr+Dx)*sn+(dr-Dx)*sh)/(2e0*sqrtK*dr);
            M(state_t::PS_X, state_t::PS_Y)  = M(state_t::PS_PX, state_t::PS_PY) =
            M(state_t::PS_Y, state_t::PS_X)  = M(state_t::PS_PY, state_t::PS_PX) = Dy*(-cs+ch)/(2e0*dr);
            M(state_t::PS_X, state_t::PS_PY) = Dy*(-sn+sh)/(2e0*sqrtK*dr);

            M(state_t::PS_PX, state_t::PS_X) = sqrtK*(-(dr+Dx)*sn+(dr-Dx)*sh)/(2e0*dr);
            M(state_t::PS_PX, state_t::PS_Y) = M(state_t::PS_PY, state_t::PS_X) =  Dy*sqrtK*(sn+sh)/(2e0*dr);

            M(state_t::PS_Y, state_t::PS_PX) = Dy*(-sn+sh)/(2e0*sqrtK*dr);
            M(state_t::PS_Y, state_t::PS_Y) = M(state_t::PS_PY, state_t::PS_PY) = ((dr-Dx)*cs+(dr+Dx)*ch)/(2e0*dr);
            M(state_t::PS_Y, state_t::PS_PY) = ((dr-Dx)*sn+(dr+Dx)*sh)/(2e0*sqrtK*dr);

            M(state_t::PS_PY, state_t::PS_Y) = sqrtK*(-(dr-Dx)*sn+(dr+Dx)*sh)/(2e0*dr);

        } else {
            M(state_t::PS_X, state_t::PS_PX) = L;
            M(state_t::PS_Y, state_t::PS_PY) = L;
        }

    } // option
    if (dstkick){
        M(state_t::PS_PX, 6) = -K3*L*(D2x-D2y);
        M(state_t::PS_PY, 6) =  2e0*K3*L*D2xy;
    }
}

void GetEdgeMatrix(const double rho, const double phi, typename MomentElementBase::value_t &M)
{
    typedef typename MomentElementBase::state_t state_t;

    M = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

    M(state_t::PS_PX, state_t::PS_X) =  tan(phi)/rho;
    M(state_t::PS_PY, state_t::PS_Y) = -tan(phi)/rho;
}


void GetEEdgeMatrix(const double fringe_x, const double fringe_y, const double kappa, typename MomentElementBase::value_t &M)
{
    // Edge focusing for electrostatic dipole.
    typedef typename MomentElementBase::state_t state_t;

    M = boost::numeric::ublas::identity_matrix<double>(state_t::maxsize);

    M(state_t::PS_PX, state_t::PS_X)  = fringe_x;
    M(state_t::PS_PX, state_t::PS_PX) = sqrt(1e0+kappa);
    M(state_t::PS_PY, state_t::PS_Y)  = fringe_y;
    M(state_t::PS_PS, state_t::PS_PS) = 1e0+kappa;
}


void GetSBendMatrix(const double L, const double phi, const double phi1, const double phi2, const double K,
                    const double IonEs, const double ref_gamma, const double qmrel,
                    const double dip_beta, const double dip_gamma, const double d, const double dip_IonK, typename MomentElementBase::value_t &M)
{
    typedef typename MomentElementBase::state_t state_t;

    MomentState::matrix_t edge1, edge2, R;

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
        dx = sqr(L)/(2e0*rho);
        sx = L/rho;
    } else if (Kx > 0e0) {
        dx = (1e0-cos(sqrt(Kx)*L))/(rho*Kx);
        sx = sin(sqrt(Kx)*L)/(rho*sqrt(Kx));
    } else {
        dx = (1e0-cosh(sqrt(-Kx)*L))/(rho*Kx);
        sx = sin(sqrt(-Kx)*L)/(rho*sqrt(-Kx));
    }

    M(state_t::PS_X,  state_t::PS_PS) = dx/(sqr(dip_beta)*dip_gamma*IonEs/MeVtoeV);
    M(state_t::PS_PX, state_t::PS_PS) = sx/(sqr(dip_beta)*dip_gamma*IonEs/MeVtoeV);

    M(state_t::PS_S,  state_t::PS_X)  = sx*dip_IonK;
    M(state_t::PS_S,  state_t::PS_PX) = dx*dip_IonK;
    // Low beta approximation.
    M(state_t::PS_S,  state_t::PS_PS) =
            ((L/rho-sx)/(Kx*rho)-L/sqr(ref_gamma))*dip_IonK
            /(sqr(dip_beta)*dip_gamma*IonEs/MeVtoeV);

    // Add dipole terms.
    M(state_t::PS_S,  6) = ((L/rho-sx)/(Kx*rho)*d-L/sqr(ref_gamma)*(d+qmrel))*dip_IonK;
    M(state_t::PS_X,  6) = dx*d;
    M(state_t::PS_PX, 6) = sx*d;

    // Edge focusing.
    GetEdgeMatrix(rho, phi1, edge1);
    GetEdgeMatrix(rho, phi2, edge2);

    M = prod(M, edge1);
    M = prod(edge2, M);

    // Longitudinal plane.
    // For total path length.
    //        M(state_t::PS_S,  state_t::PS_S) = L;
}


void GetSolMatrix(const double L, const double K, typename MomentElementBase::value_t &M)
{
    typedef typename MomentElementBase::state_t state_t;

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
                    const double dip_beta, const double dip_gamma, const double delta_KZ, const double SampleIonK, typename MomentElementBase::value_t &M)
{
    typedef typename MomentElementBase::state_t state_t;

    MomentState::matrix_t edge;

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

//    M(state_t::PS_X,  6) = dx*Nk*(delta_K+delta_KZ);
//    M(state_t::PS_PX, 6) = sx/rho*Nk*(delta_K+delta_KZ);

    double Nkz = (1e0+2e0*eta0-h)/(2e0*(1e0+eta0)*(1e0+2e0*eta0));
    M(state_t::PS_X,  6) = dx*(Nk*delta_K+Nkz*delta_KZ);
    M(state_t::PS_PX, 6) = sx/rho*(Nk*delta_K+Nkz*delta_KZ);

    // Edge focusing.
    GetEEdgeMatrix(fringe_x, fringe_y ,kappa, edge);

    M = prod(M, edge);
    M = prod(edge, M);

    // Longitudinal plane.
    // For total path length.
    //        M(state_t::PS_S,  state_t::PS_S) = L;
}

void GetCurveData(const Config &c, const unsigned ncurve, std::vector<double> &Scales,
                  std::vector<std::vector<double> > &Curves)
{
    boost::shared_ptr<Config> conf;

    std::string filename;
    std::vector<double> range;
    bool checker = c.tryGet<std::string>("CurveFile", filename),
         rngchecker = c.tryGet<std::vector<double> >("use_range", range);

    if (checker){
        std::string CurveFile =  c.get<std::string>("Eng_Data_Dir", defpath);
        CurveFile += "/" + filename;
        std::string key(SB()<<CurveFile<<"|"<<boost::filesystem::last_write_time(CurveFile));
        if ( CurveMap.find(key) == CurveMap.end() ) {
            // not found in CurveMap
            try {
                try {
                    GLPSParser P;
                    conf.reset(P.parse_file(CurveFile.c_str(), false));
                }catch(std::exception& e){
                    throw std::runtime_error(SB()<<"Parse error: "<<e.what()<<"\n");
                }

            }catch(std::exception& e){
                throw std::runtime_error(SB()<<"Error: "<<e.what()<<"\n");
            }
            CurveMap.insert(std::make_pair(key, conf));
        } else {
            // found in CurveMap
            conf=CurveMap[key];
        }
    }

    size_t prev_size = 0;
    for (unsigned n=0; n<ncurve; n++) {
        std::string num(boost::lexical_cast<std::string>(n));
        std::vector<double> cv;
        bool cvchecker;
        if (checker){
            cvchecker = conf->tryGet<std::vector<double> >("curve"+num, cv);
        } else {
            cvchecker = c.tryGet<std::vector<double> >("curve"+num, cv);
        }
        if (!cvchecker)throw std::runtime_error(SB()<<"'curve" << num << "' is missing in lattice file.\n");

        if (rngchecker){
            unsigned start, end;
            if (range.size() != 2)
                throw std::runtime_error(SB()<<"Size of 'use_range' must be 2.\n");
            start = static_cast<unsigned>(range[0]);
            end   = static_cast<unsigned>(range[1]);

            if (start > cv.size() or end > cv.size())
                throw std::runtime_error(SB()<<"'use_range' is out of curve size.\n");

            std::vector<double> part(cv.begin()+start, cv.begin()+end);
            Curves.push_back(part);
        } else {
            Curves.push_back(cv);
        }

        if (n != 0 and prev_size != cv.size())
            throw std::runtime_error(SB()<<"Size of 'curve" << n << "' (" << cv.size() <<") and 'curve" << n-1 <<
                                           "' (" << prev_size <<") are inconsistent.  All curves must have the same size.\n");
        prev_size = cv.size();

        Scales.push_back(c.get<double>("scl_fac"+num, 0.0));
    }
}