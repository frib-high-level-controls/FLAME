
#include <fstream>

#include "scsi/constants.h"
#include "scsi/moment2.h"
#include "scsi/moment2_sup.h"
#include "scsi/rf_cavity.h"

// RF Cavity beam dynamics functions.


void CavDataType::RdData(std::fstream &inf)
{
    std::string       line;
    double            s, Elong;
    std::stringstream str;

    while (getline(inf, line) && !inf.fail()) {
        str.str(line);
        str.clear();
        str >> s >> Elong;
        this->s.push_back(s), this->Elong.push_back(Elong);
        // Convert from [mm] to [m].
//        this->s[this->s.size()-1] /= MtoMM;
    }

    if (false) {
        std::cout << "\n";
        for (size_t k = 0; k < this->s.size(); k++)
            show(std::cout, k);
    }
}


void CavDataType::show(std::ostream& strm, const int k) const
{
    strm << std::scientific << std::setprecision(5)
         << std::setw(13) << s[k] << std::setw(13) << Elong[k] << "\n";
}


void CavDataType::show(std::ostream& strm) const
{
    for (unsigned int k = 0; k < s.size(); k++)
        show(strm, k);
}


void CavTLMLineType::clear(void)
{
    s.clear(); Elem.clear(); E0.clear();
    T.clear(); S.clear(); Accel.clear();
}


void CavTLMLineType::set(const double s, const std::string &Elem, const double E0,
                         const double T, const double S, const double Accel)
{
    this->s.push_back(s); this->Elem.push_back(Elem); this->E0.push_back(E0);
    this->T.push_back(T); this->S.push_back(S); this->Accel.push_back(Accel);
}


void CavTLMLineType::show(std::ostream& strm, const int k) const
{
    strm << std::fixed << std::setprecision(5)
         << std::setw(9) << s[k] << std::setw(10) << Elem[k]
         << std::scientific << std::setprecision(10)
         << std::setw(18) << T[k] << std::setw(18) << S[k]
         << std::setw(18) << Accel[k] << "\n";
}


void CavTLMLineType::show(std::ostream& strm) const
{
    for (unsigned int k = 0; k < s.size(); k++)
        show(strm, k);
}


void PrtVec(const std::vector<double> &a)
{
    for (size_t k = 0; k < a.size(); k++)
        std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << a[k];
    std::cout << "\n";
}


void PrtVec(const Moment2State::vector_t &a)
{
    for (size_t k = 0; k < a.size(); k++)
        std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << a[k];
    std::cout << "\n";
}


void PrtMat(const value_mat &M)
{
    for (size_t j = 0; j < M.size1(); j++) {
        for (size_t k = 0; k < M.size2(); k++)
            std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << M(j, k);
        std::cout << "\n";
    }
}


void rd_data(std::fstream &inf, const int column_no, std::vector<double> &Z, std::vector<double> &EM)
{
    std::string       line;
    std::stringstream str;
    double            c[8];

    inf.clear();
    inf.seekg(0, inf.beg);

    while (getline(inf, line) && !inf.fail()) {
        if (line[0] == '%') {
            // Comment.
        } else {
            str.str(line);
            str.clear();
            str >> c[0] >> c[1]>> c[2]>> c[3]>> c[4]>> c[5]>> c[6]>> c[7];
            if (false)
                std::cout << std::scientific << std::setprecision(5)
                          << std::setw(13) << c[0] << std::setw(13) << c[column_no-1] << "\n";
            Z.push_back(c[0]), EM.push_back(c[column_no-1]);
        }
    }
}


void prt_data(std::vector<double> &Z, std::vector<double> &EM)
{
    int k, n;

    n = Z.size();
    std::cout << "\nn = " << n << "\n";
    for (k = 0; k < n; k++)
        std::cout << std::fixed << std::setprecision(5)
                  << std::setw(10) << Z[k] << std::setw(15) << EM[k] << "\n";
}


int get_column(const std::string &str)
{
    int column;

    if (str == "CaviMlp_EFocus1")
        column = 3;
    else if (str == "CaviMlp_EFocus2")
        column = 4;
    else if (str == "CaviMlp_EDipole")
        column = 2;
    else if (str == "CaviMlp_EQuad")
        column = 5;
    else if (str == "CaviMlp_HMono")
        column = 7;
    else if (str == "CaviMlp_HDipole")
        column = 6;
    else if (str == "CaviMlp_HQuad")
        column = 8;
    else {
        std::cout << "get_column: undef. column: " << str << "\n";
        exit(1);
    }

    return column;
}


void calTransfac(std::fstream &inf,  int column_no, const double IonK, const bool half,
                 double &Ecenter, double &T, double &Tp, double &S, double &Sp, double &V0)
{
    // Compute electric center, amplitude, and transit time factors [T, Tp, S, Sp] for RF cavity mode.
    int                 n, k;
    double              dz, eml, em_mom;
    std::vector<double> z, EM;

    rd_data(inf, column_no, z, EM);

    n = z.size();
    if (half) n = (int)round((n-1)/2e0);

    dz = (z[n-1]-z[0])/(n-1);

//    prt_data(z, EM);

    // Start at zero.
    for (k = n-1; k >= 0; k--)
        z[k] -= z[0];

    eml = 0e0, em_mom = 0e0;
    for (k = 0; k < n-1; k++) {
        eml    += (fabs(EM[k])+fabs(EM[k+1]))/2e0*dz;
        em_mom += (z[k]+z[k+1])/2e0*(fabs(EM[k])+fabs(EM[k+1]))/2e0*dz;
    }
    Ecenter = em_mom/eml;

    for (k = 0; k < n; k++)
        z[k] -= Ecenter;

    T = 0;
    for (k = 0; k < n-1; k++)
        T += (EM[k]+EM[k+1])/2e0*cos(IonK*(z[k]+z[k+1])/2e0)*dz;
    T /= eml;

    Tp = 0;
    for (k = 0; k < n-1; k++)
        Tp -= ((z[k]+z[k+1])/2e0)*(EM[k]+EM[k+1])/2e0*sin(IonK*(z[k]+z[k+1])/2e0)*dz;
    Tp /= eml;

    S = 0e0;
    for (k = 0; k < n-1; k++)
        S += (EM[k]+EM[k+1])/2e0*sin(IonK*(z[k]+z[k+1])/2e0)*(z[k+1]-z[k]);
    S /= eml;

    Sp = 0e0;
    for (k = 0; k < n-1; k++) {
        Sp += (z[k]+z[k+1])/2e0*(EM[k]+EM[k+1])/2e0*cos(IonK*(z[k]+z[k+1])/2e0)*dz;
    }
    Sp /= eml;

    V0 = eml/MeVtoeV/MtoMM;
}


double PwrSeries(const double beta,
                 const double a0, const double a1, const double a2, const double a3,
                 const double a4, const double a5, const double a6, const double a7,
                 const double a8, const double a9)
{
    int    k;
    double f;

    const int    n   = 10;
    const double a[] = {a0, a1, a2, a3, a4, a5, a6, a7, a8, a9};

    f = a[0];
    for (k = 1; k < n; k++)
        f += a[k]*pow(beta, k);

    return f;
}


void ElementRFCavity::TransFacts(const int cavilabel, double beta, const double CaviIonK, const int gaplabel, const double EfieldScl,
                                 double &Ecen, double &T, double &Tp, double &S, double &Sp, double &V0)
{
    // Evaluate Electric field center, transit factors [T, T', S, S'] and cavity field.
    std::ostringstream  strm;

    // For debugging of TTF function.
    if (false) {
        calTransfac(inf2, 2, CaviIonK, true, Ecen, T, Tp, S, Sp, V0);
        V0 *= EfieldScl;
        printf("\n%19.12e %19.12e %19.12e %19.12e %19.12e %19.12e\n", Ecen, T, Tp, S, Sp, V0);
//        return;
    }

    switch (cavilabel) {
    case 41:
        if (beta < 0.025 || beta > 0.08) {
            std::cout << "*** TransFacts: CaviIonK out of Range " << cavilabel << "\n";
            calTransfac(inf2, 2, CaviIonK, true, Ecen, T, Tp, S, Sp, V0);
            V0 *= EfieldScl;
            return;
        }
        switch (gaplabel) {
        case 0:
            // One gap evaluation.
            Ecen = 120.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = PwrSeries(beta, -4.109, 399.9, -1.269e4, 1.991e5, -1.569e6, 4.957e6, 0.0, 0.0, 0.0, 0.0);
            Sp   = PwrSeries(beta, 61.98, -1.073e4, 4.841e5, 9.284e6, 8.379e7, -2.926e8, 0.0, 0.0, 0.0, 0.0);
            V0   = 0.98477*EfieldScl;
            break;
        case 1:
            // Two gap calculation, first gap.
            Ecen = 0.0006384*pow(beta, -1.884) + 86.69;
            T    = PwrSeries(beta, 0.9232, -123.2, 3570, -5.476e4, 4.316e5, -1.377e6, 0.0, 0.0, 0.0, 0.0);
            Tp   = PwrSeries(beta, 1.699, 924.7, -4.062e4, 7.528e5, -6.631e6, 2.277e7, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp   = PwrSeries(beta, -1.571, 25.59, 806.6, -2.98e4, 3.385e5, -1.335e6, 0.0, 0.0, 0.0, 0.0);
            V0   = 0.492385*EfieldScl;
            break;
        case 2:
            // Two gap calculation, second gap.
            Ecen = -0.0006384*pow(beta, -1.884) + 33.31;
            T    = PwrSeries(beta, -0.9232, 123.2, -3570, 5.476e4, -4.316e5, 1.377e6, 0.0, 0.0, 0.0, 0.0);
            Tp   = PwrSeries(beta, -1.699, -924.7, 4.062e4, -7.528e5, 6.631e6, -2.277e7, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp    = PwrSeries(beta, -1.571, 25.59, 806.6, -2.98e4, 3.385e5, -1.335e6, 0.0, 0.0, 0.0, 0.0);
            V0   = 0.492385*EfieldScl;
            break;
        default:
            strm << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            throw std::runtime_error(strm.str());
        }
        break;
    case 85:
        if (beta < 0.05 || beta > 0.25) {
            std::cout << "*** TransFacts: CaviIonK out of Range " << cavilabel << "\n";
            calTransfac(inf2, 2, CaviIonK, true, Ecen, T, Tp, S, Sp, V0);
            V0 *= EfieldScl;
            return;
        }
        switch (gaplabel) {
          case 0:
            Ecen = 150.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = PwrSeries(beta, -6.811, 343.9, -6385, 6.477e4, -3.914e5, 1.407e6, -2.781e6, 2.326e6, 0.0, 0.0);
            Sp   = PwrSeries(beta, 162.7, -1.631e4, 4.315e5, -5.344e6, 3.691e7, -1.462e8, 3.109e8, -2.755e8, 0.0, 0.0);
            V0   = 1.967715*EfieldScl;
            break;
        case 1:
            Ecen = 0.0002838*pow(beta, -2.13) + 76.5;
            T    = 0.0009467*pow(beta, -1.855) - 1.002;
            Tp   = PwrSeries(beta, 24.44, -334, 2468, -1.017e4, 2.195e4, -1.928e4, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp   = -0.0009751*pow(beta, -1.898) + 0.001568;
            V0   = 0.9838574*EfieldScl;
            break;
        case 2:
            Ecen = -0.0002838*pow(beta, -2.13) + 73.5;
            T    = -0.0009467*pow(beta, -1.855) + 1.002;
            Tp   = PwrSeries(beta, -24.44, 334, -2468, 1.017e4, -2.195e4, 1.928e4, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp   = -0.0009751*pow(beta, -1.898) + 0.001568;
            V0   = 0.9838574*EfieldScl;
            break;
        default:
            strm << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            throw std::runtime_error(strm.str());
        }
        break;
    case 29:
        if (beta < 0.15 || beta > 0.4) {
            std::cout << "*** TransFacts: CaviIonK out of Range " << cavilabel << "\n";
            calTransfac(inf2, 2, CaviIonK, true, Ecen, T, Tp, S, Sp, V0);
            V0 *= EfieldScl;
            return;
        }
        switch (gaplabel) {
          case 0:
            Ecen = 150.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = PwrSeries(beta, -4.285, 58.08, -248, 486, -405.6, 76.54, 0.0, 0.0, 0.0, 0.0);
            Sp   = PwrSeries(beta, 888, -2.043e4, 1.854e5, -9.127e5, 2.695e6, -4.791e6, 4.751e6, -2.025e6, 0.0, 0.0);
            V0   = 2.485036*EfieldScl;
            break;
        case 1:
            Ecen = 0.01163*pow(beta, -2.001) + 91.77;
            T    = 0.02166*pow(beta, -1.618) - 1.022;
            Tp   = PwrSeries(beta, -11.25, 534.7, -3917, 1.313e4, -2.147e4, 1.389e4, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp   = PwrSeries(beta, -0.8283, -4.409, 78.77, -343.9, 645.1, -454.4, 0.0, 0.0, 0.0, 0.0);
            V0   = 1.242518*EfieldScl;
            break;
        case 2:
            Ecen =-0.01163*pow(beta, -2.001) + 58.23;
            T    =-0.02166*pow(beta, -1.618) + 1.022;
            Tp   = PwrSeries(beta,  11.25, -534.7,  3917, -1.313e4, 2.147e4, -1.389e4, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp   = PwrSeries(beta, -0.8283, -4.409, 78.77, -343.9, 645.1, -454.4, 0.0, 0.0, 0.0, 0.0);
            V0   = 1.242518*EfieldScl;
            break;
        default:
            strm << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            throw std::runtime_error(strm.str());
        }
        break;
    case 53:
        if (beta < 0.3 || beta > 0.6) {
            std::cout << "*** TransFacts: CaviIonK out of Range " << cavilabel << "\n";
            calTransfac(inf2, 2, CaviIonK, true, Ecen, T, Tp, S, Sp, V0);
            V0 *= EfieldScl;
            return;
        }
        switch (gaplabel) {
          case 0:
            Ecen = 250.0; // [mm].
            T    = 0.0;
            Tp   = 0.0;
            S    = PwrSeries(beta, -4.222, 26.64, -38.49, -17.73, 84.12, -52.93, 0.0, 0.0, 0.0, 0.0);
            Sp   = PwrSeries(beta, -1261, -1.413e4, 5.702e4, -1.111e5, 1.075e5, -4.167e4, 0.0, 0.0 , 0.0, 0.0);
            V0   = 4.25756986*EfieldScl;
            break;
        case 1:
            Ecen = 0.01219*pow(beta, -2.348) + 137.8;
            T    = 0.04856*pow(beta, -1.68) - 1.018;
            Tp   = PwrSeries(beta, -3.612, 422.8, -1973, 4081, -4109, 1641, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp   = -0.03969*pow(beta, -1.775) +0.009034;
            V0   = 2.12878493*EfieldScl;
            break;
        case 2:
            Ecen = -0.01219*pow(beta, -2.348) + 112.2;
            T    = -0.04856*pow(beta, -1.68) + 1.018;
            Tp   = PwrSeries(beta, 3.612, -422.8, 1973, -4081, 4109, -1641, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp   = -0.03969*pow(beta, -1.775) +0.009034;
            V0   = 2.12878493*EfieldScl;
            break;
        default:
            strm << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            throw std::runtime_error(strm.str());
        }
        break;
    default:
        strm << "*** GetTransitFac: undef. cavity type" << "\n";
        throw std::runtime_error(strm.str());
    }

    // Convert from [mm] to [m].
//    Ecen /= MtoMM;

//    printf("%19.12e %19.12e %19.12e %19.12e %19.12e %19.12e\n", Ecen, T, Tp, S, Sp, V0);
}


void ElementRFCavity::TransitFacMultipole(const int cavi, const std::string &flabel, const double CaviIonK,
                                          double &T, double &S)
{
    double Ecen, Tp, Sp, V0;

    // For debugging of TTF function.
    if (false) {
        calTransfac(inf3, get_column(flabel), CaviIonK, false, Ecen, T, Tp, S, Sp, V0);
        return;
    }

    if (((cavi == 1) && (CaviIonK < 0.025 || CaviIonK > 0.055)) ||
        ((cavi == 2) && (CaviIonK < 0.006 || CaviIonK > 0.035)) ||
        ((cavi == 3) && (CaviIonK < 0.01687155 || CaviIonK > 0.0449908)) ||
        ((cavi == 4) && (CaviIonK < 0.0112477 || CaviIonK > 0.0224954))) {
        std::cout << "*** TransitFacMultipole: CaviIonK out of Range" << "\n";
        calTransfac(inf3, get_column(flabel), CaviIonK, false, Ecen, T, Tp, S, Sp, V0);
        return;
    }

    if (flabel == "CaviMlp_EFocus1") {
        switch (cavi) {
        case 1:
            T = PwrSeries(CaviIonK, 1.256386e+02, -3.108322e+04, 3.354464e+06, -2.089452e+08, 8.280687e+09, -2.165867e+11,
                          3.739846e+12, -4.112154e+13, 2.613462e14, -7.316972e14);
            S = PwrSeries(CaviIonK, 1.394183e+02, -3.299673e+04, 3.438044e+06, -2.070369e+08, 7.942886e+09, -2.013750e+11,
                         3.374738e+12, -3.605780e+13, 2.229446e+14, -6.079177e+14);
            break;
        case 2:
            T = PwrSeries(CaviIonK, -9.450041e-01, -3.641390e+01, 9.926186e+03, -1.449193e+06, 1.281752e+08, -7.150297e+09,
                          2.534164e+11, -5.535252e+12, 6.794778e+13, -3.586197e+14);
            S = PwrSeries(CaviIonK, 9.928055e-02, -5.545119e+01, 1.280168e+04, -1.636888e+06, 1.279801e+08, -6.379800e+09,
                          2.036575e+11, -4.029152e+12, 4.496323e+13, -2.161712e+14);
            break;
        case 3:
            T = PwrSeries(CaviIonK, -1.000000e+00, 2.778823e-07, 6.820327e+01, 4.235106e-03, -1.926935e+03, 1.083516e+01,
                          2.996807e+04, 6.108642e+03, -3.864554e+05, 6.094390e+05);
            S = PwrSeries(CaviIonK, -4.530303e-10, 1.625011e-07, -2.583224e-05, 2.478684e+01, -1.431967e-01, -1.545412e+03,
                          -1.569820e+02, 3.856713e+04, -3.159828e+04, -2.700076e+05);

            break;
        case 4:
            T = PwrSeries(CaviIonK, -1.000000e+00, -2.406447e-07, 9.480040e+01, -7.659927e-03, -4.926996e+03, -3.504383e+01,
                          1.712590e+05, -1.964643e+04, -4.142976e+06, 6.085390e+06);
            S = PwrSeries(CaviIonK, 3.958048e-11, -2.496811e-08, 7.027794e-06, -8.662787e+01, 1.246098e-01, 9.462491e+03,
                          4.481784e+02, -4.552412e+05, 3.026543e+05, 8.798256e+06);
            break; 
        }
    } else if (flabel == "CaviMlp_EFocus2") {
        switch (cavi) {
        case 1:
            T = PwrSeries(CaviIonK, 1.038803e+00, -9.121320e+00, 8.943931e+02, -5.619149e+04, 2.132552e+06, -5.330725e+07,
                          8.799404e+08, -9.246033e+09, 5.612073e+10, -1.499544e+11);
            S = PwrSeries(CaviIonK, 1.305154e-02, -2.585211e+00, 2.696971e+02, -1.488249e+04, 5.095765e+05, -1.154148e+07,
                          1.714580e+08, -1.604935e+09, 8.570757e+09, -1.983302e+10);
            break;
        case 2:
            T = PwrSeries(CaviIonK, 9.989307e-01, 7.299233e-01, -2.932580e+02, 3.052166e+04, -2.753614e+06, 1.570331e+08,
                          -5.677804e+09, 1.265012e+11, -1.584238e+12, 8.533351e+12);
            S = PwrSeries(CaviIonK, -3.040839e-03, 2.016667e+00, -4.313590e+02, 5.855139e+04, -4.873584e+06, 2.605444e+08,
                          -8.968899e+09, 1.923697e+11, -2.339920e+12, 1.233014e+13);
            break;
        case 3:
            T = PwrSeries(CaviIonK, 1.000000e+00, -4.410575e-06, -8.884752e+01, -7.927594e-02, 4.663277e+03, -2.515405e+02,
                         -1.797134e+05, -1.904305e+05, 8.999378e+06, -2.951362e+07);
            S = PwrSeries(CaviIonK, 6.387813e-08, -2.300899e-05, 3.676251e-03, -1.703282e+02, 2.066461e+01, 1.704569e+04,
                          2.316653e+04, -1.328926e+06, 4.853676e+06, 1.132796e+06);

            break;
        case 4:
            T = PwrSeries(CaviIonK, 1.000000e+00, -5.025186e-06, -1.468976e+02, -2.520376e-01,2.048799e+04, -2.224267e+03,
                          -2.532091e+06, -4.613480e+06, 3.611911e+08, -1.891951e+09);
            S = PwrSeries(CaviIonK, -1.801149e-08, 1.123280e-05, -3.126902e-03, 4.655245e+02, -5.431878e+01, -1.477730e+05,
                          -1.922110e+05, 2.795761e+07, -1.290046e+08, -4.656951e+08);
            break;
        }
    } else if (flabel == "CaviMlp_EDipole") {
        switch (cavi) {
        case 1:
            T = PwrSeries(CaviIonK, -1.005885e+00, 1.526489e+00, -1.047651e+02, 1.125013e+04, -4.669147e+05, 1.255841e+07,
                          -2.237287e+08, 2.535541e+09, -1.656906e+10, 4.758398e+10);
            S = PwrSeries(CaviIonK, -2.586200e-02, 5.884367e+00, -6.407538e+02, 3.888964e+04, -1.488484e+06, 3.782592e+07,
                          -6.361033e+08, 6.817810e+09, -4.227114e+10, 1.155597e+11);
            break;
        case 2:
            T = PwrSeries(CaviIonK, -9.999028e-01, -6.783669e-02, 1.415756e+02, -2.950990e+03, 2.640980e+05, -1.570742e+07,
                          5.770450e+08, -1.303686e+10, 1.654958e+11, -9.030017e+11);
            S = PwrSeries(CaviIonK, 2.108581e-04, -3.700608e-01, 2.851611e+01, -3.502994e+03, 2.983061e+05, -1.522679e+07,
                          4.958029e+08, -1.002040e+10, 1.142835e+11, -5.617061e+11);
            break;
        default:
            std::ostringstream strm;
            strm << "*** 0.29 HWR and 0.53HWR havr no dipole term\n";
            throw std::runtime_error(strm.str());
            break;
        }
    } else if (flabel == "CaviMlp_EQuad") {
        switch (cavi) {
        case 1:
            T = PwrSeries(CaviIonK, 1.038941e+00, -9.238897e+00, 9.127945e+02, -5.779110e+04, 2.206120e+06, -5.544764e+07,
                          9.192347e+08, -9.691159e+09, 5.896915e+10, -1.578312e+11);
            S = PwrSeries(CaviIonK, 1.248096e-01, -2.923507e+01, 3.069331e+03, -1.848380e+05, 7.094882e+06, -1.801113e+08,
                          3.024208e+09, -3.239241e+10, 2.008767e+11, -5.496217e+11);
            break;
        case 2:
            T = PwrSeries(CaviIonK, 1.000003e+00, -1.015639e-03, -1.215634e+02, 1.720764e+01, 3.921401e+03, 2.674841e+05,
                          -1.236263e+07, 3.128128e+08, -4.385795e+09, 2.594631e+10);
            S = PwrSeries(CaviIonK, -1.756250e-05, 2.603597e-01, -2.551122e+00, -4.840638e+01, -2.870201e+04, 1.552398e+06,
                          -5.135200e+07, 1.075958e+09, -1.277425e+10, 6.540748e+10);
            break;
        case 3:
            T = PwrSeries(CaviIonK, 1.000000e+00, 6.239107e-06, -1.697479e+02, 3.444883e-02, 1.225241e+04, -1.663533e+02,
                          -5.526645e+05, -3.593353e+05, 2.749580e+07, -9.689870e+07);
            S = PwrSeries(CaviIonK, 2.128708e-07, -7.985618e-05, 1.240259e-02, -3.211339e+02, 7.098731e+01, 3.474652e+04,
                          8.187145e+04, -3.731688e+06, 1.802053e+07, -1.819958e+07);
            break;
        case 4:
            T = PwrSeries(CaviIonK, 9.998746e-01, -2.431292e-05, -5.019138e+02, -1.176338e+00, 1.006054e+05, -9.908805e+03,
                          -1.148028e+07, -1.922707e+07, 1.432258e+09, -7.054482e+09);
            S = PwrSeries(CaviIonK, 6.003340e-08, -1.482633e-02, 1.037590e-02, -2.235440e+03, 1.790006e+02, 6.456882e+05,
                          6.261020e+05, -1.055477e+08, 4.110502e+08, 2.241301e+09);
            break;
        }
    } else if (flabel == "CaviMlp_HMono") {
        switch (cavi) {
        case 1:
            T = PwrSeries(CaviIonK, 1.703336e+00, -1.671357e+02, 1.697657e+04, -9.843253e+05, 3.518178e+07, -8.043084e+08,
                          1.165760e+10, -1.014721e+11, 4.632851e+11, -7.604796e+11);
            S = PwrSeries(CaviIonK, 1.452657e+01, -3.409550e+03, 3.524921e+05, -2.106663e+07, 8.022856e+08,
                          -2.019481e+10, 3.360597e+11, -3.565836e+12, 2.189668e+13, -5.930241e+13);
            break;
        case 2:
            T = PwrSeries(CaviIonK, 1.003228e+00, -1.783406e+00, 1.765330e+02, -5.326467e+04, 4.242623e+06, -2.139672e+08,
                          6.970488e+09, -1.411958e+11, 1.617248e+12, -8.000662e+12);
            S = PwrSeries(CaviIonK, -1.581533e-03, 1.277444e+00, -2.742508e+02, 3.966879e+04, -3.513478e+06, 1.962939e+08,
                          -6.991916e+09, 1.539708e+11, -1.910236e+12, 1.021016e+13);
            break;
        case 3:
            T = PwrSeries(CaviIonK, 9.999990e-01, 3.477993e-04, -2.717994e+02, 4.554376e+00, 3.083481e+04, 8.441315e+03,
                          -2.439843e+06, 1.322379e+06, 1.501750e+08, -6.822135e+08);
            S = PwrSeries(CaviIonK, 1.709084e-06, -6.240506e-04, 1.013278e-01, -2.649338e+02, 5.944163e+02,
                          4.588900e+04, 7.110518e+05, -2.226574e+07, 1.658524e+08, -3.976459e+08);
            break;
        case 4:
            T = PwrSeries(CaviIonK, 1.000000e+00, -4.358956e-05, -7.923870e+02, -2.472669e+00, 2.241378e+05, -2.539286e+04,
                          -3.385480e+07, -6.375134e+07, 5.652166e+09, -3.355877e+10);
            S = PwrSeries(CaviIonK, 1.163146e-07, -7.302018e-05, 2.048587e-02, -3.689694e+02, 3.632907e+02, 1.757838e+05,
                          1.327057e+06, -9.520645e+07, 9.406709e+08, -2.139562e+09);
            break;
        }
    } else if (flabel == "CaviMlp_HDipole") {
        switch (cavi) {
        case 1:
            T = PwrSeries(CaviIonK, 6.853803e-01, 7.075414e+01, -7.117391e+03, 3.985674e+05, -1.442888e+07, 3.446369e+08,
                          -5.420826e+09, 5.414689e+10, -3.116216e+11, 7.869717e+11);
            S = PwrSeries(CaviIonK, 1.021102e+00, -2.441117e+02, 2.575274e+04, -1.569273e+06, 6.090118e+07, -1.562284e+09,
                          2.649289e+10, -2.864139e+11, 1.791634e+12, -4.941947e+12);
            break;
        case 2:
            T = PwrSeries(CaviIonK, 1.014129e+00, -8.016304e+00, 1.631339e+03, -2.561826e+05, 2.115355e+07, -1.118723e+09,
                          3.821029e+10, -8.140248e+11, 9.839613e+12, -5.154137e+13);
            S = PwrSeries(CaviIonK, -4.688714e-03, 3.299051e+00, -8.101936e+02, 1.163814e+05, -1.017331e+07, 5.607330e+08,
                          -1.967300e+10, 4.261388e+11, -5.194592e+12, 2.725370e+13);
            break;
        default:
            std::ostringstream strm;
            strm << "*** 0.29 HWR and 0.53HWR have no dipole term\n";
            throw std::runtime_error(strm.str());
            break;
        }
    } else if (flabel == "CaviMlp_HQuad") {
        switch (cavi) {
        case 1:
            T = PwrSeries(CaviIonK, -1.997432e+00, 2.439177e+02, -2.613724e+04, 1.627837e+06, -6.429625e+07, 1.676173e+09,
                          -2.885455e+10, 3.163675e+11, -2.005326e+12, 5.600545e+12);
            S = PwrSeries(CaviIonK, -2.470704e+00, 5.862902e+02, -6.135071e+04, 3.711527e+06, -1.431267e+08, 3.649414e+09,
                          -6.153570e+10, 6.617859e+11, -4.119861e+12, 1.131390e+13);
            break;
        case 2:
            T = PwrSeries(CaviIonK, -1.000925e+00, 5.170302e-01, 9.311761e+01, 1.591517e+04, -1.302247e+06, 6.647808e+07,
                          -2.215417e+09, 4.603390e+10, -5.420873e+11, 2.764042e+12);
            S = PwrSeries(CaviIonK, 3.119419e-04, -4.540868e-01, 5.433028e+01, -7.571946e+03, 6.792565e+05, -3.728390e+07,
                          1.299263e+09, -2.793705e+10, 3.377097e+11, -1.755126e+12);
            break;
        case 3:
            T = PwrSeries(CaviIonK, -9.999997e-01, -1.049624e-04, 2.445420e+02, -1.288731e+00, -2.401575e+04, -1.972894e+03,
                          1.494708e+06, 2.898145e+05, -8.782506e+07, 3.566907e+08);
            S = PwrSeries(CaviIonK, -7.925695e-07, 2.884963e-04, -4.667266e-02, 2.950936e+02, -2.712131e+02, -4.260259e+04,
                          -3.199682e+05, 1.103376e+07, -7.304474e+07, 1.479036e+08);
            break;
        case 4:
            T = PwrSeries(CaviIonK, -1.000000e+00, 4.357777e-05, 7.605879e+02, 2.285787e+00, -2.009415e+05, 2.149581e+04,
                         2.773856e+07, 4.886782e+07, -4.127019e+09, 2.299278e+10);
            S = PwrSeries(CaviIonK, -1.483304e-07, 9.278457e-05, -2.592071e-02, 1.690272e+03, -4.545599e+02, -6.192487e+05,
                          -1.632321e+06, 1.664856e+08, -1.124066e+09, -3.121299e+08);
            break;            
        }
    } else {
        std::ostringstream strm;
        strm << "*** TransitFacMultipole: undef. multipole type " << flabel << "\n";
        throw std::runtime_error(strm.str());
    }
}


double GetCavPhase(const int cavi, Particle ref, const double IonFys, const double multip)
{
    /* If the cavity is not at full power, the method gives synchrotron
     * phase slightly different from the nominal value.                 */

    double IonEk, Fyc;

    IonEk = (ref.IonW-ref.IonEs)/MeVtoeV;

    switch (cavi) {
    case 1:
        Fyc = 4.394*pow(IonEk, -0.4965) - 4.731;
        break;
    case 2:
        Fyc = 5.428*pow(IonEk, -0.5008) + 1.6;
        break;
    case 3:
        Fyc = 22.35*pow(IonEk, -0.5348) + 2.026;
        break;
    case 4:
        Fyc = 41.43*pow(IonEk, -0.5775) + 2.59839;
        break;
    case 5:
        Fyc = 5.428*pow(IonEk, -0.5008) + 1.6;
        break;
    default:
        std::ostringstream strm;
        strm << "*** GetCavPhase: undef. cavity type" << "\n";
        throw std::runtime_error(strm.str());
    }

    return IonFys - Fyc - ref.phis*multip;
}


void EvalGapModel(const double dis, const double IonW0, Particle &real, const double IonFy0,
                  const double k, const double Lambda, const double Ecen,
                  const double T, const double S, const double Tp, const double Sp, const double V0,
                  double &IonW_f, double &IonFy_f)
{
    double Iongamma_f, IonBeta_f, k_f;

    IonW_f     = IonW0 + real.IonZ*V0*T*cos(IonFy0+k*Ecen)*MeVtoeV - real.IonZ*V0*S*sin(IonFy0+k*Ecen)*MeVtoeV;
    Iongamma_f = IonW_f/real.IonEs;
    IonBeta_f  = sqrt(1e0-1e0/sqr(Iongamma_f));
    k_f        = 2e0*M_PI/(IonBeta_f*Lambda);

    IonFy_f = IonFy0 + k*Ecen + k_f*(dis-Ecen)
              + real.IonZ*V0*k*(Tp*sin(IonFy0+k*Ecen)+Sp*cos(IonFy0+k*Ecen))/(2e0*(IonW0-real.IonEs)/MeVtoeV);
}

int get_MpoleLevel(const Config &conf)
{
    int MpoleLevel = 0;

    std::string str = conf.get<std::string>("MpoleLevel", "2");
    if (str == "0")
        MpoleLevel = 0;
    else if (str == "1")
        MpoleLevel = 1;
    else if (str == "2")
        MpoleLevel = 2;
    else {
        std::cout << "get_MpoleLevel: undef. MpoleLevel " << MpoleLevel << "\n";
        exit(1);
    }

    return MpoleLevel;
}


void  ElementRFCavity::GetCavMatParams(const int cavi, const double beta_tab[], const double gamma_tab[], const double CaviIonK[])
{
    // Evaluate time transit factors and acceleration.

    std::string line, Elem, Name;
    int         MpoleLevel;
    double      s, Length, Aper, E0, T, S, Accel;
    std::stringstream str;

    MpoleLevel = get_MpoleLevel(conf());

    inf1.clear();
    inf1.seekg(0, inf1.beg);

    CavTLMLineTab.clear();

    s = CavData.s[0];
    while (getline(inf1, line) && !inf1.fail()) {
        T = 0e0, S = 0e0, Accel = 0e0;
        if (line[0] == '%') {
            // Comment.
        } else {
            str.str(line);
            str.clear();
            str >> Elem >> Name >> Length >> Aper;

            s += Length;

            if ((Elem != "drift") && (Elem != "AccGap"))
                str >> E0;
            else
                E0 = 0e0;

            if (Elem == "drift") {
            } else if (Elem == "EFocus1") {
                if (s < 0e0) {
                    // First gap. By reflection 1st Gap EFocus1 is 2nd gap EFocus2.
                    ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_EFocus2", CaviIonK[0], T, S);
                    // First gap *1, transverse E field the same.
                    S = -S;
                } else {
                    // Second gap.
                    ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_EFocus1", CaviIonK[1], T, S);
                }
            } else if (Elem == "EFocus2") {
                if (s < 0e0) {
                    // First gap.
                    ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_EFocus1", CaviIonK[0], T, S);
                    S = -S;
                } else {
                    // Second gap.
                    ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_EFocus2", CaviIonK[1], T, S);
                }
            } else if (Elem == "EDipole") {
                if (MpoleLevel >= 1) {
                    if (s < 0e0) {
                        ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_EDipole", CaviIonK[0], T, S);
                        // First gap *1, transverse E field the same.
                        S = -S;
                    } else {
                        // Second gap.
                        ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_EDipole", CaviIonK[1], T, S);
                    }
                }
            } else if (Elem == "EQuad") {
                if (MpoleLevel >= 2) {
                    if (s < 0e0) {
                        // First gap.
                        ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_EQuad", CaviIonK[0], T, S);
                        S = -S;
                    } else {
                        // Second Gap
                        ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_EQuad", CaviIonK[1], T, S);
                    }
                }
            } else if (Elem == "HMono") {
                if (MpoleLevel >= 2) {
                    if (s < 0e0) {
                        // First gap.
                        ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_HMono", CaviIonK[0], T, S);
                        T = -T;
                    } else {
                        // Second Gap
                        ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_HMono", CaviIonK[1], T, S);
                    }
                }
            } else if (Elem == "HDipole") {
                if (MpoleLevel >= 1) {
                    if (s < 0e0) {
                        // First gap.
                        ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_HDipole", CaviIonK[0], T, S);
                        T = -T;
                    }  else {
                        // Second gap.
                        ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_HDipole", CaviIonK[1], T, S);
                    }
                }
            } else if (Elem == "HQuad") {
                if (MpoleLevel >= 2) {
                    if (s < 0e0) {
                        // First gap.
                        ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_HQuad", CaviIonK[0], T, S);
                        T = -T;
                    } else {
                        // Second gap.
                        ElementRFCavity::TransitFacMultipole(cavi, "CaviMlp_HQuad", CaviIonK[1], T, S);
                    }
                }
            } else if (Elem == "AccGap") {
                if (s < 0e0) {
                    // First gap.
                    Accel = (beta_tab[0]*gamma_tab[0])/((beta_tab[1]*gamma_tab[1]));
                } else {
                    // Second gap.
                    Accel = (beta_tab[1]*gamma_tab[1])/((beta_tab[2]*gamma_tab[2]));
                }
            } else {
                std::ostringstream strm;
                strm << "*** GetCavMatParams: undef. multipole element " << Elem << "\n";
                throw std::runtime_error(strm.str());
            }

            CavTLMLineTab.set(s, Elem, E0, T, S, Accel);
        }
    }

    if (false) {
        std::cout << "\n";
        CavTLMLineTab.show(std::cout);
    }
}


void ElementRFCavity::GenCavMat(const int cavi, const double dis, const double EfieldScl, const double TTF_tab[],
                                const double beta_tab[], const double gamma_tab[], const double Lambda,
                                Particle &real, const double IonFys[], const double Rm, value_mat &M)
{
    /* RF cavity model, transverse only defocusing.
     * 2-gap matrix model.                                            */

    std::string       line, Elem, Name;
    int               MpoleLevel, seg, n;
    double            Length, Aper, Efield, s, k_s[3];
    double            Ecens[2], Ts[2], Ss[2], V0s[2], ks[2], L1, L2, L3;
    double            beta, gamma, kfac, V0, T, S, kfdx, kfdy, dpy, Accel, IonFy;
    value_mat         Idmat, Mlon_L1, Mlon_K1, Mlon_L2;
    value_mat         Mlon_K2, Mlon_L3, Mlon, Mtrans, Mprob;
    std::stringstream str;

    const double IonA = 1e0;

    using boost::numeric::ublas::prod;

    Idmat = boost::numeric::ublas::identity_matrix<double>(PS_Dim);

    MpoleLevel = get_MpoleLevel(conf());

    inf1.clear();
    inf1.seekg(0, inf1.beg);

    k_s[0] = 2e0*M_PI/(beta_tab[0]*Lambda);
    k_s[1] = 2e0*M_PI/(beta_tab[1]*Lambda);
    k_s[2] = 2e0*M_PI/(beta_tab[2]*Lambda);

    // Longitudinal model:  Drift-Kick-Drift, dis: total lenghth centered at 0,
    // Ecens[0] & Ecens[1]: Electric Center position where accel kick applies, Ecens[0] < 0
    // TTFtab:              2*6 vector, Ecens, T Tp S Sp, V0;

    Ecens[0] = TTF_tab[0];
    Ts[0]    = TTF_tab[1];
    Ss[0]    = TTF_tab[3];
    V0s[0]   = TTF_tab[5];
    ks[0]    = 0.5*(k_s[0]+k_s[1]);
    L1       = dis + Ecens[0];       //try change dis/2 to dis 14/12/12

    Mlon_L1 = Idmat;
    Mlon_K1 = Idmat;
    // Pay attention, original is -
    Mlon_L1(4, 5) = -2e0*M_PI/Lambda*(1e0/cube(beta_tab[0]*gamma_tab[0])*MeVtoeV/real.IonEs*L1);
    // Pay attention, original is -k1-k2
    Mlon_K1(5, 4) = -real.IonZ*V0s[0]*Ts[0]*sin(IonFys[0]+ks[0]*L1)-real.IonZ*V0s[0]*Ss[0]*cos(IonFys[0]+ks[0]*L1);

    Ecens[1] = TTF_tab[6];
    Ts[1]    = TTF_tab[7];
    Ss[1]    = TTF_tab[9];
    V0s[1]   = TTF_tab[11];
    ks[1]    = 0.5*(k_s[1]+k_s[2]);
    L2       = Ecens[1] - Ecens[0];

    Mlon_L2 = Idmat;
    Mlon_K2 = Idmat;

    Mlon_L2(4, 5) = -2e0*M_PI/Lambda*(1e0/cube(beta_tab[1]*gamma_tab[1])*MeVtoeV/real.IonEs*L2); //Problem is Here!!
    Mlon_K2(5, 4) = -real.IonZ*V0s[1]*Ts[1]*sin(IonFys[1]+ks[1]*Ecens[1])-real.IonZ*V0s[1]*Ss[1]*cos(IonFys[1]+ks[1]*Ecens[1]);

    L3 = dis - Ecens[1]; //try change dis/2 to dis 14/12/12

    Mlon_L3       = Idmat;
    Mlon_L3(4, 5) = -2e0*M_PI/Lambda*(1e0/cube(beta_tab[2]*gamma_tab[2])*MeVtoeV/real.IonEs*L3);

    Mlon = Idmat;
    Mlon = prod(Mlon_K1, Mlon_L1);
    Mlon = prod(Mlon_L2, Mlon);
    Mlon = prod(Mlon_K2, Mlon);
    Mlon = prod(Mlon_L3, Mlon);

    // Transverse model
    // Drift-FD-Drift-LongiKick-Drift-FD-Drift-0-Drift-FD-Drift-LongiKick-Drift-FD-Drift

    seg    = 0;

    Mtrans = Idmat;
    Mprob  = Idmat;

    beta   = beta_tab[0];
    gamma  = gamma_tab[0];
    IonFy  = IonFys[0];
    kfac   = k_s[0];

    V0 = 0e0, T = 0e0, S = 0e0, kfdx = 0e0, kfdy = 0e0, dpy = 0e0;

    s = CavData.s[0];
    n = 0;
    while (getline(inf1, line) && !inf1.fail()) {
        if (line[0] == '%') {
            // Comment.
        } else {
            n++;
            str.str(line);
            str.clear();
            str >> Elem >> Name >> Length >> Aper;

            s += Length;

            if ((Elem != "drift") && (Elem != "AccGap"))
                str >> Efield;
            else
                Efield = 0e0;

            if (false)
                printf("%9.5f %8s %8s %9.5f %9.5f %9.5f\n",
                       s, Elem.c_str(), Name.c_str(), Length, Aper, Efield);

            Mprob = Idmat;
            if (Elem == "drift") {
                IonFy = IonFy + kfac*Length;

                Mprob(0, 1) = Length;
                Mprob(2, 3) = Length;
                Mtrans      = prod(Mprob, Mtrans);
            } else if (Elem == "EFocus1") {
                V0   = CavTLMLineTab.E0[n-1]*EfieldScl;
                T    = CavTLMLineTab.T[n-1];
                S    = CavTLMLineTab.S[n-1];
                kfdx = real.IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy))/Rm;
                kfdy = kfdx;

                Mprob(1, 0) = kfdx;
                Mprob(3, 2) = kfdy;
                Mtrans      = prod(Mprob, Mtrans);
            } else if (Elem == "EFocus2") {
                V0   = CavTLMLineTab.E0[n-1]*EfieldScl;
                T    = CavTLMLineTab.T[n-1];
                S    = CavTLMLineTab.S[n-1];
                kfdx = real.IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy))/Rm;
                kfdy = kfdx;

                Mprob(1, 0) = kfdx;
                Mprob(3, 2) = kfdy;
                Mtrans      = prod(Mprob, Mtrans);
            } else if (Elem == "EDipole") {
                if (MpoleLevel >= 1) {
                    V0  = CavTLMLineTab.E0[n-1]*EfieldScl;
                    T   = CavTLMLineTab.T[n-1];
                    S   = CavTLMLineTab.S[n-1];
                    dpy = real.IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy));

                    Mprob(3, 6) = dpy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "EQuad") {
                if (MpoleLevel >= 2) {
                    V0   = CavTLMLineTab.E0[n-1]*EfieldScl;
                    T    = CavTLMLineTab.T[n-1];
                    S    = CavTLMLineTab.S[n-1];
                    kfdx =  real.IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy))/Rm;
                    kfdy = -kfdx;

                    Mprob(1, 0) = kfdx;
                    Mprob(3, 2) = kfdy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "HMono") {
                if (MpoleLevel >= 2) {
                    V0   = CavTLMLineTab.E0[n-1]*EfieldScl;
                    T    = CavTLMLineTab.T[n-1];
                    S    = CavTLMLineTab.S[n-1];
                    kfdx = -MU0*C0*real.IonZ*V0/beta/gamma/IonA/AU*(T*cos(IonFy+M_PI/2e0)-S*sin(IonFy+M_PI/2e0))/Rm;
                    kfdy = kfdx;

                    Mprob(1, 0) = kfdx;
                    Mprob(3, 2) = kfdy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "HDipole") {
                if (MpoleLevel >= 1) {
                    V0  = CavTLMLineTab.E0[n-1]*EfieldScl;
                    T   = CavTLMLineTab.T[n-1];
                    S   = CavTLMLineTab.S[n-1];
                    dpy = -MU0*C0*real.IonZ*V0/beta/gamma/IonA/AU*(T*cos(IonFy+M_PI/2e0)-S*sin(IonFy+M_PI/2e0));

                    Mprob(3, 6) = dpy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "HQuad") {
                if (MpoleLevel >= 2) {
                    if (s < 0e0) {
                        // First gap.
                        beta  = (beta_tab[0]+beta_tab[1])/2e0;
                        gamma = (gamma_tab[0]+gamma_tab[1])/2e0;
                    } else {
                        beta  = (beta_tab[1]+beta_tab[2])/2e0;
                        gamma = (gamma_tab[1]+gamma_tab[2])/2e0;
                    }
                    V0   = CavTLMLineTab.E0[n-1]*EfieldScl;
                    T    = CavTLMLineTab.T[n-1];
                    S    = CavTLMLineTab.S[n-1];
                    kfdx = -MU0*C0*real.IonZ*V0/beta/gamma/IonA/AU*(T*cos(IonFy+M_PI/2e0)-S*sin(IonFy+M_PI/2e0))/Rm;
                    kfdy = -kfdx;

                    Mprob(1, 0) = kfdx;
                    Mprob(3, 2) = kfdy;
                    Mtrans      = prod(Mprob, Mtrans);
                }
            } else if (Elem == "AccGap") {
                //IonFy = IonFy + real.IonZ*V0s[0]*kfac*(TTF_tab[2]*sin(IonFy)
                //        + TTF_tab[4]*cos(IonFy))/2/((gamma-1)*real.IonEs/MeVtoeV); //TTF_tab[2]~Tp
                seg    = seg + 1;
                beta   = beta_tab[seg];
                gamma  = gamma_tab[seg];
                kfac   = 2e0*M_PI/(beta*Lambda);
                Accel  = CavTLMLineTab.Accel[n-1];

                Mprob(1, 1) = Accel;
                Mprob(3, 3) = Accel;
                Mtrans      = prod(Mprob, Mtrans);
            } else {
                std::ostringstream strm;
                strm << "*** GetCavMat: undef. multipole type " << Elem << "\n";
                throw std::runtime_error(strm.str());
            }
//            std::cout << Elem << "\n";
//            PrtMat(Mprob);
        }
    }

    M = Mtrans;

    M(4, 4) = Mlon(4, 4);
    M(4, 5) = Mlon(4, 5);
    M(5, 4) = Mlon(5, 4);
    M(5, 5) = Mlon(5, 5);
}


void ElementRFCavity::GetCavMat(const int cavi, const int cavilabel, const double Rm, Particle &real,
                                const double EfieldScl, const double IonFyi_s,
                                const double IonEk_s, const double fRF, value_mat &M)
{
    int    n;
    double CaviLambda, Ecen[2], T[2], Tp[2], S[2], Sp[2], V0[2];
    double dis, IonW_s[3], IonFy_s[3], gamma_s[3], beta_s[3], CaviIonK_s[3];
    double CaviIonK[2];

    CaviLambda  = C0/fRF*MtoMM;

    IonW_s[0]      = IonEk_s + real.IonEs;
    IonFy_s[0]     = IonFyi_s;
    gamma_s[0]     = IonW_s[0]/real.IonEs;
    beta_s[0]      = sqrt(1e0-1e0/sqr(gamma_s[0]));
    CaviIonK_s[0]  = 2e0*M_PI/(beta_s[0]*CaviLambda);

    n   = CavData.s.size();
    dis = (CavData.s[n-1]-CavData.s[0])/2e0;

    ElementRFCavity::TransFacts(cavilabel, beta_s[0], CaviIonK_s[0], 1, EfieldScl,
                                Ecen[0], T[0], Tp[0], S[0], Sp[0], V0[0]);
    EvalGapModel(dis, IonW_s[0], real, IonFy_s[0], CaviIonK_s[0], CaviLambda,
                 Ecen[0], T[0], S[0], Tp[0], Sp[0], V0[0], IonW_s[1], IonFy_s[1]);
    gamma_s[1]     = IonW_s[1]/real.IonEs;
    beta_s[1]      = sqrt(1e0-1e0/sqr(gamma_s[1]));
    CaviIonK_s[1]  = 2e0*M_PI/(beta_s[1]*CaviLambda);

    ElementRFCavity::TransFacts(cavilabel, beta_s[1], CaviIonK_s[1], 2, EfieldScl,
            Ecen[1], T[1], Tp[1], S[1], Sp[1], V0[1]);
    EvalGapModel(dis, IonW_s[1], real, IonFy_s[1], CaviIonK_s[1], CaviLambda,
                Ecen[1], T[1], S[1], Tp[1], Sp[1], V0[1], IonW_s[2], IonFy_s[2]);
    gamma_s[2]    = IonW_s[2]/real.IonEs;
    beta_s[2]     = sqrt(1e0-1e0/sqr(gamma_s[2]));
    CaviIonK_s[2] = 2e0*M_PI/(beta_s[2]*CaviLambda);

    Ecen[0] = Ecen[0] - dis;

    double TTF_tab[] = {Ecen[0], T[0], Tp[0], S[0], Sp[0], V0[0], Ecen[1], T[1], Tp[1], S[1], Sp[1], V0[1]};
    CaviIonK[0] = (CaviIonK_s[0]+CaviIonK_s[1])/2e0;
    CaviIonK[1] = (CaviIonK_s[1]+CaviIonK_s[2])/2e0;

    if (false) {
        printf("\n GetCavMat:\n");
        printf("CaviIonK: %15.10f %15.10f %15.10f\n", CaviIonK_s[0], CaviIonK_s[1], CaviIonK_s[2]);
        printf("CaviIonK: %15.10f %15.10f\n", CaviIonK[0], CaviIonK[1]);
        printf("beta:     %15.10f %15.10f %15.10f\n", beta_s[0], beta_s[1], beta_s[2]);
        printf("gamma:    %15.10f %15.10f %15.10f\n", gamma_s[0], gamma_s[1], gamma_s[2]);
        printf("Ecen:     %15.10f %15.10f\n", Ecen[0], Ecen[1]);
        printf("T:        %15.10f %15.10f\n", T[0], T[1]);
        printf("Tp:       %15.10f %15.10f\n", Tp[0], Tp[1]);
        printf("S:        %15.10f %15.10f\n", S[0], S[1]);
        printf("Sp:       %15.10f %15.10f\n", Sp[0], Sp[1]);
        printf("V0:       %15.10f %15.10f\n", V0[0], V0[1]);
    }

    ElementRFCavity::GetCavMatParams(cavi, beta_s, gamma_s, CaviIonK);
    ElementRFCavity::GenCavMat(cavi, dis, EfieldScl, TTF_tab, beta_s, gamma_s, CaviLambda, real, IonFy_s, Rm, M);
}


void ElementRFCavity::GetCavBoost(const CavDataType &CavData, Particle &state, const double IonFy0, const double fRF,
                                  const double EfieldScl, double &IonFy, double &accIonW)
{
    int    n = CavData.s.size(),
           k;

    double dis = CavData.s[n-1] - CavData.s[0],
           dz  = dis/(n-1),
           IonW0, CaviLambda, CaviIonK, IonFylast, IonGamma, IonBeta;

    CaviLambda = C0/fRF*MtoMM;

    IonW0 = state.IonW;
    IonFy = IonFy0;
    // Sample rate is different for RF Cavity; due to different RF frequencies.
//    IonK  = state.SampleIonK;
    CaviIonK = 2e0*M_PI*fRF/(state.beta*C0*MtoMM);
    for (k = 0; k < n-1; k++) {
        IonFylast = IonFy;
        IonFy += CaviIonK*dz;
        state.IonW  += state.IonZ*EfieldScl*(CavData.Elong[k]+CavData.Elong[k+1])/2e0
                       *cos((IonFylast+IonFy)/2e0)*dz/MtoMM;
        IonGamma = state.IonW/state.IonEs;
        IonBeta  = sqrt(1e0-1e0/sqr(IonGamma));
        if ((state.IonW-state.IonEs) < 0e0) {
            state.IonW = state.IonEs;
            IonBeta = 0e0;
        }
        CaviIonK = 2e0*M_PI/(IonBeta*CaviLambda);
    }

    accIonW = state.IonW - IonW0;
}


void ElementRFCavity::PropagateLongRFCav(Particle &ref)
{
    std::string CavType;
    int         cavi;
    double      fRF, multip, IonFys, EfieldScl, caviFy, IonFy_i, IonFy_o, accIonW;

    CavType = conf().get<std::string>("cavtype");
    if (CavType == "0.041QWR") {
        cavi = 1;
    } else if (conf().get<std::string>("cavtype") == "0.085QWR") {
        cavi = 2;
    } else if (conf().get<std::string>("cavtype") == "0.29HWR") {
        cavi = 3;
    } else if (conf().get<std::string>("cavtype") == "0.53HWR") {
        cavi = 4;
    } else {
        std::ostringstream strm;
        strm << "*** PropagateLongRFCav: undef. cavity type: " << CavType << "\n";
        throw std::runtime_error(strm.str());
    }

    fRF       = conf().get<double>("f");
    multip    = fRF/SampleFreq;
    IonFys    = conf().get<double>("phi")*M_PI/180e0;  // Synchrotron phase [rad].
    EfieldScl = conf().get<double>("scl_fac");         // Electric field scale factor.

    caviFy = GetCavPhase(cavi, ref, IonFys, multip);

    IonFy_i = multip*ref.phis + caviFy;
    conf().set<double>("phi_ref", caviFy);

    // For the reference particle, evaluate the change of:
    // kinetic energy, absolute phase, beta, and gamma.
    GetCavBoost(CavData, ref, IonFy_i, fRF, EfieldScl, IonFy_o, accIonW);

    ref.IonEk       = ref.IonW - ref.IonEs;
    ref.recalc();
    ref.phis       += (IonFy_o-IonFy_i)/multip;
}


void ElementRFCavity::InitRFCav(Particle &real, double &accIonW,
                                double &avebeta, double &avegamma, value_mat &M)
{
    std::string CavType;
    int         cavi, cavilabel, multip;
    double      Rm, IonFy_i, Ek_i, fRF, EfieldScl, IonFy_o;

    CavType = conf().get<std::string>("cavtype");
    if (CavType == "0.041QWR") {
        cavi       = 1;
        cavilabel  = 41;
        multip     = 1;
        Rm         = 17e0;
    } else if (conf().get<std::string>("cavtype") == "0.085QWR") {
        cavi       = 2;
        cavilabel  = 85;
        multip     = 1;
        Rm         = 17e0;
    } else if (conf().get<std::string>("cavtype") == "0.29HWR") {
        cavi       = 3;
        cavilabel  = 29;
        multip     = 4;
        Rm         = 20e0;
    } else if (conf().get<std::string>("cavtype") == "0.53HWR") {
        cavi       = 4;
        cavilabel  = 53;
        multip     = 4;
        Rm         = 20e0;
    } else if (conf().get<std::string>("cavtype") == "??EL") {
        // 5 Cell elliptical.
        cavi       = 5;
        cavilabel  = 53;
        multip     = 8;
        Rm         = 20e0;
    } else {
        std::ostringstream strm;
        strm << "*** InitRFCav: undef. cavity type: " << CavType << "\n";
        throw std::runtime_error(strm.str());
    }

    IonFy_i   = multip*real.phis + conf().get<double>("phi_ref");
    Ek_i      = real.IonEk;
    real.IonW = real.IonEk + real.IonEs;

    avebeta   = real.beta;
    avegamma  = real.gamma;
    fRF       = conf().get<double>("f");
    EfieldScl = conf().get<double>("scl_fac");         // Electric field scale factor.

    ElementRFCavity::GetCavBoost(CavData, real, IonFy_i, fRF, EfieldScl, IonFy_o, accIonW);

    real.IonEk       = real.IonW - real.IonEs;
    real.recalc();
    real.phis       += (IonFy_o-IonFy_i)/multip;
    avebeta         += real.beta;
    avebeta         /= 2e0;
    avegamma        += real.gamma;
    avegamma        /= 2e0;

    GetCavMat(cavi, cavilabel, Rm, real, EfieldScl, IonFy_i, Ek_i, fRF, M);
}
