
#include <fstream>

#include "scsi/constants.h"
#include "scsi/moment2.h"
#include "scsi/rf_cavity.h"

// RF Cavity beam dynamics functions.

static
void TransitFacMultipole(const int cavi, const std::string &flabel, const double IonK,
                         double &T, double &S);

static
void TransFacts(const int cavilabel, double beta, const int gaplabel, const double EfieldScl,
                double &Ecen, double &T, double &Tp, double &S, double &Sp, double &V0);

static
void EvalGapModel(const double dis, const double IonW0, Particle &real, const double IonFy0,
                  const double k, const double Lambda, const double Ecen,
                  const double T, const double S, const double Tp, const double Sp, const double V0,
                  double &IonW_f, double &IonFy_f);

static
double GetCavPhase(const int cavi, const Particle& ref, const double IonFys, const double multip);


void CavDataType::RdData(const std::string &FileName)
{
    std::string       line;
    double            s, Elong;
    std::stringstream str;
    std::fstream      inf;

    inf.open(FileName.c_str(), std::ifstream::in);
    if (!inf.is_open()) {
        std::ostringstream strm;
        strm << "*** RdData: failed to open " << FileName << "\n";
        throw std::runtime_error(strm.str());
    }
    while (getline(inf, line) && !inf.fail()) {
        str.str(line);
        str.clear();
        str >> s >> Elong;
        this->s.push_back(s), this->Elong.push_back(Elong);
        // Convert from [mm] to [m].
//        this->s[this->s.size()-1] /= MtoMM;
    }
    inf.close();

    if (false) {
        std::cout << "\n";
        for (size_t k = 0; k < this->s.size(); k++)
            this->show(std::cout, k);
    }
}


void CavDataType::show(std::ostream& strm, const int k) const
{
    strm << std::scientific << std::setprecision(5)
         << std::setw(13) << this->s[k] << std::setw(13) << this->Elong[k] << "\n";
}


void CavDataType::show(std::ostream& strm) const
{
    for (unsigned int k = 0; k < this->s.size(); k++)
        this->show(strm, k);
}


void CavTLMLineType::clear(void)
{
    this->s.clear(); this->Elem.clear(); this->E0.clear();
    this->T.clear(); this->S.clear(); this->Accel.clear();
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
         << std::setw(9) << this->s[k] << std::setw(10) << this->Elem[k]
         << std::setw(9) << this->T[k] << std::setw(9) << this->S[k]
         << std::setw(9) << this->Accel[k] << "\n";
}


void CavTLMLineType::show(std::ostream& strm) const
{
    for (unsigned int k = 0; k < this->s.size(); k++)
        this->show(strm, k);
}


void PrtVec(const std::vector<double> &a)
{
    for (size_t k = 0; k < a.size(); k++)
        std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << a[k];
    std::cout << "\n";
}

static
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


static
void TransFacts(const int cavilabel, double beta, const int gaplabel, const double EfieldScl,
                double &Ecen, double &T, double &Tp, double &S, double &Sp, double &V0)
{
    // Evaluate Electric field center, transit factors [T, T', S, S'] and cavity field.
    std::vector<double> vec;

    switch (cavilabel) {
    case 41:
        if (beta < 0.025 || beta > 0.08) {
            std::ostringstream strm;
            strm << "*** GetTransitFac: beta out of Range " << beta << "\n";
            throw std::runtime_error(strm.str());
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
            std::ostringstream strm;
            strm << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            throw std::runtime_error(strm.str());
        }
        break;
    case 85:
        if (beta < 0.05 || beta > 0.25) {
            std::ostringstream strm;
            strm << "*** GetTransitFac: beta out of range " << beta << "\n";
            throw std::runtime_error(strm.str());
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
            Tp   = PwrSeries(beta,  24.44, 334,  2468, 1.017e4, -2.195e4, 1.928e4, 0.0, 0.0, 0.0, 0.0);
            S    = 0.0;
            Sp   = -0.0009751*pow(beta, -1.898) + 0.001568;
            V0   = 0.9838574*EfieldScl;
            break;
        default:
            std::ostringstream strm;
            strm << "*** GetTransitFac: undef. number of gaps " << gaplabel << "\n";
            throw std::runtime_error(strm.str());
        }
        break;
    default:
        std::ostringstream strm;
        strm << "*** GetTransitFac: undef. cavity type" << "\n";
        throw std::runtime_error(strm.str());
    }

    // Convert from [mm] to [m].
//    Ecen /= MtoMM;
}


static
void TransitFacMultipole(const int cavi, const std::string &flabel, const double IonK,
                         double &T, double &S)
{

    if ((cavi == 1) && (IonK < 0.025 || IonK > 0.055)) {
        std::ostringstream strm;
        strm << "*** TransitFacMultipole: IonK out of Range 1 " << cavi << " " << IonK << "\n";
        throw std::runtime_error(strm.str());
    } else if ((cavi == 2) && (IonK < 0.006 || IonK > 0.035)) {
        std::ostringstream strm;
        strm << "*** TransitFacMultipole: IonK out of Range 2" << cavi << " " << IonK << "\n";
        throw std::runtime_error(strm.str());
    }

    if (flabel == "CaviMlp_EFocus1") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 1.256386e+02, -3.108322e+04, 3.354464e+06, -2.089452e+08, 8.280687e+09, -2.165867e+11,
                          3.739846e+12, -4.112154e+13, 2.613462e14, -7.316972e14);
            S = PwrSeries(IonK, 1.394183e+02, -3.299673e+04, 3.438044e+06, -2.070369e+08, 7.942886e+09, -2.013750e+11,
                         3.374738e+12, -3.605780e+13, 2.229446e+14, -6.079177e+14);
            break;
        case 2:
            T = PwrSeries(IonK, -9.450041e-01, -3.641390e+01, 9.926186e+03, -1.449193e+06, 1.281752e+08, -7.150297e+09,
                          2.534164e+11, -5.535252e+12, 6.794778e+13, -3.586197e+14);
            S = PwrSeries(IonK, 9.928055e-02, -5.545119e+01, 1.280168e+04, -1.636888e+06, 1.279801e+08, -6.379800e+09,
                          2.036575e+11, -4.029152e+12, 4.496323e+13, -2.161712e+14);
            break;
        }
    } else if (flabel == "CaviMlp_EFocus2") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 1.038803e+00, -9.121320e+00, 8.943931e+02, -5.619149e+04, 2.132552e+06, -5.330725e+07,
                          8.799404e+08, -9.246033e+09, 5.612073e+10, -1.499544e+11);
            S = PwrSeries(IonK, 1.305154e-02, -2.585211e+00, 2.696971e+02, -1.488249e+04, 5.095765e+05, -1.154148e+07,
                          1.714580e+08, -1.604935e+09, 8.570757e+09, -1.983302e+10);
            break;
        case 2:
            T = PwrSeries(IonK, 9.989307e-01, 7.299233e-01, -2.932580e+02, 3.052166e+04, -2.753614e+06, 1.570331e+08,
                          -5.677804e+09, 1.265012e+11, -1.584238e+12, 8.533351e+12);
            S = PwrSeries(IonK, -3.040839e-03, 2.016667e+00, -4.313590e+02, 5.855139e+04, -4.873584e+06, 2.605444e+08,
                          -8.968899e+09, 1.923697e+11, -2.339920e+12, 1.233014e+13);
            break;
        }
    } else if (flabel == "CaviMlp_EDipole") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, -1.005885e+00, 1.526489e+00, -1.047651e+02, 1.125013e+04, -4.669147e+05, 1.255841e+07,
                          -2.237287e+08, 2.535541e+09, -1.656906e+10, 4.758398e+10);
            S = PwrSeries(IonK, -2.586200e-02, 5.884367e+00, -6.407538e+02, 3.888964e+04, -1.488484e+06, 3.782592e+07,
                          -6.361033e+08, 6.817810e+09, -4.227114e+10, 1.155597e+11);
            break;
        case 2:
            T = PwrSeries(IonK, -9.999028e-01, -6.783669e-02, 1.415756e+02, -2.950990e+03, 2.640980e+05, -1.570742e+07,
                          5.770450e+08, -1.303686e+10, 1.654958e+11, -9.030017e+11);
            S = PwrSeries(IonK, 2.108581e-04, -3.700608e-01, 2.851611e+01, -3.502994e+03, 2.983061e+05, -1.522679e+07,
                          4.958029e+08, -1.002040e+10, 1.142835e+11, -5.617061e+11);
            break;
        }
    } else if (flabel == "CaviMlp_EQuad") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 1.038941e+00, -9.238897e+00, 9.127945e+02, -5.779110e+04, 2.206120e+06, -5.544764e+07,
                          9.192347e+08, -9.691159e+09, 5.896915e+10, -1.578312e+11);
            S = PwrSeries(IonK, 1.248096e-01, -2.923507e+01, 3.069331e+03, -1.848380e+05, 7.094882e+06, -1.801113e+08,
                          3.024208e+09, -3.239241e+10, 2.008767e+11, -5.496217e+11);
            break;
        case 2:
            T = PwrSeries(IonK, 1.000003e+00, -1.015639e-03, -1.215634e+02, 1.720764e+01, 3.921401e+03, 2.674841e+05,
                          -1.236263e+07, 3.128128e+08, -4.385795e+09, 2.594631e+10);
            S = PwrSeries(IonK, -1.756250e-05, 2.603597e-01, -2.551122e+00, -4.840638e+01, -2.870201e+04, 1.552398e+06,
                          -5.135200e+07, 1.075958e+09, -1.277425e+10, 6.540748e+10);
            break;
        }
    } else if (flabel == "CaviMlp_HMono") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 1.703336e+00, -1.671357e+02, 1.697657e+04, -9.843253e+05, 3.518178e+07, -8.043084e+08,
                          1.165760e+10, -1.014721e+11, 4.632851e+11, -7.604796e+11);
            S = PwrSeries(IonK, 1.452657e+01, -3.409550e+03, 3.524921e+05, -2.106663e+07, 8.022856e+08,
                          -2.019481e+10, 3.360597e+11, -3.565836e+12, 2.189668e+13, -5.930241e+13);
            break;
        case 2:
            T = PwrSeries(IonK, 1.003228e+00, -1.783406e+00, 1.765330e+02, -5.326467e+04, 4.242623e+06, -2.139672e+08,
                          6.970488e+09, -1.411958e+11, 1.617248e+12, -8.000662e+12);
            S = PwrSeries(IonK, -1.581533e-03, 1.277444e+00, -2.742508e+02, 3.966879e+04, -3.513478e+06, 1.962939e+08,
                          -6.991916e+09, 1.539708e+11, -1.910236e+12, 1.021016e+13);
            break;
        }
    } else if (flabel == "CaviMlp_HDipole") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, 6.853803e-01, 7.075414e+01, -7.117391e+03, 3.985674e+05, -1.442888e+07, 3.446369e+08,
                          -5.420826e+09, 5.414689e+10, -3.116216e+11, 7.869717e+11);
            S = PwrSeries(IonK, 1.021102e+00, -2.441117e+02, 2.575274e+04, -1.569273e+06, 6.090118e+07, -1.562284e+09,
                          2.649289e+10, -2.864139e+11, 1.791634e+12, -4.941947e+12);
            break;
        case 2:
            T = PwrSeries(IonK, 1.014129e+00, -8.016304e+00, 1.631339e+03, -2.561826e+05, 2.115355e+07, -1.118723e+09,
                          3.821029e+10, -8.140248e+11, 9.839613e+12, -5.154137e+13);
            S = PwrSeries(IonK, -4.688714e-03, 3.299051e+00, -8.101936e+02, 1.163814e+05, -1.017331e+07, 5.607330e+08,
                          -1.967300e+10, 4.261388e+11, -5.194592e+12, 2.725370e+13);
            break;
        }
    } else if (flabel == "CaviMlp_HQuad") {
        switch (cavi) {
        case 1:
            T = PwrSeries(IonK, -1.997432e+00, 2.439177e+02, -2.613724e+04, 1.627837e+06, -6.429625e+07, 1.676173e+09,
                          -2.885455e+10, 3.163675e+11, -2.005326e+12, 5.600545e+12);
            S = PwrSeries(IonK, -2.470704e+00, 5.862902e+02, -6.135071e+04, 3.711527e+06, -1.431267e+08, 3.649414e+09,
                          -6.153570e+10, 6.617859e+11, -4.119861e+12, 1.131390e+13);
            break;
        case 2:
            T = PwrSeries(IonK, -1.000925e+00, 5.170302e-01, 9.311761e+01, 1.591517e+04, -1.302247e+06, 6.647808e+07,
                          -2.215417e+09, 4.603390e+10, -5.420873e+11, 2.764042e+12);
            S = PwrSeries(IonK, 3.119419e-04, -4.540868e-01, 5.433028e+01, -7.571946e+03, 6.792565e+05, -3.728390e+07,
                          1.299263e+09, -2.793705e+10, 3.377097e+11, -1.755126e+12);
            break;
        }
    } else {
        std::ostringstream strm;
        strm << "*** TransitFacMultipole: undef. multipole type " << flabel << "\n";
        throw std::runtime_error(strm.str());
    }
}


static
double GetCavPhase(const int cavi, const Particle& ref, const double IonFys, const double multip)
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


static
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

ElementRFCavity::ElementRFCavity(const Config& c)
    :base_t(c)
    ,phi_ref(std::numeric_limits<double>::quiet_NaN())
{
    std::string cav_type = c.get<std::string>("cavtype");
    double L             = c.get<double>("L")*MtoMM;         // Convert from [m] to [mm].

    std::string CavType      = conf().get<std::string>("cavtype");
    std::string cavfile(c.get<std::string>("Eng_Data_Dir", ""));

    if (CavType == "0.041QWR") {
        CavData.RdData(cavfile+"/axisData_41.txt");
        cavfile += "/Multipole41/thinlenlon_41.txt";
    } else if (conf().get<std::string>("cavtype") == "0.085QWR") {
        CavData.RdData(cavfile+"/axisData_85.txt");
        cavfile += "/Multipole85/thinlenlon_85.txt";
    } else {
        std::ostringstream strm;
        strm << "*** InitRFCav: undef. cavity type: " << CavType << "\n";
        throw std::runtime_error(strm.str());
    }

    {
        std::ifstream fstrm(cavfile.c_str());

        std::string rawline;
        unsigned line=0;
        while(std::getline(fstrm, rawline)) {
            line++;

            size_t cpos = rawline.find_first_not_of(" \t");
            if(cpos==rawline.npos || rawline[cpos]=='%')
                continue; // skip blank and comment lines

            cpos = rawline.find_last_not_of("\r\n");
            if(cpos!=rawline.npos)
                rawline = rawline.substr(0, cpos+1);

            std::istringstream lstrm(rawline);
            RawParams params;

            lstrm >> params.type >> params.name >> params.length >> params.aperature;
            bool needE0 = params.type!="drift" && params.type!="AccGap";

            if(needE0)
                lstrm >> params.E0;
            else
                params.E0 = 0.0;

            if(lstrm.fail() && !lstrm.eof()) {
                std::ostringstream strm;
                strm<<"Error parsing line "<<line<<" in "<<cavfile;
                throw std::runtime_error(strm.str());
            }

            lattice.push_back(params);
        }

        if(fstrm.fail() && !fstrm.eof()) {
            std::ostringstream strm;
            strm<<"Error, extra chars at end of file (line "<<line<<") in "<<cavfile;
            throw std::runtime_error(strm.str());
        }
    }

    this->transfer_raw(state_t::PS_X, state_t::PS_PX) = L;
    this->transfer_raw(state_t::PS_Y, state_t::PS_PY) = L;
    // For total path length.
//        this->transfer(state_t::PS_S, state_t::PS_S)  = L;
}

void ElementRFCavity::GetCavMatParams(const int cavi, const double beta_tab[], const double gamma_tab[], const double IonK[],
                                      CavTLMLineType& lineref) const
{
    // Evaluate time transit factors and acceleration.

    lineref.clear();
//    std::cout<<"Recompute CavTLMLineTab from "<<__FUNCTION__<<"\n";

    double s = CavData.s[0];
    for(size_t i=0; i<lattice.size(); i++) {
        const RawParams& P = lattice[i];
//        std::cout<<" > "<<P.type<<" "<<P.name<<" "<<P.length<<" "<<P.aperature<<" "<<P.E0<<"\n";

        double T = 0e0, S = 0e0, Accel = 0e0;

        s += P.length;

        if (P.type == "drift") {
        } else if (P.type == "EFocus1") {
            if (s < 0e0) {
                // First gap. By reflection 1st Gap EFocus1 is 2nd gap EFocus2.
                TransitFacMultipole(cavi, "CaviMlp_EFocus2", IonK[0], T, S);
                // First gap *1, transverse E field the same.
                S = -S;
            } else {
                // Second gap.
                TransitFacMultipole(cavi, "CaviMlp_EFocus1", IonK[1], T, S);
            }
        } else if (P.type == "EFocus2") {
            if (s < 0e0) {
                // First gap.
                TransitFacMultipole(cavi, "CaviMlp_EFocus1", IonK[0], T, S);
                S = -S;
            } else {
                // Second gap.
                TransitFacMultipole(cavi, "CaviMlp_EFocus2", IonK[1], T, S);
            }
        } else if (P.type == "EDipole") {
            if (MpoleLevel >= 1) {
                if (s < 0e0) {
                    TransitFacMultipole(cavi, "CaviMlp_EDipole", IonK[0], T, S);
                    // First gap *1, transverse E field the same.
                    S = -S;
                } else {
                    // Second gap.
                    TransitFacMultipole(cavi, "CaviMlp_EDipole", IonK[1], T, S);
                }
            }
        } else if (P.type == "EQuad") {
            if (MpoleLevel >= 2) {
                if (s < 0e0) {
                    // First gap.
                    TransitFacMultipole(cavi, "CaviMlp_EQuad", IonK[0], T, S);
                    S = -S;
                } else {
                    // Second Gap
                    TransitFacMultipole(cavi, "CaviMlp_EQuad", IonK[1], T, S);
                }
            }
        } else if (P.type == "HMono") {
            if (MpoleLevel >= 2) {
                if (s < 0e0) {
                    // First gap.
                    TransitFacMultipole(cavi, "CaviMlp_HMono", IonK[0], T, S);
                    T = -T;
                } else {
                    // Second Gap
                    TransitFacMultipole(cavi, "CaviMlp_HMono", IonK[1], T, S);
                }
            }
        } else if (P.type == "HDipole") {
            if (MpoleLevel >= 1) {
                if (s < 0e0) {
                    // First gap.
                    TransitFacMultipole(cavi, "CaviMlp_HDipole", IonK[0], T, S);
                    T = -T;
                }  else {
                    // Second gap.
                    TransitFacMultipole(cavi, "CaviMlp_HDipole", IonK[1], T, S);
                }
            }
        } else if (P.type == "HQuad") {
            if (MpoleLevel >= 2) {
                if (s < 0e0) {
                    // First gap.
                    TransitFacMultipole(cavi, "CaviMlp_HQuad", IonK[0], T, S);
                    T = -T;
                } else {
                    // Second gap.
                    TransitFacMultipole(cavi, "CaviMlp_HQuad", IonK[1], T, S);
                }
            }
        } else if (P.type == "AccGap") {
            if (s < 0e0) {
                // First gap.
                Accel = (beta_tab[0]*gamma_tab[0])/((beta_tab[1]*gamma_tab[1]));
            } else {
                // Second gap.
                Accel = (beta_tab[1]*gamma_tab[1])/((beta_tab[2]*gamma_tab[2]));
            }
        } else {
            std::ostringstream strm;
            strm << "*** GetCavMatParams: undef. multipole element " << P.type << "\n";
            throw std::runtime_error(strm.str());
        }

        lineref.set(s, P.type, P.E0, T, S, Accel);
//        std::cout<<" < "<<s<<" "<<P.type<<" "<<P.E0<<" "<<T<<" "<<S<<" "<<Accel<<"\n";
    }

    if (false) {
        std::cout << "\n";
        lineref.show(std::cout);
    }
}


void ElementRFCavity::GenCavMat(const int cavi, const double dis, const double EfieldScl, const double TTF_tab[],
                                const double beta_tab[], const double gamma_tab[], const double Lambda,
                                Particle &real, const double IonFys[], const double Rm, state_t::matrix_t &M)
{
    /* RF cavity model, transverse only defocusing.
     * 2-gap matrix model.                                            */

    int               seg;
    double            Efield, k_s[3];
    double            Ecens[2], Ts[2], Ss[2], V0s[2], ks[2], L1, L2, L3;
    double            beta, gamma, kfac, V0, T, S, kfdx, kfdy, dpy, Accel, IonFy;
    state_t::matrix_t         Idmat, Mlon_L1, Mlon_K1, Mlon_L2;
    state_t::matrix_t         Mlon_K2, Mlon_L3, Mlon, Mtrans, Mprob;
    std::stringstream str;

    const double IonA = 1e0;

    using boost::numeric::ublas::prod;

    Idmat = boost::numeric::ublas::identity_matrix<double>(PS_Dim);

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
//    std::cout<<__FUNCTION__<<" Mlon "<<Mlon<<"\n";

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

    double s = CavData.s[0];
    for(size_t n=0; n<lattice.size(); n++) {
        const RawParams& P = lattice[n];

        s += P.length;

        if ((P.type != "drift") && (P.type != "AccGap"))
            str >> Efield;
        else
            Efield = 0e0;

        if (false)
            printf("%9.5f %8s %8s %9.5f %9.5f %9.5f\n",
                   s, P.type.c_str(), P.name.c_str(), P.length, P.aperature, Efield);

        Mprob = Idmat;
        if (P.type == "drift") {
            IonFy = IonFy + kfac*P.length;

            Mprob(0, 1) = P.length;
            Mprob(2, 3) = P.length;
            Mtrans      = prod(Mprob, Mtrans);
        } else if (P.type == "EFocus1") {
            V0   = CavTLMLineTab.E0[n]*EfieldScl;
            T    = CavTLMLineTab.T[n];
            S    = CavTLMLineTab.S[n];
            kfdx = real.IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy))/Rm;
            kfdy = kfdx;
//                std::cout<<" X EFocus1 kfdx="<<kfdx<<"\n";
//                std::cout<<" Y "<<CavTLMLineTab.E0[n]<<" "<<EfieldScl<<" "<<beta
//                         <<" "<<gamma<<" "<<IonFy<<" "<<Rm<<"\n Z "<<T<<" "<<S<<"\n";

            Mprob(1, 0) = kfdx;
            Mprob(3, 2) = kfdy;
            Mtrans      = prod(Mprob, Mtrans);
        } else if (P.type == "EFocus2") {
            V0   = CavTLMLineTab.E0[n]*EfieldScl;
            T    = CavTLMLineTab.T[n];
            S    = CavTLMLineTab.S[n];
            kfdx = real.IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy))/Rm;
            kfdy = kfdx;
//                std::cout<<" X EFocus2 kfdx="<<kfdx<<"\n";

            Mprob(1, 0) = kfdx;
            Mprob(3, 2) = kfdy;
            Mtrans      = prod(Mprob, Mtrans);
        } else if (P.type == "EDipole") {
            if (MpoleLevel >= 1) {
                V0  = CavTLMLineTab.E0[n]*EfieldScl;
                T   = CavTLMLineTab.T[n];
                S   = CavTLMLineTab.S[n];
                dpy = real.IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy));
//                    std::cout<<" X EDipole dpy="<<dpy<<"\n";

                Mprob(3, 6) = dpy;
                Mtrans      = prod(Mprob, Mtrans);
            }
        } else if (P.type == "EQuad") {
            if (MpoleLevel >= 2) {
                V0   = CavTLMLineTab.E0[n]*EfieldScl;
                T    = CavTLMLineTab.T[n];
                S    = CavTLMLineTab.S[n];
                kfdx =  real.IonZ*V0/sqr(beta)/gamma/IonA/AU*(T*cos(IonFy)-S*sin(IonFy))/Rm;
                kfdy = -kfdx;
//                    std::cout<<" X EQuad kfdx="<<kfdx<<"\n";

                Mprob(1, 0) = kfdx;
                Mprob(3, 2) = kfdy;
                Mtrans      = prod(Mprob, Mtrans);
            }
        } else if (P.type == "HMono") {
            if (MpoleLevel >= 2) {
                V0   = CavTLMLineTab.E0[n]*EfieldScl;
                T    = CavTLMLineTab.T[n];
                S    = CavTLMLineTab.S[n];
                kfdx = -MU0*C0*real.IonZ*V0/beta/gamma/IonA/AU*(T*cos(IonFy+M_PI/2e0)-S*sin(IonFy+M_PI/2e0))/Rm;
                kfdy = kfdx;
//                    std::cout<<" X HMono kfdx="<<kfdx<<"\n";

                Mprob(1, 0) = kfdx;
                Mprob(3, 2) = kfdy;
                Mtrans      = prod(Mprob, Mtrans);
            }
        } else if (P.type == "HDipole") {
            if (MpoleLevel >= 1) {
                V0  = CavTLMLineTab.E0[n]*EfieldScl;
                T   = CavTLMLineTab.T[n];
                S   = CavTLMLineTab.S[n];
                dpy = -MU0*C0*real.IonZ*V0/beta/gamma/IonA/AU*(T*cos(IonFy+M_PI/2e0)-S*sin(IonFy+M_PI/2e0));
//                    std::cout<<" X HDipole dpy="<<dpy<<"\n";

                Mprob(3, 6) = dpy;
                Mtrans      = prod(Mprob, Mtrans);
            }
        } else if (P.type == "HQuad") {
            if (MpoleLevel >= 2) {
                if (s < 0e0) {
                    // First gap.
                    beta  = (beta_tab[0]+beta_tab[1])/2e0;
                    gamma = (gamma_tab[0]+gamma_tab[1])/2e0;
                } else {
                    beta  = (beta_tab[1]+beta_tab[2])/2e0;
                    gamma = (gamma_tab[1]+gamma_tab[2])/2e0;
                }
                V0   = CavTLMLineTab.E0[n]*EfieldScl;
                T    = CavTLMLineTab.T[n];
                S    = CavTLMLineTab.S[n];
                kfdx = -MU0*C0*real.IonZ*V0/beta/gamma/IonA/AU*(T*cos(IonFy+M_PI/2e0)-S*sin(IonFy+M_PI/2e0))/Rm;
                kfdy = -kfdx;
//                    std::cout<<" X HQuad kfdx="<<kfdx<<"\n";

                Mprob(1, 0) = kfdx;
                Mprob(3, 2) = kfdy;
                Mtrans      = prod(Mprob, Mtrans);
            }
        } else if (P.type == "AccGap") {
            //IonFy = IonFy + real.IonZ*V0s[0]*kfac*(TTF_tab[2]*sin(IonFy)
            //        + TTF_tab[4]*cos(IonFy))/2/((gamma-1)*real.IonEs/MeVtoeV); //TTF_tab[2]~Tp
            seg    = seg + 1;
            beta   = beta_tab[seg];
            gamma  = gamma_tab[seg];
            kfac   = 2e0*M_PI/(beta*Lambda);
            Accel  = CavTLMLineTab.Accel[n];
//                std::cout<<" X AccGap Accel="<<Accel<<"\n";

            Mprob(1, 1) = Accel;
            Mprob(3, 3) = Accel;
            Mtrans      = prod(Mprob, Mtrans);
        } else {
            std::ostringstream strm;
            strm << "*** GetCavMat: undef. multipole type " << P.type << "\n";
            throw std::runtime_error(strm.str());
        }
        //            std::cout << Elem << "\n";
        //            PrtMat(Mprob);

//            std::cout<<"Elem "<<P.name<<":"<<P.type<<"\n Mtrans "<<Mtrans<<"\nMprob "<<Mprob<<"\n";
    }

//    inf.close();

    M = Mtrans;

    M(4, 4) = Mlon(4, 4);
    M(4, 5) = Mlon(4, 5);
    M(5, 4) = Mlon(5, 4);
    M(5, 5) = Mlon(5, 5);
}


void ElementRFCavity::GetCavMat(const int cavi, const int cavilabel, const double Rm, Particle &real,
                                const double EfieldScl, const double IonFyi_s,
                                const double IonEk_s, const double fRF, state_t::matrix_t &M)
{
    int    n;
    double IonLambda, Ecen[2], T[2], Tp[2], S[2], Sp[2], V0[2];
    double dis, IonW_s[3], IonFy_s[3], gamma_s[3], beta_s[3], IonK_s[3];
    double IonK[2];

    IonLambda  = C0/fRF*MtoMM;

    IonW_s[0]  = IonEk_s + real.IonEs;
    IonFy_s[0] = IonFyi_s;
    gamma_s[0] = IonW_s[0]/real.IonEs;
    beta_s[0]  = sqrt(1e0-1e0/sqr(gamma_s[0]));
    IonK_s[0]  = 2e0*M_PI/(beta_s[0]*IonLambda);

    n   = CavData.s.size();
    dis = (CavData.s[n-1]-CavData.s[0])/2e0;

    TransFacts(cavilabel, beta_s[0], 1, EfieldScl, Ecen[0], T[0], Tp[0], S[0], Sp[0], V0[0]);
    EvalGapModel(dis, IonW_s[0], real, IonFy_s[0], IonK_s[0], IonLambda,
                Ecen[0], T[0], S[0], Tp[0], Sp[0], V0[0], IonW_s[1], IonFy_s[1]);
    gamma_s[1] = IonW_s[1]/real.IonEs;
    beta_s[1]  = sqrt(1e0-1e0/sqr(gamma_s[1]));
    IonK_s[1]  = 2e0*M_PI/(beta_s[1]*IonLambda);

    TransFacts(cavilabel, beta_s[1], 2, EfieldScl, Ecen[1], T[1], Tp[1], S[1], Sp[1], V0[1]);
    EvalGapModel(dis, IonW_s[1], real, IonFy_s[1], IonK_s[1], IonLambda,
                Ecen[1], T[1], S[1], Tp[1], Sp[1], V0[1], IonW_s[2], IonFy_s[2]);
    gamma_s[2] = IonW_s[2]/real.IonEs;
    beta_s[2]  = sqrt(1e0-1e0/sqr(gamma_s[2]));
    IonK_s[2]  = 2e0*M_PI/(beta_s[2]*IonLambda);

    Ecen[0] = Ecen[0] - dis;

    double TTF_tab[] = {Ecen[0], T[0], Tp[0], S[0], Sp[0], V0[0], Ecen[1], T[1], Tp[1], S[1], Sp[1], V0[1]};
    IonK[0] = (IonK_s[0]+IonK_s[1])/2e0;
    IonK[1] = (IonK_s[1]+IonK_s[2])/2e0;

    if (false) {
        printf("\n GetCavMat:\n");
        printf("IonK  : %15.10f %15.10f %15.10f\n", IonK_s[0], IonK_s[1], IonK_s[2]);
        printf("IonK  : %15.10f %15.10f\n", IonK[0], IonK[1]);
        printf("beta  : %15.10f %15.10f %15.10f\n", beta_s[0], beta_s[1], beta_s[2]);
        printf("gamma : %15.10f %15.10f %15.10f\n", gamma_s[0], gamma_s[1], gamma_s[2]);
        printf("Ecen  : %15.10f %15.10f\n", Ecen[0], Ecen[1]);
        printf("T     : %15.10f %15.10f\n", T[0], T[1]);
        printf("Tp    : %15.10f %15.10f\n", Tp[0], Tp[1]);
        printf("S     : %15.10f %15.10f\n", S[0], S[1]);
        printf("Sp    : %15.10f %15.10f\n", Sp[0], Sp[1]);
        printf("V0    : %15.10f %15.10f\n", V0[0], V0[1]);
    }

    this->ElementRFCavity::GetCavMatParams(cavi, beta_s, gamma_s, IonK, CavTLMLineTab);
    this->ElementRFCavity::GenCavMat(cavi, dis, EfieldScl, TTF_tab, beta_s, gamma_s, IonLambda, real, IonFy_s, Rm, M);
}


void ElementRFCavity::GetCavBoost(const CavDataType &CavData, Particle &state, const double IonFy0, const double fRF,
                                  const double EfieldScl, double &IonFy)
{
    int    n = CavData.s.size(),
           k;

    double dis = CavData.s[n-1] - CavData.s[0],
           dz  = dis/(n-1),
           IonLambda, IonK, IonFylast, IonGamma, IonBeta;

    IonLambda = C0/fRF*MtoMM;

//    std::cout<<__FUNCTION__
//             <<" IonFy0="<<IonFy0
//             <<" fRF="<<fRF
//             <<" EfieldScl="<<EfieldScl
//             <<" state="<<state
//             <<"\n";

    IonFy = IonFy0;
    IonK  = state.SampleIonK;
    for (k = 0; k < n-1; k++) {
        IonFylast = IonFy;
        IonFy += IonK*dz;
        state.IonW  += state.IonZ*EfieldScl*(CavData.Elong[k]+CavData.Elong[k+1])/2e0
                       *cos((IonFylast+IonFy)/2e0)*dz/MtoMM;
        IonGamma = state.IonW/state.IonEs;
        IonBeta  = sqrt(1e0-1e0/sqr(IonGamma));
        if ((state.IonW-state.IonEs) < 0e0) {
            state.IonW = state.IonEs;
            IonBeta = 0e0;
        }
        IonK = 2e0*M_PI/(IonBeta*IonLambda);
//        std::cout<<" "<<k<<" IonK="<<IonK<<" IonW="<<state.IonW<<"\n";
    }
}


void ElementRFCavity::PropagateLongRFCav(const Config &conf, Particle &ref)
{
    std::string CavType;
    int         cavi;
    double      fRF, multip, IonFys, EfieldScl, caviFy, IonFy_i, IonFy_o;

    CavType = conf.get<std::string>("cavtype");
    if (CavType == "0.041QWR") {
        cavi = 1;
    } else if (conf.get<std::string>("cavtype") == "0.085QWR") {
        cavi = 2;
    } else {
        std::ostringstream strm;
        strm << "*** PropagateLongRFCav: undef. cavity type: " << CavType << "\n";
        throw std::runtime_error(strm.str());
    }

    fRF       = conf.get<double>("f");
    multip    = fRF/SampleFreq;
    IonFys    = conf.get<double>("phi")*M_PI/180e0;  // Synchrotron phase [rad].
    EfieldScl = conf.get<double>("scl_fac");         // Electric field scale factor.

    caviFy = GetCavPhase(cavi, ref, IonFys, multip);

    IonFy_i = multip*ref.phis + caviFy;
    phi_ref = caviFy;
//    std::cout<<"RF long phase"
//               " caviFy="<<caviFy
//             <<" multip="<<multip
//             <<" phis="<<ref.phis
//             <<" IonFy_i="<<IonFy_i
//             <<" EfieldScl="<<EfieldScl
//             <<"\n";

    // For the reference particle, evaluate the change of:
    // kinetic energy, absolute phase, beta, and gamma.
    this->GetCavBoost(CavData, ref, IonFy_i, fRF, EfieldScl, IonFy_o);

    ref.IonEk       = ref.IonW - ref.IonEs;
    ref.recalc();
    ref.phis       += (IonFy_o-IonFy_i)/multip;
}


void ElementRFCavity::InitRFCav(const Config &conf, Particle &real, state_t::matrix_t &M)
{
    std::string CavType;
    int         cavi, cavilabel, multip;
    double      Rm, IonFy_i, Ek_i, fRF, EfieldScl, IonFy_o;
    double accIonW, avebeta, avegamma;

//    std::cout<<"RF recompute start "<<real<<"\n";

    CavType = conf.get<std::string>("cavtype");
    if (CavType == "0.041QWR") {
        cavi       = 1;
        cavilabel  = 41;
        multip     = 1;
        Rm         = 17e0;
    } else if (conf.get<std::string>("cavtype") == "0.085QWR") {
        cavi       = 2;
        cavilabel  = 85;
        multip     = 1;
        Rm         = 17e0;
    } else if (conf.get<std::string>("cavtype") == "0.29HWR") {
        cavi       = 3;
        cavilabel  = 29;
        multip     = 4;
        Rm         = 20e0;
    } else if (conf.get<std::string>("cavtype") == "0.53HWR") {
        cavi       = 4;
        cavilabel  = 53;
        multip     = 4;
        Rm         = 20e0;
    } else if (conf.get<std::string>("cavtype") == "??EL") {
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

    IonFy_i   = multip*real.phis + phi_ref;
    Ek_i      = real.IonEk;
    real.IonW = real.IonEk + real.IonEs;

    fRF       = conf.get<double>("f");
    EfieldScl = conf.get<double>("scl_fac");         // Electric field scale factor.

//    J.B.: Note, this was passed:
//    CaviIonK  = 2e0*M_PI*fRF/(real.beta*C0*MtoMM);
//    vs.:
//    double SampleIonK = 2e0*M_PI/(real.beta*C0/SampleFreq*MtoMM);

    ElementRFCavity::GetCavBoost(CavData, real, IonFy_i, fRF, EfieldScl, IonFy_o);

    real.IonEk       = real.IonW - real.IonEs;
    real.recalc();
    real.phis       += (IonFy_o-IonFy_i)/multip;

//    std::cout<<"RF recompute before "<<real
//             <<" cavi="<<cavi
//             <<" cavilabel="<<cavilabel
//             <<" Rm="<<Rm
//             <<" EfieldScl="<<EfieldScl
//             <<" IonFy_i="<<IonFy_i
//             <<" Ek_i="<<Ek_i
//             <<" fRF="<<fRF
//             <<"\n";

    this->GetCavMat(cavi, cavilabel, Rm, real, EfieldScl, IonFy_i, Ek_i, fRF, M);

//    std::cout<<"RF recompute after  "<<real<<"\n"
//             <<" YY "<<M<<"\n"
//             ;
}
