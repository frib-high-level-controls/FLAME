#include "scsi/constants.h"
#include "scsi/moment2.h"
#include "scsi/moment2_sup.h"
#include "scsi/chg_stripper.h"

void GetCenofChg(const Config &conf, const Moment2State &ST,
                 Moment2State::vector_t &CenofChg, Moment2State::vector_t &BeamRMS)
{
    size_t i, j;
    const size_t n = ST.real.size(); // # of states

    CenofChg = boost::numeric::ublas::zero_vector<double>(PS_Dim);
    BeamRMS  = boost::numeric::ublas::zero_vector<double>(PS_Dim);
    Moment2State::vector_t BeamVar(boost::numeric::ublas::zero_vector<double>(PS_Dim));

    double Ntot = 0e0;
    for (i = 0; i < n; i++) {
        CenofChg += ST.moment0[i]*ST.real[i].IonQ;
        Ntot += ST.real[i].IonQ;
    }

    CenofChg /= Ntot;

    for (i = 0; i < n; i++) {
        for (j = 0; j < Moment2State::maxsize; j++) {
            BeamVar[j]  +=
                    ST.real[i].IonQ*(ST.moment1[i](j, j)
                    +(ST.moment0[i][j]-CenofChg[j])*(ST.moment0[i][j]-CenofChg[j]));
        }
    }

    for (j = 0; j < PS_Dim; j++)
        BeamRMS[j] = sqrt(BeamVar[j]/Ntot);
}


static
double Gaussian(double in, const double Q_ave, const double d)
{
    // Gaussian distribution.
    return 1e0/sqrt(2e0*M_PI)/d*exp(-0.5e0*sqr(in-Q_ave)/sqr(d));
}


static
void StripperCharge(const double beta, double &Q_ave, double &d)
{
    // Baron's formula for carbon foil.
    double Q_ave1, Y;

    Q_ave1 = Stripper_IonProton*(1e0-exp(-83.275*(beta/pow(Stripper_IonProton, 0.447))));
    Q_ave  = Q_ave1*(1e0-exp(-12.905+0.2124*Stripper_IonProton-0.00122*sqr(Stripper_IonProton)));
    Y      = Q_ave1/Stripper_IonProton;
    d      = sqrt(Q_ave1*(0.07535+0.19*Y-0.2654*sqr(Y)));
}


static
void ChargeStripper(const double beta, const std::vector<double> ChgState, double chargeAmount_Baron[])
{
    unsigned    k;
    double Q_ave, d;

    StripperCharge(beta, Q_ave, d);
    for (k = 0; k < ChgState.size(); k++)
        chargeAmount_Baron[k] = Gaussian(ChgState[k]*Stripper_IonMass, Q_ave, d);
}


static
void Stripper_Propagate_ref(const Config &conf, Particle &ref, const std::vector<double> ChgState)
{
    const int n = ChgState.size();

    double chargeAmount_Baron[n];

    // Change reference particle charge state.
    ref.IonZ = Stripper_IonZ;

    ChargeStripper(ref.beta, ChgState, chargeAmount_Baron);

    // Evaluate change in reference particle energy due to stripper model energy straggling.
    ref.IonEk      = (ref.IonEk-Stripper_Para[2])*Stripper_E0Para[1] + Stripper_E0Para[0];

    ref.IonW       = ref.IonEk + ref.IonEs;
    ref.gamma      = ref.IonW/ref.IonEs;
    ref.beta       = sqrt(1e0-1e0/sqr(ref.gamma));
    ref.SampleIonK = 2e0*M_PI/(ref.beta*SampleLambda);

    // chargeAmount = fribstripper.chargeAmount_Baron;
}

static
void PrtMat1(const state_t::matrix_t &M)
{
    for (size_t j = 0; j < M.size1(); j++) {
        for (size_t k = 0; k < M.size2(); k++)
            std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << M(j, k);
        std::cout << "\n";
    }
}


static
std::vector<double> GetStrChgState(const Config &conf)
{
    return conf.get<std::vector<double> >("Stripper_IonChargeStates");
}


void Stripper_GetMat(const Config &conf,
                     Moment2State &ST, std::vector<double> ChgState)
{
    unsigned               k, n;
    double                 tmptotCharge, Fy_abs_recomb, Ek_recomb, stdEkFoilVariation, ZpAfStr, growthRate;
    double                 stdXYp, XpAfStr, growthRateXp, YpAfStr, growthRateYp, s;
    std::vector<double>    NChg;
    Particle               ref;
    state_t                *StatePtr = &ST;
    Moment2State::vector_t CenofChg, BeamRMS;
    Moment2State::matrix_t tmpmat;

    // Evaluate beam parameter recombination.

    GetCenofChg(conf, ST, CenofChg, BeamRMS);
    assert(CenofChg==ST.moment0_env);

    std::cout<<"In "<<__FUNCTION__<<"\n";

    n = conf.get<std::vector<double> >("IonChargeStates").size();
    NChg.resize(n);
    NChg = conf.get<std::vector<double> >("NCharge");

    tmptotCharge  = 0e0;
    Fy_abs_recomb = 0e0;
    Ek_recomb     = 0e0;
    tmpmat        = boost::numeric::ublas::zero_matrix<double>(PS_Dim);
    for (k = 0; k < ST.size(); k++) {

        tmptotCharge  += NChg[k];
        Fy_abs_recomb += NChg[k]*StatePtr->real[k].phis;
        Ek_recomb     += NChg[k]*StatePtr->real[k].IonEk;
        tmpmat        +=
                NChg[k]*(StatePtr->moment1[k]+outer_prod(StatePtr->moment0[k]-CenofChg, StatePtr->moment0[k]-CenofChg));
    }

    Fy_abs_recomb /= tmptotCharge;
    Ek_recomb     /= tmptotCharge;

    // Stripper model energy straggling.
    Ek_recomb      = (Ek_recomb-Stripper_Para[2])*Stripper_E0Para[1] + Stripper_E0Para[0];
    tmpmat        /= tmptotCharge;

    // Estimate Zp stripper caused envelope increase.
    stdEkFoilVariation = sqrt(1e0/3e0)*Stripper_Para[1]/100e0*Stripper_Para[0]*Stripper_E0Para[2];
    ZpAfStr            = tmpmat(5, 5) + sqr(Stripper_E1Para) + sqr(stdEkFoilVariation);
    growthRate         = sqrt(ZpAfStr/tmpmat(5, 5));

    for (k = 0; k < 6; k++) {
        tmpmat(k, 5) = tmpmat(k, 5)*growthRate;
        // This trick allows two times growthRate for <Zp, Zp>.
        tmpmat(5, k) = tmpmat(5, k)*growthRate;
    }

    // Estimate Xp, Yp stripper caused envelope increase.
    stdXYp       = sqrt(Stripper_upara/sqr(Stripper_lambda))*1e-3;// mrad to rad
    XpAfStr      = tmpmat(1, 1) + sqr(stdXYp);
    growthRateXp = sqrt(XpAfStr/tmpmat(1, 1));
    YpAfStr      = tmpmat(3, 3) + sqr(stdXYp);
    growthRateYp = sqrt(YpAfStr/tmpmat(3, 3));

    for (k = 0; k < 6; k++) {
        tmpmat(k, 1) = tmpmat(k, 1)*growthRateXp;
        // This trick allows two times growthRate for <Zp, Zp>.
        tmpmat(1, k) = tmpmat(1, k)*growthRateXp;
        tmpmat(k, 3) = tmpmat(k, 3)*growthRateYp;
        tmpmat(3, k) = tmpmat(3, k)*growthRateYp;
    }

    // Get new charge states.
    ChgState = GetStrChgState(conf);
    NChg = conf.get<std::vector<double> >("Stripper_NCharge");
    n = ChgState.size();
    assert(NChg.size()==n);

    // Propagate reference particle.
    ref = ST.ref;
    Stripper_Propagate_ref(conf, ref, ChgState);

    s = StatePtr->pos;

    ST.real.resize(n);
    ST.moment0.resize(n);
    ST.moment1.resize(n);

    for (k = 0; k < n; k++) {

        // Length is zero.
        StatePtr->pos        = s;

        StatePtr->ref        = ref;

        StatePtr->real[k].IonZ  = ChgState[k];
        StatePtr->real[k].IonQ  = NChg[k];
        StatePtr->real[k].IonEs = ref.IonEs;
        StatePtr->real[k].IonEk = Ek_recomb;
        StatePtr->real[k].IonW  = StatePtr->real[k].IonEk + StatePtr->ref.IonEs;
        StatePtr->real[k].gamma = StatePtr->real[k].IonW/StatePtr->ref.IonEs;
        StatePtr->real[k].beta  = sqrt(1e0-1e0/sqr(StatePtr->real[k].gamma));
        StatePtr->real[k].phis  = Fy_abs_recomb;
        StatePtr->real[k].IonEk = Ek_recomb;
        StatePtr->moment0[k]    = CenofChg;
        StatePtr->moment1[k]    = tmpmat;
    }

    assert(CenofChg==StatePtr->moment0_env);
}
