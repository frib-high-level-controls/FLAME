#include "scsi/constants.h"
#include "scsi/moment2.h"
#include "scsi/chg_stripper.h"


void GetCenofChg(std::auto_ptr<Config> &conf, std::vector<boost::shared_ptr<StateBase> > &ST,
                 Moment2State::vector_t &CenofChg, Moment2State::vector_t &BeamRMS)
{
    int                    i, j, n;
    double                 Ntot;
    std::vector<double>    NChg;
    Moment2State::vector_t BeamVar;

    CenofChg = boost::numeric::ublas::zero_vector<double>(PS_Dim);
    BeamRMS  = boost::numeric::ublas::zero_vector<double>(PS_Dim);
    BeamVar  = boost::numeric::ublas::zero_vector<double>(PS_Dim);

    n = conf->get<std::vector<double> >("IonChargeStates").size();
    NChg.resize(n);
    NChg = conf->get<std::vector<double> >("NCharge");

    Ntot = 0e0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < PS_Dim; j++) {
            state_t* StatePtr = dynamic_cast<state_t*>(ST[i].get());
            CenofChg[j] += NChg[i]*StatePtr->moment0[j];
        }
        Ntot += NChg[i];
    }

    for (j = 0; j < PS_Dim; j++)
        CenofChg[j] /= Ntot;

    for (i = 0; i < n; i++) {
        for (j = 0; j < PS_Dim; j++) {
            state_t* StatePtr = dynamic_cast<state_t*>(ST[i].get());
            BeamVar[j]  +=
                    NChg[i]*(StatePtr->state(j, j)
                    +(StatePtr->moment0[j]-CenofChg[j])*(StatePtr->moment0[j]-CenofChg[j]));
        }
    }

    for (j = 0; j < PS_Dim; j++)
        BeamRMS[j] = sqrt(BeamVar[j]/Ntot);
}


double Gaussian(double in, const double Q_ave, const double d)
{
    // Gaussian distribution.
    return 1e0/sqrt(2e0*M_PI)/d*exp(-0.5e0*sqr(in-Q_ave)/sqr(d));
}


void StripperCharge(const double IonProton, const double beta, double &Q_ave, double &d)
{
    // Use Baron's formula for carbon foil.
    double Q_ave1, Y;

    Q_ave1 = IonProton*(1e0-exp(-83.275*(beta/pow(IonProton, 0.447))));
    Q_ave  = Q_ave1*(1e0-exp(-12.905+0.2124*IonProton-0.00122*sqr(IonProton)));
    Y      = Q_ave1/IonProton;
    d      = sqrt(Q_ave1*(0.07535+0.19*Y-0.2654*sqr(Y)));
}


void ChargeStripper(const double IonMass, const double IonProton, const double beta,
                    const int nChargeStates, const double IonChargeStates[],
                    double chargeAmount_Baron[])
{
    int    k;
    double Q_ave, d;

    StripperCharge(IonProton, beta, Q_ave, d);
    for (k = 0; k < nChargeStates; k++)
        chargeAmount_Baron[k] = Gaussian(IonChargeStates[k]*IonMass, Q_ave, d);
}


void PropagateLongStripper(const Config &conf, const int n, double &IonZ, const double IonEs,
                           double &IonW, double &SampleIonK, double &IonBeta)
{
    double IonEk, IonGamma, chargeAmount_Baron[Stripper_n_chg_states];

    IonZ = Stripper_IonZ;
    ChargeStripper(Stripper_IonMass, Stripper_IonProton, IonBeta,
                   Stripper_n_chg_states, Stripper_IonChargeStates,
                   chargeAmount_Baron);
    // Evaluate change in reference particle energy due to stripper model energy straggling.
//    IonEk      = (LongTab.Ek[n-2]-StripperPara[2])*Stripper_E0Para[1] + Stripper_E0Para[0];
    IonW       = IonEk + IonEs;
    IonGamma   = IonW/IonEs;
    IonBeta    = sqrt(1e0-1e0/sqr(IonGamma));
    SampleIonK = 2e0*M_PI/(IonBeta*SampleLambda);

//    LongTab.set(LongTab.s[n-2], IonEk, LongTab.FyAbs[n-2], IonBeta, IonGamma);

    // chargeAmount = fribstripper.chargeAmount_Baron;
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


void Stripper_GetMat(std::auto_ptr<Config> &conf, std::vector<boost::shared_ptr<Machine> > &sim,
                     std::vector<boost::shared_ptr<StateBase> > &ST)
{
    int                    k, n;
    double                 tmptotCharge, Fy_abs_recomb, Ek_recomb, stdEkFoilVariation, ZpAfStr, growthRate;
    double                 stdXYp, XpAfStr, growthRateXp, YpAfStr, growthRateYp;
    std::vector<double>    NChg;
    Particle               ref;
    state_t                *StatePtr;
    Moment2State::vector_t CenofChg, BeamRMS;
    Moment2State::matrix_t tmpmat;

    // Evaluate beam parameter recombination.

    GetCenofChg(conf, ST, CenofChg, BeamRMS);

    n = conf->get<std::vector<double> >("IonChargeStates").size();
    NChg.resize(n);
    NChg = conf->get<std::vector<double> >("NCharge");

    tmptotCharge  = 0e0;
    Fy_abs_recomb = 0e0;
    Ek_recomb     = 0e0;
    tmpmat        = boost::numeric::ublas::zero_matrix<double>(PS_Dim);
    for (k = 0; k < ST.size(); k++) {
        StatePtr = dynamic_cast<state_t*>(ST[k].get());

        tmptotCharge  += NChg[k];
        Fy_abs_recomb += NChg[k]*StatePtr->real.phis;
        Ek_recomb     += NChg[k]*StatePtr->real.IonEk;
        tmpmat        +=
                NChg[k]*(StatePtr->state+outer_prod(StatePtr->moment0-CenofChg, StatePtr->moment0-CenofChg));
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

    ref = dynamic_cast<state_t*>(ST[0].get())->ref;

    sim.clear();
    ST.clear();
    for (k = 0; k < Stripper_n_chg_states; k++) {
        sim.push_back(boost::shared_ptr<Machine> (new Machine(*conf)));
        ST.push_back(boost::shared_ptr<StateBase> (sim[k]->allocState()));

//        Move iterator to after charge splitter.
//        ST.pos += length;

        state_t* StatePtr = dynamic_cast<state_t*>(ST[k].get());

        StatePtr->ref        = ref;
        StatePtr->ref.IonZ   = Stripper_IonZ;

        StatePtr->real.IonZ  = Stripper_IonZ;
        StatePtr->real.IonEk = Ek_recomb;
        StatePtr->real.IonW  = StatePtr->real.IonEk + StatePtr->ref.IonEs;
        StatePtr->real.gamma = StatePtr->real.IonW/StatePtr->ref.IonEs;
        StatePtr->real.beta  = sqrt(1e0-1e0/sqr(StatePtr->real.gamma));
        StatePtr->real.phis  = Fy_abs_recomb;
        StatePtr->real.IonEk = Ek_recomb;
        StatePtr->moment0    = CenofChg;
        StatePtr->state      = tmpmat;
    }
}
