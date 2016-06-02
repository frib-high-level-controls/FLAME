#include "scsi/constants.h"
#include "scsi/moment2.h"
#include "scsi/chg_stripper.h"


void GetCenofChg(const Config &conf, std::vector<boost::shared_ptr<StateBase> > &ST,
                 Moment2State::vector_t &CenofChg, Moment2State::vector_t &BeamRMS)
{
    int                    i, j, n;
    double                 Ntot;
    std::vector<double>    NChg;
    Moment2State::vector_t BeamVar;

    CenofChg = boost::numeric::ublas::zero_vector<double>(PS_Dim);
    BeamRMS  = boost::numeric::ublas::zero_vector<double>(PS_Dim);
    BeamVar  = boost::numeric::ublas::zero_vector<double>(PS_Dim);

    n = conf.get<std::vector<double> >("IonChargeStates").size();
    NChg.resize(n);
    NChg = conf.get<std::vector<double> >("NCharge");

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


void StripperCharge(const double beta, double &Q_ave, double &d)
{
    // Baron's formula for carbon foil.
    double Q_ave1, Y;

    Q_ave1 = Stripper_IonProton*(1e0-exp(-83.275*(beta/pow(Stripper_IonProton, 0.447))));
    Q_ave  = Q_ave1*(1e0-exp(-12.905+0.2124*Stripper_IonProton-0.00122*sqr(Stripper_IonProton)));
    Y      = Q_ave1/Stripper_IonProton;
    d      = sqrt(Q_ave1*(0.07535+0.19*Y-0.2654*sqr(Y)));
}


void ChargeStripper(const double beta, const std::vector<double> ChgState, double chargeAmount_Baron[])
{
    int    k;
    double Q_ave, d;

    StripperCharge(beta, Q_ave, d);
    for (k = 0; k < ChgState.size(); k++)
        chargeAmount_Baron[k] = Gaussian(ChgState[k]*Stripper_IonMass, Q_ave, d);
}


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


std::vector<double> GetStrChgState(const Config &conf)
{
    return conf.get<std::vector<double> >("Stripper_IonChargeStates");
}


void Stripper_GetMat(const Config &conf, std::vector<boost::shared_ptr<Machine> > &sim,
                     std::vector<boost::shared_ptr<StateBase> > &ST, std::vector<double> ChgState)
{
    int                    k, n;
    double                 tmptotCharge, Fy_abs_recomb, Ek_recomb, stdEkFoilVariation, ZpAfStr, growthRate;
    double                 stdXYp, XpAfStr, growthRateXp, YpAfStr, growthRateYp, s;
    std::vector<double>    NChg;
    Particle               ref;
    state_t                *StatePtr;
    Moment2State::vector_t CenofChg, BeamRMS;
    Moment2State::matrix_t tmpmat;

    // Evaluate beam parameter recombination.

    GetCenofChg(conf, ST, CenofChg, BeamRMS);

    n = conf.get<std::vector<double> >("IonChargeStates").size();
    NChg.resize(n);
    NChg = conf.get<std::vector<double> >("NCharge");

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

    // Get new charge states.
    ChgState.clear();
    n = GetStrChgState(conf).size();
    ChgState.resize(n);
    ChgState = GetStrChgState(conf);

    // Propagate reference particle.
    ref = dynamic_cast<state_t*>(ST[0].get())->ref;
    Stripper_Propagate_ref(conf, ref, ChgState);

    StatePtr = dynamic_cast<state_t*>(ST[0].get());
    s = StatePtr->pos;

    sim.clear();
    ST.clear();
    for (k = 0; k < n; k++) {
        sim.push_back(boost::shared_ptr<Machine> (new Machine(conf)));
        ST.push_back(boost::shared_ptr<StateBase> (sim[k]->allocState()));

//        Move iterator to after charge splitter.

        state_t* StatePtr = dynamic_cast<state_t*>(ST[k].get());

        // Length is zero.
        StatePtr->pos        = s;

        StatePtr->ref        = ref;

        StatePtr->real.IonZ  = ChgState[k];
        StatePtr->real.IonEs = ref.IonEs;
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
