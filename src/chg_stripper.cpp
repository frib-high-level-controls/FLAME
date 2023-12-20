#include "flame/constants.h"
#include "flame/moment.h"
#include "flame/moment_sup.h"
#include "flame/chg_stripper.h"

#define sqr(x)  ((x)*(x))
#define cube(x) ((x)*(x)*(x))

static
double Gaussian(double in, const double Q_ave, const double d)
{
    // Gaussian distribution.
    return 1e0/sqrt(2e0*M_PI)/d*exp(-0.5e0*sqr(in-Q_ave)/sqr(d));
}


void ElementStripper::StripperCharge(const double beta, double &Q_ave, double &d)
{
    // Baron's formula for carbon foil.
    double Q_ave1, Y;

    Q_ave1 = Stripper_IonProton*(1e0-exp(-83.275*(beta/pow(Stripper_IonProton, 0.447))));
    Q_ave  = Q_ave1*(1e0-exp(-12.905+0.2124*Stripper_IonProton-0.00122*sqr(Stripper_IonProton)));
    Y      = Q_ave1/Stripper_IonProton;
    d      = sqrt(Q_ave1*(0.07535+0.19*Y-0.2654*sqr(Y)));
}


void ElementStripper::ChargeStripper(const double beta, const std::vector<double>& ChgState, std::vector<double>& chargeAmount_Baron)
{
    unsigned    k;
    double Q_ave, d;

    StripperCharge(beta, Q_ave, d);
    for (k = 0; k < ChgState.size(); k++)
        chargeAmount_Baron.push_back(Gaussian(ChgState[k]*Stripper_IonMass, Q_ave, d));
}


void ElementStripper::Stripper_Propagate_ref(Particle &ref)
{

    // Change reference particle charge state.
    ref.IonZ = Stripper_IonZ;

//    ChargeStripper(ref.beta, ChgState, chargeAmount_Baron);

    // Evaluate change in reference particle energy due to stripper model energy straggling.
    ref.IonEk      = (ref.IonEk-Stripper_Para[2])*Stripper_E0Para[1] + Stripper_E0Para[0];

    ref.recalc();
}

void ElementStripper::Stripper_GetMat(const Config &conf,
                     MomentState &ST)
{
    unsigned               k, n;
    double                 tmptotCharge, Fy_abs_recomb, Ek_recomb, stdEkFoilVariation, ZpAfStr, growthRate;
    double                 stdXYp, XpAfStr, growthRateXp, YpAfStr, growthRateYp, s;
    Particle               ref;
    state_t                *StatePtr = &ST;
    MomentState::matrix_t  tmpmat;
    std::vector<double>    chargeAmount_Set;

    //std::cout<<"In "<<__FUNCTION__<<"\n";

    // Get new charge states and stripper parameters.
    const std::vector<double>& ChgState = conf.get<std::vector<double> >("IonChargeStates");
    const std::vector<double>& NChg = conf.get<std::vector<double> >("NCharge");
    const std::string chrgmdl = conf.get<std::string>("charge_model", "baron");

    if(chrgmdl=="off" && ChgState.size()!=NChg.size())
        throw std::runtime_error("charge stripper requires that IonChargeStates[] and NCharge[] have the same length");
    if(chrgmdl!="off" && chrgmdl!="baron")
        throw std::runtime_error("charge_model key word unknown, only \"baron\" and \"off\" supported by now");

    Stripper_IonZ = conf.get<double>("Stripper_IonZ", Stripper_IonZ_default);
    Stripper_IonMass = conf.get<double>("Stripper_IonMass", Stripper_IonMass_default);
    Stripper_IonProton = conf.get<double>("Stripper_IonProton", Stripper_IonProton_default);
    Stripper_E1Para = conf.get<double>("Stripper_E1Para", Stripper_E1Para_default);
    Stripper_lambda = conf.get<double>("Stripper_lambda", Stripper_lambda_default);
    Stripper_upara = conf.get<double>("Stripper_upara", Stripper_upara_default);

    const std::vector<double> p1_default(Stripper_Para_default, Stripper_Para_default+3),
                              p2_default(Stripper_E0Para_default, Stripper_E0Para_default+3);

    Stripper_Para = conf.get<std::vector<double> >("Stripper_Para", p1_default);
    Stripper_E0Para = conf.get<std::vector<double> >("Stripper_E0Para", p2_default);

    n = ChgState.size();

    if(chrgmdl=="off" )
        if (ChgState.size()!=NChg.size())
            throw std::runtime_error("charge stripper requires that IonChargeStates[] and NCharge[] have the same length");
        else{
            for (k = 0; k < n; k++) chargeAmount_Set.push_back(NChg[k]);
        }
    else
        ChargeStripper(StatePtr->real[0].beta, ChgState, chargeAmount_Set);


    // Evaluate beam parameter recombination.

    tmptotCharge  = 0e0;
    Fy_abs_recomb = 0e0;
    Ek_recomb     = 0e0;
    tmpmat        = boost::numeric::ublas::zero_matrix<double>(PS_Dim);
    for (k = 0; k < ST.size(); k++) {
        const double Q = ST.real[k].IonQ;
        tmptotCharge  += Q;
        Fy_abs_recomb += Q*StatePtr->real[k].phis;
        Ek_recomb     += Q*StatePtr->real[k].IonEk;
        tmpmat        +=
                Q*(StatePtr->moment1[k]+outer_prod(StatePtr->moment0[k]-ST.moment0_env, StatePtr->moment0[k]-ST.moment0_env));
    }

    Fy_abs_recomb /= tmptotCharge;
    Ek_recomb     /= tmptotCharge;

    // Stripper model energy straggling.
    Ek_recomb      = (Ek_recomb-Stripper_Para[2])*Stripper_E0Para[1] + Stripper_E0Para[0];
    tmpmat        /= tmptotCharge;

    // Estimate Zp stripper caused envelope increase.
    stdEkFoilVariation = sqrt(1e0/3e0)*Stripper_Para[1]/100e0*Stripper_Para[0]*Stripper_E0Para[2];
    ZpAfStr            = tmpmat(5, 5) + sqr(Stripper_E1Para) + sqr(stdEkFoilVariation);
    if (tmpmat(5, 5) != 0e0) {
        growthRate         = sqrt(ZpAfStr/tmpmat(5, 5));
    } else {
        growthRate         = 1e0;
    }

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


    // Propagate reference particle.
    ref = ST.ref;
    Stripper_Propagate_ref(ref);

    s = StatePtr->pos;

    ST.real.resize(n);
    ST.moment0.resize(n);
    ST.moment1.resize(n);
    ST.transmat.resize(n);

    // Length is zero.
    StatePtr->pos = s;

    StatePtr->ref = ref;


    for (k = 0; k < n; k++) {
        StatePtr->real[k].IonZ  = ChgState[k];
        StatePtr->real[k].IonQ  = chargeAmount_Set[k];
        StatePtr->real[k].SampleFreq = ref.SampleFreq;
        StatePtr->real[k].IonEs = ref.IonEs;
        StatePtr->real[k].IonEk = Ek_recomb;
        StatePtr->real[k].recalc();
        StatePtr->real[k].phis  = Fy_abs_recomb;
        StatePtr->moment0[k]    = ST.moment0_env;
        StatePtr->moment1[k]    = tmpmat;
        StatePtr->transmat[k]   = boost::numeric::ublas::identity_matrix<double>(PS_Dim);
    }

    ST.calc_rms();
}

void ElementStripper::advance(StateBase &s)
{
    state_t& ST = static_cast<state_t&>(s);

    ST.recalc();
    ST.calc_rms(); // paranoia in case someone (python) has made moment0_env inconsistent

    if(ST.retreat) throw std::runtime_error(SB()<<
        "Backward propagation error: Backward propagation does not support charge stripper.");


    Stripper_GetMat(conf(), ST);
}
