

#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>

#include <boost/numeric/ublas/io.hpp>

#include <scsi/constants.h>
#include <scsi/base.h>
#include <scsi/moment2.h>
#include <scsi/state/vector.h>
#include <scsi/state/matrix.h>


std::string HomeDir = "";

typedef Moment2State state_t;

typedef boost::numeric::ublas::vector<double> value_vec;
typedef boost::numeric::ublas::matrix<double> value_mat;

extern int glps_debug;

enum {maxsize = 7};

typedef boost::numeric::ublas::vector<double,
                boost::numeric::ublas::bounded_array<double, maxsize>
> vector_t;
typedef boost::numeric::ublas::matrix<double,
                boost::numeric::ublas::row_major,
                boost::numeric::ublas::bounded_array<double, maxsize*maxsize>
> matrix_t;

// Phase space dimension; including vector for orbit/1st moment.
# define PS_Dim       7


void PrtVec(const std::vector<double> &a)
{
    for (size_t k = 0; k < a.size(); k++)
        std::cout << std::scientific << std::setprecision(10)
                      << std::setw(18) << a[k];
    std::cout << "\n";
}


void PrtVec(const vector_t &a)
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


std::vector<double> GetChgState(Config &conf, const std::string &CSstr)
{
    return conf.get<std::vector<double> >(CSstr);
}


std::vector<double> GetNChg(Config &conf, const std::string &CAstr)
{
    return conf.get<std::vector<double> >(CAstr);
}


value_vec GetBaryCenter(Config &conf, const std::string &BCstr)
{

    const std::vector<double>& BCvec = conf.get<std::vector<double> >(BCstr);
    value_vec BC(PS_Dim);
    if (BCvec.size() > BC.data().size())
        throw std::invalid_argument("Initial state size too big");
    std::copy(BCvec.begin(), BCvec.end(), BC.data().begin());

    return BC;
}


value_mat GetBeamEnvelope(Config &conf, const std::string &BEstr)
{

    const std::vector<double>& BEvec = conf.get<std::vector<double> >(BEstr);
    value_mat BE(PS_Dim, PS_Dim);
    if (BEvec.size() > BE.data().size())
        throw std::invalid_argument("Initial state size too big");
    std::copy(BEvec.begin(), BEvec.end(), BE.data().begin());

    return BE;
}


void propagate(std::auto_ptr<Config> conf)
{
    const int nChgStates = 2;

    int                                      k;
    double                                   IonEk, IonEs, IonW, Fy_absState, EkState;
    std::vector<double>                      ChgState, NChg;
    value_vec                                BC[nChgStates];
    value_mat                                BE[nChgStates];
    std::vector<boost::shared_ptr<Machine> > sims;
    Machine::p_elements_t::iterator          it;


    for (k = 0; k < nChgStates; k++) {
        BC[k].resize(PS_Dim);
        BE[k].resize(PS_Dim, PS_Dim);
    }

    ChgState.resize(nChgStates);
    ChgState = GetChgState(*conf, "IonChargeStates");

    NChg.resize(nChgStates);
    NChg = GetNChg(*conf, "NCharge");

    BC[0] = GetBaryCenter(*conf, "BaryCenter1");
    BC[1] = GetBaryCenter(*conf, "BaryCenter2");

    BE[0] = GetBeamEnvelope(*conf, "S1");
    BE[1] = GetBeamEnvelope(*conf, "S2");

    std::cout << "\nIon charge states:\n";
    for (k = 0; k < nChgStates; k++)
        std::cout << std::fixed << std::setprecision(5) << std::setw(9) << ChgState[k];
    std::cout << "\n";
    std::cout << "\nIon charge amount:\n";
    for (k = 0; k < nChgStates; k++)
        std::cout << std::fixed << std::setprecision(1) << std::setw(8) << NChg[k];
    std::cout << "\n";
    std::cout << "\nBarycenter:\n";
    PrtVec(BC[0]);
    PrtVec(BC[1]);
    std::cout << "\nBeam envelope:\n";
    PrtMat(BE[0]);
    std::cout << "\n";
    PrtMat(BE[1]);

    for (k = 0; k < nChgStates; k++) {
        sims.push_back(boost::shared_ptr<Machine> (new Machine(*conf)));
        sims[k]->set_trace(NULL);
    }

//        std::cout << "# Machine configuration\n" << *sims[0] << "\n\n";

    Config  D;
    state_t *StatePtr[nChgStates];

    std::cout << "\n";
    for (k = 0; k < nChgStates; k++) {
        std::auto_ptr<StateBase> state(sims[k]->allocState(D));

        it = sims[k]->p_elements.begin();
        // Skip over state.
        it++;

        // Propagate through first element (beam initial conditions).
        sims[k]->propagate(state.get(), 0, 1);

        IonEk = state->IonEk/MeVtoeV;
        IonEs = state->IonEs/MeVtoeV;
        IonW  = state->IonW/MeVtoeV;

        // Define initial conditions.
        Fy_absState = BC[k][state_t::PS_S];
        EkState     = IonEk + BC[k][state_t::PS_PS];

        // Initialize state.
        StatePtr[k] = dynamic_cast<state_t*>(state.get());
        StatePtr[k]->moment0 = BC[k];
        StatePtr[k]->state   = BE[k];

        StatePtr[k]->Fy_absState = StatePtr[k]->moment0[state_t::PS_S];
        StatePtr[k]->EkState     = StatePtr[k]->IonEk + StatePtr[k]->moment0[state_t::PS_PS]*MeVtoeV;

        std::cout << std::fixed << std::setprecision(5)
                  << "  IonZ = " << ChgState[k]
                  << ",  IonEs [Mev/u] = " << IonEs << ", IonEk [Mev/u] = " << IonEk
                  << ", IonW [Mev/u] = " << IonW << "\n";
    }

    //    for (k = 0; k < nChgStates; k++) {
    for (k = 0; k < 1; k++) {
        std::auto_ptr<StateBase> state(sims[k]->allocState(D));

        StatePtr[k] = dynamic_cast<state_t*>(state.get());

        std::cout << std::fixed << std::setprecision(3) << "\ns [m] = " << StatePtr[k]->pos << "\n";

        // Propagate through first element (beam initial conditions).
        sims[k]->propagate(state.get(), 0, 1);

        IonEk = state->IonEk/MeVtoeV;
        IonEs = state->IonEs/MeVtoeV;
        IonW  = state->IonW/MeVtoeV;

        // Initialize state.
        StatePtr[k]->moment0 = BC[k];
        StatePtr[k]->state   = BE[k];

        // Approximate (E_k = m0*v^2/2 vs. p*c0).
        StatePtr[k]->gamma  = (IonEs+IonEk)/IonEs;
        StatePtr[k]->beta   = sqrt(1e0-1e0/sqr(StatePtr[k]->gamma));
        StatePtr[k]->bg1    = StatePtr[k]->beta*StatePtr[k]->gamma;

        StatePtr[k]->Fy_absState = StatePtr[k]->moment0[state_t::PS_S];
        StatePtr[k]->EkState     = StatePtr[k]->IonEk + StatePtr[k]->moment0[state_t::PS_PS]*MeVtoeV;

        it = sims[k]->p_elements.begin();
        // Skip over state.
        it++;

//        ElementVoid* elem = *it;
        for (; it != sims[k]->p_elements.end(); ++it) {
//            elem->advance(*state);
//            sims[k]->propagate(state.get(), elem->index+n-1, 1);
            (*it)->advance(*state);
        }

        std::cout << "\n";
        PrtVec(StatePtr[k]->moment0);
        std::cout << "\n";
        PrtMat(StatePtr[k]->state);
        std::cout << std::fixed << std::setprecision(3) << "\ns [m] = " << StatePtr[k]->pos << "\n";
    }
}


void GetCavTLMstr(void)
{
    std::fstream inf1, inf2;

    inf1.open((HomeDir+"/data/Multipole41/thinlenlon_41.txt").c_str(), std::ifstream::in);
    CavTLMstream2[0] << inf1.rdbuf();

    inf2.open((HomeDir+"/data/Multipole85/thinlenlon_85.txt").c_str(), std::ifstream::in);
    CavTLMstream2[1] << inf2.rdbuf();
}


int main(int argc, char *argv[])
{
    try {
        std::auto_ptr<Config> conf;
        FILE                  *inf = NULL;

        if(argc > 1) {
            HomeDir = argv[2];
            inf = fopen(argv[1], "r");
            if (!inf) {
                fprintf(stderr, "Failed to open %s\n", argv[1]);
                return 2;
            }
        }

        glps_debug = 0; // 0 or 1.

        try {
            GLPSParser P;
            conf.reset(P.parse(inf));
            fprintf(stderr, "Parsing succeeds\n");
        } catch(std::exception& e) {
            fprintf(stderr, "Parse error: %s\n", e.what());
            fclose(inf);
            return 1;
        }

//        std::cout<<"# Reduced lattice\n";
//        GLPSPrint(std::cout, *conf);
//        std::cout<<"\n";

        CavData2[0].RdData(HomeDir+"/data/axisData_41.txt");
        CavData2[1].RdData(HomeDir+"/data/axisData_85.txt");

        GetCavTLMstr();

        // register state and element types
        registerLinear();
//        registerMoment();
        registerMoment2();

        propagate(conf);

        return 0;
    } catch(std::exception& e) {
        std::cerr << "Main exception: " << e.what() << "\n";
        Machine::registeryCleanup();
        return 1;
    }
}
