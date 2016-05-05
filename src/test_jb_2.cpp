

#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <ctime>

#include <boost/numeric/ublas/io.hpp>

#include <scsi/constants.h>
#include <scsi/base.h>
#include <scsi/moment2.h>
#include <scsi/state/vector.h>
#include <scsi/state/matrix.h>


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

    int                                      k;
    std::vector<value_vec>                   BC;
    std::vector<value_mat>                   BE;
    std::vector<boost::shared_ptr<Machine> > sims;
    Machine::iterator          it;

    std::vector<double>                      ChgState(GetChgState(*conf, "IonChargeStates")),
                                             NChg(GetNChg(*conf, "NCharge"));

    const int nChgStates = ChgState.size();

    BC.resize(nChgStates);
    BE.resize(nChgStates);

    for (k = 0; k < nChgStates; k++) {
        {
            std::ostringstream strm;
            strm<<"BaryCenter"<<k;
            BC[k] = GetBaryCenter(*conf, strm.str());
        }

        {
            std::ostringstream strm;
            strm<<"S"<<k;
            BE[k] = GetBeamEnvelope(*conf, strm.str());
        }
    }

    std::cout << "\nIon charge states:\n";
    for (k = 0; k < nChgStates; k++)
        std::cout << std::fixed << std::setprecision(5) << std::setw(9) << ChgState[k]
                  << std::setprecision(1) << std::setw(8) << NChg[k] << "\n";

    std::cout << "\nBarycenter:\n";
    for (k = 0; k < nChgStates; k++) PrtVec(BC[k]);
    std::cout << "\nBeam envelope:\n";
    for (k = 0; k < nChgStates; k++) {
        if(k!=0) std::cout<<"\n";
        PrtMat(BE[k]);
    }

    for (k = 0; k < nChgStates; k++) {
        conf->set<double>("cstate", k);
        sims.push_back(boost::shared_ptr<Machine> (new Machine(*conf)));
        sims[k]->set_trace(NULL);
    }

//        std::cout << "# Machine configuration\n" << *sims[0] << "\n\n";

    std::vector<state_t *> StatePtr(nChgStates);

    for (k = 0; k < nChgStates; k++) {
        std::auto_ptr<StateBase> state(sims[k]->allocState());

        StatePtr[k] = dynamic_cast<state_t*>(state.get());
        if(!StatePtr[k]) throw std::runtime_error("Only sim_type MomentMatrix2 is supported");

        std::cout << std::fixed << std::setprecision(3) << "\ns [m] = " << StatePtr[k]->pos << "\n";

        // Propagate through first element (beam initial conditions).
        sims[k]->propagate(state.get(), 0, 1);

        clock_t tStamp[2];
        tStamp[0] = clock();
        sims[k]->propagate(state.get(), 1);
        tStamp[1] = clock();

        std::cout << std::fixed << std::setprecision(5)
                  << "\npropagate: " << double(tStamp[1]-tStamp[0])/CLOCKS_PER_SEC << " sec" << "\n";

        std::cout << "\n";
        PrtVec(StatePtr[k]->moment0);
        std::cout << "\n";
        PrtMat(StatePtr[k]->state);
        std::cout << std::fixed << std::setprecision(3) << "\ns [m] = " << StatePtr[k]->pos << "\n";
    }
}


int main(int argc, char *argv[])
{
    try {
        std::auto_ptr<Config> conf;

        if(argc>2)

        glps_debug = 0; // 0 or 1.

        try {
            GLPSParser P;
            conf.reset(P.parse_file(argc>1 ? argv[1] : NULL));
            fprintf(stderr, "Parsing succeeds\n");
        } catch(std::exception& e) {
            fprintf(stderr, "Parse error: %s\n", e.what());
            return 1;
        }

//        std::cout<<"# Reduced lattice\n";
//        GLPSPrint(std::cout, *conf);
//        std::cout<<"\n";

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
