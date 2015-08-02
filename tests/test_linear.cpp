
#include <iostream>
#include <list>

#include <boost/numeric/ublas/io.hpp>

#include <scsi/linear.h>
#include <scsi/state/vector.h>
#include <scsi/state/matrix.h>

int main(int argc, char *argv[]) {

    registerLinear();
    Config conf;

    {
        std::cout << "Building description\n";
        std::list<boost::any> elements;

        conf.set<std::string>("sim-type", "Linear1D");

        {
            Config E;
            E.set<std::string>("name", "elem0");
            E.set<std::string>("type", "drift");
            E.set<double>("length", 1.0);

            elements.push_back(E);
        }

        {
            Config E;
            E.set<std::string>("name", "elem1");
            E.set<std::string>("type", "dipole");
            E.set<double>("angle", 1e-3); // 1 mrad
            E.set<double>("radius", 1.0);

            elements.push_back(E);
        }

        {
            Config E;
            E.set<std::string>("name", "elem2");
            E.set<std::string>("type", "drift");
            E.set<double>("length", 1.0);

            elements.push_back(E);
        }

        conf.setAny("elements", elements);
    }

    std::cout << "Config:\n" << conf;

    {
        std::cout << "Propogate state vector\n";

        Machine sim(conf);

        std::cout << "Machine\n"<<sim<<"\n";

        Config C;
        VectorState I(C);
        I.state[0] = 1.0;
        I.state[1] = 1e-3;

        std::cout << " Initial: " << I.state << "\n";

        sim.propogate(&I);

        std::cout << " Final: " << I.state << "\n";
    }

    conf.set<std::string>("sim-type", "Linear1DTransfer");

    {
        std::cout << "\nPropogate transfer matrix\n";

        Machine sim(conf);

        Config C;
        MatrixState I(C);

        std::cout << " Initial: " << I.state << "\n";

        sim.propogate(&I);

        std::cout << " Final: " << I.state << "\n";
    }

    return 0;
}
