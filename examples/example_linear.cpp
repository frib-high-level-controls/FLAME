
#include <iostream>
#include <list>

#include <boost/numeric/ublas/io.hpp>

#include <scsi/base.h>
#include <scsi/state/vector.h>
#include <scsi/state/matrix.h>

int main(int argc, char *argv[]) {

    registerLinear(); // register linear transfer matrix states and elements

    Config conf;

    std::cout << "Setup a pair of drifts\n";
    {
        std::cout << "Building description\n";
        std::list<boost::any> elements;

        {
            Config E;
            E.set<std::string>("name", "drift1");
            E.set<std::string>("type", "drift");
            E.set<double>("length", 1.0);

            elements.push_back(E);
        }

        {
            Config E;
            E.set<std::string>("name", "drift2");
            E.set<std::string>("type", "drift");
            E.set<double>("length", 1.0);

            elements.push_back(E);
        }

        conf.setAny("elements", elements);
    }

    // Start with 3d (size()==6) vector state
    conf.set<std::string>("sim-type", "Vector");

    std::cout << "Config:\n" << conf;

    {
        std::cout << "Propogate state vector\n";

        Machine sim(conf);

        std::cout << "Machine\n"<<sim<<"\n";

        Config C;
        std::auto_ptr<StateBase> state;
        {
            std::vector<double> I(6, 0.0);
            I[0] = 1.0;
            I[1] = 1e-3;
            I[2] = 1.0;
            I[3] = 1e-3;
            I[4] = 1.0;
            I[5] = 1e-3;
            C.setAny("initial", I);
            state.reset(sim.allocState(C));
        }

        std::cout << " Initial: " << *state << "\n";

        sim.propogate(state.get());

        std::cout << " Final: " << *state << "\n";
    }

    // Switch to 3d (6x6) transfer matrix
    conf.set<std::string>("sim-type", "TransferMatrix");

    {
        std::cout << "\nPropogate transfer matrix\n";

        Machine sim(conf);

        Config C;
        // initial value is identity matrix
        std::auto_ptr<StateBase> state(sim.allocState(C));

        std::cout << " Initial: " << *state << "\n";

        sim.propogate(state.get());

        std::cout << " Final: " << *state << "\n";
    }

    std::cout << "Setup a quad triplet\n";

    {
        std::cout << "Building description\n";
        std::list<boost::any> elements;

        {
            Config E;
            E.set<std::string>("name", "quad1");
            E.set<std::string>("type", "quad");
            E.set<double>("strength", 1e-3);
            E.set<double>("length", 1.0);

            elements.push_back(E);
        }

        {
            Config E;
            E.set<std::string>("name", "drift1");
            E.set<std::string>("type", "drift");
            E.set<double>("length", 1.0);

            elements.push_back(E);
        }

        {
            Config E;
            E.set<std::string>("name", "quad2");
            E.set<std::string>("type", "quad");
            E.set<double>("strength", -2e-3);
            E.set<double>("length", 1.0);

            elements.push_back(E);
        }

        {
            Config E;
            E.set<std::string>("name", "drift2");
            E.set<std::string>("type", "drift");
            E.set<double>("length", 1.0);

            elements.push_back(E);
        }

        {
            Config E;
            E.set<std::string>("name", "quad3");
            E.set<std::string>("type", "quad");
            E.set<double>("strength", 1e-3);
            E.set<double>("length", 1.0);

            elements.push_back(E);
        }

        conf.setAny("elements", elements);
    }

    std::cout << "Config:\n" << conf;

    {
        std::cout << "\nPropogate transfer matrix\n";

        Machine sim(conf);

        std::cout << "Machine:\n" << sim;

        Config C;
        // initial value is identity matrix
        std::auto_ptr<StateBase> state(sim.allocState(C));

        std::cout << " Initial: " << *state << "\n";

        sim.propogate(state.get());

        std::cout << " Final: " << *state << "\n";
    }

    return 0;
}
