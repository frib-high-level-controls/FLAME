
#include <iostream>
#include <list>

#include <boost/numeric/ublas/io.hpp>

#include <scsi/base.h>
#include <scsi/state/vector.h>
#include <scsi/state/matrix.h>

int main(int argc, char *argv[])
{
    FILE *in = stdin;
    if(argc>1) {
        in = fopen(argv[1], "r");
        if(!in) {
            fprintf(stderr, "Failed to open %s\n", argv[1]);
            return 2;
        }
    }

    std::auto_ptr<Config> conf;
    try {
        GLPSParser P;
        conf.reset(P.parse(in));
    }catch(std::exception& e){
        std::cerr<<"Parse error: "<<e.what()<<"\n";
        fclose(in);
        return 1;
    }

    std::cout<<"# Reduced lattice\n";
    GLPSPrint(std::cout, *conf);
    std::cout<<"\n";

    // register state and element types
    registerLinear();
    registerMoment();
    registerMoment2();

    try {
        Machine sim(*conf);
        sim.set_trace(&std::cout);

        std::cout<<"# Machine configuration\n"<<sim<<"\n\n";

        Config D;
        std::auto_ptr<StateBase> state(sim.allocState(D));
        sim.propagate(state.get());

        std::cout << "\n# Final " << *state << "\n";
    }catch(std::exception& e){
        std::cerr<<"Simulation error: "<<e.what()<<"\n";
        fclose(in);
        return 1;
    }

    fclose(in);

    return 0;
}
