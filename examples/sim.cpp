#include <stdio.h>
#include <iostream>

//! [includes]
#include "flame/config.h"
#include "flame/base.h"
//! [includes]

int main(int argc, char *argv[])
{
    try {
        if(argc<2) return 2;

        //! [Register sim_types]
        registerLinear();
        //! [Register sim_types]

        //! [Parse lattice file]
        std::unique_ptr<Config> conf;
        {
            GLPSParser parser;
            conf.reset(parser.parse_file(argv[1]));
        }
        //! [Parse lattice file]

        //! [Construct Machine]
        Machine mymachine(*conf);
        //! [Construct Machine]

        mymachine.set_trace(&std::cout); // print intermediates

        //! [Allocate State]
        std::unique_ptr<StateBase> thestate(mymachine.allocState());
        //! [Allocate State]

        std::cout<<"Initial state "<<*thestate<<"\n";

        //! [Run simulation]
        mymachine.propagate(thestate.get());
        //! [Run simulation]

        std::cout<<"Final state "<<*thestate<<"\n";

        //! [Deinit]
        Machine::registeryCleanup();
        //! [Deinit]

        return 0;
    } catch(std::exception& e) {
        std::cerr<<"Error: "<<e.what()<<"\n";
        return 1;
    }
}
