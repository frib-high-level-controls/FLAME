#include <stdio.h>
#include <iostream>

//! [includes]
#include "scsi/config.h"
#include "scsi/base.h"
//! [includes]

int main(int argc, char *argv[])
{
    try {
        if(argc<2) return 2;

        //! [Register sim_types]
        registerLinear();
        //! [Register sim_types]

        //! [Parse lattice file]
        FILE *fp = fopen(argv[1], "r");

        std::auto_ptr<Config> conf;
        try{
            GLPSParser parser;
            conf.reset(parser.parse(fp));
        }catch(...){
            fclose(fp);
            throw;
        }
        fclose(fp);
        //! [Parse lattice file]

        //! [Construct Machine]
        Machine mymachine(*conf);
        //! [Construct Machine]

        mymachine.set_trace(&std::cout); // print intermediates

        //! [Allocate State]
        std::auto_ptr<StateBase> thestate(mymachine.allocState());
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
