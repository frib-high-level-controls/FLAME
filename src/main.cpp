
#include <iostream>
#include <stdio.h>

#include "scsi/config.h"

#include <scsi/base.h>
#include <scsi/state/vector.h>
#include <scsi/state/matrix.h>

extern int glps_debug;


void prt_lat(Machine &sim)
{
    std::stringstream strm;

    std::cout << "# Machine configuration\n" << sim << "\n\n";

    strm << "sim_type: " << sim.p_info.name << "\n#Elements: " << sim.p_elements.size() << "\n";
    for(Machine::p_elements_t::const_iterator it = sim.p_elements.begin(),
        end = sim.p_elements.end(); it != end; ++it)
    {
        (*it)->show(strm);
    }
    std::cout << strm;
}


//Machine::propagate_jb(StateBase* S, size_t start, size_t max) const
//{
//    const size_t nelem = p_elements.size();

//    for(size_t i=start; S->next_elem<nelem && i<max; i++)
//    {
//        ElementVoid* E = p_elements[S->next_elem];
//        S->next_elem++;
//        E->advance(*S);
//        if(p_trace)
//            (*p_trace) << "After "<< i << " " << *S;
//    }
//}

void initialize_long(Machine &sim)
{
//    typedef ElementRFCavity<MomentElementBase> rfcav_t;

    for(Machine::p_elements_t::iterator it = sim.p_elements.begin(),
        end = sim.p_elements.end(); it != end; ++it)
    {
        ElementVoid* elem = *it;
        const Config & conf = elem->conf();
        std::string t_name = elem->type_name(); // C string -> C++ string.

        std::cout << t_name << std::endl;

        if (t_name == "marker") {
        } else if (t_name == "drift") {
            std::cout << "drift: L =" << conf.get<double>("L") << std::endl;
        } else if (t_name == "sbend") {
        } else if (t_name == "quadrupole") {
        } else if (t_name == "solenoid") {
        } else if (t_name == "rfcavity") {
            std::cout << "cavity: L =" << conf.get<double>("L")
                      << ", phi = " << conf.get<double>("phi") << std::endl;

//            rfcav_t* rfcav_ptr = dynamic_cast<rfcav_t*>(elem);
//            assert(rfcav_ptr != NULL); // Check if cavity.
//            rfcav_ptr-> = phi;
        }
    }
}


int main(int argc, char *argv[])
{
    Machine::p_elements_t::const_iterator it;

    FILE *in = stdin;
    if(argc>1) {
        in = fopen(argv[1], "r");
        if(!in) {
            fprintf(stderr, "Failed to open %s\n", argv[1]);
            return 2;
        }
    }

    glps_debug = 0; // 0 or 1.

    std::auto_ptr<Config> conf;

    try {
        GLPSParser P;
        conf.reset(P.parse(in));
        fprintf(stderr, "Parsing succeeds\n");
    } catch(std::exception& e) {
        fprintf(stderr, "Parse error: %s\n", e.what());
        fclose(in);
        return 1;
    }

    // register state and element types
    registerLinear();
    registerMoment();

    Machine sim(*conf);
    sim.set_trace(&std::cout);

    initialize_long(sim);

//    it = sim.p_elements.begin();

    std::cout << "\nLattice parameters:\n";
    std::cout << "sim_type: " << sim.p_info.name
              << "\n#Elements: " << sim.p_elements.size()
              << "\n";

//    prt_lat(sim);

//    std::cout<<"# Reduced lattice\n";
//    GLPSPrint(std::cout, *conf);
//    std::cout<<"\n";

//    std::cerr<<"Generic AST:\n";
//    std::cerr << *conf;
//    std::cerr<<"GLPS:\n";
//    GLPSPrint(std::cout, *conf);

    if (false) {
        try {
            Machine sim(*conf);
            sim.set_trace(&std::cout);

            std::cout << "# Machine configuration\n" << sim << "\n\n";

            Config D;
            std::auto_ptr<StateBase> state(sim.allocState(D));
            sim.propagate(state.get());

            std::cout << "\n# Final " << *state << "\n";
        } catch(std::exception& e) {
            std::cerr << "Simulation error: " << e.what() << "\n";
            fclose(in);
            return 1;
        }
    }

    fprintf(stderr, "Done\n");
    fclose(in);
    return 0;
}

