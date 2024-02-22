#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <typeinfo>
#include <climits>

#include <time.h>

#include <boost/numeric/ublas/io.hpp>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/call_traits.hpp>
#include <boost/lexical_cast.hpp>

#include <flame/base.h>
#include <flame/state/vector.h>
#include <flame/state/matrix.h>
#include <flame/moment.h>

#ifdef USE_HDF5
#include <flame/h5writer.h>
#endif

namespace po = boost::program_options;

namespace {

typedef std::vector<std::string> strvect;

strvect tokenize(const std::string& inp)
{
    strvect ret;
    size_t pos = 0;
    while(true) {
        size_t sep = inp.find_first_of(',', pos);
        if(sep==inp.npos) {
            ret.push_back(inp.substr(pos));
            break;
        } else if(sep!=pos) {
            ret.push_back(inp.substr(pos, sep-pos));
        } else {
            // ignore empty
        }
        pos = sep+1;
    }
    return ret;
}

static
void getargs(int argc, char *argv[], po::variables_map& args)
{
    std::ostringstream caption;
    caption<<argv[0]<<" [options] <lattice file>";
    po::options_description opts(caption.str());
    opts.add_options()
            ("help,h", "Display this message")
            ("verbose,v", po::value<int>()->default_value(0)->value_name("NUM"),
                "Make some noise.  Values 0 through 6 are useful.")
            ("define,D", po::value<std::vector<std::string> >()->composing()->value_name("name=type:val"),
                "Define variable value (\"-Dname=str:value\")")
            ("lattice", po::value<std::string>()->value_name("FILE"),
                "Input lattice file")
            ("max,M", po::value<std::string>()->value_name("NUM"),
                "Maximum number of elements propagate through. (default is all)")
#ifdef USE_HDF5
            ("format,F", po::value<std::string>()->value_name("FMT")->default_value("txt"),
            "output format (txt or hdf5)")
#else
            ("format,F", po::value<std::string>()->value_name("FMT")->default_value("txt"),
            "output format (txt)")
#endif
            ("select-all,A", "Select all elements for output")
            ("select-type,T", po::value<std::vector<std::string> >()->composing()->value_name("ETYPE"),
                "Select all elements of the given type for output")
            ("select-name,N", po::value<std::vector<std::string> >()->composing()->value_name("ENAME"),
                "Select all elements with the given name for output")
            ("select-last,L", "Select last element for output")
#ifdef CLOCK_MONOTONIC
            ("timeit", "Measure execution time")
#endif
            ;

    po::positional_options_description pos;
    pos.add("lattice", 1);

    po::store(po::command_line_parser(argc, argv).options(opts).positional(pos).run(), args);
    po::notify(args);

    if(args.count("help") || !args.count("lattice")) {
        std::cout<<opts<<"\n\n"
                   "Output formats:\n\n"
                   " txt  - Print selected outputs states to screen ('--format txt' the default default)\n"
                   "        or file ('--format txt[,file=out.txt][,verbose[=lvl#]]')\n"
                   "\n"
                   " utest - Print selected output states to screen as a python checkPropagate()\n"
                   "\n"
#ifdef USE_HDF5
                   " hdf5 - Write selected output states to an HDF5 file.\n"
                   "        eg. '--format hdf5,file=out.h5'\n"
                   "\n"
#endif
                   "Definitions:\n\n"
                   " Variable defintions made by arguments must specify a name, type, and value\n"
                   " The type may be: 'str'' or 'double', which may be abbreviated as 'S' or 'D'.\n"
                   " Definitions are overwritten by those in the lattice file.\n"
                   ;
        exit(1);
    }
}

struct ObserverFactory
{
    virtual ~ObserverFactory() {}
    virtual Observer *observe(Machine& M, ElementVoid* E) = 0;
    virtual void before_sim(Machine&) {}
    virtual void after_sim(Machine&) {}
};

struct UnitTestObserver : public Observer
{
    typedef boost::shared_ptr<std::list<std::string> > interested_t;
    interested_t interested;
    typedef std::vector<unsigned> interested_indicies_t;
    interested_indicies_t interested_indicies;

    UnitTestObserver(const interested_t& I) :interested(I) {}
    virtual ~UnitTestObserver() {}

    void lookup(StateBase *S) {
        typedef std::map<std::string, unsigned> lookup_t;
        lookup_t L;

        unsigned idx=0;
        StateBase::ArrayInfo info;
        while(S->getArray(idx++, info)) {
            bool skip = false;
            switch(info.type) {
            case StateBase::ArrayInfo::Sizet:
                if(info.ndim!=0) skip = true;
            case StateBase::ArrayInfo::Double:
                break;
            default:
                skip = true;
            }
            if(info.ndim>2) skip=true;
            if(skip) continue;

            L[info.name] = idx-1;
        }

        interested_indicies.resize(interested->size());

        interested_t::element_type::const_iterator it = interested->begin();
        for(size_t i=0; i<interested_indicies.size(); i++, ++it) {
            lookup_t::const_iterator Lit = L.find(*it);
            if(Lit!=L.end()) {
                interested_indicies[i] = Lit->second;
            }
        }
    }

    virtual void view(const ElementVoid *elem, const StateBase *state)
    {
        // hack since getArray() is non-const
        // we won't actually modify the state
        StateBase *S = const_cast<StateBase*>(state);

        if(interested_indicies.size()!=interested->size())
            lookup(S);

        std::cout<<"    def test_"<<elem->type_name()<<"(self):\n"
                   "        # "<<elem->name<<"\n"
                 <<"        self.checkPropagate(0, {}, {\n";

        for(size_t i=0; i<interested_indicies.size(); i++) {
            unsigned idx = interested_indicies[i];
            StateBase::ArrayInfo info;
            bool valid = S->getArray(idx, info);
            assert(valid);
            (void)valid;

            std::cout<<"            '"<<info.name<<"':";
            if(info.ndim==0) {
                switch(info.type) {
                case StateBase::ArrayInfo::Double:
                    std::cout<<std::scientific << std::setprecision(17)<<*(double*)info.ptr;
                    break;
                case StateBase::ArrayInfo::Sizet:
                    std::cout<<*(size_t*)info.ptr;
                }

            } else if(info.ndim==1) {
                assert(info.type==StateBase::ArrayInfo::Double);
                std::cout<<"asfarray([";
                for(size_t i=0; i<info.dim[0]; i++) {
                    std::cout<<std::scientific << std::setprecision(16)<<*info.get<double>(&i);
                    if(i!=info.dim[0]-1)
                        std::cout<<", ";
                }
                std::cout<<"])";

            } else if(info.ndim==2) {
                assert(info.type==StateBase::ArrayInfo::Double);
                std::cout<<"asfarray([\n";
                size_t idx[StateBase::ArrayInfo::maxdims];
                memset(idx, 0, sizeof(idx));
                for(idx[0]=0; idx[0]<info.dim[0]; idx[0]++) {
                    std::cout<<"                [";
                    for(idx[1]=0; idx[1]<info.dim[1]; idx[1]++) {
                        std::cout<<std::scientific << std::setprecision(16)<<*info.get<double>(idx);
                        if(idx[1]!=info.dim[1]-1)
                            std::cout<<", ";
                    }
                    std::cout<<"],\n";
                }
                std::cout<<"            ])";
            } else {
                std::cout<<"None";
            }
            std::cout<<",\n";
        }
        std::cout<<"        }, max="<<elem->index+1<<")\n";
    }

    struct Factory : public ObserverFactory
    {
        interested_t interested;
        Factory(const strvect& fmt) :interested(new interested_t::element_type)
        {
            assert(!fmt.empty() && fmt[0]=="utest");

            for(strvect::const_iterator it=fmt.begin()+1, end=fmt.end(); it!=end; ++it)
            {
                const std::string& cmd = *it;
                if(cmd.find_first_of('=')==cmd.npos) {
                    interested->push_back(cmd);
                } else {
                    std::cerr<<"Warning: -F "<<fmt[0]<<" includes unknown option "<<cmd<<"\n";
                }
            }
        }
        virtual ~Factory() {}
        virtual Observer *observe(Machine& M, ElementVoid* E)
        {
            return new UnitTestObserver(interested);
        }
    };
};

struct StreamObserver : public Observer
{
    std::ostream *strm;
    int detail;
    StreamObserver(std::ostream& strm, int detail=1) :strm(&strm), detail(detail) {}
    virtual ~StreamObserver() {}
    virtual void view(const ElementVoid* elem, const StateBase* state)
    {
        (*strm)<<"After Element ["<<elem->index<<"] "<<elem->name<<" ";
        state->show(*strm, detail);
        (*strm)<<"\n";
    }

    struct Factory : public ObserverFactory
    {
        std::unique_ptr<std::ostream> owned_strm;
        std::ostream *strm;
        int detail;
        Factory(const strvect& fmt) :strm(&std::cout), detail(1)
        {
            assert(!fmt.empty() && fmt[0]=="txt");

            for(strvect::const_iterator it=fmt.begin()+1, end=fmt.end(); it!=end; ++it)
            {
                const std::string& cmd = *it;
                if(cmd.substr(0,5)=="file=") {
                    owned_strm.reset(new std::ofstream(cmd.substr(5).c_str()));
                    strm = owned_strm.get();
                } else if(cmd=="verbose") {
                    detail=1;
                } else if(cmd.substr(0,8)=="verbose=") {
                    detail=boost::lexical_cast<int>(cmd.substr(8));
                } else {
                    std::cerr<<"Warning: -F "<<fmt[0]<<" includes unknown option "<<cmd<<"\n";
                }
            }
        }
        virtual ~Factory() {}
        virtual Observer *observe(Machine& M, ElementVoid* E)
        {
            return new StreamObserver(*strm, detail);
        }

        virtual void after_sim(Machine&)
        {
            strm->flush();
        }
    };
};

#ifdef USE_HDF5
struct H5Observer : public Observer
{
    H5StateWriter *writer;
    H5Observer(H5StateWriter *writer) : writer(writer) {}
    virtual ~H5Observer() {}

    struct Factory : public ObserverFactory
    {
        virtual ~Factory() {}
        std::unique_ptr<H5StateWriter> writer;
        Factory(const strvect& fmt)
        {
            assert(!fmt.empty() && fmt[0]=="hdf5");

            for(strvect::const_iterator it=fmt.begin()+1, end=fmt.end(); it!=end; ++it)
            {
                const std::string& cmd = *it;
                if(cmd.substr(0,5)=="file=") {
                    writer.reset(new H5StateWriter(cmd.substr(5)));
                } else {
                    std::cerr<<"Warning: -F "<<fmt[0]<<" includes unknown option "<<cmd<<"\n";
                }
            }
            if(!writer.get()) {
                std::cerr<<"Warning: hdf5 output format requires file=...\n";
            }
        }
        virtual Observer *observe(Machine &M, ElementVoid *E)
        {
            if(!writer.get()) return NULL;
            else              return new H5Observer(writer.get());
        }
        virtual void before_sim(Machine & M)
        {
            if(writer.get()) writer->setAttr("sim_type", M.simtype());
        }
        virtual void after_sim(Machine&)
        {
            if(writer.get()) writer->close();
            writer.reset();
        }
    };

    virtual void view(const ElementVoid *, const StateBase *state)
    {
        writer->append(state);
    }
};
#endif

struct Timer {
    timespec ts;
    Timer() {
#ifdef CLOCK_MONOTONIC
        clock_gettime(CLOCK_MONOTONIC, &ts);
#endif
    }
    double delta() {
#ifdef CLOCK_MONOTONIC
        timespec start = ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);

        // tv_nsec and tv_sec are signed integers
        double D = ts.tv_nsec-start.tv_nsec;
        D *= 1e-9;
        D += ts.tv_sec-start.tv_sec;
        return D;
#else
        return std::numeric_limits<double>::quiet_NaN();
#endif
    }
    void showdelta(const char *msg) {
#ifdef CLOCK_MONOTONIC
        double D = delta();
        printf("%s : %.3f ms\n", msg, D*1e3);
#endif
    }
};

} // namespace

int main(int argc, char *argv[])
{
try {
    po::variables_map args;
    getargs(argc, argv, args);

    bool showtime = args.count("timeit")>0;
    Timer timeit;

    std::unique_ptr<Config> conf;

    int verb = args["verbose"].as<int>();
    if(verb<=2)
#ifdef USE_HDF5
        H5StateWriter::dontPrint();
#endif
    if(verb>2)
        Machine::log_detail=(FLAME_ERROR-10*(verb-2));
    {
        GLPSParser P;

        if(args.count("define")) {
            const std::vector<std::string>& defs = args["define"].as<std::vector<std::string> >();

            BOOST_FOREACH(const std::string& def, defs) {
                // expected form "<name>=<type>:<value>"
                size_t equal = def.find_first_of('='),
                       colon = def.find_first_of(':', equal);
                if(equal==def.npos) {
                    std::cerr<<"-D "<<def<<" missing '='\n";
                    exit(1);
                } else if(colon==def.npos) {
                    std::cerr<<"-D "<<def<<" missing ':'\n";
                    exit(1);
                } else if(equal==0) {
                    std::cerr<<"-D "<<def<<" missing variable name\n";
                    exit(1);
                }

                std::string name(def.substr(0,equal)),
                            type(def.substr(equal+1, colon-equal-1)),
                           value(def.substr(colon+1));

                Config::value_t curval;

                if(type=="double" || type=="D") {
                    curval = boost::lexical_cast<double>(value);

                } else if(type=="str" || type=="S") {
                    curval = value;

                } else {
                    std::cerr<<"Unknown type "<<type<<" in -D "<<def<<"\n";
                    exit(1);
                }

                P.setVar(name, curval);
            }
        }

        try {
            conf.reset(P.parse_file(args["lattice"].as<std::string>().c_str()));
        }catch(std::exception& e){
            std::cerr<<"Parse error: "<<e.what()<<"\n";
            return 1;
        }
    }

    if(showtime) timeit.showdelta("Parsing");

    if(verb) {
        std::cout<<"# Reduced lattice\n";
        GLPSPrint(std::cout, *conf);
        std::cout<<"\n";
    }

    //size_t maxelem = (size_t)-1;
    //if(args.count("max")) {
    //    maxelem = boost::lexical_cast<size_t>(args["max"].as<std::string>());
    //}
    int maxelem = INT_MAX;
    if(args.count("max")) {
        maxelem = boost::lexical_cast<int>(args["max"].as<std::string>());
    }

    // register state and element types
    registerLinear();
    registerMoment();

    std::unique_ptr<ObserverFactory> ofact;

    {
        const std::string& ofactname = args["format"].as<std::string>();
        strvect fmt(tokenize(ofactname));

        if(fmt.empty()) {
            std::cerr<<"Empty output format\n";
            exit(1);
        } else if(fmt[0]=="txt") {
            ofact.reset(new StreamObserver::Factory(fmt));
#ifdef USE_HDF5
        } else if(fmt[0]=="hdf5") {
            ofact.reset(new H5Observer::Factory(fmt));
#endif
        } else if(fmt[0]=="utest") {
            ofact.reset(new UnitTestObserver::Factory(fmt));
        } else {
            std::cerr<<"Unknown output format \""<<ofactname<<"\"\n";
            exit(1);
        }
    }

    if(showtime) timeit.showdelta("Setup 1");

    Machine sim(*conf);

    if(showtime) timeit.showdelta("Create Machine");

    if(args.count("select-all")) {
        BOOST_FOREACH(ElementVoid *elem, sim) {
            assert(elem->observer()==NULL);
            elem->set_observer(ofact->observe(sim, elem));
        }
    }
    if(args.count("select-type")) {
        BOOST_FOREACH(const std::string& etype, args["select-type"].as<std::vector<std::string> >()) {

            std::pair<Machine::lookup_iterator, Machine::lookup_iterator> S(sim.equal_range_type(etype));

            if(S.first==S.second) {
                std::cerr<<"Warning: --select-type "<<etype<<" does not match any elements\n";
            } else {
                for(; S.first!=S.second; ++S.first) {
                    ElementVoid *elem = *S.first;

                    if(elem->observer()==NULL) {
                        // don't replace existing Observer
                        elem->set_observer(ofact->observe(sim, elem));
                    }
                }
            }
        }
    }
    if(args.count("select-name")) {
        BOOST_FOREACH(const std::string& ename, args["select-name"].as<std::vector<std::string> >()) {

            std::pair<Machine::lookup_iterator, Machine::lookup_iterator> S(sim.equal_range(ename));

            if(S.first==S.second) {
                std::cerr<<"Warning: --select-name "<<ename<<" does not match any elements\n";
            } else {
                for(; S.first!=S.second; ++S.first) {
                    ElementVoid *elem = *S.first;

                    if(elem->observer()==NULL) {
                        // don't replace existing Observer
                        elem->set_observer(ofact->observe(sim, elem));
                    }
                }
            }
        }
    }
    if(args.count("select-last") && sim.size()>0) {
        ElementVoid *elem = sim[sim.size()-1];

        if(elem->observer()==NULL) {
            // don't replace existing Observer
            elem->set_observer(ofact->observe(sim, elem));
        }
    }

    ofact->before_sim(sim);

    if(verb) {
        sim.set_trace(&std::cout);

        std::cout<<"# Machine configuration\n"<<sim<<"\n\n";
    }

    if(showtime) timeit.showdelta("Setup 2");

    std::unique_ptr<StateBase> state(sim.allocState());
    if(showtime) timeit.showdelta("Alloc State");
    sim.propagate(state.get(), 0, maxelem);
    if(showtime) {
        timeit.showdelta("Simulate (cache cold)");
        sim.propagate(state.get(), 0, maxelem);
        timeit.showdelta("Simulate (cache hot)");
    }

    ofact->after_sim(sim);

    if(verb) {
        std::cout << "\n# Final " << *state << "\n";
    }

    Machine::registeryCleanup();
    if(showtime) timeit.showdelta("Cleanup");

    return 0;
}catch(std::exception& e){
    std::cerr<<"Error "<<typeid(e).name()<<" : "<<e.what()<<"\n";
    return 1;
}
}
