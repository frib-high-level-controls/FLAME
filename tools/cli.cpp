
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <typeinfo>

#include <boost/numeric/ublas/io.hpp>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/call_traits.hpp>
#include <boost/lexical_cast.hpp>

#include <scsi/base.h>
#include <scsi/state/vector.h>
#include <scsi/state/matrix.h>

namespace po = boost::program_options;

namespace {

std::vector<std::string> tokenize(const std::string& inp)
{
    std::vector<std::string> ret;
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
            ("verbose,v", "Make some noise")
            ("define,D", po::value<std::vector<std::string> >()->composing()->value_name("name=val"),
                "Override variable value (\"-Dname=value\")")
            ("lattice", po::value<std::string>()->value_name("FILE"),
                "Input lattice file")
            ("max,M", po::value<std::string>()->value_name("NUM"),
                "Maximum number of elements propagate through. (default is all)")
            ("format,F", po::value<std::string>()->value_name("FMT")->default_value("csv"),
                "output format")
            ("select-all,A", "Select all elements for output")
            ("select-type,T", po::value<std::string>()->value_name("ETYPE"),
                "Select all elements of the given type for output")
            ;

    po::positional_options_description pos;
    pos.add("lattice", 1);

    po::store(po::command_line_parser(argc, argv).options(opts).positional(pos).run(), args);
    po::notify(args);

    if(args.count("help") || !args.count("lattice")) {
        std::cout<<opts<<"\n";
        exit(1);
    }
}

struct update_define : boost::static_visitor<Config::value_t>
{
    const std::string& value;
    update_define(const std::string& val) : value(val) {}

    Config::value_t operator()(double v) const
    {
        return boost::lexical_cast<double>(value);
    }
    Config::value_t operator()(const std::string& v) const
    {
        return value;
    }
    Config::value_t operator()(const std::vector<double>& v) const
    {
        throw std::runtime_error("-D can't set vector values");
    }
    Config::value_t operator()(const Config::vector_t& v) const
    {
        throw std::runtime_error("-D can't set Config values");
    }
};

struct ObserverFactory
{
    virtual ~ObserverFactory() {}
    virtual Observer *observe(Machine& M, ElementVoid* E) = 0;
    virtual void before_sim(Machine&) {}
    virtual void after_sim(Machine&) {}
};

struct StreamObserver : public Observer
{
    std::ostream *strm;
    StreamObserver(std::ostream& strm) :strm(&strm) {}
    virtual ~StreamObserver() {}
    virtual void view(const ElementVoid* elem, const StateBase* state)
    {
        (*strm)<<elem->index<<" "<<*state;
    }

    struct Factory : public ObserverFactory
    {
        std::auto_ptr<std::ostream> owned_strm;
        std::ostream *strm;
        Factory(const std::string& fmt) :strm(&std::cout)
        {
            assert(fmt.substr(0,3)=="csv");
            std::cerr<<"XXX '"<<fmt.substr(3)<<"'\n";
            BOOST_FOREACH(const std::string& cmd, tokenize(fmt.substr(3)))
            {
                std::cerr<<"YYY '"<<cmd<<"'\n";
                if(cmd.substr(0,5)=="file=") {
                    owned_strm.reset(new std::ofstream(cmd.substr(5).c_str()));
                    strm = owned_strm.get();
                } else {
                    std::cerr<<"Warning: -F "<<fmt<<" includes unknown option "<<cmd<<"\n";
                }
            }
        }
        virtual ~Factory() {}
        virtual Observer *observe(Machine& M, ElementVoid* E)
        {
            return new StreamObserver(*strm);
        }

        virtual void after_sim(Machine&)
        {
            strm->flush();
        }
    };
};

} // namespace

int main(int argc, char *argv[])
{
try {
    po::variables_map args;
    getargs(argc, argv, args);

    std::auto_ptr<Config> conf;

    try {
        GLPSParser P;
        conf.reset(P.parse_file(args["lattice"].as<std::string>().c_str()));
    }catch(std::exception& e){
        std::cerr<<"Parse error: "<<e.what()<<"\n";
        return 1;
    }

    if(args.count("define")) {
        const std::vector<std::string>& defs = args["define"].as<std::vector<std::string> >();

        BOOST_FOREACH(const std::string& def, defs) {
            size_t sep = def.find_first_of('=');
            if(sep==def.npos) {
                std::cerr<<"-D "<<def<<" missing '='\n";
                exit(1);
            } else if(sep==0) {
                std::cerr<<"-D "<<def<<" missing variable name\n";
                exit(1);
            }

            std::string name(def.substr(0,sep));

            Config::value_t curval;
            try {
                curval = conf->getAny(name);
            } catch(key_error& e) {
                std::cerr<<"-D "<<def<<" variable "<<name<<" does not exist.  -D may only redefine existing variables\n";
                exit(1);
            }

            conf->setAny(name, boost::apply_visitor(update_define(def.substr(sep+1)), curval));
        }
    }

    if(args.count("verbose")) {
        std::cout<<"# Reduced lattice\n";
        GLPSPrint(std::cout, *conf);
        std::cout<<"\n";
    }

    size_t maxelem = (size_t)-1;
    if(args.count("max")) {
        maxelem = boost::lexical_cast<size_t>(args["max"].as<std::string>());
    }

    // register state and element types
    registerLinear();
    registerMoment();

    std::auto_ptr<ObserverFactory> ofact;

    {
        const std::string& ofactname = args["format"].as<std::string>();
        if(ofactname.substr(0, 3)=="csv") {
            ofact.reset(new StreamObserver::Factory(ofactname));
        } else {
            std::cerr<<"Unknown output format \""<<ofactname<<"\"\n";
            exit(1);
        }
    }


    Machine sim(*conf);

    if(args.count("select-all")) {
        BOOST_FOREACH(ElementVoid *elem, sim) {
            assert(elem->observer()==NULL);
            elem->set_observer(ofact->observe(sim, elem));
        }
    }
    if(args.count("select-type")) {
        const std::string& etype = args["select-type"].as<std::string>();

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

    ofact->before_sim(sim);

    if(args.count("verbose")) {
        sim.set_trace(&std::cout);

        std::cout<<"# Machine configuration\n"<<sim<<"\n\n";
    }

    std::auto_ptr<StateBase> state(sim.allocState());
    sim.propagate(state.get(), 0, maxelem);

    ofact->after_sim(sim);

    if(args.count("verbose")) {
        std::cout << "\n# Final " << *state << "\n";
    }

    return 0;
}catch(std::exception& e){
    std::cerr<<"Error "<<typeid(e).name()<<" : "<<e.what()<<"\n";
    return 1;
}
}
