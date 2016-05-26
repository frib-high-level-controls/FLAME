
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
#include <scsi/h5writer.h>

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
            ("verbose,v", po::value<std::string>()->default_value("0")->value_name("NUM"),
                "Make some noise")
            ("define,D", po::value<std::vector<std::string> >()->composing()->value_name("name=val"),
                "Override variable value (\"-Dname=value\")")
            ("lattice", po::value<std::string>()->value_name("FILE"),
                "Input lattice file")
            ("max,M", po::value<std::string>()->value_name("NUM"),
                "Maximum number of elements propagate through. (default is all)")
            ("format,F", po::value<std::string>()->value_name("FMT")->default_value("txt"),
                "output format (txt or hdf5)")
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
        Factory(const strvect& fmt) :strm(&std::cout)
        {
            assert(!fmt.empty() && fmt[0]=="txt");

            for(strvect::const_iterator it=fmt.begin()+1, end=fmt.end(); it!=end; ++it)
            {
                const std::string& cmd = *it;
                if(cmd.substr(0,5)=="file=") {
                    owned_strm.reset(new std::ofstream(cmd.substr(5).c_str()));
                    strm = owned_strm.get();
                } else {
                    std::cerr<<"Warning: -F "<<fmt[0]<<" includes unknown option "<<cmd<<"\n";
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

struct H5Observer : public Observer
{
    H5StateWriter *writer;
    H5Observer(H5StateWriter *writer) : writer(writer) {}
    virtual ~H5Observer() {}

    struct Factory : public ObserverFactory
    {
        virtual ~Factory() {}
        std::auto_ptr<H5StateWriter> writer;
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

} // namespace

int main(int argc, char *argv[])
{
try {
    po::variables_map args;
    getargs(argc, argv, args);

    std::auto_ptr<Config> conf;

    size_t verb = boost::lexical_cast<size_t>(args["verbose"].as<std::string>());
    if(verb<=2)
        H5StateWriter::dontPrint();

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

    if(verb) {
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
    registerMoment2();

    std::auto_ptr<ObserverFactory> ofact;

    {
        const std::string& ofactname = args["format"].as<std::string>();
        strvect fmt(tokenize(ofactname));

        if(fmt.empty()) {
            std::cerr<<"Empty output format\n";
            exit(1);
        } else if(fmt[0]=="txt") {
            ofact.reset(new StreamObserver::Factory(fmt));
        } else if(fmt[0]=="hdf5") {
            ofact.reset(new H5Observer::Factory(fmt));
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

    if(verb) {
        sim.set_trace(&std::cout);

        std::cout<<"# Machine configuration\n"<<sim<<"\n\n";
    }

    std::auto_ptr<StateBase> state(sim.allocState());
    sim.propagate(state.get(), 0, maxelem);

    ofact->after_sim(sim);

    if(verb) {
        std::cout << "\n# Final " << *state << "\n";
    }

    Machine::registeryCleanup();

    return 0;
}catch(std::exception& e){
    std::cerr<<"Error "<<typeid(e).name()<<" : "<<e.what()<<"\n";
    return 1;
}
}
