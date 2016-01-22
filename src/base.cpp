
#include <list>
#include <sstream>

#include <boost/thread/mutex.hpp>

#include "scsi/base.h"
#include "scsi/util.h"

namespace {
// This mutex guards the global Machine::p_state_infos
typedef boost::mutex info_mutex_t;
info_mutex_t info_mutex;
}

StateBase::~StateBase() {}

StateBase::StateBase(const Config& c)
    :next_elem(0)
    ,IonZ(c.get<double >("IonZ", 0))
    ,IonEs(c.get<double>("IonEs", 0))
    ,IonEk(c.get<double>("IonEk", 0))
    ,IonW(c.get<double >("IonW", 0))
    ,pyptr(0)
{}

ElementVoid::ElementVoid(const Config& conf)
    :name(conf.get<std::string>("name"))
    ,index(0)
    ,p_conf(conf)
{}

ElementVoid::~ElementVoid() {}

void ElementVoid::show(std::ostream& strm) const
{
    strm<<"Element "<<index<<": "<<name<<" ("<<type_name()<<")\n";
}

Machine::Machine(const Config& c)
    :p_elements()
    ,p_trace(NULL)
    ,p_info()
{
    std::string type(c.get<std::string>("sim_type"));

    info_mutex_t::scoped_lock G(info_mutex);

    p_state_infos_t::iterator it = p_state_infos.find(type);
    if(it==p_state_infos.end()) {
        std::ostringstream msg;
        msg<<"Unsupport sim_type '"<<type<<"'";
        throw key_error(msg.str());
    }

    p_info = it->second;

    typedef Config::vector_t elements_t;
    elements_t Es(c.get<elements_t>("elements"));

    p_elements_t result;
    result.reserve(Es.size());

    size_t idx=0;
    for(elements_t::iterator it=Es.begin(), end=Es.end(); it!=end; ++it)
    {
        const Config& EC = *it;

        const std::string& etype(EC.get<std::string>("type"));

        state_info::elements_t::iterator eit = p_info.elements.find(etype);
        if(eit==p_info.elements.end())
            throw key_error(etype);

        element_builder_t builder = eit->second;

        ElementVoid *E;
        try{
            E = (*builder)(EC);
        }catch(key_error& e){
            std::ostringstream strm;
            strm<<"Error while initializing element "<<idx<<" '"<<EC.get<std::string>("name", "<invalid>")
               <<"' : missing required parameter '"<<e.what()<<"'";
            throw key_error(strm.str());

        }catch(std::exception& e){
            std::ostringstream strm;
            strm<<"Error while constructing element "<<idx<<" '"<<EC.get<std::string>("name", "<invalid>")
               <<"' : "<<e.what();
            throw std::runtime_error(strm.str());
        }

        *const_cast<size_t*>(&E->index) = idx++; // ugly

        result.push_back(E);
    }

    G.unlock();

    for(p_elements_t::iterator it=p_elements.begin(), end=p_elements.end(); it!=end; ++it)
    {
        (*it)->peek(result);
    }

    p_elements.swap(result);
}

Machine::~Machine()
{
    for(p_elements_t::iterator it=p_elements.begin(), end=p_elements.end(); it!=end; ++it)
    {
        delete *it;
    }
}

void
Machine::propagate(StateBase* S, size_t start, size_t max) const
{
    const size_t nelem = p_elements.size();

    S->next_elem = start;
    for(size_t i=0; S->next_elem<nelem && i<max; i++)
    {
        ElementVoid* E = p_elements[S->next_elem];
        S->next_elem++;
        E->advance(*S);
        if(p_trace)
            (*p_trace) << "After "<< i<< " " << *S;
    }
}

StateBase*
Machine::allocState(Config& c) const
{
    return (*p_info.builder)(c);
}

void Machine::reconfigure(size_t idx, const Config& c)
{
    if(idx>=p_elements.size())
        throw std::invalid_argument("element index out of range");

    const std::string& etype(c.get<std::string>("type"));

    state_info::elements_t::iterator eit = p_info.elements.find(etype);
    if(eit==p_info.elements.end())
        throw key_error(etype);

    element_builder_t builder = eit->second;

    std::auto_ptr<ElementVoid> E((*builder)(c));

    *const_cast<size_t*>(&E->index) = idx; // ugly

    delete p_elements[idx];
    p_elements[idx] = E.release();
}

Machine::p_state_infos_t Machine::p_state_infos;

void Machine::p_registerState(const char *name, state_builder_t b)
{
    info_mutex_t::scoped_lock G(info_mutex);
    if(p_state_infos.find(name)!=p_state_infos.end()) {
        std::ostringstream strm;
        strm<<"attempt to register already registered sim_type=\""<<name<<"\"";
        throw std::logic_error(strm.str());
    }
    state_info I;
    I.name = name;
    I.builder = b;
    p_state_infos[name] = I;
}

void Machine::p_registerElement(const std::string& sname, const char *ename, element_builder_t b)
{
    info_mutex_t::scoped_lock G(info_mutex);
    p_state_infos_t::iterator it = p_state_infos.find(sname);
    if(it==p_state_infos.end()) {
        std::ostringstream strm;
        strm<<"can't add element \""<<ename<<"\" for unknown sim_type=\""<<sname<<"\"";
        throw std::logic_error(strm.str());
    }
    state_info& I = it->second;
    if(I.elements.find(ename)!=I.elements.end()) {
        std::ostringstream strm;
        strm<<"element type \""<<ename<<"\" has already been registered for "
              "sim_type=\""<<sname<<"\"";
        throw std::logic_error(strm.str());
    }
    I.elements[ename] = b;
}

std::ostream& operator<<(std::ostream& strm, const Machine& m)
{
    strm<<"sim_type: "<<m.p_info.name<<"\n#Elements: "<<m.p_elements.size()<<"\n";
    for(Machine::p_elements_t::const_iterator it=m.p_elements.begin(),
        end=m.p_elements.end(); it!=end; ++it)
    {
        (*it)->show(strm);
    }
    return strm;
}
