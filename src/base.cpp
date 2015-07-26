
#include <list>

#include "scsi/base.h"
#include "scsi/util.h"

StateBase::~StateBase() {}

bool
Config::has(const std::string &s) const
{
    return p_props.find(s)!=p_props.end();
}


boost::any
Config::getAny(const std::string& s) const
{
    map_t::const_iterator it = p_props.find(s);
    if(it==p_props.end())
        throw key_error(s);
    return it->second;
}

boost::any
Config::getAny(const std::string& s, const boost::any &def) const
{
    map_t::const_iterator it = p_props.find(s);
    if(it==p_props.end())
        return def;
    else
        return it->second;
}

void
Config::setAny(const std::string& s, const boost::any& val)
{
    p_props[s] = val;
}

std::ostream& operator<<(std::ostream& strm, const Config& c)
{
    for(Config::map_t::const_iterator it=c.p_props.begin(), end=c.p_props.end();
        it!=end; ++it)
    {
        try{
            typedef std::list<boost::any> sub_t;
            const sub_t& sub = boost::any_cast<const sub_t&>(it->second);
            strm<<">>> "<<it->first<<"\n";
            for(sub_t::const_iterator it2=sub.begin(), end=sub.end(); it2!=end; ++it2)
            {
                try{
                    const Config& sc = boost::any_cast<const Config&>(*it2);
                    strm<<sc<<"\n";
                }catch(boost::bad_any_cast&){
                    strm<<it2->type().name()<<"\n";
                }
            }
            strm<<"<<<\n";
        }catch(boost::bad_any_cast&){
            strm<<"Name: '"<<it->first<<"' type "<<it->second.type().name()<<"\n";
        }
    }
    return strm;
}

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

Machine::Machine(Config& c)
    :p_elements()
    ,p_info()
{
    std::string type(c.get<std::string>("sim-type"));

    p_state_infos_t::iterator it = p_state_infos.find(type);
    if(it==p_state_infos.end())
        throw key_error(type);

    p_info = it->second;

    typedef std::list<boost::any> elements_t;
    elements_t Es(c.get<elements_t>("elements"));

    p_elements_t result;
    result.reserve(Es.size());

    size_t idx=0;
    for(elements_t::iterator it=Es.begin(), end=Es.end(); it!=end; ++it)
    {
        Config EC = boost::any_cast<Config>(*it);

        std::string etype(EC.get<std::string>("type"));

        state_info::elements_t::iterator eit = p_info.elements.find(etype);
        if(eit==p_info.elements.end())
            throw key_error(etype);

        element_builder_t builder = eit->second;

        ElementVoid *E = (*builder)(EC);

        *const_cast<size_t*>(&E->index) = idx++; // ugly

        result.push_back(E);
    }

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
Machine::propogate(StateBase* S, size_t max) const
{
    const size_t nelem = p_elements.size();

    for(size_t i=0; S->next_elem<nelem && i<max; i++)
    {
        ElementVoid* E = p_elements[S->next_elem];
        S->next_elem++;
        E->advance(*S);
    }
}

StateBase*
Machine::allocState(Config& c) const
{
    return (*p_info.builder)(c);
}

Machine::p_state_infos_t Machine::p_state_infos;

void Machine::p_registerState(const char *name, state_builder_t b)
{
    if(p_state_infos.find(name)!=p_state_infos.end())
        throw key_error(name);
    state_info I;
    I.name = name;
    I.builder = b;
    p_state_infos[name] = I;
}

void Machine::p_registerElement(const std::string& sname, const char *ename, element_builder_t b)
{
    p_state_infos_t::iterator it = p_state_infos.find(sname);
    if(it==p_state_infos.end())
        throw key_error(sname);
    state_info& I = it->second;
    if(I.elements.find(ename)!=I.elements.end())
        throw key_error(ename);
    I.elements[ename] = b;
}

std::ostream& operator<<(std::ostream& strm, const Machine& m)
{
    strm<<"sim-type: "<<m.p_info.name<<"\n#Elements: "<<m.p_elements.size()<<"\n";
    for(Machine::p_elements_t::const_iterator it=m.p_elements.begin(),
        end=m.p_elements.end(); it!=end; ++it)
    {
        (*it)->show(strm);
    }
    return strm;
}
