
#include <sstream>
#include <set>

#include <scsi/config.h>
#include <scsi/util.h>

#include "glps_parser.h"

Config::Config()
    :value_scopes(1)
{
    value_scopes[0].reset(new values_t);
}

Config::Config(const values_t& V)
    :value_scopes(1)
{
    value_scopes[0].reset(new values_t(V));
}

Config::Config(const Config& O)
    :value_scopes(O.value_scopes)
{}

Config&
Config::operator=(const Config& O)
{
    value_scopes = O.value_scopes;
    return *this;
}

//! Ensure we have a unique reference to our inner scope, making a copy if necessary
void Config::_cow()
{
    assert(!value_scopes.empty());
    if(!value_scopes.back().unique()) {
        Config::values_pointer U(new values_t(*value_scopes.back())); // copy
        value_scopes.back().swap(U);
    }
}

const Config::value_t&
Config::getAny(const std::string& name) const
{
    assert(!value_scopes.empty());
    for(values_scope_t::const_reverse_iterator it = value_scopes.rbegin(), end = value_scopes.rend()
        ; it!=end; ++it)
    {
        values_t::const_iterator S=(*it)->find(name);
        if(S!=(*it)->end()) return S->second;
    }
    throw key_error(name);
}

void
Config::setAny(const std::string& name, const value_t& val)
{
    _cow();
    (*value_scopes.back())[name] = val;
}

void
Config::swapAny(const std::string& name, value_t& val)
{
    _cow();
    {
        values_t::iterator it = value_scopes.back()->find(name);
        if(it!=value_scopes.back()->end()) {
            it->second.swap(val);
            return;
        }
    }
    std::pair<values_t::iterator, bool> ret = value_scopes.back()->insert(std::make_pair(name,value_t()));
    assert(ret.second);
    ret.first->second.swap(val);
}

void
Config::flatten()
{
    if(depth()<=1) return;

    values_scope_t replace(1);
    // copy inner scope
    replace[0].reset(new values_t(*value_scopes.back()));
    values_t& B(*replace[0]);

    // copy enclosing scopes into new base scope
    // outer most is last to take lowest priority
    for(size_t i=value_scopes.size()-1; i; i--) {
        const values_t& S(*value_scopes[i-1]);
        B.insert(S.begin(), S.end()); // note that insert() will not overwrite existing keys
    }

    value_scopes.swap(replace);
}

size_t
Config::depth() const
{
    return value_scopes.size();
}

void
Config::push_scope()
{
    values_pointer N(new values_t);
    value_scopes.push_back(N);
}

void
Config::pop_scope()
{
    if(value_scopes.size()==1) {
        // when last scope is popped, just clear
        values_pointer N(new values_t);
        value_scopes.back().swap(N);
    } else {
        value_scopes.pop_back();
    }
}

Config Config::new_scope() const
{
    values_scope_t S(value_scopes.size()+1);
    std::copy(value_scopes.begin(), value_scopes.end(), S.begin());
    S.back().reset(new values_t);
    return Config(S);
}

namespace {
struct show_value : public boost::static_visitor<void>
{
    unsigned indent;
    std::ostream& strm;
    const std::string& name;
    show_value(std::ostream& s, const std::string& n, unsigned ind=0)
        : indent(ind), strm(s), name(n) {}

    void operator()(double v) const
    {
        unsigned i=indent;
        while(i--)
            strm.put(' ');
        strm << name << " = " << v << "\n";
    }

    void operator()(const std::string& v) const
    {
        doindent();
        strm << name << " = \"" << v << "\"\n";
    }

    void operator()(const std::vector<double>& v) const
    {
        doindent();
        strm << name << " = [";
        for(size_t i=0, N=v.size(); i<N; i++)
        {
            if(i!=0)
                strm << ", ";
            strm << v[i];
        }
        strm << "]\n";
    }

    void operator()(const Config::vector_t& v) const
    {
        doindent();
        strm << name << " = [\n";
        for(size_t i=0, N=v.size(); i<N; i++)
        {
            doindent(2);
            strm << "[" << i << "] = {\n";
            v[i].show(strm, indent+4);
            doindent(2);
            strm << "},\n";
        }
        doindent();
        strm << "]\n";
    }

    void doindent(unsigned extra=0) const
    {
        unsigned i=indent+extra;
        while(i--)
            strm.put(' ');
    }
};
}

void
Config::show(std::ostream& strm, unsigned indent) const
{
    //TODO: show nested scopes?
    for(Config::values_t::const_iterator it=value_scopes.back()->begin(), end=value_scopes.back()->end();
        it!=end; ++it)
    {
        boost::apply_visitor(show_value(strm, it->first, indent), it->second);
    }
}

namespace {
// store variable definitions in parser context
struct store_ctxt_var : public boost::static_visitor<void>
{
    const std::string& name;
    parse_context& ctxt;
    store_ctxt_var(parse_context& c, const std::string& n)
        :name(n), ctxt(c)
    {}
#define VISIT(TYPE, E) \
    void operator()(TYPE v) const { \
        ctxt.vars.push_back(parse_var(name.c_str(), expr_t(E, v))); \
        ctxt.var_idx[name] = ctxt.vars.size()-1; \
    }
    VISIT(double, glps_expr_number)
    VISIT(const std::string&, glps_expr_string)
    VISIT(const std::vector<double>&, glps_expr_vector)
#undef VISIT
    void operator()(const Config::vector_t&) const
    {
        // ignore
    }
};

void assign_expr_to_Config(Config& conf, const std::string& name, const expr_t& expr)
{
    switch(expr.etype)
    {
    case glps_expr_number:
        conf.set<double>(name, boost::get<double>(expr.value));
        break;
    case glps_expr_string:
        conf.set<std::string>(name, boost::get<std::string>(expr.value));
        break;
    case glps_expr_vector:
        conf.set<std::vector<double> >(name, boost::get<std::vector<double> >(expr.value));
        break;
    default:
        throw std::logic_error("Context contained unresolved/illegal variable");
    }
}
}

struct GLPSParser::Pvt {
    typedef Config::values_t values_t;
    values_t vars;

    void fill_vars(parse_context& ctxt)
    {
        for(values_t::const_iterator it=vars.begin(), end=vars.end(); it!=end; ++it)
        {
            // fill in ctxt.vars and ctxt.var_idx
            boost::apply_visitor(store_ctxt_var(ctxt, it->first), it->second);
        }
    }

    Config* fill_context(parse_context& ctxt)
    {
        std::auto_ptr<Config> ret(new Config);
        ret->reserve(ctxt.vars.size()+2);

        // copy ctxt.vars to top level Config
        for(parse_context::vars_t::iterator it=ctxt.vars.begin(), end=ctxt.vars.end();
            it!=end; ++it)
        {
            assign_expr_to_Config(*ret, it->name, it->expr);
        }

        if(ctxt.line.size()==0)
            throw std::runtime_error("No beamlines defined by this file");

        parse_line *line = NULL;

        {
            // find the magic "USE" element.  eg "USE: linename;"
            parse_context::map_idx_t::const_iterator it=ctxt.element_idx.find("USE");
            if(it!=ctxt.element_idx.end()) {
                parse_element &elem = ctxt.elements[it->second];
                parse_context::map_idx_t::const_iterator lit = ctxt.line_idx.find(elem.etype);

                if(lit!=ctxt.line_idx.end()) {
                    line = &ctxt.line[lit->second];
                } else {
                    std::ostringstream strm;
                    strm<<"\"USE: "<<elem.etype<<";\" references undefined beamline";
                    throw std::runtime_error(strm.str());
                }
            } else {
                // no magic USE, default to last line
                line = &ctxt.line.back();
            }
        }

        assert(line);

        if(line->names.size()==0) {
            std::ostringstream strm;
            strm<<"Beamline '"<<line->label<<"' has no elements";
            throw std::runtime_error(strm.str());
        }

        Config::vector_t elements;
        elements.resize(line->names.size());

        // copy in elements
        size_t i = 0;
        for(strlist_t::list_t::const_iterator it=line->names.begin(), end=line->names.end();
            it!=end; ++it)
        {
            Config next(ret->new_scope()); // inheirt global scope
            const parse_element& elem = ctxt.elements[ctxt.element_idx[*it]];

            next.reserve(elem.props.size()+2);

            // push elements properties
            for(kvlist_t::map_t::const_iterator itx=elem.props.begin(), endx=elem.props.end();
                itx!=endx; ++itx)
            {
                assign_expr_to_Config(next, itx->first, itx->second);
            }

            // special properties
            assert(!elem.etype.empty() && !elem.label.empty());
            next.set<std::string>("type", elem.etype);
            next.set<std::string>("name", elem.label);
            elements[i++].swap(next);
        }

        ret->swap<std::string>("name", line->label);
        ret->swap<Config::vector_t>("elements", elements);

        return ret.release();
    }
};

GLPSParser::GLPSParser()
    :priv(new Pvt)
{}

GLPSParser::~GLPSParser() {}

void
GLPSParser::setVar(const std::string& name, const Config::value_t& v)
{
    priv->vars[name] = v;
}

Config*
GLPSParser::parse_file(const char *fname)
{
    boost::filesystem::path fpath(fname);
    fpath = boost::filesystem::canonical(fpath).parent_path();

    FILE *fp;
    bool closeme = fname!=NULL && strcmp(fname,"-")!=0;
    if(closeme)
        fp = fopen(fname, "r");
    else
        fp = stdin;
    if(!fp) {
        std::ostringstream strm;
        strm<<"Failed to open file for parsing '"<<fname<<"'";
        throw std::runtime_error(strm.str());
    }
    try{
        Config *ret = parse_file(fp, fpath.native().c_str());
        if(closeme) fclose(fp);
        return ret;
    }catch(...){
        if(closeme) fclose(fp);
        throw;
    }
}

Config*
GLPSParser::parse_file(FILE *fp, const char *path)
{
    parse_context ctxt(path);
    priv->fill_vars(ctxt);
    ctxt.parse(fp);
    return priv->fill_context(ctxt);
}

Config*
GLPSParser::parse_byte(const char* s, size_t len, const char *path)
{
    parse_context ctxt(path);
    priv->fill_vars(ctxt);
    ctxt.parse(s, len);
    return priv->fill_context(ctxt);
}

Config*
GLPSParser::parse_byte(const std::string& s, const char *path)
{
    parse_context ctxt(path);
    priv->fill_vars(ctxt);
    ctxt.parse(s);
    return priv->fill_context(ctxt);
}

namespace {
// show the properties of a GLPS element
struct glps_show_props : public boost::static_visitor<void>
{
    std::ostream& strm;
    const std::string& name;
    glps_show_props(std::ostream& s, const std::string& n) :strm(s), name(n) {}

    void operator()(double v) const
    {
        strm<<", "<<name<<" = "<<v;
    }

    void operator()(const std::string& v) const
    {
        strm<<", "<<name<<" = \""<<v<<"\"";
    }

    void operator()(const std::vector<double>& v) const
    {
        strm <<", " << name << " = [";
        for(size_t i=0, N=v.size(); i<N; i++)
        {
            if(i!=0)
                strm << ", ";
            strm << v[i];
        }
        strm << "]";
    }

    void operator()(const Config::vector_t& v) const
    {
        // ignore
    }
};
// show the base GLPS Config (variables and the elements array)
struct glps_show : public boost::static_visitor<void>
{
    std::ostream& strm;
    const std::string& name;
    glps_show(std::ostream& s, const std::string& n) :strm(s), name(n) {}

    void operator()(double v) const
    {
        strm<<name<<" = "<<v<<";\n";
    }

    void operator()(const std::string& v) const
    {
        strm<<name<<" = \""<<v<<"\";\n";
    }

    void operator()(const std::vector<double>& v) const
    {
        strm << name << " = [";
        for(size_t i=0, N=v.size(); i<N; i++)
        {
            if(i!=0)
                strm << ", ";
            strm << v[i];
        }
        strm << "];\n";
    }

    void operator()(const Config::vector_t& v) const
    {
        if(name!="elements") {
            // The GLPS format Only supports nested beamline definitions
            strm << "# "<<name<<" = [... skipped ...];\n";
            return;
        }
    }
};
}

void GLPSPrint(std::ostream& strm, const Config& conf)
{
    // print variables
    for(Config::const_iterator it=conf.begin(), end=conf.end();
        it!=end; ++it)
    {
        boost::apply_visitor(glps_show(strm, it->first), it->second);
    }

    const Config::vector_t *v;
    try{
        v = &conf.get<Config::vector_t>("elements");
    }catch(key_error&){
        strm<<"# Missing beamline element list\n";
        return;
    }catch(boost::bad_get&){
        strm<<"# 'elements' is not a beamline element list\n";
        return;
    }

    std::vector<std::string> line;
    line.reserve(v->size());

    std::set<std::string> eshown;

    // print element definitions
    for(Config::vector_t::const_iterator it=v->begin(), end=v->end();
        it!=end; ++it)
    {
        bool ok = true;
        try {
            const std::string& name=it->get<std::string>("name");
            const std::string& type=it->get<std::string>("type");
            if(name.empty() || type.empty())
                throw std::runtime_error("Element missing 'name' and/or 'type'");
            line.push_back(name);
            // only show element definition once
            if(eshown.find(name)!=eshown.end())
                continue;
            strm<<name<<": "<<type;
            eshown.insert(name);
        }catch(key_error&){
            ok=false;
        }catch(boost::bad_get&){
            ok=false;
        }
        if(!ok)
            strm<<"# <malformed element>";

        for(Config::const_iterator itx=it->begin(), endx=it->end();
            itx!=endx; ++itx)
        {
            if(itx->first=="name" || itx->first=="type")
                continue;
            boost::apply_visitor(glps_show_props(strm, itx->first), itx->second);
        }

        strm<<";\n";
    }

    std::string lname(conf.get<std::string>("name", "default"));
    strm<<lname<<": LINE = (";

    bool first=true;
    for(std::vector<std::string>::const_iterator it=line.begin(), end=line.end();
        it!=end; ++it)
    {
        if(!first)
            strm<<", ";
        first = false;
        strm<<*it;
    }

    strm<<");\nUSE: "<<lname<<";\n";
}
