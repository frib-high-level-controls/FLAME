
#include <iostream>
#include <sstream>
#include <set>

#include <flame/config.h>
#include <flame/util.h>

#include "glps_parser.h"

Config::Config()
    :values(new values_t)
{}

Config::Config(const Config& O)
    :values(O.values)
    ,implicit_values(O.implicit_values)
{}

Config&
Config::operator=(const Config& O)
{
    if(this!=&O) {
        values = O.values;
        implicit_values = O.implicit_values;
    }
    return *this;
}

//! Ensure we have a unique reference to our values std::map by making a copy if necessary
void Config::_cow()
{
    if(!values.unique()) {
        Config::values_pointer U(new values_t(*values)); // copy
        U.swap(values);
    }
}

bool
Config::tryGetAny(const std::string& name, value_t& ret) const
{
    values_t::const_iterator it=values->find(name);
    if(it!=values->end()) {
        ret = it->second;
        return true;
    }
    if(implicit_values) {
        it = implicit_values->find(name);
        if(it!=implicit_values->end()) {
            ret = it->second;
            return true;
        }
    }
    return false;
}

const Config::value_t&
Config::getAny(const std::string& name) const
{
    values_t::const_iterator it=values->find(name);
    if(it!=values->end()) return it->second;
    if(implicit_values) {
        it = implicit_values->find(name);
        if(it!=implicit_values->end()) return it->second;
    }
    throw key_error(SB()<<"Missing parameter '"<<name<<"'");
}
void
Config::setAny(const std::string& name, const value_t& val)
{
    _cow();
    (*values)[name] = val;
}

void
Config::swapAny(const std::string& name, value_t& val)
{
    _cow();
    {
        values_t::iterator it = values->find(name);
        if(it!=values->end()) {
            it->second.swap(val);
            return;
        }
    }
    std::pair<values_t::iterator, bool> ret = values->insert(std::make_pair(name,value_t()));
    assert(ret.second);
    ret.first->second.swap(val);
}

Config Config::new_scope() const
{
    Config ret;
    // return new Config with empty 'values' and 'implicit_values' containing our 'implicit_values'+'values'
    if(values->empty()) {
        ret.implicit_values = implicit_values;
    } else if(!implicit_values || implicit_values->empty()) {
        ret.implicit_values = values; // _cow() makes this safe
    } else {
        values_pointer ptr(new values_t(*values)); // copy
        ptr->insert(implicit_values->begin(),
                    implicit_values->end());
        // note that insert() will not overwrite existing keys
        ret.implicit_values = ptr;
    }
    return ret;
}

void Config::push_scope()
{
    flatten();
    implicit_values = values; // _cow() makes this safe
    values.reset(new values_t);
}

void Config::flatten()
{
    if(implicit_values) {
        values->insert(implicit_values->begin(),
                       implicit_values->end());
        implicit_values.reset();
        // note that insert() will not overwrite existing keys
    }
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
    for(Config::values_t::const_iterator it=values->begin(), end=values->end();
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
    case glps_expr_config: {
        boost::shared_ptr<Config> pconf(boost::get<boost::shared_ptr<Config> >(expr.value));
        std::vector<Config> cvect;
        if(pconf)
            cvect.push_back(*pconf);
        conf.set<std::vector<Config> >(name, cvect);
    }
        break;
    default:
        throw std::logic_error("Context contained unresolved/illegal variable");
    }
}
}

struct GLPSParser::Pvt {
    typedef Config::values_t values_t;
    values_t vars;
    std::ostream *printer;

    Pvt() :printer(&std::cerr) {}

    void fill_vars(parse_context& ctxt)
    {
        for(values_t::const_iterator it=vars.begin(), end=vars.end(); it!=end; ++it)
        {
            // fill in ctxt.vars and ctxt.var_idx
            boost::apply_visitor(store_ctxt_var(ctxt, it->first), it->second);
        }
    }

    Config* fill_context(parse_context& ctxt, const bool lattice=true)
    {
        std::unique_ptr<Config> ret(new Config);
        ret->reserve(ctxt.vars.size()+2);

        // copy ctxt.vars to top level Config
        for(parse_context::vars_t::iterator it=ctxt.vars.begin(), end=ctxt.vars.end();
            it!=end; ++it)
        {
            assign_expr_to_Config(*ret, it->name, it->expr);
        }

        if(lattice){
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
                Config next(ret->new_scope()); // inhiert global scope
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
        }

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

void
GLPSParser::setPrinter(std::ostream* strm)
{
    priv->printer = strm;
}

Config*
GLPSParser::parse_file(const char *fname, const bool have_lattice)
{
    boost::filesystem::path fpath;
    if(fname) {
        fpath = boost::filesystem::canonical(fname).parent_path();
    } else {
        fpath = boost::filesystem::current_path();
    }

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
        Config *ret = parse_file(have_lattice, fp, fpath.native().c_str());
        if(closeme) fclose(fp);
        return ret;
    }catch(...){
        if(closeme) fclose(fp);
        throw;
    }
}

Config*
GLPSParser::parse_file(const bool have_lattice, FILE *fp, const char *path)
{
    parse_context ctxt(path);
    ctxt.printer = priv->printer;
    priv->fill_vars(ctxt);
    ctxt.parse(fp);
    return priv->fill_context(ctxt, have_lattice);
}

Config*
GLPSParser::parse_byte(const char* s, size_t len, const char *path)
{
    parse_context ctxt(path);
    ctxt.printer = priv->printer;
    priv->fill_vars(ctxt);
    ctxt.parse(s, len);
    return priv->fill_context(ctxt);
}

Config*
GLPSParser::parse_byte(const std::string& s, const char *path)
{
    parse_context ctxt(path);
    ctxt.printer = priv->printer;
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
