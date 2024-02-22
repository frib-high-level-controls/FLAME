
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <sstream>

#include "glps_parser.h"

# define M_PI 3.14159265358979323846

#define YYSTYPE         GLPS_STYPE
extern "C" {
// the generated headers wrap some parts in extern "C" blocks, but not all...
#include "glps.par.h"
#include "glps.tab.h"
}

parse_element::parse_element(std::string L, std::string E, kvlist_t::map_t &M)
    :label(L), etype(E)
{
    assert(!label.empty() && !etype.empty());
    props.swap(M);
}

parse_line::parse_line(std::string L, std::string E, strlist_t::list_t& N)
    :label(L), etype(E)
{
    names.swap(N);
}

operation_t::operation_t(const char *name, eval_t fn, glps_expr_type R, unsigned N, va_list args)
    :name(name), fn(fn), result_type(R), arg_types(N)
{
    for(unsigned i=0; i<N; i++)
        arg_types[i] = (glps_expr_type)va_arg(args, int);
}

void glps_error(void *scan, parse_context *ctxt, const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    glps_verror(scan, ctxt, fmt, args);
    va_end(args);
}

void glps_verror(void *scan, parse_context *ctxt, const char*fmt, va_list args)
{
    if(ctxt->last_error.size())
        return; // don't clobber first error
    vsnprintf(&ctxt->error_scratch[0], ctxt->error_scratch.size()-1, fmt, args);
    ctxt->error_scratch[ctxt->error_scratch.size()-1] = '\0';
    ctxt->last_error = &ctxt->error_scratch[0];
    ctxt->last_line = glps_get_lineno(scan);
}

const char *glps_expr_type_name(glps_expr_type e)
{
    switch(e) {
    case glps_expr_number: return "Number";
    case glps_expr_vector: return "Vector";
    case glps_expr_string: return "String";
    case glps_expr_var:    return "Variable";
    case glps_expr_config: return "Config";
    case glps_expr_elem:   return "Element";
    case glps_expr_line:   return "Beamline";
    case glps_expr_invalid:return "Invalid";
    default: return "Unknown_expression_type";
    }
}

void glps_string_debug(FILE *fp, const string_t *s)
{
    fprintf(fp, "%s", s->str.c_str());
}

void glps_expr_debug(FILE *fp, const expr_t *E)
{
    fprintf(fp, "%p type %s", E, glps_expr_type_name(E->etype));
    if(E->etype==glps_expr_line) {
        try{
            const strlist_t::list_t& L(boost::get<strlist_t::list_t>(E->value));
            fprintf(fp, " [%u] (", (unsigned)L.size());
            for(size_t i=0, N=L.size(); i<N; i++)
                fprintf(fp, "%s, ", L[i].c_str());
            fprintf(fp, ")");
        }catch(std::exception& e){
            fprintf(fp, " oops %s", e.what());
        }
    }
}

string_t* glps_string_alloc(const char *s, size_t n)
{
    try{
        return new string_t(s, n);
    }catch(...){
        return NULL;
    }
}

void glps_string_cleanup(string_t* str)
{
    delete str;
}

void glps_expr_cleanup(expr_t* expr)
{
    delete expr;
}

void glps_vector_cleanup(vector_t* vec)
{
    delete vec;
}

void glps_kvlist_cleanup(kvlist_t* pval)
{
    delete pval;
}

void glps_strlist_cleanup(strlist_t* pval)
{
    delete pval;
}

kvlist_t* glps_append_kv(parse_context *ctxt, kvlist_t* L, kv_t* V)
{
    std::unique_ptr<string_t> SN(V->key);
    std::unique_ptr<expr_t> EV(V->value);
    if(!L) {
        try{
            L = new kvlist_t;
        } catch(std::bad_alloc&) {
            glps_error(ctxt->scanner, ctxt, "Allocation failure");
            return NULL;
        }
    }
    L->map[V->key->str] = *V->value;
    return L;
}

strlist_t* glps_append_expr(parse_context *ctxt, strlist_t* list, expr_t *expr)
{
    assert(expr);
    std::unique_ptr<expr_t> E(expr);
    std::unique_ptr<strlist_t> L(list);

    try{
        if(!L.get()) {
            L.reset(new strlist_t);
        }

        switch(expr->etype) {
        case glps_expr_elem:
            // prepend a single element (append, the result will be reversed in glps_add_line)
            L->list.push_back(boost::get<std::string>(expr->value));
            break;
        case glps_expr_line:
        {
            // prepend another line (append in reversed order)
            std::vector<std::string>& N(boost::get<std::vector<std::string> >(expr->value));

            size_t orig = L->list.size();
            L->list.resize(L->list.size()+N.size());
            std::copy(N.rbegin(),
                      N.rend(),
                      L->list.begin()+orig);
        }
            break;
        default:
            glps_error(ctxt->scanner, ctxt, "Lines can not be constructed from %s",
                       glps_expr_type_name(expr->etype));
            L.reset(NULL);
        }

    } catch(boost::bad_get& e) {
        glps_error(ctxt->scanner, ctxt, "Error appending to vector: incorrect type %s", expr->value.type().name());
        L.reset(NULL);
    } catch(std::exception& e) {
        glps_error(ctxt->scanner, ctxt, "Error appending to vector: %s", e.what());
        L.reset(NULL);
    }

    return L.release();
}

vector_t* glps_append_vector(parse_context *ctxt, vector_t *list, expr_t *expr)
{
    assert(expr);
    try{
        std::unique_ptr<expr_t> E(expr);
        std::unique_ptr<vector_t> V(list);
        if(!V.get())
            V.reset(new vector_t);

        if(expr->etype==glps_expr_number) {
            V->value.push_back(boost::get<double>(expr->value));
        } else {
            std::ostringstream strm;
            strm << "Vector element types must be scalar not type " << glps_expr_type_name(expr->etype);
            throw std::runtime_error(strm.str());
        }

        return V.release();
    } catch(std::exception& e) {
        glps_error(ctxt->scanner, ctxt, "Error appending to vector: %s", e.what());
        return NULL;
    }
}

expr_t *glps_add_value(parse_context *ctxt, glps_expr_type t, ...)
{
    std::unique_ptr<expr_t> ret;
    try {
        ret.reset(new expr_t);
        ret->etype = t;

        va_list args;
        va_start(args, t);
        ret->etype = t;
        switch(t) {
        case glps_expr_number:
            ret->value = va_arg(args, double);
            break;
        case glps_expr_var:
        {
            string_t *name = va_arg(args, string_t*);
            std::unique_ptr<string_t> SN(name);

            parse_context::map_idx_t::const_iterator it;

            // var may actually be a variable name, or an element/line label
            if((it=ctxt->var_idx.find(name->str))!=ctxt->var_idx.end()) {
                // expand variable
                expr_t &E = ctxt->vars[it->second].expr;

                ret->etype = E.etype;
                ret->value = E.value;

            } else if(ctxt->element_idx.find(name->str)!=ctxt->element_idx.end()) {
                // expand element
                strlist_t::list_t T(1);
                T[0] = name->str;
                ret->etype = glps_expr_line;
                ret->value = T;

            } else if((it=ctxt->line_idx.find(name->str))!=ctxt->line_idx.end()) {
                // expand beamline
                parse_line &L = ctxt->line[it->second];
                ret->etype = glps_expr_line;
                ret->value = L.names;

            } else if(name->str=="pi") {
                ret->etype = glps_expr_number;
                ret->value = M_PI;

            } else {
                /* having this check in the parser ensure that variables
                 * are defined before first use, which prevents definition
                 * loops from forming.
                 */
                glps_error(ctxt->scanner, ctxt, "Variable/Label '%s' referenced before definition", name->str.c_str());
                ret.reset();
            }
        }
            break;
        case glps_expr_string:
        {
            string_t *name = va_arg(args, string_t*);
            std::unique_ptr<string_t> SN(name);
            assert(name);
            ret->value = name->str;
        }
            break;
        case glps_expr_vector:
        {
            vector_t *vec = va_arg(args, vector_t *);
            if(!vec)
                ret->value = std::vector<double>();
            else {
                std::reverse(vec->value.begin(), vec->value.end());
                ret->value = vec->value;
                delete vec;
            }
            break;
        }
        default:
            glps_error(ctxt->scanner, ctxt, "Logic error, type=%d", (int)t);
            ret.reset();
            break;
        }
        va_end(args);

    }catch(std::exception& e){
        glps_error(ctxt->scanner, ctxt, "Exception while constructing value: %s", e.what());
        return NULL;
    }
    return ret.release();
}

static
std::string glps_describe_op(const operation_t* op)
{
    std::ostringstream strm;
    strm << glps_expr_type_name(op->result_type)
         << " '" << op->name <<"'(";
    for(unsigned i=0; i<op->arg_types.size(); i++) {
        if(i>0)
            strm<<", ";
        strm<<glps_expr_type_name(op->arg_types[i]);
    }
    strm << ")";
    return strm.str();

}

expr_t *glps_add_op(parse_context *ctxt, string_t *name, unsigned N, expr_t **args)
{
    std::unique_ptr<string_t> SN(name);
    std::unique_ptr<expr_t> ret;
    try{

        parse_context::operations_iterator_pair opit = ctxt->operations.equal_range(name->str);

        const operation_t *op=NULL;
        for(;opit.first!=opit.second; ++opit.first)
        {
            if(opit.first->second.arg_types.size()!=N)
                continue;
            bool match=true;
            for(unsigned I=0; I<N; I++)
            {
                if(opit.first->second.arg_types[I]!=args[I]->etype) {
                    match=false;
                    break;
                }
            }
            if(match) {
                op = &opit.first->second;
                break;
            }
        }

        if(!op) {
            std::ostringstream strm;
            strm << "Operation '" << name <<"'(";
            for(unsigned i=0; i<N; i++) {
                if(i>0)
                    strm<<", ";
                strm<<glps_expr_type_name(args[i]->etype);
            }
            strm << ")";
            glps_error(ctxt->scanner, ctxt, "%s", strm.str().c_str());

        } else {
            ret.reset(new expr_t);
            ret->etype = op->result_type;
            int rv = 0;
            try{
                rv = (*op->fn)(ctxt, &ret->value, args);
            }catch(std::exception& e){
                glps_error(ctxt->scanner, ctxt, "User exception evaluating '%s': %s",
                           glps_describe_op(op).c_str(), e.what());
                ret.reset(0);
            }

            if(rv)
                throw std::runtime_error(ctxt->last_error);
        }


    }catch(std::exception& e){
        glps_error(ctxt->scanner, ctxt, "Exception evaluating expression op '%s': %s",
                   name->str.c_str(), e.what());
        ret.reset(0);
    }

    while(N) glps_expr_cleanup(args[--N]);
    return ret.release();
}


void glps_assign(parse_context *ctxt, string_t *name, expr_t*value)
{
    assert(name);
    assert(value);
    try{
        std::unique_ptr<string_t> SN(name);
        std::unique_ptr<expr_t> VN(value);

        bool ok = false;
        parse_var V;
        V.name = name->str;

        switch(value->etype) {
        case glps_expr_vector:
        case glps_expr_number:
        case glps_expr_string:
        case glps_expr_config:
        {
            V.expr = *value;
            ok = true;
        }
            break;
        default:
            glps_error(ctxt->scanner, ctxt, "expression type %d may not be assigned to variable '%s'",
                       value->etype, V.name.c_str());

        }

        if(ok) {
            if(ctxt->var_idx.find(V.name)!=ctxt->var_idx.end()) {
                glps_error(ctxt->scanner, ctxt, "Name '%s' already defined", V.name.c_str());
            } else {
                ctxt->vars.push_back(V);
                ctxt->var_idx[V.name] = ctxt->vars.size()-1;
            }
        }

    } catch(std::exception& e) {
        glps_error(ctxt->scanner, ctxt, "Variable '%s' assignment error; %s", name->str.c_str(), e.what());
    }
}

void glps_add_element(parse_context *ctxt, string_t *label, string_t *etype, kvlist_t *P)
{
    std::unique_ptr<string_t> SL(label), SE(etype);
    std::unique_ptr<kvlist_t> props(P);
    try{
        if(!P)
            props.reset(new kvlist_t);
        if(label->str.empty() || etype->str.empty()) {
            glps_error(ctxt->scanner, ctxt, "Element with null label or type");

        } else if(ctxt->element_idx.find(label->str)!=ctxt->element_idx.end()) {
            glps_error(ctxt->scanner, ctxt, "Name '%s' already used", label->str.c_str());

        } else {
            ctxt->elements.push_back(parse_element(label->str, etype->str, props->map));
            ctxt->element_idx[label->str] = ctxt->elements.size()-1;
        }

    } catch(std::exception& e) {
        glps_error(ctxt->scanner, ctxt, "Element '%s' definition error; %s", label->str.c_str(), e.what());
    }
}

void glps_add_line(parse_context *ctxt, string_t *label, string_t *etype, strlist_t *N)
{
    std::unique_ptr<string_t> SL(label), SE(etype);
    try{
        std::unique_ptr<strlist_t> names(N);
        if(!N)
            names.reset(new strlist_t);

        if(strcasecmp(etype->str.c_str(),"LINE")!=0) {
            glps_error(ctxt->scanner, ctxt, "line-like definition %s with '%s' instead of 'LINE'",
                       label->str.c_str(), etype->str.c_str());

        } else if(ctxt->line_idx.find(label->str)!=ctxt->line_idx.end()) {
            glps_error(ctxt->scanner, ctxt, "Name '%s' already used", label->str.c_str());

        } else {
            // reverse order of elements
            std::reverse(names->list.begin(), names->list.end());
            ctxt->line.push_back(parse_line(label->str, etype->str, names->list));
            ctxt->line_idx[label->str] = ctxt->line.size()-1;
        }

    } catch(std::exception& e) {
        glps_error(ctxt->scanner, ctxt, "Line '%s' definition error; %s", label->str.c_str(), e.what());
    }
}

void glps_command(parse_context* ctxt, string_t *kw)
{
    std::unique_ptr<string_t> SK(kw);
    if(strcmp(kw->str.c_str(), "END")!=0) {
        glps_error(ctxt->scanner, ctxt, "Undefined command '%s'", kw->str.c_str());
    }
}

void glps_call1(parse_context *ctxt, string_t *func, expr_t *arg)
{
    std::unique_ptr<string_t> FN(func);
    std::unique_ptr<expr_t> ARG(arg);

    if(func->str!="print") {
        glps_error(ctxt->scanner, ctxt, "Undefined global function '%s'", func->str.c_str());
    } else if(ctxt->printer) {
        std::ostream& strm = *ctxt->printer;
        strm<<"On line "<<glps_get_lineno(ctxt->scanner)<<" : ";
        switch(arg->etype)
        {
        case glps_expr_number: strm<< boost::get<double>(arg->value); break;
        case glps_expr_string: strm<<"\""<< boost::get<std::string>(arg->value)<<"\""; break;
        case glps_expr_vector: {
            std::ostream_iterator<double> it(strm, ", ");
            const std::vector<double>& vect = boost::get<std::vector<double> >(arg->value);
            strm<<"[";
            std::copy(vect.begin(), vect.end(), it);
            strm<<"]";
            }
            break;
        case glps_expr_line: {
            std::ostream_iterator<std::string> it(strm, ", ");
            const std::vector<std::string>& vect = boost::get<std::vector<std::string> >(arg->value);
            strm<<"(";
            std::copy(vect.begin(), vect.end(), it);
            strm<<")";
            }
            break;
        default:
            strm<<"?? <"<<glps_expr_type_name(arg->etype)<<"> ??";
        }
        strm<<"\n";
    }
}

namespace {
void parse_common(YY_BUFFER_STATE (*setup)(yyscan_t, void*),
                  parse_context* ctxt, void *pvt)
{
    yyscan_t scanner;

    if(glps_lex_init_extra(ctxt, &scanner))
        throw std::runtime_error("Failed to allocate/init lexer");
    try{
        YY_BUFFER_STATE buf = setup(scanner, pvt);
        if(!buf)
            throw std::runtime_error("Failed to allocate lex buffer");

        ctxt->last_error.clear();
        glps_set_lineno(1, scanner);

        yacc_arg arg;
        ctxt->scanner = arg.scanner = scanner;
        arg.ctxt = ctxt;
        glps_parse(arg);
        ctxt->scanner = NULL;

        if(ctxt->last_error.size()) {
            std::ostringstream strm;
            strm << "On line " << ctxt->last_line;
            strm << ": " << ctxt->last_error;
            throw std::runtime_error(strm.str());
        }

        glps_lex_destroy(scanner);
    } catch(...) {
        glps_lex_destroy(scanner);
        throw;
    }
}

YY_BUFFER_STATE setup_file(yyscan_t scanner, void* pvt)
{
    FILE *fp = (FILE*)pvt;
    YY_BUFFER_STATE ret = glps__create_buffer(fp, 1024, scanner);
    glps__switch_to_buffer(ret, scanner);
    return ret;
}

YY_BUFFER_STATE setup_string(yyscan_t scanner, void* pvt)
{
    const std::string *str = (const std::string*)pvt;
    return glps__scan_bytes(str->c_str(), str->size(), scanner);
}

struct bytebuf {
    const char *bytes;
    size_t len;
};

YY_BUFFER_STATE setup_bytes(yyscan_t scanner, void* pvt)
{
    bytebuf *buf = static_cast<bytebuf*>(pvt);
    return glps__scan_bytes(buf->bytes, buf->len, scanner);
}
}

void parse_context::parse(FILE *fp)
{
    parse_common(&setup_file, this, (void*)fp);
}

void parse_context::parse(const char* s, size_t len)
{
    bytebuf buf = {s, len};
    parse_common(&setup_bytes, this, (void*)&buf);
}

void parse_context::parse(const std::string& s)
{
    parse_common(&setup_string, this, (void*)&s);
}

void
parse_context::addop(const char *name,
                     operation_t::eval_t fn,
                     glps_expr_type R, unsigned N, ...)
{
    va_list args;
    va_start(args, N);
    try{
        operations.insert(std::make_pair(name,
                                         operation_t(name, fn, R, N, args)
                                         ));

        va_end(args);
    }catch(...){
        va_end(args);
        throw;
    }
}
