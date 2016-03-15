#ifndef GLPS_H
#define GLPS_H

#include <stdarg.h>

#ifdef __cplusplus
#include <string>
#include <vector>
#include <map>
#include <set>

#include <boost/variant.hpp>
#include <boost/shared_ptr.hpp>
extern "C" {
#endif

typedef enum {
    glps_expr_number, // scalar value
    glps_expr_vector, // vector value
    glps_expr_string, // string value
    glps_expr_var,    // variable name
    glps_expr_line,   // beamline name
    glps_expr_elem,   // element label
    glps_expr_invalid = -1
} glps_expr_type;
const char *glps_expr_type_name(glps_expr_type);

typedef struct parse_context parse_context;

typedef struct {
    void *scanner;
    parse_context *ctxt;
} yacc_arg;

typedef struct string_t string_t;

typedef struct expr_t expr_t;

typedef struct vector_t vector_t;
typedef struct kvlist_t kvlist_t;
typedef struct strlist_t strlist_t;

typedef struct {
    string_t* key;
    expr_t *value;
} kv_t;

string_t* glps_string_alloc(const char *, size_t);

void glps_string_cleanup(string_t*);
void glps_expr_cleanup(expr_t*);
void glps_vector_cleanup(vector_t*);
void glps_kvlist_cleanup(kvlist_t*);
void glps_strlist_cleanup(strlist_t*);

void glps_expr_debug(FILE *, const expr_t *);
void glps_string_debug(FILE *, const string_t *);

void glps_error(void *scan, parse_context *ctxt, const char*, ...) __attribute__((format(printf,3,4)));
void glps_verror(void *scan, parse_context *ctxt, const char*, va_list);

void glps_assign(parse_context*, string_t *, expr_t*);
void glps_add_element(parse_context*, string_t *label, string_t *etype, kvlist_t*);
void glps_add_line(parse_context*, string_t *label, string_t *etype, strlist_t*);
void glps_command(parse_context*, string_t *kw);
void glps_call1(parse_context*, string_t*, expr_t*);

kvlist_t* glps_append_kv(parse_context *ctxt, kvlist_t*, kv_t*);
strlist_t* glps_append_expr(parse_context *ctxt, strlist_t*, expr_t *);
vector_t* glps_append_vector(parse_context *ctxt, vector_t*, expr_t *);

expr_t *glps_add_value(parse_context *ctxt, glps_expr_type t, ...);
expr_t *glps_add_op(parse_context *ctxt, string_t *, unsigned N, expr_t **);

#ifdef __cplusplus
}

typedef boost::variant<
    double, // glps_expr_number
    std::vector<double>, // glps_expr_vector
    std::string, // glps_expr_string,
    std::vector<std::string> // glps_expr_line
> expr_value_t;

struct string_t {
    std::string str;
    template<typename A>
    string_t(A a) : str(a) {}
    template<typename A, typename B>
    string_t(A a, B b) : str(a,b) {}
};

struct expr_t {
    glps_expr_type etype;
    expr_value_t value;

    expr_t(): etype(glps_expr_invalid), value(0.0) {}
    expr_t(glps_expr_type e,  expr_value_t v): etype(e), value(v) {}
};

struct vector_t {
    std::vector<double> value;
};

struct kvlist_t {
    typedef std::map<std::string, expr_t> map_t;
    map_t map;
};

struct strlist_t {
    typedef std::vector<std::string> list_t;
    list_t list;
};

struct parse_var {
    std::string name;
    expr_t expr;

    parse_var() {}
    parse_var(const char* n, expr_t E) :name(n), expr(E) {}
};

struct parse_element {
    std::string label, etype;
    kvlist_t::map_t props;

    parse_element(std::string L, std::string E, kvlist_t::map_t& M);
};

struct parse_line {
    std::string label, etype;
    strlist_t::list_t names;

    parse_line(std::string L, std::string E, strlist_t::list_t& N);
};

struct operation_t {
    typedef int (*eval_t)(parse_context*, expr_value_t *, const expr_t* const *);

    std::string name;
    eval_t fn;
    glps_expr_type result_type;
    std::vector<glps_expr_type> arg_types;

    operation_t(const char *name, eval_t, glps_expr_type R, unsigned N, va_list args);
};

struct parse_context {
    typedef std::vector<parse_var> vars_t;
    vars_t vars;

    typedef std::vector<parse_element> elements_t;
    elements_t elements;

    typedef std::vector<parse_line> line_t;
    line_t line;

    // to enforce uniquness of variable and label names
    typedef std::map<std::string, size_t> map_idx_t;
    map_idx_t var_idx, element_idx, line_idx;

    typedef std::multimap<std::string, operation_t> operations_t;
    typedef std::pair<operations_t::const_iterator, operations_t::const_iterator> operations_iterator_pair;
    operations_t operations;

    void addop(const char *name, operation_t::eval_t, glps_expr_type R, unsigned N, ...);

    std::string last_error;
    unsigned last_line;
    std::ostream *printer;

    std::vector<char> error_scratch;

    parse_context();
    ~parse_context();

    void parse(FILE *fp);
    void parse(const char* s, size_t len);
    void parse(const std::string& s);

    void *scanner;
};

#endif /* __cplusplus */

#endif /* GLPS_H */
