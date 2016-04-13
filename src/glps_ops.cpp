
#include "math.h"

#include <sstream>
#include <limits>
#include <stdexcept>

#include "glps_parser.h"

namespace {
// Numeric operations

int unary_negate(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    *R = -boost::get<double>(A[0]->value);
    return 0;
}

int unary_sin(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    *R = sin(boost::get<double>(A[0]->value));
    return 0;
}
int unary_cos(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    *R = cos(boost::get<double>(A[0]->value));
    return 0;
}
int unary_tan(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    *R = tan(boost::get<double>(A[0]->value));
    return 0;
}
int unary_asin(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    *R = asin(boost::get<double>(A[0]->value));
    return 0;
}
int unary_acos(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    *R = acos(boost::get<double>(A[0]->value));
    return 0;
}
int unary_atan(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    *R = atan(boost::get<double>(A[0]->value));
    return 0;
}

int binary_add(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    *R = boost::get<double>(A[0]->value)+boost::get<double>(A[1]->value);
    return 0;
}
int binary_sub(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    *R = boost::get<double>(A[0]->value)-boost::get<double>(A[1]->value);
    return 0;
}
int binary_mult(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    *R = boost::get<double>(A[0]->value)*boost::get<double>(A[1]->value);
    return 0;
}
int binary_div(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    double result = boost::get<double>(A[0]->value)/boost::get<double>(A[1]->value);
    if(!isfinite(result)) {
        ctxt->last_error = "division results in non-finite value";
        return 1;
    }
    *R = result;
    return 0;
}

// beamline operations

int unary_bl_negate(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    // reverse the order of the beamline
    const std::vector<std::string>& line = boost::get<std::vector<std::string> >(A[0]->value);
    std::vector<std::string> ret(line.size());
    std::copy(line.rbegin(),
              line.rend(),
              ret.begin());
    *R = ret;
    return 0;
}

template<int MULT, int LINE>
int binary_bl_mult(parse_context* ctxt, expr_value_t *R, const expr_t * const *A)
{
    // multiple of scale * beamline repeats the beamline 'scale' times
    assert(A[MULT]->etype==glps_expr_number);
    assert(A[LINE]->etype==glps_expr_line);

    double factor = boost::get<double>(A[MULT]->value);

    if(factor<0.0 || factor>std::numeric_limits<unsigned>::max()) {
        ctxt->last_error = "beamline scale by negative value or out of range value";
        return 1;
    }
    unsigned factori = (unsigned)factor;

    const std::vector<std::string>& line = boost::get<std::vector<std::string> >(A[LINE]->value);

    std::vector<std::string> ret(line.size()*factori);

    if(factori>0)
    {
        std::vector<std::string>::iterator outi = ret.begin();

        while(factori--)
            outi = std::copy(line.begin(),
                             line.end(),
                             outi);
    }
    *R = ret;
    return 0;
}

}

parse_context::parse_context(const char *path)
    :last_line(0), error_scratch(300), scanner(NULL)
{
    if(path)
        cwd = boost::filesystem::canonical(path);
    else
        cwd = boost::filesystem::current_path();
    addop("-", &unary_negate, glps_expr_number, 1, glps_expr_number);

    addop("sin", &unary_sin, glps_expr_number, 1, glps_expr_number);
    addop("cos", &unary_cos, glps_expr_number, 1, glps_expr_number);
    addop("tan", &unary_tan, glps_expr_number, 1, glps_expr_number);
    addop("asin",&unary_asin,glps_expr_number, 1, glps_expr_number);
    addop("acos",&unary_acos,glps_expr_number, 1, glps_expr_number);
    addop("atan",&unary_atan,glps_expr_number, 1, glps_expr_number);
    // aliases to capture legacy behavour :P
    addop("arcsin",&unary_asin,glps_expr_number, 1, glps_expr_number);
    addop("arccos",&unary_acos,glps_expr_number, 1, glps_expr_number);
    addop("arctan",&unary_atan,glps_expr_number, 1, glps_expr_number);

    addop("+", &binary_add, glps_expr_number, 2, glps_expr_number, glps_expr_number);
    addop("-", &binary_sub, glps_expr_number, 2, glps_expr_number, glps_expr_number);
    addop("*", &binary_mult,glps_expr_number, 2, glps_expr_number, glps_expr_number);
    addop("/", &binary_div, glps_expr_number, 2, glps_expr_number, glps_expr_number);

    addop("-", &unary_bl_negate, glps_expr_line, 1, glps_expr_line);

    addop("*", &binary_bl_mult<0,1>, glps_expr_line, 2, glps_expr_number, glps_expr_line);
    addop("*", &binary_bl_mult<1,0>, glps_expr_line, 2, glps_expr_line, glps_expr_number);
}

parse_context::~parse_context()
{
}

