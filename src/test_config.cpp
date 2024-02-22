#define BOOST_TEST_MODULE config
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <math.h>

#include "flame/config.h"

BOOST_AUTO_TEST_CASE(config_getset)
{
    Config C;

    C.set<double>("hello", 4.2);

    BOOST_CHECK_CLOSE(C.get<double>("hello"), 4.2, 0.1);

    BOOST_CHECK_THROW(C.get<double>("world"), key_error);

    BOOST_CHECK_CLOSE(C.get<double>("world", 5.2), 5.2, 0.1);
}

static const char config_print_stmt_input[] =
"X = 14;\n"
"print(X);\n"
"Y = \"test\";\n"
"print(Y);\n"
"foo : bar, x=\"y\";\n"
"baz : LINE = (foo);\n";

BOOST_AUTO_TEST_CASE(config_print_stmt)
{
    std::ostringstream strm;
    GLPSParser parse;
    parse.setPrinter(&strm);

    std::unique_ptr<Config> conf(parse.parse_byte(config_print_stmt_input, sizeof(config_print_stmt_input)-1));

    BOOST_CHECK_EQUAL(strm.str(), "On line 2 : 14\n" "On line 4 : \"test\"\n");
}
