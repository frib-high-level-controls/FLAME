#define BOOST_TEST_MODULE config
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <math.h>

#include "scsi/config.h"

BOOST_AUTO_TEST_CASE(config_getset)
{
    Config C;

    C.set<double>("hello", 4.2);

    BOOST_CHECK_CLOSE(C.get<double>("hello"), 4.2, 0.1);

    BOOST_CHECK_THROW(C.get<double>("world"), key_error);

    BOOST_CHECK_CLOSE(C.get<double>("world", 5.2), 5.2, 0.1);
}

BOOST_AUTO_TEST_CASE(config_scope)
{
    Config C;

    C.set<double>("hello", 4.2);
    C.set<double>("world", 5.2);

    // create a new inner scope
    BOOST_CHECK_EQUAL(C.depth(), 1);
    C.push_scope();
    BOOST_CHECK_EQUAL(C.depth(), 2);

    BOOST_CHECK_CLOSE(C.get<double>("hello"), 4.2, 0.1);
    BOOST_CHECK_CLOSE(C.get<double>("world"), 5.2, 0.1);

    // populate
    C.set<double>("world", 1.2);
    C.set<double>("other", 15.0);

    BOOST_CHECK_CLOSE(C.get<double>("hello"), 4.2, 0.1);  // from outer
    BOOST_CHECK_CLOSE(C.get<double>("world"), 1.2, 0.1);  // from inner
    BOOST_CHECK_CLOSE(C.get<double>("other"), 15.0, 0.1); // from inner

    C.pop_scope();

    BOOST_CHECK_CLOSE(C.get<double>("hello"), 4.2, 0.1); // same
    BOOST_CHECK_CLOSE(C.get<double>("world"), 5.2, 0.1); // from outer
    BOOST_CHECK_THROW(C.get<double>("other"), key_error); // no longer present
}

BOOST_AUTO_TEST_CASE(config_scope_mutate)
{
    Config C;

    C.set<double>("hello", 4.2);
    C.set<double>("world", 5.2);

    Config D(C);

    // create a new inner scope
    C.push_scope();
    BOOST_CHECK_EQUAL(C.depth(), 2);
    BOOST_CHECK_EQUAL(D.depth(), 1);

    // populate
    C.set<double>("world", 1.2);
    C.set<double>("other", 15.0);

    BOOST_CHECK_CLOSE(C.get<double>("hello"), 4.2, 0.1);  // from outer
    BOOST_CHECK_CLOSE(C.get<double>("world"), 1.2, 0.1);  // from inner
    BOOST_CHECK_CLOSE(C.get<double>("other"), 15.0, 0.1); // from inner

    BOOST_CHECK_CLOSE(D.get<double>("hello"), 4.2, 0.1); // same
    BOOST_CHECK_CLOSE(D.get<double>("world"), 5.2, 0.1); // from outer
    BOOST_CHECK_THROW(D.get<double>("other"), key_error); // only inner

    // Change D, which by copy on write is detached from C.
    D.set<double>("hello", 101.0);
    D.set<double>("world", 102.0);
    D.set<double>("other", 103.0);

    // C's outer scope has not changed
    BOOST_CHECK_CLOSE(C.get<double>("hello"), 4.2, 0.1);  // from outer
    BOOST_CHECK_CLOSE(C.get<double>("world"), 1.2, 0.1);  // from inner
    BOOST_CHECK_CLOSE(C.get<double>("other"), 15.0, 0.1); // from inner

    // D has new values
    BOOST_CHECK_CLOSE(D.get<double>("hello"), 101.0, 0.1);
    BOOST_CHECK_CLOSE(D.get<double>("world"), 102.0, 0.1);
    BOOST_CHECK_CLOSE(D.get<double>("other"), 103.0, 0.1);
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

    std::auto_ptr<Config> conf(parse.parse_byte(config_print_stmt_input, sizeof(config_print_stmt_input)-1));

    BOOST_CHECK_EQUAL(strm.str(), "On line 2 : 14\n" "On line 4 : \"test\"\n");
}