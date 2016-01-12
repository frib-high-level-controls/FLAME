#define BOOST_TEST_MODULE config
#include <boost/test/included/unit_test.hpp>

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
