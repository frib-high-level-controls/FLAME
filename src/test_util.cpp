#define BOOST_TEST_MODULE config
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/numeric/ublas/io.hpp>

#include "flame/util.h"

BOOST_AUTO_TEST_CASE(parse_table_empty)
{
    std::istringstream strm;

    numeric_table tbl;
    tbl.read(strm);

    BOOST_CHECK_EQUAL(tbl.table.size1(), 0);
    BOOST_CHECK_EQUAL(tbl.table.size2(), 0);
    BOOST_CHECK_EQUAL(tbl.colnames.size(), 0);
}

static const char table1[] = {
    "%A B\n"
    " 12 13\n"
    "14   15 \n"
    " 16 17\n"
};

BOOST_AUTO_TEST_CASE(parse_table1)
{
    std::istringstream strm(table1);

    numeric_table tbl;
    tbl.read(strm);

    BOOST_CHECK_EQUAL(tbl.table.size1(), 3);
    BOOST_CHECK_EQUAL(tbl.table.size2(), 2);
    BOOST_CHECK_EQUAL(tbl.colnames.size(), 2);

    BOOST_CHECK_EQUAL(tbl.colnames["A"], 0);
    BOOST_CHECK_EQUAL(tbl.colnames["B"], 1);

    std::ostringstream out;
    out<<tbl.table;

    BOOST_CHECK_EQUAL(out.str(), "[3,2]((12,13),(14,15),(16,17))");
}

static const char table2[] = {
    // no column names
    " 12 13\n"
    "14   15\n"
    " 16 17"  // no trailing newline
};

BOOST_AUTO_TEST_CASE(parse_table2)
{
    std::istringstream strm(table2);

    numeric_table tbl;
    tbl.read(strm);

    BOOST_CHECK_EQUAL(tbl.table.size1(), 3);
    BOOST_CHECK_EQUAL(tbl.table.size2(), 2);
    BOOST_CHECK_EQUAL(tbl.colnames.size(), 0);

    std::ostringstream out;
    out<<tbl.table;

    BOOST_CHECK_EQUAL(out.str(), "[3,2]((12,13),(14,15),(16,17))");
}

static const char table3[] = {
    "%A B\n"
    " 12 13\n"
    "14   15 16 17\n"  // wrong number of columns
};

BOOST_AUTO_TEST_CASE(parse_table3)
{
    std::istringstream strm(table3);

    numeric_table tbl;
    BOOST_CHECK_THROW(tbl.read(strm), std::runtime_error);
}

static const char table4[] = {
    "%A B\n"
    " 12 13\n"
    "14   hello\n" // not numeric
    "16 17\n"
};

BOOST_AUTO_TEST_CASE(parse_table4)
{
    std::istringstream strm(table4);

    numeric_table tbl;
    BOOST_CHECK_THROW(tbl.read(strm), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(test_ndindex_iterate)
{
    size_t limits[2] = {2,3};
    ndindex_iterate<2> iter(2, limits);

    BOOST_CHECK(!iter.done);
    BOOST_CHECK_EQUAL(iter.ndim, 2);
    BOOST_CHECK_EQUAL_COLLECTIONS(limits, limits+2, iter.limit, iter.limit+2);

#define CHECKME(isdone, I, J) {size_t P[2] = {I,J}; BOOST_CHECK_EQUAL(isdone, iter.done); \
    BOOST_CHECK_EQUAL_COLLECTIONS(P, P+2, iter.index, iter.index+2); }

    CHECKME(false, 0, 0);
    iter.next();
    CHECKME(false, 1, 0);
    iter.next();
    CHECKME(false, 0, 1);
    iter.next();
    CHECKME(false, 1, 1);
    iter.next();
    CHECKME(false, 0, 2);
    iter.next();
    CHECKME(false, 1, 2);
    iter.next();
    CHECKME(true,  0, 0);

#undef CHECKME
}
