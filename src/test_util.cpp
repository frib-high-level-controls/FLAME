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
    BOOST_CHECK_THROW(tbl.read(strm), std::runtime_error)
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
    BOOST_CHECK_THROW(tbl.read(strm), std::runtime_error)
}
