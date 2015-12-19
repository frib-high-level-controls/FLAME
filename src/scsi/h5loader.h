#ifndef H5LOADER_H
#define H5LOADER_H

#include <ostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>

class H5Loader
{
    struct Pvt;
    Pvt *pvt;
public:
    H5Loader();
    H5Loader(const char *);
    H5Loader(const std::string&);
    ~H5Loader();

    void open(const char *);
    void open(const std::string&);
    void close();

    typedef boost::numeric::ublas::matrix<double,
                    boost::numeric::ublas::row_major,
                    boost::numeric::ublas::unbounded_array<double>
    > matrix_t;

    matrix_t load(const char *);
    matrix_t load(const std::string&);

    static void dontPrint();
};

#endif // H5LOADER_H
