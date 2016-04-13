#ifndef H5LOADER_H
#define H5LOADER_H

#include <ostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>

/** @brief Helper to read 2d matricies
 *
 @code
 H5Loader L("file.h5:/group");
 H5Loader::matrix M(L.load("dataset"));
 @endcode
 */
class H5Loader
{
    struct Pvt;
    Pvt *pvt;
public:
    //! Construct w/o opening.  Must call open()
    H5Loader();
    //! Construct and open().
    H5Loader(const char *);
    //! Construct and open().
    H5Loader(const std::string&);
    ~H5Loader();

    void open(const char *);
    void open(const std::string&);
    void close();

    //! A 2d matrix
    typedef boost::numeric::ublas::matrix<double,
                    boost::numeric::ublas::row_major,
                    boost::numeric::ublas::unbounded_array<double>
    > matrix_t;

    matrix_t load(const char *);
    matrix_t load(const std::string&);

    static void dontPrint();
};

#endif // H5LOADER_H
