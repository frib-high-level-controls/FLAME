#ifndef UTIL_H
#define UTIL_H

#include <map>
#include <ostream>
#include <stdexcept>

#include <boost/call_traits.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#ifdef __GNUC__
#  define UNLIKELY(E) __builtin_expect(E, 0)
#else
#  define UNLIKELY(E) (E)
#endif

struct key_error : public std::runtime_error
{
    key_error(const std::string& s) : std::runtime_error(s) {}
};

//! Helper for std::map lookup
template<typename M>
bool find(const M& map,
          typename boost::call_traits<typename M::key_type>::param_type key,
          typename M::mapped_type& result)
{
    typename M::const_iterator it = map.find(key);
    if(it==map.end()) return false;
    result = it->second;
    return true;
}

//! An iterator which adapts the second member of a std::pair
//! Allows iteration over a values of a std::map<K,V> to look like
//! iteration over a std::list<V> or std::vector<V>.
template<typename I>
class value_proxy_iterator
{
    I orig;
public:
    typedef typename I::iterator_category iterator_category;
    typedef typename I::value_type::second_type value_type;
    typedef std::ptrdiff_t                   difference_type;
    typedef value_type*                  pointer;
    typedef value_type&                  reference;

    value_proxy_iterator() {}
    explicit value_proxy_iterator(const I& o) :orig(o) {}
    value_proxy_iterator(const value_proxy_iterator& o) :orig(o.orig) {}

    reference operator*() const { return orig->second; }
    pointer operator->() const { return &orig->second; }
    reference operator[](size_t i) const { return orig[i]->second; }

    bool operator==(const value_proxy_iterator& o) const { return orig==o.orig; }
    bool operator!=(const value_proxy_iterator& o) const { return orig!=o.orig; }

    value_proxy_iterator& operator++() { ++orig; return *this; }
    value_proxy_iterator  operator++(int) { return value_proxy_iterator(orig++); }
    value_proxy_iterator& operator--() { --orig; return *this; }
    value_proxy_iterator  operator--(int) { return value_proxy_iterator(orig--); }

    value_proxy_iterator& operator+=(size_t i) { orig+=i; return *this; }
    value_proxy_iterator& operator-=(size_t i) { orig-=i; return *this; }
    value_proxy_iterator  operator+(size_t i) const { return value_proxy_iterator(orig+i); }
    value_proxy_iterator  operator-(size_t i) const { return value_proxy_iterator(orig-i); }
};

/** Parse a file containing a single table of numeric values w/ optional column headers
 */
struct numeric_table {
    typedef boost::numeric::ublas::matrix<double,
        boost::numeric::ublas::row_major
    > value_t;

    typedef std::map<std::string, size_t> colnames_t;
    colnames_t colnames;

    value_t table;

    void read(std::istream&);
    void readvec(std::vector<double> vec, int numrow);
};

class numeric_table_cache {
    struct Pvt;
    std::unique_ptr<Pvt> pvt;
public:
    numeric_table_cache();
    ~numeric_table_cache();

    typedef boost::shared_ptr<const numeric_table> table_pointer;

    table_pointer fetch(const std::string& path);

    void clear();

    static numeric_table_cache* get();
};

struct SB {
    std::ostringstream strm;
    SB() {}
    operator std::string() const { return strm.str(); }
    template<typename T>
    SB& operator<<(T i) { strm<<i; return *this; }
};

//! Helper to step through the indicies of an Nd array
template<unsigned MAX>
struct ndindex_iterate {
    bool done;
    const unsigned ndim;
    size_t index[MAX], limit[MAX];
    ndindex_iterate(unsigned nd, size_t *lim) :done(false), ndim(nd) {
        std::fill(index, index+MAX, 0u);
        std::copy(lim, lim+nd, limit);
    }
    bool next() {
        if(!done) {
            for(unsigned nd=0; nd<ndim; nd++) {
                index[nd]++;
                if(index[nd]<limit[nd]) {
                    return done;
                } else {
                    index[nd] = 0;
                }
            }
            done = true;
        }
        return done;
    }
};

#ifdef __GNUC__
#define FLAME_UNUSED __attribute__((unused))
#else
#define FLAME_UNUSED
#endif

#if __cplusplus<201103L && !defined(static_assert)
#define STATIC_JOIN(x, y) STATIC_JOIN2(x, y)
#define STATIC_JOIN2(x, y) x ## y
#define static_assert(expr, msg) \
    typedef int STATIC_JOIN(static_assert_failed_at_line_, __LINE__) \
    [ (expr) ? 1 : -1 ] FLAME_UNUSED
#endif

#endif // UTIL_H
