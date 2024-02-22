#ifndef FLAME_CONFIG_H
#define FLAME_CONFIG_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
} // extern "C"

#include <ostream>
#include <vector>
#include <map>
#include <memory>
#include <stdexcept>

#include <boost/shared_ptr.hpp>
#include <boost/variant.hpp>
#include <boost/call_traits.hpp>
#include <boost/static_assert.hpp>

#include "util.h"

namespace detail {
// helper to return POD types by value and complex types by const reference
template<typename X>
struct RT { typedef const X& type; };
template<> struct RT<double> { typedef double type; };

template<typename V>
struct buildval {
    static inline V op() { return V(); }
};
template<>
struct buildval<double> {
    static inline double op() { return 0.0; }
};

// helper to ensure that attempts to call Config::get<T> for unsupported T will fail to compile.
template<typename T>
struct is_config_value {
};
#define IS_CONFIG_VALUE(TYPE) \
namespace detail {template<> struct is_config_value<TYPE> { typedef TYPE type; };}
} // namespace detail
IS_CONFIG_VALUE(double)
IS_CONFIG_VALUE(std::string)
IS_CONFIG_VALUE(std::vector<double>)

/** @brief Associative configuration container
 *
 * Typed key/value storage.
 * Value types must be:
 * # double,
 * # std::vector<double>,
 * # std::string, or
 * # std::vector<Config>
 *
 * Most common usage is get<>() and set<>() to fetch and store typed values.
 * Generic code might also use getAny() or setAny().
 *
 * Also has the notion
 */
class Config
{
public:
    //! An individual value (double, double[], string, or Config[])
    typedef boost::variant<
        double,
        std::vector<double>,
        std::string,
        std::vector<Config>
    > value_t;

    typedef std::vector<Config> vector_t;

    typedef std::map<std::string, value_t> values_t;
private:
    typedef boost::shared_ptr<values_t> values_pointer;
    typedef boost::shared_ptr<const values_t> const_values_pointer;

    values_pointer values;
    const_values_pointer implicit_values;

    void _cow();
public:
    //! New empty config
    Config();
    //! Build config from a single std::map (single scope)
    //! Copy ctor
    Config(const Config&);
    ~Config() {}
    //! Assignment
    Config& operator=(const Config&);

    /** lookup untyped.
     * @returns true if 'ret' updates, and false if no parameter with 'name'.
     */
    bool tryGetAny(const std::string& name, value_t& ret) const;
    /** lookup untyped.
     * @throws key_error if name doesn't refer to an existing parameter
     */
    const value_t& getAny(const std::string& name) const;
    /** add/replace with a new value, untyped
     */
    void setAny(const std::string& name, const value_t& val);
    //! Exchange a single parameter untyped
    void swapAny(const std::string& name, value_t& val);

    /** lookup typed.
     * @throws key_error if name doesn't refer to an existing parameter, or the parameter has other
     *                   than the requested type.
     *
     @code
     Config X;
     double v = X.get<double>("param"); // throws since "param" is not set
     @endcode
     */
    template<typename T>
    typename detail::RT<T>::type
    get(const std::string& name) const {
        try {
            return boost::get<typename detail::is_config_value<T>::type>(getAny(name));
        } catch(boost::bad_get&) {
            throw key_error(SB()<<"Wrong type for '"<<name<<"'.  should be "<<typeid(T).name());
        }
    }
    /** lookup typed with default.
     * If 'name' doesn't refer to a parameter, or it has the wrong type,
     * then 'def' is returned instead.
     *
     @code
     Config X;
     double v = X.get<double>("param", 0.0); // returns 0.0 since "param" is not set
     @endcode
     */
    template<typename T>
    typename detail::RT<T>::type
    get(const std::string& name, typename boost::call_traits<T>::param_type def) const {
        try{
            return boost::get<typename detail::is_config_value<T>::type>(getAny(name));
        } catch(boost::bad_get&) {
        } catch(key_error&) {
        }
        return def;
    }

    /** lookup where missing parameters returns false
     *
     * @returns true if val is updated, false otherwise
     */
    template<typename T>
    bool
    tryGet(const std::string& name, T& val) const {
        value_t ret;
        if(tryGetAny(name, ret)) {
            try{
                val = boost::get<typename detail::is_config_value<T>::type>(ret);
                return true;
            } catch(boost::bad_get&) {
            }
        }
        return false;
    }

    /** add/replace with a new value
     @code
     Config X;
     X.set<double>("param", 42.0);
     @endcode
     */
    template<typename T>
    void set(const std::string& name,
             typename boost::call_traits<typename detail::is_config_value<T>::type>::param_type val)
    {
        _cow();
        (*values)[name] = val;
    }

    /** Exchange a single parameter typed
     @code
     Config X;
     X.set<double>("param", 42.0);
     double v = 43.0;
     X.swap<double>("param", v);
     assert(v==42.0);
     @endcode
     */
    template<typename T>
    void swap(const std::string& name,
              typename boost::call_traits<typename detail::is_config_value<T>::type>::reference val)
    {
        value_t temp = detail::buildval<T>::op();
        val.swap(boost::get<T>(temp));
        swapAny(name, temp);
        //val.swap(boost::get<T>(temp)); // TODO: detect insert, make this a real swap
    }

    //! Exchange entire Config
    void swap(Config& c)
    {
        // no _cow() here since no value_t are modified
        values.swap(c.values);
        implicit_values.swap(c.implicit_values);
    }

    //! Print listing of inner scope
    void show(std::ostream&, unsigned indent=0) const;

    //! iterator
    typedef values_t::iterator iterator;
    //! const_iterator
    typedef values_t::const_iterator const_iterator;

    //! The first element
    inline const_iterator begin() const { return values->begin(); }
    //! one after the last element
    inline const_iterator end() const { return values->end(); }

    inline void reserve(size_t) {}
    //! Create a new inner scope
    Config new_scope() const;
    void push_scope();

    void flatten();
};

IS_CONFIG_VALUE(std::vector<Config>);

inline
std::ostream& operator<<(std::ostream& strm, const Config& c)
{
    c.show(strm);
    return strm;
}

//! @brief Interface to lattice file parser
class GLPSParser
{
    class Pvt;
    std::unique_ptr<Pvt> priv;
public:
    //! Construct an empty parser context
    GLPSParser();
    ~GLPSParser();

    /** @brief Pre-define variable
     *
     * Equivalent to inserting a variable definition at the beginning
     * of a lattice file.
     */
    void setVar(const std::string& name, const Config::value_t& v);
    //! @brief Set output for lexer/parser error messages
    void setPrinter(std::ostream*);

    /** @brief Open and parse a file
     *
     * A wrapper around parse_file(FILE*, const char*).
     *
     * @arg fname File name to open.  If NULL or '-' then parse stdin     *
     *            If stdin is parsed then it is not closed.
     */
    Config *parse_file(const char *fname, const bool have_lattice=true);
    /** @brief Parse from open FILE
     *
     * @arg fp an open file descriptor
     * @args path The directory containing the open file (used to expand relative paths in the file)
     * @post fp is not closed.
     * @returns New Config, which the caller must delete.
     * @throws std::runtime_error For various error conditions
     */
    Config *parse_file(const bool have_lattice, FILE *fp, const char *path=NULL);
    /** @brief Parse from byte buffer
     *
     * @arg s Byte array
     * @arg len Length of byte array (in bytes)
     * @args path A directory to use to expand relative paths when parsing
     * @returns New Config or NULL
     */

    Config *parse_byte(const char* s, size_t len, const char *path=NULL);
    /** @brief Parse from std::string
     *
     * @arg s String
     * @args path A directory to use to expand relative paths when parsing
     * @returns New Config or NULL
     */
    Config *parse_byte(const std::string& s, const char *path=NULL);
};

//! Print a previously parsed, or constructed, Config
void GLPSPrint(std::ostream& strm, const Config&);


#undef IS_CONFIG_VALUE

#endif /* extern "C" */

#endif // FLAME_CONFIG_H
