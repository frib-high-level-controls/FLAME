#ifndef FLAME_BASE_H
#define FLAME_BASE_H

#include <stdlib.h>

#include <ostream>
#include <string>
#include <map>
#include <vector>
#include <utility>

#include <boost/noncopyable.hpp>
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/call_traits.hpp>

#include "config.h"
#include "util.h"

#define FLAME_API_VERSION 0

struct ElementVoid;

/** @brief The abstract base class for all simulation state objects.
 *
 * Represents the state of some bunch of particles moving through a Machine
 */
struct StateBase : public boost::noncopyable
{
    virtual ~StateBase();

    //! Index of ElementVoid in Machine to follow this the current one.
    //! May be altered within ElementVoid::advance() to achieve branching or looping.
    size_t next_elem;

    double pos;        //!< absolute longitudinal position at end of Element

    //! virtual equivalent to operator=()
    //! Should only be used with another State originating from the same Machine
    //! from Machine::allocState() or clone().
    virtual void assign(const StateBase& other) =0;

    //! Print information about the state.
    //! level is a hint as to the verbosity expected by the caller.
    virtual void show(std::ostream&, int level =0) const {}

    //! Used with StateBase::getArray() to describe a single parameter
    struct ArrayInfo {
        enum {maxdims=3};
        ArrayInfo() :name(0), type(Double), ptr(NULL), ndim(0) {}
        //! The parameter name
        const char *name;
        //! The parameter type Double (double) or Sizet (size_t)
        enum Type {
            Double, Sizet
        } type;
        //! Pointer to parameter storage.
        //! Actual type depends on the type: Double (double) or Sizet (size_t)
        void *ptr;
        //! Number of dimensions.
        //! Indicates how many entries in dim[] and stride[] are valid.
        //! ndim==0 indicates a scalar.
        //! @pre ndim<=maxdims
        unsigned ndim;
        //! Array dimensions in elements.
        size_t dim[maxdims];
        //! Array strides in bytes
        size_t stride[maxdims];

        //! is the given index valid?
        bool inbounds(size_t* d) const {
            bool ret = true;
            switch(ndim) {
            case 3: ret &= d[2]<dim[2];
            case 2: ret &= d[1]<dim[1];
            case 1: ret &= d[0]<dim[0];
            }
            return ret;
        }

        void *raw(size_t* d) {
            char *ret = (char*)ptr;
            switch(ndim) {
            case 3: ret += d[2]*stride[2];
            case 2: ret += d[1]*stride[1];
            case 1: ret += d[0]*stride[0];
            }
            return ret;
        }

        //! Helper to fetch the pointer for a given index (assumed valid)
        template<typename E>
        inline E* get(size_t *d) {
            return (E*)raw(d);
        }
    };

    /** @brief Introspect named parameter of the derived class
     * @param idx The index of the parameter
     * @param Info mailbox to be filled with information about this paramter
     * @return true if 'idx' is valid and Info was filled in, otherwise false
     *
     * To introspect call this with index increasing from zero until false is returned.
     *
     * Sub-classes are required to maintain a fixed number of parameters.
     * The same ArrayInfo::name should alwaysbe returned for a given index.
     * Similarly ArrayInfo::ndim and ArrayInfo.type should not change.
     *
     * ArrayInfo::ptr, ArrayInfo::dim, ArrayInfo::stride
     * are allowed to change when the state is passed to Machine::propagate().
     *
     * Therefore, if getArray() for an instance has returned true for some index,
     * then a caller may assume that all future calls to getArray() for this instance
     * will also return true.
     * However, the caller must still re-call getArray() in future as the storage
     * and sizes may have changed.
     */
    virtual bool getArray(unsigned index, ArrayInfo& Info);

    //! Allocate a new instance which is a copy of this one.
    //! Caller is responsible to delete the returned pointer
    virtual StateBase* clone() const =0;

    //! @private
    //! Mailbox to hold the python interpreter object wrapping us.
    void *pyptr;
protected:
    /** Construct a new state
     *
     * By convention an empty Config() should be valid.
     */
    StateBase(const Config& c);
    //! Used with clone ctor
    struct clone_tag{};
    //! For use in clone()
    StateBase(const StateBase& c, clone_tag);
};

inline
std::ostream& operator<<(std::ostream& strm, const StateBase& s)
{
    s.show(strm, 0);
    return strm;
}

/**
 * @brief Allow inspection of intermediate State
 *
 * Use with ElementVoid::set_observer() to associate an observer with an Element.
 * During Machine::propagate() a call will be made to Observer::view()
 * with the element's output State.
 */
struct Observer : public boost::noncopyable
{
    virtual ~Observer() {}
    //! Called from within Machine::propagate()
    virtual void view(const ElementVoid* elem, const StateBase* state) =0;
};

/**
 * @brief Base class for all simulated elements.
 *
 * Sub-classes of ElementVoid must be registered with Machine::registerElement
 * before they will be found by Machine::Machine().
 */
struct ElementVoid : public boost::noncopyable
{
    /**
     * @brief Construct this element using the provided Config.
     *
     * Base class ctor makes use of Config parameters "name" and "length".
     * "name" is required.  "length" is option, and is 0.0 if omitted.
     *
     * Sub-classes are allowed to require certain parameters to be provided.
     *
     * @throws KeyError       If a required parameter is missing
     * @throws boost::bad_get If a parameter exists, but has the wrong value type
     */
    ElementVoid(const Config& conf);
    virtual ~ElementVoid();

    /** Sub-classes must provide an approprate short description string.
     *  Must match the type name passed to Machine::registerElement().
     */
    virtual const char* type_name() const =0;

    //! Propogate the given State through this Element
    virtual void advance(StateBase& s) =0;

    //! The Config used to construct this element.
    inline const Config& conf() const {return p_conf;}

    const std::string name; //!< Name of this element (unique in its Machine)
    size_t index; //!< Index of this element (unique in its Machine)

    double length; //!< Longitudual length of this element (added to StateBase::pos)

    //! The current observer, or NULL
    Observer *observer() const { return p_observe; }
    /** Add Observer which will inspect the output State of this Element.
     *  Observer instance musy outlive the Element.
     * @param o A new Observer or NULL, will replace any existing pointer.
     */
    void set_observer(Observer *o) { p_observe = o; }

    //! Print information about the element.
    //! level is a hint as to the verbosity expected by the caller.
    virtual void show(std::ostream&, int level) const;

    //! Used by Machine::reconfigure() to avoid re-alloc (and iterator invalidation)
    //! Assumes other has the same type.
    //! Sub-classes must call base class assign()
    //! Come c++11 this can be replaced with a move ctor.
    virtual void assign(const ElementVoid* other ) =0;
private:
    Observer *p_observe;
    Config p_conf;
    friend class Machine;
};

/**
 * @brief The core simulate Machine engine
 *
 * Provides std::vector<ElementVoid*>-like access to individual elements
 *
 * @note A Machine instance is reentrant, but not thread-safe.
 *       Any thread may create a Machine at any time.
 *       However, each instance should be accessed by a single thread.
 */
struct Machine : public boost::noncopyable
{
    /**
     * @brief Construct a new Machine
     * @param c A Config instance, such as that returned by GLPSParser::parse_file().
     */
    Machine(const Config& c);
    ~Machine();

    /** @brief Pass the given bunch State through this Machine.
     *
     * @param S The initial state, will be updated with the final state
     * @param start The index of the first Element the state will pass through
     * @param max The maximum number of elements through which the state will be passed
     * @throws std::exception sub-classes for various errors.
     *         If an exception is thrown then the state of S is undefined.
     */
    void propagate(StateBase* S,
                   size_t start=0,
                   size_t max=-1) const;

    /** @brief Allocate (with "operator new") an appropriate State object
     *
     * @param c Configuration describing the initial state
     * @return A pointer to the new state (never NULL).  The caller takes responsibility for deleteing.
     * @throws std::exception sub-classes for various errors, mostly incorrect Config.
     */
    StateBase* allocState(const Config& c) const;

    //! Allocate new State with empty Config
    //! Equivalent to allocState(Config())
    inline StateBase* allocState() const {
        Config defaults;
        return allocState(defaults);
    }

    //! Fetch Config used to construct this Machine
    inline const Config& conf() const { return p_conf; }

    /**
     * @brief Change the configuration of a single element.
     * @param idx The index of this element
     * @param c The new Config
     *
     * Triggers re-construction of a single element.
     * An optimization to avoid the overhead of reconstructing
     * the entire Machine to change a single element.
     *
     * @code
     * Machine M(...);
     * assert(M.size()>5);
     * Config e5conf(M[5]->conf());
     * @endcode
     */
    void reconfigure(size_t idx, const Config& c);

    //! Return the sim_type string found during construction.
    inline const std::string& simtype() const {return p_simtype;}

    //! The current tracing stream, or NULL.
    inline std::ostream* trace() const {return p_trace;}
    /**
     * @brief Assign new tracing stream.
     * @param v new stream or NULL to clear
     *
     * The trace stream will be written to during Machine::propagate()
     * as a debugging aid.
     */
    void set_trace(std::ostream* v) {p_trace=v;}

private:
    typedef std::vector<ElementVoid*> p_elements_t;

    struct LookupKey {
        std::string name;
        size_t index;
        LookupKey(const std::string& n, size_t i) :name(n), index(i) {}
        bool operator<(const LookupKey& o) const {
            int ord = name.compare(o.name);
            if(ord<0)      return true;
            else if(ord>0) return false;
            else           return index<o.index;
        }
    };

    typedef std::map<LookupKey, ElementVoid*> p_lookup_t;
public:

    //! @return Number of beamline elements
    inline size_t size() const { return p_elements.size(); }

    //! Access a beamline element
    inline ElementVoid* operator[](size_t i) { return p_elements[i]; }
    //! Access a beamline element
    inline const ElementVoid* operator[](size_t i) const { return p_elements[i]; }

    //! Access a beamline element
    inline ElementVoid* at(size_t i) { return p_elements.at(i); }
    //! Access a beamline element
    inline const ElementVoid* at(size_t i) const { return p_elements.at(i); }

    //! Beamline element iterator
    typedef p_elements_t::iterator iterator;
    //! Beamline element iterator (const version)
    typedef p_elements_t::const_iterator const_iterator;

    //! Points to the first element
    iterator begin() { return p_elements.begin(); }
    //! Points to the first element
    const_iterator begin() const { return p_elements.begin(); }

    //! Points just after the last element
    iterator end() { return p_elements.end(); }
    //! Points just after the last element
    const_iterator end() const { return p_elements.end(); }

    //! Find the nth element with the given name
    //! @return NULL on failure
    //! A convienence wrapper around equal_range().
    ElementVoid* find(const std::string& name, size_t nth=0) {
        p_lookup_t::const_iterator low (p_lookup.lower_bound(LookupKey(name, 0))),
                                   high(p_lookup.upper_bound(LookupKey(name, (size_t)-1)));
        size_t i=0;
        for(;low!=high;++low,++i) {
            if(i==nth)
                return low->second;
        }
        return NULL;
    }

    //! iterator for use with equal_range() and equal_range_type()
    typedef value_proxy_iterator<p_lookup_t::iterator> lookup_iterator;

    std::pair<lookup_iterator, lookup_iterator> all_range() {
        return std::make_pair(lookup_iterator(p_lookup.begin()),
                              lookup_iterator(p_lookup.end()));
    }

    //! Return a pair of iterators for the sequence [first, second) of those elements
    //! with the given name.
    std::pair<lookup_iterator, lookup_iterator> equal_range(const std::string& name) {
        p_lookup_t::iterator low (p_lookup.lower_bound(LookupKey(name, 0))),
                             high(p_lookup.upper_bound(LookupKey(name, (size_t)-1)));
        return std::make_pair(lookup_iterator(low),
                              lookup_iterator(high));
    }

    //! Return a pair of iterators for the sequence [first, second) of those elements
    //! with the given type name.
    std::pair<lookup_iterator, lookup_iterator> equal_range_type(const std::string& name) {
        p_lookup_t::iterator low (p_lookup_type.lower_bound(LookupKey(name, 0))),
                             high(p_lookup_type.upper_bound(LookupKey(name, (size_t)-1)));
        return std::make_pair(lookup_iterator(low),
                              lookup_iterator(high));
    }


private:
    p_elements_t p_elements;
    p_lookup_t p_lookup; //!< lookup by element instance name
    p_lookup_t p_lookup_type; //!< lookup by element type name
    std::string p_simtype;
    std::ostream* p_trace;
    Config p_conf;

    typedef StateBase* (*state_builder_t)(const Config& c);
    template<typename State>
    struct state_builder_impl {
        static StateBase* build(const Config& c)
        { return new State(c); }
    };
    struct element_builder_t {
        virtual ~element_builder_t() {}
        virtual ElementVoid* build(const Config& c) =0;
        virtual void rebuild(ElementVoid *o, const Config& c, const size_t idx) =0;
    };
    template<typename Element>
    struct element_builder_impl : public element_builder_t {
        virtual ~element_builder_impl() {}
        ElementVoid* build(const Config& c)
        { return new Element(c); }
        void rebuild(ElementVoid *o, const Config& c, const size_t idx)
        {
            std::auto_ptr<ElementVoid> N(build(c));
            Element *m = dynamic_cast<Element*>(o);
            if(!m)
                throw std::runtime_error("reconfigure() can't change element type");
            m->assign(N.get());
            m->index = idx; // copy index number
        }
    };

    struct state_info {
        std::string name;
        state_builder_t builder;
        typedef std::map<std::string, element_builder_t*> elements_t;
        elements_t elements;
    };

    state_info p_info;

    typedef std::map<std::string, state_info> p_state_infos_t;
    static p_state_infos_t p_state_infos;

    static void p_registerState(const char *name, state_builder_t b);

    static void p_registerElement(const std::string& sname, const char *ename, element_builder_t* b);

public:

    /**
     * @brief Register a new State with the simulation framework.
     *
     * The first step to adding a new sim_type.
     *
     * @param name The new sim_type name
     * @throws std::logic_error if name is already registered
     *
     * @note This method may be called from any thread at any time.
     *
     * @code
     * struct MyState : public StateBase { ... };
     * void myReg() {
     *   Machine::registerState<MyState>("mysimtype");
     *   ...
     * @endcode
     */
    template<typename State>
    static void registerState(const char *name)
    {
        p_registerState(name, &state_builder_impl<State>::build);
    }

    /**
     * @brief Register a new Element type with the simulation framework
     *
     * Add a new element to an existing sim_type (see registerState()).
     *
     * @param sname A sim_type name
     * @param ename The new element type name
     * @throws std::logic_error if sname has not been registered, or if ename is already registered
     *
     * @note This method may be called from any thread at any time.
     */
    template<typename Element>
    static void registerElement(const char *sname, const char *ename)
    {
        p_registerElement(sname, ename, new element_builder_impl<Element>);
    }

    /**
     * @brief Discard all registered State and Element type information.
     *
     * Clears all previous registerations made by registerState() and registerElement().
     *
     * Suggested use case is to call just before process exit
     * so that valgrind doesn't flag these as leaks
     *
     * @note This method may be called from any thread at any time.
     */
    static void registeryCleanup();

    friend std::ostream& operator<<(std::ostream&, const Machine& m);

    struct LogRecord {
        const char * fname;
        unsigned short lnum; // >=64k LoC in one file is already a bug
        unsigned short level;
        std::ostringstream strm;
        LogRecord(const char *fname, unsigned short lnum, unsigned short lvl)
            :fname(fname), lnum(lnum), level(lvl) {}
        ~LogRecord();
        template<typename T>
        LogRecord& operator<<(T v)
        { strm<<v; return *this; }
    private:
        LogRecord(const LogRecord&);
        LogRecord& operator=(const LogRecord&);
    };

    struct Logger {
        virtual ~Logger() {}
        virtual void log(const LogRecord& r)=0;
    };

    static int log_detail;

    static inline bool detail(int lvl) { return log_detail<=lvl; }

    static void set_logger(const boost::shared_ptr<Logger>& p);
private:
    static boost::shared_ptr<Logger> p_logger;
};

#define FLAME_ERROR 40
#define FLAME_WARN  30
#define FLAME_INFO  20
#define FLAME_DEBUG 10
// Super verbose logging, a la the rf cavity element
#define FLAME_FINE  0

#define FLAME_LOG_ALWAYS(LVL) Machine::LogRecord(__FILE__, __LINE__, FLAME_##LVL)

#ifndef FLAME_DISABLE_LOG
#define FLAME_LOG_CHECK(LVL) UNLIKELY(Machine::detail(FLAME_##LVL))
#else
#define FLAME_LOG_CHECK(LVL) (false)
#endif

//! Hook into redirect-able logging.  Use like an std::ostream w/ operator<<
#define FLAME_LOG(LVL) if(FLAME_LOG_CHECK(LVL)) FLAME_LOG_ALWAYS(LVL)

std::ostream& operator<<(std::ostream&, const Machine& m);

//! Register sim_types "Vector" and "TransferMatrix"
void registerLinear();
//! Register sim_type "MomentMatrix"
void registerMoment();

#endif // FLAME_BASE_H
