#ifndef SCSI_BASE_H
#define SCSI_BASE_H

#include <stdlib.h>

#include <ostream>
#include <string>
#include <map>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/call_traits.hpp>

#include "config.h"
#include "util.h"

// Macros:
#define sqr(x)  ((x)*(x))
#define cube(x) ((x)*(x)*(x))

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
    double IonZ_ref,
           IonZ,
           IonEs,
           IonEk,
           IonW_ref, // Total energy of reference particle.
           IonW;     // Total energy of ion.

    virtual void assign(const StateBase& other) =0;

    virtual void show(std::ostream&) const {}

    struct ArrayInfo {
        ArrayInfo() :name(), type(Double), ptr(NULL), ndim(0) {}
        std::string name;
        enum Type {
            Double, Sizet
        } type;
        void *ptr;
        int ndim;
        size_t dim[5];
    };

    /** @brief Introspect named parameter of the derived class
     * @param idx The index of the parameter
     * @param Info mailbox to be filled with information about this paramter
     * @return true if 'idx' is valid and Info was filled in, otherwise false
     *
     * To introspect call this with index increasing from zero until false is returned.
     *
     * @note This method requires that parameter storage be stable for the lifetime of the object
     */
    virtual bool getArray(unsigned idx, ArrayInfo& Info);

    virtual StateBase* clone() const =0;

    //! @private
    //! Mailbox to hold the python interpreter object wrapping us.
    void *pyptr;
protected:
    StateBase(const Config& c);
    struct clone_tag{};
    StateBase(const StateBase& c, clone_tag);
};

inline
std::ostream& operator<<(std::ostream& strm, const StateBase& s)
{
    s.show(strm);
    return strm;
}

struct Observer : public boost::noncopyable
{
    virtual ~Observer() {}
    virtual void view(const ElementVoid*, const StateBase*) =0;
};

struct ElementVoid : public boost::noncopyable
{
    ElementVoid(const Config& conf);
    virtual ~ElementVoid();

    virtual const char* type_name() const =0;

    //! Propogate the given State through this Element
    virtual void advance(StateBase& s) =0;

//    inline const Config& conf() const {return p_conf;}
    inline Config& conf() {return p_conf;}

    const std::string name; //!< Name of this element (unique in its Machine)
    const size_t index; //!< Index of this element (unique in its Machine)

    Observer *observer() const { return p_observe; }
    /** Add Observer which will inspect the output State of this Element.
     *  Observer instance musy outlive the Element.
     * @param o A new Observer or NULL, will replace any existing pointer.
     */
    void set_observer(Observer *o) { p_observe = o; }

    virtual void show(std::ostream&) const;

    //! @internal
    //! Used by Machine::reconfigure()
    //! Assumes other has the same type
    virtual void assign(const ElementVoid* other ) =0;
private:
    Observer *p_observe;
    Config p_conf;
    friend class Machine;
};

//std::ostream& operator<<(std::ostream& strm, const ElementVoid& s)
//{
//    s.show(strm);
//    return strm;
//}

struct Machine : public boost::noncopyable
{
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

    void reconfigure(size_t idx, const Config& c);

    inline const std::string& simtype() const {return p_simtype;}

    inline std::ostream* trace() const {return p_trace;}
    void set_trace(std::ostream* v) {p_trace=v;}

    typedef std::vector<ElementVoid*> p_elements_t;
    typedef std::map<std::string, ElementVoid*> p_lookup_t;

    inline size_t size() const { return p_elements.size(); }

    typedef p_elements_t::iterator iterator;
    typedef p_elements_t::const_iterator const_iterator;

    iterator begin() { return p_elements.begin(); }
    const_iterator begin() const { return p_elements.begin(); }

    iterator end() { return p_elements.end(); }
    const_iterator end() const { return p_elements.end(); }

    inline ElementVoid* operator[](size_t i) { return p_elements[i]; }
    inline const ElementVoid* operator[](size_t i) const { return p_elements[i]; }
private:
//    p_elements_t p_elements;
//    p_lookup_t p_lookup;
//    std::string p_simtype;
//    std::ostream* p_trace;

    typedef StateBase* (*state_builder_t)(const Config& c);
    template<typename State>
    struct state_builder_impl {
        static StateBase* build(const Config& c)
        { return new State(c); }
    };
    struct element_builder_t {
        virtual ~element_builder_t() {}
        virtual ElementVoid* build(const Config& c) =0;
        virtual void rebuild(ElementVoid *o, const Config& c) =0;
    };
    template<typename Element>
    struct element_builder_impl : public element_builder_t {
        virtual ~element_builder_impl() {}
        ElementVoid* build(const Config& c)
        { return new Element(c); }
        void rebuild(ElementVoid *o, const Config& c)
        {
            std::auto_ptr<ElementVoid> N(build(c));
            Element *m = dynamic_cast<Element*>(o);
            if(!m)
                throw std::runtime_error("reconfigure() can't change element type");
            m->assign(N.get());
        }
    };

    struct state_info {
        std::string name;
        state_builder_t builder;
        typedef std::map<std::string, element_builder_t*> elements_t;
        elements_t elements;
    };

//    state_info p_info;

    typedef std::map<std::string, state_info> p_state_infos_t;
    static p_state_infos_t p_state_infos;

    static void p_registerState(const char *name, state_builder_t b);

    static void p_registerElement(const std::string& sname, const char *ename, element_builder_t* b);

public:
    p_elements_t p_elements;
    p_lookup_t p_lookup;
    std::string p_simtype;
    std::ostream* p_trace;
    state_info p_info;

    template<typename State>
    static void registerState(const char *name)
    {
        p_registerState(name, &state_builder_impl<State>::build);
    }

    template<typename Element>
    static void registerElement(const char *sname, const char *ename)
    {
        p_registerElement(sname, ename, new element_builder_impl<Element>);
    }

    //! Discard all registered State and Element type information
    //! Suggested use case is to call just before process exit
    //! so that valgrind doesn't flag these as leaks
    static void registeryCleanup();

    friend std::ostream& operator<<(std::ostream&, const Machine& m);
};

std::ostream& operator<<(std::ostream&, const Machine& m);

void registerLinear();
void registerMoment();
void registerMoment2();

#endif // SCSI_BASE_H
