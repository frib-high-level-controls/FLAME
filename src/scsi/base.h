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

#include "util.h"

/** @brief A wrapper around map<string,any> to hold configuration information
 */
struct Config
{
    bool has(const std::string& s) const;

    boost::any getAny(const std::string& s) const;
    boost::any getAny(const std::string& s, const boost::any& def) const;
    //! Set value with implicit cast to @class boost::any
    void setAny(const std::string& s, const boost::any& val);

    /** @brief Fetch value assocated with name and cast to type T
     * @throws key_error if name is not valid
     * @throws boost::bad_any_cast if the value associated with name can't be cast to type T
     */
    template<typename T>
    T get(const std::string& name) const
    {
        return boost::any_cast<T>(getAny(name));
    }
    /** @brief @see get()
     * @returns the value assocated with name, or the provided default value
     */
    template<typename T>
    T get(const std::string& s, typename boost::call_traits<T>::param_type def) const
    {
        try{
            return boost::any_cast<T>(getAny(s));
        }catch(boost::bad_any_cast&){
        }catch(key_error&){
        }
        return def;
    }
    //! Set value with explicit type
    template<typename T>
    void set(const std::string& s, typename boost::call_traits<T>::param_type V)
    {
        setAny(s, V);
    }

    typedef std::map<std::string, boost::any> map_t;
private:
    map_t p_props;

    friend std::ostream& operator<<(std::ostream&, const Config& c);
};

std::ostream& operator<<(std::ostream&, const Config& c);

/** @brief The abstract base class for all simulation state objects.
 *
 * Represents the state of some bunch of particles moving through a @class Machine
 */
struct StateBase : public boost::noncopyable
{
    virtual ~StateBase();

    //! Index of @class ElementBase in @class Machine which will follow this one.
    //! May be altered by ElementBase::advance()
    size_t next_elem;

    virtual void show(std::ostream&) const {}

    struct ArrayInfo {
        ArrayInfo() :name(), ptr(NULL), ndim(0) {}
        std::string name;
        double *ptr;
        int ndim;
        size_t dim[5];
    };

    /** @brief Introspect named parameter of the derived class
     * @param idx The index of the parameter
     * @param Info mailbox to be filled with information about this paramter
     * @return true if the parameter and Info was filled in, otherwise false
     *
     * To introspect call this with index increasing from zero until false is returned.
     */
    virtual bool getArray(unsigned idx, ArrayInfo& Info) {return false;}

    //! Mailbox to hold the python interpreter object wrapping us.
    void *pyptr;
protected:
    StateBase(const Config& c) :next_elem(0), pyptr(0) {}
};

struct ElementVoid : public boost::noncopyable
{
    ElementVoid(const Config& conf);
    virtual ~ElementVoid();

    virtual const char* type_name() const =0;

    //! Called after all Elements are constructed
    virtual void peek(const std::vector<ElementVoid*>&) {}

    //! Propogate the given State through this Element
    virtual void advance(StateBase& s) const =0;

    inline const Config& conf() const {return p_conf;}

    const std::string name; //!< Name of this element (unique in its @class Machine)
    const size_t index; //!< Index of this element (unique in its @class Machine)

    virtual void show(std::ostream&) const;

private:
    const Config p_conf;
};

struct Machine : public boost::noncopyable
{
    Machine(const Config& c);
    ~Machine();

    void propogate(StateBase* S,
                   size_t start=0,
                   size_t max=-1) const;

    StateBase* allocState(Config& c) const;

    typedef std::vector<ElementVoid*> p_elements_t;
private:
    p_elements_t p_elements;

    typedef StateBase* (*state_builder_t)(const Config& c);
    typedef ElementVoid* (*element_builder_t)(const Config& c);
    template<typename State>
    struct state_builder_impl {
        static StateBase* build(const Config& c)
        { return new State(c); }
    };
    template<typename Element>
    struct element_builder_impl {
        static ElementVoid* build(const Config& c)
        { return new Element(c); }
    };

    struct state_info {
        std::string name;
        state_builder_t builder;
        typedef std::map<std::string, element_builder_t> elements_t;
        elements_t elements;
    };

    state_info p_info;

    typedef std::map<std::string, state_info> p_state_infos_t;
    static p_state_infos_t p_state_infos;

    static void p_registerState(const char *name, state_builder_t b);

    static void p_registerElement(const std::string& sname, const char *ename, element_builder_t b);

public:

    template<typename State>
    static void registerState()
    {
        p_registerState(State::type_name(), &state_builder_impl<State>::build);
    }

    template<typename Element>
    static void registerElement(const char *name)
    {
        typedef typename Element::state_t state_t;
        p_registerElement(state_t::type_name(), name, &element_builder_impl<Element>::build);
    }

    friend std::ostream& operator<<(std::ostream&, const Machine& m);
};

std::ostream& operator<<(std::ostream&, const Machine& m);

void registerLinear();

#endif // SCSI_BASE_H
