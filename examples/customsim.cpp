#include <iostream>

#include "flame/config.h"
#include "flame/base.h"

# define C0           2.99792458e8

// Good practice to enclose all local definitions in an anonymous namespace
// equivalent of static functions in plain C.
namespace {

//! [StateDef]
struct State1D : public StateBase
{
    double x,  // transverse position
           xv; // transverse veclocity
//! [StateDef]

    virtual void assign(const StateBase& other)
    {
        StateBase::assign(other);
        const State1D& ST = static_cast<const State1D&>(other);
        x = ST.x;
        xv= ST.xv;
    }

    virtual void show(std::ostream &strm) const
    {
        strm<<pos<<" ["<<x<<", "<<xv<<"]\n";
    }

    // allow introspection of state variable by eg. Python
    virtual bool getArray(unsigned idx, ArrayInfo& Info)
    {
        unsigned I=0;
        if(idx==I++) {
            Info.name = "x";
            Info.ptr = &x;
            // remaining defaults ok for scalar double
        } else if(idx==I++) {
            Info.name = "xv";
            Info.ptr = &xv;
        } else {
            return StateBase::getArray(idx-I, Info); // check w/ base class for any remaining
        }
        return true;
    }

    virtual StateBase* clone() const {
        return new State1D(*this, clone_tag());
    }

    virtual ~State1D() {}

    // initialize state from config
    //! [StateInit]
    State1D(const Config& c)
        :StateBase(c)
        ,x(c.get<double>("x", 0.0))   // m
        ,xv(c.get<double>("xv", 0.0)) // m/s
    {}
    //! [StateInit]
    State1D(const State1D& O, clone_tag t)
        :StateBase(O, t)
        ,x(O.x)
        ,xv(O.xv)
    {}
};

//! [ElemSrcDef]
struct Element1DSource : public ElementVoid
{
    State1D initial;
//! [ElemSrcDef]

//! [ElemSrcInit]
    Element1DSource(const Config& c)
        :ElementVoid(c)
        ,initial(c)
    {}
//! [ElemSrcInit]

//! [ElemSrcAdvance]
    virtual void advance(StateBase& s)
    {
        s.assign(initial);
        s.pos += length; // source element ususaly has zero length, but not required
    }
//! [ElemSrcAdvance]

    virtual void assign(const ElementVoid* other )
    {
        ElementVoid::assign(other);
        const Element1DSource *O = static_cast<const Element1DSource*>(other);
        initial.assign(O->initial);
    }

    virtual const char* type_name() const { return "source"; }
};

//! [ElemGenericDef]
struct Element1DGeneric : public ElementVoid
{
    double xa, // transverse acceleration
           tt; // logitudinal transit time
//! [ElemGenericDef]

//! [ElemGenericInit]
    Element1DGeneric(const Config& c)
        :ElementVoid(c)
        ,xa(c.get<double>("A", 0.0)) // m/s2
    {
        // length will never change, so only need to compute transit time once
        tt = length/C0; // sec.
    }
//! [ElemGenericInit]

//! [ElemGenericAdvance]
    virtual void advance(StateBase& s)
    {
        State1D &ST = static_cast<State1D&>(s); // safe since sim_type=Simple1D will only use State1D

        ST.pos += length;
        ST.xv  += xa*tt;
        ST.x   += xa*tt*tt;
    }
//! [ElemGenericAdvance]

    virtual void assign(const ElementVoid* other )
    {
        ElementVoid::assign(other);
        const Element1DGeneric *O = static_cast<const Element1DGeneric*>(other);
        xa = O->xa;
        tt = O->tt;
    }

    virtual const char* type_name() const { return "generic"; }
};

} // end anon namespace

//! [register]
void register1D()
{
    Machine::registerState<State1D>("Simple1D");

    Machine::registerElement<Element1DSource>("Simple1D", "source");
    Machine::registerElement<Element1DGeneric>("Simple1D", "generic");
}
//! [register]

//! [main]
static const char lattice[] = ""
"sim_type = \"Simple1D\";\n"
"src   : source, x = 1e-5, xv = 1e-6;\n"
"elem1 : generic, A = 1,  L = 10;\n"
"elem2 : generic, A = -2, L = 10;\n"
"example : LINE = (src, elem1, elem2);\n"
;

int main(int argc, char *argv[])
{
    try {
        register1D();

        std::unique_ptr<Config> conf;
        {
            GLPSParser parser;
            conf.reset(parser.parse_byte(lattice, sizeof(lattice)-1));
        }

        /* conf now contains (using python notation)
         * {"sim_type":"Simple1D",
         *  "elements": [
         *    {"type":"source",  "name":"src", "x":1e-5, "xv":1e-6},
         *    {"type":"generic", "name":"elem1", "A":1.0,  "L":10.0},
         *    {"type":"generic", "name":"elem2", "A":-2.0, "L":10.0},
         *  ],
         * }
         */

        Machine mymachine(*conf);
        /* During Machine construction
         * sim_type=Simple1D with type=source selects Element1DSource, which is constructed once
         *
         * sim_type=Simple1D with type=source selects Element1DGeneric, which is constructed twice.
         *  first with A=1, then a second time with A=-2
         *
         * mymachine now has 3 elements.  mymachine.size()==3
         */

        mymachine.set_trace(&std::cout); // print intermediates

        std::unique_ptr<StateBase> thestate(mymachine.allocState());

        // propagate through the source element to initialize the state based on
        // values of "x" and "xv" from the lattice
        mymachine.propagate(thestate.get(), 0, 1);

        std::cout<<"Source state "<<*thestate<<"\n";

        // propagate through remaining elements
        mymachine.propagate(thestate.get(), 1);

        std::cout<<"Final state "<<*thestate<<"\n";

        Machine::registeryCleanup();

        return 0;
    } catch(std::exception& e) {
        std::cerr<<"Error: "<<e.what()<<"\n";
        return 1;
    }
}
//! [main]
