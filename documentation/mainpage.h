
/**
@mainpage FLAME documentation

A discrete accelerator simulation engine

@tableofcontents

@subpage abstract

@subpage latparse

@subpage glpsexprs

@subpage typesparams

@subpage apiusage

Source doc repository https://github.com/frib-high-level-controls/FLAME

*/

// =====================================================================================

/**
@page abstract Abstraction

@section abstractoverview Overview

The core concepts of the simulation engine are the Machine, Element, State, and Config.

A Config is a container for key/value pairs.
This is the interface through which Elements are specialized (eg. length=2 vs. length=3).
The lattice file parser populates a Config, which may then be used to construct a Machine.

A Machine represents a ordered list of Elements.

@diafile elements.dia Elements in a Machine.

A State represents "the beam", which is to say, a particle or bunch of particles
which pass from Element to Element through a Machine.
The operation of the simulation is for an Element to transform a State.
This new (output) state is then passed as input to the next element.

@diafile elements-states.dia Propagation of State through a list of Elements.

In code, passing State between Elements is implemented in Machine::propagate.
Transformation of State is encapsulated in ElementVoid::advance.

The precise meanings of State and Element are governed by the "sim_type" key.
When a Machine is constructed, the "sim_type" is used to select particular
sub-classes of ElementVoid and StateBase.

For example, "sim_type=Vector" selects VectorState and a sub-class of LinearElementBase depending on the elements "type" key.
"type=drift" selects ElementDrift.

These selections are made from a global table which is populated by Machine::registerState and Machine::registerElement.

@section config Configuration

A simulation is configured by passing a Config to Machine::Machine.
A Config can be obtained in several ways:

1. Lattice parser GLPSParser::parse()
2. A Python dictionary
3. Construct an empty Config and Config::set all necessary keys.

At minimum, Machine::Machine requires that a Config has the two keys.
"sim_type" which is a string, and "elements" which is a vector<Config>.
Each element must in turn have "type" as a string.

By convention the nested Config for an element should use the
enclosing Config as an enclosing scope ( Config::new_scope() ).
*/

// =====================================================================================

/**
@page typesparams sim_types, Element types, and Parameter names

@section simtypelist sim_types

At present three "sim_type"s are defined: Vector, TransferMatrix, and MomentMatrix.

@subsection simvector sim_type=Vector

Selects use of VectorState as the state type.
Elements derive from \ref LinearElementBase.

The vector VectorState::state defaults to all zeros.

The propagation step LinearElementBase::advance multiples
the an Element's transfer matrix (LinearElementBase::transfer)
with the vector VectorState::state.

\f{eqnarray*}{ state = Transfer \cdot state \f}

@subsection simmatrix sim_type=TransferMatrix

Selects use of MatrixState as the state type.
Elements derive from \ref LinearElementBase.

The matrix MatrixState::state defaults to an identity matrix.

The propagation step LinearElementBase::advance multiples
the an Element's transfer matrix (LinearElementBase::transfer)
with the matrix MatrixState::state.

\f{eqnarray*}{ State = Transfer \cdot State \f}

@subsection simmoment sim_type=MomentMatrix

Selects use of MomentState as the state type.
Elements derive from \ref MomentElementBase

The propagation step MomentElementBase::advance multiples
the an Element's transfer matrix (MomentElementBase::transfer)
with the matrix MomentState::state, and also by the vector MomentState::moment0.

\f{eqnarray*}{
  State &=& Transfer \cdot State \cdot Transfer^t \\
  moment0 &=& Transfer \cdot moment0
\f}

@subsection simelements Element Types

This section lists all element types, and lists which "sim_type"s each is defined for.

@subsubsection elementsource source

The purpose of the "source" element is to overwrite the State.

Each instance holds an internal State.
The advance() method simply assigns this internal State to the state being propagate()d.

The parameters accepted are the same as the associated State type.

@code
sim_type = "Vector";
elemname: source, initial = [0, 0, 0, 0, 0, 0, 0];

sim_type = "TransferMatrix";
elemname: source, initial = [1, 0, 0, 0, 0, 0, 0,
                             0, 1, 0, 0, 0, 0, 0,
                             0, 0, 1, 0, 0, 0, 0,
                             0, 0, 0, 1, 0, 0, 0,
                             0, 0, 0, 0, 1, 0, 0,
                             0, 0, 0, 0, 0, 1, 0,
                             0, 0, 0, 0, 0, 0, 1];

sim_type = "MomentMatrix";
elemname: source, initial = [1, 0, 0, 0, 0, 0, 0,
                             0, 1, 0, 0, 0, 0, 0,
                             0, 0, 1, 0, 0, 0, 0,
                             0, 0, 0, 1, 0, 0, 0,
                             0, 0, 0, 0, 1, 0, 0,
                             0, 0, 0, 0, 0, 1, 0,
                             0, 0, 0, 0, 0, 0, 1],
                  moment0 = [0, 0, 0, 0, 0, 0, 0];
@endcode

@subsubsection elementdrift drift

A drift section.

@htmlonly
<table class="param">
<thead><tr><th>Name</th><th>Default</th><th>Desc.</th></tr></thead>
<tbody>
<tr><td>L</td><td>None</td><td>Length</td></tr>
</tbody>
</table>
@endhtmlonly

Supported by all "sim_type"s.

@code
elemname: drift, L = 0.1;
@endcode

@subsubsection elementmarker marker

Equivalent to a zero length drift section.  No parameters.

Supported by all "sim_type"s.

@code
elemname: marker;
@endcode

@subsubsection elementsbend sbend

Sector bend magnet.

@htmlonly
<table class="param">
<thead><tr><th>Name</th><th>Default</th><th>Desc.</th></tr></thead>
<tbody>
<tr><td>L</td><td>None</td><td>Length</td></tr>
<tr><td>phi</td><td>None</td><td>Bend angle</td></tr>
<tr><td>K</td><td>0.0</td><td>Strength</td></tr>
</tbody>
</table>
@endhtmlonly

Supported by all "sim_type"s.

@code
elemname: sbend, L = 0.1, phi = 3.14159/16, K = 10;
@endcode

@subsubsection elementquad quadrupole

Magnetic quadrupole.

@htmlonly
<table class="param">
<thead><tr><th>Name</th><th>Default</th><th>Desc.</th></tr></thead>
<tbody>
<tr><td>L</td><td>None</td><td>Length</td></tr>
<tr><td>K</td><td>0.0</td><td>Strength.  K&gt;0 focusing in horizontal. K&lt;0 focusing in vertical. </td></tr>
</tbody>
</table>
@endhtmlonly

Supported by all "sim_type"s.

@code
elemname: quadrupole, L = 0.1, K = 10;
@endcode

@subsubsection elementsol solenoid

Solenoid magnet.

@htmlonly
<table class="param">
<thead><tr><th>Name</th><th>Default</th><th>Desc.</th></tr></thead>
<tbody>
<tr><td>L</td><td>None</td><td>Length</td></tr>
<tr><td>K</td><td>0.0</td><td>Strength.</td></tr>
</tbody>
</table>
@endhtmlonly

Supported by all "sim_type"s.

@code
elemname: solenoid, L = 0.1, K = 10;
@endcode

@subsubsection elementrf rfcavity

TODO

@subsubsection elementstrip stripper

Placeholder.  Equivalent to marker.

@subsubsection elementedipole edipole

Placeholder.  Equivalent to marker.

@subsubsection elementgeneric generic

Element for which transfer matrix is directly specified.

@htmlonly
<table class="param">
<thead><tr><th>Name</th><th>Default</th><th>Desc.</th></tr></thead>
<tbody>
<tr><td>transfer</td><td>None</td><td>Flattened transfer matrix</td></tr>
</tbody>
</table>
@endhtmlonly

Supported by all "sim_type"s.

@code
elemname: generic, transfer = [1, 0, 0, 0, 0, 0, 0,
                               0, 1, 0, 0, 0, 0, 0,
                               0, 0, 1, 0, 0, 0, 0,
                               0, 0, 0, 1, 0, 0, 0,
                               0, 0, 0, 0, 1, 0, 0,
                               0, 0, 0, 0, 0, 1, 0,
                               0, 0, 0, 0, 0, 0, 1];
@endcode

*/

// =====================================================================================

/**
@page apiusage API Usage

@section apiusagecpp C++ Simulation API Usage

The example @subpage examples_sim_cpp demonstrates the basic process of running
a simulation against an existing sim_type and set of Elements.
This is a simplified version of cli.cpp, which is compiled as the run_uscsi executable.

Most user code will only need base.h and config.h.

@snippet sim.cpp includes

Before constructing a Machine it is necessary to register any "sim_type"s to be used.
This will typically be the first thing a program does.

@snippet sim.cpp Register sim_types

Running the simulation with the C++ API begins by populating a Config.
This will often by done using the lattice parser (GLPSParser).

Parsing will succeed as long as the file is syntactically correct.
For example, "sim_type" and element types are not validated.

@snippet sim.cpp Parse lattice file

Next construct a Machine.
It is at this step that "sim_type" and element types are validated.

@snippet sim.cpp Construct Machine

Now allocate an approprate StateBase based on the "sim_type".
More than one StateBase may be allocated and used with a single Machine.
The allocated StateBase must only be used with the Machine which allocated it.

This StateBase sub-class is be initialized based a Config.
In this example an empty Config is provided, which will
leave the state with some "sim_type" specific defaults
("sim_type=Vector" defaults to all zeros).
Generally this will only be done when the lattice begins with a "source" element,
which will overwrite the state based on its Config.

Alternately, a Config can be provided to Machine::allocState() to
provide non-default values.

@snippet sim.cpp Allocate State

The Machine::propagate() method is used to run the simulation.
When propagate() returns, the StateBase has been modified to
reflect the final state.

By default propagate() starts with the first element
and continues through the last element.

@snippet sim.cpp Run simulation

Before exiting some cleanup may be to clear the registry of "sim_type"s.

@snippet sim.cpp Deinit

@section apiusagepy Python Simulation API Usage

The example @subpage examples_sim_cpp demonstrates the basic process of running
a simulation using the python API.

Importing "uscsi" also registers a standard set of "sim_type"s.

@snippet sim.py includes

Read an parse a lattice file and construct a Machine.

@snippet sim.py Parse lattice file

Allocate a State using an empty Config.

@snippet sim.py Allocate State

Propagate through all elements.

@snippet sim.py Run simulation

@section apicreatesim Defining a sim_type

The example defines a new sim_type and two Elements.
The simulation itself will be trivial.

The full example can be found in @subpage examples_customsim_cpp

@subsection apicreatestate Define state struct

The first task is to define the state structure.
This contains the variables which define the state of one "bunch".
In this case "x" and "xv", which represent a transverse velocity and relative position.
In addition the absolute logitudinal StateBase::pos is inherited.

@snippet customsim.cpp StateDef

The member variables of State1D are initialized from a Config.

@snippet customsim.cpp StateInit

The remainder of State1D is boilerplate which is repeated for each member variable.

@subsection apicreatesource Define source element type

Now define the "source" element type.
By convention this is an element which simply overwrites it's input state with
predefined values (eg. from lattice file).
So it is simplest to give this struct it's own internal State1D instance.

@snippet customsim.cpp ElemSrcDef

Initialization of the element can also initialize this internal State1D from the same Config.

@snippet customsim.cpp ElemSrcInit

So the same "x" and "xv" used to initialize the state can be given to a source element.

@code
srcname : source, x=0, xv=0.1, L=0.5;
@endcode

With this lattice file entry, the "x" and "xv" are used to initialize the internal State1D.
"L" is used to initialize ElementVoid::length.

Next define the simulated action of this element with its advance() method.
This method should modify it's input State1D.
In this case, the source element simply overwrites the input state using its internal State1D.

@snippet customsim.cpp ElemSrcAdvance

@subsection apicreategeneric Define another element type

Now define another "generic" element type which transforms the state in a more interesting way.

@snippet customsim.cpp ElemGenericDef

@snippet customsim.cpp ElemGenericInit

Here we see that the name given to Config::get need not match a c++ variable name.
Also that some calculations can be done once, and the result stored in a class member variable for later re-use.

@snippet customsim.cpp ElemGenericAdvance

The advance() function incrementally updates the state.

@subsection apicreateregister Registration

So far three types have been defined: State1D, Element1DSource, and Element1DGeneric.
If nothing further is done there will be no way to make use of this code.
As these definitions appear in an anonymous C++ namespace a modern compiler is able to issue a warn
that this code is unreachable.

In order for this code to be reachable it must be tied into the FLAME framework in some way.
This is accomplished with Machine::registerState and Machine::registerElement.

By convention a function "register...()" is defined, which must be called exactly once
before Machine can make use of these definitions.

@snippet customsim.cpp register

Now to make use Simple1D we expand on @subpage examples_sim_cpp

@snippet customsim.cpp main

*/

// =====================================================================================

/**
@page examples_sim_cpp examples/sim.cpp

@include sim.cpp
*/

// =====================================================================================

/**
@page examples_customsim_cpp examples/customsim.cpp

@include customsim.cpp
*/

// =====================================================================================

/**
@page examples_sim_py examples/sim.py

@include sim.py
*/

// =====================================================================================

/**
@page latparse Lattice file syntax

@section latover Overview

The lattice file will contain one or more element definitions and one or more beamline definitions.

An element definition takes the form:

@verbatim
<instance_name> : <element_type>, <key>=<value>, ... ;
@endverbatim

For example

@verbatim
dr001 : drift, L=0.1;
@endverbatim

A beamline definition is constructed from zero or more element or beamline names.

@verbatim
<instance_name> : LINE = (<element>, ...) ;
@endverbatim

For example

@verbatim
line1 : LINE = (dr001, dr001);
@endverbatim

Entries in a beamline definition may be expressions using element and beamline names.
For example

@verbatim
d1 : drift, L=1;
d2 : drift, L=2;
line1 : LINE = (d1, d2, d1);
line2 : LINE = (d1, line1, 2*line1, -line1);
@endverbatim

A "USE" statement may be given to select which beamline is expanded and used.
If no "USE" statement is given then the last "LINE" is used.

@verbatim
USE : line1;
@endverbatim

@section latexpr Variables and Expressions

In addition to element and beamline definitinos, arbitrary variables may be defined.
String, scalar and vector floating point values may be assigned to a variable.

@verbatim
var1 = 1.0;
var2 = "hello";
var3 = [1,2,3,4];
@endverbatim

Expressions of these types may also be assigned.
These expressions may use the values of previously defined variables,
as well as operations and functions.

@verbatim
var4 = 1+var1;
var5 = cos(var4);
@endverbatim

The full list of valid expressions is given on the @subpage glpsexprs page.

@section latgram Grammar definition

@verbatim

# unquoted word
KEYWORD := [A-Za-z]([A-Za-z0-9_:]*[A-Za-z0-9_])?
# floating point number
NUM := [0-9]+(\.[0-9]*)?([eE][+-]?[0-9]+)?
# quoted string
STR := ".*"

# comment (ignored)
COMMENT := #.*\n

assignment : KEYWORD '=' expr ';'
element    : KEYWORD ':' KEYWORD properties ';'
line       : KEYWORD ':' KEYWORD '=' '(' line_list ')' ';'
func       : KEYWORD '(' expr ')' ';'
command    : KEYWORD ';'

properties :
           | ',' property properties
property : KEYWORD '=' expr

line_list :
          | expr
          | expr ',' line_list

expr : NUM
     | vector
     | STR
     | KEYWORD
     | expr '+' expr
     | expr '-' expr
     | expr '*' expr
     | expr '/' expr
     | '-' expr
     | '(' expr ')'
     | KEYWORD '(' expr ')'

expr_list :
          | expr
          | expr ',' expr_list

vector : '[' expr_list ']'

entry : assignment | element | line | func | command

file :
     | entry file

@endverbatim
*/

// =====================================================================================

/**
@page glpsexprs Expressions

@see parse_context::parse_context

@section exprops Operators

@htmlonly
<table class="param">
<thead><tr><th>Op.</th><th>Desc.</th></tr></thead>
<tbody>
<tr><td>NUM := - NUM</td><td>Floating point negation</td></tr>
<tr><td>NUM := NUM + NUM</td><td>Floating point addition</td></tr>
<tr><td>NUM := NUM - NUM</td><td>Floating point subtraction</td></tr>
<tr><td>NUM := NUM * NUM</td><td>Floating point multiplication</td></tr>
<tr><td>NUM := NUM / NUM</td><td>Floating point division.  Divide by zero will trigger a parser error.</td></tr>
</tbody>
</table>
@endhtmlonly

@htmlonly
<table class="param">
<thead><tr><th>Op.</th><th>Desc.</th></tr></thead>
<tbody>
<tr><td>NUM := sin(NUM)</td><td>Floating point trig. functions</td></tr>
<tr><td>NUM := cos(NUM)</td><td></td></tr>
<tr><td>NUM := tan(NUM)</td><td></td></tr>
<tr><td>NUM := asin(NUM) or arcsin(NUM)</td><td></td></tr>
<tr><td>NUM := acos(NUM) or arccos(NUM)</td><td></td></tr>
<tr><td>NUM := atan(NUM) or arctan(NUM)</td><td></td></tr>
</tbody>
</table>
@endhtmlonly

@htmlonly
<table class="param">
<thead><tr><th>Op.</th><th>Desc.</th></tr></thead>
<tbody>
<tr><td>NUM := rad2deg(NUM)</td><td>Convert between radians and degrees</td></tr>
<tr><td>NUM := deg2rad(NUM)</td><td></td></tr>
</tbody>
</table>
@endhtmlonly

@htmlonly
<table class="param">
<thead><tr><th>Op.</th><th>Desc.</th></tr></thead>
<tbody>
<tr><td>STR = file(STR)</td><td>Normalize file path.  Accepts a possibly relative path, returns the canonical path (no '.' or '..')</td></tr>
<tr><td>STR = h5file(STR)</td><td>Normalize h5file path.  same as file() while also ignoring any trailing "/group/name".</td></tr>
</tbody>
</table>
@endhtmlonly

*/

// =====================================================================================
