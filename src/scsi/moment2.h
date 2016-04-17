#ifndef SCSI_MOMENT_H
#define SCSI_MOMENT_H

#include <ostream>
#include <iomanip>
#include <math.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "base.h"

/** @brief Simulation state which include only a matrix
 */

class CavDataType {
// Cavity on-axis longitudinal electric field vs. s.
public:
    std::vector<double> s,     // s coordinate [m]
                        Elong; // Longitudinal Electric field [V/m].

    void RdData(const std::string&);
    void show(std::ostream&, const int) const;
    void show(std::ostream&) const;
};


extern std::stringstream   CavTLMstream2[];
extern CavDataType         CavData2[];

struct Moment2State : public StateBase
{
    enum {maxsize=7};
    enum param_t {
        PS_X, PS_PX, PS_Y, PS_PY, PS_S, PS_PS,
        PS_QQ // ???
    };

    Moment2State(const Config& c);
    virtual ~Moment2State();

    typedef boost::numeric::ublas::vector<double,
                    boost::numeric::ublas::bounded_array<double, maxsize>
    > vector_t;

    typedef boost::numeric::ublas::matrix<double,
                    boost::numeric::ublas::row_major,
                    boost::numeric::ublas::bounded_array<double, maxsize*maxsize>
    > matrix_t;

    void assign(const StateBase& other);

    virtual void show(std::ostream& strm) const;

    double pos;        // absolute longitudinal position at end of Element
    double IonZ;       // Charge state.
    double Ekinetic;   // kinetic energy of reference particle
                       // actual is Ekinetic + moment0[6]
    double SampleIonK;

    double sync_phase, // synchotron phase
           FyAbs,      // Betatron phase for reference particle.
           EkState,
           Fy_absState;

    double gamma,      // (Erest+Ekinetic)/Erest
           beta,       // sqrt(1e0-1e0/(gamma*gamma))
           bg0,        // Initial beta*gamma.
           bg1;        // Beta*gamma for reference particle.

    vector_t moment0;
    matrix_t state; // TODO: better name

    virtual bool getArray(unsigned idx, ArrayInfo& Info);

    virtual Moment2State* clone() const {
        return new Moment2State(*this, clone_tag());
    }

protected:
    Moment2State(const Moment2State& o, clone_tag);
};

/** @brief An Element which propagates the statistical moments of a bunch
 *
 */
struct Moment2ElementBase : public ElementVoid
{
    typedef Moment2State state_t;

    Moment2ElementBase(const Config& c);
    virtual ~Moment2ElementBase();

    virtual void advance(StateBase& s);

    virtual void recompute_matrix(state_t&);

    virtual void show(std::ostream& strm) const;

    typedef boost::numeric::ublas::matrix<double> value_t;

    double length;
    double FSampLength; //!< sample (rf) clock wavelength
    /** Fy += phase_factor/beta
     * phase_factor := L*2*pi*Fsamp/C
     */
    double phase_factor;
    double Erest; //!< rest energy of particle species

    double last_Kenergy_in, last_Kenergy_out;
    //! final transfer matrix
    value_t transfer;
    value_t transfer_raw, misalign, misalign_inv;

    virtual void assign(const ElementVoid *other);

protected:
    // scratch space to avoid temp. allocation in advance()
    // An Element can't be shared between multiple threads
    //TODO: non-const advance()
    mutable state_t::matrix_t scratch;
};

#endif // SCSI_MOMENT_H
