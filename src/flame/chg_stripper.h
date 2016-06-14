#ifndef CHG_STRIPPER_H
#define CHG_STRIPPER_H

#include <boost/numeric/ublas/matrix.hpp>

#include "moment2.h"

typedef Moment2State state_t;


// Phase space dimension; including vector for orbit/1st moment.
# define PS_Dim Moment2State::maxsize // Set to 7; to include orbit.


// Charge stripper parameters.

/* SRIM Fitting Model:

     f(t, E)=(AN*t^(u-1))*exp(-l*t)*exp(-0.5*((E-Eo)/E1)^2)),

     t->xdata(:,1) E->xdata(:,2)

     AN->x(1),l->x(2), u->x(3), E_0->x(4), E_1->x(5)

   where t is scattering angle, and E is particle energy.

   Stripper Thin Lens Model:

   RMS standard deviation of <> after stripper can be calculated as:

     std_e = sqrt(std_i^2+std_s^2+std_f^2)

   std_i is the initial sqrt(<>), std_s is the variance introduced by interacting with perfect stripper,
   std_f is the variance introduced by stripper foil thickness variation.
   std_s is E1Para, std_f is sqrt(1/3)*20/100*3.0*abs(-0.10681), sqrt(1/3) is introduced by assuming uniform
   distribution of stripper thickness variation, 20 is thickness variation in %,  3.0 is foil thickness in um,
   -0.10681 is thickness dependence of E0 after stripper.                                                       */
static
const double Stripper_IonZ      = 78.0/238.0,
             Stripper_IonMass   = 238.0,
             Stripper_IonProton = 92.0,
             // Thickness [microns], thickness variation [%], reference energy [eV/u].
             Stripper_Para[]    = {3.0, 20.0, 16.623e6},
             // E0 after stripper [eV/u], energy dependance, gap dependance. Unit for last parameter ???
             Stripper_E0Para[]  = {16.348e6, 1.00547, -0.10681},
             // E1 after stripper; Gaussian distribution.
             Stripper_E1Para    = 2.8874e-3,
             // Theta after stripper.
             Stripper_lambda    = 5.5740,
             // U after stripper.
             Stripper_upara     = 2.6903;


struct ElementStripper : public Moment2ElementBase
{
    // Transport (identity) matrix for a Charge Stripper.
    typedef ElementStripper          self_t;
    typedef Moment2ElementBase       base_t;
    typedef typename base_t::state_t state_t;
    ElementStripper(const Config& c)
        :base_t(c)
    {
        // Identity matrix.
        length = 0e0;

        const std::vector<double>& ChgState = c.get<std::vector<double> >("IonChargeStates");
        const std::vector<double>& NChg = c.get<std::vector<double> >("NCharge");
        if(ChgState.size()!=NChg.size())
            throw std::runtime_error("charge stripper requires that IonChargeStates[] and NCharge[] have the same length");
    }
    virtual ~ElementStripper() {}

    virtual void assign(const ElementVoid *other) { base_t::assign(other); }

    virtual void advance(StateBase &s);

    virtual const char* type_name() const {return "stripper";}
};


void Stripper_GetMat(const Config &conf,
                     Moment2State &ST);


#endif // CHG_STRIPPER_H
