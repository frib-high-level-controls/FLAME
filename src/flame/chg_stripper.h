#ifndef CHG_STRIPPER_H
#define CHG_STRIPPER_H

#include <boost/numeric/ublas/matrix.hpp>

#include "moment.h"

typedef MomentState state_t;


// Phase space dimension; including vector for orbit/1st moment.
# define PS_Dim MomentState::maxsize // Set to 7; to include orbit.


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
const double Stripper_IonZ_default      = 78.0/238.0,
             Stripper_IonMass_default   = 238.0,
             Stripper_IonProton_default = 92.0,
             // E1 after stripper; Gaussian distribution.
             Stripper_E1Para_default    = 2.8874e-3,
             // Theta after stripper.
             Stripper_lambda_default    = 5.5740,
             // U after stripper.
             Stripper_upara_default     = 2.6903,
			 // Thickness [microns], thickness variation [%], reference energy [eV/u].
             Stripper_Para_default[] = {3.0, 20.0, 16.623e6},
			 // E0 after stripper [eV/u], energy dependance, gap dependance. Unit for last parameter ???
             Stripper_E0Para_default[] = {16.348e6, 1.00547, -0.10681};


struct ElementStripper : public MomentElementBase
{
    // Transport (identity) matrix for a Charge Stripper.
    typedef ElementStripper          self_t;
    typedef MomentElementBase       base_t;
    typedef typename base_t::state_t state_t;


    double  Stripper_IonZ,
			Stripper_IonMass,
			Stripper_IonProton,
			// E1 after stripper; Gaussian distribution.
			Stripper_E1Para,
			// Theta after stripper.
			Stripper_lambda,
			// U after stripper.
			Stripper_upara;

    // Thickness [microns], thickness variation [%], reference energy [eV/u].
    // E0 after stripper [eV/u], energy dependance, gap dependance. Unit for last parameter ???
    std::vector<double> Stripper_Para, Stripper_E0Para;


    ElementStripper(const Config& c)
        :base_t(c)
    {
        length = 0e0;
    }
    virtual ~ElementStripper() {}

    virtual void assign(const ElementVoid *other) { base_t::assign(other); }

    virtual void advance(StateBase &s);

    virtual const char* type_name() const {return "stripper";}

    void StripperCharge(const double beta, double &Q_ave, double &d);
    void ChargeStripper(const double beta, const std::vector<double>& ChgState, std::vector<double>& chargeAmount_Baron);
    void Stripper_Propagate_ref(Particle &ref);
    void Stripper_GetMat(const Config &conf, MomentState &ST);
};

#endif // CHG_STRIPPER_H
