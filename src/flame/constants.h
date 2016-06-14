#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <math.h>

// Phase-space units.
#if false
    // Use [m, rad, m, rad, rad, eV/u].
    #define MtoMM   1e0
    #define MeVtoeV 1e0
#else
    // Use [mm, rad, mm, rad, rad, MeV/u].
    #define MtoMM   1e3
    #define MeVtoeV 1e6
#endif

// Physical constants

#ifndef M_PI
// normally from math.h
# define M_PI		3.14159265358979323846	/* pi */
#endif

// Speed of light [m/s].
# define C0           2.99792458e8
// Atomic mass unit [eV/c^2].
# define AU           (931.49432e6/MeVtoeV)
// Vacuum permeability.
# define MU0          4e0*M_PI*1e-7


#endif // CONSTANTS_H
