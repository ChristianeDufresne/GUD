C $Header$
C $Name$

C     *==========================================================*
C     | SEAICE_OPTIONS.h
C     | o CPP options file for sea ice package.
C     *==========================================================*
C     | Use this file for selecting options within the sea ice
C     | package.
C     *==========================================================*

#ifndef SEAICE_OPTIONS_H
#define SEAICE_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_SEAICE
C     Package-specific Options & Macros go here

C--   Write "text-plots" of certain fields in STDOUT for debugging.
#undef SEAICE_DEBUG

C--   Allow sea-ice dynamic code.
C     This option is provided to allow use of TAMC
C     on the thermodynamics component of the code only.
C     Sea-ice dynamics can also be turned off at runtime
C     using variable SEAICEuseDYNAMICS.
#undef SEAICE_ALLOW_DYNAMICS

C--   By default, the sea-ice package uses its own integrated bulk
C     formulae to compute fluxes (fu, fv, EmPmR, Qnet, and Qsw) over
C     open-ocean.  When this flag is set, these variables are computed
C     in a separate external package, for example, pkg/exf, and then
C     modified for sea-ice effects by pkg/seaice.
#define SEAICE_EXTERNAL_FLUXES

C--   By default, the sea-ice package uses 2-category thermodynamics.
C     When this flag is set, an 8-category calculation of ice
C     thermodynamics is carried out in groatb.F
C     Note the pickup_seaice.* generated by this option differ
C     from those generated with the default 2-category model.
C     Therefore it is not possible to switch between the two
C     in the middle of an integration.
#undef SEAICE_MULTICATEGORY
#undef SEAICE_MODIFY_GROWTH_ADJ
#undef SEAICE_SIMPLIFY_GROWTH_ADJ

#undef SEAICE_ADD_SUBLIMATION_TO_FWBUDGET

C--   By default cdm wind stress under sea-ice is set to the
C     same value as it would be if there was no sea-ice.
C     This is a placeholder until a more physically-sensible
C     stable solution is developed.  Apparently ocean stress
C     that is based on sea-ice velocity can cause trouble in
C     Arctic Ocean configurations.  This temporary solution
C     also corresponds to that adopted by CCSM-2.  For the
C     adventurous user who requires ice-modified stress,
C     following options are available.
#undef SEAICE_ORIGINAL_BAD_ICE_STRESS
#undef SEAICE_TEST_ICE_STRESS_1

C--   By default the freezing point of water is set to the value of
C     the parameter SEAICE_freeze (=-1.96 by default). To use a
C     simple linear dependence of the freezing point on salinity,
C     set the following flag (pressure is assumed to have no effect,
C     which is a good assumption for the top 20 meters). With this
C     option defined the parameter SEAICE_freeze has no effect.
#undef SEAICE_VARIABLE_FREEZING_POINT

C--   By default the seaice model is discretized on a B-Grid (for
C     historical reasons). Define the following flag to use a new
C     (not thoroughly) test version on a C-grid
#define SEAICE_CGRID

C--   Only for the C-grid version it is possible to enable EVP code by
C     defining the following flag
#ifdef SEAICE_CGRID
#define SEAICE_ALLOW_EVP
#endif /* SEAICE_CGRID */

C--   Seaice flooding
#define ALLOW_SEAICE_FLOODING

C--   By default sea ice is fresh.  Set following flag for salty ice.
#define SEAICE_VARIABLE_SALINITY

c--   Enable modified sea-ice thermodynamic by Ian Fenty
c--   Should improve adjoint, but doesnt work so far.
#undef SEAICE_ALLOW_TD_IF

#endif /* ALLOW_SEAICE */
#endif /* SEAICE_OPTIONS_H */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
