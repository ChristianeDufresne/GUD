C $Header$
C $Name$

#ifndef GUD_OPTIONS_H
#define GUD_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#ifdef ALLOW_GUD

#include "CPP_OPTIONS.h"

CBOP
C    !ROUTINE: GUD_OPTIONS.h
C    !INTERFACE:

C    !DESCRIPTION:
c options for gud package
CEOP

#define READ_PAR
#define MINFE
#undef  NUT_SUPPLY
#undef  CONS_SUPP
#undef  OLD_GRAZE
#undef  ALLOW_DIAZ
#undef  OLD_NSCHEME
#undef  ALLOW_MUTANTS
#define PORT_RAND
#undef OLDSEED

#undef NOTEMP
#define TEMP_VERSION 1
#define TEMP_RANGE

#undef  TWO_SPECIES_SETUP
#undef  NINE_SPECIES_SETUP

#define GEIDER
#define DYNAMIC_CHL
#undef  ALLOW_CARBON

#define CHLQUOTA

c ANNA turn wavebands on/off
#undef WAVEBANDS 
c#define WAVEBANDS

#define ALLOW_GUD_CONS
#undef  DAR_DIAG_RSTAR
#undef  DAR_DIAG_DIVER

#undef  USE_AVPAR

#define GUD_GRAZING_MONOD_STYLE

#define GUD_UNUSED 0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C dependencies
c if two or nine species setup we don't want specific temperature ranges
#ifdef  TWO_SPECIES_SETUP
#undef TEMP_RANGE
#endif
#ifdef  NINE_SPECIES_SETUP
#undef TEMP_RANGE
#endif

#ifdef DAR_DIAG_CHL
#define ALLOW_PAR_DAY
#endif
 
#endif /* ALLOW_GUD */
#endif /* GUD_OPTIONS_H */
