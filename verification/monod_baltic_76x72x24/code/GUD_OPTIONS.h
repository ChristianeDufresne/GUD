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

#undef  READ_PAR
#define USE_QSW
#define MINFE
#undef  NUT_SUPPLY
#undef  CONS_SUPP
#undef  ALLOW_DIAZ
#define PORT_RAND
#undef  OLDSEED

#define GUD_GRAZING_MONOD_STYLE

#undef NOTEMP
#define TEMP_VERSION 1
#define TEMP_RANGE

#undef TWO_SPECIES_SETUP
#undef NINE_SPECIES_SETUP

#define ALLOW_GUD_CONS
#define DAR_DIAG_RSTAR
#define DAR_DIAG_DIVER

#define GUD_VECTORIZE -1

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
