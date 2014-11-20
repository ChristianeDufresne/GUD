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
C options for gud package
CEOP

#define READ_PAR
#undef  USE_QSW
#undef  USE_AVPAR
#undef  GEIDER

#undef  WAVEBANDS

#define GUD_GRAZING_MONOD_STYLE
#undef  SER_GRAZ

#undef  NOTEMP
#define TEMP_VERSION 2
#define TEMP_RANGE
#define NOZOOTEMP

C carbon-related
#define ALLOW_CARBON
#undef DIC_NO_NEG
#define ALLOW_OLD_VIRTUALFLUX

#undef  TWO_SPECIES_SETUP
#define NINE_SPECIES_SETUP

#undef  NQUOTA
#undef  PQUOTA
#undef  FEQUOTA
#undef  SIQUOTA
#undef  CHLQUOTA

#define MINFE
#define PART_SCAV
#define IRON_SED_SOURCE_VARIABLE

#undef  UNCERTAINTY

#define ALLOW_DENIT
#undef  ALLOW_CDOM
#define ALLOW_DIAZ
#undef  ALLOW_EXUDE

C deprecated
#define PORT_RAND
#undef  OLDSEED

#define ALLOW_GUD_CONS

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

