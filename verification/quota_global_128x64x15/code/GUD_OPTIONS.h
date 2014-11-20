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

#define PORT_RAND

#define READ_PAR
#undef  USE_AVPAR

#define MINFE
#undef  NUT_SUPPLY

#undef NOTEMP
#define TEMP_VERSION 3
#define TEMP_RANGE

#undef  TWO_SPECIES_SETUP
#undef  NINE_SPECIES_SETUP

#define GEIDER
#undef  ALLOW_CARBON

c ANNA turn wavebands on/off
#undef WAVEBANDS 

#define ALLOW_GUD_CONS

#undef GUD_GRAZING_MONOD_STYLE
#define NQUOTA
#define FEQUOTA
#define CHLQUOTA

#define IRON_SED_SOURCE_VARIABLE

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
 
#endif /* ALLOW_GUD */
#endif /* GUD_OPTIONS_H */
