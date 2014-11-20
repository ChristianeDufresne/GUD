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

#undef  GUD_DEBUG

#define READ_PAR
#undef  USE_QSW
#undef  USE_AVPAR
#define GEIDER
#undef  ALLOW_RADTRANS
#undef  ALLOW_CDOM
#define GUD_GEIDER_RHO_SYNTH

#undef  GUD_GRAZING_MONOD_STYLE

#undef  NOTEMP
#define TEMP_VERSION 2
#undef  TEMP_RANGE

#undef  NINE_SPECIES_SETUP
#undef  TWO_SPECIES_SETUP

#undef  NQUOTA
#undef  PQUOTA
#undef  FEQUOTA
#undef  SIQUOTA
#define CHLQUOTA

#define MINFE
#define PART_SCAV
#define IRON_SED_SOURCE_VARIABLE

#undef  UNCERTAINTY

#define ALLOW_DENIT
#undef  ALLOW_CDOM
#define ALLOW_DIAZ
#undef  ALLOW_EXUDE

C -> ALLOW_CARBON
#define ALLOW_CARBON

C deprecated
#define PORT_RAND
#undef  OLDSEED

#define ALLOW_GUD_CONS

#define GUD_UNUSED 0

#define ALLOW_OLD_VIRTUALFLUX


#endif /* ALLOW_GUD */
#endif /* GUD_OPTIONS_H */
