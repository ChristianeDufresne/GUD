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
#undef  ALLOW_RADTRANS
#undef  ALLOW_CDOM

C grazing choices:
C   grazing a la monod:           #define GUD_GRAZING_MONOD_STYLE
C   linear grazing a la quota:    both #undef
C   quadratic grazing a la quota: #define GUD_GRAZING_SWITCH
#define GUD_GRAZING_MONOD_STYLE
#undef  GUD_GRAZING_SWITCH

#undef  NOTEMP
#define TEMP_VERSION 1
#undef  TEMP_RANGE

#undef  NINE_SPECIES_SETUP
#undef  TWO_SPECIES_SETUP

#undef  NQUOTA
#undef  PQUOTA
#undef  FEQUOTA
#undef  SIQUOTA
#undef  CHLQUOTA

#define MINFE
#undef  PART_SCAV
#undef  IRON_SED_SOURCE_VARIABLE

#undef  UNCERTAINTY

#undef  ALLOW_DENIT
#undef  ALLOW_CDOM
#undef  ALLOW_DIAZ
#undef  ALLOW_EXUDE

#undef  ALLOW_CARBON

C deprecated
#define PORT_RAND
#undef  OLDSEED

#define ALLOW_GUD_CONS

#undef GUD_DEBUG
#define iDEBUG 198
#define jDEBUG 13

#endif /* ALLOW_GUD */
#endif /* GUD_OPTIONS_H */
