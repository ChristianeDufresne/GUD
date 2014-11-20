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
#undef  ALLOW_DIAZ
#undef  ALLOW_MUTANTS
#define PORT_RAND

#define GUD_GRAZING_MONOD_STYLE

#undef NOTEMP
#define TEMP_VERSION 1
#define TEMP_RANGE

#undef TWO_SPECIES_SETUP
#undef NINE_SPECIES_SETUP

#define ALLOW_RADTRANS
#define GEIDER
#define DAR_CALC_ACDOM
#define ALLOW_CDOM

#undef  NQUOTA
#undef  PQUOTA
#undef  FEQUOTA
#undef  SIQUOTA
#undef  CHLQUOTA

#undef  RELAX_NUTS
#undef  FLUX_NUTS

#undef  CHECK_CONS
#undef  DAR_DIAG_RSTAR
#undef  DAR_DIAG_DIVER
#undef  DAR_DIAG_GROW
#define DAR_DIAG_ACDOM
#define DAR_DIAG_ABSORP
#define DAR_DIAG_SCATTER
#define DAR_DIAG_PART_SCATTER
#define DAR_DIAG_IRR
#define DAR_DIAG_IRR_AMPS

C diagnostic chlorophyll
#undef  DAR_DIAG_CHL

#define DAR_CHECK_IRR_CONT

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
