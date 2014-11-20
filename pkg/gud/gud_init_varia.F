C $Header$
C $Name$

#include "GUD_OPTIONS.h"
#include "EXF_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: GUD_INIT_VARIA

C     !INTERFACE:
      SUBROUTINE GUD_INIT_VARIA( myThid )

C     !DESCRIPTION:
C     Initialize GUD constants

C     !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#include "PARAMS.h"
#include "EXF_PARAM.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#ifdef ALLOW_DIAZ
#include "PTRACERS_PARAMS.h"
#endif
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"
#include "GUD_GENPARAMS.h"
#include "GUD_TRAITS.h"
#include "GUD_EXF.h"
#include "GUD_FIELDS.h"

C     !INPUT PARAMETERS:
      INTEGER myThid
CEOP

#ifdef ALLOW_GUD
C     !LOCAL VARIABLES:
      _RL dTsub(Nr)
      INTEGER bi,bj,jp

      CALL GUD_EXF_INIT_VARIA(myThid)

#ifdef ALLOW_CARBON
      DO bj=myByLo(mythid),myByHi(myThid)
       DO bi=myBxLo(mythid),myBxHi(myThid)
        AtmosP(:,:,bi,bj) = 1.0 _d 0
        IF (pCO2File .EQ. ' ') atmospCO2(:,:,bi,bj) = gud_atmos_pCO2
        IF (windFile .EQ. ' ') windSpeed(:,:,bi,bj) = 5.0 _d 0
       ENDDO
      ENDDO
#endif

#ifndef CHLQUOTA
      chlPrev = 0.0 _d 0
#endif

#ifdef ALLOW_DIAZ
      IF (nIter0.EQ.PTRACERS_Iter0) THEN
       DO jp=1,nplank
        IF (diazotroph(jp) .EQ. 1.0 _d 0) THEN
         Ptracer(:,:,:,:,:,ic+jp-1) = Ptracer(:,:,:,:,:,ic+jp-1) *
     &                                diaz_ini_fac
#ifdef NQUOTA
         Ptracer(:,:,:,:,:,in+jp-1) = Ptracer(:,:,:,:,:,in+jp-1) *
     &                                diaz_ini_fac
#endif
#ifdef PQUOTA
         Ptracer(:,:,:,:,:,ip+jp-1) = Ptracer(:,:,:,:,:,ip+jp-1) *
     &                                diaz_ini_fac
#endif
#ifdef FEQUOTA
         Ptracer(:,:,:,:,:,ife+jp-1) = Ptracer(:,:,:,:,:,ife+jp-1) *
     &                                 diaz_ini_fac
#endif
#ifdef SIQUOTA
         Ptracer(:,:,:,:,:,isi+jp-1) = Ptracer(:,:,:,:,:,isi+jp-1) *
     &                                 diaz_ini_fac
#endif
#ifdef CHLQUOTA
         Ptracer(:,:,:,:,:,ichl+jp-1) = Ptracer(:,:,:,:,:,ichl+jp-1) *
     &                                  diaz_ini_fac
#endif
        ENDIF
       ENDDO
      ENDIF
#endif

      IF(gud_read_phos)THEN
       Ptracer(:,:,:,:,:,ic:ec) = 120.0 _d 0*Ptracer(:,:,:,:,:,ic:ec)
      ENDIF

      IF (nIter0 .EQ. gud_chlIter0) THEN
C      Initialize Chl from balanced-growth Chl:C if requested
C      and check Chl:C bounds.
C      use current time for insol
       dTsub(:) = 0.0 _d 0
       DO bj=myByLo(mythid),myByHi(myThid)
       DO bi=myBxLo(mythid),myBxHi(myThid)
        CALL GUD_INIT_CHL(dTsub, bi, bj, startTime, nIter0, myThid)
       ENDDO
       ENDDO
       DO jp = 1,nChl
       _EXCH_XYZ_RL(pTracer(1-OLx,1-OLy,1,1,1,iChl+jp-1),myThid)
       ENDDO
      ENDIF

#ifdef ALLOW_CARBON
      CALL GUD_SURFFORCING_INIT( myThid )
#endif

#ifdef ALLOW_GUD_CONS
      IF ( myProcId.EQ.0 .AND. myThid.EQ.1 ) THEN
       CALL MDSFINDUNIT( GUD_cons_c_unit, mythid )
       open(GUD_cons_C_unit,file='gud_cons_C.txt',status='unknown')
       CALL MDSFINDUNIT( GUD_cons_p_unit, mythid )
       open(GUD_cons_P_unit,file='gud_cons_P.txt',status='unknown')
       CALL MDSFINDUNIT( GUD_cons_n_unit, mythid )
       open(GUD_cons_N_unit,file='gud_cons_N.txt',status='unknown')
       CALL MDSFINDUNIT( GUD_cons_fe_unit, mythid )
       open(GUD_cons_Fe_unit,file='gud_cons_Fe.txt',status='unknown')
       CALL MDSFINDUNIT( GUD_cons_si_unit, mythid )
       open(GUD_cons_Si_unit,file='gud_cons_Si.txt',status='unknown')
#ifdef ALLOW_CARBON
       CALL MDSFINDUNIT( GUD_cons_a_unit, mythid )
       open(GUD_cons_A_unit,file='gud_cons_A.txt',status='unknown')
       CALL MDSFINDUNIT( GUD_cons_o_unit, mythid )
       open(GUD_cons_O_unit,file='gud_cons_O.txt',status='unknown')
#endif
c     myProcId and myThid
      ENDIF
#endif

#endif /* ALLOW_GUD */

      RETURN
      END

