C $Header$
C $Name$
C TODO
C - have to call fechem after model?

#include "GUD_OPTIONS.h"

CBOP
C !ROUTINE: GUD_FORCING
C !INTERFACE: ==========================================================
      SUBROUTINE GUD_FORCING( dTsub, iMin, iMax, jMin, jMax, bi, bj,
     &                  midTime, myTime, myIter, myThid )

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#ifdef ALLOW_GUD
#include "GCHEM.h"
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"
#include "GUD_GENPARAMS.h"
#include "GUD_TRAITS.h"
#include "GUD_FIELDS.h"
#include "GUD_EXF.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      _RL dTsub(Nr)
      _RL midTime
      _RL myTime
      INTEGER iMin, iMax, jMin, jMax, bi, bj, myIter, myThid
CEOP

#ifdef ALLOW_GUD

C!LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  k                    :: vertical level
      INTEGER i,j,k,kdn,iTr
      CHARACTER*8 diagname
      _RL PAR(sNx,sNy,Nr,nlam)
      _RL diags(sNx,sNy,Nr,gud_nDiag)
      _RL flxCO2(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL flxALK(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL flxO2(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL atten
      _RL freeFe(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr)
      _RL scv,scav_poc
      _RL flx
      _RL ptr(nGud), gtr(nGud), PARl(nlam)
      _RL chlout(nPhoto)
      _RL diagsl(gud_nDiag)
      _RL photoTempFunc(nplank)
      _RL grazTempFunc(nplank)
      _RL reminTempFunc
      _RL mortTempFunc
      _RL mort2TempFunc
      _RL uptakeTempFunc
#ifdef ALLOW_DIAGNOSTICS
      INTEGER l
#endif

C === reset tendencies =================================================
      DO itr=1,nGud
      DO k=1,Nr
      DO j=1,sNy
      DO i=1,sNx
        gPtr(i,j,k,bi,bj,iTr) = 0.0 _d 0
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C === light ============================================================
C     Initialize Chl from balanced-growth Chl:C if requested
C     and check Chl:C bounds.
C     Note: myIter has already been incremented
c      IF (myIter-1 .EQ. gud_chlIter0) THEN
c        CALL GUD_INIT_CHL(bi, bj, myTime, myIter, myThid)
c      ENDIF

      CALL TIMER_START('GUD_LIGHT [GUD_FORCING]',myThid)
#ifdef ALLOW_RADTRANS
      CALL GUD_LIGHT_RADTRANS(PAR, bi, bj, myTime, myIter, myThid)
#else
      CALL GUD_LIGHT(PAR, midTime, bi, bj, myTime, myIter, myThid)
#endif
      CALL TIMER_STOP ('GUD_LIGHT [GUD_FORCING]',myThid)

C === dic ==============================================================
#ifdef ALLOW_CARBON
C carbon air-sea interaction
       CALL GUD_SURFFORCING(
     &                    bi,bj,imin,imax,jmin,jmax,
     &                    myIter,myTime,myThid)
#endif

C === iron =============================================================
      CALL TIMER_START('GUD_FE_CHEM [GUD_FORCING]',myThid)
      CALL GUD_FE_CHEM(
     U                 Ptracer(1-OLx,1-OLy,1,bi,bj,iFeT),
     O                 freeFe(1-OLx, 1-OLy, 1),
     I                 bi, bj, myThid)
      CALL TIMER_STOP ('GUD_FE_CHEM [GUD_FORCING]',myThid)

C     iron dust input
      DO j=1,sNy
      DO i=1,sNx
      IF (hFacC(i,j,1,bi,bj) .NE. 0.) THEN
        gPtr(i,j,1,bi,bj,iFeT) =
     &  gPtr(i,j,1,bi,bj,iFeT) + alpfe * inputFe(i,j,bi,bj) /
     &                           drF(1) / hFacC(i,j,1,bi,bj)
      ENDIF
      ENDDO
      ENDDO

C     scavenging
      DO k=1,Nr
      DO j=1,sNy
      DO i=1,sNx
#ifdef PART_SCAV
        scav_poc = Ptracer(i,j,k,bi,bj,iPOP)/scav_R_POPPOC
        scv = scav_rat*scav_inter*(scav_poc**scav_exp)
#else
        scv = scav
#endif
        gPtr(i,j,k,bi,bj,iFeT) = gPtr(i,j,k,bi,bj,iFeT) -
     &                           scv*freefe(i,j,k)
      ENDDO
      ENDDO
      ENDDO

C     iron sediment source (in bottom grid cell above kMaxFeSed)
      DO k = kMinFeSed, kMaxFeSed
       kdn = MIN(Nr, k+1)
       DO j=1,sNy
        DO i=1,sNx
         IF (hFacC(i,j,k,bi,bj) .GT. 0. .AND.
     &       (hFacC(i,j,kdn,bi,bj) .EQ. 0. .OR. k .EQ. Nr)) THEN
#ifdef IRON_SED_SOURCE_VARIABLE
          flx = fesedflux_pcm*wp_sink*R_CP_fesed*
     &                        Ptracer(i,j,k-1,bi,bj,iPOP)
#else
          flx = fesedflux
#endif
          gPtr(i,j,k,bi,bj,iFeT) = gPtr(i,j,k,bi,bj,iFeT) +
     &                             flx/(drF(k)*hFacC(i,j,k,bi,bj))
         ENDIF
        ENDDO
       ENDDO
      ENDDO

C === plankton =========================================================
      CALL TIMER_START('GUD_MODEL [GUD_FORCING]',myThid)
      do k=1,Nr
      do j=1,sNy
      do i=1,sNx
      if (hFacC(i,j,k,bi,bj) .GT. 0.) then
        CALL GUD_TEMPFUNC(Theta(i,j,k,bi,bj),
     &         photoTempFunc, grazTempFunc, reminTempFunc,
     &         mortTempFunc, mort2TempFunc,
     &         uptakeTempFunc, myThid)
        ptr(:) = Ptracer(i, j, k, bi, bj, 1:nGud)
        PARl(:) = PAR(i, j, k, :)
        gtr(:) = gPtr(i, j, k, bi, bj, 1:nGud)
        CALL GUD_MODEL(ptr,
     U                    gtr,
     O                    chlout, diagsl,
     I                    PARl, photoTempFunc, reminTempFunc,
     I                    uptakeTempFunc,
#ifdef GUD_DEBUG
     I                    myxgloballo+(bi-1)*sNx+i-1,
     I                    myygloballo+(bj-1)*sNy+j-1,
     I                    k,
#endif
     I                    myTime, myIter, myThid )
        diags(i, j, k, :) = diagsl
#ifndef CHLQUOTA
#ifdef ALLOW_RADTRANS
        chlPrev(i, j, k, bi, bj, :) = chlout
#else
        chlPrev(i, j, k, bi, bj) = SUM(chlout)
#endif
#endif
#ifdef GUD_TIME_GRAZING
        gPtr(i, j, k, bi, bj, 1:nGud) = gtr(:)
      endif
      enddo
      enddo
      enddo
      CALL TIMER_STOP('GUD_MODEL [GUD_FORCING]',myThid)
      CALL TIMER_START('GUD_GRAZING [GUD_FORCING]',myThid)
      do k=1,Nr
      do j=1,sNy
      do i=1,sNx
      if (hFacC(i,j,k,bi,bj) .GT. 0.) then
        ptr(:) = Ptracer(i, j, k, bi, bj, 1:nGud)
        gtr(:) = gPtr(i, j, k, bi, bj, 1:nGud)
        CALL GUD_TEMPFUNC(Theta(i,j,k,bi,bj),
     &         photoTempFunc, grazTempFunc, reminTempFunc,
     &         mortTempFunc, mort2TempFunc,
     &         uptakeTempFunc, myThid)
#endif
#ifdef GUD_GRAZING_MONOD_STYLE
        CALL GUD_GRAZING_MONOD(ptr,
     U                    gtr,
     I                    grazTempFunc, reminTempFunc,
     I                    mortTempFunc, mort2TempFunc,
     I                    myTime, myIter, myThid )
#else
        CALL GUD_GRAZING(ptr,
     U                    gtr,
     I                    grazTempFunc, reminTempFunc,
     I                    mortTempFunc, mort2TempFunc,
     I                    myTime, myIter, myThid )
#endif
        gPtr(i, j, k, bi, bj, 1:nGud) = gtr(:)
      endif
      enddo
      enddo
      enddo
#ifdef GUD_TIME_GRAZING
      CALL TIMER_STOP('GUD_GRAZING [GUD_FORCING]',myThid)
#else
      CALL TIMER_STOP('GUD_MODEL [GUD_FORCING]',myThid)
#endif

C === sinking ==========================================================
      CALL TIMER_START('GUD_SINKING [GUD_FORCING]',myThid)
      CALL GUD_SINKING( Ptracer,gPtr,bi,bj,myTime,myIter,myThid )
      CALL TIMER_STOP ('GUD_SINKING [GUD_FORCING]',myThid)

C === apply tendencies to tracers ======================================
      DO iTr=1,nGud
      DO k=1,Nr
      DO j=1,sNy
      DO i=1,sNx
        pTracer(i,j,k,bi,bj,iTr)=pTracer(i,j,k,bi,bj,iTr)
     &                          +dTsub(k)*gPtr(i,j,k,bi,bj,iTr)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C === iron =============================================================
C     re-apply free iron limit to FeT
      CALL TIMER_START('GUD_FE_CHEM [GUD_FORCING]',myThid)
      CALL GUD_FE_CHEM(
     U                 pTracer(1-OLx,1-OLy,1,bi,bj,iFeT),
     O                 freeFe(1-OLx,1-OLy,1),
     I                 bi, bj, myThid)
      CALL TIMER_STOP ('GUD_FE_CHEM [GUD_FORCING]',myThid)

C === diagnostics ======================================================
#ifdef ALLOW_DIAGNOSTICS
      IF (useDIAGNOSTICS) THEN
       CALL TIMER_START('DIAGS_FILL [GUD_FORCING]',myThid)
       DO l = 1, nlam
        WRITE(diagname, '(A,I2.2)') 'PAR', l
        CALL DIAGNOSTICS_FILL(PAR(1,1,1,l),diagname,0,Nr,3,bi,bj,myThid)
       ENDDO
       CALL DIAGNOSTICS_FILL(diags(1,1,1,iPP),   'PP      ',0,Nr,3,
     &          bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(diags(1,1,1,iNfix), 'Nfix    ',0,Nr,3,
     &          bi,bj,myThid)
       CALL DIAGNOSTICS_FILL(diags(1,1,1,iDenit),'Denit   ',0,Nr,3,
     &          bi,bj,myThid)
       CALL TIMER_STOP ('DIAGS_FILL [GUD_FORCING]',myThid)
      ENDIF
#endif

#endif /* ALLOW_GUD */

      RETURN
      END
