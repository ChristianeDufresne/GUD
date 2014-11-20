C $Header$
C $Name$

#include "GUD_OPTIONS.h"

CBOP
C !ROUTINE: GUD_SURFFORCING_INIT

C !INTERFACE: ==========================================================
      SUBROUTINE GUD_SURFFORCING_INIT(myThid)

C !DESCRIPTION:
C  Calculate first guess of pH

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#ifdef ALLOW_GUD
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"
#include "GUD_EXF.h"
#include "GUD_GENPARAMS.h"
#include "GUD_TRAITS.h"
#include "GUD_FIELDS.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER  myThid

#ifdef ALLOW_GUD
#ifdef ALLOW_CARBON

C !LOCAL VARIABLES: ====================================================
       INTEGER i,j,k,it
       _RL co3dummy
C Number of iterations for pCO2 solvers...
C Solubility relation coefficients
C local variables for carbon chem
      INTEGER iMin,iMax,jMin,jMax, bi, bj
      _RL surfdic
      _RL surfalk
      _RL surfphos
      _RL surfsi
      _RL surftemp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfsalt(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iprt,jprt
      LOGICAL pH_isLoaded
CEOP

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      k=1

#ifdef ALLOW_OFFLINE
      IF ( useOffLine ) THEN
        CALL OFFLINE_FIELDS_LOAD( startTime, nIter0, myThid )
      ENDIF
#endif

C =================================================================

      jMin=1
      jMax=sNy
      iMin=1
      iMax=sNx

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
            pH(i,j,bi,bj) = 8. _d 0
          ENDDO
         ENDDO
       ENDDO
      ENDDO

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          ak0(i,j,bi,bj)=0. _d 0
          ak1(i,j,bi,bj)=0. _d 0
          ak2(i,j,bi,bj)=0. _d 0
          akw(i,j,bi,bj)=0. _d 0
          akb(i,j,bi,bj)=0. _d 0
          akf(i,j,bi,bj)=0. _d 0
          ak1p(i,j,bi,bj)=0. _d 0
          ak2p(i,j,bi,bj)=0. _d 0
          ak3p(i,j,bi,bj)=0. _d 0
          aksi(i,j,bi,bj)=0. _d 0
          fugf(i,j,bi,bj)=0. _d 0
          ff(i,j,bi,bj)=0. _d 0
          ft(i,j,bi,bj)=0. _d 0
          st(i,j,bi,bj)=0. _d 0
          bt(i,j,bi,bj)=0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

c      pH_isLoaded = .FALSE.
c      IF ( nIter0.GT.PTRACERS_Iter0 .OR.
c     &    (nIter0.EQ.PTRACERS_Iter0 .AND. pickupSuff.NE.' ')
c     &   ) THEN
cC       Read pH from a pickup file if needed
c        CALL DIC_READ_PICKUP(
c     O                        pH_isLoaded,
c     I                        nIter0, myThid )
c      ENDIF

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

C determine inorganic carbon chem coefficients
        DO j=jMin,jMax
         DO i=iMin,iMax
          surfsalt(i,j) = 
     &            max(4. _d 0, min(50. _d 0, salt(i,j,k,bi,bj)))
          surftemp(i,j) =
     &            max(-4. _d 0, min(39. _d 0, theta(i,j,k,bi,bj)))
         ENDDO
        ENDDO

        CALL GUD_CARBON_COEFFS(
     I                       surftemp,surfsalt,
     I                       bi,bj,iMin,iMax,jMin,jMax,myThid)

C====================================================================

        IF ( .NOT.pH_isLoaded ) THEN
C set guess of pH for first step here

          WRITE(standardMessageUnit,*) 'QQ: pCO2 approximation method'
          iprt = MIN(20,sNx)
          jprt = MIN(20,sNy)
c first approximation
C$TAF LOOP = parallel
          DO j=jMin,jMax
C$TAF LOOP = parallel
           DO i=iMin,iMax
            IF ( maskC(i,j,k,bi,bj) .NE. 0. _d 0) THEN
c put bounds on tracers so pH solver doesn't blow up
             surfdic  = MAX(0.01 _d 0, MIN(4. _d 0,
     &                            ptr2mol * Ptracer(i,j,k,bi,bj,iDIC)))
     &                  * maskC(i,j,k,bi,bj)
             surfalk  = MAX(0.01 _d 0, MIN(4. _d 0,
     &                            ptr2mol * Ptracer(i,j,k,bi,bj,iALK)))
     &                  * maskC(i,j,k,bi,bj)
             surfphos = MAX(1 _d -13, MIN(1 _d -2,
     &                            ptr2mol * Ptracer(i,j,k,bi,bj,iPO4)))
     &                  * maskC(i,j,k,bi,bj)
             surfsi   = MAX(1. _d -11, MIN(0.5 _d 0,
     &                            ptr2mol * Ptracer(i,j,k,bi,bj,iSiO2)))
     &                  * maskC(i,j,k,bi,bj)
C$TAF init dic_surf = static, 10
             DO it=1,10
C$TAF STORE pH(i,j,bi,bj)                              = dic_surf
C$TAF STORE surfalk(i,j), surfphos(i,j), surfsi(i,j)   = dic_surf
              CALL GUD_CALC_PCO2_APPROX(
     I        surftemp(i,j), surfsalt(i,j),
     I        surfdic, surfphos, surfsi, surfalk,
     I        ak1(i,j,bi,bj), ak2(i,j,bi,bj),
     I        ak1p(i,j,bi,bj), ak2p(i,j,bi,bj), ak3p(i,j,bi,bj),
     I        aks(i,j,bi,bj), akb(i,j,bi,bj), akw(i,j,bi,bj),
     I        aksi(i,j,bi,bj), akf(i,j,bi,bj),
     I        ak0(i,j,bi,bj), fugf(i,j,bi,bj),
     I        ff(i,j,bi,bj),
     I        bt(i,j,bi,bj), st(i,j,bi,bj), ft(i,j,bi,bj),
     U        pH(i,j,bi,bj), pCO2(i,j,bi,bj), co3dummy,
     I        i,j,k,bi,bj, it ,  myThid )
             ENDDO
            ENDIF
            IF (i.EQ.iprt .AND. j.EQ.jprt)
     &       WRITE(standardMessageUnit,*) 'QQ first guess pH',
     &        pH(iprt,jprt,bi,bj),
     &        theta(iprt,jprt,1,bi,bj), salt(iprt,jprt,1,bi,bj),
     &        surfdic, surfphos,surfsi,surfalk
           ENDDO
          ENDDO

        ENDIF

C     end bi,bj loops
       ENDDO
      ENDDO

#endif /* ALLOW_CARBON */
#endif /* ALLOW_GUD */

      RETURN
      END
