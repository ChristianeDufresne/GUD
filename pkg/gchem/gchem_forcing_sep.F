C $Header$
C $Name$

#include "GCHEM_OPTIONS.h"
#ifdef ALLOW_DIC
#include "DIC_OPTIONS.h"
#endif
#ifdef ALLOW_DARWIN
#include "DARWIN_OPTIONS.h"
#endif
#ifdef ALLOW_GUD
#include "GUD_OPTIONS.h"
#endif


CBOP
C !ROUTINE: GCHEM_FORCING_SEP
C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_FORCING_SEP( myTime, myIter, myThid )

C !DESCRIPTION:
C     calls subroutine that will update passive tracers values
C     with a separate timestep. Since GCHEM_FORCING_SEP is now
C     called before DO_FIELDS_BLOCKING_EXCHANGES, the passive
C     tracer values in the halo regions are not up to date and
C     must not be used.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "GCHEM.h"
#ifdef ALLOW_DIC
#include "DIC_VARS.h"
#endif /* ALLOW_DIC */
#ifdef ALLOW_DARWIN
#include "DARWIN_FLUX.h"
#include "DARWIN_SIZE.h"
#endif
#ifdef ALLOW_GUD
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"
#include "GUD_GENPARAMS.h"
#ifdef ALLOW_CARBON
#include "GUD_FIELDS.h"
#endif
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      _RL myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_GCHEM
#ifdef GCHEM_SEPARATE_FORCING

C!LOCAL VARIABLES: ====================================================
C  bi,bj                :: tile indices
C  k                    :: vertical level
      INTEGER bi,bj,iMin,iMax,jMin,jMax,isub,iTracer
      INTEGER k
      _RL midTime
#ifdef ALLOW_DIAGNOSTICS
      INTEGER iTr
      CHARACTER*8 diagName
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
C--   fill-in tracer diagnostics before any GChem udate
       DO iTr = 1,PTRACERS_numInUse
        diagName = '        '
        WRITE(diagName,'(A5,A2)') 'GC_Tr', PTRACERS_ioLabel(iTr)
        CALL DIAGNOSTICS_FILL( pTracer(1-OLx,1-OLy,1,1,1,iTr), diagName,
     &                         0, Nr, 0, 1, 1, myThid )
       ENDDO
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

ccccccccccccccccccccccccc
c global calculations   c
ccccccccccccccccccccccccc
#ifdef ALLOW_OLD_VIRTUALFLUX
#ifdef ALLOW_DIC
c find global surface averages
       gsm_s = 0. _d 0
       gsm_dic = 0. _d 0
       gsm_alk = 0. _d 0
       CALL GCHEM_SURFMEAN(salt,gsm_s,myThid)
       CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,1), gsm_dic, myThid )
       print*,'mean surface dic', gsm_dic,gsm_s
       CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,2), gsm_alk, myThid )
#endif
#ifdef ALLOW_DARWIN
c find global surface averages
       gsm_s = 0. _d 0
       gsm_dic = 0. _d 0
       gsm_alk = 0. _d 0
       CALL GCHEM_SURFMEAN(salt,gsm_s,myThid)
       CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,iDIC), gsm_dic, myThid )
       print*,'mean surface dic', gsm_dic,gsm_s
       CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,iALK), gsm_alk, myThid )
#endif
ccccccccccccccccccccccccccccccccccccccccccc
#ifdef ALLOW_OFFLINE
       IF ( useOffLine ) THEN
         CALL OFFLINE_GET_SURFFORCING( myTime, myIter, myThid )
       ENDIF
#endif /* ALLOW_OFFLINE */
#endif /* ALLOW_OLD_VIRTUALFLUX */

#ifdef ALLOW_DARWIN
       IF ( useDARWIN ) THEN
         CALL DARWIN_CONS( myIter, myTime, myThid )
       ENDIF
#endif

ccccccccccccccccccccccccc
c chemical forcing      c
ccccccccccccccccccccccccc

#if defined(ALLOW_DIC) || defined(ALLOW_DARWIN)

C$taf loop = parallel
       DO bj=myByLo(myThid),myByHi(myThid)
C$taf loop = parallel
        DO bi=myBxLo(myThid),myBxHi(myThid)

        jMin=1
        jMax=sNy
        iMin=1
        iMax=sNx
c
ccccccccccccccccccccccccccc DIC cccccccccccccccccccccccccccccccc

#ifdef ALLOW_DIC
#ifdef ALLOW_FE
          CALL DIC_BIOTIC_FORCING( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,5),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,6),
     &                          bi,bj,imin,imax,jmin,jmax,
     &                          myIter,myTime,myThid)
#else
#ifdef ALLOW_O2
          CALL DIC_BIOTIC_FORCING( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,5),
     &                          bi,bj,imin,imax,jmin,jmax,
     &                          myIter,myTime,myThid)
#else
          CALL DIC_BIOTIC_Forcing( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          bi,bj,imin,imax,jmin,jmax,
     &                          myIter,myTime,myThid)
#endif
#endif

#endif
cccccccccccccccccccccccccc END DIC cccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccc DARWIN ccccccccccccccccccccccccccccccccccc
#ifdef ALLOW_DARWIN
         IF ( useDARWIN ) THEN
#ifdef NUT_SUPPLY
c articficial supply of nutrients
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,1),
     &                          bi,bj,imin,imax,jmin,jmax,
     &                          myIter,myTime,myThid)
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,2),
     &                          bi,bj,imin,imax,jmin,jmax,
     &                          myIter,myTime,myThid)
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,3),
     &                          bi,bj,imin,imax,jmin,jmax,
     &                          myIter,myTime,myThid)
          CALL DARWIN_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,4),
     &                          bi,bj,imin,imax,jmin,jmax,
     &                          myIter,myTime,myThid)
#endif
ccccccccccccccc
C darwin_forcing operates on bi,bj part only, but needs to get full
C array because of last (iPtr) index
          CALL DARWIN_FORCING(  pTracer(1-OLx,1-OLy,1,1,1,1),
     &                          bi,bj,imin,imax,jmin,jmax,
     &                          myIter,myTime,myThid)
         ENDIF
#endif /* ALLOW_DARWIN */
cccccccccccccccccccccccccc END DARWIN ccccccccccccccccccccccccccccccc

C--   end bi,bj loops.
        ENDDO
       ENDDO

#endif /* DIC or DARWIN */

#ifdef ALLOW_DARWIN
       IF ( useDARWIN ) THEN
         CALL DARWIN_CONS( myIter, myTime, myThid )
#ifdef ALLOW_CARBON
         CALL DIC_ATMOS( 1, myTime, myIter, myThid )
#endif
       ENDIF
#endif /* ALLOW_DARWIN */

#ifdef ALLOW_DIC
       CALL DIC_ATMOS( 1, myTime, myIter, myThid )
       CALL DIC_STORE_FLUXCO2( myTime, myIter, myThid )
#endif

#ifdef ALLOW_DIC
#ifdef ALLOW_COST
       CALL DIC_COST( myTime, myIter, myThid )
#endif
#endif

ccccccccccccccccccccccccccc GUD ccccccccccccccccccccccccccccc
#ifdef ALLOW_GUD
       IF ( useGUD ) THEN

         IF ( GUD_linFSConserve ) THEN
           CALL GUD_CONSERVE_SURFACE( myTime, myIter, myThid )
         ENDIF

         CALL GUD_CONS( myTime, myIter, myThid )

#ifdef ALLOW_CARBON
#ifdef ALLOW_OLD_VIRTUALFLUX
c find global surface averages
         CALL GCHEM_SURFMEAN(salt, gsm_s, myThid)
         CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,iDIC), gsm_dic, myThid )
         CALL GCHEM_SURFMEAN(
     &             pTracer(1-OLx,1-OLy,1,1,1,iALK), gsm_alk, myThid )
#endif /* ALLOW_OLD_VIRTUALFLUX */
#endif /* ALLOW_CARBON */

C$taf loop = parallel
         DO bj=myByLo(myThid),myByHi(myThid)
C$taf loop = parallel
          DO bi=myBxLo(myThid),myBxHi(myThid)

           jMin=1
           jMax=sNy
           iMin=1
           iMax=sNx

#ifdef NUT_SUPPLY
c articficial supply of nutrients
           CALL GUD_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,iPO4),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
           CALL GUD_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,iNO3),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
           CALL GUD_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,iFeT),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
           CALL GUD_NUT_SUPPLY( pTracer(1-OLx,1-OLy,1,bi,bj,iSiO2),
     &                          bi,bj,iMin,iMax,jMin,jMax,
     &                          myIter,myTime,myThid)
#endif

           CALL TIMER_START('GUD_FORCING [GCHEM_FORCING_SEP]',myThid)
C          time at middle of sub-timestep
           midTime = myTime - deltaTclock + .5*deltaTclock/nsubtime
           DO isub=1,nsubtime
            CALL GUD_FORCING( gchem_dTsub,iMin,iMax,jMin,jMax,bi,bj,
     &                        midTime,myTime,myIter,myThid)
            midTime = midTime + deltaTclock/nsubtime
           ENDDO
           CALL TIMER_STOP ('GUD_FORCING [GCHEM_FORCING_SEP]',myThid)
         
C--   end bi,bj loops.
          ENDDO
         ENDDO

         CALL GUD_CONS( myTime, myIter, myThid )

       ENDIF
#endif /* ALLOW_GUD */
cccccccccccccccccccccccccc END GUD ccccccccccccccccccccccccccccccc

#endif /* GCHEM_SEPARATE_FORCING */
#endif /* ALLOW_GCHEM */

      RETURN
      END
