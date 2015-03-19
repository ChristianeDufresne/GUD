C $Header$
C $Name$

#include "GUD_OPTIONS.h"

CBOP
C !ROUTINE: GUD_LIGHT_RADTRANS
C !INTERFACE: ==========================================================
      SUBROUTINE GUD_LIGHT_RADTRANS(
     O                         PAR,
     I                         solTime,
     I                         bi, bj,
     I                         myTime, myIter, myThid )

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
C#include "DYNVARS.h"
#include "FFIELDS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#ifdef ALLOW_GUD
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"
#include "GUD_EXF.h"
#include "GUD_PARAMS.h"
#include "GUD_GENPARAMS.h"
#include "GUD_TRAITS.h"
#include "GUD_FIELDS.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myTime :: time at end of (sub)timestep
C  myThid :: thread number
      _RL solTime
      _RL myTime
      INTEGER bi, bj, myIter, myThid

C !OUTPUT PARAMETERS: ==================================================
      _RL PAR(sNx, sNy, Nr, nlam)
CEOP

#ifdef ALLOW_GUD
#ifdef GUD_ALLOW_RADTRANS

C!LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  k                    :: vertical level
      INTEGER i,j,k,l,jp
      INTEGER iyr,imon,iday,isec,lp,wd,mydate(4)
      _RL solz(sNx,sNy)
      _RL rmud
      _RL part
      _RL aCDOM(nlam)
      _RL phychl(nPhoto)
#ifdef GUD_ALLOW_CDOM
      _RL CDOM
#else
      _RL actotref, atotref
#endif
      _RL PAR_k(Nr, nlam)
      _RL dz_k(Nr)
      _RL actot, bctot, bbctot
C
      _RL a_k(Nr,nlam)
      _RL bt_k(Nr,nlam)
      _RL bb_k(Nr,nlam)
      _RL Edwsf(nlam),Eswsf(nlam)
      _RL Edbot(nlam,Nr),Esbot(nlam,Nr),Eubot(nlam,Nr)
      _RL Estop(nlam,Nr),Eutop(nlam,Nr)
      _RL amp1(nlam,Nr), amp2(nlam,Nr)

      _RL solzlat(1-OLy:sNy+OLy)

C ======================================================================
C get current date and time of day: iyr/imon/iday+isec
      CALL CAL_GETDATE( myIter, mytime, mydate, mythid )
      CALL CAL_CONVDATE( mydate,iyr,imon,iday,isec,lp,wd,mythid )

      IF (gud_useModelCalendar) THEN
       CALL GUD_INSOL(solTime, solzlat, bj, .TRUE.)
      ENDIF

      DO j=1,sNy
       DO i=1,sNx
        IF (gud_useModelCalendar) THEN
         solz(i,j) = solzlat(j)
        ELSE
C oj: average light effective at noon?
         isec = 12*3600
         CALL SUN_SFCSOLZ(rad2deg, iyr, imon, iday, isec,
     I                    XC(i,j,bi,bj), YC(i,j,bi,bj),
     O                    solz(i,j))
        ENDIF

c compute 1/cos(zenith) for direct light below surface
        CALL SUN_SFCRMUD(rad2deg, solz(i,j),
     O                   rmud)

C ======================================================================
        DO k=1,Nr

         dz_k(k) = drF(k)*HFacC(i,j,k,bi,bj)
         part = MAX(Ptracer(i,j,k,bi,bj,iPOP), 0. _d 0)
         DO jp=1,nPhoto
#ifdef GUD_ALLOW_CHLQUOTA
          phychl(jp)=MAX(Ptracer(i,j,k,bi,bj,iChl+jp-1),0.)
#else
          phychl(jp)=MAX(chlPrev(i,j,k,bi,bj,jp), 0.)
#endif
         ENDDO

#ifdef GUD_ALLOW_CDOM
c use cdom-like tracer
         CDOM = MAX(Ptracer(i,j,k,bi,bj,iCDOM), 0.0 _d 0)
         DO l = 1,nlam
          aCDOM(l) = CDOMcoeff*CDOM*exCDOM(l)
         ENDDO
#else
         actotref = 0.0 _d 0
         atotref = 0.0 _d 0
         DO jp = 1,nPhoto
c         nb. n,k swapped from WG
          actotref = actotref +
     &          phychl(jp)*aphy_chl(iMinPhoto+jp-1,laCDOM)
         ENDDO
         atotref = aw(laCDOM) + actotref
         DO l = 1,nlam
          aCDOM(l) = gud_aCDOM_fac*atotref*exCDOM(l)
         ENDDO
#endif

         DO l = 1,nlam
c         absorption by phyto
          actot = 0.0
          bctot = 0.0
          bbctot = 0.0
          DO jp = 1, nPhoto
           actot  = actot  + phychl(jp)*aphy_chl(iMinPhoto+jp-1,l)
           bctot  = bctot  + phychl(jp)*bphy_chl(iMinPhoto+jp-1,l)
           bbctot = bbctot + phychl(jp)*bbphy_chl(iMinPhoto+jp-1,l)
          ENDDO
c         add water and CDOM
          a_k(k,l)  = aw(l) + aCDOM(l) + actot  + part*apart_P(l)
          bt_k(k,l) = bw(l)            + bctot  + part*bpart_P(l)
          bb_k(k,l) = gud_bbw*bw(l) + bbctot + part*bbpart_P(l)
          bb_k(k,l) = MAX(gud_bbmin, bb_k(k,l))
         ENDDO

C       k
        ENDDO
C ======================================================================

        Edwsf(:) = OASIM_Ed(i,j,bi,bj,:)
        Eswsf(:) = OASIM_Es(i,j,bi,bj,:)

        CALL GUD_RADTRANS_DIRECT(
     I                dz_k,rmud,Edwsf,Eswsf,a_k,bt_k,bb_k,
     I                gud_radtrans_kmax,
     O                Edbot,Esbot,Eubot,Estop,Eutop,
     O                PAR_k,
     O                amp1,amp2,
     I                myThid)

        DO l = 1,nlam
         DO k = 1,Nr
          PAR(i,j,k,l) = PAR_k(k,l)
         ENDDO
        ENDDO

C      i,j
       ENDDO
      ENDDO
C ======================================================================

#endif /* GUD_ALLOW_RADTRANS */
#endif /* ALLOW_GUD */

      RETURN
      END
