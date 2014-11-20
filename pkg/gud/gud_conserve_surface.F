C $Header$
C $Name$

#include "GUD_OPTIONS.h"

c=============================================================
c subroutine GUD_conserve_surface
c check conservation in the model
C==============================================================
      SUBROUTINE GUD_CONSERVE_SURFACE(myTime,myIter,myThid)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"
#include "GUD_GENPARAMS.h"
#include "GUD_TRAITS.h"
#include "GUD_FIELDS.h"

C     === Global variables ===
      _RL myTime
      INTEGER myIter
      INTEGER myThid

#ifdef ALLOW_GUD

C============== Local variables ============================================
      INTEGER i,j,k,bi,bj
      INTEGER np, iTr
      _RL totflxC
#ifdef PQUOTA
      _RL totflxP
#endif
#ifdef NQUOTA
      _RL totflxN
#endif
#ifdef FEQUOTA
      _RL totflxFe
#endif
#ifdef SIQUOTA
      _RL totflxSi
#endif
#ifdef ALLOW_CARBON
      _RL totflxA
      _RL totflxO
#endif
      _RL vol, voltot

      totflxC = 0. _d 0
#ifdef PQUOTA
      totflxP = 0. _d 0
#endif
#ifdef NQUOTA
      totflxN = 0. _d 0
#endif
#ifdef FEQUOTA
      totflxFe = 0. _d 0
#endif
#ifdef SIQUOTA
      totflxSi = 0. _d 0
#endif
#ifdef ALLOW_CARBON
      totflxA = 0. _d 0
      totflxO = 0. _d 0
#endif

c carbon budget
      totflxC=totflxC+totSurfCorPtr(iDIC)   !DIC
      totflxC=totflxC+totSurfCorPtr(iDOC)   !DOC
      totflxC=totflxC+totSurfCorPtr(iPOC)   !POC
      totflxC=totflxC+totSurfCorPtr(iPIC)   !PIC
      DO np=1,nplank
        iTr=ic+np-1
        totflxC=totflxC+totSurfCorPtr(iTr)  !Phy
      ENDDO
#ifdef ALLOW_CDOM
      totflxC=totflxC+totSurfCorPtr(iCDOM)*R_CP_CDOM  !CDOM
#endif

#ifdef PQUOTA
c Phosphorus budget
      totflxP=totflxP+totSurfCorPtr(iPO4)    !PO4
      totflxP=totflxP+totSurfCorPtr(iDOP)    !DOP
      totflxP=totflxP+totSurfCorPtr(iPOP)   !POP
      DO np=1,nplank
       IF (np .LE. npquota) THEN
        iTr=ip+np-1
        totflxP=totflxP+totSurfCorPtr(iTr)  !Phy 
       ELSE
        iTr=ic+np-1
        totflxP=totflxP+totSurfCorPtr(iTr)*R_PC(np)  !Phy 
       ENDIF
      ENDDO
#ifdef ALLOW_CDOM
      totflxP=totflxP+totSurfCorPtr(iCDOM)     !CDOM
#endif
#endif

#ifdef NQUOTA
c Nitrogen budget
      totflxN=totflxN+totSurfCorPtr(iNO3)    !NO3
      totflxN=totflxN+totSurfCorPtr(iNH4)   !NH4
      totflxN=totflxN+totSurfCorPtr(iNO2)   !NO2
      totflxN=totflxN+totSurfCorPtr(iDON)    !DON
      totflxN=totflxN+totSurfCorPtr(iPON)   !PON
      DO np=1,nplank
       IF (np .LE. nnquota) THEN
        iTr=in+np-1
        totflxN=totflxN+totSurfCorPtr(iTr)  !Phy
       ELSE
        iTr=ic+np-1
        totflxN=totflxN+totSurfCorPtr(iTr)*R_NC(np)  !Phy
       ENDIF
      ENDDO
#ifdef ALLOW_CDOM
      totflxN=totflxN+totSurfCorPtr(iCDOM)*R_NP_CDOM  !CDOM
#endif
#endif

#ifdef FEQUOTA
c Iron budget
      totflxFe=totflxFe+totSurfCorPtr(iFeT)    !Fet
      totflxFe=totflxFe+totSurfCorPtr(iDOFe)   !DOFe
      totflxFe=totflxFe+totSurfCorPtr(iPOFe)   !POFe
      DO np=1,nplank
       IF (np .LE. nfequota) THEN
        iTr=ife+np-1
        totflxFe=totflxFe+totSurfCorPtr(iTr)  !Phy
       ELSE
        iTr=ic+np-1
        totflxFe=totflxFe+totSurfCorPtr(iTr)*R_FeC(np)  !Phy
       ENDIF
      ENDDO
#ifdef ALLOW_CDOM
      totflxFe=totflxFe+totSurfCorPtr(iCDOM)*R_FeP_CDOM  !CDOM
#endif
#endif

#ifdef SIQUOTA
c Silica budget
      totflxSi=totflxSi+totSurfCorPtr(iSiO2)   !Si
      totflxSi=totflxSi+totSurfCorPtr(iPOSi)   !POSi
      DO np=1,nplank
       IF (np .LE. nsiquota) THEN
        iTr=isi+np-1
        totflxSi=totflxSi+totSurfCorPtr(iTr)  !Phy
       ELSE
        iTr=ic+np-1
        totflxSi=totflxSi+totSurfCorPtr(iTr)*R_SiC(np)  !Phy
       ENDIF
      ENDDO
c------------------------------------------------------------
#ifdef ALLOW_CARBON
c alkalinity budget !!!!! NOT FULL BUDGET!!!!!
      totflxA=totflxA+totSurfCorPtr(iALK)    !ALK
c oxygem budget !!!!!!!! NOT FULL BUDGET !!!!!!!!!!!
      totflxO=totflxO+totSurfCorPtr(iO2)    !O2
#endif
#endif

      voltot = 0. _d 0
      DO bj=myByLo(myThid),myByHi(myThid)
      DO bi=myBxLo(myThid),myBxHi(myThid)
        DO k=1,Nr
        DO i=1,sNx
        DO j=1,sNy
          vol=rA(i,j,bi,bj)*drF(k)*hFacC(i,j,k,bi,bj)
          voltot=voltot+vol
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO
      _GLOBAL_SUM_RL(voltot,myThid)

C--   Everyone else must wait
      _BARRIER

C     mean tracer change
      totflxC  = ptracers_dTLev(1)*totflxC /voltot
#ifdef NQUOTA
      totflxN  = ptracers_dTLev(1)*totflxN /voltot
#endif
#ifdef PQUOTA
      totflxP  = ptracers_dTLev(1)*totflxP /voltot
#endif
#ifdef FEQUOTA
      totflxFe = ptracers_dTLev(1)*totflxFe/voltot
#endif
#ifdef SIQUOTA
      totflxSi = ptracers_dTLev(1)*totflxSi/voltot
#endif
#ifdef ALLOW_CARBON
      totflxA  = ptracers_dTLev(1)*totflxA /voltot
      totflxO  = ptracers_dTLev(1)*totflxO /voltot
#endif

C     add neglected flux to tracer globally
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO k=1,Nr
         DO i=1,sNx
          DO j=1,sNy
           Ptracer(i,j,k,bi,bj,iDIC)=Ptracer(i,j,k,bi,bj,iDIC)+
     &                               totflxC*maskC(i,j,k,bi,bj)
#ifdef PQUOTA
           Ptracer(i,j,k,bi,bj,iPO4)=Ptracer(i,j,k,bi,bj,iPO4)+
     &                               totflxP *maskC(i,j,k,bi,bj)
#endif
#ifdef NQUOTA
           Ptracer(i,j,k,bi,bj,iNO3)=Ptracer(i,j,k,bi,bj,iNO3)+
     &                               totflxN *maskC(i,j,k,bi,bj)
#endif
#ifdef FEQUOTA
           Ptracer(i,j,k,bi,bj,iFeT)=Ptracer(i,j,k,bi,bj,iFeT)+
     &                               totflxFe*maskC(i,j,k,bi,bj)
#endif
#ifdef SIQUOTA
           Ptracer(i,j,k,bi,bj,iSi )=Ptracer(i,j,k,bi,bj,iSi )+
     &                               totflxSi*maskC(i,j,k,bi,bj)
#endif
#ifdef ALLOW_CARBON
           Ptracer(i,j,k,bi,bj,iALK)=Ptracer(i,j,k,bi,bj,iALK)+
     &                               totflxA*maskC(i,j,k,bi,bj)
           Ptracer(i,j,k,bi,bj,iO2 )=Ptracer(i,j,k,bi,bj,iO2 )+
     &                               totflxO*maskC(i,j,k,bi,bj)
#endif
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#endif /* ALLOW_GUD */

      RETURN
      END

