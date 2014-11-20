C $Header$
C $Name$

#include "GUD_OPTIONS.h"

c=============================================================
c subroutine GUD_cons_tend
c check conservation in the model
C==============================================================
      SUBROUTINE GUD_CONS_TEND(myTime,myIter,myThid)
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"
#include "GUD_TRAITS.h"
#include "GUD_FIELDS.h"

C     === Global variables ===
      _RL myTime
      INTEGER myIter
      INTEGER myThid

#ifdef ALLOW_GUD_CONS

C============== Local variables ============================================
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      _RL tmptotC
      _RL tmptotP
      _RL tmptotN
      _RL tmptotSi
      _RL tmptotFe
#ifdef ALLOW_CARBON
      _RL tmptotA
      _RL tmptotO
#endif
      _RL vol, voltot

      INTEGER i,j,k,bi,bj
      INTEGER np, iTr
c
      tmptotC = 0. _d 0
      tmptotP = 0. _d 0
      tmptotN = 0. _d 0
      tmptotFe = 0. _d 0
      tmptotSi = 0. _d 0
#ifdef ALLOW_CARBON
      tmptotA = 0. _d 0
      tmptotO = 0. _d 0
#endif
      voltot = 0. _d 0

      DO bj=myByLo(myThid),myByHi(myThid)
      DO bi=myBxLo(myThid),myBxHi(myThid)
 
        DO k=1,Nr
        DO i=1,sNx
        DO j=1,sNy
c-----------------------------------------------------------------
          vol=rA(i,j,bi,bj)*drF(k)*hFacC(i,j,k,bi,bj)
          voltot=voltot+vol

c carbon budget
          tmptotC=tmptotC+gPtr(i,j,k,bi,bj,iDIC)*vol   !DIC
          tmptotC=tmptotC+gPtr(i,j,k,bi,bj,iDOC)*vol   !DOC
          tmptotC=tmptotC+gPtr(i,j,k,bi,bj,iPOC)*vol   !POC
          tmptotC=tmptotC+gPtr(i,j,k,bi,bj,iPIC)*vol   !PIC
          DO np=1,nplank
            iTr=ic+np-1
            tmptotC=tmptotC+gPtr(i,j,k,bi,bj,iTr)*vol  !Phy
          ENDDO
#ifdef ALLOW_CDOM
          tmptotC=tmptotC+gPtr(i,j,k,bi,bj,iCDOM)*vol  !CDOM
#endif

c Phosphorus budget
          tmptotP=tmptotP+gPtr(i,j,k,bi,bj,iPO4)*vol    !PO4
          tmptotP=tmptotP+gPtr(i,j,k,bi,bj,iDOP)*vol    !DOP
          tmptotP=tmptotP+gPtr(i,j,k,bi,bj,iPOP)*vol   !POP
          DO np=1,nplank
#ifdef PQUOTA
            iTr=ip+np-1
            tmptotP=tmptotP+gPtr(i,j,k,bi,bj,iTr)*vol  !Phy 
#else
            iTr=ic+np-1
            tmptotP=tmptotP+gPtr(i,j,k,bi,bj,iTr)*R_PC(np)*vol  !Phy 
#endif
          ENDDO
#ifdef ALLOW_CDOM
          tmptotP=tmptotP+gPtr(i,j,k,bi,bj,iCDOM)*R_PC_CDOM*vol     !CDOM
#endif

c Nitrogen budget
          tmptotN=tmptotN+gPtr(i,j,k,bi,bj,iNO3)*vol    !NO3
          tmptotN=tmptotN+gPtr(i,j,k,bi,bj,iNH4)*vol   !NH4
          tmptotN=tmptotN+gPtr(i,j,k,bi,bj,iNO2)*vol   !NO2
          tmptotN=tmptotN+gPtr(i,j,k,bi,bj,iDON)*vol    !DON
          tmptotN=tmptotN+gPtr(i,j,k,bi,bj,iPON)*vol   !PON
          DO np=1,nplank
#ifdef NQUOTA
            iTr=in+np-1
            tmptotN=tmptotN+gPtr(i,j,k,bi,bj,iTr)*vol  !Phy
#else
            iTr=ic+np-1
            tmptotN=tmptotN+gPtr(i,j,k,bi,bj,iTr)*R_NC(np)*vol  !Phy
#endif
          ENDDO
#ifdef ALLOW_CDOM
          tmptotN=tmptotN+gPtr(i,j,k,bi,bj,iCDOM)*R_NC_CDOM*vol  !CDOM
#endif

c Iron budget
          tmptotFe=tmptotFe+gPtr(i,j,k,bi,bj,iFeT)*vol    !Fet
          tmptotFe=tmptotFe+gPtr(i,j,k,bi,bj,iDOFe)*vol   !DOFe
          tmptotFe=tmptotFe+gPtr(i,j,k,bi,bj,iPOFe)*vol   !POFe
          DO np=1,nplank
#ifdef FEQUOTA
            iTr=ife+np-1
            tmptotFe=tmptotFe+gPtr(i,j,k,bi,bj,iTr)*vol  !Phy
#else
            iTr=ic+np-1
            tmptotFe=tmptotFe+gPtr(i,j,k,bi,bj,iTr)*R_FeC(np)*vol  !Phy
#endif
          ENDDO
#ifdef ALLOW_CDOM
          tmptotFe=tmptotFe+gPtr(i,j,k,bi,bj,iCDOM)*R_FeC_CDOM*vol  !CDOM
#endif

c Silica budget
          tmptotSi=tmptotSi+gPtr(i,j,k,bi,bj,iSiO2)*vol   !Si
          tmptotSi=tmptotSi+gPtr(i,j,k,bi,bj,iPOSi)*vol   !POSi
          DO np=1,nplank
#ifdef SIQUOTA
            iTr=isi+np-1
            tmptotSi=tmptotSi+gPtr(i,j,k,bi,bj,iTr)*vol  !Phy
#else
            iTr=ic+np-1
            tmptotSi=tmptotSi+gPtr(i,j,k,bi,bj,iTr)*R_SiC(np)*vol  !Phy
#endif
          ENDDO
c------------------------------------------------------------
#ifdef ALLOW_CARBON
c alkalinity budget !!!!! NOT FULL BUDGET!!!!!
          tmptotA=tmptotA+gPtr(i,j,k,bi,bj,iALK)*vol    !ALK
c oxygem budget !!!!!!!! NOT FULL BUDGET !!!!!!!!!!!
          tmptotO=tmptotO+gPtr(i,j,k,bi,bj,iO2)*vol    !O2
#endif
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO

      _GLOBAL_SUM_RL(tmptotC,myThid)
      _GLOBAL_SUM_RL(tmptotP,myThid)
      _GLOBAL_SUM_RL(tmptotN,myThid)
      _GLOBAL_SUM_RL(tmptotFe,myThid)
      _GLOBAL_SUM_RL(tmptotSi,myThid)
      _GLOBAL_SUM_RL(voltot,myThid)
#ifdef ALLOW_CARBON
      _GLOBAL_SUM_RL(tmptotA,myThid)
      _GLOBAL_SUM_RL(tmptotO,myThid)
#endif

#define PMSGR(_MSG) CALL PRINT_MESSAGE(_MSG,standardMessageUnit,SQUEEZE_RIGHT,1)

      IF ( myProcId.EQ.0 .AND. myThid.EQ.1 ) THEN
        WRITE(msgBuf,'(A,1P2E25.16)') 'totgC ',tmptotC ,tmptotC /voltot
        PMSGR(msgBuf)
        WRITE(msgBuf,'(A,1P2E25.16)') 'totgN ',tmptotN ,tmptotN /voltot
        PMSGR(msgBuf)
        WRITE(msgBuf,'(A,1P2E25.16)') 'totgP ',tmptotP ,tmptotP /voltot
        PMSGR(msgBuf)
        WRITE(msgBuf,'(A,1P2E25.16)') 'totgFe',tmptotFe,tmptotFe/voltot
        PMSGR(msgBuf)
        WRITE(msgBuf,'(A,1P2E25.16)') 'totgSi',tmptotSi,tmptotSi/voltot
        PMSGR(msgBuf)
#ifdef ALLOW_CARBON
        WRITE(msgBuf,'(A,1P2E25.16)') 'totgA ',tmptotA ,tmptotA /voltot
        PMSGR(msgBuf)
        WRITE(msgBuf,'(A,1P2E25.16)') 'totgO ',tmptotO ,tmptotO /voltot
        PMSGR(msgBuf)
#endif
        CALL MDS_FLUSH(standardMessageUnit, myThid)
      ENDIF

C--   Everyone else must wait
      _BARRIER

#endif /* ALLOW_GUD_CONS */

      RETURN
      END

