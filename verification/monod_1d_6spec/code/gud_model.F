C $Header$
C $Name$
C TODO:
C - carbon
C - cdom
C - conservation of Chl?
C - PIC: grazers must have R_PICPOC = 0!
C   need to set all R_PICPOC to zero to emulate monod without carbon
C   (maybe not)
C   with carbon???
C - diazotrophy with quota

#include "GUD_OPTIONS.h"

CBOP
C !ROUTINE: GUD_MODEL
C !INTERFACE: ==========================================================
      SUBROUTINE GUD_MODEL(
     I     Ptr,
     U     gTr,
     O     chlout,
     I     PAR, photoTempFunc, reminTempFunc, uptakeTempFunc,
#ifdef GUD_DEBUG
     I     iG, jG, k,
#endif
     I     myTime,myIter,myThid)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"
#include "GUD_GENPARAMS.h"
#include "GUD_TRAITS.h"

C !INPUT PARAMETERS: ===================================================
C  Ptr    :: gud model tracers
C  PAR    :: PAR in uEin/s/m2
C         :: either non-spectral (tlam=1) or waveband total
C  myTime :: current time
C  myIter :: current iteration number
C  myThid :: thread number
      _RL Ptr(nGud)
      _RL PAR(nlam)
      _RL photoTempFunc(nplank)
      _RL reminTempFunc
      _RL uptakeTempFunc
      INTEGER myThid, myIter
      _RL myTime
#ifdef GUD_DEBUG
      INTEGER iG, jG
#endif

C !INPUT/OUTPUT PARAMETERS: ============================================
C  gTr    :: accumulates computed tendencies
      _RL gTr(nGud)

C !OUTPUT PARAMETERS: ==================================================
C  chlout :: computed acclimated chlorophyll if not dynamic
      _RL chlout(nPhoto)
CEOP

#ifdef ALLOW_GUD

c !LOCAL VARIABLES: ====================================================
      INTEGER j, jy, l
#ifdef ALLOW_DIAGNOSTICS
      CHARACTER*8   diagName
#endif

      _RL DIC
      _RL NH4
      _RL NO2
      _RL NO3
      _RL PO4
      _RL SiO2
      _RL FeT
      _RL DOC
      _RL DON
      _RL DOP
      _RL DOFe
      _RL POC
      _RL PON
      _RL POP
      _RL POSi
      _RL POFe
      _RL PIC
      _RL O2

      _RL X
      _RL Qc
#ifdef NQUOTA
      _RL Qn
#endif
#ifdef PQUOTA
      _RL Qp
#endif
#ifdef SIQUOTA
      _RL Qsi
#endif
#ifdef FEQUOTA
      _RL Qfe
#endif
#ifdef CHLQUOTA
      _RL QChl
#endif

      _RL regQ

      _RL limitpCO2
      _RL limitNH4
      _RL limitNO2
      _RL limitNO3
      _RL fracNH4
      _RL fracNO2
      _RL fracNO3
      _RL limitn
      _RL limitp
      _RL limitsi
      _RL limitfe
      _RL limitnut
      _RL limitI
      _RL ngrow

      _RL uptakeDIC
      _RL uptakeNH4
      _RL uptakeNO2
      _RL uptakeNO3
      _RL uptakeN
      _RL uptakePO4
      _RL uptakeSiO2
      _RL uptakeFeT
      _RL consumDIC
      _RL consumDIC_PIC
      _RL consumNH4
      _RL consumNO2
      _RL consumNO3
      _RL consumPO4
      _RL consumSiO2
      _RL consumFeT

      _RL inhibNH4

      _RL alpha_I
      _RL alpha_I_growth
      _RL PCm
      _RL PC
      _RL acclim
      _RL chl2c
      _RL growth
      _RL rhochl
      _RL Ek
      _RL EkoverE

      _RL synthChl

      _RL reminDOC
      _RL reminDON
      _RL reminDOP
      _RL reminDOFe
      _RL reminPOC
      _RL reminPON
      _RL reminPOP
      _RL reminPOSi
      _RL reminPOFe
      _RL disscPIC

      _RL prodNO2
      _RL prodNO3

      _RL PARtot

      _RL tmp

#ifdef ALLOW_CDOM
      _RL CDOM
      _RL reminPOC_CDOM
      _RL reminPON_CDOM
      _RL reminPOP_CDOM
      _RL reminPOFe_CDOM
      _RL degrCDOM_DOC
      _RL degrCDOM_DON
      _RL degrCDOM_DOP
      _RL degrCDOM_DOFe
#endif

#ifdef ALLOW_DENIT
      _RL denit, denitNH4
#endif

      _RL EPS
      PARAMETER (EPS=1D-38)

      INTEGER previter,k
      SAVE previter,k
      DATA previter/-1/

      IF (myiter .ne. previter) then
        k = 0
        previter = myiter
      ENDIF
      k = k + 1

C==== precompute a few things ==========================================

      PARtot = SUM(PAR)

C==== make all bio fields non-negative and compute quotas ==============

      DIC  = MAX(0., Ptr(iDIC))
      NH4  = MAX(0., Ptr(iNH4))
      NO2  = MAX(0., Ptr(iNO2))
      NO3  = MAX(0., Ptr(iNO3))
      PO4  = MAX(EPS, Ptr(iPO4))
      SiO2 = MAX(EPS, Ptr(iSiO2))
      FeT  = MAX(EPS, Ptr(iFeT))

      DOC  = MAX(0., Ptr(iDOC))
      DON  = MAX(0., Ptr(iDON))
      DOP  = MAX(0., Ptr(iDOP))
      DOFe = MAX(0., Ptr(iDOFe))

      PIC  = MAX(0., Ptr(iPIC))
      POC  = MAX(0., Ptr(iPOC))
      PON  = MAX(0., Ptr(iPON))
      POP  = MAX(0., Ptr(iPOP))
      POSi = MAX(0., Ptr(iPOSi))
      POFe = MAX(0., Ptr(iPOFe))
#ifdef ALLOW_CARBON
      O2   = MAX(0., Ptr(iO2))
#endif
#ifdef ALLOW_CDOM
      CDOM = MAX(0., Ptr(iCDOM))
#endif

      chlout = 0.0 _d 0

      consumDIC  = 0.0
      consumDIC_PIC = 0.0
      consumNH4  = 0.0
      consumNO2  = 0.0
      consumNO3  = 0.0
      consumPO4  = 0.0
      consumSiO2 = 0.0
      consumFeT  = 0.0

C fixed carbon quota, for now 1.0 (may change later)
C other elements: get quota from corresponding ptracer or set to fixed
c ratio if not variable.
      DO j = iMinPhoto, iMaxPhoto
        X = MAX(0.0, Ptr(ic+j-1))
        Qc = 1.0

C==== uptake and nutrient limitation ===================================
C       for quota elements, growth is limiteed by available quota,
C       for non-quota elements, by available nutrients in medium

C       to not use PO4, ..., set ksatPO4=0 and Vmax_PO4=0 (if PQUOTA)
C       or R_PC=0 (if not)
C       the result will be limitp = 1, uptakePO4 = 0

c phosphorus
         limitp = PO4/(PO4 + ksatPO4(j))
#ifdef PQUOTA
         Qp  = MAX(EPS*R_PC(j), Ptr(ip+j-1)) / MAX(EPS, X)
         regQ = MAX(0., MIN(1., (Qpmax(j)-Qp)/
     &                          (Qpmax(j)-Qpmin(j)) ))
         uptakePO4 = Vmax_PO4(j) * limitp * regQ *
     &                        uptakeTempFunc * X
c        normalized Droop limitation
         limitp = MAX(0., MIN(1., (1.0-Qpmin(j)/Qp)/
     &                            (1.0-Qpmin(j)/Qpmax(j))))
#endif

c silica
         limitsi = SiO2/(SiO2 + ksatSiO2(j))
#ifdef SIQUOTA
         Qsi = MAX(EPS*R_SiC(j), Ptr(isi+j-1)) / MAX(EPS, X)
         regQ = MAX(0., MIN(1., (Qsimax(j) - Qsi)/
     &                          (Qsimax(j) - Qsimin(j)) ))
         uptakeSiO2 = Vmax_SiO2(j) * limitsi * regQ *
     &                                  uptakeTempFunc * X

c        linear limitation
         limitsi = MAX(0., MIN(1., (Qsi - Qsimin(j))/
     &                          (Qsimax(j) - Qsimin(j)) ))
#endif

c iron
         limitfe = FeT/(FeT + ksatFeT(j))
#ifdef FEQUOTA
         Qfe = MAX(EPS*R_FeC(j), Ptr(ife+j-1)) / MAX(EPS, X)
         regQ = MAX(0., MIN(1., (Qfemax(j)-Qfe)/
     &                          (Qfemax(j)-Qfemin(j)) ))
         uptakeFeT = Vmax_FeT(j) * limitfe * regQ *
     &                        uptakeTempFunc * X

c        normalized Droop limitation
         limitfe = MAX(0., MIN(1., (1.0-Qfemin(j)/Qfe)/
     &                             (1.0-Qfemin(j)/Qfemax(j))))
#endif

c nitrogen
#ifdef NQUOTA
c       have nitrogen quota
        inhibNH4 = EXP(-amminhib(j)*NH4)
        limitNH4 = NH4/(NH4 + ksatNH4(j))
        limitNO2 = NO2/(NO2 + ksatNO2(j))*inhibNH4
        limitNO3 = NO3/(NO3 + ksatNO3(j))*inhibNH4
        Qn  = MAX(EPS*R_NC(j), Ptr(in+j-1)) / MAX(EPS, X)
        regQ = MAX(0., MIN(1., (Qnmax(j)-Qn)/
     &                         (Qnmax(j)-Qnmin(j)) ))
        uptakeNH4 = Vmax_NH4(j)*limitNH4*regQ*uptakeTempFunc*X
        uptakeNO2 = Vmax_NO2(j)*limitNO2*regQ*uptakeTempFunc*X
        uptakeNO3 = Vmax_NO3(j)*limitNO3*regQ*uptakeTempFunc*X
#if USE_FELIMIT_PS
        uptakeNO3 = uptakeNO3 * limitfe
#endif

c       linear limitation
        limitn = MAX(0., MIN(1., (Qn - Qnmin(j))/
     &                        (Qnmax(j) - Qnmin(j)) ))
#else /* not NQUOTA */
        inhibNH4 = EXP(-amminhib(j)*NH4)
        limitNH4 = useNH4(j)*NH4/(NH4 + ksatNH4(j))
        limitNO2 = useNO2(j)*NO2/
     &   (NO2 + combNO(j)*(NO3 + ksatNO3(j) - ksatNO2(j)) + ksatNO2(j))*
     &   inhibNH4
        limitNO3 = useNO3(j)*NO3/
     &   (combNO(j)*NO2 + NO3 + ksatNO3(j))*inhibNH4
        limitn = limitNH4 + limitNO2 + limitNO3
C       normalize to sum (approx) 1
        fracNH4 = limitNH4/(limitn + EPS)
        fracNO2 = limitNO2/(limitn + EPS)
        fracNO3 = limitNO3/(limitn + EPS)
C if diazo, all fracN* == 0 but want no N limitation
        limitn = MIN(1.0, limitn + diazotroph(j))
#endif /* NQUOTA */
c        IF (limitn .GT. 0.0) THEN
c          ngrow = ((10*4+2)/(10*4 + 2*limitNH4/limitn +
c     &                       8*limitNO2/limitn + 10*limitNO3/limitn))
c        ELSE
        ngrow = 1.0
c        ENDIF
      print*,'limit',myiter,k,j,limitp,limitfe,limitsi,limitn,
     &     fracNH4,fracNO2,fracNO3

        limitnut = MIN(limitn, limitp, limitsi)
#if ! USE_FELIMIT_PS
        limitnut = MIN(limitnut, limitfe)
#endif

        limitpCO2 = 1.

C==== growth ===========================================================
#ifdef GEIDER

        alpha_I = 0 _d 0
        DO l = 1, nlam
          alpha_I = alpha_I + alphachl(j,l)*PAR(l)
#ifdef GUD_DEBUG
      if (ig.eq.iDEBUG.and.jg.eq.jDEBUG) then
       print*,'alphaI',myiter,k,j,l,alpha_I,alphachl(j,l),
     &       PAR(l)
      ENDIF
#endif
        ENDDO
C       NB: for quota, PCmax(j) = Vmax_c(j)
        PCm = PCmax(j)*limitnut*photoTempFunc(j)*limitpCO2
#ifdef GUD_DEBUG
      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG)
     & print*,'pcm',myiter,k,j,PCmax(j),limitnut,phototempfunc(j),
     &      limitpco2
#endif

        IF (PCm .GT. 0.0) THEN
          acclim = MAX(chl2cmin(j), MIN(chl2cmax(j),
     &             chl2cmax(j)/(1+(chl2cmax(j)*alpha_I)/(2*PCm)) ))
        ELSE
          acclim = chl2cmin(j)
        ENDIF

#ifdef CHLQUOTA
        QChl=MAX(EPS*R_ChlC(j),Ptr(ichl+j-iMinPhoto))/MAX(EPS,X)
C       quotas are already relative to carbon
        chl2c = QChl
#ifdef GUD_DEBUG
      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG)
     & print*,'chl2c',myiter,k,j,chl2c,EPS,R_ChlC(j),
     &     Ptr(ichl+j-iMinPhoto),X

#endif
#else
        chl2c = acclim
#endif

        alpha_I_growth = alpha_I
C a la quota
#if USE_FELIMIT_PS
        alpha_I_growth = alpha_I_growth*limitfe
#endif

C       carbon-specific growth rate
C       PC = PCm*(1-EXP(-alpha_I_growth*chl2c/MAX(EPS, PCm)))
        IF (PCm .GT. 0.0 .AND. PARtot .GT. PARmin) THEN
          PC = PCm*(1-EXP(-alpha_I_growth*chl2c/PCm))
        ELSE
          PC = 0.0 _d 0
        ENDIF
#ifdef GUD_DEBUG
      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG)
     & print*,'pc',myiter,k,j,PCm,alpha_I_growth,chl2c,PC
#endif

        IF (inhibcoef_geid(j) .GT. 0.0) THEN
C         "total" PAR:
          tmp = alpha_I/alpha_mean(j)
          Ek = PCm/(chl2c*alpha_mean(j))
          EkoverE = Ek / tmp
          IF (tmp .GE. Ek) THEN
            PC = PC*EkoverE*inhibcoef_geid(j)
          ENDIF
        ENDIF

#else /* not GEIDER */

        IF (PARtot .GT. PARmin) THEN
C         only 1 waveband without GEIDER
          limitI = (1.0 _d 0 - EXP(-PARtot*ksatPAR(j)))*
     &             EXP(-PARtot*kinhPAR(j)) * normI(j)
          PC = PCmax(j)*limitnut*limitI*photoTempFunc(j)*limitpCO2
        ELSE
          PC = 0.0 _d 0
        ENDIF
        synthChl = 0.0

#endif /* GEIDER */

        growth = PC*ngrow*X

        uptakeDIC = growth

      print*,'chl2c',myiter,k,j,chl2c,PC,ngrow

C non-quota elements are taken up with fixed stoichiometry
#ifndef NQUOTA
          uptakeN = growth*R_NC(j)
          uptakeNH4 = uptakeN*fracNH4
          uptakeNO2 = uptakeN*fracNO2
          uptakeNO3 = uptakeN*fracNO3
#endif
#ifndef PQUOTA
         uptakePO4 = growth*R_PC(j)
#endif
#ifndef SIQUOTA
         uptakeSiO2 = growth*R_SiC(j)
#endif
#ifndef FEQUOTA
         uptakeFeT = growth*R_FeC(j)
#endif

C==== chlorophyll ======================================================
#ifdef GEIDER
#ifdef CHLQUOTA
#ifdef NQUOTA
C         Geider 1998
          IF (alpha_I*chl2c .GT. 0.0 _d 0) THEN
C            rhochl = Chl2Nmax/(alpha_I*chl2c)*ngrow ???
            rhochl = Chl2Nmax*PC*ngrow/(alpha_I*chl2c)
          ELSE
            rhochl = Chl2Nmax
          ENDIF
          uptakeN = uptakeNH4 + uptakeNO2 +
     &              uptakeNO3
          uptakeDIC = uptakeDIC - synthcost*uptakeN
          synthChl = rhochl*uptakeN
#else
c          synthChl = acclim*growth +
c     &                        acclimtimescl(j)*(acclim-chl2c)*X
          IF (alpha_I .GT. 0.0 _d 0 .AND. acclim .GT. 0.0 _d 0) THEN
            rhochl = Chl2Cmax(j)*PC*ngrow/(alpha_I*acclim)
          ELSE
            rhochl = 0.0 _d 0    ! should be Chl2Cmax(j) ?????
          ENDIF
          synthChl = rhochl*growth +
     &                        acclimtimescl(j)*(acclim-chl2c)*X
#ifdef GUD_DEBUG
      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
       print*,'rhochl',myiter,k,j,rhochl,chl2cmax(j),pc,alpha_I,
     &      acclim,ngrow
       print*,'synth',myiter,k,j,synthchl,rhochl,growth,
     &       acclimtimescl(j),acclim,chl2c,X
      ENDIF
#endif
#endif /* NQUOTA */
#else /* CHLQUOTA */
         chlout(j-iMinPhoto+1) = X*Qc*chl2c
         synthChl = 0.0
#endif /* CHLQUOTA */
#endif /* GEIDER */
C=======================================================================
        consumDIC_PIC = consumDIC_PIC + uptakeDIC*R_PICPOC(j)
        consumDIC  = consumDIC  + uptakeDIC 
        consumNH4  = consumNH4  + uptakeNH4 
        consumNO2  = consumNO2  + uptakeNO2 
        consumNO3  = consumNO3  + uptakeNO3 
        consumPO4  = consumPO4  + uptakePO4 
        consumSiO2 = consumSiO2 + uptakeSiO2
        consumFeT  = consumFeT  + uptakeFeT 

C=======================================================================

       gTr(ic+j-1)=gTr(ic+j-1)  + uptakeDIC
#ifdef NQUOTA
       gTr(in+j-1)=gTr(in+j-1)  + uptakeNH4 + uptakeNO2 + uptakeNO3
#endif
#ifdef PQUOTA
       gTr(ip+j-1)=gTr(ip+j-1)  + uptakePO4
#endif
#ifdef SIQUOTA
       gTr(isi+j-1)=gTr(isi+j-1) + uptakeSiO2
#endif
#ifdef FEQUOTA
       gTr(ife+j-1)=gTr(ife+j-1) + uptakeFeT
#endif
#ifdef CHLQUOTA
       gTr(iChl+j-iMinPhoto)=gTr(iChl+j-iMinPhoto) + synthChl
#endif

#ifdef GUD_DEBUG
      IF (iG.eq.iDEBUG.and.jG.eq.jDEBUG) THEN
       print*,'uptake',myiter,k,j,
     &     uptakeDIC,
     &     uptakeNH4,
     &     uptakeNO2,
     &     uptakeNO3,
     &     uptakePO4,
     &     uptakeSiO2,
     &     uptakeFeT
      ENDIF
#endif
C     j
      ENDDO

      gTr(iDIC )=gTr(iDIC ) - consumDIC - consumDIC_PIC
      print*,'consum',consumDIC,consumDIC_PIC
      gTr(iNH4 )=gTr(iNH4 ) - consumNH4
      gTr(iNO2 )=gTr(iNO2 ) - consumNO2
      gTr(iNO3 )=gTr(iNO3 ) - consumNO3
      gTr(iPO4 )=gTr(iPO4 ) - consumPO4
      gTr(iSiO2)=gTr(iSiO2) - consumSiO2
      gTr(iFeT )=gTr(iFeT ) - consumFeT

      reminDOC  = reminTempFunc*KDOC *DOC
      reminDON  = reminTempFunc*KDON *DON
      reminDOP  = reminTempFunc*KDOP *DOP
      reminDOFe = reminTempFunc*KDOFe*DOFe
      reminPOC  = reminTempFunc*KPOC *POC
      print*,'reminPOC',reminPOC,remintempfunc,kpoc,poc
      reminPON  = reminTempFunc*KPON *PON
      reminPOP  = reminTempFunc*KPOP *POP
      reminPOSi = reminTempFunc*KPOSi*POSi
      print*,'reminPOSi',reminPOSi,remintempfunc,kpoSi,poSi
      reminPOFe = reminTempFunc*KPOFe*POFe
      disscPIC = Kdissc*PIC

c nitrogen chemistry
c NH4 -> NO2 -> NO3 by bacterial action
      prodNO2 = knita*NH4
      prodNO3 = knitb*NO2
      IF (PAR_oxi .NE. 0.0 _d 0) THEN
c       print*,'prod',prodNO2,PARtot,PAR_oxi
        prodNO2 = prodNO2*MAX(0.0, 1.0 - PARtot/PAR_oxi)
        prodNO3 = prodNO3*MAX(0.0, 1.0 - PARtot/PAR_oxi)
      ENDIF

#ifdef ALLOW_CDOM
      reminPOP_CDOM = fracCDOM*reminPOP
      reminPOC_CDOM = R_CP_CDOM*reminPOP_CDOM
      reminPON_CDOM = R_NP_CDOM*reminPOP_CDOM
      reminPOFe_CDOM = R_FeP_CDOM*reminPOP_CDOM
c degradation of  CDOM - high when bleached by light
      degrCDOM_DOP = reminTempFunc*CDOM*
     &     (CDOMdegrd+CDOMbleach*MIN(1.0 _d 0, PARtot/PARCDOM))
      degrCDOM_DOC  = R_CP_CDOM  * degrCDOM_DOP
      degrCDOM_DON  = R_NP_CDOM  * degrCDOM_DOP
      degrCDOM_DOFe = R_FeP_CDOM * degrCDOM_DOP
#endif

#ifdef ALLOW_DENIT
      IF (O2 .LT. O2crit .AND. NO3 .LT. NO3crit) THEN
c no remineralization for N, P, Fe (not Si?)
        reminDOC = 0. _d 0
        reminDOP  = 0. _d 0
        reminDON  = 0. _d 0
        reminDOFe = 0. _d 0
        reminPOC = 0. _d 0
        reminPOP  = 0. _d 0
        reminPON  = 0. _d 0
        reminPOFe = 0. _d 0
#ifdef ALLOW_CDOM
        reminPOC_cdom = 0. _d 0
        reminPOP_cdom  = 0. _d 0
        reminPON_cdom  = 0. _d 0
        reminPOFe_cdom = 0. _d 0
        degrCDOM_DOP = reminTempFunc*CDOM*
     &                 CDOMbleach*MIN(1.0 _d 0, PARtot/PARCDOM)
        degrCDOM_DON  = R_NP_CDOM  * degrCDOM_DOP
        degrCDOM_DOFe = R_FeP_CDOM * degrCDOM_DOP
        degrCDOM_DOC  = R_CP_CDOM  * degrCDOM_DOP
#endif
      ENDIF
#endif /* ALLOW_DENIT */

C==== tendencies =======================================================

#ifdef ALLOW_CARBON
c production of O2 by photosynthesis
      gTr(iO2  )=gTr(iO2  ) + R_OP*consumPO4
c loss of O2 by remineralization
      IF (O2 .GT. O2crit) THEN
        gTr(iO2)=gTr(iO2) - R_OP*reminDOP
#ifndef ALLOW_CDOM
        gTr(iO2)=gTr(iO2) - R_OP*reminPOP
#endif
      ENDIF

      gTr(iALK)=gTr(iALK) - (prodNO3 - consumNO3)
     &                    - 2.0 _d 0*(consumDIC_PIC - disscPIC)
#endif /* ALLOW_CARBON */

      gTr(iDIC )=gTr(iDIC ) + reminDOC + disscPIC
      print*,'reminDOC',reminDOC,disscPIC
      gTr(iNH4 )=gTr(iNH4 ) + reminDON - prodNO2
      gTr(iNO2 )=gTr(iNO2 ) + prodNO2 - prodNO3
      gTr(iNO3 )=gTr(iNO3 )           + prodNO3
#ifdef ALLOW_DENIT
      IF (O2 .LT. O2crit) THEN
        denitNH4 = reminDON
        denit = denit_NP*reminDOP
#ifndef ALLOW_CDOM
        denitNH4 = denitNH4 + reminPON
        denit = denit + denit_NP*reminPOP
#endif
        gTr(iNH4)=gTr(iNH4) - denitNH4
        gTr(iNO3)=gTr(iNO3) - denit_NO3/denit_np*denit
        gTr(iALK)=gTr(iALK) + denit_NO3/denit_np*denit
      ENDIF
#endif /* ALLOW_DENIT */

      gTr(iPO4 )=gTr(iPO4 ) + reminDOP
      gTr(iFeT )=gTr(iFeT ) + reminDOFe
      gTr(iSiO2)=gTr(iSiO2)             + reminPOSi

      gTr(iDOC )=gTr(iDOC ) - reminDOC
      gTr(iDON )=gTr(iDON ) - reminDON
      gTr(iDOP )=gTr(iDOP ) - reminDOP
      gTr(iDOFe)=gTr(iDOFe) - reminDOFe

      gTr(iPIC )=gTr(iPIC ) - disscPIC
      gTr(iPOC )=gTr(iPOC ) - reminPOC
      gTr(iPON )=gTr(iPON ) - reminPON
      gTr(iPOP )=gTr(iPOP ) - reminPOP
      gTr(iPOSi)=gTr(iPOSi) - reminPOSi
      gTr(iPOFe)=gTr(iPOFe) - reminPOFe
#ifdef ALLOW_CDOM
      gTr(iDOC )=gTr(iDOC ) + reminPOC  - reminPOC_CDOM  + degrCDOM_DOC
      gTr(iDON )=gTr(iDON ) + reminPON  - reminPON_CDOM  + degrCDOM_DON
      gTr(iDOP )=gTr(iDOP ) + reminPOP  - reminPOP_CDOM  + degrCDOM_DOP
      gTr(iDOFe)=gTr(iDOFe) + reminPOFe - reminPOFe_CDOM + degrCDOM_DOFe

      gTr(iCDOM)=gTr(iCDOM) + reminPOP_CDOM - degrCDOM_DOP
#else
      gTr(iDIC )=gTr(iDIC ) + reminPOC
      gTr(iNH4 )=gTr(iNH4 ) + reminPON
      gTr(iPO4 )=gTr(iPO4 ) + reminPOP
      gTr(iFeT )=gTr(iFeT ) + reminPOFe
#endif /* ALLOW_CDOM */

#endif /* ALLOW_GUD */

      RETURN
      END SUBROUTINE
