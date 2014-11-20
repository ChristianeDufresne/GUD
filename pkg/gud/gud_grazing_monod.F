C $Header$
C $Name$

#include "GUD_OPTIONS.h"

CBOP
C !ROUTINE: GUD_GRAZING_MONOD
C !INTERFACE: =========================================================
      SUBROUTINE GUD_GRAZING_MONOD(
     I     Ptr,
     U     gTr,
     I     grazTempFunc, reminTempFunc, mortTempFunc, mort2TempFunc,
     I     myTime,myIter,myThid)

C !DESCRIPTION:
C     Add monod-style grazing tendencies to gPtr

C !USES: ===============================================================
      IMPLICIT NONE
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"
#include "GUD_GENPARAMS.h"
#include "GUD_TRAITS.h"

C !INPUT PARAMETERS: ===================================================
C  Ptr    :: gud model tracers
C  Temp   :: temperature field (degrees C)
C  myTime :: current time
C  myIter :: current iteration number
C  myThid :: thread number
      _RL Ptr(nGud)
      _RL grazTempFunc(nplank)
      _RL mortTempFunc
      _RL mort2TempFunc
      _RL reminTempFunc
      _RL myTime
      INTEGER myIter, myThid

C !INPUT/OUTPUT PARAMETERS: ============================================
C  gTr    :: accumulates computed tendencies
      _RL gTr(nGud)
CEOP

#ifdef ALLOW_GUD
#ifdef GUD_GRAZING_MONOD_STYLE

c !LOCAL VARIABLES: ====================================================
      INTEGER jz, jp

      _RL Qc  (nplank)
      _RL Qn  (nplank)
      _RL Qp  (nplank)
      _RL Qsi (nplank)
      _RL Qfe (nplank)
      _RL QChl(nChl)
      _RL X   (nplank)
      _RL Xi  (nplank)
      _RL preygraz (nplank)
      _RL predgrazc(nplank)
      _RL sumprey, sumpref, grazphy, expFrac
      _RL totkillc, totkilln, totkillp, totkillsi, totkillfe
#ifdef ALLOW_CARBON
      _RL totkillPIC
#endif
      _RL graz2OC, graz2ON, graz2OP, graz2OFe
      _RL graz2POC, graz2PON, graz2POP, graz2POSi, graz2POFe
      _RL graz2PIC

c exudation
      _RL Xe
      _RL mortX

      _RL exude_DOC
      _RL exude_DON
      _RL exude_DOP
      _RL exude_DOFe

      _RL exude_PIC
      _RL exude_POC
      _RL exude_PON
      _RL exude_POP
      _RL exude_POSi
      _RL exude_POFe

      _RL mort_c(nplank)

      _RL respir
      _RL respir_c

#ifdef ALLOW_CDOM
      _RL graz2CDOM, exude_CDOM
#endif

C==== make all bio fields non-negative and compute quotas ==============

C fixed carbon quota, for now 1.0 (may change later)
      Qc(:) = 1.0
      X(:) = MAX(0., Ptr(ic:ec))/Qc(:)
      WHERE (X(:) .GT. 0.0 _d 0)
        Xi(:) = 1.0 _d 0/X(:)
      ELSEWHERE
        Xi(:) = 0.
      ENDWHERE
      Qn(:) = R_NC(:)
      Qp(:) = R_PC(:)
      Qsi(:) = R_SiC(:)
      Qfe(:) = R_FeC(:)
#ifdef CHLQUOTA
      QChl(:) = MAX(0., Ptr(ichl:echl))*Xi(1:nChl)
#endif

      preygraz(:)  = 0.0
      predgrazc(:) = 0.0
      graz2POC  = 0.0
      graz2PON  = 0.0
      graz2POP  = 0.0
      graz2POSI = 0.0
      graz2POFE = 0.0
      graz2OC   = 0.0
      graz2ON   = 0.0
      graz2OP   = 0.0
      graz2OFE  = 0.0
      graz2PIC  = 0.0

      DO jz = iMinPred, iMaxPred
        sumprey = 0.0
        DO jp = iMinPrey, iMaxPrey
         sumprey = sumprey + palat(jp,jz)*X(jp)
        ENDDO
        sumprey = MAX(0.0, sumprey - phygrazmin)
        sumpref = sumprey + phygrazmin
        totkillc = grazemax(jz)*grazTempFunc(jz)*X(jz)*
     &      (sumprey**hollexp/(sumprey**hollexp+kgrazesat(jz)**hollexp))

c        totkillc  = 0.0 _d 0
        totkilln  = 0.0 _d 0
        totkillp  = 0.0 _d 0
        totkillsi = 0.0 _d 0
        totkillfe = 0.0 _d 0
#ifdef ALLOW_CARBON
        totkillPIC = 0.0 _d 0
#endif
        DO jp = iMinPrey, iMaxPrey
          grazphy = totkillc*palat(jp,jz)*X(jp)/sumpref

          preygraz(jp) = preygraz(jp) + grazphy

c          totkillc  = totkillc  + grazphy
          totkilln  = totkilln  + grazphy*Qn (jp)
          totkillp  = totkillp  + grazphy*Qp (jp)
          totkillsi = totkillsi + grazphy*Qsi(jp)
          totkillfe = totkillfe + grazphy*Qfe(jp)
#ifdef ALLOW_CARBON
          totkillPIC = totkillPIC + grazphy*R_PICPOC(jp)
#endif

          predgrazc(jz) = predgrazc(jz) + grazphy*asseff(jp,jz)
        ENDDO

        expFrac = ExportFracGraz(jz)
        graz2POC =graz2POC  + expFrac*(totkillc -predgrazc(jz))
        graz2PON =graz2PON  + expFrac*(totkilln -predgrazc(jz)*Qn (jz))
        graz2POP =graz2POP  + expFrac*(totkillp -predgrazc(jz)*Qp (jz))
        graz2POFe=graz2POFe + expFrac*(totkillfe-predgrazc(jz)*Qfe(jz))
        graz2POSi=graz2POSi +          totkillsi-predgrazc(jz)*Qsi(jz)
#ifdef ALLOW_CARBON
        graz2PIC =graz2PIC  + totkillPIC
#endif

        graz2OC  = graz2OC  + totkillc  - predgrazc(jz)
        graz2ON  = graz2ON  + totkilln  - predgrazc(jz)*Qn (jz)
        graz2OP  = graz2OP  + totkillp  - predgrazc(jz)*Qp (jz)
        graz2OFe = graz2OFe + totkillfe - predgrazc(jz)*Qfe(jz)
      ENDDO

C==== tendencies =======================================================

      gTr(iDOC )=gTr(iDOC ) + graz2OC  - graz2POC
      gTr(iDON )=gTr(iDON ) + graz2ON  - graz2PON
      gTr(iDOP )=gTr(iDOP ) + graz2OP  - graz2POP
      gTr(iDOFe)=gTr(iDOFe) + graz2OFe - graz2POFe
      gTr(iPOC )=gTr(iPOC ) + graz2POC
      gTr(iPON )=gTr(iPON ) + graz2PON
      gTr(iPOP )=gTr(iPOP ) + graz2POP
      gTr(iPOSi)=gTr(iPOSi) + graz2POSi
      gTr(iPOFe)=gTr(iPOFe) + graz2POFe
#ifdef ALLOW_CARBON
      gTr(iPIC )=gTr(iPIC ) + graz2PIC
#endif
#ifdef ALLOW_CDOM
      graz2CDOM = fracCDOM*(graz2OP - graz2POP)
      gTr(iCDOM)=gTr(iCDOM) + graz2CDOM
      gTr(iDOC )=gTr(iDOC )             - R_CP_CDOM*graz2CDOM
      gTr(iDON )=gTr(iDON )             - R_NP_CDOM*graz2CDOM
      gTr(iDOP )=gTr(iDOP ) - graz2CDOM
      gTr(iDOFe)=gTr(iDOFe)             - R_FeP_CDOM*graz2CDOM
#endif
      DO jp = iMinPrey, iMaxPrey
       gTr(ic+jp-1)=gTr(ic+jp-1)-preygraz(jp)
#ifdef CHLQUOTA
       gTr(iChl+jp-1)=gTr(iChl+jp-1) - preygraz(jp)*QChl(jp)
#endif
      ENDDO
      DO jz = iMinPred, iMaxPred
       gTr(ic+jz-1)=gTr(ic+jz-1)+predgrazc(jz)
      ENDDO


C==== mortality ========================================================
      exude_DOC  = 0.0 _d 0
      exude_POC  = 0.0 _d 0
      exude_DON  = 0.0 _d 0
      exude_PON  = 0.0 _d 0
      exude_DOFe = 0.0 _d 0
      exude_POFe = 0.0 _d 0
      exude_DOP  = 0.0 _d 0
      exude_POP  = 0.0 _d 0
      exude_POSi = 0.0 _d 0
      exude_PIC  = 0.0 _d 0
      respir     = 0.0 _d 0

      DO jp = 1, nplank
        Xe = MAX(0.0, X(jp) - Xmin(jp))
        mortX = mort(jp)*Xe*MAX(mortTempFuncMin(jp), mortTempFunc) +
     &          mort2(jp)*Xe*Xe*MAX(mort2TempFuncMin(jp), mort2TempFunc)

        mort_c(jp) = mortX

        exude_DOC = exude_DOC + (1.-ExportFrac(jp))*mortX
        exude_POC = exude_POC +     ExportFrac(jp) *mortX

        exude_DON = exude_DON + (1.-ExportFrac(jp))*mortX*Qn(jp)
        exude_PON = exude_PON +     ExportFrac(jp) *mortX*Qn(jp)

        exude_DOP = exude_DOP + (1.-ExportFrac(jp))*mortX*Qp(jp)
        exude_POP = exude_POP +     ExportFrac(jp) *mortX*Qp(jp)

        exude_DOFe = exude_DOFe + (1.-ExportFrac(jp))*mortX*Qfe(jp)
        exude_POFe = exude_POFe +     ExportFrac(jp) *mortX*Qfe(jp)

        exude_POSi = exude_POSi + mortX*Qsi(jp)

        exude_PIC = exude_PIC + mortX*R_PICPOC(jp)

        respir_c = respiration(jp)*Xe*reminTempFunc
        respir = respir + respir_c

        gTr(ic+jp-1)=gTr(ic+jp-1)  -  mortX - respir_c
#ifdef NQUOTA
        gTr(in+jp-1)=gTr(in+jp-1)  -  mortX*Qn(jp)
#endif
#ifdef PQUOTA
        gTr(ip+jp-1)=gTr(ip+jp-1)  -  mortX*Qp(jp)
#endif
#ifdef SIQUOTA
        gTr(isi+jp-1)=gTr(isi+jp-1) - mortX*Qsi(jp)
#endif
#ifdef FEQUOTA
        gTr(ife+jp-1)=gTr(ife+jp-1) - mortX*Qfe(jp)
#endif

#ifdef ALLOW_EXUDE
        exude_DOC = exude_DOC + (1.-ExportFrac(jp))*kexcC(jp)*Xe
        exude_POC = exude_POC +     ExportFrac(jp) *kexcC(jp)*Xe
        exude_DON = exude_DON + (1.-ExportFrac(jp))*kexcN(jp)*Xe*Qn(jp)
        exude_PON = exude_PON +     ExportFrac(jp) *kexcN(jp)*Xe*Qn(jp)
        exude_DOP = exude_DOP + (1.-ExportFrac(jp))*kexcP(jp)*Xe*Qp(jp)
        exude_POP = exude_POP +     ExportFrac(jp) *kexcP(jp)*Xe*Qp(jp)
        exude_DOFe=exude_DOFe+(1.-ExportFrac(jp))*kexcFe(jp)*Xe*Qfe(jp)
        exude_POFe=exude_POFe+    ExportFrac(jp) *kexcFe(jp)*Xe*Qfe(jp)
        exude_POSi = exude_POSi + kexcSi(jp)*Xe*Qsi(jp)
        gTr(ic+jp-1)=gTr(ic+jp-1)   - kexcC(jp)*Xe
#ifdef NQUOTA
        gTr(in+jp-1)=gTr(in+jp-1)   - kexcN(jp)*Xe*Qn(jp)
#endif
#ifdef PQUOTA
        gTr(ip+jp-1)=gTr(ip+jp-1)   - kexcP(jp)*Xe*Qp(jp)
#endif
#ifdef SIQUOTA
        gTr(isi+jp-1)=gTr(isi+jp-1) - kexcSi(jp)*Xe*Qsi(jp)
#endif
#ifdef FEQUOTA
        gTr(ife+jp-1)=gTr(ife+jp-1) - kexcFe(jp)*Xe*Qfe(jp)
#endif
#endif
      ENDDO

#ifdef CHLQUOTA
      DO jp = 1, nChl
        gTr(iChl+jp-1)=gTr(iChl+jp-1) - mort_c(jp)*QChl(jp)
      ENDDO
#endif

      gTr(iDIC )=gTr(iDIC ) + respir

      gTr(iDOC )=gTr(iDOC ) + exude_DOC
      gTr(iDON )=gTr(iDON ) + exude_DON
      gTr(iDOP )=gTr(iDOP ) + exude_DOP
      gTr(iDOFe)=gTr(iDOFe) + exude_DOFe

      gTr(iPIC )=gTr(iPIC ) + exude_PIC
      gTr(iPOC )=gTr(iPOC ) + exude_POC
      gTr(iPON )=gTr(iPON ) + exude_PON
      gTr(iPOP )=gTr(iPOP ) + exude_POP
      gTr(iPOSi)=gTr(iPOSi) + exude_POSi
      gTr(iPOFe)=gTr(iPOFe) + exude_POFe
#ifdef ALLOW_CDOM
      exude_CDOM = fracCDOM*exude_DOP
      gTr(iCDOM)=gTr(iCDOM) + exude_CDOM
      gTr(iDOC )=gTr(iDOC )              - R_CP_CDOM*exude_CDOM
      gTr(iDON )=gTr(iDON )              - R_NP_CDOM*exude_CDOM
      gTr(iDOP )=gTr(iDOP ) - exude_CDOM
      gTr(iDOFe)=gTr(iDOFe)              - R_FeP_CDOM*exude_CDOM
#endif

#endif /* GUD_GRAZING_MONOD_STYLE */
#endif /* ALLOW_GUD */

      RETURN
      END SUBROUTINE

