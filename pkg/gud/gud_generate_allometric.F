C $Header$
C $Name$

#include "GUD_OPTIONS.h"

CBOP
C !ROUTINE: GUD_GENERATE_ALLOMETRIC

C !INTERFACE: ==========================================================
      SUBROUTINE GUD_GENERATE_ALLOMETRIC( myThid )

C !DESCRIPTION:
C     Generate parameters for plankton types using a "Monte Carlo"
C     approach using an allometric approach.

C !USES: ===============================================================
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "GUD_SIZE.h"
#include "GUD_PARAMS.h"
#include "GUD_GENPARAMS.h"
#include "GUD_TRAITPARAMS.h"
#include "GUD_TRAITS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP

#ifdef ALLOW_GUD

C !FUNCTIONS: ==========================================================
      _RL GUD_RANDOM
      EXTERNAL GUD_RANDOM
      _RL GUD_RANDOM_NORMAL
      EXTERNAL GUD_RANDOM_NORMAL

C !LOCAL VARIABLES: ====================================================
C     msgBuf   :: Informational/error meesage buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER jp,jp2,l
      _RL RandNo
      _RL mortdays
      _RL year
      _RL month
      _RL fiveday
      _RL rtime
      _RL tmpsrt
      _RL tmpend
      _RL tmprng
      _RL prd_pry
      _RL factor
      _RL taxon_mu(nplank)
      _RL a,b,p,error
      _RL heterotrophy(nplank)
#ifdef ALLOW_RADTRANS
      INTEGER iopt
      _RL sm
#endif

C     these used to be global in quota pkg
      INTEGER use_NO3(nplank)
      INTEGER use_Si(nplank)
      INTEGER pft(nplank)
      _RL autotrophy(nplank)
      _RL biovol(nplank)
      _RL qcarbon(nplank)
      _RL pp_opt(nplank)
      _RL pp_sig(nplank)
      _RL pp_pry(nplank)

C ======================================================================

      IF (nPro+nSyn+nEuk+nDiatom .GE. nPlank .OR.
     &    ngrazer .GT. nPlank) THEN
        STOP 'Need nPlank >= ngrazer and > nPro+nSyn+nEuk+nDiatom'
      ENDIF

c..........................................................
c Generate plankton volumes and stochastic parameters
c..........................................................
      factor = 2. _d 0
c       Allocate Phytoplankton Taxa
c
c       pft != 6 are grazed
c       heterotrophy == 1 graze
c Prochloro
      do jp2=1,nPro
        jp = jp2
        biovol(jp)    = biovol0(1) * factor**(jp2-1)
        autotrophy(jp)= 1.00 _d 0
        use_NO3(jp)   = 1
        nsource(jp)   = 3
        diacoc(jp)    = 0
        taxon_mu(jp)  = 1.00 _d 0
        pft(jp)       = 1
      enddo
c Synnecho
      do jp2=1,nSyn
        jp = nPro + jp2
        biovol(jp)    = biovol0(2) * factor**(jp2-1)
        autotrophy(jp)= 1.00 _d 0
        use_NO3(jp)   = 1
        nsource(jp)   = 3
        diacoc(jp)    = 0
        taxon_mu(jp)  = 1.40 _d 0
        pft(jp)       = 2
      enddo
c Small Euk
      do jp2=1,nEuk
        jp = nPro + nSyn + jp2
        biovol(jp)    = biovol0(3) * factor**(jp2-1)
        autotrophy(jp)= 1.0 _d 0
        use_NO3(jp)   = 1
        nsource(jp)   = 3
        diacoc(jp)    = 0
        taxon_mu(jp)  = 2.10 _d 0
        pft(jp)       = 3
      enddo
c Diatoms
      do jp2=1,nDiatom
        jp = nPro + nSyn + nEuk + jp2
        biovol(jp)    = biovol0(4) * factor**(jp2-1)
        autotrophy(jp)= 1.0 _d 0
        use_NO3(jp)   = 1
        nsource(jp)   = 3
        diacoc(jp)    = 0
        taxon_mu(jp)  = 3.8 _d 0
        pft(jp)       = 4
      enddo
c Specialist grazers
      do jp2=1,ngrazer
        jp = nPlank - ngrazer + jp2
        biovol(jp)    = biovol0(5) * factor**(jp2-1)
        autotrophy(jp)= 0.00 _d 0
        use_NO3(jp)   = 0
        nsource(jp)   = 0
        diacoc(jp)    = 0
        taxon_mu(jp)  = 0.00 _d 0
        pft(jp)       = 6
      enddo
c
      do jp=1,nPlank
       heterotrophy(jp)=1.0 _d 0 - autotrophy(jp)
       diazotroph(jp) = 0
       if (diacoc(jp) .eq. 1) then
         use_si(jp) = 0
       endif
c      elemental ratios for the non-quota elements
       IF(pft(jp).EQ.6) THEN
        R_NC(jp) = val_R_NC_zoo
        R_PC(jp) = val_R_PC_zoo
        R_SiC(jp) = val_R_SiC_zoo
        R_FeC(jp) = val_R_FeC_zoo
        R_ChlC(jp) = val_R_ChlC_zoo
        R_PICPOC(jp) = val_R_PICPOC_zoo
       ELSE
        IF(diacoc(jp) .EQ. 1.0 _d 0)then
          R_SiC(jp) =  val_R_SiC_diatom
        ELSE
          R_SiC(jp) = 0.0 _d 0
        ENDIF
        IF(diacoc(jp) .EQ. 2.0 _d 0)then
          R_PICPOC(jp) =  val_R_PICPOC
        ELSE
          R_PICPOC(jp) = 0.0 _d 0
        ENDIF
c   set N:P and iron requirement according to diazotroph status
        IF(diazotroph(jp) .EQ. 1.0 _d 0)then
          R_NC(jp) = val_R_NC_diaz
          R_FeC(jp) =  val_R_FeC_diaz
        ELSE
          R_NC(jp) = val_R_NC
          R_FeC(jp) = val_R_FeC
        ENDIF
c   set C:P ratio
        R_PC(jp) = val_R_PC
        R_ChlC(jp) = val_R_ChlC
       ENDIF
      enddo
c
c ----------------------------------------------------------------------
c Allometry
#ifdef UNCERTAINTY
      error = 1.0 _d 0
#else
      error = 0.0 _d 0
C       set stdev of allometric parameters to zero
#endif
c ----------------------------------------------------------------------
      do jp=1,nplank
        ! parameters independent of nutrient element

        IF (autotrophy(jp) .EQ. 1) THEN
          Xmin(jp) = phymin
          amminhib(jp) = val_amminhib
          acclimtimescl(jp) = val_acclimtimescl
        ELSE
          Xmin(jp) = 0.0 _d 0
        ENDIF

c CARBON CONTENT
        p = gud_random(myThid)
        call invnormal(a,p,
     &         log10(a_qcarbon),log10(ae_qcarbon)*error)
        call invnormal(b,p,b_qcarbon,be_qcarbon*error)
        qcarbon(jp)    = 10. _d 0**a * biovol(jp) ** b
c INITIAL SLOPE P-I
        p = gud_random(myThid)
        call invnormal(a,p,
     &         log10(a_alphachl),log10(ae_alphachl)*error)
        call invnormal(b,p,b_alphachl,be_alphachl*error)
        do l=1,nlam
          alphachl(jp,l) = 10. _d 0**a * biovol(jp) ** b
        enddo
c RESPIRATION RATE
        p = gud_random(myThid)
        IF (a_respir.NE.0. _d 0) THEN
        call invnormal(a,p,
     &         log10(a_respir),log10(ae_respir)*error)
        call invnormal(b,p,b_respir,be_respir*error)
        respiration(jp) = 10. _d 0**a
     &                    * (12. _d 9 * qcarbon(jp)) ** b
     &                    / qcarbon(jp)
        if (pft(jp).eq.6) then
          respiration(jp) = respiration(jp) * 0.50 _d 0
        endif
        ELSE
          respiration(jp) = 0.0 _d 0
        ENDIF
c MAXIMUM GRAZING RATE
        p = gud_random(myThid)
        call invnormal(a,p,
     &         log10(a_graz),log10(ae_graz)*error)
        call invnormal(b,p,b_graz,be_graz*error)
        grazemax(jp) = 10. _d 0**a * biovol(jp) ** b
     &                   * heterotrophy(jp)
c GRAZING SIZE PREFERENCE RATIO
        if (pft(jp).eq.5) then ! dinoflagellates prey upon similar sized plankton
          pp_opt(jp)      = 1.0 _d 0
        else ! other types follow standard relationship
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_prdpry),log10(ae_prdpry)*error)
          call invnormal(b,p,b_prdpry,be_prdpry*error)
          pp_opt(jp)      = 10. _d 0**a * biovol(jp) ** b
        endif
c STANDARD DEVIATION OF SIZE PREFERENCE
c          pp_sig(jp) = 0.50 _d 0
c          pp_sig(jp) = 1.00 _d 0
          pp_sig(jp) = pp_sig_val
c FRACTION GRAZED TO POM
C            p = gud_random(myThid)
C            call invnormal(a,p,
C     &           log10(a_beta_graz(io)),log10(ae_beta_graz(io))*error)
C            call invnormal(b,p,b_beta_graz(io),be_beta_graz(io)*error)
C            beta_graz(io,jp) = 10. _d 0**a * biovol(jp) ** b
C            beta_graz(io,jp)=min(max(
C     &                       beta_graz(io,jp),0.5 _d 0),0.9 _d 0)
        if (pft(jp).lt.3) ExportFracGraz(jp)=1.0-0.8
        if (pft(jp).gt.2) ExportFracGraz(jp)=1.0-0.5
c FRACTION MORTALITY TO POM
C            p = gud_random(myThid)
C            call invnormal(a,p,
C     &           log10(a_beta_mort(io)),log10(ae_beta_mort(io))*error)
C            call invnormal(b,p,b_beta_mort(io),be_beta_mort(io)*error)
C            beta_mort(io,jp)= 10. _d 0**a * biovol(jp) ** b
C            beta_mort(io,jp)=min(max(
C     &                       beta_mort(io,jp),0.5 _d 0),0.9 _d 0)
        if (pft(jp).lt.3) ExportFrac(jp)=1.0-0.8
        if (pft(jp).gt.2) ExportFrac(jp)=1.0-0.5
c GRAZING HALF-SATURATION
        p = gud_random(myThid)
        call invnormal(a,p,
     &         log10(a_kg),log10(ae_kg)*error)
        call invnormal(b,p,b_kg,be_kg*error)
        kgrazesat(jp) = 10. _d 0**a * biovol(jp) ** b
c PHYTOPLANKTON SINKING
        p = gud_random(myThid)
        call invnormal(a,p,
     &         log10(a_biosink),log10(ae_biosink)*error)
        call invnormal(b,p,b_biosink,be_biosink*error)
        if (pft(jp).eq.6) then
          wsink(jp)    =  0.0 _d 0 ! grazers don't sink
        else
          wsink(jp)    =  (10.0 _d 0**a) * biovol(jp) ** b
        endif
c SWIMMING
        p = gud_random(myThid)
        IF (a_bioswim.NE.0. _d 0) THEN
        call invnormal(a,p,
     &         log10(a_bioswim),log10(ae_bioswim)*error)
        call invnormal(b,p,b_bioswim,be_bioswim*error)
        if (autotrophy(jp).eq.1) then
          wswim(jp)    =  0.00 _d 0 ! only grazers can swim
        else
          wswim(jp)    = (10.0 _d 0**a) * biovol(jp) ** b
        endif
        ELSE
          wswim(jp)    =  0.00 _d 0
        ENDIF
c MORTALITY
        ! constant background mortality
        mort(jp) = a_mort
        if (pft(jp).eq.6) then
          ! grazers have lower mortality
          mort(jp) = 1.00 _d 0 * a_mort
        endif
        mort2(jp) = 0.0 _d 0
        if (pft(jp).eq.6) then
          mortTempFuncMin(jp) = 0.0 _d 0
          mort2TempFuncMin(jp) = 0.0 _d 0
        else
          mortTempFuncMin(jp) = 1.0 _d 0
          mort2TempFuncMin(jp) = 1.0 _d 0
        endif
C parameters relating to inorganic nutrients
c MAXIMUM NUTRIENT UPTAKE RATE
C carbon
        p = gud_random(myThid)
        call invnormal(a,p,
     &         log10(a_vmax_DIC),log10(ae_vmax_DIC)*error)
        call invnormal(b,p,b_vmax_DIC,be_vmax_DIC*error)
        PCmax(jp)= 10. _d 0**a * biovol(jp) ** b
     &             * taxon_mu(jp)

        IF (jp .LE. nNQuota) THEN
C NH4
          p = gud_random(myThid)
          call invnormal(a,p,log10(a_vmax_NH4),log10(ae_vmax_NH4)*error)
          call invnormal(b,p,b_vmax_NH4,be_vmax_NH4*error)
C         Pro has higher NH4 affinity
          vmax_NH4(jp)= 10. _d 0**a * biovol(jp) ** b * autotrophy(jp)*
     &                  (2.0 - use_NO3(jp))
          p = gud_random(myThid)
          call invnormal(a,p,log10(a_kn_NH4),log10(ae_kn_NH4)*error)
          call invnormal(b,p,b_kn_NH4,be_kn_NH4*error)
          ksatNH4(jp) = 10. _d 0**a * biovol(jp) ** b
C     *                  * autotrophy(jp)
C NO2     
          p = gud_random(myThid)
          call invnormal(a,p,log10(a_vmax_NO2),log10(ae_vmax_NO2)*error)
          call invnormal(b,p,b_vmax_NO2,be_vmax_NO2*error)
          vmax_NO2(jp)= 10. _d 0**a * biovol(jp) ** b * autotrophy(jp)
          p = gud_random(myThid)
          call invnormal(a,p,log10(a_kn_NO2),log10(ae_kn_NO2)*error)
          call invnormal(b,p,b_kn_NO2,be_kn_NO2*error)
          ksatNO2(jp) = 10. _d 0**a * biovol(jp) ** b
C     *                  * autotrophy(jp)
C NO3     
          p = gud_random(myThid)
          call invnormal(a,p,log10(a_vmax_NO3),log10(ae_vmax_NO3)*error)
          call invnormal(b,p,b_vmax_NO3,be_vmax_NO3*error)
          vmax_NO3(jp)= 10. _d 0**a * biovol(jp) ** b * autotrophy(jp)*
     &                  use_NO3(jp)
          p = gud_random(myThid)
          call invnormal(a,p,log10(a_kn_NO3),log10(ae_kn_NO3)*error)
          call invnormal(b,p,b_kn_NO3,be_kn_NO3*error)
          ksatNO3(jp) = 10. _d 0**a * biovol(jp) ** b
C     *                  * autotrophy(jp)
        ENDIF
        IF (jp .LE. nPQuota) THEN
C PO4     
          p = gud_random(myThid)
          call invnormal(a,p,log10(a_vmax_PO4),log10(ae_vmax_PO4)*error)
          call invnormal(b,p,b_vmax_PO4,be_vmax_PO4*error)
          vmax_PO4(jp)= 10. _d 0**a * biovol(jp) ** b * autotrophy(jp)
          p = gud_random(myThid)
          call invnormal(a,p,log10(a_kn_PO4),log10(ae_kn_PO4)*error)
          call invnormal(b,p,b_kn_PO4,be_kn_PO4*error)
          ksatPO4(jp) = 10. _d 0**a * biovol(jp) ** b
C     *                  * autotrophy(jp)
        ENDIF
        IF (jp .LE. nSiQuota) THEN
C SiO2    
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_vmax_SiO2),log10(ae_vmax_SiO2)*error)
          call invnormal(b,p,b_vmax_SiO2,be_vmax_SiO2*error)
          vmax_SiO2(jp)= 10. _d 0**a * biovol(jp) ** b * autotrophy(jp)*
     &                   float(use_Si(jp))
          p = gud_random(myThid)
          call invnormal(a,p,log10(a_kn_SiO2),log10(ae_kn_SiO2)*error)
          call invnormal(b,p,b_kn_SiO2,be_kn_SiO2*error)
          ksatSiO2(jp) = 10. _d 0**a * biovol(jp) ** b*
     &                   float(use_Si(jp))
C     *                  * autotrophy(jp)
C Silicate parameters to zero for non-diatoms
        ENDIF
        IF (jp .LE. nFeQuota) THEN
C FeT     
          p = gud_random(myThid)
          call invnormal(a,p,log10(a_vmax_FeT),log10(ae_vmax_FeT)*error)
          call invnormal(b,p,b_vmax_FeT,be_vmax_FeT*error)
          vmax_FeT(jp)= 10. _d 0**a * biovol(jp) ** b * autotrophy(jp)
          p = gud_random(myThid)
          call invnormal(a,p,log10(a_kn_FeT),log10(ae_kn_FeT)*error)
          call invnormal(b,p,b_kn_FeT,be_kn_FeT*error)
          ksatFeT(jp) = 10. _d 0**a * biovol(jp) ** b
C     *                  * autotrophy(jp)
        ENDIF
c
C parameters relating to quota nutrients
c EXCRETION
        IF (a_kexc_c.NE.0. _d 0 .AND. ae_kexc_c.NE.0. _d 0) THEN
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_kexc_c),log10(ae_kexc_c)*error)
          call invnormal(b,p,b_kexc_c,be_kexc_c*error)
          kexcC(jp) = 10. _d 0**a * biovol(jp) ** b
        ENDIF

        IF (jp .LE. nNQuota .AND. a_kexc_n.NE.0. _d 0 .AND.
     &                           ae_kexc_n.NE.0. _d 0) THEN
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_kexc_n),log10(ae_kexc_n)*error)
          call invnormal(b,p,b_kexc_n,be_kexc_n*error)
          kexcN(jp) = 10. _d 0**a * biovol(jp) ** b
        ENDIF

        IF (jp .LE. nPQuota .AND. a_kexc_p.NE.0. _d 0 .AND.
     &                           ae_kexc_p.NE.0. _d 0) THEN
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_kexc_p),log10(ae_kexc_p)*error)
          call invnormal(b,p,b_kexc_p,be_kexc_p*error)
          kexcP(jp) = 10. _d 0**a * biovol(jp) ** b
        ENDIF

c        p = gud_random(myThid)
c        call invnormal(a,p,log10(a_kexc_si),log10(ae_kexc_si)*error)
c        call invnormal(b,p,b_kexc_si,be_kexc_si*error)
c        kexcSi(jp) = 10. _d 0**a * biovol(jp) ** b
        kexcSi(jp) = 0.0 _d 0

c        p = gud_random(myThid)
c        call invnormal(a,p,log10(a_kexc_fe),log10(ae_kexc_fe)*error)
c        call invnormal(b,p,b_kexc_fe,be_kexc_fe*error)
c        kexcFe(jp) = 10. _d 0**a * biovol(jp) ** b
        kexcFe(jp) = 0.0 _d 0

        IF (jp .LE. nNQuota) THEN
c MINIMUM QUOTA
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_qmin_n),log10(ae_qmin_n)*error)
          call invnormal(b,p,b_qmin_n,be_qmin_n*error)
          Qnmin(jp) = 10. _d 0**a * biovol(jp) ** b
c MAXIMUM QUOTA
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_qmax_n),log10(ae_qmax_n)*error)
          call invnormal(b,p,b_qmax_n,be_qmax_n*error)
          Qnmax(jp) = 10. _d 0**a * biovol(jp) ** b
        ENDIF

        IF (jp .LE. nPQuota) THEN
c MINIMUM QUOTA
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_qmin_p),log10(ae_qmin_p)*error)
          call invnormal(b,p,b_qmin_p,be_qmin_p*error)
          Qpmin(jp) = 10. _d 0**a * biovol(jp) ** b
c MAXIMUM QUOTA
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_qmax_p),log10(ae_qmax_p)*error)
          call invnormal(b,p,b_qmax_p,be_qmax_p*error)
          Qpmax(jp) = 10. _d 0**a * biovol(jp) ** b
        ENDIF

        IF (jp .LE. nsiQuota) THEN
c MINIMUM QUOTA
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_qmin_si),log10(ae_qmin_si)*error)
          call invnormal(b,p,b_qmin_si,be_qmin_si*error)
          Qsimin(jp) = 10. _d 0**a * biovol(jp) ** b *
     &                 *use_Si(jp)
c MAXIMUM QUOTA
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_qmax_si),log10(ae_qmax_si)*error)
          call invnormal(b,p,b_qmax_si,be_qmax_si*error)
          Qsimax(jp) = 10. _d 0**a * biovol(jp) ** b *
     &                 *use_Si(jp)
        ENDIF

        IF (jp .LE. nfeQuota) THEN
c MINIMUM QUOTA
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_qmin_fe),log10(ae_qmin_fe)*error)
          call invnormal(b,p,b_qmin_fe,be_qmin_fe*error)
          Qfemin(jp) = 10. _d 0**a * biovol(jp) ** b
c MAXIMUM QUOTA
          p = gud_random(myThid)
          call invnormal(a,p,
     &           log10(a_qmax_fe),log10(ae_qmax_fe)*error)
          call invnormal(b,p,b_qmax_fe,be_qmax_fe*error)
          Qfemax(jp) = 10. _d 0**a * biovol(jp) ** b
        ENDIF

c Zooplankton have approximately Redfieldian N:C ratio
        if (pft(jp).eq.6) then
          qnmin(jp) = 0.0755 _d 0
          qnmax(jp) = 0.1510 _d 0
        endif
c PREFERENCE FUNCTION
        ! assign grazing preference according to predator/prey radius ratio
        do jp2=1,nplank ! jp2 denotes prey
          if (heterotrophy(jp).gt.0. _d 0.and.pft(jp2).ne.6) then
            prd_pry = biovol(jp) / biovol(jp2)
            palat(jp2,jp) =
#ifdef GUD_GRAZING_MONOD_STYLE
     &        EXP(-(LOG(prd_pry/pp_opt(jp))**2) / (2*pp_sig(jp)**2))
     &                                          / pp_sig(jp)/2. _d 0
#else
            ! lower preference for larger P
     &        1.0 _d 0
C     &        biovol(jp2)**(-0.16 _d 0)
#endif
C              ! reduce diatom palatability
C              if (pft(jp2).eq.5) then
C                palat(jp2,jp) = palat(jp2,jp) * 0.8
C              endif
            if (palat(jp2,jp).lt.1. _d -4) then
              palat(jp2,jp)=0. _d 0
            endif
            asseff(jp2,jp) = val_ass_eff
          else
            palat(jp2,jp) = 0. _d 0
          endif
        enddo
c
c..........................................................
c generate phyto Temperature Function parameters
c.......................................................
        phytoTempCoeff(jp) = tempcoeff1
        phytoTempExp1(jp) = tempcoeff3
        phytoTempExp2(jp) = tempcoeff2_small
     &                      + (tempcoeff2_big-tempcoeff2_small)
     &                      * float(jp-1)/float(nplank-1)
        phytoTempOptimum(jp) = 2. _d 0
        phytoDecayPower(jp) = tempdecay

c..........................................................
      enddo

C ======================================================================
C dependent traits

      DO jp = 1, nPlank
C       should be local - leave unset
C        mQyield(jp)        = 0.0 _d 0
C        physize(jp)        = 0.0 _d 0
C        phyto_esd(jp)      = 0.0 _d 0
C        phyto_vol(jp)      = 0.0 _d 0

C for non-quota nitrogen
        IF (nsource(jp) .EQ. 3) THEN
         useNH4(jp) = 1.0
         useNO2(jp) = 1.0
         useNO3(jp) = 1.0
         combNO(jp) = 1.0
        ELSEIF (nsource(jp) .EQ. 2) THEN
         useNH4(jp) = 1.0
         useNO2(jp) = 0.0
         useNO3(jp) = 0.0
         combNO(jp) = 0.0
        ELSEIF (nsource(jp) .EQ. 1) THEN
         useNH4(jp) = 1.0
         useNO2(jp) = 1.0
         useNO3(jp) = 0.0
         combNO(jp) = 0.0
        ELSE
         useNH4(jp) = 0.0
         useNO2(jp) = 0.0
         useNO3(jp) = 0.0
         combNO(jp) = 0.0
        ENDIF
      ENDDO

      DO jp = iMinPhoto, iMaxPhoto
#ifdef GEIDER
C      should depend on size
       RandNo = gud_random(myThid)
       chl2cmax(jp) = Smallchl2cmax +RandNo*Smallchl2cmaxrange

c NB:  LL Pro should be inhibited
       inhibcoef_geid(jp) = 0. _d 0

#ifdef ALLOW_RADTRANS
       ap_type(jp) = 0

       iopt = ap_type(jp)
       IF (1 .LE. iopt .AND. iopt .LE. nOpt) THEN
        DO l = 1,nlam
         aphy_chl(jp,l) = aphy_chl_type(iopt,l)
         aphy_chl_ps(jp,l) = aphy_chl_ps_type(iopt,l)
         bphy_chl(jp,l) = bphy_chl_type(iopt,l)
         bbphy_chl(jp,l) = bbphy_chl_type(iopt,l)
        ENDDO
       ELSE
        WRITE(msgBuf,'(2A,2I4)') 'GUD_INIT_FIXED: ',
     &   "invalid optical phyto type:",jp,iopt
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R GUD_INIT_FIXED'
       ENDIF

C NB:  this overrides alphachl set from allometry above!!!
       DO l = 1,nlam
         alphachl(jp,l) = mQyield(jp) * aphy_chl_ps(jp,l)
       ENDDO

       sm = 0.0 _d 0
       DO l = 1, nlam
        sm = sm + wb_width(l)*alphachl(jp,l)
       ENDDO
       alpha_mean(jp) = sm/wb_totalWidth

       chl2cmin(jp)=chl2cmax(jp)/
     &           (1+(chl2cmax(jp)* alpha_mean(jp) *2000. _d 0)/
     &              (2*pcmax(jp)))
#else
       alpha_mean(jp) = alphachl(jp,1)
       chl2cmin(jp) = 0. _d 0
#endif /* ALLOW_RADTRANS */

#else /* not GEIDER */
C      should depend on size
       RandNo = gud_random_normal(myThid)
       ksatPAR(jp) = abs(smallksatpar+smallksatparstd*RandNo)
       RandNo = gud_random_normal(myThid)
       kinhPAR(jp) = abs(smallkinhPAR+smallkinhPARstd*RandNo)
#endif /* GEIDER */
C     jp
      ENDDO

#endif  /*ALLOW_GUD*/

      RETURN
      END
