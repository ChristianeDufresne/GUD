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
      INTEGER jp,jp2,jz,g,gz,l,jpg(nGroup),gmin
      _RL logvol, vol, volmin
      _RL RandNo
      _RL prd_pry
      _RL kappa
      _RL taxon_mu(nplank)
      _RL isphoto(nplank)
      _RL ispred(nplank)
      _RL isprey(nplank)
      _RL prey_ass_eff(nplank)
#ifdef ALLOW_RADTRANS
      INTEGER iopt
      _RL sm
#endif

C     these used to be global in quota pkg
      _RL pp_opt(nplank)
      _RL pp_sig(nplank)
      _RL pp_pry(nplank)

      INTEGER group(nplank)
      INTEGER iGroup(nplank)
      INTEGER diacoc(nplank)
      INTEGER nsource(nplank)
      _RL biovol(nplank)
      _RL qcarbon(nplank)

C ======================================================================

      DO jp=1,nplank
        IF (diazo(jp) .NE. UNINIT_I) THEN
         WRITE(msgBuf,'(2A)') 'GUD_GENERATE_ALLOMETRIC: ',
     &    'do no set diazo, set grp_diazo instead.'
         CALL PRINT_ERROR( msgBuf, myThid )
         STOP 'ABNORMAL END: S/R GUD_GENERATE_ALLOMETRIC'
        ENDIF
      ENDDO

c..........................................................
c Generate plankton volumes and stochastic parameters
c..........................................................
c       Allocate Phytoplankton Taxa
c
c       isprey == 1 are grazed
c       ispred == 1 graze
      DO g=1,nGroup
        jpg(g) = 0
        IF (logvol0ind(g) .GT. 0) THEN
          logvol = logvol0base + (logvol0ind(g)-1)*logvol0inc
          biovol0(g) = 10 _d 0 ** logvol
        ENDIF
      ENDDO
      IF (gud_sort_biovol) THEN
c      sort by volume
       DO jp=1,nplank
        volmin = 1 _d 38
        gmin = 0
C       first check phototrophs
        DO g=1,ngroup
         IF (grp_photo(g).NE.0) THEN
          vol = biovol0(g)*biovolfac(g)**jpg(g)
          IF (jpg(g).LT.grp_nplank(g) .AND. vol.LT.volmin) THEN
           gmin = g
           volmin = vol
          ENDIF
         ENDIF
        ENDDO
C       then pure heterotrophs
        IF (gmin.EQ.0) THEN
         DO g=1,ngroup
          IF (grp_photo(g).EQ.0) THEN
           vol = biovol0(g)*biovolfac(g)**jpg(g)
           IF (jpg(g).LT.grp_nplank(g) .AND. vol.LT.volmin) THEN
            gmin = g
            volmin = vol
           ENDIF
          ENDIF
         ENDDO
        ENDIF
        IF (gmin.EQ.0) THEN
         STOP 'gmin'
        ENDIF
        group(jp) = gmin
        biovol(jp) = volmin
        jpg(gmin) = jpg(gmin) + 1
        igroup(jp) = jpg(gmin)
       ENDDO
       DO g=1,ngroup
        IF (jpg(g).NE.grp_nplank(g)) THEN
         STOP 'grp_nplank'
        ENDIF
       ENDDO
      ELSE
C      sort by group
       jp = 1
       DO g=1,ngroup
        DO jp2=1,grp_nplank(g)
         IF (jp .GT. nPlank) THEN
          WRITE(msgBuf,'(2A)') 'GUD_GENERATE_ALLOMETRIC: ',
     &     'need SUM(grp_nplank) = nPlank, nPlank too small'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R GUD_GENERATE_ALLOMETRIC'
         ENDIF
         group(jp)  = g
         biovol(jp) = biovol0(g) * biovolfac(g)**(jp2-1)
         jp = jp + 1
        ENDDO
       ENDDO
       IF (jp .NE. nPlank + 1) THEN
        WRITE(msgBuf,'(2A,2I4)') 'GUD_GENERATE_ALLOMETRIC: ',
     &   'need SUM(grp_nplank) = nPlank, not ',jp-1,nPlank
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R GUD_GENERATE_ALLOMETRIC'
       ENDIF
      ENDIF

      DO jp=1,nplank
        g = group(jp)
        jp2 = igroup(jp)

        isphoto(jp)    = grp_photo(g)
        ispred(jp)     = grp_pred(g)
        isprey(jp)     = grp_prey(g)
        nsource(jp)    = grp_nsource(g)
        diacoc(jp)     = grp_diacoc(g)
        diazo(jp)      = grp_diazo(g)
        taxon_mu(jp)   = grp_taxon_mu(g)
        R_NC(jp)       = grp_R_NC(g)
        R_PC(jp)       = grp_R_PC(g)
        R_SiC(jp)      = grp_R_SiC(g)
        R_FeC(jp)      = grp_R_FeC(g)
        R_ChlC(jp)     = grp_R_ChlC(g)
        R_PICPOC(jp)   = grp_R_PICPOC(g)

        Xmin(jp) = phymin*(1.0 _d 0 - ispred(jp))
        amminhib(jp) = grp_amminhib(g)
        acclimtimescl(jp) = grp_acclimtimescl(g)

#ifdef GEIDER
        mQyield(jp) = grp_mQyield(g)
        chl2cmax(jp) = grp_chl2cmax(g)
#endif

        IF (nsource(jp) .EQ. 3) THEN
         useNH4(jp) = 1
         useNO2(jp) = 1
         useNO3(jp) = 1
         combNO(jp) = 1
        ELSEIF (nsource(jp) .EQ. 2) THEN
         useNH4(jp) = 1
         useNO2(jp) = 0
         useNO3(jp) = 0
         combNO(jp) = 0
        ELSEIF (nsource(jp) .EQ. 1) THEN
         useNH4(jp) = 1
         useNO2(jp) = 1
         useNO3(jp) = 0
         combNO(jp) = 0
        ELSE
         useNH4(jp) = 0
         useNO2(jp) = 0
         useNO3(jp) = 0
         combNO(jp) = 0
        ENDIF

c CARBON CONTENT
        qcarbon(jp) = a_qcarbon(g) * biovol(jp)**b_qcarbon(g)
c RESPIRATION RATE
        respiration(jp) = a_respir(g)
     &                    * (12. _d 9 * qcarbon(jp))**b_respir(g)
     &                    / qcarbon(jp)
c FRACTION MORTALITY TO POM
        ExportFrac(jp) = grp_ExportFrac(g)
c PHYTOPLANKTON SINKING
        wsink(jp) = a_biosink(g) * biovol(jp)**b_biosink(g)
c SWIMMING
        wswim(jp) = a_bioswim(g) * biovol(jp)**b_bioswim(g)
     &              *(1.0 _d 0 - ispred(jp))
c MORTALITY
        ! constant background mortality
        mort(jp) = grp_mort(g)
        mort2(jp) = grp_mort2(g)
        ! if 0, temperature dependent, if 1, not.
        mortTempFuncMin(jp) = MIN(1, 1 - grp_tempMort(g))
        mort2TempFuncMin(jp) = MIN(1, 1 - grp_tempMort2(g))
C parameters relating to inorganic nutrients
c MAXIMUM NUTRIENT UPTAKE RATE
C carbon
        PCmax(jp)= a_vmax_DIC(g) * biovol(jp)**b_vmax_DIC(g)
     &             * taxon_mu(jp)

        vmax_NH4(jp)  = a_vmax_NH4(g)  * biovol(jp)**b_vmax_NH4(g)
        vmax_NO2(jp)  = a_vmax_NO2(g)  * biovol(jp)**b_vmax_NO2(g)
        vmax_NO3(jp)  = a_vmax_NO3(g)  * biovol(jp)**b_vmax_NO3(g)
        vmax_N(jp)    = a_vmax_N(g)    * biovol(jp)**b_vmax_N(g)
        vmax_PO4(jp)  = a_vmax_PO4(g)  * biovol(jp)**b_vmax_PO4(g)
        vmax_SiO2(jp) = a_vmax_SiO2(g) * biovol(jp)**b_vmax_SiO2(g)
        vmax_FeT(jp)  = a_vmax_FeT(g)  * biovol(jp)**b_vmax_FeT(g)

        Qnmin(jp) = a_qmin_n(g) * biovol(jp)**b_qmin_n(g)
        Qnmax(jp) = a_qmax_n(g) * biovol(jp)**b_qmax_n(g)

        Qpmin(jp) = a_qmin_p(g) * biovol(jp)**b_qmin_p(g)
        Qpmax(jp) = a_qmax_p(g) * biovol(jp)**b_qmax_p(g)

        Qsimin(jp) = a_qmin_si(g) * biovol(jp)**b_qmin_si(g)
        Qsimax(jp) = a_qmax_si(g) * biovol(jp)**b_qmax_si(g)

        Qfemin(jp) = a_qmin_fe(g) * biovol(jp)**b_qmin_fe(g)
        Qfemax(jp) = a_qmax_fe(g) * biovol(jp)**b_qmax_fe(g)

        ksatNH4(jp) = a_kn_NH4(g) * biovol(jp)**b_kn_NH4(g)
        ksatNO2(jp) = a_kn_NO2(g) * biovol(jp)**b_kn_NO2(g)
        ksatNO3(jp) = a_kn_NO3(g) * biovol(jp)**b_kn_NO3(g)
        ksatPO4(jp) = a_kn_PO4(g) * biovol(jp)**b_kn_PO4(g)
        ksatSiO2(jp) = a_kn_SiO2(g) * biovol(jp)**b_kn_SiO2(g)
        ksatFeT(jp) = a_kn_FeT(g) * biovol(jp)**b_kn_FeT(g)

        IF (GUD_effective_ksat) THEN
c        compute effective half sat for uptake of non-quota elements
C we compute it for NO3 and scale for others
         kappa= (ksatNO3(jp)*PCmax(jp)*Qnmin(jp)*(Qnmax(jp)-Qnmin(jp)))/
     &          (vmax_NO3(jp)*Qnmax(jp) +
     &                       PCmax(jp)*Qnmin(jp)*(Qnmax(jp)-Qnmin(jp)))
#ifndef NQUOTA
         ksatNO3(jp) = kappa
         ksatNO2(jp) = kappa*grp_ksatNO2fac(g)
         ksatNH4(jp) = kappa*grp_ksatNH4fac(g)
#endif
#ifndef PQUOTA
         ksatPO4(jp) = kappa/R_NC(jp)*R_PC(jp)
#endif
#ifndef SIQUOTA
         ksatSiO2(jp) = kappa/R_NC(jp)*R_SiC(jp)
#endif
#ifndef FEQUOTA
         ksatFeT(jp) = kappa/R_NC(jp)*R_FeC(jp)
#endif
        ENDIF

        vmax_NH4(jp) = vmax_NH4(jp) * useNH4(jp)
        vmax_NO2(jp) = vmax_NO2(jp) * useNO2(jp)
        vmax_NO3(jp) = vmax_NO3(jp) * useNO3(jp)
C Silicate parameters to zero for non-diatoms
        IF (diacoc(jp) .NE. 1) THEN
          vmax_SiO2(jp) = 0.0 _d 0
          ksatSiO2(jp) = 0.0 _d 0
          R_SiC(jp) = 0.0 _d 0
        ENDIF

C parameters relating to quota nutrients
c EXCRETION
        kexcC(jp) = a_kexc_c(g) * biovol(jp)**b_kexc_c(g)
        kexcN(jp) = a_kexc_n(g) * biovol(jp)**b_kexc_n(g)
        kexcP(jp) = a_kexc_p(g) * biovol(jp)**b_kexc_p(g)
        kexcSi(jp) = a_kexc_si(g) * biovol(jp)**b_kexc_si(g)
        kexcFe(jp) = a_kexc_fe(g) * biovol(jp)**b_kexc_fe(g)

c..........................................................
c generate phyto Temperature Function parameters
c.......................................................
        phytoTempCoeff(jp)   = grp_tempcoeff1(g)
        phytoTempExp1(jp)    = grp_tempcoeff3(g)
        phytoTempExp2(jp)    = grp_tempcoeff2(g)
        phytoTempOptimum(jp) = grp_tempopt(g)
        phytoDecayPower(jp)  = grp_tempdecay(g)

c MAXIMUM GRAZING RATE
        grazemax(jp) = a_graz(g) * biovol(jp)**b_graz(g)
c GRAZING SIZE PREFERENCE RATIO
        pp_opt(jp) = a_prdpry(g) * biovol(jp)**b_prdpry(g)
c STANDARD DEVIATION OF SIZE PREFERENCE
        pp_sig(jp) = grp_pp_sig(g)
c GRAZING HALF-SATURATION
        kgrazesat(jp) = a_kg(g) * biovol(jp)**b_kg(g)

#ifdef GEIDER
        inhibcoef_geid(jp) = grp_inhibcoef_geid(g)
#else
        ksatPAR(jp) = grp_ksatPAR(g)
        kinhPAR(jp) = grp_kinhPAR(g)
#endif /* GEIDER */

C     jp
      ENDDO

c..........................................................

c PREFERENCE FUNCTION
c     assign grazing preference according to predator/prey radius ratio
      DO jz=1,nplank
        gz = group(jz)
        DO jp=1,nplank ! jp denotes prey
          g = group(jp)
          IF (ispred(jz).GT.0 .AND. isprey(jp).GT.0) THEN
            prd_pry = biovol(jz) / biovol(jp)
            palat(jp,jz) =
#ifdef GUD_GRAZING_MONOD_STYLE
     &        EXP(-(LOG(prd_pry/pp_opt(jz))**2) / (2*pp_sig(jz)**2))
     &                                          / pp_sig(jz)/2. _d 0
#else
            ! lower preference for larger P
     &        1.0 _d 0
C     &        biovol(jp)**(-0.16 _d 0)
#endif
C              ! reduce diatom palatability
C              if (pft(jp).eq.5) then
C                palat(jp,jz) = palat(jp,jz) * 0.8
C              endif
            IF (palat(jp,jz).LT.palat_min) THEN
              palat(jp,jz) = 0. _d 0
            ENDIF
c ASSIMILATION EFFICIENCY FOR PREY jp
            asseff(jp,jz) = grp_ass_eff(g)
c FRACTION GRAZED TO POM
            ExportFracPreyPred(jp,jz) = grp_ExportFracPreyPred(g,gz)
          ELSE
            palat(jp,jz) = 0. _d 0
          ENDIF
        ENDDO
      ENDDO

#endif  /*ALLOW_GUD*/

      RETURN
      END
