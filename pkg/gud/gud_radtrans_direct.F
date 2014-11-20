C $Header$
C $Name$

#include "GUD_OPTIONS.h"

CBOP
C !ROUTINE: GUD_RADTRANS_DIRECT

C !INTERFACE: ==========================================================
      subroutine GUD_RADTRANS_DIRECT(
     I                   H,rmud,Edsf,Essf,a_k,bt_k,bb_k,kmax,
     O                   Edbot,Esbot,Eubot,Estop,Eutop,
     O                   PAR,
     O                   amp1, amp2,
     I                   myThid)

C !DESCRIPTION:
c
c  Model of irradiance in the water column.  Accounts for three
c  irradiance streams [Ackleson, Balch, Holligan, JGR, 1994],
c
c  Edbot = direct downwelling irradiance in W/m2 per waveband
c  Esbot = diffuse downwelling irradiance in W/m2 per waveband
c  Eubot = diffuse upwelling irradiance in W/m2 per waveband
c
c  all defined at the bottom of each layer.  Also computed are Estop,
c  Eutop at the top of each layer which should be very close to Esbot,
c  Eubot of the layer above.
c
c  The Ed equation is integrated exactly, Es and Eu are computed by
c  solving a set of linear equation for the amplitudes in the exact
c  solution [see, e.g., Kylling, Stamnes, Tsay, JAC, 1995].
c  The boundary condition in the deepest wet layer is
c  downward-decreasing modes only (i.e., zero irradiance at infinite
c  depth, assuming the optical properties of the last layer).
c
c  Also computed are scalar radiance and PAR at the grid cell center
c  (both in uEin/m2/s).
c
C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"                 /* Nr */
#include "EEPARAMS.h"
#include "GUD_SIZE.h"          /* nlam */
#ifdef ALLOW_RADTRANS
#include "GUD_PARAMS.h"        /* WtouEin */
#include "GUD_GENPARAMS.h"  /* gud_PAR_ilamLo/Hi
                                     gud_radmodThresh
                                     gud_rmus gud_rmuu */
#endif

C !INPUT PARAMETERS: ===================================================
C     H     :: layer thickness (including hFacC!)
C     rmud  :: inv.cosine of direct (underwater solar) zenith angle
C     Edsf  :: direct downwelling irradiance below surface per waveband
C     Essf  :: diffuse downwelling irradiance below surface per waveband
C     a_k   :: absorption coefficient per level and waveband (1/m)
C     bt_k  :: total scattering coefficient per level and waveband (1/m)
C              = forward + back scattering coefficient
C     bb_k  :: backscattering coefficient per level and waveband (1/m)
C     kmax  :: maximum number of layers to compute
      _RL H(Nr)
      _RL rmud
      _RL Edsf(nlam), Essf(nlam)
      _RL a_k(Nr,nlam), bt_k(Nr,nlam), bb_k(Nr,nlam)
      INTEGER kmax
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C     Edbot  :: direct downwelling irradiance at bottom of layer
C     Esbot  :: diffuse downwelling irradiance at bottom of layer
C     Eubot  :: diffuse upwelling irradiance at bottom of layer
C     Estop  :: diffuse downwelling irradiance at top of layer
C     Eutop  :: diffuse upwelling irradiance at top of layer
C     PAR :: total scalar irradiance at cell center per waveband
C     amp1   :: amplitude of downward increasing mode
C     amp2   :: amplitude of downward decreasing mode
      _RL Edbot(nlam,Nr),Esbot(nlam,Nr),Eubot(nlam,Nr)
      _RL Estop(nlam,Nr),Eutop(nlam,Nr)
      _RL PAR(Nr,nlam)
      _RL amp1(nlam,Nr), amp2(nlam,Nr)
CEOP

#ifdef ALLOW_RADTRANS

C !LOCAL VARIABLES: ====================================================
      INTEGER k, l, kbot
      _RL Edtop(nlam,Nr)
      _RL Etopwq, Ebotwq
      _RL zd
      _RL rmus,rmuu
      _RL cd,au,Bu,Cu
      _RL as,Bs,Cs,Bd,Fd
      _RL bquad,D
      _RL kappa1,kappa2,denom
      _RL c1,c2
      _RL r2(Nr),r1(Nr),x(Nr),y(Nr)
      _RL ed(Nr),e2(Nr),e1(Nr)
      _RL a3d(2*Nr), b3d(2*Nr), c3d(2*Nr), y3d(2*Nr)
      _RL rd, ru
      data rd /1.5 _d 0/   !these are taken from Ackleson, et al. 1994 (JGR)
      data ru /3.0 _d 0/

      rmus = gud_rmus
      rmuu = gud_rmuu

c     find deepest wet layer
      kbot = kmax
      DO WHILE (H(kbot).EQ.0 .AND. kbot.GT.0)
        kbot = kbot - 1
      ENDDO

      DO l = 1,nlam
       DO k=1,Nr
        Edtop(l,k) = 0.0
        Estop(l,k) = 0.0
        Eutop(l,k) = 0.0
        Edbot(l,k) = 0.0
        Esbot(l,k) = 0.0
        Eubot(l,k) = 0.0
        amp1(l,k) = 0.0
        amp2(l,k) = 0.0
       ENDDO
      ENDDO
      IF (kbot.GT.0) THEN
       DO l = 1,nlam
        IF (Edsf(l) .GE. gud_radmodThresh .OR.
     &      Essf(l) .GE. gud_radmodThresh) THEN
         DO k=1,kbot
          zd = H(k)
          cd = (a_k(k,l)+bt_k(k,l))*rmud
          au = a_k(k,l)*rmuu
          Bu = ru*bb_k(k,l)*rmuu
          Cu = au+Bu
          as = a_k(k,l)*rmus
          Bs = rd*bb_k(k,l)*rmus
          Cs = as+Bs
          Bd = bb_k(k,l)*rmud
          Fd = (bt_k(k,l)-bb_k(k,l))*rmud
          bquad = Cs + Cu
          D = 0.5*(bquad + SQRT(bquad*bquad - 4.0*Bs*Bu))
          kappa1 = D - Cs
          kappa2 = Cs - Bs*Bu/D  ! == D - Cu
          r1(k) = Bu/D
          r2(k) = Bs/D
          denom = (cd-Cs)*(cd+Cu) + Bs*Bu
          x(k) = -((cd+Cu)*Fd+Bu*Bd)/denom
          y(k) = (-Bs*Fd+(cd-Cs)*Bd)/denom
          ed(k) = EXP(-cd*zd)
          e1(k) = EXP(-kappa1*zd)
          e2(k) = EXP(-kappa2*zd)
         ENDDO

C integrate Ed equation first
         Edtop(l,1) = Edsf(l)
         DO k=1,kbot-1
          Edbot(l,k) = Edtop(l,k)*ed(k)
          Edtop(l,k+1) = Edbot(l,k)
         ENDDO
         Edbot(l,kbot) = Edtop(l,kbot)*ed(kbot)

C setup tridiagonal matrix of continuity/boundary conditions
C variables: c2(1), c1(1), c2(2), ..., c1(kbot)
C a3d,b3d,c3d: lower, main and upper diagonal
C y3d: right-hand side
C
C top b.c.: c2(1) + e1(1)*r1(1)*c1(1) = Essf - x(1)*Edsf
         a3d(1) = 0. _d 0  ! not used
         b3d(1) = 1.           ! A(1,1)*c2(1)
         c3d(1) = e1(1)*r1(1)  ! A(1,2)*c1(1)
         y3d(1) = Essf(l) - x(1)*Edsf(l)
C continuity at layer boundaries
         DO k=1, kbot-1
           a3d(2*k) = (1. - r2(k)*r1(k+1))*e2(k)  !   A(2k,2k-1)*c2(k)
           b3d(2*k) = r1(k) - r1(k+1)             ! + A(2k,2k  )*c1(k)
           c3d(2*k) = -1. + r2(k+1)*r1(k+1)       ! + A(2k,2k+1)*c2(k+1)
           y3d(2*k)= (x(k+1) - x(k) - r1(k+1)*(y(k+1)-y(k)))*Edbot(l,k)
           a3d(2*k+1) = 1 - r1(k)*r2(k)                !   A(2k+1,2k  )*c1(k)
           b3d(2*k+1) = r2(k) - r2(k+1)                ! + A(2k+1,2k+1)*c2(k+1)
           c3d(2*k+1) = (-1. + r1(k+1)*r2(k))*e1(k+1)  ! + A(2k+1,2k+2)*c1(k+1)
           y3d(2*k+1)= (y(k+1) - y(k) - r2(k)*(x(k+1)-x(k)))*Edbot(l,k)
         ENDDO
c bottom boundary condition: c1 = 0
         a3d(2*kbot) = 0. _d 0  !   A(2*kbot,2*kbot-1)*c2(kbot)
         b3d(2*kbot) = 1. _d 0  ! + A(2*kbot,2*kbot  )*c1(kbot)
         c3d(2*kbot) = 0. _d 0  ! not used
         y3d(2*kbot) = 0. _d 0  ! = 0

         CALL SOLVE_TRIDIAGONAL_PIVOT(a3d,b3d,c3d,y3d,2*kbot,myThid)

C compute irradiances
         DO k=1,kbot
          c2 = y3d(2*k-1)
          c1 = y3d(2*k)
          Estop(l,k) = c2 + r1(k)*e1(k)*c1 + x(k)*Edtop(l,k)
          Esbot(l,k) = e2(k)*c2 + r1(k)*c1 + x(k)*Edbot(l,k)
          Eutop(l,k) = r2(k)*c2 + e1(k)*c1 + y(k)*Edtop(l,k)
          Eubot(l,k) = r2(k)*e2(k)*c2 + c1 + y(k)*Edbot(l,k)
          amp1(l,k) = c1
          amp2(l,k) = c2
         ENDDO
         IF (kbot .LT. Nr) THEN
          Estop(l,kbot+1) = Esbot(l,kbot)
          Eutop(l,kbot+1) = Eubot(l,kbot)
         ENDIF

C       endif thresh
        ENDIF

        DO k = 1,Nr
C convert to scalar irradiance in quanta
         Etopwq = (rmud*Edtop(l,k)+rmus*Estop(l,k)+rmuu*Eutop(l,k))
     &            *WtouEins(l)
         Ebotwq = (rmud*Edbot(l,k)+rmus*Esbot(l,k)+rmuu*Eubot(l,k))
     &            *WtouEins(l)
C        and interpolate
         PAR(k,l) = SQRT(Etopwq*Ebotwq)
        ENDDO

C      enddo l
       ENDDO
C     endif kbot.gt.0
      ENDIF

#endif /* ALLOW_RADTRANS */

      return
      end

