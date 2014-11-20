C $Header$
C $Name$

#include "GUD_OPTIONS.h"

CBOP
C !ROUTINE: GUD_INIT_CHL
C !INTERFACE: ==========================================================
      SUBROUTINE GUD_INIT_CHL( dTsub,
     I                         bi, bj, myTime, myIter, myThid )

C !DESCRIPTION:
C     If gud_chlInitBalanced is true, initialize Chl from c
C     assuming balanced growth at PAR.
C     In any case, make sure Chl/C is between chl2cmin and chl2cmax.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#ifdef ALLOW_GUD
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"
#include "GUD_EXF.h"
#include "GUD_TRAITS.h"
#include "GUD_FIELDS.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid :: thread number
C  PAR    :: photosynthetically available radiation for balanced growth
      _RL dTsub(Nr)
      _RL myTime
      INTEGER bi, bj, myIter, myThid
CEOP

#ifdef ALLOW_GUD

C!LOCAL VARIABLES: ====================================================
      INTEGER j,l
      _RL PAR(sNx, sNy, Nr, nlam)
      _RL alpha_I(sNx,sNy,Nr)

C ======================================================================
C     compute PAR from current Chl
C     (usually initialized to 0 if gud_chlInitBalanced)
      IF (gud_chlInitBalanced) THEN
#ifdef ALLOW_RADTRANS
       CALL GUD_LIGHT_RADTRANS(PAR, bi, bj, myTime, myIter, myThid)
#else
C!!!!! do not take ice into account for backwards compatibility
       j = -1 
       CALL GUD_LIGHT(PAR, dTsub, bi, bj, myTime, j, myThid)
#endif
      ENDIF

      DO j = 1, nChl
       IF (gud_chlInitBalanced) THEN
        alpha_I = 0.0
        DO l = 1, nlam
         alpha_I = alpha_I + alphachl(j,l) * PAR(:,:,:,l)
        ENDDO
        Ptracer(1:sNx,1:sNy,:,bi,bj,iChl+j-1) =
     &  Ptracer(1:sNx,1:sNy,:,bi,bj,ic+j-1)*
     &         chl2cmax(j)/(1+chl2cmax(j)*alpha_I/(2*PCmax(j)))
       ENDIF
       Ptracer(1:sNx,1:sNy,:,bi,bj,iChl+j-1) =
     &            MAX(Ptracer(1:sNx,1:sNy,:,bi,bj,ic+j-1)*chl2cmin(j),
     &                Ptracer(1:sNx,1:sNy,:,bi,bj,iChl+j-1))
       Ptracer(1:sNx,1:sNy,:,bi,bj,iChl+j-1) =
     &            MIN(Ptracer(1:sNx,1:sNy,:,bi,bj,ic+j-1)*chl2cmax(j),
     &                Ptracer(1:sNx,1:sNy,:,bi,bj,iChl+j-1))
      ENDDO

#endif /* ALLOW_GUD */

      RETURN
      END
