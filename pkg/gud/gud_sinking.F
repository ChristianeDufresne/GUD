C $Header$
C $Name$

#include "GUD_OPTIONS.h"

CBOP
C !ROUTINE: GUD_SINKING
C !INTERFACE: ==========================================================
      SUBROUTINE GUD_SINKING(
     I     Ptr,
     U     gTr,
     I     bi,bj,myTime,myIter,myThid)

C !DESCRIPTION:
C     compute tendencies from sinking of particulate organic matter

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
c#include "EEPARAMS.h"
#include "GRID.h"
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"
#include "GUD_GENPARAMS.h"
#include "GUD_TRAITS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid :: thread number
C  Ptr    :: gud model tracers
C  dT     :: time step (for fake tendency with useEquilibriumC)
      _RL Ptr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nGud)
      INTEGER bi,bj
      INTEGER myThid, myIter
      _RL myTime

C !OUTPUT PARAMETERS: ===================================================
C  gTr    :: computed tendencies
      _RL gTr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nGud)
CEOP

#ifdef ALLOW_GUD

c !LOCAL VARIABLES: ====================================================
      INTEGER j,k
      _RL dzup(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL dzdn(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL flux(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      DO k=1,Nr-1
       dzup = DRF(k)*hFacC(:,:,k,bi,bj)
       dzdn = DRF(k+1)*hFacC(:,:,k+1,bi,bj)
       WHERE( dzdn.GT.0.0 )
         flux = wPIC_sink*MAX(0.0, Ptr(:,:,k,bi,bj,iPIC))
         gTr(:,:,k  ,bi,bj,iPIC ) = gTr(:,:,k  ,bi,bj,iPIC ) - flux/dzup
         gTr(:,:,k+1,bi,bj,iPIC ) = gTr(:,:,k+1,bi,bj,iPIC ) + flux/dzdn
         flux = wC_sink*MAX(0.0, Ptr(:,:,k,bi,bj,iPOC))
         gTr(:,:,k  ,bi,bj,iPOC ) = gTr(:,:,k  ,bi,bj,iPOC ) - flux/dzup
         gTr(:,:,k+1,bi,bj,iPOC ) = gTr(:,:,k+1,bi,bj,iPOC ) + flux/dzdn
         flux = wN_sink*MAX(0.0, Ptr(:,:,k,bi,bj,iPON))
         gTr(:,:,k  ,bi,bj,iPON ) = gTr(:,:,k  ,bi,bj,iPON ) - flux/dzup
         gTr(:,:,k+1,bi,bj,iPON ) = gTr(:,:,k+1,bi,bj,iPON ) + flux/dzdn
         flux = wP_sink*MAX(0.0, Ptr(:,:,k,bi,bj,iPOP))
         gTr(:,:,k  ,bi,bj,iPOP ) = gTr(:,:,k  ,bi,bj,iPOP ) - flux/dzup
         gTr(:,:,k+1,bi,bj,iPOP ) = gTr(:,:,k+1,bi,bj,iPOP ) + flux/dzdn
         flux = wSi_sink*MAX(0.0, Ptr(:,:,k,bi,bj,iPOSi))
         gTr(:,:,k  ,bi,bj,iPOSi) = gTr(:,:,k  ,bi,bj,iPOSi) - flux/dzup
         gTr(:,:,k+1,bi,bj,iPOSi) = gTr(:,:,k+1,bi,bj,iPOSi) + flux/dzdn
         flux = wFe_sink*MAX(0.0, Ptr(:,:,k,bi,bj,iPOFe))
         gTr(:,:,k  ,bi,bj,iPOFe) = gTr(:,:,k  ,bi,bj,iPOFe) - flux/dzup
         gTr(:,:,k+1,bi,bj,iPOFe) = gTr(:,:,k+1,bi,bj,iPOFe) + flux/dzdn
       ENDWHERE
       DO j = 1, nplank
         WHERE( dzdn.GT.0.0 )
      flux = wsink(j)*MAX(0.0, Ptr(:,:,k,bi,bj,ic+j-1))
      gTr(:,:,k  ,bi,bj,ic+j-1 )=gTr(:,:,k  ,bi,bj,ic+j-1 ) - flux/dzup
      gTr(:,:,k+1,bi,bj,ic+j-1 )=gTr(:,:,k+1,bi,bj,ic+j-1 ) + flux/dzdn
      flux = wswim(j)*MAX(0.0, Ptr(:,:,k+1,bi,bj,ic+j-1))
      gTr(:,:,k  ,bi,bj,ic+j-1 )=gTr(:,:,k  ,bi,bj,ic+j-1 ) + flux/dzup
      gTr(:,:,k+1,bi,bj,ic+j-1 )=gTr(:,:,k+1,bi,bj,ic+j-1 ) - flux/dzdn
         ENDWHERE
         IF (j .LE. nnquota) THEN
           WHERE( dzdn.GT.0.0 )
      flux = wsink(j)*MAX(0.0, Ptr(:,:,k,bi,bj,in+j-1))
      gTr(:,:,k  ,bi,bj,in+j-1 )=gTr(:,:,k  ,bi,bj,in+j-1 ) - flux/dzup
      gTr(:,:,k+1,bi,bj,in+j-1 )=gTr(:,:,k+1,bi,bj,in+j-1 ) + flux/dzdn
      flux = wswim(j)*MAX(0.0, Ptr(:,:,k+1,bi,bj,in+j-1))
      gTr(:,:,k  ,bi,bj,in+j-1 )=gTr(:,:,k  ,bi,bj,in+j-1 ) + flux/dzup
      gTr(:,:,k+1,bi,bj,in+j-1 )=gTr(:,:,k+1,bi,bj,in+j-1 ) - flux/dzdn
           ENDWHERE
         ENDIF
         IF (j .LE. npquota) THEN
           WHERE( dzdn.GT.0.0 )
      flux = wsink(j)*MAX(0.0, Ptr(:,:,k,bi,bj,ip+j-1))
      gTr(:,:,k  ,bi,bj,ip+j-1 )=gTr(:,:,k  ,bi,bj,ip+j-1 ) - flux/dzup
      gTr(:,:,k+1,bi,bj,ip+j-1 )=gTr(:,:,k+1,bi,bj,ip+j-1 ) + flux/dzdn
      flux = wswim(j)*MAX(0.0, Ptr(:,:,k+1,bi,bj,ip+j-1))
      gTr(:,:,k  ,bi,bj,ip+j-1 )=gTr(:,:,k  ,bi,bj,ip+j-1 ) + flux/dzup
      gTr(:,:,k+1,bi,bj,ip+j-1 )=gTr(:,:,k+1,bi,bj,ip+j-1 ) - flux/dzdn
           ENDWHERE
         ENDIF
         IF (j .LE. nsiquota) THEN
           WHERE( dzdn.GT.0.0 )
      flux = wsink(j)*MAX(0.0, Ptr(:,:,k,bi,bj,isi+j-1))
      gTr(:,:,k  ,bi,bj,isi+j-1)=gTr(:,:,k  ,bi,bj,isi+j-1) - flux/dzup
      gTr(:,:,k+1,bi,bj,isi+j-1)=gTr(:,:,k+1,bi,bj,isi+j-1) + flux/dzdn
      flux = wswim(j)*MAX(0.0, Ptr(:,:,k+1,bi,bj,isi+j-1))
      gTr(:,:,k  ,bi,bj,isi+j-1)=gTr(:,:,k  ,bi,bj,isi+j-1) + flux/dzup
      gTr(:,:,k+1,bi,bj,isi+j-1)=gTr(:,:,k+1,bi,bj,isi+j-1) - flux/dzdn
           ENDWHERE
         ENDIF
         IF (j .LE. nfequota) THEN
           WHERE( dzdn.GT.0.0 )
      flux = wsink(j)*MAX(0.0, Ptr(:,:,k,bi,bj,ife+j-1))
      gTr(:,:,k  ,bi,bj,ife+j-1)=gTr(:,:,k  ,bi,bj,ife+j-1) - flux/dzup
      gTr(:,:,k+1,bi,bj,ife+j-1)=gTr(:,:,k+1,bi,bj,ife+j-1) + flux/dzdn
      flux = wswim(j)*MAX(0.0, Ptr(:,:,k+1,bi,bj,ife+j-1))
      gTr(:,:,k  ,bi,bj,ife+j-1)=gTr(:,:,k  ,bi,bj,ife+j-1) + flux/dzup
      gTr(:,:,k+1,bi,bj,ife+j-1)=gTr(:,:,k+1,bi,bj,ife+j-1) - flux/dzdn
           ENDWHERE
         ENDIF
         IF (j .LE. nChl) THEN
           WHERE( dzdn.GT.0.0 )
      flux = wsink(j)*MAX(0.0, Ptr(:,:,k,bi,bj,iChl+j-1))
      gTr(:,:,k  ,bi,bj,iChl+j-1)=gTr(:,:,k  ,bi,bj,iChl+j-1)-flux/dzup
      gTr(:,:,k+1,bi,bj,iChl+j-1)=gTr(:,:,k+1,bi,bj,iChl+j-1)+flux/dzdn
      flux = wswim(j)*MAX(0.0, Ptr(:,:,k+1,bi,bj,iChl+j-1))
      gTr(:,:,k  ,bi,bj,iChl+j-1)=gTr(:,:,k  ,bi,bj,iChl+j-1)+flux/dzup
      gTr(:,:,k+1,bi,bj,iChl+j-1)=gTr(:,:,k+1,bi,bj,iChl+j-1)-flux/dzdn
           ENDWHERE
         ENDIF
       ENDDO ! j
      ENDDO ! k

#endif /* ALLOW_GUD */

      RETURN
      END

