C $Header$
C $Name$

#include "GUD_OPTIONS.h"
#include "EXF_OPTIONS.h"

CBOP
C !ROUTINE: GUD_READPARMS

C !INTERFACE: ==========================================================
      SUBROUTINE GUD_READPARMS( myThid )

C !DESCRIPTION:
C     Initialize GUD parameters, read in data.gud

C !USES: ===============================================================
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#ifdef ALLOW_GCHEM
#include "GCHEM.h"
#endif
#ifdef ALLOW_EXF
#ifdef USE_EXF_INTERPOLATION
#include "SET_GRID.h"
#include "EXF_CONSTANTS.h"
#endif
#include "EXF_PARAM.h"
#endif
#include "GUD_SIZE.h"
#include "GUD_EXF.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none
CEOP

#ifdef ALLOW_GUD

C !LOCAL VARIABLES: ====================================================
C  iUnit   :: unit number for I/O
C  msgBuf  :: message buffer
      _RL inp_lon0, inp_lat0, inp_dLon, inp_dLat, gridNx, gridNy
      integer iUnit, oUnit, l
      CHARACTER*(MAX_LEN_MBUF) msgBuf

#include "GUD_EXF_NAMELIST.h"

C Open and read the data.gud file

      _BEGIN_MASTER(myThid)

#include "GUD_EXF_DEFAULTS.h"

      WRITE(msgBuf,'(A)') ' GUD_READPARMS: opening data.gud'
      CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
     &                   SQUEEZE_RIGHT , 1)

      CALL OPEN_COPY_DATA_FILE(
     I                   'data.gud', 'GUD_READPARMS',
     O                   iUnit,
     I                   myThid )

      READ(UNIT=iUnit, NML=gud_forcing_params)
#ifdef USE_EXF_INTERPOLATION
      READ(UNIT=iUnit, NML=gud_interp_params)
#endif

C ======================================================================
C     compute dependent parameters
      gud_loadFieldsEarly = gud_chlInitBalanced .AND.
     &                         nIter0 .EQ. gud_chlIter0

C ======================================================================
      IF ( myProcId.EQ.0 .AND. myThid.EQ.1 ) THEN
        CALL MDSFINDUNIT( oUnit, mythid )
        open(oUnit,file='gud_params.nml',status='unknown')
      ELSE
        oUnit = -1
      ENDIF

      CALL GUD_READGENPARAMS(iUnit, oUnit, myThid)
      CALL GUD_READTRAITPARAMS(iUnit, oUnit, myThid)

      IF ( oUnit .GE. 0 ) THEN
        close(oUnit)
      ENDIF
      CLOSE(iUnit)

      WRITE(msgBuf,'(A)') ' ==================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

      _END_MASTER(myThid)

C Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_GUD */

      RETURN
      END

