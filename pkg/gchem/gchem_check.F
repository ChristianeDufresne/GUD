C $Header$
C $Name$

#include "GCHEM_OPTIONS.h"

      SUBROUTINE GCHEM_CHECK( myThid )
C     *==========================================================*
C     | SUBROUTINE GCHEM_CHECK
C     | o Validate basic package setup and inter-package
C     |   dependencies.
C     *==========================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GCHEM.h"

C     === Routine arguments ===
C     myThid   ::  Number of this instance
      INTEGER myThid

#ifdef ALLOW_GCHEM

C     === Local variables ===
C     msgBuf   :: Informational/error message buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER numPkg

      _BEGIN_MASTER(myThid)

      WRITE(msgBuf,'(A)') 'GCHEM_CHECK: #define ALLOW_GCHEM'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &     SQUEEZE_RIGHT , myThid )

#ifdef GCHEM_SEPARATE_FORCING
      IF ( useCFC ) THEN
         WRITE(msgBuf,'(A)') 'GCHEM_SEPARATE_FORCING not coded for CFC'
         CALL PRINT_ERROR( msgBuf , myThid )
         STOP 'ABNORMAL END: S/R GCHEM_CHECK'
      ENDIF
#endif

C     GCHEM needs passive tracer package
      IF ( .NOT.usePTRACERS) THEN
         WRITE(msgBuf,'(A)') 'usePTRACERS has to be .TRUE. for GCHEM'
         CALL PRINT_ERROR( msgBuf , myThid )
         STOP 'ABNORMAL END: S/R GCHEM_CHECK'
      ENDIF

C check that we are not running multiple incompatible packages
      numPkg = 0
      IF ( useDIC ) numPkg = numPkg + 1
      IF ( useCFC ) numPkg = numPkg + 1
      IF ( useDARWIN ) numPkg = numPkg + 1
      IF ( useGUD ) numPkg = numPkg + 1
      IF ( numPkg .GT. 1 ) THEN
        WRITE(msgBuf,'(2A)') 'only one of useDIC, useCFC,',
     &    ' useDARWIN and useGUD can be .TRUE. at a time'
        CALL PRINT_ERROR( msgBuf , myThid )
        STOP 'ABNORMAL END: S/R GCHEM_CHECK'
      ENDIF

#ifdef ALLOW_CFC
C     run checks specific to CFC
      IF ( useCFC ) THEN
        CALL CFC_CHECK( myThid )
      ENDIF
#endif

#ifdef ALLOW_DARWIN
C     run checks specific to darwin
      IF ( useDARWIN ) THEN
        CALL DARWIN_CHECK( myThid )
      ENDIF
#endif

#ifdef ALLOW_GUD
C     run checks specific to gud
      IF ( useGUD ) THEN
        CALL GUD_CHECK( myThid )
      ENDIF
#endif

      _END_MASTER(myThid)

#endif /* ALLOW_GCHEM */

      RETURN
      END
