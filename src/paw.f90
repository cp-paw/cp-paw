!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM MAIN
!     **************************************************************************
!     **                      CP-PAW                                          **
!     **     (C) COPYRIGHT CLAUSTHAL UNIVERSITY OF TECHNOLOGY                 **
!     **     LICENSED UNDER THE GNU GENERAL PUBLIC LICENSE V3                 **
!     **************************************************************************
!     **  THIS THE CAR-PARRINELLO FIRST PRINCIPLES MOLECULAR DYNAMICS         **
!     **  PROGRAM BASED ON THE PROJECTOR-AUGMENTED PLANE WAVE METHOD          **
!     **  AND THE LOCAL DENSITY APPROXIMATION.                                **
!     **                                                                      **
!     **  PLEASE READ THE PAW.README FILE CONTAINING DISCLAIMER               **
!     **  AND USER AGGREEMENT!                                                **
!     **                                                                      **
!     **  AUTHOR: PETER E. BLOECHL                                            **
!     **************************************************************************
      USE CLOCK_MODULE
      IMPLICIT NONE
      CHARACTER(32) :: DATIME
      INTEGER(4)    :: NFILO
      INTEGER(4)    :: NTASKS,THISTASK
      LOGICAL       :: DEBUGWAIT
!     **************************************************************************
!
!     ==========================================================================
!     == INITIALIZE MPE ROUTINE FOR PARALLEL PROCESSING                       ==
!     ==========================================================================
      CALL MPE$INIT
!
!     ==========================================================================
!     == CONSISTENCY CHECKS OF THE CALLING SEQUENCY AND VERSION INFO          ==
!     == PAW_VERSION MUST BE CALLED AFTER MPI$INIT FOR PROPER ERROR HANDLING  ==
!     ==========================================================================
      CALL PAW_VERSION()
!
!     ==========================================================================
!     == ENDLESS LOOP FOR PARALLEL DEBUGGING                                  ==
!     ==========================================================================
      DEBUGWAIT=.FALSE.  !FOR PARALLEL DEBUGGING SET EQUAL TRUE
      DO WHILE (DEBUGWAIT) !SET BREAKPOINT HERE
      ENDDO
!
!     ==========================================================================
!     ==  ENTER CAR-PARRINELLO SIMULATION                                     ==
!     ==========================================================================
      CALL TRACE$PUSH('MAIN')
                              CALL TIMING$START
      CALL PAW
!
!     ==========================================================================
!     ==  END OF CARPARRINELLO CALCULATION                                    ==
!     ==========================================================================
!     ==PRINTING IS ALLOWED ONLY ON THE FIRST TASK =============================
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      IF(THISTASK.EQ.1) THEN      
        CALL FILEHANDLER$REPORT(NFILO,'USED')
      END IF
!     == TIMING MUST BE CALLED BY ALL NODES ==================================
      CALL TRACE$PASS('BEFORE TIMING')
      CALL TIMING$PRINT('MONOMER',NFILO)
      IF(THISTASK.EQ.1) THEN      
        WRITE(NFILO,*)'REMARK: THE TIMING REPORT DOES NOT CONSIDER THE FIRST ' &
     &              ,'ITERATION'
        WRITE(NFILO,*)'        TO EXCLUDE THE TIME FOR SELF-INITIALIZATION'
      END IF
      CALL TRACE$PASS('AFTER TIMING')
      CALL MPE$CLOCKREPORT(NFILO)
!
!     ==PRINTING IS ALLOWED ONLY ON THE FIRST TASK =============================
      IF(THISTASK.EQ.1) THEN
        CALL CLOCK$NOW(DATIME)        
        WRITE(NFILO,FMT='(80("="))')
        WRITE(NFILO,FMT='(80("="),T15,"  PROGRAM FINISHED ",A32,"  ")')DATIME
        WRITE(NFILO,FMT='(80("="))')
      END IF
!
!     ==========================================================================
!     ==  CLOSE DOWN                                                          ==
!     ==========================================================================
      CALL TRACE$POP
      CALL ERROR$NORMALSTOP 
      STOP
      END
