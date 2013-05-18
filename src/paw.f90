Module version_module
!uses SVN keyword substitution
character(256):: VERInf='$HeadURL: https://yap/paw/branches/pbloechl/devel/src/paw.f90 $'
character(256):: VERrev='$LastChangedRevision: 1199 $'
character(256):: VERaut='$LastChangedBy: ptpb $'
character(256):: VERdat='$LastChangedDate: 2013-01-11 11:41:15 +0100 (Fri, 11 Jan 2013) $'
end Module version_module
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM MAIN
!     **************************************************************************
!     **                      CP-PAW                                          **
!     **     (C) COPYRIGHT   CLAUSTHAL UNIVERSITY OF TECHNOLOGY               **
!     **     LICENSED MATERIALS                                               **
!     **     PROPERTY OF CLAUSTHAL UNIVERSITY OF TECHNOLOGY                   **
!     **************************************************************************
!     **                                                                      **
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
!     == consistency checks of the calling sequency and version info          ==
!     ==========================================================================
      call paw_version()
!
!     ==========================================================================
!     == INITIALIZE MPE ROUTINE FOR PARALLEL PROCESSING                       ==
!     ==========================================================================
      CALL MPE$INIT
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
      IF(THISTASK.EQ.1) THEN      
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        CALL FILEHANDLER$REPORT(NFILO,'USED')
      END IF
!     == TIMING MUST BE CALLED BY ALL NODES ====================================
      CALL TRACE$PASS('BEFORE TIMING')
      CALL TIMING$PRINT('MONOMER',NFILO)
      WRITE(NFILO,*)'REMARK: THE TIMING REPORT DOES NOT CONSIDERTHE FIRST ' &
     &              ,'ITERATION TO EXCLUDE THE TIME FOR SELF-INITIALIZATION'
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAW_VERSION()
!     **************************************************************************
!     ** CHECKS CALLING SEQUENCE AND RESOLVES VERSION INFORMATION
!     **************************************************************************
      USE VERSION_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)      :: ISVAR
      CHARACTER(512)  :: CNTLNAME
      CHARACTER(32)   :: CH32SVAR
      CHARACTER(256)  :: VERSIONTEXT
      CHARACTER(256)  :: VERSIONINFO
      COMMON/VERSION/VERSIONTEXT
!     **************************************************************************
!     ==========================================================================
!     == THESE LINES CONTAIN INFORMATION THAT CAN BE GREPPED OUT OF THE       ==
!     == EXECUTABLE. THEY HAVE NO FUNCTION IN THE EXECUTION OF THE CODE.      ==
!     ==========================================================================
      VERSIONINFO = '@(#) PAW-VERSION %R% CREATED %U% %E%'
      VERSIONTEXT = 'PROGRAM VERSION %R% CREATED %U% %E%'
!
!     ==========================================================================
!     == CHECK THAT AT LEAST ONE ARGUMENT IS GIVEN                            ==
!     ==========================================================================
      CALL LIB$NARGS(ISVAR)
      IF(ISVAR.LT.1) THEN
        CALL ERROR$MSG('THE NAME OF THE CONTROLFILE')
        CALL ERROR$MSG('MUST BE GIVEN AS ARGUMENT')
        CALL ERROR$STOP('READIN')
      END IF
!
!     ==========================================================================
!     == CHECK THAT AT LEAST ONE ARGUMENT IS GIVEN                            ==
!     ==========================================================================
      CALL LIB$GETARG(1,CNTLNAME)
      IF(+CNTLNAME.EQ.'--HELP'.OR.+CNTLNAME.EQ.'-H'.OR.+CNTLNAME.EQ.'?') THEN
        WRITE(*,FMT='("USAGE:"/"=====")')
        WRITE(*,FMT='(T5,"1ST ARGUMENT",T25,"ACTION"/T5,32("-"))')
        WRITE(*,FMT='(T5,"CONTROL FILENAME",T25,"EXECUTE ")')
        CH32SVAR=-'--HELP'
        WRITE(*,FMT='(T5,A6,T25,"PRINT HELP INFORMATION")')TRIM(CH32SVAR)
        CH32SVAR=-'-H'
        WRITE(*,FMT='(T5,A2,T25,"PRINT HELP INFORMATION")')TRIM(CH32SVAR)
        WRITE(*,FMT='(T5,"?",T25,"PRINT HELP INFORMATION")')
        CH32SVAR=-'--VERSION'
        WRITE(*,FMT='(T5,A9,T25,"PRINT VERSION INFORMATION")')TRIM(CH32SVAR)
        CH32SVAR=-'--PARMFILE'
        WRITE(*,FMT='(T5,A10,T25,"PRINT PARMFILE")')TRIM(CH32SVAR)
        WRITE(*,FMT='(T5,A10,T25,"ERRORS AFTER THIS LINE ARE irrelevant")')
        CALL ERROR$NORMALSTOP
!
      ELSE IF(+CNTLNAME.EQ.'--VERSION') THEN
        WRITE(*,FMT='()')
        WRITE(*,FMT='(72("*"))')
        WRITE(*,FMT='(72("*"),T15' &
     &              //',"  CP-PAW VERSION INFO: ")')
        WRITE(*,FMT='(72("*"))')

        VERINF=TRIM(ADJUSTL(VERINF))
        WRITE(*,FMT='(A)')VERINF(1:11)//TRIM(VERINF(43:))
        WRITE(*,FMT='(A)')TRIM(VERREV)
        WRITE(*,FMT='(A)')TRIM(VERAUT)
        WRITE(*,FMT='(A)')TRIM(VERDAT)
        WRITE(*,FMT='(72("*"))')
        WRITE(*,FMT='(T5,A10,T25,"ERRORS AFTER THIS LINE ARE irRELEVANT")')
        CALL ERROR$NORMALSTOP
      ELSE IF(+CNTLNAME.EQ.'--PARMFILE') THEN
        CALL VERSION$WRITEPARMFILE()
        WRITE(*,FMT='(T5,A10,T25,"ERRORS AFTER THIS LINE ARE irRELEVANT")')
        CALL ERROR$NORMALSTOP
      END IF
      return
      end
