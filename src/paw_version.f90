!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CPPAW_WRITEVERSION(NFIL)
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      CHARACTER(512) :: REMOTE='UNKNOWN'
      CHARACTER(512) :: HASH='UNKNOWN'
      CHARACTER(512) :: BRANCH='UNKNOWN'
      CHARACTER(512) :: SHORTREVISIONNUMBER='UNKNOWN'
      CHARACTER(512) :: AUTHOR='UNKNOWN'
      CHARACTER(512) :: COMMITDATE='UNKNOWN'
      CHARACTER(512) :: COMPILEDATE='UNKNOWN'
      CHARACTER(512) :: COMPILEPERSON='UNKNOWN'
!     **************************************************************************
      INCLUDE 'CPPAW_VERSION.INFO'  ! FILENAME MADE LOWERCASE BY F90PP.SED 
!
      WRITE(NFIL,FMT='(82("*"),T25,"  CPPAW VERSION INFO  ")')
      WRITE(NFIL,FMT='(A)')'HASH  =      '//TRIM(HASH)
      WRITE(NFIL,FMT='(A)')'BRANCH=      '//TRIM(BRANCH)
      WRITE(NFIL,FMT='(A)')'REMOTE=      '//TRIM(REMOTE)
      WRITE(NFIL,FMT='(A)')'COMMITTED ON='//TRIM(COMMITDATE)
      WRITE(NFIL,FMT='(A)')'COMMITTED BY='//TRIM(AUTHOR)
      WRITE(NFIL,FMT='(A)')'COMPILED BY  '//TRIM(COMPILEPERSON) &
     &                          //'  ON '//TRIM(COMPILEDATE)
      WRITE(NFIL,FMT='(82("*"))')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CPPAW_HELPMESSAGE(NFIL)
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      CHARACTER(512) :: CH32SVAR
!     **************************************************************************
      WRITE(NFIL,FMT='(82("*"),T25,"  HELP MESSAGE  ")')
      WRITE(NFIL,FMT='("USAGE:"/"=====")')
      WRITE(NFIL,FMT='(T5,"1ST ARGUMENT",T25,"ACTION"/T5,32("-"))')
      WRITE(NFIL,FMT='(T5,"CONTROL FILENAME",T25,"EXECUTE ")')
      CH32SVAR=-'--HELP'
      WRITE(NFIL,FMT='(T5,A,T25,"PRINT HELP INFORMATION")')TRIM(CH32SVAR)
      CH32SVAR=-'-H'
      WRITE(NFIL,FMT='(T5,A,T25,"PRINT HELP INFORMATION")')TRIM(CH32SVAR)
      WRITE(NFIL,FMT='(T5,"?",T25,"PRINT HELP INFORMATION")')
      CH32SVAR=-'--VERSION'
      WRITE(NFIL,FMT='(T5,A,T25,"PRINT VERSION INFORMATION")')TRIM(CH32SVAR)
      CH32SVAR=-'-V'
      WRITE(NFIL,FMT='(T5,A,T25,"PRINT VERSION INFORMATION")')TRIM(CH32SVAR)
      CH32SVAR=-'--PARMFILE'
      WRITE(NFIL,FMT='(T5,A,T25,"PRINT PARMFILE")')TRIM(CH32SVAR)
      CH32SVAR=-'-P'
      WRITE(NFIL,FMT='(T5,A,T25,"PRINT PARMFILE")')TRIM(CH32SVAR)
      WRITE(NFIL,FMT='(82("*"))')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAW_VERSION()
!     **************************************************************************
!     ** CHECKS CALLING SEQUENCE AND RESOLVES VERSION INFORMATION
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)      :: ISVAR
      CHARACTER(512)  :: CNTLNAME
      CHARACTER(32)   :: CH32SVAR
      INTEGER         :: ST,LENG
!     **************************************************************************
!
!     ==========================================================================
!     == CHECK THAT AT LEAST ONE ARGUMENT IS GIVEN                            ==
!     ==========================================================================
      ISVAR=COMMAND_ARGUMENT_COUNT()
      IF(ISVAR.LT.1) THEN
        CALL CPPAW_HELPMESSAGE(6)
        CALL ERROR$MSG('NO ARGUMENTS ARE GIVEN')
        CALL ERROR$STOP('PAW_VERSION')
      END IF
!
!     ==========================================================================
!     == CHECK THAT AT LEAST ONE ARGUMENT IS GIVEN                            ==
!     ==========================================================================
      CALL GET_COMMAND_ARGUMENT(1,CNTLNAME,LENG,ST)
      IF(ST.NE.0) THEN
        CALL ERROR$MSG('FAILURE COLLECTING COMMAND LINE ARGUMENT')
        CALL ERROR$I4VAL('ARGUMENT LENGTH',LENG)
        CALL ERROR$I4VAL('STATUS',ST)
        CALL ERROR$STOP('PAW_VERSION')
      END IF

     IF(+CNTLNAME.EQ.'--HELP'.OR.+CNTLNAME.EQ.'-H'.OR.+CNTLNAME.EQ.'?') THEN
        CALL CPPAW_HELPMESSAGE(6)
        WRITE(*,FMT='(T25,"ERRORS AFTER THIS LINE ARE IRRELEVANT")')
        CALL ERROR$NORMALSTOP
!
      ELSE IF(+CNTLNAME.EQ.'--VERSION'.OR.+CNTLNAME.EQ.'-V') THEN
        CALL CPPAW_WRITEVERSION(6)
        WRITE(*,FMT='(T25,"ERRORS AFTER THIS LINE ARE IRRELEVANT")')
        CALL ERROR$NORMALSTOP
!
      ELSE IF(+CNTLNAME.EQ.'--PARMFILE'.OR.+CNTLNAME.EQ.'-P') THEN
        CALL VERSION$WRITEPARMFILE()
        WRITE(*,FMT='(T25,"ERRORS AFTER THIS LINE ARE IRRELEVANT")')
        CALL ERROR$NORMALSTOP
      END IF
      RETURN
      END
