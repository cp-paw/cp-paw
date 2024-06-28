!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CPPAW_WRITEVERSION(NFIL)
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      CHARACTER(512) :: RELEASE=''
      CHARACTER(512) :: REMOTE='UNKNOWN'
      CHARACTER(512) :: HASH='UNKNOWN'
      CHARACTER(512) :: BRANCH='UNKNOWN'
      CHARACTER(512) :: SHORTREVISIONNUMBER='UNKNOWN'
      CHARACTER(512) :: AUTHOR='UNKNOWN'
      CHARACTER(512) :: COMMITDATE='UNKNOWN'
      CHARACTER(512) :: COMPILEDATE='UNKNOWN'
      CHARACTER(512) :: COMPILEPERSON='UNKNOWN'
      INTEGER        :: NUMCHANGES=0
!     **************************************************************************
      INCLUDE 'CPPAW_VERSION.INFO'  ! FILENAME MADE LOWERCASE BY F90PP.SED 
!
      WRITE(NFIL,FMT='(82("*"),T25,"  CPPAW VERSION INFO  ")')
      IF(RELEASE.EQ.'') THEN
        WRITE(NFIL,FMT='(A)')-'THIS IS A DEVELOPMENT VERSION'
        WRITE(NFIL,FMT='(A)')-'COMPILED BY  '//TRIM(COMPILEPERSON) &
     &                          //'  ON '//TRIM(COMPILEDATE)
      ELSE
        WRITE(NFIL,FMT='(A)')-'THIS IS RELEASE '//TRIM(RELEASE)
      END IF
      WRITE(NFIL,FMT='(A)')-'HASH=        '//TRIM(HASH)
      WRITE(NFIL,FMT='(A)')-'BRANCH=      '//TRIM(BRANCH)
      WRITE(NFIL,FMT='(A)')-'REMOTE=      '//TRIM(REMOTE)
      WRITE(NFIL,FMT='(A)')-'COMMITTED ON='//TRIM(COMMITDATE)
      WRITE(NFIL,FMT='(A)')-'COMMITTED BY='//TRIM(AUTHOR)
      IF(NUMCHANGES.GT.0) THEN
        WRITE(NFIL,FMT='(I6,A)')NUMCHANGES,-' CHANGES SINCE LAST COMMIT'
      ELSE
        WRITE(NFIL,FMT='(A)')-'NO CHANGES SINCE LAST COMMIT'
      END IF
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
      CHARACTER(1),PARAMETER :: DOLLARSYM=ACHAR(36)
      INTEGER(4)      :: ISVAR
      CHARACTER(512)  :: CNTLNAME
      CHARACTER(32)   :: CH32SVAR
      INTEGER         :: ST,LENG
      CHARACTER(128)  :: CMD
      CHARACTER(512)  :: CMDMSG
      INTEGER(4)      :: EXITSTAT,CMDSTAT
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
        WRITE(*,FMT='(T25,"PARMFILE IS WRITTEN TO FILE ",A)')-'PARMS.IN_USE'
        CALL VERSION$WRITEPARMFILE(-'PARMS.IN_USE')
        WRITE(*,FMT='(T25,"ERRORS AFTER THIS LINE ARE IRRELEVANT")')
        CALL ERROR$NORMALSTOP
!
      ELSE IF(+CNTLNAME.EQ.'--SRCBLOB'.OR.+CNTLNAME.EQ.'-B') THEN
!       ========================================================================
!       ==  EXTRACT SRCBLOB FROM EXECUTABLE ====================================
!       ==  PAW_GETSRC.SH IS A BASH SCRIPT OF TEH CPPAW DISTRIBUTION THAT     ==
!       ==  EXTRACTS THE SRCBLOB, WHICH IS EMBEDDED IN THE EXECUTABLE         ==
!       ==  AND PLACES IT INTO PAW_SRCBLOB.TGZ                                ==
!       ========================================================================
        CALL GET_COMMAND_ARGUMENT(0,CMD)  ! COLLECT NAME OF EXECUTABLE
        CMD=-'WHICH '//TRIM(CMD)
        CMD=DOLLARSYM//'('//TRIM(CMD)//')'
        CMD=-'PAW_GETSRC.SH -X '//TRIM(CMD)//-' -O PAW_SRCBLOB.TGZ >&2'
        CMDMSG=''
!       == REMOVE EXISTING SOURCEBLOB TO AVOID SAFETY FEATURE OF PAW_GETSRC.SH =
        CALL EXECUTE_COMMAND_LINE(-'RM -F PAW_SRCBLOB.TGZ')
        CALL EXECUTE_COMMAND_LINE(TRIM(CMD),.TRUE.,EXITSTAT,CMDSTAT,CMDMSG)
        IF(EXITSTAT.NE.0) THEN
          CALL ERROR$MSG('COULD NOT EXTRACT SRCBLOB')
          CALL ERROR$CHVAL('CMD',TRIM(CMD))
          CALL ERROR$I4VAL('EXITSTAT',EXITSTAT)
          CALL ERROR$I4VAL('CMDSTAT',CMDSTAT)
!          CALL ERROR$MSG(TRIM(CMDMSG))
          CALL ERROR$STOP('PAW_VERSION')
        END IF
        CALL ERROR$NORMALSTOP
      END IF
      RETURN
      END
