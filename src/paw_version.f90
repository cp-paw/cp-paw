!      MODULE VERSION_MODULE
!        CHARACTER(256):: VERTYP=_VERTYP
!        CHARACTER(256):: VERINF=_VERINF
!        CHARACTER(256):: VERREV=_VERREV
!        CHARACTER(256):: VERAUT=_VERAUT
!        CHARACTER(256):: VERDAT=_VERDAT
!      END MODULE VERSION_MODULE
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

        WRITE(*,FMT='(A)')TRIM(VERTYP)
        WRITE(*,FMT='(A)')TRIM(VERINF)
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
