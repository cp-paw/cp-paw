!*******************************************************************************
!*******************************************************************************
!****                                                                      *****
!****  INTERFACES-ROUTINES FOR EXTERNAL LIBRARIES                          *****
!****                                                                      *****
!****  ALL EXTERNAL LIBRARIES ARE ONLY CALLED VIA THESE INTERFACES.        *****
!****                                                                      *****
!****  SPECIFIC LIBRARIES ARE SELECTED BY THE PREPROCESSOR STATEMENTS      *****
!****                                                                      *****
!****  CONTAINS:                                                           *****
!****    1) FORTRAN UTILITY LIBRARIES                                      *****
!****    2) LAPACK/BLAS LIBRARIES                                          *****
!****    3) FOURIER TRANSFORM LIBRARIES                                    *****
!****    4) SOME OTHER                                                     *****
!****                                                                      *****
!****  THE INTERFACES SHALL BE TESTED USING LIB$TEST                       *****
!****                                                                      *****
!*******************************************************************************
!*******************************************************************************
! 
!*******************************************************************************
!*******************************************************************************
!****                                                                      *****
!****              INTERFACES TO THE UTILITY LIBRARIES                     *****
!****                                                                      *****
!****  GETARG (GETARG)                                                     *****
!****  NARGS (IARGC)                                                       *****
!****  SYSTEM (SYSTEM)                                                     *****
!****  ETIME   (ETIME)                                                     *****
!****  GETHOSTNM (HOSTNM)                                                  *****
!****  FLUSH (FLUSH)                                                       *****
!****                                                                      *****
!****                                                                      *****
!****  REMARKS:                                                            *****
!****    THE XLF SUPPORT ROUTINES CARRY ONE TRAILING UNDERSCORE            *****
!****                                                                      *****
!****    THE SUPPORT LIBRARY ROUTINES OF FORTRAN DIFFER SLIGHTLY FROM      *****
!****    COMPILER TO COMPILER (ARGUMENT TYPES)                             *****
!****                                                                      *****
!****                                                                      *****
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$GETARG(IPOS,ARG)
!     **************************************************************************
!     **  RETURNS THE VALUE OF THE I-TH COMMAND LINE ARGUMENT                 ==
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: IPOS
      CHARACTER(*),INTENT(OUT) :: ARG
      INTEGER                  :: IPOSSTD
      integer                  :: leng
      integer                  :: st
!     **************************************************************************
!     == USE FORTRAN 2003 INTRINSIC SUBROUTINE
      CALL GET_COMMAND_ARGUMENT(IPOS,ARG,LENG,ST)
      IF(ST.NE.0) THEN
        CALL ERROR$MSG('FAILURE COLLECTING COMMAND LINE ARGUMENT')
        CALL ERROR$I4VAL('STATUS',ST)
        CALL ERROR$I4VAL('POSITION',IPOS)
        CALL ERROR$I4VAL('ACTUAL LENGTH OF ARGUMENT',LENG)
        CALL ERROR$STOP('LIB$GETARG')
      END IF
!
!     == use the commented lines below if fortran 2003 is not implemented
!!$#IF DEFINED(CPPVAR_COMPILER_G95)
!!$      CALL LIB_G95_GETARG(IPOS,ARG)
!!$#ELIF DEFINED(CPPVAR_COMPILER_IFC)
!!$      CALL LIB_IFC_GETARG(IPOS,ARG)
!!$#ELIF DEFINED(CPPVAR_COMPILER_IFC7)
!!$      CALL LIB_IFC7_GETARG(IPOS,ARG)
!!$#ELIF DEFINED(CPPVAR_COMPILER_ABSOFT)
!!$      CALL LIB_ABSOFT_GETARG(IPOS,ARG)
!!$#ELIF DEFINED(CPPVAR_COMPILER_XLF)
!!$      CALL LIB_XLF_GETARG(IPOS,ARG)
!!$#ELIF DEFINED(CPPVAR_COMPILER_PGI)
!!$      CALL LIB_PGI_GETARG(IPOS,ARG)
!!$#ELIF DEFINED(CPPVAR_COMPILER_PATHSCALE)
!!$      CALL LIB_PATHSCALE_GETARG(IPOS,ARG)
!!$#ELIF DEFINED(CPPVAR_COMPILER_SUN)
!!$      CALL LIB_SUN_GETARG(IPOS,ARG)
!!$#ELSE    
!!$      ! NO EXPLICIT INTERFACE; LET US HOPE FOR THE BEST....
!!$      CALL GETARG(IPOSSTD,ARG)
!!$#ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$NARGS(NARGS)
!     **************************************************************************
!     **  RETURNS THE NUMBER OF COMMAND-LINE ARGUMENTS OF THE MAIN ROUTINE    **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: NARGS
!     **************************************************************************
!
!     == use fortran 2003 intrinsic function ===================================
      NARGS=COMMAND_ARGUMENT_COUNT()
!
!     == use commented lines if Fortran 2003 is not available
!!$#IF DEFINED(CPPVAR_COMPILER_G95)
!!$      CALL LIB_G95_NARGS(NARGS)
!!$#ELIF DEFINED(CPPVAR_COMPILER_IFC)
!!$      CALL LIB_IFC_NARGS(NARGS)
!!$#ELIF DEFINED(CPPVAR_COMPILER_IFC7)
!!$      CALL LIB_IFC7_NARGS(NARGS)
!!$#ELIF DEFINED(CPPVAR_COMPILER_ABSOFT)
!!$      CALL LIB_ABSOFT_NARGS(NARGS)
!!$#ELIF DEFINED(CPPVAR_COMPILER_XLF)
!!$      CALL LIB_XLF_NARGS(NARGS)
!!$#ELIF DEFINED(CPPVAR_COMPILER_PGI)
!!$      CALL LIB_PGI_NARGS(NARGS)
!!$#ELIF DEFINED(CPPVAR_COMPILER_PATHSCALE)
!!$      CALL LIB_PATHSCALE_NARGS(NARGS)
!!$#ELIF DEFINED(CPPVAR_COMPILER_SUN)
!!$      CALL LIB_SUN_NARGS(NARGS)
!!$#ELSE          
!!$      ! NO EXPLICIT INTERFACE; LET US HOPE FOR THE BEST....
!!$      NARGS=IARGC()
!!$#ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$GETENV(NAME,VAL)
!     **************************************************************************
!     **  RETURNS THE VALUE OF THE I-TH COMMAND LINE ARGUMENT                 ==
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: NAME
      CHARACTER(*),INTENT(OUT) :: VAL
      INTEGER                  :: LENG
      INTEGER                  :: ST
!     **************************************************************************
      CALL GET_ENVIRONMENT_VARIABLE(NAME,VAL,LENG,ST)
      IF(ST.NE.0) THEN
        CALL ERROR$MSG('FAILURE COLLECTING ENVIORONMENT VARIABLE')
        IF(ST.EQ.11) THEN
          CALL ERROR$MSG('ENVIRONMENT VARIABLE DOES NOT EXIST')
        ELSE IF(ST.EQ.-1) THEN
          CALL ERROR$MSG('ENVIRONMENT VARIABLE DOES NOT FIT INTO STRING')
        END IF
        CALL ERROR$STOP('LIB$GETENV')
      end if
!
!     == use commented lines if Fortran 2003 is not available
! #IF DEFINED(CPPVAR_COMPILER_G95)
!       CALL LIB_G95_GETENV(NAME,VAL)
! #ELIF DEFINED(CPPVAR_COMPILER_IFC)
!       CALL LIB_IFC_GETENV(NAME,VAL)
!!$#ELIF DEFINED(CPPVAR_COMPILER_IFC7)
!!$      CALL LIB_IFC7_GETENV(NAME,VAL)
! #ELIF DEFINED(CPPVAR_COMPILER_ABSOFT)
!       CALL LIB_ABSOFT_GETENV(NAME,VAL)
!!$#ELIF DEFINED(CPPVAR_COMPILER_XLF)
!!$      CALL LIB_XLF_GETENV(NAME,VAL)
!!$#ELIF DEFINED(CPPVAR_COMPILER_PGI)
!!$      CALL LIB_PGI_GETENV(NAME,VAL)
!!$#ELIF DEFINED(CPPVAR_COMPILER_PATHSCALE)
!!$      CALL LIB_PATHSCALE_GETENV(NAME,VAL)
!!$#ELIF DEFINED(CPPVAR_COMPILER_SUN)
!!$      CALL LIB_SUN_GETENV(NAME,VAL)
!!$#ELSE    
!!$      ! NO EXPLICIT INTERFACE; USE  FORTRAN 2003 INTRINSIC FUNCTION
!!$      CALL GET_ENVIRONMENT_VARIABLE(NAME,VAL,LENG,ST)
!!$      IF(ST.NE.0) THEN
!!$        CALL ERROR$MSG('FAILURE COLLECTING ENVIORONMENT VARIABLE')
!!$        IF(ST.EQ.11) THEN
!!$          CALL ERROR$MSG('ENVIRONMENT VARIABLE DOES NOT EXIST')
!!$        ELSE IF(ST.EQ.-1) THEN
!!$          CALL ERROR$MSG('ENVIRONMENT VARIABLE DOES NOT FIT INTO STRING')
!!$        END IF
!!$        CALL ERROR$STOP('LIB$GETENV')
!!$      END IF 
!!$#ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$ETIME(USRTIME,SYSTIME)
!     **************************************************************************
!     **  RETURNS THE USER AND SYSTEM TIME OF THE CURRENT PROCESS             **
!     **                                                                      **
!     **    USRTIME IS THE ELAPSED CPU TIME IN SECONDS USED BY THE USER CODE  **
!     **    SYSTIME IS THE ELAPSED CPU TIME IN SECONDS USED BY                **
!     **                                          CALLS TO THE OPERATING SYSTEM**
!     **                                                                      **
!     **    TIME FOR IO, WHERE THE PROCESSOR IS IDLE IS NOT ACCOUNTED FOR!!!  **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: USRTIME
      REAL(8)   ,INTENT(OUT) :: SYSTIME
      REAL                   :: TARRAY(2),RESULT
!     **************************************************************************
#IF DEFINED(CPPVAR_COMPILER_G95)
      CALL LIB_G95_ETIME(USRTIME,SYSTIME)
#ELIF DEFINED(CPPVAR_COMPILER_IFC)
      CALL LIB_IFC_ETIME(USRTIME,SYSTIME)
#ELIF DEFINED(CPPVAR_COMPILER_IFC7)
      CALL LIB_IFC7_ETIME(USRTIME,SYSTIME)
#ELIF DEFINED(CPPVAR_COMPILER_ABSOFT)
      CALL LIB_ABSOFT_ETIME(USRTIME,SYSTIME)
#ELIF DEFINED(CPPVAR_COMPILER_XLF)
      CALL LIB_XLF_ETIME(USRTIME,SYSTIME)
#ELIF DEFINED(CPPVAR_COMPILER_PGI)
      CALL LIB_PGI_ETIME(USRTIME,SYSTIME)
#ELIF DEFINED(CPPVAR_COMPILER_PATHSCALE)
      CALL LIB_PATHSCALE_ETIME(USRTIME,SYSTIME)
#ELIF DEFINED(CPPVAR_COMPILER_SUN)
      CALL LIB_SUN_ETIME(USRTIME,SYSTIME)
#ELSE  
      ! NO EXPLICIT INTERFACE; LET US HOPE FOR THE BEST....
      RESULT=ETIME(TARRAY)
      USRTIME=TARRAY(1)
      SYSTIME=TARRAY(2)
#ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$GETHOSTNAME(HOSTNAME)
!     **************************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                     **
!     **************************************************************************
      CHARACTER(*),INTENT(OUT)  :: HOSTNAME
      INTEGER(4)                :: RC
!     **************************************************************************
#IF DEFINED(CPPVAR_COMPILER_G95)
      CALL LIB_G95_GETHOSTNAME(HOSTNAME)
#ELIF DEFINED(CPPVAR_COMPILER_IFC)
      CALL LIB_IFC_GETHOSTNAME(HOSTNAME)
#ELIF DEFINED(CPPVAR_COMPILER_IFC7)
      CALL LIB_IFC7_GETHOSTNAME(HOSTNAME)
#ELIF DEFINED(CPPVAR_COMPILER_ABSOFT)
      CALL LIB_ABSOFT_GETHOSTNAME(HOSTNAME)
#ELIF DEFINED(CPPVAR_COMPILER_XLF)
      CALL LIB_XLF_GETHOSTNAME(HOSTNAME)
#ELIF DEFINED(CPPVAR_COMPILER_PGI)
      CALL LIB_PGI_GETHOSTNAME(HOSTNAME)
#ELIF DEFINED(CPPVAR_COMPILER_PATHSCALE)
      CALL LIB_PATHSCALE_GETHOSTNAME(HOSTNAME)
#ELIF DEFINED(CPPVAR_COMPILER_SUN)
      CALL LIB_SUN_GETHOSTNAME(HOSTNAME)
#ELSE
      ! NO EXPLICIT INTERFACE; LET US HOPE FOR THE BEST....
      RC=HOSTNM(HOSTNAME)    
      IF(RC.NE.0)HOSTNAME='UNKNOWN'
#ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$SYSTEM(COMMAND)
!     **************************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                     **
!     **************************************************************************
      CHARACTER(*),INTENT(IN)   :: COMMAND
      integer(4)                :: rc   ! return code
!     *********************************************************************
      call LIB$SYSTEMrc(COMMAND,rc)
      IF(RC.NE.0) THEN
        CALL ERROR$MSG('SYSTEM CALL FAILED')
        CALL ERROR$CHVAL('COMMAND',COMMAND)
        CALL ERROR$I4VAL('RC',RC)
        CALL ERROR$STOP('LIB$SYSTEM')
      End IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$SYSTEMrc(COMMAND,rc)
!     **************************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                     **
!     **************************************************************************
      CHARACTER(*),INTENT(IN)   :: COMMAND
      integer(4)  ,intent(out)  :: rc   ! return code
!     *********************************************************************
#IF DEFINED(CPPVAR_COMPILER_G95)
      CALL LIB_G95_SYSTEM(COMMAND,rc)
#ELIF DEFINED(CPPVAR_COMPILER_IFC)
      CALL LIB_IFC_SYSTEM(COMMAND,rc)
#ELIF DEFINED(CPPVAR_COMPILER_IFC7)
      CALL LIB_IFC7_SYSTEM(COMMAND,rc)
#ELIF DEFINED(CPPVAR_COMPILER_ABSOFT)
      CALL LIB_ABSOFT_SYSTEM(COMMAND,rc)
#ELIF DEFINED(CPPVAR_COMPILER_XLF)
      CALL LIB_XLF_SYSTEM(COMMAND,rc)
#ELIF DEFINED(CPPVAR_COMPILER_PGI)
      CALL LIB_PGI_SYSTEM(COMMAND,rc)
#ELIF DEFINED(CPPVAR_COMPILER_PATHSCALE)
      CALL LIB_PATHSCALE_SYSTEM(COMMAND,rc)
#ELIF DEFINED(CPPVAR_COMPILER_SUN)
      CALL LIB_SUN_SYSTEM(COMMAND,rc)
#ELSE
      ! NO EXPLICIT INTERFACE; LET US HOPE FOR THE BEST....
      RC=SYSTEM(COMMAND)
#ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$FLUSHFILE(N)
!     **************************************************************************
!     ** FLUSHES THE BUFFER FOR THE FILE CONNECTED TO FORTRAN UNIT N          **
!     **************************************************************************
      INTEGER(4),INTENT(IN) :: N
!     **************************************************************************
#IF DEFINED(CPPVAR_COMPILER_G95)
      CALL LIB_G95_FLUSH(N)   
#ELIF DEFINED(CPPVAR_COMPILER_IFC)
      CALL LIB_IFC_FLUSH(N)   
#ELIF DEFINED(CPPVAR_COMPILER_IFC7)
      CALL LIB_IFC7_FLUSH(N)   
#ELIF DEFINED(CPPVAR_COMPILER_ABSOFT)
      CALL LIB_ABSOFT_FLUSH(N)   
#ELIF DEFINED(CPPVAR_COMPILER_XLF)
      CALL LIB_XLF_FLUSH(N)   
#ELIF DEFINED(CPPVAR_COMPILER_PGI)
      CALL LIB_PGI_FLUSH(N)   
#ELIF DEFINED(CPPVAR_COMPILER_PATHSCALE)
      CALL LIB_PATHSCALE_FLUSH(N)   
#ELIF DEFINED(CPPVAR_COMPILER_SUN)
      CALL LIB_SUN_FLUSH(N)   
#ELSE
      ! NO EXPLICIT INTERFACE; LET US HOPE FOR THE BEST....
      FLUSH(N)  ! fortran 2003
#ENDIF
      RETURN
      END 

#IF DEFINED(CPPVAR_COMPILER_IFC)
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC_NARGS(NARGS)
!     **************************************************************************
!     **  RETURNS THE NUMBER OF COMMAND-LINE ARGUMENTS OF THE MAIN ROUTINE    **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER                             **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC10. THERE IS A CHANGE FROM **
!     **           IFC7, THEREFORE THERE IS AN EXPLICIT INTERFACE IFC7)       **
!     **  THE MODULE IFPORT IS SUPPLIED BY THE IFC COMPILER                   **
!     **************************************************************************
      USE IFPORT          !INCLUDE FORTRAN PORTABILITY FUNCTIONS
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: NARGS
!     **************************************************************************
      NARGS=IARGC()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC_GETARG(IPOS,ARG)
!     **************************************************************************
!     **  RETURNS THE VALUE OF THE I-TH COMMAND LINE ARGUMENT                 ==
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER                             **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC10. THERE IS A CHANGE FROM **
!     **           IFC7, THEREFORE THERE IS AN EXPLICIT INTERFACE IFC7)       **
!     **  THE MODULE IFPORT IS SUPPLIED BY THE IFC COMPILER                   **
!     **************************************************************************
      USE IFPORT          !INCLUDE FORTRAN PORTABILITY FUNCTIONS
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: IPOS
      CHARACTER(*),INTENT(OUT) :: ARG
      INTEGER                  :: IPOSSTD
!     **************************************************************************
      IPOSSTD=IPOS
      CALL GETARG(IPOSSTD,ARG)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC_GETENV(NAME,VAL)
!     **************************************************************************
!     **  RETURNS THE VALUE OF AN ENVIRONMENT VARIABLE                        **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER                             **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC10. THERE IS A CHANGE FROM **
!     **           IFC7, THEREFORE THERE IS AN EXPLICIT INTERFACE IFC7)       **
!     **  THE MODULE IFPORT IS SUPPLIED BY THE IFC COMPILER                   **
!     **************************************************************************
      USE IFPORT          !INCLUDE FORTRAN PORTABILITY FUNCTIONS
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: NAME
      CHARACTER(*),INTENT(OUT) :: VAL
!     **************************************************************************
      CALL GETENV(NAME,VAL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC_ETIME(USRTIME,SYSTIME)
!     **************************************************************************
!     ** RETURNS THE USER AND SYSTEM ELAPSED TIME OF THE CURRENT PROCESS      **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER                             **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC10. THERE IS A CHANGE FROM **
!     **           IFC7, THEREFORE THERE IS AN EXPLICIT INTERFACE IFC7)       **
!     **  THE MODULE IFPORT IS SUPPLIED BY THE IFC COMPILER                   **
!     **************************************************************************
      USE IFPORT
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: USRTIME
      REAL(8),INTENT(OUT) :: SYSTIME
      REAL(4)             :: TARRAY(2)
      REAL(4)             :: RESULT
!     **************************************************************************
      RESULT=ETIME(TARRAY)
      USRTIME=TARRAY(1)
      SYSTIME=TARRAY(2)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC_SYSTEM(COMMAND,RC) 
!     **************************************************************************
!     **  ISSUES A SHELL COMMAND                                              **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER                             **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC10. THERE IS A CHANGE FROM **
!     **           IFC7, THEREFORE THERE IS AN EXPLICIT INTERFACE IFC7)       **
!     **  THE MODULE IFPORT IS SUPPLIED BY THE IFC COMPILER                   **
!     **************************************************************************
      USE IFPORT
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: COMMAND
      INTEGER(4)  ,intent(out) :: RC
!     **************************************************************************
      RC=SYSTEM(COMMAND)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC_GETHOSTNAME(HOSTNAME)
!     **************************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER                             **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC10. THERE IS A CHANGE FROM **
!     **           IFC7, THEREFORE THERE IS AN EXPLICIT INTERFACE IFC7)       **
!     **  THE MODULE IFPORT IS SUPPLIED BY THE IFC COMPILER                   **
!     **************************************************************************
      USE IFPORT
      IMPLICIT NONE
      CHARACTER(*),INTENT(OUT)  :: HOSTNAME
      INTEGER(4)                :: RC
!     **************************************************************************
      RC=HOSTNAM(HOSTNAME)    
      IF(RC.NE.0)HOSTNAME='UNKNOWN'
      HOSTNAME='UNKNOWN'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC_FLUSH(NFIL)
!     **************************************************************************
!     **  FLUSHES THE FILE CONNECTED TO UNIT NFIL                             **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER                             **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC10. THERE IS A CHANGE FROM **
!     **           IFC7, THEREFORE THERE IS AN EXPLICIT INTERFACE IFC7)       **
!     **  THE MODULE IFPORT IS SUPPLIED BY THE IFC COMPILER                   **
!     **************************************************************************
      USE IFPORT
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)       :: NFIL
      INTEGER(4)                  :: LUNIT
!     **************************************************************************
      LUNIT=NFIL
      CALL FLUSH(LUNIT)
      RETURN
      END
#ELIF DEFINED(CPPVAR_COMPILER_IFC7)
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC7_ETIME(USRTIME,SYSTIME)
!     **************************************************************************
!     ** RETURNS THE USER AND SYSTEM ELAPSED TIME OF THE CURRENT PROCESS      **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER  (IFC7 AND EARLIER)         **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC7. THERE IS A CHANGE TO    **
!     **           IFC10, THEREFORE THERE IS ANOTHER EXPLICIT INTERFACE IFC)  **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        REAL(4) FUNCTION ETIME(TARRAY)
        REAL(4),OPTIONAL,INTENT(OUT) :: TARRAY(2)
        END FUNCTION ETIME
      END INTERFACE
      REAL(8),INTENT(OUT) :: USRTIME
      REAL(8),INTENT(OUT) :: SYSTIME
      REAL(4)             :: TARRAY(2)
      REAL(4)             :: RESULT      
!     **************************************************************************
      RESULT=ETIME(TARRAY)
      USRTIME=TARRAY(1)
      SYSTIME=TARRAY(2)
      RETURN
      END SUBROUTINE LIB_IFC7_ETIME
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC7_GETARG(IPOS,ARG)
!     **************************************************************************
!     **  RETURNS THE VALUE OF THE I-TH COMMAND LINE ARGUMENT                 **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER  (IFC7 AND EARLIER)         **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC7. THERE IS A CHANGE TO    **
!     **           IFC10, THEREFORE THERE IS ANOTHER EXPLICIT INTERFACE IFC)  **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE GETARG(POS,VALUE)
        INTEGER(4)     ,INTENT(IN) :: POS
        CHARACTER(*)   ,INTENT(OUT):: VALUE
        END SUBROUTINE GETARG
      END INTERFACE
      INTEGER(4)      ,INTENT(IN)  :: IPOS
      CHARACTER(*)    ,INTENT(OUT) :: ARG
!     **************************************************************************
      CALL GETARG(IPOS,ARG)
      RETURN
      END SUBROUTINE LIB_IFC7_GETARG
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC7_NARGS(NARGS)
!     **************************************************************************
!     **  RETURNS THE NUMBER OF COMMAND-LINE ARGUMENTS OF THE MAIN ROUTINE    **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER  (IFC7 AND EARLIER)         **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC7. THERE IS A CHANGE TO    **
!     **           IFC10, THEREFORE THERE IS ANOTHER EXPLICIT INTERFACE IFC)  **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        INTEGER(4) FUNCTION IARGC()
        END FUNCTION IARGC
      END INTERFACE
      INTEGER(4),INTENT(OUT) :: NARGS
!     **************************************************************************
      NARGS=IARGC()
      RETURN
      END SUBROUTINE LIB_IFC7_NARGS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC7_SYSTEM(COMMAND,rc)
!     **************************************************************************
!     **  ISSUES A SHELL COMMAND                                              **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER  (IFC7 AND EARLIER)         **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC7. THERE IS A CHANGE TO    **
!     **           IFC10, THEREFORE THERE IS ANOTHER EXPLICIT INTERFACE IFC)  **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        INTEGER(4) FUNCTION SYSTEM(STRING)
        CHARACTER(*) ,INTENT(IN) :: STRING
        END FUNCTION SYSTEM
      END INTERFACE
      CHARACTER(*),INTENT(IN) :: COMMAND
      INTEGER(4)  ,intent(out) :: RC
!     **************************************************************************
      RC=SYSTEM(COMMAND)
      RETURN
      END SUBROUTINE LIB_IFC7_SYSTEM
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC7_GETHOSTNAME(HOSTNAME)
!     **************************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER  (IFC7 AND EARLIER)         **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC7. THERE IS A CHANGE TO    **
!     **           IFC10, THEREFORE THERE IS ANOTHER EXPLICIT INTERFACE IFC)  **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE HOSTNM(STRING)
        CHARACTER(*) ,INTENT(OUT) :: STRING
        END SUBROUTINE HOSTNM
      END INTERFACE
      CHARACTER(*),INTENT(OUT)  :: HOSTNAME
      INTEGER(4)                :: RC
!     *********************************************************************
      CALL HOSTNM(HOSTNAME)    
      RETURN
      END SUBROUTINE LIB_IFC7_GETHOSTNAME
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_IFC7_FLUSH(NFIL)
!     **************************************************************************
!     **  FLUSHES THE FILE CONNECTED TO UNIT NFIL                             **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER  (IFC7 AND EARLIER)         **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC7. THERE IS A CHANGE TO    **
!     **           IFC10, THEREFORE THERE IS ANOTHER EXPLICIT INTERFACE IFC)  **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE FLUSH(LUNIT)
        INTEGER,INTENT(IN) :: LUNIT
        END SUBROUTINE FLUSH
      END INTERFACE
      INTEGER(4),INTENT(IN)       :: NFIL
      INTEGER                     :: LUNIT
!     **************************************************************************
      LUNIT=NFIL
      CALL FLUSH(LUNIT)
      RETURN
      END SUBROUTINE LIB_IFC7_FLUSH
#ELIF DEFINED(CPPVAR_COMPILER_G95)
!     **************************************************************************
!     ** G95 OFFERS THE USUAL SUPPORT FUNCTIONS AS INTRINSIC FUNCTION         **
!     ** EXTENSIONS. THESE MUST BE ADRESSED WITHOUT AN EXPLICIT INTERFACE     **
!     ** TO AVOID THAT AN EXTERNAL LIBRARY IS CALLED INSTEAD.                 **
!     **                                                                      **
!     ** USING THE EXTENSIONS IS INCOMPATIBLE WITH THE -STD95 FLAG OF THE     **
!     ** FORTRAN COMPILER!                                                    **
!     **                                                                      **
!     ** SOME OF THE ROUTINES ARE ALSO AVAILABLE AS INTRINSIC FUNCTIONS OF    **
!     ** FORTRAN 2003.                                                        **
!     **************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_G95_ETIME(USRTIME,SYSTIME)
!     **************************************************************************
!     **  RETURNS THE USER AND SYSTEM ELAPSED TIME OF THE CURRENT PROCESS     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE G95 COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: USRTIME
      REAL(8),INTENT(OUT) :: SYSTIME
      REAL                :: TARRAY(2)
      REAL                :: RESULT      
!     **************************************************************************
      call etime(tarray,result)
      usrtime=tarray(1)
      systime=tarray(2)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_G95_GETARG(IPOS,ARG)
!     **************************************************************************
!     **  RETURNS THE VALUE OF THE I-TH COMMAND LINE ARGUMENT                 **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE G95 COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
!!$    ! DO NOT USE THE INTERFACE BECAUSE THIS IS A INTRINSIC FUNCTION EXTENSION
!!$    ! IN G95
!!$      INTERFACE 
!!$        SUBROUTINE GETARG(POS,VALUE)
!!$        INTEGER         ,INTENT(IN) :: POS
!!$        CHARACTER(LEN=*),INTENT(OUT):: VALUE
!!$        END SUBROUTINE GETARG
!!$      END INTERFACE
      INTEGER(4),INTENT(IN)         :: IPOS
      CHARACTER(LEN=*) ,INTENT(OUT) :: ARG
      INTEGER                       :: POS
      INTEGER                       :: LEN,ST
!     **************************************************************************
      POS=IPOS
!      CALL GETARG(POS,ARG)
!     ==  USE  FORTRAN 2003 INTRINSIC FUNCTION
      CALL GET_COMMAND_ARGUMENT(POS,ARG,LEN,ST)
      IF(ST.NE.0) THEN
        CALL ERROR$MSG('ARGUMENT NOT PRESENT')
        CALL ERROR$STOP('LIB_G95_GETARG')
      END IF 
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_G95_GETENV(NAME,VAL)
!     **************************************************************************
!     **  RETURNS THE VALUE OF THE I-TH COMMAND LINE ARGUMENT                 **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE G95 COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
!!$    ! DO NOT USE THE INTERFACE BECAUSE THIS IS A INTRINSIC FUNCTION EXTENSION
!!$    ! IN G95
!!$      INTERFACE 
!!$        SUBROUTINE GETENV(NAME,VALUE)
!!$        CHARACTER(LEN=*),INTENT(IN) :: NAME
!!$        CHARACTER(LEN=*),INTENT(OUT):: VALUE
!!$        END SUBROUTINE GETARG
!!$      END INTERFACE
      CHARACTER(LEN=*) ,INTENT(IN) :: NAME
      CHARACTER(LEN=*) ,INTENT(OUT) :: VAL
      INTEGER                       :: LENG,ST
!     **************************************************************************
!      CALL GETENV(NAME,VAL)
!      RETURN
!     ==  USE  FORTRAN 2003 INTRINSIC FUNCTION
      CALL GET_ENVIRONMENT_VARIABLE(NAME,VAL,LENG,ST)
      IF(ST.NE.0) THEN
        CALL ERROR$MSG('FAILURE COLLECTING ENVIORONMENT VARIABLE')
        IF(ST.EQ.11) THEN
          CALL ERROR$MSG('ENVIRONMENT VARIABLE DOES NOT EXIST')
        ELSE IF(ST.EQ.-1) THEN
          CALL ERROR$MSG('ENVIRONMENT VARIABLE DOES NOT FIT INTO STRING')
        END IF
        CALL ERROR$CHVAL('NAME',TRIM(NAME))
        CALL ERROR$STOP('LIB_G95_GETENV')
      END IF 
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_G95_NARGS(NARGS)
!     **************************************************************************
!     **  RETURNS THE NUMBER OF COMMAND-LINE ARGUMENTS OF THE MAIN ROUTINE    **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE G95 COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
!!$      INTERFACE 
!!$        INTEGER FUNCTION IARGC()
!!$        END FUNCTION IARGC
!!$      END INTERFACE
      INTEGER(4),INTENT(OUT) :: NARGS
      INTEGER                :: count
!     **************************************************************************
!!$      NARGS=IARGC()
      count=COMMAND_ARGUMENT_COUNT()
      nargs=count
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_G95_SYSTEM(COMMAND,rc)
!     **************************************************************************
!     **  ISSUES A SHELL COMMAND                                              **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE G95 COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
!!$      INTERFACE 
!!$        INTEGER FUNCTION SYSTEM(STRING)
!!$        CHARACTER*(*) ,INTENT(IN) :: STRING
!!$        END FUNCTION SYSTEM
!!$      END INTERFACE
      CHARACTER(*),INTENT(IN)  :: COMMAND
      INTEGER(4)  ,intent(out) :: RC
!     **************************************************************************
      RC=SYSTEM(COMMAND)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_G95_GETHOSTNAME(HOSTNAME)
!     **************************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE G95 COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
!!$      INTERFACE 
!!$        INTEGER FUNCTION HOSTNM(STRING)
!!$        CHARACTER*(*) ,INTENT(IN) :: STRING
!!$        END FUNCTION HOSTNM
!!$      END INTERFACE
      CHARACTER(*),INTENT(OUT)  :: HOSTNAME
      INTEGER(4)                :: RC
!     *********************************************************************
      RC=HOSTNM(HOSTNAME)    
      IF(RC.NE.0)HOSTNAME='UNKNOWN'
      HOSTNAME='UNKNOWN'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_G95_FLUSH(NFIL)
!     **************************************************************************
!     **  FLUSHES THE FILE CONNECTED TO UNIT NFIL                             **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE G95 COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
!!$      INTERFACE 
!!$        SUBROUTINE FLUSH(LUNIT)
!!$        INTEGER,INTENT(IN) :: LUNIT
!!$        END SUBROUTINE FLUSH
!!$      END INTERFACE
      INTEGER(4),INTENT(IN)       :: NFIL
      INTEGER                     :: LUNIT
!     *********************************************************************
      LUNIT=NFIL
      CALL FLUSH(LUNIT)
      RETURN
      END
#ELIF DEFINED(CPPVAR_COMPILER_ABSOFT)
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_ABSOFT_ETIME(USRTIME,SYSTIME)
!     **************************************************************************
!     **  RETURNS THE USER AND SYSTEM ELAPSED TIME OF THE CURRENT PROCESS     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE ABSOFT COMPILER                          **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        REAL*4 FUNCTION ETIME(TARRAY)
        REAL*4,INTENT(OUT) :: TARRAY(2)
        END FUNCTION ETIME
      END INTERFACE
      REAL(8),INTENT(OUT) :: USRTIME
      REAL(8),INTENT(OUT) :: SYSTIME
      REAL*4              :: TARRAY(2)
      REAL*4              :: RESULT      
!     **************************************************************************
      RESULT=ETIME(TARRAY)
      USRTIME=TARRAY(1)
      SYSTIME=TARRAY(2)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_ABSOFT_GETARG(IPOS,ARG)
!     **************************************************************************
!     **  RETURNS THE VALUE OF THE I-TH COMMAND LINE ARGUMENT                 **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE ABSOFT COMPILER                          **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE GETARG(POS,VALUE)
        INTEGER*4       ,INTENT(IN) :: POS
        CHARACTER(LEN=*),INTENT(OUT):: VALUE
        END SUBROUTINE GETARG
      END INTERFACE
      INTEGER(4),INTENT(IN)         :: IPOS
      CHARACTER(LEN=*) ,INTENT(OUT) :: ARG
      INTEGER*4                     :: POS
!     **************************************************************************
      POS=IPOS
      CALL GETARG(POS,ARG)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_ABSOFT_NARGS(NARGS)
!     **************************************************************************
!     **  RETURNS THE NUMBER OF COMMAND-LINE ARGUMENTS OF THE MAIN ROUTINE    **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE ABSOFT COMPILER                          **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        INTEGER*4 FUNCTION IARGC()
        END FUNCTION IARGC
      END INTERFACE
      INTEGER(4),INTENT(OUT) :: NARGS
!     **************************************************************************
      NARGS=IARGC()
PRINT*,'NARGS ',NARGS,IARGC()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_absoft_GETENV(NAME,VAL)
!     **************************************************************************
!     **  RETURNS THE VALUE OF AN ENVIRONMENT VARIABLE                        **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE IFC COMPILER                             **
!     **  (REMARK: THIS INTERFACE IS ADAPTED TO IFC10. THERE IS A CHANGE FROM **
!     **           IFC7, THEREFORE THERE IS AN EXPLICIT INTERFACE IFC7)       **
!     **  THE MODULE IFPORT IS SUPPLIED BY THE IFC COMPILER                   **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: NAME
      CHARACTER(*),INTENT(OUT) :: VAL
!     **************************************************************************
      CALL GETENV(NAME,VAL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_ABSOFT_SYSTEM(COMMAND,rc)
!     **************************************************************************
!     **  ISSUES A SHELL COMMAND                                              **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE ABSOFT COMPILER                          **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        INTEGER*4 FUNCTION SYSTEM(STRING)
        CHARACTER*(*) ,INTENT(IN) :: STRING
        END FUNCTION SYSTEM
      END INTERFACE
      CHARACTER(*),INTENT(IN) :: COMMAND
      INTEGER(4)  ,intent(out):: RC
!     **************************************************************************
      RC=SYSTEM(COMMAND)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_ABSOFT_GETHOSTNAME(HOSTNAME)
!     **************************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE ABSOFT COMPILER                          **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        INTEGER*4 FUNCTION HOSTNM(STRING)
        CHARACTER*(*) ,INTENT(IN) :: STRING
        END FUNCTION HOSTNM
      END INTERFACE
      CHARACTER(*),INTENT(OUT)  :: HOSTNAME
      INTEGER(4)                :: RC
!     **************************************************************************
      RC=HOSTNM(HOSTNAME)    
      IF(RC.NE.0)HOSTNAME='UNKNOWN'
      HOSTNAME='UNKNOWN'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_ABSOFT_FLUSH(NFIL)
!     **************************************************************************
!     **  FLUSHES THE FILE CONNECTED TO UNIT NFIL                             **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE ABSOFT COMPILER                          **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE FLUSH(LUNIT)
        INTEGER*4,INTENT(IN) :: LUNIT
        END SUBROUTINE FLUSH
      END INTERFACE
      INTEGER(4),INTENT(IN)       :: NFIL
      INTEGER*4                   :: LUNIT
!     *********************************************************************
      LUNIT=NFIL
      CALL FLUSH(LUNIT)
      RETURN
      END
#ELIF DEFINED(CPPVAR_COMPILER_XLF)
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_XLF_ETIME(USRTIME,SYSTIME)
!     **************************************************************************
!     **  RETURNS THE USER AND SYSTEM ELAPSED TIME OF THE CURRENT PROCESS     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE XLF COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
      TYPE TB_TYPE
        SEQUENCE
        REAL(4) :: USRTIME
        REAL(4) :: SYSTIME
      END TYPE TB_TYPE
!     **  THE ETIME_ FUNCTION SETS THE USER-ELAPSED TIME AND SYSTEM-ELAPSED   **
!     **  TIME IN ETIME_STRUCT SINCE THE START OF THE EXECUTION OF A PROCESS. **
!     **  THE RETURNED VALUE, ELAPSED, IS THE SUM OF THE USER-ELAPSED TIME    **
!     **  AND THE SYSTEM-ELAPSED TIME. THE RESOLUTION FOR ALL TIMING IS 1/100 **
!     **  OF A SECOND. THE OUTPUT APPEARS IN UNITS OF SECONDS.                **
      INTERFACE 
        REAL*4 FUNCTION ETIME(TARRAY)
        REAL*4,INTENT(OUT) :: TARRAY(2)
        END FUNCTION ETIME
      END INTERFACE
      REAL(8),INTENT(OUT) :: USRTIME
      REAL(8),INTENT(OUT) :: SYSTIME
      TYPE(TB_TYPE)       :: ETIME_STRUCT
      REAL*4              :: RESULT      
!     **************************************************************************
      RESULT=ETIME_(ETIME_STRUCT)
      USRTIME=ETIME_STRUCT%USRTIME
      SYSTIME=ETIME_STRUCT%SYSTIME
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_XLF_GETARG(IPOS,ARG)
!     **************************************************************************
!     **  RETURNS THE VALUE OF THE I-TH COMMAND LINE ARGUMENT                 **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE XLF COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE GETARG(POS,VALUE)
        INTEGER(4)       ,INTENT(IN) :: POS
        CHARACTER(LEN=*),INTENT(OUT):: VALUE
        END SUBROUTINE GETARG
      END INTERFACE
      INTEGER(4),INTENT(IN)         :: IPOS
      CHARACTER(LEN=*) ,INTENT(OUT) :: ARG
!     **************************************************************************
      CALL GETARG(IPOS,ARG)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_XLF_NARGS(NARGS)
!     **************************************************************************
!     **  RETURNS THE NUMBER OF COMMAND-LINE ARGUMENTS OF THE MAIN ROUTINE    **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE XLF COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        INTEGER(4) FUNCTION IARGC()
        END FUNCTION IARGC
      END INTERFACE
      INTEGER(4),INTENT(OUT) :: NARGS
!     **************************************************************************
      NARGS=IARGC()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_XLF_SYSTEM(COMMAND,rc)
!     **************************************************************************
!     **  ISSUES A SHELL COMMAND                                              **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE XLF COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        INTEGER*4 FUNCTION SYSTEM(STRING)
        CHARACTER*(*) ,INTENT(IN) :: STRING
        END FUNCTION SYSTEM
      END INTERFACE
      CHARACTER(*),INTENT(IN)  :: COMMAND
      INTEGER(4)  ,intent(out) :: RC
!     **************************************************************************
      RC=SYSTEM(COMMAND)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_XLF_GETHOSTNAME(HOSTNAME)
!     **************************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE XLF COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        INTEGER(4) FUNCTION HOSTNM_(STRING)
        CHARACTER(32) ,INTENT(IN) :: STRING
        END FUNCTION HOSTNM
      END INTERFACE
      CHARACTER(*),INTENT(OUT)  :: HOSTNAME
      CHARACTER(32)             :: HOSTNAME1
      INTEGER(4)                :: RC
!     **************************************************************************
      RC=HOSTNM_(HOSTNAME1)    
      HOSTNAME=HOSTNAME1
      IF(RC.NE.0)HOSTNAME='UNKNOWN'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_XLF_FLUSH(NFIL)
!     **************************************************************************
!     **  FLUSHES THE FILE CONNECTED TO UNIT NFIL                             **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE XLF COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE FLUSH_(LUNIT)
        INTEGER(4),INTENT(IN) :: LUNIT
        END SUBROUTINE FLUSH_
      END INTERFACE
      INTEGER(4),INTENT(IN)       :: NFIL
      INTEGER(4)                  :: LUNIT
!     *********************************************************************
      LUNIT=NFIL
      CALL FLUSH_(LUNIT)
      RETURN
      END
#ELIF DEFINED(CPPVAR_COMPILER_PGI)
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_PGI_ETIME(USRTIME,SYSTIME)
!     **************************************************************************
!     **  RETURNS THE USER AND SYSTEM ELAPSED TIME OF THE CURRENT PROCESS     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE PGI COMPILER                             **
!     **  THE MODULE DFPORT IS SUPPLIED BY THE PGI COMPILER                   **
!     **************************************************************************
!      USE DFPORT
      IMPLICIT NONE
!     **  THE ETIME_ FUNCTION SETS THE USER-ELAPSED TIME AND SYSTEM-ELAPSED   **
!     **  TIME IN ETIME_STRUCT SINCE THE START OF THE EXECUTION OF A PROCESS. **
!     **  THE RETURNED VALUE, ELAPSED, IS THE SUM OF THE USER-ELAPSED TIME    **
!     **  AND THE SYSTEM-ELAPSED TIME. THE RESOLUTION FOR ALL TIMING IS 1/100 **
!     **  OF A SECOND. THE OUTPUT APPEARS IN UNITS OF SECONDS.                **
      INTERFACE 
        REAL FUNCTION ETIME(TARRAY)
        REAL,INTENT(OUT) :: TARRAY(2)
        END FUNCTION ETIME
      END INTERFACE
      REAL(8),INTENT(OUT) :: USRTIME
      REAL(8),INTENT(OUT) :: SYSTIME
      REAL                :: RESULT      
      REAL                :: TARRAY(2)   
!     **************************************************************************
      RESULT=ETIME(TARRAY)
      USRTIME=TARRAY(1)
      SYSTIME=TARRAY(2)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_PGI_GETARG(IPOS,ARG)
!     **************************************************************************
!     **  RETURNS THE VALUE OF THE I-TH COMMAND LINE ARGUMENT                 **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE PGI COMPILER                             **
!     **  THE MODULE DFPORT IS SUPPLIED BY THE PGI COMPILER                   **
!     **************************************************************************
!      USE DFPORT
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE GETARG(POS,VALUE)
        INTEGER(4)       ,INTENT(IN) :: POS
        CHARACTER(LEN=*),INTENT(OUT):: VALUE
        END SUBROUTINE GETARG
      END INTERFACE
      INTEGER(4),INTENT(IN)         :: IPOS
      CHARACTER(LEN=*) ,INTENT(OUT) :: ARG
!     **************************************************************************
      CALL GETARG(IPOS,ARG)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_PGI_NARGS(NARGS)
!     **************************************************************************
!     **  RETURNS THE NUMBER OF COMMAND-LINE ARGUMENTS OF THE MAIN ROUTINE    **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE PGI COMPILER                             **
!     **  THE MODULE DFPORT IS SUPPLIED BY THE PGI COMPILER                   **
!     **************************************************************************
!      USE DFPORT
      IMPLICIT NONE
      INTERFACE 
        INTEGER(4) FUNCTION IARGC()
        END FUNCTION IARGC
      END INTERFACE
      INTEGER(4),INTENT(OUT) :: NARGS
!     **************************************************************************
      NARGS=IARGC()
PRINT*,'NARGS ',NARGS,IARGC()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_PGI_SYSTEM(COMMAND,rc)
!     **************************************************************************
!     **  ISSUES A SHELL COMMAND                                              **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE PGI COMPILER                             **
!     **  THE MODULE DFPORT IS SUPPLIED BY THE PGI COMPILER                   **
!     **************************************************************************
!      USE DFPORT
      IMPLICIT NONE
      INTERFACE 
        INTEGER FUNCTION SYSTEM(STRING)
        CHARACTER(*) ,INTENT(IN) :: STRING
        END FUNCTION SYSTEM
      END INTERFACE
      CHARACTER(*),INTENT(IN) :: COMMAND
      INTEGER(4)  ,intent(out):: RC
!     **************************************************************************
      RC=SYSTEM(COMMAND)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_PGI_GETHOSTNAME(HOSTNAME)
!     **************************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE PGI COMPILER                             **
!     **  THE MODULE DFPORT IS SUPPLIED BY THE PGI COMPILER                   **
!     **************************************************************************
!      USE DFPORT
      IMPLICIT NONE
      INTERFACE 
        INTEGER FUNCTION HOSTNM(STRING)
        CHARACTER*(*) ,INTENT(IN) :: STRING
        END FUNCTION HOSTNM
      END INTERFACE
      CHARACTER(*),INTENT(OUT)  :: HOSTNAME
      INTEGER                   :: RC
!     **************************************************************************
      RC=HOSTNM(HOSTNAME)    
      IF(RC.NE.0)HOSTNAME='UNKNOWN'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_PGI_FLUSH(NFIL)
!     **************************************************************************
!     **  FLUSHES THE FILE CONNECTED TO UNIT NFIL                             **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE PGI COMPILER                             **
!     **  THE MODULE DFPORT IS SUPPLIED BY THE PGI COMPILER                   **
!     **************************************************************************
!      USE DFPORT
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE FLUSH(LUNIT)
        INTEGER,INTENT(IN) :: LUNIT
        END SUBROUTINE FLUSH
      END INTERFACE
      INTEGER(4),INTENT(IN)       :: NFIL
      INTEGER                     :: LUNIT
!     **************************************************************************
      LUNIT=NFIL
      CALL FLUSH(LUNIT)
      RETURN
      END
#ELIF DEFINED(CPPVAR_COMPILER_PATHSCALE)
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_PATHSCALE_ETIME(USRTIME,SYSTIME)
!     **************************************************************************
!     **  RETURNS THE USER AND SYSTEM ELAPSED TIME OF THE CURRENT PROCESS     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE PGI COMPILER                             **
!     **  THE MODULE DFPORT IS SUPPLIED BY THE PGI COMPILER                   **
!     **************************************************************************
!      USE DFPORT
      IMPLICIT NONE
!     **  THE ETIME_ FUNCTION SETS THE USER-ELAPSED TIME AND SYSTEM-ELAPSED   **
!     **  TIME IN ETIME_STRUCT SINCE THE START OF THE EXECUTION OF A PROCESS. **
!     **  THE RETURNED VALUE, ELAPSED, IS THE SUM OF THE USER-ELAPSED TIME    **
!     **  AND THE SYSTEM-ELAPSED TIME. THE RESOLUTION FOR ALL TIMING IS 1/100 **
!     **  OF A SECOND. THE OUTPUT APPEARS IN UNITS OF SECONDS.                **
      INTERFACE 
        REAL FUNCTION ETIME(TARRAY)
        REAL,INTENT(OUT) :: TARRAY(2)
        END FUNCTION ETIME
      END INTERFACE
      REAL(8),INTENT(OUT) :: USRTIME
      REAL(8),INTENT(OUT) :: SYSTIME
      REAL                :: RESULT      
      REAL                :: TARRAY(2)   
!     **************************************************************************
      RESULT=ETIME(TARRAY)
      USRTIME=TARRAY(1)
      SYSTIME=TARRAY(2)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_PATHSCALE_GETARG(IPOS,ARG)
!     **************************************************************************
!     **  RETURNS THE VALUE OF THE I-TH COMMAND LINE ARGUMENT                 **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE PGI COMPILER                             **
!     **  THE MODULE DFPORT IS SUPPLIED BY THE PGI COMPILER                   **
!     **************************************************************************
!      USE DFPORT
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE GETARG(POS,VALUE)
        INTEGER(4)       ,INTENT(IN) :: POS
        CHARACTER(LEN=*),INTENT(OUT):: VALUE
        END SUBROUTINE GETARG
      END INTERFACE
      INTEGER(4),INTENT(IN)         :: IPOS
      CHARACTER(LEN=*) ,INTENT(OUT) :: ARG
!     **************************************************************************
      CALL GETARG(IPOS,ARG)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_PATHSCALE_NARGS(NARGS)
!     **************************************************************************
!     **  RETURNS THE NUMBER OF COMMAND-LINE ARGUMENTS OF THE MAIN ROUTINE    **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE PGI COMPILER                             **
!     **  THE MODULE DFPORT IS SUPPLIED BY THE PGI COMPILER                   **
!     **************************************************************************
!      USE DFPORT
      IMPLICIT NONE
      INTERFACE 
        INTEGER(4) FUNCTION IARGC()
        END FUNCTION IARGC
      END INTERFACE
      INTEGER(4),INTENT(OUT) :: NARGS
!     **************************************************************************
      NARGS=IARGC()
PRINT*,'NARGS ',NARGS,IARGC()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_PATHSCALE_SYSTEM(COMMAND,rc)
!     **************************************************************************
!     **  ISSUES A SHELL COMMAND                                              **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE PGI COMPILER                             **
!     **  THE MODULE DFPORT IS SUPPLIED BY THE PGI COMPILER                   **
!     **************************************************************************
!      USE DFPORT
      IMPLICIT NONE
      INTERFACE 
        INTEGER FUNCTION SYSTEM(STRING)
        CHARACTER(*) ,INTENT(IN) :: STRING
        END FUNCTION SYSTEM
      END INTERFACE
      CHARACTER(*),INTENT(IN) :: COMMAND
      INTEGER(4)  ,intent(out):: RC
!     **************************************************************************
      RC=SYSTEM(COMMAND)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_PATHSCALE_GETHOSTNAME(HOSTNAME)
!     **************************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE PGI COMPILER                             **
!     **  THE MODULE DFPORT IS SUPPLIED BY THE PGI COMPILER                   **
!     **************************************************************************
!      USE DFPORT
      IMPLICIT NONE
      INTERFACE 
        INTEGER FUNCTION HOSTNM(STRING)
        CHARACTER*(*) ,INTENT(IN) :: STRING
        END FUNCTION HOSTNM
      END INTERFACE
      CHARACTER(*),INTENT(OUT)  :: HOSTNAME
      INTEGER                   :: RC
!     **************************************************************************
      RC=HOSTNM(HOSTNAME)    
      IF(RC.NE.0)HOSTNAME='UNKNOWN'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_PATHSCALE_FLUSH(NFIL)
!     **************************************************************************
!     **  FLUSHES THE FILE CONNECTED TO UNIT NFIL                             **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE PGI COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE FLUSH(LUNIT)
        INTEGER,INTENT(IN) :: LUNIT
        END SUBROUTINE FLUSH
      END INTERFACE
      INTEGER(4),INTENT(IN)       :: NFIL
      INTEGER                     :: LUNIT
!     **************************************************************************
      LUNIT=NFIL
      CALL FLUSH(LUNIT)
      RETURN
      END
#ELIF DEFINED(CPPVAR_COMPILER_SUN)
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_SUN_ETIME(USRTIME,SYSTIME)
!     **************************************************************************
!     **  RETURNS THE USER AND SYSTEM ELAPSED TIME OF THE CURRENT PROCESS     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE SUN COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
!     **  THE ETIME_ FUNCTION SETS THE USER-ELAPSED TIME AND SYSTEM-ELAPSED   **
!     **  TIME IN ETIME_STRUCT SINCE THE START OF THE EXECUTION OF A PROCESS. **
!     **  THE RETURNED VALUE, ELAPSED, IS THE SUM OF THE USER-ELAPSED TIME    **
!     **  AND THE SYSTEM-ELAPSED TIME. THE RESOLUTION FOR ALL TIMING IS 1/100 **
!     **  OF A SECOND. THE OUTPUT APPEARS IN UNITS OF SECONDS.                **
      INTERFACE 
        REAL FUNCTION ETIME(TARRAY)
        REAL,INTENT(OUT) :: TARRAY(2)
        END FUNCTION ETIME
      END INTERFACE
      REAL(8),INTENT(OUT) :: USRTIME
      REAL(8),INTENT(OUT) :: SYSTIME
      REAL                :: RESULT      
      REAL                :: TARRAY(2)   
!     **************************************************************************
      RESULT=ETIME(TARRAY)
      USRTIME=TARRAY(1)
      SYSTIME=TARRAY(2)
      RETURN
      END SUBROUTINE LIB_SUN_ETIME
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_SUN_GETARG(IPOS,ARG)
!     **************************************************************************
!     **  RETURNS THE VALUE OF THE I-TH COMMAND LINE ARGUMENT                 **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE SUN COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE GETARG(POS,VALUE)
        INTEGER(4)      ,INTENT(IN) :: POS
        CHARACTER(LEN=*),INTENT(OUT):: VALUE
        END SUBROUTINE GETARG
      END INTERFACE
      INTEGER(4),INTENT(IN)         :: IPOS
      CHARACTER(LEN=*) ,INTENT(OUT) :: ARG
!     **************************************************************************
      CALL GETARG(IPOS,ARG)
      RETURN
      END SUBROUTINE LIB_SUN_GETARG
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_SUN_NARGS(NARGS)
!     **************************************************************************
!     **  RETURNS THE NUMBER OF COMMAND-LINE ARGUMENTS OF THE MAIN ROUTINE    **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE SUN COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        INTEGER(4) FUNCTION IARGC()
        END FUNCTION IARGC
      END INTERFACE
      INTEGER(4),INTENT(OUT) :: NARGS
!     **************************************************************************
      NARGS=IARGC()
PRINT*,'NARGS ',NARGS,IARGC()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_SUN_SYSTEM(COMMAND,rc)
!     **************************************************************************
!     **  ISSUES A SHELL COMMAND                                              **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE SUN COMPILER                             **
!     **************************************************************************
!      USE DFPORT
      IMPLICIT NONE
      INTERFACE 
        INTEGER FUNCTION SYSTEM(STRING)
        CHARACTER(*) ,INTENT(IN) :: STRING
        END FUNCTION SYSTEM
      END INTERFACE
      CHARACTER(*),INTENT(IN) :: COMMAND
      INTEGER(4)  ,intent(out):: RC
!     **************************************************************************
      RC=SYSTEM(COMMAND)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_SUN_GETHOSTNAME(HOSTNAME)
!     **************************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                     **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE SUN COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        INTEGER FUNCTION HOSTNM(STRING)
        CHARACTER*(*) ,INTENT(IN) :: STRING
        END FUNCTION HOSTNM
      END INTERFACE
      CHARACTER(*),INTENT(OUT)  :: HOSTNAME
      INTEGER                   :: RC
!     **************************************************************************
      RC=HOSTNM(HOSTNAME)    
      IF(RC.NE.0)HOSTNAME='UNKNOWN'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_SUN_FLUSH(NFIL)
!     **************************************************************************
!     **  FLUSHES THE FILE CONNECTED TO UNIT NFIL                             **
!     **                                                                      **
!     **  SPECIFIC INTERFACE FOR THE SUN COMPILER                             **
!     **************************************************************************
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE FLUSH(LUNIT)
        INTEGER,INTENT(IN) :: LUNIT
        END SUBROUTINE FLUSH
      END INTERFACE
      INTEGER(4),INTENT(IN)       :: NFIL
      INTEGER(4)                     :: LUNIT
!     **************************************************************************
      LUNIT=NFIL
      CALL FLUSH(LUNIT)
      RETURN
      END
#ENDIF
! 
!*******************************************************************************
!*******************************************************************************
!****                                                                      *****
!****              INTERFACES TO THE LAPACK SUBROUTINES                    *****
!****                                                                      *****
!****  THE LINEAR ALGEBRA PACKAGE (LAPACK) ADRESSES THE FOLLOWING PROBLEMS *****
!****  LINEAR EQUATIONS                                                    *****
!****  LINEAR LEAST SQUARES                                                *****
!****  GENERALIZED LINEAR LEAST SQUARES                                    *****
!****  SYMMETRIC EIGENPROBLEMS                                             *****
!****  NONSYMMETRIC EIGENPROBLEMS                                          *****
!****  SINGULAR VALUE DECOMPOSITION (SVD)                                  *****
!****  GENERALIZED SYMMETRIC DEFINITE EIGENPROBLEMS                        *****
!****  GENERALIZED NONSYMMETRIC EIGENPROBLEMS                              *****
!****  GENERALIZED SINGULAR VALUE DECOMPOSITION                            *****
!****                                                                      *****
!****  HERE WE COVER ONLY A SUBSET OF THESE                                *****
!****    MATRIX INVERSION                                                  *****
!****    LINEAR EQUATIONS                                                  *****
!****    SYMMETRIC EIGENVALUE PROBLEMS                                     *****
!****                                                                      *****
!****                                                                      *****
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$INVERTR8(N,A,AINV)
!     **************************************************************************
!     **  INVERTS THE REAL, SQUARE MATRIX A                                   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N         ! DIMENSION OF THE MATRIX
      REAL(8)   ,INTENT(IN) :: A(N,N)    ! MATRIX TO BE INVERTED
      REAL(8)   ,INTENT(OUT):: AINV(N,N) ! INVERTED MATRIX
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)   ,ALLOCATABLE:: RES(:,:)
      REAL(8)               :: DEV
      INTEGER(4)            :: I
!     **************************************************************************
!      == GENERAL MATRIX INVERSE ===============================================
#IF DEFINED(CPPVAR_LAPACK_ESSL)
       CALL LIB_ESSL_DGEICD(N,A,AINV)
#ELSE 
      CALL LIB_LAPACK_DGETRI(N,A,AINV)
#ENDIF
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(RES(N,N))
        RES=MATMUL(A,AINV)
        DO I=1,N
          RES(I,I)=RES(I,I)-1.D0
        ENDDO
        DEV=MAXVAL(ABS(RES))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$INVERTR8')
        END IF
        DEALLOCATE(RES)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$INVERTC8(N,A,AINV)
!     **************************************************************************
!     **  INVERTS THE complex, SQUARE MATRIX A                                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N         ! DIMENSION OF THE MATRIX
      complex(8),INTENT(IN) :: A(N,N)    ! MATRIX TO BE INVERTED
      complex(8),INTENT(OUT):: AINV(N,N) ! INVERTED MATRIX
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      complex(8),ALLOCATABLE:: RES(:,:)
      REAL(8)               :: DEV
      INTEGER(4)            :: I
!     **************************************************************************
!      == GENERAL MATRIX INVERSE ===============================================
#IF DEFINED(CPPVAR_LAPACK_ESSL)
      CALL LIB_LAPACK_ZGETRI(N,A,AINV)
!      CALL LIB_ESSL_ZGETRI(N,A,AINV)
#ELSE 
      CALL LIB_LAPACK_ZGETRI(N,A,AINV)
#ENDIF
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(RES(N,N))
        RES=MATMUL(A,AINV)
        DO I=1,N
          RES(I,I)=RES(I,I)-(1.D0,0.d0)
        ENDDO
        DEV=MAXVAL(ABS(RES))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$INVERTR8')
        END IF
        DEALLOCATE(RES)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$PSEUDOINVERTR8(N,A,smin,AINV)
!     **************************************************************************
!     **  CONSTRUCTS PSEUDO INVERSE OF A USING SINGULAR VALUE DECOMPOSITION   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N         ! DIMENSION OF THE MATRIX
      REAL(8)   ,INTENT(IN) :: A(N,N)    ! MATRIX TO BE INVERTED
      REAL(8)   ,INTENT(IN) :: smin      ! min singular value to be considered
      REAL(8)   ,INTENT(OUT):: AINV(N,N) ! INVERTED MATRIX
      LOGICAL   ,PARAMETER  :: TTEST=.true.
      COMPLEX(8),ALLOCATABLE:: RES(:,:)
      REAL(8)               :: DEV
      INTEGER(4)            :: I
      REAL(8)               :: U(N,N)
      REAL(8)               :: V(N,N)
      REAL(8)               :: S(N)
!     **************************************************************************
#IF DEFINED(CPPVAR_LAPACK_ESSL)
      CALL ERROR$MSG('INTERFACE TO ESSL ROUTINE NOT IMPLEMENTED')
#ELSE 
      CALL LIB_LAPACK_DGESVD(N,N,A,U,S,V)
!print*,'s ',s
      DO I=1,N
        IF(ABS(S(i)).GT.SMIN) THEN
          S(i)=1.D0/S(i)
        ELSE
          S(i)=0.D0
        END IF
      ENDDO
      DO I=1,N
        U(:,I)=U(:,I)*S(I)
      ENDDO
      AINV(:,:)=transpose(MATMUL(U,V))
#ENDIF
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(RES(N,N))
        RES=MATMUL(A,AINV)
        DO I=1,N
          RES(I,I)=RES(I,I)-(1.D0,0.D0)
        ENDDO
        DEV=MAXVAL(ABS(RES))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$INVERTR8')
        END IF
        DEALLOCATE(RES)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$SVDR8(M,N,A,U,S,Vt)
!     **************************************************************************
!     **  SINGULAR VALUE DECOMPOSITION OF THE MXN MATRIX A                    **
!     **  a = u * s * vt
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: m         ! DIMENSION OF THE MATRIX
      INTEGER(4),INTENT(IN) :: N         ! DIMENSION OF THE MATRIX
      REAL(8)   ,INTENT(IN) :: A(m,n)    ! MATRIX TO BE INVERTED
      REAL(8)   ,INTENT(OUT):: u(m,m)    ! left hand orthogonal vectors
      REAL(8)   ,INTENT(OUT):: s(m)      ! singular vectors
      REAL(8)   ,INTENT(OUT):: vt(n,n)   ! transposed right-hand orthogonal vectors
!     **************************************************************************
#IF DEFINED(CPPVAR_LAPACK_ESSL)
      CALL ERROR$MSG('INTERFACE TO ESSL ROUTINE NOT IMPLEMENTED')
#ELSE 
      CALL LIB_LAPACK_DGESVD(m,N,A,U,S,Vt)
#ENDIF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$MATRIXSOLVER8(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  SOLVES THE LINEAR EQUATION SYSTEM AX=B                      **
!     **  WHERE A IS A(N,M)                                           **
!     **  IF A IS NOT SQUARE, THE EQUATION IS SOLVED IN A LEAST       **
!     **  SQUARE SENSE                                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(OUT):: X(M,NEQ)
      REAL(8)   ,INTENT(IN) :: B(N,NEQ)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: SVAR1,SVAR2
      INTEGER(4)            :: IPIV(N)
!     ******************************************************************
!
!     ==================================================================
!     == SOLVE SPECIAL CASE WITH DIMENSIONS 1 FIRST                   ==
!     ==================================================================
      IF(N.EQ.1.AND.M.EQ.1) THEN
        X(1,:)=B(1,:)/A(1,1)
        RETURN
      END IF
#IF DEFINED(CPPVAR_LAPACK_ESSL)
!
!     ==================================================================
!     == CALL ESSL DRIVER ROUTINES                                    ==
!     ==================================================================
      IF(N.EQ.M) THEN
        CALL LIB_ESSL_DGES(N,M,NEQ,A,X,B)
      ELSE IF(N.LT.M) THEN
        CALL LIB_ESSL_DGESVS(N,M,NEQ,A,X,B)
      ELSE
        CALL ERROR$MSG('SYSTEM OF EQUATIONS IS OVER DETERMINED')
        CALL ERROR$STOP(' LIB$MATRIXSOLVER8')
      END IF
#ELSE 
!
!     ==================================================================
!     == CALL LAPACK DRIVER ROUTINES                                  ==
!     ==================================================================
      IF(N.EQ.M) THEN
        CALL LIB_LAPACK_DGESV(N,M,NEQ,A,X,B)
      ELSE IF(N.LT.M) THEN
        CALL LIB_LAPACK_DGELSD(N,M,NEQ,A,X,B)
      ELSE
        CALL ERROR$MSG('SYSTEM OF EQUATIONS IS OVER DETERMINED')
        CALL ERROR$STOP(' LIB$MATRIXSOLVER8')
      END IF
#ENDIF
!
!     ==================================================================
!     ==  TEST                                                        ==
!     ==================================================================
      IF(TTEST) THEN
        SVAR1=MAXVAL(ABS(MATMUL(A,X)-B))
        SVAR2=MAXVAL(ABS(A))/MAXVAL(ABS(B))
        IF(SVAR1/SVAR2.GT.1.D-4) THEN
          CALL ERROR$R8VAL('SVAR1',SVAR1)
          CALL ERROR$R8VAL('SVAR2',SVAR2)
          CALL ERROR$R8VAL('SVAR1/SVAR2',SVAR1/SVAR2)
          CALL ERROR$STOP('LIB$MATRIXSOLVER8')
        END IF
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$MATRIXSOLVEC8(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  SOLVES THE LINEAR EQUATION SYSTEM AX=B                      **
!     **  WHERE A IS A(N,M)                                           **
!     **  IF A IS NOT SQUARE, THE EQUATION IS SOLVED IN A LEAST       **
!     **  SQUARE SENSE                                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(OUT):: X(M,NEQ)
      COMPLEX(8),INTENT(IN) :: B(N,NEQ)
!     ******************************************************************
      IF(N.EQ.1.AND.M.EQ.1) THEN
        X(1,:)=B(1,:)/A(1,1)
        RETURN
      END IF
#IF DEFINED(CPPVAR_LAPACK_ESSL)
!
!     ==================================================================
!     == CALL ESSL DRIVER ROUTINES                                    ==
!     ==================================================================
      IF(N.EQ.M) THEN
        CALL LIB_ESSL_ZGES(N,M,NEQ,A,X,B)
      ELSE IF(N.LT.M) THEN
        CALL LIB_ESSL_ZGESVS(N,M,NEQ,A,X,B) 
      ELSE
        CALL ERROR$MSG('SYSTEM OF EQUATIONS IS OVER DETERMINED')
        CALL ERROR$STOP(' LIB$MATRIXSOLVEC8')
      END IF
#ELSE 
      IF(N.EQ.M) THEN
        CALL LIB_LAPACK_ZGESV(N,M,NEQ,A,X,B)
      ELSE IF(N.LT.M) THEN
        CALL LIB_LAPACK_ZGELSD(N,M,NEQ,A,X,B)
      ELSE
        CALL ERROR$MSG('SYSTEM OF EQUATIONS IS OVER DETERMINED')
        CALL ERROR$STOP(' LIB$MATRIXSOLVEC8')
      END IF
#ENDIF
      RETURN
      END SUBROUTINE LIB$MATRIXSOLVEC8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$DIAGR8(N,H,E,U)
!     **************************************************************************
!     **                                                                      **
!     **  DIAGONALIZES THE REAL, SQUARE MATRIX H AFTER SYMMETRIZATION         **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                           **
!     **                                                                      **
!     **         U(K,I)*H(K,L)*U(L,J)=DELTA(I,J)*E(I)                         **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE              **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE                  **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES                     **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      REAL(8)   ,INTENT(OUT):: U(N,N)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: DEV
      REAL(8)  ,ALLOCATABLE :: EMAT(:,:)
      INTEGER(4)            :: I
!     **************************************************************************
!
!     ==========================================================================
!     == DIAGONALIZE                                                          ==
!     ==========================================================================
#IF DEFINED(CPPVAR_LAPACK_ESSL)
      CALL LIB_ESSL_DSPEV(N,H,E,U)
#ELSE
      CALL LIB_LAPACK_DSYEV(N,H,E,U)
#ENDIF
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(EMAT(N,N))
!       == TEST EIGENVALUE EQUATION ============================================
        EMAT(:,:)=0.D0
        DO I=1,N
          EMAT(I,I)=E(I)
        ENDDO
        DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(U,EMAT)))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$DIAGR8')
        END IF
!       == TEST ORTHONORMALITY OF EIGENVECTORS =================================
        EMAT=MATMUL(TRANSPOSE(U),U)
        DO I=1,N
          EMAT(I,I)=EMAT(I,I)-1.D0
        ENDDO
        DEV=MAXVAL(ABS(EMAT))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('ORTHONORMALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$DIAGR8')
        END IF
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$DIAGC8(N,H,E,U)
!     **************************************************************************
!     **                                                                      **
!     **  DIAGONALIZES THE HERMITEAN, SQUARE MATRIX H                         **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                           **
!     **                                                                      **
!     **      CONJG(U(K,I))*H(K,L)*U(L,J)=DELTA(I,J)*E(I)                     **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE              **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE                  **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES                     **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      COMPLEX(8),INTENT(OUT):: U(N,N)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: DEV
      COMPLEX(8),ALLOCATABLE:: EMAT(:,:)
      INTEGER(4)            :: I
!     **************************************************************************
      IF(N.EQ.1) THEN
        E(1)=real(H(1,1),kind=8)
        U(1,1)=(1.D0,0.D0)
        RETURN
      END IF
!
!     ==========================================================================
!     == DIAGONALIZE                                                          ==
!     ==========================================================================
#IF DEFINED(CPPVAR_LAPACK_ESSL)
      CALL LIB_ESSL_ZHPEV(N,H,E,U)
#ELSE
      CALL LIB_LAPACK_ZHEEV(N,H,E,U)
#ENDIF
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(EMAT(N,N))
!       == TEST EIGENVALUE EQUATION ============================================
        EMAT(:,:)=CMPLX(0.D0,0.D0)
        DO I=1,N
          EMAT(I,I)=CMPLX(E(I),0.D0)
        ENDDO
        DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(U,EMAT)))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$DIAGC8')
        END IF
!       == TEST ORTHONORMALITY OF EIGENVECTORS =================================
        EMAT=MATMUL(TRANSPOSE(CONJG(U)),U)
        DO I=1,N
          EMAT(I,I)=EMAT(I,I)-CMPLX(1.D0,0.D0)
        ENDDO
        DEV=MAXVAL(ABS(EMAT))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('ORTHONORMALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$DIAGC8')
        END IF
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$GENERALEIGENVALUER8(N,H,S,E,U)
!     **************************************************************************
!     **                                                                      **
!     ** SOLVES THE GENERALIZED, REAL, SYMMETRIC EIGENVALUE PROBLEM           **
!     **                                                                      **
!     **      H*U = S*U*E                                                     **
!     **                                                                      **
!     ** WITH EIGENVECTORS U THAT ARE ORTHOGONAL IN THE SENSE                 **
!     **                                                                      **
!     **      U^T*S*U=IDENTITY                                                **
!     **                                                                      **
!     ** REMARK: H AND S MUST BE SYMMETRIC                                    **
!     **         S MUST BE POSITIVE DEFINITE                                  **
!     **         EIGENVECTORS ARE ORTHONORMAL IN THE SENSE                    **
!     **             MATMUL(TRANSPOSE(U),MATMUL(S,U))=IDENTITY                **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)    ! HAMILTON MATRIX
      REAL(8)   ,INTENT(IN) :: S(N,N)    ! OVERLAP MATRIX
      REAL(8)   ,INTENT(OUT):: E(N)      ! EIGENVALUES
      REAL(8)   ,INTENT(OUT):: U(N,N)    ! EIGENVECTORS
      REAL(8)               :: B(N,N)    ! COPY OF OVERLAP MATRIX
      LOGICAL               :: TSYM
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE. ! IF TRUE TEST RESULT
      REAL(8)               :: DEV       ! DEVIATION
      INTEGER               :: I 
!     **************************************************************************
!
!     ==========================================================================
!     == TAKE CARE OF TRIVIAL CASES                                           ==
!     ==========================================================================
      IF(N.EQ.1) THEN
        E(1)=H(1,1)/S(1,1)
        U(1,1)=1.D0/SQRT(S(1,1))
        RETURN
      ELSE IF(N.EQ.0) THEN
        RETURN
      END IF
!
!     ==========================================================================
!     == TEST IF INPUT MATRICES ARE SYMMETRIC                                 ==
!     ==========================================================================
      DEV=MAXVAL(ABS(H-TRANSPOSE(H)))
      TSYM=(DEV.LT.1.D-5)
      IF(.NOT.TSYM) THEN
        CALL ERROR$MSG('HAMILTONIAN NOT SYMMETRIC')
        CALL ERROR$R8VAL('DEV',DEV)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
      END IF
      DEV=MAXVAL(ABS(S-TRANSPOSE(S)))
      TSYM=TSYM.AND.(DEV.LT.1.D-5)
      IF(.NOT.TSYM) THEN
        CALL ERROR$MSG('OVERLAP MATRIX NOT SYMMETRIC')
        CALL ERROR$R8VAL('DEV',DEV)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
      END IF
!
!     ==========================================================================
!     == DIAGONALIZE                                                          ==
!     ==========================================================================
#IF DEFINED(CPPVAR_LAPACK_ESSL)
      IF(TSYM) THEN
        CALL LIB_ESSL_DSYGV(N,H,S,E,U)
      ELSE
        CALL LIB_ESSL_DGEGV(N,H,S,E,U)
      END IF
#ELSE
      CALL LIB_LAPACK_DSYGV(N,H,S,E,U)
#ENDIF
!
!     ==========================================================================
!     == TEST RESULT OF THE ROUTINE                                           ==
!     ==========================================================================
      IF(TTEST) THEN
!       == CHECK ORTHONORMALITY OF EIGENVECTORS ================================
        B=MATMUL(TRANSPOSE(U),MATMUL(S,U))
        DO I=1,N
          B(I,I)=B(I,I)-1.D0
        ENDDO
        DEV=SUM(ABS(B))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('EIGENSTATES NOT ORTHONORMAL')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
        END IF
!       == CHECK EIGENVALUE PROBLEM ============================================
        DO I=1,N
          DEV=SUM(ABS(MATMUL(H-E(I)*S,U(:,I))))
          IF(DEV.GT.1.D-7) THEN
            CALL ERROR$MSG('EIGENVALUE PROBLEM FAILED')
            CALL ERROR$R8VAL('DEV',DEV)
            CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
          END IF
        ENDDO
      END IF
      RETURN
      END SUBROUTINE LIB$GENERALEIGENVALUER8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$GENERALEIGENVALUEC8(N,H,S,E,U)
!     **************************************************************************
!     **                                                                      **
!     ** SOLVES THE GENERALIZED, complex NON-SYMMETRIC EIGENVALUE PROBLEM        **
!     **      [H(:,:)-E(I)*S(:,:)]*U(:,I)=0                                   **
!     **                                                                      **
!     ** REMARK: H AND S MUST BE HERMITEANC                                   **
!     **         S MUST BE POSITIVE DEFINITE                                  **
!     **         EIGENVECTORS ARE ORTHONORMAL IN THE SENSE                    **
!     **             MATMUL(TRANSPOSE(U),MATMUL(S,U))=IDENTITY                **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)    ! HAMITON MATRIX
      COMPLEX(8),INTENT(IN) :: S(N,N)    ! OVERLAP MATRIX
      REAL(8)   ,INTENT(OUT):: E(N)      ! EIGENVALUES
      COMPLEX(8),INTENT(OUT):: U(N,N)  ! EIGENVECTORS
      COMPLEX(8)            :: S1(N,N)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: DEV
      INTEGER               :: I
!     **************************************************************************
!
!     ==========================================================================
!     == TEST IF INPUT MATRICES ARE SYMMETRIC                                 ==
!     ==========================================================================
      DEV=SUM(ABS(H-TRANSPOSE(CONJG(H))))
      IF(DEV.GT.1.D-8) THEN
        CALL ERROR$MSG('HAMILTON MATRIX NOT HERMITEAN')
        CALL ERROR$R8VAL('DEV',DEV)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
      END IF
      DEV=SUM(ABS(S-TRANSPOSE(CONJG(S))))
      IF(DEV.GT.1.D-8) THEN
        CALL ERROR$MSG('OVERLAP MATRIX NOT HERMITEAN')
        CALL ERROR$R8VAL('DEV',DEV)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
      END IF
!
!     ========================================================================
!     == CALL LAPACK ROUTINE                                                ==
!     ========================================================================
#IF DEFINED(CPPVAR_LAPACK_ESSL)
      CALL ERROR$MSG('GENERALEIGENVALUE PROBLEM NOT IMPLEMENTED FOR ESSL')
      CALL ERROR$STOP('GENERALEIGENVALUER8')  
#ELSE
      CALL LIB_LAPACK_ZHEGV(N,H,S,E,U)
#ENDIF
!
!     ========================================================================
!     == TEST RESULT OF THE ROUTINE                                         ==
!     ========================================================================
      IF(TTEST) THEN
        S1=MATMUL(TRANSPOSE(CONJG(U)),MATMUL(S,U))
        DEV=0.D0
        DO I=1,N
          S1(I,I)=S1(I,I)-(1.D0,0.D0)
          DEV=MAX(DEV,MAXVAL(ABS(MATMUL(H-E(I)*S,U(:,I)))))
        ENDDO
        IF(DEV.GT.1.D-6) THEN
          CALL ERROR$MSG('GENERAL EIGENVALUE TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
        END IF
        DEV=SUM(ABS(S1))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('EIGENSTATES NOT ORTHONORMAL')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
        END IF
      END IF
      RETURN
      END
#IF DEFINED(CPPVAR_LAPACK_ESSL)
!***********************************************************************
!***********************************************************************
!****                                                               ****
!****  DRIVER ROUTINES FOR ESSL CORRESPONDING TO LAPACK             ****
!****                                                               ****
!***********************************************************************
!***********************************************************************
!
!     ..................................................................
      SUBROUTINE LIB_ESSL_DGEICD(N,A,AINV)
!     ******************************************************************
!     **                                                              **
!     **  INVERTS THE REAL, SQUARE MATRIX A                           **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: DGEICD                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: A(N,N)
      REAL(8)   ,INTENT(OUT):: AINV(N,N)
      INTEGER(4)            :: NAUX
      REAL(8)               :: AUX(100*N)
      REAL(8)               :: RCOND
      REAL(8)               :: DET(2)
!     ******************************************************************
      NAUX=100*N
      AINV(:,:)=A(:,:)
      CALL DGEICD(AINV,N,N,0,RCOND,DET,AUX,NAUX) !ESSL
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB_ESSL_DGES(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  SOLVES THE LINEAR EQUATION SYSTEM AX=B                      **
!     **  WHERE A IS A SQUARE MATRIX                                  **
!     **                                                              **
!     **  REMARK: WHILE DGES CAN ALSO HANDLE NON-SQUARE MATRICES      **
!     **     WE USE IT HERE FOR THIS PURPOSE ONLY.                    **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(OUT):: X(M,NEQ)
      REAL(8)   ,INTENT(IN) :: B(N,NEQ)
      REAL(8)               :: AFACT(N,N)
      INTEGER(4)            :: IPVT(N)
      INTEGER(4)            :: INFO
      INTEGER(4)            :: I
      REAL(8)               :: DEV
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
!     ******************************************************************
      IF(N.NE.M) THEN  !MATRIX FACTORIZATION
        CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
        CALL ERROR$STOP('LIB_ESSL_DGES')
      END IF
      AFACT=A
!     == GENERAL MATRIX FACTORIZATION =================================
      CALL DGEF(AFACT,N,N,IPVT)
      X=B
      DO I=1,NEQ
        CALL DGES(AFACT,N,N,IPVT,X(:,I),0)
      ENDDO
!     =================================================================
!     ==  TEST                                                       ==
!     =================================================================
      IF(TTEST) THEN
        DEV=MAXVAL(ABS(MATMUL(A,X)-B))
        WRITE(*,*)'MAX ERROR OF LIB$MATRIXSOLVE ',DEV
      END IF
      RETURN
      END SUBROUTINE LIB_ESSL_DGES
!
!     ..................................................................
      SUBROUTINE LIB_ESSL_DGESVS(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  SOLVES THE LINEAR EQUATION SYSTEM AX=B                      **
!     **  WHERE A IS A(N,M)                                           **
!     **  IF A IS NOT SQUARE, THE EQUATION IS SOLVED IN A LEAST       **
!     **  SQUARE SENSE                                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(OUT):: X(M,NEQ)
      REAL(8)   ,INTENT(IN) :: B(N,NEQ)
      REAL(8)   ,ALLOCATABLE:: AFACT(:,:)
      REAL(8)   ,ALLOCATABLE:: BFACT(:,:)
      REAL(8)   ,ALLOCATABLE:: AUX(:)
      REAL(8)   ,ALLOCATABLE:: S(:)
      INTEGER(4)            :: NAUX
      INTEGER(4)            :: NM
      INTEGER(4)            :: IRANK
      REAL(8)               :: TAU=1.D-6
      INTEGER(4)            :: INFO
      INTEGER(4)            :: I
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      INTEGER(4)            :: M1,N1
!     ******************************************************************
!     ===========================================================
!     == SINGULAR VALUE DECOMPOSION WITH ESSL                  ==
!     ===========================================================
      M1=N
      N1=M
      NM=MAX(1,MAX(N,M))
      ALLOCATE(AFACT(NM,M))
      ALLOCATE(BFACT(NM,NEQ))
      AFACT(1:N,:)=A
      AFACT(N+1:,:)=0.D0
      BFACT(:N,:)=B
      NAUX=2*N1+MAX(NM,NEQ)
      ALLOCATE(S(N1))    
      ALLOCATE(AUX(NAUX))
!     ==  SINGULAR VALUES OF A GENERAL MATRIX =========================
      CALL DGESVF(2,AFACT,NM,BFACT,NM,NEQ,S,M1,N1,AUX,NAUX) !->ESSL
      DEALLOCATE(AUX)
!     ==  LINEAR LEAST SQUARES SOLUTION FOR A GENERAL MATRIX ==========
!     ==  USING THE SINGULAR VALUE DECOMPOSITION
      CALL DGESVS(AFACT,NM,BFACT,NM,NEQ,S,X,M,M1,N1,TAU) !->ESSL
      DEALLOCATE(S)
      DEALLOCATE(BFACT)
      DEALLOCATE(AFACT)
!     =================================================================
!     ==  TEST                                                       ==
!     =================================================================
      IF(TTEST) THEN
        ALLOCATE(AUX(1))
        AUX(:)=MAXVAL(ABS(MATMUL(A,X)-B))
        WRITE(*,*)'MAX ERROR OF LIB$MATRIXSOLVE ',AUX(1)
        DEALLOCATE(AUX)
      END IF
      RETURN
      END SUBROUTINE LIB_ESSL_DGESVS
!
!     ..................................................................
      SUBROUTINE LIB_ESSL_ZGES(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  DRIVER ROUTINE FOR ESSL                                     **
!     **                                                              **
!     **  SOLVES THE COMPLEX LINEAR EQUATION SYSTEM AX=B              **
!     **  WHERE A IS A SQUARE MATRIX                                  **
!     **                                                              **
!     **  REMARK: WHILE DGES CAN ALSO HANDLE NON-SQUARE MATRICES      **
!     **     WE USE IT HERE FOR THIS PURPOSE ONLY.                    **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(OUT):: X(M,NEQ)
      COMPLEX(8),INTENT(IN) :: B(N,NEQ)
      COMPLEX(8)            :: AFACT(N,N)
      INTEGER(4)            :: IPVT(N)
      INTEGER(4)            :: INFO
      INTEGER(4)            :: I
      REAL(8)               :: DEV
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
!     ******************************************************************
      IF(N.NE.M) THEN  !MATRIX FACTORIZATION
        CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
        CALL ERROR$STOP('LIB_ESSL_DGES')
      END IF
      AFACT=A
!     == GENERAL MATRIX FACTORIZATION =================================
      CALL ZGEF(AFACT,N,N,IPVT)
      X=B
      DO I=1,NEQ
        CALL ZGES(AFACT,N,N,IPVT,X(:,I),0)
      ENDDO
!     =================================================================
!     ==  TEST                                                       ==
!     =================================================================
      IF(TTEST) THEN
        DEV=MAXVAL(ABS(MATMUL(A,X)-B))
        WRITE(*,*)'MAX ERROR OF LIB$MATRIXSOLVE ',DEV
      END IF
      RETURN
      END SUBROUTINE LIB_ESSL_ZGES
!
!     ..................................................................
      SUBROUTINE LIB_ESSL_ZGESVS(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  SOLVES THE LINEAR EQUATION SYSTEM AX=B                      **
!     **  WHERE A IS A(N,M)                                           **
!     **  IF A IS NOT SQUARE, THE EQUATION IS SOLVED IN A LEAST       **
!     **  SQUARE SENSE                                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(OUT):: X(M,NEQ)
      COMPLEX(8),INTENT(IN) :: B(N,NEQ)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      COMPLEX(8)            :: APLUSA(M,M)
      COMPLEX(8)            :: APLUSB(M,NEQ)
      REAL(8)               :: AMAT(2*M,2*M)
      REAL(8)               :: BMAT(2*M,NEQ)
      REAL(8)               :: XMAT(2*M,NEQ)
      REAL(8)               :: DEV
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
!     ******************************************************************
      APLUSA(:,:)=MATMUL(TRANSPOSE(A),A)
      AMAT(1:M,1:M)=REAL(APLUSA)
      AMAT(1:M,M+1:2*M)=REAL(CI*APLUSA)
      AMAT(M+1:2*M,1:M)=REAL(-CI*APLUSA)
      AMAT(M+1:2*M,M+1:2*M)=REAL(APLUSA)
      APLUSB(:,:)=MATMUL(TRANSPOSE(A),B)
      BMAT(1:M,:)=REAL(APLUSB)
      BMAT(M+1:2*M,:)=REAL(-CI*APLUSB)
      CALL LIB_ESSL_DGESVS(2*N,2*M,NEQ,AMAT,XMAT,BMAT)
!     CALL LIB_LAPACK_DGESV(2*M,2*M,NEQ,AMAT,XMAT,BMAT)
      X=CMPLX(XMAT(1:M,:),XMAT(M+1:2*M,:))
!     =================================================================
!     ==  TEST                                                       ==
!     =================================================================
      IF(TTEST) THEN
        DEV=MAXVAL(ABS(MATMUL(AMAT,XMAT)-BMAT))
        WRITE(*,*)'MAX ERROR OF LIB$MATRIXSOLVE ----',DEV
        DEV=MAXVAL(ABS(MATMUL(APLUSA,X)-APLUSB))
        WRITE(*,*)'MAX ERROR OF LIB$MATRIXSOLVE ---+',DEV
        DEV=MAXVAL(ABS(MATMUL(A,X)-B))
        WRITE(*,*)'MAX ERROR OF LIB$MATRIXSOLVE --++',DEV
      END IF
      RETURN
      END SUBROUTINE LIB_ESSL_ZGESVS
!
!     ..................................................................
      SUBROUTINE LIB_ESSL_DSPEV(N,H,E,U)
!     ******************************************************************
!     **                                                              **
!     **  DIAGONALIZES THE REAL, SQUARE MATRIX H AFTER SYMMETRIZATION **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                   **
!     **                                                              **
!     **         U(K,I)*H(K,L)*U(L,J)=DELTA(I,J)*E(I)                 **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL DSPEV  MATRIX DIAGONALIZATION P727                   **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      REAL(8)   ,INTENT(OUT):: U(N,N)
      REAL(8)               :: WORK1((N*(N+1))/2)
      REAL(8)               :: WORK2(2*N)
      INTEGER(4)            :: K,I,J
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: DEV
      REAL(8)  ,ALLOCATABLE :: EMAT(:,:)
!     ******************************************************************
!
!     ==================================================================
!     ====  STORE IN LOWER PACKED STORAGE MODE FOR DIAGONALIZATION    ==
!     ==================================================================
      K=0
      DO J=1,N
        DO I=J,N
          K=K+1
          WORK1(K)=0.5D0*(H(I,J)+H(J,I))
        ENDDO
      ENDDO
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
      CALL DSPEV(1,WORK1,E,U,N,N,WORK2,2*N) !->ESSL
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
!
!     ==================================================================
!     == TEST                                                         ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(EMAT(N,N))
!       == TEST EIGENVALUE EQUATION ====================================
        EMAT(:,:)=0.D0
        DO I=1,N
          EMAT(I,I)=E(I)
        ENDDO
        DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(U,EMAT)))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_ESSL_DSPEV')
        END IF
!       == TEST ORTHONORMALITY OF EIGENVECTORS =========================
        EMAT=MATMUL(TRANSPOSE(U),U)
        DO I=1,N
          EMAT(I,I)=EMAT(I,I)-1.D0
        ENDDO
        DEV=MAXVAL(ABS(EMAT))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('ORTHONORMALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_ESSL_DSPEV')
        END IF
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB_ESSL_ZHPEV(N,H,E,U)
!     ******************************************************************
!     **                                                              **
!     **  DIAGONALIZES THE HERMITEAN, SQUARE MATRIX H                 **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                   **
!     **                                                              **
!     **      CONJG(U(K,I))*H(K,L)*U(L,J)=DELTA(I,J)*E(I)             **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: ZHPEV :  MATRIX DIAGONALIZATION  P727               **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTERFACE
        SUBROUTINE EINFO(ICODE,INF1,INF2) !ESSL ERROR HANDLING ROUTINE
        INTEGER                       :: ICODE
        INTEGER ,INTENT(OUT),OPTIONAL :: INF1
        INTEGER ,INTENT(OUT),OPTIONAL :: INF2
        END SUBROUTINE EINFO
      END INTERFACE
      LOGICAL(4) ,PARAMETER :: TESSLERR=.FALSE.
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      COMPLEX(8),INTENT(OUT):: U(N,N)
      COMPLEX(8)            :: WORK1((N*(N+1))/2)
      REAL(8)               :: RWORK(4*N)
      INTEGER(4)            :: K,I,J
      CHARACTER(8)          :: SAV2101
      INTEGER(4)            :: I1,I2
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      COMPLEX(8)            :: CSVAR
      LOGICAL(4)            :: TCHK
      COMPLEX(8),ALLOCATABLE:: WORK3(:,:)
      INTEGER(4)            :: INFO
!CLEMENS
      INTEGER               :: INF1
      INTEGER               :: INF2
!     ******************************************************************
!
!     ==================================================================
!     ====  STORE IN LOWER PACKED STORAGE MODE FOR DIAGONALIZATION    ==
!     ==================================================================
      K=0
      DO J=1,N
        DO I=J,N
          K=K+1
          WORK1(K)=0.5D0*(H(I,J)+CONJG(H(J,I)))
        ENDDO
      ENDDO
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
      IF(TESSLERR) THEN
        CALL EINFO(0,INF1,INF2)
        CALL ERRSAV(2101,SAV2101)
        CALL ERRSET(2101,255,0,0,0,2101)
      END IF
      CALL ZHPEV(1,WORK1,E,U,N,N,RWORK,4*N)  !ESSL
!
!     ==================================================================
!     ====  OPTIONAL TEST                                             ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(WORK3(N,N)) 
        TCHK=.TRUE.
        DO I=1,N
          DO J=1,N
            CSVAR=(0.D0,0.D0)
            DO K=1,N
              CSVAR=CSVAR+U(I,K)*E(K)*CONJG(U(J,K))
            ENDDO
            IF(ABS(CSVAR-H(I,J)).GT.1.D-8) THEN
              WRITE(*,FMT='(2I5,5F10.5)')I,J,CSVAR,H(I,J),ABS(CSVAR-H(I,J))
              TCHK=.FALSE.
            END IF
          ENDDO
        ENDDO
        DEALLOCATE(WORK3)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$STOP('LIB_ESSL_ZHPEV')
        END IF
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB_ESSL_DSYGV(N,H,S,E,U)
!     ******************************************************************
!     **                                                              **
!     **  GENERALIZED REAL SYMMETRIC EIGENSYSTEM,                     ** 
!     **        HU-SUE=0,                                             **
!     **  WHERE A IS REAL SYMMETRIC                                   **
!     **  AND B IS REAL SYMMETRIC POSITIVE DEFINITE                   **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: DSYGV                                               **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(IN) :: S(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      REAL(8)   ,INTENT(OUT):: U(N,N)
      REAL(8)               :: H1(N,N)
      REAL(8)               :: S1(N,N)
      REAL(8)               :: AUX(2*N)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      INTEGER(4)            :: I
      REAL(8)   ,ALLOCATABLE:: EMAT(:,:)
      REAL(8)               :: DEV
!     ******************************************************************
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
      H1=H
      S1=S
      CALL DSYGV(1,H1,N,S1,N,E,U,N,N,AUX,2*N)  !ESSL
!
!     ==================================================================
!     ====  OPTIONAL TEST                                             ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(EMAT(N,N)) 
        EMAT(:,:)=0.D0
        DO I=1,N
          EMAT(I,I)=E(I)
        ENDDO
        DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(S,MATMUL(U,EMAT))))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_ESSL_DSYGV')
        END IF
        DEALLOCATE(EMAT)
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB_ESSL_DGEGV(N,H,S,E,U)
!     ******************************************************************
!     **                                                              **
!     **  GENERALIZED REAL NON-SYMMETRIC EIGENSYSTEM,                 ** 
!     **        HU-SUE=0,                                             **
!     **  WHERE H AND S ARE REAL, GENERAL MATRICES                    **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: DGEGV                                               **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **   2) NOTE THAT THE EIGENVECTORS ARE NOT ORTHONORMAL!!!!      **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(IN) :: S(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      REAL(8)   ,INTENT(OUT):: U(N,N)
      REAL(8)               :: H1(N,N)
      REAL(8)               :: S1(N,N)
      REAL(8)               :: ALPHA(N)
      REAL(8)               :: BETA(N)
      REAL(8)               :: AUX(3*N)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      INTEGER(4)            :: I
      REAL(8)   ,ALLOCATABLE:: EMAT(:,:)
      REAL(8)               :: DEV
!     ******************************************************************
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
      H1=H
      S1=S
      CALL DGEGV(1,H1,N,S,N,ALPHA,BETA,U,N,N,AUX,3*N)  !ESSL
      DO I=1,N
        IF(BETA(I).EQ.0.D0) THEN
          CALL ERROR$MSG('EIGENVALUE IS INFINITE')
          CALL ERROR$I4VAL('I',I)
          CALL ERROR$STOP('LIB_ESSL_DGEGV')
        END IF
      ENDDO
      E=ALPHA/BETA
!
!     ==================================================================
!     ====  OPTIONAL TEST                                             ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(EMAT(N,N)) 
        EMAT(:,:)=0.D0
        DO I=1,N
          EMAT(I,I)=E(I)
        ENDDO
        DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(S,MATMUL(U,EMAT))))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_ESSL_DGEGV')
        END IF
        DEALLOCATE(EMAT)
      END IF
      RETURN
      END
#ELSE     
!***********************************************************************
!***********************************************************************
!****                                                               ****
!****  DRIVER ROUTINES FOR LAPACK                                   ****
!****                                                               ****
!***********************************************************************
!***********************************************************************
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_DGETRI(N,A,AINV)
!     ******************************************************************
!     **                                                              **
!     **  INVERTS THE REAL, SQUARE MATRIX A                           **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: DGEICD                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: A(N,N)
      REAL(8)   ,INTENT(OUT):: AINV(N,N)
      INTEGER(4)            :: NAUX
      REAL(8)               :: AUX(100*N)
      INTEGER(4)            :: IPIV(N)
      INTEGER(4)            :: INFO
!     ******************************************************************
      NAUX=100*N
      AINV(1:N,1:N)=A(1:N,1:N)
!
!     ==================================================================
!     == PERFORM LU FACTORIZATION OF A                                ==
!     ==================================================================
      CALL DGETRF(N,N,AINV,N,IPIV,INFO) !LAPACK
!
!     ==================================================================
!     == INVERT A USING THE LU FACTORIZATION                          ==
!     ==================================================================
      CALL DGETRI(N,AINV,N,IPIV,AUX,NAUX,INFO) !LAPACK
!
!     ==================================================================
!     == CHECK ERROR CODE                                             ==
!     ==================================================================
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB_LAPACK_DGETRI')
        ELSE
          CALL ERROR$MSG('U(I,I) IS EXACTLY ZERO')
          CALL ERROR$MSG('MATRIX IS SINGULAR AND ITS INVERSE COULD NOT BE COMPUTED.')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_LAPACK_DGETRI')
        END IF
      END IF
      RETURN
      END SUBROUTINE LIB_LAPACK_DGETRI
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_zGETRI(N,A,AINV)
!     ******************************************************************
!     **                                                              **
!     **  INVERTS THE REAL, SQUARE MATRIX A                           **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: DGEICD                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      complex(8),INTENT(IN) :: A(N,N)
      complex(8),INTENT(OUT):: AINV(N,N)
      INTEGER(4)            :: NAUX
      complex(8)            :: AUX(100*N)
      INTEGER               :: IPIV(N)
      INTEGER               :: INFO
!     ******************************************************************
      NAUX=100*N
      AINV(1:N,1:N)=A(1:N,1:N)
!
!     ==================================================================
!     == PERFORM LU FACTORIZATION OF A                                ==
!     ==================================================================
      CALL zGETRF(N,N,AINV,N,IPIV,INFO) !LAPACK
!
!     ==================================================================
!     == INVERT A USING THE LU FACTORIZATION                          ==
!     ==================================================================
      CALL zGETRI(N,AINV,N,IPIV,AUX,NAUX,INFO) !LAPACK
!
!     ==================================================================
!     == CHECK ERROR CODE                                             ==
!     ==================================================================
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB_LAPACK_DGETRI')
        ELSE
          CALL ERROR$MSG('U(I,I) IS EXACTLY ZERO')
          CALL ERROR$MSG('MATRIX IS SINGULAR AND ITS INVERSE COULD NOT BE COMPUTED.')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_LAPACK_DGETRI')
        END IF
      END IF
      RETURN
      END SUBROUTINE LIB_LAPACK_ZGETRI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_LAPACK_DGESVD(N,M,A,U,S,Vt)
!     **************************************************************************
!     ** SINGULAR VALUE DECOMPOSITION OF THE NON-SQUARE MATRIX A              **
!     ** A=U*SIGMA*Vt  
!     **  SIGMA IS AN M-TIMES-N DIAGONAL MATRIX WITH DIAGONAL ELEMENTS S(I)   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: M
      INTEGER(4),INTENT(IN)  :: N
      REAL(8)   ,INTENT(IN)  :: A(M,N)
      REAL(8)   ,INTENT(OUT) :: U(M,M)
      REAL(8)   ,INTENT(OUT) :: Vt(N,N)
      REAL(8)   ,INTENT(OUT) :: S(M)     ! SINGULAR VALUES IN DESCENDING ORDER
      REAL(8)   ,ALLOCATABLE :: WORK(:)
      REAL(8)                :: ACOPY(M,N)
      INTEGER(4)             :: LWORK
      INTEGER(4)             :: INFO
      INTEGER(4)             :: I
!     **************************************************************************
      LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
      ALLOCATE(WORK(LWORK))
      ACOPY=A  ! ACOPY WILL BE OVERWRITTEN
      S=0.D0
      U=0.D0
      Vt=0.D0
      WORK=0.D0
      CALL DGESVD('A','A',M,N,ACOPY,M,S,U,M,Vt,N,WORK,LWORK,INFO)
      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('THE I-TH ARGUMENT JHAS AN ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$I4VAL('1ST INDEX OF MATRIX',M)
        CALL ERROR$I4VAL('2ND INDEX OF MATRIX',N)
        CALL ERROR$STOP('LIB_LAPACK_DGESVD')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$MSG('DBDSQR DID NOT CONVERGE')
        CALL ERROR$I4VAL('NUMBER OF UNCONVERGED SUPERDIAGONALS',INFO)
        CALL ERROR$I4VAL('1ST INDEX OF MATRIX',M)
        CALL ERROR$I4VAL('2ND INDEX OF MATRIX',N)
        CALL ERROR$STOP('LIB_LAPACK_DGESVD')
      END IF
      DEALLOCATE(WORK)
      RETURN
      END SUBROUTINE LIB_LAPACK_DGESVD
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_DGESV(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  DRIVER ROUTINE FOR LAPACK ROUTINE DGESV                     **
!     **                                                              **
!     **  COMPUTES THE SOLUTION TO A REAL SYSTEM OF LINEAR EQUATIONS  **
!     **        A * X = B,                                            **
!     **  WHERE A IS A N-BY-N MATRIX AND X AND B ARE N-BY-NEQ MATRICES**
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(OUT):: X(M,NEQ)
      REAL(8)   ,INTENT(IN) :: B(N,NEQ)
      REAL(8)   ,allocatable:: A1(:,:)
      INTEGER               :: INFO
      INTEGER(4)            :: IPIV(N)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: SVAR1,SVAR2
!     ******************************************************************
      IF(N.NE.M) THEN
        CALL ERROR$MSG('WORKS ONLY FOR SQUARE MATRICES')
        CALL ERROR$STOP('LIB_LAPACK_DGESV')
      END IF
!
!     ==================================================================
!     == SOLVE SPECIAL CASE WITH DIMENSIONS 1 FIRST                   ==
!     ==================================================================
      IF(N.EQ.1.AND.M.EQ.1) THEN
        X(1,:)=B(1,:)/A(1,1)
        RETURN
      END IF
!
!     ==================================================================
!     == NOW CALL LAPACK ROUTINE                                      ==
!     ==================================================================
      allocate(a1(n,m))
      A1=A
      X=B
      CALL DGESV(N,NEQ,A1,N,IPIV,X,N,INFO )
      deallocate(a1)
      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('ERROR EXIT FROM DGESV')
        CALL ERROR$MSG('THE I-TH ARGUMENT HAD AN ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$STOP('LIB$MATRIXSOLVENEW')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$MSG('ERROR EXIT FROM DGESV')
        CALL ERROR$MSG('U(I,I) IS EXACTLY ZERO.')
        CALL ERROR$MSG('THE FACTORIZATION HAS BEEN COMPLETED,')
        CALL ERROR$MSG('BUT THE FACTOR U IS EXACTLY SINGULAR,')
        CALL ERROR$MSG('SO THE SOLUTION COULD NOT BE COMPUTED.')
        CALL ERROR$I4VAL('I',INFO)
        CALL ERROR$STOP('LIB$MATRIXSOLVENEW')
      END IF
!
!     ==================================================================
!     ==  TEST                                                        ==
!     ==================================================================
      IF(TTEST) THEN
        SVAR1=MAXVAL(ABS(MATMUL(A,X)-B))
        SVAR2=MAXVAL(ABS(A))/MAXVAL(ABS(B))
        IF(SVAR1/SVAR2.GT.1.D-4) THEN
          CALL ERROR$R8VAL('SVAR1',SVAR1)
          CALL ERROR$R8VAL('SVAR2',SVAR2)
          CALL ERROR$R8VAL('SVAR1/SVAR2',SVAR1/SVAR2)
          CALL ERROR$STOP('LIB$MATRIXSOLVENEW')
        END IF
      END IF

      RETURN
      END SUBROUTINE LIB_LAPACK_DGESV
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_LAPACK_DGELSD(N,M,NEQ,A,X,B)
!     **************************************************************************
!     **  DRIVER ROUTINE FOR LAPACK ROUTINE DGELSD                            **
!     **                                                                      **
!     **  COMPUTES THE MINIMUM-NORM SOLUTION TO A COMPLEX LINEAR LEAST        **
!     **  SQUARES PROBLEM:                                                    **
!     **              MINIMIZE 2-NORM(| B - A*X |)                            **
!     **  USING THE SINGULAR VALUE DECOMPOSITION (SVD) OF A.                  **
!     **  A IS AN M-BY-N MATRIX WHICH MAY BE RANK-DEFICIENT.                  **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(OUT):: X(M,NEQ)
      REAL(8)   ,INTENT(IN) :: B(N,NEQ)
      double precision            :: A1(N,M)
      double precision,ALLOCATABLE:: B1(:,:)
      double precision            :: S(N)        ! SINGULAR VALUES
      double precision,PARAMETER  :: RCOND=-1.D0 ! USE MACHINE PRECISION
      INTEGER                     :: RANK
      INTEGER                     :: LWORK
      double precision,ALLOCATABLE:: WORK(:)
      INTEGER         ,ALLOCATABLE:: IWORK(:)
      INTEGER                     :: LIWORK
      INTEGER                     :: INFO
      INTEGER                     :: nrhs
      INTEGER                     :: LDA,LDB
      INTEGER                     :: N1,M1,SMLSIZ,MINMN,MAXMN,NLVL
      INTEGER                     :: ILAENV
      EXTERNAL                    :: ILAENV
      LOGICAL         ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)         ,PARAMETER  :: TOL=1.D-5
      REAL(8)                     :: SVAR1,SVAR2
      character(6)    ,parameter  :: type='dgels'
!     **************************************************************************
      IF(N.LT.1.OR.M.LT.1) THEN
        CALL ERROR$MSG('DIMENSIONS MUST BE NONZERO AND POSITIVE')
        CALL ERROR$I4VAL('N',N)
        CALL ERROR$I4VAL('M',M)
        CALL ERROR$STOP('LIB_LAPACK_DGELSD')
      END IF
      M1=N    ! LAPACK USES M AND N OPPOSITE 
      N1=M
      nrhs=neq
      MINMN=MIN(M1,N1)
      MAXMN=MAX(M1,N1)
!
!     ==========================================================================
!     == interface to dgelsD                                                  ==
!     ==========================================================================
      if(type.eq.'dgelsd') then
        SMLSIZ=ILAENV(9,'DGELSD',' ',0,0,0,0 )
        NLVL = MAX(0,INT(LOG(REAL(MINMN,kind=8)/real(SMLSIZ+1,KIND=8))/LOG(2.D0))+1)
        LIWORK=MAX(1,3*MINMN*NLVL+11*MINMN)
        LWORK=12*MINMN+2*MINMN*SMLSIZ+8*MINMN*NLVL+MINMN*NEQ+(SMLSIZ+1)**2
        ALLOCATE(WORK(LWORK))
        ALLOCATE(IWORK(LIWORK))
        LDB=MAX(1,MAXMN)
        ALLOCATE(B1(LDB,NEQ))
        B1=0.D0
        B1(1:M1,:)=B(:,:)
        A1=A
! the call to dgelsd failed using mkl
        CALL DGELSD(M1,N1,Nrhs,A1,M1,B1,LDB,S,RCOND,RANK,WORK,LWORK,IWORK,INFO )
        deallocate(iwork)
        deallocate(work)
        X=B1(1:N1,:)
        deallocate(B1)
        IF(INFO.NE.0) THEN
          IF(INFO.LT.0) THEN
            CALL ERROR$MSG('I-TH ARGUMENT HAD AN ILLEGAL VALUE')
            CALL ERROR$I4VAL('I',-INFO)
            CALL ERROR$STOP('LIB_LAPACK_DGELSD')
          ELSE
            CALL ERROR$MSG('ALGORITHM FOR COMPUTING SVD FAILED TO CONVERGE')
            CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
            CALL ERROR$MSG('BIDIAGONAL FORM DID NOT CONVERGE TO ZERO.')
            CALL ERROR$I4VAL('I',INFO)
            CALL ERROR$STOP('LIB_LAPACK_DGELSD')
          END IF
        END IF
!
!     ==========================================================================
!     == interface to dgels                                                   ==
!     ==========================================================================
      else if(type.eq.'dgels') then
        LDA=MAX(1,M1)
        LDB=MAX(1,MAXMN)
        ALLOCATE(B1(LDB,NEQ))
        B1=0.D0
        B1(1:M1,:)=B(:,:)
        A1=A
        LWORK=MINMN+max(1,max(NRHS,maxMN))
        ALLOCATE(WORK(LWORK))
        CALL DGELS('N',M1,N1,NRHS,A1,LDA,B1,LDB,WORK,LWORK,INFO)
        DEALLOCATE(WORK)
        X=B1(1:N1,:)
        DEALLOCATE(B1)
!
!     ==========================================================================
!     == interface to dgelsS                                                  ==
!     ==========================================================================
      else if(type.eq.'dgelss') then
        LDA=MAX(1,M1)
        LDB=MAX(1,MAX(M1,N1))
        ALLOCATE(B1(LDB,NEQ))
        B1=0.D0
        B1(1:M1,:)=B(:,:)
        A1=A
        LWORK=3*MINMN+max(2*MINMN,max(NRHS,maxMN))
        ALLOCATE(WORK(LWORK))
        CALL DGELSS(M1,N1,NRHS,A1,LDA,B1,LDB,S,RCOND,RANK,WORK,LWORK,INFO)
        DEALLOCATE(WORK)
        X=B1(1:N1,:)
        DEALLOCATE(B1)
      else
        call error$msg('type not recognized') 
        call error$stop('lib_lapack_dgelsd')
      end if
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        SVAR1=MAXVAL(ABS(MATMUL(A,X)-B))
        SVAR2=TOL*MAXVAL(ABS(B))/MAXVAL(ABS(B))
        IF(SVAR1.GT.SVAR2) THEN
          CALL ERROR$MSG('TEST FAILED')
          CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
          CALL ERROR$MSG('BIDIAGONAL FORM DID NOT CONVERGE TO ZERO.')
          CALL ERROR$R8VAL('MAX DEV',SVAR1)
          CALL ERROR$R8VAL('ALLOWED DEV',SVAR2)
          CALL ERROR$STOP('LIB_LAPACK_DGELSD')
        END IF
      END IF
      RETURN
      END SUBROUTINE LIB_LAPACK_DGELSD
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_ZGESV(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  OMPLEX SYSTEM  OF  LINEAR EQUATIONS                         **
!     **              A * X = B,                                      **
!     **  WHERE A IS A(N,M) IS A SQUARE MATRIX                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(OUT):: X(M,NEQ)
      COMPLEX(8),INTENT(IN) :: B(N,NEQ)
      COMPLEX(8)            :: A1(N,M)
      COMPLEX(8)            :: B1(N,NEQ)
      INTEGER               :: INFO
      REAL(8)               :: RCOND=-1.D0
      INTEGER               :: LDWORK
      INTEGER               :: IRANK
      INTEGER               :: IPIVOT(N)
      REAL(8)               :: SING(N)
      REAL(8)   ,ALLOCATABLE:: WORK(:)
      INTEGER               :: N1,M1,NEQ1
!     ******************************************************************
      IF(N.NE.M) THEN
        CALL ERROR$MSG('ONLY SYMMETRIC MATRICES ALLOWED')
        CALL ERROR$STOP('LIB_LAPACK_ZGESV')
      END IF

      LDWORK=3*MIN(M,N)+MAX(2*MIN(M,N),MAX(M,N),NEQ)
!     -- USE 3*M+3*N+NEQ
      N1=N
      M1=M
      NEQ1=NEQ
!
!     ===================================================================
!     == SOLVE EQUATION SYSTEM                                         ==
!     ===================================================================
      A1=A 
      X=B
      ALLOCATE(WORK(LDWORK))
      CALL ZGESV(N,NEQ,A1,N,IPIVOT,X,N,INFO)
      DEALLOCATE(WORK)
!
!     ===================================================================
!     == CHECK ERROR CODE                                              ==
!     ===================================================================
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('ITH ARGUMENT HAS ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LLIB_LAPACK_ZGESV')
        ELSE IF(INFO.GT.0) THEN
          CALL ERROR$MSG('PROBLEM IS SINGULAR. NO SOLUTION CAN BE COMPUTED')
          CALL ERROR$MSG('U(I,I) IS  EXACTLY  ZERO. THE FACTORIZATION HAS BEEN') 
          CALL ERROR$MSG('COMPLETED, BUT THE FACTOR U IS EXACTLY SINGULAR,')
          CALL ERROR$MSG('SO THE SOLUTION COULD NOT  BE COMPUTED.')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_LAPACK_ZGESV')
        END IF
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_ZGELSD(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  DRIVER ROUTINE FOR LAPACK ROUTINE ZGELSD                    **
!     **                                                              **
!     **  COMPUTES THE MINIMUM-NORM SOLUTION TO A COMPLEX LINEAR LEAST**
!     **  SQUARES PROBLEM:                                            **
!     **              MINIMIZE 2-NORM(| B - A*X |)                    **
!     **  USING THE SINGULAR VALUE DECOMPOSITION (SVD) OF A.          **
!     **  A IS AN M-BY-N MATRIX WHICH MAY BE RANK-DEFICIENT.          **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(OUT):: X(M,NEQ)
      COMPLEX(8),INTENT(IN) :: B(N,NEQ)
      COMPLEX(8)            :: A1(N,M)
      COMPLEX(8),ALLOCATABLE:: B1(:,:)
      REAL(8)               :: S(N)        ! SINGULAR VALUES
      REAL(8)   ,PARAMETER  :: RCOND=-1.D0 ! USE MACHINE PRECISION
      INTEGER(4)            :: RANK
      INTEGER(4)            :: LCWORK      ! SIZE OF CWORK
      INTEGER(4)            :: LRWORK      ! SIZE OF RWORK
      INTEGER(4)            :: LIWORK      ! SIZE OF IWORK
      COMPLEX(8),ALLOCATABLE:: CWORK(:)
      REAL(8)   ,ALLOCATABLE:: RWORK(:)
      INTEGER(4),ALLOCATABLE:: IWORK(:)
      INTEGER(4)            :: INFO
      INTEGER               :: LDB
      INTEGER               :: N1,M1,SMLSIZ,MINMN,NLVL
      INTEGER               :: ILAENV
      EXTERNAL              :: ILAENV
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)   ,PARAMETER  :: TOL=1.D-5
      REAL(8)               :: SVAR1,SVAR2
!     ******************************************************************
      M1=N    ! LAPACK USES M AND N OPPOSITE 
      N1=M
      SMLSIZ=ILAENV(9,'ZGELSD',' ',0,0,0,0 )
      MINMN=MIN(M1,N1)
      NLVL = MAX(0,INT(LOG(REAL(MINMN/(SMLSIZ+1),KIND=8))/LOG(2.D0))+1)
      LCWORK=2*MINMN+MINMN*NEQ
      LRWORK=10*MINMN+2*MINMN*SMLSIZ+8*MINMN*NLVL+3*SMLSIZ*NEQ+(SMLSIZ+1)**2
      LIWORK=MAX(1,3*MINMN*NLVL+11*MINMN)
      ALLOCATE(CWORK(LCWORK))
      ALLOCATE(RWORK(LRWORK))
      ALLOCATE(IWORK(LIWORK))
      LDB=MAX(1,MAX(M1,N1))
      ALLOCATE(B1(LDB,NEQ))
      B1=0.D0
      B1(1:M1,:)=B(:,:)
      A1=A
      CALL ZGELSD(M1,N1,NEQ,A1,M1,B1,LDB,S,RCOND,RANK,CWORK,LCWORK,RWORK,IWORK,INFO )
      X=B1(1:N1,:)
      DEALLOCATE(CWORK)
      DEALLOCATE(RWORK)
      DEALLOCATE(IWORK)
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB_LAPACK_ZGELSD')
        ELSE
          CALL ERROR$MSG('ALGORITHM FOR COMPUTING SVD FAILED TO CONVERGE')
          CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
          CALL ERROR$MSG('BIDIAGONAL FORM DID NOT CONVERGE TO ZERO.')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_LAPACK_ZGELSD')
        END IF
      END IF
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      IF(TTEST) THEN
        SVAR1=MAXVAL(ABS(MATMUL(A,X)-B))
        SVAR2=TOL*MAXVAL(ABS(B))/MAXVAL(ABS(B))
        IF(SVAR1.GT.SVAR2) THEN
          CALL ERROR$MSG('TEST FAILED')
          CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
          CALL ERROR$MSG('BIDIAGONAL FORM DID NOT CONVERGE TO ZERO.')
          CALL ERROR$R8VAL('MAX DEV',SVAR1)
          CALL ERROR$R8VAL('ALLOWED DEV',SVAR2)
          CALL ERROR$STOP('LIB_LAPACK_ZGELSD')
        END IF
      END IF
      RETURN
      END SUBROUTINE LIB_LAPACK_ZGELSD
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_DSYEV(N,H,E,U)
!     ******************************************************************
!     **                                                              **
!     **  DIAGONALIZES THE REAL, SQUARE MATRIX H AFTER SYMMETRIZATION **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                   **
!     **                                                              **
!     **         H(I,K)*U(K,J)=U(I,J)*E(J)                            **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4) ,PARAMETER :: TESSLERR=.FALSE.
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      REAL(8)   ,INTENT(OUT):: U(N,N)
      REAL(8)               :: WORK(3*N)
      INTEGER(4)            :: INFO
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: DEV
      REAL(8)  ,ALLOCATABLE :: EMAT(:,:)
      INTEGER(4)            :: I
!     ******************************************************************
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
      U=0.5D0*(H+TRANSPOSE(H))
      CALL DSYEV('V','U',N,U,N,E,WORK,3*N,INFO )
!
!     ==================================================================
!     == CHECK ERROR CODE                                             ==
!     ==================================================================
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB_LAPACK_DSYEV')
        ELSE
          CALL ERROR$MSG('THE ALGORITHM FAILED TO CONVERGE')
          CALL ERROR$MSG('OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
          CALL ERROR$MSG('TRIDIAGONAL FORM DID NOT CONVERGE TO ZERO')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_LAPACK_DSYEV')
        END IF
      END IF
!
!     ==================================================================
!     == TEST                                                         ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(EMAT(N,N))
!       == TEST EIGENVALUE EQUATION ====================================
        EMAT(:,:)=0.D0
        DO I=1,N
          EMAT(I,I)=E(I)
        ENDDO
        DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(U,EMAT)))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_DSYEV')
        END IF
!       == TEST ORTHONORMALITY OF EIGENVECTORS =========================
        EMAT=MATMUL(TRANSPOSE(U),U)
        DO I=1,N
          EMAT(I,I)=EMAT(I,I)-1.D0
        ENDDO
        DEV=MAXVAL(ABS(EMAT))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('ORTHONORMALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_DSYEV')
        END IF
      END IF
      RETURN
      END SUBROUTINE LIB_LAPACK_DSYEV
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_ZHEEV(N,H,E,U)
!     ******************************************************************
!     **                                                              **
!     **  DIAGONALIZES THE HERMITEAN, SQUARE MATRIX H                 **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                   **
!     **                                                              **
!     **      CONJG(U(K,I))*H(K,L)*U(L,J)=DELTA(I,J)*E(I)             **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: ZHPEV :  MATRIX DIAGONALIZATION  P727               **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      COMPLEX(8),INTENT(OUT):: U(N,N)
      integer(4),parameter  :: lwmax=2000
      COMPLEX(8)            :: CWORK(lwmax)
      REAL(8)               :: RWORK(3*N-2)
      INTEGER(4)            :: I
      LOGICAL   ,PARAMETER  :: TTEST=.false.
      COMPLEX(8),ALLOCATABLE:: RES(:,:)
      INTEGER(4)            :: lwork
      INTEGER(4)            :: INFO
      COMPLEX(8)            :: ch(n,n)
      external zheev
!     ******************************************************************
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
      IF(LWMAX.LT.2*N) THEN
        CALL ERROR$MSG('HARDWIRED LIMIT LWMAX SMALLER THAN 2*N')
        CALL ERROR$I4VAL('LWMAX',LWMAX)
        CALL ERROR$I4VAL('2*N',2*N)
        CALL ERROR$STOP('LIB_LAPACK_ZHEEV')
      END IF
      u=0.5d0*(h+transpose(conjg(h)))
      LWORK=-1
      CALL ZHEEV('V','L',N,u,N,E,CWORK,LWORK,RWORK,INFO) !LAPACK
      LWORK=MIN(LWMAX,INT(CWORK(1)))
      CALL ZHEEV('V','L',N,u,N,E,CWORK,LWORK,RWORK,INFO) !LAPACK
!
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT TO ZHEEV HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB_LAPACK_ZHEEV')
        ELSE
          CALL ERROR$MSG('THE ALGORITHM FAILED TO CONVERGE')
          CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
          CALL ERROR$MSG('TRIDIAGONAL FORM DID NOT CONVERGE TO ZERO.')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_LAPACK_ZHEEV')
        END IF
      END IF
!
!     ==================================================================
!     ====  OPTIONAL TEST                                             ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(RES(N,N))
        RES=0.5D0*(H+TRANSPOSE(CONJG(H)))
        RES=MATMUL(RES,U)
        DO I=1,N
          RES(:,I)=RES(:,I)-U(:,I)*E(I)
        ENDDO
        IF(MAXVAL(ABS(RES)).GT.1.D-10) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV ',MAXVAL(ABS(RES)))
          CALL ERROR$STOP('LIB_LAPACK_ZHEEV')
        END IF
        RES=MATMUL(CONJG(TRANSPOSE(U)),U)
        DO I=1,N
          RES(I,I)=RES(I,I)-CMPLX(1.D0,0.D0)
        ENDDO
        IF(MAXVAL(ABS(RES)).GT.1.D-10) THEN
          CALL ERROR$MSG('ORTHONORMALITY  TEST FAILED')
          CALL ERROR$R8VAL('DEV ',MAXVAL(ABS(RES)))
          CALL ERROR$STOP('LIB_LAPACK_ZHEEV')
        END IF
        DEALLOCATE(RES)
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_DSYGV(N,H,S,E,U)
!     ******************************************************************
!     **                                                              **
!     **  DIAGONALIZES THE REAL, SQUARE MATRIX H AFTER SYMMETRIZATION **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                   **
!     **                                                              **
!     **         U(K,I)*H(K,L)*U(L,J)=DELTA(I,J)*E(I)                 **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL DSPEV  MATRIX DIAGONALIZATION P727                   **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(IN) :: S(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      REAL(8)   ,INTENT(OUT):: U(N,N)
      REAL(8)               :: B(N,N)      
      REAL(8)               :: WORK(3*N)
      INTEGER(4)            :: INFO
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)   ,ALLOCATABLE:: EMAT(:,:)
      REAL(8)               :: DEV
      INTEGER(4)            :: I
!     ******************************************************************
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
      U=H
      B=S
      CALL DSYGV(1,'V','U',N,U,N,B,N,E,WORK,3*N,INFO )
!
!     ==================================================================
!     == CHECK ERROR CODE                                             ==
!     ==================================================================
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP(' LIB_LAPACK_DSYGV')
        ELSE
          IF(INFO.LE.N) THEN
            CALL ERROR$MSG('DSYEV FAILED TO CONVERGE')
            CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
            CALL ERROR$MSG('ITRIDIAGONAL FORM DID NOT CONVERGE TO ZERO')
            CALL ERROR$I4VAL('I',INFO)
            CALL ERROR$STOP(' LIB_LAPACK_DSYGV')
          ELSE
            CALL ERROR$MSG('LEADING MINOR OF ORDER I OF B IS NOT POSITIVE DEFINITE')
            CALL ERROR$MSG('THE FACTORIZATION OF B COULD NOT BE COMPLETED')
            CALL ERROR$MSG('NO EIGENVALUES OR EIGENVECTORS WERE COMPUTED')
            CALL ERROR$I4VAL('I',INFO-N)
            CALL ERROR$STOP(' LIB_LAPACK_DSYGV')
          END IF
        END IF
      END IF
!
!     ==================================================================
!     == TEST                                                         ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(EMAT(N,N))
        EMAT=0.D0
        DO I=1,N
          EMAT(I,I)=E(I)
        ENDDO
        DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(S,MATMUL(U,EMAT))))
        IF(DEV.GT.1.D-5) THEN
          CALL ERROR$MSG('EIGENVALUE EQUATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP(' LIB_LAPACK_DSYGV')
        END IF
        DEALLOCATE(EMAT)
      END IF
      RETURN
      END SUBROUTINE LIB_LAPACK_DSYGV
!
!     ......................................................................
      SUBROUTINE LIB_LAPACK_ZHEGV(N,H,S,E,VEC)
!     **                                                                 **
!     ** SOLVES THE GENERALIZED, REAL NON-SYMMETRIC EIGENVALUE PROBLEM   **
!     **      [H(:,:)-E(I)*S(:,:)]*VEC(:,I)=0                            **
!     **                                                                 **
!     ** REMARK: H AND S MUST BE HERMITEANC                              **
!     **         S MUST BE POSITIVE DEFINITE                             **
!     **         EIGENVECTORS ARE ORTHONORMAL IN THE SENSE               **
!     **             MATMUL(TRANSPOSE(VEC),MATMUL(S,VEC))=IDENTITY       **
!     **                                                                 **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)    ! HAMITON MATRIX
      COMPLEX(8),INTENT(IN) :: S(N,N)    ! OVERLAP MATRIX
      REAL(8)   ,INTENT(OUT):: E(N)      ! EIGENVALUES
      COMPLEX(8),INTENT(OUT):: VEC(N,N)  ! EIGENVECTORS
      INTEGER               :: LDWORK
      COMPLEX(8),ALLOCATABLE:: WORK(:)
      COMPLEX(8)            :: S1(N,N)
      REAL(8)               :: RWORK(3*N-2)
      INTEGER               :: INFO
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: DEV
      INTEGER               :: I
!     *********************************************************************
!
!     ========================================================================
!     == TEST IF INPUT MATRICES ARE SYMMETRIC                               ==
!     ========================================================================
      IF(TTEST) THEN
        DEV=SUM(ABS(H-TRANSPOSE(CONJG(H))))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('HAMILTON MATRIX NOT HERMITEAN')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_ZHEGV')
        END IF
        DEV=SUM(ABS(S-TRANSPOSE(CONJG(S))))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('OVERLAP MATRIX NOT HERMITEAN')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_ZHEGV')
        END IF
      END IF
!
!     ========================================================================
!     == CALL LAPACK ROUTINE                                                ==
!     ========================================================================
      VEC=0.5D0*(H+TRANSPOSE(CONJG(H)))  ! LAPACK ROUTINE OVERWRITES HAMILTONIAN WITH EIGENVECTORS
      S1=S
      !DO WORKSPAVE QUERY
      allocate(WORK(1))
      LDWORK=-1
      CALL ZHEGV(1,'V','U',N,VEC,N,S1,N,E,WORK,LDWORK,RWORK,INFO)
      LDWORK=WORK(1)
      deallocate(WORK)
      allocate(WORK(LDWORK)) 
      CALL ZHEGV(1,'V','U',N,VEC,N,S1,N,E,WORK,LDWORK,RWORK,INFO)
      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('ITH ARGUMENT OF ZHGEV HAS ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$STOP('LIB_LAPACK_ZHEGV')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$MSG('FAILED')
        CALL ERROR$I4VAL('INFO',INFO)
        CALL ERROR$STOP('LIB_LAPACK_ZHEGV')
      END IF
!
!     ========================================================================
!     == TEST RESULT OF THE ROUTINE                                         ==
!     ========================================================================
      IF(TTEST) THEN
        S1=MATMUL(TRANSPOSE(CONJG(VEC)),MATMUL(S,VEC))
        DEV=0.D0
        DO I=1,N
          S1(I,I)=S1(I,I)-(1.D0,0.D0)
          DEV=MAX(DEV,MAXVAL(ABS(MATMUL(H-E(I)*S,VEC(:,I)))))
        ENDDO
        IF(DEV.GT.1.D-6) THEN
          CALL ERROR$MSG('GENERAL EIGENVALUE TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_ZHEGV')
        END IF
        DEV=SUM(ABS(S1))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('EIGENSTATES NOT ORTHONORMAL')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_ZHEGV')
        END IF
      END IF
      RETURN
      END

#ENDIF
! 
!*******************************************************************************
!*******************************************************************************
!****                                                                      *****
!****              INTERFACES TO THE BLAS SUBROUTINES                      *****
!****                                                                      *****
!****                                                                      *****
!*******************************************************************************
!*******************************************************************************
!
!DGEMUL: MATRIX MULTIPLICATION P441
!     ..................................................................
      SUBROUTINE LIB$MATMULR8(N,M,L,A,B,C)
!     ******************************************************************
!     **  MATRIX MULTPLICATION   A*B=C                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(IN) :: B(M,L)
      REAL(8)   ,INTENT(OUT):: C(N,L)
      CHARACTER(8),PARAMETER:: LIB='ESSL'
      REAL(8)     ,PARAMETER:: ONE=1.D0
      REAL(8)     ,PARAMETER:: ZERO=0.D0
!     ******************************************************************
#IF DEFINED(CPPVAR_BLAS_ESSL)
      CALL DGEMUL(A,N,'N',B,M,'N',C,N,N,M,L)
#ELSE
      CALL DGEMM('N','N',N,L,M,ONE,A,N,B,M,ZERO,C,N)
#ENDIF
      RETURN
      END
!
!ZGEMUL: MATRIX MULTIPLICATION P441
!     . ................................................................
      SUBROUTINE LIB$MATMULC8(N,M,L,A,B,C)
!     ******************************************************************
!     **  MATRIX MULTPLICATION   A*B=C                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: L
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(IN) :: B(M,L)
      COMPLEX(8),INTENT(OUT):: C(N,L)
      CHARACTER(8),PARAMETER:: LIB='ESSL'
      COMPLEX(8)  ,PARAMETER:: ONE=(1.D0,0.D0)
      COMPLEX(8)  ,PARAMETER:: ZERO=(0.D0,0.D0)
!     ******************************************************************
#IF DEFINED(CPPVAR_BLAS_ESSL)
      CALL ZGEMUL(A,N,'N',B,M,'N',C,N,N,M,L)
!#ELIF DEFINED(CPPVAR_BLAS_ATLAS)
!      CALL ZGEMM('N','N',N,L,M,ONE,A,N,B,M,ZERO,C,N)
#ELSE
      C(:,:)=(0.D0,0.D0)
      CALL ZGEMM('N','N',N,L,M,ONE,A,N,B,M,ZERO,C,N)
#ENDIF
      RETURN
      END
!
!ZGEMM(TRANSA,TRANSB,L,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC):
!C=BETA*C+ALPHA*A*B P456
!     ..................................................................
      SUBROUTINE LIB$ADDPRODUCTC8(TID,N,M,L,A,B,C)
!     ******************************************************************
!     **  ADD THE PRODUCT OF TWO MATRICES                             **
!     **  C=C+MATMUL(A,B)                                             **
!     **  FOR TID=.TRUE., A AND C MAY USE IDENTICAL STORAGE           **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)   :: TID
      INTEGER(4),INTENT(IN)   :: N
      INTEGER(4),INTENT(IN)   :: M
      INTEGER(4),INTENT(IN)   :: L
      COMPLEX(8),INTENT(IN)   :: A(N,M)
      COMPLEX(8),INTENT(IN)   :: B(M,L)
      COMPLEX(8),INTENT(INOUT):: C(N,L)
      COMPLEX(8)              :: SUM
      COMPLEX(8),ALLOCATABLE  :: WORK(:,:)
      INTEGER(4)              :: I,J,K
!     ******************************************************************
      IF(TID) THEN
        ALLOCATE(WORK(N,M))
        WORK=A
!       ==  C=C+MATMUL(WORK,B) 
        CALL ZGEMM('N','N',N,L,M,(1.D0,0.D0),WORK,N,B,M,(1.D0,0.D0),C,N)
        DEALLOCATE(WORK)
      ELSE
!       ==  C=C+MATMUL(A,B) 
        CALL ZGEMM('N','N',N,L,M,(1.D0,0.D0),A,N,B,M,(1.D0,0.D0),C,N)
      END IF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$DYADSUMR8(LEN1,LEN2,N,PSI1,PSI2,OPERATOR)
!     ******************************************************************
!     **   OPERATOR = SUM_I |PSI1(I)><PSI2(I)|                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LEN1
      INTEGER(4),INTENT(IN) :: LEN2
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: PSI1(LEN1,N)
      REAL(8)   ,INTENT(IN) :: PSI2(LEN2,N)
      REAL(8)   ,INTENT(OUT):: OPERATOR(LEN1,LEN2)
      INTEGER(4)            :: I,J,K
      REAL(8)               :: SUM
!     ******************************************************************
#IF DEFINED(CPPVAR_BLAS_ESSL)
      CALL DGEMUL(PSI1,LEN1,'N',PSI2,LEN2,'T',OPERATOR,LEN1,LEN1,N,LEN2)
#ELSE 
!     == OPERATOR=MATMUL(PSI1,TRANSPOSE(PSI2))
      CALL DGEMM('N','T',LEN1,LEN2,N,1.D0,PSI1,LEN1,PSI2,LEN2,0.D0 &
     &          ,OPERATOR,LEN1)
#ENDIF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$DYADSUMC8(LEN1,LEN2,N,PSI1,PSI2,OPERATOR)
!     ******************************************************************
!     **   OPERATOR = SUM_I |PSI1(I)><PSI2(I)|                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LEN1
      INTEGER(4),INTENT(IN) :: LEN2
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: PSI1(LEN1,N)
      COMPLEX(8),INTENT(IN) :: PSI2(LEN2,N)
      COMPLEX(8),INTENT(OUT):: OPERATOR(LEN1,LEN2)
      INTEGER(4)            :: I,J,K
      COMPLEX(8)            :: SUM
!     ******************************************************************
#IF DEFINED(CPPVAR_BLAS_ESSL)
      CALL ZGEMUL(PSI1,LEN1,'N',PSI2,LEN2,'C',OPERATOR,LEN1,LEN1,N,LEN2)
#ELSE 
!     == OPERATOR=MATMUL(PSI1,TRANSPOSE(PSI2))
      OPERATOR(:,:)=(0.D0,0.D0)
      CALL ZGEMM('N','C',LEN1,LEN2,N,(1.D0,0.D0) &
     &          ,PSI1(:,:),LEN1,PSI2(:,:),LEN2,(0.D0,0.D0),OPERATOR,LEN1)
#ENDIF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$SCALARPRODUCTR8(TID,LEN,N1,PSI1,N2,PSI2,OVERLAP)
!     ******************************************************************
!     ** PERFORMS THE SCALAR PRODUCTS OF TWO ARRAYS OF COMPLEX        **
!     ** STATE VECTORS                                                **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TID
      INTEGER(4),INTENT(IN) :: LEN
      INTEGER(4),INTENT(IN) :: N1
      INTEGER(4),INTENT(IN) :: N2
      REAL(8)   ,INTENT(IN) :: PSI1(LEN,N1)
      REAL(8)   ,INTENT(IN) :: PSI2(LEN,N2)
      REAL(8)   ,INTENT(OUT):: OVERLAP(N1,N2)
      INTEGER(4)            :: I,J,K
      REAL(8)               :: SUM
!     ******************************************************************
      IF(TID.AND.N1.NE.N2) THEN
        CALL ERROR$MSG('PSI2 AND PSI1 DIFFER FOR TID=.TRUE.')
        CALL ERROR$STOP('LIB$SCALARPRODUCTR8')
      END IF
      IF(TID) THEN
!       ==  OVERLAP(I,J) = 0.D0*OVERLAP+1.D0*SUM_K:PSI1(K,I)*PSI1(K,J) =
        CALL DSYRK('U','T',N1,LEN,1.D0,PSI1,LEN,0.D0,OVERLAP,N1)
        DO I=1,N1
          DO J=I+1,N2
            OVERLAP(J,I)=OVERLAP(I,J)
          ENDDO
        ENDDO
      ELSE
#IF DEFINED(CPPVAR_BLAS_ESSL)
        CALL DGEMUL(PSI1,LEN,'T',PSI2,LEN,'N',OVERLAP,N1,N1,LEN,N2)
#ELSE 
        CALL DGEMM('T','N',N1,N2,LEN,1.D0,PSI1(:,:),LEN,PSI2(:,:),LEN &
     &             ,0.D0,OVERLAP,N1)
#ENDIF
      END IF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$SCALARPRODUCTC8(TID,LEN,N1,PSI1,N2,PSI2,OVERLAP)
!     ******************************************************************
!     ** PERFORMS THE SCALAR PRODUCTS OF TWO ARRAYS OF COMPLEX        **
!     ** STATE VECTORS                                                **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TID
      INTEGER(4),INTENT(IN) :: LEN
      INTEGER(4),INTENT(IN) :: N1
      INTEGER(4),INTENT(IN) :: N2
      COMPLEX(8),INTENT(IN) :: PSI1(LEN,N1)
      COMPLEX(8),INTENT(IN) :: PSI2(LEN,N2)
      COMPLEX(8),INTENT(OUT):: OVERLAP(N1,N2)
      INTEGER(4)            :: I,J,K
      COMPLEX(8)            :: SUM
!     ******************************************************************
      IF(TID.AND.N1.NE.N2) THEN
        CALL ERROR$MSG('PSI2 AND PSI1 DIFFER FOR TID=.TRUE.')
        CALL ERROR$STOP('LIB$SCALARPRODUCTC8')
      END IF
      OVERLAP(:,:)=(0.D0,0.D0)
      IF(TID) THEN
!       == ATTENTION: SCALAR FACTORS ARE SUPPOSED TO BE REAL AS THEY ARE
        CALL ZHERK('U','C',N1,LEN,1.D0,PSI1,LEN,0.D0,OVERLAP,N1)
        DO I=1,N1
          DO J=I+1,N2
            OVERLAP(J,I)=CONJG(OVERLAP(I,J))
          ENDDO
        ENDDO
      ELSE
#IF DEFINED(CPPVAR_BLAS_ESSL)
        CALL ZGEMUL(PSI1,LEN,'C',PSI2,LEN,'N',OVERLAP,N1,N1,LEN,N2)
#ELSE 
        CALL ZGEMM('C','N',N1,N2,LEN,(1.D0,0.D0),PSI1,LEN,PSI2,LEN &
     &            ,(0.D0,0.D0),OVERLAP,N1)
#ENDIF
      END IF
      RETURN
      END
!
!ZAXPY(N,ALPHA,X,INCX,Y,INCY)  Y=Y+ALPHA*X P267
!     ..................................................................
      SUBROUTINE LIB$VECTORADDC8(N,X,FAC,Y)
!     ******************************************************************
!     **  ADD TWO VECTORS                                             **
!     **  X=X+FAC*Y  WHERE X,Y ARE VECTORS AND FAC IS A SCALAR        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N
      COMPLEX(8),INTENT(INOUT) :: X(N)
      COMPLEX(8),INTENT(IN)    :: FAC
      COMPLEX(8),INTENT(IN)    :: Y(N)
      INTEGER(4)               :: I
      CHARACTER(8),PARAMETER   :: LIB='ESSL'
!     ******************************************************************
      CALL ZAXPY(N,FAC,X,1,Y,1)
      RETURN
      END
!
!DAXPY(N,ALPHA,X,INCX,Y,INCY)  Y=Y+ALPHA*X P267
!DVES(N,X,INCX,Y,INCY,X,INCZ)  Z=X-Y       P313
!     ..................................................................
      SUBROUTINE LIB$VECTORADDR8(N,X,FAC,Y)
!     ******************************************************************
!     **  ADD TWO VECTORS                                             **
!     **  X=X+FAC*Y  WHERE X,Y ARE VECTORS AND FAC IS A SCALAR        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N
      REAL(8)   ,INTENT(INOUT) :: X(N)
      REAL(8)   ,INTENT(IN)    :: FAC
      REAL(8)   ,INTENT(IN)    :: Y(N)
      INTEGER(4)               :: I
      CHARACTER(8),PARAMETER   :: LIB='ESSL'
!     ******************************************************************
      IF(FAC.EQ.-1.D0) THEN
#IF DEFINED(CPPVAR_BLAS_ESSL)
        CALL DVES(N,X,1,Y,1,X,1)
#ELSE
        CALL DAXPY(N,FAC,X,1,Y,1)
#ENDIF
      ELSE
        CALL DAXPY(N,FAC,X,1,Y,1)
      END IF
      RETURN
      END
!
!ISORT  SORTS ELEMENTS IN A SEQUENCE P904
!     .................................................................
      SUBROUTINE LIB$SORTI4(N,X)
!     ******************************************************************
!     **  SORTS THE ARRAY X IN ASCENDING ORDER                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N
      INTEGER(4),INTENT(INOUT):: X(N)
      INTEGER(4)              :: IHELP
      INTEGER(4)              :: I,J
!     ******************************************************************
#IF DEFINED(CPPVAR_BLAS_ESSL)
      CALL ISORT(X,1,N)
#ELSE
      DO I=1,N-1
        DO J=I+1,N
          IF(X(I).GT.X(J)) THEN
            IHELP=X(I)
            X(I) =X(J)
            X(J) =IHELP
          ENDIF
        ENDDO
      ENDDO
#ENDIF
      RETURN
      END

!
!*******************************************************************************
!*******************************************************************************
!****                                                                       ****
!****  EXTERNAL INTERFACES FOR FOURIER TRANSFORM CALLS                      ****
!****  THESE ROUTINES FORK INTO THE LIBRARY SPECIFIC DRIVER ROUTINES        ****
!****                                                                       ****
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$FFTADJUSTGRD(NR)
!     **************************************************************************
!     **  THIS ROUTINE RETURNS THE ALLOWED FOURIER TRANSFORM LENGTH           ** 
!     **  THAT IS EQUAL OR LARGER THAN THE LENGTH SUPPLIED, BUT               ** 
!     **  BUT OTHERWISE AS SMALL AS POSSIBLE.                                 **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT):: NR
      INTEGER(4),PARAMETER  :: MAXI=300
      LOGICAL(4),SAVE       :: TINIT=.TRUE.
      INTEGER(4),SAVE       :: COUNT
      INTEGER(4),SAVE       :: IFR(MAXI)
      INTEGER(4)            :: H,I,J,K,M
      INTEGER(4)            :: ISVAR
      REAL(8)               :: SVAR
!     **************************************************************************
      IF (TINIT) THEN
        TINIT=.FALSE.
#IF DEFINED(CPPVAR_FFT_ESSL)
!       == ALLOWED LENGTH FOR THE ESSL FFT =====================================
        COUNT=0
        OUTER: DO H=1,25
          DO I=0,2
            DO J=0,1
              DO K=0,1
                DO M=0,1
                  IF(COUNT.GE.MAXI) EXIT OUTER
                  SVAR = 2.D0**H * 3.D0**I * 5.D0**J * 7.D0**K *11.D0**M
                  IF(SVAR.GT.37748736.D0) CYCLE
                  COUNT=COUNT+1
                  IFR(COUNT)=NINT(SVAR)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO OUTER
#ELIF DEFINED(CPPVAR_FFT_FFTW)
        COUNT=0
        OUTER: DO H=1,25
          DO I=0,2
            DO J=0,1
              DO K=0,1
                DO M=0,1
                  IF(COUNT.GE.MAXI) EXIT OUTER
                  SVAR = 2.D0**H * 3.D0**I * 5.D0**J * 7.D0**K *11.D0**M
                  IF(SVAR.GT.37748736.D0) CYCLE
                  COUNT=COUNT+1
                  IFR(COUNT)=NINT(SVAR)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO OUTER
#ELSE
!       == THE GENERAL LENGTH HAS BEEN REMOVED BECAUSE PASSF AND PASSB =========
!       == DO NOT WORK =========================================================
        COUNT=0
        OUTER: DO H=1,25
          DO I=0,4
            DO J=0,2
              IF(COUNT.GE.MAXI) EXIT OUTER
              SVAR = 2.D0**H * 3.D0**I * 5.D0**J 
              IF(SVAR.GT.37748736.D0) CYCLE
              COUNT=COUNT+1
              IFR(COUNT)=NINT(SVAR)
            ENDDO
          ENDDO
        ENDDO OUTER
#ENDIF
        DO I=1,COUNT
          DO J=I+1,COUNT
            IF(IFR(I).GT.IFR(J)) THEN
              ISVAR=IFR(I)
              IFR(I)=IFR(J)
              IFR(J)=ISVAR
            END IF
          ENDDO
        ENDDO
      ENDIF
!
      DO I=2,COUNT
        IF(IFR(I).GE.NR) THEN
          NR=IFR(I)
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('REQUESTED GRIDPOINT OUT OF RANGE')
      CALL ERROR$I4VAL('NR',NR)
      CALL ERROR$STOP('PLANEWAVE$ADJUSTFFTGRD')
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$FFTC8(DIR,LEN,NFFT,X,Y)                  
!     **************************************************************************
!     **  1-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **  PACKAGES THE ESSL ROUTINE DCFT                                      **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                             **
!     **                                                                      **
!     **  USE FFTW AS STANDARD                                                **
!     **  USE FFTESSL IF ESSL IS INSTALLED                                    **
!     **  USE FFTPACK AS BACKUP IF C-ROUTINES CANNOT BE LINKED OR             **
!     **      FFTW IS NOT AVAILABLE                                           **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR !'GTOR' OR 'RTOG'
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
!     **************************************************************************
#IF DEFINED(CPPVAR_FFT_ESSL)
      CALL LIB_FFTESSL(DIR,LEN,NFFT,X,Y)                  
#ELIF DEFINED(CPPVAR_FFT_FFTW)
      CALL LIB_FFTW(DIR,LEN,NFFT,X,Y)
#ELIF DEFINED(CPPVAR_FFT_FFTW3)
      CALL LIB_FFTW3(DIR,LEN,NFFT,X,Y)
#ELIF DEFINED(CPPVAR_FFT_ACML)
      CALL LIB_ACML_FFT1DC8(DIR,LEN,NFFT,X,Y)
#ELIF DEFINED(CPPVAR_FFT_PACK)
      CALL LIB_FFTPACK(DIR,LEN,NFFT,X,Y)
#ELSE
      CALL ERROR$MSG('NO FFT PACKAGE SELECTED DURING COMPILATION')
      CALL ERROR$STOP('LIB$FFTC8')
#ENDIF

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$3DFFTC8(DIR,N1,N2,N3,X,Y)
!     **************************************************************************
!     **  3-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **    USES THE 3D FFTW ROUTINES                                         **
!     **                                        CLEMENS FOERST, 2001          **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(4)            :: DIR !'GTOR' OR 'RTOG'
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
!     **************************************************************************
#IF DEFINED(CPPVAR_FFT_ESSL)
      CALL LIB_3DFFT_ESSL(DIR,N1,N2,N3,X,Y)
#ELIF DEFINED(CPPVAR_FFT_FFTW)
      CALL LIB_3DFFTW(DIR,N1,N2,N3,X,Y)
#ELIF DEFINED(CPPVAR_FFT_FFTW3)
      CALL LIB_3DFFTW3(DIR,N1,N2,N3,X,Y)
#ELIF DEFINED(CPPVAR_FFT_ACML)
      CALL LIB_ACML_FFT3DC8(DIR,N1,N2,N3,X,Y)
#ELIF DEFINED(CPPVAR_FFT_PACK)
      CALL LIB_3DFFTPACK(DIR,N1,N2,N3,X,Y)
#ELSE
      CALL ERROR$MSG('NO FFT PACKAGE SELECTED DURING COMPILATION')
      CALL ERROR$STOP('LIB$3DFFTC8')
#ENDIF
      RETURN
      END

#IF DEFINED(CPPVAR_FFT_ACML)
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_ACML_FFT1DC8(DIR,LEN,NFFT,X,Y)                  
!     **************************************************************************
!     **  1-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **  PACKAGES THE ACML ROUTINE ZFFT1M                                    **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                             **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      INTEGER                 :: NFFT_1
      INTEGER                 :: LEN_1
      INTEGER                 :: MODE
      INTEGER                 :: INFO
      COMPLEX(8)  ,ALLOCATABLE,SAVE :: COMM(:)
      INTEGER(4)              ,SAVE :: LENPREV=0
      REAL(8)                 ,SAVE :: SCALEFORWARD
      REAL(8)                 ,SAVE :: SCALEBACKWARD
!     **************************************************************************
      IF(LEN.EQ.0) RETURN
      NFFT_1=NFFT
      LEN_1=LEN
      IF(DIR.EQ.'RTOG') THEN
        MODE=-1
      ELSE IF(DIR.EQ.'GTOR') THEN
        MODE=1
      ELSE
        CALL ERROR$MSG('INVALID VALUE OF VARIABLE "DIR"')
        CALL ERROR$MSG('CAN BE "GTOR" OR "RTOG"')
        CALL ERROR$CHVAL('DIR',DIR)
        CALL ERROR$STOP('LIB_ACML_FFT1DC8')
      END IF
! 
!     ==========================================================================
!     == INITIALIZE PLAN                                                      ==
!     ==========================================================================
      IF(LEN.NE.LENPREV) THEN
        IF(LENPREV.NE.0) DEALLOCATE(COMM)
        ALLOCATE(COMM(3*LEN+100))
        LENPREV=LEN
        Y=X
        CALL ZFFT1M(100,NFFT,LEN,Y,COMM,INFO)
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT OF ZFFT1M HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB_ACML_FFT1DC8')
        ELSE IF(INFO.GT.0) THEN
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_ACML_FFT1DC8')
        END IF
        SCALEFORWARD=1.D0/SQRT(REAL(LEN,KIND=8))
        SCALEBACKWARD=1.D0/SCALEFORWARD
      END IF
! 
!     ==========================================================================
!     == PERFORM FOURIER TRANSFORM                                            ==
!     ==========================================================================
      Y=X
      CALL ZFFT1M(MODE,NFFT,LEN,Y,COMM,INFO)
      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('I-TH ARGUMENT OF ZFFT1M HAD AN ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$STOP('LIB_ACML_FFT1DC8')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$I4VAL('I',INFO)
        CALL ERROR$STOP('LIB_ACML_FFT1DC8')
      END IF
! 
!     ==========================================================================
!     == SCALE RESULT                                                         ==
!     ==========================================================================
      IF(DIR.EQ.'RTOG') THEN
         Y=Y*SCALEFORWARD
      ELSE 
         Y=Y*SCALEBACKWARD
      END IF
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE LIB_ACML_FFT3DC8(DIR,N1,N2,N3,X,Y)
!     **************************************************************************
!     **  3-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **    USES THE 3D FFTW ROUTINES                                         **
!     **                                        CLEMENS FOERST, 2001          **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      INTEGER                 :: N1_1,N2_1,N3_1
      INTEGER                 :: MODE
      INTEGER                 :: INFO
      COMPLEX(8)  ,ALLOCATABLE,SAVE :: COMM(:)
      INTEGER(4)              ,SAVE :: N1PREV=0
      INTEGER(4)              ,SAVE :: N2PREV=0
      INTEGER(4)              ,SAVE :: N3PREV=0
      REAL(8)                 ,SAVE :: SCALEFORWARD
      REAL(8)                 ,SAVE :: SCALEBACKWARD
!     **************************************************************************
      IF(N1*N2*N3.EQ.0) RETURN
      N1_1=N1
      N2_1=N2
      N3_1=N3
      IF(DIR.EQ.'RTOG') THEN
        MODE=-1
      ELSE IF(DIR.EQ.'GTOR') THEN
        MODE=1
      ELSE
        CALL ERROR$MSG('INVALID VALUE OF VARIABLE "DIR"')
        CALL ERROR$MSG('CAN BE "GTOR" OR "RTOG"')
        CALL ERROR$CHVAL('DIR',DIR)
        CALL ERROR$STOP('LIB_ACML_3DFFTC8')
      END IF
! 
!     ==========================================================================
!     == INITIALIZE PLAN                                                      ==
!     ==========================================================================
      IF(N1.NE.N1PREV.OR.N2.NE.N2PREV.OR.N3.NE.N3PREV) THEN
        IF(N1PREV.NE.0) DEALLOCATE(COMM)
        ALLOCATE(COMM(N1*N2*N3+3*(N1+N2+N3)))
        N1PREV=N1
        N2PREV=N2
        N3PREV=N3
        SCALEFORWARD=1.D0/SQRT(REAL(N1*N2*N3,KIND=8))
        SCALEBACKWARD=1.D0/SCALEFORWARD
      END IF
! 
!     ==========================================================================
!     == PERFORM FOURIER TRANSFORM                                            ==
!     ==========================================================================
      Y=X
      CALL ZFFT3D(MODE,N1_1,N2_1,N3_1,Y,COMM,INFO)

      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('I-TH ARGUMENT OF ZFFT1M HAD AN ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$STOP('LIB_ACML_3DFFTC8')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$I4VAL('I',INFO)
        CALL ERROR$STOP('LIB_ACML_3DFFTC8')
      END IF
! 
!     ==========================================================================
!     == SCALE RESULT                                                         ==
!     ==========================================================================
      IF(DIR.EQ.'RTOG') THEN
         Y=Y*SCALEFORWARD
      ELSE 
         Y=Y*SCALEBACKWARD
      END IF
      RETURN
      END
#ELIF DEFINED(CPPVAR_FFT_FFTW)
!
!*******************************************************************************
!*******************************************************************************
!****                                                                       ****
!****  DRIVER ROUTINES FOR THE FOURIER TRANSFORMS FROM THE FFTW PACKAGE     ****
!****                                                                       ****
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_FFTW(DIR,LEN,NFFT,X,Y)                  
!     **************************************************************************
!     **  1-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **  PACKAGES THE ESSL ROUTINE DCFT                                      **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                             **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR !'GTOR' OR 'RTOG'
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      CHARACTER(4),SAVE       :: DIRSAVE=''
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      INTEGER     ,SAVE       :: ISIGN
      REAL(8)     ,SAVE       :: SCALE
      COMPLEX(8)              :: XDUMMY(LEN,NFFT)
      INTEGER(4)              :: I
      INTEGER(4),SAVE         :: NP=0
      INTEGER(4),PARAMETER    :: NPX=10 ! #(DIFFERENT FFT PLANS)
      INTEGER(8),SAVE         :: PLANS(NPX,2),PLAN=-1
      LOGICAL                 :: DEF
      INCLUDE 'FFTW_F77.I'
!     ***********  FFTW_F77.I *******************************************
!     THIS FILE CONTAINS PARAMETER STATEMENTS FOR VARIOUS CONSTANTS
!     THAT CAN BE PASSED TO FFTW ROUTINES.  YOU SHOULD INCLUDE
!     THIS FILE IN ANY FORTRAN PROGRAM THAT CALLS THE FFTW_F77
!     ROUTINES (EITHER DIRECTLY OR WITH AN #INCLUDE STATEMENT
!     IF YOU USE THE C PREPROCESSOR).
!      INTEGER,PARAMETER :: FFTW_FORWARD=-1 ! SIGN IN THE EXPONENT OF THE FORWARD FT
!      INTEGER,PARAMETER :: FFTW_BACKWARD=1 ! SIGN IN THE EXPONENT OF THE BACKWARD FT
!      INTEGER,PARAMETER :: FFTW_REAL_TO_COMPLEX=-1
!      INTEGER,PARAMETER :: FFTW_COMPLEX_TO_REAL=1
!      INTEGER,PARAMETER :: FFTW_ESTIMATE=0
!      INTEGER,PARAMETER :: FFTW_MEASURE=1
!      INTEGER,PARAMETER :: FFTW_OUT_OF_PLACE=0
!      INTEGER,PARAMETER :: FFTW_IN_PLACE=8
!      INTEGER,PARAMETER :: FFTW_USE_WISDOM=16
!      INTEGER,PARAMETER :: FFTW_THREADSAFE=128
!     CONSTANTS FOR THE MPI WRAPPERS:
!      INTEGER,PARAMETER :: FFTW_TRANSPOSED_ORDER=1
!      INTEGER,PARAMETER :: FFTW_NORMAL_ORDER=0
!      INTEGER,PARAMETER :: FFTW_SCRAMBLED_INPUT=8192
!      INTEGER,PARAMETER :: FFTW_SCRAMBLED_OUTPUT=16384
!     **************************************************************************
!
!     ==========================================================================
!     ==  INITIALIZE FFT                                                      ==
!     ==========================================================================
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE) THEN
        IF (DIR.EQ.'GTOR') THEN
          ISIGN=1
        ELSE IF (DIR.EQ.'RTOG') THEN
          ISIGN=-1
        ELSE
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('1D-FFTW')
        END IF
!
!       == FIND PLAN IN THE LIST ===============================================
        DEF=.FALSE.
        DO I=1,NP
          IF((LEN*ISIGN).EQ.PLANS(I,1)) THEN
            DEF=.TRUE.
            PLAN=PLANS(I,2)
            EXIT
          END IF
        END DO
!
!       == CREATE NEW PLAN IF NOT IN THE LIST ==================================
        IF(.NOT.DEF) THEN
          WRITE(*,*) 'LIB_FFTW: CREATE PLAN FOR: ',TRIM(DIR),ISIGN,LEN,NP
          NP=NP+1
          IF(NP.GE.NPX) NP=NPX ! ALLOW ONLY NPX PLANS
!         CALL FFTW_F77_CREATE_PLAN(PLAN,LEN,FFTW_FORWARD,FFTW_MEASURE)
          CALL FFTW_F77_CREATE_PLAN(PLANS(NP,2),LEN,ISIGN,FFTW_MEASURE)
          PLANS(NP,1)=ISIGN*LEN
          PLAN=PLANS(NP,2)
        END IF
        LENSAVE=LEN
        DIRSAVE=DIR
        SCALE=1.D0/REAL(LEN,KIND=8)
      END IF
!
!     ==========================================================================
!     ==  NOW PERFORM FFT                                                     ==
!     ==========================================================================
      XDUMMY=X
      CALL FFTW_F77(PLAN,NFFT,XDUMMY(1,1),1,LEN,Y(1,1),1,LEN)
      IF (DIR.EQ.'RTOG') THEN
        Y(:,:)=Y(:,:)*SCALE
      END IF
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE LIB_3DFFTW(DIR,N1,N2,N3,X,Y)
!     **************************************************************************
!     **  3-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **    USES THE 3D FFTW ROUTINES                                         **
!     **                                        CLEMENS FOERST, 2001          **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER                 :: DIM(3)
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      INTEGER(4)  ,SAVE       :: NP=0
      INTEGER(4),PARAMETER    :: NPX=10
      INTEGER(8)              :: PLAN
      INTEGER(8)  ,SAVE       :: PLANS(NPX,4)
      REAL(8)     ,SAVE       :: SCALE
      INTEGER     ,SAVE       :: DIMSAVE(3)=0
      CHARACTER(4),SAVE       :: DIRSAVE=''
      LOGICAL                 :: DEF
      INTEGER(4)              :: I
      INTEGER     ,SAVE       :: ISIGN
      INCLUDE 'FFTW_F77.I'
!     **************************************************************************
      DIM(1)=N1
      DIM(2)=N2
      DIM(3)=N3
      IF(DIM(1).NE.DIMSAVE(1).OR.DIM(2).NE.DIMSAVE(2).OR. &
     &  DIM(3).NE.DIMSAVE(3).OR.DIR.NE.DIRSAVE) THEN
        IF (DIR.EQ.'GTOR') THEN
          ISIGN=1
        ELSE IF (DIR.EQ.'RTOG') THEN
          ISIGN=-1
        ELSE
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('3D-FFTW')
        END IF
!
!       == FIND PLAN IN THE LIST
        DEF=.FALSE.
        DO I=1,NP
          IF((DIM(1)*ISIGN).EQ.PLANS(I,1).AND.(DIM(2)*ISIGN).EQ.PLANS(I,2)&
     &                     .AND.(DIM(3)*ISIGN).EQ.PLANS(I,3)) THEN
            DEF=.TRUE.
            PLAN=PLANS(I,4)
            EXIT
          END IF
        END DO
!
!       == CREATE NEW PLAN IF NOT IN THE LIST ==================================
        IF(.NOT.DEF) THEN
          WRITE(*,*) '3D-FFTW CREATE PLAN FOR: ', ISIGN,DIM
          NP=NP+1
          IF(NP.GE.NPX) NP=NPX ! ALLOW ONLY NPX PLANS
          CALL FFTWND_F77_CREATE_PLAN(PLAN,3,DIM,ISIGN,FFTW_ESTIMATE)
          PLANS(NP,1:3)=ISIGN*DIM
          PLANS(NP,4)=PLAN
        END IF
        DIMSAVE=DIM
        DIRSAVE=DIR
        IF(DIR.EQ.'RTOG') SCALE=1.D0/REAL(N1*N2*N3,KIND=8)
      END IF
!
!     ==========================================================================
!     ==  NOW PERFORM FFT                                                     ==
!     ==========================================================================
      CALL FFTWND_F77_ONE(PLAN,X,Y)
      IF (DIR.EQ.'RTOG') Y=Y*SCALE
      RETURN
      END
#ENDIF
#IF DEFINED(CPPVAR_FFT_FFTW3)
!
!*******************************************************************************
!**  DRIVER ROUTINES FOR THE FOURIER TRANSFORMS FROM THE FFTW3 PACKAGE        **
!**                                                                           **
!**  NOTE THAT THE CALLS TO FFTW3 ARE NOT COMPATIBLE WITH FFTW2,              **
!**  WHICH IS CALLED HERE AS FFTW                                             **
!**                                                                           **
!**               HTTP://WWW.FFTW.ORG/                                        **
!*******************************************************************************
!     ..........................................................................
      SUBROUTINE LIB_FFTW3(DIR,LEN,NFFT,X,Y)                  
!     **************************************************************************
!     **  1-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                             **
!     **                                                                      **
!     **************************************************************************
      use, intrinsic :: iso_c_binding
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      CHARACTER(4),SAVE       :: DIRSAVE=''
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      REAL(8)     ,SAVE       :: SCALE
      INTEGER     ,SAVE       :: ISIGN
      COMPLEX(8)              :: XDUMMY(LEN)
      COMPLEX(8)              :: YDUMMY(LEN)
      INTEGER(4),SAVE         :: NP=0
      INTEGER(4),PARAMETER    :: NPX=10 ! #(DIFFERENT FFT PLANS)
      type(C_PTR),SAVE        :: PLANS2(NPX),PLAN
      INTEGER(4),SAVE         :: PLANS1(NPX)
      LOGICAL                 :: DEF
      INTEGER(4)              :: I
      include 'fftw3.f03'
!     **************************************************************************
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE) THEN
        IF (DIR.EQ.'GTOR') THEN
          ISIGN=1
        ELSE IF (DIR.EQ.'RTOG') THEN
          ISIGN=-1
        ELSE
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('1D-FFTW')
        END IF
!
!       == FIND PLAN IN THE LIST ===============================================
        DEF=.FALSE.
        DO I=1,NP
          IF((LEN*ISIGN).EQ.PLANS1(I)) THEN
            DEF=.TRUE.
            PLAN=PLANS2(I)
            EXIT
          END IF
        END DO
!
!       == CREATE NEW PLAN IF NOT IN THE LIST ==================================
        IF(.NOT.DEF) THEN
          WRITE(*,*) 'LIB_FFTW: CREATE PLAN FOR: ',TRIM(DIR),ISIGN,LEN,NP
          NP=NP+1
          IF(NP.GE.NPX) NP=NPX ! ALLOW ONLY NPX PLANS
          IF(DIR.EQ.'RTOG') THEN
            PLANS2(NP) = fftw_plan_dft_1d(LEN,XDUMMY,YDUMMY,FFTW_FORWARD,&
              &IOR(FFTW_DESTROY_INPUT,IOR(FFTW_MEASURE,FFTW_UNALIGNED)))
          ELSE IF (DIR.EQ.'GTOR') THEN
            PLANS2(NP) = fftw_plan_dft_1d(LEN,XDUMMY,YDUMMY,FFTW_BACKWARD,&
              &IOR(FFTW_DESTROY_INPUT,IOR(FFTW_MEASURE,FFTW_UNALIGNED)))
          ELSE
            CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
            CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
            CALL ERROR$CHVAL('DIR',TRIM(DIR))
            CALL ERROR$STOP('LIB_FFTW3')
          END IF 

          PLANS1(NP)=ISIGN*LEN
          PLAN=PLANS2(NP)
        END IF
        LENSAVE=LEN
        DIRSAVE=DIR
      END IF
!
!     ==========================================================================
!     ==  NOW PERFORM FFT                                                     ==
!     ==========================================================================
      DO I=1,NFFT
        XDUMMY=X(:,I)
        call fftw_execute_dft(plan, XDUMMY, Y(:,I))
      enddo
!
!     ==========================================================================
!     ==  SCALE RESULT                                                        ==
!     ==========================================================================
      IF (DIR.EQ.'RTOG') THEN
        SCALE=1.D0/REAL(LEN,KIND=8)
        Y(:,:)=Y(:,:)*SCALE
      END IF
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE LIB_3DFFTW3(DIR,N1,N2,N3,X,Y)
!     **************************************************************************
!     **  3-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **    USES THE 3D FFTW ROUTINES                                         **
!     **                                        CLEMENS FOERST, 2001          **
!     **************************************************************************
      use, intrinsic :: iso_c_binding
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      type(C_PTR) :: plan
      REAL(8)     ,SAVE       :: SCALE
      include 'fftw3.f03'
!     **************************************************************************
      print*,"3dfftw3",N1,N2,N3
      IF(DIR.EQ.'RTOG') THEN
        plan = fftw_plan_dft_3d(N3,N2,N1,X,Y,FFTW_FORWARD,FFTW_ESTIMATE)
      ELSE IF (DIR.EQ.'GTOR') THEN
        plan = fftw_plan_dft_3d(N3,N2,N1,X,Y,FFTW_BACKWARD,FFTW_ESTIMATE)
      ELSE
        CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
        CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
        CALL ERROR$CHVAL('DIR',TRIM(DIR))
        CALL ERROR$STOP('LIB_3DFFTW3')
      END IF  
!
!     ==========================================================================
!     ==  EXECUTE FOURIER TRANSFORM                                           ==
!     ==========================================================================
      call fftw_execute_dft(plan, X, Y)
      call fftw_destroy_plan(plan)
!
!     ==========================================================================
!     ==  SCALE RESULT                                                        ==
!     ==========================================================================
      IF (DIR.EQ.'RTOG') THEN
        SCALE=1.D0/REAL(N1*N2*N3,KIND=8)
        Y=Y*SCALE
      END IF
      RETURN
      END
#ENDIF
!
!*******************************************************************************
!**  DRIVER ROUTINES FOR THE FOURIER TRANSFORMS FROM THE                      **
!**     ENGINEERING AND SCIENTIFIC SUBROUTINE LIBRARY (ESSL)                  **
!**                                                                           **
!**  ESSL IS A COMMERCIAL PACKAGE FROM IBM                                    **
!*******************************************************************************
!DCFT:  1-D FFT P765
#IF DEFINED(CPPVAR_FFT_ESSL)
!     ..................................................................
      SUBROUTINE LIB_FFTESSL(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      INTEGER(4)  ,PARAMETER  :: NAUX1=20000
      INTEGER(4)  ,PARAMETER  :: NAUX2=20000
      REAL(8)     ,SAVE       :: AUX2(NAUX2)
      REAL(8)     ,SAVE       :: AUX1(NAUX1)
      CHARACTER(4),SAVE       :: DIRSAVE=''
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      INTEGER(4)  ,SAVE       :: ISIGN
      REAL(8)     ,SAVE       :: SCALE
      INTEGER(4)              :: IFFT,I
!     ******************************************************************
!  
!     ==================================================================
!     == INITIALIZATION PHASE                                         ==
!     ==================================================================
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE.OR.NFFT.NE.NFFTSAVE) THEN
        IF(LEN.GT.8192) THEN
          CALL ERROR$MSG('FFT TOO LONG')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('FFT')
        END IF
        IF(DIR.EQ.'GTOR') THEN
          ISIGN=-1
          SCALE=1.D0
        ELSE IF(DIR.EQ.'RTOG') THEN
          ISIGN=1
          SCALE=1.D0/REAL(LEN,KIND=8)
        ELSE 
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('FFT')
        END IF
        CALL DCFT(1,X,1,LEN,Y,1,LEN,LEN,NFFT,ISIGN,SCALE,AUX1,NAUX1,AUX2,NAUX2)
        DIRSAVE=DIR
        LENSAVE=LEN
        NFFTSAVE=NFFT
      END IF
!  
!     ==================================================================
!     == FOURIER TRANSFORM                                            ==
!     ==================================================================
      CALL DCFT(0,X,1,LEN,Y,1,LEN,LEN,NFFT,ISIGN,SCALE,AUX1,NAUX1,AUX2,NAUX2)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB_3DFFT_ESSL(DIR,N1,N2,N3,X,Y)
!     ******************************************************************
!     **  3-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **    USES THE 3D FFT ROUTINE OF ESSL DCFT3                     **
!     **                                        PETER BLOECHL, 2001   **
!     ******************************************************************
!     **    NOT TESTED                                                **
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: DIM(3)
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      INTEGER(4)              :: NAUX       
      REAL(8)   ,ALLOCATABLE  :: AUX(:)
      INTEGER(4)              :: ISIGN
      REAL(8)                 :: SCALE
      INTEGER(4)              :: S,PSI,LAMBDA
!     ******************************************************************
      NAUX=60000
      IF(N1.GT.2048) NAUX=NAUX+4.56*N1
      IF(N3.LT.252) THEN
        IF(N2.GE.252) THEN
          S=MIN(64,N1)
          LAMBDA=(2*N2+256)*(S+4.56)
          NAUX=NAUX+LAMBDA
        END IF
      ELSE
        IF(N2.GE.252) THEN
          S=MIN(64,N1*N2)
          PSI=(2*N3+256)*(S+4.56)
          NAUX=NAUX+PSI
        ELSE
          S=MIN(64,N1*N2)
          PSI=(2*N3+256)*(S+4.56)
          S=MIN(64,N1)
          LAMBDA=(2*N2+256)*(S+4.56)
          NAUX=NAUX+MAX(PSI,LAMBDA)
        END IF
      END IF
      IF (DIR.EQ.'GTOR') THEN
        ISIGN=1
        SCALE=1.D0
      ELSE IF (DIR.EQ.'RTOG') THEN
        ISIGN=-1
        SCALE=1.D0/REAL(N1*N2*N3,KIND=8)
      ELSE
        CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
        CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
        CALL ERROR$CHVAL('DIR',TRIM(DIR))
        CALL ERROR$STOP('LIB_3DFFT_ESSL')
      END IF
      ALLOCATE(AUX(NAUX))
      CALL DCFT3(X,N1,N1*N2,Y,N1,N1*N2,N1,N2,N3,ISIGN,SCALE,AUX,NAUX)
      DEALLOCATE(AUX)
      RETURN
      END
#ENDIF
!
!
!*******************************************************************************
!**  DRIVER ROUTINES FOR THE FOURIER TRANSFORMS FROM THE                      **
!**     FFTPACK PACKAGE                                                       **
!**                                                                           **
!**  FFTPACK IS A FORTRAN SUBROUTINE LIBRARY OF FAST FOURIER TRANSFORM        **
!**  DEVELOPED AT THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH.               **
!**                                                                           **
!**   FFTPACK5 IS DISTRIBUTED UNDER THE GNU GENERAL PUBLIC LICENSE            **
!**   FROM HTTP://WWW.CISL.UCAR.EDU/CSS/SOFTWARE/FFTPACK5/INDEX.HTML          **
!**                                                                           **
!*******************************************************************************
#IF DEFINED(CPPVAR_FFT_PACK)
!     ..................................................................
      SUBROUTINE LIB_FFTPACK(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      INTEGER(4)  ,PARAMETER  :: NAUX2=2000
      REAL(8)     ,SAVE       :: AUX2(NAUX2)
      INTEGER(4)  ,SAVE       :: IAUX(30)
      REAL(8)                 :: SEQUENCE(2*LEN)
      CHARACTER(4),SAVE       :: DIRSAVE=''
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      INTEGER(4)  ,SAVE       :: ISIGN
      REAL(8)     ,SAVE       :: SCALE
      INTEGER(4)              :: IFFT,I
!     ******************************************************************
!  
!     ==================================================================
!     == INITIALIZATION PHASE                                         ==
!     ==================================================================
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE) THEN
        IF(DIR.EQ.'GTOR') THEN
          ISIGN=-1
          SCALE=1.D0
        ELSE IF(DIR.EQ.'RTOG') THEN
          ISIGN=1
          SCALE=1.D0/REAL(LEN,KIND=8)
        ELSE 
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('LIB_FFTPACK')
        END IF
        DIRSAVE=DIR
        LENSAVE=LEN
        IF(LEN.EQ.1) RETURN
        IF(NAUX2.LT.2*LEN) THEN
          CALL ERROR$MSG('AUXILIARY ARRAY TOO SMALL: INCREASE NAUX2')
          CALL ERROR$I4VAL('NAUX2',NAUX2)
          CALL ERROR$I4VAL('2*LEN',2*LEN)
          CALL ERROR$STOP('LIB_FFTPACK')
        END IF
        CALL CFFTI(LEN,AUX2,IAUX)
      END IF
!  
!     ==================================================================
!     == FOURIER TRANSFORM                                            ==
!     ==================================================================
      IF(LEN.EQ.1) RETURN
      IF(ISIGN.EQ.1) THEN
        DO IFFT=1,NFFT
          DO I=1,LEN
            SEQUENCE(2*I-1)= REAL(X(I,IFFT),KIND=8)
            SEQUENCE(2*I  )=AIMAG(X(I,IFFT))
          ENDDO
          CALL CFFTF(LEN,SEQUENCE,AUX2,IAUX)
          DO I=1,LEN
            Y(I,IFFT)=CMPLX(SEQUENCE(2*I-1),SEQUENCE(2*I),KIND=8)*SCALE
          ENDDO
        ENDDO
      ELSE
        DO IFFT=1,NFFT
          DO I=1,LEN
            SEQUENCE(2*I-1)= REAL(X(I,IFFT),KIND=8)
            SEQUENCE(2*I  )=AIMAG(X(I,IFFT))
          ENDDO
          CALL CFFTB(LEN,SEQUENCE,AUX2,IAUX)
          DO I=1,LEN
            Y(I,IFFT)=CMPLX(SEQUENCE(2*I-1),SEQUENCE(2*I),KIND=8)*SCALE
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB_3DFFTPACK(DIR,N1,N2,N3,X,Y)
!     ******************************************************************
!     **  3-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **    USES THE 3D FFT ROUTINE OF ESSL DCFT3                     **
!     **                                        PETER BLOECHL, 2001   **
!     ******************************************************************
!     **    NOT TESTED                                                **
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: DIM(3)
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      COMPLEX(8)              :: WORK1(N1*N2*N3),WORK2(N1*N2*N3)
      INTEGER(4)              :: I,J,K,IND
!     ******************************************************************
      CALL LIB_FFTPACK(DIR,N1,N2*N3,X,WORK2)
      IND=0
      DO K=1,N3
        DO J=1,N2
          DO I=1,N1
            IND=IND+1
            WORK1(J+N2*(K-1+N1*(I-1)))=WORK2(IND)
          ENDDO
        ENDDO
      ENDDO
      CALL LIB_FFTPACK(DIR,N2,N1*N3,WORK1,WORK2)
      IND=0
      DO I=1,N1
        DO K=1,N3
          DO J=1,N2
            IND=IND+1
            WORK1(K+N3*(I-1+N1*(J-1)))=WORK2(IND)
          ENDDO
        ENDDO
      ENDDO
      CALL LIB_FFTPACK(DIR,N3,N1*N2,WORK1,WORK2)
      IND=0
      DO J=1,N2
        DO I=1,N1
          DO K=1,N3
            IND=IND+1
            Y(I,J,K)=WORK2(IND)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
#ENDIF


!
!     .................................................................
      SUBROUTINE LIB$GETUSAGE(ID,VALUE)
!     *****************************************************************
!     **                                                             **
!     **  PROVIDES INFORMATION  ON THE USAGE OF SYSTEM RESOURCES     **
!     **                                                             **
!     **  USES STANDARD C-LIBRARY ROUTINE GETRUSAGE                  **
!     **                                                             **
!     **                                                             **
!     **                                                             **
!     *****************************************************************
      IMPLICIT NONE
      TYPE USG_TYPE 
        SEQUENCE
        INTEGER(4) :: UTIME(2)   ! USER TIME    (SECONDS,MICROSECONDS)
        INTEGER(4) :: STIME(2)   ! SYSTEM TIME  (SECONDS,MICROSECONDS)
        INTEGER(4) :: MAXRSS     ! MAX SIZE [KBYTE] 
        INTEGER(4) :: IXRSS      ! INTEGRAL SHARED MEMORY SIZE [KBYTE*SEC]         
        INTEGER(4) :: IDRSS      ! INTEGRAL UNSHARED DATA [KBYTE*SEC]    
        INTEGER(4) :: ISRSS      ! INTEGRAL UNSHARED STACK [KBYTE*SEC]    
        INTEGER(4) :: MINFLT     ! #(PAGE RECLAIMS)
        INTEGER(4) :: MAJFLT     ! #(PAGE FAULTS)
        INTEGER(4) :: NSWAP      ! #(SWAPPING PROCESS OUT OF MAIN)            
        INTEGER(4) :: INBLOCK    ! #(FILE INPUTS)                             
        INTEGER(4) :: OUBLOCK    ! #(FILE OUTPUTS)                            
        INTEGER(4) :: MSGSND     ! #(IPC MESSAGES SEND)                       
        INTEGER(4) :: MSGRCV     ! #(IPC MESSAGES RECEIVED)                   
        INTEGER(4) :: NSIGNALS   ! #(SIGNALS SENT)                            
        INTEGER(4) :: NVCSW      ! #(VOLUNTARY CONTEXT SWITCHES)
        INTEGER(4) :: NIVCSW     ! #(INVOLUNTARY CONTEXT SWITCHES)
      END TYPE USG_TYPE
      CHARACTER(*),INTENT(IN)  :: ID
      REAL(8)     ,INTENT(OUT) :: VALUE
      TYPE(USG_TYPE)           :: USG
      REAL(8)                  :: KBYTE=2.D0**10
      INTEGER(4)               :: GETRUSAGE
      INTEGER(4)               :: RC
      REAL(8)                  :: CPUTIME
#IF DEFINED(CPPVAR_USAGE_EXIST)
      EXTERNAL GETRUSAGE
#ENDIF
!     ******************************************************************
      USG%UTIME=(/0,0/)
      USG%STIME=(/0,0/)
      USG%MAXRSS=0     
      USG%IXRSS=0      
      USG%IDRSS=0      
      USG%ISRSS=0      
      USG%MINFLT=0     
      USG%MAJFLT=0     
      USG%NSWAP=0      
      USG%INBLOCK=0    
      USG%OUBLOCK=0    
      USG%MSGSND=0     
      USG%MSGRCV=0     
      USG%NSIGNALS=0   
      USG%NVCSW=0      
      USG%NIVCSW=0     
!     ==================================================================
!     ==================================================================
#IF DEFINED(CPPVAR_USAGE_EXIST)
      RC=GETRUSAGE(%VAL(0),USG)    ! C-ROUTINE
      IF(RC.NE.0) THEN
        CALL ERROR$MSG('ERROR CALLING MYGETRUSAGE')
        CALL ERROR$STOP('MEMORY$GET')
      END IF
#ENDIF
!
      CPUTIME=(REAL(USG%UTIME(1)+USG%STIME(1),KIND=8) &
     &        +REAL(USG%UTIME(2)+USG%STIME(2),KIND=8))*1.D-6
      CPUTIME=MAX(CPUTIME,1.D-6)
      IF(ID.EQ.'MAXMEM') THEN
        VALUE=REAL(USG%MAXRSS,KIND=8)*KBYTE
      ELSE IF(ID.EQ.'USRTIME') THEN
        VALUE=(REAL(USG%UTIME(1),KIND=8)+REAL(USG%UTIME(2),KIND=8))*1.D-6
      ELSE IF(ID.EQ.'SYSTIME') THEN
        VALUE=(REAL(USG%STIME(1),KIND=8)+REAL(USG%STIME(2),KIND=8))*1.D-6
      ELSE IF(ID.EQ.'CPUTIME') THEN
        VALUE=CPUTIME
      ELSE IF(ID.EQ.'PAGEFAULTRATE') THEN
        VALUE=REAL(USG%MAJFLT,KIND=8)/CPUTIME
      ELSE IF(ID.EQ.'SWAPRATE') THEN
        VALUE=REAL(USG%NSWAP,KIND=8)/CPUTIME
      ELSE IF(ID.EQ.'SWITCHRATE') THEN
        VALUE=REAL(USG%NVCSW+USG%NIVCSW,KIND=8)/CPUTIME
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('MEMORY$GET')
      END IF
      RETURN
      END
!
#IF DEFINED(CPPVAR_LANGEXT_XLF)
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$ERFR8(X,Y)
!     **************************************************************************
!     **  ERROR FUNCTION                                                      **
!     **    Y=2/SQRT(PI)  INT_0^X DZ EXP(-Z**2)                               **
!     **    Y=(INFINITY)=1                                                    **
!     **************************************************************************
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: Y
!     **************************************************************************
      Y=DERF(X) 
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$ERFCR8(X,Y)
!     **************************************************************************
!     **  COMPLEMENTARY ERROR FUNCTION                                        **
!     **  Y=1-ERF(X)                                                          **
!     **************************************************************************
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: Y
!     **************************************************************************
      Y=DERFC(X)
      RETURN
      END
#ELSE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$ERFR8(X,Y)
!     **************************************************************************
!     **  COPYRIGHT(C) 1996 TAKUYA OOURA                                      **
!     **  (EMAIL: OOURA@MMM.T.U-TOKYO.AC.JP).                                 **
!     **  YOU MAY USE, COPY, MODIFY THIS CODE FOR ANY PURPOSE AND             **
!     **  WITHOUT FEE. YOU MAY DISTRIBUTE THIS ORIGINAL PACKAGE.              **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: Y
      REAL(8)            :: W
      REAL(8)            :: T
      INTEGER(4)         :: K,I
      REAL(8)            :: A(0:64)
      REAL(8)            :: B(0:64)
!     **************************************************************************
      DATA (A(I), I = 0, 12) / &
     &    0.00000000005958930743D0, -0.00000000113739022964D0, &
     &    0.00000001466005199839D0, -0.00000016350354461960D0, &
     &    0.00000164610044809620D0, -0.00001492559551950604D0, &
     &    0.00012055331122299265D0, -0.00085483269811296660D0, &
     &    0.00522397762482322257D0, -0.02686617064507733420D0, &
     &    0.11283791670954881569D0, -0.37612638903183748117D0, &
     &    1.12837916709551257377D0 / 
      DATA (A(I), I = 13, 25) /                                &
     &    0.00000000002372510631D0, -0.00000000045493253732D0, &
     &    0.00000000590362766598D0, -0.00000006642090827576D0, &
     &    0.00000067595634268133D0, -0.00000621188515924000D0, &
     &    0.00005103883009709690D0, -0.00037015410692956173D0, &
     &    0.00233307631218880978D0, -0.01254988477182192210D0, &
     &    0.05657061146827041994D0, -0.21379664776456006580D0, &
     &    0.84270079294971486929D0 / 
      DATA (A(I), I = 26, 38) /                                &
     &    0.00000000000949905026D0, -0.00000000018310229805D0, &
     &    0.00000000239463074000D0, -0.00000002721444369609D0, &
     &    0.00000028045522331686D0, -0.00000261830022482897D0, &
     &    0.00002195455056768781D0, -0.00016358986921372656D0, &
     &    0.00107052153564110318D0, -0.00608284718113590151D0, &
     &    0.02986978465246258244D0, -0.13055593046562267625D0, &
     &    0.67493323603965504676D0 / 
      DATA (A(I), I = 39, 51) /                                &
     &    0.00000000000382722073D0, -0.00000000007421598602D0, &
     &    0.00000000097930574080D0, -0.00000001126008898854D0, &
     &    0.00000011775134830784D0, -0.00000111992758382650D0, &
     &    0.00000962023443095201D0, -0.00007404402135070773D0, &
     &    0.00050689993654144881D0, -0.00307553051439272889D0, &
     &    0.01668977892553165586D0, -0.08548534594781312114D0, &
     &    0.56909076642393639985D0 / 
      DATA (A(I), I = 52, 64) /                                &
     &    0.00000000000155296588D0, -0.00000000003032205868D0, &
     &    0.00000000040424830707D0, -0.00000000471135111493D0, &
     &    0.00000005011915876293D0, -0.00000048722516178974D0, &
     &    0.00000430683284629395D0, -0.00003445026145385764D0, &
     &    0.00024879276133931664D0, -0.00162940941748079288D0, &
     &    0.00988786373932350462D0, -0.05962426839442303805D0, &
     &    0.49766113250947636708D0 / 
      DATA (B(I), I = 0, 12) /                                 &
     &   -0.00000000029734388465D0,  0.00000000269776334046D0, &
     &   -0.00000000640788827665D0, -0.00000001667820132100D0, & 
     &   -0.00000021854388148686D0,  0.00000266246030457984D0, &
     &    0.00001612722157047886D0, -0.00025616361025506629D0, &
     &    0.00015380842432375365D0,  0.00815533022524927908D0, &
     &   -0.01402283663896319337D0, -0.19746892495383021487D0, & 
     &    0.71511720328842845913D0 / 
      DATA (B(I), I = 13, 25) /                                 &
     &   -0.00000000001951073787D0, -0.00000000032302692214D0,  &
     &    0.00000000522461866919D0,  0.00000000342940918551D0,  &
     &   -0.00000035772874310272D0,  0.00000019999935792654D0,  &
     &    0.00002687044575042908D0, -0.00011843240273775776D0,  &
     &   -0.00080991728956032271D0,  0.00661062970502241174D0,  &
     &    0.00909530922354827295D0, -0.20160072778491013140D0,  &
     &    0.51169696718727644908D0 / 
      DATA (B(I), I = 26, 38) /                                 &
     &    0.00000000003147682272D0, -0.00000000048465972408D0,  &
     &    0.00000000063675740242D0,  0.00000003377623323271D0,  &
     &   -0.00000015451139637086D0, -0.00000203340624738438D0,  &
     &    0.00001947204525295057D0,  0.00002854147231653228D0,  &
     &   -0.00101565063152200272D0,  0.00271187003520095655D0,  &
     &    0.02328095035422810727D0, -0.16725021123116877197D0,  &
     &    0.32490054966649436974D0 / 
      DATA (B(I), I = 39, 51) /                                &
     &    0.00000000002319363370D0, -0.00000000006303206648D0, & 
     &   -0.00000000264888267434D0,  0.00000002050708040581D0, & 
     &    0.00000011371857327578D0, -0.00000211211337219663D0, & 
     &    0.00000368797328322935D0,  0.00009823686253424796D0, & 
     &   -0.00065860243990455368D0, -0.00075285814895230877D0, & 
     &    0.02585434424202960464D0, -0.11637092784486193258D0, & 
     &    0.18267336775296612024D0 / 
      DATA (B(I), I = 52, 64) /                                &
     &   -0.00000000000367789363D0,  0.00000000020876046746D0, & 
     &   -0.00000000193319027226D0, -0.00000000435953392472D0, & 
     &    0.00000018006992266137D0, -0.00000078441223763969D0, & 
     &   -0.00000675407647949153D0,  0.00008428418334440096D0, & 
     &   -0.00017604388937031815D0, -0.00239729611435071610D0, & 
     &    0.02064129023876022970D0, -0.06905562880005864105D0, & 
     &    0.09084526782065478489D0 / 
      W = ABS(X)
      IF (W .LT. 2.2D0) THEN
          T = W * W
          K = INT(T)
          T = T - REAL(K,KIND=8)
          K = K * 13
          Y = (((((((((((( A(K    )*T + A(K+ 1)) * T + &
     &        A(K+ 2))*T + A(K+ 3))*T + A(K+ 4)) * T + &
     &        A(K+ 5))*T + A(K+ 6))*T + A(K+ 7)) * T + &
     &        A(K+ 8))*T + A(K+ 9))*T + A(K+10)) * T + &
     &        A(K+11))*T + A(K+12))*W
      ELSE IF (W .LT. 6.9D0) THEN
          K = INT(W)
          T = W - REAL(K,KIND=8)
          K = 13 * (K - 2)
          Y = (((((((((((B(K) * T + B(K + 1)) * T + &
     &        B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T + &
     &        B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T + &
     &        B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T + &
     &        B(K + 11)) * T + B(K + 12)
          Y = Y * Y
          Y = Y * Y
          Y = Y * Y
          Y = 1 - Y * Y
      ELSE
          Y = 1
      END IF
      IF(X.LT.0) Y=-Y
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$ERFCR8(X,Y)
!     **************************************************************************
!     **  COPYRIGHT(C) 1996 TAKUYA OOURA                                      **
!     **  (EMAIL: OOURA@MMM.T.U-TOKYO.AC.JP).                                 **
!     **  YOU MAY USE, COPY, MODIFY THIS CODE FOR ANY PURPOSE AND             **
!     **  WITHOUT FEE. YOU MAY DISTRIBUTE THIS ORIGINAL PACKAGE.              **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: Y
      REAL(8)     ,PARAMETER :: PA  = 3.97886080735226000D+00
      REAL(8)     ,PARAMETER :: P0  = 2.75374741597376782D-01
      REAL(8)     ,PARAMETER :: P1  = 4.90165080585318424D-01 
      REAL(8)     ,PARAMETER :: P2  = 7.74368199119538609D-01 
      REAL(8)     ,PARAMETER :: P3  = 1.07925515155856677D+00 
      REAL(8)     ,PARAMETER :: P4  = 1.31314653831023098D+00 
      REAL(8)     ,PARAMETER :: P5  = 1.37040217682338167D+00 
      REAL(8)     ,PARAMETER :: P6  = 1.18902982909273333D+00 
      REAL(8)     ,PARAMETER :: P7  = 8.05276408752910567D-01 
      REAL(8)     ,PARAMETER :: P8  = 3.57524274449531043D-01 
      REAL(8)     ,PARAMETER :: P9  = 1.66207924969367356D-02 
      REAL(8)     ,PARAMETER :: P10 =-1.19463959964325415D-01 
      REAL(8)     ,PARAMETER :: P11 =-8.38864557023001992D-02
      REAL(8)     ,PARAMETER :: P12 = 2.49367200053503304D-03 
      REAL(8)     ,PARAMETER :: P13 = 3.90976845588484035D-02 
      REAL(8)     ,PARAMETER :: P14 = 1.61315329733252248D-02 
      REAL(8)     ,PARAMETER :: P15 =-1.33823644533460069D-02 
      REAL(8)     ,PARAMETER :: P16 =-1.27223813782122755D-02 
      REAL(8)     ,PARAMETER :: P17 = 3.83335126264887303D-03 
      REAL(8)     ,PARAMETER :: P18 = 7.73672528313526668D-03 
      REAL(8)     ,PARAMETER :: P19 =-8.70779635317295828D-04 
      REAL(8)     ,PARAMETER :: P20 =-3.96385097360513500D-03 
      REAL(8)     ,PARAMETER :: P21 = 1.19314022838340944D-04 
      REAL(8)     ,PARAMETER :: P22 = 1.27109764952614092D-03
      REAL(8)                :: T,U
!     **************************************************************************
      T = PA / (PA + ABS(X))
      U = T - 0.5D0
      Y = (((((((((P22 * U + P21) * U + P20) * U + &
     &    P19) * U + P18) * U + P17) * U + P16) * U + &
     &    P15) * U + P14) * U + P13) * U + P12 
      Y = ((((((((((((Y * U + P11) * U + P10) * U + &
     &    P9) * U + P8) * U + P7) * U + P6) * U + P5) * U + &
     &    P4) * U + P3) * U + P2) * U + P1) * U + P0) * T * &
     &    EXP(-X * X)
      IF (X .LT. 0) Y = 2 - Y
      RETURN
      END
#ENDIF
!
!.......................................................................
MODULE RANDOM_MODULE
  INTEGER(8),SAVE        :: SEED=1
  INTEGER(8),PARAMETER   :: RANFAC1=22222221
  INTEGER(8),PARAMETER   :: RANFAC2=2**24
! CHOICE EXPLICIT CPPVARIABLE_XLF RNG
! INTEGER(8),SAVE        :: SEED=1_8
! INTEGER(8),PARAMETER   :: RANFAC1=44485709377909_8
! INTEGER(8),PARAMETER   :: RANFAC2=2_8**48
END MODULE RANDOM_MODULE
!
!     ..................................................................
      SUBROUTINE LIB$RANDOM(X)
!     ****************************************************************** 
!     **  RETURNS A RANDOM NUMBER                                     **
!     ******************************************************************
      USE RANDOM_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: X
!     ******************************************************************
      SEED=MODULO(RANFAC1*SEED,RANFAC2)
      X=REAL(SEED,KIND=8)/REAL(RANFAC2,KIND=8)
!     == CHOICE EXPLICIT CPPVARIABLE_XLF RNG
!     CALL RANDOM_NUMBER(Y)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$RANDOMSEED
!     ****************************************************************** 
!     **  CHOOSE A SEED  FOR THE NADOM NUMBER GENERATOR               **
!     ******************************************************************
      USE RANDOM_MODULE
      IMPLICIT NONE
!     ******************************************************************
      SEED=1
!     == CHOICE EXPLICIT CPPVARIABLE_XLF RNG
!     SEED=1_8
!     == CHOICE ORIGINAL CPPVARIABLE_XLF RNG
!     CALL RANDOM_SEED(GENERATOR=2)
!     CALL RANDOM_SEED
      RETURN
      END
!
!***********************************************************************
!***********************************************************************
!****                                                               ****
!****  TEST ROUTINES                                                ****
!****                                                               ****
!***********************************************************************
!***********************************************************************
!
!     ................................................................
      SUBROUTINE LIB$TEST()
      CALL LIB_TEST_INVERTR8()
      CALL LIB_TEST_DIAGR8()
      CALL LIB_TEST_DIAGC8()
      CALL LIB_TEST_GENERALEIGENVALUER8()
!     == COMPLEX GENERAL EIGENVALUE PROBLEM NOT IMPLEMENTED FOR ESSL
      CALL LIB_TEST_GENERALEIGENVALUEC8()
      CALL LIB_TEST_MATRIXSOLVER8()
      CALL LIB_TEST_MATRIXSOLVEC8()
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_INVERTR8()
      INTEGER(4),PARAMETER :: N=5
      REAL(8)              :: A(N,N)
      REAL(8)              :: AINV(N,N)
      REAL(8)              :: RES1(N,N)
      REAL(8)              :: RES2(N,N)
      INTEGER(4)           :: I
      LOGICAL   ,PARAMETER :: TPR=.FALSE.
!     ****************************************************************
      WRITE(*,FMT='("TEST LIB$INVERTR8")')
!
!     == MAKE INPUT DATA =============================================
      CALL RANDOM_NUMBER(A)
!
!     == WRITE INPUT DATA ============================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("A:",10F10.3)')A(I,:)
        ENDDO
      END IF
!
!     == SOLVE PROBLEM ===============================================
      CALL LIB$INVERTR8(N,A,AINV)
!
!     == WRITE OUT DATA ==============================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("AINV:",10F10.3)')AINV(I,:)
        ENDDO
      END IF
!
!     == TEST RESULT  ================================================
      RES1=MATMUL(A,AINV)
      RES2=MATMUL(AINV,A)
      DO I=1,N
        RES1(I,I)=RES1(I,I)-1.D0
        RES2(I,I)=RES2(I,I)-1.D0
      ENDDO
!      PRINT*,'TEST ',RES1
      IF(MAXVAL(ABS(RES1)+ABS(RES2)).LT.1.D-10) THEN
        WRITE(*,FMT='(T5,"OK")')
      ELSE
        WRITE(*,FMT='(T5,"FAILED")')
      END IF        
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_DIAGR8()
      INTEGER(4),PARAMETER :: N=5
      REAL(8)              :: H(N,N)
      REAL(8)              :: E(N)
      REAL(8)              :: U(N,N)
      REAL(8)              :: EMAT(N,N)
      REAL(8)              :: RES(N,N)
      INTEGER(4)           :: I
!     ****************************************************************
      WRITE(*,FMT='("TEST LIB$DIAGR8")')
      CALL RANDOM_NUMBER(H)
      H=H+TRANSPOSE(H)
!!$      DO I=1,N
!!$        WRITE(*,FMT='("H:",10F10.3)')H(I,:)
!!$      ENDDO
      CALL LIB$DIAGR8(N,H,E,U)
      EMAT=0.D0
      DO I=1,N
        EMAT(I,I)=E(I)
      ENDDO
      RES=MATMUL(H,U)-MATMUL(U,EMAT)
!      PRINT*,'TEST ',RES1
      IF(MAXVAL(ABS(RES)).LT.1.D-10) THEN
        WRITE(*,FMT='(T5,"OK")')
      ELSE
        WRITE(*,FMT='(T5,"FAILED")')
      END IF        
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_DIAGC8()
      INTEGER(4),PARAMETER :: N=5
      COMPLEX(8)           :: H(N,N)
      REAL(8)              :: E(N)
      COMPLEX(8)           :: U(N,N)
      COMPLEX(8)           :: EMAT(N,N)
      COMPLEX(8)           :: RES(N,N)
      REAL(8)              :: RE,IM
      INTEGER(4)           :: I,J
!     ****************************************************************
      WRITE(*,FMT='("TEST LIB$DIAGC8")')
      DO I=1,N
        DO J=1,N
          CALL RANDOM_NUMBER(RE)
          CALL RANDOM_NUMBER(IM)
          H(I,J)=CMPLX(RE,IM)
        ENDDO
      ENDDO
      H=H+TRANSPOSE(CONJG(H))
!!$      DO I=1,N
!!$        WRITE(*,FMT='("H:",10F10.3)')H(I,:)
!!$      ENDDO
      CALL LIB$DIAGC8(N,H,E,U)
      EMAT=(0.D0,0.D0)
      DO I=1,N
        EMAT(I,I)=CMPLX(E(I),0.D0)
      ENDDO
      RES=MATMUL(H,U)-MATMUL(U,EMAT)
!      PRINT*,'TEST ',RES
      IF(MAXVAL(ABS(RES)).LT.1.D-7) THEN
        WRITE(*,FMT='(T5,"OK")')
      ELSE
        WRITE(*,FMT='(T5,"FAILED")')
      END IF        
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_MATRIXSOLVER8()
      INTEGER(4),PARAMETER :: N=5
      INTEGER(4),PARAMETER :: M=7
      INTEGER(4),PARAMETER :: NEQ=2
      REAL(8)              :: ASQ(N,N)
      REAL(8)              :: BSQ(N,NEQ)
      REAL(8)              :: XSQ(N,NEQ)
      REAL(8)              :: A(N,M)
      REAL(8)              :: B(N,NEQ)
      REAL(8)              :: X(M,NEQ)
      INTEGER(4)           :: I
      REAL(8)              :: DEV
      REAL(8)  ,PARAMETER  :: TOL=1.D-10
      LOGICAL(4),PARAMETER :: TPR=.FALSE.
!     ****************************************************************
      WRITE(*,FMT='("LIB_TEST_MATRIXSOLVER8")')
!
!     ================================================================
!     ================================================================
!     == TEST WITH SQUARE MATRIX                                    ==
!     ================================================================
!     ================================================================
!
!     == MAKE INPUT DATA =============================================
      CALL RANDOM_NUMBER(ASQ)
      CALL RANDOM_NUMBER(BSQ)
!
!     == WRITE INPUT DATA ============================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("A:",10F10.3)')ASQ(I,:)
        ENDDO
        DO I=1,N
          WRITE(*,FMT='("B:",10F10.3)')BSQ(I,:)
        ENDDO
      END IF
!
!     == SOLVE PROBLEM ===============================================
      CALL LIB$MATRIXSOLVER8(N,N,NEQ,ASQ,XSQ,BSQ)
!
!     == WRITE OUTPUT DATA ===========================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("X:",10F10.3)')XSQ(I,:)
        ENDDO
      END IF
!
!     == TEST RESULT =================================================
      DEV=MAXVAL(ABS(MATMUL(ASQ,XSQ)-BSQ))
      IF(DEV.LT.TOL) THEN
        WRITE(*,FMT='(T5,"OK: DEV=",E10.5)')DEV
      ELSE
        WRITE(*,FMT='(T5,"FAILED: DEV=",E10.5)')DEV
      END IF        
!
!     ================================================================
!     ================================================================
!     == TEST WITH GENERAL RECTANGULAR MATRIX                       ==
!     ================================================================
!     ================================================================
!
!     == MAKE INPUT DATA =============================================
      CALL RANDOM_NUMBER(A)
      CALL RANDOM_NUMBER(B)
!
!     == WRITE INPUT DATA ============================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("A:",10F10.3)')ASQ(I,:)
        ENDDO
        DO I=1,N
          WRITE(*,FMT='("B:",10F10.3)')BSQ(I,:)
        ENDDO
      END IF
!
!     == SOLVE PROBLEM ===============================================
      CALL LIB$MATRIXSOLVER8(N,M,NEQ,A,X,B)
!
!     == WRITE OUTPUT DATA ===========================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("X:",10F10.3)')XSQ(I,:)
        ENDDO
      END IF
!
!     == TEST RESULT =================================================
      DEV=MAXVAL(ABS(MATMUL(A,X)-B))
      IF(DEV.LT.TOL) THEN
        WRITE(*,FMT='(T5,"OK: DEV=",E10.5)')DEV
      ELSE
        WRITE(*,FMT='(T5,"FAILED: DEV=",E10.5)')DEV
      END IF        
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_MATRIXSOLVEC8()
      INTEGER(4),PARAMETER :: N=5
      INTEGER(4),PARAMETER :: M=7
      INTEGER(4),PARAMETER :: NEQ=2
      COMPLEX(8)           :: ASQ(N,N)
      COMPLEX(8)           :: XSQ(N,NEQ)
      COMPLEX(8)           :: A(N,M)
      COMPLEX(8)           :: B(N,NEQ)
      COMPLEX(8)           :: X(M,NEQ)
      REAL(8)              :: RASQ(N,N)
      REAL(8)              :: RB(N,NEQ)
      REAL(8)              :: RA(N,M)
      INTEGER(4)           :: I
      REAL(8)              :: DEV
      REAL(8)  ,PARAMETER  :: TOL=1.D-10
      COMPLEX(8),PARAMETER :: CI=(0.D0,1.D0)
      LOGICAL(4),PARAMETER :: TPR=.FALSE.
!     ****************************************************************
      WRITE(*,FMT='("LIB_TEST_MATRIXSOLVEC8")')
!
!     ================================================================
!     ================================================================
!     == TEST WITH SQUARE MATRIX                                    ==
!     ================================================================
!     ================================================================
!
!     == MAKE INPUT DATA =============================================
      CALL RANDOM_NUMBER(RASQ)
      ASQ=RASQ
      CALL RANDOM_NUMBER(RASQ)
      ASQ=ASQ+RASQ*CMPLX(0.D0,1.D0)
      CALL RANDOM_NUMBER(RB)
      B=RB
      CALL RANDOM_NUMBER(RB)
      B=B+RB*CMPLX(0.D0,1.D0)

!
!     == WRITE INPUT DATA ============================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("A:",10F10.3)')ASQ(I,:)
        ENDDO
        DO I=1,N
          WRITE(*,FMT='("B:",10F10.3)')B(I,:)
        ENDDO
      END IF
!
!     == SOLVE PROBLEM ===============================================
      CALL LIB$MATRIXSOLVEC8(N,N,NEQ,ASQ,XSQ,B)
!
!     == WRITE OUTPUT DATA ===========================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("X:",10F10.3)')XSQ(I,:)
        ENDDO
      END IF
!
!     == TEST RESULT =================================================
      DEV=MAXVAL(ABS(MATMUL(ASQ,XSQ)-B))
      IF(DEV.LT.TOL) THEN
        WRITE(*,FMT='(T5,"OK: DEV=",E10.5)')DEV
      ELSE
        WRITE(*,FMT='(T5,"FAILED: DEV=",E10.5)')DEV
      END IF        
!
!     ================================================================
!     ================================================================
!     == TEST WITH GENERAL RECTANGULAR MATRIX                       ==
!     ================================================================
!     ================================================================
!
!     == MAKE INPUT DATA =============================================
      CALL RANDOM_NUMBER(RA)
      A=RA-0.5D0
      CALL RANDOM_NUMBER(RA)
      A=A+(RA-0.5D0)*CMPLX(0.D0,1.D0)
      CALL RANDOM_NUMBER(RB)
      B=(RB-0.5D0)
      CALL RANDOM_NUMBER(RB)
      B=B+(RB-0.5D0)*CMPLX(0.D0,1.D0)
!
!     == WRITE INPUT DATA ============================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("A:",10F10.3)')A(I,:)
        ENDDO
        DO I=1,N
          WRITE(*,FMT='("B:",10F10.3)')B(I,:)
        ENDDO
      END IF
!
!     == SOLVE PROBLEM ===============================================
      CALL LIB$MATRIXSOLVEC8(N,M,NEQ,A,X,B)
!
!     == WRITE OUTPUT DATA ===========================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("X:",10F10.3)')X(I,:)
        ENDDO
      END IF
!
!     == TEST RESULT =================================================
      DEV=MAXVAL(ABS(MATMUL(A,X)-B))
!      DEV=MAX(MAXVAL(ABS(REAL(MATMUL(A,X)-B))),MAXVAL(ABS(REAL(CI*(MATMUL(A,X)-B)))))
      IF(DEV.LT.TOL) THEN
        WRITE(*,FMT='(T5,"OK: DEV=",E10.5)')DEV
      ELSE
        WRITE(*,FMT='(T5,"FAILED: DEV=",E10.5)')DEV
      END IF  
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_GENERALEIGENVALUER8()
      INTEGER(4),PARAMETER :: N=5
      REAL(8)              :: H(N,N)
      REAL(8)              :: S(N,N)
      REAL(8)              :: E(N)
      REAL(8)              :: U(N,N)
      REAL(8)              :: EMAT(N,N)
      REAL(8)              :: DEV
      INTEGER(4)           :: I
      LOGICAL   ,PARAMETER :: TPR=.FALSE.
      REAL(8)   ,PARAMETER :: TOL=1.D-8
!     ****************************************************************
      WRITE(*,FMT='("LIB_TEST_GENERALEIGENVALUEC8")')
!
!     ================================================================
!     == SET UP INPUT DATA                                          ==
!     ================================================================
      CALL RANDOM_NUMBER(H)
      H=0.5D0*(H+TRANSPOSE(H))
      CALL RANDOM_NUMBER(S)
      S=MATMUL(S,TRANSPOSE(S))
      
!!$H=0.D0
!!$H(1,2)=1.D0
!!$H(2,1)=1.D0
!!$S=0.D0
!!$DO I=1,N
!!$  S(I,I)=1.D0
!!$ENDDO
!
!     ================================================================
!     == WRITE INPUT                                                ==
!     ================================================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("H:",10F10.3)')H(I,:)
        ENDDO
        WRITE(*,*)'----------------- '
        DO I=1,N
          WRITE(*,FMT='("S:",10F10.3)')S(I,:)
        ENDDO
      END IF
!
!     ================================================================
!     == SOLVE PROBLEM                                              ==
!     ================================================================
      CALL LIB$GENERALEIGENVALUER8(N,H,S,E,U)
!
!     ================================================================
!     == WRITE OUTPUT                                               ==
!     ================================================================
      IF(TPR) THEN
        WRITE(*,*)
        WRITE(*,FMT='("E:",10F10.3)')E
        WRITE(*,*)
        DO I=1,N
          WRITE(*,FMT='("U:",10F10.3)')U(I,:)
        ENDDO
      END IF
!
!     ================================================================
!     == TEST RESULT                                                ==
!     ================================================================
      EMAT=0.D0
      DO I=1,N
        EMAT(I,I)=E(I)
      ENDDO
      DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(S,MATMUL(U,EMAT))))
      IF(DEV.LT.TOL) THEN
        WRITE(*,FMT='(T5,"OK: DEV=",E10.3)')DEV
      ELSE
        WRITE(*,FMT='(T5,"FAILED: DEV=",E10.3)')DEV
      END IF        
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_GENERALEIGENVALUEC8()
      INTEGER(4),PARAMETER :: N=5
      COMPLEX(8)           :: H(N,N)
      COMPLEX(8)           :: S(N,N)
      REAL(8)              :: E(N)
      COMPLEX(8)           :: U(N,N)
      COMPLEX(8)           :: EMAT(N,N)
      REAL(8)              :: RANMAT(N,N)
      REAL(8)              :: DEV1,DEV2
      INTEGER(4)           :: I
      LOGICAL   ,PARAMETER :: TPR=.FALSE.
      REAL(8)   ,PARAMETER :: TOL=1.D-6
!     ****************************************************************
      WRITE(*,FMT='("LIB_TEST_GENERALEIGENVALUER8")')
!
!     ================================================================
!     == SET UP INPUT DATA                                          ==
!     ================================================================
      CALL RANDOM_NUMBER(RANMAT)
      H=RANMAT
      CALL RANDOM_NUMBER(RANMAT)
      H=H+RANMAT*CMPLX(0.D0,1.D0)
      H=0.5D0*(H+TRANSPOSE(CONJG(H)))
      CALL RANDOM_NUMBER(RANMAT)
      S=RANMAT
      CALL RANDOM_NUMBER(RANMAT)
      S=S+RANMAT*CMPLX(0.D0,1.D0)
      S=MATMUL(S,TRANSPOSE(CONJG(S)))
!
!     ================================================================
!     == WRITE INPUT                                                ==
!     ================================================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("H:",10("(",F10.3,",",F10.3,")"))')H(I,:)
        ENDDO
        WRITE(*,*)
        DO I=1,N
          WRITE(*,FMT='("S:",10("(",F10.3,",",F10.3,")"))')S(I,:)
        ENDDO
        WRITE(*,*)
      END IF
!
!     ================================================================
!     == SOLVE PROBLEM                                              ==
!     ================================================================
      CALL LIB$GENERALEIGENVALUEC8(N,H,S,E,U)
!
!     ================================================================
!     == WRITE OUTPUT                                               ==
!     ================================================================
      IF(TPR) THEN
        WRITE(*,FMT='("E:",10F10.3)')E
        WRITE(*,*)
        DO I=1,N
          WRITE(*,FMT='("U:",10("(",F10.3,",",F10.3,")"))')U(I,:)
        ENDDO
      END IF
!
!     ================================================================
!     == TEST RESULT                                                ==
!     ================================================================
!     == TEST EIGENVALUE PROBLEM
      EMAT=0.D0
      DO I=1,N
        EMAT(I,I)=CMPLX(E(I),0.D0)
      ENDDO
      DEV1=MAXVAL(ABS(MATMUL(H,U)-MATMUL(S,MATMUL(U,EMAT))))
!     == TEST ORTHONORMALITY
      EMAT=0.D0
      DO I=1,N
        EMAT(I,I)=CMPLX(1.D0,0.D0)
      ENDDO
      DEV2=MAXVAL(ABS(MATMUL(TRANSPOSE(CONJG(U)),MATMUL(S,U))-EMAT))
!     ==
      IF(MAX(DEV1,DEV2).LT.TOL) THEN
        WRITE(*,FMT='(T5,"OK: DEV=",2E10.3)')DEV1,DEV2
      ELSE
        WRITE(*,FMT='(T5,"FAILED: DEV=",2E10.3)')DEV1,DEV2
      END IF        
      RETURN
      END
! 
!*******************************************************************************
!*******************************************************************************
!****                                                                      *****
!****              INTERFACES FOR SPECIAL FUNCTIONS                        *****
!****                                                                      *****
!****     INTERFACE FOR BASIC SPECIAL FUNCTIONS FROM SLATEC                *****
!****     SEE SLATEC.F FOR FURTHER DETAILS                                 *****
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$DBESJ(L,X,Y)
!     **************************************************************************
!     **  BESSEL FUNCTION OF FIRST KIND                                       **
!     **  SEE SLATEC.F FOR FURTHER DETAILS                                    **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)           :: L
      REAL(8),INTENT(IN)           :: X
      REAL(8),INTENT(OUT)          :: Y
      INTEGER(4)                   :: NZ
!     ================================================================
      CALL DBESJ(X,L,1,Y,NZ)
      IF(NZ.NE.0)THEN
        Y=0.0D0
!        CALL ERROR$MSG('UNDERFLOW OR OVERFLOW')
!        CALL ERROR$STOP('LIB$DBESJ')
      END IF
      RETURN
      END SUBROUTINE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$DBESY(L,X,Y)
!     **************************************************************************
!     **  BESSEL FUNCTION OF SECOND KIND                                      **
!     **  SEE SLATEC.F FOR FURTHER DETAILS                                    **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)           :: L
      REAL(8),INTENT(IN)           :: X
      REAL(8),INTENT(OUT)          :: Y
!     ================================================================
      CALL DBESY(X,L,1,Y)
      RETURN
      END SUBROUTINE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$DBESI(L,X,Y)
!     **************************************************************************
!     **  MODIFIED BESSEL FUNCTION OF FIRST KIND                              **
!     **  SEE SLATEC.F FOR FURTHER DETAILS                                    **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)           :: L
      REAL(8),INTENT(IN)           :: X
      REAL(8),INTENT(OUT)          :: Y
      INTEGER(4)                   :: NZ
!     ================================================================
      CALL DBESI(X,L,1,1,Y,NZ)
      IF(NZ.NE.0)THEN
        Y=0.0D0
!        CALL ERROR$MSG('UNDERFLOW OR OVERFLOW')
!        CALL ERROR$STOP('LIB$DBESI')
      END IF
      RETURN
      END SUBROUTINE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$DBESK(L,X,Y)
!     **************************************************************************
!     **  MODIFIED BESSEL FUNCTION OF THIRD KIND                              **
!     **  SEE SLATEC.F FOR FURTHER DETAILS                                    **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)           :: L
      REAL(8),INTENT(IN)           :: X
      REAL(8),INTENT(OUT)          :: Y
      INTEGER(4)                   :: NZ
!     ================================================================
      CALL DBESK(X,L,1,1,Y,NZ)
      IF(NZ.NE.0)THEN
        Y=0.0D0
!        CALL ERROR$MSG('UNDERFLOW OR OVERFLOW')
!        CALL ERROR$STOP('LIB$DBESK')
      END IF
      RETURN
      END SUBROUTINE
