!**********************************************************************
!**********************************************************************
!**                                                                  **
!**  NAME: USAGE                                                     **
!**                                                                  **
!**  PURPOSE: REPORTS THE USAGE OF THE SYSTEM RESOURCES              **
!**                                                                  **
!**  METHODS                                                         **
!**    USAGE$GET(ID,VAL)                       X                      **
!**    USAGE$REPORT(NFIL)                                            **
!**                                                                  **
!**  REMARKS: USES C-FUNCTION "GETRUSAGE"  SEE AIX INFO AND          **
!**    /USR/INCLUDE/SYS/RESOURCE.H FOR FURTHER INFORMATION           **
!**                                                                  **
!**********************************************************************
!**********************************************************************
!     .................................................................
      SUBROUTINE USAGE$GET(ID,VAL)
!     *****************************************************************
!     **                                                             **
!     **  PROVIDES INFORMATION  ON THE USAGE OF SYSTEM RESOURCES     **
!     **                                                             **
!     **                                                             **
!     *****************************************************************
      USE MPE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
!     *****************************************************************
      CALL LIB$GETUSAGE(ID,VAL)
      RETURN
      END
!     .................................................................
      SUBROUTINE USAGE$REPORT(NFIL)
!     *****************************************************************
!     **                                                             **
!     **  PROVIDES INFORMATION  ON THE USAGE OF SYSTEM RESOURCES     **
!     **                                                             **
!     **                                                             **
!     *****************************************************************
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)               :: MBYTE=2.D0**20
      INTEGER(4)            :: NTASKS,THISTASK
      INTEGER(4)            :: NITEM=6
      REAL(8)   ,ALLOCATABLE:: LOCARRAY(:)
      REAL(8)   ,ALLOCATABLE:: GLOBARRAY(:,:)
      REAL(8)               :: CPUTIME
!     *****************************************************************
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      ALLOCATE(LOCARRAY(NITEM))
!
!     =================================================================
!     == COLLECT ON EACH NODE                                        == 
!     =================================================================
      CALL LIB$GETUSAGE('CPUTIME',LOCARRAY(1))
      CALL LIB$GETUSAGE('USRTIME',LOCARRAY(2))
      CALL LIB$GETUSAGE('SYSTIME',LOCARRAY(3))
      CALL LIB$GETUSAGE('MAXMEM',LOCARRAY(4))
      CALL LIB$GETUSAGE('SWAPRATE',LOCARRAY(5))
      CALL LIB$GETUSAGE('SWITCHRATE',LOCARRAY(6))
!
!     =================================================================
!     == AVERAGE OVER ALL NODES                                      == 
!     =================================================================
      ALLOCATE(GLOBARRAY(NITEM,NTASKS))
      CALL MPE$GATHER('~',1,LOCARRAY,GLOBARRAY)
      LOCARRAY(:)=SUM(GLOBARRAY,DIM=2)/REAL(NTASKS,KIND=8)
      DEALLOCATE(GLOBARRAY)
!
!     =================================================================
!     ==  NOW REPORT                                                 == 
!     =================================================================
      IF(THISTASK.EQ.1) THEN
        CPUTIME=LOCARRAY(1)
        WRITE(NFIL,FMT='(/"USAGE OF SYSTEM RESOURCES"/25("="))')
        WRITE(NFIL,FMT='(30("."),T1,"#(PROCESSORS)",T30,I10)')NTASKS
        WRITE(NFIL,FMT='(30("."),T1,"CPU TIME PER CPU",T30,F10.2," SEC")')LOCARRAY(1)
        WRITE(NFIL,FMT='(30("."),T1,"%(USER TIME)",T30,F10.5)')LOCARRAY(2)/CPUTIME
        WRITE(NFIL,FMT='(30("."),T1,"%(SYSTEM TIME)",T30,F10.5)')LOCARRAY(3)/CPUTIME
        WRITE(NFIL,FMT='(30("."),T1,"MAX. MEMORY",T30,F10.5," MBYTE")')LOCARRAY(4)/MBYTE
        WRITE(NFIL,FMT='(30("."),T1,"#(SWAPS)/SEC",T30,F10.5)')LOCARRAY(5)
        WRITE(NFIL,FMT='(30("."),T1,"#(CONTEXT SWITCHES)/SEC",T30,F10.5)')LOCARRAY(6)
      END IF
      DEALLOCATE(LOCARRAY)
      RETURN
      END
