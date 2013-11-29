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
      INTEGER(4),parameter  :: NITEM=4
      REAL(8)               :: LOCARRAY(nitem)
      REAL(8)   ,ALLOCATABLE:: GLOBARRAY(:,:)
      REAL(8)               :: CPUTIME
!     *****************************************************************
      CALL MPE$QUERY('~',NTASKS,THISTASK)
!
!     =================================================================
!     == COLLECT ON EACH NODE                                        == 
!     =================================================================
      CALL CPU_TIME(LOCARRAY(1))
      CALL LIB$GETUSAGE('MAXMEM',LOCARRAY(2))
      CALL LIB$GETUSAGE('SWAPRATE',LOCARRAY(3))
      CALL LIB$GETUSAGE('SWITCHRATE',LOCARRAY(4))
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
        WRITE(NFIL,FMT='(30("."),T1,"MAX. MEMORY",T30,F10.5," MBYTE")')LOCARRAY(2)/MBYTE
        WRITE(NFIL,FMT='(30("."),T1,"#(SWAPS)/SEC",T30,F10.5)')LOCARRAY(3)
        WRITE(NFIL,FMT='(30("."),T1,"#(CONTEXT SWITCHES)/SEC",T30,F10.5)')LOCARRAY(4)
      END IF
      RETURN
      END
