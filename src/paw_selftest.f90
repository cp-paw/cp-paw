!
!.......................................................................
MODULE SELFTEST_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: SELFTEST                                                   **
!**                                                                   **
!**  PURPOSE: TEST WHETHER FORCES ARE CONSISTEN WITH DERIVATIVES      **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    SELFTEST$START(DENT_,LENG_,VAL_,DIFFAMP_)                      **
!**    SELFTEST$END(IDENT_,LENG_,VAL_,ENERGY_,TBACK)                  **
!**                                                                   **
!**  USAGE:                                                           **
!**    1) CALL SELFTEST$START WITH THE REFERENCE CONFIGURATION        **
!**       AND CONTINUE FORCE AND ENERGY CALCULATION WITH OUTPUT       **
!**       VALUE                                                       **
!**    2) AFTER ENERGIES AND FORCES HAVE BEEN COMPUTED,               **
!**       CALL SELFTEST$END                                           **
!**    3) IF TBACK=TRUE JUMP TO SELFTEST$BACK                         **
!**    THE PROGRAM REPEATS THE CYCLE A COUPLE OF TIME AND             **
!**    FORMS THE NUMERICAL DERIVATIVES OF THE ENERGY,                 **
!**    WHEN FINISHED THE PROGRAM PRINTS A REPORT AND RETURNS          **
!**    TBACK=.FALSE.                                                  **
!**                                                                   **
!**    EXAMPLE:                                                       **
!**      1234   SELFTEST$START('TEST1',5,POS,1.D-3)                   **
!**             .....                                                 **
!**             SELFTEST$END(TEST1,5,FORCE,ENERGY,TBACK)              **
!**             IF(TBACK) GOTO 1234                                   **
!**                                                                   **
!**  DEPENDENCIES:                                                    **
!**    ERROR                                                          **
!**    FILEHANDLER                                                    **
!**                                                                   **
!**  REMARKS:                                                         **
!**    1) AFTER SELFTEST IS FINISHED THE VARIABLES ARE UNCHANGED      **
!**    2) SELFTEST CAN ONLY DO ONE TEST AT A TIME                     **
!**    3) THE PARAMETER TON=.FALSE. HIDES ALL CALLS TO SELFTEST       **
!**       THE INTERFACE TO SWITCH TON DOES NOT YET EXIST              **
!**    4) WORKS ONLY WITH REAL*8 ARRAYS.                              **
!**       FOR COMPLEX(8) ARRAYS TRY TO RESHAPE BEFORE AND AFTER       **
!**       SELFTEST$START AND BEFORE SELFTEST$END                      **
!**                                                                   **
!***********************************************************************
INTEGER(4),PARAMETER  :: NAMP=4
INTEGER(4),PARAMETER  :: NDIS=NAMP*2+1
LOGICAL(4)            :: TON=.TRUE.
LOGICAL(4)            :: TACTIVE=.FALSE.
REAL(8)               :: ENERGY(NDIS)
REAL(8)               :: DFORCE(NAMP)
REAL(8)               :: AMPLITUDE(NDIS)
CHARACTER(64)         :: IDENT
INTEGER(4)            :: LENG=0
REAL(8)               :: DIFFAMP
REAL(8)   ,ALLOCATABLE:: XREF(:)    !(LENG)   REFERENCE COORDINATES
REAL(8)   ,ALLOCATABLE:: DIFF(:)    !(LENG)   DISTORTION PATTERN
INTEGER(4)            :: IDIS=0
END MODULE SELFTEST_MODULE
!
!     ..................................................................
      SUBROUTINE SELFTEST$START(IDENT_,LENG_,VAL_,DIFFAMP_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE SELFTEST_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)   :: IDENT_      ! IDENTIFIER OF THIS TEST
      INTEGER(4)  ,INTENT(IN)   :: LENG_       ! LENGTH OF THE COORDINATE ARRAY
      REAL(8)     ,INTENT(INOUT):: VAL_(LENG_) ! COORDINATES
      REAL(8)     ,INTENT(IN)   :: DIFFAMP_    ! SIZE OF TEH DISTORTION
      INTEGER(4)                :: I
      REAL(8)                   :: SVAR
!     ******************************************************************
      IF(.NOT.TON) RETURN
      IF(IDIS.GT.0.AND.(IDENT_.NE.IDENT.OR.LENG_.NE.LENG)) THEN
        CALL ERROR$STOP('SELFTEST')
      END IF
!     
!     ==================================================================
!     == INITIALIZE                                                   ==
!     ==================================================================
      IF(.NOT.TACTIVE) THEN
        TACTIVE=.TRUE.
!       == STORE IDENT,LENG,IDIS,XREF,DIFF,AMPLITUDE,ENERGY
        IDENT=IDENT_
        LENG=LENG_
        DIFFAMP=DIFFAMP_
!     
!       == STORE REFERENCE POSITIONS =================================
        ALLOCATE(XREF(LENG))
        XREF(:)=VAL_
!     
!       == DEFINE DISTORTION MODE = ==================================
        ALLOCATE(DIFF(LENG))
        DO I=1,LENG
          CALL RANDOM_NUMBER(SVAR)   ! INTRINSIC FORTRAN FUNCTION
          DIFF(I)=2.D0*SVAR-1.D0
        ENDDO            
!       == DEFINE DISTORTION AMPLITUDES ==============================
        DO I=1,NAMP
          AMPLITUDE(2*I-1)=DIFFAMP_*DBLE(I)
          AMPLITUDE(2*I)=-DIFFAMP_*DBLE(I)
        ENDDO
        AMPLITUDE(2*NAMP+1)=0.D0  !LAST ITERATION WITH ORIGINAL POSITIONS
      ELSE IF(IDIS.GT.NDIS) THEN
        CALL ERROR$STOP('SELFTEST')
      END IF
!     
!     ==================================================================
!     == CREATE NEW X ARRAY                                           ==
!     ==================================================================
      IDIS=IDIS+1
      VAL_(:)=XREF(:)+DIFF(:)*AMPLITUDE(IDIS)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SELFTEST$END(IDENT_,LENG_,FORCE,ENERGY_,TBACK)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE SELFTEST_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)   :: IDENT_
      INTEGER(4)  ,INTENT(IN)   :: LENG_
      REAL(8)     ,INTENT(IN)   :: FORCE(LENG_)
      REAL(8)     ,INTENT(IN)   :: ENERGY_
      LOGICAL(4)  ,INTENT(OUT)  :: TBACK
      REAL(8)                   :: DEDDIFF
      REAL(8)                   :: SVAR,SVAR2,DX
      REAL(8)                   :: X1,X2,Y1,Y2
      INTEGER(4)                :: I
      INTEGER(4)                :: NFILO
!     ******************************************************************
      TBACK=.FALSE.
      IF((.NOT.TON).OR.(.NOT.TACTIVE)) RETURN
!     
!     ==================================================================
!     ==  CHECK FOR CONSISTENCY                                       ==
!     ==================================================================
      IF(IDENT_.NE.IDENT.OR.LENG_.NE.LENG) THEN
        CALL ERROR$MSG('IDENTIFIER OR LENGTH INCONSISTENT')
        CALL ERROR$I4VAL('LENG_',LENG_)
        CALL ERROR$I4VAL('LENG',LENG)
        CALL ERROR$CHVAL('IDENT_',IDENT_)
        CALL ERROR$CHVAL('IDENT',IDENT)
        CALL ERROR$STOP('SELFTEST')
      END IF
!     
!     ==================================================================
!     ==  STORE ENERGY AND RETURN FOR THE NEXT ITERATION              ==
!     ==================================================================
      ENERGY(IDIS)=ENERGY_
      TBACK=.TRUE.
      IF(IDIS.NE.NDIS) RETURN
!     
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
      IF(IDIS.EQ.NDIS) THEN 
!       == PROJECT FORCE ON DISTORTION ARRAY ===========================
        SVAR=0.D0
        DO I=1,LENG
          SVAR=SVAR+FORCE(I)*DIFF(I)
        ENDDO
        DEDDIFF=SVAR
!     
!       == NUMERICAL DERIVATIVE ========================================   
        DO I=1,NAMP
          Y1=ENERGY(2*I-1)
          Y2=ENERGY(2*I)
          X1=AMPLITUDE(2*I-1)
          X2=AMPLITUDE(2*I)
          DFORCE(I)=(Y2-Y1)/(X2-X1)
          PRINT*,'ST-DER',X1,X2,Y1,Y2,DFORCE(I)
        ENDDO
!     
!       ================================================================
!       == PRINT OUT                                                  ==   
!       ================================================================
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,FMT='(72("="))')
        WRITE(NFILO,FMT='(72("="),T10,"  SELFTEST:",A,"  ")')IDENT
!       == ENERGY DERIVATIVE FROM ENERGIES
        WRITE(NFILO,FMT='("DX ",F10.5," DEDX(AN) ",ES20.10,"XX",ES20.10)') &
     &                0.D0,DEDDIFF,0.D0
!       == ENERGY DERIVATIVE FROM FORCES
        DO I=1,NAMP
          DX=DIFFAMP*DBLE(I)
          SVAR=DFORCE(I)
          SVAR2=(DFORCE(I)-DEDDIFF)/(DFORCE(1)-DEDDIFF)
          WRITE(NFILO,FMT='("DX ",ES10.2," DEDX(NU) ",ES20.10," XX",ES20.10)') &
     &                DX,SVAR,SVAR2
        ENDDO
        WRITE(NFILO,FMT='(72("="))')
!     
!       ================================================================
!       ==  RESET MODULE AND RETURN TBACK=.FALSE.                     ==  
!       ================================================================
        DIFFAMP=0.D0
        LENG=0
        IDENT=' '
        DEALLOCATE(DIFF)
        DEALLOCATE(XREF)
        TACTIVE=.FALSE.
        TBACK=.FALSE.
        RETURN
      END IF
      RETURN
      END
