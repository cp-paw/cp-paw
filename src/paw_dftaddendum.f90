!.......................................................................
MODULE NEWDFT_MODULE
LOGICAL(4) :: TINI=.FALSE.
CONTAINS
      SUBROUTINE NEWDFT_INITIALIZE()
      IF(TINI) RETURN
      TINI=.TRUE.
      END SUBROUTINE NEWDFT_INITIALIZE
END MODULE NEWDFT_MODULE
!     ..................................................................
      SUBROUTINE NEWDFT$EVAL1(VAL,EXC,DEXC)
!     ******************************************************************
!     **  BRIEF DESCRIPTION                                           **
!     ******************************************************************
      USE NEWDFT_MODULE
      REAL(8),INTENT(IN)  :: VAL(5)
      REAL(8),INTENT(OUT) :: EXC
      REAL(8),INTENT(OUT) :: DEXC(5)
!     ******************************************************************
!EXC=RHOT*EPSILON_XC
!VAL(1)=TOTAL DENSITY=RHOT
!VAL(2)=SPINDENSITY=RHOS= (NUP-NDOWN)
!VAL(3)= (GRAD*RHOT)**2
!VAL(4)= (GRAD*RHOS)**2
!VAL(3)= (GRAD*RHOT)*(GRAD*RHOS)
      CALL NEWDFT_INITIALIZE()
      EXC=0.D0
      DEXC(:)=0.D0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE NEWDFT$EVAL2(VAL,EXC,DEXC,D2EXC)
!     ******************************************************************
!     **  BRIEF DESCRIPTION                                           **
!     ******************************************************************
      USE NEWDFT_MODULE
      REAL(8),INTENT(IN)  :: VAL(5)
      REAL(8),INTENT(OUT) :: EXC
      REAL(8),INTENT(OUT) :: DEXC(5)
      REAL(8),INTENT(OUT) :: D2EXC(5,5)
!     ******************************************************************
      CALL NEWDFT_INITIALIZE()
      EXC=0.D0
      DEXC(:)=0.D0
      D2EXC(:,:)=0.D0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE NEWDFT$EVAL3(VAL,EXC,DEXC,D2EXC,D3EXC)
!     ******************************************************************
!     **  BRIEF DESCRIPTION                                           **
!     ******************************************************************
      USE NEWDFT_MODULE
      REAL(8),INTENT(IN)  :: VAL(5)
      REAL(8),INTENT(OUT) :: EXC
      REAL(8),INTENT(OUT) :: DEXC(5)
      REAL(8),INTENT(OUT) :: D2EXC(5,5)
      REAL(8),INTENT(OUT) :: D3EXC(5,5,5)
!     ******************************************************************
      CALL NEWDFT_INITIALIZE()
      EXC=0.D0
      DEXC(:)=0.D0
      D2EXC(:,:)=0.D0
      D3EXC(:,:,:)=0.D0
      RETURN
      END













