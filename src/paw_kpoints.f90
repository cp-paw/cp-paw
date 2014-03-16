!*******************************************************************************
!** CONTAINER FOR K-POINT-RELATED INFORMATION.                                **
!** DRIVER FOR THE PAW_BRILLOUIN OBJECT.                                      **
!** IS ALSO RESPONSIBLE FOR THE CASES                                         **
!** WHERE THE BRILLOUIN OBJECT IS NOT USED, SUCH AS ISOLATED MOLECULES        **
!*******************************************************************************
MODULE KPOINTS_MODULE
LOGICAL(4)          :: TINV           ! TIME-INVERSION SYMMETRY
INTEGER(4)          :: NKDIV(3)       ! DIVISION OF K-POINTS
INTEGER(4)          :: ISHIFT(3)      ! DISPLACEMENT OF GRID FROM GAMMA
INTEGER(4)          :: NKPT=0
REAL(8),ALLOCATABLE :: XK(:,:)  !(3,NKPT) K-POINTS IN RELATIVE COORDINATES
REAL(8),ALLOCATABLE :: WKPT(:)  ! GEOMETRIC K-POINT WEIGHTS
END MODULE KPOINTS_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$SETL4(ID,VAL)
!     **************************************************************************
      USE KPOINTS_MODULE, ONLY : TINV
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'TINV') THEN
        TINV=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('KPOINTS$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$GETL4(ID,VAL)
!     **************************************************************************
      USE KPOINTS_MODULE, ONLY : TINV
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(ID.EQ.'TINV') THEN
        VAL=TINV
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('KPOINTS$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$GETI4(ID,VAL)
!     **************************************************************************
      USE KPOINTS_MODULE, ONLY : NKPT
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(ID.EQ.'NKPT') THEN
        VAL=NKPT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('KPOINTS$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$SETI4A(ID,LEN,VAL)
!     **************************************************************************
      USE KPOINTS_MODULE, ONLY : NKDIV,ISHIFT
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: VAL(LEN)
!     **************************************************************************
      IF(ID.EQ.'DIV') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('SIZE INCONSISTENCY: LEN IS NOT 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('KPOINTS$SETI4A')
        END IF
        NKDIV=VAL
!
      ELSE IF(ID.EQ.'SHIFT') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('SIZE INCONSISTENCY: LEN IS NOT 3')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('KPOINTS$SETI4A')
        END IF
        ISHIFT=VAL
! 
     ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('KPOINTS$SETI4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$GETR8A(ID,LEN,VAL)
!     **************************************************************************
      USE KPOINTS_MODULE, ONLY : NKPT,WKPT,XK
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(OUT):: VAL(LEN)
!     **************************************************************************
      IF(ID.EQ.'XK') THEN
        IF(LEN.NE.3*NKPT) THEN
          CALL ERROR$MSG('SIZE INCONSISTENCY: LEN IS NOT 3*NKPT')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NKPT',NKPT)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('KPOINTS$GETR8A')
        END IF
        VAL=RESHAPE(XK,(/3*NKPT/))
!
      ELSE IF(ID.EQ.'WKPT') THEN
        IF(LEN.NE.NKPT) THEN
          CALL ERROR$MSG('SIZE INCONSISTENCY: LEN IS NOT NKPT')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NKPT',NKPT)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('KPOINTS$GETR8A')
        END IF
        VAL=WKPT
! 
     ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('KPOINTS$GETR8A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS_NKDIV(RBAS,RMAX,NKDIV)
!     **************************************************************************
!     **  HELPER ROUTINE                                                      **
!     **  CALCULATES THE DIVISION OF THE RECIPROCAL UNIT CELL                 **
!     **  CONSISTENT WITH A RADIUS RMAX IN REAL SPACE.                        **
!     **  THE K-POINT GRID IS CONSISTENT WITH A SUPERCELL WITH                **
!     **  THE LONGEST DIAMETER EQUAL TO RMAX. (CHECK THIS)                    **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: RMAX
      INTEGER(4),INTENT(OUT):: NKDIV(3)
      REAL(8)               :: GBAS(3,3)
      REAL(8)               :: VOL,SVAR
      INTEGER(4)            :: I
!     **************************************************************************
      CALL GBASS(RBAS,GBAS,VOL)
      DO I=1,3
        SVAR=RMAX*SQRT(DOT_PRODUCT(GBAS(:,I),GBAS(:,I))) &
     &           /DOT_PRODUCT(RBAS(:,I),GBAS(:,I))
        NKDIV(I)=INT(ABS(SVAR+1.D0))
      ENDDO
      RETURN
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE KPOINTS_KPOINTS(TINV,NKDIV,ISHIFT,NKPT,XK,WKPT)
!!$!     ******************************************************************
!!$!     **  SELECTS K-POINTS ON AN EQUISPACED GRID WITH THE PROVIDED    **
!!$!     **  DIVISIONS NKDIV. TIME INVERSION SYMMETRY IS CONSIDERED      **
!!$!     **  IF TINV=.TRUE.                                              **
!!$!     **                                                              **
!!$!     **  THE NUMBER OF K-POINTS MUST BE ESTIMATED USING KPOINTS_NKPT **
!!$!     **                                                              **
!!$!     ******************************************************************
!!$      IMPLICIT NONE
!!$      LOGICAL(4),INTENT(IN) :: TINV        ! TIME INVERSION SYMMETRY
!!$      INTEGER(4),INTENT(IN) :: NKDIV(3)    ! #(DIVISIONS)
!!$      INTEGER(4),INTENT(IN) :: ISHIFT(3)   ! DISPLACEMENT AWAY FROM GAMMA
!!$      INTEGER(4),INTENT(IN) :: NKPT        ! #(KPOINTS)
!!$      REAL(8)   ,INTENT(OUT):: XK(3,NKPT)  !KPOINT IN RELATIVE COORDINATES
!!$      REAL(8)   ,INTENT(OUT):: WKPT(NKPT)  ! GEOMETRIC K-POINT WEIGHT
!!$      REAL(8)               :: DUMMYRBAS(3,3)
!!$      INTEGER(4)            :: NKPT1
!!$!     ******************************************************************
!!$!     USE SIC UNIT CELL TO OBTAIN RELATIVE COORDINATES
!!$                               CALL TRACE$PUSH('KPOINTS_KPOINTS')
!!$      DUMMYRBAS(:,1)=(/1.D0,0.D0,0.D0/)
!!$      DUMMYRBAS(:,2)=(/0.D0,1.D0,0.D0/)
!!$      DUMMYRBAS(:,3)=(/0.D0,0.D0,1.D0/)
!!$      CALL BRILLOUIN$MSHNOSYM(TINV,DUMMYRBAS,NKDIV,ISHIFT)
!!$      CALL BRILLOUIN$GETI4('NK',NKPT1)
!!$      IF(NKPT1.NE.NKPT) THEN
!!$        CALL ERROR$STOP('KPOINTS_KPOINTS')
!!$      END IF 
!!$      CALL BRILLOUIN$GETR8A('XK',3*NKPT,XK)
!!$      CALL BRILLOUIN$GETR8A('WKPT',NKPT,WKPT)
!!$                                CALL TRACE$POP
!!$      RETURN
!!$      END
!!$!
!!$!     ..................................................................
!!$      SUBROUTINE KPOINTS_NKPT(TINV,NKDIV,ISHIFT,NKPT)
!!$!     ******************************************************************
!!$!     **                                                              **
!!$!     ******************************************************************
!!$      IMPLICIT NONE
!!$      LOGICAL(4),INTENT(IN) :: TINV
!!$      INTEGER(4),INTENT(IN) :: NKDIV(3)
!!$      INTEGER(4),INTENT(IN) :: ISHIFT(3)   ! DISPLACEMENT AWAY FROM GAMMA
!!$      INTEGER(4),INTENT(OUT):: NKPT
!!$      REAL(8)               :: DUMMYRBAS(3,3)
!!$!     ******************************************************************
!!$      DUMMYRBAS(:,1)=(/1.D0,0.D0,0.D0/)
!!$      DUMMYRBAS(:,2)=(/0.D0,1.D0,0.D0/)
!!$      DUMMYRBAS(:,3)=(/0.D0,0.D0,1.D0/)
!!$      CALL BRILLOUIN$MSHNOSYM(TINV,DUMMYRBAS,NKDIV,ISHIFT)
!!$      CALL BRILLOUIN$GETI4('NK',NKPT)
!!$      RETURN
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE KPOINTS$MAKEGRID()
!     ******************************************************************
!     **  SELECTS K-POINTS ON AN EQUISPACED GRID WITH THE PROVIDED    **
!     **  DIVISIONS NKDIV. TIME INVERSION SYMMETRY IS CONSIDERED      **
!     **  IF TINV=.TRUE.                                              **
!     **                                                              **
!     **  THE NUMBER OF K-POINTS MUST BE ESTIMATED USING KPOINTS_NKPT **
!     ******************************************************************
      USE KPOINTS_MODULE, ONLY : TINV,NKDIV,ISHIFT,NKPT,XK,WKPT
      IMPLICIT NONE
      REAL(8)               :: DUMMYRBAS(3,3)
      INTEGER(4) :: IKPT
!     ******************************************************************
                               CALL TRACE$PUSH('KPOINTS$MAKEGRID')
!     USE SIC UNIT CELL TO OBTAIN RELATIVE COORDINATES
      DUMMYRBAS(:,1)=(/1.D0,0.D0,0.D0/)
      DUMMYRBAS(:,2)=(/0.D0,1.D0,0.D0/)
      DUMMYRBAS(:,3)=(/0.D0,0.D0,1.D0/)
      CALL BRILLOUIN$MSHNOSYM(TINV,DUMMYRBAS,NKDIV,ISHIFT)
      CALL BRILLOUIN$GETI4('NK',NKPT)
      ALLOCATE(XK(3,NKPT))
      ALLOCATE(WKPT(NKPT))
      CALL BRILLOUIN$GETR8A('XK',3*NKPT,XK)
      CALL BRILLOUIN$GETR8A('WKPT',NKPT,WKPT)
                                CALL TRACE$POP
      RETURN
      END






