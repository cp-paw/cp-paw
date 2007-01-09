!
!     ..................................................................
      SUBROUTINE KPOINTS_NKDIV(RBAS,RMAX,NKDIV)
!     ******************************************************************
!     **  CALCULATES THE DIVISION OF THE RECIPROCAL UNIT CELL         **
!     **  CONSISTENT WITH A RADIUS RMAX IN REAL SPACE.                **
!     **  THE K-POINT GRID IS CONSISTENT WITH A SUPERCELL WITH        **
!     **  THE LONGEST DIAMETER EQUAL TO RMAX. (CHECK THIS)            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: RMAX
      INTEGER(4),INTENT(OUT):: NKDIV(3)
      REAL(8)               :: GBAS(3,3)
      REAL(8)               :: VOL,SVAR
      INTEGER(4)            :: I
!     ******************************************************************
      CALL GBASS(RBAS,GBAS,VOL)
      DO I=1,3
        SVAR=RMAX*SQRT(DOT_PRODUCT(GBAS(:,I),GBAS(:,I))) &
     &           /DOT_PRODUCT(RBAS(:,I),GBAS(:,I))
        NKDIV(I)=INT(ABS(SVAR+1.D0))
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE KPOINTS_KPOINTS(TINV,NKDIV,ISHIFT,NKPT,XK,WKPT)
!     ******************************************************************
!     **  SELECTS K-POINTS ON AN EQUISPACED GRID WITH THE PROVIDED    **
!     **  DIVISIONS NKDIV. TIME INVERSION SYMMETRY IS CONSIDERED      **
!     **  IF TINV=.TRUE.                                              **
!     **                                                              **
!     **  THE NUMBER OF K-POINTS MUST BE ESTIMATED USING KPOINTS_NKPT **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TINV        ! TIME INVERSION SYMMETRY
      INTEGER(4),INTENT(IN) :: NKDIV(3)    ! #(DIVISIONS)
      INTEGER(4),INTENT(IN) :: ISHIFT(3)   ! DISPLACEMENT AWAY FROM GAMMA
      INTEGER(4),INTENT(IN) :: NKPT        ! #(KPOINTS)
      REAL(8)   ,INTENT(OUT):: XK(3,NKPT)  !KPOINT IN RELATIVE COORDINATES
      REAL(8)   ,INTENT(OUT):: WKPT(NKPT)  ! GEOMETRIC K-POINT WEIGHT
      REAL(8)               :: DUMMYRBAS(3,3)
      INTEGER(4)            :: NKPT1
!     ******************************************************************
!     USE SIC UNIT CELL TO OBTAIN RELATIVE COORDINATES
                               CALL TRACE$PUSH('KPOINTS_KPOINTS')
      DUMMYRBAS(:,1)=(/1.D0,0.D0,0.D0/)
      DUMMYRBAS(:,2)=(/0.D0,1.D0,0.D0/)
      DUMMYRBAS(:,3)=(/0.D0,0.D0,1.D0/)
      CALL BRILLOUIN$MSHNOSYM(TINV,DUMMYRBAS,NKDIV,ISHIFT)
      CALL BRILLOUIN$GETI4('NK',NKPT1)
      IF(NKPT1.NE.NKPT) THEN
        CALL ERROR$STOP('KPOINTS_KPOINTS')
      END IF 
      CALL BRILLOUIN$GETR8A('XK',3*NKPT,XK)
      CALL BRILLOUIN$GETR8A('WKPT',NKPT,WKPT)
                                CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE KPOINTS_NKPT(TINV,NKDIV,ISHIFT,NKPT)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TINV
      INTEGER(4),INTENT(IN) :: NKDIV(3)
      INTEGER(4),INTENT(IN) :: ISHIFT(3)   ! DISPLACEMENT AWAY FROM GAMMA
      INTEGER(4),INTENT(OUT):: NKPT
      REAL(8)               :: DUMMYRBAS(3,3)
!     ******************************************************************
      DUMMYRBAS(:,1)=(/1.D0,0.D0,0.D0/)
      DUMMYRBAS(:,2)=(/0.D0,1.D0,0.D0/)
      DUMMYRBAS(:,3)=(/0.D0,0.D0,1.D0/)
      CALL BRILLOUIN$MSHNOSYM(TINV,DUMMYRBAS,NKDIV,ISHIFT)
      CALL BRILLOUIN$GETI4('NK',NKPT)
      RETURN
      RETURN
      END






