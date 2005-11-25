!
!     ..................................................................
      SUBROUTINE KPOINTS_NKDIV(RBAS,RMAX,NKDIV)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: RMAX
      INTEGER(4),INTENT(OUT):: NKDIV(3)
      INTEGER(4)            :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
!     ******************************************************************
      CALL BOXSPH(RBAS,0.D0,0.D0,0.D0,RMAX &
     &           ,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
      NKDIV(1)=MAX1-MIN1+1
      NKDIV(2)=MAX2-MIN2+1
      NKDIV(3)=MAX3-MIN3+1
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
      INTEGER(4)            :: I1,I2,I3
      INTEGER(4)            :: I1INV,I2INV,I3INV
      INTEGER(4)            :: N1,N2,N3
      INTEGER(4)            :: IKPT
      INTEGER(4)            :: IND,INDINV
      REAL(8)               :: WGHT,WGHT0
!     ******************************************************************
      N1=NKDIV(1)
      N2=NKDIV(2)
      N3=NKDIV(3)
      WGHT0=1.D0/REAL(N1*N2*N3,KIND=8)
      IKPT=0
      DO I1=1,N1
        DO I2=1,N2
          DO I3=1,N3
!           ===========================================================
!           ==  CHECK INVERSION SYMMETRY                             ==
!           ==  EACH POINT IN THE CELL IS DEFINED BY A SINGLE INTEGER==
!           ==  LABEL IND=1+[I1-1+N1*(I2-1+N2*I3-1)]                 ==
!           ===========================================================
            IF(TINV) THEN
              I1INV=MOD(N1-(I1+ISHIFT(1)-1),N1)+1
              I2INV=MOD(N2-(I2+ISHIFT(2)-1),N2)+1
              I3INV=MOD(N3-(I3+ISHIFT(3)-1),N3)+1
              IND   =I1   -1+N1*(I2   -1+N2*(I3   -1))
              INDINV=I1INV-1+N1*(I2INV-1+N2*(I3INV-1))
              IF(IND.EQ.INDINV) THEN  ! SPECIAL POINT (EQUIV TO ITS INVERSE)
                WGHT=WGHT0                 
              ELSE IF(IND.GT.INDINV) THEN 
                CYCLE ! INVERSE IMAGE IS NOT COUNTED
              ELSE
                WGHT=2.D0*WGHT0 ! GENERAL K-POINT GETS THE WEIGHT OF ITS IMAGE
              END IF
            ELSE
              WGHT=WGHT0
            END IF
            IKPT=IKPT+1
            IF(IKPT.GT.NKPT) CYCLE
            XK(1,IKPT)=REAL(2*(I1-1)+ISHIFT(1),KIND=8)/REAL(2*N1,KIND=8)
            XK(2,IKPT)=REAL(2*(I2-1)+ISHIFT(2),KIND=8)/REAL(2*N2,KIND=8)
            XK(3,IKPT)=REAL(2*(I3-1)+ISHIFT(3),KIND=8)/REAL(2*N3,KIND=8)
            WKPT(IKPT)=WGHT
!write(*,fmt='(i5,3f10.5,5x,f10.4)')ikpt,xk(:,ikpt),wkpt(ikpt)
          ENDDO
        ENDDO
      ENDDO   
      IF(IKPT.NE.NKPT) THEN
        CALL ERROR$MSG('INCONSISTENT #(K-POINTS)')
        CALL ERROR$I4VAL('NKPT ANTICIPATED',NKPT)
        CALL ERROR$I4VAL('NKPT ACTUAL     ',IKPT)
        CALL ERROR$STOP('KPOINTS$KPOINTS')
      END IF
!call error$stop('forced stop')
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
      INTEGER(4)            :: I1,I2,I3
      INTEGER(4)            :: I1INV,I2INV,I3INV
      INTEGER(4)            :: N1,N2,N3
      INTEGER(4)            :: IKPT
      INTEGER(4)            :: IND,INDINV
!     ******************************************************************
!     ******************************************************************
      N1=NKDIV(1)
      N2=NKDIV(2)
      N3=NKDIV(3)
      NKPT=0
      DO I1=1,N1
        DO I2=1,N2
          DO I3=1,N3
!           ===========================================================
!           ==  CHECK INVERSION SYMMETRY                             ==
!           ==  EACH POINT IN THE CELL IS DEFINED BY A SINGLE INTEGER==
!           ==  LABEL IND=1+[I1-1+N1*(I2-1+N2*I3-1)]                 ==
!           ===========================================================
            IF(TINV) THEN
              I1INV=MOD(N1-(I1+ISHIFT(1)-1),N1)+1
              I2INV=MOD(N2-(I2+ISHIFT(2)-1),N2)+1
              I3INV=MOD(N3-(I3+ISHIFT(3)-1),N3)+1
              IND   =I1   -1+N1*(I2   -1+N2*(I3   -1))
              INDINV=I1INV-1+N1*(I2INV-1+N2*(I3INV-1))
              IF(IND.GT.INDINV) THEN 
                CYCLE ! INVERSE IMAGE IS NOT COUNTED
              ENDIF
            END IF
            NKPT=NKPT+1
          ENDDO
        ENDDO
      ENDDO   
      RETURN
      END






