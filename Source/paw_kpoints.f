!
!     ..................................................................
      SUBROUTINE kpoints_nkdiv(rbas,rmax,nkdiv)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      implicit none
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: RMAX
      integer(4),INTENT(OUT):: NKDIV(3)
      INTEGER(4)            :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
!     ******************************************************************
      CALL BOXSPH(RBAS,0.D0,0.D0,0.D0,RMAX &
     &           ,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
      NKDIV(1)=MAX1-MIN1+1
      NKDIV(2)=MAX2-MIN2+1
      NKDIV(3)=MAX3-MIN3+1
      RETURN
      end
!
!     ..................................................................
      SUBROUTINE KPOINTS_KPOINTS(tinv,nkdiv,nkpt,XK,wkpt)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      implicit none
      logical(4),intent(in) :: tinv
      integer(4),INTENT(IN) :: nkdiv(3)
      integer(4),INTENT(IN) :: nkpt
      real(8)   ,INTENT(OUT):: xk(3,nkpt)
      real(8)   ,INTENT(OUT):: wkpt(nkpt)
      INTEGER(4)            :: I1,I2,I3
      INTEGER(4)            :: N1,N2,N3
      REAL(8)               :: T1FAC,T2FAC,T3FAC
      REAL(8)               :: T1,T2,T3
      REAL(8)               :: DET
      logical(4)            :: tsymm1,tsymm2,tsymm3
      logical(4)            :: tzero1,tzero2
      real(8)               :: fac
      integer(4)            :: ikpt
!     ******************************************************************
      N1=NKDIV(1)
      N2=NKDIV(2)
      N3=NKDIV(3)
      T1FAC=1.D0/REAL(N1,KIND=8)
      T2FAC=1.D0/REAL(N2,KIND=8)
      T3FAC=1.D0/REAL(N3,KIND=8)
      FAC=T1FAC*T2FAC*T3FAC      
      IKPT=0
      DO I1=1,N1
        T1=DBLE(I1-1)*T1FAC
        IF(T1.GT.0.5D0)T1=-1.D0+T1
        IF(TINV) THEN
          TSYMM1=TINV.AND.(T1.EQ.0.D0.OR.T1.EQ.0.5D0)
          IF(T1.LT.0.D0) CYCLE
          TZERO1=T1.EQ.0.D0
        END IF
        DO I2=1,N2
          T2=DBLE(I2-1)*T2FAC
          IF(T2.GT.0.5D0)T2=-1.D0+T2
          IF(TINV) THEN
            TSYMM2=TSYMM1.AND.(T2.EQ.0.D0.OR.T2.EQ.0.5D0)
            IF(TZERO1.AND.T2.LT.0.D0) CYCLE
            TZERO2=TZERO1.AND.T2.EQ.0.D0
          END IF
          DO I3=1,N3
            T3=DBLE(I3-1)*T3FAC
            IF(T3.GT.0.5D0)T3=-1.D0+T3
            IF(TINV) THEN
              TSYMM3=TSYMM2.AND.(T3.EQ.0.D0.OR.T3.EQ.0.5D0)
              IF(TZERO2.AND.T3.LT.0.D0) CYCLE
            END IF
            IKPT=IKPT+1
            IF(IKPT.GT.NKPT) THEN
              CALL ERROR$MSG('#(K-POINTS) EXCEEDED')
              CALL ERROR$I4VAL('NKPT',NKPT)
              CALL ERROR$I4VAL('N1',N1)
              CALL ERROR$I4VAL('N2',N2)
              CALL ERROR$I4VAL('N3',N3)
              CALL ERROR$STOP('KPOINTS$KPOINTS')
            END IF
            WKPT(IKPT)=FAC
            IF(TINV) THEN
               IF(.NOT.TSYMM3)WKPT(IKPT)=2.D0*WKPT(IKPT)  !avoid double couunting
            END IF
            XK(1,IKPT)=T1
            XK(2,IKPT)=T2
            XK(3,IKPT)=T3
          ENDDO
        ENDDO
      ENDDO
      IF(IKPT.NE.NKPT) THEN
        CALL ERROR$MSG('INCONSISTENT #(K-POINTS)')
        CALL ERROR$I4VAL('NKPT ANTICIPATED',NKPT)
        CALL ERROR$I4VAL('NKPT ACTUAL     ',IKPT)
        CALL ERROR$STOP('KPOINTS$KPOINTS')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE KPOINTS_nkpt(tinv,nkdiv,nkpt)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      implicit none
      logical(4),intent(in) :: tinv
      integer(4),INTENT(IN) :: nkdiv(3)
      integer(4),INTENT(out):: nkpt
      INTEGER(4)            :: I1,I2,I3
      INTEGER(4)            :: N1,N2,N3
      REAL(8)               :: T1FAC,T2FAC,T3FAC
      REAL(8)               :: T1,T2,T3
      REAL(8)               :: DET
!     ******************************************************************
      N1=NKDIV(1)
      N2=NKDIV(2)
      N3=NKDIV(3)
      T1FAC=1.D0/REAL(N1,KIND=8)
      T2FAC=1.D0/REAL(N2,KIND=8)
      T3FAC=1.D0/REAL(N3,KIND=8)
      nKPT=0
      DO I1=1,N1
        T1=DBLE(I1-1)*T1FAC
        IF(T1.GT.0.5D0)T1=-1.D0+T1
        IF(TINV.AND.T1.LT.0.D0) CYCLE 
        DO I2=1,N2
          T2=DBLE(I2-1)*T2FAC
          IF(T2.GT.0.5D0)T2=-1.D0+T2
          IF(TINV.AND.T1.EQ.0.D0.AND.T2.LT.0.D0) CYCLE 
          DO I3=1,N3
            T3=DBLE(I3-1)*T3FAC
            IF(T3.GT.0.5D0)T3=-1.D0+T3
            IF(TINV.AND.T1.EQ.0.D0.AND.T2.EQ.0.AND.T3.LT.0.D0) CYCLE 
            nKPT=nKPT+1
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END






