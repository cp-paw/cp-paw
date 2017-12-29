!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POTENTIAL_PAIRP(RBAS,NSP,NAT,ISPECIES,LMRXX,LMRX,RCSM,RCBG &
     &                          ,R0,FORCE,QLM,VQLM,EPAIR,STRESS)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES SELF ENERGY AND PAIR INTERACTION OF                      **
!     **  COMPENSATION CHARGES                                                **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1991 ************
      USE MPE_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: RBAS(3,3)         !LATTICE VECTORS
      INTEGER(4),INTENT(IN)  :: NSP               ! #(ATOM TYPES)
      INTEGER(4),INTENT(IN)  :: NAT               ! #(ATOMS)
      INTEGER(4),INTENT(IN)  :: LMRXX             !
      INTEGER(4),INTENT(IN)  :: LMRX(NSP)         !
      INTEGER(4),INTENT(IN)  :: ISPECIES(NAT)     !
      REAL(8)   ,INTENT(IN)  :: RCSM(NSP)         !(NSP)
      REAL(8)   ,INTENT(IN)  :: RCBG(NSP)         !(NSP)
      REAL(8)   ,INTENT(INOUT)  :: R0(3,NAT)         ! ATOMIC POSITIONS
      REAL(8)   ,INTENT(OUT) :: FORCE(3,NAT)      ! FORCES
      REAL(8)   ,INTENT(INOUT)  :: QLM(LMRXX,NAT)    ! MULTIPOLE MOMENTS
      REAL(8)   ,INTENT(OUT) :: VQLM(LMRXX,NAT)   ! DE/DQLM
      REAL(8)   ,INTENT(OUT) :: EPAIR
      REAL(8)   ,INTENT(OUT) :: STRESS(3,3)
      REAL(8)   ,PARAMETER   :: RMAXBYRCBG=6.D0 
      REAL(8)                :: RMAX ! CUTOFF RADIUS FOR PAIR POTENTIAL 
      INTEGER(4)             :: NFILO
      INTEGER(4)             :: NTASKS,THISTASK
      INTEGER(4)             :: LRXX     
      INTEGER(4)             :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
      INTEGER(4)             :: I1,I2,I3,I,J,LM
      INTEGER(4)             :: ISP,ISP1,ISP2,IAT,IAT1,IAT2
      REAL(8)                :: T1,T2,T3,DIS,DIS2,DX,DY,DZ
      REAL(8)                :: EPAIR1,ESELF1,ESELF0,EPAIR0
      REAL(8)                :: DR(3)
      REAL(8)                :: FORCE1(3)
      REAL(8)                :: VQLM1(LMRXX)
      REAL(8)                :: VQLM2(LMRXX)
      REAL(8)                :: SVAR
      REAL(8)   ,ALLOCATABLE :: HS(:,:)     !(0:2*LRXX+1,0:LRXX)
      REAL(8)   ,ALLOCATABLE :: HB(:,:)     !(0:2*LRXX+1,0:LRXX)
      REAL(8)   ,ALLOCATABLE :: A(:,:,:)    !(LRXX+1,2*LRXX+2,LRXX)
!     ==  VARIABLES FOR SELF-TEST ==============================================
      LOGICAL(4),PARAMETER   :: TTEST_V=.FALSE.
      LOGICAL(4),PARAMETER   :: TTEST_FORCE=.FALSE.
      LOGICAL(4)             :: TBACK
!     **************************************************************************
                             CALL TRACE$PUSH('PAIRP')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      RMAX=MAXVAL(RCBG(:))*RMAXBYRCBG
      EPAIR=0.D0
      VQLM(:,:)=0.D0
      FORCE(:,:)=0.D0
!
!     ==========================================================================
!     == ALLOCATE ARRAYS                                                      ==
!     ==========================================================================
      SVAR=SQRT(REAL(LMRXX-1,KIND=8))
      LRXX=INT(SVAR+1.D-12)
      ALLOCATE(HS(0:2*LRXX+1,0:LRXX))
      ALLOCATE(HB(0:2*LRXX+1,0:LRXX))
      ALLOCATE(A(LRXX+1,2*LRXX+2,LRXX))
!
!     ==========================================================================
!     == SELFTEST                                                             ==
!     ==========================================================================
 1000 CONTINUE
      IF(TTEST_V.AND.TTEST_FORCE) THEN
        CALL ERROR$MSG('CANNNOT TEST FORCE AND POTENTIAL AT THE SAME TIME')
        CALL ERROR$STOP('PAIRP')
      END IF
      IF(TTEST_V) THEN
        CALL SELFTEST$START('PAIRP-QLM',LMRXX*NAT,QLM,1.D-3)
        VQLM(:,:)=0.D0
      END IF
      IF(TTEST_FORCE) THEN
        CALL SELFTEST$START('PAIRP-FORCE',3*NAT,R0,1.D-4)
        FORCE(:,:)=0.D0
      END IF
!
!     ==========================================================================
!     ==   H(L=0,N)=1/R * D/DR**2N * R * H(L=0,N=0)                           ==
!     ==   H(L,N)=SUM(J): R**L * A(J,L+1,N) * (-2*R)**(2J-2)*EXP(-R**2)       ==
!     ==========================================================================
      CALL PAIRP_GETA(LRXX,A)
!
!     ==========================================================================
!     ==  CALCULATE SELF ENERGY OF PSEUDOCHARGES                              ==
!     ==========================================================================
      ESELF0=0.D0
      DO IAT=THISTASK,NAT,NTASKS
        ISP=ISPECIES(IAT)
        CALL PAIRP_ESELF(RCSM(ISP),RCBG(ISP),LMRX(ISP) &
     &                 ,QLM(1,IAT),VQLM(1,IAT),ESELF1,A,LRXX)
        ESELF0=ESELF0+ESELF1
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE PAIR TERM OF PSEUDOCHARGE (ENERGY AND FORCE)              ==
!     ==========================================================================
      EPAIR0=0.D0
      STRESS=0.D0
      DO IAT1=THISTASK,NAT,NTASKS
        ISP1=ISPECIES(IAT1)
        DO IAT2=IAT1,NAT
          ISP2=ISPECIES(IAT2)
          DX=R0(1,IAT2)-R0(1,IAT1)
          DY=R0(2,IAT2)-R0(2,IAT1)
          DZ=R0(3,IAT2)-R0(3,IAT1)
          CALL BOXSPH(RBAS,-DX,-DY,-DZ,RMAX,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
          DO I1=MIN1,MAX1
            DO I2=MIN2,MAX2
              DO I3=MIN3,MAX3
                IF(IAT1.EQ.IAT2) THEN
                  IF(I1.GT.0) THEN
                    GOTO 310
                  ELSE IF(I1.EQ.0) THEN
                    IF(I2.GT.0) THEN
                      GOTO 310
                    ELSE IF(I2.EQ.0) THEN
                      IF(I3.GE.0) GOTO 310
                    END IF
                  END IF
                END IF
                T1=DBLE(I1)
                T2=DBLE(I2)
                T3=DBLE(I3)
                DR(1)=DX+RBAS(1,1)*T1+RBAS(1,2)*T2+RBAS(1,3)*T3
                DR(2)=DY+RBAS(2,1)*T1+RBAS(2,2)*T2+RBAS(2,3)*T3
                DR(3)=DZ+RBAS(3,1)*T1+RBAS(3,2)*T2+RBAS(3,3)*T3
                DIS2 = DR(1)**2 + DR(2)**2 + DR(3)**2
                DIS=SQRT(DIS2)
                IF(DIS.GT.RMAX) GOTO 310
                IF(DIS.LT.1.D-6) GOTO 310
!               ================================================================
!               ==  RADIAL PART OF PAIR INTERACTIONS                          ==
!               ================================================================
                CALL PAIRP_GETH(DIS,RCSM(ISP1),RCSM(ISP2) &
     &                         ,RCBG(ISP1),RCBG(ISP2) &
     &                         ,LRXX,A,HS,HB)
!
!               ================================================================
!               ==  ANGULAR PART OF PAIR INTERACTIONS                         ==
!               ================================================================
                CALL PAIRP_PAIRP(LMRX(ISP1),LMRX(ISP2),DR,FORCE1 &
     &                          ,QLM(1,IAT1),QLM(1,IAT2) &
     &                          ,VQLM1,VQLM2,EPAIR1,LRXX,HS,HB)
                EPAIR0=EPAIR0+EPAIR1
                DO LM=1,LMRX(ISP1)
                  VQLM(LM,IAT1)=VQLM(LM,IAT1)+VQLM1(LM)
                ENDDO
                DO LM=1,LMRX(ISP2)
                  VQLM(LM,IAT2)=VQLM(LM,IAT2)+VQLM2(LM)
                ENDDO
                DO I=1,3
                  FORCE(I,IAT1)=FORCE(I,IAT1)-FORCE1(I)
                  FORCE(I,IAT2)=FORCE(I,IAT2)+FORCE1(I)
                ENDDO
!               == THIS STRESS TENSOR IS NOT EXACTLY SYMMETRIC !
                DO I=1,3
                  DO J=1,3
                    STRESS(I,J)=STRESS(I,J)-FORCE1(I)*DR(J)
                  ENDDO
                ENDDO
 310            CONTINUE
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  ADD UP ENERGIES                                                     ==
!     ==========================================================================
      EPAIR=ESELF0+EPAIR0
!     ==========================================================================
!     == PARALLELIZE THE RESULTS                                              ==
!     == HELP ARRAY TEMPORARY STORING THE FORCES                              ==
!     == FORCE_T IS ALLOCATED NAT+1 IN ORDER TO STORE EPAIR FOR COMM          ==
!     ==========================================================================
      CALL MPE$COMBINE('MONOMER','+',EPAIR)
      CALL MPE$COMBINE('MONOMER','+',VQLM)
      CALL MPE$COMBINE('MONOMER','+',FORCE)
      CALL MPE$COMBINE('MONOMER','+',STRESS)
!
!     ==========================================================================
!     ==  HERE OPTIONAL SELF TEST                                             ==
!     ==========================================================================
      IF(TTEST_V) THEN
        CALL SELFTEST$END('PAIRP-QLM',LMRXX*NAT,VQLM,EPAIR,TBACK)
        IF(TBACK) GOTO 1000
        CALL ERROR$MSG('NORMAL STOP AFTER SELFTEST')
        CALL ERROR$STOP('PAIRP')
      END IF
      IF(TTEST_FORCE) THEN
        ! NOTE: FORCE = -DE/DR!
        CALL SELFTEST$END('PAIRP-FORCE',3*NAT,FORCE,-EPAIR,TBACK)
        IF(TBACK) GOTO 1000
        CALL ERROR$MSG('NORMAL STOP AFTER SELFTEST')
        CALL ERROR$STOP('PAIRP')
      END IF
!
!     ==========================================================================
!     ==  CLOSE DOWN                                                          ==
!     ==========================================================================
      DEALLOCATE(A)
      DEALLOCATE(HS)
      DEALLOCATE(HB)
                          CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAIRP_GETA(LRXX,A)
!     **************************************************************************
!     **                                                                      **
!     **   H(L=0,N)=1/R * D/DR**2N * R * H(L=0,N=0)                           **
!     **   H(L,N)=SUM(J): R**L * A(J,L+1,N) * (-2*R)**(2J-2)*EXP(-R**2)       **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      INTEGER(4),INTENT(IN) :: LRXX
      REAL(8)   ,INTENT(OUT):: A(LRXX+1,2*LRXX+2,LRXX)
      INTEGER(4)            :: LUP,NUP,JUP
      INTEGER(4)            :: LR3,LR2,LR1,IL,N,J
      INTEGER(4)            :: ISVAR
!     **************************************************************************
      IF(LRXX.EQ.0) RETURN
      LUP=(2*LRXX+1)
      NUP=LRXX
      JUP=LRXX
!     == INITIALIZE A TO ZERO ==================================================
      A(:,:,:)=0.D0
!
      A(1,1,1)=-2.D0
      DO IL=2,LUP+1
        A(1,IL,1)=-2.D0*A(1,IL-1,1)
      ENDDO
!
      DO N=1,NUP-1
        A(1,1,N+1)=-6.D0*A(1,1,N)+24.D0*A(2,1,N)
        DO J=2,N+1
          A(J,1,N+1)=A(J-1,1,N)-REAL(8*J-2,KIND=8)     *A(J,1,N) &
     &                         +REAL(J*(16*J+8),KIND=8)*A(J+1,1,N)
        ENDDO
        DO IL=2,LUP+1
          DO J=1,N+1
            A(J,IL,N+1)=-2.D0*A(J,IL-1,N+1)+REAL(8*J,KIND=8)*A(J+1,IL-1,N+1)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  ROUND A TO NEAREST INTEGER                                          ==
!     ==========================================================================
      DO LR3=1,LRXX
        DO LR2=1,2*LRXX+2
          DO LR1=1,LRXX+1
!           __ ROUND TO NEXT NEAREST INTEGER____________________________________
            A(LR1,LR2,LR3)=REAL(NINT(A(LR1,LR2,LR3),KIND=8),KIND=8)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  PRINT OUT FOR TESTING                                               ==
!     ==========================================================================
      IF(TPR) THEN
        DO IL=1,LUP+1
          DO N=1,NUP
            WRITE(*,FMT='(2I5,3E20.12)')IL-1,N,(A(J,IL,N),J=1,JUP)
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAIRP_ESELF(RCSM,RCBG,LMRX,QLM,VQLM,ESELF,A,LRXX)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: LMRX
      INTEGER(4),INTENT(IN)   :: LRXX
      REAL(8)   ,INTENT(IN)   :: RCSM
      REAL(8)   ,INTENT(IN)   :: RCBG
      REAL(8)   ,INTENT(IN)   :: QLM(LMRX)
      REAL(8)   ,INTENT(IN)   :: A(LRXX+1,2*LRXX+2,LRXX)
      REAL(8)   ,INTENT(OUT)  :: VQLM(LMRX)
      REAL(8)   ,INTENT(OUT)  :: ESELF
      REAL(8)   ,PARAMETER    :: PI=4.D0*ATAN(1.D0)
      REAL(8)                 :: ROOT2
      REAL(8)                 :: RC12S,RC12B
      REAL(8)                 :: ES,EB
      REAL(8)                 :: E00
      INTEGER(4)              :: LM,L,IM
      INTEGER(4)              :: LX
!     **************************************************************************
      ESELF=0.D0
      VQLM(:)=0.D0
      ROOT2=SQRT(2.D0)
      RC12S=ROOT2*RCSM
      RC12B=ROOT2*RCBG
      ES=4.D0*SQRT(PI)/RC12S
      EB=4.D0*SQRT(PI)/RC12B
      E00=ES-EB
      VQLM(1)=2.D0*E00*QLM(1)
      ESELF = E00*QLM(1)**2
      LX=INT(SQRT(REAL(LMRX))-1.D0)
      IF((LX+1)**2.NE.LMRX) THEN   !LMX=(LX+1)**2
        CALL ERROR$MSG('LMX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMRX',LMRX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('PAIRP_ESELF')
      END IF    
      LM=1
      DO L=1,LX
        ES =-ES / (DBLE(2*L+1)*RC12S)**2
        EB =-EB / (DBLE(2*L+1)*RC12B)**2
        E00=A(1,1,L)*(ES-EB)
        DO IM=1,2*L+1
          LM=LM+1
          VQLM(LM)=2.D0*E00*QLM(LM)
          ESELF = ESELF + E00*QLM(LM)**2
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAIRP_GETH(DIS,RCSM1,RCSM2,RCBG1,RCBG2 &
     &                ,LRXX,A,HS,HB)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LRXX
      REAL(8)   ,INTENT(IN) :: DIS
      REAL(8)   ,INTENT(IN) :: RCSM1
      REAL(8)   ,INTENT(IN) :: RCSM2
      REAL(8)   ,INTENT(IN) :: RCBG1
      REAL(8)   ,INTENT(IN) :: RCBG2
      REAL(8)   ,INTENT(OUT):: HS(0:2*LRXX+1,0:LRXX)
      REAL(8)   ,INTENT(OUT):: HB(0:2*LRXX+1,0:LRXX)
      REAL(8)   ,INTENT(IN) :: A(LRXX+1,2*LRXX+2,LRXX)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)            :: LUP,NUP,JUP
      INTEGER(4)            :: L,N,J
      REAL(8)               :: RC12S,RC12B
      REAL(8)               :: XS,XB
      REAL(8)               :: GOFXS,GOFXB
      REAL(8)               :: SVAR
      REAL(8)               :: QLS,QLB
      REAL(8)               :: DQS,DQB
      REAL(8)               :: FACS,FACB
      REAL(8)               :: XSL,XBL
      REAL(8)               :: SUMS,SUMB
      REAL(8)               :: SVARS,SVARB
!     **************************************************************************
      LUP=(2*LRXX+1)
      NUP=LRXX
      JUP=LRXX
      RC12S=SQRT(RCSM1**2+RCSM2**2)
      RC12B=SQRT(RCBG1**2+RCBG2**2)
      XS    = DIS/RC12S
      XB    = DIS/RC12B
      GOFXS = EXP(-XS**2)
      GOFXB = EXP(-XB**2)
      SVAR=SQRT(PI)/2.D0
      CALL LIB$ERFR8(XS,QLS)
      QLS=SVAR*QLS
      CALL LIB$ERFR8(XB,QLB)
      QLB=SVAR*QLB
      DQS=-XS*GOFXS
      DQB=-XB*GOFXB
      FACS=1.D0/XS
      FACB=1.D0/XB
      HS(0,0)=FACS*QLS/RC12S
      HB(0,0)=FACB*QLB/RC12B
      DO L=1,LUP
        QLS=QLS+DQS
        QLB=QLB+DQB
        FACS=-DBLE(2*L-1)/XS*FACS
        FACB=-DBLE(2*L-1)/XB*FACB
        HS(L,0)=FACS*QLS/RC12S**(L+1)
        HB(L,0)=FACB*QLB/RC12B**(L+1)
        DQS=2.D0*XS**2/DBLE(2*L+1) * DQS
        DQB=2.D0*XB**2/DBLE(2*L+1) * DQB
      ENDDO
!     ==   H(L,N)=SUM(J): R**L * A(J,L+1,N) * (-2*R)**(2J-2)*EXP(-R**2)=
      SVARS=(-2.D0*XS)**2
      SVARB=(-2.D0*XB)**2
      DO L=0,LUP
        XSL=XS**L
        XBL=XB**L
        DO N=1,NUP
          FACS=GOFXS
          FACB=GOFXB
          SUMS=0.D0
          SUMB=0.D0
          DO J=1,JUP
            SUMS=SUMS+A(J,L+1,N)*FACS
            SUMB=SUMB+A(J,L+1,N)*FACB
            FACS=FACS*SVARS
            FACB=FACB*SVARB
          ENDDO
          HS(L,N)=SUMS*XSL/RC12S**(L+2*N+1)
          HB(L,N)=SUMB*XBL/RC12B**(L+2*N+1)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAIRP_PAIRP(LMRX1,LMRX2,DR,FORCE &
     &              ,QLM1,QLM2,VQLM1,VQLM2,EPAIR,LRXX,HS,HB)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: LMRX1
      INTEGER(4),INTENT(IN)   :: LMRX2
      REAL(8)   ,INTENT(OUT)  :: EPAIR
      REAL(8)   ,INTENT(IN)   :: DR(3)
      REAL(8)   ,INTENT(OUT)  :: FORCE(3)
      REAL(8)   ,INTENT(IN)   :: QLM1(LMRX1)
      REAL(8)   ,INTENT(IN)   :: QLM2(LMRX2)
      REAL(8)   ,INTENT(OUT)  :: VQLM1(LMRX1)
      REAL(8)   ,INTENT(OUT)  :: VQLM2(LMRX2)
      INTEGER(4),INTENT(IN)   :: LRXX
      REAL(8)   ,INTENT(IN)   :: HS(0:2*LRXX+1,0:LRXX)
      REAL(8)   ,INTENT(IN)   :: HB(0:2*LRXX+1,0:LRXX)
      REAL(8)   ,PARAMETER    :: PI=4.D0*ATAN(1.D0)
      REAL(8)                 :: YLM((2*LRXX+2)**2)
      REAL(8)                 :: SQ4PB3
      INTEGER(4)              :: LPPXX
      INTEGER(4)              :: LMPPXX
      INTEGER(4)              :: I,LM,LM1,LM2,LM3,LM4,L1,L2,L3,L4
      INTEGER(4)              :: M1,M2,M3,M4,N3,N4
      INTEGER(4)              :: LRX1,LRX2
      INTEGER(4)              :: L3BO,L3UP,L4MIN
      REAL(8)                 :: DFAC1,DFAC2
      REAL(8)                 :: FAC,SVAR
      REAL(8)                 :: QLMT1,QLMT2
      REAL(8)                 :: FOFR
      REAL(8)                 :: CG,CG123   ! CLEBSCH GORDAN COEFFICIENT
!     **************************************************************************
      SQ4PB3=SQRT(4.D0*PI/3.D0)
      LPPXX=2*LRXX+1
      LMPPXX=(LPPXX+1)**2
      CALL GETYLM(LMPPXX,DR,YLM)
!
!     == INITIALIZE ARRAYS =====================================================
      EPAIR=0.D0
      DO I=1,3
        FORCE(I)=0.D0
      ENDDO
      DO LM=1,LMRX1
        VQLM1(LM)=0.D0
      ENDDO
      DO LM=1,LMRX2
        VQLM2(LM)=0.D0
      ENDDO
!
      LRX1=INT(SQRT(REAL(LMRX1-1,KIND=8)+1.D-5))
      LRX2=INT(SQRT(REAL(LMRX2-1,KIND=8)+1.D-5))
      LM1=0
      DFAC1=1.D0
      DO L1=0,LRX1
        DFAC1=DFAC1*REAL(2*L1+1,KIND=8)
        DO M1=-L1,L1
          LM1=LM1+1
          LM2=0
          DFAC2=1.D0
          DO L2=0,LRX2
            DFAC2=DFAC2*DBLE(2*L2+1)
            DO M2=-L2,L2
!
              LM2=LM2+1
              FAC=32.D0*PI**1.5D0*(-1.D0)**L1/(DFAC1*DFAC2)
              QLMT1=QLM1(LM1)
              QLMT2=QLM2(LM2)
              L3BO=IABS(L1-L2)
              L3UP=L1+L2
              DO L3=L3BO,L3UP,2
                N3=(L1+L2-L3)/2
                FOFR=FAC*(HS(L3,N3)-HB(L3,N3))
                LM3=L3**2
                DO M3=-L3,L3
                  LM3=LM3+1
                  CALL CLEBSCH(LM1,LM2,LM3,CG123)
                  IF(ABS(CG123).GT.1.D-5) THEN
                    SVAR=CG123*FOFR*YLM(LM3)
                    EPAIR=EPAIR+SVAR*QLMT1*QLMT2
                    VQLM1(LM1)=VQLM1(LM1)+SVAR*QLMT2
                    VQLM2(LM2)=VQLM2(LM2)+SVAR*QLMT1
!
!                   ==  FORCES =================================================
                    L4MIN=MAX0(0,L3-1)
                    DO L4=L3+1,L4MIN,-2
                      N4=(L1+L2+1-L4)/2
                      FOFR=FAC*(HS(L4,N4)-HB(L4,N4))*QLMT1*QLMT2
                      LM4=L4**2
                      DO M4=-L4,L4
                        LM4=LM4+1
                        SVAR=SQ4PB3*FOFR*CG123*YLM(LM4)
                        CALL CLEBSCH(2,LM3,LM4,CG)
                        FORCE(1)=FORCE(1)-SVAR*CG
                        CALL CLEBSCH(4,LM3,LM4,CG)
                        FORCE(2)=FORCE(2)-SVAR*CG
                        CALL CLEBSCH(3,LM3,LM4,CG)
                        FORCE(3)=FORCE(3)-SVAR*CG
                      ENDDO
                    ENDDO
                  END IF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN 
      END
