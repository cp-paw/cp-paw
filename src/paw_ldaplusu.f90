      MODULE LDAPLUSU_MODULE
      LOGICAL(4)   :: TDO=.FALSE.
      END MODULE LDAPLUSU_MODULE
!
!     ..................................................................
      SUBROUTINE LDAPLUSU$SET
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      TDO=.FALSE.
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LDAPLUSU$SETTING(TDO_)
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(OUT) :: TDO_
      TDO_=TDO
      RETURN
      END
!
!     ...................................................AUTOPI.........
      SUBROUTINE LDAPLUSU(NRX,LMNXX,NSPIN,LOX,LNX &
     &                ,gid,NR,AEZ,AEPHI,DENMAT,DETOT,DH)
!     **                                                              **
!     **                                                              **
      USE LDAPLUSU_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NRX
      INTEGER(4),INTENT(IN) :: LMNXX
      INTEGER(4),INTENT(IN) :: NSPIN
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: gid   ! grid id
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: AEPHI(NRX,LNX)
      REAL(8)   ,INTENT(IN) :: DENMAT(LMNXX,LMNXX,NSPIN)
      REAL(8)   ,INTENT(OUT):: DH(LMNXX,LMNXX,NSPIN)
      REAL(8)   ,INTENT(OUT):: DETOT
      REAL(8)   ,ALLOCATABLE:: OCC(:,:)
      REAL(8)   ,ALLOCATABLE:: POT(:,:)
      REAL(8)               :: X(LNX,LNX)
      REAL(8)               :: DH1(LMNXX,LMNXX,NSPIN)
      REAL(8)               :: RJ(0:3)
      REAL(8)               :: RU(0:3)
      REAL(8)               :: WORK(NRX)
      REAL(8)               :: XEXP,RI
      REAL(8)               :: EV
      REAL(8)               :: RASA
      INTEGER(4)            :: LX,LMNX,LN,L,LMX,IZ,ISPIN,LM,M
      INTEGER(4)            :: LN1,LN2,L1,L2,M1,M2,LMN1,LMN2,LM1,LM2
      INTEGER(4)            :: IR
      REAL(8)               :: OCC1UP,OCC2UP,OCC1DOWN,OCC2DOWN
      REAL(8)               :: V1UP,V1DOWN,V2UP,V2DOWN
      REAL(8)               :: AVPOT,AVOCC
      REAL(8)               :: SVAR
      REAL(8)               :: r(nr)
!     *******************************************************************
!     == POINTER ARRAYS
      IF(.NOT.TDO) RETURN
      call radial$r(gid,nr,r)
      DETOT=0.D0
!     == EVALUATE LMX ===================================================
      LX=-1
      LMNX=0
      DO LN=1,LNX
        L=LOX(LN)
        LX=MAX(LX,L)
        LMNX=LMNX+2*L+1
      ENDDO
      LMX=(LX+1)**2    
!     == ALLOCATE ARRAYS ===============================================
      ALLOCATE(OCC(LMX,NSPIN))
      ALLOCATE(POT(LMX,NSPIN))
!
!     ==================================================================      
!     ==  LOOKUP U AND J                                              ==
!     ==================================================================      
      CALL CONSTANTS('EV',EV)
      IZ=NINT(AEZ)
      CALL PERIODICTABLE$GET(IZ,'R(ASA)',RASA)
      DO L=0,3
        RJ(L)=0.D0
        RU(L)=0.D0
      ENDDO
!
      IF(IZ.EQ.29) THEN
        RJ(2)=0.98D0*EV
        RU(2)=7.5D0*EV
        RU(2)=10.0D0*EV
      ELSE
        RETURN
      END IF
!
!
!     ==================================================================      
!     ==  OVERLAP                                                     ==
!     ==================================================================      
      X(:,:)=0.D0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO LN2=LN1,LNX
          L2=LOX(LN2)
          X(LN1,LN2)=0.D0
          X(LN2,LN1)=0.D0
          IF(L1.EQ.L2) THEN 
            DO IR=1,NR
              WORK(ir)=AEPHI(ir,LN1)*AEPHI(ir,LN2)*R(ir)**2
              IF(R(ir).GT.RASA) WORK(IR)=0.D0
            ENDDO
            CALL RADIAL$INTEGRAL(gid,NR,WORK,X(LN1,LN2))
            X(LN2,LN1)=X(LN1,LN2)
          END IF
        ENDDO
      ENDDO
!     PRINT*,'X'
!     DO LN1=1,LNX
!       WRITE(*,FMT='(9F10.5)')(X(LN1,LN2),LN2=1,LNX)
!     ENDDO
!
!     ==================================================================      
!     ==  CALCULATE OCCUPATIONS                                       ==
!     ==================================================================      
      DO ISPIN=1,NSPIN
        DO LM=1,LMX
          OCC(LM,ISPIN)=0.D0
        ENDDO
      ENDDO   
!
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO M1=1,2*L1+1
          LMN1=LMN1+1
          LMN2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            DO M2=1,2*L2+1
              LMN2=LMN2+1
              IF(L1.EQ.L2.AND.M1.EQ.M2) THEN
                LM=L1**2+M1
                DO ISPIN=1,NSPIN
                  OCC(LM,ISPIN)=OCC(LM,ISPIN) &
     &                         +DENMAT(LMN1,LMN2,ISPIN)*X(LN1,LN2)
                ENDDO
!               PRINT*,'LM,LMN1,LMN2,LN1,LN2',LM,LMN1,LMN2,LN1,LN2
!    &                ,DENMAT(LMN1,LMN2,1),DENMAT(LMN1,LMN2,2)
              END IF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      PRINT*,'OCC ',NSPIN,LMX
      DO ISPIN=1,NSPIN
        WRITE(*,FMT='(9F8.5)')(OCC(LM,ISPIN),LM=1,LMX)
      ENDDO
!
!     ==================================================================      
!     ==  SUBTRACT AVERAGE OCCUPATION FROM TOTAL OCCUPATION           ==
!     ==================================================================      
      DO L=0,LX
        AVOCC=0.D0
        DO M=1,2*L+1
          LM=L**2+M
          AVOCC=AVOCC+OCC(LM,1)
        ENDDO
        AVOCC=AVOCC/DBLE(1*(2*L+1))
!       PRINT*,'AVOCC ',L,AVOCC
        DO M=1,2*L+1
          LM=L**2+M
          OCC(LM,1)=OCC(LM,1)-AVOCC
        ENDDO
      ENDDO
!     PRINT*,'OCC - AVOCC'
!     DO ISPIN=1,NSPIN
!       WRITE(*,FMT='(9F8.5)')(OCC(LM,ISPIN),LM=1,LMX)
!     ENDDO
!
!     ==================================================================      
!     ==  CALCULATE TOTAL ENERGY CONTRIBUTION                         ==
!     ==================================================================      
      DETOT=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMX
          POT(LM,ISPIN)=0.D0
        ENDDO
      ENDDO
      DO L=0,LX
        DO M1=1,2*L+1
          LM1=L**2+M1
          DO M2=1,2*L+1
            LM2=L**2+M2
            OCC1UP  =0.5D0*OCC(LM1,1)
            OCC2UP  =0.5D0*OCC(LM2,1)
            OCC1DOWN=0.5D0*OCC(LM1,1)
            OCC2DOWN=0.5D0*OCC(LM2,1)
            IF(NSPIN.EQ.2) THEN
              OCC1UP  =OCC1UP  +0.5D0*OCC(LM1,2)
              OCC1DOWN=OCC1DOWN-0.5D0*OCC(LM1,2)
              OCC2UP  =OCC2UP  +0.5D0*OCC(LM2,2)
              OCC2DOWN=OCC2DOWN-0.5D0*OCC(LM2,2)
            END IF
            SVAR=0.5D0*RU(L)
            DETOT=DETOT+SVAR*(OCC1UP*OCC2DOWN+OCC1DOWN*OCC2UP)
            V1UP  =SVAR*OCC2DOWN
            V1DOWN=SVAR*OCC2UP
            V2UP  =SVAR*OCC1DOWN
            V2DOWN=SVAR*OCC1UP
            IF(M1.NE.M2) THEN
              SVAR=0.5D0*(RU(L)-RJ(L))
              DETOT =DETOT +SVAR*(OCC1UP*OCC2UP+OCC1DOWN*OCC2DOWN)
              V1UP  =V1UP  +SVAR*OCC2UP
              V1DOWN=V1DOWN+SVAR*OCC2DOWN
              V2UP  =V2UP  +SVAR*OCC1UP
              V2DOWN=V2DOWN+SVAR*OCC1DOWN
            END IF
            POT(LM1,1)=POT(LM1,1)+0.5D0*(V1UP+V1DOWN)
            POT(LM2,1)=POT(LM2,1)+0.5D0*(V2UP+V2DOWN)
            IF(NSPIN.EQ.2) THEN
              POT(LM1,2)=POT(LM1,2)+0.5D0*(V1UP-V1DOWN)
              POT(LM2,2)=POT(LM2,2)+0.5D0*(V2UP-V2DOWN)
            END IF
          ENDDO
        ENDDO
      ENDDO       
!     PRINT*,'POT-AVPOT'
!     DO ISPIN=1,NSPIN
!       WRITE(*,FMT='(9F8.5)')(POT(LM,ISPIN),LM=1,LMX)
!     ENDDO
!
!     ==================================================================      
!     ==  SUBTRACT AVERAGE POTENTIAL                                  ==
!     ==================================================================      
      DO L=0,LX
        AVPOT=0.D0
        DO M=1,2*L+1
          LM=L**2+M
          AVPOT=AVPOT+POT(LM,1)
        ENDDO
        AVPOT=AVPOT/DBLE(1*(2*L+1))
        DO M=1,2*L+1
          LM=L**2+M
          POT(LM,1)=POT(LM,1)-AVPOT
        ENDDO
      ENDDO
      PRINT*,'POT IN EV'
      DO ISPIN=1,NSPIN
        WRITE(*,FMT='(9F8.5)')(POT(LM,ISPIN)/EV,LM=1,LMX)
      ENDDO
!
!     ==================================================================      
!     ==  CALCULATE DH                                                ==
!     ==================================================================      
      DO ISPIN=1,NSPIN
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            DH1(LMN1,LMN2,ISPIN)=0.D0
          ENDDO
        ENDDO
      ENDDO   
!
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO M1=1,2*L1+1
          LMN1=LMN1+1
          LMN2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            DO M2=1,2*L2+1
              LMN2=LMN2+1
              IF(L1.EQ.L2.AND.M1.EQ.M2) THEN
                LM=L1**2+M1
                DO ISPIN=1,NSPIN
                  DH1(LMN1,LMN2,ISPIN)=DH1(LMN1,LMN2,ISPIN) &
     &                             +X(LN1,LN2)*POT(LM,ISPIN)
                  DH(LMN1,LMN2,ISPIN)=DH(LMN1,LMN2,ISPIN) &
     &                             +X(LN1,LN2)*POT(LM,ISPIN)
                ENDDO
              END IF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(POT)
      DEALLOCATE(OCC)
!
      CALL CONSTANTS('EV',EV)
      PRINT*,'LDA+U ENERGY',DETOT,'H;',DETOT/EV,' EV'
      DO ISPIN=1,NSPIN
!       PRINT*,'DENMAT  FOR SPIN ',ISPIN
        DO LMN1=1,LMNX
!          WRITE(*,FMT='(9F10.5)')
!    &         (DENMAT(LMN1,LMN2,ISPIN)/EV,LMN2=1,LMNX)
        ENDDO
      ENDDO
      DO ISPIN=1,NSPIN
!       PRINT*,'DH1 IN EV FOR SPIN ',ISPIN
        DO LMN1=1,LMNX
!         WRITE(*,FMT='(9F10.5)')(DH1(LMN1,LMN2,ISPIN)/EV,LMN2=1,LMNX)
        ENDDO
      ENDDO
      RETURN
      END
