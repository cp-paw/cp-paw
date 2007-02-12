!TODO :
! DATH IS STILL REAL AND SHOULD PROBABLY BE COMPLEX LIKE DENMAT

MODULE LDAPLUSU_MODULE
TYPE THISTYPE
LOGICAL(4)             :: TINI=.FALSE.
LOGICAL(4)             :: TON=.FALSE.
INTEGER(4)             :: GID
INTEGER(4)             :: NR
INTEGER(4)             :: LNXCHI
INTEGER(4),POINTER     :: LOXCHI(:)
INTEGER(4),POINTER     :: NORB(:)  ! X(# LOCAL FUNCTIONS PER ANGULAR MOMENTUM)
INTEGER(4)             :: LRX
INTEGER(4)             :: NCHI
REAL(8)   ,POINTER     :: CHI(:,:)   ! CORRELATED (HEAD) ORBITALS
LOGICAL(4)             :: USEDIEL=.FALSE.
LOGICAL(4)             :: USEUPAR=.FALSE.
LOGICAL(4)             :: USEJPAR=.FALSE.
REAL(8)                :: DIEL=1     ! DIELECTRIC CONSTANT
REAL(8)                :: UPAR=0.D0  ! U-PARAMETER
REAL(8)                :: JPAR=0.D0  ! J-PARAMETER
INTEGER(4)             :: MAINLN(2)=(/0,0/)  ! DEFINES SHELL TO WHICH UPAR AND JPAR REFER TO
REAL(8)                :: RCUT=0.D0  ! RADIUS WHERE NODE IN CORRELATED WAVE FUNCTION IS ENFORCED
REAL(8)   ,POINTER     :: ULITTLE(:,:,:,:,:) ! SLATER INTEGRALS (EXCEPT FACTOR)
REAL(8)   ,POINTER     :: DOWNFOLD(:,:)      ! MAPS PARTIAL WAVES TO CORRELATED ORBITALS
END TYPE THISTYPE
TYPE(THISTYPE),ALLOCATABLE,TARGET :: THISARRAY(:)
TYPE(THISTYPE),POINTER :: THIS
LOGICAL(4)             :: TON=.FALSE.
INTEGER(4)             :: NSP=0    ! #(ATOM TYPES)
INTEGER(4)             :: ISP=0    ! SELECTED ATOM TYPE
CHARACTER(8)           :: DCTYPE='FLL' ! SPECIFIES TYPE OF DOUBLE COUNTING CORRECTION
END MODULE LDAPLUSU_MODULE
!
!       .............................................................................
        SUBROUTINE LDAPLUSU(IAT,GID,NR,LNX,LMNXX,NSPIN,LOX &
       &                   ,AEZ,AEPHI,DENMAT_,DETOT,DATH)
        USE LDAPLUSU_MODULE
        USE PERIODICTABLE_MODULE
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: IAT
        INTEGER(4),INTENT(IN) :: GID
        INTEGER(4),INTENT(IN) :: NR
        INTEGER(4),INTENT(IN) :: LNX
        INTEGER(4),INTENT(IN) :: LMNXX
        INTEGER(4),INTENT(IN) :: NSPIN
        INTEGER(4),INTENT(IN) :: LOX(LNX)
        REAL(8)   ,INTENT(IN) :: AEZ
        REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX)
        COMPLEX(8),INTENT(IN) :: DENMAT_(LMNXX,LMNXX,NSPIN)
        REAL(8)   ,INTENT(OUT):: DETOT
        REAL(8)   ,INTENT(OUT):: DATH(LMNXX,LMNXX,NSPIN)
        REAL(8)               :: DENMAT(LMNXX,LMNXX,NSPIN)
        REAL(8)               :: R(NR)
        REAL(8)               :: DIEL        ! DIELECTRIC CONSTANT DIEL
        LOGICAL               :: MYUIJKL     ! CHOOSE CALCULATION OF UIJKL
        LOGICAL               :: DEFJP       ! CHOOSE DEFINITION OF JPARAMETER
        INTEGER(4)            :: CHOOSEPOT   ! CHOOSE DC POTENTIAL (1:FLL; 2:AMF)
        INTEGER(4)            :: CHOOSENUE   ! CHOOSE DEFINITION OF NUE (1:O.BENG.;2:PETER)
        REAL(8)   ,ALLOCATABLE:: UIJKL(:,:,:,:)
        INTEGER(4)            :: LN1,LMN1
        INTEGER(4)            :: L
        INTEGER(4)            :: LMSSH       !START VALUE OF S. HARMONICS
        INTEGER(4)            :: DIM
        INTEGER(4)            :: M1,M2,M3
        REAL(8)               :: UPARAMETER
        REAL(8)               :: JPARAMETER
        REAL(8)               :: JPARAMETERST
        REAL(8),ALLOCATABLE   :: NUE(:,:,:)
        REAL(8)               :: NUENN(LMNXX,LMNXX)
        REAL(8)               :: ELDAUSUMMAND
        REAL(8)               :: POTUPDOWN(LMNXX,LMNXX,2)
        REAL(8),ALLOCATABLE   :: POTAMF(:,:,:)
        REAL(8),ALLOCATABLE   :: POTFLL(:,:,:)
        REAL(8),ALLOCATABLE   :: OVER(:,:) 
        INTEGER(4)            :: NOXA        !NUMBER OF XA
        REAL(8)               :: RCUT        !RADIUS IN WHICH STRONG CORRELATION EFFECTS WERE EXPECTED.SHOULD 
                                             !BE EQUAL WITH THE AUGMENTATION RADIUS
        INTEGER(4)            :: WAN         !FOR WHICH ATOM LDA+U SHOULD BE USED
        INTEGER(4)            :: HML         !HOW MANY L NUMBERS SHOULD BE DESCRIBED WITH LDA+U? (NOT REALLY TESTED)
        INTEGER(4)            :: WHICHL(2)   !WHICH L NUMBERS HAVE THE LOCALIZED ORBITALS    
        INTEGER(4)            :: NOLDAU  !
        INTEGER(4)            :: SPIN
        REAL(8),ALLOCATABLE   :: XALPHA(:,:)
        LOGICAL,PARAMETER     :: ATOMICORBIT=.TRUE.  !SHOULD THE ATOMIC ORBITALS BE USED?
!       *****************************************************************************
CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
CALL LDAPLUSU$ETOT(ISP,LMNXX,NSPIN,DENMAT_,DETOT,DATH)
RETURN
        DETOT=0.D0        
        DATH=0.D0
        IF(.NOT.TON) RETURN
        IF(NSPIN.NE.2) THEN
           PRINT*, 'FOR AN LDA+U CALCULATION NSPIN=2 SHOULD BE USED'
           STOP
        END IF
        DENMAT(:,:,:)=REAL(DENMAT_(:,:,:))
!       CALL SETUP$AEZ(ISP,AEZ)
!       CALL SETUP$ISELECT(ISP)
!       CALL SETUP$GETI4('GID',GID)
        CALL RADIAL$R(GID,NR,R)
!
!       ===============================================================================
!       ==  READ PARAMETERS                                                          ==
!       ===============================================================================
PRINT*, 'READ LDA+U PARAMETERS'
        CALL READPARAMETERS(UPARAMETER,JPARAMETER,DIEL,MYUIJKL,DEFJP,CHOOSEPOT,CHOOSENUE,RCUT,HML,WAN,WHICHL) 
!
PRINT*, 'THE PARAMETERS HAVE THE FOLLOWING VALUES:'
WRITE(*,'(A,F10.3)') 'HUBBARD U: ',UPARAMETER
WRITE(*,'(A,F10.3)') 'EXCHANGE J: ',JPARAMETER
WRITE(*,'(A,2I5)') 'ANGULAR QUANTUM NUMBER: ', WHICHL
PRINT*, 'AEZ,LOX(LNX)', AEZ, LOX(:)
!
!       ===============================================================================
!       ==  RETURN IF WRONG ATOMIC NUMBER                                            ==
!       ===============================================================================
        IF(NINT(AEZ).NE.WAN) RETURN
!
!       ===============================================================================
!       ==                                                                           ==
!       ===============================================================================
        DO NOLDAU=1,HML
          POTUPDOWN(:,:,:)=0.D0
          L=WHICHL(NOLDAU)
          DIM=2*L+1    !DIMENSION OF UIJKL,NUE AND XALPHA 
          LMSSH=L**2   !START VALUE OF S.HARMONICS
!
!         ===========================================================================
!         ==  LOAD LOCALIZE ORBITAL X_ALPHA                                        ==
!         ===========================================================================
          ALLOCATE(XALPHA(NR,DIM))
          CALL  LOCALIZEDORBITAL(GID,NR,LNX,L,LOX,ATOMICORBIT,AEPHI,DIM,RCUT,XALPHA) 
DO SPIN=1,2
  PRINT*, 'DENMAT FOR NSPIN=',SPIN ,'AND ATOMIC-NUMBER=',AEZ  
  PRINT*, '======================================='
  LMN1=0
  DO LN1=1,LNX
    DO M1=1,2*LOX(LN1)+1
      LMN1=LMN1+1
      WRITE(*,'(I2,"   ",100(F8.3))') LOX(LN1),REAL(DENMAT(LMN1,:,SPIN))
    END DO
  END DO
  PRINT*, '======================================='
END DO
!
!         ===========================================================================
!         ==  LOAD ONE-CENTER DENSITY NUE                                          ==
!         ===========================================================================
          NOXA=DIM
          ALLOCATE(NUE(DIM,DIM,2))
          NUE(:,:,:)=0.D0
          ALLOCATE(OVER(NOXA,LMNXX))
          OVER(:,:)=0.D0
          CALL NUETENSOR(GID,NR,NOXA,OVER,DIM,DENMAT,LNX,LOX &
         &              ,LMNXX,NSPIN,AEPHI,NUE,NUENN,L,CHOOSENUE,RCUT,XALPHA)
DO SPIN=1,2
  PRINT*, 'NUETENSOR FOR NSPIN=',SPIN ,'AND ATOMIC-NUMBER=',AEZ  
  PRINT*, '======================================='
  DO M1=1,DIM
    WRITE(*,'(100(F8.3))')REAL(NUE(M1,:,SPIN))
  END DO
  PRINT*, '======================================='
END DO
!
!         ===========================================================================
!         ==  LOAD U-MATRIX                                                        ==
!         ===========================================================================
          ALLOCATE (UIJKL(DIM,DIM,DIM,DIM))
          UIJKL(:,:,:,:)=0.D0
          CALL UIJKLTENSOR(GID,NR,L,LMSSH,LNX &
         &                ,DIEL,DIM,MYUIJKL,UIJKL,UPARAMETER,JPARAMETER,XALPHA)
        
!         == CALCULATE UPARAMETER ===================================================
          UPARAMETER=0.0D0
          DO M1=1,2*L+1
            DO M3=1,2*L+1
              UPARAMETER=UPARAMETER+UIJKL(M1,M3,M1,M3)
            END DO
          END DO
        
          UPARAMETER=UPARAMETER/(REAL(2*L+1,KIND=8)**2)
          PRINT *,'UPARAMETER[EV]',UPARAMETER*27.2144D0 

!         ==  CALCULATE JPARAMETER ==========================================
          IF (.NOT.DEFJP) THEN
            JPARAMETER=0.0D0
            DO M1=1,2*L+1
              DO M3=1,2*L+1
                JPARAMETER=UIJKL(M1,M3,M1,M3)-UIJKL(M1,M3,M3,M1)+JPARAMETER
              END DO
            END DO
            JPARAMETER=UPARAMETER-1/(2*FLOAT(L)*(2*FLOAT(L)+1))*JPARAMETER
          ELSE
            JPARAMETER=0.0D0
            JPARAMETERST=0.0D0
            DO M1=1,2*L+1
              DO M3=1,2*L+1
                IF(M1.EQ.M3) THEN
                  JPARAMETERST=0.D0
                ELSE
                  JPARAMETERST=UIJKL(M1,M3,M3,M1)
                END IF
                JPARAMETER=JPARAMETERST+JPARAMETER
              ENDDO
            ENDDO
            JPARAMETER=JPARAMETER/REAL(2*L*(2*L+1),KIND=8)
            JPARAMETER=7.D0/5.D0*JPARAMETER
          END IF
        
          PRINT*,'JPARAMETER',JPARAMETER*27.2144D0 
          PRINT*,'U/J',UPARAMETER/JPARAMETER

!         ==  GET LDA-U CORRECTION FOR USED LN VALUE =========================
          ALLOCATE(POTFLL(DIM,DIM,2))
          ALLOCATE(POTAMF(DIM,DIM,2))
          POTAMF(:,:,:)=0.D0
          POTFLL(:,:,:)=0.D0
          CALL LDAUCORRECTION(LOX,LNX,POTUPDOWN,POTAMF,POTFLL,L,&
                              LMNXX,NUE,NUENN,DIM,UIJKL,UPARAMETER,&
                              JPARAMETER,NSPIN,CHOOSEPOT,ELDAUSUMMAND,&
                              DATH,OVER,NOXA,CHOOSENUE)
PRINT*,' POTUPDOWN ',POTUPDOWN 
          DETOT=DETOT+ELDAUSUMMAND
          PRINT*, '================='
          PRINT*,' TOTAL ENERGY ',DETOT,' IN EV: ',DETOT*27.21440D0
          PRINT*, '=================' 
          DEALLOCATE(UIJKL)
          DEALLOCATE(NUE)
          DEALLOCATE(POTAMF)
          DEALLOCATE(POTFLL)
          DEALLOCATE(OVER)
          DEALLOCATE(XALPHA)

!         ==  TRANSFORM INTO TOTAL AND SPIN FROM SPIN UP AND SPIN DOWN ========
          DO M1=1,LMNXX
            DO M2=1,LMNXX
              DATH(M1,M2,1)=(POTUPDOWN(M1,M2,1)+POTUPDOWN(M1,M2,2))*0.5D0
              IF (NSPIN.EQ.2) THEN
                DATH(M1,M2,2)=(POTUPDOWN(M1,M2,1)-POTUPDOWN(M1,M2,2))*0.5D0
              END IF
            ENDDO
          ENDDO
        END DO

        DO SPIN=1,2
          PRINT*, 'DATH FOR NSPIN=',SPIN ,'AND ATOMIC-NUMBER=',AEZ  
          PRINT*, '======================================='
          LMN1=0
          DO LN1=1,LNX
            DO M1=1,2*LOX(LN1)+1
              LMN1=LMN1+1
!             PRINT*, DATH(LMN1,:,SPIN) 
              WRITE(*,'(I2.1,A,100(F8.6))') LOX(LN1),'  ', DATH(LMN1,:,SPIN)
            END DO
          END DO
          PRINT*, '======================================='
        END DO
STOP 'FORCED STOP'
        RETURN
        END
!
!       .............................................................................
        SUBROUTINE NUETENSOR(GID,NR,NOXA,OVER,DIM,DENMAT,LNX,LOX,LMNXX,NSPIN,AEPHI &
       &                    ,NUE,NUENN,L,CHOOSENUE,RCUT,XALPHA)
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: GID
        INTEGER(4),INTENT(IN) :: NR
        INTEGER(4),INTENT(IN) :: DIM
        INTEGER(4),INTENT(IN) :: LNX
        INTEGER(4),INTENT(IN) :: LMNXX
        INTEGER(4),INTENT(IN) :: LOX(LNX)
        INTEGER(4),INTENT(IN) :: NSPIN
        REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX)
        REAL(8)   ,INTENT(IN) :: DENMAT(LMNXX,LMNXX,NSPIN)
        INTEGER(4),INTENT(IN) :: L
        INTEGER(4),INTENT(IN) :: CHOOSENUE
        REAL(8)   ,INTENT(OUT):: NUE(DIM,DIM,2)
        REAL(8)   ,INTENT(OUT):: NUENN(LMNXX,LMNXX)
        REAL(8)   ,INTENT(IN) :: XALPHA(NR,DIM)
        INTEGER(4),INTENT(IN) :: NOXA      !NUMBER OF XA
        REAL(8)   ,INTENT(OUT):: OVER(NOXA,LMNXX)
        REAL(8)   ,INTENT(IN) :: RCUT       !2.33D0      !20.0D0  !2.33     
        REAL(8)               :: R(NR)
        INTEGER(4)            :: ALPHA,ALPHA2
        REAL(8)               :: NU(LMNXX,LMNXX,2)
        REAL(8)   ,ALLOCATABLE:: XA(:,:)
        INTEGER(4)            :: BETA,BETA2
        INTEGER(4)            :: SHOAE(LMNXX) ! LM(LMN) FOR AEPHI
        INTEGER(4),ALLOCATABLE:: SHOXA(:)     ! LM(LMN) FOR XA        
        REAL(8)               :: WORK(NR),SVAR
        REAL(8)               :: X(LNX,LNX)
        INTEGER(4)            :: IR
        INTEGER(4)            :: L1,L2
        INTEGER(4)            :: M
        INTEGER(4)            :: LN1,LN2,LN
        INTEGER(4)            :: LM1,LM2,LM
        INTEGER(4)            :: LMN1,LMN2,LMN
        INTEGER(4)            :: SPIN
!       *****************************************************************************
        CALL RADIAL$R(GID,NR,R)
        NUE(:,:,:)=0.D0
        NUENN(:,:)=0.D0
        OVER(:,:)=0.D0
        IF(CHOOSENUE.EQ.1) THEN
!         ===========================================================================
!         ==  O. BENGONE 'S IMPLEMENTATION                                         ==
!         ==  PROJECT ON SPHERE                                                    ==
!         ===========================================================================

         DO LMN1=1,LMNXX
           DO LMN2=1,LMNXX
             NU(LMN1,LMN2,1)=0.5D0*(DENMAT(LMN1,LMN2,1)+DENMAT(LMN1,LMN2,2))
             NU(LMN1,LMN2,2)=0.5D0*(DENMAT(LMN1,LMN2,1)-DENMAT(LMN1,LMN2,2))
           END DO
         END DO
         X(:,:)=0.D0
         DO LN1=1,LNX
           L1=LOX(LN1)
           DO LN2=LN1,LNX
             L2=LOX(LN2)
             X(LN1,LN2)=0.D0
             X(LN2,LN1)=0.D0
             IF(L1.EQ.L2) THEN
               WORK(:)=AEPHI(:,LN1)*AEPHI(:,LN2)*R(:)**2
               DO IR=1,NR
                 IF(R(IR).GT.RCUT)  WORK(IR)=0.D0
               ENDDO
               CALL RADIAL$INTEGRAL(GID,NR,WORK,X(LN1,LN2))
               X(LN2,LN1)=X(LN1,LN2)
             END IF
           ENDDO
         ENDDO

         NUE(:,:,:)=0.D0
         LMN1=0
         DO LN1=1,LNX
           L1=LOX(LN1)
           DO LM1=1,2*L1+1
             LMN1=LMN1+1
             LMN2=0
             DO LN2=1,LNX
               L2=LOX(LN2)
               DO LM2=1,2*L2+1
                 LMN2=LMN2+1
                 IF ((L1 == L).AND.(L2 == L)) THEN 
                   DO SPIN=1,2  
                     NUE(LM1,LM2,SPIN)=NUE(LM1,LM2,SPIN) & 
       &                            +NU(LMN1,LMN2,SPIN)*X(LN1,LN2)  
                   ENDDO 
                 ENDIF 
               ENDDO
             ENDDO
           ENDDO
         ENDDO

         PRINT*, 'RCUT', RCUT

         LMN1=0
         DO LN1=1,LNX
           L1=LOX(LN1)
           DO LM1=1,2*L1+1
             LMN1=LMN1+1
             LMN2=0
             DO LN2=1,LNX
               L2=LOX(LN2)
               DO LM2=1,2*L2+1
                 LMN2=LMN2+1
                 NUENN(LMN1,LMN2)=X(LN1,LN2)
               END DO
             END DO
           END DO
         END DO
       ELSE IF(CHOOSENUE.EQ.2) THEN
!         ===========================================================================
!         ==  PROJECT ON LOCAL ORBITALS                                            ==
!         ===========================================================================
 
!         SPHERICAL HARMONICS OF AEPHI                  
          LMN=0
          DO LN=1,LNX
            DO M=1,2*LOX(LN)+1
              LMN=LMN+1
              SHOAE(LMN)=LOX(LN)**2+LM     ! (LM)
            END DO
          END DO

!         DEFINE XA
          ALLOCATE(XA(NR,NOXA)) 
          ALLOCATE(SHOXA(NOXA))
          DO LM=1,NOXA
            XA(:,LM)=XALPHA(:,1)     !?????
          END DO
          
!         SPHERICAL HARMONICS OF XA
          DO LM1=1,NOXA
            SHOXA(LM1)=9+LM1
          END DO
!         == TRUNCATE XA BEYOND RCUT
          DO IR=1,NR
            IF (R(IR).GT.RCUT) THEN
              XA(IR,1)=0.D0
            END IF
          END DO

!         == NORMALIZE XA
          WORK(:)=XA(:,1)*XA(:,1)*R(:)**(2)
          CALL RADIAL$INTEGRAL(GID,NR,WORK,SVAR)
          XA(:,1)=XA(:,1)/SQRT(SVAR)
          WORK(:)=XA(:,1)*XA(:,1)*R(:)**(2)
          CALL RADIAL$INTEGRAL(GID,NR,WORK,SVAR)
          PRINT*, 'RADIALINTEGRAL2',SVAR

!         ==  CALCULATE OVERLAP  <XA|AEPHI>
          DO LM1=1,NOXA
            LMN1=0
            DO LN1=1,LNX
              DO LM2=1,2*LOX(LN1)+1
                LMN1=LMN1+1
                IF(SHOAE(LMN1).EQ.SHOXA(LM1)) THEN
                  WORK(:)=AEPHI(:,LN1)*XA(:,LM1)*R(:)**2
                  CALL RADIAL$INTEGRAL(GID,NR,WORK,SVAR) 
                  OVER(LM1,LMN1)=SVAR
                ELSE 
                  OVER(LM1,LMN1)=0.D0
                END IF
              END DO
            END DO
          END DO
         
          DO LM1=1,NOXA
            LMN1=0
            DO LN1=1,LNX
              DO LM2=1,2*LOX(LN1)+1 
                LMN1=LMN1+1
                WRITE(*,'(A,I2,A,I2,A,F8.5,A,I3,A,I3,A,I3)') 'OVER   ', LM1,'  ',LMN1,'  ' &
      &           ,OVER(LM1,LMN1),'   SHOAE ',SHOAE(LMN1), '  ',SHOXA(LM1),'   ',LOX(LN1)
              END DO
            END DO
          END DO
          
!         ==  CALCULATE NUE ========================================================
          DO LMN1=1,LMNXX
            DO LMN2=1,LMNXX
              NU(LMN1,LMN2,1)=0.5D0*(DENMAT(LMN1,LMN2,1)+DENMAT(LMN1,LMN2,2))
              NU(LMN1,LMN2,2)=0.5D0*(DENMAT(LMN1,LMN2,1)-DENMAT(LMN1,LMN2,2))
            END DO
          END DO

          NUE(:,:,:)=0.D0  
          DO SPIN=1,2
            DO ALPHA=1,DIM
              DO ALPHA2=1,DIM
                DO BETA=1,17
                  DO BETA2=1,17
                    NUE(ALPHA,ALPHA2,SPIN)=NUE(ALPHA,ALPHA2,SPIN) &
       &                       +NU(BETA,BETA2,SPIN)*OVER(ALPHA,BETA)*OVER(ALPHA2,BETA2)
                  END DO
                END DO
              END DO
            END DO
          END DO

          DO SPIN=1,2  
            DO ALPHA=1,DIM
              WRITE(*,'(F6.4,A,F6.4,A,F6.4,A,F6.4,A,F6.4,A,F6.4,A,F6.4)') &
       &             NUE(ALPHA,1,SPIN),'   ',NUE(ALPHA,2,SPIN),'   ',NUE(ALPHA,3,SPIN),'   ',&
       &             NUE(ALPHA,4,SPIN),'   ',NUE(ALPHA,5,SPIN),'   ',NUE(ALPHA,6,SPIN),'   ',&
       &             NUE(ALPHA,7,SPIN)
            END DO
          END DO
          DEALLOCATE(XA)
          DEALLOCATE(SHOXA)
        END IF
        RETURN
        END 
!
!       .............................................................................
        SUBROUTINE UIJKLTENSOR(GID,NR,L,LMSSH,LNX,DIEL,DIM,&
       &                        MYUIJKL,UIJKL,UPARAMETER,JPARAMETER,XALPHA)
        IMPLICIT NONE
        INTEGER(4),INTENT(IN)   :: GID
        INTEGER(4),INTENT(IN)   :: NR
        INTEGER(4),INTENT(IN)   :: L
        INTEGER(4),INTENT(IN)   :: LMSSH
        INTEGER(4),INTENT(IN)   :: LNX
        REAL(8)   ,INTENT(IN)   :: DIEL
        INTEGER(4),INTENT(IN)   :: DIM
        LOGICAL   ,INTENT(IN)   :: MYUIJKL
        REAL(8)   ,INTENT(INOUT):: UIJKL(DIM,DIM,DIM,DIM)
        REAL(8)   ,INTENT(IN)   :: XALPHA(NR,DIM)
        REAL(8)   ,INTENT(IN)   ::UPARAMETER        !=3.0D0/27.2144D0
        REAL(8)   ,INTENT(IN)   ::JPARAMETER    !=0.D0   !0.95D0/27.2144D0  
        REAL(8)                 :: R(NR)
        INTEGER(4)              :: M1,M2,M3,M4
        INTEGER(4)              :: DM1,DM2,DM3,DM4
        INTEGER(4)              :: LM
        REAL(8)                 :: CG
        REAL(8)                 :: RHO1(NR),RHO2(NR)
        INTEGER(4)              :: LMIR
        REAL(8)                 :: V(NR)
        REAL(8)                 :: AUX(NR)
        REAL(8)                 :: SO
        ! FOR UIJKL DEFINITION  2
        REAL(8)                 ::SF(5)
        REAL(8)                 ::CON,CON2
        INTEGER(4)              ::LL
        INTEGER(4)              ::MM
        REAL(8)                 ::CG1,CG2  
!       *****************************************************************************
        CALL RADIAL$R(GID,NR,R)
        IF(MYUIJKL) THEN
          DO M1=1,2*L+1
            DO M2=1,2*L+1
              DO M3=1,2*L+1
                DO M4=1,2*L+1
                  UIJKL(M1,M3,M2,M4)=0.D0
!                 ==  A=1+(2*L)**2+2*L+2*L
                  DM1=LMSSH+M1
                  DM2=LMSSH+M2
                  DM3=LMSSH+M3
                  DM4=LMSSH+M4
                  DO LM=1,(2*L+1)**2 
                    CALL CLEBSCH(DM3,DM4,LM,CG)
                    RHO1(:)=CG*XALPHA(:,1)*XALPHA(:,1)           
                    LMIR=INT(SQRT(REAL(LM-1)+1.D-3))
                    CALL RADIAL$POISSON(GID,NR,LMIR,RHO1,V)
                    AUX(:)=V(:)*XALPHA(:,1)*XALPHA(:,1)   
                    CALL CLEBSCH(DM1,DM2,LM,CG)
                    RHO2(:)=CG*AUX(:)*R(:)**(2)
                    CALL RADIAL$INTEGRAL(GID,NR,RHO2,SO)
                    UIJKL(M1,M3,M2,M4)=SO+UIJKL(M1,M3,M2,M4)
                  END DO
                END DO
              END DO
            END DO
          END DO
          UIJKL(:,:,:,:)=UIJKL(:,:,:,:)/DIEL
        ELSE
          UIJKL(:,:,:,:)=0.D0       
          CON2=0.625D0
          SF(1)=UPARAMETER
          SF(2)=0.D0
          SF(3)=JPARAMETER*14.D0/(1.D0+CON2)
          SF(4)=0.D0
          SF(5)=SF(3)*CON2
          DO M1=1,2*L+1
            DO M2=1,2*L+1
              DO M3=1,2*L+1
                DO M4=1,2*L+1
                  DM1=LMSSH+M1
                  DM2=LMSSH+M2
                  DM3=LMSSH+M3
                  DM4=LMSSH+M4
                  UIJKL(M1,M3,M2,M4)=0.D0
                  DO LL=0,4,2
                    CON=4*3.1415D0/REAL(2*LL+1,KIND=8)
                    DO MM=-LL,LL
                      CALL CLEBSCH(DM1,DM2,1+LL**(2)+LL-MM,CG1)
                      CALL CLEBSCH(DM3,DM4,1+LL**(2)+LL-MM,CG2) 
                      UIJKL(M1,M3,M2,M4)=CON*CG1*CG2*SF(LL+1) &
       &                                   + UIJKL(M1,M3,M2,M4)
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
        RETURN        
        END
!
!       .............................................................................
        SUBROUTINE LDAUCORRECTION(LOX,LNX,POTUPDOWN,POTAMF,POTFLL,L,LMNXX,NUE,NUENN,&
                                  DIM,UIJKL,UPARAMETER,JPARAMETER,NSPIN,CHOOSEPOT,&
                                  DETOT,DATH,OVER,NOXA,CHOOSENUE)
        IMPLICIT NONE        
        INTEGER(4),INTENT(IN) :: LNX
        INTEGER(4),INTENT(IN) :: LOX(LNX)
        INTEGER(4),INTENT(IN) :: L
        INTEGER(4),INTENT(IN) :: DIM
        REAL(8)   ,INTENT(IN) :: UIJKL(DIM,DIM,DIM,DIM)
        INTEGER(4),INTENT(IN) :: LMNXX
        REAL(8)   ,INTENT(IN) :: NUE(DIM,DIM,2)
        REAL(8)   ,INTENT(IN) :: NUENN(LMNXX,LMNXX)
        REAL(8)   ,INTENT(IN) :: UPARAMETER
        REAL(8)   ,INTENT(IN) :: JPARAMETER
        INTEGER(4),INTENT(IN) :: NSPIN
        INTEGER(4),INTENT(IN) :: CHOOSEPOT
        INTEGER(4),INTENT(IN) :: CHOOSENUE
        REAL(8)   ,INTENT(OUT):: DETOT
        REAL(8)   ,INTENT(OUT):: DATH(LMNXX,LMNXX,NSPIN)
        REAL(8)               :: OCCU(2)
        REAL(8)               :: OCCUGES
        REAL(8)               :: EU
        INTEGER(4)            :: SPIN,SPINS
        INTEGER(4)            :: M1,M2,M3,M4   
        REAL(8)               :: SUMMANDEU1
        REAL(8)               :: EDCFLL
        REAL(8)               :: EDCAMF
        REAL(8)               :: EFLLTOT
        REAL(8)               :: DETOTAMF
        REAL(8),INTENT(INOUT) :: POTUPDOWN(LMNXX,LMNXX,2)
        REAL(8),INTENT(INOUT) :: POTAMF(DIM,DIM,2)
        REAL(8),INTENT(INOUT) :: POTFLL(DIM,DIM,2)
        REAL(8)               :: POTU      
        INTEGER(4)            :: LN1,LN2
        INTEGER(4)            :: LMN1,LMN2
        REAL(8),INTENT(IN)    :: OVER(NOXA,LMNXX)
        INTEGER(4)            :: NOXA
!       *****************************************************************************
        
!       CALCULATE HUBBARD ENERGY        
        EU=0.D0
        DO SPIN=1,2
          DO SPINS=1,2
            DO M1=1,2*L+1
              DO M2=1,2*L+1
                DO M3=1,2*L+1
                  DO M4=1,2*L+1
                    IF(SPIN.EQ.SPINS) THEN
                      SUMMANDEU1=UIJKL(M1,M3,M2,M4)-UIJKL(M1,M3,M4,M2)
                    ELSE
                      SUMMANDEU1=UIJKL(M1,M3,M2,M4)
                    END IF
                    EU=SUMMANDEU1*NUE(M1,M2,SPIN)*NUE(M3,M4,SPINS)+EU
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        
        EU=0.5D0*EU
        PRINT *, 'B=HUBBARD ENNERGY=',EU,' IN EV',EU*27.2144D0
        
!       CALCULATE HUBBARD POTENTIAL
        DO SPIN=1,2
          DO M1=1,2*L+1
            DO M2=1,2*L+1
              POTAMF(M1,M2,SPIN)=0.D0
              DO M3=1,2*L+1
                DO M4=1,2*L+1
                  DO SPINS=1,2
                    IF(SPIN.EQ.SPINS) THEN
                      POTU=UIJKL(M1,M3,M2,M4)-UIJKL(M1,M3,M4,M2)
                    ELSE
                      POTU=UIJKL(M1,M3,M2,M4)
                    END IF
                    POTAMF(M1,M2,SPIN)=POTU*NUE(M3,M4,SPINS)+POTAMF(M1,M2,SPIN)
                  END DO
                END DO
              END DO
              POTFLL(M1,M2,SPIN)=POTAMF(M1,M2,SPIN)
            END DO
          END DO
        END DO
!       =====================================================================
!       == DETERMINE DOUBLE COUNTING CORRECTION                            ==
!       =====================================================================
!       ==  CALCULATE OCCUPATION 
        OCCU(:)=0.0D0
        DO SPIN=1,2
          DO M1=1,2*L+1
            OCCU(SPIN)=OCCU(SPIN)+NUE(M1,M1,SPIN)
          END DO
        END DO
        OCCUGES=OCCU(1)+OCCU(2)
        PRINT*,'OCCU1',OCCU(1)
        PRINT*,'OCCU2',OCCU(2)
        
        IF(CHOOSEPOT.EQ.1) THEN
!         ===================================================================
!         ==  E DC FLL                  
!         ===================================================================
PRINT*,'OCCU ',OCCU
PRINT*,'DC PARMS ',L,UPARAMETER,JPARAMETER,OCCU
          EDCFLL=UPARAMETER/2*(OCCUGES**2-OCCUGES) &
      &         -JPARAMETER/2*(OCCU(1)**2-OCCU(1)) &
      &         -JPARAMETER/2*(OCCU(2)**2-OCCU(2))
!         PRINT*,'E DC OF FLL',EDCFLL*27.2144D0
          EFLLTOT=EU-EDCFLL
          DETOT=EFLLTOT
PRINT *, 'DOUBLE COUNTING ENERGY FLL',EDCFLL,' IN EV:',EDCFLL*27.2144D0 

!         DC POT FLL
          DO SPIN=1,2
            DO M1=1,2*L+1
              POTFLL(M1,M1,SPIN)=POTFLL(M1,M1,SPIN) &
       &             -UPARAMETER*(OCCUGES-0.5D0)+JPARAMETER*(OCCU(SPIN)-0.5D0)
             END DO
          END DO

!         NORM CORRECTION

          IF(CHOOSENUE.EQ.1) THEN          
!           ===================================================================
!           ==  BENGONE
!           ===================================================================
            DO SPIN=1,2  
              LMN1=0 
              DO LN1=1,LNX
                DO M1=1,2*LOX(LN1)+1
                  LMN1=LMN1+1
                  LMN2=0 
                  DO LN2=1,LNX
                    DO M2=1,2*LOX(LN2)+1
                      LMN2=LMN2+1
                      IF(LOX(LN1).EQ.L.AND.LOX(LN2).EQ.L) THEN
                        POTUPDOWN(LMN1,LMN2,SPIN)=POTFLL(M1,M2,SPIN)*NUENN(LMN1,LMN2)
                      END IF
                    END DO
                  END DO
                END DO
              END DO
            END DO
          ELSE IF(CHOOSENUE.EQ.2) THEN
!           ===================================================================
!           ==  PETER
!           ===================================================================
!           POTUPDOWN(:,:,:)=0.D0
            DO SPIN=1,2
              LMN1=0 
              DO LN1=1,LNX
                DO M1=1,2*LOX(LN1)+1
                  LMN1=LMN1+1
                  LMN2=0 
                  DO LN2=1,LNX
                    DO M2=1,2*LOX(LN2)+1
                      LMN2=LMN2+1
                      DO M3=1,NOXA
                        DO M4=1,NOXA
                          IF (LOX(LN1).EQ.L.AND.LOX(LN2).EQ.L) THEN
                            POTUPDOWN(LMN1,LMN2,SPIN)=POTUPDOWN(LMN1,LMN2,SPIN) &
       &                            +OVER(M3,LMN1)*OVER(M4,LMN2)*POTFLL(M3,M4,SPIN)
                          END IF
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END IF
        ELSE IF(CHOOSEPOT.EQ.2) THEN
!         ===================================================================
!         ==  E DC AMF
!         ===================================================================
          EDCAMF=UPARAMETER*OCCU(1)*OCCU(2) &
      &         +(FLOAT(L)/FLOAT(2*L+1))*(UPARAMETER-JPARAMETER)*(OCCU(1)**2+OCCU(2)**2)
          PRINT*,'E DC OF AMF',EDCAMF*27.2144D0
          DETOTAMF=EU-EDCAMF
          DETOT=DETOTAMF
PRINT *, 'DOUBLE COUNTING ENERGY AMF ',EDCAMF,' IN EV:',EDCAMF*27.2144D0 
          PRINT*,  'EAMF',DETOT*27.2144D0       
!         ==  DC POT AMF
          DO SPIN=1,2
            DO M1=1,2*L+1
              POTAMF(M1,M1,SPIN)=POTAMF(M1,M1,SPIN)&
       &           -UPARAMETER*OCCUGES+1/(2*FLOAT(L)+1) &
       &           *(UPARAMETER+2*FLOAT(L)*JPARAMETER)*OCCU(SPIN)
            END DO
          END DO

!         PETER

!         POTUPDOWN(:,:,:)=0.D0
          DO SPIN=1,2
            LMN1=0 
            DO LN1=1,LNX
              DO M1=1,2*LOX(LN1)+1
               LMN1=LMN1+1
               LMN2=0 
               DO LN2=1,LNX
                 DO M2=1,2*LOX(LN2)+1
                   LMN2=LMN2+1
                   DO M3=1,NOXA
                     DO M4=1,NOXA
                       IF(LOX(LN1).EQ.L.AND.LOX(LN2).EQ.L) THEN
                         POTUPDOWN(LMN1,LMN2,SPIN)=POTUPDOWN(LMN1,LMN2,SPIN) &
      &                          +OVER(M3,LMN1)*OVER(M4,LMN2)*POTAMF(M3,M4,SPIN)
                       END IF
                     END DO
                   END DO
                 END DO
               END DO
             END DO
           END DO
          END DO
        END IF
        RETURN
        END SUBROUTINE LDAUCORRECTION
!
!       .............................................................................
        SUBROUTINE READPARAMETERS(UPARAMETER,JPARAMETER,DIEL,MYUIJKL,DEFJP,CHOOSEPOT,CHOOSENUE,RCUT,HML,WAN,WHICHL)
        IMPLICIT NONE
        REAL(8)   ,INTENT(OUT) :: UPARAMETER
        REAL(8)   ,INTENT(OUT) :: JPARAMETER
        REAL(8)   ,INTENT(OUT) :: DIEL           ! DIELECTRIC CONST.
        LOGICAL   ,INTENT(OUT) :: MYUIJKL        ! CALCULATION OF UIJKL
        LOGICAL   ,INTENT(OUT) :: DEFJP          ! DEFINITION OF JPARAMETER
        INTEGER(4),INTENT(OUT) :: CHOOSEPOT !=1  ! CHOOSE DC POTENTIAL (1:FLL; 2:AMF)
        INTEGER(4),INTENT(OUT) :: CHOOSENUE !=2  ! CHOOSE DEFINITION OF NUE (1:O.BENG.;2:PETER)
        REAL(8)   ,INTENT(OUT) :: RCUT
        INTEGER(4),INTENT(OUT) :: HML
        INTEGER(4),INTENT(OUT) :: WAN
        CHARACTER(6)           :: ID
        INTEGER(4)             :: WHICHL(2)
        REAL(8)                :: EV
        INTEGER(4),PARAMETER   :: NFIL=102
!       *****************************************************************************
        OPEN(UNIT=NFIL,FILE='LDAUCONTROL',STATUS='OLD')
        REWIND NFIL
!
!       =============================================================================           
!       == U AND J PARAMETER                                                       ==
!       =============================================================================           
        READ(NFIL,*)ID,UPARAMETER
        IF(ID.NE.'UPARA') THEN 
          CALL ERROR$MSG('ERROR READING PARAMETER FILE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU_READPARAMETERS')
        END IF
        WRITE(*,'(A,F4.2)') 'U',UPARAMETER

        READ(NFIL,*) ID,JPARAMETER         
        IF(ID.NE.'JPARA') THEN 
          CALL ERROR$MSG('ERROR READING PARAMETER FILE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU_READPARAMETERS')
        END IF 
        WRITE(*,'(A,F4.2)') 'J',JPARAMETER
!
!       =============================================================================           
!       == RELATIVE DIELECTRIC CONSTANT                                            ==
!       =============================================================================           
        READ(NFIL,*) ID,DIEL         
        IF(ID.NE.'DIEL') THEN 
          CALL ERROR$MSG('ERROR READING PARAMETER FILE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU_READPARAMETERS')
        END IF
        WRITE(*,'(A,F4.2)') 'DIEL', DIEL
!
!       =============================================================================           
!       == CHOOSE TO OBTAIN U-MATRIX FROM ORBITALS OR NOT                          ==
!       =============================================================================           
        READ(NFIL,*) ID,MYUIJKL
        IF(ID.NE.'MYU') THEN 
          CALL ERROR$MSG('ERROR READING PARAMETER FILE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU_READPARAMETERS')
        END IF
        WRITE(*,'(A,L1)') 'MYU', MYUIJKL
!
!       =============================================================================           
!       == CHOOSE DEFINITION OF J-PARAMETER ?????                                  ==
!       =============================================================================           
        READ(NFIL,*) ID,DEFJP
        IF(ID.NE.'DEFJP') THEN 
          CALL ERROR$MSG('ERROR READING PARAMETER FILE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU_READPARAMETERS')
        END IF
        WRITE(*,'(A,L1)') 'DEFJP', DEFJP 
!
!       =============================================================================           
!       == CHOOSE TYPE OF DOUBLE COUNTING POTENTIAL (1:FLL; 2:AMF)                 ==
!       =============================================================================           
        READ(NFIL,*) ID,CHOOSEPOT
        IF(ID.NE.'CPOT') THEN 
          CALL ERROR$MSG('ERROR READING PARAMETER FILE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU_READPARAMETERS')
        END IF
        WRITE(*,'(A,I2.2)') 'CPOT', CHOOSEPOT 
!
!       =============================================================================           
!       ==  CHOOSE DEFINITION OF DENSITY MATRIC                                    ==
!       ==  (1) PROJECTION ON SPHERE (LIKE BENGONE)                                ==
!       ==  (2) PROJECTION ON ORBITALS                                             ==
!       =============================================================================           
        READ(NFIL,*) ID,CHOOSENUE
        IF(ID.NE.'CNUE') THEN 
          CALL ERROR$MSG('ERROR READING PARAMETER FILE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU_READPARAMETERS')
        END IF
        WRITE(*,'(A,I2.2)') 'CNUE',CHOOSENUE  
!
!       =============================================================================           
!       ==  DEFINE CUTOFF RADIUS FOR ORBITALS                                      ==
!       =============================================================================           
        READ(NFIL,*) ID,RCUT
        IF(ID.NE.'RCUT') THEN 
          CALL ERROR$MSG('ERROR READING PARAMETER FILE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU_READPARAMETERS')
        END IF
        WRITE(*,'(A,F4.2)') 'RCUT',RCUT  
!
!       =============================================================================           
!       ==                                                                         ==
!       =============================================================================           
        READ(NFIL,*) ID,HML
        IF(ID.NE.'HML') THEN 
          CALL ERROR$MSG('ERROR READING PARAMETER FILE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU_READPARAMETERS')
        END IF
        WRITE(*,'(A,I5)') 'HML ',HML  
!
!       =============================================================================           
!       ==  ATOMIC NUMBER  (WAN=WHICH ATOMIC NUMBER)                               ==
!       =============================================================================           
        READ(NFIL,*) ID,WAN
        IF(ID.NE.'WAN') THEN 
          CALL ERROR$MSG('ERROR READING PARAMETER FILE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU_READPARAMETERS')
        END IF
        WRITE(*,'(A,I5)') 'WAN ',WAN  
!        ALLOCATE(WHICHL(HML))
!
!       =============================================================================           
!       ==                                                                         ==
!       =============================================================================           
        READ(NFIL,*) ID,WHICHL
        IF(ID.NE.'WIL') THEN 
          CALL ERROR$MSG('ERROR READING PARAMETER FILE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU_READPARAMETERS')
        END IF
        WRITE(*,'(A,I2.2,I2.2 )') 'WICHL  ',WHICHL  
!
!       =============================================================================           
!       == TRANSFORM FROM EV INTO ATOMIC UNITS                                     ==
!       =============================================================================           
        CALL CONSTANTS('EV',EV)
        UPARAMETER=UPARAMETER*EV
        JPARAMETER=JPARAMETER*EV
        RETURN
        END SUBROUTINE READPARAMETERS
!
!       .............................................................................
        SUBROUTINE LOCALIZEDORBITAL(GID,NR,LNX,L,LOX,ATOMICORBIT,AEPHI,DIM,RCUT,XALPHA)
        IMPLICIT NONE
        LOGICAL   ,INTENT(IN)   :: ATOMICORBIT  ! USE AEPHI/READ FROM FILE
        INTEGER(4),INTENT(IN)   :: GID
        INTEGER(4),INTENT(IN)   :: NR
        INTEGER(4),INTENT(IN)   :: LNX
        INTEGER(4),INTENT(IN)   :: L
        INTEGER(4),INTENT(IN)   :: LOX(LNX)
        REAL(8)   ,INTENT(IN)   :: AEPHI(NR,LNX)
        REAL(8)   ,INTENT(IN)   :: RCUT       !ORBITAL SET TO ZERO BEYOND THE RADIUS RCUT
        INTEGER(4),INTENT(IN)   :: DIM
        REAL(8)   ,INTENT(INOUT):: XALPHA(NR,DIM) 
        INTEGER(4)              :: IR
        INTEGER(4)              :: LN,M
        REAL(8)                 :: XALPHA2(NR) 
        REAL(8)                 :: R(NR) 
!       *****************************************************************************
        IF(DIM.NE.2*L+1) THEN
          CALL ERROR$MSG('INCONSISTENT INPUT')
          CALL ERROR$STOP('LOCALIZEDORBITAL')
        END IF
        CALL RADIAL$R(GID,NR,R)
        IF(ATOMICORBIT) THEN
!         ===========================================================================
!         ==  USE AEPHI AS LOCALIZED ORBITAL                                       ==
!         ===========================================================================
!         == PICK THE FIRST PARTIAL WAVE WITH THE CORRECT L AS LOCAL ORBITAL
          DO LN=1,LNX
            IF(LOX(LN).EQ.L) THEN
              DO M=1,2*L+1
                XALPHA(:,M)=AEPHI(:,LN)    
              ENDDO
              EXIT
            END IF
          ENDDO
        ELSE
!         ===========================================================================
!         ==   READ LOCALIZED ATOMIC ORBITAL FROM FILE                             ==
!         ===========================================================================
          OPEN(UNIT=101,FILE='XALPHA.DAT',STATUS='OLD')
          REWIND 101
          DO IR=1,NR
            READ(101,*) XALPHA2(IR)
          END DO
          DO M=1,DIM
            XALPHA(:,M)=XALPHA2(:)
          END DO
        END IF
!
!       ===========================================================================
!       ==   CUT OFF THE ORBITAL BEYOND RCUT                                     ==
!       ===========================================================================
        DO IR=1,NR
          IF(R(IR).GT.RCUT) THEN
            XALPHA(IR,:)=0.D0
          END IF
        END DO
        RETURN
        END SUBROUTINE LOCALIZEDORBITAL
!
!       ............................................................................
        SUBROUTINE LDAPLUSU$NEW(NSP_)
!       **                                                                        **
!       **                                                                        **
!       **                                                                        **
        USE LDAPLUSU_MODULE
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: NSP_
!       ****************************************************************************
        IF(NSP.NE.0) THEN
          IF(NSP_.NE.NSP) THEN
            CALL ERROR$STOP('LDAPLUSU$NEW')
            RETURN
          END IF
        END IF
        NSP=NSP_
        ISP=0
        ALLOCATE(THISARRAY(NSP))
        DO ISP=1,NSP
          ALLOCATE(THISARRAY(ISP)%NORB(10))
          THISARRAY(ISP)%NORB(:)=10
          THISARRAY(ISP)%NORB(1:2)=0 ! PER DEFAULT CORRELATE ONLY D AND F-SHELLS
        ENDDO
        RETURN
        END
!
!       ............................................................................
        SUBROUTINE LDAPLUSU$SELECT(ISP_)
!       **                                                                        **
!       **  SELECT AN ATOM TYPE OR UNSELECT WITH ISP=0                            **
!       **                                                                        **
!       **                                                                        **
        USE LDAPLUSU_MODULE
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: ISP_
!       ****************************************************************************
        ISP=ISP_
        IF(ISP.GT.NSP) THEN
          IF(NSP.EQ.0) THEN
            CALL ERROR$STOP('LDAPLUSU NOT INITIALIZED')
            CALL ERROR$STOP('LDAPLUSU$SELECT')
          END IF
          CALL ERROR$STOP('ISP OUT OF RANGE')
          CALL ERROR$STOP('LDAPLUSU$SELECT')
        END IF
        IF(ISP.EQ.0) THEN
          NULLIFY(THIS)
        ELSE
          THIS=>THISARRAY(ISP)
        END IF
        RETURN
        END
!
!       .............................................................................
        SUBROUTINE LDAPLUSU$SETR8(ID,VAL)
        USE LDAPLUSU_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID
        REAL(8)     ,INTENT(IN) :: VAL
!       *****************************************************************************
        IF(ISP.EQ.0) THEN
          CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF

        IF(ID.EQ.'UPAR') THEN
          THIS%UPAR=VAL
          THIS%USEUPAR=.TRUE.
          IF(THIS%USEDIEL) THEN
            CALL ERROR$MSG('DIEL HAS ALREADY BEEN SET')
            CALL ERROR$MSG('UPAR AND DIEL CANNOT BE USED SIMULTANEOUSLY')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('LDAPLUSU$SETR8')
          END IF
        ELSE IF(ID.EQ.'JPAR') THEN
          THIS%JPAR=VAL
          THIS%USEJPAR=.TRUE.
          IF(THIS%USEDIEL) THEN
            CALL ERROR$MSG('DIEL HAS ALREADY BEEN SET')
            CALL ERROR$MSG('JPAR AND DIEL CANNOT BE USED SIMULTANEOUSLY')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('LDAPLUSU$SETR8')
          END IF
        ELSE IF(ID.EQ.'DIEL') THEN
          THIS%DIEL=VAL
          THIS%USEDIEL=.TRUE.
          IF(THIS%USEUPAR.OR.THIS%USEJPAR) THEN
            CALL ERROR$MSG('UPAR OR JPAR HAVE ALREADY BEEN SET')
            CALL ERROR$MSG('DIEL AND UPAR OR JPAR CANNOT BE USED SIMULTANEOUSLY')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('LDAPLUSU$SETR8')
          END IF
        ELSE IF(ID.EQ.'RCUT') THEN
          THIS%RCUT=VAL
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF
        RETURN
        END 
!
!       .............................................................................
        SUBROUTINE LDAPLUSU$SETL4(ID,VAL)
        USE LDAPLUSU_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID
        LOGICAL(4)  ,INTENT(IN) :: VAL
!       *****************************************************************************
        IF(ID.EQ.'ON') THEN
          TON=VAL
          RETURN
        END IF
        IF(ISP.EQ.0) THEN
          CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF
        IF(ID.EQ.'ACTIVE') THEN
          THIS%TON=VAL
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETL4')
        END IF
        RETURN
        END 
!
!       .............................................................................
        SUBROUTINE LDAPLUSU$SETI4A(ID,LENG,VAL)
        USE LDAPLUSU_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID
        INTEGER(4)  ,INTENT(IN) :: LENG
        INTEGER(4)  ,INTENT(IN) :: VAL(LENG)
!       *****************************************************************************
        IF(ISP.EQ.0) THEN
          CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETI4A')
        END IF
        IF(ID.EQ.'NCORROFL') THEN
          DEALLOCATE(THISARRAY(ISP)%NORB)
          ALLOCATE(THISARRAY(ISP)%NORB(LENG))
          THIS%NORB=VAL
!
        ELSE IF(ID.EQ.'MAINLN') THEN
          THIS%MAINLN=VAL
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETI4A')
        END IF
        RETURN
        END 
!
!       .............................................................................
        SUBROUTINE LDAPLUSU$SETCH(ID,VAL)
        USE LDAPLUSU_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID
        CHARACTER(*),INTENT(IN) :: VAL
!       *****************************************************************************
        IF(ISP.EQ.0) THEN
          CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETCH')
        END IF
        IF(ID.EQ.'DCTYPE') THEN
          DCTYPE=VAL
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETCH')
        END IF
        RETURN
        END 
!
!       ............................................................................
        SUBROUTINE LDAPLUSU$ETOT(ISP_,LMNXX,NDIMD,DENMAT,ETOT,DATH_)
!       **                                                                        **
!       **  THIS IS THE MAIN DRIVER ROUTINE FOR THE LDA+U CORRECTION              **
!       **                                                                        **
!       **                                                                        **
        USE LDAPLUSU_MODULE
        IMPLICIT NONE
        INTEGER(4),INTENT(IN)  :: ISP_
        INTEGER(4),INTENT(IN)  :: LMNXX
        INTEGER(4),INTENT(IN)  :: NDIMD
        COMPLEX(8),INTENT(IN)  :: DENMAT(LMNXX,LMNXX,NDIMD)
        REAL(8)   ,INTENT(OUT) :: ETOT
        REAL(8)   ,INTENT(OUT) :: DATH_(LMNXX,LMNXX,NDIMD)
        CHARACTER(8),PARAMETER :: CHITYPE='FROMPHI'
        INTEGER(4),PARAMETER   :: LRX=4
        INTEGER(4),PARAMETER   :: L=2
        COMPLEX(8),ALLOCATABLE :: DATH(:,:,:)
        LOGICAL(4)             :: TUFROMPARMS=.FALSE.
        INTEGER(4)             :: GID
        INTEGER(4)             :: NR
        INTEGER(4)             :: LNX
        INTEGER(4),ALLOCATABLE :: LOX(:)
        INTEGER(4)             :: NCHI
        REAL(8)   ,ALLOCATABLE :: PHITOCHI(:,:)
        INTEGER(4)             :: LNXPHI
        INTEGER(4),ALLOCATABLE :: LOXPHI(:)
        INTEGER(4)             :: NPHI
        REAL(8)   ,ALLOCATABLE :: U(:,:,:,:)
        COMPLEX(8),ALLOCATABLE :: RHO(:,:,:,:)
        COMPLEX(8),ALLOCATABLE :: HAM(:,:,:,:)
        COMPLEX(8),ALLOCATABLE :: MATSS(:,:,:,:)
        INTEGER(4)             :: IS1,IS2,I,J,LN,M
        REAL(8)                :: PI,SVAR
        REAL(8)                :: UPAR,JPAR,JPAR1
!       ****************************************************************************
        IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('LDAPLUSU$ETOT')
        CALL LDAPLUSU$SELECT(ISP_)
        IF(.NOT.THIS%TON) RETURN
        PI=4.D0*DATAN(1.D0)
        
!       ============================================================================
!       ==  CONSTRUCT LOCAL ORBITALS                                              ==
!       ============================================================================
        IF(.NOT.THIS%TINI) THEN
          CALL LDAPLUSU_CHIFROMPHI()
PRINT*,'TINI ',THIS%TINI
PRINT*,'TON ',THIS%TON
PRINT*,'LNXCHI ',THIS%LNXCHI
PRINT*,'NR     ',THIS%NR
PRINT*,'LOXCHI ',THIS%LOXCHI
PRINT*,'NORB   ',THIS%NORB
PRINT*,'NCHI   ',THIS%NCHI
PRINT*,'DIEL   ',THIS%DIEL
PRINT*,'UPAR   ',THIS%UPAR
PRINT*,'JPAR   ',THIS%JPAR

!         ==========================================================================
!         ==  CALCULATE SMALL U-TENSOR                                            ==
!         ==========================================================================
          GID=THIS%GID
          NR=THIS%NR
          LNX=THIS%LNXCHI
          ALLOCATE(LOX(LNX))
          LOX=THIS%LOXCHI
          NCHI=THIS%NCHI
          ALLOCATE(THIS%ULITTLE(LRX+1,LNX,LNX,LNX,LNX))
          CALL LDAPLUSU_ULITTLE(GID,NR,LRX,LNX,LOX,THIS%CHI,THIS%ULITTLE)
          IF(THIS%USEDIEL.AND.(THIS%USEUPAR.OR.THIS%USEJPAR)) THEN
            CALL ERROR$MSG('DIEL AND (UPAR.OR.JPAR) ARE INCOMPATIBLE')
            CALL ERROR$L4VAL('USEDIEL',THIS%USEDIEL)
            CALL ERROR$L4VAL('USEUPAR',THIS%USEUPAR)
            CALL ERROR$L4VAL('USEJPAR',THIS%USEJPAR)
            CALL ERROR$STOP('LDAPLUSU$ETOT')
          END IF
          IF(THIS%USEDIEL) THEN
            THIS%ULITTLE=THIS%ULITTLE/THIS%DIEL
          ELSE IF(THIS%USEUPAR.OR.THIS%USEJPAR) THEN
            CALL LDAPLUSU_MODULITTLEWITHPARMS(LNX,LOX,LRX,THIS%USEUPAR,THIS%UPAR &
       &                           ,THIS%USEJPAR,THIS%JPAR,THIS%MAINLN,THIS%ULITTLE)
          END IF
        ELSE
          GID=THIS%GID
          NR=THIS%NR
          LNX=THIS%LNXCHI
          ALLOCATE(LOX(LNX))
          LOX=THIS%LOXCHI
          NCHI=THIS%NCHI
        END IF
        CALL SETUP$SELECT(ISP)
        CALL SETUP$LNX(ISP,LNXPHI)
        ALLOCATE(LOXPHI(LNXPHI))
        CALL SETUP$LOFLN(ISP,LNXPHI,LOXPHI)
        NPHI=SUM(2*LOXPHI+1)

!
!       ==========================================================================
!       ==  DOWNFOLD                                                            ==
!       ==========================================================================
        ALLOCATE(PHITOCHI(NCHI,NPHI))
        ALLOCATE(RHO(NCHI,NCHI,2,2))
        ALLOCATE(HAM(NCHI,NCHI,2,2))
        ALLOCATE(MATSS(NPHI,NPHI,2,2))
        ALLOCATE(DATH(NPHI,NPHI,NDIMD))
        ALLOCATE(U(NCHI,NCHI,NCHI,NCHI))
!       == TRANSFORM FROM TOTAL SPIN TO SPIN-SPIN
        CALL LDAPLUSU_SPINDENMAT('FORWARD',NDIMD,NPHI,DENMAT(1:NPHI,1:NPHI,:),MATSS)
!       == TRANSFORM FROM PHI TO CHI
        CALL LDAPLUSU_MAPTOCHI(LNX,LOX,NCHI,LNXPHI,LOXPHI,NPHI,PHITOCHI)
PRINT*,'===================== phitochi ======================'
do i=1,nchi
  WRITE(*,FMT='(I3,100F8.3)')nchi,phitochi(I,:)
enddo

        DO IS1=1,2
          DO IS2=1,2
            RHO(:,:,IS1,IS2)=MATMUL(PHITOCHI,MATMUL(MATSS(:,:,IS1,IS2),TRANSPOSE(PHITOCHI)))
          ENDDO
        ENDDO

DO IS1=1,2
  DO IS2=1,2
    IF(SUM(ABS(RHO(:,:,IS1,IS2))).LT.1.D-3) CYCLE
    PRINT*,'===================== RHO FOR SPIN',IS1,IS2,' ======================'
    I=0
    DO LN=1,LNX
      DO M=1,2*LOX(LN)+1
        I=I+1
        WRITE(*,FMT='(I3,100F8.3)')LOX(LN),REAL(RHO(I,:,IS1,IS2))
      ENDDO
    ENDDO
  ENDDO
ENDDO
do is1=1,2
  svar=0.d0
  do ln=1,nchi
    svar=svar+rho(ln,ln,is1,is1)
  enddo
  print*,'charge= ',svar,' for spin ',is1
enddo
!
!       ==========================================================================
!       ==  CALCULATE U-TENSOR                                                  ==
!       ==========================================================================
        CALL LDAPLUSU_UTENSOR(LRX,NCHI,LNX,LOX,THIS%ULITTLE,U)
        UPAR=0.D0
        JPAR=0.D0
        DO I=1,NCHI
          DO J=1,NCHI
            UPAR=UPAR+U(I,J,I,J)
            IF(I.NE.J)JPAR=JPAR+U(I,J,J,I)
          ENDDO
        ENDDO
        UPAR=UPAR/REAL(NCHI)**2
        JPAR=JPAR/REAL(NCHI*(NCHI-1))*7.D0/5.D0
PRINT*,'UPARAMETER[EV]    ',UPAR*27.211D0 ,'UPARAMETER    ',UPAR
PRINT*,'JPARAMETER[EV](1) ',JPAR*27.211D0 ,'JPARAMETER(1) ',JPAR

!
!       ============================================================================
!       ==  CALCULATE LDA+U TOTAL ENERGY CONTRIBUTION                             ==
!       ============================================================================
        CALL LDAPLUSU_ETOT(DCTYPE,L,UPAR,JPAR,NCHI,U,RHO,ETOT,HAM)
!
!       ============================================================================
!       ==  UPFOLD                                                                ==
!       ============================================================================
!       == TRANSFORM FROM CHI TO PHI ===============================================
        DO IS1=1,2
          DO IS2=1,2
            MATSS(:,:,IS1,IS2)=MATMUL(TRANSPOSE(PHITOCHI),MATMUL(HAM(:,:,IS1,IS2),PHITOCHI))
          ENDDO
        ENDDO
!
!       == TRANSFORM FROM (SPIN,SPIN) TO (TOTAL,SPIN) ==============================
        CALL LDAPLUSU_SPINDENMAT('BACK',NDIMD,LMNXX,DATH,MATSS)
!
!       == MAKE REAL (THIS IS A FUDGE TO BE FIXED IN AUGMENTATION!)
        DATH_(:,:,:)=(0.D0,0.D0)
        DATH_(:NPHI,:NPHI,:)=REAL(DATH)
!!$DO IS1=1,NDIMD
!!$PRINT*,'===================== DATH FOR SPIN',IS1,IS2,' ======================'
!!$I=0
!!$DO LN=1,LNXPHI
!!$DO M=1,2*LOXPHI(LN)+1
!!$I=I+1
!!$WRITE(*,FMT='(I3,100F8.3)')LOXPHI(LN),REAL(DATH(I,:,IS1))
!!$ENDDO
!!$ENDDO
!!$ENDDO
!
!       ============================================================================
!       ==  UNSELECT LDAPLUSU                                                     ==
!       ============================================================================
        DEALLOCATE(U)
        DEALLOCATE(LOXPHI)
        DEALLOCATE(LOX)
        DEALLOCATE(MATSS)
        DEALLOCATE(RHO)
        DEALLOCATE(HAM)
        CALL LDAPLUSU$SELECT(0)
                                     CALL TRACE$POP()
        RETURN
        END
!
!       ............................................................................
        SUBROUTINE LDAPLUSU_CHIFROMPHI()
!       **                                                                        **
        USE LDAPLUSU_MODULE
        IMPLICIT NONE
        REAL(8)               :: RCUT
        INTEGER(4)            :: GID
        INTEGER(4)            :: NR
        INTEGER(4)            :: LNX
        INTEGER(4),ALLOCATABLE:: LOX(:)
        REAL(8)   ,ALLOCATABLE:: PHI(:,:)
        INTEGER(4)            :: LNXCHI 
        INTEGER(4),ALLOCATABLE:: LOXCHI(:)
        REAL(8)   ,ALLOCATABLE:: CHI(:,:)
        REAL(8)   ,ALLOCATABLE:: CHI1(:,:)
        REAL(8)   ,ALLOCATABLE:: A(:,:)
        REAL(8)   ,ALLOCATABLE:: A1(:,:)
        REAL(8)   ,ALLOCATABLE:: MAT(:,:)
        REAL(8)   ,ALLOCATABLE:: MATINV(:,:)
        REAL(8)   ,ALLOCATABLE:: R(:)        
        REAL(8)   ,ALLOCATABLE:: G(:)        
        REAL(8)   ,PARAMETER  :: RCG=0.3D0
        REAL(8)   ,ALLOCATABLE:: AUX(:)
        REAL(8)               :: SVAR1,SVAR2
        INTEGER(4)            :: IR
        INTEGER(4)            :: NX,N,LX,L,LN,LNCHI,LN0,NOFL,ISVAR
        INTEGER(4)            :: N1,N2,LN1,LN2,L1,L2
!       ****************************************************************************
                              CALL TRACE$PUSH('LDAPLUSU_CHIFROMPHI')
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        THIS%GID=GID
        CALL RADIAL$GETI4(GID,'NR',NR)
        THIS%NR=NR
        CALL SETUP$LNX(ISP,LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$LOFLN(ISP,LNX,LOX)
        ALLOCATE(PHI(NR,LNX))
        CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,PHI)
!
!       ============================================================================
!       ==  DIVIDE IN HEAD AN TAIL FUNCTIONS                                      ==
!       ============================================================================
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        ALLOCATE(AUX(NR))
!       == COUNT NUMBER OF PARTIAL WAVES PER L AS AS ESTIMATE FOR THE NUMBER OF CHI
!       == COUNT ONLY THOSE ANGULAR MOMENTUM CHANNELS THAT HAVE CORRELATED ORBITALS
        LX=MAXVAL(LOX)
        LNXCHI=0
        DO L=0,LX
          IF(THIS%NORB(L+1).EQ.0) CYCLE  ! FORGET CHANNELS WITHOUT CORRELATED ORBITALS
          NOFL=0
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            NOFL=NOFL+1
            LNXCHI=LNXCHI+1
          ENDDO
!         == LIMIT THE NUMBER OF CORRELATED ORBITALS TO THE NUMBER OF PARTIAL WAVES 
!         == MINUS ONE, SO THAT ONE TAIL FUNCTION IS LEFT. 
!         == HOWEVER, IF THERE IS ONLY A SINGLE PARTIAL WAVE, LEAVE THAT ONE.
          IF(THIS%NORB(L+1).GT.NOFL-1) THEN
            ISVAR=MAX(1,NOFL-1)            
            THIS%NORB(L+1)=MIN(THIS%NORB(L+1),ISVAR) 
          END IF
        ENDDO
!
!       == ORDER ACCORDING TO L ====================================================
        ALLOCATE(LOXCHI(LNXCHI))        
        ALLOCATE(CHI(NR,LNXCHI))        
!       == |CHI_I>=SUM_I |PHI_J>A(J,I)
        ALLOCATE(A(LNX,LNXCHI))        
        A(:,:)=0.D0
        LNCHI=0
        DO L=0,LX
          IF(THIS%NORB(L+1).EQ.0) CYCLE  ! SHELLS WITHIOUT LOCAL ORBITALS ARE IRRELEVANT
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE  !CONSIDER ONLY PARTIAL WAVES WITH THE CORRECT L
            LNCHI=LNCHI+1
            LOXCHI(LNCHI)=LOX(LN)
            CHI(:,LNCHI)=PHI(:,LN)
            A(LN,LNCHI)=1.D0
          ENDDO
        ENDDO
!
!       == MAKE HEAD FUNCTION ANTIBONDING WITH NODE AT RCUT ==========================
        RCUT=THIS%RCUT
        L=-1
        DO LN=1,LNXCHI
          IF(LOXCHI(LN).NE.L) NOFL=0    ! RESET ORBITAL COUNTER FOR EACH L
          L=LOXCHI(LN)
          NOFL=NOFL+1
          IF(NOFL.GT.THIS%NORB(L+1)) CYCLE  !ONLY CORRELATED ORBITALS WILL BE LOCALIZED
          IF(LN+1.GT.LNXCHI) CYCLE          !CHECK IF THERE IS A TAIL FUNCTION LEFT
          IF(LOXCHI(LN+1).NE.L) CYCLE       !CHECK IF THERE IS A TAIL FUNCTION LEFT
!         == IMPOSE NODE CONDITION======================================================
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN),RCUT,SVAR1)
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN+1),RCUT,SVAR2)
          CHI(:,LN)=CHI(:,LN)*SVAR2-CHI(:,LN+1)*SVAR1
          A(:,LN)=A(:,LN)*SVAR2-A(:,LN+1)*SVAR1
!         == CUT AT RCUT           ===================================================
          DO IR=1,NR
            IF(R(IR).GT.RCUT) CHI(IR,LN)=0.D0
          ENDDO
!         == ORTHOGONALIZE TO THE LOWER HEAD FUNCTIONS ============================
          DO LN1=LN-NOFL+1,LN-1
            AUX(:)=CHI(:,LN)*CHI(:,LN1)*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
            CHI(:,LN)=CHI(:,LN)-CHI(:,LN1)*SVAR1
            A(:,LN)=A(:,LN)-A(:,LN1)*SVAR1
          ENDDO
!         == NORMALIZE HEAD FUNCTION =============================================
          AUX(:)=CHI(:,LN)**2*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          SVAR1=1.D0/SQRT(SVAR1)
          CHI(:,LN)=CHI(:,LN)*SVAR1
          A(:,LN)=A(:,LN)*SVAR1
        END DO          
!
!       == DELOCALIZE TAIL FUNCTIONS ====================================================
        ALLOCATE(CHI1(NR,LNXCHI))
        CHI1(:,:)=CHI(:,:)
        ALLOCATE(A1(LNX,LNXCHI))        
        A1(:,:)=A(:,:)
        ALLOCATE(G(NR))
        G(:)=EXP(-(R(:)/RCG)**2)
        L=-1
        DO LN=1,LNXCHI
          IF(LOXCHI(LN).NE.L) THEN
            L=LOXCHI(LN)
            LN0=LN
            CYCLE
          END IF  
!         == MINIMIZE CONTRIBUTION NEAR THE CENTER 
          DO LN2=LN0,LN-1
            AUX(:)=G(:)*CHI1(:,LN)*CHI1(:,LN2)*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
            AUX(:)=G(:)*CHI1(:,LN2)**2*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
            SVAR1=SVAR1/SVAR2
            CHI1(:,LN)=CHI1(:,LN)-CHI1(:,LN2)*SVAR1
            A1(:,LN)=A1(:,LN)-A1(:,LN2)*SVAR1
          ENDDO
        ENDDO
!       = NOW MAP ONLY TAIL FUNCTIONS BACK AND LEAVE HEAD FUNCTIONS UNTOUCHED        
        L=-1
        DO LN=1,LNXCHI
          IF(LOXCHI(LN).NE.L) THEN
            L=LOXCHI(LN)
            NOFL=0
          END IF  
          NOFL=NOFL+1
          IF(NOFL.LE.THIS%NORB(L+1)) CYCLE ! LEAVE HEAD FUNCTIONS ALONE
          CHI(:,LN)=CHI1(:,LN)
          A(:,LN)=A1(:,LN)
        ENDDO
!
!       === CONSTRUCT TRANSFORMATION MATRIX FROM A ========================================
        LX=MAXVAL(LOXCHI)
        DO L=1,LX
!
          NX=0
          DO LN=1,LNXCHI
            IF(LOXCHI(LN).EQ.L) NX=NX+1
          ENDDO
          IF(NX.EQ.0) CYCLE   

          ALLOCATE(MAT(NX,NX))
          ALLOCATE(MATINV(NX,NX))
!          
          N1=0
          DO LN1=1,LNXCHI
            L1=LOXCHI(LN1)
            IF(L1.NE.L) CYCLE
            N1=N1+1
            N2=0
            DO LN2=1,LNX
              L2=LOX(LN2)
              IF(L2.NE.L) CYCLE
              N2=N2+1
              MAT(N2,N1)=A(LN2,LN1)
            ENDDO
          ENDDO
          CALL LIB$INVERTR8(NX,MAT,MATINV)
          N1=0
          DO LN1=1,LNXCHI
            L1=LOXCHI(LN1)
            IF(L1.NE.L) CYCLE
            N1=N1+1
            N2=0
            DO LN2=1,LNX
              L2=LOX(LN2)
              IF(L2.NE.L) CYCLE
              N2=N2+1
              A(LN2,LN1)=MATINV(N1,N2)   ! a it transposed so that the indices match
            ENDDO
          ENDDO
          DEALLOCATE(MAT)
          DEALLOCATE(MATINV)
        ENDDO
!
!       === REMOVE TAIL FUNCTIONS                         ======================================
        LN1=0
        L=-1
        DO LN=1,LNXCHI
          IF(L.NE.LOXCHI(LN)) THEN
            L=LOXCHI(LN)
            N=0
            NX=THIS%NORB(L+1)
          END IF
          N=N+1
          IF(N.LE.NX) THEN
            LN1=LN1+1
            LOXCHI(LN1)=LOXCHI(LN)
            CHI(:,LN1)=CHI(:,LN)
!           -- A HAS BEEN INVERTED. THEREFORE THE CHI-INDEX IS LEFT 
!           -- AND THE PHI-INDEX IN ON THE RIGHT HAND SIDE
            A(:,LN1)=A(:,LN)  
          END IF
        ENDDO
        LNXCHI=LN1

!PRINT*,'== WRITE CHI.DAT'
!OPEN(UNIT=109,FILE='CHI.DAT',FORM='FORMATTED')
!DO IR=1,NR
!  WRITE(109,*)R(IR),CHI(IR,1:LNXCHI)
!ENDDO
!CLOSE(109)
!STOP 'FORCED AFTER WRITING CHI'
!
!       ============================================================================
!       ==  CUT OF SUPPORT FUNCTIONS                                              ==
!       ============================================================================
        THIS%LNXCHI=LNXCHI
        ALLOCATE(THIS%LOXCHI(LNXCHI))
        THIS%LOXCHI=LOXCHI(1:LNXCHI)
        THIS%NCHI=SUM(2*LOXCHI(1:LNXCHI)+1)
        ALLOCATE(THIS%CHI(NR,LNXCHI))
        THIS%CHI=CHI(:,1:LNXCHI)
        ALLOCATE(THIS%DOWNFOLD(LNXCHI,LNX))
        do ln=1,lnxchi
           THIS%DOWNFOLD(ln,:)=A(:,ln)
        enddo
!
!       ============================================================================
!       ==  CLEAN UP                                                              ==
!       ============================================================================
        DEALLOCATE(R)
                              CALL TRACE$POP()
        RETURN
        END
!
!       ............................................................................
        SUBROUTINE LDAPLUSU_SPINDENMAT(ID,NDIMD,LMNXX,MAT1,MAT2)
!       **                                                                        **
!       ** IF="FORWARD": CONVERTS DENSITY MATRIX MAT1                             **
!       **  FROM (TOTAL,SPIN) REPRESENTATION INTO (SPIN,SPIN) REPRESENTATION MAT2 **
!       ** IF="BACK": CONVERTS HAMILTON MATRIX MAT2                               **
!       **  FROM (TOTAL,SPIN) REPRESENTATION INTO (SPIN,SPIN) REPRESENTATION MAT2 **
!       **                                                                        **
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID
        INTEGER(4),INTENT(IN)   :: NDIMD
        INTEGER(4),INTENT(IN)   :: LMNXX
        COMPLEX(8),INTENT(INOUT):: MAT1(LMNXX,LMNXX,NDIMD)
        COMPLEX(8),INTENT(INOUT):: MAT2(LMNXX,LMNXX,2,2)
        COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
!       ****************************************************************************
                              CALL TRACE$PUSH('LDAPLUSU_SPINDENMAT')
        IF(ID.EQ.'FORWARD') THEN
          MAT2(:,:,:,:)=(0.D0,0.D0)
          IF(NDIMD.EQ.1) THEN
            MAT2(:,:,1,1)=0.5D0*MAT1(:,:,1)
            MAT2(:,:,2,2)=0.5D0*MAT1(:,:,1)
          ELSE IF(NDIMD.EQ.2) THEN
            MAT2(:,:,1,1)=0.5D0*(MAT1(:,:,1)+MAT1(:,:,2))
            MAT2(:,:,2,2)=0.5D0*(MAT1(:,:,1)-MAT1(:,:,2))
          ELSE IF (NDIMD.EQ.4) THEN
            MAT2(:,:,1,1)=0.5D0*(MAT1(:,:,1)+MAT1(:,:,4))
            MAT2(:,:,2,2)=0.5D0*(MAT1(:,:,1)-MAT1(:,:,4))
            MAT2(:,:,1,2)=0.5D0*(MAT1(:,:,2)-CI*MAT1(:,:,3))
            MAT2(:,:,2,1)=0.5D0*(MAT1(:,:,2)+CI*MAT1(:,:,3))
          ELSE
            CALL ERROR$MSG('NDIMD OUT OF RANGE')
            CALL ERROR$MSG('LDAPLUSU_SPINDENMAT')
          END IF
        ELSE IF(ID.EQ.'BACK') THEN
          MAT1(:,:,:)=(0.D0,0.D0)
          IF(NDIMD.EQ.1) THEN
            MAT1(:,:,1)=0.5D0*(MAT2(:,:,1,1)+MAT2(:,:,2,2))
          ELSE IF(NDIMD.EQ.2) THEN
            MAT1(:,:,1)=0.5D0*(MAT2(:,:,1,1)+MAT2(:,:,2,2))
            MAT1(:,:,2)=0.5D0*(MAT2(:,:,1,1)-MAT2(:,:,2,2))
          ELSE IF (NDIMD.EQ.3) THEN
            MAT1(:,:,1)=0.5D0*(MAT2(:,:,1,1)+MAT2(:,:,2,2))
            MAT1(:,:,2)=0.5D0*(MAT2(:,:,1,2)+MAT2(:,:,2,1))
            MAT1(:,:,3)=-0.5D0*CI*(MAT2(:,:,1,2)-MAT2(:,:,2,1))
            MAT1(:,:,4)=0.5D0*(MAT2(:,:,1,1)-MAT2(:,:,2,2))
          ELSE
            CALL ERROR$MSG('NDIMD OUT OF RANGE')
            CALL ERROR$MSG('LDAPLUSU_SPINDENMAT')
          END IF
        ELSE
          CALL ERROR$MSG('ID MUST BE EITHER "FORWARD" OR "BACK"')
          CALL ERROR$MSG('LDAPLUSU_SPINDENMAT')
        END IF
                              CALL TRACE$POP()
        RETURN
        END

!
!       ............................................................................
        SUBROUTINE LDAPLUSU_MAPTOCHI(LNXCHI,LOXCHI,LMNXCHI,LNXPHI,LOXPHI,LMNXPHI,DOWNFOLD)
!       **                                                                        **
        USE LDAPLUSU_MODULE
        IMPLICIT NONE
        INTEGER(4),INTENT(IN)   :: LNXCHI
        INTEGER(4),INTENT(IN)   :: LOXCHI(LNXCHI)
        INTEGER(4),INTENT(IN)   :: LMNXCHI
        INTEGER(4),INTENT(IN)   :: LNXPHI
        INTEGER(4),INTENT(IN)   :: LOXPHI(LNXPHI)
        INTEGER(4),INTENT(IN)   :: LMNXPHI
        REAL(8)   ,INTENT(OUT)  :: DOWNFOLD(LMNXCHI,LMNXPHI)
        INTEGER(4)              :: LMN1,LN1,L1,LMN2,LN2,L2,M
!       ****************************************************************************
                              CALL TRACE$PUSH('LDAPLUSU_MAPTOCHI')
        DOWNFOLD=0.D0
        LMN1=0
        DO LN1=1,LNXCHI
          L1=LOXCHI(LN1)
          LMN2=0
          DO LN2=1,LNXPHI
            L2=LOXPHI(LN2)
            IF(L1.EQ.L2) THEN
              DO M=1,2*L1+1
                DOWNFOLD(LMN1+M,LMN2+M)=THIS%DOWNFOLD(LN1,LN2)
              ENDDO
            END IF
            LMN2=LMN2+2*L2+1
          ENDDO
          LMN1=LMN1+2*L1+1
        ENDDO
                              CALL TRACE$POP()
        RETURN
        END
!
!       ............................................................................
        SUBROUTINE LDAPLUSU_ULITTLE(GID,NR,LRX,LNX,LOX,CHI,ULITTLE)
!       ** CALCULATES THE INTERACTION ENERGY                                      **
!       **                                                                        **
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: GID
        INTEGER(4),INTENT(IN) :: NR
        INTEGER(4),INTENT(IN) :: LRX
        INTEGER(4),INTENT(IN) :: LNX
        INTEGER(4),INTENT(IN) :: LOX(LNX)
        REAL(8)   ,INTENT(IN) :: CHI(NR,LNX)
        REAL(8)   ,INTENT(OUT):: ULITTLE(LRX+1,LNX,LNX,LNX,LNX)
        INTEGER(4)            :: LN1,LN2,LN3,LN4
        INTEGER(4)            :: L
        INTEGER(4)            :: LMIN,LMAX,ISVAR1,ISVAR2
        REAL(8)               :: RHO(NR)
        REAL(8)               :: POT(NR)
        REAL(8)               :: AUX(NR)
        REAL(8)               :: SVAR
        REAL(8)               :: R(NR)
!       ****************************************************************************
                              CALL TRACE$PUSH('LDAPLUSU_ULITTLE')
        CALL RADIAL$R(GID,NR,R)
        ULITTLE=0.D0
        DO LN1=1,LNX
          DO LN2=LN1,LNX
            RHO(:)=CHI(:,LN1)*CHI(:,LN2)
!           == USE SELECTION RULE (NOT TO SAVE TIME HERE, BUT LATER FOR THE U-TENSOR)
            ISVAR1=ABS(LOX(LN1)+LOX(LN2))  
            ISVAR2=ABS(LOX(LN1)-LOX(LN2))
            LMIN=MIN(ISVAR1,ISVAR2)
            LMAX=MAX(ISVAR1,ISVAR2)
            LMAX=MIN(LMAX,LRX)
            DO L=LMIN,LMAX
              CALL RADIAL$POISSON(GID,NR,L,RHO,POT)
              POT(:)=POT(:)*R(:)**2
              DO LN3=1,LNX
                DO LN4=LN3,LNX
                  ISVAR1=ABS(LOX(LN3)+LOX(LN4))
                  ISVAR2=ABS(LOX(LN3)-LOX(LN4))
                  IF(L.LT.MIN(ISVAR1,ISVAR2)) CYCLE
                  IF(L.GT.MAX(ISVAR1,ISVAR2)) CYCLE
                  AUX(:)=CHI(:,LN3)*CHI(:,LN4)*POT(:)
                  CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
                  ULITTLE(L+1,LN1,LN2,LN3,LN4)=SVAR
                  ULITTLE(L+1,LN2,LN1,LN3,LN4)=SVAR
                  ULITTLE(L+1,LN1,LN2,LN4,LN3)=SVAR
                  ULITTLE(L+1,LN2,LN1,LN4,LN3)=SVAR
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
                              CALL TRACE$POP()
        RETURN
        END
!
!       ............................................................................
        SUBROUTINE LDAPLUSU_MODULITTLEWITHPARMS(LNX,LOX,LRX,USEUPAR,UPAR &
       &                                       ,USEJPAR,JPAR,MAINLN,ULITTLE)
!       ** CALCULATES THE INTERACTION ENERGY                                      **
!       **                                                                        **
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: LNX
        INTEGER(4),INTENT(IN) :: LOX(LNX)
        INTEGER(4),INTENT(IN) :: LRX
        LOGICAL(4),INTENT(IN) :: USEUPAR
        LOGICAL(4),INTENT(IN) :: USEJPAR
        INTEGER(4),INTENT(IN) :: MAINLN(2)
        REAL(8)   ,INTENT(IN) :: UPAR
        REAL(8)   ,INTENT(IN) :: JPAR
        REAL(8)   ,INTENT(INOUT):: ULITTLE(LRX+1,LNX,LNX,LNX,LNX)
        REAL(8)   ,PARAMETER  :: FIVEEIGTH=0.625D0
        REAL(8)               :: PI,FOURPI
        REAL(8)               :: FAC
        INTEGER(4)            :: L,LN,LNPROBE,N
        REAL(8)               :: RAWJPAR,RAWUPAR
        REAL(8)               :: SVAR
!       ****************************************************************************
        PI=4.D0*DATAN(1.D0)
        FOURPI=4.D0*PI
!
!       =============================================================================
!       == DETERMIN SHELL TO WHICH UPAR AND JPAR BELONG
!       =============================================================================
        LNPROBE=-1
        N=0
        DO LN=1,LNX
          L=LOX(LNX)
          IF(L.NE.MAINLN(1)) CYCLE
          N=N+1
          IF(N.EQ.MAINLN(2))THEN
            LNPROBE=LN  
            EXIT
          END IF
        ENDDO
        IF(LNX.EQ.1) LNPROBE=1
        IF(LNPROBE.EQ.-1) THEN
          CALL ERROR$MSG('LDAPLUSU_MODULITTLEWITHPARMS')
        END IF
!
!       =============================================================================
!       == SCALE UP ULITTLE TO SATISFY UPAR
!       =============================================================================
        IF(USEUPAR) THEN
          RAWUPAR=ULITTLE(1,LNPROBE,LNPROBE,LNPROBE,LNPROBE)/FOURPI
          SVAR=UPAR/RAWUPAR
          ULITTLE=ULITTLE*SVAR
        END IF
!
!       =============================================================================
!       == SCALE UP JPAR OF THE MAIN SHELL AND SET OTHERS TO ZERO
!       =============================================================================
        IF(USEJPAR) THEN
          IF(LOX(LNPROBE).EQ.2) THEN
            RAWJPAR=0.D0
            IF(LRX+1.GT.3) THEN
              SVAR=1.D0/(14.D0*FOURPI)
              RAWJPAR=RAWJPAR+ULITTLE(3,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
            END IF 
            IF(LRX+1.GT.5) THEN
              SVAR=1.D0/(14.D0*FOURPI)
              RAWJPAR=RAWJPAR+ULITTLE(5,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
            END IF
            IF(RAWJPAR.GT.0.D0) THEN
              SVAR=JPAR/RAWJPAR
              ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)=ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
            END IF
          ELSE IF(LOX(LNPROBE).EQ.3) THEN
            RAWJPAR=0.D0
            IF(LRX+1.GT.3) THEN
              SVAR=268.D0/(6435.D0*FOURPI)
              RAWJPAR=RAWJPAR+ULITTLE(3,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
            END IF 
            IF(LRX+1.GT.5) THEN
              SVAR=195.D0/(6435.D0*FOURPI)
              RAWJPAR=RAWJPAR+ULITTLE(5,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
            END IF
            IF(LRX+1.GT.7) THEN
              SVAR=250.D0/(6435.D0*FOURPI)
              RAWJPAR=RAWJPAR+ULITTLE(7,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
            END IF
            IF(RAWJPAR.GT.0.D0) THEN
              SVAR=JPAR/RAWJPAR
              ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)=ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
            END IF
          ELSE
            IF(JPAR.EQ.0) THEN
              ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)=0.D0
            ELSE
              CALL ERROR$MSG('JPAR.NE.0 CAN ONLY BE SET OF D AND F SHELLS')
              CALL ERROR$MSG('LDAPLUSU_MODULITTLEWITHPARMS')
            END IF
          END IF
        END IF
        
!
!       =============================================================================
!       == SET UP ULITTLE
!       =============================================================================
!!$
!!$        ULITTLE=0.D0
!!$        DO L=0,LRX
!!$          FAC=0.D0
!!$          IF(L.EQ.0)  FAC=UPAR
!!$          IF(LRX.GE.2)FAC=JPAR*14.D0/(1.D0+FIVEEIGTH)
!!$          IF(LRX.GE.4)FAC=FIVEEIGTH*JPAR*14.D0/(1.D0+FIVEEIGTH)
!!$          ULITTLE(L+1,:,:,:,:)=FAC*4.D0*PI/REAL(2*L+1,KIND=8)
!!$        ENDDO
        RETURN
        END
!
!       ............................................................................
        SUBROUTINE LDAPLUSU_UTENSOR(LRX,NORB,LNX,LOX,ULITTLE,U)
!       ** CALCULATES THE INTERACTION ENERGY                                      **
!       **                                                                        **
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: LRX
        INTEGER(4),INTENT(IN) :: LNX
        INTEGER(4),INTENT(IN) :: LOX(LNX)
        REAL(8)   ,INTENT(IN) :: ULITTLE(LRX+1,LNX,LNX,LNX,LNX)
        INTEGER(4),INTENT(IN) :: NORB
        REAL(8)   ,INTENT(OUT):: U(NORB,NORB,NORB,NORB)
        INTEGER(4)            :: LN1,LN2,LN3,LN4
        INTEGER(4)            :: IORB1,IORB2,IORB3,IORB4
        INTEGER(4)            :: L1,L2,L3,L4
        INTEGER(4)            :: LM1,LM2,LM3,LM4
        INTEGER(4)            :: M1,M2,M3,M4
        INTEGER(4)            :: L,M,LM,LX
        REAL(8)               :: CG1,CG2
        REAL(8)               :: SVAR
!       ****************************************************************************
                              CALL TRACE$PUSH('LDAPLUSU_UTENSOR')
        IORB1=0
        DO LN1=1,LNX
          L1=LOX(LN1)
          LM1=L1**2
          DO M1=1,2*L1+1
            IORB1=IORB1+1
            LM1=LM1+1
!
            IORB2=0
            DO LN2=1,LNX
              L2=LOX(LN2)
              LM2=L2**2
              DO M2=1,2*L2+1
                IORB2=IORB2+1
                LM2=LM2+1
!
                IORB3=0
                DO LN3=1,LNX
                  L3=LOX(LN3)
                  LM3=L3**2
                  DO M3=1,2*L3+1
                    IORB3=IORB3+1
                    LM3=LM3+1
                    IF(LM3.LT.LM1) CYCLE
!
                    IORB4=0
                    DO LN4=1,LNX
                      L4=LOX(LN4)
                      LM4=L4**2
                      DO M4=1,2*L4+1
                        IORB4=IORB4+1
                        LM4=LM4+1
                        IF(LM4.LT.LM2) CYCLE
!           
                        LX=MIN(LRX,L2+L4,L1+L3)
                        SVAR=0.D0
                        LM=0
                        DO L=0,LX
                          DO M=1,2*L+1
                            LM=LM+1
                            CALL CLEBSCH(LM2,LM4,LM,CG1)
                            CALL CLEBSCH(LM3,LM1,LM,CG2)
                            SVAR=SVAR+CG1*CG2*ULITTLE(L+1,LN2,LN4,LN3,LN1)
                          ENDDO
                        ENDDO
                        U(IORB1,IORB2,IORB3,IORB4)=SVAR
                        U(IORB1,IORB4,IORB3,IORB2)=SVAR
                        U(IORB3,IORB2,IORB1,IORB4)=SVAR
                        U(IORB3,IORB4,IORB1,IORB2)=SVAR
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
                              CALL TRACE$POP()
        RETURN
        END
!
!       ............................................................................
        SUBROUTINE LDAPLUSU_ETOT(ID,L,UPAR,JPAR,NORB,U,RHO,ETOT,HAM)
!       ** CALCULATES THE INTERACTION ENERGY                                      **
!       **                                                                        **
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID     ! ID FOR DOUBLE COUNTING CORRECTION
        INTEGER(4)  ,INTENT(IN) :: L      ! ANGULAR MOMENTUM
        REAL(8)     ,INTENT(IN) :: UPAR   ! U-PARAMETER
        REAL(8)     ,INTENT(IN) :: JPAR   ! J-PARAMETER
        INTEGER(4)  ,INTENT(IN) :: NORB   ! BASIS-SET SIZE              
        REAL(8)     ,INTENT(IN) :: U(NORB,NORB,NORB,NORB) ! U TENSOR
        COMPLEX(8)  ,INTENT(IN) :: RHO(NORB,NORB,2,2) ! DENSITY MATRIX
        REAL(8)     ,INTENT(OUT):: ETOT    ! DOUBLE COUNTINNG ENERGY
        COMPLEX(8)  ,INTENT(OUT):: HAM(NORB,NORB,2,2)  ! DE/D(RHO^*)        
        REAL(8)                 :: E
        REAL(8)                 :: F(2)
        REAL(8)                 :: V(2)
        INTEGER(4)              :: IS,I
!       ****************************************************************************
                              CALL TRACE$PUSH('LDAPLUSU_ETOT')
        ETOT=0.D0
        HAM(:,:,:,:)=(0.D0,0.D0)
!       =============================================================================
!       ==  INTERACTION ENERGY                                                     ==
!       =============================================================================
        CALL LDAPLUSU_INTERACTION(NORB,U,RHO,ETOT,HAM)
!
!       =============================================================================
!       ==  DOUBLE COUNTING CORRECTION                                             ==
!       =============================================================================
        DO IS=1,2
          F(IS)=0.D0
          DO I=1,NORB
            F(IS)=F(IS)+REAL(RHO(I,I,IS,IS))
          ENDDO
        ENDDO
        CALL LDAPLUSU_DOUBLECOUNTING(ID,L,UPAR,JPAR,F,E,V)
        ETOT=ETOT+E
        DO IS=1,2
          DO I=1,NORB
            HAM(I,I,IS,IS)=HAM(I,I,IS,IS)+CMPLX(V(IS),0.D0)
          ENDDO
        ENDDO
                              CALL TRACE$POP()
        RETURN
        END        
!
!       ............................................................................
        SUBROUTINE LDAPLUSU_INTERACTION(NORB,U,RHO,ETOT,HAM)
!       ** CALCULATES THE INTERACTION ENERGY                                      **
!       **                                                                        **
        IMPLICIT NONE
        INTEGER(4)  ,INTENT(IN) :: NORB ! BASIS-SET SIZE              
        REAL(8)     ,INTENT(IN) :: U(NORB,NORB,NORB,NORB) ! U TENSOR
        COMPLEX(8)  ,INTENT(IN) :: RHO(NORB,NORB,2,2) ! DENSITY MATRIX
        REAL(8)     ,INTENT(OUT):: ETOT    ! DOUBLE COUNTINNG ENERGY
        COMPLEX(8)  ,INTENT(OUT):: HAM(NORB,NORB,2,2)  ! DE/D(DENMAT)        
        REAL(8)                 :: UIJKL
        COMPLEX(8)              :: RHOT(NORB,NORB)
        COMPLEX(8)              :: HAMT(NORB,NORB)
        REAL(8)                 :: SVAR,EBLOCK
        INTEGER(4)              :: I,J,K,L,IS1,IS2
!       ****************************************************************************
                              CALL TRACE$PUSH('LDAPLUSU_INTERACTION')
        RHOT(:,:)=RHO(:,:,1,1)+RHO(:,:,2,2)
        ETOT=0.D0
        HAMT(:,:)=(0.D0,0.D0)
        HAM(:,:,:,:)=(0.D0,0.D0)
        DO L=1,NORB
          DO K=1,NORB
            EBLOCK=0.D0
            DO J=1,NORB
              DO I=1,NORB
                UIJKL=0.5D0*U(I,J,K,L)
                SVAR=REAL(RHOT(K,I)*RHOT(L,J))     ! HARTREE
                HAMT(K,I)=HAMT(K,I)+UIJKL*RHOT(L,J)
                HAMT(L,J)=HAMT(L,J)+UIJKL*RHOT(K,I)
                DO IS1=1,2
                  DO IS2=1,2
                    SVAR=SVAR-REAL(RHO(L,I,IS2,IS1)*RHO(K,J,IS1,IS2))  ! EXCHANGE
                    HAM(L,I,IS2,IS1)=HAM(L,I,IS2,IS1)-UIJKL*RHO(K,J,IS1,IS2)
                    HAM(K,J,IS1,IS2)=HAM(K,J,IS1,IS2)-UIJKL*RHO(L,I,IS2,IS1)
                  ENDDO
                ENDDO
                EBLOCK=EBLOCK+SVAR*UIJKL
              ENDDO
            ENDDO
            ETOT=ETOT+EBLOCK
          ENDDO
        ENDDO
        HAM(:,:,1,1)=HAM(:,:,1,1)+HAMT(:,:)
        HAM(:,:,2,2)=HAM(:,:,2,2)+HAMT(:,:)
                              CALL TRACE$POP()
        RETURN
        END
!
!       ............................................................................
        SUBROUTINE LDAPLUSU_DOUBLECOUNTING(ID,L,U,J,F,E,V)
!       ** DOUBLE COUNTING CORRECTION TO THE LDA+U TOTAL ENERGY                   **
!       ** THE ENERGY PRODUCED SHALL BE ADDED TO THE TOTAL ENERGY                 **
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID   ! SWITCH BETWEEN DIFFERENT FORMULATIONS
        INTEGER(4)  ,INTENT(IN) :: L    ! MAIN ANGULAR MOMENTUM
        REAL(8)     ,INTENT(IN) :: U    ! U-PARAMETER
        REAL(8)     ,INTENT(IN) :: J    ! J-PARAMETER
        REAL(8)     ,INTENT(IN) :: F(2) ! MEAN OCCUPATION/PER SPIN
        REAL(8)     ,INTENT(OUT):: E    ! DOUBLE COUNTINNG ENERGY
        REAL(8)     ,INTENT(OUT):: V(2)  ! DE/DF
        REAL(8)                 :: FTOT,VTOT,SVAR
!       ****************************************************************************
                              CALL TRACE$PUSH('LDAPLUSU_DOUBLECOUNTING')
!       =============================================================================
!       ==  OPTION                                                                 ==
!       =============================================================================
        IF(ID.EQ.'FLL') THEN
          FTOT=F(1)+F(2)
          E=0.5D0*U*FTOT*(FTOT-1.D0) &
      &    -0.5D0*J*F(1)*(F(1)-1.D0) &
      &    -0.5D0*J*F(2)*(F(2)-1.D0)
          VTOT=0.5D0*U*(2.D0*FTOT-1.D0)
          V(1)=VTOT-0.5D0*J*(2.D0*F(1)-1.D0)
          V(2)=VTOT-0.5D0*J*(2.D0*F(2)-1.D0)
!
!       =============================================================================
!       ==  OPTION APPROXIMATE MEAN FIELD (AMF)                                    ==
!       =============================================================================
        ELSE IF(ID.EQ.'AMF') THEN
          SVAR=REAL(L)/REAL(2*L+1)*(U-J)
          E=U*F(1)*F(2)+SVAR*(F(1)**2+F(2)**2)
          V(1)=U*F(2)+2.D0*SVAR*F(1)
          V(2)=U*F(1)+2.D0*SVAR*F(2)
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$STOP('LDAPLUSU_DOUBLECOUNTING')
        END IF
!
!       =============================================================================
!       == REVERT SIGN SO THAT THE RESULTING ENERGY ADDS WITH POSITIVE SIGN
!       =============================================================================
        E=-E
        V(:)=-V    
                              CALL TRACE$POP()
        RETURN
        END

!=====================================================================================
!=====================================================================================
!=====================================================================================
!=====================================================================================
!=====================================================================================
!=====================================================================================
!=====================================================================================
!
!       ............................................................................
        SUBROUTINE LDAPLUSU_DOWNFOLD(GID,NR,LNX,LOX,CHI,NORB,NCORR,TCORR,A)
!       ** CALCULATES THE INTERACTION ENERGY                                      **
!       **                                                                        **
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: GID
        INTEGER(4),INTENT(IN) :: NR
        INTEGER(4),INTENT(IN) :: LNX
        INTEGER(4),INTENT(IN) :: LOX(LNX)
        INTEGER(4),INTENT(IN) :: NORB
        INTEGER(4),INTENT(IN) :: NCORR       !#(CORRELATED ORBITALS)
        LOGICAL(4),INTENT(IN) :: TCORR(LNX)
        REAL(8)   ,INTENT(IN) :: CHI(NR,LNX)
        REAL(8)   ,INTENT(OUT):: A(NCORR,NORB)
        REAL(8)               :: XLITTLE(LNX,LNX)
        REAL(8)               :: R(NR)
        REAL(8)               :: AUX(NR)
        INTEGER(4)            :: N,M,L,L1,L2,LN,LN1,LN2
        INTEGER(4)            :: IORB
        INTEGER(4)            :: IA1,IB1,IA2,IB2
        REAL(8)               :: X11(NCORR,NCORR)
        REAL(8)               :: X12(NCORR,NORB-NCORR)
        REAL(8)               :: X22(NORB-NCORR,NORB-NCORR)
        REAL(8)               :: X22INV(NORB-NCORR,NORB-NCORR)
        REAL(8)               :: X12X22INV(NCORR,NORB-NCORR)
        REAL(8)               :: A11(NCORR,NCORR)
        REAL(8)               :: A12(NCORR,NORB-NCORR)
        REAL(8)               :: SVAR
        LOGICAL(4)            :: TCORR1,TCORR2
!       ****************************************************************************
        CALL RADIAL$R(GID,NR,R)
        XLITTLE(:,:)=0.D0
        DO LN1=1,LNX
          L1=LOX(LN1)
          DO LN2=LN1,LNX
            L2=LOX(LN1)
            IF(L1.NE.L2) CYCLE
            AUX(:)=CHI(:,LN1)*CHI(:,LN2)*R(:)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            XLITTLE(LN1,LN2)=SVAR
            XLITTLE(LN2,LN1)=SVAR
          ENDDO
        ENDDO
!       =============================================================================
!       == DETERMINE DIMENSIONS OF THE SUBSPACES OF THE OVERLAP MATRIX             ==
!       =============================================================================
        N=0
        DO LN=1,LNX
          L=LOX(LN)
          IF(TCORR(LN)) N=N+2*L+1
        ENDDO
        IF(N.NE.NCORR) THEN
          CALL ERROR$STOP('LDAPLUSU_DOWNFOLD')
        END IF

!       =============================================================================
!       == DETERMINE OVERLAP MATRICES                                              ==
!       =============================================================================
        IA1=0       ! POINTER FOR LEFT-HAND HEAD FUNCTIONS
        IB1=0       ! POINTER FOR LEFT-HAND TAIL FUNCTIONS
        DO LN1=1,LNX
          L1=LOX(LN1)
          TCORR1=TCORR(LN1)
          IA2=0    ! POINTER FOR RIGHT-HAND HEAD FUNCTIONS
          IB2=0    ! POINTER FOR RIGHT-HAND TAIL FUNCTIONS
          DO LN2=1,LNX
            L2=LOX(LN2)
            TCORR2=TCORR(LN2)
            IF(L1.EQ.L2) THEN
              IF(TCORR(LN1).AND.TCORR(LN2)) THEN
                DO M=1,2*L2+1
                  X11(IA1+M,IA2+M)=XLITTLE(LN1,LN2)
                ENDDO
              ELSE IF(TCORR(LN1).AND.(.NOT.TCORR(LN2))) THEN
                DO M=1,2*L2+1
                  X12(IA1+M,IB2+M)=XLITTLE(LN1,LN2)
                ENDDO
              ELSE IF((.NOT.TCORR(LN1)).AND.(.NOT.TCORR(LN2))) THEN
                DO M=1,2*L2+1
                  X22(IB1+M,IB2+M)=XLITTLE(LN1,LN2)
                ENDDO
              END IF
            END IF
            IF(TCORR2) THEN
              IA2=IA2+2*L2+1
            ELSE
              IB2=IB2+2*L2+1
            END IF
          ENDDO
          IF(TCORR1) THEN
            IA1=IA1+2*L1+1
          ELSE
            IB1=IB1+2*L1+1
          END IF
        ENDDO
!
!       =============================================================================
!       == DETERMINE OVERLAP MATRICES                                              ==
!       =============================================================================
        CALL LIB$INVERTR8(NORB-NCORR,X22,X22INV)
        X12X22INV=MATMUL(X12,X22INV)
        A11=X11-MATMUL(X12X22INV,TRANSPOSE(X12))
        CALL LIB$INVERTR8(NCORR,A11,X11)
        A11=X11
        A12=-MATMUL(A11,X12X22INV)
!
!       =============================================================================
!       == MAP ONTO ALL ORBITALS                                                   ==
!       =============================================================================
        IORB=0
        IA2=0
        IB2=0
        DO LN2=1,LNX
          L2=LOX(LN2)
          TCORR2=TCORR(LN2)
          IF(TCORR2) THEN
            DO M=1,2*L2+1
              IORB=IORB+1
              IA2=IA2+1
              A(:,IORB)=A11(:,IA2)
            ENDDO
          ELSE
            DO M=1,2*L2+1
              IORB=IORB+1
              IB2=IB2+1
              A(:,IORB)=A12(:,IB2)
            ENDDO
          END IF
        ENDDO
        RETURN
        END
