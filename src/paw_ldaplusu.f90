!Todo :
! dath is still real and should probably be complex like denmat

MODULE LDAPLUSU_MODULE
TYPE THISTYPE
LOGICAL(4)             :: TINI=.FALSE.
LOGICAL(4)             :: TON=.FALSE.
INTEGER(4)             :: GID
INTEGER(4)             :: NR
INTEGER(4)             :: LNXCHI
INTEGER(4),POINTER     :: LOXCHI(:)
INTEGER(4),POINTER     :: norb(:)  ! x(# local functions per angular momentum)
INTEGER(4)             :: LRX
INTEGER(4)             :: NCHI
REAL(8)   ,POINTER     :: CHI(:,:)
REAL(8)                :: DIEL=1     ! DIELECTRIC CONSTANT
REAL(8)                :: UPAR=0.D0  ! U-PARAMETER
REAL(8)                :: JPAR=0.D0  ! J-PARAMETER
REAL(8)                :: RCUT=0.D0  ! J-PARAMETER
REAL(8)   ,POINTER     :: ULITTLE(:,:,:,:,:)
REAL(8)   ,POINTER     :: downfold(:,:)
END TYPE THISTYPE
TYPE(THISTYPE),ALLOCATABLE,TARGET :: THISARRAY(:)
TYPE(THISTYPE),POINTER :: THIS
LOGICAL(4)             :: TON=.false.
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
return
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
print*,' potupdown ',potupdown 
          DETOT=DETOT+ELDAUSUMMAND
          PRINT*, '================='
          PRINT*,' total energy ',detot,' in ev: ',DETOT*27.21440D0
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
stop 'forced stop'
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
        INTEGER(4),intent(in) :: NOXA      !NUMBER OF XA
        REAL(8)   ,INTENT(OUT):: OVER(NOXA,LMNXX)
        REAL(8)   ,intent(in) :: RCUT       !2.33D0      !20.0D0  !2.33     
        REAL(8)               :: R(NR)
        INTEGER(4)            :: ALPHA,ALPHA2
        REAL(8)               :: NU(LMNXX,LMNXX,2)
        REAL(8)   ,ALLOCATABLE:: XA(:,:)
        INTEGER(4)            :: BETA,BETA2
        INTEGER(4)            :: SHOAE(LMNXX) ! lm(lmn) for AEPHI
        INTEGER(4),ALLOCATABLE:: SHOXA(:)     ! lm(lmn) for XA        
        REAL(8)               :: WORK(NR),svar
        REAL(8)               :: X(LNX,LNX)
        INTEGER(4)            :: IR
        INTEGER(4)            :: L1,L2
        INTEGER(4)            :: m
        INTEGER(4)            :: LN1,LN2,ln
        INTEGER(4)            :: LM1,LM2,lm
        INTEGER(4)            :: LMN1,LMN2,lmn
        INTEGER(4)            :: SPIN
!       *****************************************************************************
        CALL RADIAL$R(GID,NR,R)
        nue(:,:,:)=0.d0
        nuenn(:,:)=0.d0
        over(:,:)=0.d0
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

         nue(:,:,:)=0.d0
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
              SHOAE(LMN)=LOX(LN)**2+LM     ! (lm)
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
!         == truncate xa beyond rcut
          DO IR=1,NR
            IF (R(IR).GT.RCUT) THEN
              XA(IR,1)=0.D0
            END IF
          END DO

!         == normalize xa
          work(:)=XA(:,1)*XA(:,1)*R(:)**(2)
          CALL RADIAL$INTEGRAL(GID,NR,work,svar)
          XA(:,1)=XA(:,1)/SQRT(svar)
          work(:)=XA(:,1)*XA(:,1)*R(:)**(2)
          CALL RADIAL$INTEGRAL(GID,NR,work,svar)
          PRINT*, 'RADIALINTEGRAL2',svar

!         ==  CALCULATE OVERLAP  <xa|aephi>
          DO LM1=1,NOXA
            LMN1=0
            DO LN1=1,LNX
              DO LM2=1,2*LOX(LN1)+1
                LMN1=LMN1+1
                IF(SHOAE(LMN1).EQ.SHOXA(LM1)) THEN
                  work(:)=AEPHI(:,LN1)*XA(:,LM1)*R(:)**2
                  CALL RADIAL$INTEGRAL(GID,NR,work,svar) 
                  OVER(LM1,LMN1)=svar
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
        PRINT *, 'b=hubbard ennergy=',EU,' in ev',eu*27.2144D0
        
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
print*,'occu ',occu
print*,'dc parms ',l,uparameter,jparameter,occu
          EDCFLL=UPARAMETER/2*(OCCUGES**2-OCCUGES) &
      &         -JPARAMETER/2*(OCCU(1)**2-OCCU(1)) &
      &         -JPARAMETER/2*(OCCU(2)**2-OCCU(2))
!         PRINT*,'E DC OF FLL',EDCFLL*27.2144D0
          EFLLTOT=EU-EDCFLL
          DETOT=EFLLTOT
PRINT *, 'double counting energy FLL',edcfll,' in ev:',edcfll*27.2144D0 

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
!         == pick the first partial wave with the correct l as local orbital
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
        do isp=1,nsp
          allocate(thisarray(isp)%norb(10))
          thisarray(isp)%norb(:)=10
          thisarray(isp)%norb(1:2)=0 ! per default correlate only d and f-shells
        enddo
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
        ELSE IF(ID.EQ.'JPAR') THEN
          THIS%JPAR=VAL
        ELSE IF(ID.EQ.'DIEL') THEN
          THIS%DIEL=VAL
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
          return
        end if
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
        COMPLEX(8),allocatable :: DATH(:,:,:)
        LOGICAL(4)             :: TUFROMPARMS=.FALSE.
        INTEGER(4)             :: GID
        INTEGER(4)             :: NR
        INTEGER(4)             :: LNX
        INTEGER(4),ALLOCATABLE :: LOX(:)
        INTEGER(4)             :: Nchi
        REAL(8)   ,ALLOCATABLE :: PHITOCHI(:,:)
        INTEGER(4)             :: LNXPHI
        INTEGER(4),ALLOCATABLE :: LOXPHI(:)
        INTEGER(4)             :: Nphi
        REAL(8)   ,ALLOCATABLE :: U(:,:,:,:)
        COMPLEX(8),ALLOCATABLE :: RHO(:,:,:,:)
        COMPLEX(8),ALLOCATABLE :: HAM(:,:,:,:)
        COMPLEX(8),allocatable :: matss(:,:,:,:)
        INTEGER(4)             :: IS1,IS2,I,j,ln,m
        real(8)                :: pi,svar
        real(8)                :: upar,jpar,jpar1
!       ****************************************************************************
        if(.not.ton) return
                              call trace$push('LDAPLUSU$ETOT')
        CALL LDAPLUSU$SELECT(ISP_)
        IF(.NOT.THIS%TON) RETURN
        pi=4.d0*datan(1.d0)
        
!       ============================================================================
!       ==  CONSTRUCT LOCAL ORBITALS                                              ==
!       ============================================================================
        IF(.NOT.THIS%TINI) THEN
          CALL LDAPLUSU_CHIFROMPHI()
!!$print*,'tini ',this%tini
!!$print*,'ton ',this%ton
!!$print*,'lnxchi ',this%lnxchi
!!$print*,'nr     ',this%nr
!!$print*,'loxchi ',this%loxchi
!!$print*,'norb   ',this%norb
!!$print*,'nchi   ',this%nchi
!!$print*,'diel   ',this%diel
!!$print*,'upar   ',this%upar
!!$print*,'jpar   ',this%jpar
!!$!print*,'ulittle',this%ulittle

!         ==========================================================================
!         ==  CALCULATE SMALL U-TENSOR                                            ==
!         ==========================================================================
          GID=THIS%GID
          NR=THIS%NR
          LNX=THIS%LNXCHI
          ALLOCATE(LOX(LNX))
          LOX=THIS%LOXCHI
          nchi=THIS%NCHI
          ALLOCATE(THIS%ULITTLE(LRX+1,LNX,LNX,LNX,LNX))
          IF(TUFROMPARMS) THEN
            CALL LDAPLUSU_ULITTLEFROMPARMS(LNX,LOX,LRX,THIS%UPAR,THIS%JPAR,THIS%ULITTLE)
          ELSE
            CALL LDAPLUSU_ULITTLE(GID,NR,LRX,LNX,LOX,THIS%CHI,THIS%ULITTLE)
          END IF
          this%Ulittle=this%Ulittle/THIS%DIEL
SVAR=4.D0*PI*THIS%UPAR/THIS%ULITTLE(1,1,1,1,1)
THIS%ULITTLE=THIS%ULITTLE*SVAR
        ELSE
          GID=THIS%GID
          NR=THIS%NR
          LNX=THIS%LNXCHI
          ALLOCATE(LOX(LNX))
          LOX=THIS%LOXCHI
          Nchi=THIS%NCHI
        END IF
        CALL SETUP$SELECT(ISP)
        CALL SETUP$lnx(isp,LNXPHI)
        ALLOCATE(loxphi(lnxphi))
        CALL SETUP$lofln(isp,lnxphi,LOXPHI)
        nphi=sum(2*loxphi+1)

!
!       ==========================================================================
!       ==  DOWNFOLD                                                            ==
!       ==========================================================================
        allocate(phitochi(nchi,nphi))
        allocate(rho(nchi,nchi,2,2))
        ALLOCATE(ham(Nchi,Nchi,2,2))
        allocate(matss(nphi,nphi,2,2))
        allocate(dath(nphi,nphi,ndimd))
        ALLOCATE(U(Nchi,Nchi,Nchi,Nchi))
!       == TRANSFORM FROM TOTAL SPIN TO SPIN-SPIN
        CALL LDAPLUSU_SPINDENMAT('FORWARD',NDIMD,nphi,DENMAT(1:nphi,1:nphi,:),MATSS)
!       == TRANSFORM FROM PHI TO CHI
        call LDAPLUSU_MAPTOCHI(LNX,LOX,nchi,LNXPHI,LOXPHI,npHI,phitochi)
        DO IS1=1,2
          DO IS2=1,2
            rho(:,:,IS1,IS2)=MATMUL(PHITOCHI,MATMUL(matss(:,:,IS1,IS2),transpose(PHITOCHI)))
          ENDDO
        ENDDO

do is1=1,2
do is2=1,2
if(sum(abs(rho(:,:,is1,is2))).lt.1.d-3) cycle
PRINT*,'===================== RHO FOR SPIN',IS1,IS2,' ======================'
i=0
do ln=1,lnx
do m=1,2*lox(ln)+1
i=i+1
write(*,fmt='(i3,100f8.3)')lox(ln),real(rho(i,:,IS1,IS2))
ENDDO
ENDDO
ENDDO
ENDDO
!
!       ==========================================================================
!       ==  CALCULATE U-TENSOR                                                  ==
!       ==========================================================================
        CALL LDAPLUSU_UTENSOR(LRX,Nchi,LNX,LOX,THIS%ULITTLE,U)
        upar=0.d0
        jpar=0.d0
        jpar1=0.d0
        do i=1,nchi
          do j=1,nchi
            upar=upar+u(i,j,i,j)
            jpar=jpar+u(i,j,i,j)-u(i,j,j,i)
            if(i.ne.j)jpar1=jpar1+u(i,j,j,i)
          enddo
        enddo
        upar=upar/real(nchi)**2
        jpar=upar-jpar/real(nchi*(nchi-1))
        jpar1=jpar1/real(nchi*(nchi-1))*7.d0/5.d0
print*,'uparameter[ev]    ',upar*27.211d0
print*,'jparameter[ev](1) ',jpar*27.211d0
print*,'jparameter[ev](2) ',jpar1*27.211d0

!
!       ============================================================================
!       ==  CALCULATE LDA+U TOTAL ENERGY CONTRIBUTION                             ==
!       ============================================================================
        CALL LDAPLUSU_ETOT(DCTYPE,L,upar,JPAR,Nchi,U,RHO,ETOT,HAM)
!
!       ============================================================================
!       ==  UPFOLD                                                                ==
!       ============================================================================
!       == TRANSFORM FROM CHI TO PHI ===============================================
        DO IS1=1,2
          DO IS2=1,2
            matSS(:,:,IS1,IS2)=MATMUL(TRANSPOSE(PHITOCHI),MATMUL(HAM(:,:,IS1,IS2),PHITOCHI))
          ENDDO
        ENDDO
!
!       == TRANSFORM FROM (SPIN,SPIN) TO (TOTAL,SPIN) ==============================
        CALL LDAPLUSU_SPINDENMAT('BACK',NDIMD,LMNXX,DATH,matss)
!
!       == MAKE REAL (THIS IS A FUDGE TO BE FIXED IN AUGMENTATION!)
        dath_(:,:,:)=(0.d0,0.d0)
        DATH_(:nphi,:nphi,:)=REAL(DATH)
!!$do is1=1,ndimd
!!$PRINT*,'===================== dath FOR SPIN',IS1,IS2,' ======================'
!!$i=0
!!$do ln=1,lnxphi
!!$do m=1,2*loxphi(ln)+1
!!$i=i+1
!!$write(*,fmt='(i3,100f8.3)')loxphi(ln),real(dath(i,:,IS1))
!!$ENDDO
!!$ENDDO
!!$ENDDO
!
!       ============================================================================
!       ==  UNSELECT LDAPLUSU                                                     ==
!       ============================================================================
        deallocate(u)
        deallocate(loxphi)
        deallocate(lox)
        deallocate(matss)
        deallocate(rho)
        deallocate(ham)
        CALL LDAPLUSU$SELECT(0)
                                     call trace$pop()
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
        INTEGER(4),allocatable:: Lox(:)
        real(8)   ,allocatable:: phi(:,:)
        integer(4)            :: lnxchi 
        integer(4),allocatable:: loxchi(:)
        real(8)   ,allocatable:: chi(:,:)
        real(8)   ,allocatable:: a(:,:)
        REAL(8)   ,ALLOCATABLE:: mat(:,:)
        REAL(8)   ,ALLOCATABLE:: matinv(:,:)
        REAL(8)   ,ALLOCATABLE:: R(:)        
        REAL(8)   ,ALLOCATABLE:: g(:)        
        real(8)   ,parameter  :: rcg=0.3d0
        real(8)   ,allocatable:: aux(:)
        real(8)               :: svar1,svar2
        INTEGER(4)            :: IR
        integer(4)            :: nx,n,lx,l,ln,lnchi,ln0
        integer(4)            :: n1,n2,ln1,ln2,l1,l2
!       ****************************************************************************
                              call trace$push('LDAPLUSU_CHIFROMPHI')
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        THIS%GID=GID
        CALL RADIAL$GETI4(GID,'NR',NR)
        THIS%NR=NR
        CALL SETUP$LNX(ISP,LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$LOFLN(ISP,LNX,lox)
        allocate(phi(nr,lnx))
        CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,phi)
!
!       ============================================================================
!       ==  DIVIDE IN HEAD AN TAIL FUNCTIONS                                      ==
!       ============================================================================
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        ALLOCATE(AUX(NR))
!       == COUNT NUMBER OF CHI FUNCTIONS (PRELIMINARY)
        LX=MAXVAL(LOX)
        LNXCHI=0
        DO L=0,LX
          IF(THIS%NORB(L+1).EQ.0) CYCLE
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            LNXCHI=LNXCHI+1
           ENDDO
        ENDDO
!
!       == ORDER ACCORDING TO L ====================================================
        ALLOCATE(LOXCHI(LNXCHI))        
        ALLOCATE(CHI(NR,LNXCHI))        
        ALLOCATE(A(LNXCHI,LNX))        
        A(:,:)=0.D0
        LNCHI=0
        DO L=0,LX
          IF(THIS%NORB(L+1).EQ.0) CYCLE
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            LNCHI=LNCHI+1
            LOXCHI(LNCHI)=LOX(LN)
            CHI(:,LNCHI)=PHI(:,LN)
            A(LNCHI,LN)=1.D0
          ENDDO
        ENDDO
!
!       == MAKE HEAD FUNCTION ANTIBONDING WITH NODE AT RCUT ==========================
        rcut=this%rcut
        L=-1
        DO LN=1,LNXCHI
          IF(LOXCHI(LN).EQ.L) CYCLE
          L=LOXCHI(LN)
          IF(LN+1.GT.LNXCHI) CYCLE
          IF(LOXCHI(LN+1).NE.L) CYCLE
!         == IMPOSE NODE CONDITION======================================================
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN),RCUT,SVAR1)
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN+1),RCUT,SVAR2)
          CHI(:,LN)=CHI(:,LN)*SVAR2-CHI(:,LN+1)*SVAR1
          A(LN,:)=A(LN,:)*SVAR2-A(LN+1,:)*SVAR1
!         == CUT AT RCUT           ===================================================
          DO IR=1,NR
            IF(R(IR).GT.RCUT) CHI(IR,LN)=0.D0
          ENDDO
!         == NORMALIZE HEAD FUNCTION =============================================
          AUX(:)=CHI(:,LN)**2*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          SVAR1=1.D0/SQRT(SVAR1)
          A(LN,:)=A(LN,:)*SVAR1
          CHI(:,LN)=CHI(:,LN)*SVAR1
        END DO          
!
!       == DELOCALIZE TAIL FUNCTIONS ====================================================
        ALLOCATE(G(NR))
        G(:)=EXP(-(R(:)/RCG)**2)
        L=-1
        DO LN=1,LNXCHI
          IF(LOXCHI(LN).NE.L) THEN
            L=LOXCHI(LN)
            LN0=LN
            CYCLE
          END IF   
!         == MINIMIZE CONTRIBUTION NEAR THE CENTER FOR TAIL FUNCTIONS
          DO LN2=LN0,LN-1
            AUX(:)=G(:)*CHI(:,LN)*CHI(:,LN2)*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
            AUX(:)=G(:)*CHI(:,LN2)**2*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
            SVAR1=SVAR1/SVAR2
            CHI(:,LN)=CHI(:,LN)-CHI(:,LN2)*SVAR1
            A(LN,:)=A(LN,:)-A(LN2,:)*SVAR1
          ENDDO
        ENDDO

!
!       === construct transformation matrix from a ========================================
        lx=maxval(loxchi)
        do l=1,lx
!
          nx=0
          do ln=1,lnxchi
            if(loxchi(ln).eq.l) nx=nx+1
          enddo
          if(nx.eq.0) cycle
if(nx.gt.2) then
  call error$msg('make also tail functions antibonding')
  call error$stop('ldaplusu_chifromphi')
end if

          allocate(mat(nx,nx))
          allocate(matinv(nx,nx))
!          
          n1=0
          do ln1=1,lnxchi
            l1=loxchi(ln1)
            if(l1.ne.l) cycle
            n1=n1+1
            n2=0
            do ln2=1,lnx
              l2=lox(ln2)
              if(l2.ne.l) cycle
              n2=n2+1
              mat(n1,n2)=a(ln1,ln2)
            enddo
          enddo
          call lib$invertr8(nx,mat,matinv)
          n1=0
          do ln1=1,lnxchi
            l1=loxchi(ln1)
            if(l1.ne.l) cycle
            n1=n1+1
            n2=0
            do ln2=1,lnx
              l2=lox(ln2)
              if(l2.ne.l) cycle
              n2=n2+1
              a(ln1,ln2)=matinv(n1,n2)
            enddo
          enddo
          deallocate(mat)
          deallocate(matinv)
        enddo
!
!       === remove tail functions                         ======================================
        ln1=0
        l=-1
        do ln=1,lnxchi
          if(l.ne.loxchi(ln)) then
            l=loxchi(ln)
            n=0
            nx=0
            do ln2=ln,lnxchi
              if(loxchi(ln2).ne.l) exit
              nx=nx+1
            enddo
            nx=max(1,nx-1)   ! throw out the highest tail function in any case
            nx=min(nx,this%norb(l+1))
          end if
          n=n+1
          if(n.le.nx) then
            ln1=ln1+1
            loxchi(ln1)=loxchi(ln)
            chi(:,ln1)=chi(:,ln)
            a(ln1,:)=a(ln,:)
          end if
        enddo
        lnxchi=ln1

!print*,'== write chi.dat'
!open(unit=109,file='chi.dat',form='formatted')
!do ir=1,nr
!  write(109,*)r(ir),chi(ir,1:lnxchi)
!enddo
!close(109)
!stop 'forced after writing chi'
!
!       ============================================================================
!       ==  CUT OF SUPPORT FUNCTIONS                                              ==
!       ============================================================================
        this%lnxchi=lnxchi
        allocate(this%loxchi(lnxchi))
        this%loxchi=loxchi(1:lnxchi)
        this%nchi=sum(2*loxchi(1:lnxchi)+1)
        allocate(this%chi(nr,lnxchi))
        this%chi=chi(:,1:lnxchi)
        allocate(this%downfold(lnxchi,lnx))
        this%downfold(:,:)=a(:lnxchi,:)
!
!       ============================================================================
!       ==  CUT OF SUPPORT FUNCTIONS                                              ==
!       ============================================================================
        DEALLOCATE(R)
                              call trace$pop()
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
                              call trace$push('LDAPLUSU_spindenmat')
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
                              call trace$pop()
        RETURN
        END

!
!       ............................................................................
        SUBROUTINE LDAPLUSU_MAPTOCHI(LNXCHI,LOXCHI,LMNXCHI,LNXPHI,LOXPHI,LMNXPHI,downfold)
!       **                                                                        **
        use ldaplusu_module
        IMPLICIT NONE
        INTEGER(4),INTENT(IN)   :: LNXCHI
        INTEGER(4),INTENT(IN)   :: LOXCHI(LNXCHI)
        INTEGER(4),INTENT(IN)   :: LMNXCHI
        INTEGER(4),INTENT(IN)   :: LNXPHI
        INTEGER(4),INTENT(IN)   :: LOXPHI(LNXPHI)
        INTEGER(4),INTENT(IN)   :: LMNXPHI
        REAL(8)   ,intent(out)  :: DOWNFOLD(LMNXCHI,LMNXPHI)
        INTEGER(4)              :: LMN1,LN1,L1,LMN2,LN2,L2,M
!       ****************************************************************************
                              call trace$push('LDAPLUSU_maptochi')
        downfold=0.d0
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
                              call trace$pop()
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
                              call trace$push('LDAPLUSU_ulittle')
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
            lmax=min(lmax,lrx)
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
                              call trace$pop()
        RETURN
        END
!
!       ............................................................................
        SUBROUTINE LDAPLUSU_ULITTLEFROMPARMS(LNX,LOX,LRX,UPAR,JPAR,ULITTLE)
!       ** CALCULATES THE INTERACTION ENERGY                                      **
!       **                                                                        **
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: LNX
        INTEGER(4),INTENT(IN) :: LOX(LNX)
        INTEGER(4),INTENT(IN) :: LRX
        REAL(8)   ,INTENT(IN) :: UPAR
        REAL(8)   ,INTENT(IN) :: JPAR
        REAL(8)   ,INTENT(OUT):: ULITTLE(LRX+1,LNX,LNX,LNX,LNX)
        REAL(8)   ,PARAMETER  :: FIVEEIGTH=0.625D0
        REAL(8)               :: PI
        REAL(8)               :: FAC
        INTEGER(4)            :: L,LN
!       ****************************************************************************
        PI=4.D0*DATAN(1.D0)
        DO LN=1,LNX
          IF(LOX(LN).NE.2) THEN
            CALL ERROR$MSG('THIS ROUTINE IS ONLY IMPLEMENTED FOR THE D-SHELL')
            CALL ERROR$STOP('LDAPLUSU_ULITTLEFROMPARMS')
          END IF
        ENDDO
!
!       =============================================================================
!       == SET UP ULITTLE
!       =============================================================================
        ULITTLE=0.D0
        DO L=0,LRX
          FAC=0.D0
          IF(L.EQ.0)  FAC=UPAR
          IF(LRX.GE.2)FAC=JPAR*14.D0/(1.D0+FIVEEIGTH)
          IF(LRX.GE.4)FAC=FIVEEIGTH*JPAR*14.D0/(1.D0+FIVEEIGTH)
          ULITTLE(L+1,:,:,:,:)=FAC*4.D0*PI/REAL(2*L+1,KIND=8)
        ENDDO
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
        INTEGER(4)            :: L,M,LM,Lx
        REAL(8)               :: CG1,CG2
        REAL(8)               :: SVAR
!       ****************************************************************************
                              call trace$push('LDAPLUSU_utensor')
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
                              call trace$pop()
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
                              call trace$push('LDAPLUSU_etot')
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
                              call trace$pop()
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
                              call trace$push('LDAPLUSU_interaction')
        RHOT(:,:)=RHO(:,:,1,1)+RHO(:,:,2,2)
        ETOT=0.D0
        HAMT(:,:)=(0.D0,0.D0)
        HAM(:,:,:,:)=(0.D0,0.D0)
        DO L=1,NORB
          DO K=1,NORB
            EBLOCK=0.D0
            DO J=1,NORB
              DO I=1,NORB
                UIJKL=0.5d0*U(I,J,K,L)
                SVAR=REAL(RHOT(K,I)*RHOT(L,J))     ! hartree
                HAMT(K,I)=HAMT(K,I)+UIJKL*RHOT(L,J)
                HAMT(L,J)=HAMT(L,J)+UIJKL*RHOT(K,I)
                DO IS1=1,2
                  DO IS2=1,2
                    SVAR=SVAR-REAL(RHO(L,I,IS2,IS1)*RHO(K,J,IS1,IS2))  ! exchange
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
                              call trace$pop()
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
                              call trace$push('LDAPLUSU_DOUBLECOUNTING')
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
                              call trace$pop()
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
