!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$TEST_ATOMLIB$BOUNDSTATE()
!     **************************************************************************
!     **  TEST THE ROUTINE ATOMLIB$BOUNDSTATE AGAINST EXACT RESULT            **
!     **  USES NON-RELATIVISTIC SCHROEDINGER EQUATION                         **
!     **  EXACT RESULT IS THE HYDROGENIC SOLUTION  POT=Z/R                    **
!     **  REMAINING DIFFERENCE MAY BE DUE TIO THE FINITE NUCLEUS              **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)    ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)    ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                :: AEZ     ! ATOMIC NUMBER
      INTEGER(4)             :: NB      ! #(STATES)
      INTEGER(4),ALLOCATABLE :: LOFI(:) ! ANGULAR MOMENTUM
      INTEGER(4),ALLOCATABLE :: SOFI(:) ! SPIN-ORBIT PARAMETER
      INTEGER(4),ALLOCATABLE :: NNOFI(:)! NUMBER OF NODES
      REAL(8)   ,ALLOCATABLE :: EOFI(:) ! ENERGY EIGENVALUES
      REAL(8)                :: RBOX    ! HARD-SPHERE BOUNDARY CONDITION AT RBOX
      INTEGER(4)             :: GID     ! GRID ID
      REAL(8)                :: DMIN    ! INNERMOST SPACING OF RADIAL GRID
      REAL(8)                :: DMAX    ! OUTERMOST SPACING OF RADIAL GRID
      REAL(8)                :: RX      ! OUTERMOST GRID POINT OF RADIAL GRID
      INTEGER(4)             :: NR       ! #(RADIAL GRID POINTS)
      REAL(8)                :: R1,DEX   ! R(I)=R1*EXP(DEX*(I-1))-R1
      REAL(8)                :: SVAR
      REAL(8)   ,ALLOCATABLE :: POT(:)
      REAL(8)   ,ALLOCATABLE :: PHI(:,:)
      REAL(8)   ,ALLOCATABLE :: SPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: AUX(:),AUX1(:)
      REAL(8)   ,ALLOCATABLE :: DREL(:)
      LOGICAL(4)             :: TREL,TZORA
      INTEGER(4)             :: IB,N,L,ISO
      REAL(8)                :: JPHALF ! J+1/2
      REAL(8)                :: C      ! SPEED OF LIGHT
      REAL(8)                :: ALPHA  ! FINE STRUCTURE CONSTANT
      LOGICAL(4)             :: TFINITENUCLEUS
!     **************************************************************************
      CALL CONSTANTS('C',C)
      ALPHA=1.D0/C
!
!     ==========================================================================
!     == THE FOLLOWING PARAMETERS MUST BE SET BY HAND                         ==
!     ==========================================================================
      AEZ=92.D0   ! URANIUM
      TFINITENUCLEUS=.TRUE.
      NB=38 !50
      DMIN=1.D-6
      DMAX=1.D-1      
      RX=25.D0
!
!     ==========================================================================
!     == DETERMINE RADIAL GRID                                                ==
!     ==========================================================================
      CALL RADIAL$GRIDPARAMETERS(DMIN,DMAX,RX,R1,DEX,NR)
      CALL RADIAL$NEW('SHLOG',GID)
      CALL RADIAL$SETI4(GID,'NR',NR)
      CALL RADIAL$SETR8(GID,'DEX',DEX)
      CALL RADIAL$SETR8(GID,'R1',R1)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == DETERMINE POTENTIAL                                                  ==
!     ==========================================================================
      ALLOCATE(POT(NR))
      IF(TFINITENUCLEUS) THEN
        CALL RADIAL$NUCPOT(GID,NR,AEZ,POT)
      ELSE  ! USE THIS TO FREDUCE THE SIZE OF THE NUCLEUS 
        SVAR=1.D-3
        CALL RADIAL$NUCPOT(GID,NR,SVAR,POT)
        POT=POT*AEZ/SVAR
      END IF

      WRITE(*,*)'ATOMIC NUMBER          ',AEZ
      IF(TFINITENUCLEUS) THEN
        WRITE(*,*)'CORRECT NUCLEAR RADIUS'
      ELSE
        WRITE(*,*)'REDUCED NUCLEAR RADIUS'
      END IF
      WRITE(*,*)'INNERMOST GRID SPACING ',DMIN
      WRITE(*,*)'OUTERMOST GRID SPACING ',DMAX
      WRITE(*,*)'OUTERMOST GRID POINT   ',RX
!
!     ==========================================================================
!     == DETERMINE QUANTUM NUMBERS                                            ==
!     ==========================================================================
      ALLOCATE(LOFI(NB))
      ALLOCATE(SOFI(NB))
      ALLOCATE(NNOFI(NB))
      IB=0
      DO N=1,10
        DO L=0,N-1
          IB=IB+1
          IF(IB.GT.NB) EXIT
          LOFI(IB)=L
          NNOFI(IB)=N-L-1
          SOFI(IB)=0
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE WAVE FUNCTIONS                                             ==
!     ==========================================================================
      ALLOCATE(EOFI(NB))
      EOFI(:)=0.D0
      ALLOCATE(PHI(NR,NB))
      ALLOCATE(SPHI(NR,NB))
      RBOX=R(NR-3)
      CALL ATOMLIB$BOUNDSTATES(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI,POT &
     &                            ,.FALSE.,.FALSE.,RBOX,PHI,SPHI)
!
!     ==========================================================================
!     == WRITE RESULT                                                         ==
!     ==========================================================================
      WRITE(*,*)'TEST NON-RELATIVSISTIC BOUND STATE CALCULATION'
      WRITE(*,FMT='(4A5,4A20)')'IB','L','ISO','NN','PAW','ANALYTIC','PAW-ANALYTIC'
      DO IB=1,NB
        SVAR=-0.5D0*(AEZ/REAL(NNOFI(IB)+LOFI(IB)+1))**2
        WRITE(*,FMT='(4I5,4F20.5)')IB,LOFI(IB),SOFI(IB),NNOFI(IB),EOFI(IB),SVAR,EOFI(IB)-SVAR
      ENDDO
!
!     ==========================================================================
!     == PERFORM SCALAR RELATIVISTIC ZORA CALCULATION                         ==
!     ==========================================================================
      TREL=.TRUE.
      TZORA=.TRUE.
      CALL ATOMLIB$BOUNDSTATES(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI,POT &
     &                            ,TREL,TZORA,RBOX,PHI,SPHI)
!
!     ==========================================================================
!     == WRITE RESULT                                                         ==
!     ==========================================================================
      WRITE(*,*)'TEST SCALAR RELATIVSISTIC (ZORA) BOUND STATE CALCULATION'
      DO IB=1,NB
        SVAR=-0.5D0*(AEZ/REAL(NNOFI(IB)+LOFI(IB)+1))**2
        WRITE(*,FMT='(4I5,4F20.5)')IB,LOFI(IB),SOFI(IB),NNOFI(IB),EOFI(IB),SVAR,EOFI(IB)-SVAR
      ENDDO
!
!     ==========================================================================
!     == PERFORM SCALAR RELATIVISTIC (NON-ZORA) CALCULATION                   ==
!     ==========================================================================
      TREL=.TRUE.
      TZORA=.FALSE.
      CALL ATOMLIB$BOUNDSTATES(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI,POT &
     &                            ,TREL,TZORA,RBOX,PHI,SPHI)
!
!     ==========================================================================
!     == WRITE RESULT                                                         ==
!     ==========================================================================
      WRITE(*,*)'TEST SCALAR RELATIVISTIC BOUND STATE CALCULATION'
      DO IB=1,NB
        SVAR=-0.5D0*(AEZ/REAL(NNOFI(IB)+LOFI(IB)+1))**2
        WRITE(*,FMT='(4I5,4F20.5)')IB,LOFI(IB),SOFI(IB) &
     &                            ,NNOFI(IB),EOFI(IB),SVAR,EOFI(IB)-SVAR
      ENDDO
!
!     ==========================================================================
!     == PREPARE FOR SPIN-ORBIT CALCULATION                                   ==
!     ==========================================================================
      IB=0
      DO N=1,10
        DO L=0,N-1
          DO ISO=-1,1,2  !SIGN OF LS
            IF(L.EQ.0.AND.ISO.EQ.-1) CYCLE
            IB=IB+1
            IF(IB.GT.NB) EXIT
            LOFI(IB)=L
            NNOFI(IB)=N-L-1
            SOFI(IB)=ISO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == PERFORM ZORA CALCULATION WITH SPIN-ORBIT COUPLING                    ==
!     ==========================================================================
      TREL=.TRUE.
      TZORA=.TRUE.
      CALL ATOMLIB$BOUNDSTATES(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI,POT &
     &                            ,TREL,TZORA,RBOX,PHI,SPHI)
!
!     ==========================================================================
!     == WRITE RESULT                                                         ==
!     ==========================================================================
      WRITE(*,*)'TEST ZORA BOUND STATE CALCULATION WITH SPIN-ORBIT'
      WRITE(*,FMT='(4A5,4A20)')'IB','L','ISO','NN','PAW','ANALYTIC','PAW-ANALYTIC'
      DO IB=1,NB
        IF(SOFI(IB).EQ.1) THEN
          JPHALF=REAL(LOFI(IB)+1,KIND=8)
        ELSE
          JPHALF=REAL(LOFI(IB),KIND=8)
        END IF
        N=NNOFI(IB)+LOFI(IB)+1
        SVAR=ALPHA*AEZ/(REAL(N)-JPHALF+SQRT(JPHALF**2-(ALPHA*AEZ)**2))
        SVAR=1.D0/SQRT(1.D0+SVAR**2)-1.D0
        SVAR=SVAR*C**2
        SVAR=2.D0*C**2*SVAR/(2.D0*C**2+SVAR)
        WRITE(*,FMT='(4I5,4F20.5)')IB,LOFI(IB),SOFI(IB),NNOFI(IB),EOFI(IB),SVAR,EOFI(IB)-SVAR
      ENDDO
!
!     ==========================================================================
!     == PERFORM FULLY  RELATIVISTIC CALCULATION                              ==
!     ==========================================================================
      TREL=.TRUE.
      TZORA=.FALSE.
      CALL ATOMLIB$BOUNDSTATES(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI,POT &
     &                            ,TREL,TZORA,RBOX,PHI,SPHI)
!
!     ==========================================================================
!     == WRITE RESULT                                                         ==
!     == EXACT RESULT FROM EQ. 4.46                                           ==
!     ==  (E. VAN LENTHE, THESIS, VRIJE U. AMSTERDAM 1996)                    ==
!     ==========================================================================
      WRITE(*,*)'TEST DIRAC BOUND STATE CALCULATION'
      WRITE(*,FMT='(4A5,4A20)')'IB','L','ISO','NN','PAW','ANALYTIC','PAW-ANALYTIC'
      DO IB=1,NB
        IF(SOFI(IB).EQ.1) THEN
          JPHALF=REAL(LOFI(IB)+1,KIND=8)
        ELSE
          JPHALF=REAL(LOFI(IB),KIND=8)
        END IF
        N=NNOFI(IB)+LOFI(IB)+1
        SVAR=ALPHA*AEZ/(REAL(N)-JPHALF+SQRT(JPHALF**2-(ALPHA*AEZ)**2))
        SVAR=1.D0/SQRT(1.D0+SVAR**2)-1.D0
        SVAR=SVAR*C**2
        WRITE(*,FMT='(4I5,4F20.5)')IB,LOFI(IB),SOFI(IB),NNOFI(IB),EOFI(IB),SVAR,EOFI(IB)-SVAR
      ENDDO
      CALL ATOMLIB_WRITEPHI('NUMDIRACBIG.DAT',GID,NR,NB,PHI)
      CALL ATOMLIB_WRITEPHI('NUMDIRACSMALL.DAT',GID,NR,NB,SPHI)
!
!     == ESTIMATE ERROR OF THE FINITE NUCLEAR SIZE ON THE 1S ORBITAL ===========
      ALLOCATE(AUX(NR))
      ALLOCATE(AUX1(NR))
      AUX(:)=(R(:)**2*POT(:)*Y0+AEZ*R(:))*(PHI(:,1)**2+SPHI(:,1))
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
      WRITE(*,FMT='("CORRECTION DOE TO FINITE NUCLEAR SIZE               ",F20.5)')SVAR
!
      CALL SCHROEDINGER$DREL(GID,NR,POT,0.D0,AUX)
      AUX(:)=-1.D0/(1.D0+2.D0*C**2*R(:)/AEZ)-AUX
      AUX(:)=R(:)**2*0.5D0*AUX(:)*(EOFI(1)-POT(:)*Y0)*(PHI(:,1)**2+SPHI(:,1))
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
      WRITE(*,FMT='("KINETIC ENERGY CORRECTION DOE TO FINITE NUCLEAR SIZE ",F20.5)')SVAR
!
!     ==========================================================================
!     == OBTAIN ANALYTICAL WAVE FUNCTIONS                                     ==
!     ==========================================================================
      CALL SCHROEDINGER$HYDROGENICDIRAC(GID,NR,AEZ,NB,LOFI,NNOFI,SOFI &
     &                                       ,EOFI,PHI,SPHI)
      WRITE(*,*)'ANALYTIC RELATIVISTIC RESULT WITH POINT CHARGE'
      DO IB=1,NB
        WRITE(*,FMT='(4I5,4F20.5)')IB,LOFI(IB),SOFI(IB),NNOFI(IB),EOFI(IB)
      ENDDO
      CALL ATOMLIB_WRITEPHI('DIRACBIG.DAT',GID,NR,NB,PHI)
      CALL ATOMLIB_WRITEPHI('DIRACSMALL.DAT',GID,NR,NB,SPHI)
!
!     ==========================================================================
!     == DETERMINE SMALL COMPONENT FROM ANALYTICAL LARGE COMPONENT TO TEST =====
!     == THE ACCURACY OF THE SMALL COMPONENT CALCULATION                   =====
!     ==========================================================================
      ALLOCATE(DREL(NR))
      DO IB=1,NB
         AUX(:)=0.D0
        CALL SCHROEDINGER$DREL(GID,NR,POT,EOFI(IB),DREL)
        CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,LOFI(IB),SOFI(IB) &
     &                                     ,DREL,AUX,PHI(:,IB),SPHI(:,IB))
      ENDDO
!     == COMPARE DIRACSMALLNUM.DAT WITH DIRACSMALL.DAT =========================
      CALL ATOMLIB_WRITEPHI('DIRACSMALLNUM.DAT',GID,NR,NB,SPHI)

!
      CALL ERROR$MSG('FORCED STOP AFTER TESTING ROUTINE')
      CALL ERROR$STOP('ATOMLIB$TEST_ATOMLIB$BOUNDSTATE')
      RETURN
      END SUBROUTINE ATOMLIB$TEST_ATOMLIB$BOUNDSTATE
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE RADIALPAW_MODULE
TYPE VPAW_TYPE
  LOGICAL            :: TON=.FALSE.
  INTEGER(4)         :: GID
  INTEGER(4)         :: NR
  INTEGER(4)         :: LNX
  INTEGER(4),POINTER :: LOX(:)
  REAL(8)   ,POINTER :: DH(:,:)
  REAL(8)   ,POINTER :: PSPOT(:)
  REAL(8)   ,POINTER :: DO(:,:)
  REAL(8)   ,POINTER :: PRO(:,:)
END TYPE VPAW_TYPE
END MODULE RADIALPAW_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIALPAW$MAKEVPAW(GID,NR,LNX,LOX,PSPOT,DH,DO,PRO,VPAW)
!     **************************************************************************
!     ** SETS UP A THE INPUT FOR CALCULATING THE FOCK POTENTIAL               **
!     **                                                                      **
!     ** IT IS POSSIBLE TO SUBTRACT ALSO A LOCAL POTENTIAL VFOCK%MUX(NR)      **
!     **   AND TO SCALE THE RESULT VFOCK%SCALE.                               **
!     **   THIS CALL INITIALIZES SCALE=1. AND MUX(:)=0.                       **
!     **************************************************************************
      USE RADIALPAW_MODULE, ONLY : VPAW_TYPE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: PSPOT(NR)
      REAL(8)   ,INTENT(IN) :: DH(LNX,LNX)
      REAL(8)   ,INTENT(IN) :: DO(LNX,LNX)
      REAL(8)   ,INTENT(IN) :: PRO(NR,LNX)
      TYPE(VPAW_TYPE),INTENT(INOUT) :: VPAW
!     **************************************************************************
      IF(VPAW%TON)CALL RADIALPAW$CLEANVPAW(VPAW)
      VPAW%TON=.TRUE.
      VPAW%GID=GID
      VPAW%NR=NR
      VPAW%LNX=LNX
      ALLOCATE(VPAW%LOX(LNX))
      ALLOCATE(VPAW%DH(LNX,LNX))
      ALLOCATE(VPAW%DO(LNX,LNX))
      ALLOCATE(VPAW%PRO(NR,LNX))
      ALLOCATE(VPAW%PSPOT(NR))
      VPAW%LOX(:)=LOX(:)
      VPAW%DH(:,:)=DH(:,:)
      VPAW%DO(:,:)=DO(:,:)
      VPAW%PRO(:,:)=PRO(:,:)
      VPAW%PSPOT(:)=PSPOT(:)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIALPAW$CLEANVPAW(VPAW)
!     **************************************************************************
!     ** PREPARES THE COEFFICIENT ARRAY USED FOR DETERMINING THE 
!     **************************************************************************
      USE RADIALPAW_MODULE, ONLY : VPAW_TYPE
      IMPLICIT NONE
      TYPE(VPAW_TYPE),INTENT(INOUT) :: VPAW
!     **************************************************************************
      IF(VPAW%LNX.EQ.-1) RETURN
      VPAW%TON=.FALSE.
      VPAW%GID=-1
      VPAW%NR=0
      VPAW%LNX=0
      DEALLOCATE(VPAW%DH)
      DEALLOCATE(VPAW%DO)
      DEALLOCATE(VPAW%PRO)
      DEALLOCATE(VPAW%PSPOT)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIALPAW$VPSI(GID,NR,VPAW,L,F,VF,OF)
!     **************************************************************************
!     **  APPLIES THE FOCK-OPERATOR TO A SET OF FUNCTIONS                     **
!     **    |G_I>= V_FOCK |F_I>
!     **************************************************************************
      USE RADIALPAW_MODULE, ONLY : VPAW_TYPE
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(IN) :: GID
      INTEGER(4)      ,INTENT(IN) :: NR
      TYPE(VPAW_TYPE) ,INTENT(IN) :: VPAW
      INTEGER(4)      ,INTENT(IN) :: L
      REAL(8)         ,INTENT(IN) :: F(NR)
      REAL(8)         ,INTENT(OUT):: VF(NR)
      REAL(8)         ,INTENT(OUT):: OF(NR)
      REAL(8)         ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)         ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)         ,ALLOCATABLE:: PROJ(:)
      REAL(8)         ,ALLOCATABLE:: CPRO(:)
      REAL(8)         ,ALLOCATABLE:: DPRO(:)
      INTEGER(4)                  :: NPRO
      INTEGER(4)                  :: LNX
      INTEGER(4)                  :: LN,LN1,LN2
      INTEGER(4)                  :: IPRO,IPRO1,IPRO2
      REAL(8)                     :: AUX(NR)
      REAL(8)                     :: R(NR)
!     **************************************************************************
      LNX=VPAW%LNX
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == DETERMINE NUMBER OF RELEVANT PROJECTORS AND ALLOCATE ARRAYS          ==
!     ==========================================================================
      NPRO=0
      DO LN=1,LNX
        IF(VPAW%LOX(LN).NE.L) CYCLE
        NPRO=NPRO+1
      ENDDO
      ALLOCATE(PROJ(NPRO))
      ALLOCATE(CPRO(NPRO))
      ALLOCATE(DPRO(NPRO))
!
!     ==========================================================================
!     == DETERMINE PROJECTIONS PROJ=<PRO|F>                                   ==
!     ==========================================================================
      PROJ(:)=0.D0
      IPRO=0
      DO LN=1,LNX
        IF(VPAW%LOX(LN).NE.L) CYCLE
        IPRO=IPRO+1
        AUX(:)=R(:)**2*VPAW%PRO(:,LN)*F(:)
        CALL RADIAL$INTEGRAL(GID,NR,AUX,PROJ(IPRO))
      ENDDO
!
!     ==========================================================================
!     == DETERMINE CPRO=DH<PRO|F> AND DPRO=DO<PRO|F>                          ==
!     ==========================================================================
      CPRO=0.D0
      DPRO=0.D0
      IPRO1=0
      DO LN1=1,LNX
        IF(VPAW%LOX(LN1).NE.L) CYCLE
        IPRO1=IPRO1+1
        IPRO2=0
        DO LN2=1,LNX
          IF(VPAW%LOX(LN2).NE.L) CYCLE
          IPRO2=IPRO2+1
          CPRO(IPRO1)=CPRO(IPRO1)+VPAW%DH(LN1,LN2)*PROJ(IPRO2)
          DPRO(IPRO1)=DPRO(IPRO1)+VPAW%DO(LN1,LN2)*PROJ(IPRO2)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE VF AND OF                                                  ==
!     ==========================================================================
      VF(:)=VPAW%PSPOT(:)*Y0*F(:)
      OF(:)=F(:)
      IPRO=0
      DO LN=1,LNX
        IF(VPAW%LOX(LN).NE.L) CYCLE
        IPRO=IPRO+1
        VF(:)=VF(:)+VPAW%PRO(:,LN)*CPRO(IPRO)
        OF(:)=OF(:)+VPAW%PRO(:,LN)*DPRO(IPRO)
      ENDDO
!
!     ==========================================================================
!     == CLOSE DOWN                                                           ==
!     ==========================================================================
      DEALLOCATE(PROJ)
      DEALLOCATE(CPRO)
      DEALLOCATE(DPRO)
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE RADIALFOCK_MODULE
TYPE VFOCK_TYPE
  LOGICAL            :: TON=.FALSE.
  INTEGER(4)         :: GID
  INTEGER(4)         :: NR
  REAL(8)            :: SCALE
  INTEGER(4)         :: LRHOX
  INTEGER(4)         :: NB=-1
  INTEGER(4),POINTER :: L(:)
  INTEGER(4),POINTER :: SO(:)
  REAL(8)   ,POINTER :: F(:)
  REAL(8)   ,POINTER :: PSI(:,:)
  REAL(8)   ,POINTER :: MUX(:)
END TYPE VFOCK_TYPE
INTEGER(4)               :: LXCOEFF=-1
REAL(8)     ,ALLOCATABLE :: COEFF(:,:,:)
END MODULE RADIALFOCK_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIALFOCK$PRINTVFOCK(FILE,VFOCK)
!     **************************************************************************
!     ** SETS UP A THE INPUT FOR CALCULATING THE FOCK POTENTIAL               **
!     **                                                                      **
!     ** IT IS POSSIBLE TO SUBTRACT ALSO A LOCAL POTENTIAL VFOCK%MUX(NR)      **
!     **   AND TO SCALE THE RESULT VFOCK%SCALE.                               **
!     **   THIS CALL INITIALIZES SCALE=1. AND MUX(:)=0.                       **
!     **************************************************************************
      USE RADIALFOCK_MODULE
      IMPLICIT NONE
      CHARACTER(*)                :: FILE
      TYPE(VFOCK_TYPE),INTENT(IN) :: VFOCK
      INTEGER(4)                  :: NR
      INTEGER(4)                  :: GID
      REAL(8)         ,ALLOCATABLE:: R(:)
      INTEGER(4)                  :: IR
!     **************************************************************************
      IF(.NOT.VFOCK%TON) RETURN
      GID=VFOCK%GID
      NR=VFOCK%NR
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      REWIND(100)
      DO IR=1,NR
        WRITE(100,*)R(IR),VFOCK%MUX(IR),VFOCK%PSI(IR,:)
      ENDDO
      CLOSE(100)
      DEALLOCATE(R)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIALFOCK$MAKEVFOCK(GID,NR,NB,LOFI,SOFI,FOFI,PSI &
     &                               ,RBOX,LRHOX,VFOCK)
!     **************************************************************************
!     ** SETS UP A THE INPUT FOR CALCULATING THE FOCK POTENTIAL               **
!     **                                                                      **
!     ** IT IS POSSIBLE TO SUBTRACT ALSO A LOCAL POTENTIAL VFOCK%MUX(NR)      **
!     **   AND TO SCALE THE RESULT VFOCK%SCALE.                               **
!     **   THIS CALL INITIALIZES SCALE=1. AND MUX(:)=0.                       **
!     **************************************************************************
      USE RADIALFOCK_MODULE
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(IN) :: GID
      INTEGER(4)      ,INTENT(IN) :: NR
      INTEGER(4)      ,INTENT(IN) :: NB
      INTEGER(4)      ,INTENT(IN) :: LOFI(NB)
      INTEGER(4)      ,INTENT(IN) :: SOFI(NB)
      REAL(8)         ,INTENT(IN) :: FOFI(NB)
      INTEGER(4)      ,INTENT(IN) :: LRHOX
      REAL(8)         ,INTENT(IN) :: RBOX
      REAL(8)         ,INTENT(IN) :: PSI(NR,NB)
      TYPE(VFOCK_TYPE),INTENT(INOUT) :: VFOCK
      REAL(8)                     :: R(NR)
      REAL(8)                     :: AUX(NR),SVAR
      INTEGER(4)                  :: IR,IB
!     **************************************************************************
      IF(VFOCK%TON)CALL RADIALFOCK$CLEANVFOCK(VFOCK)
      CALL RADIAL$R(GID,NR,R)
      VFOCK%TON=.TRUE.
      VFOCK%SCALE=1.D0
      VFOCK%LRHOX=LRHOX
      VFOCK%GID=GID
      VFOCK%NB=NB
      VFOCK%NR=NR
      ALLOCATE(VFOCK%L(NB))
      ALLOCATE(VFOCK%SO(NB))
      ALLOCATE(VFOCK%F(NB))
      ALLOCATE(VFOCK%PSI(NR,NB))
      ALLOCATE(VFOCK%MUX(NR))
      VFOCK%L(:)=LOFI(:)
      VFOCK%SO(:)=SOFI(:)
      VFOCK%F(:)=FOFI(:)
      VFOCK%PSI(:,:)=PSI(:,:)
      VFOCK%MUX(:)=0.D0
! 
!     ==========================================================================
!     == CUT OF TAILS                                                         ==
!     ==========================================================================
      DO IR=1,NR
        IF(R(IR).GE.RBOX) THEN
          VFOCK%PSI(IR:,:)=0.D0
          EXIT
        END IF
      ENDDO
! 
!     ==========================================================================
!     == RE-NORMALIZE                                                         ==
!     ==========================================================================
      DO IB=1,NB
        AUX(:)=(R(:)*VFOCK%PSI(:,IB))**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        VFOCK%PSI(:,IB)=VFOCK%PSI(:,IB)/SQRT(SVAR)
      ENDDO
      RETURN
      END SUBROUTINE RADIALFOCK$MAKEVFOCK
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIALFOCK$CLEANVFOCK(VFOCK)
!     **************************************************************************
!     ** PREPARES THE COEFFICIENT ARRAY USED FOR DETERMINING THE 
!     **************************************************************************
      USE RADIALFOCK_MODULE
      IMPLICIT NONE
      TYPE(VFOCK_TYPE),INTENT(INOUT) :: VFOCK
!     **************************************************************************
      IF(VFOCK%NB.EQ.-1) RETURN
      VFOCK%TON=.FALSE.
      VFOCK%SCALE=0.D0
      VFOCK%NB=-1
      VFOCK%GID=-1
      VFOCK%NR=-1
      VFOCK%LRHOX=-1
      IF(ASSOCIATED(VFOCK%L))DEALLOCATE(VFOCK%L)
      IF(ASSOCIATED(VFOCK%SO))DEALLOCATE(VFOCK%SO)
      IF(ASSOCIATED(VFOCK%F))DEALLOCATE(VFOCK%F)
      IF(ASSOCIATED(VFOCK%PSI))DEALLOCATE(VFOCK%PSI)
      IF(ASSOCIATED(VFOCK%MUX))DEALLOCATE(VFOCK%MUX)
      RETURN
      END SUBROUTINE RADIALFOCK$CLEANVFOCK
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIALFOCK_UPDATECOEFF(LXNEW)
!     **************************************************************************
!     ** PREPARES THE COEFFICIENT ARRAY USED FOR DETERMINING THE 
!     **************************************************************************
      USE RADIALFOCK_MODULE
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(IN) :: LXNEW
      REAL(8)         ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)                  :: LF,LB,LRHO
      INTEGER(4)                  :: LMF,LMB,LMRHO
      INTEGER(4)                  :: IMB,IMRHO
      INTEGER(4)                  :: LX
      REAL(8)                     :: SVAR
      REAL(8)                     :: CG
!     **************************************************************************
      IF(LXNEW.LE.LXCOEFF) RETURN
!
      IF(LXCOEFF.NE.-1) DEALLOCATE(COEFF)
      LXCOEFF=LXNEW
      LX=LXNEW
      ALLOCATE(COEFF(LX+1,2*LX+1,LX+1))
      COEFF(:,:,:)=0.D0
      DO LF=0,LX
        DO LB=0,LX
          DO LRHO=0,2*LX
!          IF(MOD(LF+LB+LRHO,2).EQ.1) CYCLE !EMPIRICAL
            SVAR=0.D0
            LMF=LF**2+1  ! RESULT INDEPEMDENT OF IMF: PIC ANY LMF
            LMB=LB**2
            DO IMB=1,2*LB+1
              LMB=LMB+1
              LMRHO=LRHO**2
              DO IMRHO=1,2*LRHO+1
                LMRHO=LMRHO+1
                CALL SPHERICAL$GAUNT(LMF,LMB,LMRHO,CG)
                SVAR=SVAR+CG**2
              ENDDO
            ENDDO
            SVAR=4.D0*PI*SVAR/REAL((2*LB+1)*(2*LRHO+1),KIND=8) !AVERAGE OVER M_J
            COEFF(LB+1,LRHO+1,LF+1)=SVAR
          ENDDO
        ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE RADIALFOCK_UPDATECOEFF
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIALFOCK$VPSI(GID,NR,VFOCK,L,F,G)
!     **************************************************************************
!     **  APPLIES THE FOCK-OPERATOR TO A SET OF FUNCTIONS                     **
!     **    |G_I>= V_FOCK |F_I>
!     **************************************************************************
      USE RADIALFOCK_MODULE, ONLY : VFOCK_TYPE,COEFF
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(IN) :: GID
      INTEGER(4)      ,INTENT(IN) :: NR
      TYPE(VFOCK_TYPE),INTENT(IN) :: VFOCK
      INTEGER(4)      ,INTENT(IN) :: L
      REAL(8)         ,INTENT(IN) :: F(NR)
      REAL(8)         ,INTENT(OUT):: G(NR)
      REAL(8)         ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)         ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                     :: R(NR)
      REAL(8)                     :: SVAR
      INTEGER(4)                  :: LXB
      INTEGER(4)                  :: LX
      INTEGER(4)                  :: IB
      INTEGER(4)                  :: LB,LRHO
      REAL(8)                     :: AUX1(NR),AUX2(NR)
!     **************************************************************************
      G=0.D0
      IF(.NOT.VFOCK%TON.OR.VFOCK%SCALE.EQ.0.D0) RETURN
!
      LXB=MAXVAL(VFOCK%L)
      LX=MAX(LXB,L,INT((VFOCK%LRHOX+1)/2))
      CALL RADIALFOCK_UPDATECOEFF(LX)
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == APPLY FOCK POTENTIAL                                                 ==
!     ==========================================================================
      G(:)=0.D0
      DO IB=1,VFOCK%NB
        IF(VFOCK%F(IB).LT.1.D-3) CYCLE
        LB=VFOCK%L(IB)
        AUX1(:)=VFOCK%PSI(:,IB)*F(:)
        DO LRHO=0,VFOCK%LRHOX         
          IF(COEFF(LB+1,LRHO+1,L+1).LT.1.D-8) CYCLE    !COEFF.GE.0
          CALL RADIAL$POISSON(GID,NR,LRHO,AUX1,AUX2)
!         == FACTOR 0.5 CANCELS SPIN MULTIPLICITY ============================
          SVAR=-0.5D0*VFOCK%F(IB)*COEFF(LB+1,LRHO+1,L+1) &
     &                           *REAL((2*LRHO+1),KIND=8)/(4.D0*PI)
          G(:)=G(:)+SVAR*VFOCK%PSI(:,IB)*AUX2(:)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == G=SCALE*(VFOCK-MUX*F)                                                **
!     ==========================================================================
      G(:)=G(:)-VFOCK%MUX(:)*Y0*F(:)
      G(:)=VFOCK%SCALE*G(:)
      RETURN
      END SUBROUTINE RADIALFOCK$VPSI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$AESCF(GID,NR,KEY,RBOX,AEZ,NX,NB,LOFI,SOFI,FOFI,NNOFI &
     &                        ,ETOT,POT,VFOCK,EOFI,PHI,SPHI)
!     **************************************************************************
!     ** MAKES A SELF-CONSISTENT CALCULATION OF AN ATOM IN A BOX WITH         **
!     ** RADIUS RBOX (RADIUS IS LIMITED BY THE GRID RBOX<R(NR-3) )            **
!     **                                                                      **
!     ** KEY='START,REL,SO,ZORA,FOCK='                                        **
!     **                                                                      **
!     **************************************************************************
      USE RADIALFOCK_MODULE
USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID         ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR          ! #(GRID POINTS)
      CHARACTER(*),INTENT(IN)    :: KEY     
      REAL(8)    ,INTENT(IN)     :: RBOX        ! BOX RADIUS     
      REAL(8)    ,INTENT(IN)     :: AEZ         ! ATOMIC NUMBER
      INTEGER(4) ,INTENT(IN)     :: NX          ! X#(STATES)
      INTEGER(4) ,INTENT(OUT)    :: NB          ! #(STATES)
      INTEGER(4) ,INTENT(INOUT)  :: LOFI(NX)    ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(INOUT)  :: SOFI(NX)    ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(INOUT)  :: FOFI(NX)    ! OCCUPATION
      INTEGER(4) ,INTENT(INOUT)  :: NNOFI(NX)   ! #(NODES)
      REAL(8)    ,INTENT(OUT)    :: ETOT        ! TOTAL ENERGY
      REAL(8)    ,INTENT(INOUT)  :: POT(NR)     ! POTENTIAL
      TYPE(VFOCK_TYPE),INTENT(INOUT)  :: VFOCK  ! FOCK TERM
      REAL(8)   ,INTENT(OUT)     :: EOFI(NX)    ! ONE-PARTICLE ENERGIES
      REAL(8)   ,INTENT(OUT)     :: PHI(NR,NX)  ! ONE-PARTICLE WAVE FUNCTIONS
      REAL(8)   ,INTENT(OUT)     :: SPHI(NR,NX) ! SMALL COMPONENT
      REAL(8)   ,PARAMETER       :: TOL=1.D-3
      REAL(8)   ,PARAMETER       :: PRETOL=1.D-4
      REAL(8)   ,PARAMETER       :: XMAXTOL=1.D-8
      REAL(8)   ,PARAMETER       :: XAVTOL=1.D-8
      INTEGER(4),PARAMETER       :: NITER=1000
      LOGICAL(4),PARAMETER       :: TBROYDEN=.TRUE.
      LOGICAL(4),PARAMETER       :: TPR=.TRUE.
      REAL(8)   ,PARAMETER       :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER       :: Y0=1.D0/SQRT(4.D0*PI) !SPH. HARM. L=0
      REAL(8)   ,PARAMETER       :: C0LL=Y0               !GAUNT COEFF
      REAL(8)                    :: R(NR)
      REAL(8)                    :: RHO(NR)
      REAL(8)                    :: AUX(NR),AUX1(NR)   !AUXILIARY ARRAY
      REAL(8)                    :: MUX(NR)   !EXCHANGE ONLY POTENTIAL
      INTEGER(4)                 :: ITER
      REAL(8)                    :: XAV,XMAX
      REAL(8)                    :: XAVMIN,XMAXMIN,XDEMIN
      INTEGER(4)                 :: NCONV
      LOGICAL(4)                 :: CONVG
      REAL(8)                    :: EKIN,EH,EXC
      REAL(8)                    :: POTIN(NR)
      REAL(8)                    :: SVAR
      REAL(8)                    :: FMAX
      INTEGER(4)                 :: I,IB,JB,ISO,L,IR
      INTEGER(4)                 :: ISVAR,IARR(1)
      LOGICAL(4)                 :: TSTART  ! CALCULATE ANGULAR MOMENTA ETC
      LOGICAL(4)                 :: TREL    ! RELATIOVISTIC CALCULATION
      LOGICAL(4)                 :: TSO     ! CALCULATE WITH SPIN ORBIT COUPLING
      LOGICAL(4)                 :: TZORA   ! CALCULATE ZEROTH ORDER RELATIVISTIC CORRECTIONS
      LOGICAL(4)                 :: TFOCK   ! CALCULATE WITH FOCK EXCHANGE
!      INTEGER(4)          :: LMAP(19)=(/0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1/)
      INTEGER(4),PARAMETER       :: ZCORES(6)=(/2,10,18,36,54,86/)
      REAL(8)                    :: FTOT
      INTEGER(4)                 :: NBROYDENMEM
      REAL(8)                    :: BROYDENSTEP
      INTEGER(4)                 :: LRHOX=4
      INTEGER(4)                 :: IPERIOD,ISVAR1,IPOS
      CHARACTER(128)             :: STRING,STRING1
      REAL(8)                    :: SCALE
      REAL(8)                    :: EFOCK,EX
      REAL(8)                    :: EFS
      LOGICAL                    :: TSECOND
      REAL(8)                    :: RFOCK !EXTENT OF ORBITALS DEFINING FOCK TERM
      REAL(8)       ,ALLOCATABLE :: EOFI_FOCK(:)
!     **************************************************************************
!     CALL ATOMLIB$TEST_ATOMLIB$BOUNDSTATE()
!
!     ==========================================================================
!     == LMAP
!     == H  0          Z=1-2 
!     == LI 0,1,       Z=3-10
!     == NA 0,1,       Z=11-18
!     == K  0,2,1,     Z=19-36   3-D SERIES
!     == RB 0,2,1,     Z=37-54   4-D SERIES
!     == CS 0,3,2,1,   Z=55-86   4-F SERIES LANTHANIDES
!     == FR 0,3,2,1    Z=87-     5-F SERIES ACTINIDES
!     ==========================================================================
!
!     ==========================================================================
!     == RESOLVE KEY                                                          ==
!     ==========================================================================
      TREL=.FALSE.
      TSO=.FALSE.
      TZORA=.FALSE.
      TFOCK=.FALSE.
      TSTART=.FALSE.
      SCALE=0.D0
!
      STRING=KEY
      IPOS=1
      DO I=1,7
        IPOS=INDEX(STRING,',')
        IF(IPOS.NE.0) THEN
          STRING1=STRING(1:IPOS-1)
        ELSE
          STRING1=STRING
        END IF
        STRING=STRING(IPOS+1:)
        TREL=TREL.OR.TRIM(STRING1).EQ.'REL' 
        TSO=TSO.OR.STRING1.EQ.'SO' 
        TZORA=TZORA.OR.STRING1.EQ.'ZORA' 
        TSTART=TSTART.OR.STRING1.EQ.'START' 
        IF(STRING1(1:4).EQ.'FOCK') THEN
          TFOCK=.TRUE.
          IF(STRING(5:5).EQ.'=')READ(STRING1(6:),*)SCALE
        END IF
        IF(IPOS.EQ.0) EXIT
      ENDDO
      IF(TPR) THEN
        WRITE(*,FMT='("KEY: ",A)')TRIM(KEY)
        WRITE(*,FMT='("RELATIVISTIC EFFECTS SWITCHED ON? ",L5)')TREL
        WRITE(*,FMT='("SPIN-ORBIT COUPLING SWITCHED ON?  ",L5)')TSO
        WRITE(*,FMT='("ZORA SWITCHED ON?                 ",L5)')TZORA
        WRITE(*,FMT='("FOCK CONTRIBTION ON?              ",L5)')TFOCK
        IF(TFOCK) THEN
          WRITE(*,FMT='("PERCENT FOCK CONTRIBUTION: ",I5)')NINT(SCALE*100.D0)
        END IF
        WRITE(*,FMT='("START FROM SCRATCH?                ",L5)')TSTART
      END IF
      IF(.NOT.TREL.AND.INDEX(KEY,'NONREL').EQ.0) THEN
        CALL ERROR$MSG('TREL=F BUT "NONREL" NOT SPECIFIED IN KEY')
        CALL ERROR$CHVAL('KEY',KEY)
        CALL ERROR$STOP('ATOMLIB$AESCF')
      END IF
      IF(.NOT.TSO.AND.INDEX(KEY,'NONSO').EQ.0) THEN
        CALL ERROR$MSG('TSO=F BUT "NONSO" NOT SPECIFIED IN KEY')
        CALL ERROR$CHVAL('KEY',KEY)
        CALL ERROR$STOP('ATOMLIB$AESCF')
      END IF
      IF(.NOT.TZORA.AND.INDEX(KEY,'NONZORA').EQ.0) THEN
        CALL ERROR$MSG('TZORA=F BUT "NONZORA" NOT SPECIFIED IN KEY')
        CALL ERROR$CHVAL('KEY',KEY)
        CALL ERROR$STOP('ATOMLIB$AESCF')
      END IF
!
!     ==========================================================================
!     == INITIALIZATIONS                                                      ==
!     ==========================================================================
      NBROYDENMEM=2
      BROYDENSTEP=5.D-1
      CALL RADIAL$R(GID,NR,R)
!
      EOFI=0.D0
      IF(TSTART) THEN
        LOFI(:)=0
        FOFI(:)=0.D0
        SOFI(:)=0
        NB=0
        IPERIOD=0
        DO I=1,6
          IF(AEZ.LE.ZCORES(I)) EXIT
          IPERIOD=I
        ENDDO
        IB=0
        DO I=1,IPERIOD
          ISVAR=ZCORES(I)
          DO L=0,3
            IF(L.EQ.0)CALL PERIODICTABLE$GET(ISVAR,'OCC(S)',ISVAR1)
            IF(L.EQ.1)CALL PERIODICTABLE$GET(ISVAR,'OCC(P)',ISVAR1)
            IF(L.EQ.2)CALL PERIODICTABLE$GET(ISVAR,'OCC(D)',ISVAR1)
            IF(L.EQ.3)CALL PERIODICTABLE$GET(ISVAR,'OCC(F)',ISVAR1)
            IF(ISVAR1.EQ.0) CYCLE
            IB=IB+1
            LOFI(IB)=L
            FOFI(IB)=2.D0*REAL(2*L+1,KIND=8)
            SOFI(IB)=0
            IF(TSO) THEN
              IF(L.NE.0) THEN
                FOFI(IB)=REAL(2*L,KIND=8)
                SOFI(IB)=-1
                IB=IB+1
                LOFI(IB)=L
                FOFI(IB)=REAL(2*L+2,KIND=8)
                SOFI(IB)=1
              ELSE 
                SOFI(IB)=1
              END IF
            END IF
          ENDDO
        ENDDO
        DO L=0,3    
          IF(L.EQ.0)CALL PERIODICTABLE$GET(NINT(AEZ),'OCC(S)',ISVAR1)
          IF(L.EQ.1)CALL PERIODICTABLE$GET(NINT(AEZ),'OCC(P)',ISVAR1)
          IF(L.EQ.2)CALL PERIODICTABLE$GET(NINT(AEZ),'OCC(D)',ISVAR1)
          IF(L.EQ.3)CALL PERIODICTABLE$GET(NINT(AEZ),'OCC(F)',ISVAR1)
          IF(ISVAR1.EQ.0) CYCLE
          IB=IB+1
          LOFI(IB)=L
          FOFI(IB)=REAL(MIN(2*(2*L+1),ISVAR1),KIND=8)
          SOFI(IB)=0
          IF(TSO) THEN
            IF(L.NE.0) THEN
              SVAR=FOFI(IB)/REAL(2*(2*L+1),KIND=8)
              FOFI(IB)=SVAR*REAL(2*L,KIND=8)
              FOFI(IB+1)=SVAR*REAL(2*L+2,KIND=8)
              SOFI(IB+1)=1
              SOFI(IB)=-1
              LOFI(IB+1)=L
              IB=IB+1
            ELSE
              SOFI(IB)=1
            END IF
          END IF
        ENDDO
        NB=IB
!
!       ========================================================================
!       == SPECIAL TREATMENT OR EMPTY ATOMS                                   ==
!       ========================================================================
        IF(NB.EQ.0) THEN
          NB=1
          FOFI(NB)=0.D0
          LOFI(NB)=0
          SOFI(NB)=0
          IF(TSO)SOFI(NB)=1
        END IF
!
!       == CORRECT FOR NON-INTEGER ATOMIC NUMBERS ==============================
!       == OCCUPATIONS ARE FIRST FILLED ACCORDING TO NINT(AEZ). ================
!       == THIS POSSIBILY IS USED FOR DUMMY HYDROGEN ATOMS, THAT CARRY ONLY ====
!       == A FRACTIONAL NUCLEAR AND ELECTRONIC CHARGE ==========================
        FTOT=SUM(FOFI(:NB))   ! =NINT(AEZ)
        SVAR=AEZ-FTOT         ! =AEZ-NINT(AEZ)
!
!       == REMOVE ELECTRONS ====================================================
        IF(SVAR.LT.0.D0) THEN
          DO IB=NB,1,-1
            SVAR=SVAR+FOFI(IB)
            FOFI(IB)=MAX(0.D0,SVAR)
            SVAR=SVAR-FOFI(IB)
            IF(SVAR.GE.0.D0) EXIT 
          ENDDO

!       == ADD ELECTRONS =======================================================
        ELSE IF(SVAR.GT.0.D0) THEN
          DO IB=1,NB
            L=LOFI(IB)
            IF(TSO.AND.L.NE.0) THEN
              IF(SOFI(IB).EQ.-1) THEN
                FMAX=REAL(2*L,KIND=8)
              ELSE IF(SOFI(IB).EQ.1) THEN
                FMAX=REAL(2*L+2,KIND=8)
              END IF
            ELSE
              FMAX=REAL(2*(2*L+1),KIND=8)
            END IF
            SVAR=SVAR+FOFI(IB)
            FOFI(IB)=MIN(FMAX,SVAR)
            SVAR=SVAR-FOFI(IB)
            IF(SVAR.LE.0.D0) EXIT            
          ENDDO
        END IF
!
!       == CONSISTENCY CHECK  ==================================================
        FTOT=SUM(FOFI(:NB))
        IF(ABS(FTOT-AEZ).GT.1.D-8) THEN
          DO IB=1,NB
            WRITE(* &
        &       ,'("IB=",I2," L=",I1," SOFI=",I2," F=",F8.2," SUM(F)=",F8.2)') &
        &           IB,LOFI(IB),SOFI(IB),FOFI(IB),SUM(FOFI(:IB))
          ENDDO
          CALL ERROR$MSG('INCONSISTENT NUMBER OF ELECTRONS')
          CALL ERROR$R8VAL('AEZ ',AEZ)
          CALL ERROR$R8VAL('#(ELECTRONS) ',FTOT)
          CALL ERROR$R8VAL('#(ELECTRONS)-AEZ ',FTOT-AEZ)
          CALL ERROR$STOP('ATOMLIB$AESCF')
        END IF
!
        CALL RADIAL$NUCPOT(GID,NR,AEZ,POT)
!       == USE "HARD SPHERE BOUNDARY CONDITION" FOR THE POISSON EQUATION =======
!       == THAT IS CHOOSE THE POT(R)=0 FOR R>RBOX ==============================
        CALL RADIAL$VALUE(GID,NR,POT,RBOX,SVAR)
        POT(:)=POT(:)-SVAR
        DO IR=1,NR
          IF(R(IR).LT.RBOX) CYCLE
          POT(IR:)=0.D0
          EXIT
        ENDDO
!
!       ========================================================================
!       == DETERMINE NUMBER OF NODES FOR EACH SHELL                           ==
!       ========================================================================
        DO L=0,MAXVAL(LOFI)
          DO ISO=-1,1
            IF(TSO.AND.ISO.EQ.0) CYCLE
            IF(.NOT.TSO.AND.ISO.NE.0) CYCLE
            ISVAR=0
            DO IB=1,NB
              IF(LOFI(IB).NE.L.OR.SOFI(IB).NE.ISO) CYCLE
              NNOFI(IB)=ISVAR
              ISVAR=ISVAR+1
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     == SELF-CONSISTENCY LOOP                                                ==
!     ==========================================================================
      TSECOND=.FALSE.
1000  CONTINUE
      IF(.NOT.TSECOND) THEN
        WRITE(*,FMT='("DOING SCF ITERATIONS FOR ATOM WITH Z=",F10.5)')AEZ
      ELSE
        WRITE(*,FMT='("CONTINUING WITH FOCK TERM....")')
      END IF
      XAVMIN=1.D+12
      XMAXMIN=1.D+12
      XDEMIN=1.D+12
      NCONV=0
      CONVG=.FALSE.
      POTIN=POT
      CALL BROYDEN$NEW(NR,NBROYDENMEM,BROYDENSTEP)
      DO ITER=1,NITER
!
!       ========================================================================
!       == DETERMINE BOUND STATES FOR A GIVEN POTENTIAL AND ADD TO DENSITY    ==
!       ========================================================================
        IF(TFOCK.AND.TSECOND) THEN
          CALL ATOMLIB$BOUNDSTATESWITHHF(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI,POT &
     &                              ,VFOCK,TREL,TZORA,RBOX,PHI,SPHI)
        ELSE
          CALL ATOMLIB$BOUNDSTATES(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI,POT &
     &                            ,TREL,TZORA,RBOX,PHI,SPHI)
        END IF
!
!       ========================================================================
!       == ADD UP CHARGE DENSITY                                              ==
!       ========================================================================
        RHO(:)=0.D0
        DO IB=1,NB
          RHO(:)=RHO(:)+FOFI(IB)*C0LL*(PHI(:,IB)**2+SPHI(:,IB)**2)
        ENDDO
!
!       ========================================================================
!       ==  CALCULATE ENERGY                                                  ==
!       ========================================================================
        IF(TBROYDEN) POTIN=POT
        IF(CONVG) THEN
          POTIN=POT  ! SAVE INPUT POTENTIAL TO AVOID OVERWRITING
!         == DETERMINE KINETIC ENERGY ==========================================
          AUX(:)=0.D0
          DO IB=1,NB
            AUX(:)=AUX(:) &
       &          +(PHI(:,IB)**2+SPHI(:,IB)**2)*(EOFI(IB)-POT(:)*Y0)*FOFI(IB)
            IF(TFOCK.AND.TSECOND) THEN
               CALL RADIALFOCK$VPSI(GID,NR,VFOCK,LOFI(IB),PHI(:,IB),AUX1)
               AUX=AUX-PHI(:,IB)*AUX1(:)*FOFI(IB)
            END IF
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,EKIN)
          CALL ATOMLIB$BOXVOFRHO(GID,NR,RBOX,AEZ,RHO,POT,EH,EXC)
          ETOT=EKIN+EH+EXC
!
!         ======================================================================
!         == ESTIMATE THE ENERGY DUE TO THE FINITE NUCLEAR SIZE               ==
!         ======================================================================
          CALL ATOMLIB_EFINITENUCSIZE(GID,NR,RBOX,AEZ,RHO,EFS)
!
!         ======================================================================
!         == WORK OUT FOCK EXCHANGE ENERGY =====================================
!         ======================================================================
          IF(TFOCK.AND.TSECOND) THEN
            CALL DFT$SETL4('XCONLY',.TRUE.)
            CALL ATOMLIB$BOXVOFRHO(GID,NR,RBOX,AEZ,RHO,POT,EH,EX)
            CALL DFT$SETL4('XCONLY',.FALSE.)
!           == SAVE SCALE AND MUX TO RESTORE VFOCK FOR LATER USE ===============
            SCALE=VFOCK%SCALE
            MUX(:)=VFOCK%MUX(:)
!           ==  CHANGE VFOCK TO PURE FOCK TERM WITHOUT DOUBLE COUNTING =========
!           ==  AND CALCULATE FOCK TERM TO THE ENERGY ==========================
            VFOCK%SCALE=1.D0
            VFOCK%MUX=0.D0
            AUX(:)=0.D0
            DO IB=1,NB
              CALL RADIALFOCK$VPSI(GID,NR,VFOCK,LOFI(IB),PHI(:,IB),AUX1)
              AUX=AUX+0.5D0*PHI(:,IB)*AUX1(:)*FOFI(IB)
            ENDDO
            AUX(:)=AUX(:)*R(:)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,EFOCK)
!           == RESTORE VFOCK ==================================================
            VFOCK%SCALE=SCALE
            VFOCK%MUX(:)=MUX(:)
          ELSE
            EFOCK=0.D0
            EX=0.D0
          END IF
!
!         == PRINT ENERGIES ====================================================
          ETOT=EKIN+EH+EXC+SCALE*(EFOCK-EX)
!
          IF(TPR) THEN
!           __I PROVIDE RESULTS HERE WITH 10^-8 H PRECISION, BECAUSE____________
!           __ATOMIC CODES ARE TESTED IN THE MICRO-HARTREE ACCURACY_____________
            WRITE(*,FMT='(80("="),T20,"ENERGY REPORT OF ATOMLIB$AESCF")')
            WRITE(*,FMT='(30("."),T1,"TOTAL ENERGY:",T30,F15.8)')ETOT
            WRITE(*,FMT='(30("."),T1,"KINETIC ENERGY:",T30,F15.8)')EKIN
            WRITE(*,FMT='(30("."),T1,"HARTREE ENERGY:",T30,F15.8)')EH
            IF(TFOCK.AND.TSECOND) THEN
              WRITE(*,FMT='(30("."),T1,"EXACT XC MIXING FACTOR:",T30,F15.6)') &
     &                     SCALE
              WRITE(*,FMT='(30("."),T1,"MIXED XC ENERGY:",T30,F15.6)') &
     &                                                      EXC+SCALE*(EFOCK-EX)
              WRITE(*,FMT='(30("."),T1,"100% DFT XC ENERGY:",T30,F15.6)')EXC
              WRITE(*,FMT='(30("."),T1,"100% DFT EXCHANGE ENERGY:",T30' &
     &                                                      //',F15.6)') EX
            WRITE(*,FMT='(30("."),T1,"100% FOCK EXCHANGE ENERGY:",T30' &
     &                                                       //',F15.6)') EFOCK
            ELSE
              WRITE(*,FMT='(30("."),T1,"DFT XC ENERGY:",T30,F15.6)')EXC
            END IF
            WRITE(*,FMT='(30("."),T1,"FINITE NUCLEUS:",T30,F15.8)')EFS
            WRITE(*,FMT='(".... CALCULATION USES FINITE NUCLEUS")') 
            WRITE(*,FMT='(".... DO NOT ADD FINITE NUCLEUS CORRECTION")') 
          END IF

          POT=POTIN  ! RECOVER POT AS INPUT POTENTIAL
        END IF
!
!       ========================================================================
!       ==  EXIT IF CONVERGED                                                 ==
!       ========================================================================
        IF(CONVG) THEN
          IF(TFOCK.AND.(.NOT.TSECOND)) THEN
!           == LOOP MUST NOT BE EXITED BEFORE THE FOCK POTENTIAL IS CALCULATED =
!           == EXIT WILL BE DONE BELOW
          ELSE
!           == FINAL EXIT IS DONE HERE SO THAT THE POTENTIAL AND FOCK POTENTIAL
!           == CORRESPONDS TO INPUT POTENTIAL
            WRITE(*,FMT='("... SELFCONSISTENCY OBTAINED")')
            EXIT
          END IF
        END IF
!
!       ========================================================================
!       == CALCULATE OUTPUT POTENTIAL                                         ==
!       ========================================================================
        CALL ATOMLIB$BOXVOFRHO(GID,NR,RBOX,AEZ,RHO,POT,EH,EXC)
!
!       == DETERMINE WAVE FUNCTIONS FOR FOCK POTENTIAL =========================
        IF(TFOCK.AND.(TSECOND.OR.CONVG)) THEN
!!$! A TEST FOR HELIUM SHOWS THAT RESTRICTING RFOCK TO ONLY 1.1 OF THE COVALENT
!!$! RADIUS AFFECTS THE RESULTS. THE RESTRICTED CUTOFF HAS BEEN INTRODUCED 
!!$! BECAUSE OTHERWISE THE SELFCONSISTENCY OF THE PARTIAL WAVES IN 
!!$! MAKEPARTIALWAVES DOES NOT CONVERGE
!!$! A VALUE OF 3*RCOV GIVES MH ACCURACY FOR HE BUT NOT FOR THE HEAVIER ELEMENTS
          RFOCK=RBOX
          CALL PERIODICTABLE$GET(NINT(AEZ),'R(COV)',RFOCK); RFOCK=1.1D0*RFOCK   
!         CALL PERIODICTABLE$GET(NINT(AEZ),'R(COV)',RFOCK); RFOCK=3.0D0*RFOCK   
          RFOCK=MIN(RBOX,RFOCK)
          IF(.NOT.TSECOND) THEN
            ALLOCATE(EOFI_FOCK(NB))
            EOFI_FOCK(:)=EOFI(:NB)
            CALL ATOMLIB$BOUNDSTATES(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI_FOCK,POT &
     &                              ,TREL,TZORA,RFOCK,PHI,SPHI)
          ELSE
            CALL ATOMLIB$BOUNDSTATESWITHHF(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI,POT &
     &                                    ,VFOCK,TREL,TZORA,RFOCK,PHI,SPHI)
          END IF
          CALL RADIALFOCK$MAKEVFOCK(GID,NR,NB,LOFI,SOFI,FOFI,PHI,RFOCK &
     &                            ,LRHOX,VFOCK)
          RHO(:)=0.D0
          DO IB=1,NB
            RHO(:)=RHO(:)+FOFI(IB)*C0LL*(PHI(:,IB)**2+SPHI(:,IB)**2)
          ENDDO
          CALL ATOMLIB$BOXMUX(GID,NR,RFOCK,RHO,MUX)
          VFOCK%MUX(:)=MUX(:)          
          VFOCK%SCALE=SCALE
        END IF
!
!       ========================================================================
!       ==  EXIT IF CONVERGED (ONLY FOR TFOCK IN THE FIRST SEQUENCE)          ==
!       ========================================================================
        IF(CONVG) EXIT
!
!       ========================================================================
!       ==  GENERATE NEXT ITERATION USING D. G. ANDERSONS METHOD              ==
!       ========================================================================
        XAV=SQRT(SUM(R**3*(POT-POTIN)**2)/SUM(R**3))
        XMAX=MAXVAL(ABS(POT-POTIN))
        IF(TPR)PRINT*,ITER,' AV(POT-POTIN)=',XAV,' MAX:(POT-POTIN)=',XMAX,NCONV,TFOCK.AND.TSECOND
        NCONV=NCONV+1
        IF(XAV.LT.XAVMIN) THEN
          XAVMIN=XAV
          NCONV=0
        END IF
        IF(XMAX.LT.XMAXMIN) THEN
          XMAXMIN=XMAX
          NCONV=0
        END IF
!
!       == QUIT LOOP IF BOTH TOLERANCES ARE FULFILLED ==========================
        CONVG=(XMAX.LT.XMAXTOL).AND.(XAV.LT.XAVTOL)
        IF(TFOCK.AND.(.NOT.TSECOND)) THEN
          CONVG=(XMAX.LT.PRETOL).AND.(XAV.LT.PRETOL)
        END IF
!       == IF PREVIOUS CONDITIONS CANNOT BE MET DO THE BEST YOU CAN AND
!       == CECK IF MINIMUM REQUIREMENT IS FULFILLED
        CONVG=CONVG.OR.(XMAX.LT.TOL).AND.NCONV.GT.5
!
!       ========================================================================
!       ==  GENERATE NEXT ITERATION USING D. G. ANDERSONS METHOD              ==
!       ========================================================================
        CALL BROYDEN$STEP(NR,POTIN,POT-POTIN)
        POT=POTIN
      ENDDO

      CALL BROYDEN$CLEAR
      IF(.NOT.CONVG) THEN
        CALL ERROR$MSG('SELFCONSISTENCY LOOP NOT CONVERGED')
        CALL ERROR$R8VAL('AEZ',AEZ)
        CALL ERROR$R8VAL('XMAX',XMAX)
        CALL ERROR$R8VAL('TOLX',TOL)
        CALL ERROR$I4VAL('NITER',NITER)
        CALL ERROR$STOP('ATOMLIB$AESCF')
      END IF
!     == DFT CALCULATION IS CONVERGED. NOW CONVERGE WITH FOCK TERM =============
      IF(TFOCK.AND.(.NOT.TSECOND)) THEN
        TSECOND=.TRUE.
        CONVG=.FALSE.
IF(TPR) THEN
  PRINT*,'FIRST CONVERGENCE BEFORE APPLYING FOCK TERM'
  WRITE(*,FMT='(3A4,A10,A5,A20)')'IB','L','SO','F','#NODE','E'
  DO I=1,NB
    WRITE(*,FMT='(3I4,F10.2,I5,F20.3)')I,LOFI(I),SOFI(I),FOFI(I),NNOFI(I) &
 &                                   ,EOFI(I)
  ENDDO
END IF
        GOTO 1000
      END IF
IF(TPR) THEN
  CALL RADIALFOCK$PRINTVFOCK('VFOCK.DAT',VFOCK)
  WRITE(*,FMT='("FINAL ONE-PARTICLE ENERGIES FROM AESCF  ")')
  WRITE(*,FMT='(3A4,A10,A5,A20)')'IB','L','SO','F','#NODE','E'
  DO I=1,NB
    WRITE(*,FMT='(3I4,F10.2,I5,F20.3)')I,LOFI(I),SOFI(I),FOFI(I) &
 &   ,NNOFI(I),EOFI(I)
  ENDDO
  PRINT*,'#ITERATIONS ',ITER
END IF
!
!     ==========================================================================
!     == REORDER STATES ACCORDING TO ENERGY                                   ==
!     ==========================================================================
      DO IB=1,NB
        IARR=MINLOC(EOFI(IB:NB))
        JB=IARR(1)+IB-1
        IF(JB.EQ.IB) CYCLE
        ISVAR=LOFI(JB)
        LOFI(IB+1:JB)=LOFI(IB:JB-1)
        LOFI(IB)=ISVAR
        ISVAR=SOFI(JB)
        SOFI(IB+1:JB)=SOFI(IB:JB-1)
        SOFI(IB)=ISVAR
        ISVAR=NNOFI(JB)
        NNOFI(IB+1:JB)=NNOFI(IB:JB-1)
        NNOFI(IB)=ISVAR
        SVAR=EOFI(JB)
        EOFI(IB+1:JB)=EOFI(IB:JB-1)
        EOFI(IB)=SVAR
        SVAR=FOFI(JB)
        FOFI(IB+1:JB)=FOFI(IB:JB-1)
        FOFI(IB)=SVAR
        AUX=PHI(:,JB)
        PHI(:,IB+1:JB)=PHI(:,IB:JB-1)
        PHI(:,IB)=AUX
        AUX=SPHI(:,JB)
        SPHI(:,IB+1:JB)=SPHI(:,IB:JB-1)
        SPHI(:,IB)=AUX
      ENDDO
!STOP 'FORCED STOP IN AESCF'

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$BOUNDSTATES(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI,POT &
     &                              ,TREL,TZORA,RBOX,PHI,SPHI)
!     **************************************************************************
!     **  FINDS A SET OF BOUNDSTATES FOR A GIVEN POTENTIAL                    **
!     **  THE WAVE FUNCTIONS ARE NORMALIZED                                   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: NB      ! #(WAVE FUNCTION SHELLS)
      INTEGER(4) ,INTENT(IN)     :: LOFI(NB)! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SOFI(NB)! SWITCH FOR SPIN-ORBIT COUP.
      INTEGER(4) ,INTENT(IN)     :: NNOFI(NB)  !#(NODES)
      REAL(8)    ,INTENT(INOUT)  :: EOFI(NB)   !ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR)    !POTENTIAL
      LOGICAL(4) ,INTENT(IN)     :: TREL       !SWITCH FOR RELATIVISTIC CORR/
      LOGICAL(4) ,INTENT(IN)     :: TZORA      !SWITCH FOR ZORA
      REAL(8)    ,INTENT(IN)     :: RBOX       !BOX RADIUS
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,NB) !WAVE-FUNCTION (LARGE COMP.)
      REAL(8)    ,INTENT(OUT)    :: SPHI(NR,NB)!WAVE-FUNCTION (SMALL COMP.)
      REAL(8)                    :: DREL(NR)
      REAL(8)                    :: G(NR)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: SVAR
      REAL(8)                    :: AUX(NR),AUX1(NR)
      INTEGER(4)                 :: IB
      LOGICAL(4)                 :: TVARDREL
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == DEFINE HOW RELATIVISTIC EFFECTS ARE HANDLED                          ==
!     ==========================================================================
      DREL(:)=0.D0 
      TVARDREL=.FALSE.
      IF(TREL) THEN
        IF(TZORA) THEN
          CALL SCHROEDINGER$DREL(GID,NR,POT,0.D0,DREL)
        ELSE
          TVARDREL=.TRUE.
        END IF
      END IF
!     ==========================================================================
!     == LOOP OVER ALL STATES                                                 ==
!     ==========================================================================
      DO IB=1,NB
!
!       ========================================================================
!       == DETERMINE ENERGY AND LARGE COMPONENT                               ==
!       ========================================================================
        G(:)=0.D0
        CALL ATOMLIB$BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),0.D0,RBOX &
     &                        ,TVARDREL,DREL,G,NNOFI(IB),POT,EOFI(IB),PHI(:,IB))
!
!       ========================================================================
!       == DETERMINE SMALL COMPONENTS                                         ==
!       ========================================================================
        IF(TREL.AND.(.NOT.TZORA)) THEN
          G(:)=0.D0
          CALL SCHROEDINGER$DREL(GID,NR,POT,EOFI(IB),DREL)
          CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,LOFI(IB),SOFI(IB) &
     &                                       ,DREL,G,PHI(:,IB),SPHI(:,IB))
        ELSE
          SPHI(:,IB)=0.D0
        END IF
!
!       ========================================================================
!       == NORMALIZE                                                          ==
!       ========================================================================
        AUX(:)=R(:)**2*(PHI(:,IB)**2+SPHI(:,IB)**2)
        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
        CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
        PHI(:,IB)=PHI(:,IB)/SQRT(SVAR)
        SPHI(:,IB)=SPHI(:,IB)/SQRT(SVAR)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$BOUNDSTATESWITHHF(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI,POT &
     &                              ,VFOCK,TREL,TZORA,RBOX,PHI,SPHI)
!     **************************************************************************
!     **  FINDS A SET OF BOUNDSTATES FOR A GIVEN POTENTIAL                    **
!     ** THE WAVE FUNCTIONS ARE NORMALIZED                                    **
!     **************************************************************************
      USE RADIALFOCK_MODULE
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: NB      ! #(WAVE FUNCTION SHELLS)
      INTEGER(4) ,INTENT(IN)     :: LOFI(NB)! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SOFI(NB)! SWITCH FOR SPIN-ORBIT COUP.
      INTEGER(4) ,INTENT(IN)     :: NNOFI(NB)  !#(NODES)
      REAL(8)    ,INTENT(INOUT)  :: EOFI(NB)   !ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR)    !POTENTIAL
      LOGICAL(4) ,INTENT(IN)     :: TREL       !SWITCH FOR RELATIVISTIC CORR/
      LOGICAL(4) ,INTENT(IN)     :: TZORA       !SWITCH FOR ZORA
      REAL(8)    ,INTENT(IN)     :: RBOX       !BOX RADIUS
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,NB) !WAVE-FUNCTION (LARGE COMP.)
      REAL(8)    ,INTENT(OUT)    :: SPHI(NR,NB)!WAVE-FUNCTION (SMALL COMP.)
      TYPE(VFOCK_TYPE),INTENT(IN):: VFOCK  ! FOCK TERM
      REAL(8)                    :: DREL(NR)
      REAL(8)                    :: G(NR)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: SVAR
      REAL(8)                    :: AUX(NR),AUX1(NR)
      INTEGER(4)                 :: IB
      LOGICAL(4)                 :: TVARDREL
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      DREL(:)=0.D0 
      TVARDREL=TREL.AND.(.NOT.TZORA)
      IF(TREL.AND.TZORA) CALL SCHROEDINGER$DREL(GID,NR,POT,0.D0,DREL)
!
      DO IB=1,NB
!
!       ========================================================================
!       == DETERMINE ENERGY AND LARGE COMPONENT                               ==
!       ========================================================================
        G(:)=0.D0
        CALL ATOMLIB$BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),0.D0,RBOX,TVARDREL &
     &                         ,DREL,G,NNOFI(IB),POT,EOFI(IB),PHI(:,IB))
        CALL ATOMLIB$UPDATESTATEWITHHF(GID,NR,LOFI(IB),SOFI(IB),DREL,G,POT &
     &                                   ,VFOCK,RBOX,EOFI(IB),PHI(:,IB))
!
!       ========================================================================
!       == DETERMINE SMALL COMPONENTS                                         ==
!       ========================================================================
        IF(TREL.AND.(.NOT.TZORA)) THEN
          G(:)=0.D0
          CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,LOFI(IB),SOFI(IB) &
     &                                       ,DREL,G,PHI(:,IB),SPHI(:,IB))
        ELSE
          SPHI(:,IB)=0.D0
        END IF
!
!       ========================================================================
!       == NORMALIZE                                                          ==
!       ========================================================================
        AUX(:)=R(:)**2*(PHI(:,IB)**2+SPHI(:,IB)**2)
        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
        CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
        PHI(:,IB)=PHI(:,IB)/SQRT(SVAR)
        SPHI(:,IB)=SPHI(:,IB)/SQRT(SVAR)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RNS,RBOX &
     &                             ,TVARDREL,DREL,G,NN,POT,E,PHI)
!     **************************************************************************
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND         **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN.                     **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G                **
!     **                                                                      **
!     **  THE BOUNDARY CONDITION IS PHI(RBOX)=0                               **
!     **                                                                      **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE                 **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.         **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.         **
!     **                                                                      **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND           **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT           **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                         **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: RNS     ! NODES WITHIN RNS ARE IGNORED
      REAL(8)    ,INTENT(IN)     :: RBOX    ! BOX RADIUS
      LOGICAL(4) ,INTENT(IN)     :: TVARDREL! UPDATE RELATIVISTIC PARAMETER
      REAL(8)    ,INTENT(INOUT)  :: DREL(NR)! RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   ! INHOMOGENITY
      INTEGER(4) ,INTENT(IN)     :: NN      ! #(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) ! POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       ! ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) ! WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)    ,PARAMETER      :: TOL=1.D-8
      REAL(8)                    :: DER,DERO
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHIHOM(NR),PHIINHOM(NR),GHOM(NR)
      REAL(8)                    :: PHI1(NR)
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRMATCH 
      REAL(8)                    :: SVAR
      REAL(8)                    :: POT1(NR)
      REAL(8)                    :: DREL1(NR)
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20 !MAX. FACTOR IN THE WAVEFUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
      LOGICAL(4)                 :: THOM
      REAL(8)                    :: ROUT
      REAL(8)                    :: VAL1,VAL2,R1,R2
      INTEGER(4)                 :: IROUT,IRCL,IRBOX,IREND
!     **************************************************************************
!                                 CALL TRACE$PUSH('ATOMLIB$BOUNDSTATE')
      CALL RADIAL$R(GID,NR,R)
!     ==  R(IRBOX) IS THE FIRST GRIDPOINT JUST IOUTSIDE THE BOX
      IRBOX=1
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GE.RBOX) EXIT
      ENDDO
      DREL1=DREL
!          
!     ==========================================================================
!     ==========================================================================
!     ==  ITERATE TO BISECT ENERGY WITH A NODE AT RBOX                        ==
!     ==========================================================================
!     ==========================================================================
      ISTART=1
      X0=E
      DX=1.D-2
      CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      DO I=1,NITER
        E=X0
!        E=MIN(X0,EMAX)
!
!       ========================================================================
!       ==  CUT OFF THE POTENTIAL                                             ==
!       ========================================================================
!       == SPECIALRADS DETERMINES IRCL AND IROUT. 
!       == R(IRCL)=CLASS. TURNING POINT
!       == R(IROUT)= MAX RADIUS WITHOUT OVERFLOW
        CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT,E,IRCL,IROUT)
!       == BOUNDARY CONDITION PHI(ROUT)=0 ======================================
        IF(R(IROUT).LT.RBOX) THEN
          ROUT=R(IROUT)-1.D-5    !ENSURE THAT ROUT<R(IROUT)
        ELSE
          ROUT=RBOX
          IROUT=IRBOX
        END IF
!       ==  SET KINETIC ENERGY TO ZERO BEYOND ROUT TO AVOID AN OVERFLOW ========
!       == ATTENTION: THERE IS A STEP IN THE POTENTIAL =========================
        POT1(:)=POT(:)
        POT1(IROUT+3:)=E
!
!       ========================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       ========================================================================
        IDIR=1
        IF(TVARDREL)CALL SCHROEDINGER$DREL(GID,NR,POT1,E,DREL1)
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL1,SO,G,L,E,IDIR,PHI)
!       == CHECK FOR OVERFLOW
        IF(.NOT.(PHI(IROUT).GT.0.OR.PHI(IROUT).LE.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('NN',NN)
          CALL ERROR$STOP('ATOMLIB$BOUNDSTATE')
        END IF
!
!       ========================================================================
!       == ESTIMATE PHASE SHIFT                                               ==
!       ========================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RNS,ROUT,Z0)
        Z0=Z0-REAL(NN+1)
        IF(ABS(2.D0*DX).LE.TOL) EXIT
!       ========================================================================
!       ==  BISECTION                                                         ==
!       ========================================================================
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('ATOMLIB$BOUNDSTATE')
      END IF
!
!     ==========================================================================
!     ==========================================================================
!     == INTEGRATE OUTWARD AND INVARD                                         ==
!     ==========================================================================
!     ==========================================================================
!
!     ==========================================================================
!     ==  DETERMINE MATCHING POINT                                            ==
!     ==========================================================================
      IRMATCH=IRCL
      IF(R(IRCL).GT.5.D0) THEN
        CALL RADIAL$XOFR(GID,5.D0,SVAR)
        IRMATCH=INT(SVAR)
      END IF
!
!     ==========================================================================
!     ==  INTEGRATE INWARD                                                    ==
!     ==========================================================================
      IF(IRMATCH.LT.IROUT) THEN
        THOM=MAXVAL(ABS(G(:))).EQ.0.D0
        IDIR=-1
!       ==  HOMOGENEOUS SOLUTION THAT FULFILLS THE OUTER BOUNDARY CONDITION   ==
!       ==  INTEGRATE INWARD AT THE GIVEN ENERGY WITH SLIGHTLY DIFFERENT      ==
!       ==  BOUNDARY CONDITIONS AND SUPERIMPOSE THEM SO THAT OUTER BOUNDARY   ==
!       ==  CONDITION IS EXACTLY FULFILLED                                    ==
        IREND=MIN(IROUT,NR-3)
        GHOM(:)=0.D0
        GHOM(IREND+1)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL1,SO,GHOM,L,E,IDIR,PHIHOM)
        PHIHOM(:)=PHIHOM(:)/PHIHOM(IRMATCH)
        GHOM(:)=0.D0
        GHOM(IREND+2)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL1,SO,GHOM,L,E,IDIR,PHI1)
        PHI1(:)=PHI1(:)/PHI1(IRMATCH)
!       == EXTRAPOLATE =======================
        IF(IREND.LT.IROUT) THEN
          R1=R(IREND)
          R2=R(IREND+1)
          VAL1=PHIHOM(IREND)
          VAL2=PHIHOM(IREND+1)
          PHIHOM(IREND:)=VAL1+(VAL2-VAL1)/(R2-R1)*(R(IREND:)-R1)
          VAL1=PHI1(IREND)
          VAL2=PHI1(IREND+1)
          PHI1(IREND:)=VAL1+(VAL2-VAL1)/(R2-R1)*(R(IREND:)-R1)
        END IF
!       == FULFILL OUTER BOUNDARY CONDITION ====================================
        CALL RADIAL$VALUE(GID,NR,PHIHOM,ROUT,VAL1)
        CALL RADIAL$VALUE(GID,NR,PHI1,ROUT,VAL2)
        SVAR=VAL1+VAL2
        VAL1=VAL1/SVAR
        VAL2=VAL2/SVAR
        PHIHOM(:)=VAL2*PHIHOM(:)-VAL1*PHI1(:)
        PHIHOM(:)=PHIHOM(:)/PHIHOM(IRMATCH)
!       
!       == INHOMOGENEOUS SOLUTION WITH CORRECT BOUNDARY CONDITIONS =============
        IF(.NOT.THOM) THEN     
          GHOM=G   !AVOID OVERFLOW BY REMOVIN INHOMOGENEITY OUTSIDE OF THE BOX
          GHOM(IREND+2:)=0.D0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL1,SO,GHOM,L,E,IDIR &
      &                              ,PHIINHOM)
          CALL RADIAL$VALUE(GID,NR,PHIINHOM,ROUT,VAL1)
          CALL RADIAL$VALUE(GID,NR,PHI1,ROUT,VAL2)
          PHIINHOM(:)=PHIINHOM(:)-VAL1/VAL2*PHI1(:)
         ELSE
          PHIINHOM(:)=0.D0
        END IF
!
!       =======================================================================
!       ==  MATCH SOLUTION INSIDE AND OUTSIDE WITH VALUE                     ==
!       =======================================================================
        SVAR=(PHI(IRMATCH)-PHIINHOM(IRMATCH))/PHIHOM(IRMATCH)
        PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM(:)
        CALL RADIAL$DERIVATIVE(GID,NR,PHI,R(IRMATCH),DER)
        CALL RADIAL$DERIVATIVE(GID,NR,PHIINHOM,R(IRMATCH),DERO)
        SVAR=(DERO-DER)/PHI(IRMATCH)
        PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
!
!       ========================================================================
!       ==  CLEAN UP GLITCH FROM CREATING THE INWARD SOLUTION                 ==
!       ==  (IREND<IRBOX), THERE ARE STILL DATA BEYOND IREND, WHICH MIMICK    ==
!       ==  A NODE AT IREND. THIS CONFUSES NODE COUNTING LATERON.             ==
!       ========================================================================
        IF(IREND+1.LE.IRBOX) THEN
          PHI(IREND:)=0.D0
        END IF
      END IF
!!$!
! DO NOT NORMALIZE!!!! IT DOES NOT MAKE SENSE FOR THE INHOMOGENEOUS SOLUTION
!                                 CALL TRACE$POP()
      RETURN
      END SUBROUTINE ATOMLIB$BOUNDSTATE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$PAWBOUNDSTATE(GID,NR,L,NN,RBOX,PSPOT,NPRO,PRO,DH,DO,G &
     &                                ,E,PHI)
!     **************************************************************************
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND         **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN WITH R<RMINNODE.     **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G                **
!     **                                                                      **
!     **  THE BOUNDARY CONDITION IS PHI(RBOX)=0                               **
!     **                                                                      **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE                 **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.         **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.         **
!     **                                                                      **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND           **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT           **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                         **
!     **                                                                      **
!     ** WARNING! ONLY NODES AT R>RMINNODE ARE COUNTED                        **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: NN      ! #(NODES)
      REAL(8)    ,INTENT(IN)     :: RBOX    ! BOX RADIUS
      INTEGER(4) ,INTENT(IN)     :: NPRO      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: PSPOT(NR) !POTENTIAL
      REAL(8)    ,INTENT(IN)     :: PRO(NR,NPRO)
      REAL(8)    ,INTENT(IN)     :: DH(NPRO,NPRO)
      REAL(8)    ,INTENT(IN)     :: DO(NPRO,NPRO)
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      LOGICAL(4) ,PARAMETER      :: TWRITE=.FALSE.
      INTEGER(4) ,PARAMETER      :: NITER=1000
      REAL(8)    ,PARAMETER      :: TOL=1.D-12
      REAL(8)    ,PARAMETER      :: RMATCHN=4.D0 ! MIN MATCHING RADIUS
      REAL(8)    ,PARAMETER      :: RMINNODE=1.D-2 ! MIN. RADIUS FOR NODE COUNT
!     == RMINNODE MUST BE SUFFICIENTLY SMALL TO ALLOW FOR REAL NODES SUCH AS 
!     == FOR H AND LI. ON THE OTHER HAND IT SHOULD PROJECT OUT ANY SPURIOUS 
!     == NODES DUE TO THE VIOLATION OF THE NODAL THEOREM BY NON-LOCAL POTENTIALS
!     == VALUE HAS BEEN CHANGED FROM 1 TO 0.5 ON JULY 7, 2018. PEB.
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHI1(NR),PHI2(NR)
      REAL(8)                    :: DREL(NR),GHOM(NR),PHIHOM(NR)
      REAL(8)                    :: PHIINHOM(NR)
      REAL(8)                    :: VAL1,VAL2,SVAR
      INTEGER(4)                 :: IR,IRMATCH,SO,IDIR,I
      INTEGER(4)                 :: ITER
      LOGICAL                    :: THOM
      REAL(8)                    :: Y2,Y1,X2,X1,DER
      INTEGER(4)                 :: IRBOX
!     **************************************************************************
                                 CALL TRACE$PUSH('ATOMLIB$PAWBOUNDSTATE')
      IF(TWRITE) THEN
        WRITE(*,FMT='(80("+"),T20," ATOMLIB$PAWBOUNDSTATE ")')
        WRITE(*,FMT='("L=",I3," NN=",I3," RBOX=",F9.3," NPRO=",I3)') &
     &              L,NN,RBOX,NPRO
        WRITE(*,FMT='("E=",F10.5)')E
        DO I=1,NPRO
          WRITE(*,FMT='("DH=",10F10.5)')DH(I,:)
        ENDDO
        DO I=1,NPRO
          WRITE(*,FMT='("DO=",10F10.5)')DO(I,:)
        ENDDO
      END IF
!
      PHI=0.D0
      CALL RADIAL$R(GID,NR,R)
!     ==  R(IRBOX) IS THE FIRST GRIDPOINT JUST OUTSIDE THE BOX =================
      IF(RBOX.GT.R(NR-3)) THEN
        CALL ERROR$MSG('GRID TOO SMALL FOR CHOSEN BOX RADIUS')
        CALL ERROR$R8VAL('RBOX',RBOX)
        CALL ERROR$R8VAL('R(NR)',R(NR))
        CALL ERROR$R8VAL('R(NR-2)',R(NR-2))
        CALL ERROR$STOP('ATOMLIB$PAWBOUNDSTATE')
      END IF
      IRBOX=1
      DO IR=1,NR-3
        IRBOX=IR
        IF(R(IR).GE.RBOX) EXIT
      ENDDO
!          
!     ==========================================================================
!     ==========================================================================
!     ==  PRESELECT A WINDOW FOR BISECTION                                    ==
!     ==  DO SMALL STEPS TO AVOID COMING CLOSE TO A GHOST STATE               ==
!     ==========================================================================
!     ==========================================================================
      Z0=333.333D0 ! FIRST VALUE IS MEANINGLESS
      DX=1.D-2
      X0=E-DX
!!$IF(L.EQ.0.AND.NN.EQ.1) THEN
!!$X0=X0+1.D-2
!!$END IF
!!$IF(L.EQ.0.AND.NN.EQ.1) THEN
!!$OPEN(UNIT=1005,FILE='XOUT')
!!$DO ITER=1,20000
!!$  E=-2.D0+2.5D-4*REAL(ITER,KIND=8)
!!$  CALL ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO,DH,DO,G,PHI)
!!$  CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RMINNODE,RBOX,Z0)
!!$  WRITE(1005,*)E,Z0
!!$ENDDO
!!$CLOSE(1005)
!!$CALL ERROR$STOP('---')
!!$END IF
      ZM=0.D0 ! INITIALIZATION TO MAKE COMPILER HAPPY
      DO ITER=1,NITER
        E=X0
!       ========================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       ========================================================================
        CALL ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO,DH,DO,G,PHI)
!       == CHECK FOR OVERFLOW ==================================================
        IF(.NOT.(PHI(IRBOX+2).GT.0.OR.PHI(IRBOX+2).LE.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$STOP('ATOMLIB$PAWBOUNDSTATE')
        END IF
!
!       ========================================================================
!       == ESTIMATE PHASE SHIFT                                               ==
!       ========================================================================
!       == NODES WITH R<RMINNODE ARE NOT COUNTED, BECAUSE THERE IS NO  =========
!       == NODAL THEOREM FOR NON-LOCAL POTENTIALS. =============================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RMINNODE,RBOX,Z0)
        Z0=Z0-REAL(NN+1,KIND=8)
        IF(TWRITE)WRITE(*,FMT='("LOOP1 X0=",F10.5," Z0=",E12.3)')X0,Z0
        IF(Z0.GT.0.D0) THEN
          PHI1(:)=PHI(:)
        ELSE
          PHI2(:)=PHI(:)
        END IF
        IF(ITER.GT.1) THEN
          IF(Z0*ZM.LT.0.D0) EXIT
        END IF
!!$IF(L.EQ.0.AND.NN.EQ.1) THEN
!!$ CALL ATOMLIB_WRITEPHI('ERROR_PAWPSI.DAT',GID,NR,1,PHI)
!!$ CALL ERROR$I4VAL('L',L)
!!$ CALL ERROR$I4VAL('NN',NN)
!!$ CALL ERROR$R8VAL('RMINNODE',RMINNODE)
!!$ CALL ERROR$R8VAL('RBOX',RBOX)
!!$ CALL ERROR$R8VAL('X0',X0)
!!$ CALL ERROR$R8VAL('Z0',Z0)
!!$ CALL ERROR$STOP('ATOMLIB$PAWBOUNDSTATE')
!!$ENDIF
        IF(ITER.EQ.NITER) THEN
          CALL ERROR$MSG('SEARCH FOR BISECTION WINDOW FAILED')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('NN',NN)
          CALL ERROR$R8VAL('XM',XM)
          CALL ERROR$R8VAL('X0',X0)
          CALL ERROR$R8VAL('ZM',ZM)
          CALL ERROR$R8VAL('Z0',Z0)
          CALL ERROR$STOP('ATOMLIB$PAWBOUNDSTATE')
        END IF
        IF(ITER.EQ.1) DX=SIGN(DX,-Z0)
        XM=X0
        ZM=Z0
        X0=X0+DX
      ENDDO
!          
!     ==========================================================================
!     ==========================================================================
!     ==  ITERATE TO BISECT ENERGY WITH A NODE AT RBOX                        ==
!     ==========================================================================
!     ==========================================================================
      ISTART=1
      X0=E
!     == DX IS DEFINED BY THE PREVIOUS LOOP
      CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      DO ITER=1,NITER
        E=X0
!
!       ========================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       ========================================================================
        CALL ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO,DH,DO,G,PHI)
!       == CHECK FOR OVERFLOW ==================================================
        IF(.NOT.(PHI(IRBOX+2).GT.0.OR.PHI(IRBOX+2).LE.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$STOP('ATOMLIB$PAWBOUNDSTATE')
        END IF
!
!       ========================================================================
!       == ESTIMATE PHASE SHIFT                                               ==
!       ========================================================================
!       == NODES WITH R<1 A_BOHR ARE NOT COUNTED, BECAUSE THERE IS NO  =========
!       == NODAL THEOREM FOR NON-LOCAL POTENTIALS. =============================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RMINNODE,RBOX,Z0)
        Z0=Z0-REAL(NN+1,KIND=8)
        IF(TWRITE)WRITE(*,FMT='("LOOP2 X0=",F10.5," Z0=",E12.3)')X0,Z0

        IF(ABS(2.D0*DX).LE.TOL) EXIT
        IF(Z0.GT.0.D0) THEN
          PHI1(:)=PHI(:)
        ELSE
          PHI2(:)=PHI(:)
        END IF
!       ========================================================================
!       ==  BISECTION                                                         ==
!       ========================================================================
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('ATOMLIB$PAWBOUNDSTATE')
      END IF
!
!     ==========================================================================
!     ==  AVERAGE BOTH BOUNDS OF BISECTION                                    ==
!     ==========================================================================
      X1=R(IRBOX-1)
      X2=R(IRBOX)
      Y1=PHI1(IRBOX-1)
      Y2=PHI1(IRBOX)
      DER=(Y2-Y1)/(X2-X1)
      VAL1=Y1+DER*(RBOX-X1)
      Y1=PHI2(IRBOX-1)
      Y2=PHI2(IRBOX)
      DER=(Y2-Y1)/(X2-X1)
      VAL2=Y1+DER*(RBOX-X1)
      SVAR=VAL2-VAL1
      VAL1=VAL1/SVAR
      VAL2=VAL2/SVAR
      PHI=PHI1*VAL2-PHI2*VAL1
CALL TRACE$POP()
RETURN
!
!     ==========================================================================
!     ==  DETERMINE MATCHING POINT                                            ==
!     ==========================================================================
      DO IR=1,NR
        IRMATCH=IR
        IF(R(IR).LT.RMATCHN)CYCLE
        IF(SUM(ABS(PRO(IR,:))).LT.1.D-8) THEN
          EXIT
!          IF(ABS(PHI2(IR)-PHI1(IR)).GT.1.D-3*ABS(PHI1(IR))) EXIT
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  INTEGRATE INWARD                                                    ==
!     ==========================================================================
      DREL(:)=0.D0
      SO=0
      IF(IRMATCH.LT.IRBOX) THEN
        THOM=MAXVAL(ABS(G(:))).EQ.0.D0
        IDIR=-1
!       ==  HOMOGENEOUS SOLUTION THAT FULFILLS THE OUTER BOUNDARY CONDITION   ==
!       ==  INTEGRATE INWARD AT THE GIVEN ENERGY WITH SLIGHTLY DIFFERENT      ==
!       ==  BOUNDARY CONDITIONS AND SUPERIMPOSE THEM SO THAT OUTER BOUNDARY   ==
!       ==  CONDITION IS EXACTLY FULFILLED                                    ==
        GHOM(:)=0.D0
        GHOM(IRBOX+1)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
        PHIHOM(:)=PHIHOM(:)/PHIHOM(IRMATCH)
        GHOM(:)=0.D0
        GHOM(IRBOX+2)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,GHOM,L,E,IDIR,PHI1)
        PHI1(:)=PHI1(:)/PHI1(IRMATCH)
!       == FULFILL OUTER BOUNDARY CONDITION ====================================
        CALL RADIAL$VALUE(GID,NR,PHIHOM,RBOX,VAL1)
        CALL RADIAL$VALUE(GID,NR,PHI1,RBOX,VAL2)
        SVAR=VAL1+VAL2
        VAL1=VAL1/SVAR
        VAL2=VAL2/SVAR
        PHIHOM(:)=VAL2*PHIHOM(:)-VAL1*PHI1(:)
        PHIHOM(:)=PHIHOM(:)/PHIHOM(IRMATCH)
!       
!       == INHOMOGENEOUS SOLUTION WITH CORRECT BOUNDARY CONDITIONS =============
        IF(.NOT.THOM) THEN     
          CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,G,L,E,IDIR,PHIINHOM)
          CALL RADIAL$VALUE(GID,NR,PHIINHOM,RBOX,VAL1)
          CALL RADIAL$VALUE(GID,NR,PHI1,RBOX,VAL2)
          PHIINHOM(:)=PHIINHOM(:)-VAL1/VAL2*PHI1(:)
        ELSE
          PHIINHOM(:)=0.D0
        END IF
!
!       =======================================================================
!       ==  MATCH SOLUTION INSIDE AND OUTSIDE WITH VALUE                     ==
!       =======================================================================
        SVAR=(PHI(IRMATCH)-PHIINHOM(IRMATCH))/PHIHOM(IRMATCH)
        PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM(:)
        PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
!
!!$        CALL RADIAL$DERIVATIVE(GID,NR,PHI,R(IRMATCH),DER)
!!$        CALL RADIAL$DERIVATIVE(GID,NR,PHIINHOM,R(IRMATCH),DERO)
!!$        SVAR=(DERO-DER)/PHI(IRMATCH)
!!$        PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
      END IF
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$ORTHOBOUNDSTATE(GID,NR,L,NN,RBOX,PSPOT,NPRO,PRO,G &
     &                                ,E,PHI)
!     **************************************************************************
!     ** CAUTION! THIS ROUTINE SEEMS TO FAIL!!                                **
!     **                                                                      **
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND         **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN.                     **
!     **  AND THE CONDITION THAT THE SOLUTION MUST BE ORTHOGONAL TO THE       **
!     **  PROJECTOR FUNCTIONS PRO.
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G                **
!     **                                                                      **
!     **  THE BOUNDARY CONDITION IS PHI(RBOX)=0                               **
!     **                                                                      **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE                 **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.         **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.         **
!     **                                                                      **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND           **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT           **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                         **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: NN      ! #(NODES)
      REAL(8)    ,INTENT(IN)     :: RBOX    ! BOX RADIUS
      INTEGER(4) ,INTENT(IN)     :: NPRO      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: PSPOT(NR) !POTENTIAL
      REAL(8)    ,INTENT(IN)     :: PRO(NR,NPRO)
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      LOGICAL(4) ,PARAMETER      :: TWRITE=.FALSE.
      INTEGER(4) ,PARAMETER      :: NITER=100
      REAL(8)    ,PARAMETER      :: TOL=1.D-12
      REAL(8)    ,PARAMETER      :: RMATCHN=4.D0 ! MIN MATCHING RADIUS
      REAL(8)    ,PARAMETER      :: RMINNODE=1.D0 ! MIN RADIUS FOR NODE COUNT
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHI1(NR),PHI2(NR)
      REAL(8)                    :: DREL(NR),GHOM(NR),PHIHOM(NR)
      REAL(8)                    :: PHIINHOM(NR)
      REAL(8)                    :: VAL1,VAL2,SVAR
      INTEGER(4)                 :: IR,IRMATCH,SO,IDIR,I
      INTEGER(4)                 :: ITER
      LOGICAL                    :: THOM
      REAL(8)                    :: Y2,Y1,X2,X1,DER
      INTEGER(4)                 :: IRBOX
!     **************************************************************************
                                 CALL TRACE$PUSH('ATOMLIB$ORTHOBOUNDSTATE')
      IF(TWRITE) THEN
        WRITE(*,FMT='(80("+"),T20," ATOMLIB$ORTHOBOUNDSTATE ")')
        WRITE(*,FMT='("L=",I3," NN=",I3," RBOX=",F9.3," NPRO=",I3)') &
     &              L,NN,RBOX,NPRO
        WRITE(*,FMT='("E=",F10.5)')E
      END IF
!
      PHI=0.D0
      CALL RADIAL$R(GID,NR,R)
!     ==  R(IRBOX) IS THE FIRST GRIDPOINT JUST OUTSIDE THE BOX =================
      IF(RBOX.GT.R(NR-3)) THEN
        CALL ERROR$MSG('GRID TOO SMALL FOR CHOSEN BOX RADIUS')
        CALL ERROR$R8VAL('RBOX',RBOX)
        CALL ERROR$R8VAL('R(NR)',R(NR))
        CALL ERROR$R8VAL('R(NR-2)',R(NR-2))
        CALL ERROR$STOP('ATOMLIB$ORTHOBOUNDSTATE')
      END IF
      IRBOX=1
      DO IR=1,NR-3
        IRBOX=IR
        IF(R(IR).GE.RBOX) EXIT
      ENDDO
!          
!     ==========================================================================
!     ==========================================================================
!     ==  PRESELECT A WINDOW FOR BISECTION                                    ==
!     ==  DO SMALL STEPS TO AVOID COMING CLOSE TO A GHOST STATE               ==
!     ==========================================================================
!     ==========================================================================
      X0=E
      Z0=333.333D0 ! FIRST VALUE IS MEANINGLESS
      DX=1.D-2
      DO ITER=1,NITER
        E=X0
!       ========================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       ========================================================================
        CALL ATOMLIB_ORTHODER(GID,NR,L,E,PSPOT,NPRO,PRO,G,PHI)
!       == CHECK FOR OVERFLOW ==================================================
        IF(.NOT.(PHI(IRBOX+2).GT.0.OR.PHI(IRBOX+2).LE.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$STOP('ATOMLIB$ORTHOBOUNDSTATE')
        END IF
!
!       ========================================================================
!       == ESTIMATE PHASE SHIFT                                               ==
!       ========================================================================
!       == NODES WITH R<RMINNODE ARE NOT COUNTED, BECAUSE THERE IS NO  =========
!       == NODAL THEOREM FOR NON-LOCAL POTENTIALS. =============================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RMINNODE,RBOX,Z0)
!PRINT*,'E,LOGDER',E,Z0,NPRO
!!$CALL ATOMLIB_WRITEPHI('TEST_PHI',GID,NR,1,PHI)
!!$CALL ATOMLIB_ORTHODER(GID,NR,L,E,PSPOT,0,PRO,G,PHI)
!!$PRINT*,'E,LOGDER',E,Z0,NPRO
!!$CALL ATOMLIB_WRITEPHI('TEST_PHI0',GID,NR,1,PHI)
!!$CALL ATOMLIB_ORTHODER(GID,NR,L,E,PSPOT,1,PRO(:,1),G,PHI)
!!$PRINT*,'E,LOGDER',E,Z0,NPRO
!!$CALL ATOMLIB_WRITEPHI('TEST_PHI1',GID,NR,1,PHI)
!!$STOP 'FORCED'
        Z0=Z0-REAL(NN+1,KIND=8)
        IF(Z0.GT.0.D0) THEN
          PHI1(:)=PHI(:)
        ELSE
          PHI2(:)=PHI(:)
        END IF
        IF(ITER.GT.1.AND.Z0*ZM.LT.0.D0) EXIT
        IF(ITER.EQ.1) DX=SIGN(DX,-Z0)
        XM=X0
        ZM=Z0
        X0=X0+DX
      ENDDO
!          
!     ==========================================================================
!     ==========================================================================
!     ==  ITERATE TO BISECT ENERGY WITH A NODE AT RBOX                        ==
!     ==========================================================================
!     ==========================================================================
      ISTART=1
      X0=E
      Z0=333.333D0   ! FIRST VALUE IS MEANINGLESS
!     == DX IS DEFINED BY THE PREVIOUS LOOP
      CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      DO ITER=1,NITER
        E=X0
!
!       ========================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       ========================================================================
        CALL ATOMLIB_ORTHODER(GID,NR,L,E,PSPOT,NPRO,PRO,G,PHI)
!       == CHECK FOR OVERFLOW ==================================================
        IF(.NOT.(PHI(IRBOX+2).GT.0.OR.PHI(IRBOX+2).LE.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$STOP('ATOMLIB$ORTHOBOUNDSTATE')
        END IF
!
!       ========================================================================
!       == ESTIMATE PHASE SHIFT                                               ==
!       ========================================================================
!       == NODES WITH R<1 A_BOHR ARE NOT COUNTED, BECAUSE THERE IS NO  =========
!       == NODAL THEOREM FOR NON-LOCAL POTENTIALS. =============================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,1.D0,RBOX,Z0)
        Z0=Z0-REAL(NN+1,KIND=8)
        IF(TWRITE)WRITE(*,FMT='("X0=",F10.5," Z0=",E12.3)')X0,Z0

        IF(ABS(2.D0*DX).LE.TOL) EXIT
        IF(Z0.GT.0.D0) THEN
          PHI1(:)=PHI(:)
        ELSE
          PHI2(:)=PHI(:)
        END IF
!       ========================================================================
!       ==  BISECTION                                                         ==
!       ========================================================================
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('ATOMLIB$ORTHOBOUNDSTATE')
      END IF
!
!     ==========================================================================
!     ==  AVERAGE BOTH BOUNDS OF BISECTION                                    ==
!     ==========================================================================
      X1=R(IRBOX-1)
      X2=R(IRBOX)
      Y1=PHI1(IRBOX-1)
      Y2=PHI1(IRBOX)
      DER=(Y2-Y1)/(X2-X1)
      VAL1=Y1+DER*(RBOX-X1)
      Y1=PHI2(IRBOX-1)
      Y2=PHI2(IRBOX)
      DER=(Y2-Y1)/(X2-X1)
      VAL2=Y1+DER*(RBOX-X1)
      SVAR=VAL2-VAL1
      VAL1=VAL1/SVAR
      VAL2=VAL2/SVAR
      PHI=PHI1*VAL2-PHI2*VAL1
RETURN
!
!     ==========================================================================
!     ==  DETERMINE MATCHING POINT                                            ==
!     ==========================================================================
      DO IR=1,NR
        IRMATCH=IR
        IF(R(IR).LT.RMATCHN)CYCLE
        IF(SUM(ABS(PRO(IR,:))).LT.1.D-8) THEN
          EXIT
!          IF(ABS(PHI2(IR)-PHI1(IR)).GT.1.D-3*ABS(PHI1(IR))) EXIT
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  INTEGRATE INWARD                                                    ==
!     ==========================================================================
      DREL(:)=0.D0
      SO=0
      IF(IRMATCH.LT.IRBOX) THEN
        THOM=MAXVAL(ABS(G(:))).EQ.0.D0
        IDIR=-1
!       ==  HOMOGENEOUS SOLUTION THAT FULFILLS THE OUTER BOUNDARY CONDITION   ==
!       ==  INTEGRATE INWARD AT THE GIVEN ENERGY WITH SLIGHTLY DIFFERENT      ==
!       ==  BOUNDARY CONDITIONS AND SUPERIMPOSE THEM SO THAT OUTER BOUNDARY   ==
!       ==  CONDITION IS EXACTLY FULFILLED                                    ==
        GHOM(:)=0.D0
        GHOM(IRBOX+1)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
        PHIHOM(:)=PHIHOM(:)/PHIHOM(IRMATCH)
        GHOM(:)=0.D0
        GHOM(IRBOX+2)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,GHOM,L,E,IDIR,PHI1)
        PHI1(:)=PHI1(:)/PHI1(IRMATCH)
!       == FULFILL OUTER BOUNDARY CONDITION ====================================
        CALL RADIAL$VALUE(GID,NR,PHIHOM,RBOX,VAL1)
        CALL RADIAL$VALUE(GID,NR,PHI1,RBOX,VAL2)
        SVAR=VAL1+VAL2
        VAL1=VAL1/SVAR
        VAL2=VAL2/SVAR
        PHIHOM(:)=VAL2*PHIHOM(:)-VAL1*PHI1(:)
        PHIHOM(:)=PHIHOM(:)/PHIHOM(IRMATCH)
!       
!       == INHOMOGENEOUS SOLUTION WITH CORRECT BOUNDARY CONDITIONS =============
        IF(.NOT.THOM) THEN     
          CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,G,L,E,IDIR,PHIINHOM)
          CALL RADIAL$VALUE(GID,NR,PHIINHOM,RBOX,VAL1)
          CALL RADIAL$VALUE(GID,NR,PHI1,RBOX,VAL2)
          PHIINHOM(:)=PHIINHOM(:)-VAL1/VAL2*PHI1(:)
        ELSE
          PHIINHOM(:)=0.D0
        END IF
!
!       =======================================================================
!       ==  MATCH SOLUTION INSIDE AND OUTSIDE WITH VALUE                     ==
!       =======================================================================
        SVAR=(PHI(IRMATCH)-PHIINHOM(IRMATCH))/PHIHOM(IRMATCH)
        PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM(:)
        PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
!
!!$        CALL RADIAL$DERIVATIVE(GID,NR,PHI,R(IRMATCH),DER)
!!$        CALL RADIAL$DERIVATIVE(GID,NR,PHIINHOM,R(IRMATCH),DERO)
!!$        SVAR=(DERO-DER)/PHI(IRMATCH)
!!$        PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
      END IF
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB_EFINITENUCSIZE(GID,NR,RAD,AEZ,RHO,EFS)
!     **************************************************************************
!     **  ESTIMATE THE ENERGY CORRECTION DUE TO FINITE SIZE OF THE NUCLEUS    **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: RAD
      REAL(8)   ,INTENT(IN) :: RHO(NR) ! ELECTRON DENSITY
      REAL(8)   ,INTENT(OUT):: EFS
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)               :: AUX(NR),AUX1(NR)
      REAL(8)               :: R(NR)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$NUCPOT(GID,NR,AEZ,AUX)
      AUX=RHO*(R**2*AUX+AEZ*R/Y0)
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,EFS)
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$BOXVOFRHO(GID,NR,RAD,AEZ,RHO,POT,EH,EXC)
!     **************************************************************************
!     **  ELECTROSTATIC AND EXCHANGE-CORRELATION POTENTIAL            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: RAD
      REAL(8)   ,INTENT(IN) :: RHO(NR)
      REAL(8)   ,INTENT(OUT):: POT(NR)
      REAL(8)   ,INTENT(OUT):: EH
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,PARAMETER  :: RHOMIN=1.D-2
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)               :: POTH(NR)
      REAL(8)               :: POTXC(NR)
      REAL(8)               :: FOURPI
      REAL(8)               :: RHO1(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: EDEN(NR)
      REAL(8)               :: GRHO(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VGXC,VXC,EXC1,RH,GRHO2
      REAL(8)               :: DUMMY1,DUMMY2,DUMMY3
      REAL(8)               :: SVAR
      INTEGER(4)            :: IR,IRBOX
      REAL(8)               :: Q
!     **************************************************************************
      FOURPI=4.D0*PI
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRBOX=IR    ! SMALLEST GRID-POINT INDEX WITH R(IRBOX).GE.RAD
        IF(R(IR).GE.RAD) EXIT
      ENDDO
!
!     ==========================================================================
!     ==  TOTAL CHARGE                                                        ==
!     ==========================================================================
      AUX(:)=4.D0*PI*RHO(:)*R(:)**2*Y0
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,Q)
      Q=Q-AEZ
!
!     ==========================================================================
!     ==  TOTAL POTENTIAL                                                     ==
!     ==========================================================================
      RHO1=RHO
      IR=MIN(IRBOX+3,NR)
      RHO1(IR:)=0.D0        ! AVOID EXTREME NUMERICAL ROUNDING ERRORS
      CALL RADIAL$NUCPOT(GID,NR,AEZ,POTH)
      EDEN(:)=0.5D0*RHO1(:)*POTH(:)
      CALL RADIAL$POISSON(GID,NR,0,RHO1,AUX)
      POTH(:)=POTH(:)+AUX(:)
      CALL RADIAL$VALUE(GID,NR,POTH,RAD,SVAR)
!      SVAR=Q/RAD/Y0-SVAR
!      POTH(1:IRBOX-1)=POTH(1:IRBOX-1)+SVAR
!      POTH(IRBOX:NR)=Q/R(IRBOX:NR)/Y0
!CHANGE 090618 BOUNDARY CONDITIONS: HARD SPHERE IN A METAL.
      POTH(:)=POTH(:)-SVAR
      POTH(IRBOX:)=0.D0
      EDEN(:)=EDEN(:)+0.5D0*RHO1(:)*POTH(:)
!
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,EDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,EH)
!
!     ==========================================================================
!     == EXCHANGE CORRELATION                                                 ==
!     ==========================================================================
      CALL RADIAL$DERIVE(GID,NR,RHO,GRHO)
      DO IR=1,NR
        RH=RHO(IR)*Y0
        GRHO2=(Y0*GRHO(IR))**2
        CALL DFT(RH,0.D0,GRHO2,0.D0,0.D0,EXC1,VXC,DUMMY1,VGXC,DUMMY2,DUMMY3)
        EDEN(IR)=4.D0*PI*EXC1   ! ANGULAR INTEGRATION ALREADY INCLUDED
        POTXC(IR)=VXC/Y0
        GRHO(IR)=VGXC*2.D0*GRHO(IR)
      ENDDO
      AUX(:)=R(:)**2*GRHO(:)
      CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
      GRHO(2:)=AUX1(2:)/R(2:)**2
      GRHO(1:5)=GRHO(5) ! AVOID ERRORS DUE TO TERMINATION OF THE GRID
                        ! 5 POINTS OFFSET FOR 5-POINT FORMULA APPLIED TWICE...
      POTXC(:)=POTXC(:)-GRHO(:)
!
!     == EXCHANGE CORRELATION POTENTIAL IS SET TO ZERO OUTSIDE THE BOX. 
!     == OTHERWISE IT WILL HAVE A CONSTANT VALUE IN THE ZERO-DENSITY REGION
      POTXC(IRBOX:)=0.D0
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,EDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,EXC)
!
!     ==========================================================================
!     ==  CUT OF POTENTIAL FOR LOW DENSITIES                                  ==
!     ==========================================================================
      DO IR=1,NR
        IF(RHO(IR)*Y0.LT.1.D-6) POTXC(IR)=0.D0
      ENDDO
!
!     ==========================================================================
!     ==  CUT OF POTENTIAL FOTR LOW DENSITIES                                 ==
!     ==========================================================================
      POT=POTH+POTXC
!CALL ATOMLIB_WRITEPHI('POTH',GID,NR,1,POTH)
!CALL ATOMLIB_WRITEPHI('POTXC',GID,NR,1,POTXC)
!CALL ATOMLIB_WRITEPHI('RHO',GID,NR,1,RHO)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$BOXMUX(GID,NR,RAD,RHO,MUX)
!     ******************************************************************
!     **                                                              **
!     **  ELECTROSTATIC AND EXCHANGE-CORRELATION POTENTIAL            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8),   INTENT(IN) :: RAD
      REAL(8),   INTENT(IN) :: RHO(NR)
      REAL(8),   INTENT(OUT):: MUX(NR)
      REAL(8),   PARAMETER  :: RHOMIN=1.D-2
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)               :: FOURPI
      REAL(8)               :: AUX(NR),AUX1(NR)
      REAL(8)               :: GRHO(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VGXC,VXC,EXC1,RH,GRHO2
      REAL(8)               :: DUMMY1,DUMMY2,DUMMY3
      INTEGER(4)            :: IR,IRBOX
!     **************************************************************************
      FOURPI=4.D0*PI
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GE.RAD) EXIT
      ENDDO
      CALL DFT$SETL4('XCONLY',.TRUE.)
!
!     ==========================================================================
!     ==  TOTAL POTENTIAL                                                     ==
!     ==========================================================================
      CALL RADIAL$DERIVE(GID,NR,RHO,GRHO)
      DO IR=1,NR
        RH=RHO(IR)*Y0
        GRHO2=(Y0*GRHO(IR))**2
        CALL DFT(RH,0.D0,GRHO2,0.D0,0.D0,EXC1,VXC,DUMMY1,VGXC,DUMMY2,DUMMY3)
        MUX(IR)=VXC/Y0
        GRHO(IR)=VGXC*2.D0*GRHO(IR)
      ENDDO
      AUX(:)=R(:)**2*GRHO(:)
      CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
      GRHO(2:)=AUX1(2:)/R(2:)**2
      GRHO(1:5)=GRHO(5) ! AVOID ERRORS DUE TO TERMINATION OF THE GRID
                        ! 5 POINTS OFFSET FOR 5-POINT FORMULA APPLIED TWICE...
      MUX(:)=MUX(:)-GRHO
!
!     == EXCHANGE CORRELATION POTENTIAL IS SET TO ZERO OUTSIDE THE BOX. 
!     == OTHERWISE IT WILL HAVE A CONSTANT VALUE IN THE ZERO-DENSITY REGION
      MUX(IRBOX:)=0.D0
!
!     ==========================================================================
!     ==  CUT OF POTENTIAL FOR LOW DENSITIES                                  ==
!     ==========================================================================
      DO IR=1,NR
        IF(RHO(IR)*Y0.LT.1.D-6) MUX(IR)=0.D0
      ENDDO
      CALL DFT$SETL4('XCONLY',.FALSE.)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$NODELESS(GID,NR,RBOX,POT,NB,LOFI,EOFI,UOFI,TUOFI)
!     **************************************************************************
!     ** CALCULATES THE SEQUENCE OF NODELESS WAVE FUNCTIONS                   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: GID
      INTEGER(4),INTENT(IN)   :: NR
      REAL(8)   ,INTENT(IN)   :: RBOX
      REAL(8)   ,INTENT(IN)   :: POT(NR)
      INTEGER(4),INTENT(IN)   :: NB
      INTEGER(4),INTENT(IN)   :: LOFI(NB)
      REAL(8)   ,INTENT(INOUT):: EOFI(NB)
      REAL(8)   ,INTENT(OUT)  :: UOFI(NR,NB)
      REAL(8)   ,INTENT(OUT)  :: TUOFI(NR,NB)
      REAL(8)   ,PARAMETER    :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER    :: Y0=1.D0/SQRT(4.D0*PI)
      INTEGER(4)              :: L
      REAL(8)                 :: E
      REAL(8)                 :: G(NR)
      REAL(8)                 :: DREL(NR)
      INTEGER(4)              :: NN
      INTEGER(4)              :: SO
      INTEGER(4)              :: IB,I
      LOGICAL(4)              :: TVARDREL
!     **************************************************************************
      UOFI(:,:)=0.D0
      TUOFI(:,:)=0.D0
      DO IB=1,NB
        L=LOFI(IB)
        E=EOFI(IB)
        G(:)=0.D0
        DO I=IB-1,1,-1
          IF(LOFI(I).EQ.L) THEN
            G=UOFI(:,I)
            EXIT
          END IF
        ENDDO
        SO=0
        DREL(:)=0.D0
        TVARDREL=.FALSE.
        NN=0
        CALL ATOMLIB$BOUNDSTATE(GID,NR,L,SO,0.D0,RBOX,TVARDREL,DREL,G,NN,POT &
     &                         ,E,UOFI(:,IB))
        TUOFI(:,IB)=G(:)+(E-POT(:)*Y0)*UOFI(:,IB)
        EOFI(IB)=E
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO,DH,DO,G,PHI)
!     **************************************************************************
!     **                                                                      **
!     **  SOLVES THE RADIAL PAW -SCHROEDINGER EQUATION.                       **
!     **    (T+VTILDE-E+|P>(DH-E*DO<P|]|PHI>=|G>                              **
!     **  WHERE T IS THE NONRELATIVISTIC KINETIC ENERGY.                      **
!     **                                                                      **
!     **    DH=<AEPHI|T+V|AEPHI>-<PSPHI|T+VTILDE|PSPHI>                       **
!     **    DO=<AEPHI|AEPHI>-<PSPHI|PSPHI>                                    **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID     !GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR      ! #(RADIAL GRID POINTS
      INTEGER(4)  ,INTENT(IN) :: L       ! ANGULAR MOMENTUM
      REAL(8)     ,INTENT(IN) :: E       ! ENERGY
      REAL(8)     ,INTENT(IN) :: PSPOT(NR)  !VTILDE=PSPOT*Y0
      INTEGER(4)  ,INTENT(IN) :: NPRO       ! #(PROJECTOR FUNCTIONS)
      REAL(8)     ,INTENT(IN) :: PRO(NR,NPRO) ! PROJECTOR FUNCTIONS
      REAL(8)     ,INTENT(IN) :: DH(NPRO,NPRO)     
      REAL(8)     ,INTENT(IN) :: DO(NPRO,NPRO)     
      REAL(8)     ,INTENT(IN) :: G(NR)        !INHOMOGENEITY
      REAL(8)     ,INTENT(OUT):: PHI(NR)      !PAW RADIAL FUNCTION
      REAL(8)                 :: U(NR)
      REAL(8)                 :: V(NR,NPRO)
      REAL(8)                 :: AMAT(NPRO,NPRO)
      REAL(8)                 :: BMAT(NPRO,NPRO)
      REAL(8)                 :: BMATINV(NPRO,NPRO)
      REAL(8)                 :: CMAT(NPRO,NPRO)
      REAL(8)                 :: CVEC(NPRO)
      REAL(8)                 :: DVEC(NPRO)
      INTEGER(4)              :: I1,I2,I3
      REAL(8)                 :: AUX(NR)
      REAL(8)                 :: R(NR)
      REAL(8)                 :: DREL(NR)
      INTEGER(4)              :: SO
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  -1/2NABLA^2+POT-E|U>=|G>                                            ==
!     ==========================================================================
      SO=0
      DREL(:)=0.D0
      CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,G,L,E,1,U)
!
!     ==========================================================================
!     ==  -1/2NABLA^2+POT-E|V>=|PRO>                                          ==
!     ==========================================================================
      DO I1=1,NPRO
        CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,PRO(:,I1),L,E,1,V(:,I1))
      ENDDO
!!$PRINT*,'E ',E,SO,L,NR,GID
!!$OPEN(100,FILE='XXX.DAT')
!!$ DO IR=1,NR
!!$    WRITE(100,FMT='(22F30.10)')R(IR),U(IR),V(IR,:),PRO(IR,:)
!!$ ENDDO
!!$CLOSE(10)
!!$STOP 'IN PAWDER'
!
!     ==========================================================================
!     ==  AMAT=<PRO|V>  CVEC=<PRO|U>                                          ==
!     ==========================================================================
      DO I1=1,NPRO
        DO I2=1,NPRO
          AUX(:)=PRO(:,I1)*V(:,I2)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,AMAT(I1,I2))
        ENDDO
      ENDDO
      DO I1=1,NPRO
        AUX(:)=PRO(:,I1)*U(:)*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,CVEC(I1))
      ENDDO  
!
!     ==========================================================================
!     ==  BMAT=1+(DATH-EDO)<PRO|V>                                            ==
!     ==========================================================================
!     BMAT(:,:)=MATMUL(DH(:,:)-E*DO(:,:),AMAT(:,:))
!     DO I1=1,NPRO
!       BMAT(I1,I1)=BMAT(I1,I1)+1.D0
!     ENDDO
!
      DO I1=1,NPRO
        DO I2=1,NPRO
          BMAT(I1,I2)=0.D0
          DO I3=1,NPRO
            BMAT(I1,I2)=BMAT(I1,I2)+(DH(I1,I3)-E*DO(I1,I3))*AMAT(I3,I2)     
          ENDDO
        ENDDO
        BMAT(I1,I1)=BMAT(I1,I1)+1.D0
      ENDDO
!
!     ==========================================================================
!     ==  BMAT = BMAT^-1 = [1+(DATH-EDO)<PRO|V>]^-1                           ==
!     ==========================================================================
      IF(NPRO.EQ.0) THEN
        CALL ERROR$STOP('NPRO=0 NOT ALLOWED')
      END IF
      IF(NPRO.EQ.1) THEN 
        BMAT(1,1)=1.D0/BMAT(1,1)
      ELSE 
        CALL LIB$INVERTR8(NPRO,BMAT,BMATINV)
        BMAT=BMATINV
      END IF
!
!     ==========================================================================
!     ==  CMAT = -BMAT*(DATH-EDO)                                             ==
!     ==       = -[1+(DATH-EDO)*<PRO|V>]^-1 (DATH-EDO)                        ==
!     ==========================================================================
!     CMAT(:,:)=MATMUL(BMAT(:,:),DH(:,:)-E*DO(:,:))
      DO I1=1,NPRO
        DO I2=1,NPRO
          CMAT(I1,I2)=0.D0
          DO I3=1,NPRO
            CMAT(I1,I2)=CMAT(I1,I2)-BMAT(I1,I3)*(DH(I3,I2)-E*DO(I3,I2))
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  DVEC = CMAT*CVEC                                                    ==
!     ==       = -[1+(DATH-EDO)*<PRO|V>]^-1 (DATH-EDO) <PRO|U>                ==
!     ==========================================================================
!     DVEC(:)=MATMUL(CMAT(:,:),CVEC(:))
      DO I1=1,NPRO
        DVEC(I1)=0.D0
        DO I2=1,NPRO
          DVEC(I1)=DVEC(I1)+CMAT(I1,I2)*CVEC(I2)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  |PHI> = |U>+|V>DVEC                                                 ==
!     ==  = [1-|V>[1+(DATH-EDO)*<PRO|V>]^-1(DATH-EDO)<PRO|] |U>               ==
!     ==========================================================================
!     PHI(:)=U(:)+MATMUL(V(:,:),DVEC(:))
      PHI(:)=U(:)
      DO I1=1,NPRO
        PHI(:)=PHI(:)+V(:,I1)*DVEC(I1)
      ENDDO
      RETURN
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB_ORTHODER(GID,NR,L,E,PSPOT,NPRO,PRO,G,PHI)
!     **************************************************************************
!     **  SOLVES THE SCHROEDINGER EQUATION IN A SUBSPACE ORTHOGONAL TO THE    **
!     **  PROJECTOR FUNCTIONS.                                                **
!     **                                                                      **
!     **  1) THE SOLUTION PHI IS ORTHOGONAL TO ALL PROJECTOR FUNCTIONS.       **
!     **  2) THE RESIDUAL OF THE DIFFERENTIAL EQUATION LIES ENTIRELY IN THE   **
!     **     OF PROJECTOR FUNCTIONS                                           **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID     !GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR      ! #(RADIAL GRID POINTS
      INTEGER(4)  ,INTENT(IN) :: L       ! ANGULAR MOMENTUM
      REAL(8)     ,INTENT(IN) :: E       ! ENERGY
      REAL(8)     ,INTENT(IN) :: PSPOT(NR)  !VTILDE=PSPOT*Y0
      INTEGER(4)  ,INTENT(IN) :: NPRO       ! #(PROJECTOR FUNCTIONS)
      REAL(8)     ,INTENT(IN) :: PRO(NR,NPRO) ! PROJECTOR FUNCTIONS
      REAL(8)     ,INTENT(IN) :: G(NR)        !INHOMOGENEITY
      REAL(8)     ,INTENT(OUT):: PHI(NR)      !PAW RADIAL FUNCTION
      LOGICAL(4)  ,PARAMETER  :: TTEST=.TRUE.
      REAL(8)                 :: U(NR)
      REAL(8)                 :: V(NR,NPRO)
      REAL(8)                 :: AMAT(NPRO,NPRO)
      REAL(8)                 :: CVEC(NPRO)
      REAL(8)                 :: DVEC(NPRO)
      INTEGER(4)              :: I,J
      REAL(8)                 :: AUX(NR)
      REAL(8)                 :: SVAR
      REAL(8)                 :: R(NR)
      REAL(8)                 :: DREL(NR)
      INTEGER(4)              :: SO
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  -1/2NABLA^2+POT-E|U>=|G>                                            ==
!     ==========================================================================
      SO=0
      DREL(:)=0.D0
      CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,G,L,E,1,U)
!!$      SVAR=1.D0/U(NR-3)
!!$      U=U*SVAR
!
!     ==========================================================================
!     ==  -1/2NABLA^2+POT-E|V>=|PRO>                                          ==
!     ==========================================================================
      DO I=1,NPRO
        CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,PRO(:,I),L,E,1,V(:,I))
!       __AVOID EXPONENTIALLY INCREASING SOLUTION TO SOME EXTENT________________
        SVAR=V(NR-3,I)/U(NR-3)
        V(:,I)=V(:,I)-U(:)*SVAR
      ENDDO
!
!     ==========================================================================
!     ==  AMAT=<PRO|V>  CVEC=<PRO|U>                                          ==
!     ==========================================================================
      DO I=1,NPRO
        DO J=1,NPRO
          AUX(:)=PRO(:,I)*V(:,J)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,AMAT(I,J))
        ENDDO
      ENDDO
      DO I=1,NPRO
        AUX(:)=PRO(:,I)*U(:)*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,CVEC(I))
      ENDDO  
!
!     ==========================================================================
!     == DETERMINE COEFFICIENTS THAT IMPOSE ORTHOGONALITY TO PROJECTOR FUNCS  ==
!     ==========================================================================
!     == THIS EQUATION FAILS, IF A SINGULAR VALUE OF AMAT VANISHES. IN THIS   ==
!     == THE SOLUTION IS GIVEN ENTIRELY BY THE V=FUNCTIONS WITH COEFFICIENTS  ==
!     == FROM THE SINGULAR VECTOR AND THE CONTRIBUTION OF U VANISHES. ==========
!     == THIS POSSIBILITY IS NOT CONSIDERED YET. ===============================
!     ==========================================================================
      IF(NPRO.EQ.0) THEN 
      ELSE IF(NPRO.EQ.1) THEN 
        DVEC(1)=-CVEC(1)/AMAT(1,1)
      ELSE 
!       == AMAT*DVEC+C=0 =======================================================
        CALL LIB$MATRIXSOLVER8(NPRO,NPRO,1,AMAT,DVEC,-CVEC)
      END IF
!
!     ==========================================================================
!     ==  |PHI> = |U>+|V>DVEC                                                 ==
!     ==========================================================================
!     PHI(:)=U(:)+MATMUL(V(:,:),DVEC(:))
      PHI(:)=U(:)
      DO I=1,NPRO
        PHI(:)=PHI(:)+V(:,I)*DVEC(I)
      ENDDO
!
!     ==========================================================================
!     ==  TEST                                                                ==
!     ==========================================================================
      IF(TTEST.AND.NPRO.GT.0) THEN
!       == PARK RESIDUAL IN DVEC ===============================================
        DVEC=MATMUL(AMAT,DVEC)+CVEC
        DO I=1,NPRO
          AUX=PRO(:,I)*PHI(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,CVEC(I))
          AUX=PRO(:,I)**2*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          CVEC(I)=CVEC(I)/SQRT(SVAR)
        ENDDO
        AUX=R**2*PHI**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        CVEC=CVEC/SQRT(SVAR)

        IF(MAXVAL(ABS(CVEC)).GT.1.D-4) THEN
          PRINT*,'TEST MATRIX EQ ',DVEC
          PRINT*,'TEST FINAL     ',CVEC
CALL ATOMLIB_WRITEPHI('CRASH_PHI',GID,NR,1,PHI)
CALL ATOMLIB_WRITEPHI('CRASH_PRO',GID,NR,NPRO,PRO)
CALL ATOMLIB_WRITEPHI('CRASH_U',GID,NR,1,U)
CALL ATOMLIB_WRITEPHI('CRASH_V',GID,NR,NPRO,V)
          CALL ERROR$MSG('ORTHOGONALITY TEST FAILED')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$R8VAL('E',E)
          CALL ERROR$I4VAL('NPRO',NPRO)
          CALL ERROR$R8VAL('MAX COS',MAXVAL(ABS(CVEC)))
          CALL ERROR$STOP('ATOMLIB$ORTHODER')
        END IF
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$UPDATESTATEWITHHF(GID,NR,L,SO,DREL,G,POT,VFOCK &
     &                                    ,RBND,E,PSI)
!     **************************************************************************
!     ** USES AN APPROXIMATE SOLUTION TO CONSTRUCT A SOLUTION FOR             **
!     ** A HAMILTONIAN WITH A HARTEE-FOCK CONTRIBUTION MIXED IN               **
!     **                                                                      **
!     ** THE BOUNDARY CONDITIONS (LOGARITHMIC DERIVATIVE AT RBND IS PRESERVED **
!     ** IF RBND.LE.0 A FIXED ENERGY SOLUTION IS SEARCHED AND THE OUTER       **
!     ** BOUNDARY CONDITION IS RELAXED.                                       **
!     **                                                                      **
!     ** TYPE CAN BE:                                                         **
!     **   'BOUND'                                                            **
!     **   'PHASE'    LOGARITHMIC DERIVATIVE AT RBND IS FIXED                 **
!     **   'ENERGY'   ENERGY IS FIXED TO INPUT VALUE, RBND IS NOT USED        **
!     **                                                                      **
!     ** REMARK: HAS CONVERGENCE PROBLEMS IF RBND IS CHOSEN VERY LARGE,       **
!     **    PROBABLY RELATED TO LOW-LYING STATES                              **
!     **                                                                      **
!     **************************************************************************
      USE RADIALFOCK_MODULE, ONLY : VFOCK_TYPE
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(IN) :: GID
      INTEGER(4)      ,INTENT(IN) :: NR
      INTEGER(4)      ,INTENT(IN) :: L
      INTEGER(4)      ,INTENT(IN) :: SO
      REAL(8)         ,INTENT(IN) :: DREL(NR)
      REAL(8)         ,INTENT(IN) :: POT(NR)
      TYPE(VFOCK_TYPE),INTENT(IN) :: VFOCK
      REAL(8)         ,INTENT(IN) :: G(NR)
      REAL(8)         ,INTENT(IN) :: RBND
      REAL(8)         ,INTENT(INOUT) :: E
      REAL(8)         ,INTENT(INOUT) :: PSI(NR)
      REAL(8)         ,PARAMETER  :: TOL=1.D-7 !(1D-12 IS TOO SMALL)
      INTEGER(4)      ,PARAMETER  :: NITER=200
      REAL(8)         ,PARAMETER  :: XMAX=1.D+15 !MAX. FACTOR FOR WAVEFUNCTION
      LOGICAL(4)      ,PARAMETER  :: TPR=.FALSE.
      REAL(8)         ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)         ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                  :: PHASE
      REAL(8)                  :: DPSI(NR)
      REAL(8)                  :: DEVARR(NR)
      REAL(8)                  :: PHIDOT(NR)
      REAL(8)                  :: PHIPRIME(NR)
      REAL(8)                  :: POT1(NR)
      REAL(8)                  :: VALDOT,DERDOT,VAL0,DER0
      REAL(8)                  :: W0,WDOT
      REAL(8)                  :: Y1,Y2,X1,X2
      REAL(8)                  :: SOFACTOR
      REAL(8)                  :: R(NR)
      REAL(8)                  :: AUX(NR),AUX1(NR),AUX2(NR)
      REAL(8)                  :: G1(NR)
      REAL(8)                  :: RDPRIME(NR)
      REAL(8)                  :: VAL,DER,SVAR
      INTEGER(4)               :: IR,I
      INTEGER(4)               :: ITER
      LOGICAL(4)               :: CONVG
      INTEGER(4)               :: IRBND
      INTEGER(4)               :: IROUT
      INTEGER(4)               :: IRCL
      INTEGER(4)               :: IREND
      REAL(8)                  :: ROUT
      REAL(8)                  :: REND
      REAL(8)                  :: NORM
      LOGICAL(4)               :: TFIXE ! CONSTANT ENERGY CALCULATION IF TRUE
      LOGICAL(4)               :: TINHOM ! INHOMOGENEOUS EQUATION
      INTEGER(4)               :: NN
!CHARACTER(64)            :: STRING
!REAL(8)                  :: TESTARR(NR,NITER,6)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
      TFIXE=RBND.LE.0.D0
      TINHOM=MAXVAL(ABS(G)).GT.0.D0
!
!     :: REMARK: SCHROEDINGER$PHASESHIFT TAKES VALUE AND DERIVATIVE FROM A    ::
!     :: TWO-POINT FORMULA, AND THEREFORE IS NOT FULLY CONSISTENT WITH        ::
!     :: RADIAL$VALUE AND RADIAL$DERIVATIVE                                   ::
      CALL SCHROEDINGER$PHASESHIFT(GID,NR,PSI,0.D0,RBND,PHASE)
!
      DO IR=1,NR
        IRBND=IR
        IF(TFIXE.AND.R(IR).GT.3.D0) EXIT
!       == NOTE! IRBND MUST BE CONSISTENT WITH SCHROEDINGER$PHASESHIFT,  =======
!       == THUS R(IR).GE.RBND AND NOT ".GT." ===================================
        IF(.NOT.TFIXE.AND.R(IR).GE.RBND) EXIT
      ENDDO

      IF(SO.EQ.0) THEN
        SOFACTOR=0.D0
      ELSE IF(SO.EQ.1) THEN
        SOFACTOR=REAL(L,KIND=8)
      ELSE IF(SO.EQ.-1) THEN
        SOFACTOR=REAL(-L-1,KIND=8)
      ELSE
        CALL ERROR$STOP('ILLEGAL VALUE FOR SO')
        CALL ERROR$STOP('ATOMLIB$UPDATESTATEWITHHF')
      END IF
!CALL ATOMLIB_WRITEPHI('PSIBEFORE',GID,NR,1,PSI)
!
!     ==========================================================================
!     ==  NOW START LOOP                                                      ==
!     ==========================================================================
      DO ITER=1,NITER
!
!       ========================================================================
!       ==  CONSTRUCT |DPSI>=DE/DPSI                                          ==
!       ========================================================================
!       == KINETIC ENERGY OF THE WAVE FUNCTION (NONRELATIVISTIC)
        CALL RADIAL$VERLETD2(GID,NR,R*PSI,AUX1)
        AUX(:)=R(:)*AUX1(:)
        AUX(:)=AUX(:)-REAL(L*(L+1),KIND=8)*PSI(:)
        AUX(:)=-(1.D0+DREL)*AUX(:)
        CALL RADIAL$VERLETD1(GID,NR,PSI,AUX1)
        AUX(:)=AUX(:)-RDPRIME(:)*(R(:)**2*AUX1(:)-SOFACTOR*R(:)*PSI(:))
!       == POTENTIAL ENERGY ====================================================
        AUX(:)=AUX(:)+2.D0*R(:)**2*POT(:)*Y0*PSI(:)
!       == REMOVE FACTOR 2*R**2 ================================================
        AUX(:)=AUX(:)/(2.D0*R(:)**2)
!       == ADD FOCK TERM  ======================================================
        CALL RADIALFOCK$VPSI(GID,NR,VFOCK,L,PSI,AUX1)
        AUX(:)=AUX(:)+AUX1(:)
!       == SUBTRACT INHOMOGENEITY ==============================================
        AUX(:)=AUX(:)-G(:)
!       == CORRECT GLITCHES ====================================================
        AUX(1:2)=0.D0
!       == ESTIMATE ENERGY =====================================================
        IF(.NOT.TFIXE) THEN
          AUX1(:)=R(:)**2*PSI(:)*AUX(:)
          CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
          CALL RADIAL$VALUE(GID,NR,AUX2,RBND,E)
          AUX1(:)=R(:)**2*PSI(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
          CALL RADIAL$VALUE(GID,NR,AUX2,RBND,NORM)
          E=E/NORM
          IF(E.LT.-6.D+3) THEN
            CALL ERROR$MSG('ENERGY BELOW -6000 H ENCOUNTERED.')
            CALL ERROR$MSG('INDICATION OF CATASTROPHIC FAILURE.') 
            CALL ERROR$R8VAL('E',E)
            CALL ERROR$R8VAL('RBND',RBND)
            CALL ERROR$STOP('ATOMLIB$UPDATESTATEWITHHF') 
          END IF
!         == REMOVE ENERGY TERM ================================================
          AUX(:)=AUX(:)-E*PSI(:)
        ELSE
!         == REMOVE ENERGY TERM ================================================
          AUX(:)=AUX(:)-E*PSI(:)
!         = A REASONABLE NORM IS NEEDED LATER FOR JUDGING THE CONVERGENCE ======
          AUX1(:)=R(:)**2*PSI(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
          CALL RADIAL$VALUE(GID,NR,AUX2,3.D0,NORM)
        END IF
!!$IF(TINHOM.AND.(.NOT.TFIXE)) THEN
!!$  PRINT*,'L ',L,' E ',E
!!$  CALL ATOMLIB_WRITEPHI('X1.DAT',GID,NR,1,AUX)
!!$  CALL ATOMLIB_WRITEPHI('X2.DAT',GID,NR,1,PSI)
!!$END IF
!
!       == MAP INTO DPSI =======================================================
        DPSI(:)=-AUX(:)
        DPSI(1:2)=0.D0
        DEVARR=DPSI
!
!       ========================================================================
!       ==  DETERMINE PHIPRIME-PHIDOT*EPRIME                                  ==
!       ========================================================================
        CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT,E,IRCL,IROUT)
        ROUT=R(IROUT)
!       ==  SET KINETIC ENERGY TO ZERO BEYOND ROUT TO AVOID AN OVERFLOW ========
!       == ATTENTION: THERE IS A STEP IN THE POTENTIAL =========================
        POT1=POT
        POT1(IROUT:)=E/Y0
!
        IF(TFIXE) THEN
          G1=DPSI(:)
          G1(1:2)=0.D0
          G1(IROUT:)=0.D0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G1,L,E,1,PHIPRIME)
          DPSI(:)=PHIPRIME(:)
        ELSE
          G1=PSI(:)
          G1(IROUT:)=0.D0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G1,L,E,1,PHIDOT)
!         == ADJUST ENERGY UNTIL PHIPRIME HAS A NODE AT RBND. AVOIDS ===========
!         == PROBLEMS WITH EXPONENTIALLY INCREASING TAIL FOR CORE STATES =======
          SVAR=0.D0
          DO I=1,2
            G1=DPSI(:)+SVAR*PSI(:)
            G1(1:2)=0.D0
            G1(IROUT:)=0.D0
            CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G1,L,E,1,PHIPRIME)
            SVAR=SVAR-PHIPRIME(IRBND)/PHIDOT(IRBND)
          ENDDO
!
!         == DETERMINE VALUES AND DERIVATIVES CONSISTENT WITH ==================
!         == SCHROEDINGER$PHASESHIFT ===========================================
          CALL SCHROEDINGER$RESOLVEPHASESHIFT(PHASE,VAL,DER,NN)
          IREND=IRBND
          REND=RBND
          X1=R(IREND-1)
          X2=R(IREND)
          Y1=PHIDOT(IREND-1)
          Y2=PHIDOT(IREND)
          DERDOT=(Y2-Y1)/(X2-X1)
          VALDOT=Y1+DERDOT*(REND-X1)
          WDOT=VAL*DERDOT-DER*VALDOT
!
          DPSI(:)=PHIPRIME(:)
          DO I=1,4 ! LOOP TO FIX UP NUMERICAL ERRORS 
            Y1=PSI(IREND-1)+DPSI(IREND-1)
            Y2=PSI(IREND)+DPSI(IREND)
            DER0=(Y2-Y1)/(X2-X1)
            VAL0=Y1+DER0*(REND-X1)
            W0=VAL*DER0-DER*VAL0
            SVAR=-W0/WDOT  !DELTA EPSILON
            DPSI(:)=DPSI(:)+PHIDOT(:)*SVAR
          ENDDO
          IF(REND.LT.RBND) DPSI(IREND+1:)=0.D0
        END IF
!TESTARR(:,ITER,6)=DPSI
!   
!       ========================================================================
!       ==  CHECK CONVERGENCE                                                 ==
!       ========================================================================
        VAL=MAXVAL(ABS(DPSI(:IRBND-3)))/SQRT(NORM)
        CONVG=(VAL.LT.TOL)
        IF(CONVG.AND.VAL.LT.0.D0) THEN
          CALL ATOMLIB_WRITEPHI('DPSI',GID,NR,1,DPSI)
          CALL ERROR$MSG('NOT-A-NUMBER ENCOUNTERED')
          CALL ERROR$R8VAL('NORM',NORM)
          CALL ERROR$R8VAL('MAX DEVIATION',VAL)
          CALL ERROR$STOP('ATOMLIB$UPDATESTATEWITHHF')
        END IF
        IF(CONVG) EXIT
!
!       ========================================================================
!       ==  PROPAGATE                                                         ==
!       ========================================================================
        PSI(:)=PSI(:)+DPSI(:)
!
      ENDDO    ! END OF ITERATION
!CALL ATOMLIB_WRITEPHI('PSIAFTER',GID,NR,1,PSI)
!IF(TINHOM)STOP 'FORCED ATOMLIB$UPDATESTATEWITHHF'

      IF(TPR) PRINT*,'ITER=',ITER,' L=',L,' E=',E,' DEV=',VAL
!
      IF(.NOT.CONVG) THEN
        CALL ATOMLIB_WRITEPHI('PSI.DAT',GID,NR,1,PSI/SQRT(NORM))
        CALL ATOMLIB_WRITEPHI('DPSI.DAT',GID,NR,1,DPSI/SQRT(NORM))
        CALL ATOMLIB_WRITEPHI('DEVARR.DAT',GID,NR,1,DEVARR/SQRT(NORM))
!!$        DO ITER=1,NITER,10
!!$          WRITE(STRING,*)ITER
!!$          STRING='TESTARR1_'//TRIM(ADJUSTL(STRING))//'.DAT'
!!$          CALL ATOMLIB_WRITEPHI(STRING,GID,NR,10,TESTARR(:,ITER:ITER+9,1))
!!$          STRING(8:8)='2'
!!$          CALL ATOMLIB_WRITEPHI(STRING,GID,NR,10,TESTARR(:,ITER:ITER+9,2))
!!$          STRING(8:8)='3'
!!$          CALL ATOMLIB_WRITEPHI(STRING,GID,NR,10,TESTARR(:,ITER:ITER+9,3))
!!$          STRING(8:8)='4'
!!$          CALL ATOMLIB_WRITEPHI(STRING,GID,NR,10,TESTARR(:,ITER:ITER+9,4))
!!$          STRING(8:8)='5'
!!$          CALL ATOMLIB_WRITEPHI(STRING,GID,NR,10,TESTARR(:,ITER:ITER+9,5))
!!$          STRING(8:8)='6'
!!$          CALL ATOMLIB_WRITEPHI(STRING,GID,NR,10,TESTARR(:,ITER:ITER+9,6))
!!$        ENDDO
        CALL ERROR$MSG('LOOP FOR RADIAL HARTREE FOCK NOT CONVERGED')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('E',E)
        CALL ERROR$R8VAL('REND',REND)
        CALL ERROR$R8VAL('RBND',RBND)
        CALL ERROR$R8VAL('NORM',NORM)
        CALL ERROR$I4VAL('#(ITERATIONS)',ITER)
        CALL ERROR$R8VAL('MAX DEVIATION',VAL)
        CALL ERROR$R8VAL('TOLERANCE',TOL)
        CALL ERROR$STOP('ATOMLIB$UPDATESTATEWITHHF')
      END IF
      RETURN
      END SUBROUTINE ATOMLIB$UPDATESTATEWITHHF
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$ETOTWITHFOCK(GID,NR,RBOX,LRHOX,TREL,POT,VFOCK &
    &                                ,NB,LOFI,SOFI,FOFI,PSI,ETOT)
!     **************************************************************************
!     **************************************************************************
      USE RADIALFOCK_MODULE, ONLY : VFOCK_TYPE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: RBOX
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: LOFI(NB)
      INTEGER(4),INTENT(IN) :: SOFI(NB)
      LOGICAL(4),INTENT(IN) :: TREL
      INTEGER(4),INTENT(IN) :: LRHOX
      REAL(8)   ,INTENT(IN) :: FOFI(NB)
      REAL(8)   ,INTENT(IN) :: PSI(NR,NB)
      REAL(8)   ,INTENT(IN) :: POT(NR)
      TYPE(VFOCK_TYPE),INTENT(IN) :: VFOCK
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)               :: DRRPSI(NR,NB)
      REAL(8)               :: AUX(NR),AUX1(NR),AUX2(NR)
      REAL(8)               :: EKDEN(NR),EHDEN(NR),EXDEN(NR),EKIN,EH,EX
      REAL(8)               :: R(NR)
      REAL(8)               :: EOFI(NB)
      REAL(8)               :: DREL(NR),RDPRIME(NR)
      REAL(8)               :: VAL,DER
      INTEGER(4)            :: IRBOX
      REAL(8)               :: SOFACTOR
      INTEGER(4)            :: IR,IB
      INTEGER(4)            :: L
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     == INDEX OF THE FIRST GRIDPOINT OUTSIDE RBOX =============================
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GE.RBOX) EXIT
      ENDDO
!
!     ==========================================================================
!     ==  VFOCK|PSI> AND PARTIAL_R*R|\PSI>                                    ==
!     ==========================================================================
      DO IB=1,NB
        AUX(:)=R(:)*PSI(:,IB)
        VAL=AUX(IRBOX)
        DER=(AUX(IRBOX+1)-AUX(IRBOX-1))/(R(IRBOX+1)-R(IRBOX-1))
        AUX(IRBOX:)=VAL+(R(IRBOX:)-R(IRBOX))*DER
        CALL RADIAL$VERLETD1(GID,NR,AUX,AUX1)
        AUX1(1:2)=AUX1(3)  
        DRRPSI(:,IB)=AUX1(:)
      ENDDO
!
!     ==========================================================================
!     ==  ESTIMATE ENERGIES FOR RELATIVISTIC CORRECTIONS                      ==
!     ==========================================================================
      DO IB=1,NB
        EOFI(IB)=0.D0
      ENDDO
!
!     ========================================================================
!     ==  CALCULATE TOTAL ENERGY (THIS A FAKE QUANTITY)                     ==
!     ========================================================================
      AUX(:)=0.D0
      EKDEN(:)=0.D0
      EHDEN(:)=0.D0
      EXDEN(:)=0.D0
      DO IB=1,NB
        L=LOFI(IB)
        IF(TREL) THEN
          CALL SCHROEDINGER$DREL(GID,NR,POT,EOFI(IB),DREL)
          CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
        ELSE
          DREL=0.D0
          RDPRIME=0.D0
        END IF
        IF(SOFI(IB).EQ.0) THEN
          SOFACTOR=0.D0
        ELSE IF(SOFI(IB).EQ.1) THEN
          SOFACTOR=REAL(L,KIND=8)
        ELSE IF(SOFI(IB).EQ.-1) THEN
          SOFACTOR=REAL(-L-1,KIND=8)
        ELSE
          CALL ERROR$STOP('ILLEGAL VALUE FOR SOFI ')
          CALL ERROR$STOP('ATOMLIB$AESCFWITHHF')
        END IF
!       == KINETIC ENERGY ====================================================
        AUX1(:)=0.5D0*(1.D0+DREL) &
     &               *(DRRPSI(:,IB)**2+REAL(L*(L+1),KIND=8)*PSI(:,IB)**2)
        AUX1(:)=AUX1(:)+0.5D0*(SOFACTOR+1.D0)*RDPRIME*R(:)*PSI(:,IB)**2
        EKDEN(:)=EKDEN(:)+FOFI(IB)*AUX1(:)
!       == POTENTIAL ENERGY ==================================================
        AUX1(:)=0.D0
        AUX1(:)=AUX1(:)+R(:)**2*POT(:)*Y0*PSI(:,IB)**2
        EHDEN(:)=EHDEN(:)+FOFI(IB)*AUX1(:)
!       == HARTREE-FOCK CORRECTION ===========================================
        AUX1(:)=0.D0
        CALL RADIALFOCK$VPSI(GID,NR,VFOCK,LOFI(IB),PSI(:,IB),AUX2)
        AUX2(:)=PSI(:,IB)*(0.5D0*AUX2-0.5D0*VFOCK%MUX(:)*Y0*PSI(:,IB))
        AUX1(:)=AUX1(:)+R(:)**2*AUX2(:)
        EXDEN(:)=EXDEN(:)+FOFI(IB)*AUX1(:)
!       == ADD TO SUM OVER STATES ============================================
        AUX(:)=AUX(:)+FOFI(IB)*AUX1(:)
      ENDDO
      AUX(:)=EKDEN+EHDEN+EXDEN
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,ETOT)
      CALL RADIAL$INTEGRATE(GID,NR,EKDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,EKIN)
      CALL RADIAL$INTEGRATE(GID,NR,EHDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,EH)
      CALL RADIAL$INTEGRATE(GID,NR,EXDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,EX)
PRINT*,'EKIN,EH,EX ',EKIN,EH,EX
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$PHASESHIFTSTATE(GID,NR,L,SO,DREL,G,POT,RBND,PHASE &
     &                                                              ,E,PHI)
!     **************************************************************************
!     **  FINDS A STATE WITH DETERMINED PHASESHIFT PHASE AT RADIUS RC         **
!     **  OF THE RADIAL SCHROEDINGER EQUATION AND                             **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G                **
!     **                                                                      **
!     **  THE ENERGY IS DETERMINED BY BISECTION ON THE                        **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.         **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.         **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: DREL(NR)! RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   ! INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: POT(NR) ! POTENTIAL
      REAL(8)    ,INTENT(IN)     :: RBND    ! RADIUS OF PHASE SHIFT
      REAL(8)    ,INTENT(IN)     :: PHASE   ! TARGET PHASE SHIFT
      REAL(8)    ,INTENT(INOUT)  :: E       ! ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) ! WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)    ,PARAMETER      :: TOL=1.D-15
      REAL(8)                    :: R(NR)
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRBND
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20 !MAX. FACTOR IN THE WAVEFUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
      REAL(8)                    :: PHIUP(NR)
      REAL(8)                    :: PHIDOWN(NR)
      LOGICAL(4)                 :: TUP,TDOWN
      REAL(8)                    :: ZUP,ZDOWN
      REAL(8)                    :: XUP,XDOWN
!     **************************************************************************
                                 CALL TRACE$PUSH('ATOMLIB$PHASESHIFTSTATE')
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRBND=IR
        IF(R(IR).GT.RBND) EXIT
      ENDDO
!          
!     ==========================================================================
!     ==========================================================================
!     ==  ITERATE TO BISECT ENERGY WITH A NODE AT RBOX                        ==
!     ==========================================================================
!     ==========================================================================
      ISTART=1
      X0=E
      DX=1.D-2
      CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      PHIUP(NR)=0.D0
      PHIDOWN(NR)=0.D0
      TUP=.FALSE.
      TDOWN=.FALSE.
      DO I=1,NITER
        E=X0
!
!       ========================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       ========================================================================
        IDIR=1
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,E,IDIR,PHI)
!       == CHECK FOR OVERFLOW
        IF(.NOT.(PHI(IRBND).GT.0.OR.PHI(IRBND).LE.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$STOP('ATOMLIB$PHASESHIFTSTATE')
        END IF
!
!       ========================================================================
!       == ESTIMATE PHASE SHIFT                                               ==
!       ========================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,0.D0,RBND,Z0)
        Z0=Z0-PHASE
        IF(Z0.GE.0.D0) THEN
          PHIUP(:)=PHI(:)
          TUP=.TRUE.
          ZUP=Z0
          XUP=X0
        ELSE
          PHIDOWN(:)=PHI(:)
          TDOWN=.TRUE.
          ZDOWN=Z0
          XDOWN=X0
        END IF
        IF(ABS(2.D0*DX).LE.TOL) EXIT
!       ========================================================================
!       ==  BISECTION                                                         ==
!       ========================================================================
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('ATOMLIB$PHASESHIFTSTATE')
      END IF
!
!     ==========================================================================
!     == MIX RESULT LINEARLY                                                  ==
!     ==========================================================================
      IF(.NOT.TUP.AND.TDOWN) THEN
        CALL ERROR$MSG('RARE EVENT CAPTURED')
        CALL ERROR$MSG('BISECTION DOES NOT BRACKET FROM BOTH SIDES')
        CALL ERROR$STOP('ATOMLIB$PHASESHIFTSTATE')
      END IF
      PHI(:)=(ZUP*PHIDOWN-ZDOWN*PHIUP)/(ZUP-ZDOWN)
!PRINT*,'FACTORS ',ZUP/(ZUP-ZDOWN),-ZDOWN/(ZUP-ZDOWN)
!PRINT*,'XUP  =',XUP  ,' ZUP   ',ZUP  ,' DIFFZ ',ZUP-ZDOWN
!PRINT*,'XDOWN=',XDOWN,' ZDOWN ',ZDOWN,' DIFFX ',XUP-XDOWN
      CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,0.D0,RBND,Z0)
      Z0=Z0-PHASE
!PRINT*,'L=',L,' E=',E,' ZUP-ZDOWN',ZUP-ZDOWN,' Z0=',Z0

                                 CALL TRACE$POP()
      RETURN
      END SUBROUTINE ATOMLIB$PHASESHIFTSTATE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB_WRITEPHI(FILE,GID,NR,NPHI,PHI)
!     **                                                                      **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: NPHI        
      REAL(8)     ,INTENT(IN) :: PHI(NR,NPHI)
      INTEGER(4)              :: IR,I
      REAL(8)                 :: R(NR)
      REAL(8)                 :: PHI1(NPHI)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      DO IR=1,NR
!        IF(R(IR).GT.3.D0.AND.MAXVAL(ABS(PHI(IR,:))).GT.1.D+3) CYCLE
        PHI1(:)=PHI(IR,:)
        DO I=1,NPHI  ! AVOID CONFLICT WITH XMGRACE
          IF(ABS(PHI1(I)).LT.1.D-60) PHI1(I)=0.D0
          PHI1(I)=MAX(-1.D+60,MIN(1.D+60,PHI1(I)))
        ENDDO
        WRITE(100,FMT='(F15.10,2X,20(E25.10,2X))')R(IR),PHI1
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!                                GARBADGE
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE ATOMLIB$AESCFWITHHF(GID,NR,TREL,NB,SOFI,LOFI,FOFI,EOFI,PSI &
!!$     &                              ,LRHOX,RBOX,DREL,POT,MUX,SCALE,VFOCK)
!!$!     **************************************************************************
!!$!     **************************************************************************
!!$      USE RADIALFOCK_MODULE
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN)    :: GID
!!$      INTEGER(4),INTENT(IN)    :: NR
!!$      INTEGER(4),INTENT(IN)    :: NB
!!$      LOGICAL(4),INTENT(IN)    :: TREL
!!$      INTEGER(4),INTENT(IN)    :: SOFI(NB)
!!$      INTEGER(4),INTENT(IN)    :: LOFI(NB)
!!$      REAL(8)   ,INTENT(IN)    :: FOFI(NB)
!!$      INTEGER(4),INTENT(IN)    :: LRHOX
!!$      REAL(8)   ,INTENT(IN)    :: RBOX
!!$      REAL(8)   ,INTENT(IN)    :: MUX(NR)
!!$      REAL(8)   ,INTENT(INOUT) :: SCALE
!!$      REAL(8)   ,INTENT(INOUT) :: EOFI(NB)
!!$      REAL(8)   ,INTENT(INOUT) :: DREL(NR)
!!$      REAL(8)   ,INTENT(INOUT) :: POT(NR)
!!$      REAL(8)   ,INTENT(INOUT) :: PSI(NR,NB)
!!$      TYPE(VFOCK_TYPE),INTENT(INOUT) :: VFOCK
!!$      REAL(8)   ,PARAMETER     :: TOL=1.D-6
!!$      INTEGER(4),PARAMETER     :: NITER=50
!!$      REAL(8)   ,PARAMETER     :: XMAX=1.D+15 !MAX. FACTOR IN THE WAVEFUNCTION
!!$      REAL(8)                  :: DPSI(NR,NB)
!!$      REAL(8)                  :: SPSI(NR,NB)
!!$      REAL(8)                  :: DRRPSI(NR,NB)
!!$      REAL(8)                  :: PI,Y0
!!$      REAL(8)                  :: SOFACTOR
!!$      REAL(8)                  :: R(NR)
!!$      REAL(8)                  :: AUX(NR),AUX1(NR),AUX2(NR)
!!$      REAL(8)                  :: G(NR)
!!$      REAL(8)                  :: RDPRIME(NR)
!!$      REAL(8)                  :: VAL,DER,VAL1,VAL2
!!$      REAL(8)                  :: ETOT
!!$      REAL(8)                  :: SVAR
!!$      INTEGER(4)               :: L,IB,IB1,IB2,I,IR,ISO
!!$      INTEGER(4)               :: IRBOX
!!$      INTEGER(4)               :: NN
!!$      INTEGER(4)               :: ITER
!!$      LOGICAL(4)               :: CONVG
!!$      REAL(8)                  :: E
!!$      REAL(8)                  :: POT1(NR)
!!$      INTEGER(4)               :: IRCL,IROUT
!!$      REAL(8)                  :: ROUT
!!$      LOGICAL(4)               :: TPR=.TRUE.
!!$      REAL(8)                  :: OVERLAP(NB,NB)
!!$REAL(8)  :: PSITEST(NR,NB)
!!$!     **************************************************************************
!!$      CALL ERROR$MSG('ROUTINE IS MARKED FOR DELETION')
!!$      CALL ERROR$STOP('ATOMLIB$AESCFWITHHF')
!!$      PI=4.D0*ATAN(1.D0)
!!$      Y0=1.D0/SQRT(4.D0*PI)
!!$      CALL RADIAL$R(GID,NR,R)
!!$!
!!$!     == INDEX OF THE FIRST GRIDPOINT OUTSIDE RBOX =============================
!!$      DO IR=1,NR
!!$        IRBOX=IR
!!$        IF(R(IR).GE.RBOX) EXIT
!!$      ENDDO
!!$
!!$POT(IRBOX+2:)=0.D0
!!$!
!!$      IF(TPR) THEN
!!$        PRINT*,'SCALE ',SCALE,'TREL ',TREL,' RBOX ',RBOX,R(IRBOX)
!!$        PRINT*,'EOFI ',0,EOFI
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     ==  RECREATE WAVE FUNCTIONS FOR TESTING                                 ==
!!$!     ==========================================================================
!!$      DO IB=1,NB
!!$        NN=0
!!$        DO I=1,IB-1
!!$          IF(LOFI(I).EQ.LOFI(IB)) NN=NN+1
!!$        ENDDO
!!$        IF(TREL) THEN
!!$          CALL SCHROEDINGER$DREL(GID,NR,POT,EOFI(IB),DREL)
!!$        ELSE 
!!$          DREL(:)=0.D0 
!!$        END IF
!!$        G(:)=0.D0
!!$        CALL ATOMLIB$BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),0.D0,RBOX &
!!$     &                         ,DREL,G,NN,POT,EOFI(IB),PSI(:,IB))
!!$        AUX(:)=R(:)**2*PSI(:,IB)**2
!!$        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$        CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$        PSI(:,IB)=PSI(:,IB)/SQRT(VAL)
!!$      ENDDO
!!$        
!!$      IF(TPR) THEN
!!$        CALL ATOMLIB_WRITEPHI('PSI_I.DAT',GID,NR,NB,PSI)
!!$        PRINT*,'EOFI ',-1,EOFI
!!$      END IF
!!$!
!!$!     ========================================================================
!!$!     ==  ESTIMATE FIRST-ORDER CHANGE OF THE ENERGY                         ==
!!$!     ========================================================================
!!$      IF(TPR) THEN
!!$        PSITEST=PSI
!!$        VFOCK%SCALE=0.D0
!!$        CALL ATOMLIB$ETOTWITHFOCK(GID,NR,RBOX,LRHOX,TREL,POT,VFOCK &
!!$    &                      ,NB,LOFI,SOFI,FOFI,PSITEST,VAL1)
!!$        VFOCK%SCALE=SCALE
!!$        PRINT*,'ENERGY WITHOUT HF                      :',VAL1
!!$        CALL ATOMLIB$ETOTWITHFOCK(GID,NR,RBOX,LRHOX,TREL,POT,VFOCK &
!!$    &                      ,NB,LOFI,SOFI,FOFI,PSITEST,VAL2)
!!$        PRINT*,'ENERGY WITH SCALED (HF-MU_X) CORRECTION:' ,VAL2
!!$        PRINT*,'SCALED (HF-MU_X) CORRECTION            :' ,VAL2-VAL1
!!$        VFOCK%SCALE=1.D0
!!$        CALL ATOMLIB$ETOTWITHFOCK(GID,NR,RBOX,LRHOX,TREL,POT,VFOCK &
!!$    &                      ,NB,LOFI,SOFI,FOFI,PSITEST,VAL2)
!!$        VFOCK%SCALE=SCALE
!!$        PRINT*,'100% (HF-MU_X) CORRECTION              :',VAL2-VAL1
!!$        AUX(:)=VFOCK%MUX
!!$        VFOCK%MUX=0.D0
!!$        CALL ATOMLIB$ETOTWITHFOCK(GID,NR,RBOX,LRHOX,TREL,POT,VFOCK &
!!$    &                      ,NB,LOFI,SOFI,FOFI,PSITEST,VAL2)
!!$        VFOCK%MUX=AUX
!!$        PRINT*,'100% HF CORRECTION                     :' ,VAL2-VAL1
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     ==  NOW START LOOP                                                      ==
!!$!     ==========================================================================
!!$      DO ITER=1,NITER
!!$!   
!!$!       ========================================================================
!!$!       ==  ENFORCE BOUNDARY CONDITIONS AT RBOX                               ==
!!$!       ========================================================================
!!$        PSI(IRBOX+2:,:)=0.D0
!!$!
!!$!       ========================================================================
!!$!       ==  VFOCK|PSI> AND PARTIAL_R*R|\PSI>                                  ==
!!$!       ========================================================================
!!$        AUX(:)=VFOCK%MUX
!!$        SCALE=VFOCK%SCALE
!!$        CALL RADIALFOCK$MAKEVFOCK(GID,NR,NB,LOFI,SOFI,FOFI,PSI &
!!$     &                               ,RBOX,LRHOX,VFOCK)
!!$        VFOCK%SCALE=SCALE
!!$        VFOCK%MUX=AUX
!!$!
!!$!       ========================================================================
!!$!       ==  CALCULATE TOTAL ENERGY (THIS A FAKE QUANTITY)                     ==
!!$!       ========================================================================
!!$        IF(TPR) THEN
!!$          DO IB=1,NB
!!$            AUX(:)=R(:)*PSI(:,IB)
!!$            VAL=AUX(IRBOX)
!!$            DER=(AUX(IRBOX+1)-AUX(IRBOX-1))/(R(IRBOX+1)-R(IRBOX-1))
!!$            AUX(IRBOX:)=VAL+(R(IRBOX:)-R(IRBOX))*DER
!!$            CALL RADIAL$VERLETD1(GID,NR,AUX,AUX1)
!!$            AUX1(1:2)=AUX1(3)  
!!$            DRRPSI(:,IB)=AUX1(:)
!!$          ENDDO
!!$          AUX(:)=0.D0
!!$          DO IB=1,NB
!!$            L=LOFI(IB)
!!$            IF(TREL) THEN
!!$              CALL SCHROEDINGER$DREL(GID,NR,POT,EOFI(IB),DREL)
!!$              CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
!!$            ELSE
!!$              DREL=0.D0
!!$              RDPRIME=0.D0
!!$            END IF
!!$            IF(SOFI(IB).EQ.0) THEN
!!$              SOFACTOR=0.D0
!!$            ELSE IF(SOFI(IB).EQ.1) THEN
!!$              SOFACTOR=REAL(L,KIND=8)
!!$            ELSE IF(SOFI(IB).EQ.-1) THEN
!!$              SOFACTOR=REAL(-L-1,KIND=8)
!!$            ELSE
!!$              CALL ERROR$MSG('ILLEGAL VALUE FOR SOFI ')
!!$              CALL ERROR$STOP('ATOMLIB$AESCFWITHHF')
!!$            END IF
!!$!           == KINETIC ENERGY ==================================================
!!$            AUX1(:)=0.5D0*(1.D0+DREL) &
!!$         &               *(DRRPSI(:,IB)**2+REAL(L*(L+1),KIND=8)*PSI(:,IB)**2)
!!$            AUX1(:)=AUX1(:)+0.5D0*(SOFACTOR+1.D0)*RDPRIME*R(:)*PSI(:,IB)**2
!!$!           == POTENTIAL ENERGY ================================================
!!$            AUX1(:)=AUX1(:)+R(:)**2*POT(:)*Y0*PSI(:,IB)**2
!!$!           == HARTREE-FOCK CORRECTION =========================================
!!$            CALL RADIALFOCK$VPSI(GID,NR,VFOCK,LOFI(IB),PSI(:,IB),AUX2)
!!$            AUX2(:)=PSI(:,IB)*(0.5D0*AUX2-0.5D0*VFOCK%MUX(:)*Y0*PSI(:,IB))
!!$            AUX1(:)=AUX1(:)+R(:)**2*AUX2(:)
!!$!           == ADD TO SUM OVER STATES ==========================================
!!$            AUX(:)=AUX(:)+FOFI(IB)*AUX1(:)
!!$          ENDDO
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,ETOT)
!!$        END IF
!!$!
!!$!       ========================================================================
!!$!       ==  CONSTRUCT |DPSI>=DE/DPSI                                          ==
!!$!       ========================================================================
!!$        DO IB=1,NB
!!$          L=LOFI(IB)
!!$          IF(TREL) THEN
!!$            CALL SCHROEDINGER$DREL(GID,NR,POT,EOFI(IB),DREL)
!!$            CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
!!$          ELSE
!!$            DREL=0.D0
!!$            RDPRIME=0.D0
!!$          END IF
!!$          IF(SOFI(IB).EQ.0) THEN
!!$            SOFACTOR=0.D0
!!$          ELSE IF(SOFI(IB).EQ.1) THEN
!!$            SOFACTOR=REAL(L,KIND=8)
!!$          ELSE IF(SOFI(IB).EQ.-1) THEN
!!$            SOFACTOR=REAL(-L-1,KIND=8)
!!$          ELSE
!!$            CALL ERROR$MSG('ILLEGAL VALUE FOR SOFI ')
!!$            CALL ERROR$STOP('ATOMLIB$AESCFWITHHF')
!!$          END IF
!!$!         == KINETIC ENERGY OF THE WAVE FUNCTION (NONRELATIVISTIC)
!!$PSITEST(:,IB)=R*PSI(:,IB)
!!$          CALL RADIAL$VERLETD2(GID,NR,R*PSI(:,IB),AUX2)
!!$          AUX(:)=R(:)*AUX2(:)-REAL(L*(L+1),KIND=8)*PSI(:,IB)
!!$          AUX(:)=-(1.D0+DREL)*AUX(:)
!!$          CALL RADIAL$VERLETD1(GID,NR,PSI(:,IB),AUX2)
!!$          AUX(:)=AUX(:)-RDPRIME*(R(:)**2*AUX2(:)-SOFACTOR*R(:)*PSI(:,IB))
!!$!         == POTENTIAL ENERGY ==================================================
!!$          AUX(:)=AUX(:)+2.D0*R(:)**2*POT(:)*Y0*PSI(:,IB)
!!$!         == HARTREE-FOCK CORRECTION ===========================================
!!$          CALL RADIALFOCK$VPSI(GID,NR,VFOCK,LOFI(IB),PSI(:,IB),AUX2)
!!$!PSITEST(:,IB)=AUX2(:)
!!$          AUX(:)=AUX(:)+2.D0*R(:)**2*AUX2(:)
!!$!         == CORRECT GLITCHES ==================================================
!!$          AUX(1:2)=0.D0
!!$!         == ESTIMATE ENERGY ===================================================
!!$          AUX1(:)=0.5D0*AUX(:)*PSI(:,IB)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
!!$          CALL RADIAL$VALUE(GID,NR,AUX2,RBOX,EOFI(IB))
!!$          AUX1(:)=R(:)**2*PSI(:,IB)**2
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
!!$          CALL RADIAL$VALUE(GID,NR,AUX2,RBOX,SVAR)
!!$          EOFI(IB)=EOFI(IB)/SVAR
!!$!         == REMOVE ENERGY TERM ================================================
!!$          AUX(:)=AUX(:)-2.D0*R(:)**2*EOFI(IB)*PSI(:,IB)
!!$!         == COPY INTO DPSI ====================================================
!!$          DPSI(2:,IB)=AUX(2:)/(2.D0*R(2:)**2)
!!$          DPSI(1:2,IB)=0.D0
!!$        ENDDO
!!$CALL ATOMLIB_WRITEPHI('TESTFOCK.DAT',GID,NR,VFOCK%NB,VFOCK%PSI)
!!$CALL ATOMLIB_WRITEPHI('TESTMUX.DAT',GID,NR,1,VFOCK%MUX)
!!$CALL ATOMLIB_WRITEPHI('TEST1A.DAT',GID,NR,NB,PSITEST)
!!$CALL ATOMLIB_WRITEPHI('TEST1.DAT',GID,NR,1,DPSI)
!!$!
!!$!       ========================================================================
!!$!       ==  DETERMINE PHIPRIME-PHIDOT*EPRIME                                  ==
!!$!       ========================================================================
!!$        DO IB=1,NB
!!$          L=LOFI(IB)
!!$          ISO=SOFI(IB)
!!$          E=EOFI(IB)
!!$          CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT,E,IRCL,IROUT)
!!$!         == BOUNDARY CONDITION PHI(ROUT)=0 ====================================
!!$          IF(R(IROUT).LT.RBOX) THEN
!!$            ROUT=R(IROUT)
!!$          ELSE
!!$            ROUT=RBOX
!!$            IROUT=IRBOX
!!$          END IF
!!$!         ==  SET KINETIC ENERGY TO ZERO BEYOND ROUT TO AVOID AN OVERFLOW ======
!!$!         == ATTENTION: THERE IS A STEP IN THE POTENTIAL =======================
!!$          POT1(:)=POT(:)
!!$          POT1(IROUT:)=E
!!$          IF(TREL) THEN
!!$            CALL SCHROEDINGER$DREL(GID,NR,POT1,E,DREL)
!!$          ELSE 
!!$            DREL(:)=0.D0 
!!$          END IF
!!$          G=DPSI(:,IB)
!!$          G(1:2)=0.D0
!!$          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,ISO,G,L,E,1,AUX1)
!!$          G=PSI(:,IB)
!!$          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,ISO,G,L,E,1,AUX2)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,VAL1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX2,ROUT,VAL2)
!!$          DPSI(:,IB)=AUX1(:)-AUX2(:)*VAL1/VAL2
!!$          DPSI(IROUT+3:,IB)=0.D0
!!$          IF(ROUT.LT.RBOX) DPSI(IROUT:,IB)=0.D0
!!$!
!!$!         == SMALL COMPONENT ===================================================
!!$          CALL RADIAL$VALUE(GID,NR,DPSI(:,IB),RBOX,VAL)
!!$          IF(TREL) THEN
!!$            G(:)=0.D0
!!$            CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,LOFI(IB),SOFI(IB) &
!!$     &                                       ,DREL,G,PSI(:,IB),SPSI(:,IB))
!!$          ELSE
!!$            SPSI(:,IB)=0.D0
!!$          END IF
!!$
!!$        ENDDO
!!$CALL ATOMLIB_WRITEPHI('TEST2.DAT',GID,NR,NB,DPSI)
!!$CALL ATOMLIB_WRITEPHI('TEST0.DAT',GID,NR,NB,PSI)
!!$CALL ATOMLIB_WRITEPHI('POT.DAT',GID,NR,1,POT)
!!$!   
!!$!       ========================================================================
!!$!       ==  CHECK CONVERGENCE                                                 ==
!!$!       ========================================================================
!!$        VAL=MAXVAL(ABS(DPSI(:IRBOX-3,:)))
!!$        IF(TPR) PRINT*,'ETOT=',ETOT,' DEV=',VAL
!!$        CONVG=(VAL.LT.TOL)
!!$        IF(CONVG) EXIT
!!$!   
!!$!       ========================================================================
!!$!       ==  PROPAGATE                                                         ==
!!$!       ========================================================================
!!$        PSI(:,:)=PSI(:,:)-DPSI(:,:)
!!$!   
!!$!       ========================================================================
!!$!       ==  ORTHOGONALIZE                                                     ==
!!$!       ========================================================================
!!$        DO IB1=1,NB
!!$          AUX(:)=R(:)**2*(PSI(:,IB1)**2+SPSI(:,IB1)**2)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          PSI(:,IB1)=PSI(:,IB1)/SQRT(VAL)
!!$          SPSI(:,IB1)=SPSI(:,IB1)/SQRT(VAL)
!!$        ENDDO
!!$
!!$!== ORTHOGONALIZATION KILLED CONVERGENCE
!!$        IF(TPR) THEN
!!$          OVERLAP(:,:)=0.D0
!!$          DO IB1=1,NB
!!$            DO IB2=1,IB1-1
!!$              IF(LOFI(IB1).NE.LOFI(IB2)) CYCLE
!!$              IF(SOFI(IB1).NE.SOFI(IB2)) CYCLE
!!$              AUX(:)=R(:)**2*(PSI(:,IB1)*PSI(:,IB2)+SPSI(:,IB1)*SPSI(:,IB2))
!!$              CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$              CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$              OVERLAP(IB1,IB2)=VAL
!!$!              PSI(:,IB1)=PSI(:,IB1)-PSI(:,IB2)*VAL
!!$!              SPSI(:,IB1)=SPSI(:,IB1)-SPSI(:,IB2)*VAL
!!$!IF(ABS(OVERLAP(IB1,IB2)).GT.1.D-7) THEN
!!$!WRITE(*,FMT='("OVERLAP",2I5,F20.10)')IB1,IB2,OVERLAP(IB1,IB2)
!!$!END IF
!!$            ENDDO
!!$            AUX(:)=R(:)**2*(PSI(:,IB1)**2+SPSI(:,IB1)**2)
!!$            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$            OVERLAP(IB1,IB1)=VAL
!!$! ORTHOGONALIZATION
!!$!            PSI(:,IB1)=PSI(:,IB1)/SQRT(VAL)
!!$!            SPSI(:,IB1)=SPSI(:,IB1)/SQRT(VAL)
!!$IF(ABS(OVERLAP(IB1,IB1)-1.D0).GT.1.D-7) THEN
!!$!WRITE(*,FMT='("OVERLAP",2I5,F20.10)')IB1,IB1,OVERLAP(IB1,IB1)
!!$END IF
!!$          ENDDO
!!$        END IF
!!$!
!!$      ENDDO    ! END OF ITERATION
!!$
!!$!
!!$      IF(.NOT.CONVG) THEN
!!$        CALL ERROR$MSG('LOOP FOR RADIAL HARTREE FOCK NOT CONVERGED')
!!$        CALL ERROR$I4VAL('NUMBER OF ITERATIONS',ITER)
!!$        VAL=MAXVAL(ABS(DPSI(:IRBOX-3,:)))
!!$        CALL ERROR$R8VAL('MAX DEVIATION',VAL)
!!$        CALL ERROR$R8VAL('UP TO RADIUS',R(IRBOX-3))
!!$        CALL ERROR$STOP('ATOMLIB$AESCFWITHHF')
!!$      END IF
!!$      IF(TPR) PRINT*,ITER,'ITERATIONSIN LOOP FOR RADIAL HARTREE FOCK'
!!$      IF(TPR) THEN
!!$        PRINT*,'EOFI ',ITER,EOFI
!!$        CALL ATOMLIB_WRITEPHI('PSI_F.DAT',GID,NR,NB,PSI)
!!$      END IF
!!$      RETURN
!!$      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE ATOMLIB$SPHERICALWITHHF(GID,NR,POT,DREL,VFOCK,SO,G,L,E,PSI)
!!$!     **************************************************************************
!!$!     **                                                                      **
!!$!     **************************************************************************
!!$      USE RADIALFOCK_MODULE, ONLY : VFOCK_TYPE
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN)    :: GID
!!$      INTEGER(4),INTENT(IN)    :: NR
!!$      TYPE(VFOCK_TYPE),INTENT(IN) :: VFOCK
!!$      REAL(8)   ,INTENT(IN)    :: DREL(NR)
!!$      REAL(8)   ,INTENT(IN)    :: POT(NR)
!!$      INTEGER(4),INTENT(IN)    :: L
!!$      INTEGER(4),INTENT(IN)    :: SO
!!$      REAL(8)   ,INTENT(IN)    :: G(NR)
!!$      REAL(8)   ,INTENT(IN)    :: E
!!$      REAL(8)   ,INTENT(OUT)   :: PSI(NR)
!!$      REAL(8)                  :: DPSI(NR)
!!$      REAL(8)                  :: VFOCKPSI(NR)
!!$      REAL(8)                  :: PI,Y0
!!$      REAL(8)                  :: SOFACTOR
!!$      REAL(8)                  :: R(NR)
!!$      REAL(8)                  :: AUX(NR),AUX1(NR),AUX2(NR)
!!$      REAL(8)                  :: G1(NR)
!!$      REAL(8)                  :: RDPRIME(NR)
!!$      REAL(8)                  :: VAL,NORM
!!$      INTEGER(4),PARAMETER     :: NITER=20
!!$      INTEGER(4)               :: ITER
!!$      REAL(8)   ,PARAMETER     :: TOL=1.D-5
!!$      LOGICAL(4)               :: CONVG
!!$      REAL(8)   ,PARAMETER     :: XMAX=1.D+15 !MAX. FACTOR IN THE WAVEFUNCTION
!!$      LOGICAL(4)               :: TPR=.TRUE.
!!$      LOGICAL(4)               :: THOM
!!$!     **************************************************************************
!!$      CALL ERROR$MSG('ROUTINE IS MARKED FOR DELETION')
!!$      CALL ERROR$STOP('ATOMLIB$SPHERICALWITHHF')
!!$      PI=4.D0*ATAN(1.D0)
!!$      Y0=1.D0/SQRT(4.D0*PI)
!!$      CALL RADIAL$R(GID,NR,R)
!!$      THOM=MAXVAL(ABS(G)).LT.1.D-8
!!$!
!!$      IF(TPR) THEN
!!$        PRINT*,'==============================================================='
!!$        PRINT*,'L=',L,' SO=',SO,' THOM ',THOM
!!$        PRINT*,'SCALE ',VFOCK%SCALE
!!$        PRINT*,'EOFI ',-1,E
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     ==  CHECK AND PROCESS SPIN ORBIT PARAMETER                              ==
!!$!     ==========================================================================
!!$      IF(SO.EQ.0) THEN
!!$        SOFACTOR=0.D0
!!$      ELSE IF(SO.EQ.1) THEN
!!$        SOFACTOR=REAL(L,KIND=8)
!!$      ELSE IF(SO.EQ.-1) THEN
!!$        SOFACTOR=REAL(-L-1,KIND=8)
!!$      ELSE
!!$        CALL ERROR$STOP('ILLEGAL VALUE FOR SO')
!!$        CALL ERROR$STOP('ATOMLIB$BOUNDSTATEWITHHF')
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     ==  SET RELATIVISTIC CORRECTION                                         ==
!!$!     ==========================================================================
!!$      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
!!$!
!!$!     ==========================================================================
!!$!     ==  RECREATE WAVE FUNCTION                                              ==
!!$!     ==========================================================================
!!$      CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,E,1,PSI)
!!$      IF(TPR) THEN
!!$        CALL ATOMLIB_WRITEPHI('PSI_I.DAT',GID,NR,1,PSI)
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     ==  NOW START LOOP                                                      ==
!!$!     ==========================================================================
!!$      DO ITER=1,NITER
!!$!
!!$!       ========================================================================
!!$!       ==  CONSTRUCT |DPSI>=DE/DPSI                                          ==
!!$!       ========================================================================
!!$!       == KINETIC ENERGY OF THE WAVE FUNCTION (NONRELATIVISTIC)
!!$        CALL RADIAL$VERLETD2(GID,NR,R*PSI,AUX2)
!!$        AUX(:)=R(:)*AUX2(:)-REAL(L*(L+1),KIND=8)*PSI(:)
!!$        AUX(:)=-(1.D0+DREL)*AUX(:)
!!$        CALL RADIAL$VERLETD1(GID,NR,PSI,AUX2)
!!$        AUX(:)=AUX(:)-RDPRIME*(R(:)**2*AUX2(:)-SOFACTOR*R(:)*PSI(:))
!!$!       == POTENTIAL ENERGY ==================================================
!!$        AUX(:)=AUX(:)+2.D0*R(:)**2*POT(:)*Y0*PSI(:)
!!$!       == HARTREE-FOCK CORRECTION ===========================================
!!$        CALL RADIALFOCK$VPSI(GID,NR,VFOCK,L,PSI,VFOCKPSI)
!!$        AUX(:)=AUX(:)+2.D0*R(:)**2*VFOCKPSI(:)
!!$!       == CORRECT GLITCHES ==================================================
!!$        AUX(1:2)=0.D0
!!$!       == REMOVE ENERGY TERM ================================================
!!$        AUX(:)=AUX(:)-2.D0*R(:)**2*E*PSI(:)
!!$!       == WEIGHT WITH OCCUPATION ============================================
!!$        DPSI(:)=AUX(:)/(2.D0*R(:)**2)-G(:)
!!$!
!!$!       ========================================================================
!!$!       ==  DETERMINE PHIPRIME-PHIDOT*EPRIME                                  ==
!!$!       ========================================================================
!!$        CALL ATOMLIB_WRITEPHI('DPSI1.DAT',GID,NR,1,DPSI)
!!$        G1=DPSI(:)
!!$        G1(1:2)=0.D0
!!$        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G1,L,E,1,DPSI)
!!$        CALL ATOMLIB_WRITEPHI('DPSI2.DAT',GID,NR,1,DPSI)
!!$!   
!!$!       ========================================================================
!!$!       ==  CHECK CONVERGENCE                                                 ==
!!$!       ========================================================================
!!$        AUX(:)=(R*PSI)**2
!!$        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$        CALL RADIAL$VALUE(GID,NR,AUX1,3.D0,NORM)
!!$        VAL=MAXVAL(ABS(DPSI(:)))/SQRT(NORM)
!!$        IF(TPR) PRINT*,' DEV=',VAL,NORM
!!$        CONVG=(VAL.LT.TOL)
!!$        IF(CONVG) EXIT
!!$!   
!!$!       ========================================================================
!!$!       ==  PROPAGATE                                                         ==
!!$!       ========================================================================
!!$        PSI(:)=PSI(:)-DPSI(:)
!!$!
!!$      ENDDO    ! END OF ITERATION
!!$!
!!$      IF(.NOT.CONVG) THEN
!!$        CALL ATOMLIB_WRITEPHI('PSI_F.DAT',GID,NR,1,PSI)
!!$        CALL ERROR$MSG('LOOP FOR RADIAL HARTREE FOCK NOT CONVERGED')
!!$        CALL ERROR$STOP('ATOMLIB$SPHERICALWITHHF')
!!$      END IF
!!$      IF(TPR) THEN
!!$        PRINT*,ITER,'ITERATIONSIN LOOP FOR RADIAL HARTREE FOCK'
!!$        CALL ATOMLIB_WRITEPHI('PSI_F.DAT',GID,NR,1,PSI)
!!$      END IF
!!$      RETURN
!!$      END SUBROUTINE ATOMLIB$SPHERICALWITHHF
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE ATOMLIB$UPDATEPAWSTATE(GID,NR,L,VPAW,G,RBND,E,PSI)
!!$!     **************************************************************************
!!$!     ** USES AN APPROXIMATE SOLUTION TO CONSTRUCT A SOLUTION FOR             **
!!$!     ** A PAW-TYPE HAMILTONIAN                                               **
!!$!     **                                                                      **
!!$!     ** THE BOUNDARY CONDITIONS (LOGARITHMIC DERIVATIVE AT RBND IS PRESERVED **
!!$!     ** IF RBND.LE.0 A FIXED ENERGY SOLUTION IS SEARCHED AND THE OUTER       **
!!$!     ** BOUNDARY CONDITION IS RELAXED.                                       **
!!$!     **                                                                      **
!!$!     **************************************************************************
!!$      USE RADIALPAW_MODULE, ONLY: VPAW_TYPE
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN)    :: GID
!!$      INTEGER(4),INTENT(IN)    :: NR
!!$      INTEGER(4),INTENT(IN)    :: L
!!$      TYPE(VPAW_TYPE),INTENT(IN) :: VPAW
!!$      REAL(8)   ,INTENT(IN)    :: G(NR)
!!$      REAL(8)   ,INTENT(IN)    :: RBND
!!$      REAL(8)   ,INTENT(INOUT) :: E
!!$      REAL(8)   ,INTENT(INOUT) :: PSI(NR)
!!$      REAL(8)                  :: PHASE
!!$      REAL(8)                  :: DREL(NR)
!!$      REAL(8)                  :: DPSI(NR)
!!$      REAL(8)                  :: OPSI(NR)
!!$      REAL(8)                  :: POT(NR)
!!$      REAL(8)                  :: PHIDOT(NR)
!!$      REAL(8)                  :: PHIPRIME(NR)
!!$      REAL(8)                  :: VALDOT,DERDOT,VALPRIME,DERPRIME,VAL0,DER0
!!$      REAL(8)                  :: PI,Y0
!!$      REAL(8)                  :: R(NR)
!!$      REAL(8)                  :: AUX(NR),AUX1(NR),AUX2(NR)
!!$      REAL(8)                  :: G1(NR)
!!$      REAL(8)                  :: VAL,DER,SVAR
!!$      INTEGER(4)               :: IR
!!$      INTEGER(4),PARAMETER     :: NITER=20
!!$      INTEGER(4)               :: ITER
!!$      REAL(8)   ,PARAMETER     :: TOL=1.D-8
!!$      LOGICAL(4)               :: CONVG
!!$      REAL(8)   ,PARAMETER     :: XMAX=1.D+15 !MAX. FACTOR IN THE WAVEFUNCTION
!!$      INTEGER(4)               :: IRBND
!!$      REAL(8)                  :: NORM,E1
!!$      LOGICAL(4)               :: TPR=.TRUE.
!!$      LOGICAL(4)               :: TFIXE
!!$      INTEGER(4)               :: NN
!!$!     **************************************************************************
!!$      CALL ERROR$MSG('ROUTINE IS MARKED FOR DELETION')
!!$      CALL ERROR$STOP('ATOMLIB$SPHERICALWITHHF')
!!$      PI=4.D0*ATAN(1.D0)
!!$      Y0=1.D0/SQRT(4.D0*PI)
!!$      CALL RADIAL$R(GID,NR,R)
!!$      TFIXE=RBND.LE.0.D0
!!$      DREL(:)=0.D0
!!$!
!!$!     :: REMARK: SCHROEDINGER$PHASESHIFT TAKES VALUE AND DERIVATIVE FROM A    ::
!!$!     :: TWO-POINT FORMULA, AND THEREFORE IS NOT FULLY CONSISTENT WITH        ::
!!$!     :: RADIAL$VALUE AND RADIAL$DERIVATIVE                                   ::
!!$      CALL SCHROEDINGER$PHASESHIFT(GID,NR,PSI,0.D0,RBND,PHASE)
!!$!
!!$      DO IR=1,NR
!!$        IRBND=IR
!!$        IF(TFIXE.AND.R(IR).GT.3.D0) EXIT
!!$        IF(.NOT.TFIXE.AND.R(IR).GT.RBND) EXIT
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     ==  CONSTRUCT AN LOCAL, APPROXIMATE POTENTIAL
!!$!     ==========================================================================
!!$      AUX(:)=1.D0
!!$      CALL RADIALPAW$VPSI(GID,NR,VPAW,L,AUX,AUX1,OPSI)
!!$      POT(:)=AUX1-E*OPSI(:)+E*AUX(:)
!!$POT=VPAW%PSPOT
!!$!
!!$!     ==========================================================================
!!$!     ==  NOW START LOOP                                                      ==
!!$!     ==========================================================================
!!$      DO ITER=1,NITER
!!$!
!!$!       ========================================================================
!!$!       ==  CONSTRUCT |DPSI>=DE/DPSI                                          ==
!!$!       ========================================================================
!!$!       == KINETIC ENERGY OF THE WAVE FUNCTION (NONRELATIVISTIC)
!!$        CALL RADIAL$VERLETD2(GID,NR,R*PSI,AUX1)
!!$        AUX(:)=-R(:)*AUX(:)+REAL(L*(L+1),KIND=8)*PSI(:)
!!$!       == REMOVE FACTOR 2*R**2 ================================================
!!$        AUX(:)=AUX(:)/(2.D0*R(:)**2)
!!$!       == CORRECT GLITCHES ====================================================
!!$        AUX(1:2)=0.D0
!!$!       == POTENTIAL ENERGY ====================================================
!!$        CALL RADIALPAW$VPSI(GID,NR,VPAW,L,PSI,AUX1,OPSI)
!!$        AUX(:)=AUX(:)+AUX1(:)
!!$!       == ESTIMATE ENERGY =====================================================
!!$        IF(.NOT.TFIXE) THEN
!!$          AUX1(:)=R(:)**2*PSI(:)*(AUX(:)-G(:))
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
!!$          CALL RADIAL$VALUE(GID,NR,AUX2,RBND,E)
!!$          AUX1(:)=R(:)**2*PSI(:)*OPSI(:)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
!!$          CALL RADIAL$VALUE(GID,NR,AUX2,RBND,NORM)
!!$          E=E/NORM
!!$!         == REMOVE ENERGY TERM ================================================
!!$          AUX(:)=AUX(:)-E*OPSI(:)
!!$        ELSE
!!$          AUX1(:)=R(:)**2*PSI(:)*OPSI(:)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
!!$          CALL RADIAL$VALUE(GID,NR,AUX2,3.D0,NORM)
!!$        END IF
!!$!       == SUBTRACT INHOMOGENEITY ==============================================
!!$        AUX(:)=AUX(:)-G(:)
!!$!       == MAP INTO DPSI
!!$        DPSI(:)=-AUX(:)
!!$!
!!$!       ========================================================================
!!$!       ==  DETERMINE PHIPRIME-PHIDOT*EPRIME                                  ==
!!$!       ========================================================================
!!$        G1=DPSI(:)
!!$        G1(1:2)=0.D0
!!$        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,0,G1,L,E,1,PHIPRIME)
!!$        IF(.NOT.TFIXE) THEN
!!$          G1=PSI(:)
!!$          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,0,G1,L,E,1,PHIDOT)
!!$          CALL RADIAL$VALUE(GID,NR,PSI,RBND,VAL0)
!!$          CALL RADIAL$DERIVATIVE(GID,NR,PSI,RBND,DER0)
!!$          CALL RADIAL$VALUE(GID,NR,PHIPRIME,RBND,VALPRIME)
!!$          CALL RADIAL$DERIVATIVE(GID,NR,PHIPRIME,RBND,DERPRIME)
!!$          CALL RADIAL$VALUE(GID,NR,PHIDOT,RBND,VALDOT)
!!$          CALL RADIAL$DERIVATIVE(GID,NR,PHIDOT,RBND,DERDOT)
!!$          CALL SCHROEDINGER$RESOLVEPHASESHIFT(PHASE,VAL,DER,NN)
!!$          VAL0=VAL*DER0-DER*VAL0
!!$          VALPRIME=VAL*DERPRIME-DER*VALPRIME
!!$          VALDOT=VAL*DERDOT-DER*VALDOT
!!$          SVAR=(VAL0+VALPRIME)/VALDOT  !DELTA EPSILON
!!$          DPSI(:)=PHIPRIME(:)-PHIDOT(:)*SVAR
!!$        ELSE
!!$          DPSI(:)=PHIPRIME(:)
!!$        END IF
!!$!   
!!$!       ========================================================================
!!$!       ==  CHECK CONVERGENCE                                                 ==
!!$!       ========================================================================
!!$        VAL=MAXVAL(ABS(DPSI(:IRBND-3)))/SQRT(NORM)
!!$        IF(TPR) PRINT*,' DEV=',VAL
!!$        CONVG=(VAL.LT.TOL)
!!$        IF(CONVG) EXIT
!!$!   
!!$!       ========================================================================
!!$!       ==  PROPAGATE                                                         ==
!!$!       ========================================================================
!!$        PSI(:)=PSI(:)+DPSI(:)
!!$!
!!$      ENDDO    ! END OF ITERATION
!!$!
!!$      IF(.NOT.CONVG) THEN
!!$        CALL ERROR$MSG('LOOP FOR RADIAL PAW NOT CONVERGED')
!!$        CALL ERROR$STOP('ATOMLIB$UPDATEPAWSTATE')
!!$      END IF
!!$      RETURN
!!$      END SUBROUTINE ATOMLIB$UPDATEPAWSTATE
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE ATOMLIB$PHASESHIFTSTATEWITHHF(GID,NR,L,SO,DREL,G,POT &
!!$     &                             ,VFOCK,RBND,PHASE,E,PSI)
!!$!     **************************************************************************
!!$!     **                                                                      **
!!$!     **************************************************************************
!!$      USE RADIALFOCK_MODULE, ONLY : VFOCK_TYPE
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN)    :: GID
!!$      INTEGER(4),INTENT(IN)    :: NR
!!$      INTEGER(4),INTENT(IN)    :: L
!!$      INTEGER(4),INTENT(IN)    :: SO
!!$      REAL(8)   ,INTENT(IN)    :: DREL(NR)
!!$      REAL(8)   ,INTENT(IN)    :: POT(NR)
!!$      TYPE(VFOCK_TYPE),INTENT(IN) :: VFOCK
!!$      REAL(8)   ,INTENT(IN)    :: G(NR)
!!$      REAL(8)   ,INTENT(IN)    :: RBND
!!$      REAL(8)   ,INTENT(IN)    :: PHASE
!!$      REAL(8)   ,INTENT(OUT)   :: E
!!$      REAL(8)   ,INTENT(OUT)   :: PSI(NR)
!!$      REAL(8)                  :: DPSI(NR)
!!$      REAL(8)                  :: PHIDOT(NR)
!!$      REAL(8)                  :: PHIPRIME(NR)
!!$      REAL(8)                  :: VALDOT,DERDOT,VALPRIME,DERPRIME,VAL0,DER0
!!$      REAL(8)                  :: PI,Y0
!!$      REAL(8)                  :: SOFACTOR
!!$      REAL(8)                  :: R(NR)
!!$      REAL(8)                  :: AUX(NR),AUX1(NR),AUX2(NR)
!!$      REAL(8)                  :: G1(NR)
!!$      REAL(8)                  :: RDPRIME(NR)
!!$      REAL(8)                  :: VAL,DER,SVAR
!!$      INTEGER(4)               :: IR
!!$      INTEGER(4),PARAMETER     :: NITER=20
!!$      INTEGER(4)               :: ITER
!!$      REAL(8)   ,PARAMETER     :: TOL=1.D-8
!!$      LOGICAL(4)               :: CONVG
!!$      REAL(8)   ,PARAMETER     :: XMAX=1.D+15 !MAX. FACTOR IN THE WAVEFUNCTION
!!$      INTEGER(4)               :: IRBND
!!$      REAL(8)                  :: NORM
!!$      LOGICAL(4)               :: TPR=.TRUE.
!!$      LOGICAL(4)               :: THOM
!!$      INTEGER(4)               :: NN
!!$!     **************************************************************************
!!$      CALL ERROR$MSG('ROUTINE IS MARKED FOR DELETION')
!!$      CALL ERROR$STOP('ATOMLIB$PHASESHIFTSTATEWITHHF')
!!$      PI=4.D0*ATAN(1.D0)
!!$      Y0=1.D0/SQRT(4.D0*PI)
!!$      CALL RADIAL$R(GID,NR,R)
!!$      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
!!$      IF(RBND.LE.0) THEN
!!$        CALL ERROR$MSG('CONSTANT ENERGY SOLUTION NOT IMPLEMENTED')
!!$        CALL ERROR$MSG('RBND MUST BE POSITIVE')
!!$        CALL ERROR$STOP('ATOMLIB$PHASESHIFTSTATEWITHHF')
!!$      END IF
!!$!      
!!$      THOM=MAXVAL(ABS(G)).LT.1.D-8
!!$      DO IR=1,NR
!!$        IRBND=IR
!!$        IF(R(IR).GT.RBND) EXIT
!!$      ENDDO
!!$      IF(SO.EQ.0) THEN
!!$        SOFACTOR=0.D0
!!$      ELSE IF(SO.EQ.1) THEN
!!$        SOFACTOR=REAL(L,KIND=8)
!!$      ELSE IF(SO.EQ.-1) THEN
!!$        SOFACTOR=REAL(-L-1,KIND=8)
!!$      ELSE
!!$        CALL ERROR$STOP('ILLEGAL VALUE FOR SO')
!!$        CALL ERROR$STOP('ATOMLIB$BOUNDSTATEWITHHF')
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     ==  RECREATE WAVE FUNCTION                                              ==
!!$!     ==========================================================================
!!$      CALL ATOMLIB$PHASESHIFTSTATE(GID,NR,L,SO,DREL,G,POT,RBND,PHASE,E,PSI)
!!$!
!!$      IF(THOM) THEN
!!$        AUX(:)=R(:)**2*PSI(:)**2
!!$        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$        CALL RADIAL$VALUE(GID,NR,AUX1,RBND,VAL)
!!$        PSI(:)=PSI(:)/SQRT(VAL)
!!$PRINT*,'WAVE FUNCTION NORMALIZED; THOM=',THOM
!!$      END IF
!!$        
!!$      IF(TPR) THEN
!!$        CALL ATOMLIB_WRITEPHI('PSI_I.DAT',GID,NR,1,PSI)
!!$        PRINT*,'EOFI ',0,E
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     ==  NOW START LOOP                                                      ==
!!$!     ==========================================================================
!!$      DO ITER=1,NITER
!!$!
!!$!       ========================================================================
!!$!       ==  CONSTRUCT |DPSI>=DE/DPSI                                          ==
!!$!       ========================================================================
!!$!       == KINETIC ENERGY OF THE WAVE FUNCTION (NONRELATIVISTIC)
!!$        CALL RADIAL$VERLETD2(GID,NR,R*PSI,AUX1)
!!$        AUX(:)=R(:)*AUX1(:)
!!$        AUX(:)=AUX(:)-REAL(L*(L+1),KIND=8)*PSI(:)
!!$        AUX(:)=-(1.D0+DREL)*AUX(:)
!!$        CALL RADIAL$VERLETD1(GID,NR,PSI,AUX1)
!!$        AUX(:)=AUX(:)-RDPRIME(:)*(R(:)**2*AUX1(:)-SOFACTOR*R(:)*PSI(:))
!!$!       == POTENTIAL ENERGY ====================================================
!!$        AUX(:)=AUX(:)+2.D0*R(:)**2*POT(:)*Y0*PSI(:)
!!$!       == REMOVE FACTOR 2*R**2 ================================================
!!$        AUX(:)=AUX(:)/(2.D0*R(:)**2)
!!$!       == HARTREE-FOCK CORRECTION =============================================
!!$        CALL RADIALFOCK$VPSI(GID,NR,VFOCK,L,PSI,AUX1)
!!$        AUX(:)=AUX(:)+AUX1(:)
!!$!       == CORRECT GLITCHES ====================================================
!!$        AUX(1:2)=0.D0
!!$!       == ESTIMATE ENERGY =====================================================
!!$        AUX1(:)=R(:)**2*PSI(:)*(AUX(:)-G(:))
!!$        CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
!!$        CALL RADIAL$VALUE(GID,NR,AUX2,RBND,E)
!!$        AUX1(:)=R(:)**2*PSI(:)**2
!!$        CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
!!$        CALL RADIAL$VALUE(GID,NR,AUX2,RBND,NORM)
!!$        E=E/NORM
!!$PRINT*,'EE',E,NORM
!!$!       == REMOVE ENERGY TERM ================================================
!!$        AUX(:)=AUX(:)-E*PSI(:)
!!$!       == SUBTRACT INHOMOGENEITY ==============================================
!!$        AUX(:)=AUX(:)-G(:)
!!$!       == MAP INTO DPSI
!!$        DPSI(:)=-AUX(:)
!!$!
!!$!       ========================================================================
!!$!       ==  DETERMINE PHIPRIME-PHIDOT*EPRIME                                  ==
!!$!       ========================================================================
!!$        G1=DPSI(:)
!!$        G1(1:2)=0.D0
!!$        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G1,L,E,1,PHIPRIME)
!!$        G1=PSI(:)
!!$        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G1,L,E,1,PHIDOT)
!!$!
!!$!       :: REMARK: SCHROEDINGER$PHASESHIFT TAKES VALUE AND DERIVATIVE FROM A  ::
!!$!       :: TWO-POINT FORMULA, AND THEREFORE IS NOT FULLY CONSISTENT WITH      ::
!!$!       :: RADIAL$VALUE AND RADIAL$DERIVATIVE                                 ::
!!$        CALL RADIAL$VALUE(GID,NR,PSI,RBND,VAL0)
!!$        CALL RADIAL$DERIVATIVE(GID,NR,PSI,RBND,DER0)
!!$        CALL RADIAL$VALUE(GID,NR,PHIPRIME,RBND,VALPRIME)
!!$        CALL RADIAL$DERIVATIVE(GID,NR,PHIPRIME,RBND,DERPRIME)
!!$        CALL RADIAL$VALUE(GID,NR,PHIDOT,RBND,VALDOT)
!!$        CALL RADIAL$DERIVATIVE(GID,NR,PHIDOT,RBND,DERDOT)
!!$        CALL SCHROEDINGER$RESOLVEPHASESHIFT(PHASE,VAL,DER,NN)
!!$        VAL0=VAL*DER0-DER*VAL0
!!$        VALPRIME=VAL*DERPRIME-DER*VALPRIME
!!$        VALDOT=VAL*DERDOT-DER*VALDOT
!!$        SVAR=(VAL0+VALPRIME)/VALDOT  !DELTA EPSILON
!!$        DPSI(:)=PHIPRIME(:)-PHIDOT(:)*SVAR
!!$!   
!!$!       ========================================================================
!!$!       ==  CHECK CONVERGENCE                                                 ==
!!$!       ========================================================================
!!$        VAL=MAXVAL(ABS(DPSI(:IRBND-3)))/SQRT(NORM)
!!$        IF(TPR) PRINT*,' DEV=',VAL
!!$        CONVG=(VAL.LT.TOL)
!!$        IF(CONVG) EXIT
!!$!   
!!$!       ========================================================================
!!$!       ==  PROPAGATE                                                         ==
!!$!       ========================================================================
!!$        PSI(:)=PSI(:)+DPSI(:)
!!$! CALL SCHROEDINGER$PHASESHIFT(GID,NR,PSI,0.D0,RBND,SVAR)
!!$!PRINT*,'TARGET PHASE+ ',PHASE,' ACTUAL PHASE=',SVAR,' R=',RBND
!!$!   
!!$!       ========================================================================
!!$!       ==  NORMALIZE                                                         ==
!!$!       ========================================================================
!!$        IF(THOM) THEN
!!$          AUX(:)=R(:)**2*PSI(:)**2
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBND,VAL)
!!$          PSI(:)=PSI(:)/SQRT(VAL)
!!$        END IF
!!$        IF(TPR) THEN
!!$          PRINT*,'EOFI ',ITER,E
!!$        END IF
!!$!
!!$      ENDDO    ! END OF ITERATION
!!$!
!!$      IF(.NOT.CONVG) THEN
!!$        CALL ERROR$MSG('LOOP FOR RADIAL HARTREE FOCK NOT CONVERGED')
!!$        CALL ERROR$STOP('ATOMLIB$PHASESHIFTSTATEWITHHF')
!!$      END IF
!!$      IF(TPR) PRINT*,ITER,'ITERATIONSIN LOOP FOR RADIAL HARTREE FOCK'
!!$      IF(THOM) THEN
!!$        AUX(:)=R(:)**2*PSI(:)**2
!!$        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$        CALL RADIAL$VALUE(GID,NR,AUX1,RBND,VAL)
!!$        PSI(:)=PSI(:)/SQRT(VAL)
!!$      END IF
!!$      IF(TPR) THEN
!!$        PRINT*,'EOFI ',ITER,E
!!$        CALL ATOMLIB_WRITEPHI('PSI_F.DAT',GID,NR,1,PSI)
!!$      END IF
!!$!
!!$      RETURN
!!$      END SUBROUTINE ATOMLIB$PHASESHIFTSTATEWITHHF
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE ATOMLIB$BOUNDSTATEWITHHF(GID,NR,TREL,RBOX,POT &
!!$     &                             ,VFOCK,L,SO,NN,G,E,PSI,TPSI)
!!$!     **************************************************************************
!!$!     **                                                                      **
!!$!     **************************************************************************
!!$      USE RADIALFOCK_MODULE, ONLY : VFOCK_TYPE
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN)    :: GID
!!$      INTEGER(4),INTENT(IN)    :: NR
!!$      LOGICAL(4),INTENT(IN)    :: TREL
!!$      REAL(8)   ,INTENT(IN)    :: RBOX
!!$      TYPE(VFOCK_TYPE),INTENT(IN) :: VFOCK
!!$      REAL(8)   ,INTENT(IN)    :: POT(NR)
!!$      INTEGER(4),INTENT(IN)    :: NN
!!$      INTEGER(4),INTENT(IN)    :: L
!!$      INTEGER(4),INTENT(IN)    :: SO
!!$      REAL(8)   ,INTENT(IN)    :: G(NR)
!!$      REAL(8)   ,INTENT(OUT)   :: PSI(NR)
!!$      REAL(8)   ,INTENT(OUT)   :: E
!!$      REAL(8)   ,INTENT(OUT)   :: TPSI(NR)
!!$      REAL(8)   ,PARAMETER     :: TOL=1.D-6
!!$      INTEGER(4),PARAMETER     :: NITER=20
!!$      REAL(8)   ,PARAMETER     :: XMAX=1.D+15 !MAX. FACTOR IN THE WAVEFUNCTION
!!$      LOGICAL(4),PARAMETER     :: TPR=.TRUE.
!!$      REAL(8)                  :: DREL(NR)
!!$      REAL(8)                  :: DPSI(NR)
!!$      REAL(8)                  :: VFOCKPSI(NR)
!!$      REAL(8)                  :: PI,Y0
!!$      REAL(8)                  :: SOFACTOR
!!$      REAL(8)                  :: R(NR)
!!$      REAL(8)                  :: AUX(NR),AUX1(NR),AUX2(NR)
!!$      REAL(8)                  :: G1(NR)
!!$      REAL(8)                  :: RDPRIME(NR)
!!$      REAL(8)                  :: VAL,VAL1,VAL2
!!$      INTEGER(4)               :: IR,ISO
!!$      INTEGER(4)               :: IRBOX
!!$      INTEGER(4)               :: ITER
!!$      LOGICAL(4)               :: CONVG
!!$      REAL(8)                  :: POT1(NR)
!!$      INTEGER(4)               :: IRCL,IROUT
!!$      REAL(8)                  :: ROUT
!!$      REAL(8)                  :: NORM
!!$      LOGICAL(4)               :: THOM
!!$!     **************************************************************************
!!$      CALL ERROR$MSG('ROUTINE IS MARKED FOR DELETION')
!!$      CALL ERROR$STOP('ATOMLIB$BOUNDSTATEWITHHF')
!!$      PI=4.D0*ATAN(1.D0)
!!$      Y0=1.D0/SQRT(4.D0*PI)
!!$      CALL RADIAL$R(GID,NR,R)
!!$      THOM=MAXVAL(ABS(G)).LT.1.D-8
!!$!
!!$!     == INDEX OF THE FIRST GRIDPOINT OUTSIDE RBOX =============================
!!$      DO IR=1,NR
!!$        IRBOX=IR
!!$        IF(R(IR).GE.RBOX) EXIT
!!$      ENDDO
!!$!
!!$      IF(TPR) THEN
!!$        PRINT*,'==============================================================='
!!$        PRINT*,'L=',L,' SO=',SO,' NN=',NN,' THOM ',THOM
!!$        PRINT*,'SCALE ',VFOCK%SCALE,'TREL ',TREL,' RBOX ',RBOX,R(IRBOX)
!!$        PRINT*,'EOFI ',-1,E
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     ==  RECREATE WAVE FUNCTION                                              ==
!!$!     ==========================================================================
!!$      E=0.D0 
!!$      IF(TREL) THEN
!!$        CALL SCHROEDINGER$DREL(GID,NR,POT,E,DREL)
!!$      ELSE 
!!$        DREL(:)=0.D0 
!!$      END IF
!!$      CALL ATOMLIB$BOUNDSTATE(GID,NR,L,SO,0.D0,RBOX,DREL,G,NN,POT,E,PSI)
!!$!
!!$      IF(TPR.AND.THOM) THEN
!!$        AUX(:)=R(:)**2*PSI(:)**2
!!$        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$        CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$        PSI(:)=PSI(:)/SQRT(VAL)
!!$      END IF
!!$        
!!$      IF(TPR) THEN
!!$        CALL ATOMLIB_WRITEPHI('PSI_I.DAT',GID,NR,1,PSI)
!!$        PRINT*,'EOFI ',0,E
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     ==  NOW START LOOP                                                      ==
!!$!     ==========================================================================
!!$      DO ITER=1,NITER
!!$!
!!$!       ========================================================================
!!$!       ==  CONSTRUCT |DPSI>=DE/DPSI                                          ==
!!$!       ========================================================================
!!$        IF(TREL) THEN
!!$          CALL SCHROEDINGER$DREL(GID,NR,POT,E,DREL)
!!$          CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
!!$        ELSE
!!$          DREL=0.D0
!!$          RDPRIME=0.D0
!!$        END IF
!!$        IF(SO.EQ.0) THEN
!!$          SOFACTOR=0.D0
!!$        ELSE IF(SO.EQ.1) THEN
!!$          SOFACTOR=REAL(L,KIND=8)
!!$        ELSE IF(SO.EQ.-1) THEN
!!$          SOFACTOR=REAL(-L-1,KIND=8)
!!$        ELSE
!!$          CALL ERROR$STOP('ILLEGAL VALUE FOR SO')
!!$          CALL ERROR$STOP('ATOMLIB$BOUNDSTATEWITHHF')
!!$        END IF
!!$!       == KINETIC ENERGY OF THE WAVE FUNCTION (NONRELATIVISTIC)
!!$        CALL RADIAL$VERLETD2(GID,NR,R*PSI,AUX2)
!!$        AUX(:)=R(:)*AUX2(:)-REAL(L*(L+1),KIND=8)*PSI(:)
!!$        AUX(:)=-(1.D0+DREL)*AUX(:)
!!$        CALL RADIAL$VERLETD1(GID,NR,PSI,AUX2)
!!$        AUX(:)=AUX(:)-RDPRIME*(R(:)**2*AUX2(:)-SOFACTOR*R(:)*PSI(:))
!!$!       == POTENTIAL ENERGY ==================================================
!!$        AUX(:)=AUX(:)+2.D0*R(:)**2*POT(:)*Y0*PSI(:)
!!$!       == HARTREE-FOCK CORRECTION ===========================================
!!$        CALL RADIALFOCK$VPSI(GID,NR,VFOCK,L,PSI,VFOCKPSI)
!!$        AUX(:)=AUX(:)+2.D0*R(:)**2*VFOCKPSI(:)
!!$!       == CORRECT GLITCHES ==================================================
!!$        AUX(1:2)=0.D0
!!$!       == ESTIMATE ENERGY ===================================================
!!$        AUX1(:)=0.5D0*AUX(:)*PSI(:)-R(:)**2*PSI(:)*G(:)
!!$        CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
!!$        CALL RADIAL$VALUE(GID,NR,AUX2,RBOX,E)
!!$        AUX1(:)=R(:)**2*PSI(:)**2
!!$        CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
!!$        CALL RADIAL$VALUE(GID,NR,AUX2,RBOX,NORM)
!!$        E=E/NORM
!!$!       == REMOVE ENERGY TERM ================================================
!!$        AUX(:)=AUX(:)-2.D0*R(:)**2*E*PSI(:)
!!$!       == WEIGHT WITH OCCUPATION ============================================
!!$        DPSI(:)=AUX(:)/(2.D0*R(:)**2)-G(:)
!!$!
!!$!       ========================================================================
!!$!       ==  DETERMINE PHIPRIME-PHIDOT*EPRIME                                  ==
!!$!       ========================================================================
!!$        ISO=SO
!!$        CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT,E,IRCL,IROUT)
!!$!       == BOUNDARY CONDITION PHI(ROUT)=0 ====================================
!!$        IF(R(IROUT).LT.RBOX) THEN
!!$          ROUT=R(IROUT)
!!$        ELSE
!!$          ROUT=RBOX
!!$          IROUT=IRBOX
!!$        END IF
!!$!       ==  SET KINETIC ENERGY TO ZERO BEYOND ROUT TO AVOID AN OVERFLOW ======
!!$!       == ATTENTION: THERE IS A STEP IN THE POTENTIAL =======================
!!$        POT1(:)=POT(:)
!!$        POT1(IROUT:)=E
!!$        IF(TREL) THEN
!!$          CALL SCHROEDINGER$DREL(GID,NR,POT1,E,DREL)
!!$        ELSE 
!!$          DREL(:)=0.D0 
!!$        END IF
!!$        G1=DPSI(:)
!!$        G1(1:2)=0.D0
!!$        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,ISO,G1,L,E,1,AUX1)
!!$        G1=PSI(:)
!!$        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,ISO,G1,L,E,1,AUX2)
!!$        CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,VAL1)
!!$        CALL RADIAL$VALUE(GID,NR,AUX2,ROUT,VAL2)
!!$        DPSI(:)=AUX1(:)-AUX2(:)*VAL1/VAL2
!!$        DPSI(IROUT+3:)=0.D0
!!$        IF(ROUT.LT.RBOX) DPSI(IROUT:)=0.D0
!!$        CALL RADIAL$VALUE(GID,NR,DPSI(:),RBOX,VAL)
!!$!   
!!$!       ========================================================================
!!$!       ==  CHECK CONVERGENCE                                                 ==
!!$!       ========================================================================
!!$        VAL=MAXVAL(ABS(DPSI(:IRBOX-3)))/SQRT(NORM)
!!$        IF(TPR) PRINT*,' DEV=',VAL,' ITER=',ITER
!!$        CONVG=(VAL.LT.TOL)
!!$        IF(CONVG) EXIT
!!$!   
!!$!       ========================================================================
!!$!       ==  PROPAGATE                                                         ==
!!$!       ========================================================================
!!$        PSI(:)=PSI(:)-DPSI(:)
!!$!   
!!$!       ========================================================================
!!$!       ==  NORMALIZE                                                         ==
!!$!       ========================================================================
!!$        IF(TPR.AND.THOM) THEN
!!$          AUX(:)=R(:)**2*PSI(:)**2
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          PSI(:)=PSI(:)/SQRT(VAL)
!!$        END IF
!!$        IF(TPR) THEN
!!$          PRINT*,'EOFI ',ITER,E
!!$        END IF
!!$!
!!$      ENDDO    ! END OF ITERATION
!!$!
!!$      IF(.NOT.CONVG) THEN
!!$        CALL ERROR$MSG('LOOP FOR RADIAL HARTREE FOCK NOT CONVERGED')
!!$        CALL ERROR$STOP('ATOMLIB$BOUNDSTATEWITHHF')
!!$      END IF
!!$      IF(TPR) PRINT*,ITER,'ITERATIONSIN LOOP FOR RADIAL HARTREE FOCK'
!!$      IF(TPR) THEN
!!$        PRINT*,'EOFI ',ITER,E
!!$        CALL ATOMLIB_WRITEPHI('PSI_F.DAT',GID,NR,1,PSI)
!!$      END IF
!!$!   
!!$!     ==========================================================================
!!$!     ==  DETERMINE KINETIC ENERGY                                            ==
!!$!     ==========================================================================
!!$      TPSI(:)=(E-POT(:)*Y0+VFOCK%SCALE*VFOCK%MUX*Y0)*PSI(:) &
!!$     &                    -VFOCK%SCALE*VFOCKPSI(:)+G(:)
!!$      RETURN
!!$      END SUBROUTINE ATOMLIB$BOUNDSTATEWITHHF
