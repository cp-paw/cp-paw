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
      REAL(8)   ,ALLOCATABLE :: drel(:)
      LOGICAL(4)             :: TREL,TZORA
      INTEGER(4)             :: IB,N,L,ISO,IR
      REAL(8)                :: JPHALF ! J+1/2
      REAL(8)                :: C      ! SPEED OF LIGHT
      REAL(8)                :: ALPHA  ! FINE STRUCTURE CONSTANT
      REAL(8)                :: PI     ! PI
      REAL(8)                :: Y0     ! Y0 SPHERICAL HARMONICS FOR L=0
      LOGICAL(4)             :: TFINITENUCLEUS
!     **************************************************************************
      CALL CONSTANTS('C',C)
      ALPHA=1.D0/C
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
!
!     ==========================================================================
!     == THE FOLLOWING PARAMETERS MUST BE SET BY HAND                         ==
!     ==========================================================================
      AEZ=92.D0   ! URANIUM
      TFINITENUCLEUS=.true.
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
      eofi(:)=0.d0
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
        WRITE(*,FMT='(4I5,4F20.5)')IB,LOFI(IB),SOFI(IB),NNOFI(IB),EOFI(IB),SVAR,EOFI(IB)-SVAR
      ENDDO
!
!     ==========================================================================
!     == PREPARE FOR SPIN-ORBIT CALCULATION                                   ==
!     ==========================================================================
      IB=0
      DO N=1,10
        DO L=0,N-1
          DO ISO=-1,1,2  !sign of LS
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
!     == determine small component from analytical large component to test =====
!     == the accuracy of the small component calculation                   =====
!     ==========================================================================
      allocate(drel(nr))
      do ib=1,nb
         aux(:)=0.d0
        CALL SCHROEDINGER$DREL(GID,NR,POT,eofi(ib),DREL)
        CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,LOFI(IB),SOFI(IB) &
     &                                     ,DREL,aux,PHI(:,IB),SPHI(:,IB))
      enddo
!     == compare diracsmallnum.dat with diracsmall.dat =========================
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
      REAL(8)         ,ALLOCATABLE:: PROJ(:)
      REAL(8)         ,ALLOCATABLE:: CPRO(:)
      REAL(8)         ,ALLOCATABLE:: DPRO(:)
      INTEGER(4)                  :: NPRO
      INTEGER(4)                  :: LNX
      INTEGER(4)                  :: LN,LN1,LN2
      INTEGER(4)                  :: IPRO,IPRO1,IPRO2
      REAL(8)                     :: AUX(NR)
      REAL(8)                     :: R(NR)
      REAL(8)                     :: PI,Y0
!     **************************************************************************
      LNX=VPAW%LNX
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
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
  INTEGER(4)         :: NB
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
      DEALLOCATE(VFOCK%L)
      DEALLOCATE(VFOCK%SO)
      DEALLOCATE(VFOCK%F)
      DEALLOCATE(VFOCK%PSI)
      DEALLOCATE(VFOCK%MUX)
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
      INTEGER(4)                  :: LF,LB,LRHO
      INTEGER(4)                  :: LMF,LMB,LMRHO
      INTEGER(4)                  :: IMB,IMRHO
      INTEGER(4)                  :: LX
      REAL(8)                     :: SVAR
      REAL(8)                     :: CG
      REAL(8)                     :: PI
!     **************************************************************************
      IF(LXNEW.LE.LXCOEFF) RETURN
!
      IF(LXCOEFF.NE.-1) DEALLOCATE(COEFF)
      LXCOEFF=LXNEW
      LX=LXNEW
      ALLOCATE(COEFF(LX+1,2*LX+1,LX+1))
      PI=4.D0*ATAN(1.D0)
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
      REAL(8)                     :: PI,Y0
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
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
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
!     ** KEY='START,REL,SO,zora,FOCK='                                        **
!     **                                                                      **
!     **************************************************************************
      USE RADIALFOCK_MODULE
USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID       ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR        ! #(GRID POINTS)
      CHARACTER(*),INTENT(IN)    :: KEY     
      REAL(8)    ,INTENT(INOUT)  :: RBOX      ! BOX RADIUS     
      REAL(8)    ,INTENT(IN)     :: AEZ       ! ATOMIC NUMBER
      INTEGER(4) ,INTENT(IN)     :: NX        ! X#(STATES)
      INTEGER(4) ,INTENT(OUT)    :: NB        ! #(STATES)
      INTEGER(4) ,INTENT(INOUT)  :: LOFI(NX)  ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(INOUT)  :: SOFI(NX)    ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(INOUT)  :: FOFI(NX)     ! OCCUPATION
      INTEGER(4) ,INTENT(INOUT)  :: NNOFI(NX)    ! #(NODES)
      REAL(8)    ,INTENT(OUT)    :: ETOT      ! TOTAL ENERGY
      REAL(8)    ,INTENT(INOUT)  :: POT(NR)   ! POTENTIAL
      TYPE(VFOCK_TYPE),INTENT(INOUT)  :: VFOCK  ! FOCK TERM
      REAL(8)    ,INTENT(OUT)    :: EOFI(NX)  ! ONE-PARTICLE ENERGIES
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,NX)! ONE-PARTICLE WAVE FUNCTIONS
      REAL(8)    ,INTENT(OUT)    :: SPHI(NR,NX) ! SMALL COMPONENT
      REAL(8)   ,PARAMETER       :: TOL=1.D-3
      INTEGER(4),PARAMETER       :: NITER=1000
      LOGICAL(4),PARAMETER       :: TBROYDEN=.TRUE.
      LOGICAL(4),PARAMETER       :: TPR=.true.
      REAL(8)                    :: R(NR)
      REAL(8)                    :: DREL(NR)  ! RELATIVISTIC CORRECTION
      REAL(8)                    :: RHO(NR)
      REAL(8)                    :: G(NR)     ! INHOMOGENEITY
      REAL(8)                    :: AUX(NR),AUX1(NR)   !AUXILIARY ARRAY
      REAL(8)                    :: MUX(NR)   !EXCHANGE ONLY POTENTIAL
      INTEGER(4)                 :: ITER
      REAL(8)                    :: XAV,XMAX,XDE
      REAL(8)                    :: XAVMIN,XMAXMIN,XDEMIN
      INTEGER(4)                 :: NCONV
      LOGICAL(4)                 :: CONVG
      REAL(8)                    :: EKIN,EH,EXC
      REAL(8)                    :: POTIN(NR)
      REAL(8)                    :: SVAR
      INTEGER(4)                 :: I,IB,JB,ISO,L,IR
      INTEGER(4)                 :: ISVAR,IARR(1)
      LOGICAL(4)                 :: TSTART  ! CALCULATE ANGULAR MOMENTA ETC
      LOGICAL(4)                 :: TREL    ! RELATIOVISTIC CALCULATION
      LOGICAL(4)                 :: TSO     ! CALCULATE WITH SPIN ORBIT COUPLING
      LOGICAL(4)                 :: TZORA   ! CALCULATE ZEROTH ORDER RELATIVISTIC CORRECTIONS
      LOGICAL(4)                 :: TFOCK   ! CALCULATE WITH FOCK EXCHANGE
      INTEGER(4)          :: LMAP(19)=(/0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1/)
      INTEGER(4),PARAMETER       :: ZCORES(6)=(/2,10,18,36,54,86/)
      REAL(8)                    :: FTOT
      REAL(8)                    :: PI,Y0,C0LL
      INTEGER(4)                 :: NBROYDENMEM
      REAL(8)                    :: BROYDENSTEP
      INTEGER(4)                 :: LRHOX=4
      INTEGER(4)                 :: IPERIOD,ISVAR1,ipos
      CHARACTER(128)             :: STRING,string1
      REAL(8)                    :: SCALE
      LOGICAL                    :: TSECOND
      REAL(8)                    :: RFOCK !EXTENT OF ORBITALS DEFINING FOCK TERM
      REAL(8)       ,ALLOCATABLE :: EOFI_FOCK(:)
      REAL(8)                    :: POTSAVE(NR,10),POTINSAVE(NR,10)
      REAL(8)                    :: POTOUTSAVE(NR,10)
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
!     == resolve key                                                          ==
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
      DO i=1,7
        IPOS=INDEX(string,',')
        IF(IPOS.ne.0) THEN
          STRING1=STRING(1:IPOS-1)
        ELSE
          STRING1=STRING
        END IF
        STRING=STRING(IPOS+1:)
        TREL=trel.or.trim(STRING1).EQ.'REL' 
        TSO=tso.or.STRING1.EQ.'SO' 
        TZORA=tzora.or.STRING1.EQ.'ZORA' 
        TSTART=tstart.or.STRING1.EQ.'START' 
        IF(STRING1(1:4).EQ.'FOCK') THEN
          TFOCK=.TRUE.
          if(string(5:5).eq.'=')READ(STRING1(6:),*)SCALE
        END IF
        if(ipos.eq.0) exit
      ENDDO
      IF(TPR) THEN
        WRITE(*,FMT='("KEY: ",A)')TRIM(KEY)
        WRITE(*,FMT='("RELATIVISTIC EFFECTS SWITCHED ON? ",L5)')TREL
        WRITE(*,FMT='("SPIN-ORBIT COUPLING SWITCHED ON?  ",L5)')TSO
        WRITE(*,FMT='("ZORA SWITCHED ON?                 ",L5)')TZORA
        WRITE(*,FMT='("FOCK CONTRIBTION ON?              ",L5)')TFOCK
        IF(TFOCK) THEN
          WRITE(*,FMT='("PERCENT FOCK CONTRIBUTION: ",i5)')NINT(SCALE*100.D0)
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
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=Y0
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
                FOFI(IB)=REAL(2*L+2,KIND=8)
                SOFI(IB)=1
                IB=IB+1
                LOFI(IB)=L
                FOFI(IB)=REAL(2*L,KIND=8)
                SOFI(IB)=-1
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
              FOFI(IB)=REAL(MIN(2*L+2,ISVAR1),KIND=8)
              SOFI(IB)=1
              IB=IB+1
              LOFI(IB)=L
              FOFI(IB)=REAL(MAX(0,MIN(2*L,ISVAR1-2*L-2)),KIND=8)
              SOFI(IB)=-1
            ELSE
              SOFI(IB)=1
            END IF
          END IF
        ENDDO
        NB=IB
        FTOT=SUM(FOFI(:NB))
        IF(ABS(FTOT-AEZ).GT.1.D-8) THEN
          DO IB=1,NB
            WRITE(*,'("IB=",I2," L=",I1," SOFI=",I2," F=",F8.2," SUM(F)=",F8.2)') &
        &           IB,LOFI(IB),SOFI(IB),FOFI(IB),SUM(FOFI(:IB))
          ENDDO
          CALL ERROR$MSG('INCONSISTENT NUMBER OF ELECTRONS')
          CALL ERROR$R8VAL('AEZ ',AEZ)
          CALL ERROR$R8VAL('#(ELECTRONS) ',FTOT)
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
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,EKIN)
          CALL ATOMLIB$BOXVOFRHO(GID,NR,RBOX,AEZ,RHO,POT,EH,EXC)
          ETOT=EKIN+EH+EXC
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
CALL PERIODICTABLE$GET(NINT(AEZ),'R(COV)',RFOCK)
RFOCK=1.1D0*RFOCK
          IF(.NOT.TSECOND) THEN
            ALLOCATE(EOFI_FOCK(NB))
            EOFI_FOCK(:)=EOFI(:)
            CALL ATOMLIB$BOUNDSTATES(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI_FOCK,POT &
     &                            ,TREL,TZORA,RFOCK,PHI,SPHI)
          ELSE
            CALL ATOMLIB$BOUNDSTATESWITHHF(GID,NR,NB,LOFI,SOFI,NNOFI,EOFI,POT &
     &                              ,VFOCK,TREL,TZORA,RFOCK,PHI,SPHI)
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
!       ==  GENERATE NEXT ITERATION USING D. G. ANDERSON'S METHOD             ==
!       ========================================================================
        XAV=SQRT(SUM(R**3*(POT-POTIN)**2)/SUM(R**3))
        XMAX=MAXVAL(ABS(POT-POTIN))
        AUX(:)=R**2*RHO*(POT-POTIN) 
        AUX(:)=ABS(AUX)
        CALL RADIAL$INTEGRAL(GID,NR,AUX,XDE)
        NCONV=NCONV+1
        IF(XAV.LT.XAVMIN) THEN
          XAVMIN=XAV
          NCONV=0
        END IF
        IF(XMAX.LT.XMAXMIN) THEN
          XMAXMIN=XMAX
          NCONV=0
        END IF
        IF(XDE.LT.XDEMIN) THEN
          XDEMIN=XDE
          NCONV=0
        END IF
        IF(TPR)PRINT*,ITER,' AV(POT-POTIN)=',XAV,' MAX:(POT-POTIN)=',XMAX &
      &                   ,' DE ',XDE 
PRINT*,ITER,' AV(POT-POTIN)=',XAV,' MAX:(POT-POTIN)=',XMAX &
&                   ,' DE ',XDE 
!CALL ATOMLIB_WRITEPHI('AUX',GID,NR,1,AUX)
!CALL ATOMLIB_WRITEPHI('RHO',GID,NR,1,RHO)
IF(TPR) THEN
  DO I=10,2,-1
    POTSAVE(:,I)=POTSAVE(:,I-1)
    POTINSAVE(:,I)=POTINSAVE(:,I-1)
    POTOUTSAVE(:,I)=POTOUTSAVE(:,I-1)
  ENDDO  
  POTSAVE(:,1)=POT-POTIN
  POTINSAVE(:,1)=POTIN
  POTOUTSAVE(:,1)=POT
  DO I=2,10
     POTOUTSAVE(:,I)=POTOUTSAVE(:,I)-POTOUTSAVE(:,1)
     POTINSAVE(:,I)=POTINSAVE(:,I)-POTINSAVE(:,1)
  ENDDO
  DO I=1,10
     POTINSAVE(:,I)=R(:)*POTINSAVE(:,I)
     POTSAVE(:,I)=R(:)*POTSAVE(:,I)
  ENDDO
  CALL ATOMLIB_WRITEPHI('POT-POTIN',GID,NR,10,POTSAVE)
  CALL ATOMLIB_WRITEPHI('POTIN',GID,NR,9,POTINSAVE(:,2:10))
  CALL ATOMLIB_WRITEPHI('POTOUT',GID,NR,9,POTOUTSAVE(:,2:10))
END IF
        CALL BROYDEN$STEP(NR,POTIN,POT-POTIN)
        POT=POTIN
        CONVG=(XMAX.LT.TOL).AND.NCONV.GT.5
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
        GOTO 1000
      END IF
!
!CALL RADIALFOCK$PRINTVFOCK('VFOCK.DAT',VFOCK)
!DO I=1,NB
!  WRITE(*,FMT='(3I4,F10.2,I5,F20.3)')I,LOFI(I),SOFI(I),FOFI(I),NNOFI(I),EOFI(I)
!ENDDO
!PRINT*,'#ITERATIONS ',ITER
!STOP
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
        CALL ATOMLIB$BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RBOX,TVARDREL,DREL,G &
     &                         ,NNOFI(IB),POT,EOFI(IB),PHI(:,IB))
!
!       ========================================================================
!       == DETERMINE SMALL COMPONENTS                                         ==
!       ========================================================================
        IF(TREL.AND.(.NOT.TZORA)) THEN
          G(:)=0.D0
          CALL SCHROEDINGER$DREL(GID,NR,POT,eofi(ib),DREL)
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
      REAL(8)       :: TPHI(NR,NB)
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
        CALL ATOMLIB$BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RBOX,TVARDREL,DREL,G &
     &                         ,NNOFI(IB),POT,EOFI(IB),PHI(:,IB))
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
      SUBROUTINE ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RBOX,TVARDREL,DREL,G,NN,POT,E,PHI)
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
      REAL(8)    ,INTENT(IN)     :: RBOX    ! BOX RADIUS
      LOGICAL(4) ,INTENT(IN)     :: TVARDREL! UPDATE RELATIVISTIC PARAMETER
      REAL(8)    ,INTENT(INOUT)  :: DREL(NR)!RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      INTEGER(4) ,INTENT(IN)     :: NN      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
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
      REAL(8)                    :: X1,X2,Y1,Y2
      REAL(8)                    :: VAL
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
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,ROUT,Z0)
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
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL1,SO,GHOM,L,E,IDIR,PHIINHOM)
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
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)    ,PARAMETER      :: TOL=1.D-12
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHI1(NR),PHI2(NR)
      REAL(8)                    :: DREL(NR),GHOM(NR),PHIHOM(NR)
      REAL(8)                    :: PHIINHOM(NR)
      INTEGER(4) ,PARAMETER      :: NITER=100
      REAL(8)                    :: VAL1,VAL2,SVAR
      INTEGER(4)                 :: I,IR,IRMATCH,SO,IDIR
      LOGICAL                    :: THOM
      REAL(8)                    :: Y2,Y1,X2,X1,VAL,DER,DERO,R1,R2
      INTEGER(4)                 :: IRBOX
!     **************************************************************************
                                 CALL TRACE$PUSH('ATOMLIB$PAWBOUNDSTATE')
      CALL RADIAL$R(GID,NR,R)
!     ==  R(IRBOX) IS THE FIRST GRIDPOINT JUST OUTSIDE THE BOX =================
      IRBOX=1
      DO IR=1,NR-2
        IRBOX=IR
        IF(R(IR).GE.RBOX) EXIT
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
      DO I=1,NITER
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
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RBOX,Z0)
        Z0=Z0-REAL(NN+1,KIND=8)
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
!
      SVAR=VAL2-VAL1
      VAL1=VAL1/SVAR
      VAL2=VAL2/SVAR
      PHI=PHI1*VAL2-PHI2*VAL1
!
!     ==========================================================================
!     ==  DETERMINE MATCHING POINT                                            ==
!     ==========================================================================
      DO IR=1,NR
        IRMATCH=IR
        IF(SUM(ABS(PRO(IR,:))).LT.1.D-8) THEN
          IF(ABS(PHI2(IR)-PHI1(IR)).GT.1.D-3*ABS(PHI1(IR))) EXIT
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
        CALL RADIAL$DERIVATIVE(GID,NR,PHI,R(IRMATCH),DER)
        CALL RADIAL$DERIVATIVE(GID,NR,PHIINHOM,R(IRMATCH),DERO)
        SVAR=(DERO-DER)/PHI(IRMATCH)
        PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
      END IF
 

                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$BOXVOFRHO(GID,NR,RAD,AEZ,RHO,POT,EH,EXC)
!     ******************************************************************
!     **                                                              **
!     **  ELECTROSTATIC AND EXCHANGE-CORRELATION POTENTIAL            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8),   INTENT(IN) :: AEZ
      REAL(8),   INTENT(IN) :: RAD
      REAL(8),   INTENT(IN) :: RHO(NR)
      REAL(8),   INTENT(OUT):: POT(NR)
      REAL(8),   INTENT(OUT):: EH
      REAL(8),   INTENT(OUT):: EXC
      REAL(8),   PARAMETER  :: RHOMIN=1.D-2
      REAL(8)               :: POTH(NR)
      REAL(8)               :: POTXC(NR)
      REAL(8)               :: PI
      REAL(8)               :: FOURPI
      REAL(8)               :: Y0
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
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
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
      REAL(8)               :: PI
      REAL(8)               :: FOURPI
      REAL(8)               :: Y0
      REAL(8)               :: AUX(NR),AUX1(NR)
      REAL(8)               :: GRHO(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VGXC,VXC,EXC1,RH,GRHO2
      REAL(8)               :: DUMMY1,DUMMY2,DUMMY3
      INTEGER(4)            :: IR,IRBOX
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
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
      INTEGER(4)              :: L
      REAL(8)                 :: E
      REAL(8)                 :: G(NR)
      REAL(8)                 :: DREL(NR)
      INTEGER(4)              :: NN
      INTEGER(4)              :: SO
      INTEGER(4)              :: IB,I
      REAL(8)                 :: PI,Y0
      LOGICAL(4)              :: TVARDREL
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
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
        CALL ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RBOX,TVARDREL,DREL,G,NN,POT,E,UOFI(:,IB))
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
      REAL(8)         ,PARAMETER  :: TOL=1.D-7  !(1D-12 IS TOO SMALL)
!      REAL(8)         ,PARAMETER  :: TOLMAX=1.D-8  !EMERGENCY TOLERANCE
      INTEGER(4)      ,PARAMETER  :: NITER=200
      REAL(8)         ,PARAMETER  :: XMAX=1.D+15 !MAX. FACTOR FOR WAVEFUNCTION
      REAL(8)                  :: PHASE
      REAL(8)                  :: DPSI(NR)
      REAL(8)                  :: PHIDOT(NR)
      REAL(8)                  :: PHIPRIME(NR)
      REAL(8)                  :: PHIHOM(NR)
      REAL(8)                  :: POT1(NR)
      REAL(8)                  :: VALDOT,DERDOT,VALPRIME,DERPRIME,VAL0,DER0
      REAL(8)                  :: VALHOM,DERHOM,WHOM
      REAL(8)                  :: W0,WPRIME,WDOT
      REAL(8)                  :: Y1,Y2,X1,X2
      REAL(8)                  :: PI,Y0
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
      REAL(8)                  :: DEV,DEVPREV
      LOGICAL(4)               :: TPR=.false.
      LOGICAL(4)               :: TFIXE ! CONSTANT ENERGY CALCULATION IF TRUE
      LOGICAL(4)               :: TINHOM ! INHOMOGENEOUS EQUATION
      INTEGER(4)               :: NN
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
      TFIXE=RBND.LE.0.D0
      TINHOM=MAXVAL(ABS(G)).GT.0.D0
!
!     :: REMARK: SCHROEDINGER$PHASESHIFT TAKES VALUE AND DERIVATIVE FROM A    ::
!     :: TWO-POINT FORMULA, AND THEREFORE IS NOT FULLY CONSISTENT WITH        ::
!     :: RADIAL$VALUE AND RADIAL$DERIVATIVE                                   ::
      CALL SCHROEDINGER$PHASESHIFT(GID,NR,PSI,RBND,PHASE)
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
!       == HARTREE-FOCK CORRECTION =============================================
        CALL RADIALFOCK$VPSI(GID,NR,VFOCK,L,PSI,AUX1)
        AUX(:)=AUX(:)+AUX1(:)
!       == CORRECT GLITCHES ====================================================
        AUX(1:2)=0.D0
!       == ESTIMATE ENERGY =====================================================
        IF(.NOT.TFIXE) THEN
          AUX1(:)=R(:)**2*PSI(:)*(AUX(:)-G(:))
          CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
          CALL RADIAL$VALUE(GID,NR,AUX2,RBND,E)
          AUX1(:)=R(:)**2*PSI(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
          CALL RADIAL$VALUE(GID,NR,AUX2,RBND,NORM)
          E=E/NORM
          IF(E.LT.-5.D+3) THEN
            CALL ERROR$MSG('ENERGY BELOW -5000 H ENCOUNTERED.')
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
!       == SUBTRACT INHOMOGENEITY ==============================================
        AUX(:)=AUX(:)-G(:)
!       == MAP INTO DPSI
        DPSI(:)=-AUX(:)
        DPSI(1:2)=0.D0
!
!       ========================================================================
!       ==  DETERMINE PHIPRIME-PHIDOT*EPRIME                                  ==
!       ========================================================================
        CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT,E,IRCL,IROUT)

!       == BOUNDARY CONDITION PHI(ROUT)=0 ====================================
        ROUT=R(IROUT)
!       ==  SET KINETIC ENERGY TO ZERO BEYOND ROUT TO AVOID AN OVERFLOW ======
!       == ATTENTION: THERE IS A STEP IN THE POTENTIAL =======================
        POT1(:)=POT(:)
        POT1(IROUT:)=E
!
        G1=DPSI(:)
        G1(1:2)=0.D0
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G1,L,E,1,PHIPRIME)
        IF(.NOT.TFIXE) THEN
          G1=PSI(:)
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G1,L,E,1,PHIDOT)
!
!         == DETERMINE VALUES AND DERIVATIVES CONSISTENT WITH ==================
!         == SCHROEDINGER$PHASESHIFT ===========================================
          CALL SCHROEDINGER$RESOLVEPHASESHIFT(PHASE,VAL,DER,NN)!
!
          IF(IROUT.LT.IRBND) THEN
            IREND=IROUT
            VAL=0.D0
            DER=1.D0
            REND=R(IREND)
          ELSE
            IREND=IRBND
            REND=RBND
          END IF
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
!
        ELSE
          DPSI(:)=PHIPRIME(:)
        END IF
!   
!       ========================================================================
!       ==  CHECK CONVERGENCE                                                 ==
!       ========================================================================
        VAL=MAXVAL(ABS(DPSI(:IRBND-3)))/SQRT(NORM)
        IF(TPR) PRINT*,'ITER=',ITER,' L=',L,' E=',E,' DEV=',VAL
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
!
      IF(.NOT.CONVG) THEN
        CALL ATOMLIB_WRITEPHI('PSI.DAT',GID,NR,1,PSI)
        CALL ATOMLIB_WRITEPHI('DPSI.DAT',GID,NR,1,DPSI)
        CALL ERROR$MSG('LOOP FOR RADIAL HARTREE FOCK NOT CONVERGED')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('E',E)
        CALL ERROR$R8VAL('REND',REND)
        CALL ERROR$R8VAL('RBND',RBND)
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
    &                        ,NB,LOFI,SOFI,FOFI,PSI,ETOT)
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
      REAL(8)               :: DRRPSI(NR,NB)
      REAL(8)               :: AUX(NR),AUX1(NR),AUX2(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: EOFI(NB)
      REAL(8)               :: DREL(NR),RDPRIME(NR)
      REAL(8)               :: VAL,DER
      REAL(8)               :: PI,Y0
      INTEGER(4)            :: IRBOX
      REAL(8)               :: SOFACTOR
      INTEGER(4)            :: IR,IB
      INTEGER(4)            :: L
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
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
!       == POTENTIAL ENERGY ==================================================
        AUX1(:)=AUX1(:)+R(:)**2*POT(:)*Y0*PSI(:,IB)**2
!       == HARTREE-FOCK CORRECTION ===========================================
        CALL RADIALFOCK$VPSI(GID,NR,VFOCK,LOFI(IB),PSI(:,IB),AUX2)
        AUX2(:)=PSI(:,IB)*(0.5D0*AUX2-0.5D0*VFOCK%MUX(:)*Y0*PSI(:,IB))
        AUX1(:)=AUX1(:)+R(:)**2*AUX2(:)
!       == ADD TO SUM OVER STATES ============================================
        AUX(:)=AUX(:)+FOFI(IB)*AUX1(:)
      ENDDO
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,ETOT)
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
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RBND,Z0)
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
      CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RBND,Z0)
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
      REAL(8)                 :: PHI1(NR,NPHI)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      PHI1=PHI
      DO IR=1,NR
        IF(R(IR).GT.3.D0.AND.MAXVAL(ABS(PHI1(IR,:))).GT.1.D+3) EXIT
        DO I=1,NPHI  ! AVOID CONFLICT WITH XMGRACE
          IF(ABS(PHI1(IR,I)).LT.1.D-60) PHI1(IR,I)=0.D0
        ENDDO
        WRITE(100,FMT='(F15.10,2X,20(E25.10,2X))')R(IR),PHI1(IR,:)
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
!!$        CALL ATOMLIB$BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RBOX,DREL,G,NN &
!!$     &                         ,POT,EOFI(IB),PSI(:,IB))
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
!!$      CALL SCHROEDINGER$PHASESHIFT(GID,NR,PSI,RBND,PHASE)
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
!!$! CALL SCHROEDINGER$PHASESHIFT(GID,NR,PSI,RBND,SVAR)
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
!!$      CALL ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RBOX,DREL,G,NN,POT,E,PSI)
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
