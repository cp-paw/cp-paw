!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                   **
!**  INTERFACES TO THE SCHROEDINGER ROUTINES                          **
!**                                                                   **
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIAL$SCHRODINGER(GID,NR,POT,DREL,SO,G,L,E,IDIR,PHI)
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! MAINANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
                 ! SO=0: NO SO; SO=1: L/S PARALLEL; SO=-1: L,S ANTIPARALLEL
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
!                                           ! DREL= M0/MREL(R)-1
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: E       !ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: IDIR    ! IDIR=1 INTEGRATE OUTWARD
                                            ! IDIR=-1 INTEGRATE INWARD
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2007 ************
      CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,E,IDIR,PHI)
      RETURN
      END SUBROUTINE RADIAL$SCHRODINGER
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIAL$NONSPHBOUND_NONSO(GID,NR,LX,LMX,LMRX,TMAINSH &
     &                                   ,POT,DREL,G,ENU,NPHI,EB,PHI,TPHI,TOK)
!     **                                                                  **
!     **  SOLVES THE RELATIVISTIC RADIAL DIRAC EQUATION FOR THE           **
!     **  LARGE COMPONENT. DREL=1/MREL-1/M0 IS A MEASURE FOR THE          **
!     **  RELATIVISTIC EFFECTS, WHERE MREL=M0+(E-V)/2C^2 AND M0 IS THE    **
!     **  REST MASS. V=POT*Y0 IS THE POTENTIAL.                           **
!     **  SPIN ORBIT COUPLING IS MEASURED BY SO WHICH CAN HAVE THE VALUES:**
!     **    SO=0 NO-SPIN ORBIT COUPLING                                   **
!     **    SO=1 PARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM               **
!     **  G IS AN INHOMOGENEITY, WHICH MUST BE SET TO ZERO FOR THE        **
!     **  HOMOGENEOUS SOLUTION.                                           **
!     **                                                                  **
!     **  THE SOLUTIONS ARE CALCULATED TO LINEAR ORDER IN DE              **
!     **  WHERE E+DE IS THE NEW ENERGY OF THE WAVE FUNCTIONS              **
!     **                                                                  **
!     **  THE DIFFERENTIAL EQUATION IS SOLVED WITH THE VERLET ALGORITHM   **
!     **    DPHI/DX=(PHI(+)-PHI(-))/2                                     **
!     **    D2PHI/DX2=PHI(+)-2PHI(0)+PHI(-)                               **
!     **  WHERE X IS THE VARIABLE WITH R(X=I)=R_I                         **
!     **                                                                  **
!     **  THE NONRELATIVISTIC SOLUTION IS OBTAINED BY SETTINH DREL=0      **
!     **                                                                  **
!     **  ATTENTION! THE ROUTINE IS NOT GUARDED AGAINST OVERFLOW DUE      **
!     **    TO THE EXPONENTIAL INCREASE OF THE SOLUTION                   **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID-ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: LX      ! X(ANGULAR MOMENTA FOR WAVEF.)
      INTEGER(4) ,INTENT(IN)     :: LMX     ! X#(WAVE FUNCTION ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)     :: LMRX    ! X#(POTENTIAL ANGULAR MOMENTA)
      LOGICAL(4) ,INTENT(IN)     :: TMAINSH(LX+1)
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
!                                           ! DREL= M0/MREL(R)-1
      REAL(8)    ,INTENT(IN)     :: G(NR,LMX)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: ENU     ! EXPANSION ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR,LMRX) !POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: NPHI              ! #(WAVE FUNCTIONS)
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,LMX,NPHI)  ! WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: TPHI(NR,LMX,NPHI) ! P**2/(2M)*WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: EB(NPHI)          ! ONE-PARTICKE ENERGIES
      LOGICAL(4) ,INTENT(OUT)    :: TOK               ! ERROR FLAG
      INTEGER(4)                 :: L,ISVAR
!     **********************************************************************
      IF(LMX.NE.(LX+1)**2) THEN
        CALL ERROR$MSG('LMX AND LX ARE INCONSISTENT')
        CALL ERROR$STOP('RADIAL$NONSPHBOUND_NONSO')
      END IF
      ISVAR=0
      DO L=0,LX
        IF(TMAINSH(L+1)) ISVAR=ISVAR+2*L+1
      ENDDO
      IF(ISVAR.NE.NPHI) THEN
        CALL ERROR$MSG('TMAINSH AND NPHI ARE INCONSISTENT')
        CALL ERROR$STOP('RADIAL$NONSPHBOUND_NONSO')
      END IF
      CALL SCHROEDINGER$LBND_SCALREL(GID,NR,LX,LMX,LMRX,TMAINSH,POT,DREL,G,ENU &
     &                             ,NPHI,EB,PHI,TPHI,TOK)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIAL$NONSPHBOUND_NSO_SLOC_OV(GID,NR,LX,LMX,LMRX,TMAINSH &
     &                     ,POT,LXSL,POTSL,OVSL,G,ENU,NPHI,EB,PHI,TPHI,TOK)
!     **                                                                      **
!     **  SOLVES THE NONRELATIVISTIC RADIAL SCHROEDINGER EQUATION FOR THE     **
!     **  G IS AN INHOMOGENEITY, WHICH MUST BE SET TO ZERO FOR THE            **
!     **  HOMOGENEOUS SOLUTION.                                               **
!     **                                                                      **
!     **  THE SOLUTIONS ARE CALCULATED TO LINEAR ORDER IN DE                  **
!     **  WHERE E+DE IS THE NEW ENERGY OF THE WAVE FUNCTIONS                  **
!     **                                                                      **
!     **  THE DIFFERENTIAL EQUATION IS SOLVED WITH THE VERLET ALGORITHM       **
!     **    DPHI/DX=(PHI(+)-PHI(-))/2                                         **
!     **    D2PHI/DX2=PHI(+)-2PHI(0)+PHI(-)                                   **
!     **  WHERE X IS THE VARIABLE WITH R(X=I)=R_I                             **
!     **                                                                      **
!     **  ATTENTION! THE ROUTINE IS NOT GUARDED AGAINST OVERFLOW DUE          **
!     **    TO THE EXPONENTIAL INCREASE OF THE SOLUTION                       **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID-ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: LX      ! X#(MAIN ANGULAR MOMENTA OF WAVEF)
      INTEGER(4) ,INTENT(IN)     :: LMX     ! X#(WAVE FUNCTION ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)     :: LMRX    ! X#(POTENTIAL ANGULAR MOMENTA)
      LOGICAL(4) ,INTENT(IN)     :: TMAINSH(LX+1)
      INTEGER(4) ,INTENT(IN)     :: LXSL    ! #(SEMILOCAL POTENIALS, ONE PER L)
      REAL(8)    ,INTENT(IN)     :: POTSL(NR,LXSL)  ! SEMI-LOCAL POTENTIALS
      REAL(8)    ,INTENT(IN)     :: OVSL(NR,LXSL)   ! SEMI-LOCAL OVERLAP CONTRIBTION
      REAL(8)    ,INTENT(IN)     :: G(NR,LMX)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: ENU     ! EXPANSION ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR,LMRX) !POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: NPHI              ! #(WAVE FUNCTIONS)
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,LMX,NPHI)  ! WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: TPHI(NR,LMX,NPHI) ! P**2/(2M)*WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: EB(NPHI)          ! ONE-PARTICKE ENERGIES
      LOGICAL(4) ,INTENT(OUT)    :: TOK               ! ERROR FLAG
!     **************************************************************************
      CALL SCHROEDINGER$LBND_SLOC(GID,NR,LX,LMX,LMRX,TMAINSH,POT,LXSL,POTSL &
     &                                       ,OVSL,G,ENU,NPHI,EB,PHI,TPHI,TOK)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIAL$NONSPHBOUND(GID,NR,NDIMD,LMX,LMRX,POT,DREL,G,ENU &
     &                             ,NPHI,EB,PHI,TPHI,SPHI,TSPHI,TOK)
!     **                                                                  **
!     **  SOLVES THE RELATIVISTIC RADIAL DIRAC EQUATION FOR THE           **
!     **  LARGE COMPONENT. DREL=1/MREL-1/M0 IS A MEASURE FOR THE          **
!     **  RELATIVISTIC EFFECTS, WHERE MREL=M0+(E-V)/2C^2 AND M0 IS THE    **
!     **  REST MASS. V=POT*Y0 IS THE POTENTIAL.                           **
!     **  SPIN ORBIT COUPLING IS MEASURED BY SO WHICH CAN HAVE THE VALUES:**
!     **    SO=0 NO-SPIN ORBIT COUPLING                                   **
!     **    SO=1 PARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM               **
!     **  G IS AN INHOMOGENEITY, WHICH MUST BE SET TO ZERO FOR THE        **
!     **  HOMOGENEOUS SOLUTION.                                           **
!     **                                                                  **
!     **  THE SOLUTIONS ARE CALCULATED TO LINEAR ORDER IN DE              **
!     **  WHERE E+DE IS THE NEW ENERGY OF THE WAVE FUNCTIONS              **
!     **                                                                  **
!     **  THE DIFFERENTIAL EQUATION IS SOLVED WITH THE VERLET ALGORITHM   **
!     **    DPHI/DX=(PHI(+)-PHI(-))/2                                     **
!     **    D2PHI/DX2=PHI(+)-2PHI(0)+PHI(-)                               **
!     **  WHERE X IS THE VARIABLE WITH R(X=I)=R_I                         **
!     **                                                                  **
!     **  THE NONRELATIVISTIC SOLUTION IS OBTAINED BY SETTINH DREL=0      **
!     **                                                                  **
!     **  ATTENTION! THE ROUTINE IS NOT GUARDED AGAINST OVERFLOW DUE      **
!     **    TO THE EXPONENTIAL INCREASE OF THE SOLUTION                   **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: NDIMD   ! (1,2,4)#(SPINOR COMPONENTS)
      INTEGER(4) ,INTENT(IN)     :: LMX     ! X#(WAVE FUNCTION ANGULAR MOMENTA
      INTEGER(4) ,INTENT(IN)     :: LMRX    ! X#(POTENTIAL ANGULAR MOMENTA)
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
!                                           ! DREL= M0/MREL(R)-1
      REAL(8)    ,INTENT(IN)     :: G(NR,LMX,2)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: ENU       !ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR,LMRX,NDIMD) !POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: NPHI
      COMPLEX(8) ,INTENT(OUT)    :: PHI(NR,LMX,2,NPHI) ! WAVE-FUNCTION
      COMPLEX(8) ,INTENT(OUT)    :: TPHI(NR,LMX,2,NPHI) ! WAVE-FUNCTION
      COMPLEX(8) ,INTENT(OUT)    :: SPHI(NR,LMX,2,NPHI) ! WAVE-FUNCTION
      COMPLEX(8) ,INTENT(OUT)    :: TSPHI(NR,LMX,2,NPHI) ! KINETIC ENERGY * WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: EB(NPHI)
      LOGICAL(4) ,INTENT(OUT)    :: TOK
!     ************************************************************************
      CALL SCHROEDINGER$LBND_FULLYREL(GID,NR,NDIMD,LMX,LMRX,POT,DREL,G,ENU &
     &                             ,NPHI,EB,PHI,TPHI,SPHI,TSPHI,TOK)
      RETURN
      END SUBROUTINE RADIAL$NONSPHBOUND
!
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!**                                                                           **
!**  BEGINNING OF THE ACTUAL SCHROEDINGER OBJECT                              **
!**                                                                           **
!**    SCHROEDINGER$SPHERICAL                                                 **
!**                                                                           **
!**    THE LBND ROUTINES DETERMINE A SET OF BOUND STATES IN A NON-SPHERICAL   **
!**    POTENTIAL USING THE LINEAR ENERGY EXPANSION OF THE WAVE FUNCTION       **
!**    ABOUT AN ENERGY ENU                                                    **
!**                                                                           **
!**      SCHROEDINGER$LBND_SLOC     NONRELATIVASTIC WITH SEMI-LOCAL POT       **
!**      SCHROEDINGER$LBND_SCALREL  SCALAR RELATIVISTIC (REAL PHI)            **
!**      SCHROEDINGER$LBND_FULLYREL INCLUDES SPIN-ORBIT COUPLING (COMPLEX PHI)**
!**                                                                           **
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,E,IDIR,PHI)
!     **                                                                      **
!     **  SOLVES THE RELATIVISTIC RADIAL DIRAC EQUATION FOR THE               **
!     **  LARGE COMPONENT. DREL=1/MREL-1/M0 IS A MEASURE FOR THE              **
!     **  RELATIVISTIC EFFECTS, WHERE MREL=M0+(E-V)/2C^2 AND M0 IS THE        **
!     **  REST MASS. V=POT*Y0 IS THE POTENTIAL.                               **
!     **  SPIN ORBIT COUPLING IS MEASURED BY SO WHICH CAN HAVE THE VALUES:    **
!     **    SO=0 NO-SPIN ORBIT COUPLING                                       **
!     **    SO=1 PARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM                   **
!     **    SO=-1 ANTIPARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM              **
!     **  IDIR DESCRIBES THE DIRECTION OF THE INTEGRATION OF THE DIFF.EQ.     **
!     **    IDIR=1 OUTWARD INTEGRATION                                        **
!     **    IDIR=-1 INWARD INTEGRATION                                        **
!     **  G IS AN INHOMOGENEITY, WHICH MUST BE SET TO ZERO FOR THE            **
!     **  HOMOGENEOUS SOLUTION.                                               **
!     **                                                                      **
!     **  THE DIFFERENTIAL EQUATION IS SOLVED WITH THE VERLET ALGORITHM       **
!     **    DPHI/DX=(PHI(+)-PHI(-))/2                                         **
!     **    D2PHI/DX2=PHI(+)-2PHI(0)+PHI(-)                                   **
!     **  WHERE X IS THE VARIABLE WITH R(X=I)=R_I                             **
!     **                                                                      **
!     **  IN THE PRESENCE OF AN INHOMOGENEITY, THE SOLUTION STARTS            **
!     **  WITH ZERO VALUE AND DERIVATIVE.                                     **
!     **                                                                      **
!     **  IN THE ABSENCE OF AN INHOMOGENEITY, THE SOLUTION STARTS             **
!     **  WITH R**L FROM THE INSIDE AND FROM THE OUTSIDE WITH VALUE ZERO      **
!     **  AND FINITE SLOPE                                                    **
!     **                                                                      **
!     **  THE NONRELATIVISTIC SOLUTION IS OBTAINED BY SETTING DREL=0          **
!     **                                                                      **
!     **  ATTENTION! THE ROUTINE IS NOT GUARDED AGAINST OVERFLOW DUE          **
!     **    TO THE EXPONENTIAL INCREASE OF THE SOLUTION                       **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **  - POT IS ONLY THE RADIAL PART OF THE POTENTIAL.                     **
!     **    THE POTENTIAL IS POT*Y0 WHERE Y0 IS A SPHERICAL HARMONIC          **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! MAINANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
                 ! SO=0: NO SO; SO=1: L/S PARALLEL; SO=-1: L,S ANTIPARALLEL
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
!                                           ! DREL= M0/MREL(R)-1
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: E       !ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: IDIR    ! IDIR=1 INTEGRATE OUTWARD
                                            ! IDIR=-1 INTEGRATE INWARD
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      REAL(8)    ,PARAMETER      :: PI=4.D0*ATAN(1.D0)
      REAL(8)    ,PARAMETER      :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                    :: A(NR) 
      REAL(8)                    :: B(NR) 
      REAL(8)                    :: C(NR) 
      REAL(8)                    :: D(NR) 
      REAL(8)                    :: R(NR) 
      REAL(8)                    :: SOFACTOR
      REAL(8)                    :: RDPRIME(NR)
      LOGICAL(4)                 :: THOM
INTEGER(4) :: IR
!     **************************************************************************
!
!     ==========================================================================
!     == SPIN ORBIT COUPLING                                                  ==
!     ==========================================================================
      IF (SO.EQ.0) THEN
        SOFACTOR=0.D0
      ELSE IF(SO.EQ.1) THEN
        SOFACTOR=REAL(-L,KIND=8)       ! PARALLEL SPIN AND ORBIT
      ELSE IF(SO.EQ.-1) THEN
        SOFACTOR=REAL(L+1,KIND=8)        ! ANTIPARALLELSPIN AND ORBIT
      ELSE
         CALL ERROR$MSG('SO CAN ONLY HAVE VALUES -1,0,1')
         CALL ERROR$STOP('RADIAL$SCHRODINGER')
      END IF
!
!     ==========================================================================
!     == SET UP DIFFERENTIAL EQUATION                                         ==
!     ==========================================================================
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
      CALL RADIAL$R(GID,NR,R)
      A(:)=1.D0+DREL(:)
!     == AVOID DIVIDE BY ZERO IF THE FIRST GRID POINT IS THE ORIGIN.
!     == THE FORCES ON THE FIRST GRID POINT ARE NOT USED,
!     == BECAUSE RADIAL$DGL IS BASED ON THE VERLET ALGORITHM
!     == THAT CANNOT USE THE FORCES ON THE FIRST AND LAST GRID POINT
      B(2:)=2.D0*(1.D0+DREL(2:))/R(2:)+RDPRIME(2:)
      C(2:)=-(1.D0+DREL(2:))*REAL(L*(L+1),KIND=8)/R(2:)**2 &
     &    +RDPRIME(2:)*SOFACTOR/R(2:) &
     &    -2.D0*(POT(2:)*Y0-E)
      B(1)=B(2)
      C(1)=C(2)
      D(:)=-2.D0*G(:)
!
!     ==========================================================================
!     == BOUNDARY CONDITIONS                                                  ==
!     ==========================================================================
      THOM=MAXVAL(ABS(G(:))).EQ.0.D0
      IF(IDIR.GE.0) THEN
        PHI(1:2)=0.D0
        IF(THOM)PHI(1:2)=R(1:2)**L
      ELSE
        PHI(NR-1:NR)=0.D0
        IF(THOM)PHI(NR-1)=1.D-8
      END IF
!
!     ==========================================================================
!     == SOLVE DIFFERENTIAL EQUATION                                          ==
!     ==========================================================================
      CALL RADIAL$DGL(GID,IDIR,NR,A,B,C,D,PHI)
DO IR=1,NR
  IF(.NOT.(PHI(IR).GT.0.D0.OR.PHI(IR).LE.0.D0)) THEN
    PRINT*,'A(1:5) ',A(1:5)
    PRINT*,'B(1:5) ',B(1:5)
    PRINT*,'C(1:5) ',C(1:5)
    PRINT*,'D(1:5) ',D(1:5)
    PRINT*,'DREL(1:5) ',DREL(1:5)
    PRINT*,'PHI(1:5) ',PHI(1:5)
    PRINT*,'A(NR-5:NR) ',A(NR-5:NR)
    PRINT*,'E',E
    CALL ERROR$MSG('NANS ENCOUNTERED. SHUTTING DOWN...')
    CALL ERROR$MSG('NAN="NOT-A-NUMBER"')
    CALL ERROR$I4VAL('IR',IR)
    CALL ERROR$R8VAL('PHI(IR)',PHI(IR))
    CALL ERROR$I4VAL('L',L)
    CALL ERROR$I4VAL('SO',SO)
    CALL ERROR$R8VAL('E',E)
    CALL ERROR$I4VAL('IDIR',IDIR)
    CALL ERROR$STOP('SHROEDINGER$SPHERICAL')
 END IF
ENDDO
      RETURN
      END SUBROUTINE SCHROEDINGER$SPHERICAL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER$EKINPSI(GID,NR,POT,DREL,SO,L,PHI,TPHI)
!     **************************************************************************
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **  - POT IS ONLY THE RADIAL PART OF THE POTENTIAL.                     **
!     **    THE POTENTIAL IS POT*Y0 WHERE Y0 IS A SPHERICAL HARMONIC          **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2009 ************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! MAINANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
                 ! SO=0: NO SO; SO=1: L/S PARALLEL; SO=-1: L,S ANTIPARALLEL
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
!                                           ! DREL= M0/MREL(R)-1
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL (RADIAL PART ONLY)
      REAL(8)    ,INTENT(IN)     :: PHI(NR) !WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: TPHI(NR) !KINETIC ENERGY * WAVE-FUNCTION
      REAL(8)    ,PARAMETER      :: PI=4.D0*ATAN(1.D0)
      REAL(8)                    :: A(NR) 
      REAL(8)                    :: B(NR) 
      REAL(8)                    :: C(NR) 
      REAL(8)                    :: DPHI(NR) 
      REAL(8)                    :: D2PHI(NR) 
      REAL(8)                    :: R(NR) 
      REAL(8)                    :: SOFACTOR
      REAL(8)                    :: RDPRIME(NR)
!     **************************************************************************
!
!     ==========================================================================
!     == SPIN ORBIT COUPLING                                                  ==
!     ==========================================================================
      IF (SO.EQ.0) THEN
        SOFACTOR=0.D0
      ELSE IF(SO.EQ.1) THEN
        SOFACTOR=REAL(L,KIND=8)       ! PARALLEL SPIN AND ORBIT
      ELSE IF(SO.EQ.-1) THEN
        SOFACTOR=REAL(-L-1,KIND=8)    ! ANTIPARALLELSPIN AND ORBIT
      ELSE
         CALL ERROR$MSG('SO CAN ONLY HAVE VALUES -1,0,1')
         CALL ERROR$STOP('RADIAL$SCHRODINGER')
      END IF
!
!     ==========================================================================
!     == SET UP DIFFERENTIAL EQUATION                                         ==
!     ==========================================================================
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
      CALL RADIAL$R(GID,NR,R)
      A(:)=1.D0+DREL(:)
!     == AVOID DIVIDE BY ZERO IF THE FIRST GRID POINT IS THE ORIGIN.
!     == THE FORCES ON THE FIRST GRID POINT ARE NOT USED,
!     == BECAUSE RADIAL$DGL IS BASED ON THE VERLET ALGORITHM
!     == THAT CANNOT USE THE FORCES ON THE FIRST AND LAST GRID POINT
      B(2:)=2.D0*(1.D0+DREL(2:))/R(2:)+RDPRIME(2:)
      C(2:)=-(1.D0+DREL(2:))*REAL(L*(L+1),KIND=8)/R(2:)**2 &
     &    -RDPRIME(2:)*SOFACTOR/R(2:)
      B(1)=B(2)
      C(1)=C(2)
!
!     ==========================================================================
!     == BOUNDARY CONDITIONS                                                  ==
!     ==========================================================================
      CALL RADIAL$VERLETD1(GID,NR,PHI,DPHI)
      CALL RADIAL$VERLETD2(GID,NR,PHI,D2PHI)
      TPHI(:)=A(:)*D2PHI+B(:)*DPHI+C(:)*PHI
      TPHI(:)=-0.5D0*TPHI(:)
      RETURN
      END SUBROUTINE SCHROEDINGER$EKINPSI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER$LBND_SCALREL(GID,NR,LX,LMX,LMRX,TMAINSH,POT,DREL &
     &                          ,G,ENU,NPHI,EB,PHI,TPHI,TOK)
!     **                                                                      **
!     **  APPROXIMATE BOUND STATES FOR THE SCALAR-RELATIVISTIC,               **
!     **  RADIAL DIRAC EQUATION FOR THE LARGE COMPONENT.                      **
!     **                                                                      **
!     **  THE SOLUTIONS ARE CALCULATED TO LINEAR ORDER IN (E-ENU)             **
!     **                                                                      **
!     **  DREL=1/MREL-1/M0 IS A MEASURE FOR THE RELATIVISTIC EFFECTS,         **
!     **  WHERE MREL=M0+(E-V)/2C^2 AND M0 IS THE REST MASS                    **
!     **  THE NONRELATIVISTIC SOLUTION IS OBTAINED BY SETTINH DREL=0          **
!     **                                                                      **
!     **  THE SOLUTIONS ARE SELECTED BY NPHI, FROM WHICH THE DOMINATING       **
!     **  ANGULAR MOMENTUM IS EXTRACTED BY NPHI=2*L+1.                        **
!     **                                                                      **
!     **  THE DIFFERENTIAL EQUATION IS SOLVED WITH THE VERLET ALGORITHM       **
!     **    DPHI/DX=(PHI(+)-PHI(-))/2                                         **
!     **    D2PHI/DX2=PHI(+)-2PHI(0)+PHI(-)                                   **
!     **  WHERE X IS THE VARIABLE WITH R(X=I)=R_I                             **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ******    **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID-ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: LX      ! X(ANGULAR MOMENTUM FOR WAVEF.)
      INTEGER(4) ,INTENT(IN)     :: LMX     ! X#(WAVE FUNCTION ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)     :: LMRX    ! X#(POTENTIAL ANGULAR MOMENTA)
      LOGICAL(4) ,INTENT(IN)     :: TMAINSH(LX+1)
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
!                                           ! DREL= M0/MREL(R)-1
      REAL(8)    ,INTENT(IN)     :: G(NR,LMX)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: ENU     ! EXPANSION ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR,LMRX) !POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: NPHI              ! #(WAVE FUNCTIONS)
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,LMX,NPHI)  ! WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: TPHI(NR,LMX,NPHI) ! P**2/(2M)*WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: EB(NPHI)          ! ONE-PARTICKE ENERGIES
      LOGICAL(4) ,INTENT(OUT)    :: TOK               ! ERROR FLAG
      INTEGER(4)                 :: LM,LM1,LM2,LM3,L,IM,IB
      REAL(8)                    :: A(NR)
      REAL(8)                    :: B(NR)
      REAL(8)                    :: C(NR,LMX,LMX)
      REAL(8)                    :: D(NR,LMX) 
      REAL(8)                    :: R(NR)                        !
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: CPOT(NR,LMX,LMX)
      REAL(8)                    :: CG
      REAL(8)                    :: RDPRIME(NR)                  !
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20
      INTEGER(4)                 :: IRCL,IROUT,ISVAR
      LOGICAL(4)                 :: TCHK
!     **************************************************************************
      TOK=.FALSE.
      IF((LX+1)**2.NE.LMX) THEN   !LMX=(LX+1)**2
        CALL ERROR$MSG('LMX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMX',LMX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('SCHROEDINGER$LBND_SCALREL')
      END IF
      ISVAR=0
      DO L=0,LX
        IF(TMAINSH(L+1)) ISVAR=ISVAR+2*L+1
      ENDDO
      IF(ISVAR.NE.NPHI) THEN
        CALL ERROR$MSG('TMAINSH AND NPHI ARE INCONSISTENT')
        CALL ERROR$STOP('SCHROEDINGER$LBND_SCALREL')
      END IF
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
!
!     ==========================================================================
!     ==  DETERMINE CLASSICAL TURNING POINT R(IRCL)                           ==
!     ==  AND OUTMOST POINT FOR INWARD INTEGRATION R(IROUT)                   ==
!     ==========================================================================
      CALL SCHROEDINGER_SPECIALRADS(GID,NR,0,XMAX,POT(:,1),ENU,IRCL,IROUT)
!
!     ==========================================================================
!     ==  PREPARE POTENTIAL-INDEPENDENT ARRAYS                                ==
!     ==========================================================================
!     == A*D2F/DR2+B*DF/DR+C*F=D
      A(:)=1.D0+DREL(:)
!     == AVOID DIVIDE BY ZERO IF THE FIRST GRID POINT IS THE ORIGIN.
!     == THE FORCES ON THE FIRST GRID POINT ARE NOT USED,
!     == BECAUSE RADIAL$DGL IS BASED ON THE VERLET ALGORITHM
!     == THAT CANNOT USE THE FORCES ON THE FIRST AND LAST GRID POINT
      B(2:)=2.D0*(1.D0+DREL(2:))/R(2:)+RDPRIME(2:)
      B(1)=B(2)
      D(:,:)=-2.D0*G(:,:)
!
!     ==========================================================================
!     ==  COUPLING BETWEEN WAVE FUNCTION COMPONENTS VIA POTENTIAL             ==
!     ==========================================================================
      C(:,:,:)=0.D0
      DO LM1=1,LMX
        DO LM2=1,LMX
          AUX(:)=0.D0
          DO LM3=1,LMRX
            CALL CLEBSCH(LM1,LM2,LM3,CG)
            IF(CG.EQ.0.D0) CYCLE
            AUX=CG*POT(:,LM3)
            IF(LM3.EQ.1) AUX(IROUT+1:)=AUX(IROUT)  ! CONSTANT POTENTIAL BEYOND IROUT
            C(:,LM1,LM2)=C(:,LM1,LM2)+AUX
          ENDDO
        ENDDO
      ENDDO
      CPOT=C   ! KEEP POTENTIAL TO DETERMINE THE KINETIC ENERGY
      C=-2.D0*C
!
!     ==========================================================================
!     ==  ADD KINETIC ENERGY TERM TO C AND SHIFT ENERGY ZERO TO ENU           ==
!     ==========================================================================
      LM=0
      DO L=0,LX
        AUX(1)=0.D0
        AUX(2:)=-(1.D0+DREL(2:))/R(2:)**2 * REAL(L*(L+1),KIND=8)+2.D0*ENU
        DO IM=1,2*L+1
          LM=LM+1
          C(:,LM,LM)=C(:,LM,LM)+AUX(:)
        ENDDO
      ENDDO
!     ==  AVOID DIVIDE-BY-ZERO
      C(1,:,:)=C(2,:,:)
!
!     ==========================================================================
!     ==  DETERMINE BOUND STATES                                              ==
!     ==========================================================================
      CALL SCHROEDINGER_XXXR(GID,NR,LX,LMX,IRCL,IROUT,TMAINSH,A,B,C,D,NPHI,EB,PHI,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$STOP('SCHROEDINGER_XXXR FINISHED WITH ERROR')
        CALL ERROR$STOP('RADIAL$NONSPHBOUND')
      END IF
!
!     ==========================================================================
!     ==  SHIFT ENERGIES                                                      ==
!     ==========================================================================
      EB(:)=ENU+EB(:) ! CHANGE ENERGIES RELATIVE TO ENY TO ABSOLUTE ENERGIES
!
!     ==========================================================================
!     ==  DETERMINE TPHI AND TSPHI                                            ==
!     ==========================================================================
!     -- BETTER DIRECTLY WORK OUT THE KINETIC ENERGY BECAUSE  THE 
!     -- SCHROEDINGER EQUATION IS FULFILLED ONLY TO FIRST ORDER IN DE
      DO IB=1,NPHI
        TPHI(:,:,IB)=EB(IB)*PHI(:,:,IB)
        DO LM1=1,LMX
          DO LM2=1,LMX
            TPHI(:,LM1,IB)=TPHI(:,LM1,IB)-CPOT(:,LM1,LM2)*PHI(:,LM2,IB)
          ENDDO
        ENDDO
      ENDDO
      TOK=.TRUE.
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER_XXXR(GID,NR,LX,NF,IRMATCH,IROUT,TMAINSH,A,B,C,D,NPHI,DE &
     &                            ,PHI,TOK)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID             ! GRID-ID FOR RADIAL GRID
      INTEGER(4),INTENT(IN) :: NR              ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: LX              ! MAX ANGULAR MOMENTUM OF LM EXPANSION
      INTEGER(4),INTENT(IN) :: NF              ! #(ANGULAR MOMENTA)
      INTEGER(4),INTENT(IN) :: IRMATCH         ! MATCHING POINT FOR INSIDE-OUTSIDE INTEGRATION
      INTEGER(4),INTENT(IN) :: IROUT           ! OUTERMOST POINT TO BE CONSIDERED
      LOGICAL(4),INTENT(IN) :: TMAINSH(LX+1)   ! 
      REAL(8)   ,INTENT(IN) :: A(NR)
      REAL(8)   ,INTENT(IN) :: B(NR)
      REAL(8)   ,INTENT(IN) :: C(NR,NF,NF)
      REAL(8)   ,INTENT(IN) :: D(NR,NF)
      INTEGER(4),INTENT(IN) :: NPHI            ! #(WAVE FUNCTIONS)
      REAL(8)   ,INTENT(OUT):: DE(NPHI)        ! ONE-PARTICLE EIGENVALUES
      REAL(8)   ,INTENT(OUT):: PHI(NR,NF,NPHI) ! WAVE FUNCTIONS
      LOGICAL(4),INTENT(OUT):: TOK             ! ERROR FLAG
      REAL(8)               :: ALLPHIL(NR,NF,NF)
      REAL(8)               :: ALLPHIR(NR,NF,NF)
      REAL(8)               :: PHIL(NR,NF,NPHI)
      REAL(8)               :: PHIR(NR,NF,NPHI)
      REAL(8)               :: PHIL_DOT(NR,NF,NPHI)
      REAL(8)               :: PHIR_DOT(NR,NF,NPHI)
      REAL(8)               :: PHIWORK2D(NR,NF,NF)
      INTEGER(4)            :: IF,IF1,IF2
      REAL(8)               :: HA(NF-NPHI,NF-NPHI),HX(NF-NPHI,NPHI),HB(NF-NPHI,NPHI)
      REAL(8)               :: MAT(NF,NF)
      INTEGER(4)            :: IRC
      REAL(8)               :: KINK_HOM(NPHI,NPHI)
      REAL(8)               :: KINK_DOT(NPHI,NPHI)
      REAL(8)               :: KINKC(NPHI,NPHI)
      REAL(8)               :: HAM(NPHI,NPHI),OV(NPHI,NPHI)
      REAL(8)               :: R(NR)
      REAL(8)               :: DHOM(NR,NF)
      LOGICAL   ,PARAMETER  :: TWRITE=.FALSE.
      REAL(8)               :: BVECS(NF,NF),XVECS(NF,NF)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: SVAR
      INTEGER(4)            :: I,J,L,M,LM
      INTEGER(4)            :: LOX(NF)
      LOGICAL(4)            :: TMAIN(NF)  !SELECTS THE MAIN CONTRIBUTIONS
!     **************************************************************************
      TOK=.FALSE.
!PRINT*,'NEW SCHROEDINGER_XXXR STARTED',NPHI,NF
      IF(IROUT+1.GT.NR) THEN
        CALL ERROR$MSG('IROUT OUT OF RANGE')
        CALL ERROR$STOP('SCHROEDINGER_XXXR')
      END IF
      IF(IRMATCH.GT.IROUT) THEN
        CALL ERROR$MSG('IRMATCH OUT OF RANGE')
        CALL ERROR$STOP('SCHROEDINGER_XXXR')
      END IF
      IF(NF.NE.(LX+1)**2) THEN
        CALL ERROR$MSG('NF INCONSISTENT WITH LX')
        CALL ERROR$STOP('SCHROEDINGER_XXXR')
      END IF
!
!     ==========================================================================
!     ==  INITIALIZATIONS                                                     ==
!     ==========================================================================
      LM=0
      DO L=0,LX
        DO M=1,2*L+1
          LM=LM+1
          LOX(LM)=L
        ENDDO
      ENDDO
      TMAIN(:)=.FALSE.
      DO LM=1,NF
        TMAIN(LM)=TMAINSH(LOX(LM)+1)
      ENDDO
      CALL RADIAL$R(GID,NR,R)
      IRC=IRMATCH
!
!     ==========================================================================
!     ==  OBTAIN HOMOGENEOUS SOLUTION                                         ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXR: DETERMINE PHI',LOX
      ALLPHIL(:,:,:)=0.D0
      ALLPHIR(:,:,:)=0.D0
      DO IF=1,NF
        ALLPHIL(1:2,IF,IF)=R(1:2)**LOX(IF)
        DHOM(:,:)=0.D0
        CALL RADIAL$DGLGEN(GID,NR,NF,1,IRC+1,A,B,C,DHOM,ALLPHIL(:,:,IF))
        SVAR=MAXVAL(ABS(ALLPHIL(IRC,:,IF)))
        ALLPHIL(1:IRC+1,:,IF)=ALLPHIL(1:IRC+1,:,IF)/SVAR
        DHOM(IROUT,IF)=1.D-8
        CALL RADIAL$DGLGEN(GID,NR,NF,IROUT+1,IRC-1,A,B,C,DHOM,ALLPHIR(:,:,IF))
        SVAR=MAXVAL(ABS(ALLPHIR(IRC,:,IF)))
        ALLPHIR(IRC-1:IROUT+1,:,IF)=ALLPHIR(IRC-1:IROUT+1,:,IF)/SVAR
      ENDDO
!
!     ==========================================================================
!     ==  MAKE PHI_HOM SOLUTIONS CONTINUOUS                                   ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXR: MATCH PHI'
      MAT(:,:)=ALLPHIR(IRC,:,:)
!     == MATRIX IS NOT SYMMETRIC. THUS SOLVE WITH SINGULAR VALUE DECOMPOSITION
      BVECS(:,1:NF)=ALLPHIL(IRC,:,:)
      CALL LIB$MATRIXSOLVER8(NF,NF,NF,MAT,XVECS,BVECS)
      PHIWORK2D(:,:,:)=ALLPHIR(:,:,:)
      ALLPHIR(:,:,:)=0.D0
      DO IF1=1,NF
        DO IF2=1,NF
          ALLPHIR(:,:,IF1)=ALLPHIR(:,:,IF1)+PHIWORK2D(:,:,IF2)*XVECS(IF2,IF1)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == REMOVE KINKS IN NON-RELEVANT ANGULAR MOMENTA CHANNELS                ==
!     ==========================================================================
!PRINT*,'REMOVE KINKS OF OTHER ANGULAR MOMENTUM CHANNELS OF PHI'
      I=0
      DO IF1=1,NF
        IF(TMAIN(IF1)) CYCLE
        I=I+1
        J=0
        DO IF2=1,NF
          IF(TMAIN(IF2)) CYCLE
          J=J+1
          SVAR=(ALLPHIR(IRC+1,IF1,IF2)-ALLPHIR(IRC-1,IF1,IF2)) &
     &         -(ALLPHIL(IRC+1,IF1,IF2)-ALLPHIL(IRC-1,IF1,IF2))
          HA(I,J)=SVAR
        ENDDO
        J=0
        DO IF2=1,NF
          IF(.NOT.TMAIN(IF2)) CYCLE
          J=J+1
          SVAR=(ALLPHIR(IRC+1,IF1,IF2)-ALLPHIR(IRC-1,IF1,IF2)) &
     &         -(ALLPHIL(IRC+1,IF1,IF2)-ALLPHIL(IRC-1,IF1,IF2))
          HB(I,J)=-SVAR
        ENDDO
      ENDDO
      IF(NF.GT.NPHI) THEN
        CALL LIB$MATRIXSOLVER8(NF-NPHI,NF-NPHI,NPHI,HA,HX,HB)
      END IF
      PHIL(:,:,:)=0.D0
      PHIR(:,:,:)=0.D0
      I=0
      DO IF1=1,NF
        IF(.NOT.TMAIN(IF1)) CYCLE
        I=I+1
        PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF1)
        PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I) &
     &                         +ALLPHIR(IRC-1:IROUT+1,:,IF1)
        J=0
        DO IF2=1,NF
          IF(TMAIN(IF2)) CYCLE
          J=J+1
          PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF2)*HX(J,I)
          PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I) &
     &                            +ALLPHIR(IRC-1:IROUT+1,:,IF2)*HX(J,I)
        ENDDO
      ENDDO  
!
!     ==========================================================================
!     ==   ORTHOGONALIZE PHI                                                  ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXC: ORTHOGONALIZE PHI'
      DO I=1,NPHI
!       == ORTHOGONALIZE
        DO J=1,I-1
          AUX(:)=0.D0
          DO IF=1,NF
            AUX(:IRC)=AUX(:IRC)+PHIL(:IRC,IF,J)*PHIL(:IRC,IF,I)
            AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
     &                   +PHIR(IRC+1:IROUT+1,IF,J)*PHIR(IRC+1:IROUT+1,IF,I)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)-PHIL(:IRC+1,:,J)*SVAR
          PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I)-PHIR(IRC-1:IROUT+1,:,J)*SVAR
        ENDDO
!       ==  NORMALIZE
        AUX(:)=0.D0
        DO IF=1,NF
          AUX(:IRC)         =AUX(:IRC)         +PHIL(:IRC,IF,I)**2
          AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1)+PHIR(IRC+1:IROUT+1,IF,I)**2
        ENDDO
        AUX(:)=AUX(:)*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        SVAR=1.D0/SQRT(SVAR)
        PHIL(:IRC+1,:,I)       =PHIL(:IRC+1,:,I)*SVAR
        PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I)*SVAR
      ENDDO
!
!     ==========================================================================
!     ==  DETERMINE PHIDOT                                                    ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXC: DETERMINE PHIDOT'
      PHIL_DOT(:,:,:)=0.D0
      PHIR_DOT(:,:,:)=0.D0
      DO IF=1,NPHI
        CALL RADIAL$DGLGEN(GID,NR,NF,1,IRC+1,A,B,C,-2.D0*PHIL(:,:,IF) &
     &                     ,PHIL_DOT(:,:,IF))
        CALL RADIAL$DGLGEN(GID,NR,NF,IROUT+1,IRC-1,A,B,C,-2.D0*PHIR(:,:,IF) &
     &                     ,PHIR_DOT(:,:,IF))
      ENDDO
!
!     ==========================================================================
!     ==  MAKE PHI_DOT CONTINUOUS                                             ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXC: MATCH PHIDOT'
      MAT(:,:)=ALLPHIR(IRC,:,:)
      BVECS(:,:NPHI)=PHIL_DOT(IRC,:,:)-PHIR_DOT(IRC,:,:)
      CALL LIB$MATRIXSOLVER8(NF,NF,NPHI,MAT,XVECS(:,:NPHI),BVECS(:,:NPHI))
!
      DO IF1=1,NPHI
        DO IF2=1,NF
          PHIR_DOT(IRC-1:IROUT+1,:,IF1)=PHIR_DOT(IRC-1:IROUT+1,:,IF1) &
     &                            +ALLPHIR(IRC-1:IROUT+1,:,IF2)*XVECS(IF2,IF1)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == REMOVE KINKS IN NON-RELEVANT ANGULAR MOMENTUM CHANNELS               ==
!     ==========================================================================
!PRINT*,'REMOVE KINKS OF OTHER ANGULAR MOMENTUM CHANNELS OF PHIDOT'
      I=0
      DO IF1=1,NF
        IF(TMAIN(IF1)) CYCLE
        I=I+1
        DO J=1,NPHI
          SVAR=(PHIR_DOT(IRC+1,IF1,J)-PHIR_DOT(IRC-1,IF1,J)) &
     &         -(PHIL_DOT(IRC+1,IF1,J)-PHIL_DOT(IRC-1,IF1,J))
          HB(I,J)=-SVAR
        ENDDO
      ENDDO
      IF(NF.GT.NPHI) THEN
        CALL LIB$MATRIXSOLVER8(NF-NPHI,NF-NPHI,NPHI,HA,HX,HB)
      END IF
      DO I=1,NPHI
        J=0
        DO IF2=1,NF
          IF(TMAIN(IF2)) CYCLE
          J=J+1
          PHIL_DOT(:IRC+1,:,I)=PHIL_DOT(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF2)*HX(J,I)
          PHIR_DOT(IRC-1:IROUT+1,:,I)=PHIR_DOT(IRC-1:IROUT+1,:,I) &
     &                               +ALLPHIR(IRC-1:IROUT+1,:,IF2)*HX(J,I)
        ENDDO
      ENDDO  
!
!     ==========================================================================
!     ==   ORTHOGONALIZE PHIDOT TO PHI                                        ==
!     ==========================================================================
      DO I=1,NPHI
        DO J=1,NPHI
          AUX(:)=0.D0
          DO IF=1,NF
            AUX(:IRC)=AUX(:IRC)+PHIL(:IRC,IF,J)*PHIL_DOT(:IRC,IF,I)
            AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
     &             +PHIR(IRC+1:IROUT+1,IF,J)*PHIR_DOT(IRC+1:IROUT+1,IF,I)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          PHIL_DOT(:IRC+1,:,I)=PHIL_DOT(:IRC+1,:,I)-PHIL(:IRC+1,:,J)*SVAR
          PHIR_DOT(IRC-1:IROUT+1,:,I)=PHIR_DOT(IRC-1:IROUT+1,:,I)-PHIR(IRC-1:IROUT+1,:,J)*SVAR
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  REMOVE KINKS BY MIXING PHIDOT INTO PHI                              ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXC:MATCH KINKS'
      I=0
      DO IF2=1,NF
        IF(.NOT.TMAIN(IF2)) CYCLE
        I=I+1
        KINK_HOM(I,:)=(PHIR(IRC+1,IF2,:)-PHIR(IRC-1,IF2,:)) &
     &               -(PHIL(IRC+1,IF2,:)-PHIL(IRC-1,IF2,:))
        KINK_DOT(I,:)=(PHIR_DOT(IRC+1,IF2,:)-PHIR_DOT(IRC-1,IF2,:)) &
     &               -(PHIL_DOT(IRC+1,IF2,:)-PHIL_DOT(IRC-1,IF2,:))
      ENDDO
      CALL LIB$MATRIXSOLVER8(NPHI,NPHI,NPHI,-KINK_DOT,HAM,KINK_HOM)
!     == MAKE PHI DIFFERENTIABLE =========================================
      DO I=1,NPHI
        DO J=1,NPHI
          PHIL(:IRC,:,I)         =PHIL(:IRC,:,I)         +PHIL_DOT(:IRC,:,J)  *HAM(J,I)
          PHIR(IRC+1:IROUT+1,:,I)=PHIR(IRC+1:IROUT+1,:,I)+PHIR_DOT(IRC+1:IROUT+1,:,J)*HAM(J,I)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  DETERMINE OVERLAP MATRIX                                    ==
!     ==================================================================
      DO I=1,NPHI
        DO J=I,NPHI
          AUX(:)=0.D0
          DO IF=1,NF
            AUX(:IRC)         =AUX(:IRC)+PHIL(:IRC,IF,I)*PHIL(:IRC,IF,J)
            AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
     &                         +PHIR(IRC+1:IROUT+1,IF,I)*PHIR(IRC+1:IROUT+1,IF,J)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          OV(I,J)=SVAR
          OV(J,I)=OV(I,J)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  DETERMINE EIGENSTATES                                       ==
!     ==================================================================
      HAM=0.5D0*(HAM+TRANSPOSE(HAM))
      CALL LIB$GENERALEIGENVALUER8(NPHI,HAM,OV,DE,KINKC)
      PHI(:,:,:)=0.D0
      DO I=1,NPHI
        DO J=1,NPHI
          PHI(:IRC,:,I)         =PHI(:IRC,:,I)         +PHIL(:IRC,:,J)*KINKC(J,I)      
          PHI(IRC+1:IROUT+1,:,I)=PHI(IRC+1:IROUT+1,:,I)+PHIR(IRC+1:IROUT+1,:,J)*KINKC(J,I)      
        ENDDO
      ENDDO
      TOK=.TRUE.

!CALL SCHROEDINGER_WRITEPHI(GID,NR,'FINAL',NF,NPHI,IRC,PHI,PHI)
      RETURN
      END SUBROUTINE SCHROEDINGER_XXXR
!
!     ..................................................................
!SCHROEDINGER$NONSPHERICAL_BOUND_C
      SUBROUTINE SCHROEDINGER$LBND_FULLYREL(GID,NR,NDIMD,LMX,LMRX,POT,DREL,G,ENU &
     &                             ,NPHI,EB,PHI,TPHI,SPHI,TSPHI,TOK)
!     **                                                                  **
!     **  SOLVES THE RELATIVISTIC RADIAL DIRAC EQUATION FOR THE           **
!     **  LARGE COMPONENT. DREL=1/MREL-1/M0 IS A MEASURE FOR THE          **
!     **  RELATIVISTIC EFFECTS, WHERE MREL=M0+(E-V)/2C^2 AND M0 IS THE    **
!     **  REST MASS. V=POT*Y0 IS THE POTENTIAL.                           **
!     **  SPIN ORBIT COUPLING IS MEASURED BY SO WHICH CAN HAVE THE VALUES:**
!     **    SO=0 NO-SPIN ORBIT COUPLING                                   **
!     **    SO=1 PARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM               **
!     **  G IS AN INHOMOGENEITY, WHICH MUST BE SET TO ZERO FOR THE        **
!     **  HOMOGENEOUS SOLUTION.                                           **
!     **                                                                  **
!     **  THE SOLUTIONS ARE CALCULATED TO LINEAR ORDER IN DE              **
!     **  WHERE E+DE IS THE NEW ENERGY OF THE WAVE FUNCTIONS              **
!     **                                                                  **
!     **  THE DIFFERENTIAL EQUATION IS SOLVED WITH THE VERLET ALGORITHM   **
!     **    DPHI/DX=(PHI(+)-PHI(-))/2                                     **
!     **    D2PHI/DX2=PHI(+)-2PHI(0)+PHI(-)                               **
!     **  WHERE X IS THE VARIABLE WITH R(X=I)=R_I                         **
!     **                                                                  **
!     **  THE NONRELATIVISTIC SOLUTION IS OBTAINED BY SETTINH DREL=0      **
!     **                                                                  **
!     **  ATTENTION! THE ROUTINE IS NOT GUARDED AGAINST OVERFLOW DUE      **
!     **    TO THE EXPONENTIAL INCREASE OF THE SOLUTION                   **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: NDIMD   ! (1,2,4)#(SPINOR COMPONENTS)
      INTEGER(4) ,INTENT(IN)     :: LMX     ! X#(WAVE FUNCTION ANGULAR MOMENTA
      INTEGER(4) ,INTENT(IN)     :: LMRX    ! X#(POTENTIAL ANGULAR MOMENTA)
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
!                                           ! DREL= M0/MREL(R)-1
      REAL(8)    ,INTENT(IN)     :: G(NR,LMX,2)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: ENU       !ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR,LMRX,NDIMD) !POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: NPHI
      COMPLEX(8) ,INTENT(OUT)    :: PHI(NR,LMX,2,NPHI) ! WAVE-FUNCTION
      COMPLEX(8) ,INTENT(OUT)    :: TPHI(NR,LMX,2,NPHI) ! WAVE-FUNCTION
      COMPLEX(8) ,INTENT(OUT)    :: SPHI(NR,LMX,2,NPHI) ! WAVE-FUNCTION
      COMPLEX(8) ,INTENT(OUT)    :: TSPHI(NR,LMX,2,NPHI) ! KINETIC ENERGY * WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: EB(NPHI)
      LOGICAL(4) ,INTENT(OUT)    :: TOK
      REAL(8)                    :: DE(NPHI)
      INTEGER(4)                 :: LX
      INTEGER(4)                 :: LOX(LMX,2) ! ANGULAR MOMENTA
      INTEGER(4)                 :: LM,LM1,LM2,LM3,L,M,IM,IS,IB,IS1,IS2
      REAL(8)                    :: A(NR)
      REAL(8)                    :: B(NR)
      COMPLEX(8)                 :: C(NR,LMX,2,LMX,2)
      COMPLEX(8)                 :: CPOT(NR,LMX,2,LMX,2)
      COMPLEX(8)                 :: CR(LMX,2,LMX,2)
      COMPLEX(8)                 :: CA(LMX,2,LMX,2)
      COMPLEX(8)                 :: D(NR,LMX,2) 
      REAL(8)                    :: PHIR(NR,LMX,2,2*LMX)
      REAL(8)                    :: DPHI(NR,LMX,2)
      REAL(8)                    :: R(NR) 
      REAL(8)                    :: AUX(NR),AUX1(NR),AUX2(NR) 
      COMPLEX(8)                 :: CSVAR
      REAL(8)                    :: CG
      REAL(8)                    :: RDPRIME(NR)
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20
      INTEGER(4)                 :: IRCL,IROUT
      COMPLEX(8),PARAMETER       :: CI=(0.D0,1.D0)
      COMPLEX(8)                 :: CLS(LMX,2,LMX,2) ! SPIN ORBIT MATRIX (LS)
      LOGICAL(4)                 :: TCHK
!     ************************************************************************
      TOK=.FALSE.
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
      LX=INT(SQRT(REAL(LMX))-1.D0)
      IF((LX+1)**2.NE.LMX) THEN   !LMX=(LX+1)**2
        CALL ERROR$MSG('LMX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMX',LMX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('RADIAL$NONSPHBOUND')
      END IF
      LM=0
      DO L=0,LX
        DO M=1,2*L+1
          LM=LM+1
          LOX(LM,:)=L
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  DETERMINE CLASSICAL TURNING POINT R(IRCL)                   ==
!     ==  AND OUTMOST POINT FOR INWARD INTEGRATION R(IROUT)           ==
!     ==================================================================
      CALL SCHROEDINGER_SPECIALRADS(GID,NR,0,XMAX,POT(:,1,1),ENU,IRCL,IROUT)
!
!     ==================================================================
!     ==  PREPARE POTENTIAL-INDEPENDENT ARRAYS                        ==
!     ==================================================================
!     == A*D2F/DR2+B*DF/DR+C*F=D
      A(:)=1.D0+DREL(:)
!     == AVOID DIVIDE BY ZERO IF THE FIRST GRID POINT IS THE ORIGIN.
!     == THE FORCES ON THE FIRST GRID POINT ARE NOT USED,
!     == BECAUSE RADIAL$DGL IS BASED ON THE VERLET ALGORITHM
!     == THAT CANNOT USE THE FORCES ON THE FIRST AND LAST GRID POINT
      B(2:)=2.D0*(1.D0+DREL(2:))/R(2:)+RDPRIME(2:)
      B(1)=B(2)
      D(:,:,:)=-2.D0*G(:,:,:)
!
!     ==================================================================
!     ==  COUPLING BETWEEN WAVE FUNCTION COMPONENTS VIA POTENTIAL     ==
!     ==================================================================
      C(:,:,:,:,:)=(0.D0,0.D0)
      DO LM1=1,LMX
        DO LM2=1,LMX
          AUX(:)=0.D0
          DO LM3=1,LMRX
            CALL CLEBSCH(LM1,LM2,LM3,CG)
            IF(CG.EQ.0.D0) CYCLE
            AUX=CG*POT(:,LM3,1)
            IF(LM3.EQ.1) AUX(IROUT+1:)=AUX(IROUT)  ! CONSTANT POTENTIAL BEYOND IROUT
            C(:,LM1,1,LM2,1)=C(:,LM1,1,LM2,1)+AUX
            C(:,LM1,2,LM2,2)=C(:,LM1,2,LM2,2)+AUX
            IF(NDIMD.GT.1) THEN
              IF(NDIMD.EQ.2) THEN
                AUX=CG*POT(:,LM3,2)
                C(:,LM1,1,LM2,1)=C(:,LM1,1,LM2,1)+AUX
                C(:,LM1,2,LM2,2)=C(:,LM1,2,LM2,2)-AUX
              ELSE
                AUX=CG*POT(:,LM3,2)
                C(:,LM1,1,LM2,2)=C(:,LM1,1,LM2,2)+AUX
                C(:,LM1,2,LM2,1)=C(:,LM1,2,LM2,1)-AUX
                AUX=CG*POT(:,LM3,3)
                C(:,LM1,1,LM2,2)=C(:,LM1,1,LM2,2)-CI*AUX
                C(:,LM1,2,LM2,1)=C(:,LM1,2,LM2,1)+CI*AUX
                AUX=CG*POT(:,LM3,4)
                C(:,LM1,1,LM2,1)=C(:,LM1,1,LM2,1)+AUX
                C(:,LM1,2,LM2,2)=C(:,LM1,2,LM2,2)-AUX
              END IF
            END IF
          ENDDO
        ENDDO
      ENDDO
      CPOT=C    ! STORE TO EVALUATE KINETIC ENERGY
      C=-2.D0*C
!
!     ==================================================================
!     ==  KINETIC ENERGY TERM TO C                                    ==
!     ==================================================================
      LM=0
      DO L=0,LX
        AUX(1)=0.D0
        AUX(2:)=-(1.D0+DREL(2:))/R(2:)**2 * REAL(L*(L+1),KIND=8)+2.D0*ENU
        DO IM=1,2*L+1
          LM=LM+1
          DO IS=1,2
            C(:,LM,IS,LM,IS)=C(:,LM,IS,LM,IS)+AUX(:)
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  ADD SPIN ORBIT COUPLING TO E                                ==
!     ==================================================================
      CALL SCHROEDINGER_LS(LMX,CLS)
!CLS=0.D0
      AUX(2:)=-RDPRIME(2:)/R(2:)
      AUX(1)=AUX(2)
      DO LM1=1,LMX
        DO IS1=1,2
          DO LM2=1,LMX
            DO IS2=1,2
              IF(ABS(CLS(LM1,IS1,LM2,IS2)).LT.1.D-10) CYCLE
              C(:,LM1,IS1,LM2,IS2)=C(:,LM1,IS1,LM2,IS2)+AUX(:)*CLS(LM1,IS1,LM2,IS2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     ==  AVOID DIVIDE ZEROBYZERO
      C(1,:,:,:,:)=C(2,:,:,:,:)
!
!     ==================================================================
!     ==  DETERMINE BOUND STATES                                      ==
!     ==================================================================
      CALL SCHROEDINGER_XXXC(GID,NR,2*LMX,IRCL,IROUT,LOX,A,B,C,D,NPHI,DE,PHI,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$STOP('SCHROEDINGER_XXXC FINISHED WITH ERROR')
        CALL ERROR$STOP('RADIAL$NONSPHBOUND')
      END IF
!
!     ==================================================================
!     ==  DETERMINE SMALL COMPONENT                                   ==
!     ==================================================================
       SPHI=(0.D0,0.D0)
!      CALL SCHROEDINGER_SMALLCOMPONENT(GID,NR,LMX,NPHI,DREL,PHI,SPHI)
!
!!$      CALL SCHROEDINGER_SP(LMX,CR,CA)
!!$      CALL CONSTANTS$GET('C',SVAR)
!!$      A(:)=(1.D0+DREL(:))/(2.D0*SVAR)
!!$      SPHI(:,:,:,:)=(0.D0,0.D0)
!!$      DO IB=1,NPHI
!!$        DO IS=1,2
!!$          DO LM=1,LMX
!!$            CALL RADIAL$DERIVE(GID,NR,REAL(PHI(:,LM,IS,IB)),AUX1)
!!$            CALL RADIAL$DERIVE(GID,NR,AIMAG(PHI(:,LM,IS,IB)),AUX2)
!!$            DPHI(:,LM,IS)=CMPLX(AUX1,AUX2,KIND=8)
!!$          ENDDO
!!$        ENDDO
!!$!
!!$        DO IS2=1,2
!!$          DO LM2=1,LMX
!!$            DO IS1=1,2
!!$              DO LM1=1,LMX
!!$                CSVAR=CR(LM1,IS1,LM2,IS2)
!!$                IF(CSVAR.NE.0.D0) THEN
!!$                  SPHI(:,LM1,IS1,IB)=SPHI(:,LM1,IS1,IB)+CSVAR*DPHI(:,LM2,IS2)
!!$                END IF
!!$                CSVAR=CA(LM1,IS1,LM2,IS2)
!!$                IF(CSVAR.NE.0.D0) THEN
!!$                  SPHI(:,LM1,IS1,IB)=SPHI(:,LM1,IS1,IB)+CSVAR*PHI(:,LM2,IS2,IB)
!!$                END IF
!!$              ENDDO
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$!
!!$        DO LM=1,LMX
!!$          DO IS=1,2
!!$            SPHI(:,LM,IS,IB)=A(:)*SPHI(:,LM,IS,IB)                  
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!
!     ==================================================================
!     ==  SHIFT ENERGIES                                              ==
!     ==================================================================
      EB(:)=ENU+DE(:)
!
!     ==================================================================
!     ==  DETERMINE TPHI AND TSPHI                                              ==
!     ==================================================================
!     -- BETTER DIRECTLY WORK OUT THE KINETIC ENERGY BECAUSE  THE 
!     -- SCHROEDINGER EQUATION IS FULFILLED ONLY TO FIRST ORDER IN DE
      DO IB=1,NPHI
        TPHI(:,:,:,IB)=EB(IB)*PHI(:,:,:,IB)
        TSPHI(:,:,:,IB)=EB(IB)*SPHI(:,:,:,IB)
        DO IS1=1,2
          DO LM1=1,LMX
            DO IS2=1,2
              DO LM2=1,LMX
                TPHI(:,LM1,IS1,IB)=TPHI(:,LM1,IS1,IB) &
     &                         -CPOT(:,LM1,IS1,LM2,IS2)*PHI(:,LM2,IS2,IB)
                TSPHI(:,LM1,IS1,IB)=TSPHI(:,LM1,IS1,IB) &
     &                         -CPOT(:,LM1,IS1,LM2,IS2)*SPHI(:,LM2,IS2,IB)
              ENDDO
            ENDDO  
          ENDDO
        ENDDO
      ENDDO
!
      TOK=.TRUE.
      RETURN
      END SUBROUTINE SCHROEDINGER$LBND_FULLYREL
!
!     ..................................................................
      SUBROUTINE SCHROEDINGER_XXXC(GID,NR,NF,IRMATCH,IROUT,LOX,A,B,C,D,NPHI &
     &                            ,DE,PHI,TOK)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NF
      INTEGER(4),INTENT(IN) :: IRMATCH
      INTEGER(4),INTENT(IN) :: IROUT
      INTEGER(4),INTENT(IN) :: LOX(NF)
      REAL(8)   ,INTENT(IN) :: A(NR)
      REAL(8)   ,INTENT(IN) :: B(NR)
      COMPLEX(8),INTENT(IN) :: C(NR,NF,NF)
      COMPLEX(8),INTENT(IN) :: D(NR,NF)
      INTEGER(4),INTENT(IN) :: NPHI
      REAL(8)   ,INTENT(OUT):: DE(NPHI)
      COMPLEX(8),INTENT(OUT):: PHI(NR,NF,NPHI)
      LOGICAL(4),INTENT(OUT):: TOK
      COMPLEX(8)            :: ALLPHIL(NR,NF,NF)
      COMPLEX(8)            :: ALLPHIR(NR,NF,NF)
      COMPLEX(8)            :: PHIL(NR,NF,NPHI)
      COMPLEX(8)            :: PHIR(NR,NF,NPHI)
      COMPLEX(8)            :: PHIL_DOT(NR,NF,NPHI)
      COMPLEX(8)            :: PHIR_DOT(NR,NF,NPHI)
      COMPLEX(8)            :: PHIWORK2D(NR,NF,NF)
      INTEGER(4)            :: IF,IF1,IF2
      COMPLEX(8)            :: HA(NF-NPHI,NF-NPHI),HX(NF-NPHI,NPHI),HB(NF-NPHI,NPHI)
      COMPLEX(8)            :: MAT(NF,NF)
      INTEGER(4)            :: IRC
      COMPLEX(8)            :: KINK_HOM(NPHI,NPHI)
      COMPLEX(8)            :: KINK_DOT(NPHI,NPHI)
      COMPLEX(8)            :: KINKC(NPHI,NPHI)
      COMPLEX(8)            :: HAM(NPHI,NPHI),OV(NPHI,NPHI)
      REAL(8)               :: R(NR)
      COMPLEX(8)            :: DHOM(NR,NF)
      LOGICAL   ,PARAMETER  :: TWRITE=.FALSE.
      COMPLEX(8)            :: BVECS(NF,NF),XVECS(NF,NF)
      COMPLEX(8)            :: CAUX(NR)
      COMPLEX(8)            :: CSVAR
      REAL(8)               :: AUX(NR)
      REAL(8)               :: SVAR,SVAR1,SVAR2
      INTEGER(4)            :: I,J
      INTEGER(4)            :: L0
CHARACTER(32):: FILE
!     ******************************************************************
      TOK=.FALSE.
!PRINT*,'NEW SCHROEDINGER_XXXC STARTED',NPHI,NF
      IF(IROUT+1.GT.NR) THEN
        CALL ERROR$MSG('IROUT OUT OF RANGE')
        CALL ERROR$STOP('SCHROEDINGER_XXXC')
      END IF
      IF(IRMATCH.GT.IROUT) THEN
        CALL ERROR$MSG('IRMATCH OUT OF RANGE')
        CALL ERROR$STOP('SCHROEDINGER_XXXC')
      END IF
      L0=(NPHI/2-1)/2
      IF(NPHI.NE.2*(2*L0+1)) THEN
        CALL ERROR$MSG('NPHI MUST BE (2*L0+1)*2')
        CALL ERROR$STOP('SCHROEDINGER_XXXC')
      END IF
      CALL RADIAL$R(GID,NR,R)
      IRC=IRMATCH
!PRINT*,'IROUT ',IROUT,R(IROUT),IRC,R(IRC)
!
!     ==================================================================
!     ==  OBTAIN HOMOGENEOUS SOLUTION                                 ==
!     ==================================================================
!PRINT*,'SCHROEDINGER_XXXC: DETERMINE PHI',L0
      ALLPHIL(:,:,:)=(0.D0,0.D0)
      ALLPHIR(:,:,:)=(0.D0,0.D0)
      DO IF=1,NF
        ALLPHIL(1:2,IF,IF)=CMPLX(R(1:2)**LOX(IF),0.D0,KIND=8)
        DHOM(:,:)=CMPLX(0.D0,0.D0,KIND=8)
        CALL RADIAL$DGLGENC(GID,NR,NF,1,IRC+1,A,B,C,DHOM,ALLPHIL(:,:,IF))
        SVAR=MAXVAL(ABS(ALLPHIL(IRC,:,IF)))
        ALLPHIL(1:IRC+1,:,IF)=ALLPHIL(1:IRC+1,:,IF)/SVAR
!
        DHOM(IROUT,IF)=CMPLX(1.D-8,0.D0,KIND=8)
        CALL RADIAL$DGLGENC(GID,NR,NF,IROUT+1,IRC-1,A,B,C,DHOM,ALLPHIR(:,:,IF))
        SVAR=MAXVAL(ABS(ALLPHIR(IRC,:,IF)))
        ALLPHIR(IRC-1:IROUT+1,:,IF)=ALLPHIR(IRC-1:IROUT+1,:,IF)/SVAR
      ENDDO
!
!     ==================================================================
!     ==  MAKE PHI_HOM SOLUTIONS CONTINUOUS                           ==
!     ==================================================================
!PRINT*,'SCHROEDINGER_XXXC: MATCH PHI'
      MAT(:,:)=ALLPHIR(IRC,:,:)
!     == MATRIX IS NOT SYMMETRIC. THUS SOLVE WITH SINGULAR VALUE DECOMPOSITION
      BVECS(:,1:NF)=ALLPHIL(IRC,:,:)
      CALL LIB$MATRIXSOLVEC8(NF,NF,NF,MAT,XVECS,BVECS)
      PHIWORK2D(:,:,:)=ALLPHIR(:,:,:)
      ALLPHIR(:,:,:)=(0.D0,0.D0)
      DO IF1=1,NF
        DO IF2=1,NF
          ALLPHIR(:,:,IF1)=ALLPHIR(:,:,IF1)+PHIWORK2D(:,:,IF2)*XVECS(IF2,IF1)
        ENDDO
      ENDDO
!
!     ==================================================================
!     == REMOVE KINKS IN NON-RELEVANT ANGULAR MOMENTA CHANNELS        ==
!     ==================================================================
!PRINT*,'REMOVE KINKS OF OTHER ANGULAR MOMENTUM CHANNELS OF PHI'
      I=0
      DO IF1=1,NF
        IF(LOX(IF1).EQ.L0) CYCLE
        I=I+1
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).EQ.L0) CYCLE
          J=J+1
          CSVAR=(ALLPHIR(IRC+1,IF1,IF2)-ALLPHIR(IRC-1,IF1,IF2)) &
     &         -(ALLPHIL(IRC+1,IF1,IF2)-ALLPHIL(IRC-1,IF1,IF2))
          HA(I,J)=CSVAR
        ENDDO
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).NE.L0) CYCLE
          J=J+1
          CSVAR=(ALLPHIR(IRC+1,IF1,IF2)-ALLPHIR(IRC-1,IF1,IF2)) &
     &         -(ALLPHIL(IRC+1,IF1,IF2)-ALLPHIL(IRC-1,IF1,IF2))
          HB(I,J)=-CSVAR
        ENDDO
      ENDDO
      IF(NF.GT.NPHI) THEN
        CALL LIB$MATRIXSOLVEC8(NF-NPHI,NF-NPHI,NPHI,HA,HX,HB)
      END IF
      PHIL(:,:,:)=(0.D0,0.D0)
      PHIR(:,:,:)=(0.D0,0.D0)
      I=0
      DO IF1=1,NF
        IF(LOX(IF1).NE.L0) CYCLE
        I=I+1
        PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF1)
        PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I) &
     &                         +ALLPHIR(IRC-1:IROUT+1,:,IF1)
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).EQ.L0) CYCLE
          J=J+1
          PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF2)*HX(J,I)
          PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I) &
     &                            +ALLPHIR(IRC-1:IROUT+1,:,IF2)*HX(J,I)
        ENDDO
      ENDDO  
!
!     ====================================================================
!     ==   ORTHOGONALIZE PHI                                                ==
!     ====================================================================
!PRINT*,'SCHROEDINGER_XXXC: ORTHOGONALIZE PHI'
      DO I=1,NPHI
!       == ORTHOGONALIZE
        DO J=1,I-1
          CAUX(:)=(0.D0,0.D0)
          DO IF=1,NF
            CAUX(:IRC)=CAUX(:IRC)+CONJG(PHIL(:IRC,IF,J))*PHIL(:IRC,IF,I)
            CAUX(IRC+1:IROUT+1)=CAUX(IRC+1:IROUT+1) &
     &             +CONJG(PHIR(IRC+1:IROUT+1,IF,J))*PHIR(IRC+1:IROUT+1,IF,I)
          ENDDO
          CAUX(:)=CAUX(:)*R(:)**2
          AUX(:)=REAL(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          AUX(:)=AIMAG(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
          CSVAR=CMPLX(SVAR1,SVAR2,KIND=8)
          PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)-PHIL(:IRC+1,:,J)*CSVAR
          PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I)-PHIR(IRC-1:IROUT+1,:,J)*CSVAR
        ENDDO
!       ==  NORMALIZE
        AUX(:)=0.D0
        DO IF=1,NF
          AUX(:IRC)=AUX(:IRC)+ABS(PHIL(:IRC,IF,I))**2
          AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
     &                      +ABS(PHIR(IRC+1:IROUT+1,IF,I))**2
        ENDDO
        AUX(:)=AUX(:)*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
        SVAR1=1.D0/SQRT(SVAR1)
        PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)*SVAR1
        PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I)*SVAR1
      ENDDO
!
!     ====================================================================
!     ==   TEST CONTINUITY                                              ==
!     ====================================================================
!!$PRINT*,'SCHROEDINGER_XXXC: TEST CONTINUITY'
!!$      SVAR2=0.D0
!!$      DO IF1=1,NPHI
!!$        DO IF2=1,NF
!!$          SVAR1=ABS(PHIL(IRC,IF2,IF1)-PHIR(IRC,IF2,IF1))
!!$          SVAR2=MAX(SVAR2,SVAR1)
!!$        ENDDO
!!$      ENDDO
!!$PRINT*,'DEVIATION IN VALUE',SVAR2
!!$      IF(SVAR2.GT.1.D-6) THEN
!!$        CALL ERROR$MSG('MATCHING OF PHI FAILED')
!!$        CALL ERROR$STOP('SCHROEDINGER_XXXC')
!!$      END IF   
!!$!
!!$!     == TEST DIFFERENTIABILITY ===================================
!!$      SVAR2=0.D0
!!$      DO IF1=1,NPHI
!!$        DO IF2=1,NF
!!$          IF(LOX(IF2).EQ.L0) CYCLE
!!$          SVAR1=ABS( (PHIR(IRC+1,IF2,IF1)-PHIR(IRC-1,IF2,IF1)) &
!!$     &              -(PHIL(IRC+1,IF2,IF1)-PHIL(IRC-1,IF2,IF1)))
!!$          SVAR2=MAX(SVAR2,SVAR1)
!!$        ENDDO
!!$      ENDDO
!!$PRINT*,'DEVIATION IN OTHER DERIVATIVES',SVAR2
!!$      IF(SVAR2.GT.1.D-6) THEN
!!$        CALL ERROR$MSG('MATCHING OF DPHI/DR FAILED')
!!$        CALL ERROR$STOP('SCHROEDINGER_XXXC')
!!$      END IF   
!!$!     ==  TEST OVERLAP  =============================================
!!$      PRINT*,'OVERLAP'
!!$      DO I=1,NPHI
!!$        DO J=I,NPHI
!!$          CAUX(:)=(0.D0,0.D0)
!!$          DO IF=1,NF
!!$            CAUX(:IRC)=CAUX(:IRC)+CONJG(PHIL(:IRC,IF,I))*PHIL(:IRC,IF,J)
!!$            CAUX(IRC+1:)=CAUX(IRC+1:)+CONJG(PHIR(IRC+1:,IF,I))*PHIR(IRC+1:,IF,J)
!!$          ENDDO
!!$          CAUX(:)=CAUX(:)*R(:)**2
!!$          AUX(:)=REAL(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$          AUX(:)=AIMAG(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
!!$          OV(I,J)=CMPLX(SVAR1,SVAR2,KIND=8)
!!$          OV(J,I)=CONJG(OV(I,J))
!!$        ENDDO
!!$      ENDDO
!!$      DO I=1,NPHI
!!$        SVAR=1.D0/SQRT(ABS(OV(I,I)))
!!$        OV(:,I)=OV(:,I)*SVAR
!!$        OV(I,:)=OV(I,:)*SVAR
!!$      ENDDO
!!$PRINT*,'MARKE A',NPHI
!!$      DO I=1,NPHI
!!$        WRITE(*,FMT='(20("(",F10.3,",",F10.3,")"))')OV(I,:)
!!$      ENDDO
!!$PRINT*,'MARKE B',NPHI
!!$!
!!$!     ==  PLOT WAVE FUNCTIONS ======================================
!!$      DO IF1=1,NPHI
!!$        WRITE(FILE,*)IF1
!!$        FILE='PHI'//ADJUSTL(FILE)
!!$        OPEN(8,FILE=FILE,FORM='FORMATTED')
!!$        REWIND 8
!!$        DO IR=1,IRC+1
!!$          WRITE(8,FMT='(80F12.7)')R(IR),REAL(PHIL(IR,:,IF1)),AIMAG(PHIL(IR,:,IF1))
!!$        ENDDO
!!$        DO IR=IRC-1,MIN(IROUT+1,NR)
!!$          WRITE(8,FMT='(80F12.7)')R(IR),REAL(PHIR(IR,:,IF1)),AIMAG(PHIR(IR,:,IF1))
!!$        ENDDO
!!$        CLOSE(8)
!!$      ENDDO
!
!     ==================================================================
!     ==  DETERMINE PHIDOT                                            ==
!     ==================================================================
!PRINT*,'SCHROEDINGER_XXXC: DETERMINE PHIDOT'
      PHIL_DOT(:,:,:)=0.D0
      PHIR_DOT(:,:,:)=0.D0
      DO IF=1,NPHI
        CALL RADIAL$DGLGENC(GID,NR,NF,1,IRC+1,A,B,C,-2.D0*PHIL(:,:,IF) &
     &                     ,PHIL_DOT(:,:,IF))
        CALL RADIAL$DGLGENC(GID,NR,NF,IROUT+1,IRC-1,A,B,C,-2.D0*PHIR(:,:,IF) &
     &                     ,PHIR_DOT(:,:,IF))
      ENDDO
!
!     ==================================================================
!     ==  MAKE PHI_DOT CONTINUOUS                                     ==
!     ==================================================================
!PRINT*,'SCHROEDINGER_XXXC: MATCH PHIDOT'
      MAT(:,:)=ALLPHIR(IRC,:,:)
      BVECS(:,:NPHI)=PHIL_DOT(IRC,:,:)-PHIR_DOT(IRC,:,:)
      CALL LIB$MATRIXSOLVEC8(NF,NF,NPHI,MAT,XVECS(:,:NPHI),BVECS(:,:NPHI))
!
      DO IF1=1,NPHI
        DO IF2=1,NF
          PHIR_DOT(IRC-1:IROUT+1,:,IF1)=PHIR_DOT(IRC-1:IROUT+1,:,IF1) &
     &                            +ALLPHIR(IRC-1:IROUT+1,:,IF2)*XVECS(IF2,IF1)
        ENDDO
      ENDDO
!
!     ==================================================================
!     == REMOVE KINKS IN NON-RELEVANT ANGULAR MOMENTUM CHANNELS        ==
!     ==================================================================
!PRINT*,'REMOVE KINKS OF OTHER ANGULAR MOMENTUM CHANNELS OF PHIDOT'
      I=0
      DO IF1=1,NF
        IF(LOX(IF1).EQ.L0) CYCLE
        I=I+1
        DO J=1,NPHI
          CSVAR=(PHIR_DOT(IRC+1,IF1,J)-PHIR_DOT(IRC-1,IF1,J)) &
     &         -(PHIL_DOT(IRC+1,IF1,J)-PHIL_DOT(IRC-1,IF1,J))
          HB(I,J)=-CSVAR
        ENDDO
      ENDDO
      IF(NF.GT.NPHI) THEN
        CALL LIB$MATRIXSOLVEC8(NF-NPHI,NF-NPHI,NPHI,HA,HX,HB)
      END IF
      DO I=1,NPHI
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).EQ.L0) CYCLE
          J=J+1
          PHIL_DOT(:IRC+1,:,I)=PHIL_DOT(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF2)*HX(J,I)
          PHIR_DOT(IRC-1:IROUT+1,:,I)=PHIR_DOT(IRC-1:IROUT+1,:,I) &
     &                               +ALLPHIR(IRC-1:IROUT+1,:,IF2)*HX(J,I)
        ENDDO
      ENDDO  
!
!     ====================================================================
!     ==   ORTHOGONALIZE PHIDOT TO PHI                                  ==
!     ====================================================================
      DO I=1,NPHI
        DO J=1,NPHI
          CAUX(:)=0.D0
          DO IF=1,NF
            CAUX(:IRC)=CAUX(:IRC)+CONJG(PHIL(:IRC,IF,J))*PHIL_DOT(:IRC,IF,I)
            CAUX(IRC+1:IROUT+1)=CAUX(IRC+1:IROUT+1) &
     &             +CONJG(PHIR(IRC+1:IROUT+1,IF,J))*PHIR_DOT(IRC+1:IROUT+1,IF,I)
          ENDDO
          CAUX(:)=CAUX(:)*R(:)**2
          AUX(:)=REAL(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          AUX(:)=AIMAG(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
          CSVAR=CMPLX(SVAR1,SVAR2,KIND=8)
          PHIL_DOT(:IRC+1,:,I)=PHIL_DOT(:IRC+1,:,I)-PHIL(:IRC+1,:,J)*CSVAR
          PHIR_DOT(IRC-1:IROUT+1,:,I)=PHIR_DOT(IRC-1:IROUT+1,:,I)-PHIR(IRC-1:IROUT+1,:,J)*CSVAR
        ENDDO
      ENDDO
!
!     ====================================================================
!     ==   TEST CONTINUITY                                              ==
!     ====================================================================
!!$PRINT*,'SCHROEDINGER_XXXC: TEST CONTINUITY OF PHIDOT'
!!$      SVAR2=0.D0
!!$      DO IF1=1,NPHI
!!$        DO IF2=1,NF
!!$          SVAR1=ABS(PHIL_DOT(IRC,IF2,IF1)-PHIR_DOT(IRC,IF2,IF1))
!!$          SVAR2=MAX(SVAR2,SVAR1)
!!$        ENDDO
!!$      ENDDO
!!$PRINT*,'DEVIATION IN VALUE OF PHIDOT ',SVAR2
!!$      IF(SVAR2.GT.1.D-6) THEN
!!$        CALL ERROR$MSG('MATCHING FAILED')
!!$        CALL ERROR$STOP('SCHROEDINGER_XXXC')
!!$      END IF   
!!$!
!!$!     == TEST DIFFERENTIABILITY ===================================
!!$      SVAR2=0.D0
!!$      DO IF1=1,NPHI
!!$        DO IF2=1,NF
!!$          IF(LOX(IF2).EQ.L0) CYCLE
!!$          SVAR1=ABS( (PHIR_DOT(IRC+1,IF2,IF1)-PHIR_DOT(IRC-1,IF2,IF1)) &
!!$     &              -(PHIL_DOT(IRC+1,IF2,IF1)-PHIL_DOT(IRC-1,IF2,IF1)))
!!$          SVAR2=MAX(SVAR2,SVAR1)
!!$        ENDDO
!!$      ENDDO
!!$PRINT*,'DEVIATION IN OTHER DERIVATIVES OF PHIDOT',SVAR2
!!$      IF(SVAR2.GT.1.D-4) THEN
!!$        CALL ERROR$MSG('MATCHING OF DPHIDOT/DR FAILED')
!!$        CALL ERROR$STOP('SCHROEDINGER_XXXC')
!!$      END IF   
!!$!
!!$!     =============================================================
!!$      DO IF1=1,NPHI
!!$        WRITE(FILE,*)IF1
!!$        FILE='PHIDOT'//ADJUSTL(FILE)
!!$        OPEN(8,FILE=FILE,FORM='FORMATTED')
!!$        REWIND 8
!!$        DO IR=1,IRC+1
!!$          WRITE(8,FMT='(80F12.7)')R(IR),REAL(PHIL_DOT(IR,:,IF1)),AIMAG(PHIL_DOT(IR,:,IF1))
!!$        ENDDO
!!$        DO IR=IRC-1,MIN(IROUT+1,NR)
!!$          WRITE(8,FMT='(80F12.7)')R(IR),REAL(PHIR_DOT(IR,:,IF1)),AIMAG(PHIR_DOT(IR,:,IF1))
!!$        ENDDO
!!$        CLOSE(8)
!!$      ENDDO
!
!     ==================================================================
!     ==  REMOVE KINKS BY MIXING PHIDOT INTO PHI                      ==
!     ==================================================================
!PRINT*,'SCHROEDINGER_XXXC:MATCH KINKS'
      I=0
      DO IF2=1,NF
        IF(LOX(IF2).NE.L0) CYCLE
        I=I+1
        KINK_HOM(I,:)=(PHIR(IRC+1,IF2,:)-PHIR(IRC-1,IF2,:)) &
     &               -(PHIL(IRC+1,IF2,:)-PHIL(IRC-1,IF2,:))
        KINK_DOT(I,:)=(PHIR_DOT(IRC+1,IF2,:)-PHIR_DOT(IRC-1,IF2,:)) &
     &               -(PHIL_DOT(IRC+1,IF2,:)-PHIL_DOT(IRC-1,IF2,:))
      ENDDO
      CALL LIB$MATRIXSOLVEC8(NPHI,NPHI,NPHI,-KINK_DOT,HAM,KINK_HOM)
!!$      SVAR=MAXVAL(ABS(KINK_HOM+MATMUL(KINK_DOT,HAM)))
!!$PRINT*,'REMAINING KINKS ',SVAR
!     == MAKE PHI DIFFERENTIABLE =========================================
      DO I=1,NPHI
        DO J=1,NPHI
          PHIL(:IRC,:,I)         =PHIL(:IRC,:,I)         +PHIL_DOT(:IRC,:,J)  *HAM(J,I)
          PHIR(IRC+1:IROUT+1,:,I)=PHIR(IRC+1:IROUT+1,:,I)+PHIR_DOT(IRC+1:IROUT+1,:,J)*HAM(J,I)
        ENDDO
      ENDDO
!     == MAKE PHI DIFFERENTIABLE =========================================
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("HAM ",50("(",2F10.5")   "))')HAM(I,:)
!!$ENDDO
!SVAR=0.5D0*MAXVAL(ABS(HAM-TRANSPOSE(CONJG(HAM))))
!PRINT*,'DEVIATION FROM HERMITICITY',SVAR
!
!     ==================================================================
!     ==  DETERMINE OVERLAP MATRIX                                    ==
!     ==================================================================
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("H ",50("(",2F10.5")   "))')HAM(I,:)
!!$ENDDO
!!$
!!$      DO I=1,NPHI
!!$        DO J=I,NPHI
!!$          CAUX(:)=0.D0
!!$          DO IF=1,NF
!!$            CAUX(:IRC)=CAUX(:IRC)+CONJG(PHIL_DOT(:IRC,IF,I))*PHIL_DOT(:IRC,IF,J)
!!$            CAUX(IRC+1:IROUT+1)=CAUX(IRC+1:IROUT+1) &
!!$     &             +CONJG(PHIR_DOT(IRC+1:IROUT+1,IF,I))*PHIR_DOT(IRC+1:IROUT+1,IF,J)
!!$          ENDDO
!!$          CAUX(:)=CAUX(:)*R(:)**2
!!$          AUX(:)=REAL(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$          AUX(:)=AIMAG(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
!!$          OV(I,J)=CMPLX(SVAR1,SVAR2,KIND=8)
!!$          OV(J,I)=CONJG(OV(I,J))
!!$        ENDDO
!!$      ENDDO
!!$      OV=MATMUL(TRANSPOSE(CONJG(HAM)),MATMUL(OV,HAM))
!!$      DO I=1,NPHI
!!$        OV(I,I)=OV(I,I)+(1.D0,0.D0)
!!$      ENDDO
!!$!
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("O ",50("(",2F10.5")   "))')OV(I,:)
!!$ENDDO
!
      DO I=1,NPHI
        DO J=I,NPHI
          CAUX(:)=0.D0
          DO IF=1,NF
            CAUX(:IRC)=CAUX(:IRC)+CONJG(PHIL(:IRC,IF,I))*PHIL(:IRC,IF,J)
            CAUX(IRC+1:IROUT+1)=CAUX(IRC+1:IROUT+1) &
     &             +CONJG(PHIR(IRC+1:IROUT+1,IF,I))*PHIR(IRC+1:IROUT+1,IF,J)
          ENDDO
          CAUX(:)=CAUX(:)*R(:)**2
          AUX(:)=REAL(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          AUX(:)=AIMAG(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
          OV(I,J)=CMPLX(SVAR1,SVAR2,KIND=8)
          OV(J,I)=CONJG(OV(I,J))
        ENDDO
      ENDDO
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("O ",50("(",2F10.5")   "))')OV(I,:)
!!$ENDDO
!
!     ==================================================================
!     ==  DETERMINE EIGENSTATES                                       ==
!     ==================================================================
      HAM=0.5D0*(HAM+TRANSPOSE(CONJG(HAM)))
!      CALL LIB$DIAGC8(NPHI,HAM,DE,KINKC)
      CALL LIB$GENERALEIGENVALUEC8(NPHI,HAM,OV,DE,KINKC)
      PHI(:,:,:)=0.D0
      DO I=1,NPHI
        DO J=1,NPHI
          PHI(:IRC,:,I)         =PHI(:IRC,:,I)         +PHIL(:IRC,:,J)*KINKC(J,I)      
          PHI(IRC+1:IROUT+1,:,I)=PHI(IRC+1:IROUT+1,:,I)+PHIR(IRC+1:IROUT+1,:,J)*KINKC(J,I)      
        ENDDO
      ENDDO
!PRINT*,'DE',DE
!!$!
!!$!     ==================================================================
!!$!     ==  NORMALIZE SOLUTIONS                                         ==
!!$!     ==================================================================
!!$      DO I=1,NPHI
!!$        AUX(:)=0.D0
!!$        DO J=1,NF
!!$          AUX(:)=AUX(:)+ABS(PHI(:,J,I))**2
!!$        ENDDO
!!$        AUX(:)=AUX(:)*R(:)**2
!!$        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$        SVAR=1.D0/SQRT(SVAR)
!!$        PHI(:,:,I)=PHI(:,:,I)*SVAR
!!$      ENDDO
!
!     ==================================================================
!     ==  TEST OVERLAP                                                ==
!     ==================================================================
!!$      DO IF1=1,NPHI
!!$        WRITE(FILE,*)IF1
!!$        FILE='PHIFIN'//ADJUSTL(FILE)
!!$        OPEN(8,FILE=FILE,FORM='FORMATTED')
!!$        REWIND 8
!!$        DO IR=1,IROUT+1
!!$          WRITE(8,FMT='(80F20.7)')R(IR),REAL(PHI(IR,:,IF1)),AIMAG(PHI(IR,:,IF1))
!!$        ENDDO
!!$        CLOSE(8)
!!$      ENDDO
!!$!
!     == TEST ORTHONORMALITY OVERLAP MATRIX
!!$      DO I=1,NPHI
!!$        HAM(I,I)=(0.D0,0.D0)
!!$        DO J=I,NPHI
!!$          CAUX(:)=0.D0
!!$          DO IF=1,NF
!!$            CAUX(:)=CAUX(:)+CONJG(PHI(:,IF,I))*PHI(:,IF,J)
!!$          ENDDO
!!$          CAUX(:)=CAUX(:)*R(:)**2
!!$          AUX=REAL(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$          AUX=AIMAG(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
!!$          HAM(I,J)=CMPLX(SVAR1,SVAR2,KIND=8)
!!$          HAM(J,I)=CONJG(HAM(I,J))
!!$        ENDDO
!!$        HAM(I,I)=HAM(I,I)-(1.D0,0.D0)
!!$      ENDDO
!!$      DO I=1,NPHI
!!$        WRITE(*,FMT='("O ",50("(",2F10.5")   "))')HAM(I,:)
!!$      ENDDO
         SVAR=MAXVAL(ABS(HAM))
         PRINT*,'DEVIATION FROM ORTHONORMALITY',SVAR,NF
!!$      IF(SVAR.GT.0.5D0) THEN
!!$        CALL ERROR$MSG('WAVE FUNCTIONS NOT ORTHOGONAL')
!!$        CALL ERROR$STOP('SCHROEDINGER_XXXC')
!!$      END IF
      TOK=.TRUE.
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER$LBND_SLOC(GID,NR,LX,LMX,LMRX,TMAINSH,POT &
     &                             ,LXSL,POTSL,OVSL,G,ENU,NPHI,EB,PHI,TPHI,TOK)
!     **                                                                      **
!     **  SOLVES THE NONRELATIVISTIC RADIAL SCHROEDINGER EQUATION FOR THE     **
!     **  G IS AN INHOMOGENEITY, WHICH MUST BE SET TO ZERO FOR THE            **
!     **  HOMOGENEOUS SOLUTION.                                               **
!     **                                                                      **
!     **  THE SOLUTIONS ARE CALCULATED TO LINEAR ORDER IN DE                  **
!     **  WHERE E+DE IS THE NEW ENERGY OF THE WAVE FUNCTIONS                  **
!     **                                                                      **
!     **  THE DIFFERENTIAL EQUATION IS SOLVED WITH THE VERLET ALGORITHM       **
!     **    DPHI/DX=(PHI(+)-PHI(-))/2                                         **
!     **    D2PHI/DX2=PHI(+)-2PHI(0)+PHI(-)                                   **
!     **  WHERE X IS THE VARIABLE WITH R(X=I)=R_I                             **
!     **                                                                      **
!     **  ATTENTION! THE ROUTINE IS NOT GUARDED AGAINST OVERFLOW DUE          **
!     **    TO THE EXPONENTIAL INCREASE OF THE SOLUTION                       **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID-ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: LX      ! X#(ANGULAR MOMENTUM OF WAVEF)
      INTEGER(4) ,INTENT(IN)     :: LMX     ! X#(WAVE FUNCTION ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)     :: LMRX    ! X#(POTENTIAL ANGULAR MOMENTA) 
      LOGICAL(4) ,INTENT(IN)     :: TMAINSH(LX+1)
      INTEGER(4) ,INTENT(IN)     :: LXSL    ! #(SEMILOCAL POTENTIALS, ONE PER L)
      REAL(8)    ,INTENT(IN)     :: POTSL(NR,LXSL)! SEMI-LOCAL POTENTIALS
      REAL(8)    ,INTENT(IN)     :: OVSL(NR,LXSL)  ! SEMI-LOCAL OVERLAP CONTRIBTION
      REAL(8)    ,INTENT(IN)     :: G(NR,LMX)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: ENU     ! EXPANSION ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR,LMRX) !POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: NPHI              ! #(WAVE FUNCTIONS)
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,LMX,NPHI)  ! WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: TPHI(NR,LMX,NPHI) ! P**2/(2M)*WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: EB(NPHI)          ! ONE-PARTICKE ENERGIES
      LOGICAL(4) ,INTENT(OUT)    :: TOK               ! ERROR FLAG
      INTEGER(4)                 :: LOX(LMX) ! ANGULAR MOMENTA
      INTEGER(4)                 :: LM,LM1,LM2,LM3,L,M,IM,IB
      REAL(8)                    :: A(NR)
      REAL(8)                    :: B(NR)
      REAL(8)                    :: C(NR,LMX,LMX)
      REAL(8)                    :: DCDE(NR,LMX)
      REAL(8)                    :: D(NR,LMX) 
      REAL(8)                    :: R(NR)                        !
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: CPOT(NR,LMX,LMX)
      REAL(8)                    :: CG
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20
      INTEGER(4)                 :: IRCL,IROUT
      INTEGER(4)                 :: ISVAR
      LOGICAL(4)                 :: TCHK
!     **************************************************************************
      TOK=.FALSE.
      IF((LX+1)**2.NE.LMX) THEN   !LMX=(LX+1)**2
        CALL ERROR$MSG('LMX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMX',LMX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('SCHROEDINGER$LBND_SLOC')
      END IF
      ISVAR=0
      DO L=0,LX
        IF(TMAINSH(L+1)) ISVAR=ISVAR+2*L+1
      ENDDO
      IF(ISVAR.NE.NPHI) THEN
        CALL ERROR$MSG('TMAINSH AND NPHI ARE INCONSISTENT')
        CALL ERROR$I4VAL('NPHI',NPHI)
        CALL ERROR$I4VAL('LX',LX)
        PRINT*,'TMAINSH ',TMAINSH
        CALL ERROR$STOP('SCHROEDINGER$LBND_SLOC')
      END IF
      CALL RADIAL$R(GID,NR,R)
      LM=0
      DO L=0,LX
        DO M=1,2*L+1
          LM=LM+1
          LOX(LM)=L
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  DETERMINE CLASSICAL TURNING POINT R(IRCL)                           ==
!     ==  AND OUTMOST POINT FOR INWARD INTEGRATION R(IROUT)                   ==
!     ==========================================================================
      CALL SCHROEDINGER_SPECIALRADS(GID,NR,0,XMAX,POT(:,1),ENU,IRCL,IROUT)
!
!     ==========================================================================
!     ==  PREPARE POTENTIAL-INDEPENDENT ARRAYS                                ==
!     ==========================================================================
!     == A*D2F/DR2+B*DF/DR+C*F=D
      A(:)=1.D0
!     == AVOID DIVIDE BY ZERO IF THE FIRST GRID POINT IS THE ORIGIN.
!     == THE FORCES ON THE FIRST GRID POINT ARE NOT USED,
!     == BECAUSE RADIAL$DGL IS BASED ON THE VERLET ALGORITHM
!     == THAT CANNOT USE THE FORCES ON THE FIRST AND LAST GRID POINT
      B(2:)=2.D0/R(2:)
      B(1)=B(2)
      D(:,:)=-2.D0*G(:,:)
!
!     ==========================================================================
!     ==  COUPLING BETWEEN WAVE FUNCTION COMPONENTS VIA POTENTIAL             ==
!     ==========================================================================
      DCDE=0.D0
      C(:,:,:)=0.D0
      DO LM1=1,LMX
        DO LM2=1,LMX
          AUX(:)=0.D0
          DO LM3=1,LMRX
            CALL CLEBSCH(LM1,LM2,LM3,CG)
            IF(CG.EQ.0.D0) CYCLE
            AUX=CG*POT(:,LM3)
            IF(LM3.EQ.1) AUX(IROUT+1:)=AUX(IROUT) !CONSTANT POTENTIAL BEYOND IROUT
            C(:,LM1,LM2)=C(:,LM1,LM2)+AUX
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  ADD SEMI-LOCAL POTENTIAL                                            ==
!     ==========================================================================
      CPOT=C   ! KEEP POTENTIAL TO DETERMINE THE KINETIC ENERGY
!     == ATTENTION! LXSL IS THE NUMBER OF ANGULAR MOMENTA, 
!     == AND NOT THE HIGHEST ANGULAR MOMENTUM
      LM=0
      DO L=0,LXSL-1
        DO M=1,2*L+1
          LM=LM+1
          IF(LM.GT.LMX) EXIT
          C(:,LM,LM)=C(:,LM,LM)+POTSL(:,L+1)-ENU*OVSL(:,L+1)
          CPOT(:,LM,LM)=CPOT(:,LM,LM)+POTSL(:,L+1)
          DCDE(:,LM)=-OVSL(:,L+1)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  ADD KINETIC ENERGY TERM TO C AND SHIFT ENERGY ZERO TO ENU           ==
!     ==========================================================================
      LM=0
      DO L=0,LX
        AUX(1)=0.D0
        AUX(2:)=0.5D0*REAL(L*(L+1),KIND=8)/R(2:)**2 - ENU
        DO IM=1,2*L+1
          LM=LM+1
          C(:,LM,LM)=C(:,LM,LM)+AUX(:)
          DCDE(:,LM)=DCDE(:,LM)-1.D0
        ENDDO
      ENDDO
!     ==  AVOID DIVIDE-BY-ZERO
      C(1,:,:)=C(2,:,:)
      DCDE(1,:)=DCDE(2,:)
!
!     ==========================================================================
!     ==  SCALE C                                                             ==
!     ==========================================================================
      C=-2.D0*C
      DCDE=-2.D0*DCDE
!
!     ==========================================================================
!     ==  DETERMINE BOUND STATES                                              ==
!     ==========================================================================
!      CALL SCHROEDINGER_XXXR(GID,NR,LMX,IRCL,IROUT,LOX,A,B,C,D,NPHI,EB,PHI,TCHK)
      CALL SCHROEDINGER_XXXR_OV(GID,NR,LX,LMX,IRCL,IROUT,TMAINSH,A,B,C,DCDE,D,NPHI,EB,PHI,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$STOP('SCHROEDINGER_XXXR_OV FINISHED WITH ERROR')
        CALL ERROR$STOP('SCHROEDINGER$LBND_SLOC')
      END IF
!
!     ==========================================================================
!     ==  SHIFT ENERGIES                                                      ==
!     ==========================================================================
      EB(:)=ENU+EB(:) ! CHANGE ENERGIES RELATIVE TO ENY TO ABSOLUTE ENERGIES
!
!     ==========================================================================
!     ==  DETERMINE TPHI                                                      ==
!     ==========================================================================
!     -- BETTER DIRECTLY WORK OUT THE KINETIC ENERGY BECAUSE  THE 
!     -- SCHROEDINGER EQUATION IS FULFILLED ONLY TO FIRST ORDER IN DE
      TPHI(:,:,:)=0.D0
      DO IB=1,NPHI
        DO LM1=1,LMX
!         == DCDE=-(1+O)========================================================
          TPHI(:,LM1,IB)=TPHI(:,LM1,IB)+0.5D0*DCDE(:,LM)*PHI(:,LM1,IB)*EB(IB)
          DO LM2=1,LMX
            TPHI(:,LM1,IB)=TPHI(:,LM1,IB)-CPOT(:,LM1,LM2)*PHI(:,LM2,IB)
          ENDDO
        ENDDO
      ENDDO
!
      TOK=.TRUE.
      RETURN
      END 
!
!     ..........................................................................
      SUBROUTINE SCHROEDINGER_XXXR_OV(GID,NR,LX,NF,IRMATCH,IROUT,TMAINSH,A,B,C,DCDE,D,NPHI,DE,PHI,TOK)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID      ! GRID-ID FOR RADIAL GRID
      INTEGER(4),INTENT(IN) :: NR       ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: LX       ! MAX ANGULAR MOMENTUM OF LM EXPANSION
      INTEGER(4),INTENT(IN) :: NF             ! #(ANGULAR MOMENTA)
      INTEGER(4),INTENT(IN) :: IRMATCH        ! MATCHING POINT FOR INSIDE-OUTSIDE INTEGRATION
      INTEGER(4),INTENT(IN) :: IROUT          ! OUTERMOST POINT TO BE CONSIDERED
      LOGICAL(4),INTENT(IN) :: TMAINSH(LX+1)  ! 
      REAL(8)   ,INTENT(IN) :: A(NR)
      REAL(8)   ,INTENT(IN) :: B(NR)
      REAL(8)   ,INTENT(IN) :: C(NR,NF,NF)
      REAL(8)   ,INTENT(IN) :: DCDE(NR,NF)     ! ENERGY DERIVATIVE OF C
      REAL(8)   ,INTENT(IN) :: D(NR,NF)        ! INHOMOGENEITY
      INTEGER(4),INTENT(IN) :: NPHI            ! #(WAVE FUNCTIONS)
      REAL(8)   ,INTENT(OUT):: DE(NPHI)        ! ONE-PARTICLE EIGENVALUES
      REAL(8)   ,INTENT(OUT):: PHI(NR,NF,NPHI) ! WAVE FUNCTIONS
      LOGICAL(4),INTENT(OUT):: TOK             ! ERROR FLAG
      LOGICAL(4)            :: TMAIN(NF)       ! SELECTS THE MAIN CONTRIBUTIONS
      INTEGER(4)            :: LOX(NF)         ! ANGULAR MOMENTA OF WAVE FUNCTION COMPONENTS
      REAL(8)               :: G(NR,NF)
      REAL(8)               :: ALLPHIL(NR,NF,NF)
      REAL(8)               :: ALLPHIR(NR,NF,NF)
      REAL(8)               :: PHIL(NR,NF,NPHI)
      REAL(8)               :: PHIR(NR,NF,NPHI)
      REAL(8)               :: PHIL_DOT(NR,NF,NPHI)
      REAL(8)               :: PHIR_DOT(NR,NF,NPHI)
      REAL(8)               :: PHIWORK2D(NR,NF,NF)
      INTEGER(4)            :: IF,IF1,IF2
      REAL(8)               :: HA(NF-NPHI,NF-NPHI),HX(NF-NPHI,NPHI),HB(NF-NPHI,NPHI)
      REAL(8)               :: MAT(NF,NF)
      INTEGER(4)            :: IRC
      REAL(8)               :: KINK_HOM(NPHI,NPHI)
      REAL(8)               :: KINK_DOT(NPHI,NPHI)
      REAL(8)               :: KINKC(NPHI,NPHI)
      REAL(8)               :: HAM(NPHI,NPHI),OV(NPHI,NPHI)
      REAL(8)               :: R(NR)
      REAL(8)               :: DHOM(NR,NF)
      LOGICAL   ,PARAMETER  :: TWRITE=.FALSE.
      REAL(8)               :: BVECS(NF,NF),XVECS(NF,NF)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: SVAR
      INTEGER(4)            :: I,J,M,L,LM
CHARACTER(32):: FILE
!     **************************************************************************
      TOK=.FALSE.
      IF(IROUT+1.GT.NR) THEN
        CALL ERROR$MSG('IROUT OUT OF RANGE')
        CALL ERROR$STOP('SCHROEDINGER_XXXR_OV')
      END IF
      IF(IRMATCH.GT.IROUT) THEN
        CALL ERROR$MSG('IRMATCH OUT OF RANGE')
        CALL ERROR$STOP('SCHROEDINGER_XXXR_OV')
      END IF
      IF(NF.NE.(LX+1)**2) THEN
        CALL ERROR$MSG('NF INCONSISTENT WITH LX')
        CALL ERROR$STOP('SCHROEDINGER_XXXR_OV')
      END IF
!
!     ==========================================================================
!     ==  INITIALIZATIONS                                                     ==
!     ==========================================================================
      LM=0
      DO L=0,LX
        DO M=1,2*L+1
          LM=LM+1
          LOX(LM)=L
        ENDDO
      ENDDO
      TMAIN(:)=.FALSE.
      DO LM=1,NF
        TMAIN(LM)=TMAINSH(LOX(LM)+1)
      ENDDO
!
      CALL RADIAL$R(GID,NR,R)
      IRC=IRMATCH
!PRINT*,'IROUT ',IROUT,R(IROUT),IRC,R(IRC)
!
!     ==========================================================================
!     ==  OBTAIN HOMOGENEOUS SOLUTION                                         ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXR_OV: DETERMINE PHI',LOX
      ALLPHIL(:,:,:)=0.D0
      ALLPHIR(:,:,:)=0.D0
      DO IF=1,NF
        ALLPHIL(1:2,IF,IF)=R(1:2)**LOX(IF)
        DHOM(:,:)=0.D0
        CALL RADIAL$DGLGEN(GID,NR,NF,1,IRC+1,A,B,C,DHOM,ALLPHIL(:,:,IF))
        SVAR=MAXVAL(ABS(ALLPHIL(IRC,:,IF)))
        ALLPHIL(1:IRC+1,:,IF)=ALLPHIL(1:IRC+1,:,IF)/SVAR
        DHOM(IROUT,IF)=1.D-8
        CALL RADIAL$DGLGEN(GID,NR,NF,IROUT+1,IRC-1,A,B,C,DHOM,ALLPHIR(:,:,IF))
        SVAR=MAXVAL(ABS(ALLPHIR(IRC,:,IF)))
        ALLPHIR(IRC-1:IROUT+1,:,IF)=ALLPHIR(IRC-1:IROUT+1,:,IF)/SVAR
      ENDDO
!
!     ==========================================================================
!     ==  MAKE PHI_HOM SOLUTIONS CONTINUOUS                                   ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXR_OV: MATCH PHI'
      MAT(:,:)=ALLPHIR(IRC,:,:)
!     == MATRIX IS NOT SYMMETRIC. THUS SOLVE WITH SINGULAR VALUE DECOMPOSITION
      BVECS(:,1:NF)=ALLPHIL(IRC,:,:)
      CALL LIB$MATRIXSOLVER8(NF,NF,NF,MAT,XVECS,BVECS)
      PHIWORK2D(:,:,:)=ALLPHIR(:,:,:)
      ALLPHIR(:,:,:)=0.D0
      DO IF1=1,NF
        DO IF2=1,NF
          ALLPHIR(:,:,IF1)=ALLPHIR(:,:,IF1)+PHIWORK2D(:,:,IF2)*XVECS(IF2,IF1)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == REMOVE KINKS IN NON-RELEVANT ANGULAR MOMENTA CHANNELS                ==
!     ==========================================================================
!PRINT*,'REMOVE KINKS OF OTHER ANGULAR MOMENTUM CHANNELS OF PHI'
      I=0
      DO IF1=1,NF
        IF(TMAIN(IF1)) CYCLE
        I=I+1
        J=0
        DO IF2=1,NF
          IF(TMAIN(IF2)) CYCLE
          J=J+1
          SVAR=(ALLPHIR(IRC+1,IF1,IF2)-ALLPHIR(IRC-1,IF1,IF2)) &
     &         -(ALLPHIL(IRC+1,IF1,IF2)-ALLPHIL(IRC-1,IF1,IF2))
          HA(I,J)=SVAR
        ENDDO
        J=0
        DO IF2=1,NF
          IF(.NOT.TMAIN(IF2)) CYCLE
          J=J+1
          SVAR=(ALLPHIR(IRC+1,IF1,IF2)-ALLPHIR(IRC-1,IF1,IF2)) &
     &         -(ALLPHIL(IRC+1,IF1,IF2)-ALLPHIL(IRC-1,IF1,IF2))
          HB(I,J)=-SVAR
        ENDDO
      ENDDO
      IF(NF.GT.NPHI) THEN
        CALL LIB$MATRIXSOLVER8(NF-NPHI,NF-NPHI,NPHI,HA,HX,HB)
      END IF
      PHIL(:,:,:)=0.D0
      PHIR(:,:,:)=0.D0
      I=0
      DO IF1=1,NF
        IF(.NOT.TMAIN(IF1)) CYCLE
        I=I+1
        PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF1)
        PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I) &
     &                         +ALLPHIR(IRC-1:IROUT+1,:,IF1)
        J=0
        DO IF2=1,NF
          IF(TMAIN(IF2)) CYCLE
          J=J+1
          PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF2)*HX(J,I)
          PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I) &
     &                            +ALLPHIR(IRC-1:IROUT+1,:,IF2)*HX(J,I)
        ENDDO
      ENDDO  
!
!     ==========================================================================
!     ==   ORTHOGONALIZE PHI                                                  ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXC: ORTHOGONALIZE PHI'
      DO I=1,NPHI
!       == ORTHOGONALIZE
        DO J=1,I-1
          AUX(:)=0.D0
          DO IF=1,NF
            AUX(:IRC)=AUX(:IRC) &
     &               +0.5D0*DCDE(:IRC,IF)*PHIL(:IRC,IF,J)*PHIL(:IRC,IF,I)
            AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1)+0.5D0*DCDE(IRC+1:IROUT+1,IF) &
     &                        *PHIR(IRC+1:IROUT+1,IF,J)*PHIR(IRC+1:IROUT+1,IF,I)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)-PHIL(:IRC+1,:,J)*SVAR
          PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I)-PHIR(IRC-1:IROUT+1,:,J)*SVAR
        ENDDO
!       ==  NORMALIZE
        AUX(:)=0.D0
        DO IF=1,NF
          AUX(:IRC)         =AUX(:IRC)  +0.5D0*DCDE(:IRC,IF)*PHIL(:IRC,IF,I)**2
          AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
    &                 +0.5D0*DCDE(IRC+1:IROUT+1,IF)*PHIR(IRC+1:IROUT+1,IF,I)**2
        ENDDO
        AUX(:)=AUX(:)*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        SVAR=1.D0/SQRT(SVAR)
        PHIL(:IRC+1,:,I)       =PHIL(:IRC+1,:,I)*SVAR
        PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I)*SVAR
      ENDDO
!
!     ==========================================================================
!     ==  DETERMINE PHIDOT                                                    ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXC: DETERMINE PHIDOT'
      PHIL_DOT(:,:,:)=0.D0
      PHIR_DOT(:,:,:)=0.D0
      DO IF=1,NPHI
        G(:,:)=-DCDE(:,:)*PHIL(:,:,IF)
        CALL RADIAL$DGLGEN(GID,NR,NF,1,IRC+1,A,B,C,G,PHIL_DOT(:,:,IF))
        G(:,:)=-DCDE(:,:)*PHIR(:,:,IF)
        CALL RADIAL$DGLGEN(GID,NR,NF,IROUT+1,IRC-1,A,B,C,G,PHIR_DOT(:,:,IF))
      ENDDO
!
!     ==========================================================================
!     ==  MAKE PHI_DOT CONTINUOUS                                             ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXC: MATCH PHIDOT'
      MAT(:,:)=ALLPHIR(IRC,:,:)
      BVECS(:,:NPHI)=PHIL_DOT(IRC,:,:)-PHIR_DOT(IRC,:,:)
      CALL LIB$MATRIXSOLVER8(NF,NF,NPHI,MAT,XVECS(:,:NPHI),BVECS(:,:NPHI))
!
      DO IF1=1,NPHI
        DO IF2=1,NF
          PHIR_DOT(IRC-1:IROUT+1,:,IF1)=PHIR_DOT(IRC-1:IROUT+1,:,IF1) &
     &                            +ALLPHIR(IRC-1:IROUT+1,:,IF2)*XVECS(IF2,IF1)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == REMOVE KINKS IN NON-RELEVANT ANGULAR MOMENTUM CHANNELS               ==
!     ==========================================================================
!PRINT*,'REMOVE KINKS OF OTHER ANGULAR MOMENTUM CHANNELS OF PHIDOT'
      I=0
      DO IF1=1,NF
        IF(TMAIN(IF1)) CYCLE
        I=I+1
        DO J=1,NPHI
          SVAR=(PHIR_DOT(IRC+1,IF1,J)-PHIR_DOT(IRC-1,IF1,J)) &
     &         -(PHIL_DOT(IRC+1,IF1,J)-PHIL_DOT(IRC-1,IF1,J))
          HB(I,J)=-SVAR
        ENDDO
      ENDDO
      IF(NF.GT.NPHI) THEN
        CALL LIB$MATRIXSOLVER8(NF-NPHI,NF-NPHI,NPHI,HA,HX,HB)
      END IF
      DO I=1,NPHI
        J=0
        DO IF2=1,NF
          IF(TMAIN(IF2)) CYCLE
          J=J+1
          PHIL_DOT(:IRC+1,:,I)=PHIL_DOT(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF2)*HX(J,I)
          PHIR_DOT(IRC-1:IROUT+1,:,I)=PHIR_DOT(IRC-1:IROUT+1,:,I) &
     &                               +ALLPHIR(IRC-1:IROUT+1,:,IF2)*HX(J,I)
        ENDDO
      ENDDO  
!
!     ==========================================================================
!     ==   ORTHOGONALIZE PHIDOT TO PHI                                        ==
!     ==========================================================================
      DO I=1,NPHI
        DO J=1,NPHI
          AUX(:)=0.D0
          DO IF=1,NF
            AUX(:IRC)=AUX(:IRC) &
     &                  +0.5D0*DCDE(:IRC,IF)*PHIL(:IRC,IF,J)*PHIL_DOT(:IRC,IF,I)
            AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1)+0.5D0*DCDE(IRC+1:IROUT+1,IF) &
     &                    *PHIR(IRC+1:IROUT+1,IF,J)*PHIR_DOT(IRC+1:IROUT+1,IF,I)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          PHIL_DOT(:IRC+1,:,I)=PHIL_DOT(:IRC+1,:,I)-PHIL(:IRC+1,:,J)*SVAR
          PHIR_DOT(IRC-1:IROUT+1,:,I)=PHIR_DOT(IRC-1:IROUT+1,:,I) &
     &                               -PHIR(IRC-1:IROUT+1,:,J)*SVAR
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  REMOVE KINKS BY MIXING PHIDOT INTO PHI                              ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXR_OV:MATCH KINKS'
      I=0
      DO IF2=1,NF
        IF(.NOT.TMAIN(IF2)) CYCLE
        I=I+1
        KINK_HOM(I,:)=(PHIR(IRC+1,IF2,:)-PHIR(IRC-1,IF2,:)) &
     &               -(PHIL(IRC+1,IF2,:)-PHIL(IRC-1,IF2,:))
        KINK_DOT(I,:)=(PHIR_DOT(IRC+1,IF2,:)-PHIR_DOT(IRC-1,IF2,:)) &
     &               -(PHIL_DOT(IRC+1,IF2,:)-PHIL_DOT(IRC-1,IF2,:))
      ENDDO
      CALL LIB$MATRIXSOLVER8(NPHI,NPHI,NPHI,-KINK_DOT,HAM,KINK_HOM)
!!$      SVAR=MAXVAL(ABS(KINK_HOM+MATMUL(KINK_DOT,HAM)))
!!$PRINT*,'REMAINING KINKS ',SVAR
!     == MAKE PHI DIFFERENTIABLE ===============================================
      DO I=1,NPHI
        DO J=1,NPHI
          PHIL(:IRC,:,I)         =PHIL(:IRC,:,I)         +PHIL_DOT(:IRC,:,J)  *HAM(J,I)
          PHIR(IRC+1:IROUT+1,:,I)=PHIR(IRC+1:IROUT+1,:,I)+PHIR_DOT(IRC+1:IROUT+1,:,J)*HAM(J,I)
        ENDDO
      ENDDO
!     == MAKE PHI DIFFERENTIABLE ===============================================
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("HAM ",50("(",2F10.5")   "))')HAM(I,:)
!!$ENDDO
!SVAR=0.5D0*MAXVAL(ABS(HAM-TRANSPOSE(CONJG(HAM))))
!PRINT*,'DEVIATION FROM HERMITICITY',SVAR
!
!     ==========================================================================
!     ==  DETERMINE OVERLAP MATRIX                                            ==
!     ==========================================================================
      DO I=1,NPHI
        DO J=I,NPHI
          AUX(:)=0.D0
          DO IF=1,NF
            AUX(:IRC)=AUX(:IRC) &
     &              +0.5D0*DCDE(:IRC,IF)*PHIL(:IRC,IF,I)*PHIL(:IRC,IF,J)
            AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
     &          +0.5D0*DCDE(IRC+1:IROUT+1,IF) &
     &                *PHIR(IRC+1:IROUT+1,IF,I)*PHIR(IRC+1:IROUT+1,IF,J)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          OV(I,J)=SVAR
          OV(J,I)=OV(I,J)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  DETERMINE EIGENSTATES                                               ==
!     ==========================================================================
!PRINT*,'SCHROEDINGER_XXXR_OV: DETERMINE EIGENSTATES'
      HAM=0.5D0*(HAM+TRANSPOSE(HAM))
      OV=0.5D0*(OV+TRANSPOSE(OV))
      CALL LIB$GENERALEIGENVALUER8(NPHI,HAM,OV,DE,KINKC)
      PHI(:,:,:)=0.D0
      DO I=1,NPHI
        DO J=1,NPHI
          PHI(:IRC,:,I)         =PHI(:IRC,:,I)         +PHIL(:IRC,:,J)*KINKC(J,I)      
          PHI(IRC+1:IROUT+1,:,I)=PHI(IRC+1:IROUT+1,:,I)+PHIR(IRC+1:IROUT+1,:,J)*KINKC(J,I)      
        ENDDO
      ENDDO
!PRINT*,'SCHROEDINGER_XXXR_OV: DONE'
      TOK=.TRUE.
      RETURN
      END SUBROUTINE SCHROEDINGER_XXXR_OV
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER$MODPOT(GID,NR,XMAX,E,POT)
!     **************************************************************************
!     ** MODIFIES THE POTENTIAL SUCH THAT AN EXPONENTIAL INCREASE DOES DOES   **
!     ** NOT LEAD TO AN OVERFLOW, I.E. THE MAX. RATIO OF WAVE FUNCTIONS       **
!     ** SHOULD NOT EXCEED XMAX. THE OUTER POTENTIAL IS LIMITED TO A MAX.     **
!     ** VALUE                                                                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: GID      ! GRID ID
      INTEGER(4),INTENT(IN)   :: NR       ! #(RADIAL GRID POINTS)
      REAL(8)   ,INTENT(IN)   :: XMAX     ! MAXIMUM TOLERABLE RATIO OF PHI 
      REAL(8)   ,INTENT(INOUT):: POT(NR)  ! RADIAL SPHERICAL POTENTIAL
      REAL(8)   ,INTENT(IN)   :: E        ! ENERGY
      REAL(8)   ,PARAMETER    :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER    :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                 :: SVAR
      REAL(8)                 :: R(NR)
      REAL(8)                 :: F(NR)
      REAL(8)                 :: G(NR)
      INTEGER(4)              :: IR
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      F(:)=SQRT(MAX(POT(:)*Y0-E,0.D0))
      CALL RADIAL$INTEGRATE(GID,NR,F,G)
      G=G+F*(R(NR)-R(:))-LOG(XMAX)
      DO IR=2,NR-1
        IF(G(IR).GT.0.D0) THEN
          SVAR=G(IR)/(G(IR+1)-G(IR))
          SVAR=(1.D0-SVAR)*POT(IR)+SVAR*POT(IR-1)
          POT(IR:)=SVAR
          EXIT
        END IF
      ENDDO
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,V00,E,IRCL,IROUT)
!     **************************************************************************
!     **  1) ESTIMATE THE CLASSICAL TURNING POINT R(IRCL).                    **
!     **  2) ESTIMATE THE OUTERMOST GRID POINT R(IROUT) FOR INWARD INTEGRATION**
!     **     FROM A WKB SOLUTION OF THE SCHROEDINGER EQUATION. THE CRITERION  **
!     **     IROUT IS THAT THE SOLUTION GROWS FROM R(IROUT) TO R(IRCL)        **
!     **     BY A FACTOR OF LESS THAN XMAX,                                   **
!     **                                                                      **
!     **  REMARK: IT IS NOT POSSIBLE TO CHOOSE THE INNERMOST CLASSICAL        **
!     **    TURNING POINT, BECAUSE OTHERWISE WE MUST EXPLICITELY DESCRIBE     **
!     **    THE MATCHING OF THE WKB SOLUTIONS AT THE OUTER TURNING POINTS     **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID      ! GRID ID
      INTEGER(4),INTENT(IN) :: NR       ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: L        ! ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: XMAX     ! MAXIMUM TOLERABLE RATIO OF PHI 
      REAL(8)   ,INTENT(IN) :: V00(NR)  ! RADIAL SPHERICAL POTENTIAL
      REAL(8)   ,INTENT(IN) :: E        ! ENERGY
      INTEGER(4),INTENT(OUT):: IRCL     ! CLASSICAL TURNING POINT
      INTEGER(4),INTENT(OUT):: IROUT    ! OUTERMOST GRID POINT
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI) ! SPH. HARMON. FOR LM=0
      REAL(8)               :: R(NR)    ! RADIAL GRID
      REAL(8)               :: TKIN(NR)    ! RADIAL KINETIC ENERGY
      REAL(8)               :: FAC,SVAR,SUMVAL,XMAXLOG
      INTEGER(4)            :: IR
      LOGICAL               :: TCHK
      LOGICAL(4),PARAMETER  :: TINNER=.FALSE.
!     **************************************************************************
      XMAXLOG=LOG(XMAX)
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  DETERMINE THE RADIAL KINETIC ENERGY                                 ==
!     ==========================================================================
!     FAC=0.5D0*REAL(L*(L+1),KIND=8)
      FAC=0.D0   ! ZENTRIFUGAL POTENTIAL CAN LEAD TO VERY SMALL ROUT
      TKIN(2:)=E-V00(2:)*Y0-FAC/R(2:)**2  ! TKIN=E-VEFF
      TKIN(1)=TKIN(2)  ! AVOID DIVIDE-BY-ZERO FOR R(1)=0    
!
!     ==========================================================================
!     ==  IDENTIFY CLASSICAL TURNING POINT AS THE OUTERMOST POINT, WHERE      ==
!     ==  THE KINETIC ENERGY SWITCHES FROM POSITIVE TO NEGATIVE VALUES        ==
!     ==========================================================================
      IF(TINNER) THEN
        IRCL=NR
        DO IR=2,NR
          IF(TKIN(IR-1).GT.0.D0.AND.TKIN(IR).LE.0.D0) THEN
            IRCL=IR+1
            EXIT
          END IF
        ENDDO
        SUMVAL=0.D0
        DO IR=IRCL,NR
          SVAR=-TKIN(IR)
          IF(SVAR.LE.0.D0) CYCLE
          SUMVAL=SUMVAL+SQRT(2.D0*SVAR)*0.5D0*(R(IR+1)-R(IR-1))
          IROUT=IR-1
          IF(SUMVAL.GT.XMAXLOG) THEN
            TCHK=.TRUE.
            EXIT
          END IF
        ENDDO
      ELSE
        SUMVAL=0.D0
        IROUT=1
        IRCL=NR
        TCHK=.FALSE.
        DO IR=2,NR-1
          SVAR=-TKIN(IR)
          IF(SVAR.LT.0.D0) THEN   ! KINETIC ENERGY POSITIVE; DO NOTHING
            SVAR=0.D0
            SUMVAL=0.D0
            IRCL=IR+1        ! RCL WILL BE THE FIRST POINT WITH POSITIVE EKIN
            CYCLE
          END IF
!         == SUMVAL IS THE APPROXIMATE INTEGRAL OF THE MOMENTUM FROM THE =======
!         == CLASSICAL TURNING POINT OUTWARD ===================================
          SUMVAL=SUMVAL+SQRT(2.D0*SVAR)*0.5D0*(R(IR+1)-R(IR-1))
          IROUT=IR-1
          IF(SUMVAL.GT.XMAXLOG) THEN
            TCHK=.TRUE.
            EXIT
          END IF
        ENDDO
      END IF
!CALL ATOMLIB_WRITEPHI('TEST1.DAT',GID,NR,1,-TKIN)
!
!     ==========================================================================
!     == FIX UP END OF THE GRID                                               ==
!     ==========================================================================
      IF(.NOT.TCHK) THEN
        SVAR=-TKIN(NR)
        SUMVAL=SUMVAL+0.5D0*(R(NR)-R(NR-1))*SVAR
        IF(SUMVAL.GT.XMAXLOG) THEN
          IROUT=NR-1
        ELSE
!         == THE MAX VALUE OFR IROUT IS NR-1, BECAUSE IROUT IS USED TO        ==
!         == SET AN INHOMOGENEITY. THE SOLVER IS NOT SENSITIVE TO THE LAST    ==
!         == GRID POINT                                                       ==
          IROUT=NR-1  ! IROUT IS MAXIMUM AT NR-1
        END IF
      END IF
      IRCL=MIN(IRCL,IROUT-1)
      RETURN
      END
!
!     .......................................................................
      SUBROUTINE SCHROEDINGER_LSOLD(LMX,C)
!     **                                                                    **
!     **  SPIN ORBIT MATRIX ELEMENTS                                        **
!     **                                                                    **
!     **  MATRIX ELEMENTS OF L*SIGMA IN REAL SPHERICAL HARMONICS            **
!     **  WHERE SIGMA ARE THE PAULI MATRICES AND L ARE THE ANGULAR MOMENTA  **
!     **                                                                    **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX
      COMPLEX(8),INTENT(OUT):: C(LMX,2,LMX,2)
      REAL(8)               :: C1(LMX,2,LMX,2)
      INTEGER(4)            :: LX
      INTEGER(4)            :: L,M,M1,M2,LM,LM1P,LM1M,LM2P,LM2M,IS,LM1,LM2,IS1,IS2
      REAL(8)               :: SVAR
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)            :: FAC1,FAC2,FAC3,FAC4
      LOGICAL               :: TPRINT=.FALSE.
      LOGICAL               :: TTEST=.TRUE.
      REAL(8)               :: SQR2IN
      COMPLEX(8),ALLOCATABLE :: H(:,:),U(:,:)
      REAL(8)   ,ALLOCATABLE :: E(:)
!     ************************************************************************
      LX=INT(SQRT(REAL(LMX))-1.D0) ! LMX=(LX+1)**2
      IF((LX+1)**2.NE.LMX) THEN
        CALL ERROR$MSG('LMX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMX',LMX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('SCHROEDINGER_LSOLD')
      END IF
      SQR2IN=1.D0/SQRT(2.D0)
!
!     ========================================================================
!     ==  MATRIX ELEMENTS <L,M,S|\SIGMA*L|L,M,S> WHERE |L,M,S> ARE          ==
!     ==  SPHERICAL HARMONICS (ANGULAR MOMENTUM EIGENSTATES)                ==
!     ========================================================================
      C1=0.D0
      LM=0
      DO L=0,LX 
        DO M=-L,L
          LM=LM+1
          C1(LM,1,LM,1)= REAL(M)
          C1(LM,2,LM,2)=-REAL(M)
          IF(M.EQ.L) CYCLE
          SVAR=SQRT( REAL((L-M)*(L+M+1)) )
          C1(LM,1,LM+1,2)=SVAR
          C1(LM+1,2,LM,1)=SVAR
        ENDDO
      ENDDO
!
!     ==================================================================
!     == TRANSFORM TO REAL SPHERICAL HARMONICS                        ==
!     ==================================================================
      IF(TPRINT) THEN
        DO L=1,LX
          LM=L**2
          DO IS=1,2
            DO M=1,2*L+1
              WRITE(*,FMT='(20F8.3)')C1(LM+M,IS,LM+1:LM+2*L+1,:)
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==================================================================
!     == TRANSFORM TO REAL SPHERICAL HARMONICS                        ==
!     ==================================================================
      C=(0.D0,0.D0)
      DO L=1,LX    !NO S.O. VOUPLING IN S-CHANNEL
        LM1P=L**2+L+1
        LM1M=LM1P
!       == M1=M2=0
        C(LM1M,:,LM1M,:)=CMPLX(C1(LM1M,:,LM1M,:),KIND=8)
!       == (M1=0 AND M2.NEQ.0) OR (N1.NEQ.0 AND M2=0)
        LM2P=L**2+L+1
        LM2M=LM2P
        DO M2=1,L
          LM2P=LM2P+1
          LM2M=LM2M-1
          FAC1=CMPLX(SQR2IN,KIND=8)
          FAC2=CMPLX(SQR2IN*(-1.D0)**M2,KIND=8)
          C(LM1M,:,LM2P,:)=FAC1*C1(LM1M,:,LM2P,:)+FAC2*C1(LM1M,:,LM2M,:)
          C(LM2P,:,LM1M,:)=FAC1*C1(LM2P,:,LM1M,:)+FAC2*C1(LM2M,:,LM1M,:)
          FAC1=-CI*FAC1
          FAC2=CI*FAC2
          C(LM1M,:,LM2M,:)=+FAC1*C1(LM1M,:,LM2P,:)+FAC2*C1(LM1M,:,LM2M,:)
          C(LM2M,:,LM1M,:)=-FAC1*C1(LM2P,:,LM1M,:)-FAC2*C1(LM2M,:,LM1M,:)
        ENDDO
!       == M1.NEQ.0, M2.NEQ.0 ==============================================
        DO M1=1,L
          LM1P=LM1P+1
          LM1M=LM1M-1
          LM2P=L**2+L+1
          LM2M=LM2P
          DO M2=1,L
            LM2P=LM2P+1
            LM2M=LM2M-1
            FAC1=CMPLX(0.5D0,KIND=8)
            FAC2=CMPLX(0.5D0*(-1.D0)**(M1+M2),KIND=8)
            FAC3=CMPLX(0.5D0*(-1.D0)**M2,KIND=8)
            FAC4=CMPLX(0.5D0*(-1.D0)**M1,KIND=8)
            C(LM1P,:,LM2P,:)=+FAC1*C1(LM1P,:,LM2P,:)+FAC2*C1(LM1M,:,LM2M,:) &
    &                        +FAC3*C1(LM1P,:,LM2M,:)+FAC3*C1(LM1M,:,LM2P,:)
            C(LM1M,:,LM2M,:)=+FAC1*C1(LM1P,:,LM2P,:)+FAC2*C1(LM1M,:,LM2M,:) &
    &                        -FAC3*C1(LM1P,:,LM2M,:)-FAC3*C1(LM1M,:,LM2P,:)
            FAC1=FAC1*CI
            FAC2=FAC2*CI
            FAC3=FAC3*CI
            FAC4=FAC4*CI
            C(LM1P,:,LM2M,:)=-FAC1*C1(LM1P,:,LM2P,:)+FAC2*C1(LM1M,:,LM2M,:) &
    &                        +FAC3*C1(LM1P,:,LM2M,:)-FAC3*C1(LM1M,:,LM2P,:)
            C(LM1M,:,LM2P,:)=+FAC1*C1(LM1P,:,LM2P,:)-FAC2*C1(LM1M,:,LM2M,:) &
    &                        +FAC3*C1(LM1P,:,LM2M,:)-FAC3*C1(LM1M,:,LM2P,:)
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     == TRANSFORM TO REAL SPHERICAL HARMONICS                        ==
!     ==================================================================
      IF(TPRINT) THEN
        DO L=1,LX
          ALLOCATE(H(2*(2*L+1),2*(2*L+1)))
          ALLOCATE(U(2*(2*L+1),2*(2*L+1)))
          ALLOCATE(E(2*(2*L+1)))
          LM=L**2
          DO IS=1,2
            DO M=1,2*L+1
              WRITE(*,FMT='(20("(",F10.3,",",F10.3,")"))')C(LM+M,IS,LM+1:LM+2*L+1,:)
            ENDDO
          ENDDO
          H(:2*L+1,:2*L+1)=C(LM+1:LM+2*L+1,1,LM+1:LM+2*L+1,1)
          H(:2*L+1,2*L+2:)=C(LM+1:LM+2*L+1,1,LM+1:LM+2*L+1,2)
          H(2*L+2:,:2*L+1)=C(LM+1:LM+2*L+1,2,LM+1:LM+2*L+1,1)
          H(2*L+2:,2*L+2:)=C(LM+1:LM+2*L+1,2,LM+1:LM+2*L+1,2)
          CALL LIB$DIAGC8(2*(2*L+1),H,E,U)
          PRINT*,'E',L,E
          DEALLOCATE(H)
          DEALLOCATE(E)
          DEALLOCATE(U)
        ENDDO
      END IF
!
!     ==================================================================
!     == TRANSFORM TO REAL SPHERICAL HARMONICS                        ==
!     ==================================================================
      IF(TTEST) THEN
        DO LM1=1,LMX
          DO IS1=1,2
            DO LM2=1,LMX
              DO IS2=1,2
                IF(C(LM1,IS1,LM2,IS2).NE.CONJG(C(LM2,IS2,LM1,IS1))) THEN
                  CALL ERROR$MSG('SPIN ORBIT IS NOT HERMITEAN')
                  CALL ERROR$STOP('RADIAL$LS')
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END IF
      RETURN
      END SUBROUTINE SCHROEDINGER_LSOLD
!
!     .......................................................................
      SUBROUTINE SCHROEDINGER_LS(LMX,C)
!     **                                                                    **
!     **  SPIN ORBIT MATRIX ELEMENTS                                        **
!     **                                                                    **
!     **  MATRIX ELEMENTS OF L*SIGMA IN REAL SPHERICAL HARMONICS            **
!     **  WHERE SIGMA ARE THE PAULI MATRICES AND L ARE THE ANGULAR MOMENTA  **
!     **                                                                    **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX
      COMPLEX(8),INTENT(OUT):: C(LMX,2,LMX,2)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      LOGICAL               :: TTEST=.FALSE.
      COMPLEX(8)            :: LXMAT(LMX,LMX),LYMAT(LMX,LMX),LZMAT(LMX,LMX)
      REAL(8)               :: SVAR
      COMPLEX(8)            :: CTEST(LMX,2,LMX,2)
!     ************************************************************************
      CALL SPHERICAL$L(LMX,LXMAT,LYMAT,LZMAT)
!
!     ========================================================================
!     ==  MATRIX ELEMENTS <L,M,S|\SIGMA*L|L,M,S> WHERE |L,M,S> ARE          ==
!     ==  SPHERICAL HARMONICS (ANGULAR MOMENTUM EIGENSTATES)                ==
!     ========================================================================
      C(:,1,:,1)=LZMAT(:,:)
      C(:,1,:,2)=LXMAT(:,:)-CI*LYMAT(:,:)
      C(:,2,:,1)=LXMAT(:,:)+CI*LYMAT(:,:)
      C(:,2,:,2)=-LZMAT(:,:)
!
!     ==================================================================
!     == TRANSFORM TO REAL SPHERICAL HARMONICS                        ==
!     ==================================================================
      IF(TTEST) THEN
        SVAR=MAXVAL(ABS(C(:,1,:,1)-TRANSPOSE(CONJG(C(:,1,:,1)))))
        SVAR=MAX(SVAR,MAXVAL(ABS(C(:,1,:,2)-TRANSPOSE(CONJG(C(:,2,:,1))))))
        SVAR=MAX(SVAR,MAXVAL(ABS(C(:,2,:,2)-TRANSPOSE(CONJG(C(:,2,:,2))))))
        IF(SVAR.GT.1.D-10) THEN
          CALL ERROR$MSG('SPIN ORBIT IS NOT HERMITEAN')
          CALL ERROR$R8VAL('DEVIATION',SVAR)
          CALL ERROR$STOP('RADIAL$LS')
        END IF
        CALL SCHROEDINGER_LSOLD(LMX,CTEST)
        SVAR=MAXVAL(ABS(C-CTEST))
        IF(SVAR.GT.1.D-7) THEN
          CALL ERROR$MSG('SPIN ORBIT DOES NOT AGREE WITH PREVIOUS')
          CALL ERROR$R8VAL('DEVIATION',SVAR)
          CALL ERROR$STOP('SCHROEDINGER_LS')
        END IF
      END IF
      RETURN
      END SUBROUTINE SCHROEDINGER_LS
!
!     .......................................................................
      SUBROUTINE SCHROEDINGER_SP(LMX,CR,CA)
!     **                                                                   **
!     **  CALCULATES THE MATRICES REQUIRED TO EVALUATE THE SMALL COMPONENT **
!     **  (POSITRONS) FROM THE SOLUTION OF THE LARGE COMPONENT (ELECTRONS) **
!     **                                                                   **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX
      COMPLEX(8),INTENT(OUT):: CR(LMX,2,LMX,2)
      COMPLEX(8),INTENT(OUT):: CA(LMX,2,LMX,2)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)               :: X(LMX,LMX),Y(LMX,LMX),Z(LMX,LMX)
      COMPLEX(8)            :: LX(LMX,LMX),LY(LMX,LMX),LZ(LMX,LMX)
      COMPLEX(8)            :: RCROSSL(LMX,LMX,3)
!     ***********************************************************************
      CALL SPHERICAL$ER(LMX,X,Y,Z)
      CALL SPHERICAL$L(LMX,LX,LY,LZ)
!
!     =======================================================================
!     == DETERMINE THE RADIAL PART CR                                      ==
!     =======================================================================
      CR(:,1,:,1)=CMPLX(Z(:,:),0.D0,KIND=8)
      CR(:,1,:,2)=CMPLX(X(:,:),-Y(:,:),KIND=8)
      CR(:,2,:,1)=CMPLX(X(:,:),Y(:,:),KIND=8)
      CR(:,2,:,2)=CMPLX(-Z(:,:),0.D0,KIND=8)
      CR(:,:,:,:)=0.5D0*CR(:,:,:,:)
!
!     =======================================================================
!     == DETERMINE THE RADIAL PART CA                                      ==
!     =======================================================================
      RCROSSL(:,:,1)=MATMUL(Y,LZ)-MATMUL(Z,LY)
      RCROSSL(:,:,2)=MATMUL(Z,LX)-MATMUL(X,LZ)
      RCROSSL(:,:,3)=MATMUL(X,LY)-MATMUL(Y,LX)
      CA(:,1,:,1)=RCROSSL(:,:,3)
      CA(:,1,:,2)=RCROSSL(:,:,1)-CI*RCROSSL(:,:,2)
      CA(:,2,:,1)=RCROSSL(:,:,1)+CI*RCROSSL(:,:,2)
      CA(:,2,:,2)=-RCROSSL(:,:,3)
      CA(:,:,:,:)=-0.5D0*CA(:,:,:,:)
      RETURN
      END SUBROUTINE SCHROEDINGER_SP
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,L,ISO,DREL,G,PHI,SPHI)
!     **************************************************************************
!     **  DETERMINES THE RADIAL PART OF THE SMALL COMPONENT OF THE SPHERICAL  **
!     **  DIRAC EQUATION.                                                     **
!     **                                                                      **
!     **  ISO=+1 FOR     PARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM           **
!     **  ISO=-1 FOR ANTIPARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM           **
!     **  ISO= 0 IN THE ABSENCE OF SPIN-ORBIT COUPLING                        **
!     **                                                                      **
!     **  G IS THE INHOMOGENEITY IN THE SMALL COMPONENTS OF THE               **
!     **    FOUR-COMPONENT DIRAC EQUATION. FOR THE NODELESS EQUATIONS IT IS   **
!     **    EQUAL TO THE SMALL COMPONENT OF THE NODELESS WAVE FUNCTION OF THE **
!     **    NEXT LOWER BAND.                                                  **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: L
      INTEGER(4),INTENT(IN) :: ISO         ! =1 FOR LS>0; =-1 FOR LS<0
      REAL(8)   ,INTENT(IN) :: DREL(NR)
      REAL(8)   ,INTENT(IN) :: G(NR)
      REAL(8)   ,INTENT(IN) :: PHI(NR)     ! LARGE COMPONENT
      REAL(8)   ,INTENT(OUT):: SPHI(NR)    ! SMALL COMPONENT
      REAL(8)               :: LAMBDA      ! =L FOR LS>0; -L-1 FOR LS<0
      REAL(8)               :: KAPPA       ! =-L-1 FOR LS>0; L FOR LS<0
      REAL(8)               :: SGNLAMBDA   !SGN(LAMBDA)
      REAL(8)               :: SPEEDOFLIGHT! SPEED OF LIGHT
      REAL(8)               :: R(NR)       ! RADIAL GRID
      LOGICAL(4),PARAMETER  :: TON=.TRUE.  !SWITCHES SMALL COMPONENT ON AND OFF
!     **************************************************************************
      IF(.NOT.TON) THEN
        SPHI=0.D0
        RETURN
      END IF
!
      IF(ISO.EQ.1) THEN   !PARALLEL SPIN AND ORBIT
        LAMBDA=REAL(L,KIND=8)   ! LAMBDA=-(KAPPA+1)
        KAPPA=REAL(-L-1,KIND=8)
        SGNLAMBDA=1.D0          ! SGNLAMBDA=-SGNKAPPA
      ELSE IF(ISO.EQ.-1) THEN !ANTI PARALLEL SPIN AND ORBIT
        LAMBDA=REAL(-L-1,KIND=8)
        KAPPA=REAL(L,KIND=8)
        SGNLAMBDA=-1.D0
      ELSE IF(ISO.EQ.0) THEN ! SCALAR RELATIVISTIC
        KAPPA=REAL(-1,KIND=8)
        LAMBDA=0.D0
        SGNLAMBDA=1.D0
      ELSE
         CALL ERROR$MSG('ILLEGAL VALUE OF ISO (MUST BE 1,0, OR -1)')
         CALL ERROR$STOP('SCHROEDINGER$SPHSMALLCOMPONENT')
      END IF
      CALL CONSTANTS$GET('C',SPEEDOFLIGHT)
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$DERIVE(GID,NR,PHI,SPHI)
      SPHI(2:)=SPHI(2:)+(1.D0+KAPPA)*PHI(2:)/R(2:)  !LAMBDA=-(1+KAPPA)
      SPHI(1)=SPHI(2)
!      SPHI=SGNLAMBDA*SPHI-G/SPEEDOFLIGHT
!      the sgnlambda does not occur in my nodes nor in those of robert schade
      SPHI=SPHI-G/SPEEDOFLIGHT
      SPHI=0.5D0*(1.D0+DREL)/SPEEDOFLIGHT*SPHI
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER_SMALLCOMPONENT(GID,NR,LMX,NPHI,DREL,PHI,SPHI)
!     **                                                                      **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMX
      INTEGER(4),INTENT(IN) :: NPHI
      REAL(8)   ,INTENT(IN) :: DREL(NR)
      COMPLEX(8),INTENT(IN) :: PHI(NR,LMX,2,NPHI)
      COMPLEX(8),INTENT(OUT):: SPHI(NR,LMX,2,NPHI)
      COMPLEX(8)            :: CR(LMX,2,LMX,2)
      COMPLEX(8)            :: CA(LMX,2,LMX,2)
      REAL(8)               :: A(NR)
      REAL(8)               :: AUX1(NR),AUX2(NR),SVAR
      COMPLEX(8)            :: DPHI(NR,LMX,2)
      COMPLEX(8)            :: CSVAR
      INTEGER(4)            :: IB,IS,LM,IS1,IS2,LM1,LM2
!     ************************************************************************
      CALL SCHROEDINGER_SP(LMX,CR,CA)
      CALL CONSTANTS$GET('C',SVAR)
      A(:)=(1.D0+DREL(:))/(2.D0*SVAR)
      SPHI(:,:,:,:)=(0.D0,0.D0)
      DO IB=1,NPHI
        DO IS=1,2
          DO LM=1,LMX
            CALL RADIAL$DERIVE(GID,NR,REAL(PHI(:,LM,IS,IB)),AUX1)
            CALL RADIAL$DERIVE(GID,NR,AIMAG(PHI(:,LM,IS,IB)),AUX2)
            DPHI(:,LM,IS)=CMPLX(AUX1,AUX2,KIND=8)
          ENDDO
        ENDDO
!
        DO IS2=1,2
          DO LM2=1,LMX
            DO IS1=1,2
              DO LM1=1,LMX
                CSVAR=CR(LM1,IS1,LM2,IS2)
                IF(CSVAR.NE.0.D0) THEN
                  SPHI(:,LM1,IS1,IB)=SPHI(:,LM1,IS1,IB)+CSVAR*DPHI(:,LM2,IS2)
                END IF
                CSVAR=CA(LM1,IS1,LM2,IS2)
                IF(CSVAR.NE.0.D0) THEN
                  SPHI(:,LM1,IS1,IB)=SPHI(:,LM1,IS1,IB)+CSVAR*PHI(:,LM2,IS2,IB)
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!
        DO LM=1,LMX
          DO IS=1,2
            SPHI(:,LM,IS,IB)=A(:)*SPHI(:,LM,IS,IB)                  
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ....................................................................
      SUBROUTINE SCHROEDINGER_WRITEPHI(GID,NR,FILE,NF,NPHI,IRC,PHIL,PHIR)
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: NR
      INTEGER(4)  ,INTENT(IN) :: NF
      INTEGER(4)  ,INTENT(IN) :: IRC
      INTEGER(4)  ,INTENT(IN) :: NPHI
      REAL(8)     ,INTENT(IN) :: PHIL(NR,NF,NPHI)
      REAL(8)     ,INTENT(IN) :: PHIR(NR,NF,NPHI)
      INTEGER(4)              :: IR,IPHI
      CHARACTER(32)           :: FILE1
      REAL(8)                 :: R(NR)
!     **********************************************************************       
      CALL RADIAL$R(GID,NR,R)
      DO IPHI=1,NPHI
        WRITE(FILE1,FMT='(I1)')IPHI
        FILE1=TRIM(FILE)//'_'//TRIM(FILE1)//'.DAT'
        OPEN(100,FILE=FILE1)
        DO IR=1,IRC
          WRITE(100,*)R(IR),PHIL(IR,:,IPHI)
        ENDDO
        DO IR=IRC+1,NR
          WRITE(100,*)R(IR),PHIR(IR,:,IPHI)
        ENDDO
        CLOSE(100)
      ENDDO
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER$DREL(GID,NR,POT,E,DREL)
!     **                                                                      **
!     **  DREL IS A MEASURE OF THE RELATIVISTIC CORRECTIONS                   **
!     **     D:=1/MREL-1/M0                                                   **
!     **  WHERE MREL IS THE RELATIVISTIC MASS  MREL=M0+(E-POT)/(2C**2)        **
!     **  AND  M0 IS THE REST MASS                                            **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **  -  RELATIVISTIC CORRECTIONS FOR EKIN<0 ARE SWITCHED OFF!            **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID         ! GRID ID
      INTEGER(4),INTENT(IN) :: NR          ! #(RADIAL GRID POINTS
      REAL(8)   ,INTENT(IN) :: POT(NR)     ! POTENTIAL (MULTIPLY WITH Y0!)
      REAL(8)   ,INTENT(IN) :: E           ! ONE-PARTICLE ENERGY
      REAL(8)   ,INTENT(OUT):: DREL(NR)    ! RELATIVISTIC CORRECTION 
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI) ! SPH. HARMON. FOR LM=0
      REAL(8)               :: C           ! SPEED OF LIGHT
      REAL(8)               :: EKIN(NR)
!     **************************************************************************
      CALL CONSTANTS$GET('C',C)
      EKIN(:)=E-POT(:)*Y0
      EKIN(:)=MAX(EKIN,-C**2)   ! FIX TO AVOID A DIVIDE BY ZERO.
      DREL(:)=-EKIN(:)/(EKIN(:)+2.D0*C**2)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RMIN,RC,PHASE)
!     **************************************************************************
!     **  CALCULATES THE PHASE SHIFT FOR A RADIAL FUNCTION AT RADIUS RC       **
!     **                                                                      **
!     **  THE PHASE SHIFT IS DEFINED AS                                       **
!     **     0.5-1/PI * ATAN (DPHIDR/PHI)+NN                                  **
!     **  WHERE PHI AND DPHIDR ARE VALUE AND DERIVATIVE OF PHI AT RADIUS RC   **
!     **  AND NN IS THE NUMBER OF NODES INSIDE RC.                            **
!     **                                                                      **
!     **  THIS DEFINITION OF THE PHASE SHIFT DIFFERS FROM THE TERM USED       **
!     **  IN SCATTERING THEORY                                                **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: PHI(NR)
      REAL(8)   ,INTENT(IN) :: RMIN  ! NODES INSIDE RMIN ARE NOT COUNTED
      REAL(8)   ,INTENT(IN) :: RC    ! PHASESHIFT CALCULATED AT RC
      REAL(8)   ,INTENT(OUT):: PHASE
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: R(NR)
      REAL(8)               :: VAL,DER
      INTEGER(4)            :: IR
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     == THIS IS A MORE ROBUST EXPRESSION WITH LINEAR INTERPOLATION BETWEEN
!     == TWO POINTS. IT MAY BE BETTER FOR STRONG DEVIATIONS FROM POLYNOMIAL
!     == CHARACTER OF THE FUNCTION.
      IF(RC.GT.R(NR)) THEN
        CALL ERROR$MSG('RC > R(NR)')
        CALL ERROR$STOP('SCHROEDINGER$PHASESHIFT')
      END IF
      PHASE=0.D0
      VAL=0.D0
      DO IR=3,NR-1  ! DO NOT CONSIDER A ZERO AT THE ORIGIN
        IF(R(IR-1).LT.RMIN) CYCLE  ! DO NOT COUNT ZEROS WITHIN RMIN
        IF(R(IR).GE.RC) THEN
          DER=(PHI(IR)-PHI(IR-1))/(R(IR)-R(IR-1))
          VAL=PHI(IR-1)+DER*(RC-R(IR-1))
          IF(VAL*PHI(IR-1).LT.0.D0) PHASE=PHASE+1.D0
          EXIT
        END IF
        IF(PHI(IR)*PHI(IR-1).LT.0.D0) PHASE=PHASE+1.D0
        IF(PHI(IR).EQ.0.D0) THEN
          IF(PHI(IR+1)*PHI(IR-1).LT.0.D0)PHASE=PHASE+1.D0
        END IF
      ENDDO
      IF(VAL.NE.0.D0) THEN
         PHASE=0.5D0-ATAN(DER/VAL)/PI+PHASE
      ELSE
         PHASE=1.D0+PHASE
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER$RESOLVEPHASESHIFT(PHASE,VAL,DER,NN)
!     **************************************************************************
!     **                                                                      **
!     **  THE PHASE SHIFT IS DEFINED AS                                       **
!     **     0.5-1/PI * ATAN (DPHIDR/PHI)+NN                                  **
!     **  WHERE PHI AND DPHIDR ARE VALUE AND DERIVATIVE OF PHI AT RADIUS RC   **
!     **  AND NN IS THE NUMBER OF NODES INSIDE RC.                            **
!     **                                                                      **
!     **  THIS DEFINITION OF THE PHASE SHIFT DIFFERS FROM THE TERM USED       **
!     **  IN SCATTERING THEORY                                                **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: NN
      REAL(8)   ,INTENT(IN)  :: PHASE
      REAL(8)   ,INTENT(OUT) :: VAL
      REAL(8)   ,INTENT(OUT) :: DER
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)                :: SVAR
!     **************************************************************************
      SVAR=MODULO(PHASE,1.D0)
      NN=NINT(PHASE-SVAR)
      IF(ABS(SVAR-0.5D0).LT.0.25D0) THEN
        VAL=1.D0
        DER=-TAN(PI*(SVAR+0.5D0))
      ELSE
!       == COT(X)=-TAN(X+PI/2)=-TAN(X-PI/2)=====================================
        VAL=TAN(PI*SVAR)
        DER=1.D0
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER$TEST()
!     **************************************************************************
!     ** TEST ROUTINE FOR SCHROEDINGER OBJECT                                 **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)           :: GID
      INTEGER(4),PARAMETER :: NR=250
      INTEGER(4),PARAMETER :: LMRX=9   !#(ANGULAR COMPONENTS OF THE DENSITY)
      INTEGER(4),PARAMETER :: LMX=9    !#(ANGULAR MOMENTUM COMPONENTS OF PSI)
      INTEGER(4),PARAMETER :: NB=3
      INTEGER(4),PARAMETER :: NBG=5
      REAL(8)   ,PARAMETER :: AEZ=7.D0
      REAL(8)   ,PARAMETER :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER :: Y0=1.D0/SQRT(4.D0*PI) ! SPH. HARMON. FOR LM=0
      INTEGER(4)           :: LOFI(NB)
      INTEGER(4)           :: NN(NB)
      REAL(8)              :: F(NB)  ! OCCUPATION PER SHELL
      REAL(8)              :: E(NB)  
      REAL(8)              :: EBG(NBG)
      REAL(8)              :: PHI(NR,LMX,NBG)
      REAL(8)              :: TPHI(NR,LMX,NBG)
      REAL(8)              :: R(NR)
      REAL(8)              :: POT(NR,LMRX)
!     **************************************************************************
!
!     == DEFINE RADIAL GRID
      CALL RADIAL__NEW('LOG',GID)
      CALL RADIAL__SETR8(GID,'R1',1.056D-4)
      CALL RADIAL__SETR8(GID,'DEX',0.05D0)
      CALL RADIAL__SETI4(GID,'NR',NR)
!
!     == DEFINE HYDROGEN-LIKE POTENTIAL ========================================
      CALL RADIAL__R(GID,NR,R)      
      POT=0.D0
      POT(:,1)=-AEZ/R(:)/Y0
!
!     ==========================================================================
!     == TEST SCHROEDINGER$SPHERICAL                                          ==
!     ==========================================================================
!      CALL SCHROEDINGER_TESTSPHERICAL(GID,NR,AEZ)
!
!     ==========================================================================
!     == TEST SCHROEDINGER$LBND_SCALREL                                   ==
!     ==========================================================================
!     CALL SCHROEDINGER_TESTLBND_SCALREL(GID,NR,AEZ)
!
!     ==========================================================================
!     == TEST RADIAL$NONSPHBOUND_NSO_SLOC_OV                                  ==
!     ==========================================================================
      CALL SCHROEDINGER_TESTLBND_SLOC(GID,NR,AEZ)
STOP
!
!     ==  
      LOFI(:)=(/0,0,1/)
      NN(:)=(/0,1,0/)
      F(:)=(/1.D0,0.D0,0.D0/)
      E(:)=(/-0.5D0,0.125D0,-0.125D0/)
      E(:)=0.D0

!
!      CALL SF_CORESTATES_NSPH(GID,NR,LMRX,NB,LOFI,F,NN &
!     &                       ,POT,E,LMX,NBG,EBG,PHI,TPHI)

      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER_TESTLBND_SLOC(GID,NR,AEZ)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: AEZ
      INTEGER(4),PARAMETER  :: LMRX=9
      INTEGER(4),PARAMETER  :: L=1 !ANGULAR MOMENTUM USED IN THE TEST
      INTEGER(4),PARAMETER  :: N=2 !MAIN QUANTUM NUMBER USED IN THE TEST
      INTEGER(4),PARAMETER  :: LMX=(L+1)**2
      INTEGER(4),PARAMETER  :: NPHI=2*L+1
      INTEGER(4),PARAMETER  :: LXSL=2
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI) ! SPH. HARMON. FOR L=0
      REAL(8)               :: R(NR)    ! RADIAL GRID
      REAL(8)               :: DREL(NR) ! ARRAY FOR RELATIVISTIC CORRECTION
      REAL(8)               :: POT(NR,LMRX)
      REAL(8)               :: G(NR,LMX)
      REAL(8)               :: EB(NPHI)
      REAL(8)               :: ENU
      REAL(8)               :: POTSL(NR,LXSL)
      REAL(8)               :: OVSL(NR,LXSL)
      REAL(8)               :: PHI(NR,LMX,NPHI)
      REAL(8)               :: TPHI(NR,LMX,NPHI)
      LOGICAL(4)            :: TOK
!     **************************************************************************
      IF(L+1.GT.N) THEN
        CALL ERROR$MSG('L+1 MUST BE SMALLER THAN N')
        CALL ERROR$STOP('SCHROEDINGER_TESTLBND_SCALREL')
      END IF
      CALL RADIAL__R(GID,NR,R)      
!
!     ==========================================================================
!     == SET UP HYDROGEN-LIKE POTENTIAL AND BOUND STATE ENERGY                ==
!     ==========================================================================
      POT(:,:)=0.D0
      POT(:,1)=-AEZ/R(:)/Y0
      POT(:,2)=1.D-1*R(:)*EXP(-0.1D0*R(:)**2)
      POTSL=0.D0
      OVSL=0.D0
      DREL(:)=0.D0
      G(:,:)=0.D0
      ENU=-0.5D0*AEZ**2/REAL(N**2,KIND=8)   !GROUND STATE
CALL ERROR$MSG('THE FOLLOWING CALL IS INCONSISTENT WITH THE ROUTINE CALLED')
CALL ERROR$MSG('FIX THE BUG BEFORE CONTINUING')
CALL ERROR$MSG('SCHROEDINGER_TESTLBND_SLOC')
! 
!      CALL RADIAL$NONSPHBOUND_NSO_SLOC_OV(GID,NR,LMX,LMRX,POT,LXSL,POTSL,OVSL,G,ENU &
!     &                             ,NPHI,EB,PHI,TPHI,TOK)

      PRINT*,'EB ',EB
      CALL SCHROEDINGER_WRITEPHI(GID,NR,'NSOSLOCOV',LMX,NPHI,100,PHI,PHI)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER_TESTLBND_SCALREL(GID,NR,AEZ)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: AEZ
      INTEGER(4),PARAMETER  :: LMRX=9
      INTEGER(4),PARAMETER  :: L=1 !ANGULAR MOMENTUM USED IN THE TEST
      INTEGER(4),PARAMETER  :: N=2 !MAIN QUANTUM NUMBER USED IN THE TEST
      INTEGER(4),PARAMETER  :: LMX=(L+1)**2
      INTEGER(4),PARAMETER  :: NPHI=2*L+1
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI) ! SPH. HARMON. FOR L=0
      REAL(8)               :: R(NR)    ! RADIAL GRID
      REAL(8)               :: DREL(NR) ! ARRAY FOR RELATIVISTIC CORRECTION
      REAL(8)               :: POT(NR,LMRX)
      REAL(8)               :: G(NR,LMX)
      REAL(8)               :: EB(NPHI)
      REAL(8)               :: ENU
      REAL(8)               :: PHI(NR,LMX,NPHI)
      REAL(8)               :: TPHI(NR,LMX,NPHI)
      LOGICAL(4)            :: TOK
      LOGICAL(4)            :: TMAINSH(L+1)
      INTEGER(4)            :: LX=L
!     **************************************************************************
      IF(L+1.GT.N) THEN
        CALL ERROR$MSG('L+1 MUST BE SMALLER THAN N')
        CALL ERROR$STOP('SCHROEDINGER_TESTLBND_SCALREL')
      END IF
      CALL RADIAL$R(GID,NR,R)      
!
!     ==========================================================================
!     == SET UP HYDROGEN-LIKE POTENTIAL WITH AN ELECTRIC FIELD                ==
!     ==========================================================================
      POT(:,:)=0.D0
      POT(:,1)=-AEZ/R(:)/Y0
      POT(:,2)=1.D-1*R(:)*EXP(-0.1D0*R(:)**2)
      DREL(:)=0.D0
      ENU=-0.5D0*AEZ**2/REAL(N**2,KIND=8)   !GROUND STATE
      G(:,:)=0.D0
!
!     ==========================================================================
!     == DETERMINE NONPSHERICAL BOUND STATES                                  ==
!     ==========================================================================
      LX=1
      TMAINSH(1)=.FALSE.
      TMAINSH(2)=.TRUE.
      CALL SCHROEDINGER$LBND_SCALREL(GID,NR,LMX,LX,LMRX,TMAINSH,POT,DREL,G,ENU &
     &                                  ,NPHI,EB,PHI,TPHI,TOK)
!
!     ==========================================================================
!     == WRITE RESULT                                                         ==
!     ==========================================================================
      PRINT*,'EB ',EB
      CALL SCHROEDINGER_WRITEPHI(GID,NR,'LBND_SCALREL',LMX,NPHI,100,PHI,PHI)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER_TESTSPHERICAL(GID,NR,AEZ)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: AEZ
      INTEGER(4),PARAMETER  :: L=2   !ANGULAR MOMENTUM USED IN THE TEST
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI) ! SPH. HARMON. FOR L=0
      INTEGER(4)            :: N=4   !MAIN QUANTUM NUMBER USED IN THE TEST
      INTEGER(4)            :: NFIL=70
      INTEGER(4)            :: SO
      REAL(8)               :: E
      REAL(8)               :: R(NR)    ! RADIAL GRID
      REAL(8)               :: POT(NR)  ! POTENTIAL
      REAL(8)               :: G(NR)    ! INHOMOGENEITY
      REAL(8)               :: DREL(NR) ! ARRAY FOR RELATIVISTIC CORRECTION
      INTEGER(4)            :: IDIR     ! DIRECTION OF INTEGRATION
      REAL(8)               :: PHIL(NR) ! WAVE FUNCTION FROM OUTWARD INTEGRATION
      REAL(8)               :: PHIR(NR) ! WAVE FUNCTION FROM INWARD INTEGRATION
      REAL(8)               :: PHI(NR)  ! FINAL WAVE FUNCTION
      REAL(8)               :: AUX(NR)  ! AUXILIARY ARRAY 
      REAL(8)               :: XMAX=1.D+20 ! MAX RANGE FOR WAVE FUNCTION 
      INTEGER(4)            :: IRCL     ! CLASSICAL RETURN RADIUS, MAX RADIUS 
      INTEGER(4)            :: IROUT    ! MAX RADIUS WITHOUT CREATING AN OVERFLOW
      INTEGER(4)            :: IR       ! 
      REAL(8)               :: SVAR     ! AUXILIARY VARIABLE
!     **************************************************************************
      IF(L+1.GT.N) THEN
        CALL ERROR$MSG('L+1 MUST BE SMALLER THAN N')
        CALL ERROR$STOP('SCHROEDINGER_TESTSPHERICAL')
      END IF
      CALL RADIAL__R(GID,NR,R)      
!
!     ==========================================================================
!     == SET UP HYDROGEN-LIKE POTENTIAL AND BOUND STATE ENERGY                ==
!     ==========================================================================
      POT(:)=-AEZ/R(:)/Y0
      DREL(:)=0.D0
      E=-0.5D0*AEZ**2/REAL(N**2,KIND=8)   !GROUND STATE
!
!     ==========================================================================
!     == OBTAIN WAVE FUNCTION BY INTEGRATING OUTWARD AND INWARD               ==
!     ==========================================================================
      SO=0
      CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT(:),E,IRCL,IROUT)
      IDIR=1
      G(:)=0.D0
      CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,E,IDIR,PHIL)
      PHIL(IROUT+1:)=0.D0
      IDIR=-1
      G(:)=0.D0
PRINT*,'IROUT',IROUT
      G(IROUT-1)=1.D-10
      CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,E,IDIR,PHIR)
      PHIR(IROUT+1:)=0.D0
!
!     ==========================================================================
!     == GLUE WAVE FUNCTIONS TOGETHER AND NORMALIZE ALL THREE FUNCTIONS       ==
!     ==========================================================================
      PHIR=PHIR/PHIR(IRCL)
      PHIL=PHIL/PHIL(IRCL)
      PHI(1:IRCL)=PHIL(1:IRCL)
      PHI(IRCL+1:)=PHIR(IRCL+1:)
      AUX(:)=PHI(:)**2*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
PRINT*,'SVAR',SVAR,IROUT,IRCL,PHIR(IRCL-1:IRCL+1)
      PHI=PHI/SQRT(SVAR)
      PHIL=PHIL/SQRT(SVAR)
      PHIR=PHIR/SQRT(SVAR)
!
!     ==========================================================================
!     == WRITE WAVE FUNCTIONS TO FILE                                         ==
!     ==========================================================================
      OPEN(NFIL,FILE='TESTSPHERICAL.DAT')
      DO IR=1,NR
        WRITE(NFIL,*)R(IR),PHI(IR),PHIL(IR),PHIR(IR)
      ENDDO
      CLOSE(NFIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCHROEDINGER$HYDROGENICDIRAC(GID,NR,AEZ,NB,LOFI,NNOFI,SOFI &
     &                                       ,EOFI,PHI,SPHI)
!     **************************************************************************
!     **  DETERMINE ANALYTIC WAVE FUNCTIONS OF THE DIRAC EQUATION             **
!     **  FOR A HYDROGENIC ATOM, I.E. V=-Z/R                                  **
!     **  FROM W. GORDON, Z. PHYS. 48, P11 (1928)                             **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID        ! GRID ID
      INTEGER(4),INTENT(IN) :: NR         ! #(RADIAL GRID POINTS)
      REAL(8)   ,INTENT(IN) :: AEZ        ! ATOMIC NUMBER
      INTEGER(4),INTENT(IN) :: NB         ! #(STATES)
      INTEGER(4),INTENT(IN) :: LOFI(NB)   ! ANGULAR MOMENTA
      INTEGER(4),INTENT(IN) :: NNOFI(NB)  ! #(NODES)
      INTEGER(4),INTENT(IN) :: SOFI(NB)   ! SIGN OF L*S
      REAL(8)   ,INTENT(OUT):: EOFI(NB)   ! ENERGY EIGENVALUES W/O REST ENERGY
      REAL(8)   ,INTENT(OUT):: PHI(NR,NB) ! LARGE COMPONENT
      REAL(8)   ,INTENT(OUT):: SPHI(NR,NB) ! SMALL COMPONENT
      REAL(8)               :: C          ! SPEED OF LIGHT
      REAL(8)               :: ALPHA      ! FINE STRUCTURE CONSTANT
      INTEGER(4)            :: N,L,ISO,IB,IR
      REAL(8)               :: JPRIME
      REAL(8)               :: JPHALF
      REAL(8)               :: NPRIME
      REAL(8)               :: SVAR
      REAL(8)               :: E          ! ENERGY INCLUDING REST ENERGY
      REAL(8)               :: K0
      REAL(8)               :: RHO
      REAL(8)               :: C01BYC02
      REAL(8)               :: R(NR)      ! RADIAL GRID
      REAL(8)               :: SIGMA1(NR)
      REAL(8)               :: SIGMA2(NR)
!     **************************************************************************
      CALL CONSTANTS('C',C)
      ALPHA=1.D0/C
      CALL RADIAL$R(GID,NR,R)
!
      DO IB=1,NB
        L=LOFI(IB)
        ISO=SOFI(IB)
        N=NNOFI(IB)+LOFI(IB)+1
        IF(ISO.EQ.1) THEN
          JPRIME=REAL(-L-1)
          JPHALF=REAL(L+1,KIND=8)
        ELSE
          JPRIME=REAL(L)
          JPHALF=REAL(L,KIND=8)
        END IF
        NPRIME=REAL(N)-ABS(JPRIME)  ! BEFORE EQ.2 OF GORDON
!       == EQ 6 IN GORDON ======================================================
        RHO=SQRT(JPRIME**2-(AEZ*ALPHA)**2)
!       == ENERGY INCLUDING REST ENERGY (LAST EQUATION IN GORDON) ==============
        E=C**2/SQRT(1.D0+(AEZ*ALPHA/(NPRIME+RHO))**2)
!       == EQ.4 IN GORDON ======================================================
        K0=C*SQRT(1.D0-(E/C**2)**2)            
!       == EQ. 8 IN GORDON =====================================================
        C01BYC02=-NPRIME/(ALPHA/SQRT(1-(E/C**2)**2)-JPRIME)
!       == EQ.10 IN GORDON =====================================================
        DO IR=1,NR
          CALL FHYPER(1.D0-NPRIME,2.D0*RHO+1.D0,2.D0*K0*R(IR),SIGMA1(IR))
          CALL FHYPER(-NPRIME,2.D0*RHO+1.D0,2.D0*K0*R(IR),SIGMA2(IR))
        ENDDO
        SIGMA1(:)=C01BYC02*EXP(-K0*R(:))*R(:)**(RHO-1.D0)*SIGMA1(:)
        SIGMA2(:)=EXP(-K0*R(:))*R(:)**(RHO-1.D0)*SIGMA2(:)
!       ==
        SPHI(:,IB)=SQRT(1.D0-E/C**2)*(SIGMA1(:)-SIGMA2(:))
        PHI(:,IB)=SQRT(1.D0+E/C**2)*(SIGMA1(:)+SIGMA2(:))
!       == NORMALIZE ===========================================================
        CALL RADIAL$INTEGRAL(GID,NR,R**2*(PHI(:,IB)**2+PHI(:,IB)**2),SVAR)
        SVAR=1.D0/SQRT(SVAR)
        PHI(:,IB)=PHI(:,IB)*SVAR
        SPHI(:,IB)=SPHI(:,IB)*SVAR
!       == SUBTRACT REST ENERGY ================================================
        EOFI(IB)=E-C**2
      ENDDO
      RETURN
      CONTAINS
!       .1.........2.........3.........4.........5.........6.........7.........8
        SUBROUTINE FHYPER(ALPHA,BETA,X,F)
!       ************************************************************************
!       ************************************************************************
        IMPLICIT NONE
        REAL(8)   ,INTENT(IN) :: ALPHA 
        REAL(8)   ,INTENT(IN) :: BETA 
        REAL(8)   ,INTENT(IN) :: X
        REAL(8)   ,INTENT(OUT):: F
        REAL(8)               :: FAC,RI
        INTEGER(4)            :: I
!       ************************************************************************
        F=1.D0
        FAC=1.D0
        DO I=1,10
          RI=REAL(I,KIND=8)
          FAC=FAC*(ALPHA+RI-1.D0)/(BETA+RI-1.D0)/RI*X
          F=F+FAC
        ENDDO
        RETURN
        END SUBROUTINE FHYPER
      END SUBROUTINE SCHROEDINGER$HYDROGENICDIRAC
