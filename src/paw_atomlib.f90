!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$AESCF(GID,NR,KEY,RBOX,AEZ,NX,NB,LOFI,SO,F,NN &
     &                        ,ETOT,POT,EOFI,PHI,SPHI)
!     **************************************************************************
!     ** MAKES A SELF-CONSISTENT CALCULATION OF AN ATOM IN A BOX WITH         **
!     ** RADIUS RBOX (RADIUS IS LIMITED BY THE GRID RBOX<R(NR-3) )            **
!     **                                                                      **
!     ** KEY='START,REL,SO,FOCK='                                             **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID       ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR        ! #(GRID POINTS)
      CHARACTER(*),INTENT(IN)    :: KEY     
      REAL(8)    ,INTENT(INOUT)  :: RBOX      ! BOX RADIUS     
      REAL(8)    ,INTENT(IN)     :: AEZ       ! ATOMIC NUMBER
      INTEGER(4) ,INTENT(IN)     :: NX        ! X#(STATES)
      INTEGER(4) ,INTENT(OUT)    :: NB        ! #(STATES)
      INTEGER(4) ,INTENT(INOUT)  :: LOFI(NX)  ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(INOUT)  :: SO(NX)    ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(INOUT)  :: F(NX)     ! OCCUPATION
      INTEGER(4) ,INTENT(INOUT)  :: NN(NX)    ! #(NODES)
      REAL(8)    ,INTENT(OUT)    :: ETOT      ! TOTAL ENERGY
      REAL(8)    ,INTENT(INOUT)  :: POT(NR)   ! POTENTIAL
      REAL(8)    ,INTENT(OUT)    :: EOFI(NX)  ! ONE-PARTICLE ENERGIES
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,NX)! ONE-PARTICLE WAVE FUNCTIONS
      REAL(8)    ,INTENT(OUT)    :: SPHI(NR,NX) ! SMALL COMPONENT
      REAL(8)                    :: R(NR)
      REAL(8)                    :: WGHT(NR)
      REAL(8)                    :: DREL(NR)  ! RELATIVISTIC CORRECTION
      REAL(8)                    :: RHO(NR)
      REAL(8)                    :: G(NR)     ! INHOMOGENEITY
      REAL(8)                    :: AUX(NR),AUX1(NR)   !AUXILIARY ARRAY
      REAL(8)                    :: MUX(NR)   !EXCHANGE ONLY POTENTIAL
      INTEGER(4)                 :: ITER
      INTEGER(4),PARAMETER       :: NITER=1000
      REAL(8)                    :: XAV,XMAX
      LOGICAL(4)                 :: CONVG
      REAL(8)   ,PARAMETER       :: TOL=1.D-6
      REAL(8)                    :: EKIN,EH,EXC
      LOGICAL(4),PARAMETER       :: TBROYDEN=.TRUE.
      REAL(8)                    :: POTIN(NR)
      REAL(8)                    :: SVAR
      INTEGER(4)                 :: IRBOX
      INTEGER(4)                 :: IR,I,IB,JB,ISO,L
      INTEGER(4)                 :: ISVAR,IARR(1)
      CHARACTER(32)              :: ID
      LOGICAL(4)                 :: TSTART  ! CALCULATE ANGULAR MOMENTA ETC
      LOGICAL(4)                 :: TREL    ! RELATIOVISTIC CALCULATION
      LOGICAL(4)                 :: TSO     ! CALCULATE WITH SPIN ORBIT COUPLING
      LOGICAL(4)                 :: TFOCK   ! CALCULATE WITH FOCK EXCHANGE
      INTEGER(4)          :: LMAP(19)=(/0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1/)
      REAL(8)                    :: FTOT
      REAL(8)                    :: PI,Y0,C0LL
      INTEGER(4)                 :: NBROYDENMEM
      REAL(8)                    :: BROYDENSTEP
      real(8)                    :: phifock(nr,nx)
      INTEGER(4)                 :: LRHOX=4
      character(128)             :: string
      real(8)                    :: scale
      logical                    :: tsecond
!     **************************************************************************
!
!     ==========================================================================
!     == RESOLVE KEY                                                          ==
!     ==========================================================================
      TREL=INDEX(KEY,'NONREL').EQ.0
      IF(TREL.AND.INDEX(KEY,'REL').EQ.0) THEN
        CALL ERROR$STOP('ATOMLIB$AESCF')
      END IF
      TSO=INDEX(KEY,'NONSO').EQ.0
      IF(TSO.AND.INDEX(KEY,'SO').EQ.0) THEN
        CALL ERROR$STOP('ATOMLIB$AESCF')
      END IF
      TFOCK=INDEX(KEY,'FOCK').ne.0
      if(tfock) then
        isvar=index(KEY,'FOCK')+5
        string=key(isvar:)
        isvar=index(KEY,'FOCK')
        read(string,*)scale
      else
        scale=0.d0
      end if
print*,'key=',trim(key)
print*,'scale ',scale,tfock
      TSTART=INDEX(KEY,'START').NE.0
      NBROYDENMEM=2
      BROYDENSTEP=5.D-1
!
!     ==========================================================================
!     == INITIALIZATIONS                                                      ==
!     ==========================================================================
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=Y0
      CALL RADIAL$R(GID,NR,R)
      IRBOX=0
      DO IR=1,NR-1
        IRBOX=IR
        IF(R(IR).GT.RBOX) EXIT
      ENDDO
!
      EOFI=0.D0
      IF(TSTART) THEN
        LOFI(:)=0
        F(:)=0.D0
        SO(:)=0
        NB=0
        FTOT=AEZ
        DO I=1,19
          NB=NB+1
          IF(NB.GT.NX) THEN
            CALL ERROR$MSG('ACTUAL NUMBER OF BAND EXCEEDS DIMENSION')
            CALL ERROR$I4VAL('NX',NX)
            CALL ERROR$STOP('ATOMLIB$AESCF')
          END IF
          LOFI(NB)=LMAP(I)
          IF(TSO)THEN
            IF(LMAP(I).EQ.0)THEN
              SO(NB)=1
              F(NB)=MIN(FTOT,2.D0)
              FTOT=FTOT-F(NB)
            ELSE
              SO(NB)=-1
              F(NB)=MIN(FTOT,REAL(2*LMAP(I),KIND=8))
              FTOT=FTOT-F(NB)
              NB=NB+1
              IF(NB.GT.NX) THEN
                CALL ERROR$MSG('ACTUAL NUMBER OF BAND EXCEEDS DIMENSION')
                CALL ERROR$I4VAL('NX',NX)
                CALL ERROR$STOP('ATOMLIB$AESCF')
              END IF
              LOFI(NB)=LMAP(I)
              SO(NB)=1
              F(NB)=MIN(FTOT,REAL(2*LMAP(I)+2,KIND=8))
              FTOT=FTOT-F(NB)
            END IF
          ELSE
            F(NB)=MIN(FTOT,REAL(2*(2*LMAP(I)+1),KIND=8))
            FTOT=FTOT-F(NB)
            SO(NB)=0
          END IF
          IF(FTOT.LE.1.D-10) EXIT
        ENDDO
        IF(FTOT.GT.1.D-10) THEN
          CALL ERROR$MSG('SPECIFIED NUMBER OF BANDS IS TOO SMALL')
          CALL ERROR$MSG('#(ELECTRONS) REMAINING')
          CALL ERROR$I4VAL('NX',NX)
          CALL ERROR$R8VAL('AEZ',AEZ)
          CALL ERROR$STOP('ATOMLIB$AESCF')
        END IF
        CALL RADIAL$NUCPOT(GID,NR,AEZ,POT)
        DO L=0,MAXVAL(LOFI)
          DO ISO=-1,1
            IF(TSO.AND.ISO.EQ.0) CYCLE
            IF(.NOT.TSO.AND.ISO.NE.0) CYCLE
            ISVAR=0
            DO IB=1,NB
              IF(LOFI(IB).NE.L.OR.SO(IB).NE.ISO) CYCLE
              NN(IB)=ISVAR
              ISVAR=ISVAR+1
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     == SELF-CONSISTENCY LOOP                                                ==
!     ==========================================================================
      tsecond=.false.
1000  continue
      XAV=0.D0
      XMAX=0.D0
      CONVG=.FALSE.
      CALL BROYDEN$NEW(NR,NBROYDENMEM,BROYDENSTEP)
      DO ITER=1,NITER
!
!       ========================================================================
!       == DETERMINE BOUND STATES FOR A GIVEN POTENTIAL AND ADD TO DENSITY    ==
!       ========================================================================
        DO IB=1,NB
          IF(TREL) THEN
            CALL SCHROEDINGER$DREL(GID,NR,POT,EOFI(IB),DREL)
          ELSE 
            DREL(:)=0.D0 
          END IF
          G(:)=0.D0
          IF(TFOCK.AND.tsecond) THEN
!print*,'fockboundstate ',iter,ib
            CALL ATOMLIB$FOCKBOUNDSTATE(GID,NR,LOFI(IB),SO(ib),RBOX,DREL &
     &                                 ,scale,NB,LOFI,F,PHIfock,MUX,LRHOX &
     &                                 ,G,NN(IB),POT,EOFI(IB),phi(:,ib))
          ELSE
!print*,'nofockboundstate ',iter,ib
            CALL ATOMLIB$BOUNDSTATE(GID,NR,LOFI(IB),SO(IB),RBOX,DREL,G,NN(IB) &
     &                           ,POT,EOFI(IB),PHI(:,IB))
          END IF
          IF(TREL) THEN
            G(:)=0.D0
            CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,LOFI(IB),SO(IB) &
     &                                         ,DREL,G,PHI(:,IB),SPHI(:,IB))
          ELSE
            SPHI(:,IB)=0.D0
          END IF
          AUX(:)=R(:)**2*(PHI(:,IB)**2+SPHI(:,IB)**2)
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
          PHI(:,IB)=PHI(:,IB)/SQRT(SVAR)
          SPHI(:,IB)=SPHI(:,IB)/SQRT(SVAR)
        ENDDO
        if(tfock)phifock=phi
!
!       ========================================================================
!       == ADD UP CHARGE DENSITY                                              ==
!       ========================================================================
        RHO(:)=0.D0
        DO IB=1,NB
          RHO(:)=RHO(:)+F(IB)*C0LL*(PHI(:,IB)**2+SPHI(:,IB)**2)
        ENDDO
!
!       ========================================================================
!       ==  CALCULATE ENERGY                                                  ==
!       ========================================================================
        IF(TBROYDEN) POTIN=POT
        IF(CONVG) THEN
          AUX(:)=0.D0
          DO IB=1,NB
            AUX(:)=AUX(:)+(PHI(:,IB)**2+SPHI(:,IB)**2)*(EOFI(IB)-POT(:)*Y0)*F(IB)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,EKIN)
          CALL ATOMLIB$BOXVOFRHO(GID,NR,RBOX,AEZ,RHO,POT,EH,EXC)
          ETOT=EKIN+EH+EXC
!!$PRINT*,'ETOT ',ITER,ETOT
        END IF
!
!       ========================================================================
!       ==  EXIT IF CONVERGED                                                 ==
!       ========================================================================
        IF(CONVG) EXIT

!       ========================================================================
!       == CALCULATE OUTPUT POTENTIAL                                         ==
!       ========================================================================
        CALL ATOMLIB$BOXVOFRHO(GID,NR,RBOX,AEZ,RHO,POT,EH,EXC)
        IF(TFOCK) CALL ATOMLIB$BOXMUX(GID,NR,RBOX,RHO,MUX)
!
!       ========================================================================
!       ==  GENERATE NEXT ITERATION USING D. G. ANDERSON'S METHOD             ==
!       ========================================================================
!print*,'xav ',xav,xmax
        XAV=SQRT(SUM(R**3*(POT-POTIN)**2)/SUM(R**3))
        XMAX=MAXVAL(ABS(R**2*(POT-POTIN))) 
        CALL BROYDEN$STEP(NR,POTIN,(POT-POTIN))
        POT=POTIN
        CONVG=(XMAX.LT.TOL)
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
      if(tfock.and.(.not.tsecond)) then
        if(tfock)phifock=phi
        tsecond=.true.
!        goto 1000
      end if
!!$DO I=1,NB
!!$ WRITE(*,FMT='(3I4,F10.2,I5,F20.3)')I,LOFI(I),SO(I),F(I),NN(I),EOFI(I)
!!$ENDDO
!!$PRINT*,'#ITERATIONS ',ITER
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
        ISVAR=SO(JB)
        SO(IB+1:JB)=SO(IB:JB-1)
        SO(IB)=ISVAR
        ISVAR=NN(JB)
        NN(IB+1:JB)=NN(IB:JB-1)
        NN(IB)=ISVAR
        SVAR=EOFI(JB)
        EOFI(IB+1:JB)=EOFI(IB:JB-1)
        EOFI(IB)=SVAR
        SVAR=F(JB)
        F(IB+1:JB)=F(IB:JB-1)
        F(IB)=SVAR
        AUX=PHI(:,JB)
        PHI(:,IB+1:JB)=PHI(:,IB:JB-1)
        PHI(:,IB)=AUX
        AUX=SPHI(:,JB)
        SPHI(:,IB+1:JB)=SPHI(:,IB:JB-1)
        SPHI(:,IB)=AUX
      ENDDO
return
IB=1
G(:)=0.D0
PRINT*,'AEZ ',AEZ,RBOX
PRINT*,'ib ',ib,lofi(ib),eofi(ib),nn(ib),f(ib)
print*,'tfock ',tfock
!CALL ATOMLIB$FOCKSCHROEDINGER(GID,NR,POT,DREL &
!     &        ,1.D0,NB,LOFI,F,PHI,MUX,LRHOX,LOFI(IB),SO(IB),EOFI(IB),AUX)
!CALL ATOMLIB$FOCKBOUNDSTATE(GID,NR,LOFI(IB),SO,RBOX,DREL &
!     &                                 ,0.5d0,NB,LOFI,F,PHI,MUX,LRHOX &
!     &                                 ,G,NN(IB),POT,EOFI(IB),AUX)
!CALL SETUP_WRITEPHI('FOCKPHI.DAT',GID,NR,1,AUX)
!CALL SETUP_WRITEPHI('NOFOCKPHI.DAT',GID,NR,nb,PHI(:,:))
!CALL SETUP_WRITEPHI('wFOCKPHI.DAT',GID,NR,nb,PHI(:,:))
!STOP 'done'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$FOCKSCHROEDINGER(GID,NR,POT,DREL &
     &                          ,SCALE,NB,LOFI,FOFI,AEPSI,MUX,LRHOX,L,SO,E,PHI)
!     **************************************************************************
!     **  SOLVES THE RADIAL SCHRODINGER EQUATION INCLUDING                    **
!     **  A CONTRIBUTION FROM THE FOCK-OPERATOR                               **
!     **     H = P^2/(2M) + POT + SCALE*(V_FOCK-MU_X)
!     **************************************************************************
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: POT(NR)   ! POTENTIAL
      REAL(8)   ,INTENT(IN) :: DREL(NR)  ! RELATIVISTIC CORRECTION
      REAL(8)   ,INTENT(IN) :: SCALE     ! PREFACTOR FOR FOCK-CORRECTION
      INTEGER(4),INTENT(IN) :: NB        ! #(STATES IN THE FOCK OPERATOR)
      INTEGER(4),INTENT(IN) :: LOFI(NB)  ! ANGULAR MOMENTA OF STATES
      REAL(8)   ,INTENT(IN) :: FOFI(NB)  ! OCCUPATIONS OF STATES
      REAL(8)   ,INTENT(IN) :: AEPSI(NR,NB) ! STATES IN THE FOCK OPERATOR
      REAL(8)   ,INTENT(IN) :: MUX(NR)   ! EXCHANGE POTENTIAL
      INTEGER(4),INTENT(IN) :: LRHOX     ! X(ANGULAR MOMENTUM OF THE DENSITY)
      INTEGER(4),INTENT(IN) :: L         ! ANGULAR MOMENTUM OF SOLUTION
      INTEGER(4),INTENT(IN) :: SO        ! SPIN-ORBIT SWITCH (0,-1,1)
      REAL(8)   ,INTENT(IN) :: E         ! ENERGY
      REAL(8)   ,INTENT(OUT):: PHI(NR)   ! RESULTING SOLUTION
      REAL(8)               :: G(NR)
      REAL(8)               :: R(NR)     ! RADIAL GRID
      REAL(8)               :: VPSI(NR,NB)
      REAL(8)               :: ETA(NR,NB)
      REAL(8)               :: AMAT(NB,NB),BVEC(NB),XVEC(NB)
      REAL(8)               :: AUX(NR)
      INTEGER(4)            :: IB,IB1,IB2
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == VPSI=SCALE*(V_FOCK-MU_X)|PSI>
!     ==========================================================================
      CALL ATOMLIB$FOCK(GID,NR,NB,LOFI,FOFI,AEPSI,LRHOX,NB,LOFI,AEPSI,VPSI)
      DO IB=1,NB
        VPSI(:,IB)=SCALE*(VPSI(:,IB)-MUX(:)*AEPSI(:,IB))
      ENDDO
!
!     ==========================================================================
!     == (H0-E)|PHI>=0; (H0-E)|ETA>=|VPSI>
!     ==========================================================================
      G=0.D0
      CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,E,1,PHI)
      DO IB=1,NB
        G(:)=VPSI(:,IB)
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,E,1,ETA(:,IB))
      ENDDO
!
!     ==========================================================================
!     == DETERMINE COEFFICIENTS FOR |PHI> = |PHI> + SUM_I |ETA_I> C_I        ==
!     ==========================================================================
      DO IB1=1,NB
        AUX(:)=R(:)**2*VPSI(:,IB1)*PHI(:)
        CALL RADIAL$INTEGRAL(GID,NR,AUX,BVEC(IB1))
        DO IB2=1,NB
          AUX(:)=R(:)**2*VPSI(:,IB1)*(AEPSI(:,IB2)+ETA(:,IB2))
          CALL RADIAL$INTEGRAL(GID,NR,AUX,AMAT(IB1,IB2))
        ENDDO
      ENDDO
      CALL LIB$MATRIXSOLVER8(NB,NB,1,AMAT,XVEC,BVEC)
      DO IB=1,NB
        PHI(:)=PHI(:)+ETA(:,IB)*XVEC(IB)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$FOCK(GID,NR,NB,LOFI,FOFI,AEPSI,LRHOX,NF,LOF,F,G)
!     **************************************************************************
!     **  APPLIES THE FOCK-OPERATOR TO A SET OF FUNCTIONS                     **
!     **************************************************************************
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: LOFI(NB)
      REAL(8)   ,INTENT(IN) :: FOFI(NB) 
      REAL(8)   ,INTENT(IN) :: AEPSI(NR,NB)
      INTEGER(4),INTENT(IN) :: LRHOX
      INTEGER(4),INTENT(IN) :: NF
      INTEGER(4),INTENT(IN) :: LOF(NF)
      REAL(8)   ,INTENT(IN) :: F(NR,NF)
      REAL(8)   ,INTENT(OUT):: G(NR,NF)
      REAL(8)               :: PI
      REAL(8)               :: SVAR
      REAL(8)               :: CG   ! GAUNT COEFFICIENT
      INTEGER(4)            :: LXB,LXF
      INTEGER(4)            :: I,IB
      INTEGER(4)            :: LF,LB,LRHO
      INTEGER(4)            :: IMF,IMB,IMRHO
      INTEGER(4)            :: LMF,LMB,LMRHO
      REAL(8)               :: AUX1(NR),AUX2(NR)
      REAL(8)   ,ALLOCATABLE :: COEFF(:,:,:)
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      LXB=MAXVAL(LOFI)
      LXF=MAXVAL(LOF)
!
!     ==========================================================================
!     == TAKE CARE OF ANGULAR PART                                            ==
!     ==========================================================================
      ALLOCATE(COEFF(LXB+1,LRHOX+1,LXF+1))
      COEFF(:,:,:)=0.D0
      DO LF=0,LXF
        DO LB=0,LXB
          DO LRHO=0,LRHOX
!PRINT*,'========================================================='
            LMF=LF**2+1  ! RESULT INDEPEMDENT OF IMF: PIC ANY LMF
            SVAR=0.D0
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
            SVAR=4.D0*PI*SVAR/REAL((2*LB+1)*(2*LB+1),KIND=8) ! AVERAGE OVER M_J
!WRITE(6,FMT='(5I5,F20.5)')LF,LB,LRHO,IMB,IMRHO,SVAR
            COEFF(LB+1,LRHO+1,LF+1)=SVAR
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == TAKE CARE OF ANGULAR PART                                            ==
!     ==========================================================================
      DO I=1,NF
        LF=LOF(I)
        G(:,I)=0.D0
        DO IB=1,NB
          IF(FOFI(IB).LT.1.D-3) CYCLE
          LB=LOFI(IB)
          AUX1(:)=AEPSI(:,IB)*F(:,I)
          DO LRHO=0,LRHOX         
            IF(COEFF(LB+1,LRHO+1,LF+1).LT.1.D-8) CYCLE    !COEFF.GE.0
            CALL RADIAL$POISSON(GID,NR,LRHO,AUX1,AUX2)
            SVAR=REAL((2*LRHO+1)*(2*LB+1),KIND=8)*COEFF(LB+1,LRHO+1,LF+1)
            G(:,I)=G(:,I)+SVAR*AEPSI(:,IB)*AUX2(:)
          ENDDO
        ENDDO
      ENDDO
      G(:,:)=-0.5D0*G(:,:)/(4.D0*PI)
!
!CALL SETUP_WRITEPHI('FOCKF.DAT',GID,NR,NF,F)
!CALL SETUP_WRITEPHI('FOCKG.DAT',GID,NR,NF,G)
!STOP                             
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RBOX,DREL,G,NN,POT,E,PHI)
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
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      INTEGER(4) ,INTENT(IN)     :: NN      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)    ,PARAMETER      :: TOL=1.D-8
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: DER,DERO
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHIHOM(NR),PHIINHOM(NR),GHOM(NR)
      REAL(8)                    :: PHI1(NR)
      REAL(8)                    :: AUX(NR),AUX1(NR)
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRMATCH 
      REAL(8)                    :: SVAR
      REAL(8)                    :: POT1(NR)
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20 !MAX. FACTOR IN THE WAVEFUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
      LOGICAL(4)                 :: THOM
      REAL(8)                    :: ROUT
      REAL(8)                    :: VAL,VAL1,VAL2,R1,R2
      INTEGER(4)                 :: IROUT,IRCL,IRBOX,IREND
!     **************************************************************************
                                 CALL TRACE$PUSH('ATOMLIB$BOUNDSTATE')
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
!     ==  R(IRBOX) IS THE FIRST GRIDPOINT JUST IOUTSIDE THE BOX
      IRBOX=1
      DO IR=1,NR
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
        POT1(IROUT:)=E
!
!       ========================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       ========================================================================
        IDIR=1
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHI)
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
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
        PHIHOM(:)=PHIHOM(:)/PHIHOM(IRMATCH)
        GHOM(:)=0.D0
        GHOM(IREND+2)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHI1)
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
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHIINHOM)
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
      END IF
!!$!
! DO NOT NORMALIZE!!!! IT DOES NOT MAKE SENSE FOR THE INHOMOGENEOUS SOLUTION
!!$!     =======================================================================
!!$!     ==  NORMALIZE WAVE FUNCTION                                          ==
!!$!     =======================================================================
!!$      AUX(:)=R(:)**2*PHI(:)**2
!!$      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
!!$      PHI(:)=PHI(:)/SQRT(SVAR)
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$FOCKBOUNDSTATE(GID,NR,L,SO,RBOX,DREL &
     &                                 ,SCALE,NB,LOFI,FOFI,PSIOFI,MUX,LRHOX &
     &                                 ,G,NN,POT,E,PHI)
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
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: SCALE
      INTEGER(4) ,INTENT(IN)     :: NB
      INTEGER(4) ,INTENT(IN)     :: LOFI(NB)
      REAL(8)    ,INTENT(IN)     :: FOFI(NB)
      REAL(8)    ,INTENT(IN)     :: PSIOFI(NR,NB)
      REAL(8)    ,INTENT(IN)     :: MUX(NR)
      INTEGER(4) ,INTENT(IN)     :: LRHOX
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      INTEGER(4) ,INTENT(IN)     :: NN      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)    ,PARAMETER      :: TOL=1.D-8
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: DER,DERO
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHIHOM(NR),PHIINHOM(NR),GHOM(NR)
      REAL(8)                    :: PHI1(NR)
      REAL(8)                    :: AUX(NR),AUX1(NR)
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRMATCH 
      REAL(8)                    :: SVAR
      REAL(8)                    :: POT1(NR)
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20 !MAX. FACTOR IN THE WAVEFUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
      LOGICAL(4)                 :: THOM
      REAL(8)                    :: ROUT
      REAL(8)                    :: VAL,VAL1,VAL2,R1,R2
      INTEGER(4)                 :: IROUT,IRCL,IRBOX,IREND
real(8) :: test(nr,10)
!     **************************************************************************
                                 CALL TRACE$PUSH('ATOMLIB$BOUNDSTATE')
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
!     ==  R(IRBOX) IS THE FIRST GRIDPOINT JUST IOUTSIDE THE BOX
      IRBOX=1
      DO IR=1,NR
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
        POT1(IROUT:)=E
!
!       ========================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       ========================================================================
        IDIR=1
!        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHI)
         CALL ATOMLIB$FOCKSCHROEDINGER(GID,NR,POT1,DREL &
     &                               ,SCALE,NB,LOFI,FOFI,PSIOFI,MUX,LRHOX &
     &                               ,L,SO,E,PHI)
if(i.le.10)test(:,i)=phi
!
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
!!$print*,'e ',e,z0,dx
!!$z0=-real(nn)-0.5d0
!!$do ir=3,nr
!!$  if(phi(ir)*phi(ir-1).lt.0.d0) then
!!$    z0=z0+1.d0
!!$    print*,'node ',r(ir),z0
!!$  end if
!!$  if(r(ir).ge.rout) exit
!!$enddo
!!$ CALL radial$value(GID,NR,PHI,ROUT,svar)
print*,'e ',e,z0
!       ========================================================================
!       ==  BISECTION                                                         ==
!       ========================================================================
        IF(ABS(2.D0*DX).LE.TOL) EXIT
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
CALL SETUP_WRITEPHI('FOCKbound.DAT',GID,NR,10,test)
CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,ROUT,Z0)
CALL radial$value(GID,NR,PHI,ROUT,svar)
Z0=Z0-REAL(NN+1)
print*,'+-e ',e,z0,rout,svar
      if(z0.lt.0.d0) then
        e=x0+2.d0*dx
        CALL ATOMLIB$FOCKSCHROEDINGER(GID,NR,POT1,DREL &
     &                               ,SCALE,NB,LOFI,FOFI,PSIOFI,MUX,LRHOX &
     &                               ,L,SO,E,PHI)
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,ROUT,Z0)
        CALL radial$value(GID,NR,PHI,ROUT,svar)
        Z0=Z0-REAL(NN+1)
print*,'--e ',e,z0,rout,svar
      end if
return
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
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
        PHIHOM(:)=PHIHOM(:)/PHIHOM(IRMATCH)
        GHOM(:)=0.D0
        GHOM(IREND+2)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHI1)
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
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHIINHOM)
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
      END IF
!!$!
! DO NOT NORMALIZE!!!! IT DOES NOT MAKE SENSE FOR THE INHOMOGENEOUS SOLUTION
!!$!     =======================================================================
!!$!     ==  NORMALIZE WAVE FUNCTION                                          ==
!!$!     =======================================================================
!!$      AUX(:)=R(:)**2*PHI(:)**2
!!$      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
!!$      PHI(:)=PHI(:)/SQRT(SVAR)
                                 CALL TRACE$POP()
      RETURN
      END
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
      REAL(8)    ,PARAMETER      :: TOL=1.D-8
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHI1(NR),PHI2(NR)
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      REAL(8)                    :: ROUT
      REAL(8)                    :: VAL,VAL1,VAL2,R1,R2,SVAR
      INTEGER(4)                 :: IROUT,IRCL,IRBOX,IREND
!     **************************************************************************
                                 CALL TRACE$PUSH('ATOMLIB$PAWBOUNDSTATE')
      CALL RADIAL$R(GID,NR,R)
!     ==  R(IRBOX) IS THE FIRST GRIDPOINT JUST IOUTSIDE THE BOX
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
!     ==  CHOP OF TAILS WHICH MAY BE EXPONENTIALLY INCREASING                 ==
!     ==========================================================================
      CALL RADIAL$VALUE(GID,NR,PHI1,RBOX,VAL1)
      CALL RADIAL$VALUE(GID,NR,PHI2,RBOX,VAL2)
      SVAR=VAL2-VAL1
      VAL1=VAL1/SVAR
      VAL2=VAL2/SVAR
      PHI=PHI1*VAL2-PHI2*VAL1
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
      REAL(8)               :: ALPHA
      REAL(8),   PARAMETER  :: RHOMIN=1.D-2
      REAL(8)               :: POTH(NR)
      REAL(8)               :: POTXC(NR)
      REAL(8)               :: PI
      REAL(8)               :: FOURPI
      REAL(8)               :: Y0
      REAL(8)               :: RHO1(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: RHOPLUS(NR)
      REAL(8)               :: DRHOPLUSDRHO(NR)
      REAL(8)               :: EDEN(NR)
      REAL(8)               :: GRHO(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VGXC,VXC,EXC1,RH,GRHO2
      REAL(8)               :: DUMMY1,DUMMY2,DUMMY3
      REAL(8)               :: EVAL,MUVAL
      REAL(8)               :: SVAR
      REAL(8)               :: F,X,DF,YP,YM,NBYF,DFDN
      INTEGER(4)            :: IR,IRBOX
      REAL(8)               :: Q
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRBOX=IR
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
      RHO1(IR:)=0.D0
      CALL RADIAL$NUCPOT(GID,NR,AEZ,POTH)
      EDEN(:)=0.5D0*RHO1(:)*POTH(:)
      CALL RADIAL$POISSON(GID,NR,0,RHO1,AUX)
      POTH(:)=POTH(:)+AUX(:)
      CALL RADIAL$VALUE(GID,NR,POTH,RAD,SVAR)
      SVAR=Q/RAD/Y0-SVAR
      POTH(1:IRBOX-1)=POTH(1:IRBOX-1)+SVAR
      POTH(IRBOX:NR)=Q/R(IRBOX:NR)/Y0
      EDEN(:)=EDEN(:)+0.5D0*RHO1(:)*POTH(:)
!
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,EDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,EH)
!
      CALL RADIAL$DERIVE(GID,NR,RHO,GRHO)
      DO IR=1,NR
        RH=RHO(IR)*Y0
        GRHO2=(Y0*GRHO(IR))**2
        CALL DFT(RH,0.D0,GRHO2,0.D0,0.D0,EXC1,VXC,DUMMY1,VGXC,DUMMY2,DUMMY3)
        EDEN(IR)=4.D0*PI*EXC1   ! ANGULAR INTEGRATION ALREADY INCLUDED
        POTXC(IR)=VXC/Y0
        GRHO(IR)=VGXC*2.D0*GRHO(IR)
      ENDDO
      CALL RADIAL$DERIVE(GID,NR,GRHO(:),AUX)
      POTXC(:)=POTXC(:)-AUX
      IF(R(1).GT.1.D-8) THEN
         POTXC(:)=POTXC(:)-2.D0/R(:)*GRHO(:)
      ELSE
        POTXC(2:)=POTXC(2:)-2.D0/R(2:)*GRHO(2:)
!        POTXC(1)=POTXC(1)-2.D0/R(2)*GRHO(2)
        POTXC(1)=POTXC(2)
      END IF
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
      REAL(8)               :: ALPHA
      REAL(8),   PARAMETER  :: RHOMIN=1.D-2
      REAL(8)               :: PI
      REAL(8)               :: FOURPI
      REAL(8)               :: Y0
      REAL(8)               :: AUX(NR)
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: GRHO(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VGXC,VXC,EXC1,RH,GRHO2
      REAL(8)               :: DUMMY1,DUMMY2,DUMMY3
      REAL(8)               :: EVAL,MUVAL
      REAL(8)               :: SVAR
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
      CALL RADIAL$DERIVE(GID,NR,GRHO(:),AUX)
      MUX(:)=MUX(:)-AUX
      IF(R(1).GT.1.D-8) THEN
         MUX(:)=MUX(:)-2.D0/R(:)*GRHO(:)
      ELSE
        MUX(2:)=MUX(2:)-2.D0/R(2:)*GRHO(2:)
        MUX(1)=MUX(2)
      END IF
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
      INTEGER(4)              :: LXX
      INTEGER(4)              :: L
      REAL(8)                 :: E
      REAL(8)                 :: G(NR)
      REAL(8)                 :: DREL(NR)
      INTEGER(4)              :: NN
      INTEGER(4)              :: SO
      INTEGER(4)              :: IB,I
      REAL(8)                 :: PI,Y0
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      LXX=MAXVAL(LOFI(:))
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
        NN=0
        CALL ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RBOX,DREL,G,NN,POT,E,UOFI(:,IB))
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
      INTEGER(4)              :: IB,JB
      INTEGER(4)              :: I1,I2,I3
      REAL(8)                 :: SVAR
      REAL(8)                 :: AUX(NR)
      REAL(8)                 :: R(NR)
      REAL(8)                 :: DREL(NR)
      INTEGER(4)              :: SO
INTEGER(4)              :: IR
REAL(8)                 :: PI,Y0
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
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
