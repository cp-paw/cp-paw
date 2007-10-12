!     MODULE SCHGRL_INTERFACE_MODULE
!     PAWBOUNDSTATE
!        -> SCHROEDER
!        -> PAWDER
!     PAWDER
!       ->SCHROEDER
!     SCHROEDINHOM(NOT USED)f
!       ->SCHROEDER
!     FIXNODE
!     -> AEBOUNDSTATE
!        -> SCHROEDER
!           -> VTRANSFORM
!           -> GTRANSFORM
!           -> OUTWRD
!           -> INWRD
!           -> OUTBOUNDARY
!           -> INBOUNDARY
!!!!!!!!!
      MODULE SCHRGL_INTERFACE_MODULE
      INTERFACE
        SUBROUTINE SCHROEDER(TFRWRD,THOM,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI,GINH,DLG)
          LOGICAL,   INTENT(IN)   :: THOM  ! IGNORE INHOMOGENITY?   
          LOGICAL,   INTENT(IN)   :: TFRWRD! INTEGRATE OUTWARD/INWARD
          REAL(8),    INTENT(IN)   :: R1    ! RADIAL GRID R(I)=R1*EXP(DEX*(I-1))
          REAL(8),    INTENT(IN)   :: DEX   !
          INTEGER(4), INTENT(IN)   :: NR    ! NUMBER OF GRID POINTS
          INTEGER(4), INTENT(IN)   :: L     ! MAIN ANGULAR MOMENTUM
          REAL(8),    INTENT(IN)   :: AEZ   ! ATOMIC NUMBER (AEZ=0->NON-RELATIVISTIC)
          INTEGER(4), INTENT(IN)   :: MCH   ! INTEGRATION STOPS AT IR=MCH
          REAL(8),    INTENT(IN)   :: E     ! ENERGY
          REAL(8),    INTENT(IN)   :: POT(NR)    ! POTENTIAL
          REAL(8),    INTENT(IN)   :: GINH(NR)   !INHOMOGENEITY
          REAL(8),    INTENT(IN)   :: BOXRADIUS  !
          REAL(8),    INTENT(IN),OPTIONAL :: DLG ! LOGARITHMIC DERIVATIVE
          REAL(8),    INTENT(OUT)  :: PHI(NR,3)  ! PHI,DPHI/DR,NABLA**2 PHI
        END SUBROUTINE SCHROEDER     
      END INTERFACE
      INTERFACE
        SUBROUTINE AEBOUNDSTATE(R1,DEX,NR,L,E,AEZ,BOXRADIUS,POT,PHI,DLG)
          REAL(8),    INTENT(IN)   :: R1         ! RADIAL GRID R(I)=R1*EXP(DEX*(I-1))
          REAL(8),    INTENT(IN)   :: DEX        !
          INTEGER(4), INTENT(IN)   :: NR         ! NUMBER OF GRID POINTS
          INTEGER(4), INTENT(IN)   :: L          ! MAIN ANGULAR MOMENTUM
          REAL(8),    INTENT(INOUT):: E          ! ENERGY
          REAL(8),    INTENT(IN)   :: AEZ        ! ATOMIC NUMBER (AEZ=0->NON-RELATIVISTIC)
          REAL(8),    INTENT(IN)   :: POT(NR)    ! POTENTIAL
          REAL(8),    INTENT(OUT)  :: PHI(NR,3)  ! BOUND STATE, DERIVATIVE, LAPLACE
          REAL(8),    INTENT(IN)   :: BOXRADIUS  ! BOUNDARY COND.: PHI(BOXRADIUS)=0
          REAL(8),    INTENT(IN),OPTIONAL :: DLG ! LOGARITHMIC DERIVATIVE
        END SUBROUTINE AEBOUNDSTATE
      END INTERFACE
      INTERFACE
        SUBROUTINE PAWBOUNDSTATE(R1,DEX,NR,L,E,BOXRADIUS &
     &                      ,POT,NPRO,PRO,DATH,DO,PHI,DLG)
          REAL(8)   ,INTENT(IN)   :: R1    ! RADIAL GRID R(I)=R1*EXP(DEX*(I-1))
          REAL(8)   ,INTENT(IN)   :: DEX   !
          INTEGER(4),INTENT(IN)   :: NR    ! NUMBER OF GRID POINTS
          INTEGER(4),INTENT(IN)   :: L     ! MAIN ANGULAR MOMENTUM
          REAL(8)   ,INTENT(INOUT):: E     ! ENERGY
          REAL(8)   ,INTENT(IN)   :: POT(NR)      ! POTENTIAL
          INTEGER(4),INTENT(IN)   :: NPRO         ! MAIN ANGULAR MOMENTUM
          REAL(8)   ,INTENT(IN)   :: PRO(NR,NPRO) ! PROJECTOR FUNCTIONS
          REAL(8)   ,INTENT(IN)   :: DATH(NPRO,NPRO)
          REAL(8)   ,INTENT(IN)   :: DO(NPRO,NPRO)
          REAL(8)   ,INTENT(OUT)  :: PHI(NR,3)    ! BOUND STATE, DERIVATIVE, LAPLACE
          REAL(8)   ,INTENT(IN)   :: BOXRADIUS    ! BOUNDARY COND.: PHI(BOXRADIUS)=0
          REAL(8)   ,INTENT(IN),OPTIONAL  :: DLG  ! LOGARITHMIC DERIVATIVE
        END SUBROUTINE PAWBOUNDSTATE
      END INTERFACE
      INTERFACE
        SUBROUTINE FIXNODE(R1,DEX,NR,L,E,AEZ,BOXRADIUS,POT,PHI,ZEFF,NN,TOL,DLG)
          REAL(8),   INTENT(IN)   :: R1
          REAL(8),   INTENT(IN)   :: DEX
          INTEGER(4),INTENT(IN)   :: NR
          INTEGER(4),INTENT(IN)   :: L
          INTEGER(4),INTENT(IN)   :: NN
          REAL(8),   INTENT(IN)   :: AEZ
          REAL(8),   INTENT(INOUT):: E
          REAL(8),   INTENT(IN)   :: BOXRADIUS
          REAL(8),   INTENT(IN), OPTIONAL :: DLG
          REAL(8),   INTENT(IN)   :: TOL
          REAL(8),   INTENT(IN)   :: POT(NR)
          REAL(8),   INTENT(IN)   :: ZEFF
          REAL(8),   INTENT(OUT)  :: PHI(NR,3)
        END SUBROUTINE FIXNODE
      END INTERFACE
      INTERFACE
        SUBROUTINE OUTBOUNDARY(THOM,R1,DEX,NR,R,IRBND,BOXRADIUS,U,UP,UPP,CF,DLG)
          LOGICAL,  INTENT(IN)   :: THOM      ! SWITCH FOR HOMOGENEOUS SOLUTION
          INTEGER(4),INTENT(IN)   :: NR        ! #(RADIAL GRID POINTS)
          INTEGER(4),INTENT(IN)   :: IRBND     ! SOLVE FROM IRBND TO NR
          REAL(8),   INTENT(IN)   :: R1        ! FIRST POINT ON RAIDAL GRID
          REAL(8),   INTENT(IN)   :: DEX       ! 
          REAL(8),   INTENT(IN)   :: BOXRADIUS ! RADIUS FOR WHICH BOUNDARY CONDITIONS ARE SPECIFIED
          REAL(8),   INTENT(IN)   :: R(NR)     ! RADIAL GRID
          REAL(8),   INTENT(IN),OPTIONAL :: DLG  ! LOGARITHMIC DERIVATIVE AT BOXRADIUS
          REAL(8),   INTENT(INOUT):: U(NR)
          REAL(8),   INTENT(INOUT):: UP(NR)
          REAL(8),   INTENT(INOUT):: UPP(NR)
          REAL(8),   INTENT(IN)   :: CF(NR)
        END SUBROUTINE OUTBOUNDARY
      END INTERFACE
      END MODULE SCHRGL_INTERFACE_MODULE
!
!     ..................................................................
      SUBROUTINE PAWBOUNDSTATE(R1,DEX,NR,L,E,BOXRADIUS &
     &                      ,POT,NPRO,PRO,DATH,DO,PHI,DLG)
!     **                                                              **
!     **                                                              **
!     **                                                              **
      USE SCHRGL_INTERFACE_MODULE, ONLY : SCHROEDER
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: R1    ! RADIAL GRID R(I)=R1*EXP(DEX*(I-1))
      REAL(8)   ,INTENT(IN)   :: DEX   !
      INTEGER(4),INTENT(IN)   :: NR    ! NUMBER OF GRID POINTS
      INTEGER(4),INTENT(IN)   :: L     ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(INOUT):: E     ! ENERGY
      REAL(8)   ,INTENT(IN)   :: POT(NR)      ! POTENTIAL
      INTEGER(4),INTENT(IN)   :: NPRO         ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN)   :: PRO(NR,NPRO) ! PROJECTOR FUNCTIONS
      REAL(8)   ,INTENT(IN)   :: DATH(NPRO,NPRO)
      REAL(8)   ,INTENT(IN)   :: DO(NPRO,NPRO)
      REAL(8)   ,INTENT(OUT)  :: PHI(NR,3)    ! BOUND STATE, DERIVATIVE, LAPLACE
      REAL(8)   ,INTENT(IN)   :: BOXRADIUS    ! BOUNDARY COND.: PHI(BOXRADIUS)=0
      REAL(8)   ,INTENT(IN),OPTIONAL  :: DLG  ! LOGARITHMIC DERIVATIVE
      REAL(8)   ,PARAMETER    :: AEZ=0.D0     
      REAL(8)                 :: PHI1(NR,3)   ! PARTIAL WAVE OUTWARD INTEGR.
      REAL(8)                 :: PHI2(NR,3)   ! PARTIAL WAVE INWARD INTEGR
      REAL(8)                 :: PHI1DOT(NR,3)! ENERGY DERIVATIVE OF PHI1
      REAL(8)                 :: PHI2DOT(NR,3)! ENERGY DERIVATIVE OF PHI2
      REAL(8)                 :: GINH(NR)     ! INHOMOGENEITY, ALSO USED AS DUMMY
      INTEGER(4)              :: MCH          ! MATCHING GRID POINT
      REAL(8)                 :: SVAR,SVAR1,SVAR2,SVAR3,DE,STEP(3)
      INTEGER(4)              :: IR,I,I1,I2,IPRO
      REAL(8)                 :: XEXP,RI,Y0,PI
      REAL(8)                 :: AUX(NR)
      REAL(8)                 :: r(NR)
      REAL(8)                 :: VEC(NPRO)
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
      XEXP=DEXP(DEX)
      ri=r1/xexp
      do ir=1,nr
        ri=ri*xexp
        r(ir)=ri
      enddo
!
!     ==================================================================
!     == FIND CLASSICAL TURNING POINT                                 ==
!     ==================================================================
      MCH=4
      DO IPRO=1,NPRO
        DO IR=4,NR
          IF(DABS(PRO(IR,IPRO)).GT.1.D-12)MCH=MAX(MCH,IR)
        ENDDO
      ENDDO
      MCH=MIN(MCH,NR-4)
      IF(R1*DEXP(DEX*DBLE(MCH-1)).GT.BOXRADIUS) THEN
         call error$STOP('PAWBOUNDSTATE')
      END IF
!
!     ==================================================================
!     == INTEGRATE SCHROEDINGER EQUATION INWARD AND OUTWARD           ==
!     ==================================================================
!     __ INTEGRATE OUTWARD
      CALL PAWDER(.TRUE.,R1,DEX,NR,L,E,MCH &
     &                     ,POT,NPRO,PRO,DATH,DO,PHI1,GINH)
!     __ INTEGRATE INWARD
      IF(PRESENT(DLG)) THEN
        CALL SCHROEDER(.FALSE.,.TRUE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI2,GINH,DLG=DLG)
      ELSE
        CALL SCHROEDER(.FALSE.,.TRUE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI2,GINH)
      END IF
      SVAR=PHI1(MCH,1)/PHI2(MCH,1)
      PHI2(:,:)=PHI2(:,:)*SVAR
!
!     ==================================================================
!     ==  CALCULATE INHOMOGINEITY FOR PHIDOT CALCULATION              ==
!     ==  GINH=[1+|PRO>DO<PRO|] |PHI>                                 ==
!     ==  (H-E)|PHIDOT>=|GINH>                                        ==
!     ==================================================================
      PHI(1:MCH,1) =PHI1(1:MCH,1)
      PHI(MCH:NR,1)=PHI2(MCH:NR,1)
      DO I1=1,NPRO
        CALL OLDRADIAL$INTEGRAL(R1,DEX,NR,R(:)**2*PRO(:,I1)*PHI(:,1),VEC(I1))
      ENDDO
      GINH(:)=PHI(:,1)
      DO I1=1,NPRO
        SVAR=0.D0
        DO I2=1,NPRO
          SVAR=SVAR+DO(I1,I2)*VEC(I2)
        ENDDO
        GINH(:)=GINH(:)+PRO(:,I1)*SVAR
      ENDDO
!
!     == INTEGRATE PHIDOT OUTWARD ======================================
      CALL PAWDER(.FALSE.,R1,DEX,NR,L,E,MCH &
     &                   ,POT,NPRO,PRO,DATH,DO,PHI1DOT,GINH)
!
!     == INTEGRATE PHIDOT INWARD =======================================
      GINH=PHI2(:,1)
      IF(PRESENT(DLG)) THEN
        CALL SCHROEDER(.FALSE.,.FALSE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI2DOT,GINH,DLG=DLG)
      ELSE
        CALL SCHROEDER(.FALSE.,.FALSE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI2DOT,GINH)
      END IF
!
!     == MATCH OUTER AND INNER VALUE OF PHIDOT =========================
      SVAR=(PHI1DOT(MCH,1)-PHI2DOT(MCH,1))/PHI2(MCH,1)
      PHI2DOT(:,:)=PHI2DOT(:,:)+PHI2(:,:)*SVAR
!
!     ==================================================================
!     == MATCH DERIVATIVE                                             ==
!     ==================================================================
      DE=(PHI1(MCH,2)-PHI2(MCH,2))/(PHI2DOT(MCH,2)-PHI1DOT(MCH,2))
      PHI(1:MCH,:)=PHI1(1:MCH,:)+PHI1DOT(1:MCH,:)*DE
      PHI(MCH:NR,:)=PHI2(MCH:NR,:)+PHI2DOT(MCH:NR,:)*DE
      STEP(:)=PHI1(MCH,:)+PHI1DOT(MCH,:)*DE &
     &       -PHI2(MCH,:)-PHI2DOT(MCH,:)*DE
!     print*,'step',step,de
!
!     ==================================================================
!     == NORMALIZE                                                    ==
!     ==================================================================
      DO I1=1,NPRO
        CALL OLDRADIAL$INTEGRAL(R1,DEX,NR,r(:)**2*pro(:,i1)*phi(:,1),VEC(I1))
      ENDDO
      CALL OLDRADIAL$INTEGRAL(R1,DEX,NR,(phi(:,1)*r(:))**2,SVAR)
      DO I1=1,NPRO
        DO I2=1,NPRO
          SVAR=SVAR+VEC(I1)*DO(I1,I2)*VEC(I2)
        ENDDO
      ENDDO
      SVAR=1.D0/DSQRT(SVAR)
      PHI(:,:)=PHI(:,:)*SVAR
      STEP(:)=STEP(:)*SVAR
      IF(DABS(STEP(3)).GT.1.D-3) THEN
        PRINT*,'WARNING FRom PAWBOUNDSTATE:',STEP(3),DE
      END IF
!
!     ==================================================================
!     == EVALUATE NEW ENERGY                                          ==
!     ==================================================================
      E=E+DE
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PAWDER(THOM,R1,DEX,NR,L,E,MCH &
     &                     ,POT,NPRO,PRO,DATH,DO,PHI,GINH)
!     **                                                              **
!     **  SOLVES THE INHOMOGENEOUS SCHROEDINGER EQUATION              **
!     **                                                              **
      USE SCHRGL_INTERFACE_MODULE
      IMPLICIT NONE
      LOGICAL(4)  ,INTENT(IN) :: THOM
      REAL(8)     ,INTENT(IN) :: R1      
      REAL(8)     ,INTENT(IN) :: DEX      
      INTEGER(4)  ,INTENT(IN) :: NR
      INTEGER(4)  ,INTENT(IN) :: L
      REAL(8)     ,INTENT(IN) :: E
      INTEGER(4)  ,INTENT(IN) :: MCH
      REAL(8)     ,INTENT(IN) :: POT(NR)
      INTEGER(4)  ,INTENT(IN) :: NPRO
      REAL(8)     ,INTENT(IN) :: PRO(NR,NPRO)
      REAL(8)     ,INTENT(IN) :: DATH(NPRO,NPRO)     
      REAL(8)     ,INTENT(IN) :: DO(NPRO,NPRO)     
      REAL(8)     ,INTENT(IN) :: GINH(NR)
      REAL(8)     ,INTENT(OUT):: PHI(NR,3)
      REAL(8)                 :: U(NR,3)
      REAL(8)                 :: V(NR,3,NPRO)
      REAL(8)                 :: AMAT(NPRO,NPRO)
      REAL(8)                 :: BMAT(NPRO,NPRO)
      REAL(8)                 :: BMATinv(NPRO,NPRO)
      REAL(8)                 :: CMAT(NPRO,NPRO)
      REAL(8)                 :: CVEC(NPRO)
      REAL(8)                 :: DVEC(NPRO)
      INTEGER(4)              :: I1,I2,I3,IR
      REAL(8)                 :: SVAR
      INTEGER(4),PARAMETER    :: NAUX=10000
      REAL(8)                 :: DET(2),RCOND,AUX1(NAUX)
      REAL(8)                 :: AUX(NR)
      REAL(8)                 :: RI,XEXP
      XEXP=DEXP(DEX)
!
!     ==================================================================
!     ==  TEST                                                        ==
!     ==================================================================
      DO I1=1,NPRO
        DO IR=MCH+1,NR
exit
          IF(DABS(PRO(IR,I1)).GT.1.D-6) THEN
            CALL ERROR$MSG('WARNING FROM PAWDER PRO DOES NOT VANISH OUTSIDE')
            CALL ERROR$I4VAL('L',L)
            CALL ERROR$I4VAL('IPRO',I1)
            CALL ERROR$I4VAL('IR',IR)
            CALL ERROR$I4VAL('MCH',MCH)
            CALL ERROR$r8VAL('r(mcH)',r1*exp(dex*real(mch-1)))
            CALL ERROR$STOP('PAWDER')
          END IF
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  -1/2NABLA^2+POT-E|U>=|G>                                    ==
!     ==================================================================
      CALL SCHROEDER(.TRUE.,THOM,R1,DEX,NR,L,E,0.D0,MCH,0.D0,POT,U,GINH)
!
!     ==================================================================
!     ==  -1/2NABLA^2+POT-E|V>=|PRO>                                 ==
!     ==================================================================
      DO I1=1,NPRO
        CALL SCHROEDER(.TRUE.,.FALSE.,R1,DEX,NR,L,E,0.D0,MCH,0.D0 &
     &                     ,POT,V(:,:,I1),PRO(:,I1))
      ENDDO
!
!     ==================================================================
!     ==  AMAT=<PRO|V>  CVEC=<PRO|U>                                  ==
!     ==================================================================
      DO I1=1,NPRO
        DO I2=1,NPRO
          RI=R1/XEXP
          DO IR=1,NR
            RI=RI*XEXP
            AUX(IR)=RI**2*PRO(IR,I1)*V(IR,1,I2)
          ENDDO
          CALL OLDRADIAL$INTEGRAL(R1,DEX,NR,AUX,AMAT(I1,I2))
        ENDDO
      ENDDO
      DO I1=1,NPRO
        RI=R1/XEXP
        DO IR=1,NR
          RI=RI*XEXP
          AUX(IR)=RI**2*PRO(IR,I1)*U(IR,1)
        ENDDO
        CALL OLDRADIAL$INTEGRAL(R1,DEX,NR,AUX,CVEC(I1))
      ENDDO  
!
!     ==================================================================
!     ==  BMAT=1+(DATH-EDO)<PRO|V>                                    ==
!     ==================================================================
      DO I1=1,NPRO
        DO I2=1,NPRO
          BMAT(I1,I2)=0.D0
          DO I3=1,NPRO
            BMAT(I1,I2)=BMAT(I1,I2)+(DATH(I1,I3)-E*DO(I1,I3))*AMAT(I3,I2)     
          ENDDO
        ENDDO
        BMAT(I1,I1)=BMAT(I1,I1)+1.D0
      ENDDO
!
!     ==================================================================
!     ==  BMAT = BMAT^-1 = [1+(DATH-EDO)<PRO|V>]^-1                   ==
!     ==================================================================
      IF(NPRO.EQ.1) THEN
        BMAT(1,1)=1.D0/BMAT(1,1)
      else IF(NPRO.EQ.0) THEN
      ELSE 
!        CALL DGEICD(BMAT,NPRO,NPRO,0,RCOND,DET,AUX1,NAUX,*9999)
        call lib$invertr8(npro,bmat,bmatinv)
        bmat=bmatinv
      END IF
!
!     ==================================================================
!     ==  CMAT = -BMAT*(DATH-EDO)                                     ==
!     ==       = -[1+(DATH-EDO)*<PRO|V>]^-1 (DATH-EDO)                ==
!     ==================================================================
      DO I1=1,NPRO
        DO I2=1,NPRO
          CMAT(I1,I2)=0.D0
          DO I3=1,NPRO
            CMAT(I1,I2)=CMAT(I1,I2)-BMAT(I1,I3)*(DATH(I3,I2)-E*DO(I3,I2))
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  DVEC = CMAT*CVEC                                           ==
!     ==       = -[1+(DATH-EDO)*<PRO|V>]^-1 (DATH-EDO) <PRO|U>        ==
!     ==================================================================
      DO I1=1,NPRO
        DVEC(I1)=0.D0
        DO I2=1,NPRO
          DVEC(I1)=DVEC(I1)+CMAT(I1,I2)*CVEC(I2)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  |PHI> = |U>+|V>DVEC                                         ==
!     ==  = [1-|V>[1+(DATH-EDO)*<PRO|V>]^-1(DATH-EDO)<PRO|] |U>       ==
!     ==================================================================
      PHI(:,:)=0.D0
      PHI(1:MCH,:)=U(1:MCH,:)
      DO I1=1,NPRO
        PHI(1:MCH,:)=PHI(1:MCH,:)+V(1:MCH,:,I1)*DVEC(I1)
      ENDDO
      RETURN
 9999 CONTINUE
      CALL ERROR$MSG('ERROR IN ESSL MATRIX INVERSION')
      CALL ERROR$STOP('PAWDER')
      STOP
      END
!
!     ..................................................................
      SUBROUTINE SCHROEDINHOM(R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI,GINH,DLG)
!     ******************************************************************
!     **                                                              **
!     **  SOLVES THE INHOMOGENEOUS, RELATIVISTIC RADIAL SCHROEDINGER  **
!     **  EQUATION                                                    **
!     **    [-0.5*GRAD^2+V(R)*Y0-E]|PHI>=|G>                          **
!     **  BOUNDARY CONDITION FOR INWARD INTEGRATION                   **
!     **    IS PHI(BOXRADIUS)=0                                       **
!     **                                                              **
!     **  METHOD:                                                     **
!     **    TRANSFORMES ONTO NEW VARIABLES U=PHI*R; X=1+LN(R/R1)/DEX  **
!     **    [D2/DX2 -DEX*D/DX+CF(R)] |U>=|CG(R)>                      **
!     **          CF=DEX^2 [ L(L+1) + 2R^2(V(R)*Y0-E) ]               **
!     **          CG=DEX^2 R^3 G(R)                                   **
!     **                                                              **
!     **    FOR R>TWICE THE CLASSICAL TURNING POINT                   **
!     **    AND OUTSIDE THE GRID, THE SEMICLASSICAL EQUATION IS USED  **
!     **                                                              **
!     **                                                              **
!     **  OPTIONS:                                                    **
!     **    TFRWRD-> OUTWARD/INWARD INTEGRATION                       **
!     **    THOM  -> INHOMOGENEITY IGNORED/NOT IGNORED                **
!     **                                                              **
!     ******************************************************************
      USE SCHRGL_INTERFACE_MODULE
      IMPLICIT NONE
      REAL(8),    PARAMETER    :: SCL=5.D0 ! SCL*R(CLASS.TURN.P.)=> SEMICLASICAL
      REAL(8),    INTENT(IN)   :: R1    ! RADIAL GRID R(I)=R1*EXP(DEX*(I-1))
      REAL(8),    INTENT(IN)   :: DEX   !
      INTEGER(4), INTENT(IN)   :: NR    ! NUMBER OF GRID POINTS
      INTEGER(4), INTENT(IN)   :: L     ! MAIN ANGULAR MOMENTUM
      REAL(8),    INTENT(IN)   :: AEZ   ! ATOMIC NUMBER (AEZ=0->NON-RELATIVISTIC)
      INTEGER(4), INTENT(IN)   :: MCH   ! INTEGRATION STOPS AT IR=MCH
      REAL(8),    INTENT(IN)   :: E     ! ENERGY
      REAL(8),    INTENT(IN)   :: POT(NR)    ! POTENTIAL
      REAL(8),    INTENT(IN)   :: GINH(NR)   !INHOMOGENEITY
      REAL(8),    INTENT(IN)   :: BOXRADIUS  !
      REAL(8),    INTENT(IN),OPTIONAL   :: DLG  ! LOGARITHMIC DERIVATIVE
      REAL(8),    INTENT(INOUT):: PHI(NR,3)  ! PHI,DPHI/DR,NABLA**2 PHI
      REAL(8)                  :: PHI1IN(NR,3)
      REAL(8)                  :: PHI1OUT(NR,3)
      REAL(8)                  :: PHI2IN(NR,3)
      REAL(8)                  :: PHI2OUT(NR,3)
      REAL(8)                  :: A(2,2)
      REAL(8)                  :: B(2)
      REAL(8)                  :: C(2)
      REAL(8)                  :: DET
      CALL SCHROEDER(.TRUE.,.TRUE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI1IN,GINH,DLG)
      CALL SCHROEDER(.TRUE.,.FALSE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI2IN,GINH,DLG)
      IF(PRESENT(DLG)) THEN
        CALL SCHROEDER(.FALSE.,.TRUE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI1OUT,GINH,DLG=DLG)
        CALL SCHROEDER(.FALSE.,.FALSE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI2OUT,GINH,DLG=DLG)
      ELSE
        CALL SCHROEDER(.FALSE.,.TRUE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI1OUT,GINH)
        CALL SCHROEDER(.FALSE.,.FALSE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI2OUT,GINH)
      END IF
      A(1,1)=PHI2IN(MCH,1)
      A(1,2)=-PHI2OUT(MCH,1)
      A(2,1)=PHI2IN(MCH,2)
      A(2,2)=-PHI2OUT(MCH,2)
      B(1)=PHI1OUT(MCH,1)-PHI1IN(MCH,1)
      B(2)=PHI1OUT(MCH,2)-PHI1IN(MCH,2)
      DET=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      C(1)=( A(2,2)*B(1)-A(1,2)*B(2) )/DET
      C(2)=(-A(2,1)*B(1)+A(1,1)*B(2) )/DET
      PHI1IN=PHI1IN+PHI2IN*C(1)
      PHI1OUT=PHI1OUT+PHI2OUT*C(2)
      PHI(1:MCH,:)=PHI1IN(1:MCH,:)
      PHI(MCH:NR,:)=PHI1OUT(MCH:NR,:)
      RETURN 
      END
!     ............
      SUBROUTINE WRITIT(R1,DEX,NR,TEXT,FUNC)
      CHARACTER(LEN=*)  :: TEXT
      REAL(8)     :: R1  
      REAL(8)     :: DEX  
      INTEGER(4)  :: NR
      REAL(8)     :: FUNC(NR)
      REAL(8)     :: XEXP,RI
      INTEGER(4)  :: IR
      INTEGER(4)  :: ID=5
      WRITE(*,*)'WRITIT:',TEXT
      XEXP=DEXP(DEX)
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        IF(MOD(IR,ID).EQ.0)WRITE(*,FMT='(2F10.5)')RI,FUNC(IR)
      ENDDO
      STOP 'FORCED STOP'
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE FIXNODE(R1,DEX,NR,L,E,AEZ,BOXRADIUS,POT,PHI,ZEFF,NN,TOL,DLG)
!     **                                                              **
!     **                                                              **
!     **                                                              **
      USE SCHRGL_INTERFACE_MODULE,ONLY:AEBOUNDSTATE,SCHROEDER
      IMPLICIT NONE
      logical(4),parameter    :: tpr=.false.
      REAL(8),   PARAMETER    :: ITERX=200
      REAL(8),   INTENT(IN)   :: R1
      REAL(8),   INTENT(IN)   :: DEX
      INTEGER(4),INTENT(IN)   :: NR
      INTEGER(4),INTENT(IN)   :: L      ! main angular momentum QN
      INTEGER(4),INTENT(IN)   :: NN     ! target #(nodes) 
      REAL(8),   INTENT(IN)   :: AEZ    ! atomic number
      REAL(8),   INTENT(INOUT):: E      ! band energy
      REAL(8),   INTENT(IN)   :: BOXRADIUS
      REAL(8),   INTENT(IN), OPTIONAL :: DLG
      REAL(8),   INTENT(IN)   :: TOL
      REAL(8),   INTENT(IN)   :: POT(NR)
      REAL(8),   INTENT(IN)   :: ZEFF
      REAL(8),   INTENT(OUT)  :: PHI(NR,3)
      LOGICAL                 :: CONVG
      REAL(8)                 :: EOLD
      INTEGER(4)              :: NNODE
      INTEGER(4)              :: Nnext
      INTEGER(4)              :: ITER,IR
      real(8)                 :: eup,elow
      real(8)                 :: ginh(nr)
      INTEGER(4)              :: IROUT
      real(8)                 :: rout,VAL,DER,LOGDER
      real(8)                 :: pi,svar,rcl
!     ******************************************************************
      pi=4.d0*datan(1.d0)
!
!     ==================================================================
!     == ESTIMATE ENERGY OF THE BOUND STATE BY SHOOTING OUTWARD       ==
!     == AND BISECTION OF THE ENERGY WINDOW                           ==
!     == THE BOUND STATE HAS A VANISHING LOGARITHMIC DERIVATIVE       ==
!     == THE RADIUS IS CHOSEN DYNAMICALLY: FOR CORE STATES IT IS      ==
!     == RELATED  TO THE CLASSICAL TURNING POINT.                     ==
!     ==================================================================
      EUP=10.D0
      ELOW=-1.D+4
      GINH(:)=0.D0
      ITER=0
      CONVG=.FALSE.
      DO WHILE(.NOT.CONVG)
        ITER=ITER+1
        E=0.5D0*(EUP+ELOW)
!       == ESTIMATE CLASSICAL TURNING POINT FOR LOW LYING SHELLS
        IROUT=1+INT(LOG(5.D0/R1)/DEX)
        DO IR=21,IROUT
          IF(POT(IR-20).GT.E) THEN
            IROUT=IR
            EXIT
          END IF
        ENDDO    
        ROUT=R1*EXP(DEX*REAL(IROUT-1,KIND=8))
        CALL SCHROEDER(.TRUE.,.TRUE.,R1,DEX,NR,L,E,AEZ,IROUT,ROUT &
     &                     ,POT,PHI,GINH)
        NNODE=0
        DO IR=1,IROUT-1
          IF(PHI(IR,1)*PHI(IR+1,1).LT.0.D0)NNODE=NNODE+1
        ENDDO
        VAL=PHI(IROUT,1)
        DER=PHI(IROUT,2)
        LOGDER=DATAN(-ROUT*DER/VAL)/PI+REAL(NNODE)+0.5D0
        SVAR=LOGDER-(REAL(NN,KIND=8)+0.5D0)
        IF(SVAR.GT.0.D0) THEN
          EUP=E
        ELSE 
          ELOW=E
        END IF
        if(abs(eup-elow).gt.1.d-2) cycle
        if(abs(eup-elow).lt.1.d-10) exit
        CONVG=(ABS(SVAR).LT.1.D-3)
        IF(ITER.GT.ITERX) THEN
          CALL ERROR$MSG('SHOOTING LOOP NOT CONVERGED')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$R8VAL('E',E)
          CALL ERROR$I4VAL('TARGET #(NODES)',NN)
          CALL ERROR$I4VAL('ACTUAL #(NODES)',NNODE)
          CALL ERROR$STOP('FIXNODE')
        END IF
      ENDDO
!print*,'logder ',logder,e,elow,eup
!
!     ==================================================================
!     == optimize energy by energy linearization                      ==
!     ==================================================================
      ITER=0
      CONVG=.FALSE.
      DO WHILE(.NOT.CONVG)
        ITER=ITER+1
        EOLD=E
        IF(PRESENT(DLG)) THEN
          CALL AEBOUNDSTATE(R1,DEX,NR,L,E,AEZ,BOXRADIUS,POT,PHI,DLG=DLG)
        ELSE
          CALL AEBOUNDSTATE(R1,DEX,NR,L,E,AEZ,BOXRADIUS,POT,PHI)
        END IF
        NNODE=0
        DO IR=1,NR-1
          IF(PHI(IR,1)*PHI(IR+1,1).LT.0.D0)NNODE=NNODE+1
        ENDDO
!print*,'e ',e,nnode
        CONVG=ABS(EOLD-E).LT.1.D-9
        IF(ITER.GT.ITERX) THEN
          PRINT*,'PRESENT(DLG)=',PRESENT(DLG)
          CALL ERROR$MSG('LOOP NOT CONVERGED')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$R8VAL('E',E)
          CALL ERROR$I4VAL('TARGET #(NODES)',NN)
          CALL ERROR$I4VAL('ACTUAL #(NODES)',NNODE)
          CALL ERROR$STOP('FIXNODE')
        END IF
      ENDDO
!
!     ==========================================================
!     == TOLERATE INCORRECT #NODES FOR WEAKLY BOUND STATES    ==
!     ==========================================================
      IF(NNODE.NE.NN)  THEN
        IF(E.LT.-0.1D0) THEN
          CALL ERROR$MSG('#NODES INCORRECT')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$R8VAL('E',E)
          CALL ERROR$I4VAL('TARGET #(NODES)',NN)
          CALL ERROR$I4VAL('ACTUAL #(NODES)',NNODE)
          CALL ERROR$STOP('FIXNODE')
        ELSE
          PRINT*,'WARNING: #NODES INCORRECT L,NN,E ',L,NN,E
        END IF
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AEBOUNDSTATE(R1,DEX,NR,L,E,AEZ,BOXRADIUS,POT,PHI,DLG)
!     **                                                              **
!     **                                                              **
!     **                                                              **
      USE SCHRGL_INTERFACE_MODULE, ONLY : SCHROEDER
      IMPLICIT NONE
      REAL(8),    INTENT(IN)   :: R1    ! RADIAL GRID R(I)=R1*EXP(DEX*(I-1))
      REAL(8),    INTENT(IN)   :: DEX   !
      INTEGER(4), INTENT(IN)   :: NR    ! NUMBER OF GRID POINTS
      INTEGER(4), INTENT(IN)   :: L     ! MAIN ANGULAR MOMENTUM
      REAL(8),    INTENT(INOUT):: E     ! ENERGY
      REAL(8),    INTENT(IN)   :: AEZ   ! ATOMIC NUMBER (AEZ=0->NON-RELATIVISTIC)
      REAL(8),    INTENT(IN)   :: POT(NR)      ! POTENTIAL
      REAL(8),    INTENT(OUT)  :: PHI(NR,3)    ! BOUND STATE, DERIVATIVE, LAPLACE
      REAL(8),    INTENT(IN)   :: BOXRADIUS    ! BOUNDARY COND.: PHI(BOXRADIUS)=0
      REAL(8),    INTENT(IN),OPTIONAL  :: DLG  ! LOGARITHMIC DERIVATIVE
      REAL(8)                  :: PHI1(NR,3)   ! PARTIAL WAVE OUTWARD INTEGR.
      REAL(8)                  :: PHI2(NR,3)   ! PARTIAL WAVE INWARD INTEGR
      REAL(8)                  :: PHI1DOT(NR,3)! ENERGY DERIVATIVE OF PHI1
      REAL(8)                  :: PHI2DOT(NR,3)! ENERGY DERIVATIVE OF PHI2
      REAL(8)                  :: GINH(NR)     ! INHOMOGENEITY, ALSO USED AS DUMMY
      INTEGER(4)               :: MCH          ! MATCHING GRID POINT
      REAL(8)                  :: SVAR,SVAR1,SVAR2,SVAR3,DE,STEP(3)
      INTEGER(4)               :: IR,I
      REAL(8)                  :: XEXP,RI,Y0,PI
      REAL(8)                  :: R(NR)
      REAL(8)                  :: phi1x(3),phi1dotx(3)
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
      XEXP=DEXP(DEX)
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        R(IR)=RI
      ENDDO
!
!     ==================================================================
!     == FIND CLASSICAL TURNING POINT                                 ==
!     ==================================================================
      MCH=4
      DO IR=4,NR
        IF(POT(IR)*Y0.LT.E)MCH=IR
      ENDDO
      MCH=MIN(NR-4,MCH)
      if(present(dlg)) then
        MCH=MIN(nr-4,INT(1.D0+DLOG(BOXRADIUS/R1)/DEX)+2)
      ELSE
        MCH=MIN(MCH,INT(1.D0+DLOG(BOXRADIUS/R1)/DEX)-4)
      END IF
!
!     ==================================================================
!     == evaluate phi                                                 ==
!     ==================================================================
!     __ INTEGRATE PHI OUTWARD
      CALL SCHROEDER(.TRUE.,.TRUE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI1,GINH)
      SVAR=1.D0/PHI1(MCH,1)
      PHI1(:,:)=PHI1(:,:)*SVAR
!     __ INTEGRATE PHI INWARD
      IF(PRESENT(DLG)) THEN
        CALL SCHROEDER(.FALSE.,.TRUE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI2,GINH,DLG=DLG)
!       print*,'dlg   hom',dlg,phi2(mch,2)/phi2(mch,1)
      ELSE
        CALL SCHROEDER(.FALSE.,.TRUE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI2,GINH)
      END IF
      IF(DABS(PHI2(MCH,1)).LT.1.D-6) THEN
        PRINT*,'PHI2 IS ZERO: PHI2=',PHI2(MCH,:)
        STOP
      END IF
      SVAR=1.D0/PHI2(MCH,1)
      PHI2(:,:)=PHI2(:,:)*SVAR
!
!     ==================================================================
!     == Evaluate phidot through the inhomogeneous schROEDINGER EQUATION 
!     ==================================================================
!
!     __ INTEGRATE PHIDOT OUTWARD
      GINH=PHI1(:,1)
      CALL SCHROEDER(.TRUE.,.FALSE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI1DOT,GINH)
!     __ INTEGRATE INWARD
      GINH=PHI2(:,1)
      IF(PRESENT(DLG)) THEN
        CALL SCHROEDER(.FALSE.,.FALSE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI2DOT,GINH,DLG=DLG)
!       print*,'dlg inhom',dlg,phi2dot(mch,2)/phi2dot(mch,1)
      ELSE
        CALL SCHROEDER(.FALSE.,.FALSE.,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI2DOT,GINH)
      END IF
      SVAR=(PHI1DOT(MCH,1)-PHI2DOT(MCH,1))/PHI2(MCH,1)
      PHI2DOT(:,:)=PHI2DOT(:,:)+SVAR*PHI2(:,:)
!
!     ==================================================================
!     == MATCH DERIVATIVES                                            ==
!     ==================================================================
      IF(PRESENT(DLG)) THEN
        CALL OLDRADIAL$VALUE(R1,DEX,NR,PHI1(:,1),boxradius,PHI1x(1))
        CALL OLDRADIAL$VALUE(R1,DEX,NR,PHI1(:,2),boxradius,PHI1x(2))
        CALL OLDRADIAL$VALUE(R1,DEX,NR,PHI1dot(:,1),boxradius,PHI1dotx(1))
        CALL OLDRADIAL$VALUE(R1,DEX,NR,PHI1dot(:,2),boxradius,PHI1dotx(2))
        DE=-(PHI1X(2)-DLG/BOXRADIUS*PHI1X(1))/(PHI1DOTX(2)-DLG/BOXRADIUS*PHI1DOTX(1))
        PHI(1:MCH,:)=PHI1(1:MCH,:)+PHI1DOT(1:MCH,:)*DE
      ELSE
        DE=(PHI1(MCH,2)-PHI2(MCH,2))/(PHI2DOT(MCH,2)-PHI1DOT(MCH,2))
        PHI(1:MCH,:)=PHI1(1:MCH,:)+PHI1DOT(1:MCH,:)*DE
        PHI(MCH:NR,:)=PHI2(MCH:NR,:)+PHI2DOT(MCH:NR,:)*DE
        STEP(:)=PHI1(MCH,:)+PHI1DOT(MCH,:)*DE &
     &       -PHI2(MCH,:)-PHI2DOT(MCH,:)*DE
      END IF
!
!     ==================================================================
!     == NORMALIZE                                                    ==
!     ==================================================================
      CALL OLDRADIAL$INTEGRAL(R1,DEX,NR,(R(:)*PHI(:,1))**2,SVAR)
      SVAR=1.D0/DSQRT(SVAR)
      PHI(:,:)=PHI(:,:)*SVAR
      STEP(:)=STEP(:)*SVAR
!     print*,'step ',step/phi(mch,1)
!
!     ==================================================================
!     == EVALUATE NEW ENERGY                                          ==
!     ==================================================================
      GINH(:)=PHI(:,1)*(-0.5D0*PHI(:,3)+POT(:)*Y0*PHI(:,1))
      CALL OLDRADIAL$INTEGRAL(R1,DEX,NR,R(:)**2*GINH(:),SVAR)
      E=SVAR
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SCHROEDER(TFRWRD,THOM,R1,DEX,NR,L,E,AEZ,MCH,BOXRADIUS &
     &                     ,POT,PHI,GINH,DLG)
!     ******************************************************************
!     **                                                              **
!     **  SOLVES THE INHOMOGENEOUS, scalar RELATIVISTIC RADIAL        **
!     **  SCHROEDINGER EQUATION                                       **
!     **    [-0.5*GRAD^2+V(R)*Y0-E]|PHI>=|G>                          **
!     **  BOUNDARY CONDITION FOR INWARD INTEGRATION                   **
!     **    IS PHI(BOXRADIUS)=0                                       **
!     **                                                              **
!     **  METHOD:                                                     **
!     **    TRANSFORMES ONTO NEW VARIABLES U=PHI*R; X=1+LN(R/R1)/DEX  **
!     **    [D2/DX2 -DEX*D/DX+CF(R)] |U>=|CG(R)>                      **
!     **          CF=DEX^2 [ L(L+1) + 2R^2(V(R)*Y0-E) ]               **
!     **          CG=DEX^2 R^3 G(R)                                   **
!     **                                                              **
!     **    FOR R>TWICE THE CLASSICAL TURNING POINT                   **
!     **    AND OUTSIDE THE GRID, THE SEMICLASSICAL EQUATION IS USED  **
!     **                                                              **
!     **                                                              **
!     **  OPTIONS:                                                    **
!     **    TFRWRD-> OUTWARD/INWARD INTEGRATION                       **
!     **    THOM  -> INHOMOGENEITY IGNORED/NOT IGNORED                **
!     **                                                              **
!     **  REMARK: THIS ROUTINE SOMETIMES GAVE ZERO RESULT FOR INWARD  **
!     **  INTEGRATION WITHOUT INHOMOGENEITY                           **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      USE SCHRGL_INTERFACE_MODULE, ONLY : OUTBOUNDARY
      IMPLICIT NONE
      REAL(8),    PARAMETER    :: SCL=5.D0 ! SCL*R(CLASS.TURN.P.)=> SEMICLASICAL
      LOGICAL,   INTENT(IN)    :: THOM  ! IGNORE INHOMOGENITY?   
      LOGICAL,   INTENT(IN)    :: TFRWRD !INTEGRATE OUTWARD/INWARD
      REAL(8),    INTENT(IN)   :: R1    ! RADIAL GRID R(I)=R1*EXP(DEX*(I-1))
      REAL(8),    INTENT(IN)   :: DEX   !
      INTEGER(4), INTENT(IN)   :: NR    ! NUMBER OF GRID POINTS
      INTEGER(4), INTENT(IN)   :: L     ! MAIN ANGULAR MOMENTUM
      REAL(8),    INTENT(IN)   :: AEZ   ! ATOMIC NUMBER (AEZ=0->NON-RELATIVISTIC)
      INTEGER(4), INTENT(IN)   :: MCH   ! INTEGRATION STOPS AT IR=MCH
      REAL(8),    INTENT(IN)   :: E     ! ENERGY
      REAL(8),    INTENT(IN)   :: POT(NR)    ! POTENTIAL
      REAL(8),    INTENT(IN)   :: GINH(NR)   !INHOMOGENEITY
      REAL(8),    INTENT(IN)   :: BOXRADIUS  !
      REAL(8),    INTENT(IN),OPTIONAL   :: DLG  ! LOGARITHMIC DERIVATIVE
      REAL(8),    INTENT(OUT)  :: PHI(NR,3)  ! PHI,DPHI/DR,NABLA**2 PHI
      REAL(8)                  :: RI
      REAL(8)                  :: XEXP
      INTEGER(4)               :: IR
      REAL(8)                  :: R(NR)    ! RADIAL GRID
      REAL(8)                  :: CF(NR)   ! POTENTIAL COEFF. ARRAY
                               ! CF=DEX**2*( L*(L+1) + 2*(V-E)*R**2 )
      REAL(8)                  :: CG(NR)   ! INHOMOGENITY COEFF. ARRAY
                               ! CG=2*DEX**2*R**3 *GINH
      REAL(8)                  :: FR(NR)   ! RELATIVISTIC CORRECTION
      REAL(8)                  :: FRP(NR)  ! RELATIVISTIC CORRECTION
      REAL(8)                  :: U(NR)    ! R*PHI
      REAL(8)                  :: UP(NR)   ! DU/DX  X(R)=1+LN(R/R1)/DEX
      REAL(8)                  :: UPP(NR)  ! D2U/DX2
      INTEGER(4)               :: N1,N2
      REAL(8)                  :: PI,Y0
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
!
!     ==================================================================
!     ==  RADIAL GRID                                                 ==
!     ==================================================================
      XEXP=DEXP(DEX)
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        R(IR)=RI
      ENDDO
!
!     ==================================================================
!     ==  COEFFICIENT ARRAYS                                          ==
!     ==================================================================
      CALL VTRANSFORM(DEX,NR,AEZ,L,E,R,POT,CF,FR,FRP)
      IF(THOM) THEN
        DO IR=1,NR
          CG(IR)=0.D0
        ENDDO
      ELSE
        CALL GTRANSFORM(DEX,NR,R,GINH,CG)
      END IF
!
!     ==================================================================
!     ==  CALCULATE U=R*PHI                                           ==
!     ==================================================================
      u(:)=0.d0
      up(:)=0.d0
      upp(:)=0.d0
!     DO IR=1,NR
!       RI=R(IR)
!       U(IR)=RI*PHI(IR,1)
!       UP(IR)=DEX*(U(IR)+RI**2*PHI(IR,2))
!       UPP(IR)=DEX*(UP(IR)+2.D0*DEX*RI**2*PHI(IR,2)+DEX*PHI(IR,3))
!     ENDDO
!
!     ==================================================================
!     ==  RADIAL GRID                                                 ==
!     ==================================================================
      IF(TFRWRD) THEN
        N1=4
        N2=MCH
        CALL INBOUNDARY(THOM,AEZ,L,DEX,NR,R,U,UP,UPP,CF,CG,FR,FRP)
        CALL OUTWRD(DEX,NR,N1,N2,U,UP,UPP,CF,CG,FR,FRP)
      ELSE
        N1=MIN(INT(1.D0+LOG(BOXRADIUS/R1)/DEX)-4 &
     &         ,NR-4 &
     &         ,INT(MCH+LOG(SCL)/DEX))
        N2=MCH
        IF(PRESENT(DLG)) THEN
          N1=MIN(INT(1.D0+LOG(BOXRADIUS/R1)/DEX),NR-4)
          CALL OUTBOUNDARY(THOM,R1,DEX,NR,R,N1,BOXRADIUS,U,UP,UPP,CF,DLG=DLG)
          CALL INWRD(DEX,NR,N1,N2,U,UP,UPP,CF,CG,FR,FRP)
        ELSE
          N1=MIN(INT(1.D0+LOG(BOXRADIUS/R1)/DEX) &
     &         ,NR-4 &
     &         ,INT(MCH+LOG(SCL)/DEX))
          CALL OUTBOUNDARY(THOM,R1,DEX,NR,R,N1,BOXRADIUS,U,UP,UPP,CF)
          CALL INWRD(DEX,NR,N1,N2,U,UP,UPP,CF,CG,FR,FRP)
        END IF
      END IF
!     PRINT*,'SCHROEDER MARKE 2'
!
!     ==================================================================
!     ==  CALCULATE PHI(:,1)=U/R; PHI(:,2)=DPHI/DR;                   ==
!     ==================================================================
      IF(TFRWRD) THEN
        N1=1 ;N2=MCH
        PHI(MCH+1:NR,:)=0.D0
      ELSE
        IF(PRESENT(DLG)) THEN
          N1=MCH
          N2=NR
        ELSE
          N1=MCH; N2=NR
          PHI(1:MCH-1,:)=0.D0
        END IF
      END IF
      DO IR=N1,N2
        RI=R(IR)
        PHI(IR,1)=U(IR)/RI
        PHI(IR,2)=(UP(IR)/DEX-U(IR))/RI**2
        PHI(IR,3)=2.D0*(POT(IR)*Y0-E)*PHI(IR,1)
      ENDDO
      IF(.NOT.THOM) THEN
        DO IR=N1,N2
          PHI(IR,3)=PHI(IR,3)-2.D0*GINH(IR)
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE VTRANSFORM(DEX,NR,AEZ,L,E,R,POT,CF,FR,FRP)
!     **                                                              **
!     **  COEFFICENT ARRAY FOR THE POTENTIAL                          **
!     **  AND RELATIVISTIC CORRECTIONS                                **
!     **                                                              **
      IMPLICIT NONE
      REAL(8),  PARAMETER  :: FSS=(1.0D0/137.036D0)**2
      REAL(8),  INTENT(IN) :: DEX 
      INTEGER, INTENT(IN) :: NR      ! #(RADIAL GRID POINTS)
      REAL(8),  INTENT(IN) :: AEZ     ! ATOMIC NUMBER (=0 FOR NON-RELATIVISTIC)
      INTEGER, INTENT(IN) :: L       ! MAIN ANGULAR MOMENTUM
      REAL(8),  INTENT(IN) :: E       ! ONE PARTICLE ENERGY
      REAL(8),  INTENT(IN) :: R(NR)   ! RADIAL GRID
      REAL(8),  INTENT(IN) :: POT(NR) ! POTENTIAL
      REAL(8),  INTENT(OUT):: CF(NR)  ! 
      REAL(8),  INTENT(OUT):: FR(NR)  ! RELATIVISTIC CORRECTION
      REAL(8),  INTENT(OUT):: FRP(NR) ! RELATIVISTIC CORRECTION
      REAL(8)              :: SLS     ! L(L+1)
      REAL(8)              :: ALS     ! DEX**2
      REAL(8)              :: PI
      REAL(8)              :: Y0      ! SPHERICAL HARMONIC (L=0)
      REAL(8)              :: GAMMA
      REAL(8)              :: DV(2)
      REAL(8)              :: DVI
      REAL(8)              :: SVAR
      INTEGER(4)           :: IR
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
!
!     ==================================================================
!     ==  COEFFICIENT ARRAY FOR U IN DIFFERENTIAL EQ.                 ==
!     ==================================================================
      SLS=DBLE(L*(L+1))
      ALS=DEX**2
      DO IR=1,NR
        CF(IR)=ALS*( SLS+2.0D0*(Y0*POT(IR)-E)*R(IR)**2 )
      ENDDO
!
!     ==================================================================
!     ==  RELATIVISTIC - NON-RELATIVISTIC SWITCH                      ==
!     ==================================================================
      IF(AEZ.LT.1.D-5) THEN
        DO IR=1,NR
          FR(IR)=0.D0
          FRP(IR)=0.D0
        ENDDO
        RETURN
      END IF
!
!     ==================================================================
!     ==  RELATIVISTIC COEFFICIENT ARRAYS FOR U (FR) AND UP (FRP).    ==
!     ==  CALCULATE DV/DR FOR DARWIN CORRECTION                       ==
!     ==================================================================
      IF(L .EQ. 0) THEN
        GAMMA=DSQRT(1.0D0-FSS*AEZ**2)
      ELSE
        GAMMA=(L*DSQRT(L**2-FSS*AEZ**2) & 
     &   + (L+1)*DSQRT((L+1)**2-FSS*AEZ**2))/(2*L+1)
      END IF
      DV(1)=(-50.D0*Y0*POT(1)+96.D0*Y0*POT(2)-72.D0*Y0*POT(3) &
     &       +32.D0*Y0*POT(4)- 6.D0*Y0*POT(5))/(24.D0*DEX*R(1))
      DV(2)=( -6.D0*Y0*POT(1)-20.D0*Y0*POT(2)+36.D0*Y0*POT(3) &
     &       -12.D0*Y0*POT(4)+ 2.D0*Y0*POT(5))/(24.D0*DEX*R(2))
      DO IR=1,2
        SVAR=0.5D0*FSS*DV(IR)*R(IR)/(1.0D0+0.5D0*FSS*(E-Y0*POT(IR)))
        FR(IR)=ALS*( SVAR-FSS*( R(IR)*(Y0*POT(IR)-E) )**2 )
        FRP(IR)=-DEX*SVAR
      ENDDO
      DO IR=3,NR-2
!       __GRAD(POT(IR))__5-POINT FORMULA________________________________
        DVI=(2.D0*Y0*POT(IR-2)-16.D0*Y0*POT(IR-1)+16.D0*Y0*POT(IR+1) &
     &        -2.D0*Y0*POT(IR+2))/(24.D0*DEX*R(IR))
        SVAR=0.5D0*FSS*DVI*R(IR)/(1.0D0+0.5D0*FSS*(E-Y0*POT(IR)))
        FR(IR)=ALS*( SVAR-FSS*( R(IR)*(Y0*POT(IR)-E) )**2 )
        FRP(IR)=-DEX*SVAR
      ENDDO
      FR(NR-1) =0.D0
      FRP(NR-1)=0.D0
      FR(NR)   =0.D0
      FRP(NR)  =0.D0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GTRANSFORM(DEX,NR,R,GINH,CG)
!     **                                                              **
!     **  COEFFICIENT ARRAY FOR THE INHOMOGENEITY                     **
!     **    CG=-2 DEX^2 R^3 GINH                                      **
!     **                                                              **
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NR  
      REAL(8),  INTENT(IN) :: DEX 
      REAL(8),  INTENT(IN) :: GINH(NR)  ! INHOMOGENEITY
      REAL(8),  INTENT(IN) :: R(NR) ! RADIAL GRID
      REAL(8),  INTENT(OUT):: CG(NR)
      REAL(8)              :: SVAR
      INTEGER(4)           :: IR
!
!     ==================================================================
!     ==  COEFFICIENT ARRAY FOR U IN DIFFERENTIAL EQ.                 ==
!     ==================================================================
      SVAR=-2.D0*DEX**2
      DO IR=1,NR
        CG(IR)=SVAR*R(IR)**3*GINH(IR)
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE INBOUNDARY(THOM,AEZ,L,DEX,NR,R,U,UP,UPP,CF,CG,FR,FRP)
!     **                                                              **
!     **                                                              **
      IMPLICIT NONE
      REAL(8),   PARAMETER    :: FSS=(1.0D0/137.036D0)**2
      REAL(8),   INTENT(IN)   :: DEX
      REAL(8),   INTENT(IN)   :: AEZ
      INTEGER(4),INTENT(IN)   :: NR
      INTEGER(4),INTENT(IN)   :: L ! MAIN ANGULAR MOMENTUM
      REAL(8),   INTENT(IN)   :: R(NR)
      REAL(8),   INTENT(INOUT):: U(NR)
      REAL(8),   INTENT(INOUT):: UP(NR)
      REAL(8),   INTENT(INOUT):: UPP(NR)
      REAL(8),   INTENT(IN)   :: CF(NR)
      REAL(8),   INTENT(IN)   :: CG(NR)
      REAL(8),   INTENT(IN)   :: FR(NR)
      REAL(8),   INTENT(IN)   :: FRP(NR)
      LOGICAL,  INTENT(IN)   :: THOM
      REAL(8)                 :: GAMMA
      REAL(8)                 :: RIGAMMA
      INTEGER(4)              :: IR, it !cwk
      REAL(8)                 :: U0(4),UP0(4), vtest !cwk
!     ******************************************************************
      IF(.NOT.THOM) THEN
         DO IR=1,4
            U(IR)=0.D0
            UP(IR)=0.D0
            UPP(IR)=CG(IR)
         ENDDO
         RETURN
      END IF
      IF(AEZ.GT.0) THEN
        gamma = sqrt(1.0+l*(l+1)-FSS*AEZ**2)
      ELSE
         GAMMA=L+1
      END IF
!cwk my new  Attention!! U() means now u/r**gamma
      DO IR=1,4
        U(IR)=1.0
!   UP(IR)=2.0*((gamma-1.0)/(fss*aez))*(gamma+2.0)/(2.0*gamma+1.0)*dex*r(ir)
        UP(IR)=0.d0
!       print*,up(ir),-Aez*dex*r(ir) ! to jeszcze sprawdzic
        UPP(IR)=(DEX+FRP(IR)-2.0*dex*gamma)*UP(IR) &
     &         +(CF(IR)+FR(IR)+dex*gamma*(DEX+FRP(IR))-(dex*gamma)**2)*U(IR) &
     &         +CG(IR)/r(ir)**gamma
      ENDDO
      vtest = 1.0d2
      it=0
      do while (vtest.gt.1.e-8)
        it = it+1
        u0=u(1:4)
        up0 = up(1:4)
        DO IR=2,4
          U(IR)=U(ir-1)+0.5*(UP(ir)+up(ir-1))
!print*,'marke 1b',ir,up(ir-1),upp(ir),upp(ir-1)
          UP(IR)=Up(ir-1)+0.5*(UPp(ir)+upp(ir-1))
!print*,'marke 1c'
          UPP(IR)=(DEX+FRP(IR)-2.0*dex*gamma)*UP(IR) &
    &         +(CF(IR)+FR(IR)+dex*gamma*(DEX+FRP(IR))-(dex*gamma)**2)*U(IR) &
    &         +CG(IR)/r(ir)**gamma
        ENDDO
        vtest=abs((u(4)-u0(4))/u0(4))+abs((up(4)-up0(4))/up0(4))
        if(it.gt.500) then
          call error$msg('iteration not converged')
          call error$stop('inboundary')
        end if
      end do   
!     print *,'test in inbound', it, vtest !cwk
      do ir = 1,4
        RIGAMMA=R(ir)**GAMMA
!       print*,'bef',u(ir),up(ir),rigamma
!cwk here we are again at old U()
        up(ir)=(up(ir)+dex*u(ir)*gamma)*rigamma
        u(ir)=u(ir)*rigamma
        UPP(IR)=(DEX+FRP(IR))*UP(IR)+(CF(IR)+FR(IR))*U(IR)+CG(IR)
!       print '(2x, 4f20.8)',u(ir),up(ir),upp(ir)
!!$     print *,(DEX+FRP(IR))*UP(IR),(CF(IR)+FR(IR))*U(IR)
      end do   
      RETURN

      END
!
!     ..................................................................
      SUBROUTINE INBOUNDARY_old1(THOM,AEZ,L,DEX,NR,R,U,UP,UPP,CF,CG,FR,FRP)
!     **                                                              **
!     **                                                              **
      IMPLICIT NONE
      REAL(8),   PARAMETER    :: FSS=(1.0D0/137.036D0)**2
      REAL(8),   INTENT(IN)   :: DEX
      REAL(8),   INTENT(IN)   :: AEZ
      INTEGER(4),INTENT(IN)   :: NR
      INTEGER(4),INTENT(IN)   :: L ! MAIN ANGULAR MOMENTUM
      REAL(8),   INTENT(IN)   :: R(NR)
      REAL(8),   INTENT(INOUT):: U(NR)
      REAL(8),   INTENT(INOUT):: UP(NR)
      REAL(8),   INTENT(INOUT):: UPP(NR)
      REAL(8),   INTENT(IN)   :: CF(NR)
      REAL(8),   INTENT(IN)   :: CG(NR)
      REAL(8),   INTENT(IN)   :: FR(NR)
      REAL(8),   INTENT(IN)   :: FRP(NR)
      LOGICAL,  INTENT(IN)   :: THOM
      REAL(8)                 :: GAMMA
      REAL(8)                 :: RIGAMMA
      INTEGER(4)              :: IR, it !cwk
      REAL(8)                 :: U0(4),UP0(4), vtest !cwk
!     ******************************************************************
      IF(.NOT.THOM) THEN
         DO IR=1,4
            U(IR)=0.D0
            UP(IR)=0.D0
            UPP(IR)=CG(IR)
         ENDDO
         RETURN
      END IF
      IF(AEZ.GT.0) THEN
         IF(L .EQ. 0) THEN
            GAMMA=DSQRT(1.0D0-FSS*AEZ**2)
         ELSE
!$$$            GAMMA=(    L*DSQRT(   L**2 -FSS*AEZ**2)+ &
!$$$     &           (L+1)*DSQRT((L+1)**2-FSS*AEZ**2)   )/(2*L+1)
            gamma = sqrt(1.0+l*(l+1)-FSS*AEZ**2)
!cwk I have changed for gamma as defined in Koelling-Harmon
         END IF
      ELSE
         GAMMA=L+1
      END IF
! cwk your old machine
      
!!$      DO IR=1,4
!!$         RIGAMMA=R(ir)**GAMMA
!!$         U(IR)=RIGAMMA
!!$         UP(IR)=DEX*GAMMA*RIGAMMA
!!$         UPP(IR)=(DEX+FRP(IR))*UP(IR)+(CF(IR)+FR(IR))*U(IR)+CG(IR)
!!$      ENDDO
!cwk my new  Attention!! U() means now u/r**gamma
      DO IR=1,4
        U(IR)=1.0
        UP(IR)=2.0*((gamma-1.0)/(fss*aez))*(gamma+2.0)/(2.0*gamma+1.0)*dex*r(ir)
!       print*,up(ir),-Aez*dex*r(ir) ! to jeszcze sprawdzic
        UPP(IR)=(DEX+FRP(IR)-2.0*dex*gamma)*UP(IR) &
     &          +(CF(IR)+FR(IR)+dex*gamma*(DEX+FRP(IR)) &
                  -(dex*gamma)**2)*U(IR)+CG(IR)
      ENDDO
      u0=u(1:4)
      up0 = up(1:4)
      vtest = 1.0d2
      it=0
      do while (vtest .gt. 1.e-8)
        DO IR=2,4
          U(IR)=U(ir-1)+0.5*(UP(ir)+up(ir-1))
          UP(IR)=Up(ir-1)+0.5*(UPp(ir)+upp(ir-1))
          UPP(IR)=(DEX+FRP(IR)-2.0*dex*gamma)*UP(IR) &
    &            +(CF(IR)+FR(IR)+dex*gamma*(DEX+FRP(IR)) &
    &              -(dex*gamma)**2)*U(IR)+CG(IR)
        ENDDO
        vtest=abs((u(4)-u0(4))/u0(4))+abs((up(4)-up0(4))/up0(4))
        u0=u(1:4)
        up0=up(1:4)
        it = it+1
        if (it .gt. 500 ) exit
      end do   
!     print *,'test in inbound', it, vtest !cwk
      do ir = 1,4
        RIGAMMA=R(ir)**GAMMA
!       print*,'bef',u(ir),up(ir),rigamma
!cwk here we are again at old U()
        up(ir)=(up(ir)+dex*u(ir)*gamma)*rigamma
        u(ir)=u(ir)*rigamma
        UPP(IR)=(DEX+FRP(IR))*UP(IR)+(CF(IR)+FR(IR))*U(IR)+CG(IR)
!       print '(2x, 4f20.8)',u(ir),up(ir),upp(ir)
!!$     print *,(DEX+FRP(IR))*UP(IR),(CF(IR)+FR(IR))*U(IR)
      end do   
      RETURN

      END
!
!     ..................................................................
      SUBROUTINE INBOUNDARY_old(THOM,AEZ,L,DEX,NR,R,U,UP,UPP,CF,CG,FR,FRP)
!     **                                                              **
!     **                                                              **
      IMPLICIT NONE
      REAL(8),   PARAMETER    :: FSS=(1.0D0/137.036D0)**2
      REAL(8),   INTENT(IN)   :: DEX
      REAL(8),   INTENT(IN)   :: AEZ
      INTEGER(4),INTENT(IN)   :: NR
      INTEGER(4),INTENT(IN)   :: L ! MAIN ANGULAR MOMENTUM
      REAL(8),   INTENT(IN)   :: R(NR)
      REAL(8),   INTENT(INOUT):: U(NR)
      REAL(8),   INTENT(INOUT):: UP(NR)
      REAL(8),   INTENT(INOUT):: UPP(NR)
      REAL(8),   INTENT(IN)   :: CF(NR)
      REAL(8),   INTENT(IN)   :: CG(NR)
      REAL(8),   INTENT(IN)   :: FR(NR)
      REAL(8),   INTENT(IN)   :: FRP(NR)
      LOGICAL,  INTENT(IN)   :: THOM
      REAL(8)                 :: GAMMA
      REAL(8)                 :: RIGAMMA
      INTEGER(4)              :: IR
!
      IF(.NOT.THOM) THEN
        DO IR=1,4
          U(IR)=0.D0
          UP(IR)=0.D0
          UPP(IR)=CG(IR)
        ENDDO
        RETURN
      END IF
      IF(AEZ.GT.0) THEN
        IF(L .EQ. 0) THEN
          GAMMA=DSQRT(1.0D0-FSS*AEZ**2)
        ELSE
          GAMMA=(    L*DSQRT(   L**2 -FSS*AEZ**2)+ &
     &           (L+1)*DSQRT((L+1)**2-FSS*AEZ**2)   )/(2*L+1)
        END IF
      ELSE
        GAMMA=L+1
      END IF
!
      DO IR=1,4
        RIGAMMA=R(IR)**GAMMA
        U(IR)=RIGAMMA
        UP(IR)=DEX*GAMMA*RIGAMMA
        UPP(IR)=(DEX+FRP(IR))*UP(IR)+(CF(IR)+FR(IR))*U(IR)+CG(IR)
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE OUTWRD(DEX,NR,N1,N2,U,UP,UPP,CF,CG,FR,FRP)
!     ******************************************************************
!     **                                                              **
!     **  OUTWARD INTEGRATION OF THE SEMI-RELATIVISTIc                **
!     **  SCHROEDINGER EQUATION.                                      **
!     **  REQUIRES START VALUES OF U,UP,UPP FROM N1-3 TO N1           **
!     **                                                              **
!     ******************************************************************
!     **  ADAMS EXTRAPOLATION AND INTERPOLATION FORMULAS FOR          **
!     **  OUTWARD AND INWARD INTEGRATION, ABRAMOWITZ AND              **
!     **  STEGUN, P. 896                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NR
      REAL(8), INTENT(IN)    :: DEX
      INTEGER(4),INTENT(IN)  :: N1
      INTEGER(4),INTENT(IN)  :: N2
      REAL(8), INTENT(INOUT) :: U(NR)
      REAL(8), INTENT(INOUT) :: UP(NR)
      REAL(8), INTENT(INOUT) :: UPP(NR)
      REAL(8), INTENT(IN)    :: CF(NR)
      REAL(8), INTENT(IN)    :: CG(NR)
      REAL(8), INTENT(IN)    :: FR(NR)
      REAL(8), INTENT(IN)    :: FRP(NR)
      INTEGER(4)             :: I,IT
!     ******************************************************************
      DO I=N1,N2-1
        U(I+1)=U(I)+AEO(UP,I)
        UP(I+1)=UP(I)+AEO(UPP,I)
        DO IT=1,2
          UPP(I+1)= (DEX+FRP(I+1))*UP(I+1) &
     &            + (CF(I+1)+FR(I+1))*U(I+1) + CG(I+1)
          UP(I+1)=UP(I)+AIO(UPP,I)
          U(I+1)=U(I)+AIO(UP,I)
        ENDDO
      ENDDO
      CONTAINS
!       ...........................................FUNCTION AEO...........
        DOUBLE PRECISION FUNCTION AEO(Y,J)
!       **                                                            **
!       **  PREDICTOR                                                 **
!       **                                                            **
        REAL(8),   PARAMETER    :: C1=+55.D0/24.D0
        REAL(8),   PARAMETER    :: C2=-59.D0/24.D0
        REAL(8),   PARAMETER    :: C3=+37.D0/24.D0
        REAL(8),   PARAMETER    :: C4=- 9.D0/24.D0
        INTEGER(4),INTENT(IN)   :: J
        REAL(8)   ,INTENT(INOUT):: Y(J)
        AEO=C1*Y(J)+C2*Y(J-1)+C3*Y(J-2)+C4*Y(J-3)
        RETURN
        END FUNCTION AEO
!       ...........................................FUNCTION AIO...........
        DOUBLE PRECISION FUNCTION AIO(Y,J)
!       **                                                            **
!       **  CORRECTOR                                                 **
!       **                                                            **
        INTEGER(4),INTENT(IN)   :: J
        REAL(8),   PARAMETER    :: C1=+ 9.D0/24.D0
        REAL(8),   PARAMETER    :: C2=+19.D0/24.D0
        REAL(8),   PARAMETER    :: C3=- 5.D0/24.D0
        REAL(8),   PARAMETER    :: C4=+ 1.D0/24.D0
        REAL(8),   INTENT(INOUT):: Y(J+1)
        AIO=C1*Y(J+1)+C2*Y(J)+C3*Y(J-1)+C4*Y(J-2)
        RETURN
        END FUNCTION AIO
      END
!
!     ..................................................................
      SUBROUTINE INWRD(DEX,NR,N1,N2,U,UP,UPP,CF,CG,FR,FRP)
!     ******************************************************************
!     **                                                              **
!     **  INWARD INTEGRATION OF THE SEMI-RELATIVISTI!                 **
!     **  SCHROEDINGER EQUATION.                                      **
!     **  REQUIRES START VALUES OF U,UP,UPP FROM N1 TO N1+3.          **
!     **                                                              **
!     ******************************************************************
!     **  ADAMS EXTRAPOLATION AND INTERPOLATION FORMULAS FOR          **
!     **  OUTWARD AND INWARD INTEGRATION, ABRAMOWITZ AND              **
!     **  STEGUN, P. 896                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NR        ! 
      INTEGER(4), INTENT(IN) :: N1        ! INTEGRATION STARTS AT GRID POINT N1+1
      INTEGER(4), INTENT(IN) :: N2        ! N2 LAST CLACULATED GRIDPOINT FOR INTEGRATION
      REAL(8), INTENT(IN)    :: DEX
      REAL(8), INTENT(INOUT) :: U(NR)   ! R*PHI
      REAL(8), INTENT(INOUT) :: UP(NR)  ! DU/DX ; X=LN(R/R1)/DEX
      REAL(8), INTENT(INOUT) :: UPP(NR) ! D2U/DX2
      REAL(8), INTENT(IN)    :: CF(NR)
      REAL(8), INTENT(IN)    :: CG(NR)
      REAL(8), INTENT(IN)    :: FR(NR)
      REAL(8), INTENT(IN)    :: FRP(NR)
      INTEGER(4)             :: I,IT
!     ******************************************************************
      IF(NR.LE.N1+3) THEN
        STOP 'IN INWRD'
      ENDIF
      DO I=N1,N2+1,-1
        U(I-1)=U(I)+AEI(UP,I)
        UP(I-1)=UP(I)+AEI(UPP,I)
        DO IT=1,2
          UPP(I-1)=(DEX+FRP(I-1))*UP(I-1)+(CF(I-1)+FR(I-1))*U(I-1)+CG(I-1)
          UP(I-1)=UP(I)+AII(UPP,I)
          U(I-1)=U(I)+AII(UP,I)
        ENDDO
      ENDDO
      RETURN
      CONTAINS
!       ...........................................FUNCTION AII.........
        FUNCTION AII(Y,J) RESULT(RES)
        REAL(8),   PARAMETER    :: C1=- 9.D0/24.D0
        REAL(8),   PARAMETER    :: C2=-19.D0/24.D0
        REAL(8),   PARAMETER    :: C3=+ 5.D0/24.D0
        REAL(8),   PARAMETER    :: C4=- 1.D0/24.D0
        INTEGER(4),INTENT(IN)   :: J
        REAL(8)   ,INTENT(IN)   :: Y(J+2)
        REAL(8)                 :: RES  ! RESULT
        RES=C1*Y(J-1)+C2*Y(J)+C3*Y(J+1)+C4*Y(J+2)
        RETURN
        END FUNCTION AII
!       ...........................................FUNCTION AEI.........
        FUNCTION AEI(Y,J) RESULT(RES)
        REAL(8),   PARAMETER    :: C1=-55.D0/24.D0
        REAL(8),   PARAMETER    :: C2=+59.D0/24.D0
        REAL(8),   PARAMETER    :: C3=-37.D0/24.D0
        REAL(8),   PARAMETER    :: C4=+ 9.D0/24.D0
        INTEGER(4),INTENT(IN)   :: J
        REAL(8)   ,INTENT(IN)   :: Y(J+3)
        REAL(8)                 :: RES ! RESULT
        RES=C1*Y(J)+C2*Y(J+1)+C3*Y(J+2)+C4*Y(J+3)
        RETURN
        END FUNCTION AEI
      END
!
!     ..................................................................
      SUBROUTINE OUTBOUNDARY(THOM,R1,DEX,NR,R,IRBND,BOXRADIUS,U,UP,UPP,CF,DLG)
!     **                                                              **
!     **  FINDS BOUNDARY CONDITIONS FROM THE SEMICLASSICAL SOLUTION   **
!     **                                                              **
!     **   U(X)=EXP( +/- INT[DX:SQRT(CF(X))] ); R(X)=R1*EXP(DEX*(X-1))**
!     **                                                              **
!     **   IF DLG IS NOT PRESENT, DLG=INFTY, I.E. PHI=0 IS ASSUMED    **
!     **                                                              **
      IMPLICIT NONE
      LOGICAL   ,INTENT(IN)   :: THOM      ! SWITCH FOR HOMOGENEOUS SOLUTION
      INTEGER(4),INTENT(IN)   :: NR        ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN)   :: IRBND     ! SOLVE FROM IRBND TO NR
      REAL(8),   INTENT(IN)   :: R1        ! FIRST POINT ON RADIAL GRID
      REAL(8),   INTENT(IN)   :: DEX       ! 
      REAL(8),   INTENT(IN)   :: BOXRADIUS ! RADIUS FOR WHICH BOUNDARY CONDITIONS ARE SPECIFIED
      REAL(8),   INTENT(IN)   :: R(NR)     ! RADIAL GRID
      REAL(8),   INTENT(IN),OPTIONAL :: DLG  ! LOGARITHMIC DERIVATIVE AT BOXRADIUS
      REAL(8),   INTENT(INOUT):: U(NR)     ! R*PHI
      REAL(8),   INTENT(INOUT):: UP(NR)    !
      REAL(8),   INTENT(INOUT):: UPP(NR)   !
      REAL(8),   INTENT(IN)   :: CF(NR)    ! (L(L+1)+2*R**2*(POT(Y0)-E))/DEX**2
      COMPLEX(8)              :: PHASE(NR)
      COMPLEX(8)              :: DCSVAR1,DCSVAR2
      INTEGER(4)              :: IR
      INTEGER(4)              :: IR0,ir1,ir2 ! GRID POINT NEXT TO BOXRADIUS WITH IR0<NR
      REAL(8)                 :: XBOXRADIUS !
      REAL(8)                 :: SVAR,dx
      COMPLEX(8)              :: CFAC1,cfac2,csvar,dphasedx
      real(8)                 :: shift
      complex(8)              :: expphasep,expphasem
!     ******************************************************************
!
!     ==================================================================
!     ==  START WITH U==0 FOR INHOMOGENEOUS SOLUTION                  ==
!     ==================================================================
      IF((.NOT.THOM).AND.(.NOT.PRESENT(DLG))) THEN
        DO IR=IRBND,NR
          U(IR)=0.D0
          UP(IR)=0.D0
          UPP(IR)=0.D0
        ENDDO
      ENDIF
!
!     ==================================================================
!     == INTEGRATE PHASE(X0)=INT{DX,BOXRADIUS,X0) SQRT(CF(X))         ==
!     ==================================================================
      PHASE(1)=(0.D0,0.D0)      
      DO IR=2,NR
        DCSVAR1=MYSQRT(CF(IR-1))       ! SQRT(2*(V-E)) AT BOX RADIUS
        DCSVAR2=MYSQRT(CF(IR))
        PHASE(IR)=PHASE(IR-1)+0.5D0*(DCSVAR1+DCSVAR2)
      ENDDO
!
!     ==================================================================
!     == set phase at the boxradius to zero                           ==
!     ==================================================================
      XBOXRADIUS=1.D0+LOG(BOXRADIUS/R1)/DEX
      IR0=MIN(INT(XBOXRADIUS),NR-1) ! NEXT GRID POINT SMALLER THAN BOX RADIUS
      DX=MIN(1.D0,XBOXRADIUS-DBLE(IR0))
      DCSVAR1=MYSQRT(CF(IR0))       ! SQRT(2*(V-E)) AT BOX RADIUS
      DCSVAR2=MYSQRT(CF(IR0+1))
      CSVAR=PHASE(IR0)+DX*((1.d0-0.5D0*DX)*DCSVAR1+0.5D0*DX*DCSVAR2)
      IF(XBOXRADIUS.GT.DBLE(NR)) THEN
        DX=XBOXRADIUS-DBLE(NR)
        CSVAR=CSVAR+DX*MYSQRT(CF(NR))
      END IF
      DO IR=1,NR
        PHASE(IR)=PHASE(IR)-CSVAR
      ENDDO
!
!     ==================================================================
!     == SPECIFY BOUNDARY CONDITIONS                                  ==
!     == DCSVAR1=DPHASE1/DX=DPHASE2/DX, WHERE R(X)=R1*EXP(DEX*(X-1))  ==
!     == TAKEN AT R=BOXRADIUS                                         ==
!     ==================================================================
      IF(PRESENT(DLG)) THEN
!       ==  DEX*(DLG+1)=UP/U === UP=DU/DX ==============================
!       ==  DEX*(DLG+1)=UP/U ===========================================
        IF(XBOXRADIUS.GE.DBLE(NR-1)) THEN
          DPHASEDX=MYSQRT(CF(NR))   !=DPHASE1/DX
        ELSE
!         ATTENTION! OBTAINING DPHASE/DX AS (+/-)SQRT(CF)*PHASE IS INACCURATE
          IR0=INT(XBOXRADIUS)
          DPHASEDX=(PHASE(IR0+1)-PHASE(IR0))
        END IF
        CFAC2=(DPHASEDX-DEX*(DLG+1.D0))/(DPHASEDX+DEX*(DLG+1.D0))
      ELSE
        CFAC2=(-1.D0,0.D0)+(1.d-8,0.d0)
      END IF
      CFAC1=1.D0/(1.D0+CONJG(CFAC2)*CFAC2)
      CFAC2=CFAC2*CFAC1
!
!     ================================================================
!     == calculate prefactor to avoid overflow                      ==
!     ================================================================
      shift=0.d0
      DO IR=1,nr
        shift=max(shift,abs(real(phase(ir),kind=8)))
      ENDDO
      shift=shift-300.d0
!
!     ================================================================
!     == CALCULATE U(R)=EXP(PHASE1)-EXP(PHASE2)                     ==
!     ================================================================
      DO IR=2,NR-1
        EXPPHASEP=CFAC1*EXP(PHASE(IR)-SHIFT)
        EXPPHASEM=CFAC2*EXP(-PHASE(IR)-SHIFT)
        U(IR)=DBLE(EXPPHASEP+EXPPHASEM)
        DCSVAR1=0.5D0*(PHASE(IR+1)-PHASE(IR-1))
        UP(IR)=DBLE(DCSVAR1*EXPPHASEP-DCSVAR1*EXPPHASEM)
        DCSVAR1=DCSVAR1*DCSVAR1
        UPP(IR)=DBLE(DCSVAR1*EXPPHASEP+DCSVAR1*EXPPHASEM)
        DCSVAR1=PHASE(IR+1)-2.D0*PHASE(IR)+PHASE(IR-1)
        UPP(IR)=UPP(IR)+DBLE(DCSVAR1*EXPPHASEP-DCSVAR1*EXPPHASEM)
      ENDDO
!
!     ================================================================
!     == CORRECT LAST POINT IF BOX LARGER THAN GRID                   ==
!     ================================================================
      U(NR)=U(NR-1)+UP(NR-1)+0.5D0*UPP(NR-1)
      UP(NR)=UP(NR-1)+UPP(NR-1)
      UPP(NR)=UPP(NR-1)
      U(1)=U(2)-UP(2)+0.5D0*UPP(2)
      UP(1)=UP(2)-UPP(2)
      UPP(1)=UPP(2)
!
!     ================================================================
!     == RENORMALIZE                                               ==
!     ================================================================
      IR=INT(XBOXRADIUS)
      IR1=MIN(IRBND,Ir)
      IR2=MAX(IRBND,Ir+1)
      ir2=min(ir2,nr)
      SVAR=0.D0
      DO IR=IR1,IR2
        SVAR=MAX(SVAR,ABS(U(IR)))
      ENDDO
      IF(SVAR.EQ.0.D0) THEN
        DO IR=4,NR,10
          PRINT*,'IR ',IR,U(IR),EXP(PHASE(IR)-SHIFT)
        ENDDO
        CALL ERROR$MSG('NORMALIZATION IS ZERO')
        CALL ERROR$I4VAL('IR1',IR1)
        CALL ERROR$I4VAL('IR2',IR2)
        CALL ERROR$I4VAL('NR',NR)
        CALL ERROR$I4VAL('IRBND',IRBND)
        CALL ERROR$I4VAL('IRBOXRADIUS',INT(XBOXRADIUS))
        CALL ERROR$L4VAL('THOM',THOM)
        CALL ERROR$L4VAL('PRESENT(DLG)',PRESENT(DLG))
        CALL ERROR$STOP('OUTBOUNDARY')
      END IF
      SVAR=1.D0/SVAR
      DO IR=1,NR
        U(IR)=SVAR*U(IR)
        UP(IR)=SVAR*UP(IR)
        UPP(IR)=SVAR*UPP(IR)
      ENDDO
      if(.not.present(dlg)) then
        IR1=INT(XBOXRADIUS)+1
        ir1=min(ir1,nr)
        do ir=ir1,nr
          u(ir)=0.d0
          up(ir)=0.d0
          upp(ir)=0.d0
        enddo
      end if
!
!     ================================================================
!     == print for test                                             ==
!     ================================================================
      if(present(dlg)) then
         DO IR=1,NR
!          PRINT*,R(IR),U(IR)/R(IR),UP(IR),UPP(IR)
        ENDDO
!        stop
      END IF
      RETURN
      CONTAINS
        FUNCTION MYSQRT(X) RESULT(Y)
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: X
        COMPLEX(8)          :: Y
        REAL(8)             :: SVAR
        SVAR=DSQRT(ABS(X))
        IF(X.GE.0) THEN
          Y=cmplx(SVAR,0.D0)
        ELSE IF(X.LT.0) THEN
          Y=cmplx(0.D0,SVAR)
        END IF
        RETURN
        END FUNCTION MYSQRT
      END
!
!     ..................................................................
      SUBROUTINE oldOUTBOUNDARY(THOM,R1,DEX,NR,R,IRBND,BOXRADIUS,U,UP,UPP,CF,DLG)
!     **                                                              **
!     **  FINDS BOUNDARY CONDITIONS FROM THE SEMICLASSICAL SOLUTION   **
!     **                                                              **
!     **   U(X)=EXP( +/- INT[DX:SQRT(CF(X))] ); R(X)=R1*EXP(DEX*(X-1))**
!     **                                                              **
!     **   IF DLG IS NOT PRESENT, DLG=INFTY, I.E. PHI=0 IS ASSUMED    **
!     **                                                              **
      IMPLICIT NONE
      LOGICAL   ,INTENT(IN)   :: THOM      ! SWITCH FOR HOMOGENEOUS SOLUTION
      INTEGER(4),INTENT(IN)   :: NR        ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN)   :: IRBND     ! SOLVE FROM IRBND TO NR
      REAL(8),   INTENT(IN)   :: R1        ! FIRST POINT ON RADIAL GRID
      REAL(8),   INTENT(IN)   :: DEX       ! 
      REAL(8),   INTENT(IN)   :: BOXRADIUS ! RADIUS FOR WHICH BOUNDARY CONDITIONS ARE SPECIFIED
      REAL(8),   INTENT(IN)   :: R(NR)     ! RADIAL GRID
      REAL(8),   INTENT(IN),OPTIONAL :: DLG  ! LOGARITHMIC DERIVATIVE AT BOXRADIUS
      REAL(8),   INTENT(INOUT):: U(NR)     ! R*PHI
      REAL(8),   INTENT(INOUT):: UP(NR)    !
      REAL(8),   INTENT(INOUT):: UPP(NR)   !
      REAL(8),   INTENT(IN)   :: CF(NR)    ! (L(L+1)+2*R**2*(POT(Y0)-E))/DEX**2
      COMPLEX(8)              :: PHASE(NR)
      COMPLEX(8)              :: DCSVAR1,DCSVAR2
      INTEGER(4)              :: IR
      INTEGER(4)              :: IR0,ir1,ir2 ! GRID POINT NEXT TO BOXRADIUS WITH IR0<NR
      REAL(8)                 :: XBOXRADIUS !
      REAL(8)                 :: SVAR,dx
      COMPLEX(8)              :: CFAC1,cfac2,csvar,dphasedx
      real(8)                 :: phasmx,phasmn
!     ******************************************************************
      call error$msg('routine marked for deletion')
      call error$stop('oldoutboundary')
!
!     ==================================================================
!     ==  START WITH U==0 FOR INHOMOGENEOUS SOLUTION                  ==
!     ==================================================================
      IF((.NOT.THOM).AND.(.NOT.PRESENT(DLG))) THEN
        DO IR=IRBND,NR
          U(IR)=0.D0
          UP(IR)=0.D0
          UPP(IR)=0.D0
        ENDDO
!       RETURN
      ENDIF
!
!     ==================================================================
!     == INTEGRATE PHASE(X0)=INT{DX,BOXRADIUS,X0) SQRT(CF(X))         ==
!     ==================================================================
      PHASE(1)=(0.D0,0.D0)      
      DO IR=2,NR
        DCSVAR1=MYSQRT(CF(IR-1))       ! SQRT(2*(V-E)) AT BOX RADIUS
        DCSVAR2=MYSQRT(CF(IR))
        PHASE(IR)=PHASE(IR-1)+0.5D0*(DCSVAR1+DCSVAR2)
      ENDDO
!
!     ==================================================================
!     == set phase at the boxradius to zero                           ==
!     ==================================================================
      XBOXRADIUS=1.D0+LOG(BOXRADIUS/R1)/DEX
      IR0=MIN(INT(XBOXRADIUS),NR-1) ! NEXT GRID POINT SMALLER THAN BOX RADIUS
      DX=MIN(1.D0,XBOXRADIUS-DBLE(IR0))
      DCSVAR1=MYSQRT(CF(IR0))       ! SQRT(2*(V-E)) AT BOX RADIUS
      DCSVAR2=MYSQRT(CF(IR0+1))
      CSVAR=PHASE(IR0)+DX*((1.d0-0.5D0*DX)*DCSVAR1+0.5D0*DX*DCSVAR2)
      IF(XBOXRADIUS.GT.DBLE(NR)) THEN
        DX=XBOXRADIUS-DBLE(NR)
        CSVAR=CSVAR+DX*MYSQRT(CF(NR))
      END IF
      DO IR=1,NR
        PHASE(IR)=PHASE(IR)-CSVAR
      ENDDO
!
!     ==================================================================
!     == SPECIFY BOUNDARY CONDITIONS                                  ==
!     == DCSVAR1=DPHASE1/DX=DPHASE2/DX, WHERE R(X)=R1*EXP(DEX*(X-1))  ==
!     == TAKEN AT R=BOXRADIUS                                         ==
!     ==================================================================
      IF(PRESENT(DLG)) THEN
!       ==  DEX*(DLG+1)=UP/U === UP=DU/DX ==============================
!       ==  DEX*(DLG+1)=UP/U ===========================================
        IF(XBOXRADIUS.GE.DBLE(NR-1)) THEN
          DPHASEDX=MYSQRT(CF(NR))   !=DPHASE1/DX
        ELSE
!         ATTENTION! OBTAINING DPHASE/DX AS (+/-)SQRT(CF)*PHASE IS INACCURATE
          IR0=INT(XBOXRADIUS)
          DPHASEDX=(PHASE(IR0+1)-PHASE(IR0))
        END IF
        CFAC2=(DPHASEDX-DEX*(DLG+1.D0))/(DPHASEDX+DEX*(DLG+1.D0))
      ELSE
        CFAC2=(-1.D0,0.D0)
      END IF
      CFAC1=1.D0/(1.D0+CONJG(CFAC2)*CFAC2)
      CFAC2=CFAC2*CFAC1
!
!     ================================================================
!     == CALCULATE U(R)=EXP(PHASE1)-EXP(PHASE2)                     ==
!     ================================================================
      DO IR=2,NR-1
        U(IR)=DBLE(CFAC1*EXP(PHASE(IR))+CFAC2*EXP(-PHASE(IR)))
        DCSVAR1=0.5D0*(PHASE(IR+1)-PHASE(IR-1))
        UP(IR)=DBLE(DCSVAR1*CFAC1*EXP(PHASE(IR))-CFAC2*DCSVAR1*EXP(-PHASE(IR)))
        DCSVAR1=DCSVAR1*DCSVAR1
        UPP(IR)=DBLE(DCSVAR1*CFAC1*EXP(PHASE(IR))+CFAC2*DCSVAR1*EXP(-PHASE(IR)))
        DCSVAR1=PHASE(IR+1)-2.D0*PHASE(IR)+PHASE(IR-1)
        UPP(IR)=UPP(IR)+DBLE(DCSVAR1*CFAC1*EXP(PHASE(IR))-CFAC2*DCSVAR1*EXP(-PHASE(IR)))
      ENDDO
!
!     ================================================================
!     == CORRECT LAST POINT IF BOX LARGER THAN GRID                   ==
!     ================================================================
      U(NR)=U(NR-1)+UP(NR-1)+0.5D0*UPP(NR-1)
      UP(NR)=UP(NR-1)+UPP(NR-1)
      UPP(NR)=UPP(NR-1)
      U(1)=U(2)-UP(2)+0.5D0*UPP(2)
      UP(1)=UP(2)-UPP(2)
      UPP(1)=UPP(2)
!
!     ================================================================
!     == RENORMALIZE                                               ==
!     ================================================================
      IR=INT(XBOXRADIUS)
      IR1=MIN(IRBND,Ir)
      IR2=MAX(IRBND,Ir+1)
      ir2=min(ir2,nr)
      SVAR=0.D0
      DO IR=IR1,IR2
        SVAR=MAX(SVAR,ABS(U(IR)))
      ENDDO
      SVAR=1.D0/SVAR
      DO IR=1,NR
        U(IR)=SVAR*U(IR)
        UP(IR)=SVAR*UP(IR)
        UPP(IR)=SVAR*UPP(IR)
      ENDDO
      if(.not.present(dlg)) then
        IR1=INT(XBOXRADIUS)+1
        ir1=min(ir1,nr)
        do ir=ir1,nr
          u(ir)=0.d0
          up(ir)=0.d0
          upp(ir)=0.d0
        enddo
      end if
!
!     ================================================================
!     == print for test                                             ==
!     ================================================================
      if(present(dlg)) then
         DO IR=1,NR
!          PRINT*,R(IR),U(IR)/R(IR),UP(IR),UPP(IR)
        ENDDO
!        stop
      END IF
      RETURN
      CONTAINS
        FUNCTION MYSQRT(X) RESULT(Y)
        IMPLICIT NONE
        REAL(8), INTENT(IN) :: X
        COMPLEX(8)          :: Y
        REAL(8)             :: SVAR
        SVAR=DSQRT(ABS(X))
        IF(X.GE.0) THEN
          Y=cmplx(SVAR,0.D0)
        ELSE IF(X.LT.0) THEN
          Y=cmplx(0.D0,SVAR)
        END IF
        RETURN
        END FUNCTION MYSQRT
      END
