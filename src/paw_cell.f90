!........1.........2.........3.........4.........5.........6.........7.........8
!*******************************************************************************
!**  OBJECT CELL : PARRINELLO RAHMAN                                          **
!**  USAGE                                                                    **
!**    ==  INITIALIZE: ================================                       **
!**    SET MANDATORY DT,TREF,TMASS                                            **
!**    SET OPTINALLY: STRESS,PRESSURE,FRIC                                    **
!**    GET('T0',9,T0)                                                         **
!**    ==  LOOP: ======================================                       **
!**    --  CALCULATE STRESS                                                   **
!**    SETR8A('STRESS_I',9,STRESS)                                            **
!**    PROPAGATE                                                              **
!**    GET('FRICMAT',9,U)                                                     **
!**    -- ADD OTHERFRICTION TERMS TO U                                        **
!**    -- U=U*DT/2                                                            **
!**    -- PROPAGATE POSITIONS                                                 **
!**       R(+)=(1+U)**(-1)[2R(0)-(1-U)R(-)+F*DT^2/M                           **
!**    GET('MAPTOCELL',3,U)                                                   **
!**    -- RP=U*RP ; R0=U*R0 ; RM=U*RM                                         **
!**    -- RBAS=T0                                                             **
!**    GET('EKIN',EKIN)                                                       **
!**    GET('EPOT',EPOT)                                                       **
!**    SWITCH                                                                 **
!**    GET('T0',9,T0)                                                         **
!**                                                                           **
!*******************************************************************************
MODULE CELL_MODULE
IMPLICIT NONE
logical(4) :: tpARRINELLORAHMAN=.TRUE.
! CONSTRAINTTYPE CAN BE 'NONE','ISOTROPIC','NOSHEAR','FREE'
CHARACTER(32) :: CONSTRAINTTYPE='FREE'
LOGICAL(4) :: TINIT=.FALSE.
LOGICAL(4) :: TON  =.FALSE.      ! USED TO REQUEST INTERNAL STRESS
LOGICAL(4) :: TMOVE=.FALSE.      ! PROPAGATE UNIT CELL
REAL(8)    :: TREF(3,3)          ! REFERENCE CELL CONSTANTS
REAL(8)    :: STRESS(3,3)=0.D0   ! EXTERNAL STRESS TENSOR
REAL(8)    :: PRESSURE=0.D0      ! EXTERNAL PRESSURE
REAL(8)    :: DELTAT=0.D0
REAL(8)    :: FRIC=0.D0
REAL(8)    :: STRESS_I(3,3)=0.D0 ! INTERNAL STRESS TENSOR
REAL(8)    :: KINSTRESS(3,3)=0.D0 ! INTERNAL STRESS TENSOR
REAL(8)    :: TP(3,3)=0.D0
REAL(8)    :: T0(3,3)=0.D0        ! ACTUAL CELL 
REAL(8)    :: TM(3,3)=0.D0
REAL(8)    :: TMM(3,3)=0.D0
REAL(8)    :: TMASS              ! MASS FOR THE UNIT CELL DYNAMICS
LOGICAL(4) :: TPROPAGATED=.FALSE.
REAL(8)    :: EPOT=0.D0          ! POTENTIAL ENERGY
REAL(8)    :: EKIN=0.D0          ! KINETIC ENERGY
REAL(8)    :: VREF=0.D0          ! VOLUME OF THE REFERENCE UNIT CELL
REAL(8)    :: V0=0.D0            ! VOLUME OF THE ACTUAL UNIT CELL
LOGICAL(4) :: TSTOP=.FALSE.      ! RESET VELOCITIES TO ZERO
!=====================================
REAL(8)    :: SIGMA(3,3)      ! 
REAL(8)    :: TREFINV(3,3)    ! 
CONTAINS
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL_INITIALIZE()
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) :: I
      REAL(8)    :: SVAR
!     **************************************************************************
      IF(TINIT) RETURN
      TINIT=.TRUE.
!
!     ==================================================================
!     ==  CALCULATE TREFINV                                           ==
!     ==================================================================
      IF(SUM(TREF**2).EQ.0.D0) THEN
        CALL ERROR$MSG('REFERENCE CELL NOT SET')
        CALL ERROR$STOP('CELL_INITIALIZE')
      END IF
      VREF=TREF(1,1) * (TREF(2,2)*TREF(3,3)-TREF(3,2)*TREF(2,3)) &
     &    +TREF(2,1) * (TREF(3,2)*TREF(1,3)-TREF(1,2)*TREF(3,3)) &
     &    +TREF(3,1) * (TREF(1,2)*TREF(2,3)-TREF(2,2)*TREF(1,3))
      CALL LIB$INVERTR8(3,TREF,TREFINV)
!
!     ==================================================================
!     ==  SET DEFAULT CELL IF NOT DEFINED OTHERWISE                   ==
!     ==================================================================
      IF(SUM(T0**2).EQ.0.D0) THEN   ! T0 IS NOT SET
        T0=TREF
        TM=T0
        TMM=TM
      ELSE                          ! T0 IS SET
        IF(SUM(TM**2).EQ.0.D0) THEN ! T0 IS SET BUT TM IS NOT     
          TM=T0
          TMM=TM
        ELSE                        ! T0 AND TM ARE SET
          IF(SUM(TMM**2).EQ.0) THEN ! T0 AND TM ARE SET, BUT TMM IS NOT
            TMM=2.D0*TM-T0
          END IF
        END IF
      END IF
      TP=2.D0*T0-TM
      IF(TSTOP) THEN
        TMM=T0
        TM=T0
        TP=T0
      END IF
!
!     ===================================================================
!     ==  CALCULATE SIGMA MATRIX EQ 2.24                               ==
!     ===================================================================
      SIGMA=STRESS
      SVAR=0.D0
      DO I=1,3
        SVAR=SVAR+SIGMA(I,I)
      ENDDO
      SVAR=SVAR/3.D0
      DO I=1,3
        SIGMA(I,I)=SIGMA(I,I)-SVAR
      END DO
!
!     ==========================================================================
!     ==  OTHER CHECKS                                                        ==
!     ==========================================================================
      IF(DELTAT.EQ.0.D0) THEN
        CALL ERROR$MSG('TIME STEP HAS NOT BEEN SET')
        CALL ERROR$STOP('CELL_INITIALIZE')
      END IF
      IF(TMASS.EQ.0.D0) THEN
        CALL ERROR$MSG('CELL-MASS HAS NOT BEEN SET')
        CALL ERROR$STOP('CELL_INITIALIZE')
      END IF
!
      RETURN
      END SUBROUTINE CELL_INITIALIZE
END MODULE CELL_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$REPORT(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CELL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
!     **************************************************************************
      IF(.NOT.TON) RETURN
      CALL REPORT$TITLE(NFIL,'UNIT CELL')
      CALL REPORT$L4VAL(NFIL,'DYNAMICAL UNIT CELL',TMOVE) 
      CALL REPORT$R8VAL(NFIL,'MASS',TMASS,'A.U.') 
      CALL REPORT$R8VAL(NFIL,'FRICTION',FRIC,'DT/2') 
      CALL REPORT$R8VAL(NFIL,'EXTERNAL PRESSURE',PRESSURE,'A.U.') 
      CALL REPORT$CHVAL(NFIL,'CONSTRAINTYPE',CONSTRAINTTYPE) 
      CALL REPORT$R8VAL(NFIL,'REFERENCE VOLUME',VREF,'A.U.') 
      WRITE(NFIL,FMT='("REFERENCE UNIT CELL",T35," STRESS ")')
      WRITE(NFIL,FMT='(3F10.5,T35,3F10.5)')TREF(1,:),STRESS(1,:)
      WRITE(NFIL,FMT='(3F10.5,T35,3F10.5)')TREF(2,:),STRESS(2,:)
      WRITE(NFIL,FMT='(3F10.5,T35,3F10.5)')TREF(3,:),STRESS(3,:)
      END SUBROUTINE CELL$REPORT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$CONVERT(DT,B,TREF,PERIOD,MASS,FRICTION)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: DT ! TIME STEP
      REAL(8),INTENT(IN) :: B  ! BULK MODULUS
      REAL(8),INTENT(IN) :: TREF(3,3)  ! REFERENCE UNIT CELL
      REAL(8),INTENT(IN) :: PERIOD ! PERIOD OF CELL OSCILLATION
      REAL(8),INTENT(OUT):: MASS  ! MASS OF THE UNIT CELL
      REAL(8),INTENT(OUT):: FRICTION ! FRICTION PARAMETER FOR CRITICAL DAMPING
      REAL(8)            :: PI
      REAL(8)            :: OMEGA
      REAL(8)            :: VOL
      REAL(8)            :: AMAT(3,3)
!     *****************************************************************
      CALL ERROR$MSG('NOT FULLY IMPLEMENTED: DO NOT USE')
      CALL ERROR$STOP('CELL$CONVERT')
      PI=4.D0*ATAN(1.D0)
      CALL GBASS(TREF,AMAT,VOL)
      OMEGA=2.D0*PI/PERIOD
      MASS=9.D0*B/(OMEGA*VOL)**2
      FRICTION=OMEGA*DT
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$SETCH(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CELL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
!     **************************************************************************
!
!     =================================================================
!     ==  CALCULATE PRESSURE                                         ==
!     =================================================================
      IF(ID.EQ.'CONSTRAINTTYPE') THEN
        CONSTRAINTTYPE=VAL
!
!     =================================================================
!     ==  DONE                                                       ==
!     =================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('CELL$SETCH')
      END IF
      RETURN
      END SUBROUTINE CELL$SETCH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$SETL4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CELL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
!
!     =================================================================
!     ==  CALCULATE PRESSURE                                         ==
!     =================================================================
      IF(ID.EQ.'ON') THEN
        TON=VAL
!
!     =================================================================
!     ==  PROPAGATE                                                  ==
!     =================================================================
      ELSE IF(ID.EQ.'MOVE') THEN
        TMOVE=VAL
        IF(TMOVE) THEN
          TON=.TRUE.
        END IF
!
!     =================================================================
!     ==  SET VELOCITY TO ZERO IN THE NEXT TIME STEP                 ==
!     =================================================================
      ELSE IF(ID.EQ.'STOP') THEN
        TSTOP=VAL
!
!     =================================================================
!     ==  DONE                                                       ==
!     =================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('CELL$SETL4')
      END IF
      RETURN
      END SUBROUTINE CELL$SETL4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$GETL4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CELL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
!
!     =================================================================
!     ==  CALCULATE PRESSURE AND STRESS                              ==
!     =================================================================
      IF(ID.EQ.'ON') THEN
        VAL=TON
!
!     =================================================================
!     ==  PROPAGATE                                                  ==
!     =================================================================
      ELSE IF(ID.EQ.'MOVE') THEN
        VAL=TMOVE
!
!     =================================================================
!     ==  DONE                                                       ==
!     =================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('CELL$GETL4')
      END IF
      RETURN
      END SUBROUTINE CELL$GETL4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$SETR8(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CELL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
!     ******************************************************************
!
!     ==================================================================
!     ==  PRESSURE                                                    ==
!     ==================================================================
      IF(ID.EQ.'P') THEN
        PRESSURE=VAL
!
!     ==================================================================
!     ==  TIME STEP                                                   ==
!     ==================================================================
      ELSE IF(ID.EQ.'DT') THEN
        DELTAT=VAL
!
!     ==================================================================
!     ==  FRICTION                                                    ==
!     ==================================================================
      ELSE IF(ID.EQ.'FRICTION') THEN
        FRIC=VAL
!
!     ==================================================================
!     ==  MASS                                                        ==
!     ==================================================================
      ELSE IF(ID.EQ.'MASS') THEN
        TMASS=VAL
!
!     ==================================================================
!     ==  DONE                                                        ==
!     ==================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('CELL$SETR8')
      END IF
      RETURN
      END SUBROUTINE CELL$SETR8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$GETR8(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CELL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
!     **************************************************************************
!
!     =================================================================
!     ==  POTENTIAL ENERGY                                           ==
!     =================================================================
      IF(ID.EQ.'EPOT') THEN
        IF(TMOVE.AND.(.NOT.TPROPAGATED)) THEN
          CALL ERROR$MSG('DATA AVALAIBLE ONLY AFTER PROPAGATION')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('CELL$GETR8')
        END IF
        VAL=EPOT
!
!     =================================================================
!     ==  KINETIC ENERGY                                             ==
!     =================================================================
      ELSE IF(ID.EQ.'EKIN') THEN
        IF(TMOVE.AND.(.NOT.TPROPAGATED)) THEN
          CALL ERROR$MSG('DATA AVALAIBLE ONLY AFTER PROPAGATION')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('CELL$GETR8')
        END IF
        VAL=EKIN
!
!     =================================================================
!     ==  MASS FOR CELL DYNAMICS                                     ==
!     =================================================================
      ELSE IF(ID.EQ.'MASS') THEN
        VAL=TMASS
!
!     =================================================================
!     ==  DONE                                                       ==
!     =================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('CELL$GETR8')
      END IF
      RETURN
      END SUBROUTINE CELL$GETR8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$SETR8A(ID,LEN,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CELL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(IN) :: VAL(LEN)
!     ******************************************************************
!
!     ==================================================================
!     ==  REFERENCE CELL VECTORS                                   ==
!     ==================================================================
      IF(ID.EQ.'TREF') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE MISMATCH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('CELL$SETR8A')
        END IF
        TREF=RESHAPE(VAL,(/3,3/))

!
!     ==================================================================
!     ==  EXTERNAL STRESS TENSOR                                      ==
!     ==================================================================
      ELSE IF(ID.EQ.'STRESS') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE MISMATCH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('CELL$SETR8A')
        END IF
        STRESS=RESHAPE(VAL,(/3,3/))
!
!     ==================================================================
!     ==  INTERNAL STRESS TENSOR                                      ==
!     ==================================================================
      ELSE IF(ID.EQ.'STRESS_I') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE MISMATCH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('CELL$SETR8A')
        END IF
        STRESS_I=RESHAPE(VAL,(/3,3/))
!
!     ==================================================================
!     ==  INTERNAL STRESS TENSOR                                      ==
!     ==================================================================
      ELSE IF(ID.EQ.'KINSTRESS') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE MISMATCH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('CELL$SETR8A')
        END IF
        KINSTRESS=RESHAPE(VAL,(/3,3/))
!
!     ==================================================================
!     ==  CELL VECTORS                                             ==
!     ==================================================================
      ELSE IF(ID.EQ.'T0'.OR.ID.EQ.'T(0)') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE MISMATCH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('CELL$SETR8A')
        END IF
        T0=RESHAPE(VAL,(/3,3/))
!
!     ==================================================================
!     ==  CELL VECTORS OF PREVIOUS TIME STEP                       ==
!     ==================================================================
      ELSE IF(ID.EQ.'TM'.OR.ID.EQ.'T(-)') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE MISMATCH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('CELL$SETR8A')
        END IF
        TM=RESHAPE(VAL,(/3,3/))
!
!     ==================================================================
!     ==  DONE                                                        ==
!     ==================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('CELL$SETR8A')
      END IF
      RETURN
      END SUBROUTINE CELL$SETR8A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$GETR8A(ID,LEN,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CELL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(OUT):: VAL(LEN)
      REAL(8)                 :: AMAT(3,3)
!     **************************************************************************
!
!     ==========================================================================
!     ==  CELL VECTORS                                                        ==
!     ==========================================================================
      IF(ID.EQ.'T0'.OR.ID.EQ.'T(0)') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE MISMATCH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('CELL$GETR8A')
        END IF
        VAL=RESHAPE(T0,(/9/))
!
!     ==========================================================================
!     ==  CELL VECTORS FOR THE NEXT TIME STEP                                 ==
!     ==========================================================================
      ELSE IF(ID.EQ.'TP'.OR.ID.EQ.'T(+)') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE MISMATCH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('CELL$GETR8A')
        END IF
        IF(.NOT.TON) THEN
          CALL ERROR$MSG('CELL DYNAMICS IS NOT SWITCHED ON')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('CELL$GETR8A')
        END IF 
        CALL CELL_INITIALIZE()
        IF(TINIT) THEN
          VAL=RESHAPE(TP,(/9/))
        ELSE
          VAL=RESHAPE(TREF,(/9/))
        END IF
!
!     ==========================================================================
!     ==  REFERENCE CELL VECTORS                                              ==
!     ==========================================================================
      ELSE IF(ID.EQ.'TREF') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE MISMATCH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('CELL$GETR8A')
        END IF
        VAL=RESHAPE(TREF,(/9/))
!
!     ==========================================================================
!     ==  STRESS                                                              ==
!     ==========================================================================
      ELSE IF(ID.EQ.'STRESS_I') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE MISMATCH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('CELL$GETR8A')
        END IF
        VAL=RESHAPE(STRESS_I,(/9/))
!
!     ==========================================================================
!     ==  FRICTION MATRIX FOR ATOM DYNAMICS                                   ==
!     ==========================================================================
      ELSE IF(ID.EQ.'FRICMAT') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE MISMATCH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('CELL$GETR8A')
        END IF
        IF(TON) THEN
          CALL CELL_INITIALIZE()
          CALL LIB$INVERTR8(3,T0,AMAT)
          AMAT=MATMUL(TP-TM,AMAT)/(2.D0*DELTAT)
          AMAT=AMAT+TRANSPOSE(AMAT)
          AMAT=AMAT*0.5D0*DELTAT
          VAL=RESHAPE(AMAT,(/9/))
        ELSE
          VAL=0.D0
        ENDIF
!
!     ==========================================================================
!     ==  MATRIX TO MAP POSITIONS INTO NEW UNIT CELL                          ==
!     ==========================================================================
      ELSE IF(ID.EQ.'MAPTOCELL') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE MISMATCH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('CELL$GETR8A')
        END IF
        IF(TON) THEN
          CALL CELL_INITIALIZE()
          IF(.NOT.TPROPAGATED) THEN
            CALL ERROR$MSG('MAPTOCELL IS NOT AVAILABLE')
            CALL ERROR$MSG('CELL HAS NOT YET BEEN PROPAGATED')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('CELL$GETR8A')
          END IF
          CALL LIB$INVERTR8(3,T0,AMAT)
          AMAT=MATMUL(TP,AMAT)
          VAL=RESHAPE(AMAT,(/9/))
        ELSE
          VAL=(/1.D0,0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,1.D0/)
        END IF
!
!     ==========================================================================
!     ==  DONE                                                                ==
!     ==========================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('CELL$GETR8A')
      END IF
      RETURN
      END SUBROUTINE CELL$GETR8A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$STOP()
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CELL_MODULE
      IMPLICIT NONE
!     **************************************************************************
      TSTOP=.TRUE.
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$PROPAGATE()
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CELL_MODULE
      IMPLICIT NONE
      REAL(8)    :: AMAT(3,3)
      REAL(8)    :: BMAT(3,3)
      REAL(8)    :: CMAT(3,3)
      REAL(8)    :: T0INV(3,3)
      REAL(8)    :: ALPHAP(3,3),ALPHA0(3,3),ALPHAM(3,3)
      REAL(8)    :: ALPHADOT(3,3)
      REAL(8)    :: XPMAT(3,3),XMMAT(3,3)
      REAL(8)    :: ONE(3,3)
      REAL(8)    :: stress_ext(3,3)
REAL(8)    :: mat33(3,3)
      REAL(8)    :: SVAR,SVAR1,SVAR2,SVAR3
      INTEGER(4) :: I,j,ITER,IC1,IC2
      logical(4),parameter :: donothing=.false.
      integer(4) :: nconstraint !#(constraints)
      real(8)   ,ALLOCATABLE :: constraintproject(:,:,:)    ! (3,3,nc)
      real(8)   ,ALLOCATABLE :: A(:,:),X(:),B(:)
      REAL(8)                :: smat1(3,3),smat2(3,3),smat3(3,3)
!     **************************************************************************
      IF(.NOT.TON) RETURN
      IF(.NOT.TMOVE) RETURN
      CALL CELL_INITIALIZE()
if(donothing) then
  tmm=t0
  tm=t0
  tp=t0
  ekin=0.d0
  epot=0.d0
  TPROPAGATED=.TRUE.
  return
end if
      IF(TSTOP) THEN
        TMm=T0
        TM=T0
        TSTOP=.FALSE.
      END IF
      ONE=0.D0
      DO I=1,3
        ONE(I,I)=1.D0
      ENDDO
!
!     ==========================================================================
!     == ensure that stress-tensor is symmetric                               ==
!     ==========================================================================
!     == the stresses from potential_hartree and from paw_pairpotential are
!     == not exactly symmetric. The antisymmetric part resuls in a rotation
!     == therefore we remove it here
      STRESS_I=0.5D0*(STRESS_I+TRANSPOSE(STRESS_I))
      KINSTRESS=0.5D0*(KINSTRESS+TRANSPOSE(KINSTRESS))
!
!     ==========================================================================
!     == CONSTRAINTS                                                          ==
!     ==========================================================================
      IF(CONSTRAINTTYPE.EQ.'ISOTROPIC') THEN
        nconstraint=5
        allocate(constraintproject(3,3,nconstraint))
        constraintproject(:,:,:)=0.d0
        constraintproject(1,2,1)=1.d0
        constraintproject(2,1,1)=1.d0
        constraintproject(1,3,2)=1.d0
        constraintproject(3,1,2)=1.d0
        constraintproject(2,3,3)=1.d0
        constraintproject(3,2,3)=1.d0
        constraintproject(1,1,4)=1.d0
        constraintproject(2,2,4)=-1.d0
        constraintproject(2,2,5)=1.d0
        constraintproject(3,3,5)=-1.d0
!
!!$        SVAR1=(STRESS_I(1,1)+STRESS_i(2,2)+STRESS_i(3,3))/3.D0
!!$        SVAR2=(KINSTRESS(1,1)+KINSTRESS(2,2)+KINSTRESS(3,3))/3.D0
!!$        STRESS_i(:,:)=0.D0
!!$        KINSTRESS(:,:)=0.D0
!!$        DO I=1,3
!!$          STRESS_I(I,I)=SVAR1
!!$          KINSTRESS(I,I)=SVAR2
!!$        ENDDO
      ELSE IF(CONSTRAINTTYPE.EQ.'NOSHEAR') THEN
        nconstraint=3
        allocate(constraintproject(3,3,nconstraint))
        constraintproject(:,:,:)=0.d0
        constraintproject(1,2,1)=1.d0
        constraintproject(2,1,1)=1.d0
        constraintproject(1,3,2)=1.d0
        constraintproject(3,1,2)=1.d0
        constraintproject(2,3,3)=1.d0
        constraintproject(3,2,3)=1.d0
!!$        DO I=1,3
!!$          DO J=1,3
!!$            IF(I.EQ.J) CYCLE
!!$            STRESS_I(I,J)=0.D0
!!$            KINSTRESS(I,J)=0.D0
!!$          ENDDO
!!$        ENDDO
      ELSE IF(CONSTRAINTTYPE.EQ.'NOSTRESS') THEN
        nconstraint=0
        STRESS_I(:,:)=0.D0
        KINSTRESS(:,:)=0.D0
      ELSE IF(CONSTRAINTTYPE.EQ.'FREE') THEN
        nconstraint=0
      ELSE
        CALL ERROR$MSG('CONSTRAINTTYPE NOT RECOGNIZED')
        CALL ERROR$CHVAL('CONSTRAINTTYPE',CONSTRAINTTYPE)
        CALL ERROR$STOP('CELL$PROPAGATE')
      END IF
!
      V0=T0(1,1)*(T0(2,2)*T0(3,3)-T0(2,3)*T0(3,2)) &
     &  +T0(2,1)*(T0(3,2)*T0(1,3)-T0(3,3)*T0(1,2)) &
     &  +T0(3,1)*(T0(1,2)*T0(2,3)-T0(1,3)*T0(2,2)) 
      IF(ABS(V0).GT.1.D+10*ABS(VREF)) THEN
        CALL ERROR$MSG('CELL DYNAMICS UNSTABLE')
        CALL ERROR$R8VAL('CELL',T0)
        CALL ERROR$R8VAL('VOLUME ',V0)
        CALL ERROR$R8VAL('REFERENCE VOLUME ',VREF)
        CALL ERROR$STOP('CELL_PROPAGATE')
      END IF

      IF(TPARRINELLORAHMAN) THEN
!       == ENERGY AND STRESS OF THE VOLUME RESERVOIR. THE ENTHALPY IS  H=E+PV
        EPOT=PRESSURE*V0
        STRESS_EXT(:,:)=-PRESSURE*V0*ONE
!       == PROPAGATE LATTICE VECTORS ================================================
        SVAR1=2.D0/(1.D0+FRIC)
        SVAR2=1.D0-SVAR1
        SVAR3=DELTAT**2/TMASS/(1.D0+FRIC)
        CALL LIB__INVERTR8(3,T0,T0INV)
        TP=SVAR1*T0+SVAR2*TM+SVAR3*MATMUL(STRESS_I+KINSTRESS+STRESS_EXT,TRANSPOSE(T0INV))
print*,'==cell$propagate nconstraint ',nconstraint
print*,'==cell$propagate constrainttype ',constrainttype
print*,'==cell$propagate stress_I ',stress_I
print*,'==cell$propagate kinstress ',kinstress
print*,'==cell$propagate stress_ext ',stress_ext
print*,'==cell$propagate facts ',tmass,fric,deltat,svar1,svar2,svar3
print*,'==cell$propagate tm ',tm
print*,'==cell$propagate t0 ',t0
print*,'==cell$propagate tp ',tp
print*,'==cell$propagate ft ',MATMUL(STRESS_I+KINSTRESS+STRESS_EXT,TRANSPOSE(T0INV))
!       == CONSTRAINTS ========================================================
        IF(NCONSTRAINT.GT.0) THEN        
          ALLOCATE(B(NCONSTRAINT))
          ALLOCATE(X(NCONSTRAINT))
          ALLOCATE(A(NCONSTRAINT,NCONSTRAINT))
          SMAT1(:,:)=MATMUL(TP,T0INV)
          SMAT2(:,:)=MATMUL(TRANSPOSE(T0INV),T0INV)
          DO IC1=1,NCONSTRAINT
            B(IC1)=0.D0
            DO I=1,3
              DO J=1,3
                B(IC1)=B(IC1)+SMAT1(I,J)*CONSTRAINTPROJECT(I,J,IC1)
              ENDDO
            ENDDO
            SMAT3(:,:)=MATMUL(CONSTRAINTPROJECT(:,:,IC1),SMAT2(:,:))
            DO IC2=1,NCONSTRAINT
              A(IC1,IC2)=0.D0
              DO I=1,3
                DO J=1,3
                  A(IC1,IC2)=A(IC1,IC2)+SMAT3(I,J)*CONSTRAINTPROJECT(I,J,IC2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CALL LIB$MATRIXSOLVER8(NCONSTRAINT,NCONSTRAINT,1,A,X,B)
          DO IC1=1,NCONSTRAINT
            TP(:,:)=TP(:,:)-MATMUL(CONSTRAINTPROJECT(:,:,IC1),TRANSPOSE(T0INV))*X(IC1)
          ENDDO
          DEALLOCATE(B)
          DEALLOCATE(A)
          DEALLOCATE(CONSTRAINTPROJECT)
        END IF
!WRITE(*,FMT='(A10,3F10.5,5X,3F10.5,4X,3F10.5)')'TM',TM
!WRITE(*,FMT='(A10,3F10.5,5X,3F10.5,4X,3F10.5)')'T0',T0
!write(*,fmt='(a10,3f10.5,5x,3f10.5,4x,3f10.5)')'tp',tp
!write(*,fmt='(a10,3f10.5,5x,3f10.5,4x,3f10.5)')'stress',STRESS_I+KINSTRESS+STRESS_EXT
!       == CALCULATE KINETIC ENERGY 
        EKIN=0.5D0*TMASS*SUM((TP-TM)**2)/(2.D0*DELTAT)**2
      ELSE 
!       == EXTERNAL STRESS ===============================================
!        AMAT=MATMUL(SIGMA,TRANSPOSE(TREFINV))
!        AMAT=MATMUL(TREFINV,AMAT)
!   
        AMAT=MATMUL(SIGMA,TRANSPOSE(T0))
        AMAT=MATMUL(T0,AMAT)
        AMAT=AMAT*VREF       
        EPOT=0.5D0*(AMAT(1,1)+AMAT(2,2)+AMAT(3,3))
        AMAT=-AMAT
!       == EXTERNAL PRESSURE =============================================
        AMAT=AMAT-ONE*PRESSURE*V0  ! -PV - V_0 T0*SIGMA*T0^T
        EPOT=EPOT+PRESSURE*V0
!       == STRESS_I ======================================================
        AMAT=AMAT+STRESS_I+KINSTRESS  ! -DE/DALPHA -PV +V T0*SIGMA*T0^T
!       == STRESS PER VOLUME**2 ==========================================
        AMAT=AMAT/V0**2
!    
!       ==PROPAGATE ======================================================
        AMAT=AMAT/TMASS  ! ACCELERATION
        CALL LIB$INVERTR8(3,T0,T0INV)
        ALPHAM=MATMUL(TM,T0INV)-ONE
        ALPHA0=0.D0
        ALPHAP=-ALPHAM
        DO ITER=1,100
          ALPHADOT=(ALPHAP-ALPHAM)/(2.D0*DELTAT)
          BMAT=ALPHADOT-ONE*(ALPHADOT(1,1)+ALPHADOT(2,2)+ALPHADOT(3,3))
          BMAT=BMAT+TRANSPOSE(BMAT)
          BMAT=-BMAT*0.5D0*DELTAT  ! FRICTION CONSTANT
          XPMAT=(1.D0+FRIC)*ONE+BMAT         
          XMMAT=(1.D0-FRIC)*ONE-BMAT         
          CMAT=MATMUL(TRANSPOSE(ALPHADOT),ALPHADOT)
          CMAT=(CMAT-ONE*(CMAT(1,1)+CMAT(2,2)+CMAT(3,3)))  ! SIGN CHANGED. PB
          CALL LIB$INVERTR8(3,XPMAT,XPMAT)
          ALPHAP=MATMUL(2.D0*ALPHA0-MATMUL(ALPHAM,XMMAT) &
       &               +DELTAT**2*(AMAT+CMAT),XPMAT)
        END DO
        TP=MATMUL(ONE+ALPHAP,T0)
!   
!       ==  KINETIC ENERGY  ==============================================
        ALPHADOT=(ALPHAP-ALPHAM)/(2.D0*DELTAT)
        AMAT=MATMUL(ALPHADOT,TRANSPOSE(ALPHADOT))
        EKIN=0.5D0*TMASS*V0**2*(AMAT(1,1)+AMAT(2,2)+AMAT(3,3))
      end if
      TPROPAGATED=.TRUE.
      RETURN
      END SUBROUTINE CELL$PROPAGATE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$SWITCH()
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CELL_MODULE
      IMPLICIT NONE
!     ******************************************************************
      STRESS_I=0.D0
      EKIN=0.D0
      EPOT=0.D0
      IF(.NOT.TMOVE) RETURN
      TMM=TM
      TM=T0
      T0=TP
      TP=2.D0*T0-TM
      TP=3.D0*T0-3.D0*TM+TMM
!
!     == UPDATE VOUME OF THE UNIT CELL =================================
      V0=T0(1,1) * (T0(2,2)*T0(3,3)-T0(3,2)*T0(2,3)) &
     &  +T0(2,1) * (T0(3,2)*T0(1,3)-T0(1,2)*T0(3,3)) &
     &  +T0(3,1) * (T0(1,2)*T0(2,3)-T0(2,2)*T0(1,3))
      TPROPAGATED=.FALSE.
      RETURN
      END SUBROUTINE CELL$SWITCH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$READ(NFIL,NFILO,TCHK)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE RESTART_INTERFACE
      USE CELL_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)           ,INTENT(IN) :: NFIL
      INTEGER(4)           ,INTENT(IN) :: NFILO
      LOGICAL(4)           ,INTENT(OUT):: TCHK
      TYPE (SEPARATOR_TYPE),PARAMETER  :: MYSEPARATOR &
                 =SEPARATOR_TYPE(1,'CELL','NONE','MAY2000','NONE')
      TYPE (SEPARATOR_TYPE)            :: SEPARATOR
      INTEGER(4)                       :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      TCHK=.TRUE.
      SEPARATOR=MYSEPARATOR
      IF(THISTASK.EQ.1)CALL RESTART$READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      CALL MPE$BROADCAST('MONOMER',1,TCHK)
      IF(.NOT.TCHK) RETURN
!
!     ==================================================================
!     ==  READ DATA                                                   ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        IF(SEPARATOR%VERSION.NE.MYSEPARATOR%VERSION) THEN
          CALL ERROR$MSG('VERSION INCONSISTENCY')
          CALL ERROR$STOP('CELL$READ')
        END IF        
        READ(NFIL)T0,TM,TMM
      END IF
      CALL MPE$BROADCAST('MONOMER',1,T0)
      CALL MPE$BROADCAST('MONOMER',1,TM)
      CALL MPE$BROADCAST('MONOMER',1,TMM)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL$WRITE(NFIL,NFILO,TCHK)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CELL_MODULE
      USE RESTART_INTERFACE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN) :: NFIL
      INTEGER(4)            ,INTENT(IN) :: NFILO
      LOGICAL(4)            ,INTENT(OUT):: TCHK
      TYPE (SEPARATOR_TYPE),PARAMETER  :: MYSEPARATOR &
               =SEPARATOR_TYPE(1,'CELL','NONE','MAY2000','NONE')
      INTEGER(4)                        :: NTASKS,THISTASK
!     ******************************************************************
!
!     ==================================================================
!     ==  WRITE DATA                                                  ==
!     ==================================================================
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.EQ.1) THEN
        CALL RESTART$WRITESEPARATOR(MYSEPARATOR,NFIL,NFILO,TCHK)
        WRITE(NFIL)T0(:,:),TM(:,:),TMM(:,:)
      ELSE
        TCHK=.FALSE.   !WILL NOT BE USED, BUT MAKES THE COMPILER HAPPY.
      END IF
      RETURN
      END
