!**************************************************************************
!**************************************************************************
!****  OBJECT OPTFRIC                                                  ****
!****                                                                  ****
!****  1) -) CREATE INSTANCE WITH OPTFRIC$NEW                         ****
!****     -) SET STARTING FRICTION WITH OPTFRIC$SETR8('STARTFRIC',VAL) ****
!****        THE STARTIN FRICTION SHALL HAVE A VALUE BETWEEN ZERO AND  ****
!****        ONE                                                       ****
!****     -) SET RETARDATION: OPTFRIC$SETR8('RETARD',VAL)              ****
!****        THE OPTIMUM FRICTION HAS THE FORM                         ****
!****            SQRT(NUMERATOR/DENOMINATOR)                           ****
!****        THE RUNNING AVERAGE OF NUMERATOR AND DENOMINATOR IS USED. ****
!****        THE NEW VALUES ARE MIXED INTO AN AVERAGE VALUE            ****
!****        WITH WEIGHT RETARD. THUS RETARD=1 USES THE INSTANTANEOUS  ****
!****        BEST VALUE.                                               ****
!****                                                                  ****
!****  2) WHILE PROPAGATING:                                           ****
!****     -) SELECT INSTANCE WITH OPTFRIC$SELECT                      ****
!****     -) REQUEST OPTIMIZED FRICTION WITH OPTFRIC$GETR8('FRIC',VAL) ****
!****     -) USE  THIS VALUE FOR PROPAGATION OR OVERWRITE FRICTION     ****
!****        WITH OPTFRIC$SETR8('FRIC',VAL)                           ****
!****                                                                  ****
!****  3) DURING SWITCHING:                                            ****
!****     -) CALL OPTFRIC$UPDATER8 TO UPDATE FRICTIONS                ****
!****                                                                  ****
!**************************************************************************
!**************************************************************************
MODULE OPTFRIC_MODULE
TYPE OPTFRIC_TYPE
 CHARACTER(8)       :: ID
 logical(4)         :: ton
 REAL(8)            :: startFRIC
 REAL(8)            :: FRIC0
 REAL(8)            :: FRICM
 REAL(8)            :: RETARD
 REAL(8)            :: AOPTAV
 INTEGER(4)         :: LEN
 TYPE(OPTFRIC_TYPE),POINTER :: NEXT
END TYPE OPTFRIC_TYPE
LOGICAL           :: TINI=.FALSE.
TYPE(OPTFRIC_TYPE),POINTER:: FIRST_THIS
TYPE(OPTFRIC_TYPE),POINTER:: THIS
END MODULE OPTFRIC_MODULE
!
!     ..................................................................
      SUBROUTINE OPTFRIC$NEW(ID)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE OPTFRIC_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
!     ******************************************************************
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        ALLOCATE(FIRST_THIS)
        THIS=>FIRST_THIS
      ELSE
        THIS=>FIRST_THIS
        DO 
          IF(THIS%ID.EQ.ID) THEN
            CALL ERROR$MSG('CANNOT CREATE OPTFRIC INSTANCE WITH THE SAME NAME')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('OPTFRIC$NEW')
          END IF
          IF(ASSOCIATED(THIS%NEXT)) THEN
            THIS=>THIS%NEXT
          ELSE
            ALLOCATE(THIS%NEXT)
            THIS=>THIS%NEXT
            EXIT
          END IF
        ENDDO
      END IF
!   
      THIS%ID      =ID
      THIS%len=0
      THIS%aoptav=0.d0
      NULLIFY(THIS%NEXT)
      RETURN
      END 
!
!     ..................................................................
      SUBROUTINE OPTFRIC$DELETE
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE OPTFRIC_MODULE
      IMPLICIT NONE
      TYPE(OPTFRIC_TYPE),POINTER :: PREVIOUS
!     ******************************************************************
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO OPTFRIC INSTANCE SELECTED')
        CALL ERROR$STOP('OPTFRIC$DELETE')
      END IF
      IF(ASSOCIATED(THIS,FIRST_THIS)) THEN
        FIRST_THIS=THIS
        DEALLOCATE(THIS)
      ELSE
        PREVIOUS=>FIRST_THIS
        DO WHILE (.NOT.ASSOCIATED(PREVIOUS,THIS))
          PREVIOUS=>PREVIOUS%NEXT
        END DO
        PREVIOUS%NEXT=>THIS%NEXT
        DEALLOCATE(THIS)
      END IF 
      RETURN
      END 
!
!     ..................................................................
      SUBROUTINE OPTFRIC$SELECT(ID)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE OPTFRIC_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
!     ******************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('CREATE OPTFRIC INSTANCES BEFORE SELECTING')
        CALL ERROR$STOP('OPTFRIC$SELECT')
      END IF
!
!     == UNSELECT ======================================================
      IF(ID.EQ.'NONE') THEN
        NULLIFY(THIS)
      END IF
!
!     == FIND OPTFRIC INSTANCE AND SELECT
      THIS=>FIRST_THIS
      DO WHILE (THIS%ID.NE.ID)
        IF(.NOT.ASSOCIATED(THIS%NEXT)) THEN
          CALL ERROR$MSG('OPTFRIC INSTANCE DOES NOT EXIST')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('OPTFRIC$SELECT')
        END IF
        THIS=>THIS%NEXT
      ENDDO
      RETURN
      END 
!
!     .................................................................. 
      SUBROUTINE OPTFRIC$SETR8(ID,VAL)
!     ******************************************************************
!     **  CALCULATES THE ABSOLUTE VALUE OF THE FORCE                  **
!     **                                                              **
!     **  CALL AFTER CONSTRAINTS AND BEFORE SWITCH                    **
!     ****************************************************************** 
      use optfric_module
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)    ,INTENT(IN)  :: VAL
!     *******************************************************************
      IF(ID.EQ.'STARTFRIC') THEN
        this%STARTFRIC=VAL
        this%fric0=val
        this%fricm=val
      ELSE IF(ID.EQ.'RETARD') THEN
!       == VAL IS THE NUMBER OF STEPS FOR WHICH THE DEVIATION FROM THE 
!       == VALUE FALLS TO 0.5 OF THE INITIAL DEVIATION. THE INTERNAL 
!       == VALUE THIS%RETARD IS THE MIXING FACTOR                    
        THIS%RETARD=1.D0-1.D0/2.D0**(1.D0/VAL)
      ELSE IF(ID.EQ.'FRIC') THEN
        THIS%FRIC0=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$CHVAL('THIS%ID',THIS%ID)
        CALL ERROR$STOP('OPTFRIC$SETR8')
      END IF
      RETURN 
      end
!
!     .................................................................. 
      SUBROUTINE OPTFRIC$gETR8(ID,VAL)
!     ******************************************************************
!     **  CALCULATES THE ABSOLUTE VALUE OF THE FORCE                  **
!     **                                                              **
!     **  CALL AFTER CONSTRAINTS AND BEFORE SWITCH                    **
!     ****************************************************************** 
      use optfric_module
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(out):: VAL
!     *******************************************************************
      IF(ID.EQ.'FRIC') THEN
        val=this%FRIC0
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$CHVAL('THIS%ID',THIS%ID)
        CALL ERROR$STOP('OPTFRIC$gETR8')
      END IF
      RETURN 
      end
!
!     .................................................................. 
      SUBROUTINE OPTFRIC$SETl4(ID,VAL)
!     ******************************************************************
!     **  CALCULATES THE ABSOLUTE VALUE OF THE FORCE                  **
!     **                                                              **
!     **  CALL AFTER CONSTRAINTS AND BEFORE SWITCH                    **
!     ****************************************************************** 
      USE OPTFRIC_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)    ,INTENT(IN)  :: VAL
!     *******************************************************************
      IF(ID.EQ.'ON') THEN
        THIS%tON=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$CHVAL('THIS%ID',THIS%ID)
        CALL ERROR$STOP('OPTFRIC$SETL4')
      END IF
      RETURN 
      END
!
!     .................................................................. 
      SUBROUTINE OPTFRIC$gETl4(ID,VAL)
!     ******************************************************************
!     **  CALCULATES THE ABSOLUTE VALUE OF THE FORCE                  **
!     **                                                              **
!     **  CALL AFTER CONSTRAINTS AND BEFORE SWITCH                    **
!     ****************************************************************** 
      USE OPTFRIC_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(out):: VAL
!     *******************************************************************
      IF(ID.EQ.'ON') THEN
        val=THIS%tON
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$CHVAL('THIS%ID',THIS%ID)
        CALL ERROR$STOP('OPTFRIC$gETL4')
      END IF
      RETURN 
      END
!
!     .................................................................. 
      SUBROUTINE OPTFRIC$UPDATER8(ID,LEN,XP,X0,XM,X2M,MASS)
!     ******************************************************************
!     **  CALCULATES THE ABSOLUTE VALUE OF THE FORCE                  **
!     **                                                              **
!     **  CALL AFTER CONSTRAINTS AND BEFORE SWITCH                    **
!     ****************************************************************** 
      use optfric_module
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(IN) :: XP(LEN)  ! POSITIONS AT NEXT TIME STEP
      REAL(8)     ,INTENT(IN) :: X0(LEN)  ! POSITIONS AT CURRENT TIME STEP
      REAL(8)     ,INTENT(IN) :: XM(LEN)  ! POSITION AT PREVIOUS TIME STEP
      REAL(8)     ,INTENT(IN) :: X2M(LEN) ! POSITION TWO STEPS BEFORE
      REAL(8)     ,INTENT(IN) :: MASS(len)! DIAGONAL ELEMENTS OF MASS TENSOR
      REAL(8)                 :: ASUM,BSUM
      REAL(8)                 :: A0 ! EFFECTIVE FRICTION/CURRENT TIME STEP
      REAL(8)                 :: AM ! EFFECTIVE FRICTION/PREVIOUS TIME STEP
      real(8)                 :: aopt1 !numerator
      real(8)                 :: aopt2 !denominator
      real(8)                 :: aopt
      REAL(8)                 :: DF,DX,MDX
      REAL(8)                 :: mixaopt
      integer(4)              :: i
!     ****************************************************************** 
      CALL OPTFRIC$SELECT(ID)      
      IF(.NOT.THIS%TON) RETURN
      IF(THIS%LEN.EQ.0) THIS%LEN=LEN
      IF(LEN.NE.THIS%LEN) THEN
        CALL ERROR$MSG('SIZE INCONSISTENT')
        CALL ERROR$CHVAL('THIS%ID',THIS%ID)
        CALL ERROR$STOP('OPTFRIC$UPDATER8')
      END IF
!
!     =======================================================================
!     == ESTIMATE OPTIMUM FRICTION FROM FOUR TIME STEPS AND MASSES         ==
!     =======================================================================
      A0=THIS%FRIC0
      AM=THIS%FRICM
      ASUM=0.D0
      BSUM=0.D0
      DO I=1,LEN
        DF=(1.D0+A0)*XP(I)-(3.D0+AM)*X0(I)+(3.D0-A0)*XM(I)-(1.D0-AM)*X2M(I)
        DX=X0(I)-XM(I)
        MDX=MASS(I)*DX
        ASUM=ASUM-DF*MDX
        BSUM=BSUM+DX*MDX
      ENDDO
      AOPT1=ASUM
      AOPT2=BSUM
!
!     =======================================================================
!     == set new friction                                                  ==
!     =======================================================================
      this%fricm=this%fric0
      IF(aopt1.gT.0.D0) THEN
        aopt=sqrt(AOPT1/AOPT2)
        mixaopt=this%retard
        this%aoptav=MIXAOPT*aopt+(1.D0-MIXAOPT)*this%AOPTAV
        this%fric0=this%aoptav
        this%fric0=min(this%fric0,1.d0)
        this%fric0=max(this%fric0,0.d0)
      ELSE
        this%fric0=this%startfric
      END IF
!print*,'optfric ',mixaopt,this%aoptav,this%fric0,AOPT1,AOPT2
      RETURN
      END
!
!     .................................................................. 
      SUBROUTINE OPTFRIC$TESTCONV(LEN,XP,X0,XM,X2M,A0,AM,DT,MASS,DX,DE)
!     ******************************************************************
!     **  ESTIMATES THE DEVIATION IN ENERGY AND DISTANCE FROM THE     **
!     **  MINIMUM IF THE TOTAL ENERGY SURFACE                         **
!     **  this is a stand-alone routine                               **
!     ****************************************************************** 
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(IN) :: XP(LEN)  ! NEXT POSITIONS 
      REAL(8)     ,INTENT(IN) :: X0(LEN)  ! CURRENT POSITIONS 
      REAL(8)     ,INTENT(IN) :: XM(LEN)  ! PREVIOUS POSITIONS
      REAL(8)     ,INTENT(IN) :: X2M(LEN) ! POSITIONS BEFORE PREVIOUS 
      REAL(8)     ,INTENT(IN) :: A0       ! CURRENT FRICTION FACTOR 
      REAL(8)     ,INTENT(IN) :: AM       ! PREVIOUS FRICTION FACTOR
      REAL(8)     ,INTENT(IN) :: DT       ! TIME STEP
      REAL(8)     ,INTENT(IN) :: MASS(LEN)! DIAGONAL ELEMENTS OF MASS TENSOR
      REAL(8)     ,INTENT(OUT):: DX       ! ESTIMATED DISTANCE FROM MINIMUM
      REAL(8)     ,INTENT(OUT):: DE       ! ESTIMATED ENERGY DIFFERENCE
      REAL(8)                 :: F02,DXDF,DX2
!     ****************************************************************** 
!     == F02=F0**2 =========================================================
      F02=SUM((MASS(:)*( (1.D0+A0)*XP(:) - 2.D0*X0(:) + (1.D0-A0)*XM(:) ))**2)
      F02=F02/DT**4
!     == DXDF=(X0-XM)*(F0-FM) ================================================
      DXDF=DOT_PRODUCT(X0(:)-XM(:),MASS(:)*((1.D0+A0)*XP(:)-(3.D0+AM)*X0(:) &
     &                            +(3.D0-A0)*XM-(1.D0-AM)*X2M(:)))
      DXDF=DXDF/DT**2
!     == DX2=(X0-XM)**2 ======================================================
      DX2=SUM((X0(:)-XM(:))**2)
!     == avoid divide by zero for x0=xm ====================================
      if(dxdf.eq.0.d0) then
        de=huge(de)
        dx=huge(dx)
        return
      end if
!     == DE IS THE ESTIMATE OF THE ENERGY DIFFERENCE TO THE MINIMUM
      DE=-0.5D0*DX2*F02/DXDF
!     == DX IS THE ESTIMATED DISTANCE FROM THE MINIMUM
      DX=-DX2*SQRT(F02)/DXDF
      RETURN
      END
!!$!
!!$!      .......................................................................
!!$       PROGRAM TEST
!!$       IMPLICIT NONE
!!$       REAL(8)     :: C(2,2)
!!$       REAL(8)     :: MASS(2)
!!$       REAL(8)     :: DT=1.D-1
!!$       REAL(8)     :: X0(2)
!!$       REAL(8)     :: XM(2)
!!$       REAL(8)     :: XP(2)
!!$       REAL(8)     :: X2M(2)
!!$       REAL(8)     :: F(2)   ! FORCE
!!$       REAL(8)     :: V(2)   ! VELOCITY
!!$       REAL(8)     :: SVAR1,SVAR2,SVAR3
!!$       REAL(8)     :: anne
!!$       REAL(8)     :: EPOT,EKIN
!!$       INTEGER(4)  :: ITER
!!$       INTEGER(4),PARAMETER :: NITER=500
!!$!      ********************************************************************
!!$       C(:,1)=(/8.D0,0.D0/)
!!$       C(:,2)=(/0.D0,1.D0/)
!!$       X0(:)=(/5.D0,1.D0/)
!!$       XM(:)=X0(:)
!!$       X2M(:)=XM(:)
!!$       MASS(:)=1.D0
!!$!
!!$       CALL OPTFRIC$NEW('TEST')
!!$       CALL OPTFRIC$SELECT('TEST')
!!$       CALL OPTFRIC$SETR8('STARTFRIC',0.1D0)
!!$       CALL OPTFRIC$SETR8('RETARD',1.D0)  
!!$!
!!$       DO ITER=1,NITER
!!$         CALL OPTFRIC$GETR8('FRIC',ANNE)
!!$!anne=0.d0
!!$         CALL OPTFRIC$sETR8('FRIC',ANNE)
!!$         F(:)=-MATMUL(C,X0)
!!$         EPOT=-0.5D0*DOT_PRODUCT(X0,F)
!!$         SVAR1=2.D0/(1.D0+ANNE)
!!$         SVAR2=1.D0-SVAR1  
!!$         SVAR3=DT**2/(1+ANNE)
!!$         XP=SVAR1*X0+SVAR2*XM+SVAR3*F/MASS
!!$         V=0.5D0/DT*(XP-XM)
!!$         EKIN=0.5D0*SUM(MASS*V**2)
!!$         WRITE(*,FMT='(I5,4F20.10)')ITER,EKIN,EPOT,EKIN+EPOT,anne
!!$!
!!$         CALL OPTFRIC$UPDATER8('TEST',2,XP,X0,XM,X2M,MASS) 
!!$         X2M=XM
!!$         XM=X0
!!$         X0=XP       
!!$       ENDDO
!!$       STOP
!!$       END
!!$!
!!$!..............................
!!$subroutine error$msg(id)
!!$character(*),intent(in) :: id
!!$print*,'error',id
!!$return
!!$end
!!$!
!!$!..............................
!!$subroutine error$chval(id,val)
!!$character(*),intent(in) :: id
!!$character(*),intent(in) :: val
!!$print*,'error: variable ',trim(id),' has value ',trim(val)
!!$return
!!$end
!!$!
!!$!..............................
!!$subroutine error$stop(id)
!!$character(*),intent(in) :: id
!!$print*,'error stop in ',trim(id)
!!$stop
!!$end
