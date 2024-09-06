MODULE DIMER_MODULE
!USE DIMER_CONSTR_MODULE
REAL(8)                       :: D ! CONSTRAINED DISTANCE BETWEEN DIMERPOINTS
INTEGER                       :: CALCMULTIPLIERITERMAX ! EMERGENCY EXIT FOR THE ITERATION LOOP AFTER X STEPS
REAL(8)                       :: DLAMBDA !EXACTNESS OF LAMBDA (= LAGRANGE MULTIPLIER) CALCULATION
INTEGER                       :: CALCVELOCITYITERMAX !EMERGENCY EXIT FOR THE ITERATION LOOP AFTER X STEPS
REAL(8)                       :: DVELOCITY !EXACTNESS OF VELOCITY ITERATION 
REAL(8)                       :: ENERGYTRA !WRITE THE ENERGY/DISTANCE EVERY N STEPS
REAL(8)                       :: LASTTRA=0.0D0
REAL(8),ALLOCATABLE            :: RTS(:)      !COORDINATES OF THE TRANSITION STATE FROMCNTL-FILE               
REAL(8),ALLOCATABLE            :: DIRM(:)

INTEGER(4)                    :: LCS       !LENGTHCOUNTERSHOULDBE HOW MANY STEPS WITH NO CHANGE OF DIMERCENTER BEFORE SHORTEN SQD
INTEGER(4)                    :: LC =0     ! THE ACTUAK STEPS WITHOUT CHANGE
REAL(8)                       :: RCDIFFMIN !DIFFERENCE IN DIMER CENTER POSITION BELOW LC=LC+1
REAL(8)                       :: RCDIFF    !DIFFERENCE IN DIMER CENTER POSITION
REAL(8),ALLOCATABLE,SAVE            :: RC(:)        !DIMER CENTER POSITION
REAL(8),ALLOCATABLE,SAVE            :: RCOLD(:)     !DIMER CENTER POSITION IN LAST STEP
REAL(8)                       :: DSTEP     !D=D-DSTEP
INTEGER(4)                    :: NSTEPS
REAL(8)                       :: STEPFACT
LOGICAL(4)                    :: TINITSTEPFACT=.FALSE.
REAL(8)                       :: DMIN    !THE MINIMAL DIMER LENGTH (NO FURTHER REDUCEMENT OF D)
REAL(8)                       :: WDOWNFACT
REAL(8),ALLOCATABLE           :: F1W(:) !THE WEIGHTED FORCE
LOGICAL(4)                    :: FIRSTTRA=.TRUE.
LOGICAL(4)                    :: DIMER !DO WE USE DIMER PARALLELISATION? SET IN READIN
LOGICAL(4)                    :: DIMERFOLLOWDOWN !DO WE FOLLOW DOWN WITH EACH IMAGE TO ONE OF THE VICINAL GROUND STATES

LOGICAL(4)                    :: DLFLEX    !FLEXIBLE DIMER LENGTH
LOGICAL(4)                    :: KDLENGTH    !KEEP THE STARTUP LENGTH OF THE DIMER
LOGICAL(4)                    :: INHIBITUP   !INHIBIT UPWARD MOTION OF THE DIMER
LOGICAL(4)                    :: INHIBITPERP
LOGICAL(4)                    :: ONLYPERP
LOGICAL(4)                    :: ONLYROT
LOGICAL(4)                    :: WDOWN       !WEIGHT DOWN THE DIMER IMAGE 1

CHARACTER(32)                 :: CENTER_ID
REAL(8)                       :: CENTER_COORD(3)
REAL(8)                       :: DROT

INTEGER                       :: DIM ! DIMONSION =NAT*3
REAL(8),SAVE                  :: SQD     !SQUARE CONSTRAINED DISTANCE BETWEEN DIMERPOINTS
REAL(8)                       :: SQDR=0.0D0     !SQUARE CONSTRAINED DISTANCE BETWEEN DIMERPOINTS
REAL(8)                       :: ALPHA !FRICTION; GOOD:0.0156D0
LOGICAL(4)                    :: PLACEDIMER=.FALSE.

LOGICAL(4)                    :: TFIRSTSPLIT=.TRUE.
REAL(8),ALLOCATABLE           :: RKEEP(:) 
REAL(8),ALLOCATABLE           :: R2SPLIT(:) 
LOGICAL(4)                    :: STRETCH
REAL(8)                       :: STRETCHDIST


REAL(8)                       :: FMPARA,FMPERP,FMROT,FRICPARA,FRICPERP,FRICROT
LOGICAL(4)                    :: OPTFRICPARA,OPTFRICPERP,OPTFRICROT
LOGICAL(4)                    :: FAUTOPARA,FAUTOPERP,FAUTOROT
REAL(8)                       :: FRICAUTOPARA,FRICAUTOPERP,FRICAUTOROT
REAL(8),ALLOCATABLE           :: FROTM(:),FPERPM(:),FPARAM(:)
REAL(8)                       :: CONSTRSTEP !MAX STEP FOR CONSTRAINT
REAL(8)                       :: EKINPARAM
REAL(8)                       :: EKINPERPM
REAL(8)                       :: EKINROTM 
REAL(8)                       :: TPARAM
REAL(8)                       :: TPERPM
REAL(8)                       :: TROTM 
REAL(8)                       :: TFACT=2.D0/(3.166679D-6) !PER DEG. OF FREEDOM; TO DO: GET K FROM PAW CONSTANTS

!========================================================
  TYPE DIMER_CONSTR_TYPE
     CHARACTER(32)                    :: ID        !KEYWORD
     TYPE(DIMER_CONSTR_TYPE),POINTER  :: NEXT_D      !YOUNGER BROTHER
     TYPE(DIMER_CONSTR_TYPE),POINTER  :: PREV_D      !ELDER BROTHER
  END TYPE DIMER_CONSTR_TYPE

  TYPE(DIMER_CONSTR_TYPE),POINTER  :: ELDEST_D
  TYPE(DIMER_CONSTR_TYPE),POINTER  :: THIS_D
!========================================================



  !ALEX: THIS IS FOR PCLIMB-TESTING!
  LOGICAL(4)               :: CLIMBPERP,TFIRSTCLIMBSTEP
  REAL(8)                  :: FORCEDSTEP
  REAL(8), ALLOCATABLE     :: PCLIMBDIR(:)
  !ALEX: THIS IS FOR PCLIMB-TESTING!
  
  !ALEX: THIS IS FOR ANGLE MONITORING
  REAL(8),ALLOCATABLE      :: ANGLE1(:)
  REAL(8),ALLOCATABLE      :: ANGLE2(:)
  REAL(8),ALLOCATABLE      :: ANGLEDIR(:)
  LOGICAL(4)               :: TREADAM=.FALSE.
  LOGICAL(4)               :: TREADAMNOTTHERE=.FALSE.
  !ALEX: THIS IS FOR ANGLE MONITORING
  
  !ALEXP: THIS COMES FROM MERGE IN NEW CODE
  INTEGER(4)               :: NTASKS
  INTEGER(4)               :: THISTASK
  !ALEXP: THIS COMES FROM MERGE IN NEW CODE
  
  INTEGER(4)               :: DPROTFIL

  REAL(8),ALLOCATABLE      :: G1(:)
  REAL(8),ALLOCATABLE      :: G2(:)
  INTEGER(4)               :: STEPS=0

  REAL(8)                  :: ETOT,ETOT2 !FOR TS ESTIMATION
END MODULE DIMER_MODULE

MODULE DIMER_OSCILLATOR_MODULE
  INTEGER(4)                :: ODIM=3    !FOR THE DIMER: CPARA,CPERP,CROT
  REAL(8),ALLOCATABLE       :: OSCM(:)
  REAL(8),ALLOCATABLE       :: OSC0(:)
  REAL(8),ALLOCATABLE       :: OSCP(:)
  REAL(8),ALLOCATABLE       :: OSCMASS(:)  
  REAL(8),ALLOCATABLE       :: OSCANNER(:)
  REAL(8),ALLOCATABLE       :: OSCC(:) !THE HARMONIC POTENTIAL
END MODULE DIMER_OSCILLATOR_MODULE

MODULE DIMER_ESTIMATETS_MODULE
  REAL(8),ALLOCATABLE       :: Z1M(:)                  
  REAL(8),ALLOCATABLE       :: Z2M(:)                  
  REAL(8),ALLOCATABLE       :: EM(:)                  
  REAL(8),ALLOCATABLE       :: EORTHOM(:)                  
  REAL(8)                   :: FPARAM
  REAL(8)                   :: FPERPM
END MODULE DIMER_ESTIMATETS_MODULE




SUBROUTINE DIMER_INIT_FILES()
  USE STRINGS_MODULE
  USE DIMER_MODULE
  IMPLICIT NONE
  CHARACTER(32)                    :: ID


  !=== DIMERPROTOCOL FILE =================================================
  !=== EACH MONOMER HAS ITS OWN PROTOCOL FILE
  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
  ID=+'DPROT'
  IF(THISTASK.GT.1)THEN
     CALL FILEHANDLER$SETFILE(ID,.FALSE.,-'/DEV/NULL')
  ELSE
     CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.DIMER_PROT')
  END IF
  CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
  CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
  CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
  CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')

  CALL FILEHANDLER$UNIT('DPROT',DPROTFIL)



  !=== DIMERFOLLOWDON ENERGY TRAJECTORY FILE ==============================
  !=== EACH MONOMER HAS ITS OWN FILE
  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
  ID=+'ETRA'
  IF(THISTASK.GT.1)THEN
     CALL FILEHANDLER$SETFILE(ID,.FALSE.,-'/DEV/NULL')
  ELSE
     CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.DIMER_ETRA')
  END IF
  CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
  CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
  CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
  CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')



  
  RETURN
END SUBROUTINE DIMER_INIT_FILES




!     ..................................................................
      SUBROUTINE DIMER$PROPAGATE(NAT,DT,ANNER,ANNEE &
     &                 ,RMASS,EFFEMASS,FORCE,R0,RM,RP &
     &                 ,TSTRESS,CELLFRIC,CELLKIN)
!     ******************************************************************
!     **  COMMENT                         **
!     ******************************************************************
      USE MPE_MODULE
      USE DIMER_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT           ! #(ATOMS)
      REAL(8)   ,INTENT(IN) :: DT        ! TIME STEP
      REAL(8)   ,INTENT(IN) :: ANNER         ! FRICTION ON THE ATOMS
      REAL(8)   ,INTENT(IN) :: ANNEE         ! FRICTION ON THE WAVE FUNCTIONS
      REAL(8)   ,INTENT(IN) :: RMASS(NAT)    ! BARE MASS OF THE NUCLEUS
      REAL(8)   ,INTENT(IN) :: EFFEMASS(NAT) ! MASS OF TEH ELECTRONS
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)     ! R(0)
      REAL(8)   ,INTENT(IN) :: RM(3,NAT)     ! R(-)
      REAL(8)   ,INTENT(OUT):: RP(3,NAT)     ! R(+)
      REAL(8)   ,INTENT(IN) :: FORCE(3,NAT)  ! FORCE
      LOGICAL(4),INTENT(IN) :: TSTRESS
      REAL(8)   ,INTENT(IN) :: CELLFRIC(3,3) 
      REAL(8)   ,INTENT(OUT):: CELLKIN(3,3) 
      REAL(8)               :: CELLKIN1(3,3),CELLKIN2(3,3)  
      REAL(8)               :: MM(NAT*3,NAT*3) !MASSMATRIX
      REAL(8)               :: R1_(NAT*3),R2_(NAT*3) !OLD POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)               :: R1(NAT*3),R2(NAT*3) !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)               :: R1P(NAT*3),R2P(NAT*3) !NEW POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)               :: F1(NAT*3),F2(NAT*3) !FORCE FROM POTENTIAL
      INTEGER(4)            :: NVAL,WORLD_ID
      INTEGER(4)            :: IAT

      CELLKIN=0.D0
      CELLKIN1=0.0D0
      CELLKIN2=0.0D0


      ALPHA=ANNER!DO WE NEED THIS??? DONT THINK SO! CHECK IT!
      IF(.NOT.TSTRESS) THEN
         CALL MPE$QUERY('~',NTASKS,WORLD_ID)
         CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)

         CALL DIMER$GETPROCESSLEADER2(NVAL,WORLD_ID)
         !NVAL NOW EQUALS THE WORLD_ID OF THE 1ST TASK IN DIMER2
         !WORLD_ID EQUALS THE WORLD_ID OF THE CALLING PROCESS



         IF(THISTASK.EQ.1) THEN
            !GET TOTAL ENERGY FROM LAST STEP FOR TS ESTIMATE
            CALL ENERGYLIST$GET('TOTAL ENERGY',ETOT)
            
            !ONLY 1ST TASK IN DIMER HAS TO DO THE FOLLOWING
            DIM=NAT*3 !WE NEED THIS TO INITIALIZE DIM BECAUSE NAT IS READ AFTER DIMER!
            MM(:,:)=0.0D0
            DO IAT=1 , NAT
               MM((IAT*3)-2,(IAT*3)-2)=RMASS(IAT)-EFFEMASS(IAT)
               MM((IAT*3)-1,(IAT*3)-1)=RMASS(IAT)-EFFEMASS(IAT)
               MM((IAT*3),(IAT*3))    =RMASS(IAT)-EFFEMASS(IAT)
               
               R1_((IAT*3)-2)         =RM(1,IAT)
               R1_((IAT*3)-1)         =RM(2,IAT)
               R1_((IAT*3))           =RM(3,IAT)
               
               R1((IAT*3)-2)          =R0(1,IAT)
               R1((IAT*3)-1)          =R0(2,IAT)
               R1((IAT*3))            =R0(3,IAT)
               
               F1((IAT*3)-2)          =FORCE(1,IAT)
               F1((IAT*3)-1)          =FORCE(2,IAT)
               F1((IAT*3))            =FORCE(3,IAT)
            END DO




            !GIVE WORLD_ID 1 THE POSITION AND FORCES ON DIMER2
            IF(WORLD_ID.NE.1) THEN 
               !WE ARE 1ST TASK OF DIMER2
               CALL MPE$SEND('~',1,42,R1_)
               CALL MPE$SEND('~',1,43,R1)
               CALL MPE$SEND('~',1,44,F1)
               CALL MPE$SEND('~',1,45,ETOT)
               
            ELSE
               !WE ARE WORLD_ID 1
               CALL MPE$RECEIVE('~',NVAL,42,R2_)
               CALL MPE$RECEIVE('~',NVAL,43,R2)
               CALL MPE$RECEIVE('~',NVAL,44,F2)
               CALL MPE$RECEIVE('~',NVAL,45,ETOT2)
            END IF


            
            IF(WORLD_ID.EQ.1) THEN
               !ONLY THE FIRST TASK MOVES THE DIMER IMAGES
               SQD=D*D !REMEMBER: SQDIMERDIST IS THE *ACTUAL* DISTANCE**2, SQD=D*D IS THE DESIRED DISTANCE**2 
               CALL DIMER_PROPAGATE(DT,R1,R2,R1_,R2_,MM,F1,F2,R1P,R2P,CELLKIN1,CELLKIN2)
               CELLKIN=CELLKIN1
            END IF

      
         ELSE
            !ALL TASKS .NE.1 IN DIMER
         END IF




  !FROM NOW ON WE HAVE ALL TASKS BACK

         !SEND 1ST TASK IN DIMER2 THE NEW POSITION AND CELLKIN OF IMAGE 2
         IF(WORLD_ID.EQ.1) THEN 
            !WE ARE WORLD_ID 1
            CALL MPE$SEND('~',NVAL,45,R2P)
            CALL MPE$SEND('~',NVAL,46,CELLKIN2)
         ELSE IF(WORLD_ID.EQ.NVAL) THEN
            !WE ARE 1ST TASK IN DIMER2
            CALL MPE$RECEIVE('~',1,45,R1P)
            CALL MPE$RECEIVE('~',1,46,CELLKIN)
         END IF

         IF(THISTASK.EQ.1) THEN
            !WE ARE 1ST TASK IN DIMER GROUP
            
            DO IAT=1 , NAT
               RP(1,IAT)=R1P((IAT*3)-2)
               RP(2,IAT)=R1P((IAT*3)-1)
               RP(3,IAT)=R1P((IAT*3))    
            END DO
         END IF

         !SEND THE NEW POSITION, CELLKIN TO ALL TASKS IN THE GROUP
         CALL MPE$BROADCAST('MONOMER',1,RP)
         CALL MPE$BROADCAST('MONOMER',1,CELLKIN)


      ELSE 
         
         CALL ERROR$MSG('IN DIMER PROPAGATE TSTRESS=TRUE TALK ABOUT THIS WITH PETER')
         CALL ERROR$STOP('DIMER$PROPAGATE')
         
      END IF

      RETURN
      
!     ******************************************************************

      END SUBROUTINE DIMER$PROPAGATE


!     ..................................................................
      SUBROUTINE DIMER_PROPAGATE(DT,R1,R2,R1_,R2_,MM,F1,F2,R1P,R2P,CELLKIN1,CELLKIN2)
!     ******************************************************************
!     **  COMMENT                         **
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)    :: DT                    ! TIME STEP
      REAL(8)   ,INTENT(IN)    :: R1(DIM)               !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,INTENT(IN)    :: R2(DIM)               !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,INTENT(IN)    :: R1_(DIM)              !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,INTENT(IN)    :: R2_(DIM)              !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,INTENT(IN)    :: MM(DIM,DIM)           !MASSMATRIX
      REAL(8)   ,INTENT(IN)    :: F1(DIM)               !FORCE FROM POTENTIAL
      REAL(8)   ,INTENT(IN)    :: F2(DIM)               !FORCE FROM POTENTIAL

      REAL(8)   ,INTENT(OUT)   :: R1P(DIM)              !NEW POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,INTENT(OUT)   :: R2P(DIM)              !NEW POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,INTENT(OUT)   :: CELLKIN1(3,3)         ! 
      REAL(8)   ,INTENT(OUT)   :: CELLKIN2(3,3)         ! 
      REAL(8)                  :: UF1(DIM)              !M^-(1/2)*FORCE FROM POTENTIAL
      REAL(8)                  :: UF2(DIM)              !M^-(1/2)*FORCE FROM POTENTIAL
      REAL(8)                  :: MF1(DIM)              !M^(1/2)*FORCE FROM POTENTIAL
      REAL(8)                  :: MF2(DIM)              !M^(1/2)*FORCE FROM POTENTIAL
      REAL(8)                  :: SQDIMERDIST           !CURRENT (DISTANCE BETWEEN DIMERPOINTS)**2
      REAL(8)                  :: R1NC(DIM),R2NC(DIM)   !POSITION WITHOUT CONSTRAINS
      REAL(8)                  :: V1(DIM),V2(DIM)                !VELOCITY
      REAL(8)                  :: X1(DIM),X2(DIM)  !THE MASSWEIGHTED COORDINATES
      REAL(8)                  :: X1M(DIM),X2M(DIM),X1P(DIM),X2P(DIM)
      REAL(8)                  :: RMASS(DIM/3)
      REAL(8)                  :: RTS_(DIM) !FOR TS ESTIMATION
      REAL(8)                  :: RTSOLD(DIM) !FOR TS ESTIMATION (THIS IS DIRTY BECAUSE IT SHOULD BE SAVED!)
      REAL(8)                  :: FRICUSEDPERP,FRICUSEDPARA,FRICUSEDROT
      REAL(8)                  :: FROT0(DIM)  !THE ROTATIONAL PART OF THE FORCES
      REAL(8)                  :: FPERP0(DIM) !THE PERP. PART OF THE FORCES
      REAL(8)                  :: FPARA0(DIM) !THE PERP. PART OF THE FORCES
      REAL(8)                  :: SVARV(DIM),SVAR
      INTEGER(4)               :: I,J,IAT
      LOGICAL(4)               :: PERPFIRST,PARAFIRST,ROTFIRST
      REAL(8)                  :: FP1,FP2
      REAL(8)                  :: A,B,C,D_
      !******************************************************************

      CELLKIN1=0.0D0
      CELLKIN2=0.0D0
      R1P=0.0D0
      R2P=0.0D0
      IF(.NOT.ALLOCATED(RC)) ALLOCATE(RC(DIM))
      IF(.NOT.ALLOCATED(RCOLD)) ALLOCATE(RCOLD(DIM))
      IF(.NOT.ALLOCATED(F1W)) ALLOCATE(F1W(DIM))
      
      !========================================
      !== INIT THE MASSWEIGHTED COORDINATES  ==  
      !========================================
      CALL DIMER$GET_MASSWEIGHTED(DIM,R1,X1)
      CALL DIMER$GET_MASSWEIGHTED(DIM,R2,X2)
!REMEMBER: SQDIMERDIST IS THE *ACTUAL* DISTANCE**2, SQD=D*D IS THE DESIRED DISTANCE**2
      SQDIMERDIST=DOT_PRODUCT((X1-X2),(X1-X2)) 
      WRITE(DPROTFIL,*)"DIMER PROPAGATE: DIMERDISTANCE= ",SQRT(SQDIMERDIST) &
     &                ,"  UNMW= ",SQRT(DOT_PRODUCT(R1-R2,R1-R2))


      !===============================================
      !==========    PLACE THE DIMER  ================
      !===============================================
      !THIS SHOULD ONLY BE DONE IN THE FIRST ITERATION STEP
      IF(PLACEDIMER.AND..NOT.KDLENGTH) THEN !IN THE FIRST STEP  
         CALL PLACE_DIMER(D,X1,X2)!PLACE IT SO THAT LENGTH-CONSTRAIN FOR THE MASSWEIGHTED COORD.IS FULFILLED
         !THE COORECT START POSITIONS (THEY WILL BE WRITTEN TO R0 AND RM IN PAW_ATOMS.F90)
         CALL DIMER$GET_UNMASSWEIGHTED(DIM,X1,R1P)
         CALL DIMER$GET_UNMASSWEIGHTED(DIM,X2,R2P)

         !FOR OUTPUT========================================
         !INIT THE POSITION OF THE DIMER CENTER
         RC=R1P+(R2P-R1P)/2.0D0
         RCOLD=RC
         SQDIMERDIST=DOT_PRODUCT((X1-X2),(X1-X2))
         WRITE(DPROTFIL,*)"DIMER: DIMER_PROPAGATE: WE PLACED THE DIMER AND NOW HAVE A DISTANCE OF",SQRT(SQDIMERDIST)
         !FOR OUTPUT========================================

         RETURN !DO NOT PROPAGATE!!!
      END IF

      !===============================================
      !======    KEEP THE DIMERS LENGTH   ============
      !===============================================
      IF(KDLENGTH) THEN !WE KEEP THE DISTANCE 
         D=SQRT(SQDIMERDIST)   !D IS THE CONSTAINED DISTANCE, SQD THE SQUARED CONSTRAINED
         SQD=SQDIMERDIST       !AND SQDIMERDIST THE ACTUAL SQUARED DISTANCE 

         R1P=R1
         R2P=R2


         !FOR OUTPUT========================================
         SQDIMERDIST=DOT_PRODUCT((X1-X2),(X1-X2))
         WRITE(DPROTFIL,*)"DIMER: DIMER_PROPAGATE: WE KEPT THE DIMERLENGTH AND NOW HAVE A DISTANCE OF",SQRT(SQDIMERDIST)
         !FOR OUTPUT========================================


         !THIS SHOULD BE DONE IN PAW_ATOMS
         !KDLENGTH=.FALSE. !WE DONT NEED TO PLACE THE DIMER ANY MORE
         !PLACEDIMER=.FALSE.
         RETURN !DO NOT PROPAGATE!!!
      END IF


      !===============================================
      !======        PCLIMP-TESTING       ============
      !===============================================
      IF(CLIMBPERP) THEN
         CALL DIMER_PCLIMB(R1,R2,F1,F2,R1P,R2P)
         !DO NOTHING ELSE!
         RETURN
      END IF



      !==============================================
      !======    THIS IS THE NEW EQM         ========
      !==============================================
         CALL DIMER$GET_MASSWEIGHTED(DIM,R1_,X1M)
         CALL DIMER$GET_MASSWEIGHTED(DIM,R2_,X2M)

         CALL DIMER$GET_UNMASSWEIGHTED(DIM,F1,UF1)
         CALL DIMER$GET_UNMASSWEIGHTED(DIM,F2,UF2)

         CALL DIMER$GET_MASSWEIGHTED(DIM,F1,MF1)
         CALL DIMER$GET_MASSWEIGHTED(DIM,F2,MF2)

         IF(.TRUE..OR.OPTFRICPERP.OR.OPTFRICPARA.OR.OPTFRICROT) THEN
            PERPFIRST=.FALSE.
            PARAFIRST=.FALSE.
            ROTFIRST=.FALSE.
            IF(.NOT.ALLOCATED(FPERPM)) THEN
               ALLOCATE(FPERPM(DIM))
               PERPFIRST=.TRUE.
            END IF
            IF(.NOT.ALLOCATED(FPARAM)) THEN
               ALLOCATE(FPARAM(DIM))
               PARAFIRST=.TRUE.
            END IF
            IF(.NOT.ALLOCATED(FROTM)) THEN
               ALLOCATE(FROTM(DIM))
               ROTFIRST=.TRUE.
            END IF

            !======    GET THE OPTIMAL FRICTION    ========
            CALL DIMER$OPTANNER(DIM,DT,MF1,MF2,X1,X1M,X2,X2M,FMPERP,FMPARA,FMROT,&
                 &FPERPM,FPARAM,FROTM,FPERP0,FPARA0,FROT0,FRICUSEDPERP,FRICUSEDPARA,FRICUSEDROT)


            FPERPM(:)=FPERP0(:)
            FPARAM(:)=FPARA0(:)
            FROTM(:)=FROT0(:)


            !DO WE USE OPTFRIC? 
            !USE THE OPTFRIC NOT IN THE FIRST STEP (FROTM NOT INITIALIZED)
            !USE THE USERDEFINED FIXED FRICTION INSTEAD
            IF((.NOT.OPTFRICPERP).OR.PERPFIRST) FRICUSEDPERP=FRICPERP
            IF((.NOT.OPTFRICPARA).OR.PARAFIRST) FRICUSEDPARA=FRICPARA
            IF((.NOT.OPTFRICROT).OR.ROTFIRST) FRICUSEDROT=FRICROT
         ELSE
            FRICUSEDPERP=FRICPERP
            FRICUSEDPARA=FRICPARA
            FRICUSEDROT=FRICROT            
         END IF

         !UF, BECAUSE SUBR PROP NEEDS THE FORCES IN THE FORM F=M^-(1/2)F!
         !CALL DIMER_PROP(DIM,X1P,X1,X1M,UF1,X2P,X2,X2M,UF2,DT,SQD &
         !     & ,FRICUSEDPARA,FRICUSEDPERP,FRICUSEDROT,FMPARA,FMPERP,FMROT)

         CALL DIMER_PROP06NEW(DIM,X1P,X1,X1M,UF1,X2P,X2,X2M,UF2,DT,SQD &
              & ,FRICUSEDPARA,FRICUSEDPERP,FRICUSEDROT,FMPARA,FMPERP,FMROT,G1,G2)


      CALL DIMER$GET_UNMASSWEIGHTED(DIM,X1P,R1NC)
      CALL DIMER$GET_UNMASSWEIGHTED(DIM,X2P,R2NC)

      !FOR OUTPUT========================================
      !=========================================================
      !==========  FORCES IN DIMER DIRECTION         ===========
      !=========================================================
      !CALC UNIT VECTOR IN DIMER DIRECTION

        SVAR=DOT_PRODUCT((R1-R2),(R1-R2))
        SVARV=(R1-R2)/SVAR

        !PROJECT THE FORCES ON THAT DIRECTION
        FP1=DOT_PRODUCT(SVARV,F1)
        FP2=DOT_PRODUCT(SVARV,F2)
        WRITE(DPROTFIL,*)"DIMER FORCES: I1: ",FP1," I2: ",FP2


        WRITE(DPROTFIL,*)"****************************************"
        WRITE(DPROTFIL,*)"*               DIMER                  *"
        WRITE(DPROTFIL,*)"*     ESTIMATE FOR TS COORDINATES      *"
        WRITE(DPROTFIL,*)"****************************************"
        !ALEXDEBUG        RTSOLD=RTS
        !--- ESTIMATE THE SADDLEPOINT BY THE NULL OF THE PARALLEL FORCES
        !    FOR FP2>FP1 WE HAVE NO NEGATIVE SLOPE AND THE TS PREDICTES THE
        !    MIN!
        IF(FP2*FP1.GT.0.D0) THEN 
           WRITE(DPROTFIL,*)"WARNING: POSITIVE SLOPE!"
        END IF

        RTS_=R1-SVARV*(-FP1*SVAR/(FP2-FP1))
        DO I=1,(DIM-2),3
           WRITE(DPROTFIL,FMT='(F15.5,A2,F15.5,A2,F15.5)')RTS_(I),',',RTS_(I+1),',',RTS_(I+2)
        END DO
        RTSOLD=RTSOLD-RTS_
        WRITE(DPROTFIL,*)"TS ESTIMATE FROM DIMER CHANGED ABOUT ",SQRT(DOT_PRODUCT(RTSOLD,RTSOLD))
        WRITE(DPROTFIL,*)"****************************************"



        !=========================================================
        !========  ESTIMATE TS IN UNMASSWEIGHTED SPACE   =========
        !=========================================================
        CALL DIMER$ESTIMATETS(DIM,R1,R2,F1,F2,ETOT,ETOT2)




      !=========================================================
      !==========  THE DIMERS PARALLEL MOTION       ===========
      !=========================================================
 !     RC=R1NC+(R2NC-R1NC)/2.0D0
 !     RCDIFF=SQRT(DOT_PRODUCT((R2-R1),(R2-R1))) 
 !     RCDIFF=DOT_PRODUCT((R2-R1),(RC-RCOLD))/RCDIFF
 !     PRINT*,"DIMER: DIMER_PROPAGATE: DIMER PARALLEL MOTION: ",RCDIFF
 !     RCOLD=RC
        WRITE(DPROTFIL,*)'DIMER: DIMER_PROPAGATE: DIMER PARALLEL MOTION: ',&
             &DOT_PRODUCT(R1-R2,0.5D0*(R1NC+R2NC)-0.5D0*(R1+R2))/SQRT(DOT_PRODUCT(R1-R2,R1-R2))


      !=========================================================
      !==========   THE DIMERS ACTUAL DIRECTION         =======
      !==========   AND ANGULAR CHANGE IN               =======
      !=========================================================
        !PREPARE MONITORING (HARDWIRED)
        IF(.NOT.TREADAM) THEN
           OPEN(4242,FILE="DIMERANGLE.OUT",ACTION="WRITE",POSITION='REWIND')
           DO I=1,DIM/3
              WRITE(4242,*)R1NC((I-1)*3+1)
              WRITE(4242,*)R1NC((I-1)*3+2)
              WRITE(4242,*)R1NC((I-1)*3+3)
           END DO
           DO I=1,DIM/3
              WRITE(4242,*)R2NC((I-1)*3+1)
              WRITE(4242,*)R2NC((I-1)*3+2)
              WRITE(4242,*)R2NC((I-1)*3+3)
           END DO
           CLOSE(4242)
           !OPEN THE FILE AND READ THE POSITIONS
           ALLOCATE(ANGLE1(DIM))
           ALLOCATE(ANGLE2(DIM))
           ALLOCATE(ANGLEDIR(DIM))
           INQUIRE (FILE="DIMERANGLE.IN", EXIST = TREADAMNOTTHERE)
           IF (TREADAMNOTTHERE) THEN
              OPEN(4242,FILE="DIMERANGLE.IN",ACTION="READ",POSITION='REWIND')
              DO I=1,DIM
                 READ(4242,*)ANGLE1(I)
              END DO
              DO I=1,DIM
                 READ(4242,*)ANGLE2(I)
              END DO
              
              CLOSE(4242)
              !WRITE(DPROTFIL,*)ANGLE1
              !WRITE(DPROTFIL,*)"-------------------------------------------"
              !WRITE(DPROTFIL,*)ANGLE2
              !THE NORMALIZED DIRECTION FROM GROUNDSTATE 1 TO GROUNDSTATE2
              ANGLEDIR=(ANGLE2(:)-ANGLE1(:))/SQRT(DOT_PRODUCT(ANGLE2(:)-ANGLE1(:),ANGLE2(:)-ANGLE1(:)))
           ELSE
              !CHOOSE ANY DIRECTION
              ANGLEDIR(:)=0.0D0
              ANGLEDIR(1)=1.0D0
           END IF
           DEALLOCATE(ANGLE1)
           DEALLOCATE(ANGLE2)
           TREADAM=.TRUE.
        END IF

        IF(.NOT.ALLOCATED(DIRM))ALLOCATE(DIRM(DIM))
        SVARV=(R2NC-R1NC)/SQRT(DOT_PRODUCT(R2NC-R1NC,R2NC-R1NC))
        !WRITE(DPROTFIL,*)'DIMER: DIMER_PROPAGATE: DIMER ACTUAL DIRECTION: '
        !WRITE(DPROTFIL,*)SVARV

        WRITE(DPROTFIL,*) &
       &     'DIMER: DIMER_PROPAGATE: DIMER ANGLE CHANGE IN DEGREES: ', &
       &      ACOS(DOT_PRODUCT(SVARV,DIRM))*180.D0/(4.0*ATAN(1.0D0))

       DIRM=SVARV

        !MONITORING
        WRITE(DPROTFIL,*) &
      &    'DIMER: DIMER_PROPAGATE: DIMER ANGLE MONITORING IN DEGREES : ', &
      &     ACOS(DOT_PRODUCT(SVARV,ANGLEDIR))*180.D0/(4.0*ATAN(1.0D0))


        !=========================================================
        !==========   THE DIMERS LENGTH CONTROL       ===========
        !=========================================================
        WRITE(DPROTFIL,*)"DIMER: INITIALIZED STEPFACT CHECK ",TINITSTEPFACT
        IF(.NOT.TINITSTEPFACT) THEN
           !INIT STEPFACT IF NEEDED:
           IF(NSTEPS.GT.0) THEN
              STEPFACT=(DMIN/D)**(REAL(1,KIND=8)/REAL(NSTEPS,KIND=8))
              WRITE(DPROTFIL,*)"DIMER: INITIALIZED STEPFACT TO ",STEPFACT
           END IF
           TINITSTEPFACT=.TRUE.
        END IF

        RCDIFF=SQRT(DOT_PRODUCT(R1NC+(R2NC-R1NC)/2.0D0-(R1+(R2-R1)/2.0D0),&
             &R1NC+(R2NC-R1NC)/2.0D0-(R1+(R2-R1)/2.0D0)))

        IF(DLFLEX.AND.RCDIFF.LT.RCDIFFMIN) LC=LC+1 !INCREASE COUNTER
        IF(DLFLEX.AND.RCDIFF.GT.RCDIFFMIN) LC=0 !RESET COUNTER

        IF(DLFLEX.AND.LC.GE.LCS.AND.((D-DSTEP).GT.DMIN)) THEN !SHORTEN THE LENGTH
           IF(NSTEPS.EQ.0) THEN !WE USE DSTEP AS STEPSIZE
              !SHORTEN DIMER
              SQDIMERDIST=(D-DSTEP)**2
              D=D-DSTEP
              SQD=SQDIMERDIST
              WRITE(DPROTFIL,*)"DIMER: DIMER_PROPAGATE: SHORTEND DIMER LENGTH TO ",SQRT(SQDIMERDIST)
           ELSE
              !SHORTEN DIMER %-TUAL SUCH THAT AFTER NSTEPS WE HAVE DMIN 
              IF((D-D*STEPFACT).LE.DSTEP) THEN
                 SQDIMERDIST=(D*STEPFACT)**2
                 D=D*STEPFACT
                 SQD=SQDIMERDIST
                 WRITE(DPROTFIL,*)"DIMER: DIMER_PROPAGATE: % SHORTEND DIMER LENGTH TO ",SQRT(SQDIMERDIST)
              END IF
           END IF
           !RESET COUNTER
           LC=0
        END IF

      IF(DLFLEX.AND.LC.GE.LCS.AND.((D+DSTEP).LT.DMIN)) THEN !LENGTHEN THE DIMER
         SQDIMERDIST=(D+DSTEP)**2
         D=D*DSTEP
         SQD=SQDIMERDIST 
         WRITE(DPROTFIL,*)"DIMER: DIMER_PROPAGATE: LENGTHEND DIMER LENGTH TO ",SQRT(SQDIMERDIST)
         !RESET COUNTER
         LC=0
      END IF

      R1P=R1NC
      R2P=R2NC

      V1=(R1P-R1_)/2.0D0*DT
      V2=(R2P-R2_)/2.0D0*DT
      !GET RMASS
      CALL ATOMLIST$GETR8A('MASS',0,(DIM/3),RMASS)
      !CALC CELLKIN
      DO IAT=1, DIM/3
         DO I=1,3
            DO J=1,3
               CELLKIN1(I,J)=CELLKIN1(I,J)+RMASS(IAT)*V1(IAT*3-(3-I))*V1(IAT*3-(3-J))
               CELLKIN2(I,J)=CELLKIN2(I,J)+RMASS(IAT)*V2(IAT*3-(3-I))*V2(IAT*3-(3-J))
            END DO
         END DO
      END DO


      
      RETURN
    END SUBROUTINE DIMER_PROPAGATE


!##################################################################
SUBROUTINE DIMER$ESTIMATETS(N,R1,R2,F1,F2,ETOT,ETOT2)
!##################################################################
  USE DIMER_ESTIMATETS_MODULE 
  INTEGER(4),INTENT(IN)          :: N   !THE DIMENSIONALITY
  REAL(8),INTENT(IN)             :: R1(N)
  REAL(8),INTENT(IN)             :: R2(N)
  REAL(8),INTENT(IN)             :: F1(N)
  REAL(8),INTENT(IN)             :: F2(N)
  REAL(8),INTENT(IN)             :: ETOT
  REAL(8),INTENT(IN)             :: ETOT2
  REAL(8)                        :: Z10(N),Z20(N)
  REAL(8)                        :: E0(N)
  REAL(8)                        :: F1ORTHO(N)
  REAL(8)                        :: F2ORTHO(N)
  REAL(8)                        :: EORTHO0(N)
  REAL(8)                        :: CPARA,CPERP
  REAL(8)                        :: FPARA0,FPERP0
  REAL(8)                        :: SPARA,SPERP
  REAL(8)                        :: S1,S2,A,B,C,D_
  REAL(8)                        :: DIST
  REAL(8)                        :: RVAR1,RVAR2
  REAL(8)                        :: ETS
  REAL(8)                        :: XMAX



  IF(.NOT.ALLOCATED(Z1M)) ALLOCATE(Z1M(N))
  IF(.NOT.ALLOCATED(Z2M)) ALLOCATE(Z2M(N))
  IF(.NOT.ALLOCATED(EM)) ALLOCATE(EM(N))
  IF(.NOT.ALLOCATED(EORTHOM)) ALLOCATE(EORTHOM(N))

  Z10(:)=(R1(:)+R2(:))/2.D0
  Z20(:)=R1(:)-R2(:)
  DIST=SQRT(DOT_PRODUCT(Z20(:),Z20(:)))

  E0(:)=Z20(:)/SQRT(DOT_PRODUCT(Z20(:),Z20(:)))
  
  
  !=====================================================================
  !==    DETERMINE THE FORCES (DIRECTIONS REFERENCED ON IMAGE1)       ==
  !==    WE PROJECT ON EM DIRECTIONS                                  ==
  !=====================================================================
  F1ORTHO(:)=F1(:)-EM(:)*DOT_PRODUCT(EM(:),F1(:))
  F2ORTHO(:)=F2(:)-EM(:)*DOT_PRODUCT(EM(:),F2(:))
  
  EORTHO0(:)=0.5D0*(F1ORTHO+F2ORTHO)/SQRT(DOT_PRODUCT(0.5D0*(F1ORTHO+F2ORTHO),&
       &0.5D0*(F1ORTHO+F2ORTHO)))
  
  
  FPARA0=DOT_PRODUCT(EM(:),(F1(:)+F2(:)))
  FPERP0=DOT_PRODUCT(EORTHOM,0.5D0*(F1(:)+F2(:)))






  !============================================================================
  !==    ESTIMATE ETOT IN BETWEEN IMAGE1 AND 2 USING A 3RD ORDER POLYNOMIAL ===
  !============================================================================
  !--- THE SLOPE PARALLEL TO THE DIMER
  S1=-DOT_PRODUCT(E0(:),(F1(:)))
  S2=-DOT_PRODUCT(E0(:),(F2(:)))

  D_=ETOT
  C=S1
  A=-(2.D0*ETOT2+S2*DIST+C*DIST-2.D0*D_)/(-DIST)**3
  B=(S2-C-3.D0*A*DIST**2)/(-2.D0*DIST)

  RVAR1=-B/(3.D0*A)
  IF(RVAR1**2-(C/(3.D0*A)).LT.0.D0) THEN
     PRINT*,'ERROR: VALUE IN SQRT NEGATIVE!'
     RVAR2=SQRT(RVAR1**2-(C/(3.D0*A)))
  ELSE
     RVAR2=SQRT(RVAR1**2-(C/(3.D0*A)))
  END IF

  IF(S1*S2.LT.0.D0) THEN
     !FORCES WITH OPPOSITE SIGN
     IF(RVAR1+RVAR2.LT.0.D0.AND.RVAR1+RVAR2.GT.-DIST) THEN
        XMAX=RVAR1+RVAR2
     ELSE
        XMAX=RVAR1-RVAR2
     END IF
  ELSE
     !FORCES WITH SAME SIGN
     IF(S1.GT.0D0) THEN !BOTH FORCES NEGATIVE
        IF((RVAR1+RVAR2).GT.0.D0) THEN
           XMAX=RVAR1+RVAR2
        ELSE
           XMAX=RVAR1-RVAR2
        END IF
     ELSE !BOTH FORCES POSITIVE
        IF((RVAR1+RVAR2).LT.-DIST) THEN
           XMAX=RVAR1+RVAR2
        ELSE
           XMAX=RVAR1-RVAR2
        END IF
     END IF
  END IF


  ETS=A*XMAX**3+B*XMAX**2+C*XMAX+D_
  IF(ETS.LT.ETOT.OR.ETS.LT.ETOT2) THEN
     PRINT*,'ETS IS SMALLER THAN ETOT OR ETOT2!'
  END IF

  PRINT*,XMAX
  PRINT*,S1,S2
  PRINT*,'ETOT1&2',ETOT,ETOT2
  PRINT*,'ETS ESTIMATE: ',ETS



  !=====================================================================
  !==    ESTIMATE USING DELTAF - SEEMS NOT BE BE EXACT ENOUGH!!!      ==
  !=====================================================================

  !==    DETERMINE THE DISTANCES (DIRECTIONS REFERENCED ON IMAGE1)    ==
  SPARA=DOT_PRODUCT(EM,Z10-Z1M)
  SPERP=SQRT(DOT_PRODUCT((Z10-Z1M)-SPARA*EM,(Z10-Z1M)-SPARA*EM))
  
  CPARA=-(FPARA0-FPARAM)/SPARA
  CPERP=-(FPERP0-FPERPM)/SPERP
  
  PRINT*,FPARA0,FPARAM,SPARA
  PRINT*,FPERP0,FPERPM,SPERP
  PRINT*,'TS ESTIMATION: C VALUES',CPARA,CPERP
  
  !=== ESTIMATE F=0
  PRINT*,'TSDISTPARA',FPARA0/CPARA
  PRINT*,'TSDISTPERP',FPERP0/CPERP
  
  
  
  !=== KEEP SOME VALUES FOR NEXT TURN
  !=== WE NEED TO PROJECT THEM NOW ON E0 (WHICH IS EM FOR NEXT CYCLE)
  Z1M=Z10
  Z2M=Z20
  EM=E0
  EORTHOM=EORTHO0
  !=== WE NEED TO PROJECT THEM NOW ON E0 (WHICH IS EM FOR NEXT CYCLE)
  FPARAM=DOT_PRODUCT(E0(:),(F1(:)+F2(:)))
  FPERPM=DOT_PRODUCT(EORTHO0,0.5D0*(F1(:)+F2(:)))
  



END SUBROUTINE DIMER$ESTIMATETS





!##################################################################
SUBROUTINE DIMER$STRETCH(NAT,R1,R0)
!##################################################################
  USE MPE_MODULE
!CPVERSION  USE MPE_COMM_SPACE_MODULE
  USE DIMER_MODULE
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)    :: NAT 
  REAL(8),   INTENT(INOUT) :: R1(NAT*3) !INOUT BECAUSE WE SET THE POSITION!
  REAL(8),   INTENT(IN) :: R0(NAT*3) ! ONLY FOR INITIALISATION OF RKEEP
  REAL(8)              :: SQDIMERDIST
  INTEGER(4)           :: NVAL,WORLD_ID,IAT
  LOGICAL(4)           :: LOOP
  IF(PLACEDIMER) THEN
     CALL MPE$BROADCAST('~',1,R1)!BOTH IMAGES USE THE POSITION FROM 1ST DIMERS STRC
     RETURN
  ELSE
     IF(.NOT.ALLOCATED(RKEEP)) ALLOCATE(RKEEP(3*NAT))
     IF(.NOT.ALLOCATED(R2SPLIT)) ALLOCATE(R2SPLIT(3*NAT))

     CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
     CALL DIMER$GETPROCESSLEADER2(NVAL,WORLD_ID)
     !NVAL NOW EQUALS THE WORLD_ID OF THE 1ST TASK IN DIMER2

     IF(WORLD_ID.LT.NVAL) THEN
        !ONLY THE FIRST DIMER KEEPS ITS POSITION
        IF(TFIRSTSPLIT) THEN
           RKEEP=R0
           TFIRSTSPLIT=.FALSE.
        END IF
        R1=RKEEP !KEEP R1 FIXED
     END IF

     IF(THISTASK.EQ.1) THEN
        !GIVE WORLD_ID 1 THE POSITION ON DIMER2
        IF(WORLD_ID.NE.1) THEN 
           !WE ARE 1ST TASK OF DIMER2
           CALL MPE$SEND('~',1,42,R1)
        ELSE
           !WE ARE WORLD_ID 1
           CALL MPE$RECEIVE('~',NVAL,42,R2SPLIT)
        END IF
        IF(WORLD_ID.EQ.1) THEN
           !=========APPLY THE COUPLE-CONSTRAINT=========
           IF(ASSOCIATED(ELDEST_D)) THEN
              THIS_D=>ELDEST_D
              LOOP=.TRUE.
              DO WHILE (LOOP)
                 CALL ATOMLIST$INDEX(THIS_D%ID,IAT)
                 R2SPLIT((IAT*3)-3+1)=RKEEP((IAT*3)-3+1)
                 R2SPLIT((IAT*3)-3+2)=RKEEP((IAT*3)-3+2)
                 R2SPLIT((IAT*3)-3+3)=RKEEP((IAT*3)-3+3)
                 !FOR THE FOOT CONTROLLED LOOP
                 IF(ASSOCIATED(THIS_D%NEXT_D)) THEN
                    THIS_D=>THIS_D%NEXT_D
                 ELSE
                    LOOP=.FALSE.
                 END IF
              END DO
           END IF
 

           SQDIMERDIST=DOT_PRODUCT((RKEEP-R2SPLIT),(RKEEP-R2SPLIT))
           PRINT *, "DIMER STRETCH : D SHOULD BE", STRETCHDIST," AND IS", SQRT(SQDIMERDIST)
           IF(SQDIMERDIST.GE.(STRETCHDIST*STRETCHDIST)) THEN
              STRETCH=.FALSE. !COMMUNICATE THIS TO ALL OTHER TASKS
              !SET THE NEW LENGTH IN DIMER_MODULE
              D=SQRT(SQDIMERDIST)
              SQD=SQDIMERDIST
           END IF
        END IF
     ELSE
        !ALL TASKS .NE.1 IN DIMER
     END IF
     
     CALL MPE$BROADCAST('~',1,STRETCH)
     IF(.NOT.STRETCH) THEN
        DEALLOCATE(RKEEP)
        DEALLOCATE(R2SPLIT)
     END IF
     RETURN
  END IF
END SUBROUTINE DIMER$STRETCH
   

!##################################################################
SUBROUTINE PLACE_DIMER(DIST,X1_,X2_)
!##################################################################
!PLACES THE SECOND DIMER-POINT ACCORDING TO THE LENGTH CONSTRAIN
!MOVES IT THEREFORE ALONG THE DIMER AXIS  
  USE DIMER_MODULE
  IMPLICIT NONE
  REAL(8),INTENT(IN) :: DIST 
  REAL(8),INTENT(IN) :: X1_(DIM)
  REAL(8),INTENT(INOUT) :: X2_(DIM)
  INTEGER(4)            :: IAT
  LOGICAL(4)            :: LOOP

  !======DO WE HAVE TO USE COUPLE-CONSTRAINTS ?========
  IF(ASSOCIATED(ELDEST_D)) THEN
     THIS_D=>ELDEST_D
     LOOP=.TRUE.
     DO WHILE (LOOP)
        CALL ATOMLIST$INDEX(THIS_D%ID,IAT)
        X2_((IAT*3)-3+1)=X1_((IAT*3)-3+1)
        X2_((IAT*3)-3+2)=X1_((IAT*3)-3+2)
        X2_((IAT*3)-3+3)=X1_((IAT*3)-3+3)

        !FOR THE FOOT CONTROLLED LOOP
        IF(ASSOCIATED(THIS_D%NEXT_D)) THEN
           THIS_D=>THIS_D%NEXT_D
        ELSE
           LOOP=.FALSE.
        END IF
     END DO
  END IF
    !=============PACE THE SECOND IMAGE==============
  X2_(:)=X1_(:)+DIST*(X2_(:)-X1_(:))/SQRT(DOT_PRODUCT(X2_(:)-X1_(:),X2_(:)-X1_(:)))
  RETURN
END SUBROUTINE PLACE_DIMER





!      .................................................................
       SUBROUTINE DIMER_PROP06NEW(N,X1P,X10,X1M,F1,X2P,X20,X2M,F2,DT,D2 &
      &                ,ANNEPARA,ANNEPERP,ANNEROT,MPARA,MPERP,MROT,G1_,G2_)
!      ** USES MASS WEIGHTED COORDINATES                             **
!      **    X=SQRT(M)*R      F=(-DE/DR)/SQRT(M)                     **
!         USE MODEL,ONLY:CALCMULTIPLIERITERMAX,DLAMBDA,ITERMAX
         USE DIMER_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: N
       REAL(8)   ,INTENT(IN) :: DT     ! TIME STEP
       REAL(8)   ,INTENT(IN) :: D2     ! SQUARED DIMER LENGTH
       REAL(8)   ,INTENT(OUT):: X1P(N)
       REAL(8)   ,INTENT(IN) :: X10(N)
       REAL(8)   ,INTENT(IN) :: X1M(N)
       REAL(8)   ,INTENT(IN) :: F1(N)
       REAL(8)   ,INTENT(OUT):: X2P(N)
       REAL(8)   ,INTENT(IN) :: X20(N)
       REAL(8)   ,INTENT(IN) :: X2M(N)
       REAL(8)   ,INTENT(IN) :: F2(N)
       REAL(8)   ,INTENT(IN) :: ANNEPARA   
       REAL(8)   ,INTENT(IN) :: ANNEPERP
       REAL(8)   ,INTENT(IN) :: ANNEROT
       REAL(8)   ,INTENT(IN) :: MPARA
       REAL(8)   ,INTENT(IN) :: MPERP
       REAL(8)   ,INTENT(IN) :: MROT
       REAL(8)   ,INTENT(INOUT) :: G1_(N)
       REAL(8)   ,INTENT(INOUT) :: G2_(N)
       INTEGER(4),PARAMETER  :: NITERX=1000 ! #(ITERATIONS)
       REAL(8)   ,PARAMETER  :: TOL=1.D-8
       REAL(8)   ,PARAMETER  :: DU=1.D-5
       REAL(8)               :: F1USED(N),F2USED(N)
       REAL(8)               :: FSUM(N) ! F1+F2
       REAL(8)               :: FDIFF(N) ! F1-F2
       REAL(8)               :: SVAR1,SVAR2,SVAR3,MFRAC
       REAL(8)               :: VEC1(N),VEC2(N)
       REAL(8)               :: Y1P(N),Y10(N),Y1M(N)
       REAL(8)               :: Y2P(N),Y20(N),Y2M(N),Y2BAR(N)
       REAL(8)               :: Y1DOT(N),Y2DOT(N)
       REAL(8)               :: LAMBDA
       INTEGER(4)            :: I,J
       REAL(8)               :: E(3,3) !THE POINTS NEEDED FOR THE HESSIAN
       REAL(8)               :: A(2,2)
       REAL(8)               :: B(2)
       LOGICAL(4)            :: TCONV
       REAL(8)               :: EOLD
       REAL(8)               :: G1OLD(N)
       REAL(8)               :: G2OLD(N)
       REAL(8)               :: XCENTER(N)
       REAL(8)               :: XPARA(N)
       REAL(8)               :: XPERP(N)
       REAL(8)               :: XROT(N)
!      ****************************************************************
       TCONV=.FALSE.
       F1USED(:)=F1(:)
       F2USED(:)=F2(:)

       !WEIGHT DOWN IMAGE 1 (USE A POSITIVE MPARA HERE)
       IF(WDOWN) THEN
          F1USED(:)=F1(:)*WDOWNFACT
       END IF
          

       Y10(:)=0.5D0*(X10(:)+X20(:))
       Y20(:)=X10(:)-X20(:)
       Y1M(:)=0.5D0*(X1M(:)+X2M(:))
       Y2M(:)=X1M(:)-X2M(:)
       FSUM(:)=F1USED(:)+F2USED(:)
       FDIFF(:)=F1USED(:)-F2USED(:)


       DO I=1,CALCVELOCITYITERMAX
          !===================================================================
          !=== DETERMINE Y1P
          !===================================================================
          SVAR1=1.D0/(1.D0+ANNEPERP)
          SVAR2=1.D0/(1.D0+ANNEPARA)
          
          VEC1(:)=(0.5D0*FSUM(:)+G1_(:))*SVAR1*DT**2/MPERP
          VEC2(:)=(0.5D0*FSUM(:)+G1_(:))*DT**2/(MPARA*(1.D0+ANNEPARA))
          
          Y1P(:)=(2.D0*SVAR1)*(Y10(:)-Y20(:)*DOT_PRODUCT(Y20(:),Y10(:))/D2)&
               &+(1.D0-ANNEPERP)*SVAR1*(-Y1M(:)+Y20(:)*DOT_PRODUCT(Y20(:),Y1M(:))/D2)&
               &+VEC1(:)-Y20(:)*DOT_PRODUCT(Y20(:),VEC1(:))/D2&
               &+Y20(:)*DOT_PRODUCT(Y20(:),SVAR2*(2.D0*Y10(:)-(1.D0-ANNEPARA)*Y1M(:)))/D2&
               &+Y20(:)*DOT_PRODUCT(Y20(:),VEC2(:))/D2
          


          
!PRINT*,'Y10:',Y10(:)          
!PRINT*,'Y1P:',Y1P(:)

          
          !===================================================================
          !=== DETERMINE Y2P
          !===================================================================
          !--- Y2BAR
          
          SVAR1=1.D0+ANNEROT
          
          Y2BAR(:)=Y20(:)*2.D0/SVAR1-(1.D0-ANNEROT)/SVAR1*Y2M(:)&
               &+(FDIFF(:)+G2_(:))*DT**2/(MROT*SVAR1)
          



          !---------------------------------------------------------------------
          !-----------   APPLY THE COUPLE CONSTRAINT IF NECESSARY   ------------
          !---------------------------------------------------------------------
          !DO THIS IN PAW CODE, NOT IN TESTCODE
         IF(ASSOCIATED(ELDEST_D))THEN
            CALL DIMER_APPL_COUPLE_CONSTRAINT(N,Y2BAR(:))
            !Y1BAR STAYS THE SAME BECAUSE 2*(X1+X2)/2=X1+X2!
         END IF
          


          !---------------------------------------------------------------------
          !-----------   GET RID OF UNWISHED MOTION                 ------------
          !---------------------------------------------------------------------
          IF(ONLYROT) THEN
             !SUBSTRACT THE COVERED DISTANCE OF THE CENTER OF GRAV.
             !THE COG DOES NOT MOVE:
             Y1P(:)=Y10(:)
             !ALT. TEST SUBTRACT THE PARALLEL AND THE PERP. DISTANCE
          END IF
          
          IF(INHIBITUP) THEN
             !SUBSTRACT THE PARALLEL PART OF THER COVERED DISTANCE OF THE CENTER OF GRAV.
             Y1P(:)=Y1P(:)-(Y20(:)/DOT_PRODUCT(Y20(:),Y20(:)))*DOT_PRODUCT(Y20(:),(Y1P(:)-Y10(:)))
          END IF


          
          
          !---------------------------------------------------------------------
          !-- SATISFY CONSTRAINT SVAR1*LAMBDA**2+SVAR2*LAMBDA+SVAR3-0         --
          !---------------------------------------------------------------------
          MFRAC=(2.D0*DT**2)/(MROT*(1.D0+ANNEROT))
          
          SVAR1=DOT_PRODUCT(Y2BAR(:),Y20(:))/(2.D0*MFRAC*DOT_PRODUCT(Y20(:),Y20(:)))
          SVAR2=(DOT_PRODUCT(Y2BAR(:),Y2BAR(:))-D2)/(4.D0*MFRAC**2*(DOT_PRODUCT(Y20(:),Y20(:))))
          SVAR3=SVAR1**2-SVAR2
          
          IF(SVAR3.LT.0.D0) THEN
             WRITE(DPROTFIL,*)'======LAGRANGEPARAMETER: SQRT NEGATIVE!======='
             PRINT*,'D2',D2,DOT_PRODUCT(Y2BAR(:),Y2BAR(:))-D2,DOT_PRODUCT(Y20(:),Y20(:))
             PRINT*,'DPROD',DOT_PRODUCT(Y2BAR(:),Y2BAR(:))-D2,DOT_PRODUCT(Y2BAR(:),Y2BAR(:))
             PRINT*,'MFRAC',MFRAC
             PRINT*,SVAR1,SVAR1**2,SVAR2,SVAR3
             STOP 'LAGRANGEPARAMETER NEGATIVE'

          END IF
          
          ! ATTENTION IT IS IMPORTANT WHICH ROOT IS CHOSEN. 
          ! THW WRONG ROOT ROTATES THE DIMER BY ABOUT 180 DEGREE AND MESSES UP THINGS.
          ! SWITCHING FROM ONE TO THE OTHER MESSES UP THE ITERATION
          IF(SVAR1.GT.0.D0) THEN
             LAMBDA=SVAR1-SQRT(SVAR3)
          ELSE
             LAMBDA=SVAR1+SQRT(SVAR3)
          END IF
          !SHOULD NOT CHANGE!!!Y1P(:)=Y1BAR(:)+Y1C(:)*LAMBDA
          Y2P(:)=Y2BAR(:)-2.D0*MFRAC*LAMBDA*Y20(:)
          
          
          IF(ABS(DOT_PRODUCT(Y2P,Y2P)-D2).GT.DLAMBDA) THEN
             WRITE(DPROTFIL,*)' CONSTRAIN NOT FULFILLED'
             WRITE(DPROTFIL,*)'DEVIATION ',DOT_PRODUCT(Y2P,Y2P),D2,DOT_PRODUCT(Y2P,Y2P)-D2
             PRINT*,' CONSTRAIN NOT FULFILLED'
             PRINT*,'DEVIATION ',DOT_PRODUCT(Y2P,Y2P),D2,DOT_PRODUCT(Y2P,Y2P)-D2
             STOP 'ERROR CONSTRAINT NOT FULFILLED IN PAW_DIMER.F90'
          END IF
          
          
          
          !===================================================================
          !=== DETERMINE Y1DOT, Y2DOT
          !===================================================================
          
          Y1DOT(:)=(Y1P(:)-Y1M(:))/(2.D0*DT)
          Y2DOT(:)=(Y2P(:)-Y2M(:))/(2.D0*DT)
          
          
          !===================================================================
          !=== DETERMINE NEXT GUESS FOR G1 AND G2
          !===================================================================
          G1OLD(:)=G1_(:)
          G2OLD(:)=G2_(:)

          G1_(:)=(MPERP-MPARA)/D2*(Y2DOT(:)*DOT_PRODUCT(Y20(:),Y1DOT(:))&
               &+Y20(:)*DOT_PRODUCT(Y2DOT(:),Y1DOT(:)))
          G2_(:)=4.D0*(MPARA-MPERP)*Y1DOT(:)*DOT_PRODUCT(Y20(:),Y1DOT(:))/D2

!          PRINT*,'G1DIFF',G1_-G1OLD
!          PRINT*,'G2DIFF',G2_-G2OLD
          SVAR1=0.D0
          SVAR2=0.D0
          DO J=1,N
             SVAR1=SVAR1+ABS(G1_(J)-G1OLD(J))
             SVAR2=SVAR2+ABS(G2_(J)-G2OLD(J))
          END DO
!          PRINT*,'DIFFSUM',SVAR1,SVAR2
          !IF(DOT_PRODUCT(G1_-G1OLD,G1_-G1OLD).LT.TOL&
          !     &.AND.DOT_PRODUCT(G2_-G2OLD,G2_-G2OLD).LT.TOL) EXIT
          IF(SVAR1.LT.TOL.AND.SVAR2.LT.TOL) THEN
             TCONV=.TRUE.
             EXIT
          END IF
       END DO

       !STOP 'DEBUT ROUTINE'
       !ITERATION ENDS HERE
       IF(.NOT.TCONV) THEN
          PRINT*,'U,UDOT ITRATION NOT CONVERGED',E(2,2)
          STOP 'CHECK IT!'
       END IF

       
       !-------------------------------------------
       X1P(:)=(Y1P(:)+0.5D0*Y2P(:))
       X2P(:)=(Y1P(:)-0.5D0*Y2P(:))



!!$       VEC1(:)=Y2P(:)/SQRT(DOT_PRODUCT(Y2P(:),Y2P(:)))
!!$       XCENTER(:)=Y1P(:)-Y10(:)
!!$       XPARA(:)=VEC1(:)*DOT_PRODUCT(VEC1(:),XCENTER(:))
!!$       XPERP(:)=XCENTER(:)-XPARA(:)
!!$       XROT(:)=X1P(:)-XCENTER(:)-X10(:)
!!$
!!$       PRINT*,'EKIN PARA',ABS(0.5D0*MPARA*DOT_PRODUCT(XPARA(:),XPARA(:)))/DT**2
!!$       PRINT*,'EKIN PERP',ABS(0.5D0*MPERP*DOT_PRODUCT(XPERP(:),XPERP(:)))/DT**2
!!$       PRINT*,'EKIN ROT',ABS(0.5D0*MROT*DOT_PRODUCT(XROT(:),XROT(:)))/DT**2
!!$
!!$       IF(ABS(0.5D0*MPARA*DOT_PRODUCT(XPARA(:),XPARA(:)))/DT**2.GT.0.5D0) THEN
!!$          X1P(:)=X1P(:)-(XPARA(:)-XPARA(:)/SQRT(DOT_PRODUCT(XPARA(:),XPARA(:)))*SQRT(2.D0*0.5D0/ABS(MPARA))*DT)
!!$          X2P(:)=X2P(:)-(XPARA(:)-XPARA(:)/SQRT(DOT_PRODUCT(XPARA(:),XPARA(:)))*SQRT(2.D0*0.5D0/ABS(MPARA))*DT)
!!$
!!$          VEC1(:)=Y2P(:)/SQRT(DOT_PRODUCT(Y2P(:),Y2P(:)))
!!$
!!$          Y1P(:)=(X1P(:)+X2P(:))/2.D0
!!$          XCENTER(:)=Y1P(:)-Y10(:)
!!$          XPARA(:)=VEC1(:)*DOT_PRODUCT(VEC1(:),XCENTER(:))
!!$          XPERP(:)=XCENTER(:)-XPARA(:)
!!$          XROT(:)=X1P(:)-XCENTER(:)-X10(:)
!!$          
!!$          PRINT*,'EKIN PARA',ABS(0.5D0*MPARA*DOT_PRODUCT(XPARA(:),XPARA(:)))/DT**2
!!$          PRINT*,'EKIN PERP',ABS(0.5D0*MPERP*DOT_PRODUCT(XPERP(:),XPERP(:)))/DT**2
!!$          PRINT*,'EKIN ROT',ABS(0.5D0*MROT*DOT_PRODUCT(XROT(:),XROT(:)))/DT**2
!!$       END IF

       
       
       WRITE(DPROTFIL,*)'DIMER: DIMER_PROPAGATE: DIMER MASSW. PARALLEL MOTION:',&
            &DOT_PRODUCT(Y20,0.5D0*(X1P+X2P)-0.5D0*(X10+X20))/SQRT(DOT_PRODUCT(Y20,Y20))
       



       RETURN
       
     END SUBROUTINE DIMER_PROP06NEW







!!$!      .................................................................
!!$       SUBROUTINE DIMER_PROP(N,X1P,X10,X1M,F1,X2P,X20,X2M,F2,DT,D2 &
!!$      &                ,ANNEPARA,ANNEPERP,ANNEROT,MPARA,MPERP,MROT)
!!$!      ** USES MASS WEIGHTED COORDINATES                             **
!!$!      **    X=SQRT(M)*R      F=(-DE/DR)/SQRT(M)                     **
!!$  USE DIMER_MODULE
!!$       IMPLICIT NONE
!!$       INTEGER(4),INTENT(IN) :: N
!!$       REAL(8)   ,INTENT(IN) :: DT     ! TIME STEP
!!$       REAL(8)   ,INTENT(IN) :: D2     ! SQUARED DIMER LENGTH
!!$       REAL(8)   ,INTENT(OUT):: X1P(N)
!!$       REAL(8)   ,INTENT(IN) :: X10(N)
!!$       REAL(8)   ,INTENT(IN) :: X1M(N)
!!$       REAL(8)   ,INTENT(IN) :: F1(N)
!!$       REAL(8)   ,INTENT(OUT):: X2P(N)
!!$       REAL(8)   ,INTENT(IN) :: X20(N)
!!$       REAL(8)   ,INTENT(IN) :: X2M(N)
!!$       REAL(8)   ,INTENT(IN) :: F2(N)
!!$       REAL(8)   ,INTENT(IN) :: ANNEPARA   
!!$       REAL(8)   ,INTENT(IN) :: ANNEPERP
!!$       REAL(8)   ,INTENT(IN) :: ANNEROT
!!$       REAL(8)   ,INTENT(IN) :: MPARA
!!$       REAL(8)   ,INTENT(IN) :: MPERP
!!$       REAL(8)   ,INTENT(IN) :: MROT
!!$       INTEGER(4),PARAMETER  :: NITERX=1000 ! #(ITERATIONS)
!!$       REAL(8)   ,PARAMETER  :: TOL=1.D-8
!!$       REAL(8)               :: DKIN    ! X1DOT**1-X2DOT**2
!!$       REAL(8)               :: DKINOLD  
!!$       REAL(8)               :: FSUM(N) ! F1+F2
!!$       REAL(8)               :: FDIFF(N) ! F1-F2
!!$       REAL(8)               :: FSUMPARA
!!$       REAL(8)               :: FDIFFPARA
!!$       REAL(8)               :: SVAR,SVAR1,SVAR2,SVAR3,MFRAC
!!$       REAL(8)               :: U,UDOT  !U=Y2*Y1DOT
!!$       REAL(8)               :: AP(2,2),A0(2,2),AM(2,2),AF(2,2),AC(2),AINV(2,2)
!!$       REAL(8)               :: Y1P(N),Y10(N),Y1M(N),Y1BAR(N),Y1C(N)
!!$       REAL(8)               :: Y2P(N),Y20(N),Y2M(N),Y2BAR(N),Y2C(N)
!!$       REAL(8)               :: QP,Q0,QM
!!$       REAL(8)               :: LAMBDA
!!$       REAL(8)               :: DET
!!$       REAL(8)               :: P1,G1,P2,G2
!!$       INTEGER(4)            :: ITER
!!$       REAL(8)               :: R1P(N),R2P(N),R10(N),R20(N)
!!$       REAL(8)               :: CENTER1NC(3),CENTER2NC(3)
!!$       REAL(8)               :: CC1(N),CC2(N)
!!$       REAL(8)               :: RMASS(N/3)
!!$       INTEGER(4)            :: IAT
!!$!      ****************************************************************
!!$!TO PREVENT UNASSIGNED VARIABLE BEFORE 'CHECK CONVERGENCE'
!!$       P1=42
!!$       G1=42
!!$
!!$       Y10(:)=X10(:)+X20(:)
!!$       Y20(:)=X10(:)-X20(:)
!!$       Y1M(:)=X1M(:)+X2M(:)
!!$       Y2M(:)=X1M(:)-X2M(:)
!!$       FSUM(:)=F1(:)+F2(:)
!!$       FDIFF(:)=F1(:)-F2(:)
!!$       FSUMPARA=DOT_PRODUCT(Y20(:),FSUM(:))
!!$       FDIFFPARA=DOT_PRODUCT(Y20(:),FDIFF(:))
!!$       QM=DOT_PRODUCT(Y20(:),Y1M(:))
!!$       Q0=DOT_PRODUCT(Y20(:),Y10(:))
!!$       DKIN=DOT_PRODUCT(Y10(:)-Y1M(:),Y20(:)-Y2M(:))/(4.D0*DT**2)
!!$       DKINOLD=DKIN
!!$
!!$!
!!$       DO ITER=1,CALCMULTIPLIERITERMAX
!!$!        =====================================================================
!!$!        == SELECT DKIN                                                     ==
!!$!        =====================================================================
!!$         IF(ITER.GT.2) THEN
!!$           DKIN=P1-G1*(P1-P2)/(G1-G2)
!!$         ELSE IF(ITER.EQ.2) THEN
!!$           DKIN=P1+G1
!!$         END IF
!!$!        =======
!!$!DKIN=(-1.D0+REAL(ITER,8)/REAL(NITERX,8))*2.D0
!!$!        =====================================================================
!!$!        == DETERMINE U,UDOT FOR A GIVEN DKIN                               ==
!!$!        =====================================================================
!!$         SVAR=ANNEPARA/4.D0
!!$         SVAR1=2.D0/(1.D0+SVAR)
!!$         SVAR2=(1.D0-SVAR)/(1.D0+SVAR)
!!$         SVAR3=DT**2/(MPARA*(1.D0+SVAR))
!!$         QP=SVAR1*Q0-SVAR2*QM+SVAR3*(FSUMPARA-(MPARA-MPERP)*DKIN)
!!$         U=(QP-QM)/(2.D0*DT)
!!$         UDOT=DKIN+(QP-2.D0*Q0+QM)/DT**2
!!$!WRITE(DPROTFIL,*)'U,UDOT,DKIN ',U,UDOT,DKIN
!!$!         
!!$!        =====================================================================
!!$!        == SET UP LINEAR SYSTEM OF EQUATIONS FOR Y1P,Y2P                   ==
!!$!        ==       AP*YP=A0*Y0+AM*YM+AF*F+AC*Y2*LAMBDA                       ==
!!$!        =====================================================================
!!$         MFRAC=(MPARA-MPERP)/(MPERP*D2)
!!$         AP(1,1)=1.D0+ANNEPERP/2.D0
!!$         AP(1,2)=0.5D0*MFRAC*(U*DT)
!!$         AM(1,1)=-(1.D0-ANNEPERP/2.D0)
!!$         AM(1,2)=AP(1,2)
!!$         A0(1,1)=2.D0
!!$         A0(1,2)=-MFRAC*(UDOT*DT**2)-(U*DT)/(2.D0*MPERP*D2) &
!!$        &                              *(ANNEPARA*MPARA-ANNEPERP*MPERP)
!!$         AF(1,1)=DT**2/MPERP
!!$         AF(1,2)=0.D0
!!$         AC(1)=0.D0
!!$         MFRAC=(MPARA-MPERP)/(MROT*D2)
!!$         AP(2,2)=1.D0+ANNEROT/4.D0
!!$         AP(2,1)=-0.5D0*MFRAC*(U*DT)
!!$         A0(2,1)=0.D0
!!$         A0(2,2)=2.D0
!!$         AM(2,2)=-(1.D0-ANNEROT/4.D0)
!!$         AM(2,1)=AP(2,1)
!!$         AF(2,1)=0.D0
!!$         AF(2,2)=DT**2/MROT
!!$         AC(2)=4.D0*DT**2/MROT
!!$!
!!$!        =====================================================================
!!$!        == RESOLVE EQUATIONS FOR Y1BAR,Y2BAR, Y1C Y2C                      ==
!!$!        =====================================================================
!!$         DET=AP(1,1)*AP(2,2)-AP(1,2)*AP(2,1)
!!$!WRITE(DPROTFIL,*)'AP ',AP
!!$!WRITE(DPROTFIL,*)'A0 ',A0
!!$!WRITE(DPROTFIL,*)'AM ',AM
!!$!WRITE(DPROTFIL,*)'AF ',AF
!!$!WRITE(DPROTFIL,*)'AC ',AC
!!$!STOP
!!$         AINV(1,1)= AP(2,2)/DET
!!$         AINV(1,2)=-AP(1,2)/DET
!!$         AINV(2,1)=-AP(2,1)/DET
!!$         AINV(2,2)= AP(1,1)/DET
!!$         A0(:,:)=MATMUL(AINV(:,:),A0(:,:))
!!$         AM(:,:)=MATMUL(AINV(:,:),AM(:,:))
!!$         AF(:,:)=MATMUL(AINV(:,:),AF(:,:))
!!$         AC(:)  =MATMUL(AINV(:,:),AC(:))
!!$         Y1BAR(:)=A0(1,1)*Y10(:) +A0(1,2)*Y20(:) &
!!$      &          +AM(1,1)*Y1M(:) +AM(1,2)*Y2M(:) &
!!$      &          +AF(1,1)*FSUM(:)+AF(1,2)*FDIFF(:)
!!$         Y2BAR(:)=A0(2,1)*Y10(:) +A0(2,2)*Y20(:) &
!!$      &          +AM(2,1)*Y1M(:) +AM(2,2)*Y2M(:) &
!!$      &          +AF(2,1)*FSUM(:)+AF(2,2)*FDIFF(:)
!!$         Y1C(:)=AC(1)*Y20(:)
!!$         Y2C(:)=AC(2)*Y20(:)
!!$!         
!!$

!!$!        =====================================================================
!!$!        == AP: PROJECT ROTATION & TRANSLATION OUT                              ==
!!$!        =====================================================================
!!$!        DO THIS WITH BAR VALUES -> THE LENGTH CONSTRAINED POSITIONS
!!$!        DO NOT CONTAIN CELLROTATION/TRANSLATION
!!$         X1P(:)=0.5D0*(Y1BAR(:)+Y2BAR(:))
!!$         X2P(:)=0.5D0*(Y1BAR(:)-Y2BAR(:))
!!$         CALL DIMER$GET_UNMASSWEIGHTED(DIM,X1P,R1P)
!!$         CALL DIMER$GET_UNMASSWEIGHTED(DIM,X2P,R2P)
!!$         CALL DIMER$GET_UNMASSWEIGHTED(DIM,X10,R10)
!!$         CALL DIMER$GET_UNMASSWEIGHTED(DIM,X20,R20)
!!$
!!$WRITE(DPROTFIL,*)'DEBUG43R1P',R1P
!!$WRITE(DPROTFIL,*)'DEBUG43R2P',R2P
!!$WRITE(DPROTFIL,*)'DEBUG43R10',R10
!!$WRITE(DPROTFIL,*)'DEBUG43R20',R20
!!$CALL DIMER$GET_UNMASSWEIGHTED(DIM,X1M,R10)
!!$CALL DIMER$GET_UNMASSWEIGHTED(DIM,X2M,R20)
!!$WRITE(DPROTFIL,*)'DEBUG43R1M',R10
!!$WRITE(DPROTFIL,*)'DEBUG43R2M',R20
!!$
!!$         
!!$         CALL ATOMLIST$GETR8A('MASS',0,(N/3),RMASS)
!!$         !|_ THE SAME FOR BOTH IMAGES
!!$         
!!$         IF(CENTER_ID.EQ.'COG') THEN
!!$            !DO NOTHING
!!$         ELSE
!!$            !IF WE USE NOT THE CENTER OF GRAVITY FOR ROT. CENTER
!!$            CALL ATOMLIST$INDEX(CENTER_ID,IAT)
!!$            !|_ THE SAME FOR BOTH IMAGES
!!$            
!!$            !THE CENTER FOR THE ACTUAL COORDS IS M**(1/2)*CENTER_COORD
!!$            ! THE CENTER FOR NC COORDS HAVE TO BE READ
!!$            CENTER1NC(1)=R1P((IAT*3)-3+1)
!!$            CENTER1NC(2)=R1P((IAT*3)-3+2)
!!$            CENTER1NC(3)=R1P((IAT*3)-3+3)
!!$            
!!$            CENTER2NC(1)=R2P((IAT*3)-3+1)
!!$            CENTER2NC(2)=R2P((IAT*3)-3+2)
!!$            CENTER2NC(3)=R2P((IAT*3)-3+3)
!!$         END IF
!!$         
!!$         !                                                       | WE PLACED IT THERE IN THE LAST STEP!
!!$         CALL DIMER_PROPCELLCONSTR((DIM/3),R10,R1P,RMASS,RMASS,(SQRT(RMASS(IAT))*CENTER_COORD(:)),CENTER1NC,CC1)
!!$         CALL DIMER_PROPCELLCONSTR((DIM/3),R20,R2P,RMASS,RMASS,(SQRT(RMASS(IAT))*CENTER_COORD(:)),CENTER2NC,CC2)
!!$WRITE(DPROTFIL,*)'DEBUG43CC2',CC1
!!$WRITE(DPROTFIL,*)'DEBUG43CC1',CC2
!!$
!!$         CALL DIMER$GET_MASSWEIGHTED(DIM,CC1,X1P)
!!$         CALL DIMER$GET_MASSWEIGHTED(DIM,CC2,X2P)
!!$
!!$         !RECONSTRUCT THE YBAR VECTORS (NOW WITHOUT CELLROTATION/TRANSLATION)
!!$         Y1BAR(:)=X1P(:)+X2P(:)
!!$         Y2BAR(:)=X1P(:)-X2P(:)
!!$WRITE(DPROTFIL,*)'DEBUG42Y1BAR',Y1BAR
!!$WRITE(DPROTFIL,*)'DEBUG42Y2BAR',Y2BAR
!!$WRITE(DPROTFIL,*)'DEBUG42I',ITER

!!$!        =====================================================================
!!$!        ===========   APPLY THE COUPLE CONSTRAINT IF NECESSARY   ============
!!$!        =====================================================================
!!$         IF(ASSOCIATED(ELDEST_D))THEN
!!$            CALL DIMER_APPL_COUPLE_CONSTRAINT(N,Y2BAR(:))
!!$            !Y1BAR STAYS THE SAME BECAUSE 2*(X1+X2)/2=X1+X2!
!!$         END IF
!!$
!!$!        =====================================================================
!!$!        == SATISFY CONSTRAINT SVAR1*LAMBDA**2+SVAR2*LAMBDA+SVAR3=0         ==
!!$!        =====================================================================
!!$         SVAR1=DOT_PRODUCT(Y2C(:),Y2C(:))
!!$         SVAR2=2.D0*DOT_PRODUCT(Y2BAR(:),Y2C(:))
!!$         SVAR3=DOT_PRODUCT(Y2BAR(:),Y2BAR(:))-D2
!!$         SVAR=SVAR2**2-4.D0*SVAR1*SVAR3
!!$         IF(SVAR.LT.0.D0) THEN
!!$           WRITE(DPROTFIL,*)'======LAGRANGE PARAMETER COULD NOT BE FOUND======='
!!$           WRITE(DPROTFIL,*)'SVAR ',SVAR1,SVAR2,SVAR3,SVAR
!!$           WRITE(DPROTFIL,*)'Y2BAR ',Y2BAR
!!$           WRITE(DPROTFIL,*)'Y2C   ',Y2C
!!$           STOP '======LAGRANGE PARAMETER COULD NOT BE FOUND=======ERRORSTOP'
!!$         END IF
!!$! ATTENTION IT IS IMPORTANT WHICH ROOT IS CHOSEN. 
!!$! THW WRONG ROOT ROTATES THE DIMER BY ABOUT 180 DEGREE AND MESSES UP THINGS.
!!$! SWITCHING FROM ONE TO THE OTHER MESSES UP THE ITERATION
!!$         IF(SVAR2.GT.0.D0) THEN
!!$           LAMBDA=-(SVAR2-SQRT(SVAR))/(2.D0*SVAR1)
!!$         ELSE
!!$           LAMBDA=-(SVAR2+SQRT(SVAR))/(2.D0*SVAR1)
!!$         END IF
!!$         Y1P(:)=Y1BAR(:)+Y1C(:)*LAMBDA
!!$         Y2P(:)=Y2BAR(:)+Y2C(:)*LAMBDA
!!$         IF(ABS(DOT_PRODUCT(Y2P,Y2P)-D2).GT.DLAMBDA) THEN
!!$           WRITE(DPROTFIL,*)' CONSTRAIN NOT FULFILLED'
!!$           WRITE(DPROTFIL,*)'DEVIATION ',DOT_PRODUCT(Y2P,Y2P)-D2
!!$           STOP 'ERROR CONSTRAINT NOT FULFILLED IN PAW_DIMER.F90'
!!$         END IF
!!$!         
!!$!        =====================================================================
!!$!        == CHECK CONVERGENCE                                               ==
!!$!        =====================================================================
!!$         P2=P1
!!$         G2=G1
!!$         P1=DKIN
!!$         G1=DOT_PRODUCT(Y1P(:)-Y1M(:),Y2P(:)-Y2M(:))/(4.D0*DT**2)-DKIN
!!$         WRITE(DPROTFIL,*)'DIMER:USING HARDCODED TOL!!!'
!!$         WRITE(DPROTFIL,*)'DIMER:USING HARDCODED TOL!!!'
!!$         WRITE(DPROTFIL,*)'DIMER:USING HARDCODED TOL!!!'
!!$         WRITE(DPROTFIL,*)'DIMER:USING HARDCODED TOL!!!'
!!$         WRITE(DPROTFIL,*)'DIMER:USING HARDCODED TOL!!!'
!!$         IF(ABS(G1).LT.TOL) THEN
!!$           X1P(:)=0.5D0*(Y1P(:)+Y2P(:))
!!$           X2P(:)=0.5D0*(Y1P(:)-Y2P(:))
!!$!        =====================================================================
!!$!        == DO WE USE ONLY THE ROTATION?                                    ==
!!$!        =====================================================================
!!$           IF(ONLYROT) THEN
!!$              !SUBSTRACT THE COVERED DISTANCE OF THE CENTER OF GRAV.
!!$              X1P(:)=X1P(:)-0.5D0*(Y1P(:)-Y10(:))
!!$              X2P(:)=X2P(:)-0.5D0*(Y1P(:)-Y10(:))
!!$           END IF
!!$!        =====================================================================
!!$!        == DO WE INHIBIT THE PARALLEL MOTION?                              ==
!!$!        =====================================================================
!!$           IF(INHIBITUP) THEN
!!$              !SUBSTRACT THE PARALLEL PART OF THER COVERED DISTANCE OF THE CENTER OF GRAV.
!!$              X1P(:)=X1P(:)-(Y20(:)/DOT_PRODUCT(Y20(:),Y20(:)))*DOT_PRODUCT(Y20(:),0.5D0*(Y1P(:)-Y10(:)))
!!$              X2P(:)=X2P(:)-(Y20(:)/DOT_PRODUCT(Y20(:),Y20(:)))*DOT_PRODUCT(Y20(:),0.5D0*(Y1P(:)-Y10(:)))
!!$!              X1P(:)=X1P(:)-(((X10-X20)/DOT_PRODUCT(X10-X20,X10-X20))*DOT_PRODUCT(X10-X20,X1P-X10))
!!$!              X2P(:)=X2P(:)-(((X10-X20)/DOT_PRODUCT(X10-X20,X10-X20))*DOT_PRODUCT(X10-X20,X2P-X20))
!!$           END IF
!!$
!!$           WRITE(DPROTFIL,*)'DIMER: DIMER_PROPAGATE: DIMER MASSW. PARALLEL MOTION:',&
!!$                &DOT_PRODUCT(Y20,0.5D0*(X1P+X2P)-0.5D0*(X10+X20))/SQRT(DOT_PRODUCT(Y20,Y20))
!!$          RETURN
!!$         END IF
!!$       ENDDO
!!$!STOP
!!$       STOP 'PAW_DIMER: SUBROUTINE DIMER_PROP HARDCODED TOL NOT CONVERGED'
!!$       RETURN
!!$     END SUBROUTINE DIMER_PROP
!!$


!##################################################################
SUBROUTINE DIMER$WRITEENERGYTRA(DIM_,RP)
!##################################################################
  USE DIMER_MODULE
  USE MPE_MODULE
!CPVERSION  USE MPE_COMM_SPACE_MODULE
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)      :: DIM_ 
  REAL(8),INTENT(IN)         :: RP(DIM_)
  REAL(8)                    :: ENERGY
  REAL(8)                    :: DIST    !THE DISTANCE TO THE TRANSITION STATE
  INTEGER(4)                 :: NFIL
  REAL(8)                    :: ETOT_

  IF(FIRSTTRA) THEN
     FIRSTTRA=.FALSE.
     RETURN
  END IF

  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK) 
  IF(THISTASK.NE.1) RETURN
 
  CALL LIB$SCALARPRODUCTR8(.FALSE.,DIM,1,(RP-RTS),1,(RP-RTS),DIST)
  DIST=DOT_PRODUCT((RP-RTS),(RP-RTS))
  DIST=SQRT(DIST)


  CALL MPE$QUERY('~',NTASKS,THISTASK) 
  IF(THISTASK.EQ.1) THEN !DIMER IMAGE1
     WRITE(DPROTFIL,*)"ENERGYTRA: ACTUAL DISTANCE IMAGE 1= ",DIST
  ELSE
     WRITE(DPROTFIL,*)"ENERGYTRA: ACTUAL DISTANCE IMAGE 2= ",DIST
  END IF

  !--- WRITE ONLY IN ENERGYTRA SPACING
  IF(LASTTRA+ENERGYTRA.GT.DIST) THEN
    RETURN
  ELSE
     LASTTRA=DIST
     !==================================================================
     !== GET FILE UNIT                                                ==
     !==================================================================
     CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
     IF(THISTASK.EQ.1) THEN
        CALL FILEHANDLER$UNIT('ETRA',NFIL)
        CALL ENERGYLIST$GET('TOTAL ENERGY',ENERGY)     
     ELSE
        !ALL .NEQ. TASK 1 IN DIMER
        RETURN
     END IF
     
     WRITE(NFIL,FMT='(F15.5,F15.5)')DIST,ENERGY
     FLUSH(NFIL)
     
     RETURN
  END IF

  RETURN
END SUBROUTINE DIMER$WRITEENERGYTRA




       SUBROUTINE DIMER$OPTANNER(N,DT,F1,F2,X10,X1M,X20,X2M,MPERP,MPARA,MROT,&
            &FPERPM,FPARAM,FROTM,FPERP0,FPARA0,FROT0,PERPANNER,PARAANNER,ROTANNER)
!ATTENTION: DEPENDS ON NOT MASSWEIGHTED FORCES!!!         
!      ** DETERMINES THE OPTIMAL FRICTION                    **
!      ** SUPPOSING AN HARMONIC POTENTIAL                    **
         USE DIMER_MODULE ,ONLY: DPROTFIL,EKINPARAM,EKINPERPM,EKINROTM,TPARAM,TPERPM,TROTM,&
              &TFACT,STEPS,FAUTOPARA,FAUTOPERP,FAUTOROT,FRICAUTOPARA,FRICAUTOPERP,FRICAUTOROT
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: N                                  ! 3*NAT
       REAL(8)   ,INTENT(IN) :: DT                                 ! TIMESTEP
       REAL(8)   ,INTENT(IN) :: F1(N),F2(N)                        ! THE ACTUAL FORCES
       REAL(8)   ,INTENT(IN) :: X10(N),X1M(N),X20(N),X2M(N)        ! SQUARED DIMER LENGTH
       REAL(8)   ,INTENT(IN) :: MPERP,MPARA,MROT                   ! FICT. MASS
       REAL(8)   ,INTENT(IN) :: FPERPM(N),FROTM(N),FPARAM(N)       ! THE OLD FORCES
       REAL(8)   ,INTENT(OUT):: FPERP0(N),FROT0(N),FPARA0(N)       ! THE FORCES TO BE STORED OUTSIDE
       REAL(8)   ,INTENT(OUT):: PERPANNER,PARAANNER,ROTANNER       ! OPT. FRICTION 
       REAL(8)               :: F1ORTHO(N),F2ORTHO(N) ! THE ORTHOGONAL PART OF THE FORCES
       REAL(8)               :: E0(N),EM(N) !THE UNITY VECTORS IN DIMER DIRECTION
       REAL(8)               :: Y10(N),Y20(N),Y1M(N),Y2M(N)            
       REAL(8)               :: SVAR
       REAL(8)               :: XPARA(N),XPERP(N),XCENTER(N),XROT(N)   ! THE COVERED DISTANCE
       REAL(8)               :: SPARA,SPERP,SROT   ! THE SCALAR COVERED DISTANCE

       REAL(8)               :: EPARA(N),EPERP(N),EROT(N)
       REAL(8)               :: EMWPARA(N),EMWPERP(N),EMWROT(N)
       REAL(8)               :: OMEGA2PARA,OMEGA2PERP,OMEGA2ROT
       REAL(8)               :: VEC1(N)
       REAL(8)               :: EKINPARA,EKINPERP,EKINROT
       REAL(8)               :: CPARA,CPERP,CROT
       REAL(8)               :: TPARA,TPERP,TROT
       REAL(8)               :: TFPERP,TFROT
!      ****************************************************************
       Y10(:)=(X10(:)+X20(:))/2.D0
       Y20(:)=X10(:)-X20(:)
       Y2M(:)=X1M(:)-X2M(:)
       Y1M(:)=(X1M(:)+X2M(:))/2.D0

 

       E0(:)=Y20(:)/SQRT(DOT_PRODUCT(Y20(:),Y20(:)))
       EM(:)=Y2M(:)/SQRT(DOT_PRODUCT(Y2M(:),Y2M(:)))

       F1ORTHO(:)=F1(:)-E0(:)*DOT_PRODUCT(E0(:),F1(:))
       F2ORTHO(:)=F2(:)-E0(:)*DOT_PRODUCT(E0(:),F2(:))


!      =====================================================================
!      ==    DETERMINE THE FORCES (DIRECTIONS REFERENCED ON IMAGE1)       ==
!      =====================================================================
       FPARA0(:)=E0(:)*DOT_PRODUCT(E0(:),(0.5D0*(F1(:)+F2(:))))
       FPERP0(:)=0.5D0*(F1ORTHO(:)+F2ORTHO(:))


!!$       IF(DOT_PRODUCT(F1ORTHO(:),F2ORTHO(:)).LT.0.0D0) THEN
!!$          !ANTIPARALLEL
!!$          IF(DOT_PRODUCT(F1ORTHO(:),F1ORTHO(:)).LT.DOT_PRODUCT(F2ORTHO(:),F2ORTHO(:))) THEN
!!$             FROT0(:)=F1ORTHO(:)
!!$          ELSE IF(DOT_PRODUCT(F1ORTHO(:),F1ORTHO(:)).GT.DOT_PRODUCT(F2ORTHO(:),F2ORTHO(:))) THEN
!!$             FROT0(:)=-F2ORTHO(:)!- BECAUSE WE REFERENCE ALL ON IMAGE ONE!!!
!!$          ELSE
!!$             FROT0(:)=0.0D0
!!$          END IF
!!$       ELSE
          !PARALLEL
          FROT0(:)=0.5D0*(F1ORTHO(:)-F2ORTHO(:))
!        END IF

       WRITE(DPROTFIL,*)'PERPFORCES',SQRT(DOT_PRODUCT(FPERP0(:),FPERP0(:)))
       WRITE(DPROTFIL,*)'PARAFORCES',SQRT(DOT_PRODUCT(FPARA0(:),FPARA0(:)))
       WRITE(DPROTFIL,*)'ROTFORCES',SQRT(DOT_PRODUCT(FROT0(:),FROT0(:)))



!      =====================================================================
!      ==    DETERMINE THE FORCES (DIRECTIONS REFERENCED ON IMAGE1)       ==
!      =====================================================================
       FPARA0(:)=E0(:)*DOT_PRODUCT(E0(:),(0.5D0*(F1(:)+F2(:))))
       FPERP0(:)=0.5D0*(F1ORTHO(:)+F2ORTHO(:))


       IF(DOT_PRODUCT(F1ORTHO(:),F2ORTHO(:)).LT.0.0D0) THEN
          !ANTIPARALLEL
          IF(DOT_PRODUCT(F1ORTHO(:),F1ORTHO(:)).LT.DOT_PRODUCT(F2ORTHO(:),F2ORTHO(:))) THEN
             FROT0(:)=F1ORTHO(:)
          ELSE IF(DOT_PRODUCT(F1ORTHO(:),F1ORTHO(:)).GT.DOT_PRODUCT(F2ORTHO(:),F2ORTHO(:))) THEN
             FROT0(:)=-F2ORTHO(:)!- BECAUSE WE REFERENCE ALL ON IMAGE ONE!!!
          ELSE
             FROT0(:)=F1ORTHO(:)
          END IF
       ELSE
          !PARALLEL
          FROT0(:)=0.5D0*(F1ORTHO(:)-F2ORTHO(:))
       END IF

       WRITE(DPROTFIL,*)'PERPFORCES',SQRT(DOT_PRODUCT(FPERP0(:),FPERP0(:)))
       WRITE(DPROTFIL,*)'PARAFORCES',SQRT(DOT_PRODUCT(FPARA0(:),FPARA0(:)))
       WRITE(DPROTFIL,*)'ROTFORCES',SQRT(DOT_PRODUCT(FROT0(:),FROT0(:)))

!      =====================================================================
!      ==DETERMINE THE COVERED DISTANCES (DIRECTIONS REFERENCED ON IMAGE1)==
!      =====================================================================
       XCENTER(:)=Y10(:)-Y1M(:)
       XPARA(:)=EM(:)*DOT_PRODUCT(EM(:),XCENTER(:))
       XPERP(:)=XCENTER(:)-XPARA(:)
       XROT(:)=2.D0*(X10(:)-XCENTER(:)-X1M(:))


       SPARA=SQRT(DOT_PRODUCT(XPARA(:),XPARA(:)))
       IF(SPARA.LT.1.D-10) THEN
          EPARA(:)=0.D0
       ELSE
          EPARA(:)=XPARA(:)/SPARA
       END IF

       SPERP=SQRT(DOT_PRODUCT(XPERP(:),XPERP(:)))
       IF(SPERP.LT.1.D-10) THEN
          EPERP(:)=0.D0
       ELSE
          EPERP(:)=XPERP(:)/SPERP
       END IF

       SROT=SQRT(DOT_PRODUCT(XROT(:),XROT(:)))
       IF(SROT.LT.1.D-10) THEN
          EROT(:)=0.D0
       ELSE
          EROT(:)=XROT(:)/SROT
       END IF


       WRITE(DPROTFIL,*)'PARADIST MASSWEIGHTED!',SPARA
       WRITE(DPROTFIL,*)'PERPDIST MASSWEIGHTED!',SPERP
       WRITE(DPROTFIL,*)'ROTDIST MASSWEIGHTED!',SROT


!      =====================================================================
!      ==               DETERMINE THE OPTIMAL FRICTION                    ==
!      =====================================================================

       IF(SPARA.LT.1.D-10) THEN
          CPARA=0.D0
          OMEGA2PARA=0.D0
          PARAANNER=0.D0
       ELSE
          CPARA=-(DOT_PRODUCT(FPARA0(:),EPARA(:))-DOT_PRODUCT(FPARAM(:),EPARA(:)))/&
               &SPARA
          WRITE(DPROTFIL,*)'CDIRECTPARA',CPARA
          CALL DIMER_PROPOSCILLATOR(DT,1,CPARA) !CPARA OSCILLATOR USES THE DIM 1
          WRITE(DPROTFIL,*)'CMEANPARA',CPARA
          CALL DIMER$GET_MASSWEIGHTED(N,EPARA(:),EMWPARA(:))
          OMEGA2PARA=ABS(CPARA/(MPARA*DOT_PRODUCT(EMWPARA(:),EMWPARA(:))))
          PARAANNER=0.5D0*DT*SQRT(4.D0*OMEGA2PARA)
          IF(PARAANNER.GT.1.D0)PARAANNER=1.D0
       END IF


       IF(SPERP.LT.1.D-10) THEN
          CPERP=0.D0
          OMEGA2PERP=0.D0
          PERPANNER=0.D0
       ELSE
          CPERP=-(DOT_PRODUCT(FPERP0(:),EPERP(:))-DOT_PRODUCT(FPERPM(:),EPERP(:)))/&
               &SPERP
          WRITE(DPROTFIL,*)'CDIRECTPERP',CPERP
          CALL DIMER_PROPOSCILLATOR(DT,2,CPERP) !CPERP OSCILLATOR USES THE DIM 2
          WRITE(DPROTFIL,*)'CMEANPERP',CPERP
          CALL DIMER$GET_MASSWEIGHTED(N,EPERP(:),EMWPERP(:))
          OMEGA2PERP=ABS(CPERP/(MPERP*DOT_PRODUCT(EMWPERP(:),EMWPERP(:))))
          PERPANNER=0.5D0*DT*SQRT(4.D0*OMEGA2PERP)
       END IF


       IF(SROT.LT.1.D-10) THEN
          CROT=0.D0
          OMEGA2ROT=0.D0
          ROTANNER=0.D0
       ELSE
          CROT=-(DOT_PRODUCT(FROT0(:),EROT(:))-DOT_PRODUCT(FROTM(:),EROT(:)))/&
               &SROT
          WRITE(DPROTFIL,*)'CDIRECTROT',CROT
          CALL DIMER_PROPOSCILLATOR(DT,3,CROT) !CROT OSCILLATOR USES THE DIM 3
          WRITE(DPROTFIL,*)'CMEANROT',CROT
          CALL DIMER$GET_MASSWEIGHTED(N,EROT(:),EMWROT(:))
          OMEGA2ROT=ABS(CROT/(DOT_PRODUCT(EMWROT(:),EMWROT(:))))

          OMEGA2ROT=(OMEGA2ROT-OMEGA2PARA)/MROT
          IF(OMEGA2ROT.LT.0.D0) THEN
             PRINT*,'OMEGA2ROT IS NEGATIVE - CHECK THAT'
             PRINT*,OMEGA2ROT
             PRINT*,OMEGA2PERP,OMEGA2PARA,ABS(OMEGA2PARA*2.D0*MPARA)
             OMEGA2ROT=ABS(OMEGA2ROT)
          END IF
          ROTANNER=0.5D0*DT*SQRT(4.D0*OMEGA2ROT)
       END IF
       

!      =====================================================================
!      ==             DETERMINE THE FRICTION FOR MAX EKIN                 ==
!      =====================================================================
       EKINPARA=ABS(0.5D0*MPARA*DOT_PRODUCT(XPARA(:),XPARA(:)))/DT**2
       EKINPERP=ABS(0.5D0*MPERP*DOT_PRODUCT(XPERP(:),XPERP(:)))/DT**2
       EKINROT=ABS(0.5D0*MROT*DOT_PRODUCT(XROT(:),XROT(:)))/DT**2
       TFPERP=TFACT/REAL(N-1,KIND=8)
       TFROT=TFACT/REAL(N,KIND=8)
       TPARA=EKINPARA*TFACT
       TPERP=EKINPERP*TFPERP
       TROT=EKINROT*TFROT


       PRINT*,'TEMPERATURE',TPARA,TPERP,TROT
       PRINT*,'EKIN:',EKINPARA,EKINPERP,EKINROT
       PRINT*,'MAXT:',TPARAM,TPERPM,TROTM

       IF(TPARA.GT.TPARAM) THEN
          !          PARAANNER=0.5D0*DT*(EKINPARA-EKINPARAM-DOT_PRODUCT(0.5D0*FPARA0(:),XPARA(:)))&
          !               &*DT/DOT_PRODUCT(XPARA(:),XPARA(:))
          !NEW TESTCODE:WE HAVE OPTFRIC FROM ABOVE AND ADD \DELTA E
          SVAR=ABS(0.5D0*DT*DT*(EKINPARA-TPARAM/TFACT)/(DOT_PRODUCT(XPARA(:),XPARA(:))*MPARA))
          PARAANNER=PARAANNER+SVAR
          PRINT*,'CORRECTED PARALLEL FRICTION: (VALUE/DELTA)',PARAANNER,SVAR
       END IF

       IF(TPERP.GT.TPERPM) THEN
          !PERPANNER=0.5D0*DT*(EKINPERP-EKINPERPM-DOT_PRODUCT(0.5D0*FPERP0(:),XPERP(:)))&
          !     &*DT/DOT_PRODUCT(XPERP(:),XPERP(:))
          SVAR=ABS(0.5D0*DT*DT*(EKINPERP-TPERPM/TFPERP)/(DOT_PRODUCT(XPERP(:),XPERP(:))*MPERP))
          PERPANNER=PERPANNER+SVAR
          PRINT*,'CORRECTED PERP FRICTION: (VALUE/DELTA)',PERPANNER,SVAR
       END IF

       IF(TROT.GT.TROTM) THEN
          !ROTANNER=0.5D0*DT*(EKINROT-EKINROTM-DOT_PRODUCT(0.5D0*FROT0(:),XROT(:)))&
          !     &*DT/DOT_PRODUCT(XROT(:),XROT(:))
          SVAR=ABS(0.5D0*DT*DT*(EKINROT-TROTM/TFROT)/(DOT_PRODUCT(XROT(:),XROT(:))*MROT))
          ROTANNER=ROTANNER+SVAR
          PRINT*,'CORRECTED ROT FRICTION: (VALUE/DELTA)',ROTANNER,SVAR
       END IF


!      =====================================================================
!      ==                  USE FORCE AUTOPILOT                            ==
!      =====================================================================
       !CALC 10 STEPS (DO NOT QUENCH START-UP WRIGGLES)
       IF(STEPS.GT.10) THEN
          IF(FAUTOPARA) THEN
             PRINT*,'DIRPARA',DOT_PRODUCT(XPARA,FPARA0)
             IF((DOT_PRODUCT(XPARA,FPARA0).GT.0.D0).AND.(PARAANNER.LT.FRICAUTOPARA)) THEN
                PARAANNER=FRICAUTOPARA
                PRINT*,'FORCE AUTOPILOT CORRECTED PARA FRICTION'
             END IF
          END IF

          IF(FAUTOPERP) THEN
             PRINT*,'DIRPERP',DOT_PRODUCT(XPERP,FPERP0)
             IF((DOT_PRODUCT(XPERP,FPERP0).LT.0.D0).AND.(PERPANNER.LT.FRICAUTOPERP)) THEN
                PERPANNER=FRICAUTOPERP
                PRINT*,'FORCE AUTOPILOT CORRECTED PERP FRICTION'
             END IF
          END IF

          IF(FAUTOROT) THEN
             PRINT*,'DIRROT',DOT_PRODUCT(XROT,FROT0)
             IF((DOT_PRODUCT(XROT,FROT0).LT.0.D0).AND.(ROTANNER.LT.FRICAUTOROT)) THEN
                ROTANNER=FRICAUTOROT
                PRINT*,'FORCE AUTOPILOT CORRECTED ROT FRICTION'
             END IF
          END IF
       ELSE
          STEPS=STEPS+1
       END IF




!      =====================================================================
!      ==                     ESTIMATE TS                                 ==
!      =====================================================================
       CALL DIMER$GET_UNMASSWEIGHTED(N,Y10,VEC1)
       WRITE(DPROTFIL,*)'CENTERCOORD,UNMASSWEIGTED ',Y10(:)
       WRITE(DPROTFIL,*)'CENTERCOORD,UNMASSWEIGTED ',VEC1(:)

       !TO DO: ADD HERE A CORRECTION FOR A ROTATIONAL FORCE!

       VEC1(:)=0.D0
       IF(ABS(CPARA).GT.1.D-5) VEC1(:)=VEC1(:)+FPARA0(:)/CPARA
       IF(ABS(CPERP).GT.1.D-5) VEC1(:)=VEC1(:)+FPERP0(:)/CPERP
       WRITE(DPROTFIL,*)'TS ESTIMATE',Y10(:)+VEC1(:)
       WRITE(DPROTFIL,*)'TS DIST',SQRT(DOT_PRODUCT(VEC1(:),VEC1(:)))

       RETURN
     END SUBROUTINE DIMER$OPTANNER






!      .................................................................
       SUBROUTINE DIMER$OPTANNER_OLD(N,DT,F1,F2,X10,X1M,X20,X2M,MPERP,MPARA,MROT,&
            &FPERPM,FPARAM,FROTM,FPERP0,FPARA0,FROT0,PERPANNER,PARAANNER,ROTANNER)
!      ** DETERMINES THE OPTIMAL FRICTION                    **
!      ** SUPPOSING AN HARMONIC POTENTIAL                    **
         USE DIMER_MODULE ,ONLY: DPROTFIL
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: N                                  ! 3*NAT
       REAL(8)   ,INTENT(IN) :: DT                                 ! TIMESTEP
       REAL(8)   ,INTENT(IN) :: F1(N),F2(N)                        ! THE ACTUAL FORCES
       REAL(8)   ,INTENT(IN) :: X10(N),X1M(N),X20(N),X2M(N)        ! SQUARED DIMER LENGTH
       REAL(8)   ,INTENT(IN) :: MPERP,MPARA,MROT                   ! FICT. MASS
       REAL(8)   ,INTENT(IN) :: FPERPM(N),FROTM(N),FPARAM(N)       ! THE OLD FORCES
       REAL(8)   ,INTENT(OUT):: FPERP0(N),FROT0(N),FPARA0(N)       ! THE FORCES TO BE STORED OUTSIDE
       REAL(8)   ,INTENT(OUT):: PERPANNER,PARAANNER,ROTANNER       ! OPT. FRICTION 
       REAL(8)               :: F1ORTHO(N),F2ORTHO(N) ! THE ORTHOGONAL PART OF THE FORCES
       REAL(8)               :: E0(N),EM(N) !THE UNITY VECTORS IN DIMER DIRECTION
       REAL(8)               :: Y10(N),Y20(N),Y1M(N),Y2M(N)            
       REAL(8)               :: SVAR,SVAR1 
       REAL(8)               :: XPARA(N),XPERP(N),XCENTER(N),XROT(N)   ! THE COVERED DISTANCE
!      ****************************************************************
       Y10(:)=X10(:)+X20(:)
       Y20(:)=X10(:)-X20(:)
       Y2M(:)=X1M(:)-X2M(:)
       Y1M(:)=X1M(:)+X2M(:)

       E0(:)=Y20(:)/SQRT(DOT_PRODUCT(Y20(:),Y20(:)))
       EM(:)=Y2M(:)/SQRT(DOT_PRODUCT(Y2M(:),Y2M(:)))

       F1ORTHO(:)=F1(:)-E0(:)*DOT_PRODUCT(E0(:),F1(:))
       F2ORTHO(:)=F2(:)-E0(:)*DOT_PRODUCT(E0(:),F2(:))

!      =====================================================================
!      ==    DETERMINE THE FORCES (DIRECTIONS REFERENCED ON IMAGE1)       ==
!      =====================================================================
       FPARA0(:)=E0(:)*DOT_PRODUCT(E0(:),(F1(:)+F2(:))/0.5D0)
       FPERP0(:)=0.5D0*(F1ORTHO(:)+F2ORTHO(:))

       IF(DOT_PRODUCT(F1ORTHO(:),F2ORTHO(:)).LT.0.0D0) THEN
          !ANTIPARALLEL
          IF(DOT_PRODUCT(F1ORTHO(:),F1ORTHO(:)).LT.DOT_PRODUCT(F2ORTHO(:),F2ORTHO(:))) THEN
             FROT0(:)=F1ORTHO(:)
          ELSE IF(DOT_PRODUCT(F1ORTHO(:),F1ORTHO(:)).GT.DOT_PRODUCT(F2ORTHO(:),F2ORTHO(:))) THEN
             FROT0(:)=-F2ORTHO(:)!- BECAUSE WE REFERENCE ALL ON IMAGE ONE!!!
          ELSE
             FROT0(:)=0.0D0
          END IF
       ELSE
          !PARALLEL
          FROT0(:)=0.5D0*(F1ORTHO(:)-F2ORTHO(:))
       END IF

       WRITE(DPROTFIL,*)'PERPFORCES',SQRT(DOT_PRODUCT(FPERP0(:),FPERP0(:)))
       WRITE(DPROTFIL,*)'PARAFORCES',SQRT(DOT_PRODUCT(FPARA0(:),FPARA0(:)))
       WRITE(DPROTFIL,*)'ROTFORCES',SQRT(DOT_PRODUCT(FROT0(:),FROT0(:)))

WRITE(DPROTFIL,*)"SUM42-SPLITTED FORCES",2.0D0*DOT_PRODUCT(FPERP0(:),FPERP0(:))+2.0D0*DOT_PRODUCT(FPARA0(:),FPARA0(:))+&
     &2.0D0*DOT_PRODUCT(FROT0(:),FROT0(:))
WRITE(DPROTFIL,*)"SUM42",DOT_PRODUCT(F1,F1)+DOT_PRODUCT(F2,F2)
!      =====================================================================
!      ==DETERMINE THE COVERED DISTANCES (DIRECTIONS REFERENCED ON IMAGE1)==
!      =====================================================================
       XCENTER(:)=0.5D0*Y10(:)-0.5D0*Y1M(:)
       XPARA(:)=DOT_PRODUCT(EM(:),XCENTER(:))
       XPERP(:)=XCENTER(:)-XPARA(:)
       XROT(:)=X10(:)-XCENTER(:)-X1M(:)

       WRITE(DPROTFIL,*)'PERPDIST MASSWEIGHTED!',SQRT(DOT_PRODUCT(XPERP(:),XPERP(:)))
       WRITE(DPROTFIL,*)'PARADIST MASSWEIGHTED!',SQRT(DOT_PRODUCT(XPARA(:),XPARA(:)))
       WRITE(DPROTFIL,*)'ROTDIST MASSWEIGHTED!',SQRT(DOT_PRODUCT(XROT(:),XROT(:)))

!      =====================================================================
!      ==               DETERMINE THE OPTIMAL FRICTION                    ==
!      =====================================================================

       SVAR=SQRT(DOT_PRODUCT(FPERP0-FPERPM,FPERP0-FPERPM))
       SVAR1=SQRT(DOT_PRODUCT(XPERP(:),XPERP(:)))
!ORIGINAL       PERPANNER=2.0D0*SQRT(MPERP*SVAR/SVAR1)
       PERPANNER=DT*SQRT(MPERP*SVAR/SVAR1)

WRITE(DPROTFIL,*)"DEBUG42PERP","SVAR",SVAR,"SVAR1",SVAR1,"PERPANNER",PERPANNER,"END"

       SVAR=SQRT(DOT_PRODUCT(FPARA0-FPARAM,FPARA0-FPARAM))
       SVAR1=SQRT(DOT_PRODUCT(XPARA(:),XPARA(:)))
!ORIGINAL       PARAANNER=2.0D0*SQRT(ABS(MPARA)*SVAR/SVAR1)
       PARAANNER=DT*SQRT(ABS(MPARA)*SVAR/SVAR1)

WRITE(DPROTFIL,*)"DEBUG42PARA","SVAR",SVAR,"SVAR1",SVAR1,"PERPANNER",PARAANNER,"END"
WRITE(DPROTFIL,*)SVAR/SVAR1,MPARA
       SVAR=SQRT(DOT_PRODUCT(FROT0-FROTM,FROT0-FROTM))
       SVAR1=SQRT(DOT_PRODUCT(XROT(:),XROT(:)))
!ORIGINAL       ROTANNER=2.0D0*SQRT(MROT*SVAR/SVAR1)
       ROTANNER=DT*SQRT(MROT*SVAR/SVAR1)

WRITE(DPROTFIL,*)"DEBUG42ROT","SVAR",SVAR,"SVAR1",SVAR1,"PERPANNER",ROTANNER,"END"
       WRITE(DPROTFIL,*)'PERPANNER',PERPANNER
       WRITE(DPROTFIL,*)'PARAANNER',PARAANNER
       WRITE(DPROTFIL,*)'ROTANNER',ROTANNER

       RETURN
     END SUBROUTINE DIMER$OPTANNER_OLD



!     .................................................................. 
      SUBROUTINE DIMER_APPL_COUPLE_CONSTRAINT(N,Y2)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      INTEGER(4), INTENT(IN)   :: N
      REAL(8), INTENT(INOUT)   :: Y2(N)
      INTEGER                  :: I,IAT
      LOGICAL(4)               :: LOOP
!     ******************************************************************
      WRITE(DPROTFIL,*)"APPLY COUPLE CONSTRAINTS"

      THIS_D=>ELDEST_D
      LOOP=.TRUE.
      DO WHILE (LOOP)
         CALL ATOMLIST$INDEX(THIS_D%ID,IAT)
         DO I=1,3
            WRITE(DPROTFIL,*)"SEP: ",ABS(Y2((IAT*3)-3+I))
            IF(ABS(Y2((IAT*3)-3+I)).LT.CONSTRSTEP) THEN
               !IT IS OK TO SET IT 0
               Y2((IAT*3)-3+I)=0.0D0
            ELSE
               !USE THE USERDEFINED MAX STEP
               IF(Y2((IAT*3)-3+I).LT.0.0D0) Y2((IAT*3)-3+I)=Y2((IAT*3)-3+I)+CONSTRSTEP 
               IF(Y2((IAT*3)-3+I).GT.0.0D0) Y2((IAT*3)-3+I)=Y2((IAT*3)-3+I)-CONSTRSTEP 
            END IF
            WRITE(DPROTFIL,*)"SEP NEW: ",ABS(Y2((IAT*3)-3+I))

         END DO

        !FOR THE FOOT CONTROLLED LOOP
        IF(ASSOCIATED(THIS_D%NEXT_D)) THEN
           THIS_D=>THIS_D%NEXT_D
        ELSE
           LOOP=.FALSE.
        END IF
     END DO


     RETURN
    END SUBROUTINE DIMER_APPL_COUPLE_CONSTRAINT
  












!     .................................................................. 
      SUBROUTINE DIMER$SETR8(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID_
      REAL(8)     , INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'D') THEN
        D=VAL_
      ELSE IF(ID_.EQ.'STRETCHDIST') THEN
        STRETCHDIST=VAL_
      ELSE IF(ID_.EQ.'WDOWNFACT') THEN
        WDOWNFACT=VAL_
      ELSE IF(ID_.EQ.'FMPARA') THEN
        FMPARA=VAL_
      ELSE IF(ID_.EQ.'FMPERP') THEN
        FMPERP=VAL_
      ELSE IF(ID_.EQ.'FMROT') THEN
        FMROT=VAL_
      ELSE IF(ID_.EQ.'FRICPERP') THEN
        FRICPERP=VAL_
      ELSE IF(ID_.EQ.'FRICROT') THEN
        FRICROT=VAL_
      ELSE IF(ID_.EQ.'FRICPARA') THEN
        FRICPARA=VAL_
      ELSE IF(ID_.EQ.'RCDIFFMIN') THEN
        RCDIFFMIN=VAL_
      ELSE IF(ID_.EQ.'DSTEP') THEN
        DSTEP=VAL_
      ELSE IF(ID_.EQ.'DMIN') THEN
        DMIN=VAL_
      ELSE IF(ID_.EQ.'ENERGYTRA') THEN
        ENERGYTRA=VAL_
      ELSE IF(ID_.EQ.'CONSTRSTEP') THEN
        CONSTRSTEP=VAL_
      ELSE IF(ID_.EQ.'DLAMBDA') THEN
        DLAMBDA=VAL_
      ELSE IF(ID_.EQ.'TMAXPARA') THEN
        TPARAM=VAL_
      ELSE IF(ID_.EQ.'TMAXROT') THEN
        TROTM=VAL_
      ELSE IF(ID_.EQ.'TMAXPERP') THEN
        TPERPM=VAL_
      ELSE IF(ID_.EQ.'FRICAUTOPARA') THEN
        FRICAUTOPARA=VAL_
      ELSE IF(ID_.EQ.'FRICAUTOPERP') THEN
        FRICAUTOPERP=VAL_
      ELSE IF(ID_.EQ.'FRICAUTOROT') THEN
        FRICAUTOROT=VAL_

!ALEX: THIS IS FOR PCLIMB-TESTING!
      ELSE IF(ID_.EQ.'FORCEDSTEP') THEN
        FORCEDSTEP=VAL_
!ALEX: THIS IS FOR PCLIMB-TESTING!

      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIMER$SETR8')
      END IF
      RETURN
    END SUBROUTINE DIMER$SETR8
!
!     .................................................................. 
      SUBROUTINE DIMER$GETR8(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(OUT):: VAL_
!     ******************************************************************
!!$      IF(ID_.EQ.'TIMESTEP') THEN
!!$        VAL_=DELT
!!$      ELSE IF(ID_.EQ.'AMPRE') THEN
!!$        VAL_=AMPRE
!!$      ELSE IF(ID_.EQ.'FRICTION') THEN
!!$        VAL_=ANNER
!!$      ELSE
!!$        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
!!$        CALL ERROR$CHVAL('ID_',ID_)
!!$        CALL ERROR$STOP('DIMER$GETR8')
!!$      END IF
      VAL_=0.D0 !JUST TO GET RID OF COMPILER WARNINGS
      RETURN
    END SUBROUTINE DIMER$GETR8
!     .................................................................. 
      SUBROUTINE DIMER$SETV(ID_,VAL_,N_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: N_
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(IN) :: VAL_(N_)
!     ******************************************************************
      IF(ID_.EQ.'CENTERCOORD') THEN
        CENTER_COORD(:)=VAL_(:)
      ELSE IF(ID_.EQ.'RTS') THEN
         IF(ALLOCATED(RTS)) DEALLOCATE(RTS)
         ALLOCATE(RTS(N_))
         RTS=VAL_
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIMER$SETV')
      END IF
      RETURN
    END SUBROUTINE DIMER$SETV
!
!     .................................................................. 
!!$      SUBROUTINE DIMER$GETV(ID_,VAL_,N_)
!!$!     ******************************************************************
!!$!     ******************************************************************
!!$      USE DIMER_MODULE
!!$      IMPLICIT NONE
!!$      INTEGER(4)  ,INTENT(IN) :: N_
!!$      CHARACTER(*),INTENT(IN) :: ID_
!!$      REAL(8)     ,INTENT(OUT):: VAL_(N_)
!     ******************************************************************
!!$      IF(ID_.EQ.'TIMESTEP') THEN
!!$        VAL_=DELT
!!$      ELSE IF(ID_.EQ.'AMPRE') THEN
!!$        VAL_=AMPRE
!!$      ELSE IF(ID_.EQ.'FRICTION') THEN
!!$        VAL_=ANNER
!!$      ELSE
!!$        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
!!$        CALL ERROR$CHVAL('ID_',ID_)
!!$        CALL ERROR$STOP('DIMER$GETR8')
!!$      END IF
!!$      RETURN
!!$    END SUBROUTINE DIMER$GETV
!     .................................................................. 
      SUBROUTINE DIMER$SETCH(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID_
      CHARACTER(32), INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'CENTER_ID') THEN
        CENTER_ID=VAL_
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIMER$SETCH')
      END IF
      RETURN
    END SUBROUTINE DIMER$SETCH
!
!     .................................................................. 
      SUBROUTINE DIMER$GETCH(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*),  INTENT(IN) :: ID_
      CHARACTER(32) ,INTENT(OUT):: VAL_
!     ******************************************************************
!!$      IF(ID_.EQ.'TIMESTEP') THEN
!!$        VAL_=DELT
!!$      ELSE IF(ID_.EQ.'AMPRE') THEN
!!$        VAL_=AMPRE
!!$      ELSE IF(ID_.EQ.'FRICTION') THEN
!!$        VAL_=ANNER
!!$      ELSE
!!$        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
!!$        CALL ERROR$CHVAL('ID_',ID_)
!!$        CALL ERROR$STOP('DIMER$GETR8')
!!$      END IF
      VAL_=' ' !
      RETURN
    END SUBROUTINE DIMER$GETCH
!
!     .................................................................. 
      SUBROUTINE DIMER$SETL4(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL     ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'DIMER') THEN
        DIMER=VAL_
      ELSE IF(ID_.EQ.'KDLENGTH') THEN
        KDLENGTH=VAL_
      ELSE IF(ID_.EQ.'PLACEDIMER') THEN
        PLACEDIMER=VAL_
      ELSE IF(ID_.EQ.'STRETCH') THEN
        STRETCH=VAL_
      ELSE IF(ID_.EQ.'DIMERFOLLOWDOWN') THEN
        DIMERFOLLOWDOWN=VAL_
      ELSE IF(ID_.EQ.'INHIBITUP') THEN
        INHIBITUP=VAL_
      ELSE IF(ID_.EQ.'INHIBITPERP') THEN
        INHIBITPERP=VAL_
      ELSE IF(ID_.EQ.'ONLYROT') THEN
        ONLYROT=VAL_
      ELSE IF(ID_.EQ.'ONLYPERP') THEN
        ONLYPERP=VAL_
      ELSE IF(ID_.EQ.'WDOWN') THEN
        WDOWN=VAL_
      ELSE IF(ID_.EQ.'OPTFRICPARA') THEN
        OPTFRICPARA=VAL_
      ELSE IF(ID_.EQ.'OPTFRICPERP') THEN
        OPTFRICPERP=VAL_
      ELSE IF(ID_.EQ.'OPTFRICROT') THEN
        OPTFRICROT=VAL_
      ELSE IF(ID_.EQ.'DLFLEX') THEN
        DLFLEX=VAL_
      ELSE IF(ID_.EQ.'FAUTOPARA') THEN
        FAUTOPARA=VAL_
      ELSE IF(ID_.EQ.'FAUTOPERP') THEN
        FAUTOPERP=VAL_
      ELSE IF(ID_.EQ.'FAUTOROT') THEN
        FAUTOROT=VAL_

!ALEX: THIS IS FOR PCLIMB-TESTING!
      ELSE IF(ID_.EQ.'CLIMBPERP') THEN
       CLIMBPERP =VAL_
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIMER$SETL4')
      END IF
      RETURN
    END SUBROUTINE DIMER$SETL4
!
!     .................................................................. 
      SUBROUTINE DIMER$GETL4(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL     ,INTENT(OUT):: VAL_
!     ******************************************************************
      IF(ID_.EQ.'DIMER') THEN
        VAL_=DIMER
      ELSE IF(ID_.EQ.'PLACEDIMER') THEN
        VAL_=PLACEDIMER
      ELSE IF(ID_.EQ.'STRETCH') THEN
        VAL_=STRETCH
      ELSE IF(ID_.EQ.'KDLENGTH') THEN
        VAL_=KDLENGTH
      ELSE IF(ID_.EQ.'DIMERFOLLOWDOWN') THEN
        VAL_=DIMERFOLLOWDOWN
!      ELSE IF(ID_.EQ.'MOVE') THEN
!        VAL_=TDYN
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIMER$SETL4')
      END IF
      RETURN
    END SUBROUTINE DIMER$GETL4
!     .................................................................. 
      SUBROUTINE DIMER$SETI4(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'CLACMULTIPLIERITERMAX') THEN
        CALCMULTIPLIERITERMAX=VAL_
      ELSE IF(ID_.EQ.'NSTEPS') THEN
        NSTEPS=VAL_
      ELSE IF(ID_.EQ.'LCS') THEN
        LCS=VAL_
!      ELSE IF(ID_.EQ.'STRETCH') THEN
!        STRETCH=VAL_
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIMER$SETL4')
      END IF
      RETURN
    END SUBROUTINE DIMER$SETI4
!
!     .................................................................. 
      SUBROUTINE DIMER$GETI4(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(OUT):: VAL_
!     ******************************************************************
!!$      IF(ID_.EQ.'STOP') THEN
!!$        VAL_=TSTOP
!!$      ELSE IF(ID_.EQ.'RANDOMIZE') THEN
!!$        VAL_=TRANDOMIZE
!!$      ELSE IF(ID_.EQ.'MOVE') THEN
!!$        VAL_=TDYN
!!$      ELSE
!!$        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
!!$        CALL ERROR$CHVAL('ID_',ID_)
!!$        CALL ERROR$STOP('DIMER$SETL4')
!!$      END IF
      VAL_=0 !JUST TO GET RID OF COMPILER WARNINGS
      RETURN
    END SUBROUTINE DIMER$GETI4
!


!##################################################################
SUBROUTINE DIMER$GETPROCESSLEADER2(NVAL,WORLD_ID)
!##################################################################
  USE MPE_MODULE
  USE DIMER_MODULE, ONLY:THISTASK,NTASKS
  IMPLICIT NONE
  INTEGER(4),INTENT(OUT) :: NVAL !THE WORLD_ID OF PROCESSLEADER ON IMAGE2
  INTEGER(4),INTENT(OUT) :: WORLD_ID !THE WORLD_ID OF THE CALLING PROCESS
  INTEGER(4)             :: NWORLDTASKS

         NVAL=0
         CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
         CALL MPE$QUERY('~',NWORLDTASKS,WORLD_ID)
         
         IF(WORLD_ID.NE.1.AND.THISTASK.EQ.1) THEN
            !WE ARE TASK 1 IN DIMER2
            NVAL=WORLD_ID
         END IF
         
         CALL MPE$COMBINE('~','+',NVAL)
         !NVAL NOW EQUALS THE WORLD_ID OF THE 1ST TASK IN DIMER2
         !ON ALL MACHINES [THE RESULT IS KNOWN ON ALL MACHINES]
  RETURN
END SUBROUTINE DIMER$GETPROCESSLEADER2







!##################################################################
SUBROUTINE DIMER$GET_MASSWEIGHTED(N,R,X)
!##################################################################
  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: N !THE DIMENSIONALITY
  REAL(8), INTENT(IN)   :: R(N) ! THE COORDINATE VECTOR IN REAL SPACE
  REAL(8), INTENT(OUT)  :: X(N) ! THE MASSWEIGHTED COORDINATE VECTOR  
  REAL(8)               :: RMASS(N/3) ! THE MASSVECTOR
  INTEGER(4)            :: I
      !GET RMASS
      CALL ATOMLIST$GETR8A('MASS',0,(N/3),RMASS)
      DO I=1,N/3
         X((I-1)*3+1)=SQRT(RMASS(I))*R((I-1)*3+1)
         X((I-1)*3+2)=SQRT(RMASS(I))*R((I-1)*3+2)
         X((I-1)*3+3)=SQRT(RMASS(I))*R((I-1)*3+3)
      END DO
      RETURN
    END SUBROUTINE DIMER$GET_MASSWEIGHTED


!##################################################################
SUBROUTINE DIMER$GET_UNMASSWEIGHTED(N,X,R)
!##################################################################
  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: N !THE DIMENSIONALITY
  REAL(8), INTENT(OUT)   :: R(N) ! THE COORDINATE VECTOR IN REAL SPACE
  REAL(8), INTENT(IN)  :: X(N) ! THE MASSWEIGHTED COORDINATE VECTOR  
  REAL(8)               :: RMASS(N/3) ! THE MASSVECTOR
  INTEGER(4)            :: I
      !GET RMASS
      CALL ATOMLIST$GETR8A('MASS',0,(N/3),RMASS)
      DO I=1,N/3
         R((I-1)*3+1)=(1.0D0/SQRT(RMASS(I)))*X((I-1)*3+1)
         R((I-1)*3+2)=(1.0D0/SQRT(RMASS(I)))*X((I-1)*3+2)
         R((I-1)*3+3)=(1.0D0/SQRT(RMASS(I)))*X((I-1)*3+3)
      END DO
      RETURN
    END SUBROUTINE DIMER$GET_UNMASSWEIGHTED


!     ..................................................................
      SUBROUTINE DIMER$INIT()
!     ******************************************************************
!     **  INITIALIZES THE DIMER DEFAULT VALUES                        **
!     **  SHOULD BE CALLED BEFORE DIMER_READIN@PAW_IOROUTINES         **
!     ******************************************************************
      USE DIMER_MODULE
      USE DIMER_OSCILLATOR_MODULE
      IMPLICIT NONE
      !********************************
        DIMERFOLLOWDOWN       =  .FALSE.
        PLACEDIMER            =  .TRUE.
        D                     =   1.5D0   ! DIMERDISTANCE
        CALCMULTIPLIERITERMAX = 1000       ! EMERGENCY EXIT FOR THE ITERATION LOOP AFTER X STEPS
        DLAMBDA               =   1.0D-10 ! EXACTNESS OF LAMBDA (= LAGRANGE MULTIPLIER) CALCULATION
        CALCVELOCITYITERMAX   = 1000       ! EMERGENCY EXIT FOR THE ITERATION LOOP AFTER X STEPS
!!        DVELOCITY             =   1.0D-10 ! EXACTNESS OF VELOCITY ITERATION 
        LCS                   =   1       !HOW MANY STEPS WITH NO CH.OF DIMERCENT. BEF. SHORTEN SQD
        RCDIFFMIN             =   0.0001  !DIFFERENCE IN DIMER CENTER POSITION BELOW LC=LC+1 
        DSTEP                 =   0.0050   !D=D-DSTEP
        NSTEPS                =   0       !IF 0 WE USE DSTEP
        DMIN                  =   0.4D0   !THE MINIMAL DIMER LENGTH (NO FURTHER REDUCEMENT OF D)
        DLFLEX                =  .FALSE.  !FLEXIBLE DIMER LENGTH
        KDLENGTH              =  .FALSE.   !KEEP THE LENGTH OF THE STARTUP
        INHIBITUP             =  .FALSE.  !INHIBIT THE UPWARD MOTION OF THE DIMER
        INHIBITPERP           =  .FALSE.
        ONLYPERP              =  .FALSE.
        ONLYROT               =  .FALSE.
        WDOWN                 =  .FALSE.  !WEIGHT DOWNDIMER IMAGE 1 (F1=WDOWNFACT*F1)
        WDOWNFACT             =  10.0D0        
        ENERGYTRA             =   0.5       ! WRITE ENERGY/DISTANCE TO RTS ALL N STEPS 0 -> NEVER
        CENTER_ID             = 'COG'  ! CENTER OF GRAVITY IS DEFAULT
!!        DROT                  =   0.001  
        STRETCH               = .FALSE.
        STRETCHDIST           =   0.1D0
        FMPARA                =  -1.0D0
        FMPERP                =   1.0D0
        FMROT                 =   1.0D0
        FRICPARA              =   0.1D0
        FRICPERP              =   0.1D0
        FRICROT               =   0.1D0
        OPTFRICPARA           = .FALSE.
        OPTFRICPERP           = .FALSE.
        OPTFRICROT            = .FALSE.  
        CONSTRSTEP            = 0.01D0
        TPARAM                = 500.D0
        TPERPM                = 500.D0
        TROTM                 = 500.D0
        
        FAUTOPARA             =.FALSE.
        FAUTOPERP             =.FALSE.
        FAUTOROT              =.FALSE.
        FRICAUTOPARA          = 0.05D0
        FRICAUTOPERP          = 0.05D0
        FRICAUTOROT           = 0.05D0

!ALEX: THIS IS FOR PCLIMB-TESTING!
        CLIMBPERP             = .FALSE.
        FORCEDSTEP            = 0.001D0
        TFIRSTCLIMBSTEP       = .TRUE.
!ALEX: THIS IS FOR PCLIMB-TESTING!


        TREADAM=.FALSE.
        TREADAMNOTTHERE=.FALSE.


        !CONNECT THE PROTOCOL FILES ETC.
        CALL DIMER_INIT_FILES()


        !INIT THE OSCILLATOR ARRAY
        ALLOCATE(OSCM(ODIM))
        ALLOCATE(OSC0(ODIM))
        ALLOCATE(OSCP(ODIM))
        ALLOCATE(OSCMASS(ODIM))  
        ALLOCATE(OSCANNER(ODIM))
        ALLOCATE(OSCC(ODIM))

        OSCP(:)=0.D0
        OSCM(:)=0.D0
        OSC0(:)=0.D0
        
        !ACTUALLY SET BY HAND
        OSCMASS(:)=0.1D0*1822.88D0   !100. A.U.
        OSCC(:)=1.D0 
        OSCANNER(:)=0.1D0 !WE USE OPT FRICTION IN THE SUBROUTINE (WE DO NOT KNOW DT NOW) 

        !THE CODE FOR OPT. FRICTION:
        !DO I=1,ODIM
           !WE USE OPTIMIZED FRICTION HERE
        !   OSCANNER(I)=DT**2*SQRT(OSCC(I)/OSCMASS(I))
        !END DO

        RETURN
      END SUBROUTINE DIMER$INIT



!     ..................................................................
      SUBROUTINE DIMER$INITDIM()
!     ******************************************************************
!     **  INITIALIZES THE DIMER DIM VARIABLE AND THE DATA STRUCTURES  **
!     **  DEPENDING THEREOF. SHOULD BE CALLED AFTER STRCIN            **
!     ******************************************************************
        USE DIMER_MODULE
        IMPLICIT NONE
        INTEGER(4)                  :: NAT_
        !********************************

        CALL ATOMLIST$NATOM(NAT_)
        DIM=3*NAT_
        ALLOCATE(RTS(DIM))
        ALLOCATE(G1(DIM))
        ALLOCATE(G2(DIM))
        
        G1(:)=0.D0
        G2(:)=0.D0

        RETURN
      END SUBROUTINE DIMER$INITDIM


!     ..................................................................
      SUBROUTINE DIMER$REPORT_SETTINGS(NFIL)
!     ******************************************************************
!     **  REPORTS THE SETTINGS TO STDOUT                              **
!     **  SHOULD BE CALLED AFTER READIN_DIMER@PAW_IOROUTINES.F90      **
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)       :: NFIL  
      LOGICAL(4)                  :: LOOP
      !********************************
      
      WRITE(NFIL,FMT="(A)")"========================================================"
      WRITE(NFIL,FMT="(A)")"=                                                      ="
      WRITE(NFIL,FMT="(A)")"=            SETTINGS FROM !CONTROL!DIMER              ="
      WRITE(NFIL,FMT="(A)")"=                                                      ="
      WRITE(NFIL,FMT="(A)")"========================================================"
      WRITE(NFIL,*)"DIMER ",DIMER
      WRITE(NFIL,*)"DIMERFOLLOWDOWN",DIMERFOLLOWDOWN      
      WRITE(NFIL,*)"D ",D                      
      WRITE(NFIL,*) "CALCMULTIPLIERITERMAX ",CALCMULTIPLIERITERMAX  
      WRITE(NFIL,*) "DLAMBDA ",DLAMBDA                
      WRITE(NFIL,*) "CALCVELOCITYITERMAX ",CALCVELOCITYITERMAX    
      WRITE(NFIL,*) "DVELOCITY ",DVELOCITY              
      WRITE(NFIL,*) "LCS ",LCS                     
      WRITE(NFIL,*) "RCDIFFMIN ",RCDIFFMIN              
      WRITE(NFIL,*) "DSTEP ",DSTEP                   
      WRITE(NFIL,*) "NSTEPS ",NSTEPS                   
      WRITE(NFIL,*) "DMIN ",DMIN                   
      WRITE(NFIL,*) "DLFLEX ",DLFLEX                  
      WRITE(NFIL,*) "KDLENGTH ",KDLENGTH               
      WRITE(NFIL,*) "INHIBITUP ",INHIBITUP              
      WRITE(NFIL,*) "INHIBITPERP ",INHIBITPERP            
      WRITE(NFIL,*) "ONLYPERP ",ONLYPERP               
      WRITE(NFIL,*) "ONLYROT ",ONLYROT                 
      WRITE(NFIL,*) "WDOWN ",WDOWN                  
      WRITE(NFIL,*) "WDOWNFACT ",WDOWNFACT                     
      WRITE(NFIL,*) "ENERGYTRA ",ENERGYTRA
      WRITE(NFIL,*)"NATOMS ",DIM/3
      WRITE(NFIL,*)"STRETCH",STRETCH
      WRITE(NFIL,*)"STRETCHDIST",STRETCHDIST
      WRITE(NFIL,*)"FMPARA ",FMPARA
      WRITE(NFIL,*)"FMPERP ",FMPERP
      WRITE(NFIL,*)"FMROT ",FMROT
      WRITE(NFIL,*)"FRICPARA ",FRICPARA
      WRITE(NFIL,*)"FRICPERP ",FRICPERP
      WRITE(NFIL,*)"FRICROT ",FRICROT
      WRITE(NFIL,*)"OPTFRICPARA ",OPTFRICPARA
      WRITE(NFIL,*)"OPTFRICPERP ",OPTFRICPERP
      WRITE(NFIL,*)"OPTFRICROT ",OPTFRICROT
      WRITE(NFIL,*)"FAUTOPARA", FAUTOPARA   
      WRITE(NFIL,*)"FAUTOPERP",  FAUTOPERP  
      WRITE(NFIL,*)"FAUTOROT",    FAUTOROT
      WRITE(NFIL,*)"FRICAUTOPARA",FRICAUTOPARA 
      WRITE(NFIL,*)"FRICAUTOPERP",FRICAUTOPERP
      WRITE(NFIL,*)"FRICAUTOROT",FRICAUTOROT

      WRITE(NFIL,*)"TMAXPARA ",TPARAM
      WRITE(NFIL,*)"TMAXPERP ",TPERPM
      WRITE(NFIL,*)"TMAXROT ",TROTM
      WRITE(NFIL,*)"CONSTRSTEP ",CONSTRSTEP
      !REPORT THE CONSTRAINTS
      IF(ASSOCIATED(ELDEST_D)) THEN
         THIS_D=>ELDEST_D
         LOOP=.TRUE.
         DO WHILE (LOOP)
            WRITE(NFIL,*)"COUPLED ATOM=",THIS_D%ID
            IF(ASSOCIATED(THIS_D%NEXT_D)) THEN
               THIS_D=>THIS_D%NEXT_D
            ELSE
               LOOP=.FALSE.
            END IF
         END DO
      ELSE
         WRITE(NFIL,*)'NO CONSTRAINTS'
      END IF
      
      
      !ALEX: THIS IS FOR PCLIMB-TESTING!
      WRITE(NFIL,*)"CLIMBPERP ",CLIMBPERP
      WRITE(NFIL,*)"FORCEDSTEP ",FORCEDSTEP
      WRITE(NFIL,*)"TFIRSTCLIMBSTEP ",TFIRSTCLIMBSTEP
      !ALEX: THIS IS FOR PCLIMB-TESTING!

      WRITE(NFIL,*)"========================================================"
    END SUBROUTINE DIMER$REPORT_SETTINGS
    

!!$MODULE DIMER_CONSTR_MODULE
!!$  IMPLICIT NONE
!!$  TYPE DIMER_CONSTR_TYPE
!!$     CHARACTER(32)                    :: ID        !KEYWORD
!!$     TYPE(DIMER_CONSTR_TYPE),POINTER  :: NEXT      !YOUNGER BROTHER
!!$     TYPE(DIMER_CONSTR_TYPE),POINTER  :: PREV      !ELDER BROTHER
!!$     !   TYPE(DIMER_CONSTR_TYPE),POINTER  :: ELDEST    !ELDEST BROTHER
!!$ 
!!$  END TYPE DIMER_CONSTR_TYPE
!!$
!!$  TYPE(DIMER_CONSTR_TYPE),POINTER  :: ELDEST
!!$  TYPE(DIMER_CONSTR_TYPE),POINTER  :: THIS
!!$  !LOGICAL(4)                       :: CSINI
!!$  !*******************************************
!!$  PUBLIC DIMER$CONSTRLIST_INIT
!!$  PUBLIC DIMER$CONSTRLIST_ADD
!!$
!!$CONTAINS
!!$
  SUBROUTINE DIMER$CONSTRLIST_INIT(ID)
    USE DIMER_MODULE
    IMPLICIT NONE
    CHARACTER(32),INTENT(IN)         :: ID
    !*******************************************
    ALLOCATE(ELDEST_D)
    THIS_D=>ELDEST_D
    NULLIFY(THIS_D%NEXT_D)
    NULLIFY(THIS_D%PREV_D)
    THIS_D%ID=ID
    RETURN
  END SUBROUTINE DIMER$CONSTRLIST_INIT



  SUBROUTINE DIMER$CONSTRLIST_ADD(ID)
    USE DIMER_MODULE
    IMPLICIT NONE
    CHARACTER(32),INTENT(IN)        :: ID
    !******************************************
   IF(.NOT.ASSOCIATED(ELDEST_D)) THEN
      CALL DIMER$CONSTRLIST_INIT(ID)
   ELSE
      !GO TO THE TOP
      THIS_D=>ELDEST_D
      IF(THIS_D%ID.EQ.ID) THEN
         CALL ERROR$MSG('CAN NOT USE THE ID FOR DIMER CONSTRAINTS TWICE')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('DIMER$CONSTRLIST_ADD')  
      END IF
      
      !DESCEND TO THE YOUNGEST BROTHER
      !THEREBY TEST IF THE ID IS ALREADY IN USE 
      DO WHILE (ASSOCIATED(THIS_D%NEXT_D))
         IF(THIS_D%ID.EQ.ID) THEN
            CALL ERROR$MSG('CAN NOT USE THE ID FOR DIMER CONSTRAINTS TWICE')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('DIMER$CONSTRLIST_ADD')  
         END IF
         THIS_D=>THIS_D%NEXT_D
      END DO
      
      ALLOCATE(THIS_D%NEXT_D)
      NULLIFY(THIS_D%NEXT_D%NEXT_D)
      THIS_D%NEXT_D%PREV_D=>THIS_D
      THIS_D%NEXT_D%ID=ID
   END IF
   RETURN
 END SUBROUTINE DIMER$CONSTRLIST_ADD
!!$ 
!!$END MODULE DIMER_CONSTR_MODULE


      !===============================================
      !======        PCLIMB-TESTING       ============
      !===============================================


  SUBROUTINE DIMER_PCLIMB(R1,R2,F1,F2,R1P,R2P)
    USE DIMER_MODULE
    IMPLICIT NONE
      REAL(8)   ,INTENT(IN)    :: R1(DIM)               !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,INTENT(IN)    :: R2(DIM)               !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,INTENT(IN)    :: F1(DIM)               !FORCE FROM POTENTIAL
      REAL(8)   ,INTENT(IN)    :: F2(DIM)               !FORCE FROM POTENTIAL

      REAL(8)   ,INTENT(OUT)   :: R1P(DIM)              !NEW POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,INTENT(OUT)   :: R2P(DIM)              !NEW POSITION VECTORS OF DIMERPOINT 1 & 2

      REAL(8)                  :: SVARV(DIM)
      REAL(8)                  :: FSUM(DIM) !F1+F2
    !******************************************
    FSUM=F1+F2
    IF(.NOT.ALLOCATED(PCLIMBDIR)) ALLOCATE(PCLIMBDIR(DIM))

    !***** INITIALIZE THE DIRECTION ********
    IF(TFIRSTCLIMBSTEP) THEN
       !=== GET THE NORMALIZED PARALLEL DIRECTION ===
       SVARV=(R1(:)-R2(:))/SQRT(DOT_PRODUCT(R1-R2,R1-R2))

       !DETERMINE THE DIRECTION
       PCLIMBDIR(:)=FSUM(:)-SVARV(:)*DOT_PRODUCT(SVARV,FSUM)

       !NORMALIZE THE DIRECTION
       PCLIMBDIR(:)=PCLIMBDIR(:)/SQRT(DOT_PRODUCT(PCLIMBDIR,PCLIMBDIR))

       !WE INITIALIZED THE DIRECTION
       TFIRSTCLIMBSTEP=.FALSE.
       PRINT *, "PCLIMB : INITIALIZED THE FORCEDIRECTION"
       PRINT *, "PCLIMB : THE CLIMBFORCE IS",DOT_PRODUCT(PCLIMBDIR,FSUM)
       R1P(:)=R1
       R2P(:)=R2
       RETURN
    END IF


    !************* CLIMB UP : -FORCEDSTEP*PCLIMBDIR *******
    R1P(:)=R1(:)-FORCEDSTEP*PCLIMBDIR(:)
    R2P(:)=R2(:)-FORCEDSTEP*PCLIMBDIR(:)
    PRINT *, "PCLIMB : THE CLIMBFORCE IS",DOT_PRODUCT(PCLIMBDIR,FSUM)
   RETURN
 END SUBROUTINE DIMER_PCLIMB







!====================================================================
! OLD CODE TO PREVENT ROTATION/TRANSLATION ACROSS THE IMAGES
!====================================================================
!##################################################################
SUBROUTINE DIMER$CELLCONSTRAINT(NAT,RCC)
!##################################################################
  USE DIMER_MODULE
  USE MPE_MODULE
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)  :: NAT
  REAL(8)   ,INTENT(OUT) :: RCC(NAT*3)

  INTEGER(4)             :: IAT
  REAL(8),ALLOCATABLE    :: RMASS(:)
  REAL(8),ALLOCATABLE    :: R(:),R2(:)
  INTEGER(4)             :: NVAL,WORLD_ID,NWORLDTASKS
  REAL(8)                :: RCENTER(3)


!     ==================================================================
!     ==   DIMER STRUCTURE-CONSTRAINT FOR TRANSLATION AND ROTATION    ==
!     ==================================================================

         !========= CORRECT TRANSLATION ===========


         ALLOCATE(RMASS(NAT))
         ALLOCATE(R(3*NAT))

         CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)         

         IF(CENTER_ID.EQ.'COG') THEN
            !USE THE CENTER OF GRAVITY FOR ROTATION
            CALL ATOMLIST$GETR8A('MASS',0,NAT,RMASS)
            !CORRECT THE POSITIONS
            CALL PLACECENTER(NAT,R,RMASS,RCC)

PRINT *,"CENTEROFGRAVITY"
         ELSE
PRINT *,"GIVENATOM"
            !USE A GIVEN ATOM FOR ROTATION
            CALL ATOMLIST$INDEX(CENTER_ID,IAT)
            CALL ATOMLIST$GETR8A('R(0)',IAT,3,RCENTER)
            CALL PLACEATOM(.TRUE.,NAT,R,RCENTER,RCC)
         END IF

         !SET THE ACTUAL POSITIONS
         R=RCC
PRINT *,"DRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"
PRINT *,R
PRINT *,"DRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"
!DO THIS OUTSIDE
!!$         !SET THE CORRECTED POSITIONS
!!$         CALL ATOMLIST$SETR8A('R(0)',0,3*NAT,RCC)
!!$         CALL ATOMLIST$SETR8A('R(-)',0,3*NAT,RCC)
!!$         !NOTE: IF YOU USE THIS CODE ELSEWHERE YOU WILL
!!$         !HAVE TO CARE ABOUT R(-) !!!


         !========= CORRECT ROTATION  ===========
         !COMMUNICATE TO DIMER1
         ALLOCATE(R2(3*NAT))

         NVAL=0
         CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
         CALL MPE$QUERY('~',NWORLDTASKS,WORLD_ID)

         IF(WORLD_ID.NE.1.AND.THISTASK.EQ.1) THEN
            !WE ARE TASK 1 IN DIMER2
            NVAL=WORLD_ID
         END IF
         
         CALL MPE$COMBINE('~','+',NVAL)
         !NVAL NOW EQUALS THE WORLD_ID OF THE 1ST TASK IN DIMER2

         IF(WORLD_ID.EQ.NVAL) THEN 
            !WE ARE 1ST TASK OF DIMER2
            CALL MPE$SEND('~',1,42,R)
         ELSE IF(WORLD_ID.EQ.1) THEN
            !WE ARE WORLD_ID 1
            CALL MPE$RECEIVE('~',NVAL,42,R2)
         END IF
         !WE GOT BOTH VECTORS ON THE 1ST TASK

         IF(WORLD_ID.EQ.1) THEN
            CALL PLACEROT(NAT,R,R2,DROT,RCC)            
         END IF



         !COMMUNICATE THE ROT. VECTOR RCC BACK TO DIMER2  
         IF(WORLD_ID.EQ.1) THEN 
            !WE ARE WORLD_ID 1
            CALL MPE$SEND('~',NVAL,43,RCC)
         ELSE IF(WORLD_ID.EQ.NVAL) THEN !1ST IN 2ND DIMER R
            CALL MPE$RECEIVE('~',1,43,RCC)
         END IF
         
         !ALL TASKS BACK
         IF(WORLD_ID.GE.NVAL) THEN 


            !ONLY ON 2ND DIMER
            CALL MPE$BROADCAST('MONOMER',1,RCC)
            R=RCC
!!$            !SET THE CORRECTED POSITIONS
!!$            CALL ATOMLIST$SETR8A('R(0)',0,3*NAT,RCC)
!!$            CALL ATOMLIST$SETR8A('R(-)',0,3*NAT,RCC)
         END IF
         !NOTE: IF YOU USE THIS CODE ELSEWHERE YOU WILL
         !HAVE TO CARE ABOUT R(-) !!!


         !========= SET THE CENTER BACK  ===========
         CALL PLACEATOM(.FALSE.,NAT,R,CENTER_COORD,RCC)
!!$         CALL ATOMLIST$SETR8A('R(0)',0,3*NAT,RCC)
!!$         CALL ATOMLIST$SETR8A('R(-)',0,3*NAT,RCC)
PRINT *,"D2RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"
PRINT *,RCC
PRINT *,"D2RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"

         DEALLOCATE(RMASS)
         DEALLOCATE(R)
         DEALLOCATE(R2)  
  
  RETURN
END SUBROUTINE DIMER$CELLCONSTRAINT




!##################################################################
SUBROUTINE DIMER_PROPCELLCONSTR(NAT,R,R2,RMASS,RMASS2,RCENTER,RCENTER2,RCC)
!##################################################################
! TRANSLATES AND ROTATES R2 IN SUCH A WAY (EITER COG OR GIVEN ATOM AS CENTER)
! THAT THE DISTANCE \SUM{(R-R2)^2} IS MINIMIZED
! THE CENTER OF THE ROTATED VECTOR IS MOVED TO CENTER_COORD FROM .CNTL FILE
!
! THIS IS USED TO PROJECT TRANSLATION / ROTATION OUT OF THE *NOT CONSTRAINED*
! COORDS (R1NC,R2NC) BEFORE THE DIMER LENGTH CONSTRAINT IS APPLIED
! (R=R1, R2=R1NC; R=R2, R2=R2NC )


  USE DIMER_MODULE
  USE MPE_MODULE
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)  :: NAT
  REAL(8)   ,INTENT(IN)  :: R(NAT*3)
  REAL(8)   ,INTENT(IN)  :: R2(NAT*3)
  REAL(8)   ,INTENT(IN)  :: RMASS(NAT)
  REAL(8)   ,INTENT(IN)  :: RMASS2(NAT)
  REAL(8)   ,INTENT(IN)  :: RCENTER(3)
  REAL(8)   ,INTENT(IN)  :: RCENTER2(3)
  REAL(8)   ,INTENT(OUT) :: RCC(NAT*3)

  REAL(8)                :: RCC1(3*NAT),RCC2(3*NAT)





!     ==================================================================
!     ==   DIMER STRUCTURE-CONSTRAINT FOR TRANSLATION AND ROTATION    ==
!     ==================================================================

         !========= CORRECT TRANSLATION ===========
         ! SET BOTH CENTER TO THE ORIGIN
 
         IF(CENTER_ID.EQ.'COG') THEN
            !USE THE CENTER OF GRAVITY FOR ROTATION
            !CORRECT THE POSITIONS
            CALL PLACECENTER(NAT,R,RMASS,RCC1)
            CALL PLACECENTER(NAT,R2,RMASS,RCC2)
         ELSE
            !USE A GIVEN CENTER FOR ROTATION
            CALL PLACEATOM(.TRUE.,NAT,R,RCENTER,RCC1)
            CALL PLACEATOM(.TRUE.,NAT,R2,RCENTER2,RCC2)
         END IF


         !========= CORRECT ROTATION  ===========
            CALL PLACEROT(NAT,RCC1,RCC2,DROT,RCC)            

         !========= SET THE CENTER BACK  ===========
         !PLACE ONLY THE SECOND TO THE WISHED CENTER_COORD
PRINT *, "EEEEEEEEEEEEEEEEE"
PRINT *, RCC
PRINT *, "EEEEEEEEEEEEEEEEE"
         RCC2=RCC
         CALL PLACEATOM(.FALSE.,NAT,RCC2,CENTER_COORD,RCC)
  
  RETURN
END SUBROUTINE DIMER_PROPCELLCONSTR


!##################################################################
SUBROUTINE PLACECENTER(NAT,R,RMASS,RCC)
!##################################################################
  USE DIMER_MODULE
  IMPLICIT NONE
  INTEGER(4), INTENT(IN) :: NAT
  REAL(8),INTENT(IN)     :: R(3*NAT)
  REAL(8),INTENT(IN)     :: RMASS(NAT)
  REAL(8),INTENT(OUT)    :: RCC(3*NAT)

  REAL(8)                :: VAL(3,NAT)
  REAL(8)                :: RCOM(3)
  REAL(8)                :: SUMM
  INTEGER                :: I

  RCOM(:)=0.0D0
  SUMM=0.0D0

  VAL(:,:)=RESHAPE(R,(/3,NAT/))
  
  DO I=1, NAT
     RCOM(:)=RCOM(:)+RMASS(I)*VAL(:,I)
     SUMM=SUMM+RMASS(I)
  END DO
  
  RCOM=RCOM/SUMM

  DO I=1, NAT
     VAL(:,I)=VAL(:,I)-RCOM(:)
  END DO

  RCC=RESHAPE(VAL,(/3*NAT/))  
  
  RETURN
END SUBROUTINE PLACECENTER


!##################################################################
SUBROUTINE PLACEATOM(TO_ZERO,NAT,R,RCENTER,RCC)
!##################################################################
  USE DIMER_MODULE
  IMPLICIT NONE
  LOGICAL(4), INTENT(IN) :: TO_ZERO !T->RCENTER TO ZERO, F-> ZERO -RCENTER 
  INTEGER(4), INTENT(IN) :: NAT
  REAL(8),INTENT(IN)     :: R(3*NAT)
  REAL(8),INTENT(IN)     :: RCENTER(3)
  REAL(8),INTENT(OUT)    :: RCC(3*NAT)
  REAL(8)                :: VAL(3,NAT)
  INTEGER                :: I

  VAL(:,:)=RESHAPE(R,(/3,NAT/))

  IF(TO_ZERO) THEN
     DO I=1, NAT
        VAL(:,I)=VAL(:,I)-RCENTER
     END DO
  ELSE
     DO I=1, NAT
        VAL(:,I)=VAL(:,I)+RCENTER
     END DO
  END IF
  
  RCC=RESHAPE(VAL,(/3*NAT/))  
  RETURN
END SUBROUTINE PLACEATOM

!##################################################################
SUBROUTINE PLACEROT(NAT,R1,R2,ROTDIF,R2C)
!##################################################################
  USE DIMER_MODULE
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)  :: NAT
  REAL(8),INTENT(IN)     :: R1(3*NAT)
  REAL(8),INTENT(IN)     :: R2(3*NAT)
  REAL(8),INTENT(IN)     :: ROTDIF
  REAL(8),INTENT(OUT)    :: R2C(3*NAT) !ONLY DIMER2 WILL BE ROTATED
  REAL(8)     :: VAL1(3,NAT),VAL2(3,NAT),VAL2ROT(3)
  REAL(8)     :: A0
  REAL(8)     :: AACT
  REAL(8)     :: AMAX
  REAL(8)     :: ASTEP
  REAL(8)     :: B0
  REAL(8)     :: BACT
  REAL(8)     :: BMAX
  REAL(8)     :: BSTEP
  REAL(8)     :: C0
  REAL(8)     :: CACT
  REAL(8)     :: CMAX
  REAL(8)     :: CSTEP
  REAL(8)     :: ABEST
  REAL(8)     :: BBEST
  REAL(8)     :: CBEST
  REAL(8)     :: ROTMINV(3)
  REAL(8)     :: ROTMIN
  REAL(8)     :: ROTMINBEST
  REAL(8)     :: U(3,3)
  REAL(8)     :: DIFACT
  REAL(8)     :: ROTMINOLD,ROTMINSUM
  LOGICAL(4)  :: START,FIRST
  INTEGER(4)  :: I,J,K,L
  REAL(8)     :: PI
  
  START=.TRUE.
  FIRST=.TRUE.
  ROTMINBEST=1.0D42 !SOMETHING HIGH
  ROTMINOLD=0.0D0
  PI=4.0D0*ATAN(1.0D0)
  DIFACT=ROTDIF+42.D0
!WRITE(DPROTFIL,*)'ROTMINBESTDEBUG',ROTMINBEST


 DO WHILE(START.OR.(DIFACT.GT.ROTDIF))
!PRINT *,"WE ARE STILL IN ROTATIONS",DIFACT,ROTDIF
    IF(START) THEN
       A0=0.0D0
       AMAX=2.0D0*PI
       ASTEP=(2.D0*PI)/40.0D0
       B0=0.0D0
       BMAX=PI
       BSTEP=(2.D0*PI)/40.0D0
!       C0=0.0D0
!       CMAX=PI
!       CSTEP=PI/20.0D0
       START=.FALSE.
    ELSE
       A0=ABEST-ASTEP
       AMAX=ABEST+ASTEP
       ASTEP=ASTEP/10.0D0

       B0=BBEST-BSTEP
       BMAX=BBEST+BSTEP
       BSTEP=BSTEP/10.0D0

!       C0=CBEST-CSTEP
!       CMAX=CBEST+CSTEP
!       CSTEP=CSTEP/10.0D0
    END IF
!WRITE(DPROTFIL,*)'ROTMINBESTDEBUG',ROTMINBEST    
!WRITE(DPROTFIL,*)'R1DEBUG',R1
!WRITE(DPROTFIL,*)'R2',R2
    DO J=0,INT((AMAX-A0)/ASTEP)
           DO K=0,INT((BMAX-B0)/BSTEP)
!              DO L=0,((CMAX-C0)/CSTEP)

                 AACT=A0+J*ASTEP
                 BACT=B0+K*BSTEP
!                 CACT=C0+L*CSTEP
                 
                 VAL1(:,:)=RESHAPE(R1,(/3,NAT/))
                 VAL2(:,:)=RESHAPE(R2,(/3,NAT/))
                 
              U(1,1)=COS(AACT)*COS(BACT)
              U(1,2)=-SIN(AACT)
              U(1,3)=COS(AACT)*SIN(BACT)
              U(2,1)=SIN(AACT)*COS(BACT)
              U(2,2)=COS(AACT)
              U(2,3)=SIN(AACT)*SIN(BACT)
              U(3,1)=-SIN(BACT)
              U(3,2)=0.0D0
              U(3,3)=COS(BACT)
                 


!!$                 U(1,1)=COS(AACT)*COS(BACT)*COS(CACT)
!!$                 U(1,2)=-COS(AACT)*COS(BACT)*SIN(CACT)-SIN(AACT)*COS(CACT)
!!$                 U(1,3)=COS(AACT)*SIN(BACT)
!!$                 U(2,1)=SIN(AACT)*COS(BACT)*COS(CACT)+COS(AACT)*SIN(CACT)
!!$                 U(2,2)=-SIN(AACT)*COS(BACT)*SIN(CACT)+COS(AACT)*COS(CACT)
!!$                 U(2,3)=SIN(AACT)*SIN(BACT)
!!$                 U(3,1)=-SIN(BACT)*COS(CACT)
!!$                 U(3,2)=SIN(BACT)*SIN(CACT)
!!$                 U(3,3)=COS(BACT)
!WRITE(DPROTFIL,*)'ROTMINBESTDEBUG3',ROTMINBEST                 
                 ROTMINSUM=0.0D0
                 DO I=1,NAT
                    CALL LIB$MATMULR8(3,3,1,U,VAL2(:,I),VAL2ROT)
                    ROTMINV=(VAL1(:,I)-VAL2ROT)
 !                   CALL LIB$SCALARPRODUCTR8(.FALSE.,3,1,ROTMINV,1,ROTMINV,ROTMIN)
                    ROTMIN=DOT_PRODUCT(ROTMINV(:),ROTMINV(:))
                   ROTMINSUM=ROTMINSUM+ROTMIN
                 END DO
!WRITE(DPROTFIL,*)'ROTMINBESTDEBUG4',ROTMINBEST                 
!WRITE(DPROTFIL,*)'ROTMINBEST,ROTMINSUM',ROTMINBEST,ROTMINSUM
                 IF(ROTMINSUM.LT.ROTMINBEST) THEN
                    ROTMINBEST=ROTMINSUM
                    ABEST=AACT
                    BBEST=BACT
!WRITE(DPROTFIL,*)'ABEST,AACT,BBEST,BACT',ABEST,AACT,BBEST,BACT
!WRITE(DPROTFIL,*)'ABEST,AACT,BBEST,BACT',BBEST,BACT
                 END IF
                 DIFACT=ABS(ROTMINOLD-ROTMINSUM)
!WRITE(DPROTFIL,*)'DIFACT',DIFACT,ROTMINSUM                 
                 ROTMINOLD=ROTMINSUM
!              END DO
              END DO
        END DO
!PRINT *,"ROTMINBEST1,ABEST,BBEST",ROTMINBEST,ABEST,BBEST



                 !WRITE THE CORRECT VECTOR
!!$                 U(1,1)=COS(ABEST)*COS(BBEST)*COS(CBEST)
!!$                 U(1,2)=-COS(ABEST)*COS(BBEST)*SIN(CBEST)-SIN(ABEST)*COS(CBEST)
!!$                 U(1,3)=COS(ABEST)*SIN(BBEST)
!!$                 U(2,1)=SIN(ABEST)*COS(BBEST)*COS(CBEST)+COS(ABEST)*SIN(CBEST)
!!$                 U(2,2)=-SIN(ABEST)*COS(BBEST)*SIN(CBEST)+COS(ABEST)*COS(CBEST)
!!$                 U(2,3)=SIN(ABEST)*SIN(BBEST)
!!$                 U(3,1)=-SIN(BBEST)*COS(CBEST)
!!$                 U(3,2)=SIN(BBEST)*SIN(CBEST)
!!$                 U(3,3)=COS(BBEST)

                 U(1,1)=COS(ABEST)*COS(BBEST)
                 U(1,2)=-SIN(ABEST)
                 U(1,3)=COS(ABEST)*SIN(BBEST)
                 U(2,1)=SIN(ABEST)*COS(BBEST)
                 U(2,2)=COS(ABEST)
                 U(2,3)=SIN(ABEST)*SIN(BBEST)
                 U(3,1)=-SIN(BBEST)
                 U(3,2)=0.0D0
                 U(3,3)=COS(BBEST)

                 DO I=1,NAT
                    CALL LIB$MATMULR8(3,3,1,U,VAL2(:,I),VAL2ROT)
                    DO J=1,3
                       R2C(I*3-3+J)=VAL2ROT(J)
                    END DO
                 END DO
     END DO


     RETURN
END SUBROUTINE PLACEROT

!!$!THIS IS HARDWIRED AND ONLY FOR MY USE
!!$SUBROUTINE DIMER$REREAD()
!!$  USE DIMER_MODULE
!!$  IMPLICIT NONE
!!$  LOGICAL(4)                    :: EX
!!$  INTEGER(4)                    :: NFIL,EOS
!!$  CHARACTER(256)                :: STRING
!!$  LOGICAL(4)                    :: LVAR
!!$  REAL(8)                       :: RVAR
!!$
!!$  NFIL=4242
!!$  INQUIRE (FILE="CASE.CNTLNEW", EXIST = EX)
!!$  IF ( EX ) THEN
!!$     OPEN (NFIL, FILE="CASE.CNTLNEW", STATUS="OLD", &
!!$                  &ACCESS="SEQUENTIAL", ACTION="READ", POSITION="REWIND")
!!$     DO 
!!$        READ (NFIL, FMT="(A)" , ADVANCE="YES" ,IOSTAT=EOS) STRING 
!!$        IF(EOS.EQ.-1) EXIT
!!$        
!!$        IF(TRIM(ADJUSTL(STRING)).EQ.'OPTFRICPARA') THEN
!!$           READ (NFIL, FMT="(L1)" , ADVANCE="YES" ,IOSTAT= EOS) LVAR 
!!$           CALL DIMER$SETL4('OPTFRICPARA',LVAR)
!!$           WRITE(DPROTFIL,*)'OPTFRICPARA SET TO',LVAR
!!$        ELSE IF(TRIM(ADJUSTL(STRING)).EQ.'FRICPARA') THEN
!!$           READ (NFIL, FMT="(F)" , ADVANCE="YES" ,IOSTAT= EOS) RVAR
!!$           CALL DIMER$SETR8('FRICPARA',RVAR)
!!$           WRITE(DPROTFIL,*)'FRICPARA SET TO',RVAR
!!$        ELSE
!!$           WRITE(DPROTFIL,*)'REREAD CONTROLFILE: IDENTIFIER NOT RECOGNIZED:',TRIM(ADJUSTL(STRING))
!!$           !STOP
!!$        END IF
!!$     END DO
!!$     CLOSE(NFIL)
!!$     CALL SYSTEM('RM -RF CASE.CNTLNEW')
!!$  END IF
!!$
!!$  RETURN
!!$END SUBROUTINE DIMER$REREAD


SUBROUTINE DIMER_PROPOSCILLATOR(DT,I,FORCE)
  USE DIMER_OSCILLATOR_MODULE
  REAL(8),INTENT(IN)                 :: DT
  INTEGER(4),INTENT(IN)              :: I       !WHICH ARRAY
  REAL(8),INTENT(INOUT)              :: FORCE   !EXTRA FORCE
  REAL(8)                            :: SVAR1,SVAR2

  !DESRIPTION: WE USE THE INCOMING VALUE AS EXTRA FORCE OF A DRIVEN OSCILLATOR.
  !THE RETURN FORCE VALUE IS THE VALUE OF THE ACTUAL FORCE OF THE OSCILLATOR 
  !AFTER PROPAGATION.
  !THIS IS USEFUL, BECAUSE WE SUPPOSE TO HAVE AN INCOMING VALUE WHICH OSCILLATES 
  !AROUND A MEAN VALUE. WHAT WE NEED IS THE MEAN VALUE. THE ACTUAL FORCE OF THE 
  !DRIVEN OSCILLATOR (AT STEADY STATE) WILL GIVE US THIS MEAN VALUE.

  IF(I.GT.ODIM) THEN
     STOP 'DIMER OSCILLATOR ARRAY OUT OF RANGE'
  END IF

  !WE USE OPTIMIZED FRICTION HERE
  OSCANNER(I)=DT**2*SQRT(OSCC(I)/OSCMASS(I))


  SVAR1=(1.D0-OSCANNER(I))
  SVAR2=1.D0/(1.D0+OSCANNER(I))

  !PROPAGATE
  OSCP(I)=(SVAR2*DT**2/OSCMASS(I))*FORCE&
       &-(OSCC(I)*DT**2/OSCMASS(I)-2.D0)*SVAR2*OSC0(I)&
       &-SVAR1*SVAR2*OSCM(I)


  !SWITCH
  OSCM(I)=OSC0(I)
  OSC0(I)=OSCP(I)

  !USE -1 BECAUSE WE SEARCH THE STEADY STATE
  FORCE=OSCC(I)*OSC0(I)


  RETURN
END SUBROUTINE DIMER_PROPOSCILLATOR
