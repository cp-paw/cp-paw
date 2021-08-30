!........1.........2.........3.........4.........5.........6.........7.........8
MODULE ATOMS_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: ATOMS                                                              **
!**                                                                           **
!**  PURPOSE: MAINTAINS THE ATOMIC POSITIONS                                  **
!**                                                                           **
!**  FUNCTIONS:                                                               **
!**    ATOMS$SETR8(ID,VAL)                                                    **
!**    ATOMS$REPORT                                                           **
!**    ATOMS$PROPAGATE                                                        **
!**    ATOMS$INITIALIZE                                                       **
!**                                                                           **
!**  REMARKS:                                                                 **
!**    FOR HISTORICAL REASONS, PART OF THIS OBJECT IS NAMED ATOMLIST          **
!**                                                                           **
!**  STATEOFTHIS RUNS THROUGH THE FOLLOWING LIFE CYCLE:                       **
!**    'NOT INITIALIZED'                                                      **
!**    'INITIALIZED'                                                          **
!**    'PROPAGATED WITHOUT CONSTRAINTS'                                       **
!**    'PROPAGATED WITH CONSTRAINTS'                                          **
!**    'SWITCHED'                                                             **
!**    NEXT: 'PROPAGATED WITHOUT CONSTRAINTS'                                 **
!**                                                                           **
!******************PETER E. BLOECHL, IBM RESEARCH LABORATORY (1996)*************
REAL(8)          :: DELT=0.D0      ! TIME STEP
REAL(8)          :: AMPRE=0.D0     ! TARGET TEMPERATURE FOR RANDOMIZATION
REAL(8)          :: ANNER=0.D0     ! FRICTION
LOGICAL(4)       :: TSTOP=.FALSE.  ! ZERO INITIAL VELOCITIES
LOGICAL(4)       :: TDYN=.TRUE.    ! ATOMIC MOTION IS SWITCHED OFF
LOGICAL(4)       :: TRANDOMIZE=.FALSE.  ! RANDOMIZE INITIAL VELOCITIES
LOGICAL(4)       :: START=.FALSE.  ! IGNORE ATOMIC STRUCTURE ON RESTART FILE
REAL(8)          :: ANNEE=0.D0
LOGICAL(4)       :: TCONSTRAINTREFERENCE=.FALSE.
LOGICAL(4)       :: TNONEGATIVEFRICTION=.FALSE. ! REMOVES NEGATIVE FRICTIONS
                                                ! DO NOT USE WITH THERMOSTATS
LOGICAL(4)       :: TOPTFRIC=.FALSE. ! USE OPTIMUM FRICTION
                                     ! OVERWRITES FRICTION OF THERMOSTAT AND AUTOPILOT
INTEGER(4)       :: NAT=0
CHARACTER(32),ALLOCATABLE :: NAME(:)
REAL(8)      ,ALLOCATABLE :: R0(:,:)
REAL(8)      ,ALLOCATABLE :: RM(:,:)
REAL(8)      ,ALLOCATABLE :: RP(:,:)
REAL(8)      ,ALLOCATABLE :: FORCE(:,:)
REAL(8)      ,ALLOCATABLE :: RMASS(:)
REAL(8)      ,ALLOCATABLE :: PSG2(:)
REAL(8)      ,ALLOCATABLE :: PSG4(:)
INTEGER(4)   ,ALLOCATABLE :: ISPECIES(:)
REAL(8)      ,ALLOCATABLE :: Z(:)
REAL(8)      ,ALLOCATABLE :: ZV(:)
REAL(8)      ,ALLOCATABLE :: CHARGE(:)
REAL(8)      ,ALLOCATABLE :: CHARGEANDMOMENTS(:,:)  ! CHARGE AND MAGNETIC MOMENTS INTEGRATED
                                                  ! OVER ASS SPHERE (JUST STORING FOR TRAJECTORY)
CHARACTER(128)            :: STATEOFTHIS='NOT INITIALIZED'
! THE FOLLOWING IS REQUIRED TO ESTIMATE THE OPTIMUM FRICTIONS ON THE ATOMS
REAL(8)                   :: AOPT1AV       ! AOPT=SQRT(AOPT1AV/AOPT2AV)
REAL(8)                   :: AOPT2AV
REAL(8)                   :: MIXAOPT=0.3D0 ! MIXAOPT=1: NO FLOATING AVERAGE
REAL(8)      ,ALLOCATABLE :: ANNERVEC0(:)
REAL(8)      ,ALLOCATABLE :: ANNERVECM(:)
REAL(8)      ,ALLOCATABLE :: R2M(:,:)
END MODULE ATOMS_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMS$SETR8(ID_,VAL_)
!     **************************************************************************
!     **************************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(IN) :: VAL_
!     **************************************************************************
      IF(ID_.EQ.'TIMESTEP') THEN
        DELT=VAL_
      ELSE IF(ID_.EQ.'AMPRE') THEN
        AMPRE=VAL_
      ELSE IF(ID_.EQ.'FRICTION') THEN
        ANNER=VAL_
      ELSE IF(ID_.EQ.'ANNEE') THEN
        ANNEE=VAL_
      ELSE IF(ID_.EQ.'MIXAOPT') THEN
        MIXAOPT=VAL_
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMS$SETR8')
      END IF
      RETURN
      END      
!
!     .................................................................. 
      SUBROUTINE ATOMS$GETR8(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(OUT):: VAL_
!     ******************************************************************
      IF(ID_.EQ.'TIMESTEP') THEN
        VAL_=DELT
      ELSE IF(ID_.EQ.'AMPRE') THEN
        VAL_=AMPRE
      ELSE IF(ID_.EQ.'FRICTION') THEN
        VAL_=ANNER
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMS$GETR8')
      END IF
      RETURN
      END      
!
!     .................................................................. 
      SUBROUTINE ATOMS$SETL4(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL     ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'STOP') THEN
        TSTOP=VAL_
      ELSE IF(ID_.EQ.'RANDOMIZE') THEN
        TRANDOMIZE=VAL_
      ELSE IF(ID_.EQ.'MOVE') THEN
        TDYN=VAL_
      ELSE IF(ID_.EQ.'START') THEN
        START=VAL_
      ELSE IF(ID_.EQ.'USEOPTFRIC') THEN
        TOPTFRIC=VAL_
      ELSE IF(ID_.EQ.'NONEGATIVEFRICTION') THEN
        TNONEGATIVEFRICTION=VAL_
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMS$SETL4')
      END IF
      RETURN
      END      
!
!     .................................................................. 
      SUBROUTINE ATOMS$GETL4(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL     ,INTENT(OUT):: VAL_
!     ******************************************************************
      IF(ID_.EQ.'STOP') THEN
        VAL_=TSTOP
      ELSE IF(ID_.EQ.'RANDOMIZE') THEN
        VAL_=TRANDOMIZE
      ELSE IF(ID_.EQ.'MOVE') THEN
        VAL_=TDYN
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMS$SETL4')
      END IF
      RETURN
      END      
!
!     .................................................................. 
      SUBROUTINE ATOMS$REPORT(NFIL)
!     ******************************************************************
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)               :: KELVIN
      INTEGER(4)            :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK) 
      IF(THISTASK.NE.1) RETURN
      CALL CONSTANTS('KB',KELVIN)
      CALL REPORT$TITLE(NFIL,'ATOMS')
      CALL REPORT$I4VAL(NFIL,'NUMBER OF ATOMS',NAT,' ')
      IF(START) THEN
        CALL REPORT$STRING(NFIL,'INITIAL POSITIONS TAKEN FROM STRUCTURE FILE')
      ELSE
        CALL REPORT$STRING(NFIL,'INITIAL POSITIONS TAKEN FROM RESTART FILE')
      END IF
      IF(.NOT.TDYN) THEN
        CALL REPORT$STRING(NFIL,'ATOMS ARE NOT PROPAGATED')
      ELSE
        CALL REPORT$STRING(NFIL,'ATOMS ARE PROPAGATED')
        IF(TSTOP) THEN
          CALL REPORT$CHVAL(NFIL,'INITIAL VELOCITIESARE SET TO ','ZERO')
        END IF
        IF(TRANDOMIZE) THEN
          CALL REPORT$R8VAL(NFIL,'INITIAL VELOCITIES RANDOMIZED WITH' &
     &                           ,AMPRE/KELVIN,'KELVIN')
        END IF
        IF(ANNER.NE.0.D0) THEN
          CALL REPORT$R8VAL(NFIL,'FRICTION',ANNER,' ')
        END IF
        IF(TNONEGATIVEFRICTION) THEN
          CALL REPORT$STRING(NFIL,'NEGATIVE FRICTIONS REMOVED')
        END IF
        IF(TOPTFRIC) THEN
          CALL REPORT$STRING(NFIL,'OPTIMUM FRICTION USED')
        END IF
      END IF
      RETURN
      END      
!
!     ..................................................................
      SUBROUTINE ATOMLIST$REPORT(NFIL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: ISP,IAT,ISP1
      CHARACTER(LEN=100)    :: STRING
      REAL(8)               :: U         ! MASS UNIT C12/12
      REAL(8)               :: ANGSTROM  ! ANGSTROM/ABOHR
      INTEGER(4)            :: NTASKS,THISTASK
      REAL(8)               :: EFFEMASS(NAT)
      REAL(8)               :: EMASS,EMASSCG2
      REAL(8)               :: RBAS(3,3)
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK) 
      IF(THISTASK.NE.1) RETURN
      CALL CONSTANTS('U',U)
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      CALL WAVES$GETR8('EMASS',EMASS)
      CALL WAVES$GETR8('EMASSCG2',EMASSCG2)
      CALL ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
      CALL CELL$GETR8A('T(0)',9,RBAS)
!
      CALL REPORT$TITLE(NFIL,'ATOMLIST REPORT')
      WRITE(NFIL,FMT='("T1[ANGSTROM]=",3F10.6/"T2[ANGSTROM]=",3F10.6/"T3[ANGSTROM]=",3F10.6)')RBAS/ANGSTROM
      WRITE(STRING,FMT='(T1,A,T15,A,T45,A,T53,A,T65,A,T80,A)') &
     &     'NAME','POSITION[ANGSTROM]','M[U]','MPSI_EFF[U]','Q[E]' &
     &    ,'FORCE[MH/ABOHR]'
      WRITE(NFIL,FMT='(A)')TRIM(STRING)
      ISP=0
      DO IAT=1,NAT
        ISP1=ISPECIES(IAT)
        IF(ISP1.NE.ISP) THEN
          STRING=' '
          ISP=ISP1
        END IF
        STRING=' '
        WRITE(STRING(1:9),FMT='(A)')NAME(IAT)(1:9)
        WRITE(STRING(10:42),FMT='("(",F9.5,",",F9.5,",",F9.5,")")') &
     &                                   R0(:,IAT)/ANGSTROM
        WRITE(STRING(43:52),FMT='(F8.4)')RMASS(IAT)/U
        WRITE(STRING(53:64),FMT='(F8.4)')EFFEMASS(IAT)/U
        WRITE(STRING(63:72),FMT='(F8.5)')-CHARGE(IAT)
        IF(ALLOCATED(FORCE).AND.TDYN.AND..NOT.TSTOP) THEN
          WRITE(STRING(73:100),FMT='("(",F7.2,",",F7.2,",",F7.2,")")') &
     &                        FORCE(:,IAT)*1.D+3
        END IF
        WRITE(NFIL,FMT='(A)')TRIM(STRING)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMS$INITIALIZE
!     **************************************************************************
!     **  READ ATOMIC POSITIONS FROM ATOMLIST                                 **
!     **************************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      INTEGER(4)           :: IAT,I
      REAL(8)              :: EKIN
      REAL(8)              :: EFFEMASS(NAT)
      REAL(8)              :: EMASS,EMASSCG2
      REAL(8)              :: SVAR
      REAL(8)              :: PROTONM
      LOGICAL(4)           :: TERR,TCHK
      INTEGER(4)           :: NFREE
!     **************************************************************************
                              CALL TRACE$PUSH('ATOMS$INITIALIZE')
!     CALL SHADOW$INITIALIZE
!     CALL SHADOW$OPTIMIZE(100,TOL,TCHK)
      DO IAT=1,NAT
        DO I=1,3
          FORCE(I,IAT)=0.D0
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  FIX ATOM THERMOSTAT                                                 ==
!     ==========================================================================
      CALL ATOMLIST$NATOM(NAT)
      CALL CONSTRAINTS$NFREE(NAT,NFREE)
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK)CALL THERMOSTAT$SCALEGFREE(REAL(NFREE,KIND=8))
!
!     ==========================================================================
!     ==   INITIALIZE THERMOSTATS FOR THE WAVE FUNCTIONS                      ==
!     ==========================================================================
      CALL WAVES$GETR8('EMASS',EMASS)
      CALL WAVES$GETR8('EMASSCG2',EMASSCG2)
      CALL ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
!
!     ==========================================================================
!     ==   CHECK IF REDUCED MASS IS NEGATIVE                                  ==
!     ==========================================================================
      TERR=.FALSE.
      DO IAT=1,NAT
        IF(EFFEMASS(IAT).GT.0.7D0*RMASS(IAT)) THEN
          IF(.NOT.TERR) THEN
            CALL ERROR$MSG('EFFECTIVE MASS OF WAVE FUNCTIONS ...')
            CALL ERROR$MSG('... MUST ACCOUNT LESS THAN 70% OF TOTAL MASS.')
            CALL ERROR$MSG('REDUCE WAVE FUNCTION MASS OR INCREASE NUCLEAR MASS')
            TERR=.TRUE.
          END IF
          CALL CONSTANTS('U',PROTONM)
          CALL ERROR$CHVAL('ATOM',NAME(IAT))
          CALL ERROR$R8VAL('MASS[U]     ',RMASS(IAT)/PROTONM)
          CALL ERROR$R8VAL('RED. MASS[U]',(RMASS(IAT)-EFFEMASS(IAT))/PROTONM)
          CALL ERROR$R8VAL('EFF. MASS[U]',EFFEMASS(IAT)/PROTONM)
          CALL ERROR$R8VAL('EMASS       ',EMASS)
          CALL ERROR$R8VAL('EMASSCG2    ',EMASSCG2)
          CALL ERROR$R8VAL('PSG2        ',PSG2(IAT))
          CALL ERROR$R8VAL('PSG4        ',PSG4(IAT))
        END IF
      ENDDO
      IF(TERR) THEN
        CALL ERROR$MSG('SOLUTION: REDUCE WAVE FUNCTION MASS')
        CALL ERROR$MSG('OR SET PS<G2> AND PS<G4> TO ZERO')
        CALL ERROR$STOP('ATOMS$INITIALIZE')
      END IF
!
!     ==========================================================================
!     ==   ADJUST KINETIC ENERGY OF THE WAVE FUNCTION THERMOSTAT              ==
!     ==========================================================================
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK) THEN
        CALL THERMOSTAT$SELECT('ATOMS')
        CALL THERMOSTAT$GETL4('ON',TCHK)
        IF(TCHK) THEN
          CALL THERMOSTAT$GETR8('TARGET',EKIN)
          SVAR=0.D0
          DO IAT=1,NAT
            SVAR=SVAR+EFFEMASS(IAT)/RMASS(IAT)
          ENDDO
          SVAR=SVAR/REAL(NAT,KIND=8)*EKIN
          IF(SVAR.EQ.0.D0) THEN
            CALL ERROR$MSG('TEMPERATURE FOR WAVE FUNCTION THERMOSTAT IS ZERO')
            CALL ERROR$MSG('PROBABLY !STRUCTURE!SPECIES:PS<G2> IS NOT SPECIFIED')
            CALL ERROR$STOP('ATOMS$INITIALIZE')
          END IF
          CALL THERMOSTAT$SELECT('WAVES')
          CALL THERMOSTAT$GETR8('TARGET',EKIN)
          CALL THERMOSTAT$SCALEGFREE(SVAR/EKIN)
          CALL THERMOSTAT$SETR8('GFREE',0.D0)
        ELSE
!!$          CALL ERROR$MSG('WAVE FUNCTION THERMOSTAT CANNOT BE USED')
!!$          CALL ERROR$MSG('WITHOUT ATOM THERMOSTAT')
!!$          CALL ERROR$STOP('ATOMS$INITIALIZE')
        END IF
      END IF
      STATEOFTHIS='INITIALIZED'
                               CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMS$PROPAGATE()
!     **************************************************************************
!     **  PROPAGATE: PROPAGATES ATOMIC POSITIONS                              **
!     **************************************************************************
      USE ATOMS_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)     :: IAT
      REAL(8)        :: EMASS,EMASSCG2
      REAL(8)        :: EFFEMASS(NAT)
      REAL(8)        :: REDRMASS(NAT)
      REAL(8)        :: ANNER1(NAT)
      LOGICAL(4)     :: TSTRESS
      REAL(8)        :: CELLFRIC(3,3)
      REAL(8)        :: RBAS(3,3)
      REAL(8)        :: SVAR
!     **************************************************************************
                              CALL TRACE$PUSH('ATOMS$PROPAGATE')
!!$PRINT*,'FORCES FROM ATOMS$PROPAGATE',MAXVAL(ABS(FORCE))
!!$DO IAT=1,NAT
!!$  WRITE(*,FMT='("FORCE ",I3,3F15.10)')IAT,FORCE(:,IAT)
!!$ENDDO
! 
!     **************************************************************************
!     ==========================================================================
!     == CONTROL AND CHANGE STATE OF THIS                                     ==
!     ==========================================================================
      IF(STATEOFTHIS.NE.'INITIALIZED'.AND.STATEOFTHIS.NE.'SWITCHED') THEN
        CALL ERROR$MSG('ATOMS OBJECT IS IN THE INCORRECT STATE')
        CALL ERROR$MSG('FOR APLICATION OF CONSTRAINTS')
        CALL ERROR$STOP('ATOMS$PROPAGATE')
      END IF
      STATEOFTHIS='PROPAGATED WITHOUT CONSTRAINTS'
! 
!     ==========================================================================
!     == FREEZE ATOMIC POSITIONS                                              ==
!     ==========================================================================
      IF(.NOT.TDYN) THEN
        RP(:,:)=R0(:,:)
        RM(:,:)=R0(:,:)
        CALL TRACE$POP ;RETURN
      END IF
! 
!     ==========================================================================
!     == CALCULATE EFFECTIVE WAVE FUNCTION MASS AND REDUCED MASS              ==
!     ==========================================================================
      CALL WAVES$GETR8('EMASS',EMASS)
      CALL WAVES$GETR8('EMASSCG2',EMASSCG2)
      CALL ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
      DO IAT=1,NAT
        REDRMASS(IAT)=RMASS(IAT)-EFFEMASS(IAT)
      ENDDO
! 
!     ==========================================================================
!     ==  DETERMINE FRICTION VECTOR                                           ==
!     ==========================================================================
      IF(.NOT.ALLOCATED(ANNERVEC0))ALLOCATE(ANNERVEC0(NAT))
      IF(TOPTFRIC) THEN
        IF(AOPT2AV.GT.0.D0) THEN
          SVAR=SQRT(AOPT1AV/AOPT2AV)
        ELSE
          SVAR=ANNER
        END IF
        ANNER=SVAR
PRINT*,'ATOMS: OPT.FRICTION ',SVAR,AOPT1AV,AOPT2AV,MIXAOPT
!       == CORRECT FRICTION FOR ELECTRON FRICTION  =============================
        ANNERVEC0(:)=(RMASS(:)*SVAR-EFFEMASS(:)*ANNEE)/REDRMASS(:)
      ELSE
!        == CORRECT FRICTION FOR ELECTRON FRICTION  ============================
        ANNERVEC0(:)=(RMASS(:)*ANNER-EFFEMASS(:)*ANNEE)/REDRMASS(:)
      END IF
!
!     == ENSURE THAT THE FRICTION IS ALWAYS POSITIVE. THE COMPENSATION ==
!     == FOR THE DRAG BY WAVE FUNCTION CLOUD CAN LEAD TO NEGATIVE FRICTION.
!     == DO NOT USE THIS OPTION WITH A THERMOSTAT!!
      IF(TNONEGATIVEFRICTION) THEN
        ANNERVEC0(:)=MAX(ANNERVEC0,0.D0)
      END IF
! 
!     ==========================================================================
!     == STOP ATOMIC MOTION                                                   ==
!     ==========================================================================
      IF(TSTOP) THEN
        RM(:,:)=R0(:,:)
        TSTOP=.FALSE.
      END IF
!
!     ==========================================================================
!     == RANDOMIZE VELOCITIES                                                 ==
!     ==========================================================================
      IF(TRANDOMIZE) THEN
        CALL ATOMS_RANDOMIZEVELOCITY(NAT,RMASS,RM,AMPRE,DELT)
        CALL MPE$BROADCAST('MONOMER',1,RM)
        TRANDOMIZE=.FALSE.
      END IF 
! 
!     ==========================================================================
!     ==  SET REFERENCE STRUCTURE  FOR CONSTRAINTS                            ==
!     ==========================================================================
      IF(.NOT.TCONSTRAINTREFERENCE) THEN
        CALL CELL$GETR8A('T(0)',9,RBAS)
        CALL CONSTRAINTS$SETREFERENCE(RBAS,NAT,R0,RM,REDRMASS,DELT)
        TCONSTRAINTREFERENCE=.TRUE.
      END IF
! 
!     == PRECONDITIONING =======================================================
!     CALL SHADOW$PRECONDITION(NAT,RMASS,R0,FORCE)
! 
!     == SYMMETRIZE FORCE ======================================================
      CALL SYMMETRIZE$FORCE(NAT,FORCE)
! 
!     ==========================================================================
!     ==  PROPAGATE ATOMS                                                     ==
!     ==========================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
      IF(TSTRESS) CALL CELL$GETR8A('FRICMAT',9,CELLFRIC)
      CALL ATOMS_PROPAGATE(NAT,DELT,REDRMASS,ANNERVEC0 &
   &                      ,FORCE,R0,RM,RP,TSTRESS,CELLFRIC)
      IF(.NOT.TDYN) THEN
        RP=R0
      END IF
!
!     ==========================================================================
!     == TAKE CARE OF LINK BONDS AS CONSTRAINTS                               ==
!     ==========================================================================
      CALL ATOMS_FRICTIONARRAY(NAT,ANNER,ANNEE,RMASS,EFFEMASS,ANNER1)
      CALL QMMM$PROPAGATE(NAT,DELT,ANNER,REDRMASS,ANNER1,RM,R0,RP)
!
!     ==========================================================================
!     == WARMUP SYSTEM                                                        ==
!     ==========================================================================
      CALL PAW_WARMUP_APPLY
                            CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMS$PROPAGATE_DIMER()
!     ******************************************************************
!     **  PROPAGATE: PROPAGATES ATOMIC POSITIONS                      **
!     ******************************************************************
      USE ATOMS_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)     :: IAT
      REAL(8)        :: EMASS,EMASSCG2
      REAL(8)        :: EFFEMASS(NAT)
      REAL(8)        :: REDRMASS(NAT)
      REAL(8)        :: ANNER1(NAT)
      LOGICAL(4)     :: TCHK,TCHK2,TCHK3,TCHK4,TCHK5
      LOGICAL(4)     :: TSTRESS
      REAL(8)        :: CELLFRIC(3,3)
      REAL(8)        :: CELLKIN(3,3)  ! SHOULD BE REMOVED 
      REAL(8)        :: RBAS(3,3)
      REAL(8)        :: SVAR
      REAL(8)        :: X0(NAT*3),XP(NAT*3) !MASSWEIGHTED COORDINATES FOR THE DIMER
!     ******************************************************************
                              CALL TRACE$PUSH('ATOMS$PROPAGATE')
DO IAT=1,NAT
  WRITE(*,FMT='("FORCE ",I3,3F15.10)')IAT,FORCE(:,IAT)
ENDDO
     CALL DIMER$GETL4('DIMER',TCHK)
     IF(TCHK) CALL DIMER$GETL4('PLACEDIMER',TCHK)
     IF(.NOT.TDYN.AND..NOT.TCHK) THEN
        RP(:,:)=R0(:,:)
        RM(:,:)=R0(:,:)
        CALL TRACE$POP ;RETURN
     END IF
! 
!     ==================================================================
!     == CALCULATE EFFECTIVE AND REDUCED MASS
!     ==================================================================
      CALL WAVES$GETR8('EMASS',EMASS)
      CALL WAVES$GETR8('EMASSCG2',EMASSCG2)
      CALL ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
      DO IAT=1,NAT
        REDRMASS(IAT)=RMASS(IAT)-EFFEMASS(IAT)
      ENDDO
! 
!     ==================================================================
!     == STOP ATOMIC MOTION
!     ==================================================================
      IF(TSTOP) THEN
        RM(:,:)=R0(:,:)
        TSTOP=.FALSE.
      END IF
!
!     ==================================================================
!     == RANDOMIZE VELOCITIES
!     ==================================================================
      IF(TRANDOMIZE) THEN
PRINT*,'TRANDOMIZE ',TRANDOMIZE
        CALL ATOMS_RANDOMIZEVELOCITY(NAT,RMASS,RM,AMPRE,DELT)
        CALL MPE$BROADCAST('MONOMER',1,RM)
        TRANDOMIZE=.FALSE.
      END IF 
! 
!     ==================================================================
!     ==  SET REFERENCE STRUCTURE  FOR CONSTRAINTS                    ==
!     ==================================================================
      IF(.NOT.TCONSTRAINTREFERENCE) THEN
        CALL CELL$GETR8A('T(0)',9,RBAS)
        CALL CONSTRAINTS$SETREFERENCE(RBAS,NAT,R0,RM,REDRMASS,DELT)
        TCONSTRAINTREFERENCE=.TRUE.
      END IF

!     ================================================================
!     ==  DETERMINE FRICTION VECTOR                                 ==
!     ================================================================
      IF(.NOT.ALLOCATED(ANNERVEC0))ALLOCATE(ANNERVEC0(NAT))
      IF(TOPTFRIC) THEN
        IF(AOPT2AV.GT.0.D0) THEN
          SVAR=SQRT(AOPT1AV/AOPT2AV)
        ELSE
          SVAR=ANNER
        END IF
        ANNER=SVAR
PRINT*,'ATOMS: OPT.FRICTION ',SVAR,AOPT1AV,AOPT2AV,MIXAOPT
!       == CORRECT FRICTION FOR ELECTRON FRICTION  =====================
        ANNERVEC0(:)=(RMASS(:)*SVAR-EFFEMASS(:)*ANNEE)/REDRMASS(:)
      ELSE
!        == CORRECT FRICTION FOR ELECTRON FRICTION  =====================
        ANNERVEC0(:)=(RMASS(:)*ANNER-EFFEMASS(:)*ANNEE)/REDRMASS(:)
      END IF
!
!     == ENSURE THAT THE FRICTION IS ALWAYS POSITIVE. THE COMPENSATION ==
!     == FOR THE DRAG BY WAVE FUNCTION CLOUD CAN LEAD TO NEGATIVE FRICTION.
!     == DO NOT USE THIS OPTION WITH A THERMOSTAT!!
      IF(TNONEGATIVEFRICTION) THEN
        ANNERVEC0(:)=MAX(ANNERVEC0,0.D0)
      END IF



! 
!     == PRECONDITIONING =============================================
!     CALL SHADOW$PRECONDITION(NAT,RMASS,R0,FORCE)
! 
!     == SYMMETRIZE FORCE ============================================
      CALL SYMMETRIZE$FORCE(NAT,FORCE)
! 
!     ================================================================
!     ==  PROPAGATE ATOMS                                           ==
!     ================================================================
      CELLFRIC=0.D0
      TSTRESS=.FALSE.
      !     CALL CELL$GETL4('MOVE',TSTRESS)
      !     IF(TSTRESS) CALL CELL$GETR8A('FRICMAT',9,CELLFRIC)


      !ALEXP-DIMER 
      !THE OLD CALL WAS JUST  THE ONE BELOW
      !CALL ATOMS_PROPAGATE(NAT,DELT,ANNER,ANNEE,RMASS,EFFEMASS &
      !&                      ,FORCE,R0,RM,RP,TSTRESS,CELLFRIC)


      !THIS IS A HARDWIRED REREAD OF SOME CNTL PARAMETERS
      !CHECK OR DELETE THIS
      !CALL MPE$QUERY('~',NTASKS,THISTASK)
      !IF(THISTASK.EQ.1) CALL DIMER$REREAD()
      !THIS IS HARDWIRED REREAD OF SOME CNTL PARAMETERS


      CALL DIMER$GETL4('DIMER',TCHK)
      CALL DIMER$GETL4('STRETCH',TCHK2)
      CALL DIMER$GETL4('PLACEDIMER',TCHK3)
      CALL DIMER$GETL4('DIMERFOLLOWDOWN',TCHK4)
      CALL DIMER$GETL4('KDLENGTH',TCHK5)
      !   CALL DIMER$GETL4('OUTPUT?',TCHK5) INCLUDE THIS IF PETER WANTS THE OUTPUT
      PRINT*,TCHK,TCHK2,TCHK3,TCHK4,TCHK5

      IF(TCHK.AND.TCHK2.AND.TCHK5) THEN
         CALL ERROR$MSG('WE CAN NOT STRETCH THE DIMER AND KEEP ITS LENGTH')
         CALL ERROR$MSG('AT THE SAME TIME - ADJUST YOUR .CNTL-FILE')
         CALL ERROR$STOP('ATOMS$PROPAGATE')
      END IF

      IF(TCHK.AND..NOT.TCHK2.AND..NOT.TCHK3.AND..NOT.TCHK4.AND..NOT.TCHK5) THEN
         !NORMAL DIMER PROPAGATION
         CALL DIMER$PROPAGATE(NAT,DELT,ANNER,ANNEE,RMASS,EFFEMASS &
              &                   ,FORCE,R0,RM,RP,TSTRESS,CELLFRIC,CELLKIN)

      ELSE IF(TCHK.AND..NOT.TCHK2.AND..NOT.TCHK4.AND.(TCHK3.OR.TCHK5)) THEN
         CALL DIMER$PROPAGATE(NAT,DELT,ANNER,ANNEE,RMASS,EFFEMASS &
              &                   ,FORCE,R0,RM,RP,TSTRESS,CELLFRIC,CELLKIN)
         !SET THE CORRECT POSITIONS (WE GOT THEM IN THE RP POSITIONS)
         R0(:,:)=RP(:,:)
         RM(:,:)=RP(:,:)
         CALL DIMER$SETL4('PLACEDIMER',.FALSE.)
         CALL DIMER$SETL4('KDLENGTH',.FALSE.)
         !THIS WAS IT FOR THIS STEP
         CALL TRACE$POP ;RETURN

      ELSE IF (.NOT.TCHK.AND.TCHK3) THEN
         CALL ERROR$MSG('I CAN NOT PLACE DIMER IN A NON DIMER CALCULATION')
         CALL ERROR$L4VAL('PLACEDIMER',TCHK3)
         CALL ERROR$L4VAL('DIMER',TCHK)
         CALL ERROR$STOP('ATOMS$PROPAGATE')
    
      ELSE IF(TCHK.AND..NOT.TCHK4.AND.TCHK2) THEN
         !LET THE SECOND IMAGE FOLLOWDOWN UNTIL (R2-R1)^2=D^2 
         IF(TCHK3) THEN
            CALL DIMER$GET_MASSWEIGHTED(NAT*3,R0,X0)
            CALL DIMER$STRETCH(NAT,X0,XP) !THE SECOND XP IS ARBITRARY
            CALL DIMER$GET_UNMASSWEIGHTED(NAT*3,X0,R0)
            CALL DIMER$SETL4('PLACEDIMER',.FALSE.)
         ELSE 
            CALL ERROR$MSG('ARGUMENT LIST OF ATOMS_PROPAGATE INCONSISTENT')
            CALL ERROR$MSG('FIX BUG FIRST')
            CALL ERROR$STOP('ATOMS$PROPAGATE_DIMER')
! THIS LINE MUST BE CORRECTED!!!!
!            CALL ATOMS_PROPAGATE(NAT,DELT,ANNER,ANNEE,RMASS,EFFEMASS &
!                 &                   ,FORCE,R0,RM,RP,TSTRESS,CELLFRIC)
            
            CALL DIMER$GET_MASSWEIGHTED(NAT*3,R0,X0)
            CALL DIMER$GET_MASSWEIGHTED(NAT*3,RP,XP)
            CALL DIMER$STRETCH(NAT,XP,X0)
            CALL DIMER$GET_UNMASSWEIGHTED(NAT*3,XP,RP)!ONLY THIS BECAUSE R0 IS NOT CHANGED!
         END IF


      ELSE 
         !=== AS USUAL OR DIMER IMAGES DECOUPLED FOLLOW DOWN
         CALL ERROR$MSG('THIS COMES FROM THE DIMER MERGE AND SHOULD NOT HAPPEN')
         CALL ERROR$STOP('ATOMS$PROPAGATE_DIMER')
!        CALL ATOMS_PROPAGATE(NAT,DELT,ANNER,ANNEE,RMASS,EFFEMASS &
!              &                   ,FORCE,R0,RM,RP,TSTRESS,CELLFRIC)


         !WRITES THE ENERGY FOR A DIMER FOLLWODOWN CALCULATION 
         !IF PETER DOES NOT LIKE IT DELETE THE CALL
         IF(TCHK) THEN
            CALL DIMER$WRITEENERGYTRA(NAT*3,RP)
            CALL MPE$SYNC('~')
         END IF

      END IF

      !ALEXP-DIMER END 


!
!     ==================================================================
!     == TAKE CARE OF LINK BONDS AS CONSTRAINTS                       ==
!     ==================================================================
      CALL ATOMS_FRICTIONARRAY(NAT,ANNER,ANNEE &
     &                 ,RMASS,EFFEMASS,ANNER1)
      CALL QMMM$PROPAGATE(NAT,DELT,ANNER,REDRMASS,ANNER1,RM,R0,RP)
!
!     ==================================================================
!     == WARMUP SYSTEM                                                ==
!     ==================================================================
      CALL PAW_WARMUP_APPLY
DO IAT=1,NAT
  WRITE(*,FMT='("RM ",I3,3F15.10)')IAT,RM(:,IAT)
  WRITE(*,FMT='("R0 ",I3,3F15.10)')IAT,R0(:,IAT)
  WRITE(*,FMT='("RP ",I3,3F15.10)')IAT,RP(:,IAT)
ENDDO
                            CALL TRACE$POP
      RETURN
    END SUBROUTINE ATOMS$PROPAGATE_DIMER
!
!     ..................................................................
      SUBROUTINE ATOMS$CONSTRAINTS
!     ******************************************************************
!     **  PROPAGATE: PROPAGATES ATOMIC POSITIONS                      **
!     ******************************************************************
      USE ATOMS_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      REAL(8)        :: EMASS,EMASSCG2
      REAL(8)        :: EFFEMASS(NAT)
      REAL(8)        :: REDRMASS(NAT)
      LOGICAL(4)     :: TCHK
      REAL(8)        :: MAPTOCELL(3,3)
      REAL(8)        :: MAPTOCELLINV(3,3)
      INTEGER(4)     :: IAT,I,J
      LOGICAL(4)     :: TSTRESS
      REAL(8)        :: RP1(3,NAT)
      REAL(8)        :: STRESS1(3,3)
      REAL(8)        :: V(3)
      REAL(8)        :: RBAS(3,3)
!     ******************************************************************
      CALL DIMER$GETL4('DIMER',TCHK)
      IF(TCHK) THEN
         PRINT*,'********** WARNING FROM DIMER ****************'
         PRINT*,'THIS IS A DIMER CALCULATION'
         PRINT*,'IGNORING ALL CONSTRAINTS IN ATOMS$CONSTRAINTS'
         PRINT*,'********** WARNING FROM DIMER ****************'
         RETURN
      END IF
! 
!     ==================================================================
!     == CONTROL AND CHANGE STATE OF THIS                             ==
!     ==================================================================
      IF(STATEOFTHIS.NE.'PROPAGATED WITHOUT CONSTRAINTS') THEN
        CALL ERROR$MSG('ATOMS OBJECT IS IN THE INCORRECT STATE')
        CALL ERROR$MSG('FOR APLICATION OF CONSTRAINTS')
        CALL ERROR$STOP('ATOMS$CONSTRAINTS')
      END IF
      STATEOFTHIS='PROPAGATED WITH CONSTRAINTS'
                              CALL TRACE$PUSH('ATOMS$CONSTRAINTS')
! 
!     ==================================================================
!     == RETURN IF ATOMIC STRUCTURE IS STATIC                         ==
!     ==================================================================
      IF(.NOT.TDYN) THEN
        CALL TRACE$POP ;RETURN
      END IF
! 
!     ==================================================================
!     == CALCULATE EFFECTIVE AND REDUCED MASS
!     ==================================================================
      CALL WAVES$GETR8('EMASS',EMASS)
      CALL WAVES$GETR8('EMASSCG2',EMASSCG2)
      CALL ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
      DO IAT=1,NAT
        REDRMASS(IAT)=RMASS(IAT)-EFFEMASS(IAT)
      ENDDO
!
!     ==================================================================
!     == MAP NEW COORDINATES INTO THE NEXT UNIT CELL                  ==
!     ==================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
!PRINT*,'WARNING FROM ATOMS$CONSTRAINTS: TSTRESS SET TO FALSE'
!TSTRESS=.FALSE.
      IF(TSTRESS) THEN
        CALL CELL$GETR8A('MAPTOCELL',9,MAPTOCELL)
        DO IAT=1,NAT
          RP1(:,IAT)=MATMUL(MAPTOCELL,RP(:,IAT))
        ENDDO
      ELSE
        DO IAT=1,NAT
          RP1(:,IAT)=RP(:,IAT)
        ENDDO
      END IF
!
!     ==================================================================
!     == APPLY CONSTRAINTS                                            ==
!     ==================================================================
      CALL CELL$GETR8A('T(0)',9,RBAS)
      CALL CONSTRAINTS$APPLY(RBAS,NAT,R0,RP1,REDRMASS,FORCE,DELT)
      CALL SYMMETRIZE$POSITION(NAT,RP1)
!
!     ==================================================================
!     == MAP NEW COORDINATES BACK INTO THE CURRENT UNIT CELL          ==
!     ==================================================================
      IF(TSTRESS) THEN
        CALL LIB$INVERTR8(3,MAPTOCELL,MAPTOCELLINV)
        DO IAT=1,NAT
          RP(:,IAT)=MATMUL(MAPTOCELLINV,RP1(:,IAT))
        ENDDO
      ELSE
        DO IAT=1,NAT
          RP(:,IAT)=RP1(:,IAT)
        ENDDO
      END IF
                            CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMS$EKIN(EKIN_)
!     **************************************************************************
!     **  CALCULATE KINETIC ENERGY                                            **
!     **************************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: EKIN_
      REAL(8)             :: DEKIN
!     **************************************************************************
      CALL ATOMS_EKIN(NAT,DELT,RMASS,RM,RP,EKIN_)
      CALL QMMM$DEKIN(DELT,DEKIN)   !CORRECTION FOR LINK ATOMS
      EKIN_=EKIN_+DEKIN
!
!     ==================================================================
!     == ADD KINETIC FRICTION TO STRESS                               ==
!     ==================================================================
      CALL ATOMS$KINETICSTRESS()
!
      RETURN
      END
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMS$KINETICSTRESS()
!     **************************************************************************
!     **  CALCULATE IDEAL-GAS CONTRIBUTION TO THE INTERNAL STRESS             **
!     **  (=STRESS ORIGINATING FROM THE KINETIC ENERGY OF THE ATOMS)          **
!     **     
!     **  STRESS IS COMMUNICATED DIRECTLY TO THE CELL OBJECT
!     **************************************************************************
      USE ATOMS_MODULE, ONLY : NAT,RMASS,RM,RP,DELT,PSG2,PSG4
      IMPLICIT NONE
      REAL(8)             :: STRESS(3,3)
      REAL(8)             :: V(3)
      LOGICAL(4)          :: TSTRESS
      INTEGER(4)          :: IAT,I,J      
      REAL(8)             :: EMASS,EMASSCG2
      REAL(8)             :: EFFEMASS(NAT)   !EFFECTIVE MASS OF THE WAVE FUNCTIONS
      REAL(8)             :: REDRMASS(NAT)   !REDUCED MASS   
!     **************************************************************************
      CALL CELL$GETL4('MOVE',TSTRESS)
      IF(.NOT.TSTRESS) THEN
        STRESS=0.D0
        CALL CELL$SETR8A('KINSTRESS',9,STRESS)
        RETURN
      END IF
!
!     ==========================================================================
!     == CALCULATE REDUCED MASS                                               ==
!     ==========================================================================
      CALL WAVES$GETR8('EMASS',EMASS)
      CALL WAVES$GETR8('EMASSCG2',EMASSCG2)
      CALL ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
      DO IAT=1,NAT
        REDRMASS(IAT)=RMASS(IAT)-EFFEMASS(IAT)
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE STRESS                                                    ==
!     ==========================================================================
      STRESS(:,:)=0.D0
      DO IAT=1,NAT
        V(:)=(RP(:,IAT)-RM(:,IAT))/(2.D0*DELT)
        DO I=1,3
          DO J=1,3
            STRESS(I,J)=STRESS(I,J)+REDRMASS(IAT)*V(I)*V(J)
          ENDDO
        ENDDO
      ENDDO
!
      CALL CELL$SETR8A('KINSTRESS',9,STRESS)
      RETURN
      END
!
!     .................................................................. 
      SUBROUTINE ATOMS$EFFEKIN(EKIN_)
!     ******************************************************************
!     **  CALCULATE KINETIC ENERGY OF ISOLATED ATOMS ON THE           **
!     **  BORN-OPPENHEIMER SURFACE FROM THE VELOCITIES OF THE ATOMS   **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: EKIN_
      REAL(8)             :: EFFEMASS(NAT)
      REAL(8)             :: EMASS,EMASSCG2
!     ******************************************************************
      CALL WAVES$GETR8('EMASS',EMASS)
      CALL WAVES$GETR8('EMASSCG2',EMASSCG2)
      CALL ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
      CALL ATOMS_EKIN(NAT,DELT,EFFEMASS,RM,RP,EKIN_)
      RETURN
      END
!
!     .................................................................. 
      SUBROUTINE ATOMS$DEPOTDT(DEPOTDT)
!     ******************************************************************
!     **  CALCULATE KINETIC ENERGY                                    **
!     ******************************************************************!
      USE ATOMS_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: DEPOTDT
      INTEGER(4)          :: I,IAT
!     ******************************************************************
      DEPOTDT=0.D0
      DO IAT=1,NAT
        DO I=1,3
          DEPOTDT=DEPOTDT+FORCE(I,IAT)*(RP(I,IAT)-RM(I,IAT))
        ENDDO
      ENDDO
      DEPOTDT=DEPOTDT/(2.D0*DELT)
      RETURN
      END
!
!     .................................................................. 
      SUBROUTINE ATOMS$TEMPERATURE(TEMPERATURE_)
!     ******************************************************************
!     **  CALCULATE KINETIC ENERGY                                    **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: TEMPERATURE_
      REAL(8)             :: EKIN
      REAL(8)             :: DEKIN
      INTEGER(4)          :: NFREE
!     ******************************************************************
      CALL CONSTRAINTS$NFREE(NAT,NFREE)
      IF(NFREE.EQ.0) THEN
        TEMPERATURE_=0.D0
      ELSE
        CALL ATOMS_EKIN(NAT,DELT,RMASS,RM,RP,EKIN)
        CALL QMMM$DEKIN(DELT,DEKIN)   !CORRECTION FOR LINK ATOMS
        TEMPERATURE_=(EKIN+DEKIN)/(0.5D0*DBLE(NFREE))
      END IF
      RETURN
      END
!
!     .................................................................. 
      SUBROUTINE ATOMS$FORCECRITERION(FAV,FMAX)
!     ******************************************************************
!     **  CALCULATES THE ABSOLUTE VALUE OF THE FORCE                  **
!     **                                                              **
!     **  CALL AFTER CONSTRAINTS AND BEFORE SWITCH                    **
!     ****************************************************************** 
      USE ATOMS_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: FAV
      REAL(8),INTENT(OUT) :: FMAX
      REAL(8)             :: F(3,NAT)
      INTEGER(4)          :: IAT
      REAL(8)             :: SVAR
      REAL(8)             :: MRED  ! REDUCED MASS
      REAL(8)             :: A0    ! EFFECTIVE FRICTION
      REAL(8)             :: EMASS,EMASSCG2
      REAL(8)             :: EFFEMASS(NAT)
      LOGICAL(4)          :: TCHK
      REAL(8)   ,PARAMETER:: R8LARGE=1.D+20
!     ****************************************************************** 

      CALL DIMER$GETL4('DIMER',TCHK)
      IF(TCHK) THEN
         PRINT*,'********** WARNING FROM DIMER ****************'
         PRINT*,'THIS IS A DIMER CALCULATION'
         PRINT*,'IGNORING FORCECRITERION IN ATOMS$FORCECRITERION'
         PRINT*,'********** WARNING FROM DIMER ****************'
         FAV=R8LARGE
         FMAX=R8LARGE
         RETURN
      END IF

! 
!     ==================================================================
!     == CONTROL STATE OF THIS                                        ==
!     ==================================================================
      IF(STATEOFTHIS.NE.'PROPAGATED WITH CONSTRAINTS') THEN
        CALL ERROR$MSG('ATOMS OBJECT IS IN THE INCORRECT STATE')
        CALL ERROR$CHVAL('STATEOFTHIS',STATEOFTHIS)
        CALL ERROR$STOP('ATOMS$FORCECRITERION')
      END IF
!
!     == RETURN FOR STATIC CALCULATION. FORCES ARE DETERMINED FROM TRAJECTORY
      IF(.NOT.TDYN) THEN
        FAV=R8LARGE
        FMAX=R8LARGE
        RETURN
      END IF

      IF(.NOT.ALLOCATED(ANNERVEC0)) THEN
        CALL ERROR$MSG('ANNERVEC0 NOT YET ALLOCATED')
        CALL ERROR$STOP('ATOMS$FORCECRITERION')
      END IF
! 
!     ==================================================================
!     == DETERMINE MASS OF WAVE FUNCTION CLOUD                        ==
!     ==================================================================
      CALL WAVES$GETR8('EMASS',EMASS)
      CALL WAVES$GETR8('EMASSCG2',EMASSCG2)
      CALL ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
! 
!     ==================================================================
!     == DETERMINE MAXIMUM AND AVERAGE VALUE OF THE FORCE             ==
!     ==================================================================
      FAV=0.D0
      FMAX=0.D0
      DO IAT=1,NAT
        MRED=RMASS(IAT)-EFFEMASS(IAT)
        A0=ANNERVEC0(IAT)
        F(:,:)=(1.D0+A0)*RP(:,:)-2.D0*R0(:,:)+(1.D0-A0)*RM(:,:)
        F(:,IAT)=F(:,IAT)*MRED
        SVAR=DOT_PRODUCT(F(:,IAT),F(:,IAT))
        FAV=FAV+SVAR
        FMAX=MAX(FMAX,SVAR)
      ENDDO
      FAV=FAV/DELT**2
      FMAX=FMAX/DELT**2
      FAV=SQRT(FAV/REAL(NAT,KIND=8))
      FMAX=SQRT(FMAX)
      RETURN
      END
!
!     .................................................................. 
      SUBROUTINE ATOMS_OPTFRICTION(AOPT1,AOPT2,NAT,RP,R0,RM,R2M &
     &                     ,ANNERVEC0,ANNERVECM,REDUCEDMASS)
!     ******************************************************************
!     **  CALCULATES THE ABSOLUTE VALUE OF THE FORCE                  **
!     **                                                              **
!     **  CALL AFTER CONSTRAINTS AND BEFORE SWITCH                    **
!     ****************************************************************** 
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT):: AOPT1  ! OPTIMUM FRICTION =SQRT(AOPT1/AOPT2)
      REAL(8)   ,INTENT(OUT):: AOPT2  ! OPTIMUM FRICTION =SQRT(AOPT1/AOPT2)
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RP(3,NAT)  ! POSITIONS AT NEXT TIME STEP
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)  ! POSITIONS AT CURRENT TIME STEP
      REAL(8)   ,INTENT(IN) :: RM(3,NAT)  ! POSITION AT PREVIOUS TIME STEP
      REAL(8)   ,INTENT(IN) :: R2M(3,NAT) ! POSITION TWO STEPS BEFORE
      REAL(8)   ,INTENT(IN) :: ANNERVEC0(NAT)  ! CURRENT FRICTION 
      REAL(8)   ,INTENT(IN) :: ANNERVECM(NAT)  ! PREVIOUS FRICTION 
      REAL(8)   ,INTENT(IN) :: REDUCEDMASS(NAT) ! REDUCED MASS 
      INTEGER(4)            :: IAT
      REAL(8)               :: ASUM,BSUM
      REAL(8)               :: A0 ! EFFECTIVE FRICTION/CURRENT TIME STEP
      REAL(8)               :: AM ! EFFECTIVE FRICTION/PREVIOUS TIME STEP
      REAL(8)               :: DF(3),DX(3),MDX(3)
!     ****************************************************************** 
      ASUM=0.D0
      BSUM=0.D0
      DO IAT=1,NAT
        A0=ANNERVEC0(IAT)
        AM=ANNERVECM(IAT)
        DF(:)=(1.D0+A0)*RP(:,IAT)-(3.D0+AM)*R0(:,IAT) &
     &       +(3.D0-A0)*RM(:,IAT)-(1.D0-AM)*R2M(:,IAT)
        DX(:)=(R0(:,IAT)-RM(:,IAT))
        MDX(:)=REDUCEDMASS(IAT)*DX(:)
        ASUM=ASUM+DOT_PRODUCT(DF,MDX)
        BSUM=BSUM+DOT_PRODUCT(DX,MDX)
      ENDDO
      ASUM=ABS(ASUM)
      AOPT1=ASUM
      AOPT2=BSUM
      RETURN
      END
!
!     .................................................................. 
      SUBROUTINE ATOMS$SWITCH()
!     ******************************************************************
!     **  SWITCH ATOMIC POSITIONS  R(+)->R(0)->R(M)                   **
!     ****************************************************************** 
      USE ATOMS_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      REAL(8)     :: MAPTOCELL(3,3)
      LOGICAL(4)  :: TSTRESS
      INTEGER(4)  :: IAT
      REAL(8)     :: DRPM(3)
      REAL(8)     :: EFFEMASS(NAT)
      REAL(8)     :: EMASS,EMASSCG2
      REAL(8)     :: AOPT1,AOPT2
      LOGICAL(4)  :: TCHK
!     ******************************************************************
!
!     ==================================================================
!     == CONTROL STATE OF THIS                                        ==
!     ==================================================================
      CALL DIMER$GETL4('DIMER',TCHK)
      IF(TCHK) THEN
         PRINT*,'********** WARNING FROM DIMER ****************'
         PRINT*,'THIS IS A DIMER CALCULATION'
         PRINT*,'IGNORING STATEOFTHIS IN ATOMS$SWITCH'
         PRINT*,'********** WARNING FROM DIMER ****************'
      ELSE
         IF(STATEOFTHIS.NE.'PROPAGATED WITH CONSTRAINTS') THEN
            CALL ERROR$MSG('ATOMS OBJECT IS IN THE INCORRECT STATE')
            CALL ERROR$CHVAL('STATEOFTHIS',STATEOFTHIS)
            CALL ERROR$STOP('ATOMS$SWITCH')
         END IF
      END IF

      STATEOFTHIS='SWITCHED'
!
!     ==================================================================
!     == WRAP UP VARIOUS STUFF                                        ==
!     ==================================================================
      IF(ALLOCATED(R2M)) THEN
        CALL WAVES$GETR8('EMASS',EMASS)
        CALL WAVES$GETR8('EMASSCG2',EMASSCG2)
        CALL ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
        CALL ATOMS_OPTFRICTION(AOPT1,AOPT2,NAT,RP,R0,RM,R2M &
     &                     ,ANNERVEC0,ANNERVECM,RMASS-EFFEMASS)
!
!       == PERFORM FLOATING AVERAGE OF NUMERATOR AND DENOMINATOR  ======
!       == OPTIMUM FRICTION FACTOR =SQRT(AOPT1AV/AOPT2AV)
        AOPT1AV=MIXAOPT*AOPT1+(1.D0-MIXAOPT)*AOPT1AV
        AOPT2AV=MIXAOPT*AOPT2+(1.D0-MIXAOPT)*AOPT2AV
      ELSE 
!       == OPTIMUM FRICTION IS IDENTICAL TO ACTUAL FRICTION.
!       == NUMERATOR AND DENOMINATOR ARE SET TO SMALL VALUES SO THAT THEY 
!       == ARE DOMINATED QUICKLY BY THE ACTUAL VALUES.
        AOPT2AV=1.D-9
        AOPT1AV=ANNER**2*AOPT2AV
      END IF
!
!     ==================================================================
!     == NOW PERFORM ACTUAL SWITCH                                    ==
!     ==================================================================
      FORCE(:,:)=0.D0
      CALL QMMM$SWITCH
!
!     ==================================================================
!     == TRANSFORM TRAJECTORY FOR THE NEW REFERENCE CELL              ==
!     ==================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
      IF(TSTRESS) THEN
        CALL CELL$GETR8A('MAPTOCELL',9,MAPTOCELL)
        DO IAT=1,NAT
          R0(:,IAT)=MATMUL(MAPTOCELL,R0(:,IAT))
          RM(:,IAT)=MATMUL(MAPTOCELL,RM(:,IAT))
          RP(:,IAT)=MATMUL(MAPTOCELL,RP(:,IAT))
        ENDDO
      END IF
!
!     ==========================================================================
!     == SYNCHRONIZE POSITIONS AMONG PROCESSES (MPI)                          ==
!     ==========================================================================
      CALL MPE$BROADCAST('MONOMER',1,RP)
      CALL MPE$BROADCAST('MONOMER',1,R0)
      CALL MPE$BROADCAST('MONOMER',1,RM)
      IF(ALLOCATED(ANNERVEC0))CALL MPE$BROADCAST('MONOMER',1,ANNERVEC0)
!
!     ==========================================================================
!     == SWITCH VARIABLES (ALSO FOR OPTIMUM FRICTION)                         ==
!     ==========================================================================
      IF(.NOT.TDYN) RETURN
!
      IF(.NOT.ALLOCATED(R2M))ALLOCATE(R2M(3,NAT))
      IF(.NOT.ALLOCATED(ANNERVECM))ALLOCATE(ANNERVECM(NAT))
      R2M(:,:)=RM(:,:)
      RM(:,:)=R0(:,:)
      R0(:,:)=RP(:,:)
      RP(:,:)=2.D0*R0(:,:)-RM(:,:)
      ANNERVECM(:)=ANNERVEC0(:)
!
      CALL CONSTRAINTS$SWITCH()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMS_RANDOMIZEVELOCITY(NAT,RMASS,RM,EBATH,DELT)
!     **************************************************************************
!     **                                                                      **
!     **  RANDOMIZED VELOCITY ACCORDING TO A TEMPERATURE                      **
!     **                                                                      **
!     **  REMARK: ANY CONSTRAINTS MUST BE ENFORCED SUBSEQUENTLY               **
!     **                                                                      **
!     **    EBATH    K_B*T WHERE THE T IS THE TEMPERATURE OF THE              **
!     **             HEATBATH                                                 **
!     **                                                                      **
!     **    THE TARGET KINETIC ENERGY IS 1.5*NAT*EBATH                        **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NAT        ! #(ATOMS)
      REAL(8)   ,INTENT(IN)   :: DELT       ! TIME STEP
      REAL(8)   ,INTENT(IN)   :: EBATH      ! K_BT OF THE HEAT BATH
      REAL(8)   ,INTENT(IN)   :: RMASS(NAT) ! ATOMIC MASSES
      REAL(8)   ,INTENT(INOUT):: RM(3,NAT)  ! R(-)
      LOGICAL(4),PARAMETER    :: TPR=.FALSE.
      LOGICAL(4),PARAMETER    :: TSETSEED=.TRUE.
      INTEGER(4)              :: IAT,I
      REAL(8)                 :: SVAR,RAN,SUM
      INTEGER                 :: NSEED
      INTEGER   ,ALLOCATABLE  :: SEED(:),OLDSEED(:)
!     **************************************************************************
!
!     ==========================================================================
!     == INITIALIZE THE RANDOM-NUMBER GENERATION WITH A DEFINED SEED TO MAKE  ==
!     == MAKE IT REPRODUCIBLE
!     ==========================================================================
      IF(TSETSEED) THEN
        CALL RANDOM_SEED(SIZE=NSEED)
        ALLOCATE(OLDSEED(NSEED))
        CALL RANDOM_SEED(GET=OLDSEED)
        ALLOCATE(SEED(NSEED))
!       __SEE CHAPTER 7 IN NUMERICAL RECIPES: THE ART OF SCIENTIFIC COMPUTING___
        SEED(1)=314159265
        DO I=2,NSEED
          SEED(I)=MOD(8121*SEED(I-1)+28411,134456)
        ENDDO
        CALL RANDOM_SEED(PUT=SEED)
        DEALLOCATE(SEED)
      ENDIF
!
!     ==========================================================================
!     == PROVIDE A RANDOM KICK (CHANGE VELOCITIES, I.E. R(-))                 ==
!     ==========================================================================
      SUM=0.D0
      DO IAT=1,NAT
        SVAR=SQRT(2.D0*EBATH/RMASS(IAT))*DELT
        DO I=1,3
          CALL GAUSS_RANDOM_NUMBER(RAN)  !SEE PAW_GENERALPURPOSE.F90
          RM(I,IAT) = RM(I,IAT) + SVAR*RAN
          SUM=SUM+0.5D0*RMASS(IAT)*(SVAR*RAN/DELT)**2
        ENDDO
      ENDDO
      IF(TPR) THEN
        SUM=SUM/(1.5D0*DBLE(NAT))
        PRINT*,'KINETIC ENERGY OF THE KICK IN ATOM RANDOMIZATION ',SUM
      END IF
!
!     ==========================================================================
!     == RESTORE OLD SEED OF RANDOM NUMBER GENERATOR                          ==
!     ==========================================================================
      IF(TSETSEED) THEN
        CALL RANDOM_SEED(PUT=OLDSEED)
        DEALLOCATE(OLDSEED)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
!     **************************************************************************
!     **  ESTIMATES THE EFFECTIVE MASS OF THE ELECTRONS BOUND TO              **
!     **  THE ATOMS FROM A MODEL OF ISOLATED ATOMS.                           **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT           ! #(ATOMS)
      REAL(8)   ,INTENT(IN) :: EMASS         ! FICTITIOUS ELECTRON MASS 
      REAL(8)   ,INTENT(IN) :: EMASSCG2      ! G2 ENHANCEMENT FOR FICT.-E.-MASS 
      REAL(8)   ,INTENT(IN) :: PSG2(NAT)     ! 
      REAL(8)   ,INTENT(IN) :: PSG4(NAT)     ! 
      REAL(8)   ,INTENT(OUT):: EFFEMASS(NAT) ! MASS RENORMALIZATION
      INTEGER(4)            :: IAT
!     **************************************************************************
      DO IAT=1,NAT
        EFFEMASS(IAT)=2.D0/3.D0*EMASS*(PSG2(IAT)+EMASSCG2*PSG4(IAT))
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMS_EKIN(NAT,DELTAT,RMASS,RM,RP,EKIN)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATES THE KINETIC ENERGY                                **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT        ! #(ATOMS)
      REAL(8)   ,INTENT(IN) :: DELTAT     ! TIME STEP
      REAL(8)   ,INTENT(IN) :: RMASS(NAT) ! MASS OF THE NUCLEI
      REAL(8)   ,INTENT(IN) :: RP(3,NAT)  ! R(+)
      REAL(8)   ,INTENT(IN) :: RM(3,NAT)  ! R(-)
      REAL(8)   ,INTENT(OUT):: EKIN       ! KINETIC ENERGY OF THE ATOMS
      INTEGER(4)            :: IAT,I
      REAL(8)               :: SVAR
!     ******************************************************************
      EKIN=0.D0
      IF(DELTAT.EQ.0.D0) RETURN
      DO IAT=1,NAT
        SVAR=0.5D0* RMASS(IAT) / (2.D0*DELTAT)**2
        DO I=1,3
          EKIN=EKIN+SVAR*(RP(I,IAT)-RM(I,IAT))**2
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMS_FRICTIONARRAY(NAT,ANNER,ANNEE &
     &                 ,RMASS,EFFEMASS,ANNER1)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT           ! #(ATOMS)
      REAL(8)   ,INTENT(IN) :: ANNER         ! FRICTION ON THE ATOMS
      REAL(8)   ,INTENT(IN) :: ANNEE         ! FRICTION ON THE WAVE FUNCTIONS
      REAL(8)   ,INTENT(IN) :: RMASS(NAT)    ! TRUE MASS OF THE NUCLEUS
      REAL(8)   ,INTENT(IN) :: EFFEMASS(NAT) ! MASS OF TEH ELECTRONS
      REAL(8)   ,INTENT(OUT):: ANNER1(NAT)   ! FRICTION
      REAL(8)               :: RMASS0        ! BARE MASS
      INTEGER(4)            :: IAT
!     ******************************************************************
      DO IAT=1,NAT
!       == CORRECT MASS FOR EFFECTIVE MASS OF THE ELECTRONS ============
        RMASS0=RMASS(IAT)-EFFEMASS(IAT)
!       == CORRECT FRICTION FOR ELECTRON FRICTION  =====================
        ANNER1(IAT)=(RMASS(IAT)*ANNER-EFFEMASS(IAT)*ANNEE)/RMASS0
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMS_PROPAGATE(NAT,DT,RMASS,ANNERVEC,FORCE,R0,RM,RP &
     &                          ,TSTRESS,CELLFRIC)
!     **************************************************************************
!     **                                                                      **
!     **  PROPAGATES THE ATOMIC POSITIONS ACCORDING TO                        **
!     **  THE VERLET ALGORITHM WITH A FRICTION DETERMINED BY ANNER            **
!     **                                                                      **
!     **    R(+) = ( 2*R(0)-(1-ANNER)*R(-)+F*DELT**2/M ) / (1+ANNER)          **
!     **                                                                      **
!     **  FOR ANNER=0 VERLET WITHOUT FRICTION IS OBTAINED,                    **
!     **  FOR ANNER=1 STEEPEST DESCENT                                        **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT           ! #(ATOMS)
      REAL(8)   ,INTENT(IN) :: DT            ! TIME STEP
      REAL(8)   ,INTENT(IN) :: RMASS(NAT)    ! REDUCED MASS OF THE NUCLEUS
      REAL(8)   ,INTENT(IN) :: ANNERVEC(NAT) ! FRICTION ON THE NUCLEI
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)     ! R(0)
      REAL(8)   ,INTENT(IN) :: RM(3,NAT)     ! R(-)
      REAL(8)   ,INTENT(OUT):: RP(3,NAT)     ! R(+)
      REAL(8)   ,INTENT(IN) :: FORCE(3,NAT)  ! FORCE
      LOGICAL(4),INTENT(IN) :: TSTRESS
      REAL(8)   ,INTENT(IN) :: CELLFRIC(3,3) ! FORCE
      INTEGER(4)            :: IAT,I
      REAL(8)               :: SVAR1,SVAR2,SVAR3
      REAL(8)               :: MATP(3,3),MATM(3,3),MATPINV(3,3)
!     **************************************************************************
! THIS (PROPAGATING AS WITH TSTRESS=.FALSE.) IS A CHANGE DONE 200822
! WHILE SEARCHING FOR A BUG, WHICH VIOLATES ENERGY CONSERVATION WHEN
! BOTH, ATOMS AND UNIT CELL, ARE DYNAMIC. ACCORDING TO MY NOTES, THERE
! IS NO DIFFERENCE BETWEEN PROPAGATING THE ATOMS WITH AND WITHOUT
! FRICTION. THE CELL DYNAMICS BECOMES EVIDENT ONLY IN
! ATOMS$SWITCH. 
! PRINT*,'CELLFRIC IN ATOMS_PROPAGATE:',TSTRESS,CELLFRIC
!
!     ==========================================================================
!     == PROPAGATE FOR A CONSTANT UNIT CELL                                   ==
!     ==========================================================================
      IF(.NOT.TSTRESS) THEN
        DO IAT=1,NAT
          SVAR1=2.D0/(1.D0+ANNERVEC(IAT))
          SVAR2=1.D0-SVAR1
          SVAR3=DT**2/RMASS(IAT)/(1.D0+ANNERVEC(IAT))
          RP(:,IAT)=SVAR1*R0(:,IAT)+SVAR2*RM(:,IAT)+SVAR3*FORCE(:,IAT)
        ENDDO
!
!     ==========================================================================
!     == PROPAGATE FOR A DYNAMICAL UNIT CELL                                  ==
!     ==========================================================================
      ELSE 
        DO IAT=1,NAT
          MATP(:,:)=+CELLFRIC(:,:)
          MATM(:,:)=-CELLFRIC(:,:)
          DO I=1,3 
            MATP(I,I)=1.D0+ANNERVEC(IAT)+MATP(I,I)
            MATM(I,I)=1.D0-ANNERVEC(IAT)+MATM(I,I)
          ENDDO
          CALL LIB$INVERTR8(3,MATP,MATPINV)
          RP(:,IAT)=MATMUL(MATPINV,2.D0*R0(:,IAT)-MATMUL(MATM,RM(:,IAT)) &
     &                                  +FORCE(:,IAT)*DT**2/RMASS(IAT))
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMS$WRITE(NFIL,NFILO,TCHK)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE,ONLY : R0,RM,NAT,NAME
      USE RESTART_INTERFACE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN) :: NFIL
      INTEGER(4)            ,INTENT(IN) :: NFILO
      LOGICAL(4)            ,INTENT(OUT):: TCHK
      TYPE (SEPARATOR_TYPE),PARAMETER  :: MYSEPARATOR &
               =SEPARATOR_TYPE(4,'ATOMS','NONE','AUG2003','NONE')
!     TYPE (SEPARATOR_TYPE),PARAMETER  :: MYSEPARATOR &
!              =SEPARATOR_TYPE(3,'ATOMS','NONE','AUG1996','NONE')
      INTEGER(4)                        :: NTASKS,THISTASK
!     ******************************************************************
      TCHK=.FALSE.
!
!     ==================================================================
!     ==  WRITE DATA                                                  ==
!     ==================================================================
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.EQ.1) THEN
        CALL RESTART$WRITESEPARATOR(MYSEPARATOR,NFIL,NFILO,TCHK)
        WRITE(NFIL)NAT
        WRITE(NFIL)NAME(:)
        WRITE(NFIL)R0(:,:)
        WRITE(NFIL)RM(:,:)
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMS$READ(NFIL,NFILO,TCHK)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE RESTART_INTERFACE
      USE ATOMS_MODULE,ONLY : START,R0,RM,NAT,NAME
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN)  :: NFIL
      INTEGER(4)            ,INTENT(IN)  :: NFILO
      LOGICAL(4)            ,INTENT(OUT) :: TCHK
      TYPE (SEPARATOR_TYPE),PARAMETER   :: MYSEPARATOR_OLD1 &
                 =SEPARATOR_TYPE(3,'ATOMS','NONE','AUG1996','NONE')
      TYPE (SEPARATOR_TYPE),PARAMETER   :: MYSEPARATOR &
                 =SEPARATOR_TYPE(4,'ATOMS','NONE','AUG2003','NONE')
      TYPE (SEPARATOR_TYPE)             :: SEPARATOR
      INTEGER(4)                    :: NAT1 
      REAL(8)          ,ALLOCATABLE :: TMP(:,:)   !(3,NAT1)
      CHARACTER(32)    ,ALLOCATABLE :: NAME1(:)   !(NAT1)
      LOGICAL(4)                    :: TCHK1        
      INTEGER(4)                    :: NTASKS,THISTASK
      INTEGER(4)                    :: IAT,IAT1
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      TCHK=.NOT.START
      SEPARATOR=MYSEPARATOR
      IF(THISTASK.EQ.1)CALL RESTART$READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      CALL MPE$BROADCAST('MONOMER',1,TCHK)
      IF(.NOT.TCHK) RETURN
!
!
!     ==================================================================
!     ==  READ DATA                                                   ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        IF(SEPARATOR%VERSION.EQ.MYSEPARATOR%VERSION) THEN
          READ(NFIL)NAT1
          ALLOCATE(TMP(3,NAT1))
          ALLOCATE(NAME1(NAT1))
          READ(NFIL)NAME1(:)
          READ(NFIL)TMP(:,:)
!         == READ R(0)
          DO IAT=1,NAT
            TCHK1=.FALSE.
            DO IAT1=1,NAT1
              IF(NAME1(IAT1).EQ.NAME(IAT)) THEN
                R0(:,IAT)=TMP(:,IAT1)
                TCHK1=.TRUE.
                EXIT
              END IF
            ENDDO
            IF(.NOT.TCHK1) THEN
              CALL ERROR$MSG('ATOM NAMES ON RESTART FILE AND STRUCTURE FILE INCONSISTENT')
              CALL ERROR$CHVAL('NO PARTNER FOUND FOR ATOM ',NAME(IAT))
              CALL ERROR$STOP('ATOMS$READ')
            END IF
          ENDDO
!         == READ R(-)
          READ(NFIL)TMP(:,:)
          DO IAT=1,NAT
            DO IAT1=1,NAT1
              IF(NAME1(IAT1).EQ.NAME(IAT)) THEN
                RM(:,IAT)=TMP(:,IAT1)
                EXIT
              END IF
            ENDDO
          ENDDO
          DEALLOCATE(TMP) 
!       -- LEGACY VERSION -------------------------------------------------
        ELSE IF(SEPARATOR%VERSION.EQ.MYSEPARATOR_OLD1%VERSION) THEN
          READ(NFIL)NAT1
          ALLOCATE(TMP(3,NAT1))
          READ(NFIL)TMP(:,:)
          IF(NAT1.GE.NAT) THEN
            R0(:,:)=TMP(:,1:NAT)
          ELSE
            R0(:,1:NAT1)=TMP(:,:)
          END IF
          READ(NFIL)TMP(:,:)
          IF(NAT1.GE.NAT) THEN
            RM(:,:)=TMP(:,1:NAT)
          ELSE
            RM(:,1:NAT1)=TMP(:,:)
          END IF
          DEALLOCATE(TMP)
!
        ELSE
          CALL ERROR$MSG('VERSION NOT RECOGNIZED')
          CALL ERROR$CHVAL('VERSION',SEPARATOR%VERSION)
          CALL ERROR$STOP('ATOMS$READ')
        END IF
      END IF
!
!     ==================================================================
!     ==  BROADCAST RESULT                                            ==
!     ==================================================================
      CALL MPE$BROADCAST('MONOMER',1,R0)
      CALL MPE$BROADCAST('MONOMER',1,RM)
!
!     ==================================================================
!     ==  STORE DATA                                                  ==
!     ==================================================================
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$INITIALIZE(NAT_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT_
!     ******************************************************************
      NAT=NAT_
      ALLOCATE(NAME(NAT))    ;NAME=' '
      ALLOCATE(R0(3,NAT))    ;R0=0.D0
      ALLOCATE(RM(3,NAT))    ;RM=0.D0
      ALLOCATE(RP(3,NAT))    ;RP=0.D0
      ALLOCATE(FORCE(3,NAT)) ;FORCE=0.D0
      ALLOCATE(RMASS(NAT))   ;RMASS=0.D0
      ALLOCATE(ISPECIES(NAT));ISPECIES=0
      ALLOCATE(PSG2(NAT))    ;PSG2=0.D0
      ALLOCATE(PSG4(NAT))    ;PSG4=0.D0
!      ALLOCATE(EFFEMASS(NAT));EFFEMASS=0.D0
!     ALLOCATE(Z(NAT))       ;Z=0.D0
      ALLOCATE(ZV(NAT))      ;ZV=0.D0
      ALLOCATE(CHARGE(NAT))  ;CHARGE=0.D0
      ALLOCATE(CHARGEANDMOMENTS(4,NAT))  ;CHARGEANDMOMENTS=0.D0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$NATOM(NAT_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: NAT_
!     ******************************************************************
      NAT_=NAT
      RETURN
      END 
!
!     ..................................................................
      SUBROUTINE ATOMLIST$INDEX(IDENT_,IAT_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN) :: IDENT_
      INTEGER(4)   ,INTENT(OUT):: IAT_
      INTEGER(4)               :: IAT
!     ******************************************************************
      DO IAT=1,NAT
        IF(TRIM(IDENT_).EQ.TRIM(NAME(IAT))) THEN
          IAT_=IAT
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('ATOM IDENTIFIER NOT RECOGNIZED')
      CALL ERROR$CHVAL('IDENT_',IDENT_)
      CALL ERROR$STOP('ATOMLIST$INDEX')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$GETR8A(ID_,IAT_,LENG_,VAL_)
!     ******************************************************************
!     **                                                              **
!     **  RETURNS REAL(8) ARRAYS                                      **
!     **  IAT=0 RETURNS THE DATA FOR ALL ATOMS                        **
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      REAL(8)     ,INTENT(OUT):: VAL_(LENG_)
!     ******************************************************************
      IF(IAT_.GT.NAT.OR.IAT_.LT.0) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$GETR8A')
      END IF
!     ==================================================================
!     == R(0)                                                         ==
!     ==================================================================
      IF(ID_.EQ.'R(0)') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.3*NAT) GOTO 9000
          VAL_=RESHAPE(R0,(/3*NAT/))
        ELSE
          IF(LENG_.NE.3) GOTO 9000
          VAL_=R0(:,IAT_)
        END IF
!     ==================================================================
!     == R(-)                                                         ==
!     ==================================================================
      ELSE IF(ID_.EQ.'R(-)') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.3*NAT) GOTO 9000
          VAL_=RESHAPE(RM,(/3*NAT/))
        ELSE
          IF(LENG_.NE.3) GOTO 9000
          VAL_=RM(:,IAT_)
        END IF
!     ==================================================================
!     == R(+)                                                         ==
!     ==================================================================
      ELSE IF(ID_.EQ.'R(+)') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.3*NAT) GOTO 9000
          VAL_=RESHAPE(RP,(/3*NAT/))
        ELSE
          IF(LENG_.NE.3) GOTO 9000
          VAL_=RP(:,IAT_)
        END IF
!     ==================================================================
!     == FORCE                                                        ==
!     ==================================================================
      ELSE IF(ID_.EQ.'FORCE') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.3*NAT) GOTO 9000
          VAL_=RESHAPE(FORCE,(/3*NAT/))
        ELSE
          IF(LENG_.NE.3) GOTO 9000
          VAL_=FORCE(:,IAT_)
        END IF
!     ==================================================================
!     ==  MASS                                                        ==
!     ==================================================================
      ELSE IF(ID_.EQ.'MASS') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          VAL_=RMASS
        ELSE 
           CALL ERROR$MSG('GETR8A DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$GETR8A')
        END IF
!     ==================================================================
!     ==  PS<G2> PSG2                                                 ==
!     ==================================================================
      ELSE IF(ID_.EQ.'PS<G2>') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          VAL_=PSG2
        ELSE 
           CALL ERROR$MSG('GETR8A DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$GETR8A')
        END IF
!     ==================================================================
!     ==  PS<G4> PSG4                                                 ==
!     ==================================================================
      ELSE IF(ID_.EQ.'PS<G4>') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          VAL_=PSG4
        ELSE 
           CALL ERROR$MSG('GETR8A DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$GETR8A')
        END IF
!     ==================================================================
!     ==  ATOMIC NUMBER Z                                             ==
!     ==================================================================
      ELSE IF(ID_.EQ.'Z') THEN
        CALL ERROR$MSG('ID="Z" HAS BEEN DISCONNECTED')
        CALL ERROR$STOP('ATOMLIST$GETR8A')
!!$        IF(IAT_.EQ.0) THEN
!!$          IF(LENG_.NE.NAT) GOTO 9000
!!$          VAL_=Z
!!$        ELSE 
!!$           CALL ERROR$MSG('GETR8A DOES NOT HANDLE SCALARS')
!!$           CALL ERROR$STOP('ATOMLIST$GETR8A')
!!$        END IF
!     ==================================================================
!     ==  NUMBER OF VALENCE ELECTRONS                                  ==
!     ==================================================================
      ELSE IF(ID_.EQ.'ZVALENCE') THEN
        CALL ERROR$MSG('ID ZVALENCE HAS BEEN DISCONNECTED. DO NOT USE')
        CALL ERROR$STOP('ATOMLIST$GETR8A')
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          VAL_=ZV
        ELSE 
           CALL ERROR$MSG('GETR8A DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$GETR8A')
        END IF
!     ==================================================================
!     == POINT CHARGE                                                 ==
!     ==================================================================
      ELSE IF(ID_.EQ.'Q') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          VAL_=CHARGE
        ELSE 
           CALL ERROR$MSG('GETR8A DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$GETR8A')
        END IF
!     ==================================================================
!     == CHARGE AND MAGNETIC MOMENT INTEGRATED OVER THE ASA SPHERE    ==
!     ==================================================================
      ELSE IF(ID_.EQ.'CHARGEANDMOMENTS') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.4*NAT) GOTO 9000
          VAL_=RESHAPE(CHARGEANDMOMENTS,(/4*NAT/))
        ELSE 
          IF(LENG_.NE.4) GOTO 9000
          VAL_=CHARGEANDMOMENTS(:,IAT_)
        END IF
!     ==================================================================
!     == INVALID SELECTION                                            ==
!     ==================================================================
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$GETR8A')
      END IF
      RETURN
 9000 CONTINUE
      CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
      CALL ERROR$CHVAL('ID',ID_)
      CALL ERROR$I4VAL('LENG_',LENG_)
      CALL ERROR$I4VAL('IAT_',IAT_)
      CALL ERROR$I4VAL('NAT',NAT)
      CALL ERROR$STOP('ATOMLIST$GETR8A')
      STOP
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$GETR8(ID_,IAT_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      REAL(8)     ,INTENT(OUT):: VAL_
      INTEGER(4)              :: ISP
!     ******************************************************************
      IF(IAT_.LE.0.OR.IAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$GETR8')
      END IF
      IF(ID_.EQ.'MASS') THEN
        VAL_=RMASS(IAT_)
      ELSE IF(ID_.EQ.'PS<G2>') THEN
        VAL_=PSG2(IAT_)
      ELSE IF(ID_.EQ.'PS<G4>') THEN
        VAL_=PSG4(IAT_)
      ELSE IF(ID_.EQ.'Z') THEN
        CALL ATOMLIST$GETI4('ISPECIES',IAT_,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',VAL_)
        CALL SETUP$UNSELECT()
      ELSE IF(ID_.EQ.'ZVALENCE') THEN
        CALL ERROR$MSG('ID ZVALENCE HAS BEEN DISCONNECTED. DO NOT USE')
        CALL ERROR$STOP('ATOMLIST$GETR8')
!
        CALL ATOMLIST$GETI4('ISPECIES',IAT_,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('ZV',VAL_)
        CALL SETUP$UNSELECT()
      ELSE IF(ID_.EQ.'Q') THEN
        VAL_=CHARGE(IAT_)
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$GETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$GETCHA(ID_,IAT_,LENG_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      CHARACTER(*),INTENT(OUT):: VAL_(LENG_)
!     ******************************************************************
      IF(IAT_.LT.0.OR.IAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$GETCHA')
      END IF
      IF(ID_.EQ.'NAME') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          VAL_(:)=NAME(:)
        ELSE
           CALL ERROR$MSG('GETCHA DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$GETCHA')
        END IF
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$GETCHA')
      END IF
      RETURN
 9000 CONTINUE
      CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
      CALL ERROR$I4VAL('LENG_',LENG_)
      CALL ERROR$I4VAL('NAT',NAT)
      CALL ERROR$STOP('ATOMLIST$GETCHA')
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$GETCH(ID_,IAT_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      CHARACTER(*),INTENT(OUT):: VAL_
!     ******************************************************************
      IF(IAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$GET')
      END IF
      IF(ID_.EQ.'NAME') THEN
        VAL_=NAME(IAT_)
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$GETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$GETI4A(ID_,IAT_,LENG_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      INTEGER(4)  ,INTENT(OUT):: VAL_(LENG_)
!     ******************************************************************
      IF(IAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$GET')
      END IF
      IF(ID_.EQ.'ISPECIES') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) THEN
            CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
            CALL ERROR$I4VAL('LENG_',LENG_)
            CALL ERROR$I4VAL('NAT',NAT)
            CALL ERROR$STOP('ATOMLIST$GETI4A')
          END IF
          VAL_=ISPECIES
        ELSE
           CALL ERROR$MSG('GETI4A DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$GETI4A')
        END IF
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$GETI4A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$GETI4(ID_,IAT_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      INTEGER(4)  ,INTENT(OUT):: VAL_
!     ******************************************************************
      IF(IAT_.LE.0.OR.IAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$GETI4')
      END IF
      IF(ID_.EQ.'ISPECIES') THEN
        VAL_=ISPECIES(IAT_)
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$GETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$SETR8A(ID_,IAT_,LENG_,VAL_)
!     ******************************************************************
!     **                                                              **
!     **  RETURNS REAL(8) ARRAYS                                      **
!     **  IAT=0 RETURNS THE DATA FOR ALL ATOMS                        **
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      REAL(8)     ,INTENT(IN) :: VAL_(LENG_)
!     ******************************************************************
      IF(IAT_.GT.NAT.OR.IAT_.LT.0) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$SETR8A')
      END IF
!     ==================================================================
!     == R(0)                                                         ==
!     ==================================================================
      IF(ID_.EQ.'R(0)') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.3*NAT) GOTO 9000
          R0(:,:)=RESHAPE(VAL_,(/3,NAT/))
        ELSE
          IF(LENG_.NE.3) GOTO 9000
          R0(:,IAT_)=VAL_(:)
        END IF
!     ==================================================================
!     == R(-)                                                         ==
!     ==================================================================
      ELSE IF(ID_.EQ.'R(-)') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.3*NAT) GOTO 9000
          RM(:,:)=RESHAPE(VAL_,(/3,NAT/))
        ELSE
          IF(LENG_.NE.3) GOTO 9000
          RM(:,IAT_)=VAL_(:)
        END IF
!     ==================================================================
!     == R(+)                                                         ==
!     ==================================================================
      ELSE IF(ID_.EQ.'R(+)') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.3*NAT) GOTO 9000
          RP(:,:)=RESHAPE(VAL_,(/3,NAT/))
        ELSE
          IF(LENG_.NE.3) GOTO 9000
          RP(:,IAT_)=VAL_(:)
        END IF
!     ==================================================================
!     == FORCE                                                        ==
!     ==================================================================
      ELSE IF(ID_.EQ.'FORCE') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.3*NAT) GOTO 9000
          FORCE(:,:)=RESHAPE(VAL_,(/3,NAT/))
        ELSE
          IF(LENG_.NE.3) GOTO 9000
          FORCE(:,IAT_)=VAL_(:)
        END IF
!     ==================================================================
!     ==  MASS                                                        ==
!     ==================================================================
      ELSE IF(ID_.EQ.'MASS') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          RMASS=VAL_(:)
        ELSE 
          CALL ERROR$MSG('SETR8A DOES NOT HANDLE SCALARS')
          CALL ERROR$STOP('ATOMLIST$SETR8A')
        END IF
!     ==================================================================
!     ==  PSG2 'PS<G2>'                                               ==
!     ==================================================================
      ELSE IF(ID_.EQ.'PS<G2>') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          PSG2=VAL_(:)
        ELSE 
          CALL ERROR$MSG('SETR8A DOES NOT HANDLE SCALARS')
          CALL ERROR$STOP('ATOMLIST$SETR8A')
        END IF
!     ==================================================================
!     ==  PSG4 'PS<G4>'                                               ==
!     ==================================================================
      ELSE IF(ID_.EQ.'PS<G4>') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          PSG4=VAL_(:)
        ELSE 
          CALL ERROR$MSG('SETR8A DOES NOT HANDLE SCALARS')
          CALL ERROR$STOP('ATOMLIST$SETR8A')
        END IF
!     ==================================================================
!     ==  ATOMIC NUMBER Z                                             ==
!     ==================================================================
      ELSE IF(ID_.EQ.'Z') THEN 
        CALL ERROR$MSG('IF="Z" HAS BEEN DISCONNECTED')
        CALL ERROR$STOP('ATOMLIST$SETR8A')
!!$        IF(IAT_.EQ.0) THEN
!!$          IF(LENG_.NE.NAT) GOTO 9000
!!$          Z(:)=VAL_(:)
!!$        ELSE 
!!$           CALL ERROR$MSG('SETR8A DOES NOT HANDLE SCALARS')
!!$           CALL ERROR$STOP('ATOMLIST$SETR8A')
!!$        END IF
!     ==================================================================
!     ==  NUMBER OF VALENCE ELECTRONS                                  ==
!     ==================================================================
      ELSE IF(ID_.EQ.'ZVALENCE') THEN
        CALL ERROR$MSG('ID ZVALENCE HAS BEEN DISCONNECTED. DO NOT USE')
        CALL ERROR$STOP('ATOMLIST$SETR8A')
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          ZV(:)=VAL_(:)
        ELSE 
          CALL ERROR$MSG('SETR8A DOES NOT HANDLE SCALARS')
          CALL ERROR$STOP('ATOMLIST$SETR8A')
        END IF
!     ==================================================================
!     == POINT CHARGE                                                 ==
!     ==================================================================
      ELSE IF(ID_.EQ.'Q') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          CHARGE(:)=VAL_(:)
        ELSE 
           CALL ERROR$MSG('SETR8A DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$SETR8A')
        END IF
!     ==================================================================
!     == POINT CHARGE                                                 ==
!     ==================================================================
      ELSE IF(ID_.EQ.'CHARGEANDMOMENTS') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.4*NAT) GOTO 9000
          CHARGEANDMOMENTS(:,:)=RESHAPE(VAL_,(/4,NAT/))
        ELSE 
          IF(LENG_.NE.4) GOTO 9000
          CHARGEANDMOMENTS(:,IAT_)=VAL_(:)
        END IF
!     ==================================================================
!     == INVALID SELECTION                                            ==
!     ==================================================================
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$SETR8A')
      END IF
      RETURN
 9000 CONTINUE
      CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
      CALL ERROR$I4VAL('LENG_',LENG_)
      CALL ERROR$I4VAL('NAT',NAT)
      CALL ERROR$STOP('ATOMLIST$SETR8A')
      STOP
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$SETR8(ID_,IAT_,VAL_)
!     *****************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      REAL(8)     ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(IAT_.LE.0.OR.IAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$SETR8')
      END IF
      IF(ID_.EQ.'MASS') THEN
        RMASS(IAT_)=VAL_
      ELSE IF(ID_.EQ.'PS<G2>') THEN
        PSG2(IAT_)=VAL_
      ELSE IF(ID_.EQ.'PS<G4>') THEN
        PSG4(IAT_)=VAL_
      ELSE IF(ID_.EQ.'Q') THEN
        CHARGE(IAT_)=VAL_
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$SETCHA(ID_,IAT_,LENG_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      CHARACTER(*),INTENT(IN) :: VAL_(LENG_)
!     ******************************************************************
      IF(IAT_.LT.0.OR.IAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$SETCHA')
      END IF
      IF(ID_.EQ.'NAME') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          NAME(:)=VAL_(:)
        ELSE
           CALL ERROR$MSG('SETI4A DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$SETCHA')
        END IF
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$SETCHA')
      END IF
      RETURN
 9000 CONTINUE
      CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
      CALL ERROR$I4VAL('LENG_',LENG_)
      CALL ERROR$I4VAL('NAT',NAT)
      CALL ERROR$STOP('ATOMLIST$SETCHA')
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$SETCH(ID_,IAT_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      CHARACTER(*),INTENT(IN) :: VAL_
!     ******************************************************************
      IF(IAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$SETCH')
      END IF
      IF(ID_.EQ.'NAME') THEN
        NAME(IAT_)=VAL_
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$SETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$SETI4A(ID_,IAT_,LENG_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      INTEGER(4)  ,INTENT(IN) :: VAL_(LENG_)
!     ******************************************************************
      IF(IAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$SETI4A')
      END IF
      IF(ID_.EQ.'ISPECIES') THEN
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          ISPECIES=VAL_
        ELSE
           CALL ERROR$MSG('SETI4A DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$SETR8A')
        END IF
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$SETI4A')
      END IF
      RETURN
 9000 CONTINUE
      CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
      CALL ERROR$I4VAL('LENG_',LENG_)
      CALL ERROR$I4VAL('NAT',NAT)
      CALL ERROR$STOP('ATOMLIST$SETI4A')
      STOP
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$SETI4(ID_,IAT_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      INTEGER(4)  ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(IAT_.LE.0.OR.IAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$SETI4')
      END IF
      IF(ID_.EQ.'ISPECIES') THEN
        ISPECIES(IAT_)=VAL_
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$SETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$GET(ID_,NBYTE_,IAT_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: NBYTE_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      CHARACTER(1),INTENT(OUT):: VAL_(NBYTE_)
!     ******************************************************************
      IF(IAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$GET')
      END IF
      IF(ID_.EQ.'NAME') THEN
        VAL_=TRANSFER(NAME(IAT_),VAL_,32)
      ELSE IF(ID_.EQ.'ISPECIES') THEN
        VAL_=TRANSFER(ISPECIES(IAT_),VAL_,4)
      ELSE IF(ID_.EQ.'R(0)') THEN
        VAL_=TRANSFER(R0(:,IAT_),VAL_,8*3)
      ELSE IF(ID_.EQ.'R(-)') THEN
        VAL_=TRANSFER(RM(:,IAT_),VAL_,8*3)
      ELSE IF(ID_.EQ.'R(+)') THEN
        VAL_=TRANSFER(RP(:,IAT_),VAL_,8*3)
      ELSE IF(ID_.EQ.'FORCE') THEN
        VAL_=TRANSFER(FORCE(:,IAT_),VAL_,8*3)
      ELSE IF(ID_.EQ.'MASS') THEN
        VAL_=TRANSFER(RMASS(IAT_),VAL_,8)
      ELSE IF(ID_.EQ.'PS<G2>') THEN
        VAL_=TRANSFER(PSG2(IAT_),VAL_,8)
      ELSE IF(ID_.EQ.'PS<G4>') THEN
        VAL_=TRANSFER(PSG4(IAT_),VAL_,8)
!!$      ELSE IF(ID_.EQ.'Z') THEN
!!$        VAL_=TRANSFER(Z(IAT_),VAL_,8)
      ELSE IF(ID_.EQ.'ZVALENCE') THEN
        CALL ERROR$MSG('ID ZVALENCE HAS BEEN DISCONNECTED. DO NOT USE')
        CALL ERROR$STOP('ATOMLIST$GET')
        VAL_=TRANSFER(ZV(IAT_),VAL_,8)
      ELSE IF(ID_.EQ.'POINTCHARGE') THEN
!       == KEYWORD SHALL BE REPLACED BY 'Q'
        VAL_=TRANSFER(CHARGE(IAT_),VAL_,8)
      ELSE IF(ID_.EQ.'Q') THEN
        VAL_=TRANSFER(CHARGE(IAT_),VAL_,8)
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$GET')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$ORDER
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      INTEGER(4)    :: ISP,IAT,IAT1,ISP1
      CHARACTER(32) :: XNAME         
      REAL(8)       :: XR0(3)        
      REAL(8)       :: XRM(3)        
      REAL(8)       :: XRP(3)        
      REAL(8)       :: XFORCE(3)     
      REAL(8)       :: XRMASS        
      REAL(8)       :: XPSG2       
      REAL(8)       :: XPSG4       
      INTEGER(4)    :: XISPECIES     
      REAL(8)       :: XZV           
      REAL(8)       :: XCHARGE       
!     ******************************************************************
      ISP=1
      IAT=1
      DO WHILE(IAT.LT.NAT) 
        DO IAT1=IAT,NAT
          ISP1=ISPECIES(IAT1)
          IF(ISP1.LT.1) THEN
            CALL ERROR$MSG('ISPECIES MUST BE SET')
            CALL ERROR$STOP('ATOMLIST$ORDER')
          END IF
          IF(ISP1.EQ.ISP) THEN
            IF(IAT1.NE.IAT) THEN
              XNAME        =NAME(IAT)
              XR0(:)       =R0(:,IAT)    
              XRM(:)       =RM(:,IAT)    
              XRP(:)       =RP(:,IAT)    
              XFORCE(:)    =FORCE(:,IAT) 
              XRMASS       =RMASS(IAT)   
              XPSG2        =PSG2(IAT)  
              XPSG4        =PSG4(IAT)  
              XISPECIES    =ISPECIES(IAT)
              XZV          =ZV(IAT)      
              XCHARGE      =CHARGE(IAT)     
!             =========
              NAME(IAT)        =NAME(IAT1)       
              R0(:,IAT)        =R0(:,IAT1)       
              RM(:,IAT)        =RM(:,IAT1)       
              RP(:,IAT)        =RP(:,IAT1)       
              FORCE(:,IAT)     =FORCE(:,IAT1)    
              RMASS(IAT)       =RMASS(IAT1)      
              PSG2(IAT)        =PSG2(IAT1)     
              PSG4(IAT)        =PSG4(IAT1)     
              ISPECIES(IAT)    =ISPECIES(IAT1)   
              ZV(IAT)          =ZV(IAT1)         
              CHARGE(IAT)      =CHARGE(IAT1)     
!             =========
              NAME(IAT1)        =XNAME        
              R0(:,IAT1)        =XR0(:)       
              RM(:,IAT1)        =XRM(:)       
              RP(:,IAT1)        =XRP(:)       
              FORCE(:,IAT1)     =XFORCE(:)    
              RMASS(IAT1)       =XRMASS       
              PSG2(IAT1)        =XPSG2
              PSG4(IAT1)        =XPSG4
              ISPECIES(IAT1)    =XISPECIES    
              ZV(IAT1)          =XZV          
              CHARGE(IAT1)      =XCHARGE      
!             =========
            END IF
            IAT=IAT+1
            GOTO 100
          END IF
        ENDDO
        ISP=ISP+1
 100    CONTINUE
      ENDDO
      RETURN
      END      
!
!     ..................................................................
      SUBROUTINE ATOMLIST$INDEXX(NAME,IAT,IT)
!     ******************************************************************
!     **  ATOMLIST$INDEXX                                             **
!     **  RESOLVES AN ATOM NAME IN THE EXTENDED NOTATION:             **
!     **  NAME:+I+J+K WHERE NAME IS THE ATOM NAME AND +I,+J,+K ARE THE**
!     **  TRANSLATIONS INTO THE THREE LATTICE DIRECTIONS              **
!     **  THE PLUS SIGN IS OPTIONAL AND CAN BE A MINUS SIGN           **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(OUT):: IAT
      INTEGER(4)  ,INTENT(OUT):: IT(3)
      CHARACTER(32)           :: NAME1
      INTEGER(4)              :: ICH
      INTEGER(4)              :: ISIGN
      INTEGER(4)              :: I,IEND,IPOS
!     ******************************************************************
      NAME1=NAME
      IPOS=INDEX(NAME1,':',BACK=.TRUE.)
      IT=0
      IF(IPOS.GT.1) THEN   ! ':' FOUND TEST FOR EXTENDED STRUCTURE
        NAME1=NAME(1:IPOS-1)
        IEND=LEN_TRIM(NAME)
        IPOS=IPOS+1
        DO I=1,3
!         = EXTRACT SIGN
          IF(IPOS.GT.IEND) THEN
            IPOS=-1
            EXIT
          END IF
          ICH=IACHAR(NAME(IPOS:IPOS))
          IF(ICH.EQ.45) THEN        ! '-' SIGN FOUND
            ISIGN=-1
            IPOS=IPOS+1
          ELSE IF(ICH.EQ.43) THEN   ! '+' SIGN FOUND
            ISIGN=+1
            IPOS=IPOS+1
          ELSE                      ! NO SIGN ENCOUNTERED
            ISIGN=+1
          END IF
!       
          IF(IPOS.GT.IEND) THEN
            EXIT
          END IF
          ICH=IACHAR(NAME(IPOS:IPOS))-48
          IF(ICH.GE.0.AND.ICH.LE.9) THEN
            IT(I)=ICH*ISIGN
            IPOS=IPOS+1
          ELSE
            IPOS=-1
            EXIT
          END IF
        ENDDO
!       ==  IF 
        IF(IPOS.NE.IEND+1) THEN
          NAME1=NAME
          IT(:)=0
        END IF
      END IF
      CALL ATOMLIST$INDEX(NAME1,IAT)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMS$DUMMYCHECK(STRING)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: STRING
      REAL(8)               :: VEL1(3,NAT)
      REAL(8)               :: VEL2(3,NAT)
      REAL(8)               :: SUM1,SUM2
      INTEGER(4)            :: IAT
      INTEGER(4)            :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK) 
      IF(THISTASK.NE.1) RETURN
      VEL1(:,:)=R0(:,:)-RM(:,:)
      VEL2(:,:)=RP(:,:)-R0(:,:)
      SUM1=0.D0
      SUM2=0.D0
      DO IAT=1,NAT
        SUM1=SUM1+DOT_PRODUCT(VEL1(:,IAT),VEL1(:,IAT))
        SUM2=SUM2+DOT_PRODUCT(VEL2(:,IAT),VEL2(:,IAT))
      ENDDO
      WRITE(*,FMT='("CHECK ",2E10.2,A)')SUM1,SUM2,TRIM(STRING)
      RETURN
      END
