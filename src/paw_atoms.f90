!
!....................................................................... 
MODULE ATOMS_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: ATOMS                                                      **
!**                                                                   **
!**  PURPOSE: MAINTAINS THE ATOMIC POSITIONS                          **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    ATOMS$SETR8(ID,VAL)                                            **
!**    ATOMS$REPORT                                                   **
!**    ATOMS$PROPAGATE                                                **
!**    ATOMS$INITIALIZE                                               **
!**                                                                   **
!**  REMARKS:                                                         **
!**    FOR HISTORICAL REASONS, PART OF THIS OBJECT IS NAMED ATOMLIST  **
!**                                                                   **
!******************PETER E. BLOECHL, IBM RESEARCH LABORATORY (1996)*****
REAL(8)          :: DELT=0.D0      ! TIME STEP
REAL(8)          :: AMPRE=0.D0     ! TARGET TEMPERATURE FOR RANDOMIZATION
REAL(8)          :: ANNER=0.D0     ! FRICTION
LOGICAL(4)       :: TSTOP=.FALSE.  ! ZERO INITIAL VELOCITIES
LOGICAL(4)       :: TDYN=.TRUE.    ! ATOMIC MOTION IS SWITCHED OFF
LOGICAL(4)       :: TRANDOMIZE=.FALSE.  ! RANDOMIZE INITIAL VELOCITIES
LOGICAL(4)       :: START=.FALSE.  ! RANDOMIZE INITIAL VELOCITIES
REAL(8)          :: ANNEE=0.D0
LOGICAL(4)       :: TCONSTRAINTREFERENCE=.FALSE.
INTEGER(4)       :: NAT=0
CHARACTER(32),ALLOCATABLE :: NAME(:)
REAL(8)      ,ALLOCATABLE :: R0(:,:)
REAL(8)      ,ALLOCATABLE :: RM(:,:)
REAL(8)      ,ALLOCATABLE :: RP(:,:)
REAL(8)      ,ALLOCATABLE :: FORCE(:,:)
REAL(8)      ,ALLOCATABLE :: RMASS(:)
REAL(8)      ,ALLOCATABLE :: PSG2(:)
REAL(8)      ,ALLOCATABLE :: PSG4(:)
!REAL(8)     ,ALLOCATABLE :: REDUCEDMASS(:)
INTEGER(4)   ,ALLOCATABLE :: ISPECIES(:)
REAL(8)      ,ALLOCATABLE :: Z(:)
REAL(8)      ,ALLOCATABLE :: ZV(:)
REAL(8)      ,ALLOCATABLE :: CHARGE(:)
REAL(8)                   :: CELLKIN(3,3) ! SUM M*RDOT_I*RDOT_J
END MODULE ATOMS_MODULE
!
!     .................................................................. 
      SUBROUTINE ATOMS$SETR8(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'TIMESTEP') THEN
        DELT=VAL_
      ELSE IF(ID_.EQ.'AMPRE') THEN
        AMPRE=VAL_
      ELSE IF(ID_.EQ.'FRICTION') THEN
        ANNER=VAL_
      ELSE IF(ID_.EQ.'ANNEE') THEN
        ANNEE=VAL_
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
      CALL MPE$QUERY(NTASKS,THISTASK) 
      IF(THISTASK.NE.1) RETURN
      CALL CONSTANTS('KB',KELVIN)
      CALL REPORT$TITLE(NFIL,'ATOMS')
      CALL REPORT$I4VAL(NFIL,'NUMBER OF ATOMS',NAT,' ')
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
      CHARACTER(LEN=82)     :: STRING
      CHARACTER(LEN=32)     :: IDENT1
      REAL(8)               :: X(5)
      REAL(8)               :: U      ! MASS UNIT C12/12
      INTEGER(4)            :: NTASKS,THISTASK
      REAL(8)               :: EFFEMASS(NAT)
      REAL(8)               :: EMASS,EMASSCG2
      REAL(8)               :: RBAS(3,3)
!     ******************************************************************
      CALL MPE$QUERY(NTASKS,THISTASK) 
      IF(THISTASK.NE.1) RETURN
      CALL CONSTANTS('U',U)
      CALL WAVES$GETR8('EMASS',EMASS)
      CALL WAVES$GETR8('EMASSCG2',EMASSCG2)
      CALL ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
      CALL CELL$GETR8A('T(0)',9,RBAS)
!
      CALL REPORT$TITLE(NFIL,'ATOMLIST REPORT')
      WRITE(NFIL,FMT='("T1=",3F10.6/"T2=",3F10.6/"T3=",3F10.6)')RBAS
      WRITE(STRING,FMT='(T1,A,T10,A,T43,A,T53,A,T65,A)') &
     &                  'NAME','R(0)','M[U]','MPSI_EFF[U]','Q[E]'
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
        WRITE(STRING(10:42),FMT='("(",F9.5,",",F9.5,",",F9.5,")")')R0(:,IAT)
        WRITE(STRING(43:52),FMT='(F8.4)')RMASS(IAT)/U
        WRITE(STRING(53:64),FMT='(F8.4)')EFFEMASS(IAT)/U
        WRITE(STRING(63:72),FMT='(F8.5)')-CHARGE(IAT)
        WRITE(NFIL,FMT='(A)')TRIM(STRING)
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMS$INITIALIZE
!     ******************************************************************
!     **  READ ATOMIC POSITIONS FROM ATOMLIST                         **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      INTEGER(4)           :: IAT,I
      REAL(8)              :: EKIN
      REAL(8)              :: EFFEMASS(NAT)
      REAL(8)              :: EMASS,EMASSCG2
      REAL(8)              :: SVAR
      LOGICAL(4)           :: TERR,TCHK
      INTEGER(4)           :: NFREE
!     ******************************************************************
                              CALL TRACE$PUSH('ATOMS$INITIALIZE')
!     CALL SHADOW$INITIALIZE
!     CALL SHADOW$OPTIMIZE(100,TOL,TCHK)
      DO IAT=1,NAT
        DO I=1,3
          FORCE(I,IAT)=0.D0
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  FIX ATOM THERMOSTAT                                         ==
!     ==================================================================
      CALL ATOMLIST$NATOM(NAT)
      CALL CONSTRAINTS$NFREE(NAT,NFREE)
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK)CALL THERMOSTAT$SCALEGFREE(REAL(NFREE,KIND=8))
!
!     ==================================================================
!     ==   INITIALIZE THERMOSTATS FOR THE WAVE FUNCTIONS              ==
!     ==================================================================
      CALL WAVES$GETR8('EMASS',EMASS)
      CALL WAVES$GETR8('EMASSCG2',EMASSCG2)
      CALL ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
!
!     ==================================================================
!     ==   CHECK IF EFFECTIVE MASS IS NEGATIVE                        ==
!     ==================================================================
      TERR=.FALSE.
      DO IAT=1,NAT
        IF(EFFEMASS(IAT).GT.RMASS(IAT)) THEN
          IF(.NOT.TERR) THEN
            CALL ERROR$MSG('REDUCED MASS MUST NOT BE NEGATIVE')
            TERR=.TRUE.
          END IF
          CALL ERROR$CHVAL('ATOM ',NAME(IAT))
          CALL ERROR$R8VAL('MASS ',RMASS(IAT))
          CALL ERROR$R8VAL('REDUCED MASS ',RMASS(IAT)-EFFEMASS(IAT))
        END IF
      ENDDO
      IF(TERR) THEN
        CALL ERROR$MSG('SOLUTION: REDUCE WAVE FUNCTION MASS')
        CALL ERROR$MSG('OR SET PS<G2> AND PS<G4> TO ZERO')
        CALL ERROR$STOP('ATOMS$INITIALIZE')
      END IF
!
!     ==================================================================
!     ==   ADJUST KINETIC ENERGY OF THE WAVE FUNCTION THERMOSTAT      ==
!     ==================================================================
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
          SVAR=SVAR/REAL(NAT)*EKIN
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
          CALL ERROR$MSG('WAVE FUNCTION THERMOSTAT CANNOT BE USED')
          CALL ERROR$MSG('WITHOUT ATOM THERMOSTAT')
          CALL ERROR$STOP('ATOMS$INITIALIZE')
        END IF
      END IF
                               CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMS$PROPAGATE()
!     ******************************************************************
!     **  PROPAGATE: PROPAGATES ATOMIC POSITIONS                      **
!     ******************************************************************
      USE ATOMS_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)     :: IAT,I
      REAL(8)        :: EKIN
      REAL(8)        :: DEKIN
      REAL(8)        :: EMASS,EMASSCG2
      REAL(8)        :: ENOSEP
      REAL(8)        :: EFFEMASS(NAT)
      REAL(8)        :: REDRMASS(NAT)
      LOGICAL(4)     :: TCHK
      INTEGER(4)     :: ITER
      LOGICAL(4)     :: TSTRESS
      REAL(8)        :: STRESS(3,3)
      REAL(8)        :: CELLFRIC(3,3)
      REAL(8)        :: RBAS(3,3)
!     ******************************************************************
                              CALL TRACE$PUSH('ATOMS$PROPAGATE')
DO IAT=1,NAT
  WRITE(*,FMT='("FORCE ",I3,3F15.10)')IAT,FORCE(:,IAT)
ENDDO
      IF(.NOT.TDYN) THEN
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
        CALL MPE$BROADCAST(1,RM)
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
!
!     ==================================================================
!     ==  DISTRIBUTE FORCES FOR LINK ATOMS OF QM-MM COUPLING         ==
!     ==================================================================
      CALL QMMM$SCALEFORCE(NAT,FORCE)
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
      CALL ATOMS_PROPAGATE(NAT,DELT,ANNER,ANNEE,RMASS,EFFEMASS &
   &                      ,FORCE,R0,RM,RP,TSTRESS,CELLFRIC,CELLKIN)
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
      END
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
      INTEGER(4)     :: IAT,I,J
      LOGICAL(4)     :: TSTRESS
      REAL(8)        :: RP1(3,NAT)
      REAL(8)        :: STRESS(3,3)
      REAL(8)        :: STRESS1(3,3)
      REAL(8)        :: V(3)
      REAL(8)        :: RBAS(3,3)
!     ******************************************************************
                              CALL TRACE$PUSH('ATOMS$CONSTRAINTS')
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
PRINT*,'WARNING FROM ATOMS$CONSTRAINTS: TSTRESS SET TO FALSE'
TSTRESS=.FALSE.
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
        CALL lib$invertr8(3,MAPTOCELL,MAPTOCELL)
!CALL INVERT(3,MAPTOCELL,MAPTOCELL)
        DO IAT=1,NAT
          RP(:,IAT)=MATMUL(MAPTOCELL,RP1(:,IAT))
        ENDDO
      ELSE
        DO IAT=1,NAT
          RP(:,IAT)=RP1(:,IAT)
        ENDDO
      END IF
!
!     ==================================================================
!     == ADD KINETIC FRICTION TO STRESS                               ==
!     ==================================================================
      STRESS1=0.D0
      DO IAT=1,NAT
        DO I=1,3
          V(I)=(RP(I,IAT)-RM(I,IAT))/(2.D0*DELT)
        ENDDO
        DO I=1,3
          DO J=1,3
            STRESS1(I,J)=STRESS1(I,J)+REDRMASS(IAT)*V(I)*V(J)
          ENDDO
        ENDDO
      ENDDO
      CALL CELL$SETR8A('KINSTRESS',9,STRESS1)
WRITE(*,FMT='("KIN-STRESS ",3F10.5)')STRESS1(1,:)
WRITE(*,FMT='("KIN-STRESS ",3F10.5)')STRESS1(2,:)
WRITE(*,FMT='("KIN-STRESS ",3F10.5)')STRESS1(3,:)
                            CALL TRACE$POP
      RETURN
      END
!
!     .................................................................. 
      SUBROUTINE ATOMS$EKIN(EKIN_)
!     ******************************************************************
!     **  CALCULATE KINETIC ENERGY                                    **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: EKIN_
      REAL(8)             :: DEKIN
!     ******************************************************************
      CALL ATOMS_EKIN(NAT,DELT,RMASS,RM,RP,EKIN_)
      CALL QMMM$DEKIN(DELT,DEKIN)   !CORRECTION FOR LINK ATOMS
      EKIN_=EKIN_+DEKIN
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
      INTEGER(4)          :: NFREE
!     ******************************************************************
      CALL CONSTRAINTS$NFREE(NAT,NFREE)
      IF(NFREE.EQ.0) THEN
        TEMPERATURE_=0.D0
      ELSE
        CALL ATOMS_EKIN(NAT,DELT,RMASS,RM,RP,EKIN)
        TEMPERATURE_=EKIN/(0.5D0*DBLE(NFREE))
      END IF
      RETURN
      END
!
!     .................................................................. 
      SUBROUTINE ATOMS$SWITCH()
!     ******************************************************************
!     **  SWITCH ATOMIC POSITIONS  R(+)->R(0)->R(M)                   **
!     ****************************************************************** 
      USE ATOMS_MODULE
      IMPLICIT NONE
      REAL(8)     :: MAPTOCELL(3,3)
      LOGICAL(4)  :: TSTRESS
      INTEGER(4)  :: IAT
      REAL(8)     :: DRPM(3)
!     ******************************************************************
      FORCE(:,:)=0.D0
      CALL QMMM$SWITCH
      CALL CONTINUUM$SWITCH(NAT,RP,R0)
!
      CALL CELL$GETL4('MOVE',TSTRESS)
      IF(TSTRESS) THEN
        CALL CELL$GETR8A('MAPTOCELL',9,MAPTOCELL)
        DO IAT=1,NAT
          DRPM(:)=RP(:,IAT)-RM(:,IAT)
          R0(:,IAT)=MATMUL(MAPTOCELL,R0(:,IAT))
          RM(:,IAT)=MATMUL(MAPTOCELL,RM(:,IAT))
          RP(:,IAT)=RM(:,IAT)+DRPM(:)
        ENDDO
      END IF
      IF(.NOT.TDYN) RETURN
      RM(:,:)=R0(:,:)
      R0(:,:)=RP(:,:)
      RP(:,:)=2.D0*R0(:,:)-RM(:,:)
      CALL CONSTRAINTS$SWITCH()
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMS_RANDOMIZEVELOCITY(NAT,RMASS,RM,EBATH,DELT)
!     ******************************************************************
!     **                                                              **
!     **  RANDOMIZED VELOCITY ACCORDING TO A TEMPERATURE              **
!     **                                                              **
!     **  REMARK: ANY CONSTRAINTS MUST BE ENFORCED SUBSEQUENTLY       **
!     **                                                              **
!     **    EBATH    K_B*T WHERE THE T IS THE TEMPERATURE OF THE      **
!     **             HEATBATH                                         **
!     **                                                              **
!     **    THE TARGET KINETIC ENERGY IS 1.5*NAT*EBATH                **
!     **                                                              **
!     ******************************************************************
      IMPLICIT none
      INTEGER(4),INTENT(IN)   :: NAT        ! #(ATOMS)
      REAL(8)   ,INTENT(IN)   :: DELT       ! TIME STEP
      REAL(8)   ,INTENT(IN)   :: EBATH      ! K_BT OF THE HEAT BATH
      REAL(8)   ,INTENT(IN)   :: RMASS(NAT) ! ATOMIC MASSES
      REAL(8)   ,INTENT(INOUT):: RM(3,NAT)  ! R(-)
      LOGICAL(4),PARAMETER    :: TPR=.FALSE.
      INTEGER(4)              :: IAT,I
      REAL(8)                 :: SVAR,RAN,SUM
!     ******************************************************************
      SUM=0.D0
      DO IAT=1,NAT
        SVAR=DSQRT(2.D0*EBATH/RMASS(IAT))*DELT
        DO I=1,3
          CALL GAUSS_RANDOM_NUMBER(RAN)
          RM(I,IAT) = RM(I,IAT) + SVAR*RAN
          SUM=SUM+0.5D0*RMASS(IAT)*(SVAR*RAN/DELT)**2
        ENDDO
      ENDDO
      IF(TPR) THEN
        SUM=SUM/(1.5D0*DBLE(NAT))
        PRINT*,'KINETIC ENERGY OF THE KICK IN ATOM RANDOMIZATION ',SUM
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMS_EFFEMASS(NAT,EMASS,EMASSCG2,PSG2,PSG4,EFFEMASS)
!     ******************************************************************
!     **                                                              **
!     **  ESTIMATES THE EFFECTIVE MASS OF THE ELECTRONS BOUND TO      **
!     **  THE ATOMS FROM A MODEL OF ISOLATED ATOMS.                   **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT           ! #(ATOMS)
      REAL(8)   ,INTENT(IN) :: EMASS         ! FICTITIOUS ELECTRON MASS 
      REAL(8)   ,INTENT(IN) :: EMASSCG2      ! G2 ENHANCEMENT FOR FICT.-E.-MASS 
      REAL(8)   ,INTENT(IN) :: PSG2(NAT)     ! 
      REAL(8)   ,INTENT(IN) :: PSG4(NAT)     ! 
      REAL(8)   ,INTENT(OUT):: EFFEMASS(NAT) ! MASS RENORMALIZATION
      REAL(8)               :: SVAR
      INTEGER(4)            :: IAT
!     ******************************************************************
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
      SUBROUTINE ATOMS_PROPAGATE(NAT,DT,ANNER,ANNEE &
     &                 ,RMASS,EFFEMASS,FORCE,R0,RM,RP &
     &                 ,TSTRESS,CELLFRIC,CELLKIN)
!     ******************************************************************
!     **                                                              **
!     **  PROPAGATES THE ATOMIC POSITIONS ACCORDING TO                **
!     **  THE VERLET ALGORITHM WITH A FRICTION DETERMINED BY ANNER    **
!     **                                                              **
!     **    R(+) = ( 2*R(0)-(1-ANNER)*R(-)+F*DELT**2/M ) / (1+ANNER)  **
!     **                                                              **
!     **  FOR ANNER=0 VERLET WITHOUT FRICTION IS OBTAINED,            **
!     **  FOR ANNER=1 STEEPEST DESCENT                                **
!     **                                                              **
!     ******************************************************************
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
      REAL(8)   ,INTENT(IN) :: CELLFRIC(3,3) ! FORCE
      REAL(8)   ,INTENT(OUT):: CELLKIN(3,3) ! FORCE
      INTEGER(4)            :: IAT,I,J
      REAL(8)               :: ANNER1
      REAL(8)               :: RMASS0        ! BARE ATOM MASS
      REAL(8)               :: SVAR1,SVAR2,SVAR3
      REAL(8)               :: MATP(3,3),MATM(3,3),MATPINV(3,3)
      REAL(8)               :: V(3)
!     ******************************************************************
      CELLKIN=0.D0
      IF(.NOT.TSTRESS) THEN
        DO IAT=1,NAT
!         == CORRECT MASS FOR EFFECTIVE MASS OF THE ELECTRONS ============
          RMASS0=RMASS(IAT)-EFFEMASS(IAT)
!         == CORRECT FRICTION FOR ELECTRON FRICTION  =====================
          ANNER1=(RMASS(IAT)*ANNER-EFFEMASS(IAT)*ANNEE)/RMASS0
!         ==  PROPAGATE ==================================================
          SVAR1=2.D0/(1.D0+ANNER1)
          SVAR2=1.D0-SVAR1
          SVAR3=DT**2/RMASS0/(1.D0+ANNER1)
          DO I=1,3
            RP(I,IAT)=SVAR1*R0(I,IAT)+SVAR2*RM(I,IAT)+SVAR3*FORCE(I,IAT)
            V(I)=(RP(I,IAT)-RM(I,IAT))/(2.D0*DT)
          ENDDO
          DO I=1,3
            DO J=1,3
              CELLKIN(I,J)=CELLKIN(I,J)+RMASS0*V(I)*V(J)
            ENDDO
          ENDDO
        ENDDO
      ELSE 
        DO IAT=1,NAT
!         == CORRECT MASS FOR EFFECTIVE MASS OF THE ELECTRONS ============
          RMASS0=RMASS(IAT)-EFFEMASS(IAT)
!         == CORRECT FRICTION FOR ELECTRON FRICTION  =====================
          ANNER1=(RMASS(IAT)*ANNER-EFFEMASS(IAT)*ANNEE)/RMASS0
!         ==  PROPAGATE ==================================================
          DO I=1,3 
            DO J=1,3
              MATP(I,J)=CELLFRIC(I,J)
              MATM(I,J)=-CELLFRIC(I,J)
            ENDDO
            MATP(I,I)=1.D0+ANNER1+MATP(I,I)
            MATM(I,I)=1.D0-ANNER1+MATM(I,I)
          ENDDO
          CALL lib$INVERTr8(3,MATP,MATPINV)
!CALL INVERT(3,MATP,MATPINV)
          RP(:,IAT)=MATMUL(MATPINV,2.D0*R0(:,IAT)-MATMUL(MATM,RM(:,IAT)) &
     &                                  +FORCE(:,IAT)*DT**2/RMASS0)
          DO I=1,3
            V(I)=(RP(I,IAT)-RM(I,IAT))/(2.D0*DT)
          ENDDO
          DO I=1,3
            DO J=1,3
              CELLKIN(I,J)=CELLKIN(I,J)+RMASS0*V(I)*V(J)
            ENDDO
          ENDDO
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
      USE ATOMS_MODULE,ONLY : R0,RM,NAT,name
      USE RESTART_INTERFACE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN) :: NFIL
      INTEGER(4)            ,INTENT(IN) :: NFILO
      LOGICAL(4)            ,INTENT(OUT):: TCHK
      TYPE (SEPARATOR_TYPE),PARAMETER  :: MYSEPARATOR &
               =SEPARATOR_TYPE(3,'ATOMS','NONE','AUG2003','NONE')
!     TYPE (SEPARATOR_TYPE),PARAMETER  :: MYSEPARATOR &
!              =SEPARATOR_TYPE(3,'ATOMS','NONE','AUG1996','NONE')
      INTEGER(4)                        :: NTASKS,ITASK
!     ******************************************************************
      CALL WRITESEPARATOR(MYSEPARATOR,NFIL,NFILO,TCHK)
!
!     ==================================================================
!     ==  WRITE DATA                                                  ==
!     ==================================================================
      CALL MPE$QUERY(NTASKS,ITASK)
      IF(ITASK.EQ.1) THEN
        WRITE(NFIL)NAT
        WRITE(NFIL)name(:)
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
      TYPE (SEPARATOR_TYPE),PARAMETER   :: MYSEPARATOR_old1 &
                 =SEPARATOR_TYPE(3,'ATOMS','NONE','AUG1996','NONE')
      TYPE (SEPARATOR_TYPE),PARAMETER   :: MYSEPARATOR &
                 =SEPARATOR_TYPE(3,'ATOMS','NONE','AUG2003','NONE')
      TYPE (SEPARATOR_TYPE)             :: SEPARATOR
      INTEGER(4)                    :: NAT1 
      REAL(8)          ,ALLOCATABLE :: TMP(:,:)   !(3,NAT1)
      character(32)    ,ALLOCATABLE :: name1(:)   !(NAT1)
      logical(4)                    :: tchk1        
      INTEGER(4)                    :: NTASKS,ITASK
      INTEGER(4)                    :: iat,iat1
!     ******************************************************************
      TCHK=.NOT.START
      SEPARATOR=MYSEPARATOR
      CALL READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      IF(.NOT.TCHK) RETURN
!
      CALL MPE$QUERY(NTASKS,ITASK)
!
!     ==================================================================
!     ==  READ DATA                                                   ==
!     ==================================================================
      IF(ITASK.EQ.1) THEN
        IF(SEPARATOR%VERSION.EQ.MYSEPARATOR%VERSION) THEN
          READ(NFIL)NAT1
          ALLOCATE(TMP(3,NAT1))
          ALLOCATE(NAME1(NAT1))
          READ(NFIL)NAME1(:)
          READ(NFIL)TMP(:,:)
!         == read r(0)
          DO IAT=1,NAT
            tchk1=.false.
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
!         == read r(-)
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
      CALL MPE$BROADCAST(1,R0)
      CALL MPE$BROADCAST(1,RM)
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
      ALLOCATE(Z(NAT))       ;Z=0.D0
      ALLOCATE(ZV(NAT))      ;ZV=0.D0
      ALLOCATE(CHARGE(NAT))  ;CHARGE=0.D0
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
        CALL ERROR$STOP('ATOMLIST$GET')
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
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          VAL_=Z
        ELSE 
           CALL ERROR$MSG('GETR8A DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$GETR8A')
        END IF
!     ==================================================================
!     ==  NUMBER OF VALENCE ELECTRONS                                  ==
!     ==================================================================
      ELSE IF(ID_.EQ.'ZVALENCE') THEN
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
      CALL ERROR$I4VAL('LENG_',LENG_)
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
        VAL_=Z(IAT_)
      ELSE IF(ID_.EQ.'ZVALENCE') THEN
        VAL_=ZV(IAT_)
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
           CALL ERROR$MSG('GETchA DOES NOT HANDLE SCALARS')
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
          IF(LENG_.NE.NAT) GOTO 9000
          VAL_=ISPECIES
        ELSE
           CALL ERROR$MSG('GETI4A DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$GETR8A')
        END IF
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$GETI4A')
      END IF
      RETURN
 9000 CONTINUE
      CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
      CALL ERROR$I4VAL('LENG_',LENG_)
      CALL ERROR$I4VAL('NAT',NAT)
      CALL ERROR$STOP('ATOMLIST$GETI4A')
      STOP
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
        IF(IAT_.EQ.0) THEN
          IF(LENG_.NE.NAT) GOTO 9000
          Z(:)=VAL_(:)
        ELSE 
           CALL ERROR$MSG('SETR8A DOES NOT HANDLE SCALARS')
           CALL ERROR$STOP('ATOMLIST$SETR8A')
        END IF
!     ==================================================================
!     ==  NUMBER OF VALENCE ELECTRONS                                  ==
!     ==================================================================
      ELSE IF(ID_.EQ.'ZVALENCE') THEN
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
!     ******************************************************************
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
      ELSE IF(ID_.EQ.'Z') THEN
        Z(IAT_)=VAL_
      ELSE IF(ID_.EQ.'ZVALENCE') THEN
        ZV(IAT_)=VAL_
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
      INTEGER(4)  ,INTENT(IN):: VAL_(LENG_)
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
      ELSE IF(ID_.EQ.'Z') THEN
        VAL_=TRANSFER(Z(IAT_),VAL_,8)
      ELSE IF(ID_.EQ.'ZVALENCE') THEN
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
      SUBROUTINE ATOMLIST$GETALL(ID_,NBYTE_,NAT_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: NBYTE_
      INTEGER(4)  ,INTENT(IN) :: NAT_
      CHARACTER(1),INTENT(OUT):: VAL_(NBYTE_*NAT_)
!     ******************************************************************
      IF(NAT_.NE.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('NAT_',NAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$GETALL')
      END IF
      IF(ID_.EQ.'NAME') THEN
        VAL_(:)=TRANSFER(NAME,VAL_,32*NAT)
      ELSE IF(ID_.EQ.'ISPECIES') THEN
        VAL_(:)=TRANSFER(ISPECIES,VAL_,4*NAT)
      ELSE IF(ID_.EQ.'R(0)') THEN
        VAL_(:)=TRANSFER(R0,VAL_,8*3*NAT)
      ELSE IF(ID_.EQ.'R(-)') THEN
        VAL_(:)=TRANSFER(RM,VAL_,8*3*NAT)
      ELSE IF(ID_.EQ.'R(+)') THEN
        VAL_(:)=TRANSFER(RP,VAL_,8*3*NAT)
      ELSE IF(ID_.EQ.'FORCE') THEN
        VAL_(:)=TRANSFER(FORCE,VAL_,8*3*NAT)
      ELSE IF(ID_.EQ.'MASS') THEN
        VAL_(:)=TRANSFER(RMASS,VAL_,8*NAT)
      ELSE IF(ID_.EQ.'PS<G2>') THEN
        VAL_(:)=TRANSFER(PSG2,VAL_,8*NAT)
      ELSE IF(ID_.EQ.'PS<G4>') THEN
        VAL_(:)=TRANSFER(PSG4,VAL_,8*NAT)
      ELSE IF(ID_.EQ.'Z') THEN
        VAL_(:)=TRANSFER(Z,VAL_,8*NAT)
      ELSE IF(ID_.EQ.'ZVALENCE') THEN
        VAL_(:)=TRANSFER(ZV,VAL_,8*NAT)
      ELSE IF(ID_.EQ.'POINTCHARGE') THEN
!       == KEYWORD SHALL BE REPLACED BY 'Q'
        VAL_(:)=TRANSFER(CHARGE,VAL_,8*NAT)
      ELSE IF(ID_.EQ.'Q') THEN
        VAL_(:)=TRANSFER(CHARGE,VAL_,8*NAT)
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$GETALL')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$SET(ID_,NBYTE_,IAT_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: NBYTE_
      INTEGER(4)  ,INTENT(IN) :: IAT_
      CHARACTER(1),INTENT(IN) :: VAL_(NBYTE_)
      REAL(8)                 :: R8VAL
!     ******************************************************************
      IF(IAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IAT_',IAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$SET')
      END IF
      IF(ID_.EQ.'NAME') THEN
        NAME(IAT_)=TRANSFER(VAL_,NAME(1))
      ELSE IF(ID_.EQ.'ISPECIES') THEN
        ISPECIES(IAT_)=TRANSFER(VAL_,ISPECIES(1))
      ELSE IF(ID_.EQ.'R(0)') THEN
        R0(:,IAT_)=TRANSFER(VAL_,R0,3)
      ELSE IF(ID_.EQ.'R(-)') THEN
        RM(:,IAT_)=TRANSFER(VAL_,RM,3)
      ELSE IF(ID_.EQ.'R(+)') THEN
        RP(:,IAT_)=TRANSFER(VAL_,RP,3)
      ELSE IF(ID_.EQ.'FORCE') THEN
        FORCE(:,IAT_)=TRANSFER(VAL_,FORCE,3)
      ELSE IF(ID_.EQ.'MASS') THEN
        RMASS(IAT_)=TRANSFER(VAL_,R8VAL)
      ELSE IF(ID_.EQ.'PS<G2>') THEN
        PSG2(IAT_)=TRANSFER(VAL_,R8VAL)
      ELSE IF(ID_.EQ.'PS<G4>') THEN
        PSG4(IAT_)=TRANSFER(VAL_,R8VAL)
      ELSE IF(ID_.EQ.'Z') THEN
        Z(IAT_)=TRANSFER(VAL_,R8VAL)
      ELSE IF(ID_.EQ.'ZVALENCE') THEN
        ZV(IAT_)=TRANSFER(VAL_,R8VAL)
      ELSE IF(ID_.EQ.'POINTCHARGE') THEN
!       == KEYWORD SHALL BE REPLACED BY 'Q'
        CHARGE(IAT_)=TRANSFER(VAL_,R8VAL)
      ELSE IF(ID_.EQ.'Q') THEN
        CHARGE(IAT_)=TRANSFER(VAL_,R8VAL)
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$SET')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$SETALL(ID_,NBYTE_,NAT_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: NBYTE_
      INTEGER(4)  ,INTENT(IN) :: NAT_
      CHARACTER(1),INTENT(IN) :: VAL_(NBYTE_,NAT_)
      INTEGER(4)              :: IAT
!     ******************************************************************
      IF(NAT_.GT.NAT) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('NAT_',NAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('ATOMLIST$SETALL')
      END IF
      IF(ID_.EQ.'NAME') THEN
        NAME=TRANSFER(VAL_,NAME,NAT)
      ELSE IF(ID_.EQ.'ISPECIES') THEN
        ISPECIES=TRANSFER(VAL_,ISPECIES,NAT)
      ELSE IF(ID_.EQ.'R(0)') THEN
        R0=RESHAPE(TRANSFER(VAL_,R0,3*NAT),(/3,NAT/))
      ELSE IF(ID_.EQ.'R(-)') THEN
        RM=RESHAPE(TRANSFER(VAL_,RM,3*NAT),(/3,NAT/))
      ELSE IF(ID_.EQ.'R(+)') THEN
        RP=RESHAPE(TRANSFER(VAL_,RP,3*NAT),(/3,NAT/))
      ELSE IF(ID_.EQ.'FORCE') THEN
        FORCE=RESHAPE(TRANSFER(VAL_,FORCE,3*NAT),(/3,NAT/))
       ELSE IF(ID_.EQ.'MASS') THEN
         RMASS=TRANSFER(VAL_,RMASS,NAT)
       ELSE IF(ID_.EQ.'PS<G2>') THEN
         PSG2=TRANSFER(VAL_,PSG2,NAT)
       ELSE IF(ID_.EQ.'PS<G4>') THEN
         PSG4=TRANSFER(VAL_,PSG4,NAT)
       ELSE IF(ID_.EQ.'Z') THEN
         Z=TRANSFER(VAL_,Z,NAT)
       ELSE IF(ID_.EQ.'ZVALENCE') THEN
         ZV=TRANSFER(VAL_,ZV,NAT)
      ELSE IF(ID_.EQ.'POINTCHARGE') THEN
!       == KEYWORD SHALL BE REPLACED BY 'Q'
         CHARGE=TRANSFER(VAL_,CHARGE,NAT)
      ELSE IF(ID_.EQ.'Q') THEN
         CHARGE=TRANSFER(VAL_,CHARGE,NAT)
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('ATOMLIST$SETALL')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$SETLATTICE(RBAS_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: RBAS_(3,3)
      REAL(8)            :: RBAS(3,3)
!     ******************************************************************
      CALL ERROR$MSG('ROUTINE MARKED FOR DELETION')
      CALL ERROR$STOP('ATOMLIST$SETLATTICE')
      RBAS(:,:)=RBAS_(:,:)
      CALL CELL$GETR8A('T0',9,RBAS)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMLIST$LATTICE(RBAS_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE ATOMS_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: RBAS_(3,3)
      REAL(8)             :: RBAS(3,3)
!     ******************************************************************
      CALL ERROR$MSG('ROUTINE MARKED FOR DELETION')
      CALL ERROR$STOP('ATOMLIST$LATTICE')
      RBAS_(:,:)=RBAS(:,:)
      CALL CELL$GETR8A('T0',9,RBAS_)
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
      REAL(8)       :: XZ            
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
              XZ           =Z(IAT)       
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
              Z(IAT)           =Z(IAT1)          
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
              Z(IAT1)           =XZ           
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
      CHARACTER(32)           :: CHT
      INTEGER(4)              :: ICH
      LOGICAL(4)              :: TCHK
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
        TCHK=.TRUE.
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
      CALL MPE$QUERY(NTASKS,THISTASK) 
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
