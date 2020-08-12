!ATTENTION CHARGES MAY BE REVERSED BECAUSE I COUNT ELECTRONS AND 
! FORCE FIELDS COUNT ELECTRON CHARGES!
! 
!.......................................................................
MODULE QMMM_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: QMMM                                                       **
!**                                                                   **
!**  PURPOSE: COUPLES A QUANTUM MECHANICAL (QM) MOLECULE ON AN        **
!**  MOLECULAR MECHANICAL (MM) ENVIRONMENT.                           **
!**                                                                   **
!**  THREE SYSTEMS ARE TREATED:                                       **
!**  1) THE QM-MOLECULE (NOT IN THIS OBJECT)                          **
!**  2) THE MM-SHADOW OF THE QM-MOLECULE                              **
!**  3) THE MM ENVIRONMENT AND ITS EMBEDDED MOLECULE                  **
!**  THE TOTAL ENERGY OF THE COMBINED SYSTEM IS E1+(E3_E2)            **
!**                                                                   **
!**  METHODS:                                                         **
!**    QMMM$INTERFACE                                                 **
!**    QMMM$SCALEFORCE(NAT,QFORCE)                                    **
!**    QMMM$PROPAGATE                                                 **
!**    QMMM$SWITCH                                                    **
!**                                                                   **
!**  USES:                                                            **
!**    CLASSICAL                                                      **
!**                                                                   **
!**  REMARKS:                                                         **
!**  THE ENERGIES FROM THIS MODULE ARE ACCESSED VIA                   **
!**      QMMM$GETR8('EPOT',VAL)                                       **
!**      QMMM$GETR8('EKIN',VAL)                                       **
!**      QMMM$GETR8('ETHERM',VAL)                                     **
!**  THE ENERGYLIST IS NOT USED IN THIS OBJECT                        **
!**                                                                   **
!***********************************************************************

TYPE LINK_TYPE
! LINK DESCRIBES A BOND ACROSS THE QM-MM BOUNDARY. Q/M/SJOINT ARE THE 
! ATOM INDICES FOR THE COMMON ATOM OF ALL THREE SUBSYSTEMS.
! Q/S ATOM REFERS TO THE DUMMY ATOM IN THE BOND AND M ATOM IS THE 
! ATOM INDEX OF THE MM-SYSTEM CONNECTED TO THE QM-SYSTEM.
  INTEGER(4)  :: QJOINT
  INTEGER(4)  :: MJOINT
  INTEGER(4)  :: SJOINT
  INTEGER(4)  :: QATOM
  INTEGER(4)  :: MATOM
  INTEGER(4)  :: SATOM
  REAL(8)     :: ALPHA
  REAL(8)     :: QAO    !Q-ALL OUTER
  REAL(8)     :: QSO    !Q-SHADOW-OUTER
  REAL(8)     :: FAO(3) !FORCE-ALL-OUTER
  REAL(8)     :: FAI(3) !FORCE-ALL-INNER
  REAL(8)     :: FSO(3) !FORCE-SHADOW-OUTER
  REAL(8)     :: FSI(3) !FORCE-SHADOW-INNER
  INTEGER(4)  :: SHARED      ! POINTS TO ANOTHER LINK BOND WITH SHARED ATOMS
END TYPE LINK_TYPE
TYPE MAP_TYPE
! MAP PROVIDES THE INDICES FOR THE ATOMS OF THE CENTRAL CLUSTER 
! NOT PARTICIPATING IN THE LINK BONDS.
! THE ATOM HAS THE INDEX QATOM IN THE QM-SYSTEM, THE INDEX SATOM IN THE 
! SHADOW AND THE INDEX MATOM IN THE MM SYSTEM. MAP HAS ENTRIES ALSO TO 
! BOTH PARTNERS OF A BOND CROSSING THE QM-MM BOUNDARY. (IT INCLUDES ALSO 
! THE DUMMY HYDROGEN IN THE BOND.
  INTEGER(4)  :: QATOM
  INTEGER(4)  :: MATOM
  INTEGER(4)  :: SATOM
END TYPE MAP_TYPE
LOGICAL(4)                  :: TON  =.FALSE. ! ON/OFF SWITCH
LOGICAL(4)                  :: TINI =.FALSE. ! FLAG FOR INITIALIZATION
LOGICAL(4)                  :: TMOVE=.FALSE. ! MOVE ATOMS OR FREEZE MOTION
LOGICAL(4)                  :: TADIABATIC=.FALSE.  !MINIMIZATION IN EACH TIME STEP
REAL(8)                     :: DELTAT=10.D0  ! TIMESTEP
REAL(8)                     :: ANNE =0.D0    ! FRICTION
LOGICAL(4)                  :: TSTOP=.TRUE.  ! SET VELOCITIES TO ZERO
LOGICAL(4)                  :: TRANDOMIZE=.FALSE.  ! RANDOMIZE VELOCITIES
REAL(8)                     :: AMPRE=0.D0    ! K_B*T FOR RANDOMIZATION
INTEGER(4)                  :: NMULTIPLE=1   ! #(MULTIPLE TIME STEPS)
INTEGER(4)                  :: NATM=0        ! #(ATOMS IN THE QM+MM)
INTEGER(4)                  :: NATS=0        ! #(SHADOW ATOMS)
INTEGER(4)                  :: NATQ=0        ! #(QM ATOMS)
INTEGER(4)                  :: NLINK=0       ! #(LINK BONDS)
TYPE(LINK_TYPE),ALLOCATABLE :: LINK(:)       ! DEFINES LINK-BONDS
INTEGER(4)                  :: NMAP=0        ! #(QM-ATOMS EXCEPT DUMMY ATOMS)
TYPE(MAP_TYPE) ,ALLOCATABLE :: MAP(:)        ! REACTION CENTER + DUMMY ATOMS
REAL(8)                     :: EPOT_QMMM=0.D0   ! EPOT OF THE ENVIRONMENT
REAL(8)                     :: EKIN_QMMM=0.D0   ! EKIN OF THE ENVIRONMENT
REAL(8)                     :: ETHERM_QMMM=0.D0 ! EKIN+EPOT OF THE THERMOSTAT
END MODULE QMMM_MODULE      
!     ..................................................................
      SUBROUTINE QMMM_INITIALIZE
!     ******************************************************************
!     ******************************************************************
      USE QMMM_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)               :: IATQJ,IATMJ,IATQ,IATM,ILINK,IATS
      REAL(8)                  :: RS,RJ,RM
      REAL(8)     ,ALLOCATABLE :: MQ(:)
      REAL(8)     ,ALLOCATABLE :: SQ(:)
      CHARACTER(5),ALLOCATABLE :: MTYPE(:)
      CHARACTER(5),ALLOCATABLE :: STYPE(:)
      CHARACTER(32),ALLOCATABLE:: MNAME(:)
      CHARACTER(32),ALLOCATABLE:: SNAME(:)
      CHARACTER(2), ALLOCATABLE:: MELEMENT(:), SELEMENT(:)
      INTEGER(4)               :: I
      LOGICAL(4)               :: TCHK
!     ******************************************************************
      IF(TINI) RETURN
                        CALL TRACE$PUSH('QMMM_INITIALIZE')
      TINI=.TRUE.
!
!     ==================================================================
!     == COLLECT #(ATOMS)                                             ==
!     ==================================================================
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$GETI4('NAT',NATM)
      CALL CLASSICAL$SELECT('SHADOW')
      CALL CLASSICAL$GETI4('NAT',NATS)
      CALL ATOMLIST$NATOM(NATQ)
      IF(NATS.NE.NATQ) THEN
        CALL ERROR$MSG('#(QM ATOMS) NOT EQUAL #(SHADOW ATOMS)')
        CALL ERROR$STOP('QMMM_INITIALIZE')
      END IF
!
!     ==================================================================
!     ==  SET ALPHA AND M-CHARGE                                      ==
!     ==================================================================
      ALLOCATE(MQ(NATM))
      ALLOCATE(SQ(NATS))
      ALLOCATE(MTYPE(NATM))
      ALLOCATE(STYPE(NATS))
      ALLOCATE(MNAME(NATM))
      ALLOCATE(MELEMENT(NATM))
      ALLOCATE(SNAME(NATS))
      ALLOCATE(SELEMENT(NATS))
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$GETR8A('QEL',NATM,MQ)
      CALL CLASSICAL$GETCHA('TYPE',NATM,MTYPE)
      CALL CLASSICAL$GETCHA('ATOMNAME',NATM,MNAME)
      CALL CLASSICAL$GETCHA('ELEMENT',NATM,MELEMENT)
      CALL CLASSICAL$SELECT('SHADOW')
      CALL CLASSICAL$GETCHA('TYPE',NATS,STYPE)
      CALL CLASSICAL$GETCHA('ATOMNAME',NATS,SNAME)
      CALL CLASSICAL$GETCHA('ELEMENT',NATS,SELEMENT)

      DO I=1,SIZE(MNAME)
         MNAME(I)(1:2)= ADJUSTL(MNAME(I)(1:2))
         IF(MNAME(I)(2:2).EQ.' ') MNAME(I)(2:2)='_'
!NEW VERSION:
         IF(MELEMENT(I)(2:2).EQ.' ') MELEMENT(I)(2:2)='_'
!----------
      ENDDO
      DO I=1,SIZE(SNAME)
         SNAME(I)= ADJUSTL(SNAME(I)(1:2))
         IF(SNAME(I)(2:2).EQ.' ') SNAME(I)(2:2)='_'
!NEW VERSION:
         IF(SELEMENT(I)(2:2).EQ.' ') SELEMENT(I)(2:2)='_'
!----------
      ENDDO
      CALL CLASSICAL$GETR8A('QEL',NATS,SQ)
!
      DO ILINK=1,NLINK
        IATQJ=LINK(ILINK)%QJOINT
        IATMJ=LINK(ILINK)%MJOINT
        IATQ=LINK(ILINK)%QATOM
        IATS=LINK(ILINK)%SATOM
        IATM=LINK(ILINK)%MATOM
!NEW VERSION: CHANGED S/MNAME INTO S/MELEMENT
        CALL PERIODICTABLE$GET(SELEMENT(IATS)(1:2),'R(COV)',RS)
        CALL PERIODICTABLE$GET(MELEMENT(IATM)(1:2),'R(COV)',RM)
        CALL PERIODICTABLE$GET(MELEMENT(IATMJ)(1:2),'R(COV)',RJ)

        LINK(ILINK)%ALPHA=(RS+RJ)/(RM+RJ)
        LINK(ILINK)%QAO=MQ(IATM)
        LINK(ILINK)%QSO=SQ(IATS)
        LINK(ILINK)%FAO(:)=0.D0
        LINK(ILINK)%FAI(:)=0.D0
        LINK(ILINK)%FSO(:)=0.D0
        LINK(ILINK)%FSI(:)=0.D0
        LINK(ILINK)%SHARED=0
      ENDDO
      DEALLOCATE(STYPE)
      DEALLOCATE(MTYPE)
      DEALLOCATE(MNAME)
      DEALLOCATE(SNAME)
      DEALLOCATE(MELEMENT)
      DEALLOCATE(SELEMENT)
      DEALLOCATE(MQ)
      DEALLOCATE(SQ)
!
!     ===================================================================
!     ==  CHECK IF LINK BONDS OVERLAP                                  ==
!     ===================================================================
      DO ILINK=1,NLINK
        IATQ=LINK(ILINK)%QATOM
        IATQJ=LINK(ILINK)%QJOINT
        LINK(ILINK)%SHARED=0
        DO I=ILINK+1,NLINK
          IF(LINK(I)%QATOM.EQ.IATQ.OR.LINK(I)%QJOINT.EQ.IATQJ) THEN
            LINK(ILINK)%SHARED=I
            EXIT
          END IF
        ENDDO
      ENDDO
!
!     ===================================================================
!     ==  SET THERMOSTAT  FOR THE ENVIRONMENT                          ==
!     ===================================================================
      CALL THERMOSTAT$SELECT('QM-MM')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK) THEN
        CALL THERMOSTAT$SCALEGFREE(REAL(3*(NATM-NATS),KIND=8))
      END IF
                         CALL TRACE$POP
      RETURN
      END SUBROUTINE QMMM_INITIALIZE
! 
!     ..................................................................
      SUBROUTINE QMMM$SETI4A(ID,LEN,VALUE)
!     *******************************************************************
!     *******************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: VALUE(LEN)
      INTEGER(4)              :: I,IMAP,ILINK
!     *******************************************************************
      IF(ID.EQ.'MAP') THEN
        IF(ALLOCATED(MAP)) DEALLOCATE(MAP)
        NMAP=LEN/3
        IF(LEN.NE.3*NMAP) THEN
          CALL ERROR$MSG('LENGTH INCONSISTENT')
          CALL ERROR$STOP('QMMM$SETI4A')
        END IF
        ALLOCATE(MAP(NMAP))
        I=0
        DO IMAP=1,NMAP
          I=I+1; MAP(IMAP)%QATOM=VALUE(I)
          I=I+1; MAP(IMAP)%MATOM=VALUE(I)
          I=I+1; MAP(IMAP)%SATOM=VALUE(I)
        ENDDO
      ELSE IF(ID.EQ.'LINK') THEN
        IF(ALLOCATED(LINK)) DEALLOCATE(LINK)
        NLINK=LEN/6
        IF(LEN.NE.6*NLINK) THEN
          CALL ERROR$MSG('LENGTH INCONSISTENT')
          CALL ERROR$STOP('QMMM$SETI4A')
        END IF
        ALLOCATE(LINK(NLINK))
        I=0
        DO ILINK=1,NLINK
          I=I+1; LINK(ILINK)%QJOINT=VALUE(I)
          I=I+1; LINK(ILINK)%MJOINT=VALUE(I)
          I=I+1; LINK(ILINK)%SJOINT=VALUE(I)
          I=I+1; LINK(ILINK)%QATOM=VALUE(I)
          I=I+1; LINK(ILINK)%MATOM=VALUE(I)
          I=I+1; LINK(ILINK)%SATOM=VALUE(I)
          LINK(ILINK)%ALPHA=0.D0
          LINK(ILINK)%QAO=0.D0
        ENDDO
      ELSE 
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('QMMM$SETI4A')
      END IF
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE QMMM$GETI4(ID,VALUE)
!     *******************************************************************
!     *******************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT) :: VALUE
!     *******************************************************************
      IF(ID.EQ.'NAT:ENV') THEN
        VALUE=NATM-NATS
      ELSE IF(ID.EQ.'NAT:ALL') THEN
        VALUE=NATM
      ELSE IF (ID.EQ.'NAT:SHADOW') THEN 
        VALUE=NATS
      ELSE 
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('QMMM$GETI4')
      END IF
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE QMMM$SETL4(ID,VALUE)
!     *******************************************************************
!     *******************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VALUE
!     *******************************************************************
      IF(ID.EQ.'ON') THEN
        TON=VALUE
      ELSE IF(ID.EQ.'STOP') THEN
        TSTOP=VALUE
      ELSE IF(ID.EQ.'MOVE') THEN
        TMOVE=VALUE
      ELSE IF(ID.EQ.'ADIABATIC') THEN
        TADIABATIC=VALUE
      ELSE 
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('QMMM$SETL4')
      END IF
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE QMMM$GETL4(ID,VALUE)
!     *******************************************************************
!     *******************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VALUE
!     *******************************************************************
      IF(ID.EQ.'ON') THEN
        VALUE=TON
      ELSE 
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('QMMM$GETL4')
      END IF
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE QMMM$SETI4(ID,VALUE)
!     *******************************************************************
!     *******************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VALUE
!     *******************************************************************
      IF(ID.EQ.'MULTIPLE') THEN
        NMULTIPLE=VALUE
      ELSE 
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('QMMM$SETL4')
      END IF
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE QMMM$SETR8(ID,VALUE)
!     *******************************************************************
!     *******************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VALUE
!     *******************************************************************
      IF(ID.EQ.'TIMESTEP') THEN
        DELTAT=VALUE
      ELSE IF(ID.EQ.'FRICTION') THEN
        ANNE=VALUE
      ELSE IF(ID.EQ.'RANDOM') THEN
        AMPRE=VALUE
        TRANDOMIZE=.TRUE.
      ELSE 
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('QMMM$SETR8')
      END IF
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE QMMM$GETR8(ID,VALUE)
!     *******************************************************************
!     *******************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VALUE
!     *******************************************************************
      IF(ID.EQ.'EKIN') THEN
        VALUE=EKIN_QMMM
      ELSE IF(ID.EQ.'EPOT') THEN
        VALUE=EPOT_QMMM
      ELSE IF(ID.EQ.'ETHERM') THEN
        VALUE=ETHERM_QMMM
      ELSE 
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('QMMM$GETR8')
      END IF
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE QMMM$REPORT(NFIL)
!     *****************************************************************
!     *****************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      REAL(8)                :: KELVIN
      INTEGER(4)             :: NTASKS,THISTASK
!     *****************************************************************
      IF(.NOT.TON) RETURN
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
                              CALL TRACE$PUSH('QMMM$REPORT')
      IF(THISTASK.EQ.1) THEN
        CALL CONSTANTS('KB',KELVIN)
        CALL REPORT$TITLE(NFIL,'QM-MM COUPLING')
        IF(TSTOP)CALL REPORT$STRING(NFIL,'ZERO INITIAL VELOCITIES OF ENVIRONMENT')
        IF(TADIABATIC)CALL REPORT$STRING(NFIL,'MM ATOMS RELAXED IN EACH TIME STEP')
        CALL REPORT$I4VAL(NFIL,'ENVIRONMENT OVERSAMPLED BY FACTOR',NMULTIPLE,' ')
        IF(.NOT.TMOVE)CALL REPORT$STRING(NFIL,'ENVIRONMENT FROZEN')
        CALL REPORT$R8VAL(NFIL,'FRICTION',ANNE,' ')
        IF(TRANDOMIZE) THEN
          CALL REPORT$R8VAL(NFIL,'INITIAL VELOCITIES RANDOMIZED WITH',AMPRE/KELVIN,'K')
        END IF
        CALL REPORT$I4VAL(NFIL,'NUMBER OF ATOMS INCLUDING ENVIRONMENT',NATM,' ')
        CALL REPORT$I4VAL(NFIL,'NUMBER OF ATOMS IN SHADOW',NATS,' ')
!
!     ==================================================================
!     == REPORT POSITIONS, BONDS FOR THE MM-ENVIRONMENT               ==
!     ==================================================================
        CALL REPORT$STRING(NFIL,'MM-ENVIRONMENT')
      END IF
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$REPORT(NFIL)
                             CALL TRACE$POP
      RETURN
      END

! 
!     ..................................................................
      SUBROUTINE QMMM$INTERFACE(NAT,POS,CHARGE,FORCE,POT,DEPOT)
!     *****************************************************************
!     *****************************************************************
      USE QMMM_MODULE !, ONLY : MAP_TYPE,NMAP,MAP,LINK_TYPE,NLINK,LINK
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NAT          ! #(ATOMS IN QM CLUSTER)
      REAL(8)   ,INTENT(IN)  :: POS(3,NAT)   ! POSITION OF QM CLUSTER
      REAL(8)   ,INTENT(IN)  :: CHARGE(NAT)  ! CHARGES OF QM CLUSTER
      REAL(8)   ,INTENT(OUT) :: FORCE(3,NAT) ! ENVIRONMENT FORCES ON QM CLUSTER
      REAL(8)   ,INTENT(OUT) :: POT(NAT)     ! ENVIRONMENT POTENTIALS ON QM CLUSTER
      REAL(8)   ,INTENT(OUT) :: DEPOT        ! ENVIRONMENT ENERGY
      REAL(8)   ,ALLOCATABLE :: MPOS(:,:)    ! POSITIONS OF COMPLETE SYSTEM
      REAL(8)   ,ALLOCATABLE :: MPOSM(:,:)
      REAL(8)   ,ALLOCATABLE :: SPOS(:,:)
      REAL(8)   ,ALLOCATABLE :: SPOSM(:,:)
      REAL(8)   ,ALLOCATABLE :: MCHARGE(:)
      REAL(8)   ,ALLOCATABLE :: SCHARGE(:)
      REAL(8)   ,ALLOCATABLE :: MFORCE(:,:)
!     REAL(8)   ,ALLOCATABLE :: SFORCE(:,:)
      REAL(8)   ,ALLOCATABLE :: MPOT(:)
!     REAL(8)   ,ALLOCATABLE :: SPOT(:)
      REAL(8)   ,ALLOCATABLE :: FA(:,:)
      REAL(8)   ,ALLOCATABLE :: FS(:,:)
      REAL(8)   ,ALLOCATABLE :: VA(:)
      REAL(8)   ,ALLOCATABLE :: VS(:)
      REAL(8)                :: EPOTM,EPOTS
      REAL(8)                :: EKINMM
      INTEGER(4)             :: IMAP,ILINK
      INTEGER(4)             :: IATQ,IATM,IATS,IATQJ,IATMJ,IATSJ
      REAL(8)                :: ALPHA
      LOGICAL(4),PARAMETER   :: TPR=.TRUE.
      REAL(8)                :: SVAR
      LOGICAL(4)             :: TCHK
      INTEGER(4)             :: NFILO
      INTEGER(4)             :: IMULTIPLE
      LOGICAL(4)             :: TTHERMOSTAT
      REAL(8)                :: ENOSE
      INTEGER(4)             :: IND
!     *****************************************************************
      DEPOT=0.D0
      POT(:)=0.D0
      FORCE(:,:)=0.D0
      IF(.NOT.TON) RETURN
                        CALL TRACE$PUSH('QMMM$INTERFACE')
      CALL QMMM_INITIALIZE
                        CALL TIMING$CLOCKON('QM-MM')
      ALLOCATE(MPOS(3,NATM))
      ALLOCATE(MPOSM(3,NATM))
      ALLOCATE(MCHARGE(NATM))
      ALLOCATE(MFORCE(3,NATM))
      ALLOCATE(MPOT(NATM))
      ALLOCATE(SPOS(3,NATS))
      ALLOCATE(SPOSM(3,NATS))
      ALLOCATE(SCHARGE(NATS))
!     ALLOCATE(SFORCE(3,NATS))
!     ALLOCATE(SPOT(NATS))
      ALLOCATE(FA(3,NATQ))
      ALLOCATE(FS(3,NATS))
      ALLOCATE(VA(NATQ))
      ALLOCATE(VS(NATS))
!
!     ==================================================================
!     ==  APPLY CONSTRAINTS                                           ==
!     ==  AND SET VELOCITIES FOR REACTION CENTER AND LINK ATOMS       ==
!     ==  IN THE MM-PART TO ZERO                                      ==
!     ==================================================================
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$GETR8A('R(0)',3*NATM,MPOS)
      CALL CLASSICAL$GETR8A('R(-)',3*NATM,MPOSM)
      CALL CLASSICAL$GETR8A('QEL',NATM,MCHARGE)
      CALL CLASSICAL$SELECT('SHADOW')
      CALL CLASSICAL$GETR8A('R(0)',3*NATS,SPOS)
      CALL CLASSICAL$GETR8A('R(-)',3*NATS,SPOSM)
      CALL CLASSICAL$GETR8A('QEL',NATS,SCHARGE)
!
!     =========================================================================
!     ==  FOR THE ATOMS SHARED BETWEEN Q,M AND S, SET POSITIONS AND CHARGES  ==
!     ==  EQUAL TO Q-POSITIONS                                               ==
!     ==  M-VELOCITIES ARE SET TO ZERO                                       ==
!     =========================================================================
      DO IMAP=1,NMAP
        IATQ=MAP(IMAP)%QATOM
        IATM=MAP(IMAP)%MATOM
        IATS=MAP(IMAP)%SATOM
        MPOS(:,IATM)=POS(:,IATQ)
        MPOSM(:,IATM)=MPOS(:,IATM)   !SET VELOCITY TO ZERO FOR CENTRAL CLUSTER
        SPOS(:,IATS)=POS(:,IATQ)
        MCHARGE(IATM)=CHARGE(IATQ)
        SCHARGE(IATS)=CHARGE(IATQ)
!SET TO ZERO TO NEGLECT CHARGES
! MCHARGE(IATM)=0.D0
! SCHARGE(IATS)=0.D0
      ENDDO      
!
!     =========================================================================
!     ==  FOR THE LINK-BONDS, DETERMINE THE POSITION OF THE DUMMY ATOM       ==
!     ==  AND THE LINKED ENVIRONMENT ATOM IN M FROM THE Q-POSITIONS          ==
!     =========================================================================
      DO ILINK=1,NLINK
        IATQJ=LINK(ILINK)%QJOINT
        IATMJ=LINK(ILINK)%MJOINT
        IATSJ=LINK(ILINK)%SJOINT
        IATQ =LINK(ILINK)%QATOM
        IATM =LINK(ILINK)%MATOM
        IATS =LINK(ILINK)%SATOM
        ALPHA=LINK(ILINK)%ALPHA
        MPOS(:,IATM)  =POS(:,IATQJ)+(POS(:,IATQ)-POS(:,IATQJ))/ALPHA
        MPOS(:,IATMJ) =POS(:,IATQJ)   
        MPOSM(:,IATM) =MPOS(:,IATM)   !SET VELOCITIES TO ZERO
        MPOSM(:,IATMJ)=MPOS(:,IATMJ)  !SET VELOCITIES TO ZERO
        SPOS(:,IATS)=POS(:,IATQ)
        SPOS(:,IATSJ)=POS(:,IATQJ)
!       == NOTE THAT THE CHARGE OF THE M ATOM MUST BE RESET!! ========
        SVAR=CHARGE(IATQ)-LINK(ILINK)%QSO
        MCHARGE(IATM) =LINK(ILINK)%QAO+ALPHA*SVAR
        MCHARGE(IATMJ)=CHARGE(IATQJ)+(1.D0-ALPHA)*SVAR
!DOUBLE COUNTING FOR SHARED ATOMS!
        SCHARGE(IATS) =CHARGE(IATQ)
        SCHARGE(IATSJ)=CHARGE(IATQJ)

!SET TO ZERO TO NEGLECT CHARGES
!MCHARGE(IATM)=0.D0
!MCHARGE(IATMJ)=0.D0
!SCHARGE(IATS)=0.D0
!SCHARGE(IATSJ)=0.D0
      ENDDO
PRINT*,'MQ ',MCHARGE
PRINT*,'SQ ',SCHARGE
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$SETR8A('R(0)',3*NATM,MPOS)
      CALL CLASSICAL$SETR8A('R(-)',3*NATM,MPOSM)
      CALL CLASSICAL$SETR8A('QEL',NATM,MCHARGE)
      CALL CLASSICAL$SELECT('SHADOW')
      CALL CLASSICAL$SETR8A('R(0)',3*NATS,SPOS)
      CALL CLASSICAL$SETR8A('QEL',NATS,SCHARGE)
!
!     ==================================================================
!     ==  TOTAL ENERGY AND FORCES OF THE SHADOW                       ==
!     ==================================================================
      CALL CLASSICAL$SELECT('SHADOW')
      CALL CLASSICAL$NEIGHBORS
      CALL CLASSICAL$ETOT(EPOTS)  
      CALL CLASSICAL$GETR8A('FORCE',3*NATS,FS)
      CALL CLASSICAL$GETR8A('VEL',NATS,VS)
!
!     ==================================================================
!     ==  SET VELOCITIES TO ZERO IF REQUESTED                         ==
!     ==================================================================
      IF(TSTOP) THEN
        CALL CLASSICAL$SELECT('QMMM')
        CALL CLASSICAL$GETR8A('R(0)',3*NATM,MPOS)
        CALL CLASSICAL$SETR8A('R(-)',3*NATM,MPOS)
        CALL CLASSICAL$SELECT('SHADOW')
        CALL CLASSICAL$GETR8A('R(0)',3*NATS,SPOS)
        CALL CLASSICAL$SETR8A('R(-)',3*NATS,SPOS)
        TSTOP=.FALSE.
      END IF
!
!     ==================================================================
!     ==  RANDOMIZE INITIAL VELOCITIES                                ==
!     ==================================================================
      IF(TRANDOMIZE) THEN
!       == THE VELOCITIES OF THE REACTION CENTER REMAIN UNCHANGED
        CALL QMMM_RANDOMIZEVELOCITY
        TRANDOMIZE=.FALSE.
      END IF
!
!     ==================================================================
!     ==  SELECT THERMOSTAT                                           ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('QM-MM')
      CALL THERMOSTAT$GETL4('ON',TTHERMOSTAT)
!
!     ==================================================================
!     ==  MULTIPLE TIMESTEP LOOP                                      ==
!     ==================================================================
      EPOT_QMMM=0.D0
      EKIN_QMMM=0.D0
      ETHERM_QMMM=0.D0
      FA(:,:)=0.D0
      VA(:)=0.D0
      CALL CLASSICAL$SELECT('QMMM')
      DO IMULTIPLE=1,NMULTIPLE
!
!       ==================================================================
!       ==  CALCULATE FORCES                                            ==
!       ==================================================================
                        CALL TIMING$CLOCKON('QM-MM:NEIGHBORS')
        CALL CLASSICAL$NEIGHBORS
                        CALL TIMING$CLOCKOFF('QM-MM:NEIGHBORS')
        IF(TADIABATIC.AND.TMOVE) THEN
          CALL QMMM_MINIMIZE(TCHK)
        END IF
                        CALL TIMING$CLOCKON('QM-MM:ETOT')
        CALL CLASSICAL$ETOT(EPOTM)  
        EPOT_QMMM=EPOT_QMMM+EPOTM-EPOTS
                        CALL TIMING$CLOCKOFF('QM-MM:ETOT')
!
!       ==================================================================
!       ==  ACCUMULATE FORCES AND POTENTIALS                             ==
!       ==================================================================
        CALL CLASSICAL$GETR8A('FORCE',3*NATM,MFORCE)
        CALL CLASSICAL$GETR8A('VEL',NATM,MPOT)
        DO IMAP=1,NMAP
          IATQ=MAP(IMAP)%QATOM
          IATM=MAP(IMAP)%MATOM
          IATS=MAP(IMAP)%SATOM
!         THE INDEXING OF THE QM SYSTEM IS USED BECAUSE THE ARRAYS ARE 
!         ALLOCATED ONLY WITH NAT ELEMENTS
          FA(:,IATQ)=FA(:,IATQ)+MFORCE(:,IATM)
          VA(IATQ)=VA(IATQ)+MPOT(IATM)
        ENDDO        
        DO ILINK=1,NLINK
          IATQJ=LINK(ILINK)%QJOINT
          IATMJ=LINK(ILINK)%MJOINT
          IATSJ=LINK(ILINK)%SJOINT
          IATQ =LINK(ILINK)%QATOM
          IATM =LINK(ILINK)%MATOM
          IATS =LINK(ILINK)%SATOM
          ALPHA=LINK(ILINK)%ALPHA
          FA(:,IATQ) =FA(:,IATQ) +MFORCE(:,IATM)
          FA(:,IATQJ)=FA(:,IATQJ)+MFORCE(:,IATMJ)
          VA(IATQ)   =VA(IATQ)   +ALPHA*MPOT(IATM)+(1.D0-ALPHA)*MPOT(IATMJ)
          VA(IATQJ)  =VA(IATQJ)  +MPOT(IATMJ)
        ENDDO
!
!       ==================================================================
!       ==  DO NOT PROPAGATE ADIABATIC SOLUTION                         ==
!       ==================================================================
        IF(TADIABATIC.OR.(.NOT.TMOVE)) THEN
          CALL CLASSICAL$GETR8A('R(0)',3*NATM,MPOS)
          CALL CLASSICAL$SETR8A('R(+)',3*NATM,MPOS)
          CALL CLASSICAL$SETR8A('R(-)',3*NATM,MPOS)
          EXIT
        END IF
!
!       ==================================================================
!       ==  PROPAGATE SMALL TIME STEPS                                  ==
!       ==================================================================
        IF(TTHERMOSTAT) THEN
          CALL THERMOSTAT$GETR8('COOLING',ANNE)
        END IF
        CALL CLASSICAL$PROPAGATE(DELTAT,ANNE)
!
!       ==================================================================
!       ==  ENFORCE CONSTRAINTS                                         ==
!       ==================================================================
        CALL CLASSICAL$GETR8A('R(+)',3*NATM,MPOS)
        DO IMAP=1,NMAP
          IATQ=MAP(IMAP)%QATOM
          IATM=MAP(IMAP)%MATOM
          MPOS(:,IATM)=POS(:,IATQ)
        ENDDO      
        DO ILINK=1,NLINK
          IATQJ=LINK(ILINK)%QJOINT
          IATQ =LINK(ILINK)%QATOM
          IATM =LINK(ILINK)%MATOM
          ALPHA=LINK(ILINK)%ALPHA
          MPOS(:,IATM)=POS(:,IATQJ)+(POS(:,IATQ)-POS(:,IATQJ))/ALPHA
          MPOS(:,IATMJ)=POS(:,IATQJ)
        ENDDO
        CALL CLASSICAL$SETR8A('R(+)',3*NATM,MPOS)
!
!       ==================================================================
!       ==  OBTAIN KINETIC ENERGY (MM-ATOMS EXCEPT LINK ATOMS AND       ==
!       ==                        ATOMS OF THE REACTION CENTER)
!       ==================================================================
        CALL CLASSICAL$EKIN(DELTAT,EKINMM)
        EKIN_QMMM=EKIN_QMMM+EKINMM
!
!       ==================================================================
!       ==  PROPAGATE THERMOSTAT                                        ==
!       ==================================================================
        CALL THERMOSTAT$SETR8('EKIN(SYSTEM)',EKINMM)
        CALL THERMOSTAT$PROPAGATE
        CALL THERMOSTAT$GETR8('ENERGY',ENOSE)
        ETHERM_QMMM=ETHERM_QMMM+ENOSE   ! CONTAINS KINETIC AND POTENTIAL ENERGY
!
!       ==================================================================
!       ==  SWITCH EXCEPT THE LAST                                      ==
!       ==================================================================
        IF(IMULTIPLE.LT.NMULTIPLE) THEN
          CALL CLASSICAL$SWITCH
          CALL THERMOSTAT$SWITCH
        END IF
      ENDDO
      IF(.NOT.TADIABATIC) THEN
        SVAR=1.D0/DBLE(NMULTIPLE)
        EKIN_QMMM  =SVAR*EKIN_QMMM
        EPOT_QMMM  =SVAR*EPOT_QMMM
        ETHERM_QMMM=SVAR*ETHERM_QMMM
        FA(:,:)    =SVAR*FA(:,:)
        VA(:)      =SVAR*VA(:)
      ELSE
        EKIN_QMMM  = 0.D0
        ETHERM_QMMM= 0.D0
      END IF
      DEPOT      =EPOT_QMMM
PRINT*,'EKIN_QMMM',EKIN_QMMM
!
!     ==================================================================
!     ==  MAP FORCES AND POTENTIALS INTO ARGUMENTS FORCE AND POT      ==
!     ==  FORCES ON LINK ATOMS ARE DONE SPECIAL                       ==
!     ==================================================================
      POT(:)=0.D0
      DO IMAP=1,NMAP
        IATQ=MAP(IMAP)%QATOM
        IATM=MAP(IMAP)%MATOM
        IATS=MAP(IMAP)%SATOM
!       THE INDEXING OF THE QM SYSTEM IS USED BECAUSE THE ARRAYS ARE 
!       ALLOCATED ONLY WITH NAT ELEMENTS
        FORCE(:,IATQ)=FA(:,IATQ)-FS(:,IATS)
        POT(IATQ)=POT(IATQ)+VA(IATQ)-VS(IATS)
      ENDDO        
      DO ILINK=1,NLINK
        IATQJ=LINK(ILINK)%QJOINT
        IATMJ=LINK(ILINK)%MJOINT
        IATSJ=LINK(ILINK)%SJOINT
        IATQ =LINK(ILINK)%QATOM
        IATM =LINK(ILINK)%MATOM
        IATS =LINK(ILINK)%SATOM
        ALPHA=LINK(ILINK)%ALPHA
        LINK(ILINK)%FAO(:)=FA(:,IATQ)
        LINK(ILINK)%FAI(:)=FA(:,IATQJ)
        LINK(ILINK)%FSO(:)=FS(:,IATS)
        LINK(ILINK)%FSI(:)=FS(:,IATSJ)
        IF(TADIABATIC) THEN
          FORCE(:,IATQ) =FORCE(:,IATQ)+FA(:,IATQ)/ALPHA-FS(:,IATS)
          FORCE(:,IATQJ)=FORCE(:,IATQJ)+FA(:,IATQJ)+FA(:,IATQ)*(1.D0-1.D0/ALPHA)-FS(:,IATSJ)
PRINT*,'FORCE ALONG LINK',DOT_PRODUCT(FORCE(:,IATQ)-FORCE(:,IATQJ),POS(:,IATQ)-POS(:,IATQJ)) &
  /SQRT(SUM((POS(:,IATQ)-POS(:,IATQJ))**2)),SQRT(SUM((POS(:,IATQ)-POS(:,IATQJ))**2)),ILINK,ALPHA
        ELSE    ! FORCES ARE TAKEN CARE OF IN QMMM$PROPAGATE
          FORCE(:,IATQ) =0.D0
          FORCE(:,IATQJ)=0.D0
        END IF
        IND=LINK(ILINK)%SHARED
!CHECK POTENTIALS FOR SHARED BOND MORE CAREFULLY. HERE WE ASSUME
!THAT ONLY THE JOINT ATOM IS SHARED
        IF(IND.EQ.0) THEN
          POT(IATQJ)=POT(IATQJ)+VA(IATQJ)-VS(IATSJ)          
        END IF
        POT(IATQ)=POT(IATQ)+ALPHA*VA(IATQ)+(1.D0-ALPHA)*VA(IATQJ)-VS(IATS)
      ENDDO
!PRINT*,'VA',VA
!PRINT*,'VS',VS
!PRINT*,'POT',POT
!STOP
!
!     ==================================================================
!     ==  DEALLOCATE ARRAYS                                           ==
!     ==================================================================
      DEALLOCATE(MPOS)
      DEALLOCATE(MPOSM)
      DEALLOCATE(MCHARGE)
      DEALLOCATE(MFORCE)
      DEALLOCATE(MPOT)
      DEALLOCATE(SPOS)
      DEALLOCATE(SPOSM)
      DEALLOCATE(SCHARGE)
!      DEALLOCATE(SFORCE)
!      DEALLOCATE(SPOT)
      DEALLOCATE(FA)
      DEALLOCATE(FS)
      DEALLOCATE(VA)
      DEALLOCATE(VS)
!
!     ==================================================================
!     ==  PRINTOUT                                                    ==
!     ==================================================================
      IF(TPR) THEN
        WRITE(*,FMT='("QM-MM TOTAL ENERGY",F10.5)')DEPOT
        DO IATQ=1,NAT
          WRITE(*,FMT='("IAT ",I3,"|R|",F10.5," R ",3F10.5," CHA ",F10.5)') &
     &         IATQ,0.D0,POS(:,IATQ),CHARGE(IATQ)
          SVAR=SQRT(DOT_PRODUCT(FORCE(:,IATQ),FORCE(:,IATQ)))
          WRITE(*,FMT='("IAT ",I3,"|F|",F10.5," F ",3F10.5," POT ",F10.5)') &
     &         IATQ,SVAR,FORCE(:,IATQ),POT(IATQ)
        ENDDO
      ENDIF
                        CALL TIMING$CLOCKOFF('QM-MM')
                        CALL TRACE$POP
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE QMMM$PROPAGATE(NAT_,DELT,ANNER,MASS,ANNERVEC,RM,R0,RP)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NAT_
      REAL(8)     ,INTENT(IN) :: DELT
      REAL(8)     ,INTENT(IN) :: ANNER
      REAL(8)     ,INTENT(IN) :: MASS(NAT_)  ! REDUCED MASS
      REAL(8)     ,INTENT(IN) :: ANNERVEC(NAT_)  
      REAL(8)     ,INTENT(IN) :: RM(3,NAT_)  
      REAL(8)     ,INTENT(INOUT) :: RP(3,NAT_)
      REAL(8)     ,INTENT(IN) :: R0(3,NAT_)
      REAL(8)     ,ALLOCATABLE:: RAP(:,:)
      REAL(8)     ,ALLOCATABLE:: MA(:)
      REAL(8)     ,ALLOCATABLE:: RSP(:,:)
      REAL(8)     ,ALLOCATABLE:: MS(:)
      INTEGER(4)              :: ILINK,I
      INTEGER(4)              :: IATQO,IATQI,IATAO,IATAI,IATSO,IATSI
      REAL(8)                 :: RAOM(3),RAO0(3),RAOP(3)
      REAL(8)                 :: RAIM(3),RAI0(3),RAIP(3)
      REAL(8)                 :: RSOM(3),RSO0(3),RSOP(3)
      REAL(8)                 :: RSIM(3),RSI0(3),RSIP(3)
      REAL(8)                 :: FAO(3),FAI(3)
      REAL(8)                 :: FSO(3),FSI(3)
      REAL(8)                 :: CHIQO,CHIQI,CHIAO,CHIAI,CHISO,CHISI
      REAL(8)                 :: ALPHA
      REAL(8)                 :: MAT(4,4),MATINV(4,4),VEC(4),LAGR(4)
      REAL(8)                 :: DR(3)
      REAL(8)                 :: SVAR,SVAR1,SVAR2,SVAR3
!     ******************************************************************
      IF(.NOT.TON) RETURN
      IF(TADIABATIC) RETURN
                              CALL TRACE$PUSH('QMMM$PROPAGATE')
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$GETI4('NAT',NATM)
      CALL CLASSICAL$SELECT('SHADOW')
      CALL CLASSICAL$GETI4('NAT',NATS)
      CALL ATOMLIST$NATOM(NATQ)
IF(NAT_.NE.NATQ) THEN
  CALL ERROR$STOP('QMMM$PROPAGATE')
END IF
!
      ALLOCATE(MA(NATM))
      ALLOCATE(RAP(3,NATM))
      ALLOCATE(MS(NATS))
      ALLOCATE(RSP(3,NATS))
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$GETR8A('R(+)',3*NATM,RAP)
      CALL CLASSICAL$GETR8A('MASS',NATM,MA)
      CALL CLASSICAL$SELECT('SHADOW')
      CALL CLASSICAL$GETR8A('R(+)',3*NATS,RSP)
      CALL CLASSICAL$GETR8A('MASS',NATS,MS)
      DO ILINK=1,NLINK
        IATQI=LINK(ILINK)%QJOINT
        IATAI=LINK(ILINK)%MJOINT
        IATSI=LINK(ILINK)%SJOINT
        IATQO=LINK(ILINK)%QATOM
        IATAO=LINK(ILINK)%MATOM
        IATSO=LINK(ILINK)%SATOM
        ALPHA=LINK(ILINK)%ALPHA
!       == PROPAGATE MM ATOMS
        SVAR1=2.D0/(1.D0+ANNER)
        SVAR2=1.D0-SVAR1
        SVAR3=DELT**2/(1.D0+ANNER)
        RAO0(:)=R0(:,IATQI)+(R0(:,IATQO)-R0(:,IATQI))/ALPHA
        RAI0(:)=R0(:,IATQI)
        RAOM(:)=RM(:,IATQI)+(RM(:,IATQO)-RM(:,IATQI))/ALPHA
        RAIM(:)=RM(:,IATQI)
        FAO(:)=LINK(ILINK)%FAO(:)
        FAI(:)=LINK(ILINK)%FAI(:)
        RAOP(:)=RAO0(:)*SVAR1+RAOM(:)*SVAR2+FAO(:)*SVAR3/MA(IATAO)
        RAIP(:)=RAI0(:)*SVAR1+RAIM(:)*SVAR2+FAI(:)*SVAR3/MA(IATAI)
!       == PROPAGATE SHADOW ATOMS
        RSO0(:)=R0(:,IATQO)
        RSI0(:)=R0(:,IATQI)
        RSOM(:)=RM(:,IATQO)
        RSIM(:)=RM(:,IATQI)
        FSO(:)=LINK(ILINK)%FSO(:)
        FSI(:)=LINK(ILINK)%FSI(:)
        RSOP(:)=RSO0(:)*SVAR1+RSOM(:)*SVAR2+FSO(:)*SVAR3/MS(IATSO)
        RSIP(:)=RSI0(:)*SVAR1+RSIM(:)*SVAR2+FSI(:)*SVAR3/MS(IATSI)
!
!       ==  ENFORCE CONSTRAINTS ===================================
!       ==  G_?=RQI-RAI
!       ==  G_?=RSI-RAI
!       ==  G_1=RAI+ALPHA*(RAO-RAI)-RQO
!       ==  G_?=RSO-RQO
        CHIQO=1.D0/(MASS(IATQO)*(1.D0+ANNERVEC(IATQO)))
        CHIQI=1.D0/(MASS(IATQI)*(1.D0+ANNERVEC(IATQI)))
        CHIAO=1.D0/(MA(IATAO)*(1.D0+ANNER))
        CHIAI=1.D0/(MA(IATAI)*(1.D0+ANNER))
        CHISO=1.D0/(MS(IATSO)*(1.D0+ANNER))
        CHISI=1.D0/(MS(IATSI)*(1.D0+ANNER))
        MAT(:,:)=0.D0
        MAT(1,1)=(1.D0-ALPHA)*CHIAI ! DG1/DRAI
        MAT(1,2)=CHIQI+CHIAI        ! 
        MAT(1,4)=CHIQI
        MAT(2,2)=CHIQI
        MAT(2,4)=CHIQI-CHISI
        MAT(3,1)=CHIQO
        MAT(3,3)=CHIQO-CHISO
        MAT(4,1)=CHIQO+(1.D0-ALPHA)**2*CHIAI+ALPHA**2*CHIAO
        MAT(4,2)=(1.D0-ALPHA)*CHIAI
        MAT(4,3)=CHIQO
        CALL LIB$INVERTR8(4,MAT,MATINV)
        DO I=1,3
          VEC(1)=RAIP(I)-RP(I,IATQI)
          VEC(2)=RSIP(I)-RP(I,IATQI)    
          VEC(3)=RSOP(I)-RP(I,IATQO)    
          VEC(4)=(1.D0-ALPHA)*RAIP(I)+ALPHA*RAOP(I)-RP(I,IATQO)    
          LAGR(:)=MATMUL(MATINV,VEC)
          RP(I,IATQO)=RP(I,IATQO)+CHIQO*(LAGR(1)+LAGR(3))
          RP(I,IATQI)=RP(I,IATQI)+CHIQI*(LAGR(2)+LAGR(4))
          RAOP(I)=RAOP(I)+CHIAO*(-ALPHA*LAGR(1))
          RAIP(I)=RAIP(I)+CHIAI*(-(1.D0-ALPHA)*LAGR(1)-LAGR(2))
          RSOP(I)=RSOP(I)+CHISO*LAGR(3)
          RSIP(I)=RSIP(I)+CHISI*LAGR(4)
        ENDDO
        SVAR=0.D0
        DR(:)=RP(:,IATQO)-RSOP(:)
        SVAR=SVAR+SUM(DR**2)
        DR(:)=RP(:,IATQI)-RSIP(:)
        SVAR=SVAR+SUM(DR**2)
        DR(:)=RP(:,IATQO)-( RAIP(:)+(RAOP(:)-RAIP(:))*ALPHA )
        SVAR=SVAR+SUM(DR**2)
        DR(:)=RP(:,IATQI)-RAIP(:)
        SVAR=SVAR+SUM(DR**2)
        SVAR=SQRT(SVAR/12.D0)
        IF(SVAR.GT.1.D-7) THEN
PRINT*,'QOP',IATQO,RP(:,IATQO)
PRINT*,'SOP',IATSO,RSOP
PRINT*,'MOP',IATAO,RAOP
PRINT*,'QIP',IATQI,RP(:,IATQI)
PRINT*,'SIP',IATSI,RSIP
PRINT*,'MIP',IATAI,RAIP
          CALL ERROR$MSG('CONSTRAINTS NOT FULFILLED FOR LINK BONDS')
          CALL ERROR$I4VAL('ILINK',ILINK)
          CALL ERROR$R8VAL('SIGMA',SVAR)
          CALL ERROR$STOP('QMMM$PROPAGATE')
        END IF
      ENDDO
      DEALLOCATE(RAP)
      DEALLOCATE(RSP)
      DEALLOCATE(MA)
      DEALLOCATE(MS)
                              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE QMMM$DEKIN(DELT,EKIN)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE KINETIC ENERGY OF THE LINK ATOMS                  **
!     **                                                              **
!     **  THE KINETIC ENERGY OF THE ENVIRONMENT ATOMS, EXCLUDING      **
!     **  ANY LINK ATOMS, ARE CALCULATED IN $INTERFACE                **
!     **                                                              **
!     ******************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
      REAL(8)  ,INTENT(IN)  :: DELT  ! TIME STEP OF THE QM SYSTEM
      REAL(8)  ,INTENT(OUT) :: EKIN  ! KINETIC ENERGY CORRECTION
      REAL(8)               :: SVAR1,SVAR2
      REAL(8)               :: QRM(3,NATQ)
      REAL(8)               :: QRP(3,NATQ)
      REAL(8)               :: QMASS(NATQ)
      REAL(8)               :: MRP(3)
      REAL(8)               :: MRM(3)
      REAL(8)               :: MMASS(NATM)
      REAL(8)               :: SMASS(NATS)
      INTEGER(4)            :: ILINK,I
      INTEGER(4)            :: IATM,IATQ,IATQJ,IATS
      REAL(8)               :: ALPHA
!     *****************************************************************
      EKIN=0.D0
      IF (.NOT. TON) RETURN
                              CALL TRACE$PUSH('QMMM$EKIN')
!
      CALL ATOMLIST$GETR8A('R(+)',0,3*NATQ,QRP)
      CALL ATOMLIST$GETR8A('R(-)',0,3*NATQ,QRM)
      CALL ATOMLIST$GETR8A('MASS',0,NATQ,QMASS)
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$GETR8A('MASS',NATM,MMASS)
      CALL CLASSICAL$SELECT('SHADOW')
      CALL CLASSICAL$GETR8A('MASS',NATS,SMASS)
      DO ILINK=1,NLINK
        IATQ =LINK(ILINK)%QATOM
        IATQJ=LINK(ILINK)%QJOINT
        IATM =LINK(ILINK)%MATOM
        IATS =LINK(ILINK)%SATOM
        ALPHA=LINK(ILINK)%ALPHA
        MRP(:)=QRP(:,IATQJ)+(QRP(:,IATQ)-QRP(:,IATQJ))/ALPHA
        MRM(:)=QRM(:,IATQJ)+(QRM(:,IATQ)-QRM(:,IATQJ))/ALPHA
        SVAR1=0.D0
        SVAR2=0.D0
        DO I=1,3
          SVAR1=SVAR1+(MRP(I)-MRM(I))**2
          SVAR2=SVAR2+(QRP(I,IATQ)-QRM(I,IATQ))**2
        ENDDO
        SVAR1=SVAR1+0.5D0*MMASS(IATM)*SVAR1/(2.D0*DELT)**2
        SVAR2=SVAR2+0.5D0*SMASS(IATS)*SVAR2/(2.D0*DELT)**2
        EKIN=EKIN+SVAR1-SVAR2
      ENDDO
                              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE QMMM$SWITCH
!     ******************************************************************
!     **                                                              **
!     **  R(0)->R(-);   R(+)->R(0)                                    **
!     **                                                              **
!     ******************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
!     *****************************************************************
      IF(.NOT.TON) RETURN
      EKIN_QMMM=0.D0
      EPOT_QMMM=0.D0
      CALL CLASSICAL$SELECT('SHADOW')
      CALL CLASSICAL$SWITCH
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$SWITCH
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE QMMM_RANDOMIZEVELOCITY
!     ******************************************************************
!     **                                                              **
!     **  RANDOMIZE VELOCITY ACCORDING TO A GIVEN TEMPERATURE         **
!     **                                                              **
!     **    TEMP     K_B*T WHERE THE T IS THE TEMPERATURE OF THE      **
!     **             HEATBATH                                         **
!     **    THE TARGET KINETIC ENERGY IS 1.5*(NATM-NATS)*TEMP         **
!     **                                                              **
!     ******************************************************************
      USE QMMM_MODULE,ONLY  : DELTAT,AMPRE,NATM,NATS &
     &                       ,NMAP,MAP,NLINK,LINK
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
      REAL(8)               :: RM(3,NATM)
      REAL(8)               :: RMOLD(3,NATM)
      REAL(8)               :: MASS(NATM)
      INTEGER(4)            :: IAT,I,IMAP,ILINK
      REAL(8)               :: SVAR,RAN,SUM
      REAL(8)               :: KELVIN
      INTEGER(4)            :: NFILO
!     ******************************************************************
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$GETR8A('R(-)',3*NATM,RM)
      RMOLD(:,:)=RM(:,:)
      CALL CLASSICAL$GETR8A('MASS',NATM,MASS)
!
!     ==================================================================
!     ==  CHANGE VELOCITIES                                           ==
!     ==================================================================
      DO IAT=1,NATM
        SVAR=SQRT(AMPRE/MASS(IAT))*DELTAT
        DO I=1,3
          CALL GAUSS_RANDOM_NUMBER(RAN)
          RM(I,IAT) = RM(I,IAT) + SVAR*RAN
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  LEAVE THE VELOCITIES OF THE REACTION CENTER UNCHANGED       ==
!     ==================================================================
      DO IMAP=1,NMAP
        IAT=MAP(IMAP)%MATOM
        RM(:,IAT)=RMOLD(:,IAT)
      ENDDO
      DO ILINK=1,NLINK
        IAT=LINK(ILINK)%MATOM
        RM(:,IAT)=RMOLD(:,IAT)
      ENDDO
!
!     ==================================================================
!     ==  HAND THE MODIFIED R(-) BACK                                 ==
!     ==================================================================
      CALL CLASSICAL$SETR8A('R(-)',3*NATM,RM)
!
!     ==================================================================
!     ==  PRINTOUT FOR TESTS                                          ==
!     ==================================================================
      IF(TPR) THEN
        SUM=0.D0
        DO IAT=1,NATM
          DO I=1,3
            SUM=SUM+0.5D0*MASS(IAT)*((RM(I,IAT)-RMOLD(I,IAT))/DELTAT)**2
          ENDDO
        ENDDO
        CALL CONSTANTS('KB',KELVIN)
        SUM=SUM/(0.5D0*DBLE(NATM-NATS)*KELVIN)
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,*)'QMMM-RANDOMIZED: T=',SUM,' KELVIN'
      ENDIF
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE QMMM_MINIMIZE(TCONV)
!     ******************************************************************
!     **                                                              **
!     **  OPTIMISES ENVIRONMENT ATOMS, EXCLUDING ATOMS OF THE REACTION**
!     **  CENTER AND LINK ATOMS, USING A CONJUGATE GRADIENT MINIMIZER **
!     **                                                              **
!     ******************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(OUT):: TCONV
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
      INTEGER(4),PARAMETER  :: NITER=10000
      REAL(8)   ,PARAMETER  :: TOL=1.D-3
      REAL(8)               :: FORCE1(3*NATM)
      REAL(8)               :: FORCE2(3*NATM)
      REAL(8)               :: DIR1(3*NATM)
      REAL(8)               :: DIR2(3*NATM)
      REAL(8)               :: R1(3*NATM)
      REAL(8)               :: R2(3*NATM)
      INTEGER(4)            :: ITER,INNER
      REAL(8)               :: EPOT1,EPOT2,FMAX
      REAL(8)               :: ALPHA,ALPHALAST,DALPHA
      REAL(8)               :: EPOTLAST,FORCELAST(3*NATM)
      INTEGER(4)            :: NFILO
      REAL(8)               :: SVAR1,SVAR2
!     ******************************************************************
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ================================================================
!     ==  PREPARE INITIAL STEP                                      ==
!     ================================================================
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$GETR8A('R(0)',3*NATM,R1)
      CALL CLASSICAL$NEIGHBORS
      CALL ONETOT(NATM,R1,EPOT1,FORCE1)
      DIR1=FORCE1
!
!     ==================================================================
!     ==  ITERATE                                                     ==
!     ==================================================================
      DO ITER=1,NITER
        IF(MOD(ITER,100).EQ.0)CALL CLASSICAL$NEIGHBORS
!
!       ================================================================
!       ==  PRINT                                                     ==
!       ================================================================
        IF(TPR) THEN
          FMAX=SQRT(DOT_PRODUCT(FORCE1,FORCE1))
          CALL FILEHANDLER$UNIT('PROT',NFILO)
!          WRITE(NFILO,FMT='("I",I10," EPOT ",E17.10," FMAX ",E12.5)')ITER,EPOT1,FMAX
          WRITE(*,FMT='("I",I10," EPOT ",E17.10," FMAX ",E12.5)')ITER,EPOT1,FMAX
        END IF
!
!       ================================================================
!       ==  CHECK CONVERGENCE                                         ==
!       ================================================================
        FMAX=SQRT(DOT_PRODUCT(FORCE1,FORCE1))
        TCONV=(FMAX.LT.TOL)
        IF(TCONV) EXIT
!
!       ================================================================
!       ==  PERFORM CG LINE SEARCH  ALONG DIRECTION FORCE1            ==
!       ================================================================
        CALL QMMM_CG_LINESEARCH(NATM,R1,DIR1,ALPHA)
        R2=R1+ALPHA*DIR1

!!$!
!!$        ALPHA=1.D-3/SQRT(SUM(DIR1**2))
!!$        ALPHALAST=0.D0
!!$        EPOTLAST=EPOT1
!!$        FORCELAST=FORCE1
!!$        DO INNER=1,10000
!!$          R2(:)=R1(:)+ALPHA*DIR1(:)
!!$          CALL ONETOT(NATM,R2,EPOT2,FORCE2)
!!$!
!!$!         == TRY AGAIN WITH SMALLER STEP IF ENERGY WENT UP ==============
!!$          IF(EPOT2.GT.EPOTLAST) THEN
!!$            ALPHA=0.5D0*(ALPHALAST+ALPHA)
!!$PRINT*,'ENERGY GOES UP; TRY AGAIN: ',EPOT2,ALPHA,ALPHA/SQRT(SUM(DIR1**2))
!!$            CYCLE
!!$          END IF
!!$!
!!$!         == REPORT ACCEPTED MOVE ========================================
!!$          IF(TPR) THEN
!!$!IF(INNER.GT.0) THEN
!!$            WRITE(NFILO,FMT='("INNER",I6," EPOT ",E12.5," FMAX ",E12.5," ALPHA ",F10.5)') &
!!$       &                    INNER,EPOT2,DOT_PRODUCT(DIR1,FORCE2),ALPHA
!!$          END IF
!!$!
!!$!         == ENERGY WENT DOWN, CALCULATE NEW ALPHA ====================
!!$          SVAR1=DOT_PRODUCT(DIR1,FORCE2)
!!$          SVAR2=DOT_PRODUCT(DIR1,FORCE2-FORCELAST)
!!$          TCONV=(ABS(SVAR1).LT.1.D-5)
!!$          IF(TCONV) EXIT
!!$!
!!$          DALPHA=ALPHA-ALPHALAST          
!!$          IF(-DALPHA*SVAR2.GT.0.D0) THEN   ! CURVATURE POSITIVE
!!$            DALPHA=-SVAR1/SVAR2*DALPHA
!!$PRINT*,'CURVATURE POSITIVE ',DALPHA,ALPHA,ALPHALAST,EPOT2,ITER,INNER
!!$          ELSE                             ! STEEPEST DESCENT FOR NEGATIVE CURVATURE
!!$            DALPHA=SIGN(1.D-2/DOT_PRODUCT(DIR1,DIR1),SVAR2)
!!$PRINT*,'CURVATURE NEGATIVE ',DALPHA,ALPHA,ALPHALAST,EPOT2,ITER,INNER
!!$          END IF
!!$!
!!$!         == STORE SUCCESSFUL MOVE AS REFERENCE AND DETERMINE NEW ALPHA ==
!!$          ALPHALAST=ALPHA
!!$          EPOTLAST=EPOT2
!!$          FORCELAST=FORCE2
!!$          ALPHA=ALPHA+DALPHA
!!$        ENDDO
!!$        IF(.NOT.TCONV) THEN
!!$          CALL ERROR$MSG('LINE SEARCH NOT CONVERGED')
!!$          CALL ERROR$STOP('QMMM_MINIMIZE')
!!$        END IF
        CALL ONETOT(NATM,R2,EPOT2,FORCE2)
!
!       == CHOOSE NEW SEARCH DIRECTION ===================================
        CALL CG$NEWDIR(3*NATM,FORCE1,DIR1,FORCE2,DIR2)
        R1=R2
        FORCE1=FORCE2
        DIR1=DIR2
        EPOT1=EPOT2
      ENDDO
!     == END OF  CONVERGENCE LLOOP
      IF(.NOT.TCONV) THEN
        CALL ERROR$MSG('LOOP NOT CONVERGED')        
        CALL ERROR$R8VAL('FMAX',FMAX)        
        CALL ERROR$STOP('QMMM_MINIMIZE')
      END IF
      PRINT*,'QMMM$MINIMIZE FINISHED AFTER ',ITER,' ITERATIONS'
      CALL CLASSICAL$SETR8A('R(0)',3*NATM,R2)
      CALL CLASSICAL$ETOT(EPOT2)
      RETURN
      CONTAINS
!       .................................................................
        SUBROUTINE ONETOT(NAT,R,E,F)
!       USE QMMM_MODULE ,ONLY : MAP,NMAP,LINK,NLINK,MAP_TYPE,LINK_TYPE
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: NAT
        REAL(8)   ,INTENT(IN) :: R(3,NAT)
        REAL(8)   ,INTENT(OUT):: E
        REAL(8)   ,INTENT(OUT):: F(3,NAT)
        INTEGER(4)            :: IMAP,ILINK,IAT
!       *****************************************************************
        CALL CLASSICAL$SETR8A('R(0)',3*NAT,R)
CALL CLASSICAL$NEIGHBORS
        CALL CLASSICAL$ETOT(E)
        CALL CLASSICAL$GETR8A('FORCE',3*NAT,F)
        DO IMAP=1,NMAP
          IAT=MAP(IMAP)%MATOM
          F(:,IAT)=0.D0
        ENDDO
        DO ILINK=1,NLINK
          IAT=LINK(ILINK)%MATOM
          F(:,IAT)=0.D0
        ENDDO
        END SUBROUTINE ONETOT
      END SUBROUTINE QMMM_MINIMIZE

!
!     ................................................................
      SUBROUTINE QMMM_CG_LINESEARCH(NAT,RZERO,DIR,ALPHA)
!     *****************************************************************
!     **                                                             **
!     **                                                             **
!     *****************************************************************
      USE QMMM_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RZERO(3,NAT)
      REAL(8)   ,INTENT(IN) :: DIR(3,NAT)
      REAL(8)   ,INTENT(OUT):: ALPHA
      REAL(8)               :: DIRLEN
      REAL(8)               :: ALPHAM,ALPHA0,ALPHAP
      REAL(8)               :: EM,E0
      REAL(8)               :: FM,F0
      REAL(8)               :: R(3,NAT)
      REAL(8)               :: FORCE(3,NAT)
      INTEGER(4)            :: ITER
      LOGICAL(4)            :: TCONV
      LOGICAL(4)            :: TPR=.TRUE.
      INTEGER(4),PARAMETER  :: NITER=10000
      REAL(8)   ,PARAMETER  :: TOL=1.D-4
      REAL(8)   ,PARAMETER  :: FIRSTSTEPSIZE=1.D-1
      REAL(8)               :: CURVATURE,LASTCURVATURE
      LOGICAL(4)            :: THARMONIC
REAL(8),ALLOCATABLE :: EHISTORY(:) !(NITER)
REAL(8),ALLOCATABLE :: FHISTORY(:) !(NITER)
REAL(8),ALLOCATABLE :: AHISTORY(:) !(NITER)
!     *****************************************************************
ALLOCATE(EHISTORY(NITER))
ALLOCATE(FHISTORY(NITER))
ALLOCATE(AHISTORY(NITER))
EHISTORY=0.D0
FHISTORY=0.D0
AHISTORY=0.D0
      DIRLEN=SQRT(SUM(DIR**2))
      LASTCURVATURE=1.D+10
      ALPHAM=0.D0
      R=RZERO
      CALL ONETOT1(NAT,R,EM,FORCE)
      FM=SUM(DIR*FORCE)
      ALPHA0=SIGN(FIRSTSTEPSIZE/DIRLEN,FM)
      DO ITER=1,NITER
        R=RZERO+ALPHA0*DIR
        CALL ONETOT1(NAT,R,E0,FORCE)
        F0=SUM(DIR*FORCE)
EHISTORY(ITER)=E0
FHISTORY(ITER)=F0/DIRLEN
AHISTORY(ITER)=ALPHA0*DIRLEN
!
!       == EXCEPTION IF ENERGY GOES UP
        IF(E0.GT.EM) THEN
PRINT*,'ENERGY GOES UP ',E0,E0-EM,ALPHA0,ALPHA0-ALPHAM,F0,FM
          ALPHA0=0.5D0*(ALPHAM+ALPHA0)
          CYCLE
        END IF

        IF(TPR) THEN
          WRITE(*,FMT='("INNER",I6," EPOT ",E12.5," FMAX ",E12.5," ALPHA ",F10.5)') &
       &                  ITER,E0,F0/DIRLEN,ALPHA0
        END IF
!!
!       == CHECK CONVERGENCE =========================================
        TCONV=(ABS(F0/DIRLEN).LT.TOL)
        IF(TCONV) EXIT
!
!       == DETERMINE NEW MIXING FACTOR ===============================
!       == CHECK IF SYSTEM IS IN THE HARMONIC REGIME =================
!       == IF NOT OR IF CURVATURE IS NEGATIVE, SWITCH TO STEPPING ====
        CURVATURE=-(F0-FM)/(ALPHA0-ALPHAM)/DIRLEN
        THARMONIC=(ABS(CURVATURE/LASTCURVATURE-1.D0).LT.0.2)
        IF(THARMONIC.AND.CURVATURE.GT.0.D0) THEN
          ALPHAP=(ALPHA0-ALPHAM)/(F0-FM)
PRINT*,'==',E0,F0,CURVATURE
          ALPHAP=ALPHAM-FM*ALPHAP
        ELSE
PRINT*,'==',E0,F0,CURVATURE,' STEPPING'
          ALPHAP=ALPHA0+SIGN(FIRSTSTEPSIZE,F0)
        END IF
!
!       == SWITCH  ===================================================
        FM=F0
        EM=E0
        ALPHAM=ALPHA0
        ALPHA0=ALPHAP
        LASTCURVATURE=CURVATURE
      ENDDO
      IF(.NOT.TCONV) THEN
OPEN(8,FILE='DUMP')
REWIND 8
DO ITER=1,NITER
WRITE(8,*)ITER,AHISTORY(ITER),EHISTORY(ITER),FHISTORY(ITER)
ENDDO
CLOSE(8)
OPEN(8,FILE='DUMP2')
REWIND 8
DO ITER=1,NITER
WRITE(8,*)AHISTORY(ITER),EHISTORY(ITER),FHISTORY(ITER)
ENDDO
CLOSE(8)
        CALL ERROR$MSG('CG LINE SEARCH NOT CONVERGED')
        CALL ERROR$STOP('QMMM_CG_LINESEARCH')
      END IF
      ALPHA=ALPHA0
DEALLOCATE(EHISTORY)
DEALLOCATE(FHISTORY)
DEALLOCATE(AHISTORY)
      RETURN
      CONTAINS
!       .................................................................
        SUBROUTINE ONETOT1(NAT,R,E,F)
!        USE QMMM_MODULE ,ONLY : MAP,NMAP,LINK,NLINK,MAP_TYPE,LINK_TYPE
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: NAT
        REAL(8)   ,INTENT(IN) :: R(3,NAT)
        REAL(8)   ,INTENT(OUT):: E
        REAL(8)   ,INTENT(OUT):: F(3,NAT)
        INTEGER(4)            :: IMAP,ILINK,IAT
!       *****************************************************************
        CALL CLASSICAL$SETR8A('R(0)',3*NAT,R)
        CALL CLASSICAL$ETOT(E)
        CALL CLASSICAL$GETR8A('FORCE',3*NAT,F)
        DO IMAP=1,NMAP
          IAT=MAP(IMAP)%MATOM
          F(:,IAT)=0.D0
        ENDDO
        DO ILINK=1,NLINK
          IAT=LINK(ILINK)%MATOM
          F(:,IAT)=0.D0
        ENDDO
        END SUBROUTINE ONETOT1
      END SUBROUTINE QMMM_CG_LINESEARCH
!
!     ................................................................
      SUBROUTINE QMMM$WRITE(NFIL,NFILO,TCHK)
!     *****************************************************************
!     **                                                             **
!     **                                                             **
!     *****************************************************************
      USE QMMM_MODULE
      USE RESTART_INTERFACE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN) :: NFIL
      INTEGER(4)            ,INTENT(IN) :: NFILO
      LOGICAL(4)            ,INTENT(OUT):: TCHK
      TYPE (SEPARATOR_TYPE),PARAMETER  :: MYSEPARATOR &
                 =SEPARATOR_TYPE(8,'QM-MM','NONE','AUG1996','NONE')
      TYPE (SEPARATOR_TYPE)            :: SEPARATOR
      INTEGER(4)                        :: NAT
      REAL(8)               ,ALLOCATABLE:: R0(:,:)
      REAL(8)               ,ALLOCATABLE:: RM(:,:)
      REAL(8)               ,ALLOCATABLE:: QEL(:)
      INTEGER(4)                        :: NTASKS,THISTASK
!     ******************************************************************
      TCHK=.FALSE.
      IF(.NOT.TON) RETURN
!
!     == WRITE SHADOW ======
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      CALL CLASSICAL$SELECT('SHADOW')
      IF(THISTASK.EQ.1) THEN
        CALL RESTART$WRITESEPARATOR(MYSEPARATOR,NFIL,NFILO,TCHK)
        CALL CLASSICAL$GETI4('NAT',NAT)
        ALLOCATE(R0(3,NAT))
        ALLOCATE(RM(3,NAT))
        ALLOCATE(QEL(NAT))
        CALL CLASSICAL$GETR8A('R(0)',3*NAT,R0)
        CALL CLASSICAL$GETR8A('R(-)',3*NAT,RM)
        CALL CLASSICAL$GETR8A('QEL',NAT,QEL)
!
        WRITE(NFIL)NAT
        WRITE(NFIL)R0(:,:)
        WRITE(NFIL)RM(:,:)
        WRITE(NFIL)QEL(:)
!
        DEALLOCATE(R0)
        DEALLOCATE(RM)
        DEALLOCATE(QEL)
      END IF
!
!     == WRITE QMMM ======
      CALL CLASSICAL$SELECT('QMMM')
      IF(THISTASK.EQ.1) THEN
        CALL CLASSICAL$GETI4('NAT',NAT)
        ALLOCATE(R0(3,NAT))
        ALLOCATE(RM(3,NAT))
        ALLOCATE(QEL(NAT))
        CALL CLASSICAL$GETR8A('R(0)',3*NAT,R0)
        CALL CLASSICAL$GETR8A('R(-)',3*NAT,RM)
        CALL CLASSICAL$GETR8A('QEL',NAT,QEL)
        WRITE(NFIL)NAT
        WRITE(NFIL)R0(:,:)
        WRITE(NFIL)RM(:,:)
        WRITE(NFIL)QEL(:)
        DEALLOCATE(R0)
        DEALLOCATE(RM)
        DEALLOCATE(QEL)
      END IF
      RETURN
      END
!
!     ................................................................
      SUBROUTINE QMMM$READ(NFIL,NFILO,TCHK)
!     *****************************************************************
!     **                                                             **
!     **                                                             **
!     *****************************************************************
      USE QMMM_MODULE
      USE RESTART_INTERFACE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)           ,INTENT(IN) :: NFIL
      INTEGER(4)           ,INTENT(IN) :: NFILO
      LOGICAL(4)           ,INTENT(OUT):: TCHK
      TYPE (SEPARATOR_TYPE),PARAMETER  :: MYSEPARATOR &
                 =SEPARATOR_TYPE(8,'QM-MM','NONE','AUG1996','NONE')
      TYPE (SEPARATOR_TYPE)            :: SEPARATOR
      INTEGER(4)                       :: NTASKS,THISTASK
      INTEGER(4)                       :: NAT
      REAL(8)              ,ALLOCATABLE:: R0(:,:)
      REAL(8)              ,ALLOCATABLE:: RM(:,:)
      REAL(8)              ,ALLOCATABLE:: QEL(:)
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      TCHK=TON
      SEPARATOR=MYSEPARATOR
      IF(THISTASK.EQ.1)CALL RESTART$READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      CALL MPE$BROADCAST('MONOMER',1,TCHK)
      IF(.NOT.TCHK) RETURN

      IF(SEPARATOR%VERSION.NE.MYSEPARATOR%VERSION) THEN
        CALL ERROR$MSG('VERSION NOT CONSISTENT')
        CALL ERROR$STOP('QMMM$READ')
      END IF
      CALL QMMM_INITIALIZE
!
!     == READ SHADOW ======
      CALL CLASSICAL$SELECT('SHADOW')
      ALLOCATE(R0(3,NATS))
      ALLOCATE(RM(3,NATS))
      ALLOCATE(QEL(NATS))
      IF(THISTASK.EQ.1) THEN
        READ(NFIL)NAT
        IF(NAT.NE.NATS) THEN
          CALL ERROR$MSG('#(SHADOW ATOMS NOT CONSISTENT')
          CALL ERROR$I4VAL('NATS READ',NAT)
          CALL ERROR$I4VAL('NATS EXPECTED',NATS)
          CALL ERROR$STOP('QM-MM$READ')
        END IF
        READ(NFIL)R0(:,:)
        READ(NFIL)RM(:,:)
        READ(NFIL)QEL(:)
      END IF
      CALL MPE$BROADCAST('MONOMER',1,R0)
      CALL MPE$BROADCAST('MONOMER',1,RM)
      CALL MPE$BROADCAST('MONOMER',1,QEL)
      CALL CLASSICAL$SETR8A('R(0)',3*NATS,R0)
      CALL CLASSICAL$SETR8A('R(-)',3*NATS,RM)
      CALL CLASSICAL$SETR8A('QEL',NATS,QEL)
      DEALLOCATE(R0)
      DEALLOCATE(RM)
      DEALLOCATE(QEL)
!
!     == WRITE QMMM ======
      CALL CLASSICAL$SELECT('QMMM')
      ALLOCATE(R0(3,NATM))
      ALLOCATE(RM(3,NATM))
      ALLOCATE(QEL(NATM))
      IF(THISTASK.EQ.1) THEN
        READ(NFIL)NAT
        IF(NAT.NE.NATM) THEN
          CALL ERROR$MSG('#(QM-MM ATOMS NOT CONSISTENT')
          CALL ERROR$I4VAL('NATM READ',NAT)
          CALL ERROR$I4VAL('NATM EXPECTED',NATM)
          CALL ERROR$STOP('QMMM$READ')
        END IF
        READ(NFIL)R0(:,:)
        READ(NFIL)RM(:,:)
        READ(NFIL)QEL(:)   !CHARGE WILL NOT BE USED
      END IF
      CALL MPE$BROADCAST('MONOMER',1,R0)
      CALL MPE$BROADCAST('MONOMER',1,RM)
!     CALL MPE$BROADCAST('MONOMER',1,QEL)
      CALL CLASSICAL$SETR8A('R(0)',3*NATM,R0)
      CALL CLASSICAL$SETR8A('R(-)',3*NATM,RM)
!     CALL CLASSICAL$SETR8A('QEL',NATM,QEL)
      DEALLOCATE(R0)
      DEALLOCATE(RM)
      DEALLOCATE(QEL)
      RETURN
    END SUBROUTINE QMMM$READ
