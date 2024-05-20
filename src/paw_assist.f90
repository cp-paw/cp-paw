!........................................................AUTO<PI..........
MODULE DIALS_MODULE 
!***********************************************************************
!**                                                                   **
!**  NAME: DIALS                                                      **
!**                                                                   **
!**  PURPOSE: CHANGES SLOWLY THE VALUES OF AN OBJECT                   **
!**                                                                   **
!***********************************************************************
TYPE DIAL_TYPE
  CHARACTER(32) :: NAME
  CHARACTER(32) :: DESCRIPTION
  LOGICAL(4)    :: ON
  REAL(8)       :: ACTUAL    ! ACTUAL VALUE
  REAL(8)       :: FINAL     ! FINAL VALUE
  REAL(8)       :: RATE      ! SIGN IS IRRELEVANT
END TYPE DIAL_TYPE
LOGICAL(4)                  :: TINI=.FALSE.
LOGICAL(4)                  :: TSET=.FALSE.
INTEGER(4)       ,PARAMETER :: NDIALS=4
TYPE (DIAL_TYPE)            :: DIAL(NDIALS) &
     =DIAL_TYPE(" "," ",.FALSE.,0.D0,0.D0,0.D0)
INTEGER(4)                  :: ISELECT=0
REAL(8)                     :: DELTAT=0.D0  !TIME STEP
CONTAINS
!       ................................................................
        SUBROUTINE DIALS_INITIALIZE
        IMPLICIT NONE
        IF(TINI) RETURN
        TINI=.TRUE.
        DIAL(1)%NAME='TEMP(E)'
        DIAL(1)%DESCRIPTION='ELECTRONIC TEMPERATURE'
        DIAL(2)%NAME='TEMP(R)'
        DIAL(2)%DESCRIPTION='TEMPERATURE OF ATOM-THERMOSTAT'
        DIAL(3)%NAME='SPIN'
        DIAL(3)%DESCRIPTION='TOTAL SPIN'
        DIAL(4)%NAME='CHARGE'
        DIAL(4)%DESCRIPTION='TOTAL CHARGE'
        RETURN
        END SUBROUTINE DIALS_INITIALIZE
!       ...............................................................
        SUBROUTINE DIALS_CHECKOUT(NAME,DIAL_)
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN)  :: NAME
        TYPE(DIAL_TYPE) ,INTENT(OUT) :: DIAL_
        INTEGER(4)                    :: I
        DO I=1,NDIALS
          IF(TRIM(NAME).EQ.TRIM(DIAL(I)%NAME)) THEN
            DIAL_=DIAL(I)
            RETURN
          END IF
        ENDDO
        CALL ERROR$MSG('DIAL NAME NOT RECOGNIZED')
        CALL ERROR$CHVAL('NAME',NAME)
        CALL ERROR$STOP('DIALS_CHECKOUT')
        END SUBROUTINE DIALS_CHECKOUT
!       ...............................................................
        SUBROUTINE DIALS_CHECKIN(NAME,DIAL_)
        IMPLICIT NONE
        CHARACTER(*)   ,INTENT(IN)   :: NAME
        TYPE(DIAL_TYPE),INTENT(INOUT):: DIAL_
        INTEGER(4)                   :: I
        DO I=1,NDIALS
          IF(TRIM(NAME).EQ.TRIM(DIAL(I)%NAME)) THEN
            DIAL(I)=DIAL_
            RETURN
          END IF
        ENDDO
        CALL ERROR$MSG('DIAL NAME NOT RECOGNIZED')
        CALL ERROR$CHVAL('NAME',NAME)
        CALL ERROR$STOP('DIALS_CHECKIN')
        END SUBROUTINE DIALS_CHECKIN
!       ...............................................................
        SUBROUTINE DIALS_TURN(DIAL_,DT)
        IMPLICIT NONE
        TYPE(DIAL_TYPE),INTENT(INOUT):: DIAL_
        REAL(8)         ,INTENT(IN)  :: DT     ! TIME STEP
        REAL(8)                      :: CHANGE
        REAL(8)                      :: LEFT
        CHANGE=ABS(DIAL_%RATE*DT)
        LEFT=DIAL_%FINAL-DIAL_%ACTUAL
!       == AVOID OVERSHOOTING  =========================================
        CHANGE=MIN(CHANGE,ABS(LEFT))
!       == DETERMINE DIRECTION OF CHANGE================================
        CHANGE=SIGN(CHANGE,LEFT)
!       == RESET ACTUAL VALUE ==========================================
        DIAL_%ACTUAL=DIAL_%ACTUAL+CHANGE
        RETURN
        END SUBROUTINE DIALS_TURN
      END MODULE DIALS_MODULE
!
!     ...................................................AUTOPI.........
      SUBROUTINE DIALS$SELECT(ID_)
!     ******************************************************************
!     **  REPORT SETTING OF AUTOPILOT                                 **
!     ******************************************************************
      USE DIALS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)                   :: I
!     ******************************************************************
      CALL DIALS_INITIALIZE
      IF(ID_.EQ.'..'.OR.ID_.EQ.'~') THEN
        ISELECT=0
        RETURN
      END IF
      DO I=1,NDIALS
        IF(TRIM(DIAL(I)%NAME).EQ.ID_) THEN
          ISELECT=I
          RETURN
        END IF
      ENDDO
      RETURN
      END
!
!     ...................................................AUTOPI.........
      SUBROUTINE DIALS$SETR8(ID_,VAL_)
!     ******************************************************************
!     **  REPORT SETTING OF AUTOPILOT                                 **
!     ******************************************************************
      USE DIALS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'TIMESTEP') THEN
        DELTAT=VAL_
        RETURN 
      END IF
!
      CALL DIALS_INITIALIZE
      IF(ISELECT.LE.0) THEN
        CALL ERROR$MSG('A DIAL MUST BE SELECTED')
        CALL ERROR$STOP('DIALS$SETR8')
      END IF
      IF(ID_.EQ.'FINAL') THEN
        DIAL(ISELECT)%FINAL=VAL_
      ELSE IF(ID_.EQ.'RATE') THEN
        DIAL(ISELECT)%RATE=VAL_
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIALS$SETR8')
      END IF
      RETURN
      END
!
!     ...................................................AUTOPI.........
      SUBROUTINE DIALS$SETL4(ID_,VAL_)
!     ******************************************************************
!     **  REPORT SETTING OF AUTOPILOT                                 **
!     ******************************************************************
      USE DIALS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL(4)  ,INTENT(IN) :: VAL_
!     ******************************************************************
      CALL DIALS_INITIALIZE
      IF(ISELECT.LE.0) THEN
        CALL ERROR$MSG('A DIAL MUST BE SELECTED')
        CALL ERROR$STOP('DIALS$SETR8')
      END IF
      IF(ID_.EQ.'ON') THEN
        DIAL(ISELECT)%ON=VAL_
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIALS$SETL4')
      END IF
      RETURN
      END
!     
!
!     ...................................................AUTOPI.........
      SUBROUTINE DIALS$REPORT(NFIL)
!     ******************************************************************
!     **  REPORT SETTING OF AUTOPILOT                                 **
!     ******************************************************************
      USE DIALS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: IDIAL 
!     ******************************************************************
      CALL DIALS_INITIALIZE
      DO IDIAL=1,NDIALS
        IF(.NOT.DIAL(IDIAL)%ON) CYCLE
        CALL REPORT$TITLE(NFIL,'DIAL :'//TRIM(DIAL(IDIAL)%NAME))
        CALL REPORT$STRING(NFIL,TRIM(DIAL(IDIAL)%DESCRIPTION)//':')
        CALL REPORT$R8VAL(NFIL,'RATE',DIAL(IDIAL)%RATE,'A.U.')
        CALL REPORT$R8VAL(NFIL,'FINAL VALUE',DIAL(IDIAL)%FINAL,'A.U.')
      ENDDO
      RETURN
      END
!
!     ...................................................AUTOPI.........
      SUBROUTINE DIALS$APPLY
!     ******************************************************************
!     **  REPORT SETTING OF AUTOPILOT                                 **
!     ******************************************************************
      USE DIALS_MODULE 
      IMPLICIT NONE
      TYPE (DIAL_TYPE) :: DIAL1
      CHARACTER(32)    :: NAME
      REAL(8)          :: VALUE
!     ******************************************************************
      CALL DIALS_INITIALIZE
!
!     =================================================================
!     ==  SET INITIAL VALUES                                         ==
!     =================================================================
      IF(.NOT.TSET) THEN
!
!     ==  TEMPERATURE OF THE ELECTRONS ================================
        NAME='TEMP(E)'
        CALL DIALS_CHECKOUT(NAME,DIAL1)
        CALL DYNOCC$GETR8('TEMP',VALUE)
        DIAL1%ACTUAL=VALUE      
        CALL DIALS_CHECKIN(NAME,DIAL1)
        NAME='SPIN'
        CALL DIALS_CHECKOUT(NAME,DIAL1)
        CALL DYNOCC$GETR8('SPIN',VALUE)
        DIAL1%ACTUAL=VALUE      
        CALL DIALS_CHECKIN(NAME,DIAL1)
        NAME='CHARGE'
        CALL DIALS_CHECKOUT(NAME,DIAL1)
        CALL DYNOCC$GETR8('TOTCHA',VALUE)
        DIAL1%ACTUAL=VALUE      
        CALL DIALS_CHECKIN(NAME,DIAL1)
        TSET=.TRUE.
        RETURN
      END IF
!
!     =================================================================
!     ==  TEMPERATURE OF THE ELECTRONS                               ==
!     =================================================================
      NAME='TEMP(E)'
      CALL DIALS_CHECKOUT(NAME,DIAL1)
      IF(DIAL1%ON) THEN
        CALL DIALS_TURN(DIAL1,DELTAT)
        VALUE=DIAL1%ACTUAL
        CALL DIALS_CHECKIN(NAME,DIAL1)
        CALL DYNOCC$SETR8('TEMP',VALUE)
      END IF
      NAME='SPIN'
      CALL DIALS_CHECKOUT(NAME,DIAL1)
      IF(DIAL1%ON) THEN
        CALL DIALS_TURN(DIAL1,DELTAT)
        VALUE=DIAL1%ACTUAL
        CALL DIALS_CHECKIN(NAME,DIAL1)
        CALL DYNOCC$SETR8('SPIN',VALUE)
      END IF
      NAME='CHARGE'
      CALL DIALS_CHECKOUT(NAME,DIAL1)
      IF(DIAL1%ON) THEN
        CALL DIALS_TURN(DIAL1,DELTAT)
        VALUE=DIAL1%ACTUAL
        CALL DIALS_CHECKIN(NAME,DIAL1)
        CALL DYNOCC$SETR8('TOTCHA',VALUE)
      END IF
!
!     =================================================================
!     ==  NEXT                                                       ==
!     =================================================================
      RETURN
      END
!
!     ...................................................AUTOPI.........
MODULE AUTOPILOT_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: AUTOPILOT                                                  **
!**                                                                   **
!**  PURPOSE: ANNEALING SCHEME FOR FRICTION DYNAMICS.                 **
!**     TUNES THE FRICTION                                            **
!**                                                                   **
!***********************************************************************
TYPE CONTROL_TYPE
 CHARACTER(32) :: NAME
 LOGICAL(1)    :: ON
 REAL(8)       :: LOWERFRICTION 
 REAL(8)       :: LOWERFACTOR
 REAL(8)       :: UPPERFRICTION 
 REAL(8)       :: UPPERFACTOR
 REAL(8)       :: MINFRICTION
END TYPE CONTROL_TYPE
INTEGER(4),PARAMETER :: NAUTOX=5
TYPE(CONTROL_TYPE)   :: CONTROL(NAUTOX) &
          & =(/CONTROL_TYPE('WAVES',.FALSE.,0.D0,1.D0,0.D0,1.D0,0.D0) &
          &   ,CONTROL_TYPE('ATOMS',.FALSE.,0.D0,1.D0,0.D0,1.D0,0.D0) &
          &   ,CONTROL_TYPE('OCC',.FALSE.,0.D0,1.D0,0.D0,1.D0,0.D0) &
          &   ,CONTROL_TYPE('QM-MM',.FALSE.,0.D0,1.D0,0.D0,1.D0,0.D0) &
          &   ,CONTROL_TYPE('COSMO',.FALSE.,0.D0,1.D0,0.D0,1.D0,0.D0)/)
INTEGER(4)           ::NAUTO=NAUTOX
LOGICAL(4) :: TINI=.FALSE.
LOGICAL(4) :: TAUTO=.FALSE.  ! SWITCH(AUTOPILOT ON/OFF)
LOGICAL(4) :: TFIRST=.TRUE.
REAL(8)    :: EPREV=1.D+20  ! ENERGY OF PREVIOUS STEP
INTEGER(4) :: NTOL=20       ! CONVERGENCE CRITERION: #(STEPS)
REAL(8)    :: ETOL=1.D-8    ! CONVERGENCE CRITERION: ENERGY WINDOW
REAL(8)    :: EMARK         !
INTEGER(4) :: ITER          ! #(ITERATIONS WITH |ETOT-EMARK|<ETOL
INTEGER(4) :: ISELECT       ! SELECTS ONE CONTROL FOR IO INTERFACE
END MODULE AUTOPILOT_MODULE
!
!     ...................................................AUTOPI.........
      SUBROUTINE AUTO$REPORT(NFIL)
!     ******************************************************************
!     **  REPORT SETTING OF AUTOPILOT                                 **
!     ******************************************************************
      USE AUTOPILOT_MODULE
      IMPLICIT NONE
      REAL(8)   ,PARAMETER  :: DSMALL=1.D-12
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: I
      INTEGER(4)            :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      IF(.NOT.TAUTO) RETURN
      CALL REPORT$TITLE(NFIL,'AUTOMATIC MINIMIZER')
      DO I=1,NAUTO
        IF(.NOT.CONTROL(I)%ON) CYCLE
        IF(ABS(CONTROL(I)%LOWERFRICTION-CONTROL(I)%UPPERFRICTION).LE.DSMALL &
     &    .AND.ABS(CONTROL(I)%LOWERFACTOR-1.D0).GT.DSMALL &
     &    .AND.ABS(CONTROL(I)%UPPERFACTOR-1.D0).GT.DSMALL) THEN
          CALL REPORT$R8VAL(NFIL &
     &                        ,'CONSTANT FRICTION FOR '//TRIM(CONTROL(I)%NAME) &
     &                                            ,CONTROL(I)%LOWERFRICTION,' ')
        ELSE
          CALL REPORT$R8VAL(NFIL,'LOWER FRICTION FOR '//TRIM(CONTROL(I)%NAME) &
                                    ,CONTROL(I)%LOWERFRICTION,' ')
          IF(ABS(CONTROL(I)%LOWERFACTOR-1.D0).GT.DSMALL) THEN
            CALL REPORT$R8VAL(NFIL,'SCALING FACTOR OF LOWER FRICTION FOR ' &
     &                       //TRIM(CONTROL(I)%NAME),CONTROL(I)%LOWERFACTOR,' ')
          END IF
          CALL REPORT$R8VAL(NFIL,'UPPER FRICTION FOR '//TRIM(CONTROL(I)%NAME) &
                                    ,CONTROL(I)%UPPERFRICTION,' ')
          IF(ABS(CONTROL(I)%UPPERFACTOR-1.D0).GT.DSMALL) THEN
            CALL REPORT$R8VAL(NFIL &
     &         ,'SCALING FACTOR OF UPPER FRICTION FOR '//TRIM(CONTROL(I)%NAME) &
     &                                              ,CONTROL(I)%UPPERFACTOR,' ')
          END IF
          IF(CONTROL(I)%MINFRICTION.GT.DSMALL) THEN
            CALL REPORT$R8VAL(NFIL &
     &                   ,'LOWEST FRICTION VALUE FOR '//TRIM(CONTROL(I)%NAME) &
     &                                              ,CONTROL(I)%MINFRICTION,' ')
          END IF 
        END IF
      ENDDO
      CALL REPORT$R8VAL(NFIL,'ENERGY TOLERANCE FOR TERMINATION',ETOL,' ')
      CALL REPORT$I4VAL(NFIL,'#STEPS BEFORE TERMINATION IS CONSIDERED',NTOL,' ')
      CALL REPORT$STRING(NFIL,'      FRICTION=0: NEWTON DYNAMICS')
      CALL REPORT$STRING(NFIL,'      FRICTION=1: STEEPEST DESCENT')
      RETURN
      END
!
!     ...................................................AUTOPI.........
      SUBROUTINE AUTOPI(TSTOP)
!     **                                                              **
!     **  AUTOPILOT SETS THE FRICTION PARAMETERS ANNEE AND ANNER      **
!     **                                                              **
      USE AUTOPILOT_MODULE
      IMPLICIT NONE
      REAL(8)                :: ETOT 
      LOGICAL(4),INTENT(OUT) :: TSTOP
      REAL(8)                :: FRICTION
      INTEGER(4)             :: I
!     ******************************************************************
      TSTOP=.FALSE.
      IF(.NOT.TAUTO) RETURN
!
!     ==================================================================
!     == CHECK CONVERGENCE                                            ==
!     ==================================================================
      CALL ENERGYLIST$GET('TOTAL ENERGY',ETOT)
      IF(TFIRST) THEN
        TFIRST=.FALSE.
        EMARK=ETOT
        ITER=1
      ELSE
        TSTOP=.FALSE.
        IF(ABS(ETOT-EMARK).GT.ETOL) THEN
          ITER=1
          EMARK=ETOT
        ELSE 
          IF(ITER.GT.NTOL) THEN
            TSTOP=.TRUE.
          END IF
          ITER=ITER+1
        ENDIF       
      END IF
!
!     ==================================================================
!     == APPLY AND SCALE FRICTIONS                                    ==
!     ==================================================================
      DO I=1,NAUTO
        IF(.NOT.CONTROL(I)%ON) CYCLE
        IF(ETOT.GT.EPREV) THEN
          FRICTION=CONTROL(I)%UPPERFRICTION
        ELSE
          FRICTION=CONTROL(I)%LOWERFRICTION
          CONTROL(I)%LOWERFRICTION=CONTROL(I)%LOWERFRICTION*CONTROL(I)%LOWERFACTOR
          CONTROL(I)%UPPERFRICTION=CONTROL(I)%UPPERFRICTION*CONTROL(I)%UPPERFACTOR
        END IF
        FRICTION=MAX(FRICTION,CONTROL(I)%MINFRICTION)
!
        IF(TRIM(CONTROL(I)%NAME).EQ.'ATOMS') THEN
           CALL ATOMS$SETR8('FRICTION',FRICTION)
        ELSE IF(TRIM(CONTROL(I)%NAME).EQ.'OCC') THEN
          CALL DYNOCC$SETR8('FRICTION',FRICTION)
        ELSE IF(TRIM(CONTROL(I)%NAME).EQ.'WAVES') THEN
          CALL WAVES$SETR8('FRICTION',FRICTION)
        ELSE IF(TRIM(CONTROL(I)%NAME).EQ.'QM-MM') THEN
          CALL QMMM$SETR8('FRICTION',FRICTION)
        ELSE IF(CONTROL(I)%NAME.EQ.'COSMO') THEN
          CALL COSMO$SETR8('FRICTION',FRICTION)
        END IF
      ENDDO
      EPREV=ETOT
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUTO$SELECT(ID_)
!     ******************************************************************
!     ******************************************************************
      USE AUTOPILOT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)              :: I
!     ******************************************************************
      IF(ID_.EQ.'~') THEN
        ISELECT=0
        RETURN
      END IF
      DO I=1,NAUTO
        IF(ID_.EQ.TRIM(CONTROL(I)%NAME)) THEN
          ISELECT=I
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('KEYWORD NOT RECOGNIZED')
      CALL ERROR$CHVAL('ID_',ID_)
      CALL ERROR$STOP('AUTO$SELECT')
      END
!
!     ..................................................................
      SUBROUTINE AUTO$SETL4(IDENT_,VALUE)
!     ******************************************************************
!     ******************************************************************
      USE AUTOPILOT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: IDENT_
      LOGICAL(4)  ,INTENT(IN) :: VALUE
      INTEGER(4)              :: I
!     ******************************************************************
      IF(ISELECT.EQ.0) THEN
        CALL ERROR$MSG('CONTROL MUST BE SELECTED')
        CALL ERROR$STOP('AUTO$SETL4')
      END IF
      IF(IDENT_.EQ.'ON') THEN
        CONTROL(ISELECT)%ON=VALUE
        TAUTO=.FALSE.
        DO I=1,NAUTO
          TAUTO=TAUTO.OR.CONTROL(I)%ON
        ENDDO
      ELSE
        CALL ERROR$MSG('KEYWORD NOT RECOGNIZED')
        CALL ERROR$CHVAL('IDENT_',IDENT_)
        CALL ERROR$STOP('AUTO$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUTO$SETI4(IDENT_,VALUE)
!     ******************************************************************
!     ******************************************************************
      USE AUTOPILOT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: IDENT_
      INTEGER(4)  ,INTENT(IN) :: VALUE
!     ******************************************************************
!     == GLOBAL VARIABLES (NOT SPECIFIC FOR A PARTICULAR AUTOPILOT)
      IF(IDENT_.EQ.'NTOL') THEN
        NTOL=VALUE
        RETURN
      END IF
!     == CHECK IF A PARTICULAR AUTOPILOT IS SELECTED
      IF(ISELECT.EQ.0) THEN
        CALL ERROR$MSG('CONTROL MUST BE SELECTED')
        CALL ERROR$STOP('AUTO$SETI4')
      END IF
!     == SETTINGS OF A PARTICULAR AUTOPILOT 
      IF(IDENT_.EQ.'XXX') THEN
!        XXX=VALUE
      ELSE
        CALL ERROR$MSG('KEYWORD NOT RECOGNIZED')
        CALL ERROR$CHVAL('IDENT_',IDENT_)
        CALL ERROR$STOP('AUTO$SETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUTO$SETR8(IDENT_,VALUE)
!     ******************************************************************
!     ******************************************************************
      USE AUTOPILOT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: IDENT_
      REAL(8)     ,INTENT(IN) :: VALUE
!     ******************************************************************
      IF(IDENT_.EQ.'TOLERANCE') THEN
        ETOL=VALUE
        RETURN
      END IF
      IF(ISELECT.EQ.0) THEN
        CALL ERROR$MSG('CONTROL MUST BE SELECTED')
        CALL ERROR$STOP('AUTO%SETR8')
      END IF
      IF(IDENT_.EQ.'LOWERFRICTION') THEN
        CONTROL(ISELECT)%LOWERFRICTION=VALUE
      ELSE IF(IDENT_.EQ.'UPPERFRICTION') THEN
        CONTROL(ISELECT)%UPPERFRICTION=VALUE
      ELSE IF(IDENT_.EQ.'UPPERFACTOR') THEN
        CONTROL(ISELECT)%UPPERFACTOR=VALUE
      ELSE IF(IDENT_.EQ.'LOWERFACTOR') THEN
        CONTROL(ISELECT)%LOWERFACTOR=VALUE
      ELSE IF(IDENT_.EQ.'MINFRICTION') THEN
        CONTROL(ISELECT)%MINFRICTION=VALUE
      ELSE
        CALL ERROR$MSG('KEYWORD NOT RECOGNIZED')
        CALL ERROR$CHVAL('IDENT_',IDENT_)
        CALL ERROR$STOP('AUTO$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUTO$GETL4(IDENT_,VALUE)
!     ******************************************************************
!     ******************************************************************
      USE AUTOPILOT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: IDENT_
      LOGICAL(4)  ,INTENT(OUT):: VALUE
!     ******************************************************************
      IF(ISELECT.EQ.0) THEN
        CALL ERROR$MSG('CONTROL MUST BE SELECTED')
        CALL ERROR$STOP('AUTO%GETL4')
      END IF
      IF(IDENT_.EQ.'ON') THEN
        VALUE=CONTROL(ISELECT)%ON
      ELSE
        CALL ERROR$MSG('KEYWORD NOT RECOGNIZED')
        CALL ERROR$CHVAL('IDENT_',IDENT_)
        CALL ERROR$STOP('AUTO$GETL8')
      END IF
      RETURN
      END 
