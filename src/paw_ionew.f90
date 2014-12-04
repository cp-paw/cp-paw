!======================================================================
!======================================================================
!=====                                                            =====
!=====   OBJECT FOR WRITING AND READING THE RESTART FILE          =====
!=====   (ACTUALLY NOT A WELL DEFINED OBJECT,YET)                 =====
!=====                                                            =====
!======================================================================
!======================================================================
!======================================================================
!
!     ..................................................................
      SUBROUTINE READRESTART
!     ******************************************************************
!     ******************************************************************
      USE RESTART_INTERFACE
      USE MPE_MODULE
      IMPLICIT NONE
      LOGICAL(4)            :: TCHK
      INTEGER(4)            :: NFILO
      INTEGER(4)            :: NFIL
      TYPE(SEPARATOR_TYPE)  :: SEPARATOR
      LOGICAL(4)            :: TREAD
      INTEGER(4)            :: NSTEP
      REAL(8)               :: DELTAT
      INTEGER(4)            :: NTASKS,THISTASK
      INTEGER(4)            :: NFILSTDOUT=6
!     ******************************************************************
                            CALL TRACE$PUSH('READRESTART')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      NFILO=NFILSTDOUT
      WRITE(NFILO,FMT='("READING RESTART FILE")')
!
!     =================================================================
!     == CHECK FILE FORMAT                                           ==
!     =================================================================
      IF(THISTASK.EQ.1) THEN      
        CALL FILEHANDLER$UNIT('RESTART_IN',NFIL)
        TCHK=.FALSE.
        REWIND NFIL
        READ(NFIL,END=1111)SEPARATOR
        TCHK=.TRUE.
 1111   CONTINUE
        IF((.NOT.TCHK).OR.(SEPARATOR%ID.NE.HEADER%ID)) THEN
          CALL ERROR$MSG('FORMAT OF RESTART FILE NOT RECOGNIZED') 
          CALL ERROR$STOP('READRESTART')
        END IF
        CALL RESTART$CHECK(NFIL)
        REWIND NFIL
      END IF
!
!     =================================================================
!     == READ FILE                                                   ==
!     =================================================================
      SEPARATOR=HEADER
      TCHK=.TRUE.
      IF(THISTASK.EQ.1)CALL RESTART$READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      CALL MPE$BROADCAST('MONOMER',1,TCHK)
      IF(TCHK) THEN
        IF(SEPARATOR%VERSION.NE.HEADER%VERSION) THEN
          CALL ERROR$STOP('READRESTART')
        END IF
      ELSE
        CALL TRACE$POP()
        RETURN
      END IF
 100  CONTINUE
      TREAD=.FALSE.
!
!     ==================================================================
!     ==  END OF FILE                                                 ==
!     ==================================================================
      TCHK=.TRUE.
      SEPARATOR=ENDOFFILE
      IF(THISTASK.EQ.1)CALL RESTART$READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK) 
      CALL MPE$BROADCAST('MONOMER',1,TCHK)
      IF(TCHK) THEN
        IF(THISTASK.EQ.1) THEN
          CALL FILEHANDLER$CLOSE('RESTART_IN')
        END IF
        CALL MPE$BROADCAST('MONOMER',1,TCHK)
        CALL TRACE$POP()
        RETURN
      END IF
!
!     ==================================================================
!     ==  READ TIMESTEP                                               ==
!     ==================================================================
      CALL TIMESTEP$READ(NFIL,NFILO,TCHK)
      CALL TIMESTEP$GETR8('DELTAT',DELTAT)
      CALL TIMESTEP$GETI4('ISTEP',NSTEP)
      TREAD=TREAD.OR.TCHK
!
!     ==================================================================
!     ==  READ THERMOSTAT FOR THE ATOMS                               ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$READ(NFIL,NFILO,TCHK)
      TREAD=TREAD.OR.TCHK
!
!     ==================================================================
!     ==  READ THERMOSTAT FOR THE WAVE FUNCTIONS                      ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$READ(NFIL,NFILO,TCHK)
      TREAD=TREAD.OR.TCHK
!
!     ==================================================================
!     ==  READ ATOMIC POSITIONS                                       ==
!     ==================================================================
      CALL ATOMS$READ(NFIL,NFILO,TCHK)
      TREAD=TREAD.OR.TCHK
!
!     ==================================================================
!     ==  READ LATTICE VECTORS                                        ==
!     ==================================================================
      CALL CELL$READ(NFIL,NFILO,TCHK)
      TREAD=TREAD.OR.TCHK
!
!     ==================================================================
!     ==  READ OCCUPATIONS                                            ==
!     ==================================================================
      CALL DYNOCC$READ(NFIL,NFILO,TCHK)
      TREAD=TREAD.OR.TCHK
!
!     ==================================================================
!     ==  READ WAVE FUNCTIONS                                         ==
!     ==================================================================
      CALL WAVES$READ(NFIL,NFILO,TCHK)
      TREAD=TREAD.OR.TCHK
!
!     ==================================================================
!     ==  READ MOLECULAR MECHANICS ENVIRONMENT                        ==
!     ==================================================================
      CALL QMMM$READ(NFIL,NFILO,TCHK)
      TREAD=TREAD.OR.TCHK
!
!     ==================================================================
!     ==  READ COSMO ENVIRONMENT                                  ==
!     ==================================================================
      CALL COSMO$READ(NFIL,NFILO,TCHK)
      TREAD=TREAD.OR.TCHK
!
!     ==================================================================
!     ==  UNIDENTIFIED OPTION                                         ==
!     ==================================================================
      IF(.NOT.TREAD) THEN
        IF(THISTASK.EQ.1)CALL RESTART$SKIP(NFIL,NFILO)
      END IF
      GOTO 100
      END
!
!     ..................................................................
      SUBROUTINE WRITERESTART
!     ******************************************************************
!     ******************************************************************
      USE RESTART_INTERFACE
      USE MPE_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      LOGICAL(4)          :: TTEST=.FALSE.
      INTEGER(4)          :: NFIL
      INTEGER(4)          :: NFILO
      REAL(8)             :: DELTAT=1.D0
      LOGICAL(4)          :: TCHK
      INTEGER(4)          :: NTASKS,THISTASK
      INTEGER(4),PARAMETER:: NFILSTDOUT=6
      CHARACTER(128)      :: FILENAME
!     ******************************************************************
                          CALL TRACE$PUSH('WRITERESTART')
!     CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==================================================================
!     == RETURN IF WRITING TO A SCRATCH FILE                          ==
!     ==================================================================
      CALL FILEHANDLER$FILENAME('RESTART_OUT',FILENAME)
      TCHK=(FILENAME.EQ.-'/DEV/NULL')
      CALL MPE$BROADCAST('MONOMER',1,TCHK)
      IF(TCHK) THEN
        CALL TRACE$POP()
        RETURN
      END IF
!
!     ==================================================================
!     == GET FILE UNIT                                                ==
!     ==================================================================
      NFILO=NFILSTDOUT
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.EQ.1) THEN
        CALL FILEHANDLER$UNIT('RESTART_OUT',NFIL)
        REWIND NFIL
      ELSE
        NFIL=-1
      END IF
!
!     ==================================================================
!     == WRITE HEADER                                                 ==
!     ==================================================================
      IF(THISTASK.EQ.1)CALL RESTART$WRITESEPARATOR(HEADER,NFIL,NFILO,TCHK)
!
!     ==================================================================
!     == WRITE THERMOSTAT FOR THE ATOMS                               ==
!     ==================================================================
      CALL TIMESTEP$WRITE(NFIL,NFILO,TCHK)
!
!     ==================================================================
!     == WRITE THERMOSTAT FOR THE ATOMS                               ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$WRITE(NFIL,NFILO,TCHK)
!
!     ==================================================================
!     == WRITE THERMOSTAT FOR THE WAVE FUNCTIONS                      ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$WRITE(NFIL,NFILO,TCHK)
!
!     ==================================================================
!     == WRITE ATOMIC POSITIONS                                       ==
!     ==================================================================
      CALL ATOMS$WRITE(NFIL,NFILO,TCHK)
!
!     ==================================================================
!     == WRITE ATOMIC POSITIONS                                       ==
!     ==================================================================
      CALL CELL$WRITE(NFIL,NFILO,TCHK)
!
!     ==================================================================
!     == WRITE OCCUPATIONS                                            ==
!     ==================================================================
      CALL DYNOCC$WRITE(NFIL,NFILO,TCHK)
!
!     ==================================================================
!     == WRITE WAVE FUNCTIONS                                         ==
!     ==================================================================
      CALL WAVES$WRITE(NFIL,NFILO,TCHK)
!
!     ==================================================================
!     == WRITE MOLECULAR MECHANICS ENVIRONMENT                        ==
!     ==================================================================
      CALL QMMM$WRITE(NFIL,NFILO,TCHK)
!
!     ==================================================================
!     == WRITE COSMO ENVIRONMENT                                      ==
!     ==================================================================
      CALL COSMO$WRITE(NFIL,NFILO,TCHK)
!
!     ==================================================================
!     == END OF FILE                                                  ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        TCHK=.TRUE.
        CALL RESTART$WRITESEPARATOR(ENDOFFILE,NFIL,NFILO,TCHK)
      END IF
!
!     ==================================================================
!     == CHECK CONSISTENCY OF THE FILE                                ==
!     ==================================================================
      IF(TTEST.AND.THISTASK.EQ.1) THEN
        PRINT*,'CHECK CONSISTENCY OF RESTART FILE IN WRITERESTART'
        CALL RESTART$CHECK(NFIL)
      END IF
      CALL FILEHANDLER$CLOSE('RESTART_OUT')
      CALL MPE$BROADCAST('MONOMER',1,TCHK)
                                                     CALL TRACE$POP()
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP$READ(NFIL,NFILO,TCHK)
!     ******************************************************************
!     ******************************************************************
      USE RESTART_INTERFACE
      USE TIMESTEP_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN)  :: NFIL
      INTEGER(4)            ,INTENT(IN)  :: NFILO
      LOGICAL(4)            ,INTENT(OUT) :: TCHK
      TYPE (SEPARATOR_TYPE) ,PARAMETER   :: MYSEPARATOR &
             =SEPARATOR_TYPE(1,'TIMESTEP','NONE','AUG1996','NONE')
      TYPE (SEPARATOR_TYPE)              :: SEPARATOR
      INTEGER(4)                         :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      TCHK=.TRUE.
      SEPARATOR=MYSEPARATOR
      IF(THISTASK.EQ.1)CALL RESTART$READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      CALL MPE$BROADCAST('MONOMER',1,TCHK)
      IF(.NOT.TCHK) RETURN
!
      IF(SEPARATOR%VERSION.NE.MYSEPARATOR%VERSION) THEN
        CALL ERROR$MSG('VERSION NOT RECOGNIZED')
        CALL ERROR$CHVAL('ACTUALVERSION ',MYSEPARATOR%VERSION)
        CALL ERROR$CHVAL('VERSION ',SEPARATOR%VERSION)
        CALL ERROR$STOP('READTIMESTEP')
      END IF
      IF(THISTASK.EQ.1) THEN
        READ(NFIL)DELTAT,ISTEPNUMBER
      END IF
      CALL MPE$BROADCAST('MONOMER',1,ISTEPNUMBER)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP$WRITE(NFIL,NFILO,TCHK)
!     ******************************************************************
      USE RESTART_INTERFACE
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      INTEGER(4)       ,INTENT(IN)  :: NFIL
      INTEGER(4)       ,INTENT(IN)  :: NFILO
      LOGICAL(4)       ,INTENT(OUT) :: TCHK
      TYPE (SEPARATOR_TYPE)         :: MYSEPARATOR &
           =SEPARATOR_TYPE(1,'TIMESTEP','NONE','AUG1996','NONE')
      INTEGER(4)                    :: NTASKS,THISTASK
!     ******************************************************************
      TCHK=.TRUE.
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.EQ.1) THEN
        CALL RESTART$WRITESEPARATOR(MYSEPARATOR,NFIL,NFILO,TCHK)
        WRITE(NFIL)DELTAT,ISTEPNUMBER
      END IF
      RETURN
      END
