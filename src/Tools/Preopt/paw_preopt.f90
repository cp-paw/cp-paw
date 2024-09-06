!     ******************************************************************
!     **     CP-PAW PREOPTIMIZATION TOOL                              **
!     ******************************************************************
!     **                                                              **
!     **     THIS TOOL IS TO BE USED WITH THE CP-PAW CODE             **
!     **     AND SHOULD PRE-OPTIMIZE STRUCTURE INPUT FILES            **
!     **     USING THE CLASSICAL MOLECULAR DYNAMICS IMPLEMENTED       **
!     **     INTO THE CP-PAW CODE.                                    **
!     **                                                              **
!     **  AUTHOR: JOHANNES SCHIMPL  2003                              **
!     **                                                              **
!     **  COMMENTS:                                                   **
!     **   ALWAYS USES REAL MASS                                      **
!     **   BOND READING NOT YET IMPLEMENTED                           **
!     **   IT SHOULD BE CHECKED IF MD CONVERGES BETTER THAN CG LINE SEACH
!     **                                                              **
!     **   DOCUMENTATION (MANUAL)                                     **
!     **                                                              **
!     ******************************************************************
!
!     ..................................................................
MODULE POPT_MODULE
      TYPE ATOM_TYPE
        CHARACTER(32)     :: NAME
        REAL(8)           :: R(3)
        REAL(8)           :: Q
        REAL(8)           :: M
        CHARACTER(5)      :: FFTYPE
        LOGICAL(4)        :: FREEZE
        INTEGER(4)        :: QMSATOM   ! M-ATOM FOR Q ATOM
                                       ! S ATOM FOR M-ATOM
                                       ! Q ATOM FOR S ATOM
      END TYPE ATOM_TYPE
      TYPE LINK_TYPE
        INTEGER(4)        :: MJOINT
        INTEGER(4)        :: QJOINT
        INTEGER(4)        :: SJOINT
        INTEGER(4)        :: MATOM
        INTEGER(4)        :: QATOM
        INTEGER(4)        :: SATOM
      END TYPE LINK_TYPE
      TYPE BOND_TYPE
        INTEGER(4)        :: ATOM1
        INTEGER(4)        :: ATOM2
        REAL(8)           :: BO
      END TYPE BOND_TYPE
      INTEGER(4)                  :: NATQ,NATM,NATS
      INTEGER(4)                  :: NBONDM,NBONDS
      INTEGER(4)                  :: NLINK
      INTEGER(4)                  :: NMAP
      REAL(8)                     :: LUNIT
      TYPE(ATOM_TYPE),ALLOCATABLE :: MATOM(:)
      TYPE(BOND_TYPE),ALLOCATABLE :: MBOND(:)
      LOGICAL(4)                  :: TPRINTATOMS=.TRUE.
      LOGICAL(4)                  :: TPRINTBONDS=.TRUE.
      LOGICAL(4)                  :: TWARNFF=.TRUE.
      LOGICAL(4)                  :: TXYZOUT=.TRUE.
      REAL(8)                     :: TOL=1.D-4
      INTEGER(4)                  :: NITER=1000
END MODULE POPT_MODULE

!     ******************************************************************
      PROGRAM PREOPT
!     ******************************************************************
      USE POPT_MODULE
      IMPLICIT NONE
!     ******************************************************************
      CALL TRACE$SETL4('ON',.FALSE.)
      CALL TRACE$PUSH('MAIN')

!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT
!
!     ==================================================================
!     ==  RESOLVE ARGUMENTLIST AND INITIALIZE FILE HANDLER            ==
!     ==================================================================
      CALL INITIALIZEFILEHANDLER

      CALL READCNTL

      CALL READSTRC
      CALL CLASSICAL$SETCH('FF','UFF')

      CALL REPORT_START

      CALL MM_MINIMIZE_DYNAMIC
!      CALL MM_MINIMIZE_CG

      CALL REPORT_END

      DEALLOCATE(MATOM)
      DEALLOCATE(MBOND)
      CALL FILEHANDLER$CLOSEALL
      CALL ERROR$NORMALSTOP()
      END PROGRAM PREOPT
!     ==================================================================
!     ==================================================================
!     ==================================================================

!
!     ..................................................................
MODULE READ_MODULE
USE LINKEDLIST_MODULE
TYPE(LL_TYPE)   :: LL_CNTL
TYPE(LL_TYPE)   :: LL_STRC
SAVE
END MODULE READ_MODULE

!
!     ..................................................................
      SUBROUTINE READSTRC
      USE LINKEDLIST_MODULE
      USE READ_MODULE
      USE POPT_MODULE
      USE STRINGS_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE 
      INTEGER(4)                  :: NFIL
      REAL(8)                     :: PROTONMASS
      INTEGER(4)                  :: IATQ,IATM,IATS,IATM1,IATM2,IATS1,IATS2
      INTEGER(4)                  :: IBONDM,IBONDS,ILINK,I,JATM,MAXBOND
      INTEGER(4)                  :: IZ
      LOGICAL(4)                  :: TRC    ! ATOM IS PART OF THE RECTION CENTER OR DUMMY ATOM
      LOGICAL(4)                  :: TCHK
      INTEGER(4)     ,ALLOCATABLE :: LINKARRAY(:,:)
      INTEGER(4)     ,ALLOCATABLE :: MAPARRAY(:,:)
      REAL(8)        ,ALLOCATABLE :: RCOV(:) 
      CHARACTER(32)               :: CH32SVAR1
!     ******************************************************************    
!    
!     ==================================================================
!     ==  OPEN STRUCTURE FILE                                         ==
!     ==================================================================
      CALL LINKEDLIST$NEW(LL_STRC)
      CALL FILEHANDLER$UNIT('STRC',NFIL)
      CALL LINKEDLIST$READ(LL_STRC,NFIL,'~')

      CALL CONSTANTS('U',PROTONMASS)
!    
!     ==================================================================
!     ==  READ UNIT FROM BLOCK !STRUCTURE!GENERIC                     ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'GENERIC')
      CALL LINKEDLIST$GET(LL_STRC,'LUNIT',1,LUNIT)
!
!     ==================================================================
!     ==  READ BLOCK ATOMS                                            ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')

      CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NATM)
      ALLOCATE(MATOM(NATM))
      ALLOCATE(RCOV(NATM))
      DO IATM=1,NATM
         CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IATM)
!
!       ==  ATOM NAME  =================================================
         CALL LINKEDLIST$EXISTD(LL_STRC,'NAME',1,TCHK)
         IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('VARIABLE !ATOM:NAME NOT FOUND')
            CALL ERROR$STOP('READSTRC') 
         END IF
         CALL LINKEDLIST$GET(LL_STRC,'NAME',1,MATOM(IATM)%NAME)
         CALL LINKEDLIST$EXISTD(LL_STRC,'SP',1,TCHK)
         IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_STRC,'SP',1,CH32SVAR1)
            IF(+CH32SVAR1(1:2).NE.+MATOM(IATM)%NAME(1:2)) &
                 MATOM(IATM)%NAME=+CH32SVAR1(1:2)//MATOM(IATM)%NAME
         END IF
         MATOM(IATM)%QMSATOM=0
!
!       == POSITION ====================================================
         CALL LINKEDLIST$EXISTD(LL_STRC,'R',1,TCHK)
         IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('KEYWORD !ATOM:R NOT FOUND')
            CALL ERROR$STOP('READSTRC') 
         END IF
         CALL LINKEDLIST$GET(LL_STRC,'R',1,MATOM(IATM)%R)
         MATOM(IATM)%R(:)=MATOM(IATM)%R(:)*LUNIT
!
!       ==  FORCE FIELD ATOM TYPE========================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'FFTYPE',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL LINKEDLIST$SET(LL_STRC,'FFTYPE',0,'')
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'FFTYPE',1,MATOM(IATM)%FFTYPE)
!
!      ==  MASS ========================================================
         CH32SVAR1=MATOM(IATM)%NAME(1:2)
         IF(CH32SVAR1(2:2).EQ.'_') CH32SVAR1(2:2)=' '
         CALL PERIODICTABLE$GET(CH32SVAR1,'Z',IZ)
         CALL PERIODICTABLE$GET(IZ,'MASS',MATOM(IATM)%M)
!
!       == KOVALENT RADIUS =============================================
         CALL PERIODICTABLE$GET(IZ,'R(COV)',RCOV(IATM))
         MATOM(IATM)%Q=0.D0
         MATOM(IATM)%QMSATOM=0

         CALL LINKEDLIST$SELECT(LL_STRC,'..')
      END DO

      CALL CLASSICAL$SELECT('QMMM')
      MAXBOND=3*NATM
      ALLOCATE(MBOND(MAXBOND))
      NBONDM=0
!
!     == BOND INFORMATION =============================================
      CALL FINDBONDS(MAXBOND,RCOV)
      DEALLOCATE(RCOV)
      CALL FIND_FF
!
!     == READ FREEZE ==================================================
      CALL READ_FREEZE

      CALL STRCIN_SETM(NATM,MATOM,NBONDM,MBOND(1:NBONDM))
      
    END SUBROUTINE READSTRC
!       ................................................................
        SUBROUTINE READ_FREEZE
!       ****************************************************************    
!       **                                                            **    
!       ****************************************************************    
        USE POPT_MODULE
        USE READ_MODULE
        USE STRINGS_MODULE
        IMPLICIT NONE
        INTEGER(4)            :: IAT,NLIST,ILIST,IVAR
        CHARACTER(32)         :: STR
        LOGICAL(4)            :: TCHK
!       ****************************************************************    
        MATOM(:)%FREEZE=.FALSE.
        CALL LINKEDLIST$SELECT(LL_CNTL,'~')
        CALL LINKEDLIST$SELECT(LL_CNTL,'PCNTL')
!
!     == FREEZE =======================================================
        CALL LINKEDLIST$EXISTL(LL_CNTL,'FREEZEONLY',1,TCHK)
        IF(TCHK) THEN
           CALL LINKEDLIST$SELECT(LL_CNTL,'FREEZEONLY')
           MATOM(:)%FREEZE=.FALSE.
           CALL LINKEDLIST$NLISTS(LL_CNTL,'ATOM',NLIST)
           DO ILIST=1,NLIST
              CALL LINKEDLIST$SELECT(LL_CNTL,'ATOM',ILIST)
              CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
              IF(TCHK) THEN
                 CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,STR)
                 IVAR=0
                 DO IAT=1,NATM
                    IF(+MATOM(IAT)%NAME.EQ.TRIM(+STR)) THEN
                       IVAR=IAT
                       EXIT
                    END IF
                 END DO
                 IF(IVAR.EQ.0) THEN
                    CALL ERROR$MSG('ATOM NOT IN THE ATOM LIST')
                    CALL ERROR$CHVAL('ATOM',STR)
                    CALL ERROR$STOP('READ_FREEZE !PCNTL!FREEZEONLY!ATOM')
                 END IF
                 MATOM(IVAR)%FREEZE=.TRUE.
              END IF
              CALL LINKEDLIST$EXISTD(LL_CNTL,'PART',1,TCHK)
              IF(TCHK) THEN
                 CALL LINKEDLIST$GET(LL_CNTL,'PART',1,STR)
                 DO IAT=1,NATM
                    IF(INDEX(+MATOM(IAT)%NAME,TRIM(+STR)).GT.0) &
                         MATOM(IAT)%FREEZE=.TRUE.
                 END DO
              END IF
              CALL LINKEDLIST$SELECT(LL_CNTL,'..')
           END DO
           CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        END IF
!
!     == MOVE =========================================================
        CALL LINKEDLIST$EXISTL(LL_CNTL,'MOVEONLY',1,TCHK)
        IF(TCHK) THEN
           CALL LINKEDLIST$SELECT(LL_CNTL,'MOVEONLY')
           MATOM(:)%FREEZE=.TRUE.
           CALL LINKEDLIST$NLISTS(LL_CNTL,'ATOM',NLIST)
           DO ILIST=1,NLIST
              CALL LINKEDLIST$SELECT(LL_CNTL,'ATOM',ILIST)
              CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
              IF(TCHK) THEN
                 CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,STR)
                 IVAR=0
                 DO IAT=1,NATM
                    IF(+MATOM(IAT)%NAME.EQ.TRIM(+STR)) THEN
                       IVAR=IAT
                       EXIT
                    END IF
                 END DO
                 IF(IVAR.EQ.0) THEN
                    CALL ERROR$MSG('ATOM NOT IN THE ATOM LIST')
                    CALL ERROR$CHVAL('ATOM',STR)
                    CALL ERROR$STOP('READ_FREEZE !PCNTL!MOVEONLY!ATOM')
                 END IF
                 MATOM(IVAR)%FREEZE=.FALSE.
              END IF
              CALL LINKEDLIST$EXISTD(LL_CNTL,'PART',1,TCHK)
              IF(TCHK) THEN
                 CALL LINKEDLIST$GET(LL_CNTL,'PART',1,STR)
                 DO IAT=1,NATM
                    IF(INDEX(+MATOM(IAT)%NAME,TRIM(+STR)).GT.0) &
                         MATOM(IAT)%FREEZE=.FALSE.
                 END DO
              END IF
              CALL LINKEDLIST$SELECT(LL_CNTL,'..')
           END DO
           CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        END IF
!
!     == AS LAST POINT: READ FREEZE FOR EACH ATOM FROM STRC ===========
        CALL LINKEDLIST$SELECT(LL_STRC,'~')
        CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')

        CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',IAT)       
        IF(IAT.NE.NATM) THEN
           CALL ERROR$MSG('NUMBER OF ATOMS NOT CONSITENT')
           CALL ERROR$STOP('READ_FREEZE')
        END IF
        DO IAT=1,NATM
           CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IAT)
           CALL LINKEDLIST$EXISTD(LL_STRC,'FREEZE',1,TCHK)
           IF(TCHK) THEN
              CALL LINKEDLIST$GET(LL_STRC,'FREEZE',1,MATOM(IAT)%FREEZE)
           END IF
           CALL LINKEDLIST$SELECT(LL_STRC,'..')
        END DO
      END SUBROUTINE READ_FREEZE

!       ................................................................
        SUBROUTINE FINDBONDS(MAXBOND,RCOV)
!       ****************************************************************    
!       **                                                            **    
!       ****************************************************************    
        USE POPT_MODULE
        IMPLICIT NONE
        INTEGER(4), INTENT(IN)      :: MAXBOND
        REAL(8)   , INTENT(IN)      :: RCOV(NATM)
        INTEGER(4)                  :: IAT,JAT
        REAL(8)                     :: SVAR
!       ****************************************************************    
        DO IAT=1,NATM
           DO JAT=IAT+1,NATM
              SVAR=SQRT(SUM((MATOM(IAT)%R-MATOM(JAT)%R)**2))
              IF(SVAR.LT.((RCOV(IAT)+RCOV(JAT))*1.3D0)) THEN !MAKE BOND
                 NBONDM=NBONDM+1
                 IF(NBONDM.GT.MAXBOND) THEN
                    CALL ERROR$MSG('PARAMETER MAXBOND TOO SMALL')
                    CALL ERROR$STOP('FINDBONDS')
                 END IF
                 MBOND(NBONDM)%ATOM1=IAT
                 MBOND(NBONDM)%ATOM2=JAT
                 MBOND(NBONDM)%BO=1.D0  ! ONLY ORDER 1 BONDS UP TO NOW
              END IF
           END DO
        END DO
      END SUBROUTINE FINDBONDS
!       ................................................................
        SUBROUTINE FIND_FF
!       ****************************************************************    
!       **                                                            **    
!       ****************************************************************    
        USE POPT_MODULE
        IMPLICIT NONE
        INTEGER(4)             :: IAT,IBOND,NN,NFILO
        LOGICAL(4)             :: EXIST
        CHARACTER(5)           :: TYPE,STR
!       ****************************************************************    
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        DO IAT=1,NATM
           IF(TRIM(MATOM(IAT)%FFTYPE).NE.'') CYCLE ! ONLY GUESS IF NOT SPECIFIED
           NN=0
           DO IBOND=1,NBONDM
              IF(MBOND(IBOND)%ATOM1.EQ.IAT) NN=NN+1
              IF(MBOND(IBOND)%ATOM2.EQ.IAT) NN=NN+1
           END DO
           IF(NN.EQ.1.OR.NN.EQ.0) THEN
              STR=MATOM(IAT)%NAME(1:2)
           ELSE
              WRITE(STR,"(A2,I1)") MATOM(IAT)%NAME(1:2),NN-1
           END IF
           CALL UFFTABLE$EXIST(TRIM(STR),EXIST,TYPE)
           IF(EXIST) THEN
              MATOM(IAT)%FFTYPE=TYPE
           ELSE
              CALL UFFTABLE$EXIST(MATOM(IAT)%NAME(1:2),EXIST,TYPE)
              IF(TRIM(TYPE).NE.'') THEN
                 MATOM(IAT)%FFTYPE=TYPE
                 IF(TWARNFF) WRITE(NFILO,"('WARNING: UFFTYPE ',A,' USED &
                      &FOR ATOM ',A,' (',I2,' BONDS)')") &
                      TYPE,TRIM(MATOM(IAT)%NAME),NN
              ELSE
                 CALL ERROR$MSG('NO FITTING FORCE FIELD TYPE FOUND, PLEASE &
                      &SPECIFY IT MANUALLY!')
                 CALL ERROR$CHVAL('ATOM',MATOM(IAT)%NAME)
                 CALL ERROR$STOP('FIND_FF')
              END IF
           END IF
        END DO
      END SUBROUTINE FIND_FF
!
!       ................................................................
        SUBROUTINE REPORT_START
!       ****************************************************************    
!       **                                                            **    
!       ****************************************************************    
        USE POPT_MODULE
        IMPLICIT NONE
        INTEGER(4)                  :: IAT,JAT,IBOND,NFILO,LENGTH
        REAL(8)                     :: SVAR
        CHARACTER(256)              :: LINE
        CHARACTER(32)               :: STR
!       ****************************************************************    
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,*)
        IF(TPRINTATOMS) THEN
           WRITE(NFILO,"(' ATOMS AS USED BY THE MINIMIZER ')")
           WRITE(NFILO,"('================================')")
           LENGTH=1
           DO IAT=1,NATM
              LENGTH=MAX(LENGTH,LEN(TRIM(MATOM(IAT)%NAME)))
           END DO
           DO IAT=1,NATM
              STR='                              '
              LINE="!ATOM NAME='"//TRIM(MATOM(IAT)%NAME)//"'"//&
                   &STR(1:LENGTH-LEN(TRIM(MATOM(IAT)%NAME)))//" R="
              IF(MATOM(IAT)%FREEZE) THEN
                 STR='T'
              ELSE
                 STR='F'
              END IF
              WRITE(LINE,"(A,3F10.5,' FREEZE=',A1,' FFTYPE=')") &
                   TRIM(LINE),MATOM(IAT)%R/LUNIT,TRIM(STR)
              LINE=TRIM(LINE)//"'"//TRIM(MATOM(IAT)%FFTYPE)//"' !END"
              WRITE(NFILO,*) TRIM(LINE)
           END DO
           WRITE(NFILO,*)
        END IF
        IF(TPRINTBONDS) THEN
           WRITE(NFILO,"(' BONDS AS USED BY THE MINIMIZER ')")
           WRITE(NFILO,"('================================')")
           DO IBOND=1,NBONDM
              IAT=MBOND(IBOND)%ATOM1
              JAT=MBOND(IBOND)%ATOM2
              LINE="!BOND ATOM1='"//TRIM(MATOM(IAT)%NAME)//&
                   &"' ATOM2='"//TRIM(MATOM(JAT)%NAME)//"'"
              WRITE(LINE,"(A,' BO=',F3.1,' !END')") TRIM(LINE),MBOND(IBOND)%BO
              WRITE(NFILO,*) TRIM(LINE)
           END DO
           WRITE(NFILO,*)
        END IF
      END SUBROUTINE REPORT_START
!       ................................................................
        SUBROUTINE REPORT_END
!       ****************************************************************    
!       **                                                            **    
!       ****************************************************************    
        USE POPT_MODULE
        IMPLICIT NONE
        INTEGER(4)                  :: IAT,NFILO,NFIL,LENGTH
        REAL(8)                     :: SVAR,R(3*NATM),ANGSTROM
        CHARACTER(256)              :: LINE
        CHARACTER(32)               :: STR
!       ****************************************************************    
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        CALL CLASSICAL$GETR8A('R(0)',3*NATM,R)
        LENGTH=1
        DO IAT=1,NATM
           MATOM(IAT)%R=R(3*IAT-2:3*IAT)
           LENGTH=MAX(LENGTH,LEN(TRIM(MATOM(IAT)%NAME)))
        END DO

        WRITE(NFILO,"(' ATOMS AFTER OPTIMIZATION ')")
        WRITE(NFILO,"('==========================')")
        STR='                              '
        DO IAT=1,NATM
          LINE="!ATOM NAME='"//TRIM(MATOM(IAT)%NAME)//"'"//&
               &STR(1:LENGTH-LEN(TRIM(MATOM(IAT)%NAME)))//" R="
          WRITE(LINE,"(A,3F12.7,' !END')") TRIM(LINE),MATOM(IAT)%R/LUNIT
          WRITE(NFILO,*) TRIM(LINE)
        END DO
        WRITE(NFILO,*)
        IF(TXYZOUT) THEN
           CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
           CALL FILEHANDLER$UNIT('XYZ',NFIL)
           WRITE(NFIL,*) NATM
           WRITE(NFIL,*) '   PREOPTIMIZED STRUCTURE'
           DO IAT=1,NATM
              STR=MATOM(IAT)%NAME
              IF(STR(2:2).EQ.'_') STR(2:2)=' '
              WRITE(NFIL,"(A2,2X,3(F10.5,1X),2X,A)") &
                   STR(1:2),MATOM(IAT)%R/ANGSTROM,TRIM(MATOM(IAT)%NAME)
           END DO
        END IF
      END SUBROUTINE REPORT_END
!
!       ................................................................
        SUBROUTINE STRCIN_SETM(NAT,ATOM,NBOND,BOND)  
!       ****************************************************************    
!       **                                                            **    
!       ****************************************************************    
        USE POPT_MODULE
        IMPLICIT NONE
        INTEGER(4)     ,INTENT(IN) :: NAT
        INTEGER(4)     ,INTENT(IN) :: NBOND
        TYPE(ATOM_TYPE),INTENT(IN) :: ATOM(NAT)
        TYPE(BOND_TYPE),INTENT(IN) :: BOND(NBOND)
        INTEGER(4)                 :: IAT,IBOND
        REAL(8)                    :: R(3,NAT)
        LOGICAL(4)                 :: TFREEZE(NAT)
        REAL(8)                    :: MASS(NAT)
        REAL(8)                    :: CHARGE(NAT)
        CHARACTER(5)               :: FFTYPE(NAT)
        INTEGER(4)                 :: INDEX2(2,NBOND)
        REAL(8)                    :: BO(NBOND)
!       ****************************************************************    
!
!       ================================================================
!       ==  SEND ATOM AND BOND INFORMATION TO CLASSICAL OBJECT        ==
!       ================================================================
        DO IAT=1,NAT
          R(:,IAT)   =ATOM(IAT)%R(:)
          TFREEZE(IAT)=ATOM(IAT)%FREEZE
          MASS(IAT)  =ATOM(IAT)%M
          CHARGE(IAT)=ATOM(IAT)%Q
          FFTYPE(IAT)=ATOM(IAT)%FFTYPE
        ENDDO
        CALL CLASSICAL$SETI4('NAT',NAT)
        CALL CLASSICAL$SETR8A('R(0)',3*NAT,R)
        CALL CLASSICAL$SETR8A('R(-)',3*NAT,R)
        CALL CLASSICAL$SETR8A('MASS',NAT,MASS)
        CALL CLASSICAL$SETR8A('QEL',NAT,CHARGE)
        CALL CLASSICAL$SETCHA('TYPE',NAT,FFTYPE)
        CALL CLASSICAL$SETL4A('TFREEZE',NAT,TFREEZE)
        CALL CLASSICAL$SETCHA('ATOMNAME',NAT,ATOM(1:NAT)%NAME)
!
!       ================================================================
!       ==  SEND BOND INFORMATION TO CLASSICAL OBJECT                 ==
!       ================================================================
        DO IBOND=1,NBOND
          INDEX2(1,IBOND)=BOND(IBOND)%ATOM1
          INDEX2(2,IBOND)=BOND(IBOND)%ATOM2
          BO(IBOND)=BOND(IBOND)%BO
        ENDDO
        CALL CLASSICAL$SETI4('NBOND',NBOND)
        CALL CLASSICAL$SETI4A('INDEX2',2*NBOND,INDEX2)
        CALL CLASSICAL$SETR8A('BONDORDER',NBOND,BO)
        CALL CLASSICAL$SETL4('LONGRANGE',.TRUE.)
      END SUBROUTINE STRCIN_SETM
!
!     ..................................................................
      SUBROUTINE READCNTL
      USE LINKEDLIST_MODULE
      USE READ_MODULE
      USE POPT_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER :: TPR=.FALSE.
      LOGICAL(4)           :: TCHK,TCHK2
      INTEGER(4)           :: NFIL
      CHARACTER(32)        :: ID 
      CHARACTER(256)       :: FILENAME
      INTEGER(4)           :: ITH
      INTEGER(4)           :: NUM
      INTEGER(4)           :: NFILO
!     ****************************************************************
                          CALL TRACE$PUSH('READCNTL')
!
!     ==================================================================
!     ==  READ CONTROL FILE                                           ==
!     ==================================================================
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL FILEHANDLER$UNIT('PCNTL',NFIL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'~')
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO) 
        CALL LINKEDLIST$REPORT(LL_CNTL,NFILO)
      END IF
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'PCNTL')

      CALL LINKEDLIST$EXISTL(LL_CNTL,'GENERIC',1,TCHK)
      IF(TCHK) THEN
         CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')

         CALL LINKEDLIST$EXISTD(LL_CNTL,'TRACE',1,TCHK)
         IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_CNTL,'TRACE',1,TCHK2)
            CALL TRACE$SETL4('ON',TCHK2)
         END IF

         CALL LINKEDLIST$EXISTD(LL_CNTL,'TOL',1,TCHK)
         IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'TOL',1,TOL)

         CALL LINKEDLIST$EXISTD(LL_CNTL,'NSTEPS',1,TCHK)
         IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'NSTEPS',1,NITER)
         CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      END IF

      CALL LINKEDLIST$EXISTL(LL_CNTL,'FILES',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILES')
        CALL LINKEDLIST$NLISTS(LL_CNTL,'FILE',NUM)
        DO ITH=1,NUM
          CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',ITH)
          CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
          IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.FALSE.)
  !       ==  READ ACTUAL VALUES  ======================================
          CALL LINKEDLIST$GET(LL_CNTL,'ID',1,ID)
          CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILENAME)
          CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
          CALL FILEHANDLER$SETFILE(+ID,TCHK,FILENAME)
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'PCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'OUTPUT',1,TCHK)
      IF(TCHK) THEN
         CALL LINKEDLIST$SELECT(LL_CNTL,'OUTPUT')

         CALL LINKEDLIST$EXISTD(LL_CNTL,'PRINTATOMS',1,TCHK)
         IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'PRINTATOMS',1,TPRINTATOMS)

         CALL LINKEDLIST$EXISTD(LL_CNTL,'PRINTBONDS',1,TCHK)
         IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'PRINTBONDS',1,TPRINTBONDS)

         CALL LINKEDLIST$EXISTD(LL_CNTL,'WARNFF',1,TCHK)
         IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'WARNFF',1,TWARNFF)

         CALL LINKEDLIST$EXISTD(LL_CNTL,'XYZOUT',1,TCHK)
         IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'XYZOUT',1,TXYZOUT)
         CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      END IF
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE READCNTL
!      
!     ...................................................................
      SUBROUTINE INITIALIZEFILEHANDLER
      USE STRINGS_MODULE
      CHARACTER(256) :: ROOTNAME
      CHARACTER(256) :: INNAME
      INTEGER(4)     :: ISVAR,NFILO
      INTEGER        :: NARG
!     ******************************************************************
      NARG=COMMAND_ARGUMENT_COUNT()
      IF(NARG.LT.1) THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$MSG('THE CONTROL FILE OF THE PREOPTIMIZATION TOOL IS MANDATORY')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      END IF
      CALL GET_COMMAND_ARGUMENT(1,INNAME)
      ISVAR=INDEX(INNAME,-'.PCNTL',BACK=.TRUE.)
      IF(ISVAR.NE.0) THEN
        ROOTNAME=INNAME(1:ISVAR-1)
      ELSE
        ROOTNAME=' '
      END IF
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL STANDARDFILES
      CALL FILEHANDLER$SETFILE('PCNTL',.FALSE.,INNAME)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==================================================================
!     ==  WRITE HEADER                                                ==
!     ==================================================================
      WRITE(NFILO,FMT='(80("*"))')
      WRITE(NFILO,FMT='(80("*"),T15 &
     &             ,"           PREOPTIMIZATION TOOL              ")')
      WRITE(NFILO,FMT='(80("*"),T15 &
     &             ,"    FOR THE PROJECTOR-AUGMENTED WAVE METHOD  ")')
      WRITE(NFILO,FMT='(80("*"))')
      WRITE(NFILO,FMT='(T10 &
     &           ,"J. SCHIMPL, CLAUSTHAL UNIVERSITY OF TECHNOLOGY")')
      WRITE(NFILO,FMT='(T10 &
     &            ,"DISTRIBUTED UNDER THE GNU PUBLIC LICENSE V3")')
      WRITE(NFILO,*)
      RETURN
      END SUBROUTINE INITIALIZEFILEHANDLER
!
!     ..................................................................     
      SUBROUTINE STANDARDFILES
!     *****************************************************************
!     **                                                             **
!     *****************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(32)        :: CH32SVAR1
      CHARACTER(32)        :: ID
      INTEGER(4)           :: NFILO
!     *****************************************************************
                                   CALL TRACE$PUSH('STANDARDFILES')
!  
!     ==================================================================
!     == SET STANDARD FILENAMES                                       ==
!     ==================================================================
!
!     ==  ERROR FILE ===================================================
      ID=+'ERR'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.DERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PROTOCOL FILE ================================================
      ID=+'PROT'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.PPROT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  CONTROL FILE  == =============================================
      ID=+'PCNTL'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.PCNTL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  STRUCTURE FILE   =============================================
      ID=+'STRC'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.STRC')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  XYZ FILE ====================================================
      ID=+'XYZ'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.XYZ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
                                   CALL TRACE$POP
      RETURN
      END SUBROUTINE STANDARDFILES

!
!     ..................................................................
      SUBROUTINE TIMESTEP$GETI4(ID_,VAL_)
!     ******************************************************************
!     ** THIS IS A DUMMY ROUTINE IN ORDER TO STAY COMPATIBLE WITH     **
!     ** PAW.F                                                        **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(OUT):: VAL_
!     ******************************************************************
      VAL_=1
      RETURN
    END SUBROUTINE TIMESTEP$GETI4
!
!     ..................................................................
      SUBROUTINE MM_HESSIAN(NDIM,R0,HESSIAN)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE THE HESSIAN FROM THE FORCE FIELD                  **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)   :: NDIM
      REAL(8)     ,INTENT(IN)   :: R0(NDIM)
      REAL(8)     ,INTENT(OUT)  :: HESSIAN(NDIM,NDIM)
      REAL(8)                   :: E,DELTA
      REAL(8)                   :: R(NDIM),FORCE1(NDIM),FORCE2(NDIM)
      INTEGER(4)                :: IDIM,IDIM2
!     ******************************************************************
      DELTA=1.D-3
      DO IDIM=1,NDIM
         R(:)=R0(:)
         R(IDIM)=R(IDIM)-DELTA
         CALL CLASSICAL$SETR8A('R(0)',NDIM,R)
         CALL CLASSICAL$ETOT(E)
         CALL CLASSICAL$GETR8A('FORCE',NDIM,FORCE1)
         R(IDIM)=R(IDIM)+2*DELTA
         CALL CLASSICAL$SETR8A('R(0)',NDIM,R)
         CALL CLASSICAL$ETOT(E)
         CALL CLASSICAL$GETR8A('FORCE',NDIM,FORCE2)
         DO IDIM2=1,NDIM
            HESSIAN(IDIM,IDIM2)=(FORCE1(IDIM2)-FORCE2(IDIM2))/(2.D0*DELTA)
         END DO
      END DO
      !SYMMETRIZE
      DO IDIM=1,NDIM
         DO IDIM2=1,NDIM
            HESSIAN(IDIM,IDIM2)=0.5D0*(HESSIAN(IDIM,IDIM2)+HESSIAN(IDIM2,IDIM))
            HESSIAN(IDIM2,IDIM)=HESSIAN(IDIM,IDIM2)
         END DO
      END DO
      END SUBROUTINE MM_HESSIAN
!
!     ..................................................................
      SUBROUTINE MM_ETOT_FORCE(NDIM,R,E,FORCE)
!     ******************************************************************
!     **                                                              **
!     **  WRAPPER FOR THE DIFFERENT OPTIMIZATION ROUTINES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)   :: NDIM
      REAL(8)     ,INTENT(IN)   :: R(NDIM)
      REAL(8)     ,INTENT(OUT)  :: E
      REAL(8)     ,INTENT(OUT)  :: FORCE(NDIM)
!     ******************************************************************
      CALL CLASSICAL$SETR8A('R(0)',NDIM,R)
      CALL CLASSICAL$ETOT(E)
      CALL CLASSICAL$GETR8A('FORCE',NDIM,FORCE)
      END SUBROUTINE MM_ETOT_FORCE
! 
!     ..................................................................
      SUBROUTINE MM_MINIMIZE_CG
!     ******************************************************************
!     **                                                              **
!     **  OPTIMISES ENVIRONMENT ATOMS, EXCLUDING ATOMS OF THE REACTION**
!     **  CENTER AND LINK ATOMS, USING A CONJUGATE GRADIENT MINIMIZER **
!     **                                                              **
!     ******************************************************************
      USE POPT_MODULE
!      USE TESTPOT_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      REAL(8)               :: TOLI=1.D-5
      INTEGER(4)            :: NITERI=100
      REAL(8)               :: FORCE1(3*NATM)
      REAL(8)               :: FORCE2(3*NATM)
      REAL(8)               :: R1(3*NATM)
      REAL(8)               :: R2(3*NATM)
      REAL(8)               :: HESSIAN(3*NATM,3*NATM)
      REAL(8)               :: PRECON(3*NATM)
      INTEGER(4)            :: ITER,IATM
      REAL(8)               :: FMAX
      REAL(8)               :: ALPHA,SVAR,ALPHA_TEST
      REAL(8)               :: SDIR(3*NATM),SDIROLD(3*NATM)
      REAL(8)               :: ABSF1TOT,ABSF1TOTOLD
      REAL(8)               :: ABSSDIR,ABSF1,PALPHA,ABSF2,EPOT1,EPOT2
      INTEGER(4)            :: NFILO,I1,I2,INNER
      LOGICAL(4)            :: TCONV
!     ******************************************************************
      CALL CLASSICAL$SELECT('QMMM')
!     ==================================================================
!     ==  ITERATE                                                     ==
!     ==================================================================
      ALPHA=1.D-1
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,*)' OPITMIZATION '
      WRITE(NFILO,*)'=============='
      WRITE(NFILO,FMT='(55("."),": ",T1,A,T58,F10.6,A)') &
           'TOLERANCE IN MAXIMAL FORCE',TOL
      WRITE(NFILO,FMT='(55("."),": ",T1,A,T58,I10,A)') &
           'MAXIMAL NUMBER OF ITERATIONS',NITER
      WRITE(NFILO,*)
      WRITE(NFILO,*)'ITER     EPOT      ALPHA     FMAX '
!
!     ================================================================
!     ==  GET POSITIONS AND FORCES                                  ==
!     ================================================================
      CALL CLASSICAL$GETR8A('R(0)',3*NATM,R1)
!
!     ================================================================
!     ==  CALCULATE PRECONDITIONING                                 ==
!     ================================================================
      CALL MM_HESSIAN(3*NATM,R1,HESSIAN)
      DO IATM=1,3*NATM
         PRECON(IATM)=1.D0/HESSIAN(IATM,IATM)
         PRECON(IATM)=ABS(PRECON(IATM)) ! THIS LINE MAY BE DANGEROUS
      END DO

      SDIROLD=0.D0
      DO ITER=1,NITER
        IF(MOD(ITER-1,100).EQ.0) THEN
           CALL CLASSICAL$NEIGHBORS
           !RECALCULATE PRECONDITIONING
!!$           CALL MM_HESSIAN(3*NATM,R1,HESSIAN)
!!$           DO IATM=1,3*NATM
!!$              PRECON(IATM)=1.D0/HESSIAN(IATM,IATM)
!!$              PRECON(IATM)=ABS(PRECON(IATM)) ! THIS LINE MAY BE DANGEROUS
!!$           END DO
        END IF
        CALL MM_ETOT_FORCE(3*NATM,R1,EPOT1,FORCE1)
        FORCE1(:)=PRECON(:)*FORCE1(:)
        ABSF1TOT=SQRT(DOT_PRODUCT(FORCE1,FORCE1))
!
!       ================================================================
!       ==  DEFINE THE NEW SEARCH DIRECTION                           ==
!       ================================================================
        IF(ITER.EQ.1) THEN
           SDIR=FORCE1
        ELSE
           SDIR=FORCE1+SDIROLD*(ABSF1TOT/ABSF1TOTOLD)**2
        END IF
        !FORCE1 PROJECTED ONTO SDIR
        ABSSDIR=SQRT(DOT_PRODUCT(SDIR,SDIR))
        ABSF1=DOT_PRODUCT(SDIR,FORCE1)/ABSSDIR
!
!       ================================================================
!       ==  PERFORM LINE SEARCH                                       ==
!       ================================================================
        DO INNER=1,NITERI
          R2(:)=R1(:)+ALPHA*SDIR(:)
          CALL MM_ETOT_FORCE(3*NATM,R2,EPOT2,FORCE2)
          FORCE2(:)=PRECON(:)*FORCE2(:)
          ABSF2=DOT_PRODUCT(SDIR,FORCE2)/ABSSDIR
          !EXTRAPOLATION BASED ON F1F2
          ALPHA=ABSF1/(ABSF1-ABSF2)*ALPHA
          IF(ALPHA.LT.0.D0) ALPHA=1.D-3
          IF(ABS(ABSF2).LT.TOLI) EXIT
        ENDDO
        R2(:)=R1(:)+ALPHA*SDIR(:)
!
!       ================================================================
!       ==  PROPAGATE                                                 ==
!       ================================================================
        ABSF1TOTOLD=ABSF1TOT
        SDIROLD=SDIR
        R1(:)=R2(:)
!
!       ================================================================
!       ==  PRINT                                                     ==
!       ================================================================
        !FMAX=ABSF1TOT !BAD WITH PRECONDINTIONING
        FMAX=SQRT(DOT_PRODUCT(FORCE2(:)/PRECON(:),FORCE2(:)/PRECON(:)))
        IF(MOD(ITER,100).EQ.0.OR.ITER.EQ.1) THEN
           WRITE(NFILO,FMT='(I5,1X,F10.5,F10.5,F10.5)') ITER,EPOT1,ALPHA,FMAX
           FLUSH(NFILO)
        END IF
        TCONV=(FMAX.LT.TOL)
        IF(TCONV) EXIT
      ENDDO
      CALL CLASSICAL$SETR8A('R(0)',3*NATM,R1)
      WRITE(NFILO,"('MINIMIZATION FINISHED AFTER ',I10,' ITERATIONS')") ITER-1
      WRITE(*,"('MINIMIZATION FINISHED AFTER ',I10,' ITERATIONS')") ITER-1
      IF(TCONV) THEN
         WRITE(NFILO,"('MINIMIZATION CONVERGED')")
      ELSE
         WRITE(NFILO,"('MINIMIZATION NOT CONVERGED, HOWEVER THE RESULTS MAY BE USEFUL')")
      END IF
      WRITE(NFILO,*)
    END SUBROUTINE MM_MINIMIZE_CG
! 
!     ..................................................................
      SUBROUTINE MM_MINIMIZE_DYNAMIC
!     ******************************************************************
!     **                                                              **
!     **  OPTIMISES ENVIRONMENT ATOMS, EXCLUDING ATOMS OF THE REACTION**
!     **  CENTER AND LINK ATOMS, USING MOLECULAR DYNAMICS             **
!     **                                                              **
!     ******************************************************************
      USE POPT_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      REAL(8)   ,PARAMETER  :: DELT=10.D0
      REAL(8)   ,PARAMETER  :: EKINMAX=0.5D0
      REAL(8)               :: FORCE(3*NATM)
      REAL(8)               :: RM(3*NATM)
      REAL(8)               :: R0(3*NATM)
      REAL(8)               :: RP(3*NATM)
      INTEGER(4)            :: ITER,IMAP,ILINK
      INTEGER(4)            :: IATM,IATMJ
      REAL(8)               :: EPOTP,EPOT,EKIN,ECON,FMAX
      REAL(8)               :: ANNELOC
      REAL(8)               :: ALPHA,SVAR,VEC(3)
      INTEGER(4)            :: NFILO,I1,I2,INNER,IAT
      LOGICAL(4)            :: TCONV
      REAL(8)               :: FRIC,DT,LOWFRIC,HIGHFRIC,FAC,EPOT_OLD
      REAL(8)               :: MASS(3*NATM)
      REAL(8)               :: M(NATM)
!     ******************************************************************
      DT=5.D1  !1.D2
      LOWFRIC=0.3D0
      HIGHFRIC=0.3D0
      FAC=0.9D0
      CALL CLASSICAL$SELECT('QMMM')
      EPOTP=1.D+8
      ALPHA=1.D-4
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,*)' OPITMIZATION '
      WRITE(NFILO,*)'=============='
!      WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,F10.5,A)')NAME,VALUE,UNIT
      WRITE(NFILO,FMT='(55("."),": ",T1,A,T58,F10.6,A)') &
           'TOLERANCE IN MAXIMAL FORCE',TOL
      WRITE(NFILO,FMT='(55("."),": ",T1,A,T58,I10,A)') &
           'MAXIMAL NUMBER OF ITERATIONS',NITER
      WRITE(NFILO,*)
      WRITE(NFILO,*)'ITER     EPOT      FRIC      FMAX '


      CALL CLASSICAL$GETR8A('R(0)',3*NATM,R0)
      RM=R0
!
!     ================================================================
!     ==  GET MASSES AND CONVERT THEM TO 3*NATM                     ==
!     ================================================================
      CALL CLASSICAL$GETR8A('MASS',NATM,M)
      DO IAT=1,NATM
         MASS(3*(IAT-1)+1)=M(IAT)
         MASS(3*(IAT-1)+2)=M(IAT)
         MASS(3*(IAT-1)+3)=M(IAT)
      END DO
      EPOT_OLD=1.20
!
!     ==================================================================
!     ==  ITERATE                                                     ==
!     ==================================================================
      DO ITER=1,NITER
        IF(MOD(ITER-1,100).EQ.0) CALL CLASSICAL$NEIGHBORS
!
!       ================================================================
!       ==  CALCULATE FORCE                                           ==
!       ================================================================
        CALL CLASSICAL$ETOT(EPOT)
        CALL CLASSICAL$GETR8A('FORCE',3*NATM,FORCE)
!
!       ================================================================
!       ==  SET FRICTION                                              ==
!       ================================================================
        IF(EPOT.GT.EPOT_OLD) THEN
           FRIC=HIGHFRIC
        ELSE
           FRIC=LOWFRIC*FAC
           LOWFRIC=LOWFRIC*FAC
        END IF
        EPOT_OLD=EPOT
!
!       ================================================================
!       ==  CALCULATE NEXT STEP                                       ==
!       ================================================================
        RP=(2.D0*R0-(1.D0-FRIC)*RM+FORCE*DT**2/MASS)/(1.D0+FRIC)
!
!       ================================================================
!       ==  PROPAGATE                                                 ==
!       ================================================================
        RM=R0
        R0=RP
        CALL CLASSICAL$SETR8A('R(0)',3*NATM,R0)

        FMAX=SQRT(DOT_PRODUCT(FORCE,FORCE))
        IF(MOD(ITER,100).EQ.0.OR.ITER.EQ.1) THEN
           WRITE(NFILO,FMT='(I5,1X,F10.5,F10.5,F10.5)') ITER,EPOT,FRIC,FMAX
        END IF
        TCONV=(FMAX.LT.TOL)
        IF(TCONV) EXIT
      ENDDO
      WRITE(NFILO,"('MINIMIZATION FINISHED AFTER ',I10,' ITERATIONS')") ITER-1
      WRITE(*,"('MINIMIZATION FINISHED AFTER ',I10,' ITERATIONS')") ITER-1
      IF(TCONV) THEN
         WRITE(NFILO,"('MINIMIZATION CONVERGED')")
      ELSE
         WRITE(NFILO,"('MINIMIZATION NOT CONVERGED, HOWEVER THE RESULTS MAY BE USEFUL')")
      END IF
      WRITE(NFILO,*)
    END SUBROUTINE MM_MINIMIZE_DYNAMIC


