!
!............................................ENERGYLIST.................
MODULE ENERGYLIST_MODULE
!**                                                                   ** 
!** KEEPS TRACK OF ENERGIES                                           ** 
!**                                                                   ** 
!** FUNCTIONS:                                                        ** 
!**   SETUNIT                                                         ** 
!**   RESET                                                           ** 
!**   SET                                                             ** 
!**   ADD                                                             ** 
!**   RETURN                                                          ** 
!**   PRINT                                                           ** 
!**                                                                   ** 
!**                                                                   ** 
IMPLICIT NONE
INTEGER(4),PARAMETER :: NEX=100
INTEGER(4)           :: NE=0
CHARACTER(40)        :: IDENTIFIER(NEX)
REAL(8)              :: ENERGY(NEX)
CHARACTER(10)        :: UNITNAME='H'
REAL(8)              :: UNITVALUE=1.D0
END MODULE ENERGYLIST_MODULE
!
!     ..................................................................
      SUBROUTINE ENERGYLIST$SETUNIT(UNITNAME_,UNITVALUE_)
!     ==================================================================
!     == SPECIFIES THE ENERGY UNIT                                    ==
!     == (EXPRESSED IN STANDARD UNITS(HARTREE))                       ==
!     == DEFAULT IS HARTREE ATOMIC UNITS                              ==
!     ==================================================================
      USE ENERGYLIST_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: UNITNAME_
      REAL(8)     ,INTENT(IN) :: UNITVALUE_
!     ******************************************************************
      UNITNAME=UNITNAME_
      UNITVALUE=UNITVALUE_
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ENERGYLIST$RESET
!     ==================================================================
!     == SETS ALL ENERGIES TO ZERO                                    ==
!     ==================================================================
      USE ENERGYLIST_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: I
!     ******************************************************************
      NE=0
      DO I=1,NEX
        IDENTIFIER(I)=' '
        ENERGY(I)=0.D0
      ENDDO
      RETURN
      END
!     
!     ..................................................................
      SUBROUTINE ENERGYLIST$SET(STRING_,VALUE_)
!     ==================================================================
!     == SETS ENERGY                                                  ==
!     == STRING  IDENTIFIES THE ENERGY                                ==
!     == VALUE   IS THE VALUE OF THE ENERGY                           ==
!     ==================================================================
      USE ENERGYLIST_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: STRING_
      REAL(8)     ,INTENT(IN) :: VALUE_
      CHARACTER(40)           :: STRING
      LOGICAL(4)              :: TSET
      INTEGER(4)              :: I
!     ******************************************************************
      TSET=.FALSE.
      STRING=STRING_
      DO I=1,NE
        IF(STRING.EQ.IDENTIFIER(I)) THEN
          ENERGY(I)=VALUE_
          TSET=.TRUE.
          EXIT
        END IF
      ENDDO 
      IF(.NOT.TSET) THEN
        IF(NE.EQ.NEX) THEN
          CALL ERROR$MSG('ATTEMPT TO STORE MORE ENERGIES THAN ALLOWED')
          CALL ERROR$MSG('INCREASE PARAMETER NEX IN ROUTINE ENERGIES')
          CALL ERROR$STOP('ENERGYLIST$SET')
        END IF
        NE=NE+1
        IDENTIFIER(NE)=STRING
        ENERGY(NE)=VALUE_
      END IF    
      RETURN
      END
!     
!     ..................................................................
      SUBROUTINE ENERGYLIST$ADD(STRING_,VALUE_)
!     ==================================================================
!     == SETS ENERGY                                                  ==
!     == STRING  IDENTIFIES THE ENERGY                                ==
!     == VALUE   IS THE VALUE OF THE ENERGY                           ==
!     ==================================================================
      USE ENERGYLIST_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: STRING_
      REAL(8)     ,INTENT(IN) :: VALUE_
      CHARACTER(40)           :: STRING
      INTEGER(4)              :: I
      LOGICAL(4)              :: TSET
!     ******************************************************************
      TSET=.FALSE.
      STRING=STRING_
      DO I=1,NE
        IF(STRING.EQ.IDENTIFIER(I)) THEN
          ENERGY(I)=ENERGY(I)+VALUE_
          TSET=.TRUE.
          EXIT
        END IF
      ENDDO 
      IF(.NOT.TSET) THEN
        IF(NE.EQ.NEX) THEN
          CALL ERROR$MSG('ATTEMPT TO STORE MORE ENERGIES THAN ALLOWED')
          CALL ERROR$MSG('INCREASE PARAMETER NEX IN ROUTINE ENERGIES')
          CALL ERROR$STOP('ENERGYLIST$ADD')
        END IF
        NE=NE+1
        IDENTIFIER(NE)=STRING
        ENERGY(NE)=VALUE_
      END IF    
      RETURN
      END
!     
!     ..........................................................................
      SUBROUTINE ENERGYLIST$RETURN(ID,VALUE_)
!     **************************************************************************
!     ** THIS FUNCTION IS OBSOLETE AND WILL BE REPLACED BY ENERGYLIST$GET     **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN)  :: ID
      REAL(8)      ,INTENT(OUT) :: VALUE_
!     **************************************************************************
      CALL ERROR$MSG('SUBROUTINE MARKED FOR DELETION')
      CALL ERROR$STOP('ENERGYLIST$RETURN')
!
      CALL ENERGYLIST$GET(ID,VALUE_)
      RETURN
      END
!     
!     ..........................................................................
      SUBROUTINE ENERGYLIST$GET(ID,VAL)
!     **************************************************************************
!     ** RETURNS ENERGY VALUE                                                 **
!     **************************************************************************
      USE ENERGYLIST_MODULE, ONLY : NE,IDENTIFIER,ENERGY
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN)  :: ID
      REAL(8)      ,INTENT(OUT) :: VAL
      INTEGER(4)                :: I
!     **************************************************************************
      VAL=0.D0
      DO I=1,NE
        IF(ID.EQ.IDENTIFIER(I)) THEN
          VAL=ENERGY(I)
          RETURN
        END IF
      ENDDO
!!$!     __AN ITEM THAT HAS NOT BEEN SET OBTAINS THE VALUE ZERO
!!$      CALL ERROR$MSG('LIST ITEM NOT FOUND')
!!$      CALL ERROR$CHVAL('ID',ID)
!!$      CALL ERROR$STOP('ENERGYLIST$GET')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ENERGYLIST$PRINT(NFIL)
!     ==================================================================
!     == PRINTS ALL STORED ENERGIES                                   ==
!     ==================================================================
      USE ENERGYLIST_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: I
!     ******************************************************************
      CALL REPORT$TITLE(NFIL,'ENERGY REPORT')
      DO I=1,NE
        WRITE(NFIL,FMT='(A40,":",F15.7," ",A10)') &
     &        IDENTIFIER(I),ENERGY(I)/UNITVALUE,UNITNAME
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE  ENERGYLIST$PRINTHEADER(NFIL)
!     ==================================================================
!     == PRINTS ALL STORED ENERGIES                                   ==
!     ==================================================================
      USE ENERGYLIST_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
!     ******************************************************************
      WRITE(NFIL,FMT='("ENERGY REPORT"/"=============")')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ENERGYLIST$PRINTONE(NFIL,STRING_)
!     ==================================================================
!     == PRINTS ALL STORED ENERGIES                                   ==
!     ==================================================================
      USE ENERGYLIST_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NFIL
      CHARACTER(*) ,INTENT(IN) :: STRING_
      CHARACTER(40)            :: STRING
      INTEGER(4)               :: I
!     ******************************************************************
      STRING=STRING_
      DO I=1,NE
        IF(IDENTIFIER(I).EQ.STRING) THEN
          WRITE(NFIL,FMT='(A40,":",F15.7," ",A10)') &
     &        IDENTIFIER(I),ENERGY(I)/UNITVALUE,UNITNAME
        END IF
      ENDDO
      RETURN
      END
!
!     ..................................................................
      MODULE GROUPLIST_MODULE
      LOGICAL(4)                :: TINI=.FALSE.
      INTEGER(4)                :: NATOM=0
      INTEGER(4)                :: NGROUPX=0
      INTEGER(4)                :: NGROUP=0
      CHARACTER(32),ALLOCATABLE :: IDENTIFIER(:) ! NGROUPX
      LOGICAL(4)   ,ALLOCATABLE :: TMEMBER(:,:)
      END MODULE GROUPLIST_MODULE
!
!     ..................................................................
      SUBROUTINE GROUPLIST$INITIALIZE(NATOM1,NGROUPX1)
!     ******************************************************************
!     **  INITIALIZE                                                  **
!     ******************************************************************
      USE GROUPLIST_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NATOM1
      INTEGER(4),INTENT(IN)  :: NGROUPX1
      INTEGER(4)             :: IGROUP,IAT
!     ******************************************************************
      IF(TINI) THEN
        CALL ERROR$MSG('GROUPLIST  INITIALIZED')
        CALL ERROR$STOP('GROUPLIST')
      END IF
      TINI=.TRUE.
      NATOM=NATOM1
      NGROUPX=NGROUPX1
!     CALL MEMORY$ALLOCATE(32,$IDENTIFIER,NGROUPX)
!     CALL MEMORY$ALLOCATE(4,$TMEMBER,NATOM*NGROUPX)
      ALLOCATE(IDENTIFIER(NGROUPX))
      ALLOCATE(TMEMBER(NATOM,NGROUPX))
      NGROUP=0
      DO IGROUP=1,NGROUPX
        IDENTIFIER(IGROUP)=' '
        DO IAT=1,NATOM
          TMEMBER(IAT,IGROUP)=.FALSE.
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GROUPLIST$ADD(IDENT1,IAT1)
!     ******************************************************************
!     **  ADD ATOM TO THE GROUP                                       **
!     ******************************************************************
      USE GROUPLIST_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: IDENT1
      INTEGER(4),INTENT(IN)  :: IAT1
      INTEGER(4)             :: IGROUP,I
!     ******************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('GROUPLIST NOT INITIALIZED')
        CALL ERROR$STOP('GROUPLIST$ADD')
      END IF
!
!     ==================================================================
!     ==  DETERMINE GROUP INDEX                                       ==
!     ==================================================================
      IGROUP=0
      DO I=1,NGROUP
        IF(TRIM(IDENT1).EQ.TRIM(IDENTIFIER(I))) THEN
          IGROUP=I
          EXIT
        END IF
      ENDDO
      IF(IGROUP.EQ.0) THEN
        NGROUP=NGROUP+1
        IF(NGROUP.GT.NGROUPX) THEN
          CALL ERROR$MSG('MORE GROUPS THAN ALLOWED BY INITIALIZATION')
          CALL ERROR$MSG('CURRENT MAX. NUMBER OF GROUPS=NGROUPX')
          CALL ERROR$I4VAL('NGROUPX',NGROUPX)
          CALL ERROR$STOP('GROUPLIST$ADD')
        ENDIF
        IGROUP=NGROUP
        IDENTIFIER(NGROUP)=IDENT1
      END IF
      IF(IAT1.GT.NATOM) THEN
        CALL ERROR$MSG('ATOMINDEX LARGER THAN MAXIMUM')
        CALL ERROR$STOP('GROUPLIST$ADD')
      END IF
      TMEMBER(IAT1,IGROUP)=.TRUE.
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GROUPLIST$MEMBERS(IDENT1,NATOM1,TMEMBER1)
!     ******************************************************************
!     **  REQUEST FOR MEMBERS OF A GIVEN GROUP                        **
!     ******************************************************************
      USE GROUPLIST_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: IDENT1
      INTEGER(4)  ,INTENT(IN) :: NATOM1
      LOGICAL(4)  ,INTENT(OUT):: TMEMBER1(NATOM1)
      INTEGER(4)              :: IGROUP,IAT,I
!     ******************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('GROUPLIST NOT INITIALIZED')
        CALL ERROR$STOP('GROUPLIST$MEMBERS')
      END IF
      IF(NATOM1.NE.NATOM) THEN
        CALL ERROR$MSG(' NUMBER OF ATOMS INCONSISTENT')
        CALL ERROR$STOP('GROUPLIST$MEMBERS')
      END IF
!
!     ==================================================================
!     ==  DETERMINE GROUP INDEX                                       ==
!     ==================================================================
      IGROUP=0
      DO I=1,NGROUP
        IF(TRIM(IDENT1).EQ.TRIM(IDENTIFIER(I))) THEN
          IGROUP=I
          EXIT
        END IF
      ENDDO
      IF(IGROUP.EQ.0) THEN
        CALL ERROR$MSG('NO GROUP NAMED BY IDENT IN THE LIST')
        CALL ERROR$CHVAL('IDENT',IDENT1)
        CALL ERROR$STOP('GROUPLIST')
      END IF
      DO IAT=1,NATOM
        TMEMBER1(IAT)=TMEMBER(IAT,IGROUP)
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GROUPLIST$LENGTH(NGROUP1)
!     ******************************************************************
!     **  LENGTH                                                      **
!     ******************************************************************
      USE GROUPLIST_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: NGROUP1
!     ******************************************************************
      NGROUP1=NGROUP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GROUPLIST$MEMBERNUMBER(IDENT1,NMEMBER)
!     ******************************************************************
!     **  NUMBER OF MEMBERS                                           **
!     ******************************************************************
      USE GROUPLIST_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: IDENT1
      INTEGER(4)  ,INTENT(OUT) :: NMEMBER
      INTEGER(4)               :: IGROUP,IAT,I
!     ******************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('GROUPLIST NOT INITIALIZED')
        CALL ERROR$STOP('GROUPLIST')
      END IF
!
!     ==================================================================
!     ==  DETERMINE GROUP INDEX                                       ==
!     ==================================================================
      IGROUP=0
      DO I=1,NGROUP
        IF(TRIM(IDENT1).EQ.TRIM(IDENTIFIER(I))) THEN
          IGROUP=I
          EXIT
        END IF
      ENDDO
      IF(IGROUP.EQ.0) THEN
        CALL ERROR$MSG('NO GROUP NAMED BY IDENT IN THE LIST')
        CALL ERROR$CHVAL('IDENT',IDENT1)
        CALL ERROR$STOP('GROUPLIST')
      END IF
!
!     ==================================================================
!     ==  COUNT MEMBERS IN THE GROUP                                  ==
!     ==================================================================
      NMEMBER=0
      DO IAT=1,NATOM
        IF(TMEMBER(IAT,IGROUP)) NMEMBER=NMEMBER+1
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GROUPLIST$REPORT(NFIL)
!     ******************************************************************
!     **  REPORT                                                      **
!     ******************************************************************
      USE GROUPLIST_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      INTEGER(4)              :: IGROUP,IAT,NMEM,I
      CHARACTER(32)           :: NAMES(NATOM)
!     ******************************************************************
      CALL REPORT$TITLE(NFIL,'GROUPLIST REPORT')
      DO IGROUP=1,NGROUP
        NMEM=0
        DO IAT=1,NATOM
          IF(TMEMBER(IAT,IGROUP)) THEN
            NMEM=NMEM+1
            CALL ATOMLIST$GETCH('NAME',IAT,NAMES(NMEM))
          END IF
        ENDDO
        WRITE(NFIL,FMT='("GROUP:",A16,"N(MEMBERS)",I5)') &
     &         IDENTIFIER(IGROUP),NMEM
        WRITE(NFIL,FMT='(10(A,"; "))')(TRIM(NAMES(I)),I=1,NMEM)
      ENDDO
      RETURN
      END
!     
!     .................................................................
      SUBROUTINE LIST_LOOKUP(LENGTH,LIST,ITEM,INDEX)
!     **                                                             **
!     **  DETERMNIES THE POSITION OF AN ITEM IN A LIST               **
!     **  SUCH THAT                                                  **
!     **            ITEM=LIST(INDEX)                                 **
!     **                                                             **
      INTEGER(4)  ,INTENT(IN) :: LENGTH
      INTEGER(4)  ,INTENT(OUT):: INDEX
      CHARACTER(*),INTENT(IN) :: ITEM
      CHARACTER(*),INTENT(IN) :: LIST(LENGTH)
      INTEGER(4)              :: I
!     ******************************************************************
      INDEX=0
      DO I=1,LENGTH
        IF(TRIM(ITEM).EQ.TRIM(LIST(I))) THEN
          INDEX=I
          RETURN
        END IF
      ENDDO
      RETURN
      END
!     
!     .................................................................
      SUBROUTINE LIST_TESTINI(TINI,NAME)
!     ******************************************************************
!     ******************************************************************
      LOGICAL(4)  ,INTENT(IN) :: TINI
      CHARACTER(*),INTENT(IN) :: NAME
!     ******************************************************************
      IF(TINI) RETURN
!      
      CALL ERROR$MSG('OBJECT NAME IS NOT INITIALIZED')
      CALL ERROR$CHVAL('NAME',NAME)
      CALL ERROR$STOP('LIST_TESTINI')
      END

