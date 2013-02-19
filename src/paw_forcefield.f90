!     ...1.........2.........3.........4.........5.........6.........7.........8
MODULE FORCEFIELD_MODULE
!*******************************************************************************
!*  FORCEFIELD_MODULE                                                          *
!*                                                                             *
!*  PROVIDES THE INPUT/OUTPUT PROCEDURES TO READ EXTERNAL FORCEFIELD FILES     *
!*  AND STORES THEM IN DATA ARRAYS. THE FOLLOWING FORCEFIELDS ARE SUPPORTED:   *
!*                                                                             *
!*     1.)  AMBER (cornell_all files, including TIP3P water model)             *
!*                                                                             *
!*  FUNCTIONS:                                                                 *
!*    FORCEFIELD$SETCH                                                         *
!*    FORCEFIELD$GETCH                                                         *
!*    FORCEFIELD$GETI4                                                         *
!*    FORCEFIELD$GETRESNUMBER                                                  *
!*    FORCEFIELD$READ_PARMFILE                                                 *
!*    FORCEFIELD$READ_TOPFILE                                                  *
!*    FORCEFIELD$READ_MMSTRC                                                   *
!*    FORCEFIELD$WRITE_MMSTRC                                                  *
!*    FORCEFIELD$ADDATOM        (=> RESNUMBER, ATOM( TYPE(TOP_ATOM_TYPE)  ) )  *
!*    FORCEFIELD$ADDBOND        (=> RESNUMBER, BOND( TYPE(TOP_BOND_TYPE)  ) )  *
!*    FORCEFIELD$ADDIMPROPER    (=> RESNUMBER, IMPROPER( TYPE(TOP_IMPROPER_TYPE))) *
!*    FORCEFIELD$ADDIC          (=> RESNUMBER, IC( TYPE(TOP_ATOM_TYPE)  ) )    *
!*    FORCEFIELD$DELETEATOM     (=> RESNUMBER, ATOMNAME   )                    *
!*                                                                             *
!*  FORCEFIELD AMBER SPECIFIC ROUTINES:                                        *
!*    FORCEFIELD$AMBER_BONDPARMS                                               *
!*    FORCEFIELD$AMBER_ANGLEPARMS                                              *
!*    FORCEFIELD$AMBER_TORSIONPARMS                                            *
!*    FORCEFIELD$AMBER_NONBONDPARMS                                            *
!*                                                                             *
!*******************************************WRITTEN BY SASCHA HEMMEN, 2006******

!----  DATA TYPES FOR THE FORCEFIELD PARAMETERS
  TYPE FF_BOND_TYPE                ! BOND KEYWORD
     character(2)         :: ATOM1
     character(2)         :: ATOM2
     real(8)              :: f_const
     real(8)              :: r_equi
  END TYPE FF_BOND_TYPE
  TYPE FF_ANGLE_TYPE               ! THETAS KEYWORD
     character(2)         :: ATOM1
     character(2)         :: ATOM2
     character(2)         :: ATOM3
     real(8)              :: f_const
     real(8)              :: THETA_equi
  END TYPE FF_ANGLE_TYPE
  TYPE FF_TORSION_TYPE             ! PHI KEYWORD
     character(2)         :: ATOM1
     character(2)         :: ATOM2
     character(2)         :: ATOM3
     CHARACTER(2)         :: ATOM4
     real(8)              :: f_const   ! V_n/2 DIVIDED BY THE PERIODICITY
     integer(4)           :: periodicy ! # of BOND PATHS OF THE MIDDLE ATOMS
     real(8)              :: phase     ! symmetry
  END TYPE FF_TORSION_TYPE
  TYPE FF_IMPTORSION_TYPE          ! IMPHI KEYWORD
     character(2)         :: ATOM1
     character(2)         :: ATOM2
     character(2)         :: ATOM3
     CHARACTER(2)         :: ATOM4
     real(8)              :: f_const 
     integer(4)           :: periodicy
     real(8)              :: phase
  END TYPE FF_IMPTORSION_TYPE
  TYPE FF_VDW_TYPE                 ! NONBONDED KEYWORD
     character(2)         :: ATOM1
     real(8)              :: Emin
     real(8)              :: Rmin
  END TYPE FF_VDW_TYPE


!---- DATA TYPES FOR THE TOPOLGIES
  TYPE TOP_MASSES_TYPE
     integer              :: id
     CHARACTER(2)         :: atom
     real(8)              :: mass
  END TYPE TOP_MASSES_TYPE

  TYPE TOP_IC_TYPE
     character(4)         :: atom1, atom2, atom3, atom4
     real(8)              :: dist12, ang123, tor1234, ang234, dist34
  END TYPE TOP_IC_TYPE

  TYPE TOP_IMPROPER_TYPE
     character(4)         :: atom1, atom2, atom3, atom4
  END TYPE TOP_IMPROPER_TYPE

  TYPE TOP_BOND_TYPE
     character(4)         :: atom1, atom2
  END TYPE TOP_BOND_TYPE

  TYPE TOP_ATOM_TYPE
     character(4)         :: name
     character(2)         :: atom
     real(8)              :: charge
  END TYPE TOP_ATOM_TYPE

  TYPE RES_TYPE
     CHARACTER(4)                                 :: RESNAME
     REAL(8)                                      :: TOTALCHARGE
     TYPE(TOP_ATOM_TYPE),    dimension(:),pointer :: ATOM
     TYPE(TOP_BOND_TYPE),    dimension(:),pointer :: BOND
     TYPE(TOP_IMPROPER_TYPE),dimension(:),pointer :: IMPROPER
     TYPE(TOP_IC_TYPE),      dimension(:),pointer :: IC
  END TYPE RES_TYPE

!---- DATA TYPE FOR THE ATOMS AND HETATOMS
  TYPE PDB_ATOM_TYPE
      CHARACTER(6)         :: KEYWORD     ! ATOM OR HETATM. ATOM BELONGS TO AMINOACIDS, DIFFERENT AMINO ACIDS
                                          ! BECOME CONNECTED VIA PEPTIDE BONDS
     INTEGER              :: ID          ! atom serial number
     CHARACTER(5)         :: NAME        ! atom name != element name
!     CHARACTER(1)         :: altLOC      ! alternate location indicator (need?) drop this for "reduce"
     CHARACTER(4)         :: RESNAME     ! residue name. 
     CHARACTER(1)         :: CHAINID     ! chain identifier
     INTEGER              :: RESSEQ      ! Residue Sequence Number
     REAL(8)              :: R(3)        ! xyz - coordinates in angstrom !
     REAL(8)              :: OCCUPANCY
     REAL(8)              :: tempfactor 
     CHARACTER(4)         :: SEGID       ! Segment Identifier
     CHARACTER(2)         :: ELEMENT     ! ELEMENT SYMBOL     BE CAREFUL: element not provided by all PDB files
     CHARACTER(1)         :: FLAG        ! F for fixed, Q for QM ATOM, L for LINKATOM
     CHARACTER(12)        :: QMNAME      ! IF FLAG IST Q OR L THEN YOU CAN FIND THE CORRESPONDING PAW NAME IN THIS COLUMN
     CHARACTER(7)         :: LNAME
  END TYPE PDB_ATOM_TYPE

!---- DATA TYPE FOR THE CONECT DATA
  TYPE PDB_CONECT_TYPE
     CHARACTER(6)         :: KEYWORD
     INTEGER(4)           :: CATOM, ATOM1, ATOM2, ATOM3, ATOM4
  END TYPE PDB_CONECT_TYPE


  CHARACTER(LEN=*), PARAMETER           :: pdb_form='(A6,I5,1X,A5,A4,A1,I4,4X,3F8.3,2F6.2,6X,A4,A2,2X,A1,1X,A12,1X,A7)'
     
  TYPE(FF_BOND_TYPE),       ALLOCATABLE :: bond_parms(:)
  TYPE(FF_ANGLE_TYPE),      ALLOCATABLE :: angle_parms(:)
  TYPE(FF_TORSION_TYPE),    ALLOCATABLE :: torsion_parms(:)
  TYPE(FF_IMPTORSION_TYPE), ALLOCATABLE :: imptorsion_parms(:)
  TYPE(FF_VDW_TYPE),        ALLOCATABLE :: vdw_parms(:)
  INTEGER                               :: NBOND,NANGLE,NTORSION,NIMPTORSION,NVDW

  TYPE(TOP_MASSES_TYPE),    ALLOCATABLE :: masses(:)
  TYPE(RES_TYPE),           ALLOCATABLE :: TOP_RES(:)
  CHARACTER(16)                         :: FORCEFIELD
  CHARACTER(255)                        :: FF_PARMS= '/home/shemmen/PAW/forcefields/amber/cornell_all.prm'
  CHARACTER(255)                        :: FF_TOP=   '/home/shemmen/PAW/forcefields/amber/cornell_all.rtf'
  INTEGER                               :: NRES, NPATCH
  TYPE(PDB_ATOM_TYPE),      ALLOCATABLE :: MMATOM(:)
  TYPE(PDB_CONECT_TYPE),    ALLOCATABLE :: MMCONECT(:)

END MODULE FORCEFIELD_MODULE
!
!     .................................................................. 
      SUBROUTINE FORCEFIELD$SETCH(ID_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID_
        CHARACTER(*),INTENT(IN) :: VAL_
!     ******************************************************************
        IF(ID_.EQ.'FORCEFIELD') THEN
           FORCEFIELD=VAL_
        ELSE IF(ID_.EQ.'PARMFILE') THEN
           FF_PARMS=VAL_
        ELSE IF(ID_.EQ.'TOPFILE') THEN
           FF_TOP=VAL_
        ELSE
           CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
           CALL ERROR$CHVAL('ID_',ID_)
           CALL ERROR$STOP('FORCEFIELD$SETCH')
        END IF
        RETURN
      END SUBROUTINE FORCEFIELD$SETCH
!

!
!     .................................................................. 
      SUBROUTINE FORCEFIELD$GETCH(ID_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID_
        CHARACTER(*),INTENT(OUT):: VAL_
!     ******************************************************************
        IF(ID_.EQ.'FORCEFIELD') THEN
           VAL_=FORCEFIELD
        ELSE IF(ID_.EQ.'PARMFILE') THEN
           VAL_=FF_PARMS
        ELSE IF(ID_.EQ.'TOPFILE') THEN
           VAL_=FF_TOP
        ELSE
           CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
           CALL ERROR$CHVAL('ID_',ID_)
           CALL ERROR$STOP('FORCEFIELD$GETCH')
        END IF
        RETURN
      END SUBROUTINE FORCEFIELD$GETCH
!
!     .................................................................
      SUBROUTINE FORCEFIELD$GETCHA(ID_,LENG_,VAL_)
!     *****************************************************************      
!     **  FORCEFIELD$GET                                             **      
!     *****************************************************************      
        USE FORCEFIELD_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID_
        INTEGER(4)  ,INTENT(IN) :: LENG_
        CHARACTER(*),INTENT(OUT):: VAL_(LENG_)
!     *****************************************************************      
        IF(ID_.EQ.'RESNAME') THEN
           IF(LENG_.NE.SIZE(MMATOM)) THEN
              CALL ERROR$MSG('INCONSISTENT SIZE')
              CALL ERROR$CHVAL('ID_',ID_)
              CALL ERROR$I4VAL('LENG_',LENG_)
              CALL ERROR$STOP('FORCEFIELD$GETCHA')
           END IF
           VAL_=MMATOM(:)%RESNAME
        ELSE IF(ID_.EQ.'NAME') THEN
           IF(LENG_.NE.SIZE(MMATOM)) THEN
              CALL ERROR$MSG('INCONSISTENT SIZE')
              CALL ERROR$CHVAL('ID_',ID_)
              CALL ERROR$I4VAL('LENG_',LENG_)
              CALL ERROR$STOP('FORCEFIELD$GETCHA')
           END IF
           VAL_=MMATOM(:)%NAME       
        ELSE
           CALL ERROR$MSG('INVALID IDENTIFIER')
           CALL ERROR$STOP('FORCEFIELD$GETCHA')
        END IF
        RETURN
      END SUBROUTINE FORCEFIELD$GETCHA
!
!     .................................................................. 
      SUBROUTINE FORCEFIELD$GETI4(ID_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID_
        INTEGER(4),  INTENT(OUT):: VAL_
!     ******************************************************************
        IF(ID_.EQ.'NMMATOM') THEN
           VAL_= SIZE(MMATOM)
        ELSE
           CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
           CALL ERROR$CHVAL('ID_',ID_)
           CALL ERROR$STOP('FORCEFIELD$GETI4')
        END IF
        RETURN
      END SUBROUTINE FORCEFIELD$GETI4
!
!     .................................................................. 
      SUBROUTINE FORCEFIELD$GETI4A(ID_,LENG_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID_
        INTEGER(4)  ,INTENT(IN) :: LENG_
        INTEGER(4)  ,INTENT(OUT):: VAL_(LENG_)
!     ******************************************************************
        IF(ID_.EQ.'RESSEQ') THEN
           IF(LENG_.NE.SIZE(MMATOM)) THEN
              CALL ERROR$MSG('INCONSISTANT SIZE')
              CALL ERROR$CHVAL('ID_',ID_)
              CALL ERROR$I4VAL('LENG_',LENG_)
              CALL ERROR$STOP('FORCEFIELD$GETI4A')
           END IF
           VAL_=MMATOM(:)%RESSEQ
        ELSE
           CALL ERROR$MSG('INVALID IDENTIFIER')
           CALL ERROR$STOP('FORCEFIELD$GETI4A')
        END IF
        RETURN
      END SUBROUTINE FORCEFIELD$GETI4A
!
!     .................................................................. 
      SUBROUTINE FORCEFIELD$GETRESNUMBER(ID_,VAL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN)  :: ID_
        INTEGER(4),  INTENT(OUT) :: VAL_
!     ******************************************************************
        DO VAL_=1,SIZE(TOP_RES)
           IF(TOP_RES(VAL_)%RESNAME.EQ.ID_) RETURN
        ENDDO
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('FORCEFIELD$GETRESNUMBER')
      END SUBROUTINE FORCEFIELD$GETRESNUMBER

!     .................................................................. 
      SUBROUTINE FORCEFIELD$READ_PARMFILE
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        USE STRINGS_MODULE
        IMPLICIT NONE
        integer                               :: nfil,istep
        CHARACTER(255)                        :: line
        LOGICAL(4)                            :: OK1, OK2
        CHARACTER(8)                          :: ID='FF-PARMS'
!     ******************************************************************
                           CALL TRACE$PUSH('READ_PARMFILE')
        CALL FILEHANDLER$SETFILE(ID,.FALSE.,FF_PARMS)
        CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
        
        CALL FILEHANDLER$UNIT(ID,NFIL)
        IF(FORCEFIELD.EQ.'AMBER') THEN
           OK1=.TRUE.
           
           DO WHILE(OK1)
              READ(NFIL,FMT='(A255)', END=101) LINE
              IF(LINE(1:5).EQ.'BONDS') THEN
                 NBOND=0
                 DO 
                    READ(NFIL,'(A255)') LINE
                    IF(LINE(1:1).EQ.'!') CYCLE
                    IF(LINE(1:1).EQ.' ') EXIT
                    NBOND=NBOND+1
                 ENDDO
              END IF
              IF(LINE(1:6).EQ.'THETAS') THEN
                 NANGLE=0
                 DO 
                    READ(NFIL,'(A255)') LINE
                    IF(LINE(1:1).EQ.'!') CYCLE
                    IF(LINE(1:1).EQ.' ') EXIT
                    NANGLE=NANGLE+1
                 ENDDO
              END IF
              IF(LINE(1:3).EQ.'PHI') THEN
                 NTORSION=0
                 DO 
                    READ(NFIL,'(A255)') LINE
                    IF(LINE(1:1).EQ.'!') CYCLE
                    IF(LINE(1:1).EQ.' ') EXIT
                    NTORSION=NTORSION+1
                 ENDDO
              END IF
              IF(LINE(1:5).EQ.'IMPHI') THEN
                 NIMPTORSION=0
                 DO 
                    READ(NFIL,'(A255)') LINE
                    IF(LINE(1:1).EQ.'!') CYCLE
                    IF(LINE(1:1).EQ.' ') EXIT
                    NIMPTORSION=NIMPTORSION+1
                 ENDDO
              END IF
              IF(LINE(1:9).EQ.'NONBONDED') THEN
                 NVDW=0
                 READ(NFIL,'(A)') LINE
                 DO 
                    READ(NFIL,'(A255)') LINE
                    IF(+LINE(1:3).EQ.'END') EXIT
                    IF(LINE(1:1).EQ.'!') CYCLE
                    IF(LINE(1:1).EQ.' ') EXIT
                    NVDW=NVDW+1
!                    print*,nvdw,line(1:65)
                 ENDDO
              END IF
           ENDDO
101        OK1=.FALSE.
           WRITE(*,FMT='(A,I4,3XA,I4,3X,A,I4,3X,A,I4,3X,A,I4)') "#NBOND=",NBOND,"#NANGLE=",NANGLE,&
                "#NTORSION=",NTORSION,"#NIMPTORSION=",NIMPTORSION,"#NVDW=",NVDW
!STOP 'FORCED STOP in READ PARMFILE'           
           ALLOCATE(BOND_PARMS(NBOND))
           ALLOCATE(ANGLE_PARMS(NANGLE))
           ALLOCATE(TORSION_PARMS(NTORSION))
           ALLOCATE(IMPTORSION_PARMS(NIMPTORSION))
           ALLOCATE(VDW_PARMS(NVDW))
           
           REWIND(NFIL)
           DO 
              READ(NFIL,'(A)') LINE
              IF(+LINE(1:3).EQ.'END') EXIT
              
              IF(LINE(1:5).EQ.'BONDS'.AND.LINE(1:1).NE.' ') THEN
                 DO ISTEP=1,NBOND
                    READ(NFIL,FMT='(A2,3X,A2,F7.1,3X,F7.4)') BOND_PARMS(ISTEP)%ATOM1, BOND_PARMS(ISTEP)%ATOM2, &
                         & BOND_PARMS(ISTEP)%F_CONST, BOND_PARMS(ISTEP)%R_EQUI
                 ENDDO
              END IF
              
              IF(LINE(1:6).EQ.'THETAS') THEN
                 DO ISTEP=1,NANGLE
                    READ(NFIL,FMT='(A2,3X,A2,3X,A2,3X,F6.1,5X,F7.2)') ANGLE_PARMS(ISTEP)%ATOM1, ANGLE_PARMS(ISTEP)%ATOM2, &
                         & ANGLE_PARMS(ISTEP)%ATOM3, ANGLE_PARMS(ISTEP)%F_CONST, ANGLE_PARMS(ISTEP)%THETA_EQUI
                 ENDDO
              END IF
              
              IF(LINE(1:3).EQ.'PHI') THEN
                 DO ISTEP=1,NTORSION
                    READ(NFIL,FMT='(A2,2X,A2,2X,A2,2X,A2,2X,F11.8,5X,I1,3X,F6.2)') TORSION_PARMS(ISTEP)%ATOM1, &
                         & TORSION_PARMS(ISTEP)%ATOM2, TORSION_PARMS(ISTEP)%ATOM3, TORSION_PARMS(ISTEP)%ATOM4, &
                         & TORSION_PARMS(ISTEP)%F_CONST, TORSION_PARMS(ISTEP)%PERIODICY, TORSION_PARMS(ISTEP)%PHASE
                 ENDDO
              END IF
              
              IF(LINE(1:5).EQ.'IMPHI') THEN
                 DO ISTEP=1,NIMPTORSION
                    READ(NFIL,FMT='(A2,2X,A2,2X,A2,2X,A2,2X,F11.8,5X,I1,3X,F6.2)') IMPTORSION_PARMS(ISTEP)%ATOM1, &
                         & IMPTORSION_PARMS(ISTEP)%ATOM2, IMPTORSION_PARMS(ISTEP)%ATOM3, IMPTORSION_PARMS(ISTEP)%ATOM4, &
                         & IMPTORSION_PARMS(ISTEP)%F_CONST, IMPTORSION_PARMS(ISTEP)%PERIODICY, IMPTORSION_PARMS(ISTEP)%PHASE
                 ENDDO
              END IF
              
              IF(LINE(1:9).EQ.'NONBONDED') THEN
                 OK2=.TRUE.
                 READ(NFIL,'(A)') LINE
                 DO WHILE(OK2)
                    READ(NFIL,'(A)') LINE
                    print*,line(1:60)
                    IF(+LINE(1:3).EQ.'END') EXIT
                    IF(LINE(1:1).NE.' '.AND.LINE(1:1).NE.'!') OK2=.FALSE.
           !         READ(NFIL,'(A)') LINE
                 ENDDO
                 print*,line(1:60)
                 DO ISTEP=1,NVDW
                    READ(LINE,FMT='(A2,5X,F11.6,28X,F8.4)') VDW_PARMS(ISTEP)%ATOM1, VDW_PARMS(ISTEP)%EMIN, &
                         & VDW_PARMS(ISTEP)%RMIN
                    IF(ISTEP.NE.60) THEN
                       READ(NFIL,'(A)') LINE
                    ELSE
                       EXIT
                    END IF
                 ENDDO
              END IF
           END DO
        ELSE
           CALL ERROR$MSG('FORCEFIELD TYPE IS NOT KNOWN')
           CALL ERROR$CHVAL('FORCEFIELD',FORCEFIELD)
           CALL ERROR$STOP('FORCEFIELD$READ_PARMFILE')
        END IF
                         CALL TRACE$POP
!---------------------- DEBUG remove this later-------------------------
! print*,"*****************************************************"
! print*,"FLAG: TORSION"
! DO ISTEP=1,SIZE(TORSION_PARMS)
!    write(*,FMT='(4A4,F8.3,I4,F8.3)') TORSION_PARMS(ISTEP)
! END DO
! print*,"*****************************************************"
! STOP

      END SUBROUTINE FORCEFIELD$READ_PARMFILE

!     .................................................................. 
      SUBROUTINE FORCEFIELD$READ_TOPFILE
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        USE STRINGS_MODULE
        IMPLICIT NONE
!     ******************************************************************
        CHARACTER(8)                          :: ID='FF-TOP'
        integer                               :: nfil,istep
        CHARACTER(255)                        :: line, dummy
        LOGICAL(4)                            :: OK, TRESIDUE, TRESI
        INTEGER                               :: I, J, NMASS
        INTEGER                               :: NATOM_RES, NBOND_RES, NIMPROPER_RES, NIC_RES
        
        CHARACTER(4),        ALLOCATABLE      :: BONDS(:)
        INTEGER                               :: maxbonds
        CHARACTER(8)                          :: CRAP
!     ******************************************************************
                           CALL TRACE$PUSH('READ_TOPFILE')
        CALL FILEHANDLER$SETFILE(ID,.FALSE.,FF_TOP)
        CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
        CALL FILEHANDLER$UNIT(ID,NFIL)

        NRES=0
        NMASS=0
        NPATCH=0
!      --- FIRST LOOP TO GET THE NUMBER OF RESIDUES ---
        DO
           READ(NFIL,*) LINE
           IF(LINE(1:3).EQ.'RES')  NRES=NRES + 1
           IF(LINE(1:4).EQ.'MASS') NMASS= NMASS + 1
           IF(LINE(1:4).EQ.'PRES') NPATCH= NPATCH + 1
           IF(LINE(1:3).EQ.'END')  EXIT
        ENDDO

        REWIND(NFIL)
        ALLOCATE(TOP_RES(NRES+NPATCH))
        ALLOCATE(MASSES(NMASS))
        DO I=1,SIZE(TOP_RES)
           NULLIFY(TOP_RES(I)%ATOM)
           NULLIFY(TOP_RES(I)%BOND)
           NULLIFY(TOP_RES(I)%IMPROPER)
           NULLIFY(TOP_RES(I)%IC)
        ENDDO
print*,"SIZE OF MASSES: ",SIZE(MASSES)
print*,"SIZE OF RES: (including patches)",   SIZE(TOP_RES)
!      --- SECOND LOOP TO GET THE NUMBER OF ENTRIES IN EACH RESIDUE
        I=0
        DO WHILE(I.LT.SIZE(TOP_RES))
           READ(NFIL,FMT='(A255)') LINE
           IF(+LINE(1:3).EQ.'RES') THEN
              NATOM_RES=0
              NBOND_RES=0
              NIMPROPER_RES=0
              NIC_RES=0
              I = I + 1
              dummy= LINE(INDEX(LINE,' '):LEN_TRIM(LINE))
              READ(DUMMY,*) TOP_RES(I)%RESNAME, TOP_RES(I)%TOTALCHARGE
              DO
                 READ(NFIL,FMT='(A255)') LINE
                 IF(LINE(1:3).EQ.'END') EXIT
                 IF(LINE(1:4).EQ.'ATOM') NATOM_RES= NATOM_RES+1
                 IF(LINE(1:4).EQ.'BOND') NBOND_RES= NBOND_RES+1 !WARNING: NBOND_RES IS THE NUMBER OF LINES NOT OF BONDS!
                 IF(LINE(1:8).EQ.'IMPROPER') NIMPROPER_RES = NIMPROPER_RES+1 
                 IF(LINE(1:2).EQ.'IC') NIC_RES = NIC_RES+1
                 IF(LINE(1:5).EQ.'PATCH') THEN
                    ALLOCATE(TOP_RES(I)%ATOM(NATOM_RES))
                    ALLOCATE(TOP_RES(I)%BOND(NBOND_RES))
                    ALLOCATE(TOP_RES(I)%IMPROPER(NIMPROPER_RES))
                    ALLOCATE(TOP_RES(I)%IC(NIC_RES))
                    EXIT
                 END IF
              ENDDO
           END IF
           IF(LINE(1:3).EQ.'END') EXIT
        ENDDO
        REWIND(NFIL)
!      --- THIRD LOOP TO READ THE TOPOLOGIES OF THE RESIDUES
        I=0
        DO WHILE(I.LT.SIZE(TOP_RES))
!         --- WORKAROUND BECAUSE OF DIFFERENT FORMAT IN RESIDUE (NUCLEINACID) AND RES (AMINOACID) ENTRIES
           TRESI=.FALSE.
           TRESIDUE=.FALSE.
!         ----------
           READ(NFIL,FMT='(A255)') LINE

           IF(+LINE(1:5).EQ.'RESI ') TRESI=.TRUE.
           IF(+LINE(1:5).EQ.'RESID') TRESIDUE=.TRUE.

           IF(+LINE(1:3).EQ.'RES') THEN  !ENTER i-th RESIDUE
              I = I + 1
              ALLOCATE(BONDS(3* SIZE(TOP_RES(I)%BOND) * 5)) !AUXILLARY ARRAY TO SAVE THE BOND ENTRIES.
                                                            !CHECK IF THE ARRAY IS BIG ENOUGH!
              DO
                 READ(NFIL,FMT='(A255)') LINE
                 IF(LINE(1:4).EQ.'ATOM') THEN !READ ATOM DATA FROM LINE, WRITE IT TO TOP_RES AND READ NEXT LINE
                    DO J=1,SIZE(TOP_RES(I)%ATOM)
                       IF(TRESI) READ(LINE,FMT='(5X,A4,1X,A2,4X,F8.5)') TOP_RES(I)%ATOM(J)%NAME, &
                            & TOP_RES(I)%ATOM(J)%ATOM, TOP_RES(I)%ATOM(J)%CHARGE
                       IF(TRESIDUE) READ(LINE,FMT='(6X,A4,2X,A2,2X,F7.4)') TOP_RES(I)%ATOM(J)%NAME, &
                            & TOP_RES(I)%ATOM(J)%ATOM, TOP_RES(I)%ATOM(J)%CHARGE
                       READ(NFIL,FMT='(A255)') LINE
                    ENDDO
                 END IF

                 IF(LINE(1:4).EQ.'BOND') THEN !READ THE BONDS. STORE THEM IN BONDS-ARRAY. COPY THEM TO TOP_RES LATER
                    MAXBONDS=0
                    DO J=1,SIZE(TOP_RES(I)%BOND)
                       DUMMY = ADJUSTL(LINE(5:LEN_TRIM(LINE)))
                       IF(SCAN(DUMMY,'!').NE.0) DUMMY = DUMMY(1:SCAN(DUMMY,'!')-1)
                       DO 
                          MAXBONDS= MAXBONDS + 1
                          BONDS(MAXBONDS) = DUMMY(1:SCAN(DUMMY,' ')-1)
                          DUMMY=ADJUSTL(DUMMY(SCAN(DUMMY,' '):LEN_TRIM(DUMMY)))
                          IF(LEN_TRIM(ADJUSTL(DUMMY)).EQ.0) EXIT
                       ENDDO
                       READ(NFIL,FMT='(A255)') LINE
                    ENDDO
                    IF(MOD(MAXBONDS,2).NE.0) THEN
                       CALL ERROR$MSG('ODD NUMBER OF ATOM ENTRIES IN BONDS')
                       CALL ERROR$CHVAL('RESIDUE',TOP_RES(I)%RESNAME)
                       CALL ERROR$STOP('FORCEFIELD%READ_TOPFILE')
                    END IF
                 END IF
                 
                 IF(LINE(1:8).EQ.'IMPROPER') THEN !READ IMPROPER DATA FROM LINE
                    DO J=1,SIZE(TOP_RES(I)%IMPROPER)
                       READ(LINE,*) CRAP, TOP_RES(I)%IMPROPER(J)%ATOM1, & !FMT='(11X,A4,X,A4,X,A4,X,A4)'
                      &     TOP_RES(I)%IMPROPER(J)%ATOM2, TOP_RES(I)%IMPROPER(J)%ATOM3, &
                      &     TOP_RES(I)%IMPROPER(J)%ATOM4
                       READ(NFIL,FMT='(A255)') LINE
                    ENDDO
                 END IF

                 IF(LINE(1:2).EQ.'IC') THEN !READ IC DATA FROM LINE
                    DO J=1,SIZE(TOP_RES(I)%IC)
                       READ(LINE,FMT='(3X,A4,1X,A4,1X,A4,1X,A4,3X,F7.4,1X,F7.2,1X,F7.2,1X,F7.2,1X,F7.4)') &
                      &     TOP_RES(I)%IC(J)%ATOM1, TOP_RES(I)%IC(J)%ATOM2, TOP_RES(I)%IC(J)%ATOM3, &
                      &     TOP_RES(I)%IC(J)%ATOM4, TOP_RES(I)%IC(J)%DIST12, TOP_RES(I)%IC(J)%ANG123, &
                      &     TOP_RES(I)%IC(J)%TOR1234, TOP_RES(I)%IC(J)%ANG234, TOP_RES(I)%IC(J)%DIST34
                       READ(NFIL,FMT='(A255)') LINE
                    ENDDO
                 END IF
                 IF(LINE(1:5).EQ.'PATCH') THEN
!               ---- COPY BONDS TO NEW ALLOCATED TOP_RES(I)%BOND
                    DEALLOCATE(TOP_RES(I)%BOND)
                    ALLOCATE(TOP_RES(I)%BOND(MAXBONDS / 2))
                    DO J=2, MAXBONDS, 2
                       TOP_RES(I)%BOND(J/2)%ATOM1 = BONDS(J-1)
                       TOP_RES(I)%BOND(J/2)%ATOM2 = BONDS(J)
                    ENDDO
                    MAXBONDS=0
                    EXIT !LEAVE THE RESIDUE LOOP AND LOOK FOR NEW RESIDUE
                 END IF
              END DO
              DEALLOCATE(BONDS)
           END IF
           IF(LINE(1:3).EQ.'END') EXIT
        END DO
        REWIND(NFIL)
        
        CALL FORCEFIELD_READPATCHRESIDUE

! !**** PRINTOUT FOR DEBUG
! print*,"FLAG: ********** TOP_RES ***************"
! DO I=1,SIZE(TOP_RES)
!    print*,"-------------------------------------------------------"
!    print*, TOP_RES(I)%RESNAME, TOP_RES(I)%TOTALCHARGE, I

!    DO J=1, SIZE(TOP_RES(I)%ATOM)
!       WRITE(*,FMT='(2A6,F10.5)') TOP_RES(I)%ATOM(J)%NAME, TOP_RES(I)%ATOM(J)%ATOM, TOP_RES(I)%ATOM(J)%CHARGE
!    END DO
!     print*,"SIZE OF BOND ARRAY: ",SIZE(TOP_RES(I)%BOND)
!     DO J=1,SIZE(TOP_RES(I)%BOND)
!        WRITE(*,FMT='(2A5)') TOP_RES(I)%BOND(J)%ATOM1, TOP_RES(I)%BOND(J)%ATOM2
!     ENDDO

!     DO J=1,SIZE(TOP_RES(I)%IMPROPER)
!        WRITE(*,FMT='(4A5)') TOP_RES(I)%IMPROPER(J)%ATOM1, & 
!             &     TOP_RES(I)%IMPROPER(J)%ATOM2, TOP_RES(I)%IMPROPER(J)%ATOM3, &
!             &     TOP_RES(I)%IMPROPER(J)%ATOM4
!     ENDDO

!     DO J=1,SIZE(TOP_RES(I)%IC)
!        WRITE(*,FMT='(4A5,5F7.3)') TOP_RES(I)%IC(J)%ATOM1, TOP_RES(I)%IC(J)%ATOM2, TOP_RES(I)%IC(J)%ATOM3, &
!             & TOP_RES(I)%IC(J)%ATOM4, TOP_RES(I)%IC(J)%dist12, TOP_RES(I)%IC(J)%ang123, TOP_RES(I)%IC(J)%tor1234, &
!             & TOP_RES(I)%IC(J)%ang234, TOP_RES(I)%IC(J)%dist34
!     ENDDO
!    PRINT*,"-------------------------------------------------------"
!    PRINT*
! ENDDO
! print*,"SIZE OF TOP_RES: ",SIZE(TOP_RES)
! STOP
                          CALL TRACE$POP
      END SUBROUTINE FORCEFIELD$READ_TOPFILE


!     ..................................................................
      SUBROUTINE FORCEFIELD_READPATCHRESIDUE
!     ******************************************************************
!     **  COPIES THE VARIABLE TOP_RES(ISOURCE) OF TYPE RES_TYPE       **
!     **  INTO TOP_RES(ITARGET). THE TARGET WILL BE DEALLOCATED AND   **
!     **  ALLOCATED AGAIN WITH THE RIGHT NUMBER OF ENTRIES            **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        USE STRINGS_MODULE
        IMPLICIT NONE
        INTEGER                      :: NFIL, I, IRES
        CHARACTER(255)               :: line, dummy
        CHARACTER(4)                 :: NAME
        INTEGER                      :: NPATCHRES
        INTEGER                      :: J
        TYPE(TOP_ATOM_TYPE)          :: DUMMY_ATOM
        TYPE(TOP_BOND_TYPE)          :: DUMMY_BOND
        TYPE(TOP_IMPROPER_TYPE)      :: DUMMY_IMPROPER
        TYPE(TOP_IC_TYPE)            :: DUMMY_IC
        LOGICAL                      :: TCHK
!     ******************************************************************
                              CALL TRACE$PUSH('READPATCHRESIDUE')
        CALL FILEHANDLER$UNIT('FF-TOP',NFIL)
        REWIND(NFIL)
        I= 0
        READ(NFIL,FMT='(A255)') LINE
        DO WHILE(I.LT.NPATCH)
            IF(+LINE(1:4).EQ.'PRES') THEN  !ENTER i-th PATCHRESIDUE
              I = I + 1
              NAME= LINE(7:9)         !FOR AMINO ACIDS (FORMAT: NXXX OR CXXX)
              IF(SCAN(LINE(6:9),'0123456789').NE.0) NAME= LINE(6:8) !FOR NUCLEIC ACID (FORMAT: XXX5 OR XXX3)
              IF(+LINE(6:7).EQ.'DO') THEN !FOR DEOXYRIBOSES
                 IF(LINE(8:8).EQ.'A') NAME= 'ADE'
                 IF(LINE(8:8).EQ.'C') NAME= 'CYT'
                 IF(LINE(8:8).EQ.'G') NAME= 'GUA'
              END IF
              CALL FORCEFIELD$GETRESNUMBER(TRIM(ADJUSTL(NAME)),IRES)
              NPATCHRES= NRES+I
              CALL FORCEFIELD_COPYRES(IRES,NPATCHRES)
              READ(LINE,*) DUMMY, TOP_RES(NPATCHRES)%RESNAME, TOP_RES(NPATCHRES)%TOTALCHARGE
              DO
                 READ(NFIL,FMT='(A255)') LINE
                 IF(LINE(1:3).EQ.'END') EXIT
                 IF(LINE(1:4).EQ.'PRES') EXIT
!             ---- DELETE ATOMS
                 IF(LINE(1:6).EQ.'DELETE') THEN
                    READ(LINE,*) DUMMY, DUMMY, NAME
                    CALL FORCEFIELD$DELETEATOM(NPATCHRES,NAME)
                    CYCLE
                 END IF
!             ---- CHANGE VALUES: ATOMS
                 IF(LINE(1:4).EQ.'ATOM') THEN
                    TCHK=.FALSE.
                    READ(LINE,*) DUMMY, DUMMY_ATOM%NAME, DUMMY_ATOM%ATOM, DUMMY_ATOM%CHARGE
                    DO J=1, SIZE(TOP_RES(NPATCHRES)%ATOM)
                       IF(TOP_RES(NPATCHRES)%ATOM(J)%NAME.EQ.DUMMY_ATOM%NAME) THEN
                          TCHK=.TRUE.
                          TOP_RES(NPATCHRES)%ATOM(J)%ATOM = DUMMY_ATOM%ATOM
                          TOP_RES(NPATCHRES)%ATOM(J)%CHARGE = DUMMY_ATOM%CHARGE
                       END IF
                    ENDDO
                    IF(.NOT.TCHK) THEN
                       CALL FORCEFIELD$ADDATOM(NPATCHRES,DUMMY_ATOM)
                    END IF
                    CYCLE
                 END IF
!             ---- CHANGE VALUES: BONDS
                 IF(LINE(1:4).EQ.'BOND') THEN
                    READ(LINE,*) DUMMY, DUMMY_BOND%ATOM1, DUMMY_BOND%ATOM2
                    CALL FORCEFIELD$ADDBOND(NPATCHRES,DUMMY_BOND)
                    CYCLE
                 END IF
!             ---- CHANGE VALUES: IC
                 IF(LINE(1:2).EQ.'IC') THEN
                    READ(LINE,*) DUMMY, DUMMY_IC%ATOM1, DUMMY_IC%ATOM2, DUMMY_IC%ATOM3, DUMMY_IC%ATOM4,&
                         &  DUMMY_IC%dist12, DUMMY_IC%ang123, DUMMY_IC%tor1234, DUMMY_IC%ang234, DUMMY_IC%dist34
                    CALL FORCEFIELD$ADDIC(NPATCHRES,DUMMY_IC)
                    CYCLE
                 END IF                       
              END DO
           ELSE
              READ(NFIL,FMT='(A255)') LINE
           END IF
           IF(LINE(1:3).EQ.'END') EXIT
        END DO
                      CALL TRACE$POP
      END SUBROUTINE FORCEFIELD_READPATCHRESIDUE


!     ..................................................................
      SUBROUTINE FORCEFIELD_COPYRES(ISOURCE,ITARGET)
!     ******************************************************************
!     **  COPIES THE VARIABLE TOP_RES(ISOURCE) OF TYPE RES_TYPE INTO  **
!     **  TOP_RES(ITARGET). THE TARGET WILL BE DEALLOCATED AND        **
!     **  ALLOCATED AGAIN WITH THE RIGHT NUMBER OF ENTRIES            **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        IMPLICIT NONE
        INTEGER,           INTENT(IN):: ISOURCE
        INTEGER,           INTENT(IN):: ITARGET
        INTEGER                      :: I
!     ******************************************************************
                  CALL TRACE$PUSH('FORCEFIELD_COPYRES')
        IF(ASSOCIATED(TOP_RES(ITARGET)%ATOM)) DEALLOCATE(TOP_RES(ITARGET)%ATOM)
        ALLOCATE(TOP_RES(ITARGET)%ATOM(SIZE(TOP_RES(ISOURCE)%ATOM)))
        IF(ASSOCIATED(TOP_RES(ITARGET)%BOND)) DEALLOCATE(TOP_RES(ITARGET)%BOND)
        ALLOCATE(TOP_RES(ITARGET)%BOND(SIZE(TOP_RES(ISOURCE)%BOND)))
        IF(ASSOCIATED(TOP_RES(ITARGET)%IMPROPER)) DEALLOCATE(TOP_RES(ITARGET)%IMPROPER)
        ALLOCATE(TOP_RES(ITARGET)%IMPROPER(SIZE(TOP_RES(ISOURCE)%IMPROPER)))
        IF(ASSOCIATED(TOP_RES(ITARGET)%IC)) DEALLOCATE(TOP_RES(ITARGET)%IC)
        ALLOCATE(TOP_RES(ITARGET)%IC(SIZE(TOP_RES(ISOURCE)%IC)))

        DO I=1, SIZE(TOP_RES(ISOURCE)%ATOM)
           TOP_RES(ITARGET)%ATOM(I) = TOP_RES(ISOURCE)%ATOM(I)
        ENDDO
        DO I=1, SIZE(TOP_RES(ISOURCE)%BOND)
           TOP_RES(ITARGET)%BOND(I) = TOP_RES(ISOURCE)%BOND(I)
        ENDDO
        DO I=1, SIZE(TOP_RES(ISOURCE)%IMPROPER)
           TOP_RES(ITARGET)%IMPROPER(I) = TOP_RES(ISOURCE)%IMPROPER(I)
        ENDDO
        DO I=1, SIZE(TOP_RES(ISOURCE)%IC)
           TOP_RES(ITARGET)%IC(I) = TOP_RES(ISOURCE)%IC(I)
        ENDDO
                  CALL TRACE$POP

      END SUBROUTINE FORCEFIELD_COPYRES

!     ..................................................................
      SUBROUTINE FORCEFIELD$ADDATOM(RESIDUE,NEWATOM)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        IMPLICIT NONE
        INTEGER,             INTENT(IN)       :: RESIDUE
        TYPE(TOP_ATOM_TYPE), INTENT(IN)       :: NEWATOM
        TYPE(TOP_ATOM_TYPE), ALLOCATABLE      :: OLDATOMS(:)
        INTEGER                               :: IVAR
!     ******************************************************************
                  CALL TRACE$PUSH('FORCEFIELD$ADDATOM')
        ALLOCATE(OLDATOMS(SIZE(TOP_RES(RESIDUE)%ATOM)))
        DO IVAR=1, SIZE(TOP_RES(RESIDUE)%ATOM)
           OLDATOMS(IVAR) = TOP_RES(RESIDUE)%ATOM(IVAR)
        ENDDO
        DEALLOCATE(TOP_RES(RESIDUE)%ATOM)
        ALLOCATE(TOP_RES(RESIDUE)%ATOM(SIZE(OLDATOMS)+1))
        DO IVAR=1, SIZE(OLDATOMS)
           TOP_RES(RESIDUE)%ATOM(IVAR) = OLDATOMS(IVAR)
        ENDDO
        TOP_RES(RESIDUE)%ATOM(SIZE(OLDATOMS)+1) = NEWATOM
        DEALLOCATE(OLDATOMS)
                  CALL TRACE$POP
        RETURN
      END SUBROUTINE FORCEFIELD$ADDATOM

!     ..................................................................
      SUBROUTINE FORCEFIELD$ADDBOND(RESIDUE,NEWBOND)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        IMPLICIT NONE
        INTEGER,             INTENT(IN)       :: RESIDUE
        TYPE(TOP_BOND_TYPE), INTENT(IN)       :: NEWBOND
        TYPE(TOP_BOND_TYPE), ALLOCATABLE      :: OLDBONDS(:)
        INTEGER                               :: IVAR
!     ******************************************************************
                  CALL TRACE$PUSH('FORCEFIELD$ADDBOND')
        ALLOCATE(OLDBONDS(SIZE(TOP_RES(RESIDUE)%BOND)))
        DO IVAR=1, SIZE(TOP_RES(RESIDUE)%BOND)
           OLDBONDS(IVAR) = TOP_RES(RESIDUE)%BOND(IVAR)
        ENDDO
        DEALLOCATE(TOP_RES(RESIDUE)%BOND)
        ALLOCATE(TOP_RES(RESIDUE)%BOND(SIZE(OLDBONDS)+1))
        DO IVAR=1, SIZE(OLDBONDS)
           TOP_RES(RESIDUE)%BOND(IVAR) = OLDBONDS(IVAR)
        ENDDO
        TOP_RES(RESIDUE)%BOND(SIZE(OLDBONDS)+1) = NEWBOND
        DEALLOCATE(OLDBONDS)
                  CALL TRACE$POP
        RETURN
      END SUBROUTINE FORCEFIELD$ADDBOND

!     ..................................................................
      SUBROUTINE FORCEFIELD$ADDIMPROPER(RESIDUE,NEWIMPROPER)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        IMPLICIT NONE
        INTEGER,             INTENT(IN)       :: RESIDUE
        TYPE(TOP_IMPROPER_TYPE), INTENT(IN)       :: NEWIMPROPER
        TYPE(TOP_IMPROPER_TYPE), ALLOCATABLE      :: OLDIMPROPERS(:)
        INTEGER                               :: IVAR
!     ******************************************************************
                  CALL TRACE$PUSH('FORCEFIELD$ADDIMPROPER')
        ALLOCATE(OLDIMPROPERS(SIZE(TOP_RES(RESIDUE)%IMPROPER)))
        DO IVAR=1, SIZE(TOP_RES(RESIDUE)%IMPROPER)
           OLDIMPROPERS(IVAR) = TOP_RES(RESIDUE)%IMPROPER(IVAR)
        ENDDO
        DEALLOCATE(TOP_RES(RESIDUE)%IMPROPER)
        ALLOCATE(TOP_RES(RESIDUE)%IMPROPER(SIZE(OLDIMPROPERS)+1))
        DO IVAR=1, SIZE(OLDIMPROPERS)
           TOP_RES(RESIDUE)%IMPROPER(IVAR) = OLDIMPROPERS(IVAR)
        ENDDO
        TOP_RES(RESIDUE)%IMPROPER(SIZE(OLDIMPROPERS)+1) = NEWIMPROPER
        DEALLOCATE(OLDIMPROPERS)
                  CALL TRACE$POP
        RETURN
      END SUBROUTINE FORCEFIELD$ADDIMPROPER

!     ..................................................................
      SUBROUTINE FORCEFIELD$ADDIC(RESIDUE,NEWIC)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        IMPLICIT NONE
        INTEGER,             INTENT(IN)       :: RESIDUE
        TYPE(TOP_IC_TYPE), INTENT(IN)       :: NEWIC
        TYPE(TOP_IC_TYPE), ALLOCATABLE      :: OLDICS(:)
        INTEGER                               :: IVAR
!     ******************************************************************
                  CALL TRACE$PUSH('FORCEFIELD$ADDIC')
        ALLOCATE(OLDICS(SIZE(TOP_RES(RESIDUE)%IC)))
        DO IVAR=1, SIZE(TOP_RES(RESIDUE)%IC)
           OLDICS(IVAR) = TOP_RES(RESIDUE)%IC(IVAR)
        ENDDO
        DEALLOCATE(TOP_RES(RESIDUE)%IC)
        ALLOCATE(TOP_RES(RESIDUE)%IC(SIZE(OLDICS)+1))
        DO IVAR=1, SIZE(OLDICS)
           TOP_RES(RESIDUE)%IC(IVAR) = OLDICS(IVAR)
        ENDDO
        TOP_RES(RESIDUE)%IC(SIZE(OLDICS)+1) = NEWIC
        DEALLOCATE(OLDICS)
                  CALL TRACE$POP
        RETURN
      END SUBROUTINE FORCEFIELD$ADDIC

!     ..................................................................
      SUBROUTINE FORCEFIELD$DELETEATOM(RESIDUE,NAME)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        IMPLICIT NONE
        INTEGER,             INTENT(IN)       :: RESIDUE
        CHARACTER(4),        INTENT(IN)       :: NAME
        TYPE(TOP_ATOM_TYPE), ALLOCATABLE      :: OLDATOMS(:)
        INTEGER                               :: IVAR,I
!     ******************************************************************
                  CALL TRACE$PUSH('FORCEFIELD$DELETEATOM')
        ALLOCATE(OLDATOMS(SIZE(TOP_RES(RESIDUE)%ATOM)))
        DO IVAR=1,SIZE(OLDATOMS)
           OLDATOMS(IVAR) = TOP_RES(RESIDUE)%ATOM(IVAR)
        ENDDO
        DEALLOCATE(TOP_RES(RESIDUE)%ATOM)
        ALLOCATE(TOP_RES(RESIDUE)%ATOM(SIZE(OLDATOMS)-1))
        I= 0
        DO IVAR = 1, SIZE(OLDATOMS)
           IF(TRIM(ADJUSTL(OLDATOMS(IVAR)%NAME)).NE.TRIM(ADJUSTL(NAME))) THEN
              I = I + 1
              TOP_RES(RESIDUE)%ATOM(I) = OLDATOMS(IVAR)
           END IF
        ENDDO
        DEALLOCATE(OLDATOMS)
                  CALL TRACE$POP
        RETURN
      END SUBROUTINE FORCEFIELD$DELETEATOM

!     ..................................................................
      SUBROUTINE FORCEFIELD$READ_MMSTRC
!     ******************************************************************
!     **  READS THE MM-SYSTEM FROM A PDB FILE.                        **
!     **  WARNING: AT THIS STAGE NO CONECT ENTRIES ARE TREATED        **
!     ******************************************************************
        USE FORCEFIELD_MODULE
        USE STRINGS_MODULE
        IMPLICIT NONE
        INTEGER                               :: NFIL, NATOM, NHETATM, NCONECT
        LOGICAL                               :: TCHK
        CHARACTER(102)                        :: LINE, DUMMY
        CHARACTER(12)                         :: ID
        CHARACTER(LEN=*), PARAMETER           :: con_form='(A6,5I5)'  !FORMAT STRING TO READ THE CONECT ENTRIES IN THE PDB FILE
        INTEGER                               :: I, IATOM, IHETATM, ICONECT
                               CALL TRACE$PUSH('READ_MMSTRC')
!     ==  MM STRUCTURE FILE  ===========================================
        CALL FILEHANDLER$UNIT('MMSTRC',NFIL)
        REWIND(NFIL)
        TCHK=.TRUE.
        NATOM=0
        NHETATM=0
        NCONECT=0
        DO
           READ(NFIL,FMT='(A102)',END=201) LINE
           IF(LINE(1:6).EQ.'ATOM  ') NATOM = NATOM + 1
           IF(LINE(1:6).EQ.'HETATM') NHETATM = NHETATM + 1
           IF(LINE(1:6).EQ.'CONECT') NCONECT = NCONECT + 1
        END DO
201     REWIND(NFIL)
        IF(ALLOCATED(MMATOM)) THEN
           CALL ERROR$MSG('MMATOM ARRAY ALREADY ALLOCATED.')
           CALL ERROR$STOP('FORCEFIELD$READ_MMSTRC')
        END IF
        ALLOCATE(MMATOM(NATOM+NHETATM))
        ALLOCATE(MMCONECT(NCONECT))
        MMCONECT(:)%CATOM=0
        MMCONECT(:)%ATOM1=0
        MMCONECT(:)%ATOM2=0
        MMCONECT(:)%ATOM3=0
        MMCONECT(:)%ATOM4=0
        IATOM=0
        IHETATM=0
        ICONECT=0
        DO 
           READ(NFIL,FMT='(A102)',END=202) LINE
           IF(LINE(1:6).EQ.'ATOM  ') THEN 
              IATOM = IATOM + 1
              READ(LINE,pdb_form) MMATOM(IATOM)
           END IF
           IF(LINE(1:6).EQ.'HETATM') THEN
              IHETATM = IHETATM + 1
              READ(LINE,pdb_form) MMATOM(NATOM+IHETATM)
           END IF
           IF(LINE(1:6).EQ.'CONECT') THEN
              ICONECT = ICONECT + 1
              READ(LINE,CON_FORM) MMCONECT(ICONECT)
           END IF
        END DO
! print*,"--------DEBUG--------"
! DO I=1,SIZE(MMCONECT)
!    write(*,con_form) MMCONECT(I)
! END DO
! print*,"---------------------"
! STOP

202     DO I=1, SIZE(MMATOM)
!     ----- if you use 'reduce' to fix pdb-files then you get atomnames like 3HG1 instead of HG13. this loop shall fix this
           DUMMY = MMATOM(I)%NAME
           IF(SCAN(DUMMY,'1234567890').EQ.1) THEN
              MMATOM(I)%NAME = ' '//TRIM(DUMMY(2:LEN(DUMMY))) // DUMMY(1:1)
           END if
!     ----- IF ELEMENT COLUMN IS EMPTY COPY ELEMENT NAME FROM ATOM NAME
           IF(MMATOM(I)%ELEMENT.EQ.' ') THEN
              MMATOM(I)%ELEMENT = ADJUSTR(MMATOM(I)%NAME(1:2))
           END IF
!     ----- If there is only one chain in the protein, some PDB file don't specify this. Then we set the chain name to 'A'
           IF(MMATOM(I)%KEYWORD.EQ.'ATOM  '.AND.MMATOM(I)%CHAINID.EQ.' ') THEN
              MMATOM(I)%CHAINID='A'
           END IF
        ENd DO
!     --------------------------              
                     CALL TRACE$POP
        RETURN
      END SUBROUTINE FORCEFIELD$READ_MMSTRC

!     ..................................................................
      SUBROUTINE FORCEFIELD$WRITE_MMSTRC
!     ******************************************************************
!     **  WRITES THE MM-SYSTEM INTO A PDB FILE.                       **
!     **  WARNING: AT THIS STAGE NO CONECT ENTRIES ARE TREATED        **
!     ******************************************************************
        USE LINKEDLIST_MODULE
        USE FORCEFIELD_MODULE
        USE STRINGS_MODULE
        IMPLICIT NONE
        INTEGER                               :: NFIL, NATOM
        INTEGER                               :: IATOM
        REAL(8), ALLOCATABLE                  :: R(:,:)
        REAL(8)                               :: UNIT
        LOGICAL(4)                            :: TCHK
        CHARACTER(LEN=*), PARAMETER           :: con_form='(A6,5I5)'  !FORMAT STRING TO WRITE THE CONECT ENTRIES IN THE PDB FILE

                     call trace$push('FORCEFIELD$WRITE_MMSTRC')
        CALL CONSTANTS$GET('ANGSTROM',UNIT)
                     
        CALL FILEHANDLER$SETFILE('MMSTRC_OUT',.TRUE.,-'.PDB_OUT')
        CALL FILEHANDLER$SETSPECIFICATION('MMSTRC_OUT','STATUS','REPLACE')
        CALL FILEHANDLER$SETSPECIFICATION('MMSTRC_OUT','POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION('MMSTRC_OUT','ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION('MMSTRC_OUT','FORM','FORMATTED')

        CALL FILEHANDLER$UNIT('MMSTRC_OUT',NFIL)
        NATOM=SIZE(MMATOM)
        ALLOCATE(R(3,NATOM))
        CALL CLASSICAL$SELECT('QMMM')
        CALL CLASSICAL$GETR8A('R(0)',3*NATOM,R)
        DO IATOM=1,NATOM 
           MMATOM(IATOM)%R=R(:,IATOM) / UNIT
           WRITE(NFIL,pdb_form) MMATOM(IATOM)
        END DO
        IF(SIZE(MMCONECT).GT.0) THEN
           DO IATOM=1,SIZE(MMCONECT)
              WRITE(NFIL,CON_FORM) MMCONECT(IATOM)
           END DO
        END IF
        DEALLOCATE(R)
                     call trace$POP
        RETURN
      END SUBROUTINE FORCEFIELD$WRITE_MMSTRC
!     ..................................................................

!*************************************************************************************
!******   AMBER SPECIFIC ROUTINES                *************************************
!*************************************************************************************

      SUBROUTINE FORCEFIELD$AMBER_BONDPARMS(TYPE1,TYPE2,ID,R,K,D,TCHK)
        USE FORCEFIELD_MODULE, ONLY: BOND_PARMS
        IMPLICIT NONE
        CHARACTER(5),       INTENT(IN) :: TYPE1      ! FF-ATOM-TYPE
        CHARACTER(5),       INTENT(IN) :: TYPE2      ! FF-ATOM-TYPE
        CHARACTER(64),      INTENT(OUT):: ID         ! IDENTIFIER
        REAL(8),            INTENT(OUT):: R          ! BOND DISTANCE
        REAL(8),            INTENT(OUT):: K          ! FORCE CONSTANT
        REAL(8),            INTENT(OUT):: D          ! BINDING ENERGY
        LOGICAL(4),         INTENT(OUT):: TCHK       ! NON-ZERO FORCE FIELD
        INTEGER                        :: I
        REAL(8)                        :: ANGSTROM
        REAL(8)                        :: KCALBYMOL
!                        cALL TRACE$PUSH('AMBER_BONDPARMS')
!        PRINT*,"TYPES: ",TYPE1, TYPE2
        D=70.D0                   ! BOND DISSOCIATION ENERGY (ONLY FOR NON-HARMONIC)
        TCHK=.FALSE.
        DO I=1, SIZE(BOND_PARMS)
           IF((BOND_PARMS(I)%ATOM1.EQ.TYPE1.AND.BOND_PARMS(I)%ATOM2.EQ.TYPE2).OR.&
                & (BOND_PARMS(I)%ATOM1.EQ.TYPE2.AND.BOND_PARMS(I)%ATOM2.EQ.TYPE1)) THEN
              K= BOND_PARMS(I)%F_CONST
              R= BOND_PARMS(I)%R_EQUI
              TCHK=.TRUE.
              EXIT
           END IF
        END DO
        IF(.NOT.TCHK) THEN
           CALL ERROR$MSG('BOND PAIR NOT FOUND!')
           CALL ERROR$CHVAL('TYPE1: ',TYPE1)
           CALL ERROR$CHVAL('TYPE2: ',TYPE2)
           CALL ERROR$STOP('FORCEFIELD$AMBER_BONDPARMS')
        END IF
!     ==================================================================
!     == WRITE POTENTIAL ID                                          ==
!     ==================================================================
        IF(LGT(TYPE1,TYPE2)) THEN
           WRITE(ID,FMT='(A1,2A5," R=",F5.2," K=",F5.0," D=",F5.0)') &
                &        'B ',TYPE1,TYPE2,R,K,D
        ELSE
           WRITE(ID,FMT='(A1,2A5," R=",F5.2," K=",F5.0," D=",F5.0)') &
                &        'B ',TYPE2,TYPE1,R,K,D
        END IF
        CALL CONSTANTS('KCAL/MOL',KCALBYMOL)
        CALL CONSTANTS('ANGSTROM',ANGSTROM)
        R=R*ANGSTROM
        D=D*KCALBYMOL
        K=K*KCALBYMOL/ANGSTROM**2
        TCHK=(K.GT.1.D-6)
        RETURN
!                      CALL TRACE$POP
      END SUBROUTINE FORCEFIELD$AMBER_BONDPARMS

!     ******************************************************************
      SUBROUTINE FORCEFIELD$AMBER_ANGLEPARMS(ATOM1,ATOM2,ATOM3,ID,THETA,K,TCHK)
        USE FORCEFIELD_MODULE, ONLY: ANGLE_PARMS
        USE CLASSICAL_MODULE,  ONLY: POT_TYPE
        IMPLICIT NONE
        CHARACTER(5)  ,INTENT(IN)  :: ATOM1 ! FF-ATOM-TYPE TERMINAL ATOM   
        CHARACTER(5)  ,INTENT(IN)  :: ATOM2 ! FF-ATOM-TYPE CENTRAL ATOM   
        CHARACTER(5)  ,INTENT(IN)  :: ATOM3 ! FF-ATOM-TYPE TERMINAL ATOM
        CHARACTER(64) ,INTENT(OUT) :: ID    ! IDENTIFIER
        REAL(8)       ,INTENT(OUT) :: THETA ! EQUILIBRIUM ANGLE
        REAL(8)       ,INTENT(OUT) :: K     ! FORCE CONSTANT
        LOGICAL(4)    ,INTENT(OUT) :: TCHK
        INTEGER                    :: I
        REAL(8)                    :: PI,KCALBYMOL
        REAL(8)                    :: ZI,ZK
        REAL(8)                    :: RIJ,RJK,RIK,BETA
        LOGICAL(4)                 :: TCHK1
!                         CALL TRACE$PUSH('AMBER_ANGLEPARMS')
!     ******************************************************************
        PI=4.D0*DATAN(1.D0)
!     ==================================================================
!     ==  APPLY GENERAL RULE                                          ==
!     ==================================================================
        TCHK=.FALSE.
        DO I=1, SIZE(ANGLE_PARMS)
           IF(ANGLE_PARMS(I)%ATOM2.EQ.ATOM2) THEN
              IF((ANGLE_PARMS(I)%ATOM1.EQ.ATOM1.AND.ANGLE_PARMS(I)%ATOM3.EQ.ATOM3).OR.&
                   (ANGLE_PARMS(I)%ATOM1.EQ.ATOM3.AND.ANGLE_PARMS(I)%ATOM3.EQ.ATOM1)) THEN
                 K= ANGLE_PARMS(I)%F_CONST
                 THETA= ANGLE_PARMS(I)%THETA_EQUI
                 THETA= THETA / 180.d0 * PI
                 EXIT
              END IF
           END IF
        END DO
!
!     ==================================================================
!     ==  WRITE STRING IDENTIFYING POTENTIAL                          ==
!     ==================================================================
        IF(LGT(ATOM1,ATOM3)) THEN
           WRITE(ID,FMT='(A1,3A5," THETA=",F5.0," K=",F5.0)') &
                &                'A ',ATOM1,ATOM2,ATOM3,THETA/PI*180,K
        ELSE
           WRITE(ID,FMT='(A1,3A5," THETA=",F5.0," K=",F5.0)') &
                &                'A ',ATOM3,ATOM2,ATOM1,THETA/PI*180,K
        END IF
!
!     ==================================================================
!     == CONVERT TO ATOMIC UNITS                                      ==
!     ==================================================================
        CALL CONSTANTS('KCAL/MOL',KCALBYMOL)
        K=K*KCALBYMOL
!
!     ==================================================================
!     ==  TEST FOR ZERO POTENTIALS                                    ==
!     ==================================================================
        TCHK=(K.GT.1.D-6)
!                      CALL TRACE$POP
        RETURN
      END SUBROUTINE FORCEFIELD$AMBER_ANGLEPARMS



      SUBROUTINE FORCEFIELD$AMBER_ANGLEPOTA(ID,THETA0,KIJK,POT)
!     ******************************************************************
!     **  AMBER USES E_angle= K_theta * (THETA - THETA0) ** 2         **
!     ******************************************************************
        USE CLASSICAL_MODULE, ONLY: POT_TYPE
        IMPLICIT NONE
        CHARACTER(*)  ,INTENT(IN)  :: ID
        REAL(8)       ,INTENT(IN)  :: THETA0   ! OPTIMUM ANGLE
        REAL(8)       ,INTENT(IN)  :: KIJK     !FORCE CONSTANT
        TYPE(POT_TYPE),INTENT(OUT) :: POT
        INTEGER(4)    ,PARAMETER   :: NX=1001
        REAL(8)                    :: X1
        REAL(8)                    :: X2
        INTEGER(4)                 :: I
        REAL(8)                    :: THETA
        REAL(8)                    :: PI
        REAL(8)                    :: X

integer   :: nfilinfo

!                      CALL TRACE$PUSH('AMBER_ANGLEPOTA')
!     ******************************************************************
        PI=4.D0*DATAN(1.D0)
        X1 =  -PI
        X2 =  PI
!     ==================================================================
!     ==  DETERMINE ARRAY SIZE ETC                                    ==
!     ==================================================================
        POT%ID=ID
        POT%NX=NX
        POT%X1=X1
        POT%DX=(X2-X1)/DBLE(POT%NX-1)
        ALLOCATE(POT%VAL(POT%NX))
        ALLOCATE(POT%DER(POT%NX))
        DO I=1,POT%NX
           THETA=POT%X1+POT%DX*DBLE(I-1) 
           POT%VAL(I)=KIJK*(THETA-THETA0)**2
           POT%DER(I)=2*KIJK*(THETA-THETA0)
        END DO
! !---- DEBUG remove this later
!         CALL FILEHANDLER$UNIT('INFO',nfilinfo)
!         DO I=1,POT%NX
!            WRITE(nfilinfo,FMT='(I8,2F15.7)') I, POT%VAL(I), POT%DER(I)
!         END DO
!         print*,"FLAG FORCED STOP ANGLEPOT INFO"
!         STOP
! !-----------------------
!                      CALL TRACE$POP
        RETURN 
      END SUBROUTINE FORCEFIELD$AMBER_ANGLEPOTA
!

!     ......................................................................
      SUBROUTINE FORCEFIELD$AMBER_TORSIONPARMS(TYPE1,TYPE2,TYPE3,TYPE4, &
           & ID, PHI0, NJK, VBARRIER, TCHK)
!     **********************************************************************
        USE FORCEFIELD_MODULE
        USE CLASSICAL_MODULE, ONLY: POT_TYPE
        IMPLICIT NONE
        CHARACTER(5)  ,INTENT(IN) :: TYPE1
        CHARACTER(5)  ,INTENT(IN) :: TYPE2
        CHARACTER(5)  ,INTENT(IN) :: TYPE3
        CHARACTER(5)  ,INTENT(IN) :: TYPE4
        CHARACTER(64) ,INTENT(OUT):: ID
        REAL(8)       ,INTENT(OUT):: PHI0
        INTEGER       ,INTENT(OUT):: NJK
        REAL(8)       ,INTENT(OUT):: VBARRIER
        LOGICAL(4)    ,INTENT(OUT):: TCHK
        REAL(8)                   :: PI, KCALBYMOL
        INTEGER                   :: I,J,ITORS
!     **********************************************************************
!                      CALL TRACE$PUSH('AMBER_TORSIONPARMS')
        PI=4.D0*DATAN(1.D0)
        ITORS = 0
!     ---- Here I try to find the right parameters. The problem is that there are
!          more than one possibilities. E.g. one can have X-CT-CT-X or OS-CT-CT-OH.
!          I also don't know if it is possible that one gets X-CA-C-X instead of X-C-CA-X
!          so I also check the entries again. Maybe this can be optimized later.
!print*,"FLAG: TYPES= ",Type1, TYPE2, TYPE3, TYPE4
        DO I=1, SIZE(TORSION_PARMS)
           IF(TRIM(ADJUSTL(TYPE2)).EQ.TORSION_PARMS(I)%ATOM2.AND.TRIM(ADJUSTL(TYPE3))&
                &.EQ.TORSION_PARMS(I)%ATOM3) THEN
              ITORS= I
              DO J=1, SIZE(TORSION_PARMS)
                 IF(TRIM(ADJUSTL(TYPE1)).EQ.TORSION_PARMS(J)%ATOM1.AND.&
                      & TRIM(ADJUSTL(TYPE2)).EQ.TORSION_PARMS(J)%ATOM2.AND.&
                      & TRIM(ADJUSTL(TYPE3)).EQ.TORSION_PARMS(J)%ATOM3.AND.&
                      & TRIM(ADJUSTL(TYPE4)).EQ.TORSION_PARMS(J)%ATOM4) THEN
                    ITORS= J
                 END IF
              END DO
           END IF
           IF(TRIM(ADJUSTL(TYPE2)).EQ.TORSION_PARMS(I)%ATOM3.AND.TRIM(ADJUSTL(TYPE3))&
                &.EQ.TORSION_PARMS(I)%ATOM2) THEN
              ITORS= I
              DO J=1, SIZE(TORSION_PARMS)
                 IF(TRIM(ADJUSTL(TYPE1)).EQ.TORSION_PARMS(J)%ATOM4.AND.&
                      & TRIM(ADJUSTL(TYPE2)).EQ.TORSION_PARMS(J)%ATOM3.AND.&
                      & TRIM(ADJUSTL(TYPE3)).EQ.TORSION_PARMS(J)%ATOM2.AND.&
                      & TRIM(ADJUSTL(TYPE4)).EQ.TORSION_PARMS(J)%ATOM1) THEN
                    ITORS= J
                 END IF
              END DO
           END IF
           IF(ITORS.NE.0) EXIT
        END DO
        IF(ITORS.EQ.0) THEN
           CALL ERROR$MSG('TORSION PARAMETER NOT FOUND.')
           CALL ERROR$CHVAL('TYPE1 ',TYPE1)
           CALL ERROR$CHVAL('TYPE2 ',TYPE2)
           CALL ERROR$CHVAL('TYPE3 ',TYPE3)
           CALL ERROR$CHVAL('TYPE4 ',TYPE4)
           CALL ERROR$STOP('FORCEFIELD$AMBER_TORSIONPARMS')
        END IF

        PHI0     = TORSION_PARMS(ITORS)%PHASE
        PHI0     = PHI0 * PI / 180.d0
        NJK      = TORSION_PARMS(ITORS)%PERIODICY
        VBARRIER = TORSION_PARMS(ITORS)%F_CONST
        CALL CONSTANTS('KCAL/MOL',KCALBYMOL)
        VBARRIER=VBARRIER*KCALBYMOL
! !     -------------------- DEBUG remove this later----------------------
!         print*,"FLAG: CHOOSE TORSION SET #",ITORS
!         WRITE(*,FMT='(4A4,F8.3,I4,F8.3)') TORSION_PARMS(ITORS)
!         print*,"FLAG PHI0=",PHi0, PHI0/pi*180.d0
!         print*,"------------------------------------------"
! !     ------------------------------------------------------------------
!
!     ==================================================================
!     ==  READ STRING IDENTIFYING POTENTIAL                           ==
!     ==================================================================
        IF(LGT(TYPE2,TYPE3)) THEN
           WRITE(ID,FMT='(A1,2A5," PHI=",F8.4," N=",I2," V=",F3.0)') &
                &                'T ',TYPE2,TYPE3,PHI0,NJK,VBARRIER
        ELSE
           WRITE(ID,FMT='(A1,2A5," PHI=",F8.4," N=",I2," V=",F3.0)') &
                &                'T ',TYPE3,TYPE2,PHI0,NJK,VBARRIER
        END IF
!
        TCHK=(VBARRIER.GT.1.D-6)
        IF(TCHK) RETURN
        ID='ZERO POTENTIAL'
        NJK=1
        PHI0=0.D0
        VBARRIER=0.D0
        TCHK=.FALSE.
!                      CALL TRACE$POP
        RETURN
      END SUBROUTINE FORCEFIELD$AMBER_TORSIONPARMS

!     ..................................................................
      SUBROUTINE FORCEFIELD$AMBER_NONBONDPARMS(ATOM1,ATOM2,ID,RIJ,DIJ,TCHK)
!     ******************************************************************
        USE FORCEFIELD_MODULE
        USE CLASSICAL_MODULE,ONLY : POT_TYPE
        IMPLICIT NONE
        CHARACTER(5)  ,INTENT(IN)  :: ATOM1
        CHARACTER(5)  ,INTENT(IN)  :: ATOM2
        CHARACTER(64) ,INTENT(OUT) :: ID
        REAL(8)       ,INTENT(OUT) :: RIJ
        REAL(8)       ,INTENT(OUT) :: DIJ
        LOGICAL(4)    ,INTENT(OUT) :: TCHK
        REAL(8)                    :: XI,XJ,DI,DJ
        REAL(8)                    :: ANGSTROM,KCALBYMOL
        INTEGER                    :: I
        LOGICAL(4)                 :: T1, T2
!     ******************************************************************
!                      CALL TRACE$PUSH('AMBER_NONBONDPARMS')
        T1=.FALSE.
        T2=.FALSE.
        DO I=1, SIZE(VDW_PARMS)
           IF(TRIM(ADJUSTL(ATOM1)).EQ.TRIM(ADJUSTL(VDW_PARMS(I)%ATOM1))) THEN
              XI= VDW_PARMS(I)%RMIN
              DI= VDW_PARMS(I)%EMIN
              T1=.TRUE.
           END IF
           IF(TRIM(ADJUSTL(ATOM2)).EQ.TRIM(ADJUSTL(VDW_PARMS(I)%ATOM1))) THEN
              XJ= VDW_PARMS(I)%RMIN
              DJ= VDW_PARMS(I)%EMIN
              T2=.TRUE.
           END IF
        END DO
        IF((.NOT.T1).OR.(.NOT.T2)) THEN
! print*,"VANDERWAALS"
! print*,"ATOM1", ATOM1
! print*,"ATOM2", ATOM2
! print*,T1,T2
!stop
           CALL ERROR$MSG('ATOM TYPE NOT FOUND')
           CALL ERROR$CHVAL('ATOM1= ',ATOM1)
           CALL ERROR$CHVAL('ATOM2= ',ATOM2)
           CALL ERROR$STOP('FORCEFIELD$AMBER_NONBONDPARMS')
        END IF
        RIJ=XI+XJ
        DIJ=DSQRT(DI*DJ)
        WRITE(ID,FMT='(A1,2A5, " X=",F5.2," D=",F8.4)') &
             &              'N ',ATOM1,ATOM2,RIJ,DIJ
!
!     ==================================================================
!     == CONVERT TO ATOMIC UNITS                                      ==
!     ==================================================================
        CALL CONSTANTS('ANGSTROM',ANGSTROM)
        CALL CONSTANTS('KCAL/MOL',KCALBYMOL)
        RIJ=RIJ*ANGSTROM
        DIJ=DIJ*KCALBYMOL
        TCHK=(DIJ.GT.1.D-8)
!                      CALL TRACE$POP
        RETURN
      END SUBROUTINE FORCEFIELD$AMBER_NONBONDPARMS
!

!     ..................................................................
      SUBROUTINE FORCEFIELD$AMBER_NONBONDPOTA(ID,RIJ,DIJ,POT)
!     ******************************************************************
!     **                                                              **
!     **  VAN DER WAALS POTENTIAL AS FUNCTION OF 1/R                  **
!     **                                                              **
!     ******************************************************************
      USE CLASSICAL_MODULE, ONLY: POT_TYPE
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN)  :: ID
      REAL(8)      ,INTENT(IN)  :: RIJ ! EQUILIBRIUM RADIUS
      REAL(8)      ,INTENT(IN)  :: DIJ ! VAN DER WAALS ENERGY
      TYPE(POT_TYPE),INTENT(OUT):: POT
      REAL(8)      ,PARAMETER   :: X1=0.D0
      REAL(8)      ,PARAMETER   :: X2=1.D0
      INTEGER(4)   ,PARAMETER   :: NX=100
      LOGICAL(4)   ,PARAMETER   :: TLENNARDJONES=.TRUE.
      REAL(8)                   :: A,B,ALPHA
      INTEGER(4)                :: I
      REAL(8)                   :: SVAR,SVAR6,SVAR1,SVAR2
      REAL(8)                   :: R,X,DRDX
      LOGICAL(4)                :: TCHK
integer :: nfilinfo
!                      CALL TRACE$PUSH('AMBER_NONBONDPOTA')
!     ******************************************************************
!     ==  WE USE A LENNARD-JONES 6-12 POTENTIAL =======================
      POT%ID=ID
      POT%NX=NX
      POT%X1=X1
      POT%DX=(X2-X1)/DBLE(POT%NX-1)
      ALLOCATE(POT%VAL(POT%NX))
      ALLOCATE(POT%DER(POT%NX))
      SVAR1= DIJ * RIJ**12
      SVAR2= 2* DIJ * RIJ**6
      DO I=1,POT%NX
         IF(I.NE.1) THEN
            X=POT%X1+POT%DX*DBLE(I-1)
            POT%VAL(I)=SVAR1*x**12 - SVAR2*x**6  ! x=1/r
            POT%DER(I)=12.d0*SVAR1*x**11 - 6.d0*SVAR2*x**5
         ELSE
            POT%VAL(I)=0.D0
            POT%DER(I)=0.D0
         END IF
      ENDDO
!
!PETER
      POT%VAL(NX)=POT%VAL(NX-1)+0.5D0*POT%DER(NX-1)*POT%DX
      POT%DER(NX)=0.D0

!     ==================================================================
!     == CHOP OF THE INNER NEGATIVE DIVERGENCE OF THE POTENTIAL       ==
!     ==================================================================
      TCHK=.FALSE.
      DO I=1,POT%NX
        SVAR1=POT%VAL(I)
        SVAR2=POT%DER(I)
        TCHK=TCHK.OR.(SVAR1.GT.0.D0.AND.SVAR2.LT.0.D0)
        IF(TCHK) THEN
!         POT%VAL(I)=POT%VAL(I-1)
!         POT%DER(I)=0.D0
        END IF
      ENDDO

! !---- DEBUG remove this later
!         CALL FILEHANDLER$UNIT('INFO',nfilinfo)
!         DO I=1,POT%NX
!            X=POT%X1+POT%DX*DBLE(I-1)
!            SVAR6= 1.d0/X
! !           WRITE(nfilinfo,FMT='(F15.10,2D20.10)') SVAR6, POT%VAL(I), POT%DER(I)
!            WRITE(NFILINFO,FMT='(I8,2F15.6)') I, POT%VAL(I), POT%DER(I)
!         END DO
!         print*,"FLAG POT%ID=",POT%ID
!         print*,"FLAG DIJ=",DIJ, "   RIJ=",RIJ
!         print*,"FLAG FORCED STOP AMBER_NONBONDPOTA INFO" 
!         STOP
! !-----------------------------------
!                      CALL TRACE$POP
      RETURN
    END SUBROUTINE FORCEFIELD$AMBER_NONBONDPOTA
       

! !*************************************************************************************
! !******   TIP3P SPECIFIC ROUTINES                *************************************
! !*************************************************************************************
     
! !     ................................................................
!     SUBROUTINE FORCEFIELD$TIP3P(NITER_)
! !     ********************************************************************
! !     **                                                                **
! !     **    SORT THE ATOMS, WRITE THEM INTO NEW ARRAYS R0_T3P(:,:,:),   **
! !     **                                                                **
! !     **    CALCULATE THE CONSTRAINTS WITH THE LAGRANGIAN MULTIPLIERS   **
! !     **                                                                **
! !     **    AND WRITE THE COORDINATES BACK TO THE NORMAL ARRAYS         **
! !     **                                                                **
! !     ********************************************************************   
!       USE CLASSICAL_MODULE
!       USE PERIODICTABLE_MODULE
!       IMPLICIT NONE
!       INTEGER(4)  ,INTENT(IN)   :: NITER_ ! Number of iterations for the constraints
!       INTEGER(4)  ,ALLOCATABLE,  SAVE  :: b(:,:) ! intern variable to sort the atoms to groups - useful for the constraints
!       INTEGER(4)  ,SAVE         :: b_temp(3)
!       INTEGER(4)                :: i,j,k,n,m
!       REAL(8)                   :: mh       !=1837.94d0          ! hydrogen mass
!       REAL(8)                   :: mo       !=29174.74d0         ! oxygen mass
!       REAL(8)     ,PARAMETER    :: dOH=(0.9572d0*1.889726d0)     ! OH distance
!       REAL(8)     ,PARAMETER    :: dHH=(1.514207d0*1.889726d0)   ! HH distance
!       REAL(8)     ,SAVE         :: C(3)                          ! Constraints
!       REAL(8)     ,SAVE         :: LAMBDA(3)                     ! Lagrange multiplier
!       REAL(8)     ,SAVE         :: nablag10(9), nablag20(9), nablag30(9)  ! gradient(0)
!       REAL(8)     ,SAVE         :: nablag1p(9), nablag2p(9), nablag3p(9)  ! gradient(+)
!       INTEGER(4)                :: NMOL                          ! Number of molecules
!       REAL(8)     ,SAVE         :: RP9(9)                        ! intern variable for the coordinates
!       REAL(8)     ,SAVE         :: MM(9,9), MMINV(9,9)           ! mass matrix
!       REAL(8)     ,SAVE         :: A(3,3), AINV(3,3)             ! system matrix for calculating lambda
!       INTEGER(4)                :: NITER                         ! # of interations for the constraints
!       LOGICAL     ,SAVE         :: INIT=.FALSE.                  ! check if initialisation is already done
!       real(8)     ,parameter    :: tol=1.d-8
!                       CALL TRACE$PUSH('TIP3P')
!       NMOL=(MD%NAT)/3 ! one water molecule has three atoms
!       NITER=NITER_    ! set the number of constraint iterations

!       CALL PERIODICTABLE$GET(1, 'MASS', mh) ! get the hydrogen mass
!       CALL PERIODICTABLE$GET(8, 'MASS', mo) ! get the oxygen mass
      
!       IF(.NOT.INIT) THEN ! if INIT=.FALSE. do the initialisation first
         
!          ALLOCATE(MD%RP(3,MD%NAT))  ! ############################################## ACHTUNG HIERBEI !!!!!!##########
!          ALLOCATE(b(3,MD%NBOND))
!          ALLOCATE(MD%R0_T3P(3,3,(MD%NAT)/3))  ! MD%NAT/3 should be the same as MD%NBOND/2 in the case of pure water
!          ALLOCATE(MD%RP_T3P(3,3,(MD%NAT)/3))
         
!          !***********SORTROUTINE 1**************************************************
!          ! sort the atoms in groups of three: better handling for the constraints
!          !**************************************************************************
!          ! set the entries of the array b to -1
!          do i=1,3
!             do j=1, MD%NBOND
!                b(i,j) = -1
!             end do
!          end do

!          b(1,1) = MD%INDEX2(1,1)
!          b(2,1) = MD%INDEX2(2,1)
!          k=2
!          do j = 2, MD%NBOND        ! j=2, because the first entries are already made
!             do n = 1, (MD%NBOND/2) ! that's a bit dangerouse. for pure water it will work because every molecule has two bonds
!                do m = 1, 2         ! only the first 2 columns have to be checked
!                   if(MD%INDEX2(1,j).eq.b(m,n)) then
!                      b(3,n) = MD%INDEX2(2,j)
!                      goto 4321
!                   elseif(MD%INDEX2(2,j).eq.b(m,n)) then
!                      b(3,n) = MD%INDEX2(1,j)
!                      goto 4321
!                   end if
!                end do
!             end do
!             b(1,k) = MD%INDEX2(1,j)
!             b(2,k) = MD%INDEX2(2,j)
!             k = k + 1
! 4321     end do

!          !***********SORTROUTINE 2**************************************************
!          ! places the oxygen atom on the first position
!          !**************************************************************************
!          do i=1,(MD%NBOND/2)
!             b_temp(1) = b(1,i)
!             b_temp(2) = b(2,i)
!             b_temp(3) = b(3,i)
!             if(MD%TYPE(b(1,i))(1:1).eq.'H') then ! oxygen not on first position
!                if(MD%TYPE(b(2,i))(1:1).eq.'H') then ! oxygen on third position
!                   b(1,i) = b_temp(3)
!                   b(3,i) = b_temp(1)
!                elseif(MD%TYPE(b(2,i))(1:1).eq.'O') then ! oxygen on second position
!                   b(1,i) = b_temp(2)
!                   b(2,i) = b_temp(1)
!                end if
!             end if
!          end do
         
         
!          ! first run: R(+) = R(0)
!          MD%RP(:,:) = MD%R0(:,:)  ! eigentlich sollte das nicht hier geschehen, 
!                                   ! es gibt bestimmt eine andere funktion, die das macht
!          INIT=.TRUE.
!       end if ! end of the initialisation
      
!       !****************************************************************************
!       ! write the spatial coordinates to the new .._T3P arrays
!       !****************************************************************************
!       do i=1, (MD%NAT)/3
!          MD%R0_T3P(1,1,i) = MD%R0(1,b(1,i))  ! oxygen x coordinate
!          MD%R0_T3P(2,1,i) = MD%R0(2,b(1,i))  ! oxygen y coordinate
!          MD%R0_T3P(3,1,i) = MD%R0(3,b(1,i))  ! oxygen z coordinate
!          MD%R0_T3P(1,2,i) = MD%R0(1,b(2,i))  ! hydrogen_1 x coordinate
!          MD%R0_T3P(2,2,i) = MD%R0(2,b(2,i))  ! hydrogen_1 y coordinate
!          MD%R0_T3P(3,2,i) = MD%R0(3,b(2,i))  ! hydrogen_1 z coordinate
!          MD%R0_T3P(1,3,i) = MD%R0(1,b(3,i))  ! hydrogen_2 x coordinate
!          MD%R0_T3P(2,3,i) = MD%R0(2,b(3,i))  ! hydrogen_2 y coordinate
!          MD%R0_T3P(3,3,i) = MD%R0(3,b(3,i))  ! hydrogen_2 z coordinate
         
!          MD%RP_T3P(1,1,i) = MD%RP(1,b(1,i))  ! oxygen x coordinate
!          MD%RP_T3P(2,1,i) = MD%RP(2,b(1,i))  ! oxygen y coordinate
!          MD%RP_T3P(3,1,i) = MD%RP(3,b(1,i))  ! oxygen z coordinate
!          MD%RP_T3P(1,2,i) = MD%RP(1,b(2,i))  ! hydrogen_1 x coordinate
!          MD%RP_T3P(2,2,i) = MD%RP(2,b(2,i))  ! hydrogen_1 y coordinate
!          MD%RP_T3P(3,2,i) = MD%RP(3,b(2,i))  ! hydrogen_1 z coordinate
!          MD%RP_T3P(1,3,i) = MD%RP(1,b(3,i))  ! hydrogen_2 x coordinate
!          MD%RP_T3P(2,3,i) = MD%RP(2,b(3,i))  ! hydrogen_2 y coordinate
!          MD%RP_T3P(3,3,i) = MD%RP(3,b(3,i))  ! hydrogen_2 z coordinate
!       end do
      
!       !******************************************************************************
!       ! initialise the mass matrix MM(:,:)
!       !******************************************************************************
!       do i=1,3
!          MM(i,i) = mo
!       end do
!       do i=4,9
!          MM(i,i) = mh
!       end do
!       !******************************************************************************
!       ! invert the mass matrix MM(:,:) -> MMINV(:,:)
!       !******************************************************************************
!       CALL LIB$INVERTR8(9,MM,MMINV)
      
!       DO i=1,NMOL
!          DO j=1,NITER ! because the three constraints are coupled
!             !******************************************************************************
!             ! write the constraints into the vector C(:)
!             !******************************************************************************
!             CALL TIP3PG(MD%RP_T3P(:,1,i), MD%RP_T3P(:,2,i), DOH, C(1))
!             CALL TIP3PG(MD%RP_T3P(:,1,i), MD%RP_T3P(:,3,i), DOH, C(2))
!             CALL TIP3PG(MD%RP_T3P(:,2,i), MD%RP_T3P(:,3,i), DHH, C(3))
!             !******************************************************************************
!             ! write the derivations of the constraints into the 
!             ! nine-dimensional vectors nablag..(:)
!             !******************************************************************************
!             CALL TIP3PNABLAG(NMOL, i, MD%R0_T3P, nablag10, nablag20, nablag30)
!             CALL TIP3PNABLAG(NMOL, i, MD%R0_T3P, nablag1p, nablag2p, nablag3p)
!             !******************************************************************************
!             ! fill and invert the system matrix A(:,:)
!             !******************************************************************************
!             A(1,1) = dot_product(nablag10,matmul(MMINV(:,:),nablag1p))
!             A(1,2) = dot_product(nablag20,matmul(MMINV(:,:),nablag1p))
!             A(1,3) = dot_product(nablag30,matmul(MMINV(:,:),nablag1p))
!             A(2,1) = dot_product(nablag10,matmul(MMINV(:,:),nablag2p))
!             A(2,2) = dot_product(nablag20,matmul(MMINV(:,:),nablag2p))
!             A(2,3) = dot_product(nablag30,matmul(MMINV(:,:),nablag2p))
!             A(3,1) = dot_product(nablag10,matmul(MMINV(:,:),nablag3p))
!             A(3,2) = dot_product(nablag20,matmul(MMINV(:,:),nablag3p))
!             A(3,3) = dot_product(nablag30,matmul(MMINV(:,:),nablag3p))
!             CALL LIB$INVERTR8(3,A,AINV)
!             LAMBDA(:) = matmul(AINV(:,:),C(:))
!             write(*,*) 'NITER   LAMBDA', j, LAMBDA
!             !stop
!             !******************************************************************************
!             ! correct the (+) coordinates with the constraints term
!             !******************************************************************************
!             RP9(1)=MD%RP_T3P(1,1,i)
!             RP9(2)=MD%RP_T3P(2,1,i)
!             RP9(3)=MD%RP_T3P(3,1,i)
!             RP9(4)=MD%RP_T3P(1,2,i)
!             RP9(5)=MD%RP_T3P(2,2,i)
!             RP9(6)=MD%RP_T3P(3,2,i)
!             RP9(7)=MD%RP_T3P(1,3,i)
!             RP9(8)=MD%RP_T3P(2,3,i)
!             RP9(9)=MD%RP_T3P(3,3,i)
            
!             RP9(:)=RP9(:) - LAMBDA(1)*matmul(MMINV(:,:),nablag10(:)) &
!                  &- LAMBDA(2)*matmul(MMINV(:,:),nablag20(:)) &
!                  &- LAMBDA(3)*matmul(MMINV(:,:),nablag30(:))
            
!             MD%RP_T3P(1,1,i)=RP9(1)
!             MD%RP_T3P(2,1,i)=RP9(2)
!             MD%RP_T3P(3,1,i)=RP9(3)
!             MD%RP_T3P(1,2,i)=RP9(4)
!             MD%RP_T3P(2,2,i)=RP9(5)
!             MD%RP_T3P(3,2,i)=RP9(6)
!             MD%RP_T3P(1,3,i)=RP9(7)
!             MD%RP_T3P(2,3,i)=RP9(8)
!             MD%RP_T3P(3,3,i)=RP9(9)
            
!             IF(ABS(LAMBDA(1)).LT.TOL.AND.ABS(LAMBDA(2)).LT.TOL.AND.ABS(LAMBDA(3)).LT.TOL) GOTO 1010
            
            
!          END DO ! end of NITER
! 1010     CONTINUE
         
!       END DO ! end of NMOL
!                       CALL TRACE$POP
      
!     END SUBROUTINE FORCEFIELD$TIP3P

!     SUBROUTINE TIP3PNABLAG(NMOL,IMOL,R,nablag1,nablag2,nablag3)
!       IMPLICIT NONE
!       REAL(8),INTENT(IN)     :: R(3,3,nmol)
!       REAL(8),INTENT(OUT)    :: nablag1(9), nablag2(9), nablag3(9)
!       INTEGER(4)             :: i
!       INTEGER(8),INTENT(IN)  :: nmol, imol

!       DO i=1,3
!          nablag1(i) = 2.0d0 * (R(i,1,imol) - R(i,2,imol))
!          nablag1(i+3) = -2.0d0 * (R(i,1,imol) - R(i,2,imol))
!          nablag1(i+6) = 0
!       END DO

!       DO i=1,3
!          nablag2(i) = 2.0d0 * (R(i,1,imol) - R(i,3,imol))
!          nablag2(i+3) = 0
!          nablag2(i+6) = -2.0d0 * (R(i,1,imol) - R(i,3,imol))
!       END DO

!       DO i=1,3
!          nablag3(i) = 0
!          nablag3(i+3) = 2.0d0 * (R(i,2,imol) - R(i,3,imol))
!          nablag3(i+6) = -2.0d0 * (R(i,2,imol) - R(i,3,imol))
!       END DO

!       RETURN
!     END SUBROUTINE TIP3PNABLAG

!     SUBROUTINE TIP3PG(r1,r2,d12,g)
!       IMPLICIT NONE
!       REAL(8),INTENT(IN)     :: r1(3)
!       REAL(8),INTENT(IN)     :: r2(3)
!       REAL(8),INTENT(IN)     :: d12
!       REAL(8),INTENT(OUT)    :: g
!       REAL(8)                :: d(3)

!       d(:) = r1(:)-r2(:)         
!       g = d(1)**2 + d(2)**2 + d(3)*2 - d12**2

!       RETURN
!     END SUBROUTINE TIP3PG

