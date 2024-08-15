!     ..................................................................
      PROGRAM PAW_TOSTRC
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT none
      TYPE(LL_TYPE)               :: LL_STRC
      INTEGER(4)                  :: NFIL
      INteGER(4)                  :: NFILSTRC,nfilpdb
      CHARACTER(128)              :: ROOTNAME     ! COMMON ROOT OF THE FILENAMES
      CHARACTER(128)              :: OBJECTNAME,SUFFIX
      LOGICAL                     :: TCHK
      CHARACTER(32)               :: STRING
      logical(4)                  :: tcrystal
      INTEGER(4)                  :: i
      
!     ==========================================================================
!     == mpe$ init must be called also for non-parallel codes                 ==
!     ==========================================================================
      CALL MPE$INIT

      CALL LINKEDLIST$NEW(LL_STRC)
!
!     ==================================================================
!     == GET FILE NAME ROOT FROM THE ARGUMENT LIST AND CONSTRUCT      ==
!     == FILE NAMES                                                   ==
!     ==================================================================
      CALL GET_COMMAND_ARGUMENT(1,ROOTNAME)
      IF(ROOTNAME(1:1).EQ.'-') THEN
        TCRYSTAL=(+ROOTNAME(2:2).EQ.+'C')
        CALL GET_COMMAND_ARGUMENT(2,ROOTNAME)
      END IF
      IF(LEN(TRIM(ROOTNAME)).EQ.0) THEN
        STOP 'NO ROOTNAME SUPPLIED'
      END IF
      I=INDEX(ROOTNAME,'/',BACK=.TRUE.) !???
      OBJECTNAME=ROOTNAME(I+1:)
      
      I=INDEX(ROOTNAME,'.',BACK=.TRUE.)
      IF(I.EQ.0) THEN
         SUFFIX="CAR"
      ELSE
         SUFFIX=ROOTNAME(I+1:)
         ROOTNAME=ROOTNAME(1:I-1)
      END IF
     
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL FILEHANDLER$SETFILE('STRC',.TRUE.,-'.PRESTRC')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','FORM','FORMATTED')
      
      IF(+SUFFIX.EQ.'CAR'.OR.+SUFFIX.EQ.+ROOTNAME) THEN
         CALL FILEHANDLER$SETFILE('CAR',.TRUE.,-'.CAR')
         CALL FILEHANDLER$SETSPECIFICATION('CAR','STATUS','OLD')
         CALL FILEHANDLER$SETSPECIFICATION('CAR','POSITION','REWIND')
         CALL FILEHANDLER$SETSPECIFICATION('CAR','ACTION','READ')
         CALL FILEHANDLER$SETSPECIFICATION('CAR','FORM','FORMATTED')
         
         !     ==================================================================
         !     == READ CAR FILE                                                ==
         !     ==================================================================
         CALL FILEHANDLER$UNIT('CAR',NFIL)
         CALL FILEHANDLER$UNIT('STRC',NFILSTRC)
         CALL readCAR(NFIL,NFILSTRC)
      END IF

      IF(+SUFFIX.EQ.'PDB') THEN
         CALL FILEHANDLER$SETFILE('PDB',.TRUE.,-'.PDB')
         CALL FILEHANDLER$SETSPECIFICATION('PDB','STATUS','OLD')
         CALL FILEHANDLER$SETSPECIFICATION('PDB','POSITION','REWIND')
         CALL FILEHANDLER$SETSPECIFICATION('PDB','ACTION','READ')
         CALL FILEHANDLER$SETSPECIFICATION('PDB','FORM','FORMATTED')
         
         !     ==================================================================
         !     == READ PDB FILE                                                ==
         !     ==================================================================
         CALL FILEHANDLER$UNIT('PDB',NFIL)
         CALL FILEHANDLER$UNIT('STRC',NFILSTRC)
         CALL readPDB(NFIL,NFILSTRC)
      END IF

      IF(+SUFFIX.EQ.'PDB_RED') THEN
         CALL FILEHANDLER$SETFILE('PDB_RED',.TRUE.,-'.PDB_RED')
         CALL FILEHANDLER$SETSPECIFICATION('PDB_RED','STATUS','OLD')
         CALL FILEHANDLER$SETSPECIFICATION('PDB_RED','POSITION','REWIND')
         CALL FILEHANDLER$SETSPECIFICATION('PDB_RED','ACTION','READ')
         CALL FILEHANDLER$SETSPECIFICATION('PDB_RED','FORM','FORMATTED')
         
         !     ==================================================================
         !     == READ PDB FILE                                                ==
         !     ==================================================================
         CALL FILEHANDLER$UNIT('PDB_RED',NFIL)
         CALL FILEHANDLER$UNIT('STRC',NFILSTRC)
         CALL readPDB_RED(NFIL,NFILSTRC)
      END IF

      IF(+SUFFIX.EQ.'XSD') THEN
         CALL FILEHANDLER$SETFILE('XSD',.TRUE.,-'.XSD')
         CALL FILEHANDLER$SETSPECIFICATION('XSD','STATUS','OLD')
         CALL FILEHANDLER$SETSPECIFICATION('XSD','POSITION','REWIND')
         CALL FILEHANDLER$SETSPECIFICATION('XSD','ACTION','READ')
         CALL FILEHANDLER$SETSPECIFICATION('XSD','FORM','FORMATTED')
         
         !     ==================================================================
         !     == READ XSD FILE                                                ==
         !     ==================================================================
         CALL FILEHANDLER$UNIT('XSD',NFIL)
         CALL FILEHANDLER$UNIT('STRC',NFILSTRC)
         CALL readXSD(NFIL,NFILSTRC)
      END IF

      IF(+SUFFIX.EQ.'CSSR') THEN
         CALL FILEHANDLER$SETFILE('CSSR',.TRUE.,-'.CSSR')
         CALL FILEHANDLER$SETSPECIFICATION('CSSR','STATUS','OLD')
         CALL FILEHANDLER$SETSPECIFICATION('CSSR','POSITION','REWIND')
         CALL FILEHANDLER$SETSPECIFICATION('CSSR','ACTION','READ')
         CALL FILEHANDLER$SETSPECIFICATION('CSSR','FORM','FORMATTED')
         
         !     ==================================================================
         !     == READ CSSR FILE                                                ==
         !     ==================================================================
         CALL FILEHANDLER$UNIT('CSSR',NFIL)
         CALL FILEHANDLER$UNIT('STRC',NFILSTRC)
         CALL readCSSR(NFIL,NFILSTRC)
      END IF


      CALL ERROR$NORMALSTOP()
      STOP
      END PROGRAM PAW_TOSTRC



!     ..................................................................
      SUBROUTINE READCAR(NFIL,NFILSTRC)
      USE STRINGS_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      TYPE BOND_TYPE
        INTEGER(4) :: ATOM1
        INTEGER(4) :: ATOM2
        REAL(8)    :: BO
      END TYPE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      integer(4)  ,INTENT(IN) :: NFILSTRC
      INTEGER(4)              :: NAT
      REAL(8)                 :: RBAS(3,3)
      CHARACTER(8),ALLOCATABLE:: NAME(:)
      REAL(8)     ,ALLOCATABLE:: R(:,:)
      REAL(8)     ,ALLOCATABLE:: Q(:)
      INTEGER(4)  ,ALLOCATABLE:: NN(:,:)
      INTEGER(4)              :: NBOND
      TYPE(BOND_TYPE),ALLOCATABLE :: BOND(:)
      REAL(8)                 :: PI
      REAL(8)                 :: A,B,C            ! LENGTH OF LATTICE VECTORS
      REAL(8)                 :: ALPHA,BETA,GAMMA ! ANGLES BETWEEN LATTICE VECTORS
      INTEGER(4)              :: IAT,I,J,IBOND,isvar,iat1,iat2,nend,JAT
      LOGICAL(4)              :: OK1,OK2
      CHARACTER(80)           :: line
      character(LEN=*), PARAMETER  :: car_form='(A5,3F15.9,15X,A8)'
      character(5)            :: dummy

      INTEGER(4)                  :: NEL=108
      TYPE ATOMS_TYPE
         CHARACTER(2)             :: ID           ! Atomic Number
         CHARACTER(2)             :: SYMBOL       ! Atom Name
         REAL(8)                  :: RCOV         ! Covalent radius
         INTEGER(4)               :: CONFIG(4)    ! #s,p,d,f-valence electrons
         INTEGER(4)               :: ZV           ! number of total valence electrons
         INTEGER(4)               :: B            ! number of Bonds
         INTEGER(4)               :: NNB          ! number of neighbours
      END TYPE ATOMS_TYPE
      TYPE(ATOMS_TYPE),ALLOCATABLE            :: ATOMS(:)
      REAL(8)                     :: dist_cov
      REAL(8)     ,PARAMETER      :: ANGSTROM=1.8897259926D0

!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      OK1=.TRUE.
      OK2=.FALSE.
      NAT=0
      NEND=0
!
!     ==================================================================
!     == READ NAT                                                     ==
!     ==================================================================
      DO WHILE(ok1)
         READ(NFIL,FMT='(A80)',END=100)line
         IF (LINE(1:5).EQ."!DATE") THEN
            ok2=.true.
         END IF
         DO WHILE(OK2)
            READ(NFIL, FMT='(A80)', END=200)line
            line=adjustl(TRIM(line))
            IF(LINE(1:3).NE."end".AND.LINE(1:3).NE."PBC") THEN
               NAT=NAT+1
            ELSE
               IF(LINE(1:3).EQ."end") THEN
                  NEND=NEND+1
               END IF
            END IF
         ENDDO
      ENDDO
100  ok1=.false.
200  ok2=.false.

      REWIND(NFIL)
      
      READ(NFIL,*)
      READ(NFIL,*)
      READ(NFIL,*)
      READ(NFIL,*)
!     ==================================================================
!     == READ POSITIONS                                              ==
!     ==================================================================
      ALLOCATE(R(3,NAT))
      ALLOCATE(NAME(NAT))
      ALLOCATE(Q(NAT))
      ALLOCATE(ATOMS(NAT))

      Q=0.0d0
      isvar=0
      DO IAT=1,NAT+NEND
       READ(NFIL,'(A80)') line
        if (line(1:3).EQ."end") then 
           ISVAR=ISVAR+1
           cycle
        end if
         I = IAT - ISVAR
        READ(line,car_form) &
     &           DUMMY,R(:,IAT-isvar),NAME(IAT-ISVAR)
        NAME(IAT-ISVAR) = TRIM(ADJUSTL(NAME(IAT-ISVAR)))
        ATOMS(IAT-ISVAR)%ID=NAME(IAT-ISVAR)(1:2)
        J=IACHAR(NAME(iat-isvar)(2:2))
        IF(J.ge.65.AND.J.LE.90.OR.J.GE.97.AND.J.LE.122) THEN
           NAME(IAT-ISVAR) = TRIM(ADJUSTL(NAME(IAT-ISVAR)))//.ITOS.I
        else
           NAME(IAT-ISVAR) = TRIM(ADJUSTL(NAME(IAT-ISVAR)))//"_"//TRIM(ADJUSTL(.ITOS.I))
        END IF
      ENDDO
!     ==================================================================
!     == UPPERCASE                                                    ==
!     ==================================================================
      DO I=1,NAT
         J=IACHAR(ATOMS(I)%ID(2:2))
         IF (J.ge.97.AND.J.le.122) ATOMS(I)%ID(2:2)=ACHAR(J-32)
      ENDDO


!       ==================================================================
!       == GET COVALENT RADII FROM PERIODIC TABLE                       ==
!       ==================================================================
      DO I=1,NAT
         CALL PERIODICTABLE$GET(ATOMS(I)%ID,"R(COV)",ATOMS(I)%RCOV)
         CALL PERIODICTABLE$GET(ATOMS(I)%ID,"OCC(S)",ATOMS(I)%CONFIG(1))
         CALL PERIODICTABLE$GET(ATOMS(I)%ID,"OCC(P)",ATOMS(I)%CONFIG(2))
!         CALL PERIODICTABLE$GET(ATOMS(I)%ID,"OCC(D)",ATOMS(I)%CONFIG(3))
!         CALL PERIODICTABLE$GET(ATOMS(I)%ID,"OCC(F)",ATOMS(I)%CONFIG(4))
         ATOMS(I)%RCOV=ATOMS(I)%RCOV / ANGSTROM
         ATOMS(I)%ZV=SUM(ATOMS(I)%CONFIG)
         ATOMS(I)%B=0
         ATOMS(I)%NNB=0
      ENDDO

      
! !     ==================================================================
! !     == CALCULATE BONDS                                              ==
! !     ==================================================================
      ALLOCATE(BOND(NAT*8))
      NBOND=0
      DO JAT=1, NAT
         DO IAT=JAT,NAT
            dist_cov = SQRT( (R(1,IAT)-R(1,JAT))**2 + (R(2,IAT)-R(2,JAT))**2 &
           &  + (R(3,IAT)-R(3,JAT))**2 )
            IF (dist_cov.LT.1E-8) dist_cov=100.
            IF (dist_cov.LT.(ATOMS(IAT)%RCOV + ATOMS(JAT)%RCOV * 1.2)) THEN  !DRAW BONDS
               NBOND= NBOND + 1
               ATOMS(JAT)%B=ATOMS(JAT)%B + 1
               ATOMS(IAT)%B=ATOMS(IAT)%B + 1
               IF (ATOMS(JAT)%B.GT.8) THEN                          !MAX 8 BONDING NEIGHBOURS
                  CALL ERROR$MSG('MORE THAN EIGHT NEIGHBOURS')
                  CALL ERROR$I4VAL('B',ATOMS(JAT)%B)
                  CALL ERROR$STOP('PAW_TOSTRC_READCAR')
               END IF
            !   ATOMS(JAT)%NNB(ATOMS(JAT)%B) = IAT
               BOND(NBOND)%ATOM1=JAT 
               BOND(NBOND)%ATOM2=IAT
               BOND(NBOND)%BO=1
            END IF
         ENDDO
      ENDDO
 
! !     ===================================================================
! !     == CHECK BOND ORDER                                              ==
! !     ===================================================================
  DO IAT2=1,2
      DO I=1,NBOND
         IF (ATOMS(BOND(I)%ATOM1)%ZV.EQ.1.OR.ATOMS(BOND(I)%ATOM1)%ZV.EQ.7) THEN
            IF(ATOMS(BOND(I)%ATOM1)%B.EQ.1) THEN
               CYCLE
            ELSE
               print*,"WARNING, ATOM ",BOND(I)%ATOM1," HAS NO BONDS."
            END IF
         END IF
         IF (ATOMS(BOND(I)%ATOM1)%ZV.EQ.2.OR.ATOMS(BOND(I)%ATOM1)%ZV.EQ.6) THEN
            IF(ATOMS(BOND(I)%ATOM1)%B.EQ.2.OR.ATOMS(BOND(I)%ATOM1)%B.EQ.6) THEN
               CYCLE
            ELSE
               IF(ATOMS(BOND(I)%ATOM2)%B.NE.ATOMS(BOND(I)%ATOM2)%ZV.OR.ATOMS(BOND(I)%ATOM2)%B.NE.&
                    &(8-ATOMS(BOND(I)%ATOM2)%ZV)) THEN
                  BOND(I)%BO=BOND(I)%BO + 1
                  ATOMS(BOND(I)%ATOM1)%B =ATOMS(BOND(I)%ATOM1)%B+1 
                  ATOMS(BOND(I)%ATOM2)%B =ATOMS(BOND(I)%ATOM2)%B+1
               END IF
            END IF
         END IF
         IF (ATOMS(BOND(I)%ATOM1)%ZV.EQ.3.OR.ATOMS(BOND(I)%ATOM1)%ZV.EQ.5) THEN
            IF(ATOMS(BOND(I)%ATOM1)%B.EQ.3.OR.ATOMS(BOND(I)%ATOM1)%B.EQ.5) THEN
               CYCLE
            ELSE
               IF(ATOMS(BOND(I)%ATOM2)%B.NE.ATOMS(BOND(I)%ATOM2)%ZV.OR.ATOMS(BOND(I)%ATOM2)%B.NE.&
                    &(8-ATOMS(BOND(I)%ATOM2)%ZV)) THEN
                  BOND(I)%BO=BOND(I)%BO + 1
                  ATOMS(BOND(I)%ATOM1)%B =ATOMS(BOND(I)%ATOM1)%B+1 
                  ATOMS(BOND(I)%ATOM2)%B =ATOMS(BOND(I)%ATOM2)%B+1
               END IF
            END IF
         END IF
         IF (ATOMS(BOND(I)%ATOM1)%ZV.EQ.4) THEN
            IF(ATOMS(BOND(I)%ATOM1)%B.EQ.4) THEN
               CYCLE
            ELSE
               IF(ATOMS(BOND(I)%ATOM2)%B.NE.ATOMS(BOND(I)%ATOM2)%ZV.AND.ATOMS(BOND(I)%ATOM2)%B.NE.&
                    &(8-ATOMS(BOND(I)%ATOM2)%ZV)) THEN
                  BOND(I)%BO=BOND(I)%BO + 1
                  ATOMS(BOND(I)%ATOM1)%B =ATOMS(BOND(I)%ATOM1)%B+1 
                  ATOMS(BOND(I)%ATOM2)%B =ATOMS(BOND(I)%ATOM2)%B+1
               END IF
            END IF
         END IF
      ENDDO      
    ENDDO

!     ==================================================================
!     == WRITE PRESTRC FILE                                           ==
!     ==================================================================
      DO IAT=1,NAT
        WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R(:,IAT) &
     &                   ," Q=",Q(IAT)," FFTYPE='"//name(iat)(1:2)//"' !END"
      ENDDO
      DO IBOND=1,NBOND
        IAT1=BOND(IBOND)%ATOM1
        IAT2=BOND(IBOND)%ATOM2
        WRITE(NFILSTRC,fmt='(a,f5.2,a)')"!BOND ATOM1='"//TRIM(NAME(IAT1)) &
     &                  //"' ATOM2='"//TRIM(NAME(IAT2)) &
     &                  //"' BO=",BOND(IBOND)%BO," !END"
      ENDDO        
      RETURN

      END SUBROUTINE READCAR



!     ..................................................................
      SUBROUTINE READCSSR(NFIL,NFILSTRC)
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE BOND_TYPE
        INTEGER(4) :: ATOM1
        INTEGER(4) :: ATOM2
        REAL(8)    :: BO
      END TYPE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      integer(4)  ,INTENT(IN) :: NFILSTRC
      INTEGER(4)              :: NAT
      REAL(8)                 :: RBAS(3,3)
      CHARACTER(8),ALLOCATABLE:: NAME(:)
      REAL(8)     ,ALLOCATABLE:: R(:,:)
      REAL(8)     ,ALLOCATABLE:: Q(:)
      INTEGER(4)  ,ALLOCATABLE:: NN(:,:)
      INTEGER(4)              :: NBOND
      TYPE(BOND_TYPE),ALLOCATABLE :: BOND(:)
      REAL(8)                 :: PI
      REAL(8)                 :: A,B,C            ! LENGTH OF LATTICE VECTORS
      REAL(8)                 :: ALPHA,BETA,GAMMA ! ANGLES BETWEEN LATTICE VECTORS
      INTEGER(4)              :: IAT,I,J,IBOND,isvar,iat1,iat2
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
!
!     ==================================================================
!     == READ HEADER                                                  ==
!     ==================================================================
      READ(NFIL,FMT='(T39,3F8.3/T22,3F8.3)')A,B,C,ALPHA,BETA,GAMMA
      READ(NFIL,FMT='(I4)')NAT
      READ(NFIL,*)
!
!     ==================================================================
!     == CALCULATE LATTICE VECTORS                                    ==
!     ==================================================================
      RBAS(:,:)=0.D0
      RBAS(3,3)=C
      RBAS(3,2)=B*COS(ALPHA*PI/180.D0)
      RBAS(2,2)=SQRT(B**2-RBAS(3,2)**2)
      RBAS(3,1)=A*COS(BETA*PI/180.D0)
      RBAS(2,1)=(A*B*COS(GAMMA*PI/180.D0)-RBAS(3,1)*RBAS(3,2))/RBAS(2,2)
      RBAS(1,1)=SQRT(A**2-RBAS(2,1)**2-RBAS(3,1)**2)
!
!     ==================================================================
!     == READ POSITIONS                                              ==
!     ==================================================================
      ALLOCATE(R(3,NAT))
      ALLOCATE(NAME(NAT))
      ALLOCATE(Q(NAT))
      ALLOCATE(NN(8,NAT))
      DO IAT=1,NAT
        READ(NFIL,FMT='(I4,1X,A5,3F10.5,1X,8I4,F8.3)') &
     &           ISVAR,NAME(IAT),R(:,IAT),NN(:,IAT),Q(IAT)
       i=iachar(name(iat)(2:2))
       if(i.ge.48.and.i.le.57) then
         name(iat)=name(iat)(1:1)//'_'//name(iat)(2:)
       end if
      ENDDO
!
!     ==================================================================
!     == REMOVE BONDS THAT ARE COUNTED TWICE                          ==
!     ==================================================================
      NBOND=0
      DO IAT=1,NAT
        DO I=1,8
          IF(NN(I,IAT).LT.IAT) NN(I,IAT)=0
          IF(NN(i,IAT).NE.0) NBOND=NBOND+1
        ENDDO
      ENDDO
      ALLOCATE(BOND(NBOND))
      IBOND=0
      DO IAT=1,NAT
        DO I=1,8
          IF(NN(i,IAT).EQ.0) CYCLE
          IBOND=IBOND+1
          BOND(IBOND)%ATOM1=IAT
          BOND(IBOND)%ATOM2=NN(I,IAT)
          BOND(IBOND)%BO=0.D0
          BOND(IBOND)%BO=1.D0
          NN(I,IAT)=0
          DO J=I+1,8
            IF(NN(J,IAT).EQ.BOND(IBOND)%ATOM2) THEN
              BOND(IBOND)%BO=BOND(IBOND)%BO+1.D0
              NBOND=NBOND-1
              NN(J,IAT)=0
            END IF              
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     == WRITE CSSR FILE                                              ==
!     ==================================================================
      DO IAT=1,NAT
        WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R(:,IAT) &
     &                   ," Q=",Q(IAT)," FFTYPE='"//name(iat)(1:2)//"' !END"
      ENDDO
      DO IBOND=1,NBOND
        IAT1=BOND(IBOND)%ATOM1
        IAT2=BOND(IBOND)%ATOM2
        WRITE(NFILSTRC,fmt='(a,f5.2,a)')"!BOND ATOM1='"//TRIM(NAME(IAT1)) &
     &                  //"' ATOM2='"//TRIM(NAME(IAT2)) &
     &                  //"' BO=",BOND(IBOND)%BO," !END"
      ENDDO        
      RETURN
      END SUBROUTINE READCSSR






!     ..................................................................
      SUBROUTINE READPDB(NFIL,NFILSTRC)
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE ATOMS_TYPE
         INTEGER(4)           :: ID
         CHARACTER(4)         :: ATOM
         CHARACTER(3)         :: RESNAME
         CHARACTER(1)         :: CHAINID
         INTEGER(4)           :: RESID
         REAL(8)              :: R(3)
         CHARACTER(2)         :: ELEMENT
      END TYPE ATOMS_TYPE
      TYPE(ATOMS_TYPE), ALLOCATABLE :: ATOMS(:)
      TYPE(ATOMS_TYPE), ALLOCATABLE :: HETATMS(:)
      TYPE CONECT_TYPE
         INTEGER(4)           :: CATOM
         INTEGER(4)           :: BOND(4)
      END TYPE CONECT_TYPE
      TYPE(CONECT_TYPE), ALLOCATABLE :: CONECT(:),CONECT2(:)

      INTEGER(4)  ,INTENT(IN) :: NFIL
      integer(4)  ,INTENT(IN) :: NFILSTRC
      CHARACTER(80)           :: LINE
      CHARACTER(5), ALLOCATABLE :: FFTYPE(:)
      INTEGER(4)              :: IAT,NAT,JAT,NHET,NCON,IAT1,IAT2,NCON2
      LOGICAL(4)              :: OK1,UFF,TSWITCH
      CHARACTER(*), PARAMETER :: ATOM_FORM='(6X,I5,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3,22X,A2)'
      CHARACTER(*), PARAMETER :: CON_FORM='(6X,5I5)'
      INTEGER(4)              :: I,J,K,IVAR
      TYPE BOND_TYPE
        INTEGER(4) :: ATOM1
        INTEGER(4) :: ATOM2
        REAL(8)    :: BO
      END TYPE BOND_TYPE
      TYPE(BOND_TYPE), ALLOCATABLE :: BONDS(:)
      INTEGER(4)                   :: NBOND
      CHARACTER(8), ALLOCATABLE    :: NAME(:)
      REAL(8),   ALLOCATABLE       :: Q(:)
      INTEGER(4)                   :: NRES, NPROTON
      REAL(8)                      :: R_H(3)
      LOGICAL(4)                   :: AS_OK,FIRST, TBOND
      LOGICAL(4), ALLOCATABLE      :: TCHK(:)
      INTEGER(4), ALLOCATABLE      :: oxygen(:)
      INTEGER(4)                   :: NOXYGEN, con(4), atom 

      OK1=.true.
      TSWITCH=.false.
!----------------------
!WHICH FORCEFIELD? FIX THIS LATER...
      UFF=.TRUE.
!-----------------------
      REWIND(NFIL)
      NAT=0
      NHET=0
      NCON=0
      DO WHILE(OK1)
         READ(NFIL,FMT='(A80)',END=101) LINE
         IF(LINE(1:4).EQ."ATOM") NAT=NAT+1
         IF(LINE(1:6).EQ."HETATM") NHET=NHET+1
         IF(LINE(1:6).EQ."CONECT") NCON=NCON+1
      ENDDO
101 OK1=.FALSE.
      
      REWIND(NFIL)
      PRINT*,"REPORT: "
      PRINT*,"   ",NAT," ATOMS IN AMINOACIDS FOUND."
      PRINT*,"   ",NHET," HETATOMS FOUND."
      PRINT*,"   ",NCON," CONECT-ENTRIES FOUND."

      ALLOCATE(ATOMS(NAT))
      ALLOCATE(HETATMS(NHET))
      ALLOCATE(CONECT(NCON))

      DO I=1,NCON
         CONECT(I)%CATOM=0
         CONECT(I)%BOND(:)=0
      ENDDO
      I=0
      J=0
      K=0

      OK1=.TRUE.
      DO WHILE(OK1)
         READ(NFIL,FMT='(A80)',END=201) LINE
         IF(LINE(1:4).EQ."ATOM") THEN
            I= I + 1
            READ(LINE, ATOM_FORM) ATOMS(I)%ID, ATOMS(I)%ATOM, ATOMS(I)%RESNAME, ATOMS(I)%CHAINID,&
                 & ATOMS(I)%RESID, ATOMS(I)%R(:), ATOMS(I)%ELEMENT
         END IF
         IF(LINE(1:6).EQ."HETATM") THEN
            J=J+1
            READ(LINE, ATOM_FORM) HETATMS(J)%ID, HETATMS(J)%ATOM, HETATMS(J)%RESNAME, HETATMS(J)%CHAINID,&
                 & HETATMS(J)%RESID, HETATMS(J)%R(:), HETATMS(J)%ELEMENT 
         END IF
         IF(LINE(1:6).EQ."CONECT") THEN
            K= K + 1
            READ(LINE, CON_FORM) CONECT(K)%CATOM, CONECT(K)%BOND(:)
         END IF
      ENDDO
201 OK1=.FALSE.



!   ========================================================================
!   == FILL UP PROTONS AND PREPARE BONDS                                  ==
!   ========================================================================
      ALLOCATE(BONDS((NAT+NHET)*8))
      ALLOCATE(NAME(NAT*8))
      ALLOCATE(Q((NAT+NHET)*8))
      ALLOCATE(FFTYPE((NAT+NHET)*8))
      ALLOCATE(OXYGEN((NAT+NHET)*8))
      FFTYPE(:)=" "
      OXYGEN=0
      NOXYGEN=0
      NBOND=0
      NPROTON=0
      IAT=1
      Q=0.d0

      FIRST=.TRUE.

      DO WHILE(IAT.LE.NAT)
         CALL GET_NRES(ATOMS(IAT)%RESNAME,NRES)
!      ========================================================================
!      ==  CHECK IF AMINOACID HAS BOND TO HETATM                             ==
!      ========================================================================
         IF(ALLOCATED(TCHK)) DEALLOCATE(TCHK)
         ALLOCATE(TCHK(NRES))
         TCHK(:)=.FALSE.
         DO I=1,NRES
            DO J=1,NCON
               DO K=1,NHET
                  IF(HETATMS(K)%ID.EQ.CONECT(J)%CATOM.AND.ATOMS(IAT+I-1)%ID.EQ.CONECT(J)%BOND(1)) THEN
                      TCHK(I)=.TRUE.
                  END IF
                  IF(HETATMS(K)%ID.EQ.CONECT(J)%CATOM.AND.ATOMS(IAT+I-1)%ID.EQ.CONECT(J)%BOND(2)) THEN
                      TCHK(I)=.TRUE.
                   END IF
                  IF(HETATMS(K)%ID.EQ.CONECT(J)%CATOM.AND.ATOMS(IAT+I-1)%ID.EQ.CONECT(J)%BOND(3)) THEN
                      TCHK(I)=.TRUE.
                   END IF
                  IF(HETATMS(K)%ID.EQ.CONECT(J)%CATOM.AND.ATOMS(IAT+I-1)%ID.EQ.CONECT(J)%BOND(4))  THEN
                     TCHK(I)=.TRUE.
                  END IF
               ENDDO
            ENDDO
         ENDDO
!      ========================================================================
!      ==  CHECK IF AMINOACID IS OK                                          ==
!      ========================================================================
         AS_OK=.TRUE.
         DO I=1,NRES-1
            IF(SIZE(ATOMS).GT.(IAT+NRES-1)) THEN
              IF(ATOMS(IAT+I)%RESNAME.NE.ATOMS(IAT)%RESNAME) THEN
                  AS_OK=.FALSE.
                  IAT=IAT+I
                  EXIT
               END IF
            ELSE
               IF(ATOMS(IAT+I)%RESNAME.NE.ATOMS(IAT)%RESNAME) THEN
                 AS_OK=.FALSE.
                  IAT=NAT
                  EXIT 
               END IF
            END IF
         ENDDO
         IF(.NOT.AS_OK) THEN
            print*,"====================================================================="
            print*,"WARNING: AMINOACID ",ATOMS(IAT-1)%RESID,", ",ATOMS(IAT-1)%RESNAME," IS NOT COMPLETE."
            print*,"---------------------"
            print*,"INFORMATION: DAMAGED AMINOACID "
            print*,ATOMS(IAT-1)
            print*,"---------------------"
            print*,"NEXT AMINOACID: "
            print*,ATOMS(IAT)
            
            IF(IAT.NE.NAT) THEN
               print*,"         JUMP TO NEXT RESIDUE. NEW RESID: ",ATOMS(IAT)%RESID
               print*,"CURRENT POSITION: IAT ",IAT
               print*,"====================================================================="
               CYCLE
            END IF
         ENDIF

         IF(IAT.EQ.NAT) EXIT

!    ==================================================================
!    == FOR EACH AMINO ACID                                          ==
!    ==================================================================
!           ------------------------------------------------
!           -- FIRST AMINOACID = TWO MORE PROTONS         --
!           ------------------------------------------------
         FFTYPE(IAT+1)="C_3  "
         FFTYPE(IAT+2)="C_2  "
         FFTYPE(IAT+3)="O_1  "
         IF(FIRST) THEN
            FFTYPE(IAT)="N_3  "
            Q(IAT)=1.0d0
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT)%R + 0.7*(ATOMS(IAT)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)') "!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                 &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT+2)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                 &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            FIRST=.FALSE.
         ELSE 
            FFTYPE(IAT)="N_2  "
            Q(IAT)=0.0d0
            IF(ATOMS(IAT-1)%CHAINID.NE.ATOMS(IAT)%CHAINID) THEN
               print*,"FLAG: NEW CHAIN FOUND. SET TWO MORE PROTONS AND CHARGE +1 FOR NITROGEN"
               print*,"FLAG: ",ATOMS(IAT-1)%CHAINID, ATOMS(IAT-1)%ATOM,ATOMS(IAT)%CHAINID, ATOMS(IAT)%ATOM
               FFTYPE(IAT)="N_3  "
               Q(IAT) = 1.0d0
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT)%R + 0.7*(ATOMS(IAT)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)') "!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT)%R - 0.7*(ATOMS(IAT+4)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
         END IF
!           ------------------------------------------------
!           -- LAST AMINOACID = ADD -O  TO THE END        --
!           ------------------------------------------------
         IF(((IAT+NRES-1).EQ.NAT)) THEN !letzte AS  insgesamt
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("O_"//(.itos.ivar)))
            NOXYGEN=NOXYGEN+1
            R_H = ATOMS(IAT+2)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+2
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            OXYGEN(NOXYGEN)=NPROTON
            BONDS(NBOND)%BO=1
            Q(NAT+NHET+NPROTON)=-1.0d0
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                 &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"1' !END"
!            Q(IAT)=0.0d0
         ELSE
            IF(ATOMS(IAT)%CHAINID.NE.ATOMS(IAT+NRES)%CHAINID) THEN !LETZTE AS IN KETTE
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("O_"//(.itos.ivar)))
               NOXYGEN=NOXYGEN+1
               R_H = ATOMS(IAT+2)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+2
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               OXYGEN(NOXYGEN)=NPROTON
               BONDS(NBOND)%BO=1
               Q(NAT+NHET+NPROTON)=-1.0d0
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"1' !END"
!               Q(IAT)=0.0d0
            END IF
         END IF
         
         !BOND N_1 -- C_2
         NBOND = NBOND+1
         BONDS(NBOND)%ATOM1=IAT
         BONDS(NBOND)%ATOM2=IAT+1
         BONDS(NBOND)%BO=1
         !PROTON TO N_1
         IF(ATOMS(IAT)%RESNAME.NE."PRO") THEN
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT)%R - 0.7*((ATOMS(IAT)%R - ATOMS(IAT+1)%R) + (ATOMS(IAT+2)%R &
                 & - ATOMS(IAT+1)%R))
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                 &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
         END IF
         !PROTON H_2 TO C_2
         IF(.NOT.TCHK(2)) THEN
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+1)%R + 0.7*((ATOMS(IAT)%R - ATOMS(IAT+1)%R) + (ATOMS(IAT+2)%R &
                 & - ATOMS(IAT+1)%R))
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+1
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                 &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
         END IF
         !BOND C_2 -- C_3
         NBOND = NBOND + 1
         BONDS(NBOND)%ATOM1=IAT+1
         BONDS(NBOND)%ATOM2=IAT+2
         BONDS(NBOND)%BO=1
         !BOND C_3 -- O_4
         NBOND = NBOND + 1
         BONDS(NBOND)%ATOM1=IAT+2
         BONDS(NBOND)%ATOM2=IAT+3
         BONDS(NBOND)%BO=2
         !BOND C_2 -- C_5
         IF(ATOMS(IAT)%RESNAME.NE."GLY") THEN
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+1
            BONDS(NBOND)%ATOM2=IAT+4
            BONDS(NBOND)%BO=1
         END IF
!         Q(IAT)=0.0d0
         !   =======================================================================
         !   == BUILD BONDS AND PROTONS FOR EACH KIND OF AMINOACID                ==
         !   =======================================================================
         IF(AS_OK) THEN
            !      ========================================================================
            !      ==  RESIDUE SERINE OR CYSTEINE                                        ==
            !      ========================================================================
            IF((ATOMS(IAT)%RESNAME.EQ."SER").OR.(ATOMS(IAT)%RESNAME.EQ."CYS")) THEN
               FFTYPE(IAT+4)="C_3  "
               IF(ATOMS(IAT)%RESNAME.EQ."SER") FFTYPE(IAT+5)="O_2  "
               IF(ATOMS(IAT)%RESNAME.EQ."CYS") FFTYPE(IAT+5)="S_3+2"
               !           -------------------------------------------------
               !           -- REST OF AMINOACID                           --
               !           -------------------------------------------------
               !PROTON H_3 TO C_5
               IF(.NOT.TCHK(5)) THEN
                  NPROTON=NPROTON+1
                  ivar = nat+nhet+nproton
                  NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
                  R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
                  NBOND = NBOND+1
                  BONDS(NBOND)%ATOM1=IAT+4
                  BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
                  BONDS(NBOND)%BO=1
                  WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                       &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
               END IF
               !PROTON H_4 TO C_5
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
               !BOND C_5 -- O_6
               NBOND = NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=IAT+5
               BONDS(NBOND)%BO=1
               
               !PROTON H_5 TO O_6 or S_6 (CYSTEINE)
               print*,"FLAG: CYSTEINE OR SERINE: ATOM 6 BONDED TO CLUSTER?",TCHK(6)
               IF(.NOT.TCHK(6)) THEN
                  NPROTON=NPROTON+1
                  ivar = nat+nhet+nproton
                  NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
                  R_H = ATOMS(IAT+5)%R + 0.7*(ATOMS(IAT+5)%R - ATOMS(IAT+4)%R)
                  NBOND = NBOND+1
                  BONDS(NBOND)%ATOM1=IAT+5
                  BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
                  BONDS(NBOND)%BO=1
                  WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                       &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
               END IF
               
               IAT=IAT+NRES-1
            END IF !END OF SERINE OR CYSTEINE 
            
!      ========================================================================
!      ==  RESIDUE GLUTAMINEAcid                                             ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."GLU") THEN
            FFTYPE(IAT+4)="C_3  "
            FFTYPE(IAT+5)="C_3  "
            FFTYPE(IAT+6)="C_2  "
            FFTYPE(IAT+7)="O_1  "
            Q(IAT+7)=-1.0d0
            FFTYPE(IAT+8)="O_1  "            
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON H_3 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON H_4 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !PROTON H_5 TO C_6
            IF(.NOT.TCHK(6)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+5)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+5
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON H_6 TO C_6
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+5)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_6 -- C_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !BOND C_7 -- O_8
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=1
            !BOND C_7 -- O_9
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=IAT+8
            BONDS(NBOND)%BO=2
            IAT=IAT+NRES-1
         END IF !END OF GLUTAMINacid

!      ========================================================================
!      ==  RESIDUE VALINE                                                    ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."VAL") THEN
            FFTYPE(IAT+4)="C_3  "
            FFTYPE(IAT+5)="C_3  "
            FFTYPE(IAT+6)="C_3  "
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+4)%R + 0.7*((ATOMS(IAT)%R - ATOMS(IAT+1)%R) + &
                    & (ATOMS(IAT+2)%R - ATOMS(IAT+1)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 1 TO C_6
            IF(.NOT.TCHK(6)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+5)%R + 0.7*((ATOMS(IAT+5)%R - ATOMS(IAT+4)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+5
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_6
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+5)%R + 0.7*((ATOMS(IAT+4)%R - ATOMS(IAT+1)%R))
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 3 TO C_6
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+5)%R + 0.7*((ATOMS(IAT+4)%R - ATOMS(IAT+6)%R))
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
             !PROTON 1 TO C_7
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+6)%R + 0.7*((ATOMS(IAT+4)%R - ATOMS(IAT+1)%R))
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 2 TO C_7
            IF(.NOT.TCHK(7)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+6)%R + 0.7*((ATOMS(IAT+6)%R - ATOMS(IAT+4)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+6
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 3 TO C_7
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+6)%R + 0.7*((ATOMS(IAT+4)%R - ATOMS(IAT+5)%R))
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !BOND C_5 -- C_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1

            IAT= IAT+ NRES-1
         END IF !END OF VALINE

!      ========================================================================
!      ==  RESIDUE LEUCIN                                                    ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."LEU") THEN
            FFTYPE(IAT+4)="C_3"
            FFTYPE(IAT+5)="C_3"
            FFTYPE(IAT+6)="C_3"
            FFTYPE(IAT+7)="C_3"
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON1 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+4)%R + 0.7*((ATOMS(IAT)%R - ATOMS(IAT+1)%R) + &
                    & (ATOMS(IAT+2)%R - ATOMS(IAT+1)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_5
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+4)%R + 0.7*((ATOMS(IAT+1)%R - ATOMS(IAT+2)%R) )
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_6
            IF(.NOT.TCHK(6)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+5)%R + 0.7*((ATOMS(IAT+2)%R - ATOMS(IAT+1)%R) )
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !BOND C_6 -- C_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !BOND C_6 -- C_8
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_7
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+6)%R + 0.7*((ATOMS(IAT+6)%R - ATOMS(IAT+5)%R) )
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 2 TO C_7
            IF(.NOT.TCHK(7)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+6)%R + 0.7*((ATOMS(IAT+2)%R - ATOMS(IAT+1)%R) )
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+6
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 3 TO C_7
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+6)%R - 0.7*((ATOMS(IAT+2)%R - ATOMS(IAT+1)%R) )
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 1 TO C_8
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+7)%R + 0.7*((ATOMS(IAT)%R - ATOMS(IAT+1)%R) )
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 2 TO C_8
            IF(.NOT.TCHK(8)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+7)%R - 0.7*((ATOMS(IAT)%R - ATOMS(IAT+1)%R) )
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+7
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 3 TO C_8
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+7)%R + 0.7*((ATOMS(IAT+7)%R - ATOMS(IAT+5)%R) )
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"

            IAT= IAT+ NRES-1
         END IF !END OF LEUCINE

!      ========================================================================
!      ==  RESIDUE ISOLEUCIN                                                 ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."ILE") THEN
            FFTYPE(IAT+4)="C_3"
            FFTYPE(IAT+5)="C_3"
            FFTYPE(IAT+6)="C_3"
            FFTYPE(IAT+7)="C_3"
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+4)%R + 0.7*((ATOMS(IAT)%R - ATOMS(IAT+1)%R) + &
                    & (ATOMS(IAT+2)%R - ATOMS(IAT+1)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !BOND C_5 -- C_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_7
            IF(.NOT.TCHK(7)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+6)%R + 0.7*((ATOMS(IAT+4)%R - ATOMS(IAT+1)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+6
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_7
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+6)%R + 0.7*((ATOMS(IAT+6)%R - ATOMS(IAT+4)%R))
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 3 TO C_7
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+6)%R + 0.7*((ATOMS(IAT+4)%R - ATOMS(IAT+5)%R))
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 1 TO C_6
            IF(.NOT.TCHK(6)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+5)%R + 0.7*((ATOMS(IAT+5)%R - ATOMS(IAT+4)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+5
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_6
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+5)%R + 0.7*((ATOMS(IAT)%R - ATOMS(IAT+1)%R))
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_6 -- C_8
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_8
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+7)%R + 0.7*((ATOMS(IAT+7)%R - ATOMS(IAT+5)%R))
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 2 TO C_8
            IF(.NOT.TCHK(8)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+7)%R + 0.7*((ATOMS(IAT+5)%R - ATOMS(IAT+4)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+7
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 3 TO C_8
            NPROTON= NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H= ATOMS(IAT+7)%R + 0.7*((ATOMS(IAT+6)%R - ATOMS(IAT+4)%R))
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"

            IAT= IAT+ NRES-1
         END IF !END OF ISOLEUCINE

!      ========================================================================
!      ==  RESIDUE GLUTAMINE                                                 ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."GLN") THEN
            FFTYPE(IAT+4)="C_3"
            FFTYPE(IAT+5)="C_3"
            FFTYPE(IAT+6)="C_2"
            FFTYPE(IAT+7)="O_1"
            FFTYPE(IAT+8)="N_2"
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON H_3 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R - 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON H_4 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !PROTON H_5 TO C_6
            IF(.NOT.TCHK(6)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+5)%R - 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+5
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON H_6 TO C_6
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+5)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_6 -- C_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !BOND C_7 -- O_8
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=2
            !BOND C_7 -- N_9
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=IAT+8
            BONDS(NBOND)%BO=1
            !PROTON 1 TO N_9
            IF(.NOT.TCHK(9)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+8)%R + 0.7*(ATOMS(IAT+8)%R - ATOMS(IAT+6)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+8
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO N_9
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+8)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT+4)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+8
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"

            IAT=IAT+NRES-1
         END IF !END OF GLUTAMINE

!      ========================================================================
!      ==  RESIDUE Tyrosin                                                   ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."TYR") THEN
            FFTYPE(IAT+4)= "C_3"
            FFTYPE(IAT+5)= "C_R"
            FFTYPE(IAT+6)= "C_R"
            FFTYPE(IAT+7)= "C_R"
            FFTYPE(IAT+8)= "C_R"
            FFTYPE(IAT+9)= "C_R"
            FFTYPE(IAT+10)="C_R"
            FFTYPE(IAT+11)="O_2"
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON 1 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !BOND C_6 -- C_7 BO=2
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=2
            !BOND C_7 -- C_9 BO=1
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=IAT+8
            BONDS(NBOND)%BO=1
            !BOND C_9 -- C_11 BO=2
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+8
            BONDS(NBOND)%ATOM2=IAT+10
            BONDS(NBOND)%BO=2
            !BOND C_11 -- C_10 BO=1
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+10
            BONDS(NBOND)%ATOM2=IAT+9
            BONDS(NBOND)%BO=1
            !BOND C_10 -- C_8 BO=2
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+9
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=2
            !BOND C_8 -- C_6 BO=1
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !PROTON TO C_7
            IF(.NOT.TCHK(7)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+6)%R + 0.7*((ATOMS(IAT+8)%R - ATOMS(IAT+10)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+6
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON TO C_10
            IF(.NOT.TCHK(10)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+9)%R - 0.7*((ATOMS(IAT+8)%R - ATOMS(IAT+10)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+9
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON TO C_9
            IF(.NOT.TCHK(9)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+8)%R + 0.7*((ATOMS(IAT+6)%R - ATOMS(IAT+5)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+8
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON TO C_8
            IF(.NOT.TCHK(8)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+7)%R - 0.7*((ATOMS(IAT+6)%R - ATOMS(IAT+5)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+7
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !BOND C_11 -- O_12
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+10
            BONDS(NBOND)%ATOM2=IAT+11
            BONDS(NBOND)%BO=1
            !PROTON TO O_12
            IF(.NOT.TCHK(12)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+11)%R + 0.7*((ATOMS(IAT+10)%R - ATOMS(IAT+8)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+11
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF

            IAT= IAT+ NRES-1
         END IF !END OF Tyrosine

!      ========================================================================
!      ==  RESIDUE PROLINE                                                   ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."PRO") THEN
            FFTYPE(IAT+4)="C_3"
            FFTYPE(IAT+5)="C_3"
            FFTYPE(IAT+6)="C_3"
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON 1 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_6
            IF(.NOT.TCHK(6)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+5)%R - 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+5
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_6
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+5)%R + 0.7*(ATOMS(IAT+5)%R - ATOMS(IAT+6)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_6 -- C_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_7
            IF(.NOT.TCHK(7)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+6)%R + 0.7*(ATOMS(IAT+6)%R - ATOMS(IAT)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+6
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_7
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+6)%R + 0.7*(ATOMS(IAT+3)%R - ATOMS(IAT+2)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"

            !BOND C_7 -- C_1
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            FFTYPE(IAT)="N_3  "

            IAT= IAT+ NRES-1
         END IF !END OF PROLINE

!      ========================================================================
!      ==  RESIDUE LYSINE                                                    ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."LYS") THEN
            FFTYPE(IAT+4)="C_3"
            FFTYPE(IAT+5)="C_3"
            FFTYPE(IAT+6)="C_3"
            FFTYPE(IAT+7)="C_3"
            FFTYPE(IAT+8)="N_3"
            Q(IAT+8)=1.0d0
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON 1 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT+2)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_6
            IF(.NOT.TCHK(6)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+5)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+5
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_6
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+5)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT+2)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_6 -- C_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_7
            IF(.NOT.TCHK(7)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+6)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+6
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_7
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+6)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT+2)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_7 -- C_8
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_8
            IF(.NOT.TCHK(8)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+7)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+7
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_8
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+7)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT+2)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_8 -- N_9
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=IAT+8
            BONDS(NBOND)%BO=1
            !PROTON 1 TO N_9
            IF(.NOT.TCHK(9)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+8)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+8
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO N_9
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+8)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT+2)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+8
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 3 TO N_9
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+8)%R + 0.7*(ATOMS(IAT+7)%R - ATOMS(IAT+6)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+8
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"


            IAT= IAT+ NRES-1
         END IF !END OF LYSINE

!      ========================================================================
!      ==  RESIDUE ALANINE                                                   ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."ALA") THEN
            FFTYPE(IAT+4)="C_3"
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON 1 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT+2)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 3 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*((ATOMS(IAT)%R - ATOMS(IAT+1)%R) + (ATOMS(IAT+2)%R &
                 & - ATOMS(IAT+1)%R) + (ATOMS(IAT+4)%R - ATOMS(IAT+1)%R))
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"

            IAT= IAT+ NRES-1
         END IF !END OF ALANINE

!      ========================================================================
!      ==  RESIDUE ASPARAGINEACID                                            ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."ASP") THEN
            FFTYPE(IAT+4)="C_3"
            FFTYPE(IAT+5)="C_2"
            FFTYPE(IAT+6)="O_1"
            FFTYPE(IAT+7)="O_1"
            Q(IAT+6)=-1.0d0
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON 1 TO C_5
            IF (.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT+2)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*((ATOMS(IAT)%R - ATOMS(IAT+1)%R) + (ATOMS(IAT+2)%R &
                 & - ATOMS(IAT+1)%R) + (ATOMS(IAT+4)%R - ATOMS(IAT+1)%R))
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !BOND C_6 -- O_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !BOND C_6 -- O_8
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=2

            IAT= IAT+ NRES-1
         END IF !END OF ASPARAGINEACID

!      ========================================================================
!      ==  RESIDUE ASPARAGINE                                                ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."ASN") THEN
!           -------------------------------------------------
!           -- N and O can be switched in PDB file         --
!           -------------------------------------------------
            IF(ADJUSTL(ATOMS(IAT+6)%ELEMENT).EQ."O ") tswitch=.true.
!           -------------------------------------------------

            FFTYPE(IAT+4)="C_3"
            FFTYPE(IAT+5)="C_3"
            FFTYPE(IAT+6)="N_2"
            FFTYPE(IAT+7)="O_1"
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON 1 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT+2)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*((ATOMS(IAT)%R - ATOMS(IAT+1)%R) + (ATOMS(IAT+2)%R &
                 & - ATOMS(IAT+1)%R) + (ATOMS(IAT+4)%R - ATOMS(IAT+1)%R))
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1            
            !BOND C_6 -- N_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            if(tswitch) BONDS(NBOND)%BO=2
            !BOND C_6 -- O_8
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=2
            if(tswitch) BONDS(NBOND)%BO=1
            !PROTON 1 TO N_7
            IF(.NOT.TCHK(7)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+7)%R + 0.7*((ATOMS(IAT+5)%R - ATOMS(IAT+4)%R))
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+6
               if(TSWITCH) BONDS(NBOND)%ATOM1=IAT+7
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO N_7
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+7)%R + 0.7*((ATOMS(IAT+5)%R - ATOMS(IAT+6)%R))
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+6
            if(TSWITCH) BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"

            IAT= IAT+ NRES-1
         END IF !END OF ASPARAGINE

!      ========================================================================
!      ==  RESIDUE HISTIDINE                                                 ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."HIS") THEN
            FFTYPE(IAT+4)="C_3"
            FFTYPE(IAT+5)="C_2"
            FFTYPE(IAT+6)="N_2"
            FFTYPE(IAT+7)="C_2"
            FFTYPE(IAT+8)="C_2"
            FFTYPE(IAT+9)="N_2"
            Q(IAT+9)=1.0d0
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON 1 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+4)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !BOND C_6 -- N_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !BOND C_6 -- C_8
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=2
            !BOND C_8 -- N_10
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=IAT+9
            BONDS(NBOND)%BO=1
            !BOND C_9 -- N_10
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+8
            BONDS(NBOND)%ATOM2=IAT+9
            BONDS(NBOND)%BO=2
            !BOND C_9 -- N_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+8
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !PROTON 1 TO N_7
            IF(.NOT.TCHK(7)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+6)%R - 0.7*((ATOMS(IAT+5)%R - ATOMS(IAT+6)%R) + (ATOMS(IAT+8)%R - &
                    & ATOMS(IAT+6)%R))
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+6
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 1 TO C_9
            IF (.NOT.TCHK(9)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+8)%R - 0.7*((ATOMS(IAT+6)%R - ATOMS(IAT+8)%R) + (ATOMS(IAT+9)%R - &
                    & ATOMS(IAT+8)%R))
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+8
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 1 TO N_10
            IF(.NOT.TCHK(10)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+9)%R - 0.7*((ATOMS(IAT+8)%R - ATOMS(IAT+9)%R) + (ATOMS(IAT+7)%R - &
                    & ATOMS(IAT+9)%R))
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+9
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 1 TO C_8
            IF(.NOT.TCHK(8)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+7)%R - 0.7*((ATOMS(IAT+9)%R - ATOMS(IAT+7)%R) + (ATOMS(IAT+5)%R - &
                    & ATOMS(IAT+7)%R))
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+7
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF

            IAT= IAT+ NRES-1
         END IF !END OF HISTIDINE

!      ========================================================================
!      ==  RESIDUE GLYCINE                                                   ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."GLY") THEN
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON 1 TO C_2
            IF(.NOT.TCHK(2)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+1)%R - 0.7*((ATOMS(IAT)%R - ATOMS(IAT+1)%R) + (ATOMS(IAT+2)%R - &
                    & ATOMS(IAT+1)%R))
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+1
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF

            IAT= IAT+ NRES-1
         END IF !END OF GLYCINE

!      ========================================================================
!      ==  RESIDUE ARGININE                                                  ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."ARG") THEN
            FFTYPE(IAT+4)= "C_3"
            FFTYPE(IAT+5)= "C_3"
            FFTYPE(IAT+6)= "C_3"
            FFTYPE(IAT+7)= "N_2"
            FFTYPE(IAT+8)= "C_3"
            FFTYPE(IAT+9)= "N_2"
            Q(IAT+9)=1.0d0
            FFTYPE(IAT+10)="N_2"
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON 1 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_6
            IF(.NOT.TCHK(6)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+5)%R - 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+5
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_6
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+5)%R - 0.7*(ATOMS(IAT)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_6 -- C_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_7
            IF(.NOT.TCHK(7)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+6)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+6
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_7
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+6)%R + 0.7*(ATOMS(IAT)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_7 -- N_8
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=1
            !PROTON 1 TO N_8
            IF(.NOT.TCHK(8)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+7)%R + 0.7*(ATOMS(IAT)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+7
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !BOND N_8 -- C_9
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=IAT+8
            BONDS(NBOND)%BO=1
            !BOND C_9 -- N_10
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+8
            BONDS(NBOND)%ATOM2=IAT+9
            BONDS(NBOND)%BO=2
            !BOND C_9 -- N_11
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+8
            BONDS(NBOND)%ATOM2=IAT+10
            BONDS(NBOND)%BO=1
            !PROTON 1 TO N_10
            IF(.NOT.TCHK(10)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+9)%R + 0.7*(ATOMS(IAT)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+9
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO N_10
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+9)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+9
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 1 TO N_11
            IF(.NOT.TCHK(11)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+10)%R - 0.7*(ATOMS(IAT)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+10
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 1 TO N_11
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+10)%R - 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+10
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            
            IAT= IAT+ NRES-1
         END IF !END OF GLYCINE

!      ========================================================================
!      ==  RESIDUE METIONINE                                                   ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."MET") THEN
            FFTYPE(IAT+4)="C_3"
            FFTYPE(IAT+5)="C_3"
            FFTYPE(IAT+6)="S_3+2"
            FFTYPE(IAT+7)="C_3"
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON 1 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_6
            IF(.NOT.TCHK(6)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+5)%R - 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+5
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_6
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+5)%R - 0.7*(ATOMS(IAT)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_6 -- S_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !BOND S_7 -- C_8
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_8
            IF(.NOT.TCHK(8)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+7)%R - 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+7
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_8
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+7)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 3 TO C_8
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+7)%R + 0.7*(ATOMS(IAT+7)%R - ATOMS(IAT+6)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
             
            IAT= IAT+ NRES-1
         END IF !END OF METIONINE

!      ========================================================================
!      ==  RESIDUE THREONINE                                                 ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."THR") THEN
            FFTYPE(IAT+4)="C_3"
            FFTYPE(IAT+5)="O_2"
            FFTYPE(IAT+6)="C_3"
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON  TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R - 0.7*( (ATOMS(IAT+5)%R - ATOMS(IAT+4)%R) + (ATOMS(IAT+6)%R - ATOMS(IAT+4)%R) )
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !BOND C_5 -- O_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !PROTON  TO O_6
            IF(.NOT.TCHK(6)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+5)%R + 0.7*(ATOMS(IAT+5)%R - ATOMS(IAT+4)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+5
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !BOND C_5 -- C_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=1
            !PROTON 1 TO C_7
            IF (.NOT.TCHK(7)) THEN 
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+6)%R + 0.7*(ATOMS(IAT+4)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+6
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_7
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+6)%R - 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !PROTON 3 TO C_7
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+6)%R + 0.7*(ATOMS(IAT)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            
            IAT= IAT+ NRES-1
         END IF !END OF THREONINE

!      ========================================================================
!      ==  RESIDUE TRYPTOPHAN                                                ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."TRP") THEN
            FFTYPE(IAT+4)= "C_3"
            FFTYPE(IAT+5)= "C_2"
            FFTYPE(IAT+6)= "C_2"
            FFTYPE(IAT+7)= "C_R"
            FFTYPE(IAT+8)= "N_2"
            FFTYPE(IAT+9)= "C_R"
            FFTYPE(IAT+10)="C_R"
            FFTYPE(IAT+11)="C_R"
            FFTYPE(IAT+12)="C_R"
            FFTYPE(IAT+13)="C_R"

!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON 1 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON 2 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !BOND C_6 -- C_7
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=2
            !PROTON TO C_7
            IF(.NOT.TCHK(7)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+6)%R + (ATOMS(IAT+5)%R - ATOMS(IAT+7)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+6
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !BOND C_7 -- N_9
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=IAT+8
            BONDS(NBOND)%BO=1
            !PROTON TO N_9
            IF(.NOT.TCHK(9)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+8)%R -0.7* ((ATOMS(IAT+6)%R - ATOMS(IAT+8)%R) + (ATOMS(IAT+9)%R - &
                    & ATOMS(IAT+8)%R))
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+8
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !BOND N_9 -- C_10
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+8
            BONDS(NBOND)%ATOM2=IAT+9
            BONDS(NBOND)%BO=1
            !BOND C_10 -- C_8
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+9
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=2
            !BOND C_6 -- C_8
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=1
            !BOND C_8 -- C_11
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=IAT+10
            BONDS(NBOND)%BO=1
            !PROTON TO C_11
            IF(.NOT.TCHK(11)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+10)%R -0.7* ((ATOMS(IAT+7)%R - ATOMS(IAT+10)%R) + (ATOMS(IAT+12)%R - &
                    & ATOMS(IAT+10)%R))
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+10
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !BOND C_11 -- C_13
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+10
            BONDS(NBOND)%ATOM2=IAT+12
            BONDS(NBOND)%BO=2            
            !PROTON TO C_13
            IF(.NOT.TCHK(13)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+12)%R -0.7* ((ATOMS(IAT+10)%R - ATOMS(IAT+12)%R) + (ATOMS(IAT+13)%R - &
                    & ATOMS(IAT+12)%R))
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+12
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !BOND C_13 -- C_14
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+12
            BONDS(NBOND)%ATOM2=IAT+13
            BONDS(NBOND)%BO=1            
            !PROTON TO C_14
            IF(.NOT.TCHK(14)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+13)%R -0.7* ((ATOMS(IAT+12)%R - ATOMS(IAT+13)%R) + (ATOMS(IAT+11)%R - &
                    & ATOMS(IAT+13)%R))
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+13
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !BOND C_14 -- C_12
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+13
            BONDS(NBOND)%ATOM2=IAT+11
            BONDS(NBOND)%BO=2            
            !PROTON TO C_12
            IF(.NOT.TCHK(12)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+11)%R -0.7* ((ATOMS(IAT+9)%R - ATOMS(IAT+11)%R) + (ATOMS(IAT+13)%R - &
                    & ATOMS(IAT+11)%R))
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+11
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !BOND C_12 -- C_10
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+11
            BONDS(NBOND)%ATOM2=IAT+9
            BONDS(NBOND)%BO=1            


            IAT= IAT+ NRES-1
         END IF !END OF TRYPTOPHAN

!      ========================================================================
!      ==  RESIDUE PHENYLALANINE                                                ==
!      ========================================================================
         IF(ATOMS(IAT)%RESNAME.EQ."PHE") THEN
            FFTYPE(IAT+4)= "C_3"
            FFTYPE(IAT+5)= "C_R"
            FFTYPE(IAT+6)= "C_R"
            FFTYPE(IAT+7)= "C_R"
            FFTYPE(IAT+8)= "C_R"
            FFTYPE(IAT+9)= "C_R"
            FFTYPE(IAT+10)="C_R"
!           -------------------------------------------------
!           -- REST OF AMINOACID                           --
!           -------------------------------------------------
            !PROTON 1 TO C_5
            IF(.NOT.TCHK(5)) THEN
               NPROTON=NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+1)%R - ATOMS(IAT)%R)
               NBOND = NBOND+1
               BONDS(NBOND)%ATOM1=IAT+4
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
               END IF
            !PROTON 2 TO C_5
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = ATOMS(IAT+4)%R + 0.7*(ATOMS(IAT+2)%R - ATOMS(IAT+1)%R)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
              &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            !BOND C_5 -- C_6
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+4
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !BOND C_6 -- C_7 BO=2
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+5
            BONDS(NBOND)%ATOM2=IAT+6
            BONDS(NBOND)%BO=2
            !BOND C_7 -- C_9 BO=1
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+6
            BONDS(NBOND)%ATOM2=IAT+8
            BONDS(NBOND)%BO=1
            !BOND C_9 -- C_11 BO=2
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+8
            BONDS(NBOND)%ATOM2=IAT+10
            BONDS(NBOND)%BO=2
            !BOND C_11 -- C_10 BO=1
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+10
            BONDS(NBOND)%ATOM2=IAT+9
            BONDS(NBOND)%BO=1
            !BOND C_10 -- C_8 BO=2
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+9
            BONDS(NBOND)%ATOM2=IAT+7
            BONDS(NBOND)%BO=2
            !BOND C_8 -- C_6 BO=1
            NBOND = NBOND + 1
            BONDS(NBOND)%ATOM1=IAT+7
            BONDS(NBOND)%ATOM2=IAT+5
            BONDS(NBOND)%BO=1
            !PROTON TO C_7
            IF(.NOT.TCHK(7)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+6)%R + 0.7*((ATOMS(IAT+8)%R - ATOMS(IAT+10)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+6
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON TO C_10
            IF(.NOT.TCHK(10)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+9)%R - 0.7*((ATOMS(IAT+8)%R - ATOMS(IAT+10)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+9
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON TO C_9
            IF(.NOT.TCHK(9)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+8)%R + 0.7*((ATOMS(IAT+6)%R - ATOMS(IAT+5)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+8
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON TO C_8
            IF(.NOT.TCHK(8)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+7)%R - 0.7*((ATOMS(IAT+6)%R - ATOMS(IAT+5)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+7
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            !PROTON TO C_11
            IF(.NOT.TCHK(11)) THEN
               NPROTON= NPROTON+1
               ivar = nat+nhet+nproton
               NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
               R_H= ATOMS(IAT+10)%R + 0.7*((ATOMS(IAT+9)%R - ATOMS(IAT+7)%R))
               NBOND= NBOND + 1
               BONDS(NBOND)%ATOM1=IAT+10
               BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
               BONDS(NBOND)%BO=1
               WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                    &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            END IF
            
            IAT= IAT+ NRES-1
         END IF !END OF PHENYLALANIN

      END IF             
      IAT=IAT+1
      DEALLOCATE(TCHK)

      ENDDO

!   ======================================================================
!   == FILL UP HETATMS (WATER) WITH PROTONS AND MAKE FFTYPE FOR HETATMS ==
!   ======================================================================
      DO IAT=1,NHET
         IF(HETATMS(IAT)%ATOM(1:2).EQ." S") FFTYPE(NAT+IAT)="S_3+2"
         IF(HETATMS(IAT)%ATOM(1:2).EQ."FE") FFTYPE(NAT+IAT)="FE3+2"
         IF(HETATMS(IAT)%ATOM(1:2).EQ."MO") FFTYPE(NAT+IAT)="MO6+6"
         IF(HETATMS(IAT)%ATOM(1:2).EQ."CA") FFTYPE(NAT+IAT)="CA6+2"
 

         IF(HETATMS(IAT)%RESNAME.EQ."PO4") THEN
            FFTYPE(NAT+IAT)="P_3+5"

            FFTYPE(NAT+IAT+1)="O_2"
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=NAT+IAT
            BONDS(NBOND)%ATOM2=NAT+IAT+1
            BONDS(NBOND)%BO=1
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = HETATMS(IAT)%R + (/1.0, 0.0, 0.0/)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=NAT+IAT
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)') "!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                 &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"

            FFTYPE(NAT+IAT+2)="O_1"
            Q(NAT+IAT+2)=-1.0d0
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=NAT+IAT
            BONDS(NBOND)%ATOM2=NAT+IAT+2
            BONDS(NBOND)%BO=1

            FFTYPE(NAT+IAT+3)="O_1"
            Q(NAT+IAT+3)=-1.0d0
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=NAT+IAT
            BONDS(NBOND)%ATOM2=NAT+IAT+3
            BONDS(NBOND)%BO=1

            FFTYPE(NAT+IAT+4)="O_1"
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=NAT+IAT
            BONDS(NBOND)%ATOM2=NAT+IAT+4
            BONDS(NBOND)%BO=2
         END IF
            

         IF(HETATMS(IAT)%RESNAME.EQ."HOH") THEN
            FFTYPE(NAT+IAT)="O_2"
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = HETATMS(IAT)%R + (/1.0, 0.0, 0.0/)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=NAT+IAT
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)') "!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                 &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
            NPROTON=NPROTON+1
            ivar = nat+nhet+nproton
            NAME(IAT) = TRIM(ADJUSTL("H_"//(.itos.ivar)))
            R_H = HETATMS(IAT)%R + (/0.0, 1.0, 0.0/)
            NBOND = NBOND+1
            BONDS(NBOND)%ATOM1=NAT+IAT
            BONDS(NBOND)%ATOM2=NAT+NHET+NPROTON
            BONDS(NBOND)%BO=1
            WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R_H(:) &
                 &                   ," Q=",Q(NAT+NHET+NPROTON)," FFTYPE='"//name(iat)(1:2)//"' !END"
         END IF
      ENDDO

!   ======================================================================
!   == BUILD BONDS FROM CONECT-INFOS                                    ==
!   ======================================================================
      NCON2 = 0
      DO JAT=1,NCON
         DO IAT=1,4
            IF(CONECT(JAT)%BOND(IAT).NE.0) NCON2 = NCON2 + 1
         ENDDO
      ENDDO
      ALLOCATE(CONECT2(NCON2))
      
      DO I=1,NCON
         DO IAT=1,NAT
            IF(ATOMS(IAT)%ID.EQ.CONECT(I)%CATOM)    CONECT(I)%CATOM   = IAT
            IF(ATOMS(IAT)%ID.EQ.CONECT(I)%BOND(1))  CONECT(I)%BOND(1) = IAT
            IF(ATOMS(IAT)%ID.EQ.CONECT(I)%BOND(2))  CONECT(I)%BOND(2) = IAT
            IF(ATOMS(IAT)%ID.EQ.CONECT(I)%BOND(3))  CONECT(I)%BOND(3) = IAT
            IF(ATOMS(IAT)%ID.EQ.CONECT(I)%BOND(4))  CONECT(I)%BOND(4) = IAT
         ENDDO
         DO JAT=1,NHET
            IF(HETATMS(JAT)%ID.EQ.CONECT(I)%CATOM)    CONECT(I)%CATOM   = NAT + JAT
            IF(HETATMS(JAT)%ID.EQ.CONECT(I)%BOND(1))  CONECT(I)%BOND(1) = NAT + JAT
            IF(HETATMS(JAT)%ID.EQ.CONECT(I)%BOND(2))  CONECT(I)%BOND(2) = NAT + JAT
            IF(HETATMS(JAT)%ID.EQ.CONECT(I)%BOND(3))  CONECT(I)%BOND(3) = NAT + JAT
            IF(HETATMS(JAT)%ID.EQ.CONECT(I)%BOND(4))  CONECT(I)%BOND(4) = NAT + JAT
         ENDDO

      ENDDO


      DO I=1, NCON
         DO IAT=1,4
            IF(CONECT(I)%BOND(IAT).NE.0) THEN
               TBOND=.FALSE.
               DO J=1,NBOND
                  IF((BONDS(J)%ATOM1.EQ.CONECT(I)%BOND(IAT)).AND.(BONDS(J)%ATOM2 &
                       & .EQ.CONECT(I)%CATOM)) THEN
                     TBOND=.TRUE.
                  END IF
               ENDDO
               IF(TBOND) THEN
                  CYCLE
               ELSE
                  NBOND = NBOND + 1
                  BONDS(NBOND)%ATOM1=CONECT(I)%CATOM
                  BONDS(NBOND)%ATOM2=CONECT(I)%BOND(IAT)
                  BONDS(NBOND)%BO=1
               END IF
            ELSE
               CYCLE
            END IF
         ENDDO
      ENDDO
               
!   ======================================================================
!   == BUILD BONDS BETWEEN AMINOACIDS                                   ==
!   ======================================================================
!  PRODUCES ERRORS. CONNECTS THE LAST ATOM OF AN AMINOACID WITH THE NEXT AMINOACID.
!  CORRECT: THE THIRD ATOM HAS TO BE CONNECTED TO THE NEXT AMINOACID
!       DO IAT=1,NAT-1
!          IF((ATOMS(IAT)%RESNAME.NE.ATOMS(IAT+1)%RESNAME).AND.(ATOMS(IAT)%CHAINID.EQ.ATOMS(IAT+1)%CHAINID)) THEN
!             NBOND = NBOND + 1
!             BONDS(NBOND)%ATOM1=IAT
!             BONDS(NBOND)%ATOM2=IAT+1
!             BONDS(NBOND)%BO=1
!          END IF
!       ENDDO
      IAT=1
      DO WHILE(IAT.LE.NAT-1)
         CALL GET_NRES(ATOMS(IAT)%RESNAME,NRES)
         IF(ATOMS(IAT)%CHAINID.EQ.ATOMS(IAT+NRES)%CHAINID) THEN
            NBOND= NBOND + 1
            BONDS(NBOND)%ATOM1= IAT+2
            BONDS(NBOND)%ATOM2= IAT+NRES
            BONDS(NBOND)%BO= 1
         END IF
         IAT= IAT + NRES 
      ENDDO
 
!     =====================================================================
!     == WRITE ATOMS                                                     ==
!     =====================================================================
      IF(ALLOCATED(NAME)) DEALLOCATE(NAME)
      ALLOCATE(NAME(NAT))
      DO IAT=1,NAT
         NAME(IAT) = TRIM(ADJUSTL(ATOMS(IAT)%ELEMENT))
         J=IACHAR(NAME(IAT)(2:2))
         IF(J.ge.65.AND.J.LE.90.OR.J.GE.97.AND.J.LE.122) THEN
            NAME(IAT) = TRIM(ADJUSTL(NAME(IAT)))//.ITOS.IAT
         else
            NAME(IAT) = TRIM(ADJUSTL(NAME(IAT)))//"_"//TRIM(ADJUSTL(.ITOS.IAT))
         END IF
         
         WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A,A,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",ATOMS(IAT)%R(:) &
              & ," Q=",Q(IAT)," FFTYPE='",TRIM(ADJUSTL(FFTYPE(IAT))),&
              & "'   RES='"//ATOMS(IAT)%RESNAME//TRIM(ADJUSTL(.ITOS.ATOMS(IAT)%RESID))//"' !END"
      ENDDO
      DEALLOCATE(NAME)

!     =====================================================================
!     == WRITE HETATMS                                                   ==
!     =====================================================================
      ALLOCATE(NAME(NHET))
      DO IAT=1,NHET
         NAME(IAT) = TRIM(ADJUSTL(HETATMS(IAT)%ELEMENT))
         J=IACHAR(NAME(IAT)(2:2))
         ivar=nat+iat
         IF(J.ge.65.AND.J.LE.90.OR.J.GE.97.AND.J.LE.122) THEN
            NAME(IAT) = TRIM(ADJUSTL(NAME(IAT)))//TRIM(ADJUSTL(.ITOS.ivar))
         else
            NAME(IAT) = TRIM(ADJUSTL(NAME(IAT)))//"_"//TRIM(ADJUSTL(.ITOS.ivar))
         END IF
         
         WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A,A5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",HETATMS(IAT)%R(:) &
              &  ," Q=",Q(NAT+IAT)," FFTYPE='",FFTYPE(NAT+IAT),"'   RES='"//HETATMS(IAT)%RESNAME//"' !END"
      ENDDO
      DEALLOCATE(NAME)
      DEALLOCATE(Q)



!     =====================================================================
!     == WRITE BONDS                                                     ==
!     =====================================================================

print*,"==================================="
print*,"PRINTING BONDS"
         IF(ALLOCATED(NAME)) DEALLOCATE(NAME)
         ALLOCATE(NAME((NAT+NHET+NPROTON)*8))
         DO i=1,NAT
            NAME(I) = TRIM(ADJUSTL(ATOMS(I)%ELEMENT))
            J=IACHAR(NAME(I)(2:2))
            IF(J.ge.65.AND.J.LE.90.OR.J.GE.97.AND.J.LE.122) THEN
               NAME(I) = TRIM(ADJUSTL(NAME(I)))//(.ITOS.I)
            else
               NAME(I) = TRIM(ADJUSTL(NAME(I)))//"_"//TRIM(ADJUSTL(.ITOS.I))
            END IF
         ENDDO
         DO i=1,NHET
            NAME(NAT+I) = TRIM(ADJUSTL(HETATMS(I)%ELEMENT))
            J=IACHAR(NAME(I)(2:2))
            IF(J.ge.65.AND.J.LE.90.OR.J.GE.97.AND.J.LE.122) THEN
               IVAR=I+NAT
               NAME(NAT+I) = TRIM(ADJUSTL(NAME(NAT+I)))//TRIM(ADJUSTL(.ITOS.IVAR))
            else
               IVAR=I+NAT
               NAME(NAT+I) = TRIM(ADJUSTL(NAME(NAT+I)))//"_"//TRIM(ADJUSTL(.ITOS.IVAR))
            END IF
         ENDDO
         DO i=1,NPROTON
            IVAR=NAT+NHET+I
            NAME(IVAR) = "H_"//TRIM(ADJUSTL(.itos.IVAR))
         ENDDO
         DO i=1,NOXYGEN
            IVAR=NAT+NHET+OXYGEN(i)
            NAME(IVAR) = "O_"//TRIM(ADJUSTL(.itos.IVAR))
         ENDDO
            

      !WRITE BONDS
        DO I=1, NBOND
            IAT1=BONDS(I)%ATOM1
            IAT2=BONDS(I)%ATOM2
            WRITE(NFILSTRC,fmt='(a,f5.2,a)')"!BOND ATOM1='"//TRIM(NAME(IAT1)) &
     &                  //"' ATOM2='"//TRIM(NAME(IAT2)) &
     &                  //"' BO=",BONDS(I)%BO," !END"
         ENDDO
print*,"===DONE==="

      STOP

      END SUBROUTINE READPDB


      SUBROUTINE READPDB_RED(NFIL, NFILSTRC)
        USE STRINGS_MODULE
        USE FORCEFIELD_MODULE, ONLY: PDB_ATOM_TYPE
        IMPLICIT NONE
        INTEGER(4),       INTENT(IN)  :: NFIL, NFILSTRC
        CHARACTER(103)                 :: line
        CHARACTER(LEN=*), PARAMETER   :: pdb_form='(A6,I5,1X,A5,A4,A1,I4,4X,3F8.3,2F6.2,6X,A4,A2,2X,A1,1X,A12,1X,A7)'
        CHARACTER(LEN=*), PARAMETER   :: strc_form='(3X,A12,A11,A5,3F13.5,A5)'
        TYPE(PDB_ATOM_TYPE)           :: pdb_atom
        character(6)                  :: dummy
        CHARACTER(11)                 :: paw_name

        REWIND(NFIL)
        DO
           READ(NFIL,FMT='(A102)', END=301) LINE
           IF(LINE(1:6).EQ.'ATOM  '.OR.LINE(1:6).EQ.'HETATM') THEN
              READ(LINE,pdb_form)  pdb_atom
              IF(SCAN(PDB_ATOM%NAME(1:1),'1234567890').NE.0) THEN
                 paw_name = pdb_atom%name(2:2)//"_"//pdb_atom%name(3:LEN_TRIM(pdb_atom%name))//&
                      & pdb_atom%name(1:1)//"_"&
                      & //TRIM(ADJUSTL(.itos.pdb_atom%id))//"'"
              ELSE IF(PDB_ATOM%NAME(1:1).EQ.' ') THEN
                 paw_name = pdb_atom%name(2:2)//"_"//pdb_atom%name(3:LEN_TRIM(pdb_atom%name))//"_"&
                      & //TRIM(ADJUSTL(.itos.pdb_atom%id))//"'"               
              ELSE
                 paw_name = pdb_atom%name(1:2)//pdb_atom%name(3:LEN_TRIM(pdb_atom%name))//"_"&
                      &//TRIM(ADJUSTL(.itos.pdb_atom%id))//"'"
              END IF
!              WRITE(NFILSTRC,FMT='(A,3F15.5,A)')"   !ATOM NAME='"//TRIM(ADJUSTL(paw_name))//"  R=",pdb_atom%r(:),"' !END"
              WRITE(NFILSTRC,STRC_FORM) "!ATOM NAME='",TRIM(ADJUSTL(paw_name)),"   R=",pdb_atom%r(:)," !END"

           END IF
        END DO
301     RETURN
      END SUBROUTINE READPDB_RED


!     ..................................................................
      SUBROUTINE READXSD(NFIL,NFILSTRC)
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      integer(4)  ,INTENT(IN) :: NFILSTRC

      PRINT*,"NOT YET INCLUDED!"
      STOP

      END SUBROUTINE READXSD

!     ..................................................................
      SUBROUTINE GET_NRES(RESNAME,NRES)
        IMPLICIT NONE
        CHARACTER(3), INTENT(IN)  :: RESNAME
        INTEGER(4),   INTENT(OUT) :: NRES
        
        IF(RESNAME.EQ."GLY") NRES=4
        IF(RESNAME.EQ."ALA") NRES=5
        IF(RESNAME.EQ."VAL") NRES=7
        IF(RESNAME.EQ."LEU") NRES=8
        IF(RESNAME.EQ."ILE") NRES=8
        IF(RESNAME.EQ."MET") NRES=8
        IF(RESNAME.EQ."PRO") NRES=7
        IF(RESNAME.EQ."PHE") NRES=11
        IF(RESNAME.EQ."TRP") NRES=14
        IF(RESNAME.EQ."SER") NRES=6
        IF(RESNAME.EQ."THR") NRES=7
        IF(RESNAME.EQ."ASN") NRES=8
        IF(RESNAME.EQ."GLN") NRES=9
        IF(RESNAME.EQ."TYR") NRES=12
        IF(RESNAME.EQ."CYS") NRES=6
        IF(RESNAME.EQ."LYS") NRES=9
        IF(RESNAME.EQ."ARG") NRES=11
        IF(RESNAME.EQ."HIS") NRES=10
        IF(RESNAME.EQ."ASP") NRES=8
        IF(RESNAME.EQ."GLU") NRES=9

      END SUBROUTINE GET_NRES
        
