!  IS IT TRUE THAT AMBER DOES NOT USE INVERSIONS?
!  I CHANGED IF(TRIM(ADJUSTL(MD%FF)).EQ.'UFF') THEN
!  CHECK EXPRESSION FOR STRESSES (SIGMA) IN CLASSICAL_ETOTAL REGARDING  
!  LATTICE TRANSLATIONS

!........1.........2.........3.........4.........5.........6.........7.........8
MODULE CLASSICAL_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: CLASSICAL                                                          **
!**                                                                           **
!**  PURPOSE: PERFORMS DYNAMICAL SIMULATIONS WITH A CLASSICAL                 **
!**    FIELD                                                                  **
!**                                                                           **
!**  PREPARE DATA                                                             **
!**   1. SET THE FOLLOWING DATA USING CLASSICAL$SET                           **
!**      NAT      NUMBER OF ATOMS                                             **
!**      R(0)                                                                 **
!**      R(-)                                                                 **
!**      MASS     MASS OF THE ATOMS                                           **
!**      QEL      CHARGES ON THE ATOMS                                        **
!**      TYPE     ATOM SYMBOL FOR UFF                                         **
!**      NBOND    NUMBER OF BONDS                                             **
!**      BOND     BONDS                                                       **
!**                                                                           **
!**  ITERATE                                                                  **
!**   1. SELECT THE PROPER DATASET IF ANOTHER COULD BE SELECTED               **
!**      USING CLASSICAL$SELECT                                               **
!**   2. PRODUCE A NEW NEIGHBORLIST FOR NONBONDED INTERACTIONS                **
!**      IF NECCESARY                                                         **
!**      USING CLASSICAL$NEIGHBORS                                            **
!**   3. CALCULATE TOTAL ENERGY AND FORCES IF NECCESARY                       **
!**      USING CLASSICAL$ETOTAL                                               **
!**   -- CHANGE FORCES IF DESIRED                                             **
!**   4. PROPAGATE ATOMIC POSITIONS                                           **
!**      USING CLASSICAL$PROPAGATE                                            **
!**   -- CHANGE R(+) IF DESIRED                                               **
!**   5. SWITCH ATOMIC POSITIONS TO PREPARE FOR THE NEXT STEP                 **
!**      USING CLASSICAL$SWITCH                                               **
!**                                                                           **
!**   THE NONBOND NEIGHBORLIST HAS A HARDWIRED CUTOFF IN                      **
!**   CLASSICAL_NEIGHBORS                                                     **
!**                                                                           **
!**                                                                           **
!**                                                                           **
!*******************************************************************************
TYPE POT_TYPE                   ! POTENTIAL ON THE GRID
  REAL(8)            :: X1      ! FIRST GRID POINT
  REAL(8)            :: DX      ! GRID SPACING
  INTEGER(4)         :: NX      ! #(GRID POINTS)
  REAL(8)   ,POINTER :: VAL(:)  ! VALUES
  REAL(8)   ,POINTER :: DER(:)  ! DERIVATIVES
  CHARACTER(64)      :: ID      ! POTENTIAL IDENTIFIER
END TYPE POT_TYPE
!
TYPE NONBOND_TYPE               ! NONBOND NEIGBORLIST ITEM
  INTEGER(4)         :: IAT1    ! INDEX OF FIRST ATOM
  INTEGER(4)         :: IAT2    ! INDEX OF SECOND ATOM
  INTEGER(4),POINTER :: IT(:)   ! INTEGER LATTICE TRANSLATION OF SECOND ATOM
  LOGICAL(4)         :: EXCLUDE ! BOND EXCLUSION
  LOGICAL(4)         :: TONEFOUR
END TYPE NONBOND_TYPE
!
TYPE BOND_TYPE                  ! BOND ITEM
  INTEGER(4)         :: IAT1    ! INDEX OF FIRST ATOM
  INTEGER(4)         :: IAT2    ! INDEX OF SECOND ATOM
  INTEGER(4),POINTER :: IT2(:)  ! INTEGER LATTICE TRANSLATION OF SECOND ATOM
  INTEGER(4)         :: IPOT    ! POTENTIAL INDEX 
END TYPE BOND_TYPE
!
TYPE ANGLE_TYPE                 ! BOND ANGLE ITEM
  INTEGER(4)         :: IAT1    ! INDEX OF FIRST ATOM (TERMINAL)
  INTEGER(4)         :: IAT2    ! INDEX OF SECOND ATOM (CENTRAL)
  INTEGER(4)         :: IAT3    ! INDEX OF THIRD ATOM  (TERMINAL)
  INTEGER(4),POINTER :: IT1(:)  ! INTEGER LATTICE TRANSLATION OF FIRST ATOM
  INTEGER(4),POINTER :: IT3(:)  ! INTEGER LATTICE TRANSLATION OF THIRD ATOM
  INTEGER(4)         :: IPOT    ! POTENTIAL INDEX 
END TYPE ANGLE_TYPE
!
TYPE TORSION_TYPE               ! BOND-TORSION ITEM
  INTEGER(4)         :: IAT1    ! INDEX OF FIRST ATOM (TERMINAL)
  INTEGER(4)         :: IAT2    ! INDEX OF SECOND ATOM (CENTRAL)
  INTEGER(4)         :: IAT3    ! INDEX OF THIRD ATOM  (CENTRAL)
  INTEGER(4)         :: IAT4    ! INDEX OF THIRD ATOM  (TERMINAL)
  INTEGER(4),POINTER :: IT1(:)  ! INTEGER LATTICE TRANSLATION OF FIRST ATOM
  INTEGER(4),POINTER :: IT3(:)  ! INTEGER LATTICE TRANSLATION OF THIRD ATOM
  INTEGER(4),POINTER :: IT4(:)  ! INTEGER LATTICE TRANSLATION OF FOURTH ATOM
  INTEGER(4)         :: IPOT    ! POTENTIAL INDEX 
END TYPE TORSION_TYPE
!
TYPE INVERSION_TYPE               ! BOND-TORSION ITEM
  INTEGER(4)         :: IAT1    ! INDEX OF FIRST ATOM (CENTRAL)
  INTEGER(4)         :: IAT2    ! INDEX OF SECOND ATOM (TERMINAL)
  INTEGER(4)         :: IAT3    ! INDEX OF THIRD ATOM  (TERMINAL)
  INTEGER(4)         :: IAT4    ! INDEX OF THIRD ATOM  (TERMINAL)
  INTEGER(4),POINTER :: IT2(:)  ! INTEGER LATTICE TRANSLATION OF SECOND ATOM
  INTEGER(4),POINTER :: IT3(:)  ! INTEGER LATTICE TRANSLATION OF THIRD  ATOM
  INTEGER(4),POINTER :: IT4(:)  ! INTEGER LATTICE TRANSLATION OF FOURTH ATOM
  INTEGER(4)         :: IPOT    ! POTENTIAL INDEX 
END TYPE INVERSION_TYPE
!
TYPE MD_TYPE
  !== DATA THAT MUST BE SET BEFORE INITIALIZATION ==============================
  INTEGER(4)            :: NAT            !         #(ATOMS)
  REAL(8)      ,POINTER :: R0(:,:)        !(3,NAT)  ACTUAL ATOMIC POSITIONS
  REAL(8)      ,POINTER :: RMASS(:)       !(NAT)    ATOMIC MASSES
  CHARACTER(8) ,POINTER :: TYPE(:)        !(NAT)    ATOM TYPE
  CHARACTER(30),POINTER :: ATNAME(:)      !(NAT)    ATOM NAME
  CHARACTER(2), POINTER :: ELEMENT(:)     !(NAT)    ELEMENT SYMBOL
  REAL(8)      ,POINTER :: QEL(:)         !(NAT)    POINT CHARGE
  INTEGER(4)            :: NBOND               
  INTEGER(4)   ,POINTER :: INDEX2(:,:)    ! (3,NBOND)
  REAL(8)      ,POINTER :: R0_T3P(:,:,:)  !(3,3,NAT/3)
  REAL(8)      ,POINTER :: RP_T3P(:,:,:)  !(3,3,NAT/3) 
  CHARACTER(1) ,POINTER :: FIXMHK(:)      ! WHICH QM-MM ATOMS ARE FIXED?        
  TYPE(BOND_TYPE),POINTER :: BOND(:)      !(NBOND)
  LOGICAL(4)            :: MOVECELL       ! SWITCH FOR VARIABLE CELL SHAPE
  REAL(8)      ,POINTER :: RBAS0(:,:)     ! LATTICE VECTORS
  !== DATA SET DURING INITIALIZATION ===========================================
  INTEGER(4)            :: NNBX           ! DIMENSION OF NEIGHBORLIST
  INTEGER(4)            :: NNB            ! ACTUAL SIZE OF NEIGHBORLIST
  TYPE(NONBOND_TYPE),POINTER :: IBLIST(:) ! NEIGHBOR LIST
  INTEGER(4)            :: MAXDIV         ! MAX TRANSLATIONS FOR EXCLUSION FILE
  INTEGER(4)            :: NEXCL          ! #(EXCLUSIONS)
  INTEGER(4)   ,POINTER :: EXCLUSION(:,:) ! EXCLUSIONS
  INTEGER(4)   ,POINTER :: I2FIRST(:)     !(NAT)
  INTEGER(4)   ,POINTER :: ITYPE(:)       !(NAT)          
  REAL(8)      ,POINTER :: BO(:)          !(NBOND)           
  INTEGER(4)            :: NANGLE              
  TYPE(ANGLE_TYPE),POINTER :: ANGLE(:)    !(NANGLE)
  INTEGER(4)            :: NTORSION            
  TYPE(TORSION_TYPE),POINTER :: TORSION(:)        !NTORSION)
  INTEGER(4)            :: NINVERSION          
  TYPE(INVERSION_TYPE),POINTER :: INVERSION(:)      !(NINVERSION)
  INTEGER(4)            :: NTYPE               
  INTEGER(4)   ,POINTER :: NONBOND(:,:)   !(NTYPE,NTYPE)
  INTEGER(4)            :: NPOT                
  TYPE(POT_TYPE),POINTER:: POT(:)         !(NPOT)
  LOGICAL(4)            :: TLONGRANGE          
  REAL(8)               :: SCALEONEFOUR   ! SCALEFACTOR OF 1-4 VDW INTERACTIONS
                                          ! (E.G. AMBER)
  LOGICAL(4)   ,POINTER :: TFREEZE(:)     !(NAT)        
  LOGICAL(4)            :: TMMSPHERE=.FALSE.
!                                     
! == DATA ===UNNING ====UNNING =================================================
  REAL(8)      ,POINTER :: RP(:,:)         !(3,NAT)   
  REAL(8)      ,POINTER :: RM(:,:)         !(3,NAT)   
  REAL(8)      ,POINTER :: FORCE(:,:)      !(3,NAT)
  REAL(8)      ,POINTER :: VEL(:)          !(NAT)    
  REAL(8)               :: SIGMA(3,3)
  CHARACTER(32)         :: MDNAME          ! IDENTIFIER FOR SELECTING
  LOGICAL(4)            :: TINI            ! INITALIZED OR NOT
  CHARACTER(16)         :: FF='UFF'        ! TIP3P, UFF OR WHATEVER
  CHARACTER(3) ,POINTER :: RES(:)          ! (NAT)    RESIDUE NAME
  TYPE(MD_TYPE),POINTER :: NEXT            ! LINK TO THE NEXT INSTANCE
  REAL(8)      ,POINTER         ::MMO(:)   !ORIGIN OF THE MM BOX
  REAL(8)      ,POINTER         ::MMVEC(:) !SIZE OF THE MM BOX
  REAL(8)      ,POINTER         ::MMS(:)   !ORIGIN OF THE MM SPHERE
  REAL(8)      ,POINTER         ::MSR(:)   !SIZE OF THE MM SPHERE
END TYPE MD_TYPE
TYPE(MD_TYPE) ,POINTER :: MD
TYPE(MD_TYPE) ,POINTER :: MDFIRST
LOGICAL(4)             :: TINI=.FALSE.
LOGICAL(4)             :: SELECTED=.FALSE.
REAL(8)     ,PARAMETER :: RCLONGRANGE=100.D0
INTEGER(4)             :: LOD  !LEVEL OF DETAIL FOR PRINTOUT
CONTAINS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_ZERO_BOND(N,BOND)
!     **************************************************************************
!     ** INITIALIZES AN BOND ARRAY STRUCTURE                                  **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN)    :: N
      TYPE(BOND_TYPE),INTENT(INOUT) :: BOND(N)
      INTEGER(4)                    :: I
!     **************************************************************************
      DO I=1,N
        BOND(I)%IAT1=0
        BOND(I)%IAT2=0
        NULLIFY(BOND(I)%IT2)
        BOND(I)%IPOT=0
      ENDDO
      RETURN
      END SUBROUTINE CLASSICAL_ZERO_BOND
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_ZERO_ANGLE(N,ANGLE)
!     **************************************************************************
!     ** INITIALIZES AN ANGLE ARRAY STRUCTURE                                 **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(IN)    :: N
      TYPE(ANGLE_TYPE),INTENT(INOUT) :: ANGLE(N)
      INTEGER(4)                     :: I
!     **************************************************************************
      DO I=1,N
        ANGLE(I)%IAT1=0
        ANGLE(I)%IAT2=0
        ANGLE(I)%IAT3=0
        NULLIFY(ANGLE(I)%IT1)
        NULLIFY(ANGLE(I)%IT3)
        ANGLE(I)%IPOT=0
      ENDDO
      RETURN
      END SUBROUTINE CLASSICAL_ZERO_ANGLE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_ZERO_TORSION(N,TORSION)
!     **************************************************************************
!     ** INITIALIZES AN TORSION ARRAY STRUCTURE                               **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)    :: N
      TYPE(TORSION_TYPE),INTENT(INOUT) :: TORSION(N)
      INTEGER(4)                       :: I
!     **************************************************************************
      DO I=1,N
        TORSION(I)%IAT1=0
        TORSION(I)%IAT2=0
        TORSION(I)%IAT3=0
        TORSION(I)%IAT4=0
        NULLIFY(TORSION(I)%IT1)
        NULLIFY(TORSION(I)%IT3)
        NULLIFY(TORSION(I)%IT4)
        TORSION(I)%IPOT=0
      ENDDO
      RETURN
      END SUBROUTINE CLASSICAL_ZERO_TORSION
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_ZERO_INVERSION(N,INVERSION)
!     **************************************************************************
!     ** INITIALIZES AN INVERSION ARRAY STRUCTURE                             **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)          ,INTENT(IN)    :: N
      TYPE(INVERSION_TYPE),INTENT(INOUT) :: INVERSION(N)
      INTEGER(4)                         :: I
!     **************************************************************************
      DO I=1,N
        INVERSION(I)%IAT1=0
        INVERSION(I)%IAT2=0
        INVERSION(I)%IAT3=0
        INVERSION(I)%IAT4=0
        NULLIFY(INVERSION(I)%IT2)
        NULLIFY(INVERSION(I)%IT3)
        NULLIFY(INVERSION(I)%IT4)
        INVERSION(I)%IPOT=0
      ENDDO
      RETURN
      END SUBROUTINE CLASSICAL_ZERO_INVERSION
END MODULE CLASSICAL_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$SELECT(ID_)
!     **************************************************************************
!     **  SELECT A INSTANCE OF THE CLASSICAL OBJECT                           **
!     **  OR CREATE A NEW ONE IF NOT PRESENT                                  **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
!     **************************************************************************
      SELECTED=.TRUE.
!
!     ==========================================================================
!     ==  SELECT OR ALLOCATE                                                  ==
!     ==========================================================================
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        ALLOCATE(MDFIRST)
        MD=>MDFIRST
      ELSE
        MD=>MDFIRST
        DO 
          IF(MD%MDNAME.EQ.ID_) RETURN   ! INSTANCE FOUND. DONE...
          IF(ASSOCIATED(MD%NEXT)) THEN
            MD=>MD%NEXT
          ELSE
            ALLOCATE(MD%NEXT)
            MD=>MD%NEXT
            EXIT 
          END IF
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  SET NEW INSTANCE TO ZERO                                            ==
!     ==========================================================================
      MD%MDNAME=ID_
      MD%NAT=0
      NULLIFY(MD%R0)
      NULLIFY(MD%RMASS)
      NULLIFY(MD%TYPE)
      NULLIFY(MD%ATNAME)
      NULLIFY(MD%ELEMENT)
      NULLIFY(MD%QEL)
      MD%NBOND=0
      NULLIFY(MD%BOND)
      MD%MOVECELL=.FALSE.
      NULLIFY(MD%RBAS0)
      NULLIFY(MD%MMO)
      NULLIFY(MD%MMVEC)
      MD%NNBX=0
      MD%NNB=0
      NULLIFY(MD%IBLIST)
      MD%MAXDIV=0
      MD%NEXCL=0
      NULLIFY(MD%EXCLUSION)
      NULLIFY(MD%I2FIRST)
      NULLIFY(MD%ITYPE)
      NULLIFY(MD%BO)
      MD%NANGLE=0
      NULLIFY(MD%ANGLE)
      MD%NTORSION=0
      NULLIFY(MD%TORSION)
      MD%NINVERSION=0
      NULLIFY(MD%INVERSION)
      MD%NTYPE=0
      NULLIFY(MD%NONBOND)
      MD%NPOT=0
      NULLIFY(MD%POT)
      MD%TLONGRANGE=.TRUE.
      MD%TMMSPHERE=.FALSE.
      NULLIFY(MD%TFREEZE)
      NULLIFY(MD%RP)
      NULLIFY(MD%RM)
      NULLIFY(MD%FORCE)
      NULLIFY(MD%VEL)
      MD%SIGMA(:,:)=0.D0
      MD%TINI=.FALSE.
      NULLIFY(MD%NEXT)
      NULLIFY(MD%RES)
      NULLIFY(MD%MMS)
      NULLIFY(MD%MSR)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$SETR8A(ID_,LENG_,VAL_)
!     **************************************************************************
!     **  CLASSICAL$SET                                                       **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      REAL(8)     ,INTENT(IN) :: VAL_(LENG_)
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO CLASSICAL OBJECT SELECTED')
        CALL ERROR$STOP('CLASSICAL$SETR8A')
      END IF
!
!     ==========================================================================
!     == R(0) = ACTUAL ATOMIC POSITIONS                                       ==
!     ==========================================================================
      IF(ID_.EQ.'R(0)') THEN
        IF(MD%NAT.EQ.0)MD%NAT=LENG_/3
        IF(LENG_.NE.3*MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%R0))THEN
          ALLOCATE(MD%R0(3,MD%NAT))
        END IF
        MD%R0=RESHAPE(VAL_,(/3,MD%NAT/))
!
!     ==========================================================================
!     == CELL(0) = ACTUAL CELL SHAPE                                          ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'CELL(0)') THEN
        IF(.NOT.ASSOCIATED(MD%RBAS0))ALLOCATE(MD%RBAS0(3,3))
        MD%RBAS0=RESHAPE(VAL_,(/3,3/))
!
!     ==========================================================================
!     ==  R(-) = POSITIONS OF THE PREVIOUS TIME STEP                          ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'R(-)') THEN
        IF(MD%NAT.EQ.0)MD%NAT=LENG_/3
        IF(LENG_.NE.3*MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%RM))ALLOCATE(MD%RM(3,MD%NAT))
        MD%RM=RESHAPE(VAL_,(/3,MD%NAT/))
!     ==========================================================================
!     ==  MMO = ORIGIN OF THE MM BOX                                          ==
!     ==  MMVEC = SIZE OF THE MM BOX                                          ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'MMO') THEN
        IF(.NOT.ASSOCIATED(MD%MMO))ALLOCATE(MD%MMO(3))
        MD%MMO=VAL_
      ELSE IF(ID_.EQ.'MMVEC') THEN
        IF(.NOT.ASSOCIATED(MD%MMVEC)) ALLOCATE(MD%MMVEC(9))
        MD%MMVEC=VAL_

!     ==========================================================================
!     ==  MMS = ORIGIN OF MMSPHERE                                            ==
!     ==  MSR = BOUNDARY RADIUS OF MMSPHERE                                   ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'MMS') THEN
        IF(.NOT.ASSOCIATED(MD%MMS))ALLOCATE(MD%MMS(3))
        MD%MMS=VAL_
      ELSE IF(ID_.EQ.'MSR') THEN
        IF(.NOT.ASSOCIATED(MD%MSR))ALLOCATE(MD%MSR(1))
        MD%MSR=VAL_ 
!
!     ==========================================================================
!     ==  R(+) = POSITIONS OF THE NEXT TIME STEP                              ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'R(+)') THEN
        IF(MD%NAT.EQ.0)MD%NAT=LENG_/3
        IF(LENG_.NE.3*MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%RM))ALLOCATE(MD%RP(3,MD%NAT))
        MD%RP=RESHAPE(VAL_,(/3,MD%NAT/))
!
!     ==========================================================================
!     == MASS = ATOMIC MASSES                                                 ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'MASS') THEN
        IF(MD%NAT.EQ.0)MD%NAT=LENG_
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%RMASS))ALLOCATE(MD%RMASS(MD%NAT))
        MD%RMASS=VAL_
!
!     ==========================================================================
!     == QEL =  ATOMIC POINT CHARGE                                           ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'QEL') THEN
        IF(MD%NAT.EQ.0)MD%NAT=LENG_
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%QEL))ALLOCATE(MD%QEL(MD%NAT))
        MD%QEL=VAL_
!
!     ==========================================================================
!     ==  BONDORDER                                                           ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'BONDORDER') THEN
        IF(MD%NBOND.EQ.0)MD%NBOND=LENG_
        IF(LENG_.NE.MD%NBOND) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%BO))ALLOCATE(MD%BO(MD%NBOND))
        MD%BO=VAL_
!
!     ==========================================================================
!     ==  FORCE                                                               ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'FORCE') THEN
        IF(MD%NAT.EQ.0)MD%NAT=LENG_/3
        IF(LENG_.NE.3*MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%FORCE))ALLOCATE(MD%FORCE(3,MD%NAT))
        MD%FORCE=RESHAPE(VAL_,(/3,MD%NAT/))
!
!     ==========================================================================
!     == QEL =  ATOMIC POINT CHARGE                                           ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'VEL') THEN
        IF(MD%NAT.EQ.0)MD%NAT=LENG_
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%VEL))ALLOCATE(MD%VEL(MD%NAT))
        MD%VEL=VAL_
!
!     ==========================================================================
!     == NOT IDENTIFIED                                                       ==
!     ==========================================================================
      ELSE
        CALL ERROR$MSG('INVALID IDENTIFIER')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
        CALL ERROR$STOP('CLASSICAL$SETR8A')
      END IF
      RETURN
    END SUBROUTINE CLASSICAL$SETR8A
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$SETI4A(ID_,LENG_,VAL_)
!     **************************************************************************
!     **  CLASSICAL$SET                                                       **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      INTEGER(4)  ,INTENT(IN) :: VAL_(LENG_)
      INTEGER(4)              :: I,I0,IT(3)
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO CLASSICAL OBJECT SELECTED')
        CALL ERROR$STOP('CLASSICAL$SETI4')
      END IF
!
!     ==========================================================================
!     == INDEX2 = ARRAY DEFINING COVALENT BONDS                               ==
!     ==========================================================================
      IF(ID_.EQ.'INDEX2') THEN
CALL ERROR$MSG('THE INTERFACE INDEX2 SHOULD BE REPLACED BY BOND')
CALL ERROR$STOP('CLASSICAL$SETI4A')
        IF(MD%NBOND.EQ.0)MD%NBOND=LENG_/2
        IF(LENG_.NE.2*MD%NBOND) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$I4VAL('LENG',LENG_)
          CALL ERROR$I4VAL('NBOND',MD%NBOND)
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETI4A')
        END IF
!       == HERE THE NEW IMPLEMENTATION WITH BOND
        IF(.NOT.ASSOCIATED(MD%BOND))ALLOCATE(MD%BOND(MD%NBOND))
        CALL CLASSICAL_ZERO_BOND(MD%NBOND,MD%BOND)
        DO I=1,MD%NBOND
          MD%BOND(I)%IAT1=VAL_(2*(I-1)+1)
          MD%BOND(I)%IAT2=VAL_(2*(I-1)+2)
        ENDDO
!
!     ==========================================================================
!     == BOND  = ARRAY DEFINING COVALENT BONDS                                ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'BOND') THEN
        IF(MD%NBOND.EQ.0)MD%NBOND=LENG_/2
        IF(LENG_.NE.5*MD%NBOND) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$I4VAL('LENG',LENG_)
          CALL ERROR$I4VAL('NBOND',MD%NBOND)
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETI4A')
        END IF
        IF(.NOT.ASSOCIATED(MD%BOND))ALLOCATE(MD%BOND(MD%NBOND))
        CALL CLASSICAL_ZERO_BOND(MD%NBOND,MD%BOND)
        DO I=1,MD%NBOND
          I0=5*(I-1)
          MD%BOND(I)%IAT1=VAL_(I0+1)
          MD%BOND(I)%IAT2=VAL_(I0+2)
          IT(:)=VAL_(I0+3:I0+5)
          IF(IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) THEN
            ALLOCATE(MD%BOND(I)%IT2(3))
            MD%BOND(I)%IT2(:)=IT(:)
          END IF
          MD%BOND(I)%IPOT=0
        ENDDO
      ELSE
        CALL ERROR$MSG('INVALID IDENTIFIER')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
        CALL ERROR$STOP('CLASSICAL$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$SETI4(ID_,VAL_)
!     **************************************************************************
!     **  CLASSICAL$SET                                                       **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: VAL_
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO CLASSICAL OBJECT SELECTED')
        CALL ERROR$STOP('CLASSICAL$SETI4')
      END IF
!
!     ==========================================================================
!     == NAT = #(ATOMS)                                                       ==
!     ==========================================================================
      IF(ID_.EQ.'NAT') THEN
        MD%NAT=VAL_
!
!     ==========================================================================
!     == NBOND = #(BONDS)                                                     ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'NBOND') THEN
        MD%NBOND=VAL_
!
!     ==========================================================================
!     == LOD  LEVEL OF DETAIL FOR THE OUTPUT                                  ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'LOD') THEN
         LOD=VAL_
      ELSE
        CALL ERROR$MSG('INVALID IDENTIFIER')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
        CALL ERROR$STOP('CLASSICAL$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$SETL4A(ID_,LENG_,VAL_)
!     **************************************************************************
!     **  CLASSICAL$SET                                                       **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      LOGICAL(4)  ,INTENT(IN) :: VAL_(LENG_)
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO CLASSICAL OBJECT SELECTED')
        CALL ERROR$STOP('CLASSICAL$SETL4A')
      END IF
!
!     ==========================================================================
!     ==  DEFINES IF AN ATOM IS FROZEN                                        ==
!     ==========================================================================
      IF(ID_.EQ.'TFREEZE') THEN
        IF(MD%NAT.EQ.0)MD%NAT=LENG_
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETL4A')
        END IF
        IF(.NOT.ASSOCIATED(MD%TFREEZE))ALLOCATE(MD%TFREEZE(MD%NAT))
        MD%TFREEZE(:)=VAL_(:)
      ELSE
        CALL ERROR$MSG('INVALID IDENTIFIER')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
        CALL ERROR$STOP('CLASSICAL$SETL4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$SETL4(ID_,VAL_)
!     **************************************************************************
!     **  CLASSICAL$SET                                                       **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL(4)  ,INTENT(IN) :: VAL_
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO CLASSICAL OBJECT SELECTED')
        CALL ERROR$STOP('CLASSICAL$SETL4')
      END IF
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      IF(ID_.EQ.'LONGRANGE') THEN
        MD%TLONGRANGE=VAL_
      ELSE IF(ID_.EQ.'MMSPHERE') THEN 
         MD%TMMSPHERE=VAL_
      ELSE
        CALL ERROR$MSG('INVALID IDENTIFIER')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
        CALL ERROR$STOP('CLASSICAL$SETL4')
      END IF
      RETURN
    END SUBROUTINE CLASSICAL$SETL4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$SETCH(ID_,VAL_)
!     **************************************************************************
!     **  CLASSICAL$SET                                                       **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      CHARACTER(*),INTENT(IN) :: VAL_
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO CLASSICAL OBJECT SELECTED')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$STOP('CLASSICAL$SETL4')
      END IF
      IF(ID_.EQ.'FF') THEN
        IF(VAL_.NE.'UFF'.AND.VAL_.NE.'AMBER') THEN
          CALL ERROR$MSG('UNKNOWN FORCE FIELD ID')
          CALL ERROR$MSG('SUPPORTED VALUES ARE UFF AND AMBER')
          CALL ERROR$CHVAL('FF',VAL_)
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$STOP('CLASSICAL$SETL4')
        END IF   
        MD%FF=VAL_
      ELSE
        CALL ERROR$MSG('INVALID IDENTIFIER')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$STOP('CLASSICAL$SETCH')
      END IF
      RETURN
    END SUBROUTINE CLASSICAL$SETCH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$SETCHA(ID_,LENG_,VAL_)
!     **************************************************************************
!     **  CLASSICAL$SET                                                       **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      CHARACTER(*),INTENT(IN) :: VAL_(LENG_)
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO CLASSICAL OBJECT SELECTED')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$STOP('CLASSICAL$SETCHA')
      END IF
!
!     ==========================================================================
!     == TYPE = ATOM TYPE                                                     ==
!     ==========================================================================
      IF(ID_.EQ.'RES') THEN
         IF(MD%NAT.EQ.0) MD%NAT=LENG_
         IF(LENG_.NE.MD%NAT) THEN
            CALL ERROR$MSG('INCONSISTENT SIZE')
            CALL ERROR$CHVAL('ID',ID_)
            CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
            CALL ERROR$STOP('CLASSICAL$SETCHA')
         END IF
         IF(.NOT.ASSOCIATED(MD%RES))ALLOCATE(MD%RES(MD%NAT))
         MD%RES=VAL_
!     WRITE(*,*) 'RESIDUE NAME:  ', MD%RES

      ELSE IF(ID_.EQ.'TYPE') THEN
        IF(MD%NAT.EQ.0) MD%NAT=LENG_
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETCHA')
        END IF
        IF(.NOT.ASSOCIATED(MD%TYPE))ALLOCATE(MD%TYPE(MD%NAT))
        MD%TYPE=VAL_
      ELSE IF(ID_.EQ.'ATOMNAME') THEN
        IF(MD%NAT.EQ.0) MD%NAT=LENG_
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETCHA')
        END IF
        IF(.NOT.ASSOCIATED(MD%ATNAME))ALLOCATE(MD%ATNAME(MD%NAT))
        MD%ATNAME=VAL_
      ELSE IF(ID_.EQ.'ELEMENT') THEN
        IF(MD%NAT.EQ.0) MD%NAT=LENG_
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$SETCHA')
        END IF
        IF(.NOT.ASSOCIATED(MD%ELEMENT))ALLOCATE(MD%ELEMENT(MD%NAT))
        MD%ELEMENT=VAL_
      ELSE
        CALL ERROR$MSG('INVALID IDENTIFIER')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
        CALL ERROR$STOP('CLASSICAL$SETCHA')
      END IF
      RETURN
      END SUBROUTINE CLASSICAL$SETCHA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$GETR8A(ID_,LENG_,VAL_)
!     **************************************************************************
!     **  CLASSICAL$SET                                                       **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      REAL(8)     ,INTENT(OUT):: VAL_(LENG_)
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO CLASSICAL OBJECT SELECTED')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$STOP('CLASSICAL$GETR8A')
      END IF
!
!     ==========================================================================
!     == R(0) = ACTUAL ATOMIC POSITIONS                                       ==
!     ==========================================================================
      IF(ID_.EQ.'R(0)') THEN
        IF(MD%NAT.EQ.0) MD%NAT=LENG_/3
        IF(LENG_.NE.3*MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$I4VAL('LENG_',LENG_)
          CALL ERROR$I4VAL('MD%NAT',MD%NAT)
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%R0)) THEN
          CALL ERROR$MSG('NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        VAL_=RESHAPE(MD%R0,(/3*MD%NAT/))
!
!     ==========================================================================
!     == R(+) = ATOMIC POSITIONS FOR THE NEXT TIME STEP                       ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'R(+)') THEN
        IF(MD%NAT.EQ.0) MD%NAT=LENG_/3
        IF(LENG_.NE.3*MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$I4VAL('LENG_',LENG_)
          CALL ERROR$I4VAL('MD%NAT',MD%NAT)
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%RP)) THEN
          CALL ERROR$MSG('NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        VAL_=RESHAPE(MD%RP,(/3*MD%NAT/))
!
!     ==========================================================================
!     == CELL(0) = CELL SHAPE OF THE MM SYSTEM                                ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'CELL(0)') THEN
         IF(LENG_.NE.9) THEN
            CALL ERROR$MSG('INCONSISTENT NUMBER OF CELL VECTORS')
            CALL ERROR$CHVAL('ID_',ID_)
            CALL ERROR$I4VAL('LENG_',LENG_)
            CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
            CALL ERROR$STOP('CLASSICAL$GETR8A')
         END IF
         IF(.NOT.ASSOCIATED(MD%RBAS0)) THEN
            CALL ERROR$MSG('NOT ALLOCATED')
            CALL ERROR$CHVAL('ID_',ID_)
            CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
            CALL ERROR$STOP('CLASSICAL$GETR8A')
         END IF
         VAL_=RESHAPE(MD%RBAS0,(/9/))
!
!     ==========================================================================
!     == MASS = ATOMIC MASSES                                                 ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'MASS') THEN
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%RMASS)) THEN
          CALL ERROR$MSG('NOT ALLOCATED')
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        VAL_=RESHAPE(MD%RMASS,(/MD%NAT/))
!
!     ==========================================================================
!     == QEL = ATOMIC CHARGES                                                 ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'QEL') THEN
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$I4VAL('MD%NAT',MD%NAT)
          CALL ERROR$I4VAL('LENG_',LENG_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%QEL)) THEN
          CALL ERROR$MSG('NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        VAL_=RESHAPE(MD%QEL,(/MD%NAT/))
!
!     ==========================================================================
!     ==  R(-) = POSITIONS OF THE PREVIOUS TIME STEP                          ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'R(-)') THEN
        IF(LENG_.NE.3*MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%RM)) THEN
          CALL ERROR$MSG('NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        VAL_=RESHAPE(MD%RM,(/3*MD%NAT/))

!     ==========================================================================
!     ==  MMS = ORIGIN OF MMSPHERE                                            ==
!     ==  MSR = BOUNDARY RADIUS OF MMSPHERE                                   ==
!     ==========================================================================
     ELSE IF(ID_.EQ.'MMS') THEN
        IF(LENG_.NE. SIZE(MD%MMS)) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
       END IF
       IF(.NOT.ASSOCIATED(MD%MMS))THEN
          CALL ERROR$MSG('NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
       END IF
       VAL_=MD%MMS

     ELSE IF(ID_.EQ.'MSR') THEN
        IF(LENG_.NE. SIZE(MD%MSR)) THEN
          CALL ERROR$MSG('NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
       END IF
       IF(.NOT.ASSOCIATED(MD%MSR))THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
       END IF      
       VAL_=MD%MSR       
!
!     ==========================================================================
!     ==  FORCE= POSITIONS OF THE PREVIOUS TIME STEP                          ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'FORCE') THEN
        IF(LENG_.NE.3*MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%FORCE)) THEN
          CALL ERROR$MSG('NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        VAL_=RESHAPE(MD%FORCE,(/3*MD%NAT/))
!
!     ==========================================================================
!     ==  VEL = D(E)/D(VEL)                                                   ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'VEL') THEN
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        IF(.NOT.ASSOCIATED(MD%VEL)) THEN
          CALL ERROR$MSG('NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
          CALL ERROR$STOP('CLASSICAL$GETR8A')
        END IF
        VAL_=MD%VEL
      ELSE
        CALL ERROR$MSG('INVALID IDENTIFIER')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$CHVAL('SELECTION',MD%MDNAME)
        CALL ERROR$STOP('CLASSICAL$GETR8A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$GETI4A(ID_,LENG_,VAL_)
!     **************************************************************************
!     **  CLASSICAL$SET                                                       **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      INTEGER(4)  ,INTENT(OUT):: VAL_(LENG_)
      INTEGER(4)              :: I,I0
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO CLASSICAL OBJECT SELECTED')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$STOP('CLASSICAL$GETI4A')
      END IF
!
!     ==========================================================================
!     == INDEX2 = ARRAY DEFINING COVALENT BONDS                               ==
!     ==========================================================================
      IF(ID_.EQ.'INDEX2') THEN
CALL ERROR$MSG('THE INTERFACE INDEX2 SHOULD BE REPLACED BY BOND')
CALL ERROR$STOP('CLASSICAL$GETI4A')
        IF(LENG_.NE.2*MD%NBOND) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$STOP('CLASSICAL$GETI4A')
        END IF
        IF(.NOT.ASSOCIATED(MD%BOND)) THEN
          CALL ERROR$MSG('BOND NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$STOP('CLASSICAL$GETI4A')
        END IF
        DO I=1,MD%NBOND
          VAL_(2*(I-1)+1)=MD%BOND(I)%IAT1
          VAL_(2*(I-1)+2)=MD%BOND(I)%IAT2
        ENDDO
!
!     ==========================================================================
!     == BOND: ARRAY DEFINING COVALENT BONDS                                  ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'BOND') THEN
        IF(LENG_.NE.5*MD%NBOND) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$STOP('CLASSICAL$GETI4A')
        END IF
        IF(.NOT.ASSOCIATED(MD%BOND)) THEN
          CALL ERROR$MSG('BOND NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$STOP('CLASSICAL$GETI4A')
        END IF
        DO I=1,MD%NBOND
          I0=5*(I-1)
          VAL_(I0+1)=MD%BOND(I)%IAT1
          VAL_(I0+2)=MD%BOND(I)%IAT2
          VAL_(I0+3:I0+5)=MD%BOND(I)%IT2(:)
        ENDDO
      ELSE
        CALL ERROR$MSG('INVALID IDENTIFIER')
        CALL ERROR$STOP('CLASSICAL$GETI4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$GETI4(ID_,VAL_)
!     **************************************************************************
!     **  CLASSICAL$SET                                                       **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(OUT):: VAL_
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO CLASSICAL OBJECT SELECTED')
        CALL ERROR$STOP('CLASSICAL$SETI4')
      END IF
!
!     ==========================================================================
!     == NAT = #(ATOMS)                                                       ==
!     ==========================================================================
      IF(ID_.EQ.'NAT') THEN
        VAL_=MD%NAT
!
!     ==========================================================================
!     == NBOND = #(BONDS)                                                     ==
!     ==========================================================================
      ELSE IF(ID_.EQ.'NBOND') THEN
        VAL_=MD%NBOND
      ELSE
        CALL ERROR$MSG('INVALID IDENTIFIER')
        CALL ERROR$STOP('CLASSICAL$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$GETCHA(ID_,LENG_,VAL_)
!     **************************************************************************
!     **  CLASSICAL$SET                                                       **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      CHARACTER(*),INTENT(OUT):: VAL_(LENG_)
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO CLASSICAL OBJECT SELECTED')
        CALL ERROR$STOP('CLASSICAL$SETCHA')
      END IF
!
!     ==========================================================================
!     == TYPE = ATOM TYPE                                                     ==
!     ==========================================================================
      IF(ID_.EQ.'TYPE') THEN
        IF(MD%NAT.EQ.0) MD%NAT=LENG_
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$STOP('CLASSICAL$GETCHA')
        END IF
        IF(.NOT.ASSOCIATED(MD%TYPE))ALLOCATE(MD%TYPE(MD%NAT))
        VAL_=MD%TYPE
      ELSE IF(ID_.EQ.'ATOMNAME') THEN
        IF(MD%NAT.EQ.0) MD%NAT=LENG_
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$STOP('CLASSICAL$GETCHA')
        END IF
        IF(.NOT.ASSOCIATED(MD%ATNAME))ALLOCATE(MD%ATNAME(MD%NAT))
        VAL_=MD%ATNAME
      ELSE IF(ID_.EQ.'ELEMENT') THEN
        IF(MD%NAT.EQ.0) MD%NAT=LENG_
        IF(LENG_.NE.MD%NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$STOP('CLASSICAL$GETCHA')
        END IF
        IF(.NOT.ASSOCIATED(MD%ELEMENT))ALLOCATE(MD%ELEMENT(MD%NAT))
        VAL_=MD%ELEMENT
      ELSE
        CALL ERROR$MSG('INVALID IDENTIFIER')
        CALL ERROR$STOP('CLASSICAL$GETCHA')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_INITIALIZE
!     **************************************************************************
!     **  CLASSICAL$GET                                                       **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      LOGICAL(4) :: TCHK
      INTEGER(4) :: NANGLEX
      INTEGER(4) :: NTORSIONX
      INTEGER(4) :: NINVERSIONX
      INTEGER(4) :: NTYPEX
      INTEGER(4) :: NPOTX
      INTEGER(4),ALLOCATABLE :: NONBOND(:,:)
      TYPE(ANGLE_TYPE)    ,ALLOCATABLE :: ANGLE(:)
      TYPE(TORSION_TYPE)  ,ALLOCATABLE :: TORSION(:)
      TYPE(INVERSION_TYPE),ALLOCATABLE :: INVERSION(:)
      INTEGER(4)             :: ISVAR,I
      REAL(8)   ,ALLOCATABLE :: RPOS(:,:)
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO CLASSICAL OBJECT SELECTED')
        CALL ERROR$STOP('CLASSICAL_INITIALIZE')
      END IF
      IF(MD%TINI) THEN
        RETURN
      END IF
      MD%TINI=.TRUE.
                               CALL TRACE$PUSH('CLASSICAL_INITIALIZE')
!
!     ==========================================================================
!     == CHECK WHETHER ALL REQUIRED DATA HAVE BEEN SET                        ==
!     ==========================================================================
      TCHK=.TRUE.
      TCHK=TCHK.AND.(MD%NAT.GT.0)
      TCHK=TCHK.AND.(ASSOCIATED(MD%R0))
      TCHK=TCHK.AND.(ASSOCIATED(MD%RMASS))
      TCHK=TCHK.AND.(ASSOCIATED(MD%TYPE))
      TCHK=TCHK.AND.(ASSOCIATED(MD%QEL))
      TCHK=TCHK.AND.(MD%NBOND.GT.0)
      TCHK=TCHK.AND.(ASSOCIATED(MD%BOND))
      TCHK=TCHK.AND.(ASSOCIATED(MD%BO))
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('NOT ALL MANDATORY DATA ARE SET')
        CALL ERROR$I4VAL('NAT',MD%NAT)
        CALL ERROR$I4VAL('NBOND',MD%NBOND)
        CALL ERROR$L4VAL('R(0)',ASSOCIATED(MD%R0))
        CALL ERROR$L4VAL('MASS',ASSOCIATED(MD%RMASS))
        CALL ERROR$L4VAL('FFTYPE',ASSOCIATED(MD%TYPE))
        CALL ERROR$L4VAL('QEL',ASSOCIATED(MD%QEL))
        CALL ERROR$L4VAL('BOND',ASSOCIATED(MD%BOND))
        CALL ERROR$L4VAL('BOND ORDER',ASSOCIATED(MD%BO))
        CALL ERROR$STOP('CLASSICAL_INITIALIZE')
      END IF
!
!     ==========================================================================
!     == ALLOCATE ARRAYS WITH OBVIOUS DEFAULT VALUES THAT MAY BE SET          ==
!     ==========================================================================
      IF(.NOT.ASSOCIATED(MD%TFREEZE)) THEN 
        ALLOCATE(MD%TFREEZE(MD%NAT))
        MD%TFREEZE=.FALSE.
      END IF
!
      IF(.NOT.ASSOCIATED(MD%RP)) THEN
        ALLOCATE(MD%RP(3,MD%NAT))
        MD%RP=MD%R0
      END IF
!
      IF(.NOT.ASSOCIATED(MD%FORCE)) THEN
        ALLOCATE(MD%FORCE(3,MD%NAT))
        MD%FORCE=0.D0
      END IF
!
      IF(.NOT.ASSOCIATED(MD%VEL)) THEN
        ALLOCATE(MD%VEL(MD%NAT))
        MD%VEL=0.D0
      END IF
!
!     ==========================================================================
!     ==  DEFINE POTENTIALS FOR                                               ==
!     ==========================================================================
      NANGLEX=12*MD%NAT   !THE VALUE OF 12 IS PROBABLY NOT THE OPTIMUM CHOICE
      NTORSIONX=9*MD%NBOND
      NINVERSIONX=MD%NAT*3
      NTYPEX=MD%NAT
      NPOTX=500
      ALLOCATE(ANGLE(NANGLEX))
      ALLOCATE(TORSION(NTORSIONX))
      ALLOCATE(INVERSION(NINVERSIONX))
      CALL CLASSICAL_ZERO_ANGLE(NANGLEX,ANGLE)
      CALL CLASSICAL_ZERO_TORSION(NTORSIONX,TORSION)
      CALL CLASSICAL_ZERO_INVERSION(NINVERSIONX,INVERSION)
!
      ALLOCATE(NONBOND(NTYPEX,NTYPEX))
      ALLOCATE(MD%POT(NPOTX))
      ALLOCATE(MD%ITYPE(MD%NAT))
      CALL CLASSICAL_FORCEFIELDSETUP(MD%NBOND,MD%BOND,MD%BO &
     &             ,MD%NAT,MD%TYPE,MD%ITYPE &
     &             ,NANGLEX,MD%NANGLE,ANGLE &
     &             ,NTORSIONX,MD%NTORSION,TORSION &
     &             ,NINVERSIONX,MD%NINVERSION,INVERSION &
     &             ,NTYPEX,MD%NTYPE,NONBOND &
     &             ,NPOTX,MD%NPOT,MD%POT)
      ALLOCATE(MD%NONBOND(MD%NTYPE,MD%NTYPE))
      MD%NONBOND(:,:)=NONBOND(1:MD%NTYPE,1:MD%NTYPE)
      DEALLOCATE(NONBOND)
      ALLOCATE(MD%ANGLE(MD%NANGLE))
      ALLOCATE(MD%TORSION(MD%NTORSION))
      ALLOCATE(MD%INVERSION(MD%NINVERSION))
      MD%ANGLE=ANGLE(1:MD%NANGLE)
      MD%TORSION=TORSION(1:MD%NTORSION)
      MD%INVERSION=INVERSION(1:MD%NINVERSION)
      DEALLOCATE(ANGLE)
      DEALLOCATE(TORSION)
      DEALLOCATE(INVERSION)
!
!     ==========================================================================
!     ==  LATTICE VECTORS                                                     ==
!     ==========================================================================
      IF(.NOT.ASSOCIATED(MD%RBAS0)) THEN
        ALLOCATE(MD%RBAS0(3,3))
        MD%RBAS0(:,:)=0.D0
        DO I=1,3
          MD%RBAS0(I,I)=10000.D0
        ENDDO
        MD%MAXDIV=0
      ELSE
        MD%MAXDIV=1
      END IF
!   
!     ==========================================================================
!     ==  CALCULATE EXCLUSIONS                                                ==
!     ==========================================================================
      CALL TIP4P$SELECT(MD%MDNAME)
      ALLOCATE(MD%I2FIRST(MD%NAT))
!     == FIRST COUNT THE NUMBER OF EXCLUSIONS ....
      MD%NEXCL=1
      ALLOCATE(MD%EXCLUSION(2,MD%NEXCL))
      ISVAR=MD%NEXCL
      CALL CLASSICAL_EXCLUSIONS(MD%NAT,MD%NBOND,MD%BOND,MD%NANGLE,MD%ANGLE &
     &               ,MD%NTORSION,MD%TORSION,MD%MAXDIV,ISVAR,MD%NEXCL &
     &               ,MD%EXCLUSION,MD%I2FIRST)
      DEALLOCATE(MD%EXCLUSION)
!     == NOW CALCULATE EXCLUSIONS
      ALLOCATE(MD%EXCLUSION(2,MD%NEXCL))
      CALL CLASSICAL_EXCLUSIONS(MD%NAT,MD%NBOND,MD%BOND,MD%NANGLE,MD%ANGLE &
     &               ,MD%NTORSION,MD%TORSION,MD%MAXDIV,MD%NEXCL,ISVAR &
     &               ,MD%EXCLUSION,MD%I2FIRST)
      CALL TIP4P$SELECT('NONE')
!
!     ==========================================================================
!     == CALCULATE NEIGHBORLIST                                               ==
!     ==========================================================================
      MD%NNBX=1
      MD%NNB=1
      ALLOCATE(MD%IBLIST(MD%NNBX))
      DO I=1,MD%NNBX
        NULLIFY(MD%IBLIST(I)%IT)
      ENDDO    
      CALL CLASSICAL_NEIGHBORS(MD%NAT,MD%R0,MD%RBAS0,MD%NNBX,MD%NNB,MD%IBLIST &
     &                        ,MD%MAXDIV,MD%NEXCL,MD%EXCLUSION,MD%I2FIRST)
      DO I=1,MD%NNBX
        IF(ASSOCIATED(MD%IBLIST(I)%IT))DEALLOCATE(MD%IBLIST(I)%IT)
      ENDDO    
      DEALLOCATE(MD%IBLIST)
      MD%NNBX=MD%NNB
      ALLOCATE(MD%IBLIST(MD%NNBX))
      DO I=1,MD%NNBX
        NULLIFY(MD%IBLIST(I)%IT)
      ENDDO      
      CALL CLASSICAL_NEIGHBORS(MD%NAT,MD%R0,MD%RBAS0,MD%NNBX,MD%NNB,MD%IBLIST &
     &                      ,MD%MAXDIV,MD%NEXCL,MD%EXCLUSION,MD%I2FIRST)
!
!     ==========================================================================
!     == ENFORCE TIP4P CONSTRAINTS                                            ==
!     ==========================================================================
      ALLOCATE(RPOS(3,MD%NAT))
      CALL TIP4P$SELECT(MD%MDNAME)
      RPOS=MD%R0
      CALL TIP4P$CONSTRAINTS(MD%NAT,MD%RMASS,RPOS,MD%R0)
      RPOS=MD%RM
      CALL TIP4P$CONSTRAINTS(MD%NAT,MD%RMASS,RPOS,MD%RM)
      CALL TIP4P$SELECT('NONE')
                               CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$NEIGHBORS
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      INTEGER(4) :: I
!     **************************************************************************
                               CALL TRACE$PUSH('CLASSICAL$NEIGHBORS')
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('CLASSICAL OBJECT NOT SELECTED')
        CALL ERROR$STOP('CLASSICAL$NEIGHBORS')
      END IF
      CALL CLASSICAL_INITIALIZE
!
!     ==========================================================================
!     ==  NOW CALCULATE NEIGBORLIST                                           ==
!     ==========================================================================
      CALL CLASSICAL_NEIGHBORS(MD%NAT,MD%R0,MD%RBAS0,MD%NNBX,MD%NNB,MD%IBLIST &
     &                      ,MD%MAXDIV,MD%NEXCL,MD%EXCLUSION,MD%I2FIRST)
      IF(MD%NNB.GT.MD%NNBX) THEN
        DO I=1,MD%NNBX
          IF(ASSOCIATED(MD%IBLIST(I)%IT))DEALLOCATE(MD%IBLIST(I)%IT)
        ENDDO    
        MD%NNBX=MD%NNB
        DEALLOCATE(MD%IBLIST)
        ALLOCATE(MD%IBLIST(MD%NNBX))
        DO I=1,MD%NNBX
          NULLIFY(MD%IBLIST(I)%IT)
        ENDDO      
        CALL CLASSICAL_NEIGHBORS(MD%NAT,MD%R0,MD%RBAS0,MD%NNBX,MD%NNB &
     &                    ,MD%IBLIST,MD%MAXDIV,MD%NEXCL,MD%EXCLUSION,MD%I2FIRST)
      END IF
                               CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$ETOT(EPOT)  
!     **************************************************************************
!     **  CALCULATES THE TOTAL ENERGY FOR THE CURRENT STRUCTURE               **
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: EPOT
      INTEGER(4)             :: IAT,I
!     **************************************************************************
                               CALL TRACE$PUSH('CLASSICAL$ETOT')
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('CLASSICAL OBJECT NOT SELECTED')
        CALL ERROR$STOP('CLASSICAL$ETOT')
      END IF
      CALL CLASSICAL_INITIALIZE
!
!     ==========================================================================
!     ==  TAKE CARE OF FORCES ACTING ON DUMMY ATOMS                           ==
!     ==========================================================================
      CALL CLASSICAL_DUMMY_POSITION(MD%NAT,MD%TYPE,MD%R0,MD%NBOND,MD%BOND)
      CALL TIP4P$SELECT(MD%MDNAME)  ! IS USED IN CLASSICAL_DUMMY_POSITION
      CALL TIP4P$DUMMY_POSITIONS(MD%NAT,MD%R0)
!
!     ==========================================================================
!     ==  NOW CALCULATE TOTAL ENERGY AND FORCES                               ==
!     ==========================================================================
      EPOT=0.D0
      CALL CLASSICAL_ETOTAL(MD%NAT,MD%R0,MD%QEL,MD%ITYPE,EPOT &
     &             ,MD%FORCE,MD%VEL,MD%RBAS0,MD%SIGMA &
     &             ,MD%TLONGRANGE &
     &             ,MD%NBOND,MD%BOND,MD%NANGLE,MD%ANGLE &
     &             ,MD%NTORSION,MD%TORSION &
     &             ,MD%NINVERSION,MD%INVERSION,MD%NTYPE,MD%NONBOND &
     &             ,MD%NNB,MD%IBLIST,MD%NPOT,MD%POT)
!
!     ==========================================================================
!     ==  TAKE CARE OF FORCES ACTING ON DUMMY ATOMS                           ==
!     ==========================================================================
      CALL CLASSICAL_DUMMY_FORCE(MD%NAT,MD%TYPE,MD%FORCE,MD%NBOND,MD%BOND)
      CALL TIP4P$DUMMY_FORCES(MD%NAT,MD%FORCE)
      CALL TIP4P$SELECT('NONE')
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
!
!     CALL CLASSICAL_REFLECTINGBOUNDARY(NAT,R0,FORCE,RCENT,FCENT
!    &                  ,RAD,ALPHA)
      DO IAT=1,MD%NAT
        IF(MD%TFREEZE(IAT)) THEN
          DO I=1,3
            MD%FORCE(I,IAT)=0.D0
          ENDDO
        END IF
      ENDDO
                               CALL TRACE$POP
      RETURN
    END SUBROUTINE CLASSICAL$ETOT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_DUMMY_POSITION(NAT,TYPE,R,NBOND,BOND)
!     **************************************************************************
!     ** CALCULATE POSITIONS OF DUMMY ATOM                                    **
!     **                                                                      **
!     **  DUMMY ATOMS ARE INTRODUCED TO DEFINE PI BONDED SYSTEMS              **
!     **  SUCH AS CYCLOPENTADIENYL (CPR, CPR_B), ALLYL (CIR) AND              **
!     **  AN OLEFIN (PIR).                                                    **
!     **                                                                      **
!     **  THEIR POSITIONS ARE DEFINED AS THE AVERAGE POSITION                 **
!     **  OF C_2 OR C_R ATOMS BONDED TO THEM                                  **
!     **                                                                      **
!     **  THEIR MASSES MUST BE SMALL BUT NONZERO                              **
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE, ONLY : BOND_TYPE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)   :: NAT
      CHARACTER(*),INTENT(IN)   :: TYPE(NAT)
      REAL(8)     ,INTENT(INOUT):: R(3,NAT)
      INTEGER(4)  ,INTENT(IN)   :: NBOND
      TYPE(BOND_TYPE),INTENT(IN):: BOND(NBOND)
      INTEGER(4)                :: IAT,IBOND
      INTEGER(4)                :: NN,IAT1,IAT2
      REAL(8)                   :: R0(3)
!     **************************************************************************
!
!     ==========================================================================
!     == SCAN ATOMS FOR DUMMY ATOMS AND RECALCULATE THEIR POSITIONS           ==
!     ==========================================================================
      DO IAT=1,NAT
!
!       ========================================================================
!       == MIDDLE ATOM OF AROMATIC SYSTEMS, ALLYL AND OLEFIN                  ==
!       ========================================================================
        IF(TYPE(IAT).EQ.'CPR'.OR.TYPE(IAT).EQ.'CPR_B'.OR. &
     &     TYPE(IAT).EQ.'PIR'.OR.TYPE(IAT).EQ.'CIR') THEN
          R0=0.D0
          NN=0
          DO IBOND=1,NBOND
            IAT1=BOND(IBOND)%IAT1
            IAT2=BOND(IBOND)%IAT2
            IF(IAT1.EQ.IAT) THEN
              IF(TYPE(IAT2).NE.'C_2'.AND.TYPE(IAT2).NE.'C_R') CYCLE
              NN=NN+1
              R0(:)=R0(:)+R(:,IAT2)
            ELSE IF(IAT2.EQ.IAT) THEN
              IF(TYPE(IAT1).NE.'C_2'.AND.TYPE(IAT1).NE.'C_R') CYCLE
              NN=NN+1
              R0(:)=R0(:)+R(:,IAT1)
            END IF
          ENDDO
          R(:,IAT)=R0(:)/DBLE(NN)
        END IF
      ENDDO  !END OF LOOP OVER ATOMS SCANNING FOR DUMMY ATOMS
      RETURN
    END SUBROUTINE CLASSICAL_DUMMY_POSITION
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_DUMMY_FORCE(NAT,TYPE,F,NBOND,BOND)
!     **************************************************************************
!     ** CALCULATE DUMMY ATOM FORCES                                          **
!     **************************************************************************
      USE CLASSICAL_MODULE, ONLY : BOND_TYPE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)   :: NAT
      CHARACTER(*),INTENT(IN)   :: TYPE(NAT)
      REAL(8)     ,INTENT(INOUT):: F(3,NAT)
      INTEGER(4)  ,INTENT(IN)   :: NBOND
      TYPE(BOND_TYPE),INTENT(IN):: BOND(NBOND)
      INTEGER(4)                :: IAT,IBOND
      INTEGER(4)                :: NN,IAT1,IAT2
      REAL(8)                   :: F0(3)
!     **************************************************************************
!
!     ==========================================================================
!     == SCAN ATOMS FOR DUMMY ATOMS AND MAP THEIR FORCES ONTO REAL ATOMS      ==
!     ==========================================================================
      DO IAT=1,NAT
!
!       ========================================================================
!       == MIDDLE ATOM OF AROMATIC SYSTEMS, ALLYL AND OLEFIN                  ==
!       ========================================================================
        IF(TYPE(IAT).EQ.'CPR'.OR.TYPE(IAT).EQ.'CPR_B'.OR. &
     &     TYPE(IAT).EQ.'PIR'.OR.TYPE(IAT).EQ.'CIR') THEN
          NN=0
          DO IBOND=1,NBOND
            IAT1=BOND(IBOND)%IAT1
            IAT2=BOND(IBOND)%IAT2
            IF(IAT1.EQ.IAT) THEN
              IF(TYPE(IAT2).NE.'C_2'.AND.TYPE(IAT2).NE.'C_R') CYCLE
              NN=NN+1
            ELSE IF(IAT2.EQ.IAT) THEN
              IF(TYPE(IAT1).NE.'C_2'.AND.TYPE(IAT1).NE.'C_R') CYCLE
              NN=NN+1
            END IF
          ENDDO
          F0(:)=F(:,IAT)/DBLE(NN)
          DO IBOND=1,NBOND
            IAT1=BOND(IBOND)%IAT1
            IAT2=BOND(IBOND)%IAT2
            IF(IAT1.EQ.IAT) THEN
              IF(TYPE(IAT2).NE.'C_2'.AND.TYPE(IAT2).NE.'C_R') CYCLE
              F(:,IAT2)=F(:,IAT2)+F0(:)
            ELSE IF(IAT2.EQ.IAT) THEN
              IF(TYPE(IAT1).NE.'C_2'.AND.TYPE(IAT1).NE.'C_R') CYCLE
              F(:,IAT1)=F(:,IAT1)+F0(:)
            END IF
          ENDDO
          F(:,IAT)=0.D0   !SET FORCE OF DUMMY ATOM TO ZERO
        END IF
      ENDDO !END OF LOOP OVER ATOMS SCANNING FOR DUMMY ATOMS
      RETURN
    END SUBROUTINE CLASSICAL_DUMMY_FORCE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$STOP
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
!     **************************************************************************
                               CALL TRACE$PUSH('CLASSICAL$ETOT')
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('CLASSICAL OBJECT NOT SELECTED')
        CALL ERROR$STOP('CLASSICAL$ETOT')
      END IF
      CALL CLASSICAL_INITIALIZE
!
!     ==========================================================================
!     ==  NOW SET VELOCITIES TO ZERO                                          ==
!     ==========================================================================
      MD%RM=MD%R0
                               CALL TRACE$POP
      RETURN
    END SUBROUTINE CLASSICAL$STOP
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$PROPAGATE(DELT,ANNE)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE, ONLY : SELECTED,MD
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: DELT
      REAL(8)   ,INTENT(IN)   :: ANNE
      INTEGER(4)              :: IAT,I,NATM,NH2O,N
      REAL(8)                 :: SVAR1,SVAR2,SVAR3
      INTEGER(4)  ,ALLOCATABLE::NUMBER(:)
      CHARACTER(4),ALLOCATABLE::NAME(:),ATOM(:)
      REAL(8)     ,ALLOCATABLE::H2O(:,:),H2OO(:,:)
!     **************************************************************************
                               CALL TRACE$PUSH('CLASSICAL$PROPAGATE')
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('CLASSICAL OBJECT NOT SELECTED')
        CALL ERROR$STOP('CLASSICAL$ETOT')
      END IF
      CALL CLASSICAL_INITIALIZE
!
!     ==========================================================================
!     ==  PROPAGATE ATOMS WITHOUT CONSTRAINTS                                 ==
!     ==========================================================================
      DO IAT=1,MD%NAT
!IF (IAT.EQ.1) PRINT*,"FLAG: FORCE1 ",IAT,MD%FORCE(:,IAT)

        SVAR1=2.D0/(1.D0+ANNE)
        SVAR2=1.D0-SVAR1
        SVAR3=DELT**2/MD%RMASS(IAT)/(1.D0+ANNE)
        DO I=1,3
          MD%RP(I,IAT)=MD%R0(I,IAT)*SVAR1+MD%RM(I,IAT)*SVAR2 &
     &                             +MD%FORCE(I,IAT)*SVAR3
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == APPLY CONSTRAINTS                                                    ==
!     ==========================================================================
      DO IAT=1,MD%NAT
        IF(MD%TFREEZE(IAT)) THEN
          DO I=1,3
            MD%RP(I,IAT)=MD%R0(I,IAT)
          ENDDO
        END IF
      ENDDO
!
      CALL TIP4P$SELECT(MD%MDNAME)
      CALL TIP4P$CONSTRAINTS(MD%NAT,MD%RMASS,MD%R0,MD%RP)
      CALL TIP4P$SELECT('NONE')
!!$PRINT*,'  '
!!$DO IAT=1,MD%NAT
!!$PRINT*,'MD%DR(3) ',MD%RP(:,IAT)-MD%R0(:,IAT)
!!$ENDDO
!
!     =========================================================================
!     == TIP3P NEW                                                           ==
!     =========================================================================
IF(TRIM(ADJUSTL(MD%FF)).EQ.'AMBER') THEN
       CALL FORCEFIELD$GETI4('NMMATOM',NATM)
       ALLOCATE(NUMBER(NATM))
       ALLOCATE(NAME(NATM))
       ALLOCATE(ATOM(NATM))
       CALL FORCEFIELD$GETI4A('RESSEQ',NATM,NUMBER)
       CALL FORCEFIELD$GETCHA('RESNAME',NATM,NAME)
       CALL FORCEFIELD$GETCHA('NAME',NATM,ATOM)
       
       N=0
       NH2O=0
       DO I=1,NATM
          IF (NAME(I) .EQ. 'WA2') THEN
             N=N+1
             IF (N .EQ. 3) THEN
                NH2O=NH2O+1
                N=0
             END IF
          END IF
       END DO

       ALLOCATE(H2O(3,3))
       H2O(:,:)=0.D0
       ALLOCATE(H2OO(3,3))
       H2OO(:,:)=0.D0
!      =====SEND SINGLE H2O TO TIP3P=====================================
       I=1
       DO WHILE (I .LT. NATM)
          IF (.NOT. MD%TFREEZE(I)) THEN
             IF ((NUMBER(I) .EQ. NUMBER(I+1)) .AND. (NUMBER(I+1) .EQ. NUMBER(I+2))) THEN
                IF (ATOM(I) .EQ. ' OH2') THEN   ! 1ST ATOM O
                    H2O(1,:)=MD%RP(:,I) ! 1ST ATOM O
                    H2O(2,:)=MD%RP(:,I+1) ! 2ND ATOM H
                    H2O(3,:)=MD%RP(:,I+2) ! 3RD ATOM H
                    H2OO(1,:)=MD%R0(:,I) ! 1ST ATOM O
                    H2OO(2,:)=MD%R0(:,I+1) ! 2ND ATOM H
                    H2OO(3,:)=MD%R0(:,I+2) ! 3RD ATOM H
!PRINT*,"FLAG: I=",I
                    CALL FORCEFIELD$TIP3PNEW(H2O,H2OO,10)

                    MD%RP(:,I)=H2O(1,:)     ! NEW POS. 1ST ATOM
                    MD%RP(:,I+1)=H2O(2,:)   ! NEW POS. 2ND ATOM
                    MD%RP(:,I+2)=H2O(3,:)   ! NEW POS. 3RD ATOM

                ELSE IF (ATOM(I+1) .EQ. ' OH2') THEN ! 2ND ATOM O
                    H2O(1,:)=MD%RP(:,I+1)! 1ST ATOM O
                    H2O(2,:)=MD%RP(:,I) ! 2ND ATOM H
                    H2O(3,:)=MD%RP(:,I+2) ! 3RD ATOM H

                    H2OO(1,:)=MD%R0(:,I+1)! 1ST ATOM O
                    H2OO(2,:)=MD%R0(:,I) ! 2ND ATOM H
                    H2OO(3,:)=MD%R0(:,I+2) ! 3RD ATOM H

!PRINT*,"FLAG: I=",I
                    CALL FORCEFIELD$TIP3PNEW(H2O,H2OO,10)

                    MD%RP(:,I+1)=H2O(1,:)     ! NEW POS. 1ST ATOM
                    MD%RP(:,I)=H2O(2,:)   ! NEW POS. 2ND ATOM
                    MD%RP(:,I+2)=H2O(3,:)   ! NEW POS. 3RD ATOM



                ELSE IF (ATOM(I+2) .EQ. ' OH2') THEN ! 3RD ATOM O 

                    H2O(1,:)=MD%RP(:,I+2) ! 1ST ATOM O
                    H2O(2,:)=MD%RP(:,I) ! 2ND ATOM H
                    H2O(3,:)=MD%RP(:,I+1) ! 3RD ATOM H

                    H2OO(1,:)=MD%R0(:,I+2) ! 1ST ATOM O
                    H2OO(2,:)=MD%R0(:,I) ! 2ND ATOM H
                    H2OO(3,:)=MD%R0(:,I+1) ! 3RD ATOM H

!PRINT*,"FLAG: I=",I
                    CALL FORCEFIELD$TIP3PNEW(H2O,H2OO,10)

                    MD%RP(:,I+2)=H2O(1,:)     ! NEW POS. 1ST ATOM
                    MD%RP(:,I)=H2O(2,:)   ! NEW POS. 2ND ATOM
                    MD%RP(:,I+1)=H2O(3,:)   ! NEW POS. 3RD ATOM
                END IF
            
! PRINT *,'***********TIP3P********************'
! PRINT *,'I',I
! PRINT *,'NATM',NATM
! PRINT *,'NH2O',NH2O
! PRINT *,'NUMBER:',NUMBER(:)
! PRINT *,'NAME: ',NAME(:)
! PRINT *,'ATOM: ',ATOM(:)
! PRINT *,'H2O(:)',H2O(:,:)
! PRINT *,'***********TIP3P********************'
                I=I+3
              ELSE
                I=I+1
              END IF
          ELSE
             I=I+1
          END IF   

       END DO
    END IF


!STOP



                               CALL TRACE$POP
      RETURN
    END SUBROUTINE CLASSICAL$PROPAGATE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$EKIN(DELT_,EKIN_)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: DELT_
      REAL(8)   ,INTENT(OUT):: EKIN_
      INTEGER(4)            :: IAT,I
      REAL(8)               :: SVAR
!     **************************************************************************
                               CALL TRACE$PUSH('CLASSICAL$EKIN')
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('CLASSICAL OBJECT NOT SELECTED')
        CALL ERROR$STOP('CLASSICAL$EKIN')
      END IF
      CALL CLASSICAL_INITIALIZE
!
!     ==========================================================================
!     ==  NOW CALCULATE TOTAL ENERGY AND FORCES                               ==
!     ==========================================================================
      EKIN_=0.D0
      DO IAT=1,MD%NAT
        SVAR=0.5D0*MD%RMASS(IAT)/(2.D0*DELT_)**2
        DO I=1,3
          EKIN_=EKIN_+SVAR*(MD%RP(I,IAT)-MD%RM(I,IAT))**2
        ENDDO
      ENDDO
                               CALL TRACE$POP
      RETURN
    END SUBROUTINE CLASSICAL$EKIN
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$SWITCH
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: IAT,I
      LOGICAL(4)             :: TCHK
!     **************************************************************************
                             CALL TRACE$PUSH('CLASSICAL$SWITCH')           
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('CLASSICAL OBJECT NOT SELECTED')
        CALL ERROR$STOP('CLASSICAL$ETOT')
      END IF
      CALL CLASSICAL_INITIALIZE
!
!     ==========================================================================
!     ==  NOW CALCULATE TOTAL ENERGY AND FORCES                               ==
!     ==========================================================================
   IF(MD%MDNAME.EQ.'SHADOW') THEN
!**************** BEGIN OLD *************************** 
      DO IAT=1,MD%NAT
        DO I=1,3
          MD%RM(I,IAT)=MD%R0(I,IAT)
          MD%R0(I,IAT)=MD%RP(I,IAT)
          MD%RP(I,IAT)=0.D0
          MD%FORCE(I,IAT)=0.D0
        ENDDO
      ENDDO
!**************** END OLD *****************************
  ELSE
      TCHK=MD%TMMSPHERE
      IF(TCHK) THEN
         CALL FORCEFIELD$MMSPHERE
      ELSE
         DO IAT=1,MD%NAT
            DO I=1,3
               MD%RM(I,IAT)=MD%R0(I,IAT)
               MD%R0(I,IAT)=MD%RP(I,IAT)
               MD%RP(I,IAT)=0.D0
               MD%FORCE(I,IAT)=0.D0
            ENDDO
         ENDDO
      END IF
   END IF
!
      RETURN
      END SUBROUTINE CLASSICAL$SWITCH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$MAXFORCE(FMAX_)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: FMAX_
      INTEGER(4)             :: IAT,I
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('CLASSICAL OBJECT NOT SELECTED')
        CALL ERROR$STOP('CLASSICAL$ETOT')
      END IF
      CALL CLASSICAL_INITIALIZE
!
!     ==========================================================================
!     ==  NOW CALCULATE TOTAL ENERGY AND FORCES                               ==
!     ==========================================================================
      FMAX_=0.D0
      DO IAT=1,MD%NAT
        IF(MD%TFREEZE(IAT)) CYCLE ! DO NOT INCLUDE THE FORCE ON FROZEN ATOMS
        DO I=1,3
          FMAX_=FMAX_+MD%FORCE(I,IAT)**2
        ENDDO
      ENDDO
      FMAX_=SQRT(FMAX_)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL$REPORT(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      REAL(8)                :: PRTONM
      REAL(8)                :: ANGSTROM
      INTEGER(4)             :: IAT,I,IB,IA,IPOT
      INTEGER(4)             :: IAT1,IAT2,IAT3
      REAL(8)                :: D,A
      REAL(8)                :: DX1,DY1,DZ1
      REAL(8)                :: DX2,DY2,DZ2
      REAL(8)                :: D11,D22,D12
      REAL(8)                :: DR(3),DR1(3),DR2(3)
      INTEGER(4)             :: NTASKS,THISTASK
!     **************************************************************************
                               CALL TRACE$PUSH('CLASSICAL$REPORT')
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('CLASSICAL OBJECT NOT SELECTED')
        CALL ERROR$STOP('CLASSICAL$ETOT')
      END IF
      CALL CLASSICAL_INITIALIZE
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) THEN
        CALL TRACE$POP
        RETURN
      END IF
!
!     ==========================================================================
!     ==  REPORT ATOMIC POSITIONS                                             ==
!     ==========================================================================
      WRITE(NFIL,FMT='(''=== POSITIONS IN ANGSTROM ===; NAT='',I10)')MD%NAT
      IF(LOD.GE.2) THEN
        CALL CONSTANTS('U',PRTONM)
        CALL CONSTANTS('ANGSTROM',ANGSTROM)
        DO IAT=1,MD%NAT
          WRITE(NFIL,FMT='(I3,1X,A,"R=",3F10.5,2X' &
     &                 //',"A ;M=",F10.5,"U"," Q[E] ",F10.5)') &
     &                   IAT,MD%TYPE(IAT),(MD%R0(I,IAT)/ANGSTROM,I=1,3) &
     &                                   ,MD%RMASS(IAT)/PRTONM,-MD%QEL(IAT)
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  REPORT BONDS                                                        ==
!     ==========================================================================
      WRITE(NFIL,FMT='(''===  BONDS IN ANGSTROM  ===; NBOND='',I10)')MD%NBOND
      IF(LOD.GE.3) THEN
        CALL CONSTANTS('ANGSTROM',ANGSTROM)
        DO IB=1,MD%NBOND
          IAT1=MD%BOND(IB)%IAT1
          IAT2=MD%BOND(IB)%IAT2
          DR(:)=MD%R0(:,IAT1)-MD%R0(:,IAT2)
          IF(ASSOCIATED(MD%BOND(IB)%IT2)) THEN
            DR(:)=DR(:)-MD%RBAS0(:,1)*REAL(MD%BOND(IB)%IT2(1),KIND=8)
            DR(:)=DR(:)-MD%RBAS0(:,2)*REAL(MD%BOND(IB)%IT2(2),KIND=8)
            DR(:)=DR(:)-MD%RBAS0(:,3)*REAL(MD%BOND(IB)%IT2(3),KIND=8)
            D=SQRT(SUM(DR(:)**2))
            WRITE(NFIL,FMT='(2I5,F10.5," A; IT=",3I5)') &
     &                                      IAT1,IAT2,D/ANGSTROM,MD%BOND(IB)%IT2
          ELSE
            D=SQRT(SUM(DR(:)**2))
            WRITE(NFIL,FMT='(2I5,F10.5," A")')IAT1,IAT2,D/ANGSTROM
          END IF
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  REPORT BOND ANGLES                                                  ==
!     ==========================================================================
      WRITE(NFIL,FMT='(''===  ANGLES IN DEGREE  ===; NANGLE='',I10)')MD%NANGLE
      IF(LOD.GE.3) THEN
         DO IA=1,MD%NANGLE
            IAT1=MD%ANGLE(IA)%IAT1
            IAT2=MD%ANGLE(IA)%IAT2
            IAT3=MD%ANGLE(IA)%IAT3
            DR1(:)=MD%R0(:,IAT1)-MD%R0(:,IAT2)
            IF(ASSOCIATED(MD%ANGLE(IA)%IT1)) THEN
              DR1(:)=DR1(:)+MD%RBAS0(:,1)*REAL(MD%ANGLE(IA)%IT1(1),KIND=8)
              DR1(:)=DR1(:)+MD%RBAS0(:,2)*REAL(MD%ANGLE(IA)%IT1(2),KIND=8)
              DR1(:)=DR1(:)+MD%RBAS0(:,3)*REAL(MD%ANGLE(IA)%IT1(3),KIND=8)
            END IF
            DR2(:)=MD%R0(:,IAT3)-MD%R0(:,IAT2)
            IF(ASSOCIATED(MD%ANGLE(IA)%IT3)) THEN
              DR2(:)=DR2(:)+MD%RBAS0(:,1)*REAL(MD%ANGLE(IA)%IT3(1),KIND=8)
              DR2(:)=DR2(:)+MD%RBAS0(:,2)*REAL(MD%ANGLE(IA)%IT3(2),KIND=8)
              DR2(:)=DR2(:)+MD%RBAS0(:,3)*REAL(MD%ANGLE(IA)%IT3(3),KIND=8)
            END IF
            DX1=DR1(1)
            DY1=DR1(2)
            DZ1=DR1(3)
            DX2=DR2(1)
            DY2=DR2(2)
            DZ2=DR2(3)
            D11=DX1*DX1+DY1*DY1+DZ1*DZ1
            D12=DX1*DX2+DY1*DY2+DZ1*DZ2
            D22=DX2*DX2+DY2*DY2+DZ2*DZ2
            A=D12/SQRT(D11*D22)
            A=DACOS(A)/(2.D0*PI)*360.D0
            WRITE(NFIL,FMT='(3I5,F10.5)')IAT1,IAT2,IAT3,A
         ENDDO
      END IF
!
!     ==========================================================================
!     ==  REPORT BOND POTENTIALS                                              ==
!     ==========================================================================
      WRITE(NFIL,FMT='(''===  POTENTIALS  ===; NPOT='',I10)')MD%NPOT
      IF(LOD.GE.4) THEN
        DO IPOT=1,MD%NPOT
          WRITE(NFIL,FMT='(I5,1X,A)')IPOT,MD%POT(IPOT)%ID
!         WRITE(NFIL,FMT='(10F8.3)')(MD%POT(IPOT)%VAL(I),I=1,MD%POT(IPOT)%NX,10)
!         WRITE(NFIL,FMT='(10F8.3)')(MD%POT(IPOT)%DER(I),I=1,MD%POT(IPOT)%NX,10)
        ENDDO
      END IF
      CALL TRACE$POP
      RETURN
    END SUBROUTINE CLASSICAL$REPORT
! 
!     ==========================================================================
!     ==========================================================================
!     ====                                                                  ====
!     ====  ROUTINES USED BY CLASSICAL                                      ====
!     ====                                                                  ====
!     ==========================================================================
!     ==========================================================================
!     ==  CLASSICAL_ETOTAL                                                  ====
!     ==    CLASSICAL_ECOULOMB                                              ====
!     ==    CLASSICAL_EBOND                                                 ====
!     ==    CLASSICAL_EANGLE                                                ====
!     ==    CLASSICAL_ETORSION                                              ====
!     ==    CLASSICAL_EINVERSION                                            ====
!     ==  CLASSICAL_FORCEFIELDSETUP                                         ====
!     ==    CLASSICAL_BONDPOT                                               ====
!     ==    CLASSICAL_ANGLEPOT                                              ====
!     ==    CLASSICAL_TORSIONPOT                                            ====
!     ==    CLASSICAL_INVERSIONPOT                                          ====
!     ==  CLASSICAL_UFFTABLE                                                ====
!     ==  CLASSICAL_NEIGHBORS                                               ====
!     ==  CLASSICAL_EXCLUSIONS                                              ====
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_ETOTAL(NAT,R,Q,ITYPE,ETOT,F,V,RBAS,SIGMA &
     &               ,TLONGRANGE &
     &               ,NBOND,BOND,NANGLE,ANGLE,NTORSION,TORSION &
     &               ,NINVERSION,INVERSION,NTYPE,NONBOND &
     &               ,NNB,NBLIST,NPOT,POT)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE, ONLY : POT_TYPE,NONBOND_TYPE &
     &                ,BOND_TYPE,ANGLE_TYPE,TORSION_TYPE,INVERSION_TYPE
      USE MPE_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TLONGRANGE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      REAL(8)   ,INTENT(IN) :: Q(NAT)
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(OUT):: F(3,NAT)
      REAL(8)   ,INTENT(OUT):: V(NAT)
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)   ,INTENT(OUT):: SIGMA(3,3)
      INTEGER(4),INTENT(IN) :: ITYPE(NAT)
      INTEGER(4),INTENT(IN) :: NBOND
      TYPE(BOND_TYPE),INTENT(IN) :: BOND(NBOND)
      INTEGER(4),INTENT(IN) :: NANGLE
      TYPE(ANGLE_TYPE),INTENT(IN) :: ANGLE(NANGLE)
      INTEGER(4),INTENT(IN) :: NTORSION
      TYPE(TORSION_TYPE),INTENT(IN) :: TORSION(NTORSION)
      INTEGER(4),INTENT(IN) :: NINVERSION
      TYPE(INVERSION_TYPE),INTENT(IN) :: INVERSION(NINVERSION)
      INTEGER(4),INTENT(IN) :: NNB
      TYPE(NONBOND_TYPE),INTENT(IN) :: NBLIST(NNB)
      INTEGER(4),INTENT(IN) :: NTYPE
      INTEGER(4),INTENT(IN) :: NONBOND(NTYPE,NTYPE)
      INTEGER(4),INTENT(IN) :: NPOT
      TYPE(POT_TYPE),INTENT(IN) :: POT(NPOT)  ! POTENTIAL
      INTEGER(4)            :: IAT,I,J,I2BODY,IANGLE,ITORSION,IINV
      INTEGER(4)            :: IAT1,IAT2,IAT3,IAT4,IPOT
      REAL(8)               :: R1(3),R2(3),R3(3),R4(3)
      REAL(8)               :: F1(3),F2(3),F3(3),F4(3)
      REAL(8)               :: DE
      INTEGER(4)            :: THISTASK,NTASKS,ICOUNT
!     **************************************************************************
                            CALL TRACE$PUSH('CLASSICAL_ETOTAL')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
! 
!     ==========================================================================
!     ==  ZERO OUT ARRAYS                                                     ==
!     ==========================================================================
      ETOT=0.D0
      DO IAT=1,NAT
        V(IAT)=0.D0
        DO I=1,3
          F(I,IAT)=0.D0
        ENDDO
      ENDDO
      SIGMA(:,:)=0.D0
! 
!     ==========================================================================
!     ==  COULOMB INTERACTION                                                 ==
!     ==  ATTENTION! RESULT IS DISTRIBUTED OVER TASKS                         ==
!     ==========================================================================
      IF(TLONGRANGE) THEN
        CALL CLASSICAL_ECOULOMB(NAT,R,Q,ETOT,F,V,RBAS,SIGMA &
     &            ,ITYPE,NTYPE,NONBOND,NNB,NBLIST,NPOT,POT)
      END IF
 PRINT*,'COULOMB ',ETOT
! 
!     =========================================================================
!     ==  TWO BODY INTERACTION                                               ==
!     =========================================================================
!     PRINT*,'BOND',NBOND,ETOT
      ICOUNT=0


      DO I2BODY=1,NBOND
        ICOUNT=ICOUNT+1
        IF(MOD(ICOUNT-1,NTASKS).NE.THISTASK-1) CYCLE ! DISTRIBUTE ONTO TASKS
        IAT1=BOND(I2BODY)%IAT1
        IAT2=BOND(I2BODY)%IAT2
        IPOT=BOND(I2BODY)%IPOT
        R1(:)=R(:,IAT1)
        R2(:)=R(:,IAT2)
        IF(ASSOCIATED(BOND(I2BODY)%IT2)) THEN
          R2(:)=R2(:)+RBAS(:,1)*REAL(BOND(I2BODY)%IT2(1),KIND=8) &
    &                +RBAS(:,2)*REAL(BOND(I2BODY)%IT2(2),KIND=8) &
    &                +RBAS(:,3)*REAL(BOND(I2BODY)%IT2(3),KIND=8)
        END IF
        CALL CLASSICAL_EBOND(R1,R2,DE,F1,F2,POT(IPOT))
        ETOT=ETOT+DE
        F(:,IAT1)=F(:,IAT1)+F1(:)
        F(:,IAT2)=F(:,IAT2)+F2(:)
        DO I=1,3
          DO J=1,3
            SIGMA(I,J)=SIGMA(I,J)+R1(I)*F1(J)+R2(I)*F2(J)
          ENDDO
        ENDDO
      ENDDO
! PRINT*,'BOND ',ETOT
! 
!     =========================================================================
!     ==  BOND ANGLE FORCES                                                  ==
!     =========================================================================
!     PRINT*,'ANGLE',NANGLE,ETOT
      DO IANGLE=1,NANGLE
        ICOUNT=ICOUNT+1
        IF(MOD(ICOUNT-1,NTASKS).NE.THISTASK-1) CYCLE
        IAT1=ANGLE(IANGLE)%IAT1
        IAT2=ANGLE(IANGLE)%IAT2
        IAT3=ANGLE(IANGLE)%IAT3
        IPOT=ANGLE(IANGLE)%IPOT
        R1(:)=R(:,IAT1)
        R2(:)=R(:,IAT2)
        R3(:)=R(:,IAT3)
        IF(ASSOCIATED(ANGLE(IANGLE)%IT1)) THEN
          R1(:)=R1(:)+RBAS(:,1)*REAL(ANGLE(IANGLE)%IT1(1),KIND=8) &
    &                +RBAS(:,2)*REAL(ANGLE(IANGLE)%IT1(2),KIND=8) &
    &                +RBAS(:,3)*REAL(ANGLE(IANGLE)%IT1(3),KIND=8)
        END IF
        IF(ASSOCIATED(ANGLE(IANGLE)%IT3)) THEN
          R3(:)=R3(:)+RBAS(:,1)*REAL(ANGLE(IANGLE)%IT3(1),KIND=8) &
    &                +RBAS(:,2)*REAL(ANGLE(IANGLE)%IT3(2),KIND=8) &
    &                +RBAS(:,3)*REAL(ANGLE(IANGLE)%IT3(3),KIND=8)
        END IF
!       PRINT*,'ANGLE IAT1,IAT2,IAT3,IPOT ',IAT1,IAT2,IAT3,IPOT,ETOT
        CALL CLASSICAL_EANGLE(R(1,IAT1),R(1,IAT2),R(1,IAT3)  &
     &                ,DE,F1,F2,F3,POT(IPOT))
        ETOT=ETOT+DE
        F(:,IAT1)=F(:,IAT1)+F1(:)
        F(:,IAT2)=F(:,IAT2)+F2(:)
        F(:,IAT3)=F(:,IAT3)+F3(:)
        DO I=1,3
          DO J=1,3
            SIGMA(I,J)=SIGMA(I,J)+R1(I)*F1(J)+R2(I)*F2(J)+R3(I)*F3(J)
          ENDDO
        ENDDO
      ENDDO 
! PRINT*,'ANGLE ',ETOT
! 
!     =========================================================================
!     ==  TORSION ANGLE FORCES                                               ==
!     =========================================================================
!     PRINT*,'TORSION',NTORSION,ETOT
      DO ITORSION=1,NTORSION
        ICOUNT=ICOUNT+1
        IF(MOD(ICOUNT-1,NTASKS).NE.THISTASK-1) CYCLE
        IAT1=TORSION(ITORSION)%IAT1
        IAT2=TORSION(ITORSION)%IAT2
        IAT3=TORSION(ITORSION)%IAT3
        IAT4=TORSION(ITORSION)%IAT4
        IPOT=TORSION(ITORSION)%IPOT
        R1(:)=R(:,IAT1)
        R2(:)=R(:,IAT2)
        R3(:)=R(:,IAT3)
        R4(:)=R(:,IAT4)
        IF(ASSOCIATED(TORSION(ITORSION)%IT1)) THEN
          R1(:)=R1(:)+RBAS(:,1)*REAL(TORSION(ITORSION)%IT1(1),KIND=8) &
    &                +RBAS(:,2)*REAL(TORSION(ITORSION)%IT1(2),KIND=8) &
    &                +RBAS(:,3)*REAL(TORSION(ITORSION)%IT1(3),KIND=8)
        END IF
        IF(ASSOCIATED(TORSION(ITORSION)%IT3)) THEN
          R3(:)=R3(:)+RBAS(:,1)*REAL(TORSION(ITORSION)%IT3(1),KIND=8) &
    &                +RBAS(:,2)*REAL(TORSION(ITORSION)%IT3(2),KIND=8) &
    &                +RBAS(:,3)*REAL(TORSION(ITORSION)%IT3(3),KIND=8)
        END IF
        IF(ASSOCIATED(TORSION(ITORSION)%IT4)) THEN
          R4(:)=R4(:)+RBAS(:,1)*REAL(TORSION(ITORSION)%IT4(1),KIND=8) &
    &                +RBAS(:,2)*REAL(TORSION(ITORSION)%IT4(2),KIND=8) &
    &                +RBAS(:,3)*REAL(TORSION(ITORSION)%IT4(3),KIND=8)
        END IF
        CALL CLASSICAL_ETORSION(R1,R2,R3,R4,DE,F1,F2,F3,F4,POT(IPOT))
        ETOT=ETOT+DE
!       PRINT*,'TORSION ',IAT1,IAT2,IAT3,IAT4,IPOT,ETOT
        F(:,IAT1)=F(:,IAT1)+F1(:)
        F(:,IAT2)=F(:,IAT2)+F2(:)
        F(:,IAT3)=F(:,IAT3)+F3(:)
        F(:,IAT4)=F(:,IAT4)+F4(:)
        DO I=1,3
          DO J=1,3
            SIGMA(I,J)=SIGMA(I,J)+R1(I)*F1(J)+R2(I)*F2(J) &
     &                           +R3(I)*F3(J)+R4(I)*F4(J)
          ENDDO
        ENDDO
      ENDDO 
! PRINT*,'TORSION ',ETOT
! 
!     =========================================================================
!     ==  INVERSION FORCES                                                   ==
!     =========================================================================
!     PRINT*,'INVERSION',NINVERSION,ETOT
      DO IINV=1,NINVERSION
        ICOUNT=ICOUNT+1
        IF(MOD(ICOUNT-1,NTASKS).NE.THISTASK-1) CYCLE
        IAT1=INVERSION(IINV)%IAT1
        IAT2=INVERSION(IINV)%IAT2
        IAT3=INVERSION(IINV)%IAT3
        IAT4=INVERSION(IINV)%IAT4
        IPOT=INVERSION(IINV)%IPOT
        R1(:)=R(:,IAT1)
        R2(:)=R(:,IAT2)
        R3(:)=R(:,IAT3)
        R4(:)=R(:,IAT4)
        IF(ASSOCIATED(INVERSION(IINV)%IT2)) THEN
          R2(:)=R2(:)+RBAS(:,1)*REAL(INVERSION(IINV)%IT2(1),KIND=8) &
    &                +RBAS(:,2)*REAL(INVERSION(IINV)%IT2(2),KIND=8) &
    &                +RBAS(:,3)*REAL(INVERSION(IINV)%IT2(3),KIND=8)
        END IF
        IF(ASSOCIATED(INVERSION(IINV)%IT3)) THEN
          R3(:)=R3(:)+RBAS(:,1)*REAL(INVERSION(IINV)%IT3(1),KIND=8) &
    &                +RBAS(:,2)*REAL(INVERSION(IINV)%IT3(2),KIND=8) &
    &                +RBAS(:,3)*REAL(INVERSION(IINV)%IT3(3),KIND=8)
        END IF
        IF(ASSOCIATED(INVERSION(IINV)%IT4)) THEN
          R4(:)=R4(:)+RBAS(:,1)*REAL(INVERSION(IINV)%IT4(1),KIND=8) &
    &                +RBAS(:,2)*REAL(INVERSION(IINV)%IT4(2),KIND=8) &
    &                +RBAS(:,3)*REAL(INVERSION(IINV)%IT4(3),KIND=8)
        END IF
!       PRINT*,'INVERSION DETAILS',IAT1,IAT2,IAT3,IAT4,IPOT
        CALL CLASSICAL_EINVERSION(R1,R2,R3,R4,DE,F1,F2,F3,F4,POT(IPOT))
        ETOT=ETOT+DE
        F(:,IAT1)=F(:,IAT1)+F1(:)
        F(:,IAT2)=F(:,IAT2)+F2(:)
        F(:,IAT3)=F(:,IAT3)+F3(:)
        F(:,IAT4)=F(:,IAT4)+F4(:)
        DO I=1,3
          DO J=1,3
            SIGMA(I,J)=SIGMA(I,J)+R1(I)*F1(J)+R2(I)*F2(J) &
     &                           +R3(I)*F3(J)+R4(I)*F4(J)
          ENDDO
        ENDDO
      ENDDO 
!  PRINT*,'INVERSION ',ETOT
!
!     ==========================================================================
!     ==  ADD RESULTS FROM ALL TASKS                                          ==
!     ==========================================================================
      CALL MPE$COMBINE('MONOMER','+',ETOT)
      CALL MPE$COMBINE('MONOMER','+',F)
      CALL MPE$COMBINE('MONOMER','+',V)
!
!     == THIS IS AN OPTIONAL CONSISTENCY CHECK
!     CALL CLASSICAL_TESTSUMRULES(NAT,R,F)

                            CALL TRACE$POP
      RETURN
      END SUBROUTINE CLASSICAL_ETOTAL
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_TESTSUMRULES(NAT,R,F)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      REAL(8)   ,INTENT(IN) :: F(3,NAT)
      REAL(8)               :: FTOT(3)    ! GESAMTKRAFT
      REAL(8)               :: D(3)       ! DREHMOMENT
      INTEGER(4)            :: IAT
!     ******************************************************************
      FTOT(:)=0.D0
      D(:)=0.D0
      DO IAT=1,NAT
        FTOT(:)=FTOT(:)+F(:,IAT)
        D(1)=D(1)+F(2,IAT)*R(3,IAT)-F(3,IAT)*R(2,IAT)
        D(2)=D(2)+F(3,IAT)*R(1,IAT)-F(1,IAT)*R(3,IAT)
        D(3)=D(3)+F(1,IAT)*R(2,IAT)-F(2,IAT)*R(1,IAT)
      ENDDO
      WRITE(*,FMT='("SUMRULES ",6F15.7)')FTOT,D
      RETURN
      END SUBROUTINE CLASSICAL_TESTSUMRULES
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_ECOULOMB(NAT,R,Q,E,F,V,RBAS,SIGMA &
     &                  ,ITYPE,NTYPE,NONBOND &
     &                  ,NNB,NBLIST,NPOT,POT)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATE COULOMB AND VAN DER WAALS INTERACTION                      **
!     **                                                                      **
!     **  ATTENTION! THE RESULT IS DISTRIBUTED OVER ALL TASKS                 **
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE, ONLY : POT_TYPE,NONBOND_TYPE,MD
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)   :: NAT
      REAL(8)           ,INTENT(IN)   :: R(3,NAT)    ! ATOMIC POSITIONS
      REAL(8)           ,INTENT(IN)   :: Q(NAT)      ! CHARGES
      REAL(8)           ,INTENT(IN)   :: RBAS(3,3)   ! LATTICE TRANSLATIONS
      REAL(8)           ,INTENT(INOUT):: E           ! ENERGY
      REAL(8)           ,INTENT(INOUT):: F(3,NAT)    ! FORCES
      REAL(8)           ,INTENT(INOUT):: V(NAT)      ! POTENTIALS
      REAL(8)           ,INTENT(INOUT):: SIGMA(3,3)  ! R(I)*F(J)
      INTEGER(4)        ,INTENT(IN)   :: ITYPE(NAT)
      INTEGER(4)        ,INTENT(IN)   :: NTYPE
      INTEGER(4)        ,INTENT(IN)   :: NONBOND(NTYPE,NTYPE)  
      INTEGER(4)        ,INTENT(IN)   :: NNB
      TYPE(NONBOND_TYPE),INTENT(IN)   :: NBLIST(NNB)
      INTEGER(4)        ,INTENT(IN)   :: NPOT
      TYPE(POT_TYPE)    ,INTENT(IN)   :: POT(NPOT)
      REAL(8)           ,PARAMETER    :: PI=4.D0*ATAN(1.D0)
      LOGICAL(4)                      :: TEXCL
      INTEGER(4)                      :: IAT1,IAT2,IPOT,IIB
      INTEGER(4)                      :: ITYPE1,ITYPE2
      REAL(8)                         :: Q1,Q2
      REAL(8)                         :: D1,D2,D3
      REAL(8)                         :: X,FAC,SVAR
      REAL(8)                         :: RC,G,DGDX
      REAL(8)                         :: DE2,DEDX
      INTEGER(4)                      :: THISTASK,NTASKS,ICOUNT
      REAL(8)                         :: T1,T2,T3
LOGICAL(4) :: TONEFOUR   
!     **************************************************************************
                               CALL TRACE$PUSH('CLASSICAL_ECOULOMB')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      ICOUNT=0   ! USED TO DISTRIBUTE TASKS
      DO IIB=1,NNB
!       == DISTRIBUTE OVER TASKS
        ICOUNT=ICOUNT+1
        IF(MOD(ICOUNT-1,NTASKS).NE.THISTASK-1) CYCLE
!       == NOW CONTINUE ON THISTASK ...
        IAT1=NBLIST(IIB)%IAT1
        IAT2=NBLIST(IIB)%IAT2
        TEXCL=NBLIST(IIB)%EXCLUDE
        TONEFOUR=NBLIST(IIB)%TONEFOUR
        Q1=Q(IAT1)
        Q2=Q(IAT2)
        D1=R(1,IAT2)-R(1,IAT1)
        D2=R(2,IAT2)-R(2,IAT1)
        D3=R(3,IAT2)-R(3,IAT1)
        IF(ASSOCIATED(NBLIST(IIB)%IT)) THEN
          T1=REAL(NBLIST(IIB)%IT(1),KIND=8)
          T2=REAL(NBLIST(IIB)%IT(2),KIND=8)
          T3=REAL(NBLIST(IIB)%IT(3),KIND=8)
          D1=D1+RBAS(1,1)*T1+RBAS(1,2)*T2+RBAS(1,3)*T3
          D2=D2+RBAS(2,1)*T1+RBAS(2,2)*T2+RBAS(2,3)*T3
          D3=D3+RBAS(3,1)*T1+RBAS(3,2)*T2+RBAS(3,3)*T3
        END IF
        X=D1*D1+D2*D2+D3*D3
        IF(X.LT.1.D-20) X=1.D-20
        X=1.D0/SQRT(X)
        FAC=0.D0
!         
!       ========================================================================
!       ==  COULOMB INTERACTION INTERACTION                                   ==
!       ==  ASSUMES A CHARGE DISTRIBUTION WITH ATOM-CENTERED                  ==
!       ==  GAUSSIANS  RHO(R)=SUM Q_I*G(R-R_I)                                ==
!       ========================================================================
        IF(.NOT.TEXCL) THEN
          RC=1.D0       ! RC IS MEANT FOR USE WITH THE POTENTIAL OF GAUSSIANS 
                        ! INSTEAD OF POINT CHARGES TO SCALE THE WIDTH OF 
                        ! THE GAUSSIANS
          CALL VALUE(POT(1),RC*X,G,DGDX) ! POT(1) CONTAINS THE COULOMB 
                        ! POTENTIAL FOR ATOM PAIRS WITHOUT BOND EXCLUSIONS
          G=G/RC
        ELSE    ! SPECIAL TREATMENT OF THE ELECTROSTATIC OF EXCLUDED PAIRS
          IF(TONEFOUR) THEN   ! COULOMB INTERACTION FOR 1-4 EXCLUSIONS
            RC=1.D0     ! RC IS MEANT FOR USE WITH THE POTENTIAL OF GAUSSIANS 
                        ! INSTEAD OF POINT CHARGES TO SCALE THE WIDTH OF 
                        ! THE GAUSSIANS
            CALL VALUE(POT(2),RC*X,G,DGDX)!POT(1) CONTAINS THE COULOMB 
                        ! POTENTIAL FOR ATOM PAIRS WITHOUT BOND EXCLUSIONS
            G=G/RC
          ELSE ! COULOMB INTERACTION FOR 1-2 AND 1-3 EXCLUSIONS
            RC=4.D0     ! RC IS MEANT FOR USE WITH THE POTENTIAL OF GAUSSIANS 
                        ! INSTEAD OF POINT CHARGES TO SCALE THE WIDTH OF 
                        ! THE GAUSSIANS
            CALL VALUE(POT(3),RC*X,G,DGDX) ! POT(1) CONTAINS THE COULOMB 
                        ! POTENTIALFOR ATOM PAIRS WITHOUT BOND EXCLUSIONS
            G=G/RC
          END IF
        END IF
        E=E+Q1*Q2*G
        FAC=FAC+Q1*Q2*DGDX           
        V(IAT1)=V(IAT1)+Q2*G
        V(IAT2)=V(IAT2)+Q1*G
!         
!       ========================================================================
!       ==  LENNARD JONES INTERACTION                                         ==
!       ========================================================================
        IF((.NOT.TEXCL).OR.TONEFOUR) THEN
          ITYPE1=ITYPE(IAT1)
          ITYPE2=ITYPE(IAT2)
          IPOT=NONBOND(ITYPE1,ITYPE2)
          IF(IPOT.GT.0) THEN
            CALL VALUE(POT(IPOT),X,DE2,DEDX)
            IF(TONEFOUR) THEN  ! AMBER SCALES THE 1-4 VDW INTERACTION 
              DE2=DE2*MD%SCALEONEFOUR
              DEDX=DEDX*MD%SCALEONEFOUR
            END IF
            E=E+DE2
            FAC=FAC+DEDX
          END IF
        END IF
!       
!       ========================================================================
!       ==  CALCULATE FORCES                                                  ==
!       ========================================================================
        FAC=-FAC*X*X*X
        SVAR=FAC*D1
        F(1,IAT2)=F(1,IAT2)-SVAR
        F(1,IAT1)=F(1,IAT1)+SVAR
        SVAR=FAC*D2
        F(2,IAT2)=F(2,IAT2)-SVAR
        F(2,IAT1)=F(2,IAT1)+SVAR
        SVAR=FAC*D3
        F(3,IAT2)=F(3,IAT2)-SVAR
        F(3,IAT1)=F(3,IAT1)+SVAR
        SIGMA(1,:)=SIGMA(1,:)+FAC*D1**2
        SIGMA(2,:)=SIGMA(2,:)+FAC*D2**2
        SIGMA(3,:)=SIGMA(3,:)+FAC*D3**2
      ENDDO
                               CALL TRACE$POP
      RETURN
      END
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_EBOND(R1,R2,E,F1,F2,POT)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATES THE CONTRIBUTION OF BOND-LENGTH DISTORTIONS               **
!     **  TO TOTAL ENERGY AND FORCES.                                         **
!     **  THE POTENTIAL MUST BE GIVEN AS FUNCTION OF 1/|R1-R2|                **
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE, ONLY : POT_TYPE
      IMPLICIT NONE
      REAL(8)       ,INTENT(IN)  :: R1(3)
      REAL(8)       ,INTENT(IN)  :: R2(3)
      TYPE(POT_TYPE),INTENT(IN)  :: POT
      REAL(8)       ,INTENT(OUT) :: E
      REAL(8)       ,INTENT(OUT) :: F1(3)
      REAL(8)       ,INTENT(OUT) :: F2(3)
      REAL(8)                    :: D1,D2,D3
      REAL(8)                    :: X,FAC,SVAR,DEDX
!     **************************************************************************
      D1=R2(1)-R1(1)
      D2=R2(2)-R1(2)
      D3=R2(3)-R1(3)
      X=1.D0/SQRT(D1*D1+D2*D2+D3*D3)
      CALL VALUE(POT,X,E,DEDX)
      FAC=-DEDX*X*X*X
      SVAR=FAC*D1
      F2(1)=-SVAR
      F1(1)=+SVAR
      SVAR=FAC*D2
      F2(2)=-SVAR
      F1(2)=+SVAR
      SVAR=FAC*D3
      F2(3)=-SVAR
      F1(3)=+SVAR
      RETURN
      END
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_EANGLE(R1,R2,R3,E,F1,F2,F3,POT)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATES THE CONTRIBUTION OF BOND-ANGLE DISTORTIONS                **
!     **  TO TOTAL ENERGY AND FORCES.                                         **
!     **  THE POTENTIAL MUST BE GIVEN AS FUNCTION OF COS(PHI)                 **
!     **  WHERE PHI IS THE ANGLE BETWEEN THE BONDS TO R1 AND R3               **
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE, ONLY : POT_TYPE
      IMPLICIT NONE
      REAL(8)    ,INTENT(IN)   :: R1(3)
      REAL(8)    ,INTENT(IN)   :: R2(3)
      REAL(8)    ,INTENT(IN)   :: R3(3)
      REAL(8)    ,INTENT(OUT)  :: E
      REAL(8)    ,INTENT(OUT)  :: F1(3)
      REAL(8)    ,INTENT(OUT)  :: F2(3)
      REAL(8)    ,INTENT(OUT)  :: F3(3)
      TYPE(POT_TYPE),INTENT(IN):: POT
      REAL(8)                  :: D1D1,D1D2,D2D2
      INTEGER(4)               :: I
      REAL(8)                  :: D1,D2
      REAL(8)                  :: FAC,X,DEDX,SVAR
      REAL(8)                  :: A,B,C,DEDB1,DEDB2
!     **************************************************************************
      D1D1=0.D0
      D1D2=0.D0
      D2D2=0.D0
      DO I=1,3
        D1=R1(I)-R2(I)
        D2=R3(I)-R2(I)
        D1D1=D1D1+D1*D1
        D1D2=D1D2+D1*D2
        D2D2=D2D2+D2*D2
      ENDDO
      FAC=1.D0/SQRT(D1D1*D2D2)       
      X=D1D2*FAC
      CALL VALUE(POT,X,E,DEDX)
      SVAR=-DEDX*X
      A=SVAR/D1D1
      C=SVAR/D2D2
      B=DEDX*FAC
      DO I=1,3
        D1=R1(I)-R2(I)
        D2=R3(I)-R2(I)
        DEDB1=A*D1+B*D2
        DEDB2=B*D1+C*D2
        F1(I)=-DEDB1
        F2(I)=+DEDB1+DEDB2
        F3(I)=-DEDB2
      ENDDO
      RETURN
      END       
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_ETORSION(R1,R2,R3,R4,E,F1,F2,F3,F4,POT)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATES THE CONTRIBUTION OF BOND-TORSION DISTORTIONS              **
!     **  TO TOTAL ENERGY AND FORCES.                                         **
!     **  THE POTENTIAL MUST BE GIVEN AS FUNCTION OF COS(PHI)                 **
!     **                                                                      **
!     **  CALCULATES THE TORSION POTENTIAL BETWEEN THE BONDS B1=R1-R2         **
!     **  AND R43=R3-R4 ABOUT THE BOND B2=R2-R3.                              **
!     **  VARRAY CONTAINS THE POTENTIAL FOR N1*N2                             **
!     **  WHERE N1 AND N2 ARE THE NORMALIZED NORMAL VECTORS                   **
!     **  ON B1-B2 AND B2-R43 RESPECTIVELY                                    **
!     **  N1=CROSSPRODUCT(B1,B2)/SQRT(B1**2*B2**2)                            **
!     **  N2=CROSSPRODUCT(B2,R43)/SQRT(B2**2*R43**2)                          **
!     **                                                                      **
!     **  X=CROSS(R12,R32)*CROSS(R32,R43)                                     **
!     **       /SQRT(CROSS(R12,R12)**2*CROSS(R32,R43)**2)                     **
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE,ONLY : POT_TYPE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: R1(3)
      REAL(8)   ,INTENT(IN)   :: R2(3)
      REAL(8)   ,INTENT(IN)   :: R3(3)
      REAL(8)   ,INTENT(IN)   :: R4(3)
      REAL(8)   ,INTENT(OUT):: E
      REAL(8)   ,INTENT(OUT):: F1(3)
      REAL(8)   ,INTENT(OUT):: F2(3)
      REAL(8)   ,INTENT(OUT):: F3(3)
      REAL(8)   ,INTENT(OUT):: F4(3)
      TYPE(POT_TYPE),INTENT(IN):: POT
      REAL(8)                 :: R21(3),R32(3),R43(3)
      REAL(8)                 :: FR21(3),FR32(3),FR43(3)
      REAL(8)                 :: UX,UY,UZ,VX,VY,VZ
      REAL(8)                 :: UU,UV,VV
      REAL(8)                 :: ROOTUUVVINV
      REAL(8)                 :: X,DEDX
      REAL(8)                 :: A11,A12,A21,A22
      REAL(8)                 :: OX,OY,OZ,PX,PY,PZ
      REAL(8)                 :: SVAR
      INTEGER(4)              :: I
!     **************************************************************************
      DO I=1,3
         R21(I)=R2(I)-R1(I)
         R32(I)=R3(I)-R2(I)
         R43(I)=R4(I)-R3(I)
      ENDDO
!     == U IS THE NORMAL ON R1-R2-R3 ===========================================
      UX=R21(2)*R32(3)-R21(3)*R32(2)
      UY=R21(3)*R32(1)-R21(1)*R32(3)
      UZ=R21(1)*R32(2)-R21(2)*R32(1)
!     == V IS THE NORMAL ON R2-R3-R4 ===========================================
      VX=R43(2)*R32(3)-R43(3)*R32(2)
      VY=R43(3)*R32(1)-R43(1)*R32(3)
      VZ=R43(1)*R32(2)-R43(2)*R32(1)
!     ==
      UU=UX*UX+UY*UY+UZ*UZ
      VV=VX*VX+VY*VY+VZ*VZ
      UV=UX*VX+UY*VY+UZ*VZ
      IF(ABS(UU*VV).LT.1.D-8) THEN
        E=0.D0
        F1(:)=0.D0
        F2(:)=0.D0
        F3(:)=0.D0
        F4(:)=0.D0
        RETURN
      END IF
      ROOTUUVVINV=1.D0/SQRT(UU*VV)
!     == X IS THE COS(PHI) WHERE PHI IS THE ANGLE BETWEEN U AND V ==============
      X=UV*ROOTUUVVINV
!
      X=-X
      CALL VALUE(POT,X,E,DEDX)
      DEDX=-DEDX
!
!     == DE/DU=A11*U+A12*V ; DE/DV=A21*U+A22*V
      SVAR=DEDX*ROOTUUVVINV
      A21=SVAR
      A12=SVAR
      SVAR=SVAR*UV
      A11=-SVAR/UU
      A22=-SVAR/VV
!     ==  DE=O*DV+P*DU
      OX=A11*UX+A21*VX
      OY=A11*UY+A21*VY
      OZ=A11*UZ+A21*VZ
      PX=A12*UX+A22*VX
      PY=A12*UY+A22*VY
      PZ=A12*UZ+A22*VZ
      FR21(1)=  R32(2)*OZ-R32(3)*OY 
      FR21(2)=  R32(3)*OX-R32(1)*OZ 
      FR21(3)=  R32(1)*OY-R32(2)*OX 
      FR32(1)=-(R21(2)*OZ-R21(3)*OY+R43(2)*PZ-R43(3)*PY)
      FR32(2)=-(R21(3)*OX-R21(1)*OZ+R43(3)*PX-R43(1)*PZ)
      FR32(3)=-(R21(1)*OY-R21(2)*OX+R43(1)*PY-R43(2)*PX)
      FR43(1)=  R32(2)*PZ-R32(3)*PY 
      FR43(2)=  R32(3)*PX-R32(1)*PZ 
      FR43(3)=  R32(1)*PY-R32(2)*PX 
      DO I=1,3
        F1(I)=        +FR21(I)
        F2(I)=-FR21(I)+FR32(I)
        F3(I)=-FR32(I)+FR43(I)
        F4(I)=-FR43(I)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_EINVERSION(R1,R2,R3,R4,E,F1,F2,F3,F4,POT)
!     **************************************************************************
!     **                                                                      **
!     **  R1 IS THE CENTRAL ATOM                                              **
!     **  THE POTENTIAL DEPENDS ON THE SCALAR PRODUCT OF THE                  **
!     **  OF THE NORMAL ON THE PLANE SPANNED BY R1-R2-R3                      **
!     **  WITH THE NORMALIZED DIRECTION  R1-R4                                **
!     **                                                                      **
!     **  X=CROSSPRODUCT(R21,R31)*R41                                         **
!     **       /|CROSSPRODUCT(R21*R31)|*|R41|                                 **
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE,ONLY : POT_TYPE
      IMPLICIT NONE
      REAL(8)       ,INTENT(IN)    :: R1(3),R2(3),R3(3),R4(3)
      TYPE(POT_TYPE),INTENT(IN)    :: POT
      REAL(8)       ,INTENT(OUT)   :: E
      REAL(8)       ,INTENT(OUT)   :: F1(3),F2(3),F3(3),F4(3)
      REAL(8)                      :: R21X,R21Y,R21Z
      REAL(8)                      :: R31X,R31Y,R31Z
      REAL(8)                      :: R41X,R41Y,R41Z,R41L
      REAL(8)                      :: UX,UY,UZ,UL
      REAL(8)                      :: FUX,FUY,FUZ
      REAL(8)                      :: F31X,F31Y,F31Z
      REAL(8)                      :: F21X,F21Y,F21Z
      REAL(8)                      :: F41X,F41Y,F41Z
      REAL(8)                      :: UR41,FAC,X,DEDX,SVAR,A,B
!     **************************************************************************
      R21X=R2(1)-R1(1)
      R21Y=R2(2)-R1(2)
      R21Z=R2(3)-R1(3)
      R31X=R3(1)-R1(1)
      R31Y=R3(2)-R1(2)
      R31Z=R3(3)-R1(3)
      R41X=R4(1)-R1(1)
      R41Y=R4(2)-R1(2)
      R41Z=R4(3)-R1(3)
      UX=R21Y*R31Z-R21Z*R31Y    ! U=(R2-E1)X(R3-R1)
      UY=R21Z*R31X-R21X*R31Z
      UZ=R21X*R31Y-R21Y*R31X
      UL=UX*UX+UY*UY+UZ*UZ      ! UL=U**2
      IF(UL.LT.1.D-8) THEN
!       == U=0 I.E. R21 AND R31 ARE COLLINEAR
        E=0.D0
        F1(:)=0.D0
        F2(:)=0.D0
        F3(:)=0.D0
        F4(:)=0.D0
        RETURN
      END IF
      R41L=R41X*R41X+R41Y*R41Y+R41Z*R41Z  !R41L=(R4-R1)**2
      UR41=UX*R41X+UY*R41Y+UZ*R41Z
      FAC=1.D0/SQRT(UL*R41L)
      X=UR41*FAC
!
      CALL VALUE(POT,X,E,DEDX)
!
      SVAR=DEDX*FAC
      A=SVAR
      B=-SVAR*UR41/R41L
      F41X=A*UX+B*R41X
      F41Y=A*UY+B*R41Y
      F41Z=A*UZ+B*R41Z
      A=-SVAR*UR41/UL
      B=SVAR
      FUX=A*UX+B*R41X
      FUY=A*UY+B*R41Y
      FUZ=A*UZ+B*R41Z
      F31X=FUY*R21Z-FUZ*R21Y
      F31Y=FUZ*R21X-FUX*R21Z
      F31Z=FUX*R21Y-FUY*R21X
      F21X=FUZ*R31Y-FUY*R31Z 
      F21Y=FUX*R31Z-FUZ*R31X
      F21Z=FUY*R31X-FUX*R31Y
      F1(1)=+F21X+F31X+F41X
      F1(2)=+F21Y+F31Y+F41Y
      F1(3)=+F21Z+F31Z+F41Z
      F2(1)=-F21X
      F2(2)=-F21Y
      F2(3)=-F21Z
      F3(1)=-F31X
      F3(2)=-F31Y
      F3(3)=-F31Z
      F4(1)=-F41X
      F4(2)=-F41Y
      F4(3)=-F41Z
      RETURN
    END SUBROUTINE CLASSICAL_EINVERSION
! 
!     ==========================================================================
!     ==========================================================================
!     ====                                                                  ====
!     ====  ROUTINES THAT DEFINE THE POTENTIALS                             ====
!     ====                                                                  ====
!     ==========================================================================
!     ==========================================================================
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_FORCEFIELDSETUP(NBOND,BOND,BO &
     &                 ,NAT,TYPE,ITYPE &
     &                 ,NANGLEX,NANGLE,ANGLE &
     &                 ,NTORSIONX,NTORSION,TORSION &
     &                 ,NINVERSIONX,NINVERSION,INVERSION &
     &                 ,NTYPEX,NTYPE,NONBOND &
     &                 ,NPOTX,NPOT,POT)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE, ONLY : POT_TYPE,MD &
     &                           ,BOND_TYPE,ANGLE_TYPE,TORSION_TYPE &
     &                           ,INVERSION_TYPE
      IMPLICIT NONE
      INTEGER(4)   ,PARAMETER    :: NNEIGHX=30
      INTEGER(4)   ,INTENT(IN)   :: NBOND
      TYPE(BOND_TYPE),INTENT(INOUT) :: BOND(NBOND)
      REAL(8)      ,INTENT(IN)   :: BO(NBOND)
      INTEGER(4)   ,INTENT(IN)   :: NAT
      CHARACTER(*) ,INTENT(IN)   :: TYPE(NAT)
      INTEGER(4)   ,INTENT(OUT)  :: ITYPE(NAT)
      INTEGER(4)   ,INTENT(IN)   :: NANGLEX
      INTEGER(4)   ,INTENT(OUT)  :: NANGLE
      TYPE(ANGLE_TYPE),INTENT(OUT) :: ANGLE(NANGLEX)
      INTEGER(4)   ,INTENT(IN)   :: NTORSIONX
      INTEGER(4)   ,INTENT(OUT)  :: NTORSION
      TYPE(TORSION_TYPE),INTENT(OUT) :: TORSION(NTORSIONX)
      INTEGER(4)   ,INTENT(IN)   :: NINVERSIONX
      INTEGER(4)   ,INTENT(OUT)  :: NINVERSION
      TYPE(INVERSION_TYPE),INTENT(OUT) :: INVERSION(NINVERSIONX)
      INTEGER(4)   ,INTENT(IN)   :: NTYPEX
      INTEGER(4)   ,INTENT(OUT)  :: NTYPE
      INTEGER(4)   ,INTENT(OUT)  :: NONBOND(NTYPEX,NTYPEX)
      INTEGER(4)   ,INTENT(IN)   :: NPOTX
      INTEGER(4)   ,INTENT(OUT)  :: NPOT
      TYPE(POT_TYPE),INTENT(INOUT) :: POT(NPOTX)
      LOGICAL(4)                 :: TCHK
      CHARACTER(64)              :: ID
      CHARACTER(8)               :: TYPE1,TYPE2,TYPE3,TYPE4
!     INTEGER(4)   ,ALLOCATABLE  :: IWORK(:)
      INTEGER(4)                 :: NNEIGH(NAT)
      INTEGER(4)                 :: INEIGH(NNEIGHX,NAT)
      INTEGER(4)                 :: ITNEIGH(3,NNEIGHX,NAT) !TRANSLATION
      REAL(8)                    :: BONEIGH(NNEIGHX,NAT)
      INTEGER(4)                 :: IB,IPOT,IAT,I
      INTEGER(4)                 :: IAT1,IAT2,IAT3,IAT4,NN1,NN2,NN3
      INTEGER(4)                 :: IN2,IN3
      INTEGER(4)                 :: IT1(3),IT2(3),IT3(3),IT4(3)
      REAL(8)                    :: BO1,BO2,BO3,BO12,BO23,BO34,BO13,BO14
      INTEGER(4)                 :: NN,II,II1,II3,I1,I2,NIJ
      REAL(8)                    :: X,K,D
      LOGICAL(4),PARAMETER       :: TFOURIER=.FALSE.
!     **************************************************************************
                                  CALL TRACE$PUSH('CLASSICAL_FORCEFIELDSETUP')
      NPOT=0
!     ==========================================================================
!     ==  SET UP ELECTROSTATIC POTENTIAL FOR LONGRANGED INTERACTIONS          ==
!     ==========================================================================
      NPOT=NPOT+1
      IF(MD%FF.EQ.'UFF') THEN
        CALL CLASSICAL_COULOMBPOTB(POT(NPOT))   ! INTERACTING GAUSSIANS
      ELSE IF(MD%FF.EQ.'AMBER') THEN
        CALL CLASSICAL_COULOMBPOTA(POT(NPOT))   ! PURE COULOMB
      END IF
!
!     ==========================================================================
!     ==  SET UP ELECTROSTATIC POTENTIAL FOR 1-4 INTERACTIONS                 ==
!     ==========================================================================
      NPOT=NPOT+1
      IF(MD%FF.EQ.'UFF') THEN
        CALL CLASSICAL_COULOMBPOTB(POT(NPOT))   ! INTERACTING GAUSSIANS
      ELSE IF(MD%FF.EQ.'AMBER') THEN
        CALL CLASSICAL_COULOMBPOTA(POT(NPOT))   ! PURE COULOMB
        POT(NPOT)%VAL(:)=POT(NPOT)%VAL(:)/1.2D0
        POT(NPOT)%DER(:)=POT(NPOT)%DER(:)/1.2D0
        POT(NPOT)%ID=TRIM(POT(NPOT)%ID)//'(SCALED BY 5/6)'
      END IF
!
!     ==========================================================================
!     ==  SET UP ELECTROSTATIC POTENTIAL FOR 1-2 AND 1-3 INTERACTIONS         ==
!     ==========================================================================
      NPOT=NPOT+1
      IF(MD%FF.EQ.'UFF') THEN
        CALL CLASSICAL_COULOMBPOTB(POT(NPOT))   ! INTERACTING GAUSSIANS
      ELSE IF(MD%FF.EQ.'AMBER') THEN
        CALL CLASSICAL_COULOMBPOTB(POT(NPOT))   ! INTERACTING GAUSSIANS
!       POT(NPOT)%ID='CONSTANT POTENTIAL WITH VALUE -1 HARTREE '
!       POT*NPOT)%VAL(:)=-1.D0
!       POT*NPOT)%DER(:)=0.D0
      END IF
!
!     ==========================================================================
!     ==  SET SCALING OF 1-4 VDW INTERACTIONS                                 ==
!     ==========================================================================
      IF(MD%FF.EQ.'UFF') THEN
        MD%SCALEONEFOUR=1.D0
      ELSE IF(MD%FF.EQ.'AMBER') THEN
        MD%SCALEONEFOUR=0.5D0
      END IF
!     
!     ==========================================================================
!     ==  PREPARE BOND INDEX ARRAYS AND SET BOND DISTANCE POTENTIALS          ==
!     ==========================================================================
      DO IB=1,NBOND
!
!       ========================================================================
!       == WRITE IDENTIFIER FOR THE POTENTIAL                                 ==
!       ========================================================================
        IAT1=BOND(IB)%IAT1
        IAT2=BOND(IB)%IAT2
!
        IF(IAT1.GT.NAT.OR.IAT2.GT.NAT) THEN
          CALL ERROR$MSG('ATOM INDEX IN BOND OUT OF RANGE')
          CALL ERROR$I4VAL('IAT1',IAT1)
          CALL ERROR$I4VAL('IAT2',IAT2)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('UFFINITIALIZE')
        END IF
        TYPE1=TYPE(IAT1)
        TYPE2=TYPE(IAT2)
        IF(TRIM(ADJUSTL(MD%FF)).EQ.'UFF') THEN
           CALL UFFTABLE$BONDPARMS(TYPE1,TYPE2,BO(IB),ID,X,K,D,TCHK)
        ELSE IF(TRIM(ADJUSTL(MD%FF)).EQ.'AMBER') THEN
           CALL FORCEFIELD$AMBER_BONDPARMS(TYPE1,TYPE2,ID,X,K,D,TCHK)
!    -------- UFF USES K/2 WHERE AMBER USES K
           K= 2.D0 * K
!    ----------------------------------------
        ELSE
           CALL ERROR$MSG('FORCEFIELD NOT RECOGNIZED')
           CALL ERROR$CHVAL('MD%FF',MD%FF)
           CALL ERROR$STOP('FORCEFIELDSETUP')
        END IF
!
!       =================================================================
!       == SEARCH WHETHER POTENTIAL ALREADY EXISTS                     ==
!       == AND CALCULATE A NEW ONE IF IT DOES NOT                      ==
!       =================================================================
        TCHK=.FALSE.
        DO IPOT=1,NPOT
          IF(ID.EQ.POT(IPOT)%ID) THEN
            BOND(IB)%IPOT=IPOT
            TCHK=.TRUE.
            EXIT
          END IF
        ENDDO
        IF(TCHK) CYCLE
!
!       == CREATE NEW POTENTIAL IF NECCESARY ===================================
        NPOT=NPOT+1
        IF(NPOT.GT.NPOTX) THEN
          CALL ERROR$MSG('NPOT TOO LARGE')
          CALL ERROR$STOP('CLASSICAL_UFFINITIALIZE')
        END IF
        CALL CLASSICAL_BONDPOTA(ID,X,K,D,POT(NPOT))
        BOND(IB)%IPOT=NPOT
      ENDDO
!
!     ==========================================================================
!     ==  DEFINE NEIGHBOR LIST                                                ==
!     ==========================================================================
      NNEIGH(:)=0
      INEIGH(:,:)=0
      DO IB=1,NBOND
        IAT1=BOND(IB)%IAT1
        IAT2=BOND(IB)%IAT2
        IF(ASSOCIATED(BOND(IB)%IT2)) THEN
          IT2=BOND(IB)%IT2
        ELSE
          IT2(:)=0
        END IF
        BO1=BO(IB)
        NN1=NNEIGH(IAT1)+1
        NN2=NNEIGH(IAT2)+1
        IF(NN1.GT.NNEIGHX.OR.NN2.GT.NNEIGHX) THEN
          CALL ERROR$STOP('CLASSICAL_UFFINITIALIZE')
        END IF
        NNEIGH(IAT1)=NN1
        INEIGH(NN1,IAT1)=IAT2
        ITNEIGH(:,NN1,IAT1)=IT2
        BONEIGH(NN1,IAT1)=BO1
        NNEIGH(IAT2)=NN2
        INEIGH(NN2,IAT2)=IAT1
        ITNEIGH(:,NN2,IAT2)=-IT2(:)
        BONEIGH(NN2,IAT2)=BO1
      ENDDO
      DO IAT=1,NAT
!       PRINT*,'TYPE   ',TYPE(IAT)
        PRINT*,'IAT',IAT,'INEIGH ',(INEIGH(I,IAT),I=1,NNEIGH(IAT))
!       PRINT*,'TYPE   ',(TYPE(INEIGH(I,IAT)),I=1,NNEIGH(IAT))
!       PRINT*,'BO     ',(BONEIGH(I,IAT),I=1,NNEIGH(IAT))
      ENDDO
!     STOP
!
!     ==========================================================================
!     ==  SET BOND ANGLE POTENTIALS                                           ==
!     ==========================================================================
!     PRINT*,'BEFORE ANGLE-LOOP',NBOND
      NANGLE=0
      DO IAT2=1,NAT
        NN2=NNEIGH(IAT2)
        DO II1=1,NN2
          DO II3=II1+1,NN2
            IAT1=INEIGH(II1,IAT2)
            IT1(:)=ITNEIGH(:,II1,IAT2)
            BO1=BONEIGH(II1,IAT2)
            IAT3=INEIGH(II3,IAT2)
            IT3(:)=ITNEIGH(:,II3,IAT2)
            BO3=BONEIGH(II3,IAT2)
            TYPE1=TYPE(IAT1)
            TYPE2=TYPE(IAT2)
            TYPE3=TYPE(IAT3)
            IF(TRIM(ADJUSTL(MD%FF)).EQ.'UFF') THEN
               CALL UFFTABLE$ANGLEPARMS(TYPE1,TYPE2,TYPE3,BO1,BO3,ID,X,K,TCHK)
            ELSE IF(TRIM(ADJUSTL(MD%FF)).EQ.'AMBER') THEN
               CALL FORCEFIELD$AMBER_ANGLEPARMS(TYPE1,TYPE2,TYPE3,ID,X,K,TCHK)
            ELSE
               CALL ERROR$MSG('FORCEFIELD NOT RECOGNIZED')
               CALL ERROR$CHVAL('MD%FF',MD%FF)
               CALL ERROR$STOP('FORCEFIELDSETUP')
            END IF
            IF(.NOT.TCHK) CYCLE
!
!           == ADD NEW ANGLE ===================================================
            NANGLE=NANGLE+1
            IF(NANGLE.GT.NANGLEX) THEN
              CALL ERROR$MSG('NUMBER OF ANGLES LARGER THAN MAXIMUM')
              CALL ERROR$OVERFLOW('NANGLE',NANGLE,NANGLEX)
              CALL ERROR$STOP('CLASSICAL_UFFINITIALIZE')
            END IF
            ANGLE(NANGLE)%IAT1=IAT1
            ANGLE(NANGLE)%IAT2=IAT2
            ANGLE(NANGLE)%IAT3=IAT3
            ANGLE(NANGLE)%IPOT=0
            IF(IT1(1).NE.0.OR.IT1(2).NE.0.OR.IT1(3).NE.0) THEN
              ALLOCATE(ANGLE(NANGLE)%IT1(3))
              ANGLE(NANGLE)%IT1(:)=IT1(:)
            END IF
            IF(IT3(1).NE.0.OR.IT3(2).NE.0.OR.IT3(3).NE.0) THEN
              ALLOCATE(ANGLE(NANGLE)%IT3(3))
              ANGLE(NANGLE)%IT3(:)=IT3(:)
            END IF
!
!           == SEARCH WHETHER POTENTIAL ALREADY EXISTS =========================
            TCHK=.FALSE.
            DO IPOT=1,NPOT
              IF(ID.EQ.POT(IPOT)%ID) THEN
                ANGLE(NANGLE)%IPOT=IPOT
                TCHK=.TRUE.
                EXIT
              END IF
            ENDDO
            IF(TCHK) CYCLE
!
!           == CREATE NEW POTENTIAL IF IT DOES NOTE YET EXIST  =================
            NPOT=NPOT+1
            IF(NPOT.GT.NPOTX) THEN
              CALL ERROR$MSG('NUMBER OF POTENTIALS TOO LARGE')
              CALL ERROR$STOP('CLASSICAL_UFFINITIALIZE')
            END IF
            IF(TRIM(ADJUSTL(MD%FF)).EQ.'UFF') THEN
               CALL CLASSICAL_ANGLEPOTA(ID,X,K,POT(NPOT))
            ELSE IF(TRIM(ADJUSTL(MD%FF)).EQ.'AMBER') THEN
               IF(TFOURIER) THEN
                  K=2.D0*K  ! HERE WE USE K^F = 2*K^H/N^2 (N=1)
                  CALL CLASSICAL_ANGLEPOTA(ID,X,K,POT(NPOT))
               ELSE
                  CALL FORCEFIELD$AMBER_ANGLEPOTA(ID,X,K,POT(NPOT))
               END IF
            END IF
            ANGLE(NANGLE)%IPOT=NPOT
          ENDDO
        ENDDO
      ENDDO
!
!     == TORSIONS ==============================================================
!     PRINT*,'BEFORE TORSION-LOOP: NANGLE=',NANGLE
      NTORSION=0
      DO IB=1,NBOND
        IAT2=BOND(IB)%IAT1
        IAT3=BOND(IB)%IAT2
        IT2(:)=0
        IF(ASSOCIATED(BOND(IB)%IT2)) THEN
          IT3(:)=BOND(IB)%IT2(:)
        ELSE
          IT3(:)=0
        END IF
        BO23=BO(IB)
        NN2=NNEIGH(IAT2)
        NN3=NNEIGH(IAT3)
        DO IN2=1,NN2
          IAT1=INEIGH(IN2,IAT2)
          IT1(:)=ITNEIGH(:,IN2,IAT2)
!         ____AVOID USING CENTRAL BOND AS TERMINAL BOND_________________________
          IF(IAT1.EQ.IAT3) THEN
            IF(IT1(1)-IT3(1).EQ.0) THEN
              IF(IT1(2)-IT3(2).EQ.0) THEN
                IF(IT1(3)-IT3(3).EQ.0) THEN
                  CYCLE
                END IF
              END IF
            END IF
          END IF
          DO IN3=1,NN3
            IAT4=INEIGH(IN3,IAT3)
            IT4(:)=ITNEIGH(:,IN3,IAT3)+IT3(:)
!           ==  AVOID USING CENTRAL BOND AS TERMINAL BOND ======================
            IF(IAT4.EQ.IAT2) THEN
              IF(IT4(1)-IT2(1).EQ.0) THEN
                IF(IT4(2)-IT2(2).EQ.0) THEN
                  IF(IT4(3)-IT2(3).EQ.0) THEN
                    CYCLE
                  END IF
                END IF
              END IF
            END IF
!
!           == WRITE IDENTIFIER FOR THE POTENTIAL ==============================
            TYPE1=TYPE(IAT1)
            TYPE2=TYPE(IAT2)
            TYPE3=TYPE(IAT3)
            TYPE4=TYPE(IAT4)
            BO12=0.D0
            BO34=0.D0
            IF(TRIM(ADJUSTL(MD%FF)).EQ.'UFF') THEN
               CALL UFFTABLE$TORSIONPARMS(TYPE1,TYPE2,TYPE3,TYPE4 &
                    &                       ,BO12,BO23,BO34,ID,X,NIJ,K,TCHK)
            ELSE IF(TRIM(ADJUSTL(MD%FF)).EQ.'AMBER') THEN
               CALL FORCEFIELD$AMBER_TORSIONPARMS(TYPE1,TYPE2,TYPE3,TYPE4 &
         &                                       ,ID,X, NIJ, K, TCHK)
            ELSE
               CALL ERROR$MSG('FORCEFIELD NOT RECOGNIZED')
               CALL ERROR$CHVAL('MD%FF',MD%FF)
               CALL ERROR$STOP('FORCEFIELDSETUP')
            END IF
               
            IF(.NOT.TCHK) CYCLE
!
!           == ADD NEW TORSION =================================================
            NTORSION=NTORSION+1
            IF(NTORSION.GT.NTORSIONX) THEN
              CALL ERROR$MSG('NUMBER OF TORSIONS LARGER THAN MAXIMUM')
              CALL ERROR$I4VAL('NTORSION',NTORSION)
              CALL ERROR$I4VAL('NTORSIONX',NTORSIONX)
              CALL ERROR$STOP('CLASSICAL_UFFINITIALIZE')
            END IF
            TORSION(NTORSION)%IAT1=IAT1
            TORSION(NTORSION)%IAT2=IAT2
            TORSION(NTORSION)%IAT3=IAT3
            TORSION(NTORSION)%IAT4=IAT4
            IF(IT1(1).NE.0.OR.IT1(2).NE.0.OR.IT1(3).NE.0) THEN
              ALLOCATE(TORSION(NTORSION)%IT1(3))
              TORSION(NTORSION)%IT1(:)=IT1(:)
            END IF
            IF(IT3(1).NE.0.OR.IT3(2).NE.0.OR.IT3(3).NE.0) THEN
              ALLOCATE(TORSION(NTORSION)%IT3(3))
              TORSION(NTORSION)%IT3(:)=IT3(:)
            END IF
            IF(IT4(1).NE.0.OR.IT4(2).NE.0.OR.IT4(3).NE.0) THEN
              ALLOCATE(TORSION(NTORSION)%IT4(3))
              TORSION(NTORSION)%IT4(:)=IT4(:)
            END IF
              
!
!           == SEARCH WHETHER POTENTIAL ALREADY EXISTS =========================
            TCHK=.FALSE.
            DO IPOT=1,NPOT
              IF(ID.EQ.POT(IPOT)%ID) THEN
                TORSION(NTORSION)%IPOT=IPOT
                TCHK=.TRUE.
                EXIT
              END IF
            ENDDO
            IF(TCHK) CYCLE
!
!           == CREATE NEW POTENTIAL IF IT DOES NOTE YET EXIST  =================
            NPOT=NPOT+1
            IF(NPOT.GT.NPOTX) THEN
              CALL ERROR$MSG('NUMBER OF POTENTIALS TOO LARGE')
              CALL ERROR$STOP('CLASSICAL_UFFINITIALIZE')
            END IF
            IF(TRIM(ADJUSTL(MD%FF)).EQ.'UFF') THEN
               CALL CLASSICAL_TORSIONPOTA(ID,X,NIJ,K,POT(NPOT))
            ELSE IF(TRIM(ADJUSTL(MD%FF)).EQ.'AMBER') THEN
               CALL FORCEFIELD$AMBER_TORSIONPOTA(ID,X,NIJ,K,POT(NPOT))
            END IF
            TORSION(NTORSION)%IPOT=NPOT
          ENDDO
        ENDDO
      ENDDO
!
!     == INVERSION =============================================================
!     == AMBER APPARENTLY DOES NOT USE INVERSIONS....
      IF(TRIM(ADJUSTL(MD%FF)).NE.'AMBER') THEN
        NINVERSION=0
        DO IAT=1,NAT
          NN=NNEIGH(IAT)
          IT1(:)=0
          IF(NN.EQ.3) THEN
            DO II=1,3
              IAT2=INEIGH(1+MOD(II-1,3),IAT)
              IAT3=INEIGH(1+MOD(II  ,3),IAT)
              IAT4=INEIGH(1+MOD(II+1,3),IAT)
              IT2(:)=ITNEIGH(:,1+MOD(II-1,3),IAT)
              IT3(:)=ITNEIGH(:,1+MOD(II  ,3),IAT)
              IT4(:)=ITNEIGH(:,1+MOD(II+1,3),IAT)
              TYPE1=TYPE(IAT)
              TYPE2=TYPE(IAT2)
              TYPE3=TYPE(IAT3)
              TYPE4=TYPE(IAT4)
!             NO BO NEEDED FOR INVERSION                               
!             BO12=BONEIGH(1+MOD(II-1,3),IAT)
!             BO13=BONEIGH(1+MOD(II  ,3),IAT)
!             BO14=BONEIGH(1+MOD(II+1,3),IAT)
              CALL UFFTABLE$INVERSIONPARMS(TYPE1,TYPE2,TYPE3,TYPE4 &
     &              ,ID,X,K,TCHK)
              IF(.NOT.TCHK) CYCLE
!
              NINVERSION=NINVERSION+1
!     PRINT*,'NINVERSION',NINVERSION
 
              INVERSION%IAT1=IAT
              INVERSION%IAT2=IAT2
              INVERSION%IAT3=IAT3
              INVERSION%IAT4=IAT4
              IF(IT2(1).NE.0.OR.IT2(2).NE.0.OR.IT2(3).NE.0) THEN
                ALLOCATE(INVERSION(NINVERSION)%IT2(3))
                INVERSION(NINVERSION)%IT2(:)=IT2(:)
              END IF
              IF(IT3(1).NE.0.OR.IT3(2).NE.0.OR.IT3(3).NE.0) THEN
                ALLOCATE(INVERSION(NINVERSION)%IT3(3))
                INVERSION(NINVERSION)%IT3(:)=IT3(:)
              END IF
              IF(IT4(1).NE.0.OR.IT4(2).NE.0.OR.IT4(3).NE.0) THEN
                ALLOCATE(INVERSION(NINVERSION)%IT4(3))
                INVERSION(NINVERSION)%IT4(:)=IT4(:)
              END IF
!
!             == SEARCH WHETHER POTENTIAL ALREADY EXISTS =======================
              TCHK=.FALSE.
              DO IPOT=1,NPOT
                IF(ID.EQ.POT(IPOT)%ID) THEN
                  INVERSION%IPOT=IPOT
                  TCHK=.TRUE.
                  EXIT
                END IF
              ENDDO
              IF(TCHK) CYCLE
!        
!             == CREATE NEW POTENTIAL IF IT DOES NOTE YET EXIST  ===============
              NPOT=NPOT+1
              IF(NPOT.GT.NPOTX) THEN
                CALL ERROR$MSG('NUMBER OF POTENTIALS TOO LARGE')
                CALL ERROR$STOP('CLASSICAL_UFFINITIALIZE')
              END IF
              INVERSION%IPOT=NPOT
              CALL CLASSICAL_INVERSIONPOTA(ID,X,K,POT(NPOT))
            ENDDO
          END IF
        ENDDO
      END IF
!
!     ==========================================================================
!     ==   DEFINE POINTER TO ATOM TYPES USED                                  ==
!     ==========================================================================
      NTYPE=0
      DO IAT1=1,NAT
        TYPE1=TYPE(IAT1)
        TCHK=.FALSE.
        DO IAT2=1,IAT1-1
          TYPE2=TYPE(IAT2)
          IF(TYPE1.EQ.TYPE2) THEN
            ITYPE(IAT1)=ITYPE(IAT2)
            TCHK=.TRUE.
          END IF
        ENDDO
        IF(.NOT.TCHK) THEN
          NTYPE=NTYPE+1
          ITYPE(IAT1)=NTYPE
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  NONBOND                                                             ==
!     ==========================================================================
!     PRINT*,'BEFORE NONBOND LOOP '
      NONBOND(:,:)=0
      DO IAT1=1,NAT
        DO IAT2=IAT1,NAT
          I1=ITYPE(IAT1)
          I2=ITYPE(IAT2)
          IF(NONBOND(I1,I2).EQ.0) THEN
            TYPE1=TYPE(IAT1)
            TYPE2=TYPE(IAT2)
            IF(TRIM(ADJUSTL(MD%FF)).EQ.'UFF') THEN
              CALL UFFTABLE$NONBONDPARMS(TYPE1,TYPE2,ID,X,D,TCHK)
            ELSE IF(TRIM(ADJUSTL(MD%FF)).EQ.'AMBER') THEN
              CALL FORCEFIELD$AMBER_NONBONDPARMS(TYPE1,TYPE2,ID,X,D,TCHK)
            END IF
            IF(TCHK) THEN
              NPOT=NPOT+1
              IF(NPOT.GT.NPOTX) THEN
                 CALL ERROR$MSG('NUMBER OF POTENTIALS TOO LARGE AT NONBOND')
                 CALL ERROR$STOP('CLASSICAL_UFFINITIALIZE')
              END IF
              IF(TRIM(ADJUSTL(MD%FF)).EQ.'UFF') THEN
                 CALL CLASSICAL_NONBONDPOTA(ID,X,D,POT(NPOT))
              ELSE IF(TRIM(ADJUSTL(MD%FF)).EQ.'AMBER') THEN
                 CALL FORCEFIELD$AMBER_NONBONDPOTA(ID,X,D,POT(NPOT))
              END IF
              NONBOND(I1,I2)=NPOT
              NONBOND(I2,I1)=NPOT
            ELSE
              NONBOND(I1,I2)=-1
              NONBOND(I2,I1)=-1
            END IF
          END IF
        ENDDO
      ENDDO
                                  CALL TRACE$POP
      RETURN
    END SUBROUTINE CLASSICAL_FORCEFIELDSETUP
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_COULOMBPOTA(POT)
!     ******************************************************************
!     **  RETURNS RC*F(RC/R) AS FUNCTION OF X=RC/R                    **
!     **  WHERE RC*F(RC/R) APPROACHES 1/R FOR R->INFTY                **
!     **                                                              **
!     ******************************************************************
      USE CLASSICAL_MODULE, ONLY: POT_TYPE
      IMPLICIT NONE
      TYPE(POT_TYPE),INTENT(OUT) :: POT
      INTEGER(4)    ,PARAMETER   :: NX=1000
      REAL(8)       ,PARAMETER   :: X1=0.D0
      REAL(8)       ,PARAMETER   :: X2=5.D0
      REAL(8)       ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)                 :: I
      REAL(8)                    :: X
      REAL(8)                    :: ERFX
!     ******************************************************************
      POT%ID='COULOMB WITH 1/R'
      POT%NX=NX
      POT%X1=X1
      POT%DX=(X2-X1)/DBLE(POT%NX-1)
      ALLOCATE(POT%VAL(POT%NX))
      ALLOCATE(POT%DER(POT%NX))
      DO I=1,POT%NX
        X=POT%X1+POT%DX*DBLE(I-1)
        IF(X.NE.0.D0) THEN
          CALL LIB$ERFR8(1.D0/X,ERFX)
!           POT%VAL(I)=X*ERFX    
!           POT%DER(I)=ERFX-2.D0/(X*SQRT(PI))*EXP(-1.D0/X**2) 
          POT%VAL(I)=X
          POT%DER(I)=1.D0
       ELSE
          POT%VAL(I)=0.D0
          POT%DER(I)=1.D0
        END IF
      ENDDO
!      
!       CALL FILEHANDLER$UNIT('INFO',NFILINFO)
!       DO I=1,POT%NX
!          X=POT%X1+POT%DX*DBLE(I-1)
!          IF(X.NE.0.D0) THEN
!             WRITE(NFILINFO,FMT='(3F15.6)') 1/X,POT%VAL(I),POT%DER(I)
!          END IF
!       END DO
!       STOP 'FORCED STOP IN CLASSICAL_COULOMBPOTA'
!         
      RETURN
    END SUBROUTINE CLASSICAL_COULOMBPOTA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_COULOMBPOTB(POT)
!     ******************************************************************
!     **                                                              **
!     **  RETURNS RC*F(RC/R) AS FUNCTION OF X=RC/R                    **
!     **  WHERE RC*F(RC/R) APPROACHES 1/R FOR R->INFTY                **
!     **                                                              **
!     ******************************************************************
      USE CLASSICAL_MODULE, ONLY: POT_TYPE
      IMPLICIT NONE
      TYPE(POT_TYPE),INTENT(OUT) :: POT
      INTEGER(4)    ,PARAMETER   :: NX=1000
      REAL(8)       ,PARAMETER   :: X1=0.D0
      REAL(8)       ,PARAMETER   :: X2=5.D0
      REAL(8)       ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)                 :: I
      REAL(8)                    :: X
      REAL(8)                    :: ERFX
!     ******************************************************************
      POT%ID='COULOMB WITH GAUSS CUTOFF'
      POT%NX=NX
      POT%X1=X1
      POT%DX=(X2-X1)/DBLE(POT%NX-1)
      ALLOCATE(POT%VAL(POT%NX))
      ALLOCATE(POT%DER(POT%NX))
      DO I=1,POT%NX
        X=POT%X1+POT%DX*DBLE(I-1)
        IF(X.NE.0.D0) THEN
          CALL LIB$ERFR8(1.D0/X,ERFX)
          POT%VAL(I)=X*ERFX    
          POT%DER(I)=ERFX-2.D0/(X*SQRT(PI))*EXP(-1.D0/X**2) 
        ELSE
          POT%VAL(I)=0.D0
          POT%DER(I)=1.D0
        END IF
      ENDDO
      RETURN
    END SUBROUTINE CLASSICAL_COULOMBPOTB
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_BONDPOTA(ID,RIJ,KIJ,DIJ,POT)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE CLASSICAL_MODULE, ONLY: POT_TYPE
      IMPLICIT NONE
      CHARACTER(*)  ,INTENT(IN)  :: ID
      REAL(8)       ,INTENT(IN)  :: RIJ     ! OPTIMUM BONDLENGTH
      REAL(8)       ,INTENT(IN)  :: KIJ     !FORCE CONSTANT
      REAL(8)       ,INTENT(IN)  :: DIJ     !FORCE CONSTANT
      TYPE(POT_TYPE),INTENT(OUT) :: POT
      INTEGER(4)    ,PARAMETER   :: NX=1000
      REAL(8)       ,PARAMETER   :: X1=1.D-4 
      REAL(8)       ,PARAMETER   :: X2=5.D0
      LOGICAL(4)    ,PARAMETER   :: THARMONIC=.TRUE.
      INTEGER(4)                 :: I
      REAL(8)                    :: X,R,DRDX,SVAR,ALPHA
!     ******************************************************************
      POT%ID=ID
      POT%NX=NX
      POT%X1=X1
      POT%DX=(X2-X1)/DBLE(POT%NX-1)
      ALLOCATE(POT%VAL(POT%NX))
      ALLOCATE(POT%DER(POT%NX))
      DO I=1,POT%NX
        X=POT%X1+POT%DX*DBLE(I-1)
!       == POTENTIAL DIVIDE-BY-ZERO: SEE X1
        R=1/X
        DRDX=-R**2
        IF(THARMONIC) THEN
          POT%VAL(I)=0.5D0*KIJ*(R-RIJ)**2
          POT%DER(I)=KIJ*(R-RIJ)*DRDX
        ELSE   !MORSE TYPE POTENTIAL
          ALPHA=SQRT(KIJ/(2.D0*DIJ))
          SVAR=EXP(-ALPHA*(R-RIJ))
          POT%VAL(I)=DIJ*(SVAR-1.D0)**2
          POT%DER(I)=DIJ*2.D0*(SVAR-1.D0)*(-ALPHA*SVAR)*DRDX
        END IF
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_ANGLEPOTA(ID,THETA0,KIJK,POT)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE CLASSICAL_MODULE, ONLY: POT_TYPE
      IMPLICIT NONE
      CHARACTER(*)  ,INTENT(IN)  :: ID
      REAL(8)       ,INTENT(IN)  :: THETA0   ! OPTIMUM ANGLE
      REAL(8)       ,INTENT(IN)  :: KIJK     !FORCE CONSTANT
      TYPE(POT_TYPE),INTENT(OUT) :: POT
      INTEGER(4)    ,PARAMETER   :: NX=101
      REAL(8)       ,PARAMETER   :: X1=-1.D0
      REAL(8)       ,PARAMETER   :: X2=+1.D0
      REAL(8)       ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)                 :: I,N
      REAL(8)                    :: C(0:4)
      REAL(8)                    :: SVAR1,SVAR2,THETA,THETAN,A,B
      REAL(8)                    :: X
!     ******************************************************************
!
!     ==========================================================================
!     ==  DETERMINE FOURIER COMPONENTS                                        ==
!     ==========================================================================
      C(:)=0.D0
      IF(ABS(THETA0-0.5D0*PI).LT.1.D-5) THEN
!     ==  SPECIAL CASE 90 DEGREE: OCTAHEDRAL OR SQUARE-PLANAR ==========
        C(0)=KIJK/4.D0**2
        C(4)=-C(0)
!     == SPECIAL CASE 180 DEGREE (LINEAR) TO AVOID DIVIDE-BY-ZERO ======
      ELSE IF(ABS(THETA0-PI).LT.1.D-5) THEN
        C(0)=-KIJK   ! EQ. UFF10 HAS A SIGN ERROR FOR THIS SPECIAL CASE
        C(1)=-C(0)
!     == SPECIAL CASE 120 DEGREE (TRIGONAL-PLANAR) =====================
      ELSE IF(ABS(THETA0-2.D0*PI/3.D0).LT.1.D-5) THEN
        C(0)=KIJK/3.D0**2
        C(3)=-C(0)
      ELSE
!     == GENERAL NONLINEAR CASE  =======================================
!     == USE COS(2X)=2*COS(X)**2-1 TO GET A POLYNOMIAL IN COS(X) =======
!     == V=(C0-C2) + C1*COS(X) + 2*C2*COS(X)**2  =======================
!     == THE MINIMUM IS AT COS(THETA) AND THE VALUE AT THE MINIMUM IS ZERO
        C(2)=KIJK/(4.D0*SIN(THETA0)**2)
        C(1)=-4.D0*C(2)*COS(THETA0)
        C(0)=C(2)*(2.D0*COS(THETA0)**2+1.D0)
      END IF
!
!     ==========================================================================
!     ==  DETERMINE ARRAY SIZE ETC                                            ==
!     ==========================================================================
      POT%ID=ID
      POT%NX=NX
      POT%X1=X1
      POT%DX=(X2-X1)/DBLE(POT%NX-1)
      ALLOCATE(POT%VAL(POT%NX))
      ALLOCATE(POT%DER(POT%NX))
      DO I=1,POT%NX
        X=POT%X1+POT%DX*DBLE(I-1) 
        THETA=DACOS(X)
        SVAR1=0.D0
        SVAR2=0.D0
        IF(I.NE.1.AND.I.NE.POT%NX) THEN
          DO N=0,4
            THETAN=THETA*DBLE(N)
            SVAR1=SVAR1+C(N)*COS(THETAN)
            SVAR2=SVAR2+DBLE(N)*C(N)*SIN(THETAN)/SIN(THETA)
          ENDDO
        ELSE
!         == AVOID A DIVIDE BY ZERO AT THETA=0 AND THETA=PI
          DO N=0,4
            THETAN=THETA*DBLE(N)
            SVAR1=SVAR1+C(N)*COS(THETAN)
            SVAR2=SVAR2+DBLE(N)**2*C(N)*COS(THETAN)/COS(THETA)
          ENDDO
        END IF
        POT%VAL(I)=SVAR1
        POT%DER(I)=SVAR2
      ENDDO
!
!     ==================================================================
!     ==  REMOVE THE POTENTIAL-WELL AT THETA=0                        ==
!     ==  (NECCESARY TO AVOID MINIMA FOR BOND ANGLES AT THETA=0)      ==
!     ==================================================================
      IF(POT%DER(POT%NX).GT.0) RETURN
      DO I=POT%NX-1,1,-1
        N=I
        IF(POT%VAL(I).LT.POT%VAL(I+1))EXIT
      ENDDO
      X=POT%X1+POT%DX*DBLE(N-1)
      B=POT%DER(N)/(2.D0*(X-1.D0))
      A=POT%VAL(N)-B*(X-1.D0)**2
      IF(B.GT.0.D0) THEN
        CALL ERROR$MSG('ATTEMPT TO REMOVE MINIMUM AT ZERO BOND-ANGLE FAILED')
        CALL ERROR$I4VAL('N',N)
        CALL ERROR$STOP('CLASSICAL%ANGLEPOTA')
      END IF
      DO I=N,POT%NX
        X=POT%X1+POT%DX*DBLE(I-1)
        POT%VAL(I)=A+B*(X-1.D0)**2
        POT%DER(I)=2.D0*B*(X-1.D0)
      ENDDO
      RETURN 
    END SUBROUTINE CLASSICAL_ANGLEPOTA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_TORSIONPOTA(ID,PHI0,NJK,V,POT)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES A TORSION POTENTIAL FROM SPECIFICATIONS          **
!     **                                                              **
!     **  THE POTENTIAL IS A FUNCTION OF COS(PHI)                     **
!     **                                                              **
!     ******************************************************************
      USE CLASSICAL_MODULE,ONLY: POT_TYPE
      IMPLICIT NONE
      CHARACTER(*)  ,INTENT(IN) :: ID   !POTENTIAL IDENTIFIER
      REAL(8)       ,INTENT(IN) :: PHI0 !IDEAL TORSION ANGLE
      INTEGER(4)    ,INTENT(IN) :: NJK  !MULTIPLICITY OF MINIMA
      REAL(8)       ,INTENT(IN) :: V    !TORSIONAL BARRIER HEIGHT
      TYPE(POT_TYPE),INTENT(OUT):: POT
      INTEGER(4)    ,PARAMETER  :: NX=101
      REAL(8)       ,PARAMETER  :: X1=-1.D0
      REAL(8)       ,PARAMETER  :: X2=+1.D0
      INTEGER(4)                :: I
      REAL(8)                   :: PHI
      REAL(8)                   :: FAC
      REAL(8)                   :: X
!     ******************************************************************
      POT%ID=ID
      POT%NX=NX
      POT%X1=X1
      POT%DX=(X2-X1)/DBLE(POT%NX-1) 
      ALLOCATE(POT%VAL(POT%NX))
      ALLOCATE(POT%DER(POT%NX))
!
!     ==========================================================================
!     == CALCULATE POTENTIAL                                                  ==
!     ==========================================================================
      FAC=COS(NJK*PHI0)   ! CAN BE EITHER 1 OR -1
      DO I=1,POT%NX
        X=POT%X1+DBLE(I-1)*POT%DX
        X=MIN(1.D0,X)
        PHI=ACOS(X)
        POT%VAL(I)=0.5D0*V*(1.D0-FAC*COS(DBLE(NJK)*PHI))
        IF(I.NE.1.AND.I.NE.POT%NX) THEN
          POT%DER(I)=-0.5D0*V*FAC*DBLE(NJK)*SIN(DBLE(NJK)*PHI)/SIN(PHI)
        ELSE   ! AVOID DIVIDE BY ZERO AT THE END POINTS
          POT%DER(I)=-0.5D0*V*FAC*DBLE(NJK)**2*COS(DBLE(NJK)*PHI)/COS(PHI)
        END IF
      ENDDO
      RETURN
    END SUBROUTINE CLASSICAL_TORSIONPOTA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_INVERSIONPOTA(ID,GAMMA0,K,POT)
!     ******************************************************************
!     **                                                              **
!     **  INVERSION POTENTIAL                                         **
!     **                                                              **
!     **  GAMMA0 IS THE COSINUS OF THE EQUILIBRIUM ANGLE BETWEEN THE  **
!     **  IL BOND AND THE NORMAL OF THE IJK PLANE.                    **
!     **                                                              **
!     **  K IS THE FORCE CONSTANT D2E/DGAMMA2                         **
!     **                                                              **
!     **  THE PARAMETER IS THE COSINUS OF THE ANGLE BETWEEN THE       **
!     **  IL BOND AND THE NORMAL OF THE IJK PLANE.                    **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      USE CLASSICAL_MODULE, ONLY: POT_TYPE
      IMPLICIT NONE
      CHARACTER(*)  ,INTENT(IN)  :: ID
      REAL(8)       ,INTENT(IN)  :: GAMMA0  ! OPTIMUM ANGLE
      REAL(8)       ,INTENT(IN)  :: K     !FORCE CONSTANT
      TYPE(POT_TYPE),INTENT(OUT) :: POT
      INTEGER(4)    ,PARAMETER   :: NX=20
      REAL(8)       ,PARAMETER   :: X1=-1.D0
      REAL(8)       ,PARAMETER   :: X2=1.D0
      REAL(8)       ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)                 :: I
      REAL(8)                    :: X0,X,FAC
      REAL(8)                    :: C0,C1,C2
!     ******************************************************************
      POT%ID=ID
      POT%NX=NX
      POT%X1=X1
      POT%DX=(X2-X1)/DBLE(POT%NX-1) 
      ALLOCATE(POT%VAL(POT%NX))
      ALLOCATE(POT%DER(POT%NX))
!
!     ==================================================================
!     ==================================================================
      X0=COS(GAMMA0)
      IF(1.D0-X0**2.GT.1.D-5) THEN
        FAC=0.5D0*K/(1.D0-X0**2)  !SCALE UP TO PROPER FORCE CONSTANT
        C0=FAC*X0**2
        C1=-2.D0*FAC*X0
        C2=FAC
      ELSE
        FAC=K/X0
        C0=FAC
        C1=-FAC
        C2=0
      END IF
!
!     ==================================================================
!     ==================================================================
      DO I=1,POT%NX
        X=POT%X1+DBLE(I-1)*POT%DX    !X=COS(GAMMA)
        POT%VAL(I)=C0+C1*X+C2*X**2
        POT%DER(I)=C1+2.D0*C2*X
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_NONBONDPOTA(ID,RIJ,DIJ,POT)
!     **************************************************************************
!     **                                                                      **
!     **  VAN DER WAALS POTENTIAL AS FUNCTION OF 1/R                          **
!     **                                                                      **
!     **************************************************************************
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
      INTEGER(4)                :: I
      REAL(8)                   :: SVAR,SVAR6,SVAR1,SVAR2
      REAL(8)                   :: R,X,DRDX
      LOGICAL(4)                :: TCHK
!     ******************************************************************
!     ==  WE USES A LENNARD-JONES 6-12 POTENTIAL =======================
      POT%ID=ID
      POT%NX=NX
      POT%X1=X1
      POT%DX=(X2-X1)/DBLE(POT%NX-1)
      ALLOCATE(POT%VAL(POT%NX))
      ALLOCATE(POT%DER(POT%NX))
      DO I=1,POT%NX
        X=POT%X1+POT%DX*DBLE(I-1)
        IF(I.NE.1) THEN
          R=1/X
          DRDX=-R**2
          IF(TLENNARDJONES) THEN
            SVAR=RIJ*X
            SVAR6=SVAR**6
            POT%VAL(I)=DIJ*(-2.D0*SVAR6+SVAR6**2)
            POT%DER(I)=DIJ*12.D0*(-SVAR6+SVAR6**2)*R
          ELSE
            CALL ERROR$MSG('INVALID NONBOND POTENTIAL')
            CALL ERROR$STOP('CLASSICAL_NONBONDPOTA')
          END IF
        ELSE
          POT%VAL(I)=0.D0
          POT%DER(I)=0.D0
        END IF
      ENDDO
!
!     ==========================================================================
!     == CHOP OF THE INNER NEGATIVE DIVERGENCE OF THE POTENTIAL               ==
!     ==========================================================================
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
      RETURN
    END SUBROUTINE CLASSICAL_NONBONDPOTA
! 
!     =================================================================
!     =================================================================
!     ====                                                         ====
!     ====  GENERAL PURPOSE ROUTINES                               ====
!     ====                                                         ====
!     =================================================================
!     =================================================================
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE VALUE(POT,X,F,DFDX)
!     ******************************************************************
!     **                                                             **
!     **  OBTAIN VALUE F AND DERIVATIVE DFDX OF A FUNCTION           **
!     **  WHOSE VALUE ARRAY(1,I) AND DERIVATIVE ARRAY(2,I)           **
!     **  ARE GIVEN ON A LINEAR GRID  X(I)=X1+DX*(I-1)               **
!     **                                                             **
!     ******************************************************************
      USE CLASSICAL_MODULE, ONLY: POT_TYPE
      IMPLICIT NONE
      TYPE(POT_TYPE),INTENT(IN) :: POT
      REAL(8)       ,INTENT(IN) :: X
      REAL(8)       ,INTENT(OUT):: F
      REAL(8)       ,INTENT(OUT):: DFDX
      REAL(8)                   :: D2FDX
      REAL(8)                   :: XI
      INTEGER(4)                :: IX
      REAL(8)                   :: F1V,F1D,F2V,F2D,DY,SVAR1
      INTEGER(4)                :: NX
      REAL(8)                   :: X1 ! FIRST GRID POINT
      REAL(8)                   :: XN ! LAST GRID POINT
!     ******************************************************************
!
!     ==========================================================================
!     == LINEARLY EXTRAPOLATE IF POINT FALLS OUT OF THE GRID RANGE            ==
!     ==========================================================================
      NX=POT%NX
      X1=POT%X1
      XN=POT%X1+REAL(NX-1,KIND=8)*POT%DX
!     == LINEAR EXTRAPOLATION ON THE LEFT SIDE
      IF(X.LT.X1) THEN
        F=POT%VAL(1)+(X-X1)*POT%DER(1)        
        DFDX=POT%DER(1)
        RETURN
!     == LINEAR EXTRAPOLATION ON THE RIGHT SIDE
      ELSE IF(X.GT.XN) THEN
        F=POT%VAL(NX)+(X-XN)*POT%DER(NX)        
        DFDX=POT%DER(NX)
        RETURN
      END IF
!
!     ==========================================================================
!     == POINT IS WITHIN RANGE:  SPLINE INTERPOLATION                         ==
!     ==========================================================================
      XI=(X-POT%X1)/POT%DX+1.D0
      IX=INT(XI)
      IX=MAX(1,IX)
      IX=MIN(IX,POT%NX-1)
      F1V=POT%VAL(IX)
      F2V=POT%VAL(IX+1)
      F1D=POT%DER(IX)*POT%DX
      F2D=POT%DER(IX+1)*POT%DX
      XI=XI-REAL(IX,KIND=8)
      DY=1.D0-XI
!
!     == STRAIGHT LINE THROUGH END POINTS ======================================
      F=DY*F1V+XI*F2V
      DFDX=F2V-F1V
      F1D=F1D-DFDX
      F2D=F2D-DFDX
! 
!     == SECOND ORDER ==========================================================
      D2FDX=F2D-F1D
      F=F-0.5D0*XI*DY*D2FDX
      DFDX=DFDX-0.5D0*(DY-XI)*D2FDX
      F1D=F1D+0.5D0*D2FDX
      F2D=F2D-0.5D0*D2FDX
!
!     == THIRD ORDER ===================================================
      SVAR1=-2.D0*F1D
      F=F+XI*DY*(XI-0.5D0)*SVAR1
      DFDX=DFDX+(-3.D0*XI*(XI-1)-0.5D0)*SVAR1
!
!     ==   DF/DX=DF/DI / DX/DI ========================================
      DFDX=DFDX/POT%DX
      RETURN
      END 
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE VALUE3(POT,X,F,DFDX,D2FDX2,D3FDX3)
!     ******************************************************************
!     **                                                             **
!     **  OBTAIN VALUE F AND FIRST 3 DERIVATIVES OF A FUNCTION       **
!     **  WHOSE VALUE ARRAY(1,I) AND DERIVATIVE ARRAY(2,I)           **
!     **  ARE GIVEN ON A LINEAR GRID  X(I)=X1+DX*(I-1)               **
!     **                                                             **
!     ******************************************************************
      USE CLASSICAL_MODULE, ONLY: POT_TYPE
      IMPLICIT NONE
      TYPE(POT_TYPE),INTENT(IN) :: POT
      REAL(8)       ,INTENT(IN) :: X
      REAL(8)       ,INTENT(OUT):: F
      REAL(8)       ,INTENT(OUT):: DFDX
      REAL(8)       ,INTENT(OUT):: D2FDX2
      REAL(8)       ,INTENT(OUT):: D3FDX3
      REAL(8)                   :: XI
      INTEGER(4)                :: IX
      REAL(8)                   :: F1V,F1D,F2V,F2D,DY,SVAR1
      INTEGER(4)                :: NX
      REAL(8)                   :: X1 ! FIRST GRID POINT
      REAL(8)                   :: XN ! LAST GRID POINT
!     ******************************************************************
CALL ERROR$MSG('ROUTINE NEEDS TO BE CHECKED')
CALL ERROR$STOP('VALUE3')
!
!     ===================================================================
!     == LINEARLY EXTRAPOLATE IF POINT FALLS OUT OF THE GRID RANGE     ==
!     ===================================================================
      NX=POT%NX
      X1=POT%X1
      XN=POT%X1+REAL(NX-1,KIND=8)*POT%DX
!     == LINEAR EXTRAPOLATION ON THE LEFT SIDE
      IF(X.LT.X1) THEN
        F=POT%VAL(1)+(X-X1)*POT%DER(1)        
        DFDX=POT%DER(1)
        RETURN
!     == LINEAR EXTRAPOLATION ON THE RIGHT SIDE
      ELSE IF(X.GT.XN) THEN
        F=POT%VAL(NX)+(X-XN)*POT%DER(NX)        
        DFDX=POT%DER(NX)
        RETURN
      END IF
!
!     ===================================================================
!     == POINT IS WITHIN RANGE:  SPLINE INTERPOLATION                  ==
!     ===================================================================
      XI=(X-POT%X1)/POT%DX+1.D0
      IX=INT(XI)
      IX=MAX(1,IX)
      IX=MIN(IX,POT%NX-1)
      F1V=POT%VAL(IX)
      F2V=POT%VAL(IX+1)
      F1D=POT%DER(IX)*POT%DX
      F2D=POT%DER(IX+1)*POT%DX
      XI=XI-REAL(IX,KIND=8)
      DY=1.D0-XI
!
!     == STRAIGHT LINE THROUGH END POINTS =============================
      F=DY*F1V+XI*F2V
      DFDX=F2V-F1V
      F1D=F1D-DFDX
      F2D=F2D-DFDX
! 
!     == SECOND ORDER =================================================
      D2FDX2=F2D-F1D
      F=F-0.5D0*XI*DY*D2FDX2
      DFDX=DFDX-0.5D0*(DY-XI)*D2FDX2
      F1D=F1D+0.5D0*D2FDX2
      F2D=F2D-0.5D0*D2FDX2
!
!     == THIRD ORDER ===================================================
      SVAR1=-2.D0*F1D
      F=F+XI*DY*(XI-0.5D0)*SVAR1
      DFDX=DFDX+(-3.D0*XI*(XI-1)-0.5D0)*SVAR1
      D2FDX2=D2FDX2+(-6.D0*XI+3.D0)*SVAR1
      D3FDX3=-6.D0*SVAR1
!
!     ==   DF/DX=DF/DI / DX/DI ========================================
      DFDX=DFDX/POT%DX
      D2FDX2=D2FDX2/POT%DX**2
      D3FDX3=D3FDX3/POT%DX**3
      RETURN
      END 
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE UFFTABLE_MODULE
INTEGER(4),PARAMETER :: NTYPEX=145
TYPE UFFPAR_TYPE
  CHARACTER(8) :: NAME
  REAL(8)      :: R        ! COVALENT BOND LENGTH
  REAL(8)      :: THETA    ! COVALENT BOND ANGLE
  REAL(8)      :: X        ! VDW-RADIUS
  REAL(8)      :: D        ! VDW-ENERGY
  REAL(8)      :: ZETA     ! VDW-SCALE
  REAL(8)      :: Z        ! EFFECTIVE CHARGE
END TYPE UFFPAR_TYPE
LOGICAL(4)        :: TINI=.FALSE.
TYPE(UFFPAR_TYPE) :: PAR(NTYPEX)
INTEGER(4)        :: NTYPE
CONTAINS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE UFFTABLE_INI
!     **************************************************************************
!     **  TABLE I FROM RAPPE, ET AL. JACS 114, 10024 (1992)                   **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)     :: SVAR
      REAL(8)     :: ANGSTROM
      REAL(8)     :: KCALBYMOL
!     **************************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
!     ==  (NAME,RCOV,BOND ANGLE,VDW RADIUS, VDW ENERGY, VDW-SCALE, CHARGE) =====
PAR(  1)=UFFPAR_TYPE('H_   ',0.354,180.000,2.886,0.044,12.000,0.712)
PAR(  2)=UFFPAR_TYPE('H___B',0.460, 83.500,2.886,0.044,12.000,0.712)
PAR(  3)=UFFPAR_TYPE('HE4+4',0.849, 90.000,2.362,0.056,12.000,0.098)
PAR(  4)=UFFPAR_TYPE('LI   ',1.336,180.000,2.451,0.025,15.240,1.026)
PAR(  5)=UFFPAR_TYPE('BE3+2',1.074,109.471,2.745,0.085,12.000,1.565)
PAR(  6)=UFFPAR_TYPE('B_3  ',0.838,109.471,4.083,0.180,12.052,1.755)
PAR(  7)=UFFPAR_TYPE('B_2  ',0.828,120.000,4.083,0.180,12.052,1.755)
PAR(  8)=UFFPAR_TYPE('C_3  ',0.757,109.471,3.851,0.105,12.730,1.912)
PAR(  9)=UFFPAR_TYPE('C_R  ',0.729,120.000,3.851,0.105,12.730,1.912)
PAR( 10)=UFFPAR_TYPE('C_2  ',0.732,120.000,3.851,0.105,12.730,1.912)
PAR( 11)=UFFPAR_TYPE('C_1  ',0.706,180.000,3.851,0.105,12.730,1.912)
PAR( 12)=UFFPAR_TYPE('N_3  ',0.700,106.700,3.660,0.069,13.407,2.544)
PAR( 13)=UFFPAR_TYPE('N_R  ',0.699,120.000,3.660,0.069,13.407,2.544)
PAR( 14)=UFFPAR_TYPE('N_2  ',0.685,111.300,3.660,0.069,13.407,2.544)
PAR( 15)=UFFPAR_TYPE('N_1  ',0.656,180.000,3.660,0.069,13.407,2.544)
PAR( 16)=UFFPAR_TYPE('O_3  ',0.658,104.510,3.500,0.060,14.085,2.300)
PAR( 17)=UFFPAR_TYPE('O_3_Z',0.528,145.500,3.500,0.060,14.085,2.300)
PAR( 18)=UFFPAR_TYPE('O_R  ',0.680,110.300,3.500,0.060,14.085,2.300)
PAR( 19)=UFFPAR_TYPE('O_2  ',0.634,120.000,3.500,0.060,14.085,2.300)
PAR( 20)=UFFPAR_TYPE('O_1  ',0.639,180.000,3.500,0.060,14.085,2.300)
PAR( 21)=UFFPAR_TYPE('F_   ',0.668,180.000,3.364,0.050,14.762,1.735)
PAR( 22)=UFFPAR_TYPE('NE4+4',0.920, 90.000,3.243,0.042,15.440,0.194)
PAR( 23)=UFFPAR_TYPE('NA   ',1.539,180.000,2.983,0.030,12.000,1.081)
PAR( 24)=UFFPAR_TYPE('MG3+2',1.421,109.471,3.021,0.111,12.000,1.787)
PAR( 25)=UFFPAR_TYPE('AL3  ',1.244,109.471,4.499,0.505,12.000,1.792)
PAR( 26)=UFFPAR_TYPE('SI3  ',1.117,109.471,4.295,0.402,12.000,2.323)
PAR( 27)=UFFPAR_TYPE('P_3+3',1.101, 93.800,4.147,0.305,12.000,2.863)
PAR( 28)=UFFPAR_TYPE('P_3+5',1.056,109.471,4.147,0.305,12.000,2.863)
PAR( 29)=UFFPAR_TYPE('P_3+Q',1.056,109.471,4.147,0.305,12.000,2.863)
PAR( 30)=UFFPAR_TYPE('S_3+2',1.064, 92.100,4.035,0.274,12.000,2.703)
PAR( 31)=UFFPAR_TYPE('S_3+4',1.049,103.200,4.035,0.274,12.000,2.703)
PAR( 32)=UFFPAR_TYPE('S_3+6',1.027,109.471,4.035,0.274,12.000,2.703)
PAR( 33)=UFFPAR_TYPE('S_R  ',1.077, 92.200,4.035,0.274,12.000,2.703)
PAR( 34)=UFFPAR_TYPE('S_2  ',0.854,120.000,4.035,0.274,12.000,2.703)
PAR( 35)=UFFPAR_TYPE('CL   ',1.044,180.000,3.947,0.227,12.000,2.348)
PAR( 36)=UFFPAR_TYPE('AR4+4',1.032, 90.000,3.868,0.185,12.000,0.300)
PAR( 37)=UFFPAR_TYPE('K_   ',1.953,180.000,3.812,0.035,12.000,1.165)
PAR( 38)=UFFPAR_TYPE('CA6+2',1.761, 90.000,3.399,0.238,12.000,2.141)
PAR( 39)=UFFPAR_TYPE('SC3+3',1.513,109.471,3.295,0.019,12.000,2.592)
PAR( 40)=UFFPAR_TYPE('TI3+4',1.412,109.471,3.175,0.017,12.000,2.659)
PAR( 41)=UFFPAR_TYPE('TI6+4',1.412, 90.000,3.175,0.017,12.000,2.659)
PAR( 42)=UFFPAR_TYPE('V_3+5',1.402,109.471,3.144,0.016,12.000,2.679)
PAR( 43)=UFFPAR_TYPE('CR6+3',1.345, 90.000,3.023,0.015,12.000,2.463)
PAR( 44)=UFFPAR_TYPE('MN6+2',1.382, 90.000,2.961,0.013,12.000,2.430)
PAR( 45)=UFFPAR_TYPE('FE3+2',1.412,109.470,2.912,0.013,12.000,2.430)
PAR( 46)=UFFPAR_TYPE('FE6+2',1.335, 90.000,2.912,0.013,12.000,2.430)
PAR( 47)=UFFPAR_TYPE('CO6+3',1.241, 90.000,2.872,0.014,12.000,2.430)
PAR( 48)=UFFPAR_TYPE('NI4+2',1.164, 90.000,2.834,0.015,12.000,2.430)
PAR( 49)=UFFPAR_TYPE('CU3+1',1.302,109.471,3.495,0.005,12.000,1.756)
PAR( 50)=UFFPAR_TYPE('ZN3+2',1.193,109.471,2.763,0.124,12.000,1.308)
PAR( 51)=UFFPAR_TYPE('GA3+3',1.260,109.471,4.383,0.415,12.000,1.821)
PAR( 52)=UFFPAR_TYPE('GE3  ',1.197,109.471,4.280,0.379,12.000,2.789)
PAR( 53)=UFFPAR_TYPE('AS3+3',1.211, 92.100,4.230,0.309,12.000,2.864)
PAR( 54)=UFFPAR_TYPE('SE3+2',1.190, 90.600,4.205,0.291,12.000,2.764)
PAR( 55)=UFFPAR_TYPE('BR   ',1.192,180.000,4.189,0.251,12.000,2.519)
PAR( 56)=UFFPAR_TYPE('KR4+4',1.147, 90.000,4.141,0.220,12.000,0.452)
PAR( 57)=UFFPAR_TYPE('RB   ',2.260,180.000,4.114,0.040,12.000,1.592)
PAR( 58)=UFFPAR_TYPE('SR6+2',2.052, 90.000,3.641,0.235,12.000,2.449)
PAR( 59)=UFFPAR_TYPE('Y_3+3',1.698,109.471,3.345,0.072,12.000,3.257)
PAR( 60)=UFFPAR_TYPE('ZR3+4',1.564,109.471,3.124,0.069,12.000,3.667)
PAR( 61)=UFFPAR_TYPE('NB3+5',1.473,109.471,3.165,0.059,12.000,3.618)
PAR( 62)=UFFPAR_TYPE('MO6+6',1.467, 90.000,3.052,0.056,12.000,3.400)
PAR( 63)=UFFPAR_TYPE('MO3+6',1.484,109.471,3.052,0.056,12.000,3.400)
PAR( 64)=UFFPAR_TYPE('TC6+5',1.322, 90.000,2.998,0.048,12.000,3.400)
PAR( 65)=UFFPAR_TYPE('RU6+2',1.478, 90.000,2.963,0.056,12.000,3.400)
PAR( 66)=UFFPAR_TYPE('RH6+3',1.332, 90.000,2.929,0.053,12.000,3.508)
PAR( 67)=UFFPAR_TYPE('PD4+2',1.338, 90.000,2.899,0.048,12.000,3.210)
PAR( 68)=UFFPAR_TYPE('AG1+1',1.386,180.000,3.148,0.036,12.000,1.956)
PAR( 69)=UFFPAR_TYPE('CD3+2',1.403,109.471,2.848,0.228,12.000,1.650)
PAR( 70)=UFFPAR_TYPE('IN3+3',1.459,109.471,4.463,0.599,12.000,2.070)
PAR( 71)=UFFPAR_TYPE('SN3  ',1.398,109.471,4.392,0.567,12.000,2.961)
PAR( 72)=UFFPAR_TYPE('SB3+3',1.407, 91.600,4.420,0.449,12.000,2.704)
PAR( 73)=UFFPAR_TYPE('TE3+2',1.386, 90.250,4.470,0.398,12.000,2.882)
PAR( 74)=UFFPAR_TYPE('I_   ',1.382,180.000,4.500,0.339,12.000,2.650)
PAR( 75)=UFFPAR_TYPE('XE4+4',1.267, 90.000,4.404,0.332,12.000,0.556)
PAR( 76)=UFFPAR_TYPE('CS   ',2.570,180.000,4.517,0.045,12.000,1.573)
PAR( 77)=UFFPAR_TYPE('BA6+2',2.277, 90.000,3.703,0.364,12.000,2.727)
PAR( 78)=UFFPAR_TYPE('LA3+3',1.943,109.471,3.522,0.017,12.000,3.300)
PAR( 79)=UFFPAR_TYPE('CE6+3',1.841, 90.000,3.556,0.013,12.000,3.300)
PAR( 80)=UFFPAR_TYPE('PR6+3',1.823, 90.000,3.606,0.010,12.000,3.300)
PAR( 81)=UFFPAR_TYPE('ND6+3',1.816, 90.000,3.575,0.010,12.000,3.300)
PAR( 82)=UFFPAR_TYPE('PM6+3',1.801, 90.000,3.547,0.009,12.000,3.300)
PAR( 83)=UFFPAR_TYPE('SM6+3',1.780, 90.000,3.520,0.008,12.000,3.300)
PAR( 84)=UFFPAR_TYPE('EU6+3',1.771, 90.000,3.493,0.008,12.000,3.300)
PAR( 85)=UFFPAR_TYPE('GD6+3',1.735, 90.000,3.368,0.009,12.000,3.300)
PAR( 86)=UFFPAR_TYPE('TB6+3',1.732, 90.000,3.451,0.007,12.000,3.300)
PAR( 87)=UFFPAR_TYPE('DY6+3',1.710, 90.000,3.428,0.007,12.000,3.300)
PAR( 88)=UFFPAR_TYPE('HO6+3',1.696, 90.000,3.409,0.007,12.000,3.416)
PAR( 89)=UFFPAR_TYPE('ER6+3',1.673, 90.000,3.391,0.007,12.000,3.300)
PAR( 90)=UFFPAR_TYPE('TM6+3',1.660, 90.000,3.374,0.006,12.000,3.300)
PAR( 91)=UFFPAR_TYPE('YB6+3',1.637, 90.000,3.355,0.228,12.000,2.618)
PAR( 92)=UFFPAR_TYPE('LU6+3',1.671, 90.000,3.640,0.041,12.000,3.271)
PAR( 93)=UFFPAR_TYPE('HF3+4',1.611,109.471,3.141,0.072,12.000,3.921)
PAR( 94)=UFFPAR_TYPE('TA3+5',1.511,109.471,3.170,0.081,12.000,4.075)
PAR( 95)=UFFPAR_TYPE('W_6+6',1.392, 90.000,3.069,0.067,12.000,3.700)
PAR( 96)=UFFPAR_TYPE('W_3+4',1.526,109.471,3.069,0.067,12.000,3.700)
PAR( 97)=UFFPAR_TYPE('W_3+6',1.380,109.471,3.069,0.067,12.000,3.700)
PAR( 98)=UFFPAR_TYPE('RE6+5',1.372, 90.000,2.954,0.066,12.000,3.700)
PAR( 99)=UFFPAR_TYPE('RE3+7',1.314,109.471,2.954,0.066,12.000,3.700)
PAR(100)=UFFPAR_TYPE('OS6+6',1.372, 90.000,3.120,0.037,12.000,3.700)
PAR(101)=UFFPAR_TYPE('IR6+3',1.371, 90.000,2.840,0.073,12.000,3.731)
PAR(102)=UFFPAR_TYPE('PT4+2',1.364, 90.000,2.754,0.080,12.000,3.382)
PAR(103)=UFFPAR_TYPE('AU4+3',1.262, 90.000,3.293,0.039,12.000,2.625)
PAR(104)=UFFPAR_TYPE('HG1+2',1.340,180.000,2.705,0.385,12.000,1.750)
PAR(105)=UFFPAR_TYPE('TL3+3',1.518,120.000,4.347,0.680,12.000,2.068)
PAR(106)=UFFPAR_TYPE('PB3  ',1.459,109.471,4.297,0.663,12.000,2.846)
PAR(107)=UFFPAR_TYPE('BI3+3',1.512, 90.000,4.370,0.518,12.000,2.470)
PAR(108)=UFFPAR_TYPE('PO3+2',1.500, 90.000,4.709,0.325,12.000,2.330)
PAR(119)=UFFPAR_TYPE('AT   ',1.545,180.000,4.750,0.284,12.000,2.240)
PAR(120)=UFFPAR_TYPE('RN4+4',1.420, 90.000,4.765,0.248,12.000,0.583)
PAR(121)=UFFPAR_TYPE('FR   ',2.880,180.000,4.900,0.050,12.000,1.847)
PAR(122)=UFFPAR_TYPE('RA6+2',2.512, 90.000,3.677,0.404,12.000,2.920)
PAR(123)=UFFPAR_TYPE('AC6+3',1.983, 90.000,3.478,0.033,12.000,3.900)
PAR(124)=UFFPAR_TYPE('TH6+4',1.721, 90.000,3.396,0.026,12.000,4.202)
PAR(125)=UFFPAR_TYPE('PA6+4',1.711, 90.000,3.424,0.022,12.000,3.900)
PAR(126)=UFFPAR_TYPE('U_6+4',1.684, 90.000,3.395,0.022,12.000,3.900)
PAR(127)=UFFPAR_TYPE('NP6+4',1.666, 90.000,3.424,0.019,12.000,3.900)
PAR(128)=UFFPAR_TYPE('PU6+4',1.657, 90.000,3.424,0.016,12.000,3.900)
PAR(129)=UFFPAR_TYPE('AM6+4',1.660, 90.000,3.381,0.014,12.000,3.900)
PAR(130)=UFFPAR_TYPE('CM6+3',1.801, 90.000,3.326,0.013,12.000,3.900)
PAR(131)=UFFPAR_TYPE('BK6+3',1.761, 90.000,3.339,0.013,12.000,3.900)
PAR(132)=UFFPAR_TYPE('CF6+3',1.750, 90.000,3.313,0.013,12.000,3.900)
PAR(133)=UFFPAR_TYPE('ES6+3',1.724, 90.000,3.299,0.012,12.000,3.900)
PAR(134)=UFFPAR_TYPE('FM6+3',1.712, 90.000,3.286,0.012,12.000,3.900)
PAR(135)=UFFPAR_TYPE('MD6+3',1.689, 90.000,3.274,0.011,12.000,3.900)
PAR(136)=UFFPAR_TYPE('NO6+3',1.679, 90.000,3.248,0.011,12.000,3.900)
PAR(137)=UFFPAR_TYPE('LW6+3',1.698, 90.000,3.236,0.011,12.000,3.900)
! CYCLOPENTADIENYL CENTER
PAR(138)=UFFPAR_TYPE('CPR  ',0.551, 90.000,3.851,0.000,12.000,1.912)
! CYCLOPENTADIENYL CENTER WITH BACK BONDING SUCH AS FERROCENE
PAR(139)=UFFPAR_TYPE('CPR_B',0.340, 90.000,3.851,0.000,12.000,1.912)
! DOUBLE BOND CENTER
PAR(140)=UFFPAR_TYPE('CIR  ',0.616, 90.000,3.851,0.000,12.000,1.912)
! ALLYL CENTER
PAR(141)=UFFPAR_TYPE('PIR  ',0.616, 90.000,3.851,0.000,12.000,1.912)
!
!     ==========================================================================
!     ==  TIP4P WATER MODEL INHERITS O_3,H_,H_                                ==
!     ==========================================================================
!     == 1) THE EFFECTIVE CHARGES DETERMINE THE FORCE CONSTANTS FOR BOND-     ==
!     ==    AND ANGLE FORCES. THEY ARE SET TO ZERO BECAUSE ALL INTERNAL       ==
!     ==    MOLECULAR FORCES ARE ACCOUNTED FOR BY CONSTRAINTS.                ==
!        2) THE LENNARD JONES POTENTIAL IS ATTACHED ONLY TO THE O-ATOM        ==
!     == 3) M_TIP4P IS A DUMMY ATOM, CARRIES A CHARGE                         ==
      PAR(142)=UFFPAR_TYPE('O_TIP4P',0.658,104.510,3.500,0.060,14.085,0.0)
      PAR(143)=UFFPAR_TYPE('H_TIP4P',0.354,180.000,2.886,0.000,12.000,0.0)
      PAR(144)=UFFPAR_TYPE('M_TIP4P',0.354,180.000,2.886,0.000,12.000,0.0)
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
      CALL CONSTANTS$GET('KCAL/MOL',KCALBYMOL)
      CALL TIP4P$GETR8('RLJ',SVAR) ; SVAR=SVAR/ANGSTROM
      PAR(142)%X=SVAR
      CALL TIP4P$GETR8('DLJ',SVAR) ; SVAR=SVAR/KCALBYMOL
      PAR(142)%D=SVAR
      CALL TIP4P$GETR8('ROH',SVAR) ; SVAR=SVAR/ANGSTROM
      PAR(143)%R=SVAR-PAR(142)%R   !R_COV,H = DISTANCE - R_COV,O
!
      NTYPE=144
      TINI=.TRUE.
      RETURN
      END SUBROUTINE UFFTABLE_INI
END MODULE UFFTABLE_MODULE
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE UFFTABLE$EXIST(TYPE_,EXIST,SUGGESTION)
!     **************************************************************************
!     **  CHECK IF FORCE FIELD TYPE EXIST.                                    **
!     **  USED IN PAW_PREOPT TOOL                                             **
!     **                                                                      **
!     **  REMARK                                                              **
!     **  IF THE STRING TYPE IS CONTAINED IN AN EXISTING ATOM TYPE            **
!     **  BUT NOT IDENTICAL TO AN EXISTING ATOM TYPE                          **
!     **  ONE OF THE EXISTING ATOM TYPES IS RETURNED AS SUGGESTION            **
!     **                                                                      **
!     ************************WRITTEN BY JOHANNES KAESTNER 2002*****************
      USE UFFTABLE_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: TYPE_
      LOGICAL(4)  ,INTENT(OUT):: EXIST
      CHARACTER(*),INTENT(OUT):: SUGGESTION
      INTEGER(4)              :: ITYPE
!     **************************************************************************
      CALL UFFTABLE_INI
!
!     ==========================================================================
!     ==  FIND ATOM TYPE                                                      ==
!     ==========================================================================
      ITYPE=1
      SUGGESTION=' '
      EXIST=.FALSE.
      DO ITYPE=1,NTYPE
         IF(TYPE_.EQ.PAR(ITYPE)%NAME) THEN
            EXIST=.TRUE.
            SUGGESTION=PAR(ITYPE)%NAME
            EXIT
         ELSE
            IF(INDEX(PAR(ITYPE)%NAME,TYPE_).GT.0) SUGGESTION=PAR(ITYPE)%NAME
         END IF
      END DO
    END SUBROUTINE UFFTABLE$EXIST
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE UFFTABLE$GET(TYPE_,ID,VAL_)
!     **************************************************************************
!     **************************************************************************
      USE UFFTABLE_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: TYPE_
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL_
      INTEGER(4)              :: ITYPE
!     ******************************************************************
      CALL UFFTABLE_INI
!
!     ==================================================================
!     ==  FIND ATOM TYPE                                              ==
!     ==================================================================
      ITYPE=1
      DO WHILE(TYPE_.NE.PAR(ITYPE)%NAME) 
        ITYPE=ITYPE+1
        IF(ITYPE.GT.NTYPE) THEN
          CALL ERROR$MSG('FFTYPE NOT RECOGNIZED')
          CALL ERROR$CHVAL('TYPE ',TYPE_)
          CALL ERROR$STOP('UFFTABLE$GET')
        END IF
      ENDDO
!
!     ==================================================================
!     ==  RESOLVE ID                                                  ==
!     ==================================================================
      IF(ID.EQ.'R') THEN          ! BOND LENGTH
        VAL_=PAR(ITYPE)%R
      ELSE IF(ID.EQ.'THETA') THEN ! COVALENT ANGLE
        VAL_=PAR(ITYPE)%THETA
      ELSE IF(ID.EQ.'X') THEN     ! VAN DER WAALS RADIUS
        VAL_=PAR(ITYPE)%X
      ELSE IF(ID.EQ.'D') THEN     ! VAN DER WAALS ENERGY
        VAL_=PAR(ITYPE)%D
      ELSE IF(ID.EQ.'ZETA') THEN
        VAL_=PAR(ITYPE)%ZETA
      ELSE IF(ID.EQ.'Z') THEN      ! EFFECTIVE CHARGE (FOR FORCE CONSTANTS)
        VAL_=PAR(ITYPE)%Z
      ELSE IF(ID.EQ.'ENC') THEN    !PAULING ELECTRONEGATIVITY
        CALL PERIODICTABLE$GET(TYPE_(1:2),'EN',VAL_)
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$STOP('UFFTABLE$GET')
      END IF
      RETURN
      END        
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE UFFTABLE_BONDLENGTH(TYPE1,TYPE2,BO,R)
!     ******************************************************************
!     **                                                              **
!     **  EQULIBRIUM COVALENT BOND LENGTH                             **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: TYPE1      ! UFF ATOM TYPE FIRST ATOM
      CHARACTER(*),INTENT(IN) :: TYPE2      ! UFF ATOM TYPE SECOND ATOM
      REAL(8)     ,INTENT(IN) :: BO         ! BOND ORDER
      REAL(8)     ,INTENT(OUT):: R          ! EQUILIBRIUM BOND LENGTH
      REAL(8)                 :: RI,RJ      ! COVALENT RADII
      REAL(8)                 :: ENCI,ENCJ  ! ELECTRONEGATIVITIES
      REAL(8)                 :: ZI,ZJ      ! EFFECTIVE BOND CHARGES
      REAL(8)                 :: RBO        ! BOND ORDER CORRECTION
      REAL(8)                 :: REN        ! ELECTRONEGATIVITY CORRECTION
!     ******************************************************************
!
!     ==================================================================
!     == COLLECT UFF DATA                                            ==
!     ==================================================================
      CALL UFFTABLE$GET(TYPE1,'R',RI)
      CALL UFFTABLE$GET(TYPE2,'R',RJ)
      CALL UFFTABLE$GET(TYPE1,'ENC',ENCI)
      CALL UFFTABLE$GET(TYPE2,'ENC',ENCJ)
      CALL UFFTABLE$GET(TYPE1,'Z',ZI)
      CALL UFFTABLE$GET(TYPE2,'Z',ZJ)
!
!     ==================================================================
!     == BOND ORDER CORRECTION : EQ.3                                 ==
!     ==================================================================
      RBO=-0.1332D0*(RI+RJ)*LOG(BO)
!
!     ==================================================================
!     == ELECTRONEGATIVITY CORRECTION EQ.4                            ==
!     == ATTENTION! WRONG SIGN IN THE UFF PAPER (CORRECTED HERE)      ==
!     == ( SEE ALSO COMMENTS FOR UFF FROM MCCCS TOWHEE)
!     ==================================================================
      REN=-RI*RJ*(SQRT(ENCI)-SQRT(ENCJ))**2/(ENCI*RI+ENCJ*RJ)
!     
      R=RI+RJ+RBO+REN               ! BOND DISTANCE
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE UFFTABLE$BONDPARMS(TYPE1,TYPE2,BO,ID,R,K,D,TCHK)
!     ******************************************************************
!     **                                                              **
!     ** UFF BONDING POTENTIAL AS FUNCTION OF 1/DISTANCE              **
!     **                                                              **
!     **   TYPE1,TYPE2 ARE MNEMOTIC NAMES FOR THE ATOMTYPE            **
!     **   BOND ORDERS ARE SINGLE BOND: BO=1; DOUBLE BOND: BO=2;      **
!     **   TRIPLE BOND:BO=3; AROMATIC RING: BO=1.5                    **
!     **                                                              **
!     **  REMARK: UFFTABLE GIVES DATA IN KCAL/MOL AND ANGSTROM.       **
!     **          CLASSICAL_BONDPOT PRODUCES DATA IN A.U.             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*)  ,INTENT(IN)  :: TYPE1 ! FF-ATOM-TYPE
      CHARACTER(*)  ,INTENT(IN)  :: TYPE2 ! FF-ATOM-TYPE
      REAL(8)       ,INTENT(IN)  :: BO    ! BOND ORDER
      CHARACTER(64) ,INTENT(OUT) :: ID    ! IDENTIFIER
      REAL(8)       ,INTENT(OUT) :: R     ! BOND DISTANCE
      REAL(8)       ,INTENT(OUT) :: K     ! FORCE CONSTANT
      REAL(8)       ,INTENT(OUT) :: D     ! BINDING ENERGY
      LOGICAL(4)    ,INTENT(OUT) :: TCHK  ! NON-ZERO FORCE FIELD
      REAL(8)                    :: ZI,ZJ
      REAL(8)                    :: ANGSTROM
      REAL(8)                    :: KCALBYMOL
      LOGICAL(4)                 :: TCHK1
!     ******************************************************************
!
!     ==================================================================
!     == THE FOLLOWING PARAMETERS DETERMINE THE BOND-POTENTIAL        ==
!     ==================================================================
      CALL UFFTABLE_BONDSPECIAL(TYPE1,TYPE2,BO,R,K,D,TCHK1)
      IF(.NOT.TCHK1) THEN
        CALL UFFTABLE_BONDLENGTH(TYPE1,TYPE2,BO,R)
        CALL UFFTABLE$GET(TYPE1,'Z',ZI)
        CALL UFFTABLE$GET(TYPE2,'Z',ZJ)
        D=BO*70.D0            ! BOND DISSOCIATION ENERGY (ONLY FOR NON-HARMONIC)
        K=2.D0*332.06D0*ZI*ZJ/R**3 ! FORCE CONSTANT
      END IF
!
!     ==================================================================
!     == WRITE POTENTIAL ID                                          ==
!     ==================================================================
      IF(LGT(TYPE1,TYPE2)) THEN
        WRITE(ID,FMT='(A1,2A8," BO=",F5.1," R=",F5.2," K=",F5.0," D=",F5.0)') &
     &        'B ',TYPE1,TYPE2,BO,R,K,D
      ELSE
        WRITE(ID,FMT='(A1,2A8," BO=",F5.1," R=",F5.2," K=",F5.0," D=",F5.0)') &
     &        'B ',TYPE2,TYPE1,BO,R,K,D
      END IF
!
!     ==================================================================
!     == CONVERT TO ATOMIC UNITS                                      ==
!     ==================================================================
      CALL CONSTANTS('KCAL/MOL',KCALBYMOL)
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      R=R*ANGSTROM
      D=D*KCALBYMOL
      K=K*KCALBYMOL/ANGSTROM**2
!
!     ==================================================================
!     == RETURN IF FORCE CONSTANT IS ZERO                             ==
!     == NOT POSSIBLE IN THIS VERSION BECAUSE BONDS ARE PREASSIGNED   ==
!     ==================================================================
      TCHK=(K.GT.1.D-6)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE UFFTABLE_BONDSPECIAL(ATOM1,ATOM2,BO,R,K,D,TCHK)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      TYPE BONDRULE_TYPE 
        CHARACTER(8) :: ATOM1
        CHARACTER(8) :: ATOM2
        REAL(8)      :: R       ! BOND LENGTH
        REAL(8)      :: K       ! FORCE CONSTANT
        REAL(8)      :: D       ! BINDING ENERGY, USED WITH MORSE
      END TYPE BONDRULE_TYPE 
      CHARACTER(*),INTENT(IN) :: ATOM1   ! FIRST ATOM TYPE
      CHARACTER(*),INTENT(IN) :: ATOM2   ! SECOND ATOM TYPE
      REAL(8)     ,INTENT(IN) :: BO      ! BOND ORDER
      REAL(8)     ,INTENT(OUT):: R       ! BOND DISTANCE
      REAL(8)     ,INTENT(OUT):: D       ! DISSOCIATION ENERGY (MORSE)
      REAL(8)     ,INTENT(OUT):: K       ! FORCE CONSTANT
      LOGICAL(4)  ,INTENT(OUT):: TCHK    ! SPECIAL RULE APPLIED 
      INTEGER(4)  ,PARAMETER  :: NRULE=4 ! # OF SPECIAL RULES
      TYPE(BONDRULE_TYPE)     :: RULE(NRULE)
      INTEGER(4)              :: I
!     ******************************************************************
      RULE(1)=BONDRULE_TYPE('CPR  ','C_R  ',1.1733D0,700.D0,0.D0)
      RULE(2)=BONDRULE_TYPE('CPR_B','C_R  ',1.1733D0,700.D0,0.D0)
      RULE(3)=BONDRULE_TYPE('CIR  ','C_2  ',0.6644D0,700.D0,0.D0)
      RULE(4)=BONDRULE_TYPE('PIR  ','C_2  ',0.9951D0,0.D0,700.D0)
      TCHK=.TRUE.
      DO I=1,NRULE
        CALL COMPARE(ATOM1,ATOM2,RULE(I),TCHK)
        IF(TCHK) THEN
          R=RULE(I)%R
          D=RULE(I)%D
          K=RULE(I)%K
          RETURN
        END IF
      ENDDO
      TCHK=.FALSE.
      RETURN
      CONTAINS
!       ................................................................
        SUBROUTINE COMPARE(ATOM1,ATOM2,RULE,TCHK)
!       ****************************************************************
        CHARACTER(*),INTENT(IN)        :: ATOM1
        CHARACTER(*),INTENT(IN)        :: ATOM2
        TYPE(BONDRULE_TYPE),INTENT(IN) :: RULE
        LOGICAL(4)         ,INTENT(OUT):: TCHK
!       ****************************************************************
        TCHK=.TRUE.
        IF(ATOM1.EQ.RULE%ATOM1) THEN
          IF(ATOM2.EQ.RULE%ATOM2) RETURN
        ELSE IF(ATOM1.EQ.RULE%ATOM2) THEN
          IF(ATOM2.EQ.RULE%ATOM1) RETURN
        END IF
        TCHK=.FALSE.
        RETURN
        END SUBROUTINE COMPARE
      END SUBROUTINE UFFTABLE_BONDSPECIAL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE UFFTABLE$ANGLEPARMS(ATOM1,ATOM2,ATOM3,BO1,BO2,ID,THETA,K,TCHK)
!     ******************************************************************
!     **                                                              **
!     ** UFF BOND BENDING POTENTIAL AS FUNCTION OF COS(PHI)           **
!     **                                                              **
!     **   ATOM1,ATOM2,ATOM3 MNEMOTIC NAME FOR THE ATOMTYPE           **
!     **   BO BOND ORDER (SINGLE BOND: BO=1; DOUBLE BOND: BO=2;       **
!     **   TRIPLE BOND:BO=3; AROMATIC RING: BO=1.5                    **
!     **   ON INPUT VARRAY CONTAINS THE ANGLE THETA(X)                **
!     **  0.664                                                       **
!     **                                                              **
!     ******************************************************************
      USE CLASSICAL_MODULE, ONLY: POT_TYPE
      IMPLICIT NONE
      CHARACTER(*)  ,INTENT(IN)  :: ATOM1 ! FF-ATOM-TYPE TERMINAL ATOM   
      CHARACTER(*)  ,INTENT(IN)  :: ATOM2 ! FF-ATOM-TYPE CENTRAL ATOM   
      CHARACTER(*)  ,INTENT(IN)  :: ATOM3 ! FF-ATOM-TYPE TERMINAL ATOM
      REAL(8)       ,INTENT(IN)  :: BO1   ! BOND ORDER (ATOM1-ATOM2)
      REAL(8)       ,INTENT(IN)  :: BO2   ! BOND ORDER (ATOM2-ATOM3)
      CHARACTER(64) ,INTENT(OUT) :: ID    ! IDENTIFIER
      REAL(8)       ,INTENT(OUT) :: THETA ! EQUILIBRIUM ANGLE
      REAL(8)       ,INTENT(OUT) :: K     ! FORCE CONSTANT
      LOGICAL(4)    ,INTENT(OUT) :: TCHK
      REAL(8)       ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      REAL(8)                    :: KCALBYMOL
      REAL(8)                    :: ZI,ZK
      REAL(8)                    :: RIJ,RJK,RIK,BETA
      LOGICAL(4)                 :: TCHK1
!     ******************************************************************
!
!     ==================================================================
!     ==  SCAN SPECIAL RULES                                          ==
!     ==================================================================
      CALL UFFTABLE_ANGLESPECIAL(ATOM1,ATOM2,ATOM3,BO1,BO2,THETA,K,TCHK1)
!
!     ==================================================================
!     ==  APPLY GENERAL RULE                                          ==
!     ==================================================================
      IF(.NOT.TCHK1) THEN
        CALL UFFTABLE$GET(ATOM1,'Z',ZI)
        CALL UFFTABLE$GET(ATOM3,'Z',ZK)
        CALL UFFTABLE$GET(ATOM2,'THETA',THETA); THETA=THETA/180.D0*PI
        CALL UFFTABLE_BONDLENGTH(ATOM1,ATOM2,BO1,RIJ)
        CALL UFFTABLE_BONDLENGTH(ATOM2,ATOM3,BO2,RJK)
        RIK=SQRT(RIJ**2+RJK**2-2.D0*RIJ*RJK*COS(THETA))
        BETA=664.12/(RIJ*RJK)
        K=BETA*ZI*ZK*RIJ*RJK/RIK**5 &
     &                      *(3.D0*RIJ*RJK*SIN(THETA)**2-RIK**2*COS(THETA))
      END IF
!
!     ==================================================================
!     ==  WRITE STRING IDENTIFYING POTENTIAL                          ==
!     ==================================================================
      IF(LGT(ATOM1,ATOM3)) THEN
        WRITE(ID,FMT='(A1,3A8," BO=",F5.1," BO=",F5.1' &
     &              //'," THETA=",F5.0," K=",F5.0)') &
     &                'A ',ATOM1,ATOM2,ATOM3,BO1,BO2,THETA/PI*180,K
      ELSE
        WRITE(ID,FMT='(A1,3A8," BO=",F5.1," BO=",F5.1' &
     &              //'," THETA=",F5.0," K=",F5.0)') &
     &                'A ',ATOM3,ATOM2,ATOM1,BO2,BO1,THETA/PI*180,K
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
      TCHK=(K.GT.1.D-6) !DO NOT EXCLUDE ANGLES HERE, BECAUSE THEY DEFINE 
                         ! BOND EXCLUSIONS
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE UFFTABLE_ANGLESPECIAL(ATOM1,ATOM2,ATOM3,BO1,BO2,THETA,K,TCHK)
!     ******************************************************************
!     **                                                              **
!     **  SPECIAL RULES MAY CONTAIN JOKERS 'X' WHICH APPLY TO ANY     **
!     **  ATOM TYPE                                                   **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      TYPE ANGLERULE_TYPE 
        CHARACTER(8) :: ATOM1
        CHARACTER(8) :: ATOM2
        CHARACTER(8) :: ATOM3
        REAL(8)      :: THETA   ! EQUILIBRIUM ANGLE
        REAL(8)      :: K       ! FORCE CONSTANT
      END TYPE ANGLERULE_TYPE 
      CHARACTER(*),INTENT(IN) :: ATOM1   ! FIRST ATOM TYPE
      CHARACTER(*),INTENT(IN) :: ATOM2   ! SECOND (CENTRAL) ATOM TYPE
      CHARACTER(*),INTENT(IN) :: ATOM3   ! THIRD ATOM TYPE
      REAL(8)     ,INTENT(IN) :: BO1     ! BOND ORDER (ATOM1,ATOM2)
      REAL(8)     ,INTENT(IN) :: BO2     ! BOND ORDER (ATOM2-ATOM3)
      REAL(8)     ,INTENT(OUT):: THETA   ! EQUILIBRIUM ANGLE
      REAL(8)     ,INTENT(OUT):: K       ! FORCE CONSTANT
      LOGICAL(4)  ,INTENT(OUT):: TCHK    ! SPECIAL RULE APPLIED 
      REAL(8)     ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)  ,PARAMETER  :: NRULE=17 ! # OF SPECIAL RULES
      TYPE(ANGLERULE_TYPE)    :: RULE(NRULE)
      INTEGER(4)              :: I
!     ******************************************************************
      RULE(1)=ANGLERULE_TYPE('C_R  ','CPR  ','X    ', 90.D0,100.D0)
      RULE(2)=ANGLERULE_TYPE('CPR  ','C_R  ','X    ',180.D0,100.D0)
      RULE(3)=ANGLERULE_TYPE('C_R  ','CPR  ','C_R  ', 72.D0,  0.D0)
      RULE(4)=ANGLERULE_TYPE('CPR  ','C_R  ','C_R  ', 54.D0,  0.D0)
      RULE(5)=ANGLERULE_TYPE('CPR_B','C_R  ','X    ',180.D0,100.D0)
      RULE(6)=ANGLERULE_TYPE('CPR_B','C_R  ','C_R  ', 54.D0,  0.D0)
      RULE(7)=ANGLERULE_TYPE('PIR  ','C_R  ','X    ',180.D0,100.D0)
      RULE(8)=ANGLERULE_TYPE('PIR  ','C_R  ','C_R  ', 54.D0,  0.D0)
      RULE(9)=ANGLERULE_TYPE('CIR  ','C_2  ','C_2  ',  0.D0,  0.D0)
      RULE(10)=ANGLERULE_TYPE('CPR_B','MN6+2','X    ',135.D0,100.D0)
      RULE(11)=ANGLERULE_TYPE('CPR_B','FE6+2','X    ',135.D0,100.D0)
      RULE(12)=ANGLERULE_TYPE('C_R  ','CPR_B','X    ', 90.D0,100.D0)
      RULE(13)=ANGLERULE_TYPE('C_R  ','CPR_B','C_R  ', 72.D0,  0.D0)
      RULE(14)=ANGLERULE_TYPE('C_2  ','CIR  ','X    ', 90.D0,100.D0)
      RULE(15)=ANGLERULE_TYPE('C_2  ','CIR  ','C_2  ',180.D0,  0.D0)
      RULE(16)=ANGLERULE_TYPE('C_2  ','PIR  ','X    ', 90.D0,100.D0)
      RULE(17)=ANGLERULE_TYPE('C_2  ','PIR  ','C_R  ', 72.D0,  0.D0)
      DO I=1,NRULE
        CALL UFFTABLE_ANGLESPECIAL_COMPARE(ATOM1,ATOM2,ATOM3,RULE(I),TCHK)
        IF(TCHK) THEN
          THETA=RULE(I)%THETA/180.D0*PI
          K=RULE(I)%K
          RETURN
        END IF
      ENDDO
      TCHK=.FALSE.
      RETURN
      CONTAINS
!
!       .1.........2.........3.........4.........5.........6.........7.........8
        SUBROUTINE UFFTABLE_ANGLESPECIAL_COMPARE(ATOM1,ATOM2,ATOM3,RULE,TCHK)
!       ****************************************************************
        CHARACTER(*)        ,INTENT(IN) :: ATOM1
        CHARACTER(*)        ,INTENT(IN) :: ATOM2
        CHARACTER(*)        ,INTENT(IN) :: ATOM3
        TYPE(ANGLERULE_TYPE),INTENT(IN) :: RULE
        LOGICAL(4)          ,INTENT(OUT):: TCHK
!       ****************************************************************
        IF(ATOM2.NE.RULE%ATOM2) THEN
          TCHK=.FALSE.
          RETURN
        END IF
        TCHK=.TRUE.
        IF(ATOM1.EQ.RULE%ATOM1) THEN
          IF(ATOM3.EQ.RULE%ATOM3) RETURN
        ELSE IF(ATOM1.EQ.RULE%ATOM3) THEN
          IF(ATOM3.EQ.RULE%ATOM1) RETURN
        END IF
!       == GENERIC RULES 
        IF(RULE%ATOM3.EQ.'X') THEN
          IF(ATOM1.EQ.RULE%ATOM1.OR.ATOM3.EQ.RULE%ATOM1) RETURN
        END IF
        TCHK=.FALSE.
        RETURN
        END SUBROUTINE UFFTABLE_ANGLESPECIAL_COMPARE
      END SUBROUTINE UFFTABLE_ANGLESPECIAL
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE UFFTABLE$TORSIONPARMS(TYPE1,TYPE2,TYPE3,TYPE4 &
     &            ,BO12,BO23,BO34 &
     &            ,ID,PHI0,NJK,VBARRIER,TCHK)
!     ******************************************************************
!     **  THE TORSION ANGLE IS DEFINED SUCH THAT A TRIANGLE HAS       **
!     **  TORSION ANGLES OF ZERO                                      **
!     **                                                              **
!     ******************************************************************
      USE CLASSICAL_MODULE,ONLY: POT_TYPE
      IMPLICIT NONE
      CHARACTER(*)  ,INTENT(IN) :: TYPE1
      CHARACTER(*)  ,INTENT(IN) :: TYPE2
      CHARACTER(*)  ,INTENT(IN) :: TYPE3
      CHARACTER(*)  ,INTENT(IN) :: TYPE4
      REAL(8)       ,INTENT(IN) :: BO12
      REAL(8)       ,INTENT(IN) :: BO23
      REAL(8)       ,INTENT(IN) :: BO34
      CHARACTER(64) ,INTENT(OUT):: ID
      REAL(8)       ,INTENT(OUT):: PHI0
      INTEGER(4)    ,INTENT(OUT):: NJK
      REAL(8)       ,INTENT(OUT):: VBARRIER
      LOGICAL(4)    ,INTENT(OUT):: TCHK
      REAL(8)       ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      LOGICAL(4)                :: TPR
      INTEGER(4)    ,PARAMETER  :: NDATA=20 ! # OF DATA ENTRIES
      CHARACTER(3)              :: CHEL(NDATA)
      REAL(8)                   :: USP3(NDATA),USP2(5)
      INTEGER(4)                :: IPER(NDATA),IGRP(NDATA)
      INTEGER(4)                :: NTORSION
      REAL(8)                   :: SVAR1,SVAR2
      INTEGER(4)                :: IGRP2,IGRP3,IPER2,IPER3
      INTEGER(4)                :: I
      REAL(8)                   :: VSP3,VSP2
      INTEGER(4)                :: IHYB2,IHYB3
      REAL(8)                   :: KCALBYMOL,SVAR11,SVAR12
!     ******************************************************************
      TCHK=.TRUE.
      TPR=.FALSE.
!
!     ==================================================================
!     == EVALUATE STANDARD TORSION POTENTIALS    =======================
!     ==================================================================
      CHEL( 1)='C_' ;  USP3( 1)=2.119 ;  IGRP( 1)=4 ;  IPER( 1)=2
      CHEL( 2)='N_' ;  USP3( 2)=0.450 ;  IGRP( 2)=5 ;  IPER( 2)=2
      CHEL( 3)='O_' ;  USP3( 3)=0.018 ;  IGRP( 3)=6 ;  IPER( 3)=2
      CHEL( 4)='SI' ;  USP3( 4)=1.225 ;  IGRP( 4)=4 ;  IPER( 4)=3
      CHEL( 5)='P_' ;  USP3( 5)=2.400 ;  IGRP( 5)=5 ;  IPER( 5)=3
      CHEL( 6)='S_' ;  USP3( 6)=0.484 ;  IGRP( 6)=6 ;  IPER( 6)=3
      CHEL( 7)='GE' ;  USP3( 7)=0.701 ;  IGRP( 7)=4 ;  IPER( 7)=4
      CHEL( 8)='AS' ;  USP3( 8)=1.500 ;  IGRP( 8)=5 ;  IPER( 8)=4
      CHEL( 9)='SE' ;  USP3( 9)=0.335 ;  IGRP( 9)=6 ;  IPER( 9)=4
      CHEL(10)='SN' ;  USP3(10)=0.199 ;  IGRP(10)=4 ;  IPER(10)=5
      CHEL(11)='SB' ;  USP3(11)=1.100 ;  IGRP(11)=5 ;  IPER(11)=5
      CHEL(12)='TE' ;  USP3(12)=0.300 ;  IGRP(12)=6 ;  IPER(12)=5
      CHEL(13)='PB' ;  USP3(13)=0.100 ;  IGRP(13)=4 ;  IPER(13)=6
      CHEL(14)='BI' ;  USP3(14)=1.000 ;  IGRP(14)=5 ;  IPER(14)=6
      CHEL(15)='PO' ;  USP3(15)=0.300 ;  IGRP(15)=6 ;  IPER(15)=6
      CHEL(16)='B_' ;  USP3(16)=0.000 ;  IGRP(16)=3 ;  IPER(16)=2
      CHEL(17)='AL' ;  USP3(17)=0.000 ;  IGRP(17)=3 ;  IPER(17)=3
      CHEL(18)='GA' ;  USP3(18)=0.000 ;  IGRP(18)=3 ;  IPER(18)=4
      CHEL(19)='IN' ;  USP3(19)=0.000 ;  IGRP(19)=3 ;  IPER(19)=5
      CHEL(20)='TL' ;  USP3(20)=0.000 ;  IGRP(20)=3 ;  IPER(20)=6
      USP2=(/2.D0,1.25D0,0.7D0,0.2D0,0.1D0/)
!     == ASSIGN GROUP AND PERIOD AND DETERMINE SP3-SP3 TORSIONAL BARRIERS
      SVAR1=0.D0
      SVAR2=0.D0
      IGRP2=0
      IGRP3=0
      IPER2=0
      IPER3=0
      DO I=1,15
        IF(TYPE2(1:2).EQ.CHEL(I)) THEN
          SVAR1=USP3(I)
          IGRP2=IGRP(I)
          IPER2=IPER(I)
        ENDIF
        IF(TYPE3(1:2).EQ.CHEL(I)) THEN
          SVAR2=USP3(I)
          IGRP3=IGRP(I)
          IPER3=IPER(I)
        END IF
      ENDDO
      VSP3=SQRT(SVAR1*SVAR2)
!
!     == DETERMINE SP2-SP2 TORSIONAL BARRIERS
      SVAR1=0.D0
      SVAR2=0.D0
      IF(IPER2-1.GE.1) SVAR1=USP2(IPER2-1)
      IF(IPER3-1.GE.1) SVAR2=USP2(IPER3-1)
      VSP2=5.D0*SQRT(SVAR1*SVAR2)*(1.D0+4.18D0*LOG(BO23)) !RAPPE EQ. 17
!
!     ==================================================================
!     == DETERMINE HYBRIDIZATION                                      ==
!     ==================================================================
      IF(TYPE2(3:3).EQ.'1') THEN
        IHYB2=1
      ELSE IF(TYPE2(3:3).EQ.'2'.OR.TYPE2(3:3).EQ.'R') THEN
        IHYB2=2
      ELSE IF(TYPE2(3:3).EQ.'3') THEN
        IHYB2=3
      ELSE IF(TYPE2(3:3).EQ.'4') THEN
        IHYB2=4
      ELSE IF(TYPE2(3:3).EQ.'6') THEN
        IHYB2=6
      ELSE
        CALL ERROR$MSG('HYBRIDIZATION NOT RECOGNIZED')
        CALL ERROR$CHVAL('TYPE2',TYPE2)
        CALL ERROR$STOP('CLASSICAL_TORSIONPOT')
      END IF
!
      IF(TYPE3(3:3).EQ.'1') THEN
        IHYB3=1
      ELSE IF(TYPE3(3:3).EQ.'2'.OR.TYPE3(3:3).EQ.'R') THEN
        IHYB3=2
      ELSE IF(TYPE3(3:3).EQ.'3') THEN
        IHYB3=3
      ELSE IF(TYPE3(3:3).EQ.'4') THEN
        IHYB3=4
      ELSE IF(TYPE3(3:3).EQ.'6') THEN
        IHYB3=6
      ELSE
        CALL ERROR$MSG('HYBRIDIZATION NOT RECOGNIZED')
        CALL ERROR$CHVAL('TYPE3',TYPE3)
        CALL ERROR$STOP('CLASSICAL_TORSIONPOT')
      END IF
!     ==================================================================
!     ==  SPECIAL CASES                                               ==
!     ==================================================================
      IF((TYPE2.EQ.'ZR3+4'.AND.TYPE3.EQ.'C_3').OR. &
     &   (TYPE3.EQ.'ZR3+4'.AND.TYPE2.EQ.'C_3')) THEN
        VBARRIER=2.D0
        PHI0=0.D0
        NJK=3
        TCHK=.TRUE.
        GOTO 1000
      ENDIF
      IF((TYPE2.EQ.'CIR'.AND.TYPE3.EQ.'C_2').OR. &
     &   (TYPE3.EQ.'CIR'.AND.TYPE2.EQ.'C_2')) THEN
        VBARRIER=1.D0
        PHI0=0.D0
        NJK=4
        TCHK=.TRUE.
        GOTO 1000
      ENDIF
!
!
!     ==================================================================
!     ==  NO TORSIONS WITH DUMMY ATOMS                                ==
!     ==================================================================
      SELECTCASE(TYPE1)
         CASE('CPR','CPR_B','CIR','PIR')
           GOTO 9999
      ENDSELECT
      SELECTCASE(TYPE2)
         CASE('CPR','CPR_B','CIR','PIR')
           GOTO 9999
      ENDSELECT
      SELECTCASE(TYPE3)
         CASE('CPR','CPR_B','CIR','PIR')
           GOTO 9999
      ENDSELECT
      SELECTCASE(TYPE4)
         CASE('CPR','CPR_B','CIR','PIR')
           GOTO 9999
      ENDSELECT
!
!
!     ==================================================================
!     == EXCLUDE CENTRAL BONDS CONTAINING NON-MAINGROUP ELEMENTS =======
!     ==================================================================
      TCHK=.TRUE.
      IF(IGRP2.EQ.0.OR.IGRP3.EQ.0) THEN
        GOTO 9999
!
!     == BOND CONTAINING TWO SP3 CENTERS ===============================
      ELSE IF(IHYB2.EQ.3.AND.IHYB3.EQ.3) THEN     
        NJK=3
        PHI0=PI
        VBARRIER=VSP3
!       == EXCEPTION: BOTH CENTRAL ATOMS ARE GROUP 6 (O IS SPECIAL AGAIN)
        IF(IGRP2.EQ.6.AND.IGRP3.EQ.6) THEN
          IF(TYPE2.EQ.'O_3') THEN
            SVAR11=2.D0
          ELSE
            SVAR11=6.8D0
          END IF
          IF(TYPE3.EQ.'O_3') THEN
            SVAR12=2.D0
          ELSE
            SVAR12=6.8D0
          END IF
          VBARRIER=SQRT(SVAR11*SVAR12)
          NJK=2
          PHI0=0.5D0*PI
        END IF
!
!     == CENTRAL BOND WITH ONE SP2 AND ONE SP3 CENTER    ==============
      ELSE IF((IHYB2.EQ.2.AND.IHYB3.EQ.3).OR. &
     &        (IHYB2.EQ.3.AND.IHYB3.EQ.2)) THEN
        VBARRIER=1.D0
        NJK=6
        PHI0=0.D0
!       == EXCEPTIONS =================================================
        IF(ABS(BO23-1.D0).LE.1.D-1) THEN
!       == CENTRAL BOND IS SINGLE
            IF(((TYPE1(3:3).EQ.'2'.OR.TYPE1(3:3).EQ.'R').AND.IHYB2.EQ.2).OR. &
     &         ((TYPE4(3:3).EQ.'2'.OR.TYPE4(3:3).EQ.'R').AND.IHYB3.EQ.2)) THEN
!           == SP2 OF CENTRAL BOND IS BONDED TO ANOTHER SP2
                VBARRIER=2.D0
                PHI0=PI
                NJK=3
            ENDIF
            IF((IHYB2.EQ.3.AND.IGRP2.EQ.6).OR.(IHYB3.EQ.3.AND.IGRP3.EQ.6)) THEN
!           == CENTRAL BOND WITH SP2 AND GROUP 6 SP3
                VBARRIER=VSP2
                PHI0=0.5D0*PI
                NJK=2
            ENDIF
        ENDIF
!
!     == BOND OF TWO SP2 HYBRIDIZED ATOMS ==============================
      ELSE IF(IHYB2.EQ.2.AND.IHYB3.EQ.2) THEN
        VBARRIER=VSP2
        NJK=2
        PHI0=PI
!     == BOND INVOLVING SP CENTRE ======================================        
      ELSE IF(IHYB2.EQ.1.OR.IHYB3.EQ.1) THEN
        GOTO 9999
      ELSE
        GOTO 9999
      ENDIF
!
!     ==================================================================
!     ==  READ STRING IDENTIFYING POTENTIAL                           ==
!     ==================================================================
      IF(LGT(TYPE2,TYPE3)) THEN
        WRITE(ID,FMT='(A1,2A8," BO=",F5.1," PHI=",F4.0," N=",I2," V=",F3.0)') &
     &                'T ',TYPE2,TYPE3,BO23,PHI0/PI*180.D0,NJK,VBARRIER
      ELSE
        WRITE(ID,FMT='(A1,2A8," BO=",F5.1," PHI",F4.0," N=",I2," V=",F3.0)') &
     &                'T ',TYPE3,TYPE2,BO23,PHI0/PI*180.D0,NJK,VBARRIER
      END IF
!
!     ==================================================================
!     ==  DIVIDE BY THE NUMBER OF TORSIONS ABOUT THIS BOND            ==
!     ==  AND CONVERT INTO ATOMIC UNITS                               ==
!     ==================================================================
 1000 CONTINUE     
      NTORSION=IHYB2*IHYB3
      VBARRIER=VBARRIER/DBLE(NTORSION)  
      CALL CONSTANTS('KCAL/MOL',KCALBYMOL)
      VBARRIER=VBARRIER*KCALBYMOL
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      TCHK=(VBARRIER.GT.1.D-6)
      IF(TCHK) RETURN
 9999 CONTINUE
      ID='ZERO POTENTIAL'
      NTORSION=1
      NJK=1
      PHI0=0.D0
      VBARRIER=0.D0
      TCHK=.FALSE.
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE UFFTABLE$INVERSIONPARMS(TYPE1,TYPE2,TYPE3,TYPE4 &
     &              ,ID,GAMMA0,K,TCHK)
!     ******************************************************************
!     ******************************************************************
      USE CLASSICAL_MODULE,ONLY : POT_TYPE
      IMPLICIT NONE
      CHARACTER(*)  ,INTENT(IN)  :: TYPE1
      CHARACTER(*)  ,INTENT(IN)  :: TYPE2
      CHARACTER(*)  ,INTENT(IN)  :: TYPE3
      CHARACTER(*)  ,INTENT(IN)  :: TYPE4
      CHARACTER(64) ,INTENT(OUT) :: ID
      REAL(8)       ,INTENT(OUT) :: GAMMA0
      REAL(8)       ,INTENT(OUT) :: K
      LOGICAL(4)    ,INTENT(OUT) :: TCHK
      REAL(8)       ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      REAL(8)                    :: KCALBYMOL
!     ******************************************************************
      TCHK=.TRUE.
!
!     ==================================================================
!     ==  NO INVERSIONS WITH DUMMY ATOMS                              ==
!     ==================================================================
      SELECTCASE(TYPE1)
         CASE('CPR','CPR_B','CIR','PIR')
           GOTO 9999
      ENDSELECT
      SELECTCASE(TYPE2)
         CASE('CPR','CPR_B','CIR','PIR')
           GOTO 9999
      ENDSELECT
      SELECTCASE(TYPE3)
         CASE('CPR','CPR_B','CIR','PIR')
           GOTO 9999
      ENDSELECT
      SELECTCASE(TYPE4)
         CASE('CPR','CPR_B','CIR','PIR')
           GOTO 9999
      ENDSELECT
!
!     ==================================================================
!     ==  EVALUATE PARAMETERS                                         ==
!     ==================================================================
      SELECTCASE(TYPE1)
         CASE('C_2','C_R')
            GAMMA0=PI/2.D0   ! PLANAR
            K=6.D0
            IF(TYPE2(1:3).EQ.'O_2'.OR.TYPE3(1:3).EQ.'O_2' &
           &                     .OR.TYPE4(1:3).EQ.'O_2') K=50.D0
         CASE('B_2','N_2','N_R','O_2','O_R','S_2','S_R')
            GAMMA0=PI/2.D0   ! PLANAR
            K=6.D0
         CASE('N_3')
            GAMMA0=5.022534D-1
            K=0.D0
         CASE('P_3+3')
            GAMMA0=9.714677D-2
            K=11.D0          ! K=E(BARRIER)/2
         CASE('AS3+3')
            GAMMA0=5.282239D-2
            K=11.D0          ! K=E(BARRIER)/2
         CASE('SB3+3')
            GAMMA0=4.00605D-2
            K=11.D0          ! K=E(BARRIER)/2
         CASE('BI3+3')
            GAMMA0=0.D0
            K=11.D0          ! K=E(BARRIER)/2
         CASE DEFAULT
            GOTO 9999
      ENDSELECT
!
!     ==================================================================
!     ==  DIVIDE BY THE NUMBER OF INVERSIONS AT THIS SITE             ==
!     ==================================================================
      K=K/3.D0
!
!     ==================================================================
!     ==  READ STRING IDENTIFYING POTENTIAL                           ==
!     ==================================================================
      WRITE(ID,FMT='(A1,2A8," GAMMA0=",F5.0," K=",F5.0)') &
     &              'I ',TYPE1,TYPE2,GAMMA0/PI*180,K
!
!     ==================================================================
!     ==  EVALUATE POTENTIAL ON THE GRID                              ==
!     ==================================================================
      CALL CONSTANTS('KCAL/MOL',KCALBYMOL)
      K=K*KCALBYMOL
!
      TCHK=(K.GT.1.D-6)
      IF(TCHK) RETURN
 9999 CONTINUE
      TCHK=.FALSE.
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE UFFTABLE$NONBONDPARMS(ATOM1,ATOM2,ID,RIJ,DIJ,TCHK)
!     ******************************************************************
!     **                                                              **
!     ** UFF NON BONDING POTENTIAL                                    **
!     **   ATOM1,ATOM2 MNEMOTIC NAME FOR THE ATOMTYPE                 **
!     **                                                              **
!     ******************************************************************
      USE CLASSICAL_MODULE,ONLY : POT_TYPE
      IMPLICIT NONE
      CHARACTER(*)  ,INTENT(IN)  :: ATOM1
      CHARACTER(*)  ,INTENT(IN)  :: ATOM2
      CHARACTER(64) ,INTENT(OUT) :: ID
      REAL(8)       ,INTENT(OUT) :: RIJ
      REAL(8)       ,INTENT(OUT) :: DIJ
      LOGICAL(4)    ,INTENT(OUT) :: TCHK
      REAL(8)                    :: XI,XJ,DI,DJ
      REAL(8)                    :: ANGSTROM,KCALBYMOL
!     ******************************************************************
      CALL UFFTABLE$GET(ATOM1,'X',XI)
      CALL UFFTABLE$GET(ATOM2,'X',XJ)
      CALL UFFTABLE$GET(ATOM1,'D',DI)
      CALL UFFTABLE$GET(ATOM2,'D',DJ)
!
      RIJ=SQRT(XI*XJ)
      DIJ=SQRT(DI*DJ)
      WRITE(ID,FMT='(A1,2A8, " X=",F5.2," D=",F6.4)') &
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
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE TIP4P_MODULE
!*******************************************************************************
!**  TIP4P MODEL OF WATER FROM JOERGENSEN ET AL. J.CHEM.PHYS. 79, 926 (1983)  **
!**                                                                           **
!**  THE MOLECULAR STRUCTURE IS FROZEN.                                       **
!**  A LENNARD-JONES POTENTIAL IS ATTACHED TO THE OXYGEN ION                  **
!**  ELECTROSTATIC CHARGES ARE ATTACHED TO THE HYDROGEN ATOM AND A DUMMY ATOM **
!**  THE DUMMY ATOM IS EXPRESSED BY THE POSITIONS OF THE REAL ATOMS           **
!**                                                                           **
!**  THE TIP4P OBJECT HAS MULTIPLE INSTANCES WHICH NEED TO BE SELECTED AND    **
!**  UNSELECTED WITH THE FUNCTION TIP4P$SELECT(ID) OR TIP4P$SELECT('NONE')    **
!**                                                                           **
!*******************************************************************************
LOGICAL(4)        :: TINI=.FALSE.
REAL(8),PARAMETER :: QO=0.D0
REAL(8),PARAMETER :: QH=0.52D0
REAL(8),PARAMETER :: QM=-2.D0*QH
REAL(8)           :: ROH  ! O-H DISTANCE
REAL(8)           :: RHH  ! H1-H2 DISTANCE
REAL(8)           :: RLJ  ! EQUILIBRIUM DISTANCE LENNARD-JONES DISTANCE
REAL(8)           :: DLJ  ! BINDING ENERGY LENNARD-JONES POTENTIAL
REAL(8)           :: WO   ! RMIDDLE=WO*RO+WH*RH1+WH*RH2
REAL(8)           :: WH
TYPE TIP4P_TYPE
  CHARACTER(32)      :: ID
  INTEGER(4)         :: NWATER
  INTEGER(4),POINTER :: IATS(:,:) => NULL()
  TYPE(TIP4P_TYPE),POINTER :: NEXT
END TYPE TIP4P_TYPE
TYPE(TIP4P_TYPE),TARGET,SAVE :: FIRST
TYPE(TIP4P_TYPE),POINTER     :: THIS =>NULL()
LOGICAL(4)                   :: SELECTED=.FALSE.
REAL(8)                      :: NWATER
INTEGER(4),ALLOCATABLE       :: IATS(:,:)
END MODULE TIP4P_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIP4P$GETR8(ID,VAL)
!     **************************************************************************
!     **  GET REAL(8)-DATA FROM TIP4P OBJECT                                  **
!     **************************************************************************
      USE TIP4P_MODULE, ONLY : TINI,RLJ,DLJ,ROH
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(.NOT.TINI)CALL TIP4P_INI()  
      IF(ID.EQ.'RLJ') THEN
         VAL=RLJ            !VAN DER WAALS RADIUS
      ELSE IF(ID.EQ.'DLJ') THEN
         VAL=DLJ            !VAN DER WAALS BINDING ENERGY
      ELSE IF(ID.EQ.'ROH') THEN
         VAL=ROH            ! O-H DISTANCE
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('TIP4P$GETR8')
      END IF
      RETURN
      END      
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIP4P$SELECT(ID_)
!     **************************************************************************
!     **  SELECT A INSTANCE OF THE CLASSICAL OBJECT                           **
!     **  OR CREATE A NEW ONE IF NOT PRESENT                                  **
!     **************************************************************************
      USE TIP4P_MODULE, ONLY : TINI,SELECTED,FIRST,THIS
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
!     **************************************************************************
      IF(ID_.EQ.'NONE') THEN
        IF(SELECTED) THEN
          SELECTED=.FALSE.
        ELSE
          CALL ERROR$MSG('ATTEMPT TO UNSELECT, WHILE NONE IS SELECTED')
          CALL ERROR$STOP('TIP4P$SELECT')
        END IF
      ELSE
        IF(SELECTED) THEN
          CALL ERROR$MSG('ATTEMPT TO SELECT, WHILE ONE IS ALREADY SELECTED')
          CALL ERROR$CHVAL('ID REQUESTED',ID_)
          CALL ERROR$CHVAL('ID SELECTED ',THIS%ID)
          CALL ERROR$STOP('TIP4P$SELECT')
        ELSE
          SELECTED=.TRUE.
        END IF
      END IF
!
!     ==========================================================================
!     ==  SELECT OR ALLOCATE                                                  ==
!     ==========================================================================
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        CALL TIP4P_INI()  !CALCULATES COMMON PARAMETER 
        THIS=>FIRST
      ELSE
        THIS=>FIRST
        DO 
          IF(THIS%ID.EQ.ID_) RETURN   ! INSTANCE FOUND. DONE...
          IF(ASSOCIATED(THIS%NEXT)) THEN
            THIS=>THIS%NEXT
          ELSE
            ALLOCATE(THIS%NEXT)
            THIS=>THIS%NEXT
            EXIT 
          END IF
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  SET NEW INSTANCE TO ZERO                                            ==
!     ==========================================================================
      THIS%ID=ID_
      THIS%NWATER=-1
      NULLIFY(THIS%IATS)
      NULLIFY(THIS%NEXT)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIP4P_INI()
!     **************************************************************************
!     ** DETERMINE PARAMETERS FOR THE TIP4P MODEL OF WATER                    **
!     ** SOURCE: JOERGENSEN ET AL. J. CHEM. PHYS. 79, 926 (1983)              **
!     **************************************************************************
      USE TIP4P_MODULE, ONLY : ROH,RHH,RLJ,DLJ,WO,WH
      IMPLICIT NONE 
      REAL(8),PARAMETER :: PI=4.D0*ATAN(1.D0)
      REAL(8)        :: ANGSTROM
      REAL(8)        :: KCALBYMOL
      REAL(8)        :: PHIHOH
      REAL(8)        :: ROM
      REAL(8)        :: SVAR
      REAL(8)        :: C12,C6
!     **************************************************************************
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
      CALL CONSTANTS$GET('KCAL/MOL',KCALBYMOL)
!
!     ==========================================================================
!     == DISTANCES BETWEEN THE THREE REAL ATOMS                               ==
!     ==========================================================================
      ROH=0.9572D0*ANGSTROM
      ROM=0.15D0*ANGSTROM
      PHIHOH=104.52D0/180.D0*PI
      RHH=2.D0*ROH*SIN(PHIHOH/2.D0)
!
!     ==========================================================================
!     == WEIGHTS THAT EXPRESS THE DUMMY ATOM M IN TERMS OF REAL ATOMS         ==
!     ==========================================================================
      WH=0.5D0*ROM/(ROH*COS(PHIHOH/2.D0))
      WO=1.D0-2.D0*WH
!
!     ==========================================================================
!     == LENNARD-JONES PARAMETERS                                             ==
!     ==========================================================================
      C12=6.D+5*KCALBYMOL*ANGSTROM**12
      C6=610.D0*KCALBYMOL*ANGSTROM**6
      RLJ=(2.D0*C12/C6)**(1.D0/6.D0) !EQUILIBRIUM DISTANCE
      SVAR=1.D0/RLJ**6
      DLJ=-C6**2/(4.D0*C12)        !EQUILIBRIUM ENERGY D>0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIP4P$EXCLUSIONS(NAT,MAXDIV,NEXCLUSIONX,NEXCLUSION,IEXCLUSION)
!     **************************************************************************
!     ** ADD INTRAMOLECULAR PAIRS TO EXCLUSION FILE                           **
!     **                                                                      **
!     ** IF THE ARRAY-SIZE  IS TOO SMALL THE EXLCUSIONS ARE ONLY COUNTED      **
!     ** TO PROVIDE GUIDANCE FOR ARRAY ALLOCATION                             **
!     **                                                                      **
!     ** IEXCLUSION(:,2) ARE THE EXCLUSION TYPES. FOR TIP4P THE EXCLUSION     **
!     ** TYPE IS SET TO  2, INDICATING COMPLETE EXCLUSION OF LONG-RANGED      **
!     ** FORCES                                                               **
!     **************************************************************************
      USE TIP4P_MODULE    , ONLY : SELECTED,THIS
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NAT
      INTEGER(4),INTENT(IN)    :: MAXDIV
      INTEGER(4),INTENT(IN)    :: NEXCLUSIONX
      INTEGER(4),INTENT(INOUT) :: NEXCLUSION
      INTEGER(4),INTENT(INOUT) :: IEXCLUSION(2,NEXCLUSIONX)
      INTEGER(4)               :: IT(3) ! LATTICE TRANSLATION
      INTEGER(4)               :: IA(4) ! ATOM INDICES OF THE MOLECULE
      INTEGER(4)               :: IEXCL
      INTEGER(4)               :: IW,I,J
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO INSTANCE SELECTED')
        CALL ERROR$STOP('TIP4P$DUMMYPOSITIONS')
      END IF
!
!     ==========================================================================
!     == CONSTRUCT MAPPING ARRAY TO SPEED UP OPERATIONS                       ==
!     ==========================================================================
      IF(THIS%NWATER.LE.0) THEN
        CALL TIP4P$COLLECT()
      END IF
!
!     ==========================================================================
!     == COUNT AND ADD EXCLUSIONS                                             ==
!     ==========================================================================
      IT=0
      DO IW=1,THIS%NWATER
        IA(1)=THIS%IATS(1,IW)
        IA(2)=THIS%IATS(2,IW)
        IA(3)=THIS%IATS(3,IW)
        IA(4)=THIS%IATS(4,IW)
        DO I=1,3
          DO J=I+1,4
            CALL CLASSICAL_ENCODEEXCLUSION(NAT,MAXDIV,IA(I),IA(J),IT,IEXCL)
            NEXCLUSION=NEXCLUSION+1
            IF(NEXCLUSION.LE.NEXCLUSIONX) THEN
              IEXCLUSION(1,NEXCLUSION)=IEXCL
              IEXCLUSION(2,NEXCLUSION)=2
            END IF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIP4P$DUMMY_POSITIONS(NAT,R)
!     **************************************************************************
      USE TIP4P_MODULE    , ONLY : SELECTED,THIS,WO,WH
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NAT
      REAL(8)   ,INTENT(INOUT) :: R(3,NAT)
      INTEGER(4)               :: IO,IM,IH1,IH2
      INTEGER(4)               :: IW
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO INSTANCE SELECTED')
        CALL ERROR$STOP('TIP4P$DUMMYPOSITIONS')
      END IF
!
!     ==========================================================================
!     == CONSTRUCT MAPPING ARRAY TO SPEED UP OPERATIONS                       ==
!     ==========================================================================
      IF(THIS%NWATER.LE.0) THEN
        CALL TIP4P$COLLECT()
      END IF
!
!     ==========================================================================
!     == CONSTRUCT MAPPING ARRAY TO SPEED UP OPERATIONS                       ==
!     ==========================================================================
      DO IW=1,THIS%NWATER
        IO =THIS%IATS(1,IW)
        IM =THIS%IATS(2,IW)
        IH1=THIS%IATS(3,IW)
        IH2=THIS%IATS(4,IW)
        R(:,IM)=WO*R(:,IO)+WH*(R(:,IH1)+R(:,IH2))
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIP4P$DUMMY_FORCES(NAT,F)
!     **************************************************************************
      USE TIP4P_MODULE    , ONLY : SELECTED,THIS,WO,WH
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NAT
      REAL(8)   ,INTENT(INOUT) :: F(3,NAT)
      INTEGER(4)               :: IO,IM,IH1,IH2
      INTEGER(4)               :: IW
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO INSTANCE SELECTED')
        CALL ERROR$STOP('TIP4P$DUMMY_FORCES')
      END IF
!
!     ==========================================================================
!     == CONSTRUCT MAPPING ARRAY TO SPEED UP OPERATIONS                       ==
!     ==========================================================================
      IF(THIS%NWATER.LE.0) THEN
        CALL TIP4P$COLLECT()
      END IF
!
!     ==========================================================================
!     == CONSTRUCT MAPPING ARRAY TO SPEED UP OPERATIONS                       ==
!     ==========================================================================
      DO IW=1,THIS%NWATER
        IO =THIS%IATS(1,IW)
        IM =THIS%IATS(2,IW)
        IH1=THIS%IATS(3,IW)
        IH2=THIS%IATS(4,IW)
        F(:,IO) =F(:,IO) +WO*F(:,IM)
        F(:,IH1)=F(:,IH1)+WH*F(:,IM)
        F(:,IH2)=F(:,IH2)+WH*F(:,IM)
        F(:,IM)=0.D0
      ENDDO
      RETURN
      END
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIP4P$COLLECT()
!     **************************************************************************
      USE CLASSICAL_MODULE, ONLY : MD
      USE TIP4P_MODULE    , ONLY : QM,QH,QO,SELECTED,THIS
      IMPLICIT NONE
      INTEGER(4)              :: IAT,IW,IBOND,IAT1,IAT2
      INTEGER(4)              :: NWATER
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO INSTANCE SELECTED')
        CALL ERROR$STOP('TIP4P$COLLECT')
      END IF
!
!     ==========================================================================
!     == COUNT NUMBER OF TIP4P WATER MOLECULES                                ==
!     ==========================================================================
      NWATER=0
      DO IAT=1,MD%NAT
        IF(MD%TYPE(IAT).EQ.'O_TIP4P') NWATER=NWATER+1
      ENDDO
      THIS%NWATER=NWATER
!
!     ==========================================================================
!     == FILL ATOM POINTERS IATS TO TIP4P WATER MOLECULES                     ==
!     == IATS(:,IW)=(IAT_O,IAT_M,IAT_H2,IAT_H2)                               ==
!     ==========================================================================
      ALLOCATE(THIS%IATS(4,NWATER))
      THIS%IATS(:,:)=0
      IW=0
      DO IAT=1,MD%NAT
        IF(MD%TYPE(IAT).NE.'O_TIP4P') CYCLE
        IW=IW+1
        THIS%IATS(1,IW)=IAT
        DO IBOND=1,MD%NBOND
          IAT1=MD%BOND(IBOND)%IAT1
          IAT2=MD%BOND(IBOND)%IAT2
          IF(IAT1.NE.IAT.AND.IAT2.NE.IAT) CYCLE
          IF(MD%TYPE(IAT1)(3:).NE.'TIP4P') CYCLE
          IF(MD%TYPE(IAT2)(3:).NE.'TIP4P') CYCLE
          IF(IAT2.EQ.IAT)IAT2=IAT1
          IF(MD%TYPE(IAT2).EQ.'M_TIP4P') THEN
            THIS%IATS(2,IW)=IAT2
          ELSE IF(MD%TYPE(IAT2).EQ.'H_TIP4P') THEN
            IF(THIS%IATS(3,IW).EQ.0) THEN
              THIS%IATS(3,IW)=IAT2
            ELSE
              THIS%IATS(4,IW)=IAT2
              EXIT
            END IF
          END IF
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == FILL ATOM POINTERS IATS TO TIP4P WATER MOLECULES                     ==
!     ==========================================================================
      DO IW=1,NWATER
        MD%QEL(THIS%IATS(1,IW))=-QO
        MD%QEL(THIS%IATS(2,IW))=-QM
        MD%QEL(THIS%IATS(3,IW))=-QH
        MD%QEL(THIS%IATS(4,IW))=-QH
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIP4P$CONSTRAINTS(NAT,MASS,R0,RP)
!     **************************************************************************
!     **************************************************************************
      USE TIP4P_MODULE, ONLY : SELECTED,THIS,WO,WH
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NAT
      REAL(8)   ,INTENT(IN)    :: MASS(NAT)
      REAL(8)   ,INTENT(IN)    :: R0(3,NAT)
      REAL(8)   ,INTENT(INOUT) :: RP(3,NAT)
      REAL(8)   ,PARAMETER     :: TOL=1.D-8
      REAL(8)                  :: MO,MH1,MH2   !MASSES
      REAL(8)                  :: RO(3),RH1(3),RH2(3) 
      REAL(8)                  :: ROP(3),RH1P(3),RH2P(3) 
      INTEGER(4)               :: IW
      INTEGER(4)               :: IO,IM,IH1,IH2
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('NO INSTANCE SELECTED')
        CALL ERROR$STOP('TIP4P$CONSTRAINTS')
      END IF
      IF(THIS%NWATER.EQ.0) RETURN
!
!     ==========================================================================
!     == CONSTRUCT MAPPING ARRAY TO SPEED UP OPERATIONS                       ==
!     ==========================================================================
      IF(THIS%NWATER.LE.0) THEN
        CALL TIP4P_INI
        CALL TIP4P$COLLECT()
      END IF
!
!     ==========================================================================
!     == ENFORCE STRUCTURE CONSTRAINTS FOR ALL WATER MOLECULES                ==
!     ==========================================================================
      DO IW=1,THIS%NWATER
        IO =THIS%IATS(1,IW)
        IM =THIS%IATS(2,IW)
        IH1=THIS%IATS(3,IW)
        IH2=THIS%IATS(4,IW)
        MO  =MASS(IO)
        RO  =R0(:,IO)
        ROP =RP(:,IO)
        MH1 =MASS(IH1)
        RH1 =R0(:,IH1)
        RH1P=RP(:,IH1)
        MH2 =MASS(IH2)
        RH2 =R0(:,IH2)
        RH2P=RP(:,IH2)
        CALL TIP4P_CONSTRAINT(MO,MH1,MH2,RO,RH1,RH2,ROP,RH1P,RH2P)
        RP(:,IO)=ROP
        RP(:,IM)=WO*ROP+WH*(RH1P+RH2P)
        RP(:,IH1)=RH1P
        RP(:,IH2)=RH2P
      ENDDO

      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIP4P_CONSTRAINT(MO,MH1,MH2,RO,RH1,RH2,ROP,RH1P,RH2P)
!     **************************************************************************
!     ** GIVEN THE COORDINATES OF THE ACTUAL POSITIONS AND THOSE PROPAGATED   **
!     ** WITHOUT CONSTRAINTS, THE NEXT POSITIONS PROPAGATED WITH CONSTRAINTS  **
!     ** ARE RETURNED                                                         **
!     **************************************************************************
      USE TIP4P_MODULE, ONLY : ROH,RHH
      IMPLICIT NONE 
      INTEGER(4),PARAMETER  :: NITER=200
      REAL(8)   ,PARAMETER  :: TOL=1.D-8
      REAL(8),INTENT(IN)    :: MO      ! OXYGEN MASS
      REAL(8),INTENT(IN)    :: MH1     ! HYDROGEN MASS
      REAL(8),INTENT(IN)    :: MH2     ! HYDROGEN MASS
      REAL(8),INTENT(IN)    :: RO(3)   ! OXYGEN POSITION
      REAL(8),INTENT(IN)    :: RH1(3)  ! 1ST HYDROGEN POSITION
      REAL(8),INTENT(IN)    :: RH2(3)  ! 2ND HYDROGEN POSITION
      REAL(8),INTENT(INOUT) :: ROP(3)  ! PROPAGATED OXYGEN POSITION
      REAL(8),INTENT(INOUT) :: RH1P(3) ! PROPAGATED 1ST HYDROGEN POSITION
      REAL(8),INTENT(INOUT) :: RH2P(3) ! PROPAGATED 2ND HYDROGEN POSITION
      REAL(8)               :: MT      ! MOLECULAR MASS
      REAL(8)               :: G(3)    ! CONSTRAINTS
      REAL(8)               :: LAMBDA(3)  ! LAGRANGE MULTIPLIERS
      REAL(8)               :: DGDR(9,3) ! CONSTRAINT DERIVATIVES
      REAL(8)               :: RP(9)     ! CONTAINER FOR NEXT POSITIONS
      REAL(8)               :: A(3,3)
      REAL(8)               :: AINV(3,3)
      INTEGER(4)            :: ITER
      LOGICAL(4)            :: CONVG
!     **************************************************************************
!     ==========================================================================
!     == CALCULATE THE DERIVATIVE OF THE CONSTRAINT AT THE ORIGINAL POSITIONS ==
!     == G(1)=(RO-RH1)^2-ROH^2
!     == G(2)=(RO-RH2)^2-ROH^2
!     == G(3)=(RH1-RH2)^2-RHH^2
!     == LAMBDA=(ROH1,ROH2,RHH)
!     == R=(RO,RH1,RH2)
!     == DGDR(I,J)=DG_J/DR_I
!     ==========================================================================
      DGDR=0.D0
      DGDR(1:3,1)=2.D0*(RO-RH1)
      DGDR(4:6,1)=2.D0*(RH1-RO)
      DGDR(1:3,2)=2.D0*(RO-RH2)
      DGDR(7:9,2)=2.D0*(RH2-RO)
      DGDR(4:6,3)=2.D0*(RH1-RH2)
      DGDR(7:9,3)=2.D0*(RH2-RH1)
!
!     ==========================================================================
!     == LINEAR MATRIX A_IJ=SUM_K: DG_I/DR_K 1/M_K DGJ/DR_K
!     ==========================================================================
      MT=MO+MH1+MH2  ! TOTAL MASS USED TO VOID LARGE NUMBERS
      DGDR(1:3,:)  =DGDR(1:3,:)/SQRT(MO/MT)
      DGDR(4:6,:)  =DGDR(4:6,:)/SQRT(MH1/MT)
      DGDR(7:9,:)  =DGDR(7:9,:)/SQRT(MH2/MT)
      A=MATMUL(TRANSPOSE(DGDR),DGDR)
      CALL LIB$INVERTR8(3,A,AINV)
!
!     == DGDR=1/M*DGDR =========================================================
      DGDR(1:3,:)  =DGDR(1:3,:)/SQRT(MO/MT)
      DGDR(4:6,:)  =DGDR(4:6,:)/SQRT(MH1/MT)
      DGDR(7:9,:)  =DGDR(7:9,:)/SQRT(MH2/MT)
!
!     ==========================================================================
!     == ITERATE TO FULFILL CONSTRAINTS                                       ==
!     ==========================================================================
      RP(1:3)=ROP
      RP(4:6)=RH1P
      RP(7:9)=RH2P
      DO ITER=1,NITER
        G(1)=SUM((ROP-RH1P)**2)-ROH**2
        G(2)=SUM((ROP-RH2P)**2)-ROH**2
        G(3)=SUM((RH1P-RH2P)**2)-RHH**2
        CONVG=MAXVAL(ABS(G)).LT.TOL
        IF(CONVG) EXIT
        LAMBDA=MATMUL(AINV,G)
        RP=RP-MATMUL(DGDR(:,:),LAMBDA)
        ROP =RP(1:3)
        RH1P=RP(4:6)
        RH2P=RP(7:9)
      ENDDO
      IF(.NOT.CONVG) THEN
        CALL ERROR$MSG('LOOP NOT CONVERGED')
        CALL ERROR$STOP('TIP4P_CONSTRAINT')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_NEIGHBORS(NAT,R,RBAS,NNBX,NNB,NBLIST &
     &                              ,MAXDIV,NEXCL,IEXCLUSION,I2FIRST)
!     **************************************************************************
!     **  CALCULATE A NEIGHBORLIST                                            **
!     **                                                                      **
!     **************************************************************************
      USE CLASSICAL_MODULE, ONLY : NONBOND_TYPE,RCLONGRANGE
      IMPLICIT NONE
      LOGICAL(4)        ,PARAMETER  :: TPR=.FALSE.
      INTEGER(4)        ,INTENT(IN) :: NAT
      REAL(8)           ,INTENT(IN) :: R(3,NAT)
      REAL(8)           ,INTENT(IN) :: RBAS(3,3)
      INTEGER(4)        ,INTENT(IN) :: NNBX
      INTEGER(4)        ,INTENT(OUT):: NNB
      TYPE(NONBOND_TYPE),INTENT(OUT):: NBLIST(NNBX)
      INTEGER(4)        ,INTENT(IN) :: MAXDIV
      INTEGER(4)        ,INTENT(IN) :: NEXCL
      INTEGER(4)        ,INTENT(IN) :: IEXCLUSION(2,NEXCL)
      INTEGER(4)        ,INTENT(IN) :: I2FIRST(NAT)
      INTEGER(4)                    :: EXCLUSION
      REAL(8)                       :: RMAX2
      INTEGER(4)                    :: IAT,IAT1,IAT2,IAT2A,IAT2B,NN,I
      INTEGER(4)                    :: THISTASK,NTASKS
      INTEGER(4)                    :: IEX
      LOGICAL(4)                    :: TCHK
      REAL(8)                       :: D(3),D2
      INTEGER(4)                    :: IT1,IT2,IT3
      INTEGER(4)                    :: I1,I2,I3
      INTEGER(4)                    :: J1,J2,J3
      INTEGER(4)                    :: ITVEC(3)
      INTEGER(4)        ,PARAMETER  :: NDIV(3)=(/3,3,3/) !RECALCULATED LATER
      REAL(8)                       :: X(3)
      INTEGER(4)                    :: IDIV(3,NAT),ITR(3,NAT)
      INTEGER(4)                    :: IATPNT(NAT)
      INTEGER(4)        ,ALLOCATABLE:: NATINBOX(:,:,:)
      INTEGER(4)        ,ALLOCATABLE:: CURRENT(:,:,:)
      INTEGER(4)        ,ALLOCATABLE:: FIRST(:,:,:)
      REAL(8)                       :: RBASINV(3,3)
      REAL(8)                       :: TBOX(3,3)
      REAL(8)                       :: X0,Y0,Z0
      INTEGER(4)                    :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      NBLIST(:)%TONEFOUR=.FALSE.
!
!     ==========================================================================
!     == DIVIDE THE UNIT CELL INTO BOXES AND ATTRIBUTE EACH ATOM TO A BOX     ==
!     ==========================================================================
      CALL LIB$INVERTR8(3,RBAS,RBASINV)
      ALLOCATE(NATINBOX(NDIV(1),NDIV(2),NDIV(3)))
      NATINBOX(:,:,:)=0
      DO IAT=1,NAT
!       == TRANSFORM COORDINATES TO RELATIVE COORDINATES
        X(:)=MATMUL(RBASINV,R(:,IAT))
!       == MAP ATOMS INTO FIRST UNIT CELL AND KEEP THE LATTICE TRANSLATION =====
!       == R-RBAS*IT LIES IN THE FIRST UNIT CELL ===============================
        DO I=1,3
          IF(X(I).GE.0.D0) THEN  ! ROUND DOWN
            ITR(I,IAT)=INT(X(I))
          ELSE
            ITR(I,IAT)=INT(X(I))-1
          END IF
        ENDDO
        X(:)=X(:)-REAL(ITR(:,IAT),KIND=8)
!       == IDENTIFY SUB-BOX ====================================================
!       == IDIV ARE THE NUMBERS OF SUB-BOX TRANSLATIONS POINTING TO THE ORIGIN==
!       == OF THE RESPECTIVE SUB=BOX ===========================================
        IDIV(:,IAT)=INT(X(:)*REAL(NDIV(:),KIND=8))
!       == AVOID ROUNDING ERROR IF X IS ROUNDED TO 1.D0 ========================
        IDIV(:,IAT)=MIN(NDIV(:)-1,IDIV(:,IAT))
        I1=1+IDIV(1,IAT)
        I2=1+IDIV(2,IAT)
        I3=1+IDIV(3,IAT)
!       == INCREASE THE ATOM COUNTER OF THE RESPECTIVE SUB-BOX =================
        NATINBOX(I1,I2,I3)=NATINBOX(I1,I2,I3)+1
      ENDDO
!
!     ==========================================================================
!     ==  FIRST(I,J,K) POINTS TO THE FIRST ATOM IN THE GIVEN BOX              ==
!     ==========================================================================
      ALLOCATE(FIRST(NDIV(1),NDIV(2),NDIV(3)))
      I=0
      DO I3=1,NDIV(3)
        DO I2=1,NDIV(2)
          DO I1=1,NDIV(1)
            FIRST(I1,I2,I3)=I+1
            I=I+NATINBOX(I1,I2,I3)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  IATPNT POINTS FROM THE ATOMS IN A GIVEN BOX TO THE BARE ATOM INDEX  ==
!     ==========================================================================
      IATPNT(:)=0
      ALLOCATE(CURRENT(NDIV(1),NDIV(2),NDIV(3)))
      CURRENT(:,:,:)=FIRST(:,:,:)  ! A SUB-BOX-SPECIFIC ATOM POINTER
      DO IAT=1,NAT
        I1=IDIV(1,IAT)+1
        I2=IDIV(2,IAT)+1
        I3=IDIV(3,IAT)+1
        IATPNT(CURRENT(I1,I2,I3))=IAT
        CURRENT(I1,I2,I3)=CURRENT(I1,I2,I3)+1
      ENDDO
      DEALLOCATE(CURRENT)
!
!     ==========================================================================
!     == SET UP NEIGHBORLIST                                                  ==
!     ==========================================================================
      RMAX2=RCLONGRANGE**2
      NNB=0
      DO I=1,3
        TBOX(:,I)=RBAS(:,I)/REAL(NDIV(I),KIND=8) ! SUB-BOX TRANSLATION VECTORS
      ENDDO
      DO IAT1=1,NAT
        X0=R(1,IAT1)
        Y0=R(2,IAT1)
        Z0=R(3,IAT1)
        CALL BOXSPH(TBOX,X0,Y0,Z0,RCLONGRANGE,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
!       == LOOP OVER BOXES IN THE NEIGHBORHOOD ================================
        DO I1=MIN1-1,MAX1
          J1=MODULO(I1,NDIV(1))
          IT1=(I1-J1)/NDIV(1)
          J1=J1+1
          DO I2=MIN2-1,MAX2
            J2=MODULO(I2,NDIV(2))
            IT2=(I2-J2)/NDIV(2)
            J2=J2+1
            DO I3=MIN3-1,MAX3
              J3=MODULO(I3,NDIV(3))
              IT3=(I3-J3)/NDIV(3)
              J3=J3+1
              IAT2A=FIRST(J1,J2,J3)
              IAT2B=IAT2A-1+NATINBOX(J1,J2,J3)
              DO I=IAT2A,IAT2B
                IAT2=IATPNT(I)
                ITVEC(1)=IT1-ITR(1,IAT2)
                ITVEC(2)=IT2-ITR(2,IAT2)
                ITVEC(3)=IT3-ITR(3,IAT2)
!               == AVOID DOUBLE COUNTING=======================================
                IF(IAT1.GT.IAT2) CYCLE
                IF(IAT2.EQ.IAT1) THEN
                  IF(ITVEC(3).GT.0) CYCLE
                  IF(ITVEC(3).EQ.0) THEN
                    IF(ITVEC(2).GT.0) CYCLE
                    IF(ITVEC(2).EQ.0) THEN
                      IF(ITVEC(1).GE.0) CYCLE
                    END IF
                  END IF
                END IF   
!               == DISTANCE CRITERION =========================================
                D(:)=R(:,IAT2)+MATMUL(RBAS,REAL(ITVEC,KIND=8))-R(:,IAT1)
                D2=SUM(D(:)**2)
                IF(D2.GT.RMAX2) CYCLE
                NNB=NNB+1
                IF(NNB.LE.NNBX) THEN
                  NBLIST(NNB)%IAT1=IAT1
                  NBLIST(NNB)%IAT2=IAT2
                  TCHK=(ITVEC(1).NE.0).OR.(ITVEC(2).NE.0).OR.(ITVEC(3).NE.0)
                  IF(ASSOCIATED(NBLIST(NNB)%IT)) THEN
                    IF(TCHK) THEN
                      NBLIST(NNB)%IT(:)=ITVEC(:)
                    ELSE
                      DEALLOCATE(NBLIST(NNB)%IT)
                    END IF
                  ELSE
                    IF(TCHK) THEN
                      ALLOCATE(NBLIST(NNB)%IT(3))
                      NBLIST(NNB)%IT(:)=ITVEC(:)
                    END IF
                  END IF                               
!         == IF ISVAR=IAT1-1+NAT*(IAT2-1)+...
!         == THEN 1) IT IS THE TRANSLATION OF IAT1
!         ==      2) IAT1.GE.IAT2
                  NBLIST(NNB)%EXCLUDE=.FALSE.
                  IF(ABS(ITVEC(1)).GT.MAXDIV) CYCLE                  
                  IF(ABS(ITVEC(2)).GT.MAXDIV) CYCLE                  
                  IF(ABS(ITVEC(3)).GT.MAXDIV) CYCLE                  
                  EXCLUSION=1+(ITVEC(1)+MAXDIV) &
               &             +(1+2*MAXDIV)*((ITVEC(2)+MAXDIV) &
               &             +(1+2*MAXDIV)*((ITVEC(3)+MAXDIV) &
               &             +(1+2*MAXDIV)*(IAT2-1+NAT*(IAT1-1))))
                  IAT2A=I2FIRST(IAT1)
                  NBLIST(NNB)%EXCLUDE=.FALSE.
                  NBLIST(NNB)%TONEFOUR=.FALSE.
                  DO IEX=IAT2A,NEXCL
                     IF(IEXCLUSION(1,IEX).LT.EXCLUSION) CYCLE
                     IF(IEXCLUSION(1,IEX).GT.EXCLUSION) EXIT
!                    == IEXCLUSION(1,IEX).EQ.EXCLUSION: EXCLUSION ENCOUNTERED ==
                     IF(IEXCLUSION(2,IEX).EQ.0) THEN
!                      __ 1,2 OR 1,3 EXCLUSION__________________________________
                       NBLIST(NNB)%EXCLUDE=.TRUE. 
                       NBLIST(NNB)%TONEFOUR=.FALSE.
                     ELSE IF(IEXCLUSION(2,IEX).EQ.1) THEN
!                      __1,4 EXCLUSION__________________________________________
                       NBLIST(NNB)%EXCLUDE=.TRUE. 
                       NBLIST(NNB)%TONEFOUR=.TRUE.
                     ELSE IF(IEXCLUSION(2,IEX).EQ.2) THEN
!                      __TOTAL EXCLUSION. REMOVE FROM NEIGHBORLIST______________
                       NNB=NNB-1
                     ELSE
                       CALL ERROR$MSG('EXCLUSION TYPE NOT RECOGNIZED')
                       CALL ERROR$I4VAL('EXCLUSION TYPE',IEXCLUSION(2,IEX))
                       CALL ERROR$I4VAL('IAT1',IAT1)
                       CALL ERROR$I4VAL('IAT2',IAT2)
                       CALL ERROR$I4VAL('IT(1)',ITVEC(1))
                       CALL ERROR$I4VAL('IT(2)',ITVEC(2))
                       CALL ERROR$I4VAL('IT(3)',ITVEC(3))
                       CALL ERROR$STOP('CLASSICAL_NEIGHBORS')
                     END IF
                     EXIT
                  END DO
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  IF ARRAY TOO SMALL RETURN WITH AN UNUSABLE NEIGHBORLIST     ==
!     ==================================================================
      IF(NNB.GT.NNBX) THEN
        NBLIST(:)%IAT1=0
        NBLIST(:)%IAT2=0
        RETURN
      END IF
!
!     ==================================================================
!     ==================================================================
!     ==================================================================
      IF(TPR) THEN
        PRINT*,' NUMBER OF NEIGHBORS ',NNB
        DO NN=1,MIN(NNB,NNBX)
          IF(ASSOCIATED(NBLIST(NN)%IT)) THEN         
            WRITE(*,FMT='("NN: IAT1=",I5," IAT2=",I5," T=",3I5)') &
     &           NBLIST(NN)%IAT1,NBLIST(NN)%IAT2,NBLIST(NN)%IT
          ELSE
            WRITE(*,FMT='("NN: IAT1=",I5," IAT2=",I5,"EXCL:",L4)') &
                 NBLIST(NN)%IAT1,NBLIST(NN)%IAT2,NBLIST(NN)%EXCLUDE
          END IF
        ENDDO        
      END IF
      RETURN
      END SUBROUTINE CLASSICAL_NEIGHBORS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_ENCODEEXCLUSION(NAT,MAXDIV,IAT1,IAT2,IT,IEXCL)
!     **************************************************************************
!     **  EXCLUSIONS DEFINE A PAIR OF ATOMS THAT FOR WHICH NO NON-BOND FORCES **
!     **  SHALL BE CALCULATED. AN EXCLUSION IS SPECIFIED BY TWO ATOM INDICES  **
!     **  AND A LATTICE TRANSLATION OF THE SECOND ATOM. THESE FIVE INTEGERS   **
!     **  ARE MAPPEN ONTO A SINGLE NUMBER.                                    **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: MAXDIV
      INTEGER(4),INTENT(IN) :: IAT1
      INTEGER(4),INTENT(IN) :: IAT2
      INTEGER(4),INTENT(IN) :: IT(3)
      INTEGER(4),INTENT(OUT):: IEXCL
!     **************************************************************************
      IEXCL=1+(IT(1)+MAXDIV)+(1+2*MAXDIV)*((IT(2)+MAXDIV) &
     &                      +(1+2*MAXDIV)*((IT(3)+MAXDIV) &
     &                      +(1+2*MAXDIV)*(IAT2-1+NAT*(IAT1-1))))
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_DECODEEXCLUSION(NAT,MAXDIV,IAT1,IAT2,IT,IEXCL)
!     **************************************************************************
!     **  EXCLUSIONS DEFINE A PAIR OF ATOMS THAT FOR WHICH NO NON-BOND FORCES **
!     **  SHALL BE CALCULATED. AN EXCLUSION IS SPECIFIED BY TWO ATOM INDICES  **
!     **  AND A LATTICE TRANSLATION OF THE SECOND ATOM. THESE FIVE INTEGERS   **
!     **  ARE MAPPEN ONTO A SINGLE NUMBER.                                    **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NAT
      INTEGER(4),INTENT(IN)  :: MAXDIV
      INTEGER(4),INTENT(OUT) :: IAT1
      INTEGER(4),INTENT(OUT) :: IAT2
      INTEGER(4),INTENT(OUT) :: IT(3)
      INTEGER(4),INTENT(IN)  :: IEXCL
      INTEGER(4)             :: ID
!     **************************************************************************
!     == IEXCL=1+(IT(1)+MAXDIV)+(1+2*MAXDIV)*((IT(2)+MAXDIV) ===================
!     ==                       +(1+2*MAXDIV)*((IT(3)+MAXDIV) &
!     ==                      +(1+2*MAXDIV)*(IAT2-1+NAT*(IAT1-1))))
      ID=IEXCL-1
      IT(1)=MOD(ID,(1+2*MAXDIV))-MAXDIV
      ID=INT(ID/(1+2*MAXDIV))
      IT(2)=MOD(ID,(1+2*MAXDIV))-MAXDIV
      ID=INT(ID/(1+2*MAXDIV))
      IT(3)=MOD(ID,(1+2*MAXDIV))-MAXDIV
      ID=INT(ID/(1+2*MAXDIV))
      IAT2=MOD(ID,NAT)+1
      ID=INT(ID/NAT)
      IAT1=ID+1
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_PRINTEXCLUSION(NAT,MAXDIV,NEXCLX,NEXCL,IEXCL)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NAT
      INTEGER(4),INTENT(IN)  :: MAXDIV
      INTEGER(4),INTENT(IN)  :: NEXCLX
      INTEGER(4),INTENT(IN)  :: NEXCL
      INTEGER(4),INTENT(IN)  :: IEXCL(2,NEXCLX)
      INTEGER(4)             :: I
      INTEGER(4)             :: IAT1,IAT2,IT(3),IEXCL1
!     **************************************************************************
      DO I=1,NEXCL
        IF(IEXCL(1,I).EQ.-1) EXIT
        CALL CLASSICAL_DECODEEXCLUSION(NAT,MAXDIV,IAT1,IAT2,IT,IEXCL(1,I))
        CALL CLASSICAL_ENCODEEXCLUSION(NAT,MAXDIV,IAT1,IAT2,IT,IEXCL1)
        WRITE(*,FMT='("IAT1=",I5," IAT2=",I5," IT=",3I4," TEST=",2I10)') &
     &                     IAT1,IAT2,IT,IEXCL(1,I),IEXCL1
      ENDDO        
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_EXCLUSIONS(NAT,NBOND,BOND,NANGLE,ANGLE,NTORSION &
     &               ,TORSION,MAXDIV,NEXCLUSIONX,NEXCLUSION,IEXCLUSION,I2FIRST)
!     **************************************************************************
!     **  NON-BOND INTERACTIONS SUCH AS VAN DER WAALS AND COULOMB INTERACTION **
!     **  MUST NOT BE CALCULATED THAT ARE CONNECTED VIA TWO BONDS OR LESS.    **
!     **  THE EXCLUSIONS ARE A LIST OF THOSE NEIGHBORS.                       **
!     **                                                                      **
!     **  ONE EXCLUSION CONTAINS THE INFORMATION OF THE TWO ATOMS IN THE PAIR **
!     **  AND THE INTEGER TRANSLATIONS OF THE SECOND ATOM. IN ORDER TO        **
!     **  REDUCE MEMORY, ALL FIVE NUMBERS ARE MAPPED ONTO A SINGLE INTEGER    **
!     **  NUMBER.                                                             **
!     **                                                                      **
!     **  NONBOND-EXCLUSIONS STORED AS 1-D ARRAY                              **
!     **  1+(IT1+2)+5*((IT2+2)+5*(IT3+2)+5*((IAT1-1)+NAT*(IAT2-1)             **
!     **  SO THAT IAT1>IAT2                                                   **
!     **************************************************************************
      USE CLASSICAL_MODULE, ONLY : BOND_TYPE,ANGLE_TYPE,TORSION_TYPE,MD
      IMPLICIT NONE
      INTEGER(4)       ,INTENT(IN) :: NAT
      INTEGER(4)       ,INTENT(IN) :: NBOND
      TYPE(BOND_TYPE)  ,INTENT(IN) :: BOND(NBOND)
      INTEGER(4)       ,INTENT(IN) :: NANGLE
      TYPE(ANGLE_TYPE) ,INTENT(IN) :: ANGLE(NANGLE)
      INTEGER(4)       ,INTENT(IN) :: NTORSION
      TYPE(TORSION_TYPE) ,INTENT(IN) :: TORSION(NTORSION)
      INTEGER(4)       ,INTENT(IN) :: MAXDIV
      INTEGER(4)       ,INTENT(IN) :: NEXCLUSIONX
      INTEGER(4)       ,INTENT(OUT):: NEXCLUSION
      INTEGER(4)       ,INTENT(OUT):: IEXCLUSION(2,NEXCLUSIONX)
      INTEGER(4)       ,INTENT(OUT):: I2FIRST(NAT)
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
      INTEGER(4)            :: LARGEST    ! LARGEST NUMBER ON EXLCUSION FILE
      INTEGER(4)            :: IB,IANGLE,IAT1,IAT2,I,I1,IMAX,ISVAR,ITORSION
      INTEGER(4)            :: IT(3),IT1(3),IT3(3),IT4(3)
      REAL(8)               :: WORK(NEXCLUSIONX)
      INTEGER(4)            :: FROM, TO, ISVAR1, ISVAR2
!     **************************************************************************
      LARGEST=1+(1+2*MAXDIV)**3+NAT**2+1
      IF(LARGEST.GT.HUGE(IEXCLUSION)) THEN
        CALL ERROR$MSG('NUMBERS LARGER THAN NUMBER MODEL')
        CALL ERROR$STOP('CLASSICAL_EXCLUSIONS')
      END IF
      NEXCLUSION=0
      IEXCLUSION(1,:)=-1
      IEXCLUSION(2,:)=-1  
                          ! 0 STANDS FOR A 1-2 OR 1-3 EXCLUSION
                          ! 1 STANDS FOR 1-4 EXCLUSION
                          ! 2 STANDS FOR DEFINITELY NO INTERACTION
!
!     ==========================================================================
!     ==  TIP4P-INTRAMOLECULAR EXCLUSIONS                                     ==
!     ==========================================================================
      CALL TIP4P$EXCLUSIONS(NAT,MAXDIV,NEXCLUSIONX,NEXCLUSION,IEXCLUSION)
!
!     ==========================================================================
!     ==  FIND BOND EXCLUSIONS                                               ==
!     ==========================================================================
      DO IB=1,NBOND
        IAT1=BOND(IB)%IAT1
        IAT2=BOND(IB)%IAT2
        NEXCLUSION=NEXCLUSION+1
        IF(NEXCLUSION.LE.NEXCLUSIONX) THEN
          IT(:)=0
          IF(ASSOCIATED(BOND(IB)%IT2)) THEN
            IT(:)=BOND(IB)%IT2(:)
          END IF
          IF(IAT2.GE.IAT1) THEN 
            ISVAR=IAT2-1+NAT*(IAT1-1)
            IF(IAT2.EQ.IAT1) THEN
              IF(IT(3).GT.0) THEN
                IT(:)=-IT(:)
              ELSE IF(IT(3).EQ.0)  THEN
                IF(IT(2).GT.0) THEN
                  IT(:)=-IT(:)
                ELSE IF(IT(2).EQ.0)  THEN
                  IF(IT(1).GT.0) THEN
                    IT(:)=-IT(:)
                  END IF
                END IF
              END IF
            END IF
          ELSE
            ISVAR=IAT1-1+NAT*(IAT2-1)
            IT(:)=-IT(:)
          END IF
          IF(ABS(IT(1)).GT.MAXDIV.OR.ABS(IT(2)).GT.MAXDIV &
     &                           .OR.ABS(IT(3)).GT.MAXDIV) THEN
            CALL ERROR$MSG('BOND TRANSLATIONS TOO LARGE')
            CALL ERROR$MSG('CANNOT HANDLE AT PRESENT')
            CALL ERROR$I4VAL('MAXDIV',MAXDIV)
            CALL ERROR$I4VAL('IT(1)',IT(1))
            CALL ERROR$I4VAL('IT(2)',IT(2))
            CALL ERROR$I4VAL('IT(3)',IT(3))
            CALL ERROR$STOP('CLASSICAL_EXCLUSIONS')
          END IF
!         == IT IS THE TRANSLATION OF THE  ATOM IAT1 ===========================
!         == IF ISVAR=IAT1-1+NAT*(IAT2-1) ======================================
          ISVAR=(IT(1)+MAXDIV)+(1+2*MAXDIV)*((IT(2)+MAXDIV) &
     &                        +(1+2*MAXDIV)*((IT(3)+MAXDIV) &
     &                        +(1+2*MAXDIV)*ISVAR))
          IEXCLUSION(1,NEXCLUSION)=1+ISVAR
          IEXCLUSION(2,NEXCLUSION)=0   ! 0 STANDS FOR A 1-2 OR 1-3 EXCLUSION
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  FIND ANGLE EXLCUSIONS                                               ==
!     ==========================================================================
      DO IANGLE=1,NANGLE
        IAT1=ANGLE(IANGLE)%IAT1
        IAT2=ANGLE(IANGLE)%IAT3
        NEXCLUSION =NEXCLUSION+1
        IF(NEXCLUSION.LE.NEXCLUSIONX) THEN
          IT1(:)=0
          IT3(:)=0
          IF(ASSOCIATED(ANGLE(IANGLE)%IT1)) IT1(:)=ANGLE(IANGLE)%IT1(:)
          IF(ASSOCIATED(ANGLE(IANGLE)%IT3)) IT3(:)=ANGLE(IANGLE)%IT3(:)
          IT(:)=IT3(:)-IT1(:)
          IF(IAT2.GE.IAT1) THEN
            ISVAR=IAT2-1+NAT*(IAT1-1)
            IF(IAT2.EQ.IAT1) THEN
              IF(IT(3).GT.0) THEN
                IT(:)=-IT(:)
              ELSE IF(IT(3).EQ.0)  THEN
                IF(IT(2).GT.0) THEN
                  IT(:)=-IT(:)
                ELSE IF(IT(2).EQ.0)  THEN
                  IF(IT(1).GT.0) THEN
                    IT(:)=-IT(:)
                  END IF
                END IF
              END IF
            END IF
          ELSE
            ISVAR=IAT1-1+NAT*(IAT2-1)
            IT(:)=-IT(:)
          END IF
          IF(ABS(IT(1)).GT.MAXDIV.OR.ABS(IT(2)).GT.MAXDIV &
     &                           .OR.ABS(IT(3)).GT.MAXDIV) THEN
            CALL ERROR$MSG('BOND TRANSLATIONS TOO LARGE')
            CALL ERROR$MSG('CANNOT HANDLE AT PRESENT')
            CALL ERROR$I4VAL('MAXDIV',MAXDIV)
            CALL ERROR$I4VAL('IT(1)',IT(1))
            CALL ERROR$I4VAL('IT(2)',IT(2))
            CALL ERROR$I4VAL('IT(3)',IT(3))
            CALL ERROR$STOP('CLASSICAL_EXCLUSIONS')
          END IF
!         == IF ISVAR=IAT1-1+NAT*(IAT2-1)+...
!         == THEN 1) IT IS THE TRANSLATION OF IAT1
!         ==      2) IAT1.GE.IAT2
          ISVAR=(IT(1)+MAXDIV)+(1+2*MAXDIV)*((IT(2)+MAXDIV) &
     &                        +(1+2*MAXDIV)*((IT(3)+MAXDIV) &
     &                        +(1+2*MAXDIV)*ISVAR))
          IEXCLUSION(1,NEXCLUSION)=1+ISVAR
          IEXCLUSION(2,NEXCLUSION)=0   ! 0 STANDS FOR A 1-2 OR 1-3 EXCLUSION
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  FIND TORSION EXCLUSIONS                                             ==
!     ==========================================================================
PRINT*,"FLAG: NTORSION=",NTORSION
      IF(MD%FF.NE.'UFF') THEN
         DO ITORSION=1,NTORSION
            IAT1=TORSION(ITORSION)%IAT1
            IAT2=TORSION(ITORSION)%IAT4
            NEXCLUSION =NEXCLUSION+1
            IF(NEXCLUSION.LE.NEXCLUSIONX) THEN
               IT1(:)=0
               IT4(:)=0
               IF(ASSOCIATED(TORSION(ITORSION)%IT1)) &
    &                        IT1(:)=TORSION(ITORSION)%IT1(:)
               IF(ASSOCIATED(TORSION(ITORSION)%IT4)) &
    &                        IT4(:)=TORSION(ITORSION)%IT4(:)
               IT(:)=IT4(:)-IT1(:)
               IF(IAT2.GE.IAT1) THEN
                  ISVAR=IAT2-1+NAT*(IAT1-1)
                  IF(IAT2.EQ.IAT1) THEN
                     IF(IT(3).GT.0) THEN
                        IT(:)=-IT(:)
                     ELSE IF(IT(3).EQ.0)  THEN
                        IF(IT(2).GT.0) THEN
                           IT(:)=-IT(:)
                        ELSE IF(IT(2).EQ.0)  THEN
                           IF(IT(1).GT.0) THEN
                              IT(:)=-IT(:)
                           END IF
                        END IF
                     END IF
                  END IF
               ELSE
                  ISVAR=IAT1-1+NAT*(IAT2-1)
                  IT(:)=-IT(:)
               END IF
               IF(ABS(IT(1)).GT.MAXDIV.OR.ABS(IT(2)).GT.MAXDIV &
       &                              .OR.ABS(IT(3)).GT.MAXDIV) THEN
                  CALL ERROR$MSG('BOND TRANSLATIONS TOO LARGE')
                  CALL ERROR$MSG('CANNOT HANDLE AT PRESENT')
                  CALL ERROR$I4VAL('MAXDIV',MAXDIV)
                  CALL ERROR$I4VAL('IT(1)',IT(1))
                  CALL ERROR$I4VAL('IT(2)',IT(2))
                  CALL ERROR$I4VAL('IT(3)',IT(3))
                  CALL ERROR$STOP('CLASSICAL_EXCLUSIONS')
               END IF
!         == IF ISVAR=IAT1-1+NAT*(IAT2-1)+...
!         == THEN 1) IT IS THE TRANSLATION OF IAT1
!         ==      2) IAT1.GE.IAT2
               ISVAR=(IT(1)+MAXDIV)+(1+2*MAXDIV)*((IT(2)+MAXDIV) &
                    &                        +(1+2*MAXDIV)*((IT(3)+MAXDIV) &
                    &                        +(1+2*MAXDIV)*ISVAR))
               IEXCLUSION(1,NEXCLUSION)=1+ISVAR
               IEXCLUSION(2,NEXCLUSION)=1
! DO NOT SET IEXCLUSION (2,NEXCL...) TO AVOID OVERWRITING 1-2 AND 1-3 EXCLUSIONS
! IT IS SET TO 1 BEFORE THE 1-2 AND 1-3 EXCLUSIONS ARE EVALUATED

! PRINT*,"FLAG: ANGLE EXCLUSION FOUND: ",IAT1, IAT2
! C=IEXCLUSION(NEXCLUSION)
! A= MAXDIV
! B= (1+2*MAXDIV)
! D= INT( (REAL(C - 1 - A - B*A - A*B**2))/REAL(B**3) +1 - NAT*(IAT1-1))
! E= INT( (REAL(C-1-A-B*A-A*B**2)) / REAL(NAT*B**3) + REAL(1)/REAL(NAT) 
!                                                   - REAL(IAT2)/REAL(NAT) + 1)
! PRINT*,"FLAG: ANGLEEXCL: ", IAT1, IAT2, C, D, E 
! PRINT*,"FLAG: ANGLEEXCL: ",NEXCLUSION,IAT1, IAT2, IT(:)
! PRINT*,"FLAG: ANGLEEXCL: ",IEXCLUSION(NEXCLUSION),NAT,MAXDIV
            END IF
         ENDDO
      END IF
!
!     ==========================================================================
!     ==  RETURN THE NUMBER OF EXCLUSIONS IF GREATER THAN MAXIMUM             ==
!     ==========================================================================
      IF(NEXCLUSION+1.GT.NEXCLUSIONX) THEN
!       == ADD ONE ELEMENT FOR A FLAG INDICATING THE LAST ELEMENT
        NEXCLUSION=NEXCLUSION+1
        RETURN
      END IF
!
!     ==========================================================================
!     ==  SORT IN INCREASING MAGNITUDE                                        ==
!     ==========================================================================
      IF(NEXCLUSION.GT.1) THEN
       WORK(:)=REAL(IEXCLUSION(1,:))
       CALL SORT$SET(NEXCLUSION,WORK)
       CALL SORT$RESTART()
       CALL SORT$FLIP(FROM,TO)
       DO WHILE (FROM.NE.0.OR.TO.NE.0)
         IF(TO.EQ.0) THEN
           ISVAR1=IEXCLUSION(1,FROM)
           ISVAR2=IEXCLUSION(2,FROM)
         ELSE IF(FROM.EQ.0) THEN
           IEXCLUSION(1,TO)=ISVAR1
           IEXCLUSION(2,TO)=ISVAR2
         ELSE
           IEXCLUSION(1,TO)=IEXCLUSION(1,FROM)
           IEXCLUSION(2,TO)=IEXCLUSION(2,FROM)
         END IF
         CALL SORT$FLIP(FROM,TO)
       ENDDO
       CALL SORT$UNSET
!
!       ==  EXCLUDE DOUBLE EXCLUSIONS ==========================================
        I1=1
        DO I=2,NEXCLUSION
          IF(IEXCLUSION(1,I1).NE.IEXCLUSION(1,I)) THEN
            I1=I1+1
            IEXCLUSION(1,I1)=IEXCLUSION(1,I)        
            IEXCLUSION(2,I1)=IEXCLUSION(2,I)        
          ELSE
!           __REPORT DISCARDED EXCLUSIONS_______________________________________
            IAT2=(IEXCLUSION(1,I)-1)/NAT+1
            IAT1=IEXCLUSION(1,I)-NAT*(IAT2-1)
            PRINT*,'DOUBLE EXCLUSION BETWEEN ATOMS:',IAT1,IAT2
!           == UNRAVEL PREVALENCE OF EXCLUSION TYPES ===========================
            IF(IEXCLUSION(2,I).EQ.2) THEN
!             __TOTAL EXCLUSION HAS PREVALENCE OVER ALL OTHERS__________________
              IEXCLUSION(2,I1)=IEXCLUSION(2,I)
            ELSE IF(IEXCLUSION(2,I).EQ.0) THEN
!             __1,2 AND 1,3 EXCLUSIONS HAVE PREVALENCE OVER 1,4 EXCLUSIONS______
              IF(IEXCLUSION(2,I1).EQ.1) THEN
                IEXCLUSION(2,I1)=IEXCLUSION(2,I)
              END IF
            ELSE IF(IEXCLUSION(2,I).EQ.0) THEN
            ELSE
              CALL ERROR$MSG('EXCLUSION TYPE NOT RECOGNIZED')
              CALL ERROR$I4VAL('IAT1',IAT1)
              CALL ERROR$I4VAL('IAT2',IAT2)
              CALL ERROR$STOP('CLASSICAL_EXCLUSIONS')
            END IF
          END IF
        ENDDO
        IF(NEXCLUSION.NE.I1) THEN
          PRINT*,'DOUBLE COUNTING IN CLASSICAL_EXCLUSIONS'
        END IF
        NEXCLUSION=I1
      END IF
!
!     ==========================================================================
!     ==  ADD ONE NUMBER (=-1) TO CHECK FOR LAST EXCLUSION                    ==
!     ==========================================================================
      NEXCLUSION=NEXCLUSION+1
      IF(NEXCLUSION.LE.NEXCLUSIONX) THEN
        IEXCLUSION(1,NEXCLUSION)=-1
      END IF
!
!     ==================================================================
!     ==  MAKE AN ARRAY....                                           ==
!     ==================================================================
!         == IF ISVAR=IAT1-1+NAT*(IAT2-1)+...
!         == THEN 1) IT IS THE TRANSLATION OF IAT1
!         ==      2) IAT1.GE.IAT2
      I2FIRST(:)=0
      DO I=1,NEXCLUSION
        IAT2=1+INT(REAL(IEXCLUSION(1,I))/REAL(((1+2*MAXDIV)**3*NAT)))
        IF(I2FIRST(IAT2).EQ.0) I2FIRST(IAT2)=I
      ENDDO
      IF(I2FIRST(NAT).EQ.0)I2FIRST(NAT)=NEXCLUSION
      DO IAT2=NAT-1,1,-1
        IF(I2FIRST(IAT2).EQ.0)I2FIRST(IAT2)=I2FIRST(IAT2+1)
      ENDDO
! PRINT*,'I2FIRST',I2FIRST
! PRINT*,"I   I2FIRST(I)"
! DO I=1,NAT
!    PRINT*,I,"   ",I2FIRST(I)
! END DO
! DO  I=1,NEXCLUSION-1
!   C=IEXCLUSION(I)-1
!   IAT1=1+C/((1+2*MAXDIV)**3*NAT)
!   C=C-(IAT1-1)*(1+2*MAXDIV)**3*NAT
!   IAT2=1+C/(1+2*MAXDIV)**3
!   C=C-(IAT2-1)*(1+2*MAXDIV)**3
!   IT(3)=C/(1+2*MAXDIV)**2-MAXDIV
!   C=C-(IT(3)+MAXDIV)*(1+2*MAXDIV)**2
!   IT(2)=C/(1+2*MAXDIV)-MAXDIV
!   C=C-(IT(2)+MAXDIV)*(1+2*MAXDIV)
!   IT(1)=C-MAXDIV
!   ISVAR=IAT2-1+NAT*(IAT1-1)
!   ISVAR=(IT(1)+MAXDIV)+(1+2*MAXDIV)*((IT(2)+MAXDIV) &
!  &                        +(1+2*MAXDIV)*((IT(3)+MAXDIV) &
!  &                        +(1+2*MAXDIV)*ISVAR))
!   ISVAR=1+ISVAR
!   PRINT*,"EXCLUSION: ",I,IAT1,IAT2,IT(:),IEXCLUSION(I)
!   IF(IAT2.EQ.1+INT(REAL(IEXCLUSION(I)-1)/REAL(((1+2*MAXDIV)**3*NAT)))) THEN
!      CALL ERROR$STOP('--------')
!   END IF
!   IF(ISVAR.NE.IEXCLUSION(I)) THEN
!     CALL ERROR$STOP('.......')
!   END IF
! ENDDO
! !STOP 'FORCED'

!
!     ==================================================================
!     ==  PRINTOUT IF REQUESTED                                       ==
!     ==================================================================
      IF(TPR) THEN
        IMAX=MIN(NEXCLUSION,NEXCLUSIONX)
        PRINT*,'IEXCLUSION',IEXCLUSION(1,1:IMAX)
      END IF
      RETURN
      END                
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLASSICAL_REFLECTINGBOUNDARY(NAT,R,F,RCENT,FCENT &
     &                    ,RAD,ALPHA)
!     ******************************************************************
!     **                                                            **
!     **  ADDS THE FORCES OF A REFLECTING BOUNDARY                  **
!     **  THE BOUNDARY IS A SPHERE CENTERED AT RCENT WITH RADIUS R  **
!     **  THE WALLS HAVE A HARDNESS SPECIFIED BY ALPHA              **
!     **  THE TOTAL FORCES ACTING ON THE WALLS IS RETURNED IN FCENT **
!     **                                                            **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NAT
      REAL(8)   ,INTENT(IN)   :: R(3,NAT)
      REAL(8)   ,INTENT(INOUT):: F(3,NAT)
      REAL(8)   ,INTENT(IN)   :: RCENT(3)
      REAL(8)   ,INTENT(IN)   :: RAD
      REAL(8)   ,INTENT(IN)   :: ALPHA
      REAL(8)   ,INTENT(OUT)  :: FCENT(3)
      REAL(8)                 :: RAD2
      INTEGER(4)              :: IAT
      REAL(8)                 :: D2,D,FAC
      REAL(8)                 :: DR(3)
      REAL(8)                 :: DF(3)
!     ******************************************************************
      RAD2=RAD**2
      FCENT(:)=0.D0
      DO IAT=1,NAT
        DR(:)=R(:,IAT)-RCENT(:)
        D2=DOT_PRODUCT(DR,DR)
        IF(D2.GT.RAD2) THEN
          D=SQRT(D2)
          FAC=(D-RAD)/D*ALPHA
          DF(:)=-FAC*DR(:)
          FCENT(:)=FCENT(:)-DF(:)
          F(:,IAT)=F(:,IAT)+DF(:)
        END IF
      ENDDO
      RETURN
      END
!     ==================================================================
!     ==================================================================
!     ====  HIGHER LEVEL TOOLS TO INTERACT WITH CLASSICAL           ====
!     ==================================================================
!     ==================================================================
!     ==  CLASSICAL_BUILDSOLVENT                                      ==
!     ==  CLASSICAL_MINIMIZE                                          ==
! 
!     ..................................................................
      SUBROUTINE CLASSICAL$BUILDSOLVENT(DENSITY,RSPHERE &
     &                           ,IDENTEXCL,IDENTMOL,IDENTSOLV)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: IDENTEXCL
      CHARACTER(*),INTENT(IN)  :: IDENTSOLV
      CHARACTER(*),INTENT(IN)  :: IDENTMOL
      REAL(8)     ,INTENT(IN)  :: RSPHERE
      REAL(8)     ,INTENT(IN)  :: DENSITY
      REAL(8)     ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      LOGICAL(4)               :: TCHK
      REAL(8)                  :: RBAS(3,3)
      REAL(8)                  :: RBASINV(3,3)
      LOGICAL(4)  ,ALLOCATABLE :: TINDEX(:,:,:) !(N1,N2,N3)
      REAL(8)     ,ALLOCATABLE :: R0S(:,:)       !(3,NATS)
      REAL(8)     ,ALLOCATABLE :: QELS(:)        !(NATS)
      REAL(8)     ,ALLOCATABLE :: RMASSS(:)      !(NATS)
      CHARACTER(8),ALLOCATABLE :: TYPES(:)       !(NATS)
      INTEGER(4)  ,ALLOCATABLE :: IBONDSNEW(:,:)    !(6,NBONDS)
      REAL(8)     ,ALLOCATABLE :: BOS(:)         !(NBONDS)
      REAL(8)     ,ALLOCATABLE :: R0(:,:)        !(3,NAT)
      REAL(8)     ,ALLOCATABLE :: QEL(:)         !(NAT)
      REAL(8)     ,ALLOCATABLE :: RMASS(:)       !(NAT)
      CHARACTER(8),ALLOCATABLE :: TYPE(:)        !(NAT)
      INTEGER(4)  ,ALLOCATABLE :: IBONDNEW(:,:)!(6,NBOND)(IAT1,IAT2,IT2(3),IPOT)
      REAL(8)     ,ALLOCATABLE :: BO(:)          !(NBOND)
      REAL(8)     ,ALLOCATABLE :: REXCL(:,:)     !(3,NATEXCL)
      INTEGER(4)               :: NAT,NBOND
      INTEGER(4)               :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
      INTEGER(4)               :: N1,N2,N3
      INTEGER(4)               :: I1,I2,I3
      REAL(8)                  :: T1,T2,T3
      REAL(8)                  :: X,Y,Z,X1,Y1,Z1
      REAL(8)                  :: SVAR
      REAL(8)                  :: TWOPI
      INTEGER(4)               :: NMOL
      INTEGER(4)               :: NATS,NBONDS
      INTEGER(4)               :: IBS,IMOL,IAT,IATS,IB,I
      INTEGER(4)               :: NATEXCL,NBONDEXCL
      REAL(8)                  :: ALAT,DET
!     ******************************************************************
      TWOPI=2.D0*PI
!
!     ==================================================================
!     ==  COLLECT DATA FOR THE SOLVENT MOLECULE                       ==
!     ==================================================================
      CALL CLASSICAL$SELECT(IDENTMOL)
      CALL CLASSICAL$GETI4('NAT',NAT)
      CALL CLASSICAL$GETI4('NBOND',NBOND)
!
      ALLOCATE(R0(3,NAT))
      ALLOCATE(TYPE(NAT))
      ALLOCATE(QEL(NAT))
      ALLOCATE(RMASS(NAT))
      ALLOCATE(IBONDNEW(6,NBOND))
      ALLOCATE(BO(NBOND))
      CALL CLASSICAL$GETR8A('R(0)',3*NAT,R0)
      CALL CLASSICAL$GETCHA('TYPE',NAT,TYPE)
      CALL CLASSICAL$GETR8A('QEL',NAT,QEL)
      CALL CLASSICAL$GETR8A('MASS',NAT,RMASS)
      CALL CLASSICAL$GETI4A('BOND',5*NBOND,IBONDNEW)
      CALL CLASSICAL$GETR8A('BONDORDER',NBOND,BO)
!
!     ==================================================================
!     ==  DEFINE GRID                                                 ==
!     ==================================================================
      ALAT=(4.D0/DENSITY)**(1.D0/3.D0)
      RBAS(:,:)=0.5D0*ALAT
      DO I=1,3
        RBAS(I,I)=0.D0
      ENDDO
      CALL BOXSPH(RBAS,0.D0,0.D0,0.D0,RSPHERE &
     &           ,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
      N1=MAX1-MIN1+1
      N2=MAX2-MIN2+1
      N3=MAX3-MIN3+1
!
!     ==================================================================
!     ==  SPECIFY GRID POSITIONS IN THE SPHERE                        ==
!     ==================================================================
      ALLOCATE(TINDEX(N1,N2,N3))
      DO I1=MIN1,MAX1
        T1=DBLE(I1)
        DO I2=MIN2,MAX2 
          T2=DBLE(I2)
          X1=RBAS(1,1)*T1+RBAS(1,2)*T2
          Y1=RBAS(2,1)*T1+RBAS(2,2)*T2
          Z1=RBAS(3,1)*T1+RBAS(3,2)*T2
          DO I3=MIN3,MAX3
            T3=DBLE(I3)
            X=X1+RBAS(1,3)*T3
            Y=Y1+RBAS(2,3)*T3
            Z=Z1+RBAS(3,3)*T3
            SVAR=X**2+Y**2+Z**2
            TINDEX(I1-MIN1+1,I2-MIN2+1,I3-MIN3+1)=(SVAR.LT.RSPHERE**2)
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  EXCLUDE CENTRAL MOLECULE                                    ==
!     ==================================================================
      CALL CLASSICAL$SELECT(IDENTEXCL)
      CALL CLASSICAL$GETI4('NAT',NATEXCL)
      ALLOCATE(REXCL(3,NATEXCL))
      CALL CLASSICAL$GETR8A('R(0)',3*NATEXCL,REXCL)
      CALL GBASS(RBAS,RBASINV,DET)
      RBASINV(:,:)=TRANSPOSE(RBASINV)/TWOPI
      DO IAT=1,NATEXCL
        X=RBASINV(1,1)*REXCL(1,IAT) &
     &   +RBASINV(1,2)*REXCL(2,IAT) &
     &   +RBASINV(1,3)*REXCL(3,IAT)
        Y=RBASINV(2,1)*REXCL(1,IAT) &
     &   +RBASINV(2,2)*REXCL(2,IAT) &
     &   +RBASINV(2,3)*REXCL(3,IAT)
        Z=RBASINV(3,1)*REXCL(1,IAT) &
     &   +RBASINV(3,2)*REXCL(2,IAT) &
     &   +RBASINV(3,3)*REXCL(3,IAT)
        I1=NINT(X)
        I2=NINT(Y)
        I3=NINT(Z)
        TCHK=(I1.LT.MIN1).OR.(I2.LT.MIN2).OR.(I3.LT.MIN3) &
     &   .OR.(I1.GT.MAX1).OR.(I2.GT.MAX2).OR.(I3.GT.MAX3)
        IF(TCHK) THEN
          CALL ERROR$MSG('MOLECULE NOT FULLY SOLVATED')
          CALL ERROR$STOP('SOLVENT')
        ELSE
          TINDEX(I1-MIN1+1,I2-MIN2+1,I3-MIN3+1)=.FALSE.
        END IF
      ENDDO
      DEALLOCATE(REXCL)
!
!     ==================================================================
!     ==  COUNT NUMBER OF SOLVENT MOLECULES                           ==
!     ==================================================================
      NMOL=0
      DO I1=1,N1
        DO I2=1,N2
          DO I3=1,N3
            IF(TINDEX(I1,I2,I3))NMOL=NMOL+1
          ENDDO
        ENDDO
      ENDDO
      PRINT*,'NUMBER OF SOLVENT MOLECULES ',NMOL
!
!     ==================================================================
!     ==  ALLOCATE ARRAYS FOR THE COMBINED SYSTEM                     == 
!     == AND SET CENTRAL MOLECULESET BOND ARRAY FOR CLASSICAL         ==
!     ==================================================================
      CALL CLASSICAL$GETI4('NBOND',NBONDEXCL)
      NATS=NMOL*NAT+NATEXCL
      NBONDS=NMOL*NBOND+NBONDEXCL
!
      ALLOCATE(R0S(3,NATS))
      ALLOCATE(TYPES(NATS))
      ALLOCATE(RMASSS(NATS))
      ALLOCATE(QELS(NATS))
      ALLOCATE(IBONDSNEW(6,NBONDS))
      ALLOCATE(BOS(NBONDS))
!
      CALL CLASSICAL$GETR8A('R(0)',8*3*NATEXCL,R0S)
      CALL CLASSICAL$GETCHA('TYPE',NATEXCL,TYPES)
      CALL CLASSICAL$GETR8A('MASS',NATEXCL,RMASSS)
      CALL CLASSICAL$GETR8A('QEL',NATEXCL,QELS)
      CALL CLASSICAL$GETI4A('BOND',5*NBONDEXCL,IBONDSNEW)
      CALL CLASSICAL$GETR8A('BONDORDER',NBONDEXCL,BOS)
!
!     ==================================================================
!     ==  SET BOND ARRAY FOR CLASSICAL                                ==
!     ==================================================================
      IBS=NBONDEXCL
      DO IMOL=1,NMOL
        DO IB=1,NBOND
          IBS=IBS+1
          IBONDSNEW(1,IBS)=IBONDNEW(1,IB)+NAT*(IMOL-1)+NATEXCL
          IBONDSNEW(2,IBS)=IBONDNEW(2,IB)+NAT*(IMOL-1)+NATEXCL
          IBONDSNEW(3:5,IBS)=IBONDNEW(3:5,IB)
          BOS(IBS)=BO(IB)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  SET BOND POSITIONS AND ATOMIC QUANTITIES                    ==
!     ==================================================================
      IATS=NATEXCL
      DO I1=MIN1,MAX1
        T1=DBLE(I1)
        DO I2=MIN2,MAX2 
          T2=DBLE(I2)
          DO I3=MIN3,MAX3
            T3=DBLE(I3)
            X=RBAS(1,1)*T1+RBAS(1,2)*T2+RBAS(1,3)*T3
            Y=RBAS(2,1)*T1+RBAS(2,2)*T2+RBAS(2,3)*T3
            Z=RBAS(3,1)*T1+RBAS(3,2)*T2+RBAS(3,3)*T3
            IF(TINDEX(I1-MIN1+1,I2-MIN2+1,I3-MIN3+1)) THEN
              DO IAT=1,NAT
                IATS=IATS+1
                R0S(1,IATS)=R0(1,IAT)+X
                R0S(2,IATS)=R0(2,IAT)+Y
                R0S(3,IATS)=R0(3,IAT)+Z
                QELS(IATS)=QEL(IAT)
                RMASSS(IATS)=RMASS(IAT)
                TYPES(IATS)=TYPE(IAT)
              ENDDO
            END IF
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  SET DATA IN THE CLASSICAL OBJECT FOR THE COMBINED SYSTEM    ==
!     ==================================================================
      CALL CLASSICAL$SELECT(IDENTSOLV)
      CALL CLASSICAL$SETI4('NAT',NATS)
      CALL CLASSICAL$SETR8A('R(0)',3*NATS,R0S)
      CALL CLASSICAL$SETR8A('R(-)',3*NATS,R0S)
      CALL CLASSICAL$SETR8A('QEL',NATS,QELS)
      CALL CLASSICAL$SETCHA('TYPE',NATS,TYPES)
      CALL CLASSICAL$SETR8A('MASS',NATS,RMASSS)
      CALL CLASSICAL$SETI4('NBOND',NBONDS)
      CALL CLASSICAL$SETI4A('BOND',5*NBONDS,IBONDSNEW)
      CALL CLASSICAL$SETR8A('BONDORDER',NBONDS,BOS)
      DEALLOCATE(IBONDSNEW)
      DEALLOCATE(BOS)
      DEALLOCATE(R0S)
      DEALLOCATE(QELS)
      DEALLOCATE(RMASSS)
      DEALLOCATE(TYPES)
      DEALLOCATE(TINDEX)
      DEALLOCATE(R0)
      DEALLOCATE(TYPE)
      DEALLOCATE(QEL)
      DEALLOCATE(RMASS)
      DEALLOCATE(IBONDNEW)
      DEALLOCATE(BO)
      RETURN
      END       
! 
!     ..................................................................
      SUBROUTINE CLASSICAL$MINIMIZE(TOL_,MAXSTEP,TCONV)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
      REAL(8)   ,INTENT(IN) :: TOL_
      INTEGER(4),INTENT(IN) :: MAXSTEP
      LOGICAL(4),INTENT(OUT):: TCONV
      INTEGER(4)            :: NITER
      REAL(8)               :: TOL
      INTEGER(4)            :: NAT
      INTEGER(4)            :: NFILO
      REAL(8)   ,ALLOCATABLE:: RMASS(:)  !(NAT)
      REAL(8)   ,ALLOCATABLE:: RMASS1(:)  !(NAT)
      REAL(8)               :: EPOTP,ETOTP
      REAL(8)               :: ANNE
      REAL(8)               :: FMAX
      REAL(8)               :: ECON,EKIN,EPOT
      REAL(8)               :: DELT
      INTEGER(4)            :: IAT,ITER
!     ******************************************************************
!     == MAXIMUM NUMBER OF ITERATIONS ==================================
      IF(MAXSTEP.EQ.0) THEN
        NITER=10000
      ELSE
        NITER=MAXSTEP
      END IF
!     == ITERATION IS STOPED IF FMAX IS SMALLER THAN TOL ===============
      IF(TOL_.EQ.0.D0) THEN
        TOL=1.D-4
      ELSE 
        TOL=TOL_
      END IF
      DELT=10.D0
!     ==  SET MASSES ABOUT EQUAL TO THE PROTON MASS ====================
      CALL CLASSICAL$GETI4('NAT',NAT)
      ALLOCATE(RMASS(NAT))
      CALL CLASSICAL$GETR8A('MASS',NAT,RMASS)
      ALLOCATE(RMASS1(NAT))
      DO IAT=1,NAT
!       RMASS1(IAT)=2.D0*1899.D0
      ENDDO
!     CALL CLASSICAL$SETR8A('MASS',NAT,RMASS1)
!
      ETOTP=1.D+8
      DO ITER=1,NITER
!       ================================================================
!       ==  PROPAGATE                                                 ==
!       ================================================================
        CALL CLASSICAL$ETOT(EPOT)
        ANNE=0.D0
        IF(EPOT.GT.EPOTP) ANNE=1.0D0
        EPOTP=EPOT
        CALL CLASSICAL$PROPAGATE(DELT,ANNE)
!
!       ================================================================
!       ==  CHECK CONVERGENCE                                         ==
!       ================================================================
        CALL CLASSICAL$MAXFORCE(FMAX)
!
!       ================================================================
!       ==  PRINT                                                     ==
!       ================================================================
        IF(TPR) THEN
          CALL FILEHANDLER$UNIT('PROT',NFILO)
          CALL CLASSICAL$EKIN(DELT,EKIN)
          ECON=EPOT+EKIN
          WRITE(NFILO,FMT='("I",I10,"EPOT",E12.5," EKIN",E12.5' &
     &                 //'," ECONS",E12.5," ANNE",F3.1,"FMAX",E12.5)') &
     &         ITER,EPOT,EKIN,ECON,ANNE,FMAX
        END IF
!
!       ================================================================
!       ==  CHECK CONVERGENCE                                         ==
!       ================================================================
        TCONV=(FMAX.LT.TOL)
        IF(TCONV) EXIT
!
!       ================================================================
!       == SWITCH                                                     ==
!       ================================================================
        CALL CLASSICAL$SWITCH
      ENDDO
      PRINT*,'CLASSICAL$MINIMIZE CONVERGED AFTER ',ITER,' ITERATIONS'
      CALL CLASSICAL$SETR8A('MASS',NAT,RMASS)
      DEALLOCATE(RMASS)
      DEALLOCATE(RMASS1)
      RETURN
      END   
! 
!     ..................................................................
      SUBROUTINE CLASSICAL$D2EDX2(NAT_,DYNMAT)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT_
      REAL(8)   ,INTENT(OUT):: DYNMAT(3,NAT_,3,NAT_)
      INTEGER(4)            :: NAT
      REAL(8)   ,ALLOCATABLE:: RSAVE(:,:)  !(3,NAT)
      REAL(8)   ,ALLOCATABLE:: R0(:,:)     !(3,NAT)
      REAL(8)   ,ALLOCATABLE:: FORCE1(:,:) !(3,NAT)
      REAL(8)   ,ALLOCATABLE:: FORCE2(:,:) !(3,NAT)
      REAL(8)   ,PARAMETER  :: DELTA=1.D-5
      INTEGER(4)            :: I,IAT
      INTEGER(4)            :: IND1,IND2,IAT1,IAT2,I1,I2
      REAL(8)               :: EPOT
      REAL(8)               :: FAC,SVAR
!     ******************************************************************
      CALL CLASSICAL$GETI4('NAT',NAT)
      IF(NAT.NE.NAT_) THEN
        CALL ERROR$MSG('NAT OF CLASSICAL NOT CONSISTENT WITH INPUT')
        CALL ERROR$STOP('CLASSICAL$DYNAMICALMATRIX')
      END IF
      ALLOCATE(RSAVE(3,NAT))
!
!     ==================================================================
!     == CALCULATED D2E/DX2                                           ==
!     ==================================================================
      ALLOCATE(R0(3,NAT))
      ALLOCATE(FORCE1(3,NAT))
      ALLOCATE(FORCE2(3,NAT))
      CALL CLASSICAL$GETR8A('R(0)',3*NAT,RSAVE)
      DO IAT=1,NAT
        DO I=1,3
          R0(:,:)=RSAVE(:,:)
          R0(I,IAT)=R0(I,IAT)-DELTA
          CALL CLASSICAL$SETR8A('R(0)',3*NAT,R0)
          CALL CLASSICAL$ETOT(EPOT)
          CALL CLASSICAL$GETR8A('FORCE',3*NAT,FORCE1)
          R0(:,:)=RSAVE(:,:)
          R0(I,IAT)=R0(I,IAT)+DELTA
          CALL CLASSICAL$SETR8A('R(0)',3*NAT,R0)
          CALL CLASSICAL$ETOT(EPOT)
          CALL CLASSICAL$GETR8A('FORCE',3*NAT,FORCE2)
          FAC=-1.D0/(2.D0*DELTA)
          DYNMAT(:,:,I,IAT)=FAC*(FORCE2(:,:)-FORCE1(:,:))
        ENDDO
      ENDDO
!
!     == SYMMETRIZE DYNMAT =============================================
      DO IND1=1,3*NAT
        IAT1=(IND1-1)/3+1
        I1=IND1-3*(IAT1-1)
        DO IND2=I1+1,3*NAT
          IAT2=(IND2-1)/3+1
          I2=IND2-3*(IAT2-1)
          SVAR=0.5D0*(DYNMAT(I1,IAT1,I2,IAT2)+DYNMAT(I2,IAT2,I1,IAT1))
          DYNMAT(I1,IAT1,I2,IAT2)=SVAR 
          DYNMAT(I2,IAT2,I1,IAT1)=SVAR
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  RESTORE ORIGINAL STATE                                      ==
!     ==================================================================
      CALL CLASSICAL$SETR8A('R(0)',3*NAT,RSAVE)
      CALL CLASSICAL$ETOT(EPOT)
      DEALLOCATE(FORCE2)
      DEALLOCATE(FORCE1)
      DEALLOCATE(R0)
      DEALLOCATE(RSAVE)
      RETURN
      END   
!     ==================================================================
!     ==================================================================
!     ====    ROUTINES TO TEST CLASSICAL   =============================
!     ==================================================================
!     ==================================================================
!
!     ..................................................................
      SUBROUTINE TESTVALUE(POT)
!     ******************************************************************
!     ******************************************************************
      USE CLASSICAL_MODULE,ONLY: POT_TYPE
      IMPLICIT NONE
      TYPE(POT_TYPE),INTENT(IN)  :: POT
      INTEGER(4),PARAMETER :: NDIFF=10
      REAL(8)              :: X(-NDIFF:NDIFF)
      REAL(8)              :: E(-NDIFF:NDIFF)
      REAL(8)              :: DEDX(-NDIFF:NDIFF)
      INTEGER(4)           :: I0  
      INTEGER(4)           :: I
      REAL(8)              :: X0
      REAL(8)              :: DX1
      REAL(8)              :: SVAR
!     ******************************************************************
      I0=POT%NX/2
      X0=POT%X1+POT%DX*(I0-1)
      DX1=POT%DX/DBLE(NDIFF)
      DO I=-NDIFF,NDIFF
        X(I)=X0+DX1*DBLE(I)
        CALL VALUE(POT,X(I),E(I),DEDX(I))
      ENDDO
      PRINT*,'X1,DX,NX,E(0) ',POT%X1,POT%DX,POT%NX,X0,E(0)
      DO I=-NDIFF,NDIFF
        E(I)=E(I)-E(0)
      ENDDO
      DO I=-NDIFF,NDIFF
        IF(I.GT.-NDIFF.AND.I.LT.NDIFF) THEN
          SVAR=(E(I+1)-E(I-1))/(2.D0*DX1)
        ELSE
          SVAR=0.D0
        END IF
        WRITE(*,FMT='("I=",I5,"X=",E12.5,"E=",E12.5,"DEDX=",E12.5' &
     &          //',"TEST=",E12.5)')I,X(I),E(I),DEDX(I),SVAR
      ENDDO
      CALL ERROR$MSG('NORMAL STOP IN TESTVALUE')
      CALL ERROR$STOP('TESTVALUE')
      END
!
!     ..................................................................
      SUBROUTINE TESTFORCES
!     ****************************************************************** 
!     ****************************************************************** 
      IMPLICIT NONE
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R0SAVE(:,:)
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      REAL(8)   ,ALLOCATABLE :: FSAVE(:,:)
      REAL(8)   ,PARAMETER   :: DEL=1.D-2
      REAL(8)                :: EPOTSAVE,EPOTM2,EPOTM1,EPOTP1,EPOTP2
      REAL(8)                :: D0,D1,D2
      REAL(8)                :: DPRED
      INTEGER(4)             :: IAT,I
!     ****************************************************************** 
      PRINT*,'=======FORCETEST======'
      CALL CLASSICAL$GETI4('NAT',NAT)
      ALLOCATE(R0SAVE(3,NAT))
      ALLOCATE(R0(3,NAT))
      ALLOCATE(FSAVE(3,NAT))
!
      CALL CLASSICAL$GETR8A('R(0)',3*NAT,R0SAVE)
      CALL CLASSICAL$ETOT(EPOTSAVE)
      CALL CLASSICAL$GETR8A('FORCE',3*NAT,FSAVE)
      DO IAT=1,NAT 
        DO I=1,3
!
          R0=R0SAVE
          R0(I,IAT)=R0(I,IAT)-2*DEL
          CALL CLASSICAL$SETR8A('R(0)',3*NAT,R0)
          CALL CLASSICAL$ETOT(EPOTM2)
!
          R0=R0SAVE
          R0(I,IAT)=R0(I,IAT)-DEL
          CALL CLASSICAL$SETR8A('R(0)',3*NAT,R0)
          CALL CLASSICAL$ETOT(EPOTM1)
!
          R0=R0SAVE
          R0(I,IAT)=R0(I,IAT)+DEL
          CALL CLASSICAL$SETR8A('R(0)',3*NAT,R0)
          CALL CLASSICAL$ETOT(EPOTP1)
!
          R0=R0SAVE
          R0(I,IAT)=R0(I,IAT)+2*DEL
          CALL CLASSICAL$SETR8A('R(0)',3*NAT,R0)
          CALL CLASSICAL$ETOT(EPOTP2)
!
          EPOTM2=EPOTM2-EPOTSAVE
          EPOTM1=EPOTM1-EPOTSAVE
          EPOTP1=EPOTP1-EPOTSAVE
          EPOTP2=EPOTP2-EPOTSAVE
          D2=(EPOTP2-EPOTM2)/(4*DEL)
          D1=(EPOTP1-EPOTM1)/(2*DEL)
          D0=-FSAVE(I,IAT)
          DPRED=(4.D0*D1-D2)/3.D0
          WRITE(*,FMT='("==",2I3,5E15.5)')IAT,I,D0,DPRED,D1,D2
        ENDDO
      ENDDO      
!
!     ==================================================================
!     == RESTORE ORIGINAL STATE                                       == 
!     ==================================================================
      CALL CLASSICAL$SETR8A('R(0)',3*NAT,R0SAVE)
      CALL CLASSICAL$ETOT(EPOTSAVE)
      DEALLOCATE(R0SAVE)
      DEALLOCATE(R0)
      DEALLOCATE(FSAVE)
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE LENNARDJONES(EPSLJ,SIGLJ,POT)
!     ******************************************************************
!     **                                                              **
!     ** LENNARD JONES POTENTIAL IN TERMS OF INVERSE DISTANCE         **
!     **                                                              **
!     ******************************************************************
      USE CLASSICAL_MODULE,ONLY: POT_TYPE
      IMPLICIT NONE
      REAL(8)       ,INTENT(IN)  :: EPSLJ
      REAL(8)       ,INTENT(IN)  :: SIGLJ
      TYPE(POT_TYPE),INTENT(OUT) :: POT
      REAL(8)                    :: SIGLJ6,FACLJ1,FACLJ2,X,X5,SVAR
      INTEGER(4)                 :: I
!     ******************************************************************
      POT%NX=1000
      ALLOCATE(POT%VAL(POT%NX))
      ALLOCATE(POT%DER(POT%NX))
      POT%X1=0.D0
      POT%DX=1.3D0*SIGLJ/DBLE(POT%NX-1)
      SIGLJ6=SIGLJ**6
      FACLJ1=4.D0*EPSLJ
      FACLJ2=FACLJ1*SIGLJ6*6.D0
      DO I=1,POT%NX
        X=POT%X1+POT%DX*DBLE(I-1)
        X5=X**5
        SVAR=SIGLJ6*X5*X
        POT%VAL(I)=FACLJ1*(SVAR**2-SVAR)
        POT%DER(I)=FACLJ2*(2*SVAR-1.D0)*X5
      ENDDO
      RETURN
      END       
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ANGLEARRAY(STRENGTH,PHI0,POT)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE ANGLE ARRAY IN TERMS OF THE COSINE OF THE ANGLE   **
!     **                                                              **
!     **    VARRAY=0.5D0*STRENGTH*(PHI-PHI0)**2;  PHI=ACOS(X)         **
!     **                                                              **
!     ******************************************************************
      USE CLASSICAL_MODULE,ONLY: POT_TYPE
      IMPLICIT NONE
      REAL(8)       ,INTENT(IN) :: STRENGTH
      REAL(8)       ,INTENT(IN) :: PHI0
      TYPE(POT_TYPE),INTENT(OUT):: POT
      REAL(8)       ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)                   :: XI,PHI,DPHIDX
      INTEGER(4)                :: IX
!     ******************************************************************
      POT%NX=1000
      ALLOCATE(POT%VAL(POT%NX))
      ALLOCATE(POT%DER(POT%NX))
      POT%X1=-1.D0
      POT%DX=2.D0/DBLE(POT%NX-1)
      DO IX=2,POT%NX-1
        XI=POT%X1+POT%DX*DBLE(IX-1)
        PHI=DACOS(XI)
        DPHIDX=-1.D0/SIN(PHI)
        POT%VAL(IX)=0.5D0*STRENGTH*(PHI-PHI0)**2
        POT%DER(IX)=STRENGTH*(PHI-PHI0)*DPHIDX
      ENDDO
      POT%VAL(1)=POT%VAL(2)+0.5D0*POT%DER(2)*POT%DX
      POT%DER(1)=0.D0
      POT%VAL(POT%NX)=POT%VAL(POT%NX)+0.5D0*POT%DER(POT%NX-1)*POT%DX
      POT%DER(POT%NX)=0.D0
      RETURN 
      END
