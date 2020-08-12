!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE BRILLOUIN_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: BRILLOUIN                                                          **
!**                                                                           **
!**  PURPOSE: BRILLOUIN ZONE INTEGRATION USING THE IMPROVED VERSION           **
!**    OF THE TETRAHEDRON METHOD                                              **
!**                                                                           **
!**  DESCRIPTION:                                                             **
!**    BRILLOUIN$MSH CREATE INFORMATION WHICH IS ONLY DEPENDENT               **
!**    ON THE STRUCTURE OF THE CRYSTAL. THIS FUNCTION IS EXECUTED             **
!**    BEFORE THE SELFCONSITENCY ITERATIONS. THE IRREDUCIBLE K-POINTS         **
!**    ARE THEN AVAILABLE THROUGH THE $GET.. INTERFACE                        **
!**                                                                           **
!**    BRILLOUIN$DOS IS EXECUTED IN EACH SELFCONSISTENCY ITERATION.           **
!**    THEY NEED AS INPUT THE ENERGY EIGENVALUES EB AND THE NUMBER OF         **
!**    OCCUPIED STATES RNTOT. FOR SPIN POLARIZED CALCULATIONS TREAT           **
!**    SPIN UP AND SPIN DOWN AS SEPARATE BANDS. OUTPUT ARE THE                **
!**    WEIGHTS WGHT. THE INTEGRATION OF AN ARBITRARY FUNCTION A(N,K)          **
!**    OVER THE OCCUPIED STATES IS PERFORMED BY SUMMATION                     **
!**            <A>=SUM OVER K AND N OF WGHT(N,K)*A(N,K)                       **
!**    WHERE THE SUM RUNS OVER OCCUPIED AND! UNOCCUPIED STATES.               **
!**    THE WEIGHTS CONTAIN BOTH THE GEOMETRICAL WEIGHT AND THE                **
!**    INFLUENCE OF THE FERMI FUNCTION.                                       **
!**                                                                           **
!**    A) FOR INSULATORS THE METHOD IS IDENTICAL TO THE SPECIAL POINT         **
!**       SCHEME OF MONKHORST AND PACK.                                       **
!**    B) FOR METALS IT IS IDENTICAL TO THE "TRADITIONAL" TETRAHEDRON         **
!**       METHOD OF ANDERSEN AND JEPSEN, IF THE CORRECTION IS SWITCHED        **
!**       OFF.                                                                **
!**    C) WITH THE CORRECTION ALSO THE RESULTS FOR METALS ARE                 **
!**       COMPARABLE TO THAT OBTAINED FOR INSULATORS                          **
!**                                                                           **
!**    THE CORRECTION FORMULA FOR LINEAR INTERPOLATION CAN BE SWITCHED        **
!**    ON AND OFF BY THE PARAMETER "ICOR" IN THE SUBROUTINE "WEIGHTS"         **
!**                                                                           **
!**    THE SYMMETRY OPERATIONS USED AS INPUT ARE TABULATED IN:                **
!**      C.J.BRADLEY AND A.P.CRACKNELL,                                       **
!**      THE MATHEMATICAL THEORY OF SYMMETRY IN SOLIDS,                       **
!**      OXFORD 1972                                                          **
!**                                                                           **
!**    SOMETIMES THE ROUTINE HAS PROBLEMS FINDING THE FERMI LEVEL,            **
!**    IF THE FERMI LEVEL IS PINNED A TETRAHEDRON WITH IDENTICAL              **
!**    ENERGIES ON ALL 4 CORNERS, WHICH LEADS TO A DELTA PEAK IN THE          **
!**    DENSITIES OF STATES. IN THIS CASE ONE SHOULD INCREASE THE K-MESH,      **
!**    OR CHANGE THE K-MESH FROM AN EVEN NUMBER OF DIVISIONS FOR THE          **
!**    RECIPROCAL LATTICE VECTORS TO AN ODD NUMBER OR VICE VERSA.             **
!**                                                                           **
!**    THE FERMI LEVEL IS DETERMINED TO AN ACCURACY OF 1.E-5 IN THE           **
!**    NUMBER OF STATES. THIS TOLERANCE CAN BE MODIFIED BY CHANGING           **
!**    THE SETTING OF THE PARAMETER "TOLMAX" IN SUBROUTINE DOS.               **
!**                                                                           **
!**    THE GRID OF K-POINTS MAY OR MAY NOT INCLUDE THE GAMMA POINT            **
!**    SWITCH USING TSHIFT IN BRILLOUIN_REDUZ                                 **
!**                                                                           **
!**  FUNCTIONS:                                                               **
!**    BRILLOUIN$MSH                                                          **
!**    BRILLOUIN$DOS                                                          **
!**    BRILLOUIN$GETI4                                                        **
!**    BRILLOUIN$GETR8A                                                       **
!**                                                                           **
!**  USAGE:                                                                   **
!**    1) DEFINE LATTICE AND SYMMETRY OPERATIONS USING BRILLOUIN$MSH.         **
!**    2) OBTAIN THE IRREDUCIBLE K-POINTS USING                               **
!**           BRILLOUIN$GETI4('NK',NK)                                        **
!**           BRILLOIIN$GETR8A('K',3*NK,K)                                    **
!**    3) CALCULATE ONE-PARTICLE ENERGIES AT THE K-POINTS                     **
!**    4) OBTAIN SAMPLING WEIGHTS USING BRILLOUIN$DOS                         **
!**                                                                           **
!*******************************************************************************
TYPE THIS_TYPE
  INTEGER(4)             :: NKP       ! #(IRREDUCIBLE KPOINTS)
  INTEGER(4)             :: NTET      ! #(IRREDUCIBLE TETRAHEDRA)
  REAL(8)                :: VOL       ! 1/#(GENERAL TETRAHEDRA)
  REAL(8)                :: RBAS(3,3) ! REAL SPACE LATTICE VECTORS
  INTEGER(4)             :: NKDIV(3)
  INTEGER(4)             :: ISHIFT(3)
  REAL(8)       ,POINTER :: XK(:,:)   !(3,NKP) IRR. K-POINTS IN RELATIVE COORDINATES
  INTEGER(4)    ,POINTER :: IKP(:,:)  !(4,NTET) TETRAHEDRON CORNERS
  INTEGER(4)    ,POINTER :: MULT(:)   !(NTET) MULTIPLICITY OF THE TETRAHEDRON
  INTEGER(4)    ,POINTER :: IRRKP(:)  !(NMSHP) POINTER TO IRR. K-POINT
END TYPE THIS_TYPE
!
! THE EWGHT STRUCTURE KEEPS DATA ON THE RELEVANT PROTION OF THE ENERGY GRID.
! THE FIRST AND LAST GRID POINTS ARE SPECIFIED BY I1 AND I2.
TYPE EWGHT_TYPE
  INTEGER(4)             :: I1
  INTEGER(4)             :: I2
  REAL(8),POINTER        :: WGHT(:)
END TYPE EWGHT_TYPE
TYPE(THIS_TYPE) :: THIS
END MODULE BRILLOUIN_MODULE
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!----  BLOCK TEST                                                           ----
!----  TEST SUBROUTINE                                                      ----
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$TESTING
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NSYMX=50
      INTEGER(4),PARAMETER :: NBX=2
      REAL(8)              :: RBAS(3,3)
      INTEGER(4)           :: IARB(3)
      INTEGER(4)           :: IIO(3,3,NSYMX)
      REAL(8)   ,ALLOCATABLE:: BK(:,:)
      REAL(8)   ,ALLOCATABLE:: EB(:,:)
      REAL(8)   ,ALLOCATABLE:: WGHT(:,:)
      REAL(8)   ,ALLOCATABLE:: A(:,:)
      INTEGER(4)           :: NFIL1=1001
      REAL(8)              :: RNTOT
      REAL(8)              :: EF
      REAL(8)              :: SUMA
      INTEGER(4)           :: NSYM ! #(SYMMETRY OPERATIONS)
      INTEGER(4)           :: NKP ! #(K-POINTS)
      INTEGER(4)           :: ISYM,I,J,NB,IKP,IB
!     **************************************************************************
                                    CALL TRACE$PUSH('BRILLOUIN$TESTING')
!     ==========================================================================
!     ==  READ LATTICE VECTORS                                                ==
!     ==========================================================================
      OPEN(UNIT=NFIL1,FILE='TEST.IN')
      REWIND NFIL1
!     ==== RBAS(I,J) = REAL SPACE LATTICE VECTORS (I=X,Y,Z)
      DO I=1,3
        READ(NFIL1,*)RBAS(:,I)
      ENDDO
!     ==========================================================================
!     ==  INPUT OF SYMMETRY OPERATIONS                                        ==
!     ==========================================================================
      READ(NFIL1,*)IARB(:)
!     ==  NSYM = NUMBER OF SYMMETRY OPERATIONS
      READ(NFIL1,*)NSYM
!     ==  SYMMETRY MATRIX
      DO ISYM=1,NSYM
        READ(NFIL1,*)(IIO(:,J,ISYM),J=1,3)
        DO I=1,3
          DO J=1,3
            IIO(I,J,ISYM+NSYM)=-IIO(I,J,ISYM)
          ENDDO
        ENDDO
      ENDDO
      NSYM=NSYM*2
!     == NKP = NUMBER OF K-POINTS IN THE WHOLE UNIT CELL
      READ(NFIL1,*)NKP
!
!     ==========================================================================
!     ==  FIND IRREDUCIBLE K-POINTS AND TETRAHEDRA                            ==
!     ==========================================================================
!      CALL BRILLOUIN$MSH(RBAS,NKP,NSYM,IIO,IARB)

      CALL BRILLOUIN$MSHNOSYM(.TRUE.,RBAS,(/5,5,5/),(/1,1,1/))
 
!     ==========================================================================
!     ==  CALCULATE ENERGIES AT THE IRREDUCIBLE K-POINTS                      ==
!     ==========================================================================
      CALL BRILLOUIN$GETI4('NK',NKP)
      ALLOCATE(BK(3,NKP))
      CALL BRILLOUIN$GETR8A('K',3*NKP,BK)
!     ====  NB  = NUMBER OF BANDS
      NB=1
!     ====  RNTOT = NUMBER OF OCCUPIED STATES
      RNTOT=0.5D0
      PRINT*,'NKP ',NKP
!     ====  EB(IB,IKP) = ENERGY OF BAND IB AT K-POINT IKP
!     ====  F(IB,IKP) = MATRIX ELEMENT OF BAND IB AT K-POINT IKP
      ALLOCATE(EB(NBX,NKP))
      ALLOCATE(WGHT(NBX,NKP))
      ALLOCATE(A(NBX,NKP))
      DO IKP=1,NKP
        DO IB=1,NB
          EB(IB,IKP)=COS(BK(1,IKP))+COS(BK(2,IKP))+COS(BK(3,IKP))
          A(IB,IKP)=EB(IB,IKP)
          WRITE(6,FMT="(2I5,' BK ',3F10.4,' EB ',F10.4,' A ',F10.4)") &
     &            IKP,IB,(BK(I,IKP),I=1,3),EB(IB,IKP),A(IB,IKP)
        ENDDO
      ENDDO
      DO IKP=1,NKP
        DO IB=NB+1,NBX
          EB(IB,IKP)=5.D0
          A(IB,IKP)=EB(IB,IKP)
          WRITE(6,FMT="(2I5,' BK ',3F10.4,' EB ',F10.4,' A ',F10.4)") &
     &            IKP,IB,(BK(I,IKP),I=1,3),EB(IB,IKP),A(IB,IKP)
        ENDDO
      ENDDO
 
!     ==========================================================================
!     ==  CALCULATE WEIGHTS                                                   ==
!     ==========================================================================
      CALL BRILLOUIN$DOS(NB,NKP,EB,WGHT,RNTOT,EF)
 
!     ==========================================================================
!     ==  PERFORM BRILLOUIN ZONE INTEGRATION OF F(K)                          ==
!     ==========================================================================
      SUMA=0.D0
      DO IKP=1,NKP
        DO IB=1,NB
          SUMA=SUMA+WGHT(IB,IKP)*A(IB,IKP)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  PRINT OUT                                                           ==
!     ==========================================================================
      PRINT*,'INTEGRAL OF A : ',SUMA 
!
!     ==========================================================================
!     ==  TEST DERIVATIVES                                                    ==
!     ==========================================================================
      CALL BRILLOUIN_TESTWEIGHTANDDER

                                    CALL TRACE$POP
      END SUBROUTINE BRILLOUIN$TESTING
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_TESTWEIGHTANDDER
!     **************************************************************************
!     ** TEST ROUTINE FOR BRILLOUIN$WEIGHTANDER                               **
!     ** CHECKS DERIVATIVES AGAINST NUMERICAL DERIVATIVES                     **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8) :: VOL=1.D0
      REAL(8) :: E0(4)=(/2.D0,1.D0,0.D0,3.D0/)
      REAL(8) :: DELTA=1.D-8
      REAL(8) :: EF
      REAL(8) :: E(4)
      REAL(8) :: WGHT(4)
      REAL(8) :: WGHT0(4)
      REAL(8) :: DWGHT(4,4)
      REAL(8) :: DWGHT0(4,4)
      REAL(8) :: WGHTARR(4,2,4)
      REAL(8) :: DWGHTTEST(4,4)
      INTEGER(4) :: I,J,K
!     **************************************************************************
      DO K=1,3
        EF=0.3+REAL(K-1,KIND=8)
        PRINT*,'NOW TEST WEIGHTANDDER CASE ',K,EF
        CALL BRILLOUIN_WEIGHTANDDER(VOL,E0,EF,WGHT0,DWGHT0)
        DO I=1,4
          DO J=1,2
            E=E0
            E(I)=E(I)+REAL(2*J-3,KIND=8)*DELTA
            CALL BRILLOUIN_WEIGHTANDDER(VOL,E,EF,WGHT,DWGHT)
            WGHTARR(:,J,I)=WGHT
          ENDDO
          DWGHTTEST(:,I)=(WGHTARR(:,2,I)-WGHTARR(:,1,I))/(2.D0*DELTA)
          DO J=1,4
            PRINT*,'TEST ',DWGHT(J,I),DWGHTTEST(J,I),DWGHTTEST(J,I)-DWGHT(J,I)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_TESTCOMPLETE()
!     **************************************************************************
!     ** TESTROUTINE FOR THE ROUTINE BRILLOUIN_COMPLETE                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NOPX=48
      INTEGER(4)           :: NOP
      INTEGER(4)           :: OP(3,3,NOPX)
      INTEGER(4),DIMENSION(3,3) :: C2Z,C2X,C2A,C31P,INV
      INTEGER(4)           :: I
!     **************************************************************************
!     == LIST OF GENERATORS ====================================================
!     == GROUP 228 ==GAMMA^F_C==================================================
      C2Z(:,1)=(/0,1,0/)   ; C2Z(:,2)=(/1,0,0/)  ; C2Z(:,3)=(/-1,-1,-1/) 
      C2X(:,1)=(/-1,-1,-1/); C2X(:,2)=(/0,0,1/)  ; C2X(:,3)=(/0,1,0/)
      C2A(:,1)=(/-1,0,0/)  ; C2A(:,2)=(/0,-1,0/) ; C2A(:,3)=(/1,1,1/)
      C31P(:,1)=(/0,1,0/)  ; C31P(:,2)=(/0,0,1/) ; C31P(:,3)=(/1,0,0/)
      INV(:,1)=(/-1,0,0/)  ; INV(:,2)=(/0,-1,0/) ; INV(:,3)=(/0,0,-1/)
!
!     ==========================================================================
!     == PUT GENERATORS ONTO INPUT ARRAY                                      ==
!     ==========================================================================
      NOP=5
      OP(:,:,1)=C2Z
      OP(:,:,2)=C2X
      OP(:,:,3)=C2A
      OP(:,:,4)=C31P
      OP(:,:,5)=INV
      WRITE(*,FMT='(82("="),T10," GENERATORS OF THE GROUP ")')
      DO I=1,NOP
        WRITE(*,FMT='(I5,T20,3("|",3I5,"|"))')I,OP(:,:,I)
      ENDDO
!
!     ==========================================================================
!     ==  RUN BRILLOUIN_COMPLETE
!     ==========================================================================
      CALL BRILLOUIN_COMPLETE(NOPX,NOP,OP)
!
!     ==========================================================================
!     ==  RUN BRILLOUIN_COMPLETE
!     ==========================================================================
      WRITE(*,FMT='(82("="),T10," FULL POINT GROUP ")')
      DO I=1,NOP
        WRITE(*,FMT='(I5,T20,3("|",3I5,"|"))')I,OP(:,:,I)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$SETR8A(ID,LENG,VAL)
!     **************************************************************************
!     **************************************************************************
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LENG
      REAL(8)     ,INTENT(IN) :: VAL(LENG)
!     **************************************************************************
      IF(ID.EQ.'RBAS')THEN
        IF(LENG.NE.9) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BRILLOUIN$SETR8A')
        END IF
        THIS%RBAS=RESHAPE(VAL,(/3,3/))
      ELSE
        CALL ERROR$MSG('UNKNOWN ID')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BRILLOUIN$SETR8A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$GETR8A(ID,LENG,VAL)
!     **************************************************************************
!     **************************************************************************
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LENG
      REAL(8)     ,INTENT(OUT):: VAL(LENG)
      INTEGER(4)              :: I
      REAL(8)                 :: GBAS(3,3)
      REAL(8)                 :: SVAR
!     **************************************************************************
      IF(ID.EQ.'K') THEN
        IF(3*THIS%NKP.NE.LENG) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BRILLOUIN$GETR8A')
        END IF
        CALL GBASS(THIS%RBAS,GBAS,SVAR)
        DO I=1,THIS%NKP
          VAL(1+3*(I-1):3*I)=MATMUL(GBAS,THIS%XK(:,I))
        ENDDO
!
      ELSE IF(ID.EQ.'XK') THEN
        IF(3*THIS%NKP.NE.LENG) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BRILLOUIN$GETR8A')
        END IF
        VAL(:)=RESHAPE(THIS%XK,(/LENG/))
!
      ELSE IF(ID.EQ.'WKPT') THEN
        IF(THIS%NKP.NE.LENG) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BRILLOUIN$GETR8A')
        END IF
!       == CALCULATES WEIGHTS FOR A COMPLETELY FILLED BAND. ====================
!       == ATTENTION DEPENDS ON ACTUAL RBAS ====================================
        CALL BRILLOUIN_SAMFAC(1,THIS%NKP,(/(0.D0,I=1,THIS%NKP)/),1.D0,VAL,THIS)
        IF(ABS(SUM(VAL)-1.D0).GT.1.D-9) THEN
          CALL ERROR$MSG('K-POINT WEIGHTS DO NOT SUM UP TO ONE')
          CALL ERROR$R8VAL('TOTAL WEIGHT',SUM(VAL))
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BRILLOUIN$GETR8A')
        END IF
      ELSE
        CALL ERROR$MSG('UNKNOWN ID')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BRILLOUIN$GETR8A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$GETI4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(ID.EQ.'NK') THEN
        VAL=THIS%NKP
      ELSE
        CALL ERROR$MSG('UNKNOWN ID')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BRILLOUIN$GETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$GETI4A(ID,LEN,VAL)
!     **************************************************************************
!     **************************************************************************
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(OUT):: VAL(LEN)
!     **************************************************************************
      IF(ID.EQ.'NKDIV') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BRILLOUIN$GETI4A')
        END IF
        VAL=THIS%NKDIV
      ELSE
        CALL ERROR$MSG('UNKNOWN ID')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BRILLOUIN$GETI4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$IKP(XK,IKPT)
!     **************************************************************************
!     **  TRANSLATES A K-POINT IN RELATIVE COORDINATES INTO AN INDEX          **
!     **  OF THE CORRESPONDING IRREDUCIBLE K-POINT                            **
!     **                                                                      **
!     **************************************************************************
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: XK(3)
      INTEGER(4),INTENT(OUT) :: IKPT
      REAL(8)                :: XK1(3)
      INTEGER(4)             :: IP(3)
      INTEGER(4)             :: I
!     **************************************************************************
      XK1(:)=MODULO(XK(:),1.D0)
      XK1(:)=XK1(:)*REAL(THIS%NKDIV(:),KIND=8)-0.5D0*REAL(THIS%ISHIFT(:),KIND=8)
      IP(:)=NINT(XK1(:))
      IF(SUM((XK1(:)-REAL(IP(:),KIND=8))**2).GT.1.D-5) THEN
        CALL ERROR$MSG('INVALID VALUE FOR K-POINT POSITION IN REL. COORD.')
        CALL ERROR$STOP('BRILLOUIN$IKP')
      END IF
!     __THIS COORDINATE CONTAINS ALL BOUNDARIES OF THE BOX!!___________________
      I=1+IP(3)+(THIS%NKDIV(3)+1)*(IP(2)+(THIS%NKDIV(2)+1)*IP(1))  ! COORDINATE
      IKPT=THIS%IRRKP(I)
      RETURN
      END
!
!===============================================================================
!===============================================================================
!====                                                                       ====
!====  BLOCK ARBMSH:                                                        ====
!====                                                                       ====
!====  READS SYMMETRYELEMENTS AND FINDS IRR. K-POINTS                       ====
!====  AND TETRAHEDRA, USED IN THE BRILLOUIN ZONE INTEGRATION               ====
!====                                                                       ====
!===============================================================================
!===============================================================================
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$GRID(RBAS,NSYM,IIO,NKTARGET,TSHIFT,NKDIV,ISHIFT)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      INTEGER(4),INTENT(IN) :: NSYM       ! #(SYMMETRY OPERATIONS)
      INTEGER(4),INTENT(IN) :: IIO(3,3,NSYM) !SYMMETRY OPERATIONS
      LOGICAL(4),INTENT(IN) :: TSHIFT   ! ATTEMPT GRID SHIFT
      INTEGER(4),INTENT(IN) :: NKTARGET
      INTEGER(4),INTENT(OUT):: NKDIV(3)
      INTEGER(4),INTENT(OUT):: ISHIFT(3)
      REAL(8)               :: GBAS(3,3)
      REAL(8)               :: SVAR
      INTEGER(4)            :: IARB(3)
      INTEGER(4)            :: I,NMSHP
      LOGICAL(4)            :: TCHK
!     **************************************************************************
!
!     ==========================================================================
!     == CALCULATE IARB                                                       ==
!     ==========================================================================
      IARB=(/0,0,0/)
      DO I=1,NSYM
        IF(IIO(2,3,I).NE.0.OR.IIO(3,2,I).NE.0) IARB(1)=1
        IF(IIO(1,3,I).NE.0.OR.IIO(3,1,I).NE.0) IARB(2)=1
        IF(IIO(1,2,I).NE.0.OR.IIO(2,1,I).NE.0) IARB(3)=1
      ENDDO       
!
!     ==========================================================================
!     == DEFINE GRID                                                          ==
!     ==========================================================================
      CALL GBASS(RBAS,GBAS,SVAR)                                             
      NMSHP=NKTARGET
      CALL BRILLOUIN_BASDIV(NKDIV,NMSHP,GBAS,IARB)
!
!     ==========================================================================
!     == ATTEMPT TO SHIFT GRID                                                ==
!     ==========================================================================
      IF(TSHIFT) THEN
        ISHIFT(:)=1
        CALL BRILLOUIN_CHECKSHIFT(ISHIFT,NSYM,IIO,TCHK)
        IF(.NOT.TCHK) ISHIFT(:)=0
      ELSE
        ISHIFT(:)=0
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$MSH(RBAS,NGKP,NSYM,IIO,IARB,TSHIFT)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATE IRREDUCIBLE K-POINTS                                      **
!     **  AND FINDS INEQUIVALENT TETRAHEDRA                                   **
!     **  FOR BRILLOUIN ZONE INTEGRATION                                      **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **    1) IARB: DEPENDENCIES FOR DIVISIONS OF RECIPROCAL                 **
!     **       LATTICE VECTORS.                                               **
!     **                ( IF IARB(1)=1 THEN 1ST AND 2ND LATTICE               **
!     **                  VECTORS ARE DIVIDED BY AN EQUAL NUMBER;             **
!     **                  IF IARB(2)=1 THEN SAME FOR 2ND AND 3RD;             **
!     **                  IF IARB(3)=1 THEN SAME FOR 3RD AND 1ST)             **
!     **                                                                      **
!     **   AUTHOR: PETER E. BLOECHL                                           **
!     **                                                                      **
!     **   SUBROUTINES USED:                                                  **
!     **     GBASS,BASDIV,REDUZ,ZUORD,TETDIV,TETCNT,ORD1                      **
!     **                                                                      **
!     **************************************************************************
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)  ! LATTICE VECTORS (REAL SPACE)
      INTEGER(4),INTENT(IN) :: NGKP       ! TARGET FOR #(GENERAL K-POINTS)
      INTEGER(4),INTENT(IN) :: NSYM       ! #(SYMMETRY OPERATIONS)
      INTEGER(4),INTENT(IN) :: IIO(3,3,NSYM) !SYMMETRY OPERATIONS
      INTEGER(4),INTENT(IN) :: IARB(3)    ! DEPENDENCE
      LOGICAL(4),INTENT(INOUT) :: TSHIFT  ! ATTEMPTS A SHIFT FROM GAMMA POINT
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      REAL(8)   ,ALLOCATABLE:: XK(:,:)   !(3,NGKP) K-POINTS
      INTEGER(4),PARAMETER  :: NOPX=48    ! MAX #(POINT GROUP OPERATIONS)
      INTEGER(4)            :: NOP        ! #(POINT GROUP OPERATIONS)
      INTEGER(4)            :: OP(3,3,NOPX)  ! POINT GROUP OPERATIONS
      REAL(8)               :: GBAS(3,3)
      INTEGER(4)            :: N(3)       !DIVISION OF REC.LATT. VECT.
      INTEGER(4)            :: IWBK=0
      INTEGER(4)            :: IWNUM=0
      INTEGER(4)            :: NMSHP
      INTEGER(4),ALLOCATABLE:: NUM(:)   ! POINTS FROM GENERAL TO IRR. K-POINT
      INTEGER(4)            :: ISHIF(3)
      INTEGER(4)            :: ISYM,I,J,K
      INTEGER(4)            :: INV
      INTEGER(4)            :: NKP
      INTEGER(4)            :: INDEX1,INDEX2,NDIMM,NDI1,NDI2
      INTEGER(4)            :: TET0(3,4,6)
      REAL(8)               :: DUMMY
      LOGICAL(4)            :: TCHK
!     **************************************************************************
                                   CALL TRACE$PUSH('BRILLOUIN$MSH') 
      CALL ERROR$MSG('OBSOLETE ROUTINE')
      CALL ERROR$MSG('USE BRILLOUIN$GRID AND BRILLOUIN$MSHSYM INSTEAD')
      CALL ERROR$STOP('BRILLOUIN$MSH')
      THIS%RBAS=RBAS
      INV=0
!     ==========================================================================
!     == DEFINE MESH                                                          ==
!     ==========================================================================
      NMSHP=NGKP                                                         
      CALL GBASS(RBAS,GBAS,DUMMY)                                             
      CALL BRILLOUIN_BASDIV(N,NMSHP,GBAS,IARB)
      THIS%NKDIV(:)=N(:)
!
!     ==========================================================================
!     == EXPAND GENERATORS TO COMPLETE FULL POINT GROUP                       ==
!     ==========================================================================
      NOP=NSYM
      OP(:,:,1:NSYM)=IIO(:,:,1:NSYM)
      CALL BRILLOUIN_COMPLETE(NOPX,NOP,OP)
!
!     ==========================================================================
!     == ATTEMPT TO SHIFT GRID                                                ==
!     ==========================================================================
      IF(TSHIFT) THEN
        ISHIF(:)=1
        CALL BRILLOUIN_CHECKSHIFT(ISHIF,NOP,OP,TCHK)
        IF(.NOT.TCHK) ISHIF(:)=0
        TSHIFT=TCHK  ! RETURN FEEDBACK WHETHER SHIFT WAS SUCCESSFUL
      ELSE
        ISHIF(:)=0
      END IF
      THIS%ISHIFT(:)=ISHIF
!     ==========================================================================
!     == FIND IRREDUCIBLE K-POINTS                                            ==
!     ==========================================================================
      ALLOCATE(NUM(NMSHP))
      CALL BRILLOUIN_REDUZ(N,NMSHP,ISHIF,NOP,OP,NUM)
      IF(ISHIF(1).NE.0.OR.ISHIF(2).NE.0.OR.ISHIF(3).NE.0)INV=0   
      ALLOCATE(XK(3,NGKP))
      CALL BRILLOUIN_ZUORD(NMSHP,NUM,N,ISHIF,NGKP,NKP,XK)            
      THIS%NKP=NKP
      ALLOCATE(THIS%XK(3,NKP))
      THIS%XK(:,:)=XK(:,1:NKP)
      DEALLOCATE(XK)
!
!     ==========================================================================
!     == CONSTRUCT MAPPING ONTO IRREDUCIBLE K-POINTS                          ==
!     ==========================================================================
!     __ THE GENERAL POSITION IS ENCODED ACCORDING TO __________________________
!     __IND=1+IP(3)+(THIS%NKDIV(3)+1)*(IP(2)+(THIS%NKDIV(2)+1)*IP(1))___________
!     __ THIS MAPPING ALSO INCLUDES ALL BOUNDARIES OF THE UNIT CELL! ___________
      ALLOCATE(THIS%IRRKP(NMSHP))
      THIS%IRRKP(:)=NUM(:)
!
!     ==========================================================================
!     -- PRINTOUT OF SYMMETRY OPERATIONS                                      ==
!     ==========================================================================
      IF(TPR) THEN
        DO ISYM=1,NOP                                             
          WRITE(*,FMT='("SYMMETRYMATRIX NR. : ",I5/3(" ",3I10/))') &
     &            ISYM,((OP(I,J,ISYM),J=1,3),I=1,3)                 
        ENDDO
      END IF                                                            
!     ==========================================================================
!     ==  PRINTOUT OF MAPPING TO IRREDUCIBLE K-POINTS                         ==
!     ==========================================================================
      IF(IWNUM.EQ.1.AND.N(3)+1.LE.25) THEN                              
        DO I=0,N(1)                                                  
          DO J=N(2),0,-1                                               
            INDEX1=I*(N(2)+1)*(N(3)+1)+J*(N(3)+1)+1                         
            INDEX2=INDEX1+N(3)                                              
            WRITE(*,FMT='(25I5)')(NUM(K),K=INDEX1,INDEX2)                      
          ENDDO
          WRITE(*,FMT='(/)')
        ENDDO
      END IF                                                            
!     ==========================================================================
!     ==  PRINT IREDUCIBLE K-POINTS                                           ==
!     ==========================================================================
      WRITE(*,FMT='("NO. OF INEQUIVALENT K-POINTS ",I5)')NKP
      IF(IWBK.NE.0) THEN                                                
        WRITE(*,FMT='("INEQUIVALENT BLOCH VECTORS")')                       
        NDIMM=NKP/2+1                                                   
        DO I=1,NDIMM                                                 
          NDI1=(I-1)*2+1                                                  
          NDI2=NDI1+1                                                     
          IF(I.EQ.NDIMM)NDI2=NKP                                          
          WRITE(6,FMT='(I4,"(",3F10.6,")",I4," (",3F10.6,")")') &
     &          (K,(XK(J,K),J=1,3),K=NDI1,NDI2)                    
        ENDDO
      END IF                                                            
!     ==========================================================================
!     ==  CHOOSE TETRAHEDRA                                                   ==
!     ==========================================================================
      CALL BRILLOUIN_TETDIV(N,GBAS,TET0)                                      
      IF(TPR) THEN                                                 
        WRITE(*,FMT='(3(I2,"-TES NORMTETRAHEDRON "))')1,2,3
        WRITE(*,FMT='(4I5,5X,4I5,5X,4I5)') &
     &         (((TET0(J,K,I),K=1,4),I=1,3),J=1,3)    
        WRITE(*,FMT='(3(I2,"-TES NORMTETRAHEDRON "))')4,5,6
        WRITE(*,FMT='(4I5,5X,4I5,5X,4I5)') &
     &         (((TET0(J,K,I),K=1,4),I=4,6),J=1,3)    
      END IF                                                            
!     ==========================================================================
!     ==  FIND INEQUIVALENT TETRAHEDRA                                        ==
!     ==========================================================================
      CALL BRILLOUIN_TETCNT(NMSHP,NUM,TET0,N,INV,THIS)
      DEALLOCATE(NUM)
                                   CALL TRACE$POP
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$MSHSYM(RBAS,NKDIV,ISHIFT,TINV,NSYM,IIO)
!     **************************************************************************
!     **  CALCULATE IRREDUCIBLE K-POINTS AND INEQUIVALENT TETRAHEDRA          **
!     **  FOR BRILLOUIN ZONE INTEGRATION                                      **
!     **                                                                      **
!     **  IIO ARE THE GENERATORS OF THE SPATIAL SYMMETRY OPERATIONS.          **
!     **  ADD TIME INVERSION SYMMETRY USING TINV=.TRUE.                       **
!     **                                                                      **
!     **   SUBROUTINES USED:                                                  **
!     **     GBASS,REDUZ,ZUORD,TETDIV,TETCNT,ORD1                             **
!     **                                                                      **
!     ********************************PETER E. BLOECHL, STUTTGART, GOSLAR*******
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)  ! LATTICE VECTORS (REAL SPACE)
      INTEGER(4),INTENT(IN) :: NKDIV(3)   ! DIV. OF THE REC. LATTICE VECTORS
      INTEGER(4),INTENT(IN) :: ISHIFT(3)  ! GRID DOSPLACEMENT FROM GAMMA
      LOGICAL(4),INTENT(IN) :: TINV       ! TIME INVERSION SYMMETRY?
      INTEGER(4),INTENT(IN) :: NSYM       ! #(SYMMETRY OPERATIONS (GENERATORES))
      INTEGER(4),INTENT(IN) :: IIO(3,3,NSYM) !GENERATORS OF THE SYMMETRY GROUP
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      REAL(8)   ,ALLOCATABLE:: XK(:,:)    ! (3,NKP) K-POINTS
      INTEGER(4),PARAMETER  :: NOPX=48    ! MAX #(POINT GROUP OPERATIONS)
      INTEGER(4)            :: NOP        ! #(POINT GROUP OPERATIONS)
      INTEGER(4)            :: OP(3,3,NOPX)  ! POINT GROUP OPERATIONS
      REAL(8)               :: GBAS(3,3)
      INTEGER(4)            :: IWBK=0
      INTEGER(4)            :: IWNUM=0
      INTEGER(4)            :: NMSHP    ! #(KPOINTS IN CELL+BOUNDARIES)
      INTEGER(4),ALLOCATABLE:: NUM(:)   ! POINTS FROM GENERAL TO IRR. K-POINT
      INTEGER(4)            :: ISYM,I,J,K
      INTEGER(4)            :: NKP      !#(IRRED. K-POINTS)
      INTEGER(4)            :: INDEX1,INDEX2,NDIMM,NDI1,NDI2
      INTEGER(4)            :: TET0(3,4,6)
      INTEGER(4)            :: INV      ! ?? IS ALWAYS 0. IS USED IN TETCNT
      REAL(8)               :: DUMMY
      LOGICAL(4)            :: TCHK
!     **************************************************************************
                                   CALL TRACE$PUSH('BRILLOUIN$MSH') 
      THIS%RBAS=RBAS
      CALL GBASS(RBAS,GBAS,DUMMY)                                             
      NMSHP=(NKDIV(1)+1)*(NKDIV(2)+1)*(NKDIV(3)+1)
      THIS%NKDIV(:)=NKDIV(:)
      THIS%ISHIFT(:)=ISHIFT
!
!     ==========================================================================
!     == EXPAND GENERATORS TO COMPLETE FULL POINT GROUP                       ==
!     ==========================================================================
      NOP=NSYM
      OP(:,:,1:NSYM)=IIO(:,:,:)
      IF(TINV) THEN
!       == ADD TIME INVERSION UNLESS THERE IS ALREADY A SPATIAL INVERSION.
        TCHK=.FALSE.  ! SPATIAL INVERSION PRESENT?
        DO I=1,NSYM
          IF(IIO(1,1,I).NE.-1) CYCLE
          IF(IIO(2,2,I).NE.-1) CYCLE
          IF(IIO(3,3,I).NE.-1) CYCLE
          IF(IIO(1,2,I).NE.0) CYCLE
          IF(IIO(2,3,I).NE.0) CYCLE
          IF(IIO(3,1,I).NE.0) CYCLE
          IF(IIO(1,3,I).NE.0) CYCLE
          IF(IIO(3,2,I).NE.0) CYCLE
          IF(IIO(2,1,I).NE.0) CYCLE
          TCHK=.TRUE. ! SPATIAL INVERSION FOUND!
          EXIT
        ENDDO
        IF(.NOT.TCHK) THEN
          NOP=NOP+1
          OP(:,1,NOP)=(/-1,0,0/)
          OP(:,2,NOP)=(/0,-1,0/)
          OP(:,3,NOP)=(/0,0,-1/)
        END IF
      ENDIF
      CALL BRILLOUIN_COMPLETE(NOPX,NOP,OP)
!
!     ==========================================================================
!     == FIND IRREDUCIBLE K-POINTS                                            ==
!     ==========================================================================
      ALLOCATE(NUM(NMSHP))
      CALL BRILLOUIN_REDUZ(NKDIV,NMSHP,ISHIFT,NOP,OP,NUM)
      ALLOCATE(XK(3,NMSHP))
      CALL BRILLOUIN_ZUORD(NMSHP,NUM,NKDIV,ISHIFT,NMSHP,NKP,XK)            
      THIS%NKP=NKP
      ALLOCATE(THIS%XK(3,NKP))
      THIS%XK(:,:)=XK(:,1:NKP)
      DEALLOCATE(XK)
!
!     ==========================================================================
!     == CONSTRUCT MAPPING ONTO IRREDUCIBLE K-POINTS                          ==
!     ==========================================================================
!     __ THE GENERAL POSITION IS ENCODED ACCORDING TO __________________________
!     __IND=1+IP(3)+(THIS%NKDIV(3)+1)*(IP(2)+(THIS%NKDIV(2)+1)*IP(1))___________
!     __ THIS MAPPING ALSO INCLUDES ALL BOUNDARIES OF THE UNIT CELL! ___________
      ALLOCATE(THIS%IRRKP(NMSHP))
      THIS%IRRKP(:)=NUM(:)
!
!     ==========================================================================
!     -- PRINTOUT OF SYMMETRY OPERATIONS                                      ==
!     ==========================================================================
      IF(TPR) THEN
        DO ISYM=1,NOP                                             
          WRITE(*,FMT='("SYMMETRYMATRIX NR. : ",I5/3(" ",3I10/))') &
     &            ISYM,((OP(I,J,ISYM),J=1,3),I=1,3)                 
        ENDDO
      END IF                                                            
!
!     ==========================================================================
!     ==  PRINTOUT OF MAPPING TO IRREDUCIBLE K-POINTS                         ==
!     ==========================================================================
      IF(IWNUM.EQ.1.AND.NKDIV(3)+1.LE.25) THEN                              
        DO I=0,NKDIV(1)
          DO J=NKDIV(2),0,-1
            INDEX1=I*(NKDIV(2)+1)*(NKDIV(3)+1)+J*(NKDIV(3)+1)+1
            INDEX2=INDEX1+NKDIV(3)                                              
            WRITE(*,FMT='(25I5)')(NUM(K),K=INDEX1,INDEX2)                      
          ENDDO
          WRITE(*,FMT='(/)')
        ENDDO
      END IF                                                            
!
!     ==========================================================================
!     ==  PRINT IREDUCIBLE K-POINTS                                           ==
!     ==========================================================================
      WRITE(*,FMT='("NO. OF INEQUIVALENT K-POINTS ",I5)')NKP
      IF(IWBK.NE.0) THEN                                                
        WRITE(*,FMT='("INEQUIVALENT BLOCH VECTORS")')                       
        NDIMM=NKP/2+1                                                   
        DO I=1,NDIMM                                                 
          NDI1=(I-1)*2+1                                                  
          NDI2=NDI1+1                                                     
          IF(I.EQ.NDIMM)NDI2=NKP                                          
          WRITE(6,FMT='(I4,"(",3F10.6,")",I4," (",3F10.6,")")') &
     &          (K,(XK(J,K),J=1,3),K=NDI1,NDI2)                    
        ENDDO
      END IF                                                            
!
!     ==========================================================================
!     ==  CHOOSE TETRAHEDRA                                                   ==
!     ==========================================================================
      CALL BRILLOUIN_TETDIV(NKDIV,GBAS,TET0)
      IF(TPR) THEN                                                 
        WRITE(*,FMT='(3(I2,"-TES NORMTETRAHEDRON "))')1,2,3
        WRITE(*,FMT='(4I5,5X,4I5,5X,4I5)') &
     &         (((TET0(J,K,I),K=1,4),I=1,3),J=1,3)    
        WRITE(*,FMT='(3(I2,"-TES NORMTETRAHEDRON "))')4,5,6
        WRITE(*,FMT='(4I5,5X,4I5,5X,4I5)') &
     &         (((TET0(J,K,I),K=1,4),I=4,6),J=1,3)    
      END IF                                                            
!
!     ==========================================================================
!     ==  FIND INEQUIVALENT TETRAHEDRA                                        ==
!     ==========================================================================
      INV=0
      CALL BRILLOUIN_TETCNT(NMSHP,NUM,TET0,NKDIV,INV,THIS)
      DEALLOCATE(NUM)
                                   CALL TRACE$POP
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$MSHNOSYM(TINV,RBAS,NKDIV,ISHIFT)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATE IRREDUCIBLE K-POINTS                                      **
!     **  AND FINDS INEQUIVALENT TETRAHEDRA                                   **
!     **  FOR BRILLOUIN ZONE INTEGRATION                                      **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **    1) IARB: DEPENDENCIES FOR DIVISIONS OF RECIPROCAL                 **
!     **       LATTICE VECTORS.                                               **
!     **                ( IF IARB(1)=1 THEN 1ST AND 2ND LATTICE               **
!     **                  VECTORS ARE DIVIDED BY AN EQUAL NUMBER;             **
!     **                  IF IARB(2)=1 THEN SAME FOR 2ND AND 3RD;             **
!     **                  IF IARB(3)=1 THEN SAME FOR 3RD AND 1ST)             **
!     **                                                                      **
!     **   AUTHOR: PETER E. BLOECHL                                           **
!     **                                                                      **
!     **   SUBROUTINES USED:                                                  **
!     **     GBASS,BASDIV,REDUZ,ZUORD,TETDIV,TETCNT,ORD1                      **
!     **                                                                      **
!     **************************************************************************
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TINV       ! FLAG FOR TIME INVERSION SYMMETRY
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)  ! LATTICE VECTORS (REAL SPACE)
      INTEGER(4),INTENT(IN) :: NKDIV(3)   ! DIVISIONS OF RECIPR. LATTICE VECTORS
      INTEGER(4),INTENT(IN) :: ISHIFT(3)  ! GRID IS SHIFTED 
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      REAL(8)   ,ALLOCATABLE:: XK(:,:)   !(3,NGKP) K-POINTS
      REAL(8)               :: GBAS(3,3)
      INTEGER(4)            :: IWBK=0
      INTEGER(4)            :: IWNUM=0
      INTEGER(4)            :: NGKP
      INTEGER(4)            :: NMSHP
      INTEGER(4),ALLOCATABLE:: NUM(:)   ! POINTS FROM GENERAL TO IRR. K-POINT
      INTEGER(4)            :: ISYM,I,J,K
      INTEGER(4)            :: INV
      INTEGER(4)            :: NKP
      INTEGER(4)            :: INDEX1,INDEX2,NDIMM,NDI1,NDI2
      INTEGER(4)            :: TET0(3,4,6)
      REAL(8)               :: DUMMY
      INTEGER(4)            :: IIO(3,3,2)
      INTEGER(4)            :: NSYM  !#(SYMMETRY OPERATIONS) 1 OR 2
!     **************************************************************************
                                   CALL TRACE$PUSH('BRILLOUIN$MSHNOSYM')
      THIS%RBAS=RBAS
!
      INV=0
      IF(TINV) INV=1
      IIO(:,:,:)=0
!     == IDENTITY
      IIO(1,1,1)=1
      IIO(2,2,1)=1
      IIO(3,3,1)=1
!     == INVERSION
      IIO(:,:,2)=-IIO(:,:,1)
      NSYM=1
      IF(TINV) NSYM=2
      NGKP=NKDIV(1)*NKDIV(2)*NKDIV(3)
      THIS%ISHIFT(:)=ISHIFT(:)
      THIS%NKDIV(:)=NKDIV(:)
!     ==========================================================================
!     == DEFINE MESH                                                          ==
!     ==========================================================================
      NMSHP=(NKDIV(1)+1)*(NKDIV(2)+1)*(NKDIV(3)+1)
      CALL GBASS(RBAS,GBAS,DUMMY)                                             
!     ==========================================================================
!     == FIND IRREDUCIBLE K-POINTS                                            ==
!     ==========================================================================
      ALLOCATE(NUM(NMSHP))
      CALL BRILLOUIN_REDUZ(NKDIV,NMSHP,ISHIFT,NSYM,IIO,NUM)
!     == INV=1 EXPLOITS TIME INVERSION SYMMETRY FOR UNSHIFTED GRIDS  ==
!     == WHEN SEARCHING IRREDUCIBLE TETRAHEDRA. PROBABLY IRRELEVANT  ==
      IF(ISHIFT(1).NE.0.OR.ISHIFT(2).NE.0.OR.ISHIFT(3).NE.0)INV=0   
      IF(.NOT.TINV) INV=0
!
!     ==========================================================================
!     ==  DETERMINE IRREDUCIBLE K-POINTS                                      ==
!     ==========================================================================
      ALLOCATE(XK(3,NGKP))
      CALL BRILLOUIN_ZUORD(NMSHP,NUM,NKDIV,ISHIFT,NGKP,NKP,XK)            
      THIS%NKP=NKP
      ALLOCATE(THIS%XK(3,NKP))
      THIS%XK(:,:)=XK(:,1:NKP)
      DEALLOCATE(XK)
!
!     ==========================================================================
!     == CONSTRUCT MAPPING ONTO IRREDUCIBLE K-POINTS                          ==
!     ==========================================================================
!     __ THE GENERAL POSITION IS ENCODED ACCORDING TO __________________________
!     __IND=1+IP(3)+(THIS%NKDIV(3)+1)*(IP(2)+(THIS%NKDIV(2)+1)*IP(1))___________
!     __ THIS MAPPING ALSO INCLUDES ALL BOUNDARIES OF THE UNIT CELL! ___________
      ALLOCATE(THIS%IRRKP(NMSHP))
      THIS%IRRKP(:)=NUM(:)
!
!     ==========================================================================
!     == PRINTOUT OF SYMMETRY OPERATIONS                                      ==
!     ==========================================================================
      IF(TPR) THEN
        DO ISYM=1,NSYM                                             
          WRITE(*,FMT='("SYMMETRYMATRIX NR. : ",I5/3(" ",3I10/))') &
     &            ISYM,((IIO(I,J,ISYM),J=1,3),I=1,3)                 
        ENDDO
      END IF                                                            
!     ==========================================================================
!     ==  PRINTOUT OF MAPPING TO IRREDUCIBLE K-POINTS                         ==
!     ==========================================================================
      IF(IWNUM.EQ.1.AND.NKDIV(3)+1.LE.25) THEN                              
        DO I=0,NKDIV(1)                                                  
          DO J=NKDIV(2),0,-1                                               
            INDEX1=I*(NKDIV(2)+1)*(NKDIV(3)+1)+J*(NKDIV(3)+1)+1
            INDEX2=INDEX1+NKDIV(3)                                              
            WRITE(*,FMT='(25I5)')(NUM(K),K=INDEX1,INDEX2)                      
          ENDDO
          WRITE(*,FMT='(/)')
        ENDDO
      END IF      
!     ==========================================================================
!     --  PRINT IREDUCIBLE K-POINTS                                  -- 
!     ==========================================================================
      WRITE(*,FMT='("NO. OF INEQUIVALENT K-POINTS ",I5)')NKP
      IF(IWBK.NE.0) THEN                                                
        WRITE(*,FMT='("INEQUIVALENT BLOCH VECTORS")')                       
        NDIMM=NKP/2+1                                                   
        DO I=1,NDIMM                                                 
          NDI1=(I-1)*2+1                                                  
          NDI2=NDI1+1                                                     
          IF(I.EQ.NDIMM)NDI2=NKP                                          
          WRITE(*,FMT='(I4,"(",3F10.6,")",I4," (",3F10.6,")")') &
     &          (K,(XK(J,K),J=1,3),K=NDI1,NDI2)                    
        ENDDO
      END IF                                                            
!
!     ==========================================================================
!     ==  CHOOSE TETRAHEDRA                                                   ==
!     ==========================================================================
      CALL BRILLOUIN_TETDIV(NKDIV,GBAS,TET0)
      IF(TPR) THEN                                                 
        WRITE(*,FMT='(3(I2,"-TES NORMTETRAHEDRON "))')1,2,3
        WRITE(*,FMT='(4I5,5X,4I5,5X,4I5)') &
     &         (((TET0(J,K,I),K=1,4),I=1,3),J=1,3)    
        WRITE(*,FMT='(3(I2,"-TES NORMTETRAHEDRON "))')4,5,6
        WRITE(*,FMT='(4I5,5X,4I5,5X,4I5)') &
     &         (((TET0(J,K,I),K=1,4),I=4,6),J=1,3)    
      END IF                                                            
!
!     ==========================================================================
!     ==  FIND INEQUIVALENT TETRAHEDRA                                        ==
!     ==========================================================================
      CALL BRILLOUIN_TETCNT(NMSHP,NUM,TET0,NKDIV,INV,THIS)
      DEALLOCATE(NUM)
                                   CALL TRACE$POP
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_BASDIV(N,NMSHP,GBAS,IARB)
!     **************************************************************************
!     **  BASDIV DETERMINES DIVISION OF BASEVECTORS OF REC. LATT.             **
!     **  SO THAT THE NUMBER OF MESHPOINTS IS JUST BELOW NMSHP AND            **
!     **  TAKES INTO ACCOUNT THE DEPENDENCY BY POINTSYMMETRY                  **
!     **  INPUT :                                                             **
!     **    GBAS        RECIPROCAL LATTICE VECTORS                            **
!     **    IARB        DEPENDENCIES FOR DIVISIONS                            **
!     **                OF RECIPROCAL LATTICE VECTORS                         **
!     **                ( IF IARB(1)=1 THEN 1ST AND 2ND LATTICE               **
!     **                  VECTORS ARE DIVIDED BY AN EQUAL NUMBER;             **
!     **                  IF IARB(2)=1 THEN SAME FOR 2ND AND 3RD;             **
!     **                  IF IARB(3)=1 THEN SAME FOR 3RD AND 1ST)             **
!     **    NMSHP       TOTAL NUMBER OF GRID POINTS                           **
!     **  OUTPUT:                                                             **
!     **    N           NUMBER OF DIVISIONS OF RECIPROCAL LATTICE             **
!     **                VECTORS FOR DEFINITION OF SUBLATTICE                  **
!     **    NMSHP       NUMBER OF SUBLATTICE POINTS IN A REC. UNIT CEL        L*
!     **                AND ON ALL ITS FACES                                  **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT)   :: N(3)
      INTEGER(4),INTENT(INOUT) :: NMSHP
      REAL(8)   ,INTENT(IN)    :: GBAS(3,3)
      INTEGER(4),INTENT(IN)    :: IARB(3)
      REAL(8)                  :: BETR(3)
      REAL(8)                  :: RN(3)
      REAL(8)   ,PARAMETER     :: OPAR=1.D-6
      INTEGER(4)               :: I
      REAL(8)                  :: SVAR
!     **************************************************************************
!     ==========================================================================
!     == FIND UPPER LIMIT FOR THE LENGTH OF SUBLATTICE VECTORS                ==
!     ==========================================================================
      DO I=1,3                                                       
        BETR(I)=SQRT(GBAS(I,1)**2+GBAS(I,2)**2+GBAS(I,3)**2)             
      ENDDO
      SVAR=(REAL(NMSHP,KIND=8)/(BETR(1)*BETR(2)*BETR(3)))**(1.D0/3.D0)         
      DO I=1,3                                                       
        RN(I)=BETR(I)*SVAR                                                
      ENDDO
!
!     ==========================================================================
!     == FIND DIVISIONS OF LATTICE VECTORS                                    ==
!     ==========================================================================
      IF(IARB(1).EQ.1.AND.IARB(2).EQ.1)THEN                             
        N(1)=INT((RN(1)*RN(2)*RN(3))**(1.D0/3.D0)+OPAR)                 
        N(2)=N(1)                                                       
        N(3)=N(1)                                                       
      ELSE IF(IARB(1).EQ.1) THEN                                        
        N(1)=INT(SQRT(RN(1)*RN(2))+OPAR)                               
        N(2)=N(1)                                                       
        N(3)=INT(RN(3))                                                 
      ELSE IF(IARB(2).EQ.1) THEN                                        
        N(1)=INT(SQRT(RN(1)*RN(3))+OPAR)                               
        N(2)=INT(RN(2))                                                 
        N(3)=N(1)                                                       
      ELSE IF (IARB(3).EQ.1) THEN                                       
        N(1)=INT(RN(1))                                                 
        N(2)=INT(SQRT(RN(2)*RN(3))+OPAR)                               
        N(3)=N(2)                                                       
      ELSE                                                              
        N(1)=INT(RN(1)+OPAR)                                            
        N(2)=INT(RN(2)+OPAR)                                            
        N(3)=INT(RN(3)+OPAR)                                            
      END IF                                                            
      N(1)=MAX(1,N(1))                                                 
      N(2)=MAX(1,N(2))                                                 
      N(3)=MAX(1,N(3))                                                 
!     ==========================================================================
!     == USE THIS TO FIX THE K-MESH PY HAND                                   ==
!     ==========================================================================
!     PRINT*,' K-MESH PUT IN BY HAND!!!!!!!!!!'                         
!     N(1)=12                                                           
!     N(2)=12                                                           
!     N(3)=2                                                            
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      NMSHP=(N(1)+1)*(N(2)+1)*(N(3)+1)                                  
      WRITE(*,FMT='("NO. OF MESH POINTS IN THE BRILLOUIN ZONE =",I6)') &
     &            N(1)*N(2)*N(3)      
      WRITE(*,FMT= &
     &        '("DIVISION OF RECIPROCAL LATTICE VECTORS (INTERVALS)=",3I5)')N(:)
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_REDUZ(N,NMSHP,ISHIFT,NSYM,IO,NUM)
!     **************************************************************************
!     **  REDUZ CREATES THE RELATION BETWEEN THE MESHPOINTS AND               **
!     **  THE POINTS IN THE "IRREDUCIBLE ZONE"                                **
!     **  INPUT :                                                             **
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS           **
!     **    NMSHP       NUMBER OF SUBLATTICE POINTS INSIDE AND                **
!     **                ON ALL FACES OF A REC. UNIT CELL                      **
!     **    NSYM        NUMBER OF SYMMETRY OPERATIONS                         **
!     **    IO          SYMMETRY OPERATIONS (BASIS ARE REC. LATT. VEC.)       **
!     **  OUTPUT :                                                            **
!     **    NUM(I)      MAPPING FROM A GENERAL POINT (I) TO THE               **
!     **                CORRESPONDING IRREDUCIBLE POINT (NUM)                 **
!     **  REMARKS :                                                           **
!     **    THE MAPPING FROM COORDINATES TO NUMBERS IS GIVEN BY :             **
!     **   (X,Y,Z)=RBAS*(I,J,K)                                               **
!     **   (I,J,K)  <->  NUM = I*(N(2)+1)*(N(3)+1)+J*(N(3)+1)+K+1             **
!     **                                                                      **
!     **  REMARK: A CHANGE FROM 5B51787 TO B86AE50 RESULTED IN A CHANGE OF THE**
!     **    ORDER OF IRREDUCIBLE K-POINTS, WHICH RENDERED PREVIOUS RESTART    **
!     **    FILES UNREADABLE                                                  **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N(3)
      INTEGER(4),INTENT(IN) :: NMSHP
      INTEGER(4),INTENT(IN) :: ISHIFT(3)
      INTEGER(4),INTENT(IN) :: NSYM
      INTEGER(4),INTENT(IN) :: IO(3,3,NSYM)
      INTEGER(4),INTENT(OUT):: NUM(NMSHP)
      LOGICAL(4),PARAMETER  :: TSHIFT=.TRUE.
      INTEGER(4)            :: ISYM(NMSHP)  ! SYMMETRY FROM IRR. P. OR ZERO
      INTEGER(4)            :: I,J,JPRIME
      INTEGER(4)            :: I1,I2,I3,J1,J2,J3
      INTEGER(4)            :: IX1,IX2,IX3,IND,IND1
      LOGICAL(4)            :: TCHK
!     **************************************************************************
                           CALL TRACE$PUSH('BRILLOUIN_REDUZ')
      CALL BRILLOUIN_CHECKSHIFT(ISHIFT,NSYM,IO,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('SHIFT INCOMPATIBLE WITH SYMMETRY')
        CALL ERROR$STOP('BRILLOUIN_REDUZ')
      END IF
!
!     ==========================================================================
!     ==  REDUCTION OF REDUCIBLE K-POINTS                                     ==
!     ==  APPLY THE FULL POINT GROUP ON EACH K-POINT. THE K-POINT WITH THE    ==
!     ==  MINIMUM INDEX IN THE RESULTING STAR IS THE IRRDUCIBLE POINT         ==
!     ==                                                                      ==
!     ==  POSITION IS ENCODED AS NUM(I1,I2,I3)=1+I3+(N(3)+1)*(I2+(N(2)+1)*I1) ==
!     ==  WHERE UM RUNS FROM 1 TO NMSHP                                       ==
!     ==========================================================================
!     == SCAN ALL POINTS INSIDE THE FIRST UNIT CELL (EXCLUDING OUTER BOUNDARIES)
!     == THE ORDER OF THE POINTS DETERMINES THE ORDER AND THE CHOICE OF 
!     == IRREDUCIBLE K-POINTS.
      NUM(:)=0
      DO IND1=1,N(1)*N(2)*N(3)
!       == INVERT IND1=1+I1+N(1)*(I2+N(2)*I3) ==================================
        IND=IND1-1
        I3=IND/(N(1)*N(2))
        IND=IND-I3*N(1)*N(2)
        I2=IND/N(1)
        I1=IND-I2*N(1)
!       == ENCODE POSITION AND CHECK IF POINT IS ALREADY REDUCIBLE ============
        J=1+I3+(N(3)+1)*(I2+(N(2)+1)*I1)  ! ENCODES POSITION
        IF(NUM(J).NE.0) CYCLE  ! NO IRREDUCIBLE POINT
!
!       == PLACE ON A GRID THAT IS TWICE AS FINE TO ACCOUNT FOR THE SHIFT=======
        I1=2*I1+ISHIFT(1)                                            
        I2=2*I2+ISHIFT(2)                                            
        I3=2*I3+ISHIFT(3)                                            
!
!       == CONSTRUCT THE STAR OF EQUIVALENT K-POINTS ===========================
        DO I=1,NSYM                                                    
!         == APPLY SYMMETRY ====================================================
          J1=IO(1,1,I)*I1+IO(1,2,I)*I2+IO(1,3,I)*I3             
          J2=IO(2,1,I)*I1+IO(2,2,I)*I2+IO(2,3,I)*I3             
          J3=IO(3,1,I)*I1+IO(3,2,I)*I2+IO(3,3,I)*I3             
!         == CHECK IF THE POINT FALLS ONTO THE ORIGINAL SHIFTED GRID ===========
          IX1=MODULO(J1-ISHIFT(1),2)                                           
          IX2=MODULO(J2-ISHIFT(2),2)                                           
          IX3=MODULO(J3-ISHIFT(3),2)                                           
          IF(IX1.NE.0.OR.IX2.NE.0.OR.IX3.NE.0) CYCLE  
!         == COARSEN GRID AGAIN ================================================
          J1=(J1-ISHIFT(1))/2                                               
          J2=(J2-ISHIFT(2))/2                                               
          J3=(J3-ISHIFT(3))/2                                               
!         == TRANSLATE TARGET POINT INTO FIRST UNIT CELL =======================
          J1=MODULO(J1,N(1))
          J2=MODULO(J2,N(2))
          J3=MODULO(J3,N(3))
!         == NEW POSITION ======================================================
          JPRIME=1+J3+(N(3)+1)*(J2+(N(2)+1)*J1)
!         == POINT TO IRREDUCIBLE POINT AND SPECIFY SYMMETRY  ==================
          NUM(JPRIME)=J   !POINTER TOWARDS IRREDICIBLE POINT
          ISYM(JPRIME)=I  !SYMMETRY OPERATION FROM IRREDUCEBLE POINT TO JPRIME
        ENDDO
!       == ENFORCE THAT THE IDENTITY IS USED TO MAP THE POINT ON ITSELF ========
        NUM(J)=J
        ISYM(J)=0
      ENDDO
!
!     ==========================================================================
!     == COMNPLETE OUTER BOUNDARIES                                           ==
!     ==========================================================================
      DO I1=0,N(1)-1
        DO I2=0,N(2)-1
          I3=0
          J3=N(3)
          J     =1+I3+(N(3)+1)*(I2+(N(2)+1)*I1)
          JPRIME=1+J3+(N(3)+1)*(I2+(N(2)+1)*I1)
          NUM(JPRIME)=NUM(J)
          ISYM(JPRIME)=ISYM(J)
        ENDDO
      ENDDO
      DO I1=0,N(1)-1
        DO I3=0,N(3)
          I2=0
          J2=N(2)
          J     =1+I3+(N(3)+1)*(I2+(N(2)+1)*I1)
          JPRIME=1+I3+(N(3)+1)*(J2+(N(2)+1)*I1)
          NUM(JPRIME)=NUM(J)
          ISYM(JPRIME)=ISYM(J)
        ENDDO
      ENDDO
      DO I2=0,N(2)
        DO I3=0,N(3)
          I1=0
          J1=N(1)
          J     =1+I3+(N(3)+1)*(I2+(N(2)+1)*I1)
          JPRIME=1+I3+(N(3)+1)*(I2+(N(2)+1)*J1)
          NUM(JPRIME)=NUM(J)
          ISYM(JPRIME)=ISYM(J)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == THE FOLLOWING IS NEEDED BECAUSE THE OPERATIONS MAY NOT COMMUTE?
!     == A STAR IS CHARACTERIZED BY HAVING THE SAME VALUE OF NUM(I)    
!     ==========================================================================
      DO I=1,NMSHP
        IF(NUM(I).NE.NUM(NUM(I))) THEN
          CALL ERROR$MSG('MAPPING FAILED')
          CALL ERROR$I4VAL('I',I)
          CALL ERROR$I4VAL('NUM(I)',NUM(I))
          CALL ERROR$I4VAL('NUM(NUM(I))',NUM(NUM(I)))
          CALL ERROR$STOP('BRILLOUIN_REDUZ')
        END IF
      ENDDO
                             CALL TRACE$POP()
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_CHECKSHIFT(ISHIFT,NSYM,IO,TCHK)
!     **************************************************************************
!     **                                                                      **
!     **  CHECKSHIFT CHECKS THE CONSISTENCY OF A K-GRID SHIFT WITH            **
!     **  THE SYMMETRY OPERATIONS                                             **
!     **  INPUT :                                                             **
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS           **
!     **    NSYM        NUMBER OF SYMMETRY OPERATIONS                         **
!     **    IO          SYMMETRY OPERATIONS (BASIS ARE REC. LATT. VEC.)       **
!     **  OUTPUT :                                                            **
!     **    TCHK        FLAG FOR CONSISTENCY OF THE SHIFT                     **
!     **  REMARKS :                                                           **
!     **    THE MAPPING FROM COORDINATES TO NUMBERS IS GIVEN BY :             **
!     **   (X,Y,Z)=RBAS*(I,J,K)                                               **
!     **   (I,J,K)  <->  NUM = I*(N(2)+1)*(N(3)+1)+J*(N(3)+1)+K+1             **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISHIFT(3)
      INTEGER(4),INTENT(IN) :: NSYM
      INTEGER(4),INTENT(IN) :: IO(3,3,NSYM)
      LOGICAL(4),INTENT(OUT):: TCHK
      INTEGER(4)            :: IN(3,8)      
      INTEGER(4)            :: I,J,K,L
      INTEGER(4)            :: I1,I2,I3,J1,J2,J3
!     **************************************************************************
!
!     ==================================================================
!     ==  TEST IF SUBMESH SHIFT IS COMPATIBLE WITH THE POINT GROUP    ==
!     ==================================================================
      DO I=1,8                                                       
        IN(1,I)=(I-1)/4                                                   
        IN(2,I)=(I-IN(1,I)*4-1)/2                                         
        IN(3,I)=I-IN(1,I)*4-IN(2,I)*2-1                                   
      ENDDO
      DO I=1,NSYM                                                    
        DO J=1,8                                                       
          I1=2*IN(1,J)+ISHIFT(1)                                            
          I2=2*IN(2,J)+ISHIFT(2)                                            
          I3=2*IN(3,J)+ISHIFT(3)                                            
          J1=IO(1,1,I)*I1+IO(1,2,I)*I2+IO(1,3,I)*I3                         
          J2=IO(2,1,I)*I1+IO(2,2,I)*I2+IO(2,3,I)*I3                         
          J3=IO(3,1,I)*I1+IO(3,2,I)*I2+IO(3,3,I)*I3                         
          IF(MOD(REAL(J1-ISHIFT(1),KIND=8),2.D0).NE.0.D0.OR. &
     &       MOD(REAL(J2-ISHIFT(2),KIND=8),2.D0).NE.0.D0.OR. &
     &       MOD(REAL(J3-ISHIFT(3),KIND=8),2.D0).NE.0.D0) THEN
             TCHK=.FALSE.
             PRINT*,'SUBMESH SHIFT CONFLICTS WITH POINTGROUP'                
             WRITE(6,FMT="('SHIFT=',3I5)")ISHIFT
             WRITE(6,FMT="('SYMMETRYMATRIX NR. : ',I5/3(' ',3I10/))") &
     &               I,((IO(K,L,I),L=1,3),K=1,3)                  
             RETURN 
          END IF                                                            
        ENDDO
      ENDDO
      TCHK=.TRUE.
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_ZUORD(NMSHP,NUM,N,ISHIFT,IDKP,NKP,XK)             
!     **                                                              **
!     **  ZUORD CREATES A RELATION BETWEEN "IRREDUCIBLE POINTS" AND   **
!     **  THEIR PLACES ON THE FILE, COUNTS THEM, AND CALCULATES THEIR **
!     **  COORDINATES IN K-SPACE                                      **
!     **                                                              **
!     **  INPUT:                                                      **
!     **    NMSHP       NUMBER OF SUBLATTICE POINTS INSIDE AND        **
!     **                ON ALL FACES OF A REC. UNIT CELL              **
!     **    NUM(I)      MAPPING FROM A GENERAL POINT (I) TO THE       **
!     **                CORRESPONDING IRREDUCIBLE POINT (NUM)         **
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    ISHIFT                                                    **
!     **    IDKP        MAXIMUM NUMBER OF IRREDUCIBLE K-POINTS        **
!     **    GBAS        RECIPROCAL LATTICE VECTORS                    **
!     **  OUTPUT :                                                    **
!     **    NKP         ACTUAL NUMBER OF IRREDUCIBLE K-POINTS         **
!     **    NUM(I)      MAPPING FROM A GENERAL POINT (I) TO THE       **
!     **                CORRESPONDING IRREDUCIBLE POINT (NUM)         **
!     **                (SEE REMARKS)                                 **
!     **    XK          IRREDUCIBLE K-POINT IN RELATIVE COORDINATES  **
!     **                                                              **
!     **  REMARKS :                                                   **
!     **    NUM(I) REFERS TO THE ACTUAL POSITION OF THE IRREDUCIBLE   **
!     **    K-POINT ON INPUT. ON OUTPUT IT ONLY REFERS TO THE         **
!     **    SECOND INDEX OF XK(I,IKP)                                 **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NMSHP
      INTEGER(4),INTENT(OUT):: NUM(NMSHP)
      INTEGER(4),INTENT(IN) :: N(3)
      INTEGER(4),INTENT(IN) :: ISHIFT(3)
      INTEGER(4),INTENT(IN) :: IDKP
      INTEGER(4),INTENT(OUT):: NKP
      REAL(8)   ,INTENT(OUT):: XK(3,IDKP)
      INTEGER(4)            :: NDIM
      INTEGER(4)            :: I1,I2,I3,I
      INTEGER(4)            :: MAP(NMSHP)
!     ******************************************************************
                           CALL TRACE$PUSH('BRILLOUIN_ZUORD')
      MAP(:)=0
      NDIM=0                                                            
      DO I1=1,N(1)                                                
        DO I2=1,N(2)
          DO I3=1,N(3)
!           == HAD TO CHANGE THE ORDER OF THE DO LOOPS BECAUSE OF     ==
!           == COMPATIBILITY WITH PREVIOUS IMPLEMENTATION IN CP-PAW.  ==
!           == THUS THE INDEX I BELOW DOES NOT INCREASE MONOTONICALLY!==
!            I=I1+(N(1)+1)*(I2-1+(N(2)+1)*(I3-1))  
            I=I3+(N(3)+1)*(I2-1+(N(2)+1)*(I1-1))  ! COORDINATE
            IF(I.GT.NMSHP) THEN                                               
              CALL ERROR$MSG('I.GT.NMSHP,STOP')
              CALL ERROR$I4VAL('I',I)
              CALL ERROR$I4VAL('I1',I1)
              CALL ERROR$I4VAL('I2',I2)
              CALL ERROR$I4VAL('I3',I3)
              CALL ERROR$STOP('BRILLOUIN_ZUORD')
            END IF                                                            
            IF(I.EQ.NUM(I))THEN                                               
              NDIM=NDIM+1                                                     
              IF(NDIM.GT.IDKP) THEN                                           
                CALL ERROR$MSG('NUMBER OF INEQUIVALENT POINTS ECCEEDS IDKP')
                CALL ERROR$STOP('BRILLOUIN_ZUORD')
              END IF                                                          
              MAP(I)=NDIM    
              XK(1,NDIM)=(REAL(I1-1,KIND=8)+REAL(ISHIFT(1),KIND=8)/2.D0) &
     &                                     /REAL(N(1),KIND=8)              
              XK(2,NDIM)=(REAL(I2-1,KIND=8)+REAL(ISHIFT(2),KIND=8)/2.D0) &
     &                                     /REAL(N(2),KIND=8)              
              XK(3,NDIM)=(REAL(I3-1,KIND=8)+REAL(ISHIFT(3),KIND=8)/2.D0) &
     &                                     /REAL(N(3),KIND=8)              
            END IF
          ENDDO
        ENDDO
      ENDDO
      NKP=NDIM                                                          
      DO I=1,NMSHP
        NUM(I)=MAP(NUM(I))
        IF(NUM(I).EQ.0) THEN
          CALL ERROR$MSG('MAPPING ERROR')      
          CALL ERROR$STOP('BRILLOUIN_ZUORD')
        END IF
      ENDDO
                             CALL TRACE$POP()
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_TETDIV(N,GBAS,TET0)
!     **************************************************************************
!     **  TETDIV DETERMINES THE DIVISION OF THE PARALLELEPIPEDS               **
!     **  IN TETRAHEDRONS ACCORDING TO THE SHORTEST DIAGONAL                  **
!     **                                                                      **
!     **  INPUT:                                                      **
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    GBAS        RECIPROCAL LATTICE VECTORS                    **
!     **  OUTPUT:                                                    ** 
!     **    TET0(I,J,K) COORDINATES (I) IN THE BASIS OF SUBLATTICE   ** 
!     **                VECTORS, OF THE 4 CORNERS (J) OF A           ** 
!     **                TETRAHEDRON (K) IN A SUBLATTICE UNIT CELL    ** 
!     **                                                             ** 
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N(3)
      REAL(8)   ,INTENT(IN) :: GBAS(3,3)
      INTEGER(4),INTENT(OUT):: TET0(3,4,6)
      REAL(8)               :: P(8,3)
      REAL(8)               :: DIAG(4)
      INTEGER(4)            :: IACHT(8)
      INTEGER(4)            :: TET(4,6)
      INTEGER(4)            :: I,J,K,L
      INTEGER(4)            :: MNDG
      INTEGER(4)            :: ISVAR
!     **************************************************************************
                           CALL TRACE$PUSH('BRILLOUIN_TETDIV')
!     ------------------------------------------------------------------
!     -- SEARCH FOR THE SHORTEST DIAGONAL                             --
!     ------------------------------------------------------------------
      DO I=0,1                                                       
        DO J=0,1                                                       
          DO K=0,1                                                       
            ISVAR=4*I+2*J+K+1                                                 
            DO L=1,3                                                       
              P(ISVAR,L)=GBAS(1,L)*REAL(I,KIND=8)/REAL(N(1),KIND=8) &
     &                  +GBAS(2,L)*REAL(J,KIND=8)/REAL(N(2),KIND=8) &
     &                  +GBAS(3,L)*REAL(K,KIND=8)/REAL(N(3),KIND=8)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!*                                                                      
      DO I=1,4                                                      
        DIAG(I)=0.D0                                                     
        DO J=1,3                                                      
          DIAG(I)=DIAG(I)+(P(I,J)-P(9-I,J))**2                             
        ENDDO
      ENDDO
!*                                                                     
      MNDG=1                                                           
      DO I=2,4                                                      
        IF(DIAG(I).LT.DIAG(MNDG)) THEN                                   
          MNDG=I                                                         
        END IF                                                           
      ENDDO
!     ------------------------------------------------------------------
!     -- ROTATE PARALLELEPIPED                                        --
!     ------------------------------------------------------------------
      IF(MNDG.EQ.1)THEN                                                 
        DO I=1,8                                                     
          IACHT(I)=I                                                      
        ENDDO
      ELSE IF(MNDG.EQ.2) THEN                                           
        DO I=1,4                                                     
          IACHT(2*I-1)=2*I                                                
          IACHT(2*I)=2*I-1                                                
        ENDDO
      ELSE IF(MNDG.EQ.3) THEN                                           
        DO I=0,1                                                     
          DO J=1,2                                                     
            IACHT(4*I+J)=4*I+J+2                                            
            IACHT(4*I+J+2)=4*I+J                                            
          ENDDO
        ENDDO
      ELSE IF(MNDG.EQ.4) THEN                                           
        DO I=1,4                                                     
          IACHT(I)=I+4                                                    
          IACHT(I+4)=I                                                    
        ENDDO
      END IF                                                            
!      **  CREATION OF TETRAHEDRA  **                                   
!      **  (1248);(1438);(1378);(1758);(1568);(1628)                    
       DO I=1,6                                                      
         TET(1,I)=IACHT(1)                                                
         TET(4,I)=IACHT(8)                                                
       ENDDO
       TET(2,1)=IACHT(2)                                                
       TET(3,1)=IACHT(4)                                                
       TET(2,2)=IACHT(4)                                                
       TET(3,2)=IACHT(3)                                                
       TET(2,3)=IACHT(3)                                                
       TET(3,3)=IACHT(7)                                                
       TET(2,4)=IACHT(7)                                                
       TET(3,4)=IACHT(5)                                                
       TET(2,5)=IACHT(5)                                                
       TET(3,5)=IACHT(6)                                                
       TET(2,6)=IACHT(6)                                                
       TET(3,6)=IACHT(2)                                                
!*                                                                      
      DO I=1,4                                                       
        DO J=1,6                                                       
          TET0(1,I,J)=(TET(I,J)-1)/4                                        
          TET0(2,I,J)=(TET(I,J)-TET0(1,I,J)*4-1)/2                          
          TET0(3,I,J)=TET(I,J)-TET0(1,I,J)*4-TET0(2,I,J)*2-1                
        ENDDO
      ENDDO
                             CALL TRACE$POP()
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_TETCNT(NMSHP,NUM,TET0,N,INV,THIS)
!     **************************************************************************
!     **  TETCNT CALCULATES ALL DIFFERENT TETRAHEDRA AND COUNTS THEM          **
!     **  INPUT :                                                             **
!     **    NMSHP       NUMBER OF SUBLATTICE POINTS INSIDE AND                **
!     **                ON ALL FACES OF A REC. UNIT CELL                      **
!     **    NUM(I)      MAPPING FROM A GENERAL POINT (I) TO THE               **
!     **                CORRESPONDING IRREDUCIBLE POINT (NUM)                 **
!     **    TET0(I,J,K) COORDINATES (I) IN THE BASIS OF SUBLATTICE            **
!     **                VECTORS, OF THE 4 CORNERS (J) OF A                    **
!     **                TETRAHEDRON (K) IN A SUBLATTICE UNIT CELL             **
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS           **
!     **    INV         DUMMY NUMBER (MUST BE 0)                              **
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                        **
!     **    MWRIT       INFORMATION FOR MWRIT TETRAHEDRA ARE WRITTEN          **
!     **                AT ONE TIME.                                          **
!     **                                                                      **
!     **************************************************************************
      USE BRILLOUIN_MODULE, ONLY : THIS_TYPE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NMSHP        ! #(GENERAL K-POINTS +FACES)
      INTEGER(4),INTENT(IN) :: NUM(NMSHP)
      INTEGER(4),INTENT(IN) :: N(3)         ! DIV.OF LATTICE VECTORS
      INTEGER(4),INTENT(IN) :: TET0(3,4,6)
      INTEGER(4),INTENT(IN) :: INV
      TYPE(THIS_TYPE),INTENT(INOUT) :: THIS
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
      INTEGER(4),ALLOCATABLE:: ITET(:,:)
      INTEGER(4),ALLOCATABLE:: IY(:,:)
      INTEGER(4)            :: NTT    !#(IRREDUCIBLE TETRAHEDRA)
      INTEGER(4)            :: NTMAX
      INTEGER(4)            :: IPP
      INTEGER(4)            :: NTET
      REAL(8)               :: SUMA
      INTEGER(4)            :: IND,IP
      INTEGER(4)            :: ISVAR1,ISVAR2
      INTEGER(4)            :: IXX
      INTEGER(4)            :: I,J,K,L,M,K1,K2,K3
!     ******************************************************************
                           CALL TRACE$PUSH('BRILLOUIN_TETCNT')
      ALLOCATE(ITET(4,6*N(3)))
      ALLOCATE(IY(4,6*NMSHP))
      NTMAX=N(1)*N(2)*N(3)*6                                            
      NTMAX=NTMAX/(1+INV)                                               
      IPP=0                                                             
      DO K1=1,N(1)                                                   
        DO K2=1,N(2)                                                   
          IND=0                                                             
          DO K3=1,N(3)                                                   
            IP=K3+(N(3)+1)*((K2-1)+(N(2)+1)*(K1-1))                           
            DO I=1,6      ! SUM OVER TETRAHEDRA IN THIS PARALLELEPIPED
              IND=IND+1                                                         
              DO J=1,4    ! SUM OVER TETRAHEDRON EDGES
                IXX=TET0(1,J,I)*(N(2)+1)*(N(3)+1)  &
     &             +TET0(2,J,I)*(N(3)+1) &
     &             +TET0(3,J,I)
                ITET(J,IND)=IP+IXX
              ENDDO
            ENDDO
          ENDDO
                                                                        
!     --  TRANSFORM THE EDGEPOINTS ONTO THE IRREDUCIBLE POINTS          
          DO M=1,4                                                       
            DO J=1,N(3)*6                                                  
              ITET(M,J)=NUM(ITET(M,J))                                          
            ENDDO
          ENDDO
                                                                        
!     --  ORDER THE POINTS OF EACH TETR. ACC. TO INREASING NUMBER       
          DO K=1,3                                                       
            DO J=K+1,4                                                     
              DO L=1,N(3)*6                                                  
                ISVAR1=ITET(K,L)
                ISVAR2=ITET(J,L)
                ITET(K,L)=MIN(ISVAR1,ISVAR2)
                ITET(J,L)=MAX(ISVAR1,ISVAR2)
              ENDDO
            ENDDO
          ENDDO
                                                                        
!     --  IDENTIFY THE TETRAHEDRA WITH INTEGERS                         
          DO M=1,N(3)*6                                                  
            IPP=IPP+1                                                         
            IY(:,IPP)=ITET(:,M)
          ENDDO
          IF(IPP.GE.NTMAX) THEN                                             
            IPP=NTMAX                                                       
            GOTO 100                                                        
          END IF                                                            
!*                                                                      
        ENDDO
      ENDDO
      CALL ERROR$MSG('UNNORMAL END OF LOOP')
      CALL ERROR$STOP('BRILLOUIN_TETCNT')
                                                                        
100   CONTINUE                                                          
      DEALLOCATE(ITET)
!     ------------------------------------------------------------------
!     --  ORDER TETRAHEDRA                                            --
!     ------------------------------------------------------------------
      NTET=IPP                                                          
      CALL BRILLOUIN_ORD1(NTET,IY)
!
!     ==================================================================
!     == DETERMINE NUMBER OF IRREDUCIBLE TETRAHEDRA                   ==
!     ==================================================================
      NTT=1
      DO I=1,NTET-1
        IF(IY(4,I+1).NE.IY(4,I)) THEN
          NTT=NTT+1
        ELSE
          IF(IY(3,I+1).NE.IY(3,I)) THEN
            NTT=NTT+1
          ELSE
            IF(IY(2,I+1).NE.IY(2,I)) THEN
              NTT=NTT+1
            ELSE
              IF(IY(1,I+1).NE.IY(1,I)) THEN
                NTT=NTT+1
              END IF
            END IF
          END IF
        END IF
      ENDDO     
      WRITE(*,FMT='("NUMBER OF DIFFERENT TETRAHEDRA :",I5)')NTT                 
      THIS%NTET=NTT
!
!     ==================================================================
!     == COLLECT IRREDUCIBLE TETRAHEDRA                               ==
!     ==================================================================
      THIS%VOL =1.D0/REAL(6*N(1)*N(2)*N(3),KIND=8)
      ALLOCATE(THIS%MULT(NTT))
      ALLOCATE(THIS%IKP(4,NTT))
      THIS%MULT(:)=0
      NTT=1
      THIS%MULT(1)=1+INV
      THIS%IKP(:,1)=IY(:,1)
      DO I=1,NTET-1                                                 
        IF(IY(1,I+1).EQ.IY(1,I)) THEN
          IF(IY(2,I+1).EQ.IY(2,I)) THEN
            IF(IY(3,I+1).EQ.IY(3,I)) THEN
              IF(IY(4,I+1).EQ.IY(4,I)) THEN
                THIS%MULT(NTT)=THIS%MULT(NTT)+1+INV
              ELSE IF(IY(4,I+1).GT.IY(4,I)) THEN
                NTT=NTT+1
                THIS%MULT(NTT)=1+INV
                THIS%IKP(:,NTT)=IY(:,I+1)
              ELSE 
                CALL ERROR$MSG('ORDERING FAILED')
                CALL ERROR$STOP('BRILLOUIN_TETCNT')
              END IF
            ELSE IF(IY(3,I+1).GT.IY(3,I)) THEN
              NTT=NTT+1
              THIS%MULT(NTT)=1+INV
              THIS%IKP(:,NTT)=IY(:,I+1)
            ELSE 
              CALL ERROR$MSG('ORDERING FAILED')
              CALL ERROR$STOP('BRILLOUIN_TETCNT')
            END IF
          ELSE IF(IY(2,I+1).GT.IY(2,I)) THEN
            NTT=NTT+1
            THIS%MULT(NTT)=1+INV
            THIS%IKP(:,NTT)=IY(:,I+1)
          ELSE 
            CALL ERROR$MSG('ORDERING FAILED')
            CALL ERROR$STOP('BRILLOUIN_TETCNT')
          END IF
        ELSE IF(IY(1,I+1).GT.IY(1,I)) THEN
          NTT=NTT+1
          THIS%MULT(NTT)=1+INV
          THIS%IKP(:,NTT)=IY(:,I+1)
        ELSE 
          CALL ERROR$MSG('ORDERING FAILED')
          CALL ERROR$STOP('BRILLOUIN_TETCNT')
        END IF
      ENDDO
      DEALLOCATE(IY)
                                                                        
!     ------------------------------------------------------------------
!     --  CHECK SUMRULE                                               --
!     ------------------------------------------------------------------
      SUMA=0.D0
      DO I=1,THIS%NTET
        SUMA=SUMA+THIS%VOL*REAL(THIS%MULT(I),KIND=8)
      ENDDO
      IF(ABS(SUMA-1.D0).GT.1.D-5) THEN                                  
        CALL ERROR$MSG('SUMRULE NOT FULLFILLED')
        CALL ERROR$MSG('SUM SHOULD BE EQUAL TO 1')
        CALL ERROR$R8VAL('SUM',SUMA)
        CALL ERROR$I4VAL('THIS%NTET',THIS%NTET)
        CALL ERROR$R8VAL('THIS%VOL',THIS%VOL)
        CALL ERROR$I4VAL('N(1)',N(1))
        CALL ERROR$I4VAL('N(2)',N(2))
        CALL ERROR$I4VAL('N(3)',N(3))
        CALL ERROR$STOP('BRILLOUIN_TETCNT')
      END IF                                                            
                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_ORD1(NMAX,IX)
!     ******************************************************************
!     **                                                              **
!     **  ORD1 ORDERES THE ARRAY IX WITH SIZE NMAX ACCORDING TO       **
!     **  INCREASING NUMBER                                           **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NMAX
      INTEGER(4),INTENT(INOUT):: IX(4,NMAX)
      INTEGER(4)              :: I,J,ISVAR1,ISVAR2
      INTEGER(4)              :: IX1(4)
      INTEGER(4)              :: I1A,I2A,I1B,I2B,I1C,I2C
!     ******************************************************************
      DO I=1,NMAX-1                                               
        DO J=I+1,NMAX                                               
          ISVAR1=IX(1,I)                                                    
          ISVAR2=IX(1,J)                                                    
          IF(ISVAR2.LT.ISVAR1) THEN
            IX1(:)=IX(:,I)
            IX(:,I)=IX(:,J)
            IX(:,J)=IX1(:)
          END IF
        ENDDO
      ENDDO
      I1A=1
      DO
        DO I=I1A,NMAX
          IF(IX(1,I).NE.IX(1,I1A)) EXIT
          I2A=I
        ENDDO
        DO I=I1A,I2A-1
          DO J=I+1,I2A
            ISVAR1=IX(2,I)                                                    
            ISVAR2=IX(2,J)                                                    
            IF(ISVAR2.LT.ISVAR1) THEN
              IX1(:)=IX(:,I)
              IX(:,I)=IX(:,J)
              IX(:,J)=IX1(:)
            END IF
          ENDDO
        ENDDO
        I1B=I1A
        DO
          DO I=I1B,I2A
            IF(IX(2,I).NE.IX(2,I1B)) EXIT
            I2B=I
          ENDDO
          DO I=I1B,I2B-1
            DO J=I+1,I2B
              ISVAR1=IX(3,I)                                                    
              ISVAR2=IX(3,J)                                                    
              IF(ISVAR2.LT.ISVAR1) THEN
                IX1(:)=IX(:,I)
                IX(:,I)=IX(:,J)
                IX(:,J)=IX1(:)
              END IF
            ENDDO
          ENDDO
          I1C=I1B
          DO
            DO I=I1C,I2B
              IF(IX(3,I).NE.IX(3,I1C)) EXIT
              I2C=I
            ENDDO
            DO I=I1C,I2C-1
              DO J=I+1,I2C
                ISVAR1=IX(4,I)
                ISVAR2=IX(4,J)
                IF(ISVAR2.LT.ISVAR1) THEN
                  IX1(:)=IX(:,I)
                  IX(:,I)=IX(:,J)
                  IX(:,J)=IX1(:)
                END IF
              ENDDO
            ENDDO
            IF(I2C.EQ.I2B) EXIT
            I1C=I2C+1
          ENDDO     
          IF(I2B.EQ.I2A) EXIT
          I1B=I2B+1
        ENDDO
        IF(I2A.EQ.NMAX) EXIT
        I1A=I2A+1
      ENDDO
      RETURN                                                            
      END                                                               
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ----                                                         -----
!     ----  BLOCK KINTR                                            -----
!     ----  CALCULATION OF SAMPLING WEIGHTS                        -----
!     ----                                                         -----
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$DOS(NB,NKP,EB,WGHT,RNTOT,EF)
!     **************************************************************************
!     **  CALCULATES THE SAMPLING WEIGHTS FROM TETRAHEDRON INTEGRATION        **
!     **                                                                      **
!     **  INPUT :                                                             **
!     **    NB          NUMBER OF BANDS                                       **
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                        **
!     **    EB          ENERGIES ( DETERMINE FERMI SURFACE )                  **
!     **    RNTOT       NUMBER OF OCCUPIED STATES                             **
!     **  OUTPUT :                                                            **
!     **    WGHT        SAMPLING WEIGHTS                                      **
!     **    EF          FERMI LEVEL                                           **
!     **                                                                      **
!     **  USES THE TETRAHEDRON INFORMATION ON THE BRILLOUIN_MODULE            **
!     **  PRODUCED BY BRILLOUIN$MSG OR BRILLOUIN$MSHNOSYM                     **
!     **                                                                      **
!     **  AUTHOR : PETER E. BLOECHL                                           **
!     **                                                                      **
!     **  SUBROUTINES USED:                                                   **
!     **  EFERMI,SAMFAC,TOTNOS,EFI,WEIGHT                                     **
!     **                                                                      **
!     **************************************************************************
      USE BRILLOUIN_MODULE, ONLY : THIS
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NB          ! #(BANDS)
      INTEGER(4),INTENT(IN) :: NKP         ! #(IRREDUCIBLE K-POINTS)
      REAL(8)   ,INTENT(IN) :: EB(NB*NKP)  ! ENERGY BANDS
      REAL(8)   ,INTENT(OUT):: WGHT(NB*NKP)! INTEGRATION WEIGHTS
      REAL(8)   ,INTENT(IN) :: RNTOT       ! #(ELECTRONS)
      REAL(8)   ,INTENT(OUT):: EF          ! FERMI LEVEL
      REAL(8)   ,PARAMETER  :: TOLMAX=1.D-5
      LOGICAL(4)            :: TCHECK=.TRUE.      
!     **************************************************************************
                       CALL TRACE$PUSH('BRILLOUIN$DOS')
!
!     ==========================================================================
!     --  CALCULATE FERMI LEVEL                                       --
!     ==========================================================================
      CALL BRILLOUIN_EFERMI(RNTOT,EF,TOLMAX,NB,NKP,EB,THIS)
      PRINT*,'FERMI ENERGY AT ',EF                                      
!
!     ==========================================================================
!     ==  CALCULATE WEIGHTS                                                   ==
!     ==========================================================================
      CALL BRILLOUIN_SAMFAC(NB,NKP,EB,EF,WGHT,THIS)
!      CALL BRILLOUIN$WGHT(NKP,NB,EF,EB,WGHT)  ! NEW VERSION. NEEDS TO BE TESTED
!                                                                       
!     ==========================================================================
!     --  CHECK WHETHER SUMRULE IS FULLFILLED                         --
!     ==========================================================================
      IF(TCHECK) THEN 
        IF(ABS(SUM(WGHT)-RNTOT).GT.TOLMAX) THEN                                
          CALL ERROR$MSG('INTEGRATION FAILED')
          CALL ERROR$R8VAL('RESULT OF INTEGRATION: ',SUM(WGHT))
          CALL ERROR$R8VAL('SHOULD BE: ',RNTOT)
          CALL ERROR$R8VAL('TOLMAX',TOLMAX)
          CALL ERROR$STOP('BRILLOUIN_DOS')
        END IF                                                            
      END IF
                                CALL TRACE$POP()
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$WDOS(NB,NKP,EB,EF,WGHT)
!     **************************************************************************
!     **  CALCULATES THE SAMPLING WEIGHTS FOR THE DENSITY OF STATES AT        **
!     **  ENERGY EF                                                           **
!     **                                                                      **
!     **  INPUT :                                                             **
!     **    NB          NUMBER OF BANDS                                       **
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                        **
!     **    EB          ENERGIES ( DETERMINE FERMI SURFACE )                  **
!     **    NE          #(POINTS IN ENERGY GRID)                              **
!     **    EI          GRID ENERGIES                                         **
!     **  OUTPUT :                                                            **
!     **    WGHT        SAMPLING WEIGHTS                                      **
!     **                                                                      **
!     **  AUTHOR : PETER E. BLOECHL                                           **
!     **                                                                      **
!     **  SUBROUTINES USED:                                                   **
!     **  EFERMI,SAMFAC,TOTNOS,EFI,WEIGHT                                     **
!     **                                                                      **
!     **************************************************************************
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NB          ! #(BANDS)
      INTEGER(4),INTENT(IN) :: NKP         ! #(IRREDUCIBLE K-POINTS)
      REAL(8)   ,INTENT(IN) :: EB(NB,NKP)  ! ENERGY BANDS
      REAL(8)   ,INTENT(IN) :: EF
      REAL(8)   ,INTENT(OUT):: WGHT(NB,NKP)! DOS WEIGHTS
!     **************************************************************************
                       CALL TRACE$PUSH('BRILLOUIN$WDOS')
!
!     ==========================================================================
!     ==  CALCULATE WEIGHTS                                                   ==
!     ==========================================================================
      CALL BRILLOUIN_SAMDOS(NB,NKP,EB,EF,WGHT,THIS)
                               CALL TRACE$POP()
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_EFERMI(RNTOT,EF,TOLMAX,NB,NKP,EB,THIS)
!     **************************************************************************
!     **  CALCUALTES THE FERMILEVEL BY INTEGRATION OF THE TOTAL               **
!     **  DENSITY OF STATES                                                   **
!     **                                                                      **
!     **  INPUT :                                                             **
!     **    RNTOT       NUMBER OF OCCUPIED STATES                             **
!     **    NB          NUMBER OF BANDS                                       **
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                        **
!     **    EB          ENERGIES ( DETERMINE FERMI SURFACE )                  **
!     **    TOLMAX      TOLERANCE IN THE NUMBER OF STATES AT EF               **
!     **  OUTPUT :                                                            **
!     **    EF          FERMI LEVEL                                           **
!     **                                                                      **
!     **************************************************************************
      USE BRILLOUIN_MODULE,ONLY : THIS_TYPE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: NKP
      REAL(8)   ,INTENT(IN) :: RNTOT
      REAL(8)   ,INTENT(OUT):: EF
      REAL(8)   ,INTENT(IN) :: TOLMAX
      REAL(8)   ,INTENT(IN) :: EB(NB,NKP)      
      TYPE(THIS_TYPE),INTENT(IN) :: THIS
      REAL(8)               :: E(4)
      INTEGER(4)            :: IKP(4)
      INTEGER(4),PARAMETER  :: NP=1000
      INTEGER(4)            :: IB,I,ILOOP,IP
      REAL(8)               :: DE,EMIN,EMAX
      REAL(8)   ,ALLOCATABLE:: NOS(:)
      REAL(8)   ,ALLOCATABLE:: SOS(:)
      INTEGER(4)            :: ITET
      INTEGER(4)            :: NTET
      REAL(8)               :: SUMA
      REAL(8)               :: TOL
      REAL(8)               :: VOL
      REAL(8)               :: ESTEP
!     **************************************************************************
                       CALL TRACE$PUSH('BRILLOUIN_EFERMI')
!     ------------------------------------------------------------------
!     --  FIND EMIN EMAX     (ENERGYBANDS ARE ASSUMED                 --
!     --                      TO BE ORDERED WITH RESPECT TO ENERGY    --
!     ------------------------------------------------------------------
      EMIN=MINVAL(EB(:,:))
      EMAX=MAXVAL(EB(:,:))
      DE=(EMAX-EMIN)/REAL(NP-1,KIND=8)  !NOT THE REAL GRID-SPACING!
      EMAX=EMAX+DE                                                      
      EMIN=EMIN-DE                                                      
      ALLOCATE(NOS(NP))
      ALLOCATE(SOS(NP))
      NTET=THIS%NTET
      ILOOP=0                                                           
1000  CONTINUE                                                          
      ILOOP=ILOOP+1                                                     
      NOS(:)=0.D0
      SOS(:)=0.D0
      SUMA=0.D0                                                          
      DO ITET=1,NTET                                                
        VOL=THIS%VOL*REAL(THIS%MULT(ITET),KIND=8)
        IKP(:)=THIS%IKP(:,ITET)
        SUMA=SUMA+VOL                                                       
        DO IB=1,NB                                                    
          DO I=1,4                                                      
            E(I)=EB(IB,IKP(I))                                                
          ENDDO
          CALL BRILLOUIN_TOTNOS(VOL,E,EMIN,EMAX,NP,NOS,SOS)                   
        ENDDO
      ENDDO
      IF(ABS(SUMA-1.D0).GT.1.D-5) THEN                                  
        CALL ERROR$MSG('TETRAHEDRA DO NOT FILL VOLUME')
        CALL ERROR$MSG('SUM SHOULD BE 1.')
        CALL ERROR$R8VAL('SUM ',SUMA)
        CALL ERROR$STOP('BRILLOUIN_EFERMI')
      END IF                                                            
                                                                        
!     ------------------------------------------------------------------
!     --  GET FERMI LEVEL                                             --
!     ------------------------------------------------------------------
      TOL=TOLMAX
      CALL BRILLOUIN_EFI(RNTOT,TOL,EF,EMIN,EMAX,NP,NOS,SOS)               
!                                                                        
!     ------------------------------------------------------------------
!     --  CHECK ACCURACY AND RESTART IF NECCESARY                     --
!     ------------------------------------------------------------------
      IF(TOL.GT.TOLMAX) THEN                                            
        ESTEP=(EMAX-EMIN)/REAL(NP-1,KIND=8)                                    
        IP=1+INT((EF-EMIN)/ESTEP)                                            
        EMIN=EMIN+ESTEP*REAL(IP-1,KIND=8)                                      
        EMAX=EMIN+ESTEP                                                 
        IF(RNTOT-NOS(IP).LE.TOLMAX) THEN                      
          EF=EMIN                                                       
          GOTO 2000                                                     
        ELSE IF(NOS(IP+1)-RNTOT.LE.TOLMAX) THEN               
          EF=EMAX                                                       
          GOTO 2000                                                     
        END IF                                                          
        IF(ILOOP.GT.5) THEN                                             
          CALL ERROR$MSG('CANNOT FIND FERMI LEVEL')
          CALL ERROR$R8VAL('TOL',TOL)
          CALL ERROR$R8VAL('EMIN',EMIN)
          CALL ERROR$R8VAL('EMAX',EMAX)
          CALL ERROR$R8VAL('NOS(EMIN)',NOS(IP))
          CALL ERROR$R8VAL('NOS(EMAX)',NOS(IP+1))
          CALL ERROR$STOP('BRILLOUIN_EFERMI')
        END IF                                                          
!        PRINT*,'ILOOP ',ILOOP                                           
        GOTO 1000                                                       
      END IF                                                            
2000  CONTINUE                                                          
      DEALLOCATE(NOS)
      DEALLOCATE(SOS)
                                  CALL TRACE$POP()
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_EFI(RNTOT,TOL,EFERMI,EMIN,EMAX,NP,NOS,SOS)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NP     ! #(ENERGY GRID POINTS)
      REAL(8)   ,INTENT(IN)   :: RNTOT  ! #(ELECTRONS)
      REAL(8)   ,INTENT(INOUT):: TOL
      REAL(8)   ,INTENT(IN)   :: SOS(NP)!
      REAL(8)   ,INTENT(INOUT):: NOS(NP)
      REAL(8)   ,INTENT(OUT)  :: EFERMI
      REAL(8)                 :: NOSUP,NOSLOW,NOSIP
      REAL(8)                 :: DNOS
      REAL(8)                 :: ELOW
      REAL(8)                 :: ADD
      REAL(8)                 :: ESTEP,EMIN,EMAX
      INTEGER(4)              :: IP,IPUP,IPLOW
      INTEGER(4)              :: IFIND,I
!     **************************************************************************
                       CALL TRACE$PUSH('BRILLOUIN_EFI')
      ADD=0.D0                                                          
      DO I=1,NP                                                     
        ADD=ADD+SOS(I)                                                    
        NOS(I)=NOS(I)+ADD                                                 
      ENDDO
!++++  NOS(I) IS THE INTEGRATED DENSITY OF STATES
      IF(NOS(1).GT.RNTOT+.5D0*TOL.OR.NOS(NP).LT.RNTOT-.5D0*TOL) THEN    
        CALL ERROR$MSG('EFERMI OUT OF ENERGY RANGE')
        CALL ERROR$R8VAL('ENERGY OF LOWER BOUND',EMIN)
        CALL ERROR$R8VAL('NUMBER OF STATES AT THE LOWER BOUND',NOS(1))         
        CALL ERROR$R8VAL('ENERGY OF UPPER BOUND',EMAX)           
        CALL ERROR$R8VAL('NUMBER OF STATES AT THE UPPER BOUND',NOS(NP))        
        CALL ERROR$R8VAL('NUMBER OF STATES REQUESTED',RNTOT)          
        CALL ERROR$R8VAL('ADD ',ADD)
        CALL ERROR$MSG('YOU PROBABLY HAVE TO INCREASE THE NUMBER OF BANDS NB')
        CALL ERROR$STOP('BRILLOUIN_EFI')
      END IF                                                            
      IPUP=NP                                                           
      IPLOW=1                                                           
      NOSLOW=NOS(1)
      NOSUP=NOS(NP)
      DO IFIND=1,NP                                                 
        IP=IPLOW+INT(0.5*REAL(IPUP-IPLOW))
        NOSIP=NOS(IP)                                                     
        IF(RNTOT.GT.NOSIP) THEN                                         
          IPLOW=IP                                                        
          NOSLOW=NOSIP                                                    
        ELSE                                                              
          IPUP=IP                                                         
          NOSUP=NOSIP                                                     
        END IF                                                            
        IF(IPUP.EQ.IPLOW+1) GOTO 300                                      
      ENDDO
      CALL ERROR$MSG('EFERMI NOT FOUND')
      CALL ERROR$STOP('BRILLOUIN_EFI')
300   CONTINUE                                                          
      TOL=NOSUP-NOSLOW                                                  
      ESTEP=(EMAX-EMIN)/REAL(NP-1,KIND=8)
      ELOW=EMIN+REAL(IPLOW-1,KIND=8)*ESTEP                                     
      DNOS=NOSUP-NOSLOW                                                 
      IF(DNOS.NE.0.D0) THEN                                             
        EFERMI=ELOW+(RNTOT-NOSLOW)/(NOSUP-NOSLOW)*ESTEP                 
      ELSE                                                              
        EFERMI=ELOW                                                     
      END IF                                                            
      IF(EFERMI.LT.ELOW) PRINT*,'ERROR IN EFI '                       
                            CALL TRACE$POP()
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_TOTNOS(VOL,E_,EMIN,EMAX,NP,NOS,SOS)
!     **************************************************************************
!     **  CALCULATES THE INTEGRATED DOS                                       **
!     **  FOR ONE TETRAHEDRON ON THE ENERGYMESH                               **
!     **  INPUT :                                                             **
!     **    VOL         WEIGHT OF THIS TETRAHEDRON                            **
!     **    E           ENERGIES AT THE EDGEPOINTS                            **
!     **    EMIN        MINIMUM VALUE OF ENERGY MESH FOR NOS AND SOS          **
!     **    EMAX        MAXIMUM VALUE OF ENERGY MESH FOR NOS AND SOS          **
!     **    NP          NUMBER OF POINTS ON THE ENERGY MESH                   **
!     **  OUTPUT:                                                             **
!     **    SOS         > NOS(E)+ SUM OVER E: SOS(E) =                        **
!     **    NOS         > NUMBER OF STATES BELOW E                            **
!     **                                                                      **
!     **    SOS(I) CONTAINS THE WEIGHTS OFF ALL TETRAHEDRA WITH THE           **
!     **    HIGHEST ENERGY EX IN THE WINDOW E(I-1)<EX<E(I).                   **
!     **    SOS(1) CONTAINS THE WEIGHTS OF ALL TETRAHEDRA WITH THEIR          **
!     **    HIGHEST ENERGY EX BELOW E(I)                                      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: VOL
      REAL(8)   ,INTENT(IN)   :: E_(4)
      REAL(8)   ,INTENT(IN)   :: EMIN
      REAL(8)   ,INTENT(IN)   :: EMAX
      INTEGER(4),INTENT(IN)   :: NP
      REAL(8)   ,INTENT(INOUT):: NOS(NP)
      REAL(8)   ,INTENT(INOUT):: SOS(NP)
      REAL(8)                 :: E(4)
      REAL(8)                 :: SVAR1,SVAR2
      REAL(8)                 :: E21,E31,E41,E32,E42,E43
      REAL(8)                 :: ESTEP,DE,EN
      REAL(8)                 :: A,B,C,D
      INTEGER(4)              :: IMIN,IMAX
      INTEGER(4)              :: I,J
      REAL(8)                 :: EIND(4)
      INTEGER(4)              :: IND(4)
!     **************************************************************************
      E(:)=E_(:)
!
!     ------------------------------------------------------------------
!     --  INTEGRATION WITHOUT FERMISURFACE                            --
!     ------------------------------------------------------------------
      IF(MINVAL(E).GE.EMAX) RETURN
      IF(MAXVAL(E).LE.EMIN) THEN                                                
        SOS(1)=SOS(1)+VOL                                               
        RETURN                                                          
      END IF                                                            
!     ------------------------------------------------------------------
!     --  ORDER ENERGIES                                              --
!     ------------------------------------------------------------------
      DO I=1,3                                                      
        DO J=I+1,4                                                    
          SVAR1=MIN(E(I),E(J))                                            
          SVAR2=MAX(E(I),E(J))                                            
          E(I)=SVAR1                                                        
          E(J)=SVAR2                                                        
        ENDDO
      ENDDO
!
!     ==================================================================
!     == CONSTRUCT THE INDICES OF THE GRID POINTS JUST BELOW THE ENERGIES
!     == THIS CONSTRUCTION IS NECESSARY, BECAUSE AT LEAST IN THE       ==
!     == G95 COMPILER LARGE NEGATIVE INTEGERS RESULT FROM              ==
!     == THE INT FUNCTION APPLIED TO A VERY LARGE POSITIVE REAL NUMBER ==
!     ==================================================================
      ESTEP=(EMAX-EMIN)/REAL(NP-1,KIND=8)                                      
      EIND(:)=MIN(E(:),EMAX+ESTEP)
      EIND(:)=MAX(EIND(:),EMIN-ESTEP)
      IND(:)=INT(2.D0+(EIND(:)-EMIN)/ESTEP)-1
      IND(:)=MAX(0,IND(:))
      IND(:)=MIN(NP,IND(:))
!     ------------------------------------------------------------------
!     --  CALCULATE UNCORRECTED INTEGRAL AS MEANVALUE                 --
!     ------------------------------------------------------------------
      E21=E(2)-E(1)                                                     
      E31=E(3)-E(1)                                                     
      E41=E(4)-E(1)                                                     
      E32=E(3)-E(2)                                                     
      E42=E(4)-E(2)                                                     
      E43=E(4)-E(3)                                                     
!     == ENERGY WINDOW FROM E1 TO E2 ====================================
      IMIN=IND(1)+1                !E(IMIN) IS THE FIRST POINT IN [E1,E2]
      IMAX=IND(2)                  !E(IMAX) IS THE LAST POINT IN [E1,E2]
      EN=EMIN+ESTEP*REAL(IMIN-1,KIND=8)
      IF(IMAX.GE.IMIN) THEN                                             
        A=VOL/(E21*E31*E41)                                             
        DO I=IMIN,IMAX                                              
          NOS(I)=NOS(I)+A*(EN-E(1))**3                                    
          EN=EN+ESTEP                                                     
        ENDDO
      END IF                                                            
!     == ENERGY WINDOW FROM E2 TO E3 ====================================
      IMIN=IND(2)+1              !E(IMIN) IS THE FIRST POINT IN [E2,E3]
      IMAX=IND(3)                !E(IMAX) IS THE LAST POINT IN [E2,E3]
      IF(IMAX.GE.IMIN) THEN                                             
        A=VOL*E21**2/(E31*E41)                                          
        B=3.D0*VOL*E21/(E31*E41)                                        
        C=3.D0*VOL/(E31*E41)                                            
        D=-VOL/(E32*E41*E31*E42)*(E31+E42)                              
        DO I=IMIN,IMAX                                              
          DE=EN-E(2)                                                      
!         NOS(I)=NOS(I)+A+B*DE+C*DE**2+D*DE**3                            
          NOS(I)=NOS(I)+A+DE*(B+DE*(C+D*DE))                              
          EN=EN+ESTEP                                                     
        ENDDO
      END IF                                                            
!     == ENERGY WINDOW FROM E3 TO E4 ====================================
      IMIN=IND(3)+1                !E(IMIN) IS THE FIRST POINT IN [E3,E4]
      IMAX=IND(4)                  !E(IMAX) IS THE LAST POINT IN [E3,E4]
      IF(E43.GT.0.D0) THEN                                              
        A=VOL                                                           
        D=VOL/(E41*E42*E43)                                             
        DO I=IMIN,IMAX                                              
          NOS(I)=NOS(I)+A+D*(EN-E(4))**3                                  
          EN=EN+ESTEP                                                     
        ENDDO
      END IF                                                            
      IMIN=IND(4)+1 !E(IMIN) IS FIRST POINT IN [E4,INFTY]
      IF(IMIN.GT.NP) RETURN                                             
      SOS(IMIN)=SOS(IMIN)+VOL                                           
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_SAMFAC(NB,NKP,EB,EF,WGHT,THIS)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES SAMPLING WEIGHTS                                 **
!     **  INPUT :                                                     **
!     **    NB          NUMBER OF BANDS                               **
!     **    NKP         NUMBER OF K-POINTS                            **
!     **    EF          FERMI LEVEL                                   **
!     **  OUTPUT :                                                    **
!     **    WGHT        SAMPLING WEIGHTS                              **
!     **                                                              **
!     ******************************************************************
      USE BRILLOUIN_MODULE, ONLY : THIS_TYPE
      IMPLICIT NONE
      INTEGER(4)    ,INTENT(IN)  :: NB
      INTEGER(4)    ,INTENT(IN)  :: NKP
      REAL(8)       ,INTENT(IN)  :: EF
      REAL(8)       ,INTENT(IN)  :: EB(NB,NKP)
      TYPE(THIS_TYPE),INTENT(IN) :: THIS
      REAL(8)       ,INTENT(OUT) :: WGHT(NB,NKP)
      REAL(8)                    :: E(4)
      REAL(8)                    :: WGHT0(4)
      REAL(8)                    :: VOL
      INTEGER(4)                 :: IKP(4)
      INTEGER(4)                 :: ITET,NTET
      INTEGER(4)                 :: IB,I
!     ******************************************************************
                       CALL TRACE$PUSH('BRILLOUIN_SAMFAC')
      WGHT(:,:)=0.D0
      NTET=THIS%NTET
      DO ITET=1,NTET                                                
        VOL=THIS%VOL*THIS%MULT(ITET)
        IKP(:)=THIS%IKP(:,ITET)
        DO IB=1,NB                                                    
          DO I=1,4                                                      
            E(I)=EB(IB,IKP(I))                                                
            WGHT0(I)=0.D0
          ENDDO
          CALL BRILLOUIN_WEIGHT(VOL,E,EF,WGHT0)
          DO I=1,4                                                      
            WGHT(IB,IKP(I))=WGHT(IB,IKP(I))+WGHT0(I)                          
          ENDDO
        ENDDO
      ENDDO
                           CALL TRACE$POP()
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_SAMDOS(NB,NKP,EB,EF,WGHT,THIS)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES SAMPLING WEIGHTS                                 **
!     **  INPUT :                                                     **
!     **    NB          NUMBER OF BANDS                               **
!     **    NKP         NUMBER OF K-POINTS                            **
!     **    EF          FERMI LEVEL                                   **
!     **  OUTPUT :                                                    **
!     **    WGHT        SAMPLING WEIGHTS                              **
!     **                                                              **
!     ******************************************************************
      USE BRILLOUIN_MODULE, ONLY : THIS_TYPE
      IMPLICIT NONE
      INTEGER(4)    ,INTENT(IN)  :: NB
      INTEGER(4)    ,INTENT(IN)  :: NKP
      REAL(8)       ,INTENT(IN)  :: EF
      REAL(8)       ,INTENT(IN)  :: EB(NB,NKP)
      TYPE(THIS_TYPE),INTENT(IN) :: THIS
      REAL(8)       ,INTENT(OUT) :: WGHT(NB,NKP)
      REAL(8)                    :: E(4)
      REAL(8)                    :: WGHT0(4),DWGHT0(4)
      REAL(8)                    :: VOL
      INTEGER(4)                 :: IKP(4)
      INTEGER(4)                 :: ITET,NTET
      INTEGER(4)                 :: IB,I
!     ******************************************************************
                       CALL TRACE$PUSH('BRILLOUIN_SAMFAC')
      WGHT(:,:)=0.D0
      NTET=THIS%NTET
      DO ITET=1,NTET                                                
        VOL=THIS%VOL*THIS%MULT(ITET)
        IKP(:)=THIS%IKP(:,ITET)
        DO IB=1,NB                                                    
          DO I=1,4                                                      
            E(I)=EB(IB,IKP(I))                                                
            WGHT0(I)=0.D0
          ENDDO
          CALL BRILLOUIN_DWEIGHT(VOL,E,EF,WGHT0,DWGHT0)
          DO I=1,4                                                      
            WGHT(IB,IKP(I))=WGHT(IB,IKP(I))+DWGHT0(I)                          
          ENDDO
        ENDDO
      ENDDO
                           CALL TRACE$POP()
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_SAMFACANDDER(NB,NKP,EB,EF,WGHT,DWGHT,THIS)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES SAMPLING WEIGHTS                                 **
!     **  INPUT :                                                     **
!     **    NB          NUMBER OF BANDS                               **
!     **    NKP         NUMBER OF K-POINTS                            **
!     **    EF          FERMI LEVEL                                   **
!     **  OUTPUT :                                                    **
!     **    WGHT        SAMPLING WEIGHTS                              **
!     **                                                              **
!     ******************************************************************
      USE BRILLOUIN_MODULE, ONLY : THIS_TYPE
      IMPLICIT NONE
      INTEGER(4)    ,INTENT(IN)  :: NB
      INTEGER(4)    ,INTENT(IN)  :: NKP
      REAL(8)       ,INTENT(IN)  :: EF
      REAL(8)       ,INTENT(IN)  :: EB(NB,NKP)
      TYPE(THIS_TYPE),INTENT(IN)  :: THIS
      REAL(8)       ,INTENT(OUT) :: WGHT(NB,NKP)
      REAL(8)       ,INTENT(OUT) :: DWGHT(NB,NKP,NKP)
      REAL(8)                    :: E(4)
      REAL(8)                    :: WGHT0(4)
      REAL(8)                    :: DWGHT0(4,4)
      REAL(8)                    :: VOL
      INTEGER(4)                 :: IKP(4)
      INTEGER(4)                 :: ITET,NTET
      INTEGER(4)                 :: IB,I,J
!     ******************************************************************
      WGHT(:,:)=0.D0
      DWGHT(:,:,:)=0.D0
      NTET=THIS%NTET
      DO ITET=1,NTET                                                
        VOL=THIS%VOL*THIS%MULT(ITET)
        IKP(:)=THIS%IKP(:,ITET)
        DO IB=1,NB                                                    
          DO I=1,4                                                      
            E(I)=EB(IB,IKP(I))                                                
            WGHT0(I)=0.D0
            DWGHT0(I,:)=0.D0
          ENDDO
          CALL BRILLOUIN_WEIGHTANDDER(VOL,E,EF,WGHT0,DWGHT0)
          DO I=1,4                                                      
            WGHT(IB,IKP(I))=WGHT(IB,IKP(I))+WGHT0(I)                          
            DO J=1,4
              DWGHT(IB,IKP(I),IKP(J))=DWGHT(IB,IKP(I),IKP(J))+DWGHT0(I,J)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_WEIGHT(VOL,E_,EF,WGHT)
!     **************************************************************************
!     **  CALCULATES THE WEIGHTS FOR TETRAHEDRON-SAMPLING                     **
!     **  CORRESPONDING TO INTEGRATION OVER ONE TETRAHEDRON                   **
!     **                                                                      **
!     **  CORRECTION FOR THE NONLINEAR SHAPE INCLUDED IF ICOR=1               **
!     **                                                                      **
!     **  AUTHOR : P.BLOECHL                                                  **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)  ,INTENT(IN) :: VOL      ! WEIGHT OF THIS TETRAHEDRON
      REAL(8)  ,INTENT(IN) :: EF       ! FERMI LEVEL
      REAL(8)  ,INTENT(IN) :: E_(4)    ! ENERGY BANDS AT THE CORNERS
      REAL(8)  ,INTENT(OUT):: WGHT(4)  ! INTEGRATION WEIGHTS
      INTEGER(4),PARAMETER :: ICOR=1   ! ON/OFF SWITCH FOR CORRECTION 
      REAL(8)              :: E(4)     ! ENERGY BANDS AT THE CORNERS
      REAL(8)              :: FA(4)
      REAL(8)              :: FB(4)
      INTEGER(4)           :: INDEX(4)
      REAL(8)              :: X
      REAL(8)              :: VPRIME
      REAL(8)              :: DE
      REAL(8)              :: DOS
      REAL(8)              :: E21,E31,E41,E32,E42,E43
      REAL(8)              :: DE1,DE2,DE3,DE4
      INTEGER(4)           :: I,J,IP,N,M,K
      REAL(8)              :: DA,DB,DC
      REAL(8)              :: VOL14
!     **************************************************************************
      E(:)=E_(:)
      WGHT(:)=0.D0
!     ------------------------------------------------------------------
!     --  INTEGRATION WITHOUT FERMISURFACE                            --
!     ------------------------------------------------------------------
      IF(MINVAL(E).GE.EF) RETURN
      IF(MAXVAL(E).LE.EF) THEN                                                  
        VPRIME=.25D0*VOL                                                
        WGHT(:)=VPRIME                                                  
        RETURN                                                          
      END IF                                                            
!     ------------------------------------------------------------------
!     --  ORDER ENERGIES                                              --
!     ------------------------------------------------------------------
!     -- INDEX HOLDS THE ORIGINAL POSITION OF THE ENERGIES AND WEIGHTS  
      DO I=1,4                                                      
        INDEX(I)=I                                                        
      ENDDO
      DO I=1,3                                                      
        IP=I                                                              
        DO J=I+1,4                                                    
          IF(E(IP).GT.E(J)) IP=J                                            
        ENDDO
        IF(IP.GT.I) THEN                                                  
          X=E(IP)                                                         
          E(IP)=E(I)                                                      
          E(I)=X                                                          
          K=INDEX(IP)                                                     
          INDEX(IP)=INDEX(I)                                              
          INDEX(I)=K                                                      
        END IF                                                            
      ENDDO
!
!     ------------------------------------------------------------------
!     --  CALCULATE UNCORRECTED INTEGRAL AS MEANVALUE                 --
!     ------------------------------------------------------------------
      E21=E(2)-E(1)                                                     
      E31=E(3)-E(1)                                                     
      E41=E(4)-E(1)                                                     
      E32=E(3)-E(2)                                                     
      E42=E(4)-E(2)                                                     
      E43=E(4)-E(3)                                                     
      IF(EF.GT.E(1).AND.EF.LE.E(2)) THEN                                
        DE=EF-E(1)                                                      
        VPRIME=.25D0*VOL*DE**3/(E21*E31*E41)                            
        WGHT(1)=VPRIME*(4.D0-DE/E21-DE/E31-DE/E41)                      
        WGHT(2)=VPRIME*DE/E21                                           
        WGHT(3)=VPRIME*DE/E31                                           
        WGHT(4)=VPRIME*DE/E41                                           
!       ------  PARAMETERS FOR CORRECION                                
        DOS=3.D0*VPRIME*4.D0/(EF-E(1))                                  
      ELSE IF(EF.GT.E(2).AND.EF.LT.E(3)) THEN                           
        DE1=EF-E(1)                                                     
        DE2=EF-E(2)                                                     
        DE3=E(3)-EF                                                     
        DE4=E(4)-EF                                                     
!       ------  TETRAHEDRON X1,X2,X13P,X14P                             
        VPRIME=VOL*DE1**2/(E41*E31)*.25D0                               
        WGHT(2)=VPRIME                                                  
        WGHT(3)=VPRIME*(DE1/E31)                                        
        WGHT(4)=VPRIME*(DE1/E41)                                        
        WGHT(1)=VPRIME*(3.D0-DE1/E41-DE1/E31)                           
!       ------  TETRAHEDRON X2,X13P,X23P,X14P                           
        VPRIME=.25D0*VOL*DE2*DE3*DE1/(E32*E31*E41)                      
        WGHT(1)=WGHT(1)+VPRIME*(2.D0-DE1/E31-DE1/E41)                   
        WGHT(2)=WGHT(2)+VPRIME*(2.D0-DE2/E32)                           
        WGHT(3)=WGHT(3)+VPRIME*(DE2/E32+DE1/E31)                        
        WGHT(4)=WGHT(4)+VPRIME*(DE1/E41)                                
!       ------  TETRAHEDRON X2,X23P,X24P,X14P                           
        VPRIME=.25D0*VOL*DE2**2*DE4/(E42*E32*E41)                       
        WGHT(1)=WGHT(1)+VPRIME*(1.D0-DE1/E41)                           
        WGHT(2)=WGHT(2)+VPRIME*(3.D0-DE2/E32-DE2/E42)                   
        WGHT(3)=WGHT(3)+VPRIME*(DE2/E32)                                
        WGHT(4)=WGHT(4)+VPRIME*(DE2/E42+DE1/E41)                        
!       ------  DOS=A+B*(EF-E2)+C*(EF-E2)**2                            
        DA=3.D0*VOL*E21/(E31*E41)                                       
        DB=6.D0*VOL/(E31*E41)                                           
        DC=-3.D0*VOL/(E32*E41*E31*E42)*(E31+E42)                        
        DOS=DA+DB*DE2+DC*DE2**2                                         
      ELSE IF(EF.GE.E(3).AND.EF.LT.E(4)) THEN                           
        DE=E(4)-EF                                                      
        VPRIME=.25D0*VOL*DE**3/(E41*E42*E43)                            
        VOL14=.25D0*VOL                                                 
        WGHT(1)=VOL14-VPRIME*DE/E41                                     
        WGHT(2)=VOL14-VPRIME*DE/E42                                     
        WGHT(3)=VOL14-VPRIME*DE/E43                                     
        WGHT(4)=VOL14-VPRIME*(4.D0-DE/E41-DE/E42-DE/E43)                
!       ------  PARAMETERS FOR CORRECION                                
        DOS=3.D0*VPRIME*4.D0/(E(4)-EF)                                  
      ELSE                                                              
        CALL ERROR$MSG('ENERGIES OUT OF ORDER OR FERMI LEVEL OUT OF RANGE')
        CALL ERROR$R8VAL('EF',EF)
        CALL ERROR$R8VAL('E1',E(1))
        CALL ERROR$R8VAL('E2',E(2))
        CALL ERROR$R8VAL('E3',E(3))
        CALL ERROR$R8VAL('E4',E(4))
        CALL ERROR$STOP('BRILLOUIN_WEIGHT')
      END IF                                                            
!     ------------------------------------------------------------------
!     --  ADD CORRECTION FOR QUADRATIC DEVIATION                      --
!     ------------------------------------------------------------------
      IF(ICOR.EQ.1) THEN                                                
        DO M=1,4                                                    
          DO N=1,4                                                    
            WGHT(M)=WGHT(M)+.25D0*(E(N)-E(M))*DOS*.1D0                      
          ENDDO
        ENDDO
      END IF                                                            
!     ------------------------------------------------------------------
!     --  REORDER WEIGHTS                                             --
!     ------------------------------------------------------------------
      DO I=1,4                                                      
        FA(INDEX(I))=WGHT(I)                                              
        FB(INDEX(I))=E(I)                                                 
      ENDDO
      DO I=1,4                                                      
        WGHT(I)=FA(I)                                                     
        E(I)=FB(I)                                                        
      ENDDO
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_DWEIGHT(VOL,E_,EF,WGHT,DWGHT)
!     **************************************************************************
!     **  CALCULATES THE WEIGHTS FOR TETRAHEDRON-SAMPLING                     **
!     **  CORRESPONDING TO INTEGRATION OVER ONE TETRAHEDRON                   **
!     **  AND THEIR DERIVATIVES WITH RESPECT TO THE FERMILEVEL                **
!     **                                                                      **
!     **  CORRECTION FOR THE NONLINEAR SHAPE INCLUDED IF ICOR=1               **
!     **                                                                      **
!     **  AUTHOR : P.BLOECHL, GOSLAR 2011                                     **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)  ,INTENT(IN) :: VOL      ! WEIGHT OF THIS TETRAHEDRON
      REAL(8)  ,INTENT(IN) :: EF       ! FERMI LEVEL
      REAL(8)  ,INTENT(IN) :: E_(4)    ! ENERGY BANDS AT THE CORNERS
      REAL(8)  ,INTENT(OUT):: WGHT(4)  ! INTEGRATION WEIGHTS
      REAL(8)  ,INTENT(OUT):: DWGHT(4)  ! ENERGY DERIVATIVE OF I.W.
      INTEGER(4),PARAMETER :: ICOR=0   ! ON/OFF SWITCH FOR CORRECTION 
      REAL(8)              :: E(4)     ! ENERGY BANDS AT THE CORNERS
      REAL(8)              :: FA(4),DFA(4)
      REAL(8)              :: FB(4)
      INTEGER(4)           :: INDEX(4)
      REAL(8)              :: X
      REAL(8)              :: VPRIME,DVPRIME
      REAL(8)              :: DE
      REAL(8)              :: DOS,DDOS
      REAL(8)              :: E21,E31,E41,E32,E42,E43
      REAL(8)              :: DE1,DE2,DE3,DE4
      INTEGER(4)           :: I,J,IP,N,M,K
      REAL(8)              :: DA,DB,DC
      REAL(8)              :: VOL14
!     **************************************************************************
      E(:)=E_(:)
      WGHT(:)=0.D0
      DWGHT(:)=0.D0
!     ==========================================================================
!     ==  INTEGRATION WITHOUT FERMISURFACE                                    ==
!     ==========================================================================
      IF(MINVAL(E).GE.EF) RETURN
      IF(MAXVAL(E).LE.EF) THEN                                                  
        VPRIME=.25D0*VOL                                                
        WGHT(:)=VPRIME                                                  
        DWGHT(:)=0.D0
        RETURN                                                          
      END IF                                                            
!     ==========================================================================
!     ==  ORDER ENERGIES                                                      ==
!     ==========================================================================
!     -- INDEX HOLDS THE ORIGINAL POSITION OF THE ENERGIES AND WEIGHTS  
      DO I=1,4                                                      
        INDEX(I)=I                                                        
      ENDDO
      DO I=1,3                                                      
        IP=I                                                              
        DO J=I+1,4                                                    
          IF(E(IP).GT.E(J)) IP=J                                            
        ENDDO
        IF(IP.GT.I) THEN                                                  
          X=E(IP)                                                         
          E(IP)=E(I)                                                      
          E(I)=X                                                          
          K=INDEX(IP)                                                     
          INDEX(IP)=INDEX(I)                                              
          INDEX(I)=K                                                      
        END IF                                                            
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE UNCORRECTED INTEGRAL AS MEANVALUE                         ==
!     ==========================================================================
      E21=E(2)-E(1)                                                     
      E31=E(3)-E(1)                                                     
      E41=E(4)-E(1)                                                     
      E32=E(3)-E(2)                                                     
      E42=E(4)-E(2)                                                     
      E43=E(4)-E(3)                                                     
      IF(EF.GT.E(1).AND.EF.LE.E(2)) THEN                                
        DE=EF-E(1)                                                      
        VPRIME=.25D0*VOL*DE**3/(E21*E31*E41)                            
        DVPRIME=.25D0*VOL*3.D0*DE**2/(E21*E31*E41)                            
        WGHT(1)=VPRIME*(4.D0-DE/E21-DE/E31-DE/E41)                      
        WGHT(2)=VPRIME*DE/E21                                           
        WGHT(3)=VPRIME*DE/E31                                           
        WGHT(4)=VPRIME*DE/E41                                           
        DWGHT(1)=DVPRIME*(4.D0-DE/E21-DE/E31-DE/E41) &
       &                     +VPRIME*(-1.D0/E21-1.D0/E31-1.D0/E41)
        DWGHT(2)=DVPRIME*DE/E21+VPRIME/E21
        DWGHT(3)=DVPRIME*DE/E31+VPRIME/E31
        DWGHT(4)=DVPRIME*DE/E41+VPRIME/E41
!       ------  PARAMETERS FOR CORRECION                                
        DOS=3.D0*VPRIME*4.D0/(EF-E(1))
        DDOS=-DOS/(EF-E(1))
      ELSE IF(EF.GT.E(2).AND.EF.LT.E(3)) THEN                           
        DE1=EF-E(1)                                                     
        DE2=EF-E(2)                                                     
        DE3=E(3)-EF                                                     
        DE4=E(4)-EF                                                     
!       ------  TETRAHEDRON X1,X2,X13P,X14P                             
        VPRIME=VOL*DE1**2/(E41*E31)*.25D0                               
        DVPRIME=VOL*2.D0*DE1/(E41*E31)*.25D0                               
        WGHT(2)=VPRIME                                                  
        WGHT(3)=VPRIME*(DE1/E31)                                        
        WGHT(4)=VPRIME*(DE1/E41)                                        
        WGHT(1)=VPRIME*(3.D0-DE1/E41-DE1/E31)                           
        DWGHT(2)=DVPRIME                                                  
        DWGHT(3)=DVPRIME*(DE1/E31)+VPRIME/E31
        DWGHT(4)=DVPRIME*(DE1/E41)+VPRIME/E41
        DWGHT(1)=DVPRIME*(3.D0-DE1/E41-DE1/E31)+VPRIME*(-1.D0/E41-1.D0/E31)
!       ------  TETRAHEDRON X2,X13P,X23P,X14P                           
        VPRIME=.25D0*VOL*DE2*DE3*DE1/(E32*E31*E41)                      
        DVPRIME=.25D0*VOL*(DE3*DE1-DE2*DE1+DE2*DE3)/(E32*E31*E41)
        WGHT(1)=WGHT(1)+VPRIME*(2.D0-DE1/E31-DE1/E41)                   
        WGHT(2)=WGHT(2)+VPRIME*(2.D0-DE2/E32)                           
        WGHT(3)=WGHT(3)+VPRIME*(DE2/E32+DE1/E31)                        
        WGHT(4)=WGHT(4)+VPRIME*(DE1/E41)                                
        DWGHT(1)=DWGHT(1)+DVPRIME*(2.D0-DE1/E31-DE1/E41) &
       &                                           +VPRIME*(-1.D0/E31-1.D0/E41)
        DWGHT(2)=DWGHT(2)+DVPRIME*(2.D0-DE2/E32)   +VPRIME*(-1.D0/E32)
        DWGHT(3)=DWGHT(3)+DVPRIME*(DE2/E32+DE1/E31)+VPRIME*(1.D0/E32+1.D0/E31)
        DWGHT(4)=DWGHT(4)+DVPRIME*(DE1/E41)        +VPRIME*(1.D0/E41)
!       ------  TETRAHEDRON X2,X23P,X24P,X14P  ---------------------------------
        VPRIME=.25D0*VOL*DE2**2*DE4/(E42*E32*E41)                       
        DVPRIME=.25D0*VOL*(2.D0*DE2*DE4-DE2**2)/(E42*E32*E41)
        WGHT(1)=WGHT(1)+VPRIME*(1.D0-DE1/E41)                           
        WGHT(2)=WGHT(2)+VPRIME*(3.D0-DE2/E32-DE2/E42)                   
        WGHT(3)=WGHT(3)+VPRIME*(DE2/E32)                                
        WGHT(4)=WGHT(4)+VPRIME*(DE2/E42+DE1/E41)                        
        DWGHT(1)=DWGHT(1)+DVPRIME*(1.D0-DE1/E41)   +VPRIME*(-1.D0/E41)
        DWGHT(2)=DWGHT(2)+DVPRIME*(3.D0-DE2/E32-DE2/E42) &
       &                                           +VPRIME*(-1.D0/E32-1.D0/E42)
        DWGHT(3)=DWGHT(3)+DVPRIME*(DE2/E32)        +VPRIME*(1.D0/E32)
        DWGHT(4)=DWGHT(4)+DVPRIME*(DE2/E42+DE1/E41)+VPRIME*(1.D0/E42+1.D0/E41)
!       ------  DOS=A+B*(EF-E2)+C*(EF-E2)**2  ----------------------------------
        DA=3.D0*VOL*E21/(E31*E41)                                       
        DB=6.D0*VOL/(E31*E41)                                           
        DC=-3.D0*VOL/(E32*E41*E31*E42)*(E31+E42)                        
        DOS=DA+DB*DE2+DC*DE2**2
        DDOS=DB+2.D0*DC*DE2
      ELSE IF(EF.GE.E(3).AND.EF.LT.E(4)) THEN                           
        DE=E(4)-EF                                                      
        VPRIME=.25D0*VOL*DE**3/(E41*E42*E43)                            
        DVPRIME=.25D0*VOL*(-3.D0*DE**2)/(E41*E42*E43)
        VOL14=.25D0*VOL                                                 
        WGHT(1)=VOL14-VPRIME*DE/E41                                     
        WGHT(2)=VOL14-VPRIME*DE/E42                                     
        WGHT(3)=VOL14-VPRIME*DE/E43                                     
        WGHT(4)=VOL14-VPRIME*(4.D0-DE/E41-DE/E42-DE/E43)                
        DWGHT(1)=-DVPRIME*DE/E41+VPRIME/E41
        DWGHT(2)=-DVPRIME*DE/E42+VPRIME/E42
        DWGHT(3)=-DVPRIME*DE/E43+VPRIME/E43                                     
        DWGHT(4)=-DVPRIME*(4.D0-DE/E41-DE/E42-DE/E43) &
       &                      -VPRIME*(1.D0/E41+1.D0/E42+1.D0/E43)
!       ------  PARAMETERS FOR CORRECION                                
        DOS=3.D0*VPRIME*4.D0/(E(4)-EF)
        DDOS=DOS/(E(4)-EF)                                  
      ELSE                                                              
!        CALL ERROR$MSG('ERROR')
!        CALL ERROR$STOP('BRILLOUIN_WEIGHT')
      END IF                                                            
!     ==========================================================================
!     ==  ADD CORRECTION FOR QUADRATIC DEVIATION                              ==
!     ==========================================================================
      IF(ICOR.EQ.1) THEN                                                
        DO M=1,4                                                    
          DO N=1,4                                                    
            WGHT(M) =WGHT(M) +.25D0*(E(N)-E(M))*DOS*.1D0                      
            DWGHT(M)=DWGHT(M)+.25D0*(E(N)-E(M))*DDOS*.1D0                      
          ENDDO
        ENDDO
      END IF                                                            
!     ==========================================================================
!     ==  REORDER WEIGHTS                                                     ==
!     ==========================================================================
      DO I=1,4                                                      
        FA(INDEX(I))=WGHT(I)                                              
        DFA(INDEX(I))=DWGHT(I)                                              
        FB(INDEX(I))=E(I)                                                 
      ENDDO
      DO I=1,4                                                      
        WGHT(I)=FA(I)                                                     
        DWGHT(I)=DFA(I)                                                     
        E(I)=FB(I)                                                        
      ENDDO
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_WEIGHTANDDER(VOL,E_,EF,WGHT,DWGHT)
!     **************************************************************************
!     **  CALCULATES THE WEIGHTS FOR TETRAHEDRON-SAMPLING                     **
!     **  CORRESPONDING TO INTEGRATION OVER ONE TETRAHEDRON                   **
!     **  AND THEIR DERIVATIVES WITH RESPECT TO THE ENERGIES AT THE CORNERS   **
!     **                                                                      **
!     **  CORRECTION FOR THE NONLINEAR SHAPE INCLUDED IF ICOR=1               **
!     **                                                                      **
!     **  AUTHOR : P.BLOECHL                                                  **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)  ,INTENT(IN) :: VOL      ! WEIGHT OF THIS TETRAHEDRON
      REAL(8)  ,INTENT(IN) :: EF       ! FERMI LEVEL
      REAL(8)  ,INTENT(IN) :: E_(4)    ! ENERGY BANDS AT THE CORNERS
      REAL(8)  ,INTENT(OUT):: WGHT(4)  ! INTEGRATION WEIGHTS
      REAL(8)  ,INTENT(OUT):: DWGHT(4,4)! DERIVATIVE OFINTEGRATION WEIGHTS
                                       ! WITH RESPECT TO ENERGIES
      INTEGER(4),PARAMETER :: ICOR=1   ! ON/OFF SWITCH FOR CORRECTION 
      REAL(8)              :: E(4)     ! ENERGY BANDS AT THE CORNERS
      REAL(8)              :: FA(4),DFA(4,4)
      REAL(8)              :: FB(4)
      INTEGER(4)           :: INDEX(4)
      REAL(8)              :: X
      REAL(8)              :: WGHTFAC(4),DWGHTFAC(4,4)
      REAL(8)              :: VPRIME,DVPRIME(4)
      REAL(8)              :: DE
      REAL(8)              :: DOS,DDOS(4)
      REAL(8)              :: E21,E31,E41,E32,E42,E43
      REAL(8)              :: DE1,DE2,DE3,DE4
      INTEGER(4)           :: I,J,IP,N,M,K
      REAL(8)              :: DA,DB,DC
      REAL(8)              :: VOL14
!     **************************************************************************
      E(:)=E_(:)
      WGHT(:)=0.D0
!
!     ==========================================================================
!     ==  INTEGRATION WITHOUT FERMISURFACE                            ==========
!     ==========================================================================
      X=MIN(E(1),E(2),E(3),E(4))                                      
      IF(X.GE.EF) RETURN
      X=MAX(E(1),E(2),E(3),E(4))                                      
      IF(X.LE.EF) THEN                                                  
        VPRIME=.25D0*VOL                                                
        DO I=1,4                                                     
          WGHT(I)=VPRIME                                                  
        ENDDO
        DWGHT(:,:)=0.D0
        RETURN                                                          
      END IF                                                            
!
!     ==========================================================================
!     ==  ORDER ENERGIES                                                      ==
!     ==========================================================================
!     -- INDEX HOLDS THE ORIGINAL POSITION OF THE ENERGIES AND WEIGHTS  
      DO I=1,4                                                      
        INDEX(I)=I                                                        
      ENDDO
      DO I=1,3                                                      
        IP=I                                                              
        DO J=I+1,4                                                    
          IF(E(IP).GT.E(J)) IP=J                                            
        ENDDO
        IF(IP.GT.I) THEN                                                  
          X=E(IP)                                                         
          E(IP)=E(I)                                                      
          E(I)=X                                                          
          K=INDEX(IP)                                                     
          INDEX(IP)=INDEX(I)                                              
          INDEX(I)=K                                                      
        END IF                                                            
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE UNCORRECTED INTEGRAL AS MEANVALUE                         ==
!     ==========================================================================
      E21=E(2)-E(1)                                                     
      E31=E(3)-E(1)                                                     
      E41=E(4)-E(1)                                                     
      E32=E(3)-E(2)                                                     
      E42=E(4)-E(2)                                                     
      E43=E(4)-E(3)                                                     
      IF(EF.GT.E(1).AND.EF.LE.E(2)) THEN                                
        DE=EF-E(1)                                                      
!
        VPRIME=.25D0*VOL*DE**3/(E21*E31*E41)                            
        DVPRIME(1)=.25D0*VOL*(-3.D0*DE**2)/(E21*E31*E41) &
    &              +VPRIME*(1.D0/E21+1.D0/E31+1.D0/E41)
        DVPRIME(2)=-VPRIME/E21
        DVPRIME(3)=-VPRIME/E31
        DVPRIME(4)=-VPRIME/E41
!
        WGHTFAC(1)=(4.D0-DE/E21-DE/E31-DE/E41)
        DWGHTFAC(1,1)=(1.D0-DE/E21)/E21+(1.D0-DE/E31)/E31+(1.D0-DE/E41)/E41
        DWGHTFAC(1,2)=+DE/E21**2
        DWGHTFAC(1,3)=+DE/E31**2
        DWGHTFAC(1,4)=+DE/E41**2
!
        WGHTFAC(2)=DE/E21                                           
        DWGHTFAC(2,1)=-(1.D0-DE/E21)/E21
        DWGHTFAC(2,2)=-DE/E21**2
        DWGHTFAC(2,3)=0.D0
        DWGHTFAC(2,4)=0.D0
!
        WGHTFAC(3)=DE/E31                                           
        DWGHTFAC(3,1)=-(1.D0-DE/E31)/E31
        DWGHTFAC(3,2)=0.D0
        DWGHTFAC(3,3)=-DE/E31**2
        DWGHTFAC(3,4)=0.D0
!
        WGHTFAC(4)=DE/E41                                           
        DWGHTFAC(4,1)=-(1.D0-DE/E41)/E41
        DWGHTFAC(4,2)=0.D0
        DWGHTFAC(4,3)=0.D0
        DWGHTFAC(4,4)=-DE/E41**2
!
        WGHT(:)=VPRIME*WGHTFAC(:)
        DO I=1,4        
          DWGHT(:,I)=DVPRIME(I)*WGHTFAC(:)+VPRIME*DWGHTFAC(:,I)
        ENDDO
!       ------  PARAMETERS FOR CORRECION                                
        DOS=12.D0*VPRIME/DE
        DDOS(:)=12.D0*DVPRIME(:)/DE
        DDOS(1)=DDOS(1)+DOS/DE
!
!     ==========================================================================
!     == SECOND CASE E1,E2<EF<E3,E4                                           ==
!     ==========================================================================
      ELSE IF(EF.GT.E(2).AND.EF.LT.E(3)) THEN                           
        DE1=EF-E(1)                                                     
        DE2=EF-E(2)                                                     
        DE3=E(3)-EF                                                     
        DE4=E(4)-EF                                                     
!       ------  TETRAHEDRON X1,X2,X13P,X14P                             
        VPRIME=0.25D0*VOL*DE1**2/(E41*E31)
        DVPRIME(1)=-0.5D0*VOL*DE1/(E41*E31)+VPRIME*(1.D0/E41+1.D0/E31)
        DVPRIME(2)=0.D0
        DVPRIME(3)=-VPRIME/E31
        DVPRIME(4)=-VPRIME/E41
!
        WGHTFAC(1)=3.D0-DE1/E41-DE1/E31
        DWGHTFAC(1,1)=(1.D0-DE1/E41)/E41+(1.D0-DE1/E31)/E31
        DWGHTFAC(1,2)=0.D0
        DWGHTFAC(1,3)=+DE1/E31**2
        DWGHTFAC(1,4)=+DE1/E41**2
!
        WGHTFAC(2)=1.D0
        DWGHTFAC(2,:)=0.D0
!
        WGHTFAC(3)=DE1/E31
        DWGHTFAC(3,1)=-(1.D0-DE1/E31)/E31
        DWGHTFAC(3,2)=0.D0
        DWGHTFAC(3,3)=-DE1/E31**2
        DWGHTFAC(3,4)=0.D0

        WGHTFAC(4)=DE1/E41
        DWGHTFAC(4,1)=-(1.D0-DE1/E41)/E41
        DWGHTFAC(4,2)=0.D0
        DWGHTFAC(4,3)=0.D0
        DWGHTFAC(4,4)=-DE1/E41**2
!
        WGHT(:)=VPRIME*WGHTFAC(:)
        DO I=1,4        
          DWGHT(:,I)=DVPRIME(I)*WGHTFAC(:)+VPRIME*DWGHTFAC(:,I)
        ENDDO
!       ------  TETRAHEDRON X2,X13P,X23P,X14P                           
        VPRIME=.25D0*VOL*DE2*DE3*DE1/(E32*E31*E41)                      
        DVPRIME(1)=VPRIME*(-1.D0/DE1+1.D0/E31+1.D0/E41)
        DVPRIME(2)=VPRIME*(-1.D0/DE2+1.D0/E32)
        DVPRIME(3)=VPRIME*(1.D0/DE3-1.D0/E32-1.D0/E31)
        DVPRIME(4)=VPRIME*(-1.D0/E41)
!
        WGHTFAC(1)=2.D0-DE1/E31-DE1/E41
        DWGHTFAC(1,1)=(1.D0-DE1/E31)/E31+(1.D0-DE1/E41)/E41
        DWGHTFAC(1,2)=0.D0
        DWGHTFAC(1,3)=DE1/E31**2
        DWGHTFAC(1,4)=DE1/E41**2
!
        WGHTFAC(2)=2.D0-DE2/E32
        DWGHTFAC(2,1)=0.D0
        DWGHTFAC(2,2)=(1.D0-DE2/E32)/E32
        DWGHTFAC(2,3)=DE2/E32**2
        DWGHTFAC(2,4)=0.D0
!
        WGHTFAC(3)=DE2/E32+DE1/E31
        DWGHTFAC(3,1)=-(1.D0-DE1/E31)/E31
        DWGHTFAC(3,2)=-(1.D0-DE2/E32)/E32
        DWGHTFAC(3,3)=-DE2/E32**2-DE1/E31**2
        DWGHTFAC(3,4)=0.D0
!
        WGHTFAC(4)=DE1/E41
        DWGHTFAC(4,1)=-(1.D0-DE1/E41)/E41
        DWGHTFAC(4,2)=0.D0
        DWGHTFAC(4,3)=0.D0
        DWGHTFAC(4,4)=-DE1/E41**2
!
        WGHT(:)=WGHT(:)+VPRIME*WGHTFAC(:) 
        DO I=1,4        
          DWGHT(:,I)=DWGHT(:,I)+DVPRIME(I)*WGHTFAC(:)+VPRIME*DWGHTFAC(:,I)
        ENDDO
!
!       ------  TETRAHEDRON X2,X23P,X24P,X14P  ---------------------------------
        VPRIME=.25D0*VOL*DE2**2*DE4/(E42*E32*E41)                       
        DVPRIME(1)=VPRIME/E41
        DVPRIME(2)=VPRIME*(-2.D0/DE2+1.D0/E42+1.D0/E32)
        DVPRIME(3)=VPRIME*(-1.D0/E32)
        DVPRIME(4)=VPRIME*(1.D0/DE4-1.D0/E42-1.D0/E41)

        WGHTFAC(1)=1.D0-DE1/E41
        DWGHTFAC(1,1)=(1.D0-DE1/E41)/E41
        DWGHTFAC(1,2)=0.D0
        DWGHTFAC(1,3)=0.D0 
        DWGHTFAC(1,4)=DE1/E41**2
!
        WGHTFAC(2)=3.D0-DE2/E32-DE2/E42
        DWGHTFAC(2,1)=0.D0
        DWGHTFAC(2,2)=(1.D0-DE2/E32)/E32+(1.D0-DE2/E42)/E42
        DWGHTFAC(2,3)=DE2/E32**2
        DWGHTFAC(2,4)=DE2/E42**2

        WGHTFAC(3)=DE2/E32
        DWGHTFAC(3,1)=0.D0
        DWGHTFAC(3,2)=-(1.D0-DE2/E32)/E32
        DWGHTFAC(3,3)=-DE2/E32**2
        DWGHTFAC(3,4)=0.D0
!
        WGHTFAC(4)=DE2/E42+DE1/E41
        DWGHTFAC(4,1)=-(1.D0-DE1/E41)/E41
        DWGHTFAC(4,2)=-(1.D0-DE2/E42)/E42
        DWGHTFAC(4,3)=0.D0
        DWGHTFAC(4,4)=-DE2/E42**2-DE1/E41**2
!
        WGHT(:)=WGHT(:)+VPRIME*WGHTFAC(:) 
        DO I=1,4        
          DWGHT(:,I)=DWGHT(:,I)+DVPRIME(I)*WGHTFAC(:)+VPRIME*DWGHTFAC(:,I)
        ENDDO
!
!       ------  DOS=A+B*(EF-E2)+C*(EF-E2)**2                            
        DA=3.D0*VOL*E21/(E31*E41)                                       
        DB=6.D0*VOL/(E31*E41)                                           
        DC=-3.D0*VOL/(E32*E41*E31*E42)*(E31+E42)                        
        DOS=DA+DB*DE2+DC*DE2**2                                         
        DDOS(1)=DA*(-1.D0/E21+1.D0/E31+1.D0/E41) &
       &       +DB*(1.D0/E31+1.D0/E41)*DE2 &
       &       +DC*(1.D0/E41+1.D0/E31-1.D0/(E31+E42))*DE2**2
        DDOS(2)=DA*(1.D0/E21) &
       &       +DC*(1.D0/E32+1.D0/E42-1.D0/(E31+E42))*DE2**2 &
       &       -DB-2.D0*DC*DE2
        DDOS(3)=DA*(-1.D0/E31) &
       &       +DB*(-1.D0/E31)*DE2 &
       &       +DC*(-1.D0/E32-1.D0/E31+1.D0/(E31+E42))*DE2**2
        DDOS(4)=DA*(-1.D0/E41) &
       &       +DB*(-1.D0/E41)*DE2 &
       &       +DC*(-1.D0/E41-1.D0/E42+1.D0/(E31+E42))*DE2**2 
!
!     ==========================================================================
!     == LAST CASE                                                            ==
!     ==========================================================================
      ELSE IF(EF.GE.E(3).AND.EF.LT.E(4)) THEN                           
        DE=E(4)-EF                                                      
        VPRIME=.25D0*VOL*DE**3/(E41*E42*E43)                            
        DVPRIME(1)=VPRIME*(1.D0/E41)
        DVPRIME(2)=VPRIME*(1.D0/E42)
        DVPRIME(3)=VPRIME*(1.D0/E43)
        DVPRIME(4)=VPRIME*3.D0/DE-VPRIME*(1.D0/E41+1.D0/E42+1.D0/E43)
!
        WGHTFAC(1)=-DE/E41
        DWGHTFAC(1,1)=-DE/E41**2
        DWGHTFAC(1,2)=0.D0
        DWGHTFAC(1,3)=0.D0
        DWGHTFAC(1,4)=-(1.D0-DE/E41)/E41
!
        WGHTFAC(2)=-DE/E42
        DWGHTFAC(2,1)=0.D0
        DWGHTFAC(2,2)=-DE/E42**2
        DWGHTFAC(2,3)=0.D0
        DWGHTFAC(2,4)=-(1.D0-DE/E42)/E42

        WGHTFAC(3)=-DE/E43
        DWGHTFAC(3,1)=0.D0
        DWGHTFAC(3,2)=0.D0
        DWGHTFAC(3,3)=-DE/E43**2
        DWGHTFAC(3,4)=-(1.D0-DE/E43)/E43

        WGHTFAC(4)=-(4.D0-DE/E41-DE/E42-DE/E43)
        DWGHTFAC(4,1)=DE/E41**2
        DWGHTFAC(4,2)=DE/E42**2
        DWGHTFAC(4,3)=DE/E43**2
        DWGHTFAC(4,4)=(1.D0-DE/E41)/E41+(1.D0-DE/E42)/E42+(1.D0-DE/E43)/E43
        VOL14=.25D0*VOL                                                 
        WGHT(:)=VOL14+VPRIME*WGHTFAC(:) 
        DO I=1,4        
          DWGHT(:,I)=DVPRIME(I)*WGHTFAC(:)+VPRIME*DWGHTFAC(:,I)
        ENDDO
!       ------  PARAMETERS FOR CORRECION                                
        DOS=12.D0*VPRIME/DE                                
        DDOS(1)=12.D0*DVPRIME(1)/DE
        DDOS(2)=12.D0*DVPRIME(2)/DE
        DDOS(3)=12.D0*DVPRIME(3)/DE
        DDOS(4)=12.D0*DVPRIME(4)/DE-DOS/DE
      ELSE                                                              
        CALL ERROR$MSG('ERROR')
        CALL ERROR$STOP('BRILLOUIN_WEIGHT')
      END IF                                                            
!
!     ==========================================================================
!     ==  ADD CORRECTION FOR QUADRATIC DEVIATION                              ==
!     ==========================================================================
      IF(ICOR.EQ.1) THEN                                                
        DO M=1,4                                                    
          DO N=1,4                                                    
            WGHT(M)=WGHT(M)+.25D0*(E(N)-E(M))*DOS*.1D0                      
            DWGHT(M,:)=DWGHT(M,:)+0.25*(E(N)-E(M))*DDOS(:)*0.1D0
            DWGHT(M,N)=DWGHT(M,N)+0.25*DOS*0.1D0
            DWGHT(M,M)=DWGHT(M,M)-0.25*DOS*0.1D0
          ENDDO
        ENDDO
      END IF                                                            
!
!     ==========================================================================
!     ==  REORDER WEIGHTS                                                     ==
!     ==========================================================================
      DO I=1,4                                                      
        FA(INDEX(I))=WGHT(I)                                              
        FB(INDEX(I))=E(I)                                                 
        DO J=1,4
          DFA(INDEX(I),INDEX(J))=DWGHT(I,J)
        ENDDO
      ENDDO
      DO I=1,4                                                      
        WGHT(I)=FA(I)                                                     
        E(I)=FB(I)                                                        
        DO J=1,4
          DWGHT(I,J)=DFA(I,J)
        ENDDO
      ENDDO
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$WGHT(NKP,NB,EF,EB,WGHT)
!     **************************************************************************
!     ** CALCULATE BRILLOUIN INTEGRATION WEIGHTS FOR SPECIFIED BANDS AND      **
!     ** A GIVEN FERMI LEVEL                                                  **
!     **                                                                      **
!     ** IS BASICALLY IDENTICAL TO BRILLOUIN_SAMFAC                           **
!     **************************************************************************
      USE BRILLOUIN_MODULE, ONLY : THIS
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NKP
      INTEGER(4),INTENT(IN)  :: NB
      REAL(8)   ,INTENT(IN)  :: EF
      REAL(8)   ,INTENT(IN)  :: EB(NB,NKP)
      REAL(8)   ,INTENT(OUT) :: WGHT(NB,NKP)
      REAL(8)                :: E(4)
      REAL(8)                :: WGHT0(4)
      REAL(8)                :: VOL
      INTEGER(4)             :: IKP(4)
      INTEGER(4)             :: ITET,NTET
      INTEGER(4)             :: IB,I
!     **************************************************************************
                                        CALL TRACE$PUSH('BRILLOUIN$WGHT')
      WGHT(:,:)=0.D0
      NTET=THIS%NTET
      DO ITET=1,NTET                                                
        VOL=THIS%VOL*THIS%MULT(ITET)
        IKP(:)=THIS%IKP(:,ITET)
        DO IB=1,NB                                                    
          DO I=1,4                                                      
            E(I)=EB(IB,IKP(I))                                                
            WGHT0(I)=0.D0
          ENDDO
          CALL BRILLOUIN_WEIGHT(VOL,E,EF,WGHT0)
          DO I=1,4                                                      
            WGHT(IB,IKP(I))=WGHT(IB,IKP(I))+WGHT0(I)                          
          ENDDO
        ENDDO
      ENDDO
                                        CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$EWGHT(NKP,NB,EB,EMIN,EMAX,NE,EWGHT)
!     **************************************************************************
!     ** GENERATES ENERGY DEPENDENT WEIGHTS FOR DOS AND PDOS                  **
!     **                                                                      **
!     ** THE DENSITY OF STATES FOR A MATRIX ELEMENT A IS GIVEN                **
!     **      D_A(E)=SUM_N SUM_K V_N(K,E)*A_N(K)                              **
!     ** WHERE V_N(K)=EWGHT(N,K)%WGHT(E) ARE THE ENERGY DERIVATIVES           **
!     ** OF THE INTEGRATION WEIGHTS.                                          **
!     **                                                                      **
!     ** THE NUMBER OF STATES IS OBTAINED FROM WEIGHTS AT THE MINIMUM ENERGY  **
!     ** AND THE ENERGY INTEGRAL OF THE DENSITY OF STATES                     **
!     **     N_A(E)=SUM_N SUM_K W_N(K,E)*A_N(K)                               **
!     ** WHERE W_N(K,E)=W_N(K,EMIN)+SUM_{E<EMIN} V_N(K,E)*DE                  **
!     ** AND W_N(K,EMIN) FOR THE MINIMUM ENERGY IS CALCULATED WITH            **
!     ** BRILLOUIN_SAMFAC                                                     **
!     **************************************************************************
      USE BRILLOUIN_MODULE, ONLY : THIS,EWGHT_TYPE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NKP
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: NE
      REAL(8)   ,INTENT(IN) :: EB(NB,NKP)
      REAL(8)   ,INTENT(IN) :: EMIN
      REAL(8)   ,INTENT(IN) :: EMAX
      TYPE(EWGHT_TYPE),INTENT(INOUT) :: EWGHT(NB,NKP)
      INTEGER(4)            :: IKP(4)
      REAL(8)               :: VOL
      REAL(8)               :: E(4)
      REAL(8)               :: EF
      REAL(8)               :: WGHT0(4)
      REAL(8)               :: AE,BE,DE
      REAL(8)               :: SVARM,SVARP
      INTEGER(4)            :: NTET  !#(TETRAHEDRA)
      INTEGER(4)            :: ITET,IK,IB,I,IE
      INTEGER(4)            :: I1,I2
!     **************************************************************************
      NTET=THIS%NTET
!
!     ==========================================================================
!     == PARAMETERS FOR THE ENERGY GRID                                       ==
!     ==  E=EMIN+(EMAX-EMIN)/(NE-1)*(IE-1)  ====================================
!     ==  IE=1+(NE-1)*(E-EMIN)/(EMAX-EMIN) =====================================
!     ==  IE=[(NE-1)/(EMAX-EMIN)] * E + [1-(NE-1)*EMIN/(EMAX-EMIN)]
!     ==  IE=AE * E + BE
!     ==========================================================================
      AE=REAL(NE-1,KIND=8)/(EMAX-EMIN)
      BE=1.D0-REAL(NE-1,KIND=8)*EMIN/(EMAX-EMIN)
      DE=(EMAX-EMIN)/REAL(NE-1,KIND=8)
!
!     ==========================================================================
!     == DETERMINE GRID RANGE FOR EACH BAND AND K-POINT                       ==
!     ==========================================================================
      EWGHT(:,:)%I1=NE
      EWGHT(:,:)%I2=1
      DO ITET=1,NTET
        IKP(:)=THIS%IKP(:,ITET)
        DO IB=1,NB
!         == WGHT(I2) CONTAINS THE COMPLETE INTEGRAL ===========================
          I1=INT(AE*MINVAL(EB(IB,IKP(:)))+BE)
          I2=1+INT(AE*MAXVAL(EB(IB,IKP(:)))+BE)
          EWGHT(IB,IKP)%I1=MIN(EWGHT(IB,IKP)%I1,I1)
          EWGHT(IB,IKP)%I2=MAX(EWGHT(IB,IKP)%I2,I2)
        ENDDO
      ENDDO
      DO IK=1,NKP
        DO IB=1,NB
          I1=EWGHT(IB,IK)%I1
          I2=EWGHT(IB,IK)%I2
          I1=MAX(I1,1)
          I2=MIN(I2,NE)
          EWGHT(IB,IK)%I1=I1
          EWGHT(IB,IK)%I2=I2
          ALLOCATE(EWGHT(IB,IK)%WGHT(I1:I2))
          EWGHT(IB,IK)%WGHT(:)=0.D0
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == ADD UP INTEGRATION WEIGHTS                                           ==
!     ==========================================================================
      DO ITET=1,NTET
        VOL=THIS%VOL*THIS%MULT(ITET)
        IKP(:)=THIS%IKP(:,ITET)
        DO IB=1,NB
          I1=NE
          I2=1
          DO I=1,4                                                      
            E(I)=EB(IB,IKP(I))                                                
            I1=MIN(I1,1+INT(AE*E(I)+BE))
            I2=MAX(I2,INT(AE*E(I)+BE))
          ENDDO
          I1=MAX(I1,1)
          I2=MIN(I2,NE)
          DO IE=I1,I2          
            EF=EMIN+REAL(IE-1,KIND=8)*DE
            WGHT0(:)=0.D0
            CALL BRILLOUIN_WEIGHT(VOL,E,EF,WGHT0)
            DO I=1,4                                                      
              EWGHT(IB,IKP(I))%WGHT(IE)=EWGHT(IB,IKP(I))%WGHT(IE)+WGHT0(I)
            ENDDO
          ENDDO
          DO I=1,4                                                      
            EWGHT(IB,IKP(I))%WGHT(I2+1:)=EWGHT(IB,IKP(I))%WGHT(I2+1:)+0.25D0*VOL
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DIFFERENTIATE INTEGRATION WEIGHTS                                    ==
!     ==========================================================================
      DO IK=1,NKP
        DO IB=1,NB
          I1=EWGHT(IB,IK)%I1
          I2=EWGHT(IB,IK)%I2
          SVARM=0.5D0*EWGHT(IB,IK)%WGHT(I1)/DE
          DO IE=I1,I2-1
            SVARP=0.5D0*(EWGHT(IB,IK)%WGHT(IE+1)-EWGHT(IB,IK)%WGHT(IE))/DE
            EWGHT(IB,IK)%WGHT(IE)=SVARM+SVARP
            SVARM=SVARP
          ENDDO
          EWGHT(IB,IK)%WGHT(I2)=SVARM
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE BRILLOUIN$EWGHT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_COMPLETE(NOPX,NOP,OP)
!     **************************************************************************
!     ** PRODUCES THE COMPLETE SET OF POINT GROUP OPERATIONS FROM A           **
!     ** SET OF GENERATORS OF THE GROUP.                                      **
!     **                                                                      **
!     ** REMARK: THE ROUTINE TOLERATES IF MORE THAN A MINIMUM SET OF          **
!     **         GENERATORS IS PROVIDED                                       **
!     ** REMARK: THE IDENTITY NEED NOT BE INCLUDED IN THE SET OF GENERATORS   **
!     **         GENERATORS IS PROVIDED                                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NOPX ! MAX #(POINT GROUP OPERATIONS)
      INTEGER(4),INTENT(INOUT) :: NOP  ! IN: #(GENERATORS) OUT:  #(OPERATIONS)
      INTEGER(4),INTENT(INOUT) :: OP(3,3,NOPX) ! IN: GENERATORS OUT: OPERATIONS
      INTEGER(4)               :: GENERATOR(3,3,NOPX) ! COPY OF THE GENERATORS
      INTEGER(4)               :: OPNEW(3,3)
      INTEGER(4)               :: NGENERATOR
      INTEGER(4)               :: I,IOP,IOP2
      LOGICAL(4)               :: TOLD
!     **************************************************************************
      NGENERATOR=NOP
      GENERATOR(:,:,1:NOP)=OP(:,:,1:NOP)
!
!     == IDENTITY ==============================================================
      OP(:,:,:)=0
      OP(1,1,1)=1
      OP(2,2,1)=1
      OP(3,3,1)=1
      NOP=1
!
!     ==========================================================================
!     == COMPLETE SYMMETRY GROUP                                              ==
!     == ENSURES THAT NO NEW OPERATION IS OBTAINED BY APPLYING ANY OF THE     ==
!     == GENERATORS TO ANT OF THE OPERATIONS IN THE GROUP                     ==
!     ==========================================================================
      NOP=1
      IOP=1
      DO WHILE(IOP.LE.NOP) !LOOP OVER OPERATIONS
!WRITE(*,FMT='("OLD OPERATION ",I5,T20,3("|",3I5,"|"))')IOP,OP(:,:,IOP)
        DO I=1,NGENERATOR
          OPNEW=MATMUL(GENERATOR(:,:,I),OP(:,:,IOP))
!WRITE(*,FMT='("OPNEW ",T20,3("|",3I5,"|"))') OPNEW
!         == CHECK IF OPNEW ALREADY PRESENT IN OP ============================
          TOLD=.FALSE.
          DO IOP2=1,NOP
            TOLD=SUM(ABS(OPNEW-OP(:,:,IOP2))).EQ.0
            IF(TOLD) EXIT
          ENDDO
          IF(TOLD) CYCLE
          NOP=NOP+1
          IF(NOP.GT.NOPX) THEN
            CALL ERROR$MSG('NUMBER OF OPERATIONS EXCEEDS MAXIMUM')
            CALL ERROR$I4VAL('NOP',NOP)
            CALL ERROR$I4VAL('NOPX',NOPX)
            CALL ERROR$STOP('BRILLOUIN_COMPLETE')
          END IF
          OP(:,:,NOP)=OPNEW
!WRITE(*,FMT='("NEW OP: ",I5,T20,3("|",3I5,"|"))')NOP,OP(:,:,NOP)
        ENDDO ! END OF LOOP OVER GENERATORS
        IOP=IOP+1
      ENDDO !END OF LOOP OVER OPERATIONS
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE SPACEGROUP_MODULE
!*******************************************************************************
!
!** THE TRANSLATION T0 OF SG_TYPE IS THE TRANSLATION OF THE ORIGIN THAT,      **
!** SUBSTITUTED IN EQ.32.5.11 OF BRADLEY+CRACKNELL CHANGES THE                **
!** ELEMENTS {R|V} THAT WE USE INTO THE ELEMENTS {R'|V'} OF THE               **
!** INTERNATIONAL TABLES OF CRYSTALLOGRAPHY. (SEE FOOTNOTE V OF TABLE         **
!** 3.7 OF BRADLEY CRACKNELL.)                                                **
!*******************************************************************************
TYPE SG_TYPE  !SPACE-GROUP TYPE
  CHARACTER(3) :: BRAVAIS
  CHARACTER(16):: INTERNATIONALSYMBOL
  REAL(8)      :: T0(3)   ! T0 IN TABLE 3.7 OF BRADLEY CRACKNELL
END TYPE SG_TYPE
LOGICAL(4)    :: TINI=.FALSE.
INTEGER(4)    :: ISPACEGROUP=0
TYPE(SG_TYPE) :: SG(230)
END MODULE SPACEGROUP_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPACEGROUP_INI()
!     **************************************************************************
!

!     **************************************************************************
      USE SPACEGROUP_MODULE, ONLY : TINI,SG,SG_TYPE
      IMPLICIT NONE
!     **************************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
      SG(  1)=SG_TYPE('G1','P1',(/0.D0,0.D0,0.D0/))
      SG(  2)=SG_TYPE('G1','P\BAR{1}',(/0.D0,0.D0,0.D0/))
      SG(  3)=SG_TYPE('GM','P2',(/0.D0,0.D0,0.D0/))  
      SG(  4)=SG_TYPE('GM','P2_1',(/0.D0,0.D0,0.D0/))
      SG(  5)=SG_TYPE('GMB','B2',(/0.D0,0.D0,0.D0/))
      SG(  6)=SG_TYPE('GM','PM',(/0.D0,0.D0,0.D0/))
      SG(  7)=SG_TYPE('GM','PB',(/0.D0,0.D0,0.D0/))
      SG(  8)=SG_TYPE('GMB','BM',(/0.D0,0.D0,0.D0/))
      SG(  9)=SG_TYPE('GMB','BP',(/0.D0,0.D0,0.D0/))  
      SG( 10)=SG_TYPE('GM','P2/M',(/0.D0,0.D0,0.D0/))
      SG( 11)=SG_TYPE('GM','P2_1/M',(/0.D0,0.D0,0.25D0/))
      SG( 12)=SG_TYPE('GBM','B2/M',(/0.D0,0.D0,0.D0/))
      SG( 13)=SG_TYPE('GM','P2/B',(/-0.25D0,0.D0,0.25D0/))
      SG( 14)=SG_TYPE('GM','P2_1/B',(/0.D0,0.D0,0.D0/))
      SG( 15)=SG_TYPE('GMB','B2/B',(/0.25D0,0.D0,0.D0/))  
      SG( 16)=SG_TYPE('GO','P222',(/0.D0,0.D0,0.D0/))
      SG( 17)=SG_TYPE('GO','P222_1',(/0.D0,0.D0,0.D0/))
      SG( 18)=SG_TYPE('GO','P2_12_12',(/0.D0,0.D0,0.D0/))
      SG( 19)=SG_TYPE('GO','P2_12_12_1',(/0.D0,0.D0,0.D0/))
      SG( 20)=SG_TYPE('GOB','C222_1',(/0.D0,0.D0,0.D0/))
      SG( 21)=SG_TYPE('GOB','C222',(/0.D0,0.D0,0.D0/))
      SG( 22)=SG_TYPE('GOF','F222',(/0.D0,0.D0,0.D0/))
      SG( 23)=SG_TYPE('GOV','I222',(/0.D0,0.D0,0.D0/))
      SG( 24)=SG_TYPE('GOV','I2_12_12_1',(/0.5D0,0.D0,0.D0/))
      SG( 25)=SG_TYPE('GO','PMM2',(/0.D0,0.D0,0.D0/))
      SG( 26)=SG_TYPE('GO','PMC2_1',(/0.D0,0.D0,0.D0/))
      SG( 27)=SG_TYPE('GO','PCC2',(/0.D0,0.D0,0.D0/))
      SG( 28)=SG_TYPE('GO','PMA2',(/-0.25D0,0.D0,0.D0/))
      SG( 29)=SG_TYPE('GO','PCA2_1',(/-0.25D0,0.D0,0.D0/))
      SG( 30)=SG_TYPE('GO','PNC2',(/-0.25D0,0.D0,0.D0/))
      SG( 31)=SG_TYPE('GO','PMN2_1',(/0.D0,0.D0,0.D0/))
      SG( 32)=SG_TYPE('GO','PBA2',(/-0.25D0,0.25D0,0.D0/))
      SG( 33)=SG_TYPE('GO','PNA2_1',(/-0.25D0,0.25D0,0.D0/))
      SG( 34)=SG_TYPE('GO','PNN2',(/-0.25D0,0.25D0,0.D0/))
      SG( 35)=SG_TYPE('GOB','CMM2',(/0.D0,0.D0,0.D0/))
      SG( 36)=SG_TYPE('GOB','CMC2_1',(/0.D0,0.D0,0.D0/))
      SG( 37)=SG_TYPE('GOB','CCC2',(/0.D0,0.D0,0.D0/))
      SG( 38)=SG_TYPE('GOB','AMM2',(/0.D0,0.D0,0.D0/))
      SG( 39)=SG_TYPE('GOB','ABM2',(/0.D0,0.D0,0.D0/))
      SG( 40)=SG_TYPE('GOB','AMA2',(/0.D0,0.D0,0.D0/))
      SG( 41)=SG_TYPE('GOB','ABA2',(/0.D0,0.D0,0.D0/))
      SG( 42)=SG_TYPE('GOF','FMM2',(/0.D0,0.D0,0.D0/))
      SG( 43)=SG_TYPE('GOF','FDD2',(/0.D0,0.D0,0.25D0/))
      SG( 44)=SG_TYPE('GOV','IMM2',(/0.D0,0.D0,0.D0/))
      SG( 45)=SG_TYPE('GOV','IBA2',(/0.D0,0.D0,0.D0/))
      SG( 46)=SG_TYPE('GOV','IMA2',(/0.D0,0.D0,0.D0/))
      SG( 47)=SG_TYPE('GO','PMMM',(/0.D0,0.D0,0.D0/))
      SG( 48)=SG_TYPE('GO','PNNN',(/0.D0,0.D0,0.D0/))
      SG( 49)=SG_TYPE('GO','PCCM',(/0.D0,0.D0,0.25D0/))
      SG( 50)=SG_TYPE('GO','PBAN',(/0.D0,0.D0,0.D0/))
      SG( 51)=SG_TYPE('GO','PMMA',(/0.D0,0.D0,0.D0/))
      SG( 52)=SG_TYPE('GO','PNNA',(/0.25D0,0.25D0,0.D0/))
      SG( 53)=SG_TYPE('GO','PMNA',(/-0.25D0,0.D0,0.D0/))
      SG( 54)=SG_TYPE('GO','PCCA',(/0.D0,0.25D0,0.D0/))
      SG( 55)=SG_TYPE('GO','PBAM',(/0.D0,0.D0,0.D0/))
      SG( 56)=SG_TYPE('GO','PCCN',(/-0.25D0,0.25D0,0.25D0/))
      SG( 57)=SG_TYPE('GO','PBCM',(/0.D0,0.25D0,0.D0/))
      SG( 58)=SG_TYPE('GO','PNNM',(/0.D0,0.D0,-0.25D0/))
      SG( 59)=SG_TYPE('GO','PMMN',(/0.D0,0.D0,0.D0/))
      SG( 60)=SG_TYPE('GO','PBCN',(/-0.25D0,0.D0,0.25D0/))
      SG( 61)=SG_TYPE('GO','PBCA',(/0.D0,0.D0,0.D0/))
      SG( 62)=SG_TYPE('GO','PNMA',(/0.25D0,0.25D0,0.D0/))
      SG( 63)=SG_TYPE('GOB','CMCM',(/0.D0,0.D0,0.D0/))
      SG( 64)=SG_TYPE('GOB','CMCA',(/0.25D0,0.25D0,0.D0/))
      SG( 65)=SG_TYPE('GOB','CMMM',(/0.D0,0.D0,0.D0/))
      SG( 66)=SG_TYPE('GOB','CCCM',(/0.D0,0.D0,0.25D0/))
      SG( 67)=SG_TYPE('GOB','CMMA',(/0.25D0,0.25D0,0.D0/))
      SG( 68)=SG_TYPE('GOB','CCCA',(/-0.25D0,0.25D0,0.25D0/))
      SG( 69)=SG_TYPE('GOF','FMMM',(/0.D0,0.D0,0.D0/))
      SG( 70)=SG_TYPE('GOF','FDDD',(/0.D0,0.D0,0.D0/))
      SG( 71)=SG_TYPE('GOV','IMMM',(/0.D0,0.D0,0.D0/))
      SG( 72)=SG_TYPE('GOV','IBAM',(/0.25D0,0.25D0,0.D0/))
      SG( 73)=SG_TYPE('GOV','IBAC',(/0.25D0,0.D0,0.D0/))
      SG( 74)=SG_TYPE('GOV','IMMA',(/0.75D0,0.25D0,0.D0/))
      SG( 75)=SG_TYPE('GQ','P4',(/0.D0,0.D0,0.D0/))
      SG( 76)=SG_TYPE('GQ','P4_1',(/0.D0,0.D0,0.D0/))
      SG( 77)=SG_TYPE('GQ','P4_2',(/0.D0,0.D0,0.D0/))
      SG( 78)=SG_TYPE('GQ','P4_3',(/0.D0,0.D0,0.D0/))
      SG( 79)=SG_TYPE('GQV','I4',(/0.D0,0.D0,0.D0/))
      SG( 80)=SG_TYPE('GQV','I4_1',(/0.D0,0.D0,0.D0/))
      SG( 81)=SG_TYPE('GQ','P\BAR{4}',(/0.D0,0.D0,0.D0/))
      SG( 82)=SG_TYPE('GQV','I\BAR{4}',(/0.D0,0.D0,0.D0/))
      SG( 83)=SG_TYPE('GQ','P4/M',(/0.D0,0.D0,0.D0/))
      SG( 84)=SG_TYPE('GQ','P4_2/M',(/0.D0,0.D0,0.25D0/))
      SG( 85)=SG_TYPE('GQ','P4/N',(/0.D0,0.D0,0.D0/))
      SG( 86)=SG_TYPE('GQ','P4_2/N',(/0.D0,0.D0,0.D0/))
      SG( 87)=SG_TYPE('GQV','I4/M',(/0.D0,0.D0,0.D0/))
      SG( 88)=SG_TYPE('GQV','I4_1/M',(/0.D0,0.D0,0.D0/))
      SG( 89)=SG_TYPE('GQ','P422',(/0.D0,0.D0,0.D0/))
      SG( 90)=SG_TYPE('GQ','P42_12',(/0.5D0,0.D0,0.D0/))
      SG( 91)=SG_TYPE('GQ','P4_122',(/0.D0,0.D0,0.25D0/))
      SG( 92)=SG_TYPE('GQ','P4_12_12',(/0.5D0,0.D0,-0.375D0/))
      SG( 93)=SG_TYPE('GQ','P4_222',(/0.D0,0.D0,0.D0/))
      SG( 94)=SG_TYPE('GQ','P4_22_12',(/0.5D0,0.D0,0.25D0/))
      SG( 95)=SG_TYPE('GQ','P4_322',(/0.D0,0.D0,0.25D0/))
      SG( 96)=SG_TYPE('GQ','P4_32_12',(/0.5D0,0.D0,-0.125D0/))
      SG( 97)=SG_TYPE('GQV','I422',(/0.D0,0.D0,0.D0/))
      SG( 98)=SG_TYPE('GQV','I4_122',(/0.125D0,0.125D0,0.D0/))
      SG( 99)=SG_TYPE('GQ','P4MM',(/0.D0,0.D0,0.D0/))
      SG(100)=SG_TYPE('GQ','P4BM',(/0.D0,0.D0,0.D0/))
      SG(101)=SG_TYPE('GQ','P4_2CM',(/0.D0,0.D0,0.D0/))
      SG(102)=SG_TYPE('GQ','P4_2NM',(/0.5D0,0.D0,0.D0/))
      SG(103)=SG_TYPE('GQ','P4CC',(/0.D0,0.D0,0.D0/))
      SG(104)=SG_TYPE('GQ','P4NC',(/0.D0,0.D0,0.D0/))
      SG(105)=SG_TYPE('GQ','P4_2MC',(/0.D0,0.D0,0.D0/))
      SG(106)=SG_TYPE('GQ','P4_2BC',(/0.D0,0.D0,0.D0/))
      SG(107)=SG_TYPE('GQV','I4MM',(/0.D0,0.D0,0.D0/))
      SG(108)=SG_TYPE('GQV','I4CM',(/0.D0,0.D0,0.D0/))
      SG(109)=SG_TYPE('GQV','I4_1MD',(/0.D0,0.D0,0.D0/))
      SG(110)=SG_TYPE('GQV','I4_1CD',(/0.D0,0.D0,0.D0/))
      SG(111)=SG_TYPE('GQ','P\BAR{4}2M',(/0.D0,0.D0,0.D0/))
      SG(112)=SG_TYPE('GQ','P\BAR{4}2C',(/0.D0,0.D0,0.D0/))
      SG(113)=SG_TYPE('GQ','P\BAR{4}2_1M',(/0.D0,0.D0,0.D0/))
      SG(114)=SG_TYPE('GQ','P\BAR{4}2_1C',(/0.D0,0.D0,0.D0/))
      SG(115)=SG_TYPE('GQ','P\BAR{4}M2',(/0.D0,0.D0,0.D0/))
      SG(116)=SG_TYPE('GQ','P\BAR{4}C2',(/0.D0,0.D0,0.D0/)) 
      SG(117)=SG_TYPE('GQ','P\BAR{4}B2',(/0.D0,0.D0,0.D0/)) 
      SG(118)=SG_TYPE('GQ','P\BAR{4}N2',(/0.D0,0.D0,0.D0/))
      SG(119)=SG_TYPE('GQV','I\BAR{4}M2',(/0.D0,0.D0,0.D0/))
      SG(120)=SG_TYPE('GQV','I\BAR{4}C2',(/0.D0,0.D0,0.D0/))
      SG(121)=SG_TYPE('GQV','I\BAR{4}2M',(/0.D0,0.D0,0.D0/))
      SG(122)=SG_TYPE('GQV','I\BAR{4}2D',(/0.D0,0.D0,0.D0/))
      SG(123)=SG_TYPE('GQ','P4/MMM',(/0.D0,0.D0,0.D0/))
      SG(124)=SG_TYPE('GQ','P4/MCC',(/0.D0,0.D0,0.25D0/))
      SG(125)=SG_TYPE('GQ','P4/NBM',(/0.25D0,-0.25D0,0.D0/))
      SG(126)=SG_TYPE('GQ','P4/NNC',(/0.5D0,0.D0,0.D0/))
      SG(127)=SG_TYPE('GQ','P4/MBM',(/0.5D0,0.D0,0.D0/))
      SG(128)=SG_TYPE('GQ','P4/MNC',(/0.5D0,0.D0,0.25D0/))
      SG(129)=SG_TYPE('GQ','P4/NMM',(/0.5D0,0.D0,0.D0/))
      SG(130)=SG_TYPE('GQ','P4/NCC',(/0.5D0,0.D0,0.25D0/))
      SG(131)=SG_TYPE('GQ','P4_2/MMC',(/0.D0,0.D0,0.D0/))
      SG(132)=SG_TYPE('GQ','P4_2/MCM',(/0.D0,0.D0,0.25D0/))
      SG(133)=SG_TYPE('GQ','P4_2/NBC',(/0.D0,0.D0,0.25D0/))
      SG(134)=SG_TYPE('GQ','P4_2/NNM',(/0.D0,0.D0,0.D0/))
      SG(135)=SG_TYPE('GQ','P4_2/MBC',(/0.5D0,0.D0,0.D0/))
      SG(136)=SG_TYPE('GQ','P4_2/MNM',(/0.D0,0.D0,0.25D0/))
      SG(137)=SG_TYPE('GQ','P4_2/NMC',(/0.5D0,0.D0,0.25D0/))
      SG(138)=SG_TYPE('GQ','P4_2/NCM',(/0.5D0,0.D0,0.D0/))
      SG(139)=SG_TYPE('GQV','I4/MMM',(/0.D0,0.D0,0.D0/))
      SG(140)=SG_TYPE('GQV','I4/MCM',(/0.75D0,0.25D0,0.5D0/))
      SG(141)=SG_TYPE('GQV','I4_1/AMD',(/0.25D0,0.25D0,0.D0/))
      SG(142)=SG_TYPE('GQV','I4_1/ACD',(/0.D0,0.D0,0.D0/))
      SG(143)=SG_TYPE('GH','P3',(/0.D0,0.D0,0.D0/))
      SG(144)=SG_TYPE('GH','P3_1',(/0.D0,0.D0,0.D0/))
      SG(145)=SG_TYPE('GH','P3_2',(/0.D0,0.D0,0.D0/))
      SG(146)=SG_TYPE('GRH','R3',(/0.D0,0.D0,0.D0/))
      SG(147)=SG_TYPE('GH','P\BAR{3}',(/0.D0,0.D0,0.D0/)) 
      SG(148)=SG_TYPE('GRH','R\BAR{3}',(/0.D0,0.D0,0.D0/)) 
      SG(149)=SG_TYPE('GH','P312',(/0.D0,0.D0,0.D0/)) 
      SG(150)=SG_TYPE('GH','P321',(/0.D0,0.D0,0.D0/)) 
      SG(151)=SG_TYPE('GH','P3_112',(/0.D0,0.D0,(1.D0/6.D0)/))
      SG(152)=SG_TYPE('GH','P3_121',(/0.D0,0.D0,0.D0/)) 
      SG(153)=SG_TYPE('GH','P3_212',(/0.D0,0.D0,-(1.D0/6.D0)/))
      SG(154)=SG_TYPE('GH','P3_221',(/0.D0,0.D0,0.D0/))
      SG(155)=SG_TYPE('GRH','R32',(/0.D0,0.D0,0.D0/))     
      SG(156)=SG_TYPE('GH','P3M1',(/0.D0,0.D0,0.D0/))     
      SG(157)=SG_TYPE('GH','P31M',(/0.D0,0.D0,0.D0/))     
      SG(158)=SG_TYPE('GH','P3C1',(/0.D0,0.D0,0.D0/))     
      SG(159)=SG_TYPE('GH','P31C',(/0.D0,0.D0,0.D0/))  
      SG(160)=SG_TYPE('GRH','R3M',(/0.D0,0.D0,0.D0/))
      SG(161)=SG_TYPE('GRH','R3C',(/0.D0,0.D0,0.D0/))
      SG(162)=SG_TYPE('GH','P\BAR{3}1M',(/0.D0,0.D0,0.D0/)) 
      SG(163)=SG_TYPE('GH','P\BAR{3}1C',(/0.D0,0.D0,0.D0/)) 
      SG(164)=SG_TYPE('GH','P\BAR{3}M1',(/0.D0,0.D0,0.D0/)) 
      SG(165)=SG_TYPE('GH','P\BAR{3}C1',(/0.D0,0.D0,0.D0/))    
      SG(166)=SG_TYPE('GRH','R\BAR{3}M',(/0.D0,0.D0,0.D0/)) 
      SG(167)=SG_TYPE('GRH','R\BAR{3}C',(/0.D0,0.D0,0.D0/)) 
      SG(168)=SG_TYPE('GH','P6',(/0.D0,0.D0,0.D0/))
      SG(169)=SG_TYPE('GH','P6_1',(/0.D0,0.D0,0.D0/))
      SG(170)=SG_TYPE('GH','P6_5',(/0.D0,0.D0,0.D0/))
      SG(171)=SG_TYPE('GH','P6_2',(/0.D0,0.D0,0.D0/))
      SG(172)=SG_TYPE('GH','P6_4',(/0.D0,0.D0,0.D0/))
      SG(173)=SG_TYPE('GH','P6_3',(/0.D0,0.D0,0.D0/))
      SG(174)=SG_TYPE('GH','P\BAR{6}',(/0.D0,0.D0,0.D0/))
      SG(175)=SG_TYPE('GH','P6/M',(/0.D0,0.D0,0.D0/))
      SG(176)=SG_TYPE('GH','P6_3/M',(/0.D0,0.D0,0.25D0/))
      SG(177)=SG_TYPE('GH','P622',(/0.D0,0.D0,0.D0/))
      SG(178)=SG_TYPE('GH','P6_122',(/0.D0,0.D0,0.D0/))
      SG(179)=SG_TYPE('GH','P6_522',(/0.D0,0.D0,0.D0/))     
      SG(180)=SG_TYPE('GH','P6_222',(/0.D0,0.D0,0.D0/))   
      SG(181)=SG_TYPE('GH','P6_422',(/0.D0,0.D0,0.D0/))   
      SG(182)=SG_TYPE('GH','P6_322',(/0.D0,0.D0,0.D0/)) 
      SG(183)=SG_TYPE('GH','P6MM',(/0.D0,0.D0,0.D0/)) 
      SG(184)=SG_TYPE('GH','P6CC',(/0.D0,0.D0,0.D0/)) 
      SG(185)=SG_TYPE('GH','P6_3CM',(/0.D0,0.D0,0.D0/)) 
      SG(186)=SG_TYPE('GH','P6_3MC',(/0.D0,0.D0,0.D0/)) 
      SG(187)=SG_TYPE('GH','P\BAR{6}M2',(/0.D0,0.D0,0.D0/)) 
      SG(188)=SG_TYPE('GH','P\BAR{6}C2',(/0.D0,0.D0,0.25D0/)) 
      SG(189)=SG_TYPE('GH','P\BAR{6}2M',(/0.D0,0.D0,0.D0/)) 
      SG(190)=SG_TYPE('GH','P\BAR{6}2C',(/0.D0,0.D0,0.25D0/)) 
      SG(191)=SG_TYPE('GH','P6/MMM',(/0.D0,0.D0,0.D0/)) 
      SG(192)=SG_TYPE('GH','P6/MMC',(/0.D0,0.D0,0.D0/)) 
      SG(193)=SG_TYPE('GH','P6/MCM',(/0.D0,0.D0,0.D0/)) 
      SG(194)=SG_TYPE('GH','P6_3/MMC',(/0.D0,0.D0,0.D0/)) 
      SG(195)=SG_TYPE('GC','P23',(/0.D0,0.D0,0.D0/)) 
      SG(196)=SG_TYPE('GCF','F23',(/0.D0,0.D0,0.D0/))
      SG(197)=SG_TYPE('GCV','I23',(/0.D0,0.D0,0.D0/))
      SG(198)=SG_TYPE('GC','P2_13',(/0.D0,0.D0,0.D0/))
      SG(199)=SG_TYPE('GCV','I2_13',(/0.D0,0.D0,0.D0/))
      SG(200)=SG_TYPE('GC','PM3',(/0.D0,0.D0,0.D0/))
      SG(201)=SG_TYPE('GC','PN3',(/0.D0,0.D0,0.D0/))
      SG(202)=SG_TYPE('GCF','FM3',(/0.D0,0.D0,0.D0/))
      SG(203)=SG_TYPE('GCF','FD3',(/0.D0,0.D0,0.D0/))
      SG(204)=SG_TYPE('GCV','IM3',(/0.D0,0.D0,0.D0/))
      SG(205)=SG_TYPE('GC','PA3',(/0.D0,0.D0,0.D0/))
      SG(206)=SG_TYPE('GCV','IA3',(/0.D0,0.D0,0.D0/))
      SG(207)=SG_TYPE('GC','P432',(/0.D0,0.D0,0.D0/))
      SG(208)=SG_TYPE('GC','P4_232',(/0.D0,0.D0,0.D0/))
      SG(209)=SG_TYPE('GCF','F432',(/0.D0,0.D0,0.D0/))
      SG(210)=SG_TYPE('GCF','F4_132',(/0.D0,0.D0,0.D0/))
      SG(211)=SG_TYPE('GCV','I432',(/0.D0,0.D0,0.D0/))
      SG(212)=SG_TYPE('GC','P4_332',(/0.D0,0.D0,0.D0/))
      SG(213)=SG_TYPE('GC','P4_132',(/0.D0,0.D0,0.D0/))
      SG(214)=SG_TYPE('GCV','I4_132',(/0.D0,0.D0,0.D0/))
      SG(215)=SG_TYPE('GC','P\BAR{4}3M',(/0.D0,0.D0,0.D0/))
      SG(216)=SG_TYPE('GCF','F\BAR{4}3M',(/0.D0,0.D0,0.D0/))
      SG(217)=SG_TYPE('GCV','I\BAR{4}3M',(/0.D0,0.D0,0.D0/))
      SG(218)=SG_TYPE('GC','P\BAR{4}3N',(/0.D0,0.D0,0.D0/))
      SG(219)=SG_TYPE('GCF','F\BAR{4}3C',(/0.D0,0.D0,0.D0/))
      SG(220)=SG_TYPE('GCV','I\BAR{4}3D',(/0.D0,0.D0,0.D0/))
      SG(221)=SG_TYPE('GC','PM3M',(/0.D0,0.D0,0.D0/))
      SG(222)=SG_TYPE('GC','PN3N',(/0.D0,0.D0,0.D0/))
      SG(223)=SG_TYPE('GC','PM3N',(/0.D0,0.D0,0.D0/))
      SG(224)=SG_TYPE('GC','PN3M',(/0.D0,0.D0,0.D0/))
      SG(225)=SG_TYPE('GCF','FM3M',(/0.D0,0.D0,0.D0/))
      SG(226)=SG_TYPE('GCF','FM3C',(/0.25D0,0.25D0,0.25D0/))
      SG(227)=SG_TYPE('GCF','FD3M',(/0.D0,0.D0,0.D0/))
      SG(228)=SG_TYPE('GCF','FD3C',(/0.D0,0.D0,0.D0/))
      SG(229)=SG_TYPE('GCV','IM3M',(/0.D0,0.D0,0.D0/))
      SG(230)=SG_TYPE('GCV','IA3D',(/0.D0,0.D0,0.D0/))
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPACEGROUP$SETI4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE SPACEGROUP_MODULE, ONLY : ISPACEGROUP
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      CALL SPACEGROUP_INI()
      IF(ID.EQ.'SPACEGROUP') THEN
        ISPACEGROUP=VAL
        IF(ISPACEGROUP.LT.1.OR.ISPACEGROUP.GT.230) THEN
          CALL ERROR$MSG('SPACE GROUP NUMBER OUT OF RANGE')
          CALL ERROR$I4VAL('ISPACEGROUP',ISPACEGROUP)
          CALL ERROR$STOP('SPACEGROUP$SETI4')
        END IF
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SPACEGROUP$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPACEGROUP$GETCH(ID,VAL)
      USE SPACEGROUP_MODULE, ONLY : ISPACEGROUP,SG
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(OUT):: VAL
!     **************************************************************************
      CALL SPACEGROUP_INI()
      IF(ISPACEGROUP.EQ.0) THEN
        CALL ERROR$MSG('NO SPACE GROUP SELECTED')
        CALL ERROR$MSG('CALL SPACEGROUP$SETCH WITH ID="SPACEGROUP:"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SPACEGROUP$GETCH')
      END IF
!
      IF(ID.EQ.'BRAVAIS') THEN
        VAL=SG(ISPACEGROUP)%BRAVAIS
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SPACEGROUP$GETCH')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPACEGROUP$GENERATORS(ID,NOPX,NOP,OPERATION,C)
!     **************************************************************************
!     **   IMPLEMENTS TABLE 3.4 OF BRADLEY CRACKNELL                          **
!     ** SYMMETRY ELEMENTS EXPLAINED IN TABLE 1.2 OF BRADLEY CRACKNELL        **
!     **                                                                      **
!     ** THE TRANSOFRMATION OF RELATIVE COORDINATES IS                        **
!     **    RPRIME = OPERATION*R + C                                          **
!     **                                                                      **
!     **************************************************************************
      USE SPACEGROUP_MODULE, ONLY : ISPACEGROUP
      IMPLICIT NONE
      CHARACTER(4),INTENT(IN) :: ID ! CAN BE 'REAL' OR 'RECI'
      INTEGER(4)  ,INTENT(IN) :: NOPX
      INTEGER(4)  ,INTENT(OUT):: NOP
      INTEGER(4)  ,INTENT(OUT):: OPERATION(3,3,NOPX)
      REAL(8)     ,INTENT(OUT):: C(3,NOPX)  ! CENTER OF OPERATION
      REAL(8)     ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)              :: E(3,3)      ! IDENTITY
      INTEGER(4)              :: INV(3,3)    ! INVERSION
      INTEGER(4)              :: C4XP(3,3),  C4XM(3,3)
      INTEGER(4)              :: C4YP(3,3),  C4YM(3,3)
      INTEGER(4)              :: C4ZP(3,3),  C4ZM(3,3)
      INTEGER(4)              :: C2Z(3,3)    ! TWO-FOLD ROTATION ABOUT Z-AXIS
      INTEGER(4)              :: C2X(3,3)    ! TWO-FOLD ROTATION ABOUT X-AXIS
      INTEGER(4)              :: C2Y(3,3)    ! TWO-FOLD ROTATION ABOUT Y-AXIS
      INTEGER(4)              :: C2A(3,3)
      INTEGER(4)              :: C2B(3,3)
      INTEGER(4)              :: C2C(3,3)
      INTEGER(4)              :: C2D(3,3)
      INTEGER(4)              :: C2E(3,3)
      INTEGER(4)              :: C2F(3,3)
      INTEGER(4)              :: C6P(3,3), C6M(3,3)
      INTEGER(4)              :: C3P(3,3), C3M(3,3)
      INTEGER(4)              :: C2(3,3)
      INTEGER(4)              :: C21S(3,3), C22S(3,3), C23S(3,3)
      INTEGER(4)              :: C21SS(3,3), C22SS(3,3), C23SS(3,3)
      INTEGER(4)              :: C31P(3,3),C32P(3,3),C33P(3,3),C34P(3,3) 
      INTEGER(4)              :: C31M(3,3),C32M(3,3),C33M(3,3),C34M(3,3)
      REAL(8)                 :: MAT(3,3)
      REAL(8)                 :: PIH ! PIH=PI/2
      REAL(8)                 :: RBAS(3,3),TTAUT(3,3),TTAUTINV(3,3)
      INTEGER(4)              :: I
      CHARACTER(3)            :: BRAVAIS
!     **************************************************************************
      CALL SPACEGROUP_INI()
      IF(ISPACEGROUP.EQ.0) THEN 
        CALL ERROR$MSG('NO SPACE GROUP SELECTED')
        CALL ERROR$MSG('CALL SPACEGROUP$SETCH WITH ID="SPACEGROUP:"')
        CALL ERROR$STOP('SPACEGROUP$GENERATORS')
      END IF
!
!     ==========================================================================
!     == SELECT THE BRAVAIS LATTICE FOR THE SPECIFIED POINT GROUP
!     ==========================================================================
      CALL SPACEGROUP$GETCH('BRAVAIS',BRAVAIS)
!
!     ==========================================================================
!     == GENERATORS IN RECIPROCAL SPACE FOR THE BRAVAIS LATTICES              ==
!     ==========================================================================
      E=0
      INV=0
      C4XP=0;  C4XM=0; C4YP=0;  C4YM=0; C4ZP=0;  C4ZM=0
      C2Z=0; C2X=0; C2Y=0; C2A=0; C2B=0; C2C=0; C2D=0; C2E=0; C2F=0
      C6P=0; C6M=0
      C3P=0; C3M=0; C2=0; C21S=0; C22S=0; C23S=0; C21SS=0; C22SS=0; C23SS=0
      C31P=0; C32P=0;C33P=0;C34P=0; C31M=0;C32M=0;C33M=0;C34M=0

      E(:,1)=(/1,0,0/); E(:,3)=(/0,1,0/); E(:,3)=(/0,0,1/); 
      INV(:,1)=(/-1,0,0/); INV(:,2)=(/0,-1,0/); INV(:,3)=(/0,0,-1/); 
!
      IF(BRAVAIS.EQ.'G1') THEN ! MONOCLINIC
        CONTINUE

      ELSE IF(BRAVAIS.EQ.'GM ') THEN ! MONOCLINIC
        C2Z(:,1)=(/-1,0,0/); C2Z(:,2)=(/0,-1,0/); C2Z(:,3)=(/0,0,1/)
        !
      ELSE IF(BRAVAIS.EQ.'GMB') THEN
        C2Z(:,1)=(/-1,0,0/); C2Z(:,2)=(/0,0,-1/); C2Z(:,3)=(/0,-1,0/)
        !
      ELSE IF(BRAVAIS.EQ.'GO') THEN
        C2X(:,1)=(/-1,0,0/); C2X(:,2)=(/0,1,0/); C2X(:,3)=(/0,0,-1/)
        C2Y(:,1)=(/1,0,0/); C2Y(:,2)=(/0,-1,0/); C2Y(:,3)=(/0,0,-1/)
        C2Z(:,1)=(/-1,0,0/); C2Z(:,2)=(/0,-1,0/); C2Z(:,3)=(/0,0,1/)
        !
      ELSE IF(BRAVAIS.EQ.'GOB') THEN
        C2X(:,1)=(/0,1,0/); C2X(:,2)=(/1,0,0/); C2X(:,3)=(/0,0,-1/)
        C2Y(:,1)=(/0,-1,0/); C2Y(:,2)=(/-1,0,0/); C2Y(:,3)=(/0,0,-1/)
        C2Z(:,1)=(/-1,0,0/); C2Z(:,2)=(/0,-1,0/); C2Z(:,3)=(/0,0,1/) 
        !
      ELSE IF(BRAVAIS.EQ.'GOV') THEN
        C2X(:,1)=(/0,-1,1/); C2X(:,2)=(/0,-1,0/); C2X(:,3)=(/1,-1,0/)
        C2Y(:,1)=(/-1,0,0/); C2Y(:,2)=(/-1,0,1/); C2Y(:,3)=(/-1,1,0/)
        C2Z(:,1)=(/0,1,-1/); C2Z(:,2)=(/1,0,-1/); C2Z(:,3)=(/0,0,-1/) 
        !
      ELSE IF(BRAVAIS.EQ.'GOF') THEN
        C2X(:,1)=(/0,0,1/); C2X(:,2)=(/-1,-1,-1/); C2X(:,3)=(/1,0,0/)
        C2Y(:,1)=(/-1,-1,-1/); C2Y(:,2)=(/0,0,1/); C2Y(:,3)=(/0,1,0/)
        C2Z(:,1)=(/0,1,0/); C2Z(:,2)=(/1,0,0/); C2Z(:,3)=(/-1,-1,-1/) 
        !
      ELSE IF(BRAVAIS.EQ.'GQ') THEN
        C4ZP(:,1)=(/0,1,0/);  C4ZP(:,2)=(/-1,0,0/);  C4ZP(:,3)=(/0,0,1/)
        C2Z=MATMUL(C4ZP,C4ZP)
        C4ZM=MATMUL(C4ZP,C2Z)
        C2X(:,1)=(/1,0,0/); C2X(:,2)=(/0,-1,0/); C2X(:,3)=(/0,0,-1/)
        C2Y(:,1)=(/-1,0,0/); C2Y(:,2)=(/0,1,0/); C2Y(:,3)=(/0,0,-1/)
        C2A(:,1)=(/0,1,0/); C2A(:,2)=(/1,0,0/); C2A(:,3)=(/0,0,-1/) 
        C2B(:,1)=(/0,-1,0/); C2B(:,2)=(/-1,0,0/); C2B(:,3)=(/0,0,-1/) 
        !
      ELSE IF(BRAVAIS.EQ.'GQV') THEN
        C4ZP(:,1)=(/1,0,-1/);  C4ZP(:,2)=(/1,0,0/);  C4ZP(:,3)=(/1,-1,0/)
        C2Z=MATMUL(C4ZP,C4ZP)
        C4ZM=MATMUL(C4ZP,C2Z)
        C2X(:,1)=(/-1,0,0/); C2X(:,2)=(/-1,0,1/); C2X(:,3)=(/-1,1,0/)
        C2Y(:,1)=(/0,-1,1/); C2Y(:,2)=(/0,-1,0/); C2Y(:,3)=(/1,-1,0/)
        C2A(:,1)=(/-1,0,1/); C2A(:,2)=(/0,-1,1/); C2A(:,3)=(/0,0,1/) 
        C2B(:,1)=(/0,-1,0/); C2B(:,2)=(/-1,0,0/); C2B(:,3)=(/0,0,-1/) 
        !
      ELSE IF(BRAVAIS.EQ.'GRH') THEN
        C3P(:,1)=(/0,1,0/); C3P(:,2)=(/0,0,1/); C3P(:,3)=(/1,0,0/)
        C3M=MATMUL(C3P,C3P) 
        C21S(:,1)=(/-1,0,0/); C21S(:,2)=(/0,0,-1/); C21S(:,3)=(/0,-1,0/) 
        C22S(:,1)=(/0,0,-1/); C22S(:,2)=(/0,-1,0/); C22S(:,3)=(/-1,0,0/)
        C23S(:,1)=(/0,-1,0/); C23S(:,2)=(/-1,0,0/); C23S(:,3)=(/0,0,-1/)
        !
      ELSE IF(BRAVAIS.EQ.'GRH') THEN
        C6P(:,1)=(/0,1,0/); C6P(:,2)=(/-1,1,0/); C6P(:,3)=(/0,0,1/)
        C3P=MATMUL(C6P,C6P)
        C2(:,1)=(/-1,0,0/); C2(:,2)=(/0,-1,0/); C2(:,3)=(/0,0,1/) 
        C3M(:,1)=(/0,-1,0/); C3M(:,2)=(/1,-1,0/); C3M(:,3)=(/0,0,1/)
        C6M=MATMUL(C3M,C6P) 
        C21S(:,1)=(/-1,1,0/); C21S(:,2)=(/0,1,0/); C21S(:,3)=(/0,0,-1/) 
        C22S(:,1)=(/1,0,0/); C22S(:,2)=(/1,-1,0/); C22S(:,3)=(/0,0,-1/)
        C23S(:,1)=(/0,-1,0/); C23S(:,2)=(/-1,0,0/); C23S(:,3)=(/0,0,-1/)
        C21SS(:,1)=(/1,-1,0/); C21SS(:,2)=(/0,-1,0/); C21SS(:,3)=(/0,0,-1/) 
        C22SS(:,1)=(/-1,0,0/); C22SS(:,2)=(/-1,1,0/); C22SS(:,3)=(/0,0,-1/)
        C23SS(:,1)=(/0,1,0/); C23SS(:,2)=(/1,0,0/); C23SS(:,3)=(/0,0,-1/)  
        !
      ELSE IF(BRAVAIS.EQ.'GC') THEN
        C2X(:,1)=(/1,0,0/); C2X(:,2)=(/0,-1,0/); C2X(:,3)=(/0,0,-1/)
        C2Y(:,1)=(/-1,0,0/); C2Y(:,2)=(/0,1,0/); C2Y(:,3)=(/0,0,-1/)
        C2Z(:,1)=(/-1,0,0/); C2Z(:,2)=(/0,-1,0/); C2Z(:,3)=(/0,0,1/)         
        C31P(:,1)=(/0,1,0/); C31P(:,2)=(/0,0,1/); C31P(:,3)=(/1,0,0/)
        C32P(:,1)=(/0,1,0/); C32P(:,2)=(/0,0,-1/); C32P(:,3)=(/-1,0,0/)
        C33P(:,1)=(/0,-1,0/); C33P(:,2)=(/0,0,1/); C33P(:,3)=(/-1,0,0/)
        C34P(:,1)=(/0,-1,0/); C34P(:,2)=(/0,0,-1/); C34P(:,3)=(/1,0,0/)
        C31M(:,1)=(/0,0,1/); C31M(:,2)=(/1,0,0/); C31M(:,3)=(/0,1,0/)
        C32M(:,1)=(/0,0,-1/); C32M(:,2)=(/1,0,0/); C32M(:,3)=(/0,-1,0/)
        C33M(:,1)=(/0,0,-1/); C33M(:,2)=(/-1,0,0/); C33M(:,3)=(/0,1,0/)
        C34M(:,1)=(/0,0,1/); C34M(:,2)=(/-1,0,0/); C34M(:,3)=(/0,-1,0/)
        C4XP(:,1)=(/1,0,0/); C4XP(:,2)=(/0,0,1/); C4XP(:,3)=(/0,-1,0/)
        C4YP(:,1)=(/0,0,-1/); C4YP(:,2)=(/0,1,0/); C4YP(:,3)=(/1,0,0/)
        C4ZP(:,1)=(/0,1,0/); C4ZP(:,2)=(/-1,0,0/); C4ZP(:,3)=(/0,0,1/)
        C4XM(:,1)=(/1,0,0/); C4XM(:,2)=(/0,0,-1/); C4XM(:,3)=(/0,1,0/)
        C4YM(:,1)=(/0,0,1/); C4YM(:,2)=(/0,1,0/); C4YM(:,3)=(/-1,0,0/)
        C4ZM(:,1)=(/0,-1,0/); C4ZM(:,2)=(/1,0,0/); C4ZM(:,3)=(/0,0,1/) 
        C2A(:,1)=(/0,1,0/); C2A(:,2)=(/1,0,0/); C2A(:,3)=(/0,0,-1/) 
        C2B(:,1)=(/0,-1,0/); C2B(:,2)=(/-1,0,0/); C2B(:,3)=(/0,0,-1/)
        C2C(:,1)=(/0,0,1/); C2C(:,2)=(/0,-1,0/); C2C(:,3)=(/1,0,0/) 
        C2D(:,1)=(/-1,0,0/); C2D(:,2)=(/0,0,1/); C2D(:,3)=(/0,1,0/)
        C2E(:,1)=(/0,0,-1/); C2E(:,2)=(/0,-1,0/); C2E(:,3)=(/-1,0,0/) 
        C2F(:,1)=(/-1,0,0/); C2F(:,2)=(/0,0,-1/); C2F(:,3)=(/0,-1,0/)

      ELSE IF(BRAVAIS.EQ.'GCF') THEN
        C2X(:,1)=(/-1,-1,-1/); C2X(:,2)=(/0,0,1/); C2X(:,3)=(/0,1,0/)
        C2Y(:,1)=(/0,0,1/); C2Y(:,2)=(/-1,-1,-1/); C2Y(:,3)=(/1,0,0/)
        C2Z(:,1)=(/0,1,0/); C2Z(:,2)=(/1,0,0/); C2Z(:,3)=(/-1,-1,-1/)         
        C31P(:,1)=(/0,1,0/); C31P(:,2)=(/0,0,1/); C31P(:,3)=(/1,0,0/)
        C32P(:,1)=(/-1,-1,-1/); C32P(:,2)=(/1,0,0/); C32P(:,3)=(/0,0,1/)
        C33P(:,1)=(/1,0,0/); C33P(:,2)=(/-1,-1,-1/); C33P(:,3)=(/0,1,0/)
        C34P(:,1)=(/0,0,1/); C34P(:,2)=(/0,1,0/); C34P(:,3)=(/-1,-1,-1/)
        C31M(:,1)=(/0,0,1/); C31M(:,2)=(/1,0,0/); C31M(:,3)=(/0,1,0/)
        C32M(:,1)=(/0,1,0/); C32M(:,2)=(/-1,-1,-1/); C32M(:,3)=(/0,0,1/)
        C33M(:,1)=(/1,0,0/); C33M(:,2)=(/0,0,1/); C33M(:,3)=(/-1,-1,-1/)
        C34M(:,1)=(/-1,-1,-1/); C34M(:,2)=(/0,1,0/); C34M(:,3)=(/1,0,0/)
        C4XP(:,1)=(/0,0,-1/); C4XP(:,2)=(/-1,0,0/); C4XP(:,3)=(/1,1,1/)
        C4YP(:,1)=(/1,1,1/); C4YP(:,2)=(/-1,0,0/); C4YP(:,3)=(/0,-1,0/)
        C4ZP(:,1)=(/0,0,-1/); C4ZP(:,2)=(/1,1,1/); C4ZP(:,3)=(/0,-1,0/)
        C4XM(:,1)=(/0,-1,0/); C4XM(:,2)=(/1,1,1/); C4XM(:,3)=(/-1,0,0/)
        C4YM(:,1)=(/0,-1,0/); C4YM(:,2)=(/0,0,-1/); C4YM(:,3)=(/1,1,1/)
        C4ZM(:,1)=(/1,1,1/); C4ZM(:,2)=(/0,0,-1/); C4ZM(:,3)=(/-1,0,0/) 
        C2A(:,1)=(/-1,0,0/); C2A(:,2)=(/0,-1,0/); C2A(:,3)=(/1,1,1/) 
        C2B(:,1)=(/0,-1,0/); C2B(:,2)=(/-1,0,0/); C2B(:,3)=(/0,0,-1/)
        C2C(:,1)=(/-1,0,0/); C2C(:,2)=(/1,1,1/); C2C(:,3)=(/0,0,-1/) 
        C2D(:,1)=(/1,1,1/); C2D(:,2)=(/0,-1,0/); C2D(:,3)=(/0,0,-1/)
        C2E(:,1)=(/0,0,-1/); C2E(:,2)=(/0,-1,0/); C2E(:,3)=(/-1,0,0/) 
        C2F(:,1)=(/-1,0,0/); C2F(:,2)=(/0,0,-1/); C2F(:,3)=(/0,-1,0/)      
        !
      ELSE IF(BRAVAIS.EQ.'GCV') THEN
        C2X(:,1)=(/-1,0,0/); C2X(:,2)=(/-1,0,1/); C2X(:,3)=(/-1,1,0/)
        C2Y(:,1)=(/0,-1,1/); C2Y(:,2)=(/0,-1,0/); C2Y(:,3)=(/1,-1,0/)
        C2Z(:,1)=(/0,1,-1/); C2Z(:,2)=(/1,0,-1/); C2Z(:,3)=(/0,0,-1/)         
        C31P(:,1)=(/0,1,0/); C31P(:,2)=(/0,0,1/); C31P(:,3)=(/1,0,0/)
        C32P(:,1)=(/0,-1,0/); C32P(:,2)=(/1,-1,0/); C32P(:,3)=(/0,-1,1/)
        C33P(:,1)=(/1,0,-1/); C33P(:,2)=(/0,0,-1/); C33P(:,3)=(/0,1,-1/)
        C34P(:,1)=(/-1,0,1/); C34P(:,2)=(/-1,1,0/); C34P(:,3)=(/-1,0,0/)
        C31M(:,1)=(/0,0,1/); C31M(:,2)=(/1,0,0/); C31M(:,3)=(/0,1,0/)
        C32M(:,1)=(/-1,1,0/); C32M(:,2)=(/-1,0,0/); C32M(:,3)=(/-1,0,1/)
        C33M(:,1)=(/1,-1,0/); C33M(:,2)=(/0,-1,1/); C33M(:,3)=(/0,-1,0/)
        C34M(:,1)=(/0,0,-1/); C34M(:,2)=(/0,1,-1/); C34M(:,3)=(/1,0,-1/)
        C4XP(:,1)=(/0,1,-1/); C4XP(:,2)=(/-1,1,0/); C4XP(:,3)=(/0,1,0/)
        C4YP(:,1)=(/0,0,1/); C4YP(:,2)=(/-1,0,1/); C4YP(:,3)=(/0,-1,1/)
        C4ZP(:,1)=(/1,0,-1/); C4ZP(:,2)=(/1,0,0/); C4ZP(:,3)=(/1,-1,0/)
        C4XM(:,1)=(/0,-1,1/); C4YM(:,2)=(/0,0,1/); C4YM(:,3)=(/-1,0,1/)
        C4YM(:,1)=(/1,-1,0/); C4YM(:,2)=(/1,0,-1/); C4YM(:,3)=(/1,0,0/)
        C4ZM(:,1)=(/0,1,0/); C4ZM(:,2)=(/0,1,-1/); C4ZM(:,3)=(/-1,1,0/) 
        C2A(:,1)=(/-1,0,1/); C2A(:,2)=(/0,-1,1/); C2A(:,3)=(/0,0,1/) 
        C2B(:,1)=(/0,-1,0/); C2B(:,2)=(/-1,0,0/); C2B(:,3)=(/0,0,-1/)
        C2C(:,1)=(/-1,1,0/); C2C(:,2)=(/0,1,0/); C2C(:,3)=(/0,1,-1/) 
        C2D(:,1)=(/1,0,0/); C2D(:,2)=(/1,-1,0/); C2D(:,3)=(/1,0,-1/)
        C2E(:,1)=(/0,0,-1/); C2E(:,2)=(/0,-1,0/); C2E(:,3)=(/-1,0,0/) 
        C2F(:,1)=(/-1,0,0/); C2F(:,2)=(/0,0,-1/); C2F(:,3)=(/0,-1,0/)
      ELSE
        CALL ERROR$MSG('LATTICE SYMBOL NOT RECOGNIZED')
        CALL ERROR$STOP('SPACEGROUP$GENERATORS')
      END IF
!
!     ==========================================================================
!     == BUILD UP ALL POINT GROUPS                                            ==
!     ==========================================================================
      OPERATION=0
      C=0.D0
      NOP=0
      IF(ISPACEGROUP.EQ.1) THEN
      ELSE IF(ISPACEGROUP.EQ.2) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.3) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
      ELSE IF(ISPACEGROUP.EQ.4) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.5) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
      ELSE IF(ISPACEGROUP.EQ.6) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Z
      ELSE IF(ISPACEGROUP.EQ.7) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Z ; C(:,NOP)=(/0.5D0,0.D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.8) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Z 
      ELSE IF(ISPACEGROUP.EQ.9) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Z ; C(:,NOP)=(/0.5D0,0.D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.10) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.11) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.12) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.13) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV ; C(:,NOP)=(/0.5D0,0.D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.5D0,0.D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.14) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV ; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.5D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.15) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=INV ; C(:,NOP)=(/0.5D0,0.D0,0.0D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.5D0,0.D0,0.0D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.16) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
      ELSE IF(ISPACEGROUP.EQ.17) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.18) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.19) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.20) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.21) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X         
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
      ELSE IF(ISPACEGROUP.EQ.22) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X         
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
      ELSE IF(ISPACEGROUP.EQ.23) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X         
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
      ELSE IF(ISPACEGROUP.EQ.24) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.25) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y 
      ELSE IF(ISPACEGROUP.EQ.26) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
      ELSE IF(ISPACEGROUP.EQ.27) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
      ELSE IF(ISPACEGROUP.EQ.28) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.29) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.D0,0.5D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.30) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.31) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.32) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.33) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.34) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.35) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y
      ELSE IF(ISPACEGROUP.EQ.36) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.37) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.38) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y
      ELSE IF(ISPACEGROUP.EQ.39) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Z; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.5D0,0.D0/) 
      ELSE IF(ISPACEGROUP.EQ.40) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Z; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.D0,05.D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.41) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Z; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/) 
      ELSE IF(ISPACEGROUP.EQ.42) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y
      ELSE IF(ISPACEGROUP.EQ.43) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X;C(:,NOP)=(/0.D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y;C(:,NOP)=(/0.5D0,0.D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X;C(:,NOP)=(/0.75D0,0.75D0,0.75D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y;C(:,NOP)=(/0.75D0,0.75D0,0.75D0/)
      ELSE IF(ISPACEGROUP.EQ.44) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y
      ELSE IF(ISPACEGROUP.EQ.45) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.46) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.D0,0.5D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2Y; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.47) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.48) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/) 
      ELSE IF(ISPACEGROUP.EQ.49) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.50) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/) 
      ELSE IF(ISPACEGROUP.EQ.51) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.52) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.53) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.D0,0.5D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.54) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.55) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.56) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.57) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.58) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.59) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.60) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.61) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.5D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.62) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.5D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.D0,0.D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.63) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.64) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.65) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.66) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.67) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.68) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.69) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.70) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.25D0,0.25D0,0.25D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.75D0,0.75D0,0.75D0/)
      ELSE IF(ISPACEGROUP.EQ.71) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.72) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.73) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.74) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Y; C(:,NOP)=(/0.D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.75) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
      ELSE IF(ISPACEGROUP.EQ.76) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.25D0/)
      ELSE IF(ISPACEGROUP.EQ.77) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.78) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.75D0/)
      ELSE IF(ISPACEGROUP.EQ.79) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
      ELSE IF(ISPACEGROUP.EQ.80) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.75D0,0.25D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.81) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
      ELSE IF(ISPACEGROUP.EQ.82) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
      ELSE IF(ISPACEGROUP.EQ.83) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.84) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.85) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.86) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.87) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.88) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.75D0,0.25D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.75D0,0.25D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.89) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
      ELSE IF(ISPACEGROUP.EQ.90) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.91) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.25D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.25D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.92) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.25D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.25D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.75D0/)
      ELSE IF(ISPACEGROUP.EQ.93) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
      ELSE IF(ISPACEGROUP.EQ.94) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.95) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.75D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.75D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.96) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.75D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.75D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.25D0/)
      ELSE IF(ISPACEGROUP.EQ.97) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
      ELSE IF(ISPACEGROUP.EQ.98) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.75D0,0.25D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.75D0,0.25D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.75D0,0.25D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.99) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
      ELSE IF(ISPACEGROUP.EQ.100) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.101) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.102) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.103) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.104) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.105) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
      ELSE IF(ISPACEGROUP.EQ.106) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.107) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
      ELSE IF(ISPACEGROUP.EQ.108) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.109) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.75D0,0.25D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X
      ELSE IF(ISPACEGROUP.EQ.110) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.75D0,0.25D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.111) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
      ELSE IF(ISPACEGROUP.EQ.112) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.113) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.114) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.115) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A
      ELSE IF(ISPACEGROUP.EQ.116) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.117) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.118) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.119) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A
      ELSE IF(ISPACEGROUP.EQ.120) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.121) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
      ELSE IF(ISPACEGROUP.EQ.122) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.25D0,0.75D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C4ZM
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.75D0,0.25D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.123) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.124) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.125) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.126) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.127) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.128) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.129) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.130) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.131) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.132) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.133) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.134) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP;  C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.135) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.136) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.137) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.138) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.139) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.140) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.141) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.142) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.5D0,0.D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
        NOP=NOP+1; OPERATION(:,:,NOP)=C4ZP; C(:,NOP)=(/0.D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.143) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
      ELSE IF(ISPACEGROUP.EQ.144) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P; C(:,NOP)=(/0.D0,0.D0,(1.D0/3.D0)/)
      ELSE IF(ISPACEGROUP.EQ.145) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P; C(:,NOP)=(/0.D0,0.D0,(2.D0/3.D0)/)
      ELSE IF(ISPACEGROUP.EQ.146) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
      ELSE IF(ISPACEGROUP.EQ.147) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C3M
      ELSE IF(ISPACEGROUP.EQ.148) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C3M
      ELSE IF(ISPACEGROUP.EQ.149) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S
      ELSE IF(ISPACEGROUP.EQ.150) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
        NOP=NOP+1; OPERATION(:,:,NOP)=C21SS
      ELSE IF(ISPACEGROUP.EQ.151) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P; C(:,NOP)=(/0.D0,0.D0,(1.D0/3.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S; C(:,NOP)=(/0.D0,0.D0,(2.D0/3.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P; C(:,NOP)=(/0.D0,0.D0,(1.D0/3.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S; C(:,NOP)=(/0.D0,0.D0,(1.D0/3.D0)/)
      ELSE IF(ISPACEGROUP.EQ.152) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P; C(:,NOP)=(/0.D0,0.D0,(1.D0/3.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21SS; C(:,NOP)=(/0.D0,0.D0,(2.D0/3.D0)/)
      ELSE IF(ISPACEGROUP.EQ.153) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P; C(:,NOP)=(/0.D0,0.D0,(2.D0/3.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S; C(:,NOP)=(/0.D0,0.D0,(1.D0/3.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P; C(:,NOP)=(/0.D0,0.D0,(2.D0/3.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S; C(:,NOP)=(/0.D0,0.D0,(2.D0/3.D0)/)
      ELSE IF(ISPACEGROUP.EQ.154) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P; C(:,NOP)=(/0.D0,0.D0,(2.D0/3.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21SS; C(:,NOP)=(/0.D0,0.D0,(1.D0/3.D0)/)
      ELSE IF(ISPACEGROUP.EQ.155) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
        NOP=NOP+1; OPERATION(:,:,NOP)=C21SS
      ELSE IF(ISPACEGROUP.EQ.156) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS
      ELSE IF(ISPACEGROUP.EQ.157) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21S
      ELSE IF(ISPACEGROUP.EQ.158) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.159) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21S; C(:,NOP)=(/0.D0,0.5D0,0.D0/)
      ELSE IF(ISPACEGROUP.EQ.160) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21S
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS
      ELSE IF(ISPACEGROUP.EQ.161) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21S; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/) 
        NOP=NOP+1; OPERATION(:,:,NOP)=C3P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.162) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21S
      ELSE IF(ISPACEGROUP.EQ.163) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21S; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.164) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS
      ELSE IF(ISPACEGROUP.EQ.165) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.166) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21S
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS
      ELSE IF(ISPACEGROUP.EQ.167) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21S; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.168) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P
      ELSE IF(ISPACEGROUP.EQ.169) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,(1.D0/6.D0)/)
      ELSE IF(ISPACEGROUP.EQ.170) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,(5.D0/6.D0)/)
      ELSE IF(ISPACEGROUP.EQ.171) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,(1.D0/3.D0)/)
      ELSE IF(ISPACEGROUP.EQ.172) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,(2.D0/3.D0)/)
      ELSE IF(ISPACEGROUP.EQ.173) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.174) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6M
      ELSE IF(ISPACEGROUP.EQ.175) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2
      ELSE IF(ISPACEGROUP.EQ.176) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.177) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S
      ELSE IF(ISPACEGROUP.EQ.178) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,(1.D0/6.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,(1.D0/6.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21SS
      ELSE IF(ISPACEGROUP.EQ.179) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,(5.D0/6.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,(5.D0/6.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21SS
      ELSE IF(ISPACEGROUP.EQ.180) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,(1.D0/3.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S
      ELSE IF(ISPACEGROUP.EQ.181) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,(2.D0/3.D0)/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S
      ELSE IF(ISPACEGROUP.EQ.182) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21SS
      ELSE IF(ISPACEGROUP.EQ.183) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS
      ELSE IF(ISPACEGROUP.EQ.184) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.185) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.186) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS
      ELSE IF(ISPACEGROUP.EQ.187) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6M
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS
      ELSE IF(ISPACEGROUP.EQ.188) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6M
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C3P; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21SS; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.189) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6M
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21S
      ELSE IF(ISPACEGROUP.EQ.190) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6M
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21S; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C6M; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C21S; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.191) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.192) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.193) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.194) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C6P; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C21S; C(:,NOP)=(/0.D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.195) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.196) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.197) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.198) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.199) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.200) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.201) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.202) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.203) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.25D0,0.25D0,0.25D0/)
      ELSE IF(ISPACEGROUP.EQ.204) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.205) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.206) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.207) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.208) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.209) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.210) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.25D0,0.25D0,0.25D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.211) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.212) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.25D0,0.75D0,0.75D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.213) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.75D0,0.25D0,0.25D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.214) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.5D0,0.D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.215) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2A
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.216) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2A
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.217) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2A
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.218) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2A; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.219) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2A; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.220) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=-C2A; C(:,NOP)=(/0.5D0,0.D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
      ELSE IF(ISPACEGROUP.EQ.221) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.222) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.223) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.224) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
      ELSE IF(ISPACEGROUP.EQ.225) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.226) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.5D0,0.5D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.227) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.25D0,0.25D0,0.25D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.25D0,0.25D0,0.25D0/)
      ELSE IF(ISPACEGROUP.EQ.228) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.25D0,0.25D0,0.25D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV; C(:,NOP)=(/0.75D0,0.75D0,0.75D0/)
      ELSE IF(ISPACEGROUP.EQ.229) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      ELSE IF(ISPACEGROUP.EQ.221) THEN
        NOP=NOP+1; OPERATION(:,:,NOP)=C2Z; C(:,NOP)=(/0.5D0,0.D0,0.5D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2X; C(:,NOP)=(/0.5D0,0.5D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C2A; C(:,NOP)=(/0.5D0,0.D0,0.D0/)
        NOP=NOP+1; OPERATION(:,:,NOP)=C31P
        NOP=NOP+1; OPERATION(:,:,NOP)=INV
      END IF
!
!     ==========================================================================
!     ==  TRANSFORM RECIPROCAL SPACE OPERATIONS INTO REAL SPACE OPERATIONS
!     ==========================================================================
      IF(ID.EQ.'REAL') THEN
        CALL ERROR$MSG('REAL SPACE GENERATORES NOT YET IMPLEMENTED')
        CALL ERROR$STOP('BRILLOUIN$GENERATORS')
        PIH=PI/2.D0
        CALL SPACEGROUP$RBAS(BRAVAIS,1.D0,1.D0,1.D0,PIH,PIH,PIH,RBAS)
        TTAUT=MATMUL(TRANSPOSE(RBAS),RBAS)
        CALL LIB$INVERTR8(3,TTAUT,TTAUTINV)
        DO I=1,NOP
          MAT=REAL(OPERATION(:,:,I))
          MAT=MATMUL(TTAUTINV,MATMUL(MAT,TTAUT))
          OPERATION(:,:,I)=NINT(MAT)
          MAT=MAT-NINT(MAT)
          IF(MAXVAL(ABS(MAT)).GT.1.D-6) THEN
            CALL ERROR$MSG('SYMMETRY OP. INCONSISTENT WITH BRAVAIS LATTICE')
            CALL ERROR$STOP('BRILLOUIN$GENERATORS')
          END IF
        ENDDO
      ELSE IF(ID.EQ.'RECI') THEN
        C(:,:)=0.D0
      ELSE
        CALL ERROR$MSG('INVALID VALUE OF ID (MUST BE "REAL" OR "RECI")')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BRILLOUIN$GENERATORS')
      END IF
!
!     ==========================================================================
!     ==  CHECK CONSISTENCY OF POINT GROUP OPERATIONS                         ==
!     ==========================================================================
      DO I=1,NOP
        IF(SUM(ABS(OPERATION(:,:,I))).EQ.0) THEN
          CALL ERROR$MSG('INTERNAL ERROR: INVALID SYMMETRY OPERATION')
          CALL ERROR$MSG('SPACE GROUP NUMBER INCONSISTENT WTH BRAVAIS LATT.')
          CALL ERROR$MSG('SPACE GROUP NUMBER POINT GROUP OPERATIONS')
          CALL ERROR$I4VAL('SPACE GROUP NUMBER',ISPACEGROUP)
          CALL ERROR$CHVAL('BRAVAIS LATTICE',BRAVAIS)
          CALL ERROR$I4VAL('NUMBER OF GENERATORS',NOP)
          CALL ERROR$I4VAL('PROBLEM WITH GENERATOR',I)
          CALL ERROR$STOP('BRILLOUIN$GENERATORS')
        END IF
      ENDDO
      RETURN
      END SUBROUTINE SPACEGROUP$GENERATORS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPACEGROUP$RBAS(BRAVAIS,A0,B0,C0,ALPHA,BETA,GAMMA,RBAS)
!     **************************************************************************
!     ** LATTICE VECTORS RBAS FOR THE NAMED BRAVAIS LATTICE                   **
!     ** FOLLOWING TABLE 3.1 OF BRADLEY CRACKNELL                             **
!     **                                                                      **
!     ** LATTICE CONSTANTS AND ANGLES THAT ARE NOT NEEDED ARE NOT USED        **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(3),INTENT(IN)   :: BRAVAIS
      REAL(8)     ,INTENT(IN)   :: A0
      REAL(8)     ,INTENT(IN)   :: B0
      REAL(8)     ,INTENT(IN)   :: C0
      REAL(8)     ,INTENT(IN)   :: ALPHA
      REAL(8)     ,INTENT(IN)   :: BETA
      REAL(8)     ,INTENT(IN)   :: GAMMA
      REAL(8)     ,INTENT(OUT)  :: RBAS(3,3)
      REAL(8)                   :: A,B,C,ALPHA1,BETA1,GAMMA1
!     **************************************************************************
      IF(BRAVAIS.EQ.'G1') THEN ! TRICLINIC PRIMITIVE
!      __REMAP VARIABLES TO AVOID COMPLAINT REGARDING CONFLICTING INTENT(IN/OUT)
       A=A0
       B=B0
       C=C0
       ALPHA1=ALPHA
       BETA1=BETA
       GAMMA1=GAMMA
       CALL SPACEGROUP$ABCALPHABETAGAMMA(.TRUE.,A,B,C,ALPHA1,BETA1,GAMMA1,RBAS)

!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'G1') THEN ! MONOCLIN PRIMITIVE 
        RBAS(1,1)=0.0D0
        RBAS(2,1)=-1.0*B0
        RBAS(3,1)=0.0D0
        RBAS(1,2)=A0*DSIN(GAMMA)
        RBAS(2,2)=-A0*DCOS(GAMMA)  
        RBAS(3,2)=0.D0
        RBAS(1,3)=0.0D0
        RBAS(2,3)=0.0D0
        RBAS(3,3)=C0 
!
!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'GMB') THEN ! MONOCLIN BASE-CENTRED
        RBAS(1,1)=0.0D0
        RBAS(2,1)=-1.0D0*B0
        RBAS(3,1)=0.0D0
        RBAS(1,2)=0.5D0*A0*DSIN(GAMMA)
        RBAS(2,2)=-0.5D0*A0*DCOS(GAMMA)  
        RBAS(3,2)=-0.5D0*C0
        RBAS(1,3)=0.5D0*A0*DSIN(GAMMA)
        RBAS(2,3)=-0.5D0*A0*DCOS(GAMMA) 
        RBAS(3,3)=0.5D0*C0  
!
!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'GO') THEN ! ORTHORHOMBIC PRIMITIVE
        RBAS(1,1)=0.0D0
        RBAS(2,1)=-1.0D0*B0
        RBAS(3,1)=0.0D0
        RBAS(1,2)=A0
        RBAS(2,2)=0.0D0
        RBAS(3,2)=0.0D0
        RBAS(1,3)=0.0D0
        RBAS(2,3)=0.0D0
        RBAS(3,3)=C0
!
!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'GOB') THEN ! ORTHORHOMBIC BASE-CENTRED
        RBAS(1,1)=0.5D0*A0
        RBAS(2,1)=-0.5D0*B0
        RBAS(3,1)=0.0D0
        RBAS(1,2)=0.5D0*A0
        RBAS(2,2)=0.5D0*B0
        RBAS(3,2)=0.0D0
        RBAS(1,3)=0.0D0
        RBAS(2,3)=0.0D0
        RBAS(3,3)=C0
!
!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'GOV') THEN ! ORTHORHOMBIC BODY-CENTRED
        RBAS(1,1)=0.5D0*A0
        RBAS(2,1)=0.5D0*B0
        RBAS(3,1)=0.5D0*C0
        RBAS(1,2)=-0.5D0*A0
        RBAS(2,2)=-0.5D0*B0
        RBAS(3,2)=0.5D0*C0
        RBAS(1,3)=0.5D0*A0
        RBAS(2,3)=-0.5D0*B0
        RBAS(3,3)=-0.5D0*C0
!
!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'GOF') THEN ! ORTHORHOMBIC FACE-CENTRED
        RBAS(1,1)=0.5D0*A0
        RBAS(2,1)=0.0D0
        RBAS(3,1)=0.5D0*C0
        RBAS(1,2)=0.0D0
        RBAS(2,2)=-0.5D0*B0
        RBAS(3,2)=0.5D0*C0
        RBAS(1,3)=0.5D0*A0
        RBAS(2,3)=-0.5D0*B0
        RBAS(3,3)=0.0D0
!
!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'GQ') THEN ! TETRAGONAL PRIMITIVE
        RBAS(1,1)=A0
        RBAS(2,1)=0.0D0
        RBAS(3,1)=0.0D0
        RBAS(1,2)=0.0D0
        RBAS(2,2)=A0
        RBAS(3,2)=0.0D0
        RBAS(1,3)=0.0D0
        RBAS(2,3)=0.0D0
        RBAS(3,3)=C0
!
!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'GQV') THEN ! TETRAGONAL BODY-CENTRED
        RBAS(1,1)=-0.5D0*A0
        RBAS(2,1)=0.5D0*A0
        RBAS(3,1)=0.5D0*C0
        RBAS(1,2)=0.5D0*A0
        RBAS(2,2)=-0.5D0*A0
        RBAS(3,2)=0.5D0*C0
        RBAS(1,3)=0.5D0*A0
        RBAS(2,3)=0.5D0*A0
        RBAS(3,3)=-0.5D0*C0
!
!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'GRH') THEN ! TRIGONAL PRIMITIVE
        RBAS(1,1)=0. 
        RBAS(2,1)=-1.D0*A0
        RBAS(3,1)=1.D0*C0
        RBAS(1,2)=0.5D0*DSQRT(3.D0)*A0 
        RBAS(2,2)=0.5D0*A0
        RBAS(3,2)=1.0D0*C0
        RBAS(1,3)=-0.5D0*DSQRT(3.D0)*A0 
        RBAS(2,3)=0.5D0*A0
        RBAS(3,3)=1.0*C0
!
!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'GH') THEN ! HEXAGONAL PRIMITIVE
        RBAS(1,1)=0.0D0
        RBAS(2,1)=-1.D0*A0
        RBAS(3,1)=0.0D0
        RBAS(1,2)=0.5D0*DSQRT(3.D0)*A0 
        RBAS(2,2)=0.5D0*A0
        RBAS(3,2)=0.0D0
        RBAS(1,3)=0.0D0
        RBAS(2,3)=0.0D0
        RBAS(3,3)=1.D0*C0
!
!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'GC') THEN ! SIMPLE CUBIC
        RBAS(1,1)=A0  
        RBAS(2,1)=0.0D0
        RBAS(3,1)=0.0D0
        RBAS(1,2)=0.0D0
        RBAS(2,2)=A0    
        RBAS(3,2)=0.0D0
        RBAS(1,3)=0.0D0
        RBAS(2,3)=0.0D0
        RBAS(3,3)=A0
!
!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'GCF') THEN !FACE-CENTRED CUBIC
        RBAS(1,1)=0.D0  
        RBAS(2,1)=0.5D0*A0
        RBAS(3,1)=0.5D0*A0
        RBAS(1,2)=0.5D0*A0
        RBAS(2,2)=0.0D0    
        RBAS(3,2)=0.5D0*A0
        RBAS(1,3)=0.5D0*A0
        RBAS(2,3)=0.5D0*A0
        RBAS(3,3)=0.0D0
!
!     ==========================================================================
      ELSE IF(BRAVAIS.EQ.'GCV') THEN !BODY-CENTRED CUBIC
        RBAS(1,1)=-0.5D0*A0  
        RBAS(2,1)=0.5D0*A0
        RBAS(3,1)=0.5D0*A0
        RBAS(1,2)=0.5D0*A0
        RBAS(2,2)=-0.5D0*A0    
        RBAS(3,2)=0.5D0*A0
        RBAS(1,3)=0.5D0*A0
        RBAS(2,3)=0.5D0*A0
        RBAS(3,3)=-0.5D0*A0
      ELSE
        CALL ERROR$MSG('ID FOR BRAVAIS LATTICE NOT RECOGNIZED')
        CALL ERROR$CHVAL('BRAVAIS',BRAVAIS)
        CALL ERROR$STOP('SPACEGROUP$RBAS')
      ENDIF
      
      RETURN
      END SUBROUTINE SPACEGROUP$RBAS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPACEGROUP$ABCALPHABETAGAMMA(SWITCH,A,B,C,ALPHA,BETA,GAMMA,T)
!     **************************************************************************
!     **  CONVERTS THE LATTICE REPRESENTATION OF LENTHS OF AND ANGLES         **
!     **  BETWEEN LATTICE VECTORS                                             **
!     **                                                                      **
!     **  DERIVED FROM ABCALPHABETAGAMMA IN PAW_GENERALPURPOSE.F90            **
!     **************************************************************************
      IMPLICIT NONE
      LOGICAL,  INTENT(IN)    :: SWITCH  ! ABC...-> T / T->ABC...
      REAL(8),  INTENT(INOUT) :: A,B,C   ! LENGTH OF LATTICE VECTORS
      REAL(8),  INTENT(INOUT) :: ALPHA,BETA,GAMMA  ! ANGLES BTWN LATTICE VECTORS
      REAL(8),  INTENT(INOUT) :: T(3,3)  !LATTICE VECTORS  
      REAL(8)                 :: COSA,COSB,COSG,SING
!     **************************************************************************
      IF(SWITCH) THEN
        COSA=COS(ALPHA)
        COSB=COS(BETA )
        COSG=COS(GAMMA)
        SING=SQRT(1.D0-COSG**2)
        T(1,1)=A
        T(2,1)=0.D0   
        T(3,1)=0.D0   
        T(1,2)=B*COSG
        T(2,2)=B*SING
        T(3,2)=0.D0   
        T(1,3)=C*COSB
        T(2,3)=C*(COSA-COSB*COSG)/SING
        T(3,3)=C*SQRT(SING**2+2.D0*COSA*COSB*COSG-COSA**2-COSB**2)/SING
      ELSE
        A=SQRT(T(1,1)**2+T(2,1)**2+T(3,1)**2)      
        B=SQRT(T(1,2)**2+T(2,2)**2+T(3,2)**2)      
        C=SQRT(T(1,3)**2+T(2,3)**2+T(3,3)**2)      
        COSA=(T(1,2)*T(1,3)+T(2,2)*T(2,3)+T(3,2)*T(3,3))/(B*C)
        COSB=(T(1,1)*T(1,3)+T(2,1)*T(2,3)+T(3,1)*T(3,3))/(A*C)
        COSG=(T(1,1)*T(1,2)+T(2,1)*T(2,2)+T(3,1)*T(3,2))/(A*B)
        ALPHA=ACOS(COSA)
        BETA =ACOS(COSB)
        GAMMA=ACOS(COSG)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPACEGROUP$COMPLETE(NOPX,NOP,OP,C)
!     **************************************************************************
!     ** PRODUCES THE COMPLETE SET OF POINT GROUP OPERATIONS FROM A           **
!     ** SET OF GENERATORS OF THE GROUP.                                      **
!     **                                                                      **
!     ** REMARK: THE ROUTINE TOLERATES IF MORE THAN A MINIMUM SET OF          **
!     **         GENERATORS IS PROVIDED                                       **
!     ** REMARK: THE IDENTITY NEED NOT BE INCLUDED IN THE SET OF GENERATORS   **
!     **         GENERATORS IS PROVIDED                                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NOPX ! MAX #(POINT GROUP OPERATIONS)
      INTEGER(4),INTENT(INOUT) :: NOP  ! IN: #(GENERATORS) OUT:  #(OPERATIONS)
      INTEGER(4),INTENT(INOUT) :: OP(3,3,NOPX) ! IN: GENERATORS OUT: OPERATIONS
      REAL(8)   ,INTENT(INOUT) :: C(3,NOPX)    ! IN: GENERATORS OUT: OPERATIONS
      INTEGER(4)               :: GENERATOR(3,3,NOPX) ! COPY OF THE GENERATORS
      REAL(8)                  :: CGENERATOR(3,NOPX)
      INTEGER(4)               :: OPNEW(3,3)
      REAL(8)                  :: CNEW(3)
      INTEGER(4)               :: NGENERATOR
      INTEGER(4)               :: I,IOP,IOP2
      LOGICAL(4)               :: TOLD
!     **************************************************************************
      NGENERATOR=NOP
      GENERATOR(:,:,:NOP)=OP(:,:,:NOP)
      CGENERATOR(:,:NOP)=C(:,:NOP)
!
!     == IDENTITY ==============================================================
      OP(:,:,:)=0
      OP(1,1,1)=1
      OP(2,2,1)=1
      OP(3,3,1)=1
      C(:,1)=0.D0
      NOP=1
!
!     ==========================================================================
!     == COMPLETE SYMMETRY GROUP                                              ==
!     == ENSURES THAT NO NEW OPERATION IS OBTAINED BY APPLYING ANY OF THE     ==
!     == GENERATORS TO ANT OF THE OPERATIONS IN THE GROUP                     ==
!     ==========================================================================
      NOP=1
      IOP=1
      DO WHILE(IOP.LE.NOP) !LOOP OVER OPERATIONS
!WRITE(*,FMT='("OLD OPERATION ",I5,T20,3("|",3I5,"|"))')IOP,OP(:,:,IOP)
        DO I=1,NGENERATOR
          OPNEW=MATMUL(GENERATOR(:,:,I),OP(:,:,IOP))
          CNEW=MATMUL(REAL(GENERATOR(:,:,I)),C(:,IOP))+CGENERATOR(:,I)
!WRITE(*,FMT='("OPNEW ",T20,3("|",3I5,"|"))') OPNEW
!         == CHECK IF OPNEW ALREADY PRESENT IN OP ============================
          TOLD=.FALSE.
          DO IOP2=1,NOP
            TOLD=(SUM(ABS(OPNEW-OP(:,:,IOP2))).EQ.0) &
     &           .AND.(MAXVAL(ABS(CNEW-C(:,IOP2))).GT.1.D-6)
            IF(TOLD) EXIT
          ENDDO
          IF(TOLD) CYCLE
          NOP=NOP+1
          IF(NOP.GT.NOPX) THEN
            CALL ERROR$MSG('NUMBER OF OPERATIONS EXCEEDS MAXIMUM')
            CALL ERROR$I4VAL('NOP',NOP)
            CALL ERROR$I4VAL('NOPX',NOPX)
            CALL ERROR$STOP('SPACEGROUP$COMPLETE')
          END IF
          OP(:,:,NOP)=OPNEW
          C(:,NOP)=CNEW
!WRITE(*,FMT='("NEW OP: ",I5,T20,3("|",3I5,"|"))')NOP,OP(:,:,NOP)
        ENDDO ! END OF LOOP OVER GENERATORS
        IOP=IOP+1
      ENDDO !END OF LOOP OVER OPERATIONS
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$CHECKRBAS(BRAVAIS,A0,B0,C0,ALPHA,BETA,GAMMA,RBAS)
!     **************************************************************************
!     ** THIS FUNCTION CHECKS IF THE GIVEN LATTICE VECTORS RBAS ARE           **
!     ** COMPATIBLE WITH THE GIVEN TYPE OF BRAVAIS LATTICE WITHIN A TOLERANCE **
!     ** TOL AND RETURNS THE LATTICE PARAMETERS A,B,C,ALPHA,BETA,GAMMA        **
!     ** ACCORDING TO BRADLEY AND CRACKNELL                                   **
!     ******************************************************ROBERT SCHADE 2013**
      IMPLICIT NONE
      INTEGER(4),PARAMETER         :: NSYMX=48
      CHARACTER(3),INTENT(IN)      :: BRAVAIS
      REAL(8),INTENT(IN)           :: RBAS(3,3)
      REAL(8),INTENT(OUT)          :: A0,B0,C0,ALPHA,BETA,GAMMA
      INTEGER(4)                   :: ISYM,I,J,NB,IKP,IB,K
      REAL(8)                      :: RBASBRAVAIS(3,3)
      REAL(8),PARAMETER            :: TOL=1.D-6
!     **************************************************************************
      !CUBIC
      IF(BRAVAIS.EQ.'GC') THEN
        A0=RBAS(1,1)
      ENDIF
      IF(BRAVAIS.EQ.'GCF') THEN
        A0=2.0D0*RBAS(2,1)
      ENDIF
      IF(BRAVAIS.EQ.'GCV') THEN
        A0=-2.0D0*RBAS(1,1)
      ENDIF
      !MONOCLIN PRIMITIVE
      IF(BRAVAIS.EQ.'GM') THEN
        A0=SQRT(RBAS(1,2)**2+RBAS(2,2)**2)
        B0=-RBAS(2,1)
        C0=RBAS(3,3)
        GAMMA=ASIN(RBAS(2,1)/A0)
      ENDIF
      !MONOCLIN BASE-CENTRED
      IF(BRAVAIS.EQ.'GMB') THEN
        A0=SQRT((2.0D0*RBAS(2,1))**2+(2.0D0*RBAS(2,2))**2)
        B0=-RBAS(1,2)
        C0=2.0D0*RBAS(3,3)
        GAMMA=ASIN(RBAS(2,1)/(0.5D0*A0))
      ENDIF
      !ORTHORHOMBIC PRIMITIVE
      IF(BRAVAIS.EQ.'GO') THEN
        A0=RBAS(1,2)
        B0=-RBAS(2,1)
        C0=RBAS(3,3)
      ENDIF
      !ORTHORHOMBIC BASE-CENTRED
      IF(BRAVAIS.EQ.'GOB') THEN
        A0=2.0D0*RBAS(1,1)
        B0=-2.0D0*RBAS(2,1)
        C0=RBAS(3,3)
      ENDIF
      !ORTHORHOMBIC BODY-CENTRED
      IF(BRAVAIS.EQ.'GOV') THEN
        A0=2.0D0*RBAS(1,1)
        B0=2.0D0*RBAS(2,1)
        C0=2.0D0*RBAS(3,1)
      ENDIF
      !ORTHORHOMBIC FACE-CENTRED
      IF(BRAVAIS.EQ.'GOF') THEN
        A0=2.0D0*RBAS(1,1)
        B0=-2.0D0*RBAS(2,2)
        C0=2.0D0*RBAS(3,1)
      ENDIF
      !TETRAGONAL PRIMITIVE
      IF(BRAVAIS.EQ.'GQ') THEN
        A0=RBAS(1,1)
        C0=RBAS(3,3)
      ENDIF
      !TETRAGONAL BODY-CENTRED
      IF(BRAVAIS.EQ.'GQV') THEN
        A0=-2.0D0*RBAS(1,1)
        C0=2.0D0*RBAS(3,1)
      ENDIF
      !TRIGONAL PRIMITIVE
      IF(BRAVAIS.EQ.'GRH') THEN
        A0=-RBAS(2,1)
        C0=RBAS(3,1)
      ENDIF
      !HEXAGONAL PRIMITIVE
      IF(BRAVAIS.EQ.'GH') THEN
        A0=-RBAS(2,1)
        C0=RBAS(3,3)
      ENDIF
      !TRICLINIC
      IF(BRAVAIS.EQ.'G1') THEN
        RETURN !DO NOT CHECK IF TRICLINIC
      ENDIF

      CALL SPACEGROUP$RBAS(BRAVAIS,A0,B0,C0,ALPHA,BETA,GAMMA,RBASBRAVAIS)
      IF(MAXVAL(ABS(RBAS(:,:)-RBASBRAVAIS(:,:))).GT.TOL)THEN
        CALL ERROR$MSG('LATTICE VECTORS FROM PAW CALCULATION AND')
        CALL ERROR$MSG('GIVEN BRAVAIS LATTICE ARE NOT CONSISTENT')
        CALL ERROR$CHVAL("BRAVAIS TYPE",BRAVAIS)
        CALL ERROR$R8VAL("MAXIMAL ALLOWED TOLERANCE",TOL)
        CALL ERROR$R8VAL("MAX(RBAS-RBASBRAVAIS)",MAXVAL(ABS(RBAS-RBASBRAVAIS)))
        DO I=1,3
          WRITE(*,*)"RBAS       ",I,RBAS(I,:)
          WRITE(*,*)"RBASBRAVAIS",I,RBASBRAVAIS(I,:)
        ENDDO
        CALL ERROR$STOP('BRILLOUIN$CHECKRBAS')
      ENDIF
      RETURN
      END SUBROUTINE BRILLOUIN$CHECKRBAS
