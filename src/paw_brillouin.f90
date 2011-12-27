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
  integer(4)             :: NKDIV(3)
  integer(4)             :: ishift(3)
  REAL(8)       ,POINTER :: XK(:,:)   !(3,NKP) IRR. K-POINTS IN RELATIVE COORDINATES
  INTEGER(4)    ,POINTER :: IKP(:,:)  !(4,NTET) TETRAHEDRON CORNERS
  INTEGER(4)    ,POINTER :: MULT(:)   !(NTET) MULTIPLICITY OF THE TETRAHEDRON
  INTEGER(4)    ,POINTER :: irrkp(:)  !(nmshp) pointer to irr. k-point
END TYPE THIS_TYPE
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
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE BRILLOUIN_OPTIMIZE(NKP_,EKP_,WGHT_)
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: NKP_
!!$      REAL(8)    :: EKP_(NKP_)
!!$      REAL(8)    :: WGHT_(NKP_)
!!$!     **************************************************************************
!!$      CALL ERROR$MSG('BRILLOUIN_OPTIMIZE IS NOT YET TO BE USED')
!!$      CALL ERROR$STOP('BRILLOUIN_OPTIMIZE')
!!$
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN$SETR8A(ID,LENG,VAL)
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
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LENG
      REAL(8)     ,INTENT(OUT):: VAL(LENG)
      INTEGER(4)              :: I
      REAL(8)                 :: GBAS(3,3)
      REAL(8)                 :: SVAR
!     ******************************************************************
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
!       == CALCULATES WEIGHTS FOR A COMPLETELY FILLED BAND.            ==
!       == ATTENTION DEPENDS ON ACTUAL RBAS                            ==
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
      SUBROUTINE BRILLOUIN$GETI4a(ID,len,VAL)
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      integer(4)  ,intent(in) :: len
      INTEGER(4)  ,INTENT(OUT):: VAL(len)
!     **************************************************************************
      IF(ID.EQ.'NKDIV') THEN
        if(len.ne.3) then
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BRILLOUIN$GETi4A')
        END IF
        VAL=THIS%nkdiv
      ELSE
        CALL ERROR$MSG('UNKNOWN ID')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BRILLOUIN$GETI4a')
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
      real(8)                :: xk1(3)
      INTEGER(4)             :: IP(3)
      INTEGER(4)             :: I
!     **************************************************************************
      XK1(:)=MODULO(XK(:),1.D0)
      XK1(:)=XK1(:)*REAL(this%NKDIV(:),KIND=8)-0.5D0*REAL(THIS%ISHIFT(:),KIND=8)
      IP(:)=NINT(XK1(:))
      IF(SUM((XK1(:)-REAL(IP(:),KIND=8))**2).GT.1.D-5) THEN
        CALL ERROR$MSG('INVALID VALUE FOR K-POINT POSITION IN REL. COORD.')
        CALL ERROR$STOP('BRILLOUIN$IKP')
      END IF
!     __this coordinate contains all boundaries of the box!!___________________
      I=1+IP(3)+(this%nkdiv(3)+1)*(IP(2)+(this%nkdiv(2)+1)*IP(1))  ! COORDINATE
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
      SUBROUTINE BRILLOUIN$MSH(RBAS,NGKP,NSYM,IIO,IARB,tshift)
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
      logical(4),intent(inout) :: tshift     ! attempts a shift from the Gamma point
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      REAL(8)   ,ALLOCATABLE:: XK(:,:)   !(3,NGKP) K-POINTS
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
      THIS%RBAS=RBAS
      INV=0
!     ------------------------------------------------------------------
!     -- DEFINE MESH                                                  --
!     ------------------------------------------------------------------
      NMSHP=NGKP                                                         
      CALL GBASS(RBAS,GBAS,DUMMY)                                             
      CALL BRILLOUIN_BASDIV(N,NMSHP,GBAS,IARB)
      this%nkdiv(:)=n(:)
!
!     ==========================================================================
!     == ATTEMPT TO SHIFT GRID                                                ==
!     ==========================================================================
      IF(TSHIFT) THEN
        ISHIF(:)=1
        CALL BRILLOUIN_CHECKSHIFT(ISHIF,NSYM,IIO,TCHK)
        IF(.NOT.TCHK) ISHIF(:)=0
        TSHIFT=TCHK  ! return feedback whether shift was successful
      ELSE
        ISHIF(:)=0
      END IF
      THIS%ISHIFT(:)=ISHIF
!     ==========================================================================
!     == FIND IRREDUCIBLE K-POINTS                                            ==
!     ==========================================================================
      ALLOCATE(NUM(NMSHP))
      CALL BRILLOUIN_REDUZ(N,NMSHP,ISHIF,NSYM,IIO,NUM)
      IF(ISHIF(1).NE.0.OR.ISHIF(2).NE.0.OR.ISHIF(3).NE.0)INV=0   
      ALLOCATE(XK(3,NGKP))
      CALL BRILLOUIN_ZUORD(NMSHP,NUM,N,ISHIF,NGKP,NKP,XK)            
      THIS%NKP=NKP
      ALLOCATE(THIS%XK(3,NKP))
      THIS%XK(:,:)=XK(:,1:NKP)
      DEALLOCATE(XK)
!
!     ==========================================================================
!     == construct mapping onto irreducible k-points                          ==
!     ==========================================================================
!     __ the general position is encoded according to __________________________
!     __Ind=1+IP(3)+(this%nkdiv(3)+1)*(IP(2)+(this%nkdiv(2)+1)*IP(1))___________
!     __ this mapping also includes all boundaries of the unit cell! ___________
      allocate(this%irrkp(nmshp))
      this%irrkp(:)=num(:)
!
!     ==========================================================================
!     -- PRINTOUT OF SYMMETRY OPERATIONS                                      ==
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
      INTEGER(4),INTENT(IN) :: NKDIV(3)   ! DIVISIONS OF RECIPROCAL LATTICE VECTORS
      INTEGER(4),INTENT(IN) :: ISHIFT(3)  ! GRID IS SHIFTED 
      LOGICAL(4),PARAMETER  :: TPR=.false.
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
!     ******************************************************************
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
      this%ishift(:)=ishift(:)
      this%nkdiv(:)=nkdiv(:)
!     ------------------------------------------------------------------
!     -- DEFINE MESH                                                  --
!     ------------------------------------------------------------------
      NMSHP=(NKDIV(1)+1)*(NKDIV(2)+1)*(NKDIV(3)+1)
      CALL GBASS(RBAS,GBAS,DUMMY)                                             
!     ------------------------------------------------------------------
!     -- FIND IRREDUCIBLE K-POINTS                                    --
!     ------------------------------------------------------------------
      ALLOCATE(NUM(NMSHP))
      CALL BRILLOUIN_REDUZ(NKDIV,NMSHP,ISHIFT,NSYM,IIO,NUM)
!     == INV=1 EXPLOITS TIME INVERSION SYMMETRY FOR UNSHIFTED GRIDS  ==
!     == WHEN SEARCHING IRREDUCIBLE TETRAHEDRA. PROBABLY IRRELEVANT  ==
      IF(ISHIFT(1).NE.0.OR.ISHIFT(2).NE.0.OR.ISHIFT(3).NE.0)INV=0   
      IF(.NOT.TINV) INV=0
!
!     ==================================================================
!     ==  DETERMINE IRREDUCIBLE K-POINTS                              ==
!     ==================================================================
      ALLOCATE(XK(3,NGKP))
      CALL BRILLOUIN_ZUORD(NMSHP,NUM,NKDIV,ISHIFT,NGKP,NKP,XK)            
      THIS%NKP=NKP
      ALLOCATE(THIS%XK(3,NKP))
      THIS%XK(:,:)=XK(:,1:NKP)
      DEALLOCATE(XK)
!
!     ==========================================================================
!     == construct mapping onto irreducible k-points                          ==
!     ==========================================================================
!     __ the general position is encoded according to __________________________
!     __Ind=1+IP(3)+(this%nkdiv(3)+1)*(IP(2)+(this%nkdiv(2)+1)*IP(1))___________  
!     __ this mapping also includes all boundaries of the unit cell! ___________
      allocate(this%irrkp(nmshp))
      this%irrkp(:)=num(:)
!
!     ------------------------------------------------------------------
!     -- PRINTOUT OF SYMMETRY OPERATIONS                              --
!     ------------------------------------------------------------------
      IF(TPR) THEN
        DO ISYM=1,NSYM                                             
          WRITE(*,FMT='("SYMMETRYMATRIX NR. : ",I5/3(" ",3I10/))') &
     &            ISYM,((IIO(I,J,ISYM),J=1,3),I=1,3)                 
        ENDDO
      END IF                                                            
!     ----------------------------------------------------------------- 
!     --  PRINTOUT OF MAPPING TO IRREDUCIBLE K-POINTS                -- 
!     ----------------------------------------------------------------- 
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
!     ----------------------------------------------------------------- 
!     --  PRINT IREDUCIBLE K-POINTS                                  -- 
!     ----------------------------------------------------------------- 
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
!     ----------------------------------------------------------------- 
!     --  CHOOSE TETRAHEDRA                                          -- 
!     ----------------------------------------------------------------- 
      CALL BRILLOUIN_TETDIV(NKDIV,GBAS,TET0)                                      
      IF(TPR) THEN                                                 
        WRITE(*,FMT='(3(I2,"-TES NORMTETRAHEDRON "))')1,2,3                      
        WRITE(*,FMT='(4I5,5X,4I5,5X,4I5)') &
     &         (((TET0(J,K,I),K=1,4),I=1,3),J=1,3)    
        WRITE(*,FMT='(3(I2,"-TES NORMTETRAHEDRON "))')4,5,6
        WRITE(*,FMT='(4I5,5X,4I5,5X,4I5)') &
     &         (((TET0(J,K,I),K=1,4),I=4,6),J=1,3)    
      END IF                                                            
!     ----------------------------------------------------------------- 
!     --  FIND INEQUIVALENT TETRAHEDRA                               -- 
!     ----------------------------------------------------------------- 
      CALL BRILLOUIN_TETCNT(NMSHP,NUM,TET0,NKDIV,INV,THIS)
      DEALLOCATE(NUM)
                                   CALL TRACE$POP
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_BASDIV(N,NMSHP,GBAS,IARB)                              
!     **                                                              **
!     **  BASDIV DETERMINES DIVISION OF BASEVECTORS OF REC. LATT.     **
!     **  SO THAT THE NUMBER OF MESHPOINTS IS JUST BELOW NMSHP AND    **
!     **  TAKES INTO ACCOUNT THE DEPENDENCY BY POINTSYMMETRY          **
!     **  INPUT :                                                     **
!     **    GBAS        RECIPROCAL LATTICE VECTORS                    **
!     **    IARB        DEPENDENCIES FOR DIVISIONS                    **
!     **                OF RECIPROCAL LATTICE VECTORS                 **
!     **                ( IF IARB(1)=1 THEN 1ST AND 2ND LATTICE       **
!     **                  VECTORS ARE DIVIDED BY AN EQUAL NUMBER;     **
!     **                  IF IARB(2)=1 THEN SAME FOR 2ND AND 3RD;     **
!     **                  IF IARB(3)=1 THEN SAME FOR 3RD AND 1ST)     **
!     **    NMSHP       TOTAL NUMBER OF GRID POINTS                   **
!     **  OUTPUT:                                                     **
!     **    N           NUMBER OF DIVISIONS OF RECIPROCAL LATTICE     **
!     **                VECTORS FOR DEFINITION OF SUBLATTICE          **
!     **    NMSHP       NUMBER OF SUBLATTICE POINTS IN A REC. UNIT CELL*
!     **                AND ON ALL ITS FACES                          **
!     **                                                              **
!     ******************************************************************
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
!     ******************************************************************
!     ==================================================================
!     == FIND UPPER LIMIT FOR THE LENGTH OF SUBLATTICE VECTORS        ==
!     ==================================================================
      DO I=1,3                                                       
        BETR(I)=SQRT(GBAS(I,1)**2+GBAS(I,2)**2+GBAS(I,3)**2)             
      ENDDO
      SVAR=(REAL(NMSHP,KIND=8)/(BETR(1)*BETR(2)*BETR(3)))**(1.D0/3.D0)         
      DO I=1,3                                                       
        RN(I)=BETR(I)*SVAR                                                
      ENDDO
!     ==================================================================
!     == FIND DIVISIONS OF LATTICE VECTORS                            ==
!     ==================================================================
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
!     ==================================================================
!     == USE THIS TO FIX THE K-MESH PY HAND                           ==
!     ==================================================================
!     PRINT*,' K-MESH PUT IN BY HAND!!!!!!!!!!'                         
!     N(1)=12                                                           
!     N(2)=12                                                           
!     N(3)=2                                                            
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      NMSHP=(N(1)+1)*(N(2)+1)*(N(3)+1)                                  
      WRITE(*,FMT='("NO. OF MESH POINTS IN THE BRILLOUIN ZONE =",I6)') &
     &            N(1)*N(2)*N(3)      
      WRITE(*,FMT='("DIVISION OF RECIPROCAL LATTICE VECTORS (INTERVALS)=" &
     &        ,3I5)')N(:)
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_REDUZ(N,NMSHP,ISHIFT,NSYM,IO,NUM)
!     **                                                              **
!     **  REDUZ CREATES THE RELATION BETWEEN THE MESHPOINTS AND       **
!     **  THE POINTS IN THE "IRREDUCIBLE ZONE"                        **
!     **  INPUT :                                                     **
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    NMSHP       NUMBER OF SUBLATTICE POINTS INSIDE AND        **
!     **                ON ALL FACES OF A REC. UNIT CELL              **
!     **    NSYM        NUMBER OF SYMMETRY OPERATIONS                 **
!     **    IO          SYMMETRY OPERATIONS (BASIS ARE REC. LATT. VEC.)*
!     **  OUTPUT :                                                    **
!     **    NUM(I)      MAPPING FROM A GENERAL POINT (I) TO THE       **
!     **                CORRESPONDING IRREDUCIBLE POINT (NUM)         **
!     **  REMARKS :                                                   **
!     **    THE MAPPING FROM COORDINATES TO NUMBERS IS GIVEN BY :     **
!     **   (X,Y,Z)=RBAS*(I,J,K)                                       **
!     **   (I,J,K)  <->  NUM = I*(N(2)+1)*(N(3)+1)+J*(N(3)+1)+K+1     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N(3)
      INTEGER(4),INTENT(IN) :: NMSHP
      INTEGER(4),INTENT(IN) :: ISHIFT(3)
      INTEGER(4),INTENT(IN) :: NSYM
      INTEGER(4),INTENT(IN) :: IO(3,3,NSYM)
      INTEGER(4),INTENT(OUT):: NUM(NMSHP)
      LOGICAL(4),PARAMETER  :: TSHIFT=.TRUE.
      INTEGER(4)            :: I,J
      INTEGER(4)            :: I1,I2,I3,J1,J2,J3
      INTEGER(4)            :: IX1,IX2,IX3,IND,IND1,IND2
      LOGICAL(4)            :: TCHK
!     ******************************************************************
                           CALL TRACE$PUSH('BRILLOUIN_REDUZ')
      CALL BRILLOUIN_CHECKSHIFT(ISHIFT,NSYM,IO,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('SHIFT INCOMPATIBLE WITH SYMMETRY')
        CALL ERROR$STOP('BRILLOUIN_REDUZ')
      END IF
!     ==================================================================
!     ==  INITIALIZE                                                  ==
!     ==================================================================
      DO I=1,NMSHP                                                   
        IND=I-1    ! I1 RUNS FASTEST
        I1=IND/((N(3)+1)*(N(2)+1))      
        IND=IND-I1*(N(3)+1)*(N(2)+1)
        I2=IND/(N(3)+1)
        I3=IND-I2*(N(3)+1)
        J1=MODULO(I1,N(1))
        J2=MODULO(I2,N(2))
        J3=MODULO(I3,N(3))
        NUM(I)=1+J3+(N(3)+1)*(J2+(N(2)+1)*J1)
      ENDDO
!     ==================================================================
!     ==  REDUCTION OF REDUCIBLE K-POINTS                             ==
!     ==================================================================
      DO I=1,NSYM                                                    
        DO J=1,NMSHP                                                   
          IND=NUM(J)-1
          I1=IND/((N(3)+1)*(N(2)+1))
          IND=IND-I1*(N(3)+1)*(N(2)+1)
          I2=IND/(N(3)+1)
          I3=IND-I2*(N(3)+1)
          IND1=1+I1+N(1)*(I2+N(2)*I3)
!         == PLACE ON A GRID THAT IS DOUBLE AS FINE
          I1=2*I1+ISHIFT(1)                                            
          I2=2*I2+ISHIFT(2)                                            
          I3=2*I3+ISHIFT(3)                                            
!         == APPLY SYMMETRY ============================================
          J1=IO(1,1,I)*I1+IO(1,2,I)*I2+IO(1,3,I)*I3             
          J2=IO(2,1,I)*I1+IO(2,2,I)*I2+IO(2,3,I)*I3             
          J3=IO(3,1,I)*I1+IO(3,2,I)*I2+IO(3,3,I)*I3             
!         == CHECK IF THE POINT FALLS ONTO THE ORIGINAL SHIFTED GRID ===
          IX1=MODULO(J1-ISHIFT(1),2)                                           
          IX2=MODULO(J2-ISHIFT(2),2)                                           
          IX3=MODULO(J3-ISHIFT(3),2)                                           
          IF(IX1.NE.0.OR.IX2.NE.0.OR.IX3.NE.0) CYCLE  
!         == COARSEN GRID AGAIN ======================================
          J1=(J1-ISHIFT(1))/2                                               
          J2=(J2-ISHIFT(2))/2                                               
          J3=(J3-ISHIFT(3))/2                                               
!         == DO NOT CONSIDER POINTS AT THE OUTER BOUNDARY AS 
          J1=MODULO(J1,N(1))
          J2=MODULO(J2,N(2))
          J3=MODULO(J3,N(3))
          IND2=1+J1+N(1)*(J2+N(2)*J3)
!         == IDENTIFY WITH ONE OF THE TWO SYMMETRY RELATED K-POINTS
          IF(IND2.LT.IND1) THEN
            NUM(J)=1+J3+(N(3)+1)*(J2+(N(2)+1)*J1)
          END IF
        ENDDO
      ENDDO
!
!     ==================================================================
!     == THE FOLLOWING IS NEEDED BECAUSE THE OPERATIONS MAY NOT COMMUTE?
!     == A STAR IS CHARACTERIZED BY HAVING THE SAME VALUE OF NUM(I)    
!     ==================================================================
      DO I=1,NMSHP
        IF(NUM(I).NE.NUM(NUM(I))) THEN
          CALL ERROR$MSG('MAPPING FAILED')
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
!     ******************************************************************
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
            I=I3+(N(3)+1)*(I2-1+(N(2)+1)*(I1-1))  ! coordinate
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
              XK(1,NDIM)=(REAL(I1-1,KIND=8)+REAL(ISHIFT(1),KIND=8)/2.D0)/REAL(N(1),KIND=8)              
              XK(2,NDIM)=(REAL(I2-1,KIND=8)+REAL(ISHIFT(2),KIND=8)/2.D0)/REAL(N(2),KIND=8)              
              XK(3,NDIM)=(REAL(I3-1,KIND=8)+REAL(ISHIFT(3),KIND=8)/2.D0)/REAL(N(3),KIND=8)              
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
!     ******************************************************************
!     **                                                             ** 
!     **  TETDIV DETERMINES THE DIVISION OF THE PARALLELEPIPEDS      ** 
!     **  IN TETRAHEDRONS ACCORDING TO THE SHORTEST DIAGONAL         ** 
!     **                                                             ** 
!     **  INPUT:                                                      **
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    GBAS        RECIPROCAL LATTICE VECTORS                    **
!     **  OUTPUT:                                                    ** 
!     **    TET0(I,J,K) COORDINATES (I) IN THE BASIS OF SUBLATTICE   ** 
!     **                VECTORS, OF THE 4 CORNERS (J) OF A           ** 
!     **                TETRAHEDRON (K) IN A SUBLATTICE UNIT CELL    ** 
!     **                                                             ** 
!     ******************************************************************
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
!     ******************************************************************
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
!     **                                                              **
!     **  TETCNT CALCULATES ALL DIFFERENT TETRAHEDRA AND COUNTS THEM  **
!     **  INPUT :                                                     **
!     **    NMSHP       NUMBER OF SUBLATTICE POINTS INSIDE AND        **
!     **                ON ALL FACES OF A REC. UNIT CELL              **
!     **    NUM(I)      MAPPING FROM A GENERAL POINT (I) TO THE       **
!     **                CORRESPONDING IRREDUCIBLE POINT (NUM)         **
!     **    TET0(I,J,K) COORDINATES (I) IN THE BASIS OF SUBLATTICE    ** 
!     **                VECTORS, OF THE 4 CORNERS (J) OF A            ** 
!     **                TETRAHEDRON (K) IN A SUBLATTICE UNIT CELL     ** 
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    INV         DUMMY NUMBER (MUST BE 0)                      **
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                **
!     **    MWRIT       INFORMATION FOR MWRIT TETRAHEDRA ARE WRITTEN  **
!     **                AT ONE TIME.                                  **
!     **                                                              **
!     ******************************************************************
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
!     **                                                              **
!     **  CALCULATES THE SAMPLING WEIGHTS FROM TETRAHEDRON INTEGRATION**
!     **                                                              **
!     **  INPUT :                                                     **
!     **    NB          NUMBER OF BANDS                               **
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                **
!     **    EB          ENERGIES ( DETERMINE FERMI SURFACE )          **
!     **    RNTOT       NUMBER OF OCCUPIED STATES                     **
!     **  OUTPUT :                                                    **
!     **    WGHT        SAMPLING WEIGHTS                              **
!     **    EF          FERMI LEVEL                                   **
!     **                                                              **
!     **  AUTHOR : PETER E. BLOECHL                                   **
!     **                                                              **
!     **  SUBROUTINES USED:                                           **
!     **  EFERMI,SAMFAC,TOTNOS,EFI,WEIGHT                             **
!     **                                                              **
!     ******************************************************************
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NB          ! #(BANDS)
      INTEGER(4),INTENT(IN) :: NKP         ! #(IRREDUCIBLE K-POINTS)
      REAL(8)   ,INTENT(IN) :: EB(NB*NKP)  ! ENERGY BANDS
      REAL(8)   ,INTENT(OUT):: WGHT(NB*NKP)! INTEGRATION WEIGHTS
      REAL(8)   ,INTENT(IN) :: RNTOT       ! #(ELECTRONS)
      REAL(8)   ,INTENT(OUT):: EF          ! FERMI LEVEL
      REAL(8)   ,PARAMETER  :: TOLMAX=1.D-5
      LOGICAL(4)            :: TCHECK=.TRUE.      
!     ******************************************************************
                       CALL TRACE$PUSH('BRILLOUIN$DOS')
!
!     ------------------------------------------------------------------
!     --  CALCULATE FERMI LEVEL                                       --
!     ------------------------------------------------------------------
      CALL BRILLOUIN_EFERMI(RNTOT,EF,TOLMAX,NB,NKP,EB,THIS)
                                                                        
      PRINT*,'FERMI ENERGY AT ',EF                                      
!     ------------------------------------------------------------------
!     --  CALCULATE WEIGHTS                                           --
!     ------------------------------------------------------------------
      CALL BRILLOUIN_SAMFAC(NB,NKP,EB,EF,WGHT,THIS)
                                                                        
!     ------------------------------------------------------------------
!     --  CHECK WHETHER SUMRULE IS FULLFILLED                         --
!     ------------------------------------------------------------------
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
      SUBROUTINE BRILLOUIN$wDOS(NB,NKP,EB,ef,WGHT)
!     **************************************************************************
!     **  CALCULATES THE SAMPLING WEIGHTS for the density of states at        **
!     **  energy ef                                                           **
!     **                                                                      **
!     **  INPUT :                                                             **
!     **    NB          NUMBER OF BANDS                                       **
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                        **
!     **    EB          ENERGIES ( DETERMINE FERMI SURFACE )                  **
!     **    ne          #(points in energy grid)                              **
!     **    ei          grid energies                                         **
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
      REAL(8)   ,INTENT(IN) :: ef
      REAL(8)   ,INTENT(OUT):: WGHT(NB,NKP)! dos WEIGHTS
!     **************************************************************************
                       CALL TRACE$PUSH('BRILLOUIN$wDOS')
!
!     ------------------------------------------------------------------
!     --  CALCULATE WEIGHTS                                           --
!     ------------------------------------------------------------------
      CALL BRILLOUIN_SAMdos(NB,NKP,EB,Ef,WGHT,THIS)
                               CALL TRACE$POP()
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_EFERMI(RNTOT,EF,TOLMAX,NB,NKP,EB,THIS)
!     **                                                              **
!     **  CALCUALTES THE FERMILEVEL BY INTEGRATION OF THE TOTAL       **
!     **  DENSITY OF STATES                                           **
!     **                                                              **
!     **  INPUT :                                                     **
!     **    RNTOT       NUMBER OF OCCUPIED STATES                     **
!     **    NB          NUMBER OF BANDS                               **
!     **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                **
!     **    EB          ENERGIES ( DETERMINE FERMI SURFACE )          **
!     **    TOLMAX      TOLERANCE IN THE NUMBER OF STATES AT EF       **
!     **  OUTPUT :                                                    **
!     **    EF          FERMI LEVEL                                   **
!     **                                                              **
!     ******************************************************************
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
!     ******************************************************************
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
        PRINT*,'ILOOP ',ILOOP                                           
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
!     ******************************************************************
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
!     **                                                              **
!     **  CALCULATES THE INTEGRATED DOS                               **
!     **  FOR ONE TETRAHEDRON ON THE ENERGYMESH                       **
!     **  INPUT :                                                     **
!     **    VOL         WEIGHT OF THIS TETRAHEDRON                    **
!     **    E           ENERGIES AT THE EDGEPOINTS                    **
!     **    EMIN        MINIMUM VALUE OF ENERGY MESH FOR NOS AND SOS  **
!     **    EMAX        MAXIMUM VALUE OF ENERGY MESH FOR NOS AND SOS  **
!     **    NP          NUMBER OF POINTS ON THE ENERGY MESH           **
!     **  OUTPUT:                                                     **
!     **    SOS         > NOS(E)+ SUM OVER E: SOS(E) =                **
!     **    NOS         > NUMBER OF STATES BELOW E                    **
!     **                                                              **
!     **    SOS(I) CONTAINS THE WEIGHTS OFF ALL TETRAHEDRA WITH THE   **
!     **    HIGHEST ENERGY EX IN THE WINDOW E(I-1)<EX<E(I).           **
!     **    SOS(1) CONTAINS THE WEIGHTS OF ALL TETRAHEDRA WITH THEIR  **
!     **    HIGHEST ENERGY EX BELOW E(I)                              **
!     **                                                              **
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
!     ******************************************************************
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
      REAL(8)                    :: WGHT0(4),dwght0(4)
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
          call BRILLOUIN_DWEIGHT(VOL,E,EF,WGHT0,DWGHT0)
          DO I=1,4                                                      
            WGHT(IB,IKP(I))=WGHT(IB,IKP(I))+dWGHT0(I)                          
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
!     **                                                                      **
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
!       ------  TETRAHEDRON X1,X2,X13',X14'                             
        VPRIME=VOL*DE1**2/(E41*E31)*.25D0                               
        WGHT(2)=VPRIME                                                  
        WGHT(3)=VPRIME*(DE1/E31)                                        
        WGHT(4)=VPRIME*(DE1/E41)                                        
        WGHT(1)=VPRIME*(3.D0-DE1/E41-DE1/E31)                           
!       ------  TETRAHEDRON X2,X13',X23',X14'                           
        VPRIME=.25D0*VOL*DE2*DE3*DE1/(E32*E31*E41)                      
        WGHT(1)=WGHT(1)+VPRIME*(2.D0-DE1/E31-DE1/E41)                   
        WGHT(2)=WGHT(2)+VPRIME*(2.D0-DE2/E32)                           
        WGHT(3)=WGHT(3)+VPRIME*(DE2/E32+DE1/E31)                        
        WGHT(4)=WGHT(4)+VPRIME*(DE1/E41)                                
!       ------  TETRAHEDRON X2,X23',X24',X14'                           
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
      REAL(8)  ,INTENT(OUT):: dWGHT(4)  ! energy derivative of i.W.
      INTEGER(4),PARAMETER :: ICOR=0   ! ON/OFF SWITCH FOR CORRECTION 
      REAL(8)              :: E(4)     ! ENERGY BANDS AT THE CORNERS
      REAL(8)              :: FA(4),dFA(4)
      REAL(8)              :: FB(4)
      INTEGER(4)           :: INDEX(4)
      REAL(8)              :: X
      REAL(8)              :: VPRIME,dvprime
      REAL(8)              :: DE
      REAL(8)              :: DOS,ddos
      REAL(8)              :: E21,E31,E41,E32,E42,E43
      REAL(8)              :: DE1,DE2,DE3,DE4
      INTEGER(4)           :: I,J,IP,N,M,K
      REAL(8)              :: DA,DB,DC
      REAL(8)              :: VOL14
!     **************************************************************************
      E(:)=E_(:)
      WGHT(:)=0.D0
      dWGHT(:)=0.D0
!     ==========================================================================
!     ==  INTEGRATION WITHOUT FERMISURFACE                                    ==
!     ==========================================================================
      IF(MINVAL(E).GE.EF) RETURN
      IF(MAXVAL(E).LE.EF) THEN                                                  
        VPRIME=.25D0*VOL                                                
        WGHT(:)=VPRIME                                                  
        dWGHT(:)=0.d0
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
        dVPRIME=.25D0*VOL*3.d0*DE**2/(E21*E31*E41)                            
        WGHT(1)=VPRIME*(4.D0-DE/E21-DE/E31-DE/E41)                      
        WGHT(2)=VPRIME*DE/E21                                           
        WGHT(3)=VPRIME*DE/E31                                           
        WGHT(4)=VPRIME*DE/E41                                           
        dWGHT(1)=dVPRIME*(4.D0-DE/E21-DE/E31-DE/E41) &
       &                     +VPRIME*(-1.d0/E21-1.d0/E31-1.d0/E41)
        dWGHT(2)=dVPRIME*DE/E21+VPRIME/E21
        dWGHT(3)=dVPRIME*DE/E31+VPRIME/E31
        dWGHT(4)=dVPRIME*DE/E41+VPRIME/E41
!       ------  PARAMETERS FOR CORRECION                                
        DOS=3.D0*VPRIME*4.D0/(EF-E(1))
        ddos=-dos/(EF-E(1))
      ELSE IF(EF.GT.E(2).AND.EF.LT.E(3)) THEN                           
        DE1=EF-E(1)                                                     
        DE2=EF-E(2)                                                     
        DE3=E(3)-EF                                                     
        DE4=E(4)-EF                                                     
!       ------  TETRAHEDRON X1,X2,X13',X14'                             
        VPRIME=VOL*DE1**2/(E41*E31)*.25D0                               
        dVPRIME=VOL*2.d0*DE1/(E41*E31)*.25D0                               
        WGHT(2)=VPRIME                                                  
        WGHT(3)=VPRIME*(DE1/E31)                                        
        WGHT(4)=VPRIME*(DE1/E41)                                        
        WGHT(1)=VPRIME*(3.D0-DE1/E41-DE1/E31)                           
        dWGHT(2)=dVPRIME                                                  
        dWGHT(3)=dVPRIME*(DE1/E31)+vprime/e31
        dWGHT(4)=dVPRIME*(DE1/E41)+vprime/e41
        dWGHT(1)=dVPRIME*(3.D0-DE1/E41-DE1/E31)+VPRIME*(-1.d0/E41-1.d0/E31)
!       ------  TETRAHEDRON X2,X13',X23',X14'                           
        VPRIME=.25D0*VOL*DE2*DE3*DE1/(E32*E31*E41)                      
        dVPRIME=.25D0*VOL*(DE3*DE1-DE2*DE1+DE2*DE3)/(E32*E31*E41)
        WGHT(1)=WGHT(1)+VPRIME*(2.D0-DE1/E31-DE1/E41)                   
        WGHT(2)=WGHT(2)+VPRIME*(2.D0-DE2/E32)                           
        WGHT(3)=WGHT(3)+VPRIME*(DE2/E32+DE1/E31)                        
        WGHT(4)=WGHT(4)+VPRIME*(DE1/E41)                                
        dWGHT(1)=dWGHT(1)+dVPRIME*(2.D0-DE1/E31-DE1/E41) &
       &                                           +VPRIME*(-1.d0/E31-1.d0/E41)
        dWGHT(2)=dWGHT(2)+dVPRIME*(2.D0-DE2/E32)   +VPRIME*(-1.d0/E32)
        dWGHT(3)=dWGHT(3)+dVPRIME*(DE2/E32+DE1/E31)+VPRIME*(1.d0/E32+1.d0/E31)
        dWGHT(4)=dWGHT(4)+dVPRIME*(DE1/E41)        +VPRIME*(1.d0/E41)
!       ------  TETRAHEDRON X2,X23',X24',X14'                           
        VPRIME=.25D0*VOL*DE2**2*DE4/(E42*E32*E41)                       
        dVPRIME=.25D0*VOL*(2.d0*DE2*DE4-DE2**2)/(E42*E32*E41)                       
        WGHT(1)=WGHT(1)+VPRIME*(1.D0-DE1/E41)                           
        WGHT(2)=WGHT(2)+VPRIME*(3.D0-DE2/E32-DE2/E42)                   
        WGHT(3)=WGHT(3)+VPRIME*(DE2/E32)                                
        WGHT(4)=WGHT(4)+VPRIME*(DE2/E42+DE1/E41)                        
        dWGHT(1)=dWGHT(1)+dVPRIME*(1.D0-DE1/E41)   +VPRIME*(-1.d0/E41)
        dWGHT(2)=dWGHT(2)+dVPRIME*(3.D0-DE2/E32-DE2/E42) &
       &                                           +VPRIME*(-1.d0/E32-1.d0/E42)
        dWGHT(3)=dWGHT(3)+dVPRIME*(DE2/E32)        +VPRIME*(1.d0/E32)
        dWGHT(4)=dWGHT(4)+dVPRIME*(DE2/E42+DE1/E41)+VPRIME*(1.d0/E42+1.d0/E41)
!       ------  DOS=A+B*(EF-E2)+C*(EF-E2)**2                            
        DA=3.D0*VOL*E21/(E31*E41)                                       
        DB=6.D0*VOL/(E31*E41)                                           
        DC=-3.D0*VOL/(E32*E41*E31*E42)*(E31+E42)                        
        DOS=Da+DB*de2+DC*DE2**2
        dDOS=DB+2.d0*DC*DE2
      ELSE IF(EF.GE.E(3).AND.EF.LT.E(4)) THEN                           
        DE=E(4)-EF                                                      
        VPRIME=.25D0*VOL*DE**3/(E41*E42*E43)                            
        dVPRIME=.25D0*VOL*(-3.d0*DE**2)/(E41*E42*E43)                            
        VOL14=.25D0*VOL                                                 
        WGHT(1)=VOL14-VPRIME*DE/E41                                     
        WGHT(2)=VOL14-VPRIME*DE/E42                                     
        WGHT(3)=VOL14-VPRIME*DE/E43                                     
        WGHT(4)=VOL14-VPRIME*(4.D0-DE/E41-DE/E42-DE/E43)                
        dWGHT(1)=-dVPRIME*DE/E41+VPRIME/E41
        dWGHT(2)=-dVPRIME*DE/E42+VPRIME/E42
        dWGHT(3)=-dVPRIME*DE/E43+VPRIME/E43                                     
        dWGHT(4)=-dVPRIME*(4.D0-DE/E41-DE/E42-DE/E43) &
       &                      -VPRIME*(1.d0/E41+1.d0/E42+1.d0/E43)
!       ------  PARAMETERS FOR CORRECION                                
        DOS=3.D0*VPRIME*4.D0/(E(4)-EF)                                  
        dDOS=dos/(E(4)-EF)                                  
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
            dWGHT(M)=dWGHT(M)+.25D0*(E(N)-E(M))*dDOS*.1D0                      
          ENDDO
        ENDDO
      END IF                                                            
!     ==========================================================================
!     ==  REORDER WEIGHTS                                                     ==
!     ==========================================================================
      DO I=1,4                                                      
        FA(INDEX(I))=WGHT(I)                                              
        dFA(INDEX(I))=dWGHT(I)                                              
        FB(INDEX(I))=E(I)                                                 
      ENDDO
      DO I=1,4                                                      
        WGHT(I)=FA(I)                                                     
        dWGHT(I)=dFA(I)                                                     
        E(I)=FB(I)                                                        
      ENDDO
      RETURN                                                            
      END                                                               
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BRILLOUIN_WEIGHTANDDER(VOL,E_,EF,WGHT,DWGHT)                                  
!     **************************************************************************
!     **                                                                      **
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
!       ------  TETRAHEDRON X1,X2,X13',X14'                             
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
!       ------  TETRAHEDRON X2,X13',X23',X14'                           
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
!       ------  TETRAHEDRON X2,X23',X24',X14'                           
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
