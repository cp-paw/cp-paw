!
!.......................................................................
MODULE BRILLOUIN_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: BRILLOUIN                                                  **
!**                                                                   **
!**  PURPOSE: BRILLOUIN ZONE INTEGRATION USING THE IMPROVED VERSION   **
!**    OF THE TETRAHEDRON METHOD                                      **
!**                                                                   **
!**  DESCRIPTION:                                                     **
!**    BRILLOUIN$MSH CREATE INFORMATION WHICH IS ONLY DEPENDENT       **
!**    ON THE STRUCTURE OF THE CRYSTAL. THIS FUNCTION IS EXECUTED     **
!**    BEFORE THE SELFCONSITENCY ITERATIONS. THE IRREDUCIBLE K-POINTS **
!**    ARE THEN AVAILABLE THROUGH THE $GET.. INTERFACE                ** 
!**                                                                   **
!**    BRILLOUIN$DOS IS EXECUTED IN EACH SELFCONSISTENCY ITERATION.   **
!**    THEY NEED AS INPUT THE ENERGY EIGENVALUES EB AND THE NUMBER OF **
!**    OCCUPIED STATES RNTOT. FOR SPIN POLARIZED CALCULATIONS TREAT   **
!**    SPIN UP AND SPIN DOWN AS SEPARATE BANDS. OUTPUT ARE THE        **
!**    WEIGHTS WGHT. THE INTEGRATION OF AN ARBITRARY FUNCTION A(N,K)  **
!**    OVER THE OCCUPIED STATES IS PERFORMED BY SUMMATION             **
!**            <A>=SUM OVER K AND N OF WGHT(N,K)*A(N,K)               **
!**    WHERE THE SUM RUNS OVER OCCUPIED AND! UNOCCUPIED STATES.       **
!**    THE WEIGHTS CONTAIN BOTH THE GEOMETRICAL WEIGHT AND THE        **
!**    INFLUENCE OF THE FERMI FUNCTION.                               **
!**                                                                   **
!**    A) FOR INSULATORS THE METHOD IS IDENTICAL TO THE SPECIAL POINT **
!**       SCHEME OF MONKHORST AND PACK.                               **
!**    B) FOR METALS IT IS IDENTICAL TO THE "TRADITIONAL" TETRAHEDRON **
!**       METHOD OF ANDERSEN AND JEPSEN, IF THE CORRECTION IS SWITCHED**
!**       OFF.                                                        **
!**    C) WITH THE CORRECTION ALSO THE RESULTS FOR METALS ARE         **
!**       COMPARABLE TO THAT OBTAINED FOR INSULATORS                  **
!**                                                                   **
!**    THE CORRECTION FORMULA FOR LINEAR INTERPOLATION CAN BE SWITCHED**
!**    ON AND OFF BY THE PARAMETER "ICOR" IN THE SUBROUTINE "WEIGHTS" **
!**                                                                   **
!**    THE SYMMETRY OPERATIONS USED AS INPUT ARE TABULATED IN:        **
!**      C.J.BRADLEY AND A.P.CRACKNELL,                               **
!**      THE MATHEMATICAL THEORY OF SYMMETRY IN SOLIDS,               **
!**      OXFORD 1972                                                  **
!**                                                                   **
!**    SOMETIMES THE ROUTINE HAS PROBLEMS FINDING THE FERMI LEVEL,    **
!**    IF THE FERMI LEVEL IS PINNED A TETRAHEDRON WITH IDENTICAL      **
!**    ENERGIES ON ALL 4 CORNERS, WHICH LEADS TO A DELTA PEAK IN THE  **
!**    DENSITIES OF STATES. IN THIS CASE ONE SHOULD INCREASE THE K-MESH,
!**    OR CHANGE THE K-MESH FROM AN EVEN NUMBER OF DIVISIONS FOR THE  ** 
!**    RECIPROCAL LATTICE VECTORS TO AN ODD NUMBER OR VICE VERSA.     **
!**                                                                   **
!**    THE FERMI LEVEL IS DETERMINED TO AN ACCURACY OF 1.E-5 IN THE   **
!**    NUMBER OF STATES. THIS TOLERANCE CAN BE MODIFIED BY CHANGING   **
!**    THE SETTING OF THE PARAMETER "TOLMAX" IN SUBROUTINE DOS.       **
!**                                                                   **
!**    THE GRID OF K-POINTS MAY OR MAY NOT INCLUDE THE GAMMA POINT    **
!**    SWITCH USING TSHIFT IN BRILLOUIN_REDUZ                         **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    BRILLOUIN$MSH                                                  **
!**    BRILLOUIN$DOS                                                  **
!**    BRILLOUIN$GETI4                                                **
!**    BRILLOUIN$GETR8A                                               **
!**                                                                   **
!**  USAGE:                                                           **
!**    1) DEFINE LATTICE AND SYMMETRY OPERATIONS USING BRILLOUIN$MSH. **
!**    2) OBTAIN THE IRREDUCIBLE K-POINTS USING                       **
!**           BRILLOUIN$GETI4('NK',NK)                                **
!**           BRILLOIIN$GETR8A('K',3*NK,K)                            **
!**    3) CALCULATE ONE-PARTICLE ENERGIES AT THE K-POINTS             **
!**    4) OBTAIN SAMPLING WEIGHTS USING BRILLOUIN$DOS                 **
!**                                                                   **
!***********************************************************************
TYPE THIS_TYPE
  INTEGER(4)             :: NKP      ! #(IRREDUCIBLE KPOINTS)
  INTEGER(4)             :: NTET     ! #(IRREDUCIBLE TETRAHEDRA)
  REAL(8)                :: VOL      ! 1/#(GENERAL TETRAHEDRA)
  REAL(8)                :: RBAS(3,3)! REAL SPACE LATTICE VECTORS
  REAL(8)                :: NKDIV(3)
  REAL(8)       ,POINTER :: XK(:,:)   !(3,NKP) IRR. K-POINTS IN RELATIVE COORDINATES
  INTEGER(4)    ,POINTER :: IKP(:,:) !(4,NTET) TETRAHEDRON CORNERS
  INTEGER(4)    ,POINTER :: MULT(:)  !(NTET) MULTIPLICITY OF THE TETRAHEDRON
END TYPE THIS_TYPE
TYPE(THIS_TYPE) :: THIS
END MODULE BRILLOUIN_MODULE
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!----  BLOCK TEST                                                   ----
!----  TEST SUBROUTINE                                              ----
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE BRILLOUIN$TESTING
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)                               
      INTEGER(4),PARAMETER :: NSYMX=50
      INTEGER(4),PARAMETER :: NBX=2
      REAL(8)              :: RBAS(3,3)
      INTEGER(4)           :: N(3)
      INTEGER(4)           :: IARB(3)
      INTEGER(4)           :: IIO(3,3,NSYMX)
      REAL(8)   ,ALLOCATABLE:: BK(:,:)
      REAL(8)   ,ALLOCATABLE:: EB(:,:)
      REAL(8)   ,ALLOCATABLE:: WGHT(:,:)
      REAL(8)   ,ALLOCATABLE:: A(:,:)
      INTEGER(4)           :: NFIL1=1001
!     ******************************************************************
                                    call trace$push('BRILLOUIN$TESTING')
!     ==================================================================
!     ==  READ LATTICE VECTORS                                        ==
!     ==================================================================
      OPEN(UNIT=NFIL1,FILE='TEST.IN')
      rewind nfil1
!     ==== RBAS(I,J) = REAL SPACE LATTICE VECTORS (I=X,Y,Z)
      DO I=1,3
        READ(NFIL1,*)RBAS(:,I)
      ENDDO
!     ==================================================================
!     ==  INPUT OF SYMMETRY OPERATIONS                                ==
!     ==================================================================
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
!     ==================================================================
!     ==  FIND IRREDUCIBLE K-POINTS AND TETRAHEDRA                    ==
!     ==================================================================
!      CALL BRILLOUIN$MSH(RBAS,NKP,NSYM,IIO,IARB)

      call BRILLOUIN$MSHNOSYM(.true.,RBAS,(/5,5,5/),(/1,1,1/))
 
!     ==================================================================
!     ==  CALCULATE ENERGIES AT THE IRREDUCIBLE K-POINTS              ==
!     ==================================================================
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
 
!     ==================================================================
!     ==  CALCULATE WEIGHTS                                           ==
!     ==================================================================
      CALL BRILLOUIN$DOS(NB,NKP,EB,WGHT,RNTOT,EF)
 
!     ==================================================================
!     ==  PERFORM BRILLOUIN ZONE INTEGRATION OF F(K)                  ==
!     ==================================================================
      SUMA=0.D0
      DO IKP=1,NKP
        DO IB=1,NB
          SUMA=SUMA+WGHT(IB,IKP)*A(IB,IKP)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  PRINT OUT                                                   ==
!     ==================================================================
      PRINT*,'INTEGRAL OF A : ',SUMA 
!
!     ==================================================================
!     ==  test derivatives                                            ==
!     ==================================================================
      CALL brillouin_testweightandder

                                    call trace$pop
      END SUBROUTINE BRILLOUIN$TESTING
!
!     ..................................................................
      subroutine brillouin_testweightandder
!     ** test routine for brillouin$weightander                       **
!     ** checks derivatives against numerical derivatives             **
      implicit none
      real(8) :: vol=1.d0
      real(8) :: e0(4)=(/2.d0,1.d0,0.d0,3.d0/)
      real(8) :: delta=1.d-8
      real(8) :: ef
      real(8) :: e(4)
      real(8) :: wght(4)
      real(8) :: wght0(4)
      real(8) :: dwght(4,4)
      real(8) :: dwght0(4,4)
      real(8) :: wghtarr(4,2,4)
      real(8) :: dwghttest(4,4)
      integer(4) :: I,J,k
!     ******************************************************************
      do k=1,3
        ef=0.3+real(k-1,kind=8)
        print*,'now test weightandder case ',k,ef
        call BRILLOUIN_WEIGHTandder(VOL,E0,EF,WGHT0,dwght0)
        do i=1,4
          do j=1,2
            e=e0
            e(i)=e(i)+real(2*j-3,kind=8)*delta
            call BRILLOUIN_WEIGHTandder(VOL,E,EF,WGHT,dwght)
            wghtarr(:,j,i)=wght
          enddo
          dwghttest(:,i)=(wghtarr(:,2,i)-wghtarr(:,1,i))/(2.D0*delta)
          do j=1,4
            print*,'test ',dwght(j,i),dwghttest(j,i),dwghttest(j,i)-dwght(j,i)
          enddo
        enddo
      enddo
      return
      end
!
!     ..................................................................
      subroutine brillouin_OPTIMIZE(NKP_,EKP_,WGHT_)
!      implicit none
!     ******************************************************************


      return
      end
!
!     ..................................................................
      SUBROUTINE BRILLOUIN$SETR8A(ID,LEN,VAL)
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(OUT):: VAL(LEN)
      INTEGER(4)              :: I
      REAL(8)                 :: GBAS(3,3)
      REAL(8)                 :: SVAR
!     ******************************************************************
      IF(ID.EQ.'RBAS')THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
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
!     ..................................................................
      SUBROUTINE BRILLOUIN$GETR8A(ID,LEN,VAL)
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(OUT):: VAL(LEN)
      INTEGER(4)              :: I
      REAL(8)                 :: GBAS(3,3)
      REAL(8)                 :: SVAR
!     ******************************************************************
      IF(ID.EQ.'K') THEN
        IF(3*THIS%NKP.NE.LEN) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$STOP('BRILLOUIN$GETR8A')
        END IF
        CALL GBASS(THIS%RBAS,GBAS,SVAR)
        DO I=1,THIS%NKP
          VAL(1+3*(I-1):3*I)=MATMUL(GBAS,THIS%XK(:,I))
        ENDDO
      ELSE IF(ID.EQ.'XK') THEN
        IF(3*THIS%NKP.NE.LEN) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$STOP('BRILLOUIN$GETR8A')
        END IF
        VAL(:)=RESHAPE(THIS%XK,(/LEN/))
      ELSE IF(ID.EQ.'WKPT') THEN
        IF(THIS%NKP.NE.LEN) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$STOP('BRILLOUIN$GETR8A')
        END IF
!       == CALCULATES WEIGHTS FOR A COMPLETELY FILLED BAND.            ==
!       == ATTENTION DEPENDS ON ACTUAL RBAS                            ==
        CALL BRILLOUIN_SAMFAC(1,THIS%NKP,(/(0.D0,I=1,THIS%NKP)/),1.D0,VAL,THIS)
      ELSE
        CALL ERROR$MSG('UNKNOWN ID')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BRILLOUIN$GETR8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE BRILLOUIN$GETI4(ID,VAL)
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     ******************************************************************
      IF(ID.EQ.'NK') THEN
        VAL=THIS%NKP
      ELSE
        CALL ERROR$MSG('UNKNOWN ID')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BRILLOUIN$GETI4')
      END IF
      RETURN
      END
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ----                                                         -----
!     ----  BLOCK ARBMSH:                                          -----
!     ----                                                         -----
!     ----  READS SYMMETRYELEMENTS AND FINDS IRR. K-POINTS         -----
!     ----  AND TETRAHEDRA, USED IN THE BRILLOUIN ZONE INTEGRATION -----
!     ----                                                         -----
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     .....................................................ARBMSH.......
      SUBROUTINE BRILLOUIN$MSH(RBAS,NGKP,NSYM,IIO,IARB)
!     **                                                              **
!     **  CALCULATE IRREDUCIBLE K-POINTS                              **
!     **  AND FINDS INEQUIVALENT TETRAHEDRA                           **
!     **  FOR BRILLOUIN ZONE INTEGRATION                              **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    1) IARB: DEPENDENCIES FOR DIVISIONS OF RECIPROCAL         **
!     **       LATTICE VECTORS.                                       **
!     **                ( IF IARB(1)=1 THEN 1ST AND 2ND LATTICE       **
!     **                  VECTORS ARE DIVIDED BY AN EQUAL NUMBER;     **
!     **                  IF IARB(2)=1 THEN SAME FOR 2ND AND 3RD;     **
!     **                  IF IARB(3)=1 THEN SAME FOR 3RD AND 1ST)     **
!     **                                                              **
!     **   AUTHOR: PETER E. BLOECHL                                   **
!     **                                                              **
!     **   SUBROUTINES USED:                                          **
!     **     GBASS,BASDIV,REDUZ,ZUORD,TETDIV,TETCNT,ORD1              **
!     **                                                              **
!     ******************************************************************
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)  ! LATTICE VECTORS (REAL SPACE)
      INTEGER(4),INTENT(IN) :: NGKP       ! TARGET FOR #(GENERAL K-POINTS)
      INTEGER(4),INTENT(IN) :: NSYM       ! #(SYMMETRY OPERATIONS)
      INTEGER(4),INTENT(IN) :: IIO(3,3,NSYM) !SYMMETRY OPERATIONS
      INTEGER(4),INTENT(IN) :: IARB(3)    ! DEPENDENCE
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
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
      INTEGER(4)            :: MWRIT
      INTEGER(4)            :: INDEX1,INDEX2,NDIMM,NDI1,NDI2
      INTEGER(4)            :: TET0(3,4,6)
      INTEGER(4)            :: N3
      REAL(8)               :: DUMMY
      INTEGER(4)            :: NTT !#(IRR. TETRAHEDRA)
      REAL(8)               :: V   ! UNIT CELL VOLUME
      LOGICAL(4)            :: TCHK
!     ******************************************************************
                                   call trace$push('BRILLOUIN$MSH')
      THIS%RBAS=RBAS
      INV=0
!     ------------------------------------------------------------------
!     -- DEFINE MESH                                                  --
!     ------------------------------------------------------------------
      NMSHP=NGKP                                                         
      CALL GBASS(RBAS,GBAS,DUMMY)                                             
      CALL BRILLOUIN_BASDIV(N,NMSHP,GBAS,IARB)                                    
!     ------------------------------------------------------------------
!     -- ATTEMPT TO SHIFT GRID                                        --
!     ------------------------------------------------------------------
      ISHIF(:)=1
      CALL BRILLOUIN_CHECKSHIFT(ISHIF,NSYM,IIO,TCHK)
      IF(.NOT.TCHK) ISHIF=0
!     ------------------------------------------------------------------
!     -- FIND IRREDUCIBLE K-POINTS                                    --
!     ------------------------------------------------------------------
      ALLOCATE(NUM(NMSHP))
      CALL BRILLOUIN_REDUZ(N,NMSHP,ISHIF,NSYM,IIO,NUM)
      IF(ISHIF(1).NE.0.OR.ISHIF(2).NE.0.OR.ISHIF(3).NE.0)INV=0   
      ALLOCATE(XK(3,NGKP))
      CALL BRILLOUIN_ZUORD(NMSHP,NUM,N,ISHIF,GBAS,NGKP,NKP,XK)            
      THIS%NKP=NKP
      ALLOCATE(THIS%XK(3,NKP))
      THIS%XK(:,:)=XK(:,1:NKP)
      DEALLOCATE(XK)
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
          WRITE(6,FMT='(I4,"(",3F10.6,")",I4," (",3F10.6,")")') &
     &          (K,(XK(J,K),J=1,3),K=NDI1,NDI2)                    
        ENDDO
      END IF                                                            
!     ----------------------------------------------------------------- 
!     --  CHOOSE TETRAHEDRA                                          -- 
!     ----------------------------------------------------------------- 
      CALL BRILLOUIN_TETDIV(N,GBAS,TET0)                                      
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
      MWRIT=100                                                         
      CALL BRILLOUIN_TETCNT(NMSHP,NUM,TET0,N,INV,THIS)
      DEALLOCATE(NUM)
                                   call trace$pop
      RETURN                                                            
      END                                                               
!     .....................................................ARBMSH.......
      SUBROUTINE BRILLOUIN$MSHNOSYM(TINV,RBAS,NKDIV,ISHIFT)
!     **                                                              **
!     **  CALCULATE IRREDUCIBLE K-POINTS                              **
!     **  AND FINDS INEQUIVALENT TETRAHEDRA                           **
!     **  FOR BRILLOUIN ZONE INTEGRATION                              **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    1) IARB: DEPENDENCIES FOR DIVISIONS OF RECIPROCAL         **
!     **       LATTICE VECTORS.                                       **
!     **                ( IF IARB(1)=1 THEN 1ST AND 2ND LATTICE       **
!     **                  VECTORS ARE DIVIDED BY AN EQUAL NUMBER;     **
!     **                  IF IARB(2)=1 THEN SAME FOR 2ND AND 3RD;     **
!     **                  IF IARB(3)=1 THEN SAME FOR 3RD AND 1ST)     **
!     **                                                              **
!     **   AUTHOR: PETER E. BLOECHL                                   **
!     **                                                              **
!     **   SUBROUTINES USED:                                          **
!     **     GBASS,BASDIV,REDUZ,ZUORD,TETDIV,TETCNT,ORD1              **
!     **                                                              **
!     ******************************************************************
      USE BRILLOUIN_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TINV       ! FLAG FOR TIME INVERSION SYMMETRY
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)  ! LATTICE VECTORS (REAL SPACE)
      INTEGER(4),INTENT(IN) :: NKDIV(3)   ! DIVISIONS OF RECIPROCAL LATTICE VECTORS
      INTEGER(4),INTENT(IN) :: ISHIFT(3)  ! GRID IS SHIFTED 
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
      REAL(8)   ,ALLOCATABLE:: XK(:,:)   !(3,NGKP) K-POINTS
      REAL(8)               :: GBAS(3,3)
      INTEGER(4)            :: N(3)       !DIVISION OF REC.LATT. VECT.
      INTEGER(4)            :: IWBK=0
      INTEGER(4)            :: IWNUM=0
      INTEGER(4)            :: NGKP
      INTEGER(4)            :: NMSHP
      INTEGER(4),ALLOCATABLE:: NUM(:)   ! POINTS FROM GENERAL TO IRR. K-POINT
      INTEGER(4)            :: ISHIF(3)
      INTEGER(4)            :: ISYM,I,J,K
      INTEGER(4)            :: INV
      INTEGER(4)            :: NKP
      INTEGER(4)            :: MWRIT
      INTEGER(4)            :: INDEX1,INDEX2,NDIMM,NDI1,NDI2
      INTEGER(4)            :: TET0(3,4,6)
      INTEGER(4)            :: N3
      REAL(8)               :: DUMMY
      INTEGER(4)            :: NTT !#(IRR. TETRAHEDRA)
      REAL(8)               :: V   ! UNIT CELL VOLUME
      INTEGER(4)            :: IIO(3,3,2)
      INTEGER(4)            :: NSYM  !#(SYMMETRY OPERATIONS) 1 OR 2
!     ******************************************************************
                                   call trace$push('BRILLOUIN$MSHNOSYM')
      this%rbas=rbas
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
!     ------------------------------------------------------------------
!     -- DEFINE MESH                                                  --
!     ------------------------------------------------------------------
      NMSHP=(NKDIV(1)+1)*(NKDIV(2)+1)*(NKDIV(3)+1)
      CALL GBASS(RBAS,GBAS,DUMMY)                                             
!     ------------------------------------------------------------------
!     -- FIND IRREDUCIBLE K-POINTS                                    --
!     ------------------------------------------------------------------
      ALLOCATE(NUM(NMSHP))
      CALL BRILLOUIN_REDUZ(NKDIV,NMSHP,ISHIF,NSYM,IIO,NUM)
!     == INV=1 EXPLOITS TIME INVERSION SYMMETRY FOR UNSHIFTED GRIDS  ==
!     == WHEN SEARCHING IRREDUCIBLE TETRAHEDRA. PROBABLY IRRELEVANT  ==
      IF(ISHIF(1).NE.0.OR.ISHIF(2).NE.0.OR.ISHIF(3).NE.0)INV=0   
      IF(.NOT.TINV) INV=0
!
!     ==================================================================
!     ==  DETERMINE IRREDUCIBLE K-POINTS                              ==
!     ==================================================================
      ALLOCATE(XK(3,NGKP))
      CALL BRILLOUIN_ZUORD(NMSHP,NUM,NKDIV,ISHIF,GBAS,NGKP,NKP,XK)            
      THIS%NKP=NKP
      ALLOCATE(THIS%XK(3,NKP))
      THIS%XK(:,:)=XK(:,1:NKP)
      DEALLOCATE(XK)
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
          WRITE(6,FMT='(I4,"(",3F10.6,")",I4," (",3F10.6,")")') &
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
      MWRIT=100                                                         
      CALL BRILLOUIN_TETCNT(NMSHP,NUM,TET0,NKDIV,INV,THIS)
      DEALLOCATE(NUM)
                                   call trace$pop
      RETURN                                                            
      END                                                               
!
!     ..................................................................
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
      INTEGER(4),PARAMETER     :: OPAR=1.D-6
      INTEGER(4)               :: I
      REAL(8)                  :: SVAR
!     ******************************************************************
!     ==================================================================
!     == FIND UPPER LIMIT FOR THE LENGTH OF SUBLATTICE VECTORS        ==
!     ==================================================================
      DO I=1,3                                                       
        BETR(I)=DSQRT(GBAS(I,1)**2+GBAS(I,2)**2+GBAS(I,3)**2)             
      ENDDO
      SVAR=(DBLE(NMSHP)/(BETR(1)*BETR(2)*BETR(3)))**(1.D0/3.D0)         
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
        N(1)=INT(DSQRT(RN(1)*RN(2))+OPAR)                               
        N(2)=N(1)                                                       
        N(3)=INT(RN(3))                                                 
      ELSE IF(IARB(2).EQ.1) THEN                                        
        N(1)=INT(DSQRT(RN(1)*RN(3))+OPAR)                               
        N(2)=INT(RN(2))                                                 
        N(3)=N(1)                                                       
      ELSE IF (IARB(3).EQ.1) THEN                                       
        N(1)=INT(RN(1))                                                 
        N(2)=INT(DSQRT(RN(2)*RN(3))+OPAR)                               
        N(3)=N(2)                                                       
      ELSE                                                              
        N(1)=INT(RN(1)+OPAR)                                            
        N(2)=INT(RN(2)+OPAR)                                            
        N(3)=INT(RN(3)+OPAR)                                            
      END IF                                                            
      N(1)=MAX0(1,N(1))                                                 
      N(2)=MAX0(1,N(2))                                                 
      N(3)=MAX0(1,N(3))                                                 
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
!     ..................................................................
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
      INTEGER(4)            :: IN(3,NMSHP)      
      INTEGER(4)            :: I,J,K,L,ISYM
      INTEGER(4)            :: I1,I2,I3,J1,J2,J3
      INTEGER(4)            :: IX1,IX2,IX3
      LOGICAL(4)            :: TCHK
!     ******************************************************************
      CALL BRILLOUIN_CHECKSHIFT(ISHIFT,NSYM,IO,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('SHIFT INCOMPATIBLE WITH SYMMETRY')
        CALL ERROR$STOP('BRILLOUIN_REDUZ')
      END IF
!     ==================================================================
!     ==  INITIALIZE                                                  ==
!     ==================================================================
      DO I=1,NMSHP                                                   
        NUM(I)=I                                                          
        IN(1,I)=(I-1)/((N(3)+1)*(N(2)+1))                                 
        IN(2,I)=(I-IN(1,I)*(N(2)+1)*(N(3)+1)-1)/(N(3)+1)                  
        IN(3,I)=I-IN(1,I)*(N(2)+1)*(N(3)+1)-IN(2,I)*(N(3)+1)-1            
      ENDDO
!     ==================================================================
!     ==  REDUCTION OF REDUCIBLE K-POINTS                             ==
!     ==================================================================
      DO I=1,NSYM                                                    
        DO J=1,NMSHP                                                   
          I1=2*IN(1,J)+ISHIFT(1)                                            
          I2=2*IN(2,J)+ISHIFT(2)                                            
          I3=2*IN(3,J)+ISHIFT(3)                                            
          J1=MOD(IO(1,1,I)*I1+IO(1,2,I)*I2+IO(1,3,I)*I3,2*N(1))             
          J2=MOD(IO(2,1,I)*I1+IO(2,2,I)*I2+IO(2,3,I)*I3,2*N(2))             
          J3=MOD(IO(3,1,I)*I1+IO(3,2,I)*I2+IO(3,3,I)*I3,2*N(3))             
          J1=J1+(1-ISIGN(1,J1))*N(1)                                        
          J2=J2+(1-ISIGN(1,J2))*N(2)                                        
          J3=J3+(1-ISIGN(1,J3))*N(3)                                        
          IX1=MOD(J1-ISHIFT(1),2)                                           
          IX2=MOD(J2-ISHIFT(2),2)                                           
          IX3=MOD(J3-ISHIFT(3),2)                                           
          IF(IX1.NE.0.OR.IX2.NE.0.OR.IX3.NE.0) CYCLE
          J1=(J1-ISHIFT(1))/2                                               
          J2=(J2-ISHIFT(2))/2                                               
          J3=(J3-ISHIFT(3))/2                                               
          NUM(J)=MIN0(NUM(J),J1*(N(2)+1)*(N(3)+1)+J2*(N(3)+1)+J3+1)                    
        ENDDO
      ENDDO
      RETURN                                                            
      END                                                               
!
!     ..................................................................
      SUBROUTINE BRILLOUIN_CHECKSHIFT(ISHIFT,NSYM,IO,TCHK)
!     **                                                              **
!     **  CHECKSHIFT CHECKS THE CONSISTENCY OF A K-GRID SHIFT WITH    **
!     **  THE SYMMETRY OPERATIONS                                      **
!     **  INPUT :                                                     **
!     **    N           NUMBER OF DIVISIONS OF REC. LATTICE VECTORS   **
!     **    NSYM        NUMBER OF SYMMETRY OPERATIONS                 **
!     **    IO          SYMMETRY OPERATIONS (BASIS ARE REC. LATT. VEC.)*
!     **  OUTPUT :                                                    **
!     **    TCHK        FLAG FOR CONSISTENCY OF THE SHIFT             **
!     **  REMARKS :                                                   **
!     **    THE MAPPING FROM COORDINATES TO NUMBERS IS GIVEN BY :     **
!     **   (X,Y,Z)=RBAS*(I,J,K)                                       **
!     **   (I,J,K)  <->  NUM = I*(N(2)+1)*(N(3)+1)+J*(N(3)+1)+K+1     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISHIFT(3)
      INTEGER(4),INTENT(IN) :: NSYM
      INTEGER(4),INTENT(IN) :: IO(3,3,NSYM)
      LOGICAL(4),INTENT(OUT):: TCHK
      INTEGER(4)            :: IN(3,8)      
      INTEGER(4)            :: I,J,K,L,ISYM
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
          IF(DMOD(DBLE(J1-ISHIFT(1)),2.D0).NE.0.D0.OR. &
     &       DMOD(DBLE(J2-ISHIFT(2)),2.D0).NE.0.D0.OR. &
     &       DMOD(DBLE(J3-ISHIFT(3)),2.D0).NE.0.D0) THEN                    
             TCHK=.FALSE.
             PRINT*,'SUBMESH SHIFT CONFLICTS WITH POINTGROUP'                
             WRITE(6,FMT="('SHIFT=',3I5)")ISHIFT
             WRITE(6,FMT="('SYMMETRYMATRIX NR. : ',I5/3(1H ,3I10/))") &
     &               I,((IO(K,L,I),L=1,3),K=1,3)                  
             RETURN 
          END IF                                                            
        ENDDO
      ENDDO
      TCHK=.TRUE.
      RETURN                                                            
      END                                                               
!
!     ..................................................................
      SUBROUTINE BRILLOUIN_ZUORD(NMSHP,NUM,N,ISHIFT,GBAS,IDKP,NKP,XK)             
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
      REAL(8)   ,INTENT(IN) :: GBAS(3,3)
      INTEGER(4),INTENT(OUT):: NKP
      REAL(8)   ,INTENT(OUT):: XK(3,IDKP)
      INTEGER(4)            :: NDIM
      INTEGER(4)            :: I1,I2,I3,I
!     ******************************************************************
      NDIM=0                                                            
      I=0                                                               
      DO I1=1,N(1)+1                                                 
        DO I2=1,N(2)+1                                                 
          DO I3=1,N(3)+1                                                 
            I=I+1                                                             
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
              NUM(I)=NDIM                                                     
              XK(1,NDIM)=(DBLE(I1-1)+DBLE(ISHIFT(1))/2.D0)/DBLE(N(1))              
              XK(2,NDIM)=(DBLE(I2-1)+DBLE(ISHIFT(2))/2.D0)/DBLE(N(2))              
              XK(3,NDIM)=(DBLE(I3-1)+DBLE(ISHIFT(3))/2.D0)/DBLE(N(3))              
            ELSE                                                              
              IF(NUM(I).GT.NMSHP) PRINT*,'ERROR'                              
              NUM(I)=NUM(NUM(I))                                              
            END IF                                                            
          ENDDO
        ENDDO
      ENDDO
      NKP=NDIM                                                          
      RETURN                                                            
      END                                                               
!
!     .....................................................TETDIV...... 
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
!     ------------------------------------------------------------------
!     -- SEARCH FOR THE SHORTEST DIAGONAL                             --
!     ------------------------------------------------------------------
      DO I=0,1                                                       
        DO J=0,1                                                       
          DO K=0,1                                                       
            ISVAR=4*I+2*J+K+1                                                 
            DO L=1,3                                                       
              P(ISVAR,L)=GBAS(1,L)*DBLE(I)/DBLE(N(1)) &
     &                  +GBAS(2,L)*DBLE(J)/DBLE(N(2)) &
     &                  +GBAS(3,L)*DBLE(K)/DBLE(N(3))
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
      RETURN                                                            
      END                                                               
!
!     .....................................................TETCNT.......
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
      INTEGER(4),PARAMETER  :: MWRIT=100
      LOGICAL(4),PARAMETER  :: TCHK=.TRUE.
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
      INTEGER(4),ALLOCATABLE:: ITET(:,:)
      INTEGER(4),ALLOCATABLE:: IY(:)
      INTEGER(4)            :: NKP    !#(IRR.K-POINTS)
      INTEGER(4)            :: NTT    !#(IRREDUCIBLE TETRAHEDRA)
      REAL(8)               :: V      ! VOLUME OF ONE TETRAHEDRON
      INTEGER(4)            :: IBGGST
      INTEGER(4)            :: NTMAX
      INTEGER(4)            :: IPP
      INTEGER(4)            :: NTET
      REAL(8)               :: SUM
      INTEGER(4)            :: IND,IP,NI
      INTEGER(4)            :: ISVAR1,ISVAR2
      INTEGER(4)            :: IXX,ISVAR
      INTEGER(4)            :: I,J,K,L,M,K1,K2,K3
!     ******************************************************************
      IBGGST=2*(2**30-1)+1
      NKP=THIS%NKP
      IF((NKP+1)**4.GT.HUGE(IY)) THEN
        CALL ERROR$MSG('TEST FOR BIGGEST INTEGER FAILED')
        CALL ERROR$MSG('CHECK WHETHER IBGGST IS THE BIGGEST INTEGER')
        CALL ERROR$I4VAL('IBGGST=',IBGGST)
        CALL ERROR$MSG('IF NOT, CORRECT IT IN ROUTINE TETCNT')
        CALL ERROR$MSG('ELSE REDUCE NUMBER OF K-POINTS')
        CALL ERROR$MSG('OR CHANGE TETCNT')
        CALL ERROR$STOP('BRILLOUIN_TETCNT')
      END IF
      ALLOCATE(ITET(4,6*N(3)))
      ALLOCATE(IY(6*NMSHP))
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
                ITET(K,L)=MIN0(ISVAR1,ISVAR2)                                     
                ITET(J,L)=MAX0(ISVAR1,ISVAR2)                                     
              ENDDO
            ENDDO
          ENDDO
                                                                        
!     --  IDENTIFY THE TETRAHEDRA WITH INTEGERS                         
          DO M=1,N(3)*6                                                  
            IPP=IPP+1                                                         
            IY(IPP)=ITET(1,M)*(NKP+1)**3 &
     &             +ITET(2,M)*(NKP+1)**2 &
     &             +ITET(3,M)*(NKP+1)    &
     &             +ITET(4,M)                                                 
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
!     ==  CHECK ORDERING AND CALCULATE NUMBER OF INEQUIVALENT TETRAHEDRA
      NTT=1
      DO I=1,NTET-1                                                 
        IF(IY(I+1).GT.IY(I)) THEN                                         
          NTT=NTT+1                                                       
        END IF                                                            
        IF(IY(I+1).LT.IY(I)) THEN                                         
          PRINT*,' I ',I,' IY ',IY(I)                                     
        END IF                                                            
      ENDDO
      WRITE(*,FMT='("NUMBER OF DIFFERENT TETRAHEDRA :",I5)')NTT                 
      THIS%NTET=NTT
      THIS%VOL =1.D0/DBLE(6*N(1)*N(2)*N(3))
!     ------------------------------------------------------------------
!     --  WRITE ON FILE                                               --
!     ------------------------------------------------------------------
      ALLOCATE(THIS%MULT(NTT))
      ALLOCATE(THIS%IKP(4,NTT))
      NTT=0
      DO I=1,NTET                                                   
        IF(I.NE.1.AND.IY(I).EQ.IY(I-1)) THEN                              
          THIS%MULT(NTT)=THIS%MULT(NTT)+1+INV
        ELSE                                                              
          NTT=NTT+1                                                       
          THIS%MULT(NTT)=1+INV
          ISVAR=IY(I)
          THIS%IKP(1,NTT)=ISVAR/(NKP+1)**3
          ISVAR=ISVAR-THIS%IKP(1,NTT)*(NKP+1)**3
          THIS%IKP(2,NTT)=ISVAR/(NKP+1)**2
          ISVAR=ISVAR-THIS%IKP(2,NTT)*(NKP+1)**2
          THIS%IKP(3,NTT)=ISVAR/(NKP+1)
          ISVAR=ISVAR-THIS%IKP(3,NTT)*(NKP+1)
          THIS%IKP(4,NTT)=ISVAR
        END IF                                                            
      ENDDO
      DEALLOCATE(IY)
                                                                        
!     ------------------------------------------------------------------
!     --  CHECK SUMRULE                                               --
!     ------------------------------------------------------------------
      sum=0.d0
      DO I=1,THIS%NTET
        SUM=SUM+THIS%VOL*DBLE(THIS%MULT(I))
      ENDDO
      IF(DABS(SUM-1.D0).GT.1.D-5) THEN                                  
        CALL ERROR$MSG('SUMRULE NOT FULLFILLED')
        CALL ERROR$MSG('SUM SHOULD BE EQUAL TO 1')
        CALL ERROR$R8VAL('SUM',SUM)
        CALL ERROR$I4VAL('THIS%NTET',THIS%NTET)
        CALL ERROR$R8VAL('THIS%VOL',THIS%VOL)
        CALL ERROR$I4VAL('N(1)',N(1))
        CALL ERROR$I4VAL('N(2)',N(2))
        CALL ERROR$I4VAL('N(3)',N(3))
        CALL ERROR$STOP('BRILLOUIN_TETCNT')
      END IF                                                            
      RETURN
      END
!
!     .....................................................ORD1.........
      SUBROUTINE BRILLOUIN_ORD1(NMAX,IX)
!     ******************************************************************
!     **                                                              **
!     **  ORD1 ORDERES THE ARRAY IX WITH SIZE NMAX ACCORDING TO       **
!     **  INCREASING NUMBER                                           **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NMAX
      INTEGER(4),INTENT(INOUT):: IX(NMAX)
      INTEGER(4)              :: WORK(NMAX)
      INTEGER(4)              :: IWORK(257)
      INTEGER(4)              :: I,J,ISVAR1,ISVAR2
!     ******************************************************************
      DO I=1,NMAX-1                                               
        DO J=I+1,NMAX                                               
          ISVAR1=IX(I)                                                    
          ISVAR2=IX(J)                                                    
          IX(I)=MIN0(ISVAR1,ISVAR2)                                       
          IX(J)=MAX0(ISVAR1,ISVAR2)                                       
        ENDDO
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
!     .....................................................DOS..........
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
      INTEGER(4),PARAMETER  :: MWRIT=100
      INTEGER(4)            :: I
      REAL(8)               :: SUM
      INTEGER(4)            :: ICHECK
!     ******************************************************************
!     ------------------------------------------------------------------
!     --  CALCULATE FERMI LEVEL                                       --
!     ------------------------------------------------------------------
      CALL BRILLOUIN_EFERMI(RNTOT,EF,TOLMAX,NB,NKP,EB,THIS)
                                                                        
      PRINT*,'FERMI ENERGY AT ',EF                                      
!     ------------------------------------------------------------------
!     --  CALCULATE WEIGHTS                                           --
!     ------------------------------------------------------------------
      CALL BRILLOUIN_SAMFAC(NB,NKP,EB,EF,WGHT,THIS)
                                                                        
      IF(ICHECK.EQ.0) RETURN                                            
!     ------------------------------------------------------------------
!     --  CHECK WHETHER SUMRULE IS FULLFILLED                         --
!     ------------------------------------------------------------------
      SUM=0.D0                                                          
      DO I=1,NB*NKP                                                 
        SUM=SUM+WGHT(I)                                                   
      ENDDO
      IF(DABS(SUM-RNTOT).GT.TOLMAX) THEN                                
        CALL ERROR$MSG('INTEGRATION FAILED')
        CALL ERROR$R8VAL('RESULT OF INTEGRATION: ',SUM)
        CALL ERROR$R8VAL('SHOULD BE: ',RNTOT)
        CALL ERROR$STOP('BRILLOUIN_DOS')
      END IF                                                            
      RETURN                                                            
      END                                                               
!
!     .....................................................EFERMI.......
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
      INTEGER(4)            :: IK,IB,I,ILOOP,IP
      REAL(8)               :: DE,EMIN,EMAX
      REAL(8)   ,ALLOCATABLE:: NOS(:)
      REAL(8)   ,ALLOCATABLE:: SOS(:)
      INTEGER(4)            :: INIT
      INTEGER(4)            :: IWGHT
      INTEGER(4)            :: ITET
      INTEGER(4)            :: MWRIT
      INTEGER(4)            :: NTET
      REAL(8)               :: SUM
      REAL(8)               :: TOL
      REAL(8)               :: VOL1
      REAL(8)               :: VOL
      REAL(8)               :: ESTEP
!     ******************************************************************
!     ------------------------------------------------------------------
!     --  FIND EMIN EMAX     (ENERGYBANDS ARE ASSUMED                 --
!     --                      TO BE ORDERED WITH RESPECT TO SYMMETRY  --
!     ------------------------------------------------------------------
      EMIN=EB(1,1)                                                      
      EMAX=EB(NB,1)                                                     
      DO IK=1,NKP                                                   
        DO IB=1,NB
          EMIN=DMIN1(EMIN,EB(IB,IK))                                        
          EMAX=DMAX1(EMAX,EB(IB,IK))                                        
        ENDDO
      ENDDO

!++++  CHRIS  FIX EMIN AND EMAX  AS MINIMUM OF THE ENERGY GRID  ++++
      DE=(EMAX-EMIN)/DBLE(NP-1)                                         
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
      SUM=0.D0                                                          
      DO ITET=1,NTET                                                
        VOL=THIS%VOL*DBLE(THIS%MULT(ITET))                                              
        IKP(:)=THIS%IKP(:,ITET)
        SUM=SUM+VOL                                                       
        DO IB=1,NB                                                    
          DO I=1,4                                                      
            E(I)=EB(IB,IKP(I))                                                
          ENDDO
          CALL BRILLOUIN_TOTNOS(VOL,E,EMIN,EMAX,NP,NOS,SOS)                   
        ENDDO
      ENDDO
      IF(DABS(SUM-1.D0).GT.1.D-5) THEN                                  
        CALL ERROR$MSG('TETRAHEDRA DO NOT FILL VOLUME')
        CALL ERROR$MSG('SUM SHOULD BE 1.')
        CALL ERROR$R8VAL('SUM ',SUM)
        CALL ERROR$STOP('BRILLOUIN_EFERMI')
      END IF                                                            
                                                                        
!     ------------------------------------------------------------------
!     --  GET FERMI LEVEL                                             --
!     ------------------------------------------------------------------
      TOL=TOLMAX
      CALL BRILLOUIN_EFI(RNTOT,TOL,EF,EMIN,EMAX,NP,NOS,SOS)               
                                                                        
!     ------------------------------------------------------------------
!     --  CHECK ACCURACY AND RESTART IF NECCESARY                     --
!     ------------------------------------------------------------------
      IF(TOL.GT.TOLMAX) THEN                                            
        ESTEP=(EMAX-EMIN)/DBLE(NP-1)                                    
        IP=1+(EF-EMIN)/ESTEP                                            
        EMIN=EMIN+ESTEP*DBLE(IP-1)                                      
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
      RETURN                                                            
      END                                                               
!     .....................................................EFI..........
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
      DO IFIND=1,NP                                                 
        IP=IPLOW+0.5*(IPUP-IPLOW)                                         
        NOSIP=NOS(IP)                                                     
        IF(RNTOT-NOSIP.GT.0) THEN                                         
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
      ESTEP=(EMAX-EMIN)/DBLE(NP-1)                                      
      ELOW=EMIN+DBLE(IPLOW-1)*ESTEP                                     
      DNOS=NOSUP-NOSLOW                                                 
      IF(DNOS.NE.0.D0) THEN                                             
        EFERMI=ELOW+(RNTOT-NOSLOW)/(NOSUP-NOSLOW)*ESTEP                 
      ELSE                                                              
        EFERMI=ELOW                                                     
      END IF                                                            
      IF(EFERMI-ELOW.LT.0) PRINT*,'ERROR IN EFI '                       
      RETURN                                                            
      END                                                               
!     .....................................................TOTNOS.......
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
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: VOL
      REAL(8)   ,INTENT(IN)  :: E_(4)
      REAL(8)   ,INTENT(IN)  :: EMIN
      REAL(8)   ,INTENT(IN)  :: EMAX
      INTEGER(4),INTENT(IN)  :: NP
      REAL(8)   ,INTENT(OUT) :: NOS(NP)
      REAL(8)   ,INTENT(OUT) :: SOS(NP)
      REAL(8)                :: E(4)
      REAL(8)                :: X
      REAL(8)                :: SVAR1,SVAR2
      REAL(8)                :: E21,E31,E41,E32,E42,E43
      REAL(8)                :: ESTEP,DE,EN
      REAL(8)                :: A,B,C,D
      INTEGER(4)             :: IMIN,IMAX
      INTEGER(4)             :: I,J
!     ******************************************************************
      E(:)=E_(:)
!
!     ------------------------------------------------------------------
!     --  INTEGRATION WITHOUT FERMISURFACE                            --
!     ------------------------------------------------------------------
      X=DMIN1(E(1),E(2),E(3),E(4))                                      
      IF(X.GE.EMAX) RETURN
      X=DMAX1(E(1),E(2),E(3),E(4))                                      
      IF(X.LE.EMIN) THEN                                                
        SOS(1)=SOS(1)+VOL                                               
        RETURN                                                          
      END IF                                                            
!     ------------------------------------------------------------------
!     --  ORDER ENERGIES                                              --
!     ------------------------------------------------------------------
      DO I=1,3                                                      
        DO J=I+1,4                                                    
          SVAR1=DMIN1(E(I),E(J))                                            
          SVAR2=DMAX1(E(I),E(J))                                            
          E(I)=SVAR1                                                        
          E(J)=SVAR2                                                        
        ENDDO
      ENDDO
!     ------------------------------------------------------------------
!     --  CALCULATE UNCORRECTED INTEGRAL AS MEANVALUE                 --
!     ------------------------------------------------------------------
      E21=E(2)-E(1)                                                     
      E31=E(3)-E(1)                                                     
      E41=E(4)-E(1)                                                     
      E32=E(3)-E(2)                                                     
      E42=E(4)-E(2)                                                     
      E43=E(4)-E(3)                                                     
      ESTEP=(EMAX-EMIN)/DBLE(NP-1)                                      
      IMIN=IDINT(2.D0+(E(1)-EMIN)/ESTEP)                                
      IMIN=MAX0(1,IMIN)                                                 
      IMAX=IDINT(1.D0+(E(2)-EMIN)/ESTEP)                                
      IMAX=MIN0(NP,IMAX)                                                
      EN=EMIN+ESTEP*(IMIN-1)                                            
      IF(IMAX.GE.IMIN) THEN                                             
        A=VOL/(E21*E31*E41)                                             
        DO I=IMIN,IMAX                                              
          NOS(I)=NOS(I)+A*(EN-E(1))**3                                    
          EN=EN+ESTEP                                                     
        ENDDO
      END IF                                                            
      IMIN=MAX0(1,IMAX+1)                                               
      IMAX=INT(1.D0+(E(3)-EMIN)/ESTEP)                                  
      IMAX=MIN0(NP,IMAX)                                                
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
      IMIN=MAX0(1,IMAX+1)                                               
      IMAX=INT(1.D0+(E(4)-EMIN)/ESTEP)                                  
      IMAX=MIN0(NP,IMAX)                                                
      IF(E43.GT.0.D0) THEN                                              
        A=VOL                                                           
        D=VOL/(E41*E42*E43)                                             
        DO I=IMIN,IMAX                                              
          NOS(I)=NOS(I)+A+D*(EN-E(4))**3                                  
          EN=EN+ESTEP                                                     
        ENDDO
      END IF                                                            
      IMIN=MAX0(1,IMAX+1)                                               
      IF(IMIN.GT.NP) RETURN                                             
      SOS(IMIN)=SOS(IMIN)+VOL                                           
      RETURN                                                            
      END                                                               
!
!     .....................................................SAMFAC.......
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
      TYPE(THIS_TYPE),INTENT(IN)  :: THIS
      REAL(8)       ,INTENT(OUT) :: WGHT(NB,NKP)
      REAL(8)                    :: E(4)
      REAL(8)                    :: WGHT0(4)
      REAL(8)                    :: VOL0,VOL
      INTEGER(4)                 :: IKP(4)
      INTEGER(4)                 :: ITET,IWGHT,MWRIT,NTET,INIT
      INTEGER(4)    ,ALLOCATABLE :: WORK(:)
      INTEGER(4)                 :: IB,I
!     ******************************************************************
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
      RETURN                                                            
      END                                                               
!
!     .....................................................SAMFAC.......
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
      REAL(8)       ,INTENT(OUT) :: dWGHT(NB,NKP,nkp)
      REAL(8)                    :: E(4)
      REAL(8)                    :: WGHT0(4)
      REAL(8)                    :: dWGHT0(4,4)
      REAL(8)                    :: VOL0,VOL
      INTEGER(4)                 :: IKP(4)
      INTEGER(4)                 :: ITET,IWGHT,MWRIT,NTET,INIT
      INTEGER(4)    ,ALLOCATABLE :: WORK(:)
      INTEGER(4)                 :: IB,I,j
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
            do j=1,4
              DWGHT(IB,IKP(I),IKP(J))=DWGHT(IB,IKP(I),IKP(J))+DWGHT0(I,J)                          
            enddo
          ENDDO
        ENDDO
      ENDDO
      RETURN                                                            
      END                                                               
!
!     .....................................................WHEIGT.......
      SUBROUTINE BRILLOUIN_WEIGHT(VOL,E_,EF,WGHT)                                  
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE WEIGHTS FOR TETRAHEDRON-SAMPLING             **
!     **  CORRESPONDING TO INTEGRATION OVER ONE TETRAHEDRON           **
!     **                                                              **
!     **  CORRECTION FOR THE NONLINEAR SHAPE INCLUDED IF ICOR=1       **
!     **                                                              **
!     **  AUTHOR : P.BLOECHL                                          **
!     ******************************************************************
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
!     ******************************************************************
      E(:)=E_(:)
      WGHT(:)=0.D0
!     ------------------------------------------------------------------
!     --  INTEGRATION WITHOUT FERMISURFACE                            --
!     ------------------------------------------------------------------
      X=MIN(E(1),E(2),E(3),E(4))                                      
      IF(X.GE.EF) RETURN
      X=MAX(E(1),E(2),E(3),E(4))                                      
      IF(X.LE.EF) THEN                                                  
        VPRIME=.25D0*VOL                                                
        DO I=1,4                                                     
          WGHT(I)=VPRIME                                                  
        ENDDO
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
        CALL ERROR$MSG('ERROR')
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
!     .....................................................WHEIGT.......
      SUBROUTINE BRILLOUIN_WEIGHTandder(VOL,E_,EF,WGHT,dwght)                                  
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE WEIGHTS FOR TETRAHEDRON-SAMPLING             **
!     **  CORRESPONDING TO INTEGRATION OVER ONE TETRAHEDRON           **
!     **                                                              **
!     **  CORRECTION FOR THE NONLINEAR SHAPE INCLUDED IF ICOR=1       **
!     **                                                              **
!     **  AUTHOR : P.BLOECHL                                          **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)  ,INTENT(IN) :: VOL      ! WEIGHT OF THIS TETRAHEDRON
      REAL(8)  ,INTENT(IN) :: EF       ! FERMI LEVEL
      REAL(8)  ,INTENT(IN) :: E_(4)    ! ENERGY BANDS AT THE CORNERS
      REAL(8)  ,INTENT(OUT):: WGHT(4)  ! INTEGRATION WEIGHTS
      REAL(8)  ,INTENT(OUT):: dwght(4,4)! derivative ofINTEGRATION WEIGHTS
                                       ! with respect to energies
      INTEGER(4),PARAMETER :: ICOR=1   ! ON/OFF SWITCH FOR CORRECTION 
      REAL(8)              :: E(4)     ! ENERGY BANDS AT THE CORNERS
      REAL(8)              :: FA(4),dfa(4,4)
      REAL(8)              :: FB(4)
      INTEGER(4)           :: INDEX(4)
      REAL(8)              :: X
      REAL(8)              :: wghtfac(4),dwghtfac(4,4)
      REAL(8)              :: VPRIME,dvprime(4)
      REAL(8)              :: DE
      REAL(8)              :: DOS,ddos(4)
      REAL(8)              :: E21,E31,E41,E32,E42,E43
      REAL(8)              :: DE1,DE2,DE3,DE4
      INTEGER(4)           :: I,J,IP,N,M,K
      REAL(8)              :: DA,DB,DC
      REAL(8)              :: VOL14
!     ******************************************************************
      E(:)=E_(:)
      WGHT(:)=0.D0
!     ------------------------------------------------------------------
!     --  INTEGRATION WITHOUT FERMISURFACE                            --
!     ------------------------------------------------------------------
      X=MIN(E(1),E(2),E(3),E(4))                                      
      IF(X.GE.EF) RETURN
      X=MAX(E(1),E(2),E(3),E(4))                                      
      IF(X.LE.EF) THEN                                                  
        VPRIME=.25D0*VOL                                                
        DO I=1,4                                                     
          WGHT(I)=VPRIME                                                  
        ENDDO
        dwght(:,:)=0.d0
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
        DOS=12.d0*VPRIME/de
        DDOS(:)=12.d0*DVPRIME(:)/de
        DDOS(1)=ddos(1)+dos/de
!
!     ======================================================================
!     == second case e1,e2<ef<e3,e4                                       ==
!     ======================================================================
      ELSE IF(EF.GT.E(2).AND.EF.LT.E(3)) THEN                           
        DE1=EF-E(1)                                                     
        DE2=EF-E(2)                                                     
        DE3=E(3)-EF                                                     
        DE4=E(4)-EF                                                     
!       ------  TETRAHEDRON X1,X2,X13',X14'                             
        VPRIME=0.25D0*VOL*DE1**2/(E41*E31)
        dvprime(1)=-0.5d0*vol*de1/(E41*E31)+vprime*(1.d0/e41+1.d0/e31)
        dvprime(2)=0.d0
        dvprime(3)=-vprime/e31
        dvprime(4)=-vprime/e41
!
        wghtfac(1)=3.D0-DE1/E41-DE1/E31
        dwghtfac(1,1)=(1.d0-de1/e41)/e41+(1.d0-de1/e31)/e31
        dwghtfac(1,2)=0.d0
        dwghtfac(1,3)=+de1/e31**2
        dwghtfac(1,4)=+de1/e41**2
!
        wghtfac(2)=1.d0
        dwghtfac(2,:)=0.d0
!
        wghtfac(3)=DE1/E31
        dwghtfac(3,1)=-(1.d0-de1/e31)/e31
        dwghtfac(3,2)=0.d0
        dwghtfac(3,3)=-de1/e31**2
        dwghtfac(3,4)=0.d0

        wghtfac(4)=DE1/E41
        dwghtfac(4,1)=-(1.d0-de1/e41)/e41
        dwghtfac(4,2)=0.d0
        dwghtfac(4,3)=0.d0
        dwghtfac(4,4)=-de1/e41**2
!
        wght(:)=vprime*wghtfac(:)
        DO I=1,4        
          DWGHT(:,I)=DVPRIME(I)*WGHTFAC(:)+VPRIME*DWGHTFAC(:,I)
        ENDDO
!       ------  TETRAHEDRON X2,X13',X23',X14'                           
        VPRIME=.25D0*VOL*DE2*DE3*DE1/(E32*E31*E41)                      
        dvprime(1)=vprime*(-1.d0/de1+1.d0/e31+1.d0/e41)
        dvprime(2)=vprime*(-1.d0/de2+1.d0/e32)
        dvprime(3)=vprime*(1.d0/de3-1.d0/e32-1.d0/e31)
        dvprime(4)=vprime*(-1.d0/e41)
!
        WGHTfac(1)=2.D0-DE1/E31-DE1/E41
        dwghtfac(1,1)=(1.d0-de1/e31)/e31+(1.d0-de1/e41)/e41
        dwghtfac(1,2)=0.d0
        dwghtfac(1,3)=de1/e31**2
        dwghtfac(1,4)=de1/e41**2
!
        WGHTfac(2)=2.D0-DE2/E32
        dwghtfac(2,1)=0.d0
        dwghtfac(2,2)=(1.d0-de2/e32)/e32
        dwghtfac(2,3)=de2/e32**2
        dwghtfac(2,4)=0.d0
!
        WGHTfac(3)=DE2/E32+DE1/E31
        dwghtfac(3,1)=-(1.d0-de1/e31)/e31
        dwghtfac(3,2)=-(1.d0-de2/e32)/e32
        dwghtfac(3,3)=-de2/e32**2-de1/e31**2
        dwghtfac(3,4)=0.d0
!
        WGHTfac(4)=DE1/E41
        dwghtfac(4,1)=-(1.d0-de1/e41)/e41
        dwghtfac(4,2)=0.d0
        dwghtfac(4,3)=0.d0
        dwghtfac(4,4)=-de1/e41**2
!
        wght(:)=wght(:)+vprime*wghtfac(:) 
        DO I=1,4        
          DWGHT(:,I)=dwght(:,i)+DVPRIME(I)*WGHTFAC(:)+VPRIME*DWGHTFAC(:,I)
        ENDDO
!
!       ------  TETRAHEDRON X2,X23',X24',X14'                           
        VPRIME=.25D0*VOL*DE2**2*DE4/(E42*E32*E41)                       
        dvprime(1)=vprime/e41
        dvprime(2)=vprime*(-2.d0/de2+1.d0/e42+1.d0/e32)
        dvprime(3)=vprime*(-1.d0/e32)
        dvprime(4)=vprime*(1.d0/de4-1.d0/e42-1.d0/e41)

        WGHTfac(1)=1.D0-DE1/E41
        dwghtfac(1,1)=(1.d0-de1/e41)/e41
        dwghtfac(1,2)=0.d0
        dwghtfac(1,3)=0.d0 
        dwghtfac(1,4)=de1/e41**2
!
        WGHTfac(2)=3.D0-DE2/E32-DE2/E42
        dwghtfac(2,1)=0.d0
        dwghtfac(2,2)=(1.d0-de2/e32)/e32+(1.d0-de2/e42)/e42
        dwghtfac(2,3)=de2/e32**2
        dwghtfac(2,4)=de2/e42**2

        WGHTfac(3)=DE2/E32
        dwghtfac(3,1)=0.d0
        dwghtfac(3,2)=-(1.d0-de2/e32)/e32
        dwghtfac(3,3)=-de2/e32**2
        dwghtfac(3,4)=0.d0
!
        WGHTfac(4)=DE2/E42+DE1/E41
        dwghtfac(4,1)=-(1.d0-de1/e41)/e41
        dwghtfac(4,2)=-(1.d0-de2/e42)/e42
        dwghtfac(4,3)=0.d0
        dwghtfac(4,4)=-de2/e42**2-de1/e41**2
!
        wght(:)=wght(:)+vprime*wghtfac(:) 
        DO I=1,4        
          DWGHT(:,I)=dwght(:,i)+DVPRIME(I)*WGHTFAC(:)+VPRIME*DWGHTFAC(:,I)
        ENDDO
!
!       ------  DOS=A+B*(EF-E2)+C*(EF-E2)**2                            
        DA=3.D0*VOL*E21/(E31*E41)                                       
        DB=6.D0*VOL/(E31*E41)                                           
        DC=-3.D0*VOL/(E32*E41*E31*E42)*(E31+E42)                        
        DOS=DA+DB*DE2+DC*DE2**2                                         
        ddos(1)=da*(-1.d0/e21+1.d0/e31+1.d0/e41) &
       &       +db*(1.d0/e31+1.d0/e41)*de2 &
       &       +dc*(1.d0/e41+1.d0/e31-1.d0/(e31+e42))*de2**2
        ddos(2)=da*(1.d0/e21) &
       &       +dc*(1.d0/e32+1.d0/e42-1.d0/(e31+e42))*DE2**2 &
       &       -db-2.d0*dc*de2
        ddos(3)=da*(-1.d0/e31) &
       &       +db*(-1.d0/e31)*de2 &
       &       +dc*(-1.d0/e32-1.d0/e31+1.d0/(e31+e42))*de2**2
        ddos(4)=da*(-1.d0/e41) &
       &       +db*(-1.d0/e41)*de2 &
       &       +dc*(-1.d0/e41-1.d0/e42+1.d0/(e31+e42))*de2**2 
!
!     ======================================================================
!     == last case                                                        ==
!     ======================================================================
      ELSE IF(EF.GE.E(3).AND.EF.LT.E(4)) THEN                           
        DE=E(4)-EF                                                      
        VPRIME=.25D0*VOL*DE**3/(E41*E42*E43)                            
        dvprime(1)=vprime*(1.d0/e41)
        dvprime(2)=vprime*(1.d0/e42)
        dvprime(3)=vprime*(1.d0/e43)
        dvprime(4)=vprime*3.d0/de-vprime*(1.d0/e41+1.d0/e42+1.d0/e43)
!
        wghtfac(1)=-DE/E41
        dwghtfac(1,1)=-de/e41**2
        dwghtfac(1,2)=0.d0
        dwghtfac(1,3)=0.d0
        dwghtfac(1,4)=-(1.d0-de/e41)/e41
!
        wghtfac(2)=-DE/E42
        dwghtfac(2,1)=0.d0
        dwghtfac(2,2)=-de/e42**2
        dwghtfac(2,3)=0.d0
        dwghtfac(2,4)=-(1.d0-de/e42)/e42

        wghtfac(3)=-DE/E43
        dwghtfac(3,1)=0.d0
        dwghtfac(3,2)=0.d0
        dwghtfac(3,3)=-de/e43**2
        dwghtfac(3,4)=-(1.d0-de/e43)/e43

        wghtfac(4)=-(4.D0-DE/E41-DE/E42-DE/E43)
        dwghtfac(4,1)=de/e41**2
        dwghtfac(4,2)=de/e42**2
        dwghtfac(4,3)=de/e43**2
        dwghtfac(4,4)=(1.d0-de/e41)/e41+(1.d0-de/e42)/e42+(1.d0-de/e43)/e43
        VOL14=.25D0*VOL                                                 
        wght(:)=vol14+vprime*wghtfac(:) 
        DO I=1,4        
          DWGHT(:,I)=DVPRIME(I)*WGHTFAC(:)+VPRIME*DWGHTFAC(:,I)
        ENDDO
!       ------  PARAMETERS FOR CORRECION                                
        DOS=12.D0*VPRIME/de                                
        ddos(1)=12.d0*dvprime(1)/de
        ddos(2)=12.d0*dvprime(2)/de
        ddos(3)=12.d0*dvprime(3)/de
        ddos(4)=12.d0*dvprime(4)/de-dos/de
      ELSE                                                              
        CALL ERROR$MSG('ERROR')
        CALL ERROR$STOP('BRILLOUIN_WEIGHT')
      END IF                                                            
!     ------------------------------------------------------------------
!     --  ADD CORRECTION FOR QUADRATIC DEVIATION                      --
!     ------------------------------------------------------------------
      IF(ICOR.EQ.1) THEN                                                
        DO M=1,4                                                    
          DO N=1,4                                                    
            WGHT(M)=WGHT(M)+.25D0*(E(N)-E(M))*DOS*.1D0                      
            dwght(m,:)=dwght(m,:)+0.25*(e(n)-e(m))*ddos(:)*0.1d0
            dwght(m,n)=dwght(m,n)+0.25*dos*0.1d0
            dwght(m,m)=dwght(m,m)-0.25*dos*0.1d0
          ENDDO
        ENDDO
      END IF                                                            
!     ------------------------------------------------------------------
!     --  REORDER WEIGHTS                                             --
!     ------------------------------------------------------------------
      DO I=1,4                                                      
        FA(INDEX(I))=WGHT(I)                                              
        FB(INDEX(I))=E(I)                                                 
        do j=1,4
          dfa(index(i),index(j))=dwght(i,j)
        enddo
      ENDDO
      DO I=1,4                                                      
        WGHT(I)=FA(I)                                                     
        E(I)=FB(I)                                                        
        do j=1,4
          dwght(i,j)=dfa(i,j)
        enddo
      ENDDO
      RETURN                                                            
      END                                                               
