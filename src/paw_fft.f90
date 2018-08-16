!======================================================================
!======================================================================
!====                                                              ====
!====  PLANEWAVE: GENERIC PLANE WAVE OPERATION                     ====
!====                                                              ====
!======================================================================
!======================================================================
!======================================================================
!**********************************************************************
!**********************************************************************
!**                                                                  **
!**  FUNCTIONS:                                                      **
!**    DIVIDERGRIDONTASKS(EPWRHO,RBAS,NR1START,NR1L,NR2,NR3)         **        
!**    INITIALIZE(KEY,RBAS,K,EPW,NR1START,NR1L,NR2,NR3)              **
!**    SELECT(KEY)                                                   **
!**    SUPFFT(ID,NFFT,NGL,FOFG,NRL,FOFR)                             **
!**    FFT(ID,NFFT,NGL,FOFG,NRL,FOFR)                                **
!**    STRUCTUREFACTOR(R,NGL,EIGR)                                   **
!**    INVERTG(NGL,FOFG,FOFMG)                                       **
!**    SCALARPRODUCT(ID,NGL,NDIM,N1,F1,N2,F2,MAT)                    **
!**    COLLECT(NGL,FL,NGG,FG)    : R8; I4, C8                        **
!**    DISTRIBUTE(NGL,FL,NGG,FG) : R8, I4, C8                        **
!**    GETL4(ID,VAL) : SUPER; TINV                                   **
!**    GETR8A(ID,LEN,VAL): GVEC; G2                                  **
!**    GETC8A(ID,LEN,VAL): EIKR                                      **
!**    GETI4(ID,VAL) : NGL,NRL,NGAMMA                                **
!**    SETL4(ID,VAL) : SUPER; TINV                                   **
!**                                                                  **
!**  DESCRIPTION:                                                    **
!**    WAVE FUNCTIONS PSI(R) ARE DECOMPOSED INTO A PERIODIC PHI(R)   **
!**    AND A PHASE FACTOR EIKR                                       **
!**    THE PROGRAM WORKS WITH PSI IN G-SPACE AND PHI IN R-SPACE      **
!**    (EXCEPT FOR SUPER WAVE FUNCTIONS)                             **
!**                                                                  **
!**    SUPERWAVE FUNCTIONS WITH SUPERWAVE FUNCTION INTERFACE:        **
!**    THE OPTION OF SUPERWAVE FUNCTIONS IS POSSIBLE FOR THE GAMMA   **
!**    POINT ONLY, I.E. FOR THE CHARGE DENSITY. THE WAVE FUNCTIONS   **
!**    IN R-SPACE ARE REAL (PSI, NOT PHI) AND ONLY ONE-HALF OF THE   **
!**    RECIPROCAL SPACE VECTORS NEED TO BE KEPT, NAMELY THOSE WITH   **
!**    MINUSG(IG).GE.IG                                              **
!**                                                                  **
!**    SUPERWAVE FUNCTIONS WITHOUT SUPERWAVE FUNCTION INTERFACE:     **
!**    SUPER WAVEFUNCTIONS CAN ALSO BE HANDLED EXTERNALLY. IN THIS   **
!**    CASE PSI=PSI1+I*PSI2, WHEREW PSI1,AND PSI2 ARE REAL WAVE      **
!**    FUNCTIONS. NOTE THAT                                          **
!**      PHI1=RE[PHI*EIKR]/EIKR AND                                  **
!**      PHI2=IM[PHI*EIKR]/EIKR, THAT IS THE WAVE FUNCTIONS ARE NOT  **
!**    AUTOMATICALLY UNRAVELED IN R-SPACE. THIS ALLOWS TO TREAT      **
!**    A PAIR OF REAL WAVE FUNCTIONS IN ONE COMPLEX ARRAY.           **
!**    THE FUNCTION "INVERTG" ALLOWS TO OBTAIN CONJG(PSI(-G))        **
!**    AS NEW WAVE FUNCTION, WHICH IS NEEDED TO UNRAVEL THE WAVE     **
!**    FUNCTIONS IN G-SPACE:                                         **
!**       PSI1(G)=0.5   *(PSI(G)+CONJG(PSI(-G))                      **
!**       PSI2(G)=0.5*CI*(PSI(G)-CONJG(PSI(-G))                      **
!**                                                                  **
!**    A "STRIPE" IS A COLUMN OF G-VECTORS ALONG THE FIRST REC.      **
!**      LATTICE VECTOR                                              **
!**    A "PLANE" IS A PLANE OF G-VECTORS THAT DIFFER ONLY BY         **
!**      THE SECOND AND THIRD LATTICE VECTORS                        **
!**    THE PARALLELIZATION IS SUCH THAT A NUMBER OF G-VECTORS        **
!**    ARE KEPT IN G-SPACE SO THAT NO STRIPE IS DIVIDED BETWEEN      **
!**    TWO TASKS. IN REAL SPACE A NUMBER OF PLANES IS KEPT ON EACH   **
!**    TASK.                                                         **
!**                                                                  **
!**    THE FOURIER TRANSFORM FROM G-TO R-SPACE, MAPS G-VECTORS       **
!**    ONTO AND ARRAY OF STRIPES USING THE ARRAY "IGTOSTRIPE".       **
!**    THE STRIPES ARE DIVIDED ACCORING TO SETS OF PLANES            **
!**    FOR EACH TASK. THE CORRESPONDING PORTIONS ARE EXCHANGED,      **
!**    SO THAT EACH TASK HOLDS THE INFORMATION OF ALL STRIPES,       **
!**    BUT ONLY THE PART CORRESPONDING TO ITS PLANES.                **
!**    THE PARTIAL STRIPES ARE DISTRIBUTED ONTO PLANES USING THE     **
!**    ARRAY "ISTRIPEOFYZ". WE CONSIDER                              **
!**    THAT SOME LINE-FFTS NEED NOT BE DONE, BECAUSE THERE ARE       **
!**    NO NON-ZERO VECTORS. FINALLY THE THIRD DIRECTION IS PROCESSED.**
!**                                                                  **
!**                                                                  **
!**********************************************************************
!**********************************************************************
!.......................................................................
MODULE PLANEWAVE_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: PLANEWAVE                                                  **
!**                                                                   **
!**  PURPOSE:                                                         **
!**                                                                   **
!**  REMARKS:                                                         **
!**  EVERY DATA OBJECT IS IDENTIFIED BY A KEYWORD AND A NUMBER        **
!**  EVERY DATA OBJECT HAS AN ENTRY IN A LINKEDLIST ($LISTSTART)      **
!**    WHICH STORES ALL NECCESARY INFORMATION                         **
!**                                                                   **
!***********************************************************************
TYPE PWPARALLEL_TYPE
  CHARACTER(32)      :: ID         ! PRODUCED BY "KEYCRYPT"
  CHARACTER(32)      :: CID        ! COMMUNICATOR ID FOR PARALLELIZATION
  REAL(8)            :: KVEC(3)    ! K-POINT IN RELATIVE COORDINATES
  REAL(8)            :: GBAS(3,3)  ! RECIPROCAL LATTICE VECTORS
  INTEGER(4)         :: NR1        ! 1ST DIMENSION FO FFT GRID (GLOBAL)
  INTEGER(4)         :: NR2        ! 2ND DIMENSION FO FFT GRID
  INTEGER(4)         :: NR3        ! 3RD DIMENSION FO FFT GRID
  INTEGER(4)         :: NR1L       ! #(PLANES ON THIS TASK)
  INTEGER(4)         :: NR1START   ! INDEX OF THE FIRST LOCAL PLANE
  INTEGER(4),POINTER :: NR1LARR(:) !(NTASKS) #(PLANES ON EACH TASK)
  INTEGER(4),POINTER :: NGLARR(:)  !(NTASKS) #(GVECTORS ON EACH TASK)
  INTEGER(4),POINTER :: NSTRIPELARR(:)  !(NTASKS) #(STRIPES ON EACH TASK)
  INTEGER(4),POINTER :: IGVEC(:,:) !(3,NGL) G=GBAS*(KVEC+IGVEC)
  REAL(8)            :: RWEIGHT    !R-INTEGRATION WEIGHT
  REAL(8)            :: GWEIGHT    !G-INTEGRATION WEIGHT
  !== ARRAYS USED FOR FFT ==============================================
  INTEGER(4),POINTER :: IGTOSTRIPE(:)    !(NGL)  IG -> IR1+NR1*(ISTRIPE-1)
  INTEGER(4),POINTER :: ISTRIPETOYZ(:,:) !(NSTRIPELX,NTASKS) & 
                                         !STRIPE -> (IR2+NR2*(IR3-1))
  !== ARRAYS USED FOR COLLECT AND DISTRIBUTE ===========================
  INTEGER(4),POINTER :: TASKOFIG(:)      !(NG) (IG -> TASK)
  !== USED FOR SUPERWAVE FUNCTIONS (TIME INVERSION SYMMETRY)
  LOGICAL(4)         :: TINV       ! (G <-> -G)?
  LOGICAL(4)         :: TSUPER     ! SUPER WAVE FUNCTION INTERFACE SELECTED
  INTEGER(4),POINTER :: NGHARR(:)  ! #(INDEPENDENT G-VECTORS PER TASK)
  INTEGER(4),POINTER :: MINUSG(:)  ! (NGL) IG(G) -> IG(-G)
  INTEGER(4)         :: NGAMMA     ! G(NGAMMA) IS THE GAMMA POINT 
                                   ! NGAMMA=0 IF NOT ON THIS TASK 
  INTEGER(4),POINTER :: IGVECH(:,:) ! IGVEC OFR MINUSG(IG)>IG
  ! POINTER TO NEXT THIS
  TYPE(PWPARALLEL_TYPE),POINTER :: NEXT
END TYPE PWPARALLEL_TYPE
LOGICAL(4)                     :: TINI=.FALSE.
TYPE (PWPARALLEL_TYPE),POINTER :: THIS
END MODULE PLANEWAVE_MODULE
!*******************************************************************************
!     
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PLANEWAVE$REPORT(NFIL)
!     **************************************************************************
!     **  REPORT SETTINGS OF THE PLANEWAVE INSTANCES                          **
!     **                                                                      **
!     **  NOTE, THAT INSTANCES MAY EXIST ON ONLY A PORTION OF ALL TASKS       **
!     **                                                                      **
!     **************************************************************************
      USE MPE_MODULE
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)  :: NFIL
      INTEGER(4)              :: ITASK
      INTEGER(4)              :: NTASKS,THISTASK
      INTEGER(4)              :: NTASKS1,THISTASK1
      TYPE (PWPARALLEL_TYPE),POINTER :: THIS1
      CHARACTER(64)           :: ID
      LOGICAL(4)              :: TINV
      LOGICAL(4)              :: TSUPER
      REAL(8)                 :: KVEC(3)
      INTEGER(4)              :: NR1,NR2,NR3
      INTEGER(4)              :: I,J
      INTEGER(4),ALLOCATABLE  :: ICOUNT(:)
      INTEGER(4),ALLOCATABLE  :: IWORK(:,:)
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     ==========================================================================
!     == DETERMINE THE MAIN TASKS FOR ALL INSTANCES.                          ==
!     == TASK ITASK IS RESPONSIBLE (IS MAIN TASK) FOR ICOUNT(ITASK) INSTANCES.==
!     == THE NUMBER OF INSTANCES IS SUM(ICOUNT).                              ==
!     ==========================================================================
      ALLOCATE(ICOUNT(NTASKS))
      ICOUNT(:)=0
      THIS1=>THIS
      DO 
        CALL MPE$QUERY(THIS1%CID,NTASKS1,THISTASK1)
        IF(THISTASK1.EQ.1) ICOUNT(THISTASK)=ICOUNT(THISTASK)+1
        THIS1=>THIS1%NEXT
        IF(ASSOCIATED(THIS,THIS1)) EXIT
      ENDDO
      CALL MPE$COMBINE('MONOMER','+',ICOUNT)
!
!     ==========================================================================
!     == NOW LOOP OVER ALL INSTANCES                                          ==
!     ==========================================================================
      DO ITASK=1,NTASKS
!       __ONLY MAIN TASKS OF INSTANCES REPORT INFORMATION_______________________
        IF(ICOUNT(ITASK).EQ.0) CYCLE  
!       __ONLY MAIN TASK OF AN INSTANCE AND FIRST TASK OF MONOMER COMMUNICATE___
        IF(THISTASK.NE.ITASK.AND.THISTASK.NE.1) CYCLE 
        THIS1=>THIS
!       __LOOP OVER ALL INSTANCES THAT HAVE ITASK AS MAIN TASK__________________
        DO I=1,ICOUNT(ITASK)
          IF(THISTASK.EQ.ITASK) THEN
!           __LOOP UNTIL AN INSTANCE WITH MAIN TASK ITASK HAS BEEN IDENTIFIED___
            DO 
              CALL MPE$QUERY(THIS1%CID,NTASKS1,THISTASK1)
              IF(THISTASK1.EQ.1) EXIT
              THIS1=>THIS1%NEXT
            ENDDO
            ID=THIS1%ID
            TINV=THIS1%TINV
            TSUPER=THIS1%TSUPER
            KVEC(:)=THIS1%KVEC
            NR1=THIS1%NR1
            NR2=THIS1%NR2
            NR3=THIS1%NR3
            CALL MPE$QUERY(THIS1%CID,NTASKS1,THISTASK1)
          END IF
!         __COMMUNICATE INFORMATION TO THE FIRST TASK OF MONOMER________________
          CALL MPE$SENDRECEIVE('MONOMER',ITASK,1,ID)
          CALL MPE$SENDRECEIVE('MONOMER',ITASK,1,TINV)
          CALL MPE$SENDRECEIVE('MONOMER',ITASK,1,TSUPER)
          CALL MPE$SENDRECEIVE('MONOMER',ITASK,1,KVEC)
          CALL MPE$SENDRECEIVE('MONOMER',ITASK,1,NR1)
          CALL MPE$SENDRECEIVE('MONOMER',ITASK,1,NR2)
          CALL MPE$SENDRECEIVE('MONOMER',ITASK,1,NR3)
          CALL MPE$SENDRECEIVE('MONOMER',ITASK,1,NTASKS1)
          ALLOCATE(IWORK(NTASKS1,3))
          IF(THISTASK.EQ.ITASK) THEN
            IWORK(:,1)=THIS1%NR1LARR(:)          
            IWORK(:,2)=THIS1%NGLARR(:)          
            IWORK(:,3)=THIS1%NSTRIPELARR(:)          
          END IF
          CALL MPE$SENDRECEIVE('MONOMER',ITASK,1,IWORK)
!         __WRITE INFORMATION___________________________________________________
          IF(THISTASK.EQ.1) THEN
            CALL REPORT$TITLE(NFIL,'PLANEWAVE')
            CALL REPORT$CHVAL(NFIL,'ID',ID)
            CALL REPORT$L4VAL(NFIL,'TINV',TINV)
            CALL REPORT$L4VAL(NFIL,'TSUPER',TSUPER)
            WRITE(NFIL,FMT='("GBASINV*K=(",3F10.5,")")')KVEC
            CALL REPORT$I4VAL(NFIL,'NR1',NR1,' ')
            CALL REPORT$I4VAL(NFIL,'NR2',NR2,' ')
            CALL REPORT$I4VAL(NFIL,'NR3',NR3,' ')
            CALL REPORT$I4VAL(NFIL,'BELONGS TO TASK GROUP '//TRIM(THIS1%CID)//' WITH MAIN TASK',ITASK,' ')
            DO J=1,NTASKS1
              WRITE(NFIL,FMT='("SUB-TASK ",I4," NR1L ",I5," NGL ",I10," NSTRIPEL ",I5)') &
              J,IWORK(J,1),IWORK(J,2),IWORK(J,3)
            ENDDO
          END IF
          DEALLOCATE(IWORK)
          IF(THISTASK.EQ.ITASK)THIS1=>THIS1%NEXT
        ENDDO
      ENDDO
      DEALLOCATE(ICOUNT)
      RETURN
      END
!     
!     .................................................................
      SUBROUTINE PLANEWAVE$SELECT(ID)
!     ******************************************************************
!     **  PLANEWAVE$SELECT                                            **
!     **                                                              **
!     **  SELECTS                                                     **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      TYPE(PWPARALLEL_TYPE),POINTER :: THIS1
!     ******************************************************************
      IF(.NOT.ASSOCIATED(THIS))THEN
        CALL ERROR$MSG('LIST NOT INITIALIZED')
        CALL ERROR$STOP('PLANEWAVE$SELECT')
      ENDIF
!
!     ==================================================================
!     ==  CHECK IF ALREADY SELECTED                                   ==
!     ==================================================================
      IF(THIS%ID.EQ.ID) RETURN
!
!     ==================================================================
!     ==  SCROLL TO THE PROPER OBJECT                                 ==
!     ==================================================================
      THIS1=>THIS
      DO WHILE(THIS%ID.NE.ID) 
        THIS=>THIS%NEXT 
        IF(ASSOCIATED(THIS,THIS1)) THEN
          CALL ERROR$MSG('PWPARALLELOBJECT NOT FOUND')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PLANEWAVE$SELECT')
        END IF
      ENDDO      
      RETURN
      END
!     
!     .................................................................
      SUBROUTINE PLANEWAVE_NEW(ID)
!     ******************************************************************
!     **  PLANEWAVE_NEW                                               **
!     **                                                              **
!     **  INSTANCES FORM A RING, THAT IS THE FIRST INSTANCE IS        **
!     **  THIS%NEXT OF THE LAST INSTANCE                              **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*)      ,INTENT(IN) :: ID
      TYPE(PWPARALLEL_TYPE),POINTER :: THIS1
!     ***************UCH***************************************************
      IF(.NOT.TINI)THEN
        TINI=.TRUE.
        ALLOCATE(THIS)
        THIS%NEXT=>THIS
      ELSE
        ALLOCATE(THIS1)
        THIS1%NEXT=>THIS%NEXT
        THIS%NEXT=>THIS1
        THIS=>THIS%NEXT
      ENDIF
!
!     ==================================================================
!     ==  DO SOME INITIALIZATIONS                                     ==
!     ==================================================================
      THIS%ID=ID
      THIS%CID='NONE'
      THIS%TSUPER=.FALSE.
      THIS%NGAMMA=0
      NULLIFY(THIS%IGVECH)
!
!     ==================================================================
!     ==  CHECK IF ID ALREADY EXISTS                                  ==
!     ==================================================================
      THIS1=>THIS%NEXT
      DO WHILE(.NOT.ASSOCIATED(THIS,THIS1))
        IF(THIS%ID.EQ.THIS1%ID) THEN
          CALL ERROR$MSG('ID ALREADY EXISTS')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PLANEWAVE_NEW')
        END IF
        THIS1=>THIS1%NEXT 
      ENDDO      
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$GETL4(ID,VAL)
!     ******************************************************************
!     **  HANDLES REQUESTS FOR LOGICAL VARIALES                       **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VAL
!     ******************************************************************
      IF(ID.EQ.'SUPER') THEN ! SELECTION OF SUPER WAVE FUNCTION INTERFACE
        VAL=THIS%TSUPER
      ELSE IF(ID.EQ.'TINV') THEN ! TIME INVERSION SYMMETRY
        VAL=THIS%TINV
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PLANEWAVE$GETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$SETL4(ID,VAL)
!     ******************************************************************
!     **  ACCEPTS LOGICAL NUMBERS                                     **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'SUPER') THEN ! SELECT SUPER WAVE FUNCTION INTERFACE
        THIS%TSUPER=VAL
        IF(VAL.AND.DOT_PRODUCT(THIS%KVEC,THIS%KVEC).GT.1.D-6) THEN
          CALL ERROR$MSG('SUPER WAVE FUNCTION INTERFACE CAN BE SELECTED')
          CALL ERROR$MSG('ONLY FOR THE GAMMA POINT')
          CALL ERROR$CHVAL(ID,ID)
          CALL ERROR$STOP('PLANEWAVE$SETL4')
        END IF
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PLANEWAVE$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$GETI4(ID,VAL)
!     ******************************************************************
!     **  HANDLES REQUESTS FOR INTEGER NUMBERS                        **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
      LOGICAL(4)              :: TSUPER
      INTEGER(4)              :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
      TSUPER=THIS%TSUPER
      IF(ID.EQ.'NGL') THEN         ! LOCAL #(GVECTORS)
        IF(TSUPER) THEN
          VAL=THIS%NGHARR(THISTASK)
        ELSE
          VAL=THIS%NGLARR(THISTASK)
        END IF
      ELSE IF(ID.EQ.'NGG') THEN    ! GLOBAL #(GVECTORS)
        IF(TSUPER) THEN
          VAL=SUM(THIS%NGHARR)
        ELSE
          VAL=SUM(THIS%NGLARR)
        END IF
      ELSE IF(ID.EQ.'NRL') THEN    ! LOCAL #(R-POINTS)
        VAL=THIS%NR1L*THIS%NR2*THIS%NR3
      ELSE IF(ID.EQ.'NGAMMA') THEN ! INDEX OF THE GAMMA POINT OR ZERO
        IF(TSUPER) THEN
          VAL=THIS%NGAMMA
        ELSE
          CALL ERROR$MSG('NO SUPER WAVE FUNCTION INTERFACE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PLANEWAVE$GETI4')
        END IF
      ELSE IF(ID.EQ.'NR1') THEN    ! #(R-POINTS ALONG 1ST LATTICE VECTOR) 
        VAL=THIS%NR1
      ELSE IF(ID.EQ.'NR1L') THEN   ! #(LOCAL R-POINTS ALONG 1ST LATTICE VECTOR)
        VAL=THIS%NR1L
      ELSE IF(ID.EQ.'NR2') THEN    ! #(R-POINTS ALONG 2ND LATTICE VECTOR) 
        VAL=THIS%NR2
      ELSE IF(ID.EQ.'NR3') THEN    ! #(R-POINTS ALONG 3RD LATTICE VECTOR) 
        VAL=THIS%NR3
      ELSE IF(ID.EQ.'NR1START') THEN  !#(1ST LOCAL R-POINT ALONG 1ST LATTICE VECTOR)
        VAL=THIS%NR1START
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PLANEWAVE$GETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$GETI4A(ID,LEN,VAL)
!     ******************************************************************
!     **  HANDLES REQUESTS FOR REAL ARRAYS                            **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(OUT):: VAL(LEN)
      INTEGER(4)              :: NGL
      LOGICAL(4)              :: SUPER
      INTEGER(4)              :: IG,IGH
      INTEGER(4)              :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
      SUPER=THIS%TSUPER
!
!     ==================================================================
!     ==  GET G-VECTORS                                               ==
!     ==================================================================
      IF(ID.EQ.'IGVEC') THEN   ! G-VECTORS
        NGL=THIS%NGLARR(THISTASK)
        IF(.NOT.SUPER) THEN
          IF(LEN.NE.3*NGL) THEN
            CALL ERROR$MSG('SIZE INCONSISTENT')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$I4VAL('LEN',LEN)
            CALL ERROR$I4VAL('3*NGL',3*NGL)
            CALL ERROR$STOP('PLANEWAVE$GETI4A')
          END IF
          VAL=RESHAPE(THIS%IGVEC,(/3*NGL/))
        ELSE
          IGH=0
          DO IG=1,NGL
            IF(THIS%MINUSG(IG).LT.IG) CYCLE
            IF(IGH+3.GT.LEN) THEN
              CALL ERROR$MSG('SIZE INCONSISTENT')
              CALL ERROR$CHVAL('ID',ID)
              CALL ERROR$I4VAL('LEN',LEN)
              CALL ERROR$I4VAL('1+3*(NGL-1)/2',1+3*(NGL-1)/2)
              CALL ERROR$STOP('PLANEWAVE$GETI4A')
            END IF
            VAL(IGH+1)=THIS%IGVEC(1,IG)
            VAL(IGH+2)=THIS%IGVEC(2,IG)
            VAL(IGH+3)=THIS%IGVEC(3,IG)
            IGH=IGH+3
          ENDDO
          IF(IGH.NE.LEN) THEN
            CALL ERROR$MSG('SIZE INCONSISTENT')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$I4VAL('LEN',LEN)
            CALL ERROR$I4VAL('IGH',IGH)
            CALL ERROR$STOP('PLANEWAVE$GETI4A')
          END IF
        END IF
!
!     ==================================================================
!     ==  NR1LARR                                                     ==
!     ==================================================================
      ELSE IF(ID.EQ.'NR1LARR') THEN
        IF(LEN.NE.NTASKS) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NTASKS',NTASKS)
          CALL ERROR$STOP('PLANEWAVE$GETI4A')
        END IF
        VAL(:)=THIS%NR1LARR(:)
!
!     ==================================================================
!     ==  UNKNOWN ID                                                  ==
!     ==================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PLANEWAVE$GETI4A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$GETR8(ID,VAL)
!     ******************************************************************
!     **  HANDLES REQUESTS FOR REAL NUMBERS                           **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
!     ******************************************************************
      IF(ID.EQ.'RWEIGHT') THEN   
        VAL=THIS%RWEIGHT
      ELSE IF(ID.EQ.'GWEIGHT') THEN   
        VAL=THIS%GWEIGHT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PLANEWAVE$GETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$GETR8A(ID,LEN,VAL)
!     ******************************************************************
!     **  HANDLES REQUESTS FOR REAL ARRAYS                            **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(OUT):: VAL(LEN)
      INTEGER(4)              :: NGL
      REAL(8)     ,ALLOCATABLE:: GVEC(:,:)
      LOGICAL(4)              :: SUPER
      INTEGER(4)              :: IG,IGH
      INTEGER(4)              :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
      SUPER=THIS%TSUPER
!
!     ==================================================================
!     ==  GET G-VECTORS                                               ==
!     ==================================================================
      IF(ID.EQ.'GVEC') THEN   ! G-VECTORS
        NGL=THIS%NGLARR(THISTASK)
        IF(.NOT.SUPER) THEN
          IF(LEN.NE.3*NGL) THEN
            CALL ERROR$MSG('SIZE INCONSISTENT')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$I4VAL('LEN',LEN)
            CALL ERROR$I4VAL('3*NGL',3*NGL)
            CALL ERROR$STOP('PLANEWAVE$GETR8A')
          END IF
          CALL PLANEWAVE_GVECTORS(NGL,VAL)
        ELSE
          ALLOCATE(GVEC(3,NGL))
          CALL PLANEWAVE_GVECTORS(NGL,GVEC)
          IGH=0
          DO IG=1,NGL
            IF(THIS%MINUSG(IG).LT.IG) CYCLE
            IF(IGH+3.GT.LEN) THEN
              CALL ERROR$MSG('SIZE INCONSISTENT')
              CALL ERROR$CHVAL('ID',ID)
              CALL ERROR$I4VAL('LEN',LEN)
              CALL ERROR$I4VAL('1+3*(NGL-1)/2',1+3*(NGL-1)/2)
              CALL ERROR$STOP('PLANEWAVE$GETR8A')
            END IF
            VAL(IGH+1)=GVEC(1,IG)
            VAL(IGH+2)=GVEC(2,IG)
            VAL(IGH+3)=GVEC(3,IG)
            IGH=IGH+3
          ENDDO
          DEALLOCATE(GVEC)
          IF(IGH.NE.LEN) THEN
            CALL ERROR$MSG('SIZE INCONSISTENT')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$I4VAL('LEN',LEN)
            CALL ERROR$I4VAL('IGH',IGH)
            CALL ERROR$STOP('PLANEWAVE$GETR8A')
          END IF
        END IF
!
!     ==================================================================
!     ==  GET G**2                                                    ==
!     ==================================================================
      ELSE IF(ID.EQ.'G2') THEN   ! G-VECTORS
        NGL=THIS%NGLARR(THISTASK)
        ALLOCATE(GVEC(3,NGL))
        CALL PLANEWAVE_GVECTORS(NGL,GVEC)
        IF(.NOT.SUPER) THEN
          IF(LEN.NE.NGL) THEN
            CALL ERROR$MSG('SIZE INCONSISTENT')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('PLANEWAVE$GETR8A')
          END IF
          DO IG=1,NGL
            VAL(IG)=GVEC(1,IG)**2+GVEC(2,IG)**2+GVEC(3,IG)**2
          ENDDO
        ELSE
          IGH=0
          DO IG=1,NGL
            IF(THIS%MINUSG(IG).LT.IG) CYCLE
            IGH=IGH+1
            IF(IGH.GT.LEN) THEN
              CALL ERROR$MSG('SIZE INCONSISTENT')
              CALL ERROR$CHVAL('ID',ID)
              CALL ERROR$STOP('PLANEWAVE$GETR8A')
            END IF
            VAL(IGH)=GVEC(1,IG)**2+GVEC(2,IG)**2+GVEC(3,IG)**2
          ENDDO
          IF(IGH.NE.LEN) THEN
            CALL ERROR$MSG('SIZE INCONSISTENT')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('PLANEWAVE$GETR8A')
          END IF
        END IF
        DEALLOCATE(GVEC)
!
!     ==================================================================
!     ==  K-POINT                                                     ==
!     ==================================================================
      ELSE IF(ID.EQ.'K') THEN
        IF(3.NE.LEN) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PLANEWAVE$GETR8A')
        END IF
        VAL(:)=MATMUL(THIS%GBAS,THIS%KVEC)
!
!     ==================================================================
!     ==  K-POINT IN RELATIVE COORDINATES                             ==
!     ==================================================================
      ELSE IF(ID.EQ.'XK') THEN
        IF(3.NE.LEN) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PLANEWAVE$GETR8A')
        END IF
        VAL(:)=THIS%KVEC
!
!     ==================================================================
!     ==  RECIPROCAL LATTICE VECTORS                                  ==
!     ==================================================================
      ELSE IF(ID.EQ.'GBAS') THEN
        IF(9.NE.LEN) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PLANEWAVE$GETR8A')
        END IF
        VAL(:)=RESHAPE(THIS%GBAS,(/9/))
!
!     ==================================================================
!     ==  UNKNOWN ID                                                  ==
!     ==================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PLANEWAVE$GETR8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$SETR8A(ID,LEN,VAL)
!     ******************************************************************
!     **  ACCEPTS REAL ARRAYS                                         **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(IN) :: VAL(LEN)
      REAL(8)                 :: RBAS(3,3)
!     ******************************************************************
!
!     ==================================================================
!     ==  SET G-VECTORS                                               ==
!     ==================================================================
      IF(ID.EQ.'RBAS') THEN   ! LATTICE VECTORS
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('PLANEWAVE$SETR8A')
        END IF
        CALL GBASS(VAL,THIS%GBAS,THIS%GWEIGHT)
        THIS%RWEIGHT=THIS%GWEIGHT/(THIS%NR1*THIS%NR2*THIS%NR3)
!
!     ==================================================================
!     ==  SET RECIPROCAL LATTICE VECTORS                              ==
!     ==================================================================
      ELSE IF(ID.EQ.'GBAS') THEN   ! LATTICE VECTORS
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('PLANEWAVE$SETR8A')
        END IF
        CALL GBASS(VAL,RBAS,THIS%GWEIGHT)
        CALL GBASS(RBAS,THIS%GBAS,THIS%GWEIGHT)
        THIS%RWEIGHT=THIS%GWEIGHT/(THIS%NR1*THIS%NR2*THIS%NR3)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PLANEWAVE$SETR8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$GETC8A(ID,LEN,VAL)
!     ******************************************************************
!     **  HANDLES REQUESTS FOR COMPLEX ARRAYS                         **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      COMPLEX(8)  ,INTENT(OUT):: VAL(LEN)
!     ******************************************************************
!
!     ==================================================================
!     ==  GET PHASE FACTOR EIKR                                       ==
!     ==================================================================
      IF(ID.EQ.'EIKR') THEN   ! PHASE FACTOR
        IF(LEN.NE.THIS%NR1L*THIS%NR2*THIS%NR3) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PLANEWAVE$GETC8A')
        END IF
        CALL PLANEWAVE_EIKR(THIS%NR1,THIS%NR2,THIS%NR3 &
     &                     ,THIS%NR1START,THIS%NR1START+THIS%NR1L-1 &
     &                     ,THIS%KVEC,VAL)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PLANEWAVE$GETC8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$GETCH(ID,VAL)
!     ******************************************************************
!     **  HANDLES REQUESTS FOR COMPLEX ARRAYS                         **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(OUT):: VAL
!     ******************************************************************
!
!     ==================================================================
!     ==  GET CURRENT ID                                              ==
!     ==================================================================
      IF(ID.EQ.'ID') THEN   ! PHASE FACTOR
        VAL=THIS%ID
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PLANEWAVE$GETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$SETCH(ID,VAL)
!     ******************************************************************
!     **  HANDLES REQUESTS FOR COMPLEX ARRAYS                         **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN):: VAL
!     ******************************************************************
!
!     ==================================================================
!     ==  GET CURRENT ID                                              ==
!     ==================================================================
      IF(ID.EQ.'CID') THEN   ! PHASE FACTOR
        THIS%CID=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PLANEWAVE$GETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$CHECKINVERSION(KVEC,TINV)
!     ******************************************************************
!     **  THIS TOOL IS USED TO DETERMINE IF THE RECIPROCAL SPACE GRID **
!     **  IS INVERSION SYMETRIC                                       **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: KVEC(3)! K-POINT IN RELATIVE COORDINATES
      LOGICAL(4),INTENT(OUT) :: TINV
      REAL(8)   ,PARAMETER   :: DSMALL=1.D-12
      REAL(8)                :: SVAR
      INTEGER(4)             :: I
!     ******************************************************************
      TINV=.TRUE.
      DO I=1,3
        SVAR=MODULO(KVEC(I),0.5D0)
        TINV=TINV.AND.(ABS(SVAR).LT.DSMALL)
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$DIVIDERGRIDONTASKS(EPW,RBAS,CID,NR1START,NR1L,NR2,NR3)
!     ******************************************************************
!     **  THIS TOOL IS USED TO DETERMINE THE REAL SPACE GRID,         **
!     **  WHICH IS THE INPUT FOR THE PLANEWAVE OBJECT.                **
!     **  THE REAL SPACE GRID MUST BE DEFINED AHEAD OF TIME BECAUSE   **
!     **  DIFFERENT PLANEWAVE INSTANCES COMMUNICATE VIA THE REAL SPACE**
!     **  GRID                                                        **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: CID      ! COMMUNICATOR ID (SEE MPE)
      REAL(8)   ,INTENT(IN)  :: EPW      ! PLANE WAVE CUTOFF
      REAL(8)   ,INTENT(IN)  :: RBAS(3,3)! R-SPACE LATTICE VECTORS
      INTEGER(4),INTENT(OUT) :: NR1START ! INDEX OF FIRST R-GRID PLANE
      INTEGER(4),INTENT(OUT) :: NR1L     ! #(R-GRID PLANES ON THIS TASK)
      INTEGER(4),INTENT(OUT) :: NR2      ! #(R-POINTS IN SECOND DIRECTION)
      INTEGER(4),INTENT(OUT) :: NR3      ! #(R-POINTS IN THIRD DIRECTION)
      REAL(8)                :: GBAS(3,3)
      REAL(8)                :: VOL
      REAL(8)                :: RAD
      INTEGER(4)             :: NR1G
      INTEGER(4)             :: NTASKS,THISTASK
      INTEGER(4)             :: ISVAR
      INTEGER(4)             :: MIN1,MAX1
      INTEGER(4)             :: MIN2,MAX2
      INTEGER(4)             :: MIN3,MAX3
!     ******************************************************************
      CALL MPE$QUERY(CID,NTASKS,THISTASK)
      CALL GBASS(RBAS,GBAS,VOL)
      RAD=SQRT(2.D0*EPW)+1.D-6
      CALL BOXSPH(GBAS,0.D0,0.D0,0.D0,RAD &
     &           ,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
      NR1G=MAX1-MIN1+1
      NR2=MAX2-MIN2+1
      NR3=MAX3-MIN3+1
      CALL LIB$FFTADJUSTGRD(NR1G)
      CALL LIB$FFTADJUSTGRD(NR2)
      CALL LIB$FFTADJUSTGRD(NR3)
      NR1L=NR1G/NTASKS
      ISVAR=NR1G-NTASKS*NR1L
      IF(THISTASK.LE.ISVAR) THEN
        NR1L=NR1L+1
        NR1START=(THISTASK-1)*NR1L+1
      ELSE
        NR1START=(THISTASK-1)*NR1L+ISVAR+1
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$INITIALIZE(ID,CID,RBAS,KVEC,TINV,EPW &
     &                               ,NR1START,NR1L,NR2,NR3)
!     ******************************************************************
!     **                                                              **
!     **  REMARK:                                                     **
!     **  - THE DISTRIBUTION OF PLANES ON PARALLEL TASKS HAS TO BE    **
!     **    DONE AHEAD OF TIME                                        **
!     **    USE TOOL PLANEWAVE$DIVIDERGRIDONTASKS                     **
!     **                                                              **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID        ! IDENTIFIER
      CHARACTER(*),INTENT(IN) :: CID       ! COMMUNICATOR ID (SEE MPE)
      REAL(8)     ,INTENT(IN) :: RBAS(3,3) ! LATTICE VECTORS 
      REAL(8)     ,INTENT(IN) :: KVEC(3)   ! K-VECTOR IN RELATIVE! COORDINATES
      LOGICAL(4)  ,INTENT(IN) :: TINV      ! SELECT SUPERWAVEFUNCTIONS
      REAL(8)     ,INTENT(IN) :: EPW       ! PLANE WAVE CUTOFF
      INTEGER(4)  ,INTENT(IN) :: NR1START  ! FIRST PLANE OF R-GRID POINTS
      INTEGER(4)  ,INTENT(IN) :: NR1L,NR2,NR3  ! REAL SPACE GRID
      REAL(8)                 :: GBAS(3,3)
      REAL(8)                 :: VOL
      INTEGER(4)              :: NG
      INTEGER(4)              :: NGH
      INTEGER(4)              :: NGL
      INTEGER(4)              :: I1,I2,I3,IG,IGL,ITASK,ISTRIPE,IGH
      INTEGER(4)              :: NSTRIPELX
      INTEGER(4)              :: NR1
      INTEGER(4) ,ALLOCATABLE :: IGVEC(:,:)
      INTEGER(4) ,ALLOCATABLE :: MINUSG(:)
      INTEGER(4) ,ALLOCATABLE :: ITASKOFYZ(:,:)
      INTEGER(4) ,ALLOCATABLE :: MAP(:)
      INTEGER(4)              :: NTASKS,THISTASK
!     ******************************************************************
                           CALL TRACE$PUSH('PLANEWAVE$INITIALIZE')
      IF(NR1L.EQ.0) THEN
        CALL ERROR$MSG('ZERO GRID POINTS IN FIRST DIRECTION OF THE REAL SPACE GRID')
        CALL ERROR$MSG('STOPPING, BECAUSE PARALLELIZATION WILL FAIL')
        CALL ERROR$STOP('PLANEWAVE$INITIALIZE')
      ENDIF 
      CALL MPE$QUERY(CID,NTASKS,THISTASK) 
      CALL PLANEWAVE_NEW(ID)
      CALL GBASS(RBAS,GBAS,VOL)
      THIS%CID=CID
      THIS%KVEC=KVEC
      THIS%GBAS=GBAS
      THIS%NR2=NR2
      THIS%NR3=NR3
      THIS%NR1L=NR1L
      THIS%NR1START=NR1START
      THIS%TINV=TINV
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
      ALLOCATE(THIS%NR1LARR(NTASKS))
      CALL MPE$GATHER(THIS%CID,1,NR1L,THIS%NR1LARR)
      CALL MPE$BROADCAST(THIS%CID,1,THIS%NR1LARR)
      THIS%NR1=SUM(THIS%NR1LARR)
      NR1=THIS%NR1
      THIS%GWEIGHT=VOL
      THIS%RWEIGHT=VOL/REAL(THIS%NR1*THIS%NR2*THIS%NR3)
!
!     ===================================================================
!     == CALCULATE G-VECTORS                                           ==
!     ===================================================================
      CALL PLANEWAVE_COUNTG(GBAS,KVEC,EPW,NG)
      ALLOCATE(IGVEC(3,NG))
      CALL PLANEWAVE_GETIG(GBAS,KVEC,EPW,NG,IGVEC)
      ALLOCATE(MINUSG(NG))
      IF(TINV) THEN
        CALL PLANEWAVE$MINUSG(KVEC,NG,IGVEC,MINUSG)
      ELSE
        MINUSG(:)=0
      END IF
!
!     ===================================================================
!     == DISTRIBUTE G-VECTORS ONTO TASKS                               ==
!     ===================================================================
      ALLOCATE(ITASKOFYZ(NR2,NR3))
      CALL PLANEWAVE_MAPSTRIPESTOTASKS(NR1,NR2,NR3 &
     &            ,NG,IGVEC,TINV,MINUSG,NTASKS,ITASKOFYZ)
!
!     ===================================================================
!     == DETERMINE NUMBER OF STRIPES                                   ==
!     ===================================================================
      ALLOCATE(THIS%NSTRIPELARR(NTASKS))
      THIS%NSTRIPELARR(:)=0
      DO I2=1,NR2
        DO I3=1,NR3
          ITASK=ITASKOFYZ(I2,I3)
          IF(ITASK.EQ.0) CYCLE
          THIS%NSTRIPELARR(ITASK)=THIS%NSTRIPELARR(ITASK)+1
        ENDDO
      ENDDO
!
!     ===================================================================
!     == DETERMINE #(LOCAL G-VECTORS)                                  ==
!     ===================================================================
      ALLOCATE(THIS%NGLARR(NTASKS))
      THIS%NGLARR(:)=0
      DO IG=1,NG
        I2=MODULO(IGVEC(2,IG),NR2)+1   !<
        I3=MODULO(IGVEC(3,IG),NR3)+1   !<
        ITASK=ITASKOFYZ(I2,I3)
        THIS%NGLARR(ITASK)=THIS%NGLARR(ITASK)+1
      ENDDO
!
!     ===================================================================
!     == FOLD DOWN ARRAYS                                              ==
!     ===================================================================
      NGL=THIS%NGLARR(THISTASK)
      IF(TINV) THEN
        ALLOCATE(THIS%MINUSG(NGL))
      ELSE
        NULLIFY(THIS%MINUSG)
      END IF
      ALLOCATE(THIS%IGVEC(3,NGL))
      ALLOCATE(THIS%TASKOFIG(NG))
      ALLOCATE(MAP(NG))
      MAP(:)=0
      IGL=0
      IGH=0
      DO IG=1,NG
        I2=MODULO(IGVEC(2,IG),NR2)+1
        I3=MODULO(IGVEC(3,IG),NR3)+1
        ITASK=ITASKOFYZ(I2,I3)
        THIS%TASKOFIG(IG)=ITASK
        IF(ITASK.EQ.0) THEN
          CALL ERROR$MSG('FAILURE IN TASK ASSIGNMENT')
          CALL ERROR$CHVAL('ID',THIS%ID)
          CALL ERROR$I4VAL('IG',IG)
          CALL ERROR$I4VAL('I2',I2)
          CALL ERROR$I4VAL('I3',I3)
          CALL ERROR$STOP('PLANEWAVE$INITIALIZE')
        END IF
        IF(ITASK.NE.THISTASK) CYCLE
        IGL=IGL+1
        THIS%IGVEC(:,IGL)=IGVEC(:,IG)
        IF(TINV) THEN
          THIS%MINUSG(IGL)=MINUSG(IG)
          MAP(IG)=IGL                   ! USED TO TRANSLATE VALUE OF MINUSG
          IF(MINUSG(IG).GE.IG) IGH=IGH+1
        END IF
      ENDDO
      DEALLOCATE(MINUSG)
      DEALLOCATE(IGVEC)
      IF(THISTASK.NE.1) DEALLOCATE(THIS%TASKOFIG)
      IF(THIS%TINV) THEN
        DO IG=1,NGL
          IF(MAP(THIS%MINUSG(IG)).EQ.0) THEN
            CALL ERROR$MSG('MINUSG MUST BE ON THE SAME TASK AS G')
            CALL ERROR$STOP('PLANEWAVE$INITIALIZE')
          END IF
          THIS%MINUSG(IG)=MAP(THIS%MINUSG(IG))
        ENDDO
      END IF
      DEALLOCATE(MAP)
!
!     ===================================================================
!     == MAPPING FROM ISTRIPE TO YZ                                    ==
!     ===================================================================
      NSTRIPELX=MAXVAL(THIS%NSTRIPELARR)
      ALLOCATE(THIS%ISTRIPETOYZ(NSTRIPELX,NTASKS))
      DO ITASK=1,NTASKS
        ISTRIPE=0
        DO I2=1,NR2
          DO I3=1,NR3
            IF(ITASK.NE.ITASKOFYZ(I2,I3)) CYCLE
            ISTRIPE=ISTRIPE+1
            THIS%ISTRIPETOYZ(ISTRIPE,ITASK)=I2+NR2*(I3-1)
          ENDDO
        ENDDO
      ENDDO
!
!     ===================================================================
!     == MAPPING FROM IG TO IPLANE*ISTRIPE                             ==
!     ===================================================================
      ALLOCATE(THIS%IGTOSTRIPE(NGL))
!     == CONVERT ITASKOFYZ TO YZ->ISTRIPE
      ISTRIPE=0
      DO I2=1,NR2
        DO I3=1,NR3
          ITASK=ITASKOFYZ(I2,I3)
          IF(ITASK.EQ.THISTASK) THEN
            ISTRIPE=ISTRIPE+1
            ITASKOFYZ(I2,I3)=ISTRIPE
          ELSE
            ITASKOFYZ(I2,I3)=0
          END IF
        ENDDO
      ENDDO
!
      DO IG=1,NGL
        I1=MODULO(THIS%IGVEC(1,IG),NR1)+1
        I2=MODULO(THIS%IGVEC(2,IG),NR2)+1
        I3=MODULO(THIS%IGVEC(3,IG),NR3)+1
        ISTRIPE=ITASKOFYZ(I2,I3)
        IF(ISTRIPE.EQ.0) THEN
          CALL ERROR$MSG('FAILED CONSISTENCY CHECK "ISTRIPE=0"')
          CALL ERROR$STOP('PLANEWAVE$INITIALIZE')
        END IF
        THIS%IGTOSTRIPE(IG)=I1+NR1*(ISTRIPE-1)
      ENDDO
      DEALLOCATE(ITASKOFYZ)
!
!     ===================================================================
!     ==  SOME SPECIAL THINGS FOR SUPER WAVE FUNCTIONS                 ==
!     ===================================================================
      THIS%NGAMMA=0
      IF(TINV) THEN
        NGH=0
        DO IG=1,NGL
          IF(THIS%MINUSG(IG).LT.IG) CYCLE
          NGH=NGH+1
          IF(THIS%MINUSG(IG).EQ.IG) THEN
            IF(THIS%NGAMMA.EQ.0) THEN
              THIS%NGAMMA=NGH
            ELSE
              CALL ERROR$MSG(' MORE THAN ONE GAMMA POINT?')
              CALL ERROR$STOP('PLANEWAVE$INITIALIZE')
            END IF
          END IF
        ENDDO
        ALLOCATE(THIS%NGHARR(NTASKS))
        CALL MPE$GATHER(THIS%CID,1,NGH,THIS%NGHARR)
        CALL MPE$BROADCAST(THIS%CID,1,THIS%NGHARR)
!
        ALLOCATE(THIS%IGVECH(3,NGL))
        IGH=0
        DO IG=1,NGL
          IF(THIS%MINUSG(IG).LT.IG) CYCLE
          IGH=IGH+1
          THIS%IGVECH(:,IGH)=THIS%IGVEC(:,IG)
        ENDDO
      ELSE
        NULLIFY(THIS%NGHARR)
      END IF
                             CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE_COUNTG(GBAS,KVEC,EPW,NG)
!     ******************************************************************
!     **  COUNTS NUMBER OF G-VECTORS                                  **
!     **  G=KVEC+GBAS*IG                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: GBAS(3,3)
      REAL(8)   ,INTENT(IN)  :: KVEC(3)
      REAL(8)   ,INTENT(IN)  :: EPW
      INTEGER(4),INTENT(OUT) :: NG
      REAL(8)                :: RAD
      REAL(8)                :: G1,G2,G3
      REAL(8)                :: T1,T2,T3
      INTEGER(4)             :: MIN1,MIN2,MIN3,MAX1,MAX2,MAX3
      INTEGER(4)             :: I,J,K
      REAL(8)                :: KARTK(3)
!     ******************************************************************
      RAD=SQRT(2.D0*EPW)+1.D-6
      KARTK(:)=MATMUL(GBAS,KVEC)
      CALL BOXSPH(GBAS,-KARTK(1),-KARTK(2),-KARTK(3),RAD &
     &           ,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
!     ==================================================================
!     ==  FIND ALL LATTICE VECTORS WITH 0.5*G**2 < EPW                ==
!     ==================================================================
      NG=0
      DO I=MIN1,MAX1
        T1=REAL(I,KIND=8)
        DO J=MIN2,MAX2
          T2=REAL(J,KIND=8)
          DO K=MIN3,MAX3
            T3=REAL(K,KIND=8)
            G1=KARTK(1)+GBAS(1,1)*T1+GBAS(1,2)*T2+GBAS(1,3)*T3
            G2=KARTK(2)+GBAS(2,1)*T1+GBAS(2,2)*T2+GBAS(2,3)*T3
            G3=KARTK(3)+GBAS(3,1)*T1+GBAS(3,2)*T2+GBAS(3,3)*T3
            G2=G1**2+G2**2+G3**2
            IF(0.5D0*G2.GT.EPW) CYCLE
            NG=NG+1
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$MINUSG(KVEC,NG,IG,MINUSG)
!     ******************************************************************
!     **  EVALUATES THE POINTERS BETWEEN INVERSION SYMMETRIC PARTNERS **
!     **  IN G-SPACE                                                  **
!     **  G=KVEC+GBAS*IG                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: KVEC(3)
      INTEGER(4),INTENT(IN)  :: NG
      INTEGER(4),INTENT(IN)  :: IG(3,NG)
      INTEGER(4),INTENT(OUT) :: MINUSG(NG)
      REAL(8)   ,PARAMETER   :: SMALL=1.D-5
      LOGICAL(4)             :: TCHK
      INTEGER(4)             :: IK(3)
      INTEGER(4)             :: I,IG1,IG2
      INTEGER(4)             :: IGM(3)
!     ******************************************************************
!     ==================================================================
!     ==  CHECK IF LATTICE IS INVERSION SYMMETRIC                     ==
!     ==================================================================
      DO I=1,3
        IK(I)=NINT(2.D0*KVEC(I))
        TCHK=(MODULO(2.D0*KVEC(I),1.D0).LT.SMALL)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('SUPER WAVE FUNCTIONS SELECTED')
          CALL ERROR$MSG('FOR A NON-INVERSION SYMMETRIC RECIPROCAL GRID')
          CALL ERROR$STOP('PLANEWAVE$MINUSG')
        END IF
      ENDDO
!
!     ==================================================================
!     ==   MAP G-VECTORS ONTO INVERSION SYMMETRIC PARTNERS            ==
!     ==================================================================
      IK(:)=-IK(:)
      MINUSG(:)=0
      DO IG1=1,NG
        IF(MINUSG(IG1).NE.0) CYCLE
        IGM(:)=IK(:)-IG(:,IG1)        
        DO IG2=IG1,NG
          IF(IGM(1).NE.IG(1,IG2)) CYCLE
          IF(IGM(2).NE.IG(2,IG2)) CYCLE
          IF(IGM(3).NE.IG(3,IG2)) CYCLE
          MINUSG(IG1)=IG2
          MINUSG(IG2)=IG1
          GOTO 100
        ENDDO 
        CALL ERROR$MSG('INVERSION SYMMETRIC PARTNER NOT FOUND')
        CALL ERROR$STOP('PLANEWAVE$MINUSG')
 100    CONTINUE
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE_GETIG(GBAS,KVEC,EPW,NG,IGVEC)
!     ******************************************************************
!     **  COUNTS NUMBER OF G-VECTORS                                  **
!     **  IF(NGX=0) ONLY THE NUMBER OF G-VECTORS IS CALCULATED        **
!     **  G=KVEC+GBAS*IG                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: GBAS(3,3)
      REAL(8)   ,INTENT(IN)  :: KVEC(3)
      REAL(8)   ,INTENT(IN)  :: EPW
      INTEGER(4),INTENT(IN)  :: NG
      INTEGER(4),INTENT(OUT) :: IGVEC(3,NG)
      REAL(8)                :: RAD
      REAL(8)                :: G1,G2,G3
      REAL(8)                :: T1,T2,T3
      INTEGER(4)             :: MIN1,MIN2,MIN3,MAX1,MAX2,MAX3
      INTEGER(4)             :: I,J,K,IG
      REAL(8)                :: G2A(NG)
      REAL(8)                :: GSQUARE
      INTEGER(4)             :: ICOUNT,IGARR(3),ISTART,ISTOP
      REAL(8)                :: KARTK(3)
      INTEGER(4)             :: FROM,TO
      REAL(8)                :: SVAR
!     ******************************************************************
      KARTK=MATMUL(GBAS,KVEC)
!
!     ==================================================================
!     ==  FIND ALL LATTICE VECTORS WITH 0.5*G**2 < EPW                ==
!     ==================================================================
      RAD=SQRT(2.D0*EPW)+1.D-6
      CALL BOXSPH(GBAS,-KARTK(1),-KARTK(2),-KARTK(3),RAD &
     &           ,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
!     ==================================================================
!     ==  FIND ALL LATTICE VECTORS WITH 0.5*G**2 < EPW                ==
!     ==================================================================
      ICOUNT=0
      DO I=MIN1,MAX1
        T1=REAL(I,KIND=8)
        DO J=MIN2,MAX2
          T2=REAL(J,KIND=8)
          DO K=MIN3,MAX3
            T3=REAL(K,KIND=8)
            G1=KARTK(1)+GBAS(1,1)*T1+GBAS(1,2)*T2+GBAS(1,3)*T3
            G2=KARTK(2)+GBAS(2,1)*T1+GBAS(2,2)*T2+GBAS(2,3)*T3
            G3=KARTK(3)+GBAS(3,1)*T1+GBAS(3,2)*T2+GBAS(3,3)*T3
            GSQUARE=G1**2+G2**2+G3**2
            IF(0.5D0*GSQUARE.GT.EPW) CYCLE
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.NG) THEN
              CALL ERROR$MSG('INCONSISTENT #(GVECTORS)')
              CALL ERROR$STOP('PLANEWAVE$GETIG')
            END IF
            IGVEC(1,ICOUNT)=I
            IGVEC(2,ICOUNT)=J
            IGVEC(3,ICOUNT)=K
            G2A(ICOUNT)=GSQUARE
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  ORDER G-VECTORS ACCORDING TO INCREASING G**2                ==
!     ==================================================================
      CALL SORT$SET(NG,G2A)
      CALL SORT$RESTART
      CALL SORT$FLIP(FROM,TO)
      DO WHILE (FROM.NE.0.OR.TO.NE.0)
        IF(TO.EQ.0) THEN
          IGARR(:)=IGVEC(:,FROM)
          SVAR=G2A(FROM)
        ELSE IF(FROM.EQ.0) THEN
          IGVEC(:,TO)=IGARR(:)
          G2A(TO)=SVAR
        ELSE
          IGVEC(:,TO)=IGVEC(:,FROM)
          G2A(TO)=G2A(FROM)
        END IF
        CALL SORT$FLIP(FROM,TO)
      ENDDO                
      CALL SORT$UNSET
!
!     =============================================================================
!     ==  ORDER G-VECTORS WITHIN A START OF G-VECTORS OF SAME LENGTH             ==
!     =============================================================================
      ISTART=1
      SVAR=G2A(ISTART)
      DO IG=2,NG
        IF(ABS(G2A(IG)-SVAR).LT.1.D-5.AND.IG.NE.NG) CYCLE
        ISTOP=IG-1
        IF(IG.EQ.NG)ISTOP=NG
        DO I=ISTART,ISTOP
          DO J=I+1,ISTOP
            IF(IGVEC(1,I).GT.IGVEC(1,J)) THEN
              IGARR(:)=IGVEC(:,I)
              IGVEC(:,I)=IGVEC(:,J)
              IGVEC(:,J)=IGARR(:)
            ELSE IF(IGVEC(1,I).EQ.IGVEC(1,J)) THEN
              IF (IGVEC(2,I).GT.IGVEC(2,J)) THEN
                IGARR(:)=IGVEC(:,I)
                IGVEC(:,I)=IGVEC(:,J)
                IGVEC(:,J)=IGARR(:)
              ELSE IF(IGVEC(2,I).EQ.IGVEC(2,J)) THEN
                IF(IGVEC(3,I).GT.IGVEC(3,J)) THEN
                  IGARR(:)=IGVEC(:,I)
                  IGVEC(:,I)=IGVEC(:,J)
                  IGVEC(:,J)=IGARR(:)
                END IF
              END IF
            END IF
          END DO
        END DO
        IF(ISTOP.EQ.NG) EXIT
        ISTART=ISTOP+1
        SVAR=G2A(ISTART)
      ENDDO
      RETURN
      END 
!MEBAREK: PLEASE TEST THIS ALSO WITH NR1.NE.NR2.NE.NR3
!
!     ..................................................................
      SUBROUTINE PLANEWAVE_MAPSTRIPESTOTASKS(NR1,NR2,NR3 &
     &            ,NG,IGVEC,TINV,MINUSG,NTASKS,ITASKOFYZ)
!     ******************************************************************
!     ***  MAP STRIPES (IG2,IG3) ONTO TASKS                          ***
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NR1
      INTEGER(4),INTENT(IN) :: NR2
      INTEGER(4),INTENT(IN) :: NR3
      INTEGER(4),INTENT(IN) :: NG
      LOGICAL(4),INTENT(IN) :: TINV
      INTEGER(4),INTENT(IN) :: IGVEC(3,NG)     
      INTEGER(4),INTENT(IN) :: NTASKS
      INTEGER(4),INTENT(IN) :: MINUSG(NG)         ! IGGLOB -> IXYZ
      INTEGER(4),INTENT(OUT):: ITASKOFYZ(NR2,NR3)  ! IYZ    -> ITASK
      INTEGER(4)            :: IGROUPOFYZ(NR2,NR3)
      INTEGER(4)            :: NGROUP
      INTEGER(4),ALLOCATABLE:: NGOFGROUP(:) !(NGROUP)
      INTEGER(4),ALLOCATABLE:: ITASKOFGROUP(:) !(NGROUP)
      INTEGER(4)            :: IG,I2,I3,IGROUP,IGROUP1,IGROUP2,I
      INTEGER(4)            :: ITASK
      INTEGER(4)            :: CHOICE
      INTEGER(4)            :: ONEARR(1)
      REAL(8)   ,ALLOCATABLE:: COSTOFTASK(:)  !(NGROUP)
      REAL(8)               :: COSTOFG
      REAL(8)               :: COSTOFSTRIPE
!     ******************************************************************
             CALL TRACE$PUSH('PLANEWAVE_MAPSTRIPESTOTASKS')
!
!     ==================================================================
!     ==  DISTRIBUTE STRIPES CHARACTERIZED BY (I2,I3) INTO GROUPS     ==
!     ==  ASSUMING AN INFINITE NUMBER OF TASKS.                       ==
!     ==  EVERY PAIR (+G AND -G) MUST BE IN THE SAME GROUP            ==
!     ==                                                              ==
!     ==  NGROUP IS THE NUMBER OF INDEPENDENT GROUPS POSSIBLE         ==
!     ==  IGROUPOFYZ MAPS EACH STRIPE TO A GROUP                      ==
!     ==================================================================
      IGROUPOFYZ(:,:)=0
      NGROUP=0
      DO IG=1,NG
        I2=MODULO(IGVEC(2,IG),NR2)+1  !<
        I3=MODULO(IGVEC(3,IG),NR3)+1  !<
        IF(IGROUPOFYZ(I2,I3).EQ.0) THEN
          NGROUP=NGROUP+1
          IGROUPOFYZ(I2,I3)=NGROUP
        END IF
        IGROUP1=IGROUPOFYZ(I2,I3)
        IF(TINV) THEN
          I2=MODULO(IGVEC(2,MINUSG(IG)),NR2)+1  !<
          I3=MODULO(IGVEC(3,MINUSG(IG)),NR3)+1  !<
          IGROUP2=IGROUPOFYZ(I2,I3)
          IF(IGROUP1.NE.IGROUP2) THEN
            IF(IGROUP2.EQ.0) THEN
              IGROUPOFYZ(I2,I3)=IGROUP1
            ELSE
              CALL ERROR$MSG('STRIPES RELATED BY INVERSION ENDED UP')
              CALL ERROR$MSG('IN DIFFERENT GROUPS')
              CALL ERROR$STOP('PLANEWAVE_MAPSTRIPESTOTASKS')
            END IF
          END IF
        END IF
      ENDDO         
!
!     ==================================================================
!     ==  DETERMINE THE NUMBER OF G-VECTORS IN EACH GROUP             ==
!     ==================================================================
      ALLOCATE(NGOFGROUP(NGROUP))
      NGOFGROUP=0
      DO IG=1,NG
        I2=MODULO(IGVEC(2,IG),NR2)+1   !<
        I3=MODULO(IGVEC(3,IG),NR3)+1   !<
        IGROUP=IGROUPOFYZ(I2,I3)
        NGOFGROUP(IGROUP)=NGOFGROUP(IGROUP)+1
      ENDDO
!PRINT*,'NGOFGROUP ',NGOFGROUP
!
!     ==================================================================
!     ==  DISTRIBUTE GROUPS ONTO TASKS                                ==
!     ==================================================================
      CHOICE=2
      IF(CHOICE.EQ.1) THEN
        ALLOCATE(ITASKOFGROUP(NGROUP))
        ITASKOFGROUP(:)=0
        DO IGROUP=1,NGROUP
!         == ASSUME INFINITE NUMBER OF TASKS AND ATTRIBUTE WITH DECREASING 
!         == NUMBERS OF G-VECTORS
!         == SELECT THE LARGEST REMAINING GROUP ==========================
          ITASKOFGROUP(IGROUP:IGROUP)=MAXLOC(NGOFGROUP,MASK=ITASKOFGROUP.EQ.0)
          ITASK=ITASKOFGROUP(IGROUP)
          IF(ITASK.EQ.0) THEN
            CALL ERROR$MSG('ERROR')
            CALL ERROR$STOP('PLANEWAVE_MAPSTRIPESTOTASKS')
          END IF
!         == MAP TO FINITE NUMBER OF TASKS ACCORDING 1,2,..N,N,N-1,...1,1,2..
          ITASK=MODULO(ITASK-1,2*NTASKS)+1
          IF(ITASK.GT.NTASKS) ITASK=2*NTASKS-(ITASK-1)
          ITASKOFGROUP(IGROUP)=ITASK
        ENDDO
      ELSE IF(CHOICE.EQ.2) THEN
!       == USE THESE PARAMETERS TO GUIDE LOADBALANCE ===================
        COSTOFG=1.D0   
        COSTOFSTRIPE=1.D-1*REAL(NR1,KIND=8)*COSTOFG
!
        ALLOCATE(COSTOFTASK(NTASKS))
        COSTOFTASK=0
        ALLOCATE(ITASKOFGROUP(NGROUP))
        ITASKOFGROUP(:)=0
        DO I=1,NGROUP
!         == SELECT THE LARGEST REMAINING GROUP ==========================
          ONEARR=MAXLOC(NGOFGROUP,MASK=ITASKOFGROUP.EQ.0)
          IGROUP=ONEARR(1)
!         == SELECT THE TASK WITH THE LEAST G-VECTORS ====================
          ONEARR=MINLOC(COSTOFTASK)
          ITASK=ONEARR(1)
          ITASKOFGROUP(IGROUP)=ITASK
          COSTOFTASK(ITASK)=COSTOFTASK(ITASK) &
     &                     +COSTOFG*REAL(NGOFGROUP(IGROUP),KIND=8) &
     &                     +COSTOFSTRIPE
        ENDDO
        DEALLOCATE(COSTOFTASK)
      END IF
      DEALLOCATE(NGOFGROUP)
!PRINT*,'ITASKOFGROUP ',ITASKOFGROUP
!
!     ==================================================================
!     ==  CREATE MAPPING OF STRIPES TO TASKS                          ==
!     ==================================================================
      DO I2=1,NR2
        DO I3=1,NR3
          IGROUP=IGROUPOFYZ(I2,I3)
          IF(IGROUP.NE.0) THEN
            ITASKOFYZ(I2,I3)=ITASKOFGROUP(IGROUP)
          ELSE
            ITASKOFYZ(I2,I3)=0.D0
          END IF
        ENDDO
      ENDDO
      DEALLOCATE(ITASKOFGROUP)
!PRINT*,'ITASKOFYZ ',ITASKOFYZ
                                            CALL TRACE$POP
      RETURN
      END
!
!     .....................................................GEGENEW......
      SUBROUTINE PLANEWAVE_GVECTORS(NGL,GK)
!     **                                                              **
!     **  GENERATES THE RECIPROCAL LATTICE VECTORS                    **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NGL
      REAL(8)   ,INTENT(OUT):: GK(3,NGL)
      REAL(8)               :: GBAS(3,3)
      REAL(8)               :: KVEC(3)
      INTEGER(4)            :: IG
      REAL(8)               :: T1,T2,T3
      REAL(8)               :: KARTK(3)
      INTEGER(4)            :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
      IF(NGL.NE.THIS%NGLARR(THISTASK)) THEN
        CALL ERROR$MSG('SIZE INCONSISTENT')
        CALL ERROR$STOP('PLANEWAVE_GVECTORS')
      END IF
      GBAS(:,:)=THIS%GBAS(:,:)
      KVEC(:)=THIS%KVEC(:)
      KARTK=MATMUL(GBAS,KVEC)
      DO IG=1,NGL
        T1=REAL(THIS%IGVEC(1,IG),KIND=8)
        T2=REAL(THIS%IGVEC(2,IG),KIND=8)
        T3=REAL(THIS%IGVEC(3,IG),KIND=8)
        GK(1,IG) = KARTK(1)+GBAS(1,1)*T1+GBAS(1,2)*T2+GBAS(1,3)*T3
        GK(2,IG) = KARTK(2)+GBAS(2,1)*T1+GBAS(2,2)*T2+GBAS(2,3)*T3
        GK(3,IG) = KARTK(3)+GBAS(3,1)*T1+GBAS(3,2)*T2+GBAS(3,3)*T3
      ENDDO
      RETURN
      END
!
!     ................................................WAVES_EIKR........
      SUBROUTINE PLANEWAVE_EIKR(NR1,NR2,NR3,NR1START,NR1END,KVEC,EIKR)
!     ******************************************************************
!     **                                                              **
!     ** CALCULATES PHASE  FACTOR EXP( I*K*R ) ON THE REAL SPACE GRID **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NR1,NR2,NR3
      INTEGER(4),INTENT(IN) :: NR1START,NR1END
      REAL(8)   ,INTENT(IN) :: KVEC(3)
      COMPLEX(8),INTENT(OUT):: EIKR((NR1END-NR1START+1)*NR2*NR3)
      COMPLEX(8)            :: DEIK1,DEIK2,DEIK3
      COMPLEX(8)            :: EIKR1,EIKR2,EIKR3
      REAL(8)               :: PI
      REAL(8)               :: PHAS1,PHAS2,PHAS3
      INTEGER(4)            :: IRR,IR1,IR2,IR3
!     ******************************************************************
!     == SHORTCUT FOR K=0 ==============================================
      IF(KVEC(1)**2+KVEC(2)**2+KVEC(3)**2.LT.1.D-8) THEN
        EIKR(:)=(1.D0,0.D0)
        RETURN
      END IF        
!     == NOW OTHER K-POINTS ============================================
      PI=4.D0*ATAN(1.D0)
      PHAS1=2.D0*PI * KVEC(1) / DBLE(NR1)
      PHAS2=2.D0*PI * KVEC(2) / DBLE(NR2)
      PHAS3=2.D0*PI * KVEC(3) / DBLE(NR3)
      DEIK1=CMPLX(COS(PHAS1),SIN(PHAS1),8)
      DEIK2=CMPLX(COS(PHAS2),SIN(PHAS2),8)
      DEIK3=CMPLX(COS(PHAS3),SIN(PHAS3),8)
      EIKR3=(1.D0,0.D0)
      IRR=0
      DO IR3=1,NR3
        EIKR2=EIKR3
        DO IR2=1,NR2
          EIKR1=EIKR2*DEIK1**(NR1START-1)
          DO IR1=NR1START,NR1END
            IRR=IRR+1
            EIKR(IRR)=EIKR1    ! HERE EIKR IS EVALUATED
            EIKR1=EIKR1*DEIK1
          ENDDO
          EIKR2=EIKR2*DEIK2
        ENDDO
        EIKR3=EIKR3*DEIK3
      ENDDO
      RETURN
      END
!
!     .....................................................PHFAC........
      SUBROUTINE PLANEWAVE$STRUCTUREFACTOR(R,NGL,EIGR)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE STRUCTURE FACTORS EXP(-I*G*R)                     **
!     **                                                              **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: R(3)
      INTEGER(4),INTENT(IN)  :: NGL
      COMPLEX(8),INTENT(OUT) :: EIGR(NGL)
      INTEGER(4)             :: NGH
      INTEGER(4)             :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
      IF(.NOT.THIS%TSUPER) THEN
        IF(NGL.NE.THIS%NGLARR(THISTASK)) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$STOP('PLANEWAVE$STRUCTUREFACTOR')
        END IF
        CALL PLANEWAVE_STRUCTUREFACTOR(R,THIS%GBAS,THIS%KVEC &
     &                         ,THIS%NR1,THIS%NR2,THIS%NR3,NGL &
     &                         ,THIS%IGVEC,EIGR)
      ELSE
        NGH=THIS%NGHARR(THISTASK)
        IF(NGL.NE.NGH) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$STOP('PLANEWAVE$STRUCTUREFACTOR')
        END IF
        CALL PLANEWAVE_STRUCTUREFACTOR(R,THIS%GBAS,THIS%KVEC &
     &                         ,THIS%NR1,THIS%NR2,THIS%NR3,NGL &
     &                         ,THIS%IGVECH,EIGR)
      END IF
      RETURN
      END
!
!     .....................................................PHFAC........
      SUBROUTINE PLANEWAVE_STRUCTUREFACTOR(R,GBAS,KVEC,NR1,NR2,NR3,NG &
     &           ,IGVEC,EIGR)
!     **                                                              **
!     **  EVALUATES  EXP(-I*(G+K)*R)                                  **
!     **                                                              **
!     **  INPUT:                                                      **
!     **    R       ATOMIC POSITION                                   **
!     **    XK      POSITION OF THE K-VECTOR                          **
!     **    GBAS    RECIPROCAL LATTICE VECTORS                        **
!     **    NR1,NR2,NR3                                               **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R(3)      ! ATOMIC POSITION
      REAL(8)   ,INTENT(IN) :: KVEC(3)   ! K-VECTOR
      REAL(8)   ,INTENT(IN) :: GBAS(3,3) ! REC. LATTICE VECTORS
      INTEGER(4),INTENT(IN) :: NR1       ! #(GRID-POINTS IN R-SPACE || FIRST LATTICE VECTOR)
      INTEGER(4),INTENT(IN) :: NR2
      INTEGER(4),INTENT(IN) :: NR3
      INTEGER(4),INTENT(IN) :: NG        ! #(LOCAL PLANE WAVES)
      INTEGER(4),INTENT(IN) :: IGVEC(3,NG)
      COMPLEX(8),INTENT(OUT):: EIGR(NG)  ! E**(-I*G*R)
      INTEGER(4)            :: NRMAX
      COMPLEX(8),POINTER,SAVE :: DCWORK(:,:)
      COMPLEX(8)            :: CFACP,CSVARP,CFACM,CSVARM
      REAL(8)               :: GR1,GR2,GR3
      REAL(8)               :: RK1,RK2,RK3
      INTEGER(4)            :: I,J,K,IG 
      REAL(8)               :: KARTK(3)
      INTEGER(4),SAVE       :: NRMAXSAVE=-1
!     ******************************************************************
      NRMAX=MAX(NR1,NR2,NR3)
      IF(NRMAXSAVE.LT.NRMAX) THEN
        IF(NRMAXSAVE.NE.-1)DEALLOCATE(DCWORK)
        ALLOCATE(DCWORK(-NRMAX:NRMAX,3))
        NRMAXSAVE=NRMAX
      END IF
      KARTK=MATMUL(GBAS,KVEC)
!
!     ==================================================================
!     ==  FIRST                                                       ==
!     ==================================================================
      GR1=GBAS(1,1)*R(1)+GBAS(2,1)*R(2)+GBAS(3,1)*R(3)
      RK1=R(1)*KARTK(1)
      CFACP =CMPLX(COS(GR1),-SIN(GR1),8)
      CSVARP=CMPLX(COS(RK1),-SIN(RK1),8)
      CFACM=CONJG(CFACP)
      CSVARM=CSVARP
      DO I=1,NR1
        DCWORK(I-1,1)=CSVARP
        CSVARP=CSVARP*CFACP
        DCWORK(-I+1,1)=CSVARM
        CSVARM=CSVARM*CFACM
      ENDDO
!
!     ==================================================================
!     ==  SECOND DIRECTION                                            ==
!     ==================================================================
      GR2=GBAS(1,2)*R(1)+GBAS(2,2)*R(2)+GBAS(3,2)*R(3)
      RK2=R(2)*KARTK(2)
      CFACP =CMPLX(COS(GR2),-SIN(GR2),8)
      CSVARP=CMPLX(COS(RK2),-SIN(RK2),8)
      CFACM =CONJG(CFACP)
      CSVARM=CSVARP
      DO I=1,NR2
        DCWORK(I-1,2)=CSVARP
        CSVARP=CSVARP*CFACP
        DCWORK(-I+1,2)=CSVARM
        CSVARM=CSVARM*CFACM
      ENDDO
!
!     ==================================================================
!     ==  THIRD DIRECTION                                             ==
!     ==================================================================
      GR3=GBAS(1,3)*R(1)+GBAS(2,3)*R(2)+GBAS(3,3)*R(3)
      RK3=R(3)*KARTK(3)
      CFACP =CMPLX(COS(GR3),-SIN(GR3),8)
      CSVARP=CMPLX(COS(RK3),-SIN(RK3),8)
      CFACM =CONJG(CFACP)
      CSVARM=CSVARP
      DO I=1,NR3
        DCWORK(I-1,3)=CSVARP
        CSVARP=CSVARP*CFACP
        DCWORK(-I+1,3)=CSVARM
        CSVARM=CSVARM*CFACM
      ENDDO
!
!     ==================================================================
!     ==  COMBINE THREE DIRECTIONS TO STRUCTURE FACTOR                ==
!     ==================================================================
      DO IG=1,NG
        I=IGVEC(1,IG) 
        J=IGVEC(2,IG) 
        K=IGVEC(3,IG) 
        EIGR(IG)=DCWORK(I,1)*DCWORK(J,2)*DCWORK(K,3)
      ENDDO
      RETURN
      END
!
!     .................................................................
      SUBROUTINE PLANEWAVE$UNSUPER(NGL,FOFG,FOFGUNSUPER)
!     ******************************************************************
!     **  PLANEWAVE$UNSUPER                                           **
!     **  RECONSTRUCTS THE REGULAR WAVE FUNCTIONS FROM THE            **
!     **  SUPER WAVE FUNCTIONS                                        **
!     **  FOFG(1)=PSI+   FOFGUNSUPER(1)=PSI(1)                        **
!     **                 FOFGUNSUPER(2)=PSI(2)                        **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NGL
      COMPLEX(8),INTENT(IN) :: FOFG(NGL)
      COMPLEX(8),INTENT(OUT):: FOFGUNSUPER(NGL,2)
      INTEGER(4)            :: IG
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)            :: CSVAR1,CSVAR2
      INTEGER(4)            :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
      IF(.NOT.THIS%TINV) THEN
        CALL ERROR$MSG('NO TIME INVERSION SYMMETRY')
        CALL ERROR$STOP('PLANEWAVE$UNSUPER')
      END IF
      IF(THIS%TSUPER) THEN
        CALL ERROR$MSG('INVERSION DOES NOT MAKE SENSE')
        CALL ERROR$MSG('IF SUPER WAVE FUNCTION INTERFACE IS SELECTED')
        CALL ERROR$STOP('PLANEWAVE$UNSUPER')
      END IF
      IF(NGL.NE.THIS%NGLARR(THISTASK)) THEN
        CALL ERROR$MSG('INCONSISTENT SIZE')
        CALL ERROR$STOP('PLANEWAVE$UNSUPER')
      END IF
      DO IG=1,NGL
        CSVAR1=FOFG(IG)
        CSVAR2=CONJG(FOFG(THIS%MINUSG(IG)))
        FOFGUNSUPER(IG,1)=0.5D0*(CSVAR1+CSVAR2)
        FOFGUNSUPER(IG,2)=-0.5D0*CI*(CSVAR1-CSVAR2)
      ENDDO
      RETURN
      END
!
!     .................................................................
      SUBROUTINE PLANEWAVE$INVERTG(NGL,FOFG,FOFMG)
!     ******************************************************************
!     **  PLANEWAVE$INVERTG                                           **
!     **  RETURNS THE COMPLEX CONJUGATE OF F(-G)                      **
!     **  IF FOFG=PSI+   THEN FOFMG=PSI-                              **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NGL
      COMPLEX(8),INTENT(IN) :: FOFG(NGL)
      COMPLEX(8),INTENT(OUT):: FOFMG(NGL)
      INTEGER(4)            :: IG
      INTEGER(4)            :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
      IF(.NOT.THIS%TINV) THEN
        CALL ERROR$MSG('NO TIME INVERSION SYMMETRY')
        CALL ERROR$STOP('PLANEWAVE$INVERTG')
      END IF
      IF(THIS%TSUPER) THEN
        CALL ERROR$MSG('INVERSION DOES NOT MAKE SENSE')
        CALL ERROR$MSG('IF SUPER WAVE FUNCTION INTERFACE IS SELECTED')
        CALL ERROR$STOP('PLANEWAVE$INVERTG')
      END IF
      IF(NGL.NE.THIS%NGLARR(THISTASK)) THEN
        CALL ERROR$MSG('INCONSISTENT SIZE')
        CALL ERROR$STOP('PLANEWAVE$INVERTG')
      END IF
      DO IG=1,NGL
        FOFMG(IG)=CONJG(FOFG(THIS%MINUSG(IG)))
      ENDDO
      RETURN
      END
!
!     .................................................................
      SUBROUTINE PLANEWAVE$SUPFFT(ID,NFFT,NGLH,FOFG,NRL,FOFR)
!     ******************************************************************
!     **                                                              **
!     **  SUPER WAVE FUNCTION KEEP ONLY ON G-VECTOR FOR EACH PAIR     **
!     **  OF INVERSION SYMMETRIC G-VECTORS IN G-SPACE, AND ONLY       ** 
!     **  THE REAL PART OF THE WAVE FUNCTION IN R-SPACE               **
!     **                                                              **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4)  ,INTENT(IN)    :: NFFT
      INTEGER(4)  ,INTENT(IN)    :: NGLH
      INTEGER(4)  ,INTENT(IN)    :: NRL
      COMPLEX(8)  ,INTENT(INOUT) :: FOFG(NGLH,NFFT)
      REAL(8)     ,INTENT(INOUT) :: FOFR(NRL,NFFT)
      INTEGER(4)                 :: NGL
      INTEGER(4)                 :: NR1,NR2,NR3
      INTEGER(4)                 :: NSTRIPELX
      INTEGER(4)                 :: IGH,IG1,IG2,IFFT,IR,IG
      COMPLEX(8),PARAMETER       :: CI=(0.D0,1.D0)
      COMPLEX(8),ALLOCATABLE     :: FOFG1(:)
      COMPLEX(8),ALLOCATABLE     :: FOFR1(:)
      COMPLEX(8)                 :: F1,F2
      INTEGER(4)                 :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
                                 CALL TRACE$PUSH('PLANEWAVE$SUPFFT')
                                 CALL TIMING$CLOCKON('PLANEWAVE$SUPFFT')
!
!     ==================================================================
!     == SET UP INITIAL DATA                                          ==
!     ==================================================================
      NR1=THIS%NR1
      NR2=THIS%NR2
      NR3=THIS%NR3
      NGL=THIS%NGLARR(THISTASK)
      NSTRIPELX=MAXVAL(THIS%NSTRIPELARR)
!
!     ==================================================================
!     == CONSISTENCY CHECKS                                           ==
!     ==================================================================
      IF(.NOT.THIS%TSUPER) THEN
        CALL ERROR$MSG('SUPERWAVE FUNCTIONS ONLY WORK AT THE GAMMA POINT')
        CALL ERROR$STOP('PLANEWAVE$SUPFFT')
      END IF
      IF(NRL.NE.THIS%NR1L*THIS%NR2*THIS%NR3) THEN
        CALL ERROR$MSG('SIZE INCONSISTENT')
        CALL ERROR$I4VAL('NRL',NRL)
        CALL ERROR$I4VAL('THIS%(NR1L*NR2*%NR3)',THIS%NR1L*THIS%NR2*THIS%NR3)
        CALL ERROR$I4VAL('THIS%NR1L',THIS%NR1L)
        CALL ERROR$I4VAL('THIS%NR2',THIS%NR2)
        CALL ERROR$I4VAL('THIS%NR3',THIS%NR3)
        CALL ERROR$STOP('PLANEWAVE$SUPFFT')
      END IF
      IF(.NOT.THIS%TSUPER) THEN
        CALL ERROR$MSG('SUPER WAVE FUNCTION ONLY ACCESSIBLE FOR GAMMA POINT')
        CALL ERROR$STOP('PLANEWAVE$SUPFFT')
      END IF
      IGH=0
      DO IG=1,NGL
        IF(THIS%MINUSG(IG).GE.IG) IGH=IGH+1
      ENDDO
      IF(IGH.NE.NGLH) THEN
        CALL ERROR$MSG('SIZE INCONSISTENT')
        CALL ERROR$I4VAL('IGH',IGH)
        CALL ERROR$I4VAL('NGLH',NGLH)
        CALL ERROR$STOP('PLANEWAVE$SUPFFT')
      END IF
!
!     ==================================================================
!     == FOURIER TRANSFORM VIA SUPER WAVE FUNCTIONS                   ==
!     ==================================================================
      ALLOCATE(FOFG1(NGL))
      ALLOCATE(FOFR1(NRL))
      IF(ID.EQ.'GTOR') THEN
        DO IFFT=1,NFFT,2
!
!         ==============================================================
!         == SET UP SUPER WAVE FUNCTIONS IN G-SPACE                   ==
!         ==============================================================
          IGH=0
          DO IG1=1,NGL
            IG2=THIS%MINUSG(IG1)
            IF(IG2.LT.IG1) CYCLE
            IGH=IGH+1
            IF(IFFT+1.LE.NFFT) THEN
              FOFG1(IG1)=      FOFG(IGH,IFFT) +CI*      FOFG(IGH,IFFT+1)
              FOFG1(IG2)=CONJG(FOFG(IGH,IFFT))+CI*CONJG(FOFG(IGH,IFFT+1))
            ELSE
              FOFG1(IG1)=FOFG(IGH,IFFT)
              FOFG1(IG2)=CONJG(FOFG(IGH,IFFT))
            END IF
          ENDDO
!
!         ==============================================================
!         == NOW TRANSFORM SUPER WAVE FUNCTIONS TO R-SPACE            ==
!         ==============================================================
          CALL PLANEWAVE_FFTGTOR(THIS%CID,NTASKS,NGL,NRL,NR1,NR2,NR3 &
     &           ,NSTRIPELX,THIS%NSTRIPELARR,THIS%NR1LARR &
     &           ,THIS%IGTOSTRIPE,THIS%ISTRIPETOYZ &
     &           ,FOFG1,FOFR1)
!
!         ==============================================================
!         == UNRAVEL SUPER WAVE FUNCTIONS IN R-SPACE                  ==
!         ==============================================================
          IF(IFFT+1.LE.NFFT) THEN
            DO IR=1,NRL
              FOFR(IR,IFFT)  =REAL(FOFR1(IR),KIND=8)
              FOFR(IR,IFFT+1)=AIMAG(FOFR1(IR))
            ENDDO
          ELSE
            DO IR=1,NRL
              FOFR(IR,IFFT)=REAL(FOFR1(IR),KIND=8)
            ENDDO
          END IF
        ENDDO
      ELSE IF(ID.EQ.'RTOG') THEN
        DO IFFT=1,NFFT,2
!
!         ==============================================================
!         == SET UP SUPER WAVE FUNCTIONS IN R-SPACE                   ==
!         ==============================================================
          IF(IFFT+1.LE.NFFT) THEN
            DO IR=1,NRL
              FOFR1(IR)=CMPLX(FOFR(IR,IFFT),FOFR(IR,IFFT+1),8)
            ENDDO
          ELSE
            DO IR=1,NRL
              FOFR1(IR)=CMPLX(FOFR(IR,IFFT),0.D0,8)
            ENDDO
          END IF
! 
!         ==============================================================
!         == NOW TRANSFORM SUPER WAVE FUNCTIONS TO G-SPACE            ==
!         ==============================================================
          CALL PLANEWAVE_FFTRTOG(THIS%CID,NTASKS,NGL,NRL,NR1,NR2,NR3 &
     &           ,NSTRIPELX,THIS%NSTRIPELARR,THIS%NR1LARR &
     &           ,THIS%IGTOSTRIPE,THIS%ISTRIPETOYZ &
     &           ,FOFG1,FOFR1)
!
!         ==============================================================
!         == UNRAVEL SUPER WAVE FUNCTIONS IN R-SPACE                  ==
!         ==============================================================
          IGH=0
          DO IG1=1,NGL
            IG2=THIS%MINUSG(IG1)
            IF(IG2.LT.IG1) CYCLE
            IGH=IGH+1
            F1=FOFG1(IG1)
            F2=CONJG(FOFG1(IG2))
            IF(IFFT+1.LE.NFFT) THEN
              FOFG(IGH,IFFT)  =0.5D0*(F1+F2)
              FOFG(IGH,IFFT+1)=0.5D0*CI*(F2-F1)
            ELSE
              FOFG(IGH,IFFT)=0.5D0*(F1+F2)
            END IF
          ENDDO
          
        ENDDO
      ELSE
        CALL ERROR$MSG('ID MUST BE WITHER GTOR OR RTOG')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PLANEWAVE$SUPFFT')
      END IF
      DEALLOCATE(FOFG1)
      DEALLOCATE(FOFR1)
                                 CALL TIMING$CLOCKOFF('PLANEWAVE$SUPFFT')
                                 CALL TRACE$POP
      RETURN
      END      
!
!     .................................................................
      SUBROUTINE PLANEWAVE$FFT(ID,NFFT,NGL,FOFG,NRL,FOFR)
!     ******************************************************************
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID 
      INTEGER(4)  ,INTENT(IN)    :: NFFT
      INTEGER(4)  ,INTENT(IN)    :: NGL
      INTEGER(4)  ,INTENT(IN)    :: NRL
      COMPLEX(8)  ,INTENT(INOUT) :: FOFG(NGL,NFFT)
      COMPLEX(8)  ,INTENT(INOUT) :: FOFR(NRL,NFFT)
      INTEGER(4)                 :: NR1,NR2,NR3
      INTEGER(4)                 :: NSTRIPELX
      INTEGER(4)                 :: IFFT
      INTEGER(4)                 :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
                                 CALL TIMING$CLOCKON('PLANEWAVE$FFT')
      IF(NGL.NE.THIS%NGLARR(THISTASK)) THEN
        CALL ERROR$MSG('SIZE OF F(G) INCONSISTENT')
        CALL ERROR$I4VAL('NGL',NGL)
        CALL ERROR$I4VAL('THIS%NGLARR(THISTASK)',THIS%NGLARR(THISTASK))
        CALL ERROR$STOP('PLANEWAVE$FFT')
      END IF
      IF(NRL.NE.THIS%NR1L*THIS%NR2*THIS%NR3) THEN
        CALL ERROR$MSG('SIZE OF F(R) INCONSISTENT')
        CALL ERROR$I4VAL('NRL',NRL)
        CALL ERROR$I4VAL('THIS%NR1L',THIS%NR1L)
        CALL ERROR$I4VAL('THIS%NR2',THIS%NR2)
        CALL ERROR$I4VAL('THIS%NR3',THIS%NR3)
        CALL ERROR$STOP('PLANEWAVE$FFT')
      END IF
      IF(THIS%TSUPER) THEN
        CALL ERROR$MSG('WRONG INTERFACE FOR SUPER WAVE FUNCTIONS')
        CALL ERROR$STOP('PLANEWAVE$FFT')
      END IF
      NR1=THIS%NR1
      NR2=THIS%NR2
      NR3=THIS%NR3
      NSTRIPELX=MAXVAL(THIS%NSTRIPELARR)
      IF(ID.EQ.'GTOR') THEN
        DO IFFT=1,NFFT
          CALL PLANEWAVE_FFTGTOR(THIS%CID,NTASKS,NGL,NRL,NR1,NR2,NR3 &
     &           ,NSTRIPELX,THIS%NSTRIPELARR,THIS%NR1LARR &
     &           ,THIS%IGTOSTRIPE,THIS%ISTRIPETOYZ &
     &           ,FOFG(1,IFFT),FOFR(1,IFFT))
        ENDDO
      ELSE IF(ID.EQ.'RTOG') THEN
        DO IFFT=1,NFFT
          CALL PLANEWAVE_FFTRTOG(THIS%CID,NTASKS,NGL,NRL,NR1,NR2,NR3 &
     &           ,NSTRIPELX,THIS%NSTRIPELARR,THIS%NR1LARR &
     &           ,THIS%IGTOSTRIPE,THIS%ISTRIPETOYZ &
     &           ,FOFG(1,IFFT),FOFR(1,IFFT))
        ENDDO
      ELSE
        CALL ERROR$MSG('ID MUST BE WITHER GTOR OR RTOG')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PLANEWAVE$FFT')
      END IF
                                 CALL TIMING$CLOCKOFF('PLANEWAVE$FFT')
      RETURN
      END      
!
!     ..................................................................
      SUBROUTINE PLANEWAVE_FFTGTOR(CID,NTASKS_,NGL,NRL,NR1,NR2,NR3 &
     &           ,NSTRIPELX,NSTRIPELARR,NR1LARR,IGTOSTRIPE,ISTRIPETOYZ &
     &           ,FOFG,FOFR)
!     ******************************************************************
!     ***                        PLANEWAVE_FFTGTOR                    **
!     **                                                              **
!     **  PARALLEL FOURIER TRANSFORM FROM RECIPROCAL TO REAL SPACE    **
!     **                                                              **
!     ******************************************************************
      USE MPE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)   :: CID   ! COMMUNICATOR ID (SEE MPELIB)
      INTEGER(4)  ,INTENT(IN)   :: NTASKS_
      INTEGER(4)  ,INTENT(IN)   :: NGL
      INTEGER(4)  ,INTENT(IN)   :: NRL
      INTEGER(4)  ,INTENT(IN)   :: NR1,NR2,NR3    !GLOBAL
      INTEGER(4)  ,INTENT(IN)   :: IGTOSTRIPE(NGL)
      INTEGER(4)  ,INTENT(IN)   :: NSTRIPELARR(NTASKS_)
      INTEGER(4)  ,INTENT(IN)   :: NR1LARR(NTASKS_)
      INTEGER(4)  ,INTENT(IN)   :: NSTRIPELX
      INTEGER(4)  ,INTENT(IN)   :: ISTRIPETOYZ(NSTRIPELX,NTASKS_)
      COMPLEX(8)  ,INTENT(IN)   :: FOFG(NGL)
      COMPLEX(8)  ,INTENT(OUT)  :: FOFR(NRL)
      COMPLEX(8)  ,ALLOCATABLE  :: UIU(:,:,:)  !(NR1LX,NSTRIPEX,NTASKS)
      COMPLEX(8)  ,ALLOCATABLE  :: FOFR1(:)
      COMPLEX(8)  ,ALLOCATABLE  :: FOFG1(:)
      INTEGER(4)                :: NSTRIPEL
      INTEGER(4)                :: NR1LX
      INTEGER(4)                :: NR1L
      INTEGER(4)                :: ISTRIPEL,ITASK
      INTEGER(4)                :: IR1,IR2,IR3,IR1L
      INTEGER(4)                :: IMID,I23
      INTEGER(4)                :: NR3END1,NR3END2
      INTEGER(4)                :: IND1,IND2,STRIDE2
      INTEGER(4)                :: IG
      INTEGER(4)                 :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(CID,NTASKS,THISTASK)
      IF(NTASKS.NE.NTASKS_) THEN
        CALL ERROR$MSG('INCONSISTENT NTASKS')
        CALL ERROR$STOP('PLANEWAVE_FFTRTOG')
      END IF
      IF(MAX(NR1,NR2,NR3).GT.2048) THEN
        CALL ERROR$MSG('INSUFFICIENT SIZE FOR AUXILIARY ARRAYS') 
        CALL ERROR$STOP('PLANEWAVE_FFTGTOR')
      END IF
      NSTRIPEL=NSTRIPELARR(THISTASK)
      NR1LX=MAXVAL(NR1LARR)
      NR1L=NR1LARR(THISTASK)
!
!     ==================================================================
!     == MAP TO STRIPES                                               ==
!     ==================================================================
      ALLOCATE(FOFG1(NR1*NSTRIPEL))
      FOFG1(:)=(0.D0,0.D0)
      DO IG=1,NGL
        FOFG1(IGTOSTRIPE(IG))=FOFG(IG)
      ENDDO
!
!     ==================================================================
!     == TRANSFORM 1 ST DIMENSION (NR1)    FOFG(NR1,NSTRIPE)          ==
!     ==================================================================
      CALL LIB$FFTC8('GTOR',NR1,NSTRIPEL,FOFG1,FOFG1)
!
!     ==================================================================
!     == COMMUNICATE AMONG TASKS                                      ==
!     ==================================================================
      ALLOCATE(UIU(NR1LX,NSTRIPELX,NTASKS))
      UIU(:,:,:)=CMPLX(0.D0,0.D0,8)
!
      DO ISTRIPEL=1,NSTRIPEL
        IR1=NR1*(ISTRIPEL-1)
        DO ITASK=1,NTASKS
          DO IR1L=1,NR1LARR(ITASK)
            IR1=IR1+1
            UIU(IR1L,ISTRIPEL,ITASK)=FOFG1(IR1)
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(FOFG1)
!
      CALL MPE$TRANSPOSE(CID,UIU)
!
      FOFR(:)=CMPLX(0.D0,0.D0,8)
      IMID=NR3/2
      NR3END1=1
      NR3END2=NR3
      STRIDE2=NR2*NR3
      DO ITASK=1,NTASKS
        DO ISTRIPEL=1,NSTRIPELARR(ITASK)
          I23=ISTRIPETOYZ(ISTRIPEL,ITASK)   ! I23=IR2+NR2*(IR3-1)
          IND1=I23-STRIDE2
          DO IR1L=1,NR1L
!           == FOFR(NR2,NR3,NR1)==========================================
            IND1=IND1+STRIDE2
            FOFR(IND1)=UIU(IR1L,ISTRIPEL,ITASK)
!           FOFR(I23+NR2*NR3*(IR1L-1))=UIU(IR1L,ISTRIPEL,ITASK)
          ENDDO
          IR3=1+INT((I23-1)/NR2)
          IF(IR3.LE.IMID) THEN        ! OBTAIN I23END1/2 TO OBTAIN NR3END1/2
            NR3END1=MAX(NR3END1,IR3)
          ELSE
            NR3END2=MIN(NR3END2,IR3)
          END IF          
        ENDDO
      ENDDO
      DEALLOCATE(UIU)
!
!     ==================================================================
!     == FFT 2ND DIMENSION (NR2):      FOFR(NR2,NR3,NR1L)             ==
!     ==================================================================
      DO IR1L=1,NR1L
        CALL LIB$FFTC8('GTOR',NR2,NR3END1,FOFR(1+(IR1L-1)*NR2*NR3) &
     &                                   ,FOFR(1+(IR1L-1)*NR2*NR3))
      ENDDO
      DO IR1L=1,NR1L
        CALL LIB$FFTC8('GTOR',NR2,NR3-NR3END2+1 &
     &                ,FOFR(1+NR2*(NR3END2-1)+(IR1L-1)*NR2*NR3) &
     &                ,FOFR(1+NR2*(NR3END2-1)+(IR1L-1)*NR2*NR3))
      ENDDO
!
!     ==================================================================
!     == TRANSPOSE TO OPTIMIZE STRIDE                                 ==
!     == FOFR1(IR3,IR2,IR1L)=FOFR(IR2,IR3,IR1L)                       ==
!     ==================================================================
      ALLOCATE(FOFR1(NR3*NR2*NR1L))
      DO IR1L=1,NR1L
        DO IR2=1,NR2
         IND1=IR2+NR2*(-1+NR3*(IR1L-1))
         IND2=NR3*(IR2-1+NR2*(IR1L-1))
         DO IR3=1,NR3
!           IND1=IR2+NR2*(IR3-1+NR3*(IR1L-1))
!           IND2=IR3+NR3*(IR2-1+NR2*(IR1L-1))
            IND2=IND2+1
            IND1=IND1+NR2
            FOFR1(IND2)=FOFR(IND1)
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     == FFT 3RD DIMENSION (NR3)                                      ==
!     ==================================================================
      CALL LIB$FFTC8('GTOR',NR3,NR2*NR1L,FOFR1,FOFR1)
!
!     ==================================================================
!     == TRANSPOSE TO OPTIMIZE STRIDE                                 ==
!     == FOFR(IR1,IR2,IR3)=FOFR(IR2,IR3,IR1)                         ==
!     ==================================================================
      STRIDE2=NR1L*NR2
      DO IR1L=1,NR1L
        DO IR2=1,NR2
          IND1=NR3 *(IR2-1+NR2*(IR1L-1))
          IND2=IR1L+NR1L*(IR2-1-NR2)
          DO IR3=1,NR3
!           IND1=IR3 +NR3 *(IR2-1+NR2*(IR1L-1))
!           IND2=IR1L+NR1L*(IR2-1+NR2*(IR3-1))
            IND1=IND1+1
            IND2=IND2+STRIDE2
            FOFR(IND2)=FOFR1(IND1)
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(FOFR1)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE_FFTRTOG(CID,NTASKS_,NGL,NRL,NR1,NR2,NR3 &
     &          ,NSTRIPELX,NSTRIPELARR,NR1LARR,IGTOSTRIPE,ISTRIPETOYZ &
     &          ,FOFG,FOFR)
!     ******************************************************************
!     **                     PLANEWAVE_FFTRTOG                        **
!     **                                                              **
!     **  PARALLEL FOURIER TRANSFORM FROM RECIPROCAL TO REAL SPACE    **
!     **                                                              **
!     ******************************************************************
      USE MPE_MODULE
      IMPLICIT NONE 
      CHARACTER(*),INTENT(IN)   :: CID
      INTEGER(4)  ,INTENT(IN)   :: NTASKS_
      INTEGER(4)  ,INTENT(IN)   :: NGL
      INTEGER(4)  ,INTENT(IN)   :: NRL
      INTEGER(4)  ,INTENT(IN)   :: NR1,NR2,NR3    !GLOBAL
      INTEGER(4)  ,INTENT(IN)   :: NSTRIPELARR(NTASKS_)
      INTEGER(4)  ,INTENT(IN)   :: NR1LARR(NTASKS_)
      INTEGER(4)  ,INTENT(IN)   :: IGTOSTRIPE(NGL)
      INTEGER(4)  ,INTENT(IN)   :: NSTRIPELX
      INTEGER(4)  ,INTENT(IN)   :: ISTRIPETOYZ(NSTRIPELX,NTASKS_)
      COMPLEX(8)  ,INTENT(OUT)  :: FOFG(NGL)
      COMPLEX(8)  ,INTENT(IN)   :: FOFR(NRL)
      COMPLEX(8)  ,ALLOCATABLE  :: UIU(:,:,:)  !(NR1LX,NSTRIPELX,NTASKS)
      COMPLEX(8)  ,ALLOCATABLE  :: FOFR1(:)
      COMPLEX(8)  ,ALLOCATABLE  :: FOFR2(:)
      COMPLEX(8)  ,ALLOCATABLE  :: FOFG1(:)
      INTEGER(4)                :: NR1L
      INTEGER(4)                :: NR1LX
      INTEGER(4)                :: ISTRIPEL,ITASK
      INTEGER(4)                :: IR1,IR2,IR3,IR1L
      INTEGER(4)                :: IMID,I23
      INTEGER(4)                :: NR3END1,NR3END2
      INTEGER(4)                :: IND1,IND2
      INTEGER(4)                :: NSTRIPEL
      INTEGER(4)                :: IG
      INTEGER(4)                 :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(CID,NTASKS,THISTASK)
      IF(NTASKS.NE.NTASKS_) THEN
        CALL ERROR$MSG('INCONSISTENT NTASKS')
        CALL ERROR$STOP('PLANEWAVE_FFTRTOG')
      END IF
      NR1L=NR1LARR(THISTASK)
      NR1LX=MAXVAL(NR1LARR)
      NSTRIPEL=NSTRIPELARR(THISTASK)
!
!     ==================================================================
!     == TRANSPOSE TO OPTIMIZE STRIDE                                 ==
!     ==================================================================
      ALLOCATE(FOFR1(NR3*NR2*NR1L))
      DO IR1L=1,NR1L
        DO IR2=1,NR2
          DO IR3=1,NR3
            IND1=IR1L+NR1L*(IR2-1+NR2*(IR3-1))
            IND2=IR3 +NR3 *(IR2-1+NR2*(IR1L-1))
            FOFR1(IND2)=FOFR(IND1)
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     == TRANSFORM 3 DIMENSION (NR3)                                  ==
!     ==================================================================
      CALL LIB$FFTC8('RTOG',NR3,NR2*NR1L,FOFR1,FOFR1)
!
!     ==================================================================
!     == EVALUATE NR3END1,NR3END2                                     ==
!     ==================================================================
      IMID=NR3/2
      NR3END1=1
      NR3END2=NR3
      DO ITASK=1,NTASKS
        DO ISTRIPEL=1,NSTRIPELARR(ITASK)
          I23=ISTRIPETOYZ(ISTRIPEL,ITASK)    ! I23=IR2+NR2*(IR3-1)
          IR3=1+INT((I23-1)/NR2)         
          IF(IR3.LE.IMID) THEN       ! OBTAIN I23END1/2 TO OBTAIN NR3END1/2
            NR3END1=MAX(NR3END1,IR3)
          ELSE
            NR3END2=MIN(NR3END2,IR3)
          END IF          
        ENDDO
      ENDDO
      IF(NR3END1.EQ.NR3END2) NR3END2=NR3END2+1
!
!     ==================================================================
!     == TRANSPOSE TO OPTIMIZE STRIDE                                 ==
!     ==================================================================
      ALLOCATE(FOFR2(NR2*NR3*NR1L))
      FOFR2=CMPLX(0.D0,0.D0,8)
      DO IR1L=1,NR1L
        DO IR2=1,NR2
          DO IR3=1,NR3END1
            IND1=IR3+NR3*(IR2-1+NR2*(IR1L-1))
            IND2=IR2+NR2*(IR3-1+NR3*(IR1L-1))
            FOFR2(IND2)=FOFR1(IND1)
          ENDDO
          DO IR3=NR3END2,NR3
            IND1=IR3+NR3*(IR2-1+NR2*(IR1L-1))
            IND2=IR2+NR2*(IR3-1+NR3*(IR1L-1))
            FOFR2(IND2)=FOFR1(IND1)
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(FOFR1)
!
!     ==================================================================
!     == TRANSFORM 2ND DIMENSION (NR2)                                ==
!     ==================================================================
      DO IR1L=1,NR1L
        CALL LIB$FFTC8('RTOG',NR2,NR3END1,FOFR2(1+(IR1L-1)*NR2*NR3) &
     &                                   ,FOFR2(1+(IR1L-1)*NR2*NR3))
      ENDDO
      DO IR1L=1,NR1L
        CALL LIB$FFTC8('RTOG',NR2,NR3-NR3END2+1 &
     &                  ,FOFR2(1+NR2*(NR3END2-1)+(IR1L-1)*NR2*NR3) &
     &                  ,FOFR2(1+NR2*(NR3END2-1)+(IR1L-1)*NR2*NR3))
      ENDDO
!
!     ==================================================================
!     == COMMUNICATE AMONG TASKS                                      ==
!     ==================================================================
      ALLOCATE(UIU(NR1LX,NSTRIPELX,NTASKS))
      UIU(:,:,:)=CMPLX(0.D0,0.D0,8)
      DO ITASK=1,NTASKS
        DO ISTRIPEL=1,NSTRIPELARR(ITASK)
          I23=ISTRIPETOYZ(ISTRIPEL,ITASK)    ! I23=IR2+NR2*(IR3-1)
          DO IR1L=1,NR1L
!           == FOFR(NR2,NR3,NR1)==========================================
            UIU(IR1L,ISTRIPEL,ITASK)=FOFR2(I23+NR2*NR3*(IR1L-1))
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(FOFR2)
      CALL MPE$TRANSPOSE(CID,UIU)
      ALLOCATE(FOFG1(NR1*NSTRIPEL))
      FOFG1(:)=CMPLX(0.D0,0.D0,8)
      DO ISTRIPEL=1,NSTRIPEL
        IR1=NR1*(ISTRIPEL-1)
        DO ITASK=1,NTASKS
          DO IR1L=1,NR1LARR(ITASK)
            IR1=IR1+1
            FOFG1(IR1)=UIU(IR1L,ISTRIPEL,ITASK)
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(UIU)
!
!     ==================================================================
!     == TRANSFORM 1 ST DIMENSION (NR1)                               ==
!     ==================================================================
      CALL LIB$FFTC8('RTOG',NR1,NSTRIPEL,FOFG1,FOFG1)
!
!     ==================================================================
!     == MAP STRIPES TO IG                                            ==
!     ==================================================================
      DO IG=1,NGL
        FOFG(IG)=FOFG1(IGTOSTRIPE(IG))
      ENDDO
      DEALLOCATE(FOFG1)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PLANEWAVE$SCALARPRODUCT(ID,NGL,NDIM,N1,F1,N2,F2,MAT)
!     **************************************************************************
!     **  SCALAR PRODUCT OF SETS OF FUNCTIONS IN RECIPROCAL SPACE             **
!     **                                                                      **
!     **  THE SUBROUTINE IS NOT PARALLELIZED ADD THE RESULT EXTERNALLY        **
!     **  USING MPE$COMBINE('NONE','+',MAT)                                   **
!     **                                                                      **
!     **  THREE OPTIONS SPECIFIED BY ID                                       **
!     **    ID='=' INDICATES THAT F1=F2                                       **
!     **    ID='-' <F1(-)|F2(+)>                                              **
!     **    ID=' '  <F1|F2>                                                   **
!     **                                                                      **
!     **************************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE 
      CHARACTER(*),INTENT(IN)   :: ID    ! ' ','=','-'
      INTEGER(4)  ,INTENT(IN)   :: NDIM
      INTEGER(4)  ,INTENT(IN)   :: NGL
      INTEGER(4)  ,INTENT(IN)   :: N1
      INTEGER(4)  ,INTENT(IN)   :: N2
      COMPLEX(8)  ,INTENT(IN)   :: F1(NGL,NDIM,N1)
      COMPLEX(8)  ,INTENT(IN)   :: F2(NGL,NDIM,N2)
      COMPLEX(8)  ,INTENT(OUT)  :: MAT(N1,N2)
      INTEGER(4)                :: NGAMMA
      INTEGER(4)                :: I1,I2,IDIM
      COMPLEX(8)  ,ALLOCATABLE  :: F2M(:,:)
      LOGICAL(4)  ,PARAMETER    :: TESSL=.TRUE.
!     **************************************************************************
                          CALL TIMING$CLOCKON('PLANEWAVE$SCALARPRODUCT')
      IF(ID.EQ.' '.OR.ID.EQ.'=') THEN
        CALL LIB$SCALARPRODUCTC8((ID.EQ.'='),NGL*NDIM,N1,F1,N2,F2,MAT)
      ELSE IF(ID.EQ.'-') THEN
        IF(.NOT.THIS%TINV.OR.THIS%TSUPER) THEN
          CALL ERROR$MSG('OPTION - CAN ONLY USED INVERSION SYMMETRY')
          CALL ERROR$MSG('AND NOT FOR SUPER WAVE FUNCTIONS')
          CALL ERROR$STOP('PLANEWAVE$SCALARPRODUCT')
        END IF
        ALLOCATE(F2M(NGL,NDIM))
!       == SUM_G F1(-G) F2(G) =CONJG(SUM_G F1(G)^* F2(-G)^*) ============
        DO I2=1,N2
          DO IDIM=1,NDIM
!           __F2M(G)=F2(-G)^*____________________________________________
            CALL PLANEWAVE$INVERTG(NGL,F2(1,IDIM,I2),F2M(1,IDIM)) 
          ENDDO
          CALL LIB$SCALARPRODUCTC8(.FALSE.,NGL*NDIM,N1,F1,1,F2M,MAT(1,I2))
        ENDDO
        DO I1=1,N1
          DO I2=1,N2
            MAT(I1,I2)=CONJG(MAT(I1,I2))
          ENDDO
        ENDDO
        DEALLOCATE(F2M)
      END IF
!
!     ==================================================================
!     ==  COMPLETE FOR SUPERWAVE FUNCTION INTERFACE                   ==
!     ==================================================================
      IF(THIS%TSUPER) THEN
        MAT(:,:)=2.D0*MAT(:,:)
        NGAMMA=THIS%NGAMMA
        IF(NGAMMA.NE.0) THEN
          DO I1=1,N1
            DO I2=1,N2
              DO IDIM=1,NDIM
                MAT(I1,I2)=MAT(I1,I2) &
     &                    -CONJG(F1(NGAMMA,IDIM,I1))*F2(NGAMMA,IDIM,I2)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      END IF
!
!     ==================================================================
!     ==  MULTIPLY WITH GWEIGHT                                       ==
!     ==================================================================
      MAT=MAT*THIS%GWEIGHT
                         CALL TIMING$CLOCKOFF('PLANEWAVE$SCALARPRODUCT')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PLANEWAVE$ADDPRODUCT(ID,NGL,NDIM,N1,F1,N2,F2,MAT)
!     ******************************************************************
!     **                                                              **
!     **  F1=F1+F2*MAT                                                **
!     **  THREE OPTIONS SPECIFIED BY ID                               **
!     **    ID='=' INDICATES THAT F1=F2                               **
!     **    ID='-' <F1(-)|F2(+)>                                      **
!     **    ID=' '  <F1|F2>                                            **
!     **                                                              **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      IMPLICIT NONE 
      CHARACTER(*),INTENT(IN)   :: ID    ! ' ','=','-'
      INTEGER(4)  ,INTENT(IN)   :: NDIM
      INTEGER(4)  ,INTENT(IN)   :: NGL
      INTEGER(4)  ,INTENT(IN)   :: N1
      INTEGER(4)  ,INTENT(IN)   :: N2
      COMPLEX(8)  ,INTENT(INOUT):: F1(NGL,NDIM,N1)
      COMPLEX(8)  ,INTENT(IN)   :: F2(NGL,NDIM,N2)
      COMPLEX(8)  ,INTENT(IN)   :: MAT(N2,N1)
      INTEGER(4)                :: I2,IDIM
      COMPLEX(8)  ,ALLOCATABLE  :: F2M(:,:,:)
      LOGICAL(4)  ,PARAMETER    :: TESSL=.TRUE.
!     ******************************************************************
                             CALL TIMING$CLOCKON('PLANEWAVE$ADDPRODUCT')
      IF(ID.EQ.' ') THEN
        CALL LIB$ADDPRODUCTC8(.FALSE.,NGL*NDIM,N2,N1,F2,MAT,F1)
!
!     ==================================================================
!     == ALLOW F1=F2                                                  ==
!     ==================================================================
      ELSE IF(ID.EQ.'=') THEN
        IF(N1.NE.N2) THEN
          CALL ERROR$MSG('N1 MUST BE EQUAL TO N2 IF ID="=" IS SELECTED')
          CALL ERROR$STOP('PLANEWAVE$ADDPRODUCT')
        END IF
        CALL LIB$ADDPRODUCTC8(.TRUE.,NGL*NDIM,N2,N1,F2,MAT,F1)
!
!     ==================================================================
!     == BUILD SCALARPRODUCT <F1-|F2+>                                ==
!     ==================================================================
      ELSE IF(ID.EQ.'-') THEN
        IF(.NOT.THIS%TINV.OR.THIS%TSUPER) THEN
          CALL ERROR$MSG('OPTION - CAN ONLY USED INVERSION SYMMETRY')
          CALL ERROR$MSG('AND NOT FOR SUPER WAVE FUNCTIONS')
        END IF
        ALLOCATE(F2M(NGL,NDIM,N2))
        DO I2=1,N2
          DO IDIM=1,NDIM
            CALL PLANEWAVE$INVERTG(NGL,F2(1,IDIM,I2),F2M(1,IDIM,I2))
          ENDDO
        ENDDO
        CALL LIB$ADDPRODUCTC8(.FALSE.,NGL*NDIM,N2,N1,F2M,MAT,F1)
        DEALLOCATE(F2M)
      END IF
                            CALL TIMING$CLOCKOFF('PLANEWAVE$ADDPRODUCT')
      RETURN
      END
!
!.......................................................................
#TEMPLATE PLANEWAVE$DISTRIBUTE
(<TYPEID><TYPE>)=([R8][REAL(8)])
                 ([I4][INTEGER(4)])
                 ([C8][COMPLEX(8)])
#BODY
!
!     .................................................................
      SUBROUTINE PLANEWAVE$DISTRIBUTE<TYPEID>(NDIM,NGGLOB,XGLOB,NGLOC,XLOC)
!     ******************************************************************
!     **  PLANEWAVE$DISTRIBUTE                                        **
!     ******************************************************************
      USE MPE_MODULE
      USE PLANEWAVE_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NDIM
      INTEGER(4)  ,INTENT(IN) :: NGGLOB
      INTEGER(4)  ,INTENT(IN) :: NGLOC
      <TYPE>      ,INTENT(IN) :: XGLOB(NDIM,NGGLOB)
      <TYPE>      ,INTENT(OUT):: XLOC(NDIM,NGLOC)
      INTEGER(4)              :: NGSEND
      INTEGER(4)              :: TOTASK
      <TYPE>      ,ALLOCATABLE:: XTMP(:,:)
      INTEGER(4)              :: I
      INTEGER(4)              :: NGLX
      INTEGER(4)              :: NG,NGL
      LOGICAL(4)              :: TSUPER
      INTEGER(4)                 :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
!                              CALL TRACE$PUSH('PLANEWAVE$DISTRIBUTE')
      IF (.NOT.TINI) THEN
        CALL ERROR$MSG('LIST NOT INITIALIZED')
        CALL ERROR$STOP('PLANEWAVE$DISTRIBUTE<TYPEID>')
      ENDIF
      TSUPER=THIS%TSUPER
!     IF (THIS%TSUPER) THEN
!       CALL ERROR$MSG('NOT PERMITTED FOR SUPERWAVEFUNCTION INTERFACE')
!       CALL ERROR$STOP('PLANEWAVE%DISTRIBUTE<TYPEID>')
!     ENDIF
      IF(TSUPER) THEN
        NG=SUM(THIS%NGHARR)
        NGL=THIS%NGHARR(THISTASK)
        NGLX=MAXVAL(THIS%NGHARR(:))
      ELSE
        NG=SUM(THIS%NGLARR)
        NGL=THIS%NGLARR(THISTASK)
        NGLX=MAXVAL(THIS%NGLARR(:))
      END IF
      IF(NGGLOB.NE.NG.AND.THISTASK.EQ.1) THEN
        CALL ERROR$MSG('GLOBAL %(G-VECTORS) INCONSISTENT WITH OBJECT')
        CALL ERROR$I4VAL('NG (ON INPUT)',NGGLOB)
        CALL ERROR$I4VAL('NG (IN OBJECT)',NG)
        CALL ERROR$STOP('PLANEWAVE%DISTRIBUTE<TYPEID>')
      END IF
      IF(NGLOC.NE.NGL) THEN
        CALL ERROR$MSG('LOCAL %(G-VECTORS) INCONSISTENT WITH OBJECT')
        CALL ERROR$STOP('PLANEWAVE%DISTRIBUTE<TYPEID>')
      END IF
!     ******************************************************************
      IF(THISTASK.EQ.1) THEN 
        ALLOCATE(XTMP(NDIM,NGLX))
        DO TOTASK=NTASKS,1,-1
          NGSEND=0
          DO I=1,NG
            IF(THIS%TASKOFIG(I).EQ.TOTASK) THEN
              NGSEND=NGSEND+1
              IF(NGSEND.GT.NGLX) THEN
                CALL ERROR$MSG('PACKAGE LARGER THAN ALLOCATED ARRAY')
                CALL ERROR$STOP('PLANEWAVE_DISTRIBUTE')
              END IF
              XTMP(:,NGSEND)=XGLOB(:,I)
            ENDIF
          ENDDO
          IF(TOTASK.NE.1) THEN
            CALL MPE$SEND(THIS%CID,TOTASK,10,XTMP(:,1:NGSEND))
          ELSE
            XLOC(:,:)=XTMP(:,1:NGSEND)
          END IF
        ENDDO
        DEALLOCATE(XTMP)
      ELSE
        CALL MPE$RECEIVE(THIS%CID,1,10,XLOC)
      ENDIF
!                     CALL TRACE$POP
      RETURN
      END
#END TEMPLATE PLANEWAVE$DISTRIBUTE
!
!.......................................................................
#TEMPLATE PLANEWAVE$COLLECT
(<TYPEID><TYPE>)=([R8][REAL(8)])
                 ([I4][INTEGER(4)])
                 ([C8][COMPLEX(8)])
#BODY
!
!     .................................................................
      SUBROUTINE PLANEWAVE$COLLECT<TYPEID>(NDIM,NGLOC,XLOC,NGGLOB,XGLOB)
!     ******************************************************************
!     **  PLANEWAVE$COLLECT                                           **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NDIM
      INTEGER(4)  ,INTENT(IN) :: NGLOC
      <TYPE>      ,INTENT(IN) :: XLOC(NDIM,NGLOC)
      INTEGER(4)  ,INTENT(IN) :: NGGLOB
      <TYPE>      ,INTENT(OUT):: XGLOB(NDIM,NGGLOB)
      <TYPE>      ,ALLOCATABLE :: X_LOCTMP(:,:)   !(NGLOCX)
      INTEGER(4)               :: NGFROM
      INTEGER(4)               :: FROMTASK
      INTEGER(4)               :: IG,IGL
      INTEGER(4)               :: NGL,NGLX,NG
      LOGICAL(4)               :: TSUPER
      INTEGER(4)                 :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
!                     CALL TRACE$PUSH('PLANEWAVE$COLLECT')
      IF (.NOT.TINI) THEN
        CALL ERROR$MSG('LIST NOT INITIALIZED')
        CALL ERROR$STOP('PLANEWAVE$COLLECT')
      ENDIF
      TSUPER=THIS%TSUPER
!     IF (THIS%TSUPER) THEN
!       CALL ERROR$MSG('NOT PERMITTED FOR SUPERWAVEFUNCTION INTERFACE')
!       CALL ERROR$STOP('PLANEWAVE$COLLECT')
!     ENDIF
      IF(TSUPER) THEN
        NG=SUM(THIS%NGHARR)
        NGL=THIS%NGHARR(THISTASK)
        NGLX=MAXVAL(THIS%NGHARR(:))
      ELSE
        NG=SUM(THIS%NGLARR)
        NGL=THIS%NGLARR(THISTASK)
        NGLX=MAXVAL(THIS%NGLARR(:))
      END IF
      IF(NGGLOB.NE.NG.AND.THISTASK.EQ.1) THEN
        CALL ERROR$MSG('GLOBAL %(G-VECTORS) INCONSISTENT WITH OBJECT')
        CALL ERROR$I4VAL('NG (ON INPUT)',NGGLOB)
        CALL ERROR$I4VAL('NG (IN OBJECT)',NG)
        CALL ERROR$STOP('PLANEWAVE$COLLECT')
      END IF
      IF(NGLOC.NE.NGL) THEN
        CALL ERROR$MSG('LOCAL %(G-VECTORS) INCONSISTENT WITH OBJECT')
        CALL ERROR$STOP('PLANEWAVE$COLLECT')
      END IF
      ALLOCATE(X_LOCTMP(NDIM,NGLX))
      CALL MPE$SYNC(THIS%CID) !THIS IS ONLY USED TO AVOID FILLING UP THE BUFFER
      IF(THISTASK.EQ.1)THEN
!
!       ================================================================
!       == RECEIVE FROM ALL TASKS AND MAP INTO GLOBAL ARRAY           ==
!       ================================================================
        DO FROMTASK=1,NTASKS
          IF(FROMTASK.NE.1) THEN
            CALL MPE$RECEIVE(THIS%CID,FROMTASK,FROMTASK*2,NGFROM)
            IF(NGFROM.GT.NGLX) THEN
              CALL ERROR$MSG('PACKAGE TOO LARGE')
              CALL ERROR$I4VAL('FROMTASK',FROMTASK)
              CALL ERROR$I4VAL('NGFROM',NGFROM)
              CALL ERROR$I4VAL('NGLX',NGLX)
              CALL ERROR$STOP('PLANEWAVE_COLLECT')
            END IF
            CALL MPE$RECEIVE(THIS%CID,FROMTASK,FROMTASK*2+1,X_LOCTMP)
            IGL=0
            DO IG=1,NG
              IF(THIS%TASKOFIG(IG).NE.FROMTASK) CYCLE
              IGL=IGL+1
              XGLOB(:,IG)=X_LOCTMP(:,IGL)
            ENDDO
          ELSE
            IGL=0
            DO IG=1,NG
              IF(THIS%TASKOFIG(IG).NE.1) CYCLE
              IGL=IGL+1
              XGLOB(:,IG)=XLOC(:,IGL)
            ENDDO
          END IF
        ENDDO
      ELSE
!       ================================================================
!       == OTHER TASKS SEND TO TASK ONE                               ==
!       ================================================================
        CALL MPE$SEND(THIS%CID,1,THISTASK*2,NGL)
        X_LOCTMP(:,1:NGL)=XLOC(:,:)
        CALL MPE$SEND(THIS%CID,1,THISTASK*2+1,X_LOCTMP)
      ENDIF
      DEALLOCATE(X_LOCTMP)
!                     CALL TRACE$POP
      RETURN
      END
#END TEMPLATE PLANEWAVE$COLLECT
!
!     .................................................................
      SUBROUTINE PLANEWAVE$RSPACECOLLECTR8(NRL,XL,NRG,XG)
!     ******************************************************************
!     **  PLANEWAVE$COLLECT                                           **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NRL
      INTEGER(4)  ,INTENT(IN) :: NRG
      REAL(8)     ,INTENT(IN) :: XL(NRL)
      REAL(8)     ,INTENT(OUT):: XG(NRG)
      INTEGER(4)              :: NR1L,NR1G,NR2,NR3,NRLX,NRLFROM,NR1START
      REAL(8)     ,ALLOCATABLE:: X_LOCTMP(:)
      INTEGER(4)              :: I23,IRG,IRL,IR1L,IR2,IR3
      INTEGER(4)              :: FROMTASK
      INTEGER(4)                 :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
                     CALL TRACE$PUSH('PLANEWAVE$RSPACECOLLECTR8')
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('LIST NOT INITIALIZED')
        CALL ERROR$STOP('PLANEWAVE$RSPACECOLLECTR8')
      ENDIF
      NR1L=THIS%NR1L
      NR1G=SUM(THIS%NR1LARR)
      NR2=THIS%NR2
      NR3=THIS%NR3
      IF(NRG.NE.NR1G*NR2*NR3) THEN
        CALL ERROR$MSG('GLOBAL %(R-POINTS) INCONSISTENT WITH OBJECT')
        CALL ERROR$STOP('PLANEWAVE$RSPACECOLLECTR8')
      END IF
      IF(NRL.NE.NR1L*NR2*NR3) THEN
        CALL ERROR$MSG('LOCAL %(R-POINTS) INCONSISTENT WITH OBJECT')
        CALL ERROR$STOP('PLANEWAVE$RSPACECOLLECTR8')
      END IF
      XG(:)=0.D0
      NRLX=MAXVAL(THIS%NR1LARR(:))*NR2*NR3
      ALLOCATE(X_LOCTMP(NRLX))
      CALL MPE$SYNC(THIS%CID) !THIS IS ONLY USED TO AVOID FILLING UP THE BUFFER
      IF(THISTASK.EQ.1)THEN
!
!       ================================================================
!       == RECEIVE FROM ALL TASKS AND MAP INTO GLOBAL ARRAY           ==
!       ================================================================
        NR1START=1
        DO FROMTASK=1,NTASKS
          NR1L=THIS%NR1LARR(FROMTASK)
          IF(FROMTASK.NE.1) THEN
            CALL MPE$RECEIVE(THIS%CID,FROMTASK,FROMTASK*2,NRLFROM)
            IF(NRLFROM.GT.NRLX) THEN
              CALL ERROR$MSG('PACKAGE TO LARGE')
              CALL ERROR$STOP('PLANEWAVE$RSPACECOLLECTR8')
            END IF
            CALL MPE$RECEIVE(THIS%CID,FROMTASK,FROMTASK*2+1,X_LOCTMP)
            I23=-1
            DO IR3=1,NR3
              DO IR2=1,NR2
                I23=I23+1   !I23=IR2-1+NR2*(IR3-1)
                IRG=NR1START-1+NR1G*I23
                IRL=NR1L*I23
                DO IR1L=1,NR1L
                  IRL=IRL+1
                  IRG=IRG+1
                  XG(IRG)=X_LOCTMP(IRL)
                ENDDO
              ENDDO
            ENDDO        
          ELSE
            I23=-1
            DO IR3=1,NR3
              DO IR2=1,NR2
                I23=I23+1   !I23=IR2-1+NR2*(IR3-1)
                IRG=NR1START-1+NR1G*I23
                IRL=NR1L*I23
                DO IR1L=1,NR1L
                  IRL=IRL+1
                  IRG=IRG+1
                  XG(IRG)=XL(IRL)
                ENDDO
              ENDDO
            ENDDO        
          END IF
          NR1START=NR1START+NR1L
        ENDDO
      ELSE
!       ================================================================
!       == OTHER TASKS SEND TO TASK ONE                               ==
!       ================================================================
        CALL MPE$SEND(THIS%CID,1,THISTASK*2,NRL)
        X_LOCTMP(1:NRL)=XL(:)
        CALL MPE$SEND(THIS%CID,1,THISTASK*2+1,X_LOCTMP)
      ENDIF
      DEALLOCATE(X_LOCTMP)
                     CALL TRACE$POP
      
      RETURN
      END
!
!     .................................................................
      SUBROUTINE PLANEWAVE$RSPACECOLLECTC8(NRL,XL,NRG,XG)
!     ******************************************************************
!     **  PLANEWAVE$COLLECT                                           **
!     ******************************************************************
      USE PLANEWAVE_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NRL
      INTEGER(4)  ,INTENT(IN) :: NRG
      COMPLEX(8)  ,INTENT(IN) :: XL(NRL)
      COMPLEX(8)  ,INTENT(OUT):: XG(NRG)
      INTEGER(4)              :: NR1L,NR1G,NR2,NR3,NRLX,NRLFROM,NR1START
      REAL(8)     ,ALLOCATABLE:: X_LOCTMP(:)
      INTEGER(4)              :: I23,IRG,IRL,IR1L,IR2,IR3
      INTEGER(4)              :: FROMTASK
      INTEGER(4)                 :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(THIS%CID,NTASKS,THISTASK)
                     CALL TRACE$PUSH('PLANEWAVE$RSPACECOLLECTR8')
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('LIST NOT INITIALIZED')
        CALL ERROR$STOP('PLANEWAVE$RSPACECOLLECTR8')
      ENDIF
      NR1L=THIS%NR1L
      NR1G=SUM(THIS%NR1LARR)
      NR2=THIS%NR2
      NR3=THIS%NR3
      IF(NRG.NE.NR1G*NR2*NR3) THEN
        CALL ERROR$MSG('GLOBAL %(R-POINTS) INCONSISTENT WITH OBJECT')
        CALL ERROR$STOP('PLANEWAVE$RSPACECOLLECTR8')
      END IF
      IF(NRL.NE.NR1L*NR2*NR3) THEN
        CALL ERROR$MSG('LOCAL %(R-POINTS) INCONSISTENT WITH OBJECT')
        CALL ERROR$STOP('PLANEWAVE$RSPACECOLLECTR8')
      END IF
      XG(:)=0.D0
      NRLX=MAXVAL(THIS%NR1LARR(:))*NR2*NR3
      ALLOCATE(X_LOCTMP(NRLX))
      CALL MPE$SYNC(THIS%CID) !THIS IS ONLY USED TO AVOID FILLING UP THE BUFFER
      IF(THISTASK.EQ.1)THEN
!
!       ================================================================
!       == RECEIVE FROM ALL TASKS AND MAP INTO GLOBAL ARRAY           ==
!       ================================================================
        NR1START=1
        DO FROMTASK=1,NTASKS
          NR1L=THIS%NR1LARR(FROMTASK)
          IF(FROMTASK.NE.1) THEN
            CALL MPE$RECEIVE(THIS%CID,FROMTASK,FROMTASK*2,NRLFROM)
            IF(NRLFROM.GT.NRLX) THEN
              CALL ERROR$MSG('PACKAGE TO LARGE')
              CALL ERROR$STOP('PLANEWAVE$RSPACECOLLECTR8')
            END IF
            CALL MPE$RECEIVE(THIS%CID,FROMTASK,FROMTASK*2+1,X_LOCTMP)
            I23=-1
            DO IR3=1,NR3
              DO IR2=1,NR2
                I23=I23+1   !I23=IR2-1+NR2*(IR3-1)
                IRG=NR1START-1+NR1G*I23
                IRL=NR1L*I23
                DO IR1L=1,NR1L
                  IRL=IRL+1
                  IRG=IRG+1
                  XG(IRG)=X_LOCTMP(IRL)
                ENDDO
              ENDDO
            ENDDO        
          ELSE
            I23=-1
            DO IR3=1,NR3
              DO IR2=1,NR2
                I23=I23+1   !I23=IR2-1+NR2*(IR3-1)
                IRG=NR1START-1+NR1G*I23
                IRL=NR1L*I23
                DO IR1L=1,NR1L
                  IRL=IRL+1
                  IRG=IRG+1
                  XG(IRG)=XL(IRL)
                ENDDO
              ENDDO
            ENDDO        
          END IF
          NR1START=NR1START+NR1L
        ENDDO
      ELSE
!       ================================================================
!       == OTHER TASKS SEND TO TASK ONE                               ==
!       ================================================================
        CALL MPE$SEND(THIS%CID,1,THISTASK*2,NRL)
        X_LOCTMP(1:NRL)=XL(:)
        CALL MPE$SEND(THIS%CID,1,THISTASK*2+1,X_LOCTMP)
      ENDIF
      DEALLOCATE(X_LOCTMP)
                     CALL TRACE$POP
      
      RETURN
      END
