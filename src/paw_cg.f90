!     ******************************************************************
!     **                                                              **
!     **                   INTERFACE TO THE                           **
!     **                                                              **
!     **       STATE BY STATE CONJUGATE GRADIENT LIBRARY              **
!     **                                                              **
!     **             IT ACCESSES THE WAVES MODULE                     **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
MODULE CG_INTERFACE_MODULE
INTEGER (4)    :: IKPT
INTEGER (4)    :: ISPIN
INTEGER(4)     :: NAT
INTEGER(4)     :: LMNXX
REAL(8) ,ALLOCATABLE :: DH(:,:,:,:) ! (LMNXX,LMNXX,NDIMD,NAT)
!INTEGER(4)     :: NPRO
INTEGER(4)               :: NRL
REAL(8) ,ALLOCATABLE :: V(:,:) ! (NRL,NDIMD)
LOGICAL(4)               :: TSAVEMEM=.FALSE.
!LOGICAL(4)               :: TSAVEMEM=.TRUE.
! AVOID RECALCULATION OF PROJECTIONS
!COMPLEX(8)  ,ALLOCATABLE :: PROJ(:,:)
LOGICAL(4)               :: TNEWPRO=.TRUE.
! AVOID RECALCULATION OF OPSI(OLD)
COMPLEX(8)  ,ALLOCATABLE :: OPSIOLD(:,:) ! (NB,NGL)
END MODULE CG_INTERFACE_MODULE


!
!     ..................................................................
      SUBROUTINE CG$STATE_BY_STATE(NRL_,NDIMD_,V_,CONV,NAT_,LMNXX_,DH_)
!     ******************************************************************
!     ** INTERFACE SUBROUTINE                                         **
!     ******************************************************************
      USE CG_INTERFACE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)   , INTENT(IN)   :: NRL_
      INTEGER(4)   , INTENT(IN)   :: NDIMD_
      REAL(8)      , INTENT(IN)   :: V_(NRL_,NDIMD_)
      REAL(8)      , INTENT(INOUT)   :: CONV 
      INTEGER(4)   , INTENT(IN)   :: NAT_
      INTEGER(4)   , INTENT(IN)   :: LMNXX_
      REAL(8)      , INTENT(IN)   :: DH_(LMNXX_,LMNXX_,NDIMD_,NAT_)
      LOGICAL(4)                  :: TCONV
      INTEGER(4)                  :: NGL,NB,NBH,IB,NPRO
      INTEGER(4)                  :: NITERX
      REAL(8)     ,ALLOCATABLE    :: G2(:) ! (NGL)
      REAL(8)     ,ALLOCATABLE    :: R0(:,:)
!     ******************************************************************
      IF(NDIM.NE.1) THEN
         CALL ERROR$I4VAL('NDIMD',NDIMD)
         CALL ERROR$MSG('NON-COLLINEAR CG NOT IMPLEMENTED')
         CALL ERROR$STOP('CG$STATE_BY_STATE')
      END IF
      NAT=NAT_
      LMNXX=LMNXX_
      ALLOCATE(DH(LMNXX,LMNXX,NDIMD,NAT))
      DH=DH_
      NPRO=MAP%NPRO
      
      NRL=NRL_
      ALLOCATE(V(NRL,NDIMD_))
      V=V_
WRITE(*,"('CG-MIXER: POT   :',3F10.7)") V(1:3,1)
      NITERX=30 !400
      DO IKPT=1,NKPT
         DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            CALL PLANEWAVE$SELECT(GSET%ID) 
            NGL=GSET%NGL
            NBH=THIS%NBH !TAKE CARE OF SUPER WAVE FUNCTIONS!!
            NB=THIS%NB
            !ALLOCATE(THIS%EIGVAL(NB))
            !ALLOCATE(PROJ(NDIM,NPRO))
            IF((.NOT.TSAVEMEM).AND.(NB.GT.1)) THEN
               PRINT*,'SAVEMEM=FALSE, ALLOCATING ',NB*NGL*16.D0/1024.D0,' KB'
               ALLOCATE(OPSIOLD(NB-1,NGL))
            END IF
            !IF(.NOT.ASSOCIATED(THIS%HPSI))ALLOCATE(THIS%HPSI(NGL,NDIM,NBH))
            IF(NB.NE.NBH) THEN
               CALL ERROR$MSG('SUPERWAVEFUNCTIONS NOT YET IMPLEMENTED')
               CALL ERROR$STOP('CG$STATE_BY_STATE')
            END IF
            IF(ALLOCATED(G2)) DEALLOCATE(G2) !NEW GSET?
            ALLOCATE(G2(NGL))
            CALL PLANEWAVE$GETR8A('G2',NGL,G2)
            !G2(:)=1.D0 ! IF THIS LINE IS USED, PRECONDITIONING IS AVOIDED
!!$            ! IF THE FOLLOWING LOOP IS USED, START FROM RANDOM WF
!!$            DO IB=1,NB
!!$               CONV=1.D-8
!!$               CALL WAVES_RANDOMIZE(NGL,1,1,1.D4,G2,THIS%PSI0(:,1,IB))
!!$            END DO

            CALL CG_STATE_BY_STATE(NGL,NB,NDIM,NPRO,CONV,NITERX,G2,&
                 THIS%EIGVAL, &
                 THIS%PSI0,TCONV)
            CALL TIMING$CLOCKON('CG-PRO2')
            ALLOCATE(R0(3,NAT))
            CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R0,NGL,NDIM,NBH,MAP%NPRO &
                ,THIS%PSI0,THIS%PROJ) 
            DEALLOCATE(R0)
            !DEALLOCATE(PROJ)
            CALL TIMING$CLOCKOFF('CG-PRO2')
            IF(.NOT.TSAVEMEM) THEN
               DEALLOCATE(OPSIOLD)
            END IF
         END DO
      END DO
      DEALLOCATE(DH)
      DEALLOCATE(V)
      PRINT*,'CG FINISHED !!'
      PRINT*,'CG EPS:',THIS%EIGVAL*27.211396
      !CALL ERROR$STOP(' BREAK IN CG$STATE_BY_STATE')
      IF(.NOT.TCONV) THEN
         PRINT*,'CG CYCLE NOT CONVERGED'
         CALL ERROR$MSG('CG CYCLE NOT CONVERGED')
         CALL ERROR$STOP('CG$STATE_BY_STATE')
      END IF
      END SUBROUTINE CG$STATE_BY_STATE
!
!     ..................................................................
      SUBROUTINE CG_OPERATOR_DOT_VEC(ID,NGL,NDIM,NPRO,PSI,OPSI,PROJ)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE CG_INTERFACE_MODULE
      USE WAVES_MODULE, ONLY: MAP,GSET,NDIMD,THIS 
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: NGL
      INTEGER(4)  ,INTENT(IN)  :: NDIM
      INTEGER(4)  ,INTENT(IN)  :: NPRO
      COMPLEX(8)  ,INTENT(IN)  :: PSI(NGL)
      COMPLEX(8)  ,INTENT(OUT) :: OPSI(NGL)
      COMPLEX(8)  ,INTENT(INOUT):: PROJ(NDIM,NPRO)
      REAL(8)     ,ALLOCATABLE :: R0(:,:)
!     ******************************************************************
      IF(ID.EQ.'HAMILTON') THEN
         CALL TIMING$CLOCKON('CG-HPSI')
         ALLOCATE(R0(3,NAT))
         CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!WRITE(*,"('HPSI PROJ 1',3F10.5)") REAL(PROJ(1,1:3))
         IF(TNEWPRO) THEN
            CALL TIMING$CLOCKON('CG-PROJ')
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R0,NGL,1,1,NPRO &
                 ,PSI,PROJ)
            CALL TIMING$CLOCKOFF('CG-PROJ')
         END IF
!WRITE(*,"('HPSI PROJ 2',3F10.5)") REAL(PROJ(1,1:3))
         DEALLOCATE(R0)
         CALL WAVES_HPSI(MAP,GSET,ISPIN,NGL,NDIM,NDIMD,1,NPRO,LMNXX,NAT,NRL &
     &                  ,PSI,V(1,ISPIN),R0,PROJ,DH,OPSI)  !OPSI IST EIGENTLICH HPSI
!        IF(NDIM.EQ.1) THEN
!            CALL WAVES_HPSI(IKPT,ISPIN,NRL,V(:,ISPIN),NAT,LMNXX,DH,&
!                NGL,1,1,NPRO,PROJ,PSI,OPSI)
!         ELSE
!            CALL WAVES_HPSI(IKPT,ISPIN,NRL,V(:,1:NDIMD),NAT,LMNXX,DH,&
!                 NGL,1,1,NPRO,PROJ,PSI,OPSI)
!         END IF
         CALL TIMING$CLOCKOFF('CG-HPSI')
      ELSE IF (ID.EQ.'OVERLAP') THEN
         CALL TIMING$CLOCKON('CG-OPSI')
         ALLOCATE(R0(3,NAT))
         CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
         IF(TNEWPRO) THEN
            CALL TIMING$CLOCKON('CG-PROJ')
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R0,NGL,1,1,NPRO &
                 ,PSI,PROJ)
            CALL TIMING$CLOCKOFF('CG-PROJ')
         END IF
!WRITE(*,"('OPSI PROJ  ',3F10.5)") REAL(PROJ(1,1:3))
         OPSI=PSI
         CALL WAVES_OPSI(1,1,NPRO,NAT,NGL,R0,PROJ,OPSI)
         DEALLOCATE(R0)
         CALL TIMING$CLOCKOFF('CG-OPSI')
      ELSE
         CALL ERROR$MSG('ID IN CG_OPERATOR_DOT_VEC MUST BE HAMILTON OR OVERLAP')
         CALL ERROR$STOP('CG_OPERATOR_DOT_VEC')
      END IF
      END SUBROUTINE CG_OPERATOR_DOT_VEC

!
!     ..................................................................
      SUBROUTINE CG_DOT_PRODUCT(NGL,PSI1,PSI2,CVAR)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NGL
      COMPLEX(8),INTENT(IN)  :: PSI1(NGL)
      COMPLEX(8),INTENT(IN)  :: PSI2(NGL)
      COMPLEX(8),INTENT(OUT) :: CVAR
!     ******************************************************************
      CALL TIMING$CLOCKON('CG-PSIPSI')
      CALL PLANEWAVE$SCALARPRODUCT(' ',NGL,1,1,PSI1,1,PSI2,CVAR)
      CALL MPE$COMBINE('NONE','+',CVAR)
      CALL TIMING$CLOCKOFF('CG-PSIPSI')
      END SUBROUTINE CG_DOT_PRODUCT


!     ******************************************************************
!     ******************************************************************
!     ******************************************************************
!     **                                                              **
!     **       STATE BY STATE CONJUGATE GRADIENT LIBRARY              **
!     **                                                              **
!     ******************************************************************
!     ******************************************************************
!     ******************************************************************
!     **                                                              **
!     ** ON ENTRY THIS SUBROUTINE REQUIRES THE FOLLOWING DATA:        **
!     **   NGL           ... MATRIX SIZE                              **
!     **   NB            ... NUMBER OF EIGENVALUES/VECTORS NEEDED     **
!     **   CONV          ... CONVERGENCE CRITERION FOR EIGENVALUE     **
!     **   NITERX        ... MAX. NUMBER OF CG ITERATIONS             **
!     **   G2(NGL)       ... THE KINETIC ENERGY OF THE BASIS FUNCTIONS**
!     **                     (NEEDED FOR PRECONDITIONING)             **
!     **   EIGVEC(NGL,NB)... INITIAL GUESS FOR EIGENVECTORS           **
!     **                                                              **
!     ** ON EXIT IT PROVIDES THE FOLLOWING DATA:                      **
!     **   EPS(NB)       ... THE FIRST NB EIGENVALUES                 **
!     **   EIGVEC(NGL,NB)... THE FIRST NB EIGENVECTORS                **
!     **   TCONV         ... TRUE IF ALL EIGENVECTORS ARE CONVERGED   **
!     **                                                              **
!     **                                                              **
!     ** MIND! THAT THIS LIBRARY REQUIRES TWO ADDITIONAL SUBROUTINES  **
!     **       FROM THE USER:                                         **
!     **                                                              **
!     **   CG_OPERATOR_DOT_VEC(ID,NGL,PSI(:),OPSI(:))                 **
!     **   WHERE 'ID' CAN BE 'HAMILTON' OR 'OVERLAP'                  **
!     **   EITHER H|PSI> OR O|PSI> MUST BE WRITTEN INTO OPSI(:)       **
!     **                                                              **
!     **   CG_DOT_PRODUCT(NGL,PSI1(:),PSI2(:),CVAR)                   **
!     **   THE SCALAR PRODUCT <PSI1|PSI2> MUST BE WRITTEN INTO CVAR   **
!     **                                                              **
!     **                                                              **
!     ** (C) CLEMENS FOERST, 2004                                     **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
!     ******************************************************************
      SUBROUTINE CG_STATE_BY_STATE(NGL,NB,NDIM,NPRO,CONV,NITERX,G2,EPS,&
           EIGVEC,TCONV)
      USE MPE_MODULE
      USE CG_INTERFACE_MODULE, ONLY: TNEWPRO,TSAVEMEM,OPSIOLD
      IMPLICIT NONE
      INTEGER(4)   , INTENT(IN)   :: NGL
      INTEGER(4)   , INTENT(IN)   :: NB
      INTEGER(4)   , INTENT(IN)   :: NDIM
      INTEGER(4)   , INTENT(IN)   :: NPRO
      REAL(8)      , INTENT(IN)   :: CONV
      INTEGER(4)   , INTENT(IN)   :: NITERX
      REAL(8)      , INTENT(IN)   :: G2(NGL)
      REAL(8)      , INTENT(INOUT):: EPS(NB)
      COMPLEX(8)   , INTENT(INOUT):: EIGVEC(NGL,NB)
      LOGICAL(4)   , INTENT(OUT)  :: TCONV
      INTEGER(4)                  :: ITER
      INTEGER(4)                  :: IB,IB1,IB2
      INTEGER(4)                  :: IG,IG1,IG2
      COMPLEX(8)                  :: HEIG(NGL) 
      COMPLEX(8)                  :: GRAD(NGL) 
      COMPLEX(8)                  :: GRADOLD(NGL) 
      COMPLEX(8)                  :: PRECOND(NGL) 
      COMPLEX(8)                  :: SEARCH(NGL) 
      COMPLEX(8)                  :: SEARCHOLD(NGL) 
      COMPLEX(8)                  :: GAMMA
      REAL(8)                     :: EPSOLD
      REAL(8)                     :: H11,H12,H21,H22
      COMPLEX(8)                  :: CVAR
      REAL(8)                     :: PI
      REAL(8)                     :: THETA
      REAL(8)                     :: SVAR
      COMPLEX(8)                  :: OPSI(NGL)
      COMPLEX(8)                  :: GTIMESGC
      COMPLEX(8)                  :: GTIMESGCOLD
      REAL(8)                     :: K(NGL)
      INTEGER(4)                  :: ITERSUMMARY(NB)
      LOGICAL(4)                  :: TRESTART
      LOGICAL(4)                  :: TLASTCHECK
      LOGICAL(4)                  :: TCONVERGED
!
!     SUBSPACE ROTATION
      COMPLEX(8)                  :: SUBH(NB,NB)
      COMPLEX(8)                  :: SUBO(NB,NB)
      COMPLEX(8)                  :: SUBVEC(NB,NB)
      REAL(8)                     :: SUBEPS(NB)
!
!     DEBUG STUFF
      REAL(8)                     :: EIGVAL
      INTEGER(4)                  :: ISTEP
      INTEGER(4)   , PARAMETER    :: NSTEP=30
      REAL(8)                     :: WORK2D(NB,NB)
      INTEGER(4)                  :: THISTASK,NTASKS

      !JOHANNES
      COMPLEX(8)                  :: PROJPSI(NDIM,NPRO)
      COMPLEX(8)                  :: PROJSEARCH(NDIM,NPRO)
      
!     .................................................................
!
      CALL MPE$QUERY('NONE',NTASKS,THISTASK)
      PI = 4.D0*ATAN(1.D0)
!
!     =================================================================
!     =================================================================
!     =================================================================
!     == LOOP OVER ALL BANDS TO BE OPTIMIZED ==========================
!     =================================================================
!     =================================================================
!     =================================================================
      DO IB = 1, NB
!
!       == ORTHOGONALIZE TO THE ALREADY OPTIMIZED ONES ================
!!$   ! DAS SOLLTE EH UNTEN GEMACHT WERDEN!!!!!----------------
!!$        DO IB1 = 1, IB-1
!!$          CALL CG_INTERNAL_ORTHOGONALIZE(NGL,EIGVEC(:,IB),EIGVEC(:,IB1))
!!$        END DO
!!$        CALL CG_INTERNAL_NORMALIZE(NGL,EIGVEC(:,IB))
!CALL CG_INTERNAL_OVERLAP('OVERLAP',NGL,EIGVEC(:,IB),EIGVEC(:,IB),CVAR)
!PRINT*,'<PSI|PSI>=',CVAR

!
        EPSOLD = EPS(IB) !1.D100 
        TCONVERGED = .FALSE.
        TRESTART = .TRUE.
        TLASTCHECK = .FALSE.   ! WHEN CONVERGENCE IS REACHED, ANOTHER
                               ! STEEPEST DESCENT STEP IS DONE AND CONVERGENCE IS RECHECKED
!
!       ===============================================================
!       ===============================================================
!       == CG LOOP FOR ONE BAND =======================================
!       ===============================================================
!       ===============================================================
        DO ITER = 1, NITERX
!
!         =============================================================
!         == ORTHOGONALIZE WAVE-FUNCTION TO ALL LOWER BANDS          ==
!         == STRICTLY NECESSARY IN FIRST ITERATION, THEN EVERY 10    ==
!         == TO AVOID ACCUMULATION OF ROUNDING ERRORS                ==
!         =============================================================
          IF ((MOD(ITER,10).EQ.1).OR.TLASTCHECK) THEN
            CALL TIMING$CLOCKON('CG-M1-ORTHOGONALIZE')
            DO IB1 = 1, IB-1
              IF(TSAVEMEM) THEN
!IN THE FOLLOWING LINE THE ARGUMENT OPSI WAS MISSING IN THE ORIGINAL VERSION
                CALL CG_INTERNAL_ORTHOGONALIZE(NGL,EIGVEC(:,IB),EIGVEC(:,IB1),OPSI)
              ELSE
                CALL CG_INTERNAL_ORTHO_PSI(NGL,EIGVEC(:,IB),&
                     EIGVEC(:,IB1),IB1)
              END IF
            END DO
            ! CALCULATE P(PSI) HERE
            TNEWPRO=.TRUE.
            CALL CG_INTERNAL_NORMALIZE(NGL,NDIM,NPRO,EIGVEC(:,IB),PROJPSI)
            CALL TIMING$CLOCKOFF('CG-M1-ORTHOGONALIZE')
          END IF
!
          CALL TIMING$CLOCKON('CG-M2-HPSI') !TIME CRITICAL
          TNEWPRO=.FALSE.
          !CALL CG_INTERNAL_MULTIPLY('HAMILTON',NGL,EIGVEC(:,IB),HEIG(:))
          CALL CG_OPERATOR_DOT_VEC('HAMILTON',NGL,NDIM,NPRO,EIGVEC(:,IB),&
               HEIG(:),PROJPSI)
          TNEWPRO=.TRUE.
          CALL CG_INTERNAL_SCALARPRODUCT(NGL,EIGVEC(:,IB),HEIG(:),CVAR)
          EPS(IB) = REAL(CVAR,KIND=8)
          CALL TIMING$CLOCKOFF('CG-M2-HPSI')
!
!         ==============================================================
!         == EXIT WHEN CONVERGED                                      ==
!         == CONVERGENCE CRITERION: ABSOLUTE CONVERGENCE OF EIGENVALUE==
!         == AFTERWARDS, A STEEPEST DESCENT STEP IS PERFORMED TO      ==
!         == RECHECK, IF ENERGY IS STILL WITHIN THE RANGE  -> EXIT    ==
!         ==============================================================
PRINT*,IB,ITER,EPS(IB)

!PRINT*,'TASK ',THISTASK,' RES:',ABS(EPS(IB) - EPSOLD),'CONV:',CONV
!PRINT*,'TASK ',THISTASK,(ABS(EPS(IB) - EPSOLD).LT.CONV)
          IF (ABS(EPS(IB) - EPSOLD).LT.CONV) THEN
            IF (.NOT.TLASTCHECK) THEN
              TLASTCHECK = .TRUE.
PRINT*,'TRIGGERING A LAST STEEPEST DECENT CHECK : TASK ',THISTASK
            ELSE
              PRINT*,'CONVERGED', ITER
              TCONVERGED = .TRUE.
              ITERSUMMARY(IB) = ITER
              EXIT
            END IF  
          ELSE
            TLASTCHECK = .FALSE.
          END IF
          EPSOLD = EPS(IB)
!
!         == COMPUTE GRADIENT =========================================
          CALL TIMING$CLOCKON('CG-M3-GRAD')
          TNEWPRO=.FALSE.
          !CALL CG_INTERNAL_MULTIPLY('OVERLAP',NGL,EIGVEC(:,IB),OPSI(:))
          CALL CG_OPERATOR_DOT_VEC('OVERLAP',NGL,NDIM,NPRO,EIGVEC(:,IB),&
               OPSI(:),PROJPSI)
          TNEWPRO=.TRUE.
          GRAD(:) = HEIG(:) - CMPLX(EPS(IB), KIND=8)*OPSI(:)
          CALL TIMING$CLOCKOFF('CG-M3-GRAD')
!PRINT*,'TASK',THISTASK,'AFTER GRADIENT'
!
!         == PRECONDITIONING ==========================================
          CALL TIMING$CLOCKON('CG-M4-PRECON')
          CALL CG_PRECONDITION_FULL(NGL,G2(:),EIGVEC(:,IB),GRAD(:),PRECOND(:)) 
          CALL TIMING$CLOCKOFF('CG-M4-PRECON')
!PRINT*,'TASK',THISTASK,'AFTER PRECON'
!
!         == ORTHOGONALIZE GRADIENT TO ALREADY OPTIMIZED WAVE-FUNCTIONS
          CALL TIMING$CLOCKON('CG-M5-ORTHOGONALIZE G')
          DO IB1 = 1, IB - 1
            IF(TSAVEMEM) THEN
!PB THE FOLLOWING LINE DOES NOT MAKE SENSE
!               CALL CG_INTERNAL_ORTHOGONALIZE(NGL,PRECOND(:),EIGVEC(:,IB1))
CALL ERROR$MSG('INTERNAL CODE ERROR. STOPPING')
CALL ERROR$STOP('CG_STATE_BY_STATE')
            ELSE
               CALL CG_INTERNAL_ORTHO_PSI(NGL,PRECOND(:),EIGVEC(:,IB1),IB1)
            END IF
          END DO
          CALL TIMING$CLOCKOFF('CG-M5-ORTHOGONALIZE G')
!PRINT*,'TASK',THISTASK,'AFTER ORTHGRAD'
!

          CALL TIMING$CLOCKON('CG-M6-CG')
!         == CONJUGATE GRADIENT =======================================
          CALL CG_DOT_PRODUCT(NGL,GRAD(:),PRECOND(:),GTIMESGC)
          CALL CG_DOT_PRODUCT(NGL,GRADOLD(:),PRECOND(:),CVAR)!TRICK TO IMPROVE CONVERGENCE (E.G. KRESSE DISS)
          IF (TRESTART.OR.(MOD(ITER,20).EQ.1).OR.TLASTCHECK) THEN 
            SEARCH(:) = PRECOND(:)
            TRESTART=.FALSE.
          ELSE
            SEARCH(:) = PRECOND(:) + SEARCHOLD(:)*(GTIMESGC-CVAR)/GTIMESGCOLD
          END IF
!          CALL CG_DOT_PRODUCT(NGL,GRAD(:),PRECOND(:),GTIMESGCOLD)
          GTIMESGCOLD = GTIMESGC !CANGE ON AUG 25 2004
          SEARCHOLD(:) = SEARCH(:)
          GRADOLD(:) = GRAD(:)
          CALL TIMING$CLOCKOFF('CG-M6-CG')
!        
!         == ORTHOGONALIZE TO ACTUAL WAVE-FUNCTION ====================
          CALL TIMING$CLOCKON('CG-M7-ORTHSEARCH1')
          CALL CG_INTERNAL_ORTHOGONALIZE(NGL,SEARCH(:),EIGVEC(:,IB),OPSI)
          ! HERE, P(SEARCH) HAVE TO BE CALCULATED
          CALL TIMING$CLOCKOFF('CG-M7-ORTHSEARCH1')
          TNEWPRO=.TRUE.
          CALL TIMING$CLOCKON('CG-M7-ORTHSEARCH2') ! TIME CRITICAL
          CALL CG_INTERNAL_NORMALIZE(NGL,NDIM,NPRO,SEARCH(:),PROJSEARCH)
          CALL TIMING$CLOCKOFF('CG-M7-ORTHSEARCH2')
!
!         == CALCULATE TRANSFORMATION ANGLE THETA -====================
!         == MIND: ONE SHOULD ALSO CONSIDER OVERLAP!!! ================
          CALL TIMING$CLOCKON('CG-M8-THETA') !TIME CRITICAL
          H11 = EPS(IB)
!PRINT*,'CALL3'
          CALL CG_INTERNAL_SCALARPRODUCT(NGL,SEARCH(:),HEIG(:),CVAR)
          H12 = REAL(CVAR, KIND=8)  
          H21 = H12
!PRINT*,'CALL H|SEARCH>'
          TNEWPRO=.FALSE.
          CALL CG_INTERNAL_OVERLAP('HAMILTON',NGL,NDIM,NPRO,SEARCH(:),&
               SEARCH(:),PROJSEARCH,CVAR)
          TNEWPRO=.TRUE.
          H22 = REAL(CVAR, KIND=8)
!
          THETA = 0.5D0*ATAN(((H21 + H12)/(H11 - H22)))
          IF (H11.GT.H22) THETA = THETA + PI/2.D0
          CALL TIMING$CLOCKOFF('CG-M8-THETA')
!         == PROPAGATE WAVE-FUNCTION ==================================
          CALL TIMING$CLOCKON('CG-M9-PROPAGATE')
          EIGVEC(:,IB) = CMPLX(COS(THETA), KIND=8)*EIGVEC(:,IB) + &
                         CMPLX(SIN(THETA), KIND=8)*SEARCH(:)
          PROJPSI(:,:) = CMPLX(COS(THETA), KIND=8)*PROJPSI(:,:) + &
                         CMPLX(SIN(THETA), KIND=8)*PROJSEARCH(:,:)
          CALL TIMING$CLOCKOFF('CG-M9-PROPAGATE')
!
        END DO
        IF((.NOT.TSAVEMEM).AND.(IB.LT.NB)) THEN
           ! SAVE OPSI(IB) IN OPSIOLD
           !CALL CG_INTERNAL_MULTIPLY('OVERLAP',NGL,EIGVEC(:,IB),OPSIOLD(IB,:))
           CALL CG_OPERATOR_DOT_VEC('OVERLAP',NGL,NDIM,NPRO,EIGVEC(:,IB),&
                OPSIOLD(IB,:),PROJPSI)
        END IF
           
        IF (.NOT.TCONVERGED) THEN
          ITERSUMMARY(IB) = ITER-1
          TCONV = .FALSE.
          !RETURN
        END IF
PRINT*,"------------------------NEXT-BAND-----------------------------"
      END DO
!
!     ==================================================================
!     == ANALYSIS OF RESULT: PRINT SUBSPACE MATRICES AFTER CG ==========
!     ==================================================================
!      PRINT*,'SUBSPACE OVERLAP'
!      DO IB1 = 1, NB
!        DO IB2 = 1, NB
!          CALL CG_INTERNAL_OVERLAP('OVERLAP',NGL,EIGVEC(:,IB1),EIGVEC(:,IB2),CVAR)
!          WRITE(*,'(2E10.3,2X)',ADVANCE='NO') CVAR
!        END DO
!        WRITE(*,*)
!      END DO
!      PRINT*
!      PRINT*,'SUBSPACE HAMILTON'
!      DO IB1 = 1, NB
!        DO IB2 = 1, NB
!          CALL CG_INTERNAL_OVERLAP('HAMILTON',NGL,EIGVEC(:,IB1),EIGVEC(:,IB2),CVAR)
!          WRITE(*,'(2E10.3,2X)',ADVANCE='NO') CVAR
!        END DO
!        WRITE(*,*)
!      END DO
!
WRITE(*,ADVANCE='NO',FMT="(A)") 'CG: ITERSUMMARY: '
DO IB=1,NB
   WRITE(*,ADVANCE='NO',FMT="(I3)") ITERSUMMARY(IB)
END DO
WRITE(*,*)
!PRINT*,'CG: ITERSUMMARY: ',ITERSUMMARY(:)
!
      TCONV = .TRUE.
      RETURN
      END SUBROUTINE
!
!     ..................................................................
      SUBROUTINE CG_PRECONDITION_FULL(NGL,G2,EIGVEC,GRAD,PRECOND)
!     ******************************************************************
!     ** PRECONDITIONS VECTOR GRAD AND RETURNS VECTOR PRECOND         **
!     **                 AFTER KRESSE ET AL., PRB 54, 11169 (1996)    **
!     **                 OR    TETER ET AL.,  PRB 40, 12255 (1989)    **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)       , INTENT(IN)    :: NGL 
      REAL(8)          , INTENT(IN)    :: G2(NGL)
      COMPLEX(8)       , INTENT(IN)    :: EIGVEC(NGL)
      COMPLEX(8)       , INTENT(IN)    :: GRAD(NGL)
      COMPLEX(8)       , INTENT(OUT)   :: PRECOND(NGL)
      REAL(8)                          :: EKIN
      INTEGER(4)                       :: IG
      REAL(8)                          :: SVAR
      REAL(8)                          :: X
      REAL(8)                          :: FACT
      REAL(8)                          :: PREFACT
      COMPLEX(8)                       :: CVAR
!     ******************************************************************
      EKIN = 0.D0
      DO IG = 1, NGL
        CVAR = EIGVEC(IG)
        EKIN = EKIN + 0.5D0*G2(IG)*REAL(CONJG(CVAR)*CVAR,KIND=8)
      END DO
! 
!     == IMPROVED CONVERGENCE THIS WAY - EMPIRICAL PARAMETER ===========
      IF (EKIN.LT.0.5D0) EKIN = 0.5D0
!
!A-LA KRESSE
!      PREFACT = 2.D0/(1.5D0*EKIN)
!      FACT = 0.5D0/(1.5D0*EKIN)
!
!A-LA TETER
      PREFACT = 1.D0
      FACT = 0.5D0/EKIN
!
      DO IG = 1, NGL
        X = G2(IG)*FACT
        SVAR = 27.D0 + 18.D0*X + 12.D0*X**2 + 8.D0*X**3
        SVAR = PREFACT * SVAR / (SVAR + 16.D0*X**4)
        PRECOND(IG) = CMPLX(SVAR, KIND=8)*GRAD(IG)
      END DO
      RETURN
      END SUBROUTINE
!
!     ..................................................................
      SUBROUTINE CG_INTERNAL_MULTIPLY_DEL(ID,NGL,PSI,OPSI)
!     ******************************************************************
!     ** CALCULATES  H|PSI>                                           **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*)     , INTENT(IN)    :: ID
      INTEGER(4)       , INTENT(IN)    :: NGL 
      COMPLEX(8)       , INTENT(IN)    :: PSI(NGL)
      COMPLEX(8)       , INTENT(OUT)   :: OPSI(NGL)
!     ******************************************************************
!      CALL CG_OPERATOR_DOT_VEC(ID,NGL,PSI(:),OPSI(:))
      CALL ERROR$MSG('CODE ERROR')
      CALL ERROR$MSG('INCONSISTENCE CALLING SEQUENCE FOR CG_OPERATOR_DOT_VEC')
      CALL ERROR$STOP('CG_INTERNAL_MULTIPLY_DEL')
      RETURN
      END SUBROUTINE
!
!     ..................................................................
      SUBROUTINE CG_INTERNAL_SCALARPRODUCT(NGL,PSI1,PSI2,CVAR)
!     ******************************************************************
!     ** CALCULATES THE SCALAR PRODUCT    <PSI1|PSI2>                 **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)       , INTENT(IN)    :: NGL 
      COMPLEX(8)       , INTENT(IN)    :: PSI1(NGL)
      COMPLEX(8)       , INTENT(IN)    :: PSI2(NGL)
      COMPLEX(8)       , INTENT(OUT)   :: CVAR
!     ******************************************************************
      CALL CG_DOT_PRODUCT(NGL,PSI1(:),PSI2(:),CVAR)
      RETURN
      END SUBROUTINE
!
!     ..................................................................
      SUBROUTINE CG_INTERNAL_OVERLAP(ID,NGL,NDIM,NPRO,PSI1,PSI2,PROJ,CVAR)
!     ******************************************************************
!     ** CALCULATES THE EXPECTATION VALUE <PSI|O|PSI>                 **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*)     , INTENT(IN)    :: ID
      INTEGER(4)       , INTENT(IN)    :: NGL 
      INTEGER(4)       , INTENT(IN)    :: NDIM 
      INTEGER(4)       , INTENT(IN)    :: NPRO
      COMPLEX(8)       , INTENT(IN)    :: PSI1(NGL)
      COMPLEX(8)       , INTENT(IN)    :: PSI2(NGL)
      COMPLEX(8)       , INTENT(INOUT) :: PROJ(NDIM,NPRO)
      COMPLEX(8)       , INTENT(OUT)   :: CVAR
      COMPLEX(8)                       :: OPSI(NGL)
!     ******************************************************************
      CALL CG_OPERATOR_DOT_VEC(ID,NGL,NDIM,NPRO,PSI2(:),OPSI(:),PROJ)
      CALL CG_DOT_PRODUCT(NGL,PSI1(:),OPSI(:),CVAR)
      RETURN
      END SUBROUTINE
!
!     ..................................................................
      SUBROUTINE CG_INTERNAL_EXPECTVAL_WEG(ID,NGL,PSI,EVAL)
!     ******************************************************************
!     ** CALCULATES THE EXPECTATION VALUE <PSI|H|PSI>                 **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*)     , INTENT(IN)    :: ID
      INTEGER(4)       , INTENT(IN)    :: NGL 
      COMPLEX(8)       , INTENT(IN)    :: PSI(NGL)
      REAL(8)          , INTENT(OUT)   :: EVAL
      COMPLEX(8)                       :: CVAR
!     ******************************************************************
!     CALL CG_INTERNAL_OVERLAP(ID,NGL,PSI(:),PSI(:),CVAR)
      CALL ERROR$MSG('CODE ERROR')
      CALL ERROR$MSG('INCONSISTENT ARGUMENT LIST')
      CALL ERROR$STOP('CG_INTERNAL_EXPECTVAL_WEG')
      IF (ABS(AIMAG(CVAR)).GT.1.D-10) THEN
        PRINT*,'WARNING: EXPECTATION VALUE IS COMPLEX!!!!!'
      END IF
      EVAL = REAL(CVAR, KIND=8)
      RETURN
      END SUBROUTINE
!
!     ..................................................................
      SUBROUTINE CG_INTERNAL_ORTHOGONALIZE(NGL,SEARCH,PSI,OPSI)
!     ******************************************************************
!     ** ORTHOGONALIZES PSI1 TO PSI2                                  **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)       , INTENT(IN)    :: NGL 
      COMPLEX(8)       , INTENT(INOUT) :: SEARCH(NGL)
      COMPLEX(8)       , INTENT(IN)    :: PSI(NGL)
      COMPLEX(8)       , INTENT(IN)    :: OPSI(NGL)
      COMPLEX(8)                       :: CVAR
!     ******************************************************************
      !CALL CG_INTERNAL_OVERLAP('OVERLAP',NGL,PSI2(:),PSI1(:),CVAR)
      CALL CG_DOT_PRODUCT(NGL,SEARCH,OPSI,CVAR)
      SEARCH(:) = SEARCH(:) - CVAR*PSI(:)
      RETURN
      END SUBROUTINE
!
!     ..................................................................
      SUBROUTINE CG_INTERNAL_ORTHO_PSI(NGL,PSI1,PSI2,IB2)
!     ******************************************************************
!     ** ORTHOGONALIZES PSI1 TO PSI2                                  **
!     ******************************************************************
      USE CG_INTERFACE_MODULE, ONLY: OPSIOLD
      IMPLICIT NONE
      INTEGER(4)       , INTENT(IN)    :: NGL 
      COMPLEX(8)       , INTENT(INOUT) :: PSI1(NGL)
      COMPLEX(8)       , INTENT(IN)    :: PSI2(NGL)
      INTEGER(4)       , INTENT(IN)    :: IB2
      COMPLEX(8)                       :: CVAR   ,CVAR2
!     ******************************************************************
      !CALL CG_INTERNAL_OVERLAP('OVERLAP',NGL,PSI2(:),PSI1(:),CVAR)
      CALL CG_DOT_PRODUCT(NGL,PSI1(:),OPSIOLD(IB2,:),CVAR)
      PSI1(:) = PSI1(:) - CONJG(CVAR)*PSI2(:)
      RETURN
      END SUBROUTINE
!
!     ..................................................................
      SUBROUTINE CG_INTERNAL_NORMALIZE(NGL,NDIM,NPRO,PSI,PROJ)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)       , INTENT(IN)    :: NGL 
      INTEGER(4)       , INTENT(IN)    :: NDIM
      INTEGER(4)       , INTENT(IN)    :: NPRO
      COMPLEX(8)       , INTENT(INOUT) :: PSI(NGL)
      COMPLEX(8)       , INTENT(INOUT) :: PROJ(NDIM,NPRO)
      COMPLEX(8)                       :: CVAR
!     ******************************************************************
      CALL CG_INTERNAL_OVERLAP('OVERLAP',NGL,NDIM,NPRO,PSI(:),PSI(:),PROJ,CVAR)
!PRINT*,'PAR CG_INTERNAL_NORMALIZE',CVAR
!WRITE(*,"('NORM PROJ 1',3F10.5)") REAL(PROJ(1,1:3))
      PSI(:) = PSI(:) / SQRT(CVAR)
      PROJ=PROJ / SQRT(CVAR)
!WRITE(*,"('NORM PROJ 2',3F10.5)") REAL(PROJ(1,1:3))
      RETURN
      END SUBROUTINE CG_INTERNAL_NORMALIZE




