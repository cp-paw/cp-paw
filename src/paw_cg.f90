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
logical(4)               :: tsavemem=.false.
!logical(4)               :: tsavemem=.true.
! avoid recalculation of projections
!COMPLEX(8)  ,ALLOCATABLE :: PROJ(:,:)
logical(4)               :: tnewpro=.true.
! avoid recalculation of opsi(old)
COMPLEX(8)  ,ALLOCATABLE :: opsiold(:,:) ! (nb,ngl)
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
      REAL(8)      , INTENT(INout)   :: CONV 
      INTEGER(4)   , INTENT(IN)   :: NAT_
      INTEGER(4)   , INTENT(IN)   :: LMNXX_
      REAL(8)      , INTENT(IN)   :: DH_(LMNXX_,LMNXX_,NDIMD_,NAT_)
      LOGICAL(4)                  :: TCONV
      INTEGER(4)                  :: NGL,NB,NBH,ib,npro
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
      ALLOCATE(DH(LMNXX,LMNXX,ndimd,NAT))
      DH=DH_
      NPRO=MAP%NPRO
      
      NRL=NRL_
      ALLOCATE(V(NRL,NDIMD_))
      V=V_
write(*,"('CG-Mixer: pot   :',3f10.7)") v(1:3,1)
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
            if((.not.tsavemem).and.(nb.gt.1)) then
               print*,'savemem=false, allocating ',nb*ngl*16.D0/1024.D0,' KB'
               allocate(opsiold(nb-1,ngl))
            end if
            !IF(.NOT.ASSOCIATED(THIS%HPSI))ALLOCATE(THIS%HPSI(NGL,NDIM,NBH))
            IF(NB.NE.NBH) THEN
               CALL ERROR$MSG('SUPERWAVEFUNCTIONS NOT YET IMPLEMENTED')
               CALL ERROR$STOP('CG$STATE_BY_STATE')
            END IF
            if(allocated(g2)) deallocate(g2) !new gset?
            ALLOCATE(G2(NGL))
            CALL PLANEWAVE$GETR8A('G2',NGL,G2)
            !G2(:)=1.D0 ! IF THIS LINE IS USED, PRECONDITIONING IS AVOIDED
!!$            ! if the following loop is used, start from random WF
!!$            do ib=1,nb
!!$               conv=1.D-8
!!$               CALL WAVES_RANDOMIZE(NGL,1,1,1.D4,G2,THIS%PSI0(:,1,ib))
!!$            end do

            CALL CG_STATE_BY_STATE(NGL,NB,ndim,npro,CONV,NITERX,G2,&
                 THIS%EIGVAL, &
                 THIS%PSI0,TCONV)
            call timing$clockon('cg-pro2')
            ALLOCATE(R0(3,NAT))
            CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R0,NGL,NDIM,NBH,map%NPRO &
                ,THIS%PSI0,THIS%PROJ) 
            deallocate(r0)
            !deALLOCATE(PROJ)
            call timing$clockoff('cg-pro2')
            if(.not.tsavemem) then
               deallocate(opsiold)
            end if
         END DO
      END DO
      DEALLOCATE(DH)
      DEALLOCATE(V)
      PRINT*,'CG FINISHED !!'
      PRINT*,'CG EPS:',THIS%EIGVAL*27.211396
      !CALL ERROR$STOP(' BREAK IN CG$STATE_BY_STATE')
      IF(.NOT.TCONV) THEN
         print*,'CG CYCLE NOT CONVERGED'
         CALL ERROR$MSG('CG CYCLE NOT CONVERGED')
         CALL ERROR$STOP('CG$STATE_BY_STATE')
      END IF
      END SUBROUTINE CG$STATE_BY_STATE
!
!     ..................................................................
      SUBROUTINE CG_OPERATOR_DOT_VEC(ID,NGL,ndim,npro,PSI,OPSI,proj)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE CG_INTERFACE_MODULE
      USE WAVES_MODULE, ONLY: MAP,GSET,NDIMD,THIS 
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: NGL
      INTEGER(4)  ,INTENT(IN)  :: Ndim
      INTEGER(4)  ,INTENT(IN)  :: Npro
      COMPLEX(8)  ,INTENT(IN)  :: PSI(NGL)
      COMPLEX(8)  ,INTENT(OUT) :: OPSI(NGL)
      COMPLEX(8)  ,INTENT(inOUT):: proj(ndim,npro)
      REAL(8)     ,ALLOCATABLE :: R0(:,:)
!     ******************************************************************
      IF(ID.EQ.'HAMILTON') THEN
         call timing$clockon('cg-hpsi')
         ALLOCATE(R0(3,NAT))
         CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!write(*,"('Hpsi proj 1',3f10.5)") real(proj(1,1:3))
         if(tnewpro) then
            call timing$clockon('cg-proj')
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R0,NGL,1,1,NPRO &
                 ,PSI,PROJ)
            call timing$clockoff('cg-proj')
         end if
!write(*,"('Hpsi proj 2',3f10.5)") real(proj(1,1:3))
         DEALLOCATE(R0)
         CALL WAVES_HPSI(MAP,GSET,ISPIN,NGL,NDIM,NDIMD,1,NPRO,LMNXX,NAT,NRL &
     &                  ,PSI,v(1,ISPIN),R0,PROJ,DH,OPSI)  !opsi ist eigentlich hpsi
!        IF(NDIM.EQ.1) THEN
!            CALL WAVES_HPSI(IKPT,ISPIN,NRL,V(:,ISPIN),NAT,LMNXX,DH,&
!                NGL,1,1,NPRO,PROJ,PSI,OPSI)
!         ELSE
!            CALL WAVES_HPSI(IKPT,ISPIN,NRL,V(:,1:NDIMD),NAT,LMNXX,DH,&
!                 NGL,1,1,NPRO,PROJ,PSI,OPSI)
!         END IF
         call timing$clockoff('cg-hpsi')
      ELSE IF (ID.EQ.'OVERLAP') THEN
         call timing$clockon('cg-opsi')
         ALLOCATE(R0(3,NAT))
         CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
         if(tnewpro) then
            call timing$clockon('cg-proj')
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R0,NGL,1,1,NPRO &
                 ,PSI,PROJ)
            call timing$clockoff('cg-proj')
         end if
!write(*,"('opsi proj  ',3f10.5)") real(proj(1,1:3))
         OPSI=PSI
         CALL WAVES_OPSI(1,1,NPRO,NAT,NGL,R0,PROJ,OPSI)
         DEALLOCATE(R0)
         call timing$clockoff('cg-opsi')
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
      call timing$clockon('cg-psipsi')
      CALL PLANEWAVE$SCALARPRODUCT(' ',NGL,1,1,PSI1,1,PSI2,CVAR)
      CALL MPE$COMBINE('+',cvar)
      call timing$clockoff('cg-psipsi')
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
      SUBROUTINE CG_STATE_BY_STATE(NGL,NB,ndim,npro,CONV,NITERX,G2,EPS,&
           EIGVEC,TCONV)
      use mpe_module
      use CG_INTERFACE_MODULE, only: tnewpro,tsavemem,opsiold
      IMPLICIT NONE
      INTEGER(4)   , INTENT(IN)   :: NGL
      INTEGER(4)   , INTENT(IN)   :: NB
      INTEGER(4)   , INTENT(IN)   :: Ndim
      INTEGER(4)   , INTENT(IN)   :: Npro
      REAL(8)      , INTENT(IN)   :: CONV
      INTEGER(4)   , INTENT(IN)   :: NITERX
      REAL(8)      , INTENT(IN)   :: G2(NGL)
      REAL(8)      , INTENT(inOUT):: EPS(NB)
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
      integer(4)                  :: thistask,ntasks

      !johannes
      COMPLEX(8)                  :: projpsi(ndim,npro)
      COMPLEX(8)                  :: projsearch(ndim,npro)
      
!     .................................................................
!
      call mpe$query(ntasks,thistask)
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
!!$   ! das sollte eh unten gemacht werden!!!!!----------------
!!$        DO IB1 = 1, IB-1
!!$          CALL CG_INTERNAL_ORTHOGONALIZE(NGL,EIGVEC(:,IB),EIGVEC(:,IB1))
!!$        END DO
!!$        CALL CG_INTERNAL_NORMALIZE(NGL,EIGVEC(:,IB))
!CALL CG_INTERNAL_OVERLAP('OVERLAP',NGL,EIGVEC(:,IB),EIGVEC(:,IB),CVAR)
!print*,'<Psi|Psi>=',cvar

!
        EPSOLD = eps(ib) !1.D100 
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
            call timing$clockon('cg-m1-orthogonalize')
            DO IB1 = 1, IB-1
              if(tsavemem) then
                CALL CG_INTERNAL_ORTHOGONALIZE(NGL,EIGVEC(:,IB),EIGVEC(:,IB1))
              else
                CALL CG_INTERNAL_ORTHO_psi(NGL,EIGVEC(:,IB),&
                     EIGVEC(:,IB1),IB1)
              end if
            END DO
            ! calculate P(PSI) here
            tnewpro=.true.
            CALL CG_INTERNAL_NORMALIZE(NGL,ndim,npro,EIGVEC(:,IB),projpsi)
            call timing$clockoff('cg-m1-orthogonalize')
          END IF
!
          call timing$clockon('cg-m2-hpsi') !time critical
          tnewpro=.false.
          !CALL CG_INTERNAL_MULTIPLY('HAMILTON',NGL,EIGVEC(:,IB),HEIG(:))
          CALL CG_OPERATOR_DOT_VEC('HAMILTON',NGL,ndim,npro,EIGVEC(:,IB),&
               HEIG(:),projpsi)
          tnewpro=.true.
          CALL CG_INTERNAL_SCALARPRODUCT(NGL,EIGVEC(:,IB),HEIG(:),CVAR)
          EPS(IB) = REAL(CVAR,KIND=8)
          call timing$clockoff('cg-m2-hpsi')
!
!         ==============================================================
!         == EXIT WHEN CONVERGED                                      ==
!         == CONVERGENCE CRITERION: ABSOLUTE CONVERGENCE OF EIGENVALUE==
!         == AFTERWARDS, A STEEPEST DESCENT STEP IS PERFORMED TO      ==
!         == RECHECK, IF ENERGY IS STILL WITHIN THE RANGE  -> EXIT    ==
!         ==============================================================
PRINT*,IB,ITER,EPS(IB)

!print*,'task ',thistask,' res:',ABS(EPS(IB) - EPSOLD),'conv:',conv
!print*,'task ',thistask,(ABS(EPS(IB) - EPSOLD).LT.CONV)
          IF (ABS(EPS(IB) - EPSOLD).LT.CONV) THEN
            IF (.NOT.TLASTCHECK) THEN
              TLASTCHECK = .TRUE.
PRINT*,'TRIGGERING A LAST STEEPEST DECENT CHECK : task ',thistask
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
          call timing$clockon('cg-m3-grad')
          tnewpro=.false.
          !CALL CG_INTERNAL_MULTIPLY('OVERLAP',NGL,EIGVEC(:,IB),OPSI(:))
          call CG_OPERATOR_DOT_VEC('OVERLAP',NGL,ndim,npro,EIGVEC(:,IB),&
               OPSI(:),projpsi)
          tnewpro=.true.
          GRAD(:) = HEIG(:) - CMPLX(EPS(IB), KIND=8)*OPSI(:)
          call timing$clockoff('cg-m3-grad')
!print*,'task',thistask,'after gradient'
!
!         == PRECONDITIONING ==========================================
          call timing$clockon('cg-m4-precon')
          CALL CG_PRECONDITION_FULL(NGL,G2(:),EIGVEC(:,IB),GRAD(:),PRECOND(:)) 
          call timing$clockoff('cg-m4-precon')
!print*,'task',thistask,'after precon'
!
!         == ORTHOGONALIZE GRADIENT TO ALREADY OPTIMIZED WAVE-FUNCTIONS
          call timing$clockon('cg-m5-orthogonalize G')
          DO IB1 = 1, IB - 1
            if(tsavemem) then
               CALL CG_INTERNAL_ORTHOGONALIZE(NGL,PRECOND(:),EIGVEC(:,IB1))
            else
               CALL CG_INTERNAL_ORTHO_psi(NGL,PRECOND(:),&
                    EIGVEC(:,IB1),IB1)
            end if
          END DO
          call timing$clockoff('cg-m5-orthogonalize G')
!print*,'task',thistask,'after orthgrad'
!

          call timing$clockon('cg-m6-cg')
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
          gtimesgcold = gtimesgc !cange on Aug 25 2004
          SEARCHOLD(:) = SEARCH(:)
          GRADOLD(:) = GRAD(:)
          call timing$clockoff('cg-m6-cg')
!        
!         == ORTHOGONALIZE TO ACTUAL WAVE-FUNCTION ====================
          call timing$clockon('cg-m7-orthsearch1')
          CALL CG_INTERNAL_ORTHOGONALIZE(NGL,SEARCH(:),EIGVEC(:,IB),opsi)
          ! here, P(search) have to be calculated
          call timing$clockoff('cg-m7-orthsearch1')
          tnewpro=.true.
          call timing$clockon('cg-m7-orthsearch2') ! time critical
          CALL CG_INTERNAL_NORMALIZE(NGL,ndim,npro,SEARCH(:),projsearch)
          call timing$clockoff('cg-m7-orthsearch2')
!
!         == CALCULATE TRANSFORMATION ANGLE THETA -====================
!         == MIND: ONE SHOULD ALSO CONSIDER OVERLAP!!! ================
          call timing$clockon('cg-m8-theta') !time critical
          H11 = EPS(IB)
!PRINT*,'CALL3'
          CALL CG_INTERNAL_SCALARPRODUCT(NGL,SEARCH(:),HEIG(:),CVAR)
          H12 = REAL(CVAR, KIND=8)  
          H21 = H12
!PRINT*,'CALL H|SEARCH>'
          tnewpro=.false.
          CALL CG_INTERNAL_OVERLAP('HAMILTON',NGL,ndim,npro,SEARCH(:),&
               SEARCH(:),projsearch,CVAR)
          tnewpro=.true.
          H22 = REAL(CVAR, KIND=8)
!
          THETA = 0.5D0*ATAN(((H21 + H12)/(H11 - H22)))
          IF (H11.GT.H22) THETA = THETA + PI/2.D0
          call timing$clockoff('cg-m8-theta')
!         == PROPAGATE WAVE-FUNCTION ==================================
          call timing$clockon('cg-m9-propagate')
          EIGVEC(:,IB) = CMPLX(COS(THETA), KIND=8)*EIGVEC(:,IB) + &
                         CMPLX(SIN(THETA), KIND=8)*SEARCH(:)
          projpsi(:,:) = CMPLX(COS(THETA), KIND=8)*projpsi(:,:) + &
                         CMPLX(SIN(THETA), KIND=8)*projsearch(:,:)
          call timing$clockoff('cg-m9-propagate')
!
        END DO
        if((.not.tsavemem).and.(ib.lt.nb)) then
           ! save opsi(ib) in opsiold
           !CALL CG_INTERNAL_MULTIPLY('OVERLAP',NGL,EIGVEC(:,IB),opsiold(ib,:))
           call CG_OPERATOR_DOT_VEC('OVERLAP',NGL,ndim,npro,EIGVEC(:,IB),&
                opsiold(ib,:),projpsi)
        end if
           
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
write(*,ADVANCE='NO',FMT="(a)") 'CG: ITERSUMMARY: '
do ib=1,nb
   write(*,ADVANCE='NO',FMT="(i3)") ITERSUMMARY(ib)
end do
write(*,*)
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
      SUBROUTINE CG_INTERNAL_MULTIPLY_del(ID,NGL,PSI,OPSI)
!     ******************************************************************
!     ** CALCULATES  H|PSI>                                           **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*)     , INTENT(IN)    :: ID
      INTEGER(4)       , INTENT(IN)    :: NGL 
      COMPLEX(8)       , INTENT(IN)    :: PSI(NGL)
      COMPLEX(8)       , INTENT(OUT)   :: OPSI(NGL)
!     ******************************************************************
      CALL CG_OPERATOR_DOT_VEC(ID,NGL,PSI(:),OPSI(:))
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
      SUBROUTINE CG_INTERNAL_OVERLAP(ID,NGL,ndim,npro,PSI1,PSI2,proj,CVAR)
!     ******************************************************************
!     ** CALCULATES THE EXPECTATION VALUE <PSI|O|PSI>                 **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*)     , INTENT(IN)    :: ID
      INTEGER(4)       , INTENT(IN)    :: NGL 
      INTEGER(4)       , INTENT(IN)    :: Ndim 
      INTEGER(4)       , INTENT(IN)    :: Npro
      COMPLEX(8)       , INTENT(IN)    :: PSI1(NGL)
      COMPLEX(8)       , INTENT(IN)    :: PSI2(NGL)
      COMPLEX(8)       , INTENT(INout) :: Proj(ndim,npro)
      COMPLEX(8)       , INTENT(OUT)   :: CVAR
      COMPLEX(8)                       :: OPSI(NGL)
!     ******************************************************************
      CALL CG_OPERATOR_DOT_VEC(ID,NGL,ndim,npro,PSI2(:),OPSI(:),proj)
      CALL CG_DOT_PRODUCT(NGL,PSI1(:),OPSI(:),CVAR)
      RETURN
      END SUBROUTINE
!
!     ..................................................................
      SUBROUTINE CG_INTERNAL_EXPECTVAL_weg(ID,NGL,PSI,EVAL)
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
      CALL CG_INTERNAL_OVERLAP(ID,NGL,PSI(:),PSI(:),CVAR)
      IF (ABS(AIMAG(CVAR)).GT.1.D-10) THEN
        PRINT*,'WARNING: EXPECTATION VALUE IS COMPLEX!!!!!'
      END IF
      EVAL = REAL(CVAR, KIND=8)
      RETURN
      END SUBROUTINE
!
!     ..................................................................
      SUBROUTINE CG_INTERNAL_ORTHOGONALIZE(NGL,search,PSI,opsi)
!     ******************************************************************
!     ** ORTHOGONALIZES PSI1 TO PSI2                                  **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)       , INTENT(IN)    :: NGL 
      COMPLEX(8)       , INTENT(INOUT) :: search(NGL)
      COMPLEX(8)       , INTENT(IN)    :: PSI(NGL)
      complex(8)       , intent(in)    :: opsi(ngl)
      COMPLEX(8)                       :: CVAR
!     ******************************************************************
      !CALL CG_INTERNAL_OVERLAP('OVERLAP',NGL,PSI2(:),PSI1(:),CVAR)
      CALL CG_DOT_PRODUCT(NGL,search,OPSI,CVAR)
      search(:) = search(:) - CVAR*PSI(:)
      RETURN
      END SUBROUTINE
!
!     ..................................................................
      SUBROUTINE CG_INTERNAL_ORTHO_psi(NGL,PSI1,psi2,ib2)
!     ******************************************************************
!     ** ORTHOGONALIZES PSI1 TO PSI2                                  **
!     ******************************************************************
      use CG_INTERFACE_MODULE, only: opsiold
      IMPLICIT NONE
      INTEGER(4)       , INTENT(IN)    :: NGL 
      COMPLEX(8)       , INTENT(INOUT) :: PSI1(NGL)
      COMPLEX(8)       , INTENT(IN)    :: PSI2(NGL)
      integer(4)       , INTENT(IN)    :: ib2
      COMPLEX(8)                       :: CVAR   ,cvar2
!     ******************************************************************
      !CALL CG_INTERNAL_OVERLAP('OVERLAP',NGL,PSI2(:),PSI1(:),CVAR)
      CALL CG_DOT_PRODUCT(NGL,PSI1(:),OPSIold(ib2,:),CVAR)
      PSI1(:) = PSI1(:) - conjg(CVAR)*PSI2(:)
      RETURN
      END SUBROUTINE
!
!     ..................................................................
      SUBROUTINE CG_INTERNAL_NORMALIZE(NGL,ndim,npro,PSI,proj)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)       , INTENT(IN)    :: NGL 
      INTEGER(4)       , INTENT(IN)    :: ndim
      INTEGER(4)       , INTENT(IN)    :: npro
      COMPLEX(8)       , INTENT(INOUT) :: PSI(NGL)
      COMPLEX(8)       , INTENT(INOUT) :: proj(ndim,npro)
      COMPLEX(8)                       :: CVAR
!     ******************************************************************
      CALL CG_INTERNAL_OVERLAP('OVERLAP',NGL,ndim,npro,PSI(:),PSI(:),proj,CVAR)
!print*,'par CG_INTERNAL_NORMALIZE',cvar
!write(*,"('norm proj 1',3f10.5)") real(proj(1,1:3))
      PSI(:) = PSI(:) / SQRT(CVAR)
      proj=proj / SQRT(CVAR)
!write(*,"('norm proj 2',3f10.5)") real(proj(1,1:3))
      RETURN
      end subroutine CG_INTERNAL_NORMALIZE




