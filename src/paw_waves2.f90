!***********************************************************************
!***********************************************************************
!**                                                                   **
!** continuation of waves object: see paw_waves1.f90 for header       **
!**                                                                   **
!***********************************************************************
!***********************************************************************
!     ..................................................................
      SUBROUTINE WAVES$ORTHOGONALIZE()
!     ******************************************************************
!     ** ENFORCES THE ORTHONORMALITY CONDITION OF THE WAVE FUNCTIONS  **
!     ******************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      COMPLEX(8),ALLOCATABLE :: OPROJ(:,:,:)
      REAL(8)   ,ALLOCATABLE :: MARR(:)
      REAL(8)   ,ALLOCATABLE :: R0(:,:)      !CURRENT ATOMIC POSITIONS
      REAL(8)   ,ALLOCATABLE :: RP(:,:)      ! NEXT ATOMIC POSITIONS
      REAL(8)   ,ALLOCATABLE :: RM(:,:)      ! LAST ATOMIC POSITIONS
      COMPLEX(8)             :: CSUM
      COMPLEX(8),ALLOCATABLE :: MAT(:,:)
      COMPLEX(8),ALLOCATABLE :: OMAT(:,:)
      COMPLEX(8),ALLOCATABLE :: OOMAT(:,:)
      COMPLEX(8),ALLOCATABLE :: AUXMAT(:,:)
      COMPLEX(8),ALLOCATABLE :: LAMBDA(:,:)
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DO(:,:)      ! 1-CENTER OVERLAP
      REAL(8)   ,ALLOCATABLE :: G2(:)        ! SQUARE OF REC. LATTICE VECTORS
      REAL(8)   ,ALLOCATABLE :: GVEC(:,:)    ! RECIPROCAL LATTICE VECTORS
      REAL(8)                :: OCCI,OCCJ
      REAL(8)                :: FAC
      INTEGER(4)             :: NBX
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: IKPT,ISPIN,IPRO,IAT,ISP,M
      INTEGER(4)             :: LN1,L1,M1,LMN1,IPRO1
      INTEGER(4)             :: LN2,L2,M2,LMN2,IPRO2
      INTEGER(4)             :: IBH,IB,IDIM,IG,I,J
      INTEGER(4)             :: NGL,NBH,NB,LMNX,LNX,LMX,LN
      INTEGER(4)             :: NAT,IND
      REAL(8)   ,PARAMETER   :: DSMALL=1.D-12
      REAL(8)                :: SVAR,DOVER1
      REAL(8)                :: RBAS(3,3)      ! UNIT CELL
      REAL(8)                :: GBAS(3,3)      ! RECIPROCAL UNIT CELL
      REAL(8)                :: CELLVOL        ! UNIT CELL VOLUME
      REAL(8)                :: MAPTOCELL(3,3) ! R(+)=MAPTOCELL*X(+)
      REAL(8)                :: CELLSCALE      ! PSI(+)=CELLSCALE*PSI(+)
      REAL(8)   ,ALLOCATABLE :: YLM(:)
      LOGICAL(4)             :: TSTRESS
      LOGICAL(4)             :: TINV
      LOGICAL(4),PARAMETER   :: TTEST=.FALSE.
      COMPLEX(8)             :: CSVAR
      REAL(8)   ,ALLOCATABLE :: NORM(:)
      REAL(8)   ,ALLOCATABLE :: RMAT(:,:),ROMAT(:,:),ROOMAT(:,:),RLAMBDA(:,:)
      INTEGER(4),ALLOCATABLE :: SMAP(:)
      INTEGER(4)             :: I1,J1,K
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$ORTHOGONALIZE')
                             CALL TIMING$CLOCKON('WAVES$ORTHOGONALIZE')
      NPRO=MAP%NPRO
      NAT=MAP%NAT
      CALL CELL$GETL4('MOVE',TSTRESS)
!
!     ==================================================================
!     == COLLECT OCCUPATIONS                                          ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPT,NSPIN))
      CALL DYNOCC$GETR8A('OCC',NBX*NKPT*NSPIN,OCC)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==================================================================
!     ==  CALCULATE FORCE OF CONSTRAINT                               ==
!     ==================================================================
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NBH=THIS%NBH
          NB=THIS%NB
!
!         ==============================================================
!         ==  EVALUATE  DO<P|PSI>  ASSUMING <PRO|PSI> IS STILL VALID  ==
!         ==============================================================
          ALLOCATE(OPROJ(NDIM,NBH,NPRO))
          OPROJ(:,:,:)=(0.D0,0.D0)
          IPRO=1
          DO IAT=1,NAT
            ISP=MAP%ISP(IAT)
            LNX=MAP%LNX(ISP)
            LMNX=MAP%LMNX(ISP)
            ALLOCATE(DO(LNX,LNX))
            CALL SETUP$1COVERLAP(ISP,LNX,DO)
            CALL WAVES_OPROJ(LNX,MAP%LOX(1:lnx,ISP),DO,NDIM,LMNX,NBH &
      &          ,THIS%PROJ(:,:,ipro:ipro+lmnx-1),OPROJ(:,:,ipro:ipro+lmnx-1))
            DEALLOCATE(DO)
            IPRO=IPRO+LMNX
          ENDDO
!
!         ==============================================================
!         ==  ADD  |PSI>+|P>DO<P|PSI>                                 ==
!         ==============================================================
          ALLOCATE(THIS%OPSI(NGL,NDIM,NBH))
          THIS%OPSI(:,:,:)=THIS%PSI0(:,:,:)
          CALL WAVES_ADDPRO(MAP,GSET,NAT,R0,NGL,NDIM,NBH,NPRO,THIS%OPSI,OPROJ)
          DEALLOCATE(OPROJ)
!
!         ==============================================================
!         ==  DIVIDE BY WAVE FUNCTION MASS                            ==
!         ==============================================================
          ALLOCATE(MARR(NGL))
          CALL PLANEWAVE$GETR8A('G2',NGL,MARR)
          IF(ASSOCIATED(GSET%DMPSI)) THEN
            DO IG=1,NGL
!             SVAR=1.D0+ANNEE+GSET%DMPSI(IG)
              SVAR=1.D0+ANNEE 
              MARR(IG)=DELT**2/(SVAR*GSET%MPSI(IG))
            ENDDO
          ELSE
            SVAR=DELT**2/(1.D0+ANNEE)
            DO IG=1,NGL
!MARR(IG)=EMASS*(1.D0+EMASSCG2*MARR(IG))  !OLD VERSION
!MARR(IG)=DELT**2/(MARR(IG)*(1.D0+ANNEE))
              MARR(IG)=SVAR/GSET%MPSI(IG)
            ENDDO
          END IF
          DO IB=1,NBH
            DO IDIM=1,NDIM
              DO IG=1,NGL
                THIS%OPSI(IG,IDIM,IB)=MARR(IG)*THIS%OPSI(IG,IDIM,IB)
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(MARR)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  UPDATE STRUCTUREFACTOR, ETC                                 ==
!     ==================================================================
      ALLOCATE(RP(3,NAT))
      ALLOCATE(RM(3,NAT))
      CALL ATOMLIST$GETR8A('R(+)',0,3*NAT,RP)
      IF(TSTRESS) THEN
!       == PREDICT NEW POSITIONS =======================================
        CALL CELL$GETR8A('TP',9,RBAS)
        CALL CELL$GETR8A('MAPTOCELL',9,MAPTOCELL)
        CALL ATOMLIST$GETR8A('R(-)',0,3*NAT,RM)
        DO IAT=1,NAT
          RP(:,IAT)=RP(:,IAT)-RM(:,IAT)+MATMUL(MAPTOCELL,RM(:,IAT))
        ENDDO
!       ==  SCALING FACTOR FOR WAVE FUNCTIONS ==========================
        CALL GBASS(RBAS,GBAS,CELLVOL)
        CELLSCALE=MAPTOCELL(1,1)*(MAPTOCELL(2,2)*MAPTOCELL(3,3)  &
     &                           -MAPTOCELL(3,2)*MAPTOCELL(2,3)) &
     &           +MAPTOCELL(2,1)*(MAPTOCELL(3,2)*MAPTOCELL(1,3)  &
     &                           -MAPTOCELL(1,2)*MAPTOCELL(3,3)) &
     &           +MAPTOCELL(3,1)*(MAPTOCELL(1,2)*MAPTOCELL(2,3)  &
     &                           -MAPTOCELL(2,2)*MAPTOCELL(1,3))
        CELLSCALE=1.D0/SQRT(CELLSCALE)
!       == NOW K-POINT DEPENDENT DATA ==================================
        DO IKPT=1,NKPT
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
!
!         ==  UPDATE PROJECTOR FUNCTIONS ===============================
          ALLOCATE(G2(NGL))
          CALL PLANEWAVE$GETR8A('G2',NGL,G2)
          IND=0
          DO ISP=1,MAP%NSP
            CALL SETUP$ISELECT(ISP)
            DO LN=1,MAP%LNX(ISP)
              IND=IND+1
              CALL SETUP$GETFOFG('PRO',.FALSE.,LN,NGL,G2,CELLVOL,GSET%PRO(1,IND))
              CALL SETUP$GETFOFG('PRO',.TRUE.,LN,NGL,G2,CELLVOL,GSET%DPRO(1,IND))
            ENDDO
          ENDDO
          DEALLOCATE(G2)
!
!         ==  UPDATE SPHERICAL HARMONICS  ==============================
          ALLOCATE(GVEC(3,NGL))
          CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
          LMX=MAP%LMX
          ALLOCATE(YLM(LMX))
          DO IG=1,NGL
            CALL GETYLM(LMX,GVEC(1,IG),YLM)
            GSET%YLM(IG,:)=YLM(:) 
          ENDDO        
          DEALLOCATE(YLM)
!
!         == NOW THE STRAINED SPHERICAL HARMONICS ======================
          CALL WAVES_STRAINEDYLM(NGL,LMX,GVEC,GSET%YLM,GSET%SYLM)
          DEALLOCATE(GVEC)
!
!         == RESCALE WAVE FUNCTIONS ====================================
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            THIS%PSIM(:,:,:)=THIS%PSIM(:,:,:)*CELLSCALE
            THIS%OPSI(:,:,:)=THIS%OPSI(:,:,:)*CELLSCALE 
          ENDDO
        ENDDO
      ELSE
        CELLSCALE=1.D0
      END IF
!
!     ==================================================================
!     ==  NOW ORTHOGONALIZE                                           ==
!     ==================================================================
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NBH=THIS%NBH
          NB=THIS%NB
!
!         ==============================================================
!         ==  CALCULATE PROJECTIONS FOR THE NEW POSITIONS             ==
!         ==============================================================
          CALL WAVES_PROJECTIONS(MAP,GSET,NAT,RP,NGL,NDIM,NBH,NPRO &
     &                                             ,THIS%PSIM,THIS%PROJ)
          ALLOCATE(OPROJ(NDIM,NBH,NPRO))
          CALL WAVES_PROJECTIONS(MAP,GSET,NAT,RP,NGL,NDIM,NBH,NPRO &
     &                                                 ,THIS%OPSI,OPROJ)
!
!         ==============================================================
!         ==  1C-OVERLAP OF <PSI0|PSI0>, <OPSI|PSI0> AND <OPSI|OPSI>  ==
!         ==============================================================
          ALLOCATE(MAT(NB,NB))
          CALL WAVES_1COVERLAP(.TRUE.,MAP,NDIM,NBH,NB,NPRO &
         &                    ,THIS%PROJ,THIS%PROJ,MAT)
          ALLOCATE(OMAT(NB,NB))
          CALL WAVES_1COVERLAP(.FALSE.,MAP,NDIM,NBH,NB,NPRO &
       &                      ,OPROJ,THIS%PROJ,OMAT)
          ALLOCATE(OOMAT(NB,NB))
          CALL WAVES_1COVERLAP(.TRUE.,MAP,NDIM,NBH,NB,NPRO &
       &                      ,OPROJ,OPROJ,OOMAT)
!
!         ==============================================================
!         ==  NOW ADD OVERLAP OF PSEUDO WAVE FUNCTIONS                ==
!         ==============================================================
          ALLOCATE(AUXMAT(NB,NB))
          CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB &
      &                     ,THIS%PSIM,THIS%PSIM,AUXMAT)
          DO I=1,NB
            DO J=1,NB
              MAT(I,J)=MAT(I,J)+AUXMAT(I,J)
            ENDDO
          ENDDO
          CALL WAVES_OVERLAP(.FALSE.,NGL,NDIM,NBH,NB &
      &                     ,THIS%OPSI,THIS%PSIM,AUXMAT)
          DO I=1,NB
            DO J=1,NB
              OMAT(I,J)=OMAT(I,J)+AUXMAT(I,J)
            ENDDO
          ENDDO
          CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB &
       &                    ,THIS%OPSI,THIS%OPSI,AUXMAT)
          DO I=1,NB
            DO J=1,NB
              OOMAT(I,J)=OOMAT(I,J)+AUXMAT(I,J)
            ENDDO
          ENDDO
          DEALLOCATE(AUXMAT)
!
!         ==================================================================
!         ==  CALCULATE LAGRANGE PARAMETERS                               ==
!         ==================================================================
          ALLOCATE(LAMBDA(NB,NB))
          LAMBDA(:,:)=THIS%RLAM0(:,:)
          DO I=1,NB
            DO J=I+1,NB
              CSVAR=0.5D0*(LAMBDA(I,J)+CONJG(LAMBDA(J,I)))
              LAMBDA(I,J)=CSVAR
              LAMBDA(J,I)=CONJG(CSVAR)
            ENDDO
          ENDDO
!===============================================================          
          IF(.NOT.TSAFEORTHO) THEN
            ALLOCATE(SMAP(NB))
            DO I=1,NB
              SMAP(I)=I
            ENDDO
            IF(TSWAPSTATES) THEN
              SVAR=0.D0
              DO I=1,NB
                DO J=I+1,NB
                  SVAR=MAX(SVAR,ABS(LAMBDA(I,J)))
                  SVAR=MAX(SVAR,ABS(LAMBDA(J,I)))
                  I1=SMAP(I)
                  J1=SMAP(J)
                  IF(REAL(LAMBDA(J1,J1)).LT.REAL(LAMBDA(I1,I1))) THEN
                    K=SMAP(I)
                    SMAP(I)=SMAP(J)
                    SMAP(J)=K
                  END IF
                ENDDO
              ENDDO  
!PRINT*,'SMAP ',SMAP
!PRINT*,'LAMBDA ',(REAL(LAMBDA(I,I))*27.211D0,I=1,NB)
!PRINT*,'MAX ',SVAR*27.211D0
              DO I=1,NB-1
                I1=SMAP(I)
                J1=SMAP(I+1)
                IF(REAL(LAMBDA(I1,I1)).GT.REAL(LAMBDA(J1,J1))) THEN
                  CALL ERROR$MSG('STATE ORDERING FAILED')
                  CALL ERROR$STOP('WAVES$ORTHOGONALIZE')
                END IF
              ENDDO
!!$PRINT*,'---------- LAMBDA SPIN: ',ISPIN,' ------------------'
!!$WRITE(*,FMT='("LAMBDA SMAP:",20I4)')SMAP
!!$WRITE(*,FMT='("LAMBDA IMAP:",20I4)') (I,I=1,NB)
!!$PRINT*,'MAX NON-DIAGNAL LAMBDA',SVAR*27.211D0
!!$WRITE(*,FMT='("LAMBDA BEFORE:",10F8.3)') (REAL(LAMBDA(I,I))*27.211D0,I=1,NB)
!!$WRITE(*,FMT='("LAMBDA OCC:   ",10F8.3)')OCC(:,IKPT,ISPIN)
!!$ELSE
!!$   PRINT*,'LAMBDA SMAP SWITCHED OFF!!'
            END IF
          END IF
!===============================================================          
          DO I=1,NB
            DO J=1,NB
              OCCI=OCC(I,IKPT,ISPIN)
              OCCJ=OCC(J,IKPT,ISPIN)
              IF(OCCI+OCCJ.LT.DSMALL)OCCJ=DSMALL
              LAMBDA(I,J)=LAMBDA(I,J)*0.5D0*OCCI/(OCCI+OCCJ)
            ENDDO
          ENDDO
!
          IF(TSAFEORTHO) THEN
            CALL PLANEWAVE$GETL4('TINV',TINV)
            IF(TINV) THEN
              ALLOCATE(RMAT(NB,NB))
              RMAT=REAL(MAT)
              ALLOCATE(ROMAT(NB,NB))
              ROMAT=REAL(OMAT)
              ALLOCATE(ROOMAT(NB,NB))
              ROOMAT=REAL(OOMAT)
              ALLOCATE(RLAMBDA(NB,NB))
              RLAMBDA=REAL(LAMBDA)
              CALL WAVES_ORTHO_X(NB,OCC(1,IKPT,ISPIN) &
       &                        ,ROOMAT,RMAT,ROMAT,RLAMBDA)
              LAMBDA=CMPLX(RLAMBDA,0.D0,8)
              DEALLOCATE(RLAMBDA)
              DEALLOCATE(RMAT)
              DEALLOCATE(ROMAT)
              DEALLOCATE(ROOMAT)
            ELSE
              CALL WAVES_ORTHO_X_C(NB,OCC(1,IKPT,ISPIN),OOMAT,MAT,OMAT,LAMBDA)
            END IF
          ELSE
            CALL WAVES_ORTHO_Y_C(NB,MAT,OMAT,OOMAT,LAMBDA,SMAP)
          END IF
          DEALLOCATE(MAT)
          DEALLOCATE(OMAT)
          DEALLOCATE(OOMAT)
          IF(.NOT.TSAFEORTHO)DEALLOCATE(SMAP)
!
!         ==================================================================
!         ==  CALCULATE |PSI(+)>=|PSI>+|CHI>LAMBDA                        ==
!         ==================================================================
          CALL WAVES_ADDOPSI(NGL,NDIM,NBH,NB,THIS%PSIM,THIS%OPSI,LAMBDA)
          DEALLOCATE(THIS%OPSI)
PRINT*,'WARNING FROM WAVES$ORTHOGONALIZE:'
PRINT*,'MAKE SURE THAT PDOS AND GRAPHICS PICK UP A CONSISTENT SET OF '
PRINT*,'WAVE FUNCTIONS  AND PROJECTOR FUNCTIONS'
          CALL WAVES_ADDOPROJ(NPRO,NDIM,NBH,NB,THIS%PROJ,OPROJ,LAMBDA)
          DEALLOCATE(OPROJ)
!
!         ==================================================================
!         ==  RESCALE GAMMA                                               ==
!         ==================================================================
!         == MASS IS NOT INCLUDED HERE (MASS TENSOR IS TAKEN CARE OF IN 
!         == PSIBAR AND (1/M)O|PSI> 
          DO I=1,NB
            DO J=I,NB
              LAMBDA(I,J)=0.5D0*(LAMBDA(I,J)+CONJG(LAMBDA(J,I)))
              LAMBDA(J,I)=CONJG(LAMBDA(I,J))
            ENDDO
          ENDDO
!PRINT*,'LAMBDA[EV] ',(REAL(LAMBDA(I,I))*27.211D0,I=1,NB)
          IF(.NOT.ASSOCIATED(THIS%RLAM0))ALLOCATE(THIS%RLAM0(NB,NB))
          THIS%RLAM0(:,:)=LAMBDA(:,:)
          DEALLOCATE(LAMBDA)
!
!         ==================================================================
!         ==  TEST ORTHONORMALITY                                         ==
!         ==================================================================
          IF(TTEST) THEN
            ALLOCATE(OPROJ(NDIM,NBH,NPRO))
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,RP,NGL,NDIM,NBH,NPRO,THIS%PSIM,OPROJ)
            ALLOCATE(AUXMAT(NB,NB))
            CALL WAVES_1COVERLAP(.TRUE.,MAP,NDIM,NBH,NB,NPRO,OPROJ,OPROJ,AUXMAT)
            DEALLOCATE(OPROJ)
            CSUM=(0.D0,0.D0)
            DO I=1,NB
              CSUM=CSUM+AUXMAT(I,I)*OCC(I,IKPT,ISPIN)
            ENDDO
!PRINT*,'1C-CHARGE AFTER ORTHOGONALIZATION ',CSUM
            ALLOCATE(MAT(NB,NB))
            CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB,THIS%PSIM,THIS%PSIM,MAT)
            CSUM=(0.D0,0.D0)
            DO I=1,NB
              CSUM=CSUM+MAT(I,I)*OCC(I,IKPT,ISPIN)
            ENDDO
!PRINT*,'PS-CHARGE AFTER ORTHOGONALIZATION ',CSUM
            DO J=1,NB
              DO I=1,NB
                MAT(I,J)=MAT(I,J)+AUXMAT(I,J)
              ENDDO
            ENDDO
            ALLOCATE(NORM(NB))
            CALL LIB$DIAGC8(NB,MAT,NORM,AUXMAT)
!CALL CDIAG(NB,NB,MAT,NORM,AUXMAT)
            DEALLOCATE(MAT)
            DEALLOCATE(AUXMAT)
!PRINT*,'NORM ',NORM
            DO I=1,NB
              IF(ABS(NORM(I)-1.D0).GT.1.D-4) THEN
                CALL ERROR$MSG('ORTHOGONALIZATION FAILED')
                CALL ERROR$I4VAL('STATE',I)
                CALL ERROR$R8VAL('NORM',NORM(I))
                CALL ERROR$STOP('WAVES$ORTHOGONALIZE')
              END IF
            ENDDO
            DEALLOCATE(NORM)
          END IF
        ENDDO
      ENDDO
      DEALLOCATE(OCC)
!
!     ==================================================================
!     ==  NOW TRANSFORM MACK INTO ORIGINAL CELL                       ==
!     ==================================================================
      IF(TSTRESS) THEN
        CALL CELL$GETR8A('T0',9,RBAS)
        DO IKPT=1,NKPT
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            THIS%PSIM(:,:,:)=THIS%PSIM(:,:,:)/CELLSCALE
          ENDDO   
        ENDDO
      END IF
      DEALLOCATE(RP)
      DEALLOCATE(R0)
      DEALLOCATE(RM)
!
!     ==================================================================
!     ==  CALCULATE SECOND PART OF WAVE FUNCTION KINETIC ENERGY       ==
!     ==================================================================
      CALL WAVES_WAVEKINETIC(WAVEEKIN2)
                             CALL TIMING$CLOCKOFF('WAVES$ORTHOGONALIZE')
                                    CALL TRACE$POP
      RETURN
      END
!
!      .................................................................
       SUBROUTINE WAVES_OPROJ(LNX,LOX,DO,NDIM,LMNX,NB,PROJ,OPROJ)
!      *****************************************************************
!      **                                                             **
!      **  DO<PRO|PSI>                                                **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: LNX
       INTEGER(4),INTENT(IN) :: LOX(LNX)
       REAL(8)   ,INTENT(IN) :: DO(LNX,LNX)
       INTEGER(4),INTENT(IN) :: LMNX
       INTEGER(4),INTENT(IN) :: NDIM
       INTEGER(4),INTENT(IN) :: NB
       COMPLEX(8),INTENT(IN) :: PROJ(NDIM,NB,LMNX)
       COMPLEX(8),INTENT(OUT):: OPROJ(NDIM,NB,LMNX)
       INTEGER(4)            :: LMN1,LMN2
       INTEGER(4)            :: LMN10,LMN20
       INTEGER(4)            :: LN1,LN2
       INTEGER(4)            :: L1,L2
       INTEGER(4)            :: IB,IDIM,M
       REAL(8)               :: DOVER1
!      *****************************************************************
       OPROJ(:,:,:)=(0.D0,0.D0)
       LMN10=0
       DO LN1=1,LNX
         L1=LOX(LN1)
         LMN20=0
         DO LN2=1,LNX
           L2=LOX(LN2)
           IF(L1.EQ.L2) THEN
             DOVER1=DO(LN1,LN2)
             DO M=1,2*L1+1
               LMN1=LMN10+M
               LMN2=LMN20+M
               DO IDIM=1,NDIM
                 DO IB=1,NB
                   OPROJ(IDIM,IB,LMN1)=OPROJ(IDIM,IB,LMN1) &
     &                          +DOVER1*PROJ(IDIM,IB,LMN2)
                 ENDDO
               ENDDO
             END DO
           END IF
           LMN20=LMN20+2*L2+1
         ENDDO
         LMN10=LMN10+2*L1+1
       ENDDO
       RETURN
       END
!
!      .................................................................
       SUBROUTINE WAVES_ADDOPSI(NGL,NDIM,NBH,NB,PSIBAR,OPSI,LAMBDA)
!      *****************************************************************
!      **                                                             **
!      **  |PSI(+)>=|PSIBAR>+O|PSI(0)>*LAMBDA                         **
!      **                                                             **
!      **                                                             **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)   :: NGL
       INTEGER(4),INTENT(IN)   :: NDIM
       INTEGER(4),INTENT(IN)   :: NBH
       INTEGER(4),INTENT(IN)   :: NB
       COMPLEX(8),INTENT(INOUT):: PSIBAR(NGL*NDIM,NBH)
       COMPLEX(8),INTENT(INOUT):: OPSI(NGL*NDIM,NBH)
       COMPLEX(8),INTENT(IN)   :: LAMBDA(NB,NB)
       LOGICAL(4),PARAMETER    :: TESSL=.TRUE.
       LOGICAL(4)              :: TINV
       INTEGER(4)              :: IBH1,IBH2,I,IDIM,IBH
       INTEGER(4)              :: IB1A,IB1B,IB2A,IB2B
       COMPLEX(8),ALLOCATABLE  :: LAMBDA1(:,:)
       COMPLEX(8),ALLOCATABLE  :: LAMBDA2(:,:)
       COMPLEX(8),ALLOCATABLE  :: TPSI(:)
       COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
       COMPLEX(8)              :: CSVAR,CSVAR1,CSVAR2
       INTEGER(4)              :: NGLNDIM
!      *****************************************************************
                               CALL TIMING$CLOCKON('WAVES_ADDOPSI')
       TINV=NBH.NE.NB
       NGLNDIM=NGL*NDIM
       IF(.NOT.TINV) THEN
         CALL LIB$ADDPRODUCTC8(.FALSE.,NGLNDIM,NBH,NBH,OPSI,LAMBDA,PSIBAR)
!        IF(TESSL) THEN
!          CALL ZGEMM('N','N',NGLNDIM,NBH,NBH,(1.D0,0.D0),OPSI,NGLNDIM &
!                     ,LAMBDA,NBH,(1.D0,0.D0),PSIBAR,NGLNDIM)
!        ELSE
!          DO IBH1=1,NBH
!            DO IBH2=1,NBH
!              CSVAR=LAMBDA(IBH2,IBH1)
!              DO I=1,NGLNDIM
!                PSIBAR(I,IBH1)=PSIBAR(I,IBH1)+OPSI(I,IBH2)*CSVAR
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDIF
       ELSE
         ALLOCATE(LAMBDA1(NBH,NBH))
         ALLOCATE(LAMBDA2(NBH,NBH))
         DO IBH1=1,NBH
           IB1A=2*IBH1-1
           IB1B=2*IBH1
           DO IBH2=1,NBH
             IB2A=2*IBH2-1
             IB2B=2*IBH2
             CSVAR1=   LAMBDA(IB1A,IB2A)+CI*LAMBDA(IB1A,IB2B)
             CSVAR2=CI*LAMBDA(IB1B,IB2A)-   LAMBDA(IB1B,IB2B)
             LAMBDA1(IBH1,IBH2)=0.5D0*(CSVAR1-CSVAR2)
             LAMBDA2(IBH1,IBH2)=0.5D0*(CSVAR1+CSVAR2)
           ENDDO
         ENDDO
!        == ADD O|PSI_+>LAMBDA1 =======================================
         CALL LIB$ADDPRODUCTC8(.FALSE.,NGLNDIM,NBH,NBH,OPSI,LAMBDA1,PSIBAR)
!        IF(TESSL) THEN
!          CALL ZGEMM('N','N',NGLNDIM,NBH,NBH,(1.D0,0.D0),OPSI,NGLNDIM &
!                     ,LAMBDA1,NBH,(1.D0,0.D0),PSIBAR,NGLNDIM)
!        ELSE 
!          DO IBH1=1,NBH
!            DO IBH2=1,NBH
!              CSVAR=LAMBDA1(IBH2,IBH1)
!              DO I=1,NGLNDIM
!                PSIBAR(I,IBH1)=PSIBAR(I,IBH1)+OPSI(I,IBH2)*CSVAR
!              ENDDO
!            ENDDO
!          ENDDO
!        END IF
         DEALLOCATE(LAMBDA1)
!        == INVERT OPSI ================================================
         ALLOCATE(TPSI(NGLNDIM))
         DO IBH1=1,NBH
           I=1
           DO IDIM=1,NDIM
             CALL PLANEWAVE$INVERTG(NGL,OPSI(I,IBH1),TPSI(I))
             I=I+NGL
           ENDDO
           OPSI(:,IBH1)=TPSI(:)
         ENDDO
         DEALLOCATE(TPSI)
!
!        == ADD O|PSI_+>LAMBDA1 =======================================
         CALL LIB$ADDPRODUCTC8(.FALSE.,NGLNDIM,NBH,NBH,OPSI,LAMBDA2,PSIBAR)
!        IF(TESSL) THEN
!          CALL ZGEMM('N','N',NGLNDIM,NBH,NBH,(1.D0,0.D0),OPSI,NGLNDIM &
!                     ,LAMBDA2,NBH,(1.D0,0.D0),PSIBAR,NGLNDIM)
!        ELSE 
!          DO IBH1=1,NBH
!            DO IBH2=1,NBH
!              CSVAR=LAMBDA2(IBH2,IBH1)
!              DO I=1,NGLNDIM
!                PSIBAR(I,IBH1)=PSIBAR(I,IBH1)+OPSI(I,IBH2)*CSVAR
!              ENDDO
!            ENDDO
!          ENDDO
!        END IF
         DEALLOCATE(LAMBDA2)
       END IF
                               CALL TIMING$CLOCKOFF('WAVES_ADDOPSI')
       RETURN
       END
!
!      .................................................................
       SUBROUTINE WAVES_ADDOPROJ(NPRO,NDIM,NBH,NB,PROJ,OPROJ,LAMBDA)
!      *****************************************************************
!      **                                                             **
!      **  |PSI(+)>=|PSIBAR>+O|PSI(0)>*LAMBDA                         **
!      **                                                             **
!      **                                                             **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)   :: NDIM
       INTEGER(4),INTENT(IN)   :: NBH
       INTEGER(4),INTENT(IN)   :: NB
       INTEGER(4),INTENT(IN)   :: NPRO
       COMPLEX(8),INTENT(INOUT):: PROJ(NDIM,NBH,NPRO)
       COMPLEX(8),INTENT(IN)   :: OPROJ(NDIM,NBH,NPRO)
       COMPLEX(8),INTENT(IN)   :: LAMBDA(NB,NB)
       LOGICAL(4),PARAMETER    :: TESSL=.TRUE.
       LOGICAL(4)              :: TINV
       INTEGER(4)              :: IBH1,IBH2,I,IDIM,IBH
       INTEGER(4)              :: IB1A,IB1B,IB2A,IB2B
       COMPLEX(8),ALLOCATABLE  :: LAMBDA1(:,:)
       COMPLEX(8),ALLOCATABLE  :: LAMBDA2(:,:)
       COMPLEX(8),ALLOCATABLE  :: TPSI(:)
       COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
       COMPLEX(8)              :: CSVAR,CSVAR1,CSVAR2
       INTEGER(4)              :: IPRO
!      *****************************************************************
       TINV=NBH.NE.NB
       IF(.NOT.TINV) THEN
         DO IBH1=1,NBH
           DO IBH2=1,NBH
             CSVAR=LAMBDA(IBH2,IBH1)
             DO IPRO=1,NPRO
               DO IDIM=1,NDIM
                 PROJ(IDIM,IBH1,IPRO)=PROJ(IDIM,IBH1,IPRO) &
     &                              +OPROJ(IDIM,IBH2,IPRO)*CSVAR
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ELSE
         ALLOCATE(LAMBDA1(NBH,NBH))
         ALLOCATE(LAMBDA2(NBH,NBH))
         DO IBH1=1,NBH
           IB1A=2*IBH1-1
           IB1B=2*IBH1
           DO IBH2=1,NBH
             IB2A=2*IBH2-1
             IB2B=2*IBH2
             CSVAR1=   LAMBDA(IB1A,IB2A)+CI*LAMBDA(IB1A,IB2B)
             CSVAR2=CI*LAMBDA(IB1B,IB2A)-   LAMBDA(IB1B,IB2B)
             LAMBDA1(IBH1,IBH2)=0.5D0*(CSVAR1-CSVAR2)
             LAMBDA2(IBH1,IBH2)=0.5D0*(CSVAR1+CSVAR2)
           ENDDO
         ENDDO
!        == ADD O|PSI_+>LAMBDA1 =======================================
         DO IBH1=1,NBH
           DO IBH2=1,NBH
             CSVAR1=LAMBDA1(IBH2,IBH1)
             CSVAR2=LAMBDA2(IBH2,IBH1)
             DO IPRO=1,NPRO
               DO IDIM=1,NDIM
                 PROJ(IDIM,IBH1,IPRO)=PROJ(IDIM,IBH1,IPRO) &
     &                             +OPROJ(IDIM,IBH2,IPRO) *CSVAR1 &
     &                       +CONJG(OPROJ(IDIM,IBH2,IPRO))*CSVAR2
               ENDDO
             ENDDO
           ENDDO
         ENDDO
         DEALLOCATE(LAMBDA1)
         DEALLOCATE(LAMBDA2)
       END IF
       RETURN
       END
!
!      ..............................................................
       SUBROUTINE WAVES_OVERLAP(TID,NGL,NDIM,NBH,NB,PSI1,PSI2,MAT)
!      **                                                          **
!      **  CALCULATES <PSI1|PSI2>                                  **
!      **                                                          **
       USE MPE_MODULE
       IMPLICIT NONE
       LOGICAL(4),INTENT(IN) :: TID !INDICATES THAT PSI1=PSI2
       INTEGER(4),INTENT(IN) :: NGL
       INTEGER(4),INTENT(IN) :: NDIM
       INTEGER(4),INTENT(IN) :: NBH
       INTEGER(4),INTENT(IN) :: NB
       COMPLEX(8),INTENT(IN) :: PSI1(NGL,NDIM,NBH)
       COMPLEX(8),INTENT(IN) :: PSI2(NGL,NDIM,NBH)
       COMPLEX(8),INTENT(OUT):: MAT(NB,NB)
       INTEGER(4)            :: IBH1,IBH2,IB1,IB2,IDIM,IG
       INTEGER(4)            :: IB1A,IB1B,IB2A,IB2B,I,J
       COMPLEX(8)            :: CSVARPP,CSVARPM,CSVARP,CSVARM
       COMPLEX(8)            :: CSVAR
       COMPLEX(8),ALLOCATABLE:: TMAT(:,:)
       REAL(8)               :: RE,IM
       REAL(8)               :: MAT2(2,2)
       LOGICAL(4)            :: TINV
!      **************************************************************
                             CALL TIMING$CLOCKON('WAVES_OVERLAP')
       TINV=(NB.NE.NBH)
       IF(.NOT.TINV) THEN
         IF(TID) THEN
           CALL PLANEWAVE$SCALARPRODUCT('=',NGL,NDIM,NB,PSI1,NB,PSI2,MAT)
         ELSE
           CALL PLANEWAVE$SCALARPRODUCT(' ',NGL,NDIM,NB,PSI1,NB,PSI2,MAT)
         END IF
         CALL MPE$COMBINE('+',MAT)
                             CALL TIMING$CLOCKOFF('WAVES_OVERLAP')
         RETURN
       ENDIF
!
!      =================================================================
!      == NOW DEAL WITH SUPER WAVE FUNCTIONS                          ==
!      =================================================================
!      == <PSI_+|PSI_+> ================================================
       ALLOCATE(TMAT(NBH,NBH))
       IF(TID) THEN
         CALL PLANEWAVE$SCALARPRODUCT('=',NGL,NDIM,NBH,PSI1,NBH,PSI2,TMAT)
       ELSE
         CALL PLANEWAVE$SCALARPRODUCT(' ',NGL,NDIM,NBH,PSI1,NBH,PSI2,TMAT)
       END IF
       DO IBH1=1,NBH
         IB1A=2*IBH1-1
         IB1B=MIN(2*IBH1,NB)
         DO IBH2=1,NBH
           IB2A=2*IBH2-1
           IB2B=MIN(2*IBH2,NB)
           RE=0.5D0* REAL(TMAT(IBH1,IBH2))
           IM=0.5D0*AIMAG(TMAT(IBH1,IBH2))
!          == WATCH THE ORDER OF THE FOLLOWING 4 LINES!!! =============
!          == BECAUSE IB1B=IB1A FOR 2*IBH1>NB =========================
!          == AND     IB2B=IB2A FOR 2*IBH2>NB =========================
           MAT2(1,1)=+RE
           MAT2(1,2)=+IM
           MAT2(2,1)=-IM
           MAT2(2,2)=+RE
           DO I=IB1A,IB1B
             DO J=IB2A,IB2B
               MAT(I,J)=CMPLX(MAT2(I-IB1A+1,J-IB2A+1),0.D0,8)
             ENDDO
           ENDDO
         ENDDO
       ENDDO
!      == NOW <PSI-|PSI+> ==============================================
       CALL PLANEWAVE$SCALARPRODUCT('-',NGL,NDIM,NBH,PSI1,NBH,PSI2,TMAT)
       DO IBH1=1,NBH
         IB1A=2*IBH1-1
         IB1B=MIN(2*IBH1,NB)
         DO IBH2=1,NBH
!          == <PSI-(I)|PSI+(J)> ========================================
           IB2A=2*IBH2-1
           IB2B=MIN(2*IBH2,NB)
           RE=0.5D0* REAL(TMAT(IBH1,IBH2))
           IM=0.5D0*AIMAG(TMAT(IBH1,IBH2))
           MAT2(1,1)=+RE
           MAT2(1,2)=+IM
           MAT2(2,1)=+IM
           MAT2(2,2)=-RE
           DO I=IB1A,IB1B
             DO J=IB2A,IB2B
               MAT(I,J)=MAT(I,J)+CMPLX(MAT2(I-IB1A+1,J-IB2A+1),0.D0,8)
             ENDDO
           ENDDO
         ENDDO
       ENDDO
       DEALLOCATE(TMAT)
       CALL MPE$COMBINE('+',MAT)
                             CALL TIMING$CLOCKOFF('WAVES_OVERLAP')
       RETURN
       END
!
!     ..................................................................
      SUBROUTINE WAVES_1COVERLAP(TID,MAP,NDIM,NBH,NB,NPRO,PROJ1,PROJ2,MAT)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE 1C CONTRIBUTION TO THE OVERLAP MATRIX        **
!     **  ( RESULT IS ADDED TO "MAT"! )                               **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY: MAP_TYPE
      IMPLICIT NONE
      LOGICAL(4)  ,INTENT(IN)   :: TID
      TYPE(MAP_TYPE),INTENT(IN) :: MAP
      INTEGER(4)  ,INTENT(IN)   :: NDIM
      INTEGER(4)  ,INTENT(IN)   :: NBH
      INTEGER(4)  ,INTENT(IN)   :: NB
      INTEGER(4)  ,INTENT(IN)   :: NPRO
      COMPLEX(8)  ,INTENT(IN)   :: PROJ1(NDIM,NBH,NPRO) ! <P|PSI1>
      COMPLEX(8)  ,INTENT(IN)   :: PROJ2(NDIM,NBH,NPRO) ! <P|PSI2>
      COMPLEX(8)  ,INTENT(OUT)  :: MAT(NB,NB)
      REAL(8)     ,ALLOCATABLE  :: DOVER(:,:)
      INTEGER(4)                :: NTASKS,THISTASK
      INTEGER(4)                :: ICOUNT,IAT,ISP
      INTEGER(4)                :: LMN1,LN1,L1
      INTEGER(4)                :: LMN2,LN2,L2
      INTEGER(4)                :: IB1,IB2,IPRO,IPRO1,IPRO2,M
      INTEGER(4)                :: IBH1,IBH2,IDIM,LNX
      INTEGER(4)                :: IB1A,IB1B,IB2A,IB2B,I,J
      COMPLEX(8)                :: CSVAR,CSVAR1,CSVAR2
      REAL(8)                   :: RE1,IM1,RE2,IM2
      REAL(8)                   :: RMAT(2,2)
      LOGICAL(4)                :: TINV
!     ******************************************************************
      CALL MPE$QUERY(NTASKS,THISTASK)
      MAT(:,:)=(0.D0,0.D0)
      TINV=NB.NE.NBH
      ICOUNT=0
      IPRO=0
      DO IAT=1,MAP%NAT
        ISP=MAP%ISP(IAT)
!       __ SELECTION FOR PARALLEL PROCESSING____________________________
        ICOUNT=ICOUNT+1
        IF(MOD(ICOUNT-1,NTASKS).NE.THISTASK-1) THEN
          IPRO=IPRO+MAP%LMNX(ISP)
          CYCLE
        END IF
!       __ NOW CONTINUE_________________________________________________
        LNX=MAP%LNX(ISP)
        ALLOCATE(DOVER(LNX,LNX))
        CALL SETUP$1COVERLAP(ISP,LNX,DOVER)
!
        DO IBH1=1,NBH
          DO IBH2=1,NBH
            CSVAR1=(0.D0,0.D0)
            CSVAR2=(0.D0,0.D0)
            IPRO1=IPRO
            DO LN1=1,LNX
              L1=MAP%LOX(LN1,ISP)
              IPRO2=IPRO
              DO LN2=1,LNX
                L2=MAP%LOX(LN2,ISP)
                IF(L1.EQ.L2) THEN
                  CSVAR=(0.D0,0.D0)
                  DO IDIM=1,NDIM
                    DO M=1,2*L1+1
                      CSVAR=CSVAR + CONJG(PROJ1(IDIM,IBH1,IPRO1+M)) &
     &                            *       PROJ2(IDIM,IBH2,IPRO2+M)
                    ENDDO
                  ENDDO
                  CSVAR1=CSVAR1+CSVAR*DOVER(LN1,LN2)
                  IF(TINV) THEN
                    CSVAR=(0.D0,0.D0)
                    DO IDIM=1,NDIM
                      DO M=1,2*L1+1
                        CSVAR=CSVAR + PROJ1(IDIM,IBH1,IPRO1+M) &
     &                              * PROJ2(IDIM,IBH2,IPRO2+M)
                      ENDDO
                    ENDDO
                    CSVAR2=CSVAR2+CSVAR*DOVER(LN1,LN2)
                  END IF
                END IF
                IPRO2=IPRO2+2*L2+1
              ENDDO
              IPRO1=IPRO1+2*L1+1
            ENDDO

            IF(.NOT.TINV) THEN
              MAT(IBH1,IBH2)=MAT(IBH1,IBH2)+CSVAR1
            ELSE
              RE1=  REAL(CSVAR1)
              IM1= AIMAG(CSVAR1)
              RE2=  REAL(CSVAR2)
              IM2=-AIMAG(CSVAR2)  ! COMPLEX CONJUGATE OF CSVAR2 USED
              RMAT(1,1)=0.5D0*( RE1+RE2)
              RMAT(1,2)=0.5D0*( IM1-IM2)
              RMAT(2,1)=0.5D0*(-IM1-IM2)
              RMAT(2,2)=0.5D0*( RE1-RE2)
              IB1A=2*IBH1-1
              IB1B=MIN(2*IBH1,NB)
              IB2A=2*IBH2-1
              IB2B=MIN(2*IBH2,NB)
              DO I=IB1A,IB1B
                DO J=IB2A,IB2B
                  MAT(I,J)=MAT(I,J)+CMPLX(RMAT(I-IB1A+1,J-IB2A+1),0.D0,8)
                ENDDO
              ENDDO
            END IF
          ENDDO
        ENDDO
        DEALLOCATE(DOVER)
        IPRO=IPRO+MAP%LMNX(ISP)
      ENDDO
      CALL MPE$COMBINE('+',MAT)
      RETURN
      END
!
!      .................................................................
       SUBROUTINE WAVES_ORTHO_Y_C(NB,PHIPHI,CHIPHI,CHICHI,X,MAP)
!      **                                                             **
!      **  CALCULATE LAGRANGE MULTIPLIERS FOR ORTHOGONALIZATION       **
!      **    |PHI(I)>=|PHI(I)>+SUM_J |CHI(J)>X(J,I)                   **
!      **  WITH                                                       **
!      **    X(I>J)=0                                                 **
!      **                                                             **
!      **  ATTENTION!! CHIPHI0=<CHI|O|PHI>                            **
!      **        AND   PHIPHI0=<PHI|O|PHI>                            **
!      **        CONVERSION IS DONE IN THE INITIALIZATION             **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NB
       COMPLEX(8),INTENT(IN) :: PHIPHI(NB,NB) !<PHI_I|O|PHI_J>
       COMPLEX(8),INTENT(IN) :: CHIPHI(NB,NB) !<PHI_I|O|CHI>
       COMPLEX(8),INTENT(IN) :: CHICHI(NB,NB) !<CHI_I|O|CHI_J>
       COMPLEX(8),INTENT(OUT):: X(NB,NB)      ! X(I>J)=0
       INTEGER(4),INTENT(IN) :: MAP(NB)
       COMPLEX(8)            :: A(NB,NB)
       COMPLEX(8)            :: B(NB,NB)
       COMPLEX(8)            :: C(NB,NB)
       COMPLEX(8)            :: ALPHA(NB,NB)   ! (I>J)=0
       COMPLEX(8)            :: WORK(NB,NB)    
       COMPLEX(8)            :: Z(NB)
       INTEGER(4)            :: I,J,K,L,N,M
       INTEGER(4)            :: N0
       COMPLEX(8)            :: CSVAR
       REAL(8)               :: SVAR,SVAR1,SVAR2
       REAL(8)               :: MAXDEV
       LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
       REAL(8)   ,PARAMETER  :: TOL=1.D-10
       INTEGER(4)            :: NU,NU0,IU,JU
!      *****************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_Y_C')
!
!      =================================================================
!      ==  INITIALIZE                                                 ==
!      =================================================================
       DO I=1,NB
         DO J=1,NB
           A(I,J)=PHIPHI(I,J)
           B(I,J)=CHIPHI(I,J)
           C(I,J)=CHICHI(I,J)
           X(I,J)=(0.D0,0.D0)
           ALPHA(I,J)=(0.D0,0.D0)
         ENDDO
         A(I,I)=A(I,I)-(1.D0,0.D0)
         ALPHA(I,I)=(1.D0,0.D0)
       ENDDO
!
!      =================================================================
!      ==  ORTHOGONALIZATION LOOP                                     ==
!      =================================================================
!                            CALL TRACE$PASS('BEFORE ORTHOGONALIZATION LOOP')
       DO NU=1,NB
         N=MAP(NU)
!
!        ===============================================================
!        == NORMALIZE PHI(N)                                          ==
!        == PHI(N)=PHI(N)+CHI(N)*Z(N)                                 ==
!        ===============================================================
!                            CALL TRACE$PASS('NORMALIZE')
         SVAR=CONJG(B(N,N))*B(N,N)-A(N,N)*C(N,N)-AIMAG(B(N,N))**2
         SVAR1=-REAL(B(N,N))
         IF(SVAR.GE.0.D0) THEN
           SVAR2=SQRT(SVAR)
           IF(SVAR1*SVAR2.GE.0.D0) THEN
             Z(N)=SVAR1-SVAR2
           ELSE
             Z(N)=SVAR1+SVAR2
           END IF
         ELSE
           PRINT*,'ORTHOGONALIZATION FAILED! TRYING BEST APPROXIMATION...'
           Z(N)=SVAR1
         END IF
         Z(N)=Z(N)/REAL(C(N,N))
!
!        ===============================================================
!        == NOW UPDATE MATRICES                                       ==
!        ===============================================================
         NU0=NU   !SET N0=N FOR FAST CALCULATION AND N0=1 FOR TESTS
!N0=1
         DO I=1,NB
!          == A(N,M)+B(N,N)*DELTA(M)=0  ======================
           X(I,N)=X(I,N)+ALPHA(I,N)*Z(N)
         ENDDO           
         DO J=1,NB
           A(N,J)=A(N,J)+CONJG(Z(N))*B(N,J)
         ENDDO
         DO I=1,NB
           A(I,N)=A(I,N)+CONJG(B(N,I))*Z(N) 
         ENDDO
         A(N,N)=A(N,N)+CONJG(Z(N))*C(N,N)*Z(N)
         DO I=1,NB
           B(I,N)=B(I,N)+C(I,N)*Z(N)
         ENDDO
         IF(TTEST) THEN
           CALL TESTA(N,N,CSVAR)
           IF(ABS(CSVAR).GT.TOL) THEN
             WRITE(*,FMT='("NORMALIZATION OF PHI(",I4,")")')N
             PRINT*,N,N,CSVAR
           END IF
         END IF
!
!        ===============================================================
!        == ORTHOGONALIZE HIGHER PHI'S TO THIS PHI                    ==
!        == PHI(J)=PHI(J)+CHI(N)*Z(J)       J>N                       ==
!        ===============================================================
         DO IU=1,NU
           I=MAP(IU)
           Z(I)=(0.D0,0.D0)
         ENDDO
         DO IU=NU+1,NB
           I=MAP(IU)
           Z(I)=-CONJG(A(I,N)/B(N,N))
         ENDDO
!               CALL TRACE$PASS('ORTHOGONALIZE HIGHER PHIS TO THIS PHI')
!
!        ===============================================================
!        == NOW UPDATE MATRICES                                       ==
!        ===============================================================
         NU0=NU+1   !SET N0=N FOR FAST CALCULATION AND N0=1 FOR TESTS
!N0=1
         DO I=1,NB
!          == A(N,M)+B(N,N)*DELTA(M)=0  ======================
           DO JU=NU0,NB
             J=MAP(JU)
             X(I,J)=X(I,J)+ALPHA(I,N)*Z(J)
           ENDDO
         ENDDO           
         DO IU=NU0,NB
           I=MAP(IU)
           DO J=1,NB
             A(I,J)=A(I,J)+CONJG(Z(I))*B(N,J)
           ENDDO
         ENDDO
         DO I=1,NB
           DO JU=NU0,NB
             J=MAP(JU)
             A(I,J)=A(I,J)+CONJG(B(N,I))*Z(J) 
           ENDDO
         ENDDO
         DO IU=NU0,NB
           I=MAP(IU)
           DO JU=NU0,NB
             J=MAP(JU)
             A(I,J)=A(I,J)+CONJG(Z(I))*C(N,N)*Z(J)
           ENDDO
         ENDDO
         DO I=1,NB
           DO JU=NU0,NB
             J=MAP(JU)
             B(I,J)=B(I,J)+C(I,N)*Z(J)
           ENDDO
         ENDDO
         IF(TTEST) THEN
           DO IU=1,NU
             I=MAP(IU)
             DO JU=NU,NB
               J=MAP(JU)
               CALL TESTA(I,J,CSVAR)
               IF(ABS(CSVAR).GT.TOL) THEN
                 WRITE(*,FMT='("HIGHER PHIS ORTHOGONALIZED TO PHI(",I4,")")')N
                 PRINT*,I,J,CSVAR
               END IF
             ENDDO
           ENDDO
         END IF
!
!        ===============================================================
!        == ORTHOGONALIZE HIGHER CHI'S TO THIS PHI                    ==
!        == CHI(M)=CHI(M)+CHI(N)*DELTA(M)   M>N                       ==
!        ===============================================================
!               CALL TRACE$PASS('ORTHOGONALIZE HIGHER CHIS TO THIS PHI')
         DO IU=1,NU
           I=MAP(IU)
           Z(I)=(0.D0,0.D0)
         ENDDO
         DO IU=NU+1,NB
           I=MAP(IU)
!          == |CHI(J)>=|CHI(J)>+|CHI(N)>*Z(J) ==========================
!          == B(M,N)+B(N,N)*DELTA(M)=0
           Z(I)=-CONJG(B(I,N)/B(N,N))
         ENDDO
         NU0=NU+1   !SET N0=N+1 FOR FAST CALCULATION AND N0=1 FOR TESTS
!N0=1
         DO I=1,NB
           DO JU=NU0,NB
             J=MAP(JU)
             ALPHA(I,J)=ALPHA(I,J)+ALPHA(I,N)*Z(J)
           ENDDO
         ENDDO
         DO IU=NU0,NB
           I=MAP(IU)
           DO J=1,NB
             B(I,J)=B(I,J)+CONJG(Z(I))*B(N,J)
           ENDDO 
         ENDDO
         WORK(:,:)=0.D0
         DO IU=NU0,NB
           I=MAP(IU)
           DO J=1,NB
             WORK(I,J)=WORK(I,J)+CONJG(Z(I))*C(N,J) 
           ENDDO
         ENDDO
         DO I=1,NB
           DO JU=NU0,NB
             J=MAP(JU)
             WORK(I,J)=WORK(I,J)+C(I,N)*Z(J)
           ENDDO
         ENDDO
         DO IU=NU0,NB
           I=MAP(IU)
           DO JU=NU0,NB
             J=MAP(JU)
             WORK(I,J)=WORK(I,J)+CONJG(Z(I))*C(N,N)*Z(J)
           ENDDO
         ENDDO
         C(:,:)=C(:,:)+WORK(:,:)
         IF(TTEST) THEN
           DO IU=NU+1,NB
             I=MAP(IU)
             DO JU=1,NU
               J=MAP(JU)
               CALL TESTB(I,J,CSVAR)
               IF(ABS(CSVAR).GT.TOL) THEN
!                WRITE(*,FMT='("HIGHER CHIS ORTHOGONALIZED TO PHI(",I4,")")')N
!                PRINT*,I,J,CSVAR
               END IF
             ENDDO
            ENDDO
         END IF
       ENDDO
!
!      =================================================================
!      == TEST ORTHOGONALITY                                          ==
!      =================================================================
       IF(TTEST) THEN
         MAXDEV=0.D0
         DO I=1,NB
           DO J=1,NB
             CSVAR=PHIPHI(I,J)
             DO K=1,NB
               CSVAR=CSVAR+CONJG(X(K,I))*CHIPHI(K,J) &
      &                   +CONJG(CHIPHI(K,I))*X(K,J)
               DO L=1,NB
                 CSVAR=CSVAR+CONJG(X(K,I))*CHICHI(K,L)*X(L,J)
               ENDDO
             ENDDO
             IF(I.EQ.J) CSVAR=CSVAR-(1.D0,0.D0)
             MAXDEV=MAX(MAXDEV,ABS(CSVAR))
             IF(ABS(CSVAR).GT.TOL) THEN
               CALL ERROR$MSG('ORTHOGONALIZATION FAILED')
               CALL ERROR$I4VAL('I',I)
               CALL ERROR$I4VAL('J',J)
               CALL ERROR$C8VAL('<PHI(+)|O|PHI(+)>-1',CSVAR)
               CALL ERROR$STOP('WAVES_ORTHO_Y')
             END IF
           ENDDO
         ENDDO
         PRINT*,'MAX. DEVIATION IN ORTHO_Y_C',MAXDEV
       END IF
                             CALL TRACE$POP
       RETURN
       CONTAINS
!      .................................................................
       SUBROUTINE TESTA(I,J,CSVAR)
       INTEGER(4),INTENT(IN) :: I
       INTEGER(4),INTENT(IN) :: J
       COMPLEX(8),INTENT(OUT):: CSVAR
!      *****************************************************************
       CSVAR=PHIPHI(I,J)
       DO K=1,NB
         CSVAR=CSVAR+CONJG(X(K,I))*CHIPHI(K,J)+CONJG(CHIPHI(K,I))*X(K,J)
         DO L=1,NB
           CSVAR=CSVAR+CONJG(X(K,I))*CHICHI(K,L)*X(L,J)
         ENDDO
       ENDDO
       IF(I.EQ.J) CSVAR=CSVAR-(1.D0,0.D0)
       RETURN
       END SUBROUTINE TESTA
!      .................................................................
       SUBROUTINE TESTB(I,J,CSVAR)
       INTEGER(4),INTENT(IN) :: I
       INTEGER(4),INTENT(IN) :: J
       COMPLEX(8),INTENT(OUT):: CSVAR
!      *****************************************************************
       CSVAR=CHIPHI(I,J)
       DO K=1,NB
         CSVAR=CSVAR+CONJG(ALPHA(K,I))*CHIPHI(K,J)+CHIPHI(I,K)*X(K,J)
         DO L=1,NB
           CSVAR=CSVAR+CONJG(ALPHA(K,I))*CHICHI(K,L)*X(L,J)
         ENDDO
       ENDDO
       RETURN
       END SUBROUTINE TESTB
      END
!
!      .................................................................
       SUBROUTINE WAVES_ORTHO_Y(NB,PHIPHI0,CHIPHI0,CHICHI0,RLAMBDA)
!      **                                                             **
!      **  CALCULATE LAGRANGE MULTIPLIERS FOR ORTHOGONALIZATION       **
!      **    |PHI_N(+)>=|PHIBAR_N>+SUM_M |CHI_M(0)>LAMBDA(M,N)        **
!      **  WITH                                                       **
!      **    |CHI>=O|PHI>                                             **
!      **    |PHIBAR>=2|PHI(0)>-PHI(0)+H|PSI(0)>DT^2/M_PSI            **
!      **    LAMDA(I>J)=0                                             **
!      **                                                             **
!      **  ATTENTION!! CHIPHI0=<CHI|O|PHI>                            **
!      **        AND   PHIPHI0=<PHIBAR|O|PHIBAR>                      **
!      **        CONVERSION IS DONE IN THE INITIALIZATION             **
!      **                                                             **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NB
       REAL(8)   ,INTENT(IN) :: PHIPHI0(NB,NB) !<PHIBAR_I|O|PHIBAR_J>
       REAL(8)   ,INTENT(IN) :: CHIPHI0(NB,NB) !<PHIBAR_I|O|CHI>
       REAL(8)   ,INTENT(IN) :: CHICHI0(NB,NB) !<CHI_I|O|CHI_J>
       REAL(8)   ,INTENT(OUT):: RLAMBDA(NB,NB) ! (I>J)=0
       REAL(8)               :: PHIPHI(NB,NB)
       REAL(8)               :: CHIPHI(NB,NB)
       REAL(8)               :: CHICHI(NB,NB)
       REAL(8)               :: ALPHA(NB,NB)   ! (I>J)=0
       REAL(8)               :: WORK(NB,NB)   ! (I>J)=0
       REAL(8)               :: DELTA(NB)
       INTEGER(4)            :: I,J,K,L,N,M,M1,M2
       INTEGER(4)            :: NU,MU,M1U,IU
       REAL(8)               :: SVAR
       LOGICAL   ,PARAMETER  :: TPR=.FALSE.
       LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
       REAL(8)   ,PARAMETER  :: TOL=1.D-10
       REAL(8)               :: TEST(NB,NB)
       REAL(8)               :: TEST1(NB,NB)
       REAL(8)               :: TEST2(NB,NB)
       INTEGER(4)            :: MAP(NB)
!      *****************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_Y')
!
!      =================================================================
!      ==  INITIALIZE                                                 ==
!      =================================================================
       DO I=1,NB
         MAP(I)=I
         DO J=1,NB
           PHIPHI(I,J)=PHIPHI0(I,J)
           CHIPHI(I,J)=CHIPHI0(I,J)
           CHICHI(I,J)=CHICHI0(I,J)
           RLAMBDA(I,J)=0.D0
           ALPHA(I,J)=0.D0
         ENDDO
         PHIPHI(I,I)=PHIPHI(I,I)-1.D0
         ALPHA(I,I)=1.D0
       ENDDO
                             CALL TRACE$PASS('BEFORE TESTME')
       IF(TPR)CALL TESTME
!
!      =================================================================
!      ==  ORTHOGONALIZATION LOOP                                     ==
!      =================================================================
                             CALL TRACE$PASS('BEFORE ORTHOGONALIZATION LOOP')
       DO NU=1,NB
         N=MAP(NU)
!
!        ===============================================================
!        == NORMALIZE PHI(N)                                          ==
!        == PHI(N)=PHI(N)+CHI(N)*SVAR                                 ==
!        ===============================================================
                             CALL TRACE$PASS('NORMALIZE')
         SVAR=1.D0-PHIPHI(N,N)*CHICHI(N,N)/CHIPHI(N,N)**2
         IF(SVAR.GT.0.D0) THEN
           SVAR=CHIPHI(N,N)/CHICHI(N,N)*(DSQRT(SVAR)-1.D0)             
         ELSE
           PRINT*,'ORTHOGONALIZATION FAILED! TRYING BEST APPROXIMATION...'
           SVAR=-CHIPHI(N,N)/CHICHI(N,N)
         END IF       
         DO I=1,NB
           RLAMBDA(I,N)=RLAMBDA(I,N)+ALPHA(I,N)*SVAR
         ENDDO
         PHIPHI(N,N)=PHIPHI(N,N)+CHICHI(N,N)*SVAR**2
         DO I=1,NB
           PHIPHI(I,N)=PHIPHI(I,N)+CHIPHI(N,I)*SVAR
           PHIPHI(N,I)=PHIPHI(N,I)+CHIPHI(N,I)*SVAR
         ENDDO
         DO I=1,NB
           CHIPHI(I,N)=CHIPHI(I,N)+SVAR*CHICHI(N,I)
         ENDDO
         IF(TPR) THEN
           PRINT*,'AFTER NORMALIZATION'
           CALL TESTME
         END IF
!
!        ===============================================================
!        == ORTHOGONALIZE HIGHER PHI'S TO THIS PHI                    ==
!        == PHI(M)=PHI(M)+CHI(N)*DELTA(M)   M>N                       ==
!        ===============================================================
                CALL TRACE$PASS('ORTHOGONALIZE HIGHER PHIS TO THIS PHI')
         DELTA=0.D0
         DO MU=NU+1,NB
           M=MAP(MU)
!          == PHIPHI(N,M)+CHIPHI(N,N)*DELTA(M)=0  ======================
           DELTA(M)=-PHIPHI(N,M)/CHIPHI(N,N)
           DO IU=1,NU
             I=MAP(IU)
             RLAMBDA(I,M)=RLAMBDA(I,M)+ALPHA(I,N)*DELTA(M)
           ENDDO
         ENDDO           
         DO M1=1,NB
           DO M2=1,NB
             PHIPHI(M1,M2)=PHIPHI(M1,M2)+DELTA(M1)*CHIPHI(N,M2) &
     &                    +CHIPHI(N,M1)*DELTA(M2) &
     &                    +DELTA(M1)*CHICHI(N,N)*DELTA(M2)
           ENDDO
         ENDDO
         DO M1U=NU+1,NB
           M1=MAP(M1U)
           DO M2=1,NB
             CHIPHI(M2,M1)=CHIPHI(M2,M1)+DELTA(M1)*CHICHI(N,M2)
           ENDDO
         ENDDO
         IF(TPR) THEN
           WRITE(*,FMT='("HIGHER PHIS ORTHOGONALIZED TO PHI(",I4,")")')N
           CALL TESTME
         END IF
!
!        ===============================================================
!        == ORTHOGONALIZE HIGHER CHI'S TO THIS PHI                    ==
!        == CHI(M)=CHI(M)+CHI(N)*DELTA(M)   M>N                       ==
!        ===============================================================
                CALL TRACE$PASS('ORTHOGONALIZE HIGHER CHIS TO THIS PHI')
         DELTA=0.D0
         DO MU=NU+1,NB
           M=MAP(MU)
!          == |CHI(M)>=|CHI(M)>+|CHI(N)>*DELTA(M) ============================
!          == CHIPHI(M,N)+CHIPHI(N,N)*DELTA(M)=0
           DELTA(M)=-CHIPHI(M,N)/CHIPHI(N,N)
         ENDDO
         DO MU=NU+1,NB
           M=MAP(MU)
           DO I=1,NB
             ALPHA(I,M)=ALPHA(I,M)+ALPHA(I,N)*DELTA(M)
           ENDDO
           DO I=1,NB
             CHIPHI(M,I)=CHIPHI(M,I)+DELTA(M)*CHIPHI(N,I)
           ENDDO 
         ENDDO
         WORK(:,:)=0.D0
         DO M1=1,NB
           DO M2=1,NB
             WORK(M1,M2)=WORK(M1,M2)+DELTA(M1)*CHICHI(N,M2) &
        &                                    +CHICHI(M1,N)*DELTA(M2) &
        &                          +DELTA(M1)*CHICHI(N,N) *DELTA(M2)
           ENDDO
         ENDDO
         CHICHI(:,:)=CHICHI(:,:)+WORK(:,:)
         IF(TPR) THEN
           WRITE(*,FMT='("HIGHER CHIS ORTHOGONALIZED TO PHI(",I4,")")')
           CALL TESTME
         END IF
       ENDDO
!
!      =================================================================
!      == TEST ORTHOGONALITY                                          ==
!      =================================================================
       IF(TTEST) THEN
         CALL TESTME
         DO I=1,NB
           DO J=1,NB
             SVAR=PHIPHI0(I,J)
             DO K=1,NB
               SVAR=SVAR+CHIPHI0(K,I)*RLAMBDA(K,J)+RLAMBDA(K,I)*CHIPHI0(K,J)
               DO L=1,NB
                 SVAR=SVAR+RLAMBDA(K,I)*CHICHI0(K,L)*RLAMBDA(L,J)
               ENDDO
             ENDDO
             IF(I.EQ.J) SVAR=SVAR-1.D0
             IF(ABS(SVAR).GT.TOL) THEN
               CALL ERROR$MSG('ORTHOGONALIZATION FAILED')
               CALL ERROR$I4VAL('I',I)
               CALL ERROR$I4VAL('J',J)
               CALL ERROR$R8VAL('<PHI(+)|O|PHI(+)>-1',SVAR)
               CALL ERROR$STOP('WAVES_ORTHO_Y')
             END IF
           ENDDO
         ENDDO
       END IF
                             CALL TRACE$POP
       RETURN
       CONTAINS
!      ............................................................
       SUBROUTINE TESTME
       REAL(8)              :: SVAR,SVAR1,SVAR2
       WRITE(*,FMT='("STARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTART")')
       DO I=1,NB
         DO J=1,NB
           SVAR=PHIPHI0(I,J)
           SVAR1=0.D0
           SVAR2=0.D0
           DO K=1,NB
             SVAR=SVAR+RLAMBDA(K,I)*CHIPHI0(K,J)+CHIPHI0(K,I)*RLAMBDA(K,J)
             SVAR1=SVAR1+CHIPHI0(K,I)*ALPHA(K,J)
             DO L=1,NB
               SVAR=SVAR+RLAMBDA(K,I)*CHICHI0(K,L)*RLAMBDA(L,J)
               SVAR1=SVAR1+RLAMBDA(K,I)*CHICHI0(K,L)*ALPHA(L,J)
               SVAR2=SVAR2+ALPHA(K,I)*CHICHI0(K,L)*ALPHA(L,J)
             ENDDO
           ENDDO
           TEST(I,J)=SVAR
           TEST1(I,J)=SVAR1
           TEST2(I,J)=SVAR2
         ENDDO
         TEST(I,I)=TEST(I,I)-1.D0
       ENDDO
       DO I=1,NB
         WRITE(*,FMT='("TEST  =",8E10.2)')TEST(I,:)
       ENDDO
       DO I=1,NB
         WRITE(*,FMT='("PHIPHI=",8E10.2)')PHIPHI(I,:)
       ENDDO
       DO I=1,NB
         WRITE(*,FMT='("CHIPHI=",8E10.2)')CHIPHI(I,:)
       ENDDO
       DO I=1,NB
         WRITE(*,FMT='("CHICHI=",8E10.2)')CHICHI(I,:)
       ENDDO
       DO I=1,NB
         WRITE(*,FMT='("LAMBDA=",8F10.3)')RLAMBDA(I,:)
       ENDDO
       DO I=1,NB
         WRITE(*,FMT='("ALPHA =",8F10.3)')ALPHA(I,:)
       ENDDO
       DO I=1,NB
         DO J=1,NB
           IF(DABS(TEST(I,J)-PHIPHI(I,J)).GT.1.D-7) THEN
!             PRINT*,'ERROR IN PHIPHI FOR ',I,J,PHIPHI(I,J),TEST(I,J)
           END IF
           IF(DABS(TEST1(I,J)-CHIPHI(J,I)).GT.1.D-7) THEN
!             PRINT*,'ERROR IN CHIPHI FOR ',I,J,CHIPHI(J,I),TEST1(I,J)
           END IF
           IF(DABS(TEST2(I,J)-CHICHI(I,J)).GT.1.D-7) THEN
!             PRINT*,'ERROR IN CHICHI FOR ',I,J,CHICHI(I,J),TEST2(I,J)
           END IF
         ENDDO
       ENDDO
       WRITE(*,FMT='("ENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDEND")')
       RETURN
       END SUBROUTINE TESTME
      END
!
!     .....................................................ORTHO........
      SUBROUTINE WAVES_ORTHO_X_C(NB,OCC,CHICHI,PSIPSI,CHIPSI,LAMBDA)
!     ******************************************************************
!     **                                                              **
!     **  IMPOSES THE ORTHOGONALITY CONSTRAINT ONTO THE ELECTRONS     **
!     **                                                              **
!     **  NEW VERSION WITH DIAGONALIZATION FOR CHIPSI                    **
!     **                                                              **
!     **                                                              **
!     **  THE METHOD IS DESCRIBED IN :                                **
!     **    R.CAR AND M.PARRINELLO, IN "SIMPLE MOLECULAR SYSTEMS      **
!     **    AT VERY HIGH DENSITY", PAGE 455                           **
!     **    ED. A.POLIAN, PLOUBEYRE AND N.BOCCARA                     **
!     **    (PLENUM PUBLISHING CORPORATION,1989)                      **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
!     IMPLICIT NONE
      REAL(8)   ,PARAMETER     :: EPS    = 1.D-8
      REAL(8)   ,PARAMETER     :: DSMALL = 1.D-12
      INTEGER(4),PARAMETER     :: ITERX    = 200
      INTEGER(4),INTENT(IN)    :: NB
      REAL(8)   ,INTENT(IN)    :: OCC(NB)
      COMPLEX(8),INTENT(INOUT) :: LAMBDA(NB,NB)
      COMPLEX(8),INTENT(IN)    :: PSIPSI(NB,NB)
      COMPLEX(8),INTENT(IN)    :: CHIPSI(NB,NB)   !
      COMPLEX(8),INTENT(IN)    :: CHICHI(NB,NB)
      COMPLEX(8),ALLOCATABLE   :: GAMN(:,:)  
      REAL(8)                  :: EIG(NB)
      INTEGER(4)               :: IND,ITER,I,J,K ! RUNNING VARIABLES
      INTEGER(4)               :: IMAX,I0,J0   ! AUXILARY VARIABLES
      REAL(8)                  :: DIGAM,SVAR,FI,FJ,EIGI ! AUXILARY VARIABLES
      COMPLEX(8)               :: HAUX(NB,NB)    
      COMPLEX(8)               :: U(NB,NB)       
      LOGICAL(4)               :: TCONVERGED
      LOGICAL(4)               :: TESSL=.TRUE.
      COMPLEX(8)               :: CSVAR
      REAL(8)                  :: OCCI,OCCJ
      INTEGER(4),EXTERNAL      :: IDAMAX
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_X_C')
      ALLOCATE(GAMN(NB,NB))
!
!     ==================================================================
!     ==  CALCULATE  PSIPSI(I,J)= <PSIBAR(I)|PSIBAR(J)>-1             ==
!     ==        AND  CHIPSI(I,J)   = <PSI0(I)|PSIBAR(J)>              ==
!     ==================================================================
!
!     ==================================================================
!     ==  DIAGONALIZE 0.5*(CHIPSI(I,J)+CHIPSI(J,I))                   ==
!     ==================================================================
      DO I=1,NB
        DO J=1,NB
           GAMN(I,J)=0.5D0*(CHIPSI(I,J)+CONJG(CHIPSI(J,I)))
        ENDDO
      ENDDO
      CALL LIB$DIAGC8(NB,GAMN,EIG,U)
!CALL CDIAG(NB,NB,GAMN,EIG,U)
!WRITE(*,FMT='("EIG",20E10.3)')EIG
!DO I=1,NB
!  WRITE(*,FMT='("U",I2,20E10.3)')I,U(I,:)
!ENDDO

!
!     ==================================================================
!     ==================================================================
!     ==  ITERATIVE CALCULATION OF GAMMA                              ==
!     ==================================================================
!     ==================================================================
      DO ITER=1,ITERX
!PRINTA*,'==================',ITER,'==========================='
!       ================================================================
!       ==  CALCULATE <PHI(+)|PHI(+)>-1 WITH PRESENT LAMBDA           ==
!       ==  GAMN(I,J)=PSIPSI(I,J)+LAMBDA(K,I)*CHIPSI(K,J)             ==
!       ==                       +CHIPSI(K,I)*LAMBDA(K,J)             == 
!       ==           +LAMBDA(K,I)*CHICHI(K,L)*LAMBDA(L,J)-1(I,J)      ==
!       ================================================================
        IF(TESSL) THEN
!         __GAMN(I,J) = CHICHI(I,K)*LAMBDA(K,J)___________________________
          CALL LIB$MATMULC8(NB,NB,NB,CHICHI,LAMBDA,HAUX)
!CALL ZGEMUL(CHICHI,NB,'N',LAMBDA,NB,'N',HAUX,NB,NB,NB,NB)
!         __HAUX(I,J) = HAUX(I,J)+2*CHIPSI(I,J)___________________________
          HAUX=HAUX+2.D0*CHIPSI
!CALL ZAXPY(NB*NB,(2.D0,0.D0),CHIPSI,1,HAUX,1)
!         __GAMN(I,J) = LAMBDA(K,I)*HAUX(K,J)_____________________________
          CALL LIB$SCALARPRODUCTC8(.FALSE.,NB,NB,LAMBDA,NB,HAUX,GAMN)
!CALL ZGEMUL(LAMBDA,NB,'C',HAUX,NB,'N',GAMN,NB,NB,NB,NB)
!         __HAUX(I,J) = HAUX(I,J)+2*CHIPSI(I,J)___________________________
          GAMN(:,:)=GAMN(:,:)+PSIPSI(:,:)
!         __GAMN(I,J) = GAMN(I,J)-1_______________________________________
          DO I=1,NB
            DO J=I,NB
              CSVAR=0.5D0*(GAMN(I,J)+CONJG(GAMN(J,I)))
              GAMN(I,J)=CSVAR
              GAMN(J,I)=CONJG(CSVAR)
            ENDDO
            GAMN(I,I)=GAMN(I,I)-1.D0
          ENDDO
        ELSE
          DO I=1,NB
            DO J=1,NB
              GAMN(I,J)=PSIPSI(I,J)
              IF(I.EQ.J) GAMN(I,J)=GAMN(I,J)-(1.D0,0.D0)
              HAUX(I,J)=CHIPSI(I,J)
              DO K=1,NB
                GAMN(I,J)=GAMN(I,J)+CONJG(CHIPSI(K,I))*LAMBDA(K,J)
                HAUX(I,J)=HAUX(I,J)+CHICHI(I,K)*LAMBDA(K,J)
              ENDDO
            ENDDO
          ENDDO
          DO I=1,NB
            DO J=1,NB
              DO K=1,NB
                GAMN(I,J)=GAMN(I,J)+CONJG(LAMBDA(K,I))*HAUX(K,J)
              ENDDO
            ENDDO
          ENDDO
        END IF
!DO I=1,NB
!WRITE(*,FMT='("A",I2,20E10.3)')I,GAMN(I,:)
!ENDDO
!
!       ================================================================
!       == FIND LARGEST ELEMENT OF THE OVERLAP MATRIX                 ==
!       ================================================================
        DIGAM=0.D0
        TCONVERGED=.TRUE.
        DO I=1,NB
          DO J=I,NB
            IF(ABS(GAMN(I,J)).GT.EPS) THEN
              TCONVERGED=.FALSE.
              DIGAM=MAX(DIGAM,ABS(GAMN(I,J)))
            END IF
          ENDDO
        ENDDO
        IF(TCONVERGED) GOTO 9000
PRINT*,'ITER ',ITER,DIGAM
!
!       ==================================================================
!       ==  OBTAIN CHANGE OF THE LAMBDA MATRIX                          ==
!       ==================================================================
!       == TRANSFORM OVERLAP MATRIX GAMN
        IF(TESSL) THEN
!         ----  HAUX(I,L)=U(K,I)*H0(K,L)
          CALL LIB$SCALARPRODUCTC8(.FALSE.,NB,NB,U,NB,GAMN,HAUX)
!CALL ZGEMUL(U,NB,'C',GAMN,NB,'N',HAUX,NB,NB,NB,NB)
!         ----  GAMN(I,J)=HAUX(I,L)*U(L,I)
          CALL LIB$MATMULC8(NB,NB,NB,HAUX,U,GAMN)
!CALL ZGEMUL(HAUX,NB,'N',U,NB,'N',GAMN,NB,NB,NB,NB)
!
!         ==  MULTIPLY WITH 1/(EIG(I)+EIG(J))
          DO I=1,NB
            EIGI=EIG(I)
            DO J=1,NB
              GAMN(I,J)=GAMN(I,J)/(EIGI+EIG(J))
            ENDDO
          ENDDO
!
!         == TRANSFORM OVERLAP MATRIX GAMN BACK
!         ----  HAUX(I,L)=U(K,I)*H0(K,L)
          CALL LIB$MATMULC8(NB,NB,NB,U,GAMN,HAUX)
!CALL ZGEMUL(U,NB,'N',GAMN,NB,'N',HAUX,NB,NB,NB,NB)
!         ----  GAMN(I,J)=HAUX(I,L)*U(L,I)
          CALL LIB$DYADSUMC8(NB,NB,NB,HAUX,U,GAMN)
!CALL ZGEMUL(HAUX,NB,'N',U,NB,'C',GAMN,NB,NB,NB,NB)
        ELSE
          DO I=1,NB
            DO J=1,NB
              HAUX(I,J)=(0.D0,0.D0)
              DO K=1,NB
                HAUX(I,J)=HAUX(I,J)+CONJG(U(K,I))*GAMN(K,J)
              ENDDO
            ENDDO
          ENDDO
          DO I=1,NB
            DO J=1,NB
              GAMN(I,J)=(0.D0,0.D0)
              DO K=1,NB
                GAMN(I,J)=GAMN(I,J)+HAUX(I,K)*U(K,J)
              ENDDO
              GAMN(I,J)=GAMN(I,J)/(EIG(I)+EIG(J))
            ENDDO
          ENDDO
          DO I=1,NB
            DO J=1,NB
              HAUX(I,J)=(0.D0,0.D0)
              DO K=1,NB
                HAUX(I,J)=HAUX(I,J)+U(I,K)*GAMN(K,J)
              ENDDO
            ENDDO
          ENDDO
          DO I=1,NB
            DO J=1,NB
              GAMN(I,J)=(0.D0,0.D0)
              DO K=1,NB
                GAMN(I,J)=GAMN(I,J)+HAUX(I,K)*CONJG(U(J,K))
              ENDDO
            ENDDO
          ENDDO
        END IF
!
!       ================================================================
!       ==  PROPAGATE GAMMA                                           ==
!       ================================================================
        DO I=1,NB
          OCCI=OCC(I)+DSMALL
          DO J=1,NB
            OCCJ=OCC(J)+DSMALL
            SVAR=2.D0*OCCI/(OCCI+OCCJ)
            LAMBDA(I,J)=LAMBDA(I,J)-SVAR*GAMN(I,J)
          ENDDO
        ENDDO
!DO I=1,NB
!WRITE(*,FMT='("DL",I2,20E10.3)')I,GAMN(I,:)
!ENDDO
!
!       ================================================================
!       == SYMMETRIZE LAMBDA                                          ==
!       ================================================================
        DO I=1,NB
          DO J=1,NB
            IF(OCC(I).LE.OCC(J)) THEN
              LAMBDA(I,J)=CONJG(LAMBDA(J,I))*(OCC(I)+DSMALL)/(OCC(J)+DSMALL)
            END IF
          ENDDO
        ENDDO
!DO I=1,NB
!WRITE(*,FMT='("L",I2,20E10.3)')I,GAMN(I,:)
!ENDDO
!
!       ================================================================
!       ==  ALTERNATIVE                                               ==
!       ================================================================
!       DO I=1,NB
!         FI=F(I)
!         DO J=I+1,NB
!           FJ=F(J)
!           IF(FI+FJ.GT.1.D-6) THEN
!             FI=1.D0
!             FJ=0.D0
!             LAMBDA(J,I)=0.D0
!           ELSE
!             IF(FI.LE.FJ) THEN
!               LAMBDA(I,J)=LAMBDA(J,I)*FI/FJ
!             ELSE IF(FI.GT.FJ) THEN
!               LAMBDA(J,I)=LAMBDA(I,J)*FJ/FI
!             END IF
!           ENDIF
!           SVAR=1.D0/(FI*EIG(I)+FJ*EIG(J))
!           LAMBDA(I,J)=LAMBDA(I,J)-SVAR*GAMN(I,J)*FI
!           LAMBDA(J,I)=LAMBDA(J,I)-SVAR*GAMN(J,I)*FJ
!         ENDDO
!       ENDDO
!
!       DO I=1,NB
!         LAMBDA(I,I)=LAMBDA(I,I)-GAMN(I,I)/(2.D0*EIG(I))
!       ENDDO
!
      ENDDO
      CALL ERROR$MSG('LOOP FOR ORTHOGONALIZATION IS NOT CONVERGED')
      CALL ERROR$STOP('WAVES_ORTHO_X_C')
!
9000  CONTINUE
      DEALLOCATE(GAMN)
                             CALL TRACE$POP
      RETURN
      END
!
!     .....................................................ORTHO........
      SUBROUTINE WAVES_ORTHO_X(NB,OCC,CHICHI,PSIPSI,CHIPSI,LAMBDA)
!     ******************************************************************
!     **                                                              **
!     **  IMPOSES THE ORTHOGONALITY CONSTRAINT ONTO THE ELECTRONS     **
!     **                                                              **
!     **  NEW VERSION WITH DIAGONALIZATION FOR CHIPSI                    **
!     **                                                              **
!     **                                                              **
!     **  THE METHOD IS DESCRIBED IN :                                **
!     **    R.CAR AND M.PARRINELLO, IN "SIMPLE MOLECULAR SYSTEMS      **
!     **    AT VERY HIGH DENSITY", PAGE 455                           **
!     **    ED. A.POLIAN, PLOUBEYRE AND N.BOCCARA                     **
!     **    (PLENUM PUBLISHING CORPORATION,1989)                      **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
!     IMPLICIT NONE
      REAL(8)   ,PARAMETER     :: EPS    = 1.D-8
      REAL(8)   ,PARAMETER     :: DSMALL = 1.D-12
      INTEGER(4),PARAMETER     :: MAX    = 200
      INTEGER(4),INTENT(IN)    :: NB
      REAL(8)   ,INTENT(IN)    :: OCC(NB)
      REAL(8)   ,INTENT(INOUT) :: LAMBDA(NB,NB)
      REAL(8)   ,INTENT(IN)    :: PSIPSI(NB,NB)
      REAL(8)   ,INTENT(IN)    :: CHIPSI(NB,NB)   !
      REAL(8)   ,INTENT(IN)    :: CHICHI(NB,NB)
      REAL(8)   ,ALLOCATABLE   :: GAMN(:,:)  
      REAL(8)                  :: EIG(NB)
      INTEGER(4)               :: IND,ITER,I,J ! RUNNING VARIABLES
      INTEGER(4)               :: IMAX,I0,J0   ! AUXILARY VARIABLES
      REAL(8)                  :: DIGAM,SVAR,FI,FJ,EIGI ! AUXILARY VARIABLES
      REAL(8)                  :: HAUX(NB,NB)    
      REAL(8)                  :: U(NB,NB)       
      REAL(8)                  :: OCCI,OCCJ
      INTEGER(4),EXTERNAL      :: IDAMAX
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_X')
      ALLOCATE(GAMN(NB,NB))
!
!     ==================================================================
!     ==  CALCULATE  PSIPSI(I,J)= <PSIBAR(I)|PSIBAR(J)>-1             ==
!     ==        AND  CHIPSI(I,J)   = <PSI0(I)|PSIBAR(J)>              ==
!     ==================================================================
!
!     ==================================================================
!     ==  DIAGONALIZE 0.5*(CHIPSI(I,J)+CHIPSI(J,I))                   ==
!     ==================================================================
      CALL LIB$DIAGR8(NB,CHIPSI,EIG,U)
!CALL DIAG(NB,NB,CHIPSI,EIG,U)
!WRITE(*,FMT='("EIG",20E10.3)')EIG
!DO I=1,NB
!  WRITE(*,FMT='("U",I2,20E10.3)')I,U(I,:)
!ENDDO
!
!     ==================================================================
!     ==================================================================
!     ==  ITERATIVE CALCULATION OF GAMMA                              ==
!     ==================================================================
!     ==================================================================
      DO ITER=1,MAX
!PRINT*,'==================',ITER,'==========================='
!       ================================================================
!       ==  CALCULATE <PHI(+)|PHI(+)>-1 WITH PRESENT LAMBDA           ==
!       ==  GAMN(I,J)=PSIPSI(I,J)+LAMBDA(K,I)*CHIPSI(K,J)             ==
!       ==                       +CHIPSI(K,I)*LAMBDA(K,J)             == 
!       ==           +LAMBDA(K,I)*CHICHI(K,L)*LAMBDA(L,J)-1(I,J)      ==
!       ================================================================
!       __GAMN(I,J) = CHICHI(I,K)*LAMBDA(K,J)___________________________
        CALL LIB$MATMULR8(NB,NB,NB,CHICHI,LAMBDA,HAUX)
!CALL DGEMUL(CHICHI,NB,'N',LAMBDA,NB,'N',HAUX,NB,NB,NB,NB)
!       __HAUX(I,J) = HAUX(I,J)+2*CHIPSI(I,J)___________________________
        HAUX=HAUX+2.D0*CHIPSI
!CALL DAXPY(NB*NB,2.D0,CHIPSI,1,HAUX,1)
!       __GAMN(I,J) = LAMBDA(K,I)*HAUX(K,J)_____________________________
        CALL LIB$SCALARPRODUCTR8(.FALSE.,NB,NB,LAMBDA,NB,HAUX,GAMN)
!CALL DGEMUL(LAMBDA,NB,'T',HAUX,NB,'N',GAMN,NB,NB,NB,NB)
!       __GAMN(I,J) = GAMN(I,J)-1_______________________________________
        DO I=1,NB
          DO J=I,NB
            SVAR=0.5D0*(GAMN(I,J)+GAMN(J,I))+PSIPSI(I,J)
            GAMN(J,I)=SVAR
            GAMN(I,J)=SVAR
          ENDDO
          GAMN(I,I)=GAMN(I,I)-1.D0
        ENDDO
!DO I=1,NB
!WRITE(*,FMT='("A",I2,20E10.3)')I,GAMN(I,:)
!ENDDO
!
!       ================================================================
!       == FIND LARGEST ELEMENT OF THE OVERLAP MATRIX                 ==
!       ================================================================
        IMAX=IDAMAX(NB*NB,GAMN(1,1),1)
        J0=(IMAX-1)/NB+1
        I0=IMAX-(J0-1)*NB
        DIGAM=DABS(GAMN(I0,J0))
!       PRINT*,'ITER ',ITER,I0,J0,DIGAM,NCON
        IF(DIGAM.LT.EPS) GOTO 9000
PRINT*,'ITER ',ITER,DIGAM
!
!       ==================================================================
!       ==  OBTAIN CHANGE OF THE LAMBDA MATRIX                          ==
!       ==================================================================
!       == TRANSFORM OVERLAP MATRIX GAMN
!       ----  HAUX(I,L)=U(K,I)*H0(K,L)
        CALL LIB$SCALARPRODUCTR8(.FALSE.,NB,NB,U,NB,GAMN,HAUX)
!CALL DGEMUL(U,NB,'T',GAMN,NB,'N',HAUX,NB,NB,NB,NB)
!       ----  GAMN(I,J)=HAUX(I,L)*U(L,I)
        CALL LIB$MATMULR8(NB,NB,NB,HAUX,U,GAMN)
!CALL DGEMUL(HAUX,NB,'N',U,NB,'N',GAMN,NB,NB,NB,NB)
!
!       ==  MULTIPLY WITH 1/(EIG(I)+EIG(J))
        DO I=1,NB
          EIGI=EIG(I)
          DO J=1,NB
            GAMN(I,J)=GAMN(I,J)/(EIGI+EIG(J))
          ENDDO
        ENDDO
!
!       == TRANSFORM OVERLAP MATRIX GAMN BACK
!       ----  HAUX(I,L)=U(K,I)*H0(K,L)
        CALL LIB$MATMULR8(NB,NB,NB,U,GAMN,HAUX)
!CALL DGEMUL(U,NB,'N',GAMN,NB,'N',HAUX,NB,NB,NB,NB)
!       ----  GAMN(I,J)=HAUX(I,L)*U(L,I)
        CALL LIB$DYADSUMR8(NB,NB,NB,HAUX,U,GAMN)
!CALL DGEMUL(HAUX,NB,'N',U,NB,'T',GAMN,NB,NB,NB,NB)
!DO I=1,NB
!WRITE(*,FMT='("DL",I2,20E10.3)')I,GAMN(I,:)
!ENDDO
!
!       ================================================================
!       ==  PROPAGATE GAMMA                                           ==
!       ================================================================
        DO I=1,NB
          OCCI=OCC(I)+DSMALL
          DO J=1,NB
            OCCJ=OCC(J)+DSMALL
            SVAR=2.D0*OCCI/(OCCI+OCCJ)
            LAMBDA(I,J)=LAMBDA(I,J)-SVAR*GAMN(I,J)
          ENDDO
        ENDDO
!DO I=1,NB
!WRITE(*,FMT='("L",I2,20E10.3)')I,GAMN(I,:)
!ENDDO
!
!       ================================================================
!       == SYMMETRIZE LAMBDA                                          ==
!       ================================================================
        DO I=1,NB
          DO J=1,NB
            IF(OCC(I).LT.OCC(J)) THEN
              LAMBDA(I,J)=LAMBDA(J,I)*(OCC(I)+DSMALL)/(OCC(J)+DSMALL)
            END IF
          ENDDO
        ENDDO
!
!       ================================================================
!       ==  ALTERNATIVE                                               ==
!       ================================================================
!       DO I=1,NB
!         FI=F(I)
!         DO J=I+1,NB
!           FJ=F(J)
!           IF(FI+FJ.GT.1.D-6) THEN
!             FI=1.D0
!             FJ=0.D0
!             LAMBDA(J,I)=0.D0
!           ELSE
!             IF(FI.LE.FJ) THEN
!               LAMBDA(I,J)=LAMBDA(J,I)*FI/FJ
!             ELSE IF(FI.GT.FJ) THEN
!               LAMBDA(J,I)=LAMBDA(I,J)*FJ/FI
!             END IF
!           ENDIF
!           SVAR=1.D0/(FI*EIG(I)+FJ*EIG(J))
!           LAMBDA(I,J)=LAMBDA(I,J)-SVAR*GAMN(I,J)*FI
!           LAMBDA(J,I)=LAMBDA(J,I)-SVAR*GAMN(J,I)*FJ
!         ENDDO
!       ENDDO
!
!       DO I=1,NB
!         LAMBDA(I,I)=LAMBDA(I,I)-GAMN(I,I)/(2.D0*EIG(I))
!       ENDDO
!
      ENDDO
!      PRINT*,'EIG ',EIG
!      PRINT*,'OCC ',OCC
!      DO I=1,NB
!        DO J=1,NB
!          WRITE(*,FMT='(2I3,7F10.5)')I,J,LAMBDA(I,J),U(I,J),PSIPSI(I,J) &
!    &                               ,CHIPSI(I,J),CHICHI(I,J)
!        ENDDO
!      ENDDO
      CALL ERROR$MSG('LOOP FOR ORTHOGONALIZATION IS NOT CONVERGED')
      CALL ERROR$STOP('WAVES_ORTHO_X')
      
!
9000  CONTINUE
      DEALLOCATE(GAMN)
                             CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_GRAMSCHMIDT(MAP,GSET,NAT,R,NGL,NDIM,NBH,NB,PSI)
!     ******************************************************************
!     **                                                              **
!     **  GRAM-SCHMIDT ORTHOGONALIZATION OF A SET OF WAVE FUNCTIONS   **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1999)***
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY: MAP_TYPE,GSET_TYPE
      IMPLICIT NONE
      TYPE(MAP_TYPE)  ,INTENT(IN) :: MAP
      TYPE(GSET_TYPE) ,INTENT(IN) :: GSET
      INTEGER(4)      ,INTENT(IN) :: NAT
      REAL(8)         ,INTENT(IN) :: R(3,NAT)
      INTEGER(4)      ,INTENT(IN) :: NGL
      INTEGER(4)      ,INTENT(IN) :: NDIM
      INTEGER(4)      ,INTENT(IN) :: NBH
      INTEGER(4)      ,INTENT(IN) :: NB
      COMPLEX(8)      ,INTENT(INOUT):: PSI(NGL,NDIM,NBH)
      INTEGER(4)                  :: NPRO
      INTEGER(4)                  :: IDIM,IG,I,J
      INTEGER(4)                  :: IBH1,IBH2,IB1A,IB1B,IB2A,IB2B
      COMPLEX(8)                  :: CSVAR
      COMPLEX(8)      ,ALLOCATABLE:: X(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: PROJ(:,:,:)
      COMPLEX(8)      ,ALLOCATABLE:: OVERLAP(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: AUXMAT(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: X1(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: X2(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: PSIINV(:,:,:)
      COMPLEX(8)      ,ALLOCATABLE:: CWORK(:,:,:)
      COMPLEX(8)                  :: CSVAR1,CSVAR2
      COMPLEX(8)      ,PARAMETER  :: CI=(0.D0,1.D0)
      LOGICAL(4)                  :: TINV
      REAL(8)                     :: NORM(NB),SVAR
      COMPLEX(8)                  :: XTWOBYTWO(2,2)
      LOGICAL(4)      ,PARAMETER  :: TTEST=.false.
      INTEGER(4)      ,ALLOCATABLE:: SMAP(:)
!     ******************************************************************
!
!     ==================================================================
!     ==  CALCULATE PROJECTIONS FOR THE NEW POSITIONS                 ==
!     ==================================================================
      CALL PLANEWAVE$SELECT(GSET%ID)
      TINV=GSET%TINV
      NPRO=MAP%NPRO
      ALLOCATE(PROJ(NDIM,NBH,NPRO))
      CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NBH,NPRO,PSI,PROJ)
!     
!     ==================================================================
!     ==  OVERLAP OF <PSI0|PSI0>,                                     ==
!     ==================================================================
      ALLOCATE(OVERLAP(NB,NB))
      ALLOCATE(AUXMAT(NB,NB))
      CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB,PSI,PSI,OVERLAP)
      CALL WAVES_1COVERLAP(.TRUE.,MAP,NDIM,NBH,NB,NPRO,PROJ,PROJ,AUXMAT)
      DO J=1,NB
        DO I=1,NB
          OVERLAP(I,J)=OVERLAP(I,J)+AUXMAT(I,J)
        ENDDO
      ENDDO
      DEALLOCATE(AUXMAT)
      DEALLOCATE(PROJ)
!     
!     =================================================================
!     ==  OBTAIN ORTHOGONALIZATION TRANSFORM                         ==
!     =================================================================
      ALLOCATE(SMAP(NB))
      DO I=1,NB
        SMAP(I)=I
      ENDDO
      ALLOCATE(X(NB,NB))
      CALL WAVES_ORTHO_Y_C(NB,OVERLAP,OVERLAP,OVERLAP,X,SMAP)
      DEALLOCATE(OVERLAP)
      DEALLOCATE(SMAP)
!     
!     =================================================================
!     ==  TRANSFORM WAVE FUNCTIONS  |PSI>=|PSI>X                     ==
!     =================================================================
      IF(TINV) THEN
        ALLOCATE(X1(NBH,NBH))
        ALLOCATE(X2(NBH,NBH))
!
!       == FIRST FOLD DOWN OVERLAP MATRIX TO SUPER WAVE FUNCTIONS ======
        DO IBH1=1,NBH
          IB1A=2*IBH1-1
          IB1B=MIN(2*IBH1,NB)
          DO IBH2=1,NBH
            IB2A=2*IBH2-1
            IB2B=MIN(2*IBH2,NB)
            XTWOBYTWO(:,:)=(0.D0,0.D0)
            DO I=IB1A,IB1B
              DO J=IB2A,IB2B
                XTWOBYTWO(I-IB1A+1,J-IB2A+1)=X(I,J)
              ENDDO
            ENDDO
            CSVAR1=    XTWOBYTWO(1,1)+CI*XTWOBYTWO(1,2)
            CSVAR2=-CI*XTWOBYTWO(2,1)   +XTWOBYTWO(2,2)
            X1(IBH1,IBH2)=0.5D0*(CSVAR1+CSVAR2)
            X2(IBH1,IBH2)=0.5D0*(CSVAR1-CSVAR2)
          ENDDO
        ENDDO
        DEALLOCATE(X)
        ALLOCATE(PSIINV(NGL,NDIM,NBH))
        PSIINV=PSI
        CALL PLANEWAVE$ADDPRODUCT(' ',NGL,NDIM,NBH,PSI,NBH,PSIINV,X1)
        CALL PLANEWAVE$ADDPRODUCT('-',NGL,NDIM,NBH,PSI,NBH,PSIINV,X2)
        DEALLOCATE(PSIINV)
!        DO I=NBH,1,-1
!          DO J=1,I       !WORKS ONLY FOR TRIANGULAR X
!            DO IDIM=1,NDIM
!              DO IG=1,NGL
!                PSI(IG,IDIM,I)=PSI(IG,IDIM,I) &
!     &                        +PSI(IG,IDIM,J)   *X1(J,I) &
!     &                        +PSIINV(IG,IDIM,J)*X2(J,I)
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDDO         
        DEALLOCATE(X1)
        DEALLOCATE(X2)
      ELSE
        ALLOCATE(PSIINV(NGL,NDIM,NB))
        PSIINV=PSI
        CALL PLANEWAVE$ADDPRODUCT(' ',NGL,NDIM,NB,PSI,NB,PSIINV,X)
        DEALLOCATE(PSIINV)
!        DO I=NB,1,-1
!          DO J=1,I
!            DO IDIM=1,NDIM
!              DO IG=1,NGL
!                PSI(IG,IDIM,I)=PSI(IG,IDIM,I)+PSI(IG,IDIM,J)*X(J,I)
!              ENDDO
!            ENDDO
!          ENDDO
!       ENDDO         
        DEALLOCATE(X)
      ENDIF
! 
!     =================================================================
!     ==                                                             ==
!     =================================================================
      IF(TTEST) THEN
        ALLOCATE(PROJ(NDIM,NBH,NPRO))
        CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NBH,NPRO,PSI,PROJ)
        ALLOCATE(AUXMAT(NB,NB))
        CALL WAVES_1COVERLAP(.TRUE.,MAP,NDIM,NBH,NB,NPRO,PROJ,PROJ,AUXMAT)
        DEALLOCATE(PROJ)
        ALLOCATE(OVERLAP(NB,NB))
        CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB,PSI,PSI,OVERLAP)
        DO J=1,NB
          DO I=1,NB
            OVERLAP(I,J)=OVERLAP(I,J)+AUXMAT(I,J)
          ENDDO
        ENDDO
        DEALLOCATE(AUXMAT)
        ALLOCATE(X(NB,NB))
        CALL LIB$DIAGC8(NB,OVERLAP,NORM,X)
!CALL CDIAG(NB,NB,OVERLAP,NORM,X)
        DEALLOCATE(OVERLAP)
        DO I=1,NB
          IF(ABS(NORM(I)-1.D0).GT.1.D-4) THEN
            CALL ERROR$MSG('GRAM-SCHMIDT ORTHOGONALIZATION FAILED')
            CALL ERROR$I4VAL('STATE',I)
            CALL ERROR$R8VAL('NORM',NORM(I))
            CALL ERROR$STOP('WAVES_GRAMSCHMIDT')
          END IF
        ENDDO
        DEALLOCATE(X)
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SWITCH()
!     ******************************************************************
!     **  STEP FORWARD                                                **
!     **                                                              **
!     **  REMARKS: THE PROPAGATION MUST BE CONSISTENT WITH THAT IN    **
!     **   ORTHOGONALIZE                                              **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: IKPT
      INTEGER(4)             :: ISPIN
      COMPLEX(8),POINTER     :: TPSI(:,:,:)
      LOGICAL(4)             :: TSTRESS
      REAL(8)                :: CELLSCALE
      REAL(8)                :: MAPTOCELL(3,3)
      INTEGER(4)             :: NLAMBDA
      INTEGER(4),PARAMETER   :: NLAMBDAX=2
      INTEGER(4)             :: NB
      COMPLEX(8),ALLOCATABLE :: RLAMP(:,:)
      COMPLEX(8),POINTER     :: RLAMPTR(:,:)
!     ******************************************************************
      CALL CELL$GETL4('MOVE',TSTRESS)
      WAVEEKIN1=0.D0
      WAVEEKIN2=0.D0
      IF(TSTRESS) THEN
        CALL CELL$GETR8A('MAPTOCELL',9,MAPTOCELL)
        CELLSCALE=MAPTOCELL(1,1)*(MAPTOCELL(2,2)*MAPTOCELL(3,3)  &
     &                           -MAPTOCELL(3,2)*MAPTOCELL(2,3)) &
     &           +MAPTOCELL(2,1)*(MAPTOCELL(3,2)*MAPTOCELL(1,3)  &
     &                           -MAPTOCELL(1,2)*MAPTOCELL(3,3)) &
     &           +MAPTOCELL(3,1)*(MAPTOCELL(1,2)*MAPTOCELL(2,3)  &
     &                           -MAPTOCELL(2,2)*MAPTOCELL(1,3))
        CELLSCALE=1.D0/SQRT(CELLSCALE)
      ENDIF
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          TPSI=>THIS%PSI0
          THIS%PSI0=>THIS%PSIM
          THIS%PSIM=>TPSI
          IF(TSTRESS) THEN
            THIS%PSI0=THIS%PSI0*CELLSCALE
            THIS%PSIM=THIS%PSIM*CELLSCALE
          END IF
!
!         ==============================================================
!         == EXTRAPOLATE LAGRANGE MULTIPLIERS                         ==
!         ==============================================================
          NLAMBDA=4
          IF(.NOT.ASSOCIATED(THIS%RLAM3M)) NLAMBDA=3
          IF(.NOT.ASSOCIATED(THIS%RLAM2M)) NLAMBDA=2
          IF(.NOT.ASSOCIATED(THIS%RLAMM)) NLAMBDA=1
          IF(.NOT.ASSOCIATED(THIS%RLAM0)) NLAMBDA=0
          NLAMBDA=MIN(NLAMBDA,NLAMBDAX)
          NB=THIS%NB
          ALLOCATE(RLAMP(NB,NB))
          IF(NLAMBDA.EQ.0) THEN
            RLAMP=0.D0
          ELSE IF(NLAMBDA.EQ.1) THEN
            RLAMP=THIS%RLAM0
            THIS%RLAMM=>THIS%RLAM0
            NULLIFY(THIS%RLAM0)
            IF(NLAMBDA.EQ.NLAMBDAX) THEN
              THIS%RLAM0=>THIS%RLAMM
              NULLIFY(THIS%RLAMM)
            END IF 
          ELSE IF(NLAMBDA.EQ.2) THEN
            RLAMP=2.D0*THIS%RLAM0-THIS%RLAMM
            THIS%RLAM2M=>THIS%RLAMM
            THIS%RLAMM=>THIS%RLAM0
            NULLIFY(THIS%RLAM0)
            IF(NLAMBDA.EQ.NLAMBDAX) THEN  ! REUSE MEMORY 
              THIS%RLAM0=>THIS%RLAM2M     
              NULLIFY(THIS%RLAM2M)        
            END IF 
          ELSE IF(NLAMBDA.EQ.3) THEN
            RLAMP=3.D0*THIS%RLAM0-3.D0*THIS%RLAMM+THIS%RLAM2M
            THIS%RLAM3M=>THIS%RLAM2M
            THIS%RLAM2M=>THIS%RLAMM
            THIS%RLAMM=>THIS%RLAM0
            NULLIFY(THIS%RLAM0)
            IF(NLAMBDA.EQ.NLAMBDAX) THEN
              THIS%RLAM0=>THIS%RLAM3M
              NULLIFY(THIS%RLAM3M)
            END IF 
          ELSE IF(NLAMBDA.EQ.4) THEN
            RLAMP=4.D0*THIS%RLAM0-6.D0*THIS%RLAMM+4.D0*THIS%RLAM2M-THIS%RLAM3M
            RLAMPTR=>THIS%RLAM3M
            THIS%RLAM3M=>THIS%RLAM2M
            THIS%RLAM2M=>THIS%RLAMM
            THIS%RLAMM=>THIS%RLAM0
            THIS%RLAM0=>RLAMPTR
            NULLIFY(RLAMPTR)
          END IF
          IF(.NOT.ASSOCIATED(THIS%RLAM0))ALLOCATE(THIS%RLAM0(NB,NB))
          THIS%RLAM0=RLAMP             
          DEALLOCATE(RLAMP)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  DELETE HAMILTONIAN AND EIGENVALUES                          ==
!     ==================================================================
      IF(THAMILTON) THEN
        DO IKPT=1,NKPT
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            CALL PLANEWAVE$SELECT(GSET%ID)
            DEALLOCATE(THIS%EIGVAL)         
            DEALLOCATE(THIS%EIGVEC)         
            DEALLOCATE(THIS%EXPECTVAL)
          ENDDO
        ENDDO
        THAMILTON=.FALSE.
      END IF
      RETURN
      END
!
!     ......................................................RANWAV......
      SUBROUTINE WAVES_RANDOMIZE(NG,NDIM,NB,AMPLITUDE,G2,PSI)
!     ******************************************************************
!     **                                                              **
!     **  CREATE RANDOM WAVE FUNCTIONS                                **
!     **                                                              **
!     **  THE MAXIMUM WEIGHT OF THE WAVE FUNCTIONS AT EPW[RY]=GC2     **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB              ! #(BANDS)
      INTEGER(4),INTENT(IN)    :: NG              ! #(PLANE WAVES),MAX
      INTEGER(4),INTENT(IN)    :: NDIM            ! #(PLANE WAVES),MAX
      REAL(8)   ,INTENT(IN)    :: AMPLITUDE       ! SCALE FACTOR
      REAL(8)   ,INTENT(IN)    :: G2(NG)          ! G**2
      COMPLEX(8),INTENT(INOUT) :: PSI(NG,NDIM,NB) ! PS-WAVE FUNCTION
      INTEGER(4)               :: IB,IG,IDIM
      REAL(8)                  :: PI,GC,FAC
      REAL(8)   ,PARAMETER     :: GC2=10.D0
      REAL(8)                  :: SCALE(NG)
      REAL(8)                  :: REC,RIM
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      FAC=2.D0*SQRT(PI*GC2)
      FAC=FAC**3/REAL(NDIM,8)*(2.D0/3.D0)
      FAC=AMPLITUDE/FAC
      DO IG=1,NG
        SCALE(IG)=FAC*EXP(-0.5D0*G2(IG)/GC2)
      ENDDO
      CALL LIB$RANDOMSEED
      DO IB=1,NB
        DO IDIM=1,NDIM
          DO IG=1,NG
            CALL LIB$RANDOM(REC)
            CALL LIB$RANDOM(RIM)
            REC=2.D0*REC-1.D0
            RIM=2.D0*RIM-1.D0
            PSI(IG,IDIM,IB)=PSI(IG,IDIM,IB)+CMPLX(REC,RIM)*SCALE(IG)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ......................................................RANWAV......
      SUBROUTINE WAVES_INITIALIZERANDOM(NG,NDIM,NB,G2,PSI)
!     ******************************************************************
!     **                                                              **
!     **  CREATE RANDOM WAVE FUNCTIONS                                **
!     **                                                              **
!     **  THE MAXIMUM WEIGHT OF THE WAVE FUNCTIONS AT EPW[RY]=GC2     **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB              ! #(BANDS)
      INTEGER(4),INTENT(IN)    :: NG              ! #(PLANE WAVES),MAX
      INTEGER(4),INTENT(IN)    :: NDIM            ! #(PLANE WAVES),MAX
      REAL(8)   ,INTENT(IN)    :: G2(NG)          ! G**2
      COMPLEX(8),INTENT(OUT)   :: PSI(NG,NDIM,NB) ! PS-WAVE FUNCTION
      INTEGER(4)               :: IB,IG,IDIM
      REAL(8)                  :: PI,GC,FAC
      REAL(8)   ,PARAMETER     :: GC2=10.D0
      REAL(8)                  :: SCALE(NG)
      REAL(8)                  :: REC,RIM
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      FAC=2.D0*SQRT(PI*GC2)
      FAC=FAC**3/REAL(NDIM,8)*(2.D0/3.D0)
      FAC=1.D0/FAC
      DO IG=1,NG
        SCALE(IG)=FAC*EXP(-0.5D0*G2(IG)/GC2)
      ENDDO
      CALL LIB$RANDOMSEED
      DO IB=1,NB
        DO IDIM=1,NDIM
          DO IG=1,NG
            CALL LIB$RANDOM(REC)
            CALL LIB$RANDOM(RIM)
            REC=2.D0*REC-1.D0
            RIM=2.D0*RIM-1.D0
            PSI(IG,IDIM,IB)=CMPLX(REC,RIM)*SCALE(IG)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!======================================================================
!======================================================================
!======================================================================
!======================================================================
!======================================================================
!======================================================================
!
!     ..................................................................
      SUBROUTINE WAVES$REPORT(NFIL)
!     ******************************************************************
!     **  WAVES$GET                                                   **
!     **  GET                                                         **
!     ******************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)   ,PARAMETER  :: MBYTE=2.D0**20
      INTEGER(4)            :: IKPT
      REAL(8)               :: RY
      REAL(8)               :: SVAR
      REAL(8)               :: MEMORY
      INTEGER(4)            :: NG
!     ******************************************************************
      CALL CONSTANTS('RY',RY)
      CALL REPORT$TITLE(NFIL,'WAVE FUNCTIONS')
      CALL WAVES_SELECTWV(1,1)
      CALL REPORT$I4VAL(NFIL,'NUMBER OF BANDS',THIS%NB,' ')

      CALL REPORT$I4VAL(NFIL,'NUMBER OF K-POINTS',NKPT,' ')
      CALL REPORT$I4VAL(NFIL,'NUMBER OF SPINS',NSPIN,' ')
      CALL REPORT$I4VAL(NFIL,'NUMBER OF SPINOR COMPONENTS',NDIM,' ')
      CALL REPORT$R8VAL(NFIL,'PLANE WAVE CUTOFF',EPWPSI/RY,'RY')
      CALL REPORT$R8VAL(NFIL,'WAVE FUNCTION MASS',EMASS,'A.U.')
      CALL REPORT$R8VAL(NFIL,'G**2 ENHANCEMENT OF FUNCTION MASS',EMASSCG2,' ')
      CALL REPORT$R8VAL(NFIL,'BUCKET POTENTIAL STARTS AT 0.5G^2=',EPWPSI0/RY,'RY')
      CALL REPORT$R8VAL(NFIL,'BUCKET POTENTIAL PREFACTOR',D2EPWPSI,'H')
      IF(ANNEE.NE.0) THEN
        CALL REPORT$R8VAL(NFIL,'FRICTION',ANNEE,' ')
      END IF
      IF(TSTOP) THEN
        CALL REPORT$CHVAL(NFIL,'INITIAL VELOCITY IS SET TO','ZERO')
      END IF
!     IF(TRANDOMIZE) THEN
!       CALL REPORT$R8VAL(NFIL &
!    &      ,'INITIAL VELOCITIES ARE RANDOMIZED WITH ENERGY',AMPRE,'H')
!     END IF
      IF(.NOT.TSAFEORTHO) THEN
        WRITE(NFIL,FMT='("EIGENSTATES ARE CALCULATED."' &
     &                 //'," (NO STRICT ENERGY CONSERVATION)")')
      END IF
!     
!     ================================================================
!     ==  REPORT INFORMATION ABOUT G-VECTORS                        ==
!     ================================================================
      DO IKPT=1,NKPT
        CALL WAVES_SELECTWV(IKPT,1)
        NG=GSET%NGL
        CALL MPE$COMBINE('+',NG)
        CALL REPORT$I4VAL(NFIL &
     &       ,'NUMBER OF (GLOBAL) PLANE WAVES FOR WAVE FUNCTION',NG,' ')
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$WRITEPDOS
!     ******************************************************************
!     **  WAVES$GET                                                   **
!     **  GET                                                         **
!     ******************************************************************
      USE WAVES_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: NFIL
      REAL(8)   ,ALLOCATABLE :: VAL(:,:)
      REAL(8)   ,ALLOCATABLE :: DER(:,:)
      REAL(8)   ,ALLOCATABLE :: OV(:,:,:)
      REAL(8)   ,ALLOCATABLE :: R(:,:)
      REAL(8)   ,ALLOCATABLE :: RAD(:)
      REAL(8)   ,ALLOCATABLE :: WORK(:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)                :: AEZ
      INTEGER(4),ALLOCATABLE :: IZ(:)
      INTEGER(4)             :: NAT,NSP
      REAL(8)                :: R1,DEX,XEXP,RI
      INTEGER(4)             :: NR
      COMPLEX(8)             :: CSVAR
      INTEGER(4)             :: NB,NBH
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      COMPLEX(8),ALLOCATABLE :: VEC(:,:)      
      COMPLEX(8),ALLOCATABLE :: VECTOR1(:,:,:)
      INTEGER(4)             :: LNXX,LNX,NPRO
      INTEGER(4)             :: NGL
      INTEGER(4),ALLOCATABLE :: LOX(:)
      LOGICAL(4)             :: TINV
      REAL(8)                :: RBAS(3,3)
      INTEGER(4)             :: ISP,IR
      INTEGER(4)             :: NTASKS,THISTASK
      INTEGER(4)             :: IB1,IB2,IBH,LN1,LN2,IDIM,IKPT,ISPIN,IPRO
      COMPLEX(8),ALLOCATABLE :: PROJ(:,:,:)
      CHARACTER(16),ALLOCATABLE :: ATOMID(:)
      REAL(8)                :: SVAR      
      CHARACTER(32)          :: FLAG='011004'
      INTEGER(4)             :: NBX
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$WRITEPDOS')
      IF(.NOT.THAMILTON) THEN
        CALL ERROR$MSG('EIGENVALUES NOT PRESENT')
        CALL ERROR$STOP('WAVES$WRITEPDOS')
      END IF
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     == GENERAL QUANTITIES                                           ==
!     ==================================================================
      NSP=MAP%NSP
      NAT=MAP%NAT
      NPRO=MAP%NPRO
      LNXX=MAP%LNXX
      CALL PDOS$SETI4('NAT',NAT)
      CALL PDOS$SETI4('NSP',NSP)
      CALL PDOS$SETI4('NKPT',NKPT)
      CALL PDOS$SETI4('NSPIN',NSPIN)
      CALL PDOS$SETI4('NDIM',NDIM)
      CALL PDOS$SETI4('NPRO',NPRO)
      CALL PDOS$SETI4('LNXX',LNXX)
      CALL PDOS$SETI4A('LNX',NSP,MAP%LNX)
      CALL PDOS$SETI4A('LOX',LNXX*NSP,MAP%LOX)
      CALL PDOS$SETI4A('ISPECIES',NAT,MAP%ISP)
      IF(THISTASK.EQ.1) THEN
        CALL FILEHANDLER$UNIT('PDOS',NFIL)
        REWIND NFIL
        WRITE(NFIL)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,FLAG
        WRITE(NFIL)MAP%LNX(:),MAP%LOX(:,:),MAP%ISP(:)
      END IF
!
!     ==================================================================
!     == ATOMIC STRUCTURE                                             ==
!     ==================================================================
      ALLOCATE(R(3,NAT))
      ALLOCATE(ATOMID(NAT))
      CALL CELL$GETR8A('T(0)',9,RBAS)
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
      CALL ATOMLIST$GETCHA('NAME',0,NAT,ATOMID)
      CALL PDOS$SETR8A('RBAS',9,RBAS)
      CALL PDOS$SETR8A('R',3*NAT,R)
!     CALL PDOS$SETCHA('ATOMID',NAT,ATOMID)
      IF(THISTASK.EQ.1) THEN
        WRITE(NFIL)RBAS,R,ATOMID
      END IF
      DEALLOCATE(ATOMID)
!
!     ==================================================================
!     == ELEMENT SPECIFIC QUANTITIES                                  ==
!     ==================================================================
      ALLOCATE(VAL(LNXX,NSP))
      ALLOCATE(DER(LNXX,NSP))
      ALLOCATE(OV(LNXX,LNXX,NSP))
      ALLOCATE(LOX(LNXX))
      ALLOCATE(IZ(NSP))
      ALLOCATE(RAD(NSP))
      CALL SETUP$GETI4('NR',NR)
      ALLOCATE(AEPHI(NR,LNXX))
      ALLOCATE(WORK(NR))
      CALL SETUP$GETR8('R1',R1)
      CALL SETUP$GETR8('DEX',DEX)
      XEXP=EXP(DEX)
      DO ISP=1,NSP
        OV(:,:,:)=0.D0
        CALL SETUP$ISELECT(ISP)
        LNX=MAP%LNX(ISP)
        LOX=MAP%LOX(:,ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        IZ(ISP)=NINT(AEZ)
        CALL PERIODICTABLE$GET(IZ(ISP),'R(ASA)',RAD(ISP))
        CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,AEPHI)
        DO LN1=1,LNX
          CALL RADIAL$VALUE(R1,DEX,NR,AEPHI(1,LN1),RAD(ISP),VAL(LN1,ISP))
          CALL RADIAL$DERIVATIVE(R1,DEX,NR,AEPHI(1,LN1),RAD(ISP),DER(LN1,ISP))
          DO LN2=LN1,LNX
            IF(LOX(LN1).NE.LOX(LN2)) CYCLE
            RI=R1/XEXP
            DO IR=1,NR
              RI=RI*XEXP
              WORK(IR)=RI**2*AEPHI(IR,LN1)*AEPHI(IR,LN2)
            ENDDO
            CALL RADIAL$INTEGRAL1(R1,DEX,NR,WORK,RAD(ISP),OV(LN1,LN2,ISP))
            OV(LN2,LN1,ISP)=OV(LN1,LN2,ISP)
          ENDDO
        ENDDO
        IF(THISTASK.EQ.1) THEN
          WRITE(NFIL)IZ(ISP),RAD(ISP),VAL(1:LNX,ISP),DER(1:LNX,ISP) &
     &              ,OV(1:LNX,1:LNX,ISP)
        END IF
      ENDDO
      CALL PDOS$SETI4A('IZ',NSP,IZ)
      CALL PDOS$SETR8A('RAD',NSP,RAD)
      CALL PDOS$SETR8A('PHI',LNXX*NSP,VAL)
      CALL PDOS$SETR8A('DPHIDR',LNXX*NSP,DER)
      CALL PDOS$SETR8A('OVERLAP',LNXX*LNXX*NSP,OV)
      DEALLOCATE(VAL)
      DEALLOCATE(IZ)
      DEALLOCATE(RAD)
      DEALLOCATE(DER)
      DEALLOCATE(OV)
      DEALLOCATE(WORK)
      DEALLOCATE(AEPHI)
      DEALLOCATE(LOX)
!
!     ==================================================================
!     ==  NOW WRITE PROJECTIONS                                       ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPT,NSPIN))
      CALL DYNOCC$GETR8A('OCC',NBX*NKPT*NSPIN,OCC)
      ALLOCATE(XK(3,NKPT))
      CALL DYNOCC$GETR8A('XK',3*NKPT,XK)
      CALL PDOS$SETR8A('XK',3*NKPT,XK)
      ALLOCATE(VEC(NDIM,NPRO))
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(.NOT.ASSOCIATED(THIS%EIGVEC)) THEN
            CALL ERROR$MSG('EIGENVALUES NOT PRESENT')
            CALL ERROR$STOP('WAVES$WRITEPDOS')
          END IF
          NB=THIS%NB
          NBH=THIS%NBH
          NGL=GSET%NGL
          ALLOCATE(VECTOR1(NDIM,NPRO,NB))
          TINV=GSET%TINV
          IF(THISTASK.EQ.1) THEN
            WRITE(NFIL)XK(:,IKPT),NB
          END IF
!
!         =============================================================
!         ==  CALCULATE PROJECTIONS                                  ==
!         =============================================================
          ALLOCATE(PROJ(NDIM,NBH,NPRO))
          CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NBH,NPRO &
     &                            ,THIS%PSI0,PROJ)
!
!         ==============================================================
!         == UNRAVEL SUPER WAVE FUNCTIONS                             ==
!         ==============================================================
          IF(TINV) THEN
            DO IBH=1,NBH
              IB1=1+2*(IBH-1)
              IB2=2+2*(IBH-1)
              DO IPRO=1,NPRO
                DO IDIM=1,NDIM
                  VECTOR1(IDIM,IPRO,IB1)=REAL(PROJ(IDIM,IBH,IPRO))
                  VECTOR1(IDIM,IPRO,IB2)=AIMAG(PROJ(IDIM,IBH,IPRO))
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO IB1=1,NB
              DO IPRO=1,NPRO
                DO IDIM=1,NDIM
                  VECTOR1(IDIM,IPRO,IB1)=PROJ(IDIM,IB1,IPRO)
                ENDDO
              ENDDO
            ENDDO
          END IF
          DEALLOCATE(PROJ)
!
!         ==============================================================
!         ==  TRANSFORM TO EIGENSTATES IF SAFEORTHO=.TRUE.,           ==
!         ==  AND WRITE                                               ==
!         ==============================================================
          IF(TSAFEORTHO) THEN
            DO IB1=1,NB
              VEC=0.D0
              DO IB2=1,NB
                CSVAR=THIS%EIGVEC(IB2,IB1)
                DO IPRO=1,NPRO
                  DO IDIM=1,NDIM
                    VEC(IDIM,IPRO)=VEC(IDIM,IPRO)+VECTOR1(IDIM,IPRO,IB2)*CSVAR
                  ENDDO
                ENDDO
              ENDDO
              IF(THISTASK.EQ.1) THEN
                IF(FLAG.EQ.'011004') THEN
                  WRITE(NFIL)THIS%EIGVAL(IB1),OCC(IB1,IKPT,ISPIN),VEC
                ELSE
                  WRITE(NFIL)THIS%EIGVAL(IB1),VEC
                END IF
!PRINT*,'E ',THIS%EIGVAL(IB1)
!WRITE(*,FMT='(10F10.5)')VEC
              END IF
            ENDDO
          ELSE
            DO IB1=1,NB
              IF(THISTASK.EQ.1) THEN
                IF(FLAG.EQ.'011004') THEN
                  WRITE(NFIL)THIS%EXPECTVAL(IB1),OCC(IB1,IKPT,ISPIN) &
    &                       ,VECTOR1(:,:,IB1)
                ELSE
                  WRITE(NFIL)THIS%EXPECTVAL(IB1),VECTOR1(:,:,IB1)
                END IF
              END IF
            ENDDO
          END IF
          DEALLOCATE(VECTOR1)
        ENDDO
      ENDDO
      DEALLOCATE(XK)
      DEALLOCATE(VEC)
      DEALLOCATE(R)
      DEALLOCATE(OCC)
      IF(THISTASK.EQ.1) THEN
        CALL LIB$FLUSHFILE(NFIL)
        CALL FILEHANDLER$CLOSE('PDOS')
      END IF
                             CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$REPORTEIG(NFIL)
!     ******************************************************************
!     **  WAVES$GET                                                   **
!     **  GET                                                         **
!     ******************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: IKPT
      INTEGER(4)            :: ISPIN
      INTEGER(4)            :: IB
      INTEGER(4)            :: ITEN
      REAL(8)               :: EV
      CHARACTER(64)         :: STRING
      INTEGER(4)            :: NB
!     ******************************************************************
      CALL CONSTANTS('EV',EV)
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(.NOT.ASSOCIATED(THIS%EIGVAL)) THEN
            CALL ERROR$MSG('EIGENVALUES NOT PRESENT')
            CALL ERROR$STOP('WAVES$REPORTEIG')
          END IF
          NB=THIS%NB
          IF(NSPIN.EQ.1) THEN
            WRITE(STRING,FMT='("EIGENVALUES [EV] FOR K-POINT ",I4)')IKPT
          ELSE
            WRITE(STRING,FMT='("EIGENVALUES [EV] FOR K-POINT ",I4' &
     &                        //'," AND SPIN ",I1)')IKPT,ISPIN
          END IF
          CALL REPORT$TITLE(NFIL,STRING)
          ITEN=0
          DO WHILE (NB.GT.ITEN)
            WRITE(NFIL,FMT='(I3,":",10F8.3)') &
     &           ITEN,(THIS%EIGVAL(IB)/EV,IB=ITEN+1,MIN(ITEN+10,NB))
            ITEN=ITEN+10
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$WRITE(NFIL,NFILO,TCHK)
!     ******************************************************************
!     **                                                              **
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                        **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      USE RESTART_INTERFACE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN) :: NFIL
      INTEGER(4)            ,INTENT(IN) :: NFILO
      LOGICAL(4)            ,INTENT(OUT):: TCHK
      TYPE (SEPARATOR_TYPE),PARAMETER  :: SEP_WAVES &
           =SEPARATOR_TYPE(0,'WAVES','NONE','AUG1996','NONE')
      TYPE (SEPARATOR_TYPE)            :: SEPARATOR
      INTEGER(4)                       :: THISTASK,NTASKS
      INTEGER(4)                       :: IKPT,ISPIN
      INTEGER(4)                       :: NB,NBH
!     ******************************************************************
              CALL TRACE$PUSH('WAVES$WRITE')
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     == WRITE WAVE FUNCTION FOR THE THIS TIME STEP ==================
      TCHK=.TRUE.
      SEPARATOR=SEP_WAVES
      SEPARATOR%NREC=-1
      CALL WAVES_WRITEPSI(NFIL,SEPARATOR%NREC)
      CALL WRITESEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      CALL WAVES_WRITEPSI(NFIL,SEPARATOR%NREC)
              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$READ(NFIL,NFILO,TCHK)
!     ******************************************************************
!     **                                                              **
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                        **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **  --TCHK RETURNS IF ANYTHING HAS BEEN READ. AS RETURN CODE    **
!     **    IT SHOULD BETTER BE THAT ALL RELEVANT QUANTITIES HAVE     **
!     **    BEEN READ                                                 **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      USE RESTART_INTERFACE
      IMPLICIT NONE
      INTEGER(4)           ,INTENT(IN) :: NFIL
      INTEGER(4)           ,INTENT(IN) :: NFILO
      LOGICAL(4)           ,INTENT(OUT):: TCHK  ! SOMETHING HAS BEEN READ
      TYPE (SEPARATOR_TYPE),PARAMETER  :: SEP_WAVES &
           =SEPARATOR_TYPE(0,'WAVES','NONE','AUG1996','NONE')
      TYPE (SEPARATOR_TYPE)            :: SEPARATOR
      LOGICAL(4)                       :: TREAD
      REAL(8)              ,ALLOCATABLE:: EIG(:,:,:)
      INTEGER(4)                       :: IKPT,ISPIN,IB
      INTEGER(4)                       :: NKPT1,NSPIN1,NB1,NBX
      INTEGER(4)                       :: NB
      COMPLEX(8)           ,ALLOCATABLE:: TMP(:,:)
      REAL(8)              ,ALLOCATABLE:: TMPR8(:,:)
      INTEGER(4)                       :: THISTASK,NTASKS
!     ******************************************************************
                                  CALL TRACE$PUSH('WAVES$READ')
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     ==  READ WAVE FUNCTIONS FOR THIS TIME STEP                      ==
!     ==================================================================
      TCHK=.TRUE.
      SEPARATOR=SEP_WAVES
      CALL READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      IF(.NOT.TCHK) RETURN
      CALL WAVES_READPSI(NFIL)
!
!     ==================================================================
!     == SET OCCUPATIONS FOR DYNOCC OBJECT                            ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(EIG(NBX,NKPT,NSPIN))
      EIG(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          DO IB=1,THIS%NB
            EIG(IB,IKPT,ISPIN)=REAL(THIS%RLAM0(IB,IB))
          ENDDO
        ENDDO
      ENDDO
      CALL DYNOCC$SETR8A('EPSILON',NBX*NKPT*NSPIN,EIG)
      DEALLOCATE(EIG)
!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
                                  CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_WRITEPSI(NFIL,NREC)
!     ******************************************************************
!     **                                                              **
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                        **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL 
      INTEGER(4)  ,INTENT(INOUT):: NREC
      COMPLEX(8)  ,ALLOCATABLE:: PSIG(:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSI1(:,:)
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)              :: NTASKS,THISTASK
      INTEGER(4)              :: IOS
      INTEGER(4)              :: IKPT,ISPIN,IB,IDIM,IWAVE
      INTEGER(4)              :: LEN
      INTEGER(4)              :: NG1,NGG,NGL,NBH,NB
      LOGICAL(4)              :: TSUPER
      CHARACTER(64)           :: IOSTATMSG
      REAL(8)                 :: XK(3)
      REAL(8)     ,ALLOCATABLE:: GVECL(:,:)
      REAL(8)     ,ALLOCATABLE:: GVECG(:,:)
      REAL(8)                 :: GBAS(3,3)
      INTEGER(4)              :: NREC1,ISVAR
      CHARACTER(8)            :: KEY
      REAL(8)                 :: RBAS(3,3)
      INTEGER(4)  ,ALLOCATABLE:: IGVECG(:,:)
      INTEGER(4)  ,ALLOCATABLE:: IGVECL(:,:)
      INTEGER(4)              :: NWAVE=2
!     ******************************************************************
              CALL TRACE$PUSH('WAVES_WRITEPSI')
      CALL MPE$QUERY(NTASKS,THISTASK)
      IF(RSTRTTYPE.EQ.'STATIC') THEN
        NWAVE=1
      ELSE IF(RSTRTTYPE.EQ.'DYNAMIC') THEN
        NWAVE=2
      ELSE
        CALL ERROR$STOP('WAVES_WRITEPSI')
      END IF
!
!     ==================================================================
!     ==  COUNT #(RECORDS)                                            ==
!     ==================================================================
      IF(NREC.EQ.-1) THEN
        NREC=1
        DO IKPT=1,NKPT
          NREC=NREC+2
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          TSUPER=GSET%TINV
          NREC=NREC+NWAVE*THIS%NBH*NSPIN
          ISVAR=-1
          CALL WAVES_WRITELAMBDA(NFIL,IKPT,ISVAR)
          NREC=NREC+ISVAR
        ENDDO  
        CALL TRACE$POP
        RETURN
      END IF
      NREC1=0
!
!     ==================================================================
!     ==  GET DIMENSIONS ETC.                                         ==
!     ==================================================================
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     ==  WRITE SIZES                                                 ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
!       WRITE(NFIL)NKPT,NSPIN !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        CALL CELL$GETR8A('T(0)',9,RBAS)
        WRITE(NFIL)NKPT,NSPIN,RBAS,NWAVE !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        NREC1=NREC1+1
      END IF
!
!     ==================================================================
!     ==  LOOP OVER K-POINTS AND SPINS                                ==
!     ==================================================================
              CALL TRACE$PASS('MARKE 5')
      DO IKPT=1,NKPT
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NGL=GSET%NGL
        NBH=THIS%NBH
        NB=THIS%NB
        CALL PLANEWAVE$GETI4('NGG',NGG)
!       
!       ================================================================
!       == WRITE SIZE AND TYPE OF WAVE FUNCTION                       ==
!       ================================================================
        CALL PLANEWAVE$GETL4('TINV',TSUPER)
        IF(THISTASK.EQ.1) THEN
          KEY='PSI'
          WRITE(NFIL)KEY,NGG,NDIM,NB,TSUPER  !<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC1=NREC1+1
        END IF
!       
!       ================================================================
!       == WRITE K-POINTS AND G-VECTORS                               ==
!       ================================================================
        CALL PLANEWAVE$GETR8A('GBAS',9,GBAS)
        CALL PLANEWAVE$GETR8A('XK',3,XK)
        ALLOCATE(IGVECL(3,NGL))
        CALL PLANEWAVE$GETI4A('IGVEC',3*NGL,IGVECL)
        ALLOCATE(IGVECG(3,NGG))
        CALL PLANEWAVE$COLLECTI4(3,NGL,IGVECL,NGG,IGVECG)
        IF(THISTASK.EQ.1) THEN
          WRITE(NFIL)XK,IGVECG !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC1=NREC1+1
        END IF
        DEALLOCATE(IGVECG)
        DEALLOCATE(IGVECL)
        ALLOCATE(PSIG(NGG,NDIM))
        ALLOCATE(PSI1(NGL,NDIM))
        DO IWAVE=1,NWAVE
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            CALL PLANEWAVE$SELECT(GSET%ID)
            DO IB=1,NBH
              DO IDIM=1,NDIM
                IF(IWAVE.EQ.1) THEN
                  PSI1(:,IDIM)=THIS%PSI0(:,IDIM,IB)
                ELSE IF(IWAVE.EQ.2) THEN
                  PSI1(:,IDIM)=THIS%PSIM(:,IDIM,IB)
                END IF
              ENDDO
              DO IDIM=1,NDIM
                CALL PLANEWAVE$COLLECTC8(1,NGL,PSI1(1,IDIM),NGG,PSIG(1,IDIM))
              ENDDO
              IF(THISTASK.EQ.1) THEN
                WRITE(NFIL,ERR=9999,IOSTAT=IOS)PSIG !<<<<<<<<<<<<<<<<<<<
                NREC1=NREC1+1
              END IF
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(PSI1)
        DEALLOCATE(PSIG)
        ISVAR=0
        CALL WAVES_WRITELAMBDA(NFIL,IKPT,ISVAR)
        NREC1=NREC1+ISVAR
      ENDDO 
 !
!     ==================================================================
!     ==  DEALLOCATE                                                  ==
!     ==================================================================
      IF(THISTASK.EQ.1.AND.NREC1.NE.NREC) THEN
        CALL ERROR$MSG('#(RECORDS WRITTEN DIFFERENT FROM TARGET')
        CALL ERROR$I4VAL('TARGET',NREC)
        CALL ERROR$I4VAL('ACTUAL',NREC1)
        CALL ERROR$STOP('WAVES_WRITEPSI')
      END IF
              CALL TRACE$POP
      RETURN
 9999 CONTINUE
      CALL FILEHANDLER$IOSTATMESSAGE(IOS,IOSTATMSG)
      CALL ERROR$MSG('ERROR WRITING WAVE FUNCTION TO RESTART FILE')
      CALL ERROR$I4VAL('IOS',IOS)
      CALL ERROR$CHVAL('IOSTATMSG',IOSTATMSG)
      CALL ERROR$I4VAL('IB',IB)
      CALL ERROR$I4VAL('IKPT',IKPT)
      CALL ERROR$I4VAL('ISPIN',ISPIN)
      CALL ERROR$I4VAL('NGG',NGG)
      CALL ERROR$STOP('WAVES_WRITEPSI')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_READPSI(NFIL)
!     ******************************************************************
!     **                                                              **
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                        **
!     **                                                              **
!     ******************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL 
      INTEGER(4)              :: NTASKS,THISTASK
      INTEGER(4)              :: NKPT_,NSPIN_
      INTEGER(4)              :: NGG_,NDIM_,NB_,NBH_
      LOGICAL(4)              :: TSUPER_,TSUPER
      INTEGER(4)              :: IKPT_,ISPIN_
      INTEGER(4)              :: IKPT1,IKPT,IB,IDIM,IG,IB1,IB2,ISPIN,I
      INTEGER(4)              :: IKPT0
      INTEGER(4)              :: IWAVE
      LOGICAL(4)              :: TREAD(NKPT)
      REAL(8)                 :: KREAD(3,NKPT)
      REAL(8)                 :: K(3)        ! ACTUAL K-POINT 
      REAL(8)                 :: K_(3)       ! K-POINT ON FILE
      REAL(8)                 :: GBAS_(3,3)  ! REC. LATT. VECT. ON FILE
      REAL(8)     ,ALLOCATABLE:: GVECG_(:,:) ! G-VECTORS ON FILE
      REAL(8)     ,ALLOCATABLE:: GVECG(:,:)  ! G-VECTORS (GLOBAL)
      REAL(8)     ,ALLOCATABLE:: GVECL(:,:)  ! G-VECTORS (LOCAL)
      REAL(8)                 :: RBAS(3,3)
      REAL(8)                 :: SVAR,SVAR1
      INTEGER(4)              :: NGG,NGL,NBH,NB
      INTEGER(4)  ,ALLOCATABLE:: MAPG(:)
      COMPLEX(8)  ,ALLOCATABLE:: PSI(:,:,:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSIL(:,:,:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSIIN(:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSIG(:,:)
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)              :: ISVAR1
      INTEGER(4)              :: IOS
      CHARACTER(8)            :: KEY
      LOGICAL(4)              :: TCYCLE
      LOGICAL(4)              :: GBASFIX
      INTEGER(4)              :: IFORMAT
      INTEGER(4) ,ALLOCATABLE :: IGVECG_(:,:)
      REAL(8)                 :: XG1,XG2,XG3
      INTEGER(4)              :: NWAVE
      INTEGER(4) ,ALLOCATABLE :: MINUSG(:)
      LOGICAL(4)              :: TCHK
      COMPLEX(8)              :: F1,F2
      COMPLEX(8),ALLOCATABLE  :: PSIINSUPER(:,:)
      COMPLEX(8)              :: csvar,cmat(1,1)
!     ******************************************************************
                               CALL TRACE$PUSH('WAVES_READPSI')
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     ==  READ SIZES                                                  ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
!       IFORMAT=1
!        READ(NFIL,ERR=100)NKPT_,NSPIN_,RBAS
!        IFORMAT=2
! 100    CONTINUE
        NWAVE=2
        IFORMAT=0
        READ(NFIL,ERR=100)NKPT_,NSPIN_
        IFORMAT=1
        BACKSPACE(NFIL)
        READ(NFIL,ERR=100)NKPT_,NSPIN_,RBAS
        IFORMAT=2
        BACKSPACE(NFIL)
        READ(NFIL,ERR=100)NKPT_,NSPIN_,RBAS,NWAVE
 100    CONTINUE
        IF(IFORMAT.EQ.0) THEN
          CALL ERROR$MSG('FORMAT NOT RECOGNIZED')
          CALL ERROR$STOP('WAVES_READPSI')
        END IF
      END IF
      CALL MPE$BROADCAST(1,NSPIN_)
      CALL MPE$BROADCAST(1,NKPT_)
      CALL MPE$BROADCAST(1,IFORMAT)
      CALL MPE$BROADCAST(1,NWAVE)
      IF(IFORMAT.EQ.2) THEN
        CALL MPE$BROADCAST(1,RBAS)
        DO IKPT=1,NKPT
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
        ENDDO
        CALL GBASS(RBAS,GBAS_,SVAR)
      END IF
!
!     ==================================================================
!     ==  LOOP OVER K-POINTS AND SPINS                                ==
!     ==================================================================
              CALL TRACE$PASS('MARKE 5')
      GBASFIX=.FALSE.
      TREAD(:)=.FALSE.
      DO IKPT_=1,NKPT_
!
!       ================================================================
!       ==  READ COORDINATES OF THE WAVE FUNCTIONS                    ==
!       ================================================================
        IF(THISTASK.EQ.1) THEN
          READ(NFIL)KEY,NGG_,NDIM_,NB_,TSUPER_   !<<<<<<<<<<<<<<<<<<<<<<
        END IF
        CALL MPE$BROADCAST(1,KEY)
        CALL MPE$BROADCAST(1,NGG_)
        CALL MPE$BROADCAST(1,NDIM_)
        CALL MPE$BROADCAST(1,NB_)
        CALL MPE$BROADCAST(1,TSUPER_)
        IF(KEY.NE.'PSI') THEN
          CALL ERROR$MSG('ID IS NOT "PSI"')
          CALL ERROR$MSG('FILE IS CORRUPTED')
          CALL ERROR$I4VAL('IKPT_',IKPT_)
          CALL ERROR$STOP('WAVES_READPSI')
        END IF
        ALLOCATE(GVECG_(3,NGG_))
        IF(THISTASK.EQ.1) THEN
          IF(IFORMAT.EQ.1) THEN
            READ(NFIL)K_,GBAS_,GVECG_ !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          ELSE
            ALLOCATE(IGVECG_(3,NGG_))
            READ(NFIL)K_,IGVECG_ !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            DO IG=1,NGG_
              XG1=REAL(IGVECG_(1,IG),8)+K_(1)
              XG2=REAL(IGVECG_(2,IG),8)+K_(2)
              XG3=REAL(IGVECG_(3,IG),8)+K_(3)
              DO I=1,3
                GVECG_(I,IG)=GBAS_(I,1)*XG1+GBAS_(I,2)*XG2+GBAS_(I,3)*XG3
              ENDDO
            ENDDO
            IF(TSUPER_) THEN
              ALLOCATE(MINUSG(NGG_))
              CALL PLANEWAVE$MINUSG(K_,NGG_,IGVECG_,MINUSG)
            END IF
            DEALLOCATE(IGVECG_)
          END IF
        END IF
        CALL MPE$BROADCAST(1,K_)
        CALL MPE$BROADCAST(1,GVECG_)
!
        IF(IFORMAT.EQ.1.AND.(.NOT.GBASFIX)) THEN
          GBASFIX=.TRUE.
          CALL MPE$BROADCAST(1,GBAS_)
          CALL GBASS(GBAS_,RBAS,SVAR)
          DO IKPT1=1,NKPT
            CALL WAVES_SELECTWV(IKPT1,1)
            CALL PLANEWAVE$SELECT(GSET%ID)
            CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
          ENDDO
        END IF
!       
!       ==================================================================
!       ==  FIND NEAREST K-POINT                                        ==
!       ==================================================================
        SVAR=1.D+10
        IKPT=1
        DO IKPT1=1,NKPT
          CALL WAVES_SELECTWV(IKPT1,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETR8A('XK',3,K)
          SVAR1=(K(1)-K_(1))**2+(K(2)-K_(2))**2+(K(3)-K_(3))**2
          IF(SVAR1.LT.SVAR) THEN
            IKPT=IKPT1
            SVAR=SVAR1
          END IF
        ENDDO
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
!       
!       ==================================================================
!       ==  FIND NEAREST K-POINT                                        ==
!       ==================================================================
        IF(TREAD(IKPT)) THEN
          CALL PLANEWAVE$GETR8A('XK',3,K)
          SVAR1=(K(1)-K_(1))**2-(K(1)-KREAD(1,IKPT))**2 &
     &         +(K(2)-K_(2))**2-(K(2)-KREAD(2,IKPT))**2 &
     &         +(K(3)-K_(3))**2-(K(3)-KREAD(3,IKPT))**2 
          IF(SVAR1.GT.0.D0) THEN  ! PREVIOUS CHOICE WAS BETTER
            NBH_=NB_
            IF(TSUPER_)NBH_=INT(0.5D0*REAL(NB+1,8))
            DO IWAVE=1,NWAVE
              DO ISPIN_=1,NSPIN_
                DO IB=1,NBH_
                  IF(THISTASK.EQ.1)READ(NFIL)
                ENDDO
              ENDDO
            ENDDO
            CALL WAVES_READLAMBDA(NFIL,IKPT,.TRUE.)
            DEALLOCATE(GVECG_)
            IF(ALLOCATED(MINUSG))DEALLOCATE(MINUSG)
            CYCLE
          END IF
        END IF
        TREAD(IKPT)=.TRUE.
        KREAD(:,IKPT)=K_(:)
        CALL TRACE$PASS('MARKE 6')
!       == SET DATA FOR FURTHER USE ==================================
        NGL=GSET%NGL
        CALL PLANEWAVE$GETI4('NGG',NGG)
!       
!       ==============================================================
!       ==  DEFINE MAPPING OF THE ARRAYS                            ==
!       ==============================================================
        ALLOCATE(GVECL(3,NGL))
        ALLOCATE(GVECG(3,NGG))
        CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVECL)
        CALL PLANEWAVE$COLLECTR8(3,NGL,GVECL,NGG,GVECG)
        DEALLOCATE(GVECL)
        ALLOCATE(MAPG(NGG))
        IF(THISTASK.EQ.1) THEN
          CALL WAVES_MAPG(NGG_,GBAS_,GVECG_,NGG,GVECG,MAPG)
        END IF
        CALL MPE$BROADCAST(1,MAPG)
!!$DO IG=1,NGG
!!$  IF(MAPG(IG).NE.IG) THEN
!!$     PRINT*,'MAPG ',IKPT_,IG,MAPG(IG),NGG,NGG_
!!$  END IF
!!$ENDDO
        DEALLOCATE(GVECG_)
        DEALLOCATE(GVECG)
              CALL TRACE$PASS('MARKE 7')
!       
!       ==============================================================
!       ==  COLLECT DATA                                            ==
!       ==============================================================
IF(TSUPER_.AND.NDIM_.EQ.2) THEN
  CALL ERROR$MSG('NONCOLLINEAR WAVE FUNCTIONS AND SUPER WAVE FUNCTIONS')
  CALL ERROR$MSG('ARE NOT PROPERLY IMPLEMENTED BELOW')
  CALL ERROR$STOP('WAVES_READPSI')
END IF
        NBH=THIS%NBH
        NB=THIS%NB
        ALLOCATE(PSIIN(NGG_,NDIM_))
        ALLOCATE(PSIG(NGG,NDIM_))
        ALLOCATE(PSIL(NGL,NDIM_,NB_,NSPIN_))
        ALLOCATE(PSI(NGL,NDIM,NB,NSPIN))
        IF(TSUPER_) ALLOCATE(PSIINSUPER(NGG_,NDIM_))
        DO IWAVE=1,NWAVE
          DO ISPIN=1,NSPIN_
            TCHK=.FALSE.
            DO IB=1,NB_
              IF(THISTASK.EQ.1) THEN 
                IF(.NOT.TSUPER_) THEN
                  READ(NFIL,ERR=9999,IOSTAT=IOS)PSIIN !<<<<<<<<<<<<<<<<<<<
                ELSE  ! THIS FOR READING SUPER WAVE FUNCTIONS
                  IF(.NOT.TCHK) THEN   
                    READ(NFIL,ERR=9999,IOSTAT=IOS)PSIINSUPER !<<<<<<<<<<<<
                    DO IG=1,NGG_
                      IF(MINUSG(IG).LT.IG) CYCLE
                      DO IDIM=1,NDIM_
                        F1=PSIINSUPER(IG,IDIM)
                        F2=CONJG(PSIINSUPER(MINUSG(IG),IDIM))
                        PSIIN(IG,IDIM)=0.5D0*(F1+F2)
                        PSIIN(MINUSG(IG),IDIM)=CONJG(PSIIN(IG,IDIM))
                      ENDDO
                    ENDDO
                    TCHK=.TRUE.
                  ELSE
                    DO IG=1,NGG_
                      IF(MINUSG(IG).LT.IG) CYCLE
                      DO IDIM=1,NDIM_
                        F1=PSIINSUPER(IG,IDIM)
                        F2=CONJG(PSIINSUPER(MINUSG(IG),IDIM))
!this has been corrected 9.Nov2003: (- sign included). PEB
                        PSIIN(IG,IDIM)=-0.5D0*CI*(F1-F2)
                        PSIIN(MINUSG(IG),IDIM)=CONJG(PSIIN(IG,IDIM))
                      END DO
                    ENDDO
                    TCHK=.FALSE.
                  END IF
                END IF              
                DO IDIM=1,NDIM_
                  DO IG=1,NGG
                    IF(MAPG(IG).NE.0) THEN
                      PSIG(IG,IDIM)=PSIIN(MAPG(IG),IDIM)
                    ELSE
                      PSIG(IG,IDIM)=(0.D0,0.D0)
                    END IF
                  ENDDO
                ENDDO
              END IF
              DO IDIM=1,NDIM_
                CALL PLANEWAVE$DISTRIBUTEC8(1,NGG,PSIG(1,IDIM),NGL,PSIL(1,IDIM,IB,ISPIN))
              ENDDO
            ENDDO
          ENDDO
!       
!         ==============================================================
!         ==  TRANSFORM INTO THE CORRECT FORMAT                       ==
!         ==============================================================
!         == SAME REPRESENTATION =========================================
          CALL TRACE$PASS('MARKE 9')
          IF(NDIM.EQ.NDIM_.AND.NSPIN.EQ.NSPIN_) THEN
            DO IB=1,NB
              IF(IB.GT.NB_) EXIT
              PSI(:,:,IB,:)=PSIL(:,:,IB,:)
            ENDDO
!         == FROM NON-SPIN POLARIZED CALCULATION =========================
          ELSE IF(NSPIN_.EQ.1.AND.NDIM_.EQ.1) THEN
            IF(NSPIN.EQ.2) THEN
              DO ISPIN=1,NSPIN
                DO IB=1,NB
                  IF(IB.GT.NB_) EXIT
                  PSI(:,1,IB,ISPIN)=PSIL(:,1,IB,1)
                ENDDO
              ENDDO
            ELSE IF(NDIM.EQ.2) THEN
              DO IB=1,NB
                IB1=INT(0.5*REAL(IB+1)) ! IB1 =1,1,2,2,3,3,...
                IF(IB1.GT.NB_) EXIT
                IDIM=2+IB-2*IB1         ! IDIM=1,2,1,2,1,2,...
                PSI(:,IDIM,IB,1)=PSIL(:,1,IB1,1)
              ENDDO
            ELSE
              CALL ERROR$MSG('TRANSFORMATION NOT IMPLEMENTED')
              CALL ERROR$STOP('WAVES$READPSI')
            END IF
!         == FROM SPIN POLARIZED CALCULATION ===========================
          ELSE IF(NSPIN_.EQ.2.AND.NDIM_.EQ.1) THEN
            IF(NSPIN.EQ.1.AND.NDIM.EQ.1) THEN
              DO IB=1,NB
                PSI(:,1,IB,1)=PSIL(:,1,IB,1)
              ENDDO
            ELSE IF(NSPIN.EQ.1.AND.NDIM.EQ.2) THEN
              DO IB=1,(NB+1)/2
                IB1=2*IB-1
                IB2=2*IB
                PSI(:,1,IB1,1)=PSIL(:,1,IB,1)
                PSI(:,2,IB1,1)=(0.D0,0.D0)
                IF(IB2.GT.NB) EXIT
                PSI(:,1,IB2,1)=(0.D0,0.D0)
                PSI(:,2,IB2,1)=PSIL(:,1,IB,2)
              ENDDO
            ELSE
              CALL ERROR$MSG('TRANSFORMATION NOT IMPLEMENTED')
              CALL ERROR$STOP('WAVES$READPSI')
            END IF
!         == FROM NONCOLLINEAR CALCULATION  ==============================
          ELSE IF(NSPIN_.EQ.2.AND.NDIM_.EQ.1) THEN
            IF(NSPIN.EQ.1.AND.NDIM.EQ.1) THEN
              DO IB=1,NB 
                PSI(:,1,IB,1)=PSIL(:,1,2*IB-1,1)
              ENDDO
            ELSE IF(NSPIN.EQ.1.AND.NDIM.EQ.1) THEN
              DO IB=1,NB 
                IB1=2*IB-1
                IB2=2*IB
                PSI(:,1,IB,1)=PSIL(:,1,IB1,1)
                PSI(:,1,IB,2)=PSIL(:,1,IB2,1)
              ENDDO
            ELSE
              CALL ERROR$MSG('TRANSFORMATION NOT IMPLEMENTED')
              CALL ERROR$STOP('WAVES$READPSI')
            END IF
          END IF           
          CALL TRACE$PASS('MARKE 10')
!       
!         ==============================================================
!         ==  MAP ONTO SUPER WAVE FUNCTIONS                           ==
!         ==============================================================
          CALL PLANEWAVE$GETL4('TINV',TSUPER)
          IF(TSUPER) THEN
            DO ISPIN=1,NSPIN
!             == remove the factor exp(i*phi) from the wave function
              do ib=1,nb
                call PLANEWAVE$SCALARPRODUCT('-',NGL,1,1 &
        &                   ,psi(:,1,ib,ispin),1,psi(:,1,ib,ispin),cmat)
                csvar=cmat(1,1)
                csvar=conjg(sqrt(csvar/sqrt(csvar*conjg(csvar))))
                psi(:,1,ib,ispin)=psi(:,1,ib,ispin)*csvar
              enddo
              DO IB=1,NBH
                IB1=2*IB-1
                IB2=2*IB
                PSI(:,1,IB,ISPIN)=PSI(:,1,IB1,ISPIN)+CI*PSI(:,1,IB2,ISPIN)
              ENDDO
            ENDDO
          END IF
!       
!         ==============================================================
!         ==  MAP BACK                                                ==
!         ==============================================================
          CALL TRACE$PASS('MARKE 11')
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            IF(IWAVE.EQ.1) THEN
              DO IB=1,NBH
                DO IDIM=1,NDIM
                  DO IG=1,NGL
                    THIS%PSI0(IG,IDIM,IB)=PSI(IG,IDIM,IB,ISPIN)
                  ENDDO
                ENDDO
              ENDDO
              IF(NWAVE.EQ.1) THEN
                THIS%PSIM=THIS%PSI0
              END IF
            ELSE IF(IWAVE.EQ.2) THEN
              DO IB=1,NBH
                DO IDIM=1,NDIM
                  DO IG=1,NGL
                    THIS%PSIM(IG,IDIM,IB)=PSI(IG,IDIM,IB,ISPIN)
                  ENDDO
                ENDDO
              ENDDO
            END IF
          ENDDO
        ENDDO 
        DEALLOCATE(PSIL)
        DEALLOCATE(PSIG)
        DEALLOCATE(MAPG)
        DEALLOCATE(PSI)
        DEALLOCATE(PSIIN)
        IF(ALLOCATED(PSIINSUPER))DEALLOCATE(PSIINSUPER)
        IF(ALLOCATED(MINUSG)) DEALLOCATE(MINUSG)
        CALL WAVES_READLAMBDA(NFIL,IKPT,.FALSE.)
      ENDDO
!
!     ==================================================================
!     ==  COMPLETE K-POINTS                                           ==
!     ==================================================================
      CALL TRACE$PASS('MARKE 14')
      DO IKPT=1,NKPT
        IF(TREAD(IKPT)) CYCLE
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETR8A('K',3,K)
        SVAR=1.D+10
        DO IKPT_=1,NKPT_
          IF(.NOT.TREAD(IKPT_)) CYCLE
          SVAR1=(K(1)-KREAD(1,IKPT_))**2+(K(2)-KREAD(2,IKPT_))**2 &
     &                                  +(K(3)-KREAD(3,IKPT_))**2
          IF(SVAR1.LT.SVAR) THEN
            IKPT0=IKPT_
            SVAR=SVAR1
          END IF
        ENDDO       
        CALL WAVES_COPYPSI(IKPT,IKPT0)
        CALL WAVES_COPYLAMBDA(IKPT0,IKPT)
      ENDDO 

!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
                               CALL TRACE$POP
      RETURN
 9999 CONTINUE
      CALL ERROR$MSG('ERROR WHILE READING FROM FILE')
      CALL ERROR$I4VAL('IOS',IOS)
      CALL ERROR$STOP('WAVES_READPSI')
      END
!
!     ..................................................................
      SUBROUTINE WAVES_WRITELAMBDA(NFIL,IKPT,NREC)
!     ******************************************************************
!     **                                                              **
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                        **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NFIL
      INTEGER(4)   ,INTENT(IN) :: IKPT
      INTEGER(4)   ,INTENT(INOUT) :: NREC
      INTEGER(4)               :: NTASKS,THISTASK
      INTEGER(4)               :: ISPIN,I
      INTEGER(4)               :: NB
      CHARACTER(8)             :: KEY
      COMPLEX(8)   ,ALLOCATABLE:: LAMBDA(:,:,:)
      INTEGER(4)               :: NLAMBDA
!     ******************************************************************
      CALL MPE$QUERY(NTASKS,THISTASK)
      CALL WAVES_SELECTWV(IKPT,1)
!
!     ==================================================================
!     == COUNT DEPTH OF HISTORY FOR LAGRANGE PARAMETERS               ==
!     ==================================================================
      NLAMBDA=4
      IF(.NOT.ASSOCIATED(THIS%RLAM3M)) NLAMBDA=3
      IF(.NOT.ASSOCIATED(THIS%RLAM2M)) NLAMBDA=2
      IF(.NOT.ASSOCIATED(THIS%RLAMM)) NLAMBDA=1
      IF(.NOT.ASSOCIATED(THIS%RLAM0)) NLAMBDA=0
      IF(RSTRTTYPE.EQ.'STATIC') NLAMBDA=1
!
!     ==================================================================
!     == RETURN #(RECORDS) OR IF ON TASK OTHER THAN THE FIRST         ==
!     ==================================================================
      IF(NREC.EQ.-1) THEN
        NREC=1+NSPIN*NLAMBDA
        RETURN
      END IF
      IF(THISTASK.NE.1) RETURN
      NREC=0
!
!     ==================================================================
!     == NOW WRITE TO FILE                                            ==
!     ==================================================================
      KEY='LAMBDA'
      NB=THIS%NB
      WRITE(NFIL)KEY,NB,NDIM,NSPIN,NLAMBDA  !<<<<<<<<<<<<<<<<<<<<<<<<<<<
      NREC=NREC+1
      IF(NLAMBDA.GE.1) THEN
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          WRITE(NFIL)THIS%RLAM0 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC=NREC+1
        ENDDO
      END IF
      IF(NLAMBDA.GE.2) THEN
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          WRITE(NFIL)THIS%RLAMM !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC=NREC+1
        ENDDO
      END IF
      IF(NLAMBDA.GE.3) THEN
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          WRITE(NFIL)THIS%RLAM2M !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC=NREC+1
        ENDDO
      END IF
      IF(NLAMBDA.GE.4) THEN
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          WRITE(NFIL)THIS%RLAM3M !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC=NREC+1
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_READLAMBDA(NFIL,IKPT,TSKIP)
!     ******************************************************************
!     **                                                              **
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                        **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NFIL
      INTEGER(4)   ,INTENT(IN) :: IKPT
      LOGICAL(4)   ,INTENT(IN) :: TSKIP
      INTEGER(4)               :: NTASKS,THISTASK
      INTEGER(4)               :: ISPIN,I,IB1,IB2
      INTEGER(4)               :: NB_,NDIM_,NSPIN_
      INTEGER(4)               :: NB,NBA
      INTEGER(4)               :: NLAMBDA
      CHARACTER(8)             :: KEY
      COMPLEX(8)   ,ALLOCATABLE:: LAMBDA(:,:,:)
      COMPLEX(8)   ,ALLOCATABLE:: LAMBDA2(:,:,:)
!     ******************************************************************
                           CALL TRACE$PUSH('WAVES_READLAMBDA')
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     == READ AND BROADCAST HISTORY-DEPTH OF LAGRANGE PARAMETERS      ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        READ(NFIL)KEY,NB_,NDIM_,NSPIN_,NLAMBDA  !<<<<<<<<<<<<<<<<<<<<<<<<<<
        IF(KEY.NE.'LAMBDA') THEN
          CALL ERROR$MSG('ID IS NOT "LAMBDA"')
          CALL ERROR$MSG('FILE IS CORRUPTED')
          CALL ERROR$STOP('WAVES_READLAMBDA')
        END IF
      END IF
      IF(TSKIP) THEN
        DO I=1,NLAMBDA
          DO ISPIN=1,NSPIN_
           IF(THISTASK.EQ.1)READ(NFIL)
          ENDDO
        ENDDO
        CALL TRACE$POP
        RETURN
      END IF
      CALL MPE$BROADCAST(1,NLAMBDA)
!
!     ==================================================================
!     ==  READ AND FOLD DOWN                                          ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        CALL WAVES_SELECTWV(IKPT,1)
        NB=THIS%NB
        ALLOCATE(LAMBDA(NB_,NB_,NSPIN_))
        ALLOCATE(LAMBDA2(NB,NB,NSPIN))
        NBA=MIN(NB,NB_)
        DO I=1,NLAMBDA
          DO ISPIN=1,NSPIN_
            READ(NFIL)LAMBDA(:,:,ISPIN)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          ENDDO
!
!         ==============================================================
!         ==  TRANSFORM BETWEEN DATA MODELS                           ==
!         ==============================================================
          LAMBDA2=(0.D0,0.D0)
          IF(NSPIN.EQ.NSPIN_.AND.NDIM.EQ.NDIM_) THEN
            LAMBDA2(1:NBA,1:NBA,:)=LAMBDA(1:NBA,1:NBA,:)
          ELSE
!           == FROM NON-SPIN POLARIZED CALCULATION =====================
            IF(NDIM_.EQ.1.AND.NSPIN_.EQ.1) THEN
              IF(NDIM.EQ.1.AND.NSPIN.EQ.2) THEN
                LAMBDA2(1:NBA,1:NBA,1)=LAMBDA(1:NBA,1:NBA,1)
                LAMBDA2(1:NBA,1:NBA,2)=LAMBDA(1:NBA,1:NBA,1)
              ELSE IF(NDIM.EQ.2.AND.NSPIN.EQ.1) THEN
                DO IB1=1,NB_
                  DO IB2=1,NB_
                    IF(2*IB1-1.GT.NB.OR.2*IB2-1.GT.NB)CYCLE
                    LAMBDA2(2*IB1-1,2*IB2-1,1)=LAMBDA(IB1,IB2,1)
                    IF(2*IB1.GT.NB.OR.2*IB2.GT.NB)CYCLE
                    LAMBDA2(2*IB1  ,2*IB2  ,1)=LAMBDA(IB1,IB2,1)
                  ENDDO
                ENDDO
              ELSE
                CALL ERROR$MSG('INVALID OPTION TRANSFORMING FROM')
                CALL ERROR$MSG('NON-SPINPOLARIZED DATA')
                CALL ERROR$I4VAL('NDIM',NDIM)
                CALL ERROR$I4VAL('NSPIN',NSPIN)
                CALL ERROR$I4VAL('NDIM_',NDIM_)
                CALL ERROR$I4VAL('NSPIN_',NSPIN_)
                CALL ERROR$STOP('WAVES_READLAMBDA')
              END IF
!           == FROM SPIN POLARIZED CALCULATION =========================
            ELSE IF(NDIM_.EQ.1.AND.NSPIN_.EQ.2) THEN
              IF(NDIM.EQ.1.AND.NSPIN.EQ.1) THEN
                LAMBDA2(1:NBA,1:NBA,1)=LAMBDA(1:NBA,1:NBA,1)
              ELSE IF(NDIM.EQ.2.AND.NSPIN.EQ.1) THEN
                DO IB1=1,NB_
                  DO IB2=1,NB_
                    IF(2*IB1-1.GT.NB.OR.2*IB2-1.GT.NB)CYCLE
                    LAMBDA2(2*IB1-1,2*IB2-1,1)=LAMBDA(IB1,IB2,1)
                    IF(2*IB1.GT.NB.OR.2*IB2.GT.NB)CYCLE
                    LAMBDA2(2*IB1  ,2*IB2  ,1)=LAMBDA(IB1,IB2,2)
                  ENDDO
                ENDDO
              ELSE
                CALL ERROR$MSG('INVALID OPTION TRANSFORMING FROM')
                CALL ERROR$MSG('SPINPOLARIZED DATA')
                CALL ERROR$I4VAL('NDIM',NDIM)
                CALL ERROR$I4VAL('NSPIN',NSPIN)
                CALL ERROR$I4VAL('NDIM_',NDIM_)
                CALL ERROR$I4VAL('NSPIN_',NSPIN_)
                CALL ERROR$STOP('WAVES_READLAMBDA')
              END IF
!           == FROM NONCOLLINEAR WAVE FUNCTION =========================
            ELSE IF(NDIM_.EQ.2.AND.NSPIN_.EQ.1) THEN
              IF(NDIM.EQ.1.AND.NSPIN.EQ.1) THEN
                DO IB1=1,NB
                  DO IB2=1,NB
                    IF(2*IB1-1.GT.NB_.OR.2*IB2-1.GT.NB_)CYCLE
                    LAMBDA2(IB1,IB2,1)=LAMBDA(2*IB1-1,2*IB2-1,1)
                  ENDDO
                ENDDO
              ELSE IF(NDIM.EQ.2.AND.NSPIN.EQ.1) THEN
                DO IB1=1,NB
                  DO IB2=1,NB
                    IF(2*IB1-1.GT.NB_.OR.2*IB2-1.GT.NB_)CYCLE
                    LAMBDA2(IB1,IB2,1)=LAMBDA(2*IB1-1,2*IB2-1,1)
                    IF(2*IB1.GT.NB_.OR.2*IB2.GT.NB_)CYCLE
                    LAMBDA2(IB1,IB2,2)=LAMBDA(2*IB1,2*IB2,1)
                  ENDDO
                ENDDO
              ELSE
                CALL ERROR$MSG('INVALID OPTION TRANSFORMING FROM')
                CALL ERROR$MSG('NON-COLLINEAR DATA')
                CALL ERROR$I4VAL('NDIM',NDIM)
                CALL ERROR$I4VAL('NSPIN',NSPIN)
                CALL ERROR$I4VAL('NDIM_',NDIM_)
                CALL ERROR$I4VAL('NSPIN_',NSPIN_)
                CALL ERROR$STOP('WAVES_READLAMBDA')
              END IF
            ELSE
              CALL ERROR$MSG('INVALID OPTION')
              CALL ERROR$I4VAL('NDIM',NDIM)
              CALL ERROR$I4VAL('NSPIN',NSPIN)
              CALL ERROR$I4VAL('NDIM_',NDIM_)
              CALL ERROR$I4VAL('NSPIN_',NSPIN_)
              CALL ERROR$STOP('WAVES_READLAMBDA')
            END IF
          END IF
!
!         ==============================================================
!         ==  MAP ONTO THIS                                           ==
!         ==============================================================
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            NB=THIS%NB
            IF(I.EQ.1) THEN
              IF(.NOT.ASSOCIATED(THIS%RLAM0))ALLOCATE(THIS%RLAM0(NB,NB))
              THIS%RLAM0=LAMBDA2(:,:,ISPIN)
            ELSE IF(I.EQ.2) THEN
              IF(.NOT.ASSOCIATED(THIS%RLAMM))ALLOCATE(THIS%RLAMM(NB,NB))
              THIS%RLAMM=LAMBDA2(:,:,ISPIN)
            ELSE IF(I.EQ.3) THEN 
              IF(.NOT.ASSOCIATED(THIS%RLAM2M))ALLOCATE(THIS%RLAM2M(NB,NB))
              THIS%RLAM2M=LAMBDA2(:,:,ISPIN)
            ELSE IF(I.EQ.4) THEN
              IF(.NOT.ASSOCIATED(THIS%RLAM3M))ALLOCATE(THIS%RLAM3M(NB,NB))
              THIS%RLAM3M=LAMBDA2(:,:,ISPIN)
            END IF
          ENDDO
        ENDDO
        DEALLOCATE(LAMBDA)
        DEALLOCATE(LAMBDA2)
      END IF
!
!     ==================================================================
!     ==  BROADCAST LAMBDA                                            ==
!     ==================================================================
      DO ISPIN=1,NSPIN
        CALL WAVES_SELECTWV(IKPT,ISPIN)
        NB=THIS%NB
        IF(NLAMBDA.GE.1)THEN
          IF(.NOT.ASSOCIATED(THIS%RLAM0))ALLOCATE(THIS%RLAM0(NB,NB))
          CALL MPE$BROADCAST(1,THIS%RLAM0)
        END IF
        IF(NLAMBDA.GE.2)THEN
          IF(.NOT.ASSOCIATED(THIS%RLAMM))ALLOCATE(THIS%RLAMM(NB,NB))
          CALL MPE$BROADCAST(1,THIS%RLAMM)
        END IF
        IF(NLAMBDA.GE.3)THEN
          IF(.NOT.ASSOCIATED(THIS%RLAM2M))ALLOCATE(THIS%RLAM2M(NB,NB))
          CALL MPE$BROADCAST(1,THIS%RLAM2M)
        END IF
        IF(NLAMBDA.GE.4)THEN
          IF(.NOT.ASSOCIATED(THIS%RLAM3M))ALLOCATE(THIS%RLAM3M(NB,NB))
          CALL MPE$BROADCAST(1,THIS%RLAM3M)
        END IF
      ENDDO
                           CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_MAPG(NGA,GBASA,GVECA,NGB,GVECB,MAP)
!     ******************************************************************
!     **  ARRAY TO MAP WAVE FUNCTIONS READ FROM FILE INTO             **
!     **  AN EXISTING GRID                                            **
!     **  NGA,GBASA,GVECA DESCRIBE THE ARRAY READ FROM FILE           **
!     **  NGB,GVECB DESCRIBE THE EXISTING G-ARRAY                     **
!     **                                                              **
!     **  DO IG=1,NG                                                  **
!     **    IF(MAP(IG).EQ.0) THEN                                     **
!     **      PSI(IG)=(0.D0,0.D0)                                     **
!     **    ELSE                                                      **
!     **      PSI(IG)=PSIIN(MAP(IG))                                  **
!     **    END IF                                                    **
!     **  ENDDO                                                       **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NGA
      REAL(8)   ,INTENT(IN) :: GBASA(3,3)
      REAL(8)   ,INTENT(IN) :: GVECA(3,NGA)
      INTEGER(4),INTENT(IN) :: NGB
      REAL(8)   ,INTENT(IN) :: GVECB(3,NGB)
      INTEGER(4),INTENT(OUT):: MAP(NGB)
      INTEGER(4)            :: IG,I
      REAL(8)   ,ALLOCATABLE:: MAP3D(:,:,:)
      REAL(8)               :: GBASIN(3,3)
      INTEGER(4)            :: IVEC1(3),IVECA(3,NGA)
      REAL(8)               :: KA(3)
      REAL(8)               :: G(3)
      INTEGER(4)            :: MING(3),MAXG(3)
logical(4):: tchk
!     ******************************************************************
      CALL LIB$INVERTR8(3,GBASA,GBASIN)
!
!     ==================================================================
!     == DETERMINE K-POINT                                            ==
!     ==================================================================
      G(:)=MATMUL(GBASIN,GVECA(:,1))
      IVEC1=NINT(G)
      DO I=1,3
        IVEC1(I)=NINT(G(I))
      ENDDO
      G=REAL(NINT(G))
      G(:)=MATMUL(GBASA,REAL(IVEC1,KIND=8))
      KA(:)=GVECA(:,1)-G(:)
!
!     ==================================================================
!     == DIMENSION AND ALLOCATE 3-D GRID                              ==
!     ==================================================================
      MING=0
      MAXG=0
      DO IG=1,NGA
        IVECA(:,IG)=NINT(MATMUL(GBASIN,GVECA(:,IG)-KA(:)))
        DO I=1,3
          MING(I)=MIN(MING(I),IVECA(I,IG))
          MAXG(I)=MAX(MAXG(I),IVECA(I,IG))
        ENDDO
      ENDDO
      ALLOCATE(MAP3D(MING(1):MAXG(1),MING(2):MAXG(2),MING(3):MAXG(3)))
!
!     ==================================================================
!     == MAP FIRST ARRAY ONTO 3-D GRID                                ==
!     ==================================================================
      MAP3D(:,:,:)=0
      DO IG=1,NGA
        IF(MAP3D(IVECA(1,IG),IVECA(2,IG),IVECA(3,IG)).EQ.0) THEN
          MAP3D(IVECA(1,IG),IVECA(2,IG),IVECA(3,IG))=IG
        ELSE
          CALL ERROR$MSG('TWO G-VECTORS ARE MAPPED ONTO THE SAME POINT')
          CALL ERROR$STOP('WAVES_MAPG')
        END IF
      ENDDO
!
!     ==================================================================
!     ==  FOR EACH VECTOR (B) PICK OUT THE CLOSEST GRID (A) POINT     ==
!     ==================================================================
      LOOP1:DO IG=1,NGB
        MAP(IG)=0
        G(:)=MATMUL(GBASIN,GVECB(:,IG)-KA(:))
        DO I=1,3
          IVEC1(I)=NINT(G(I))
          IF(IVEC1(I).LT.MING(I).OR.IVEC1(I).GT.MAXG(I)) CYCLE LOOP1
        ENDDO
        MAP(IG)=MAP3D(IVEC1(1),IVEC1(2),IVEC1(3))
      ENDDO LOOP1
      DEALLOCATE(MAP3D)
!
!     ==================================================================
!     ==  print if remapping has been done                            ==
!     ==================================================================
      tchk=.false.
      do ig=1,ngb
        tchk=tchk.or.(map(ig).ne.ig)
      end do
      IF(TCHK)PRINT*,'ORDER OF G-VECTORS FROM RESTART FILE HAS BEEN CHANGED'
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_COPYPSI(IKPT1,IKPT2)
!     ******************************************************************
!     **                                                              **
!     **  MAP WAVE FUNCTIONS FROM ONE K-POINT IKPT2 ONTO IKPT1        **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IKPT1
      INTEGER(4),INTENT(IN) :: IKPT2
      INTEGER(4)            :: NGL1,NGL2,NGG1,NGG2
      REAL(8)   ,ALLOCATABLE:: GVECL(:,:)
      REAL(8)   ,ALLOCATABLE:: GVECG1(:,:)
      REAL(8)   ,ALLOCATABLE:: GVECG2(:,:)
      INTEGER(4),ALLOCATABLE:: MAPG(:)
      LOGICAL(4)            :: TSUPER1,TSUPER2
      INTEGER(4)            :: NB1,NB2,NBH1,NBH2,NBN,NBHN
      COMPLEX(8),ALLOCATABLE:: PSI1(:,:,:)
      COMPLEX(8),ALLOCATABLE:: PSI2(:,:,:)
      COMPLEX(8),ALLOCATABLE:: PSIG1(:)
      COMPLEX(8),ALLOCATABLE:: PSIG2(:)
      COMPLEX(8),ALLOCATABLE:: PSITMP2(:,:)
      INTEGER(4)            :: IWAVE,ISPIN,IB,IB1,IB2,IDIM,IG
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)               :: GBAS(3,3)
      INTEGER(4)            :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     ==  MAPPING FROM IKPT2 TO IKPT1                                 ==
!     ==================================================================
!     == GET GLOBAL GVECTORS FROM IKPT1 ================================
      CALL WAVES_SELECTWV(IKPT1,1)
      CALL PLANEWAVE$SELECT(GSET%ID)
      NGL1=GSET%NGL
      ALLOCATE(GVECL(3,NGL1))
      CALL PLANEWAVE$GETR8A('GVEC',3*NGL1,GVECL)
      CALL PLANEWAVE$GETI4('NGG',NGG1)
      IF(THISTASK.NE.1) NGG1=1 
      ALLOCATE(GVECG1(3,NGG1))
      CALL PLANEWAVE$COLLECTR8(3,NGL1,GVECL,NGG1,GVECG1)
      DEALLOCATE(GVECL)
!
!     == GET GLOBAL GVECTORS FROM IKPT2 ================================
      CALL WAVES_SELECTWV(IKPT2,1)
      CALL PLANEWAVE$SELECT(GSET%ID)
      NGL2=GSET%NGL
      ALLOCATE(GVECL(3,NGL2))
      CALL PLANEWAVE$GETR8A('GVEC',3*NGL2,GVECL)
      CALL PLANEWAVE$GETI4('NGG',NGG2)
      IF(THISTASK.NE.1) NGG2=1 
      ALLOCATE(GVECG2(3,NGG2))
      CALL PLANEWAVE$COLLECTR8(3,NGL2,GVECL,NGG2,GVECG2)
      DEALLOCATE(GVECL)
!
!     == COPY GBAS      ================================================
      CALL WAVES_SELECTWV(IKPT2,1)
      CALL PLANEWAVE$SELECT(GSET%ID)
      CALL PLANEWAVE$GETR8A('GBAS',9,GBAS)
      CALL WAVES_SELECTWV(IKPT2,1)
      CALL PLANEWAVE$SELECT(GSET%ID)
      CALL PLANEWAVE$SETR8A('GBAS',9,GBAS)
!
!     == DEFINE MAPPING ================================================
      IF(THISTASK.EQ.1) THEN
        ALLOCATE(MAPG(NGG1))
        CALL WAVES_MAPG(NGG2,GBAS,GVECG2,NGG1,GVECG1,MAPG)
      END IF
      DEALLOCATE(GVECG1)
      DEALLOCATE(GVECG2)
!
!     ==================================================================
!     ==  COPY                                                        ==
!     ==================================================================
      CALL WAVES_SELECTWV(IKPT1,1)
      CALL PLANEWAVE$SELECT(GSET%ID)
      CALL PLANEWAVE$GETL4('TINV',TSUPER1)
      NB1=THIS%NB
      NBH1=THIS%NBH
      CALL WAVES_SELECTWV(IKPT2,1)
      CALL PLANEWAVE$SELECT(GSET%ID)
      CALL PLANEWAVE$GETL4('TINV',TSUPER2)
      NB2=THIS%NB
      NBH2=THIS%NBH

      NBN=MIN(NB1,NB2)
      IF(TSUPER1) THEN
        NBHN=MIN(NBN/2,NBH1)
      ELSE
        NBHN=NBN
      END IF
      ALLOCATE(PSI1(NGL1,NDIM,NB1))
      ALLOCATE(PSI2(NGL2,NDIM,NB2))
      ALLOCATE(PSIG1(NGG1))
      ALLOCATE(PSIG2(NGG2))
      ALLOCATE(PSITMP2(NGL2,2))
!
!     ==================================================================
!     == MAP DATA ONTO TEMP ARRAY                                     ==
!     ==================================================================
      DO IWAVE=1,2
        DO ISPIN=1,NSPIN
!
!         ==============================================================
!         == MAP DATA ONTO TEMP ARRAY                                 ==
!         ==============================================================
          CALL WAVES_SELECTWV(IKPT2,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(IWAVE.EQ.1) THEN
            DO IB=1,NBH2
              PSI2(:,:,IB)=THIS%PSI0(:,:,IB)
            ENDDO
          ELSE IF(IWAVE.EQ.2) THEN
            DO IB=1,NBH2
              PSI2(:,:,IB)=THIS%PSIM(:,:,IB)
            ENDDO
          END IF
!
!         ==============================================================
!         ==  EXPAND TO REGULAR WAVE FUNCTIONS                        ==
!         ==============================================================
          IF(TSUPER2) THEN
            DO IB=NBH2,1,-1
              IB1=2*IB-1
              IB2=2*IB
              DO IDIM=1,NDIM
                PSITMP2(:,1)=PSI2(:,IDIM,IB)
                CALL PLANEWAVE$INVERTG(NGL2,PSI2(1,IDIM,IB),PSITMP2(1,2))
                PSI2(:,IDIM,IB1)= 0.5D0   *(PSITMP2(:,1)+PSITMP2(:,2))
                PSI2(:,IDIM,IB2)=-0.5D0*CI*(PSITMP2(:,1)-PSITMP2(:,2))
              END DO
            ENDDO
          END IF
!
!         ==============================================================
!         ==  MAP ONTO IKPT1                                          ==
!         ==============================================================
          DO IB=1,NBN
            DO IDIM=1,NDIM
              CALL WAVES_SELECTWV(IKPT2,ISPIN)
              CALL PLANEWAVE$SELECT(GSET%ID)
              CALL PLANEWAVE$COLLECTC8(1,NGL2,PSI2(1,IDIM,IB),NGG2,PSIG2)
              IF(THISTASK.EQ.1) THEN
                DO IG=1,NGG1
                  IF(MAPG(IG).NE.0) THEN
                    PSIG1(IG)=PSIG2(MAPG(IG))
                  ELSE
                    PSIG1(IG)=(0.D0,0.D0)
                  END IF
                ENDDO
              END IF
              CALL WAVES_SELECTWV(IKPT1,ISPIN)
              CALL PLANEWAVE$SELECT(GSET%ID)
              CALL PLANEWAVE$DISTRIBUTEC8(1,NGG1,PSIG1,NGL1,PSI1(1,IDIM,IB))
            ENDDO
          ENDDO
!
!         ==============================================================
!         == MAP ONTO SUPER WAVE FUNCTIONS                            ==
!         ==============================================================
          IF(TSUPER1) THEN
            DO IB=1,NBHN
              IB1=2*IB-1
              IB2=2*IB
              PSI1(:,1,IB)=PSI1(:,1,IB1)+CI*PSI1(:,1,IB2)
            ENDDO
          END IF
!
!         ==============================================================
!         == SAVE INTO THIS                                           ==
!         ==============================================================
          CALL WAVES_SELECTWV(IKPT1,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(IWAVE.EQ.1) THEN
            DO IB=1,NBHN
              THIS%PSI0(:,:,IB)=PSI1(:,:,IB)
            ENDDO
          ELSE IF(IWAVE.EQ.2) THEN
            DO IB=1,NBHN
              THIS%PSIM(:,:,IB)=PSI1(:,:,IB)
            ENDDO
          END IF
        ENDDO
      ENDDO
      IF(THISTASK.EQ.1)DEALLOCATE(MAPG)
      DEALLOCATE(PSI1)
      DEALLOCATE(PSI2)
      DEALLOCATE(PSIG1)
      DEALLOCATE(PSIG2)
      DEALLOCATE(PSITMP2)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_COPYLAMBDA(IKPT1,IKPT2)
!     ******************************************************************
!     **                                                              **
!     **  COPY LAGRANGE MULTIPLIERS FROM K-POINT IKPT1 TO IKPT2       **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: IKPT1
      INTEGER(4)   ,INTENT(IN) :: IKPT2
      INTEGER(4)               :: NB1,NB2,NBA
      INTEGER(4)               :: ISPIN
      COMPLEX(8)   ,ALLOCATABLE:: LAMBDA(:,:)
!     ******************************************************************
      CALL WAVES_SELECTWV(IKPT1,1)
      NB1=THIS%NB
      CALL WAVES_SELECTWV(IKPT2,1)
      NB2=THIS%NB
      NBA=MIN(NB1,NB2)
      ALLOCATE(LAMBDA(NBA,NBA))
      DO ISPIN=1,NSPIN
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAM0)) RETURN
        LAMBDA=THIS%RLAM0(1:NBA,1:NBA)
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAM0)) ALLOCATE(THIS%RLAM0(NB2,NB2))
        THIS%RLAM0=(0.D0,0.D0)        
        THIS%RLAM0(1:NBA,1:NBA)=LAMBDA
      ENDDO
      DO ISPIN=1,NSPIN
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAMM)) RETURN
        LAMBDA=THIS%RLAMM(1:NBA,1:NBA)
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAMM))ALLOCATE(THIS%RLAMM(NB2,NB2))
        THIS%RLAMM=(0.D0,0.D0)        
        THIS%RLAMM(1:NBA,1:NBA)=LAMBDA
      ENDDO
      DO ISPIN=1,NSPIN
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAM2M)) RETURN
        LAMBDA=THIS%RLAM2M(1:NBA,1:NBA)
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAM2M)) ALLOCATE(THIS%RLAM2M(NB2,NB2))
        THIS%RLAM2M=(0.D0,0.D0)        
        THIS%RLAM2M(1:NBA,1:NBA)=LAMBDA
      ENDDO
      DO ISPIN=1,NSPIN
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAM3M)) RETURN
        LAMBDA=THIS%RLAM3M(1:NBA,1:NBA)
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAM3M)) ALLOCATE(THIS%RLAM3M(NB2,NB2))
        THIS%RLAM2M=(0.D0,0.D0)        
        THIS%RLAM3M(1:NBA,1:NBA)=LAMBDA
      ENDDO
      DEALLOCATE(LAMBDA)
      RETURN
      END SUBROUTINE WAVES_COPYLAMBDA!
!      .................................................................
       SUBROUTINE WAVES_SPINOROVERLAP(NBH,NB,IKPT,QMAT)
!      *****************************************************************
!      **  CALCULATES Q_K,I,L,J=1/2*<PSI_I,K|PSI_J,L>                 **
!      **  I AND J ARE BAND-INDICES                                   **
!      **  K AND L ARE 1 OR 2 AND ARE THE SPINOR PART OF THE WAVEFUNCTION
!      *****************************************************************
       USE WAVES_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NBH,NB,IKPT
       COMPLEX(8),INTENT(OUT):: QMAT(2*NB*NSPIN,2*NB*NSPIN)
       COMPLEX(8),ALLOCATABLE :: AUXMAT(:,:)
       COMPLEX(8),ALLOCATABLE :: AUXMAT2(:,:)
       COMPLEX(8),ALLOCATABLE :: OPROJ(:,:,:)
       COMPLEX(8),ALLOCATABLE :: PSI0(:,:,:)
       INTEGER(4)            :: NGL,NBD,I,J,NDIMHALF,ISPIN
       IF (NDIM.EQ.1.AND.NSPIN.EQ.1) THEN
          CALL ERROR$MSG('S^2 ONLY POSSIBLE FOR NOT SPIN RESTRICTED CALCULATION')
          CALL ERROR$I4VAL('NDIM',NDIM)
          CALL ERROR$I4VAL('NSPIN',NSPIN)
          CALL ERROR$STOP('WAVES_SPINOROVERLAP')
       END IF
       IF(NDIM.EQ.2) THEN   !NON-COLLINEAR
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)      
          NGL=GSET%NGL   
          NBD=2*NB
          NDIMHALF=1
          ALLOCATE(AUXMAT(NBD,NBD))
          CALL WAVES_1COVERLAP(.TRUE.,MAP,NDIMHALF,NBD,NBD,MAP%NPRO,THIS%PROJ,THIS%PROJ,AUXMAT)
          CALL WAVES_OVERLAP(.TRUE.,NGL,NDIMHALF,NBD,NBD,THIS%PSI0,THIS%PSI0,QMAT)
          DO J=1,NBD
             DO I=1,NBD
                QMAT(I,J)=(QMAT(I,J)+AUXMAT(I,J))*0.5D0
             ENDDO
          ENDDO
          DEALLOCATE(AUXMAT)
       ELSE  ! COLLINEAR SPIN POLARIZED
          ALLOCATE(AUXMAT(NB*2,NB*2))
          DO ISPIN=1,NSPIN
             CALL WAVES_SELECTWV(IKPT,ISPIN)
             CALL PLANEWAVE$SELECT(GSET%ID)      
             NGL=GSET%NGL              
             NBD=2*NB
             IF(ISPIN.EQ.1) ALLOCATE(OPROJ(1,NBH*2,MAP%NPRO))
             IF(ISPIN.EQ.1) THEN
                OPROJ(:,1:NBH,:)=THIS%PROJ
             ELSE
                OPROJ(:,NBH+1:2*NBH,:)=THIS%PROJ
             END IF
          END DO
          CALL WAVES_1COVERLAP(.TRUE.,MAP,1,NBH*2,NB*2,MAP%NPRO,OPROJ,OPROJ,AUXMAT)
          DEALLOCATE(OPROJ)
          ALLOCATE(AUXMAT2(NB*2,NB*2))
          DO ISPIN=1,NSPIN
             CALL WAVES_SELECTWV(IKPT,ISPIN)
             CALL PLANEWAVE$SELECT(GSET%ID)      
             NGL=GSET%NGL
             IF(ISPIN.EQ.1) ALLOCATE(PSI0(NGL,1,2*NBH))
             IF(ISPIN.EQ.1) THEN
                PSI0(:,:,1:NBH)=THIS%PSI0
             ELSE
                PSI0(:,:,NBH+1:2*NBH)=THIS%PSI0
             END IF
          END DO
          CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH*2,NB*2,PSI0,PSI0,AUXMAT2)
          DEALLOCATE(PSI0)

          QMAT=(0.D0,0.D0)
          DO I=1,NB  ! UP UP
             DO J=1,NB
                QMAT(2*I-1,2*J-1)=(AUXMAT(I,J)+AUXMAT2(I,J))*0.5D0
             ENDDO
          ENDDO
          DO I=NB+1,2*NB  ! DOWN DOWN
             DO J=NB+1,2*NB
                QMAT(2*I,2*J)=(AUXMAT(I,J)+AUXMAT2(I,J))*0.5D0
             ENDDO
          ENDDO
          DO I=1,NB  ! UP DOWN
             DO J=NB+1,2*NB
                QMAT(2*I-1,2*J)=(AUXMAT(I,J)+AUXMAT2(I,J))*0.5D0
             ENDDO
          ENDDO         
          DO I=NB+1,2*NB  ! DOWN UP
             DO J=1,NB
                QMAT(2*I,2*J-1)=(AUXMAT(I,J)+AUXMAT2(I,J))*0.5D0
             ENDDO
          ENDDO
          DEALLOCATE(AUXMAT)
          DEALLOCATE(AUXMAT2)
       END IF
       RETURN
     END SUBROUTINE WAVES_SPINOROVERLAP
!
!.....................................................................
module wavesfixrho_module
! this module is used to keep the density fixed
  real(8), allocatable    :: qlm(:,:)
  real(8), allocatable    :: rho(:,:)
  real(8), allocatable    :: denmat(:,:,:,:)
end module wavesfixrho_module
!
!      .................................................................
       subroutine waves_fixrhoget(nrl,ndimd,lmrxx,nat,qlm_,rho_,denmat_)
       use wavesfixrho_module 
       implicit none
       integer(4) ,intent(in)  :: nrl,ndimd,lmrxx,nat
       real(8)    ,intent(out) :: qlm_(lmrxx,nat),rho_(nrl,ndimd)
       real(8)    ,intent(out) :: denmat_(lmrxx,lmrxx,ndimd,nat)
!      ******************************************************************
       if (.not.allocated(qlm)) then
          allocate(qlm(lmrxx,nat))
          allocate(rho(nrl,ndimd))
          allocate(denmat(lmrxx,lmrxx,ndimd,nat))
          call waves$fixrhoread()
       end if
       qlm_(:,:)=qlm(:,:)
       rho_(:,:)=rho(:,:)
       denmat_(:,:,:,:)=denmat(:,:,:,:)
       end
!
!      .................................................................
       subroutine waves_fixrhoset(nrl,ndimd,lmrxx,nat,qlm_,rho_,denmat_)
       use wavesfixrho_module 
       implicit none
       integer(4) ,intent(in)  :: nrl,ndimd,lmrxx,nat
       real(8)    ,intent(in)  :: qlm_(lmrxx,nat),rho_(nrl,ndimd)
       real(8)    ,intent(in)  :: denmat_(lmrxx,lmrxx,ndimd,nat)
!      ******************************************************************
       if (.not.allocated(qlm)) then
          allocate(qlm(lmrxx,nat))
          allocate(rho(nrl,ndimd))
          allocate(denmat(lmrxx,lmrxx,ndimd,nat))
       end if
       qlm(:,:)=qlm_(:,:)
       rho(:,:)=rho_(:,:)
       denmat(:,:,:,:)=denmat_(:,:,:,:)
       end
!
!      .................................................................
       subroutine waves$fixrhoread()
       use wavesfixrho_module
       implicit none
       integer(4)                 :: nfil
!      ******************************************************************
       CALL FILEHANDLER$SETFILE('FIXRHO',.FALSE.,'fixrho.bin')
       CALL FILEHANDLER$SETSPECIFICATION('FIXRHO','FORM','UNFORMATTED')
       CALL FILEHANDLER$SETSPECIFICATION('FIXRHO','POSITION','REWIND')
       CALL FILEHANDLER$UNIT('FIXRHO',NFIL)
       read(nfil)qlm(:,:)
       read(nfil)rho(:,:)
       read(nfil)denmat(:,:,:,:)
       call filehandler$close('FIXRHO')
       return
       end
!
!      .................................................................
       subroutine waves$fixrhowrite()
       use wavesfixrho_module
       implicit none
       integer(4)                 :: nfil
!      ******************************************************************
       CALL FILEHANDLER$SETFILE('FIXRHO',.FALSE.,'fixrho.bin')
       CALL FILEHANDLER$SETSPECIFICATION('FIXRHO','FORM','UNFORMATTED')
       CALL FILEHANDLER$SETSPECIFICATION('FIXRHO','POSITION','REWIND')
       CALL FILEHANDLER$UNIT('FIXRHO',NFIL)
       write(nfil)qlm(:,:)
       write(nfil)rho(:,:)
       write(nfil)denmat(:,:,:,:)
       call filehandler$close('FIXRHO')
       return
       end
       
       

!
!.....................................................................
MODULE TOTALSPIN_MODULE
REAL(8)               :: TOTSPIN(4) ! DIFFERS FROM JOHANNES' VERSION
END MODULE TOTALSPIN_MODULE
!
!      .................................................................
       SUBROUTINE WAVES_TOTALSPIN(NB,NKPT,IKPT,NSPIN,OCC,QMAT)
!      *****************************************************************
!      **  CALCULATES <S^2> IN UNITS OF HBAR^2                        **
!      **                                                             **
!      **                                                             **
!      *****************************************************************
       USE MPE_MODULE
       USE TOTALSPIN_MODULE, ONLY: TOTSPIN ! DIFFERS FROM JOHANNES' VERSION
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NB,NKPT,IKPT,NSPIN
       REAL(8)   ,INTENT(IN) :: OCC(NB,NKPT,NSPIN) 
       COMPLEX(8),INTENT(IN) :: QMAT(2,NB*NSPIN,2,NB*NSPIN) ! Q_I,J,K,L=1/2*<PSI_I,K|PSI_J,K>  
       COMPLEX(8)            :: SUM,PART,SPINX,SPINY,SPINZ
       REAL(8)               :: IOCC(NB*NSPIN),SVAR
       INTEGER(4)            :: I,J,NTASKS,THISTASK,NBD
       
       CALL MPE$QUERY(NTASKS,THISTASK)

       DO I=1,NB
          IOCC(I)=OCC(I,IKPT,1)
       END DO
       IF (NSPIN.EQ.2) THEN
          DO I=NB+1,2*NB
             IOCC(I)=OCC(I-NB,IKPT,2)
          END DO
       END IF
       NBD=NSPIN*NB !ONE BAND FOR ONE ELECTRON

       SUM=0.D0
       DO I=1,NBD
             PART=QMAT(1,I,2,I)+QMAT(2,I,1,I)
             SUM=SUM+PART*IOCC(I)
       END DO
       SPINX=SUM

       SUM=0.D0
       DO I=1,NBD
             PART=QMAT(1,I,2,I)-QMAT(2,I,1,I)
             SUM=SUM+PART*IOCC(I)
       END DO
       SPINY=SUM/(0,1)

       SUM=0.D0
       DO I=1,NBD
             PART=QMAT(1,I,1,I)-QMAT(2,I,2,I)
             SUM=SUM+PART*IOCC(I)
       END DO
       SPINZ=SUM


       SUM=0.D0
       DO I=1,NBD
          SUM=SUM+IOCC(I)*0.75D0
       END DO
       
       SUM=SUM+SPINX**2+SPINY**2+SPINZ**2
       PRINT*,'TOTALSPIN: PART 3/4+X+Y+Z: ',SUM
       
       SVAR=SUM
       DO I=1,NBD
          DO J=1,NBD
             PART=-(QMAT(1,I,1,J)-QMAT(2,I,2,J))*(QMAT(1,J,1,I)-QMAT(2,J,2,I))
             SUM=SUM+PART*SQRT(IOCC(I)*IOCC(J))
          END DO
       END DO
       PRINT*,'TOTALSPIN: PART UP UP    : ',SUM-SVAR

       SVAR=SUM
       DO I=1,NBD
          DO J=1,NBD
             PART=-4.D0*QMAT(1,I,2,J)*QMAT(2,J,1,I)
             SUM=SUM+PART*SQRT(IOCC(I)*IOCC(J))
          END DO
       END DO       
       PRINT*,'TOTALSPIN: PART DOWN UP  : ',SUM-SVAR

       IF (THISTASK.EQ.1) THEN
          PRINT*,'SPIN: 2*<S_X>: ',2.D0*REAL(SPINX),' HBAR'
          PRINT*,'SPIN: 2*<S_Y>: ',2.D0*REAL(SPINY),' HBAR'
          PRINT*,'SPIN: 2*<S_Z>: ',2.D0*REAL(SPINZ),' HBAR'
          PRINT*,'TOTALSPIN: <S^2>: ',REAL(SUM),' HBAR^2'
          PRINT*,'SPIN QUATUM NUMBER S=',-0.5+SQRT(0.25+REAL(SUM))
          IF (IKPT.EQ.1) TOTSPIN=0.D0 ! SUM UP TOTAL SPIN OVER K-POINTS
          TOTSPIN(1)=TOTSPIN(1)+DBLE(SUM)
          TOTSPIN(2)=TOTSPIN(2)+DBLE(SPINX)
          TOTSPIN(3)=TOTSPIN(3)+DBLE(SPINY)
          TOTSPIN(4)=TOTSPIN(4)+DBLE(SPINZ)
       END IF
       RETURN
     END SUBROUTINE WAVES_TOTALSPIN
!
!     ..................................................................
      SUBROUTINE WAVES$REPORTSPIN(NFIL)
!     ******************************************************************
!     **  WAVES$GET                                                   **
!     **  GET                                                         **
!     ******************************************************************
      USE TOTALSPIN_MODULE, ONLY: TOTSPIN ! DIFFERS FROM JOHANNES' VERSION
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)               :: SVAR
!     ******************************************************************
      IF (NSPIN.EQ.1.AND.NDIM.EQ.1) THEN
         RETURN
      ELSE
         CALL REPORT$TITLE(NFIL,'TOTAL SPIN ANALYSIS')
         IF (NDIM.EQ.2) CALL REPORT$R8VAL(NFIL,'<S_X>',TOTSPIN(2),'HBAR')
         IF (NDIM.EQ.2) CALL REPORT$R8VAL(NFIL,'<S_Y>',TOTSPIN(3),'HBAR')
         CALL REPORT$R8VAL(NFIL,'<S_Z>',TOTSPIN(4),'HBAR')
         SVAR=DSQRT(TOTSPIN(2)**2+TOTSPIN(3)**2+TOTSPIN(4)**2)
         IF (NDIM.EQ.2) CALL REPORT$R8VAL(NFIL,'|S|',SVAR,'HBAR')
!        CALL REPORT$R8VAL(NFIL,'TOTAL SPIN <S^2>',TOTSPIN(1),'HBAR^2')
!        SVAR=-0.5+SQRT(0.25+TOTSPIN(1))
!        CALL REPORT$R8VAL(NFIL,'SPIN QUANTUM NUMBER S',SVAR,'')
         RETURN
      END IF
    END SUBROUTINE WAVES$REPORTSPIN
!
!     ....................................................................
      subroutine waves_comparepsi(id,ngl,ndim,nbh,psi1,psi2)
      USE WAVES_MODULE, only : gset,delt
      implicit none
      character(*)    ,intent(in) :: id
      INTEGER(4)      ,INTENT(IN) :: NGL
      INTEGER(4)      ,INTENT(IN) :: NDIM
      INTEGER(4)      ,INTENT(IN) :: NBH
      COMPLEX(8)      ,INTENT(IN) :: PSI1(NGL,NDIM,NBH)
      COMPLEX(8)      ,INTENT(IN) :: PSI2(NGL,NDIM,NBH)
      complex(8)                  :: csvar
      integer(4)                  :: ig,idim,ib
      real(8)                     :: sum,sumtot
      real(8)                     :: rbas(3,3),gbas(3,3),cellvol
!     *********************************************************************
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL GBASS(RBAS,GBAS,CELLVOL)
      sumtot=0.d0
      do ib=1,nbh
        sum=0.d0
        do ig=1,ngl
          do idim=1,ndim
            csvar=psi1(ig,idim,ib)-psi2(ig,idim,ib)
            sum=sum+GSET%MPSI(IG)*(real(csvar)**2+aimag(csvar)**2)
          enddo
        enddo
        sumtot=sumtot+sum
        write(*,*)' kineticenergytest ',ib,sum*CELLVOL/DELT**2,id
      enddo
        write(*,*)' kineticenergytest total',sumtot*CELLVOL/DELT**2,id
      return
      end
