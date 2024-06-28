!
!*******************************************************************************
!*******************************************************************************
!**                                                                           **
!** CONTINUATION OF WAVES OBJECT: SEE PAW_WAVES1.F90 FOR HEADER               **
!**                                                                           **
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$ORTHOGONALIZE()
!     **************************************************************************
!     ** ENFORCES THE ORTHONORMALITY CONDITION OF THE WAVE FUNCTIONS          **
!     **************************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : THIS &
     &                        ,GSET &
     &                        ,MAP &
     &                        ,TSAFEORTHO &
     &                        ,TSWAPSTATES &
     &                        ,NKPTL &
     &                        ,NSPIN &
     &                        ,NDIM &
     &                        ,WAVEEKIN2 &
     &                        ,DELT,ANNEE &
     &                        ,WAVES_SELECTWV !SUBROUTINE
      IMPLICIT NONE
      COMPLEX(8),ALLOCATABLE :: OPROJ(:,:,:)
      REAL(8)   ,ALLOCATABLE :: MARR(:)
      REAL(8)   ,ALLOCATABLE :: R0(:,:)      !CURRENT ATOMIC POSITIONS
      REAL(8)   ,ALLOCATABLE :: RP(:,:)      ! NEXT ATOMIC POSITIONS
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
      INTEGER(4)             :: NBX
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: IKPT,ISPIN,IPRO,IAT,ISP
      INTEGER(4)             :: IB,IDIM,IG,I,J
      INTEGER(4)             :: NGL,NBH,NB,LMNX,LNX,LMX,LN
      INTEGER(4)             :: NAT,IND
      REAL(8)   ,PARAMETER   :: DSMALL=1.D-12
      REAL(8)                :: SVAR
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
!     **************************************************************************
                             CALL TRACE$PUSH('WAVES$ORTHOGONALIZE')
                             CALL TIMING$CLOCKON('WAVES$ORTHOGONALIZE')
      NPRO=MAP%NPRO
      NAT=MAP%NAT
      CALL CELL$GETL4('MOVE',TSTRESS)
!
!     ==========================================================================
!     == COLLECT OCCUPATIONS                                                  ==
!     ==========================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     ==  CALCULATE FORCE OF CONSTRAINT                                       ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NBH=THIS%NBH
          NB=THIS%NB
IF(1.EQ.0) THEN ! CHANGE FOR KAESTNERS CONJUGATE GRADIENT
!
!         ======================================================================
!         ==  EVALUATE  DO<P|PSI>  ASSUMING <PRO|PSI> IS STILL VALID          ==
!         ======================================================================
          ALLOCATE(OPROJ(NDIM,NBH,NPRO))
          OPROJ(:,:,:)=(0.D0,0.D0)
          IPRO=1
          DO IAT=1,NAT
            ISP=MAP%ISP(IAT)
            CALL SETUP$ISELECT(ISP)
            LNX=MAP%LNX(ISP)
            LMNX=MAP%LMNX(ISP)
            ALLOCATE(DO(LNX,LNX))
            CALL SETUP$GETR8A('DO',LNX*LNX,DO) 
                                        !FORMER CALL SETUP$1COVERLAP(ISP,LNX,DO)
            CALL WAVES_OPROJ(LNX,MAP%LOX(1:LNX,ISP),DO,NDIM,LMNX,NBH &
      &          ,THIS%PROJ(:,:,IPRO:IPRO+LMNX-1),OPROJ(:,:,IPRO:IPRO+LMNX-1))
            DEALLOCATE(DO)
            IPRO=IPRO+LMNX
            CALL SETUP$UNSELECT()
          ENDDO
!
!         ======================================================================
!         ==  ADD  |PSI>+|P>DO<P|PSI>                                         ==
!         ======================================================================
          ALLOCATE(THIS%OPSI(NGL,NDIM,NBH))
!         == THIS LOOP IS DONE EXPLICIT, BECAUSE THE IFC10 COULD NOT HANDLE IT
          DO IB=1,NBH
            DO IDIM=1,NDIM
              DO IG=1,NGL
                THIS%OPSI(IG,IDIM,IB)=THIS%PSI0(IG,IDIM,IB)
              ENDDO
            ENDDO
          ENDDO
          CALL WAVES_ADDPRO(MAP,GSET,NAT,R0,NGL,NDIM,NBH,NPRO,THIS%OPSI,OPROJ)
          DEALLOCATE(OPROJ)
ELSE
          ALLOCATE(THIS%OPSI(NGL,NDIM,NBH))
          THIS%OPSI(:,:,:)=THIS%PSI0(:,:,:)
!++++++++++++++++++++++++ FROM HERE +++++++++++++++++++++++++++++++++++++
!         __ THIS$PROJ=<PTILDE|THIS%PSI0>_______________________________________
          CALL WAVES_OPSI(NB,NBH,NPRO,NAT,NGL,R0,THIS%PROJ,THIS%OPSI)
!++++++++++++++++++++++++ TO HERE +++++++++++++++++++++++++++++++++++++++
END IF
!
!         ======================================================================
!         ==  DIVIDE BY WAVE FUNCTION MASS                                    ==
!         ======================================================================
          ALLOCATE(MARR(NGL))
          CALL PLANEWAVE$GETR8A('G2',NGL,MARR)
!PB070802          IF(ASSOCIATED(GSET%DMPSI)) THEN
!PB070802            DO IG=1,NGL
!PB070802!             SVAR=1.D0+ANNEE+GSET%DMPSI(IG)
!PB070802              SVAR=1.D0+ANNEE 
!PB070802              MARR(IG)=DELT**2/(SVAR*GSET%MPSI(IG))
!PB070802            ENDDO
!PB070802          ELSE
            SVAR=DELT**2/(1.D0+ANNEE)
            DO IG=1,NGL
!MARR(IG)=EMASS*(1.D0+EMASSCG2*MARR(IG))  !OLD VERSION
!MARR(IG)=DELT**2/(MARR(IG)*(1.D0+ANNEE))
              MARR(IG)=SVAR/GSET%MPSI(IG)
            ENDDO
!PB070802          END IF
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
      CALL ATOMLIST$GETR8A('R(+)',0,3*NAT,RP)
      IF(TSTRESS) THEN
!       == PREDICT NEW POSITIONS =======================================
        CALL CELL$GETR8A('TP',9,RBAS)
        CALL CELL$GETR8A('MAPTOCELL',9,MAPTOCELL)
        DO IAT=1,NAT
           RP(:,IAT)=MATMUL(MAPTOCELL,RP(:,IAT))
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
        DO IKPT=1,NKPTL
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
!
!         ==  UPDATE PROJECTOR FUNCTIONS =======================================
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
            CALL SETUP$UNSELECT()
          ENDDO
          DEALLOCATE(G2)
!
!         ==  UPDATE SPHERICAL HARMONICS  ======================================
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
!         == NOW THE STRAINED SPHERICAL HARMONICS ==============================
          CALL WAVES_STRAINEDYLM(NGL,LMX,GVEC,GSET%YLM,GSET%SYLM)
          DEALLOCATE(GVEC)
!
!         == RESCALE WAVE FUNCTIONS ============================================
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
!     ==========================================================================
!     ==  NOW ORTHOGONALIZE                                                   ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NBH=THIS%NBH
          NB=THIS%NB
!
!         ======================================================================
!         ==  CALCULATE PROJECTIONS FOR THE NEW POSITIONS                     ==
!         ======================================================================
          CALL WAVES_PROJECTIONS(MAP,GSET,NAT,RP,NGL,NDIM,NBH,NPRO &
     &                                                     ,THIS%PSIM,THIS%PROJ)
          CALL MPE$COMBINE('K','+',THIS%PROJ)
          ALLOCATE(OPROJ(NDIM,NBH,NPRO))
          CALL WAVES_PROJECTIONS(MAP,GSET,NAT,RP,NGL,NDIM,NBH,NPRO &
     &                                                         ,THIS%OPSI,OPROJ)
          CALL MPE$COMBINE('K','+',OPROJ)
!
!         ======================================================================
!         ==  1C-OVERLAP OF <PSI0|PSI0>, <OPSI|PSI0> AND <OPSI|OPSI>          ==
!         ======================================================================
          ALLOCATE(MAT(NB,NB))
          CALL WAVES_1COVERLAP(MAP,NDIM,NBH,NB,NPRO,THIS%PROJ,THIS%PROJ,MAT)
          ALLOCATE(OMAT(NB,NB))
          CALL WAVES_1COVERLAP(MAP,NDIM,NBH,NB,NPRO,OPROJ,THIS%PROJ,OMAT)
          ALLOCATE(OOMAT(NB,NB))
          CALL WAVES_1COVERLAP(MAP,NDIM,NBH,NB,NPRO,OPROJ,OPROJ,OOMAT)
!
!         ======================================================================
!         ==  NOW ADD OVERLAP OF PSEUDO WAVE FUNCTIONS                        ==
!         ======================================================================
          ALLOCATE(AUXMAT(NB,NB))
          CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB,THIS%PSIM,THIS%PSIM,AUXMAT)
          DO I=1,NB
            DO J=1,NB
              MAT(I,J)=MAT(I,J)+AUXMAT(I,J)
            ENDDO
          ENDDO
          CALL WAVES_OVERLAP(.FALSE.,NGL,NDIM,NBH,NB,THIS%OPSI,THIS%PSIM,AUXMAT)
          DO I=1,NB
            DO J=1,NB
              OMAT(I,J)=OMAT(I,J)+AUXMAT(I,J)
            ENDDO
          ENDDO
          CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB,THIS%OPSI,THIS%OPSI,AUXMAT)
          DO I=1,NB
            DO J=1,NB
              OOMAT(I,J)=OOMAT(I,J)+AUXMAT(I,J)
            ENDDO
          ENDDO
          DEALLOCATE(AUXMAT)
!
!         ======================================================================
!         ==  CALCULATE LAGRANGE PARAMETERS                                   ==
!         ======================================================================
          ALLOCATE(LAMBDA(NB,NB))
          LAMBDA(:,:)=THIS%RLAM0(:,:)
          DO I=1,NB
            DO J=I+1,NB
              CSVAR=0.5D0*(LAMBDA(I,J)+CONJG(LAMBDA(J,I)))
              LAMBDA(I,J)=CSVAR
              LAMBDA(J,I)=CONJG(CSVAR)
            ENDDO
          ENDDO
!
!===============================================================          
          IF(.NOT.TSAFEORTHO) THEN
            ALLOCATE(SMAP(NB))
            DO I=1,NB
              SMAP(I)=I
            ENDDO
            IF(TSWAPSTATES) THEN
!             ==================================================================
!             == CONSTRUCT SMAP SO THAT |LAMBDA(SMAP(I),SMAP(I))| INCREASES   ==
!             ==================================================================
              SVAR=0.D0   ! WILL BECOME MAX OFFDIAGONAL ELEMENTS OF LAMBDA 
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
!
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
              RMAT=REAL(MAT,KIND=8)
              ALLOCATE(ROMAT(NB,NB))
              ROMAT=REAL(OMAT,KIND=8)
              ALLOCATE(ROOMAT(NB,NB))
              ROOMAT=REAL(OOMAT,KIND=8)
              ALLOCATE(RLAMBDA(NB,NB))
              RLAMBDA=REAL(LAMBDA,KIND=8)
              CALL WAVES_ORTHO_X(NB,OCC(1,IKPT,ISPIN) &
       &                        ,ROOMAT,RMAT,ROMAT,RLAMBDA)
              LAMBDA=CMPLX(RLAMBDA,0.D0,KIND=8)
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
!         ======================================================================
!         ==  CALCULATE |PSI(+)>=|PSI>+|CHI>LAMBDA                            ==
!         ======================================================================
          CALL WAVES_ADDOPSI(NGL,NDIM,NBH,NB,THIS%PSIM,THIS%OPSI,LAMBDA)
          DEALLOCATE(THIS%OPSI)
!PRINT*,'WARNING FROM WAVES$ORTHOGONALIZE:'
!PRINT*,'MAKE SURE THAT PDOS AND GRAPHICS PICK UP A CONSISTENT SET OF '
!PRINT*,'WAVE FUNCTIONS  AND PROJECTOR FUNCTIONS'
          CALL WAVES_ADDOPROJ(NPRO,NDIM,NBH,NB,THIS%PROJ,OPROJ,LAMBDA)
          DEALLOCATE(OPROJ)
!
!         ======================================================================
!         ======================================================================
!         ==  MAP LAMBDA ONTO RLAM0                                           ==
!         ======================================================================
!         == MASS IS NOT INCLUDED HERE (MASS TENSOR IS TAKEN CARE OF IN       ==
!         == PSIBAR AND (1/M)O|PSI>                                           ==
!         ======================================================================
!         == RLAM0 IS DEFINED VIA THE EQUATION OF MOTION                      ==
!         ==         M|PSIDOTDOT>=-H|PSI>+O|PSI>*RLAM0                        ==
!         == RLAM0(I,J)*OCC(J) IS THE LAGRANGE PARAMETER, WHICH IS            ==
!         == HERMITEAN. ITS MATRIX ELEMENTS INVOLVING UNOCCUPIED STATES       ==
!         == VANISHES.                                                        ==
!         == RLAM0 IS NOT HERMITEAN! FOR A SYSTEM AT REST THE RELATION        ==
!         ==        RLAM0=<PSI|H|PSI>                                         ==
!         == HOLDS. THEN THE MATRIX FOR STATES WITH DIFFERENT OCCUPATIONS     ==
!         == VANISH.                                                          ==
!         ======================================================================
          IF(.NOT.ASSOCIATED(THIS%RLAM0))ALLOCATE(THIS%RLAM0(NB,NB))
          THIS%RLAM0(:,:)=LAMBDA(:,:)
          DEALLOCATE(LAMBDA)
!
!         ======================================================================
!         ==  TEST ORTHONORMALITY                                             ==
!         ======================================================================
          IF(TTEST) THEN
            ALLOCATE(OPROJ(NDIM,NBH,NPRO))
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,RP,NGL,NDIM,NBH,NPRO,THIS%PSIM,OPROJ)
            CALL MPE$COMBINE('K','+',OPROJ)
            ALLOCATE(AUXMAT(NB,NB))
            CALL WAVES_1COVERLAP(MAP,NDIM,NBH,NB,NPRO,OPROJ,OPROJ,AUXMAT)
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
!     ==  NOW TRANSFORM BACK INTO ORIGINAL CELL                       ==
!     ==================================================================
      IF(TSTRESS) THEN
        CALL CELL$GETR8A('T0',9,RBAS)
        DO IKPT=1,NKPTL
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE WAVES_ADDOPSI(NGL,NDIM,NBH,NB,PSIBAR,OPSI,LAMBDA)
!      *************************************************************************
!      **  |PSI(+)>=|PSIBAR>+O|PSI(0)>*LAMBDA                                 **
!      **                                                                     **
!      **  CAUTION! OVERWRITES OPSI! USE WAVES_ADDPSI WITHOUT "O"             **
!      **   FOR A ROUTINE THAT PRESERVES OPSI                                 **
!      **                                                                     **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)   :: NGL
       INTEGER(4),INTENT(IN)   :: NDIM
       INTEGER(4),INTENT(IN)   :: NBH
       INTEGER(4),INTENT(IN)   :: NB
       COMPLEX(8),INTENT(INOUT):: PSIBAR(NGL*NDIM,NBH)
       COMPLEX(8),INTENT(INOUT):: OPSI(NGL*NDIM,NBH)
       COMPLEX(8),INTENT(IN)   :: LAMBDA(NB,NB)
       LOGICAL(4)              :: TINV
       INTEGER(4)              :: IBH1,IBH2,I,IDIM
       INTEGER(4)              :: IB1A,IB1B,IB2A,IB2B
       COMPLEX(8),ALLOCATABLE  :: LAMBDA1(:,:)
       COMPLEX(8),ALLOCATABLE  :: LAMBDA2(:,:)
       COMPLEX(8),ALLOCATABLE  :: TPSI(:)
       COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
       COMPLEX(8)              :: CSVAR1,CSVAR2
       INTEGER(4)              :: NGLNDIM
!      *************************************************************************
                               CALL TIMING$CLOCKON('WAVES_ADDOPSI')
       TINV=NBH.NE.NB
       NGLNDIM=NGL*NDIM
       IF(.NOT.TINV) THEN
!        == PSIBAR=PSIBAR+OPSI*LAMBDA ==========================================
         CALL LIB$ADDPRODUCTC8(.FALSE.,NGLNDIM,NBH,NBH,OPSI,LAMBDA,PSIBAR)
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
!        == ADD O|PSI_+>LAMBDA1 ================================================
         CALL LIB$ADDPRODUCTC8(.FALSE.,NGLNDIM,NBH,NBH,OPSI,LAMBDA1,PSIBAR)
         DEALLOCATE(LAMBDA1)
!        == INVERT OPSI ========================================================
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
!        == ADD O|PSI_+>LAMBDA1 ================================================
         CALL LIB$ADDPRODUCTC8(.FALSE.,NGLNDIM,NBH,NBH,OPSI,LAMBDA2,PSIBAR)
         DEALLOCATE(LAMBDA2)
       END IF
                               CALL TIMING$CLOCKOFF('WAVES_ADDOPSI')
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE WAVES_ADDPSI(NGL,NDIM,NBH,NB,PSIBAR,OPSI,LAMBDA)
!      *************************************************************************
!      **                                                                     **
!      **  |PSI(+)>=|PSIBAR>+|OPSI>*LAMBDA                                    **
!      **                                                                     **
!      **  CAUTION: PRESEVERS OPSI, BUT PRODUCES AN ADDITIONAL ARRAY TPSI     **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)   :: NGL
       INTEGER(4),INTENT(IN)   :: NDIM
       INTEGER(4),INTENT(IN)   :: NBH
       INTEGER(4),INTENT(IN)   :: NB
       COMPLEX(8),INTENT(INOUT):: PSIBAR(NGL*NDIM,NBH)
       COMPLEX(8),INTENT(IN)   :: OPSI(NGL*NDIM,NBH)
       COMPLEX(8),INTENT(IN)   :: LAMBDA(NB,NB)
       LOGICAL(4)              :: TINV
       INTEGER(4)              :: IBH1,IBH2,I,IDIM
       INTEGER(4)              :: IB1A,IB1B,IB2A,IB2B
       COMPLEX(8),ALLOCATABLE  :: LAMBDA1(:,:)
       COMPLEX(8),ALLOCATABLE  :: LAMBDA2(:,:)
       COMPLEX(8),ALLOCATABLE  :: TPSI(:,:)
       COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
       COMPLEX(8)              :: CSVAR1,CSVAR2
       INTEGER(4)              :: NGLNDIM
!      *************************************************************************
                               CALL TIMING$CLOCKON('WAVES_ADDPSI')
       TINV=NBH.NE.NB
       NGLNDIM=NGL*NDIM
       IF(.NOT.TINV) THEN
!        == PSIBAR=PSIBAR+OPSI*LAMBDA ==========================================
         CALL LIB$ADDPRODUCTC8(.FALSE.,NGLNDIM,NBH,NBH,OPSI,LAMBDA,PSIBAR)
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
!        == ADD O|PSI_+>LAMBDA1 ================================================
         CALL LIB$ADDPRODUCTC8(.FALSE.,NGLNDIM,NBH,NBH,OPSI,LAMBDA1,PSIBAR)
         DEALLOCATE(LAMBDA1)
!        == INVERT OPSI ========================================================
         ALLOCATE(TPSI(NGLNDIM,NBH))
         DO IBH1=1,NBH
           I=1
           DO IDIM=1,NDIM
             CALL PLANEWAVE$INVERTG(NGL,OPSI(I,IBH1),TPSI(I,IBH1))
             I=I+NGL
           ENDDO
         ENDDO
!
!        == ADD O|PSI_+>LAMBDA1 ================================================
         CALL LIB$ADDPRODUCTC8(.FALSE.,NGLNDIM,NBH,NBH,TPSI,LAMBDA2,PSIBAR)
         DEALLOCATE(TPSI)
         DEALLOCATE(LAMBDA2)
       END IF
                               CALL TIMING$CLOCKOFF('WAVES_ADDPSI')
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE WAVES_ADDOPROJ(NPRO,NDIM,NBH,NB,PROJ,OPROJ,LAMBDA)
!      *************************************************************************
!      **                                                                     **
!      **  |PSI(+)>=|PSIBAR>+O|PSI(0)>*LAMBDA                                 **
!      **                                                                     **
!      **                                                                     **
!      **                                                                     **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)   :: NDIM
       INTEGER(4),INTENT(IN)   :: NBH
       INTEGER(4),INTENT(IN)   :: NB
       INTEGER(4),INTENT(IN)   :: NPRO
       COMPLEX(8),INTENT(INOUT):: PROJ(NDIM,NBH,NPRO)
       COMPLEX(8),INTENT(IN)   :: OPROJ(NDIM,NBH,NPRO)
       COMPLEX(8),INTENT(IN)   :: LAMBDA(NB,NB)
       LOGICAL(4)              :: TINV
       INTEGER(4)              :: IBH1,IBH2,IDIM
       INTEGER(4)              :: IB1A,IB1B,IB2A,IB2B
       COMPLEX(8),ALLOCATABLE  :: LAMBDA1(:,:)
       COMPLEX(8),ALLOCATABLE  :: LAMBDA2(:,:)
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
       INTEGER(4)            :: IBH1,IBH2
       INTEGER(4)            :: IB1A,IB1B,IB2A,IB2B,I,J
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
         CALL MPE$COMBINE('K','+',MAT)
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
           RE=0.5D0* REAL(TMAT(IBH1,IBH2),KIND=8)
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
               MAT(I,J)=CMPLX(MAT2(I-IB1A+1,J-IB2A+1),0.D0,KIND=8)
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
           RE=0.5D0* REAL(TMAT(IBH1,IBH2),KIND=8)
           IM=0.5D0*AIMAG(TMAT(IBH1,IBH2))
           MAT2(1,1)=+RE
           MAT2(1,2)=+IM
           MAT2(2,1)=+IM
           MAT2(2,2)=-RE
           DO I=IB1A,IB1B
             DO J=IB2A,IB2B
               MAT(I,J)=MAT(I,J)+CMPLX(MAT2(I-IB1A+1,J-IB2A+1),0.D0,KIND=8)
             ENDDO
           ENDDO
         ENDDO
       ENDDO
       DEALLOCATE(TMAT)
       CALL MPE$COMBINE('K','+',MAT)
                             CALL TIMING$CLOCKOFF('WAVES_OVERLAP')
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_1COVERLAP(MAP,NDIM,NBH,NB,NPRO,PROJ1,PROJ2,MAT)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE 1C CONTRIBUTION TO THE OVERLAP MATRIX        **
!     **  ( RESULT IS ADDED TO "MAT"! )                               **
!     **                                                              **
!     *******************************************P.E. BLOECHL, (1992)***
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY: MAP_TYPE
      IMPLICIT NONE
      TYPE(MAP_TYPE),INTENT(IN) :: MAP
      INTEGER(4)  ,INTENT(IN)   :: NDIM
      INTEGER(4)  ,INTENT(IN)   :: NBH
      INTEGER(4)  ,INTENT(IN)   :: NB
      INTEGER(4)  ,INTENT(IN)   :: NPRO
      COMPLEX(8)  ,INTENT(IN)   :: PROJ1(NDIM,NBH,NPRO) ! <P|PSI1>
      COMPLEX(8)  ,INTENT(IN)   :: PROJ2(NDIM,NBH,NPRO) ! <P|PSI2>
      COMPLEX(8)  ,INTENT(OUT)  :: MAT(NB,NB)
      REAL(8)     ,ALLOCATABLE  :: DOVER(:,:)
      COMPLEX(8)  ,ALLOCATABLE  :: PROJ1A(:,:,:)  !(NDIM,NBH,NPROX)
      COMPLEX(8)  ,ALLOCATABLE  :: PROJ2A(:,:,:)  !(NDIM,NBH,NPROX)
      COMPLEX(8)                :: CMAT1(NBH,NBH),CMAT2(NBH,NBH)
      INTEGER(4)                :: NTASKS,THISTASK
      INTEGER(4)                :: NPROX
      INTEGER(4)                :: LMNX
      INTEGER(4)                :: IAT,ISP
      INTEGER(4)                :: LN1,L1
      INTEGER(4)                :: LN2,L2
      INTEGER(4)                :: IPRO,IPROX,IPRO1,IPRO2
      INTEGER(4)                :: IBH1,IBH2,IB1,IB2,IDIM,LNX
      COMPLEX(8)                :: CSVAR1,CSVAR2
      LOGICAL(4)                :: TINV
!     **************************************************************************
      CALL MPE$QUERY('K',NTASKS,THISTASK)
      MAT(:,:)=(0.D0,0.D0)
      TINV=NB.NE.NBH
!
!     ==========================================================================     
!     == DETERMINE, HOW MANY PROJECTIONS ARE TREATED ON THE CURRENT PROCESSOR ==
!     ==========================================================================
      IF(NTASKS.EQ.1) THEN
        NPROX=NPRO
      ELSE 
        NPROX=0
        DO IAT=1,MAP%NAT
          IF(MOD(IAT-1,NTASKS).NE.THISTASK-1) CYCLE
          ISP=MAP%ISP(IAT)
!         __ SELECTION FOR PARALLEL PROCESSING____________________________
          NPROX=NPROX+MAP%LMNX(ISP)
        ENDDO
      END IF
!
!     ==========================================================================     
!     == LIMIT PROJ1 AND PROJ2 TO THE TERMS REQUIRED ON THE CURRENT PROCESSOR ==
!     == PROJ2A=DO*<P|PSI>                                                    ==
!     ==========================================================================
      ALLOCATE(PROJ1A(NDIM,NBH,NPROX))
      ALLOCATE(PROJ2A(NDIM,NBH,NPROX))
      PROJ1A(:,:,:)=(0.D0,0.D0)
      PROJ2A(:,:,:)=(0.D0,0.D0)
      IPRO=0
      IPROX=0
      DO IAT=1,MAP%NAT
        ISP=MAP%ISP(IAT)
!       __ SELECTION FOR PARALLEL PROCESSING____________________________________
        IF(MOD(IAT-1,NTASKS).NE.THISTASK-1) THEN
          IPRO=IPRO+MAP%LMNX(ISP)
          CYCLE
        END IF
        LMNX=MAP%LMNX(ISP)
!
!       == FOLD DOWN INDEX ARRAY ===============================================
!
        PROJ1A(:,:,IPROX+1:IPROX+LMNX)=PROJ1(:,:,IPRO+1:IPRO+LMNX)
!
!       == FOLD DOWN INDEX ARRAY FOR PROJ2 AND MULTIPLY WITH DOVER =============
        LNX=MAP%LNX(ISP)
        ALLOCATE(DOVER(LNX,LNX))
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8A('DO',LNX*LNX,DOVER) 
        CALL SETUP$ISELECT(0)
                                     !FORMER CALL SETUP$1COVERLAP(ISP,LNX,DOVER)
!
        IPRO1=IPROX
        DO LN1=1,LNX
          L1=MAP%LOX(LN1,ISP)
          IPRO2=IPRO
          DO LN2=1,LNX
            L2=MAP%LOX(LN2,ISP)
            IF(L1.EQ.L2) THEN
             PROJ2A(:,:,IPRO1+1:IPRO1+2*L1+1)=PROJ2A(:,:,IPRO1+1:IPRO1+2*L1+1) &
      &                        +DOVER(LN1,LN2)*PROJ2(:,:,IPRO2+1:IPRO2+2*L1+1)
            END IF
            IPRO2=IPRO2+2*L2+1
          ENDDO
          IPRO1=IPRO1+2*L1+1
        ENDDO
        DEALLOCATE(DOVER)
        IPRO=IPRO+LMNX
        IPROX=IPROX+LMNX
      ENDDO
!
!     ==========================================================================     
!     ==  MAT(I,J)= SUM_{K,L} <PSI_I|P_K>   * [ DO(K,L)*<P_L|PSI_J>]          ==
!     ==          = SUM_{K}   PROJ1A^*(K,I) *   PROJ2A(K,J)                   ==
!     ==========================================================================
      IF(.NOT.TINV) THEN
        MAT(:,:)=(0.D0,0.D0)
        DO IBH1=1,NBH
          DO IBH2=1,NBH
            DO IPRO1=1,NPROX
              DO IDIM=1,NDIM
                MAT(IBH1,IBH2)=MAT(IBH1,IBH2) &
     &                   +CONJG(PROJ1A(IDIM,IBH1,IPRO1))*PROJ2A(IDIM,IBH2,IPRO1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CALL MPE$COMBINE('K','+',MAT)
      ELSE
        CMAT1(:,:)=(0.D0,0.D0)
        CMAT2(:,:)=(0.D0,0.D0)
        DO IBH1=1,NBH
          DO IBH2=1,NBH
            DO IPRO1=1,NPROX
              DO IDIM=1,NDIM
                CMAT1(IBH1,IBH2)=CMAT1(IBH1,IBH2) &
    &                   +CONJG(PROJ1A(IDIM,IBH1,IPRO1))*PROJ2A(IDIM,IBH2,IPRO1)
                IF(TINV) THEN
                  CMAT2(IBH1,IBH2)=CMAT2(IBH1,IBH2) &
    &                          +PROJ1A(IDIM,IBH1,IPRO1)*PROJ2A(IDIM,IBH2,IPRO1)
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CALL MPE$COMBINE('K','+',CMAT1)
        CALL MPE$COMBINE('K','+',CMAT2)
      END IF
      DEALLOCATE(PROJ1A)
      DEALLOCATE(PROJ2A)
!
!     ==========================================================================     
!     ==  UNRAVEL SUPER WAVE FUNCTIONS  PSISUP=PSI1+CI*PSI2 WITH REAL PSI1/2  ==
!     == CMAT1=CONJG(PSISUP)*PSISUP=PSI1*PSI1+PSI2*PSI2+I[PSI1*PSI2-PSI2*PSI1]==
!     == CMAT2=      PSISUP *PSISUP=PSI1*PSI1-PSI2*PSI2+I[PSI1*PSI2+PSI2*PSI1]==
!     == => 0.5*[CMAT1+CMAT2]=PSI1*PSI1+I*PSI1*PSI2                           ==
!     ==    0.5*[CMAT1-CMAT2]=PSI2*PSI2-I*PSI2*PSI1                           ==
!     ==========================================================================
      IF(TINV) THEN
        MAT(:,:)=(0.D0,0.D0)
        DO IBH2=THISTASK,NBH,NTASKS   ! PARALLELIZED LOOP_______________________
          IB2=2*(IBH2-1)
          DO IBH1=1,NBH
            IB1=2*(IBH1-1)
            CSVAR1=0.5D0*(CMAT1(IBH1,IBH2)+CMAT2(IBH1,IBH2))
            CSVAR2=0.5D0*(CMAT1(IBH1,IBH2)-CMAT2(IBH1,IBH2))
            MAT(IB1+1,IB2+1)=REAL(CSVAR1,KIND=8)
            MAT(IB1+1,IB2+2)=AIMAG(CSVAR1)
            MAT(IB1+2,IB2+1)=-AIMAG(CSVAR2)
            MAT(IB1+2,IB2+2)=REAL(CSVAR2,KIND=8)
          ENDDO
        ENDDO
        CALL MPE$COMBINE('K','+',MAT)
      ENDIF  
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
       INTEGER(4)            :: I,J,K,L,N
       COMPLEX(8)            :: CSVAR
       REAL(8)               :: SVAR,SVAR1,SVAR2
       REAL(8)               :: MAXDEV
       LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
       REAL(8)   ,PARAMETER  :: TOL=1.D-10
       INTEGER(4)            :: NU,NU0,IU,JU
       REAL(8)               :: R8SMALL
!      *****************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_Y_C')
       R8SMALL=TINY(R8SMALL)*1.D+5
       R8SMALL=1.D-10
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
!CHANGED THIS TO MAKE IT SIMPLER PB 070127
         SVAR = REAL( REAL(B(N,N))**2-A(N,N)*C(N,N) )
         SVAR1=-REAL(B(N,N),KIND=8)
         IF(SVAR.GE.0.D0) THEN
           SVAR2=SQRT(SVAR)
!          == CHOOSE THE SIGN OF THE ROOT SUCH RESULT IS SMALL
           IF(SVAR1*SVAR2.GE.0.D0) THEN
             Z(N)=SVAR1-SVAR2
           ELSE
             Z(N)=SVAR1+SVAR2
           END IF
         ELSE
           PRINT*,'ORTHOGONALIZATION FAILED! TRYING BEST APPROXIMATION...',SVAR
           Z(N)=SVAR1
         END IF
         Z(N)=Z(N)/REAL(C(N,N)+R8SMALL,KIND=8)
!        == TEST FOR NANS (NAN = NOT A NUMBER)
         IF(SVAR.NE.SVAR) THEN
PRINT*,'MARKE 1',N,B(N,N),A(N,N),C(N,N)
PRINT*,'MARKE 2',N,B(N,N)**2-A(N,N)*C(N,N)
PRINT*,'PHIPHI',(PHIPHI(I,I),I=1,NB)
PRINT*,'A     ',(A(I,I),I=1,NB)
           CALL ERROR$MSG('NAN DETECTED')
           CALL ERROR$STOP('WAVES_ORTHO_Y_C')
         END IF
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
!        == ORTHOGONALIZE HIGHER PHIS TO THIS PHI                     ==
!        == PHI(J)=PHI(J)+CHI(N)*Z(J)       J>N                       ==
!        ===============================================================
         DO IU=1,NU
           I=MAP(IU)
           Z(I)=(0.D0,0.D0)
         ENDDO
         DO IU=NU+1,NB
           I=MAP(IU)
           Z(I)=-CONJG(A(I,N)/(B(N,N)+R8SMALL))
         ENDDO
!
!        ===============================================================
!        == NOW UPDATE MATRICES                                       ==
!        ===============================================================
         NU0=NU+1   !SET N0=N FOR FAST CALCULATION AND N0=1 FOR TESTS
!N0=1
!          == A(N,M)+B(N,N)*DELTA(M)=0  ======================
         DO JU=NU0,NB
           J=MAP(JU)
           DO I=1,NB
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
!
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
!        == ORTHOGONALIZE HIGHER CHIS TO THIS PHI                     ==
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
           Z(I)=-CONJG(B(I,N)/(B(N,N)+R8SMALL))
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
!IF(ABS(B(N,N)).GT.HUGE(B)) THEN
! PRINT*,B(N,N),Z
! CALL ERROR$STOP('ERROR')
!END IF
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
       REAL(8)               :: R8SMALL
!      *****************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_Y')
       R8SMALL=100.D0*TINY(R8SMALL)
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
         SVAR=1.D0-PHIPHI(N,N)*CHICHI(N,N)/(CHIPHI(N,N)+R8SMALL)**2
         IF(SVAR.GT.0.D0) THEN
           SVAR=CHIPHI(N,N)/(CHICHI(N,N)+R8SMALL)*(SQRT(SVAR)-1.D0)             
         ELSE
           PRINT*,'ORTHOGONALIZATION FAILED! TRYING BEST APPROXIMATION...'
           SVAR=-CHIPHI(N,N)/(CHICHI(N,N)+R8SMALL)
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
!        == ORTHOGONALIZE HIGHER PHIS TO THIS PHI                     ==
!        == PHI(M)=PHI(M)+CHI(N)*DELTA(M)   M>N                       ==
!        ===============================================================
                CALL TRACE$PASS('ORTHOGONALIZE HIGHER PHIS TO THIS PHI')
         DELTA=0.D0
         DO MU=NU+1,NB
           M=MAP(MU)
!          == PHIPHI(N,M)+CHIPHI(N,N)*DELTA(M)=0  ======================
           DELTA(M)=-PHIPHI(N,M)/(CHIPHI(N,N)+R8SMALL)
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
!        == ORTHOGONALIZE HIGHER CHIS TO THIS PHI                     ==
!        == CHI(M)=CHI(M)+CHI(N)*DELTA(M)   M>N                       ==
!        ===============================================================
                CALL TRACE$PASS('ORTHOGONALIZE HIGHER CHIS TO THIS PHI')
         DELTA=0.D0
         DO MU=NU+1,NB
           M=MAP(MU)
!          == |CHI(M)>=|CHI(M)>+|CHI(N)>*DELTA(M) ============================
!          == CHIPHI(M,N)+CHIPHI(N,N)*DELTA(M)=0
           DELTA(M)=-CHIPHI(M,N)/(CHIPHI(N,N)+R8SMALL)
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
           IF(ABS(TEST(I,J)-PHIPHI(I,J)).GT.1.D-7) THEN
!             PRINT*,'ERROR IN PHIPHI FOR ',I,J,PHIPHI(I,J),TEST(I,J)
           END IF
           IF(ABS(TEST1(I,J)-CHIPHI(J,I)).GT.1.D-7) THEN
!             PRINT*,'ERROR IN CHIPHI FOR ',I,J,CHIPHI(J,I),TEST1(I,J)
           END IF
           IF(ABS(TEST2(I,J)-CHICHI(I,J)).GT.1.D-7) THEN
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
!     **  NEW VERSION WITH DIAGONALIZATION FOR CHIPSI                 **
!     **                                                              **
!     **  THE METHOD IS DESCRIBED IN :                                **
!     **    R.CAR AND M.PARRINELLO,                                   **
!     **    IN "SIMPLE MOLECULAR SYSTEMSAT VERY HIGH DENSITY", P.455  **
!     **    ED. A.POLIAN, PLOUBEYRE AND N.BOCCARA                     **
!     **    (PLENUM PUBLISHING CORPORATION,1989)                      **
!     **                                                              **
!     **  THE MATRIX LAMBDA IS NOT THE LAGRANGE MULTIPLIER BUT THE    **
!     **  NON-HERMITAEAN FORM. THE LAGRANGE MULTIPLIERS, WHICH ARE    **
!     **  HERMITEAN, ARE OBTAINED FROM THIS LAMBDA MATRIX BY          **
!     **  MULTIPLICATION WITH THE OCCUPATION TO THE RIGHT.            **
!     **                                                              **
!     *******************************************P.E. BLOECHL, (1992)***
!     IMPLICIT NONE
      REAL(8)   ,PARAMETER     :: EPS    = 1.D-8
      REAL(8)   ,PARAMETER     :: DSMALL = 1.D-12
      INTEGER(4),PARAMETER     :: ITERX    = 200
      LOGICAL(4),PARAMETER     :: TPR    = .FALSE.
      INTEGER(4),INTENT(IN)    :: NB
      REAL(8)   ,INTENT(IN)    :: OCC(NB)
      COMPLEX(8),INTENT(INOUT) :: LAMBDA(NB,NB)
      COMPLEX(8),INTENT(IN)    :: PSIPSI(NB,NB)
      COMPLEX(8),INTENT(IN)    :: CHIPSI(NB,NB)   !
      COMPLEX(8),INTENT(IN)    :: CHICHI(NB,NB)
      COMPLEX(8),ALLOCATABLE   :: GAMN(:,:)  
      REAL(8)                  :: EIG(NB)
      INTEGER(4)               :: ITER,I,J,K ! RUNNING VARIABLES
      REAL(8)                  :: DIGAM,SVAR,EIGI ! AUXILARY VARIABLES
      COMPLEX(8)               :: HAUX(NB,NB)    
      COMPLEX(8)               :: U(NB,NB)       
      LOGICAL(4)               :: TCONVERGED
      LOGICAL(4)               :: TESSL=.TRUE.
      COMPLEX(8)               :: CSVAR
      REAL(8)                  :: OCCI,OCCJ
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
!         __HAUX(I,J) = HAUX(I,J)+2*CHIPSI(I,J)___________________________
          HAUX=HAUX+2.D0*CHIPSI
!         __GAMN(I,J) = LAMBDA(K,I)*HAUX(K,J)_____________________________
          CALL LIB$SCALARPRODUCTC8(.FALSE.,NB,NB,LAMBDA,NB,HAUX,GAMN)
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
!       ========================================================================
!       == CHECK CONVERGENCE MAXVAL(ABS(OVERLAP-1))<EPS ; GAMN=OVERLAP-1      ==
!       ========================================================================
!       __IN CASE THE MAXVAL EVALUATION OVER THE ENTIRE ARRAY IS EXPENSIVE______
!       __CONSIDER AN EXPLICIT LOOP THAT SETS TCONVERGED=FALSE AND EXITS, WHEN__
!       __THE FIRST ELEMENT WITH ABS(GAMN)>EPS IS ENCOUNTERES___________________
        DIGAM=MAXVAL(ABS(GAMN))
        TCONVERGED=DIGAM.LT.EPS
!       __CHECK CONVERGENCE_____________________________________________________
        IF(TCONVERGED) GOTO 9000
!
!       __CHECK WHETHER LOOP DIVERGES___________________________________________
        IF(TPR.OR.DIGAM.GT.1.D+5)  THEN
          WRITE(*,*)'WAVES_ORTHO_X_C: ITER ',ITER,' MAX(OVERLAP-1)=',DIGAM
          IF(DIGAM.GT.1.D+200) THEN
            CALL ERROR$MSG('ITERATION LOOP FOR ORTHOGONALIZATION')
            CALL ERROR$MSG('IS ABOUT TO DIVERGE')
            CALL ERROR$I4VAL('ITER',ITER)
            CALL ERROR$R8VAL('(MAX(OVERLAP-1)',DIGAM)
            CALL ERROR$STOP('WAVES_ORTHO_X_C')
          END IF
        END IF
!
!       ==================================================================
!       ==  OBTAIN CHANGE OF THE LAMBDA MATRIX                          ==
!       ==================================================================
!       == TRANSFORM OVERLAP MATRIX GAMN
        IF(TESSL) THEN
!         ----  HAUX(I,L)=U(K,I)*H0(K,L)
          CALL LIB$SCALARPRODUCTC8(.FALSE.,NB,NB,U,NB,GAMN,HAUX)
!         ----  GAMN(I,J)=HAUX(I,L)*U(L,I)
          CALL LIB$MATMULC8(NB,NB,NB,HAUX,U,GAMN)
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
!         ----  GAMN(I,J)=HAUX(I,L)*U(L,I)
          CALL LIB$DYADSUMC8(NB,NB,NB,HAUX,U,GAMN)
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
      ENDDO
      CALL ERROR$MSG('LOOP FOR ORTHOGONALIZATION IS NOT CONVERGED')
      CALL ERROR$MSG('THIS IS NOT AN UNUSUAL PROBLEM DURING STARTUP')
      CALL ERROR$MSG('1) SPECIFY SAFEORTHO=F, WHICH IS MORE ROBUST')
      CALL ERROR$MSG('2) REDUCE PLANE WAVE CUTOFF')
      CALL ERROR$MSG('3) REDUCE TIME STEP WHILE KEEPING WAVE FUNCTION MASS')
      CALL ERROR$MSG('4) CHECK FOR GHOST STATES')
      CALL ERROR$MSG('5) POTENTIAL PROBLEM WITH LINEARLY DEPENDENT INITIAL')
      CALL ERROR$MSG('   WAVE FUNCTIONS. CHECK AFTER GRAM-SCHMIDT ORTHO.')
      CALL ERROR$MSG('AFTER THE FIRST 5-10 ITERATIONS YOU CAN GO BACK')
      CALL ERROR$MSG('TO THE ORIGINAL SETTINGS')
      CALL ERROR$STOP('WAVES_ORTHO_X_C')
!
9000  CONTINUE
      DEALLOCATE(GAMN)
                             CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_ORTHO_X(NB,OCC,CHICHI,PSIPSI,CHIPSI,LAMBDA)
!     **************************************************************************
!     **                                                                      **
!     **  IMPOSES THE ORTHOGONALITY CONSTRAINT ONTO THE ELECTRONS             **
!     **                                                                      **
!     **  NEW VERSION WITH DIAGONALIZATION FOR CHIPSI                         **
!     **                                                                      **
!     **                                                                      **
!     **  THE METHOD IS DESCRIBED IN :                                        **
!     **    R.CAR AND M.PARRINELLO, IN                                        **
!     **    "SIMPLE MOLECULAR SYSTEMS AT VERY HIGH DENSITY", P. 455           **
!     **    ED. A.POLIAN, PLOUBEYRE AND N.BOCCARA                             **
!     **    (PLENUM PUBLISHING CORPORATION,1989)                              **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1992)***********
      IMPLICIT NONE
      REAL(8)   ,PARAMETER     :: EPS    = 1.D-8
      REAL(8)   ,PARAMETER     :: DSMALL = 1.D-12
      INTEGER(4),PARAMETER     :: MAX    = 200
      LOGICAL(4),PARAMETER     :: TPR    = .FALSE.
      INTEGER(4),INTENT(IN)    :: NB
      REAL(8)   ,INTENT(IN)    :: OCC(NB)
      REAL(8)   ,INTENT(INOUT) :: LAMBDA(NB,NB)
      REAL(8)   ,INTENT(IN)    :: PSIPSI(NB,NB)
      REAL(8)   ,INTENT(IN)    :: CHIPSI(NB,NB)   !
      REAL(8)   ,INTENT(IN)    :: CHICHI(NB,NB)
      REAL(8)   ,ALLOCATABLE   :: GAMN(:,:)  
      REAL(8)                  :: EIG(NB)
      INTEGER(4)               :: ITER,I,J ! RUNNING VARIABLES
      REAL(8)                  :: DIGAM,SVAR,EIGI ! AUXILARY VARIABLES
      REAL(8)                  :: HAUX(NB,NB)    
      REAL(8)                  :: U(NB,NB)       
      REAL(8)                  :: OCCI,OCCJ
!     **************************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_X')
      ALLOCATE(GAMN(NB,NB))
!
!     ==========================================================================
!     ==  CALCULATE  PSIPSI(I,J)= <PSIBAR(I)|PSIBAR(J)>-1                     ==
!     ==        AND  CHIPSI(I,J)   = <PSI0(I)|PSIBAR(J)>                      ==
!     ==========================================================================
!
!     ==========================================================================
!     ==  DIAGONALIZE 0.5*(CHIPSI(I,J)+CHIPSI(J,I))                           ==
!     ==========================================================================
      CALL LIB$DIAGR8(NB,CHIPSI,EIG,U)
!CALL DIAG(NB,NB,CHIPSI,EIG,U)
!WRITE(*,FMT='("EIG",20E10.3)')EIG
!DO I=1,NB
!  WRITE(*,FMT='("U",I2,20E10.3)')I,U(I,:)
!ENDDO
!
!     ==========================================================================
!     ==========================================================================
!     ==  ITERATIVE CALCULATION OF GAMMA                                      ==
!     ==========================================================================
!     ==========================================================================
      DO ITER=1,MAX
!PRINT*,'==================',ITER,'==========================='
!       ========================================================================
!       ==  CALCULATE <PHI(+)|PHI(+)>-1 WITH PRESENT LAMBDA                   ==
!       ==  GAMN(I,J)=PSIPSI(I,J)+LAMBDA(K,I)*CHIPSI(K,J)                     ==
!       ==                       +CHIPSI(K,I)*LAMBDA(K,J)                     ==
!       ==           +LAMBDA(K,I)*CHICHI(K,L)*LAMBDA(L,J)-1(I,J)              ==
!       ========================================================================
!       __GAMN(I,J) = CHICHI(I,K)*LAMBDA(K,J)___________________________________
        CALL LIB$MATMULR8(NB,NB,NB,CHICHI,LAMBDA,HAUX)
!CALL DGEMUL(CHICHI,NB,'N',LAMBDA,NB,'N',HAUX,NB,NB,NB,NB)
!       __HAUX(I,J) = HAUX(I,J)+2*CHIPSI(I,J)___________________________________
        HAUX=HAUX+2.D0*CHIPSI
!CALL DAXPY(NB*NB,2.D0,CHIPSI,1,HAUX,1)
!       __GAMN(I,J) = LAMBDA(K,I)*HAUX(K,J)_____________________________________
        CALL LIB$SCALARPRODUCTR8(.FALSE.,NB,NB,LAMBDA,NB,HAUX,GAMN)
!CALL DGEMUL(LAMBDA,NB,'T',HAUX,NB,'N',GAMN,NB,NB,NB,NB)
!       __GAMN(I,J) = GAMN(I,J)-1_______________________________________________
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
!
!       ========================================================================
!       == CHECK CONVERGENCE MAXVAL(ABS(OVERLAP-1))<EPS ; GAMN=OVERLAP-1      ==
!       ========================================================================
        DIGAM=MAXVAL(ABS(GAMN))
!       __CHECK CONVERGENCE_____________________________________________________
        IF(DIGAM.LT.EPS) GOTO 9000
!
!       __CHECK WHETHER LOOP DIVERGES___________________________________________
        IF(TPR.OR.DIGAM.GT.1.D+5)  THEN
          WRITE(*,*)'WAVES_ORTHO_X: ITER ',ITER,' MAX(OVERLAP-1)=',DIGAM
          IF(DIGAM.GT.1.D+200) THEN
            CALL ERROR$MSG('ITERATION LOOP FOR ORTHOGONALIZATION')
            CALL ERROR$MSG('IS ABOUT TO DIVERGE')
            CALL ERROR$I4VAL('ITER',ITER)
            CALL ERROR$R8VAL('(MAX(OVERLAP-1)',DIGAM)
            CALL ERROR$STOP('WAVES_ORTHO_X')
          END IF
        END IF
!
!       ========================================================================
!       ==  OBTAIN CHANGE OF THE LAMBDA MATRIX                                ==
!       ========================================================================
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
!       ========================================================================
!       ==  PROPAGATE GAMMA                                                   ==
!       ========================================================================
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
!       ========================================================================
!       == SYMMETRIZE LAMBDA                                                  ==
!       ========================================================================
        DO I=1,NB
          DO J=1,NB
            IF(OCC(I).LT.OCC(J)) THEN
              LAMBDA(I,J)=LAMBDA(J,I)*(OCC(I)+DSMALL)/(OCC(J)+DSMALL)
            END IF
          ENDDO
        ENDDO
      ENDDO
      CALL ERROR$MSG('LOOP FOR ORTHOGONALIZATION IS NOT CONVERGED')
      CALL ERROR$MSG('THIS IS NOT AN UNUSUAL PROBLEM DURING STARTUP')
      CALL ERROR$MSG('1) SPECIFY SAFEORTHO=F, WHICH IS MORE ROBUST')
      CALL ERROR$MSG('2) REDUCE PLANE WAVE CUTOFF')
      CALL ERROR$MSG('3) REDUCE TIME STEP WHILE KEEPING WAVE FUNCTION MASS')
      CALL ERROR$MSG('4) CHECK FOR GHOST STATES')
      CALL ERROR$MSG('AFTER THE FIRST 5-10 ITERATIONS YOU CAN GO BACK')
      CALL ERROR$MSG('TO THE ORIGINAL SETTINGS')
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
!     *******************************************P.E. BLOECHL, (1999)***
      USE WAVES_MODULE, ONLY: MAP_TYPE,GSET_TYPE
      USE MPE_MODULE
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
      LOGICAL(4)      ,PARAMETER  :: TTEST=.FALSE.
      COMPLEX(8)      ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)                  :: NPRO
      INTEGER(4)                  :: I,J
      INTEGER(4)                  :: IBH1,IBH2,IB1A,IB1B,IB2A,IB2B
      COMPLEX(8)      ,ALLOCATABLE:: X(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: PROJ(:,:,:)
      COMPLEX(8)      ,ALLOCATABLE:: OVERLAP(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: AUXMAT(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: X1(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: X2(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: PSIINV(:,:,:)
      COMPLEX(8)                  :: CSVAR1,CSVAR2
      LOGICAL(4)                  :: TINV
      REAL(8)                     :: NORM(NB)
      COMPLEX(8)                  :: XTWOBYTWO(2,2)
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
      CALL MPE$COMBINE('K','+',PROJ)
!     
!     ==================================================================
!     ==  OVERLAP OF <PSI0|PSI0>,                                     ==
!     ==================================================================
      ALLOCATE(OVERLAP(NB,NB))
      ALLOCATE(AUXMAT(NB,NB))
      CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB,PSI,PSI,OVERLAP)
      CALL WAVES_1COVERLAP(MAP,NDIM,NBH,NB,NPRO,PROJ,PROJ,AUXMAT)
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
!== SPEICHERZUGRIFFSFEHLER IFC10
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
!       == PSI_I=PSI_I+PSI_J*X1_JI+PSIINV_J*X2_JI
        CALL PLANEWAVE$ADDPRODUCT(' ',NGL,NDIM,NBH,PSI,NBH,PSIINV,X1)
        CALL PLANEWAVE$ADDPRODUCT('-',NGL,NDIM,NBH,PSI,NBH,PSIINV,X2)
        DEALLOCATE(PSIINV)
        DEALLOCATE(X1)
        DEALLOCATE(X2)
      ELSE
        ALLOCATE(PSIINV(NGL,NDIM,NB))
        PSIINV=PSI
!       == PSI_I=PSI_I+PSI_J*X1_JI
        CALL PLANEWAVE$ADDPRODUCT(' ',NGL,NDIM,NB,PSI,NB,PSIINV,X)
        DEALLOCATE(PSIINV)
        DEALLOCATE(X)
      ENDIF
! 
!     =================================================================
!     ==                                                             ==
!     =================================================================
      IF(TTEST) THEN
        ALLOCATE(PROJ(NDIM,NBH,NPRO))
        CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NBH,NPRO,PSI,PROJ)
        CALL MPE$COMBINE('K','+',PROJ)
        ALLOCATE(AUXMAT(NB,NB))
        CALL WAVES_1COVERLAP(MAP,NDIM,NBH,NB,NPRO,PROJ,PROJ,AUXMAT)
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
      USE WAVES_MODULE, ONLY : THIS &
     &                        ,GSET &
     &                        ,NKPTL &
     &                        ,NSPIN &
     &                        ,THAMILTON &
     &                        ,WAVEEKIN1 &
     &                        ,WAVEEKIN2 &
     &                        ,WAVES_SELECTWV & !SUBROUTINE
     &                        ,OPTIMIZERTYPE
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
      IF(OPTIMIZERTYPE.EQ.'CG') RETURN   !KAESTNERCG
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
WRITE(*,FMT='("MAPTOCELL=",9F8.4)')MAPTOCELL
PRINT*,'CELLSCALE ',CELLSCALE
      ELSE 
        CELLSCALE=1.D0
      ENDIF
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          TPSI=>THIS%PSI0
          THIS%PSI0=>THIS%PSIM
          THIS%PSIM=>TPSI
!===============================================================================
!I DO NOT UNDERSTAND WHY THE FOLLOWING HAD BEEN COMMENTED OUT.  I AM
!CONCERNED, BECAUSE THE CELL DYNAMICS WITH MOVING ATOMS IS NOT ENERGY
!CONSERVING, WHILE THAT WITHOUT MOVING ATOMS AND THAT WITHOUT MOVING
!CELL WORK FINE.
!===============================================================================
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
        DO IKPT=1,NKPTL
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_RANDOMIZE(NG,NDIM,NB,AMPLITUDE,G2,PSI)
!     **************************************************************************
!     **                                                                      **
!     **  CREATE RANDOM WAVE FUNCTIONS                                        **
!     **                                                                      **
!     **  THE MAXIMUM WEIGHT OF THE WAVE FUNCTIONS AT EPW[RY]=GC2             **
!     **                                                                      **
!     **  IF THE PARAMETER TDETERMINISTIC IS SET TO TRUE, THE RESULT          **
!     **  IS INDEPENDENT OF THE PARALLELIZATION                               **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1991,2007)******
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB              ! #(BANDS)
      INTEGER(4),INTENT(IN)    :: NG              ! #(PLANE WAVES),MAX
      INTEGER(4),INTENT(IN)    :: NDIM            ! #(PLANE WAVES),MAX
      REAL(8)   ,INTENT(IN)    :: AMPLITUDE       ! SCALE FACTOR
      REAL(8)   ,INTENT(IN)    :: G2(NG)          ! G**2
      COMPLEX(8),INTENT(INOUT) :: PSI(NG,NDIM,NB) ! PS-WAVE FUNCTION
      INTEGER(4)               :: IB,IG,IDIM
      REAL(8)   ,PARAMETER     :: PI=4.D0*ATAN(1.D0)
      REAL(8)                  :: FAC
      REAL(8)   ,PARAMETER     :: GC2=10.D0
      REAL(8)                  :: SCALE(NG)
      REAL(8)                  :: REC,RIM
!     **************************************************************************
!
!     ==========================================================================
!     == DETERMINE ENVELOPE FUNCTION OF THE RANDOM NUMBERS TO AVOID LARGE     ==
!     == CONTRIBUTIONS WITH LARGE KINETIC ENERGY                              ==
!     ==========================================================================
      FAC=2.D0*SQRT(PI*GC2)
      FAC=FAC**3/REAL(NDIM,KIND=8)*(2.D0/3.D0)
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
            PSI(IG,IDIM,IB)=PSI(IG,IDIM,IB)+CMPLX(REC,RIM,KIND=8)*SCALE(IG)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_INITIALIZERANDOM(NG,NDIM,NB,G,PSI)
!     **************************************************************************
!     **  CREATE RANDOM WAVE FUNCTIONS                                        **
!     **                                                                      **
!     **  THE MAXIMUM WEIGHT OF THE WAVE FUNCTIONS AT EPW[RY]=GC2             **
!     **                                                                      **
!     **  ATTENTION: THE CHOICE TDETERMINISTIC MAY LEAD TO LINEARLY DEPENDENT **
!     **     WAVE FUNCTIONS THAT WILL FAIL DURING THE ORTHOGONALIZATION       **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1991)***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB              ! #(BANDS)
      INTEGER(4),INTENT(IN)    :: NG              ! #(PLANE WAVES),MAX
      INTEGER(4),INTENT(IN)    :: NDIM            ! #(PLANE WAVES),MAX
      REAL(8)   ,INTENT(IN)    :: G(3,NG)         ! G
      COMPLEX(8),INTENT(OUT)   :: PSI(NG,NDIM,NB) ! PS-WAVE FUNCTION
      LOGICAL(4),SAVE          :: TINI=.FALSE.
      LOGICAL(4),PARAMETER     :: TDETERMINISTIC=.TRUE.
      INTEGER(4),PARAMETER     :: NRAN=10000
      REAL(8)   ,ALLOCATABLE,SAVE     :: XRAN(:)   !(NRAN)
      REAL(8)   ,PARAMETER     :: GC2=3.D0
      REAL(8)   ,PARAMETER     :: PI=4.D0*ATAN(1.D0)
      REAL(8)                  :: G2
      INTEGER(4)               :: IB,IG,IDIM,I
      REAL(8)                  :: FAC
      REAL(8)                  :: SCALE(NG)
      REAL(8)                  :: RANDARR(NG)
      REAL(8)                  :: REC,RIM
      INTEGER(4)               :: ISVAR,IND
      REAL(8)                  :: SVAR,RAN,GX,GY,GZ
      REAL(8)                  :: RSMALL
!     **************************************************************************
      RSMALL=100.D0*TINY(RSMALL)
!     ==========================================================================
!     == CREATE A SERIES OF NRAN PSEUDO-RANDOM NUMBERS                        ==
!     ==========================================================================
      IF(.NOT.TINI) THEN
        IF(TDETERMINISTIC) THEN
          ALLOCATE(XRAN(NRAN))
          DO I=1,NRAN
!           == CONSTRUCT A FIXED SERIES OF RANDOM NUMBERS USING
!           == THE WELL-DEFINED MINIMAL STANDARD LINEAR CONGRUENTIAL
!           == RANDOM NUMBER GENERATOR (IN PAW_GENERALPURPOSE.F90)
            CALL RANDOM_MSLNG(XRAN(I))
          ENDDO
        ELSE
!         __WITH TDETERMINISTIC=.FALSE., XRAN IS USED BEFORE ALLOCATION AND 
!         __BEFORE ASSIGNING A CALUE
          CALL ERROR$MSG('HARDWIRED TDETERMINISTIC=.FALSE. IS NOT ALLOWED')
          CALL ERROR$L4VAL('TDETERMINISTIC',TDETERMINISTIC)
          CALL ERROR$STOP('WAVES_INITIALIZERANDOM')
        END IF
        TINI=.TRUE.
      END IF
!
!     ==========================================================================
!     == DETERMINE ENVELOPE FUNCTION OF THE RANDOM NUMBERS TO AVOID LARGE     ==
!     == CONTRIBUTIONS WITH LARGE KINETIC ENERGY                              ==
!     ==========================================================================
      FAC=2.D0*SQRT(PI*GC2)
      FAC=FAC**3/REAL(NDIM,KIND=8)*(2.D0/3.D0)
      FAC=1.D0/FAC
      DO IG=1,NG
!       == ENVELOPE FUNCTION THAT ATTENUATES HIGH-G COMPONENTS
        G2=G(1,IG)**2+G(2,IG)**2+G(3,IG)**2
        SCALE(IG)=FAC*EXP(-0.5D0*G2/GC2)
!       == CREATE A RANDOM-LIKE NUMBER DEPENDING ON THE ANGLE
        SVAR=1.D0/SQRT(G2+RSMALL)
!       == MAP G-COMPONENTS ONTO INTERVAL [-1,1]
        GX=SVAR*G(1,IG)
        GY=SVAR*G(2,IG)
        GZ=SVAR*G(3,IG)
!       == MAP ONTO NUMBERS IN THE INTERVAL [0.1]
        GX=0.5D0*(1.D0+GX)
        GY=0.5D0*(1.D0+GY)
        GZ=0.5D0*(1.D0+GZ)
        RAN=XRAN(1+INT(GX*REAL(NRAN-1)))
        RAN=XRAN(1+INT(RAN*GY*REAL(NRAN-1)))
        RAN=XRAN(1+INT(RAN*GZ*REAL(NRAN-1)))
        RANDARR(IG)=RAN
      ENDDO
!
!     ==========================================================================
!     == SELECT RANDOM NUMBERS                                                ==
!     ==========================================================================
      IF(TDETERMINISTIC) THEN
        IND=0
        DO IB=1,NB
          DO IDIM=1,NDIM
            IND=IND+2
            IF(IND+1.GT.NRAN) IND=1
            DO IG=1,NG
              SVAR=0.5D0*(XRAN(IND)*RANDARR(IG))
              ISVAR=1+INT(SVAR*REAL(NRAN,KIND=8))
              REC=XRAN(ISVAR)
!             == THE FACTOR PI IS ONLY TO MAKE RIM DIFFERENT FROM REC ==========
              SVAR=0.5D0*(XRAN(IND+1)*RANDARR(IG))
              ISVAR=1+INT(SVAR*REAL(NRAN,KIND=8))
              RIM=XRAN(ISVAR)
              REC=2.D0*REC-1.D0
              RIM=2.D0*RIM-1.D0
              PSI(IG,IDIM,IB)=CMPLX(REC,RIM,KIND=8)*SCALE(IG)
            ENDDO
          ENDDO
        ENDDO
      ELSE 
        CALL LIB$RANDOMSEED
        DO IB=1,NB
          DO IDIM=1,NDIM
            DO IG=1,NG
              CALL LIB$RANDOM(REC)
              CALL LIB$RANDOM(RIM)
              REC=2.D0*REC-1.D0
              RIM=2.D0*RIM-1.D0
              PSI(IG,IDIM,IB)=CMPLX(REC,RIM,KIND=8)*SCALE(IG)
            ENDDO
          ENDDO
        ENDDO
      END IF
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
      INTEGER(4)            :: IKPT
      REAL(8)               :: RY
      INTEGER(4)            :: NG
      INTEGER(4)            :: NTASKS,THISTASK
      INTEGER(4)            :: IKPTL
      CHARACTER(64)         :: STRING
!     ******************************************************************
      CALL CONSTANTS('RY',RY)
      CALL REPORT$TITLE(NFIL,'WAVE FUNCTIONS')
      CALL WAVES_SELECTWV(1,1)
      CALL REPORT$I4VAL(NFIL,'NUMBER OF BANDS',THIS%NB,' ')

      CALL REPORT$I4VAL(NFIL,'NUMBER OF K-POINTS',NKPT,' ')
      CALL REPORT$I4VAL(NFIL,'NUMBER OF SPINS',NSPIN,' ')
      CALL REPORT$I4VAL(NFIL,'NUMBER OF SPINOR COMPONENTS',NDIM,' ')
      CALL REPORT$R8VAL(NFIL,'WAVEFUNCTION PLANE WAVE CUTOFF',EPWPSI/RY,'RY')
      CALL REPORT$R8VAL(NFIL,'WAVEFUNCTION MASS',EMASS,'A.U.')
      CALL REPORT$R8VAL(NFIL,'G**2 ENHANCEMENT OF FUNCTION MASS',EMASSCG2,' ')
      CALL REPORT$R8VAL(NFIL,'BUCKET POTENTIAL STARTS AT 0.5G^2=' &
     &                       ,EPWPSI0/RY,'RY')
      CALL REPORT$R8VAL(NFIL,'BUCKET POTENTIAL PREFACTOR',D2EPWPSI,'H')
      CALL REPORT$R8VAL(NFIL,'CUTOFF SCALE-FACTOR FOR NEIGHBORLIST' & 
     &                      ,SCALERCUT,' ')
      CALL REPORT$CHVAL(NFIL,'|R1-R2|<(RCOV1+RCOV2)*SCALERCUT',' ')
      IF(ANNEE.NE.0.D0) THEN
        CALL REPORT$R8VAL(NFIL,'FRICTION',ANNEE,' ')
      END IF
      IF(TSTOP) THEN
        CALL REPORT$CHVAL(NFIL,'INITIAL VELOCITY IS SET TO','ZERO')
      END IF
!     IF(TRANDOMIZE) THEN
!       CALL REPORT$R8VAL(NFIL &
!    &      ,'INITIAL VELOCITIES ARE RANDOMIZED WITH ENERGY',AMPRE,'H')
!     END IF
      CALL REPORT$L4VAL(NFIL,'SAFEORTHO',TSAFEORTHO)
      IF(.NOT.TSAFEORTHO) THEN
        WRITE(NFIL,FMT='("I.E. ORTHOGONALIZATION CONVERGES TO EIGENSTATES."' &
     &                 //'," (NO STRICT ENERGY CONSERVATION)")')
      ELSE
        WRITE(NFIL,FMT='("I.E. ORTHOGONALIZATION CONSERVES ENERGY"' &
     &                 //'," (DO NOT USE WITH VARIABLE OCCUPATIONS)")')
      END IF
!     
!     ================================================================
!     ==  REPORT INFORMATION ABOUT G-VECTORS                        ==
!     ================================================================
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IKPTL=0
      DO IKPT=1,NKPT
        IF(KMAP(IKPT).EQ.THISTASK) THEN
          IKPTL=IKPTL+1
          CALL WAVES_SELECTWV(IKPTL,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETI4('NGG',NG)
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT),1,NG)
        CALL REPORT$I4VAL(NFIL &
     &       ,'NUMBER OF (GLOBAL) PLANE WAVES FOR WAVE FUNCTION',NG,' ')
      ENDDO
!
!!$      DO IKPT=1,NKPT
!!$        WRITE(STRING,FMT='("KPOINT NR.",I5," RESIDES ON TASK ")')IKPT
!!$        CALL REPORT$I4VAL(NFIL,STRING,KMAP(IKPT),' ')
!!$      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$WRITEPDOS
!     ******************************************************************
!     **  THIS FUNCTION WRITES THE PDOS FILE                          **
!     ******************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: NFIL
      REAL(8)   ,ALLOCATABLE :: VAL(:,:)
      REAL(8)   ,ALLOCATABLE :: DER(:,:)
      REAL(8)   ,ALLOCATABLE :: OV(:,:,:)
      REAL(8)   ,ALLOCATABLE :: R(:,:)
      REAL(8)   ,ALLOCATABLE :: RAD(:)
      REAL(8)   ,ALLOCATABLE :: WORK(:)
      REAL(8)   ,ALLOCATABLE :: WORK1(:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)                :: AEZ
      INTEGER(4),ALLOCATABLE :: IZ(:)
      INTEGER(4)             :: NAT,NSP
      INTEGER(4)             :: NR
      INTEGER(4)             :: NB,NBH
      REAL(8)   ,ALLOCATABLE :: XK(:,:)   ! K-POINT IN RELATIVE COORDINATES
      REAL(8)   ,ALLOCATABLE :: WKPT(:)   ! K-POINT WEIGHTS
      COMPLEX(8),ALLOCATABLE :: VEC(:,:,:)      
      COMPLEX(8),ALLOCATABLE :: VECTOR1(:,:,:)
      REAL(8)   ,ALLOCATABLE :: EIG(:)
      INTEGER(4)             :: LNXX,LNX,NPRO
      INTEGER(4)             :: NGL
      INTEGER(4),ALLOCATABLE :: LOX(:)
      LOGICAL(4)             :: TINV
      REAL(8)                :: RBAS(3,3)
      INTEGER(4)             :: ISP
      INTEGER(4)             :: NTASKS,THISTASK
      INTEGER(4)             :: IB1,IB2,IBH,LN1,LN2,IDIM,IKPT,ISPIN,IPRO
      INTEGER(4)             :: IKPTG
      COMPLEX(8),ALLOCATABLE :: PROJ(:,:,:)
      CHARACTER(16),ALLOCATABLE :: ATOMID(:)
      INTEGER(4)             :: NBX
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      LOGICAL(4)             :: TKGROUP
      INTEGER(4)             :: GID
      REAL(8)   ,ALLOCATABLE :: RI(:)
      INTEGER(4)             :: NKDIV(3)
      INTEGER(4)             :: ISHIFT(3)
      REAL(8)                :: NEL
      REAL(8)                :: RNTOT
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$WRITEPDOS')
      IF(.NOT.THAMILTON) THEN
        CALL ERROR$MSG('EIGENVALUES NOT PRESENT')
        CALL ERROR$STOP('WAVES$WRITEPDOS')
      END IF
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
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
      CALL DYNOCC$GETI4A('NKDIV',3,NKDIV)
      CALL PDOS$SETI4A('NKDIV',3,NKDIV)
      CALL DYNOCC$GETI4A('ISHIFT',3,ISHIFT)
      CALL PDOS$SETI4A('ISHIFT',3,ISHIFT)
      CALL DYNOCC$GETR8('NEL',NEL)
      CALL PDOS$SETR8('NEL',NEL)
      IF(NSPIN.EQ.1.AND.NDIM.EQ.1)THEN
        RNTOT=0.5D0*NEL
      ELSE
        RNTOT=NEL
      ENDIF
      CALL PDOS$SETR8('RNTOT',RNTOT)
      CALL PDOS$SETI4('NSPIN',NSPIN)
      CALL PDOS$SETI4('NDIM',NDIM)
      CALL PDOS$SETI4('NPRO',NPRO)
      CALL PDOS$SETI4('LNXX',LNXX)
      CALL PDOS$SETI4A('LNX',NSP,MAP%LNX)
      CALL PDOS$SETI4A('LOX',LNXX*NSP,MAP%LOX)
      CALL PDOS$SETI4A('ISPECIES',NAT,MAP%ISP)
      CALL KPOINTS$GETL4('TINV',TINV)
      IF(TINV)THEN
        !TRICLINIC WITH INVERSION SYMMETRY
        CALL PDOS$SETI4('SPACEGROUP',2)
      ELSE
        !TRICLINIC WITHOUT INVERSION SYMMETRY
        CALL PDOS$SETI4('SPACEGROUP',1)
      ENDIF
      CALL PDOS$SETL4('TSHIFT',.FALSE.)
      CALL PDOS$SETL4('TINV',TINV)
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
      CALL PDOS$SETCHA('ATOMID',NAT,ATOMID)
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
      OV(:,:,:)=0.D0
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL RADIAL$GETI4(GID,'NR',NR)
        ALLOCATE(RI(NR))
        CALL RADIAL$R(GID,NR,RI)
        LNX=MAP%LNX(ISP)
        LOX=MAP%LOX(:,ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        IZ(ISP)=NINT(AEZ)
        CALL SETUP$GETR8('RAD',RAD(ISP))
        ALLOCATE(AEPHI(NR,LNXX))
        CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
                            !FORMER CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,AEPHI)
        CALL SETUP$UNSELECT()
        ALLOCATE(WORK(NR))
        ALLOCATE(WORK1(NR))
        DO LN1=1,LNX
          CALL RADIAL$VALUE(GID,NR,AEPHI(1,LN1),RAD(ISP),VAL(LN1,ISP))
          CALL RADIAL$DERIVATIVE(GID,NR,AEPHI(1,LN1),RAD(ISP),DER(LN1,ISP))
          DO LN2=LN1,LNX
            IF(LOX(LN1).NE.LOX(LN2)) CYCLE
            WORK(:)=AEPHI(:,LN1)*AEPHI(:,LN2)*RI(:)**2
            CALL RADIAL$INTEGRATE(GID,NR,WORK,WORK1)
            CALL RADIAL$VALUE(GID,NR,WORK1,RAD(ISP),OV(LN1,LN2,ISP))
            OV(LN2,LN1,ISP)=OV(LN1,LN2,ISP)
          ENDDO
        ENDDO
        DEALLOCATE(RI)
        DEALLOCATE(WORK)
        DEALLOCATE(WORK1)
        DEALLOCATE(AEPHI)
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
      DEALLOCATE(LOX)

      IF(THISTASK.EQ.1) THEN
        CALL FILEHANDLER$UNIT('PDOS',NFIL)
        REWIND NFIL
        CALL PDOS$WRITE(NFIL,'181213')
      ENDIF
!
!     ==================================================================
!     ==  NOW WRITE PROJECTIONS                                       ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPT,NSPIN))
      CALL DYNOCC$GETR8A('OCC',NBX*NKPT*NSPIN,OCC)
      ALLOCATE(XK(3,NKPT))
      CALL DYNOCC$GETR8A('XK',3*NKPT,XK)
      ALLOCATE(WKPT(NKPT))
      CALL DYNOCC$GETR8A('WKPT',NKPT,WKPT)

      CALL PDOS$SETR8A('XK',3*NKPT,XK)

      IKPT=0    ! IKPT IS THE K-POINT INDEX LOCAL TO MPE-GROUP K
      DO IKPTG=1,NKPT
        TKGROUP=THISTASK.EQ.KMAP(IKPTG)
        CALL MPE$BROADCAST('K',1,TKGROUP)
        IF(.NOT.(THISTASK.EQ.1.OR.TKGROUP)) CYCLE
        IF(TKGROUP)IKPT=IKPT+1
        DO ISPIN=1,NSPIN
          IF(TKGROUP) THEN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            CALL PLANEWAVE$SELECT(GSET%ID)
            IF(.NOT.ASSOCIATED(THIS%EIGVEC)) THEN
              CALL ERROR$MSG('EIGENVALUES NOT PRESENT')
              CALL ERROR$STOP('WAVES$WRITEPDOS')
            END IF
            NB=THIS%NB
            NBH=THIS%NBH
            NGL=GSET%NGL
            TINV=GSET%TINV
!  
!           =============================================================
!           ==  CALCULATE PROJECTIONS                                  ==
!           =============================================================
            ALLOCATE(PROJ(NDIM,NBH,NPRO))
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NBH,NPRO &
       &                          ,THIS%PSI0,PROJ)
            CALL MPE$COMBINE('K','+',PROJ)
          END IF
          IF(KMAP(IKPTG).EQ.THISTASK) THEN
            ALLOCATE(VECTOR1(NDIM,NPRO,NB))
!
!           ==============================================================
!           == UNRAVEL SUPER WAVE FUNCTIONS                             ==
!           ==============================================================
            IF(TINV) THEN
              DO IBH=1,NBH
                IB1=1+2*(IBH-1)
                IB2=2+2*(IBH-1)
                DO IPRO=1,NPRO
                  DO IDIM=1,NDIM
                    VECTOR1(IDIM,IPRO,IB1)=REAL(PROJ(IDIM,IBH,IPRO),KIND=8)
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
!
!           ==============================================================
!           ==  TRANSFORM TO EIGENSTATES IF SAFEORTHO=.TRUE.,           ==
!           ==============================================================
            ALLOCATE(EIG(NB))
            IF(TSAFEORTHO) THEN
!             == CONSTRUCT EIGENSTATES ===================================
              ALLOCATE(VEC(NDIM,NPRO,NB))
              VEC(:,:,:)=0.D0
              DO IB1=1,NB
                DO IB2=1,NB
                  VEC(:,:,IB1)=VEC(:,:,IB1) &
      &                       +VECTOR1(:,:,IB2)*THIS%EIGVEC(IB2,IB1)
                ENDDO
              ENDDO
              VECTOR1(:,:,:)=VEC(:,:,:)
              DEALLOCATE(VEC)
              EIG(:)=THIS%EIGVAL(:)
            ELSE
              EIG(:)=THIS%EXPECTVAL(:)
            END IF
          END IF
          IF(TKGROUP) DEALLOCATE(PROJ)
!
!         ================================================================
!         == COMMUNICATE RESULT TO FIRST NODE OF MONOMER                ==
!         == NB MAY IN FUTURE DEPEND ON THE K-POINT                     ==
!         == NDIM AND NPRO ARE INDEPENDENT OF THE K-POINTS              ==
!         ================================================================
          IF(KMAP(IKPTG).NE.1) THEN  !ONLY COMMUNICATE WHEN NECESSARY
            CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,NB)
            IF(THISTASK.EQ.1) THEN
              ALLOCATE(EIG(NB))
              ALLOCATE(VECTOR1(NDIM,NPRO,NB))
            END IF
            CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,EIG)
            CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,VECTOR1)
            IF(THISTASK.EQ.KMAP(IKPTG)) THEN
              DEALLOCATE(EIG)
              DEALLOCATE(VECTOR1)
            END IF
          END IF
!
!         ======================================================================
!         == NOW THE DATA IS AVAILABLE ON TASK 1 OF MONOMER; NEXT WRITE       ==
!         == TO FILE. ATTENTION: DATA ARE ONLY AVAILABLE ON TASK 1.           ==
!         == THIS BLOCK CANNOT BE INTEGRATED INTO THE COMMUNICATION           ==
!         == NBLOCK ABOVE, BECAUSE THERE IS NO COMMUNICATION IF THE DATA      ==
!         == ARE ALREADY IN TASK 1.                                           ==
!         ======================================================================
          IF(THISTASK.EQ.1) THEN
            CALL PDOS$WRITEK(NFIL,XK(:,IKPTG),NB,NDIM,NPRO,&
      &              WKPT(IKPTG),EIG,OCC(:,IKPTG,ISPIN),VECTOR1(:,:,:))
            DEALLOCATE(EIG)
            DEALLOCATE(VECTOR1)
          END IF
        ENDDO  ! END LOOP OVER ISPIN
      ENDDO  ! END LOOP OVER IKPTG
      DEALLOCATE(XK)
      DEALLOCATE(WKPT)
      DEALLOCATE(R)
      DEALLOCATE(OCC)
      IF(THISTASK.EQ.1) THEN
        FLUSH(NFIL)
        CALL FILEHANDLER$CLOSE('PDOS')
      END IF
                             CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$REPORTEIG(NFIL)
!     **************************************************************************
!     **  PRODUCES A PRINTOUT OF ONE-PARTICLE ENERGIES AS WELL AS INFORMATION **
!     **  ON HOMO/LUMO ETC.                                                   **
!     **                                                                      **
!     **  REMARK: THIS ROUTINE IS ONLY USED WITH FIXED OCCUPATIONS            **
!     **         , I.E. WHEN DYNOCC$GETL4('DYN',TCHK) RETURNS TCHK=.FALSE.    **
!     **************************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: IKPT,IKPTL
      INTEGER(4)            :: ISPIN
      INTEGER(4)            :: IB
      INTEGER(4)            :: ITEN
      REAL(8)               :: EV
      CHARACTER(64)         :: STRING
      INTEGER(4)            :: NB
      REAL(8)  ,ALLOCATABLE :: EIG(:)
      INTEGER(4)            :: NTASKS,THISTASK
      REAL(8)  ,ALLOCATABLE :: OCC(:,:,:)
      REAL(8)  ,ALLOCATABLE :: WKPT(:)
      REAL(8)               :: ESPINHOMO(NSPIN)
      REAL(8)               :: ESPINLUMO(NSPIN)
      INTEGER(4)            :: IBSPINHOMO(NSPIN)
      INTEGER(4)            :: IBSPINLUMO(NSPIN)
      REAL(8)               :: EHOMO,ELUMO,EGDIRECT
      INTEGER(4)            :: IBHOMO,IKHOMO,IKLUMO,IKDIRECT
      INTEGER(4)            :: IBLUMO,ISHOMO,ISLUMO,ISDIRECT
      CHARACTER(256)        :: FORMAT
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      CALL CONSTANTS('EV',EV)

!
!     ==========================================================================
!     == DETERMINE TOTAL NUMBER OF ELECTRONS TO ESTIMATE BAND GAP             ==
!     ==========================================================================
      CALL DYNOCC$GETI4('NB',NB)
      ALLOCATE(OCC(NB,NKPT,NSPIN))
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
      ALLOCATE(WKPT(NKPT))
      CALL DYNOCC$GETR8A('WKPT',NKPT,WKPT) !K-POINT WEIGHT
      EGDIRECT    =+1.D+10
      EHOMO       =-1.D+10
      ELUMO       =+1.D+10
      ESPINHOMO(:)=-1.D+10
      ESPINLUMO(:)=+1.D+10
      IBHOMO       =-1
      IBLUMO       =-1
      IBSPINHOMO(:)=-1
      IBSPINLUMO(:)=-1
      IKDIRECT=0
      IKHOMO=0
      IKLUMO=0
      ISDIRECT=0
      ISHOMO=0
      ISLUMO=0
!
!     ==========================================================================
!     == WRITE ENERGY BANDS                                                   ==
!     ==========================================================================
      IKPTL=0
      DO IKPT=1,NKPT
        IF(KMAP(IKPT).EQ.THISTASK) IKPTL=IKPTL+1
        IF(THISTASK.NE.1.AND.THISTASK.NE.KMAP(IKPT)) CYCLE
        DO ISPIN=1,NSPIN
          IF(THISTASK.EQ.KMAP(IKPT)) THEN
            CALL WAVES_SELECTWV(IKPTL,ISPIN)
            CALL PLANEWAVE$SELECT(GSET%ID)
            IF(.NOT.ASSOCIATED(THIS%EIGVAL)) THEN
              CALL ERROR$MSG('EIGENVALUES NOT PRESENT')
              CALL ERROR$STOP('WAVES$REPORTEIG')
            END IF
            NB=THIS%NB
          END IF
          CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT),1,NB)
          ALLOCATE(EIG(NB))
          IF(THISTASK.EQ.KMAP(IKPT)) THEN
            EIG(:)=THIS%EIGVAL(:)
          END IF
          CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT),1,EIG)
          IF(THISTASK.EQ.1) THEN
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
     &                               ITEN,(EIG(IB)/EV,IB=ITEN+1,MIN(ITEN+10,NB))
              ITEN=ITEN+10
            ENDDO
!
!           == SCAN EIGENVALUES FOR BAND GAPS ==================================
            DO IB=1,NB
!!$              IF(OCC(IB,IKPT,ISPIN).GT.1.D-6.AND.IB.GT.IBHOMO(ISPIN)) THEN
!!$                IBHOMO(ISPIN)=IB
!!$              END IF
!             __LOOK FOR HOMO___________________________________________________
              IF(OCC(IB,IKPT,ISPIN).GT.1.D-6.AND.EIG(IB).GT.EHOMO) THEN
                EHOMO=EIG(IB)
                IKHOMO=IKPT
                ISHOMO=ISPIN
                IBHOMO=IB
              END IF
!             __LOOK FOR HOMO___________________________________________________
              IF(OCC(IB,IKPT,ISPIN).GT.1.D-6 &
         &                                .AND.EIG(IB).GT.ESPINHOMO(ISPIN)) THEN
                ESPINHOMO(ISPIN)=EIG(IB)
                IBSPINHOMO(ISPIN)=IB
              END IF
!             __LOOK FOR LUMO___________________________________________________
              IF(OCC(IB,IKPT,ISPIN).LT.1.D-6.AND.EIG(IB).LT.ELUMO) THEN
                ELUMO=EIG(IB)
                IKLUMO=IKPT
                ISLUMO=ISPIN
                IBLUMO=IB
              END IF
!             __LOOK FOR LUMO___________________________________________________
              IF(OCC(IB,IKPT,ISPIN).LT.1.D-6 &
        &                                 .AND.EIG(IB).LT.ESPINLUMO(ISPIN)) THEN
                ESPINLUMO(ISPIN)=EIG(IB)
                IBSPINLUMO(ISPIN)=IB
              END IF
!             __LOOK FOR DIRECT GAP_____________________________________________
              IF(IB.GT.1.AND.OCC(IB,IKPT,ISPIN).LT.1.D-6) THEN
                IF(OCC(IB-1,IKPT,ISPIN).GT.1.D-6) THEN
                  IF(EIG(IB)-EIG(IB-1).LT.EGDIRECT) THEN
                    EGDIRECT=EIG(IB)-EIG(IB-1)
                    IKDIRECT=IKPT
                    ISDIRECT=ISPIN
                  END IF
                END IF
              END IF
            ENDDO   !END OF LOOP OVER BANDS
          END IF
          DEALLOCATE(EIG)
        ENDDO ! END OF LOOP OVER SPINS
      ENDDO   ! END OF LOOP OVER K-POINTS
!
!     ==========================================================================
!     == REPORT ENERGY GAPS                                                   ==
!     ==========================================================================
      IF(NSPIN.EQ.1) THEN
        CALL REPORT$I4VAL(NFIL,'BAND INDEX OF HOMO',IBSPINHOMO(1),' ')
      ELSE     
        CALL REPORT$I4VAL(NFIL,'BAND INDEX OF HOMO FOR SPIN 1' &
     &                                                       ,IBSPINHOMO(1),' ')
        CALL REPORT$I4VAL(NFIL,'BAND INDEX OF HOMO FOR SPIN 2' &
     &                                                       ,IBSPINHOMO(2),' ')
      END IF
      FORMAT='(55("."),": ",T1,A,T58,F10.4," EV AT IK=",I5," AND ISPIN=",I1)'
      WRITE(NFIL,FMT=FORMAT)'SMALLEST DIRECT GAP',EGDIRECT/EV,IKDIRECT,ISDIRECT 
      FORMAT='(55("."),": ",T1,A,T58,F10.4," EV ")'
      WRITE(NFIL,FMT=FORMAT)'ABSOLUTE GAP',(ELUMO-EHOMO)/EV
      FORMAT=       '(T10," FROM BAND=",I5," AT IK=",I5," WITH ISPIN=",I1'
      FORMAT=TRIM(FORMAT)//'" TO BAND ",I5," AT IK=",I5," WITH ISPIN=",I1)'
      WRITE(NFIL,FMT=FORMAT)IBHOMO,IKHOMO,ISHOMO,IBLUMO,IKLUMO,ISLUMO
      CALL REPORT$R8VAL(NFIL,'HOMO-ENERGY',EHOMO/EV,'EV')
      CALL REPORT$R8VAL(NFIL,'LUMO-ENERGY',ELUMO/EV,'EV')
      IF(ELUMO-EHOMO.LT.0.D0) THEN
        CALL REPORT$STRING(NFIL,'MATERIAL IS A METAL: USE VARIABLE OCCUPATIONS')
      END IF
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
!     ******************************************************************
              CALL TRACE$PUSH('WAVES$WRITE')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     == WRITE WAVE FUNCTION FOR THE THIS TIME STEP ==================
      TCHK=.TRUE.
      SEPARATOR=SEP_WAVES
      SEPARATOR%NREC=-1
      CALL WAVES_WRITEPSI(NFIL,SEPARATOR%NREC)
      IF(THISTASK.EQ.1)CALL RESTART$WRITESEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
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
      REAL(8)              ,ALLOCATABLE:: EIG(:,:,:)
      INTEGER(4)                       :: IKPT,ISPIN,IB
      INTEGER(4)                       :: NBX
      INTEGER(4)                       :: THISTASK,NTASKS
!     ******************************************************************
                                  CALL TRACE$PUSH('WAVES$READ')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     ==================================================================
!     ==  READ WAVE FUNCTIONS FOR THIS TIME STEP                      ==
!     ==================================================================
      TCHK=.TRUE.
      SEPARATOR=SEP_WAVES
      IF(THISTASK.EQ.1)CALL RESTART$READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      CALL MPE$BROADCAST('MONOMER',1,TCHK)
      IF(.NOT.TCHK) THEN
                                  CALL TRACE$POP
        RETURN
      END IF
!
!     ==================================================================
!     == UPDATE GSET TO ACCOUNT FOR A CHANGE OF THE LATTICE VECTORS 
!     == BY A CELL$READ. ENSURE THAT CELL$READ IS CALLED BEFORE WAVES$READ!
!     ==================================================================
      CALL WAVES_UPDATEGSET()
!
!     ==================================================================
!     == READ WAVE FUNCTIONS AND LAGRANGE PARAMETERS                  ==
!     ==================================================================
      CALL WAVES_READPSI(NFIL)
!      CALL WAVES_READPSI_OLD(NFIL)
!
!     ==================================================================
!     == SET OCCUPATIONS FOR DYNOCC OBJECT                            ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(EIG(NBX,NKPTL,NSPIN))
      EIG(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPTL
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          DO IB=1,THIS%NB
            EIG(IB,IKPT,ISPIN)=REAL(THIS%RLAM0(IB,IB),KIND=8)
          ENDDO
        ENDDO
      ENDDO
      CALL WAVES_DYNOCCSETR8A('EPSILON',NBX*NKPTL*NSPIN,EIG)
      DEALLOCATE(EIG)
!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
                                  CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_WRITEPSI(NFIL,NREC)
!     **************************************************************************
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                                **
!     **                                                                      **
!     **  CALL FIRST WITH NREC=-1. ON OUTPUT NREC IS THE NUMBER OF RECORDS.    **
!     **  CALL WITH CORRECT NREC TO WRITE                                     **
!     **                                                                      **
!     **************************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL 
      INTEGER(4)  ,INTENT(INOUT):: NREC
      COMPLEX(8)  ,ALLOCATABLE:: PSIG(:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSI1(:,:)
      INTEGER(4)              :: NTASKS,THISTASK
      INTEGER(4)              :: IOS
      INTEGER(4)              :: IKPT,IKPTG,ISPIN,IB,IDIM,IWAVE
      INTEGER(4)              :: NGG,NGL,NBH,NB
      LOGICAL(4)              :: TSUPER
      CHARACTER(64)           :: IOSTATMSG
      REAL(8)                 :: XK(3)
      REAL(8)                 :: GBAS(3,3)
      INTEGER(4)              :: NREC1,ISVAR
      CHARACTER(8)            :: KEY
      REAL(8)                 :: RBAS(3,3)
      INTEGER(4)  ,ALLOCATABLE:: IGVECG(:,:)
      INTEGER(4)  ,ALLOCATABLE:: IGVECL(:,:)
      INTEGER(4)              :: NWAVE=2
      INTEGER(4)              :: ISUM
      LOGICAL(4)              :: TKGROUP
      INTEGER(4)              :: ILOGICAL
!     **************************************************************************
              CALL TRACE$PUSH('WAVES_WRITEPSI')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(RSTRTTYPE.EQ.'STATIC') THEN
        NWAVE=1
      ELSE IF(RSTRTTYPE.EQ.'DYNAMIC') THEN
        NWAVE=2
      ELSE
        CALL ERROR$STOP('WAVES_WRITEPSI')
      END IF
!
!     ==========================================================================
!     ==  COUNT #(RECORDS)                                                    ==
!     ==========================================================================
      IF(NREC.EQ.-1) THEN
        NREC=1
        ISUM=0
        IKPT=0
        DO IKPTG=1,NKPT
          IF(THISTASK.EQ.KMAP(IKPTG)) THEN
            IKPT=IKPT+1
            CALL WAVES_SELECTWV(IKPT,1)
            CALL PLANEWAVE$SELECT(GSET%ID)
            ISUM=ISUM+2+NWAVE*THIS%NBH*NSPIN
!           == ISVAR=-1: LET WAVES$WRITELAMBDA WORK OUT HOW MANY RECORDS IT NEEDS
            ISVAR=-1
            CALL WAVES_WRITELAMBDA(NFIL,IKPTG,ISVAR)
            ISUM=ISUM+ISVAR
          END IF
        ENDDO  
        CALL MPE$COMBINE('MONOMER','+',ISUM)
        NREC=NREC+ISUM
        CALL TRACE$POP
        RETURN
      END IF
      NREC1=0
!
!     ==========================================================================
!     ==  WRITE SIZES                                                         ==
!     ==========================================================================
      IF(THISTASK.EQ.1) THEN
        CALL CELL$GETR8A('T(0)',9,RBAS)
        WRITE(NFIL)NKPT,NSPIN,RBAS,NWAVE !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        NREC1=NREC1+1
      END IF
!
!     ==========================================================================
!     ==  LOOP OVER K-POINTS AND SPINS                                        ==
!     ==========================================================================
      IKPT=0
      DO IKPTG=1,NKPT
        TKGROUP=(THISTASK.EQ.KMAP(IKPTG))
        CALL MPE$BROADCAST('K',1,TKGROUP) 
        IF(.NOT.(THISTASK.EQ.1.OR.TKGROUP)) CYCLE
        IF(TKGROUP) THEN
          IKPT=IKPT+1
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NBH=THIS%NBH
          NB=THIS%NB
          CALL PLANEWAVE$GETI4('NGG',NGG)
          CALL PLANEWAVE$GETL4('TINV',TSUPER)
        END IF
!       
!       ========================================================================
!       == WRITE SIZE AND TYPE OF WAVE FUNCTION                               ==
!       ========================================================================
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,NGG)
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,NB)
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,NBH)
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,TSUPER)
        IF(THISTASK.EQ.1) THEN
          KEY='PSI'
!         == RSTRT FILE WITH GNU STANDARD OF LOGICAL BIT REPRESENTATION ==
          ILOGICAL=1
          IF(.NOT.TSUPER) ILOGICAL=0
          WRITE(NFIL)KEY,NGG,NDIM,NB,ILOGICAL  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC1=NREC1+1
        END IF
!       
!       ========================================================================
!       == WRITE K-POINTS AND G-VECTORS                                       ==
!       ========================================================================
        ALLOCATE(IGVECG(3,NGG))
        IF(TKGROUP) THEN
          CALL PLANEWAVE$GETR8A('GBAS',9,GBAS)
          CALL PLANEWAVE$GETR8A('XK',3,XK)
          ALLOCATE(IGVECL(3,NGL))
          CALL PLANEWAVE$GETI4A('IGVEC',3*NGL,IGVECL)
          CALL PLANEWAVE$COLLECTI4(3,NGL,IGVECL,NGG,IGVECG)
          DEALLOCATE(IGVECL)
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,XK)
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,IGVECG)
        IF(THISTASK.EQ.1) THEN
          WRITE(NFIL)XK,IGVECG !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC1=NREC1+1
        END IF
        DEALLOCATE(IGVECG)
!
        ALLOCATE(PSIG(NGG,NDIM))
        IF(TKGROUP) ALLOCATE(PSI1(NGL,NDIM))
        DO IWAVE=1,NWAVE
          DO ISPIN=1,NSPIN
            DO IB=1,NBH
              IF(TKGROUP) THEN
                CALL WAVES_SELECTWV(IKPT,ISPIN)
                CALL PLANEWAVE$SELECT(GSET%ID)
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
              END IF
              CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,PSIG)
              IF(THISTASK.EQ.1) THEN
                WRITE(NFIL,ERR=9999,IOSTAT=IOS)PSIG !<<<<<<<<<<<<<<<<<<<<<<<<<<<
                NREC1=NREC1+1
              END IF
            ENDDO
          ENDDO
        ENDDO
        IF(TKGROUP) DEALLOCATE(PSI1)
        DEALLOCATE(PSIG)
        ISVAR=0
        CALL WAVES_WRITELAMBDA(NFIL,IKPTG,ISVAR)
        NREC1=NREC1+ISVAR
      ENDDO 
!
!     ==========================================================================
!     ==  DEALLOCATE                                                          ==
!     ==========================================================================
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
      SUBROUTINE WAVES_READPSI_OLD(NFIL)
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
      INTEGER(4)              :: IKPT0,IKPTG,IKPTL
      INTEGER(4)              :: IWAVE
      LOGICAL(4)              :: TREAD(NKPT)
      REAL(8)                 :: KREAD(3,NKPT)
      REAL(8)     ,ALLOCATABLE:: XK(:,:)     ! K-POINTS IN RELATIVE COORDINATES
      REAL(8)                 :: DK(3)       ! K-POINT DIFFERENCE IN A.U.
      REAL(8)                 :: K_(3)  !K-POINT ON FILE IN RELATIVE COORDINATES
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
      COMPLEX(8)  ,ALLOCATABLE:: PSITMP(:,:)
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)              :: IOS
      CHARACTER(8)            :: KEY
      LOGICAL(4)              :: GBASFIX
      INTEGER(4)              :: IFORMAT
      INTEGER(4) ,ALLOCATABLE :: IGVECG_(:,:)
      REAL(8)                 :: XG1,XG2,XG3
      INTEGER(4)              :: NWAVE,NFILO
      INTEGER(4) ,ALLOCATABLE :: MINUSG(:)
      LOGICAL(4)              :: TCHK
      COMPLEX(8)              :: F1,F2
      COMPLEX(8),ALLOCATABLE  :: PSIINSUPER(:,:)
      COMPLEX(8)              :: CSVAR,CMAT(1,1)
      LOGICAL(4)              :: TKGROUP
      LOGICAL(4)              :: TSKIP
      INTEGER(4)              :: ILOGICAL
!     ******************************************************************
                               CALL TRACE$PUSH('WAVES_READPSI')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==================================================================
!     ==  READ SIZES                                                  ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
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
        IF(IFORMAT.NE.2) THEN
          CALL ERROR$MSG('OLD FORMATS ARE NO MORE SUPPORTED')
          CALL ERROR$STOP('WAVES_READPSI')
        END IF
      END IF
      CALL MPE$BROADCAST('MONOMER',1,NSPIN_)
      CALL MPE$BROADCAST('MONOMER',1,NKPT_)
      CALL MPE$BROADCAST('MONOMER',1,IFORMAT)
      CALL MPE$BROADCAST('MONOMER',1,NWAVE)
      IF(IFORMAT.EQ.2) THEN
        CALL MPE$BROADCAST('MONOMER',1,RBAS)
        DO IKPT=1,NKPTL
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
      ALLOCATE(XK(3,NKPT))
      CALL DYNOCC$GETR8A('XK',3*NKPT,XK) !K-POINTS IN RELATIVE COORDINATES
      GBASFIX=.FALSE.
      TREAD(:)=.FALSE.
      KREAD(:,:)=0.D0
      DO IKPT_=1,NKPT_   ! LOOP OVER ALL K-POINTS ON THE FILE
!       == FIND OUT IF DATA ARE TO BE READ AND TO WHICH K-GROUP THEY BELONG
!
!       ================================================================
!       ==  READ COORDINATES OF THE WAVE FUNCTIONS                    ==
!       ================================================================
        IF(THISTASK.EQ.1) THEN
!         == GFORTRAN LOGICAL REPRESENTATION DEFINED WITH TRUE=1, FALSE=0     ==
!         https://gcc.gnu.org/onlinedocs/gfortran/compiler-characteristics/
!         internal-representation-of-logical-variables.html
!         == IFORT LOGICAL REPRESENTATION DEFINED WITH VALUE OF LAST BIT      ==
!         https://www.intel.com/content/www/us/en/docs/fortran-compiler/
!         developer-guide-reference/2024-2/logical-data-representations.html
!         == BOTH SHARE MEANING OF LAST BIT 1=TRUE, 0=FALSE                   ==
!         == ENSURES BACKWARDS COMPATIBILITY WITH OLD RESTART FILES           ==
          READ(NFIL)KEY,NGG_,NDIM_,NB_,ILOGICAL   !<<<<<<<<<<<<<<<<<<<<<<
          TSUPER_=BTEST(ILOGICAL,0)
          IF(KEY.NE.'PSI') THEN
            CALL ERROR$MSG('ID IS NOT "PSI"')
            CALL ERROR$MSG('FILE IS CORRUPTED')
            CALL ERROR$I4VAL('IKPT_',IKPT_)
            CALL ERROR$STOP('WAVES_READPSI')
          END IF
          ALLOCATE(GVECG_(3,NGG_))
          IF(IFORMAT.EQ.1) THEN
            READ(NFIL)K_,GBAS_,GVECG_ !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          ELSE
            ALLOCATE(MINUSG(NGG_))
            ALLOCATE(IGVECG_(3,NGG_))
!           GBAS_ HAS BEEN CALCULATED EARLIER FOR IFROMAT=2
            READ(NFIL)K_,IGVECG_ !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            IF(TSUPER_) THEN
!             == MINUSG IS ONLY USED ON THE READING TASK. 
!             == PLANEWAVE$MINUSG IS INDEPENDENT FROM THE PLANEWAVE_MODULE
              CALL PLANEWAVE$MINUSG(K_,NGG_,IGVECG_,MINUSG)
            END IF
            DO IG=1,NGG_
              XG1=REAL(IGVECG_(1,IG),KIND=8)+K_(1)
              XG2=REAL(IGVECG_(2,IG),KIND=8)+K_(2)
              XG3=REAL(IGVECG_(3,IG),KIND=8)+K_(3)
              DO I=1,3
                GVECG_(I,IG)=GBAS_(I,1)*XG1+GBAS_(I,2)*XG2+GBAS_(I,3)*XG3
              ENDDO
            ENDDO
            DEALLOCATE(IGVECG_)
          END IF
        END IF
!
!       ==================================================================
!       == SET THE NEW LATTICE VECTORS ON ALL TASKS ======================
!       == THIS SHOULD BETTER BE DONE VIA THE ATOM OR THE CELL OBJECT   ==
!       ==================================================================
        IF(IFORMAT.EQ.1.AND.(.NOT.GBASFIX)) THEN
          GBASFIX=.TRUE. 
          CALL MPE$BROADCAST('MONOMER',1,GBAS_)
          CALL GBASS(GBAS_,RBAS,SVAR)
          DO IKPT1=1,NKPTL
            CALL WAVES_SELECTWV(IKPT1,1)
            CALL PLANEWAVE$SELECT(GSET%ID)
            CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
          ENDDO
        END IF
!       
!       ==================================================================
!       ==  FIND NEAREST K-POINT                                        ==
!       ==================================================================
        IF(THISTASK.EQ.1) THEN
          SVAR=1.D+10
          DO IKPT1=1,NKPT
            DK(:)=MATMUL(GBAS_,(XK(:,IKPT1)-K_(:)))
            SVAR1=DOT_PRODUCT(DK,DK)
            IF(SVAR1.LT.SVAR) THEN
              IKPTG=IKPT1
              SVAR=SVAR1
            END IF
          ENDDO
!       
!         ==================================================================
!         ==  CHECK IF THERE IS NOT ANOTHER K-POINT, THAT HAS ALREADY BEEN==
!         ==  READ IN INTO THE SELECTED ARRAY AND THAT IS CLOSER          ==
!         ==================================================================
          TSKIP=.FALSE.
          IF(TREAD(IKPTG)) THEN
            DK(:)=MATMUL(GBAS_,(XK(:,IKPTG)-K_(:)))
            SVAR=DOT_PRODUCT(DK,DK)
            DK(:)=MATMUL(GBAS_,(XK(:,IKPTG)-KREAD(:,IKPTG)))
            SVAR1=SVAR-DOT_PRODUCT(DK,DK)
            IF(SVAR1.GT.0.D0) THEN  ! PREVIOUS CHOICE WAS BETTER
              NBH_=NB_
              IF(TSUPER_)NBH_=INT(0.5D0*REAL(NB_+1,KIND=8))
              DO IWAVE=1,NWAVE
                DO ISPIN_=1,NSPIN_
                  DO IB=1,NBH_
                    IF(THISTASK.EQ.1)READ(NFIL)
                  ENDDO
                ENDDO
              ENDDO
              CALL WAVES_READLAMBDA(NFIL,IKPTG,.TRUE.)
              TSKIP=.TRUE.
            END IF
          END IF
          TREAD(IKPTG)=.TRUE.
          KREAD(:,IKPTG)=K_(:)
        END IF
        CALL MPE$BROADCAST('MONOMER',1,TSKIP)
!
!       == DETERMINE THE K-GROUP WHERE THE K-POINT IKPTG RESIDES       ==
        IF(.NOT.TSKIP) THEN
          CALL MPE$BROADCAST('MONOMER',1,IKPTG)
          TKGROUP=THISTASK.EQ.KMAP(IKPTG)
          CALL MPE$BROADCAST('K',1,TKGROUP)
        ELSE
          TKGROUP=.FALSE.   ! NOT NEEDED
          IKPTG=0
        END IF
!
!       ===================================================================
!       ==  SKIP REST OF THE LOOP,                                       ==
!       ==  FOR ALL TASKS IF THERE ARE ALREADY BETTER DATA FOR THIS K-POINT
!       ==  OR, IF DATA ARE READ, ON ALL TASKS THAT ARE NOT INVOLVED     ==
!       ===================================================================
!CALL MPE$REPORT(NFILO)
CALL MPE$SYNC('MONOMER')
        IF(TSKIP) THEN
          IF(ALLOCATED(GVECG_))DEALLOCATE(GVECG_)
          IF(ALLOCATED(MINUSG))DEALLOCATE(MINUSG)
          CYCLE   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<CYCLE! HERE<<<<
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',1,KMAP(IKPTG),NDIM_)
        CALL MPE$SENDRECEIVE('MONOMER',1,KMAP(IKPTG),NB_)
        IF(TKGROUP) THEN
          CALL MPE$BROADCAST('K',1,NDIM_)
          CALL MPE$BROADCAST('K',1,NB_)
          CALL MPE$BROADCAST('K',1,NSPIN_)
          CALL MPE$BROADCAST('K',1,NWAVE)
        END IF
!
!       == DETERMINE LOKAL K-POINT ===================================        
        IF(TKGROUP) THEN
          IKPTL=0
          DO I=1,IKPTG
            IF(KMAP(I).EQ.KMAP(IKPTG)) IKPTL=IKPTL+1
          ENDDO            
        ELSE
          IKPTL=0 ! MAKE SURE THE FIRST MONOMER-TASK CANNOT USE IT
        END IF
!
!       == SET DATA FOR FURTHER USE ==================================
        NGG=1
        IF(TKGROUP) THEN
          CALL WAVES_SELECTWV(IKPTL,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          CALL PLANEWAVE$GETI4('NGG',NGG)
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,NGG)
!       
!       ==============================================================
!       ==  DEFINE MAPPING OF THE ARRAYS                            ==
!       ==============================================================
        ALLOCATE(GVECG(3,NGG))
        IF(TKGROUP) THEN
          ALLOCATE(GVECL(3,NGL))
          CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVECL)
          CALL PLANEWAVE$COLLECTR8(3,NGL,GVECL,NGG,GVECG)
          DEALLOCATE(GVECL)
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,GVECG)
        IF(THISTASK.EQ.1) THEN
          ALLOCATE(MAPG(NGG))
          CALL WAVES_MAPG(NGG_,GBAS_,GVECG_,NGG,GVECG,MAPG)
        END IF
        IF(THISTASK.EQ.1)DEALLOCATE(GVECG_)
        DEALLOCATE(GVECG)
!       --------------------------------------------------------------
!       == THE FIRST TASK KNOWS:
!       ==   NGG,MAPG,MINUSG
!       ==   NGG_,NDIM_,NB_,TSUPER_
!       --------------------------------------------------------------
!       
!       ==============================================================
!       ==  COLLECT DATA                                            ==
!       ==============================================================
IF(THISTASK.EQ.1) THEN
  IF(TSUPER_.AND.NDIM_.EQ.2) THEN
    CALL ERROR$MSG('NONCOLLINEAR WAVE FUNCTIONS AND SUPER WAVE FUNCTIONS')
    CALL ERROR$MSG('ARE NOT PROPERLY IMPLEMENTED BELOW')
    CALL ERROR$STOP('WAVES_READPSI')
  END IF
END IF
              CALL TRACE$PASS('READPSI MARKE 9')
        ALLOCATE(PSIG(NGG,NDIM_))
        IF(TKGROUP) THEN
          NBH=THIS%NBH
          NB=THIS%NB
          ALLOCATE(PSIL(NGL,NDIM_,NB_,NSPIN_))
          ALLOCATE(PSI(NGL,NDIM,NB,NSPIN))
        END IF
        IF(THISTASK.EQ.1) THEN
          ALLOCATE(PSIIN(NGG_,NDIM_))
          IF(TSUPER_) ALLOCATE(PSIINSUPER(NGG_,NDIM_))
        END IF
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
!THIS HAS BEEN CORRECTED 9.NOV2003: (- SIGN INCLUDED). PEB
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
              CALL MPE$SENDRECEIVE('MONOMER',1,KMAP(IKPTG),PSIG)
              IF(TKGROUP) THEN
                DO IDIM=1,NDIM_
                  CALL PLANEWAVE$DISTRIBUTEC8(1,NGG,PSIG(1,IDIM),NGL,PSIL(1,IDIM,IB,ISPIN))
                ENDDO
              END IF
            ENDDO
          ENDDO
!       
!         ==============================================================
!         ==  TRANSFORM INTO THE CORRECT FORMAT                       ==
!         ==============================================================
!         == SAME REPRESENTATION =========================================
          CALL TRACE$PASS('MARKE 9')
          IF(TKGROUP) THEN
!           ===== FIRST MAP RANDOM WAVE FUNCTION INTO NEW FILE           
            CALL PLANEWAVE$GETL4('TINV',TSUPER)
            IF(TSUPER) THEN
              ALLOCATE(PSITMP(NGL,2))
              DO ISPIN=1,NSPIN
                CALL WAVES_SELECTWV(IKPTL,ISPIN)
                CALL PLANEWAVE$SELECT(GSET%ID)
                DO IB=1,NBH
                  IB1=2*IB-1
                  IB2=2*IB
                  DO IDIM=1,NDIM
                    PSITMP(:,1)=THIS%PSI0(:,IDIM,IB)
                    CALL PLANEWAVE$INVERTG(NGL,THIS%PSI0(:,IDIM,IB),PSITMP(:,2))
                    PSI(:,IDIM,IB1,ISPIN)= 0.5D0   *(PSITMP(:,1)+PSITMP(:,2))
                    PSI(:,IDIM,IB2,ISPIN)=-0.5D0*CI*(PSITMP(:,1)-PSITMP(:,2))
                  END DO
                ENDDO
              ENDDO
              DEALLOCATE(PSITMP)
            ELSE
              DO ISPIN=1,NSPIN
                CALL WAVES_SELECTWV(IKPTL,ISPIN)
                CALL PLANEWAVE$SELECT(GSET%ID)
                PSI(:,:,:,ISPIN)=THIS%PSI0(:,:,:)
              ENDDO
            END IF
!           ===================================================================
!           ==
!           ===================================================================
            IF(NDIM.EQ.NDIM_.AND.NSPIN.EQ.NSPIN_) THEN
              DO IB=1,NB
                IF(IB.GT.NB_) EXIT
                PSI(:,:,IB,:)=PSIL(:,:,IB,:)
              ENDDO
!           == FROM NON-SPIN POLARIZED CALCULATION =========================
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
!           == FROM SPIN POLARIZED CALCULATION ===========================
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
!           == FROM NONCOLLINEAR CALCULATION  ==============================
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
!           ==============================================================
!           ==  MAP ONTO SUPER WAVE FUNCTIONS                           ==
!           ==============================================================
            CALL PLANEWAVE$GETL4('TINV',TSUPER)
            IF(TSUPER) THEN
              DO ISPIN=1,NSPIN
!               == REMOVE THE FACTOR EXP(I*PHI) FROM THE WAVE FUNCTION
                DO IB=1,NB
                  CALL PLANEWAVE$SCALARPRODUCT('-',NGL,1,1 &
          &                   ,PSI(:,1,IB,ISPIN),1,PSI(:,1,IB,ISPIN),CMAT)
                  CSVAR=CMAT(1,1)
                  CSVAR=CONJG(SQRT(CSVAR/SQRT(CSVAR*CONJG(CSVAR))))
                  PSI(:,1,IB,ISPIN)=PSI(:,1,IB,ISPIN)*CSVAR
                ENDDO
                DO IB=1,NBH
                  IB1=2*IB-1
                  IB2=2*IB
                  PSI(:,1,IB,ISPIN)=PSI(:,1,IB1,ISPIN)+CI*PSI(:,1,IB2,ISPIN)
                ENDDO
              ENDDO
            END IF
!         
!           ==============================================================
!           ==  MAP BACK                                                ==
!           ==============================================================
            CALL TRACE$PASS('MARKE 11')
            DO ISPIN=1,NSPIN
              CALL WAVES_SELECTWV(IKPTL,ISPIN)
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
          END IF
        ENDDO 
        DEALLOCATE(PSIG)
        IF(TKGROUP) THEN
          DEALLOCATE(PSIL)
          DEALLOCATE(PSI)
        END IF
        IF(THISTASK.EQ.1) THEN
          DEALLOCATE(MAPG)
          DEALLOCATE(PSIIN)
          IF(ALLOCATED(PSIINSUPER))DEALLOCATE(PSIINSUPER)
          IF(ALLOCATED(MINUSG)) DEALLOCATE(MINUSG)
        END IF
        CALL WAVES_READLAMBDA(NFIL,IKPTG,.FALSE.)
      ENDDO
!
!     ==================================================================
!     ==  COMPLETE K-POINTS                                           ==
!     ==================================================================
      CALL MPE$BROADCAST('MONOMER',1,TREAD)
      CALL MPE$BROADCAST('MONOMER',1,KREAD)
      CALL TRACE$PASS('MARKE 14')
      DO IKPT=1,NKPT
        IF(TREAD(IKPT)) CYCLE
        SVAR=1.D+10
        IKPT0=0
        DO IKPT_=1,NKPT_
          IF(.NOT.TREAD(IKPT_)) CYCLE
          DK(:)=XK(:,IKPT)-KREAD(:,IKPT_)
          SVAR1=DOT_PRODUCT(DK,DK)
          IF(SVAR1.LT.SVAR) THEN
            IKPT0=IKPT_
            SVAR=SVAR1
          END IF
        ENDDO       
        IF(IKPT0.NE.0) THEN
          CALL WAVES_COPYPSI(IKPT0,IKPT)
          CALL WAVES_COPYLAMBDA(IKPT0,IKPT)
        END IF
      ENDDO 
!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
      DEALLOCATE(XK)
                               CALL TRACE$POP
      RETURN
 9999 CONTINUE
      CALL ERROR$MSG('ERROR WHILE READING FROM FILE')
      CALL ERROR$I4VAL('IOS',IOS)
      CALL ERROR$STOP('WAVES_READPSI')
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_WRITELAMBDA(NFIL,IKPTG,NREC)
!     **************************************************************************
!     **  WRITE LAGRANGE MULTIPLIERS  TO RESTART FILE                         **
!     **                                                                      **
!     **************************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)   :: NFIL
      INTEGER(4)  ,INTENT(IN)   :: IKPTG
      INTEGER(4)  ,INTENT(INOUT):: NREC
      INTEGER(4)                :: NTASKS,THISTASK
      INTEGER(4)                :: ISPIN,IKPTL,IKPT
      INTEGER(4)                :: NB
      CHARACTER(8)              :: KEY
      COMPLEX(8)  ,ALLOCATABLE  :: LAMBDA(:,:)
      INTEGER(4)                :: NLAMBDA
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!     == ONLY THE FORST TASK OF THE K-GROUP AND THE WRITING TASK IS INVOLVED
      IF(.NOT.(THISTASK.EQ.1.OR.THISTASK.EQ.KMAP(IKPTG))) RETURN
      IKPTL=0
      DO IKPT=1,IKPTG
        IF(KMAP(IKPT).EQ.KMAP(IKPTG))IKPTL=IKPTL+1
      ENDDO
!
!     ==========================================================================
!     == COUNT DEPTH OF HISTORY FOR LAGRANGE PARAMETERS                       ==
!     ==========================================================================
      IF(THISTASK.EQ.KMAP(IKPTG)) THEN      
        CALL WAVES_SELECTWV(IKPTL,1)
        NLAMBDA=4
        IF(.NOT.ASSOCIATED(THIS%RLAM3M)) NLAMBDA=3
        IF(.NOT.ASSOCIATED(THIS%RLAM2M)) NLAMBDA=2
        IF(.NOT.ASSOCIATED(THIS%RLAMM)) NLAMBDA=1
        IF(.NOT.ASSOCIATED(THIS%RLAM0)) NLAMBDA=0
        IF(RSTRTTYPE.EQ.'STATIC') NLAMBDA=1
      END IF
!
!     ==========================================================================
!     == RETURN #(RECORDS) OR IF ON TASK OTHER THAN THE FIRST                 ==
!     ==========================================================================
      IF(NREC.EQ.-1) THEN
!       == THE RESULT IS COMMUNICATED IN WAVES_WRITEPSI
        NREC=0
        IF(THISTASK.EQ.KMAP(IKPTG)) NREC=1+NSPIN*NLAMBDA
        RETURN
      END IF
      CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,NLAMBDA)
!
!     ==========================================================================
!     == NOW WRITE TO FILE                                                    ==
!     ==========================================================================
      NREC=0
      KEY='LAMBDA'
      IF(THISTASK.EQ.KMAP(IKPTG)) THEN      
        NB=THIS%NB
      END IF
      CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,NB)
      IF(THISTASK.EQ.1) THEN
        WRITE(NFIL)KEY,NB,NDIM,NSPIN,NLAMBDA  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        NREC=NREC+1
      END IF
      ALLOCATE(LAMBDA(NB,NB))
!
      IF(NLAMBDA.GE.1) THEN
        DO ISPIN=1,NSPIN
          IF(THISTASK.EQ.KMAP(IKPTG)) THEN      
            CALL WAVES_SELECTWV(IKPTL,ISPIN)
            LAMBDA(:,:)=THIS%RLAM0(:,:)
          END IF
          CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,LAMBDA)
          IF(THISTASK.EQ.1) THEN
            WRITE(NFIL)LAMBDA !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            NREC=NREC+1
          END IF
        ENDDO
      END IF
      IF(NLAMBDA.GE.2) THEN
        DO ISPIN=1,NSPIN
          IF(THISTASK.EQ.KMAP(IKPTG)) THEN      
            CALL WAVES_SELECTWV(IKPTL,ISPIN)
            LAMBDA(:,:)=THIS%RLAMM(:,:)
          END IF
          CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,LAMBDA)
          IF(THISTASK.EQ.1) THEN
            WRITE(NFIL)LAMBDA  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            NREC=NREC+1
          END IF
        ENDDO
      END IF
      IF(NLAMBDA.GE.3) THEN
        DO ISPIN=1,NSPIN
          IF(THISTASK.EQ.KMAP(IKPTG)) THEN      
            CALL WAVES_SELECTWV(IKPTL,ISPIN)
            LAMBDA(:,:)=THIS%RLAM2M(:,:)
          END IF
          CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,LAMBDA)
          IF(THISTASK.EQ.1) THEN
            WRITE(NFIL)LAMBDA !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            NREC=NREC+1
          END IF
        ENDDO
      END IF
      IF(NLAMBDA.GE.4) THEN
        DO ISPIN=1,NSPIN
          IF(THISTASK.EQ.KMAP(IKPTG)) THEN      
            CALL WAVES_SELECTWV(IKPTL,ISPIN)
          END IF
          LAMBDA(:,:)=THIS%RLAM3M(:,:)
          CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,LAMBDA)
          IF(THISTASK.EQ.1) THEN
            WRITE(NFIL)LAMBDA !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            NREC=NREC+1
          END IF
        ENDDO
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_READLAMBDA(NFIL,IKPTG,TSKIP)
!     **************************************************************************
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                                **
!     **                                                                      **
!     **************************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NFIL
      INTEGER(4)   ,INTENT(IN) :: IKPTG
      LOGICAL(4)   ,INTENT(IN) :: TSKIP
      INTEGER(4)               :: NTASKS,THISTASK
      INTEGER(4)               :: ISPIN,I,IB1,IB2,IKPTL
      INTEGER(4)               :: NB_,NDIM_,NSPIN_
      INTEGER(4)               :: NB,NBA
      INTEGER(4)               :: NLAMBDA
      CHARACTER(8)             :: KEY
      COMPLEX(8)   ,ALLOCATABLE:: LAMBDA(:,:,:)
      COMPLEX(8)   ,ALLOCATABLE:: LAMBDA2(:,:,:)
      LOGICAL(4)               :: TKGROUP
!     **************************************************************************
                           CALL TRACE$PUSH('WAVES_READLAMBDA')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     ==========================================================================
!     == READ AND BROADCAST HISTORY-DEPTH OF LAGRANGE PARAMETERS              ==
!     ==========================================================================
      IF(THISTASK.EQ.1) THEN
        READ(NFIL)KEY,NB_,NDIM_,NSPIN_,NLAMBDA  !<<<<<<<<<<<<<<<<<<<<<<<<<<
        IF(KEY.NE.'LAMBDA') THEN
          CALL ERROR$MSG('ID IS NOT "LAMBDA"')
          CALL ERROR$MSG('FILE IS CORRUPTED')
          CALL ERROR$STOP('WAVES_READLAMBDA')
        END IF
      END IF
!
      IF(TSKIP) THEN
        IF(THISTASK.EQ.1) THEN !SCROLL FORWARD ON TASK 1.
          DO I=1,NLAMBDA
            DO ISPIN=1,NSPIN_
              READ(NFIL)
            ENDDO
          ENDDO
        END IF
        CALL TRACE$POP
        RETURN
      END IF
!
!     ==========================================================================
!     == SET TKGROUP=TRUE FOR ALL TASKS FOR WHICH THIS K-POINT IS LOCAL       ==
!     == KMAP CONTAINS THE TASK ID OF THE MASTER TASK IN THE KGROUP           ==
!     == DISTRIBUTE WITHIN THE K-GROUP                                        ==
!     ==========================================================================
      TKGROUP=THISTASK.EQ.KMAP(IKPTG)
      CALL MPE$BROADCAST('K',1,TKGROUP)
!
!     == ONLY THE READING TASK AND THE RECEIVING TASKS NEED TO BE PRESENT ======
      IF(THISTASK.NE.1.AND.(.NOT.TKGROUP)) THEN
        CALL TRACE$POP
        RETURN
      END IF
!
!     ==========================================================================
!     ==  COMMUNICATE                                                         ==
!     ==========================================================================
!     == SEND DATE FROM THE FIRST TASK OF MONOMER-GROUP 
!     ==             TO THE FIRST TASAK OF THE RELEVANT K-GROUP
      CALL MPE$SENDRECEIVE('MONOMER',1,KMAP(IKPTG),NLAMBDA)
      CALL MPE$SENDRECEIVE('MONOMER',1,KMAP(IKPTG),NB_)
      CALL MPE$SENDRECEIVE('MONOMER',1,KMAP(IKPTG),NDIM_)
      CALL MPE$SENDRECEIVE('MONOMER',1,KMAP(IKPTG),NSPIN_)
      IF(TKGROUP) THEN
        CALL MPE$BROADCAST('K',1,NLAMBDA)
        CALL MPE$BROADCAST('K',1,NDIM_)
        CALL MPE$BROADCAST('K',1,NB_)
        CALL MPE$BROADCAST('K',1,NSPIN_)
!       __ DETERMINE LOCAL K-POINT INDEX FOR THIS K-POINT_______________________
        IKPTL=0
        DO I=1,IKPTG
          IF(KMAP(I).EQ.KMAP(IKPTG))IKPTL=IKPTL+1
        ENDDO
        CALL WAVES_SELECTWV(IKPTL,1)
        NB=THIS%NB
        NBA=MIN(NB,NB_)
      END IF
!
!     ==========================================================================
!     ==  READ AND FOLD DOWN                                                  ==
!     ==========================================================================
      IF(THISTASK.EQ.1.OR.THISTASK.EQ.KMAP(IKPTG)) THEN
        ALLOCATE(LAMBDA(NB_,NB_,NSPIN_))
      ELSE
        ALLOCATE(LAMBDA(1,1,1))
      END IF
      IF(TKGROUP) ALLOCATE(LAMBDA2(NB,NB,NSPIN))
      DO I=1,NLAMBDA
        IF(THISTASK.EQ.1) THEN
          DO ISPIN=1,NSPIN_
            READ(NFIL)LAMBDA(:,:,ISPIN)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          ENDDO
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',1,KMAP(IKPTG),LAMBDA)
!
!       ========================================================================
!       ==  TRANSFORM BETWEEN DATA MODELS (COPY INTO LAMBDA2)                 ==
!       ========================================================================
        IF(THISTASK.EQ.KMAP(IKPTG)) THEN  ! RESULT WILL BE COMMUNICATED LATER
          LAMBDA2=(0.D0,0.D0)
          IF(NSPIN.EQ.NSPIN_.AND.NDIM.EQ.NDIM_) THEN
            LAMBDA2(1:NBA,1:NBA,:)=LAMBDA(1:NBA,1:NBA,:)
          ELSE
!!$CALL ERROR$MSG('OPTION TEMPORARILY DISABLED')
!!$CALL ERROR$STOP('WAVES_READLAMBDA')
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
        END IF
!
!       ==============================================================
!       ==  MAP ONTO THIS                                           ==
!       ==============================================================
        IF(TKGROUP) THEN
          CALL MPE$BROADCAST('K',1,LAMBDA2)
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPTL,ISPIN)
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
        ENDIF
      ENDDO
      IF(ALLOCATED(LAMBDA))DEALLOCATE(LAMBDA)
      IF(ALLOCATED(LAMBDA))DEALLOCATE(LAMBDA2)
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
      INTEGER(4),ALLOCATABLE:: MAP3D(:,:,:)
      REAL(8)               :: GBASIN(3,3)
      INTEGER(4)            :: IVEC1(3)
      INTEGER(4),ALLOCATABLE:: IVECA(:,:)  !(3,NGA)
      REAL(8)               :: XK(3)  !KPOINT IN RELATIVE COORDINATES
      REAL(8)               :: XG(3)  !G-VECTOR IN RECIPROCAL COORDINATES
      INTEGER(4)            :: MING(3),MAXG(3)
      REAL(8)               :: SVAR
!     ******************************************************************
      CALL LIB$INVERTR8(3,GBASA,GBASIN)
!!$      IF(NGA.NE.NGB) THEN
!!$        CALL ERROR$MSG('INCONSISTENT NUMBER OF PLANE WAVES')
!!$        CALL ERROR$MSG('THIS ERROR STATEMENT IS NOT NECESSARY')
!!$        CALL ERROR$MSG('BECAUSE WAVES_MAPG CAN DEAL WITH NGA -NE NGB')
!!$        CALL ERROR$MSG('IT HAS BEEN INTRODUCED FOR TESTING PURPOSES')
!!$        CALL ERROR$I4VAL('NGA',NGA)
!!$        CALL ERROR$I4VAL('NGB',NGB)
!!$        CALL ERROR$STOP('WAVES_MAPG')
!!$      END IF
!
!     ==================================================================
!     == DETERMINE K-POINT                                            ==
!     ==================================================================
      XG(:)=MATMUL(GBASIN,GVECA(:,1))
      XK(:)=XG(:)-REAL(NINT(XG(:)),KIND=8) !KPOINT IN RELATIVE COORDINATES
PRINT*,'XK ',XK
!
!     ==================================================================
!     == DIMENSION AND ALLOCATE 3-D GRID                              ==
!     ==================================================================
      ALLOCATE(IVECA(3,NGA))
      MING=+10000000
      MAXG=-10000000
      DO IG=1,NGA
        XG(:)=MATMUL(GBASIN,GVECA(:,IG))
! CAUTION! WHAT IS THE RESULT OF NINT FOR A HALF INTEGER?
!          NINT=INT(A+0.5) FOR A>0 AND INT(A-0.5) FOR A<0 OR A=0
!          INT(A) STRIPS THE FRACTIONAL PART OF A
        IVECA(:,IG)=NINT(XG(:)-XK(:))
!       == CATCH PROBLEM WHEN MAPPING FAILS ====================================
        IF(SUM((XG(:)-XK(:)-REAL(IVECA(:,IG),KIND=8))**2).GT.1.D-6) THEN
           CALL ERROR$MSG('INTERNAL ERROR')
           CALL ERROR$MSG('G-VECTORS CANNOT BE REPRESENTED AS SUM OF')
           CALL ERROR$MSG('K-VECTOR AND RECIPROCAL LATTICE VECTORS')
           CALL ERROR$MSG('FAILED TO MAP G-VECTORS')
           CALL ERROR$I4VAL('NGA',NGA)
           CALL ERROR$I4VAL('IG',IG)
           CALL ERROR$R8VAL('XG(1)',XG(1))
           CALL ERROR$R8VAL('XG(2)',XG(2))
           CALL ERROR$R8VAL('XG(3)',XG(3))
           CALL ERROR$R8VAL('XK(1)',XK(1))
           CALL ERROR$R8VAL('XK(2)',XK(2))
           CALL ERROR$R8VAL('XK(3)',XK(3))
           CALL ERROR$STOP('WAVES_MAPG')
        END IF
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
          CALL ERROR$I4VAL('NGA',NGA)
          CALL ERROR$I4VAL('IG',IG)
          CALL ERROR$I4VAL('MAP3D',MAP3D(IVECA(1,IG),IVECA(2,IG),IVECA(3,IG)))
          CALL ERROR$I4VAL('IVECA(1,IG)',IVECA(1,IG))
          CALL ERROR$I4VAL('IVECA(2,IG)',IVECA(2,IG))
          CALL ERROR$I4VAL('IVECA(3,IG)',IVECA(3,IG))
          CALL ERROR$R8VAL('XK_1',XK(1))
          CALL ERROR$R8VAL('XK_2',XK(2))
          CALL ERROR$R8VAL('XK_3',XK(3))
          CALL ERROR$STOP('WAVES_MAPG')
        END IF
      ENDDO
      DEALLOCATE(IVECA)
!
!     ==================================================================
!     ==  FOR EACH VECTOR (B) PICK OUT THE CLOSEST GRID (A) POINT     ==
!     ==================================================================
      LOOP1:DO IG=1,NGB
        MAP(IG)=0
        XG(:)=MATMUL(GBASIN,GVECB(:,IG))-XK(:)
        DO I=1,3
          IVEC1(I)=NINT(XG(I))
          IF(IVEC1(I).LT.MING(I).OR.IVEC1(I).GT.MAXG(I)) CYCLE LOOP1
        ENDDO
        MAP(IG)=MAP3D(IVEC1(1),IVEC1(2),IVEC1(3))
      ENDDO LOOP1
      DEALLOCATE(MAP3D)
!
!     ==================================================================
!     ==  PRINT IF REMAPPING HAS BEEN DONE                            ==
!     ==================================================================
      DO IG=1,NGB
        IF(MAP(IG).NE.IG) THEN
          PRINT*,'ORDER OF G-VECTORS FROM RESTART FILE HAS CHANGED'
          EXIT
        END IF
      END DO
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
      USE MPE_MODULE
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
      LOGICAL               :: TKGROUP1,TKGROUP2
      INTEGER(4)            :: IKPTL1,IKPTL2,IKPT
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      TKGROUP1=(THISTASK.EQ.KMAP(IKPT1))
      CALL MPE$BROADCAST('K',1,TKGROUP1) 
      TKGROUP2=(THISTASK.EQ.KMAP(IKPT2))
      IF(.NOT.(TKGROUP1.OR.TKGROUP2)) RETURN
      IKPTL1=0
      DO IKPT=1,IKPT1
        IF(KMAP(IKPT).EQ.KMAP(IKPT1))IKPTL1=IKPTL1+1
      ENDDO
      IKPTL2=0
      DO IKPT=1,IKPT2
        IF(KMAP(IKPT).EQ.KMAP(IKPT2))IKPTL2=IKPTL2+1
      ENDDO
!
!     ==================================================================
!     ==  MAPPING FROM IKPT2 TO IKPT1                                 ==
!     ==================================================================
!     == GET GLOBAL GVECTORS FROM IKPT1 ================================
      IF(TKGROUP1) THEN
        CALL WAVES_SELECTWV(IKPTL1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NGL1=GSET%NGL
        ALLOCATE(GVECL(3,NGL1))
        CALL PLANEWAVE$GETR8A('GVEC',3*NGL1,GVECL)
        CALL PLANEWAVE$GETI4('NGG',NGG1)
        ALLOCATE(GVECG1(3,NGG1))
        CALL PLANEWAVE$COLLECTR8(3,NGL1,GVECL,NGG1,GVECG1)
        DEALLOCATE(GVECL)
      END IF
!
!     == GET GLOBAL GVECTORS FROM IKPT2 ================================
      IF(TKGROUP2) THEN
        CALL WAVES_SELECTWV(IKPTL2,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NGL2=GSET%NGL
        ALLOCATE(GVECL(3,NGL2))
        CALL PLANEWAVE$GETR8A('GVEC',3*NGL2,GVECL)
        CALL PLANEWAVE$GETI4('NGG',NGG2)
        ALLOCATE(GVECG2(3,NGG2))
        CALL PLANEWAVE$COLLECTR8(3,NGL2,GVECL,NGG2,GVECG2)
        DEALLOCATE(GVECL)
      END IF
!
!     == COPY GBAS      ================================================
      IF(THISTASK.EQ.KMAP(IKPT2)) THEN
        CALL WAVES_SELECTWV(IKPTL2,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETR8A('GBAS',9,GBAS)
      END IF
      CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),GBAS)
      IF(TKGROUP1) THEN
        CALL MPE$BROADCAST('K',1,GBAS)        
        CALL WAVES_SELECTWV(IKPTL1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$SETR8A('GBAS',9,GBAS)
      END IF
!
!     == DEFINE MAPPING ================================================
      CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),NGG2)
      IF(THISTASK.EQ.KMAP(IKPT1).AND..NOT.TKGROUP2) ALLOCATE(GVECG2(3,NGG2))
      CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),GVECG2)
      IF(THISTASK.EQ.KMAP(IKPT1)) THEN
        ALLOCATE(MAPG(NGG1))
        CALL WAVES_MAPG(NGG2,GBAS,GVECG2,NGG1,GVECG1,MAPG)
        DEALLOCATE(GVECG2)
      END IF
      IF(ALLOCATED(GVECG1)) DEALLOCATE(GVECG1)
      IF(ALLOCATED(GVECG2)) DEALLOCATE(GVECG2)
!
!     ==================================================================
!     ==  COPY                                                        ==
!     ==================================================================
      IF(TKGROUP1) THEN
        CALL WAVES_SELECTWV(IKPTL1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETL4('TINV',TSUPER1)
        NB1=THIS%NB
        NBH1=THIS%NBH
      END IF
!CAUTION: SHOULD THIS BE TKGROUP2???
      IF(TKGROUP1) THEN
        CALL WAVES_SELECTWV(IKPTL2,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETL4('TINV',TSUPER2)
        NB2=THIS%NB
        NBH2=THIS%NBH
      END IF
      CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),NB2)
      CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),NBH2)

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
          IF(TKGROUP2) THEN
!
!           ==============================================================
!           == MAP DATA ONTO TEMP ARRAY                                 ==
!           ==============================================================
            CALL WAVES_SELECTWV(IKPTL2,ISPIN)
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
!           ==============================================================
!           ==  EXPAND TO REGULAR WAVE FUNCTIONS                        ==
!           ==============================================================
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
          END IF
!
!         ==============================================================
!         ==  MAP ONTO IKPT1                                          ==
!         ==============================================================
          DO IB=1,NBN
            DO IDIM=1,NDIM
              IF(TKGROUP2) THEN
                CALL WAVES_SELECTWV(IKPTL2,ISPIN)
                CALL PLANEWAVE$SELECT(GSET%ID)
                CALL PLANEWAVE$COLLECTC8(1,NGL2,PSI2(1,IDIM,IB),NGG2,PSIG2)
              END IF
              CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),PSIG2)
              IF(THISTASK.EQ.KMAP(IKPT1)) THEN
                DO IG=1,NGG1
                  IF(MAPG(IG).NE.0) THEN
                    PSIG1(IG)=PSIG2(MAPG(IG))
                  ELSE
                    PSIG1(IG)=(0.D0,0.D0)
                  END IF
                ENDDO
              END IF
              IF(TKGROUP1) THEN
                CALL WAVES_SELECTWV(IKPTL1,ISPIN)
                CALL PLANEWAVE$SELECT(GSET%ID)
                CALL PLANEWAVE$DISTRIBUTEC8(1,NGG1,PSIG1,NGL1,PSI1(1,IDIM,IB))
              END IF
            ENDDO
          ENDDO
!
!         ==============================================================
!         == MAP ONTO SUPER WAVE FUNCTIONS                            ==
!         ==============================================================
          IF(TKGROUP1.AND.TSUPER1) THEN
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
          IF(TKGROUP1) THEN
            CALL WAVES_SELECTWV(IKPTL1,ISPIN)
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
          END IF
        ENDDO
      ENDDO
      IF(THISTASK.EQ.KMAP(IKPT1))DEALLOCATE(MAPG)
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
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: IKPT1
      INTEGER(4)   ,INTENT(IN) :: IKPT2
      INTEGER(4)               :: NB1,NB2,NBA
      INTEGER(4)               :: ISPIN
      COMPLEX(8)   ,ALLOCATABLE:: LAMBDA(:,:)
      LOGICAL(4)               :: TKGROUP1,TKGROUP2
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: IKPTL1,IKPTL2,IKPT
      INTEGER(4)               :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      TKGROUP1=(THISTASK.EQ.KMAP(IKPT1))
      CALL MPE$BROADCAST('K',1,TKGROUP1) 
      TKGROUP2=(THISTASK.EQ.KMAP(IKPT2))
      IF(.NOT.(TKGROUP1.OR.TKGROUP2)) RETURN
      IKPTL1=0
      DO IKPT=1,IKPT1
        IF(KMAP(IKPT).EQ.KMAP(IKPT1))IKPTL1=IKPTL1+1
      ENDDO
      IKPTL2=0
      DO IKPT=1,IKPT2
        IF(KMAP(IKPT).EQ.KMAP(IKPT2))IKPTL2=IKPTL2+1
      ENDDO
!
      IF(TKGROUP1) THEN
        CALL WAVES_SELECTWV(IKPTL1,1)
        NB1=THIS%NB
      END IF
      IF(TKGROUP2) THEN
        CALL WAVES_SELECTWV(IKPTL2,1)
        NB2=THIS%NB
      END IF
      CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),NB2)
      CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT1),KMAP(IKPT2),NB1)
      NBA=MIN(NB1,NB2)
      ALLOCATE(LAMBDA(NBA,NBA))
!
!     ====================================================================
!     == COPY RLAM(0)                                                   ==
!     ====================================================================
      DO ISPIN=1,NSPIN
!       -- CHECK IF A LAMBDA MATRIX IS PRESENT ----------
        IF(THISTASK.EQ.KMAP(IKPT2)) THEN
          CALL WAVES_SELECTWV(IKPTL2,ISPIN)
          TCHK=(.NOT.ASSOCIATED(THIS%RLAM0))
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),TCHK)
        CALL MPE$BROADCAST('K',1,TCHK)
        IF(TCHK) RETURN
!
!       -- NOW COPY
        IF(THISTASK.EQ.KMAP(IKPT2)) THEN
          LAMBDA=THIS%RLAM0(1:NBA,1:NBA)
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),LAMBDA)
        IF(THISTASK.EQ.KMAP(IKPT1)) THEN
          CALL WAVES_SELECTWV(IKPTL1,ISPIN)
          IF(.NOT.ASSOCIATED(THIS%RLAM0)) ALLOCATE(THIS%RLAM0(NB1,NB1))
          THIS%RLAM0=(0.D0,0.D0)        
          THIS%RLAM0(1:NBA,1:NBA)=LAMBDA
        END IF
        IF(TKGROUP1) THEN
          CALL MPE$BROADCAST('K',1,THIS%RLAM0)
        END IF
      ENDDO
!
!     ====================================================================
!     == COPY RLAM(-)                                                   ==
!     ====================================================================
      DO ISPIN=1,NSPIN
!       -- CHECK IF A LAMBDA MATRIX IS PRESENT ----------
        IF(THISTASK.EQ.KMAP(IKPT2)) THEN
          CALL WAVES_SELECTWV(IKPTL2,ISPIN)
          TCHK=(.NOT.ASSOCIATED(THIS%RLAMM))
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),TCHK)
        CALL MPE$BROADCAST('K',1,TCHK)
        IF(TCHK) RETURN
!
!       -- NOW COPY
        IF(THISTASK.EQ.KMAP(IKPT2)) THEN
          LAMBDA=THIS%RLAMM(1:NBA,1:NBA)
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),LAMBDA)
        IF(THISTASK.EQ.KMAP(IKPT1)) THEN
          CALL WAVES_SELECTWV(IKPTL1,ISPIN)
          IF(.NOT.ASSOCIATED(THIS%RLAMM)) ALLOCATE(THIS%RLAMM(NB1,NB1))
          THIS%RLAMM=(0.D0,0.D0)        
          THIS%RLAMM(1:NBA,1:NBA)=LAMBDA
        END IF
        IF(TKGROUP1) THEN
          CALL MPE$BROADCAST('K',1,THIS%RLAMM)
        END IF
      ENDDO
!
!     ====================================================================
!     == COPY RLAM(2-)                                                  ==
!     ====================================================================
      DO ISPIN=1,NSPIN
!       -- CHECK IF A LAMBDA MATRIX IS PRESENT ----------
        IF(THISTASK.EQ.KMAP(IKPT2)) THEN
          CALL WAVES_SELECTWV(IKPTL2,ISPIN)
          TCHK=(.NOT.ASSOCIATED(THIS%RLAM2M))
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),TCHK)
        CALL MPE$BROADCAST('K',1,TCHK)
        IF(TCHK) RETURN
!
!       -- NOW COPY
        IF(THISTASK.EQ.KMAP(IKPT2)) THEN
          LAMBDA=THIS%RLAM2M(1:NBA,1:NBA)
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),LAMBDA)
        IF(THISTASK.EQ.KMAP(IKPT1)) THEN
          CALL WAVES_SELECTWV(IKPTL1,ISPIN)
          IF(.NOT.ASSOCIATED(THIS%RLAM2M)) ALLOCATE(THIS%RLAM2M(NB1,NB1))
          THIS%RLAM2M=(0.D0,0.D0)        
          THIS%RLAM2M(1:NBA,1:NBA)=LAMBDA
        END IF
        IF(TKGROUP1) THEN
          CALL MPE$BROADCAST('K',1,THIS%RLAM2M)
        END IF
      ENDDO
!
!     ====================================================================
!     == COPY RLAM(3-)                                                  ==
!     ====================================================================
      DO ISPIN=1,NSPIN
!       -- CHECK IF A LAMBDA MATRIX IS PRESENT ----------
        IF(THISTASK.EQ.KMAP(IKPT2)) THEN
          CALL WAVES_SELECTWV(IKPTL2,ISPIN)
          TCHK=(.NOT.ASSOCIATED(THIS%RLAM3M))
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),TCHK)
        CALL MPE$BROADCAST('K',1,TCHK)
        IF(TCHK) RETURN
!
!       -- NOW COPY
        IF(THISTASK.EQ.KMAP(IKPT2)) THEN
          LAMBDA=THIS%RLAM3M(1:NBA,1:NBA)
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPT2),KMAP(IKPT1),LAMBDA)
        IF(THISTASK.EQ.KMAP(IKPT1)) THEN
          CALL WAVES_SELECTWV(IKPTL1,ISPIN)
          IF(.NOT.ASSOCIATED(THIS%RLAM3M)) ALLOCATE(THIS%RLAM3M(NB1,NB1))
          THIS%RLAM3M=(0.D0,0.D0)        
          THIS%RLAM3M(1:NBA,1:NBA)=LAMBDA
        END IF
        IF(TKGROUP1) THEN
          CALL MPE$BROADCAST('K',1,THIS%RLAM3M)
        END IF
      ENDDO
      DEALLOCATE(LAMBDA)
      RETURN
      END SUBROUTINE WAVES_COPYLAMBDA
!
!     ....................................................................
      SUBROUTINE WAVES_SPINOROVERLAP(NBH,NB,IKPT,QMAT)
!     ********************************************************************
!     **  CALCULATES THE OVERLAP BETWEEN THE SPIN CONTRIBUTION          **
!     **  OF THE SPINORS NEEDED TO EVALUATE THE TOTAL SPIN.             **
!     **  CALCULATES Q_K,I,L,J=1/2*<PSI_I,K|PSI_J,L>                    **
!     **  I AND J ARE BAND-INDICES                                      **
!     **  K AND L ARE 1 OR 2 AND ARE THE SPINOR PART OF THE WAVEFUNCTION *
!     ********************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NBH,NB,IKPT
      COMPLEX(8),INTENT(OUT) :: QMAT(2*NB*NSPIN,2*NB*NSPIN)
      COMPLEX(8),ALLOCATABLE :: AUXMAT(:,:)
      COMPLEX(8),ALLOCATABLE :: AUXMAT2(:,:)
      COMPLEX(8),ALLOCATABLE :: OPROJ(:,:,:)
      COMPLEX(8),ALLOCATABLE :: PSI0(:,:,:)
      INTEGER(4)             :: NGL,NBD,I,J,NDIMHALF,ISPIN
!     ******************************************************************
      IF(NDIM.EQ.1.AND.NSPIN.EQ.1) THEN
        CALL ERROR$MSG('S^2 ONLY POSSIBLE FOR NOT SPIN RESTRICTED CALCULATION')
        CALL ERROR$I4VAL('NDIM',NDIM)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('WAVES_SPINOROVERLAP')
      END IF
!
!     ==================================================================
!     ==  NONCOLLINEAR CASE                                           ==
!     ==================================================================
      IF(NDIM.EQ.2) THEN   !NON-COLLINEAR
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)      
        NGL=GSET%NGL   
        NBD=2*NB
        NDIMHALF=1
        ALLOCATE(AUXMAT(NBD,NBD))
        CALL WAVES_1COVERLAP(MAP,NDIMHALF,NBD,NBD,MAP%NPRO,THIS%PROJ,THIS%PROJ,AUXMAT)
        CALL WAVES_OVERLAP(.TRUE.,NGL,NDIMHALF,NBD,NBD,THIS%PSI0,THIS%PSI0,QMAT)
        DO J=1,NBD
          DO I=1,NBD
            QMAT(I,J)=(QMAT(I,J)+AUXMAT(I,J))*0.5D0
          ENDDO
        ENDDO
        DEALLOCATE(AUXMAT)
!
!     ==================================================================
!     ==  COLLINEAR SPIN-POLARIZED CASE                                           ==
!     ==================================================================
      ELSE  
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
        CALL WAVES_1COVERLAP(MAP,1,NBH*2,NB*2,MAP%NPRO,OPROJ,OPROJ,AUXMAT)
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
!.......................................................................
MODULE WAVESFIXRHO_MODULE
! THIS MODULE IS USED TO KEEP THE DENSITY FIXED
  REAL(8)   ,ALLOCATABLE    :: QLM(:,:)
  REAL(8)   ,ALLOCATABLE    :: RHO(:,:)
  COMPLEX(8),ALLOCATABLE    :: DENMAT(:,:,:,:)
END MODULE WAVESFIXRHO_MODULE
!
!      ...................................................................
       SUBROUTINE WAVES_FIXRHOGET(NRL,NDIMD,LMRXX,NAT,QLM_,RHO_,DENMAT_)
       USE WAVESFIXRHO_MODULE 
       IMPLICIT NONE
       INTEGER(4) ,INTENT(IN)  :: NRL,NDIMD,LMRXX,NAT
       REAL(8)    ,INTENT(OUT) :: QLM_(LMRXX,NAT),RHO_(NRL,NDIMD)
       COMPLEX(8) ,INTENT(OUT) :: DENMAT_(LMRXX,LMRXX,NDIMD,NAT)
!     ********************************************************************

       IF (.NOT.ALLOCATED(QLM)) THEN
          ALLOCATE(QLM(LMRXX,NAT))
          ALLOCATE(RHO(NRL,NDIMD))
          ALLOCATE(DENMAT(LMRXX,LMRXX,NDIMD,NAT))
          CALL WAVES$FIXRHOREAD()
       END IF
       QLM_(:,:)=QLM(:,:)
       RHO_(:,:)=RHO(:,:)
       DENMAT_(:,:,:,:)=DENMAT(:,:,:,:)
       END
!
!      .................................................................
       SUBROUTINE WAVES_FIXRHOSET(NRL,NDIMD,LMRXX,NAT,QLM_,RHO_,DENMAT_)
       USE WAVESFIXRHO_MODULE 
       IMPLICIT NONE
       INTEGER(4) ,INTENT(IN)  :: NRL,NDIMD,LMRXX,NAT
       REAL(8)    ,INTENT(IN)  :: QLM_(LMRXX,NAT),RHO_(NRL,NDIMD)
       COMPLEX(8) ,INTENT(IN)  :: DENMAT_(LMRXX,LMRXX,NDIMD,NAT)
!      ******************************************************************
       IF (.NOT.ALLOCATED(QLM)) THEN
          ALLOCATE(QLM(LMRXX,NAT))
          ALLOCATE(RHO(NRL,NDIMD))
          ALLOCATE(DENMAT(LMRXX,LMRXX,NDIMD,NAT))
       END IF
       QLM(:,:)=QLM_(:,:)
       RHO(:,:)=RHO_(:,:)
       DENMAT(:,:,:,:)=DENMAT_(:,:,:,:)
       END
!
!      .................................................................
       SUBROUTINE WAVES$FIXRHOREAD()
       USE WAVESFIXRHO_MODULE
       IMPLICIT NONE
       INTEGER(4)                 :: NFIL
!      ******************************************************************
       CALL FILEHANDLER$SETFILE('FIXRHO',.FALSE.,'FIXRHO.BIN')
       CALL FILEHANDLER$SETSPECIFICATION('FIXRHO','FORM','UNFORMATTED')
       CALL FILEHANDLER$SETSPECIFICATION('FIXRHO','POSITION','REWIND')
       CALL FILEHANDLER$UNIT('FIXRHO',NFIL)
       READ(NFIL)QLM(:,:)
       READ(NFIL)RHO(:,:)
       READ(NFIL)DENMAT(:,:,:,:)
       CALL FILEHANDLER$CLOSE('FIXRHO')
       RETURN
       END
!
!      .................................................................
       SUBROUTINE WAVES$FIXRHOWRITE()
       USE WAVESFIXRHO_MODULE
       IMPLICIT NONE
       INTEGER(4)                 :: NFIL
!      ******************************************************************
       CALL FILEHANDLER$SETFILE('FIXRHO',.FALSE.,'FIXRHO.BIN')
       CALL FILEHANDLER$SETSPECIFICATION('FIXRHO','FORM','UNFORMATTED')
       CALL FILEHANDLER$SETSPECIFICATION('FIXRHO','POSITION','REWIND')
       CALL FILEHANDLER$UNIT('FIXRHO',NFIL)
       WRITE(NFIL)QLM(:,:)
       WRITE(NFIL)RHO(:,:)
       WRITE(NFIL)DENMAT(:,:,:,:)
       CALL FILEHANDLER$CLOSE('FIXRHO')
       RETURN
       END
!
!.....................................................................
MODULE TOTALSPIN_MODULE
!** TOTALSPIN CONTAINS <S^2>,<S_X>,<S_Y>,<S_X>                      **
REAL(8)               :: TOTSPIN(4) ! DIFFERS FROM JOHANNES VERSION
END MODULE TOTALSPIN_MODULE
!
!     ..................................................................
      SUBROUTINE WAVES_TOTALSPINRESET()
      USE TOTALSPIN_MODULE, ONLY: TOTSPIN ! DIFFERS FROM JOHANNES VERSION
      TOTSPIN=0.D0
      END SUBROUTINE WAVES_TOTALSPINRESET
!
!     ..................................................................
      SUBROUTINE WAVES_TOTALSPIN(NBX,NB,NKPT,IKPT,NSPIN,OCC,QMAT)
!     ******************************************************************
!     **  CALCULATES <S^2>,<S_X>,<S_Y>,<S_X> IN UNITS OF HBAR^2       **
!     **                                    1                          **
!     **                                                              **
!     ******************************************************************
      USE TOTALSPIN_MODULE, ONLY: TOTSPIN ! DIFFERS FROM JOHANNES VERSION
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NBX,NB,NKPT,IKPT,NSPIN
      REAL(8)   ,INTENT(IN) :: OCC(NBX,NKPT,NSPIN) 
      COMPLEX(8),INTENT(IN) :: QMAT(2,NB*NSPIN,2,NB*NSPIN) ! Q_I,J,K,L=1/2*<PSI_I,K|PSI_J,K>  
      COMPLEX(8)            :: CSUM,PART,SPINX,SPINY,SPINZ
      COMPLEX(8)            :: EXPECTS2
      REAL(8)               :: IOCC(NB*NSPIN)
      INTEGER(4)            :: I,J,NTASKS,THISTASK,NBD
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)

      DO I=1,NB
         IOCC(I)=OCC(I,IKPT,1)
      END DO
      IF(NSPIN.EQ.2) THEN
        DO I=1,NB
          IOCC(NB+I)=OCC(I,IKPT,2)
        END DO
      END IF
      NBD=NSPIN*NB !ONE BAND FOR ONE ELECTRON
!
!     ==================================================================
!     == SPIN EXPECTATION VALUE OF THE SLATER DETERMINANT             ==
!     ==================================================================
      SPINX=(0.D0,0.D0)
      SPINY=(0.D0,0.D0)
      SPINZ=(0.D0,0.D0)
      DO I=1,NBD
        SPINX=SPINX+(QMAT(1,I,2,I)+QMAT(2,I,1,I))*IOCC(I)
        SPINY=SPINY-CI*(QMAT(1,I,2,I)-QMAT(2,I,1,I))*IOCC(I)
        SPINZ=SPINZ+(QMAT(1,I,1,I)-QMAT(2,I,2,I))*IOCC(I)
      END DO
!
!     ==================================================================
!     == EXPECTATION VALUE OF S**2 OF THE SLATER DETERMINANT          ==
!     ==================================================================
      EXPECTS2=0.75D0*CMPLX(SUM(IOCC),0.D0,KIND=8)
!
!     ==================================================================
!     == NOW ADD THE EXCHANGE-LIKE CONTRIBUTION TO THE SPIN           ==
!     ==================================================================
      CSUM=(0.D0,0.D0)
      DO I=1,NBD
        DO J=1,NBD
          PART=(QMAT(1,I,1,J)-QMAT(2,I,2,J))*(QMAT(1,J,1,I)-QMAT(2,J,2,I)) &
     &        +4.D0*QMAT(1,I,2,J)*QMAT(2,J,1,I)
          CSUM=CSUM-PART*CMPLX(SQRT(ABS(IOCC(I)*IOCC(J))),0.D0,KIND=8)
        END DO
      END DO
      PRINT*,'TOTALSPIN: EXCHANGE LIKE CONTRIBUTION : ',CSUM
      PRINT*,'TOTALSPIN: PART 3/4+X+Y+Z: ',EXPECTS2+SPINX**2+SPINY**2+SPINZ**2
      CSUM=CSUM+EXPECTS2+SPINX**2+SPINY**2+SPINZ**2
!
!     ==================================================================
!     == PRINT RESULT                                                 ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        PRINT*,'SPIN: 2*<S_X>: ',2.D0*REAL(SPINX,KIND=8),' HBAR'
        PRINT*,'SPIN: 2*<S_Y>: ',2.D0*REAL(SPINY,KIND=8),' HBAR'
        PRINT*,'SPIN: 2*<S_Z>: ',2.D0*REAL(SPINZ,KIND=8),' HBAR'
        PRINT*,'TOTALSPIN: <S^2>: ',REAL(CSUM,KIND=8),' HBAR^2'
        PRINT*,'SPIN QUANTUM NUMBER S=',-0.5+SQRT(0.25+REAL(CSUM,KIND=8))
        IF (IKPT.EQ.1) TOTSPIN=0.D0 ! SUM UP TOTAL SPIN OVER K-POINTS
        TOTSPIN(1)=TOTSPIN(1)+REAL(CSUM,KIND=8)
        TOTSPIN(2)=TOTSPIN(2)+REAL(SPINX,KIND=8)
        TOTSPIN(3)=TOTSPIN(3)+REAL(SPINY,KIND=8)
        TOTSPIN(4)=TOTSPIN(4)+REAL(SPINZ,KIND=8)
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
      USE TOTALSPIN_MODULE, ONLY: TOTSPIN ! DIFFERS FROM JOHANNES VERSION
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
         SVAR=SQRT(TOTSPIN(2)**2+TOTSPIN(3)**2+TOTSPIN(4)**2)
         IF (NDIM.EQ.2) CALL REPORT$R8VAL(NFIL,'|S|',SVAR,'HBAR')
!        CALL REPORT$R8VAL(NFIL,'TOTAL SPIN <S^2>',TOTSPIN(1),'HBAR^2')
!        SVAR=-0.5+SQRT(0.25+TOTSPIN(1))
!        CALL REPORT$R8VAL(NFIL,'SPIN QUANTUM NUMBER S',SVAR,' ')
         RETURN
      END IF
    END SUBROUTINE WAVES$REPORTSPIN
!
!     ....................................................................
      SUBROUTINE WAVES_COMPAREPSI(ID,NGL,NDIM,NBH,PSI1,PSI2)
      USE WAVES_MODULE, ONLY : GSET,DELT
      IMPLICIT NONE
      CHARACTER(*)    ,INTENT(IN) :: ID
      INTEGER(4)      ,INTENT(IN) :: NGL
      INTEGER(4)      ,INTENT(IN) :: NDIM
      INTEGER(4)      ,INTENT(IN) :: NBH
      COMPLEX(8)      ,INTENT(IN) :: PSI1(NGL,NDIM,NBH)
      COMPLEX(8)      ,INTENT(IN) :: PSI2(NGL,NDIM,NBH)
      COMPLEX(8)                  :: CSVAR
      INTEGER(4)                  :: IG,IDIM,IB
      REAL(8)                     :: RSUM,SUMTOT
      REAL(8)                     :: RBAS(3,3),GBAS(3,3),CELLVOL
!     *********************************************************************
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL GBASS(RBAS,GBAS,CELLVOL)
      SUMTOT=0.D0
      DO IB=1,NBH
        RSUM=0.D0
        DO IG=1,NGL
          DO IDIM=1,NDIM
            CSVAR=PSI1(IG,IDIM,IB)-PSI2(IG,IDIM,IB)
            RSUM=RSUM+GSET%MPSI(IG)*(REAL(CSVAR,KIND=8)**2+AIMAG(CSVAR)**2)
          ENDDO
        ENDDO
        SUMTOT=SUMTOT+RSUM
        WRITE(*,*)' KINETICENERGYTEST ',IB,RSUM*CELLVOL/DELT**2,ID
      ENDDO
        WRITE(*,*)' KINETICENERGYTEST TOTAL',SUMTOT*CELLVOL/DELT**2,ID
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_READPSI(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **  READ WAVE FUNCTIONS FROM RESTART FILE                               **
!     **                                                                      **
!     **************************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL 
      INTEGER(4)              :: NTASKS,THISTASK
      INTEGER(4)              :: NKPT_,NSPIN_
      INTEGER(4)              :: NGG_,NDIM_,NB_,NBH_
      INTEGER(4)              :: NGG,NGL,NBH
      LOGICAL(4)              :: TSUPER_,TSUPER
      INTEGER(4)              :: IBH
      INTEGER(4)              :: IB,IDIM,IG,ISPIN,I
      INTEGER(4)              :: IKPTG,IKPTL
      INTEGER(4)              :: IWAVE
      REAL(8)                 :: K_(3) ! K-POINT ON FILE IN RELATIVE COORDINATES
      REAL(8)                 :: GBAS_(3,3)  ! REC. LATT. VECT. ON FILE
      REAL(8)     ,ALLOCATABLE:: GVECG_(:,:) ! G-VECTORS ON FILE
      REAL(8)     ,ALLOCATABLE:: GVECG(:,:)  ! G-VECTORS (GLOBAL)
      REAL(8)     ,ALLOCATABLE:: GVECL(:,:)  ! G-VECTORS (LOCAL)
      REAL(8)                 :: RBAS(3,3)
      REAL(8)                 :: SVAR
      INTEGER(4)  ,ALLOCATABLE:: MAPG(:)
      COMPLEX(8)  ,ALLOCATABLE:: PSITMP(:,:,:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSIL(:,:,:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSIG(:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSIIN(:,:)
      INTEGER(4)              :: IOS
      CHARACTER(8)            :: KEY
      INTEGER(4) ,ALLOCATABLE :: IGVECG_(:,:)
      REAL(8)                 :: XG(3)
      INTEGER(4)              :: NWAVE
      INTEGER(4)              :: NFILO
      LOGICAL(4)              :: TKGROUP
      REAL(8)    ,ALLOCATABLE :: XK(:,:) !K-POINTS IN RELATIVE COORDINATES
      REAL(8)                 :: GBASIN(3,3)
      INTEGER(4)              :: ILOGICAL
 REAL(8)     ,ALLOCATABLE:: TEST(:,:)
!     **************************************************************************
                               CALL TRACE$PUSH('WAVES_READPSI')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
      ALLOCATE(XK(3,NKPT))
      CALL DYNOCC$GETR8A('XK',3*NKPT,XK)
!
!     ==================================================================
!     ==  READ SIZES                                                  ==
!     ==   NWAVE=2 PSI AND PSIDOT; NWAVE=1: ONLY PSI                  ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        READ(NFIL,ERR=100)NKPT_,NSPIN_,RBAS,NWAVE
 100    CONTINUE
      END IF
      CALL MPE$BROADCAST('MONOMER',1,NSPIN_)
      CALL MPE$BROADCAST('MONOMER',1,NKPT_)
      CALL MPE$BROADCAST('MONOMER',1,RBAS)
      CALL MPE$BROADCAST('MONOMER',1,NWAVE)
      IF(NWAVE.EQ.1)TSTOP=.TRUE.  ! SET VELOCITY TO ZERO IF RESTART=STATIC
      IF(NKPT.NE.NKPT_) THEN
        CALL ERROR$MSG('KPOINT SETS ARE FIXED')
        CALL ERROR$I4VAL('NKPT ON FILE',NKPT_)
        CALL ERROR$I4VAL('NKPT EXPECTED',NKPT)
        CALL ERROR$I4VAL('NSPIN ON FILE',NSPIN_)
        CALL ERROR$I4VAL('NSPIN EXPECTED',NSPIN)
        CALL ERROR$STOP('WAVES_READPSI')
      END IF
!
!     ==========================================================================
!     ==  LOOP OVER K-POINTS AND SPINS                                        ==
!     ==========================================================================
      DO IKPTG=1,NKPT_   ! LOOP OVER ALL K-POINTS ON THE FILE
        TKGROUP=(THISTASK.EQ.KMAP(IKPTG))
        CALL MPE$BROADCAST('K',1,TKGROUP)
!       == NOW, TKGROUP=TRUE FOR ALL TASKS FOR WHICH THIS K-POINT IS LOCAL =====
!
!       ========================================================================
!       == COLLECT INFORMATION OF THE REQUIRED FORMAT                         ==
!       ========================================================================
        IF(TKGROUP) THEN
          IKPTL=0
          DO I=1,IKPTG
            IF(KMAP(I).EQ.KMAP(IKPTG)) IKPTL=IKPTL+1
          ENDDO            
          CALL WAVES_SELECTWV(IKPTL,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NBH=THIS%NBH
          CALL PLANEWAVE$GETI4('NGG',NGG)
          CALL PLANEWAVE$GETL4('TINV',TSUPER)
        ELSE
          IKPTL=0
          NGG=1
          NGL=0
          NBH=0
          TSUPER=.FALSE.
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,NGG)
!
!       ================================================================
!       ==  READ COORDINATES OF THE WAVE FUNCTIONS                    ==
!       ================================================================
        IF(THISTASK.EQ.1) THEN
!         == GFORTRAN LOGICAL REPRESENTATION DEFINED WITH TRUE=1, FALSE=0     ==
!         https://gcc.gnu.org/onlinedocs/gfortran/compiler-characteristics/
!         internal-representation-of-logical-variables.html
!         == IFORT LOGICAL REPRESENTATION DEFINED WITH VALUE OF LAST BIT      ==
!         https://www.intel.com/content/www/us/en/docs/fortran-compiler/
!         developer-guide-reference/2024-2/logical-data-representations.html
!         == BOTH SHARE MEANING OF LAST BIT 1=TRUE, 0=FALSE                   ==
!         == ENSURES BACKWARDS COMPATIBILITY WITH OLD RESTART FILES           ==
          READ(NFIL)KEY,NGG_,NDIM_,NB_,ILOGICAL   !<<<<<<<<<<<<<<<<<<<<<<
          TSUPER_=BTEST(ILOGICAL,0)
          IF(KEY.NE.'PSI') THEN
            CALL ERROR$MSG('ID IS NOT "PSI"')
            CALL ERROR$MSG('FILE IS CORRUPTED')
            CALL ERROR$I4VAL('IKPTG',IKPTG)
            CALL ERROR$STOP('WAVES_READPSI')
          END IF
!
          ALLOCATE(IGVECG_(3,NGG_))
          READ(NFIL)K_,IGVECG_ !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!$PRINT*,'READPSI (K ON FILE):  ',IKPTG,K_
!!$PRINT*,'READPSI (INTERNAL K): ',IKPTG,XK(:,IKPTG),GSET%ID
IF(SUM(ABS(K_(:)-XK(:,IKPTG))).GT.1.D-8) THEN
  CALL ERROR$MSG('INCONSISTENT K-POINTS')
  CALL ERROR$CHVAL('GSET%ID',GSET%ID)
  CALL ERROR$L4VAL('TKGROUP',TKGROUP)
  CALL ERROR$I4VAL('IKPTG',IKPTG)
  CALL ERROR$I4VAL('IKPTL',IKPTL)
  CALL ERROR$R8VAL('XK(1) ON FILE',K_(1))
  CALL ERROR$R8VAL('XK(2) ON FILE',K_(2))
  CALL ERROR$R8VAL('XK(3) ON FILE',K_(3))
  CALL ERROR$R8VAL('INTERNAL XK(1) ',XK(1,IKPTG))
  CALL ERROR$R8VAL('INTERNAL XK(2) ',XK(2,IKPTG))
  CALL ERROR$R8VAL('INTERNAL XK(3) ',XK(3,IKPTG))
  CALL ERROR$STOP('WAVES_READPSI')
END IF
!
          CALL GBASS(RBAS,GBAS_,SVAR)
          ALLOCATE(GVECG_(3,NGG_))
          DO IG=1,NGG_
            XG(:)=REAL(IGVECG_(:,IG),KIND=8)+K_(:)
            GVECG_(:,IG)=MATMUL(GBAS_,XG)
          ENDDO
          DEALLOCATE(IGVECG_)
!
! THE FOLLOWING ARE STATEMENTS THAT SEEM TO PREVENT A BUG.
! THE ERROR IS THAT GVECG_ CAN NO MORE BE REPRESENTED AS GBAS*(INTEGER+K)
! WHICH INDICATES THAT GVECG_ HAS CHANGED ITS VALUES. THERE IS A TEST IN
! MAPG, WHICH TESTS FOR THIS. AS SOON AS I TEST FOR IT IN READPSI, 
! THE ERROR DISAPPEARS.
!
ALLOCATE(TEST(3,NGG_))  ! WILL BE USED BEFORE WAVES_MAPG
TEST=GVECG_
IF(.FALSE.) THEN
  PRINT*,'-----------------TEST START'
  PRINT*,'K_',K_
  CALL LIB$INVERTR8(3,GBAS_,GBASIN)
  ALLOCATE(GVECG(3,NGG_))
  GVECG=MATMUL(GBASIN,GVECG_)
  DO I=1,3
    GVECG(I,:)=GVECG(I,:)-K_(I)
  ENDDO
  GVECG=GVECG-REAL(NINT(GVECG,KIND=8),KIND=8)  !!!!!
  DO I=1,NGG_
    IF(SUM(GVECG(:,I)**2).GT.1.D-6) THEN
      CALL ERROR$MSG('FAILURE')
      CALL ERROR$I4VAL('IG',IG)
      CALL ERROR$R8VAL('GVECG(1,1)',GVECG(1,1))
      CALL ERROR$R8VAL('GVECG(2,1)',GVECG(2,1))
      CALL ERROR$R8VAL('GVECG(3,1)',GVECG(3,1))
      CALL ERROR$R8VAL('GVECG(1,I)',GVECG(1,I))
      CALL ERROR$R8VAL('GVECG(2,I)',GVECG(2,I))
      CALL ERROR$R8VAL('GVECG(3,I)',GVECG(3,I))
      CALL ERROR$STOP('WAVES_READPSI')
    END IF
  ENDDO
  DEALLOCATE(GVECG)
  PRINT*,'-----------------TEST END'
END IF


        ELSE  !IF(THISTASK.NE.1) THEN
          KEY=' '
          NGG_=0
          NDIM_=0
          NB_=0
          TSUPER_=.FALSE.
          K_(:)=-999999999.D0
          GBAS_=0.D0
          SVAR=0.D0
          ! GVECG_ NOT ALLOCATED
          XG(:)=0.D0
        END IF
!
!       == NOW COMMUNICATE TO K-GROUP
        CALL MPE$SENDRECEIVE('MONOMER',1,KMAP(IKPTG),NDIM_)
        CALL MPE$SENDRECEIVE('MONOMER',1,KMAP(IKPTG),NB_)
        CALL MPE$SENDRECEIVE('MONOMER',1,KMAP(IKPTG),TSUPER_)
        IF(TKGROUP) THEN
          CALL MPE$BROADCAST('K',1,NDIM_)
          CALL MPE$BROADCAST('K',1,NB_)
          CALL MPE$BROADCAST('K',1,TSUPER_)
        END IF
        NBH_=NB_
        IF(TSUPER_)NBH_=(NB_+1)/2
!       
!       ==============================================================
!       ==  DEFINE G-VECTOR MAPPING OF THE ARRAYS                   ==
!       ==============================================================
!       == COLLECT G-VECTORS ON THE FIRST TASK OF K-GROUP
        ALLOCATE(GVECG(3,NGG))
        IF(TKGROUP) THEN
          ALLOCATE(GVECL(3,NGL))
          CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVECL)
          CALL PLANEWAVE$COLLECTR8(3,NGL,GVECL,NGG,GVECG)
          DEALLOCATE(GVECL)
        END IF
!       == SEND FROM TASK1 OF K-GROUP TO TASK1 OF MONOMER
        CALL MPE$SENDRECEIVE('MONOMER',KMAP(IKPTG),1,GVECG)
!
!       == CALCULATE MAPPING ARRAY FOR G-VECTOR MAPPING ==============
        IF(THISTASK.EQ.1) THEN
IF(SUM((TEST-GVECG_)**2).GT.1.D-6) THEN
  CALL ERROR$MSG('TEST AND GVECG_ DISAGREE:ERROR')
  CALL ERROR$STOP('WAVES_READPSI')
END IF
DEALLOCATE(TEST)
          ALLOCATE(MAPG(NGG))
          CALL WAVES_MAPG(NGG_,GBAS_,GVECG_,NGG,GVECG,MAPG)
          DEALLOCATE(GVECG_)
        END IF
        DEALLOCATE(GVECG)
!       
!       ==============================================================
!       ==  COLLECT DATA                                            ==
!       ==============================================================
        IF(THISTASK.EQ.1) ALLOCATE(PSIIN(NGG_,NDIM_))
        ALLOCATE(PSIG(NGG,NDIM_))
        IF(TKGROUP) ALLOCATE(PSIL(NGL,NDIM_,NBH_,NSPIN_))
!
        DO IWAVE=1,NWAVE  ! LOOP OVER PSI (AND PSIDOT )
!       
!         ==============================================================
!         ==  READ WAVE FUNCTIONS FROM FILE AND DISTRIBUTE TO KGROUP  ==
!         ==============================================================
          DO ISPIN=1,NSPIN_  ! SUM OVER SPIN COMPONENTS 
            DO IB=1,NBH_
              IF(THISTASK.EQ.1) THEN 
                READ(NFIL,ERR=9999,IOSTAT=IOS)PSIIN !<<<<<<<<<<<<<<<<<<<
!
!               == MAP G-VECTORS ONTO THEIR FINAL DESTINATION  ========
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
!
!             == SEND TO TKGROUP AND DISTRIBUTE
              CALL MPE$SENDRECEIVE('MONOMER',1,KMAP(IKPTG),PSIG)
!
              IF(TKGROUP) THEN
                DO IDIM=1,NDIM_
                  CALL PLANEWAVE$DISTRIBUTEC8(1,NGG,PSIG(1,IDIM),NGL,PSIL(1,IDIM,IB,ISPIN))
                ENDDO
              END IF
            ENDDO
          ENDDO
!       
!         ==============================================================
!         ==  CHANGE FORMAT                                           ==
!         ==============================================================
          IF(TKGROUP) THEN
!!$PRINT*,'TRANSFER(TSUPER,IBH)',TRANSFER(TSUPER,IBH)
!!$PRINT*,'TRANSFER(TSUPER_,IBH)',TRANSFER(TSUPER_,IBH)
!!$PRINT*,'TSUPER ',TSUPER,'TSUPER_ ',TSUPER_,'TSUPER.EQV.TSUPER_',TSUPER.EQV.TSUPER_,'TSUPER.NEQV.TSUPER_',TSUPER.NEQV.TSUPER_
            IF((TSUPER.EQV.TSUPER_).AND.(TSUPER.NEQV.TSUPER_)) THEN
              CALL ERROR$MSG('COMPILER BUG! SHUTTING DOWN..')
              CALL ERROR$MSG('.EQV. OR .NEQ. DO NOT TEST ONLY THE RELEVANT BIT')
              CALL ERROR$MSG('A.EQV.B AND A.NEQV.B CAN BE  BOTH TRUE')
              CALL ERROR$STOP('WAVES_READPSI')
            END IF

!            IF(TSUPER.NEQV.TSUPER_) THEN
            IF(.NOT.((TSUPER.AND.TSUPER_).OR.(.NOT.TSUPER.AND..NOT.TSUPER_))) THEN
!!$PRINT*,'(TSUPER.AND.TSUPER_)=',(TSUPER.AND.TSUPER_),'  (.NOT.TSUPER.AND..NOT.TSUPER_)=',(.NOT.TSUPER.AND..NOT.TSUPER_)
              CALL ERROR$MSG('TRANSFORMATION BETWEEN REGULAR ..')
              CALL ERROR$MSG('... AND SUPER WAVE FUNCTIONS NOT IMPLEMENTED')
              CALL ERROR$L4VAL('TSUPER_ ON FILE',TSUPER_)
              CALL ERROR$L4VAL('TSUPER  EXPECTED',TSUPER)
              CALL ERROR$I4VAL('IKPTG',IKPTG)
              CALL ERROR$I4VAL('IKPTL',IKPTL)
              CALL ERROR$I4VAL('IWAVE',IWAVE)
              CALL ERROR$STOP('WAVES_READPSI')
            END IF


            IF(NDIM_.NE.NDIM) THEN
              CALL ERROR$MSG('TRANSFORMATION BETWEEN SCALAR AND SPINOR ..')
              CALL ERROR$MSG('... WAVE FUNCTIONS NOT IMPLEMENTED')
              CALL ERROR$I4VAL('NDIM__ ON FILE',NDIM_)
              CALL ERROR$I4VAL('NDIM  EXPECTED',NDIM)
              CALL ERROR$STOP('WAVES_READPSI')
            END IF
            ALLOCATE(PSITMP(NGL,NDIM,NBH,NSPIN))
!           == COPY RANDOM WAVE FUNCTIONS INTO HOLDING SPACE
            DO ISPIN=1,NSPIN
              CALL WAVES_SELECTWV(IKPTL,ISPIN)
              CALL PLANEWAVE$SELECT(GSET%ID)
              IF(IWAVE.EQ.1) THEN
                PSITMP(:,:,:,ISPIN)=THIS%PSI0(:,:,:)
              ELSE
                PSITMP(:,:,:,ISPIN)=THIS%PSIM(:,:,:)
              END IF
            ENDDO
!
!           == COPY WAVE FUNCTIONS FROM THE FILE INTO THE HOLDING SPACE
            ISPIN=MIN(NSPIN,NSPIN_)
            IBH=MIN(NBH,NBH_)
            IDIM=MIN(NDIM,NDIM_)
            PSITMP(:,1:IDIM,1:IBH,1:ISPIN)=PSIL(:,1:IDIM,1:IBH,1:ISPIN)
            IF(NSPIN_.EQ.1.AND.NSPIN.EQ.2) THEN
              PSITMP(:,:,:,2)=PSITMP(:,:,:,1)
            END IF
!
!           == COPY HOLDING SPACE INTO THE MEMORY ===========================
            DO ISPIN=1,NSPIN
              CALL WAVES_SELECTWV(IKPTL,ISPIN)
              CALL PLANEWAVE$SELECT(GSET%ID)
              IF(IWAVE.EQ.1) THEN
                THIS%PSI0(:,:,:)=PSITMP(:,:,:,ISPIN)
              END IF
              IF(IWAVE.EQ.2.OR.NWAVE.EQ.1) THEN
                THIS%PSIM(:,:,:)=PSITMP(:,:,:,ISPIN)
              END IF
            ENDDO
            DEALLOCATE(PSITMP)
          END IF
        ENDDO  ! END LOOP OVER WAVES
        IF(THISTASK.EQ.1) DEALLOCATE(MAPG)
        IF(THISTASK.EQ.1) DEALLOCATE(PSIIN)
        DEALLOCATE(PSIG)
        IF(TKGROUP) DEALLOCATE(PSIL)
!
        CALL WAVES_READLAMBDA(NFIL,IKPTG,.FALSE.)
      ENDDO    ! END LOOK OVER K-POINTS
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
      END SUBROUTINE WAVES_READPSI
!!$!
!!$!     ...............................................................................
!!$      SUBROUTINE DUMPPSI(STRING)
!!$      USE WAVES_MODULE
!!$      IMPLICIT NONE
!!$      INTEGER(4)              :: ISPIN,IKPT
!!$      CHARACTER(*),INTENT(IN) :: STRING
!!$      INTEGER(4),PARAMETER :: NFIL=11
!!$!     **********************************************************************
!!$      OPEN(UNIT=NFIL,FILE=STRING,FORM='FORMATTED')
!!$      DO IKPT=1,NKPTL
!!$        DO ISPIN=1,NSPIN
!!$          CALL WAVES_SELECTWV(IKPT,ISPIN)
!!$          CALL PLANEWAVE$SELECT(GSET%ID)
!!$          WRITE(NFIL,*)THIS%PSI0
!!$          WRITE(NFIL,*)THIS%PSIM
!!$        ENDDO
!!$      ENDDO
!!$      CLOSE(NFIL)
!!$      RETURN
!!$      END




!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$TESTORTHO(ID)
!     **************************************************************************
!     **  DETERMINES THE OVERLAP MATRIX OF THE WAVE FUNCTIONS                 **
!     **************************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(1),INTENT(IN) :: ID
      REAL(8)   ,ALLOCATABLE :: R0(:,:)      !CURRENT ATOMIC POSITIONS
      COMPLEX(8),ALLOCATABLE :: MAT(:,:)
      COMPLEX(8),ALLOCATABLE :: AUXMAT(:,:)
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: IKPT,ISPIN
      INTEGER(4)             :: I,J
      INTEGER(4)             :: NGL,NBH,NB
      INTEGER(4)             :: NAT
      COMPLEX(8),ALLOCATABLE :: THISPROJ(:,:,:)
!     **************************************************************************
                             CALL TRACE$PUSH('WAVES$TESTORTHO')
      NPRO=MAP%NPRO
      NAT=MAP%NAT
!
!     ==========================================================================
!     == COLLECT OCCUPATIONS                                                  ==
!     ==========================================================================
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     ==  NOW ORTHOGONALIZE                                                   ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NBH=THIS%NBH
          NB=THIS%NB
!
!         ======================================================================
!         ==  CALCULATE PROJECTIONS FOR THE NEW POSITIONS                     ==
!         ======================================================================
          ALLOCATE(THISPROJ(NDIM,NBH,NPRO))
          IF(ID.EQ.'0') THEN
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R0,NGL,NDIM,NBH,NPRO &
     &                                                     ,THIS%PSI0,THISPROJ)
          ELSE IF(ID.EQ.'-') THEN
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R0,NGL,NDIM,NBH,NPRO &
     &                                                     ,THIS%PSIM,THISPROJ)
          ELSE
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('WAVES$TESTORTHO')
          END IF
          CALL MPE$COMBINE('K','+',THISPROJ)
!
!         ======================================================================
!         ==  1C-OVERLAP OF <PSI0|PSI0>, <OPSI|PSI0> AND <OPSI|OPSI>          ==
!         ======================================================================
          ALLOCATE(MAT(NB,NB))
          CALL WAVES_1COVERLAP(MAP,NDIM,NBH,NB,NPRO,THISPROJ,THISPROJ,MAT)
!
!         ======================================================================
!         ==  NOW ADD OVERLAP OF PSEUDO WAVE FUNCTIONS                        ==
!         ======================================================================
          ALLOCATE(AUXMAT(NB,NB))
          IF(ID.EQ.'0') THEN
            CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB,THIS%PSI0,THIS%PSI0,AUXMAT)
          ELSE IF(ID.EQ.'-') THEN
            CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB,THIS%PSIM,THIS%PSIM,AUXMAT)
          ELSE
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('WAVES$TESTORTHO')
          END IF
          DO I=1,NB
            DO J=1,NB
              MAT(I,J)=MAT(I,J)+AUXMAT(I,J)
            ENDDO
          ENDDO
WRITE(*,FMT='("%%%",80("="))')
DO J=1,NB
  WRITE(*,FMT='("%%%",20F10.5)')MAT(:,J)
ENDDO
WRITE(*,FMT='(80("-"))')
        ENDDO
      ENDDO
                                    CALL TRACE$POP
      RETURN
      END
