      MODULE ISOLATE_MODULE
!     ******************************************************************
!     **                                                              **
!     **  THE ISOLATE OBJECT                                          **
!     **  1) PROVIDES AN INTERFACE FOR COUPLING AN "EXTERNAL" SYSTEM  **
!     **     TO POSITIONS AND CHARGES                                 **
!     **  2) USES THIS INTERFACE TO SUBTRACT THE ELECTROSTATIC        **
!     **     INTERACTION OF PERIODIC IMAGES                           **
!     **                                                              **
!     **  REFERENCE:                                                  **
!     **    "ELECTROSTATIC DECOUPLING OF PERIODIC IMAGES OF           **
!     **     OF PLANE WAVE EXPANDED DENSITIES                         **
!     **     AND DERIVED POINT CHARGES";                          .   **
!     **     P.E. BLOECHL, J. CHEM.PHYS.103, P7422 (1995)             **
!     **                                                              **
!     ******************************************************************
      LOGICAL(4) :: TON=.FALSE.      ! ON/OFF SWITCH FOR ISOLATE OBJECT
      LOGICAL(4) :: TISOLATE=.FALSE. ! ON/OFF SWITCH FOR ISOLATION
      INTEGER(4) :: NFCT=3           ! #(GAUSSIANS) PER SITE
      REAL(8)    :: RC0=0.5D0        ! SMALLEST DECAY OF GAUSSIANS
      REAL(8)    :: RCFAC=1.5        ! FACTOR SEPARATING THE DECAY LENGTHS
      REAL(8)    :: G2MAX=3.D0       ! PLANE WAVE CUTOFF IN THE FIT
      END MODULE ISOLATE_MODULE
!
!     ..................................................................
      SUBROUTINE ISOLATE$ONOFF(STRING,NFCT_,RC0_,RCFAC_,G2MAX_,DECOUPLE_)
!     ******************************************************************
!     **  SWITCH ISOLATE ON AND OFF                                   **
!     ******************************************************************
      USE ISOLATE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: STRING
      INTEGER(4)  ,INTENT(IN) :: NFCT_
      REAL(8)     ,INTENT(IN) :: RC0_ 
      REAL(8)     ,INTENT(IN) :: RCFAC_
      REAL(8)     ,INTENT(IN) :: G2MAX_
      LOGICAL(4)  ,INTENT(IN) :: DECOUPLE_
!     ******************************************************************
      NFCT=NFCT_
      RC0=RC0_
      RCFAC=RCFAC_
      G2MAX=G2MAX_
      TISOLATE=DECOUPLE_
      IF(STRING.EQ.'ON') THEN
        TON=.TRUE.
!       __ TISOLATE DECIDES WHETHER ONLY CHARGES SHALL BE CALCULATED____
      ELSE IF(STRING.EQ.'OFF') THEN
        TON=.FALSE.
      ELSE
        CALL ERROR$MSG('STRING INCORRECT IN ISOLATE$ONOFF')
        CALL ERROR$CHVAL('STRING',STRING)
        CALL ERROR$STOP('ISOLATE$ONOFF')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ISOLATE$REPORT(NFIL)
!     ******************************************************************
!     **  REPORT SETING OF ISOLATE OBJECT                             **
!     ******************************************************************
      USE ISOLATE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL 
!     ******************************************************************
      IF(.NOT.TON) RETURN
      WRITE(NFIL,FMT='(/"ISOLATE OBJECT"/15("="))')
      IF(TON) THEN
        WRITE(NFIL,FMT='(55("."),": ",T1,"NUMBER OF GAUSSIANS"' &
     &           //',T58,I10)')NFCT
        WRITE(NFIL,FMT='(55("."),": ",T1,"SMALLEST GAUSSIAN CUTOFF"' &
     &           //',T58,F10.5," A.U.")')RC0
        WRITE(NFIL,FMT='(55("."),": ",T1,"FACTOR FOR GAUSSIAN SPACING"' &
     &           //',T58,F10.5)')RCFAC
        WRITE(NFIL,FMT='(55("."),": ",T1,"PLANE WAVE CUTOFF"' &
     &           //',T58,F10.5," RY")')G2MAX
        IF(TISOLATE) THEN
          WRITE(NFIL,FMT='("PERIODIC IMAGES SUBTRACTED")')
        ELSE 
          WRITE(NFIL,FMT='("PERIODIC IMAGES NOT SUBTRACTED")')
        END IF
      ELSE
        WRITE(NFIL,FMT='("ISOLATE OBJECT IS SWITCHED OFF ")')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ISOLATE(NSP,NAT,ISPECIES,POS,RBAS,FORCE,E &
     &                  ,LMRXX,LMRX,QLM,VQLM,RHOB,NG,RHO,POT)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: RHOB IS SET TO ZERO IF ISOLATION ON. RHOB HAS NO    **
!     **          OTHER FUNCTION IN THIS ROUTINE                      **
!     **                                                              **
!     ******************************************************************
      USE ISOLATE_MODULE
      IMPLICIT NONE
      LOGICAL    ,PARAMETER    :: TPR=.False.
      REAL(8)    ,PARAMETER    :: XK(3)=(/0.D0,0.D0,0.D0/)
      INTEGER(4) ,INTENT(IN)   :: NSP           ! #(ATOM TYPES)
      INTEGER(4) ,INTENT(IN)   :: NAT           ! #(ATOMS)
      INTEGER(4) ,INTENT(IN)   :: ISPECIES(NAT) ! ATOMS->ATOM TYPES
      REAL(8)    ,INTENT(IN)   :: POS(3,NAT)    ! ATOMIC POSITIONS
      REAL(8)    ,INTENT(IN)   :: RBAS(3,3)     ! UNIT CELL
      REAL(8)    ,INTENT(INOUT):: FORCE(3,NAT)  ! FORCE ON POSITIONS
      REAL(8)    ,INTENT(INOUT):: E             ! ENERGY
      INTEGER(4) ,INTENT(IN)   :: LMRXX         ! MAX(ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)   :: LMRX(NSP)
      REAL(8)    ,INTENT(IN)   :: QLM(LMRXX,NAT)
      REAL(8)    ,INTENT(INOUT):: VQLM(LMRXX,NAT)
      REAL(8)    ,INTENT(INOUT):: RHOB          ! COMPENSATING CHARGE BACKGROUND
      INTEGER(4) ,INTENT(IN)   :: NG            ! #(G-VECTORS)
      COMPLEX(8) ,INTENT(IN)   :: RHO(NG)
      COMPLEX(8) ,INTENT(INOUT):: POT(NG)
      REAL(8)    ,allocatable  :: G2(:)
      real(8)    ,allocatable  :: gvec(:,:)
      complex(8) ,allocatable  :: eigr(:)
      INTEGER(4)               :: NGS
      INTEGER(4) ,ALLOCATABLE  :: MAP(:)     !(NGS)       
      COMPLEX(8) ,ALLOCATABLE  :: RHOS(:)    !(NGS)   
      COMPLEX(8) ,ALLOCATABLE  :: F(:,:,:)   !(NGS,NFCT,NAT) 
      REAL(8)    ,ALLOCATABLE  :: GvecS(:,:) !(3,NGS)
      REAL(8)                  :: RC(NFCT,NAT)   
      REAL(8)                  :: RCSM(NAT)         
      real(8)                  :: decay(nfct)
      REAL(8)                  :: GBAS(3,3)
      REAL(8)                  :: OMEGA
      REAL(8)                  :: FORCET(3,NAT)     
      REAL(8)                  :: VQLMT(LMRXX,NAT)  
      REAL(8)                  :: ET
      INTEGER(4)               :: IAT,ISP,IFCT,IG,I,LM,igs
      INTEGER(4)               :: ngamma
!     ******************************************************************
      IF(.NOT.TON)RETURN
                              CALL TRACE$PUSH('ISOLATE')
!
!     == RMAX DETERMINES THE MAXIMUM PLANE WAVE CONSIDERED
!     == RC DETERMINES THE DECAYB OF THE GAUSSIAN
      IF(TPR) THEN
        WRITE(*,FMT='("NF=",I2," RC0=",F5.3," RCFAC=",F5.3' &
     &              //'," G2MAX ",F10.3)')NFCT,RC0,RCFAC,G2MAX
      END IF
      CALL GBASS(RBAS,GBAS,OMEGA)
      E=0.D0
!
!     ==================================================================
!     == SET CUTOFF FOR FIT FUNCTIONS                                 ==
!     ==================================================================
      DO IAT=1,NAT
        RC(1,IAT)=RC0
        DO IFCT=2,NFCT
          RC(IFCT,IAT)=RC(IFCT-1,IAT)*RCFAC
        ENDDO
      ENDDO
!      
!     == GET DECAY CONSTANT FOR THE GAUSSION OF THE COMPENSATION DENSITY
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        CALL SETUP$RCSM(ISP,RCSM(IAT))
      ENDDO
!
!     ==================================================================
!     == OBTAIN MAPPING TO LOWER PLANE WAVES                          ==
!     ==================================================================
      ALLOCATE(G2(NG))
      CALL PLANEWAVE$GETR8A('G2',NG,G2)
      CALL ISOLATE_LOWERNUMBER(G2MAX,NG,G2,NGS)
      ALLOCATE(MAP(NGS))
      CALL ISOLATE_MAPLOWER(G2MAX,NG,G2,NGS,MAP)
!
!     ==================================================================
!     == CALCULATE FIT FUNCTIONS                                      ==
!     ==================================================================
      ALLOCATE(F(NGS,NFCT,NAT))
      ALLOCATE(EIGR(NG))
      DO IAT=1,NAT
        CALL PLANEWAVE$STRUCTUREFACTOR(POS(1,IAT),NG,EIGR)
        DO IFCT=1,NFCT
          DECAY(IFCT)=-(0.5D0*RC(IFCT,IAT))**2
        ENDDO
        DO IGS=1,NGS
          IG=MAP(IGS)
          DO IFCT=1,NFCT
            F(IGS,IFCT,IAT)=EIGR(IG)*DEXP(DECAY(IFCT)*G2(IG))/OMEGA
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(EIGR)
      DEALLOCATE(G2)
!
!     ==================================================================
!     == MAP ONTO LOWER PLANE WAVES                                   ==
!     ==================================================================
      ALLOCATE(GVEC(3,NG))
      CALL PLANEWAVE$GETR8A('GVEC',3*NG,GVEC)
      CALL PLANEWAVE$GETI4('NGAMMA',NGAMMA)
      ALLOCATE(RHOS(NGS))
      ALLOCATE(GVECS(3,NGS))
      DO IGS=1,NGS
        IG=MAP(IGS)
        RHOS(IGS)=RHO(IG)
        GVECS(:,IGS)=GVEC(:,IG)
        if(ig.eq.ngamma) ngamma=igs
      ENDDO
      DEALLOCATE(GVEC)
!
!     ==================================================================
!     == CONSTRUCT CORRECTIONONS                                      ==
!     ==================================================================
      CALL ISOLATE_A(NGS,NAT,RBAS,ET,POS,FORCET,GvecS &
     &              ,ngamma,RHOS,LMRXX,QLM,VQLMT,NFCT,F,RC,RCSM)
      DEALLOCATE(F)
      DEALLOCATE(GvecS)

!     ==================================================================
!     ==  ADD CORRECTION TO POTENTIAL                                 ==
!     ==================================================================
      IF(TISOLATE) THEN
        RHOB=0.D0
        E=ET
        DO IAT=1,NAT
          DO I=1,3
            FORCE(I,IAT)=FORCE(I,IAT)+FORCET(I,IAT)
          ENDDO
          DO LM=1,LMRX(ISPECIES(IAT))
            VQLM(LM,IAT)=VQLM(LM,IAT)+VQLMT(LM,IAT)
          ENDDO
        ENDDO
        DO IG=1,NGS
          POT(MAP(IG))=POT(MAP(IG))+RHOS(IG)
        ENDDO
      END IF
      DEALLOCATE(MAP)
      DEALLOCATE(RHOS)
                              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ISOLATE_A(NGS,NAT,RBAS,E,BAS,FORCE,G &
     &                    ,ngamma,RHO,LMRXX,QLM,VQLM,NFCT,F,RC,RCSM)
!     ******************************************************************
!     ******************************************************************
      USE MPE_MODULE
      implicit none
      LOGICAL(4),PARAMETER    :: TPR=.FALSE.
      INTEGER(4),INTENT(IN)   :: NGS
      INTEGER(4),INTENT(IN)   :: NAT
      REAL(8)   ,INTENT(IN)   :: RBAS(3,3)
      REAL(8)   ,INTENT(OUT)  :: E
      REAL(8)   ,INTENT(IN)   :: BAS(3,NAT)
      REAL(8)   ,INTENT(OUT)  :: FORCE(3,NAT)
      REAL(8)   ,INTENT(IN)   :: G(3,NGS)
      integer(4),INTENT(IN)   :: ngamma
      COMPLEX(8),INTENT(INOUT):: RHO(NGS)
      INTEGER(4),INTENT(IN)   :: LMRXX
      REAL(8)   ,INTENT(IN)   :: QLM(LMRXX,NAT)
      REAL(8)   ,INTENT(OUT)  :: VQLM(LMRXX,NAT)
      INTEGER(4),INTENT(IN)   :: NFCT
      COMPLEX(8),INTENT(IN)   :: F(NGS,NFCT*NAT)
      REAL(8)   ,INTENT(IN)   :: RC(NFCT*NAT)
      REAL(8)   ,INTENT(IN)   :: RCSM(NAT)
      REAL(8)                 :: GBAS(3,3)
      REAL(8)                 :: RHOGAMMA
      REAL(8)                 :: EBACKGROUND
      REAL(8)                 :: VOL         ! CELL SIZE VOLUME   
      REAL(8)                 :: PI,Y0       ! PI AND S-SPHERICAL HARMONICS
      INTEGER(4)              :: I,J              ! LOOP INDICES
      INTEGER(4)              :: IAT,IAT1,IAT2    ! LOOP INDICES
      INTEGER(4)              :: IFCT,IFCT1,IFCT2 ! LOOP INDICES
      INTEGER(4)              :: IG,I1,I2,I3      ! LOOP INDICES
      INTEGER(4)              :: ICOUNT      ! INDEX FOR LOOP PARALLELIZATION
      INTEGER(4)              :: NTASKS      ! #(PARALLEL TASKS)
      INTEGER(4)              :: THISTASK    ! NR. OF THIS TASK
      REAL(8)                 :: SVAR,SUM    ! AUXILIARY VARIABLE
      REAL(8)                 :: FAC         ! AUXILIARY VARIABLE
      COMPLEX(8)              :: CSVAR       ! AUXILIARY VARIABLE
      INTEGER(4)              :: NFILO       ! PROTOCOLL FILE NUMBER
      REAL(8)   ,ALLOCATABLE  :: WORK1(:)    ! WORKSPACE FOR INVERSION
      REAL(8)                 :: XQI(NFCT+1,NAT)
      REAL(8)                 :: XRC(NFCT+1,NAT)
      REAL(8)                 :: XVI(NFCT+1,NAT)
      REAL(8)     :: GGAMMA(NFCT*NAT)
      REAL(8)     :: GRHO(NFCT*NAT)
      REAL(8)     :: GRADGRHO(3,NFCT*NAT)
      REAL(8)     :: QMAD(NAT)
      REAL(8)     :: VMAD(NAT)
      REAL(8)     :: GGINV(NFCT*NAT,NFCT*NAT)
      REAL(8)     :: GG(NFCT*NAT,NFCT*NAT)
      REAL(8)     :: GRADGG(3,NFCT*NAT,NFCT*NAT)
      REAL(8)     :: QI(NFCT*NAT)
      REAL(8)     :: VI(NFCT*NAT)
      REAL(8)     :: GRADQI(3,NFCT*NAT,NFCT*NAT)
      REAL(8)     :: XI(NFCT*NAT)
      REAL(8)     :: D(NFCT*NAT)
      REAL(8)     :: TDIP,DIP(3)
      real(8)     :: debye
      real(8)     :: g2,w
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
      CALL GBASS(RBAS,GBAS,VOL)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL MPE$QUERY('~',NTASKS,THISTASK)
!
!     ==================================================================
!     == INITIALIZE                                                   ==
!     ==================================================================
      DO IAT=1,NAT
        DO I=1,3
          FORCE(I,IAT)=0.D0
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  DEFINE GGAMMA,RHOGAMMA                                      ==
!     ==================================================================
      IF(NGAMMA.NE.0) THEN
        RHOGAMMA=REAL(RHO(NGAMMA))
        DO I=1,NAT*NFCT
           GGAMMA(I)=REAL(F(NGAMMA,I))
        ENDDO
      ELSE
        RHOGAMMA=0.D0
        DO I=1,NAT*NFCT
          GGAMMA(I)=0.D0
        ENDdo
      END IF
      CALL MPE$COMBINE('+',GGAMMA)
      CALL MPE$COMBINE('+',RHOGAMMA)
!
!     ==================================================================
!     ==      GRHO(I)=      <F(I)|RHO>      GG(I,J)=      <F(I)|F(J)> ==
!     ==  GRADGRHO(I)=<NABLA*F(I)|RHO>  GRADGG(I,J)=<NABLA*F(I)|F(J)> ==
!     ==================================================================
                              CALL TIMING$CLOCKON('ISOLATE_A_DOT')
      DO I=1,NAT*NFCT
!       ____GRHO(I)=<F(I)|w|RHO>________________________________________
        CALL ISOLATE_DOT(NGS,G,F(1,I),RHO,VOL,GRHO(I),GRADGRHO(1,I))
        DO J=I,NAT*NFCT
!         ____GG(I,J)=<F(I)|w|F(J)>_____________________________________
          CALL ISOLATE_DOT(NGS,G,F(1,I),F(1,J),VOL,GG(I,J),GRADGG(1,I,J))
          GG(J,I)=GG(I,J)
          GRADGG(1,J,I)=-GRADGG(1,I,J)
          GRADGG(2,J,I)=-GRADGG(2,I,J)
          GRADGG(3,J,I)=-GRADGG(3,I,J)
        ENDDO
      ENDDO
      CALL MPE$COMBINE('+',GRHO)
      CALL MPE$COMBINE('+',GRADGRHO)
      CALL MPE$COMBINE('+',GG)
      CALL MPE$COMBINE('+',GRADGG)
                              CALL TIMING$CLOCKOFF('ISOLATE_A_DOT')
!
!     ==================================================================
!     ==  COEFFICIENTS FOR THE FIT FUNCTIONS WITHOUT CONSTRAINT       ==
!     ==  C= [<G|G>]**(-1)*<G|RHO>       |RHO-BAR>=|G>QI              ==
!     ==================================================================
!
!     ____GGINV =[<G|G>]**(-1)__________________________________________
      call lib$invertr8(nat*nfct,gg,gginv)
!
!     ____QI= [<G|G>]**(-1)*<G|RHO>______________________________________
      CALL DGEMV('N',NAT*NFCT,NAT*NFCT,1.D0,GGINV,NAT*NFCT &
     &          ,GRHO,1,0.D0,QI,1)
!
!     ==================================================================
!     == CALCULATE QUALITY OF THE FIT                                 ==
!     ==================================================================
      IF(TPR) THEN
        SUM=0.D0
        DO IG=1,NGS
          CSVAR=RHO(IG)
          DO I=1,NAT*NFCT
            CSVAR=CSVAR-QI(I)*F(IG,I)
          ENDDO
          SUM=SUM+CSVAR*CONJG(CSVAR)
!         WRITE(*,FMT='(I3,("(",3F10.5,")"),2("(",2F12.7,")"))')
!    &         IG,(G(I,IG),I=1,3),RHO(IG),CSVAR
        ENDDO
        SUM=DSQRT(SUM)/DBLE(NGS)
        PRINT*,'LEAST SQUARE ',SUM
      END IF
!
!     ==================================================================
!     ==  ADD CONSTRAINT OF CHARGE CONSERVATION                       ==
!     ==================================================================
!     == XI=FORCE OF CONSTRAINT NORMALIZED TO ONE UNIT CHARGE ===========      
!     == XI=GGINV*(GGAMMA*VOL)/ (GGAMMA*VOL*GGINV*GGAMMA*VOL)
      SUM=0.D0
      DO I=1,NAT*NFCT
        SVAR=0.D0
        DO J=1,NAT*NFCT
          SVAR=SVAR+GGINV(I,J)*GGAMMA(J)*VOL
        ENDDO
        XI(I)=SVAR
        SUM=SUM+XI(I)*GGAMMA(I)*VOL
      ENDDO
!     ____XI(I)=XI(I)/SUM_________________________________________________
      do i=1,nat*nfct
        xi(i)=xi(i)/sum
      enddo
!
!     == APPLY CONSTRAINT TO CHARGES ===================================
      SUM=RHOGAMMA*VOL
      DO I=1,NAT*NFCT
        SUM=SUM-QI(I)*GGAMMA(I)*VOL
      ENDDO
!     __QI(I)=QI(I)+XI(I)*SUM______________________________________________
      do i=1,nat*nfct
        qi(i)=qi(i)+xi(i)*sum
      enddo
!
!     ==================================================================
!     ==  GRADIENT OF QI                                               ==
!     ==  REMARK: THIS IS NOT THE DERIVATIVE OF THE UNCONSTRAINED     ==
!     ==  QI, BUT ALSO CONTAINS PART OF THE LAGRANGE MULTIPLYER TERM!! ==
!     ==  (LOOKS WRONG, BUT IS CORRECT)                               ==
!     ==================================================================
!     ____GRAD QI=[<F|F>]**(-1)*GRAD<F|RHO> _____________________________
!     __________-[<F|F>]**(-1)*GRAD<F|F>[<F|F>]**(-1)*<F|RHO>___________
!     __GRADQI(I,I2,I1)=GRAD(I,I2) QI_I1__________________________________
      DO I1=1,NAT*NFCT
        DO I2=1,NAT*NFCT
          DO I=1,3
            GRADQI(I,I2,I1)=0.D0
          ENDDO
        ENDDO
      ENDDO
      ICOUNT=0
      DO I1=1,NAT*NFCT
        DO I2=1,NAT*NFCT
          ICOUNT=ICOUNT+1
          IF(MOD(ICOUNT-1,NTASKS).EQ.THISTASK-1) THEN
            DO I=1,3
              SVAR=GGINV(I1,I2)*GRADGRHO(I,I2)
              DO I3=1,NAT*NFCT
                SVAR=SVAR-GGINV(I1,I2)*GRADGG(I,I2,I3)*QI(I3) &
     &                   -GGINV(I1,I3)*GRADGG(I,I2,I3)*QI(I2)
              ENDDO
              GRADQI(I,I2,I1)=SVAR
            ENDDO
          END IF
        ENDDO
      ENDDO
      CALL MPE$COMBINE('+',GRADQI)
!
!     ==================================================================
!     ==  QI= [<F|F>]**(-1)*<F|RHO>       |RHO-BAR>=|F>QI               ==
!     ==================================================================
!     == APPLY CONSTRAINT TO GRADIENT OF CHARGES ==== 
      DO I=1,3
        DO I1=1,NAT*NFCT
          SUM=0.D0
          DO I2=1,NAT*NFCT
            SUM=SUM+GRADQI(I,I1,I2)*GGAMMA(I2)
          ENDDO
          SUM=SUM*VOL
          DO I2=1,NAT*NFCT
            GRADQI(I,I1,I2)=GRADQI(I,I1,I2)-XI(I2)*SUM
          ENDDO
        ENDDO
      ENDDO 
!
!     ==================================================================
!     == INTERFACE FOR COUPLING                                       ==
!     ==================================================================
      I=0
      DO IAT=1,NAT
        XQI(1,IAT)=QLM(1,IAT)/Y0
        XRC(1,IAT)=RCSM(IAT)
        DO IFCT=2,NFCT+1
          I=I+1
          XQI(IFCT,IAT)=QI(I)*GGAMMA(I)*VOL
          XRC(IFCT,IAT)=RC(I)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  SET CHARGES IN ATOMLIST AND GET EXTERNAL POTENTIALS        ==
!     ==================================================================
      DO IAT=1,NAT
        QMAD(IAT)=0.D0
        DO IFCT=1,NFCT+1
          QMAD(IAT)=QMAD(IAT)+XQI(IFCT,IAT)
        ENDDO
      ENDDO
      CALL ATOMLIST$SETR8A('Q',0,NAT,QMAD)
!
!     ==================================================================
!     ==  SET CHARGES IN ATOMLIST AND GET EXTERNAL POTENTIALS        ==
!     ==================================================================
      CALL ISOLATE_INTERFACE(RBAS,NAT,BAS,NFCT+1,XQI,XRC,E,XVI,FORCE)
!
!     ==================================================================
!     == CONTRIBUTION FROM THE POSITIVE BACKGROUND                    ==
!     ==================================================================
      FAC=-0.5D0*PI/VOL
      SUM=0.D0
      DO IAT1=1,NAT
        DO IFCT1=1,NFCT+1
          SVAR=0.D0
          DO IAT2=1,NAT
            DO IFCT2=1,NFCT+1
              SVAR=SVAR+FAC*(XRC(IFCT1,IAT1)**2+XRC(IFCT2,IAT2)**2) &
     &             *XQI(IFCT2,IAT2)
            ENDDO
          ENDDO
          SUM=SUM+SVAR*XQI(IFCT1,IAT1)
          XVI(IFCT1,IAT1)=XVI(IFCT1,IAT1)+2.D0*SVAR
        ENDDO
      ENDDO
      EBACKGROUND=SUM
      E=E+EBACKGROUND
!
!     ==================================================================
!     ==  BACK TRANSFORM                                              ==
!     ==================================================================
      I=0
      DO IAT=1,NAT
        VQLM(:,IAT)=0.D0
        VQLM(1,IAT)=XVI(1,IAT)/Y0
        DO IFCT=1,NFCT
          I=I+1
          VI(I)=XVI(1+IFCT,IAT)*GGAMMA(I)*VOL
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  DIPOLE MOMENT                                               ==
!     ==================================================================
      IF(TPR) THEN
        CALL CONSTANTS('DEBYE',DEBYE)
        DO I=1,3
          DIP(I)=0.D0
          DO IAT=1,NAT
            DIP(I)=DIP(I)+QMAD(IAT)*BAS(I,IAT)
          ENDDO
        ENDDO
        TDIP=DSQRT(DIP(1)**2+DIP(2)**2+DIP(3)**2)
        IF(TDIP.GT.1.D-5) THEN
          WRITE(*,FMT='("DIPOLE ",F10.5," DIR= ",3F3.0)') &
     &          TDIP/DEBYE,(DIP(I)/TDIP,I=1,3)
        ELSE
          WRITE(*,FMT='("DIPOLE =0")')
        END IF
        DO IAT=1,NAT
          if(lmrxx.lt.4) cycle   !attention when lmrx(isp)<lmrxx
          DIP(1)=DIP(1)+QLM(2,IAT)
          DIP(2)=DIP(2)+QLM(4,IAT)
          DIP(3)=DIP(3)+QLM(3,IAT)
        ENDDO
        TDIP=DSQRT(DIP(1)**2+DIP(2)**2+DIP(3)**2)
        IF(TDIP.GT.1.D-5) THEN
          WRITE(*,FMT='("DIPOLE ",F10.5," DIR= ",3F3.0)') &
     &            TDIP/DEBYE,(DIP(I)/TDIP,I=1,3)
        ELSE
          WRITE(*,FMT='("DIPOLE =0")')
        END IF
      END IF
!
      IF(TPR) THEN
        PRINT*,'E    ',E
        PRINT*,'BACKGROUND CORRECTION: ',EBACKGROUND
        WRITE(*,FMT='("QMAD",T10,10F10.5)')(QMAD(I),I=1,NAT)
        WRITE(*,FMT='("VMAD",T10,10F10.5)')(VMAD(I),I=1,NAT)
        WRITE(*,FMT='("QI",T10,10F10.5)') &
     &     (QI(I)*GGAMMA(I)*VOL,I=1,NAT*NFCT)
        WRITE(*,FMT='("VI",T10,10F10.5)')(VI(I),I=1,NAT*NFCT)
        WRITE(*,FMT='("QLM/Y0",T10,10F10.5)')(QLM(1,IAT)/Y0,IAT=1,NAT)
        WRITE(*,FMT='("VQLM",T10,10F10.5)')(VQLM(1,IAT),IAT=1,NAT)
      END IF
!
!     ==================================================================
!     == POTENTIAL                                                    ==
!     ==================================================================
      DO IG=1,NGS
        RHO(IG)=(0.D0,0.D0)
      ENDDO
!
!     == G=0 TERM ======================================================
      SUM=0.D0
      DO I=1,NAT*NFCT
        SUM=SUM+XI(I)*VI(I)
      ENDDO
      IF(NGAMMA.NE.0) RHO(NGAMMA)=SUM
!
!     ==  G.NE.0 TERM ===================================================
!     __ D= SUM_J V_J*DQ_J/D(GRHO_I)
!     ==  D=[<G|G>]**(-1) (V-CI<XI*VI>) =================================
      DO I=1,NAT*NFCT
        D(I)=VI(I)-SUM*GGAMMA(I)*VOL
      ENDDO
      ALLOCATE(WORK1(NAT*NFCT))
      CALL DCOPY(NAT*NFCT,D,1,WORK1,1)
      CALL DGEMV('N',NAT*NFCT,NAT*NFCT,1.D0,GGINV,NAT*NFCT &
     &          ,WORK1,1,0.D0,D,1)
      DEALLOCATE(WORK1)
!
      DO I=1,NAT*NFCT
        DO IG=1,NGS
          if(ig.eq.ngamma) cycle
          G2=G(1,IG)**2+G(2,IG)**2+G(3,IG)**2
          CALL ISOLATE_WEIGHTFUNCTION(G2,W)
          RHO(IG)=RHO(IG)+F(IG,I)*W*D(I)
        ENDDO                   
      ENDDO
!
!     ==================================================================
!     == FORCES                                                       ==
!     ==================================================================
      I1=0
      DO IAT=1,NAT
        DO IFCT=1,NFCT
          I1=I1+1
          DO I=1,3
            DO I2=1,NAT*NFCT
              FORCE(I,IAT)=FORCE(I,IAT)-GRADQI(I,I1,I2)*VI(I2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ISOLATE_INTERFACE(RBAS,NAT,POS,NG,QI,RC,ENERGY,VI,FORCE)
!     ******************************************************************
!     **                                                              **
!     **  THIS INTERFACE ALLOWS TO COUPLE ONTO                        **
!     **  POSITIONS AND POINTCHARGES OF THE ATOMS                     **
!     **                                                              **
!     **  CAUTION! ASSURE THAT FORCES POTENTIALS AND ENERGIES         **
!     **  ARE PROPERLY ADDED AND NOT RESET BY EXTERNAL ROUTINES       **
!     **                                                              **
!     ******************************************************************
      USE ISOLATE_MODULE, ONLY : TISOLATE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT           ! #(ATOMS)
      INTEGER(4),INTENT(IN) :: NG            ! #(GAUSSIANS PER SITE)
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)     ! CELL SIZE
      REAL(8)   ,INTENT(IN) :: POS(3,NAT)    ! ATOMIC POSITIONS
      REAL(8)   ,INTENT(IN) :: QI(NG,NAT)    ! CHARGES OF THE GAUSSIANS
      REAL(8)   ,INTENT(IN) :: RC(NG,NAT)    ! DECAY LENGTH OF GAUSSIANS
      REAL(8)   ,INTENT(OUT):: ENERGY        ! ENERGY 
      REAL(8)   ,INTENT(OUT):: FORCE(3,NAT)  ! FORCE
      REAL(8)   ,INTENT(OUT):: VI(NG,NAT)    ! POTENTIAL OF THE GAUSSIAN CHARGES
      REAL(8)               :: VI1(NG,NAT)   
      REAL(8)               :: ENERGY1       ! ENERGY 
      REAL(8)               :: FORCE1(3,NAT) ! FORCE
      REAL(8)               :: POT1(NAT)     ! POT
      REAL(8)               :: CHARGE(NAT)   ! ATOMIC POINT CHARGES
      INTEGER(4)            :: IAT,IG
!     ******************************************************************
      ENERGY=0.D0
      FORCE(:,:)=0.D0
      VI(:,:)=0.D0
!
!     ==================================================================
!     ==  ADD ALL GAUSSIAN CHARGES ON EACH ATOM TOGETHER
!     ==================================================================
      DO IAT=1,NAT
        CHARGE(IAT)=0.D0
        DO IG=1,NG
          CHARGE(IAT)=CHARGE(IAT)+QI(IG,IAT)
        ENDDO
      ENDDO
!
!     ==================================================================
!     == ISOLATE PERIODIC IMAGES                                      ==
!     ==================================================================
      IF(TISOLATE) THEN
        CALL ISOLATE_DECOUPLEPERIODIC(RBAS,NAT,POS,CHARGE,ENERGY1,POT1,FORCE1)
        ENERGY=ENERGY+ENERGY1
        DO IAT=1,NAT
          FORCE(:,IAT)=FORCE(:,IAT)+FORCE1(:,IAT)
          DO IG=1,NG
            VI(IG,IAT)=VI(IG,IAT)+POT1(IAT)
          ENDDO
        ENDDO
      END IF
!
!     ==================================================================
!     ==  COUPLE CLASSICAL ENVIRONMENT                                ==
!     ==================================================================
      CALL QMMM$INTERFACE(NAT,POS,CHARGE,FORCE1,POT1,ENERGY1)
      ENERGY=ENERGY+ENERGY1
      DO IAT=1,NAT
        FORCE(:,IAT)=FORCE(:,IAT)+FORCE1(:,IAT)
        DO IG=1,NG
          VI(IG,IAT)=VI(IG,IAT)+POT1(IAT)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  Calgary QM/MM interface                                     ==
!     ==================================================================
!     call mm_paw_md_interface(NG,RC,QI,VI1,FORCE1)
!     VI(:,:) = VI(:,:) + VI1(:,:)
!     FORCE(:,:)=FORCE(:,:)+FORCE1(:,:)
!
!     ==================================================================
!     ==  COUPLE CONTINUUM                                            ==
!     ==================================================================
      CALL CONTINUUM$PROPAGATE(NAT,POS,NG,RC,QI,ENERGY1,VI1,FORCE1)
      ENERGY=ENERGY+ENERGY1
      VI(:,:)=VI(:,:)+VI1(:,:)
      FORCE(:,:)=FORCE(:,:)+FORCE1(:,:)
!     THE FOLLOWING IS THE OLD CALL TO THE CONTINUUM, THE ONE
!     WHICH COUPLED TO THE POINT CHARGES ONLY. THE SOCKET WITHIN
!     PAW_CONTINUUM.F STILL EXISTS AND CAN BE HOOKED UP IF NECESSARY
!     CALL CONTINUUM$PROPAGATE(NAT,POS,CHARGE,ENERGY1,POT1,FORCE1)

      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ISOLATE_DECOUPLEPERIODIC(RBAS,NAT,POS,CHARGE,ENERGY,POT,FORCE)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: POS(3,NAT)
      REAL(8)   ,INTENT(IN) :: CHARGE(NAT)
      REAL(8)   ,INTENT(OUT):: ENERGY
      REAL(8)   ,INTENT(OUT):: FORCE(3,NAT)
      REAL(8)   ,INTENT(OUT):: POT(NAT)
      REAL(8)               :: ENERGY1
      REAL(8)               :: POT1(NAT)
      REAL(8)               :: FORCE1(3,NAT)
      REAL(8)               :: ENERGY2
      REAL(8)               :: POT2(NAT)
      REAL(8)               :: FORCE2(3,NAT)
      integer(4)            :: i,iat
!     ******************************************************************
!
!     ==================================================================
!     == INTERACTION OF CHARGES IN THE CLUSTER                        ==
!     ==================================================================
      CALL ISOLATE_CELLMADELUNG(NAT,POS,CHARGE,ENERGY1,POT1,FORCE1)
!
!     ==================================================================
!     == INTERACTION OF CHARGES IN THE CRYSTAL                        ==
!     ==================================================================
                              CALL TIMING$CLOCKON('ISOLATE_MADELUNG')
      CALL MADELUNG(NAT,RBAS,POS,CHARGE,ENERGY2,POT2,FORCE2)
                              CALL TIMING$CLOCKOFF('ISOLATE_MADELUNG')
!
!     ==================================================================
!     == SUBTRACT CRYSTAL FROM CLUSTER                                ==
!     ==================================================================
!     ____POT=POT1-POT2_________________________________________________
      ENERGY=ENERGY1-ENERGY2
      DO IAT=1,NAT
        POT(IAT)=POT1(IAT)-POT2(IAT)
      ENDDO
      DO IAT=1,NAT
        DO I=1,3
          FORCE(I,IAT)=FORCE1(I,IAT)-FORCE2(I,IAT)
        ENDDO
      ENDDO
!     CALL DVES(NAT,POT1,1,POT2,1,POT,1)
!     ____FORCE=FORCE1-FORCE2___________________________________________
!     CALL DVES(3*NAT,FORCE1,1,FORCE2,1,FORCE,1)
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE ISOLATE_LOWERNUMBER(G2MAX,NG,G2,NGS)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: G2MAX
      INTEGER(4),INTENT(IN) :: NG
      REAL(8)   ,INTENT(IN) :: G2(NG)
      INTEGER(4),INTENT(OUT):: NGS
      INTEGER(4)            :: IG
!     ******************************************************************
      NGS=0
      DO IG=1,NG
        IF(G2(IG).LT.G2MAX) NGS=NGS+1
      ENDDO  
      RETURN
      END     
!
!     ...............................................................
      SUBROUTINE ISOLATE_MAPLOWER(G2MAX,NG,G2,NGS,MAP)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: G2MAX
      INTEGER(4),INTENT(IN) :: NG
      REAL(8)   ,INTENT(IN) :: G2(NG)
      INTEGER(4),INTENT(IN) :: NGS
      INTEGER(4),INTENT(OUT):: MAP(NGS)
      INTEGER(4)            :: IG,IGS
!     ******************************************************************
      IGS=0
      DO IG=1,NG
        IF(G2(IG).LT.G2MAX) THEN
          IGS=IGS+1
          IF(IGS.GT.NGS) THEN
            CALL ERROR$STOP('ISOLATE_MAPLOWER')
          ENDIF
          MAP(IGS)=IG
        END IF
      ENDDO  
      IF(IGS.NE.NGS) THEN
        CALL ERROR$MSG('INCONSISTENT NUMBER OF SHORT VECTORS')
        CALL ERROR$STOP('ISOLATE_MAPLOWER')
      END IF
      RETURN
      END     
!
!     ..................................................................
      SUBROUTINE ISOLATE_WEIGHTFUNCTION(G2,W)
      USE ISOLATE_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: G2
      REAL(8)   ,INTENT(OUT):: W
!     ******************************************************************
      IF(G2.GT.1.D-8) THEN
        W=(G2-G2MAX)**2/G2
      ELSE
        W=0.D0
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ISOLATE_DOT(NG,G,F1,F2,OMEGA,RES,GRAD)
!     ******************************************************************
!     **                                                              **
!     **  res=<f1|w(G)|f2>                                            **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: NG
      REAL(8)    ,INTENT(IN) :: G(3,NG)
      REAL(8)    ,INTENT(IN) :: OMEGA           ! UNIT-CELL VOLUME
      COMPLEX(8),INTENT(IN)  :: F1(NG),F2(NG)
      REAL(8)   ,INTENT(OUT) :: RES
      REAL(8)   ,INTENT(OUT) :: GRAD(3)
      COMPLEX(8)             :: CSVAR
      REAL(8)                :: G2,SVAR,w
      INTEGER(4)             :: IG
!     ******************************************************************
      RES=0.D0
      GRAD(1)=0.D0
      GRAD(2)=0.D0
      GRAD(3)=0.D0
      DO IG=1,NG
        G2=G(1,IG)**2+G(2,IG)**2+G(3,IG)**2
        CALL ISOLATE_WEIGHTFUNCTION(G2,W)
        CSVAR=F1(IG)*CONJG(F2(IG))*W
        RES=RES+REAL(CSVAR)
        SVAR=AIMAG(CSVAR) 
        GRAD(1)=GRAD(1)+G(1,IG)*SVAR
        GRAD(2)=GRAD(2)+G(2,IG)*SVAR
        GRAD(3)=GRAD(3)+G(3,IG)*SVAR
      ENDDO
      RES=RES*2.D0*OMEGA
      GRAD(1)=GRAD(1)*2.D0*OMEGA
      GRAD(2)=GRAD(2)*2.D0*OMEGA
      GRAD(3)=GRAD(3)*2.D0*OMEGA
      RETURN
      END
!      
!     ......................................................MADELUNG...
      SUBROUTINE ISOLATE_CELLMADELUNG(NBAS,BAS,Q,EMAD,VMAD,FMAD)
!     ******************************************************************
!     **                                                              **
!     ** INTERACTION OF POINT CHARGES WITHIN ONE UNIT CELL            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NBAS
      REAL(8)   ,INTENT(IN) :: BAS(3,NBAS)
      REAL(8)   ,INTENT(IN) :: Q(NBAS)
      REAL(8)   ,INTENT(OUT):: EMAD
      REAL(8)   ,INTENT(OUT):: FMAD(3,NBAS)
      REAL(8)   ,INTENT(OUT):: VMAD(NBAS)
      INTEGER(4)            :: I,J,IAT,IAT1,IAT2
      REAL(8)               :: Q12,DX,DY,DZ,DLEN,RFAC1,RFAC2
!     ******************************************************************
!
      DO I=1,NBAS
        VMAD(I)=0.D0
        DO J=1,3
          FMAD(J,I)=0.D0
        ENDDO
      ENDDO
!
!     ==================================================================
!     == R-SPACE SUM                                                  ==
!     ==================================================================
      DO IAT1=1,NBAS
        DO IAT2=1,NBAS
          Q12=Q(IAT1)*Q(IAT2)
          DX=BAS(1,IAT2)-BAS(1,IAT1)
          DY=BAS(2,IAT2)-BAS(2,IAT1)
          DZ=BAS(3,IAT2)-BAS(3,IAT1)
          DLEN=DSQRT(DX*DX+DY*DY+DZ*DZ)
          IF(IAT1.NE.IAT2) THEN
            RFAC1=1.D0/DLEN
            RFAC2=-1.D0/DLEN**3
            RFAC1=0.5D0*RFAC1
            RFAC2=0.5D0*RFAC2*Q12
            VMAD(IAT1)=VMAD(IAT1)+RFAC1*Q(IAT2)
            VMAD(IAT2)=VMAD(IAT2)+RFAC1*Q(IAT1)
            FMAD(1,IAT1)=FMAD(1,IAT1)-RFAC2*DX     
            FMAD(2,IAT1)=FMAD(2,IAT1)-RFAC2*DY     
            FMAD(3,IAT1)=FMAD(3,IAT1)-RFAC2*DZ     
            FMAD(1,IAT2)=FMAD(1,IAT2)+RFAC2*DX     
            FMAD(2,IAT2)=FMAD(2,IAT2)+RFAC2*DY     
            FMAD(3,IAT2)=FMAD(3,IAT2)+RFAC2*DZ     
          END IF              
        ENDDO
      ENDDO
!
!     PRINT*,'VMAD IN ISOLATE_CELLMADELUNG',VMAD
!     PRINT*,'Q    IN ISOLATE_CELLMADELUNG',Q
!     PRINT*,'BAS  IN ISOLATE_CELLMADELUNG',BAS
!     PRINT*,'FMAD IN ISOLATE_CELLMADELUNG',FMAD
!
!     ==================================================================
!     == MADELUNG ENERGY                                              ==
!     ==================================================================
      EMAD=0.D0
      DO IAT=1,NBAS
        EMAD=EMAD+0.5D0*Q(IAT)*VMAD(IAT)
      ENDDO
!
!     ==================================================================
!     == TURN DERIVATIVES INTO FORCES                                 ==
!     ==================================================================
      DO IAT=1,NBAS
        DO I=1,3
          FMAD(I,IAT)=-FMAD(I,IAT)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE DIPOLE(RBAS,NGX,NG,G,RHO &
     &                ,NSP,NAT,ISPECIES,POS,LMRXX,LMRX,QLM)
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
      IMPLICIT none
      REAL(8)   ,intent(in) :: RBAS(3,3)
      integer(4),intent(in) :: ngx
      integer(4),intent(in) :: ng
      REAL(8)   ,intent(in) :: G(3,NGX)
      COMPLEX(8),intent(in) :: RHO(NG)
      integer(4),intent(in) :: nsp
      integer(4),intent(in) :: nat
      integer(4),intent(in) :: ispecies(nat)
      real(8)   ,intent(in) :: POS(3,NAT)
      integer(4),intent(in) :: lmrxx
      integer(4),intent(in) :: lmrx(nsp)
      REAL(8)   ,intent(in) :: QLM(LMRXX,NAT)
      REAL(8)               :: D(3)
      real(8)               :: pi,fourpi
      real(8)               :: y0
      real(8)               :: debye
      REAL(8)               :: GBAS(3,3)
      real(8)               :: vol
      real(8)               :: rmax
      real(8)               :: qtot,qttot
      real(8)               :: dtot
      real(8)               :: fac
      real(8)               :: f
      real(8)               :: grmax
      real(8)               :: g1
      integer(4)            :: i,ig,iat
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/DSQRT(FOURPI)
      CALL CONSTANTS('DEBYE',DEBYE)
      CALL GBASS(RBAS,GBAS,VOL)
      WRITE(*,FMT='("DIPOLE MOMENT EVALUATION")')
      WRITE(*,FMT='("========================")')
      WRITE(*,FMT='("INTEGRATES DENSITY WITHIN A SPHERE CENTERED AT THE ORIGIN")')
      WRITE(*,FMT='("THE VOLUME OF THE SPHERE EQUALS THE VOLUME OF THE UNIT CELL")')
      WRITE(*,FMT='("ASSUMES THAT ALL ATOMIC SITES ARE LOCATED WITHIN THIS SPHERE")')
      WRITE(*,FMT='("AND NOT ITS PERIODIC IMAGES!")')
!
!     ==================================================================
!     ==   DIPOLE MOMENT FROM THE PSEUDO DENSITY AND THE PSEUDO CORE  ==
!     ==================================================================
      RMAX=(3.D0*VOL/FOURPI)**(1.D0/3.D0)
      QTOT=0.D0
      DO I=1,3
        D(I)=0.D0
      ENDDO
      QTOT=QTOT+REAL(RHO(1))*FOURPI/3.D0*RMAX**3
      QTTOT=REAL(RHO(1))*VOL
      DO IG=2,NG
        G1=DSQRT(G(1,IG)**2+G(2,IG)**2+G(3,IG)**2)
        GRMAX=G1*RMAX
        F=(DSIN(GRMAX)-GRMAX*DCOS(GRMAX))/GRMAX**3
        QTOT=QTOT+FOURPI*F*RMAX**3*(2.D0*REAL(RHO(IG)))
!
        F=(3.D0*DSIN(GRMAX)-3.D0*GRMAX*DCOS(GRMAX)-GRMAX**2*DSIN(GRMAX)) &
     &      /GRMAX**5
        FAC=FOURPI*F*RMAX**5*(-2.D0*AIMAG(RHO(IG)))
        D(1)=D(1)+G(1,IG)*FAC
        D(2)=D(2)+G(2,IG)*FAC
        D(3)=D(3)+G(3,IG)*FAC
      ENDDO
      DTOT=DSQRT(D(1)**2+D(2)**2+D(3)**2)
!
      write(*,fmt='("plane wave contributrion to dipole moment:")')
      PRINT*,' QTOT ',QTOT,' QTTOT ',QTTOT
      WRITE(*,FMT='("dipole moment:",F10.5,"DEBYE; DIRECTION: (",3F3.0,")")') &
     &      DTOT/DEBYE,(D(I)/DTOT,I=1,3)
!
!     ==================================================================
!     ==  DIPOLE MOMENT                                               ==
!     ==================================================================
      DO IAT=1,NAT
        if(lmrx(ISPECIES(IAT)).le.4) cycle
        FAC=QLM(1,IAT)/Y0
        QTOT=QTOT+FAC
        QTTOT=QTTOT+FAC
        D(1)=D(1)+QLM(2,IAT)
        D(2)=D(2)+QLM(4,IAT)
        D(3)=D(3)+QLM(3,IAT)
        DO I=1,3
          D(I)=D(I)+POS(I,IAT)*FAC
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  printout                                                    ==
!     ==================================================================
      DTOT=DSQRT(D(1)**2+D(2)**2+D(3)**2)
      D(1)=D(1)/DTOT
      D(2)=D(2)/DTOT
      D(3)=D(3)/DTOT
!
      write(*,fmt='("complete dipole moment:")')
      PRINT*,' QTOT ',QTOT,'QTTOT ',QTTOT
      WRITE(*,FMT='("Dipole moment",F10.5,"DEBYE; DIRECTION: (",3F3.0,")")') &
     &      DTOT/DEBYE,D
      RETURN
      END
