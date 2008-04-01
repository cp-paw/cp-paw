!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOIN,EXC,VXC,excarr)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE EXCHANGE AND CORRELATION ENERGY                      **
!     **  FOR A DENSITY GIVEN ON A RADIAL LOGARITHMIC GRID                    **
!     **  TIMES REAL SPHERICAL HARMONICS                                      **
!     **                                                                      **
!     **  THE TOTAL ENERGY IS AN EXPANSION ABOUT THE                          **
!     **  SPHERICAL CONTRIBUTION OF THE DENSITY UP TO QUADRATIC               **
!     **  ORDER IN THE NON-SPHERICAL CONTRIBUTIONS                            **
!     **                                                                      **
!     **  EXC = EXC(XVAL(L=0)*Y0)                                             **
!     **      + 0.5 * D2[EXC]/D[XVAL(L=0)*Y0]**2 * XVAL(L>0)**2               **
!     **                                                                      **
!     **  WHERE XVAL=(/RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS/)              **
!     **  IS AN SPHERICAL HARMONICS EXPANSION ON THE RADIAL GRID.             **
!     **                                                                      **
!     **  DEPENDECIES:                                                        **
!     **    DFT                                                               **
!     **    TIMING                                                            **
!     **    TRACE                                                             **
!     **                                                                      **
!     **  REMARKS: THE GRADIENTS ARE CORRECT ONLY IF DFT SUPPORTS             **
!     **    THIRD DERIVATIVES OF THE XC-ENERGY                                **
!     **   - WHEN USING SELFTEST ON THIS ROUTINE, THEN                        **
!     **     D(EXC)/DRHO(I)=POT(I)*DEX*R(I)**3                                **
!     **     AND THE VALUES AT LARGE RADII MUST BE SURPRESSED                 **
!     **                                                                      **
!     **  REMARK: FOR A COLLINEAR DENSITY THE ROUTINE GIVES DIFFERENT RESULTS **
!     **          WITH NDIMD=2 AND NDIMD=4 DUE TO THE TAYLOR EXPANSION IN     **
!     **          ANGULAR MOMENTUM EXPANSIONS                                 **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1996 ************
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TNS=.TRUE. ! NON-SPHERICAL CONTRIBUTIONS ON
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD      ! CAN BE 1,2,4
      REAL(8)   ,INTENT(IN) :: RHOIN(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: VXC(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(OUT):: excarr(NR)
      LOGICAL(4)            :: TGRA   ! SWITCH FOR GRADIENT CORRECTION
      INTEGER(4)            :: NSPIN
      REAL(8)               :: EXC1
      REAL(8)               :: R(NR)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: GRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VGRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: B(:,:)    
      REAL(8)               :: B0INV(NR) 
      REAL(8)               :: ROOTB0(NR) 
      REAL(8)   ,ALLOCATABLE:: C(:,:)    
      REAL(8)   ,ALLOCATABLE:: VB(:,:)    
      REAL(8)   ,ALLOCATABLE:: VC(:,:)    
      REAL(8)               :: VAL5(5),VXC5(5),V2XC5(5,5),V3XC5(5,5,5)
      REAL(8)               :: XVAL(NR,5,LMRX)
      REAL(8)               :: XDER(NR,5,LMRX)
      REAL(8)               :: PI,FOURPI
      REAL(8)               :: Y0
      INTEGER(4)            :: IR,L,II,ISPIN,ISPIN1,ISPIN2,I,J
      INTEGER(4)            :: LM,LM1,LM2,LM3
      INTEGER(4)            :: IMAX
      REAL(8)               :: RI,FAC
      REAL(8)               :: CG0LL
      REAL(8)               :: CG
      REAL(8)               :: SVAR
      REAL(8)               :: WORK(NR)
      REAL(8)               :: WORK1(NR)
      REAL(8)               :: WORK2(NR)
      REAL(8)               :: WORK3(NR)
!     **************************************************************************
      CALL TRACE$PUSH('AUGMENTATION_XC')
      EXC=0.D0
      VXC(:,:,:)=0.D0
!
!     ==========================================================================
!     ==   CALCULATE SOME CONSTANTS NEEDED LATER                              ==
!     ==========================================================================
      CALL DFT$GETL4('GC',TGRA)
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      CG0LL=Y0
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  OBTAIN SPIN DENSITY                                                 ==
!     ==========================================================================
      NSPIN=1
      IF(NDIMD.GT.1) NSPIN=2
      ALLOCATE(RHO(NR,LMRX,NSPIN))
      ALLOCATE(GRHO(NR,LMRX,NSPIN))
      ALLOCATE(VRHO(NR,LMRX,NSPIN))
      RHO(:,:,1)=RHOIN(:,:,1)
      IF(NDIMD.EQ.2) THEN
        RHO(:,:,2)=RHOIN(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
!       == HERE WE NEED TO CALCULATE THE ABSOLUTE VALUE OF THE SPIN DENSITY ====
!       == IN AN ANGULAR MOMENTUM EXPANSION. THIS IS ONLY POSSIBLE APPROXIMATELY
!       == USING A TAYLOR EXPANSION ABOUT THE SPHERICAL PART OF THE SQUARE =====
!       == OF THE SPIN DENSITY. ================================================
        VRHO(:,:,:)=0.D0
        CALL AUGMENTATION_NCOLLTRANS('RHO',NR,LMRX,RHOIN,RHO,VRHO,VXC)
      END IF
!
!     == IMAX ALLOWS TO RESTRICT SOME LOOPS (1:5) TO (1:IMAX)
      IF(TGRA) THEN
        IF(NSPIN.EQ.2) THEN; IMAX=5; ELSE; IMAX=3; END IF
      ELSE 
        IF(NSPIN.EQ.2) THEN; IMAX=2; ELSE; IMAX=1; END IF
      END IF
!
!     ==========================================================================
!     ==  CALCULATE RADIAL GRADIENT OF THE DENSITY                            ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE GRADIENTS')
      IF(TGRA) THEN
        GRHO(:,:,:)=0.D0
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            CALL RADIAL$DERIVE(GID,NR,RHO(:,LM,ISPIN),GRHO(:,LM,ISPIN))
          ENDDO
        ENDDO
      ELSE
        GRHO(:,:,:)=0.D0
      END IF
!
!     ==========================================================================
!     ==  DEFINE VECTOR (RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS)             ==
!     ==========================================================================
      XVAL(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
          XVAL(:,ISPIN,LM)=RHO(:,LM,ISPIN)
        ENDDO
      ENDDO
      IF(TGRA) THEN
        II=2
        DO ISPIN1=1,NSPIN          ! THIS LOOP PUTS T,T->3; S,S->4 ;T,S->5
          DO ISPIN2=ISPIN1,1,-1    ! AND ASSURES CONSISTENCY WITH NSPIN
            II=II+1
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=DBLE(L*(L+1))
              XVAL(:,II,1)=XVAL(:,II,1) &
        &         +CG0LL*(GRHO(:,LM,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                +FAC*RHO(:,LM,ISPIN1)*RHO(:,LM,ISPIN2)/R(:)**2)
            ENDDO
            DO LM=2,LMRX 
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              XVAL(:,II,LM)=XVAL(:,II,LM) &
        &         +0.5D0*CG0LL*(GRHO(:,1,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                      +GRHO(:,LM,ISPIN1)*GRHO(:,1,ISPIN2))
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  CALCULATE EXCHANGE ENERGY FOR THE SPHERICAL DENSITY                 ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE DFT')
      WORK1(:)=0.D0
      XDER(:,:,:)=0.D0
      DO IR=1,NR
!       ==  CYCLE IF THE TOTAL DENSITY VANISHES ================================
        IF(XVAL(IR,1,1).LE.0.D0) CYCLE
!       == NOW CALL DFT ROUTINE ================================================
        VAL5(:)=XVAL(IR,:,1)*Y0
        CALL DFT3(VAL5,EXC1,VXC5,V2XC5,V3XC5)
!       == NOW CALCULATE ENERGY DENSITY AND DERIAVTIVES ========================
        WORK1(IR)=FOURPI*EXC1
        XDER(IR,:,1)  =FOURPI*VXC5(:)*Y0
!PRINT*,IR,R(IR),VAL5(1),VXC5(1),V2XC5(1,1),V3XC5(1,1,1)
        DO LM=2,LMRX
          DO I=1,IMAX        ! IMAX=<5 
            DO J=1,IMAX
              WORK1(IR)=WORK1(IR) &
       &              +0.5D0*XVAL(IR,I,LM)*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,:,1)=XDER(IR,:,1) &
       &              +0.5D0*Y0*XVAL(IR,I,LM)*V3XC5(:,I,J)*XVAL(IR,J,LM)
              XDER(IR,I,LM)=XDER(IR,I,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,J,LM)=XDER(IR,J,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,I,LM)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      excarr(:)=work1(:)/(fourpi*y0)
      CALL RADIAL$INTEGRAL(GID,NR,WORK1(:)*R(:)**2,EXC)
!
!     ==========================================================================
!     ==  TRANSFORM POTENTIALS FOR SPHERICAL PART                             ==
!     ==========================================================================
      ALLOCATE(VGRHO(NR,LMRX,NSPIN))
      VRHO(:,:,:)=0.D0
      VGRHO(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          VRHO(:,LM,ISPIN)=XDER(:,ISPIN,LM)
        ENDDO
      ENDDO
      IF(TGRA) THEN
        II=2
        DO ISPIN1=1,NSPIN
          DO ISPIN2=ISPIN1,1,-1
            II=II+1
!           == FIRST RESOLVE XVAL(:,II,1) ======================================
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=DBLE(L*(L+1))
              VRHO(:,LM,ISPIN1)  =VRHO(:,LM,ISPIN1) &
      &                 +CG0LL*FAC/R(:)**2*XDER(:,II,1)*RHO(:,LM,ISPIN2)
              VRHO(:,LM,ISPIN2)  =VRHO(:,LM,ISPIN2) &
      &                 +CG0LL*FAC/R(:)**2*XDER(:,II,1)*RHO(:,LM,ISPIN1)
              VGRHO(:,LM,ISPIN1) =VGRHO(:,LM,ISPIN1) &
      &                 +CG0LL*XDER(:,II,1)*GRHO(:,LM,ISPIN2)
              VGRHO(:,LM,ISPIN2) =VGRHO(:,LM,ISPIN2) &
      &                 +CG0LL*XDER(:,II,1)*GRHO(:,LM,ISPIN1)
            ENDDO
!           == NOW RESOLVE XVAL(:,II,LM) =======================================
            DO LM=2,LMRX
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              VGRHO(:,1,ISPIN1) =VGRHO(:,1,ISPIN1) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,LM,ISPIN2)
              VGRHO(:,1,ISPIN2) =VGRHO(:,1,ISPIN2) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,LM,ISPIN1)
              VGRHO(:,LM,ISPIN2)=VGRHO(:,LM,ISPIN2) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,1,ISPIN1)
              VGRHO(:,LM,ISPIN1)=VGRHO(:,LM,ISPIN1) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,1,ISPIN2)
            ENDDO
          ENDDO
        ENDDO               
      END IF
!
!     ==========================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS                     ==
!     ==  V = V -1/R**2 D/DR [ R**2 VGRHO ]                                   ==
!     ==  V = V -[2/R VGRHO+ D/DR VGRHO ]                                     ==
!     ==========================================================================
      IF(TGRA) THEN
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
!           == FIRST ALTERNATIVE
!           CALL RADIAL$DERIVE(GID,NR,VGRHO(:,LM,ISPIN),WORK2)   !NOT SO 
!           WORK1(:)=2.D0/R(:)*VGRHO(:,LM,ISPIN)+WORK2(:)  !GOOD
!           ==  SECOND ALTERNATIVE APPEARS TO BE MORE ACCURATE
            WORK2(:)=VGRHO(:,LM,ISPIN)*R(:)**2
            CALL RADIAL$DERIVE(GID,NR,WORK2,WORK1)
            WORK1(:)=WORK1(:)/R(:)**2
!           == ALTERNATIVES FINISHED
            VRHO(:,LM,ISPIN)=VRHO(:,LM,ISPIN)-WORK1(:)
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(VGRHO)
!
!     ==========================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS                     ==
!     ==========================================================================
      VXC(:,:,1)=VRHO(:,:,1)
      IF(NDIMD.EQ.2) THEN
        VXC(:,:,2)=VRHO(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
        CALL AUGMENTATION_NCOLLTRANS('POT',NR,LMRX,RHOIN,RHO,VRHO,VXC)
      END IF     
      DEALLOCATE(RHO)
      DEALLOCATE(GRHO)
      DEALLOCATE(VRHO)
!
!     ==========================================================================
!     ==   CORRECT FOR DIVERGENCE AT THE ORIGIN:                              ==
!     ==   IF A SHIFTED LOGARITHMIC GRID IS USED THE FIRST GRID POINT         ==
!     ==   IS MESSED UP BECAUSE OF FACTORS 1/R                                ==
!     ==========================================================================
      IF(R(1).LT.1.D-5) THEN
        VXC(1,:,:)=VXC(2,:,:)
      END IF
                      CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_NCOLLTRANS(ID,NR,LMRX,RHO4,RHO2,POT2,POT4)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4)  ,INTENT(IN)    :: NR
      INTEGER(4)  ,INTENT(IN)    :: LMRX
      REAL(8)     ,INTENT(IN)    :: RHO4(NR,LMRX,4)
      REAL(8)     ,INTENT(OUT)   :: RHO2(NR,LMRX,2)
      REAL(8)     ,INTENT(IN)    :: POT2(NR,LMRX,2)
      REAL(8)     ,INTENT(OUT)   :: POT4(NR,LMRX,4)
      REAL(8)                    :: A(NR,LMRX,3)
      REAL(8)                    :: VA(NR,LMRX,3)
      REAL(8)                    :: Q(NR)
      REAL(8)                    :: VQ(NR)
      REAL(8)                    :: PI,Y0
      REAL(8)     ,PARAMETER     :: R8SMALL=1.D-20
      INTEGER(4)                 :: ISIG,LM
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      RHO2(:,:,:)=0.D0
      POT4(:,:,:)=0.D0
!
!     ==========================================================================
!     == WORK OUT AUXILIARY VARIABLES                                         ==
!     ==========================================================================
      Q(:)=SQRT(RHO4(:,1,2)**2+RHO4(:,1,3)**2+RHO4(:,1,4)**2)
!
      DO ISIG=1,3
        DO LM=1,LMRX
          A(:,LM,ISIG)=RHO4(:,LM,ISIG+1)/(Q(:)+R8SMALL)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == WORK OUT COLLINEAR DENSITY                                           ==
!     ==========================================================================
      RHO2(:,:,1)=RHO4(:,:,1)
      RHO2(:,1,2)=Q(:)
      DO LM=2,LMRX
        RHO2(:,LM,2)=0.D0
        DO ISIG=1,3
          RHO2(:,LM,2)=RHO2(:,LM,2)+A(:,1,ISIG)*A(:,LM,ISIG)
        ENDDO
        RHO2(:,LM,2)=RHO2(:,LM,2)*Q(:)
      ENDDO
!
!     ==========================================================================
!     == RETURN IF ONLY DENSITY IS REQUIRED                                   ==
!     ==========================================================================
      IF(ID.EQ.'RHO') RETURN
!
!     ==========================================================================
!     == WORK OUT POTENTIAL                                                   ==
!     ==========================================================================
      IF(ID.NE.'POT') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AUGMENTATION_NCTRANS1')
      END IF
!     -- CALCULATE VQ ----------------------------------------------------------
      VQ(:)=0.D0
      DO LM=1,LMRX
        VQ(:)=VQ(:)+POT2(:,LM,2)*RHO2(:,LM,2)/(Q(:)+R8SMALL)
      ENDDO
!     -- CALCULATE VA ----------------------------------------------------------
      DO ISIG=1,3
        VA(:,1,ISIG)=0.D0
        DO LM=2,LMRX
          VA(:,1,ISIG)=VA(:,1,ISIG)+POT2(:,LM,2)*A(:,LM,ISIG)
          VA(:,LM,ISIG)=POT2(:,LM,2)*A(:,1,ISIG)
        ENDDO
      ENDDO
!     -- NEW VQ ----------------------------------------------------------------
      DO ISIG=1,3
        DO LM=1,LMRX
          VQ(:)=VQ(:)-VA(:,LM,ISIG)*A(:,LM,ISIG)
        ENDDO
      ENDDO
!     --------------------------------------------------------------------------
      POT4(:,:,1)=POT2(:,:,1)
      DO ISIG=1,3
        DO LM=1,LMRX
          POT4(:,LM,ISIG+1)=VA(:,LM,ISIG)
        ENDDO
        POT4(:,1,ISIG+1)=POT4(:,1,ISIG+1)+VQ(:)*A(:,1,ISIG)
      ENDDO
      RETURN
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      subroutine spherical$transform(gid,nr,ro,lmx,lmx3,lx1,f,ft)
!!$!     **************************************************************************    
!!$!     **  latex spherical-harmonics-reexp.tex
!!$!     **  reexpands a spherical harmonics expansion centered @ ro around 
!!$!     **  the origin --> the shift happens from ro to O
!!$!     **  rotates the offsite expansion situated @ ro in such a way that 
!!$!     **  the new local z axis goes through the global origin.
!!$!     **  in a second step, we shift (reexpand) the expansion to the origin
!!$!     **  the third step rotates the local z axis back to the global one
!!$!     **************************************************************************    
!!$      implicit none
!!$      integer(4),intent(in) :: gid               !grid identifier
!!$      integer(4),intent(in) :: nr                !#(radial grid points)
!!$      real(8)   ,intent(in) :: ro(3)             !the vector to the offsite center
!!$      integer(4),intent(in) :: lmx               !(il+1)**2 lmx for the original expansion
!!$      integer(4),intent(in) :: lmx3              !(il3+1)**2 lmx for the new expansion @ O
!!$      integer(4),intent(in) :: lx1               !l'max for the radial expansion (note:m'=0 per definition)
!!$      real(8)   ,intent(in) :: f(nr,lmx)         !the original radial parts 
!!$      real(8)   ,intent(out):: ft(nr,lmx3)       !the transformed radial parts 
!!$      integer(4)            :: lx                !the max values for l
!!$      integer(4)            :: lx2
!!$      integer(4)            :: lx3
!!$      real(8)               :: dtheta
!!$      integer(4),parameter  :: ntheta=300
!!$      real(8)               :: fr(nr,lmx)        !the rotated radial parts
!!$      real(8)               :: frt(nr,lmx3)      !the transformed rad. parts bef. backrotation
!!$      real(8),allocatable   :: ylm1(:)
!!$      real(8),allocatable   :: intexp(:,:)       !the original radial parts 
!!$      real(8)               :: dist
!!$      real(8)               :: cg
!!$      real(8)               :: rot(3,3)          !rotates the z axis onto the connecting line \bar{Oro}
!!$      real(8)               :: rotinv(3,3)       !rotates back
!!$      real(8)               :: ylmrot(lmx,lmx)   !rotates the radial part
!!$      real(8)               :: ylmrot3(lmx3,lmx3)!rotates the reexp. radial parts back
!!$      real(8)               :: dummy(nr)
!!$      real(8)               :: pi
!!$      real(8)               :: rvar,rvar2,rvar3
!!$      real(8)               :: xexp
!!$      real(8)               :: ri 
!!$      real(8)               :: r(nr)
!!$      integer(4)            :: il,il1,il2,il3
!!$      integer(4)            :: im,im1,im2,im3
!!$      integer(4)            :: lm,lm1,lm2,lm3
!!$      integer(4)            :: ilm,jlm
!!$      integer(4)            :: ival
!!$      integer(4)            :: i,j
!!$      integer(8)            :: isvar
!!$      real(8)               :: SQfact(0:20)  
!!$!     **************************************************************************    
!!$      pi=4.d0*datan(1.d0)
!!$      call radial$r(gid,nr,r)
!!$      dtheta=pi/real(ntheta,kind=8)
!!$! dtheta=1.d-2
!!$!
!!$!     ==========================================================================
!!$!     == just copy the density for a '0' shift transformation ==================
!!$!     ==========================================================================
!!$      if(abs(ro(1)).lt.1.d-8.and.abs(ro(2)).lt.1.d-8.and.abs(ro(3)).lt.1.d-8) then
!!$        if(lmx.gt.lmx3) then
!!$          print*,'PROBLEMS IN SPHERICAL$TRANSFORM: we have a 0 transformtion'
!!$          print*,'and the original expansion lmx is larger than the new lmx3'
!!$          print*,'which leads to precission loss'
!!$          stop'SPHERICAL$TRANSFORM'
!!$        end if
!!$        ft(:,:)=0.d0
!!$        ft(:,1:lmx)=f(:,1:lmx)
!!$        return
!!$      end if
!!$!
!!$!     ==========================================================================
!!$!     == set the max l values                                                 ==
!!$!     ==========================================================================
!!$      if(int(sqrt(real(lmx,kind=8)))**2.lt.lmx) then
!!$        print*,'YOU HAVE TO USE A FULL SHELL IN TRANSFORM_SPH'
!!$        stop
!!$      end if
!!$      if(int(sqrt(real(lmx3,kind=8)))**2.lt.lmx3) then
!!$        print*,'YOU HAVE TO USE A FULL SHELL IN TRANSFORM_SPH'
!!$        stop
!!$      end if
!!$      lx=int(sqrt(real(lmx,kind=8)))-1
!!$      lx2=lx  !?fixed from l-l''' coupling
!!$      lx3=int(sqrt(real(lmx3,kind=8)))-1 !we use the incoming value
!!$
!!$      if(lmx3.ne.(lx3+1)**2) then
!!$        print*,'THE DIMENSION OF ft DOES NOT FIT TO LMX3'
!!$        print*,(lx3+1)**2,lmx3
!!$        stop
!!$      end if
!!$      dist=sqrt(dot_product(ro(:),ro(:)))
!!$  !print*,'warning using debug dist!!!'
!!$  !dist=-ro(1)
!!$!
!!$!     ==========================================================================
!!$!     == calculate the SQUARE ROOT OF THE factorial                           ==
!!$!     == Factorial(20) is the highest representable as INTEGER(8)             ==
!!$!     ==========================================================================
!!$      sqfact(0)=1.d0
!!$      isvar=1
!!$      do i=1,20
!!$        isvar=isvar*i
!!$        SQfact(i)=SQRT(real(isvar))
!!$      enddo
!!$!
!!$!     ==========================================================================
!!$!     ==  rotate coordinate systes so that distance vector points along z     ==
!!$!     ==========================================================================
!!$      call rotate_base(ro(:),rot)
!!$!     ==  for debug: the e-14 scattering influences the potofdensity 
!!$!     ==  (on the same scale -> no problem for production, but for debugging)
!!$      if(.true.) then
!!$         do i=1,3
!!$            do j=1,3
!!$              if(abs(rot(i,j)).lt.1e-13)rot(i,j)=0.d0
!!$           end do
!!$        end do
!!$      end if
!!$!
!!$!     ==  obtain transformation of real spherical harmonics
!!$      call ROTATEYLM(lmx,rot,ylmrot(:,:))
!!$
!!$!     ==  see peters methods 
!!$!     ==  note: needs a fix: errorstop for sqrt(lmx) not integer because it 
!!$!     ==  returns only 0 for partly filled shells!
!!$!     ==  clac the new radial functions
!!$!
!!$!     ==  fr is the rotated  function
!!$      fr(:,:)=0.d0
!!$      do ilm=1,lmx    !loop over the components of the rotated vector
!!$        do jlm=1,lmx !the ilm component of the lmx vector for a given r after rotation
!!$          fr(:,ilm)=fr(:,ilm)+ylmrot(ilm,jlm)*f(:,jlm)
!!$        enddo
!!$      enddo
!!$!
!!$!     ==========================================================================
!!$!     ==                                                                      ==
!!$!     ==========================================================================
!!$      allocate(ylm1(maxval((/lmx,lmx3,(lx1+1)**2/))))
!!$      allocate(intexp(nr,(lx+1)**2))
!!$
!!$  !get all possible ylm's for Yl''m''(R)
!!$      call getylm((lx2+1)**2,(/0.d0,0.d0,1.d0/),ylm1(:)) !we have R at the z axis ->theta,phi=0 
!!$
!!$
!!$  !==================== SHIFT PARALLEL TO NEW Z AXIS ====================
!!$      frt(:,:)=0.d0
!!$      do il1=0,lx1
!!$!       == remember: im1=0
!!$        im1=0    !because of delta m',0 (cylinder symmetry)!
!!$        lm1=il1**2+il1+1+im1
!!$!       == call radialintegralpart() --> intexp(nr,lx) for every nr of |r|
!!$        call radialintegralpart(gid,nr,dtheta,(lx+1)**2,il1,dist,fr(:,:),intexp)
!!$
!!$        do il=0,lx
!!$          do im=-il,il
!!$            lm=il**2+il+1+im
!!$            do il2=0,lx2
!!$!             ==  CLEBSCH routine gives an error for 0 argument (1st part)
!!$              if(il-il2.lt.0) cycle
!!$              rvar2=dist**il2*(-1.d0)**il2
!!$              do im2=-il2,il2
!!$                lm2=il2**2+il2+1+im2
!!$!               ==  CLEBSCH routine gives an error for 0 argument (2nd part)
!!$                if((il-il2)**2+(il-il2)+1+(im-im2).lt.1 &
!!$     &             .or.abs(im-im2).gt.(il-il2)) cycle
!!$!               ==         l',m'=0,    do not substract combined L,L''!!!
!!$!               ==  the factorial factor
!!$                rvar=sqrt(4.d0*pi*(2.d0*real(il,kind=8)+1.d0) &
!!$     &                   /real((2*il-2*il2+1)*(2*il2+1),kind=8)) &
!!$     &               *sqfact(il-im)*SQfact(il+im) &
!!$     &               /(sqfact(il2-im2)*sqfact(il2+im2) &
!!$     &                    *sqfact(il-il2-im+im2)*sqfact(il-il2+im-im2))
!!$                do il3=0,lx3
!!$                  do im3=-il3,il3
!!$                    lm3=il3**2+il3+1+im3
!!$                    call CLEBSCH(lm1,(il-il2)**2+(il-il2)+1+(im-im2),lm3,CG)
!!$                    if(cg.eq.0.d0) cycle
!!$!                   == get the m'=0 term
!!$!                   == ylm1((il1)**2+il1+1)    !the same for every |r|
!!$!
!!$!                   ==  bring it togehter:
!!$                    dummy(:)=intexp(:,lm)*r(:)**(il-il2)
!!$                    dummy(:)=dummy(:)*ylm1(lm2)*rvar*cg*rvar2
!!$                    frt(:,lm3)=frt(:,lm3)+dummy(:)
!!$                  end do
!!$                end do
!!$              end do
!!$            end do
!!$          end do
!!$        end do
!!$      end do
!!$!
!!$!     ==========================================================================
!!$!     ==  rotate coordinate system back to original orientation               ==
!!$!     ==========================================================================
!!$!     == a transpose should replace the inversion (rot is unitary!)
!!$      call LIB$INVERTR8(3,rot(:,:),rotinv(:,:))
!!$!     == is rotateylm unitary as well?
!!$      call ROTATEYLM(lmx3,rotinv,ylmrot3(:,:))
!!$!     ==  see peters methods
!!$!     ==  note: needs a fix: errorstop for sqrt(lmx) not integer because it returns only 0
!!$!     ==  for partly filled shells!
!!$!     ==clac the new radial functions
!!$      ft(:,:)=0.d0
!!$      do ilm=1,lmx3    !loop over the components of the rotated vector
!!$        do jlm=1,lmx3 !the ilm component of the lmx vector for a given r after rotation
!!$          ft(:,ilm)=ft(:,ilm)+ylmrot3(ilm,jlm)*frt(:,jlm)
!!$        enddo
!!$      enddo
!!$      return
!!$    end subroutine spherical$transform
!
module sph_module
  interface operator(.fact.)
     module procedure factorial
  end interface

  contains
    function factorial(n) result(ival)
      integer(4),intent(in)      :: n
      integer(4)                 :: ival
      integer(4)                 :: i
      ival=1
      do i=2,n
         ival=ival*i
      end do
    end function factorial
end module sph_module



subroutine spherical$transform(gid,nr,ro,lmx,lmx3,f,ft)
  !latex spherical-harmonics-reexp.tex
  !reexpands a spherical harmonics expansion centered @ ro around the origin
  !--> the shift happens from ro to O
 !rotates the offsite expansion situated @ ro in such a way that the new local z axis 
  !goes through the global origin.
  !in a second step, we shift (reexpand) the expansion to the origin
  !the third step rotates the local z axis back to the global one
  use sph_module
  implicit none
  integer(4),intent(in)     :: gid
  integer(4),intent(in)     :: nr                !grid data
  real(8),intent(in)        :: ro(3)             !the vector to the offsite center
  integer(4),intent(in)     :: lmx               !(il+1)**2 lmx for the original expansion
  integer(4),intent(in)     :: lmx3              !(il3+1)**2 lmx for the new expansion @ O
  integer(4)                :: lx1               !l'max for the radial expansion (note:m'=0 per definition)
  real(8),intent(in)        :: f(nr,lmx)         !the original radial parts 
  real(8),intent(out)       :: ft(nr,lmx3)       !the transformed radial parts 
  integer(4)                :: lx                !the max values for l
  integer(4)                :: lx2
  integer(4)                :: lx3
  real(8)                   :: dtheta
  integer(4)                :: ntheta=300
  real(8)                   :: fr(nr,lmx)        !the rotated radial parts
  real(8)                   :: frt(nr,lmx3)      !the transformed rad. parts bef. backrotation
  real(8),allocatable       :: ylm1(:)
  real(8),allocatable       :: intexp(:,:)       !the original radial parts 
  real(8)                   :: dist
  real(8)                   :: cg
  real(8)                   :: rot(3,3)          !rotates the z axis onto the connecting line \bar{Oro}
  real(8)                   :: rotinv(3,3)       !rotates back
  real(8)                   :: ylmrot(lmx,lmx)   !rotates the radial part
  real(8)                   :: ylmrot3(lmx3,lmx3)!rotates the reexp. radial parts back
  real(8)                   :: dummy(nr)
  real(8)                   :: pi
  real(8)                   :: rvar,rvar2,rvar3
  real(8)                   :: xexp
  real(8)                   :: ri 
  integer(4)                :: il,il1,il2,il3
  integer(4)                :: im,im1,im2,im3
  integer(4)                :: ilm,jlm
  integer(4)                :: ival
  integer(4)                :: ir,i,j
  real(8)                   :: r(nr)
#ifdef TRACEINTO
  print*,'INTO spherical$transform'
  print*,'    with lmx,lmx3    =',lmx,lmx3
#endif
call radial$r(gid,nr,r)


  !--- just copy the density for a '0' shift transformation
  if(abs(ro(1)).lt.1.d-8.and.abs(ro(2)).lt.1.d-8.and.abs(ro(3)).lt.1.d-8) then
     if(lmx.gt.lmx3) then
        print*,'PROBLEMS IN SPHERICAL$TRANSFORM: we have a 0 transformtion'
        print*,'and the original expansion lmx is larger than the new lmx3'
        print*,'which leads to precission loss'
        stop'SPHERICAL$TRANSFORM'
     end if
     ft(:,:)=0.d0
     ft(:,1:lmx)=f(:,1:lmx)
     return
  end if

  !set the max l values
  if(int(sqrt(real(lmx,kind=8)))**2.lt.lmx) then
     print*,'YOU HAVE TO USE A FULL SHELL IN TRANSFORM_SPH'
     stop
  end if
  if(int(sqrt(real(lmx3,kind=8)))**2.lt.lmx3) then
     print*,'YOU HAVE TO USE A FULL SHELL IN TRANSFORM_SPH'
     stop
  end if

  ft(:,:)=0.d0
  frt(:,:)=0.d0
  fr(:,:)=0.d0

  lx=int(sqrt(real(lmx,kind=8)))-1
  lx2=lx  !?fixed from l-l''' coupling
  lx3=int(sqrt(real(lmx3,kind=8)))-1 !we use the incoming value

  if(lmx3.ne.(lx3+1)**2) then
     print*,'THE DIMENSION OF ft DOES NOT FIT TO LMX3'
     print*,(lx3+1)**2,lmx3
     stop
  end if



  !======================================
  !=== fix lx1 to lx+lx3 (Auswahlregel)
  !======================================
  print*,'NOTE: WE FIX lx1 to lx+lx3 (Auswahlregel) in transform'
  lx1=lx+lx3



  pi=4.d0*datan(1.d0)
  dist=sqrt(dot_product(ro(:),ro(:)))
  !print*,'warning using debug dist!!!'
  !dist=-ro(1)


  !=== get dtheta
  dtheta=pi/real(ntheta,kind=8)



  !==================== ROTATION ====================
  call rotate_base(ro(:),rot)

 !for debug: the e-14 scattering influences the potofdensity (on the same scale ->
  ! no problem for production, but for debugging)
  if(.true.) then
     do i=1,3
        do j=1,3
           if(abs(rot(i,j)).lt.1e-13)rot(i,j)=0.d0
        end do
     end do
  end if

  call ROTATEYLM(lmx,rot,ylmrot(:,:))
  !see peters methods
  !note: needs a fix: errorstop for sqrt(lmx) not integer because it returns only 0
  !for partly filled shells!
  !clac the new radial functions
  do ir=1,nr         !do this for every radial gridpoint
     do ilm=1,lmx    !loop over the components of the rotated vector
        do jlm=1,lmx !the ilm component of the lmx vector for a given r after rotation
           fr(ir,ilm)=fr(ir,ilm)+ylmrot(ilm,jlm)*f(ir,jlm)
        enddo
     enddo
  enddo





!!$ri=r1/xexp
!!$do ir=1,nr
!!$   ri=ri*xexp
!!$   print*,ri,f(ir,2),f(ir,3),f(ir,4)
!!$end do
!!$print*,'-------------------------'
!!$
!!$ri=r1/xexp
!!$do ir=1,nr
!!$   ri=ri*xexp
!!$   print*,ri,fr(ir,2),fr(ir,3),fr(ir,4)
!!$end do
!!$stop

!debug: no rotation
!fr(:,:)=f(:,:)
  allocate(ylm1(maxval((/lmx,lmx3,(lx1+1)**2/))))
  allocate(intexp(nr,(lx+1)**2))

  !get all possible ylm's for Yl''m''(R)
  call getylm((lx2+1)**2,(/0.d0,0.d0,1.d0/),ylm1(:)) !we have R at the z axis ->theta,phi=0 


  !==================== SHIFT PARALLEL TO NEW Z AXIS ====================
  im1=0    !because of delta m',0 (cylinder symmetry)
  do il1=0,lx1
     !remeber: im1=0
     !call radialintegralpart() --> intexp(nr,lx) for every nr of |r|
     call radialintegralpart(gid,nr,dtheta,(lx+1)**2,il1,dist,fr(:,:),intexp)

     do il=0,lx
        do im=-il,il

           do il2=0,lx2
              !CLEBSCH routine gives an error for 0 argument (1st part)
              if(il-il2.lt.0) cycle
              rvar2=dist**il2*(-1.d0)**il2

              do im2=-il2,il2

                 !CLEBSCH routine gives an error for 0 argument (2nd part)
                 if((il-il2)**2+(il-il2)+1+(im-im2).lt.1&
                      &.or.abs(im-im2).gt.(il-il2)) cycle
                 !             l',m'=0,    do not substract combined L,L''!!!
                 
                 
                 !the factorial factor
                 rvar=sqrt(4.d0*pi*(2.d0*real(il,kind=8)+1.d0)/&
                      &real((2*il-2*il2+1)*(2*il2+1),kind=8))*&
                      &sqrt(real(.fact.(il-im)*.fact.(il+im),kind=8)/real(.fact.(il2-im2)*.fact.(il2+im2),kind=8))&
                      &/sqrt(real(.fact.(il-il2-im+im2)*.fact.(il-il2+im-im2),kind=8))
                 
                 do il3=0,lx3
                    do im3=-il3,il3

                       call CLEBSCH(il1**2+il1+1,(il-il2)**2+(il-il2)+1+(im-im2),il3**2+il3+1+im3,CG)

                       if(cg.eq.0.d0) cycle

                       !get the m'=0 term
                       !ylm1((il1)**2+il1+1)    !the same for every |r|


                       !bring it togehter:
                       do ir=1,nr
                          dummy(ir)=intexp(ir,il**2+il+1+im)*r(ir)**(il-il2)
                       end do
                       dummy(:)=dummy(:)*ylm1(il2**2+il2+1+im2)*rvar*cg*rvar2
                       
frt(:,il3**2+il3+1+im3)=frt(:,il3**2+il3+1+im3)+dummy(:)
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do



  !==================== BACK ROTATION ====================
  call LIB$INVERTR8(3,rot(:,:),rotinv(:,:))
  call ROTATEYLM(lmx3,rotinv,ylmrot3(:,:))
  !see peters methods
  !note: needs a fix: errorstop for sqrt(lmx) not integer because it returns only 0
  !for partly filled shells!
  !clac the new radial functions
  do ir=1,nr         !do this for every radial gridpoint
     do ilm=1,lmx3    !loop over the components of the rotated vector
        do jlm=1,lmx3 !the ilm component of the lmx vector for a given r after rotation
           ft(ir,ilm)=ft(ir,ilm)+ylmrot3(ilm,jlm)*frt(ir,jlm)
        enddo
     enddo
  enddo




!debug: no rotation
!ft(:,:)=frt(:,:)
  return
end subroutine spherical$transform

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine radialintegralpart(gid,nr,dtheta,lmx,l1,RC,f,val)
!     **************************************************************************    
!     **************************************************************************    
      implicit none
      integer(4),intent(in)           :: gid
      integer(4),intent(in)           :: nr
      real(8),intent(in)              :: dtheta       !spwdidth for the integration
      integer(4),intent(in)           :: lmx          !we do the integrals for all lmx
      integer(4),intent(in)           :: l1           !l'
      real(8),intent(in)              :: RC           !dist to old center
      real(8),intent(in)              :: f(nr,lmx)    !old radial part
      real(8),intent(out)             :: val(nr,lmx)  !value of the integral for each |r| of the radial grid
                                                  !and all lmx
      integer(4)                      :: itheta
      integer(4)                      :: ir
      integer(4)                      :: ilm
      real(8)                         :: pi
      real(8)                         :: thetai
      real(8)                         :: r(nr)
      real(8),allocatable             :: yl10(:)
      real(8),allocatable             :: ylmtemp(:)
      real(8)                         :: x,fx
      real(8)                         :: dummy(nr)
!     **************************************************************************    
      call radial$r(gid,nr,r)
      pi=4.d0*atan(1.d0)
  !check the shell
      if(int(sqrt(real(lmx,kind=8)))**2.lt.lmx) then
        print*,'YOU HAVE TO USE A FULL SHELL IN RADIALINTEGRALPART'
        stop
      end if

  !======= prepare Yl',0 THIS HAS TO BE OPTIMIZED! =======
      allocate(yl10(int(pi/dtheta)))
      allocate(ylmtemp((l1+1)**2))
      thetai=-dtheta/2.d0
      do itheta=1,int(pi/dtheta)
        thetai=thetai+dtheta
        call getylm((l1+1)**2,(/sin(thetai),0.d0,cos(thetai)/),ylmtemp(:))
        yl10(itheta)=ylmtemp((l1)**2+l1+1) !we need only the m=0 term
      end do

      val(:,:)=0.d0
      do ilm=1,lmx !we loop over the l values
        thetai=-dtheta/2.d0
        do itheta=1,int(pi/dtheta)
          thetai=thetai+dtheta
          do ir=1,nr !for all possible radial values 
            x=sqrt(RC**2+r(ir)**2-2.d0*r(ir)*RC*cos(thetai))
            !interpolate the value
            call RADIAL$VALUE(gid,NR,f(:,ilm),x,fx)
            dummy(ir)=fx*1.d0/x**(ceiling(sqrt(real(ilm,kind=8)))-1)
          end do
          val(:,ilm)=val(:,ilm)+dtheta*yl10(itheta)*sin(thetai)*dummy(:)
        end do
      end do
      val(:,:)=val(:,:)*2.d0*pi
      return
    end subroutine radialintegralpart
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine rotate_base(vec,rot)
!     **************************************************************************    
!     **************************************************************************    
  !=====================================================
  !== determines the rotation matrix such, vec'=rot*vec
  !== with vec'=(/0,0,|vec|/)
  !== this corresponds to a rotation of the basis with theta
  !== and phi coming from the angles of the vector
  !=====================================================
      implicit none
      real(8),intent(in)         :: vec(3)
      real(8),intent(out)        :: rot(3,3)
      integer(4)                 :: i,ivec(1)
      real(8)                    :: dis,dr1(3),DR2(3),DR3(3)
      real(8),parameter          :: ex(3)=(/1.d0,0.d0,0.d0/)
      real(8),parameter          :: ey(3)=(/0.d0,1.d0,0.d0/)
      real(8),parameter          :: ez(3)=(/0.d0,0.d0,1.d0/)
!     **************************************************************************    
      dis=sqrt(sum(vec(:)**2))
      if(dis.eq.0.d0) then
        rot(:,:)=0.d0
        do i=1,3
          rot(i,i)=1.d0
        enddo
      end if
      dr1(:)=vec(:)/dis
      ivec=minloc(abs(dr1))
      i=ivec(1)
      if(i.eq.1) Then
        dr2(1)=dr1(2)*ex(3)-dr1(3)*ex(2)
        dr2(2)=dr1(3)*ex(1)-dr1(1)*ex(3)
        dr2(3)=dr1(1)*ex(2)-dr1(2)*ex(1)
      else if(i.eq.2) Then
        dr2(1)=dr1(2)*ey(3)-dr1(3)*ey(2)
        dr2(2)=dr1(3)*ey(1)-dr1(1)*ey(3)
        dr2(3)=dr1(1)*ey(2)-dr1(2)*ey(1)
      else if(i.eq.3) Then
        dr2(1)=dr1(2)*ez(3)-dr1(3)*ez(2)
        dr2(2)=dr1(3)*ez(1)-dr1(1)*ez(3)
        dr2(3)=dr1(1)*ez(2)-dr1(2)*ez(1)
      end if
      dr2(:)=dr2(:)/sqrt(sum(dr2(:)**2))
      dr3(1)=dr1(2)*dr2(3)-dr1(3)*dr2(2)
      dr3(2)=dr1(3)*dr2(1)-dr1(1)*dr2(3)
      dr3(3)=dr1(1)*dr2(2)-dr1(2)*dr2(1)
!     ==  the strange mapping is for consistency reasons
      ROT(3,:)=DR1(:)
      ROT(2,:)=-DR2(:)
      ROT(1,:)=DR3(:)
      return
      end subroutine rotate_base
