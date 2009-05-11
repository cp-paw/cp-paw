!........1.........2.........3.........4.........5.........6.........7.........8
module radialPAW_module
type vpaw_type
  logical            :: ton=.false.
  integer(4)         :: gid
  integer(4)         :: nr
  integer(4)         :: lnx
  integer(4),pointer :: lox(:)
  real(8)   ,pointer :: dh(:,:)
  real(8)   ,pointer :: pspot(:)
  real(8)   ,pointer :: do(:,:)
  real(8)   ,pointer :: pro(:,:)
end type vpaw_type
end module radialPAW_module
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE radialpaw$makevpaw(gid,nr,lnx,lox,pspot,dh,do,pro,vpaw)
!     **************************************************************************
!     ** sets up a the input for calculating the fock potential               **
!     **                                                                      **
!     ** it is possible to subtract also a local potential vfock%mux(nr)      **
!     **   and to scale the result vfock%scale.                               **
!     **   this call initializes scale=1. and mux(:)=0.                       **
!     **************************************************************************
      use radialpaw_module, only : vpaw_type
      integer(4),intent(in) :: gid
      integer(4),intent(in) :: nr
      integer(4),intent(in) :: lnx
      integer(4),intent(in) :: lox(lnx)
      real(8)   ,intent(in) :: pspot(nr)
      real(8)   ,intent(in) :: dh(lnx,lnx)
      real(8)   ,intent(in) :: do(lnx,lnx)
      real(8)   ,intent(in) :: pro(nr,lnx)
      type(vpaw_type),intent(inout) :: vpaw
!     **************************************************************************
      if(vpaw%ton)call radialpaw$cleanvpaw(vpaw)
      vpaw%ton=.true.
      vpaw%gid=gid
      vpaw%nr=nr
      vpaw%lnx=lnx
      allocate(vpaw%lox(lnx))
      allocate(vpaw%dh(lnx,lnx))
      allocate(vpaw%do(lnx,lnx))
      allocate(vpaw%pro(nr,lnx))
      allocate(vpaw%pspot(nr))
      vpaw%lox(:)=lox(:)
      vpaw%dh(:,:)=dh(:,:)
      vpaw%do(:,:)=do(:,:)
      vpaw%pro(:,:)=pro(:,:)
      vpaw%pspot(:)=pspot(:)
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE radialpaw$cleanvpaw(vpaw)
!     **************************************************************************
!     ** prepares the coefficient array used for determining the 
!     **************************************************************************
      use radialpaw_module, only : vpaw_type
      implicit none
      type(vpaw_type),intent(inout) :: vpaw
!     **************************************************************************
      if(vpaw%lnx.eq.-1) return
      vpaw%ton=.false.
      vpaw%gid=-1
      vpaw%nr=0
      vpaw%lnx=0
      deallocate(vpaw%dh)
      deallocate(vpaw%do)
      deallocate(vpaw%pro)
      deallocate(vpaw%pspot)
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE radialpaw$vpsi(GID,NR,vpaw,l,F,vf,of)
!     **************************************************************************
!     **  APPLIES THE FOCK-OPERATOR TO A SET OF FUNCTIONS                     **
!     **    |g_i>= v_fock |f_i>
!     **************************************************************************
      use radialpaw_module, only : vpaw_type
      implicit none
      INTEGER(4)      ,INTENT(IN) :: GID
      INTEGER(4)      ,INTENT(IN) :: NR
      type(vpaw_type) ,intent(in) :: vpaw
      INTEGER(4)      ,INTENT(IN) :: L
      REAL(8)         ,INTENT(IN) :: F(NR)
      REAL(8)         ,INTENT(OUT):: vf(NR)
      REAL(8)         ,INTENT(OUT):: of(NR)
      real(8)         ,allocatable:: proj(:)
      real(8)         ,allocatable:: cpro(:)
      real(8)         ,allocatable:: dpro(:)
      integer(4)                  :: npro
      integer(4)                  :: lnx
      integer(4)                  :: ln,ln1,ln2
      integer(4)                  :: ipro,ipro1,ipro2
      real(8)                     :: aux(nr)
      real(8)                     :: r(nr)
      real(8)                     :: pi,y0
!     **************************************************************************
      lnx=vpaw%lnx
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      call radial$r(gid,nr,r)
!
!     ==========================================================================
!     == determine number of relevant projectors and allocate arrays          ==
!     ==========================================================================
      npro=0
      do ln=1,lnx
        if(vpaw%lox(ln).ne.l) cycle
        npro=npro+1
      enddo
      allocate(proj(npro))
      allocate(cpro(npro))
      allocate(dpro(npro))
!
!     ==========================================================================
!     == determine projections proj=<pro|f>                                   ==
!     ==========================================================================
      proj(:)=0.d0
      ipro=0
      do ln=1,lnx
        if(vpaw%lox(ln).ne.l) cycle
        ipro=ipro+1
        aux(:)=r(:)**2*vpaw%pro(:,ln)*f(:)
        call radial$integral(gid,nr,aux,proj(ipro))
      enddo
!
!     ==========================================================================
!     == determine cpro=dh<pro|f> and dpro=do<pro|f>                          ==
!     ==========================================================================
      cpro=0.d0
      dpro=0.d0
      ipro1=0
      do ln1=1,lnx
        if(vpaw%lox(ln1).ne.l) cycle
        ipro1=ipro1+1
        ipro2=0
        do ln2=1,lnx
          if(vpaw%lox(ln2).ne.l) cycle
          ipro2=ipro2+1
          cpro(ipro1)=cpro(ipro1)+vpaw%dh(ln1,ln2)*proj(ipro2)
          dpro(ipro1)=dpro(ipro1)+vpaw%do(ln1,ln2)*proj(ipro2)
        enddo
      enddo
!
!     ==========================================================================
!     == determine vf and of                                                  ==
!     ==========================================================================
      vf(:)=vpaw%pspot(:)*y0*f(:)
      of(:)=f(:)
      ipro=0
      do ln=1,lnx
        if(vpaw%lox(ln).ne.l) cycle
        ipro=ipro+1
        vf(:)=vf(:)+vpaw%pro(:,ln)*cpro(ipro)
        of(:)=of(:)+vpaw%pro(:,ln)*dpro(ipro)
      enddo
!
!     ==========================================================================
!     == close down                                                           ==
!     ==========================================================================
      deallocate(proj)
      deallocate(cpro)
      deallocate(dpro)
      return
      end
!
!........1.........2.........3.........4.........5.........6.........7.........8
module radialfock_module
type vfock_type
  logical            :: ton=.false.
  integer(4)         :: gid
  integer(4)         :: nr
  real(8)            :: scale
  integer(4)         :: lrhox
  integer(4)         :: nb
  integer(4),pointer :: l(:)
  integer(4),pointer :: so(:)
  real(8)   ,pointer :: f(:)
  real(8)   ,pointer :: psi(:,:)
  real(8)   ,pointer :: mux(:)
end type vfock_type
integer(4)               :: lxcoeff=-1
real(8)     ,allocatable :: coeff(:,:,:)
end module radialfock_module
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE radialfock$makevfock(gid,nr,nb,lofi,sofi,fofi,psi &
     &                               ,rbox,lrhox,vfock)
!     **************************************************************************
!     ** sets up a the input for calculating the fock potential               **
!     **                                                                      **
!     ** it is possible to subtract also a local potential vfock%mux(nr)      **
!     **   and to scale the result vfock%scale.                               **
!     **   this call initializes scale=1. and mux(:)=0.                       **
!     **************************************************************************
      use radialfock_module
      implicit none
      integer(4)      ,intent(in) :: gid
      integer(4)      ,intent(in) :: nr
      integer(4)      ,intent(in) :: nb
      integer(4)      ,intent(in) :: lofi(nb)
      integer(4)      ,intent(in) :: sofi(nb)
      real(8)         ,intent(in) :: fofi(nb)
      integer(4)      ,intent(in) :: lrhox
      real(8)         ,intent(in) :: rbox
      real(8)         ,intent(in) :: psi(nr,nb)
      type(vfock_type),intent(inout) :: vfock
      real(8)                     :: r(nr)
      real(8)                     :: aux(nr),svar
      integer(4)                  :: ir,ib
!     **************************************************************************
      if(vfock%ton)call radialfock$cleanvfock(vfock)
      call radial$r(gid,nr,r)
      vfock%ton=.true.
      vfock%scale=1.d0
      vfock%lrhox=lrhox
      vfock%gid=gid
      vfock%nb=nb
      vfock%nr=nr
      allocate(vfock%l(nb))
      allocate(vfock%so(nb))
      allocate(vfock%f(nb))
      allocate(vfock%psi(nr,nb))
      allocate(vfock%mux(nr))
      vfock%l(:)=lofi(:)
      vfock%so(:)=sofi(:)
      vfock%f(:)=fofi(:)
      vfock%psi(:,:)=psi(:,:)
      vfock%mux(:)=0.d0
! 
!     ==========================================================================
!     == cut of tails                                                         ==
!     ==========================================================================
      do ir=1,nr
        if(r(ir).gt.rbox) then
          vfock%psi(ir:,:)=0.d0
          exit
        end if
      enddo
! 
!     ==========================================================================
!     == re-normalize                                                         ==
!     ==========================================================================
      do ib=1,nb
        aux(:)=(r(:)*vfock%psi(:,ib))**2
        call radial$integral(gid,nr,aux,svar)
        vfock%psi(:,ib)=psi(:,ib)/sqrt(svar)
      enddo
      return
      end SUBROUTINE radialfock$makevfock
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE radialfock$cleanvfock(vfock)
!     **************************************************************************
!     ** prepares the coefficient array used for determining the 
!     **************************************************************************
      use radialfock_module
      implicit none
      type(vfock_type),intent(inout) :: vfock
!     **************************************************************************
      if(vfock%nb.eq.-1) return
      vfock%ton=.false.
      vfock%scale=0.d0
      vfock%nb=-1
      vfock%gid=-1
      vfock%nr=-1
      vfock%lrhox=-1
      deallocate(vfock%l)
      deallocate(vfock%so)
      deallocate(vfock%f)
      deallocate(vfock%psi)
      deallocate(vfock%mux)
      return
      end SUBROUTINE radialfock$cleanvfock
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE radialfock_updatecoeff(lxnew)
!     **************************************************************************
!     ** prepares the coefficient array used for determining the 
!     **************************************************************************
      use radialfock_module
      implicit none
      integer(4)      ,intent(in) :: lxnew
      integer(4)                  :: lf,lb,lrho
      integer(4)                  :: lmf,lmb,lmrho
      integer(4)                  :: imb,imrho
      integer(4)                  :: lx
      real(8)                     :: svar
      real(8)                     :: cg
      real(8)                     :: pi
!     **************************************************************************
      if(lxnew.le.lxcoeff) return
!
      if(lxcoeff.ne.-1) deallocate(coeff)
      lxcoeff=lxnew
      lx=lxnew
      allocate(coeff(lx+1,2*lx+1,lx+1))
      pi=4.d0*atan(1.d0)
      COEFF(:,:,:)=0.D0
      DO LF=0,LX
        DO LB=0,LX
          DO LRHO=0,2*Lx
!          if(mod(lf+lb+lrho,2).eq.1) cycle !empirical
!PRINT*,'========================================================='
            LMF=LF**2+1  ! RESULT INDEPEMDENT OF IMF: PIC ANY LMF
            SVAR=0.D0
            LMB=LB**2
            DO IMB=1,2*LB+1
              LMB=LMB+1
              LMRHO=LRHO**2
              DO IMRHO=1,2*LRHO+1
                LMRHO=LMRHO+1
                CALL SPHERICAL$GAUNT(LMF,LMB,LMRHO,CG)
                SVAR=SVAR+CG**2
              ENDDO
            ENDDO
            SVAR=4.D0*PI*SVAR/REAL((2*LB+1)*(2*Lrho+1),KIND=8) !AVERAGE OVER M_J
!WRITE(6,FMT='(5I5,F20.5)')LF,LB,LRHO,IMB,IMRHO,SVAR
            COEFF(LB+1,LRHO+1,LF+1)=SVAR
          ENDDO
        ENDDO
      ENDDO
!
      return
      end SUBROUTINE radialfock_updatecoeff
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE radialfock$vpsi(GID,NR,vfock,l,F,G)
!     **************************************************************************
!     **  APPLIES THE FOCK-OPERATOR TO A SET OF FUNCTIONS                     **
!     **    |g_i>= v_fock |f_i>
!     **************************************************************************
      use radialfock_module, only : vfock_type,coeff
      implicit none
      INTEGER(4)      ,INTENT(IN) :: GID
      INTEGER(4)      ,INTENT(IN) :: NR
      type(vfock_type),intent(in) :: vfock
      INTEGER(4)      ,INTENT(IN) :: L
      REAL(8)         ,INTENT(IN) :: F(NR)
      REAL(8)         ,INTENT(OUT):: G(NR)
      REAL(8)                     :: PI,y0
      REAL(8)                     :: r(nr)
      REAL(8)                     :: SVAR
      INTEGER(4)                  :: LXB
      INTEGER(4)                  :: LX
      INTEGER(4)                  :: IB
      INTEGER(4)                  :: LB,LRHO
      REAL(8)                     :: AUX1(NR),AUX2(NR)
!     **************************************************************************
      g=0.d0
      if(.not.vfock%ton.or.vfock%scale.eq.0.d0) return
!
      PI=4.D0*ATAN(1.D0)
      y0=1.d0/sqrt(4.d0*pi)
      LXB=MAXVAL(vfock%l)
      lx=max(lxb,l,int((vfock%lrhox+1)/2))
      call radialfock_updatecoeff(lx)
      call radial$r(gid,nr,r)
!
!     ==========================================================================
!     == apply fock potential                                                 ==
!     ==========================================================================
      G(:)=0.D0
      DO IB=1,vfock%NB
        IF(vfock%f(IB).LT.1.D-3) CYCLE
        LB=vfock%L(IB)
        AUX1(:)=vfock%PSI(:,IB)*F(:)
        DO LRHO=0,vfock%LRHOX         
          IF(COEFF(LB+1,LRHO+1,L+1).LT.1.D-8) CYCLE    !COEFF.GE.0
          CALL RADIAL$POISSON(GID,NR,LRHO,AUX1,AUX2)
!         == factor 0.5 cancels spin multiplicity ============================
          SVAR=-0.5d0*vfock%f(ib)*COEFF(LB+1,LRHO+1,L+1) &
     &                           *REAL((2*LRHO+1),KIND=8)/(4.d0*pi)
          G(:)=G(:)+SVAR*vfock%psi(:,IB)*AUX2(:)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == g=scale*(vfock-mux*f)                                                **
!     ==========================================================================
      g(:)=g(:)-vfock%mux(:)*y0*f(:)
      g(:)=vfock%scale*g(:)
      RETURN
      END SUBROUTINE radialfock$vpsi
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$AESCF(GID,NR,KEY,RBOX,AEZ,NX,NB,LOFI,SOfi,Fofi,NNofi &
     &                        ,ETOT,POT,vfock,EOFI,PHI,SPHI)
!     **************************************************************************
!     ** MAKES A SELF-CONSISTENT CALCULATION OF AN ATOM IN A BOX WITH         **
!     ** RADIUS RBOX (RADIUS IS LIMITED BY THE GRID RBOX<R(NR-3) )            **
!     **                                                                      **
!     ** KEY='START,REL,SO,FOCK='                                             **
!     **                                                                      **
!     **************************************************************************
      use radialfock_module
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID       ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR        ! #(GRID POINTS)
      CHARACTER(*),INTENT(IN)    :: KEY     
      REAL(8)    ,INTENT(INOUT)  :: RBOX      ! BOX RADIUS     
      REAL(8)    ,INTENT(IN)     :: AEZ       ! ATOMIC NUMBER
      INTEGER(4) ,INTENT(IN)     :: NX        ! X#(STATES)
      INTEGER(4) ,INTENT(OUT)    :: NB        ! #(STATES)
      INTEGER(4) ,INTENT(INOUT)  :: LOFI(NX)  ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(INOUT)  :: SOfi(NX)    ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(INOUT)  :: Fofi(NX)     ! OCCUPATION
      INTEGER(4) ,INTENT(INOUT)  :: NNOFI(NX)    ! #(NODES)
      REAL(8)    ,INTENT(OUT)    :: ETOT      ! TOTAL ENERGY
      REAL(8)    ,INTENT(INOUT)  :: POT(NR)   ! POTENTIAL
      type(vfock_type),INTENT(INOUT)  :: vfock  ! fock term
      REAL(8)    ,INTENT(OUT)    :: EOFI(NX)  ! ONE-PARTICLE ENERGIES
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,NX)! ONE-PARTICLE WAVE FUNCTIONS
      REAL(8)    ,INTENT(OUT)    :: SPHI(NR,NX) ! SMALL COMPONENT
      REAL(8)                    :: R(NR)
      REAL(8)                    :: DREL(NR)  ! RELATIVISTIC CORRECTION
      REAL(8)                    :: RHO(NR)
      REAL(8)                    :: G(NR)     ! INHOMOGENEITY
      REAL(8)                    :: AUX(NR),AUX1(NR)   !AUXILIARY ARRAY
      REAL(8)                    :: MUX(NR)   !EXCHANGE ONLY POTENTIAL
      INTEGER(4)                 :: ITER
      INTEGER(4),PARAMETER       :: NITER=1000
      REAL(8)                    :: XAV,XMAX
      LOGICAL(4)                 :: CONVG
      REAL(8)   ,PARAMETER       :: TOL=1.D-6
      REAL(8)                    :: EKIN,EH,EXC
      LOGICAL(4),PARAMETER       :: TBROYDEN=.TRUE.
      REAL(8)                    :: POTIN(NR)
      REAL(8)                    :: SVAR
      INTEGER(4)                 :: I,IB,JB,ISO,L
      INTEGER(4)                 :: ISVAR,IARR(1)
      LOGICAL(4)                 :: TSTART  ! CALCULATE ANGULAR MOMENTA ETC
      LOGICAL(4)                 :: TREL    ! RELATIOVISTIC CALCULATION
      LOGICAL(4)                 :: TSO     ! CALCULATE WITH SPIN ORBIT COUPLING
      LOGICAL(4)                 :: TFOCK   ! CALCULATE WITH FOCK EXCHANGE
      INTEGER(4)          :: LMAP(19)=(/0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1/)
      REAL(8)                    :: FTOT
      REAL(8)                    :: PI,Y0,C0LL
      INTEGER(4)                 :: NBROYDENMEM
      REAL(8)                    :: BROYDENSTEP
      INTEGER(4)                 :: LRHOX=4
      character(128)             :: string
      real(8)                    :: scale
      logical                    :: tsecond
!     **************************************************************************
!
!     ==========================================================================
!     == RESOLVE KEY                                                          ==
!     ==========================================================================
      TREL=INDEX(KEY,'NONREL').EQ.0
      IF(TREL.AND.INDEX(KEY,'REL').EQ.0) THEN
        CALL ERROR$MSG('TREL=T BUT "REL" NOT SPECIFIED IN KEY')
        CALL ERROR$CHVAL('KEY',KEY)
        CALL ERROR$STOP('ATOMLIB$AESCF')
      END IF
      TSO=INDEX(KEY,'NONSO').EQ.0
      IF(TSO.AND.INDEX(KEY,'SO').EQ.0) THEN
        CALL ERROR$MSG('TSO=T BUT "SO" NOT SPECIFIED IN KEY')
        CALL ERROR$CHVAL('KEY',KEY)
        CALL ERROR$STOP('ATOMLIB$AESCF')
      END IF
      TFOCK=INDEX(KEY,'FOCK').NE.0
      IF(TFOCK) THEN
        ISVAR=INDEX(KEY,'FOCK')+5
        STRING=KEY(ISVAR:)
        ISVAR=INDEX(KEY,'FOCK')
        READ(STRING,*)SCALE
      ELSE
        scale=0.d0
      end if
print*,'key=',trim(key)
print*,'scale=',scale,' tfock=',tfock
      TSTART=INDEX(KEY,'START').NE.0
      NBROYDENMEM=2
      BROYDENSTEP=5.D-1
!
!     ==========================================================================
!     == INITIALIZATIONS                                                      ==
!     ==========================================================================
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=Y0
      CALL RADIAL$R(GID,NR,R)
!
      EOFI=0.D0
      IF(TSTART) THEN
        LOFI(:)=0
        Fofi(:)=0.D0
        SOFI(:)=0
        NB=0
        FTOT=AEZ
        DO I=1,19
          NB=NB+1
          IF(NB.GT.NX) THEN
            CALL ERROR$MSG('ACTUAL NUMBER OF BAND EXCEEDS DIMENSION')
            CALL ERROR$I4VAL('NX',NX)
            CALL ERROR$STOP('ATOMLIB$AESCF')
          END IF
          LOFI(NB)=LMAP(I)
          IF(TSO)THEN
            IF(LMAP(I).EQ.0)THEN
              SOFI(NB)=1
              Fofi(NB)=MIN(FTOT,2.D0)
              FTOT=FTOT-Fofi(NB)
            ELSE
              SOFI(NB)=-1
              Fofi(NB)=MIN(FTOT,REAL(2*LMAP(I),KIND=8))
              FTOT=FTOT-Fofi(NB)
              NB=NB+1
              IF(NB.GT.NX) THEN
                CALL ERROR$MSG('ACTUAL NUMBER OF BAND EXCEEDS DIMENSION')
                CALL ERROR$I4VAL('NX',NX)
                CALL ERROR$STOP('ATOMLIB$AESCF')
              END IF
              LOFI(NB)=LMAP(I)
              SOFI(NB)=1
              Fofi(NB)=MIN(FTOT,REAL(2*LMAP(I)+2,KIND=8))
              FTOT=FTOT-Fofi(NB)
            END IF
          ELSE
            Fofi(NB)=MIN(FTOT,REAL(2*(2*LMAP(I)+1),KIND=8))
            FTOT=FTOT-Fofi(NB)
            SOFI(NB)=0
          END IF
          IF(FTOT.LE.1.D-10) EXIT
        ENDDO
        IF(FTOT.GT.1.D-10) THEN
          CALL ERROR$MSG('SPECIFIED NUMBER OF BANDS IS TOO SMALL')
          CALL ERROR$MSG('#(ELECTRONS) REMAINING')
          CALL ERROR$I4VAL('NX',NX)
          CALL ERROR$R8VAL('AEZ',AEZ)
          CALL ERROR$STOP('ATOMLIB$AESCF')
        END IF
        CALL RADIAL$NUCPOT(GID,NR,AEZ,POT)
        DO L=0,MAXVAL(LOFI)
          DO ISO=-1,1
            IF(TSO.AND.ISO.EQ.0) CYCLE
            IF(.NOT.TSO.AND.ISO.NE.0) CYCLE
            ISVAR=0
            DO IB=1,NB
              IF(LOFI(IB).NE.L.OR.SOFI(IB).NE.ISO) CYCLE
              NNOFI(IB)=ISVAR
              ISVAR=ISVAR+1
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     == SELF-CONSISTENCY LOOP                                                ==
!     ==========================================================================
      tsecond=.false.
1000  continue
      XAV=0.D0
      XMAX=0.D0
      CONVG=.FALSE.
      CALL BROYDEN$NEW(NR,NBROYDENMEM,BROYDENSTEP)
      DO ITER=1,NITER
!
!       ========================================================================
!       == DETERMINE BOUND STATES FOR A GIVEN POTENTIAL AND ADD TO DENSITY    ==
!       ========================================================================
        IF(TFOCK.AND.tsecond) THEN
print*,'fockboundstate ',iter
          call radialfock$makevfock(gid,nr,nb,lofi,sofi,fofi,phi,rbox &
     &                            ,lrhox,vfock)
          vfock%scale=scale
          vfock%mux(:)=mux(:)          
          call atomlib$aescfwithHF(gid,nr,trel,nb,sofi,lofi,fofi,eofi,phi &
     &                            ,lrhox,rbox,drel,pot,mux,scale,vfock)
        else
print*,'nofockboundstate ',iter
          DO IB=1,NB
            IF(TREL) THEN
              CALL SCHROEDINGER$DREL(GID,NR,POT,EOFI(IB),DREL)
            ELSE 
              DREL(:)=0.D0 
            END IF
            G(:)=0.D0
            CALL ATOMLIB$BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RBOX,DREL,G &
     &                             ,NNOFI(IB),POT,EOFI(IB),PHI(:,IB))
          enddo
        end if
        do ib=1,nb
          IF(TREL) THEN
            G(:)=0.D0
            CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,LOFI(IB),SOFI(IB) &
     &                                       ,DREL,G,PHI(:,IB),SPHI(:,IB))
          ELSE
            SPHI(:,IB)=0.D0
          END IF
          AUX(:)=R(:)**2*(PHI(:,IB)**2+SPHI(:,IB)**2)
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
          PHI(:,IB)=PHI(:,IB)/SQRT(SVAR)
          SPHI(:,IB)=SPHI(:,IB)/SQRT(SVAR)
        ENDDO
!
!       ========================================================================
!       == ADD UP CHARGE DENSITY                                              ==
!       ========================================================================
        RHO(:)=0.D0
        DO IB=1,NB
          RHO(:)=RHO(:)+Fofi(IB)*C0LL*(PHI(:,IB)**2+SPHI(:,IB)**2)
        ENDDO
!
!       ========================================================================
!       ==  CALCULATE ENERGY                                                  ==
!       ========================================================================
        IF(TBROYDEN) POTIN=POT
        IF(CONVG) THEN
          potin=pot
          AUX(:)=0.D0
          DO IB=1,NB
            AUX(:)=AUX(:) &
       &          +(PHI(:,IB)**2+SPHI(:,IB)**2)*(EOFI(IB)-POT(:)*Y0)*Fofi(IB)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,EKIN)
          CALL ATOMLIB$BOXVOFRHO(GID,NR,RBOX,AEZ,RHO,POT,EH,EXC)
          ETOT=EKIN+EH+EXC
          pot=potin
!!$PRINT*,'ETOT ',ITER,ETOT
        END IF
!
!       ========================================================================
!       ==  EXIT IF CONVERGED                                                 ==
!       ========================================================================
        IF(CONVG) EXIT

!       ========================================================================
!       == CALCULATE OUTPUT POTENTIAL                                         ==
!       ========================================================================
        CALL ATOMLIB$BOXVOFRHO(GID,NR,RBOX,AEZ,RHO,POT,EH,EXC)
        IF(TFOCK) CALL ATOMLIB$BOXMUX(GID,NR,RBOX,RHO,MUX)
!
!       ========================================================================
!       ==  GENERATE NEXT ITERATION USING D. G. ANDERSON'S METHOD             ==
!       ========================================================================
print*,'xav ',xav,xmax
        XAV=SQRT(SUM(R**3*(POT-POTIN)**2)/SUM(R**3))
        XMAX=MAXVAL(ABS(R**2*(POT-POTIN))) 
        CALL BROYDEN$STEP(NR,POTIN,(POT-POTIN))
        POT=POTIN
        CONVG=(XMAX.LT.TOL)
      ENDDO
      CALL BROYDEN$CLEAR
      IF(.NOT.CONVG) THEN
        CALL ERROR$MSG('SELFCONSISTENCY LOOP NOT CONVERGED')
        CALL ERROR$R8VAL('AEZ',AEZ)
        CALL ERROR$R8VAL('XMAX',XMAX)
        CALL ERROR$R8VAL('TOLX',TOL)
        CALL ERROR$I4VAL('NITER',NITER)
        CALL ERROR$STOP('ATOMLIB$AESCF')
      END IF
!     == dft calculation is converged. now converge with fock term =============
      if(tfock.and.(.not.tsecond)) then
        tsecond=.true.
        goto 1000
      end if
!!$DO I=1,NB
!!$ WRITE(*,FMT='(3I4,F10.2,I5,F20.3)')I,LOFI(I),SOFI(I),F(I),NNOFI(I),EOFI(I)
!!$ENDDO
!!$PRINT*,'#ITERATIONS ',ITER
!
!     ==========================================================================
!     == REORDER STATES ACCORDING TO ENERGY                                   ==
!     ==========================================================================
      DO IB=1,NB
        IARR=MINLOC(EOFI(IB:NB))
        JB=IARR(1)+IB-1
        IF(JB.EQ.IB) CYCLE
        ISVAR=LOFI(JB)
        LOFI(IB+1:JB)=LOFI(IB:JB-1)
        LOFI(IB)=ISVAR
        ISVAR=SOFI(JB)
        SOFI(IB+1:JB)=SOFI(IB:JB-1)
        SOFI(IB)=ISVAR
        ISVAR=NNOFI(JB)
        NNOFI(IB+1:JB)=NNOFI(IB:JB-1)
        NNOFI(IB)=ISVAR
        SVAR=EOFI(JB)
        EOFI(IB+1:JB)=EOFI(IB:JB-1)
        EOFI(IB)=SVAR
        SVAR=Fofi(JB)
        Fofi(IB+1:JB)=Fofi(IB:JB-1)
        Fofi(IB)=SVAR
        AUX=PHI(:,JB)
        PHI(:,IB+1:JB)=PHI(:,IB:JB-1)
        PHI(:,IB)=AUX
        AUX=SPHI(:,JB)
        SPHI(:,IB+1:JB)=SPHI(:,IB:JB-1)
        SPHI(:,IB)=AUX
      ENDDO
!CALL ATOMLIB$FOCKSCHROEDINGER(GID,NR,POT,DREL &
!     &        ,1.D0,NB,LOFI,F,PHI,MUX,LRHOX,LOFI(IB),SOFI(IB),EOFI(IB),AUX)
!CALL ATOMLIB$FOCKBOUNDSTATE(GID,NR,LOFI(IB),SOFI,RBOX,DREL &
!     &                                 ,0.5d0,NB,LOFI,F,PHI,MUX,LRHOX &
!     &                                 ,G,NNOFI(IB),POT,EOFI(IB),AUX)
!CALL ATOMLIB_WRITEPHI('FOCKPHI.DAT',GID,NR,1,AUX)
!CALL ATOMLIB_WRITEPHI('NOFOCKPHI.DAT',GID,NR,nb,PHI(:,:))
!CALL ATOMLIB_WRITEPHI('wFOCKPHI.DAT',GID,NR,nb,PHI(:,:))
!STOP 'done'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RBOX,DREL,G,NN,POT,E,PHI)
!     **************************************************************************
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND         **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN.                     **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G                **
!     **                                                                      **
!     **  THE BOUNDARY CONDITION IS PHI(RBOX)=0                               **
!     **                                                                      **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE                 **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.         **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.         **
!     **                                                                      **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND           **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT           **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                         **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: RBOX    ! BOX RADIUS
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      INTEGER(4) ,INTENT(IN)     :: NN      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)    ,PARAMETER      :: TOL=1.D-8
      REAL(8)                    :: DER,DERO
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHIHOM(NR),PHIINHOM(NR),GHOM(NR)
      REAL(8)                    :: PHI1(NR)
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRMATCH 
      REAL(8)                    :: SVAR
      REAL(8)                    :: POT1(NR)
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20 !MAX. FACTOR IN THE WAVEFUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
      LOGICAL(4)                 :: THOM
      REAL(8)                    :: ROUT
      REAL(8)                    :: VAL1,VAL2,R1,R2
      INTEGER(4)                 :: IROUT,IRCL,IRBOX,IREND
!     **************************************************************************
                                 CALL TRACE$PUSH('ATOMLIB$BOUNDSTATE')
      CALL RADIAL$R(GID,NR,R)
!     ==  R(IRBOX) IS THE FIRST GRIDPOINT JUST IOUTSIDE THE BOX
      IRBOX=1
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GE.RBOX) EXIT
      ENDDO
!          
!     ==========================================================================
!     ==========================================================================
!     ==  ITERATE TO BISECT ENERGY WITH A NODE AT RBOX                        ==
!     ==========================================================================
!     ==========================================================================
      ISTART=1
      X0=E
      DX=1.D-2
      CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      DO I=1,NITER
        E=X0
!        E=MIN(X0,EMAX)
!
!       ========================================================================
!       ==  CUT OFF THE POTENTIAL                                             ==
!       ========================================================================
        CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT,E,IRCL,IROUT)
!       == BOUNDARY CONDITION PHI(ROUT)=0 ======================================
        IF(R(IROUT).LT.RBOX) THEN
          ROUT=R(IROUT)-1.D-5    !ENSURE THAT ROUT<R(IROUT)
        ELSE
          ROUT=RBOX
          IROUT=IRBOX
        END IF
!       ==  SET KINETIC ENERGY TO ZERO BEYOND ROUT TO AVOID AN OVERFLOW ========
!       == ATTENTION: THERE IS A STEP IN THE POTENTIAL =========================
        POT1(:)=POT(:)
        POT1(IROUT:)=E
!
!       ========================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       ========================================================================
        IDIR=1
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHI)
!       == CHECK FOR OVERFLOW
        IF(.NOT.(PHI(IROUT).GT.0.OR.PHI(IROUT).LE.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('NN',NN)
          CALL ERROR$STOP('ATOMLIB$BOUNDSTATE')
        END IF
!
!       ========================================================================
!       == ESTIMATE PHASE SHIFT                                               ==
!       ========================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,ROUT,Z0)
        Z0=Z0-REAL(NN+1)
        IF(ABS(2.D0*DX).LE.TOL) EXIT
!       ========================================================================
!       ==  BISECTION                                                         ==
!       ========================================================================
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('ATOMLIB$BOUNDSTATE')
      END IF
!
!     ==========================================================================
!     ==========================================================================
!     == INTEGRATE OUTWARD AND INVARD                                         ==
!     ==========================================================================
!     ==========================================================================
!
!     ==========================================================================
!     ==  DETERMINE MATCHING POINT                                            ==
!     ==========================================================================
      IRMATCH=IRCL
      IF(R(IRCL).GT.5.D0) THEN
        CALL RADIAL$XOFR(GID,5.D0,SVAR)
        IRMATCH=INT(SVAR)
      END IF
!
!     ==========================================================================
!     ==  INTEGRATE INWARD                                                    ==
!     ==========================================================================
      IF(IRMATCH.LT.IROUT) THEN
        THOM=MAXVAL(ABS(G(:))).EQ.0.D0
        IDIR=-1
!       ==  HOMOGENEOUS SOLUTION THAT FULFILLS THE OUTER BOUNDARY CONDITION   ==
!       ==  INTEGRATE INWARD AT THE GIVEN ENERGY WITH SLIGHTLY DIFFERENT      ==
!       ==  BOUNDARY CONDITIONS AND SUPERIMPOSE THEM SO THAT OUTER BOUNDARY   ==
!       ==  CONDITION IS EXACTLY FULFILLED                                    ==
        IREND=MIN(IROUT,NR-3)
        GHOM(:)=0.D0
        GHOM(IREND+1)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
        PHIHOM(:)=PHIHOM(:)/PHIHOM(IRMATCH)
        GHOM(:)=0.D0
        GHOM(IREND+2)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHI1)
        PHI1(:)=PHI1(:)/PHI1(IRMATCH)
!       == EXTRAPOLATE =======================
        IF(IREND.LT.IROUT) THEN
          R1=R(IREND)
          R2=R(IREND+1)
          VAL1=PHIHOM(IREND)
          VAL2=PHIHOM(IREND+1)
          PHIHOM(IREND:)=VAL1+(VAL2-VAL1)/(R2-R1)*(R(IREND:)-R1)
          VAL1=PHI1(IREND)
          VAL2=PHI1(IREND+1)
          PHI1(IREND:)=VAL1+(VAL2-VAL1)/(R2-R1)*(R(IREND:)-R1)
        END IF
!       == FULFILL OUTER BOUNDARY CONDITION ====================================
        CALL RADIAL$VALUE(GID,NR,PHIHOM,ROUT,VAL1)
        CALL RADIAL$VALUE(GID,NR,PHI1,ROUT,VAL2)
        SVAR=VAL1+VAL2
        VAL1=VAL1/SVAR
        VAL2=VAL2/SVAR
        PHIHOM(:)=VAL2*PHIHOM(:)-VAL1*PHI1(:)
        PHIHOM(:)=PHIHOM(:)/PHIHOM(IRMATCH)
!       
!       == INHOMOGENEOUS SOLUTION WITH CORRECT BOUNDARY CONDITIONS =============
        IF(.NOT.THOM) THEN     
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHIINHOM)
          CALL RADIAL$VALUE(GID,NR,PHIINHOM,ROUT,VAL1)
          CALL RADIAL$VALUE(GID,NR,PHI1,ROUT,VAL2)
          PHIINHOM(:)=PHIINHOM(:)-VAL1/VAL2*PHI1(:)
         ELSE
          PHIINHOM(:)=0.D0
        END IF
!
!       =======================================================================
!       ==  MATCH SOLUTION INSIDE AND OUTSIDE WITH VALUE                     ==
!       =======================================================================
        SVAR=(PHI(IRMATCH)-PHIINHOM(IRMATCH))/PHIHOM(IRMATCH)
        PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM(:)
        CALL RADIAL$DERIVATIVE(GID,NR,PHI,R(IRMATCH),DER)
        CALL RADIAL$DERIVATIVE(GID,NR,PHIINHOM,R(IRMATCH),DERO)
        SVAR=(DERO-DER)/PHI(IRMATCH)
        PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
      END IF
!!$!
! DO NOT NORMALIZE!!!! IT DOES NOT MAKE SENSE FOR THE INHOMOGENEOUS SOLUTION
                                 CALL TRACE$POP()
      RETURN
      END SUBROUTINE ATOMLIB$BOUNDSTATE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$phaseshiftstate(GID,NR,L,SO,DREL,G,pot,rbnd,phase &
     &                                                              ,E,PHI)
!     **************************************************************************
!     **  FINDS A state with determined phaseshift phase at radius rc         **
!     **  OF THE RADIAL SCHROEDINGER EQUATION AND                             **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G                **
!     **                                                                      **
!     **  THE ENERGY IS DETERMINED BY BISECTION ON THE                        **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.         **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.         **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: DREL(NR)! RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   ! INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: POT(NR) ! POTENTIAL
      REAL(8)    ,INTENT(IN)     :: rbnd    ! radius of phase shift
      REAL(8)    ,INTENT(IN)     :: phase   ! target phase shift
      REAL(8)    ,INTENT(INOUT)  :: E       ! ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) ! WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)    ,PARAMETER      :: TOL=1.D-8
      REAL(8)                    :: R(NR)
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRbnd
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20 !MAX. FACTOR IN THE WAVEFUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
!     **************************************************************************
                                 CALL TRACE$PUSH('ATOMLIB$BOUNDSTATE')
      CALL RADIAL$R(GID,NR,R)
      do ir=1,nr
        irbnd=ir
        if(r(ir).gt.rbnd) exit
      enddo
!          
!     ==========================================================================
!     ==========================================================================
!     ==  ITERATE TO BISECT ENERGY WITH A NODE AT RBOX                        ==
!     ==========================================================================
!     ==========================================================================
      ISTART=1
      X0=E
      DX=1.D-2
      CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      DO I=1,NITER
        E=X0
!
!       ========================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       ========================================================================
        IDIR=1
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,E,IDIR,PHI)
!       == CHECK FOR OVERFLOW
        IF(.NOT.(PHI(IRbnd).GT.0.OR.PHI(IRbnd).LE.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$STOP('ATOMLIB$BOUNDSTATE')
        END IF
!
!       ========================================================================
!       == ESTIMATE PHASE SHIFT                                               ==
!       ========================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,Rbnd,Z0)
        Z0=Z0-Phase
        IF(ABS(2.D0*DX).LE.TOL) EXIT
!       ========================================================================
!       ==  BISECTION                                                         ==
!       ========================================================================
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('ATOMLIB$BOUNDSTATE')
      END IF
                                 CALL TRACE$POP()
      RETURN
      END SUBROUTINE ATOMLIB$phaseshiftstate
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$PAWBOUNDSTATE(GID,NR,L,NN,RBOX,PSPOT,NPRO,PRO,DH,DO,G &
     &                                ,E,PHI)
!     **************************************************************************
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND         **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN.                     **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G                **
!     **                                                                      **
!     **  THE BOUNDARY CONDITION IS PHI(RBOX)=0                               **
!     **                                                                      **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE                 **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.         **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.         **
!     **                                                                      **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND           **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT           **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                         **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: NN      ! #(NODES)
      REAL(8)    ,INTENT(IN)     :: RBOX    ! BOX RADIUS
      INTEGER(4) ,INTENT(IN)     :: NPRO      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: PSPOT(NR) !POTENTIAL
      REAL(8)    ,INTENT(IN)     :: PRO(NR,NPRO)
      REAL(8)    ,INTENT(IN)     :: DH(NPRO,NPRO)
      REAL(8)    ,INTENT(IN)     :: DO(NPRO,NPRO)
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)    ,PARAMETER      :: TOL=1.D-8
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHI1(NR),PHI2(NR)
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      REAL(8)                    :: VAL1,VAL2,SVAR
      INTEGER(4)                 :: IRBOX
!     **************************************************************************
                                 CALL TRACE$PUSH('ATOMLIB$PAWBOUNDSTATE')
      CALL RADIAL$R(GID,NR,R)
!     ==  R(IRBOX) IS THE FIRST GRIDPOINT JUST IOUTSIDE THE BOX
      IRBOX=1
      DO IR=1,NR-2
        IRBOX=IR
        IF(R(IR).GE.RBOX) EXIT
      ENDDO
!          
!     ==========================================================================
!     ==========================================================================
!     ==  ITERATE TO BISECT ENERGY WITH A NODE AT RBOX                        ==
!     ==========================================================================
!     ==========================================================================
      ISTART=1
      X0=E
      DX=1.D-2
      CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      DO I=1,NITER
        E=X0
!
!       ========================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       ========================================================================
        CALL ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO,DH,DO,G,PHI)
!       == CHECK FOR OVERFLOW ==================================================
        IF(.NOT.(PHI(IRBOX+2).GT.0.OR.PHI(IRBOX+2).LE.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$STOP('ATOMLIB$PAWBOUNDSTATE')
        END IF
!
!       ========================================================================
!       == ESTIMATE PHASE SHIFT                                               ==
!       ========================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RBOX,Z0)
        Z0=Z0-REAL(NN+1,KIND=8)
        IF(ABS(2.D0*DX).LE.TOL) EXIT
        IF(Z0.GT.0.D0) THEN
          PHI1(:)=PHI(:)
        ELSE
          PHI2(:)=PHI(:)
        END IF
!       ========================================================================
!       ==  BISECTION                                                         ==
!       ========================================================================
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('ATOMLIB$PAWBOUNDSTATE')
      END IF
!     ==========================================================================
!     ==  CHOP OF TAILS WHICH MAY BE EXPONENTIALLY INCREASING                 ==
!     ==========================================================================
      CALL RADIAL$VALUE(GID,NR,PHI1,RBOX,VAL1)
      CALL RADIAL$VALUE(GID,NR,PHI2,RBOX,VAL2)
      SVAR=VAL2-VAL1
      VAL1=VAL1/SVAR
      VAL2=VAL2/SVAR
      PHI=PHI1*VAL2-PHI2*VAL1
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine atomlib$updatePAWstate(gid,nr,l,vpaw,g,rbnd,e,psi)
!     **************************************************************************
!     ** uses an approximate solution to construct a solution for             **
!     ** a PAW-type hamiltonian                                               **
!     **                                                                      **
!     ** The boundary conditions (logarithmic derivative at rbnd is preserved **
!     ** If rbnd.le.0 a fixed energy solution is searched and the outer       **
!     ** boundary condition is relaxed.                                       **
!     **                                                                      **
!     **************************************************************************
      use radialpaw_module, only: vpaw_type
      implicit none
      integer(4),intent(in)    :: gid
      integer(4),intent(in)    :: nr
      integer(4),intent(in)    :: l
      type(vpaw_type),intent(in) :: vpaw
      real(8)   ,intent(in)    :: g(nr)
      real(8)   ,intent(in)    :: rbnd
      real(8)   ,intent(inout) :: e
      real(8)   ,intent(inout) :: psi(nr)
      real(8)                  :: phase
      real(8)                  :: drel(nr)
      real(8)                  :: dpsi(nr)
      real(8)                  :: opsi(nr)
      real(8)                  :: pot(nr)
      real(8)                  :: phidot(nr)
      real(8)                  :: phiprime(nr)
      real(8)                  :: valdot,derdot,valprime,derprime,val0,der0
      real(8)                  :: pi,y0
      real(8)                  :: r(nr)
      real(8)                  :: aux(nr),aux1(nr),aux2(nr)
      real(8)                  :: g1(nr)
      real(8)                  :: val,der,svar
      integer(4)               :: ir
      integer(4),parameter     :: niter=20
      integer(4)               :: iter
      real(8)   ,parameter     :: tol=1.d-8
      logical(4)               :: convg
      REAL(8)   ,PARAMETER     :: XMAX=1.D+15 !MAX. FACTOR IN THE WAVEFUNCTION
      integer(4)               :: irbnd
      real(8)                  :: norm,e1
      logical(4)               :: tpr=.true.
      logical(4)               :: tfixe
      integer(4)               :: nn
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      call radial$r(gid,nr,r)
      tfixe=rbnd.le.0.d0
      drel(:)=0.d0
!
!     :: remark: schroedinger$phaseshift takes value and derivative from a    ::
!     :: two-point formula, and therefore is not fully consistent with        ::
!     :: radial$value and radial$derivative                                   ::
      CALL SCHROEDINGER$PHASESHIFT(GID,NR,Psi,Rbnd,phase)
print*,'phase',phase
!
      do ir=1,nr
        irbnd=ir
        if(tfixe.and.r(ir).gt.3.d0) exit
        if(.not.tfixe.and.r(ir).gt.rbnd) exit
      enddo
!
!     ==========================================================================
!     ==  construct an local, approximate potential
!     ==========================================================================
      aux(:)=1.d0
      call radialpaw$vpsi(GID,NR,vpaw,l,aux,aux1,opsi)
      pot(:)=aux1-e*opsi(:)+e*aux(:)
pot=vpaw%pspot
!
!     ==========================================================================
!     ==  now start loop                                                      ==
!     ==========================================================================
      do iter=1,niter
!
!       ========================================================================
!       ==  construct |dpsi>=de/dpsi                                          ==
!       ========================================================================
!       == kinetic energy of the wave function (nonrelativistic)
        call radial$verletd2(gid,nr,r*psi,aux1)
        aux(:)=-r(:)*aux(:)+real(l*(l+1),kind=8)*psi(:)
!       == remove factor 2*r**2 ================================================
        aux(:)=aux(:)/(2.d0*r(:)**2)
!       == correct glitches ====================================================
        aux(1:2)=0.d0
!       == potential energy ====================================================
        call radialpaw$vpsi(GID,NR,vpaw,l,psi,aux1,opsi)
        aux(:)=aux(:)+aux1(:)
!       == estimate energy =====================================================
        if(.not.tfixe) then
          aux1(:)=r(:)**2*psi(:)*(aux(:)-g(:))
          call radial$integrate(gid,nr,aux1,aux2)
          call radial$value(gid,nr,aux2,rbnd,e)
          aux1(:)=r(:)**2*psi(:)*opsi(:)
          call radial$integrate(gid,nr,aux1,aux2)
          call radial$value(gid,nr,aux2,rbnd,norm)
          e=e/norm
!         == remove energy term ================================================
          aux(:)=aux(:)-e*opsi(:)
        else
          aux1(:)=r(:)**2*psi(:)*opsi(:)
          call radial$integrate(gid,nr,aux1,aux2)
          call radial$value(gid,nr,aux2,3.d0,norm)
        end if
!       == subtract inhomogeneity ==============================================
        aux(:)=aux(:)-g(:)
!       == map into dpsi
        dpsi(:)=-aux(:)
!
!       ========================================================================
!       ==  Determine phiprime-phidot*eprime                                  ==
!       ========================================================================
        g1=dpsi(:)
        g1(1:2)=0.d0
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,0,G1,L,E,1,phiprime)
        if(.not.tfixe) then
          g1=psi(:)
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,0,G1,L,E,1,phidot)
          call radial$value(gid,nr,psi,rbnd,val0)
          call radial$derivative(gid,nr,PSI,rbnd,der0)
          call radial$value(gid,nr,phiprime,rbnd,valprime)
          call radial$derivative(gid,nr,phiprime,rbnd,derprime)
          call radial$value(gid,nr,phidot,rbnd,valdot)
          call radial$derivative(gid,nr,phidot,rbnd,derdot)
          call SCHROEDINGER$resolvePHASESHIFT(PHASE,val,der,nn)
          val0=val*der0-der*val0
          valprime=val*derprime-der*valprime
          valdot=val*derdot-der*valdot
          svar=(val0+valprime)/valdot  !delta epsilon
          DPSI(:)=phiprime(:)-phidot(:)*svar
        else
          DPSI(:)=phiprime(:)
        end if
!   
!       ========================================================================
!       ==  check convergence                                                 ==
!       ========================================================================
        val=maxval(abs(dpsi(:irbnd-3)))/sqrt(norm)
        if(tpr) print*,' dev=',val
        convg=(val.lt.tol)
        if(convg) exit
!   
!       ========================================================================
!       ==  propagate                                                         ==
!       ========================================================================
        psi(:)=psi(:)+dpsi(:)
!
      ENDDO    ! end of iteration
!
      IF(.NOT.CONVG) THEN
        call error$msg('LOOP FOR RADIAL paw NOT CONVERGED')
        call error$stop('atomlib$updatepawstate')
      end if
      RETURN
      end subroutine atomlib$updatePAWstate
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$BOXVOFRHO(GID,NR,RAD,AEZ,RHO,POT,EH,EXC)
!     ******************************************************************
!     **                                                              **
!     **  ELECTROSTATIC AND EXCHANGE-CORRELATION POTENTIAL            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8),   INTENT(IN) :: AEZ
      REAL(8),   INTENT(IN) :: RAD
      REAL(8),   INTENT(IN) :: RHO(NR)
      REAL(8),   INTENT(OUT):: POT(NR)
      REAL(8),   INTENT(OUT):: EH
      REAL(8),   INTENT(OUT):: EXC
      REAL(8),   PARAMETER  :: RHOMIN=1.D-2
      REAL(8)               :: POTH(NR)
      REAL(8)               :: POTXC(NR)
      REAL(8)               :: PI
      REAL(8)               :: FOURPI
      REAL(8)               :: Y0
      REAL(8)               :: RHO1(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: EDEN(NR)
      REAL(8)               :: GRHO(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VGXC,VXC,EXC1,RH,GRHO2
      REAL(8)               :: DUMMY1,DUMMY2,DUMMY3
      REAL(8)               :: SVAR
      INTEGER(4)            :: IR,IRBOX
      REAL(8)               :: Q
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GE.RAD) EXIT
      ENDDO
!
!     ==========================================================================
!     ==  TOTAL CHARGE                                                        ==
!     ==========================================================================
      AUX(:)=4.D0*PI*RHO(:)*R(:)**2*Y0
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,Q)
      Q=Q-AEZ
!
!     ==========================================================================
!     ==  TOTAL POTENTIAL                                                     ==
!     ==========================================================================
      RHO1=RHO
      IR=MIN(IRBOX+3,NR)
      RHO1(IR:)=0.D0
      CALL RADIAL$NUCPOT(GID,NR,AEZ,POTH)
      EDEN(:)=0.5D0*RHO1(:)*POTH(:)
      CALL RADIAL$POISSON(GID,NR,0,RHO1,AUX)
      POTH(:)=POTH(:)+AUX(:)
      CALL RADIAL$VALUE(GID,NR,POTH,RAD,SVAR)
      SVAR=Q/RAD/Y0-SVAR
      POTH(1:IRBOX-1)=POTH(1:IRBOX-1)+SVAR
      POTH(IRBOX:NR)=Q/R(IRBOX:NR)/Y0
      EDEN(:)=EDEN(:)+0.5D0*RHO1(:)*POTH(:)
!
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,EDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,EH)
!
      CALL RADIAL$DERIVE(GID,NR,RHO,GRHO)
      DO IR=1,NR
        RH=RHO(IR)*Y0
        GRHO2=(Y0*GRHO(IR))**2
        CALL DFT(RH,0.D0,GRHO2,0.D0,0.D0,EXC1,VXC,DUMMY1,VGXC,DUMMY2,DUMMY3)
!write(6,fmt='("XCTEST",i5,10f30.10)')ir,r(ir),rh,grho2,exc1,vxc,vgxc
        EDEN(IR)=4.D0*PI*EXC1   ! ANGULAR INTEGRATION ALREADY INCLUDED
        POTXC(IR)=VXC/Y0
        GRHO(IR)=VGXC*2.D0*GRHO(IR)
      ENDDO
      CALL RADIAL$DERIVE(GID,NR,GRHO(:),AUX)
      POTXC(:)=POTXC(:)-AUX
      IF(R(1).GT.1.D-8) THEN
         POTXC(:)=POTXC(:)-2.D0/R(:)*GRHO(:)
      ELSE
        POTXC(2:)=POTXC(2:)-2.D0/R(2:)*GRHO(2:)
!        POTXC(1)=POTXC(1)-2.D0/R(2)*GRHO(2)
        POTXC(1)=POTXC(2)
      END IF
!
!     == EXCHANGE CORRELATION POTENTIAL IS SET TO ZERO OUTSIDE THE BOX. 
!     == OTHERWISE IT WILL HAVE A CONSTANT VALUE IN THE ZERO-DENSITY REGION
      POTXC(IRBOX:)=0.D0
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,EDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,EXC)
!
!     ==========================================================================
!     ==  CUT OF POTENTIAL FOR LOW DENSITIES                                  ==
!     ==========================================================================
      DO IR=1,NR
        IF(RHO(IR)*Y0.LT.1.D-6) POTXC(IR)=0.D0
      ENDDO
!
!     ==========================================================================
!     ==  CUT OF POTENTIAL FOTR LOW DENSITIES                                 ==
!     ==========================================================================
      POT=POTH+POTXC
!do ir=1,nr
!write(6,fmt='("potTEST",i5,10f30.10)')ir,rho(ir),pot(ir),poth(ir),potxc(ir)
!enddo
!stop
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$BOXMUX(GID,NR,RAD,RHO,MUX)
!     ******************************************************************
!     **                                                              **
!     **  ELECTROSTATIC AND EXCHANGE-CORRELATION POTENTIAL            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8),   INTENT(IN) :: RAD
      REAL(8),   INTENT(IN) :: RHO(NR)
      REAL(8),   INTENT(OUT):: MUX(NR)
      REAL(8),   PARAMETER  :: RHOMIN=1.D-2
      REAL(8)               :: PI
      REAL(8)               :: FOURPI
      REAL(8)               :: Y0
      REAL(8)               :: AUX(NR)
      REAL(8)               :: GRHO(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VGXC,VXC,EXC1,RH,GRHO2
      REAL(8)               :: DUMMY1,DUMMY2,DUMMY3
      INTEGER(4)            :: IR,IRBOX
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GE.RAD) EXIT
      ENDDO
      CALL DFT$SETL4('XCONLY',.TRUE.)
!
!     ==========================================================================
!     ==  TOTAL POTENTIAL                                                     ==
!     ==========================================================================
      CALL RADIAL$DERIVE(GID,NR,RHO,GRHO)
      DO IR=1,NR
        RH=RHO(IR)*Y0
        GRHO2=(Y0*GRHO(IR))**2
        CALL DFT(RH,0.D0,GRHO2,0.D0,0.D0,EXC1,VXC,DUMMY1,VGXC,DUMMY2,DUMMY3)
        MUX(IR)=VXC/Y0
        GRHO(IR)=VGXC*2.D0*GRHO(IR)
      ENDDO
      CALL RADIAL$DERIVE(GID,NR,GRHO(:),AUX)
      MUX(:)=MUX(:)-AUX
      IF(R(1).GT.1.D-8) THEN
         MUX(:)=MUX(:)-2.D0/R(:)*GRHO(:)
      ELSE
        MUX(2:)=MUX(2:)-2.D0/R(2:)*GRHO(2:)
        MUX(1)=MUX(2)
      END IF
!
!     == EXCHANGE CORRELATION POTENTIAL IS SET TO ZERO OUTSIDE THE BOX. 
!     == OTHERWISE IT WILL HAVE A CONSTANT VALUE IN THE ZERO-DENSITY REGION
      MUX(IRBOX:)=0.D0
!
!     ==========================================================================
!     ==  CUT OF POTENTIAL FOR LOW DENSITIES                                  ==
!     ==========================================================================
      DO IR=1,NR
        IF(RHO(IR)*Y0.LT.1.D-6) MUX(IR)=0.D0
      ENDDO
      CALL DFT$SETL4('XCONLY',.FALSE.)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$NODELESS(GID,NR,RBOX,POT,NB,LOFI,EOFI,UOFI,TUOFI)
!     **************************************************************************
!     ** CALCULATES THE SEQUENCE OF NODELESS WAVE FUNCTIONS                   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: GID
      INTEGER(4),INTENT(IN)   :: NR
      REAL(8)   ,INTENT(IN)   :: RBOX
      REAL(8)   ,INTENT(IN)   :: POT(NR)
      INTEGER(4),INTENT(IN)   :: NB
      INTEGER(4),INTENT(IN)   :: LOFI(NB)
      REAL(8)   ,INTENT(INOUT):: EOFI(NB)
      REAL(8)   ,INTENT(OUT)  :: UOFI(NR,NB)
      REAL(8)   ,INTENT(OUT)  :: TUOFI(NR,NB)
      INTEGER(4)              :: L
      REAL(8)                 :: E
      REAL(8)                 :: G(NR)
      REAL(8)                 :: DREL(NR)
      INTEGER(4)              :: NN
      INTEGER(4)              :: SO
      INTEGER(4)              :: IB,I
      REAL(8)                 :: PI,Y0
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      DO IB=1,NB
        L=LOFI(IB)
        E=EOFI(IB)
        G(:)=0.D0
        DO I=IB-1,1,-1
          IF(LOFI(I).EQ.L) THEN
            G=UOFI(:,I)
            EXIT
          END IF
        ENDDO
        SO=0
        DREL(:)=0.D0
        NN=0
        CALL ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RBOX,DREL,G,NN,POT,E,UOFI(:,IB))
        TUOFI(:,IB)=G(:)+(E-POT(:)*Y0)*UOFI(:,IB)
        EOFI(IB)=E
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO,DH,DO,G,PHI)
!     **************************************************************************
!     **                                                                      **
!     **  SOLVES THE RADIAL PAW -SCHROEDINGER EQUATION.                       **
!     **    (T+VTILDE-E+|P>(DH-E*DO<P|]|PHI>=|G>                              **
!     **  WHERE T IS THE NONRELATIVISTIC KINETIC ENERGY.                      **
!     **                                                                      **
!     **    DH=<AEPHI|T+V|AEPHI>-<PSPHI|T+VTILDE|PSPHI>                       **
!     **    DO=<AEPHI|AEPHI>-<PSPHI|PSPHI>                                    **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID     !GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR      ! #(RADIAL GRID POINTS
      INTEGER(4)  ,INTENT(IN) :: L       ! ANGULAR MOMENTUM
      REAL(8)     ,INTENT(IN) :: E       ! ENERGY
      REAL(8)     ,INTENT(IN) :: PSPOT(NR)  !VTILDE=PSPOT*Y0
      INTEGER(4)  ,INTENT(IN) :: NPRO       ! #(PROJECTOR FUNCTIONS)
      REAL(8)     ,INTENT(IN) :: PRO(NR,NPRO) ! PROJECTOR FUNCTIONS
      REAL(8)     ,INTENT(IN) :: DH(NPRO,NPRO)     
      REAL(8)     ,INTENT(IN) :: DO(NPRO,NPRO)     
      REAL(8)     ,INTENT(IN) :: G(NR)        !INHOMOGENEITY
      REAL(8)     ,INTENT(OUT):: PHI(NR)      !PAW RADIAL FUNCTION
      REAL(8)                 :: U(NR)
      REAL(8)                 :: V(NR,NPRO)
      REAL(8)                 :: AMAT(NPRO,NPRO)
      REAL(8)                 :: BMAT(NPRO,NPRO)
      REAL(8)                 :: BMATINV(NPRO,NPRO)
      REAL(8)                 :: CMAT(NPRO,NPRO)
      REAL(8)                 :: CVEC(NPRO)
      REAL(8)                 :: DVEC(NPRO)
      INTEGER(4)              :: I1,I2,I3
      REAL(8)                 :: AUX(NR)
      REAL(8)                 :: R(NR)
      REAL(8)                 :: DREL(NR)
      INTEGER(4)              :: SO
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  -1/2NABLA^2+POT-E|U>=|G>                                            ==
!     ==========================================================================
      SO=0
      DREL(:)=0.D0
      CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,G,L,E,1,U)
!
!     ==========================================================================
!     ==  -1/2NABLA^2+POT-E|V>=|PRO>                                          ==
!     ==========================================================================
      DO I1=1,NPRO
        CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,PRO(:,I1),L,E,1,V(:,I1))
      ENDDO
!!$PRINT*,'E ',E,SO,L,NR,GID
!!$OPEN(100,FILE='XXX.DAT')
!!$ DO IR=1,NR
!!$    WRITE(100,FMT='(22F30.10)')R(IR),U(IR),V(IR,:),PRO(IR,:)
!!$ ENDDO
!!$CLOSE(10)
!!$STOP 'IN PAWDER'
!
!     ==========================================================================
!     ==  AMAT=<PRO|V>  CVEC=<PRO|U>                                          ==
!     ==========================================================================
      DO I1=1,NPRO
        DO I2=1,NPRO
          AUX(:)=PRO(:,I1)*V(:,I2)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,AMAT(I1,I2))
        ENDDO
      ENDDO
      DO I1=1,NPRO
        AUX(:)=PRO(:,I1)*U(:)*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,CVEC(I1))
      ENDDO  
!
!     ==========================================================================
!     ==  BMAT=1+(DATH-EDO)<PRO|V>                                            ==
!     ==========================================================================
!     BMAT(:,:)=MATMUL(DH(:,:)-E*DO(:,:),AMAT(:,:))
!     DO I1=1,NPRO
!       BMAT(I1,I1)=BMAT(I1,I1)+1.D0
!     ENDDO
!
      DO I1=1,NPRO
        DO I2=1,NPRO
          BMAT(I1,I2)=0.D0
          DO I3=1,NPRO
            BMAT(I1,I2)=BMAT(I1,I2)+(DH(I1,I3)-E*DO(I1,I3))*AMAT(I3,I2)     
          ENDDO
        ENDDO
        BMAT(I1,I1)=BMAT(I1,I1)+1.D0
      ENDDO
!
!     ==========================================================================
!     ==  BMAT = BMAT^-1 = [1+(DATH-EDO)<PRO|V>]^-1                           ==
!     ==========================================================================
      IF(NPRO.EQ.0) THEN
        CALL ERROR$STOP('NPRO=0 NOT ALLOWED')
      END IF
      IF(NPRO.EQ.1) THEN 
        BMAT(1,1)=1.D0/BMAT(1,1)
      ELSE 
        CALL LIB$INVERTR8(NPRO,BMAT,BMATINV)
        BMAT=BMATINV
      END IF
!
!     ==========================================================================
!     ==  CMAT = -BMAT*(DATH-EDO)                                             ==
!     ==       = -[1+(DATH-EDO)*<PRO|V>]^-1 (DATH-EDO)                        ==
!     ==========================================================================
!     CMAT(:,:)=MATMUL(BMAT(:,:),DH(:,:)-E*DO(:,:))
      DO I1=1,NPRO
        DO I2=1,NPRO
          CMAT(I1,I2)=0.D0
          DO I3=1,NPRO
            CMAT(I1,I2)=CMAT(I1,I2)-BMAT(I1,I3)*(DH(I3,I2)-E*DO(I3,I2))
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  DVEC = CMAT*CVEC                                                    ==
!     ==       = -[1+(DATH-EDO)*<PRO|V>]^-1 (DATH-EDO) <PRO|U>                ==
!     ==========================================================================
!     DVEC(:)=MATMUL(CMAT(:,:),CVEC(:))
      DO I1=1,NPRO
        DVEC(I1)=0.D0
        DO I2=1,NPRO
          DVEC(I1)=DVEC(I1)+CMAT(I1,I2)*CVEC(I2)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  |PHI> = |U>+|V>DVEC                                                 ==
!     ==  = [1-|V>[1+(DATH-EDO)*<PRO|V>]^-1(DATH-EDO)<PRO|] |U>               ==
!     ==========================================================================
!     PHI(:)=U(:)+MATMUL(V(:,:),DVEC(:))
      PHI(:)=U(:)
      DO I1=1,NPRO
        PHI(:)=PHI(:)+V(:,I1)*DVEC(I1)
      ENDDO
      RETURN
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine atomlib$aescfwithHF(gid,nr,trel,nb,sofi,lofi,fofi,eofi,psi &
     &                              ,lrhox,rbox,drel,pot,mux,scale,vfock)
!     **************************************************************************
!     **************************************************************************
      use radialfock_module
      implicit none
      integer(4),intent(in)    :: gid
      integer(4),intent(in)    :: nr
      integer(4),intent(in)    :: nb
      logical(4),intent(in)    :: trel
      integer(4),intent(in)    :: sofi(nb)
      integer(4),intent(in)    :: lofi(nb)
      real(8)   ,intent(in)    :: fofi(nb)
      integer(4),intent(in)    :: lrhox
      real(8)   ,intent(in)    :: rbox
      real(8)   ,intent(in)    :: mux(nr)
      real(8)   ,intent(inout) :: scale
      real(8)   ,intent(inout) :: eofi(nb)
      real(8)   ,intent(inout) :: drel(nr)
      real(8)   ,intent(inout) :: pot(nr)
      real(8)   ,intent(inout) :: psi(nr,nb)
      type(vfock_type),intent(inout) :: vfock
      real(8)                  :: dpsi(nr,nb)
      real(8)                  :: drrpsi(nr,nb)
      real(8)                  :: pi,y0
      real(8)                  :: sofactor
      real(8)                  :: r(nr)
      real(8)                  :: aux(nr),aux1(nr),aux2(nr)
      real(8)                  :: g(nr)
      real(8)                  :: rdprime(nr)
      real(8)                  :: val,der,val1,val2
      real(8)                  :: etot
      integer(4)               :: l,ib,ib1,ib2,i,ir,iso
      integer(4)               :: irbox
      integer(4)               :: nn
      integer(4),parameter     :: niter=20
      integer(4)               :: iter
      real(8)   ,parameter     :: tol=1.d-3
      logical(4)               :: convg
      real(8)                  :: e
      real(8)                  :: pot1(nr)
      REAL(8)   ,PARAMETER     :: XMAX=1.D+15 !MAX. FACTOR IN THE WAVEFUNCTION
      integer(4)               :: ircl,irout
      real(8)                  :: rout
      logical(4)               :: tpr=.false.
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      call radial$r(gid,nr,r)
!
!     == index of the first gridpoint outside rbox =============================
      do ir=1,nr
        irbox=ir
        if(r(ir).ge.rbox) exit
      enddo

POT(IRBOX+2:)=0.D0
!
      IF(TPR) THEN
        PRINT*,'SCALE ',SCALE,'TREL ',TREL,' RBOX ',RBOX,R(IRBOX)
        PRINT*,'EOFI ',0,EOFI
      END IF
!
!     ==========================================================================
!     ==  RECREATE WAVE FUNCTIONS FOR TESTING                                 ==
!     ==========================================================================
      DO IB=1,NB
        NN=0
        DO I=1,IB-1
          IF(LOFI(I).EQ.LOFI(IB)) NN=NN+1
        ENDDO
        IF(TREL) THEN
          CALL SCHROEDINGER$DREL(GID,NR,POT,EOFI(IB),DREL)
        ELSE 
          DREL(:)=0.D0 
        END IF
        G(:)=0.D0
        CALL ATOMLIB$BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RBOX,DREL,G,NN &
     &                           ,POT,EOFI(IB),PSI(:,IB))
        AUX(:)=R(:)**2*PSI(:,IB)**2
        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
        CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
        PSI(:,IB)=PSI(:,IB)/SQRT(VAL)
      ENDDO
        
      if(tpr) then
        CALL ATOMLIB_WRITEPHI('PSI_I.DAT',GID,NR,NB,PSI)
        PRINT*,'EOFI ',-1,EOFI
      end if
!
!     ========================================================================
!     ==  estimate first-order change of the energy                         ==
!     ========================================================================
      if(tpr) then
        vfock%scale=0.d0
        call ATOMLIB$etotwithfoCK(GID,NR,rbox,lrhox,trel,pot,vfock &
    &                      ,NB,LOFI,sofi,FOFI,PSI,val1)
        vfock%scale=scale
        print*,'energy without HF                      :',val1
        call ATOMLIB$etotwithfoCK(GID,NR,rbox,lrhox,trel,pot,vfock &
    &                      ,NB,LOFI,sofi,FOFI,PSI,val2)
        print*,'energy with scaled (HF-mu_X) correction:' ,val2
        print*,'scaled (HF-mu_X) correction            :' ,val2-val1
        vfock%scale=1.d0
        call ATOMLIB$etotwithfoCK(GID,NR,rbox,lrhox,trel,pot,vfock &
    &                      ,NB,LOFI,sofi,FOFI,PSI,val2)
        vfock%scale=scale
        print*,'100% (HF-Mu_X) correction              :',val2-val1
        aux(:)=vfock%mux
        vfock%mux=0.d0
        call ATOMLIB$etotwithfoCK(GID,NR,rbox,lrhox,trel,pot,vfock &
    &                      ,NB,LOFI,sofi,FOFI,PSI,val2)
        vfock%mux=aux
        print*,'100% HF correction                     :' ,val2-val1
      end if
!
!     ==========================================================================
!     ==  now start loop                                                      ==
!     ==========================================================================
      do iter=1,niter
!   
!       ========================================================================
!       ==  enforce boundary conditions at rbox                               ==
!       ========================================================================
        psi(irbox+2:,:)=0.d0
!
!       ========================================================================
!       ==  vfock|psi> and partial_r*r|\psi>                                  ==
!       ========================================================================
        call radialfock$makevfock(gid,nr,nb,lofi,sofi,fofi,psi &
     &                               ,rbox,lrhox,vfock)
        vfock%scale=scale
        vfock%mux=mux
!
!       ========================================================================
!       ==  calculate total energy (this a fake quantity)                     ==
!       ========================================================================
        if(tpr) then
          do ib=1,nb
            aux(:)=r(:)*psi(:,ib)
            val=aux(irbox)
            der=(aux(irbox+1)-aux(irbox-1))/(r(irbox+1)-r(irbox-1))
            aux(irbox:)=val+(r(irbox:)-r(irbox))*der
            call radial$verletd1(gid,nr,aux,aux1)
            aux1(1:2)=aux1(3)  
            drrpsi(:,ib)=aux1(:)
          enddo
          aux(:)=0.d0
          do ib=1,nb
            l=lofi(ib)
            if(trel) then
              CALL SCHROEDINGER$DREL(GID,NR,POT,EOFI(IB),DREL)
              CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
            else
              drel=0.d0
              rdprime=0.d0
            end if
            if(sofi(ib).eq.0) then
              sofactor=0.d0
            else if(sofi(ib).eq.1) then
              sofactor=real(l,kind=8)
            else if(sofi(ib).eq.-1) then
              sofactor=real(-l-1,kind=8)
            else
              call error$msg('illegal value for sofi ')
              call error$stop('atomlib$aescfwithHF')
            end if
!           == kinetic energy ==================================================
            aux1(:)=0.5d0*(1.d0+drel) &
         &               *(drrpsi(:,ib)**2+real(l*(l+1),kind=8)*psi(:,ib)**2)
            aux1(:)=aux1(:)+0.5d0*(sofactor+1.d0)*rdprime*r(:)*psi(:,ib)**2
!           == potential energy ================================================
            aux1(:)=aux1(:)+r(:)**2*pot(:)*y0*psi(:,ib)**2
!           == Hartree-Fock correction =========================================
            call radialfock$vpsi(GID,NR,vfock,lofi(ib),psi(:,ib),aux2)
            aux2(:)=psi(:,ib)*(0.5d0*aux2-0.5d0*vfock%mux(:)*y0*psi(:,ib))
            aux1(:)=aux1(:)+r(:)**2*aux2(:)
!           == add to sum over states ==========================================
            aux(:)=aux(:)+fofi(ib)*aux1(:)
          enddo
          call radial$integrate(gid,nr,aux,aux1)
          call radial$value(gid,nr,aux1,rbox,etot)
        end if
!
!       ========================================================================
!       ==  construct |dpsi>=de/dpsi                                          ==
!       ========================================================================
        do ib=1,nb
          l=lofi(ib)
          if(trel) then
            CALL SCHROEDINGER$DREL(GID,NR,POT,EOFI(IB),DREL)
            CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
          else
            drel=0.d0
            rdprime=0.d0
          end if
          if(sofi(ib).eq.0) then
            sofactor=0.d0
          else if(sofi(ib).eq.1) then
            sofactor=real(l,kind=8)
          else if(sofi(ib).eq.-1) then
            sofactor=real(-l-1,kind=8)
          else
            call error$msg('illegal value for sofi ')
            call error$stop('atomlib$aescfwithHF')
          end if
!         == kinetic energy of the wave function (nonrelativistic)
          call radial$verletd2(gid,nr,r*psi(:,ib),aux2)
          aux(:)=r(:)*aux2(:)-real(l*(l+1),kind=8)*psi(:,ib)
          aux(:)=-(1.d0+drel)*aux(:)
          call radial$verletd1(gid,nr,psi(:,ib),aux2)
          aux(:)=aux(:)-rdprime*(r(:)**2*aux2(:)-sofactor*r(:)*psi(:,ib))
!         == potential energy ==================================================
          aux(:)=aux(:)+2.d0*r(:)**2*pot(:)*y0*psi(:,ib)
!         == Hartree-Fock correction ===========================================
          call radialfock$vpsi(GID,NR,vfock,lofi(ib),psi(:,ib),aux2)
          aux(:)=aux(:)+2.d0*r(:)**2*aux2(:)
!         == correct glitches ==================================================
          aux(1:2)=0.d0
!         == estimate energy ===================================================
          aux1(:)=0.5d0*aux(:)*psi(:,ib)
          call radial$integrate(gid,nr,aux1,aux2)
          call radial$value(gid,nr,aux2,rbox,eofi(ib))
!         == remove energy term ================================================
          aux(:)=aux(:)-2.d0*r(:)**2*eofi(ib)*psi(:,ib)
!         == weight with occupation ============================================
          dpsi(:,ib)=aux(:)/(2.d0*r(:)**2)
        enddo
!
!       ========================================================================
!       ==  Determine phiprime-phidot*eprime                                  ==
!       ========================================================================
        do ib=1,nb
          l=lofi(ib)
          iso=sofi(ib)
          e=eofi(ib)
          CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT,E,IRCL,IROUT)
!         == BOUNDARY CONDITION PHI(ROUT)=0 ====================================
          IF(R(IROUT).lt.RBOX) THEN
            rout=r(irout)
          else
            rout=rbox
            IROUT=IRBOX
          END IF
!         ==  SET KINETIC ENERGY TO ZERO BEYOND ROUT TO AVOID AN OVERFLOW ======
!         == ATTENTION: THERE IS A STEP IN THE POTENTIAL =======================
          POT1(:)=POT(:)
          POT1(IROUT:)=E
          IF(TREL) THEN
            CALL SCHROEDINGER$DREL(GID,NR,POT1,E,DREL)
          ELSE 
            DREL(:)=0.D0 
          END IF
          g=dpsi(:,ib)
          g(1:2)=0.d0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,iso,G,L,E,1,aux1)
          g=psi(:,ib)
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,iso,G,L,E,1,aux2)
          call radial$value(gid,nr,aux1,rout,val1)
          call radial$value(gid,nr,aux2,rout,val2)
          dpsi(:,ib)=aux1(:)-aux2(:)*val1/val2
          dpsi(irout+3:,ib)=0.d0
          if(rout.lt.rbox) dpsi(irout:,ib)=0.d0
          call radial$value(gid,nr,dpsi(:,ib),rbox,val)
        enddo
!   
!       ========================================================================
!       ==  check convergence                                                 ==
!       ========================================================================
        val=maxval(abs(dpsi(:irbox-3,:)))
        if(tpr) print*,'etot=',etot,' dev=',val
        convg=(val.lt.tol)
        if(convg) exit
!   
!       ========================================================================
!       ==  propagate                                                         ==
!       ========================================================================
        psi(:,:)=psi(:,:)-dpsi(:,:)
!   
!       ========================================================================
!       ==  orthogonalize                                                     ==
!       ========================================================================
        do ib1=1,nb
          do ib2=1,ib1-1
            if(lofi(ib1).ne.lofi(ib2)) cycle
            if(sofi(ib1).ne.sofi(ib2)) cycle
            aux(:)=r(:)**2*psi(:,ib1)*psi(:,ib2)
            call radial$integrate(gid,nr,aux,aux1)
            call radial$value(gid,nr,aux1,rbox,val)
            psi(:,ib1)=psi(:,ib1)-psi(:,ib2)*val
          enddo
          aux(:)=r(:)**2*psi(:,ib1)**2
          call radial$integrate(gid,nr,aux,aux1)
          call radial$value(gid,nr,aux1,rbox,val)
          psi(:,ib1)=psi(:,ib1)/sqrt(val)
        enddo
      ENDDO    ! end of iteration
!
      IF(.NOT.CONVG) THEN
        call error$msg('LOOP FOR RADIAL HARTREE FOCK NOT CONVERGED')
        call error$stop('atomlib$aescfwithHF')
      end if
      if(tpr) PRINT*,ITER,'ITERATIONSIN LOOP FOR RADIAL HARTREE FOCK'
      if(tpr) then
        PRINT*,'EOFI ',ITER,EOFI
        CALL ATOMLIB_WRITEPHI('PSI_F.DAT',GID,NR,NB,PSI)
      end if
      RETURN
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine atomlib$boundstatewithHF(gid,nr,trel,rbox,pot &
     &                             ,vfock,l,so,nn,g,e,psi,tpsi)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      use radialfock_module, only : vfock_type
      implicit none
      integer(4),intent(in)    :: gid
      integer(4),intent(in)    :: nr
      logical(4),intent(in)    :: trel
      real(8)   ,intent(in)    :: rbox
      type(vfock_type),intent(in) :: vfock
      real(8)   ,intent(in)    :: pot(nr)
      integer(4),intent(in)    :: nn
      integer(4),intent(in)    :: l
      integer(4),intent(in)    :: so
      real(8)   ,intent(in)    :: g(nr)
      real(8)   ,intent(out)   :: psi(nr)
      real(8)   ,intent(out)   :: e
      real(8)   ,intent(out)   :: tpsi(nr)
      real(8)                  :: drel(nr)
      real(8)                  :: dpsi(nr)
      real(8)                  :: vfockpsi(nr)
      real(8)                  :: pi,y0
      real(8)                  :: sofactor
      real(8)                  :: r(nr)
      real(8)                  :: aux(nr),aux1(nr),aux2(nr)
      real(8)                  :: g1(nr)
      real(8)                  :: rdprime(nr)
      real(8)                  :: val,val1,val2
      integer(4)               :: ir,iso
      integer(4)               :: irbox
      integer(4),parameter     :: niter=20
      integer(4)               :: iter
      real(8)   ,parameter     :: tol=1.d-8
      logical(4)               :: convg
      real(8)                  :: pot1(nr)
      REAL(8)   ,PARAMETER     :: XMAX=1.D+15 !MAX. FACTOR IN THE WAVEFUNCTION
      integer(4)               :: ircl,irout
      real(8)                  :: rout
      real(8)                  :: norm
      logical(4)               :: tpr=.false.
      logical(4)               :: thom
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      call radial$r(gid,nr,r)
      thom=maxval(abs(g)).lt.1.d-8
!
!     == index of the first gridpoint outside rbox =============================
      do ir=1,nr
        irbox=ir
        if(r(ir).ge.rbox) exit
      enddo
!
      IF(TPR) THEN
        print*,'==============================================================='
        PRINT*,'l=',l,' so=',so,' nn=',nn,' thom ',thom
        PRINT*,'SCALE ',vfock%SCALE,'TREL ',TREL,' RBOX ',RBOX,R(IRBOX)
        PRINT*,'EOFI ',-1,E
      END IF
!
!     ==========================================================================
!     ==  RECREATE WAVE FUNCTION                                              ==
!     ==========================================================================
      e=0.d0 
      IF(TREL) THEN
        CALL SCHROEDINGER$DREL(GID,NR,POT,E,DREL)
      ELSE 
        DREL(:)=0.D0 
      END IF
      CALL ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RBOX,DREL,G,NN,POT,E,PSI)
!
      if(tpr.and.thom) then
        AUX(:)=R(:)**2*PSI(:)**2
        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
        CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
        PSI(:)=PSI(:)/SQRT(VAL)
      end if
        
      if(tpr) then
        CALL ATOMLIB_WRITEPHI('PSI_I.DAT',GID,NR,1,PSI)
        PRINT*,'EOFI ',0,E
      end if
!
!     ==========================================================================
!     ==  now start loop                                                      ==
!     ==========================================================================
      do iter=1,niter
!
!       ========================================================================
!       ==  construct |dpsi>=de/dpsi                                          ==
!       ========================================================================
        if(trel) then
          CALL SCHROEDINGER$DREL(GID,NR,POT,E,DREL)
          CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
        else
          drel=0.d0
          rdprime=0.d0
        end if
        if(so.eq.0) then
          sofactor=0.d0
        else if(so.eq.1) then
          sofactor=real(l,kind=8)
        else if(so.eq.-1) then
          sofactor=real(-l-1,kind=8)
        else
          call error$stop('illegal value for so')
          call error$stop('atomlib$boundstatewithHF')
        end if
!       == kinetic energy of the wave function (nonrelativistic)
        call radial$verletd2(gid,nr,r*psi,aux2)
        aux(:)=r(:)*aux2(:)-real(l*(l+1),kind=8)*psi(:)
        aux(:)=-(1.d0+drel)*aux(:)
        call radial$verletd1(gid,nr,psi,aux2)
        aux(:)=aux(:)-rdprime*(r(:)**2*aux2(:)-sofactor*r(:)*psi(:))
!       == potential energy ==================================================
        aux(:)=aux(:)+2.d0*r(:)**2*pot(:)*y0*psi(:)
!       == Hartree-Fock correction ===========================================
        call radialfock$vpsi(GID,NR,vfock,l,psi,vfockpsi)
        aux(:)=aux(:)+2.d0*r(:)**2*vfockpsi(:)
!       == correct glitches ==================================================
        aux(1:2)=0.d0
!       == estimate energy ===================================================
        aux1(:)=0.5d0*aux(:)*psi(:)-r(:)**2*psi(:)*g(:)
        call radial$integrate(gid,nr,aux1,aux2)
        call radial$value(gid,nr,aux2,rbox,e)
        aux1(:)=r(:)**2*psi(:)**2
        call radial$integrate(gid,nr,aux1,aux2)
        call radial$value(gid,nr,aux2,rbox,norm)
        e=e/norm
!       == remove energy term ================================================
        aux(:)=aux(:)-2.d0*r(:)**2*e*psi(:)
!       == weight with occupation ============================================
        dpsi(:)=aux(:)/(2.d0*r(:)**2)-g(:)
!
!       ========================================================================
!       ==  Determine phiprime-phidot*eprime                                  ==
!       ========================================================================
        iso=so
        CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT,E,IRCL,IROUT)
!       == BOUNDARY CONDITION PHI(ROUT)=0 ====================================
        IF(R(IROUT).lt.RBOX) THEN
          rout=r(irout)
        else
          rout=rbox
          IROUT=IRBOX
        END IF
!       ==  SET KINETIC ENERGY TO ZERO BEYOND ROUT TO AVOID AN OVERFLOW ======
!       == ATTENTION: THERE IS A STEP IN THE POTENTIAL =======================
        POT1(:)=POT(:)
        POT1(IROUT:)=E
        IF(TREL) THEN
          CALL SCHROEDINGER$DREL(GID,NR,POT1,E,DREL)
        ELSE 
          DREL(:)=0.D0 
        END IF
        g1=dpsi(:)
        g1(1:2)=0.d0
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,iso,G1,L,E,1,aux1)
        g1=psi(:)
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,iso,G1,L,E,1,aux2)
        call radial$value(gid,nr,aux1,rout,val1)
        call radial$value(gid,nr,aux2,rout,val2)
        dpsi(:)=aux1(:)-aux2(:)*val1/val2
        dpsi(irout+3:)=0.d0
        if(rout.lt.rbox) dpsi(irout:)=0.d0
        call radial$value(gid,nr,dpsi(:),rbox,val)
!   
!       ========================================================================
!       ==  check convergence                                                 ==
!       ========================================================================
        val=maxval(abs(dpsi(:irbox-3)))/sqrt(norm)
        if(tpr) print*,' dev=',val,' iter=',iter
        convg=(val.lt.tol)
        if(convg) exit
!   
!       ========================================================================
!       ==  propagate                                                         ==
!       ========================================================================
        psi(:)=psi(:)-dpsi(:)
!   
!       ========================================================================
!       ==  normalize                                                         ==
!       ========================================================================
        if(tpr.and.thom) then
          aux(:)=r(:)**2*psi(:)**2
          call radial$integrate(gid,nr,aux,aux1)
          call radial$value(gid,nr,aux1,rbox,val)
          psi(:)=psi(:)/sqrt(val)
        end if
        IF(TPR) THEN
          PRINT*,'EOFI ',iter,E
        END IF
!
      ENDDO    ! end of iteration
!
      IF(.NOT.CONVG) THEN
        call error$msg('LOOP FOR RADIAL HARTREE FOCK NOT CONVERGED')
        call error$stop('atomlib$boundstatewithHF')
      end if
      if(tpr) PRINT*,ITER,'ITERATIONSIN LOOP FOR RADIAL HARTREE FOCK'
      if(tpr) then
        PRINT*,'EOFI ',ITER,E
        CALL ATOMLIB_WRITEPHI('PSI_F.DAT',GID,NR,1,PSI)
      end if
!   
!     ==========================================================================
!     ==  determine kinetic energy                                            ==
!     ==========================================================================
      tpsi(:)=(e-pot(:)*y0+vfock%scale*vfock%mux*y0)*psi(:) &
     &                    -vfock%scale*vfockpsi(:)+g(:)
      RETURN
      end subroutine atomlib$boundstatewithHF
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine atomlib$phaseshiftstatewithHF(gid,nr,l,so,drel,g,pot &
     &                             ,vfock,rbnd,phase,e,psi)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      use radialfock_module, only : vfock_type
      implicit none
      integer(4),intent(in)    :: gid
      integer(4),intent(in)    :: nr
      integer(4),intent(in)    :: l
      integer(4),intent(in)    :: so
      real(8)   ,intent(in)    :: drel(nr)
      real(8)   ,intent(in)    :: pot(nr)
      type(vfock_type),intent(in) :: vfock
      real(8)   ,intent(in)    :: g(nr)
      real(8)   ,intent(in)    :: rbnd
      real(8)   ,intent(in)    :: phase
      real(8)   ,intent(out)   :: e
      real(8)   ,intent(out)   :: psi(nr)
      real(8)                  :: dpsi(nr)
      real(8)                  :: phidot(nr)
      real(8)                  :: phiprime(nr)
      real(8)                  :: valdot,derdot,valprime,derprime,val0,der0
      real(8)                  :: pi,y0
      real(8)                  :: sofactor
      real(8)                  :: r(nr)
      real(8)                  :: aux(nr),aux1(nr),aux2(nr)
      real(8)                  :: g1(nr)
      real(8)                  :: rdprime(nr)
      real(8)                  :: val,der,svar
      integer(4)               :: ir
      integer(4),parameter     :: niter=20
      integer(4)               :: iter
      real(8)   ,parameter     :: tol=1.d-8
      logical(4)               :: convg
      REAL(8)   ,PARAMETER     :: XMAX=1.D+15 !MAX. FACTOR IN THE WAVEFUNCTION
      integer(4)               :: irbnd
      real(8)                  :: norm
      logical(4)               :: tpr=.true.
      logical(4)               :: thom
      integer(4)               :: nn
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      call radial$r(gid,nr,r)
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
!      
      thom=maxval(abs(g)).lt.1.d-8
thom=.false.
      do ir=1,nr
        irbnd=ir
        if(r(ir).gt.rbnd) exit
      enddo
      if(so.eq.0) then
        sofactor=0.d0
      else if(so.eq.1) then
        sofactor=real(l,kind=8)
      else if(so.eq.-1) then
        sofactor=real(-l-1,kind=8)
      else
        call error$stop('illegal value for so')
        call error$stop('atomlib$boundstatewithHF')
      end if
!
!     ==========================================================================
!     ==  RECREATE WAVE FUNCTION                                              ==
!     ==========================================================================
      CALL ATOMLIB$phaseshiftstate(GID,NR,L,SO,DREL,G,POT,rbnd,phase,E,PSI)
!
      if(thom) then
        AUX(:)=R(:)**2*PSI(:)**2
        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
        CALL RADIAL$VALUE(GID,NR,AUX1,rbnd,VAL)
        PSI(:)=PSI(:)/SQRT(VAL)
print*,'wave function normalized; thom=',thom
      end if
        
      if(tpr) then
        CALL ATOMLIB_WRITEPHI('PSI_I.DAT',GID,NR,1,PSI)
        PRINT*,'EOFI ',0,E
      end if
!
!     ==========================================================================
!     ==  now start loop                                                      ==
!     ==========================================================================
      do iter=1,niter
!
!       ========================================================================
!       ==  construct |dpsi>=de/dpsi                                          ==
!       ========================================================================
!       == kinetic energy of the wave function (nonrelativistic)
        call radial$verletd2(gid,nr,r*psi,aux1)
        aux(:)=r(:)*aux1(:)
        aux(:)=aux(:)-real(l*(l+1),kind=8)*psi(:)
        aux(:)=-(1.d0+drel)*aux(:)
        call radial$verletd1(gid,nr,psi,aux1)
        aux(:)=aux(:)-rdprime(:)*(r(:)**2*aux1(:)-sofactor*r(:)*psi(:))
!       == potential energy ====================================================
        aux(:)=aux(:)+2.d0*r(:)**2*pot(:)*y0*psi(:)
!       == remove factor 2*r**2 ================================================
        aux(:)=aux(:)/(2.d0*r(:)**2)
!       == Hartree-Fock correction =============================================
        call radialfock$vpsi(GID,NR,vfock,l,psi,aux1)
        aux(:)=aux(:)+aux1(:)
!       == correct glitches ====================================================
        aux(1:2)=0.d0
!       == estimate energy =====================================================
        aux1(:)=r(:)**2*psi(:)*(aux(:)-g(:))
        call radial$integrate(gid,nr,aux1,aux2)
        call radial$value(gid,nr,aux2,rbnd,e)
        aux1(:)=r(:)**2*psi(:)**2
        call radial$integrate(gid,nr,aux1,aux2)
        call radial$value(gid,nr,aux2,rbnd,norm)
        e=e/norm
print*,'ee',e,norm
!       == remove energy term ================================================
        aux(:)=aux(:)-e*psi(:)
!       == subtract inhomogeneity ==============================================
        aux(:)=aux(:)-g(:)
!       == map into dpsi
        dpsi(:)=-aux(:)
!
!       ========================================================================
!       ==  Determine phiprime-phidot*eprime                                  ==
!       ========================================================================
        g1=dpsi(:)
        g1(1:2)=0.d0
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,so,G1,L,E,1,phiprime)
        g1=psi(:)
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,so,G1,L,E,1,phidot)
!
!       :: remark: schroedinger$phaseshift takes value and derivative from a  ::
!       :: two-point formula, and therefore is not fully consistent with      ::
!       :: radial$value and radial$derivative                                 ::
        call radial$value(gid,nr,psi,rbnd,val0)
        call radial$derivative(gid,nr,PSI,rbnd,der0)
        call radial$value(gid,nr,phiprime,rbnd,valprime)
        call radial$derivative(gid,nr,phiprime,rbnd,derprime)
        call radial$value(gid,nr,phidot,rbnd,valdot)
        call radial$derivative(gid,nr,phidot,rbnd,derdot)
        call SCHROEDINGER$resolvePHASESHIFT(PHASE,val,der,nn)
        val0=val*der0-der*val0
        valprime=val*derprime-der*valprime
        valdot=val*derdot-der*valdot
        svar=(val0+valprime)/valdot  !delta epsilon
        DPSI(:)=phiprime(:)-phidot(:)*svar
!   
!       ========================================================================
!       ==  check convergence                                                 ==
!       ========================================================================
        val=maxval(abs(dpsi(:irbnd-3)))/sqrt(norm)
        if(tpr) print*,' dev=',val
        convg=(val.lt.tol)
        if(convg) exit
!   
!       ========================================================================
!       ==  propagate                                                         ==
!       ========================================================================
        psi(:)=psi(:)+dpsi(:)
! CALL SCHROEDINGER$PHASESHIFT(GID,NR,Psi,Rbnd,svar)
!print*,'target phase+ ',phase,' actual phase=',svar,' r=',rbnd
!   
!       ========================================================================
!       ==  normalize                                                         ==
!       ========================================================================
        if(thom) then
          aux(:)=r(:)**2*psi(:)**2
          call radial$integrate(gid,nr,aux,aux1)
          call radial$value(gid,nr,aux1,rbnd,val)
          psi(:)=psi(:)/sqrt(val)
        end if
        IF(TPR) THEN
          PRINT*,'EOFI ',iter,E
        END IF
!
      ENDDO    ! end of iteration
!
      IF(.NOT.CONVG) THEN
        call error$msg('LOOP FOR RADIAL HARTREE FOCK NOT CONVERGED')
        call error$stop('atomlib$phaseshiftstatewithHF')
      end if
      if(tpr) PRINT*,ITER,'ITERATIONSIN LOOP FOR RADIAL HARTREE FOCK'
      if(thom) then
        aux(:)=r(:)**2*psi(:)**2
        call radial$integrate(gid,nr,aux,aux1)
        call radial$value(gid,nr,aux1,rbnd,val)
        psi(:)=psi(:)/sqrt(val)
      end if
      if(tpr) then
        PRINT*,'EOFI ',ITER,E
        CALL ATOMLIB_WRITEPHI('PSI_F.DAT',GID,NR,1,PSI)
      end if
!
      RETURN
      end subroutine atomlib$phaseshiftstatewithHF
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine atomlib$updatestatewithHF(gid,nr,l,so,drel,g,pot,vfock &
     &                                    ,rbnd,e,psi)
!     **************************************************************************
!     ** uses an approximate solution to construct a solution for             **
!     ** a hamiltonian with a Hartee-fock contribution mixed in               **
!     **                                                                      **
!     ** The boundary conditions (logarithmic derivative at rbnd is preserved **
!     ** If rbnd.le.0 a fixed energy solution is searched and the outer       **
!     ** boundary condition is relaxed.                                       **
!     **                                                                      **
!     **************************************************************************
      use radialfock_module, only : vfock_type
      implicit none
      integer(4)      ,intent(in) :: gid
      integer(4)      ,intent(in) :: nr
      integer(4)      ,intent(in) :: l
      integer(4)      ,intent(in) :: so
      real(8)         ,intent(in) :: drel(nr)
      real(8)         ,intent(in) :: pot(nr)
      type(vfock_type),intent(in) :: vfock
      real(8)   ,intent(in)    :: g(nr)
      real(8)   ,intent(in)    :: rbnd
      real(8)   ,intent(inout) :: e
      real(8)   ,intent(inout) :: psi(nr)
      real(8)                  :: phase(nr)
      real(8)                  :: dpsi(nr)
      real(8)                  :: phidot(nr)
      real(8)                  :: phiprime(nr)
      real(8)                  :: valdot,derdot,valprime,derprime,val0,der0
      real(8)                  :: pi,y0
      real(8)                  :: sofactor
      real(8)                  :: r(nr)
      real(8)                  :: aux(nr),aux1(nr),aux2(nr)
      real(8)                  :: g1(nr)
      real(8)                  :: rdprime(nr)
      real(8)                  :: val,der,svar
      integer(4)               :: ir
      integer(4),parameter     :: niter=50
      integer(4)               :: iter
      real(8)   ,parameter     :: tol=1.d-8
      logical(4)               :: convg
      REAL(8)   ,PARAMETER     :: XMAX=1.D+15 !MAX. FACTOR IN THE WAVEFUNCTION
      integer(4)               :: irbnd
      real(8)                  :: norm
      logical(4)               :: tpr=.false.
      logical(4)               :: tfixe
      integer(4)               :: nn
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      call radial$r(gid,nr,r)
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
      tfixe=rbnd.le.0.d0
!
!     :: remark: schroedinger$phaseshift takes value and derivative from a    ::
!     :: two-point formula, and therefore is not fully consistent with        ::
!     :: radial$value and radial$derivative                                   ::
      CALL SCHROEDINGER$PHASESHIFT(GID,NR,Psi,Rbnd,phase)
!
      do ir=1,nr
        irbnd=ir
        if(tfixe.and.r(ir).gt.3.d0) exit
        if(.not.tfixe.and.r(ir).gt.rbnd) exit
      enddo

      if(so.eq.0) then
        sofactor=0.d0
      else if(so.eq.1) then
        sofactor=real(l,kind=8)
      else if(so.eq.-1) then
        sofactor=real(-l-1,kind=8)
      else
        call error$stop('illegal value for so')
        call error$stop('atomlib$boundstatewithHF')
      end if
!
!     ==========================================================================
!     ==  now start loop                                                      ==
!     ==========================================================================
      do iter=1,niter
!
!       ========================================================================
!       ==  construct |dpsi>=de/dpsi                                          ==
!       ========================================================================
!       == kinetic energy of the wave function (nonrelativistic)
        call radial$verletd2(gid,nr,r*psi,aux1)
        aux(:)=r(:)*aux1(:)
        aux(:)=aux(:)-real(l*(l+1),kind=8)*psi(:)
        aux(:)=-(1.d0+drel)*aux(:)
        call radial$verletd1(gid,nr,psi,aux1)
        aux(:)=aux(:)-rdprime(:)*(r(:)**2*aux1(:)-sofactor*r(:)*psi(:))
!       == potential energy ====================================================
        aux(:)=aux(:)+2.d0*r(:)**2*pot(:)*y0*psi(:)
!       == remove factor 2*r**2 ================================================
        aux(:)=aux(:)/(2.d0*r(:)**2)
!       == Hartree-Fock correction =============================================
        call radialfock$vpsi(GID,NR,vfock,l,psi,aux1)
        aux(:)=aux(:)+aux1(:)
!       == correct glitches ====================================================
        aux(1:2)=0.d0
!       == estimate energy =====================================================
        if(.not.tfixe) then
          aux1(:)=r(:)**2*psi(:)*(aux(:)-g(:))
          call radial$integrate(gid,nr,aux1,aux2)
          call radial$value(gid,nr,aux2,rbnd,e)
          aux1(:)=r(:)**2*psi(:)**2
          call radial$integrate(gid,nr,aux1,aux2)
          call radial$value(gid,nr,aux2,rbnd,norm)
          e=e/norm
!         == remove energy term ================================================
          aux(:)=aux(:)-e*psi(:)
        else
!         = a reasonable norm is needed later for judging the convergence ======
          aux1(:)=r(:)**2*psi(:)**2
          call radial$integrate(gid,nr,aux1,aux2)
          call radial$value(gid,nr,aux2,3.d0,norm)
        end if
!       == subtract inhomogeneity ==============================================
        aux(:)=aux(:)-g(:)
!       == map into dpsi
        dpsi(:)=-aux(:)
!
!       ========================================================================
!       ==  Determine phiprime-phidot*eprime                                  ==
!       ========================================================================
        g1=dpsi(:)
        g1(1:2)=0.d0
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,so,G1,L,E,1,phiprime)
        if(.not.tfixe) then
          g1=psi(:)
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,so,G1,L,E,1,phidot)
          call radial$value(gid,nr,psi,rbnd,val0)
          call radial$derivative(gid,nr,PSI,rbnd,der0)
          call radial$value(gid,nr,phiprime,rbnd,valprime)
          call radial$derivative(gid,nr,phiprime,rbnd,derprime)
          call radial$value(gid,nr,phidot,rbnd,valdot)
          call radial$derivative(gid,nr,phidot,rbnd,derdot)
          call SCHROEDINGER$resolvePHASESHIFT(PHASE,val,der,nn)
          val0=val*der0-der*val0
          valprime=val*derprime-der*valprime
          valdot=val*derdot-der*valdot
          svar=(val0+valprime)/valdot  !delta epsilon
          DPSI(:)=phiprime(:)-phidot(:)*svar
        else
          CALL ATOMLIB_WRITEPHI('g.dat',GID,NR,1,g1)
          CALL ATOMLIB_WRITEPHI('phiprime.dat',GID,NR,1,phiprime)
          DPSI(:)=phiprime(:)
        end if
!   
!       ========================================================================
!       ==  check convergence                                                 ==
!       ========================================================================
        val=maxval(abs(dpsi(:irbnd-3)))/sqrt(norm)
        if(tpr) print*,' dev=',val
        convg=(val.lt.tol)
        if(convg) exit
!
!       ========================================================================
!       ==  propagate                                                         ==
!       ========================================================================
        psi(:)=psi(:)+dpsi(:)
!
      ENDDO    ! end of iteration
!
      IF(.NOT.CONVG) THEN
        call error$msg('LOOP FOR RADIAL HARTREE FOCK NOT CONVERGED')
        call error$stop('atomlib$updatestatewithHF')
      end if
      RETURN
      end subroutine atomlib$updatestatewithHF
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine atomlib$sphericalwithHF(gid,nr,pot,drel,vfock,so,g,l,e,psi)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      use radialfock_module, only : vfock_type
      implicit none
      integer(4),intent(in)    :: gid
      integer(4),intent(in)    :: nr
      type(vfock_type),intent(in) :: vfock
      real(8)   ,intent(in)    :: drel(nr)
      real(8)   ,intent(in)    :: pot(nr)
      integer(4),intent(in)    :: l
      integer(4),intent(in)    :: so
      real(8)   ,intent(in)    :: g(nr)
      real(8)   ,intent(in)    :: e
      real(8)   ,intent(out)   :: psi(nr)
      real(8)                  :: dpsi(nr)
      real(8)                  :: vfockpsi(nr)
      real(8)                  :: pi,y0
      real(8)                  :: sofactor
      real(8)                  :: r(nr)
      real(8)                  :: aux(nr),aux1(nr),aux2(nr)
      real(8)                  :: g1(nr)
      real(8)                  :: rdprime(nr)
      real(8)                  :: val,norm
      integer(4),parameter     :: niter=20
      integer(4)               :: iter
      real(8)   ,parameter     :: tol=1.d-5
      logical(4)               :: convg
      REAL(8)   ,PARAMETER     :: XMAX=1.D+15 !MAX. FACTOR IN THE WAVEFUNCTION
      logical(4)               :: tpr=.true.
      logical(4)               :: thom
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      call radial$r(gid,nr,r)
      thom=maxval(abs(g)).lt.1.d-8
!
      IF(TPR) THEN
        print*,'==============================================================='
        PRINT*,'l=',l,' so=',so,' thom ',thom
        PRINT*,'SCALE ',vfock%SCALE
        PRINT*,'EOFI ',-1,E
      END IF
!
!     ==========================================================================
!     ==  check and process spin orbit parameter                              ==
!     ==========================================================================
      if(so.eq.0) then
        sofactor=0.d0
      else if(so.eq.1) then
        sofactor=real(l,kind=8)
      else if(so.eq.-1) then
        sofactor=real(-l-1,kind=8)
      else
        call error$stop('illegal value for so')
        call error$stop('atomlib$boundstatewithHF')
      end if
!
!     ==========================================================================
!     ==  set relativistic correction                                         ==
!     ==========================================================================
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
!
!     ==========================================================================
!     ==  RECREATE WAVE FUNCTION                                              ==
!     ==========================================================================
      CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,so,G,L,E,1,psi)
      if(tpr) then
        CALL ATOMLIB_WRITEPHI('PSI_I.DAT',GID,NR,1,PSI)
      end if
!
!     ==========================================================================
!     ==  now start loop                                                      ==
!     ==========================================================================
      do iter=1,niter
!
!       ========================================================================
!       ==  construct |dpsi>=de/dpsi                                          ==
!       ========================================================================
!       == kinetic energy of the wave function (nonrelativistic)
        call radial$verletd2(gid,nr,r*psi,aux2)
        aux(:)=r(:)*aux2(:)-real(l*(l+1),kind=8)*psi(:)
        aux(:)=-(1.d0+drel)*aux(:)
        call radial$verletd1(gid,nr,psi,aux2)
        aux(:)=aux(:)-rdprime*(r(:)**2*aux2(:)-sofactor*r(:)*psi(:))
!       == potential energy ==================================================
        aux(:)=aux(:)+2.d0*r(:)**2*pot(:)*y0*psi(:)
!       == Hartree-Fock correction ===========================================
        call radialfock$vpsi(GID,NR,vfock,l,psi,vfockpsi)
        aux(:)=aux(:)+2.d0*r(:)**2*vfockpsi(:)
!       == correct glitches ==================================================
        aux(1:2)=0.d0
!       == remove energy term ================================================
        aux(:)=aux(:)-2.d0*r(:)**2*e*psi(:)
!       == weight with occupation ============================================
        dpsi(:)=aux(:)/(2.d0*r(:)**2)-g(:)
!
!       ========================================================================
!       ==  Determine phiprime-phidot*eprime                                  ==
!       ========================================================================
        CALL ATOMLIB_WRITEPHI('dpsi1.dat',GID,NR,1,dPSI)
        g1=dpsi(:)
        g1(1:2)=0.d0
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,so,G1,L,E,1,dpsi)
        CALL ATOMLIB_WRITEPHI('dpsi2.dat',GID,NR,1,dPSI)
!   
!       ========================================================================
!       ==  check convergence                                                 ==
!       ========================================================================
        aux(:)=(r*psi)**2
        call radial$integrate(gid,nr,aux,aux1)
        call radial$value(gid,nr,aux1,3.d0,norm)
        val=maxval(abs(dpsi(:)))/sqrt(norm)
        if(tpr) print*,' dev=',val,norm
        convg=(val.lt.tol)
        if(convg) exit
!   
!       ========================================================================
!       ==  propagate                                                         ==
!       ========================================================================
        psi(:)=psi(:)-dpsi(:)
!
      ENDDO    ! end of iteration
!
      IF(.NOT.CONVG) THEN
        CALL ATOMLIB_WRITEPHI('PSI_F.DAT',GID,NR,1,PSI)
        call error$msg('LOOP FOR RADIAL HARTREE FOCK NOT CONVERGED')
        call error$stop('atomlib$sphericalwithHF')
      end if
      if(tpr) then
        PRINT*,ITER,'ITERATIONSIN LOOP FOR RADIAL HARTREE FOCK'
        CALL ATOMLIB_WRITEPHI('PSI_F.DAT',GID,NR,1,PSI)
      end if
      RETURN
      end subroutine atomlib$sphericalwithHF
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$etotwithfoCK(GID,NR,rbox,lrhox,trel,pot,vfock &
    &                        ,NB,LOFI,sofi,FOFI,PSI,etot)
!     **************************************************************************
!     **************************************************************************
      use radialfock_module, only : vfock_type
      implicit none
      integer(4),intent(in) :: gid
      integer(4),intent(in) :: nr
      real(8)   ,intent(in) :: rbox
      integer(4),intent(in) :: nb
      integer(4),intent(in) :: lofi(nb)
      integer(4),intent(in) :: sofi(nb)
      logical(4),intent(in) :: trel
      integer(4),intent(in) :: lrhox
      real(8)   ,intent(in) :: fofi(nb)
      real(8)   ,intent(in) :: psi(nr,nb)
      real(8)   ,intent(in) :: pot(nr)
      type(vfock_type),intent(in) :: vfock
      real(8)   ,intent(out):: etot
      real(8)               :: drrpsi(nr,nb)
      real(8)               :: aux(nr),aux1(nr),aux2(nr)
      real(8)               :: r(nr)
      real(8)               :: eofi(nb)
      real(8)               :: drel(nr),rdprime(nr)
      real(8)               :: val,der
      real(8)               :: pi,y0
      integer(4)            :: irbox
      real(8)               :: sofactor
      integer(4)            :: ir,ib
      integer(4)            :: l
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      call radial$r(gid,nr,r)
!
!     == index of the first gridpoint outside rbox =============================
      do ir=1,nr
        irbox=ir
        if(r(ir).ge.rbox) exit
      enddo
!
!     ==========================================================================
!     ==  vfock|psi> and partial_r*r|\psi>                                    ==
!     ==========================================================================
      do ib=1,nb
        aux(:)=r(:)*psi(:,ib)
        val=aux(irbox)
        der=(aux(irbox+1)-aux(irbox-1))/(r(irbox+1)-r(irbox-1))
        aux(irbox:)=val+(r(irbox:)-r(irbox))*der
        call radial$verletd1(gid,nr,aux,aux1)
        aux1(1:2)=aux1(3)  
        drrpsi(:,ib)=aux1(:)
      enddo
!
!     ==========================================================================
!     ==  estimate energies for relativistic corrections                      ==
!     ==========================================================================
      do ib=1,nb
        eofi(ib)=0.d0
      enddo
!
!     ========================================================================
!     ==  calculate total energy (this a fake quantity)                     ==
!     ========================================================================
      aux(:)=0.d0
      do ib=1,nb
        l=lofi(ib)
        if(trel) then
          CALL SCHROEDINGER$DREL(GID,NR,POT,EOFI(IB),DREL)
          CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
        else
          drel=0.d0
          rdprime=0.d0
        end if
        if(sofi(ib).eq.0) then
          sofactor=0.d0
        else if(sofi(ib).eq.1) then
          sofactor=real(l,kind=8)
        else if(sofi(ib).eq.-1) then
          sofactor=real(-l-1,kind=8)
        else
          call error$stop('illegal value for sofi ')
          call error$stop('atomlib$aescfwithHF')
        end if
!       == kinetic energy ====================================================
        aux1(:)=0.5d0*(1.d0+drel) &
     &               *(drrpsi(:,ib)**2+real(l*(l+1),kind=8)*psi(:,ib)**2)
        aux1(:)=aux1(:)+0.5d0*(sofactor+1.d0)*rdprime*r(:)*psi(:,ib)**2
!       == potential energy ==================================================
        aux1(:)=aux1(:)+r(:)**2*pot(:)*y0*psi(:,ib)**2
!       == Hartree-Fock correction ===========================================
        call radialfock$vpsi(GID,NR,vfock,lofi(ib),psi(:,ib),aux2)
        aux2(:)=psi(:,ib)*(0.5d0*aux2-0.5d0*vfock%mux(:)*y0*psi(:,ib))
        aux1(:)=aux1(:)+r(:)**2*aux2(:)
!       == add to sum over states ============================================
        aux(:)=aux(:)+fofi(ib)*aux1(:)
      enddo
      call radial$integrate(gid,nr,aux,aux1)
      call radial$value(gid,nr,aux1,rbox,etot)
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE atomlib_WRITEPHI(FILE,GID,NR,NPHI,PHI)
!     **                                                                      **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: NPHI        
      REAL(8)     ,INTENT(IN) :: PHI(NR,NPHI)
      INTEGER(4)              :: IR
      REAL(8)                 :: R(NR)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      DO IR=1,NR
        IF(R(IR).GT.3.D0.AND.MAXVAL(ABS(PHI(IR,:))).GT.1.D+3) EXIT
        WRITE(100,FMT='(F15.10,2X,20(F25.10,2X))')R(IR),PHI(IR,:)
      ENDDO
      CLOSE(100)
      RETURN
      END
