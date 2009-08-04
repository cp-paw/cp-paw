!........1.........2.........3.........4.........5.........6.........7.........8
module localize_module
!*******************************************************************************
!*******************************************************************************
!****  dlist contains the expansion coefficients of the correlated orbitals ****
!****        centered at iat1 in partial waves centered at iat2 translated  ****
!****        by the translations it
!*******************************************************************************
!*******************************************************************************
!
!      call LOCALIZE$NBLIST()
!      do ik=1,nk
!        call localize$chidenmat
!      enddo
!
!
type nblist_type
  integer(4)         :: iat1
  integer(4)         :: iat2
  integer(4)         :: it(3)
end type nblist_type
type matlist_type
  integer(4)         :: iat1
  integer(4)         :: iat2
  integer(4)         :: it(3)
  complex(8),pointer :: mat(:,:) 
  complex(8),pointer :: mat1(:,:,:,:) 
end type matlist_type
logical(4)        ,parameter   :: tonsite=.true.
real(8)           ,parameter   :: rmax=6.d0
integer(4)        ,save        :: nnbx=0   ! x#(neighbors)
integer(4)        ,save        :: nnb      ! #(neighbors)
type(nblist_type) ,allocatable :: nblist(:)! neighbor list
integer(4)        ,save        :: npaird=0   ! length of dlist
type(matlist_type),allocatable :: dlist(:) ! expansion coefficients of 
                                           ! local orbitals into partial waves
type(matlist_type),allocatable :: chidenmat(:) ! expansion coefficients of 
                                           ! local orbitals into partial waves
logical(4),save        :: tini=.false.
integer(4),save        :: nat              ! #(atoms)
integer(4),save        :: ndim
integer(4),save        :: nphi
integer(4),save        :: nchi
integer(4),allocatable :: lnxchi(:)
integer(4),allocatable :: loxchi(:,:)
integer(4),allocatable :: lnxphi(:)
integer(4),allocatable :: loxphi(:,:)
integer(4),allocatable :: iposchi(:)
integer(4),allocatable :: iposphi(:)
contains
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LOCALIZE_CLEANCHIDENMAT()
      IMPLICIT NONE
      INTEGER(4) :: LENG,I
!     **************************************************************************
      IF(.NOT.Allocated(CHIDENMAT)) RETURN
      LENG=SIZE(CHIDENMAT)
      DO I=1,LENG
        DEALLOCATE(CHIDENMAT(I)%MAT1)
      ENDDO
      DEALLOCATE(CHIDENMAT)
      RETURN
      END SUBROUTINE LOCALIZE_CLEANCHIDENMAT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LOCALIZE_INICHIDENMAT()
      IMPLICIT NONE
      INTEGER(4) :: I,IAT1,IAT2,LMNX1,LMNX2
!     **************************************************************************
      IF(Allocated(CHIDENMAT)) CALL LOCALIZE_CLEANCHIDENMAT()
      ALLOCATE(CHIDENMAT(NNB))
      DO I=1,NNB
        IAT1=NBLIST(I)%IAT1
        IAT2=NBLIST(I)%IAT2
        CHIDENMAT(I)%IAT1=IAT1
        CHIDENMAT(I)%IAT2=IAT2
        CHIDENMAT(I)%IT(:)=NBLIST(I)%IT(:)
        LMNX1=SUM(2*LOXCHI(:LNXCHI(IAT1),IAT1)+1)             
        LMNX2=SUM(2*LOXCHI(:LNXCHI(IAT2),IAT2)+1)             
        ALLOCATE(CHIDENMAT(I)%MAT1(LMNX1,LMNX2,NDIM,NDIM))
        CHIDENMAT(I)%MAT1(:,:,:,:)=(0.D0,0.D0)
      ENDDO
      RETURN
      END SUBROUTINE LOCALIZE_INICHIDENMAT
END MODULE LOCALIZE_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LOCALIZE_initialize()
!     **************************************************************************
!     ** SET UP A NEIGHBORLIST NBLIST WITH RADIUS RMAX                        **
!     **************************************************************************
      use localize_module
      implicit none
      integer(4)             :: iat,isp,lnxx,ln,l
      integer(4)             :: NSP           !#(ATOM TYPES (SPECIES))
      integer(4),ALLOCATABLE :: ISPECIES(:)
      integer(4),ALLOCATABLE :: LOXCHI1(:)
      integer(4),ALLOCATABLE :: LOXPHI1(:)
      INTEGER(4)             :: LNXPHI1
      INTEGER(4)             :: LNXCHI1
      INTEGER(4)             :: LmNXpHI
      INTEGER(4)             :: LmNXCHI
      real(8)                :: rcut
      integer(4)             :: nofl(4)
!     **************************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
      CALL WAVES$GETI4('SPINORDIM',NDIM)
!
!     ==========================================================================
!     == communicate data from ldaplusu object to setup object                ==
!     ==========================================================================
      CALL SETUP$NSPECIES(NSP)
      do isp=1,nsp
        CALL LDAPLUSU$SELECT(ISP)
        CALL LDAPLUSU$GETR8('RCUT',RCUT)
        CALL LDAPLUSU$GETI4A('NCORROFL',4,NOFL)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$SETR8('RADCHI',RCUT)
        CALL SETUP$SETI4A('NOFLCHI',4,NOFL)
      enddo
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      CALL SETUP$LNXX(LNXX)
      ALLOCATE(LNXPHI(NAT))
      ALLOCATE(LNXCHI(NAT))
      ALLOCATE(LOXPHI(LNXX,NAT))
      ALLOCATE(LOXCHI(LNXX,NAT))
      ALLOCATE(LOXPHI1(LNXX))
      ALLOCATE(LOXCHI1(LNXX))
      nchi=0
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNXPHI1) 
        CALL SETUP$GETI4('LNXCHI',LNXCHI1) 
        LOXPHI1(:)=0
        LOXCHI1(:)=0
        CALL SETUP$GETI4A('LOX',LNXPHI1,LOXPHI1(1:LNXPHI1))
        CALL SETUP$GETI4A('LOXCHI',LNXCHI1,LOXCHI1(1:LNXCHI1))
        lmnxphi=sum(2*loxphi1(1:lnxphi1)+1)
        lmnxchi=sum(2*loxchi1(1:lnxchi1)+1)
        DO IAT=1,NAT
          IF(ISPECIES(IAT).NE.ISP) CYCLE
          LNXPHI(IAT)=LNXPHI1
          LNXCHI(IAT)=LNXCHI1
          LOXPHI(:,IAT)=LOXPHI1(:)
          LOXCHI(:,IAT)=LOXCHI1(:)
          nchi=nchi+lmnxchi
          nphi=nphi+lmnxphi
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == determine helper array for the placement of submatrices              ==
!     ==========================================================================
      allocate(iposchi(nat+1))
      allocate(iposphi(nat+1))
      iposchi(1)=1
      iposphi(1)=1
      do iat=1,nat
        iposchi(iat+1)=iposchi(iat)
        do ln=1,lnxchi(iat)
          l=loxchi(ln,iat)
          iposchi(iat+1)=iposchi(iat+1)+(2*l+1)
        enddo
        iposphi(iat+1)=iposphi(iat)
        do ln=1,lnxphi(iat)
          l=loxphi(ln,iat)
          iposphi(iat+1)=iposphi(iat+1)+(2*l+1)
        enddo
      enddo         

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine localize_dlistsimple()
!     **************************************************************************
!     **  this is a simple version for the transformation from partial waves  **
!     **  to local orbitals.                                                  **
!     **  here we use the first partial wave for each l as local orbital      **
!     **                                                                      **
!     **                                                                      **
!     ** this routine must be replaced!!!!!!                                  **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      use localize_module
      implicit none
      integer(4) :: iat,l,ln,lx,lnx
!     **************************************************************************
      npaird=nat
      allocate(dlist(npaird))
      do iat=1,nat
        dlist(iat)%iat1=iat       
        dlist(iat)%iat2=iat       
        dlist(iat)%it(:)=0
!
!       == define lnxchi and loxchi =============================================
        lx=maxval(loxphi(:,iat))
        lnxchi(iat)=lx+1
        do l=0,lx
          loxchi(l+1,iat)=l
        enddo
!
        lnx=lnxphi(iat)
        allocate(dlist(iat)%mat(lx+1,lnx))
        dlist(iat)%mat(:,:)=(0.d0,0.d0)
        do l=0,lx
          do ln=1,lnx
            if(loxphi(ln,iat).eq.l) then
              dlist(iat)%mat(l+1,ln)=(1.d0,0.d0)
              exit
            end if
          enddo
        enddo
      enddo
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LOCALIZE$NBLIST()
!     **************************************************************************
!     ** SET UP A NEIGHBORLIST NBLIST WITH RADIUS RMAX                        **
!     **                                                                      **
!     ** USES: RBAS (FROM CELL OBJECT)                                        **
!     **       NAT  (FROM ATOM OBJECT, DATA HELD LOCALLY IN LOCALIZE MODULE   **
!     **       R(3,NAT)  (FROM ATOM OBJECT)                                   **
!     **       Rmax  (parameter from localize module)                         **
!     ** PRODUCES:                                                            **
!     **       NNBX                                                           **
!     **       NNB                                                            **
!     **       NBLIST                                                         **
!     **                                                                      **
!     **************************************************************************
      use localize_module, only : nat,rmax,nnbx,nnb,nblist
      implicit none
      real(8)                :: rbas(3,3)    ! Lattice translations
      real(8)                :: rbasinv(3,3) ! inverse of rbas
      real(8)   ,allocatable :: r(:,:)       ! atomic positions
      real(8)                :: xr(3)
      real(8)                :: x0,y0,z0
      integer(4),allocatable :: itr(:,:)
      integer(4)             :: itvec(3)
      real(8)                :: d(3)         ! distance vector
      real(8)                :: d2           ! squared distance
      real(8)                :: rmax2        ! squared maximum distance
      integer(4)             :: iat,iat1,iat2,it1,it2,it3
      integer(4)             :: min1,max1,min2,max2,min3,max3
      integer(4)             :: i
      integer(4)             :: lmnx1,lmnx2
!     **************************************************************************
      CALL LOCALIZE_INITIALIZE()
!
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL LIB$INVERTR8(3,RBAS,RBASINV)
!     == the array r must be allocatable because the dimension nat may change
!     == within this routine (see above, in localize_initialize)
      ALLOCATE(R(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
!
!     ==========================================================================
!     == obtain translation relative to first unit cell                       ==
!     == r-rbas*itr lies in the first unit cell                               ==
!     ==========================================================================
      allocate(itr(3,nat))
      do iat=1,nat
        xr(:)=matmul(rbasinv(:,:),r(:,iat))
        DO I=1,3
          IF(Xr(I).GE.0.D0) THEN  ! ROUND DOWN
            ITR(I,IAT)=INT(Xr(I))
          ELSE
            ITR(I,IAT)=INT(Xr(I))-1
          END IF
        ENDDO
      enddo
!
!     ==========================================================================
!     == SET UP NEIGHBORLIST                                                  ==
!     ==========================================================================
1000  continue
      rmax2=rmax**2
      NNB=0
      DO IAT1=1,NAT
        X0=R(1,IAT1)
        Y0=R(2,IAT1)
        Z0=R(3,IAT1)
        CALL BOXSPH(Rbas,X0,Y0,Z0,Rmax,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
!       == LOOP OVER BOXES IN THE NEIGHBORHOOD ================================
        DO IT1=MIN1-1,MAX1
          DO IT2=MIN2-1,MAX2
            DO IT3=MIN3-1,MAX3
              DO IAT2=1,NAT
                ITVEC(1)=IT1-ITR(1,IAT2)
                ITVEC(2)=IT2-ITR(2,IAT2)
                ITVEC(3)=IT3-ITR(3,IAT2)
!               == AVOID DOUBLE COUNTING=======================================
                IF(IAT1.GT.IAT2) CYCLE
                IF(IAT2.EQ.IAT1) THEN
                  IF(ITVEC(3).GT.0) CYCLE
                  IF(ITVEC(3).EQ.0) THEN
                    IF(ITVEC(2).GT.0) CYCLE
                    IF(ITVEC(2).EQ.0) THEN
                      IF(ITVEC(1).GE.0) CYCLE
                    END IF
                  END IF
                END IF   
!               == DISTANCE CRITERION =========================================
                D(:)=R(:,IAT2)+MATMUL(RBAS,REAL(ITVEC,KIND=8))-R(:,IAT1)
                D2=SUM(D(:)**2)
                IF(D2.GT.RMAX2) CYCLE
                NNB=NNB+1
                IF(NNB.LE.NNBX) THEN
                  NBLIST(NNB)%IAT1=IAT1
                  NBLIST(NNB)%IAT2=IAT2
                  NBLIST(NNB)%IT(:)=ITVEC(:)
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     == RESIZE NBLIST AND REDO IF NECESSARY ===================================
      IF(NNB.gt.NNBX) THEN
        IF(NNBX.GT.0)DEALLOCATE(NBLIST)
        ALLOCATE(NBLIST(NNB))
        NNBX=NNB
        GOTO 1000
      END IF
!
!     =============================================================================
!     == report for test                                                         ==
!     =============================================================================
!!$      do i=1,nnb
!!$        iat1=nblist(i)%iat1
!!$        iat2=nblist(i)%iat2
!!$        itvec=nblist(i)%it
!!$        D(:)=R(:,IAT2)+MATMUL(RBAS,REAL(ITVEC,KIND=8))-R(:,IAT1)
!!$        D2=sqrt(SUM(D(:)**2))
!!$        write(*,fmt='(i5,2i5,f10.3,3f10.3)')i,iat1,iat2,d2,d/d2
!!$!        write(*,fmt='(i5,2i3,3i2)')i,nblist(i)%iat1,nblist(i)%iat2,nblist(i)%it
!!$      enddo
      RETURN
      END

!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      subroutine localize_chidenmat(xk,tinv,ndim,nbh,nb,OCC,npro_,proj,nchi_,cchi)
!!$!     **************************************************************************
!!$!     ** convert the projections <ptilde|psitilde(k)> into contribions to the **
!!$!     ** density matrix chidenmat in terms of local orbitals.                 **
!!$!     ** The transformation from partial waves to local orbitals is given by  **
!!$!     ** dlist.                                                               **
!!$!     **                                                                      **
!!$!     ** remarks:                                                             **
!!$!     **   This routine adds the contribution from this k-point to chidenmat. **
!!$!     **   Hence chidenmat must be first set to zero and this routine must be **
!!$!     **   called for all k-points                                            **
!!$!     **                                                                      **
!!$!     **************************************************************************
!!$      implicit none
!!$      real(8)   ,intent(in) :: xk(3)
!!$      logical(4),intent(in) :: tinv
!!$      integer(4),intent(in) :: ndim
!!$      integer(4),intent(in) :: nbh
!!$      integer(4),intent(in) :: nb
!!$      REAL(8)   ,INTENT(IN) :: OCC(NB)
!!$      integer(4),intent(in) :: npro_
!!$      integer(4),intent(in) :: nchi_
!!$      complex(8),intent(in) :: proj(ndim,nbh,npro_)      
!!$      complex(8)            :: cchi(nchi_,ndim,nbh)
!!$      real(8)               :: svar2
!!$      integer(4)            :: ipair,idim,j
!!$!     **************************************************************************
!!$!
!!$!     ==========================================================================
!!$!     ==  CHECK IF SUPERWAVEFUNCTIONS ARE USED AND IF #(BANDS) CORRECT        ==
!!$!     ==========================================================================
!!$      IF(TINV) THEN
!!$        IF(NBH.NE.(NB+1)/2) THEN
!!$          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
!!$          CALL ERROR$STOP('WAVES_DENMAT')
!!$        END IF
!!$      ELSE 
!!$        IF(NBH.NE.NB) THEN
!!$          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
!!$          CALL ERROR$STOP('WAVES_DENMAT')
!!$        END IF
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     ==                                                                      ==
!!$!     ==========================================================================
!!$      do ipair=1,npaird
!!$        iat1=chidenmat(ipair)%iat1
!!$        iat2=chidenmat(ipair)%iat2
!!$        it(:)=chidenmat(ipair)%it(:)
!!$        i1=ipos(iat1)
!!$        i2=ipos(iat1+1)-1
!!$        j1=jpos(iat2)
!!$        j2=jpos(iat2+1)-1
!!$        kt=xk(1)*real(it(1),kind=8) &
!!$     &    +xk(2)*real(it(2),kind=8) &
!!$     &    +xk(3)*real(it(3),kind=8) 
!!$        kt=twopi*kt
!!$        eikt=(cos(kt),-sin(kt))
!!$        ib=0
!!$        do ibh=1,nbh
!!$          IF(TINV) THEN
!!$            SVAR1=0.5D0*(OCC(2*IB-1)+OCC(2*IB))
!!$            SVAR2=0.5D0*(OCC(2*IB-1)-OCC(2*IB))
!!$          ELSE
!!$            SVAR1=OCC(IB)
!!$          END IF
!!$!
!!$          do j=1,j2-j1+1
!!$            DO IDIM1=1,NDIM
!!$              IF(TINV) THEN
!!$!               == FUNC=(F1+F2)/2<PSI1-I*PSI2|P>+(F1-F2)/2<PSI1+I*PSI2|P>
!!$                FUNC(j,IDIM1)=SVAR1*CONJG(cchi(IDIM1,j1-1+j,ibh)) &
!!$     &                       +SVAR2*cchi(IDIM1,j1-1+j,ibh)
!!$              ELSE 
!!$                FUNC(j,IDIM1)=SVAR1*CONJG(cchiI(IDIM1,j1-1+j,IB))
!!$              END IF
!!$            ENDDO
!!$          ENDDO
!!$          FUNC(:)=FUNC(:)*EIKT
!!$          DO j=1,j2-j1+1
!!$            DO IDIM=1,NDIM
!!$              CHIDENMAT(IPAIR)%MAT(:,j,:,IDIM)=CHIDENMAT(IPAIR)%MAT(:,j,:,IDIM) &
!!$     &                                        +cchi(i1:i2,:,IBh)*FUNC(j,IDIM)
!!$              ENDDO
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$      return
!!$      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine localize_expand(xk,ndim,nbh,npro,proj,nchi,cchi)
!     **************************************************************************
!     ** converts the projections into coefficients for local orbitals        **
!     **  |psi(k)> approx sum_{i,j} |chi_i> dlist_{i,j} <ptilde_j|psitilde(k)>**
!     **************************************************************************
      implicit none
      real(8)   ,intent(in) :: xk(3)
      integer(4),intent(in) :: ndim
      integer(4),intent(in) :: nbh
      integer(4),intent(in) :: npro
      integer(4),intent(in) :: nchi
      complex(8),intent(in) :: proj(ndim,nbh,npro)      
      complex(8)            :: propsi(npro,ndim,nbh)      
      complex(8)            :: cchi(nchi,ndim,nbh)
      complex(8)            :: dofk(npro,nchi)
      integer(4)            :: i,j,k
!     **************************************************************************
!
!     ==========================================================================
!     == transpose projections so that they fit into the library routine      ==
!     ==========================================================================
      do i=1,npro
        do j=1,nbh
          do k=1,ndim
            propsi(i,j,k)=proj(k,j,i)
          enddo
        enddo
      enddo
!
!     ==========================================================================
!     == determine expansion coefficients for the given k-point               ==
!     ==========================================================================
      call localize_dofk(xk,nchi,npro,dofk)
!
!     ==========================================================================
!     == solve equation system to obtain expansion coefficients in corr. orb  ==
!     ==========================================================================
      call LIB$MATRIXSOLVEC8(nchi,npro,ndim*nbh,dofk,cchi,propsi)
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine localize_dofk(xk,nchi_,nphi_,dofk)
!     **************************************************************************
!     ** determine the expansion coefficients of a bloch wave of partial waves**
!     ** into localized orbitals                                              **
!     **  |chi>=sum |phi>d
!     **************************************************************************
      use localize_module
      implicit none
      real(8)   ,intent(in) :: xk(3) ! k in relative coordinates k=gbas*xk
      integer(4),intent(in) :: nchi_ 
      integer(4),intent(in) :: nphi_ 
      complex(8),intent(out):: dofk(nphi_,nchi_)
      integer(4)            :: ipos(nat+1)
      integer(4)            :: jpos(nat+1)
      complex(8),parameter  :: ci=(0.d0,1.d0)
      real(8)               :: kt
      complex(8)            :: eikt
      real(8)               :: pi
      real(8)               :: twopi
      integer(4)            :: it(3)
      integer(4)            :: iat,ln,l,ipair,iat1,iat2,i1,i2,j1,j2
!     **************************************************************************
      pi=atan(1.d0)
      twopi=2.d0*pi
      if(nchi_.ne.nchi.or.nphi_.ne.nphi) then
        call error$msg('array bounds inconsistent') 
        call error$msg('array bounds inconsistent') 
         call error$stop('localize_dofk')
     end if
!
!     ==========================================================================
!     == determine helper array for the placement of submatrices              ==
!     ==========================================================================
      ipos(1)=1
      jpos(1)=1
      do iat=1,nat
        ipos(iat+1)=ipos(iat)
        do ln=1,lnxchi(iat)
          l=loxchi(ln,iat)
          ipos(iat+1)=ipos(iat+1)+(2*l+1)
        enddo
        jpos(iat+1)=jpos(iat)
        do ln=1,lnxphi(iat)
          l=loxphi(ln,iat)
          jpos(iat+1)=jpos(iat+1)+(2*l+1)
        enddo
      enddo         
!
!     ==========================================================================
!     == calculate dofk                                                       ==
!     ==========================================================================
      dofk(:,:)=(0.d0,0.d0)
      do ipair=1,npaird
        iat1=dlist(ipair)%iat1
        iat2=dlist(ipair)%iat2
        it(:)=dlist(ipair)%it(:)
        i1=ipos(iat1)
        i2=ipos(iat1+1)-1
        j1=jpos(iat2)
        j2=jpos(iat2+1)-1
        kt=xk(1)*real(it(1),kind=8) &
     &    +xk(2)*real(it(2),kind=8) &
     &    +xk(3)*real(it(3),kind=8) 
        kt=twopi*kt
        eikt=cmplx(cos(kt),sin(kt))
        dofk(i1:i2,j1:j2)=dofk(i1:i2,j1:j2)+dlist(ipair)%mat(:,:)*eikt
      enddo
      return
      end      
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LOCALIZE$nlexchange()
!     **************************************************************************
!     ** determines the approximate offsite exchange energy                   **
!     **************************************************************************
      USE LOCALIZE_MODULE
      IMPLICIT NONE
      real(8)    :: ecorr
      real(8)    :: d(3),dis
      real(8)    :: r(3,nat)
      real(8)    :: rbas(3,3)
      real(8)    :: svar
      integer(4) :: iat1,iat2
      integer(4) :: it(3)
      integer(4) :: lmnx1,lmnx2
      integer(4) :: i,idim1,idim2,lmn1,lmn2
!     **************************************************************************
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
      ECORR=0.D0
      DO I=1,NNB
        IAT1=CHIDENMAT(I)%IAT1
        IAT2=CHIDENMAT(I)%IAT2
        IT=CHIDENMAT(I)%IT
        LMNX1=SUM(2*LOXCHI(:LNXCHI(IAT1),IAT1)+1)
        LMNX2=SUM(2*LOXCHI(:LNXCHI(IAT2),IAT2)+1)
        D(:)=R(:,IAT2)+MATMUL(RBAS,REAL(IT,KIND=8))-R(:,IAT1)
        DIS=SQRT(SUM(D(:)**2))
        SVAR=0.D0
        DO IDIM1=1,NDIM
          DO IDIM2=1,NDIM
            DO LMN2=1,LMNX2
              DO LMN1=1,LMNX1
                SVAR=SVAR+CHIDENMAT(I)%MAT1(LMN1,LMN2,IDIM1,IDIM2)**2
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        ECORR=ECORR-0.5D0*SVAR/DIS
      ENDDO
PRINT*,'ECORR ',ECORR
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LOCALIZE$switchCHIDENMAT(on)
!     **************************************************************************
!     **************************************************************************
      USE LOCALIZE_MODULE
      IMPLICIT NONE
      logical(4),intent(in) :: on
!     **************************************************************************
      if(on) then
        call LOCALIZE$NBLIST()
        call localize_inichidenmat
      else
        call localize_CLEANCHIDENMAT()
      end if
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LOCALIZE$reportCHIDENMAT(nfil)
!     **************************************************************************
!     **************************************************************************
      USE LOCALIZE_MODULE
      IMPLICIT NONE
      integer(4),intent(in) :: nfil
      real(8)               :: rbas(3,3)
      real(8)               :: r(3,nat)
      real(8)               :: d(3)
      real(8)               :: dis
      integer(4)            :: iat1,iat2,lmnx1,lmnx2
      integer(4)            :: i,lmn1
      integer(4)            :: it(3)
!     **************************************************************************
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
!
      do i=1,nnb
        iat1=chidenmat(i)%iat1
        iat2=chidenmat(i)%iat2
        it=chidenmat(i)%it
        lmnx1=sum(2*loxchi(:lnxchi(iat1),iat1)+1)
        lmnx2=sum(2*loxchi(:lnxchi(iat2),iat2)+1)
        D(:)=R(:,IAT2)+MATMUL(RBAS,REAL(IT,KIND=8))-R(:,IAT1)
        Dis=sqrt(sum(d(:)**2))
        write(nfil,fmt='(82("="),t2,"iat1=",i3," iat2=",i3," dis=",f7.3," dir=",3f7.3)')iat1,iat2,dis,d(:)/dis
        do lmn1=1,lmnx1
          write(nfil,fmt='(20("(",2f7.3,") "))')chidenmat(i)%mat1(lmn1,:,1,1)
!          write(nfil,fmt='(20f10.3)')chidenmat(i)%mat1(lmn1,:,1,1)+chidenmat(i)%mat1(lmn1,:,2,2)
        enddo
      enddo
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LOCALIZE$CHIDENMAT(XK,TINV,NB,OCC,NDIM_,NBH,NPRO,PROJ)
!     **************************************************************************
!     ** CONVERT THE PROJECTIONS <PTILDE|PSITILDE(K)> INTO PROJECTIONS ONTO   **
!     ** LOCAL ORBITALS
!     **************************************************************************
      USE LOCALIZE_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: XK(3)
      LOGICAL(4),INTENT(IN) :: TINV
      INTEGER(4),INTENT(IN) :: NDIM_
      INTEGER(4),INTENT(IN) :: NBH
      INTEGER(4),INTENT(IN) :: NB
      REAL(8)   ,INTENT(IN) :: OCC(NB)
      INTEGER(4),INTENT(IN) :: NPRO
      COMPLEX(8),INTENT(IN) :: PROJ(NDIM_,NBH,NPRO)      
      COMPLEX(8),ALLOCATABLE:: PROJCHI(:,:,:)
      INTEGER(4)            :: IPAIR
      INTEGER(4)            :: IAT1,IAT2
      INTEGER(4)            :: I1,I2,j1,j2
      INTEGER(4)            :: IT(3)
      REAL(8)               :: KT
      COMPLEX(8)            :: EIKT
      REAL(8)               :: PI,TWOPI
      REAL(8)               :: SVAR1,SVAR2
      INTEGER(4)            :: IB,IBH,IDIM,jdim,J,iat
      INTEGER(4)            :: lmnxx
      COMPLEX(8),ALLOCATABLE:: FUNC(:,:)
!     **************************************************************************
      CALL LOCALIZE_INITIALIZE()
      PI=4.D0*ATAN(1.D0)
      TWOPI=2.D0*PI
!
      if(ndim_.ne.ndim) then 
        call error$stop('LOCALIZE$CHIDENMAT')
      end if
!
!     ==========================================================================
!     ==  TRANSFORM PROJECTIONS FROM PARTIAL WAVES TO LOCAL ORBITALS          ==
!     ==========================================================================
      ALLOCATE(PROJCHI(NDIM,NBH,NCHI))
      CALL LOCALIZE_MAPPROJTOCHI(NDIM,NBH,NPRO,PROJ,NCHI,PROJCHI)
!
!     ==========================================================================
!     ==  ADD UP DENSITY MATRIX                                               ==
!     ==========================================================================
      LMNXX=0
      DO IAT=1,NAT
        LMNXX=MAX(LMNXX,IPOSCHI(IAT+1)-IPOSCHI(IAT))
      ENDDO
      ALLOCATE(FUNC(LMNXX,NDIM))
      DO IPAIR=1,Nnb
        IAT1=CHIDENMAT(IPAIR)%IAT1
        IAT2=CHIDENMAT(IPAIR)%IAT2
        IT(:)=CHIDENMAT(IPAIR)%IT(:)
        I1=IPOSCHI(IAT1)
        I2=IPOSCHI(IAT1+1)-1
        j1=IPOSCHI(IAT2)
        j2=IPOSCHI(IAT2+1)-1
        KT=XK(1)*REAL(IT(1),KIND=8) &
     &    +XK(2)*REAL(IT(2),KIND=8) &
     &    +XK(3)*REAL(IT(3),KIND=8) 
        KT=TWOPI*KT
        EIKT=CMPLX(COS(KT),SIN(KT))
        IB=0
        DO IBH=1,NBH
          IF(TINV) THEN
            SVAR1=0.5D0*(OCC(2*IBh-1)+OCC(2*IBh))
            SVAR2=0.5D0*(OCC(2*IBh-1)-OCC(2*IBh))
          ELSE
            SVAR1=OCC(IBh)
          END IF
!
          DO j=j1,j2
            DO JDIM=1,NDIM
              IF(TINV) THEN
!               == FUNC=(F1+F2)/2<PSI1-I*PSI2|P>+(F1-F2)/2<PSI1+I*PSI2|P>
                FUNC(j-j1+1,JDIM)=SVAR1*CONJG(PROJCHI(JDIM,ibh,j)) &
     &                           +SVAR2*PROJCHI(JDIM,ibh,j)
              ELSE 
                FUNC(j-j1+1,JDIM)=SVAR1*CONJG(PROJCHI(JDIM,ibh,j))
              END IF
            ENDDO
          ENDDO
          FUNC(:,:)=FUNC(:,:)*EIKT
          DO j=1,j2-j1+1
            DO IDIM=1,NDIM
              DO JDIM=1,NDIM
                IF(TINV) THEN
                  CHIDENMAT(IPAIR)%MAT1(:,J,IDIM,JDIM)=CHIDENMAT(IPAIR)%MAT1(:,J,IDIM,JDIM) &
     &                                        +REAL(PROJCHI(IDIM,IBH,I1:I2)*FUNC(J,JDIM))
                ELSE
                  CHIDENMAT(IPAIR)%MAT1(:,J,IDIM,JDIM)=CHIDENMAT(IPAIR)%MAT1(:,J,IDIM,JDIM) &
     &                                        +PROJCHI(IDIM,IBH,I1:I2)*FUNC(J,JDIM)
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LOCALIZE_MAPPROJTOCHI(NDIM,NBH,NPRO,PROJ,NCHI,PROJchi)
!     **************************************************************************
!     ** CONVERT THE PROJECTIONS <PTILDE|PSITILDE(K)> INTO PROJECTIONS ONTO   **
!     ** LOCAL ORBITALS
!     **************************************************************************
      IMPLICIT NONE
      TYPE SPECIES_LOCALTYPE
        INTEGER(4)         :: LNX
        INTEGER(4)         :: LNXCHI
        INTEGER(4),POINTER :: LOX(:)
        INTEGER(4),POINTER :: LOXCHI(:)
        real(8)   ,POINTER :: BMAT(:,:)
      END TYPE SPECIES_LOCALTYPE
      INTEGER(4),INTENT(IN) :: NDIM
      INTEGER(4),INTENT(IN) :: NBH
      INTEGER(4),INTENT(IN) :: NPRO
      INTEGER(4),INTENT(IN) :: NCHI
      COMPLEX(8),INTENT(IN) :: PROJ(NDIM,NBH,NPRO)      
      COMPLEX(8),INTENT(OUT):: projchi(NDIM,NBH,NCHI)
      INTEGER(4),ALLOCATABLE:: ISPECIES(:)
      TYPE(SPECIES_LOCALTYPE),ALLOCATABLE :: SPECIES(:)
      INTEGER(4)            :: nsp   ! #(atom types)
      INTEGER(4)            :: nat   ! #(atoms)
      INTEGER(4)            :: isp,iat,ipro,ichi
      INTEGER(4)            :: LNX,LNXCHI
      INTEGER(4)            :: L     ! angular momentu,
      INTEGER(4)            :: LNchi,ln,lmnchi,lmn
      INTEGER(4)            :: i1chi,i2chi
      INTEGER(4)            :: i1pro,i2pro
      real(8)               :: svar
!     **************************************************************************
!
!     ==========================================================================
!     ==  COLLECT SPECIES INFOMATION                                          ==
!     ==========================================================================
      CALL SETUP$NSPECIES(NSP)
      ALLOCATE(SPECIES(NSP))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        CALL SETUP$GETI4('LNXCHI',LNXCHI)
        SPECIES(ISP)%LNX=LNX
        SPECIES(ISP)%LNXCHI=LNXCHI
        ALLOCATE(SPECIES(ISP)%LOX(LNX))
        ALLOCATE(SPECIES(ISP)%LOXCHI(LNXCHI))
        CALL SETUP$GETI4A('LOX',LNX,SPECIES(ISP)%LOX)
        CALL SETUP$GETI4A('LOXCHI',LNXCHI,SPECIES(ISP)%LOXCHI)
        ALLOCATE(SPECIES(ISP)%BMAT(LNXCHI,LNX))
        CALL SETUP$GETr8A('BMATCHI',LNXCHI*LNX,SPECIES(ISP)%BMAT)
      ENDDO
!
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
!
!     ==========================================================================
!     ==  TRANSFORM PROJECTIONS INTO THOSE ON LOCAL ORBITALS                  ==
!     ==========================================================================
      projchi(:,:,:)=(0.d0,0.d0)
      IPRO=0
      ICHI=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        LNX=SPECIES(ISP)%LNX
        LNXCHI=SPECIES(ISP)%LNXCHI
        LMNCHI=0
        DO LNCHI=1,LNXCHI
          L=SPECIEs(ISP)%LOXCHI(LNCHI)
          LMN=0
          DO LN=1,LNX
            IF(SPECIES(ISP)%LOX(LN).EQ.L) THEN
              I1CHI=ICHI+LMNCHI+1
              I2CHI=I1CHI-1+2*L+1
              I1PRO=IPRO+LMN+1
              I2PRO=I1PRO-1+2*L+1
              SVAR=species(isp)%BMAT(LNCHI,LN)
              projCHI(:,:,I1CHI:I2CHI)=projCHI(:,:,I1CHI:I2CHI) &
      &                               +SVAR*PROJ(:,:,I1PRO:I2PRO)
            END IF
            LMN=LMN+2*SPECIES(ISP)%LOX(LN)+1
          ENDDO
          LMNCHI=LMNCHI+2*L+1
        ENDDO        
        ICHI=ICHI+LMnchi
        IPRO=IPRO+LMN
      ENDDO
!
!     ==========================================================================
!     ==   CLEAN UP SPECIES                                                   ==
!     ==========================================================================
      DO ISP=1,NSP
        DEALLOCATE(SPECIES(ISP)%LOX)
        DEALLOCATE(SPECIES(ISP)%LOXCHI)
        DEALLOCATE(SPECIES(ISP)%BMAT)
      ENDDO
      DEALLOCATE(SPECIES)
      RETURN
      END
