!........1.........2.........3.........4.........5.........6.........7.........8
module localize_module
!*******************************************************************************
!*******************************************************************************
!****  dlist contains the expansion coefficients of the correlated orbitals ****
!****        centered at iat1 in partial waves centered at iat2 translated  ****
!****        by the translations it
!*******************************************************************************
!*******************************************************************************
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
end type matlist_type
integer(4)        ,save        :: npaird
real(8)           ,parameter   :: rmax=6.d0
integer(4)        ,save        :: nnbx=0
type(nblist_type) ,allocatable :: nblist(:) ! neighbor list
type(matlist_type),allocatable :: dlist(:) ! expansion coefficients of 
                                           ! local orbitals into partial waves
type(matlist_type),allocatable :: chidenmat(:) ! expansion coefficients of 
                                           ! local orbitals into partial waves
integer(4),save        :: nat
integer(4),save        :: nphi
integer(4),save        :: nchi
real(8)   ,allocatable :: lnxchi(:)
real(8)   ,allocatable :: loxchi(:,:)
real(8)   ,allocatable :: lnxphi(:)
real(8)   ,allocatable :: loxphi(:,:)
end module localize_module
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine localize_nblist()
!     **************************************************************************
!     **************************************************************************
      use localize_module
      implicit none
      real(8)                :: rbas(3,3)    ! Lattice translations
      real(8)                :: rbasinv(3,3) ! inverse of rbas
      integer(4)             :: nat          ! #(atoms)
      real(8)   ,allocatable :: r(:,:)       ! atomic positions
      real(8)                :: xr(3)
      real(8)                :: x0,y0,z0
      integer(4),allocatable :: itr(3,nat)
      integer(4)             :: itvec(3)
      real(8)                :: d(3)         ! distance vector
      real(8)                :: d2           ! squared distance
      real(8)                :: rmax2        ! squared maximum distance
      integer(4)             :: iat,iat1,iat2,it1,it2,it3
      integer(4)             :: min1,max1,min2,max2,min3,max3
!     **************************************************************************
      CALL CELL$GETR8A('T0',9,RBAS)
      call lib$invertr8(3,rbas,rbasinv)
      CALL ATOMLIST$NATOM(NAT)
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
        CALL BOXSPH(RBOX,X0,Y0,Z0,Rmax,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
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
      IF(NNB.LT.NNBX) THEN
        IF(NNBX.GT.0)DEALLOCATE(NBLIST)
        ALLOCATE(NBLIST(NNB))
        NNBX=NNB
        GOTO 1000
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine localize_chidenmat(xk,tinv,ndim,nbh,nb,OCC,npro_,proj,nchi_,cchi)
!     **************************************************************************
!     **************************************************************************
      implicit none
      real(8)   ,intent(in) :: xk(3)
      logical(4),intent(in) :: tinv
      integer(4),intent(in) :: ndim
      integer(4),intent(in) :: nbh
      integer(4),intent(in) :: nb
      REAL(8)   ,INTENT(IN) :: OCC(NB)
      integer(4),intent(in) :: npro_
      integer(4),intent(in) :: nchi_
      complex(8),intent(in) :: proj(ndim,nbh,npro_)      
      complex(8)            :: cchi(nchi_,ndim,nbh)
!     **************************************************************************
!
!     ==========================================================================
!     ==  CHECK IF SUPERWAVEFUNCTIONS ARE USED AND IF #(BANDS) CORRECT        ==
!     ==========================================================================
      IF(TINV) THEN
        IF(NBH.NE.(NB+1)/2) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
          CALL ERROR$STOP('WAVES_DENMAT')
        END IF
      ELSE 
        IF(NBH.NE.NB) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
          CALL ERROR$STOP('WAVES_DENMAT')
        END IF
      END IF
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      do ipair=1,npaird
        iat1=chidenmat(ipair)%iat1
        iat2=chidenmat(ipair)%iat2
        it(:)=chidenmat(ipair)%it(:)
        i1=ipos(iat1)
        i2=ipos(iat1+1)
        j1=jpos(iat2)
        j2=jpos(iat2+1)
        kt=xk(1)*real(it(1),kind=8) &
     &    +xk(2)*real(it(2),kind=8) &
     &    +xk(3)*real(it(3),kind=8) 
        kt=twopi*kt
        eikt=(cos(kt),-sin(kt))
        ib=0
        do ibh=1,nbh
          IF(TINV) THEN
            SVAR1=0.5D0*(OCC(2*IB-1)+OCC(2*IB))
            SVAR2=0.5D0*(OCC(2*IB-1)-OCC(2*IB))
          ELSE
            SVAR1=OCC(IB)
          END IF
!
          do j=1,j2-j1+1
            DO IDIM1=1,NDIM
              IF(TINV) THEN
!               == FUNC=(F1+F2)/2<PSI1-I*PSI2|P>+(F1-F2)/2<PSI1+I*PSI2|P>
                FUNC(j,IDIM1)=SVAR1*CONJG(cchi(IDIM1,j1-1+j,ibh)) &
     &                       +SVAR2*cchi(IDIM1,j1-1+j,ibh)
              ELSE 
                FUNC(j,IDIM1)=SVAR1*CONJG(cchiI(IDIM1,j1-1+j,IB))
              END IF
            ENDDO
          ENDDO
          FUNC(:)=FUNC(:)*EIKT
          DO j=1,j2-j1+1
            DO IDIM=1,NDIM
              CHIDENMAT(IPAIR)%MAT(:,j,:,IDIM)=CHIDENMAT(IPAIR)%MAT(:,j,:,IDIM)
     &                                        +cchi(i1:i2,:,IBh)*FUNC(j,IDIM)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine localize_expand(xk,ndim,nbh,npro_,proj,nchi_,cchi)
!     **************************************************************************
!     **************************************************************************
      implicit none
      real(8)   ,intent(in) :: xk(3)
      integer(4),intent(in) :: ndim
      integer(4),intent(in) :: nbh
      integer(4),intent(in) :: npro_
      integer(4),intent(in) :: nchi_
      complex(8),intent(in) :: proj(ndim,nbh,npro_)      
      complex(8)            :: propsi(npro,ndim,nbh)      
      complex(8)            :: cchi(nchi_,ndim,nbh)
      real(8)               :: dofk(npro,nchi)
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
      call localize_dofk(xk,nchi,nphi,dofk)
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
      use localize module
      implicit none
      real(8)   ,intent(in) :: xk(3) ! k in relative coordinates k=gbas*xk
      integer(4),intent(in) :: nchi 
      integer(4),intent(in) :: nphi 
      complex(8),intent(out):: dofk(nphi,nchi)
      integer(4)            :: ipos(nat+1)
      integer(4)            :: jpos(nat+1)
      complex(8),parameter  :: ci=(0.d0,1.d0)
      real(8)               :: kt
      complex(8)            :: eikt
      real(8)               :: pi
      real(8)               :: twopi
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
      ipos1(1)=1
      jpos1(1)=1
      do iat=1,nat
        do ln=1,lnxchi(iat)
          l=loxchi(ln,iat)
          ipos2(iat)=ipos1(iat)-1+2*l+1
          ipos1(iat+1)=ipos2(iat)+1
        enddo
        do ln=1,lnxphi(iat)
          l=loxphi(ln,iat)
          jpos2(iat)=jpos1(iat)-1+2*l+1
          jpos1(iat+1)=jpos2(iat)+1
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
        i2=ipos(iat1+1)
        j1=jpos(iat2)
        j2=jpos(iat2+1)
        kt=xk(1)*real(it(1),kind=8) &
     &    +xk(2)*real(it(2),kind=8) &
     &    +xk(3)*real(it(3),kind=8) 
        kt=twopi*kt
        eikt=(cos(kt),sin(kt))
        dofk(i1:i2,j1,j2)=dofk(i1:i2,j1,j2)+dlist(ipair)%mat(:,:)*eikt
      enddo
      return
      end      
