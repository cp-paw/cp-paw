      program main
      implicit none
      integer(4),parameter :: nb=1
      integer(4),parameter :: ngrid=5
      integer(4),parameter :: ndim=1
      real(8)   ,parameter :: t=1.d0
      real(8)   ,parameter :: u=1.d0
      real(8)              :: pi
      complex(8)           :: hamilton(nb,nb)
      complex(8)           :: proj(ndim,ngrid,nb)
      real(8)              :: k
      integer(4)           :: ib,j,idim
      real(8)              :: svar
!***************************************************************
      pi=4.d0*datan(1.d0)
      call dmft__setl4('on',.true.)
      call dmft__seti4('npro',ngrid)
      call dmft__seti4('nspin',1)
      call dmft__seti4('nkpt',1)
      call dmft__seti4('nat',1)
      call dmft__seti4('ndim',1)
      call dmft__seti4('nb',nb)
      call dmft__seti4a('natorb',1,(/1/))
!
!     == evaluate states =======================================
      hamilton(:,:)=(0.d0,0.d0)
      do ib=1,nb
!       == wave function
        k=pi*real(ib)/real(ngrid+1)
        do j=1,ngrid
          proj(1,j,ib)=sin(k*real(j))
        enddo
!       ==  normalize
        svar=0.d0
        do idim=1,ndim
          svar=svar+dot_product(proj(idim,:,ib),proj(idim,:,ib))
        enddo
        proj(:,:,ib)=proj(:,:,ib)/sqrt(svar)
!       == energy expectation value
        svar=0.d0
        do j=1,ngrid-1
          do idim=1,ndim
            svar=svar-2.d0*t*proj(idim,j,ib)*proj(idim,j+1,ib)
          enddo
        enddo
        hamilton(ib,ib)=svar
      enddo
!
!     ==========================================================
      do ib=1,nb
        write(*,fmt='("e=",f5.1," p=",10f5.2)') &
     &              real(hamilton(ib,ib)),real(proj(:,:,ib))
      enddo
!
!     == set projectors            ======================================
      call dmft__setproj(1,1,1,1,nb,ngrid,proj,hamilton)
!
!     == evaluate weiss field      ======================================
      call dmft__etot
      stop 'finished'
      end

!************************************************************************
!  demonstrator code for dmft
!  code written for a single site 
!************************************************************************
module dmft_module
type atomblock_type
  integer(4)         :: norb
  integer(4)         :: norbs
  real(8)   ,pointer :: orbpro(:,:)
  complex(8),pointer :: hybridf(:,:,:,:)  ! (norbs,norbs,NOMEGA,NSPIN)
  complex(8),pointer :: dsigma(:,:,:,:) ! (norbs,norbs,NOMEGA,NSPIN)
  complex(8),pointer :: sigma0(:,:,:)   ! (norbs,norbs,Nspin)
end type atomblock_type
logical(4)             :: ton=.false.
logical(4)             :: tini=.false.
integer(4)             :: nomega=0
real(8)   ,allocatable :: omega(:)
integer(4)             :: nat=0
integer(4)             :: nkpt=0
integer(4)             :: nb=0
integer(4)             :: nspin=0
integer(4)             :: ndim=0
integer(4)             :: norb=0
integer(4)             :: norbs=0
integer(4)             :: npro=0
complex(8),allocatable :: u(:,:,:,:)        !(nb,nb,nspin,nkpt)
real(8)   ,allocatable :: epsilon(:,:,:)    !(nb,nspin,nkpt)
real(8)                :: mu                ! chemical potential
complex(8),allocatable :: orbcoeff(:,:,:,:) !(norbs,nb,nspin,nkpt)
complex(8),allocatable :: green0(:,:,:,:,:) !(norbs,norbs,NOMEGA,nspin,ikpt)
complex(8),allocatable :: green(:,:,:,:,:)  !(norbs,norbs,nomega,NSPIN,ikpt)
complex(8),allocatable :: dmat(:,:,:,:)        !(nb,nb,nspin,ikpt)
type(atomblock_type),pointer :: atomblock(:)         ! (nat)
end module dmft_module
!
!      .................................................................
       subroutine dmft_initialize
!      *****************************************************************
!      **  evaluates the coefficients of the local orbitals           **
!      *****************************************************************
       use dmft_module
       implicit none
       real(8)    :: omegamin,omegamax
       real(8)    :: svar
       integer(4) :: i,iat
       integer(4) :: norb1,norbs1
!      *****************************************************************
       tini=.true.
       if(nomega.eq.0) then        
         !call error__msg('nomega not defined')
         !call error__stop('dmft__seti4')
       end if
       if(nat.eq.0) then        
         !call error__msg('nomega not defined')
         !call error__stop('dmft__seti4')
STOP 'DMFT_INITIALIZE,nat'
       end if
       if(NSPIN.eq.0) then        
         !call error__msg('nSPIN not defined')
         !call error__stop('dmft__seti4')
STOP 'DMFT_INITIALIZE,nspin'
       end if
       if(nb.eq.0) then        
         !call error__msg('nb not defined')
         !call error__stop('dmft__seti4')
STOP 'DMFT_INITIALIZE,nb'
       end if
!
!      =================================================================
!      ==  evaluate omega grid                                        ==
!      =================================================================
       nomega=1000
       omegamin=-5.d0
       omegamax=5.d0
       allocate(omega(nomega))
       svar=(omegamax-omegamin)/real(nomega-1,kind=8)
       do i=1,nomega
         omega(i)=omegamin+real(i-1,kind=8)*svar
       enddo
!
!      =================================================================
!      ==  define orbital projetors                                   ==
!      =================================================================
       allocate(atomblock(nat))
       do iat=1,nat
         atomblock(iat)%norb=0
         atomblock(iat)%norbs=0
       enddo
       iat=1
       norb=1
       atomblock(iat)%norb=norb
       allocate(atomblock(iat)%orbpro(npro,norb))
       atomblock(iat)%orbpro(:,:)=0.d0
       atomblock(iat)%orbpro(2,1)=1.d0
!
!      =================================================================
!      ==  determine total number of localized orbitals
!      =================================================================
       norb=0
       norbs=0
       do iat=1,nat
         norb1=atomblock(iat)%norb
         norbs1=norb1*ndim
         norb=norb+norb1
         norbs=norbs+norbs1
         atomblock(iat)%norbs=norb1*ndim
         allocate(atomblock(iat)%sigma0(norbs1,norbs1,nspin))
         atomblock(iat)%sigma0(:,:,:)=(0.d0,0.d0)
         allocate(atomblock(iat)%dsigma(norbs1,norbs1,nspin,nomega))
         atomblock(iat)%dsigma(:,:,:,:)=(0.d0,0.d0)
         allocate(atomblock(iat)%hybridf(norbs1,norbs1,nspin,nomega))
         atomblock(iat)%hybridf(:,:,:,:)=(0.d0,0.d0)
       enddo
       return
       end
!
!      .................................................................
       subroutine dmft__seti4(id,val)
       use dmft_module
       implicit none
       character(*),intent(in) :: id
       integer(4)  ,intent(in) :: val
!      *****************************************************************
       if(id.eq.'nomega') then
         nomega=val
         !call error__msg('nomega is still hard wired; do not set')
         !call error__stop('dmft__seti4')
stop 'error setting nomega'
       else if(id.eq.'npro') then
         npro=val
       else if(id.eq.'nspin') then
         nspin=val
       else if(id.eq.'ndim') then
         ndim=val
       else if(id.eq.'nkpt') then
         nkpt=val
       else if(id.eq.'nb') then
         nb=val
       else if(id.eq.'nat') then
         if(nat.ne.0.and.nat.ne.val) then
           !call error__msg('nat has been specified before and cannot be changed')
           !call error__stop('dmft__seti4')
stop 'error setting nat'
         end if
         nat=val
       else
         !call error__msg('id not recognized')
         !call error__chval('id',id)
         !call error__stop('dmft__seti4')
stop 'error id not recognized in seti4'
       end if
       return
       end
!
!      .................................................................
       subroutine dmft__seti4a(id,len,val)
       use dmft_module
       implicit none
       character(*),intent(in) :: id
       integer(4)  ,intent(in) :: len
       integer(4)  ,intent(in) :: val(len)
       integer(4)              :: i
!      *****************************************************************
       if(id.eq.'natorb') then
         if(.not.associated(atomblock).and.nat.eq.0) nat=len
         if(len.ne.nat) then
           !call error__msg('size inconsistent')
           !call error__chval('id',id)
           !call error__stop('dmft__seti4a')
stop 'error setting natorb'
         end if
         if(.not.associated(atomblock))allocate(atomblock(nat))
         do i=1,nat
           atomblock(i)%norb=val(i)
         enddo
       else
         !call error__msg('id not recognized')
         !call error__chval('id',id)
         !call error__stop('dmft__seti4a')
stop 'error id not recognized in seti4a'
       end if
       return
       end
!
!      .................................................................
       subroutine dmft__setl4(id,val)
       use dmft_module
       implicit none
       character(*),intent(in) :: id
       logical(4)  ,intent(in) :: val
!      *****************************************************************
       if(id.eq.'on') then
         ton=val
       else
         !call error__msg('id not recognized')
         !call error__chval('id',id)
         !call error__stop('dmft__seti4')
stop 'error id not recognized in setl4'
       end if
       return
       end
!
!      .................................................................
       subroutine dmft__setproj(ispin,ikpt,nspin_,ndim_,nb_,npro_,proj,hamilton)
!      *****************************************************************
!      **  projects onto wave functions onto local orbitals           **
!      **  and adds result to weiss function                          **
!      **  nspin=1 non-spin polarized ndim=1
!      **  nspin=2 spin polarized     ndim=1
!      **  nspin=3 non-collinear      ndim=2
!      *****************************************************************
       use dmft_module
       implicit none
       integer(4),intent(in) :: ispin      ! spin direction (1,2)
       integer(4),intent(in) :: ikpt       ! k-point index
       integer(4),intent(in) :: nspin_     ! (1,2,3) 
       integer(4),intent(in) :: ndim_      ! #(spinor components)
       integer(4),intent(in) :: nb_        ! #(bands)
       integer(4),intent(in) :: npro_      ! #(partial waves)
       complex(8),intent(in) :: proj(ndim,npro,nb_)
       complex(8),intent(in) :: hamilton(nb_,nb_)
       real(8)   ,allocatable:: orbpro(:,:)     ! (npro,norb)
       integer(4)            :: iorb,idim,ib,io,io1,iat
!      *****************************************************************
       if(.not.ton) return
       if(.not.tini) call dmft_initialize
       if(npro_.ne.npro) then
         !call error__msg('npro inconsistent')
         !call error__stop('dmft__seproj')
       end if
!
!      =================================================================
!      ==  project wave function onto local orbitals                  ==
!      =================================================================
       allocate(orbcoeff(norbs,nb,NSPIN,NKPT))
       allocate(orbpro(npro,norb))
       io=0
       do iat=1,nat
         do io1=1,atomblock(iat)%norb
           io=io+1
           orbpro(:,io)=atomblock(iat)%orbpro(:,io1)
         enddo
       enddo
       call dmft_projections(nb,npro,ndim,norb,proj,orbpro &
      &                     ,orbcoeff(:,:,ispin,ikpt))
       deallocate(orbpro)
!
!      =================================================================
!      ==  transform to eigenstates                                   ==
!      =================================================================
       if(.not.allocated(epsilon))allocate(epsilon(nb,ispin,ikpt))
       if(.not.allocated(u))allocate(u(nb,nb,ispin,ikpt))
!      call lib__diagc8(nb,hamilton,epsilon(:,ispin,ikpt),u(:,:,ispin,ikpt)
u=(0.d0,0.d0)
do ib=1,nb
  epsilon(ib,ispin,ikpt)=real(hamilton(ib,ib))
  u(ib,ib,ispin,ikpt)=(1.d0,0.d0)
enddo
       orbcoeff(:,:,ispin,ikpt)=matmul(orbcoeff(:,:,ispin,ikpt),u(:,:,ispin,ikpt))
       return
       end
!
!      .................................................................
       subroutine dmft_projections(nb,npro,ndim,norb,proj,orbpro &
      &                           ,orbcoeff)
!      *****************************************************************
!      **  evaluates the coefficients of the local orbitals           **
!      **  for each one-particle state                                **
!      **                                                             **
!      **      |psi> = |chi>*orbpro*<ptilde|psitilde>                 **
!      *****************************************************************
       implicit none
       integer(4),intent(in) :: nb         ! #(bands)
       integer(4),intent(in) :: npro       ! #(partial waves)
       integer(4),intent(in) :: ndim       ! 
       integer(4),intent(in) :: norb       ! #(localized orbitals)
       complex(8),intent(in) :: proj(ndim,npro,nb)  ! <p|psitilde>
       real(8)   ,intent(in) :: orbpro(npro,norb)    !
       complex(8),intent(out):: orbcoeff(norb*ndim,nb)
       integer(4)            :: iorb,iorbs,ipro,ib,idim
       real(8)               :: svar
!      *****************************************************************
       orbcoeff(:,:)=(0.d0,0.d0)
       do iorb=1,norb
         do ipro=1,npro
           svar=orbpro(ipro,iorb)
           if(svar.eq.0.d0) cycle
           do ib=1,nb
             do idim=1,ndim
               IORBS=IORB*(NDIM-1)+IDIM
               orbcoeff(iorbs,ib)=orbcoeff(iorbs,ib) &
      &                          +proj(idim,ipro,ib)*svar
             enddo
           enddo
         enddo
       enddo
       return
       end
!
!      .................................................................
       subroutine dmft__etot
!      *****************************************************************
!      **  projects onto wave functions onto local orbitals           **
!      **  and adds result to weiss function                          **
!      *****************************************************************
       use dmft_module
       implicit none
       complex(8),allocatable :: smat(:,:)
       complex(8),allocatable :: smatin(:,:)
       complex(8),allocatable :: dsigma(:,:) !(norbs,norbs,ispin,nomega)
       complex(8),allocatable :: sigma0(:,:)  !(norbs,norbs,ispin)
       complex(8),allocatable :: gonsite(:,:)  !(norbs,norbs,ispin)
       complex(8),allocatable :: dmat1(:,:)  !(norbs,norbs,ispin)
       complex(8),allocatable :: dsigmatilde(:,:) !(norbs,norbs,ispin)
       complex(8),parameter   :: ci=(0.d0,1.d0)
       integer(4)             :: iomega
       integer(4)             :: nfil
       integer(4)             :: io1,io2
       integer(4)             :: iat,ispin,ikpt
       integer(4)             :: istart,norbs1
!      *****************************************************************
       if(.not.ton) return
!
!      =================================================================
!      ==  determine chemical potential from one-particle energies    ==
!      =================================================================
 mu=0.d0
!
!      =================================================================
!      ==  evaluate greens function of the reference system           ==
!      =================================================================
       if(.not.allocated(green0))allocate(green0(norbs,norbs,nomega,nspin,nkpt))
       green0(:,:,:,:,:)=(0.d0,0.d0)
       do ikpt=1,nkpt
         do ispin=1,nspin
            call dmft_addgreen(nb,norbs,epsilon(:,ISPIN,IKPT),mu &
      &        ,orbcoeff(:,:,ISPIN,IKPT),nomega,omega,green0(:,:,:,:,ikpt))
         enddo
       enddo
!
!      =================================================================
!      ==  print greens function for test                             ==
!      =================================================================
!      call writeone(norbs,nomega,omega,green0(:,:,:,1,1))
!      stop 'after writing greens function'
!
!      =================================================================
!      ==  determine full green's function and sigmatilde
!      =================================================================
       if(.not.allocated(green)) &
      &         allocate(green(norbs,norbs,nomega,nspin,nkpt))
       green(:,:,:,:,:)=(0.d0,0.d0)
!
       allocate(dsigma(norbs,norbs))
       do ikpt=1,nkpt
         do ispin=1,nspin
           do iomega=1,nomega
             dsigma(:,:)=(0.d0,0.d0)
             io2=0
             do iat=1,nat
               io1=io2+1
               io2=io1-1+atomblock(iat)%norb*ndim
               dsigma(io1:io2,io1:io2)=atomblock(iat)%dsigma(:,:,ispin,iomega)
             enddo
             call dmft_fullgreen(norbs,green0(:,:,iomega,ispin,ikpt) &
        &                        ,dsigma(:,:),green(:,:,iomega,ispin,ikpt))
           enddo
         enddo
       enddo
       deallocate(dsigma)
!call writeone(atomblock(1)%norbs,nspin,nomega,omega,green(:,:,:,:,1))
!
!      =================================================================
!      ==  calculate hybridization function                           ==
!      =================================================================
       do iat=1,nat
         norbs1=atomblock(iat)%norbs
         allocate(atomblock(iat)%hybridf(norbs1,norbs1,iomega,nspin))
       enddo
         
       do ispin=1,nspin
         do iomega=1,nomega
           io2=0
           do iat=1,nat
             norbs1=atomblock(iat)%norbs
             if(norbs1.eq.0) cycle
             io1=io2+1
             io2=io1-1+norbs1
             allocate(gonsite(norbs1,norbs1))
             allocate(smat(norbs1,norbs1))
!            == k-sum of the onsite greens function =======================
             gonsite(:,:)=0.d0
             do ikpt=1,nkpt
               gonsite(:,:)=gonsite(:,:)+green(io1:io2,io1:io2,iomega,ispin,ikpt)
             enddo
             call lib__invertc8(norbs1,gonsite,smat)
!            == calculate INVERSE weiss field ==============================
             smat=smat+atomblock(iat)%dsigma(:,:,iomega,ispin) &
            &         +atomblock(iat)%sigma0(:,:,ispin)
!            == CALCULATE hybridizATION FUNCTION ===========================
             SMAT=-SMAT
             DO Io1=1,NORBS
               SMAT(Io1,Io1)=CI*OMEGA(IOMEGA)+MU+SMAT(io1,io1)
             ENDDO
             atomblock(iat)%hybridf(:,:,iomega,ispin)=smat(:,:)
             deallocate(smat)
             deallocate(gonsite)
           enddo
         enddo
       enddo
!
!      =================================================================
!      ==  call impurity solver                                       ==
!      =================================================================
       do iat=1,nat
         norbs1=atomblock(iat)%norbs
         call writeone(norbs1,nomega,omega,atomblock(IAT)%hybridf(:,:,:,1))
       enddo
!
!      =================================================================
!      ==  evaluate new density matrix                                ==
!      =================================================================
       allocate(dsigma(norbs,norbs))
       allocate(dsigmatilde(norbs,norbs))
       allocate(dmat(nb,nb,nspin,nkpt))
       allocate(dmat1(nb,nb))
       dmat(:,:,:,:)=0.d0
       do ikpt=1,nkpt
         do ispin=1,nspin
           do iomega=1,nomega
             dsigma(:,:)=0.d0
             io2=0
             do iat=1,nat
               io1=io2+1
               io2=io2+atomblock(iat)%norb*ndim
               dsigma(io1:io2,io1:io2)=atomblock(iat)%dsigma(:,:,iomega,ispin)
             enddo
             call dmft_dsigmatilde(norbs,dsigma,green0(:,:,iomega,ispin,ikpt) &
      &                                         ,dsigmatilde)
             call dmft_doccup(nb,norb,omega(iomega),epsilon,mu &
      &                      ,dsigmatilde,orbcoeff(:,:,ispin,ikpt),dmat1)
             dmat(:,:,ispin,ikpt)=dmat1(:,:)
           enddo
         enddo
       enddo
       deallocate(dsigma)
       deallocate(dmat1)
       deallocate(dsigmatilde)
!
!      =================================================================
!      ==  apply unitary transformation
!      =================================================================
       do ikpt=1,nkpt
         do ispin=1,nspin
            dmat(:,:,ikpt,ispin)=matmul(transpose(conjg(u(:,:,ispin,ikpt))) &
      &                     ,matmul(dmat(:,:,ispin,ikpt),u(:,:,ispin,ikpt)))       
         enddo
       enddo
       return
       end
!
!      .................................................................
       subroutine dmft_addgreen(nb,norb,epsilon,mu &
      &                        ,orbcoeff,nomega,omega,green)
!      *****************************************************************
!      **  evaluates the matrix elements of the greens function       **
!      **  and its inverse between the local orbitals                 **
!      **                                                             **
!      *****************************************************************
       implicit none
       integer(4),intent(in) :: nb         ! #(bands)
       integer(4),intent(in) :: norb       ! #(localized orbitals)
       real(8)   ,intent(in) :: epsilon(nb)
       real(8)   ,intent(in) :: mu
       complex(8),intent(in) :: orbcoeff(norb,nb)
       integer(4),intent(in) :: nomega
       real(8)   ,intent(in) :: omega(nomega)
       complex(8),intent(out):: green(norb,norb,nomega)
       complex(8),parameter  :: ci=(0.d0,1.d0)
       integer(4)            :: io1,io2
       integer(4)            :: ib,iomega
       complex(8)            :: csvar,cfac
!      *****************************************************************
print*,'epsilon',epsilon
print*,'orbcoeff ',orbcoeff
       do iomega=1,nomega
         do ib=1,nb
           cfac=1.d0/(ci*omega(iomega)+mu-epsilon(ib))
           do io2=1,norb
             do io1=1,norb
               csvar=orbcoeff(io1,ib)*conjg(orbcoeff(io2,ib))
               green(io1,io2,iomega)=green(io1,io2,iomega)+csvar*cfac
             enddo
           enddo
         enddo
       enddo
       return
       end
!
!      .....................................................................
       subroutine dmft_fullgreen(norb,g0,dsigma,g)
       implicit none
       integer(4),intent(in) :: norb
       complex(8),intent(in) :: dsigma(norb,norb)
       complex(8),intent(in) :: g0(norb,norb)
       complex(8),intent(out):: g(norb,norb)
       complex(8)            :: smat(norb,norb)
       complex(8)            :: smatinv(norb,norb)
       integer(4)            :: i
!      ********************************************************************
       smat=-matmul(g0,dsigma)         
       do i=1,norb
         smat(i,i)=(1.d0,0.d0)+smat(i,i)
       enddo
       call lib__invertc8(norb,smat,smatinv)  !
       g=matmul(smatinv,g0)           
       return
       end
!
!      .....................................................................
       subroutine dmft_dsigmatilde(norb,dsigma,g0,dsigmatilde)
       implicit none
       integer(4),intent(in) :: norb
       complex(8),intent(in) :: dsigma(norb,norb)
       complex(8),intent(in) :: g0(norb,norb)
       complex(8),intent(out):: dsigmatilde(norb,norb)
       complex(8)            :: smat(norb,norb)
       complex(8)            :: smatinv(norb,norb)
       integer(4)            :: i
!      ********************************************************************
       smat=-matmul(dsigma,g0)         
       do i=1,norb
         smat(i,i)=(1.d0,0.d0)+smat(i,i)
       enddo
       call lib__invertc8(norb,smat,smatinv)  !
       dsigmatilde=matmul(smatinv,dsigma)           !smat= d-sigmatilde
       return
       end
!
!      .....................................................................
       subroutine dmft_doccup(nb,norb,omega,epsilon,mu,dsigmatilde,proj,dmat)
       implicit none
       integer(4),intent(in) :: norb
       integer(4),intent(in) :: nb
       real(8)   ,intent(in) :: omega
       complex(8),intent(in) :: dsigmatilde(norb,norb)     
       complex(8),intent(in) :: proj(norb,nb)
       real(8)   ,intent(in) :: epsilon(nb)
       real(8)   ,intent(in) :: mu
       complex(8),intent(out):: dmat(nb,nb)
       complex(8),parameter  :: ci=(0.d0,1.d0)
       integer(4)            :: ib1,ib2
       complex(8)            :: cvec(nb)
!      ********************************************************************
       dmat=matmul(transpose(conjg(proj)),matmul(dsigmatilde,proj))
       cvec=1.d0/(ci*omega-epsilon(ib1)+mu)
       do ib1=1,nb
         do ib2=1,nb
           dmat(ib1,ib2)=cvec(ib1)*dmat(ib1,ib2)*cvec(ib2)
         enddo
       enddo
       return
       end
!
!     ..................................................................
      subroutine lib__invertc8(n,a,ainv)
!     ******************************************************************
!     **                                                              **
!     **  inverts the real, square matrix a                           **
!     **                                                              **
!     **  dependencies:                                               **
!     **    essl: dgeicd                                              **
!     **                                                              **
!     ******************************************************************
      implicit none
      integer(4),intent(in) :: n
      complex(8),intent(in) :: a(n,n)
      complex(8),intent(out):: ainv(n,n)
      integer(4)            :: naux
      complex(8)            :: aux(100*n)
      integer(4)            :: i,j
      integer(4)            :: ipiv(n)
      integer(4)            :: info
!     ******************************************************************
      if(n.eq.1) then
        ainv(1,1)=(1.d0,0.d0)/a(1,1)
        return
      end if
      naux=100*n
      ainv(1:n,1:n)=a(1:n,1:n)
!      call zgetrf(n,n,ainv,n,ipiv,info) !lapack
!      call zgetri(n,ainv,n,ipiv,aux,naux,info) !lapack
      return
      end
!
!      .................................................................
       subroutine writeone(norbs,nomega,omega,arr)
       implicit none
       integer(4),intent(in) :: norbs
       integer(4),intent(in) :: nomega
       real(8)   ,intent(in) :: omega(nomega)
       complex(8),intent(in) :: arr(norbs,norbs,nomega)
       integer(4)            :: nfil,iomega
       nfil=10
       open(nfil,file='greens.dat',form='formatted')
       rewind nfil
!nfil=6
       do iomega=1,nomega
         write(nfil,*)omega(iomega),real(arr(:,:,iomega)) &
      &                           ,aimag(arr(:,:,iomega))
!         write(nfil,*)omega(iomega),arr(:,:,iomega)
       enddo
       return
       end
