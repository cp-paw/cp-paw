!
!.......................................................................
module mixpulay_module
  integer(4)           ,save :: nhistory=3
  REAL(8)  ,ALLOCATABLE,save :: RHOold(:,:)  ! CHARGE DENSITY
  REAL(8)  ,ALLOCATABLE,save :: resold(:,:) ! 1CENTER DENSITY MATRIX
  REAL(8)  ,ALLOCATABLE,save :: rhoout(:) ! 1CENTER DENSITY MATRIX
  real(8)  ,allocatable,save :: a(:,:)
  real(8)  ,allocatable,save :: ainv(:,:)
  real(8)  ,allocatable,save :: test(:,:)
  real(8)  ,allocatable,save :: b(:)
  real(8)  ,allocatable,save :: alphabar(:)
  logical(4)           ,save :: tfirst=.true.
  integer(4)           ,save :: iiter
end module mixpulay_module

!     ..................................................................
      SUBROUTINE waves_mixrho_pul(rhodim,dendim,rho,denmat,tconv,&
           convpsi) ! using R and only A
!     ******************************************************************
!     ** allpy a pulay mixing scheme,                                 **
!     ** Ref: Kresse, Furtmueller Com. Mater. Sci. 6, 15 (1996)       **
!     **     Equations 83 to 87                                       **
!     ******************************************************************
      ! todo: deallocate-mechanism
      use mixpulay_module
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)    :: rhodim
      INTEGER(4) ,INTENT(IN)    :: dendim
      REAL(8)    ,INTENT(INout) :: rho(rhodim) !rho(NRL,NDIMD)
      REAL(8)    ,INTENT(INout) :: DENMAT(dendim) !(LMNXX,LMNXX,NDIMD,NAT)
      logical(4) ,intent(out)   :: tconv
      real(8)    ,intent(out)   :: convpsi
      real(8)                   :: conv=1.D-8
      real(8)                   :: g=0.5D0 !mixing factor
      real(8)                   :: svar
      real(8)                   :: r22,r12,r11
      integer(4)                :: ih1,ih2,nback
      integer(4)                :: rhox
!     ******************************************************************
      !first step
      if(tfirst) then
         iiter=1
         if(nhistory.lt.1) then
            call error$msg('nhistory must be at least 1 (for linear mixing)')
            call error$stop('waves_mixrho')
         end if
         ! allocate arrays
         allocate(rhoold(rhodim+dendim,nhistory))
         rhoold=0.D0
         allocate(resold(rhodim+dendim,nhistory))
         resold=0.D0
         print*,'mixer, allocating ',nhistory*(rhodim+dendim)*2.D0* &
              8.D0/1024.D0,' KB'
         allocate(a(nhistory,nhistory))
         a=0.D0
         allocate(ainv(nhistory,nhistory))
         ainv=0.D0
         allocate(alphabar(nhistory))
         alphabar=0.D0
         ! assign first density
         rhoold(1:rhodim,1)=rho(:)
         rhoold(rhodim+1:rhodim+dendim,1)=denmat(:)
         tfirst=.false.
         tconv=.false.
         convpsi=1.D-3 ! 3 auch gut
         return
      end if
      ! after first step: do mixing

!!$! do nothing but returning the start density
!!$rho(:)=rhoold(1:rhodim,1)
!!$denmat(:)=rhoold(rhodim+1:rhodim+dendim,1)
!!$convpsi=1.D-6 
!!$return
      allocate(rhoout(rhodim+dendim))
      iiter=iiter+1

      print*,' -------------- mixer iter: ',iiter,' ----------'
      call LIB$SCALARPRODUCTR8(.true.,rhodim,1,&
           rho(1:rhodim)-rhoold(1:rhodim,1),1,&
           rho(1:rhodim)-rhoold(1:rhodim,1),svar)
      !svar=maxval(abs(rho(1:rhodim)-rhoold(1:rhodim,1)))
      CALL MPE$COMBINE('+',svar)
write(*,"('CG-Mixer: rho   :',3f10.7)") rho(1:3)
write(*,"('CG-Mixer: rhoold:',3f10.7)") rhoold(1:3,1)
print*,'CG-Mixer iter ',iiter,' rho-diff (res) =',svar
      if(svar.gt.1.D-1) then
         convpsi=1.D-4 ! 4 auch gut
      else
         convpsi=1.D-6
      end if
      if(svar.lt.1.D-4) convpsi=1.D-8
      tconv=.false. !(svar.lt.conv) ! do not use desity criterion

      ! resold(:,1) contains the newest 
      resold(1:rhodim,1)=rho(:)-rhoold(1:rhodim,1)
      resold(rhodim+1:rhodim+dendim,1)=denmat(:) &
           -rhoold(rhodim+1:rhodim+dendim,1)
      !rhoout=rhoold(:,1)+g*resold(:,1)  ! linear part of eq. 92
      nback=min(nhistory,iiter-1)
print*,'mixer nback', nback

do ih1=1,nhistory 
   write(*,"('    mix resold(',i1,')    ',4f10.5)") ih1,resold(1:4,ih1)
end do
print*,'  mix      ======='
do ih1=1,nhistory 
   write(*,"('    mix rhoold(',i1,')    ',4f10.5)") ih1,rhoold(1:4,ih1)
end do
print*,'  mix      ======='

      rhox=rhodim ! use only density for adjusting fitting
      !rhox=rhodim+dendim ! also use denamt for fitting
!!$      do ih1=2,nback ! eq. 91
!!$         do ih2=ih1,nback
!!$            call LIB$SCALARPRODUCTR8(.false.,rhox,1,&
!!$                 resold(1:rhox,ih1),1,&
!!$                 resold(1:rhox,ih2),&
!!$                 a(ih1,ih2))
!!$            a(ih2,ih1)=a(ih1,ih2)
!!$         end do
!!$      end do
      ! eq. 87
      call LIB$SCALARPRODUCTR8(.true.,rhox,nback,&
           resold(1:rhox,1:nback),nback,resold(1:rhox,1:nback),&
           a(1:nback,1:nback))

do ih1=1,nhistory 
  write(*,"('mix a        ',5f10.5)") a(:,ih1)
end do
      if(nback.ge.1) then
         call LIB$INVERTR8(nback,A(1:nback,1:nback),AINV(1:nback,1:nback))
      end if
do ih1=1,nhistory 
  write(*,"('mix ainv     ',5f15.5)") ainv(ih1,:)
end do
!!$call lib$matmulr8(nback-1,nback-1,nback-1,a,ainv,test)
!!$do ih1=2,nhistory 
!!$  write(*,"('mix ainv*a   ',4f15.5)") test(ih1,:)
!!$end do
      svar=0.D0
      do ih1=1,nback ! eq 87
         do ih2=1,nback
            svar=svar+ainv(ih1,ih2)
         end do
      end do
print*,'mix sum ainv:',svar
      alphabar=0.D0
      do ih1=1,nback ! eq 87
         do ih2=1,nback
            alphabar(ih1)=alphabar(ih1)+ainv(ih1,ih2)
         end do
      end do
      alphabar=alphabar/svar
write(*,"('mix alphabar ',5f10.5)") alphabar(:)
      rhoout=0.D0
      do ih1=1,nback ! eq 83
         rhoout=rhoout+alphabar(ih1)*(rhoold(:,ih1)+g*resold(:,ih1))
      end do
write(*,"('CG-Mixer: rhoout:',3f10.7)") rhoout(1:3)
      ! map rhoold on rho and denmat
      rho(:)=rhoout(1:rhodim)
      denmat(:)=rhoout(rhodim+1:rhodim+dendim)
      ! propagate rhoold and resold
      do ih1=nhistory,2,-1
         rhoold(:,ih1)=rhoold(:,ih1-1)
         resold(:,ih1)=resold(:,ih1-1)
      end do
      rhoold(:,1)=rhoout
      deallocate(rhoout)
      return
      end 
!
!     ..................................................................
      SUBROUTINE waves_mixrho_pul_dR(rhodim,dendim,rho,denmat,tconv,&
           convpsi) ! using dR first large, then small
!     ******************************************************************
!     ** allpy a pulay mixing scheme,                                 **
!     ** Ref: Kresse, Furtmueller Com. Mater. Sci. 6, 15 (1996)       **
!     ** using the delta-R scheme of equations 88 to 92               **
!     ******************************************************************
      ! todo: deallocate-mechanism
      use mixpulay_module
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)    :: rhodim
      INTEGER(4) ,INTENT(IN)    :: dendim
      REAL(8)    ,INTENT(INout) :: rho(rhodim) !rho(NRL,NDIMD)
      REAL(8)    ,INTENT(INout) :: DENMAT(dendim) !(LMNXX,LMNXX,NDIMD,NAT)
      logical(4) ,intent(out)   :: tconv
      real(8)    ,intent(out)   :: convpsi
      real(8)                   :: conv=1.D-6
      real(8)                   :: g=0.5D0 !mixing factor
      real(8)                   :: svar
      real(8)                   :: r22,r12,r11
      integer(4)                :: ih1,ih2,nback
      integer(4)                :: rhox
!     ******************************************************************
      !first step
      if(tfirst) then
         iiter=1
         if(nhistory.lt.1) then
            call error$msg('nhistory must be at least 1 (for linear mixing)')
            call error$stop('waves_mixrho')
         end if
         allocate(rhoold(rhodim+dendim,nhistory))
         allocate(resold(rhodim+dendim,nhistory))
         print*,'mixer, allocating ',nhistory*(rhodim+dendim)*2.D0* &
              8.D0/1024.D0,' KB'
         allocate(a(2:nhistory,2:nhistory))
         a=0.D0
         allocate(ainv(2:nhistory,2:nhistory))
         ainv=0.D0
         allocate(test(2:nhistory,2:nhistory))
         test=0.D0
         allocate(b(2:nhistory))
         b=0.D0
         allocate(alphabar(2:nhistory))
         alphabar=0.D0
         rhoold=0.D0
         resold=0.D0
         rhoold(1:rhodim,1)=rho(:)
         rhoold(rhodim+1:rhodim+dendim,1)=denmat(:)
         tfirst=.false.
         tconv=.false.
         convpsi=1.D-3 ! 3 auch gut
         return
      end if
      ! steps after first step: do mixing
      allocate(rhoout(rhodim+dendim))
      iiter=iiter+1

print*,' -------------- mixer iter: ',iiter,' ----------'
svar=maxval(abs(rho(1:rhodim)-rhoold(1:rhodim,1)))
write(*,"('CG-Mixer: rho   :',3f10.7)") rho(1:3)
write(*,"('CG-Mixer: rhoold:',3f10.7)") rhoold(1:3,1)
print*,'CG-Mixer iter ',iiter,' rho-diff (res) =',svar
      if(svar.gt.1.D-1) then
         convpsi=1.D-4 ! 4 auch gut
      else
         convpsi=1.D-6
      end if
      if(svar.lt.1.D-4) convpsi=1.D-8
      tconv=(svar.lt.conv)

      ! resold(:,1) contains the newest 
      resold(1:rhodim,1)=rho(:)-rhoold(1:rhodim,1)
      resold(rhodim+1:rhodim+dendim,1)=denmat(:) &
           -rhoold(rhodim+1:rhodim+dendim,1)
      if(nhistory.ge.2) then
         resold(:,2)=resold(:,1)-resold(:,2) ! store delta R in resold
      end if
      rhoout=rhoold(:,1)+g*resold(:,1)  ! linear part of eq. 92
      nback=min(nhistory,iiter-1)
print*,'mixer nback', nback

do ih1=1,nhistory 
   write(*,"('    mix resold(',i1,')    ',4f10.5)") ih1,resold(1:4,ih1)
end do
print*,'        ======='
do ih1=1,nhistory 
   write(*,"('    mix rhoold(',i1,')    ',4f10.5)") ih1,rhoold(1:4,ih1)
end do
print*,'        ======='

      rhox=rhodim ! use only density for adjusting fitting
      !rhox=rhodim+dendim ! also use denamt for fitting
      ! the following loop may be done by only one call

!!$      do ih1=2,nback ! eq. 91
!!$         do ih2=ih1,nback
!!$            call LIB$SCALARPRODUCTR8(.false.,rhox,1,&
!!$                 resold(1:rhox,ih1),1,&
!!$                 resold(1:rhox,ih2),&
!!$                 a(ih1,ih2))
!!$            a(ih2,ih1)=a(ih1,ih2)
!!$         end do
!!$      end do
      call LIB$SCALARPRODUCTR8(.true.,rhox,nback-1,&
           resold(1:rhox,2:nback),nback-1,resold(1:rhox,2:nback),&
           a(2:nback,2:nback))

do ih1=2,nhistory ! eq. 91
  write(*,"('mix a        ',4f10.5)") a(:,ih1)
end do
      do ih1=2,nback ! eq 90 b=<dR_j|R1>
         call LIB$SCALARPRODUCTR8(.false.,rhox,1, &
              resold(1:rhox,ih1),1,resold(1:rhox,1),b(ih1))
      end do
write(*,"('mix b        ',4f10.5)") b
      if(nback.ge.2) then
         call LIB$INVERTR8(nback-1,A(2:nback,2:nback),AINV(2:nback,2:nback))
      end if
do ih1=2,nhistory 
  write(*,"('mix ainv     ',4f15.5)") ainv(ih1,:)
end do
!!$call lib$matmulr8(nback-1,nback-1,nback-1,a,ainv,test)
!!$do ih1=2,nhistory 
!!$  write(*,"('mix ainv*a   ',4f15.5)") test(ih1,:)
!!$end do
      alphabar=0.D0
      do ih1=2,nback ! eq 90 
         do ih2=ih1,nback
            alphabar(ih1)=alphabar(ih1)-ainv(ih1,ih2)*b(ih2)
         end do
      end do
write(*,"('mix alphabar ',4f10.5)") alphabar(:)
      do ih1=2,nback ! sum in eq 92  
         rhoout=rhoout+alphabar(ih1)*(rhoold(:,ih1-1)-rhoold(:,ih1)+ &
              g*(resold(:,ih1)))
      end do
      ! map rhoold on rho and denmat
!write(*,"('CG-Mixer: rhoout:',3f10.7)") rhoout(1:3)
      rho(:)=rhoout(1:rhodim)
      denmat(:)=rhoout(rhodim+1:rhodim+dendim)
      ! propagate rhoold and resold
      do ih1=nhistory,2,-1
         rhoold(:,ih1)=rhoold(:,ih1-1)
         resold(:,ih1)=resold(:,ih1-1)
      end do
      rhoold(:,1)=rhoout
      deallocate(rhoout)
      return
      end 


