module debug_module
character(32)        :: filename="debugout"
integer(8)           :: dunit=4243         
logical(4)           :: alltasks=.true.
end module debug_module
!      ........................................................................
       subroutine debug$nan_R8(len,output,comment)
       USE MPE_MODULE
       USE debug_module
       implicit none
       integer(4),   intent(in)           ::  len
       real(8),      intent(in)           ::  output(len)
       character(*), intent(in)           ::  comment
       integer(4)                         ::  Ntasks,Thistask
       integer(4)                         :: i
 !     **************************************************************************
       CALL MPE$QUERY(NTASKS,THISTASK)
       if(.not.alltasks) then
         OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
         &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")
         WRITE (dunit, FMT="(3A)", ADVANCE="YES" )"==============",comment,"==============" 
         if(thistask.eq.1) then
           do i=1,len
             if(output(i).ne.output(i)) then
               WRITE (dunit,FMT="(A,I10,F25.14)")"the following element is dubious: ",i,output(i) 
             end if
           end do
           close(dunit)
         end if
       else
         OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
      &     ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")
       do i=1,len
         if(output(i).ne.output(i)) then
           WRITE (dunit,FMT="(2A,I10,A,I4,F25.14)")comment," dubious element: ",i," task: ",thistask,output(i) 
         end if
       end do
       close(dunit)
       end if
       end subroutine debug$nan_R8
!
!.....................................................................................
subroutine debug$nan_C8(len,output,comment)
  USE MPE_MODULE
  USE debug_module
  implicit none
  integer(4),   intent(in)           ::  len
  complex(8),   intent(in)           ::  output(len)
  character(*), intent(in)           ::  comment
  integer(4)                         ::  Ntasks,Thistask
  integer(4)                         :: i
 !...................................................
  CALL MPE$QUERY(NTASKS,THISTASK)
if(.not.alltasks) then
  OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
       &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")
  WRITE (dunit, FMT="(3A)", ADVANCE="YES" )"==============",comment,"==============" 

  if(thistask.eq.1) then
     do i=1,len
        if(output(i).ne.output(i)) then
           WRITE (dunit,FMT="(A,I10,2F25.14)")"the following element is dubious: ",i,output(i) 
        end if
     end do
     close(dunit)
  end if
else
  OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
       &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")
     do i=1,len
        if(output(i).ne.output(i)) then
           WRITE (dunit,FMT="(2A,I10,A,I4,2F25.14)")comment," dubious element: ",i," task: ",thistask,output(i) 
        end if
     end do
     close(dunit)

end if
end subroutine debug$nan_C8


subroutine debug$nan_I4(len,output,comment)
  USE MPE_MODULE
  USE debug_module
  implicit none
  integer(4),   intent(in)           ::  len
  integer(4),   intent(in)           ::  output(len)
  character(*), intent(in)           ::  comment
  integer(4)                         ::  Ntasks,Thistask
  integer(4)                         :: i
 !...................................................

  CALL MPE$QUERY(NTASKS,THISTASK)
if(.not.alltasks) then
  OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
       &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")
  WRITE (dunit, FMT="(3A)", ADVANCE="YES" )"==============",comment,"==============" 

  if(thistask.eq.1) then
     do i=1,len
        if(output(i).ne.output(i)) then
           WRITE (dunit,FMT="(A,I10,I10)")"the following element is dubious: ",i,output(i) 
        end if
     end do
     close(dunit)
  end if
else
  OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
       &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")

     do i=1,len
        if(output(i).ne.output(i)) then
           WRITE (dunit,FMT="(2A,I10,A,I4,I10)")comment," dubious element: ",i," task: ",thistask,output(i) 
        end if
     end do
     close(dunit)

end if
end subroutine debug$nan_I4






subroutine debug$write_R8(len,output,comment)
  USE MPE_MODULE
  USE debug_module
  implicit none
  integer(4),   intent(in)           ::  len
  real(8),      intent(in)           ::  output(len)
  character(*), intent(in)           ::  comment
  integer(4)                         ::  Ntasks,Thistask
  !...................................................

  CALL MPE$QUERY(NTASKS,THISTASK)
  if(thistask.eq.1) then
     OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")
     WRITE (dunit, FMT="(3A)", ADVANCE="YES" )"==============",comment,"==============" 
!     WRITE (dunit,*)output 
     WRITE (dunit,FMT="(F25.14)")output 
     close(dunit)
  end if
  
end subroutine debug$write_R8





subroutine debug$write_I4(len,output,comment)
  USE MPE_MODULE
  USE debug_module
  implicit none
  integer(4),   intent(in)           ::  len
  integer(4),   intent(in)           ::  output(len)
  character(*), intent(in)           ::  comment
  integer(4)                         ::  Ntasks,Thistask
  !...................................................

  CALL MPE$QUERY(NTASKS,THISTASK)
  if(thistask.eq.1) then
     OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")
     WRITE (UNIT=dunit, FMT="(3A)", ADVANCE="YES" )"==============",comment,"==============" 
     WRITE (dunit,*)output 
     close(dunit)
  end if
  
end subroutine debug$write_I4
!
!  
subroutine debug$write_l4(len,output,comment)
  USE MPE_MODULE
  USE debug_module
  implicit none
  integer(4),   intent(in)           ::  len
  logical(4),   intent(in)           ::  output(len)
  character(*), intent(in)           ::  comment
  integer(4)                         ::  Ntasks,Thistask
  !...................................................

  CALL MPE$QUERY(NTASKS,THISTASK)
  if(thistask.eq.1) then
     OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")
     WRITE (UNIT=dunit, FMT="(3A)", ADVANCE="YES" )"==============",comment,"==============" 
     WRITE (dunit,*)output 
     close(dunit)
  end if
  
end subroutine debug$write_l4


subroutine debug$write_C8(len,output,comment)
  USE MPE_MODULE
  USE debug_module
  implicit none
  integer(4),   intent(in)           ::  len
  complex(8),   intent(in)           ::  output(len)
  character(*), intent(in)           ::  comment
  integer(4)                         ::  Ntasks,Thistask
  !...................................................

  CALL MPE$QUERY(NTASKS,THISTASK)
  if(thistask.eq.1) then
     OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")
     WRITE (UNIT=dunit, FMT="(3A)", ADVANCE="YES" )"==============",comment,"==============" 
!     WRITE (dunit,*)output 
     WRITE (dunit,FMT="(2f25.14)")output 
     close(dunit)
  end if
  
end subroutine debug$write_C8







!the following code is here for historical reasons and will be deleted!


subroutine paw$debug_R8(len,output,comment)
  USE MPE_MODULE
  USE debug_module
  implicit none
  integer(4),   intent(in)           ::  len
  real(8),      intent(in)           ::  output(len)
  character(*), intent(in)           ::  comment
  integer(4)                         ::  Ntasks,Thistask
  !...................................................

  CALL MPE$QUERY(NTASKS,THISTASK)
  if(thistask.eq.1) then
     OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")
     WRITE (dunit, FMT="(3A)", ADVANCE="YES" )"==============",comment,"==============" 
!     WRITE (dunit,*)output 
     WRITE (dunit,FMT="(F25.14)")output 
     close(dunit)
  end if
  
end subroutine paw$debug_R8





subroutine paw$debug_I4(len,output,comment)
  USE MPE_MODULE
  USE debug_module
  implicit none
  integer(4),   intent(in)           ::  len
  integer(4),   intent(in)           ::  output(len)
  character(*), intent(in)           ::  comment
  integer(4)                         ::  Ntasks,Thistask
  !...................................................

  CALL MPE$QUERY(NTASKS,THISTASK)
  if(thistask.eq.1) then
     OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")
     WRITE (UNIT=dunit, FMT="(3A)", ADVANCE="YES" )"==============",comment,"==============" 
     WRITE (dunit,*)output 
     close(dunit)
  end if
  
end subroutine paw$debug_I4

subroutine paw$debug_C8(len,output,comment)
  USE MPE_MODULE
  USE debug_module
  implicit none
  integer(4),   intent(in)           ::  len
  complex(8),   intent(in)           ::  output(len)
  character(*), intent(in)           ::  comment
  integer(4)                         ::  Ntasks,Thistask
  !...................................................

  CALL MPE$QUERY(NTASKS,THISTASK)
  if(thistask.eq.1) then
     OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")
     WRITE (UNIT=dunit, FMT="(3A)", ADVANCE="YES" )"==============",comment,"==============" 
!     WRITE (dunit,*)output 
     WRITE (dunit,FMT="(2f25.14)")output 
     close(dunit)
  end if
end subroutine paw$debug_C8
!
!     ................................................................................
      subroutine debug$round_r8(len,val,tol)
      USE MPE_MODULE
      USE debug_module
      implicit none
      integer(4),intent(in)    ::  len
      real(8)   ,intent(inout) ::  val(len)
      real(8)   ,intent(in)    ::  tol
!     **********************************************************************************
      val(:)=nint(val(:)/tol)*tol
      end subroutine debug$round_r8
!
!     ................................................................................
      subroutine debug$wavefunctions()
      use waves_module
      USE debug_module
      IMPLICIT NONE
      COMPLEX(8)  ,ALLOCATABLE:: PSIG(:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSI1(:,:)
      INTEGER(4)              :: NTASKS,THISTASK
      INTEGER(4)              :: IKPT,ISPIN,IB,IDIM,IWAVE
      INTEGER(4)              :: NG1,NGG,NGL,NBH,NB
      REAL(8)                 :: XK(3)
      REAL(8)     ,ALLOCATABLE:: GVECL(:,:)
      REAL(8)     ,ALLOCATABLE:: GVECG(:,:)
      REAL(8)                 :: GBAS(3,3)
      INTEGER(4)              :: NREC1,ISVAR
      CHARACTER(8)            :: KEY
      REAL(8)                 :: RBAS(3,3)
      REAL(8)                 :: a,b
      INTEGER(4)  ,ALLOCATABLE:: IGVECG(:,:)
      INTEGER(4)  ,ALLOCATABLE:: IGVECL(:,:)
      INTEGER(4)              :: ig
      character(64)           :: comment='wave function printout'
      real(8)      ,parameter :: tol=1.d-12
!     **********************************************************************************
      CALL MPE$QUERY(NTASKS,THISTASK)
print*,'Anfang debug$wavefunctions',thistask
  if(thistask.eq.1) then
     OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")
     WRITE (UNIT=dunit, FMT="(3A)", ADVANCE="YES" )"==============",comment,"==============" 
  end if
      DO IKPT=1,NKPT
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NGL=GSET%NGL
        NBH=THIS%NBH
        NB=THIS%NB
        CALL PLANEWAVE$GETI4('NGG',NGG)
        ALLOCATE(IGVECL(3,NGL))
        CALL PLANEWAVE$GETI4A('IGVEC',3*NGL,IGVECL)
        ALLOCATE(IGVECG(3,NGG))
        CALL PLANEWAVE$COLLECTI4(3,NGL,IGVECL,NGG,IGVECG)
        DEALLOCATE(IGVECL)
        ALLOCATE(PSIG(NGG,NDIM))
        ALLOCATE(PSI1(NGL,NDIM))
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          DO IB=1,NBH
            DO IDIM=1,NDIM
              PSI1(:,IDIM)=THIS%Hpsi(:,IDIM,IB)  !<<<<=================================
              do ig=1,ngl
                a=real(psi1(ig,idim))
                b=aimag(psi1(ig,idim))
                a=tol*real(nint(a/tol),kind=8)
                b=tol*real(nint(b/tol),kind=8)
                psi1(ig,idim)=cmplx(a,b,kind=8)
              enddo
            ENDDO
            DO IDIM=1,NDIM
              CALL PLANEWAVE$COLLECTC8(1,NGL,PSI1(1,IDIM),NGG,PSIG(1,IDIM))
            ENDDO
            if(thistask.eq.1) then
              do ig=1,ngg
                write(dunit,*)ib,igvecg(:,ig),psig(ig,:)
              ENDDO
            end if
          ENDDO
          DEALLOCATE(PSI1)
          DEALLOCATE(PSIG)
        ENDDO 
        DEALLOCATE(IGVECG)
      enddo
      if(thistask.eq.1) then
        close(dunit)
      end if
print*,'Ende debug$wavefunctions',thistask
      call mpe$sync
      call error$normalstop
      end subroutine debug$wavefunctions

