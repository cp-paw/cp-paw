module debug_module
character(32)        :: filename="debugout"
integer(8)           :: dunit=4243         
logical(4)           :: alltasks=.true.
end module debug_module

subroutine debug$nan_R8(len,output,comment)
  USE MPE_MODULE
  USE debug_module
  implicit none
  integer(4),   intent(in)           ::  len
  real(8),      intent(in)           ::  output(len)
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
           WRITE (dunit,FMT="(A,I10,F25.14)")"the following element is dubious: ",i,output(i) 
        end if
     end do
     close(dunit)
  end if
else
  OPEN (dunit, FILE=trim(adjustl(filename)), STATUS="UNKNOWN", &
       &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")

      do i=1,len
        if(output(i).ne.output(i)) then
           WRITE (dunit,FMT="(2A,I10,A,I4,F25.14)")comment," dubious element: ",i," task: ",thistask,output(i) 
        end if
     end do
     close(dunit)
  end if
end subroutine debug$nan_R8

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
