       program main
!      ** program reads an input file and produces a fortran subroutine      **
!      ** which prints the input file to a specified fortran channel         **
       integer,parameter  :: maxlen=77
       integer,parameter  :: linelen=512
       integer,parameter  :: overheadlen=25
       character(linelen):: line
       character(maxlen) :: lineout
       integer           :: len
!      **************************************************************************
       write(*,fmt='(A)')'subroutine version$writeparmfile()'
       do
         read(*,fmt='(A)',err=100,end=100)line
         if(len_trim(line).eq.0) cycle
         do i=1,256
           if(iachar(line(i:i)).eq.39)line(i:i)='"'
         enddo
         len=len_trim(line)
         do while (len.gt.0)
           lineout="write(*,fmt='(A)')'"//trim(line(1:maxlen-overheadlen))//"'"
           write(*,fmt='(A)')trim(lineout)
           line=line(maxlen-overheadlen+1:len)
           len=len_trim(line)
         enddo
       enddo  
100    continue
       write(*,fmt='(A)')'return'
       write(*,fmt='(A)')'end'
       stop
       end
