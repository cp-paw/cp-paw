       program main
!      ** program reads an input file and produces a fortran subroutine      **
!      ** which prints the input file to a specified fortran channel         **
       character(256) :: line
!      **************************************************************************
       line='subroutine version$writeparmfile(nfil)'
       write(*,*)trim(line)
       line='integer(4),intent(in) :: nfil'
       write(*,*)trim(line)
       do
         read(*,fmt='(A)',err=100,end=100)line
         if(len_trim(line).eq.0) cycle
         do i=1,256
           if(iachar(line(i:i)).eq.39)line(i:i)='"'
         enddo
         line="write(nfil,*)'"//trim(line)//"'"
         write(*,*)trim(line)
       enddo  
100    continue
       line='return'
       write(*,*)trim(line)
       line='end'
       write(*,*)trim(line)
       stop
       end
