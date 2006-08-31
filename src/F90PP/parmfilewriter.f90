       program main
!      ** program reads an input file and produces a fortran subroutine      **
!      ** which prints the input file to a specified fortran channel         **
       character(256) :: line
!      **************************************************************************
       write(*,fmt='(A)')'subroutine version$writeparmfile(nfil)'
       write(*,fmt='(A)')'integer(4),intent(in) :: nfil'
       do
         read(*,fmt='(A)',err=100,end=100)line
         if(len_trim(line).eq.0) cycle
         do i=1,256
           if(iachar(line(i:i)).eq.39)line(i:i)='"'
         enddo
         if(len_trim(line).gt.256-24) then
           stop 'error in writeparmfile: string too large'
         end if
         line="write(*,fmt='(A)')'"//trim(line)//"'"
         write(*,fmt='(A)')trim(line)
       enddo  
100    continue
       write(*,fmt='(A)')'return'
       write(*,fmt='(A)')'end'
       stop
       end
