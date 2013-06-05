       program main
!      *************************************************************************
!      ** program reads an input file and produces a fortran subroutine       **
!      ** which prints the input file to a specified fortran channel          **
!      **                                                                     **
!      ** it is used to convert the parmfile for the PAW installation into    **
!      ** a fortran routine, that writes the parmfile to a file               **
!      **                                                                     **
!      ********************************peter Bloechl, Goslar Nov. 9, 2010*******
       implicit none
       integer,parameter  :: maxlen=77
       integer,parameter  :: linelen=512
       integer,parameter  :: overheadlen=25
       character(linelen):: line
       character(maxlen) :: lineout
       integer           :: len,i
       integer,parameter :: nfilo=6    ! standard in
       integer,parameter :: nfilin=5   ! standard out
!      *************************************************************************
!
!      =========================================================================
!      ==  produce subroutine version$writeparmfile                           ==
!      =========================================================================
       write(nfilo,fmt='(A)')'subroutine version$writeparmfile()'
       do
         read(nfilin,fmt='(A)',err=100,end=100)line
         if(len_trim(line).eq.0) cycle
         do i=1,256
           if(iachar(line(i:i)).eq.39)line(i:i)='"'
         enddo
         len=len_trim(line)
         do while (len.gt.0)
           lineout="write(*,fmt='(A)')'"//trim(line(1:maxlen-overheadlen))//"'"
           write(nfilo,fmt='(A)')trim(lineout)
           if(len.gt.maxlen-overheadlen+1) then
              line=line(maxlen-overheadlen+1:len)
           else
              line=''
           end if
           len=len_trim(line)
         enddo
       enddo  
100    continue
       write(nfilo,fmt='(A)')'return'
       write(nfilo,fmt='(A)')'end'
!
!      =========================================================================
!      ==  PRODUCE SUBROUTINE VERSION$GETPAWDIR                               ==
!      =========================================================================
       WRITE(NFILO,FMT='(A)')'SUBROUTINE VERSION$GETPAWDIR(PAWDIR)'
       WRITE(NFILO,FMT='(A)')'CHARACTER(*),INTENT(OUT):: PAWDIR'
       WRITE(NFILO,FMT='(A)')"PAWDIR='XXPAWDIRXX'" ! TO BE MODIFIED
       WRITE(NFILO,FMT='(A)')'RETURN'
       WRITE(NFILO,FMT='(A)')'END'
!
!      =========================================================================
!      ==  PRODUCE SUBROUTINE VERSION$GETPAWVERSION                           ==
!      =========================================================================
       WRITE(NFILO,FMT='(A)')'SUBROUTINE VERSION$GETPAWVERSIONR(PAWVERSION)'
       WRITE(NFILO,FMT='(A)')'CHARACTER(*),INTENT(OUT):: PAWVERSION'
       WRITE(NFILO,FMT='(A)')"PAWVERSION='XXPAWVERSIONXX'" ! TO BE MODIFIED
       WRITE(NFILO,FMT='(A)')'RETURN'
       WRITE(NFILO,FMT='(A)')'END'
!
       stop
       end
