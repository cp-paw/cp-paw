     PROGRAM MAIN
       IMPLICIT NONE
       
       character(80)   :: line
       integer(4)      :: i, nframes, nat
       logical(4)      :: ok
       
       
       open(UNIT=101,FILE='case.info',STATUS="OLD")
       OPEN(UNIT=102,FILE='charges.movie.xyz',STATUS="UNKNOWN",POSITION="REWIND")
       ok=.true.
       
       nframes= 0
       DO WHILE(OK)
          read(101,FMT='(A80)',END=500) line
          if(line(11:20).eq.'    101010') THEN
             nframes= nframes + 1
             write(102,FMT='(a10,I10)') "    NONAME",nframes
             print*,"FOUND FRAME NUMBER ",nframes
          else
             write(102,FMT='(A80)') line
          end if
       end DO
500 ok=.false.
       print*,"PROGRAM FINISHED"
       close(101)
       close(102)
       STOP
     end PROGRAM MAIN
