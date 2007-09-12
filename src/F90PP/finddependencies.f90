      program main
!     **
!     ** parses a fortran file from standard input 
!     ** writes dependencies to standard output
!     **
!     ** usage:  finddependencies.x filename >>dependencyfile
!     **
!     **     mode can take the values 'p' and 's'. if the value is p, the dependency 
!     **     for paw_mpelib.f90 will go to paw_mpelib_p.o, the parallel version
!     **     ands similar with s for the scalar version
!     **
!     **     filename is a the name of the fortran file. it must end by .f*, where
!     **     '*' is a wildcard for any string. preceding directories and extensions are stripped
!     **
      implicit none
      character(256)       :: line
      character(512)       :: file
      integer(4)           :: nfilin=100
      integer(4),parameter :: nmux=100
      integer(4),parameter :: nmdx=100
      integer(4)           :: nmu=0
      integer(4)           :: nmd=0
      character(128)       :: moduleused(nmux)
      character(128)       :: moduledefined(nmdx)
      logical              :: tchk
      integer(4)           :: i,j,k,ipos
      integer(4)           :: ipos_use,ipos_module
      character(128)       :: dependency
      character(4)         :: dotmod,srcslash
      character(2)         :: dotf,doto
!     *******************************************************************
      DOTF='.F'
      DOTO='.O'
      DOTMOD='.MOD'
      SRCSLASH='SRC/'
      CALL DOWNCASE(DOTF)
      CALL DOWNCASE(DOTO)
      CALL DOWNCASE(DOTMOD)
      CALL DOWNCASE(SRCSLASH)
!     ======================================================================
      call getarg(1,file)
!     == consider only fortan files
      ipos=index(file,dotf)
      if(ipos.eq.0) stop
!     ==  remove backup files of emacs
      ipos=index(file,'~')
      if(ipos.ne.0) stop
      ipos=index(file,'#')
      if(ipos.ne.0) stop
!
!     ================================================================
!     == write filename relative to the PAW main directory                 
!     ================================================================
      ipos=index(file,srcslash)
      print*,'## file:',trim(file(ipos:))

!     ================================================================
!     == parse input file                                          ==
!     ================================================================
      nmu=0
      nmu=0
      open(nfilin,file=trim(file))
      DO
        READ(nfilin,FMT='(A)',END=1000)LINE
!       == REMOVE PRECEEDING BLANKS
        LINE=ADJUSTL(LINE)
!       == MAKE LINE UPPERCASE
        CALL UPCASE(LINE)
!       == identify "use" commands
        IPOS_USE=0
        IF(LINE(1:3).EQ.'USE') IPOS_USE=1
!       == identify module definitions
        IPOS_MODULE=0
        IF(LINE(1:6).EQ.'MODULE') IPOS_MODULE=1
!       == CYCLE IF NOTHING THERE
        IF(IPOS_USE.EQ.0.AND.IPOS_MODULE.EQ.0) CYCLE
!       == REMOVE EVERYTHING PAST A COMMENT
        IPOS=INDEX(LINE,'!')
        IF(IPOS.NE.0)LINE(IPOS:)=''

        if(ipos_use.ne.0) then
          line=trim(adjustl(line(4:)))
!         == strip everything past a comma
          ipos=index(line,',')
          if(ipos.ne.0)line=line(:ipos-1)
!         == strip everything past the first blank, and get length
          ipos=index(line,' ')
          line=line(:ipos)
!         == check if already there
          tchk=.false.
          do i=1,nmu
            tchk=(line.eq.moduleused(i))
            if(tchk) exit
          enddo
          if(.not.tchk) then
            nmu=nmu+1
            moduleused(nmu)=trim(line(:ipos))
          end if
        end if
!
!       ================================================================
!       == analyse module definitions                                 ==
!       ================================================================
        if(ipos_module.ne.0) then
          line=trim(adjustl(line(7:)))
!         == strip everything past the first blank, and get length
          ipos=index(line,' ')
          line=line(:ipos)
!         == get rid of module procedures
          IF(TRIM(LINE).NE.'PROCEDURE') then
            NMD=NMD+1
            moduledefined(nmd)=trim(line)
          end if
        end if
      enddo
1000  continue

!!$print*,nmu,' modules used'
!!$do i=1,nmu
!!$  print*,moduleused(i)
!!$enddo
!!$print*,nmd,' modules defined'
!!$do i=1,nmd
!!$  print*,moduledefined(i)
!!$enddo
!
!     ===================================================================
!     == exclude dependency on modules defined in the same file        ==
!     ===================================================================
      do i=1,nmd
        k=0
        do j=1,nmu
          k=k+1
          moduleused(k)=moduleused(j)
          if(moduledefined(i).eq.moduleused(k)) then
            k=k-1
          end if
        enddo
        moduleused(k+1:nmu)=''
        nmu=k
      enddo
!
!!$print*,nmu,' modules used'
!!$do i=1,nmu
!!$  print*,moduleused(i)
!!$enddo
!!$print*,nmd,' modules defined'
!!$do i=1,nmd
!!$  print*,moduledefined(i)
!!$enddo
!
!     ===================================================================
!     == convert into module file names                                ==
!     ===================================================================
      do i=1,nmd
        call upcase(moduledefined(i))
        moduledefined(i)=trim(moduledefined(i))//dotmod
      enddo
      do i=1,nmu
        call upcase(moduleused(i))
        moduleused(i)=trim(moduleused(i))//dotmod
      enddo
!
!     ===================================================================
!     == convert into dependencies                                     ==
!     ===================================================================
!     == strip directories
      ipos=index(file,'/',.true.)
      if(ipos.ne.0)file=adjustl(file(ipos+1:))
!     == strip extension
      ipos=index(file,dotf)
      if(ipos.eq.0) stop 'error 1'
      file=file(1:ipos-1)
!     == attach mode to mpelib
      line='paw_mpelib'
      call downcase(line)
      if(file.eq.trim(line)) then
        file=trim(file)//'${MODE}'
      end if
!
      file='${OBJDIR}/'//trim(file)//doto
      do i=1,nmu
        dependency=trim(file)//' : '//trim(moduleused(i))
        write(*,*) dependency
      enddo
      do i=1,nmd
        dependency=trim(moduledefined(i))//' : '//trim(file)
        write(*,*) dependency
      enddo
      stop
      end
!
!     .....................................................................
      subroutine upcase(string)
      character(*),intent(inout):: string      
      integer(4)                :: length
      integer(4)                :: i
!     *********************************************************************
      length=len(string)
      do i=1,length
        iascii=iachar(string(i:i))
        if(iascii.lt.97) cycle
        if(iascii.gt.122) cycle
        iascii=iascii-97+65
        string(i:i)=char(iascii)
      enddo
      return
      end
!
!     .....................................................................
      subroutine downcase(string)
      character(*),intent(inout):: string      
      integer(4)                :: length
      integer(4)                :: i
!     *********************************************************************
      length=len(string)
      do i=1,length
        iascii=iachar(string(i:i))
        if(iascii.lt.65) cycle
        if(iascii.gt.90) cycle
        iascii=iascii-65+97
        string(i:i)=char(iascii)
      enddo
      return
      end
