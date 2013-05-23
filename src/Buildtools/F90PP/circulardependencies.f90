      program main
!     **
!     ** parses a dependency file and looks for duplicated and circular dependencies
!     **
!     ** usage:  circulardependencies <dependencyfile
!     **
!     ** a dependency connects two items.
!     **
      implicit none
      character(256)          :: line
      integer(4)              :: nfilin=5
      integer(4)    ,parameter:: ndepx=500 ! x(number of dependencies)
      integer(4)    ,parameter:: nitemx=2*ndepx ! x(number of items)
      character(128)          :: dep(2,ndepx)   ! dependency
      integer(4)              :: idep(2,ndepx)  ! dependency pointers to items
      character(128)          :: item(nitemx)   ! item list
      integer(4)              :: map(nitemx)    ! mapping used to shrink item list
      integer(4)              :: ndep           ! actual number of dependencies
      integer(4)              :: nitem          ! actual number of items
      integer(4)              :: nitemlast
      integer(4)              :: ndeplast
      integer(4)              :: i,j,k,ipos,iter
      logical(4)              :: tchk
!     *******************************************************************
!
!     =================================================================
!     ==  read dependencies from standard input and count them       ==
!     =================================================================
      ndep=0
      DO
        READ(nfilin,FMT='(A)',END=1000)LINE
        ipos=index(line,':')
        ndep=ndep+1
        if(ndep.gt.ndepx) stop 'error: ndep>ndepx'
        dep(1,ndep)=adjustl(trim(line(:ipos-1)))
        dep(2,ndep)=adjustl(trim(line(ipos+1:)))
      enddo
1000  continue
!
!     =================================================================
!     ==  extract and count items                                    ==
!     =================================================================
      nitem=0
      do i=1,ndep
        do j=1,2
          tchk=.false.
          do k=1,nitem
            tchk=dep(j,i).eq.item(k)
            if(tchk) then
              idep(j,i)=k
              exit
            end if
          enddo
          if(.not.tchk) then
            nitem=nitem+1
            if(nitem.gt.nitemx) stop 'error: nitem>nitemx'
            item(nitem)=dep(j,i)
            idep(j,i)=nitem
          end if
        enddo
      enddo
!
!     =================================================================
!     ==  warn regardong duplicated dependencies                     ==
!     =================================================================
      print*,'looking for duplicated dependencies:'
      do i=1,ndep
        do j=1,i-1
          if(idep(1,i).eq.idep(1,j).and.idep(2,i).eq.idep(2,j)) then
            print*,trim(dep(1,i))//' : '//trim(dep(2,i))
            exit
          end if
        enddo
      enddo
      print*,'duplicated dependencies done'
!
!     =================================================================
!     == big loop to get rid of non-circular dependencies            ==
!     == exclude all items occuring only on the right hand side      ==
!     == and remove all dependencies related to then                 ==
!     =================================================================
      print*,'looking for circular dependencies:'
      do iter=1, 100
        nitemlast=nitem
        ndeplast=ndep
        print*,'step,',iter,'   ',nitem,' items and ',ndep,' dependencies'
!  
!       == check for items that only occur on the right side.
        do i=1,nitem
          tchk=.false.
          do j=1,ndep
            tchk=idep(1,j).eq.i        
            if(tchk) exit
          enddo
          if(.not.tchk) then
!           == remove all dependencies containing this as right side
            do j=1,ndep
              if(idep(2,j).eq.i) then
                idep(1,j)=0
                idep(2,j)=0
              end if
            enddo
            item(i)=''
          end if
        enddo
!
!       == clean up dependencies. ===========================================
        j=0
        do i=1,ndep
          j=j+1
          idep(:,j)=idep(:,i)
          dep(:,j)=dep(:,i)
          if(idep(2,j).eq.0)j=j-1
        enddo
        ndep=j          
!
!       == clean up items ================================================
        j=0
        do i=1,nitem
          j=j+1
          item(j)=item(i)
          map(i)=j
          if(item(j).eq.'')j=j-1
        enddo
        nitem=j
!
!       == update dependency pointer to items  ===========================
        do i=1,ndep
          idep(1,i)=map(idep(1,i))
          idep(2,i)=map(idep(2,i))
        enddo
!  
!       == exit from loop if the remaining set is irreducible ============
        if(nitem.eq.nitemlast.and.ndep.eq.ndeplast) exit
      enddo
!
!     ====================================================================
!     == print result on circular dependencies                          ==
!     ====================================================================
      print*,nitem,' items left'
      do i=1,nitem
        print*,'  ',item(i)
      enddo
      print*,ndep,' dependencies left'
      do i=1,ndep
        print*,trim(dep(1,i))//' : '//trim(dep(2,i))
      enddo
      print*,'if no items and dependencies are left, there are no circular dependencies'
      stop
      end
