!
!CAUTION! -- CAUTION! -- CAUTION! -- CAUTION! -- CAUTION! -- CAUTION! -- CAUTION
!
! CAUTION! THIS FILES CONTAINS DOLLARS THAT MUST NOT BE CONVERTED TO "__" 
! AS DONE BY THE PREPROCESSOR. FURTHERMORE IT CONTAINS STRING WITH
! UPPERCASE AND LOWER CASE LETTERS THAT MUST NOT BE CONVERTED.
!
!CAUTION! -- CAUTION! -- CAUTION! -- CAUTION! -- CAUTION! -- CAUTION! -- CAUTION
!
module stpreport_module
character(256),save :: precommand
character(256),save :: tmpdir
end module stpreport_module
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       Program main
       use stpreport_module
       implicit none
       character(512) :: filename
       character(512) :: rootname
       integer        :: nfil=1321
       character(512) :: command
       character(256)  :: texfile
       integer        :: nargs,i
       integer        :: rc
       real(8)        :: svar
!      *************************************************************************
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT
!
!      =========================================================================
!      == obtain filename and rootname                                        ==
!      =========================================================================
       NARGS=COMMAND_ARGUMENT_COUNT()
       if(nargs.ne.1) then
         print*,'filename is not given'
         stop 'command line error'
       end if
       call get_command_argument(1,filename)
       i=index(filename,'.')
       rootname=filename(1:i-1)
!
!      =========================================================================
!      == define a random directory name as temporary workspace ================
!      =========================================================================
       call init_random_seed() ! old: call random_seed()
       call random_number(svar)
       tmpdir=''
       write(tmpdir,fmt='(i12)')abs(nint(svar*1.d+6))
       tmpdir='tmpdir'//trim(adjustl(tmpdir))//'/'
       command='mkdir '//trim(adjustl(tmpdir))
!       call lib__system(command) !       rc=system(command)
       call execute_command_line(trim(command),exitstat=rc)
       if(rc.ne.0) then
         call error__msg('command line argument failed')
         call error__chval('command',trim(command))
         call error__i4val('return code',rc)
         call error__stop('main')
       end if
       command='cp '//trim(filename)//' '//trim(tmpdir)
!       call lib__system(command) ! rc=system(command)
       call execute_command_line(trim(command),exitstat=rc)
       if(rc.ne.0) then
         call error__msg('command line argument failed')
         call error__chval('command',trim(command))
         call error__i4val('return code',rc)
         call error__stop('main')
       end if
       precommand='cd '//trim(tmpdir)//'; '
!
!      =========================================================================
!      == read output file of setup construction and create intermediate files==
!      =========================================================================
       call executestpa(filename)
!
!      =========================================================================
!      == prepare figures                                                     ==
!      =========================================================================
       call gnuplotpotentials()
       call gnuplotscattering()
       call gnuplotorbitals()
       call gnuplotprojectors()
       call gnuplotwf()
!
!      =========================================================================
!      == write tex file                                                      ==
!      =========================================================================
       texfile='tmp.tex'
       open(unit=nfil,file=trim(tmpdir)//texfile,form='formatted')
       call writetexfile(nfil)
       close(nfil)
!
!      =========================================================================
!      == execute gnuplot to produce postscript file
!      =========================================================================
       command='pdflatex '//adjustl(texfile)
       command=trim(precommand)//command
!      ==== call lib__system(command)
!       call lib__system(trim(command)) ! rc=system(trim(command))
       call execute_command_line(trim(command),exitstat=rc)
       if(rc.ne.0) then
         call error__msg('command line argument failed')
         call error__chval('command',trim(command))
         call error__i4val('return code',rc)
         call error__stop('main')
       end if
       command='cp '//trim(tmpdir)//'tmp.pdf '//trim(rootname)//'.pdf'
!       call lib__system(trim(command)) ! rc=system(trim(command))
       call execute_command_line(trim(command),exitstat=rc)
       if(rc.ne.0) then
         call error__msg('command line argument failed')
         call error__chval('command',trim(command))
         call error__i4val('return code',rc)
         call error__stop('main')
       end if
!
       CALL ERROR$NORMALSTOP()
       stop
       end
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE init_random_seed()
!      *************************************************************************
!      ** initialize a random number sequence using the system clock          **
!      ** from http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html    **
!      *************************************************************************
       implicit none
       INTEGER :: i, n, clock
       INTEGER, DIMENSION(:), ALLOCATABLE :: seed
!      *************************************************************************
       CALL RANDOM_SEED(size = n)
       ALLOCATE(seed(n))
       CALL SYSTEM_CLOCK(COUNT=clock)
       seed = clock + 37 * (/ (i - 1, i = 1, n) /)
       CALL RANDOM_SEED(PUT = seed)
       DEALLOCATE(seed)
       return
       END 
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine executestpa(filename)
       use stpreport_module
       implicit none
       character(*),intent(in) :: filename
       integer      ,parameter :: stpa_nargs=(29) 
       character(50)           :: stpa_args(stpa_nargs)
       character(2048)         :: command
       integer                 :: i
       integer                 :: rc
!      *************************************************************************
       stpa_args( 1)='-s nb               -o nb.dat'
       stpa_args( 2)='-s id               -o id.dat'
       stpa_args( 3)='-s parms.rcsm       -o rcsm.dat'                                   
       stpa_args( 4)='-s parms.rbox       -o rbox.dat'                                   
       stpa_args( 5)='-s lpro             -o lpro.dat'
       stpa_args( 6)='-s npro             -o npro.dat'
       stpa_args( 7)='-s zv               -o zv.dat'
       stpa_args( 8)='-s parms.psi.type   -o psi.type.dat'
       stpa_args( 9)='-s parms.psi.rcl    -o psi.rcl.dat'
       stpa_args(10)='-s parms.psi.lambda -o psi.lambda.dat'
       stpa_args(11)='-s parms.core.rc    -o core.rc.dat'
       stpa_args(12)='-s parms.core.pow   -o core.pow.dat'
       stpa_args(13)='-s parms.pot.rc     -o pot.rc.dat'
       stpa_args(14)='-s parms.pot.pow    -o pot.pow.dat'
       stpa_args(15)='-s atom.f           -o atom.f.dat'
       stpa_args(16)='-s atom.e           -o atom.e.dat'
       stpa_args(17)='-s atom.l           -o atom.l.dat'
       stpa_args(18)='-s AEPHI            -o AEPHI.dat'
       stpa_args(19)='-s PSPHI            -o PSPHI.dat'
       stpa_args(20)='-s NLPHI            -o NLPHI.dat'
       stpa_args(21)='-s QPHI             -o QPHI.dat'
       stpa_args(22)='-s PRO              -o PRO.dat'
       stpa_args(23)='-s AEPHIDOT         -o AEPHIDOT.dat'
       stpa_args(24)='-s PSPHIDOT         -o PSPHIDOT.dat'
       stpa_args(25)='-s NLPHIDOT         -o NLPHIDOT.dat'
       stpa_args(26)='-s AEPSI            -o AEPSI.dat'
       stpa_args(27)='-s UPSI             -o UPSI.dat'
       stpa_args(28)='-s POT              -o POT.dat'
       stpa_args(29)='-s SCATTERING       -o SCATTERING.dat'
       command='paw_stpa.x'
       do i=1,stpa_nargs
         command=trim(adjustl(command))//' '//trim(adjustl(stpa_args(i)))
       enddo
       command=trim(adjustl(command))//' '//trim(adjustl(filename))
       command=trim(precommand)//command
!       call lib__system(command) ! call system(command)
       call execute_command_line(trim(command),exitstat=rc)
       if(rc.ne.0) then
         call error__msg('command line argument failed')
         call error__chval('command',trim(command))
         call error__i4val('return code',rc)
         call error__stop('executestpa')
       end if
       return
       end
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine gnuplotpotentials()
!      *************************************************************************
!      ** converts the input file "pot.dat" into an eps file                  **
!      *************************************************************************
       use stpreport_module
       implicit none
       integer      ,parameter :: nline=20
       character(80)           :: line(nline)
       character(1024)         :: filename
       integer                 :: nfilgp=1001
       character(1024)         :: command
       integer                 :: i
       integer                 :: rc
!      *************************************************************************
!
!      =========================================================================
!      == write gnuplot script                                                ==
!      =========================================================================
       line(:)=' '
       line( 1)='set terminal postscript eps enhanced color lw 4.0 24'
       line( 1)='set terminal pdf enhanced color lw 4'
       line( 2)='set xrange [0.0:3.2]'
       line( 3)='set yrange [-100:0]'
       line( 4)='set xlabel "r"'
       line( 5)='set ylabel ""'
!       line( 6)='set data style lines'
       line( 7)='set xzeroaxis'
       line( 8)='set border 3'
       line( 9)='set xtics nomirror'
       line(10)='set ytics nomirror'
       line(11)='set output "potential.pdf"'
       line(12)='set title "Potentials"'
       line(13)='plot "POT.dat" using 1:2 ls 1 title "AEPOT" with lines \ '    
       line(14)='    ,"POT.dat" using 1:3 ls 2 title "PSPOT" with lines \ '    
       line(15)='    ,"POT.dat" using 1:4 ls 3 title "VPOT" with lines'
!
       filename='pot.gp'
       open(unit=nfilgp,file=trim(tmpdir)//filename,form='formatted')
       do i=1,15
         write(nfilgp,fmt='(A)')trim(line(i))
       enddo
       close(nfilgp)
!
!      =========================================================================
!      == execute gnuplot to produce postscript file
!      =========================================================================
       command='gnuplot '//trim(adjustl(filename))
       command=trim(precommand)//command
!      ==== call lib__system(command)
!       call lib__system(command) ! rc=system(command)
       call execute_command_line(trim(command),exitstat=rc)
       if(rc.ne.0) then
         call error__msg('command line argument failed')
         call error__chval('command',trim(command))
         call error__i4val('return code',rc)
         call error__stop('gnuplotpotentials')
       end if
       return 
       end
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine gnuplotscattering()
!      *************************************************************************
!      ** converts the input file "pot.dat" into an eps file                  **
!      *************************************************************************
       use stpreport_module
       implicit none
       integer      ,parameter :: nline=20
       character(124)          :: line(nline)
       character(1024)         :: filename
       integer                 :: nfilgp=1001
       character(1024)         :: command
       integer(4),parameter    :: nvalx=10
       integer(4)              :: nval
       integer(4)              :: ival(nvalx)
       integer(4)              :: lx
       integer(4)              :: il,i
       integer(4)              :: iline
       integer(4)              :: rc
!      *************************************************************************
!
!      =========================================================================
!      == write gnuplot script                                                ==
!      =========================================================================
       line(:)=' '
       line( 1)='set terminal postscript eps enhanced color lw 3.0 18'
       line( 1)='set terminal pdf enhanced color lw 3'
       line( 2)='set output "scatt.pdf"'
       line( 3)='set xrange [-1.1:1.1]'
       line( 4)='set yrange [0:*]'
       line( 5)='set xlabel "Energy[H]"'
       line( 6)='set ylabel ""'
!       line( 7)='set data style lines'
       line( 8)='set border 3'
       line( 9)='set xtics nomirror'
       line(10)='set ytics nomirror'
       line(11)='set title "Scattering"'
       call readlinei4('lpro.dat',nvalx,nval,ival)
       lx=maxval(ival(1:nval))
       write(line(12),fmt='(a,i1,a)')'plot "SCATTERING.dat" using 1:',2 &
      &                                ,' ls 1  title "AEPHASE" \'
       write(line(13),fmt='(a,i1,a)')',"SCATTERING.dat" using 1:',lx+3 &
      &                               ,' ls 2 title "PAW-PHASE" \'
       iline=13
       do il=2,lx+1
         iline=iline+1
         write(line(iline),fmt='(a,i1,a)')',"SCATTERING.dat" using 1:',1+il &
      &                                   ,' ls 1  title "" \'    
         iline=iline+1
         write(line(iline),fmt='(a,i1,a)')',"SCATTERING.dat" using 1:',lx+2+il &
      &                              ,' ls 2  title "" \'
       enddo
       iline=iline+1
       line(iline:)=''
!
       filename='scatt.gp'
       open(unit=nfilgp,file=trim(tmpdir)//filename,form='formatted')
       do i=1,iline
         write(nfilgp,fmt='(A)')trim(line(i))
       enddo
       close(nfilgp)
!
!      =========================================================================
!      == execute gnuplot to produce postscript file
!      =========================================================================
       command='gnuplot '//adjustl(filename)
       command=trim(precommand)//command
!      ==== call lib__system(command)
!       call lib__system(command) ! rc=system(command)
       call execute_command_line(trim(command),exitstat=rc)
       if(rc.ne.0) then
         call error__msg('command line argument failed')
         call error__chval('command',trim(command))
         call error__i4val('return code',rc)
         call error__stop('gnuplotscattering')
       end if
       return 
       end
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine gnuplotorbitals()
!      *************************************************************************
!      ** converts the input file "pot.dat" into an eps file                  **
!      *************************************************************************
       use stpreport_module
       implicit none
       integer      ,parameter :: nline=20
       character(124)          :: line(nline)
       character(1024)         :: filename='phi.gp'
       integer                 :: nfilgp=1001
       character(1024)         :: command
       integer(4),parameter    :: nvalx=10
       integer(4)              :: nval
       integer(4)              :: ival(nvalx)
       integer(4)              :: lx
       integer(4)              :: l,ipro
       integer(4)              :: iline,i
       integer(4)              :: rc
       logical                 :: tfirst
!      *************************************************************************
       call readlinei4('lpro.dat',nvalx,nval,ival)
       lx=maxval(ival(1:nval))
!
!      =========================================================================
!      == write gnuplot script                                                ==
!      =========================================================================
       do l=0,lx
         line(:)=' '
         line( 1)='set terminal pdf enhanced color lw 3'
         write(line( 2),fmt='(a,i1,a)')'set output "phi_',l,'.pdf"'
         line( 3)='set xrange [0:3.2]'
         line( 4)='set yrange [-2.5:2.5]'
         line( 5)='set xlabel "r[a0]"'
         line( 6)='set ylabel ""'
!         line( 7)='set data style lines'
         line( 8)='set border 3'
         line( 9)='set xtics nomirror'
         line(10)='set ytics nomirror'
         line(11)='set title "Phi"'
         tfirst=.true.
         iline=11
         do ipro=1,nval
           if(ival(ipro).ne.l) cycle
            iline=iline+1
            if(tfirst) then
              write(line(iline),fmt='(a,i1,a)')'plot "AEPHI.dat" using 1:',1+ipro &
      &                                ,' ls 1 with lines title "AE-phi" \ '    
              tfirst=.false.
            else      
              write(line(iline),fmt='(a,i1,a)')', "AEPHI.dat" using 1:',1+ipro &
      &                                ,' ls 2 with lines title "AE-phi" \ '    
            end if
            iline=iline+1
            write(line(iline),fmt='(a,i1,a)')',"PSPHI.dat" using 1:',1+ipro &
      &                                    ,' ls 3 with lines title "PS-phi" \ '
            iline=iline+1
            write(line(iline),fmt='(a,i1,a)')',"NLPHI.dat" using 1:',1+ipro &
      &                          ,' ls 4 with lines title "Nodeless-phi" \ '
            iline=iline+1
            write(line(iline),fmt='(a,i1,a)')',"QPHI.dat" using 1:',1+ipro &
      &                                ,' ls 5 with lines title "qphi" \ '    
         enddo
         iline=iline+1
         line(iline:)=''
         open(unit=nfilgp,file=trim(tmpdir)//filename,form='formatted')
         do i=1,iline
           write(nfilgp,fmt='(A)')trim(line(i))
         enddo
         close(nfilgp)
!
!        =======================================================================
!        == execute gnuplot to produce postscript file
!        =======================================================================
         command='gnuplot '//adjustl(filename)
         command=trim(precommand)//command
!        ==== call lib__system(command)
!         call lib__system(command) ! rc=system(command)
         call execute_command_line(trim(command),exitstat=rc)
         if(rc.ne.0) then
           call error__msg('command line argument failed')
           call error__chval('command',trim(command))
           call error__i4val('return code',rc)
           call error__stop('gnuplotorbitals')
         end if
       enddo

       return 
       end
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine gnuplotprojectors()
!      *************************************************************************
!      ** converts the input file "pot.dat" into an eps file                  **
!      *************************************************************************
       use stpreport_module
       implicit none
       integer      ,parameter :: nline=20
       character(124)          :: line(nline)
       character(1024)         :: filename='phi.gp'
       integer                 :: nfilgp=1001
       character(1024)         :: command
       integer                 :: i
       integer                 :: rc
       integer(4),parameter    :: nvalx=10
       integer(4)              :: nval
       integer(4)              :: ival(nvalx)
       integer(4)              :: lx
       integer(4)              :: l,ipro
       integer(4)              :: iline
       logical                 :: tfirst
!      *************************************************************************
       call readlinei4('lpro.dat',nvalx,nval,ival)
       lx=maxval(ival(1:nval))
!
!      =========================================================================
!      == write gnuplot script                                                ==
!      =========================================================================
       do l=0,lx
         line(:)=' '
         line( 1)='set terminal pdf enhanced color lw 3'
         write(line( 2),fmt='(a,i1,a)')'set output "pro_',l,'.pdf"'
         line( 3)='set xrange [0:3.2]'
         line( 4)='set yrange [*:*]'
         line( 5)='set xlabel "r[a0]"'
         line( 6)='set ylabel ""'
!         line( 7)='set data style lines'
         line( 8)='set border 3'
         line( 9)='set xtics nomirror'
         line(10)='set ytics nomirror'
         line(11)='set title "Projectors"'
         tfirst=.true.
         iline=11
         do ipro=1,nval
           if(ival(ipro).ne.l) cycle
            iline=iline+1
            if(tfirst) then
              write(line(iline),fmt='(a,i1,a)')'plot "PRO.dat" using 1:',1+ipro &
      &                                ,' ls 1 with lines title "PRO" \ '    
              tfirst=.false.
            else      
              write(line(iline),fmt='(a,i1,a)')', "PRO.dat" using 1:',1+ipro &
      &                                ,' ls 2 with lines title "Pro" \ '    
            end if
         enddo
         iline=iline+1
         line(iline:)=''
         open(unit=nfilgp,file=trim(tmpdir)//filename,form='formatted')
         do i=1,iline
           write(nfilgp,fmt='(A)')trim(line(i))
         enddo
         close(nfilgp)
!
!        =======================================================================
!        == execute gnuplot to produce postscript file
!        =======================================================================
         command='gnuplot '//adjustl(filename)
         command=trim(precommand)//command
!        ==== call lib__system(command)
!         call lib__system(command) ! rc=system(command)
         call execute_command_line(trim(command),exitstat=rc)
         if(rc.ne.0) then
           call error__msg('command line argument failed')
           call error__chval('command',trim(command))
           call error__i4val('return code',rc)
           call error__stop('gnuplotprojectors')
         end if
       enddo

       return 
       end
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine gnuplotwf()
!      *************************************************************************
!      ** converts the input file "pot.dat" into an eps file                  **
!      *************************************************************************
       use stpreport_module
       implicit none
       integer      ,parameter :: nline=20
       character(124)          :: line(nline)
       character(1024)         :: filename='psi.gp'
       character(64)           :: string
       integer                 :: nfilgp=1001
       integer                 :: nfil1=11
       integer                 :: nfil2=12
       character(1024)         :: command
       integer(4),parameter    :: nvalx=30
       integer(4)              :: nval
       integer(4)              :: ival(nvalx)
       integer(4)              :: i,ib
       integer(4)              :: rc
       integer(4)              :: nb
       integer(4)              :: nr
       real(8)   ,allocatable  :: r(:)
       real(8)   ,allocatable  :: psi(:,:)
       real(8)   ,allocatable  :: upsi(:,:)
       integer(4),allocatable  :: lb(:)
       real(8)                 :: svar
!      *************************************************************************
       call readlinei4('nb.dat',nvalx,nval,ival)
       nb=ival(1)
       allocate(lb(nb))
       call readlinei4('atom.l.dat',nvalx,nval,lb)
!
!      =========================================================================
!      == renormalize so that max(uphi)=1                                     ==
!      =========================================================================
       command='wc -l UPSI.dat > nr.dat'
       command=trim(precommand)//command
!       call lib__system(command) ! rc=system(command)
       call execute_command_line(trim(command),exitstat=rc)
       if(rc.ne.0) then
         call error__msg('command line argument failed')
         call error__chval('command',trim(command))
         call error__i4val('return code',rc)
         call error__stop('gnuplotwf')
       end if
       call readline('nr.dat',line(1))
       read(line(1),*)nr
       allocate(r(nr))
       allocate(psi(nr,nb))
       allocate(upsi(nr,nb))
       open(nfil1,file=trim(tmpdir)//'AEPSI.dat')
       open(nfil2,file=trim(tmpdir)//'UPSI.dat')
       do i=1,nr
         read(nfil1,*)r(i),psi(i,:)
         read(nfil2,*)r(i),upsi(i,:)
       enddo
       do i=1,nb
         svar=maxval(abs(upsi(:,i)))
         psi(:,i)=psi(:,i)/svar
         upsi(:,i)=upsi(:,i)/svar
       enddo
       rewind(nfil1)
       rewind(nfil2)
       do i=1,nr
         write(nfil1,*)r(i),psi(i,:)
         write(nfil2,*)r(i),upsi(i,:)
       enddo
       close(nfil1)
       close(nfil2)
       deallocate(r)
       deallocate(psi)
       deallocate(upsi)
!
!      =========================================================================
!      == write gnuplot script                                                ==
!      =========================================================================
       do ib=1,nb
         line(:)=' '
         line( 1)='set terminal pdf enhanced color lw 3'
         write(string,*)ib
         string=adjustl(string)
         line(2)='set output "psi_'//trim(string)//'.pdf"'
!         write(line( 2),fmt='(a,i1,a)')'set output "psi_',ib,'.pdf"'
         line( 3)='set xrange [0:3.2]'
         line( 4)='set yrange [-2:2]'
         line( 5)='set xlabel "r[a0]"'
         line( 6)='set ylabel ""'
!         line( 7)='set data style lines'
         line( 8)='set border 3'
         line( 9)='set xtics nomirror'
         line(10)='set ytics nomirror'
         write(line(11),fmt='(a,i3,a,i1,a)')'set title "wave functions ib=' &
      &                                     ,ib,' and l=',lb(ib),'"'
         write(line(12),fmt='(a,i3,a)')'plot "AEPSI.dat" using 1:',1+ib &
      &                                ,' ls 1 with lines title "AE-psi" \ '    
         write(line(13),fmt='(a,i3,a)')',"UPSI.dat" using 1:',1+ib &
      &                                ,' ls 2 with lines title "upsi" \ '    
         line(14)=''
         open(unit=nfilgp,file=trim(tmpdir)//filename,form='formatted')
         do i=1,14
           write(nfilgp,fmt='(A)')trim(line(i))
         enddo
         close(nfilgp)
!
!        =======================================================================
!        == execute gnuplot to produce postscript file
!        =======================================================================
         command='gnuplot '//adjustl(filename)
         command=trim(precommand)//command
!        ==== call lib__system(command)
!         call lib__system(command) ! rc=system(command)
         call execute_command_line(trim(command),exitstat=rc)
         if(rc.ne.0) then
           call error__msg('command line argument failed')
           call error__chval('command',trim(command))
           call error__i4val('return code',rc)
           call error__stop('gnuplotwf')
         end if
       enddo

       return 
       end
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine writetexfile(nfil)
!      *************************************************************************
!      ** converts the input file "pot.dat" into an eps file                  **
!      *************************************************************************
       implicit none
       integer      ,intent(in) :: nfil
       integer      ,parameter  :: nline=200
       character(512)           :: line(nline)
       character(512)           :: inline(5)
       character(4),parameter   :: dbs='\\'  !double backslash
       character(64)           :: string
       integer(4)               :: nb
       integer(4)               :: lx
       integer(4)  ,allocatable :: lb(:)
       real(8)     ,allocatable :: eb(:)
       real(8)     ,allocatable :: fb(:)
       integer(4)  ,parameter   :: nvalx=10
       integer(4)               :: nval
       integer(4)               :: i,iline,l,ib
       real(8)                  :: r8arr(nvalx)
       integer(4)               :: i4arr(nvalx)
!      *************************************************************************
!
!      =========================================================================
!      ==   
!      =========================================================================
       call readlinei4('lpro.dat',nvalx,nval,i4arr)
       lx=maxval(i4arr(:nval))
!
!      =========================================================================
!      ==   
!      =========================================================================
!    replace \t \b \v
       write(nfil,*)'\documentclass{article}'
       write(nfil,*)'\usepackage{geometry}'
       write(nfil,*)'\geometry{a4paper,left=2cm,right=3cm, top=3cm, bottom=2cm}'
       write(nfil,*)'\usepackage[ansinew]{inputenc}'
       write(nfil,*)'\usepackage{amssymb}'
       write(nfil,*)'\usepackage{ngerman}'
       write(nfil,*)'\usepackage{color,graphics,graphicx,epsfig}'
       write(nfil,*)'\usepackage{fancyhdr}'
       write(nfil,*)'\pagestyle{fancy}'
       write(nfil,*)'\begin{document}'
       write(nfil,*)'\title{PAW-Report}'
       write(nfil,*)'\author{PAW mkreport-Tool}'
       write(nfil,*)'\rhead{\today}'
       write(nfil,*)'\cfoot{\thepage}'
       write(nfil,*)''


       ! line( 1)='\documentclass{article}'
       ! line( 2)='\usepackage{geometry}'
       ! line( 3)='\geometry{a4paper,left=2cm,right=3cm, top=3cm, bottom=2cm}'
       ! line( 4)='\usepackage[ansinew]{inputenc}'
       ! line( 5)='\usepackage{amssymb}'
       ! line( 6)='\usepackage{ngerman}'
       ! line( 7)='\usepackage{color,graphics,graphicx,epsfig}'
       ! line( 8)='\usepackage{fancyhdr}'
       ! line( 9)='\pagestyle{fancy}'
       ! line(10)='\begin{document}'
       ! line(11)='\title{PAW-Report}'
       ! line(12)='\author{PAW mkreport-Tool}'
       ! line(13)='\rhead{\today}'
       ! line(14)='\cfoot{\thepage}'
       ! line(15)=''
       ! do i=1,15
       !   write(nfil,fmt='(a)')trim(line(i))
       ! enddo
!
!      =========================================================================
!      ==   Report of input data                                              ==
!      =========================================================================
       line( 1)='\section{Input Data}'
       line( 2)='\begin{center}'
       line( 3)='\begin{tabular}{|l|l|}'
       line( 4)='\hline'
       call readline('id.dat',inline(1))  ! setup id
       line( 5)='Setup-id &\verb|'//trim(inline(1))//'|'//dbs
       if(inline(1)(2:2).eq.'_') inline(1)(2:2)=' '
       line( 6)='Element &'//trim(inline(1)(1:2))//dbs
       call readline('zv.dat',inline(1))
       line( 7)='Atomic number &'//trim(inline(1))//dbs
       call readline('psi.type.dat',inline(1))   
       line( 8)='Pseudization methiod &'//trim(inline(1))//dbs
       call readline('rbox.dat',inline(1))  
       line( 9)='box radius &'//trim(inline(1))//dbs
       call readliner8('rcsm.dat',nvalx,nval,r8arr)
       write(inline(1),fmt='(10f10.3)')r8arr(1:nval)
       line(10)='Rcsm/$r_{cov}$ &'//trim(inline(1))//dbs
       call readliner8('psi.rcl.dat',nvalx,nval,r8arr)  
       write(inline(1),fmt='(10f10.3)')r8arr(1:nval)
       line(11)='RCL/$r_{cov}$ &'//trim(inline(1))//dbs
       call readliner8('psi.lambda.dat',nvalx,nval,r8arr)
       write(inline(1),fmt='(10f10.3)')r8arr(1:nval)
       line(12)='$\lambda_\ell$ &'//trim(inline(1))//dbs
       call readline('pot.pow.dat',inline(1))  
       line(13)='power from potential pseudization &'//trim(inline(1))//dbs
       call readline('pot.rc.dat',inline(1))  
       line(14)='Cutoff radius for potential pseudization &'//trim(inline(1))//dbs
       call readline('core.pow.dat',inline(1))  
       line(15)='power from core pseudization &'//trim(inline(1))//dbs
       call readline('core.rc.dat',inline(1))  
       line(16)='Cutoff radius for core pseudization &'//trim(inline(1))//dbs
       line(17)='\hline'
       line(18)='\end{tabular}'
       line(19)='\end{center}'
       do i=1,19
         write(nfil,fmt='(a)')trim(line(i))
       enddo
!
!      =========================================================================
!      == number of projector functions                                       ==
!      =========================================================================
       line( 1)="\subsection{Configuration}"
       line( 2)='\begin{tabular}{ll}'
       call readline('npro.dat',inline(1))  
       line( 3)='\textbf{NPRO} & '//trim(inline(1))//dbs
       call readline('lpro.dat',inline(1))  
       line( 4)='\textbf{LPRO} & '//trim(inline(1))//dbs
       line( 5)='\end{tabular}'
       line( 6)=''
       line( 7)=''
       do i=1,7
         write(nfil,fmt='(a)')trim(line(i))
       enddo
!
!      =========================================================================
!      == atomic energy levels                                                ==
!      =========================================================================
       line( 1)='\subsection{Atomic energy levels}'
       line( 2)='\begin{center}'
       line( 3)='\begin{tabular}{|l|r|r|}'
       line( 4)='\hline'
       line( 5)='$\ell$ & $\epsilon[H]$ & Occupation'//dbs
       line( 6)='\hline'
       do i=1,6
         write(nfil,fmt='(a)')trim(line(i))
       enddo

       call readline('nb.dat',inline(1))  
       call readline('atom.f.dat',inline(2))  
       call readline('atom.e.dat',inline(3))  
       call readline('atom.l.dat',inline(4))  
       read(inline(1),*)nb
       allocate(fb(nb))
       allocate(eb(nb))
       allocate(lb(nb))
       read(inline(2),*)fb
       read(inline(3),*)eb
       read(inline(4),*)lb
       do i=1,nb
         write(nfil,'(i3,"&",f15.5,"&",f10.2,a)')lb(i),eb(i),fb(i),dbs
       enddo
!
       line( 1)='\hline'
       line( 2)='\end{tabular}'
       line( 3)='\end{center}'
       line( 4)=''
       do i=1,4
         write(nfil,fmt='(a)')trim(line(i))
       enddo
!
!      =========================================================================
!      == potential and scattering properties                                 ==
!      =========================================================================
       line(1)='\subsection{Potentials and Scattering Properties}'  
       line(2)='\begin{figure}[h!]'
       line(3)='  \centering'
       line(4)='  \begin{minipage}[b]{7.5cm}'
       line(5)='    \includegraphics[width=7.5cm]{potential}'
       line(6)='  \end{minipage}'
       line(7)='  \begin{minipage}[b]{7.5cm}'
       line(8)='    \includegraphics[width=7.5cm]{scatt}'
       line(9)='  \end{minipage}'
       line(10)='\end{figure}'
       line(11)=''
       do i=1,11
         write(nfil,fmt='(a)')trim(line(i))
       enddo
!
!      =========================================================================
!      == partial waves and projector functions                               ==
!      =========================================================================
       line( 1)='\newpage'
       line( 2)='\section{Partial waves and projector functions}'
       line( 3)='\begin{figure}[h!]'
       line( 4)='\centering'
       iline=4
       do l=0,lx
         iline=iline+1
         write(line(iline),fmt='(a,i1,a)') &
      &                  '\centerline{partial waves and projectors for l=',l,'}'
         iline=iline+1
         line( iline)='\begin{minipage}[b]{7.5cm}'
         iline=iline+1
         write(line(iline),fmt='(a,i1,a)') &
      &                               '\includegraphics[width=7.5cm]{phi_',l,'}'
         iline=iline+1
         line(iline)='\end{minipage}'
         iline=iline+1
         line(iline)='\begin{minipage}[b]{7.5cm}'
         iline=iline+1
         write(line(iline),fmt='(a,i1,a)') &
      &                               '\includegraphics[width=7.5cm]{pro_',l,'}'
         iline=iline+1
         line(iline)='\end{minipage}'
       enddo
       iline=iline+1
       line(iline)='\end{figure}'
       do i=1,iline
         write(nfil,fmt='(a)')trim(line(i))
       enddo
!
!      =========================================================================
!      == atomic wave functions                                               ==
!      =========================================================================
       line( 1)='\newpage'
       line( 2)='\section{Atomic wave functions}'
       line( 3)='\begin{figure}[h!]'
       line( 4)='\centering'
       iline=4
       do ib=1,nb
         iline=iline+1
         line( iline)='\begin{minipage}[b]{7.5cm}'
         iline=iline+1
         write(line(iline),fmt='(a,i2,a,i2,a)') &
     &        '\centerline{Atomic wave functions for shell ',ib  &
     &       ,' and  l=',lb(ib),'}'
         iline=iline+1
         write(string,*)ib
         string=adjustl(string)
         line(iline)='\includegraphics[width=7.5cm]{psi_'//trim(string)//'}'
         iline=iline+1
         line(iline)='\end{minipage}'
         if(modulo(ib,8).eq.0.and.ib.ne.nb) then
            iline=iline+1
            line(iline)='\end{figure}'
            iline=iline+1
            line(iline)='\begin{figure}[h!]'
            iline=iline+1
            line(iline)='\centering'
         end if
       enddo
       iline=iline+1
       line(iline)='\end{figure}'
       do i=1,iline
         write(nfil,fmt='(a)')trim(line(i))
       enddo
!
!      =========================================================================
!      =========================================================================
!      =========================================================================
       line(1)='\end{document}'
       do i=1,1
         write(nfil,fmt='(a)')trim(line(i))
       enddo
       deallocate(lb)
       deallocate(eb)
       deallocate(fb)
       return
       end       
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine readline(file,line)
       use stpreport_module
       implicit none
!      *************************************************************************
!      ** reads a single line from a specified file                           **
!      *************************************************************************
       character(*),intent(in) :: file
       character(*),intent(out):: line
       integer     ,parameter  :: nfil=1211
!      *************************************************************************
       open(nfil,file=trim(tmpdir)//file)
       read(nfil,fmt='(a)')line
       close(nfil)
       line=adjustl(line)
       return
       end
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine readlinei4(file,nvalx,nval,ival)
!      *************************************************************************
!      ** reads a single line from a specified file                           **
!      *************************************************************************
       use stpreport_module
       implicit none
       character(*),intent(in) :: file
       integer(4)  ,intent(in) :: nvalx
       integer(4)  ,intent(out):: nval
       integer(4)  ,intent(out):: ival(nvalx)
       integer     ,parameter  :: nfil=1211
       character(1024)         :: line
       integer(4)              :: i,ip
!      *************************************************************************
       open(nfil,file=trim(tmpdir)//file)
       read(nfil,fmt='(a)')line
       close(nfil)
!
       nval=0
       do i=1,nvalx
         line=adjustl(line)
         ip=index(line,' ')
         if(ip.eq.1) exit
         nval=nval+1
         if(nval.gt.nvalx) then
           stop 'in readlinei4'
         end if
         read(line(1:ip),*)ival(nval)
         line=line(ip:)
       enddo
      return
       end
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       subroutine readliner8(file,nvalx,nval,rval)
!      *************************************************************************
!      ** reads a single line from a specified file                           **
!      *************************************************************************
       use stpreport_module
       implicit none
       character(*),intent(in) :: file
       integer(4)  ,intent(in):: nvalx
       integer(4)  ,intent(out):: nval
       real(8)     ,intent(out):: rval(nvalx)
       integer     ,parameter  :: nfil=1211
       character(1024)         :: line
       integer(4)              :: i,ip
!      *************************************************************************
       open(nfil,file=trim(tmpdir)//file)
       read(nfil,fmt='(a)')line
       close(nfil)
!
       nval=0
       do i=1,nvalx
         line=adjustl(line)
         ip=index(line,' ')
         if(ip.eq.1) exit
         nval=nval+1
         if(nval.gt.nvalx) then
           stop 'in readlinei4'
         end if
         read(line(1:ip),*)rval(nval)
         line=line(ip:)
       enddo
      return
       end
!!$!
!!$!      ..1.........2.........3.........4.........5.........6.........7.........8
!!$       subroutine writegpaw(file)
!!$!      *************************************************************************
!!$!      ** reads a single line from a specified file                           **
!!$!      *************************************************************************
!!$       use strings_module
!!$       implicit none
!!$       character(2)   :: sy !element symbol
!!$       real(8)        :: aez  ! atomic number
!!$       integer(4)     :: nc   ! #(core states)
!!$       integer(4)     :: nv   ! #(valence states)
!!$       character(8)   :: xctype !lda or gga
!!$       character(64)  :: xcname ! identifier for the xc functional
!!$       character(64)  :: reltype !type of relativistic treatment
!!$       real(8)        :: aeekin     ! all electron kinetic energy
!!$       real(8)        :: aexc       ! all electron xc energy
!!$       real(8)        :: aeehartree ! all electron coulomb energy
!!$       real(8)        :: aeekin_core ! core kinetic energy
!!$       integer(4)     :: npro ! #(projector functions)
!!$!      *************************************************************************
!!$       write(nfil,fmt=*)-'<?xml version="1.0"?>'
!!$       write(nfil,fmt=*)-'<paw_setup version="0.6">'
!!$       write(nfil,fmt='("<!--- ",a60," --->")')-'test setup writing'
!!$       write(nfil,fmt='("<!--- ",a60," --->")')-'Units: Hartree atomic units'
!!$!
!!$!      =========================================================================
!!$!      ==  atom element                                                       ==
!!$!      =========================================================================
!!$       sy=
!!$       aez=
!!$       nc=
!!$       nv=aez-nc
!!$!      == write string =========================================================
!!$       string=-'<atom'
!!$       call xml$addvarch(string,-'symbol',sy)
!!$       call xml$addvari4(string,-'Z',nint(aez))
!!$       call xml$addvari4(string,-'core',nc)
!!$       call xml$addvari4(string,-'valence',nv)
!!$       string=trim(string)//'/>'
!!$       write(nfil,a)string
!!$!
!!$!      =========================================================================
!!$!      ==  exchange correlation                                               ==
!!$!      == LDA VWN: Vosko, Wilk, Nusair                                        ==
!!$!      == LDA PZ: Perdew-Zunger                                               ==
!!$!      == LDA PW: Perdew-Wang                                                 ==
!!$!      == GGA PBE: Perdew-Burke-Ernzerhof                                     ==
!!$!      == GGA RPBE: Hammer, Hansen, Norskov                                   ==
!!$!      == GGA PW91: Perdew Wang 91                                            ==
!!$!      =========================================================================
!!$       xctype='GGA'
!!$       xcname='PBE'
!!$!      == write string =========================================================
!!$       string=-'<xc_functional'
!!$       call xml$addvarch(string,-'type',xctype)
!!$       call xml$addvarch(string,-'name',xcname)
!!$       string=trim(string)//'/>'
!!$       write(nfil,a)string
!!$!
!!$!      =========================================================================
!!$!      ==  Generator                                                          ==
!!$!      ==  non-relativistic, scalar-relativistic, relativistic                ==
!!$!      =========================================================================
!!$       reltype=-'scalar-relativistic'
!!$!      == write string =========================================================
!!$       string=-'<generator'
!!$       call xml$addvarch(string,-'type',reltype)
!!$       call xml$addvarch(string,-'name',-'CP-PAW 2010')
!!$       string=trim(string)//'/>'
!!$       write(nfil,a)string
!!$       if(nc.eq.2) then
!!$         write(nfil,a)-'frozen core: [He]'
!!$       else if(nc.eq.10) then
!!$         write(nfil,a)-'frozen core: [He]'
!!$       else if(nc.eq.18) then
!!$         write(nfil,a)-'frozen core: [Ar]'
!!$       else if(nc.eq.36) then
!!$         write(nfil,a)-'frozen core: [Kr]'
!!$       else if(nc.eq.54) then
!!$         write(nfil,a)-'frozen core: [Xe]'
!!$       else if(nc.eq.86) then
!!$         write(nfil,a)-'frozen core: [Rn]'
!!$! here the other cases
!!$       else
!!$         call error$msg('core not identified')
!!$         call error$stop('writegpaw')
!!$       end if
!!$       write(nfil,a)-'</generator>'
!!$!
!!$!      =========================================================================
!!$!      ==  Energies                                                           ==
!!$!      =========================================================================
!!$       aeekin=
!!$       aeexc=
!!$       aeehartree=
!!$       aeekin_core=
!!$!      == write string =========================================================
!!$       string=-'<aeenergy'
!!$       call xml$addvarr8(string,-'kinetic',aeekin)
!!$       call xml$addvarr8(string,-'xc',aeexc)
!!$       call xml$addvarr8(string,-'electrostatic',aeehartree)
!!$       call xml$addvarr8(string,-'total',aeekin+aeexc+aeehartree)
!!$       string=trim(string)//'/>'
!!$       write(nfil,a)string
!!$       string=-'<core_energy'
!!$       call xml$addvarr8(string,-'kinetic',aeekin_core)
!!$       string=trim(string)//'/>'
!!$       write(nfil,a)string
!!$!
!!$!      =========================================================================
!!$!      ==  valence states                                                     ==
!!$!      =========================================================================
!!$!unclear
!!$       npro=
!!$       idpro=
!!$       lpro=
!!$       epro=
!!$       fpro=
!!$!
!!$!      == write string =========================================================
!!$       string=-'<valence_states>' 
!!$       write(nfil,a)string
!!$       do ipro=1,npro
!!$         string=-'<state' 
!!$         call xml$addvari4(string,-'n',nnofi(ipro)+llofi(ipro)+1)
!!$         if(fofi(ib).ne.0.d0) then
!!$           call xml$addvari4(string,-'f',nint(fofi(ib)))
!!$         end if
!!$         call xml$addvari4(string,-'e',nint(eofi(ib)))
!!$         write(idpro(ipro),*)ipro
!!$         call xml$addvarch(string,-'id',idpro)
!!$         string=trim(string)//'/>'
!!$         write(nfil,a)string
!!$       enddo
!!$       string=-'</valence_states>' 
!!$       write(nfil,a)string
!!$!
!!$!      =========================================================================
!!$!      ==  radial gridses                                                     ==
!!$!      =========================================================================
!!$       dex=
!!$       r1=
!!$       nr=
!!$!
!!$!      == write string =========================================================
!!$       string=-'<radial_grid'
!!$       call xml$addvarch(string,-'eq',-'r=a*(exp(d*i)-1)')
!!$       call xml$addvarr8(string,-'a',r1)
!!$       call xml$addvarr8(string,-'d',dex)
!!$       call xml$addvari4(string,-'n',nr)
!!$       call xml$addvari4(string,-'istart',0)
!!$       call xml$addvari4(string,-'iend',nr-1)
!!$       call xml$addvarch(string,-'id','shlog')
!!$       string=trim(string)//'/>'
!!$       write(nfil,a)string
!!$!
!!$!      =========================================================================
!!$!      ==  shape function for compensation charge                             ==
!!$!      =========================================================================
!!$       rcsm=
!!$!
!!$!      == write string =========================================================
!!$       string=-'<shape_function'
!!$       call xml$addvarch(string,-'type',-'gauss')
!!$       call xml$addvarr8(string,-'rc',rcsm)
!!$       string=trim(string)//'/>'
!!$       write(nfil,a)string
!!$!
!!$!      =========================================================================
!!$!      ==  radial functions                                                   ==
!!$!      =========================================================================
!!$       aecore=
!!$       pscore=
!!$       aevalencerho=
!!$       psvalencerho=
!!$       vzero=
!!$!
!!$!      == write string =========================================================
!!$       string=-'<ae_core_density'
!!$       call xml$addvarch(string,-'grid',-'shlog')
!!$       string=trim(string)//'/>'
!!$       write(nfil,a)string
!!$       write(nfil,fmt='(5e20.10)')aecore(:)                  
!!$       string=-'</ae_core_density>'
!!$       write(nfil,a)string
!!$!
!!$       string=-'<ps_core_density'
!!$       call xml$addvarch(string,-'grid',-'shlog')
!!$       string=trim(string)//'/>'
!!$       write(nfil,a)string
!!$       write(nfil,fmt='(5e20.10)')pscore(:)                  
!!$       string=-'</ps_core_density>'
!!$       write(nfil,a)string
!!$!
!!$       string=-'<ae_valence_density'
!!$       call xml$addvarch(string,-'grid',-'shlog')
!!$       string=trim(string)//'/>'
!!$       write(nfil,a)string
!!$       write(nfil,fmt='(5e20.10)')aevalencerho(:)                  
!!$       string=-'</ae_valence_density>'
!!$       write(nfil,a)string
!!$!
!!$       string=-'<ps_valence_density'
!!$       call xml$addvarch(string,-'grid',-'shlog')
!!$       string=trim(string)//'/>'
!!$       write(nfil,a)string
!!$       write(nfil,fmt='(5e20.10)')psvalencerho(:)                  
!!$       string=-'</ps_valence_density>'
!!$       write(nfil,a)string
!!$!
!!$       string=-'<zero_potential'
!!$       call xml$addvarch(string,-'grid',-'shlog')
!!$       string=trim(string)//'/>'
!!$       write(nfil,a)string
!!$       write(nfil,fmt='(5e20.10)')vzero(:)                  
!!$       string=-'</zero_potential>'
!!$       write(nfil,a)string
!!$!
!!$!      =========================================================================
!!$!      ==  partial waves and projector functions                              ==
!!$!      =========================================================================
!!$       vstateid
!!$       aephi 
!!$       psphi
!!$       pro
!!$!
!!$       do ipro=1,npro
!!$         string=-'<ae_partial_wave'
!!$         call xml$addvarch(string,-'state',idpro(ipro))
!!$         call xml$addvarch(string,-'grid',-'shlog')
!!$         string=trim(string)//'/>'
!!$         write(nfil,a)string
!!$         write(nfil,fmt='(5e20.10)')aephi(:,ipro)                  
!!$         string=-'</ae_partial_wave>'
!!$         write(nfil,a)string
!!$!
!!$         string=-'<pseudo_partial_wave'
!!$         call xml$addvarch(string,-'state',idpro(ipro))
!!$         call xml$addvarch(string,-'grid',-'shlog')
!!$         string=trim(string)//'/>'
!!$         write(nfil,a)string
!!$         write(nfil,fmt='(5e20.10)')psphi(:,ipro)                  
!!$         string=-'</pseudo_partial_wave>'
!!$         write(nfil,a)string
!!$!
!!$         string=-'<projector function'
!!$         call xml$addvarch(string,-'state',idpro(ipro))
!!$         call xml$addvarch(string,-'grid',-'shlog')
!!$         string=trim(string)//'/>'
!!$         write(nfil,a)string
!!$         write(nfil,fmt='(5e20.10)')psphi(:,ipro)                  
!!$         string=-'</projector function'
!!$         write(nfil,a)string
!!$       enddo
!!$!
!!$!      =========================================================================
!!$!      ==  kinetic energy differences                                         ==
!!$!      =========================================================================
!!$       npro=
!!$       dtkin(:,:)=
!!$!
!!$!      == write string =========================================================
!!$       string=-'<kinetic_energy_differences>'
!!$       write(nfil,a)string
!!$       write(nfil,fmt='(5e20.10)')dtkin(:,:)                  
!!$       string=-'</kinetic_energy_differences>'
!!$       write(nfil,a)string
!!$!
!!$!      =========================================================================
!!$!      ==  close down                                                         ==
!!$!      =========================================================================
!!$       string=-'</paw_setup>'
!!$       write(nfil,a)string
!!$       return
!!$       end
!!$
!!$
!!$!
!!$!      ..1.........2.........3.........4.........5.........6.........7.........8
!!$       subroutine xml$addvari4(string,id,val)
!!$!      *************************************************************************
!!$!      ** reads a single line from a specified file                           **
!!$!      *************************************************************************
!!$       implicit none
!!$       character(*),intent(in) :: id
!!$       integer(4)  ,intent(in) :: val
!!$       character(*),intent(inout):: string
!!$       character(64)           :: numstring
!!$       integer                 :: length
!!$!      *************************************************************************
!!$       pos=len_trim(string)
!!$       write(numstring,*)val
!!$       if(pos+len_trim(id)+4+len_trim(numstring).gt.len(string)) then
!!$         call error$stop('string is too small')
!!$         call error$stop('xml$advari4')
!!$       end if
!!$       string=trim(string)//' '//trim(adjustl(id))//'="'
!!$       string=trim(string)//trim(adjustl(numstring))//'"'
!!$       return
!!$       end
