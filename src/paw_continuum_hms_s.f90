!blo todo
!blo 1)there are calles from linkedlist with a double underscore instead of a dollar
!blo 2) there is a remark in rforces saying that the implentation may not be correct
! there was an error in cont_calc_ener when calculating the face-self energy
!..........................................................................
Module continuum_module
!***************************************************************************
!
!   Continuum Module
!
!  original code by Peter Margl
!  rewritten by R. Schmid  (Feb. 1999), Calgary
!  modified by Hans Martin Senn (???), Calgary
!
!***************************************************************************
!
!  changes:
!     - implement everyhting as one MODULE.
!     - all faces kept in straight arrays instead of list-structure
!      to allow parallelization
!     - calculate forces for nuclei/gaussians only when necessary in MTS
!     - only MTS implementation!
!     - calculate the position of the faces on the fly (no call of __switch
!        paw_atoms necessary
!     - allow to keep surfac-function weighted face-face-distance matrix 
!         which means a lot of memory but a significant speedup 
!
!***************************************************************************
!
!  there are NO external variables
!  only a couple of function calls are made public
!
!  on the other hand: there are no subroutines outside the module
!  so every calling scope MUST include a USE statement
! 
!***************************************************************************
!hms 
!hms Changes by H.M. Senn, May 2001
!hms 
!hms - The surface-tension term is switched off (simply by setting the 
!hms   coefficients to 0).
!hms - The face self-interaction (intra-face) term is multiplied by the surface
!hms   function. Hence, there is no self-interaction from buried faces.
!hms - The restart file is now read and written correctly
!hms 
!hms - Sep 2001: Surface tension ON again (we should provide a keyword for this!)
!hms - Oct 2001: New values for surface-tension parameters
!hms             Splining for close-lying charges
!***************************************************************************
!blo   Continuum$propagate is called from isolate object
!blo   
!blo   continuum$propagate
!blo   ->continuum_init
!blo     ->atomlist$natom(nat)
!blo     ->ISOLATE$GETI4('NFCT',ival)
!blo   ->continuum_setup
!blo     ->fermi\_function
!blo     ->dfermi\_function
!blo     ->continuum_filter
!blo     ->calculate_wdff (only if contv_keep_dff=T)
!blo   ->continuum_gforces
!blo   ->continuum_rforces
!blo   ->continuum_charge_pot
!blo   ->propagate_charge
!blo   ->continuum_register_energies
!blo   ->energylist$set('SURFACE Q EKIN',val)
!blo   ->energylist$add('constant energy',val)
!blo   ->continuum_refilter
!blo   
!***************************************************************************
implicit none
!private
!===========================================================================
!== this is the interface to the module                                   ==
!===========================================================================
public continuum$setl4
public continuum$setr8
public continuum$setr8a
public continuum$seti4
public continuum$setcha
public continuum$read
public continuum$write
public continuum$returnprotocol
public continuum$propagate
!===========================================================================
! these are the internal variables of the continuum                       ==
! they are all kept SAVEed                                                ==
! as a convention all "control" variables which have been in the          ==
! original code in teh control module get a contv_                        ==
!===========================================================================
logical,save                  :: contv_restart
logical,save                  :: contv_move_netcharges
logical,save                  :: contv_on=.false.
logical,save                  :: contv_stop
logical,save                  :: contv_move_atoms
logical,save                  :: contv_keep_dff=.true.
logical,save                  :: contv_read_old_restart=.false.
! debug                                
logical,save                  :: contv_debug_print = .true.
logical,save                  :: contv_init = .false.
logical,save                  :: contv_filter = .false.
logical,save                  :: contv_filter_init = .false.
integer,save                  :: contv_multiple
real(kind=8),save             :: contv_mass
real(kind=8),save             :: contv_timestep
real(kind=8),save             :: contv_charge_friction
real(kind=8),save             :: contv_scalingfactor
real(kind=8),save,allocatable :: contv_rsolv(:)
character(8),save,allocatable :: contv_tesselation_name(:)
!===========================================================================
!== these are the working arrays for the various routines                 ==
!== we define them here and try to save them throughout                   ==
!==                                                                       ==
!==  since this is all private it is only accesible from within this      ==
!==  module therefore I do not define the whole stuff inside the functions== 
!==  and allocate/deallocate it all the time.                             ==
!==                                                                       ==
!==  instead all arrays are defined here as allocatable and SAVEed and    ==
!==  they get allocated once in continuum_initialize                      ==
!==  They don't get deallocated and this is bad programming style, I know :-)
!==                                                                       ==
!==   FACE and related arrays                                             ==
!==                                                                       ==
!==   comment:                                                            ==
!==   the major point of this code review is to get rid of the            ==
!==   list-structure for the faces. it is kinda cool, but not practical   ==
!==   and at most points in the old code all values get remaped to        ==
!==   staright arrays. this is good for performance and necessary for     ==
!==   parallelization. since no faces are deleted or inserted during      ==
!==   runtime (then lists would be good) but this is handeled by the      ==
!==   surface function it is much easier to do it that way.               ==
!==   we keep one full integer array to maintain pointers to the          ==
!==   corresponding atom to which this face belongs this is memory intense==
!==    (could be done by "magic" formulas etc.) but easy and fast.        ==
!==   we do not need to keep the area, because it is the same for all     ==
!==   faces of one atom. we also don't need to keep the edges of the      ==
!==   vertices since they were used only for the matlab printout.         ==
!==   BUT: in order to make things fast (and easy to program :-) I will   ==
!==   keep a complete array of the original face positions (as if all     ==
!==   atoms would sit at the origin). Every time the continuum is called  ==
!==   the actual face positions are calculated by adding the vector of the==
!==   corresponding atomposition. Therfore the old call to                ==
!==   continuum$switch in paw_atoms is no longer necessary to move the    ==
!==   spheres. They always sit were the atoms are (I guess, this is just  ==
!==   more safe for future changes etc.)                                  ==
!===========================================================================
integer,save ::  natoms  !  this is a copy of the value in atomlist copied in init.
integer,save ::  nfaces  !  this is the magic number!!
integer,save ::  nfaces_full ! number of faces kept for filtering
integer,save ::  ngauss  !  copy of the external value
!===========================================================================
! FACE charges
!===========================================================================
real(kind=8),save,allocatable :: face_q(:)      ! charge at t0
real(kind=8),save,allocatable :: face_qm(:)     ! charge at t-
real(kind=8),save,allocatable :: face_q_keep(:) ! charge at t0
real(kind=8),save,allocatable :: face_qm_keep(:)! charge at t-
real(kind=8),save,allocatable :: face_qp(:)     ! charge at t+
real(kind=8),save,allocatable :: face_fq(:)     ! force on charge
real(kind=8),save,allocatable :: face_fmem(:)! force on charge memorized for MTS
integer     ,save,allocatable :: face_atom(:)! pointer to the corresponding atom
real(kind=8),save,allocatable :: face_fmem_keep(:)! force on charge memorized for MTS
integer     ,save,allocatable :: face_atom_keep(:)! pointer to the corresponding atom
integer     ,save,allocatable :: face_find(:) ! pointer of the filtered face to the original
!===========================================================================
! FACE positions
!===========================================================================
real(kind=8),save,allocatable ::  face_pos(:,:)    !(3,nfaces) actual positions
real(kind=8),save,allocatable ::  face_relpos(:,:) !(nfaces) relative position to nuclei
!===========================================================================
! per atom arrays
!===========================================================================
real(kind=8),save,allocatable ::  face_area(:)    !(nat) area of one face (1-natoms) 
integer     ,save,allocatable ::  face_peratom(:) !(nat) number of faces per atom (sum is nfaces)
!===========================================================================
! ENERGY contributions and other things
!===========================================================================
real(kind=8),save             :: cont_ener
real(kind=8),save             :: face_self_ener
real(kind=8),save             :: inter_face_ener
real(kind=8),save             :: hardness_ener
real(kind=8),save             :: face_atom_ener
real(kind=8),save             :: surface_tension_ener
real(kind=8),save             :: total_charge
real(kind=8),save             :: total_weighed_charge
real(kind=8),save,allocatable :: atom_surf(:)
real(kind=8),save,allocatable :: atom_charge(:)
real(kind=8),save,allocatable :: atom_weighed_charge(:)
!===========================================================================
! working arrays and other stuff
!===========================================================================
real(kind=8),save,allocatable :: afv(:,:,:)  ! distance matrix nuclei <-> face
real(kind=8),save,allocatable :: afd  (:,:)  ! distance matrix nuclei <-> face
real(kind=8),save,allocatable :: wdff(:,:)   ! weighed distance matrix face <-> face
real(kind=8),save,allocatable :: bigf(:)     ! product of fermi functions
real(kind=8),save,allocatable :: fermi(:,:)  ! fermi functions per atom/face
real(kind=8),save,allocatable :: dbigf(:,:,:)! derivative of bigf wrt nuc. coords.
real(kind=8),save,allocatable :: png(:)      ! potential of nuclear gaussians on the charge
real(kind=8),save,allocatable :: hfng(:,:,:) ! "hamiltonian" between faces and indiv. nuc. gaussians
!***************************************************************************
! the module continuum
contains
!
!     ......................................................................
      subroutine continuum_init
!     **                                                                  **
!     **  this is the setup of the continuum called from $propagate       **
!     **  in the first timestep.                                          **
!     **                                                                  **
!     **  it determines from the tesselation files and the solvation      **
!     **  radius                                                          **
!     **    face_relpos, the position of the face center relative to      **
!     **                 the atom position                                **
!     **    face area    The area of the face                             **
!     **  the faces are treated as a one dimensional vector.              **
!     **  and allocates a number of big arrays                            **
!     **                                                                  **
!     **  this is taken from the original code and heavily modified       **
!     **  (means simplified :-)                                           **
!     **  we read the tesselations right here and allocate everything     **
!     **  and report what we've done                                      **
!     **  (read_tesselations is no longer necessary)                      **
!     **                                                                  **
!     **  we do not need the atompositions here cause we just read        **
!     **  into face_relpos                                                **
!     **                                                                  **
!     **  as opposed to the old code we read everything in anyways,       **
!     **  because if a restart file is read all values will be            **
!     **  overwritten anyways (the stopping when reading the restart      **
!     **  is handeled by continuum$read, because we don't know the        **
!     **  face_q before that anyway)                                      **
!     **                                                                  **
!blo  **  remark(blo). It is yet unclear to me why the tesselation file   **
!blo  **    provides the corners of the faces instead of their centers    **
!     **                                                                  **
      implicit none
      integer(4)           :: nfct
      integer              :: nat, ng
      integer              :: nfil, nfilprot, ia, ifg, ip, i, iv
      integer              :: nmem, dummy
      integer, allocatable :: npoints(:)
      real(8), allocatable :: points(:,:)
      integer              :: vertex(3)
      real(8), parameter   :: MB = 1.D0/1024.D0/1024.D0
      real(8)              :: pi
!     ***********************************************************************
                                   call trace$push('continuum_init')
      pi=4.d0*atan(1.d0)
!
!     =======================================================================
!     ==  figure out the number of atoms and gaussians                     ==
!     =======================================================================
      call atomlist$natom(nat)
      CALL ISOLATE$GETI4('NFCT',NFCT)
      ng = nfct + 1

      call filehandler$unit('PROT', nfilprot)
      write (nfilprot, *)
      WRITE(nfilprot,fmt='("CONTINUUM INITIALIZATION" &
     &                    /"========================")')
!
!     =======================================================================
!     ==  set our private number of atoms                                  ==
!     =======================================================================
      natoms = nat
      ngauss = ng
      allocate(npoints(natoms))
      allocate(face_peratom(natoms))
      allocate(face_area(natoms)) 
!
!     =======================================================================
!     ==  read the first two lines of all tesselation files for all atoms  ==
!     ==  its a bit clumsy to close and reopen all the files, but we have  ==
!     ==  to know the number of faces before we can start (and it is only  ==
!     ==  for startup :-)                                                  ==
!     =======================================================================
      do ia=1,natoms
        call filehandler$unit(contv_tesselation_name(ia), nfil)
        rewind(nfil)
        read(nfil,*)npoints(ia)
        read(nfil,*)face_peratom(ia)  ! #(faces on this atom)
        call filehandler$close(contv_tesselation_name(ia))
      end do      
!
!     =======================================================================
!     ==  calculate the number of faces                                    ==
!     =======================================================================
      nfaces=sum(face_peratom)
      WRITE(nfilprot,fmt='(A30, I8)')"The total number faces is   : ", nfaces
!
!     =======================================================================
!     ==  now allocate all the permanent arrays                            ==
!     ==  note: for simplicity I allocate all these arrays once in full    ==
!     ==        size here even if filtering is switched on                 ==
!     =======================================================================
      allocate(face_q(nfaces))
      write (nfilprot, FMT='(A30, F8.4)')"Mem for q [MB]           : " &
     &                                  ,DBLE(nfaces*8*MB)
      allocate(face_qm(nfaces))
      write (nfilprot, FMT='(A30, F8.4)')"Mem for qm [MB]          : " &
     &                                  ,DBLE(nfaces*8*MB)
      if(.not.contv_filter) then
        allocate(face_qp(nfaces))
        write (nfilprot, FMT='(A30, F8.4)')"Mem for qp [MB]          : " &
     &                                    ,DBLE(nfaces*8*MB)
        allocate(face_fq(nfaces))
        write (nfilprot, FMT='(A30, F8.4)')"Mem for fq [MB]          : " &
     &                                    ,DBLE(nfaces*8*MB)
      end if
!
      allocate(face_fmem(nfaces))
      write (nfilprot, FMT='(A30, F8.4)')"Mem for fmem [MB]        : " &
     &                                  ,DBLE(nfaces*8*MB)
!
      allocate(face_atom(nfaces))
      write (nfilprot, FMT='(A30, F8.4)')"Mem for atom [MB]        : " &
     &                                  ,DBLE(nfaces*4*MB)
!
      allocate(face_pos(3,nfaces))
      write (nfilprot, FMT='(A30, F8.4)')"Mem for pos [MB]         : " &
     &                                  ,DBLE(nfaces*3*8*MB)
      allocate(face_relpos(3,nfaces))
      write (nfilprot, FMT='(A30, F8.4)')"Mem for relpos [MB]      : " &
     &                                   ,DBLE(nfaces*3*8*MB)
      write (nfilprot, FMT='(A30, F8.4)')"MEM TOTAL [MB]           : " &
     &                                   , DBLE(nfaces*(11*8+4)*MB)
!
!     =======================================================================
!     ==  allocate little stuff                                            ==
!     =======================================================================
      allocate(atom_surf(natoms))
!<<blo
      atom_surf(:)=0.d0 !used in net_report but assigned only later in setup
!>>blo
      allocate(atom_charge(natoms))
      allocate(atom_weighed_charge(natoms))
!      
!     =======================================================================
!     ==  allocate big working storage                                     ==
!     =======================================================================
      nmem = 0
      allocate(afv(3,natoms,nfaces))    ! atom-face vector
      nmem = nmem + 3*natoms*nfaces
      allocate(afd(natoms,nfaces))      ! atom-face distance
      nmem = nmem + natoms*nfaces
      allocate(bigf(nfaces))            ! product of fermis per face
      nmem = nmem + nfaces
      allocate(fermi(natoms,nfaces))    ! fermi functions per atom/face
      nmem = nmem + natoms*nfaces
      if(.not.contv_filter) then
!       __  deriv. of prod. of fermis per face wrt to nuc. coords ___________
        allocate(dbigf(3,natoms,nfaces))  
        nmem = nmem + 3*natoms*nfaces
        allocate(png(nfaces))             ! potential of nuclear gaussians
        nmem = nmem + nfaces
!       __   "hamiltonian" face <-> ng (for derivs.) _______________________
        allocate(hfng(ng, natoms,nfaces)) 
        nmem = nmem + ng*natoms*nfaces
        if(contv_keep_dff) then
          allocate(wdff(nfaces,nfaces))
          nmem = nmem + nfaces*nfaces
        end if
      else
!       __ Full INDex to "find" the charge when refiltering_________________
        allocate(face_find(nfaces))
        nmem = nmem + nfaces
      endif
      write (nfilprot, *)
      write (nfilprot, FMT='(A30, F8.4)')"Mem for working arrays [MB]: " &
     &                                  ,DBLE(nmem*8*MB)
!      
!     =======================================================================
!     ==  read the rest of the tesselation files and initialize the values ==
!     =======================================================================
      ifg = 1
      do ia = 1, natoms
        call filehandler$unit(contv_tesselation_name(ia), nfil)
        rewind(nfil)
!       __advance the two lines for nfaces and npoints_______________________
        read(nfil, *) dummy
        read(nfil, *) dummy
!       __load the points____________________________________________________
        allocate(points(3, npoints(ia)))
        do ip = 1, npoints(ia)
          read(nfil,*) (points(i, ip), i=1,3)
        end do
!       __ cycle over all faces for this atom________________________________
        do i=1,face_peratom(ia)
          read(nfil, *) iv, vertex  !identify the corner points for this face
          if(iv /= i) then
            call error$msg("error reading tesselation file")
            call error$stop("continuum_init")
            return
          end if
!         __calculate center of vertexpoints_________________________________
          face_relpos(:,ifg)=points(:,vertex(1)) &
      &                     +points(:,vertex(2)) &
      &                     +points(:,vertex(3))
!         __normalize and scale by rsolv______________________________________
          face_relpos(:,ifg)=face_relpos(:,ifg) &
      &            /dsqrt(sum(face_relpos(:,ifg)**2)) * contv_rsolv(ia)
          face_q(ifg)        = 0.0D0
          face_qm(ifg)       = 0.0D0
          face_fmem(ifg)     = 0.0D0
!         __index of corresponding atom________________________________________
          face_atom(ifg)     = ia
!         __increment main pointer___________________________________________
          ifg = ifg + 1
        end do
!       __area is 4*pi*rsolv**2 divided by the number of faces________________
        face_area(ia)=4.D0*pi*contv_rsolv(ia)**2/DBLE(face_peratom(ia))
PRINT*,'In continuum_init: ',contv_rsolv(ia),face_peratom(ia),face_area(ia)
!       __deallocate the points______________________________________________
        deallocate(points)
!       __close the tesselation file fo0r this atom (cause we are done)______
        call filehandler$close(contv_tesselation_name(ia))
      end do
      deallocate(npoints)
!
!     =======================================================================
!     ==  calculate correct timestep                                       ==
!     ==  gets overwritten from in paw.f90  !!!!!!!! grrrrrr!!!!!!         ==
!     ==  contv_timestep = contv_timestep/contv_multiple                   ==
!     =======================================================================
      call report_continuum_controls(nfilprot)

!blo      call net_report  ! uses many data that are not yet assigned
      contv_init = .true.
                                         call trace$pop
      end subroutine continuum_init
!
!     ..................................................................
      SUBROUTINE CONTINUUM$READ(NFIL,NFILO,TCHK)
!     ******************************************************************
!     **  input/output                                                **       
!     **                                                              **       
!     **  this is a modified restart read and write                   **      
!     **  because of the array structure this can be done at once, so **   
!     **  the single underscore routine s continuum_write             **
!     **  and continuum_read are gone.                                **
!     **                                                              **       
!     **  no attempts to do anyhting in parallel here!                **      
!     **  (hms: We have to fix some bugs here: The separator doesn't  **      
!     **  give the correct no. of records, and hence the restart file **      
!     **  is not handled properly.)                                   **      
!     ******************************************************************
      USE RESTART_INTERFACE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN)  :: NFIL
      INTEGER(4)            ,INTENT(IN)  :: NFILO
      LOGICAL(4)            ,INTENT(OUT) :: TCHK
!RS       new separator type to avoid reading old restart info
!hms  The no. of records in the separator has to be right!!!!!!
      TYPE (SEPARATOR_TYPE),PARAMETER    :: MYSEPARATOR &
             =SEPARATOR_TYPE(3,'CONTINUUM','NONE','FEB1999','NONE')
      TYPE (SEPARATOR_TYPE)              :: SEPARATOR, xseparator
!RS       new separator type to avoid reading old restart info
      TYPE (SEPARATOR_TYPE),PARAMETER    :: OLD_SEPARATOR &
             =SEPARATOR_TYPE(2,'CONTINUUM','NONE','JUN1997','NONE')
!
      INTEGER                            :: nfaces_, natoms_, ifc, iat, i
!hms -->
      INTEGER,ALLOCATABLE,DIMENSION(:)       :: face_peratom_ ! natoms
      INTEGER,ALLOCATABLE,DIMENSION(:)       :: face_atom_    ! nfaces
      REAL(8),ALLOCATABLE,DIMENSION(:)       :: face_area_    ! natoms
      REAL(8),ALLOCATABLE,DIMENSION(:,:)     :: face_relpos_  ! (3,nfaces)
!hms <--
! dummys for old restart
      INTEGER,ALLOCATABLE,DIMENSION(:) :: FACE_IDENT_old
      INTEGER,ALLOCATABLE,DIMENSION(:) :: FACE_ATOM_IDENT_old
      REAL(8),ALLOCATABLE,DIMENSION(:,:) :: FACE_MIDPOINT_old
      REAL(8),ALLOCATABLE,DIMENSION(:,:) :: FACE_TRI1_old
      REAL(8),ALLOCATABLE,DIMENSION(:,:) :: FACE_TRI2_old
      REAL(8),ALLOCATABLE,DIMENSION(:,:) :: FACE_TRI3_old
      REAL(8),ALLOCATABLE,DIMENSION(:) :: FACE_AREA_old
      REAL(8),ALLOCATABLE,DIMENSION(:) :: FACE_Q_old
      REAL(8),ALLOCATABLE,DIMENSION(:) :: FACE_QM_old
      REAL(8),ALLOCATABLE,DIMENSION(:) :: FACE_fmem_old
!     ******************************************************************
                                 call trace$push('CONTINUUM$READ')
      TCHK= contv_restart .and. contv_on
      ! contv_restart=T means that we read from rstrt (keyword START=F)

      if (contv_read_old_restart) then
           separator=old_separator
           xseparator=old_separator
      else
           SEPARATOR=MYSEPARATOR
           XSEPARATOR=MYSEPARATOR
      end if

      CALL restart$READSEPARATOR(XSEPARATOR,NFIL,NFILO,TCHK)

      IF(.not.TCHK) then
         call trace$pop
         RETURN
      endif

      IF(XSEPARATOR%VERSION.NE.SEPARATOR%VERSION) THEN
        CALL ERROR$MSG('VERSION NOT RECOGNIZED')
        CALL ERROR$CHVAL('ACTUAL VERSION ',XSEPARATOR%VERSION)
        CALL ERROR$CHVAL('EXPECTED VERSION ',SEPARATOR%VERSION)
        CALL ERROR$STOP('CONTINUUM$READ')
      END IF
!       
!     =======================================================================
!     == the continuum has not been initialized yet (all is read from the  ==
!     == cntl and strc to this point but the arrays haven't been allocated ==
!     == etc.) for a non restart run this can be done the first time       ==
!     == continuum$propagate is called but for a restrt we have to do      ==
!     == it right here                                                     ==
!     =======================================================================
      call continuum_init

      IF (contv_read_old_restart) THEN
! we assume that everything has been read and setup properly from STRC/CNTL
! and read only the necessary part (q, qm, fmem)

        ALLOCATE(FACE_IDENT_old(NFACES))
        ALLOCATE(FACE_ATOM_IDENT_old(NFACES))
        ALLOCATE(FACE_MIDPOINT_old(3,NFACES))
        ALLOCATE(FACE_TRI1_old(3,NFACES))
        ALLOCATE(FACE_TRI2_old(3,NFACES))
        ALLOCATE(FACE_TRI3_old(3,NFACES))
        ALLOCATE(FACE_AREA_old(NFACES))
        ALLOCATE(FACE_Q_old(NFACES))
        ALLOCATE(FACE_QM_old(NFACES))
        ALLOCATE(FACE_fmem_old(NFACES))

        read (nfil) nfaces_
        if (nfaces_ /= nfaces) then
          call error$msg('RESTART FILE AND INPUT ARE INCONSISTENT')
          call error$i4val('nfaces_',nfaces_)
          call error$i4val('nfaces',nfaces)
          call error$stop("CONTINUUM$READ OLD RESTART")
          return 
        end if

print *, "cont readin"
print *, "nfaces from file is  ", nfaces_
print *, "nfaces from setup is ", nfaces
print *, "current size of face_q ", size(face_q)
        READ(NFIL)FACE_IDENT_old,FACE_ATOM_IDENT_old,FACE_MIDPOINT_old,&
         & FACE_TRI1_old,FACE_TRI2_old,FACE_TRI3_old,FACE_AREA_old,FACE_Q_old,FACE_QM_old,face_fmem_old

        do ifc = 1,nfaces
            if (face_atom_ident_old(ifc) /= face_atom(ifc)) then
               call error$msg("face_atom in old restart and input inconsistent")
               call error$stop("continuum$read old restart")
               return
            end if
            face_q(ifc) = face_q_old(ifc)
            face_qm(ifc)= face_qm_old(ifc)
            face_fmem(ifc) = face_fmem_old(ifc) 
        end do

        DEALLOCATE(FACE_IDENT_old)
        DEALLOCATE(FACE_ATOM_IDENT_old)
        DEALLOCATE(FACE_MIDPOINT_old)
        DEALLOCATE(FACE_TRI1_old)
        DEALLOCATE(FACE_TRI2_old)
        DEALLOCATE(FACE_TRI3_old)
        DEALLOCATE(FACE_AREA_old)
        DEALLOCATE(FACE_Q_old)
        DEALLOCATE(FACE_QM_old)
        DEALLOCATE(FACE_fmeM_old)

      else  ! This is for the RS version of the continuum rstrt

        read (nfil) nfaces_
        read (nfil) natoms_
        if ((nfaces_ /= nfaces) .or. (natoms_ /= natoms)) then
          call error$msg('RESTART FILE AND INPUT ARE INCONSISTENT')
          call error$i4val('nfaces_',nfaces_)
          call error$i4val('nfaces',nfaces)
          call error$stop('CONTINUUM$READ')
          return 
        end if      
!hms -->
        !hms  We check consistency between values from input (cntl, strc) and rstrt
        ALLOCATE(face_atom_(nfaces))
        ALLOCATE(face_relpos_(3,nfaces))
        ALLOCATE(face_peratom_(natoms))
        ALLOCATE(face_area_(natoms))
        
        READ(nfil) face_q, face_qm, face_fmem, face_atom_, face_relpos_, face_peratom_, face_area_

        DO iat = 1,natoms
          IF ((face_peratom_(iat) /= face_peratom(iat)).OR.(face_area_(iat) /= face_area(iat))) THEN
            CALL error$msg('RESTART FILE AND INPUT ARE INCONSISTENT')
            call error$i4val('iat',iat)
            call error$i4val('face_peratom_',face_peratom_(iat))
            call error$i4val('face_peratom',face_peratom(iat))
            call error$r8val('face_area_',face_area_(iat))
            call error$r8val('face_area',face_area(iat))
            CALL error$stop('CONTINUUM$READ')
          ENDIF
        END DO
        DO ifc = 1,nfaces
          IF (face_atom_(ifc) /= face_atom(ifc)) THEN
            CALL error$msg('RESTART FILE AND INPUT ARE INCONSISTENT')
            call error$i4val('face_atom(ifc)',face_atom(ifc))
            call error$i4val('face_atom_(ifc)',face_atom_(ifc))
            CALL error$stop('CONTINUUM$READ')
          ENDIF
          DO i = 1,3
            IF(abs(face_relpos_(i,ifc)-face_relpos(i,ifc)).gt.1.d-8) THEN
              CALL error$msg('RESTART FILE AND INPUT ARE INCONSISTENT')
              call error$i4val('ifc',ifc)
              call error$i4val('i',i)
              call error$r8val('face_relpos(i,ifc)',face_relpos(i,ifc))
              call error$r8val('face_relpos_(i,ifc)',face_relpos_(i,ifc))
              CALL error$stop('CONTINUUM$READ')
            ENDIF
          END DO
        END DO
        DEALLOCATE(face_atom_)
        DEALLOCATE(face_relpos_)
        DEALLOCATE(face_peratom_)
        DEALLOCATE(face_area_)
!hms <--
      end if
!     
      if (contv_stop) then
            face_qm(:) = face_q(:)
      endif

      call trace$pop
      RETURN
      END subroutine continuum$read
!
!     ..................................................................
      SUBROUTINE CONTINUUM$WRITE(NFIL,NFILO,TCHK)
!     ******************************************************************
      USE RESTART_INTERFACE
      IMPLICIT NONE
      INTEGER(4)       ,INTENT(IN)  :: NFIL
      INTEGER(4)       ,INTENT(IN)  :: NFILO
      LOGICAL(4)       ,INTENT(OUT) :: TCHK
!hms  Again, the no. of records is 3, not 2 !!!!!!!!
      TYPE (SEPARATOR_TYPE),PARAMETER   :: MYSEPARATOR &
             =SEPARATOR_TYPE(3,'CONTINUUM','NONE','FEB1999','NONE')
!     ******************************************************************
      call trace$push('CONTINUUM$WRITE')
      IF(.NOT. contv_ON) RETURN
      TCHK=.TRUE.
      CALL restart$WRITESEPARATOR(MYSEPARATOR,NFIL,NFILO,TCHK)
!
!     write it right here

      write(nfil) nfaces
      write(nfil) natoms

      if (contv_filter) then
          write(nfil) face_q_keep, face_qm_keep, face_fmem_keep, face_atom_keep, face_relpos, face_peratom, face_area
      else
          write(nfil) face_q, face_qm, face_fmem, face_atom, face_relpos, face_peratom, face_area
      end if
!
      call net_report    
!
!
      call trace$pop
      RETURN
      END subroutine continuum$write

! ##########################################################################
!
!  utility functions and calculation stuff
!
! ##########################################################################
!
!     ........................................................................
      FUNCTION FERMI_FUNCTION(R,R0) RESULT(FERMI)
      IMPLICIT NONE
      REAL(8) :: R,R0,FERMI
!     ***********************************************************************
      IF (R > (R0-0.9D0)) THEN
       FERMI= 1.0D0 - DEXP(-(R-R0+0.9D0)**48.D0)
      ELSE
       FERMI= 0.D0
      ENDIF
      END FUNCTION FERMI_FUNCTION
!
!     ........................................................................
      FUNCTION DFERMI_function(R,R0) RESULT(DFERMI)
      IMPLICIT NONE
      REAL(8) :: R,R0,DFERMI, fermi
!     ***********************************************************************
      IF (R > (R0-0.9D0)) THEN
          DFERMI = 48.D0*(0.9D0-R0+R)**47.D0*(DEXP(-(R-R0+0.9D0)**48.D0))
      ELSE
       DFERMI = 0.D0
      ENDIF
      END FUNCTION DFERMI_function
!
!     .......................................................................
      SUBROUTINE PROPAGATE_CHARGE(KINETIC_ENERGY)
      implicit none
      real(kind=8), intent(out)   :: kinetic_energy
      integer                     :: i
      real(kind=8)                :: timestep
      real(8)                     :: svar1,svar2,svar3
!     ***********************************************************************
                                          call trace$push('PROPAGATE_CHARGE')
      KINETIC_ENERGY=0.D0
!     __has to be calculated here____________________________________________
!     __(for some magic reasons: paw.f90 overwrites it every step)___________
      timestep = contv_timestep/contv_multiple
!     ======================================================================= 
!     == vc0 = 2.D0/(1.D0+contv_charge_friction)
!     == vcm = -(1.D0-contv_charge_friction)/(1.D0+contv_charge_friction)
!     == vcf = (contv_timestep**2)/contv_mass/(1.D0-contv_charge_friction)
!     == ve  = contv_mass/(8.0D0 * contv_timestep**2)
!     ======================================================================= 
      if(contv_move_netcharges) then
!RS debug
if (contv_debug_print) then
  print *, "force, q, qm, qp"
end if
!RS debug
        svar1=2.d0/(1.D0+contv_CHARGE_FRICTION)
        svar2=1.d0-svar1
        svar3=timestep**2/contv_mass/(1.D0+contv_CHARGE_FRICTION)
        do i=1,nfaces 
!         __face_qp(i) = (vc0*face_q(i)) + (vcm*face_qm(i)) + (vcf*face_fq(i)) 
          FACE_QP(i)=svar1*FACE_Q(i)+svar2*FACE_QM(i)+svar3*FACE_FQ(i)
!
!blo FACE_QP(i)=(2.D0*FACE_Q(i)-(1.D0-contv_CHARGE_FRICTION)*FACE_QM(i) &
!blo&               +(FACE_FQ(i))/(contv_mass)*TIMESTEP**2) &
!blo&               /(1.D0+contv_CHARGE_FRICTION)
!RS debug
if (contv_debug_print) then
 print '(4F20.16)', face_fq(i), face_q(i), face_qm(i), face_qp(i)
end if
!RS debug
!         __ kinetic_energy = kinetic_energy + (face_qp(i)-face_qm(i))*ve
          KINETIC_ENERGY = KINETIC_ENERGY+(FACE_QP(i)-FACE_QM(i))**2
!blo          KINETIC_ENERGY = KINETIC_ENERGY+((FACE_QP(i)-FACE_QM(i))/&
!blo                &(2.D0*TIMESTEP))**2*contv_MASS*0.5D0
        end do
        KINETIC_ENERGY=0.5D0*CONTV_MASS*KINETIC_ENERGY/(2.D0*TIMESTEP)**2
      else
        face_qp = face_q
      endif
      call trace$pop
      RETURN
      END SUBROUTINE PROPAGATE_CHARGE
!
!     ........................................................................
      subroutine continuum_setup(nat,r0,ng,rc,qmad)
!     **                                                                    **
!     **  is called in every time step from continuum$propagate             **
!     **                                                                    **
!     **  1) determines the product of the fermi functions for each face    **
!     **  and its derivative with respect to atomic positions.              **
!     **  2) evaluates the surface tension term, wich is independent of     **
!     **  face-charges                                                      **
!     **  3) calculates the electrostatic potential "png" from the          **
!     **  atom-centeredgaussians from isolate with the face charges and its **
!     **  derivative "hfng" with respect to the charges qmad of the Gaussians*
!     **                                                                    **
!     **                                                                    **
      implicit none
      INTEGER,INTENT(IN)  :: NAT          !#(ATOMS)
      INTEGER,INTENT(IN)  :: NG           !#(GAUSSIANS PER ATOM)
      REAL(8),INTENT(IN)  :: R0(3,NAT)    !ATOM POSITIONS
      REAL(8),INTENT(IN)  :: RC(NG,NAT)   !DECAY LENGTHS OF ATOM CENTERED GAUSSIANS      
      REAL(8),INTENT(IN)  :: QMAD(NG,NAT) !CHARGES ATTRIBUTED TO A GIVEN GAUSSIAN
      REAL(8)             :: POS(3)
      real(8)             :: cafv(3)
      REAL(8)             :: cafd, cfermi, cbigf, erf
      integer             :: i, ifc, ia, ig, ia2, ifc_eff, iaf
      real(8)             :: k_tension, factor, dbigf_scal
      logical             :: face_skip
!     ***********************************************************************
                                        call trace$push('continuum_setup')
      if(nat /= natoms) then
        CALL error$msg('number of atoms inconsistent')
        CALL error$stop('continuum_setup')
      end if
!
!     =======================================================================
!     == calculate the actual face positions "pos" this step               ==
!     =======================================================================
                                            call timing$clockon('cont: setup')
      ifc_eff = 1                                                    
      do ifc=1,nfaces
        face_skip = .false.
!       __this is very clumsy, I know! but I want to keep the nonfiltered 
!       __algorithm
        if(contv_filter .and. contv_filter_init) then
          iaf = face_atom_keep(ifc)
        else
          iaf = face_atom(ifc)
        end if
        pos(:) = face_relpos(:,ifc) + r0(:,iaf)
!
!       =====================================================================
!       == calculate the vectors and distance to the nuclei                ==
!       == and the fermi functions                                         ==
!       == convention: atom - face => vecor points from face to atom       ==
!       ==                                                                 ==
!       ==blo this is smart and includes only those faces that are not in  ==
!       ==blo overlap region of another atom. It may be possible to exclude== 
!       ==blo faces that are far from entering an overlap region           ==
!       =====================================================================
!       __calculate bigf (product over the fermi functions)
        cbigf=1.D0
        do ia=1,natoms
          if(iaf == ia) then
!           __face belongs to atom
            cafv(:) = -face_relpos(:,ifc)
            cafd    = contv_rsolv(ia)
            cfermi  = 1.0D0
          else
!           __face does not belong to the atom
            cafv(:) = r0(:,ia)-pos(:)
            cafd    = dsqrt(sum(cafv*cafv))
            cfermi  = fermi_function(cafd, contv_rsolv(ia))
          end if
!         == 
          if(contv_filter.and.(cfermi<1.D-20)) then
            face_skip = .true.
            exit     ! out of the natoms-loop
          end if
          cbigf = cbigf * cfermi      ! product of fermi functions of other atoms
          afv(:,ia,ifc_eff) = cafv(:) ! atom position relative to face center
          afd(ia,ifc_eff)   = cafd    ! atom distance from face center
          fermi(ia,ifc_eff) = cfermi  ! fermi function tied to the atom
        end do
         
        if(.not.face_skip) then
          bigf(ifc_eff)       = cbigf  ! product of fermi functions
          face_pos(:,ifc_eff) = pos(:)
          if(contv_filter) then
!           __"inverted" logic: find points to the subarray 
!           __(is of size nafces_full)
!           __skiped charges get a zero
            face_find(ifc) = ifc_eff   ! lookup table
          end if
          ifc_eff = ifc_eff + 1
        else
!         this must be a contv_filter run (otherwise skip could not be true)
          face_find(ifc) = 0  
        end if
      end do

      if(contv_filter) then
        call continuum_filter(ifc_eff)
      end if
!
!     =======================================================================
!     == if we move the atoms: calculate the derivatives dbigf of bigf     ==
!     == with respect to the atomic positions                              ==
!     =======================================================================
      if (contv_move_atoms) then
        dbigf(:,:,:) = 0.D0
        do ifc=1,nfaces
          do ia=1,natoms
            if(face_atom(ifc) /= ia) then
!             __only if face does not belong to the atom grad on this atom is 
              dbigf_scal=1.D0
              do ia2=1,natoms
                if(ia == ia2) then
!                 __this is the atom itself: we need the derivative
                  dbigf_scal = dbigf_scal &
     &                       * dfermi_function(afd(ia,ifc), contv_rsolv(ia))
                else
                  ! this is one of the other atoms : multiply by fermi
                  dbigf_scal = dbigf_scal * fermi(ia2,ifc)
                endif
              end do
              dbigf(:,ia,ifc) = dbigf(:,ia,ifc) &
     &                        + dbigf_scal *(afv(:,ia,ifc)/afd(ia,ifc))
!             __the neg. gradient acts on the atom to which this face belongs
              dbigf(:,face_atom(ifc),ifc) = dbigf(:,face_atom(ifc),ifc)  &
     &                        - dbigf(:,ia,ifc)
            end if
          end do
        end do
      endif
!
!     =======================================================================
!     == calculate the accessible surface of each atom                     ==
!     == and the surface energy initial value from old code (Peter Margl)  ==
!     =======================================================================
!     ----- THIS IS FITTED TO R(C)=3.32AU AND R(H)=2.17AU
!     __hms   SURFACE_TENSION_ENER=0.0018276D0
!     __hms   K_TENSION=4.6647D-6
      SURFACE_TENSION_ENER=0.D0
      K_TENSION=0.D0
!
      do ia = 1,natoms
        atom_surf(ia) = 0.0D0
      end do
!
      do ifc = 1,nfaces
        atom_surf(face_atom(ifc)) = atom_surf(face_atom(ifc)) + bigf(ifc)
      end do
!     __hms  atom_surf(:) now contains the number of exposed faces for 
!     __each atom
      do ia = 1,natoms
        atom_surf(ia)  = atom_surf(ia) * face_area(ia)
!     __ hms  Now we have the surface areas in atom_surf(:)
! debug
        if(contv_debug_print) then
          print *,'atom_surface of atom ', ia, ' : ', atom_surf(ia)
        end if
! debug
        surface_tension_ener = surface_tension_ener + k_tension * atom_surf(ia)
      end do
!
!     =======================================================================
!     == calculate the potential (png) on the charges due to the           ==
!     == nuclei-gaussians. this stays fixed during the MTS propagation     ==
!     == of the charges (do the hfng here too for the derivs.)             ==
!     =======================================================================
      do ifc = 1,nfaces
        png(ifc) = 0.D0
        do ia = 1,natoms
          factor = bigf(ifc)/afd(ia,ifc)
          do ig = 1,ng
            CALL lib$erfr8(afd(ia,ifc)/rc(ig,ia), erf)
            hfng(ig, ia, ifc) = factor * erf
            png(ifc) = png(ifc) + hfng(ig,ia,ifc) * qmad(ig,ia)
          end do
        end do
      end do
!
!     =======================================================================
!     == if contv_keep_dff is on calculate the wdff matrix right now       ==
!     == it should be allocated in continuum_init already                  ==
!     =======================================================================
      if(contv_keep_dff) then
!hms    print *, "calculating wdff"
        call calculate_wdff
      end if
                                            call timing$clockoff('cont: setup')
                                            call trace$pop
      return
      end subroutine continuum_setup
!
!     .......................................................................
      subroutine continuum_filter(ifc_eff)
!     **                                                                   **
      implicit none
      integer, intent(in)      :: ifc_eff
      integer                  :: ifc, ifc_sub
      integer                  :: nfilo
!     ***********************************************************************
      call filehandler$unit('PROT', nfilo)

!    __ this is the number of effective faces considered
      if(.not. contv_filter_init) then
!       __save all the "full" values
print *, "filtering faces"
print *, " all faces :", nfaces
print *, "filtered faces :", ifc_eff -1
!print *, "filterlist"
!print *, face_find
        nfaces_full = nfaces
        allocate(face_q_keep(nfaces))
        allocate(face_qm_keep(nfaces))
        allocate(face_fmem_keep(nfaces))
        allocate(face_atom_keep(nfaces))
        face_q_keep    = face_q
        face_qm_keep   = face_qm
        face_fmem_keep = face_fmem
        face_atom_keep = face_atom
        deallocate(face_q)
        deallocate(face_qm)
        deallocate(face_fmem)
        deallocate(face_atom)
        contv_filter_init = .true.
      end if
!
!     =======================================================================
!     == here we begin with the filtering                                  ==
!     ==   ->  arrays have to be alocated first                            ==
!     ==   ->  and some values have to be copied                           ==
!     == switch face number to nfaces_eff                                  ==
!     =======================================================================
      nfaces = ifc_eff-1
      allocate(face_q(nfaces))
      allocate(face_qm(nfaces))
      allocate(face_qp(nfaces))
      allocate(face_fq(nfaces))
      allocate(face_fmem(nfaces))
      allocate(face_atom(nfaces))
      allocate(dbigf(3,natoms,nfaces))
      allocate(png(nfaces))
      allocate(hfng(ngauss,natoms,nfaces))
      if(contv_keep_dff) then
        allocate(wdff(nfaces,nfaces))
      end if

      do ifc = 1,nfaces_full
        ifc_sub = face_find(ifc)
        if (ifc_sub .ne. 0) then
!print *, "face " , ifc, "stored as face ", ifc_sub
          face_q(ifc_sub)    = face_q_keep(ifc)
          face_qm(ifc_sub)   = face_qm_keep(ifc)
          face_fmem(ifc_sub) = face_fmem_keep(ifc)
          face_atom(ifc_sub) = face_atom_keep(ifc)
        else
!print *, "face ", ifc, " skiped"
          if (face_q_keep(ifc) > 1.D-10) then
            write(nfilo, FMT='("charge # ", I5 &
    &            , " is skiped (q / qm / deltaq) ", 3F15.10)')  &
    &            ifc, face_q_keep(ifc), face_qm_keep(ifc) &
    &            ,face_q_keep(ifc)-face_qm_keep(ifc)
          end if
        end if
      end do
      return
      end subroutine continuum_filter
!
!     .......................................................................
      subroutine continuum_refilter
!     **                                                                   **
      implicit none
      integer              :: ifc, ifc_sub
!     ***********************************************************************
!     __ switch back to full nfaces
      nfaces = nfaces_full
!     __cleanup filtering by storing values back into the keep-list
      do ifc = 1,nfaces
        ifc_sub = face_find(ifc)
        if (ifc_sub .ne. 0) then
          face_q_keep(ifc)    = face_q(ifc_sub)
          face_qm_keep(ifc)   = face_qm(ifc_sub)
          face_fmem_keep(ifc) = face_fmem(ifc_sub)
        endif
      end do
!     __free the memory of the arrays with variing size
      deallocate(face_q)
      deallocate(face_qm)
      deallocate(face_qp)
      deallocate(face_fq)
      deallocate(face_fmem)
      deallocate(face_atom)
      deallocate(dbigf)
      deallocate(png)
      deallocate(hfng)
      if (contv_keep_dff) then
        deallocate(wdff)
      endif
      return
      end subroutine continuum_refilter
!
!     .......................................................................
      SUBROUTINE calculate_wdff
!     **                                                                   **
!     **  calculate weighed face2face distance matrix                      **
!     **                                                                   **
!     **  comment:                                                         **
!     **    a lot of work in this routine (actualy half of it) is spared   **
!     **    by using the symmetrical nature of the matrix,                 **
!     **    but the symmetry is broken if we work parallel, so for         **
!     **    the parallel implement. this will have to be changed           **
!     **    to facilitate this I keep the whole matrix                     **
!     **    (which is not very efficient in terms of the memory reqs.      **
!     **                                                                   **
      implicit none
      integer                    :: ifc1, ifc2
      real(kind=8)               :: self_int_fact 
      real(kind=8)               :: vect(3)
      REAL(kind=8)               :: dist,rspline
!     ***********************************************************************
!hms   self_int_fact = 2.0D0 * 1.9D0 * contv_scalingfactor
      self_int_fact = 1.9D0 * contv_scalingfactor

!     =======================================================================
!     ==hms The expression for the self-energy is                          ==
!     ==hms    E_sff = 1.9 * f(eps) * Sum_i{q_i^2 / Sqrt(s_i)}             ==
!     ==hms  Note that the expression for the diagonal elements A_ii in    ==
!     ==hms  Klamt's paper                                                 ==
!     ==hms    A_ii = 3.8 * f(eps) / Sqrt(s_i)                             ==
!     ==hms  represents TWICE the self-energy of one face. This is fine    ==
!     ==hms  if one evaluates the whole face-face energy as a double sum   ==
!     ==hms  and multiplies by 1/2 to account for the double counting:     ==
!     ==hms    E_ff = 1/2 * Sum_ij{q_i * A_ij * q_j}                       ==
!     ==hms  Here, however, we evaluate E_sff separately, so there is no   ==
!     ==hms  double counting. Hence, we have a factor 1.9 (as in PM's code)==
!     =======================================================================
      DO ifc1 = 1,nfaces
!       =====================================================================
!       == do self interaction first. this is what peters paper says       ==
!       == wdff(i,i) = (f(eps) * 1.07 * 4pi * bigf(i)**2)/sqrt(area(i))    ==
!       == but in the code he deviates from this (especially he does not   ==
!       == use bigf**2 and a formfactor = 1.9D0                            ==
!       == he also multiplies the potential here by two (consistently the  ==
!       == energy is evaluated without the factor 2.0, because             ==
!       ==   E(face-face) = 0.5 * sum(pot*q)                               ==
!       == I stay close to the original code to check the numbers          ==
!       =====================================================================
!hms      wdff(ifc1,ifc1) = self_int_fact/dsqrt(face_area(face_atom(ifc1)))
        wdff(ifc1,ifc1) = self_int_fact * bigf(ifc1)**2 &
     &                  / dsqrt(face_area(face_atom(ifc1)))
!
!       =====================================================================
!       ==hms  We include the surface fctunction here also for the         ==
!       ==hms  self-interaction term: No self-interaction from buried faces.=
!       ==hms  P.S.: This is a minor contribution, in the order of 0.1%    ==
!       ==hms  of the self-interaction energy.                             ==
!       =====================================================================
        DO ifc2 = ifc1+1, nfaces
!         ===================================================================
!         ==  now calculate the rest of the matrix                         ==
!         ==  in the parallel approach this ifc2 has to run from 1 to nfaces,
!         ==  but there are two different nfaces (one is the full length and 
!         ==  the other just a portion of the whole length of the matrix   ==
!         ===================================================================
!         ==hms Put in the splining for close-lying faces                  ==
!         ==hms (cf. parallel version). The idea is the following:         ==
!         ==hms The true 1/r potential is replaced by a linear function,   ==
!         ==hms if it becomes too big. The criterion is the average of the ==
!         ==hms respective self-interaction potentials of the two charges  ==
!         ==hms considered, expressed as a distance Rspline. The linear    ==
!         ==hms function is fit such that it matches the true potential in ==
!         ==hms value and derivativeat Rspline. We deviate here from RS by ==
!         ==hms including f(eps) in the self-int. potential.               ==
!         ===================================================================
          vect(:) = face_pos(:,ifc1)-face_pos(:,ifc2)
          dist    = dsqrt(sum(vect*vect))
          rspline=2.D0/(self_int_fact*(1.d0/SQRT(face_area(face_atom(ifc1))) &
     &                                +1.d0/SQRT(face_area(face_atom(ifc2)))))
          IF(dist.GE.rspline) THEN
            wdff(ifc2,ifc1) = (contv_scalingfactor*bigf(ifc1)*bigf(ifc2))/dist
          ELSE
            wdff(ifc2,ifc1) = contv_scalingfactor * bigf(ifc1) * bigf(ifc2) &
     &                     * (2.0D0/rspline - dist/rspline**2)
!PRINT *, "## Splining applied for pair: ", ifc1, ifc2
!print *, "## Rspline: ", rspline, " dist: ", dist
!print *, "## The normal potential would be: ", (contv_scalingfactor*bigf(ifc1)*bigf(ifc2))/dist
!PRINT *, "## The potential used is        : ", wdff(ifc2,ifc1)
          ENDIF
!         __ copy to symmetric position (skip in parallel)
          wdff(ifc1, ifc2) = wdff(ifc2, ifc1) 
        END DO
      END DO
      RETURN
      END SUBROUTINE calculate_wdff
!
!     .......................................................................
      subroutine continuum_gforces(nat,ng,pcont)
!     **                                                                   **
!     **  calculate  the force                                             **
!     **      pcont(ig,iat)=dE/dqmad(ig,iat)                               **
!     **  on the gaussians (on "qmad")                                     **
!     **  (assumes pcont to be initialized and adds to it!!!)              **
!     **                                                                   **
!     **  There seems to be some confusion about "potential" and "force".  **
!     **  What we need here is actually an (electrostatic) potential,      **
!     **  not a force. So, this is the potential acting on a nuclear       **
!     **  Gaussian due to the surface charges,  phi_I,nu = E_af/Q_I,nu     **
!     **                                                                   **
      implicit none
      integer,intent(in)    :: nat
      integer,intent(in)    :: ng
      real(8),intent(inout) :: pcont(ng,nat)
      integer               :: ifc, ia, ig
!     ***********************************************************************
                                            call timing$clockon('cont: gforces')
      do ifc = 1,nfaces
        do ia = 1,natoms
          do ig = 1,ng
            pcont(ig,ia) = pcont(ig,ia) + face_q(ifc)*hfng(ig, ia, ifc)
          end do
        end do
      end do
                                            call timing$clockoff('cont: gforces')
      return
      end subroutine continuum_gforces
!
!     .......................................................................
      subroutine continuum_rforces(nat,r0,ng,rc,qmad,cforce)
!     **                                                                   **
!     **  calculate the forces on the nuclei (this is a hard one!!!)       **
!     **                                                                   **
!     **  RS  debug                                                        **
!     **    for debugging calculate all contributions into a seperate      **
!     **    array and add them at the end                                  **
!     **                                                                   **
      implicit none
      integer     ,intent(in)    :: nat
      integer     ,intent(in)    :: ng
      real(kind=8),intent(in)    :: r0(3,nat)
      real(kind=8),intent(in)    :: rc(ng,nat)
      real(kind=8),intent(in)    :: qmad(ng,nat)
      real(kind=8),intent(inout) :: cforce(3,nat)
      integer                    :: ifc1, ifc2, ia, ia1, ia2, ia3, ig
      real(kind=8)               :: vect(3), force(3)
      real(kind=8)               :: dist, wq1, wq2, factor, bigfactor, expon
      real(kind=8)               :: two_o_spi, errfct 
      REAL(kind=8)               :: rspline, eff_dist
      real(kind=8), parameter    :: K1 = 1.0D1 ! hardness constant (from peters code)
!     ----- THIS IS FITTED TO R(C)=3.32AU AND R(H)=2.17AU
!hms   real(kind=8), parameter    :: K_TENSION=4.6647D-6
      real(kind=8), parameter    :: K_TENSION=0.D0
!RS debug
      real(kind=8), dimension(3,nat) :: cf_ffd, cf_ffsf1, cf_ffsf2, cf_afd, cf_afsf, cf_fh, cf_st
      real(8)                    :: pi
!     ***********************************************************************
                                                call timing$clockon('cont: rforces')  
      pi=4.d0*datan(1.d0)
      two_o_spi = 2.0D0/dsqrt(pi)
!
      cf_ffd(:,:)   = 0.0D0   !   face-face distance
      cf_ffsf1(:,:) = 0.0D0   !   face-face surfacefct. face1
      cf_ffsf2(:,:) = 0.0D0   !   face-face surfacefct. face2
      cf_afd(:,:)   = 0.0D0   !   face-atom distance
      cf_afsf(:,:)  = 0.0D0   !   face-atom surfacefct. face
      cf_fh(:,:)    = 0.0D0   !   face hardness
      cf_st(:,:)    = 0.0D0   !   surface-tension
!
!     =======================================================================
!     == (1):  face-face interaction                                       ==
!     ==      involves a double loop ofer the faces.                       ==
!     ==      we can use the symmetric nature of the matrix to save time   ==
!     ==      for a parallel implementation this will have to be changed   ==
!     ==                                                                   ==
!     ==  note: the term for the intra-charge interaction is not quite right=
!     ==    (we follow here the original implementation)                   ==
!     ==    the intra-charge interaction is NOT weighed by the surface fct.==
!     ==    which means something like an additional hardness term         ==
!     ==    therefore, the intra-face term is NOT geometry dependent and   == 
!     ==    does not show up here (could in principle be handled according ==
!     ==    to the hardness)                                               ==
!     =======================================================================
!     ==hms  There are no forces from the self-interaction term, even if it     
!     ==hms  is multiplied by the surface function. In the COSMO 
!     ==hms  approximation, the surface area is always taken to be 
!     ==hms  geometry-independent, and hence, there is no contribution 
!     ==hms  to the force acting on the nuclei. In PAW/COSMO, the surface 
!     ==hms  area is constant anyway, because we don't recalculate the 
!     ==hms  surface. Therefore, the force contribution vanishes both if 
!     ==hms  the face is switched on and if it is off.                                                            
!     =======================================================================
                                                     call timing$clockon('cont: ff-rforces')
      do ifc1 = 1,nfaces
        do ifc2 = ifc1+1,nfaces
          ia1 = face_atom(ifc1)
          ia2 = face_atom(ifc2)
          wq1     = bigf(ifc1)*face_q(ifc1)
          wq2     = bigf(ifc2)*face_q(ifc2)
!         ===================================================================
!         == (1a)  deriv wrt to face-face dist                             ==
!         == -> same convention as dbigf to get the right derivative       ==
!         ==    vectors (sign)                                             ==
!         ==    vect = pos(2) - pos(1)  : points from 1 to 2               ==
!         ===================================================================
          vect(:) = face_pos(:,ifc2) - face_pos(:,ifc1)
          dist    = dsqrt(sum(vect*vect))
          if (ia1 /= ia2) then
!           __hms Put in splining for short face-face distances
            rspline=2.0D0/(1.9*contv_scalingfactor*(1/SQRT(face_area(face_atom(ifc1))) &
     &                                            + 1/SQRT(face_area(face_atom(ifc2))))) 
            eff_dist= max(dist, rspline)
            force(:)=contv_scalingfactor * (wq1 * wq2 * (vect(:)/eff_dist/eff_dist/eff_dist))
!
!           __the gradients go with oposite signs on the corresponding atoms___ 
!           __(to get forces)__________________________________________________
!cforce(:,ia1) = cforce(:,ia1) - force(:)
!cforce(:,ia2) = cforce(:,ia2) + force(:)
            cf_ffd(:,ia1) = cf_ffd(:,ia1) - force(:)
            cf_ffd(:,ia2) = cf_ffd(:,ia2) + force(:)
          end if
!
!         ===================================================================
!         == (1b) now the derivs wrt to the surface functions              ==
!         ==      each face-face interaction gives one deriv for each atom ==
!         ===================================================================
          do ia = 1,natoms
            if(ia /= ia1) then
!             __derv wrt to bigf of face 1___________________________________
              force(:) = contv_scalingfactor * (wq2*face_q(ifc1)/dist) * dbigf(:,ia,ifc1)
!             __this acts on the current atom and the atom where the ________
!             __charge 1 sits
!cforce(:,ia) = cforce(:,ia) - force(:)
!cforce(:,ia1)= cforce(:,ia1)+ force(:)
              cf_ffsf1(:,ia) = cf_ffsf1(:,ia) - force(:)
              cf_ffsf1(:,ia1)= cf_ffsf1(:,ia1)+ force(:)
            end if

           if(ia /= ia2) then
!             __derv wrt to bigf of face 2____________________________________
              force(:) = contv_scalingfactor * (wq1*face_q(ifc2)/dist) * dbigf(:,ia,ifc2)
!             == this acts on the current atom and the atom where the ========
!             == charge 2 sits ===============================================
!cforce(:,ia) = cforce(:,ia) - force(:)
!cforce(:,ia2)= cforce(:,ia2)+ force(:)
              cf_ffsf2(:,ia) = cf_ffsf2(:,ia) - force(:)
              cf_ffsf2(:,ia2)= cf_ffsf2(:,ia2)+ force(:)
            end if
          end do
        end do
      end do

                                                   call timing$clockoff('cont: ff-rforces')

!     =======================================================================
!     == (2) :  derivatives of the face-atom interaction (single loop)     ==
!     =======================================================================
      do ifc1 = 1,nfaces
        ia1 = face_atom(ifc1)
        bigfactor =  0.0D0
        do ia2 = 1, natoms
!         =======================================================================
!         == (2a): derivatives of the interaction itself                       ==
!         ==  contribs from each gaussian, only if charge NOT on the atom itself=
!         =======================================================================
          factor = 0.0D0
          do ig = 1,ng
            expon = afd(ia2, ifc1)/ rc(ig,ia2)
            CALL lib$erfr8(expon, errfct)
            if (ia1 /= ia2) then
              force(:) = (bigf(ifc1) * face_q(ifc1) * qmad(ig, ia2) &
     &                 * afv(:,ia2,ifc1)/afd(ia2,ifc1))  &
     &                 * (-(errfct/(afd(ia2,ifc1)**2)) &
     &                    + (two_o_spi*exp(-expon**2))/(afd(ia2,ifc1)*rc(ig,ia2))) 
!cforce(:,ia2) = cforce(:,ia2) - force(:)
!cforce(:,ia1) = cforce(:,ia1) + force(:)
              cf_afd(:,ia2) = cf_afd(:,ia2) - force(:)
              cf_afd(:,ia1) = cf_afd(:,ia1) + force(:)
            end if
            factor = factor + qmad(ig, ia2) * errfct 
          end do  ! ends loop over the gaussians on atom site ia2

          factor = factor * face_q(ifc1)/afd(ia2,ifc1)
!
!         =======================================================================
!         == (2b): derivative wrt to the bigf (acts on all atoms ia3)          ==
!         ==       only if the third atom (the one that influences bigf) is NOT==
!         ==       the same as the one the charge sits on                      ==
!         =======================================================================
          do ia3=1,natoms
            if(ia1 /= ia3) then
              force(:) = factor * dbigf(:,ia3,ifc1)
!             cforce(:,ia3) = cforce(:,ia3) - force(:)
!             cforce(:,ia1) = cforce(:,ia1) + force(:)
              cf_afsf(:,ia3) = cf_afsf(:,ia3) - force(:)
              cf_afsf(:,ia1) = cf_afsf(:,ia1) + force(:)
            end if
          end do  ! ends loop over atom site 3 
!
!         =======================================================================
!         == (3): and add  derivatives of surface-tension and hardness energy  ==
!         =======================================================================
          if(ia1 /= ia2) then
            force(:) = k_tension * face_area(ia1) * dbigf(:,ia2,ifc1)
            cf_st(:,ia2) = cf_st(:,ia2) - force(:)
            cf_st(:,ia1) = cf_st(:,ia1) + force(:)
            force(:) = - K1*(face_q(ifc1)**2) * dbigf(:,ia2,ifc1)
            cf_fh(:,ia2) = cf_fh(:,ia2) - force(:)
            cf_fh(:,ia1) = cf_fh(:,ia1) + force(:)
          end if
        end do  ! ends loop over the atoms ia2 
      end do   ! ends loop over all the faces ifc1 (belonging to atom ia1)

      cforce(:,:) = cforce(:,:) + cf_ffd(:,:) + cf_ffsf1(:,:) + cf_ffsf2(:,:) + &
     &                       cf_afd(:,:) + cf_afsf(:,:) + cf_st(:,:) + cf_fh(:,:)
      cf_ffsf1(:,:) = cf_ffsf1(:,:) + cf_ffsf2(:,:)
      if (contv_debug_print) then
        print *, "(Final) Nuclear forces "
        print '(3F20.15)', cforce
        print *, " face-face dist contribution"
        print '(3F20.15)', cf_ffd
        print *, " face-face sf cont."
        print '(3F20.15)', cf_ffsf1
        print *, " face-atom dist contribution"
        print '(3F20.15)', cf_afd
        print *, " face-atom sf contribution"
        print '(3F20.15)', cf_afsf
        print *, " face hardness contribution"
        print '(3F20.15)', cf_fh
        print *, " surface tension contribution"
        print '(3F20.15)', cf_st
        call cont_num_deriv(natoms,nfaces,ng,r0,rc,qmad)
      end if
                                         call timing$clockoff('cont: rforces')   
      return
      end subroutine continuum_rforces
!
!     .......................................................................
      subroutine continuum_charge_pot
!     **                                                                   **
!     ** calculates the forces on the charges (means the negative of       **
!     ** the potential) (and all the charge dependent energies)            **
!     **                                                                   **
!     ** watch out for the signs (force or potential)                      **
!     ** we take the force for fq (neg. pot)                               **
!     ** but for the energy we need the potential  }8-<                    **
!     **                                                                   **
!     ** hms Here we need the forces acting on the surface charges         **
!     ** hms (which are dynamical variables), F_q_i = -del(E)/del(q_i)     **
!     **                                                                   **
      implicit none
      integer                    :: ifc1, ifc2, ia
      real(kind=8), parameter    :: K1 = 1.0D1 ! hardness constant (from peters code)
      real(kind=8)               :: cwdff, cpot
      real(kind=8)               :: self_int_fact
      real(kind=8), dimension(3) :: vect
      REAL(kind=8)               :: dist, rspline
!     ***********************************************************************
                                            call timing$clockon('cont: charge_pot')
!     ==hms   self_int_fact = 2.0D0 * 1.9D0 * contv_scalingfactor
      self_int_fact = 1.9D0 * contv_scalingfactor
!     ==hms  No factor 2 here (see above)

!     == zero things
      face_self_ener   = 0.0D0
      inter_face_ener  = 0.0D0
      hardness_ener    = 0.0D0
      face_atom_ener   = 0.0D0
      total_charge     = 0.0D0
      total_weighed_charge = 0.0D0
      atom_charge(:)   = 0.0D0
      atom_weighed_charge(:) = 0.0D0

!     ===========================================================================
!     ==  single loop (pot from gaussians and hardness)                        ==
!     ===========================================================================
      do ifc1 = 1,nfaces
!       __ set constant png contribution and energy
        face_fq(ifc1) = -png(ifc1)
        face_atom_ener = face_atom_ener + png(ifc1)*face_q(ifc1) !hms  png contains bigf
!       __add hardness potential and energy
        face_fq(ifc1) = face_fq(ifc1) &
     &                - (2.0D0 * K1 * (1.D0 - bigf(ifc1)) *face_q(ifc1))
        hardness_ener = hardness_ener &
     &                + (K1 * (1.D0 - bigf(ifc1)) * face_q(ifc1)**2 )
!       __ integrate the charges and the weiged charges first atomwise
!       __hms  "weig(ht)ed" means multiplied by Fermi function bigf
        atom_charge(face_atom(ifc1)) = atom_charge(face_atom(ifc1)) + face_q(ifc1)
        atom_weighed_charge(face_atom(ifc1))= atom_weighed_charge(face_atom(ifc1))  &
     &                                      + face_q(ifc1)*bigf(ifc1)
      end do

!     ===========================================================================
!     ==  double loop over all the charges                                     ==
!     ==  we use again the symmetry of the problem (has to be changed for a    ==
!     ==  parallel implementation for efficiancy we have two pieces of code,   ==
!     ==  one calculates cwdff on the fly, the other uses the matrix wdffe     ==
!     ===========================================================================
      IF(contv_keep_dff) THEN
!       __we use the corresponding entries of wdff
        do ifc1 = 1, nfaces
!         __do self interaction here first
!hms cpot = wdff(ifc1,ifc1) * face_q(ifc1)
          cpot = 2 * wdff(ifc1,ifc1) * face_q(ifc1)
!         __hms   Factor 2 from del(E_sff)/del(q_i)
          face_fq(ifc1) = face_fq(ifc1) - cpot
!         __this has to be multiplied by 0.5 at the end
          face_self_ener = face_self_ener + cpot * face_q(ifc1)
!         =======================================================================
!         ==  now do the "inner" loop
!         =======================================================================
          do ifc2 = ifc1+1, nfaces
            cpot = wdff(ifc1, ifc2) * face_q(ifc2)
!           __pot from ifc2 on ifc1
            face_fq(ifc1) = face_fq(ifc1) - cpot
!           __pot from ifc1 on ifc2
            face_fq(ifc2) = face_fq(ifc2) - wdff(ifc1, ifc2) * face_q(ifc1)
!           __ and the energy
            inter_face_ener = inter_face_ener + (cpot * face_q(ifc1)) 
          end do
        end do
      else
!       =========================================================================
!       == we calculate cwdff on the fly
!       =========================================================================
        DO ifc1 = 1, nfaces
!         =======================================================================
!         ==  calculate cwdff for self interaction                             ==
!         ==  this is what peters paper says                                   ==
!         ==  wdff(i,i) = (f(eps) * 1.07 * 4pi * bigf(i)**2)/sqrt(area(i))     ==
!         ==  but in the code he deviates from this (especially he does not use==
!         ==  bigf**2 and a formfactor = 1.9D0                                 ==
!         ==  he also multiplies the potential here by two (consistently the   ==
!         ==  energy is evaluated without the factor 2.0, because              ==
!         ==     E(face-face) = 0.5 * sum(pot*q)                               ==
!         ==  I stay close to the original code to check the numbers           ==
!         =======================================================================
!hms cwdff = self_int_fact/dsqrt(face_area(face_atom(ifc1)))
          cwdff = self_int_fact * bigf(ifc1)**2 / dsqrt(face_area(face_atom(ifc1)))
!         __hms  We include the multiplication by the surface fct. here as well
!         __ do self interaction here first
!hms cpot = cwdff * face_q(ifc1)
          cpot = 2 * cwdff * face_q(ifc1)
!         __hms  Again a factor 2 for the potential
          face_fq(ifc1) = face_fq(ifc1) - cpot
!         __this has to be multiplied by 0.5 at the end
          face_self_ener = face_self_ener + cpot * face_q(ifc1)
!         == now do the "inner" loop
          DO ifc2 = ifc1+1, nfaces
!           __calculate the cwdff first
            vect(:) = face_pos(:,ifc1)-face_pos(:,ifc2)
            dist    = dsqrt(sum(vect*vect))
!           __hms We include the splining
            rspline = 2.0D0 / (self_int_fact*(1/SQRT(face_area(face_atom(ifc1))) &
     &                                      + 1/SQRT(face_area(face_atom(ifc2)))))            
            IF(dist.GE.rspline) THEN
              cwdff = (contv_scalingfactor*bigf(ifc1)*bigf(ifc2))/dist
            ELSE
              cwdff = contv_scalingfactor * bigf(ifc1) * bigf(ifc2) &
     &                                     * (2.0D0/rspline - dist/rspline**2)
PRINT *, "## Splining applied for pair: ", ifc2, ifc1
print *, "## Rspline: ", rspline, " dist: ", dist
print *, "## The normal potential would be: ", (contv_scalingfactor*bigf(ifc1)*bigf(ifc2))/dist
PRINT *, "## The potential used is        : ", cwdff
            ENDIF
            cpot = cwdff * face_q(ifc2)
!           __pot from ifc2 on ifc1
            face_fq(ifc1) = face_fq(ifc1) - cpot
!           __pot from ifc1 on ifc2
            face_fq(ifc2) = face_fq(ifc2) - cwdff * face_q(ifc1)
!           __and the energy
            inter_face_ener = inter_face_ener + (cpot * face_q(ifc1)) 
          END DO
        END DO
      END IF
!
!     =======================================================================
!     == multiply by 0.5 (face-face energy double counting)
!     == because of the fact that we ran only over half of the
!     == off-diagonal elements, there is no double counting in the inter_face_ener
!     == in the parallel implmentation this will no longer be the case and we have
!     == to multiply here
!     =======================================================================
!     ==hms To wrap the factor-1/2 story up:
!     ==hms   E_sff = 1.9 * f(eps) * Sum_i{q_i^2 / Sqrt(s_i)}
!     ==hms -> self_int_fact is 1.9 * f(eps)
!     ==hms   del(E_sff)/del(q_i) = 2 * 1.9 * f(eps) * q_i / Sqrt(s_i)
!     ==hms -> This is cpot (which is not actually the potential, but the neg. force)
!     ==hms E_sff is evaluated from cpot
!     ==hms  E_sff = 1/2 * Sum_i{cpot_i * q_i}
!     ==hms -> So we need the factor 0.5 again!
!     =======================================================================
      face_self_ener = face_self_ener * 0.5D0
      cont_ener =   surface_tension_ener   &
     &          + face_self_ener         &
     &          + inter_face_ener        &
     &          + hardness_ener          &
     &          + face_atom_ener
!
!     =======================================================================
!     ==  loop over atoms to integrate the total charges                   ==
!     =======================================================================
      do ia = 1, natoms
        total_charge = total_charge + atom_charge(ia)
        total_weighed_charge = total_weighed_charge + atom_weighed_charge(ia)
      end do
                                            call timing$clockoff('cont: charge_pot')
      return
      end subroutine continuum_charge_pot
!
!     .......................................................................
      subroutine continuum_register_energies
!     **                                                                   **
!     ** registers the energies of the cont. in the energylists            **
!     **                                                                   **
      implicit none
!     ***********************************************************************
!     ----- SET ENERGY LISTS FOR PROTOCOL -----
      CALL ENERGYLIST$SET('INTER-TRIANGLE-POTENTIAL',inter_face_ener)
      CALL ENERGYLIST$SET('INTRA-TRIANGLE-POTENTIAL',face_self_ener)
      CALL ENERGYLIST$SET('ION-TRIANGLE-POTENTIAL',face_atom_ener)
      CALL ENERGYLIST$SET('CAVITY-HARDNESS-POTENTIAL',hardness_ener )
      CALL ENERGYLIST$SET('SURFACETENSION-POTENTIAL',surface_tension_ener )
      CALL ENERGYLIST$SET('TOTAL-SURFACE-POTENTIAL',cont_ener )
      call energylist$set('TOTAL-EFFECTIVE_CHARGE',total_weighed_charge )
      return
      end subroutine continuum_register_energies
!
!     .......................................................................
      SUBROUTINE CONTINUUM$PROPAGATE(NAT,R0,NG,RC,QMAD,ECONT,PCONT,CFORCE)
!     **                                                                   **
!     ** PROPAGATE                                                         **
!     **                                                                   **
!     ** this is the heart of the code: it gets called from paw_isolate    **
!     ** and does all the runtime calculations!  (it is PUBLIC)            **
!     **                                                                   **
!     ** we take it from the original code but modify it significantly     **
!     ** especially we break up the charge propagation and the calculation **
!     ** of the forces on the nuclei and the gaussian charges              **
!     ** ONLY multiple timestep is possible                                **
!     ** -> this is ensured by the fact that MULTIPLE must be an even      **
!     ** number (at least 2 :-)                                            **
!     **                                                                   **
      IMPLICIT NONE
      INTEGER,INTENT(IN)    :: NAT
      INTEGER,INTENT(IN)    :: NG
      REAL(8),INTENT(OUT)   :: ECONT
      REAL(8),intent(out)   :: CFORCE(3,nat)
      REAL(8),intent(out)   :: PCONT(ng,nat)
      REAL(8),intent(in)    :: R0(3,nat)
      REAL(8),intent(in)    :: QMAD(ng,nat)
      REAL(8),intent(in)    :: RC(ng,nat)
      real(kind=8)          :: cforcemem(3,nat)
      real(kind=8)          :: pcontmem(ng,nat)
      real(kind=8)          :: kinetic_energy
      integer               :: i, ideb
      logical(kind=4),save  :: tfirst = .true.
!RS DEBUG ONLY
      real(kind=8)          :: e_ff, e_sff, e_af, e_fh, e_st
!RS
!     ***********************************************************************
      IF (.NOT. contv_on) RETURN
                               CALL TRACE$PUSH('CONTINUUM$PROPAGATE')
      IF (TFIRST .and. (.not. contv_init)) THEN
        CALL CONTINUUM_INIT
      ENDIF
                                      call timing$clockon('cont: propagate')
      if((nat /= natoms) .or. (ng /= ngauss)) then
        call error$msg("number of atoms/gaussians in continuum are inconsistent")
        call error$stop("continuum$propagate")
      end if

      ECONT=0.D0
      PCONT(:,:)=0.D0
      PCONTMEM(:,:)=0.D0
      CFORCE(:,:)=0.D0
      CFORCEMEM(:,:)=0.D0
!
!     =======================================================================
!     == include here cosmo for comparison
!     =======================================================================
      call continuum_setup(nat,r0,ng,rc,qmad)
!
      CALL COSMO$READTESSELATION(CONTV_TESSELATION_NAME(1))
      CALL COSMO$SETR8A('RSOLV',NAT,CONTV_RSOLV)
      CALL COSMO$SETR8A('Q0',NFACES,FACE_Q)
      CALL COSMO$SETL4('ON',.TRUE.)
      CALL COSMO$SETR8('EPSILONR',80.1D0)
      CALL COSMO$INTERFACE(NAT,R0,NG,RC,QMAD,ECONT,PCONT,CFORCE)
!
      call cont_calc_ener(nat,nfaces,ng,face_q,face_pos,face_atom,bigf &
    &              ,face_area,r0,rc,qmad,contv_scalingfactor,e_ff,e_sff,e_af,e_fh,e_st)
      print *,"  e_ff   :", e_ff     !face-face interaction
      print *,"  e_sff  :", e_sff    ! face-electrostatic self interaction
      print *,"  e_af   :", e_af     ! atom face interaction
      print *,"  e_fh   :", e_fh     ! overlap region
      print *,"  e_st   :", e_st     ! surface tension
print*,'hms etot ',e_ff+e_sff+e_af+e_fh+e_st
stop
!
!     =======================================================================
!     == do multiple timesteps here (First half of cyle??)                 ==
!     =======================================================================
      call continuum_setup(nat,r0,ng,rc,qmad)
      call continuum_gforces(nat, ng, pcontmem)
      if (contv_move_atoms) then
        call continuum_rforces(nat, r0, ng, rc, qmad, cforcemem)
      end if
!RS  BIG DEBUG
if (contv_debug_print ) then
   print *, "pcont"
   print *, pcontmem
   print *, "cforce"
   print *, cforcemem
endif
!RS
!     =======================================================================
!     ==            first HALF OF CYCLE                                    ==
!     =======================================================================
      FIRST_HALF: DO I=1,contv_MULTIPLE/2
        if(tfirst .and. (.not. contv_restart)) exit
        call continuum_charge_pot
!RS DEBUG
        if(contv_debug_print) then
          call cont_calc_ener(nat,nfaces,ng,face_q,face_pos,face_atom,bigf &
    &              ,face_area,r0,rc,qmad,contv_scalingfactor,e_ff,e_sff,e_af,e_fh,e_st)
          print *,"  e_ff   :", e_ff
          print *,"  e_sff  :", e_sff
          print *,"  e_af   :", e_af
          print *,"  e_fh   :", e_fh
          print *,"  e_st   :", e_st
        end if
!RS DEBUG
        if(i == 1) then
          face_fq = (face_fq + face_fmem)/2.0D0
        endif
!
        CALL PROPAGATE_CHARGE(kinetic_energy)
!       == SWITCH CHARGES =================================================== 
        FACE_QM = FACE_Q
        FACE_Q  = FACE_QP
      END DO FIRST_HALF
!
!     =======================================================================
!     ==            SECOND HALF OF CYCLE                                   ==
!     =======================================================================
      SECOND_HALF: DO I=1,contv_MULTIPLE/2
        call continuum_charge_pot
!RS DEBUG
        if (contv_debug_print) then
          call cont_calc_ener(nat, nfaces, ng, face_q, face_pos, face_atom, bigf, face_area, &
      &                             r0, rc, qmad, contv_scalingfactor, e_ff, e_sff, e_af, e_fh, e_st)
          print *,"  e_ff   :", e_ff
          print *,"  e_sff  :", e_sff
          print *,"  e_af   :", e_af
          print *,"  e_fh   :", e_fh
          print *,"  e_st   :", e_st
        end if
!RS DEBUG
!
        if (i == 1) then
           call continuum_register_energies
        endif
!
        CALL PROPAGATE_CHARGE(kinetic_energy)
!
        IF (I == 1) THEN
            CALL ENERGYLIST$SET('SURFACE Q EKIN',KINETIC_ENERGY)
        ENDIF
!
!       ===== SWITCH CHARGES
        face_qm = face_q
        face_q  = face_qp
      END DO SECOND_HALF
!
!     =======================================================================
!     ===== GET FINAL FORCES AND POTENTIAL ON NUCLEI
!     =======================================================================
      call continuum_gforces(nat, ng, pcont)
      if (contv_move_atoms) then
        call continuum_rforces(nat, r0, ng, rc, qmad, cforce)
      endif
!
!     ===== STORE MEMORY FORCE ON CHARGES ----
!
!     =======================================================================
!     ==  to be consistent with peters code we have to call                ==
!     ==  continuum_charge_pot again because the charges have been         ==
!     ==  propagated and he memorizes the force from NOW                   ==
!     =======================================================================
      call continuum_charge_pot
!RS DEBUG
      if (contv_debug_print) then
        call cont_calc_ener(nat, nfaces, ng, face_q, face_pos, face_atom, bigf, face_area, &
      &                             r0, rc, qmad, contv_scalingfactor, e_ff, e_sff, e_af, e_fh, e_st)
        print *,"  e_ff   :", e_ff
        print *,"  e_sff  :", e_sff
        print *,"  e_af   :", e_af
        print *,"  e_fh   :", e_fh
        print *,"  e_st   :", e_st
      end if
!RS DEBUG
      FACE_FMEM=FACE_FQ
      IF ((.NOT. TFIRST) .OR. contv_RESTART) THEN
        PCONT=(PCONT+PCONTMEM)/2.D0
        CFORCE=(CFORCE+CFORCEMEM)/2.D0
      ENDIF
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',KINETIC_ENERGY)
      TFIRST=.FALSE.
      econt = cont_ener
!
!RS  BIG DEBUG
if (contv_debug_print ) then
   print *, "pcont"
   print *, pcontmem
   print *, "cforce"
   print *, cforcemem
endif
!RS
      if (contv_filter) then
        call continuum_refilter
      endif
                          call timing$clockoff('cont: propagate')
                          CALL TRACE$POP
      RETURN 
      END SUBROUTINE CONTINUUM$PROPAGATE
!
!****************************************************************************
!**                                                                        **
!**                              INTERFACE                                 **
!**                                                                        **
!** routines to set private variables                                      **
!** these routines are taken from the original code and modified           **
!**                                                                        **
!** all varibles are private to the module but SAVEed and start with contv_**
!**                                                                        **
!****************************************************************************
!     .......................................................................
      SUBROUTINE CONTINUUM$SETL4(ID_,VALUE_)
!     **                                                                   **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: ID_
      LOGICAL(4),intent(in)   :: VALUE_
!     ***********************************************************************
      IF (ID_ == 'MOVE') THEN
        contv_move_atoms = value_
      ELSE IF (ID_ == 'ON') THEN
        contv_on         = value_
      ELSE IF (ID_ == 'START') THEN
        contv_restart    = .not. value_
      ELSE IF (ID_ == 'FREEZE') THEN
        contv_move_netcharges = .not. value_
      ELSE IF (ID_ == 'STOP') THEN
        contv_stop = value_
      ELSE if (id_ == 'KEEP_DFF') then
        contv_keep_dff = value_
        print *,'contv_keep_dff set to value : ', value_
      else if (id_ == 'DEBUG_PRINT') then
        contv_debug_print = value_
      else if (id_ == 'READ_OLD_RESTART') then
        contv_read_old_restart = value_
        print *,'contv_read_old_restart set to value : ', value_
      else if (id_ == 'EXCLUDEOVERLAP') then
         print *, 'this option is no longer supported in the new continuum'
      else if (id_ == 'FILTER_FACES') then
        contv_filter = value_
      else
        call error$msg('id not recognized')
        call error$chval('id',id_)
        call error$stop('continuum$setl4')
      ENDIF
      return
      END SUBROUTINE CONTINUUM$SETL4
!
!     .......................................................................
      SUBROUTINE CONTINUUM$SETR8(ID_,VALUE_)
!     **                                                                   **
      IMPLICIT NONE
      CHARACTER(*) :: ID_
      REAL(8)      :: VALUE_
!     ***********************************************************************
      IF (ID_ == 'TIMESTEP') THEN
        contv_timestep = value_
      ELSE IF (ID_ == 'MASS') THEN
        contv_mass = value_
      ELSE IF (ID_ == 'EPSILON') THEN
        contv_scalingfactor = VALUE_/(VALUE_-1.D0)
!hms The scaling factor is set to eps/(eps-1). This is the inverse of the usual definition!!!
      ELSE IF (ID_ == 'FRICTION') THEN
        contv_charge_friction = value_
      ELSE 
        call error$msg('id not recognized')
        call error$chval('id',id_)
        call error$stop('continuum$setr8')
      ENDIF
      END SUBROUTINE CONTINUUM$SETR8
!
!     .......................................................................
      SUBROUTINE CONTINUUM$SETR8A(ID_,LENGTH_,VALUE_)
!     **                                                                   **
      IMPLICIT NONE
      CHARACTER(*) :: ID_
      INTEGER(4)   :: LENGTH_
      REAL(8),DIMENSION(LENGTH_)    :: VALUE_
      INTEGER(4)               :: COUNTER
!     ***********************************************************************
      IF (ID_ == 'RAD') THEN
        allocate(contv_rsolv(length_))
        contv_rsolv(:) = value_(:)
      ELSE
        call error$msg('id not recognized')
        call error$chval('id',id_)
        call error$stop('continuum$setr8a')
      ENDIF
      return
      END SUBROUTINE CONTINUUM$SETR8A
!
!     .......................................................................
      SUBROUTINE CONTINUUM$SETCHA(ID_,LENGTH_,VALUE_)
!     **                                                                   **
      IMPLICIT NONE
      CHARACTER(*)  :: ID_
      INTEGER(4)    :: LENGTH_
      CHARACTER(*)  :: VALUE_(LENGTH_)
      INTEGER(4)    :: I,NAT
!     ***********************************************************************
      IF (ID_ == 'TESS') THEN
        ALLOCATE(contv_tesselation_name(LENGTH_))
        CALL ATOMLIST$NATOM(NAT)
        IF (NAT /= LENGTH_) THEN
            CALL ERROR$MSG('SIZE INCONSISTENT')
            CALL ERROR$CHVAL('ID_',ID_)
            CALL ERROR$STOP('CONTINUUM$SETCHA')
        END IF
        DO I=1,LENGTH_
             WRITE(contv_tesselation_name(I),FMT='("TESS",I4)')I
             CALL FILEHANDLER$SETFILE(contv_tesselation_name(I),.FALSE.,VALUE_(I))
             CALL FILEHANDLER$SETSPECIFICATION(contv_tesselation_name(I),'FORM','FORMATTED')
        ENDDO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('CONTINUUM$SETCHA')
      END IF
      return
      END SUBROUTINE CONTINUUM$SETCHA
!
!     .......................................................................
      SUBROUTINE CONTINUUM$SETi4(ID_,VALUE_)
!     **                                                                   **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      integer(4)  ,intent(in) :: VALUE_
!     ***********************************************************************
      IF(ID_ == 'MULTIPLE') THEN
        contv_multiple = value_
        if(modulo(contv_multiple, 2) /= 0) then
          call error$msg("Odd multiple in Continuum not allowed!")
          call error$stop("continuum$seti4")
        endif
      ELSE
        call error$msg('id not recognized')
        call error$chval('id',id_)
        call error$stop('continuum$seti4')
      ENDIF
      return
      END SUBROUTINE CONTINUUM$SETi4
!RS
!RS  
!
!     .......................................................................
      SUBROUTINE REPORT_CONTINUUM_CONTROLS(NFILO)
!     **                                                                   **
!     **  report cont controls (called from initialize)                    **
!     **  modified from orig. code                                         **
!     **                                                                   **
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NFILO
!     ***********************************************************************
      WRITE(NFILO,*)
      WRITE(NFILO,fmt='("CONTINUUM CONTROLS"/"==================")')
      WRITE(NFILO,fmt='(" New version (R. Schmid, Feb 1999, Calgary)")')
      WRITE(NFILO,fmt='(" Modified    (H.M. Senn, May 2001, Calgary)")')
      IF (contv_on) THEN
         WRITE(NFILO,fmt='("CONTINUUM MODELING OF ELECTROSTATIC SOLVATION IS -> ON")')
      else
         WRITE(NFILO,fmt='("CONTINUUM MODELING OF ELECTROSTATIC SOLVATION IS -> OFF")')
      endif
!
      IF (contv_restart) THEN
         WRITE(NFILO,fmt='("MESH DATA READ FROM RESTART FILE")')
      ELSE
         WRITE(NFILO,fmt='("MESH GENERATED AROUND IONS READ FROM INPUT")')
      ENDIF
!
      IF (contv_MOVE_NETCHARGES) THEN
         WRITE(NFILO,fmt='("CHARGES ON MESH ARE MOBILE")')
      ELSE
         WRITE(NFILO,fmt='("CHARGES ON MESH ARE FIXED")')
      ENDIF
!
      WRITE(NFILO,fmt='("MULTIPLE TIMESTEP INTEGRATION FOR SURFACE CHARGES IS ON")')
      WRITE(NFILO,fmt='(A25,I5)')"OVERSAMPLING FACTOR : ", contv_MULTIPLE
      WRITE(NFILO,fmt='(A25,F10.5)')"EFFECTIVE TIMESTEP  : ",contv_timestep/dble(contv_multiple)
!      IF (FORCE_ZERO_OVERLAP) THEN
!       WRITE(NFILO,*)'ZERO FACE CHARGE ENFORCED IN OVERLAP REGION'
!      ENDIF
!      IF (contv_AUTO) THEN
!       WRITE(NFILO,*)'CONVERGENCE SPED UP BY VARIABLE FRICTION'
!      ENDIF
      WRITE(NFILO,FMT='("DIELECTRIC SCALING (EPS-1)/EPS    :",F12.5)')&
      & 1/contv_SCALINGFACTOR
      IF (contv_MOVE_ATOMS) THEN
         WRITE(NFILO,fmt='("FORCES ON NUCLEI DUE TO SOLVATION WILL BE COMPUTED")')
      ENDIF
      WRITE(NFILO,FMT='("FICTITIOUS MASS OF CHARGE DENSITY :",F12.5)')&
      & contv_MASS
!      IF (NOSE_THERMOSTAT) THEN
!       WRITE(NFILO,FMT='(" NOSE THERMOSTAT FOR CHARGES IS ACTIVE")')
!       IF (STOP_NOSE_THERMOSTAT) THEN
!        WRITE(NFILO,FMT='(" NOSE VARIABLE IS STOPPED")')
!       ENDIF
!       WRITE(NFILO,FMT='(" AVERAGE KINETIC ENERGY OF CHARGES:",F12.6)&
!        &')AVERAGE_KINETIC
!       WRITE(NFILO,FMT='(" FREQUENCY OF NOSE VARIABLE       :",F12.6)&
!        &')NOSE_FREQUENCY
!      ENDIF
      if (contv_keep_dff) then
         WRITE(NFILO,fmt='("THE WEIGHTED FACE-FACE DISTANCE MATRIX WILL BE KEPT")')
      else
         WRITE(NFILO,fmt='("THE WEIGHTED FACE-FACE DISTANCE MATRIX IS RECALCULATED EACH TIME")')
      endif
      if (contv_filter) then
         WRITE(NFILO,fmt='("UNEXPOSED FACES ARE FILTERED OUT")')
      else
         WRITE(NFILO,fmt='("ALL FACES ARE CONSIDERED")')
      endif
      WRITE(NFILO,*)
!
      END SUBROUTINE REPORT_CONTINUUM_CONTROLS
!
! ######################################################################
!     .......................................................................
      subroutine net_report
!     **                                                                   **
      implicit none
      integer      :: nfilo, i
      real(8)      :: total_surf, total_surf_perc, surf_perc
!hms -->
       character(16):: name(natoms)
!hms <--
!     ***********************************************************************
      call filehandler$unit('PROT',nfilo)
      write (nfilo, *)
      WRITE(NFILO,fmt='("CONTINUUM REPORT"/"================")')
      write (nfilo, FMT='(A16,1X, A12,1X, A12,1X, A17,1X, A12,1X, A12)') &
     &                   "Name", "Rsolv [au]", "Surface", "Eff. Surf. [%]" &
     &                   , "Eff. chrg.", "Charge"

      total_surf=0.0D0
      total_surf_perc = 0.0D0
      CALL atomlist$getcha('NAME',0,natoms,name(:))
print*,'marke1 net_report'
      do i=1,natoms
print*,'iata ',i,atom_surf(i)
print*,'iatb ',i,face_area(i)
print*,'iatc ',i,face_peratom(i)
        total_surf = total_surf + atom_surf(i)
        surf_perc  = atom_surf(i)/(face_area(i)*face_peratom(i))*100.D0
        total_surf_perc = total_surf_perc + surf_perc
!hms     write (nfilo, FMT='(I4,1X, F12.3,1X, F12.5,3X, F14.1,1x, F12.5,1X, F12.5)') &
!hms        &    i, contv_rsolv(i), atom_surf(i), surf_perc, atom_weighed_charge(i), atom_charge(i)
print*,'iatd ',i,contv_rsolv(i)
print*,'iatf ',i,atom_weighed_charge(i)
print*,'iatg ',i,atom_charge(i)
print*,'marke1a net_report'
        write (nfilo, FMT='(A16,1X,F12.3,1X,F12.5,3X,F14.1,1x,F12.5,1X,F12.5)') &
     &    name(i),contv_rsolv(i),atom_surf(i),surf_perc,atom_weighed_charge(i), atom_charge(i)
print*,'marke1b net_report'
      end do
      write (nfilo,*) ' ' 
print*,'marke2 net_report'
print*,'total_surf',total_surf
      write (nfilo,FMT='(A18, F12.5,1X, F12.1,1X, F12.5,1X, F12.5)') &
     &    "SUM:", total_surf, total_surf_perc/natoms, total_weighed_charge, total_charge
      write (nfilo, *) ' '
      return
      end subroutine net_report 

! #######################################################################

! called during energy printout
! taken from original code and slightly modified
! this routine is PUBLIC!

!     .......................................................................
      SUBROUTINE CONTINUUM$RETURNPROTOCOL(ON,ESOLV,EKINQ,QFRIC,QTOT)
!     **                                                                   **
      IMPLICIT NONE
      LOGICAL ,INTENT(OUT) :: ON
      REAL(8),INTENT(OUT) :: QFRIC,ESOLV,EKINQ,QTOT
!      REAL(8) :: F2,F3,F4
!      REAL(8) :: NETCHARGE,INSIDE_CHARGE,POSITIVE_INTEGRAL,&
!       &NEGATIVE_INTEGRAL,INSIDE_PLUS,INSIDE_MINUS
      INTEGER :: NFILO
!     ***********************************************************************
!
      ON = contv_on
      IF (.NOT. ON) THEN
       ESOLV=0.D0
       EKINQ=0.D0
       QFRIC=0.D0
       RETURN
      ENDIF
!      CALL FILEHANDLER$UNIT('CONTINUUM_PROTOCOL',NFILO)
      QFRIC = contv_charge_friction
      CALL ENERGYLIST$RETURN('TOTAL-SURFACE-POTENTIAL',ESOLV)
      CALL ENERGYLIST$RETURN('SURFACE Q EKIN',EKINQ)
!      CALL GET_CHARGE(NETCHARGE,INSIDE_CHARGE,POSITIVE_INTEGRAL,&
!       &NEGATIVE_INTEGRAL,INSIDE_PLUS,INSIDE_MINUS)
!      QTOT=NETCHARGE
       QTOT = total_weighed_charge
!      WRITE(NFILO,FMT='(6(1X,F12.6))')NETCHARGE,INSIDE_CHARGE,&
!       &POSITIVE_INTEGRAL,NEGATIVE_INTEGRAL,INSIDE_PLUS,INSIDE_MINUS
!
      END SUBROUTINE CONTINUUM$RETURNPROTOCOL
!

! ######################################################################
! thats it

! these are sume utility functions for debugging purposes
!
!     .......................................................................
      subroutine cont_calc_bigf(na, nf, ra, rf, fba, bigf)
!     **                                                                   **
  implicit none
  integer, intent(in)                        :: na, nf
  integer, dimension(nf), intent(in)         :: fba
  real(kind=8), dimension(3,na), intent(in)  :: ra 
  real(kind=8), dimension(3,nf), intent(in)  :: rf 
  real(kind=8), dimension(nf), intent(out)   :: bigf
  real(kind=8), dimension(3)                 :: afv
  real(kind=8)                               :: afd
  integer                                    :: ifc, ia
!     ***********************************************************************

!    calculate bigf (product over the fermi functions)
!    bigf is either 1 or 0, depending on whether the face is exposed or buried within
!    another atom's sphere
  do ifc = 1,nf
     bigf(ifc) = 1.0D0
     do ia = 1,na
        if (.not.(fba(ifc) == ia)) then
             ! face belongs not to atom
             afv(:)        = ra(:,ia)-rf(:,ifc)
             afd           = dsqrt(sum(afv*afv))
             bigf(ifc) = bigf(ifc) * fermi_function(afd, contv_rsolv(ia))
        end if 
     end do
  end do
 
      end subroutine cont_calc_bigf
!
!     .......................................................................
      subroutine cont_calc_ener(na,nf,ng,q,qpos,fba,bigf,fa,r0,rc,qmad,feps &
     &                         ,e_ff,e_sff,e_af,e_fh,e_st)
!     **                                                                   **
!     **  this subroutine calculates all energy contributions of the       **
!     **  current implementation in a straightforward way (i.e. NO         **
!     **  attempts to make it efficient) as simple as possible             **
!     **   => this is to verify things not to calculate things             **
!     **  we get ALL informations as parameters!                           **
!     **                                                                   **
      implicit none
      integer, intent(in)                             :: na, nf, ng
      real(kind=8), dimension(nf), intent(in)         :: q, bigf
      integer, dimension(nf), intent(in)              :: fba
      real(kind=8), dimension(na), intent(in)         :: fa
      real(kind=8), dimension(3,nf), intent(in)       :: qpos
      real(kind=8), dimension(3,na), intent(in)       :: r0
      real(kind=8), dimension(ng,na), intent(in)      :: rc, qmad
      real(kind=8), intent(in)                        :: feps
      real(kind=8), intent(out)                       :: e_ff,e_sff,e_af,e_fh,e_st
    
      real(kind=8)                                    :: k_tension, k1, dist, errfct
      real(kind=8), dimension(3)                      :: vect
      integer                                         :: ifc, ifc1, ifc2, ia, ig
!     ***********************************************************************
!
!     =======================================================================
!     ==  FACE FACE  INTERACTION (INTER)
!     =======================================================================
      e_ff = 0.0D0
      do ifc1 = 1,nf
        do ifc2 = 1,nf
          if (ifc1 /= ifc2) then
            vect(:) = qpos(:,ifc1) - qpos(:,ifc2)
            dist = dsqrt(sum(vect*vect))
            e_ff = e_ff + (q(ifc1) * bigf(ifc1) * q(ifc2) * bigf(ifc2))/dist
          end if
        end do
      end do
      e_ff = e_ff * feps * 0.5D0
!
!     =======================================================================
!     ==  FACE FACE INTERACTION (INTRA)
!     =======================================================================
      e_sff = 0.0D0
      do ifc = 1,nf
!blo     e_sff = e_sff + ((q(ifc)**2)/dsqrt(fa(fba(ifc))))
        e_sff = e_sff + ((bigf(ifc)*q(ifc))**2/dsqrt(fa(fba(ifc))))
      end do
!blo      e_sff = e_sff * 1.90D0 * feps
      e_sff = e_sff * 1.07d0*sqrt(4.d0*datan(1.d0)) * feps
!  
!     =======================================================================
!     ==  ATOM FACE INTERACTION
!     =======================================================================
      e_af = 0.0D0
      do ifc = 1,nf
        do ia = 1,na
          vect(:) = qpos(:,ifc) - r0(:,ia)
          dist = dsqrt(sum(vect*vect))
          do ig = 1,ng
            CALL lib$erfr8(dist/rc(ig,ia), errfct)
            e_af = e_af + (bigf(ifc) * q(ifc) * qmad(ig, ia) * errfct / dist) 
          end do
        end do 
      end do
!
!     =======================================================================
!     ==  FACE HARDNESS
!     =======================================================================
      K1 = 1.0D1
      e_fh = 0.0D0
      do ifc = 1,nf
        e_fh = e_fh + k1 * (1.0D0 - bigf(ifc)) * (q(ifc)**2)
      end do
!
!     =======================================================================
  !   SURFACE TENSION
!     =======================================================================
  !   ----- THIS IS FITTED TO R(C)=3.32AU AND R(H)=2.17AU
!hms   e_st = 0.0018276D0
!hms   K_TENSION=4.6647D-6
      e_st = 0.D0
      K_TENSION = 0.D0
      do ifc = 1,nf
        e_st = e_st + k_tension * bigf(ifc) * fa(fba(ifc))
      end do
      return
      end subroutine cont_calc_ener
!
!     .......................................................................
      subroutine cont_num_deriv(na, nf, ng, r0, rc, qmad)
!     **                                                                   **
      implicit none
      integer     ,intent(in)      :: na, nf, ng
      real(kind=8),intent(in)      :: r0(3,na)
      real(kind=8),intent(in)      :: rc(ng,na)
      real(kind=8),intent(in)      :: qmad(ng,na)
      integer                      :: ia, id, ifc
      real(kind=8)                 :: delta
! local stuff 
      real(kind=8), dimension(3,na):: r0_d
      real(kind=8), dimension(3,nf):: face_pos_d
      real(kind=8), dimension(nf)  :: bigf_d
!
      real(kind=8)                 :: e_0, e_ff_0, e_sff_0, e_fa_0, e_fh_0, e_st_0
      real(kind=8)                 :: e_1, e_ff_1, e_sff_1, e_fa_1, e_fh_1, e_st_1
      real(kind=8), dimension(3,na):: e_d, e_ff_d, e_sff_d, e_fa_d, e_fh_d, e_st_d
      real(kind=8), dimension(3,na):: e_dr, e_ff_dr, e_sff_dr, e_fa_dr, e_fh_dr, e_st_dr
      real(kind=8), dimension(3,na):: e_ds, e_ff_ds, e_sff_ds, e_fa_ds, e_fh_ds, e_st_ds
!     ***********************************************************************
      delta = 1.0D-7

! calculate iinitial values
  call cont_calc_ener(na, nf, ng, face_q, face_pos, face_atom, bigf, face_area, &
      &               r0, rc, qmad, contv_scalingfactor, e_ff_0, e_sff_0, e_fa_0, e_fh_0, e_st_0)
  e_0 = e_ff_0 + e_sff_0 + e_fa_0 + e_fh_0 + e_st_0

  do ia = 1,na
     do id = 1,3
        ! move it over
        r0_d(:,:)   = r0(:,:)
        r0_d(id,ia) = r0_d(id,ia) + delta
        ! calculate new charge positions
        do ifc = 1,nf
            face_pos_d(:,ifc) = face_relpos(:,ifc) + r0_d(:,face_atom(ifc))
        end do
        ! calculate new bigf
        call cont_calc_bigf(na, nf, r0_d, face_pos_d, face_atom, bigf_d)

        ! calculate fully displaced values
        call cont_calc_ener(na, nf, ng, face_q, face_pos_d, face_atom, bigf_d, face_area, &
      &               r0_d, rc, qmad, contv_scalingfactor, e_ff_1, e_sff_1, e_fa_1, e_fh_1, e_st_1)
        e_1 = e_ff_1 + e_sff_1 + e_fa_1 + e_fh_1 + e_st_1
        e_d(id,ia) = (e_1 - e_0)/delta  
        e_ff_d(id,ia)  = (e_ff_1  - e_ff_0)/delta  
        e_sff_d(id,ia) = (e_sff_1 - e_sff_0)/delta  
        e_fa_d(id,ia)  = (e_fa_1 - e_fa_0)/delta  
        e_fh_d(id,ia)  = (e_fh_1 - e_fh_0)/delta  
        e_st_d(id,ia)  = (e_st_1 - e_st_0)/delta  

        ! calculate with displaced points but bigf retained
        call cont_calc_ener(na, nf, ng, face_q, face_pos_d, face_atom, bigf, face_area, &
      &               r0_d, rc, qmad, contv_scalingfactor, e_ff_1, e_sff_1, e_fa_1, e_fh_1, e_st_1)
        e_1 = e_ff_1 + e_sff_1 + e_fa_1 + e_fh_1 + e_st_1
        e_dr(id,ia) = (e_1 - e_0)/delta  
        e_ff_dr(id,ia)  = (e_ff_1  - e_ff_0)/delta  
        e_sff_dr(id,ia) = (e_sff_1 - e_sff_0)/delta  
        e_fa_dr(id,ia)  = (e_fa_1 - e_fa_0)/delta  
        e_fh_dr(id,ia)  = (e_fh_1 - e_fh_0)/delta  
        e_st_dr(id,ia)  = (e_st_1 - e_st_0)/delta  

        ! calculate with new bigf but original positions
        call cont_calc_ener(na, nf, ng, face_q, face_pos, face_atom, bigf_d, face_area, &
      &               r0, rc, qmad, contv_scalingfactor, e_ff_1, e_sff_1, e_fa_1, e_fh_1, e_st_1)
        e_1 = e_ff_1 + e_sff_1 + e_fa_1 + e_fh_1 + e_st_1
        e_ds(id,ia) = (e_1 - e_0)/delta  
        e_ff_ds(id,ia)  = (e_ff_1  - e_ff_0)/delta  
        e_sff_ds(id,ia) = (e_sff_1 - e_sff_0)/delta  
        e_fa_ds(id,ia)  = (e_fa_1 - e_fa_0)/delta  
        e_fh_ds(id,ia)  = (e_fh_1 - e_fh_0)/delta  
        e_st_ds(id,ia)  = (e_st_1 - e_st_0)/delta  

     end do
  end do

print *, "+++++++++++++++++++++++++++++++++++++++++++++++++++++"
print *, "          Numeric Derivatives                        "
print *, " "
print *, " total "
print *, e_d
print *, "  face-face  "
print *, e_ff_d
print *, " self face-face  "
print *, e_sff_d
print *, " face-atom      "
print *, e_fa_d 
print *, " hardness  "
print *, e_fh_d
print *, " surface_tension  "
print *, e_st_d
print *, "+++++++++++++++++++++++++++++++++++++++++++++++++++++"
print *, "special contributions                                "
print *, "   face-face  dist only  "
print *, e_ff_dr
print *, "   face-face  surface fct. only"
print *, e_ff_ds
print *, "   face-atom  dist only  "
print *, e_fa_dr
print *, "   face_atom  surfface fct. only"
print *, e_fa_ds
end subroutine cont_num_deriv


! ######################################################################

end module continuum_module







