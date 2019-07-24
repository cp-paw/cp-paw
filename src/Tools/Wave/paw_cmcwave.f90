
!     ...1.........2.........3.........4.........5.........6.........7.........8
      Program main
!     **************************************************************************
!     ** the program takes a .wv file (wave function,density or potential)    **
!     ** from the CP-PAW  simulation code and produces a .cmcv input file     **
!     ** for crymolcad                                                        **
!     **************************************************************************
      implicit none
      integer(4)   ,parameter   :: nfilin=10
      integer(4)   ,parameter   :: nfilout=11
      character(128)            :: title
      real(8)                   :: rbas(3,3)
      integer(4)                :: nat=0
      real(8)      ,allocatable :: pos(:,:)
      real(8)      ,allocatable :: z(:)
      character(32),allocatable :: name(:)
      integer(4)                :: nr1
      integer(4)                :: nr2
      integer(4)                :: nr3
      real(8)      ,allocatable :: wave(:,:,:)
      character(256)            :: infile
      character(256)            :: outfile
      integer(4)                :: i
!     **************************************************************************
      call GET_COMMAND_ARGUMENT(1,infile)
      i=index(infile,'.',.true.)
      outfile=infile(1:i)//'cmcv'
!
!     ==========================================================================
!     == obtain dimensions                                                    ==
!     ==========================================================================
print*,'reading from file ',infile,'....'
      open(nfilin,file=infile,form='unformatted')
      nat=0
      allocate(pos(3,1))
      allocate(z(1))
      allocate(name(1))
      allocate(wave(1,1,1))
      call READWAVEPLOT(NFILin,TITLE,RBAS,NAT,POS,Z,NAME,NR1,NR2,NR3,WAVE)
      deallocate(pos)
      deallocate(z)
      deallocate(name)
      deallocate(wave)
!
!     ==========================================================================
!     == read wv file                                                         ==
!     ==========================================================================
      allocate(pos(3,nat))
      allocate(z(nat))
      allocate(name(nat))
      allocate(wave(nr1,nr2,nr3))
      call READWAVEPLOT(NFILin,TITLE,RBAS,NAT,POS,Z,NAME,NR1,NR2,NR3,WAVE)
      close(nfilin)
!
!     ==========================================================================
!     == write cmcv file                                                      ==
!     ==========================================================================
print*,'wrting to file ',outfile,' ....'
      open(nfilout,file=outfile,form='formatted')
      call writecmcv(nfilout,title,rbas,nr1,nr2,nr3,wave)
      close(nfilout)
      stop
      end

!
!     ..................................................................
      SUBROUTINE READWAVEPLOT(NFIL,TITLE,RBAS,NAT,POS,Z,NAME,NR1,NR2,NR3,WAVE)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)   :: NFIL
      CHARACTER(*) ,INTENT(OUT ) :: TITLE
      INTEGER(4)   ,INTENT(INOUT):: NAT
      REAL(8)      ,INTENT(OUT)  :: RBAS(3,3)
      REAL(8)      ,INTENT(OUT)  :: POS(3,NAT)
      REAL(8)      ,INTENT(OUT)  :: Z(NAT)
      REAL(8)                    :: Q(NAT) !POINT CHARGES (NOT YET USED)
      CHARACTER(32),INTENT(OUT)  :: NAME(NAT)
      INTEGER(4)   ,INTENT(INOUT):: NR1
      INTEGER(4)   ,INTENT(INOUT):: NR2
      INTEGER(4)   ,INTENT(INOUT):: NR3
      REAL(8)      ,INTENT(OUT)  :: WAVE(NR1,NR2,NR3)
      CHARACTER(8)               :: BEGINID
      INTEGER(4)                 :: LENTITLE
      CHARACTER(11)              :: ENDID
!     ******************************************************************
      REWIND NFIL
      READ(NFIL)BEGINID,LENTITLE
      LENTITLE=MIN(LENTITLE,LEN(TITLE))
      READ(NFIL)TITLE(1:LENTITLE)
      READ(NFIL)RBAS,NAT
      READ(NFIL)NR1,NR2,NR3
      IF(NAT.EQ.0) THEN
        RETURN
      END IF
      READ(NFIL)NAME
      READ(NFIL)Z
      READ(NFIL)POS
      READ(NFIL)Q
      READ(NFIL)WAVE
      READ(NFIL)ENDID 
      IF(ENDID.NE.'END OF FILE') THEN
      END IF
      RETURN
      END
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine writecmcv(nfil,title,rbas,n1,n2,n3,field)
      implicit none
      integer(4)  ,intent(in) :: nfil
      character(*),intent(in) :: title
      real(8)     ,intent(in) :: rbas(3,3)
      integer(4)  ,intent(in) :: n1,n2,n3
      real(8)     ,intent(in) :: field(n1,n2,n3)
      integer(4)              :: i,j,k
!     **************************************************************************
      rewind(nfil)
      write(nfil,*)title
      write(nfil,*)n1,n2,n3
      write(nfil,*)rbas
      do k=1,n3
        do j=1,n2
          do i=1,n1
            write(nfil,*)field(i,j,k)
          enddo
        enddo
      enddo
      return
      end
