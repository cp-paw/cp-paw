      PROGRAM MAIN
!     ******************************************************************
!     ** use as                                                       **
!     **  converttra infile outfile deltat                            **
!     ******************************************************************
      implicit none
      REAL(8),ALLOCATABLE :: R(:)
      REAL(8),ALLOCATABLE :: q(:)
      REAL(8)             :: DELTAT
      INTEGER(4)          :: ISTEP
      INTEGER(4)          :: NAT
      INTEGER(4)          :: NFILIN=1000
      INTEGER(4)          :: NFILOUT=1001
      REAL(8)             :: T
      CHARACTER(256)      :: FILEIN
      CHARACTER(256)      :: FILEOUT
      CHARACTER(256)      :: TEXT
!     ******************************************************************      
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT

      CALL GET_COMMAND_ARGUMENT(1,FILEIN)
      CALL GET_COMMAND_ARGUMENT(2,FILEOUT)
      CALL GET_COMMAND_ARGUMENT(3,TEXT)
      READ(TEXT,*)DELTAT
print*,'filein ',trim(filein)
print*,'fileout ',trim(fileout)
print*,'deltat ',deltat
      OPEN(NFILIN,FILE=FILeiN,STATUS='OLD',FORM='UNFORMATTED',ACTION='READ')
      OPEN(NFILOUT,FILE=FILeOUT,STATUS='NEW',FORM='UNFORMATTED',ACTION='WRITE')
      REWIND NFILIN
      READ(NFILIN)ISTEP,NAT
print*,'nat ',nat
      ALLOCATE(R(3*NAT))
      ALLOCATE(Q(NAT))
      REWIND NFILIN
print*,'convertingg....'
      DO
        READ(NFILIN,END=1000)ISTEP,NAT,R
        T=DELTAT*REAL(ISTEP-1)
        WRITE(NFILOUT)ISTEP,T,4*NAT,R,Q
      ENDDO
 1000 CONTINUE
print*,'conversion finished'
      DEALLOCATE(R)
      DEALLOCATE(Q)
      CALL ERROR$NORMALSTOP()
      STOP
      END
