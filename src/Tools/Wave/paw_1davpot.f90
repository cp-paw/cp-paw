!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM MAIN
!     **************************************************************************
!     **  CALCULATES THE POTENTIAL AVERAGED OVER LATTICE PLANES.              **
!     **  THIS POTENTIAL IS USED TO DETERMINE WORK FUNCTIONS OF BAND OFFSETS. **
!     **  THE INPUT FOR THIS CODE IS THE .WV FILE PRODUCED BY THE SIMULATION  **
!     **  CODE USING !CONTROL!ANALYSE!POTENTIAL                               **
!     **                                                                      **
!     **  THE INPUT ARGUMENT IS THE NAME OF THE .WV FILE.                     **
!     **  THE PROGRAM PRODUCES 3 FILES, ONE FOR EACH DIRECTION.               **
!     **                                                                      **
!     **  NOTE THAT THE DIRECTIONS ARE NOT THE LATTICE VECTORS, BUT THE       **
!     **  LATTICE PLANE NORMALS. THE DISTANCE IS ALSO MEASURED ALONG THE      **
!     **  LATTICE PLANE NORMAL.                                               **
!     **                                                                      **
!     **  EACH POTENTIAL IS CONSTRUCTED OVER TWO UNIT CELLS, SO THAT AN       **
!     **  APPROPRIATE SHEET CAN BE SELECTED                                   **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(256)             :: ISONAME
      INTEGER(4)    ,PARAMETER   :: NFILIN=11
      INTEGER(4)    ,PARAMETER   :: NFILOUT=10
      REAL(8)       ,ALLOCATABLE :: POTENTIAL(:,:,:)  !potential
      INTEGER(4)                 :: NR1,NR2,NR3  ! #(grid points)
      REAL(8)                    :: RBAS(3,3)    ! real space lattice vectors
      CHARACTER(8)               :: DUMMY
      INTEGER(4)                 :: LENTITLE
      CHARACTER(1)  ,ALLOCATABLE :: TITLE(:)
      INTEGER(4)                 :: NATOMS       ! #(atoms)
      REAL(8)       ,ALLOCATABLE :: POS(:,:)     ! atomic positions
      REAL(8)       ,ALLOCATABLE :: Z(:)         ! atomic numbers
      REAL(8)       ,ALLOCATABLE :: Q(:)         ! charges
      CHARACTER(32) ,ALLOCATABLE :: NAME(:)      ! atom names
      REAL(8)                    :: POT
      REAL(8)                    :: DELTA100,DELTA010,DELTA001 ! spacings
      REAL(8)                    :: VEC(3)       ! normal vector
      INTEGER(4)                 :: I
!     **************************************************************************
      CALL MPE$INIT()
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      CALL GET_COMMAND_ARGUMENT(1,ISONAME)
!
!     ==========================================================================
!     ==  READ WV FILE                                                        ==
!     ==========================================================================
      OPEN(NFILIN,FILE=TRIM(ISONAME),STATUS='OLD',FORM='UNFORMATTED')
      WRITE(*,*) 'ISONAME=',ISONAME
      READ(NFILIN) DUMMY,LENTITLE
      WRITE(*,*) 'DUMMY=',DUMMY
      WRITE(*,*) 'LENTITLE=',LENTITLE
      ALLOCATE(TITLE(LENTITLE))
      READ(NFILIN)TITLE(:)
      WRITE(*,*) 'TITLE=',TITLE
      DEALLOCATE(TITLE)
      READ(NFILIN)RBAS,NATOMS
      WRITE(*,*) 'RBAS=',RBAS(:,1)
      WRITE(*,*) RBAS(:,2)
      WRITE(*,*) RBAS(:,3)
      WRITE(*,*) 'NATOMS=',NATOMS
      READ(NFILIN) NR1, NR2, NR3
      WRITE(*,*) 'NR1=',NR1
      WRITE(*,*) 'NR2=',NR2
      WRITE(*,*) 'NR3=',NR3
      ALLOCATE(POTENTIAL(NR1,NR2,NR3))
      POTENTIAL = 0.D0
      WRITE(*,*) 'INITIAL POTENTIAL(1,1,1)=',POTENTIAL(1,1,1)
      WRITE(*,*) 'INITIAL POTENTIAL(NR1,NR2,NR3)=',POTENTIAL(NR1,NR2,NR3)
      ALLOCATE(NAME(NATOMS))
      ALLOCATE(Z(NATOMS))
      ALLOCATE(POS(3,NATOMS))
      ALLOCATE(Q(NATOMS))
      READ(NFILIN)NAME
      WRITE(*,*) 'NAME=',NAME
      READ(NFILIN)Z
      WRITE(*,*) 'Z=',Z
      READ(NFILIN)POS
      WRITE(*,*) 'POS(1)=',POS(:,1)
      WRITE(*,*) 'POS(NATOMS)=',POS(:,NATOMS)
      READ(NFILIN)Q
      WRITE(*,*) 'Q=',Q
      READ(NFILIN)POTENTIAL(:,:,:)
      WRITE(*,*) 'POTENTIAL(1,1,1)=',POTENTIAL(1,1,1)
      WRITE(*,*) 'POTENTIAL(NR1,NR2,NR3)=',POTENTIAL(NR1,NR2,NR3)
      READ(NFILIN)DUMMY
      WRITE(*,*) 'DUMMY=',DUMMY
      IF (DUMMY(1:3).EQ.'END') PRINT*,'READING OF ISO-FILE SUCCESSFUL'
      DEALLOCATE(NAME)
      DEALLOCATE(Z)
      DEALLOCATE(POS)
      DEALLOCATE(Q)
      CLOSE(NFILIN)      
!
!     ==========================================================================
!     ==  PREPARE NORMAL VECTORS AND SPACING                                  ==
!     ==========================================================================
      VEC(1)=RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)
      VEC(2)=RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)
      VEC(3)=RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)
      VEC=VEC/SQRT(SUM(VEC**2))
      DELTA100=DOT_PRODUCT(RBAS(:,1),VEC)/REAL(NR1,KIND=8)

      VEC(1)=RBAS(2,3)*RBAS(3,1)-RBAS(3,3)*RBAS(2,1)
      VEC(2)=RBAS(3,3)*RBAS(1,1)-RBAS(1,3)*RBAS(3,1)
      VEC(3)=RBAS(1,3)*RBAS(2,1)-RBAS(2,3)*RBAS(1,1)
      VEC=VEC/SQRT(SUM(VEC**2))
      DELTA010=DOT_PRODUCT(RBAS(:,2),VEC)/REAL(NR2,KIND=8)

      VEC(1)=RBAS(2,1)*RBAS(3,2)-RBAS(3,1)*RBAS(2,2)
      VEC(2)=RBAS(3,1)*RBAS(1,2)-RBAS(1,1)*RBAS(3,2)
      VEC(3)=RBAS(1,1)*RBAS(2,2)-RBAS(2,1)*RBAS(1,2)
      VEC=VEC/SQRT(SUM(VEC**2))
      DELTA001=DOT_PRODUCT(RBAS(:,3),VEC)/REAL(NR3,KIND=8)
!
!     ==========================================================================
!     ==  WRITE POTENTIAL AVERAGED OVER (100) PLANES                          ==
!     ==========================================================================
      OPEN(NFILOUT,FILE="AVPOT100.DAT")      
      DO I=1,NR1
         POT=SUM(POTENTIAL(I,:,:))/REAL(NR2*NR3,KIND=8)
         WRITE(NFILOUT,*) DELTA100*REAL(I-1,KIND=8),POT
      END DO
      DO I=1,NR1
         POT=SUM(POTENTIAL(I,:,:))/REAL(NR2*NR3,KIND=8)
         WRITE(NFILOUT,*) DELTA100*REAL(NR1+I-1,KIND=8),POT
      END DO
      CLOSE(NFILOUT)      
!
!     ==========================================================================
!     ==  WRITE POTENTIAL AVERAGED OVER (010) PLANES                          ==
!     ==========================================================================
      OPEN(NFILOUT,FILE="AVPOT010.DAT")      
      DO I=1,NR2
         POT=SUM(POTENTIAL(:,I,:))/REAL(NR1*NR3,KIND=8)
         WRITE(NFILOUT,*) DELTA010*REAL(I-1,KIND=8),POT
      END DO
      DO I=1,NR2
         POT=SUM(POTENTIAL(:,I,:))/REAL(NR1*NR3,KIND=8)
         WRITE(NFILOUT,*) DELTA010*REAL(NR2+I-1,KIND=8),POT
      END DO
      CLOSE(NFILOUT)      
!
!     ==========================================================================
!     ==  WRITE POTENTIAL AVERAGED OVER (001) PLANES                          ==
!     ==========================================================================
      OPEN(NFILOUT,FILE="AVPOT001.DAT")      
      DO I=1,NR3
         POT=SUM(POTENTIAL(:,:,I))/REAL(NR1*NR2, KIND=8)
         WRITE(NFILOUT,*) DELTA001*REAL(I-1,KIND=8),POT
      END DO
      DO I=1,NR3
         POT=SUM(POTENTIAL(:,:,I))/REAL(NR1*NR2, KIND=8)
         WRITE(NFILOUT,*) DELTA001*REAL(NR3+I-1,KIND=8),POT
      END DO
      CLOSE(NFILOUT)      
!
      CALL ERROR$NORMALSTOP()
      STOP
      END PROGRAM MAIN

