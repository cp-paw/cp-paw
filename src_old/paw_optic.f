!
!..........................................................OPTIC........
MODULE OPTIC_MODULE
LOGICAL(4) :: ON=.FALSE.    !TURNS SUBROUTINE OPTIC OFF AND ON
LOGICAL(4) :: TINI=.FALSE.  
INTEGER(4) :: DHIAT(2000)  
INTEGER(4) :: DOIAT(2000) 
END MODULE OPTIC_MODULE
!
!     ..................................................................
      SUBROUTINE OPTIC3$WAVES(NGWAVES,NBANDS,NKPOINTS,NSPIN,C0)
      USE OPTIC_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NGWAVES
      INTEGER(4),INTENT(IN) :: NSPIN
      INTEGER(4),INTENT(IN) :: NKPOINTS   
      INTEGER(4),INTENT(IN) :: NBANDS   
      COMPLEX(8),INTENT(IN) :: C0(NGWAVES,NBANDS,NKPOINTS,NSPIN)  
      INTEGER(4)            :: IG   
      INTEGER(4)            :: IBAND  
      INTEGER(4)            :: IKPT   
      INTEGER(4)            :: NFIL
!     ******************************************************************
      IF(.NOT.ON) RETURN
      if(.not.tini)call OPTIC3_INIT
!     CALL FILEHANDLER$SETFILE(+'OPTICS_WAVES',.TRUE.-'.OPTICS_WAVES')
!     WRITE(NFIL,*)NGWAVES,NBANDS,NKPOINTS,NSPIN
!     DO IKPT=1,NKPOINTS  
!       DO IBAND=1,NBANDS  
!         DO IG=1,NGWAVES   
!           WRITE(NFIL,*)(C0(IG,IBAND,IKPT,ISPIN),ISPIN=1,NSPIN)  
!         ENDDO 
!       ENDDO  
!     ENDDO  
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE OPTIC3$VOFG(NR1,NR2,NR3,NVOFG,RHOE,NSPINV,INDV1,INDV2,INDV3)
      USE OPTIC_MODULE
      IMPLICIT NONE
      INTEGER(4) :: NR1   !DIMENSION OF FFT IN X-DIRECTION  
      INTEGER(4) :: NR2   !DIMENSION OF FFT IN Y-DIRECTION
      INTEGER(4) :: NR3   !DIMENSION OF FFT IN Z-DIRECTION
      INTEGER(4) :: NVOFG !NUMBER OF PLANE WAVES IN REPRESENTATION OF 
                    !THE THE PLANE WAVE POTENTIAL V~.  
                    !COMPLEX CONJUGATES ARE NOT INCLUDED.  
      REAL(8) :: RHOE(NR1*NR2*NR3,NSPINV) !REAL SPACE REPRESENTATION OF THE 
                                      !ELECTRON DENSITY DUE TO THE
                                      !PLANE WAVE PART OF THE WAVE
                                      !FUNCTIONS FOR SPIN UP/DOWN.      
      INTEGER(4) :: NSPINV !SPIN INDEX USED TO DIMENSION ARRAYS ASSOCIATED
                     !WITH V~      
      INTEGER(4) :: INDV1(NVOFG) !INTEGER(4) :: NUMBER OF DISPLACEMENTS 
                           !IN X-DIRECTION FOR V~ 
      INTEGER(4) :: INDV2(NVOFG) !INTEGER(4) :: NUMBER OF DISPLACEMENTS 
                           !IN Y-DIRECTION FOR V~ 
      INTEGER(4) :: INDV3(NVOFG) !INTEGER(4) :: NUMBER OF DISPLACEMENTS 
                           !IN Z-DIRECTION FOR V~ 
      REAL(8) :: DUMMY1(NR1*NR2*NR3) !DUMMY ARRAY USED FOR THE ELECTRON DENSITY
                                 !OF SPIN 1
      REAL(8) :: DUMMY2(NR1*NR2*NR3) !DUMMY ARRAY USED FOR THE ELECTRON DENSITY
                                 !OF SPIN 2
      INTEGER(4) :: NFFT  !NUMBER OF FFT'S TO BE PERFORMED.
                    !THE PROGRAM USES A GENERAL COMPLEX FFT
                    !ROUTINE WHICH CAN BE USED TO 
                    !  (1) PERFORM 1 STANDARD FFT ON 1 COMPLEX ARRAY
                    !  (2) PERFORM 2 STANDARD FFT'S ON 2 REAL ARRAYS 
                    !  (3) PERFORM 1 STANDARD FFT ON 1 REAL ARRAY
                    !  WHEN DOING (3) => NFFT=1   
      INTEGER(4) :: NVOFR !NUMBER OF POINTS IN REAL SPACE REPRESENTATION
                    !OF THE PLANE WAVE POTENTIAL V~. 
                    !NVOFR=NR1*NR2*NR3
      COMPLEX(8) :: VOFG(NVOFG)         !FOURIER COEFFICIENTS OF V~
      COMPLEX(8) :: VOFG1               !VARIABLE NEEDED FOR FFT BUT
      INTEGER(4) :: I
      INTEGER(4) :: NFIL
!     ******************************************************************
      IF(.NOT.ON)RETURN  
      if(.not.tini)call OPTIC3_INIT
      DO I=1,NR1*NR2*NR3
        DUMMY1(I)=RHOE(I,1)
        DUMMY2(I)=RHOE(I,NSPINV)
      ENDDO   
      NFFT=1
      NVOFR=NR1*NR2*NR3   
      CALL PLANEWAVE$SELECT('DENSITY',1)
      CALL PLANEWAVE$SUPFFT('RTOG',' ',NFFT,NVOFR,NR1,NR2,NR3,NVOFG &
     &                   ,DUMMY1,DUMMY2,VOFG,VOFG1)
      CALL FILEHANDLER$UNIT('OPTICS_VOFG',NFIL)
      WRITE(NFIL,*)NVOFG
      DO I=1,NVOFG  
        WRITE(NFIL,*)INDV1(I),INDV2(I),INDV3(I),VOFG(I)
      ENDDO  
      CALL FILEHANDLER$CLOSE('OPTICS_VOFG')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE OPTIC3$ATOMS(NAT,LNXX,NSP,NRX,R1,DEX &
     &                  ,LOX,LNX,LMNX,PRO,AEPHI,PSPHI)
!     ******************************************************************
!     ******************************************************************
      USE OPTIC_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT   !NUMBER OF ATOMS
      INTEGER(4),INTENT(IN) :: LNXX  ! MAX NUMBER OF PROJECTORS (NOT COUNTING M 
                    ! MULTIPLICITY) OVER ALL SPECIES.
                    ! EXAMPLE: IF A CALCULATION IS PERFORMED
                    ! USING THE FOLLOWING SETUPS...
                    ! 
                    ! SPECIES 1: GA
                    !     QUANTUM NUMBER L    NUMBER OF PROJECTORS
                    !     ----------------    --------------------
                    !            0                      2
                    !            1                      2
                    !            2                      2
                    ! 
                    ! SPECIES 2: AS
                    !     QUANTUM NUMBER L    NUMBER OF PROJECTORS
                    !     ----------------    --------------------
                    !            0                      2
                    !            1                      2
                    !            2                      1
                    ! 
                    ! IN THIS CASE, LNXX=6.  HOWEVER, IF THE CALCULATION
                    ! INVOLVED AS ONLY, LNXX WOULD HAVE BEEN 5.  
      INTEGER(4),INTENT(IN) :: NSP   !NUMBER OF SPECIES
      INTEGER(4),INTENT(IN) :: NRX   !NUMBER OF POINTS ON RADIAL GRID  
      REAL(8)   ,INTENT(IN) :: R1     !FIRST POINT ON RADIAL GRID  
      REAL(8)   ,INTENT(IN) :: DEX    !DELTA X ON RADIAL GRID      
      INTEGER(4),INTENT(IN) :: LOX(LNXX,NSP)          !QUANTUM NUMBER L FOR EACH
                                     !PROJECTOR AND EACH SPECIES.
                                     !FOR THE EXAMPLE ABOVE 
                                     !(SEE EXPLANATION OF LNXX)
                                     ! LOX(1..6,1)=(0,0,1,1,2,2)
                                     ! LOX(1..6,2)=(0,0,1,1,2,0)
      INTEGER(4),INTENT(IN) :: LNX(NSP)               !NUMBER OF PROJECTORS FOR EACH
                                     !SPECIES.  FOR THE EXAMPLE ABOVE
                                     !(SEE EXPLANATION OF LNXX)  
                                     !     LNX(1)=6
                                     !     LNX(2)=5   
      INTEGER(4),INTENT(IN) :: LMNX(NSP)              ! IF N IS THE NUMBER OF PROJECTORS
                                     ! AT L, THEN LMNX(NSP) IS THE 
                                     ! TOTAL NUMBER OF LMN COMBOS
                                     ! FOR EACH SPECIES:
                                     !
                                     !   L=L_MAX
                                     !    -----
                                     !     \
                                     !      )   (2 L + 1)(# PROS FOR L)
                                     !     /
                                     !    -----
                                     !     L=0
      REAL(8) ,INTENT(IN) :: PRO(NRX,LNXX,NSP)!PROJECTOR FUNCTIONS ON RADIAL
                                    !LOGARITHMIC GRID.  
      REAL(8) ,INTENT(IN) :: AEPHI(NRX,LNXX,NSP) !ALL-ELECTRON 
                                       !ATOMIC WAVE FUNCTIONS 
                                       !ON RADIAL LOGARITHMIC GRID.  
      REAL(8) ,INTENT(IN) :: PSPHI(NRX,LNXX,NSP) !PSEUDO-ATOMIC WAVE FUNCTIONS ON 
                                       !RADIAL LOGARITHMIC GRID.  


      CHARACTER(32)       ::  NAME !VARIABLE USE TO TEMPORARILY STORE NAMES OF ATOMS 
      REAL(8)             :: RPOS(3)       ! DUMMY ATOMIC POSITION VECTOR   
      INTEGER(4)          :: ISP   !SPECIES INDEX    
      INTEGER(4)          :: GRIDPT !RADIAL GRID POINT INDEX      
      INTEGER(4)          :: LN    !DUMMY LN INDEX      
      INTEGER(4)          :: IXYZ  !ATOM INDEX    
      INTEGER(4)            :: IAT   !ATOM INDEX  
      INTEGER(4)            :: NFIL  !ATOM INDEX  
!     ******************************************************************
      IF(.NOT.ON)RETURN
      if(.not.tini)call OPTIC3_INIT
      CALL FILEHANDLER$UNIT('OPTICS_ATOMS',NFIL)
      WRITE(NFIL,*)'ATOMS   '  
      WRITE(NFIL,*)NAT  
      DO IAT=1,NAT
        CALL ATOMLIST$GETCH('NAME',IAT,NAME)
        CALL ATOMLIST$GETR8A('R(0)',IAT,3,RPOS)
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        WRITE(NFIL,*)IAT  
        WRITE(NFIL,*)ISP
        WRITE(NFIL,*)NAME  
        WRITE(NFIL,*)RPOS  
        WRITE(NFIL,*)LNX(ISP)
        WRITE(NFIL,*)LMNX(ISP)
        WRITE(NFIL,*)(LOX(LN,ISP),LN=1,LNX(ISP))
        WRITE(NFIL,*)NRX,R1,DEX
        DO GRIDPT=1,NRX
            WRITE(NFIL,*)(PRO(GRIDPT,LN,ISP),LN=1,LNX(ISP))  
        ENDDO  
        DO GRIDPT=1,NRX
            WRITE(NFIL,*)(AEPHI(GRIDPT,LN,ISP),LN=1,LNX(ISP))  
        ENDDO  
        DO GRIDPT=1,NRX
            WRITE(NFIL,*)(PSPHI(GRIDPT,LN,ISP),LN=1,LNX(ISP))  
        ENDDO  
      ENDDO  
      RETURN   
      END
!
!     ..................................................................
      SUBROUTINE OPTIC3$LATTICE(RBAS,GBAS,OMEGA &
     &                    ,NG,NGMAX,NPAWKPT,D_INDG1,D_INDG2,D_INDG3)
      USE OPTIC_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)     !LATTICE BASIS VECTORS.  
      REAL(8)   ,INTENT(IN) :: GBAS(3,3)     !RECIPROCAL LATTICE BASIS VECTORS.  
      REAL(8)   ,INTENT(IN) :: OMEGA         !UNIT CELL VOLUME   
      INTEGER(4),INTENT(IN) :: NG    !NUMBER OF PLANE WAVES IN WAVE FUNCTION
                                     !(ALSO THE ORDER OF THE HAMILTONIAN MATRIX)  
      INTEGER(4),INTENT(IN) :: NGMAX !MAXIMUM NUMBER OF PLANE WAVES IN WAVE FUNCTION
      INTEGER(4),INTENT(IN) :: NPAWKPT !NUMBER OF K-POINTS USED IN PAW CODE  
      INTEGER(4),INTENT(IN) :: D_INDG1(NGMAX,NPAWKPT) ! NUMBER OF DISPLACEMENTS 
                        !IN X-DIRECTION FOR PLANE WAVE
                        !PART OF WAVE FUNCTION FOR EACK K-POINT.   
      INTEGER(4),INTENT(IN) :: D_INDG2(NGMAX,NPAWKPT) ! NUMBER OF DISPLACEMENTS 
                        !IN Y-DIRECTION FOR PLANE WAVE
                        !PART OF WAVE FUNCTION FOR EACK K-POINT. 
      INTEGER(4),INTENT(IN) :: D_INDG3(NGMAX,NPAWKPT) ! NUMBER OF DISPLACEMENTS 
                        !IN Z-DIRECTION FOR PLANE WAVE
                        !PART OF WAVE FUNCTION FOR EACK K-POINT. 
      INTEGER(4)            :: INDG1(NG) !INTEGER(4) :: NUMBER OF DISPLACEMENTS 
                        !IN X-DIRECTION FOR PLANE WAVE
                        !PART OF WAVE FUNCTION.   
      INTEGER(4)            :: INDG2(NG) !INTEGER(4) :: NUMBER OF DISPLACEMENTS 
                        !IN Y-DIRECTION FOR PLANE WAVE
                        !PART OF WAVE FUNCTION.   
      INTEGER(4)            :: INDG3(NG) !INTEGER(4) :: NUMBER OF DISPLACEMENTS 
                        !IN Z-DIRECTION FOR PLANE WAVE
                        !PART OF WAVE FUNCTION.   
      INTEGER(4)            :: I
      INTEGER(4)            :: NFIL
!     ******************************************************************
      IF(.NOT.ON)RETURN
      if(.not.tini)call OPTIC3_INIT
      DO I=1,NG
        INDG1(I)=D_INDG1(I,1)
        INDG2(I)=D_INDG2(I,1)
        INDG3(I)=D_INDG3(I,1)
      ENDDO  
      CALL FILEHANDLER$UNIT('OPTICS_LATTICE',NFIL)

      WRITE(NFIL,*)RBAS
      WRITE(NFIL,*)GBAS
      WRITE(NFIL,*)OMEGA
      WRITE(NFIL,*)NG
      DO I=1,NG    
        WRITE(NFIL,*)INDG1(I),INDG2(I),INDG3(I)
      ENDDO  
      CALL FILEHANDLER$CLOSE('OPTICS_LATTICE')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE OPTIC3$DH(NSPINDH,LMNXXDH,NATDH,IATDH,DH)
!     ******************************************************************
!     **  DH(I,J)=DH_IJ EQN. 99 IN PAW PAPER) 
!     **  DH(I,J)=<PHI_I| H_AT  |PHI_J>  - <PHI_I~| H_AT~ |PHI_J~>
!     **  WHERE H_AT IS THE RADIAL ATOMIC HAMILTONIAN AND 
!     **  H_AT~ IS THE RADIAL PSEUDO-ATOMIC HAMILTONIAN.  
!     **  PHI_I =ATOMIC WAVE FUNCTION
!     ** PHI_1~=PSEUDO-ATOMIC WAVE FUNCTION
!     ** I,J ARE INDICES FOR THE (LMN)A COMBINATIONS.  FOR THE EXAMPLE ABOVE...

         ! SPECIES 1: GA
         ! QUANTUM NUMBER L  QUANTUM NUMBER M  PROJECTOR INDEX   I  
         ! ----------------  ----------------  ---------------  ---
         !        0                0                 1           1
         !        0                0                 2           2
         !        1               -1                 1           3
         !        1                0                 1           4
         !        1               +1                 1           5
         !        1               -1                 2           6
         !        1                0                 2           7
         !        1               +1                 2           8
         !        2               -2                 1           9
         !        2               -1                 1          10
         !        2                0                 1          11
         !        2               +1                 1          12
         !        2               +2                 1          13
         !        2               -2                 2          14
         !        2               -1                 2          15
         !        2                0                 2          16
         !        2               +1                 2          17
         !        2               +2                 2          18
         ! 
         ! SPECIES 2: AS
         ! QUANTUM NUMBER L  QUANTUM NUMBER M  PROJECTOR INDEX   I
         ! ----------------  ----------------  ---------------  ---
         !        0                0                 1           1
         !        0                0                 2           2
         !        1               -1                 1           3
         !        1                0                 1           4
         !        1               +1                 1           5
         !        1               -1                 2           6
         !        1                0                 2           7
         !        1               +1                 2           8
         !        2               -2                 1           9
         !        2               -1                 1          10
         !        2                0                 1          11
         !        2               +1                 1          12
         !        2               +2                 1          13

!     ******************************************************************
      USE OPTIC_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSPINDH!NUMBER OF SPINS 
      INTEGER(4),INTENT(IN) :: NATDH  !NUMBER OF ATOMS   
      INTEGER(4),INTENT(IN) :: IATDH   !ATOM INDEX USED ONLY IN $DH    
      REAL(8)   ,INTENT(IN) :: DH(LMNXXDH,LMNXXDH,NSPINDH)   
      INTEGER(4) :: ISUM  
      INTEGER(4) :: IAT2  !ATOM INDEX    
      INTEGER(4) :: LMNXXDH                ! MAX VALUE OF LMNX.  IN EXAMPLE 
                                     ! ABOVE LMNXX=18.  
                                     ! USED ONLY IN $DH
      INTEGER(4) :: LMNXXDH2                ! MAX VALUE OF LMNX.  IN EXAMPLE 
      INTEGER(4) :: ISPIN !DUMMY SPIN INDEX      
      INTEGER(4) :: LMN1  !DUMMY LMN INDEX      
      INTEGER(4) :: NFIL
      INTEGER(4) :: LMN2  !DUMMY LMN INDEX      
!     ******************************************************************
      IF(.NOT.ON)RETURN
      if(.not.tini)call OPTIC3_INIT
      ISUM=0
      DO IAT2=1,NATDH
       ISUM=ISUM+DHIAT(IAT2)
      ENDDO 
      IF(ISUM.EQ.NATDH)RETURN  
      CALL FILEHANDLER$UNIT('OPTICS_ATOMS',NFIL)
      IF(DHIAT(IATDH).EQ.0)THEN  
        DHIAT(IATDH)=1
        WRITE(NFIL,*)'DH      '
        WRITE(NFIL,*)NATDH
        WRITE(NFIL,*)IATDH
        WRITE(NFIL,*)NSPINDH,LMNXXDH
          DO ISPIN=1,NSPINDH
            DO LMN1=1,LMNXXDH
              DO LMN2=1,LMNXXDH
                WRITE(NFIL,*)LMN1,LMN2,ISPIN,IATDH &
     &                     ,DH(LMN1,LMN2,ISPIN)
              ENDDO 
            ENDDO 
          ENDDO 
      ENDIF  
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE OPTIC3$DO(NSPINDO,LMNXXDO,NATDO,IATDO,DO)
! DO(I,J)=DO_IJ (EQN. 100 IN PAW PAPER) DH(I,J)=<PHI_I|PHI_J>-<PHI_I~|PHI_J~>
! PHI_I =ATOMIC WAVE FUNCTION;  PHI_1~=PSEUDO-ATOMIC WAVE FUNCTION
      USE OPTIC_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSPINDO ! NUMBER OF SPINS (NSPIN=1 FOR NO SPIN POLARIZATION)
      INTEGER(4),INTENT(IN) :: LMNXXDO ! MAX VALUE OF LMNX.  IN EXAMPLE ABOVE LMNXX=18. 
      INTEGER(4),INTENT(IN) :: NATDO   !NUMBER OF ATOMS   
      INTEGER(4),INTENT(IN) :: IATDO   !ATOM INDEX USED ONLY IN $DO
      REAL(8)   ,INTENT(IN) :: DO(LMNXXDO,LMNXXDO,NSPINDO)
      INTEGER(4)            :: LMNXXDO2 ! MAX VALUE OF LMNX.  IN EXAMPLE 
      INTEGER(4)            :: ISUM,IAT2
      INTEGER(4)            :: NFIL
      INTEGER(4)            :: ISPIN,LMN1,LMN2
!     ******************************************************************
      IF(.NOT.ON)RETURN
      if(.not.tini)call OPTIC3_INIT
      ISUM=0
      DO IAT2=1,NATDO
       ISUM=ISUM+DOIAT(IAT2)
      ENDDO 
      IF(ISUM.EQ.NATDO) RETURN  
      CALL FILEHANDLER$UNIT('OPTICS_ATOMS',NFIL)
      IF(DOIAT(IATDO).EQ.0)THEN  
        DOIAT(IATDO)=1
        WRITE(NFIL,*)'DO      '
        WRITE(NFIL,*)NATDO
        WRITE(NFIL,*)IATDO
        WRITE(NFIL,*)NSPINDO,LMNXXDO
          DO ISPIN=1,NSPINDO
            DO LMN1=1,LMNXXDO
              DO LMN2=1,LMNXXDO
                WRITE(NFIL,*)LMN1,LMN2,ISPIN,IATDO,DO(LMN1,LMN2,ISPIN)
              ENDDO 
            ENDDO 
          ENDDO 
      ENDIF  
      RETURN  
      END   
!
!     ..................................................................
      SUBROUTINE OPTIC3$SETL4(ID,VAL)
      USE OPTIC_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'ON') THEN
        ON=VAL
      ELSE
        CALL ERROR$MSG('UNKNOWN ID')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('OPTIC3$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE OPTIC3_INIT
      USE OPTIC_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
!     ******************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
      DHIAT(:)=0  
      DOIAT(:)=0  
!
!     ==================================================================
!     ==  DEFINE FILES                                                ==
!     ==================================================================
      CALL FILEHANDLER$SETFILE(+'OPTICS_VOFG',.TRUE.,-'.OPTICS_VOFG')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_VOFG','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_VOFG','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_VOFG','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_VOFG','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE(+'OPTICS_ATOMS',.TRUE.,-'.OPTICS_ATOMS')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_ATOMS','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_ATOMS','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_ATOMS','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_ATOMS','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE(+'OPTICS_LATTICE',.TRUE.,-'.OPTICS_LATTICE')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_LATTICE','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_LATTICE','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_LATTICE','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_LATTICE','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE(+'OPTICS_WAVES',.TRUE.,-'.OPTICS_WAVES')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_WAVES','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_WAVES','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_WAVES','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(+'OPTICS_WAVES','FORM','FORMATTED')
      RETURN 
      END
