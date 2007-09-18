!     ..................................................................
      PROGRAM PAW_STRC
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT none
      TYPE(LL_TYPE)               :: LL_STRC
      INTEGER(4)                  :: NFIL
      LOGICAL                     :: TCRYSTAL=.FALSE.
      REAL(8)                     :: ANGSTROM
      CHARACTER(128)              :: ROOTNAME     ! COMMON ROOT OF THE FILENAMES
      CHARACTER(128)              :: OBJECTNAME
      LOGICAL                     :: TCHK
      REAL(8)                     :: RUNIT        ! LENGTH UNIT ON STRUCTURE FILE
      REAL(8)                     :: RBAS(3,3)    ! LATTICE VECTORS
      CHARACTER(32),ALLOCATABLE   :: NAME(:)      ! ATOM NAMES
      integer(4)                  :: nat          ! #(atoms in the qm part)
      REAL(8),      ALLOCATABLE   :: R(:,:)       ! ATOMIC POSITIONS
      REAL(8),      ALLOCATABLE   :: Q(:)         ! CHARGES
      integer(4)                  :: natmm        ! #(atoms of the MM part)
      CHARACTER(32),ALLOCATABLE   :: MMNAME(:)    ! ATOM NAMES
      REAL(8),      ALLOCATABLE   :: MMR(:,:)     ! ATOMIC POSITIONS FOR MM PART OF QM-MM
      REAL(8),      ALLOCATABLE   :: MMQ(:)       ! ATOMIC CHARGES MM PART OF QM-MM
      REAL(8),      ALLOCATABLE   :: RSH(:,:)     ! ATOMIC POSITIONS FOR SHADOW PART OF QM-MM
      INTEGER(4),   ALLOCATABLE   :: NEIGH(:,:)   ! NEIGHBOR LIST
      LOGICAL(4)                  :: TQMMM
      INTEGER(4)                  :: IAT,IAT1
      CHARACTER(32)               :: STRING
      integer(4)                  :: i    
!     ******************************************************************
      CALL LINKEDLIST$NEW(LL_STRC)
!
!     ==================================================================
!     == GET FILE NAME ROOT FROM THE ARGUMENT LIST AND CONSTRUCT      ==
!     == FILE NAMES                                                   ==
!     ==================================================================
      CALL GETARG(1,ROOTNAME)
      IF(ROOTNAME(1:1).EQ.'-') THEN
        TCRYSTAL=(+ROOTNAME(2:2).EQ.+'C')
        CALL GETARG(2,ROOTNAME)
      END IF
      IF(LEN(TRIM(ROOTNAME)).EQ.0) THEN
        STOP 'NO ROOTNAME SUPPLIED'
      END IF
      I=INDEX(ROOTNAME,'/',BACK=.TRUE.)
      OBJECTNAME=ROOTNAME(I+1:)
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL FILEHANDLER$SETFILE('STRC',.TRUE.,-'.STRC_OUT')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('CSSR',.TRUE.,-'.CSSR')
      CALL FILEHANDLER$SETSPECIFICATION('CSSR','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('CSSR','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('CSSR','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('CSSR','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('MMCSSR',.TRUE.,-'_MM.CSSR')
      CALL FILEHANDLER$SETSPECIFICATION('MMCSSR','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('MMCSSR','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('MMCSSR','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('MMCSSR','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('SHCSSR',.TRUE.,-'_SH.CSSR')
      CALL FILEHANDLER$SETSPECIFICATION('SHCSSR','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('SHCSSR','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('SHCSSR','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('SHCSSR','FORM','FORMATTED')
!
!     ==================================================================
!     == READ STRUCTURE FILE                                          ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('STRC',NFIL)
      CALL LINKEDLIST$READ(LL_STRC,NFIL,'~')
!
!     ==================================================================
!     == GET LATTICE VECTORS                                          ==
!     ==================================================================
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'GENERIC')
      CALL LINKEDLIST$GET(LL_STRC,'LUNIT',1,RUNIT)
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'LATTICE')
      CALL LINKEDLIST$GET(LL_STRC,'T',1,RBAS)
      RBAS=RBAS*RUNIT
!
!     ==================================================================
!     ==  READ ATOM DATA FROM STRC FILE                               ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NAT)
      ALLOCATE(NAME(NAT))
      ALLOCATE(R(3,NAT))
      ALLOCATE(Q(NAT)) ; Q(:)=0.D0
      ALLOCATE(NEIGH(8,NAT))
      NEIGH=0
      DO IAT=1,NAT
        CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IAT)
        CALL LINKEDLIST$GET(LL_STRC,'R',1,R(:,IAT))
        R(:,IAT)=R(:,IAT)*RUNIT
        CALL LINKEDLIST$EXISTD(LL_STRC,'Q',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'Q',1,Q(IAT))
        ELSE
          Q(IAT)=0.D0
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'NAME',1,NAME(IAT))
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
!
!     ==================================================================
!     ==  READ MM ATOM DATA FROM STRC FILE                            ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$EXISTL(LL_STRC,'QM-MM',1,TQMMM)
      IF(TQMMM) THEN
        CALL LINKEDLIST$SELECT(LL_STRC,'QM-MM')
        CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NATMM)
        ALLOCATE(MMNAME(NATMM))
        ALLOCATE(MMR(3,NATMM)) ;MMR(:,:)=0.D0
        ALLOCATE(MMQ(NATMM)) ;MMQ(:)=0.D0
        DO IAT=1,NATMM
          CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IAT)
          CALL LINKEDLIST$GET(LL_STRC,'NAME',1,MMNAME(IAT))
          CALL LINKEDLIST$EXISTD(LL_STRC,'R',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_STRC,'R',1,MMR(:,IAT))
            MMR(:,IAT)=MMR(:,IAT)*RUNIT
          ELSE
            CALL LINKEDLIST$GET(LL_STRC,'QMATOM',1,STRING)
            DO IAT1=1,NAT
              IF(STRING.EQ.NAME(IAT1)) THEN
                MMR(:,IAT)=R(:,IAT1)
                EXIT
              END IF
            ENDDO
          END IF
          CALL LINKEDLIST$SELECT(LL_STRC,'..')
        ENDDO  
      END IF
!
!     ==================================================================
!     ==  READ SH ATOM DATA FROM STRC FILE                            ==
!     ==================================================================
!
!     ==================================================================
!     == CONVERT DATA TO ANGSTROM AND ELECTRON CHARGES                ==
!     ==================================================================
      RBAS=RBAS/ANGSTROM
      R=R/ANGSTROM
      Q=-Q
      IF(TQMMM) THEN
        MMR=MMR/ANGSTROM
      END IF
!
!     ==================================================================
!     == WRITE CSSR FILE                                              ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('CSSR',NFIL)
      CALL WRITECSSR(NFIL,OBJECTNAME,RBAS,NAT,NAME,R,Q,TCRYSTAL)
      IF(TQMMM) THEN
        CALL FILEHANDLER$UNIT('MMCSSR',NFIL)
        CALL WRITECSSR(NFIL,TRIM(OBJECTNAME)//-'MM',RBAS,NATMM,MMNAME,MMR,MMQ,TCRYSTAL)
      END IF

      IF(TCRYSTAL) THEN
        WRITE(*,FMT='("CRYSTAL OUTPUT PRODUCED")')
      ELSE
        WRITE(*,FMT='("MOLECULAR OUTPUT PRODUCED")')
      END IF
      call filehandler$closeall
      WRITE(*,FMT='("======= TASK FINISHED ========")')
      STOP
      END
!     ..................................................................
      SUBROUTINE WRITECSSR(NFIL,OBJECTNAME,RBAS,NAT,NAME,R,Q,TCRYSTAL)
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTERFACE OPERATOR (.DYAD.)
        FUNCTION DYADISCHES_PRODUCT(R1,R2) RESULT(R3)
          REAL(8), INTENT(IN) :: R1(3)
          REAL(8), INTENT(IN) :: R2(3)
          REAL(8)             :: R3(3)
        END FUNCTION DYADISCHES_PRODUCT
      END INTERFACE 
      INTEGER(4)  ,INTENT(IN) :: NFIL
      LOGICAL(4)  ,INTENT(IN) :: TCRYSTAL
      CHARACTER(*),INTENT(IN) :: OBJECTNAME
      INTEGER(4)  ,INTENT(IN) :: NAT
!     REAL(8)     ,INTENT(IN) :: RBAS(3,3)
      REAL(8)                 :: RBAS(3,3)
      CHARACTER(*),INTENT(IN) :: NAME(NAT)
      REAL(8)                 :: R(3,NAT)
!     REAL(8)     ,INTENT(IN) :: R(3,NAT)
      REAL(8)     ,INTENT(IN) :: Q(NAT)
      REAL(8)                 :: PI
      INTEGER(4)              :: NEIGH(8,NAT)
      REAL(8)                 :: RBASINV(3,3) ! RBAS**(-1)
      REAL(8)                 :: A,B,C            ! LENGTH OF LATTICE VECTORS
      REAL(8)                 :: ALPHA,BETA,GAMMA ! ANGLES BETWEEN LATTICE VECTORS
      REAL(8)                 :: DET
      INTEGER(4)              :: IAT,I,J
      REAL(8)                 :: VEC(3),SVAR,RBASNEU(3,3)
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      NEIGH(:,:)=0
!
!     ==================================================================
!     == CONVERT DATA TO UNITS OF LATTICE VECTORS                     ==
!     ==================================================================
      IF(TCRYSTAL) THEN
        RBASINV(1,:)=RBAS(:,2).DYAD.RBAS(:,3)
        RBASINV(2,:)=RBAS(:,3).DYAD.RBAS(:,1)
        RBASINV(3,:)=RBAS(:,1).DYAD.RBAS(:,2)
        DET=DOT_PRODUCT(RBAS(:,1),RBASINV(1,:))
        RBASINV=RBASINV/DET
        A=DSQRT(DOT_PRODUCT(RBAS(:,1),RBAS(:,1)))
        B=DSQRT(DOT_PRODUCT(RBAS(:,2),RBAS(:,2)))
        C=DSQRT(DOT_PRODUCT(RBAS(:,3),RBAS(:,3)))
        GAMMA =ACOS(DOT_PRODUCT(RBAS(:,1),RBAS(:,2))/(A*B))
        ALPHA =ACOS(DOT_PRODUCT(RBAS(:,2),RBAS(:,3))/(B*C))
        BETA  =ACOS(DOT_PRODUCT(RBAS(:,3),RBAS(:,1))/(C*A))
        ALPHA=180.D0/PI*ALPHA
        BETA =180.D0/PI*BETA
        GAMMA=180.D0/PI*GAMMA
        RBASNEU(:,:)=0.D0
        RBASNEU(3,3)=C
        RBASNEU(3,2)=B*COS(ALPHA*PI/180.D0)
        RBASNEU(2,2)=SQRT(B**2-RBASNEU(3,2)**2)
        RBASNEU(3,1)=A*COS(BETA*PI/180.D0)
        RBASNEU(2,1)=(A*B*COS(GAMMA*PI/180.D0)-RBASNEU(3,1)*RBASNEU(3,2))/RBASNEU(2,2)
        RBASNEU(1,1)=SQRT(A**2-RBASNEU(2,1)**2-RBASNEU(3,1)**2)
      END IF
!
!     ==================================================================
!     == WRITE CSSR FILE                                              ==
!     ==================================================================
      IF(TCRYSTAL) THEN
        WRITE(NFIL,FMT='(T39,3F8.3 &
     &   /T22,3F8.3,T50,"SPGR = 1 P 1",T72,"OPT = 1" &
      &   /I4,''   1 CREATED BY PAW    '' &
     &   /"     0 ",A4,": ",A4)') &
     &   A,B,C,ALPHA,BETA,GAMMA,NAT,OBJECTNAME(1:4),OBJECTNAME(1:4)
      ELSE
        WRITE(NFIL,FMT='(//I4,''   1 CREATED BY PAW    '' &
     &   /"     0 ",A4,": ",A4)')NAT,OBJECTNAME(1:4),OBJECTNAME(1:4)
      END IF
      DO IAT=1,NAT
        IF(TCRYSTAL) THEN
          VEC=R(:,IAT)
          VEC=MATMUL(RBASINV,VEC)+100.D0
          VEC=MOD(VEC,1.D0)
          VEC=MATMUL(RBASNEU,VEC)
          WRITE(NFIL,FMT='(I4,1X,A5,3F10.5,1X,8I4,F8.3)') &
     &           IAT,NAME(IAT),VEC(:),NEIGH(:,IAT),Q(IAT)

        ELSE
           WRITE(NFIL,FMT='(I4,1X,A5,3F10.5,1X,8I4,F8.3)') &
      &           IAT,NAME(IAT),R(:,IAT),NEIGH(:,IAT),Q(IAT)
         END IF
      ENDDO

      RETURN
      END
!     ..................................................................
      FUNCTION DYADISCHES_PRODUCT(R1,R2) RESULT(R3)
        REAL(8), INTENT(IN) :: R1(3)
        REAL(8), INTENT(IN) :: R2(3)
        REAL(8)             :: R3(3)
        R3(1)=R1(2)*R2(3)-R1(3)*R2(2)
        R3(2)=R1(3)*R2(1)-R1(1)*R2(3)
        R3(3)=R1(1)*R2(2)-R1(2)*R2(1)
      END FUNCTION DYADISCHES_PRODUCT




