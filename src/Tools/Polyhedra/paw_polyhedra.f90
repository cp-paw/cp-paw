      PROGRAM MAIN
!     **************************************************************************
!     ** PAW TOOL: PAW_POLYHEDRA                                              **
!     **                                                                      **
!     ** THIS TOOL SHALL DIVIDE THE STRUCTURE INTO POLYHEDRA                  **
!     ** AND REPORT THEIR MAIN PARAMETERS                                     **
!     **                                                                      **
!     ** CAUTION: CURRENTLY DEDICATED TO ANALYZE OCTAHEDRA IN MANGANITES      **
!     **                                                                      **
!     ** 1) exctracts cluster of atoms surrounding atoms named 'Mn...'        **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE)             :: LL_STRC
      INTEGER(4)                :: NFIL 
      INTEGER(4)                :: NFILO    !PROTOCOL-FILE UNIT 
      INTEGER(4)                :: NARGS 
      LOGICAL(4)                :: TCHK
      CHARACTER(64)             :: ROOTNAME
      REAL(8)                   :: RBAS(3,3)  !LATTICE CONSTANTS
      INTEGER(4)                :: NAT  
      REAL(8)                   :: LUNIT  !LENGTH UNIT OF INPUT STRUCTURE FILE
      REAL(8)                   :: ANGSTROM
      CHARACTER(32),ALLOCATABLE :: NAME(:) ! ATOM NAMES
      REAL(8)      ,ALLOCATABLE :: R(:,:)  ! ATOMIC POSITIONS
      INTEGER(4)                :: IAT
      REAL(8)                   :: TRANSFORM(3,3)
!     **************************************************************************
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT
!
!     ==========================================================================
!     ==  READ ATOMIC STRUCTURE                                               ==
!     ==========================================================================
      CALL LINKEDLIST$NEW(LL_STRC)
!
!     ==========================================================================
!     == GET FILE NAME ROOT FROM THE ARGUMENT LIST AND CONSTRUCT              ==
!     == FILE NAMES                                                           ==
!     ==========================================================================
      NARGS=COMMAND_ARGUMENT_COUNT()
      IF(NARGS.NE.1) THEN
        CALL ERROR$MSG('INCORRECT NUMBER OF ARGUMENTS GIVEN')
        CALL ERROR$STOP('MAIN')
      END IF
      CALL GET_COMMAND_ARGUMENT(NARGS,ROOTNAME) !LAST ARGUMENT IS THE ROOT NAME
      WRITE(*,FMT='("ROOTNAME: ",A)')TRIM(ROOTNAME)
      IF(LEN(TRIM(ROOTNAME)).EQ.0) THEN
        STOP 'NO ROOTNAME SUPPLIED'
      END IF
!
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL FILEHANDLER$SETFILE('PROT',.TRUE.,-'.POLYPROT')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','FORM','FORMATTED')
!
      CALL FILEHANDLER$SETFILE('STRC',.TRUE.,-'.STRC_OUT')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','FORM','FORMATTED')
!
!     ==========================================================================
!     == WARNING: DO NOT USE THIS TOOL WITHOUT READING THE CODE!
!     == IT CONTAINS A NUMBER OF IMPLICIT ASSUMPTIONS, SUCH AS ALIGNMENT OF 
!     == OCTAHEDRA ALONG CARTESIAN COORDINATES!
!     ==========================================================================
      WRITE(*,FMT='(80("="))')
      WRITE(*,*)'WARNING FROM PAW_POLYHEDRA:'
      WRITE(*,*)'DO NOT USE THIS TOOL WITHOUT HAVING READ THE CODE!'
      WRITE(*,*)'IT CONTAINS A NUMBER OF IMPLICIT ASSUMPTIONS'
      WRITE(*,*)'SUCH AS ALIGNMENT OF OCTAHEDRA ALONG CARTESIAN COORDINATES!'
      WRITE(*,FMT='(80("="))')
      CALL FILEHANDLER$UNIT('PROT',NFIL)
      WRITE(NFIL,FMT='(80("="))')
      WRITE(NFIL,*)'WARNING FROM PAW_POLYHEDRA:'
      WRITE(NFIL,*)'DO NOT USE THIS TOOL WITHOUT HAVING READ THE CODE!'
      WRITE(NFIL,*)'IT CONTAINS A NUMBER OF IMPLICIT ASSUMPTIONS'
      WRITE(NFIL,*)'SUCH AS ALIGNMENT OF OCTAHEDRA ALONG CARTESIAN COORDINATES!'
      WRITE(NFIL,FMT='(80("="))')
!
!     ==========================================================================
!     == READ STRUCTURE FILE                                                  ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('STRC',NFIL)
      CALL LINKEDLIST$READ(LL_STRC,NFIL,'~')
!
!     ==========================================================================
!     == GET LENGTH UNIT                                                      ==
!     ==========================================================================
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'GENERIC')

      CALL LINKEDLIST$EXISTD(LL_STRC,'LUNIT[AA]',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_STRC,'LUNIT[AA]',1,LUNIT)
        LUNIT=LUNIT*ANGSTROM
      ELSE
        CALL LINKEDLIST$GET(LL_STRC,'LUNIT',1,LUNIT)
      END IF
!
!     ==========================================================================
!     == GET LATTICE VECTORS                                                  ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'LATTICE')
      CALL LINKEDLIST$GET(LL_STRC,'T',1,RBAS)
      RBAS=RBAS*LUNIT
!
!     ==========================================================================
!     ==  READ ATOM DATA FROM STRC FILE                                       ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NAT)
      ALLOCATE(NAME(NAT))
      ALLOCATE(R(3,NAT))
      DO IAT=1,NAT
        CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IAT)
        CALL LINKEDLIST$GET(LL_STRC,'R',1,R(:,IAT))
        R(:,IAT)=R(:,IAT)*LUNIT
        CALL LINKEDLIST$GET(LL_STRC,'NAME',1,NAME(IAT))
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO

!!$! THIS IS TO ROTATE THE AXES
!!$TRANSFORM(:,:)=0.D0
!!$TRANSFORM(1:2,1:2)=1.D0/SQRT(2.D0)
!!$TRANSFORM(1,2)=-TRANSFORM(1,2)
!!$TRANSFORM(3,3)=1.D0
!!$RBAS=MATMUL(TRANSFORM,RBAS)
!!$R=MATMUL(TRANSFORM,R)
!
!     =========================================================================
!     ==  WRITE ATOMIC STRUCTURE                                             ==
!     =========================================================================
      WRITE(*,FMT='(82("="),T10,"  ATOMIC STRUCTURE IN ANGSTROM")')
      WRITE(*,FMT='("T1 ",3F10.5)')RBAS(:,1)/ANGSTROM
      WRITE(*,FMT='("T2 ",3F10.5)')RBAS(:,2)/ANGSTROM
      WRITE(*,FMT='("T3 ",3F10.5)')RBAS(:,3)/ANGSTROM
      DO IAT=1,NAT
        WRITE(*,FMT='("TYPE ",A6," R=",3F10.5)')NAME(IAT),R(:,IAT)/ANGSTROM
      ENDDO
!
!     =========================================================================
!     ==  extract clusters of nearest neigbbors for all atoms named 'MN...'  ==
!     =========================================================================
      DO IAT=1,NAT
        IF(NAME(IAT)(1:2).EQ.'MN') THEN
          CALL EXTRACTCLUSTERS(IAT,RBAS,NAT,R)
        END IF
      ENDDO
      CALL ERROR$NORMALSTOP()
      STOP
      END
!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXTRACTCLUSTERS(IATC,RBAS,NAT,R)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IATC
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      REAL(8)   ,PARAMETER  :: ANGSTROM=1.D0/.529177D0
      REAL(8)               :: RC(3,NAT)
      INTEGER(4)            :: IAT
      INTEGER(4)            :: IT1,IT2,IT3
      INTEGER(4)            :: IC
      INTEGER(4),PARAMETER  :: NCLUSTERX=7
      REAL(8)               :: RCLUSTER(3,NCLUSTERX)
      REAL(8)               :: DCLUSTER(NCLUSTERX)
      INTEGER(4)            :: NCLUSTER
      REAL(8)               :: RT(3)  ! LATTICE TRANSLATION
      REAL(8)               :: R2(3)  ! RELATIVE POSITION OF NEIGHBOR
      REAL(8)               :: D2     ! DISTANCE OF NEIGHBOR
      REAL(8)               :: PI
!     *************************************************************************
      PI=4.D0*ATAN(1.D0)
!
!     =========================================================================
!     ==  DETERMINE CLUSTER OF NEAREST NEIGHBORS                             ==
!     =========================================================================
      DO IAT=1,NAT
        RC(:,IAT)=R(:,IAT)-R(:,IATC)
      ENDDO
      NCLUSTER=1
      RCLUSTER(:,1)=(/0.D0,0.D0,0.D0/)
      DCLUSTER(1)=0.D0
      DO IT1=-1,1
        DO IT2=-1,1
          DO IT3=-1,1
            RT=MATMUL(RBAS(:,:),REAL((/IT1,IT2,IT3/)))
            DO IAT=1,NAT
              R2(:)=RC(:,IAT)+RT(:)
              D2=SQRT(SUM(R2**2))
              IF(D2.EQ.0.D0) CYCLE
              DO IC=NCLUSTER,1,-1
                IF(D2.LT.DCLUSTER(IC)) THEN
                  IF(IC+1.LE.NCLUSTERX) THEN
                    NCLUSTER=MAX(NCLUSTER,IC+1)
                    RCLUSTER(:,IC+1)=RCLUSTER(:,IC)
                    DCLUSTER(IC+1)=DCLUSTER(IC)
                  END IF
                ELSE ! D2>DCLUSTER(IC)
                  IF(IC+1.LE.NCLUSTERX) THEN
                    RCLUSTER(:,IC+1)=R2(:)
                    DCLUSTER(IC+1)=D2
                    NCLUSTER=MAX(NCLUSTER,IC+1)
                  END IF
                  EXIT
                END IF  
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     =========================================================================
!     ==  WRITE ATOMIC STRUCTURE  OF THE CLUSTER                             ==
!     =========================================================================
      WRITE(*,FMT='(82("="),T10,"  CLUSTER FOR ATOM ",I4,"  ")')IATC
      DO IAT=1,NCLUSTER
        WRITE(*,FMT='(" R=",3F10.5," AA D=",F10.5," AA")') &
     &                        RCLUSTER(:,IAT)/ANGSTROM,DCLUSTER(IAT)/ANGSTROM
      ENDDO
!
!     =========================================================================
!     ==  EXTRACT NORM OCTAHEDRON                                            ==
!     =========================================================================
      CALL OCTAPARAMETERS(RCLUSTER(:,2:7))
      RETURN
      END
!!$!
!!$!     ..1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE OCTAcanonicalorder(Rin,rout)
!!$!     *************************************************************************
!!$!     **  draft only! not finished!
!!$      implicit none
!!$      integer(4),parameter   :: nat=6
!!$      real(8)   ,intent(in)  :: rin(3,nat)
!!$      real(8)   ,intent(out) :: rout(3,nat)
!!$      real(8)                :: cosang(nat*(nat-1)/2)
!!$      integer(4)             :: ij(2)
!!$
!!$!     == normalized distance vectors ===========================================
!!$      do iat=1,nat
!!$        rout(:,iat)=rin(:,iat)/sqrt(dot_product(rin(:,iat),rin(:,iat)))
!!$      enddo
!!$!     == angles
!!$      cosang=0.d0
!!$      ind=0
!!$      do i=1,6
!!$        do j=i+1,6
!!$          ind=ind+1
!!$          cosang(ind)=dot_product(rout(:,iat1),rout(:,iat2))
!!$        enddo
!!$      enddo
!!$      iax1=minloc(cosang,dim=1)
!!$      cosang(iax1)=1.d+10
!!$      iax2=minloc(cosang,dim=1)
!!$      cosang(iax2)=1.d+10
!!$      iax3=minloc(cosang,dim=1)
!!$      cosang(iax3)=1.d+10
!!$      ind=0
!!$      do i=1,6
!!$        do j=i+1,6
!!$          ind=ind+1
!!$          if(ind.eq.iax1) then
!!$            dir(:,1)=rin(:,i)-rin(:,j)
!!$            if(sum(dir(:,1).ge.0.d0) then
!!$              rout(:,1)=rin(:,i)            
!!$              rout(:,2)=rin(:,j)            
!!$            else
!!$              rout(:,2)=rin(:,i)            
!!$              rout(:,1)=rin(:,j)            
!!$              dir(:,1)=-dir(:,1)
!!$            end if
!!$          end if
!!$          if(ind.eq.iax2) then
!!$            dir(:,2)=rin(:,i)-rout(:,j)
!!$            if(sum(dir(:,2).ge.0.d0) then
!!$              rout(:,3)=rin(:,i)            
!!$              rout(:,4)=rin(:,j)            
!!$            else
!!$              rout(:,4)=rin(:,i)            
!!$              rout(:,3)=rin(:,j)            
!!$              dir(:,2)=-dir(:,2)
!!$            end if
!!$          end if
!!$          if(ind.eq.iax3) then
!!$            dir(:,3)=rin(:,i)-rout(:,j)
!!$            if(sum(dir(:,3).ge.0.d0) then
!!$              rout(:,5)=rin(:,i)            
!!$              rout(:,6)=rin(:,j)            
!!$            else
!!$              rout(:,6)=rin(:,i)            
!!$              rout(:,5)=rin(:,j)            
!!$              dir(:,3)=-dir(:,3)
!!$            end if
!!$          end if
!!$        enddo
!!$      enddo
!!$      return
!!$      end

!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE OCTAPARAMETERS(R0)
!     *************************************************************************
!     **  EXTRACT THE PARAMETERS FOR OCTAHEDRAL DISTORTIONS FOR A GIVEN SET  **
!     **  OF SIX LIGANDS                                                     **
!     **  RESULT IS PRINTED TO THE SCREEN                                    **
!     **  PAY ATTENTION TO THE ORDER OF THE OPERATIONS!!                     **
!     *************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R0(3,6)
      REAL(8)   ,PARAMETER  :: ANGSTROM=1.D0/.529177D0
      REAL(8)               :: R(3,6)
      REAL(8)               :: FEDIS(3) ! FERRO-ELECTRIC DISPLACEMENT
      REAL(8)               :: DAV      ! AVERAGE BOND LENGTH
      REAL(8)               :: PHI(3)   ! ROTATION ANGLE
      REAL(8)               :: DNON(3,3)   ! 
      REAL(8)               :: DI(6)
      REAL(8)               :: VEC(3)
      REAL(8)               :: SVAR,AB
      INTEGER(4)            :: I,J
      REAL(8)               :: RIDEAL(3,6)
      REAL(8)               :: DIS(3,6)
      REAL(8)               :: ANGLE
      REAL(8)               :: Q2,Q3  ! JAHN TELLER MODES
      REAL(8)               :: PI
!     *************************************************************************
      PI=4.D0*ATAN(1.D0)
!     __IDEAL OCTAHEDRAL OXYGEN POSITIONS (LEFT,RIGHT,BACK,FRONT,BOTTOM,TOP)___
      RIDEAL=0.D0   
      DO I=1,3
        RIDEAL(I,2*I-1)=-1.D0
        RIDEAL(I,2*I)  =1.D0
      ENDDO
!
!     =========================================================================
!     ==  BRING ATOMS INTO CANONICAL ORDER                                   ==
!     =========================================================================
!     __LEFT___________________________________________________________________
      I=MINLOC(R0(1,:),DIM=1)  
      R(:,1)=R0(:,I)    
!     __RIGHT__________________________________________________________________
      I=MAXLOC(R0(1,:),DIM=1)    
      R(:,2)=R0(:,I)    
!     __BACK___________________________________________________________________
      I=MINLOC(R0(2,:),DIM=1)  
      R(:,3)=R0(:,I)    
!     __FRONT__________________________________________________________________
      I=MAXLOC(R0(2,:),DIM=1)    
      R(:,4)=R0(:,I)   
!     __BOTTOM_________________________________________________________________
      I=MINLOC(R0(3,:),DIM=1)  
      R(:,5)=R0(:,I)   
!     __TOP____________________________________________________________________
      I=MAXLOC(R0(3,:),DIM=1)    
      R(:,6)=R0(:,I)   
!
!     =========================================================================
!     == FERROELECTRIC DISPLACEMENT OF CENTRAL ATOM                          ==
!     =========================================================================
      VEC=0.D0
      DO I=1,6
        VEC=VEC+R(:,I)
      ENDDO
      FEDIS(:)=VEC(:)/6.D0
!     __ REMOVE FERROELECTRIC DISPLACEMENT ____________________________________
      DO I=1,6
        R(:,I)=R(:,I)-FEDIS(:)
      ENDDO
!
!     =========================================================================
!     == REMOVE AXIS ASYMMETRY.                                              ==
!     == AFTER THIS STEP THE POSITION ARE OBTAINED FROM THE IDEAL OCTAHEDRON ==
!     == BY LINEAR TRANSFORMATION (ROTATION,STRETCH,SHEAR)                   ==
!     =========================================================================
      DO I=1,3
        VEC=(R(:,2*I-1)+R(:,2*I))/2.D0
        R(:,2*I-1)=R(:,2*I-1)-VEC
        R(:,2*I)  =R(:,2*I)  -VEC
        DNON(:,I)=VEC
      ENDDO
!
!     =========================================================================
!     == EXTRACT AVERAGE BOND LENGTH                                         ==
!     =========================================================================
      DAV=0.D0
      DO I=1,6
        DI(I)=SQRT(SUM(R(:,I)**2))
        DAV=DAV+DI(I)
      ENDDO
      DAV=DAV/6.D0
!
!     =========================================================================
!     == JAHN-TELLER DISTORTIONS                                             ==
!     =========================================================================
      DO I=1,3
        DI(I)=SQRT(SUM((R(:,2*I)-R(:,2*I-1))**2))/2.D0
      ENDDO
!     __ USE CONVENTION OF DAGOTTO, "NANOSCALE PHASE SEPARATION...",2003
      Q2=(DI(1)-DI(2))/SQRT(2.D0)
      Q3=(-DI(1)-DI(2)+2.D0*DI(3))/SQRT(6.D0)
!
!     =========================================================================
!     == NORMALIZE BOND LENGTHS (THEY ARE ENCODED IN DAV,Q2,Q3)              ==
!     =========================================================================
      DO I=1,6
        R(:,I)=R(:,I)/SQRT(SUM(R(:,I)**2))
      ENDDO
!
!     =========================================================================
!     == SHEAR                                                               ==
!     =========================================================================
!     == BONDS SHOULD BE NORMALIZED AFTER EXTRACTING JAHN TELLER DISTORTIONS ==
!     == NORMALIZATION IS DONE ANY TO MAKE SURE
      DO I=1,6
        SVAR=DSQRT(SUM(R(:,I)**2))
        R(:,I)=R(:,I)/SVAR
      ENDDO       
      DO I=1,4
!       __ PROJECT OUT AXIAL COMPONENT FROM EQUATORIAL LIGANDS ________________
        R(:,I)=R(:,I)-R(:,6)*DOT_PRODUCT(R(:,6),R(:,I))
!       __ RENORMALIZE LENGTH AGAIN. PROCEDURE AMOUNTS TO BOND ROTATION________
        R(:,I)=R(:,I)/SQRT(SUM(R(:,I)**2))
      ENDDO
      R(:,1)=R(:,1)/SQRT(SUM(R(:,1)**2))
      R(:,2)=R(:,2)/SQRT(SUM(R(:,2)**2))
      R(:,3)=R(:,3)/SQRT(SUM(R(:,3)**2))
      R(:,4)=R(:,4)/SQRT(SUM(R(:,4)**2))
!     ==  CORRECT ANGLE BETWEEN AXIAL LIGANDS =================================
!     == APRIME=A+SVAR*B; BPRIME=B+SVAR*A
!     == USE A^2=B^2=1: AB+2*SVAR+AB*SVAR**2=0  (BINOMIAL FORMULA IS UNSTABLE)
!     == ITERATE SVAR=-AB/2 * (1+SVAR**2)
      AB=DOT_PRODUCT(R(:,4),R(:,2))
      SVAR=0.D0
      DO I=1,5
        SVAR=-0.5D0*AB*(1+SVAR**2)
      ENDDO
      VEC=R(:,2)
      R(:,1)=R(:,1)-R(:,4)*SVAR
      R(:,2)=R(:,2)+R(:,4)*SVAR
      R(:,3)=R(:,3)-VEC*SVAR
      R(:,4)=R(:,4)+VEC*SVAR
      R(:,1)=R(:,1)/SQRT(SUM(R(:,1)**2))
      R(:,2)=R(:,2)/SQRT(SUM(R(:,2)**2))
      R(:,3)=R(:,3)/SQRT(SUM(R(:,3)**2))
      R(:,4)=R(:,4)/SQRT(SUM(R(:,4)**2))
!
!     =========================================================================
!     == OBTAIN TILT ANGLE                                                   ==
!     =========================================================================
      DIS=R-RIDEAL
      PHI=0.D0
      DO I=1,3
        DO J=I+1,3
          CALL VECTORPRODUCT(DIS(:,2*I),DIS(:,2*J),VEC)
          IF(DOT_PRODUCT(VEC,PHI).LT.0.D0) VEC=-VEC
          PHI(:)=PHI(:)+VEC(:)
        ENDDO
      ENDDO
      ANGLE=SQRT(SUM(PHI**2))
      IF(ANGLE.NE.0.D0) THEN
        PHI(:)=PHI(:)/ANGLE
      ELSE
        PHI(3)=1.D0
      END IF
      ANGLE=0.D0
      DO I=1,6
        CALL GETANGLE(PHI,RIDEAL(:,I),R(:,I),SVAR)
        ANGLE=ANGLE+SVAR
      ENDDO
      ANGLE=ANGLE/6.D0
      PHI=PHI*ANGLE
!
!     == UNDO ROTATION ========================================================
      DO I=1,6
        CALL ROTATE(-PHI,R(:,I),VEC)
        R(:,I)=VEC
      ENDDO     
!
!     =========================================================================
!     ==  WRITE ATOMIC STRUCTURE  OF THE OCTAHEDRON                          ==
!     =========================================================================
      WRITE(*,FMT='(82("="),T20,"  REPORT ON OCTAHEDRAL PARAMETERS  ")')
      WRITE(*,FMT='("FERROELECTRIC DISPLACEMENT ",T30,3F10.5," AA")') &
     &      FEDIS/ANGSTROM
      DO I=1,3
        WRITE(*,FMT='("DNON=",T30,3F10.5," AA")')DNON(:,I)/ANGSTROM
      ENDDO
      WRITE(*,FMT='("AVERAGE BOND LENGTH",T30,F10.5," AA")')DAV/ANGSTROM
!     __ USE CONVENTION OF DAGOTTO, "NANOSCALE PHASE SEPARATION...",2003
      WRITE(*,FMT='("Q2=[D(+X)-D(+Y)]/SQRT(2)")')
      WRITE(*,FMT='("Q3=[-D(+X)-D(+Y)+2D(+Z)]/SQRT(6)")')
      WRITE(*,FMT='("JAHN-TELLER DISTORTION (Q2,Q3)",T30,2F10.5," AA")') &
     &                                                  Q2/ANGSTROM,Q3/ANGSTROM
      WRITE(*,FMT='("JAHN-TELLER AMPL. |(Q2,Q3)|",T30,F10.5," AA")') &
     &                                               SQRT(Q2**2+Q3**2)/ANGSTROM
      IF(Q2**2+Q3**2.GT.1.D-8) THEN
        SVAR=ACOS(Q2/SQRT(Q2**2+Q3**2))
      ELSE
        SVAR=0.D0
      ENDIF
      IF(Q3.LT.0.D0)SVAR=-SVAR
      SVAR=SVAR/(PI/6.D0)
      I=INT(SVAR+12)-12
      J=NINT(100*(SVAR-REAL(I,KIND=8)))
      IF(I.EQ.-6) THEN
         WRITE(*,FMT='(I3,"% (-1,1,0)  AND ",I3,"% (-1,2,-1) JT-TYPE ")')100-J,J
      ELSE IF(I.EQ.-5) THEN
         WRITE(*,FMT='(I3,"% (-1,2,-1) AND ",I3,"% (0,1,-1)  JT-TYPE ")')100-J,J
      ELSE IF(I.EQ.-4) THEN
         WRITE(*,FMT='(I3,"% (0,1,-1)  AND ",I3,"% (1,1,-2)  JT-TYPE ")')100-J,J
      ELSE IF(I.EQ.-3) THEN
         WRITE(*,FMT='(I3,"% (1,1,-2)  AND ",I3,"% (1,0,-1)  JT-TYPE ")')100-J,J
      ELSE IF(I.EQ.-2) THEN
         WRITE(*,FMT='(I3,"% (1,0,-1)  AND ",I3,"% (2,-1,-1) JT-TYPE ")')100-J,J
      ELSE IF(I.EQ.-1) THEN
         WRITE(*,FMT='(I3,"% (2,-1,-1) AND ",I3,"% (1,-1,0)  JT-TYPE ")')100-J,J
      ELSE IF(I.EQ.0) THEN
         WRITE(*,FMT='(I3,"% (1,-1,0)  AND ",I3,"% (1,-2,1)  JT-TYPE ")')100-J,J
      ELSE IF(I.EQ.1) THEN
         WRITE(*,FMT='(I3,"% (1,-2,1)  AND ",I3,"% (0,-1,1)  JT-TYPE ")')100-J,J
      ELSE IF(I.EQ.2) THEN
         WRITE(*,FMT='(I3,"% (1,-1,1)  AND ",I3,"% (-1,-1,2) JT-TYPE ")')100-J,J
      ELSE IF(I.EQ.3) THEN
         WRITE(*,FMT='(I3,"% (-1,-1,2) AND ",I3,"% (-1,0,1)  JT-TYPE ")')100-J,J
      ELSE IF(I.EQ.4) THEN
         WRITE(*,FMT='(I3,"% (-1,0,1)  AND ",I3,"% (-2,1,1)  JT-TYPE ")')100-J,J
      ELSE IF(I.EQ.5) THEN
         WRITE(*,FMT='(I3,"% (-2,1,1)  AND ",I3,"% (-1,1,0)  JT-TYPE ")')100-J,J
      ELSE
        STOP 'SELECTION ERROR'
      END IF

      IF(SQRT(SUM(PHI**2)).GT.0.D0)  THEN
        WRITE(*,FMT='("TILT ANGLE=",T30,F10.5," DEGREE")') &
     &                     SQRT(SUM(PHI**2))/PI*180.D0
        WRITE(*,FMT='("TILT AXIS=",T30,3F10.5)')PHI/SQRT(SUM(PHI**2))
      ELSE
        WRITE(*,FMT='("TILT ANGLE=",T30,F10.5," DEGREE")')0.D0
      END IF
      RETURN
      END
!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE VECTORPRODUCT(A,B,C)
!     *************************************************************************
!     ** CONSTRUCT THE VECTOR PRODUCT C = A VECTORPRODUCT B                  **
!     *************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: A(3)
      REAL(8),INTENT(IN)  :: B(3)
      REAL(8),INTENT(OUT) :: C(3)
!     ************************************************************************
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      RETURN
      END
!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ROTATE(PHI,R,RNEW)
!     *************************************************************************
!     ** ROTATE R WITH THE ANGLE VECTOR PHI (RIGHT-HAND RULE) INTO RNEW      **
!     *************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: PHI(3)
      REAL(8),INTENT(IN) :: R(3)
      REAL(8),INTENT(OUT):: RNEW(3)
      REAL(8)            :: ANGLE
      REAL(8)            :: EPHI(3)
!     **************************************************************************
      ANGLE=SQRT(SUM(PHI**2))
      EPHI=PHI/ANGLE
      CALL VECTORPRODUCT(EPHI,R,RNEW)
      RNEW=R+RNEW*SIN(ANGLE)+(R-EPHI*DOT_PRODUCT(EPHI,R))*(COS(ANGLE)-1.D0)
      RETURN
      END
!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GETANGLE(AXIS,R,RNEW,ANGLE)
!     *************************************************************************
!     ** SPECIFY THE ROTATION ANGLE FOR A GIVEN ROTATION AXIS, THAT BRINGS   **
!     ** POINT R INTO THE TRANSFORMED VECTOR RNEW
!     *************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: AXIS(3)
      REAL(8),INTENT(IN) :: R(3) 
      REAL(8),INTENT(IN) :: RNEW(3) 
      REAL(8),INTENT(OUT):: ANGLE
      REAL(8)            :: EPHI(3) 
      REAL(8)            :: S(3)
      REAL(8)            :: SNEW(3)
      REAL(8)            :: U(3),V(3),W(3)
      REAL(8)            :: UU,UV,UW,VV,VW
      REAL(8)            :: DET
      REAL(8)            :: COSPHI,SINPHI
!     *************************************************************************
      EPHI=AXIS/SQRT(SUM(AXIS**2))
      S=R-EPHI*DOT_PRODUCT(EPHI,R)
      SNEW=RNEW-EPHI*DOT_PRODUCT(EPHI,RNEW)
      S=S/SQRT(SUM(S**2))
      SNEW=SNEW/SQRT(SUM(SNEW**2))

      U=EPHI*DOT_PRODUCT(EPHI,S)-S
      CALL VECTORPRODUCT(EPHI,S,V)
      CALL VECTORPRODUCT(EPHI,SNEW,W)
      UU=DOT_PRODUCT(U,U)
      UV=DOT_PRODUCT(U,V)
      UW=DOT_PRODUCT(U,W)
      VV=DOT_PRODUCT(V,V)
      VW=DOT_PRODUCT(V,W)
      DET=UU*VV-UV*UV
      IF(DET.EQ.0.D0) THEN
        ANGLE=0.D0
        RETURN
      END IF
      COSPHI=(UU*VW-UV*UW)/DET
      SINPHI=(VV*UW-UV*VW)/DET
      ANGLE=ACOS(COSPHI)
      IF(SINPHI.LT.0.D0) ANGLE=-ANGLE
      RETURN
      END
!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GETANGLEVEC(R,RNEW,PHI)
!     *************************************************************************
!     ** SPECIFY THE ROTATION ANGLE PHI (DIRECTION=AXIS,LENGTH=ANGLE)        **
!     ** FOR A ROTATION WITH AN AXIS PASSING THROUGH THE ORIGIN, WHICH MAPS  **
!     ** R TO RNEW. THE CHANGE IN LENGTH FROM R TO RNEW IS IGNORED           **
!     *************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: R(3) 
      REAL(8),INTENT(IN) :: RNEW(3) 
      REAL(8),INTENT(OUT):: PHI(3)
      REAL(8)            :: X(3),XNEW(3)
      REAL(8)            :: EPHI(3)
      REAL(8)            :: U(3),V(3),W(3)
      REAL(8)            :: UU,UV,UW,VV,VW
      REAL(8)            :: DET
      REAL(8)            :: COSPHI,SINPHI
      REAL(8)            :: ANGLE
!     *************************************************************************
      X=R/DSQRT(SUM(R**2))
      XNEW=RNEW/DSQRT(SUM(RNEW**2))
      CALL VECTORPRODUCT(X,XNEW,EPHI)
      EPHI=EPHI/DSQRT(SUM(EPHI**2))
      U=EPHI*DOT_PRODUCT(EPHI,X)-X
      CALL VECTORPRODUCT(EPHI,X,V)
      CALL VECTORPRODUCT(EPHI,XNEW,W)
      UU=DOT_PRODUCT(U,U)
      UV=DOT_PRODUCT(U,V)
      UW=DOT_PRODUCT(U,W)
      VV=DOT_PRODUCT(V,V)
      VW=DOT_PRODUCT(V,W)
      DET=(UU*VV-UV*UV)
      IF(DET.EQ.0.D0) THEN
        PHI=0.D0
        RETURN
      END IF
      COSPHI=(UU*VW-UV*UW)/DET
      SINPHI=(VV*UW-UV*VW)/DET
      ANGLE=ACOS(COSPHI)
      IF(SINPHI.LT.0.D0) ANGLE=-ANGLE
      PHI=EPHI*ANGLE
      RETURN
      END
