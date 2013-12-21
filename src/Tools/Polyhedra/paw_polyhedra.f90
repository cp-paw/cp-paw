      Program main
!     **************************************************************************
!     ** paw tool: paw_polyhedra                                              **
!     **                                                                      **
!     ** This tool shall divide the structure into polyhedra                  **
!     ** and report their main parameters                                     **
!     **                                                                      **
!     ** Caution: Currently dedicated to analyze octahedra in manganites      **
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      implicit none
      type(ll_type)             :: ll_strc
      integer(4)                :: nfil 
      integer(4)                :: nfilo    !protocoll-file unit 
      integer(4)                :: nargs 
      logical(4)                :: tchk
      character(64)             :: rootname
      real(8)                   :: rbas(3,3)  !lattice constants
      integer(4)                :: nat  
      real(8)                   :: lunit  !length unit of input structure file
      real(8)                   :: angstrom
      character(32),allocatable :: name(:) ! atom names
      real(8)      ,allocatable :: r(:,:)  ! atomic positions
      integer(4)                :: iat
      real(8)                   :: transform(3,3)
!     **************************************************************************
!
!     ==========================================================================
!     ==  read atomic structure                                               ==
!     ==========================================================================
      CALL LINKEDLIST$NEW(LL_STRC)
!
!     ==========================================================================
!     == GET FILE NAME ROOT FROM THE ARGUMENT LIST AND CONSTRUCT              ==
!     == FILE NAMES                                                           ==
!     ==========================================================================
      CALL LIB$NARGS(NARGS)
      if(nargs.ne.1) then
        call error$msg('incorrect number of arguments given')
        call error$stop('main')
      end if
      CALL LIB$GETARG(NARGS,ROOTNAME) !LAST ARGUMENT IS THE ROOT NAME
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
!     == READ STRUCTURE FILE                                                  ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('STRC',NFIL)
      CALL LINKEDLIST$READ(LL_STRC,NFIL,'~')
!
!     ==========================================================================
!     == GET Length unit                                                      ==
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

!!$! this is to rotate the axes
!!$transform(:,:)=0.d0
!!$transform(1:2,1:2)=1.d0/sqrt(2.d0)
!!$transform(1,2)=-transform(1,2)
!!$transform(3,3)=1.d0
!!$rbas=matmul(transform,rbas)
!!$r=matmul(transform,r)
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
!     ==  write atomic structure                                             ==
!     =========================================================================
      do iat=1,nat
        if(name(iat)(1:2).eq.'MN') Then
          call extractclusters(iat,rbas,nat,r)
        end if
      enddo
      stop
      end
!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      subroutine extractclusters(iatc,rbas,nat,r)
      implicit none
      integer(4),intent(in) :: iatc
      integer(4),intent(in) :: nat
      real(8)   ,intent(in) :: rbas(3,3)
      real(8)   ,intent(in) :: r(3,nat)
      real(8)   ,parameter  :: angstrom=1.d0/.529177d0
      real(8)               :: rc(3,nat)
      integer(4)            :: iat
      integer(4)            :: it1,it2,it3
      integer(4)            :: ic
      integer(4),parameter  :: nclusterx=7
      real(8)               :: rcluster(3,nclusterx)
      real(8)               :: dcluster(nclusterx)
      integer(4)            :: ncluster
      REAL(8)               :: RT(3)  ! LATTICE TRANSLATION
      REAL(8)               :: r2(3)  ! relative position of neighbor
      REAL(8)               :: d2     ! distance of neighbor
      real(8)               :: pi
!     *************************************************************************
      pi=4.d0*atan(1.d0)
!
!     =========================================================================
!     ==  determine cluster of nearest neighbors                             ==
!     =========================================================================
      DO IAT=1,NAT
        RC(:,IAT)=R(:,IAT)-R(:,IATC)
      ENDDO
      ncluster=1
      rcluster(:,1)=(/0.d0,0.d0,0.d0/)
      dcluster(1)=0.d0
      DO IT1=-1,1
        DO IT2=-1,1
          DO IT3=-1,1
            RT=MATMUL(RBAS(:,:),REAL((/IT1,IT2,IT3/)))
            DO IAT=1,NAT
              R2(:)=RC(:,IAT)+RT(:)
              D2=SQRT(SUM(R2**2))
              if(d2.eq.0.d0) cycle
              DO IC=NCLUSTER,1,-1
                IF(D2.LT.DCLUSTER(IC)) THEN
                  IF(IC+1.Le.NCLUSTERX) THEN
                    NCLUSTER=MAX(NCLUSTER,IC+1)
                    RCLUSTER(:,IC+1)=RCLUSTER(:,IC)
                    DCLUSTER(IC+1)=DCLUSTER(IC)
                  END IF
                ELSE ! D2>DCLUSTER(IC)
                  IF(IC+1.Le.NCLUSTERX) THEN
                    RCLUSTER(:,IC+1)=R2(:)
                    DCLUSTER(IC+1)=D2
                    ncluster=max(ncluster,ic+1)
                  END IF
                  exit
                END IF  
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     =========================================================================
!     ==  write atomic structure  of the cluster                             ==
!     =========================================================================
      write(*,fmt='(82("="),t10,"  cluster for atom ",i4,"  ")')iatc
      do iat=1,ncluster
        write(*,fmt='(" r=",3f10.5," AA d=",f10.5," AA")') &
     &                        rcluster(:,iat)/angstrom,dcluster(iat)/angstrom
      enddo
!
!     =========================================================================
!     ==  extract norm octahedron                                            ==
!     =========================================================================
      call octaparameters(rcluster(:,2:7))
      RETURN
      END
!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      subroutine octaparameters(r0)
!     *************************************************************************
!     **  extract the parameters for octahedral distortions for a given set  **
!     **  of six ligands                                                     **
!     **  result is printed to the screen                                    **
!     **  pay attention to the order of the operations!!                     **
!     *************************************************************************
      implicit none
      real(8)   ,intent(in) :: r0(3,6)
      real(8)   ,parameter  :: angstrom=1.d0/.529177d0
      real(8)               :: r(3,6)
      real(8)               :: fedis(3) ! ferro-electric displacement
      real(8)               :: dav      ! average bond length
      real(8)               :: phi(3)   ! rotation angle
      real(8)               :: dnon(3,3)   ! 
      real(8)               :: di(6)
      real(8)               :: vec(3)
      real(8)               :: svar,ab
      integer(4)            :: i,j
      real(8)               :: rideal(3,6)
      real(8)               :: dis(3,6)
      real(8)               :: angle
      real(8)               :: Q2,Q3  ! Jahn Teller modes
      real(8)               :: pi
!     *************************************************************************
      pi=4.d0*atan(1.d0)
!     __ideal octahedral oxygen positions (left,right,back,front,bottom,top)___
      rideal=0.d0   
      do i=1,3
        rideal(i,2*i-1)=-1.d0
        rideal(i,2*i)  =1.d0
      enddo
!
!     =========================================================================
!     ==  bring atoms into canonical order                                   ==
!     =========================================================================
!     __left___________________________________________________________________
      i=minloc(r0(1,:),dim=1)  
      r(:,1)=r0(:,i)    
!     __right__________________________________________________________________
      i=maxloc(r0(1,:),dim=1)    
      r(:,2)=r0(:,i)    
!     __back___________________________________________________________________
      i=minloc(r0(2,:),dim=1)  
      r(:,3)=r0(:,i)    
!     __front__________________________________________________________________
      i=maxloc(r0(2,:),dim=1)    
      r(:,4)=r0(:,i)   
!     __bottom_________________________________________________________________
      i=minloc(r0(3,:),dim=1)  
      r(:,5)=r0(:,i)   
!     __top____________________________________________________________________
      i=maxloc(r0(3,:),dim=1)    
      r(:,6)=r0(:,i)   
!
!     =========================================================================
!     == FERROELECTRIC DISPLACEMENT OF CENTRAL ATOM                          ==
!     =========================================================================
      vec=0.D0
      DO I=1,6
        vec=vec+R(:,I)
      ENDDO
      fedis(:)=vec(:)/6.D0
!     __ remove ferroelectric displacement ____________________________________
      DO I=1,6
        R(:,I)=R(:,I)-fedis(:)
      ENDDO
!
!     =========================================================================
!     == remove axis ASYMMETRY.                                              ==
!     == after this step the position are obtained from the ideal octahedron ==
!     == by linear transformation (rotation,stretch,shear)                   ==
!     =========================================================================
      DO I=1,3
        vec=(R(:,2*I-1)+R(:,2*I))/2.D0
        R(:,2*I-1)=R(:,2*I-1)-vec
        R(:,2*I)  =R(:,2*I)  -vec
        dnon(:,i)=vec
      ENDDO
!
!     =========================================================================
!     == extract average bond length                                         ==
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
!     __ use convention of dagotto, "nanoscale phase separation...",2003
      Q2=(DI(1)-DI(2))/sqrt(2.d0)
      Q3=(-DI(1)-DI(2)+2.D0*DI(3))/sqrt(6.d0)
!
!     =========================================================================
!     == NORMALIZE BOND LENGTHS (THEY ARE ENCODED IN DAV,Q2,Q3)              ==
!     =========================================================================
      DO I=1,6
        R(:,i)=R(:,I)/sqrt(sum(r(:,i)**2))
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
      ab=DOT_PRODUCT(R(:,4),R(:,2))
      svar=0.d0
      do i=1,5
        svar=-0.5d0*ab*(1+svar**2)
      enddo
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
      angle=SQRT(SUM(PHI**2))
      if(angle.ne.0.d0) then
        PHI(:)=PHI(:)/angle
      else
        phi(3)=1.d0
      end if
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
!     __ use convention of dagotto, "nanoscale phase separation...",2003
      WRITE(*,FMT='("Q2=[D(+X)-D(+Y)]/SQRT(2)")')
      WRITE(*,FMT='("Q3=[-D(+X)-D(+Y)+2D(+Z)]/sqrt(6)")')
      WRITE(*,FMT='("JAHN-TELLER DISTORTION (Q2,Q3)",T30,2F10.5," AA")') &
     &                                                  Q2/ANGSTROM,Q3/ANGSTROM
      WRITE(*,FMT='("JAHN-TELLER AMPL. |(Q2,Q3)|",T30,F10.5," AA")') &
     &                                               SQRT(Q2**2+Q3**2)/ANGSTROM
      SVAR=ACOS(Q2/SQRT(Q2**2+Q3**2))
      IF(Q3.LT.0.D0)SVAR=-SVAR
      SVAR=SVAR/(PI/6.D0)
      I=INT(SVAR+12)-12
      J=NINT(100*(SVAR-REAL(I,KIND=8)))
      IF(I.EQ.-6) THEN
         WRITE(*,FMT='(i3,"% (-1,1,0)  and ",i3,"% (-1,2,-1) JT-type ")')100-j,j
      else if(i.eq.-5) then
         WRITE(*,FMT='(i3,"% (-1,2,-1) and ",i3,"% (0,1,-1)  JT-type ")')100-j,j
      else if(i.eq.-4) then
         WRITE(*,FMT='(i3,"% (0,1,-1)  and ",i3,"% (1,1,-2)  JT-type ")')100-j,j
      else if(i.eq.-3) then
         WRITE(*,FMT='(i3,"% (1,1,-2)  and ",i3,"% (1,0,-1)  JT-type ")')100-j,j
      else if(i.eq.-2) then
         WRITE(*,FMT='(i3,"% (1,0,-1)  and ",i3,"% (2,-1,-1) JT-type ")')100-j,j
      else if(i.eq.-1) then
         WRITE(*,FMT='(i3,"% (2,-1,-1) and ",i3,"% (1,-1,0)  JT-type ")')100-j,j
      else if(i.eq.0) then
         WRITE(*,FMT='(i3,"% (1,-1,0)  and ",i3,"% (1,-2,1)  JT-type ")')100-j,j
      else if(i.eq.1) then
         WRITE(*,FMT='(i3,"% (1,-2,1)  and ",i3,"% (0,-1,1)  JT-type ")')100-j,j
      else if(i.eq.2) then
         WRITE(*,FMT='(i3,"% (1,-1,1)  and ",i3,"% (-1,-1,2) JT-type ")')100-j,j
      else if(i.eq.3) then
         WRITE(*,FMT='(i3,"% (-1,-1,2) and ",i3,"% (-1,0,1)  JT-type ")')100-j,j
      else if(i.eq.4) then
         WRITE(*,FMT='(i3,"% (-1,0,1)  and ",i3,"% (-2,1,1)  JT-type ")')100-j,j
      else if(i.eq.5) then
         WRITE(*,FMT='(i3,"% (-2,1,1)  and ",i3,"% (-1,1,0)  JT-type ")')100-j,j
      else
        stop 'selection error'
      end if

      WRITE(*,FMT='("TILT ANGLE=",T30,F10.5," DEGREE")')sqrt(sum(phi**2))/PI*180.D0
      if(SQRT(SUM(PHI**2)).gt.0.d0)  then
        WRITE(*,FMT='("TILT AXIS=",T30,3F10.5)')PHI/SQRT(SUM(PHI**2))
      end if
      RETURN
      END
!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE VECTORPRODUCT(A,B,C)
!     *************************************************************************
!     ** construct the vector product C = a vectorproduct B                  **
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
      subroutine rotate(phi,r,rnew)
!     *************************************************************************
!     ** rotate r with the angle vector phi (right-hand rule) into rnew      **
!     *************************************************************************
      implicit none
      real(8),intent(in) :: phi(3)
      real(8),intent(in) :: r(3)
      real(8),intent(out):: rnew(3)
      real(8)            :: angle
      real(8)            :: ephi(3)
!     **************************************************************************
      angle=sqrt(sum(phi**2))
      ephi=phi/angle
      call vectorproduct(ephi,r,rnew)
      rnew=r+rnew*sin(angle)+(r-ephi*dot_product(ephi,r))*(cos(angle)-1.d0)
      return
      end
!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GETANGLE(AXIS,R,RNEW,ANGLE)
!     *************************************************************************
!     ** specify the rotation angle for a given rotation axis, that brings   **
!     ** point r into the transformed vector rnew
!     *************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: AXIS(3)
      REAL(8),INTENT(IN) :: R(3) 
      REAL(8),INTENT(IN) :: RNEW(3) 
      REAL(8),INTENT(OUT):: ANGLE
      REAL(8)            :: EPHI(3) 
      REAL(8)            :: s(3)
      REAL(8)            :: snew(3)
      real(8)            :: u(3),v(3),w(3)
      real(8)            :: uu,uv,uw,vv,vw
      real(8)            :: det
      real(8)            :: cosphi,sinphi
!     *************************************************************************
      EPHI=AXIS/SQRT(SUM(AXIS**2))
      s=r-ephi*dot_product(ephi,r)
      snew=rnew-ephi*dot_product(ephi,rnew)
      s=s/sqrt(sum(s**2))
      snew=snew/sqrt(sum(snew**2))

      U=EPHI*DOT_PRODUCT(EPHI,s)-s
      call VECTORPRODUCT(EPHI,s,v)
      call VECTORPRODUCT(EPHI,sNEW,w)
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
      REAL(8)            :: x(3),xnew(3)
      REAL(8)            :: EPHI(3)
      REAL(8)            :: U(3),V(3),W(3)
      REAL(8)            :: UU,UV,UW,VV,VW
      REAL(8)            :: DET
      REAL(8)            :: COSPHI,SINPHI
      REAL(8)            :: angle
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
        phi=0.d0
        RETURN
      END IF
      COSPHI=(UU*VW-UV*UW)/DET
      SINPHI=(VV*UW-UV*VW)/DET
      ANGLE=ACOS(COSPHI)
      IF(SINPHI.LT.0.D0) ANGLE=-ANGLE
      phi=ephi*angle
      RETURN
      END
