MODULE DIMER_MODULE
!USE DIMER_CONSTR_MODULE
REAL(8)                       :: D ! CONSTRAINED DISTANCE BETWEEN DIMERPOINTS
INTEGER                       :: CALCMULTIPLIERITERMAX ! EMERGENCY EXIT FOR THE ITERATION LOOP AFTER X STEPS
REAL(8)                       :: DLAMBDA !EXACTNESS OF LAMBDA (= LAGRANGE MULTIPLIER) CALCULATION
INTEGER                       :: CALCVELOCITYITERMAX !EMERGENCY EXIT FOR THE ITERATION LOOP AFTER X STEPS
REAL(8)                       :: DVELOCITY !EXACTNESS OF VELOCITY iteration 
real(8)                       :: energytra !write the energy/distance every n steps
real(8)                       :: lasttra=0.0d0
real(8),allocatable            :: RTS(:)      !coordinates of the transition state fromcntl-file               
real(8),allocatable            :: DIRM(:)

integer(4)                    :: lcs       !lengthcountershouldbe how many steps with no change of dimercenter before shorten sqd
integer(4)                    :: lc =0     ! the actuak steps without change
real(8)                       :: RCDIFFMIN !difference in dimer center position below lc=lc+1
real(8)                       :: RCDIFF    !difference in dimer center position
real(8),allocatable,save            :: RC(:)        !dimer center position
real(8),allocatable,save            :: RCOLD(:)     !dimer center position in last step
real(8)                       :: DSTEP     !d=d-dstep
integer(4)                    :: NSTEPS
real(8)                       :: stepfact
logical(4)                    :: tinitstepfact=.false.
real(8)                       :: DMIN    !the minimal dimer length (no further reducement of d)
real(8)                       :: wdownfact
real(8),allocatable           :: F1W(:) !the weighted force
logical(4)                    :: FIRSTTRA=.true.
logical(4)                    :: dimer !do we use dimer parallelisation? set in readin
logical(4)                    :: dimerfollowdown !do we follow down with each image to one of the vicinal ground states

logical(4)                    :: DLFLEX    !flexible dimer length
logical(4)                    :: KDLENGTH    !keep the startup length of the dimer
logical(4)                    :: INHIBITUP   !inhibit upward motion of the dimer
logical(4)                    :: inhibitperp
logical(4)                    :: onlyperp
logical(4)                    :: onlyrot
logical(4)                    :: wdown       !weight down the dimer image 1

CHARACTER(32)                 :: CENTER_ID
real(8)                       :: CENTER_COORD(3)
real(8)                       :: DROT

INTEGER                       :: DIM ! DIMONSION =NAT*3
REAL(8),save                  :: SQD     !SQUARE CONSTRAINED DISTANCE BETWEEN DIMERPOINTS
REAL(8)                       :: SQDR=0.0d0     !SQUARE CONSTRAINED DISTANCE BETWEEN DIMERPOINTS
REAL(8)                       :: ALPHA !FRICTION; good:0.0156D0
LOGICAL(4)                    :: PLACEDIMER=.false.

logical(4)                    :: TFIRSTSPLIT=.true.
real(8),allocatable           :: RKEEP(:) 
real(8),allocatable           :: R2split(:) 
logical(4)                    :: STRETCH
real(8)                       :: STRETCHDIST


real(8)                       :: fmpara,fmperp,fmrot,fricpara,fricperp,fricrot
logical(4)                    :: optfricpara,optfricperp,optfricrot
logical(4)                    :: fautopara,fautoperp,fautorot
real(8)                       :: fricautopara,fricautoperp,fricautorot
real(8),allocatable           :: frotm(:),fperpm(:),fparam(:)
real(8)                       :: CONSTRSTEP !max step for constraint
real(8)                       :: ekinparam
real(8)                       :: ekinperpm
real(8)                       :: ekinrotm 
real(8)                       :: Tparam
real(8)                       :: Tperpm
real(8)                       :: Trotm 
real(8)                       :: tfact=2.d0/(3.166679d-6) !per deg. of freedom; to do: get k from paw constants

!========================================================
  TYPE DIMER_CONSTR_TYPE
     CHARACTER(32)                    :: ID        !KEYWORD
     TYPE(DIMER_CONSTR_TYPE),POINTER  :: NEXT_D      !YOUNGER BROTHER
     TYPE(DIMER_CONSTR_TYPE),POINTER  :: PREV_D      !ELDER BROTHER
  END TYPE DIMER_CONSTR_TYPE

  type(DIMER_CONSTR_TYPE),POINTER  :: ELDEST_D
  type(DIMER_CONSTR_TYPE),POINTER  :: THIS_D
!========================================================



  !alex: this is for pclimb-testing!
  logical(4)               :: climbperp,TFIRSTCLIMBSTEP
  real(8)                  :: FORCEDSTEP
  real(8), allocatable     :: PCLIMBDIR(:)
  !alex: this is for pclimb-testing!
  
  !alex: this is for angle monitoring
  real(8),allocatable      :: angle1(:)
  real(8),allocatable      :: angle2(:)
  real(8),allocatable      :: angledir(:)
  logical(4)               :: treadam=.false.
  logical(4)               :: treadamnotthere=.false.
  !alex: this is for angle monitoring
  
  !alexp: this comes from merge in new code
  integer(4)               :: ntasks
  integer(4)               :: thistask
  !alexp: this comes from merge in new code
  
  integer(4)               :: dprotfil

  real(8),allocatable      :: g1(:)
  real(8),allocatable      :: g2(:)
  integer(4)               :: steps=0

  real(8)                  :: etot,etot2 !for ts estimation
end MODULE DIMER_MODULE

module dimer_oscillator_module
  integer(4)                :: odim=3    !for the dimer: cpara,cperp,crot
  real(8),allocatable       :: oscm(:)
  real(8),allocatable       :: osc0(:)
  real(8),allocatable       :: oscp(:)
  real(8),allocatable       :: oscmass(:)  
  real(8),allocatable       :: oscanner(:)
  real(8),allocatable       :: oscc(:) !the harmonic potential
end module dimer_oscillator_module

module dimer_estimatets_module
  real(8),allocatable       :: z1m(:)                  
  real(8),allocatable       :: z2m(:)                  
  real(8),allocatable       :: em(:)                  
  real(8),allocatable       :: eorthom(:)                  
  real(8)                   :: fparam
  real(8)                   :: fperpm
end module dimer_estimatets_module




subroutine dimer_init_files()
  USE STRINGS_MODULE
  USE DIMER_MODULE
  implicit none
  character(32)                    :: id
  LOGICAL(4),PARAMETER             :: T=.TRUE.
  LOGICAL(4),PARAMETER             :: F=.FALSE.


  !=== DIMERPROTOCOLL FILE ================================================
  !=== EACH MONOMER HAS ITS OWN PROTOCOLL FILE
  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
  ID=+'DPROT'
  IF(THISTASK.GT.1)THEN
     CALL FILEHANDLER$SETFILE(ID,F,-'/DEV/NULL')
  ELSE
     CALL FILEHANDLER$SETFILE(ID,T,-'.DIMER_PROT')
  END IF
  CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
  CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
  CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
  CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')

  CALL FILEHANDLER$UNIT('DPROT',dprotfil)



  !=== DIMERFOLLOWDON ENERGY TRAJECTORY FILE ==============================
  !=== EACH MONOMER HAS ITS OWN FILE
  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
  ID=+'ETRA'
  IF(THISTASK.GT.1)THEN
     CALL FILEHANDLER$SETFILE(ID,F,-'/DEV/NULL')
  ELSE
     CALL FILEHANDLER$SETFILE(ID,T,-'.DIMER_ETRA')
  END IF
  CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
  CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
  CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
  CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')



  
  return
end subroutine dimer_init_files




!     ..................................................................
      SUBROUTINE DIMER$PROPAGATE(NAT,DT,ANNER,ANNEE &
     &                 ,RMASS,EFFEMASS,FORCE,R0,RM,RP &
     &                 ,TSTRESS,CELLFRIC,CELLKIN)
!     ******************************************************************
!     **  COMMENT                         **
!     ******************************************************************
      USE MPE_MODULE
      USE DIMER_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT           ! #(ATOMS)
      REAL(8)   ,INTENT(IN) :: DT        ! TIME STEP
      REAL(8)   ,INTENT(IN) :: ANNER         ! FRICTION ON THE ATOMS
      REAL(8)   ,INTENT(IN) :: ANNEE         ! FRICTION ON THE WAVE FUNCTIONS
      REAL(8)   ,INTENT(IN) :: RMASS(NAT)    ! BARE MASS OF THE NUCLEUS
      REAL(8)   ,INTENT(IN) :: EFFEMASS(NAT) ! MASS OF TEH ELECTRONS
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)     ! R(0)
      REAL(8)   ,INTENT(IN) :: RM(3,NAT)     ! R(-)
      REAL(8)   ,INTENT(OUT):: RP(3,NAT)     ! R(+)
      REAL(8)   ,INTENT(IN) :: FORCE(3,NAT)  ! FORCE
      LOGICAL(4),INTENT(IN) :: TSTRESS
      REAL(8)   ,INTENT(IN) :: CELLFRIC(3,3) 
      REAL(8)   ,INTENT(OUT):: CELLKIN(3,3) 
      REAL(8)               :: CELLKIN1(3,3),CELLKIN2(3,3)  
      REAL(8)               :: MM(nat*3,nat*3) !MASSMATRIX
      REAL(8)               :: R1_(nat*3),R2_(nat*3) !OLD POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)               :: R1(nat*3),R2(nat*3) !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)               :: R1P(nat*3),R2P(nat*3) !NEW POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)               :: F1(nat*3),F2(nat*3) !FORCE FROM POTENTIAL
      integer(4)            :: NVAL,WORLD_ID,NWORLDTASKS
      INTEGER(4)            :: IAT

      CELLKIN=0.D0
      CELLKIN1=0.0d0
      CELLKIN2=0.0d0


      ALPHA=ANNER!DO WE NEED THIS??? DON'T THINK SO! check it!
      IF(.NOT.TSTRESS) THEN
         call MPE$QUERY('~',NTASKS,world_id)
         call MPE$QUERY('MONOMER',NTASKS,THISTASK)

         CALL DIMER$GETPROCESSLEADER2(NVAL,WORLD_ID)
         !NVAL now equals the WORLD_ID of the 1st task in dimer2
         !WORLD_ID EQUALS THE WORLD_ID OF THE CALLING PROCESS



         if(thistask.eq.1) then
            !get total energy from last step for TS estimate
            call ENERGYLIST$RETURN('TOTAL ENERGY',etot)
            
            !only 1st task in dimer has to do the following
            DIM=nat*3 !we need this TO INITIALIZE DIM because nat is read after dimer!
            MM(:,:)=0.0d0
            do iat=1 , NAT
               MM((IAT*3)-2,(IAT*3)-2)=RMASS(IAT)-EFFEMASS(IAT)
               MM((IAT*3)-1,(IAT*3)-1)=RMASS(IAT)-EFFEMASS(IAT)
               MM((IAT*3),(IAT*3))    =RMASS(IAT)-EFFEMASS(IAT)
               
               R1_((IAT*3)-2)         =RM(1,IAT)
               R1_((IAT*3)-1)         =RM(2,IAT)
               R1_((IAT*3))           =RM(3,IAT)
               
               R1((IAT*3)-2)          =R0(1,IAT)
               R1((IAT*3)-1)          =R0(2,IAT)
               R1((IAT*3))            =R0(3,IAT)
               
               F1((IAT*3)-2)          =FORCE(1,IAT)
               F1((IAT*3)-1)          =FORCE(2,IAT)
               F1((IAT*3))            =FORCE(3,IAT)
            end do




            !give WORLD_ID 1 the position and forces on dimer2
            if(WORLD_ID.ne.1) then 
               !we are 1st task of dimer2
               call MPE$SEND('~',1,42,R1_)
               call MPE$SEND('~',1,43,R1)
               call MPE$SEND('~',1,44,F1)
               call MPE$SEND('~',1,45,etot)
               
            else
               !we are WORLD_ID 1
               call MPE$RECEIVE('~',NVAL,42,R2_)
               call MPE$RECEIVE('~',NVAL,43,R2)
               call MPE$RECEIVE('~',NVAL,44,F2)
               call MPE$RECEIVE('~',NVAL,45,etot2)
            end if


            
            if(WORLD_ID.eq.1) then
               !ONLY THE FIRST TASK MOVES THE DIMER IMAGES
               SQD=D*D !REMEMBER: SQDIMERDIST IS THE *ACTUAL* DISTANCE**2, SQD=D*D IS THE DESIRED DISTANCE**2 
               call DIMER_PROPAGATE(DT,R1,R2,R1_,R2_,MM,F1,F2,R1P,R2P,CELLKIN1,CELLKIN2)
               CELLKIN=CELLKIN1
            end if

      
         else
            !all tasks .ne.1 in dimer
         end if




  !from now on we have all tasks back

         !send 1st task in dimer2 the new position and cellkin of image 2
         if(WORLD_ID.eq.1) then 
            !we are WORLD_ID 1
            call MPE$SEND('~',NVAL,45,R2p)
            call MPE$SEND('~',NVAL,46,CELLKIN2)
         else if(WORLD_ID.eq.NVAL) then
            !we are 1st task in dimer2
            call MPE$RECEIVE('~',1,45,R1p)
            call MPE$RECEIVE('~',1,46,CELLKIN)
         end if

         if(thistask.eq.1) then
            !we are 1st task in dimer group
            
            do iat=1 , NAT
               RP(1,IAT)=R1p((IAT*3)-2)
               RP(2,IAT)=R1p((IAT*3)-1)
               RP(3,IAT)=R1p((IAT*3))    
            end do
         end if

         !send the new position, cellkin to all tasks in the group
         call MPE$BROADCAST('MONOMER',1,RP)
         call MPE$BROADCAST('MONOMER',1,CELLKIN)


      ELSE 
         
         CALL ERROR$MSG('IN DIMER PROPAGATE TSTRESS=TRUE TALK ABOUT THIS WITH PETER')
         !CALL ERROR$CHVAL('',NAME(IAT))
         CALL ERROR$STOP('DIMER$PROPAGATE')
         
      END IF

      RETURN
      
!     ******************************************************************

      end SUBROUTINE DIMER$PROPAGATE


!     ..................................................................
      SUBROUTINE DIMER_PROPAGATE(DT,R1,R2,R1_,R2_,MM,F1,F2,R1P,R2P,CELLKIN1,CELLKIN2)
!     ******************************************************************
!     **  COMMENT                         **
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)    :: DT                    ! TIME STEP
      REAL(8)   ,intent(in)    :: R1(dim)               !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,intent(in)    :: R2(dim)               !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,intent(in)    :: R1_(dim)              !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,intent(in)    :: R2_(dim)              !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,intent(in)    :: MM(dim,dim)           !MASSMATRIX
      REAL(8)   ,intent(in)    :: F1(dim)               !FORCE FROM POTENTIAL
      REAL(8)   ,intent(in)    :: F2(dim)               !FORCE FROM POTENTIAL

      REAL(8)   ,intent(out)   :: R1P(dim)              !NEW POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,intent(out)   :: R2P(dim)              !NEW POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,INTENT(OUT)   :: CELLKIN1(3,3)         ! 
      REAL(8)   ,INTENT(OUT)   :: CELLKIN2(3,3)         ! 
      REAL(8)                  :: UF1(dim)              !m^-(1/2)*FORCE FROM POTENTIAL
      REAL(8)                  :: UF2(dim)              !m^-(1/2)*FORCE FROM POTENTIAL
      REAL(8)                  :: MF1(dim)              !m^(1/2)*FORCE FROM POTENTIAL
      REAL(8)                  :: MF2(dim)              !m^(1/2)*FORCE FROM POTENTIAL
      real(8)                  :: sqdimerdist           !current (distance between dimerpoints)**2
      REAL(8)                  :: R1NC(dim),R2NC(dim)   !POSITION WITHOUT CONSTRAINS
      REAL(8)                  :: V1(DIM),V2(DIM)                !VELOCITY
      real(8)                  :: X1(dim),X2(dim)  !the massweighted coordinates
      real(8)                  :: x1m(dim),x2m(dim),x1p(dim),x2p(dim)
      real(8)                  :: RMASS(dim/3)
      Real(8)                  :: RTS_(dim) !FOR TS ESTIMATION
      Real(8)                  :: RTSOLD(dim) !FOR TS ESTIMATION (THIS IS DIRTY BECAUSE IT SHOULD BE SAVED!)
      real(8)                  :: fricusedperp,fricusedpara,fricusedrot
      real(8)                  :: frot0(dim)  !the rotational part of the forces
      real(8)                  :: fperp0(dim) !the perp. part of the forces
      real(8)                  :: fpara0(dim) !the perp. part of the forces
      real(8)                  :: SVARV(DIM),SVAR
      integer(4)               :: i,j,IAT
      logical(4)               :: perpfirst,parafirst,rotfirst
      real(8)                  :: fp1,fp2
      real(8)                  :: s1,s2,a,b,c,d_
      !******************************************************************

      CELLKIN1=0.0d0
      CELLKIN2=0.0d0
      R1P=0.0d0
      R2P=0.0d0
      IF(.not.ALLOCATED(RC)) allocate(RC(dim))
      IF(.not.ALLOCATED(RCOLD)) allocate(RCOLD(dim))
      IF(.not.ALLOCATED(F1W)) allocate(F1W(dim))
      
      !========================================
      !== init the massweighted coordinates  ==  
      !========================================
      call DIMER$GET_massweighted(dim,R1,X1)
      call DIMER$GET_massweighted(dim,R2,X2)
!REMEMBER: SQDIMERDIST IS THE *ACTUAL* DISTANCE**2, SQD=D*D IS THE DESIRED DISTANCE**2
      SQDIMERDIST=dot_product((X1-X2),(X1-X2)) 
      write(dprotfil,*)"DIMER PROPAGATE: DIMERDISTANCE= ",sqrt(sqdimerdist) &
     &                ,"  UNMW= ",sqrt(dot_product(R1-R2,R1-R2))


      !===============================================
      !==========    PLACE THE DIMER  ================
      !===============================================
      !THIS SHOULD ONLY BE DONE IN THE FIRST ITERATION STEP
      if(placedimer.and..not.KDLENGTH) then !in the first step  
         CALL PLACE_DIMER(D,X1,X2)!place it so that length-constrain for the massweighted coord.is fulfilled
         !the coorect start positions (they will be written to R0 and RM in paw_atoms.f90)
         call DIMER$GET_unmassweighted(dim,X1,R1P)
         call DIMER$GET_unmassweighted(dim,X2,R2P)

         !FOR OUTPUT========================================
         !init the position of the dimer center
         RC=R1P+(R2P-R1P)/2.0d0
         RCOLD=RC
         SQDIMERDIST=dot_product((x1-x2),(x1-x2))
         write(dprotfil,*)"DIMER: DIMER_PROPAGATE: WE PLACED THE DIMER AND NOW HAVE A DISTANCE OF",SQRT(SQDIMERDIST)
         !FOR OUTPUT========================================

         return !do not propagate!!!
      end if

      !===============================================
      !======    KEEP THE DIMERS LENGTH   ============
      !===============================================
      if(KDLENGTH) then !we keep the distance 
         D=sqrt(SQDIMERDIST)   !d is the constained distance, sqd the squared constrained
         SQD=SQDIMERDIST       !and sqdimerdist the actual squared distance 

         R1P=R1
         R2P=R2


         !FOR OUTPUT========================================
         sqdimerdist=dot_product((x1-x2),(x1-x2))
         write(dprotfil,*)"DIMER: DIMER_PROPAGATE: WE KEPT THE DIMERLENGTH AND NOW HAVE A DISTANCE OF",SQRT(SQDIMERDIST)
         !FOR OUTPUT========================================


         !this should be done in paw_atoms
         !KDLENGTH=.false. !we don't need to place the dimer any more
         !placedimer=.false.
         return !do not propagate!!!
      end if


      !===============================================
      !======        PCLIMP-TESTING       ============
      !===============================================
      if(CLIMBPERP) then
         call dimer_pclimb(R1,R2,F1,F2,R1P,R2P)
         !do nothing else!
         return
      end if



      !==============================================
      !======    THIS IS THE NEW EQM         ========
      !==============================================
         call DIMER$GET_massweighted(dim,R1_,X1m)
         call DIMER$GET_massweighted(dim,R2_,X2m)

         call DIMER$GET_unmassweighted(dim,f1,Uf1)
         call DIMER$GET_unmassweighted(dim,f2,Uf2)

         call DIMER$GET_massweighted(dim,f1,mf1)
         call DIMER$GET_massweighted(dim,f2,mf2)

         if(.true..or.optfricperp.or.optfricpara.or.optfricrot) then
            perpfirst=.false.
            parafirst=.false.
            rotfirst=.false.
            if(.not.allocated(fperpm)) then
               allocate(fperpm(dim))
               perpfirst=.true.
            end if
            if(.not.allocated(fparam)) then
               allocate(fparam(dim))
               parafirst=.true.
            end if
            if(.not.allocated(frotm)) then
               allocate(frotm(dim))
               rotfirst=.true.
            end if

            !======    get the optimal friction    ========
            call dimer$optanner(dim,dt,mf1,mf2,x1,X1m,x2,x2m,fmperp,fmpara,fmrot,&
                 &fperpm,fparam,frotm,fperp0,fpara0,frot0,fricusedperp,fricusedpara,fricusedrot)


            fperpm(:)=fperp0(:)
            fparam(:)=fpara0(:)
            frotm(:)=frot0(:)


            !do we use optfric? 
            !use the optfric not in the first step (frotm not initialized)
            !use the userdefined fixed friction instead
            if((.not.optfricperp).or.perpfirst) fricusedperp=fricperp
            if((.not.optfricpara).or.parafirst) fricusedpara=fricpara
            if((.not.optfricrot).or.rotfirst) fricusedrot=fricrot
         else
            fricusedperp=fricperp
            fricusedpara=fricpara
            fricusedrot=fricrot            
         end if

         !UF, because subr prop needs the forces in the form f=m^-(1/2)f!
         !call dimer_prop(dim,x1p,x1,x1m,Uf1,x2p,x2,x2m,Uf2,dt,SQD &
         !     & ,fricusedpara,fricusedperp,fricusedrot,fmpara,fmperp,fmrot)

         call dimer_prop06new(dim,x1p,x1,x1m,Uf1,x2p,x2,x2m,Uf2,dt,SQD &
              & ,fricusedpara,fricusedperp,fricusedrot,fmpara,fmperp,fmrot,g1,g2)


      call DIMER$GET_unmassweighted(dim,X1P,R1NC)
      call DIMER$GET_unmassweighted(dim,X2P,R2NC)

      !FOR OUTPUT========================================
      !=========================================================
      !==========  Forces in Dimer Direction         ===========
      !=========================================================
      !calc unit vector in dimer direction

        SVAR=dot_product((R1-R2),(R1-R2))
        SVARV=(R1-R2)/SVAR

        !project the forces on that direction
        fp1=dot_product(SVARV,F1)
        fp2=dot_product(SVARV,F2)
        write(dprotfil,*)"DIMER FORCES: I1: ",fp1," I2: ",fp2


        write(dprotfil,*)"****************************************"
        write(dprotfil,*)"*               DIMER                  *"
        write(dprotfil,*)"*     Estimate for TS Coordinates      *"
        write(dprotfil,*)"****************************************"
        !alexdebug        RTSOLD=RTS
        !--- estimate the saddlepoint by the null of the parallel forces
        !    for fp2>fp1 we have no negative slope and the TS predictes the
        !    min!
        if(fp2*fp1.gt.0.d0) then 
           write(dprotfil,*)"WARNING: POSITIVE SLOPE!"
        end if

        RTS_=R1-SVARV*(-fp1*SVAR/(fp2-fp1))
        do i=1,(dim-2),3
           write(dprotfil,FMT='(F15.5,A2,F15.5,A2,F15.5)')RTS_(i),',',RTS_(i+1),',',RTS_(i+2)
        end do
        RTSOLD=RTSOLD-RTS_
        write(dprotfil,*)"TS ESTIMATE FROM DIMER CHANGED ABOUT ",sqrt(dot_product(RTSOLD,RTSOLD))
        write(dprotfil,*)"****************************************"



        !=========================================================
        !========  ESTIMATE TS IN UNMASSWEIGHTED SPACE   =========
        !=========================================================
        call DIMER$ESTIMATETS(DIM,R1,R2,F1,F2,etot,etot2)




      !=========================================================
      !==========  THE DIMER'S PARALLEL MOTION       ===========
      !=========================================================
 !     RC=R1NC+(R2NC-R1NC)/2.0d0
 !     RCDIFF=sqrt(dot_product((R2-R1),(R2-R1))) 
 !     RCDIFF=dot_product((R2-R1),(RC-RCOLD))/RCDIFF
 !     print*,"DIMER: DIMER_PROPAGATE: DIMER PARALLEL MOTION: ",RCDIFF
 !     RCOLD=RC
        write(dprotfil,*)'DIMER: DIMER_PROPAGATE: DIMER PARALLEL MOTION: ',&
             &dot_product(r1-r2,0.5d0*(r1nc+r2nc)-0.5d0*(r1+r2))/sqrt(dot_product(r1-r2,r1-r2))


      !=========================================================
      !==========   THE DIMER'S actual DIRECTION         =======
      !==========   AND ANGULAR CHANGE IN °              =======
      !=========================================================
        !prepare monitoring (hardwired)
        if(.not.treadam) then
           open(4242,FILE="dimerangle.out",ACTION="WRITE",POSITION='REWIND')
           do i=1,dim/3
              write(4242,*)R1NC((i-1)*3+1)
              write(4242,*)R1NC((i-1)*3+2)
              write(4242,*)R1NC((i-1)*3+3)
           end do
           do i=1,dim/3
              write(4242,*)R2NC((i-1)*3+1)
              write(4242,*)R2NC((i-1)*3+2)
              write(4242,*)R2NC((i-1)*3+3)
           end do
           close(4242)
           !open the file and read the positions
           allocate(angle1(dim))
           allocate(angle2(dim))
           allocate(angledir(dim))
           INQUIRE (FILE="dimerangle.in", EXIST = treadamnotthere)
           IF (treadamnotthere) THEN
              open(4242,FILE="dimerangle.in",ACTION="READ",POSITION='REWIND')
              do i=1,dim
                 read(4242,*)angle1(i)
              end do
              do i=1,dim
                 read(4242,*)angle2(i)
              end do
              
              close(4242)
              !write(dprotfil,*)angle1
              !write(dprotfil,*)"-------------------------------------------"
              !write(dprotfil,*)angle2
              !the normalized direction from groundstate 1 to groundstate2
              angledir=(angle2(:)-angle1(:))/sqrt(dot_product(angle2(:)-angle1(:),angle2(:)-angle1(:)))
           else
              !choose any direction
              angledir(:)=0.0d0
              angledir(1)=1.0d0
           end IF
           deallocate(angle1)
           deallocate(angle2)
           treadam=.true.
        end if

        if(.not.allocated(dirm))allocate(dirm(dim))
        svarv=(r2nc-r1nc)/sqrt(dot_product(r2nc-r1nc,r2nc-r1nc))
        !write(dprotfil,*)'DIMER: DIMER_PROPAGATE: DIMER ACTUAL DIRECTION: '
        !write(dprotfil,*)SVARV

        write(dprotfil,*)'DIMER: DIMER_PROPAGATE: DIMER ANGLE CHANGE IN °: ',&
             &acos(dot_product(svarv,dirm))&
             &*180.d0/(4.0*atan(1.0d0))

       dirm=svarv

        !monitoring
        write(dprotfil,*)'DIMER: DIMER_PROPAGATE: DIMER ANGLE MONITORING °: ',&
             &acos(dot_product(svarv,angledir))&
             &*180.d0/(4.0*atan(1.0d0))


        !=========================================================
        !==========   THE DIMER'S LENGTH CONTROL       ===========
        !=========================================================
        write(dprotfil,*)"DIMER: INITIALIZED STEPFACT check ",tinitstepfact
        if(.not.tinitstepfact) then
           !init stepfact if needed:
           if(nsteps.gt.0) then
              stepfact=(dmin/D)**(real(1,kind=8)/real(nsteps,kind=8))
              write(dprotfil,*)"DIMER: INITIALIZED STEPFACT TO ",stepfact
           end if
           tinitstepfact=.true.
        end if

        RCDIFF=sqrt(dot_product(R1NC+(R2NC-R1NC)/2.0d0-(R1+(R2-R1)/2.0d0),&
             &R1NC+(R2NC-R1NC)/2.0d0-(R1+(R2-R1)/2.0d0)))

        if(DLFLEX.and.RCDIFF.lt.RCDIFFMIN) lc=lc+1 !increase counter
        if(DLFLEX.and.RCDIFF.gt.RCDIFFMIN) lc=0 !reset counter

        if(DLFLEX.and.lc.ge.lcs.and.((D-DSTEP).gt.DMIN)) then !shorten the length
           if(nsteps.eq.0) then !we use dstep as stepsize
              !shorten dimer
              sqdimerdist=(D-DSTEP)**2
              D=D-DSTEP
              SQD=SQDIMERDIST
              write(dprotfil,*)"DIMER: DIMER_PROPAGATE: SHORTEND DIMER LENGTH TO ",SQRT(SQDIMERDIST)
           else
              !shorten dimer %-tual such that after nsteps we have DMIN 
              if((d-d*stepfact).le.dstep) then
                 sqdimerdist=(D*stepfact)**2
                 D=D*stepfact
                 SQD=SQDIMERDIST
                 write(dprotfil,*)"DIMER: DIMER_PROPAGATE: % SHORTEND DIMER LENGTH TO ",SQRT(SQDIMERDIST)
              end if
           end if
           !reset counter
           lc=0
        end if

      if(DLFLEX.and.lc.ge.lcs.and.((D+DSTEP).lt.DMIN)) then !lengthen the dimer
         sqdimerdist=(D+DSTEP)**2
         D=D*DSTEP
         SQD=SQDIMERDIST 
         write(dprotfil,*)"DIMER: DIMER_PROPAGATE: LENGTHEND DIMER LENGTH TO ",SQRT(SQDIMERDIST)
         !reset counter
         lc=0
      end if

      R1P=R1NC
      R2P=R2NC

      V1=(R1P-R1_)/2.0d0*DT
      V2=(R2P-R2_)/2.0d0*DT
      !get RMASS
      CALL ATOMLIST$GETR8A('MASS',0,(dim/3),RMASS)
      !CALC CELLKIN
      do IAT=1, dim/3
         do i=1,3
            do j=1,3
               CELLKIN1(i,j)=CELLKIN1(I,J)+RMASS(IAT)*V1(IAT*3-(3-i))*V1(IAT*3-(3-j))
               CELLKIN2(i,j)=CELLKIN2(I,J)+RMASS(IAT)*V2(IAT*3-(3-i))*V2(IAT*3-(3-j))
            end do
         end do
      end do


      
      return
    end SUBROUTINE DIMER_PROPAGATE


!##################################################################
SUBROUTINE DIMER$ESTIMATETS(N,R1,R2,F1,F2,etot,etot2)
!##################################################################
  use dimer_estimatets_module 
  integer(4),intent(in)          :: N   !the dimensionality
  real(8),intent(in)             :: r1(n)
  real(8),intent(in)             :: r2(n)
  real(8),intent(in)             :: f1(n)
  real(8),intent(in)             :: f2(n)
  real(8),intent(in)             :: etot
  real(8),intent(in)             :: etot2
  real(8)                        :: z10(n),z20(n)
  real(8)                        :: e0(n)
  real(8)                        :: f1ortho(n)
  real(8)                        :: f2ortho(n)
  real(8)                        :: eortho0(n)
  real(8)                        :: cpara,cperp
  real(8)                        :: fpara0,fperp0
  real(8)                        :: spara,sperp
  real(8)                        :: s1,s2,a,b,c,d_
  real(8)                        :: dist
  real(8)                        :: rvar1,rvar2
  real(8)                        :: ets
  real(8)                        :: xmax



  if(.not.allocated(z1m)) allocate(z1m(n))
  if(.not.allocated(z2m)) allocate(z2m(n))
  if(.not.allocated(em)) allocate(em(n))
  if(.not.allocated(eorthom)) allocate(eorthom(n))

  z10(:)=(r1(:)+r2(:))/2.d0
  z20(:)=r1(:)-r2(:)
  dist=sqrt(dot_product(z20(:),z20(:)))

  e0(:)=z20(:)/sqrt(dot_product(z20(:),z20(:)))
  
  
  !=====================================================================
  !==    determine the forces (directions referenced on image1)       ==
  !==    we project on em directions                                  ==
  !=====================================================================
  f1ortho(:)=f1(:)-em(:)*dot_product(em(:),f1(:))
  f2ortho(:)=f2(:)-em(:)*dot_product(em(:),f2(:))
  
  eortho0(:)=0.5d0*(f1ortho+f2ortho)/sqrt(dot_product(0.5d0*(f1ortho+f2ortho),&
       &0.5d0*(f1ortho+f2ortho)))
  
  
  fpara0=dot_product(em(:),(f1(:)+f2(:)))
  fperp0=dot_product(eorthom,0.5d0*(f1(:)+f2(:)))






  !============================================================================
  !==    estimate etot in between image1 and 2 using a 3rd order polynomial ===
  !============================================================================
  !--- the slope parallel to the dimer
  s1=-dot_product(e0(:),(f1(:)))
  s2=-dot_product(e0(:),(f2(:)))

  d_=etot
  c=s1
  a=-(2.d0*etot2+s2*dist+c*dist-2.d0*d_)/(-dist)**3
  b=(s2-c-3.d0*a*dist**2)/(-2.d0*dist)

  rvar1=-b/(3.d0*a)
  if(rvar1**2-(c/(3.d0*a)).lt.0.d0) then
     print*,'Error: value in sqrt negative!'
     rvar2=sqrt(rvar1**2-(c/(3.d0*a)))
  else
     rvar2=sqrt(rvar1**2-(c/(3.d0*a)))
  end if

  if(s1*s2.lt.0.d0) then
     !forces with opposite sign
     if(rvar1+rvar2.lt.0.d0.and.rvar1+rvar2.gt.-dist) then
        xmax=rvar1+rvar2
     else
        xmax=rvar1-rvar2
     end if
  else
     !forces with same sign
     if(s1.gt.0d0) then !both forces negative
        if((rvar1+rvar2).gt.0.d0) then
           xmax=rvar1+rvar2
        else
           xmax=rvar1-rvar2
        end if
     else !both forces positive
        if((rvar1+rvar2).lt.-dist) then
           xmax=rvar1+rvar2
        else
           xmax=rvar1-rvar2
        end if
     end if
  end if


  ets=a*xmax**3+b*xmax**2+c*xmax+d_
  if(ets.lt.etot.or.ets.lt.etot2) then
     print*,'Ets is smaller than etot or etot2!'
  end if

  print*,xmax
  print*,s1,s2
  print*,'ETOT1&2',etot,etot2
  print*,'ETS ESTIMATE: ',ets



  !=====================================================================
  !==    estimate using deltaF - seems not be be exact enough!!!      ==
  !=====================================================================

  !==    determine the distances (directions referenced on image1)    ==
  spara=dot_product(em,z10-z1m)
  sperp=sqrt(dot_product((z10-z1m)-spara*em,(z10-z1m)-spara*em))
  
  cpara=-(fpara0-fparam)/spara
  cperp=-(fperp0-fperpm)/sperp
  
  print*,fpara0,fparam,spara
  print*,fperp0,fperpm,sperp
  print*,'TS ESTIMATION: C VALUES',cpara,cperp
  
  !=== estimate F=0
  print*,'TSDISTPARA',fpara0/cpara
  print*,'TSDISTPERP',fperp0/cperp
  
  
  
  !=== keep some values for next turn
  !=== we need to project them now on e0 (which is em for next cycle)
  z1m=z10
  z2m=z20
  em=e0
  eorthom=eortho0
  !=== we need to project them now on e0 (which is em for next cycle)
  fparam=dot_product(e0(:),(f1(:)+f2(:)))
  fperpm=dot_product(eortho0,0.5d0*(f1(:)+f2(:)))
  



end SUBROUTINE DIMER$ESTIMATETS





!##################################################################
SUBROUTINE DIMER$STRETCH(NAT,R1,R0)
!##################################################################
  USE MPE_MODULE
!cpversion  USE MPE_COMM_SPACE_MODULE
  USE DIMER_MODULE
  IMPLICIT NONE
  integer(4),intent(in)    :: NAT 
  real(8),   intent(inout) :: R1(NAT*3) !inout because we set the position!
  real(8),   intent(in) :: R0(NAT*3) ! only for initialisation of rkeep
  real(8)              :: SQDIMERDIST
  integer(4)           :: NVAL,WORLD_ID,NWORLDTASKS,IAT
  LOGICAL(4)           :: LOOP
  if(placedimer) then
     call MPE$BROADCAST('~',1,R1)!both images use the position from 1st dimers strc
     return
  else
     IF(.not.ALLOCATED(RKEEP)) allocate(RKEEP(3*NAT))
     IF(.not.ALLOCATED(R2SPLIT)) allocate(R2SPLIT(3*NAT))

     call MPE$QUERY('MONOMER',NTASKS,THISTASK)
     CALL DIMER$GETPROCESSLEADER2(NVAL,WORLD_ID)
     !NVAL now equals the WORLD_ID of the 1st task in dimer2

     if(World_id.lt.NVAL) then
        !ONLY THE FIRST DIMER keeps its Position
        if(TFIRSTSPLIT) then
           Rkeep=R0
           TFIRSTSPLIT=.false.
        end if
        R1=Rkeep !keep R1 fixed
     end if

     if(thistask.eq.1) then
        !give WORLD_ID 1 the position on dimer2
        if(WORLD_ID.ne.1) then 
           !we are 1st task of dimer2
           call MPE$SEND('~',1,42,R1)
        else
           !we are WORLD_ID 1
           call MPE$RECEIVE('~',NVAL,42,R2SPLIT)
        end if
        if(WORLD_ID.eq.1) then
           !=========APPLY THE COUPLE-CONSTRAINT=========
           IF(ASSOCIATED(ELDEST_D)) THEN
              THIS_D=>ELDEST_D
              LOOP=.TRUE.
              do while (LOOP)
                 call ATOMLIST$INDEX(THIS_D%ID,IAT)
                 R2SPLIT((IAT*3)-3+1)=RKEEP((IAT*3)-3+1)
                 R2SPLIT((IAT*3)-3+2)=RKEEP((IAT*3)-3+2)
                 R2SPLIT((IAT*3)-3+3)=RKEEP((IAT*3)-3+3)
                 !FOR THE FOOT CONTROLLED LOOP
                 IF(associated(THIS_D%NEXT_D)) THEN
                    THIS_D=>THIS_D%NEXT_D
                 ELSE
                    LOOP=.FALSE.
                 END IF
              end do
           END IF
 

           sqdimerdist=dot_product((Rkeep-R2SPLIT),(Rkeep-R2SPLIT))
           print *, "DIMER STRETCH : d should be", stretchdist," and is", sqrt(sqdimerdist)
           if(sqdimerdist.ge.(stretchdist*stretchdist)) then
              STRETCH=.false. !communicate this to all other tasks
              !set the new length in dimer_module
              d=sqrt(sqdimerdist)
              sqd=sqdimerdist
           end if
        end if
     else
        !all tasks .ne.1 in dimer
     end if
     
     call MPE$BROADCAST('~',1,stretch)
     if(.not.stretch) then
        deallocate(RKEEP)
        deallocate(R2SPLIT)
     end if
     return
  end if
end SUBROUTINE DIMER$STRETCH
   

!##################################################################
SUBROUTINE PLACE_DIMER(DIST,X1_,X2_)
!##################################################################
!places the second dimer-point according to the length constrain
!moves it therefore along the dimer axis  
  USE DIMER_MODULE
  IMPLICIT NONE
  real(8),intent(in) :: dist 
  real(8),intent(in) :: X1_(dim)
  real(8),intent(inout) :: X2_(dim)
  integer(4)            :: IAT
  LOGICAL(4)            :: LOOP

  !======DO WE HAVE TO USE COUPLE-CONSTRAINTS ?========
  IF(ASSOCIATED(ELDEST_D)) THEN
     THIS_D=>ELDEST_D
     LOOP=.TRUE.
     do while (LOOP)
        call ATOMLIST$INDEX(THIS_D%ID,IAT)
        X2_((IAT*3)-3+1)=X1_((IAT*3)-3+1)
        X2_((IAT*3)-3+2)=X1_((IAT*3)-3+2)
        X2_((IAT*3)-3+3)=X1_((IAT*3)-3+3)

        !FOR THE FOOT CONTROLLED LOOP
        IF(associated(THIS_D%NEXT_D)) THEN
           THIS_D=>THIS_D%NEXT_D
        ELSE
           LOOP=.FALSE.
        END IF
     end do
  END IF
    !=============PACE THE SECOND IMAGE==============
  X2_(:)=X1_(:)+DIST*(X2_(:)-X1_(:))/sqrt(dot_product(X2_(:)-X1_(:),X2_(:)-X1_(:)))
  return
end SUBROUTINE PLACE_DIMER





!      .................................................................
       subroutine dimer_prop06new(n,x1p,x10,x1m,f1,x2p,x20,x2m,f2,dt,d2 &
      &                ,annepara,anneperp,annerot,mpara,mperp,mrot,g1_,g2_)
!      ** uses mass weighted coordinates                             **
!      **    x=sqrt(m)*r      f=(-de/dr)/sqrt(m)                     **
!         use model,only:calcmultiplieritermax,dlambda,itermax
         USE DIMER_MODULE
       implicit none
       integer(4),intent(in) :: n
       real(8)   ,intent(in) :: dt     ! time step
       real(8)   ,intent(in) :: d2     ! squared dimer length
       real(8)   ,intent(out):: x1p(n)
       real(8)   ,intent(in) :: x10(n)
       real(8)   ,intent(in) :: x1m(n)
       real(8)   ,intent(in) :: f1(n)
       real(8)   ,intent(out):: x2p(n)
       real(8)   ,intent(in) :: x20(n)
       real(8)   ,intent(in) :: x2m(n)
       real(8)   ,intent(in) :: f2(n)
       real(8)   ,intent(in) :: annepara   
       real(8)   ,intent(in) :: anneperp
       real(8)   ,intent(in) :: annerot
       real(8)   ,intent(in) :: mpara
       real(8)   ,intent(in) :: mperp
       real(8)   ,intent(in) :: mrot
       real(8)   ,intent(inout) :: g1_(n)
       real(8)   ,intent(inout) :: g2_(n)



       integer(4),parameter  :: niterx=1000 ! #(iterations)
       real(8)   ,parameter  :: tol=1.d-8
       real(8)   ,parameter  :: du=1.d-5
       real(8)               :: f1used(n),f2used(n)
       real(8)               :: fsum(n) ! f1+f2
       real(8)               :: fdiff(n) ! f1-f2
       real(8)               :: fsumpara
       real(8)               :: fdiffpara
       real(8)               :: svar1,svar2,svar3,mfrac
       real(8)               :: vec1(n),vec2(n)
       real(8)               :: uold
       real(8)               :: ap(2,2),a0(2,2),am(2,2),af(2,2),ac(2),ainv(2,2)
       real(8)               :: y1p(n),y10(n),y1m(n),y1bar(n),y1c(n)
       real(8)               :: y2p(n),y20(n),y2m(n),y2bar(n),y2c(n)
       real(8)               :: y1dot(n),y2dot(n)
       real(8)               :: lambda
       real(8)               :: det
       integer(4)            :: i,j,k
       real(8)               :: du1,du2,uold_,u0
       real(8)               :: udotold,ucenter,udotcenter
       real(8)               :: e(3,3) !the points needed for the hessian
       real(8)               :: a(2,2)
       real(8)               :: b(2)
       real(8)               :: x0(2),x0old(2),y1p_old(n),y2p_old(n)

       real(8)               :: di(2),ri(2),rip(2),xi(2),xip(2)
       real(8)               :: alphai,betai
       logical(4)            :: tconv
       real(8)               :: eold

       real(8)               :: g1old(n)
       real(8)               :: g2old(n)


       real(8)               :: xcenter(n)
       real(8)               :: xpara(n)
       real(8)               :: xperp(n)
       real(8)               :: xrot(n)
!      ****************************************************************
       tconv=.false.
       f1used(:)=f1(:)
       f2used(:)=f2(:)

       !weight down image 1 (use a positive mpara here)
       if(wdown) then
          f1used(:)=f1(:)*wdownfact
       end if
          

       y10(:)=0.5d0*(x10(:)+x20(:))
       y20(:)=x10(:)-x20(:)
       y1m(:)=0.5d0*(x1m(:)+x2m(:))
       y2m(:)=x1m(:)-x2m(:)
       fsum(:)=f1used(:)+f2used(:)
       fdiff(:)=f1used(:)-f2used(:)


       do i=1,CALCVELOCITYITERMAX
          !===================================================================
          !=== determine y1p
          !===================================================================
          svar1=1.d0/(1.d0+anneperp)
          svar2=1.d0/(1.d0+annepara)
          
          vec1(:)=(0.5d0*fsum(:)+g1_(:))*svar1*dt**2/mperp
          vec2(:)=(0.5d0*fsum(:)+g1_(:))*dt**2/(mpara*(1.d0+annepara))
          
          y1p(:)=(2.d0*svar1)*(y10(:)-y20(:)*dot_product(y20(:),y10(:))/d2)&
               &+(1.d0-anneperp)*svar1*(-y1m(:)+y20(:)*dot_product(y20(:),y1m(:))/d2)&
               &+vec1(:)-y20(:)*dot_product(y20(:),vec1(:))/d2&
               &+y20(:)*dot_product(y20(:),svar2*(2.d0*y10(:)-(1.d0-annepara)*y1m(:)))/d2&
               &+y20(:)*dot_product(y20(:),vec2(:))/d2
          


          
!print*,'y10:',y10(:)          
!print*,'y1p:',y1p(:)

          
          !===================================================================
          !=== determine y2p
          !===================================================================
          !--- y2bar
          
          svar1=1.d0+annerot
          
          y2bar(:)=y20(:)*2.d0/svar1-(1.d0-annerot)/svar1*y2m(:)&
               &+(fdiff(:)+g2_(:))*dt**2/(mrot*svar1)
          



          !---------------------------------------------------------------------
          !-----------   APPLY THE COUPLE CONSTRAINT IF NECESSARY   ------------
          !---------------------------------------------------------------------
          !do this in paw code, not in testcode
         IF(ASSOCIATED(ELDEST_D))THEN
            CALL DIMER_APPL_COUPLE_CONSTRAINT(N,Y2BAR(:))
            !Y1BAR STAYS THE SAME BECAUSE 2*(X1+X2)/2=X1+X2!
         END IF
          


          !---------------------------------------------------------------------
          !-----------   GET RID OF UNWISHED MOTION                 ------------
          !---------------------------------------------------------------------
          if(onlyrot) then
             !substract the covered distance of the center of grav.
             !the cog does not move:
             y1p(:)=y10(:)
             !alt. test subtract the parallel and the perp. distance
          end if
          
          if(inhibitup) then
             !substract the parallel part of ther covered distance of the center of grav.
             y1p(:)=y1p(:)-(y20(:)/dot_product(y20(:),y20(:)))*dot_product(y20(:),(y1p(:)-y10(:)))
          end if


          
          
          !---------------------------------------------------------------------
          !-- satisfy constraint svar1*lambda**2+svar2*lambda+svar3-0         --
          !---------------------------------------------------------------------
          mfrac=(2.d0*dt**2)/(mrot*(1.d0+annerot))
          
          svar1=dot_product(y2bar(:),y20(:))/(2.d0*mfrac*dot_product(y20(:),y20(:)))
          svar2=(dot_product(y2bar(:),y2bar(:))-d2)/(4.d0*mfrac**2*(dot_product(y20(:),y20(:))))
          svar3=svar1**2-svar2
          
          if(svar3.lt.0.d0) then
             write(dprotfil,*)'======lagrangeparameter: sqrt negative!======='
             print*,'d2',d2,dot_product(y2bar(:),y2bar(:))-d2,dot_product(y20(:),y20(:))
             print*,'dprod',dot_product(y2bar(:),y2bar(:))-d2,dot_product(y2bar(:),y2bar(:))
             print*,'mfrac',mfrac
             print*,svar1,svar1**2,svar2,svar3
             stop 'lagrangeparameter negative'

          end if
          
          ! ATTENTION IT IS IMPORTANT WHICH ROOT IS CHOSEN. 
          ! THW WRONG ROOT ROTATES THE DIMER BY ABOUT 180 DEGREE AND MESSES UP THINGS.
          ! SWITCHING FROM ONE TO THE OTHER MESSES UP THE ITERATION
          if(svar1.gt.0.d0) then
             lambda=svar1-sqrt(svar3)
          else
             lambda=svar1+sqrt(svar3)
          end if
          !should not change!!!y1p(:)=y1bar(:)+y1c(:)*lambda
          y2p(:)=y2bar(:)-2.d0*mfrac*lambda*y20(:)
          
          
          if(abs(dot_product(y2p,y2p)-d2).gt.DLAMBDA) then
             write(dprotfil,*)' constrain not fulfilled'
             write(dprotfil,*)'deviation ',dot_product(y2p,y2p),d2,dot_product(y2p,y2p)-d2
             print*,' constrain not fulfilled'
             print*,'deviation ',dot_product(y2p,y2p),d2,dot_product(y2p,y2p)-d2
             stop 'error constraint not fulfilled in paw_dimer.f90'
          end if
          
          
          
          !===================================================================
          !=== determine y1dot, y2dot
          !===================================================================
          
          y1dot(:)=(y1p(:)-y1m(:))/(2.d0*dt)
          y2dot(:)=(y2p(:)-y2m(:))/(2.d0*dt)
          
          
          !===================================================================
          !=== determine next guess for g1 and g2
          !===================================================================
          g1old(:)=g1_(:)
          g2old(:)=g2_(:)

          g1_(:)=(mperp-mpara)/d2*(y2dot(:)*dot_product(y20(:),y1dot(:))&
               &+y20(:)*dot_product(y2dot(:),y1dot(:)))
          g2_(:)=4.d0*(mpara-mperp)*y1dot(:)*dot_product(y20(:),y1dot(:))/d2

!          print*,'g1diff',g1_-g1old
!          print*,'g2diff',g2_-g2old
          svar1=0.d0
          svar2=0.d0
          do j=1,n
             svar1=svar1+abs(g1_(j)-g1old(j))
             svar2=svar2+abs(g2_(j)-g2old(j))
          end do
!          print*,'diffsum',svar1,svar2
          !if(dot_product(g1_-g1old,g1_-g1old).lt.tol&
          !     &.and.dot_product(g2_-g2old,g2_-g2old).lt.tol) exit
          if(svar1.lt.tol.and.svar2.lt.tol) then
             tconv=.true.
             exit
          end if
       end do

       !stop 'debut routine'
       !iteration ends here
       if(.not.tconv) then
          print*,'u,udot itration not converged',e(2,2)
          stop 'check it!'
       end if

       
       !-------------------------------------------
       x1p(:)=(y1p(:)+0.5d0*y2p(:))
       x2p(:)=(y1p(:)-0.5d0*y2p(:))



!!$       vec1(:)=y2p(:)/sqrt(dot_product(y2p(:),y2p(:)))
!!$       xcenter(:)=y1p(:)-y10(:)
!!$       xpara(:)=vec1(:)*dot_product(vec1(:),xcenter(:))
!!$       xperp(:)=xcenter(:)-xpara(:)
!!$       xrot(:)=x1p(:)-xcenter(:)-x10(:)
!!$
!!$       print*,'Ekin para',abs(0.5d0*mpara*dot_product(xpara(:),xpara(:)))/dt**2
!!$       print*,'Ekin perp',abs(0.5d0*mperp*dot_product(xperp(:),xperp(:)))/dt**2
!!$       print*,'Ekin rot',abs(0.5d0*mrot*dot_product(xrot(:),xrot(:)))/dt**2
!!$
!!$       if(abs(0.5d0*mpara*dot_product(xpara(:),xpara(:)))/dt**2.gt.0.5d0) then
!!$          x1p(:)=x1p(:)-(xpara(:)-xpara(:)/sqrt(dot_product(xpara(:),xpara(:)))*sqrt(2.d0*0.5d0/abs(mpara))*dt)
!!$          x2p(:)=x2p(:)-(xpara(:)-xpara(:)/sqrt(dot_product(xpara(:),xpara(:)))*sqrt(2.d0*0.5d0/abs(mpara))*dt)
!!$
!!$          vec1(:)=y2p(:)/sqrt(dot_product(y2p(:),y2p(:)))
!!$
!!$          y1p(:)=(x1p(:)+x2p(:))/2.d0
!!$          xcenter(:)=y1p(:)-y10(:)
!!$          xpara(:)=vec1(:)*dot_product(vec1(:),xcenter(:))
!!$          xperp(:)=xcenter(:)-xpara(:)
!!$          xrot(:)=x1p(:)-xcenter(:)-x10(:)
!!$          
!!$          print*,'Ekin para',abs(0.5d0*mpara*dot_product(xpara(:),xpara(:)))/dt**2
!!$          print*,'Ekin perp',abs(0.5d0*mperp*dot_product(xperp(:),xperp(:)))/dt**2
!!$          print*,'Ekin rot',abs(0.5d0*mrot*dot_product(xrot(:),xrot(:)))/dt**2
!!$       end if

       
       
       write(dprotfil,*)'DIMER: DIMER_PROPAGATE: DIMER MASSW. PARALLEL MOTION:',&
            &dot_product(y20,0.5d0*(x1p+x2p)-0.5d0*(x10+x20))/sqrt(dot_product(y20,y20))
       



       return
       
     end subroutine dimer_prop06new







!!$!      .................................................................
!!$       subroutine dimer_prop(n,x1p,x10,x1m,f1,x2p,x20,x2m,f2,dt,d2 &
!!$      &                ,annepara,anneperp,annerot,mpara,mperp,mrot)
!!$!      ** uses mass weighted coordinates                             **
!!$!      **    x=sqrt(m)*r      f=(-de/dr)/sqrt(m)                     **
!!$  USE DIMER_MODULE
!!$       implicit none
!!$       integer(4),intent(in) :: n
!!$       real(8)   ,intent(in) :: dt     ! time step
!!$       real(8)   ,intent(in) :: d2     ! squared dimer length
!!$       real(8)   ,intent(out):: x1p(n)
!!$       real(8)   ,intent(in) :: x10(n)
!!$       real(8)   ,intent(in) :: x1m(n)
!!$       real(8)   ,intent(in) :: f1(n)
!!$       real(8)   ,intent(out):: x2p(n)
!!$       real(8)   ,intent(in) :: x20(n)
!!$       real(8)   ,intent(in) :: x2m(n)
!!$       real(8)   ,intent(in) :: f2(n)
!!$       real(8)   ,intent(in) :: annepara   
!!$       real(8)   ,intent(in) :: anneperp
!!$       real(8)   ,intent(in) :: annerot
!!$       real(8)   ,intent(in) :: mpara
!!$       real(8)   ,intent(in) :: mperp
!!$       real(8)   ,intent(in) :: mrot
!!$       integer(4),parameter  :: niterx=1000 ! #(iterations)
!!$       real(8)   ,parameter  :: tol=1.d-8
!!$       real(8)               :: dkin    ! x1dot**1-x2dot**2
!!$       real(8)               :: dkinold  
!!$       real(8)               :: fsum(n) ! f1+f2
!!$       real(8)               :: fdiff(n) ! f1-f2
!!$       real(8)               :: fsumpara
!!$       real(8)               :: fdiffpara
!!$       real(8)               :: svar,svar1,svar2,svar3,mfrac
!!$       real(8)               :: u,udot  !u=y2*y1dot
!!$       real(8)               :: ap(2,2),a0(2,2),am(2,2),af(2,2),ac(2),ainv(2,2)
!!$       real(8)               :: y1p(n),y10(n),y1m(n),y1bar(n),y1c(n)
!!$       real(8)               :: y2p(n),y20(n),y2m(n),y2bar(n),y2c(n)
!!$       real(8)               :: qp,q0,qm
!!$       real(8)               :: lambda
!!$       real(8)               :: det
!!$       real(8)               :: p1,g1,p2,g2
!!$       integer(4)            :: iter
!!$       real(8)               :: R1p(n),R2p(n),R10(n),R20(n)
!!$       real(8)               :: CENTER1NC(3),CENTER2NC(3)
!!$       real(8)               :: CC1(n),CC2(n)
!!$       real(8)               :: RMASS(n/3)
!!$       integer(4)            :: IAT
!!$!      ****************************************************************
!!$!to prevent unassigned variable before 'check convergence'
!!$       p1=42
!!$       g1=42
!!$
!!$       y10(:)=x10(:)+x20(:)
!!$       y20(:)=x10(:)-x20(:)
!!$       y1m(:)=x1m(:)+x2m(:)
!!$       y2m(:)=x1m(:)-x2m(:)
!!$       fsum(:)=f1(:)+f2(:)
!!$       fdiff(:)=f1(:)-f2(:)
!!$       fsumpara=dot_product(y20(:),fsum(:))
!!$       fdiffpara=dot_product(y20(:),fdiff(:))
!!$       qm=dot_product(y20(:),y1m(:))
!!$       q0=dot_product(y20(:),y10(:))
!!$       dkin=dot_product(y10(:)-y1m(:),y20(:)-y2m(:))/(4.d0*dt**2)
!!$       dkinold=dkin
!!$
!!$!
!!$       do iter=1,CALCMULTIPLIERITERMAX
!!$!        =====================================================================
!!$!        == select dkin                                                     ==
!!$!        =====================================================================
!!$         If(iter.gt.2) then
!!$           dkin=p1-g1*(p1-p2)/(g1-g2)
!!$         else if(iter.eq.2) then
!!$           dkin=p1+g1
!!$         end if
!!$!        =======
!!$!dkin=(-1.d0+real(iter,8)/real(niterx,8))*2.d0
!!$!        =====================================================================
!!$!        == determine u,udot for a given dkin                               ==
!!$!        =====================================================================
!!$         svar=annepara/4.d0
!!$         svar1=2.d0/(1.d0+svar)
!!$         svar2=(1.d0-svar)/(1.d0+svar)
!!$         svar3=dt**2/(mpara*(1.d0+svar))
!!$         qp=svar1*q0-svar2*qm+svar3*(fsumpara-(mpara-mperp)*dkin)
!!$         u=(qp-qm)/(2.d0*dt)
!!$         udot=dkin+(qp-2.d0*q0+qm)/dt**2
!!$!write(dprotfil,*)'u,udot,dkin ',u,udot,dkin
!!$!         
!!$!        =====================================================================
!!$!        == set up linear system of equations for y1p,y2p                   ==
!!$!        ==       ap*yp=a0*y0+am*ym+af*f+ac*y2*lambda                       ==
!!$!        =====================================================================
!!$         mfrac=(mpara-mperp)/(mperp*d2)
!!$         ap(1,1)=1.d0+anneperp/2.d0
!!$         ap(1,2)=0.5d0*mfrac*(u*dt)
!!$         am(1,1)=-(1.d0-anneperp/2.d0)
!!$         am(1,2)=ap(1,2)
!!$         a0(1,1)=2.d0
!!$         a0(1,2)=-mfrac*(udot*dt**2)-(u*dt)/(2.d0*mperp*d2) &
!!$        &                              *(annepara*mpara-anneperp*mperp)
!!$         af(1,1)=dt**2/mperp
!!$         af(1,2)=0.d0
!!$         ac(1)=0.d0
!!$         mfrac=(mpara-mperp)/(mrot*d2)
!!$         ap(2,2)=1.d0+annerot/4.d0
!!$         ap(2,1)=-0.5d0*mfrac*(u*dt)
!!$         a0(2,1)=0.d0
!!$         a0(2,2)=2.d0
!!$         am(2,2)=-(1.d0-annerot/4.d0)
!!$         am(2,1)=ap(2,1)
!!$         af(2,1)=0.d0
!!$         af(2,2)=dt**2/mrot
!!$         ac(2)=4.d0*dt**2/mrot
!!$!
!!$!        =====================================================================
!!$!        == resolve equations for y1bar,y2bar, y1c y2c                      ==
!!$!        =====================================================================
!!$         det=ap(1,1)*ap(2,2)-ap(1,2)*ap(2,1)
!!$!write(dprotfil,*)'ap ',ap
!!$!write(dprotfil,*)'a0 ',a0
!!$!write(dprotfil,*)'am ',am
!!$!write(dprotfil,*)'af ',af
!!$!write(dprotfil,*)'ac ',ac
!!$!stop
!!$         ainv(1,1)= ap(2,2)/det
!!$         ainv(1,2)=-ap(1,2)/det
!!$         ainv(2,1)=-ap(2,1)/det
!!$         ainv(2,2)= ap(1,1)/det
!!$         a0(:,:)=matmul(ainv(:,:),a0(:,:))
!!$         am(:,:)=matmul(ainv(:,:),am(:,:))
!!$         af(:,:)=matmul(ainv(:,:),af(:,:))
!!$         ac(:)  =matmul(ainv(:,:),ac(:))
!!$         Y1bar(:)=a0(1,1)*y10(:) +a0(1,2)*y20(:) &
!!$      &          +am(1,1)*y1m(:) +am(1,2)*y2m(:) &
!!$      &          +af(1,1)*fsum(:)+af(1,2)*fdiff(:)
!!$         Y2bar(:)=a0(2,1)*y10(:) +a0(2,2)*y20(:) &
!!$      &          +am(2,1)*y1m(:) +am(2,2)*y2m(:) &
!!$      &          +af(2,1)*fsum(:)+af(2,2)*fdiff(:)
!!$         y1c(:)=ac(1)*y20(:)
!!$         y2c(:)=ac(2)*y20(:)
!!$!         
!!$

!!$!        =====================================================================
!!$!        == ap: project rotation & translation out                              ==
!!$!        =====================================================================
!!$!        do this with bar values -> the length constrained positions
!!$!        do not contain cellrotation/translation
!!$         x1p(:)=0.5d0*(y1bar(:)+y2bar(:))
!!$         x2p(:)=0.5d0*(y1bar(:)-y2bar(:))
!!$         call DIMER$GET_unmassweighted(dim,X1p,R1p)
!!$         call DIMER$GET_unmassweighted(dim,X2p,R2p)
!!$         call DIMER$GET_unmassweighted(dim,X10,R10)
!!$         call DIMER$GET_unmassweighted(dim,X20,R20)
!!$
!!$write(dprotfil,*)'debug43R1P',R1P
!!$write(dprotfil,*)'debug43R2P',R2P
!!$write(dprotfil,*)'debug43R10',R10
!!$write(dprotfil,*)'debug43R20',R20
!!$call DIMER$GET_unmassweighted(dim,X1m,R10)
!!$call DIMER$GET_unmassweighted(dim,X2m,R20)
!!$write(dprotfil,*)'debug43R1m',R10
!!$write(dprotfil,*)'debug43R2m',R20
!!$
!!$         
!!$         CALL ATOMLIST$GETR8A('MASS',0,(n/3),RMASS)
!!$         !|_ the same for both images
!!$         
!!$         if(CENTER_ID.eq.'COG') then
!!$            !do nothing
!!$         else
!!$            !if we use not the center of gravity for rot. center
!!$            call ATOMLIST$INDEX(CENTER_ID,IAT)
!!$            !|_ the same for both images
!!$            
!!$            !the center for the actual coords is m**(1/2)*CENTER_COORD
!!$            ! the center for NC coords have to be read
!!$            CENTER1NC(1)=R1p((IAT*3)-3+1)
!!$            CENTER1NC(2)=R1p((IAT*3)-3+2)
!!$            CENTER1NC(3)=R1p((IAT*3)-3+3)
!!$            
!!$            CENTER2NC(1)=R2p((IAT*3)-3+1)
!!$            CENTER2NC(2)=R2p((IAT*3)-3+2)
!!$            CENTER2NC(3)=R2p((IAT*3)-3+3)
!!$         end if
!!$         
!!$         !                                                       | we placed it there in the last step!
!!$         call DIMER_PROPCELLCONSTR((dim/3),R10,R1p,RMASS,RMASS,(sqrt(RMASS(IAT))*CENTER_COORD(:)),CENTER1NC,CC1)
!!$         call DIMER_PROPCELLCONSTR((dim/3),R20,R2p,RMASS,RMASS,(sqrt(RMASS(IAT))*CENTER_COORD(:)),CENTER2NC,CC2)
!!$write(dprotfil,*)'debug43cc2',CC1
!!$write(dprotfil,*)'debug43cc1',CC2
!!$
!!$         call DIMER$GET_massweighted(dim,CC1,X1p)
!!$         call DIMER$GET_massweighted(dim,CC2,X2p)
!!$
!!$         !reconstruct the ybar vectors (now without cellrotation/translation)
!!$         y1bar(:)=x1p(:)+x2p(:)
!!$         y2bar(:)=x1p(:)-x2p(:)
!!$write(dprotfil,*)'debug42y1bar',y1bar
!!$write(dprotfil,*)'debug42y2bar',y2bar
!!$write(dprotfil,*)'debug42i',iter

!!$!        =====================================================================
!!$!        ===========   APPLY THE COUPLE CONSTRAINT IF NECESSARY   ============
!!$!        =====================================================================
!!$         IF(ASSOCIATED(ELDEST_D))THEN
!!$            CALL DIMER_APPL_COUPLE_CONSTRAINT(N,Y2BAR(:))
!!$            !Y1BAR STAYS THE SAME BECAUSE 2*(X1+X2)/2=X1+X2!
!!$         END IF
!!$
!!$!        =====================================================================
!!$!        == satisfy constraint svar1*lambda**2+svar2*lambda+svar3=0         ==
!!$!        =====================================================================
!!$         svar1=dot_product(y2c(:),y2c(:))
!!$         svar2=2.d0*dot_product(y2bar(:),y2c(:))
!!$         svar3=dot_product(y2bar(:),y2bar(:))-d2
!!$         svar=svar2**2-4.d0*svar1*svar3
!!$         if(svar.lt.0.d0) then
!!$           write(dprotfil,*)'======lagrange parameter could not be found======='
!!$           write(dprotfil,*)'svar ',svar1,svar2,svar3,svar
!!$           write(dprotfil,*)'y2bar ',y2bar
!!$           write(dprotfil,*)'y2c   ',y2c
!!$           stop '======lagrange parameter could not be found=======errorstop'
!!$         end if
!!$! ATTENTION IT IS IMPORTANT WHICH ROOT IS CHOSEN. 
!!$! THW WRONG ROOT ROTATES THE DIMER BY ABOUT 180 DEGREE AND MESSES UP THINGS.
!!$! SWITCHING FROM ONE TO THE OTHER MESSES UP THE ITERATION
!!$         if(svar2.gt.0.d0) then
!!$           lambda=-(svar2-sqrt(svar))/(2.d0*svar1)
!!$         else
!!$           lambda=-(svar2+sqrt(svar))/(2.d0*svar1)
!!$         end if
!!$         y1p(:)=y1bar(:)+y1c(:)*lambda
!!$         y2p(:)=y2bar(:)+y2c(:)*lambda
!!$         if(abs(dot_product(y2p,y2p)-d2).gt.DLAMBDA) then
!!$           write(dprotfil,*)' constrain not fulfilled'
!!$           write(dprotfil,*)'deviation ',dot_product(y2p,y2p)-d2
!!$           stop 'error constraint not fulfilled in paw_dimer.f90'
!!$         end if
!!$!         
!!$!        =====================================================================
!!$!        == Check convergence                                               ==
!!$!        =====================================================================
!!$         p2=p1
!!$         g2=g1
!!$         p1=dkin
!!$         g1=dot_product(y1p(:)-y1m(:),y2p(:)-y2m(:))/(4.d0*dt**2)-dkin
!!$         write(dprotfil,*)'DIMER:USING HARDCODED TOL!!!'
!!$         write(dprotfil,*)'DIMER:USING HARDCODED TOL!!!'
!!$         write(dprotfil,*)'DIMER:USING HARDCODED TOL!!!'
!!$         write(dprotfil,*)'DIMER:USING HARDCODED TOL!!!'
!!$         write(dprotfil,*)'DIMER:USING HARDCODED TOL!!!'
!!$         if(abs(g1).lt.tol) then
!!$           x1p(:)=0.5d0*(y1p(:)+y2p(:))
!!$           x2p(:)=0.5d0*(y1p(:)-y2p(:))
!!$!        =====================================================================
!!$!        == do we use only the rotation?                                    ==
!!$!        =====================================================================
!!$           if(onlyrot) then
!!$              !substract the covered distance of the center of grav.
!!$              x1p(:)=x1p(:)-0.5d0*(y1p(:)-y10(:))
!!$              x2p(:)=x2p(:)-0.5d0*(y1p(:)-y10(:))
!!$           end if
!!$!        =====================================================================
!!$!        == do we inhibit the parallel motion?                              ==
!!$!        =====================================================================
!!$           if(inhibitup) then
!!$              !substract the parallel part of ther covered distance of the center of grav.
!!$              x1p(:)=x1p(:)-(y20(:)/dot_product(y20(:),y20(:)))*dot_product(y20(:),0.5d0*(y1p(:)-y10(:)))
!!$              x2p(:)=x2p(:)-(y20(:)/dot_product(y20(:),y20(:)))*dot_product(y20(:),0.5d0*(y1p(:)-y10(:)))
!!$!              x1p(:)=x1p(:)-(((x10-x20)/dot_product(x10-x20,x10-x20))*dot_product(x10-x20,x1p-x10))
!!$!              x2p(:)=x2p(:)-(((x10-x20)/dot_product(x10-x20,x10-x20))*dot_product(x10-x20,x2p-x20))
!!$           end if
!!$
!!$           write(dprotfil,*)'DIMER: DIMER_PROPAGATE: DIMER MASSW. PARALLEL MOTION:',&
!!$                &dot_product(y20,0.5d0*(x1p+x2p)-0.5d0*(x10+x20))/sqrt(dot_product(y20,y20))
!!$          return
!!$         end if
!!$       enddo
!!$!stop
!!$       stop 'PAW_DIMER: SUBROUTINE DIMER_PROP HARDCODED TOL NOT CONVERGED'
!!$       return
!!$     end subroutine dimer_prop
!!$


!##################################################################
SUBROUTINE DIMER$WRITEENERGYTRA(dim_,RP)
!##################################################################
  USE DIMER_MODULE
  USE MPE_MODULE
!cpversion  USE MPE_COMM_SPACE_MODULE
  IMPLICIT NONE
  integer(4),intent(in)      :: dim_ 
  real(8),intent(in)         :: RP(dim_)
  real(8)                    :: ENERGY
  real(8)                    :: dist    !the distance to the transition state
  integer(4)                 :: nfil
  real(8)                    :: etot_

  if(firsttra) then
     firsttra=.false.
     return
  end if

  call MPE$QUERY('MONOMER',NTASKS,THISTASK) 
  if(thistask.ne.1) return
 
  CALL LIB$SCALARPRODUCTR8(.FALSE.,DIM,1,(RP-RTS),1,(RP-RTS),dist)
  dist=dot_product((RP-RTS),(RP-RTS))
  dist=sqrt(dist)


  call MPE$QUERY('~',NTASKS,THISTASK) 
  if(thistask.eq.1) then !dimer image1
     write(dprotfil,*)"ENERGYTRA: ACTUAL DISTANCE IMAGE 1= ",dist
  else
     write(dprotfil,*)"ENERGYTRA: ACTUAL DISTANCE IMAGE 2= ",dist
  end if

  !--- write only in energytra spacing
  if(lasttra+energytra.gt.dist) then
    return
  else
     lasttra=dist
     !==================================================================
     !== GET FILE UNIT                                                ==
     !==================================================================
     CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
     IF(THISTASK.EQ.1) THEN
        CALL FILEHANDLER$UNIT('ETRA',NFIL)
        CALL ENERGYLIST$RETURN('TOTAL ENERGY',ENERGY)     
     ELSE
        !all .neq. task 1 in dimer
        return
     END IF
     
     write(NFIL,FMT='(F15.5,F15.5)')dist,ENERGY
     CALL LIB$FLUSHFILE(NFIL)
     
     return
  end if

  return
end SUBROUTINE DIMER$WRITEENERGYTRA




       subroutine dimer$optanner(n,dt,f1,f2,x10,x1m,x20,x2m,mperp,mpara,mrot,&
            &fperpm,fparam,frotm,fperp0,fpara0,frot0,perpanner,paraanner,rotanner)
!ATTENTION: depends on not massweighted forces!!!         
!      ** determines the optimal friction                    **
!      ** supposing an harmonic potential                    **
         use dimer_module ,only: dprotfil,ekinparam,ekinperpm,ekinrotm,Tparam,Tperpm,Trotm,&
              &tfact,steps,fautopara,fautoperp,fautorot,fricautopara,fricautoperp,fricautorot
       implicit none
       integer(4),intent(in) :: n                                  ! 3*nat
       real(8)   ,intent(in) :: dt                                 ! timestep
       real(8)   ,intent(in) :: f1(n),f2(n)                        ! the actual forces
       real(8)   ,intent(in) :: X10(n),X1m(n),x20(n),X2m(n)        ! squared dimer length
       real(8)   ,intent(in) :: mperp,mpara,mrot                   ! fict. mass
       real(8)   ,intent(in) :: fperpm(n),frotm(n),fparam(n)       ! the old forces
       real(8)   ,intent(out):: fperp0(n),frot0(n),fpara0(n)       ! the forces to be stored outside
       real(8)   ,intent(out):: perpanner,paraanner,rotanner       ! opt. friction 
       real(8)               :: f1ortho(n),f2ortho(n) ! the orthogonal part of the forces
       real(8)               :: e0(n),em(n) !the unity vectors in dimer direction
       real(8)               :: y10(n),y20(n),y1m(n),y2m(n)            
       real(8)               :: svar,svar1 
       real(8)               :: xpara(n),xperp(n),xcenter(n),xrot(n)   ! the covered distance
       real(8)               :: spara,sperp,srot   ! the scalar covered distance

       real(8)               :: epara(n),eperp(n),erot(n)
       real(8)               :: emwpara(n),emwperp(n),emwrot(n)
       real(8)               :: omega2para,omega2perp,omega2rot
       real(8)               :: vec1(n)
       real(8)               :: ekinpara,ekinperp,ekinrot
       real(8)               :: cpara,cperp,crot
       real(8)               :: tpara,tperp,trot
       real(8)               :: tfperp,tfrot
!      ****************************************************************
       y10(:)=(x10(:)+x20(:))/2.d0
       y20(:)=x10(:)-x20(:)
       y2m(:)=x1m(:)-x2m(:)
       y1m(:)=(x1m(:)+x2m(:))/2.d0

 

       e0(:)=y20(:)/sqrt(dot_product(y20(:),y20(:)))
       em(:)=y2m(:)/sqrt(dot_product(y2m(:),y2m(:)))

       f1ortho(:)=f1(:)-e0(:)*dot_product(e0(:),f1(:))
       f2ortho(:)=f2(:)-e0(:)*dot_product(e0(:),f2(:))


!      =====================================================================
!      ==    determine the forces (directions referenced on image1)       ==
!      =====================================================================
       fpara0(:)=e0(:)*dot_product(e0(:),(0.5d0*(f1(:)+f2(:))))
       fperp0(:)=0.5d0*(f1ortho(:)+f2ortho(:))


!!$       if(dot_product(f1ortho(:),f2ortho(:)).lt.0.0d0) then
!!$          !antiparallel
!!$          if(dot_product(f1ortho(:),f1ortho(:)).lt.dot_product(f2ortho(:),f2ortho(:))) then
!!$             frot0(:)=f1ortho(:)
!!$          else if(dot_product(f1ortho(:),f1ortho(:)).gt.dot_product(f2ortho(:),f2ortho(:))) then
!!$             frot0(:)=-f2ortho(:)!- because we reference all on image one!!!
!!$          else
!!$             frot0(:)=0.0d0
!!$          end if
!!$       else
          !parallel
          frot0(:)=0.5d0*(f1ortho(:)-f2ortho(:))
!        end if

       write(dprotfil,*)'perpforces',sqrt(dot_product(fperp0(:),fperp0(:)))
       write(dprotfil,*)'paraforces',sqrt(dot_product(fpara0(:),fpara0(:)))
       write(dprotfil,*)'rotforces',sqrt(dot_product(frot0(:),frot0(:)))



!      =====================================================================
!      ==    determine the forces (directions referenced on image1)       ==
!      =====================================================================
       fpara0(:)=e0(:)*dot_product(e0(:),(0.5d0*(f1(:)+f2(:))))
       fperp0(:)=0.5d0*(f1ortho(:)+f2ortho(:))


       if(dot_product(f1ortho(:),f2ortho(:)).lt.0.0d0) then
          !antiparallel
          if(dot_product(f1ortho(:),f1ortho(:)).lt.dot_product(f2ortho(:),f2ortho(:))) then
             frot0(:)=f1ortho(:)
          else if(dot_product(f1ortho(:),f1ortho(:)).gt.dot_product(f2ortho(:),f2ortho(:))) then
             frot0(:)=-f2ortho(:)!- because we reference all on image one!!!
          else
             frot0(:)=f1ortho(:)
          end if
       else
          !parallel
          frot0(:)=0.5d0*(f1ortho(:)-f2ortho(:))
       end if

       write(dprotfil,*)'perpforces',sqrt(dot_product(fperp0(:),fperp0(:)))
       write(dprotfil,*)'paraforces',sqrt(dot_product(fpara0(:),fpara0(:)))
       write(dprotfil,*)'rotforces',sqrt(dot_product(frot0(:),frot0(:)))

!      =====================================================================
!      ==determine the covered distances (directions referenced on image1)==
!      =====================================================================
       xcenter(:)=y10(:)-y1m(:)
       xpara(:)=em(:)*dot_product(em(:),xcenter(:))
       xperp(:)=xcenter(:)-xpara(:)
       xrot(:)=2.d0*(x10(:)-xcenter(:)-x1m(:))


       spara=sqrt(dot_product(xpara(:),xpara(:)))
       if(spara.lt.1.d-10) then
          epara(:)=0.d0
       else
          epara(:)=xpara(:)/spara
       end if

       sperp=sqrt(dot_product(xperp(:),xperp(:)))
       if(sperp.lt.1.d-10) then
          eperp(:)=0.d0
       else
          eperp(:)=xperp(:)/sperp
       end if

       srot=sqrt(dot_product(xrot(:),xrot(:)))
       if(srot.lt.1.d-10) then
          erot(:)=0.d0
       else
          erot(:)=xrot(:)/srot
       end if


       write(dprotfil,*)'paradist massweighted!',spara
       write(dprotfil,*)'perpdist massweighted!',sperp
       write(dprotfil,*)'rotdist massweighted!',srot


!      =====================================================================
!      ==               determine the optimal friction                    ==
!      =====================================================================

       if(spara.lt.1.d-10) then
          cpara=0.d0
          omega2para=0.d0
          paraanner=0.d0
       else
          cpara=-(dot_product(fpara0(:),epara(:))-dot_product(fparam(:),epara(:)))/&
               &spara
          write(dprotfil,*)'cdirectpara',cpara
          call dimer_proposcillator(dt,1,cpara) !cpara oscillator uses the dim 1
          write(dprotfil,*)'cmeanpara',cpara
          call DIMER$GET_massweighted(n,epara(:),emwpara(:))
          omega2para=abs(cpara/(mpara*dot_product(emwpara(:),emwpara(:))))
          paraanner=0.5d0*dt*sqrt(4.d0*omega2para)
          if(paraanner.gt.1.d0)paraanner=1.d0
       end if


       if(sperp.lt.1.d-10) then
          cperp=0.d0
          omega2perp=0.d0
          perpanner=0.d0
       else
          cperp=-(dot_product(fperp0(:),eperp(:))-dot_product(fperpm(:),eperp(:)))/&
               &sperp
          write(dprotfil,*)'cdirectperp',cperp
          call dimer_proposcillator(dt,2,cperp) !cperp oscillator uses the dim 2
          write(dprotfil,*)'cmeanperp',cperp
          call DIMER$GET_massweighted(n,eperp(:),emwperp(:))
          omega2perp=abs(cperp/(mperp*dot_product(emwperp(:),emwperp(:))))
          perpanner=0.5d0*dt*sqrt(4.d0*omega2perp)
       end if


       if(srot.lt.1.d-10) then
          crot=0.d0
          omega2rot=0.d0
          rotanner=0.d0
       else
          crot=-(dot_product(frot0(:),erot(:))-dot_product(frotm(:),erot(:)))/&
               &srot
          write(dprotfil,*)'cdirectrot',crot
          call dimer_proposcillator(dt,3,crot) !crot oscillator uses the dim 3
          write(dprotfil,*)'cmeanrot',crot
          call DIMER$GET_massweighted(n,erot(:),emwrot(:))
          omega2rot=abs(crot/(dot_product(emwrot(:),emwrot(:))))

          omega2rot=(omega2rot-omega2para)/mrot
          if(omega2rot.lt.0.d0) then
             print*,'omega2rot is negative - check that'
             print*,omega2rot
             print*,omega2perp,omega2para,abs(omega2para*2.d0*mpara)
             omega2rot=abs(omega2rot)
          end if
          rotanner=0.5d0*dt*sqrt(4.d0*omega2rot)
       end if
       

!      =====================================================================
!      ==             determine the friction for max ekin                 ==
!      =====================================================================
       Ekinpara=abs(0.5d0*mpara*dot_product(xpara(:),xpara(:)))/dt**2
       Ekinperp=abs(0.5d0*mperp*dot_product(xperp(:),xperp(:)))/dt**2
       Ekinrot=abs(0.5d0*mrot*dot_product(xrot(:),xrot(:)))/dt**2
       tfperp=tfact/real(n-1,kind=8)
       tfrot=tfact/real(n,kind=8)
       Tpara=ekinpara*tfact
       Tperp=ekinperp*tfperp
       Trot=ekinrot*tfrot


       print*,'Temperature',tpara,tperp,trot
       print*,'Ekin:',ekinpara,ekinperp,ekinrot
       print*,'maxt:',tparam,tperpm,trotm

       if(tpara.gt.tparam) then
          !          paraanner=0.5d0*dt*(ekinpara-ekinparam-dot_product(0.5d0*fpara0(:),xpara(:)))&
          !               &*dt/dot_product(xpara(:),xpara(:))
          !new testcode:we have optfric from above and add \Delta E
          svar=abs(0.5d0*dt*dt*(ekinpara-tparam/tfact)/(dot_product(xpara(:),xpara(:))*mpara))
          paraanner=paraanner+svar
          print*,'corrected parallel friction: (value/delta)',paraanner,svar
       end if

       if(tperp.gt.tperpm) then
          !perpanner=0.5d0*dt*(ekinperp-ekinperpm-dot_product(0.5d0*fperp0(:),xperp(:)))&
          !     &*dt/dot_product(xperp(:),xperp(:))
          svar=abs(0.5d0*dt*dt*(ekinperp-tperpm/tfperp)/(dot_product(xperp(:),xperp(:))*mperp))
          perpanner=perpanner+svar
          print*,'corrected perp friction: (value/delta)',perpanner,svar
       end if

       if(trot.gt.trotm) then
          !rotanner=0.5d0*dt*(ekinrot-ekinrotm-dot_product(0.5d0*frot0(:),xrot(:)))&
          !     &*dt/dot_product(xrot(:),xrot(:))
          svar=abs(0.5d0*dt*dt*(ekinrot-trotm/tfrot)/(dot_product(xrot(:),xrot(:))*mrot))
          rotanner=rotanner+svar
          print*,'corrected rot friction: (value/delta)',rotanner,svar
       end if


!      =====================================================================
!      ==                  USE FORCE AUTOPILOT                            ==
!      =====================================================================
       !calc 10 steps (do not quench start-up wriggles)
       if(steps.gt.10) then
          if(fautopara) then
             print*,'dirpara',dot_product(xpara,fpara0)
             if((dot_product(xpara,fpara0).gt.0.d0).and.(paraanner.lt.fricautopara)) then
                paraanner=fricautopara
                print*,'force autopilot corrected para friction'
             end if
          end if

          if(fautoperp) then
             print*,'dirperp',dot_product(xperp,fperp0)
             if((dot_product(xperp,fperp0).lt.0.d0).and.(perpanner.lt.fricautoperp)) then
                perpanner=fricautoperp
                print*,'force autopilot corrected perp friction'
             end if
          end if

          if(fautorot) then
             print*,'dirrot',dot_product(xrot,frot0)
             if((dot_product(xrot,frot0).lt.0.d0).and.(rotanner.lt.fricautorot)) then
                rotanner=fricautorot
                print*,'force autopilot corrected rot friction'
             end if
          end if
       else
          steps=steps+1
       end if




!      =====================================================================
!      ==                     estimate TS                                 ==
!      =====================================================================
       call dimer$get_unmassweighted(n,y10,vec1)
       write(dprotfil,*)'CENTERCOORD,unmassweigted ',y10(:)
       write(dprotfil,*)'CENTERCOORD,unmassweigted ',vec1(:)

       !TO DO: ADD HERE A CORRECTION FOR A ROTATIONAL FORCE!

       vec1(:)=0.d0
       if(abs(cpara).gt.1.d-5) vec1(:)=vec1(:)+fpara0(:)/cpara
       if(abs(cperp).gt.1.d-5) vec1(:)=vec1(:)+fperp0(:)/cperp
       write(dprotfil,*)'TS estimate',y10(:)+vec1(:)
       write(dprotfil,*)'TS dist',sqrt(dot_product(vec1(:),vec1(:)))

       return
     end subroutine dimer$optanner






!      .................................................................
       subroutine dimer$optanner_old(n,dt,f1,f2,x10,x1m,x20,x2m,mperp,mpara,mrot,&
            &fperpm,fparam,frotm,fperp0,fpara0,frot0,perpanner,paraanner,rotanner)
!      ** determines the optimal friction                    **
!      ** supposing an harmonic potential                    **
         use dimer_module ,only: dprotfil
       implicit none
       integer(4),intent(in) :: n                                  ! 3*nat
       real(8)   ,intent(in) :: dt                                 ! timestep
       real(8)   ,intent(in) :: f1(n),f2(n)                        ! the actual forces
       real(8)   ,intent(in) :: X10(n),X1m(n),x20(n),X2m(n)        ! squared dimer length
       real(8)   ,intent(in) :: mperp,mpara,mrot                   ! fict. mass
       real(8)   ,intent(in) :: fperpm(n),frotm(n),fparam(n)       ! the old forces
       real(8)   ,intent(out):: fperp0(n),frot0(n),fpara0(n)       ! the forces to be stored outside
       real(8)   ,intent(out):: perpanner,paraanner,rotanner       ! opt. friction 
       real(8)               :: f1ortho(n),f2ortho(n) ! the orthogonal part of the forces
       real(8)               :: e0(n),em(n) !the unity vectors in dimer direction
       real(8)               :: y10(n),y20(n),y1m(n),y2m(n)            
       real(8)               :: svar,svar1 
       real(8)               :: xpara(n),xperp(n),xcenter(n),xrot(n)   ! the covered distance
!      ****************************************************************
       y10(:)=x10(:)+x20(:)
       y20(:)=x10(:)-x20(:)
       y2m(:)=x1m(:)-x2m(:)
       y1m(:)=x1m(:)+x2m(:)

       e0(:)=y20(:)/sqrt(dot_product(y20(:),y20(:)))
       em(:)=y2m(:)/sqrt(dot_product(y2m(:),y2m(:)))

       f1ortho(:)=f1(:)-e0(:)*dot_product(e0(:),f1(:))
       f2ortho(:)=f2(:)-e0(:)*dot_product(e0(:),f2(:))

!      =====================================================================
!      ==    determine the forces (directions referenced on image1)       ==
!      =====================================================================
       fpara0(:)=e0(:)*dot_product(e0(:),(f1(:)+f2(:))/0.5d0)
       fperp0(:)=0.5d0*(f1ortho(:)+f2ortho(:))

       if(dot_product(f1ortho(:),f2ortho(:)).lt.0.0d0) then
          !antiparallel
          if(dot_product(f1ortho(:),f1ortho(:)).lt.dot_product(f2ortho(:),f2ortho(:))) then
             frot0(:)=f1ortho(:)
          else if(dot_product(f1ortho(:),f1ortho(:)).gt.dot_product(f2ortho(:),f2ortho(:))) then
             frot0(:)=-f2ortho(:)!- because we reference all on image one!!!
          else
             frot0(:)=0.0d0
          end if
       else
          !parallel
          frot0(:)=0.5d0*(f1ortho(:)-f2ortho(:))
       end if

       write(dprotfil,*)'perpforces',sqrt(dot_product(fperp0(:),fperp0(:)))
       write(dprotfil,*)'paraforces',sqrt(dot_product(fpara0(:),fpara0(:)))
       write(dprotfil,*)'rotforces',sqrt(dot_product(frot0(:),frot0(:)))

write(dprotfil,*)"sum42-splitted forces",2.0d0*dot_product(fperp0(:),fperp0(:))+2.0d0*dot_product(fpara0(:),fpara0(:))+&
     &2.0d0*dot_product(frot0(:),frot0(:))
write(dprotfil,*)"sum42",dot_product(f1,f1)+dot_product(f2,f2)
!      =====================================================================
!      ==determine the covered distances (directions referenced on image1)==
!      =====================================================================
       xcenter(:)=0.5d0*y10(:)-0.5d0*y1m(:)
       xpara(:)=dot_product(em(:),xcenter(:))
       xperp(:)=xcenter(:)-xpara(:)
       xrot(:)=x10(:)-xcenter(:)-x1m(:)

       write(dprotfil,*)'perpdist massweighted!',sqrt(dot_product(xperp(:),xperp(:)))
       write(dprotfil,*)'paradist massweighted!',sqrt(dot_product(xpara(:),xpara(:)))
       write(dprotfil,*)'rotdist massweighted!',sqrt(dot_product(xrot(:),xrot(:)))

!      =====================================================================
!      ==               determine the optimal friction                    ==
!      =====================================================================

       svar=sqrt(dot_product(fperp0-fperpm,fperp0-fperpm))
       svar1=sqrt(dot_product(xperp(:),xperp(:)))
!original       perpanner=2.0d0*sqrt(mperp*svar/svar1)
       perpanner=dt*sqrt(mperp*svar/svar1)

write(dprotfil,*)"debug42perp","svar",svar,"svar1",svar1,"perpanner",perpanner,"end"

       svar=sqrt(dot_product(fpara0-fparam,fpara0-fparam))
       svar1=sqrt(dot_product(xpara(:),xpara(:)))
!original       paraanner=2.0d0*sqrt(abs(mpara)*svar/svar1)
       paraanner=dt*sqrt(abs(mpara)*svar/svar1)

write(dprotfil,*)"debug42para","svar",svar,"svar1",svar1,"perpanner",paraanner,"end"
write(dprotfil,*)svar/svar1,mpara
       svar=sqrt(dot_product(frot0-frotm,frot0-frotm))
       svar1=sqrt(dot_product(xrot(:),xrot(:)))
!original       rotanner=2.0d0*sqrt(mrot*svar/svar1)
       rotanner=dt*sqrt(mrot*svar/svar1)

write(dprotfil,*)"debug42rot","svar",svar,"svar1",svar1,"perpanner",rotanner,"end"
       write(dprotfil,*)'perpanner',perpanner
       write(dprotfil,*)'paraanner',paraanner
       write(dprotfil,*)'rotanner',rotanner

       return
     end subroutine dimer$optanner_old



!     .................................................................. 
      SUBROUTINE DIMER_APPL_COUPLE_CONSTRAINT(N,Y2)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      INTEGER(4), INTENT(IN)   :: N
      REAL(8), INTENT(INOUT)   :: Y2(N)
      INTEGER                  :: I,IAT
      LOGICAL(4)               :: LOOP
!     ******************************************************************
      write(dprotfil,*)"APPLY COUPLE CONSTRAINTS"

      THIS_D=>ELDEST_D
      LOOP=.TRUE.
      do while (LOOP)
         call ATOMLIST$INDEX(THIS_D%ID,IAT)
         DO I=1,3
            write(dprotfil,*)"SEP: ",ABS(Y2((IAT*3)-3+I))
            IF(ABS(Y2((IAT*3)-3+I)).LT.CONSTRSTEP) THEN
               !IT'S OK TO SET IT 0
               Y2((IAT*3)-3+I)=0.0D0
            ELSE
               !USE THE USERDEFINED MAX STEP
               IF(Y2((IAT*3)-3+I).LT.0.0D0) Y2((IAT*3)-3+I)=Y2((IAT*3)-3+I)+CONSTRSTEP 
               IF(Y2((IAT*3)-3+I).GT.0.0D0) Y2((IAT*3)-3+I)=Y2((IAT*3)-3+I)-CONSTRSTEP 
            END IF
            write(dprotfil,*)"SEP NEW: ",ABS(Y2((IAT*3)-3+I))

         END DO

        !FOR THE FOOT CONTROLLED LOOP
        IF(associated(THIS_D%NEXT_D)) THEN
           THIS_D=>THIS_D%NEXT_D
        ELSE
           LOOP=.FALSE.
        END IF
     end do


     RETURN
    END SUBROUTINE DIMER_APPL_COUPLE_CONSTRAINT
  












!     .................................................................. 
      SUBROUTINE DIMER$SETR8(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID_
      REAL(8)     , INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'D') THEN
        D=VAL_
      ELSE IF(ID_.EQ.'STRETCHDIST') THEN
        STRETCHDIST=VAL_
      ELSE IF(ID_.EQ.'WDOWNFACT') THEN
        WDOWNFACT=VAL_
      ELSE IF(ID_.EQ.'FMPARA') THEN
        FMPARA=VAL_
      ELSE IF(ID_.EQ.'FMPERP') THEN
        FMPERP=VAL_
      ELSE IF(ID_.EQ.'FMROT') THEN
        FMROT=VAL_
      ELSE IF(ID_.EQ.'FRICPERP') THEN
        FRICPERP=VAL_
      ELSE IF(ID_.EQ.'FRICROT') THEN
        FRICROT=VAL_
      ELSE IF(ID_.EQ.'FRICPARA') THEN
        FRICPARA=VAL_
      ELSE IF(ID_.EQ.'RCDIFFMIN') THEN
        RCDIFFMIN=VAL_
      ELSE IF(ID_.EQ.'DSTEP') THEN
        DSTEP=VAL_
      ELSE IF(ID_.EQ.'DMIN') THEN
        DMIN=VAL_
      ELSE IF(ID_.EQ.'ENERGYTRA') THEN
        ENERGYTRA=VAL_
      ELSE IF(ID_.EQ.'CONSTRSTEP') THEN
        CONSTRSTEP=VAL_
      ELSE IF(ID_.EQ.'DLAMBDA') THEN
        DLAMBDA=VAL_
      ELSE IF(ID_.EQ.'TMAXPARA') THEN
        Tparam=VAL_
      ELSE IF(ID_.EQ.'TMAXROT') THEN
        Trotm=VAL_
      ELSE IF(ID_.EQ.'TMAXPERP') THEN
        Tperpm=VAL_
      ELSE IF(ID_.EQ.'FRICAUTOPARA') THEN
        fricautopara=VAL_
      ELSE IF(ID_.EQ.'FRICAUTOPERP') THEN
        fricautoperp=VAL_
      ELSE IF(ID_.EQ.'FRICAUTOROT') THEN
        fricautorot=VAL_

!alex: this is for pclimb-testing!
      ELSE IF(ID_.EQ.'FORCEDSTEP') THEN
        FORCEDSTEP=VAL_
!alex: this is for pclimb-testing!

      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIMER$SETR8')
      END IF
      RETURN
    END SUBROUTINE DIMER$SETR8
!
!     .................................................................. 
      SUBROUTINE DImer$GETR8(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(OUT):: VAL_
!     ******************************************************************
!!$      IF(ID_.EQ.'TIMESTEP') THEN
!!$        VAL_=DELT
!!$      ELSE IF(ID_.EQ.'AMPRE') THEN
!!$        VAL_=AMPRE
!!$      ELSE IF(ID_.EQ.'FRICTION') THEN
!!$        VAL_=ANNER
!!$      ELSE
!!$        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
!!$        CALL ERROR$CHVAL('ID_',ID_)
!!$        CALL ERROR$STOP('DIMER$GETR8')
!!$      END IF
      val_=0.d0 !just to get rid of compiler warnings
      RETURN
    END SUBROUTINE DIMER$GETR8
!     .................................................................. 
      SUBROUTINE DIMER$SETV(ID_,VAL_,N_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: N_
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(IN) :: VAL_(N_)
!     ******************************************************************
      IF(ID_.EQ.'CENTERCOORD') THEN
        CENTER_COORD(:)=VAL_(:)
      ELSE IF(ID_.EQ.'RTS') THEN
         IF(ALLOCATED(RTS)) DEALLOCATE(RTS)
         ALLOCATE(RTS(N_))
         RTS=VAL_
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIMER$SETV')
      END IF
      RETURN
    END SUBROUTINE DIMER$SETV
!
!     .................................................................. 
!!$      SUBROUTINE DImer$GETV(ID_,VAL_,N_)
!!$!     ******************************************************************
!!$!     ******************************************************************
!!$      USE DIMER_MODULE
!!$      IMPLICIT NONE
!!$      INTEGER(4)  ,INTENT(IN) :: N_
!!$      CHARACTER(*),INTENT(IN) :: ID_
!!$      REAL(8)     ,INTENT(OUT):: VAL_(N_)
!     ******************************************************************
!!$      IF(ID_.EQ.'TIMESTEP') THEN
!!$        VAL_=DELT
!!$      ELSE IF(ID_.EQ.'AMPRE') THEN
!!$        VAL_=AMPRE
!!$      ELSE IF(ID_.EQ.'FRICTION') THEN
!!$        VAL_=ANNER
!!$      ELSE
!!$        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
!!$        CALL ERROR$CHVAL('ID_',ID_)
!!$        CALL ERROR$STOP('DIMER$GETR8')
!!$      END IF
!!$      RETURN
!!$    END SUBROUTINE DIMER$GETV
!     .................................................................. 
      SUBROUTINE DIMER$SETCH(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: ID_
      CHARACTER(32), INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'CENTER_ID') THEN
        CENTER_ID=VAL_
!      ELSE IF(ID_.EQ.'') THEN
!        =VAL_
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIMER$SETCH')
      END IF
      RETURN
    END SUBROUTINE DIMER$SETCH
!
!     .................................................................. 
      SUBROUTINE DImer$GETCH(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*),  INTENT(IN) :: ID_
      CHARACTER(32) ,INTENT(OUT):: VAL_
!     ******************************************************************
!!$      IF(ID_.EQ.'TIMESTEP') THEN
!!$        VAL_=DELT
!!$      ELSE IF(ID_.EQ.'AMPRE') THEN
!!$        VAL_=AMPRE
!!$      ELSE IF(ID_.EQ.'FRICTION') THEN
!!$        VAL_=ANNER
!!$      ELSE
!!$        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
!!$        CALL ERROR$CHVAL('ID_',ID_)
!!$        CALL ERROR$STOP('DIMER$GETR8')
!!$      END IF
      val_=' ' !
      RETURN
    END SUBROUTINE DIMER$GETCH
!
!     .................................................................. 
      SUBROUTINE DIMER$SETL4(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL     ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'DIMER') THEN
        DIMER=VAL_
      ELSE IF(ID_.EQ.'KDLENGTH') THEN
        KDLENGTH=VAL_
      ELSE IF(ID_.EQ.'PLACEDIMER') THEN
        PLACEDIMER=VAL_
      ELSE IF(ID_.EQ.'STRETCH') THEN
        STRETCH=VAL_
      ELSE IF(ID_.EQ.'DIMERFOLLOWDOWN') THEN
        DIMERFOLLOWDOWN=VAL_
      ELSE IF(ID_.EQ.'INHIBITUP') THEN
        INHIBITUP=VAL_
      ELSE IF(ID_.EQ.'INHIBITPERP') THEN
        INHIBITPERP=VAL_
      ELSE IF(ID_.EQ.'ONLYROT') THEN
        ONLYROT=VAL_
      ELSE IF(ID_.EQ.'ONLYPERP') THEN
        ONLYPERP=VAL_
      ELSE IF(ID_.EQ.'WDOWN') THEN
        WDOWN=VAL_
      ELSE IF(ID_.EQ.'OPTFRICPARA') THEN
        OPTFRICPARA=VAL_
      ELSE IF(ID_.EQ.'OPTFRICPERP') THEN
        OPTFRICPERP=VAL_
      ELSE IF(ID_.EQ.'OPTFRICROT') THEN
        OPTFRICROT=VAL_
      ELSE IF(ID_.EQ.'DLFLEX') THEN
        DLFLEX=VAL_
      ELSE IF(ID_.EQ.'FAUTOPARA') THEN
        FAUTOPARA=VAL_
      ELSE IF(ID_.EQ.'FAUTOPERP') THEN
        FAUTOPERP=VAL_
      ELSE IF(ID_.EQ.'FAUTOROT') THEN
        FAUTOROT=VAL_

!alex: this is for pclimb-testing!
      ELSE IF(ID_.EQ.'CLIMBPERP') THEN
       CLIMBPERP =VAL_
!alex: this is for pclimb-testing!
!      ELSE IF(ID_.EQ.'') THEN
!        =VAL_
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIMER$SETL4')
      END IF
      RETURN
    END SUBROUTINE DIMER$SETL4
!
!     .................................................................. 
      SUBROUTINE DIMER$GETL4(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL     ,INTENT(OUT):: VAL_
!     ******************************************************************
      IF(ID_.EQ.'DIMER') THEN
        VAL_=DIMER
      ELSE IF(ID_.EQ.'PLACEDIMER') THEN
        VAL_=PLACEDIMER
      ELSE IF(ID_.EQ.'STRETCH') THEN
        VAL_=STRETCH
      ELSE IF(ID_.EQ.'KDLENGTH') THEN
        VAL_=KDLENGTH
      ELSE IF(ID_.EQ.'DIMERFOLLOWDOWN') THEN
        VAL_=DIMERFOLLOWDOWN
!      ELSE IF(ID_.EQ.'MOVE') THEN
!        VAL_=TDYN
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIMER$SETL4')
      END IF
      RETURN
    END SUBROUTINE DIMER$GETL4
!     .................................................................. 
      SUBROUTINE DIMER$SETI4(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'CLACMULTIPLIERITERMAX') THEN
        CALCMULTIPLIERITERMAX=VAL_
      ELSE IF(ID_.EQ.'NSTEPS') THEN
        NSTEPS=VAL_
      ELSE IF(ID_.EQ.'LCS') THEN
        LCS=VAL_
!      ELSE IF(ID_.EQ.'STRETCH') THEN
!        STRETCH=VAL_
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('DIMER$SETL4')
      END IF
      RETURN
    END SUBROUTINE DIMER$SETI4
!
!     .................................................................. 
      SUBROUTINE DIMER$GETI4(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(OUT):: VAL_
!     ******************************************************************
!!$      IF(ID_.EQ.'STOP') THEN
!!$        VAL_=TSTOP
!!$      ELSE IF(ID_.EQ.'RANDOMIZE') THEN
!!$        VAL_=TRANDOMIZE
!!$      ELSE IF(ID_.EQ.'MOVE') THEN
!!$        VAL_=TDYN
!!$      ELSE
!!$        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
!!$        CALL ERROR$CHVAL('ID_',ID_)
!!$        CALL ERROR$STOP('DIMER$SETL4')
!!$      END IF
      val_=0 !just to get rid of compiler warnings
      RETURN
    END SUBROUTINE DIMER$GETI4
!


!##################################################################
SUBROUTINE DIMER$GETPROCESSLEADER2(NVAL,WORLD_ID)
!##################################################################
  USE MPE_MODULE
  use dimer_module, only:thistask,ntasks
  IMPLICIT NONE
  integer(4),intent(out) :: NVAL !the world_id of processleader on image2
  integer(4),intent(out) :: WORLD_ID !the world_id of the calling process
  integer(4)             :: NWORLDTASKS

         NVAL=0
         call MPE$QUERY('MONOMER',NTASKS,THISTASK)
         call MPE$QUERY('~',NWORLDTASKS,WORLD_ID)
         
         if(WORLD_ID.ne.1.and.THISTASK.eq.1) then
            !we are task 1 in dimer2
            NVAL=WORLD_ID
         end if
         
         CALL MPE$COMBINE('~','+',NVAL)
         !NVAL now equals the WORLD_ID of the 1st task in dimer2
         !on all machines [the result is known on all machines]
  return
end SUBROUTINE DIMER$GETPROCESSLEADER2







!##################################################################
SUBROUTINE DIMER$GET_massweighted(n,R,X)
!##################################################################
  IMPLICIT NONE
  integer(4),intent(in) :: N !the dimensionality
  real(8), intent(in)   :: R(N) ! the coordinate vector in real space
  real(8), intent(out)  :: X(N) ! the massweighted coordinate vector  
  real(8)               :: Rmass(n/3) ! the massvector
  integer(4)            :: i
      !get RMASS
      CALL ATOMLIST$GETR8A('MASS',0,(n/3),RMASS)
      do i=1,n/3
         X((i-1)*3+1)=sqrt(RMASS(i))*R((i-1)*3+1)
         X((i-1)*3+2)=sqrt(RMASS(i))*R((i-1)*3+2)
         X((i-1)*3+3)=sqrt(RMASS(i))*R((i-1)*3+3)
      end do
      return
    end SUBROUTINE DIMER$GET_massweighted


!##################################################################
SUBROUTINE DIMER$GET_unmassweighted(n,X,R)
!##################################################################
  IMPLICIT NONE
  integer(4),intent(in) :: N !the dimensionality
  real(8), intent(out)   :: R(N) ! the coordinate vector in real space
  real(8), intent(in)  :: X(N) ! the massweighted coordinate vector  
  real(8)               :: Rmass(n/3) ! the massvector
  integer(4)            :: i
      !get RMASS
      CALL ATOMLIST$GETR8A('MASS',0,(n/3),RMASS)
      do i=1,n/3
         R((i-1)*3+1)=(1.0d0/sqrt(RMASS(i)))*X((i-1)*3+1)
         R((i-1)*3+2)=(1.0d0/sqrt(RMASS(i)))*X((i-1)*3+2)
         R((i-1)*3+3)=(1.0d0/sqrt(RMASS(i)))*X((i-1)*3+3)
      end do
      return
    end SUBROUTINE DIMER$GET_unmassweighted


!     ..................................................................
      SUBROUTINE DIMER$INIT()
!     ******************************************************************
!     **  INITIALIZES THE DIMER DEFAULT VALUES                        **
!     **  SHOULD BE CALLED BEFORE DIMER_READIN@PAW_IOROUTINES         **
!     ******************************************************************
      USE DIMER_MODULE
      use dimer_oscillator_module
      IMPLICIT NONE
      !********************************
        dimerfollowdown       =  .false.
        PLACEDIMER            =  .TRUE.
        D                     =   1.5d0   ! dimerdistance
        CALCMULTIPLIERITERMAX = 1000       ! EMERGENCY EXIT FOR THE ITERATION LOOP AFTER X STEPS
        DLAMBDA               =   1.0D-10 ! EXACTNESS OF LAMBDA (= LAGRANGE MULTIPLIER) CALCULATION
        CALCVELOCITYITERMAX   = 1000       ! EMERGENCY EXIT FOR THE ITERATION LOOP AFTER X STEPS
!!        DVELOCITY             =   1.0D-10 ! EXACTNESS OF VELOCITY iteration 
        lcs                   =   1       !how many steps with no ch.of dimercent. bef. shorten sqd
        RCDIFFMIN             =   0.0001  !difference in dimer center position below lc=lc+1 
        DSTEP                 =   0.0050   !d=d-dstep
        NSTEPS                =   0       !if 0 we use dstep
        DMIN                  =   0.4d0   !the minimal dimer length (no further reducement of d)
        DLFLEX                =  .false.  !flexible dimer length
        KDLENGTH              =  .false.   !keep the length of the startup
        INHIBITUP             =  .false.  !inhibit the upward motion of the dimer
        INHIBITPERP           =  .false.
        ONLYPERP              =  .false.
        ONLYROT               =  .false.
        WDOWN                 =  .false.  !weight downdimer image 1 (F1=WDOWNFACT*F1)
        WDOWNFACT             =  10.0d0        
        ENERGYTRA             =   0.5       ! write energy/distance to RTS all n steps 0 -> never
        CENTER_ID             = 'COG'  ! center of gravity is default
!!        DROT                  =   0.001  
        STRETCH               = .false.
        STRETCHDIST           =   0.1d0
        FMPARA                =  -1.0d0
        FMPERP                =   1.0d0
        FMROT                 =   1.0d0
        fricpara              =   0.1d0
        fricperp              =   0.1d0
        fricrot               =   0.1d0
        optfricpara           = .false.
        optfricperp           = .false.
        optfricrot            = .false.  
        constrstep            = 0.01d0
        tparam                = 500.d0
        tperpm                = 500.d0
        trotm                 = 500.d0
        
        FAUTOPARA             =.false.
        FAUTOPERP             =.false.
        FAUTOROT              =.false.
        FRICAUTOPARA          = 0.05d0
        FRICAUTOPERP          = 0.05d0
        FRICAUTOROT           = 0.05d0

!alex: this is for pclimb-testing!
        CLIMBPERP             = .false.
        FORCEDSTEP            = 0.001d0
        TFIRSTCLIMBSTEP       = .true.
!alex: this is for pclimb-testing!


        treadam=.false.
        treadamnotthere=.false.


        !connect the protocoll files etc.
        call dimer_init_files()


        !init the oscillator array
        allocate(oscm(odim))
        allocate(osc0(odim))
        allocate(oscp(odim))
        allocate(oscmass(odim))  
        allocate(oscanner(odim))
        allocate(oscc(odim))

        oscp(:)=0.d0
        oscm(:)=0.d0
        osc0(:)=0.d0
        
        !actually set by hand
        oscmass(:)=0.1d0*1822.88d0   !100. a.u.
        oscc(:)=1.d0 
        oscanner(:)=0.1d0 !we use opt friction in the subroutine (we do not know dt now) 

        !the code for opt. friction:
        !do i=1,odim
           !we use optimized friction here
        !   oscanner(i)=dt**2*sqrt(oscc(i)/oscmass(i))
        !end do

        return
      end SUBROUTINE DIMER$INIT



!     ..................................................................
      SUBROUTINE DIMER$INITDIM()
!     ******************************************************************
!     **  INITIALIZES THE DIMER DIM Variable and the data structures  **
!     **  depending thereof. should be called after strcin            **
!     ******************************************************************
        USE DIMER_MODULE
        IMPLICIT NONE
        integer(4)                  :: nat_
        !********************************

        call ATOMLIST$NATOM(NAT_)
        dim=3*nat_
        allocate(rts(dim))
        allocate(g1(dim))
        allocate(g2(dim))
        
        g1(:)=0.d0
        g2(:)=0.d0

        return
      end SUBROUTINE DIMER$INITDIM


!     ..................................................................
      SUBROUTINE DIMER$REPORT_SETTINGS(NFIL)
!     ******************************************************************
!     **  REPORTS THE SETTINGS TO STDOUT                              **
!     **  SHOULD BE CALLED AFTER READIN_DIMER@PAW_IOROUTINES.F90      **
!     ******************************************************************
      USE DIMER_MODULE
      IMPLICIT NONE
      INTEGER(4),intent(in)       :: NFIL  
      LOGICAL(4)                  :: LOOP
      !********************************
      
      write(nfil,FMT="(A)")"========================================================"
      write(nfil,FMT="(A)")"=                                                      ="
      write(nfil,FMT="(A)")"=            SETTINGS FROM !CONTROL!DIMER              ="
      write(nfil,FMT="(A)")"=                                                      ="
      write(nfil,FMT="(A)")"========================================================"
      write(nfil,*)"dimer ",dimer
      write(nfil,*)"dimerfollowdown",dimerfollowdown      
      write(nfil,*)"D ",D                      
      write(nfil,*) "CALCMULTIPLIERITERMAX ",CALCMULTIPLIERITERMAX  
      write(nfil,*) "DLAMBDA ",DLAMBDA                
      write(nfil,*) "CALCVELOCITYITERMAX ",CALCVELOCITYITERMAX    
      write(nfil,*) "DVELOCITY ",DVELOCITY              
      write(nfil,*) "lcs ",lcs                     
      write(nfil,*) "RCDIFFMIN ",RCDIFFMIN              
      write(nfil,*) "DSTEP ",DSTEP                   
      write(nfil,*) "NSTEPS ",NSTEPS                   
      write(nfil,*) "DMIN ",DMIN                   
      write(nfil,*) "DLFLEX ",DLFLEX                  
      write(nfil,*) "KDLENGTH ",KDLENGTH               
      write(nfil,*) "INHIBITUP ",INHIBITUP              
      write(nfil,*) "INHIBITPERP ",INHIBITPERP            
      write(nfil,*) "ONLYPERP ",ONLYPERP               
      write(nfil,*) "ONLYROT ",ONLYROT                 
      write(nfil,*) "WDOWN ",WDOWN                  
      write(nfil,*) "WDOWNFACT ",WDOWNFACT                     
      write(nfil,*) "ENERGYTRA ",ENERGYTRA
      write(nfil,*)"NATOMS ",dim/3
      write(nfil,*)"STRETCH",stretch
      write(nfil,*)"STRETCHDIST",stretchdist
      write(nfil,*)"FMPARA ",fmpara
      write(nfil,*)"FMPERP ",fMPERP
      write(nfil,*)"FMROT ",fmrot
      write(nfil,*)"fricpara ",fricpara
      write(nfil,*)"fricperp ",fricperp
      write(nfil,*)"fricrot ",fricrot
      write(nfil,*)"optfricpara ",optfricpara
      write(nfil,*)"optfricperp ",optfricperp
      write(nfil,*)"optfricrot ",optfricrot
      write(nfil,*)"FAUTOPARA", FAUTOPARA   
      write(nfil,*)"FAUTOPERP",  FAUTOPERP  
      write(nfil,*)"FAUTOROT",    FAUTOROT
      write(nfil,*)"FRICAUTOPARA",FRICAUTOPARA 
      write(nfil,*)"FRICAUTOPERP",FRICAUTOPERP
      write(nfil,*)"FRICAUTOROT",FRICAUTOROT

      write(nfil,*)"TMAXPARA ",tparam
      write(nfil,*)"TMAXPERP ",tperpm
      write(nfil,*)"TMAXROT ",trotm
      write(nfil,*)"CONSTRSTEP ",CONSTRSTEP
      !REPORT THE CONSTRAINTS
      IF(ASSOCIATED(ELDEST_D)) THEN
         THIS_D=>ELDEST_D
         LOOP=.TRUE.
         do while (LOOP)
            write(nfil,*)"COUPLED ATOM=",THIS_D%ID
            IF(associated(THIS_D%NEXT_D)) THEN
               THIS_D=>THIS_D%NEXT_D
            ELSE
               LOOP=.FALSE.
            END IF
         end do
      else
         write(nfil,*)'NO CONSTRAINTS'
      END IF
      
      
      !alex: this is for pclimb-testing!
      write(nfil,*)"CLIMBPERP ",CLIMBPERP
      write(nfil,*)"FORCEDSTEP ",FORCEDSTEP
      write(nfil,*)"TFIRSTCLIMBSTEP ",TFIRSTCLIMBSTEP
      !alex: this is for pclimb-testing!

      write(nfil,*)"========================================================"
    END SUBROUTINE DIMER$REPORT_SETTINGS
    

!!$MODULE DIMER_CONSTR_MODULE
!!$  implicit none
!!$  TYPE DIMER_CONSTR_TYPE
!!$     CHARACTER(32)                    :: ID        !KEYWORD
!!$     TYPE(DIMER_CONSTR_TYPE),POINTER  :: NEXT      !YOUNGER BROTHER
!!$     TYPE(DIMER_CONSTR_TYPE),POINTER  :: PREV      !ELDER BROTHER
!!$     !   TYPE(DIMER_CONSTR_TYPE),POINTER  :: ELDEST    !ELDEST BROTHER
!!$ 
!!$  END TYPE DIMER_CONSTR_TYPE
!!$
!!$  type(DIMER_CONSTR_TYPE),POINTER  :: ELDEST
!!$  type(DIMER_CONSTR_TYPE),POINTER  :: THIS
!!$  !logical(4)                       :: CSINI
!!$  !*******************************************
!!$  PUBLIC DIMER$CONSTRLIST_INIT
!!$  PUBLIC DIMER$CONSTRLIST_ADD
!!$
!!$CONTAINS
!!$
  subroutine DIMER$CONSTRLIST_INIT(ID)
    USE DIMER_MODULE
    implicit none
    CHARACTER(32),intent(in)         :: ID
    !*******************************************
    ALLOCATE(ELDEST_D)
    this_D=>eldest_D
    NULLIFY(THIS_D%NEXT_D)
    NULLIFY(THIS_D%PREV_D)
    this_D%id=id
    return
  end subroutine DIMER$CONSTRLIST_INIT



  subroutine DIMER$CONSTRLIST_ADD(ID)
    USE DIMER_MODULE
    implicit none
    character(32),intent(in)        :: ID
    !******************************************
   IF(.NOT.ASSOCIATED(ELDEST_D)) THEN
      CALL DIMER$CONSTRLIST_INIT(ID)
   ELSE
      !go to the top
      THIS_D=>ELDEST_D
      if(THIS_D%ID.eq.ID) then
         CALL ERROR$MSG('CAN NOT USE THE ID FOR DIMER CONSTRAINTS TWICE')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('DIMER$CONSTRLIST_ADD')  
      end if
      
      !descend to the youngest brother
      !thereby test if the id is already in use 
      do while (associated(THIS_D%NEXT_D))
         if(THIS_D%ID.eq.ID) then
            CALL ERROR$MSG('CAN NOT USE THE ID FOR DIMER CONSTRAINTS TWICE')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('DIMER$CONSTRLIST_ADD')  
         end if
         THIS_D=>THIS_D%NEXT_D
      end do
      
      ALLOCATE(THIS_D%NEXT_D)
      NULLIFY(THIS_D%NEXT_D%NEXT_D)
      THIS_D%NEXT_D%PREV_D=>THIS_D
      THIS_D%NEXT_D%ID=ID
   END IF
   return
 end subroutine DIMER$CONSTRLIST_ADD
!!$ 
!!$END MODULE DIMER_CONSTR_MODULE


      !===============================================
      !======        PCLIMB-TESTING       ============
      !===============================================


  subroutine DIMER_PCLIMB(R1,R2,F1,F2,R1P,R2P)
    USE DIMER_MODULE
    implicit none
      REAL(8)   ,intent(in)    :: R1(dim)               !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,intent(in)    :: R2(dim)               !CURRRENT POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,intent(in)    :: F1(dim)               !FORCE FROM POTENTIAL
      REAL(8)   ,intent(in)    :: F2(dim)               !FORCE FROM POTENTIAL

      REAL(8)   ,intent(out)   :: R1P(dim)              !NEW POSITION VECTORS OF DIMERPOINT 1 & 2
      REAL(8)   ,intent(out)   :: R2P(dim)              !NEW POSITION VECTORS OF DIMERPOINT 1 & 2

      REAL(8)                  :: SVARV(dim)
      REAL(8)                  :: FSUM(DIM) !F1+F2
    !******************************************
    FSUM=F1+F2
    if(.not.allocated(PCLIMBDIR)) allocate(PCLIMBDIR(dim))

    !***** initialize the direction ********
    if(TFIRSTCLIMBSTEP) then
       !=== get the normalized parallel direction ===
       SVARV=(R1(:)-R2(:))/sqrt(dot_product(R1-R2,R1-R2))

       !determine the direction
       PCLIMBDIR(:)=FSUM(:)-SVARV(:)*dot_product(SVARV,FSUM)

       !normalize the direction
       PCLIMBDIR(:)=PCLIMBDIR(:)/sqrt(dot_product(PCLIMBDIR,PCLIMBDIR))

       !we initialized the direction
       TFIRSTCLIMBSTEP=.false.
       print *, "PCLIMB : INITIALIZED THE FORCEDIRECTION"
       print *, "PCLIMB : THE CLIMBFORCE IS",dot_product(PCLIMBDIR,FSUM)
       R1p(:)=R1
       R2p(:)=R2
       return
    end if


    !************* climb up : -FORCEDSTEP*PCLIMBDIR *******
    R1P(:)=R1(:)-FORCEDSTEP*PCLIMBDIR(:)
    R2P(:)=R2(:)-FORCEDSTEP*PCLIMBDIR(:)
    print *, "PCLIMB : THE CLIMBFORCE IS",dot_product(PCLIMBDIR,FSUM)
   return
 end subroutine DIMER_PCLIMB







!====================================================================
! OLD CODE TO PREVENT ROTATION/TRANSLATION ACROSS THE IMAGES
!====================================================================
!##################################################################
SUBROUTINE DIMER$CELLCONSTRAINT(NAT,RCC)
!##################################################################
  USE DIMER_MODULE
  USE MPE_MODULE
  IMPLICIT NONE
  integer(4),intent(in)  :: NAT
  real(8)   ,intent(out) :: RCC(NAT*3)

  integer(4)             :: IAT
  real(8),allocatable    :: RMASS(:)
  real(8),allocatable    :: R(:),R2(:)
  integer(4)             :: NVAL,WORLD_ID,NWORLDTASKS
  REAL(8)                :: RCENTER(3)


!     ==================================================================
!     ==   DIMER STRUCTURE-CONSTRAINT FOR TRANSLATION AND ROTATION    ==
!     ==================================================================

         !========= CORRECT TRANSLATION ===========


         allocate(RMASS(NAT))
         allocate(R(3*NAT))

         CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)         

         if(CENTER_ID.eq.'COG') then
            !use the center of gravity for rotation
            CALL ATOMLIST$GETR8A('MASS',0,NAT,RMASS)
            !correct the positions
            call PLACECENTER(NAT,R,RMASS,RCC)

print *,"centerofgravity"
         else
print *,"givenatom"
            !use a given atom for rotation
            call ATOMLIST$INDEX(CENTER_ID,IAT)
            CALL ATOMLIST$GETR8A('R(0)',IAT,3,RCENTER)
            CALL PLACEATOM(.TRUE.,NAT,R,RCENTER,RCC)
         end if

         !set the actual Positions
         R=RCC
print *,"drrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr"
print *,R
print *,"drrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr"
!DO THIS OUTSIDE
!!$         !set the corrected positions
!!$         CALL ATOMLIST$SETR8A('R(0)',0,3*NAT,RCC)
!!$         CALL ATOMLIST$SETR8A('R(-)',0,3*NAT,RCC)
!!$         !NOTE: IF YOU USE THIS CODE ELSEWHERE YOU'LL
!!$         !HAVE TO CARE ABOUT R(-) !!!


         !========= CORRECT ROTATION  ===========
         !COMMUNICATE TO DIMER1
         allocate(R2(3*NAT))

         NVAL=0
         call MPE$QUERY('MONOMER',NTASKS,THISTASK)
         call MPE$QUERY('~',NWORLDTASKS,WORLD_ID)

         if(WORLD_ID.ne.1.and.THISTASK.eq.1) then
            !we are task 1 in dimer2
            NVAL=WORLD_ID
         end if
         
         CALL MPE$COMBINE('~','+',NVAL)
         !NVAL now equals the WORLD_ID of the 1st task in dimer2

         if(WORLD_ID.eq.NVAL) then 
            !we are 1st task of dimer2
            call MPE$SEND('~',1,42,R)
         else if(WORLD_ID.eq.1) then
            !we are WORLD_ID 1
            call MPE$RECEIVE('~',NVAL,42,R2)
         end if
         !we got both vectors on the 1st task

         if(WORLD_ID.eq.1) then
            call PLACEROT(NAT,R,R2,DROT,RCC)            
         end if



         !Communicate the rot. vector RCC back to dimer2  
         if(WORLD_ID.eq.1) then 
            !we are WORLD_ID 1
            call MPE$SEND('~',NVAL,43,RCC)
         else if(WORLD_ID.eq.NVAL) then !1st in 2nd dimer r
            call MPE$RECEIVE('~',1,43,RCC)
         end if
         
         !all tasks back
         if(WORLD_ID.ge.NVAL) then 


            !only on 2nd dimer
            call MPE$BROADCAST('MONOMER',1,RCC)
            R=RCC
!!$            !set the corrected positions
!!$            CALL ATOMLIST$SETR8A('R(0)',0,3*NAT,RCC)
!!$            CALL ATOMLIST$SETR8A('R(-)',0,3*NAT,RCC)
         end if
         !NOTE: IF YOU USE THIS CODE ELSEWHERE YOU'LL
         !HAVE TO CARE ABOUT R(-) !!!


         !========= SET THE CENTER BACK  ===========
         call PLACEATOM(.false.,NAT,R,CENTER_COORD,RCC)
!!$         CALL ATOMLIST$SETR8A('R(0)',0,3*NAT,RCC)
!!$         CALL ATOMLIST$SETR8A('R(-)',0,3*NAT,RCC)
print *,"d2rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr"
print *,RCC
print *,"d2rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr"

         deallocate(RMASS)
         deallocate(R)
         deallocate(R2)  
  
  return
end SUBROUTINE DIMER$CELLCONSTRAINT




!##################################################################
SUBROUTINE DIMER_PROPCELLCONSTR(NAT,R,R2,RMASS,RMASS2,RCENTER,RCENTER2,RCC)
!##################################################################
! TRANSLATES and rotates R2 in such a way (eiter COG or given atom as center)
! that the distance \sum{(R-R2)^2} is minimized
! the center of the rotated vector is moved to center_coord from .cntl file
!
! this is used to project translation / rotation out of the *not constrained*
! coords (R1NC,R2NC) before the dimer length constraint is applied
! (R=R1, R2=R1nc; R=R2, R2=R2NC )


  USE DIMER_MODULE
  USE MPE_MODULE
  IMPLICIT NONE
  integer(4),intent(in)  :: NAT
  real(8)   ,intent(in)  :: R(NAT*3)
  real(8)   ,intent(in)  :: R2(NAT*3)
  real(8)   ,intent(in)  :: RMASS(NAT)
  real(8)   ,intent(in)  :: RMASS2(NAT)
  REAL(8)   ,intent(in)  :: RCENTER(3)
  REAL(8)   ,intent(in)  :: RCENTER2(3)
  real(8)   ,intent(out) :: RCC(NAT*3)

  real(8)                :: RCC1(3*nat),RCC2(3*NAT)





!     ==================================================================
!     ==   DIMER STRUCTURE-CONSTRAINT FOR TRANSLATION AND ROTATION    ==
!     ==================================================================

         !========= CORRECT TRANSLATION ===========
         ! set both center to the origin
 
         if(CENTER_ID.eq.'COG') then
            !use the center of gravity for rotation
            !correct the positions
            call PLACECENTER(NAT,R,RMASS,RCC1)
            call PLACECENTER(NAT,R2,RMASS,RCC2)
         else
            !use a given center for rotation
            CALL PLACEATOM(.TRUE.,NAT,R,RCENTER,RCC1)
            CALL PLACEATOM(.TRUE.,NAT,R2,RCENTER2,RCC2)
         end if


         !========= CORRECT ROTATION  ===========
            call PLACEROT(NAT,RCC1,RCC2,DROT,RCC)            

         !========= SET THE CENTER BACK  ===========
         !PLACE ONLY THE SECOND TO THE WISHED CENTER_COORD
print *, "eeeeeeeeeeeeeeeee"
print *, RCC
print *, "eeeeeeeeeeeeeeeee"
         RCC2=RCC
         call PLACEATOM(.false.,NAT,RCC2,CENTER_COORD,RCC)
  
  return
end SUBROUTINE DIMER_PROPCELLCONSTR


!##################################################################
SUBROUTINE PLACECENTER(NAT,R,RMASS,RCC)
!##################################################################
  USE DIMER_MODULE
  IMPLICIT NONE
  integer(4), intent(in) :: NAT
  real(8),intent(in)     :: R(3*NAT)
  REAL(8),intent(in)     :: RMASS(NAT)
  real(8),intent(out)    :: RCC(3*NAT)

  real(8)                :: VAL(3,NAT)
  real(8)                :: RCOM(3)
  real(8)                :: SUMM
  integer                :: i

  RCOM(:)=0.0d0
  SUMM=0.0d0

  VAL(:,:)=RESHAPE(R,(/3,NAT/))
  
  do i=1, NAT
     RCOM(:)=RCOM(:)+RMASS(i)*VAL(:,i)
     SUMM=SUMM+RMASS(i)
  end do
  
  RCOM=RCOM/SUMM

  do i=1, NAT
     VAL(:,i)=VAL(:,i)-RCOM(:)
  end do

  RCC=RESHAPE(VAL,(/3*NAT/))  
  
  return
end SUBROUTINE PLACECENTER


!##################################################################
SUBROUTINE PLACEATOM(TO_ZERO,NAT,R,RCENTER,RCC)
!##################################################################
  USE DIMER_MODULE
  IMPLICIT NONE
  logical(4), intent(in) :: TO_ZERO !T->RCENTER to zero, F-> zero -RCENTER 
  integer(4), intent(in) :: NAT
  real(8),intent(in)     :: R(3*NAT)
  REAL(8),intent(in)     :: RCENTER(3)
  real(8),intent(out)    :: RCC(3*NAT)

  real(8)                :: VAL(3,NAT)
  real(8)                :: RCOM(3)
  real(8)                :: SUMM
  integer                :: i


  VAL(:,:)=RESHAPE(R,(/3,NAT/))

  if(to_ZERO) then
     do i=1, NAT
        VAL(:,i)=VAL(:,i)-RCENTER
     end do
  else
     do i=1, NAT
        VAL(:,i)=VAL(:,i)+RCENTER
     end do
  end if
  
  RCC=RESHAPE(VAL,(/3*NAT/))  
  return
end SUBROUTINE PLACEATOM

!##################################################################
SUBROUTINE PLACEROT(NAT,R1,R2,rotdif,R2C)
!##################################################################
  USE DIMER_MODULE
  IMPLICIT NONE
  integer(4),intent(in)  :: NAT
  real(8),intent(in)     :: R1(3*NAT)
  real(8),intent(in)     :: R2(3*NAT)
  real(8),intent(in)     :: rotdif
  real(8),intent(out)    :: R2C(3*NAT) !only dimer2 will be rotated



  real(8)     :: VAL1(3,NAT),VAL2(3,NAT),VAL2ROT(3)
  real(8)     :: a0
  real(8)     :: aact
  real(8)     :: amax
  real(8)     :: astep
  real(8)     :: b0
  real(8)     :: bact
  real(8)     :: bmax
  real(8)     :: bstep
  real(8)     :: c0
  real(8)     :: cact
  real(8)     :: cmax
  real(8)     :: cstep
  real(8)     :: abest
  real(8)     :: bbest
  real(8)     :: cbest
  real(8)     :: rotminv(3)
  real(8)     :: rotmin
  real(8)     :: rotminbest
  real(8)     :: U(3,3)
  real(8)     :: difact
  real(8)     :: rotminold,rotminsum
  logical(4)  :: start,first
  integer(4)  :: i,j,k,l
  real(8)     :: pi
  
  start=.true.
  first=.true.
  rotminbest=1.0d42 !something high
  rotminold=0.0d0
  pi=4.0d0*atan(1.0d0)
  difact=rotdif+42.d0
!write(dprotfil,*)'rotminbestdebug',rotminbest


 DO WHILE(start.or.(difact.gt.rotdif))
!print *,"We are still in rotations",difact,rotdif
    if(start) then
       a0=0.0d0
       amax=2.0d0*pi
       astep=(2.d0*pi)/40.0d0
       b0=0.0d0
       bmax=pi
       bstep=(2.d0*pi)/40.0d0
!       c0=0.0d0
!       cmax=pi
!       cstep=pi/20.0d0
       start=.false.
    else
       a0=abest-astep
       amax=abest+astep
       astep=astep/10.0d0

       b0=bbest-bstep
       bmax=bbest+bstep
       bstep=bstep/10.0d0

!       c0=cbest-cstep
!       cmax=cbest+cstep
!       cstep=cstep/10.0d0
    end if
!write(dprotfil,*)'rotminbestdebug',rotminbest    
!write(dprotfil,*)'R1debug',R1
!write(dprotfil,*)'R2',R2
    do j=0,int((amax-a0)/astep)
           do k=0,int((bmax-b0)/bstep)
!              do l=0,((cmax-c0)/cstep)

                 aact=a0+j*astep
                 bact=b0+k*bstep
!                 cact=c0+l*cstep
                 
                 VAL1(:,:)=RESHAPE(R1,(/3,NAT/))
                 VAL2(:,:)=RESHAPE(R2,(/3,NAT/))
                 
              U(1,1)=cos(aact)*cos(bact)
              U(1,2)=-sin(aact)
              U(1,3)=cos(aact)*sin(bact)
              U(2,1)=sin(aact)*cos(bact)
              U(2,2)=cos(aact)
              U(2,3)=sin(aact)*sin(bact)
              U(3,1)=-sin(bact)
              U(3,2)=0.0d0
              U(3,3)=cos(bact)
                 


!!$                 U(1,1)=cos(aact)*cos(bact)*cos(cact)
!!$                 U(1,2)=-cos(aact)*cos(bact)*sin(cact)-sin(aact)*cos(cact)
!!$                 U(1,3)=cos(aact)*sin(bact)
!!$                 U(2,1)=sin(aact)*cos(bact)*cos(cact)+cos(aact)*sin(cact)
!!$                 U(2,2)=-sin(aact)*cos(bact)*sin(cact)+cos(aact)*cos(cact)
!!$                 U(2,3)=sin(aact)*sin(bact)
!!$                 U(3,1)=-sin(bact)*cos(cact)
!!$                 U(3,2)=sin(bact)*sin(cact)
!!$                 U(3,3)=cos(bact)
!write(dprotfil,*)'rotminbestdebug3',rotminbest                 
                 rotminsum=0.0d0
                 do i=1,Nat
                    call LIB$MATMULR8(3,3,1,U,VAL2(:,i),VAL2rot)
                    rotminv=(VAL1(:,i)-VAL2rot)
 !                   CALL LIB$SCALARPRODUCTR8(.FALSE.,3,1,rotminv,1,rotminv,rotmin)
                    rotmin=dot_product(rotminv(:),rotminv(:))
                   rotminsum=rotminsum+rotmin
                 end do
!write(dprotfil,*)'rotminbestdebug4',rotminbest                 
!write(dprotfil,*)'rotminbest,rotminsum',rotminbest,rotminsum
                 if(rotminsum.lt.rotminbest) then
                    rotminbest=rotminsum
                    abest=aact
                    bbest=bact
!write(dprotfil,*)'abest,aact,bbest,bact',abest,aact,bbest,bact
!write(dprotfil,*)'abest,aact,bbest,bact',bbest,bact
                 end if
                 difact=abs(rotminold-rotminsum)
!write(dprotfil,*)'difact',difact,rotminsum                 
                 rotminold=rotminsum
!              end DO
              end do
        end do
!print *,"rotminbest1,abest,bbest",rotminbest,abest,bbest



                 !write the correct vector
!!$                 U(1,1)=cos(abest)*cos(bbest)*cos(cbest)
!!$                 U(1,2)=-cos(abest)*cos(bbest)*sin(cbest)-sin(abest)*cos(cbest)
!!$                 U(1,3)=cos(abest)*sin(bbest)
!!$                 U(2,1)=sin(abest)*cos(bbest)*cos(cbest)+cos(abest)*sin(cbest)
!!$                 U(2,2)=-sin(abest)*cos(bbest)*sin(cbest)+cos(abest)*cos(cbest)
!!$                 U(2,3)=sin(abest)*sin(bbest)
!!$                 U(3,1)=-sin(bbest)*cos(cbest)
!!$                 U(3,2)=sin(bbest)*sin(cbest)
!!$                 U(3,3)=cos(bbest)

                 U(1,1)=cos(abest)*cos(bbest)
                 U(1,2)=-sin(abest)
                 U(1,3)=cos(abest)*sin(bbest)
                 U(2,1)=sin(abest)*cos(bbest)
                 U(2,2)=cos(abest)
                 U(2,3)=sin(abest)*sin(bbest)
                 U(3,1)=-sin(bbest)
                 U(3,2)=0.0d0
                 U(3,3)=cos(bbest)

                 do i=1,Nat
                    call LIB$MATMULR8(3,3,1,U,VAL2(:,i),VAL2rot)
                    do j=1,3
                       R2C(i*3-3+j)=VAL2ROT(j)
                    end do
                 end do
     end DO


     return
end SUBROUTINE PLACEROT

!!$!this is hardwired and only for my use
!!$subroutine dimer$reread()
!!$  use DIMER_MODULE
!!$  implicit none
!!$  logical(4)                    :: ex
!!$  integer(4)                    :: nfil,eos
!!$  character(256)                :: string
!!$  logical(4)                    :: lvar
!!$  real(8)                       :: rvar
!!$
!!$  nfil=4242
!!$  INQUIRE (FILE="case.cntlnew", EXIST = EX)
!!$  IF ( EX ) THEN
!!$     OPEN (nfil, FILE="case.cntlnew", STATUS="OLD", &
!!$                  &ACCESS="SEQUENTIAL", ACTION="READ", POSITION="REWIND")
!!$     do 
!!$        READ (nfil, FMT="(A)" , ADVANCE="YES" ,IOSTAT=eos) string 
!!$        if(eos.eq.-1) exit
!!$        
!!$        if(trim(adjustl(string)).eq.'optfricpara') then
!!$           READ (nfil, FMT="(L1)" , ADVANCE="YES" ,IOSTAT= eos) lvar 
!!$           call dimer$setl4('OPTFRICPARA',lvar)
!!$           write(dprotfil,*)'OPTFRICPARA SET TO',lvar
!!$        else if(trim(adjustl(string)).eq.'fricpara') then
!!$           READ (nfil, FMT="(F)" , ADVANCE="YES" ,IOSTAT= eos) rvar
!!$           call dimer$setr8('FRICPARA',rvar)
!!$           write(dprotfil,*)'FRICPARA SET TO',rvar
!!$        else
!!$           write(dprotfil,*)'REREAD CONTROLFILE: IDENTIFIER NOT RECOGNIZED:',trim(adjustl(string))
!!$           !stop
!!$        end if
!!$     end do
!!$     close(nfil)
!!$     call system('rm -rf case.cntlnew')
!!$  end IF
!!$
!!$  return
!!$end subroutine dimer$reread


subroutine dimer_proposcillator(dt,i,force)
  use dimer_oscillator_module
  real(8),intent(in)                 :: dt
  integer(4),intent(in)              :: i       !which array
  real(8),intent(inout)              :: force   !extra force
  real(8)                            :: svar1,svar2

  !desription: we use the incoming value as extra force of a driven oscillator.
  !the return force value is the value of the actual force of the oscillator 
  !after propagation.
  !this is useful, because we suppose to have an incoming value which oscillates 
  !around a mean value. what we need is the mean value. the actual force of the 
  !driven oscillator (at steady state) will give us this mean value.

  if(i.gt.odim) then
     stop 'dimer oscillator array out of range'
  end if

  !we use optimized friction here
  oscanner(i)=dt**2*sqrt(oscc(i)/oscmass(i))


  svar1=(1.d0-oscanner(i))
  svar2=1.d0/(1.d0+oscanner(i))

  !propagate
  oscp(i)=(svar2*dt**2/oscmass(i))*force&
       &-(oscc(i)*dt**2/oscmass(i)-2.d0)*svar2*osc0(i)&
       &-svar1*svar2*oscm(i)


  !switch
  oscm(i)=osc0(i)
  osc0(i)=oscp(i)

  !use -1 because we search the steady state
  force=oscc(i)*osc0(i)


  return
end subroutine dimer_proposcillator
