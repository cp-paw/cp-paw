!     ..................................................................
!     ******************************************************************
!     **  WRITTEN BY JOHANNES SCHIMPL AND CLEMENS FOERST IN 2002/03   **
!     ******************************************************************
      PROGRAM STRC2XYZ
      USE LINKEDLIST_MODULE
      USE CONSTANTS_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE)               :: LL_STRC
      INTEGER(4)                  :: NFIL,I,NAT
      LOGICAL                     :: TSTRCOUT=.FALSE.
      LOGICAL                     :: TCELL=.FALSE.
      LOGICAL                     :: PAWOK=.TRUE.
      REAL(8)                     :: ANGSTROM
      CHARACTER(128)              :: ROOTNAME     ! COMMON ROOT OF THE FILENAMES
      CHARACTER(128)              :: OBJECTNAME,LINE
      LOGICAL                     :: TCHK,TCHK1
      REAL(8)                     :: RUNIT        ! LENGTH UNIT ON STRUCTURE FILE
      REAL(8)                     :: RBAS(3,3)    ! LATTICE VECTORS
      CHARACTER(32),ALLOCATABLE   :: NAME(:)      ! ATOM NAMES
      CHARACTER(32),ALLOCATABLE   :: SPECIES(:)   ! SPECIES NAME
      REAL(8),      ALLOCATABLE   :: R(:,:)       ! ATOMIC POSITIONS
      REAL(8),      ALLOCATABLE   :: Q(:)         ! CHARGES
      CHARACTER(32),ALLOCATABLE   :: MMNAME(:)      ! ATOM NAMES
      REAL(8),      ALLOCATABLE   :: MMR(:,:)       ! ATOMIC POSITIONS FOR MM PART OF QM-MM
      REAL(8),      ALLOCATABLE   :: MMQ(:)       ! ATOMIC CHARGES MM PART OF QM-MM
      REAL(8),      ALLOCATABLE   :: RSH(:,:)       ! ATOMIC POSITIONS FOR SHADOW PART OF QM-MM
      INTEGER(4),   ALLOCATABLE   :: NEIGH(:,:)   ! NEIGHBOR LIST
      LOGICAL(4)                  :: TQMMM
      REAL(8)                     :: SVAR
      INTEGER(4)                  :: IAT,IAT1,JAT,DIR(3)
      CHARACTER(32)               :: STRING,STR
      INTEGER(4)                  :: MX,MY,MZ,IX,IY,IZ,IARG,NARG

      LOGICAL(4)                  :: TFORMAT=.FALSE.
      REAL(8)                     :: CELLDISTANCE,DIST
!     ******************************************************************
!     ==========================================================================
!     == MPE$ INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                 ==
!     ==========================================================================
      CALL MPE$INIT

      CALL LINKEDLIST$NEW(LL_STRC)
!
!     ==================================================================
!     == GET FILE NAME ROOT FROM THE ARGUMENT LIST AND CONSTRUCT      ==
!     == FILE NAMES                                                   ==
!     ==================================================================
      IARG=1
      MX=1
      MY=1
      MZ=1
      ROOTNAME=''
      CELLDISTANCE=6.D0*1.889726878D0   ! DEFAULT CELL DISTANCE IN A.U.
      NARG=COMMAND_ARGUMENT_COUNT()
      DO WHILE (IARG.LE.NARG)
         CALL GET_COMMAND_ARGUMENT(IARG,LINE)
         IARG=IARG+1
         IF(LINE(1:1).EQ.'-') THEN
            IF(+LINE(2:2).EQ.'O') THEN
               TSTRCOUT=.TRUE.
            ELSE IF(+LINE(2:2).EQ.'F') THEN
               TFORMAT=.TRUE.
            ELSE IF(+LINE(2:2).EQ.'C') THEN
               TCELL=.TRUE.
               PRINT*,'  +--------------------------------------------------------+'            
               PRINT*,' /|                                                       /|'            
               PRINT*,'/--------------------------------------------------------/ |'             
               PRINT*,'| |     MULECULE_RBAS                                    | |'     
               PRINT*,'| |      OPTIMIZATION OF THE UNIT-CELL                   | |'       
               PRINT*,'| /------------------------------------------------------|-/'
               PRINT*,'|/            WRITTEN IN DEZ. 2001 BY JOHANNES SCHIMPL   |/ '
               PRINT*,'+--------------------------------------------------------+  '
            ELSE IF(+LINE(2:2).EQ.'D') THEN
               CALL GET_COMMAND_ARGUMENT(IARG,STRING)
               IARG=IARG+1
               CALL CHAR_TO_REAL(STRING,DIST,TCHK)
               IF(TCHK) THEN
                  CELLDISTANCE=DIST*1.889726878D0
                  WRITE(*,"('USING ',F10.5,' ANG AS TARGET DISTANCE')") DIST
               END IF
            ELSE IF (+LINE(2:2).EQ.'M') THEN
               CALL GET_COMMAND_ARGUMENT(IARG,STRING)
               IARG=IARG+1
               READ(STRING,'(I4)') MX
               IF (MX.LT.1) MX=1
               CALL GET_COMMAND_ARGUMENT(IARG,STRING)
               IARG=IARG+1
               READ(STRING,'(I4)') MY
               IF (MY.LT.1) MY=1
               CALL GET_COMMAND_ARGUMENT(IARG,STRING)
               IARG=IARG+1
               READ(STRING,'(I4)') MZ
               IF (MZ.LT.1) MZ=1
               WRITE(*,"(' UNIT CELL WILL BE REPEATIED INTO DIRECTIONS X,Y,Z: ',3I3,' TIMES')") &
                    MX,MY,MZ
            ELSE IF (+LINE(2:2).EQ.'H') THEN
               WRITE(0,'("STRC2XYZ CONVERTS THE STRUCTURE FILE OF THE PAW-INPUT (STRC) TO AN XYZ FILE")')
               WRITE(0,*)
               WRITE(0,'("USAGE: STRC2XYZ [OPTIONS] ROOTNAME")')
               WRITE(0,*)
               WRITE(0,'("  OPTIONS: ")')
               WRITE(0,'("   -O : USE STRC_OUT FOR GENERATING THE XYZ FILE")')
               WRITE(0,'("   -F : PRINT FORMATTED OUTPUT TO BE USED AS A STRUCUTRE INPUT")')
               WRITE(0,'("   -M MX MY MZ : MULTIPLY THE UNIT CELL INTO DIRECTIONS X,Y AND Z")')
               WRITE(0,'("   -C : OPTIMIZE THE UNIT-CELL FOR ISOLATED MOLECULES")')
               WRITE(0,'("   -D NUMBER : TARGET DISTANCE FOR UNIT-CELL OPTIMIZATIONS")')
               WRITE(0,'("        THAT DISTANCE IS NOT ALWAYS FULFILLED")')
               STOP
            ELSE
               PRINT*,'UNRECOGNIZED OPTION: ',TRIM(LINE)
            END IF
         ELSE
            ROOTNAME=TRIM(LINE)
         END IF
      END DO

      IF(LEN(TRIM(ROOTNAME)).EQ.0) THEN
         ROOTNAME=-'CASE'
      END IF
      I=INDEX(ROOTNAME,'/',BACK=.TRUE.)
      OBJECTNAME=ROOTNAME(I+1:)
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      IF(TSTRCOUT) THEN
         PRINT*,'USING THE .STRC_OUT FILE'
         CALL FILEHANDLER$SETFILE('STRC',.TRUE.,-'.STRC_OUT')
      ELSE
         PRINT*,'USING THE STRC FILE'
         CALL FILEHANDLER$SETFILE('STRC',.TRUE.,-'.STRC')
      END IF
      CALL FILEHANDLER$SETSPECIFICATION('STRC','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('XYZ',.TRUE.,-'.XYZ')
      CALL FILEHANDLER$SETSPECIFICATION('XYZ','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('XYZ','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('XYZ','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('XYZ','FORM','FORMATTED')
!
!     ==================================================================
!     == READ STRUCTURE FILE                                          ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('STRC',NFIL)
      CALL LINKEDLIST$READ(LL_STRC,NFIL,'~')
!
!     ==================================================================
!     == GET CONVERSION FACOTOR                                       ==
!     ==================================================================
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'GENERIC')
!
!     ==  READ LUNIT ===========================================================
      CALL LINKEDLIST$EXISTD(LL_STRC,'LUNIT[AA]',1,TCHK1)
      IF(TCHK1) THEN 
        CALL LINKEDLIST$EXISTD(LL_STRC,'LUNIT',1,TCHK)
        IF(TCHK.AND.TCHK1) THEN
          CALL ERROR$MSG('LUNIT AND LUNIT[AA] ARE MUTUALLY EXCLUSIVE')
          CALL ERROR$STOP('STRC2XYZ (MAIN)')
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'LUNIT[AA]',0,SVAR)
        CALL CONSTANTS$GET('ANGSTROM ',ANGSTROM)
        SVAR=SVAR*ANGSTROM
        CALL LINKEDLIST$SET(LL_STRC,'LUNIT',0,SVAR)
        CALL LINKEDLIST$RMDATA(LL_STRC,'LUNIT[AA]',1)
      END IF
      CALL LINKEDLIST$EXISTD(LL_STRC,'LUNIT',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'LUNIT',0,1.D0)
!
      CALL LINKEDLIST$GET(LL_STRC,'LUNIT',1,RUNIT)
!
!     ==================================================================
!     == GET LATTICE VECTORS                                          ==
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
      ALLOCATE(SPECIES(NAT))
      ALLOCATE(R(3,NAT))
      ALLOCATE(NEIGH(8,NAT))
      NEIGH=0
      DO IAT=1,NAT
        CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IAT)
        CALL LINKEDLIST$GET(LL_STRC,'R',1,R(:,IAT))
        R(:,IAT)=R(:,IAT)*RUNIT
        CALL LINKEDLIST$GET(LL_STRC,'NAME',1,NAME(IAT))
        CALL LINKEDLIST$EXISTD(LL_STRC,'SP',1,TCHK)
        IF(.NOT.TCHK) THEN
          SPECIES(IAT)=NAME(IAT)(1:2)
        ELSE
          CALL LINKEDLIST$GET(LL_STRC,'SP',1,SPECIES(IAT))
        END IF
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
!
!     ==================================================================
!     == OPTIMIZE UNIT-CELL IF DESIRED                                ==
!     ==================================================================
      IF(TCELL) THEN
         CALL OPTIMIZE_RBAS(NAT,R,RBAS,NAME,CELLDISTANCE)
         STOP
      END IF
!
!     ==================================================================
!     == CONVERT DATA TO ANGSTROM                                     ==
!     ==================================================================
      RBAS=RBAS/ANGSTROM
      R=R/ANGSTROM
!
!     ==================================================================
!     == CHECK THE DISTANCE BETWEEN ATOMS IN ONE UNIT-CELL            ==
!     == CHECK THE NAMES OF THE ATOMS                                 ==
!     ==================================================================
      DO IAT=1,NAT
         DO JAT=IAT+1,NAT
            SVAR=SQRT(SUM((R(:,IAT)-R(:,JAT))**2))
            IF(SVAR.LT.0.5D0) THEN 
               WRITE(*,"('WARNING: ATOM  ',A,'  AND  ',A, '  ARE &
                    &TOO CLOSE. DISTANCE: ',F10.5,' ANG.')") &   
                    TRIM(NAME(IAT)),TRIM(NAME(JAT)),SVAR
               PAWOK=.FALSE.
            END IF
            IF(TRIM(NAME(IAT)).EQ.TRIM(NAME(JAT))) THEN
               WRITE(*,"('WARNING: THE ATOM& 
                    &NAME ',A,' IS USED FOR TWO ATOMS, WICH IS NOT ALLOWED!')")& 
                    TRIM(NAME(IAT))
               PAWOK=.FALSE.
            END IF
         END DO
      END DO
      CALL NEXTDIST(NAT,R,RBAS,0,DIST,IAT,JAT,DIR)
      IF(IAT.EQ.-2) PRINT*,'THIS IS DUE TO STRANGE COMPILER BEHAVIOUR!'
      WRITE(*,"('CLOSEST ATOMS BETWEEN THE PERIODIC IMAGES ARE ATOMS ',A,' AND ',A,&
           &' INTO DIRECTION ',3I3,' WITH DISTANCE: ',F8.3,' ANG')")&
           TRIM(NAME(IAT)),TRIM(NAME(JAT)),DIR(1:3),DIST
!
!     ==================================================================
!     == WRITE FORMATTED BLOCK WITH ATOM POSITIONS                    ==
!     ==================================================================
      IF(TFORMAT) THEN
         DO IAT=1,NAT
            WRITE(*,'(A,A,T22,A,3F13.7,A)') "!ATOM NAME='"//TRIM(NAME(IAT)),&
                 "'"," R=",R(:,IAT)," !END"
         END DO
      END IF
!
!     ==================================================================
!     == WRITE XYZ FILE                                               ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('XYZ',NFIL)
      WRITE(NFIL,*)NAT*MX*MY*MZ
      WRITE(NFIL,FMT='(3X,A30)')ROOTNAME
      DO IX=0,MX-1
         DO IY=0,MY-1
            DO IZ=0,MZ-1
               DO IAT=1,NAT
!ORIGINAL VERSION:
!THE NAME OF THE ATOM WAS USED
!                  STR=NAME(IAT)(1:2)
!UPDATED VERSION:
!THE SPECIES IS USED
                  STR=SPECIES(IAT)(1:2)
                  IF(STR(2:2).EQ.'_') STR(2:2)=' '
                  WRITE(NFIL,FMT='(A2,2X,3(F10.5,1X),2X,A)')STR,R(:,IAT)+ &
                       &      IX*RBAS(:,1)+IY*RBAS(:,2)+IZ*RBAS(:,3),TRIM(NAME(IAT))
               ENDDO
            END DO
         END DO
      END DO
      
      CALL FILEHANDLER$FILENAME('XYZ',OBJECTNAME)
      WRITE(*,FMT='("FILE WRITTEN: ",A)') TRIM(ADJUSTL(OBJECTNAME))
      CALL FILEHANDLER$CLOSEALL
      IF(PAWOK.AND.DIST.GT.6.D0) WRITE(*,*) ' THE STRC FILE IS OK FOR A PAW SIMULATION'
      IF(PAWOK.AND.DIST.LE.6.D0) WRITE(*,*) ' THE STRC FILE IS OK FOR A PAW SIMULATION &
           &BUT THE CELL IS TOO SMALL FOR ISOLTED MOLECULES'
      IF(.NOT.PAWOK) WRITE(*,*) ' THE STRC FILE CONTAIONS ERRORS - A PAW SIMULATION WILL NOT WORK'
      CALL ERROR$NORMALSTOP()
      STOP
      END PROGRAM STRC2XYZ
!
!     ..................................................................
      SUBROUTINE OPTIMIZE_RBAS(NAT,R,RBAS_ORIG,NAME,CELLDISTANCE)
!     ******************************************************************
!     ** THIS PROGRAM OPTIMIZES THE LATTICE VECTORS OF AN ISOLATED     
!     ** MOLECULE. THE ACTIONS OF THE PROGRAM AND UNDERLYING IDEAS ARE  
!     ** DESCRIBED IN MOLECULE_RBAS.TEX
!     ** IMPLEMENTED INTO STRC2XYZ ON 2003 01 10 BY JOHANNES SCHIMPL
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NAT
      REAL(8)      ,INTENT(IN) :: R(3,NAT)
      REAL(8)      ,INTENT(IN) :: RBAS_ORIG(3,3)
      CHARACTER(32),INTENT(IN) :: NAME(NAT)
      REAL(8)      ,INTENT(IN) :: CELLDISTANCE
      REAL(8)                  :: SVAR,DIST,SVAR1,SVAR2
      INTEGER(4)               :: AT1,AT2,DIROUT(3),I,J,K,ISTEP
      REAL(8)                  :: RBASX,RBASY,RBASZ,A
      REAL(8)                  :: RBAS_P(3,3),RBAS_F(3,3),RBAS_B(3,3),RBAS_OLD(3,3)
      REAL(8)     ,ALLOCATABLE :: U(:,:),KOV(:,:),E(:),COMPONENT(:)
!     ******************************************************************
      CALL TIMING$START ! NNEDED FOR THE PAW-LIBS
!
!     ==================================================================
!     == REPORT PROPERTIES OF THE OLD LATTICE VECTORS
!     ================================================================== 
      PRINT*,'OLD LATTICE VECTORS (ANG):'
      WRITE(*,"(3F9.5)") RBAS_ORIG/1.889726878D0
      SVAR=RBAS_ORIG(1,1)*(RBAS_ORIG(2,2)*RBAS_ORIG(3,3)-RBAS_ORIG(2,3)*RBAS_ORIG(3,2)) &
           -RBAS_ORIG(1,2)*(RBAS_ORIG(2,1)*RBAS_ORIG(3,3)-RBAS_ORIG(2,3)*RBAS_ORIG(3,1)) &
           +RBAS_ORIG(1,3)*(RBAS_ORIG(2,1)*RBAS_ORIG(3,2)-RBAS_ORIG(2,2)*RBAS_ORIG(3,1))
      WRITE(*,"('OLD CELL VOLUME: ',F15.5)") ABS(SVAR)

      CALL NEXTDIST(NAT,R,RBAS_ORIG,0,DIST,AT1,AT2,DIROUT)
      WRITE(*,"('CLOSEST ARE ATOMS ',A,', ',A,' INTO DIRECTION ',3I2.1,' WITH DISTANCE: ', &
           &F8.3,' ANG')") TRIM(NAME(AT1)),TRIM(NAME(AT2)),DIROUT,DIST/1.889726878D0 
!
!     ==================================================================
!     == CALCULATE COVARIANCE AND DIAGONALIZE IT
!     ================================================================== 
      ALLOCATE(KOV(3,3))
      DO I=1,3
         DO J=I,3
            SVAR=0.D0
            DO K=1,NAT
               SVAR=SVAR+R(I,K)*R(J,K)
            END DO
            SVAR1=0.D0
            DO K=1,NAT
               SVAR1=SVAR1+R(I,K)
            END DO
            SVAR2=0.D0
            DO K=1,NAT
               SVAR2=SVAR2+R(J,K)
            END DO
            KOV(I,J)=(SVAR-SVAR1*SVAR2)/NAT
            KOV(J,I)=KOV(I,J)
         END DO
      END DO
      ALLOCATE(U(3,3))
      ALLOCATE(E(3))
      CALL LIB$DIAGR8(3,KOV,E,U)
      DEALLOCATE(KOV)

!!$       PRINT*,'PRINCIPAL COMPONENTS WITH EIGENVALUES (4. COL)'
!!$       WRITE(*,"('1.:',4F9.3)") U(:,1),E(1)       
!!$       WRITE(*,"('2.:',4F9.3)") U(:,2),E(2)         
!!$       WRITE(*,"('3.:',4F9.3)") U(:,3),E(3)        

       DEALLOCATE(E)
!
!     ==================================================================
!     == FIRST ESTIMATE OF THE LENGTH
!     ==================================================================
       PRINT*
       PRINT*,'---------------- PRIMITIVE -------------------------' 
       ALLOCATE(COMPONENT(NAT))
       DO I=1,3  !3 LATTICE VECTORS
          DO K=1,NAT
             COMPONENT(K)=DOT_PRODUCT(U(:,I),R(:,K))
          END DO
          SVAR=MAXVAL(COMPONENT)-MINVAL(COMPONENT)+CELLDISTANCE
          RBAS_P(:,I)=SVAR*U(:,I)
       END DO
!
!     ==================================================================
!     == SCALE LENGTH OF PRIMITIVE LATTICE
!     ================================================================== 
       RBAS_OLD=RBAS_P
       DO ISTEP=1,10
          DO K=1,3
             CALL NEXTDIST(NAT,R,RBAS_P,K,DIST,AT1,AT2,DIROUT)
             I=AT1
             J=AT2
             RBASX=RBAS_P(1,K)
             RBASY=RBAS_P(2,K)
             RBASZ=RBAS_P(3,K)
             A=(SQRT(-R(1,I)**2*(RBASY**2+RBASZ**2)+2*R(1,I)*(R(1,J)*(RBASY**2+RBASZ**2)+&
                  RBASX*(R(2,I)*RBASY-R(2,J)*RBASY+R(3,I)*RBASZ-R(3,J)*RBASZ)) &
                  -R(1,J)**2*(RBASY**2+RBASZ**2)-2*R(1,J)*RBASX*(R(2,I)*RBASY  &
                  -R(2,J)*RBASY+RBASZ*(R(3,I)-R(3,J)))-RBASX**2*(R(2,I)**2-2*R(2,I)*R(2,J)+R(2,J)**2 &
                  +R(3,I)**2-2*R(3,I)*R(3,J)+R(3,J)**2-CELLDISTANCE**2)-R(2,I)**2*RBASZ**2+ &
                  2*R(2,I)*RBASZ*(R(2,J)*RBASZ+RBASY*(R(3,I)-R(3,J)))-R(2,J)**2*RBASZ**2+2*R(2,J) &
                  *RBASY*RBASZ*(R(3,J)-R(3,I))-RBASY**2*(R(3,I)**2-2*R(3,I)*R(3,J)+R(3,J)**2 &
                  -CELLDISTANCE**2)+RBASZ**2*CELLDISTANCE**2)+R(1,I)*RBASX-R(1,J)*RBASX+ &
                  R(2,I)*RBASY-R(2,J)*RBASY+RBASZ*(R(3,I)-R(3,J)))/(RBASX**2+RBASY**2+RBASZ**2)
             RBAS_P(:,K)=RBAS_P(:,K)*A
          END DO
          IF (MAXVAL(ABS(RBAS_OLD-RBAS_P)).LT.0.1D0) EXIT
          RBAS_OLD=RBAS_P
       END DO

       ! MAKE FIRST LATTICE VECTOR LONGEST (FOR PARALLELISATION)
       DO I=2,3
          SVAR=RBAS_P(1,1)**2+RBAS_P(2,1)**2+RBAS_P(3,1)**2
          SVAR1=RBAS_P(1,I)**2+RBAS_P(2,I)**2+RBAS_P(3,I)**2
          IF(SVAR1.GT.SVAR) THEN !EXCHANGE THEM
             RBAS_OLD(:,1)=RBAS_P(:,1)
             RBAS_P(:,1)=RBAS_P(:,I)
             RBAS_P(:,I)=RBAS_OLD(:,1)
          END IF
       END DO

       PRINT*,'PRIMITIVE LATTICE (ANG):'
       WRITE(*,"(3F15.9)") RBAS_P/1.889726878D0
       SVAR=RBAS_P(1,1)*(RBAS_P(2,2)*RBAS_P(3,3)-RBAS_P(2,3)*RBAS_P(3,2)) &
           -RBAS_P(1,2)*(RBAS_P(2,1)*RBAS_P(3,3)-RBAS_P(2,3)*RBAS_P(3,1)) &
           +RBAS_P(1,3)*(RBAS_P(2,1)*RBAS_P(3,2)-RBAS_P(2,2)*RBAS_P(3,1))
       WRITE(*,"('CELL VOLUME: ',F15.5)") ABS(SVAR)       

       CALL NEXTDIST(NAT,R,RBAS_P,0,DIST,AT1,AT2,DIROUT)
       WRITE(*,"('CLOSEST ARE ATOMS ',A,', ',A,' INTO DIRECTION ',3I2.1,' WITH DISTANCE: ', &
            &F8.3,' ANG')") TRIM(NAME(AT1)),TRIM(NAME(AT2)),DIROUT,DIST/1.889726878D0 

!
!     ==================================================================
!     == FACE-CENTERED LATTICE
!     ================================================================== 
       PRINT*
       PRINT*,'---------------- FACE CENTERED -------------------------' 
       RBAS_F(:,1)=(RBAS_P(:,1)+RBAS_P(:,2))!/SQRT(2.)
       RBAS_F(:,2)=(RBAS_P(:,2)+RBAS_P(:,3))!/SQRT(2.)
       RBAS_F(:,3)=(RBAS_P(:,1)+RBAS_P(:,3))!/SQRT(2.)

!!$         PRINT*,'FIRST ESTIMATE  OF FACE-CENTERED LATTICE (ANG):'
!!$       WRITE(*,"(3F9.5)") RBAS_F/1.889726878D0
!!$       SVAR=RBAS_F(1,1)*(RBAS_F(2,2)*RBAS_F(3,3)-RBAS_F(2,3)*RBAS_F(3,2)) &
!!$           -RBAS_F(1,2)*(RBAS_F(2,1)*RBAS_F(3,3)-RBAS_F(2,3)*RBAS_F(3,1)) &
!!$           +RBAS_F(1,3)*(RBAS_F(2,1)*RBAS_F(3,2)-RBAS_F(2,2)*RBAS_F(3,1))
!!$       WRITE(*,"('CELL VOLUME: ',F15.5)") ABS(SVAR)       
!!$
!!$       CALL NEXTDIST(NAT,R,RBAS_F,0,DIST,AT1,AT2,DIROUT)
!!$       WRITE(*,"('CLOSEST ARE ATOMS ',I3.0,', ',I3.0,' INTO DIRECTION ',I3.0,' WITH DISTANCE: ', &
!!$            &F8.3,' ANG')") AT1,AT2,DIROUT,DIST/1.889726878D0     

       RBAS_OLD=RBAS_F
       DO ISTEP=1,10
          DO K=1,3
             CALL NEXTDIST(NAT,R,RBAS_F,K,DIST,AT1,AT2,DIROUT)
             I=AT1
             J=AT2
             RBASX=RBAS_F(1,K)
             RBASY=RBAS_F(2,K)
             RBASZ=RBAS_F(3,K)
             A=(SQRT(-R(1,I)**2*(RBASY**2+RBASZ**2)+2*R(1,I)*(R(1,J)*(RBASY**2+RBASZ**2)+&
                  RBASX*(R(2,I)*RBASY-R(2,J)*RBASY+R(3,I)*RBASZ-R(3,J)*RBASZ)) &
                  -R(1,J)**2*(RBASY**2+RBASZ**2)-2*R(1,J)*RBASX*(R(2,I)*RBASY  &
                  -R(2,J)*RBASY+RBASZ*(R(3,I)-R(3,J)))-RBASX**2*(R(2,I)**2-2*R(2,I)*R(2,J)+R(2,J)**2 &
                  +R(3,I)**2-2*R(3,I)*R(3,J)+R(3,J)**2-CELLDISTANCE**2)-R(2,I)**2*RBASZ**2+ &
                  2*R(2,I)*RBASZ*(R(2,J)*RBASZ+RBASY*(R(3,I)-R(3,J)))-R(2,J)**2*RBASZ**2+2*R(2,J) &
                  *RBASY*RBASZ*(R(3,J)-R(3,I))-RBASY**2*(R(3,I)**2-2*R(3,I)*R(3,J)+R(3,J)**2 &
                  -CELLDISTANCE**2)+RBASZ**2*CELLDISTANCE**2)+R(1,I)*RBASX-R(1,J)*RBASX+ &
                  R(2,I)*RBASY-R(2,J)*RBASY+RBASZ*(R(3,I)-R(3,J)))/(RBASX**2+RBASY**2+RBASZ**2)

!!$             RBAS_F(:,K)=RBAS_F(:,K)*A
!!$             CALL NEXTDIST(NAT,R,RBAS_F,0,DIST,AT1,AT2,DIROUT)
!!$             IF (DIST.LT.CELLDISTANCE) RBAS_F(:,K)=RBAS_F(:,K)/A
!!$          END DO
!!$          IF (MAXVAL(ABS(RBAS_OLD-RBAS_F)).LT.0.1D0.AND.DIST.GE.CELLDISTANCE) EXIT
!!$          RBAS_OLD=RBAS_F

             RBAS_F(:,K)=RBAS_F(:,K)*A
          END DO
          IF (MAXVAL(ABS(RBAS_OLD-RBAS_F)).LT.0.1D0) EXIT
          RBAS_OLD=RBAS_F
       END DO

       ! MAKE FIRST LATTICE VECTOR LONGEST (FOR PARALLELISATION)
       DO I=2,3
          SVAR=RBAS_F(1,1)**2+RBAS_F(2,1)**2+RBAS_F(3,1)**2
          SVAR1=RBAS_F(1,I)**2+RBAS_F(2,I)**2+RBAS_F(3,I)**2
          IF(SVAR1.GT.SVAR) THEN !EXCHANGE THEM
             RBAS_OLD(:,1)=RBAS_F(:,1)
             RBAS_F(:,1)=RBAS_F(:,I)
             RBAS_F(:,I)=RBAS_OLD(:,1)
          END IF
       END DO

       PRINT*,'FACE-CENTERED LATTICE (ANG):'
       WRITE(*,"(3F15.9)") RBAS_F/1.889726878D0
       SVAR=RBAS_F(1,1)*(RBAS_F(2,2)*RBAS_F(3,3)-RBAS_F(2,3)*RBAS_F(3,2)) &
           -RBAS_F(1,2)*(RBAS_F(2,1)*RBAS_F(3,3)-RBAS_F(2,3)*RBAS_F(3,1)) &
           +RBAS_F(1,3)*(RBAS_F(2,1)*RBAS_F(3,2)-RBAS_F(2,2)*RBAS_F(3,1))
       WRITE(*,"('CELL VOLUME: ',F15.5)") ABS(SVAR)       

       CALL NEXTDIST(NAT,R,RBAS_F,0,DIST,AT1,AT2,DIROUT)
       WRITE(*,"('CLOSEST ARE ATOMS ',A,', ',A,' INTO DIRECTION ',3I2.1,' WITH DISTANCE: ', &
            &F8.3,' ANG')") TRIM(NAME(AT1)),TRIM(NAME(AT2)),DIROUT,DIST/1.889726878D0 
!
!     ==================================================================
!     == BODY-CENTERED LATTICE
!     ================================================================== 
       PRINT*
       PRINT*,'---------------- BODY CENTERED -------------------------'  
       RBAS_B(:,1)=(RBAS_P(:,1)+RBAS_P(:,2)+RBAS_P(:,3))*0.6D0!/SQRT(2.)
       RBAS_B(:,2)=(-RBAS_P(:,1)-RBAS_P(:,2)+RBAS_P(:,3))*0.6D0!/SQRT(2.)
       RBAS_B(:,3)=(RBAS_P(:,1)-RBAS_P(:,2)-RBAS_P(:,3))*0.6D0!/SQRT(2.)

!!$         PRINT*,'FIRST ESTIMATE  OF BODY-CENTERED LATTICE (ANG):'
!!$       WRITE(*,"(3F9.5)") RBAS_B/1.889726878D0
!!$       SVAR=RBAS_B(1,1)*(RBAS_B(2,2)*RBAS_B(3,3)-RBAS_B(2,3)*RBAS_B(3,2)) &
!!$           -RBAS_B(1,2)*(RBAS_B(2,1)*RBAS_B(3,3)-RBAS_B(2,3)*RBAS_B(3,1)) &
!!$           +RBAS_B(1,3)*(RBAS_B(2,1)*RBAS_B(3,2)-RBAS_B(2,2)*RBAS_B(3,1))
!!$       WRITE(*,"('CELL VOLUME: ',F15.5)") ABS(SVAR)       
!!$
!!$       CALL NEXTDIST(NAT,R,RBAS_B,0,DIST,AT1,AT2,DIROUT)
!!$       WRITE(*,"('CLOSEST ARE ATOMS ',I3.0,', ',I3.0,' INTO DIRECTION ',I3.0,' WITH DISTANCE: ', &
!!$            &F8.3,' ANG')") AT1,AT2,DIROUT,DIST/1.889726878D0     


       RBAS_OLD=RBAS_B
       DO ISTEP=1,10
          DO K=1,3
             CALL NEXTDIST(NAT,R,RBAS_B,K,DIST,AT1,AT2,DIROUT)
             I=AT1
             J=AT2
             RBASX=RBAS_B(1,K)
             RBASY=RBAS_B(2,K)
             RBASZ=RBAS_B(3,K)
             A=(SQRT(-R(1,I)**2*(RBASY**2+RBASZ**2)+2*R(1,I)*(R(1,J)*(RBASY**2+RBASZ**2)+&
                  RBASX*(R(2,I)*RBASY-R(2,J)*RBASY+R(3,I)*RBASZ-R(3,J)*RBASZ)) &
                  -R(1,J)**2*(RBASY**2+RBASZ**2)-2*R(1,J)*RBASX*(R(2,I)*RBASY  &
                  -R(2,J)*RBASY+RBASZ*(R(3,I)-R(3,J)))-RBASX**2*(R(2,I)**2-2*R(2,I)*R(2,J)+R(2,J)**2 &
                  +R(3,I)**2-2*R(3,I)*R(3,J)+R(3,J)**2-CELLDISTANCE**2)-R(2,I)**2*RBASZ**2+ &
                  2*R(2,I)*RBASZ*(R(2,J)*RBASZ+RBASY*(R(3,I)-R(3,J)))-R(2,J)**2*RBASZ**2+2*R(2,J) &
                  *RBASY*RBASZ*(R(3,J)-R(3,I))-RBASY**2*(R(3,I)**2-2*R(3,I)*R(3,J)+R(3,J)**2 &
                  -CELLDISTANCE**2)+RBASZ**2*CELLDISTANCE**2)+R(1,I)*RBASX-R(1,J)*RBASX+ &
                  R(2,I)*RBASY-R(2,J)*RBASY+RBASZ*(R(3,I)-R(3,J)))/(RBASX**2+RBASY**2+RBASZ**2)
             RBAS_B(:,K)=RBAS_B(:,K)*A
             CALL NEXTDIST(NAT,R,RBAS_B,0,DIST,AT1,AT2,DIROUT)
             IF (DIST.LT.CELLDISTANCE) RBAS_B(:,K)=RBAS_B(:,K)/A
          END DO
          IF (MAXVAL(ABS(RBAS_OLD-RBAS_B)).LT.0.1D0.AND.DIST.GE.CELLDISTANCE) EXIT
          RBAS_OLD=RBAS_B
       END DO

       ! MAKE FIRST LATTICE VECTOR LONGEST (FOR PARALLELISATION)
       DO I=2,3
          SVAR=RBAS_B(1,1)**2+RBAS_B(2,1)**2+RBAS_B(3,1)**2
          SVAR1=RBAS_B(1,I)**2+RBAS_B(2,I)**2+RBAS_B(3,I)**2
          IF(SVAR1.GT.SVAR) THEN !EXCHANGE THEM
             RBAS_OLD(:,1)=RBAS_B(:,1)
             RBAS_B(:,1)=RBAS_B(:,I)
             RBAS_B(:,I)=RBAS_OLD(:,1)
          END IF
       END DO

       PRINT*,'BODY-CENTERED LATTICE (ANG):'
       WRITE(*,"(3F15.9)") RBAS_B/1.889726878D0
       SVAR=RBAS_B(1,1)*(RBAS_B(2,2)*RBAS_B(3,3)-RBAS_B(2,3)*RBAS_B(3,2)) &
           -RBAS_B(1,2)*(RBAS_B(2,1)*RBAS_B(3,3)-RBAS_B(2,3)*RBAS_B(3,1)) &
           +RBAS_B(1,3)*(RBAS_B(2,1)*RBAS_B(3,2)-RBAS_B(2,2)*RBAS_B(3,1))
       WRITE(*,"('CELL VOLUME: ',F15.5)") ABS(SVAR)       

       CALL NEXTDIST(NAT,R,RBAS_B,0,DIST,AT1,AT2,DIROUT)
       WRITE(*,"('CLOSEST ARE ATOMS ',A,', ',A,' INTO DIRECTION ',3I2.1,' WITH DISTANCE: ', &
            &F8.3,' ANG')") TRIM(NAME(AT1)),TRIM(NAME(AT2)),DIROUT,DIST/1.889726878D0 
     END SUBROUTINE OPTIMIZE_RBAS
!
!     ..................................................................
      SUBROUTINE NEXTDIST(NAT,R,RBAS,DIRIN,DIST,AT1,AT2,DIROUT)
!     ******************************************************************
!     ** CALCULATES THE CLOSEST DISTANCE BETWEEN ATOMS OF 2 PERIODIC
!     ** IMAGES AND GIVE BACK THE NUMBERS OF THESE 2 ATOMS, THE DISTANCE
!     ** AND THE DIRECTION (NUMBER OF LATICE VECTOR OF THE CLOSEST 
!     ** DISTANCE
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      INTEGER(4),INTENT(IN) :: DIRIN      
      REAL(8)   ,INTENT(OUT):: DIST
      INTEGER(4),INTENT(OUT):: AT1,AT2,DIROUT(3) 
      INTEGER(4)            :: I,IAT1,IAT2,FROM(3),TO(3),A1,A2,A3
      REAL(8)               :: BESTDIST,SVAR
!     ******************************************************************
      IF (DIRIN.EQ.0) THEN
         FROM=1
         TO=-1
      ELSE
         FROM=0
         TO=0
         FROM(DIRIN)=1
         TO(DIRIN)=1
      END IF
      DIST=5000000.D0
      DO A1=FROM(1),TO(1),-1
         DO A2=FROM(2),TO(2),-1
            DO A3=FROM(3),TO(3),-1
               IF (A1.EQ.0.AND.A2.EQ.0.AND.A3.EQ.0) CYCLE
               DO IAT1=1,NAT
                  DO IAT2=1,NAT
                     SVAR=SQRT((R(1,IAT1)-(R(1,IAT2)+A1*RBAS(1,1)+A2*RBAS(1,2)+A3*RBAS(1,3) ))**2+ &
                          (R(2,IAT1)-(R(2,IAT2)+A1*RBAS(2,1)+A2*RBAS(2,2)+A3*RBAS(2,3) ))**2+ &
                          (R(3,IAT1)-(R(3,IAT2)+A1*RBAS(3,1)+A2*RBAS(3,2)+A3*RBAS(3,3) ))**2)
                     IF (SVAR.LT.DIST) THEN
                        DIST=SVAR
                        DIROUT(1)=A1
                        DIROUT(2)=A2
                        DIROUT(3)=A3
                        AT1=IAT1
                        AT2=IAT2
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
    END SUBROUTINE NEXTDIST
!
!     ....................................................................
    SUBROUTINE CHAR_TO_REAL(STR,VAL,TCHK)
!     ********************************************************************
!     ** TRANSFORMS A STRING INTO A REAL VARIABLE                       **
!     ** STR MAY HAVE THE SHAPES: 2.3234 324587 2345. .32455            **
!     ** ON EXIT, TCHK=.TRUE. IF SUCCESSFUL                             **
!     ********************************************************************
      IMPLICIT NONE
      CHARACTER(*)  ,INTENT(IN)   :: STR
      REAL(8)       ,INTENT(OUT)  :: VAL
      LOGICAL(4)    ,INTENT(OUT)  :: TCHK
!     ********************************************************************
      READ(STR,FMT=*,ERR=123) VAL
      TCHK=.TRUE.
      RETURN
123   TCHK=.FALSE.
      VAL=0.D0
    END SUBROUTINE CHAR_TO_REAL
