MODULE CONSTRAINTS_MODULE
!***********************************************************************
!***********************************************************************
!**                                                                   **
!**  NAME: CONSTRAINTS_MODULE                                          **
!**                                                                   **
!**  SETORIENTATION                                                   **
!**  SETANGULARMOMENTUM                                               **
!**  SETTRANSLATIONALMOMENTUM                                         **
!**  SETFIXGROUP                                                      **
!**  SETFIXATOM                                                       **
!**  SETRIGID                                                         **
!**  SETSADDLE                                                        **
!**                                                                   **
!**  SETSADDLES                                                       **
!**  SADDLEDADS                                                       **
!**  SETREFERENCE                                                     **
!**  NFREE                                                            **
!**                                                                   **
!**  SETCONSTR                                                        **
!**                                                                   **
!**  REPORT                                                           **
!***********************************************************************
!***********************************************************************
USE LINKEDLIST_MODULE
TYPE(LL_TYPE)      :: LL_CNSTR
LOGICAL(4)         :: TINI=.FALSE.
LOGICAL(4)         :: TSET=.FALSE.
LOGICAL(4)         :: TOPEN=.FALSE.
LOGICAL(4)         :: TON=.FALSE.   ! ON/OFF SWITCH
LOGICAL(4)         :: TREF=.FALSE.  ! REFERENCE STRUCTURE SET
INTEGER(4)         :: NC=0          ! #(CONSTRAINTS)
TYPE CONSTRAINT_TYPE 
  REAL(8)            :: C0
  INTEGER(4)         :: N1
  INTEGER(4),POINTER :: I1(:)
  REAL(8)   ,POINTER :: C1(:)
  INTEGER(4)         :: N2
  INTEGER(4),POINTER :: I2(:)
  INTEGER(4),POINTER :: J2(:)
  REAL(8)   ,POINTER :: C2(:)
  REAL(8)            :: LAMBDAP
  REAL(8)            :: LAMBDA0
  REAL(8)            :: LAMBDAM
  REAL(8)            :: LAMBDAMM
  LOGICAL(4)         :: SHOW
  CHARACTER(64)      :: DESCRIPTION
  INTEGER(4)         :: NSTEP
  REAL(8)            :: C0FINAL
END TYPE CONSTRAINT_TYPE
END MODULE CONSTRAINTS_MODULE
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_INITIALIZE
!     ******************************************************************
!     **  CONSTRAINTS$SET                                             **
!     ******************************************************************
      USE CONSTRAINTS_MODULE
!     ******************************************************************
      IF(TINI) RETURN
      CALL LINKEDLIST$NEW(LL_CNSTR)
      TINI=.TRUE.
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS$OPEN(TYPE_,TITLE_)
!     ******************************************************************
!     **  CONSTRAINTS$OPEN                                            **
!     ******************************************************************
      USE CONSTRAINTS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: TYPE_
      CHARACTER(*),INTENT(IN) :: TITLE_
      CHARACTER(32)           :: TYPE
      CHARACTER(256)          :: STRING
!     ******************************************************************
      CALL CONSTRAINTS_INITIALIZE
!
      IF(TSET) THEN 
        CALL ERROR$MSG('ADDING OF FURTHER CONSTRAINTS NOT ALLOWED')
        CALL ERROR$STOP('CONSTRAINTS$OPEN')
      END IF
!
      IF(TOPEN) THEN
        CALL ERROR$MSG('CANNOT ADD CONSTRAINT WHILE ADDING ANOTHER')
        CALL ERROR$STOP('CONSTRAINTS$OPEN')
      END IF
!
      CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
      CALL LINKEDLIST$SELECT(LL_CNSTR,'CONSTRAINT',-1)
      TYPE=TYPE_
      CALL LINKEDLIST$SET(LL_CNSTR,'TYPE',0,TYPE)
      STRING=TITLE_
      CALL LINKEDLIST$SET(LL_CNSTR,'TITLE',0,STRING)
      TOPEN=.TRUE.
      TON=.TRUE.
      RETURN 
      END
#TEMPLATE CONSTRAINTS$DEFINE
(<TYPEID><TYPE>)=([R8],[REAL(8)])
                 ([I4],[INTEGER(4)])
                 ([L4][LOGICAL(4)])
                 ([CH][CHARACTER(*)])
(<RANKID><,NBYTE_><(NBYTE_)><!>)=([][][][!])
                            ([A][,NBYTE_][(NBYTE_)][])
#BODY
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS$DEFINE<TYPEID><RANKID>(IDENT_<,NBYTE_>,VAL_)
!     ******************************************************************
!     **  CONSTRAINTS$DEFINE                                          **
!     ******************************************************************
      USE CONSTRAINTS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: IDENT_
   <!>INTEGER(4)  ,INTENT(IN) :: NBYTE_
      <TYPE>      ,INTENT(IN) :: VAL_<(NBYTE_)>
      LOGICAL(4)              :: TCHK
      CHARACTER(32)           :: TYPE
      INTEGER(4)              :: NOPTION
      CHARACTER(32)           :: OPTION(32)
      CHARACTER(256)          :: TITLE
      INTEGER(4)              :: I
!     ******************************************************************
      CALL CONSTRAINTS_INITIALIZE
      IF(TSET) THEN
        CALL ERROR$MSG('ADDING OF FURTHER CONSTRAINTS NOT ALLOWED')
        CALL ERROR$STOP('CONSTRAINTS$DEFINE')
      ELSE IF(.NOT.TOPEN) THEN
        CALL ERROR$MSG('CANNOT SET DATA UNLESS CONSTRAINT IS OPEN')
        CALL ERROR$STOP('CONSTRAINTS$DEFINE')
      END IF
!
!     ===================================================================
!     ==  CHECK WHETHER OPTION IS PERMITTED                            ==
!     ===================================================================
      IF(TOPEN) THEN
        CALL LINKEDLIST$GET(LL_CNSTR,'TYPE',1,TYPE)
        IF(TYPE.EQ.'ORIENTATION') THEN        
          OPTION(1)='GROUP'
          OPTION(2)='PHI'
          NOPTION=2
        ELSE IF(TYPE.EQ.'BONDLENGTH') THEN        
          OPTION(1)='ATOM1'
          OPTION(2)='ATOM2'
          OPTION(3)='NSTEP'
          OPTION(4)='FLOAT'
          OPTION(5)='SMASS'
          OPTION(6)='ANNES'
          OPTION(7)='MOVE'
          OPTION(8)='NSTEP'
          OPTION(9)='SFINAL'
          OPTION(10)='VFINAL'
          OPTION(11)='SHOW'
          NOPTION=11
        ELSE IF(TYPE.EQ.'ANGULARMOMENTUM') THEN        
          OPTION(1)='GROUP'
          OPTION(2)='PHI'
          OPTION(3)='SHOW'
          NOPTION=3
        ELSE IF(TYPE.EQ.'TRANSLATIONALMOMENTUM') THEN        
          OPTION(1)='GROUP'
          OPTION(2)='DIR'
          OPTION(3)='SHOW'
          NOPTION=3
        ELSE IF(TYPE.EQ.'MIDPLANE') THEN        
          OPTION(1) ='ATOM1'
          OPTION(2) ='ATOM2'
          OPTION(3) ='ATOM3'
          OPTION(4) ='FLOAT'
          OPTION(5) ='SMASS'
          OPTION(6) ='ANNES'
          OPTION(7) ='MOVE'
          OPTION(8) ='NSTEP'
          OPTION(9) ='SFINAL'
          OPTION(10)='VFINAL'
          OPTION(11)='SHOW'
          NOPTION=11
        ELSE IF(TYPE.EQ.'GENERALLINEAR') THEN        
          OPTION(1)='VEC'
          OPTION(2)='FLOAT'
          OPTION(3)='SMASS'
          OPTION(4)='ANNES'
          OPTION(5)='MOVE'
          OPTION(6)='NSTEP'
          OPTION(7)='SFINAL'
          OPTION(8)='VFINAL'
          OPTION(9)='SHOW'
          NOPTION=9
        ELSE IF(TYPE.EQ.'RIGID') THEN        
          OPTION(1)='GROUP'
          OPTION(2)='SHOW'
          NOPTION=1
        ELSE IF(TYPE.EQ.'FIXGROUP') THEN        
          OPTION(1)='GROUP'
          OPTION(2)='SHOW'
          NOPTION=2
        ELSE IF(TYPE.EQ.'FIXATOM') THEN        
          OPTION(1)='ATOM'
          OPTION(2)='SHOW'
          NOPTION=2
        ELSE IF(TYPE.EQ.'COGSEP') THEN   !FROM PMARGL
          OPTION(1)='GROUP1'
          OPTION(2)='GROUP2'
          OPTION(3)='NSTEP'
          OPTION(4)='FLOAT'
          OPTION(5)='SMASS'
          OPTION(6)='ANNES'
          OPTION(7)='MOVE'
          OPTION(8)='NSTEP'
          OPTION(9)='SFINAL'
          OPTION(10)='VFINAL'
          OPTION(11)='SHOW'
          NOPTION=11
        ELSE IF(TYPE.EQ.'BONDANGLE') THEN   !FROM PMARGL
          OPTION(1)='ATOM1'
          OPTION(2)='ATOM2'
          OPTION(3)='ATOM3'
          OPTION(4)='NSTEP'
          OPTION(5)='FLOAT'
          OPTION(6)='SMASS'
          OPTION(7)='ANNES'
          OPTION(8)='MOVE'
          OPTION(9)='NSTEP'
          OPTION(10)='SFINAL'
          OPTION(11)='VFINAL'
          OPTION(12)='SHOW'
          NOPTION=12
        ELSE IF(TYPE.EQ.'TORSION') THEN     ! FROM PMARGL
          OPTION(1)='ATOM1'
          OPTION(2)='ATOM2'
          OPTION(3)='ATOM3'
          OPTION(4)='ATOM4'
          OPTION(5)='NSTEP'
          OPTION(6)='FLOAT'
          OPTION(7)='SMASS'
          OPTION(8)='ANNES'
          OPTION(9)='MOVE'
          OPTION(10)='NSTEP'
          OPTION(11)='SFINAL'
          OPTION(12)='VFINAL'
          OPTION(13)='SHOW'
          NOPTION=13
        ELSE
          CALL ERROR$MSG('CONSTRAINT TYPE NOT RECOGNIZED')
          CALL ERROR$CHVAL('TYPE',TYPE)
          CALL LINKEDLIST$GET(LL_CNSTR,'TITLE',1,TITLE)
          CALL ERROR$CHVAL('TITLE',TITLE)
          CALL ERROR$STOP('CONSTRAINTS$OPEN')
        END IF
!
        TCHK=.FALSE.
        DO I=1,NOPTION
          IF(IDENT_.EQ.TRIM(OPTION(I)))TCHK=.TRUE.
        ENDDO
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('DATA IDENTIFIER NOT RECOGNIZED')
          CALL ERROR$CHVAL('IDENT_',IDENT_)
          CALL ERROR$STOP('CONSTRAINTS$DEFINE')
        END IF
      END IF
!
!     ===================================================================
!     ==   SET OPTION                                                  ==
!     ===================================================================
      CALL LINKEDLIST$SET(LL_CNSTR,IDENT_,0,VAL_)
      RETURN
      END
#END TEMPLATE CONSTRAINTS$DEFINE
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS$CLOSE
!     ******************************************************************
!     **  CONSTRAINTS$CLOSE                                           **
!     ******************************************************************
      USE CONSTRAINTS_MODULE
      IMPLICIT NONE
      LOGICAL(4)    :: TCHK
!     ******************************************************************
      CALL CONSTRAINTS_INITIALIZE
      IF(.NOT.TOPEN) THEN
        CALL ERROR$MSG('CANNOT CLOSE CONSTRAINT THAT IS NOT OPEN')
        CALL ERROR$STOP('CONSTRAINTS$CLOSE')
      END IF
!     == CHECK WHETHER PRINTOUT IS REQUESTED =======================
      CALL LINKEDLIST$EXISTD(LL_CNSTR,'SHOW',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNSTR,'SHOW',0,.FALSE.)
!     
      TOPEN=.FALSE.
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS$SETREFERENCE(RBAS,NAT_,R0_,RM_,RMASS_,DELT)
!     ******************************************************************
!     **  SET REFERENCE CONFIGURATION                                 **
!     ******************************************************************
      USE CONSTRAINTS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT_    ! #(ATOMS)
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: R0_(3,NAT_)
      REAL(8)   ,INTENT(IN) :: RM_(3,NAT_)
      REAL(8)   ,INTENT(IN) :: RMASS_(NAT_)
      REAL(8)   ,INTENT(IN) :: DELT
      INTEGER(4)            :: NAT      ! #(ATOMS)
      INTEGER(4)            :: NUM      ! #(CONSTRAINTS)
      INTEGER(4)            :: I,II,III,IC
      INTEGER(4)            :: NBYTE
      INTEGER(4)            :: NC1
      INTEGER(4)            :: IAT1,IAT2,IAT3,IAT4,IAT
      INTEGER(4)            :: IT1(3),IT2(3),IT3(3)
      INTEGER(4)            :: NMEMBER
      CHARACTER(32)         :: TYPE
      CHARACTER(256)        :: STRING,STRING1
      LOGICAL(4)            :: TCHK,TCHK1,SETS,SETV
      REAL(8)   ,ALLOCATABLE:: S0(:)    !(NC1)
      REAL(8)   ,ALLOCATABLE:: SM(:)    !(NC1)
      REAL(8)   ,ALLOCATABLE:: VEC(:,:) !(3,NAT)
      LOGICAL(4),ALLOCATABLE:: TMEMBER(:) !(NAT)
      REAL(8)   ,ALLOCATABLE:: VF(:)      !(NC1)
      REAL(8)   ,ALLOCATABLE:: RREF(:,:)  !(3*NAT,3*NAT,NC)
      REAL(8)               :: ANNES
      REAL(8)               :: SMASS
      REAL(8)               :: PHI(3)
!     ******************************************************************
      IF(.NOT.TON) RETURN
      NAT=NAT_
      CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
      IF(TSET) THEN
        CALL ERROR$MSG('CONSTRAINTS ALREADY INITIALIZED')
        CALL ERROR$STOP('CONSTRAINTS$SETREFERENCE')
      END IF
      TSET=.TRUE.
      CALL LINKEDLIST$SET(LL_CNSTR,'REFERENCE',0,R0_)
      CALL LINKEDLIST$SET(LL_CNSTR,'ICOUNT',0,1)
      ALLOCATE(RREF(3,NAT))
      CALL LINKEDLIST$GET(LL_CNSTR,'REFERENCE',1,RREF)
      TSET=.TRUE.
      TREF=.TRUE.
      IC=0
      CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
      CALL LINKEDLIST$NLISTS(LL_CNSTR,'CONSTRAINT',NUM)
      DO I=1,NUM
        CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
        CALL LINKEDLIST$SELECT(LL_CNSTR,'CONSTRAINT',I)
        CALL LINKEDLIST$GET(LL_CNSTR,'TYPE',1,TYPE)
!
!       ==  BONDLENGTH =================================================
        STRING='BONDLENGTH'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM1',1,STRING)
          CALL ATOMLIST$INDEXX(STRING,IAT1,IT1)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM2',1,STRING)
          CALL ATOMLIST$INDEXX(STRING,IAT2,IT2)
          NC1=1
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
!
          ALLOCATE(S0(NC1))
          ALLOCATE(SM(NC1))
          CALL CONSTRAINTS_BONDLENGTHVALUE(NAT,IAT1,IT1,IAT2,IT2,RBAS,R0_,S0(1))
          CALL CONSTRAINTS_BONDLENGTHVALUE(NAT,IAT1,IT1,IAT2,IT2,RBAS,RM_,SM(1))
        END IF
!
!       ==  GENERAL LINEAR =============================================
        STRING='GENERALLINEAR'
        IF(TYPE.EQ.STRING(1:32)) THEN
          ALLOCATE(VEC(3,NAT))
          CALL LINKEDLIST$GET(LL_CNSTR,'VEC',1,VEC)
          NC1=1
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
          ALLOCATE(S0(NC1))
          ALLOCATE(SM(NC1))
          CALL CONSTRAINTS_GLCVALUE(NAT,VEC,R0_,S0(1))
          CALL CONSTRAINTS_GLCVALUE(NAT,VEC,RM_,SM(1))
          DEALLOCATE(VEC)
        END IF
!
!         ==  FIXATOM   ================================================
        STRING='FIXATOM'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT)
          NC1=3
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
!
          ALLOCATE(S0(NC1))
          ALLOCATE(SM(NC1))
          CALL CONSTRAINTS_FIXATOMVALUE(NAT,IAT,R0_,S0)
          CALL CONSTRAINTS_FIXATOMVALUE(NAT,IAT,RM_,SM)
        END IF
!
!       ==  FIXGROUP  ==================================================
        STRING='FIXGROUP'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP',1,STRING)
          ALLOCATE(TMEMBER(NAT))

          CALL GROUPLIST$MEMBERS(STRING,NAT,TMEMBER)
          NC1=0
          DO IAT=1,NAT
            IF(TMEMBER(IAT))NC1=NC1+3
          ENDDO
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
!
          ALLOCATE(S0(NC1))
          ALLOCATE(SM(NC1))
          II=1
          DO IAT=1,NAT
            IF(TMEMBER(IAT)) THEN
              CALL CONSTRAINTS_FIXATOMVALUE(NAT,IAT,R0_,S0(II))
              CALL CONSTRAINTS_FIXATOMVALUE(NAT,IAT,RM_,SM(II))
              II=II+3
            ENDIF
          ENDDO
          DEALLOCATE(TMEMBER)
        ENDIF
!
!         == RIGID =====================================================
        STRING='RIGID'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP',1,STRING)
          ALLOCATE(TMEMBER(NAT_))
          CALL GROUPLIST$MEMBERS(STRING,NAT_,TMEMBER)
          NMEMBER=0
          DO IAT=1,NAT_
            IF(TMEMBER(IAT)) NMEMBER=NMEMBER+1
          ENDDO
          NC1=3*NMEMBER-6
          IF(NMEMBER.EQ.2) THEN
            NC1=1
            CALL ERROR$MSG('USE BONDLENGTH CONSTRAINT, IF TWO ATOMS') 
            CALL ERROR$MSG('SHALL BE KEPT RIGID') 
            CALL ERROR$STOP('CONSTRAINTS$SETREFERENCE')
          END IF
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
          CALL CONSTRAINTS_RIGID$BASIS(NAT,TMEMBER,R0_,IAT1,IAT2,IAT3)
          CALL LINKEDLIST$SET(LL_CNSTR,'IAT1',0,IAT1)
          CALL LINKEDLIST$SET(LL_CNSTR,'IAT2',0,IAT2)
          CALL LINKEDLIST$SET(LL_CNSTR,'IAT3',0,IAT3)
!
          ALLOCATE(S0(NC1))
          ALLOCATE(SM(NC1))
          CALL CONSTRAINTS_RIGID$VALUE(NAT,TMEMBER,R0_,IAT1,IAT2,IAT3 &
     &                          ,NC1,S0)
          PRINT*,'IAT1,IAT2,IAT3 ',IAT1,IAT2,IAT3
          PRINT*,'NC1 ',NC1
          PRINT*,'S0 ',(S0(III),III=1,NC1)
          CALL CONSTRAINTS_RIGID$VALUE(NAT,TMEMBER,RM_,IAT1,IAT2,IAT3 &
     &                          ,NC1,SM)
          DEALLOCATE(TMEMBER)
        ENDIF
!
!         ==  TRANSLATIONALMOMENTUM  =====================================
        STRING='TRANSLATIONALMOMENTUM'
        IF(TYPE.EQ.STRING(1:32)) THEN
          NC1=1
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
          ALLOCATE(TMEMBER(NAT_))
          ALLOCATE(VEC(3,1))
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP',1,STRING)
          CALL GROUPLIST$MEMBERS(STRING,NAT,TMEMBER)
          CALL LINKEDLIST$GET(LL_CNSTR,'DIR',1,VEC)
!
          ALLOCATE(S0(NC1))
          ALLOCATE(SM(NC1))
          CALL CONSTRAINTS_COGVALUE(NAT,VEC(1,1),VEC(2,1),VEC(3,1) &
     &                 ,TMEMBER,RMASS_,R0_,S0(1))
          CALL CONSTRAINTS_COGVALUE(NAT,VEC(1,1),VEC(2,1),VEC(3,1) &
     &                 ,TMEMBER,RMASS_,RM_,SM(1))
          DEALLOCATE(TMEMBER)
          DEALLOCATE(VEC)
        ENDIF
!
!       ==  ANGULARMOMENTUM ==========================================
        STRING='ANGULARMOMENTUM'
        IF(TYPE.EQ.STRING(1:32)) THEN
          NC1=1
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
          ALLOCATE(TMEMBER(NAT_))
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP',1,STRING)
          CALL GROUPLIST$MEMBERS(STRING,NAT,TMEMBER)
          ALLOCATE(VEC(3,1))
          CALL LINKEDLIST$GET(LL_CNSTR,'PHI',1,VEC)
          ALLOCATE(S0(NC1))
          ALLOCATE(SM(NC1))
          CALL CONSTRAINTS_ROT$VALUE(NAT,VEC(1,1),VEC(2,1),VEC(3,1) &
     &                                ,TMEMBER,R0_,RM_,DELT,RMASS_,S0(1))
          SM(1)=S0(1)
          DEALLOCATE(TMEMBER)
          DEALLOCATE(VEC)
        ENDIF
!
!         ==  ORIENTATION  =============================================
        STRING='ORIENTATION'
        IF(TYPE.EQ.STRING(1:32)) THEN
          NC1=1
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
          ALLOCATE(TMEMBER(NAT_))
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP',1,STRING)
          CALL GROUPLIST$MEMBERS(STRING,NAT,TMEMBER)
          CALL LINKEDLIST$GET(LL_CNSTR,'PHI',1,PHI)
          ALLOCATE(VEC(3,NAT))
          CALL CONSTRAINTS_ORIENTATIONVECTOR(NAT,TMEMBER,PHI,R0_,RMASS_,VEC)
          DEALLOCATE(TMEMBER)
!         == CONVERT INTO GENERALLINEAR ==============================
          TYPE='GENERALLINEAR'
          CALL LINKEDLIST$SET(LL_CNSTR,'TYPE',0,TYPE)
          CALL LINKEDLIST$SET(LL_CNSTR,'VEC',0,VEC)
          NC1=1
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
          ALLOCATE(S0(NC1))
          ALLOCATE(SM(NC1))
          CALL CONSTRAINTS_GLCVALUE(NAT,VEC,R0_,S0(1))
          CALL CONSTRAINTS_GLCVALUE(NAT,VEC,RM_,SM(1))
          DEALLOCATE(VEC)
        ENDIF
!
!       ==  MIDPLANE==================================================
        STRING='MIDPLANE'
        IF(TYPE.EQ.STRING(1:32)) THEN
          NC1=1
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM1',1,STRING)
          CALL ATOMLIST$INDEXX(STRING,IAT1,IT1)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM2',1,STRING)
          CALL ATOMLIST$INDEXX(STRING,IAT2,IT2)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM3',1,STRING)
          CALL ATOMLIST$INDEXX(STRING,IAT3,IT3)
!
          ALLOCATE(S0(NC1))
          ALLOCATE(SM(NC1))
          CALL CONSTRAINTS_MIDPLANEVALUE(NAT,IAT1,IT1,IAT2,IT2,IAT3,IT3 &
      &                                 ,RBAS,R0_,S0(1))
          CALL CONSTRAINTS_MIDPLANEVALUE(NAT,IAT1,IT1,IAT2,IT2,IAT3,IT3 &
      &                                 ,RBAS,RM_,SM(1))
        END IF

!       ==  COGSEP    ===============================================
        STRING='COGSEP'
        IF(TYPE.EQ.STRING(1:32)) THEN
          NC1=1
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP1',1,STRING)
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP2',1,STRING1)
!  
          ALLOCATE(S0(NC1))
          ALLOCATE(SM(NC1))
          CALL CONSTRAINTS_COGSEPVALUE(NAT,STRING,STRING1,R0_,S0(1))
          CALL CONSTRAINTS_COGSEPVALUE(NAT,STRING,STRING1,RM_,SM(1))
        END IF
!  
!       ==  BONDANGLE  ===============================================
        STRING='BONDANGLE'
        IF(TYPE.EQ.STRING(1:32)) THEN
          NC1=1
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM1',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT1)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM2',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT2)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM3',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT3)
!  
          ALLOCATE(S0(NC1))
          ALLOCATE(SM(NC1))
          CALL CONSTRAINTS_BONDANGLEVALUE(NAT,IAT1,IAT2,IAT3,R0_,S0(1))
          CALL CONSTRAINTS_BONDANGLEVALUE(NAT,IAT1,IAT2,IAT3,RM_,SM(1))
        END IF
   
!       ==  TORSION =================================================
        STRING='TORSION'
        IF(TYPE.EQ.STRING(1:32)) THEN
          NC1=1
          CALL LINKEDLIST$SET(LL_CNSTR,'NC1',0,NC1)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM1',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT1)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM2',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT2)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM3',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT3)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM4',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT4)
!  
          ALLOCATE(S0(NC1))
          ALLOCATE(SM(NC1))
          CALL CONSTRAINTS_TORSIONVALUE(NAT,IAT1,IAT2,IAT3,IAT4 &
   &                                   ,R0_,S0(1))
          CALL CONSTRAINTS_TORSIONVALUE(NAT,IAT1,IAT2,IAT3,IAT4 &
   &                                   ,RM_,SM(1))
        END IF
!
!       ==  GENERAL PROPERTIES  ======================================
        CALL LINKEDLIST$EXISTD(LL_CNSTR,'S0',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$EXISTD(LL_CNSTR,'SM',1,TCHK1)
          IF(.NOT.TCHK1) THEN
            CALL LINKEDLIST$SET(LL_CNSTR,'SM',0,S0)
          END IF
        ELSE
          CALL LINKEDLIST$SET(LL_CNSTR,'S0',0,S0)
          CALL LINKEDLIST$SET(LL_CNSTR,'SM',0,SM)
        END IF         
!
        CALL LINKEDLIST$EXISTD(LL_CNSTR,'FLOAT',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$EXISTD(LL_CNSTR,'ANNES',1,TCHK1)
          IF(.NOT.TCHK1)THEN
            ANNES=0.D0
            CALL LINKEDLIST$SET(LL_CNSTR,'ANNES',0,ANNES)
          END IF
          CALL LINKEDLIST$EXISTD(LL_CNSTR,'SMASS',1,TCHK1)
          IF(.NOT.TCHK1)THEN
            CALL CONSTANTS('U',SMASS)
            CALL LINKEDLIST$SET(LL_CNSTR,'SMASS',0,SMASS)
          END IF
        ELSE
          CALL LINKEDLIST$SET(LL_CNSTR,'FLOAT',0,.FALSE.)
        END IF
!
!       ==============================================================
!       ==  CHECK PARAMETERS FOR MOVING CONSTRAINTS                 ==
!       ==============================================================
        CALL LINKEDLIST$EXISTD(LL_CNSTR,'MOVE',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNSTR,'MOVE',0,.FALSE.)
        CALL LINKEDLIST$GET(LL_CNSTR,'MOVE',1,TCHK)
        IF(TCHK) THEN
!         == CHECK FINAL VALUE OF THE CONSTRAINT =====================
          CALL LINKEDLIST$EXISTD(LL_CNSTR,'SFINAL',1,TCHK1)
          SETS=TCHK1
          CALL LINKEDLIST$SET(LL_CNSTR,'SETS',0,SETS)
          IF(SETS) THEN
            CALL LINKEDLIST$SIZE(LL_CNSTR,'SFINAL',1,NBYTE)
            IF(NBYTE.NE.NC1) THEN
              CALL ERROR$MSG('SIZE OF SFINAL INCORRECT')
              CALL ERROR$CHVAL('TYPE',TYPE)
              CALL ERROR$STOP('CONSTRAINTS$SETREFERENCE')
            END IF
          ELSE
            CALL LINKEDLIST$SET(LL_CNSTR,'SFINAL',0,S0)
          END IF
!           == CHECK FINAL VELOCITY OF THE CONSTRAINT ==================
          CALL LINKEDLIST$EXISTD(LL_CNSTR,'VFINAL',1,TCHK1)
          SETV=TCHK1
          CALL LINKEDLIST$SET(LL_CNSTR,'SETV',0,SETV)
          IF(SETV) THEN
            CALL LINKEDLIST$SIZE(LL_CNSTR,'VFINAL',1,NBYTE)
            IF(NBYTE.NE.NC1) THEN
              CALL ERROR$MSG('SIZE OF VFINAL INCORRECT')
              CALL ERROR$CHVAL('TYPE',TYPE)
              CALL ERROR$STOP('CONSTRAINTS$SETREFERENCE')
            END IF
          ELSE
            ALLOCATE(VF(NC1))
            DO II=1,NC1
              VF(II)=0.D0
            ENDDO
            CALL LINKEDLIST$SET(LL_CNSTR,'VFINAL',0,VF)
            DEALLOCATE(VF)
          END IF
          CALL LINKEDLIST$EXISTD(LL_CNSTR,'NSTEP',1,TCHK1)
          IF(.NOT.TCHK1) THEN
            CALL ERROR$MSG('NSTEP MUST BE SET WHEN CONSTRAINTS MOVE')
            CALL ERROR$CHVAL('TYPE',TYPE)
            CALL ERROR$STOP('CONSTRAINTS$SETREFERENCE')
          END IF
        END IF
        DEALLOCATE(S0)
        DEALLOCATE(SM)
        IC=IC+NC1
      ENDDO
      DEALLOCATE(RREF)
      NC=IC
      CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS$NFREE(NAT_,NFREE_)
!     ******************************************************************
!     **  RETURN NUMBER OF CONSTRAINTS                                **
!     ******************************************************************
      USE CONSTRAINTS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT_
      INTEGER(4),INTENT(OUT):: NFREE_
      INTEGER(4)            :: NUM
      INTEGER(4)            :: I,IAT
      CHARACTER(32)         :: TYPE
      CHARACTER(32)         :: STRING
      INTEGER(4)            :: NC1
      INTEGER(4)            :: NAT
      INTEGER(4)            :: NMEMBER
      LOGICAL(4)            :: TMEMBER(NAT_)
!     ******************************************************************
      IF(.NOT.TON) THEN
        NFREE_=3*NAT_
        RETURN
      END IF
      IF(TSET) THEN
        NFREE_=3*NAT_-NC
        RETURN
      END IF
      NAT=NAT_
      CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
      CALL LINKEDLIST$NLISTS(LL_CNSTR,'CONSTRAINT',NUM)
      NC=0
      DO I=1,NUM
        NC1=0
        CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
        CALL LINKEDLIST$SELECT(LL_CNSTR,'CONSTRAINT',I)
        CALL LINKEDLIST$GET(LL_CNSTR,'TYPE',1,TYPE)
        STRING='BONDLENGTH'
        IF(TYPE.EQ.STRING(1:32)) NC1=1
        STRING='GENERALLINEAR'
        IF(TYPE.EQ.STRING(1:32)) NC1=1
        STRING='FIXATOM'
        IF(TYPE.EQ.STRING(1:32)) NC1=3
        STRING='FIXGROUP'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP',1,STRING)
          CALL GROUPLIST$MEMBERS(STRING,NAT,TMEMBER)
          NC1=0
          DO IAT=1,NAT
            IF(TMEMBER(IAT))NC1=NC1+3
          ENDDO
        ENDIF
        STRING='RIGID'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP',1,STRING)
          CALL GROUPLIST$MEMBERS(STRING,NAT_,TMEMBER)
          NMEMBER=0
          DO IAT=1,NAT_
            IF(TMEMBER(IAT)) NMEMBER=NMEMBER+1
          ENDDO
          NC1=3*NMEMBER-6
        ENDIF
        STRING='TRANSLATIONALMOMENTUM'
        IF(TYPE.EQ.STRING(1:32)) NC1=1
        STRING='ANGULARMOMENTUM'
        IF(TYPE.EQ.STRING(1:32)) NC1=1
        STRING='ORIENTATION'
        IF(TYPE.EQ.STRING(1:32)) NC1=1
        STRING='MIDPLANE'
        IF(TYPE.EQ.STRING(1:32)) NC1=1
        STRING='COGSEP'
        IF(TYPE.EQ.STRING(1:32)) NC1=1
        NC=NC+NC1
      ENDDO
      NFREE_=3*NAT-NC
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS$SWITCH
!     ******************************************************************
!     **  STEP FORWARD                                                **
!     ******************************************************************
      USE CONSTRAINTS_MODULE
      IMPLICIT NONE
      INTEGER(4)            :: ICOUNT
      INTEGER(4)            :: I
      INTEGER(4)            :: NUM
      INTEGER(4)            :: NC1
      REAL(8)   ,ALLOCATABLE:: S0(:)  !(NC1)
!     ******************************************************************
      IF(.NOT.TON) RETURN
      CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
      CALL LINKEDLIST$GET(LL_CNSTR,'ICOUNT',1,ICOUNT)
      ICOUNT=ICOUNT+1
      CALL LINKEDLIST$SET(LL_CNSTR,'ICOUNT',0,ICOUNT)
!
      CALL LINKEDLIST$NLISTS(LL_CNSTR,'CONSTRAINT',NUM)
      DO I=1,NUM
        CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
        CALL LINKEDLIST$SELECT(LL_CNSTR,'CONSTRAINT',I)
        CALL LINKEDLIST$GET(LL_CNSTR,'NC1',1,NC1)
        ALLOCATE(S0(NC1))
        CALL LINKEDLIST$GET(LL_CNSTR,'S0',1,S0)
        CALL LINKEDLIST$SET(LL_CNSTR,'SM',1,S0)
        CALL LINKEDLIST$GET(LL_CNSTR,'SP',1,S0)
        CALL LINKEDLIST$SET(LL_CNSTR,'S0',1,S0)
        DEALLOCATE(S0)
      ENDDO
      RETURN
      END
!
!     .................................................................
      SUBROUTINE CONSTRAINTS$APPLY(RBAS,NAT,R0,RP,RMASS,FORCE,DELT)
!     ******************************************************************
!     **  APPLY CONSTRAINTS                                           **
!     ******************************************************************
      USE CONSTRAINTS_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: RBAS(3,3)
      INTEGER(4),INTENT(IN)   :: NAT
      REAL(8)   ,INTENT(IN)   :: RMASS(NAT)
      REAL(8)   ,INTENT(IN)   :: R0(3,NAT)
      REAL(8)   ,INTENT(IN)   :: DELT
      REAL(8)   ,INTENT(INOUT):: RP(3,NAT)
      REAL(8)   ,INTENT(INOUT):: FORCE(3,NAT)
      REAL(8)   ,ALLOCATABLE  :: RLAM(:)  !(NC)
      REAL(8)   ,ALLOCATABLE  :: A(:)     !(NC)
      REAL(8)   ,ALLOCATABLE  :: B(:,:)   !(3*NAT,NC)
      REAL(8)   ,ALLOCATABLE  :: C(:,:,:) !(3*NAT,3*NAT,NC)
      REAL(8)                 :: VEC(3,NAT)
      INTEGER(4)              :: I,IAT
      REAL(8)                 :: EKINS
      LOGICAL(4)              :: TSTRESS
!     ******************************************************************
      IF(.NOT.TON) RETURN
                            CALL TIMING$CLOCKON('CONSTRAINTS')
                            CALL TRACE$PUSH('CONSTRAINTS$APPLY')
      IF(.NOT.TSET) THEN
        CALL ERROR$MSG('CONSTRAINTS OBJECT NOT INITIALIZED')
        CALL ERROR$STOP('CONSTRAINTS$APPLY')
      END IF

!     ==========================================================================
!     == THE STRESSES \SUM_I R_I\OTIMES F_I FROM THE CONSTRAINT FORCES ARE NOT==
!     == CONSIDERED FOR THE UNIT-CELL PROPAGATION. ATOMS$CONSTRAINTS NEED     ==
!     == TO BE CALLED BEFORE CELL$PROPAGATE                                   ==
!     ==========================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
      IF(TSTRESS) THEN
        CALL ERROR$MSG('ATOM-CONSTRAINTS AND CELL DYNAMICS ARE INCOMPATIBLE')
        CALL ERROR$MSG('STRESS FROM CONSTRAINT FORCES ARE NOT YET CONSIDERED')
        CALL ERROR$MSG('THIS IS AN IMPLEMENTATION ISSUE')
        CALL ERROR$STOP('ATOMS$CONSTRAINTS')
     END IF
!     
!     ==========================================================================
!     == COLLECT CONSTRAINTS                                                  ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
      IF(NC.EQ.0) RETURN
      ALLOCATE(A(NC))
      ALLOCATE(B(3*NAT,NC))
      ALLOCATE(C(3*NAT,3*NAT,NC))
      CALL CONSTRAINTS_APPLY_COLLECT(LL_CNSTR,NC,NAT,DELT,RBAS,RMASS,R0,A,B,C)
!     
!     ================================================================
!     == APPLY CONSTRAINTS                                          ==
!     ================================================================
      ALLOCATE(RLAM(NC))
      DO IAT=1,NAT
        DO I=1,3
          VEC(I,IAT)=RMASS(IAT)
        ENDDO
      ENDDO
      CALL CONSTRAINTS_APPLY(3*NAT,R0,RP,VEC,DELT,NC,A,B,C,RLAM)
      DEALLOCATE(A)
      DEALLOCATE(B)
      DEALLOCATE(C)
!     
!     ================================================================
!     == CALCULATE FORCE ON S AND PROPAGATE S                       ==
!     ================================================================
      CALL CONSTRAINTS_APPLY_MOVE(LL_CNSTR,DELT,EKINS,NC,RLAM)
      CALL ENERGYLIST$SET('CONSTRAINT KINETIC ENERGY',EKINS)
      DEALLOCATE(RLAM)
                            CALL TIMING$CLOCKOFF('CONSTRAINTS')
                            CALL TRACE$POP
      RETURN
      END
!     ..................................................................
      SUBROUTINE CONSTRAINTS_APPLY_MOVE(LL_CNSTR_,DELT,EKINS,NC,RLAM)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNSTR_
      INTEGER(4)   ,INTENT(IN) :: NC
      REAL(8)      ,INTENT(IN) :: DELT
      REAL(8)      ,INTENT(IN) :: RLAM(NC)
      REAL(8)      ,INTENT(OUT):: EKINS
      TYPE(LL_TYPE)            :: LL_CNSTR
      INTEGER(4)               :: NUM
      INTEGER(4)               :: I,IC,II
      INTEGER(4)               :: NC1
      REAL(8)     ,ALLOCATABLE :: S0(:) !NC1
      REAL(8)     ,ALLOCATABLE :: SM(:) !NC1
      REAL(8)     ,ALLOCATABLE :: SP(:) !NC1
      REAL(8)     ,ALLOCATABLE :: SF(:) !NC1
      REAL(8)     ,ALLOCATABLE :: VF(:) !NC1
      LOGICAL(4)               :: TFLOAT
      LOGICAL(4)               :: TMOVE
      REAL(8)                  :: ANNES
      REAL(8)                  :: SMASS
      REAL(8)                  :: DEKIN
      LOGICAL(4)               :: SETS,SETV
      INTEGER(4)               :: NSTEP      
      INTEGER(4)               :: ICOUNT
!     ******************************************************************
      EKINS=0.D0
      LL_CNSTR=LL_CNSTR_
      CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
      CALL LINKEDLIST$GET(LL_CNSTR,'ICOUNT',1,ICOUNT)
      CALL LINKEDLIST$NLISTS(LL_CNSTR,'CONSTRAINT',NUM)
      DO I=1,NUM
        CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
        CALL LINKEDLIST$SELECT(LL_CNSTR,'CONSTRAINT',I)
        CALL LINKEDLIST$GET(LL_CNSTR,'IC',1,IC)
        CALL LINKEDLIST$GET(LL_CNSTR,'NC1',1,NC1)
        CALL LINKEDLIST$SET(LL_CNSTR,'FC',0,RLAM(IC:IC+NC1-1))
        ALLOCATE(S0(NC1))
        ALLOCATE(SM(NC1))
        ALLOCATE(SP(NC1))
        CALL LINKEDLIST$GET(LL_CNSTR,'S0',1,S0)
        CALL LINKEDLIST$GET(LL_CNSTR,'SM',1,SM)
        CALL LINKEDLIST$GET(LL_CNSTR,'FLOAT',1,TFLOAT)
        CALL LINKEDLIST$GET(LL_CNSTR,'MOVE',1,TMOVE)
        IF(TFLOAT) THEN          
          CALL LINKEDLIST$GET(LL_CNSTR,'ANNES',1,ANNES)
          CALL LINKEDLIST$GET(LL_CNSTR,'SMASS',1,SMASS)
          DO II=1,NC1
            CALL CONSTRAINTS_FLOATS(S0(II),SM(II),SP(II),DELT &
     &                             ,SMASS,ANNES,-RLAM(IC+II-1),DEKIN)
          ENDDO
          EKINS=EKINS+DEKIN
          CALL LINKEDLIST$SET(LL_CNSTR,'SP',0,SP)
        ELSE IF(TMOVE) THEN          
          CALL LINKEDLIST$GET(LL_CNSTR,'NSTEP',1,NSTEP)
          CALL LINKEDLIST$GET(LL_CNSTR,'SETS',1,SETS)
          CALL LINKEDLIST$GET(LL_CNSTR,'SETV',1,SETV)
          ALLOCATE(SF(NC1))
          ALLOCATE(VF(NC1))
          CALL LINKEDLIST$GET(LL_CNSTR,'SFINAL',1,SF)
          CALL LINKEDLIST$GET(LL_CNSTR,'VFINAL',1,VF)
          DO II=1,NC1
            CALL CONSTRAINTS_MOVES(ICOUNT,S0(II),SM(II),SP(II) &
     &                      ,DELT,NSTEP,SETS,SF(II),SETV,VF(II))
          ENDDO
          CALL LINKEDLIST$SET(LL_CNSTR,'SP',0,SP)
          DEALLOCATE(SF)
          DEALLOCATE(VF)
        ELSE
          CALL LINKEDLIST$SET(LL_CNSTR,'SP',0,S0)
        END IF
        DEALLOCATE(S0)
        DEALLOCATE(SM)
        DEALLOCATE(SP)
      ENDDO
      RETURN
      END
!
!     .................................................................
      SUBROUTINE CONSTRAINTS_APPLY_COLLECT(LL_CNSTR_,NC,NAT,DELT,RBAS &
     &                                    ,RMASS,R0,A,B,C)
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN)  :: LL_CNSTR_
      INTEGER(4)   ,INTENT(IN)  :: NC
      INTEGER(4)   ,INTENT(IN)  :: NAT
      REAL(8)      ,INTENT(IN)  :: DELT
      REAL(8)      ,INTENT(IN)  :: RBAS(3,3)
      REAL(8)      ,INTENT(IN)  :: R0(3,NAT)
      REAL(8)      ,INTENT(IN)  :: RMASS(NAT)
      REAL(8)      ,INTENT(OUT) :: A(NC)  
      REAL(8)      ,INTENT(OUT) :: B(3*NAT,NC)  
      REAL(8)      ,INTENT(OUT) :: C(3*NAT,3*NAT,NC)  
      TYPE(LL_TYPE)             :: LL_CNSTR
      LOGICAL(4)                :: TMEMBER(NAT)
      INTEGER(4)                :: NC1
      INTEGER(4)                :: NUM
      CHARACTER(32)             :: TYPE
      CHARACTER(32)             :: STRING,STRING1
      REAL(8)                   :: DIR(3)
      REAL(8)                   :: VEC(3,NAT)
      REAL(8)      ,ALLOCATABLE :: S0(:)            !(NC1)
      REAL(8)                   :: RREF(3,NAT)
      INTEGER(4)                :: IC,I,IAT1,IAT2,IAT3,IAT4,IAT,II,ICI
      INTEGER(4)                :: IT1(3),IT2(3),IT3(3)
!     ******************************************************************
      A(:)=0.D0
      B(:,:)=0.D0
      C(:,:,:)=0.D0
      LL_CNSTR=LL_CNSTR_
      CALL LINKEDLIST$GET(LL_CNSTR,'REFERENCE',1,RREF)
!     
      CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
      CALL LINKEDLIST$NLISTS(LL_CNSTR,'CONSTRAINT',NUM)
      IC=1
      DO I=1,NUM
        CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
        CALL LINKEDLIST$SELECT(LL_CNSTR,'CONSTRAINT',I)
        CALL LINKEDLIST$GET(LL_CNSTR,'TYPE',1,TYPE)
        CALL LINKEDLIST$SET(LL_CNSTR,'IC',0,IC)
        CALL LINKEDLIST$GET(LL_CNSTR,'NC1',1,NC1)
!     
!     ==  BONDLENGTH =================================================
        STRING='BONDLENGTH'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM1',1,STRING)
          CALL ATOMLIST$INDEXX(STRING,IAT1,IT1)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM2',1,STRING)
          CALL ATOMLIST$INDEXX(STRING,IAT2,IT2)
          CALL CONSTRAINTS_SETBONDLENGTH(NAT,IAT1,IT1,IAT2,IT2 &
     &                  ,A(IC),B(1,IC),C(1,1,IC),RBAS,R0)
        ENDIF
!     
!       ==  TRANSLATIONALMOMENTUM  ===================================
        STRING='TRANSLATIONALMOMENTUM'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP',1,STRING)
          CALL GROUPLIST$MEMBERS(STRING,NAT,TMEMBER)
          CALL LINKEDLIST$GET(LL_CNSTR,'DIR',1,DIR)
          CALL CONSTRAINTS_SETCOG(NAT,DIR(1),DIR(2),DIR(3) &
     &                    ,TMEMBER,RMASS,A(IC),B(1,IC),C(1,1,IC))
        ENDIF
!     
!       ==  ANGULARMOMENTUM ==========================================
        STRING='ANGULARMOMENTUM'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP',1,STRING)
          CALL GROUPLIST$MEMBERS(STRING,NAT,TMEMBER)
          CALL LINKEDLIST$GET(LL_CNSTR,'PHI',1,DIR)
          CALL CONSTRAINTS_ROT$SET &
     &         (NAT,DIR(1),DIR(2),DIR(3),TMEMBER &
     &         ,R0,DELT,RMASS,A(IC),B(1,IC),C(1,1,IC))
        ENDIF
!     
!       ==  RIGID  ===================================================
        STRING='RIGID'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP',1,STRING)
          CALL GROUPLIST$MEMBERS(STRING,NAT,TMEMBER)
          CALL LINKEDLIST$GET(LL_CNSTR,'IAT1',1,IAT1)
          CALL LINKEDLIST$GET(LL_CNSTR,'IAT2',1,IAT2)
          CALL LINKEDLIST$GET(LL_CNSTR,'IAT3',1,IAT3)
          CALL CONSTRAINTS_RIGID$SET(NAT,TMEMBER,R0,IAT1,IAT2,IAT3 &
     &                     ,NC1,A(IC),B(1,IC),C(1,1,IC))
        ENDIF
!     
!       ==  FIXGROUP  ================================================
        STRING='FIXGROUP'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP',1,STRING)
          CALL GROUPLIST$MEMBERS(STRING,NAT,TMEMBER)
          II=IC
          DO IAT=1,NAT
            IF(TMEMBER(IAT)) THEN
              CALL CONSTRAINTS_SETFIXATOM(NAT,IAT &
     &                       ,A(II),B(1,II),C(1,1,II))
              II=II+3
            END IF
          ENDDO
        ENDIF
!     
!       ==  FIXATOM   ================================================
        STRING='FIXATOM'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT)
          CALL CONSTRAINTS_SETFIXATOM(NAT,IAT,A(IC),B(1,IC),C(1,1,IC))
        ENDIF
!     
!       ==  GENERAL LINEAR ===========================================
        STRING='GENERALLINEAR'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'VEC',1,VEC)
          CALL CONSTRAINTS_SETGLC(NAT,VEC,A(IC),B(1,IC),C(1,1,IC))
        ENDIF
!     
!       ==  MIDPLANE==================================================
        STRING='MIDPLANE'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM1',1,STRING)
          CALL ATOMLIST$INDEXX(STRING,IAT1,IT1)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM2',1,STRING)
          CALL ATOMLIST$INDEXX(STRING,IAT2,IT2)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM3',1,STRING)
          CALL ATOMLIST$INDEXX(STRING,IAT3,IT3)
          CALL CONSTRAINTS_SETMIDPLANE(NAT,IAT1,IT1,IAT2,IT2,IAT3,IT3 &
     &                              ,A(IC),B(1,IC),C(1,1,IC),RBAS,R0)
        ENDIF
!
!       ==  COG SEPARATION===============================================
        STRING='COGSEP'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP1',1,STRING)
          CALL LINKEDLIST$GET(LL_CNSTR,'GROUP2',1,STRING1)
          CALL CONSTRAINTS_SETCOGSEP(NAT,STRING,STRING1 &
     &                  ,A(IC),B(1,IC),C(1,1,IC),R0)
        ENDIF
!
!       ==  BONDANGLE   =================================================
        STRING='BONDANGLE'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM1',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT1)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM2',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT2)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM3',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT3)
          CALL CONSTRAINTS_SETBONDANGLE(NAT,IAT1,IAT2,IAT3 &
     &                  ,A(IC),B(1,IC),C(1,1,IC),R0)
        ENDIF
!
!       ==  TORSION   ===================================================
        STRING='TORSION'
        IF(TYPE.EQ.STRING(1:32)) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM1',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT1)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM2',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT2)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM3',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT3)
          CALL LINKEDLIST$GET(LL_CNSTR,'ATOM4',1,STRING)
          CALL ATOMLIST$INDEX(STRING,IAT4)
!       
          CALL CONSTRAINTS_SETTORSION(NAT,IAT1,IAT2,IAT3,IAT4 &
     &                  ,A(IC),B(1,IC),C(1,1,IC),R0)
        ENDIF

!     
        ALLOCATE(S0(NC1))
        CALL LINKEDLIST$GET(LL_CNSTR,'S0',1,S0)
        DO II=1,NC1
          ICI=IC+II-1
          A(ICI)=A(ICI)-S0(II)
          CALL CONSTRAINTS_TESTSET(3*NAT,R0 &
     &                    ,A(ICI),B(1,ICI),C(1,1,ICI),1.D-3)
        ENDDO
        DEALLOCATE(S0)
        IC=IC+NC1
      ENDDO
!     
      IF(IC.NE.NC+1) THEN
        CALL ERROR$MSG('NUMBER OF CONSTRAINTS INCONSISTENT')
        CALL ERROR$I4VAL('IC',IC)
        CALL ERROR$I4VAL('NC',NC)
        CALL ERROR$STOP('CONSTRAINTS$APPLY')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS$REPORT(NFIL,IDENT_)
!     ******************************************************************
!     **  RETURN NUMBER OF CONSTRAINTS                                **
!     ******************************************************************
      USE CONSTRAINTS_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NFIL
      CHARACTER(*) ,INTENT(IN) :: IDENT_
      INTEGER(4)               :: NUM
      INTEGER(4)               :: I
      INTEGER(4)               :: II
      INTEGER(4)               :: NC1
      CHARACTER(32)            :: TYPE
      CHARACTER(256)           :: STRING
      LOGICAL(4)               :: TCHK
      LOGICAL(4)               :: TLONG
      REAL(8)      ,ALLOCATABLE:: S0(:)
      REAL(8)      ,ALLOCATABLE:: SM(:)
!     ******************************************************************
      IF(.NOT.TON) RETURN
      IF(.NOT.TSET) RETURN
      TLONG=(IDENT_.NE.'SHORT')
      IF(TLONG) THEN
        CALL REPORT$TITLE(NFIL,'CONSTRAINTS')
      END IF
      CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
      CALL LINKEDLIST$NLISTS(LL_CNSTR,'CONSTRAINT',NUM)
      IF(NUM.EQ.0) THEN
        WRITE(NFIL,FMT='("NO CONSTRAINTS SPECIFIED")')
        RETURN
      END IF
      DO I=1,NUM
        CALL LINKEDLIST$SELECT(LL_CNSTR,'~')
        CALL LINKEDLIST$SELECT(LL_CNSTR,'CONSTRAINT',I)
        CALL LINKEDLIST$GET(LL_CNSTR,'TYPE',1,TYPE)
        CALL LINKEDLIST$GET(LL_CNSTR,'TITLE',1,STRING)
        IF(.NOT.TLONG) THEN
          CALL LINKEDLIST$GET(LL_CNSTR,'SHOW',1,TCHK)
          CALL LINKEDLIST$GET(LL_CNSTR,'NC1',1,NC1)
          IF(TSET.AND.TCHK) THEN
            ALLOCATE(S0(NC1))
            ALLOCATE(SM(NC1))
            CALL LINKEDLIST$GET(LL_CNSTR,'S0',1,S0)
            CALL LINKEDLIST$GET(LL_CNSTR,'FC',1,SM)
            WRITE(NFIL,ADVANCE='NO',FMT='("!>",I5)')I
            DO II=1,NC1
              IF(II.NE.NC1) THEN
                WRITE(NFIL,ADVANCE='NO',FMT='("(",2ES10.4E1,")")') &
                         S0(II),-SM(II)
              ELSE
                WRITE(NFIL,ADVANCE='YES',FMT='("(",2ES10.4E1,")")') &
                         S0(II),-SM(II)
              END IF
            ENDDO
            DEALLOCATE(SM)
            DEALLOCATE(S0)
          END IF          
        ELSE
          CALL LINKEDLIST$GET(LL_CNSTR,'SHOW',1,TCHK)
          WRITE(NFIL,FMT='(A)')TRIM(STRING)
          IF(TSET.AND.TCHK) THEN
            CALL LINKEDLIST$GET(LL_CNSTR,'NC1',1,NC1)
            ALLOCATE(S0(NC1))
            ALLOCATE(SM(NC1))
            CALL LINKEDLIST$GET(LL_CNSTR,'S0',1,S0)
            CALL LINKEDLIST$GET(LL_CNSTR,'FC',1,SM)
            WRITE(NFIL,FMT='("(VAL/FORCE):"' &
     &                    //',5("(",ES10.4E1,"/",E10.4E1,")"))') &
     &              (S0(II),-SM(II),II=1,NC1)
            DEALLOCATE(SM)
            DEALLOCATE(S0)
          END IF
        END IF
      ENDDO        
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_FLOATS(S0,SM,SP,DELT,SMASS,ANNES,FS,EKINS)
      IMPLICIT NONE
      REAL(8)    ,INTENT(IN) :: S0
      REAL(8)    ,INTENT(IN) :: SM
      REAL(8)    ,INTENT(OUT):: SP
      REAL(8)    ,INTENT(OUT):: EKINS
      REAL(8)    ,INTENT(IN) :: ANNES
      REAL(8)    ,INTENT(IN) :: FS
      REAL(8)    ,INTENT(IN) :: DELT
      REAL(8)    ,INTENT(IN) :: SMASS
      INTEGER(4)             :: NFILO
      REAL(8)                :: SVAR1,SVAR2,SVAR3
!     ******************************************************************
      SVAR1=2.D0/(1.D0+ANNES)
      SVAR2=1.D0-SVAR1
      SVAR3=DELT**2/SMASS/(1.D0+ANNES)
      SP=SVAR1*S0+SVAR2*SM+SVAR3*FS
      EKINS=0.5D0*SMASS*((SP-SM)/(2.D0*DELT))**2
!
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,*)'S0,SM,SP,FS',S0,SM,SP,FS
      WRITE(NFILO,*)'SMASS,ANNES,DELT',SMASS,ANNES,DELT
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_MOVES(ICOUNT,S0,SM,SP,DELTA,NSTEP &
     &                              ,SETS,SF_,SETV,VF_)
!     ** 
!     **  MOVES CONSTRAINT TO A GIVEN VALUE OR VELOCITY               **
!     ** 
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT) :: ICOUNT
      REAL(8)   ,INTENT(IN)    :: S0
      REAL(8)   ,INTENT(IN)    :: SM
      REAL(8)   ,INTENT(INOUT) :: SP
      REAL(8)   ,INTENT(IN)    :: DELTA
      REAL(8)   ,INTENT(IN)    :: SF_
      REAL(8)   ,INTENT(IN)    :: VF_
      INTEGER(4),INTENT(IN)    :: NSTEP
      LOGICAL(4),INTENT(IN)    :: SETS
      LOGICAL(4),INTENT(IN)    :: SETV
      REAL(8)                  :: TBYDELTA
      REAL(8)                  :: T
      REAL(8)                  :: VF,SF
!     ******************************************************************
      TBYDELTA=DBLE(ICOUNT-NSTEP)
!
!     ==================================================================
!     ==  AFTER ACCELERATION PHASE                                    ==
!     ==================================================================
      IF(TBYDELTA.GE.0.D0) THEN
        IF(SETS) THEN
          IF(SETV) THEN
            SP=SF_+DELTA*VF_*DBLE(ICOUNT-NSTEP)
          ELSE
            SP=SF_
          END IF
        ELSE
          IF(SETV) THEN
            SP=S0+VF_*DELTA
          ELSE
            SP=2.D0*S0-SM
          END IF
        END IF
        RETURN
      END IF
!
!     ==================================================================
!     ==  DURING ACCELERATION PHASE                                   ==
!     ==  THE METHOD IS TO ADD A CONTRIBUTION THAT DOES NOT CONFLICT  ==
!     ==  WITH CONSTRAINTS IMPOSED FORMERLY AND TO CORRECT THE FINAL  ==   
!     ==  VALUES RELATIVE TO THAT TRAJECTORY                          ==
!     ==================================================================
      T=TBYDELTA*DELTA
      SF=SF_
      VF=VF_
      IF(SETS.AND.(.NOT.SETV)) VF=0.D0
!
!     ==================================================================
!     ==  CONSTANT PROPAGATION  WITHOUT PREDEFINED FINAL VALUES       ==
!     ==================================================================
!     == SP=S0+(S0-SM)/DELTA*(TAU-T)
      SP=2.D0*S0-SM
      SF=SF-(S0-(S0-SM)/DELTA*T)
      VF=VF-(S0-SM)/DELTA
!     
!     ==================================================================
!     ==  ENFORCE FINAL VALUE OF THE CONSTRAINT                       ==
!     ==================================================================
      IF(SETS) THEN
!       ==  DELTA(SP)=SF/(-T)*(TAU-T)    TAU(0)=T 
        SP=SP-SF/T*DELTA
        SF=SF-SF
        VF=VF+SF/T
      END IF
!     
!     ==================================================================
!     ==  ENFORCE FINAL VELOCITY OF THE CONSTRAINT                    ==
!     ==================================================================
      IF(SETV) THEN
        IF(SETS) THEN
!         == DELTA(SP)=VF*TAU*(TAU-T)/(-T)  =============================
          SP=SP-VF*(T+DELTA)*DELTA/T
        ELSE
!         == DELTA(SP)=VF*(TAU-T)**2/(-2*T)  =============================
          SP=SP-VF*DELTA**2/(2.D0*T)
        END IF      
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_BONDLENGTHVALUE(NAT,IAT1,IT1,IAT2,IT2,RBAS,R0,VALUE)
!     ******************************************************************
!     **                                                              **
!     **  BOND-LENGTH CONSTRAINT : SET VALUE                          **   
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: IAT1
      INTEGER(4),INTENT(IN) :: IAT2
      INTEGER(4),INTENT(IN) :: IT1(3)
      INTEGER(4),INTENT(IN) :: IT2(3)
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)
      REAL(8)   ,INTENT(OUT):: VALUE
      REAL(8)               :: T2(3)
      INTEGER(4)            :: I
!     ******************************************************************
      T2=MATMUL(RBAS,REAL(IT2-IT1,KIND=8))
      VALUE=0.D0
      DO I=1,3
        VALUE=VALUE+(R0(I,IAT1)-(R0(I,IAT2)+T2(I)))**2
      ENDDO
      VALUE=SQRT(VALUE)
      RETURN      
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_SETBONDLENGTH(NAT,IAT1,IT1,IAT2,IT2,A,B,C,RBAS,R0)
!     ******************************************************************
!     **                                                              **
!     **  BOND-LENGTH CONSTRAINT : SET UP COORDINATES                 **   
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: IAT1
      INTEGER(4),INTENT(IN) :: IT1(3)
      INTEGER(4),INTENT(IN) :: IAT2
      INTEGER(4),INTENT(IN) :: IT2(3)
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)
      REAL(8)   ,INTENT(OUT):: A
      REAL(8)   ,INTENT(OUT):: B(3,NAT)
      REAL(8)   ,INTENT(OUT):: C(3,NAT,3,NAT)
      INTEGER(4)            :: I,J
      REAL(8)               :: DIS,DI,DJ
      REAL(8)               :: DELTAK
      REAL(8)               :: SVAR
      REAL(8)               :: DR(3)
!     ******************************************************************
      CALL CONSTRAINTS_INIT(3*NAT,A,B,C) 
      DR(:)=R0(:,IAT2)-R0(:,IAT1)+MATMUL(RBAS,REAL(IT2(:)-IT1(:),KIND=8))
!
!     ==================================================================
!     == CALCULATE VALUE AND FIRST TWO DERIVATIVES OF D=|R_2-R_2|     ==
!     ==================================================================
      DIS=SQRT(DOT_PRODUCT(DR,DR))
      A=DIS
      DO I=1,3
        DI=DR(I)/DIS
        B(I,IAT2)=DI
        B(I,IAT1)=-DI
        DO J=1,3
          DELTAK=0.D0
          IF(I.EQ.J) DELTAK=1.D0
          DJ=DR(J)/DIS
          SVAR=(DELTAK-DI*DJ)/DIS
          C(I,IAT1,J,IAT1)=+SVAR
          C(I,IAT1,J,IAT2)=-SVAR
          C(I,IAT2,J,IAT1)=-SVAR
          C(I,IAT2,J,IAT2)=+SVAR
        ENDDO
      ENDDO
!
!     ==================================================================
!     == TRANSFORM INTO QUADRATIC EQUATION FOR R: A+B*R+0.5*R*C*R=0   ==
!     ==================================================================
      DO I=1,3
        A=A-R0(I,IAT1)*B(I,IAT1)-R0(I,IAT2)*B(I,IAT2)
        DO J=1,3
          A=A+0.5D0*(R0(I,IAT1)*(C(I,IAT1,J,IAT1)*R0(J,IAT1) &
     &                          +C(I,IAT1,J,IAT2)*R0(J,IAT2)) &
     &              +R0(I,IAT2)*(C(I,IAT2,J,IAT1)*R0(J,IAT1) &
     &                          +C(I,IAT2,J,IAT2)*R0(J,IAT2)))
        ENDDO
      ENDDO
      DO I=1,3
        DO J=1,3
          B(I,IAT1)=B(I,IAT1)-C(I,IAT1,J,IAT1)*R0(J,IAT1) &
     &                       -C(I,IAT1,J,IAT2)*R0(J,IAT2)
          B(I,IAT2)=B(I,IAT2)-C(I,IAT2,J,IAT1)*R0(J,IAT1) &
     &                       -C(I,IAT2,J,IAT2)*R0(J,IAT2)
        ENDDO
      ENDDO
!
!     CALL FILEHANDLER$UNIT('PROT',NFILO)
!     WRITE(NFILO,FMT='("BONDLENGTH CONSTRAINT",F10.5)')DIS
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_MIDPLANEVALUE(NAT,IAT1,IT1,IAT2,IT2,IAT3,IT3 &
     &                                    ,RBAS,R0,VALUE)
!     ******************************************************************
!     **  CONSTRAINTS_MIDPLANEVALUE                                   **
!     **  CALCULATE VALUE OF CONSTRAINT FROM THE REFERENCE STRUCTURE  **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT            ! #(ATOMS)
      INTEGER(4),INTENT(IN) :: IAT1           ! INDEX OF FIRST ATOM
      INTEGER(4),INTENT(IN) :: IT1(3)         ! #(LATTICE TRANSLATIONS FOR IAT1)
      INTEGER(4),INTENT(IN) :: IAT2           ! INDEX OF CENTRAL ATOM
      INTEGER(4),INTENT(IN) :: IT2(3)         ! #(LATTICE TRANSLATIONS FOR IAT2)
      INTEGER(4),INTENT(IN) :: IAT3           ! INDEX OF THIRD ATOM
      INTEGER(4),INTENT(IN) :: IT3(3)         ! #(LATTICE TRANSLATIONS FOR IAT3)
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)      ! LATTICE TRANSLATION VECTORS
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)      ! ATOMIC POSITIONS
      REAL(8)   ,INTENT(OUT):: VALUE          ! VALUE OF THE CONSTRAINT
      INTEGER(4)            :: I
      REAL(8)               :: R31R21,R31R31
      REAL(8)               :: R31(3),R21(3)
      REAL(8)               :: T31(3),T21(3)
!     ******************************************************************
      T21=MATMUL(RBAS,REAL(IT2-IT1,KIND=8))
      T31=MATMUL(RBAS,REAL(IT3-IT1,KIND=8))
      R31R21=0.D0
      R31R31=0.D0
      DO I=1,3
        R31(I)=R0(I,IAT3)-R0(I,IAT1)+T31(I)
        R21(I)=R0(I,IAT2)-R0(I,IAT1)+T21(I)
        R31R21=R31R21+R31(I)*R21(I)
        R31R31=R31R31+R31(I)**2
      ENDDO
      VALUE=R31R21/R31R31
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_SETMIDPLANE(NAT,IAT1,IT1,IAT2,IT2,IAT3,IT3 &
     &                                  ,A,B,C,RBAS,R0)
!     ******************************************************************
!     **  CONSTRAINTS_SETMIDPLANE                                     **
!     **  APPROXIMATES CONSTRAINT FUNCTION (WITHOUT VALUE)            **
!     **  BY AN QUADRATIC POLYNOM IN A+B*R+0.5*R*C*R                  **
!     **  THE CONSTRAINT IS THE POLYNOMIAL MINUS THE VALUE FROM THE   **
!     **  REFERENCE STRUCTURE                                         **
!     **                                                              **
!     **    X=(R2-R1)*(R3-R1)/|R3-R1|**2                               **
!     **                                                              **
!     **  A VALUE OF ZERO IMPLIES THAT R2 LIES IN A PLANE NORMAL TO   **
!     **  R3-R1 THROUGH R1                                            **
!     **  A VALUE OF ONE IMPLIES THAT R2 LIES IN A PLANE NORMAL TO    **
!     **  R3-R1 THROUGH R3                                            **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: IAT1
      INTEGER(4),INTENT(IN) :: IT1(3)
      INTEGER(4),INTENT(IN) :: IAT2
      INTEGER(4),INTENT(IN) :: IT2(3)
      INTEGER(4),INTENT(IN) :: IAT3
      INTEGER(4),INTENT(IN) :: IT3(3)
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)
      REAL(8)   ,INTENT(OUT):: A
      REAL(8)   ,INTENT(OUT):: B(3,NAT)
      REAL(8)   ,INTENT(OUT):: C(3,NAT,3,NAT)
      INTEGER(4)            :: I,J,I1,I2,IAT10,IAT20
      REAL(8)               :: R31R21,R31R31
      REAL(8)               :: R31(3),R21(3)
      REAL(8)               :: SVAR,SVAR1,SVAR2
      REAL(8)               :: DELTAK
      REAL(8)               :: T21(3),T31(3)
!     ******************************************************************
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      CALL CONSTRAINTS_INIT(3*NAT,A,B,C) 
      T21=MATMUL(RBAS,REAL(IT2-IT1,KIND=8))
      T31=MATMUL(RBAS,REAL(IT3-IT1,KIND=8))
!
!     ==================================================================
!     == CALCULATE TAYLOR EXPANSION COEFFICIENTS                      ==
!     ==================================================================
      R31R21=0.D0
      R31R31=0.D0
      DO I=1,3
        R31(I)=R0(I,IAT3)-R0(I,IAT1)+T31(I)
        R21(I)=R0(I,IAT2)-R0(I,IAT1)+T21(I)
        R31R21=R31R21+R31(I)*R21(I)
        R31R31=R31R31+R31(I)**2
      ENDDO
!
      SVAR1=1.D0/R31R31
      SVAR2=2.D0*R31R21/R31R31**2
      A=R31R21/R31R31
      DO I=1,3
        B(I,IAT1)=-SVAR1*(R21(I)+R31(I))+SVAR2*R31(I)
        B(I,IAT2)=+SVAR1        *R31(I)          
        B(I,IAT3)=+SVAR1 *R21(I)        -SVAR2*R31(I)
        DO J=1,3
          DELTAK=0.D0
          IF(I.EQ.J) DELTAK=1.D0
          SVAR=-DELTAK/R31R31+2.D0/R31R31**2*R31(I)*R31(J)
          C(I,IAT1,J,IAT1)=-2.D0*SVAR
          C(I,IAT1,J,IAT2)=SVAR
          C(I,IAT1,J,IAT3)=SVAR
          C(I,IAT2,J,IAT1)=SVAR
          C(I,IAT2,J,IAT2)=0.D0
          C(I,IAT2,J,IAT3)=-SVAR
          C(I,IAT3,J,IAT1)=SVAR
          C(I,IAT3,J,IAT2)=-SVAR
          C(I,IAT3,J,IAT3)=+SVAR2*DELTAK &
     &                    -2.D0*SVAR1**2*(R21(I)*R31(J)+R31(I)*R21(J)) &
     &                    +4.D0*R31(I)*R31(J)*SVAR1*SVAR2     
        ENDDO
      ENDDO
!
!     ==================================================================
!     == TRANSFORM INTO QUADRATIC EQUATION FOR R: A+B*R+0.5*R*C*R=0   ==
!     ==================================================================
      DO I=1,3
        A=A-R0(I,IAT1)*B(I,IAT1)-R0(I,IAT2)*B(I,IAT2) &
     &     -R0(I,IAT3)*B(I,IAT3)
        DO J=1,3
          A=A+0.5D0*(R0(I,IAT1)*(C(I,IAT1,J,IAT1)*R0(J,IAT1) &
     &                          +C(I,IAT1,J,IAT2)*R0(J,IAT2) &
     &                          +C(I,IAT1,J,IAT3)*R0(J,IAT3)) &
     &              +R0(I,IAT2)*(C(I,IAT2,J,IAT1)*R0(J,IAT1) &
     &                          +C(I,IAT2,J,IAT2)*R0(J,IAT2) &
     &                          +C(I,IAT2,J,IAT3)*R0(J,IAT3)) &
     &              +R0(I,IAT3)*(C(I,IAT3,J,IAT1)*R0(J,IAT1) &
     &                          +C(I,IAT3,J,IAT2)*R0(J,IAT2) &
     &                          +C(I,IAT3,J,IAT3)*R0(J,IAT3)))
        ENDDO
      ENDDO
      DO I=1,3
        DO J=1,3
          B(I,IAT1)=B(I,IAT1)-C(I,IAT1,J,IAT1)*R0(J,IAT1) &
     &                       -C(I,IAT1,J,IAT2)*R0(J,IAT2) &
     &                       -C(I,IAT1,J,IAT3)*R0(J,IAT3)
          B(I,IAT2)=B(I,IAT2)-C(I,IAT2,J,IAT1)*R0(J,IAT1) &
     &                       -C(I,IAT2,J,IAT2)*R0(J,IAT2) &
     &                       -C(I,IAT2,J,IAT3)*R0(J,IAT3)
          B(I,IAT3)=B(I,IAT3)-C(I,IAT3,J,IAT1)*R0(J,IAT1) &
     &                       -C(I,IAT3,J,IAT2)*R0(J,IAT2) &
     &                       -C(I,IAT3,J,IAT3)*R0(J,IAT3)
        ENDDO
      ENDDO
!
!     ==================================================================
!     == TEST                                                         ==
!     ==================================================================
      DO IAT10=1,NAT
        DO I1=1,3
          DO IAT20=IAT1,NAT
            DO I2=1,3
              SVAR=C(I1,IAT10,I2,IAT20)-C(I2,IAT20,I1,IAT10)
              IF(ABS(SVAR).GT.1.D-3) THEN
                WRITE(*,FMT='("CTEST",4I3,3F10.5)')I1,IAT10,I2,IAT20 &
     &               ,C(I1,IAT10,I2,IAT20),C(I2,IAT20,I1,IAT10),SVAR
              END IF
            ENDDO
          ENDDO
        ENDDO
      ENDDO  
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_GLCVALUE(NAT,VEC,R,VALUE)
!     ******************************************************************
!     **  CONSTRAINTS_GLCVALUE                                        **
!     **  CALCULATE VALUE OF CONSTRAINT FROM THE REFERENCE STRUCTURE  **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: VEC(3,NAT)
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      REAL(8)   ,INTENT(OUT):: VALUE
      INTEGER(4)            :: IAT,I
!     ******************************************************************
      VALUE=0.D0
      DO IAT=1,NAT
        DO I=1,3
          VALUE=VALUE+VEC(I,IAT)*R(I,IAT)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_SETGLC(NAT,VEC,A,B,C)
!     ******************************************************************
!     **  CONSTRAINTS_SETGLC                                          **
!     **  SET GENERAL LINEAR CONSTRAINT (GLC) CONST+VEC*R=0           **
!     **                                                              **
!     **  APPROXIMATES CONSTRAINT FUNCTION (WITHOUT VALUE)            **
!     **  BY AN QUADRATIC POLYNOM IN A+B*R+0.5*R*C*R                  **
!     **  THE CONSTRAINT IS THE POLYNOMIAL MINUS THE VALUE FROM THE   **
!     **  REFERENCE STRUCTURE                                         **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: VEC(3,NAT)
      REAL(8)   ,INTENT(OUT):: A
      REAL(8)   ,INTENT(OUT):: B(3,NAT)
      REAL(8)   ,INTENT(OUT):: C(3,NAT,3,NAT)
      INTEGER(4)            :: IAT,I
!     ******************************************************************
      CALL CONSTRAINTS_INIT(3*NAT,A,B,C) 
      DO IAT=1,NAT
        DO I=1,3
          B(I,IAT)=VEC(I,IAT)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_COGVALUE(NAT,X_,Y_,Z_,TMEMBER,RMASS,R0,VALUE_)
!     ******************************************************************
!     **  VALUE OF CENTER OF GRAVITY                                  **
!     ****************************************************************** 
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: X_,Y_,Z_
      LOGICAL(4),INTENT(IN) :: TMEMBER(NAT)
      REAL(8)   ,INTENT(IN) :: RMASS(NAT)
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)
      REAL(8)   ,INTENT(OUT):: VALUE_
      INTEGER(4)            :: IAT
      REAL(8)               :: DLENG,TMASS,SVAR
      REAL(8)               :: X,Y,Z
!     ******************************************************************
      DLENG=SQRT(X_**2+Y_**2+Z_**2)
      X=X_/DLENG
      Y=Y_/DLENG
      Z=Z_/DLENG
      TMASS=0.D0
      SVAR=0.D0
      DO IAT=1,NAT
        IF(TMEMBER(IAT)) THEN
          TMASS=TMASS+RMASS(IAT)
          SVAR=SVAR+RMASS(IAT)*X*R0(1,IAT) &
     &             +RMASS(IAT)*Y*R0(2,IAT) &
     &             +RMASS(IAT)*Z*R0(3,IAT)
        END IF
      ENDDO
      VALUE_=SVAR/TMASS
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_SETCOG(NAT,X_,Y_,Z_,TMEMBER,RMASS,A,B,C)
!     ******************************************************************
!     **  CONSTRAINT: FIX CENTER OF GRAVITY ALONG ONE DIRECTION       **
!     **                                                              **
!     **  SETS THE TRANSLATIONAL MOMENTUM OF A FRAGMENT               **
!     **  SPECIFIED BY IND ALONG (X,Y,X) TO ZERO                      **
!     **                                                              **
!     ****************************************************************** 
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: X_,Y_,Z_
      LOGICAL(4),INTENT(IN) :: TMEMBER(NAT)
      REAL(8)   ,INTENT(IN) :: RMASS(NAT)
      REAL(8)   ,INTENT(OUT):: A
      REAL(8)   ,INTENT(OUT):: B(3,NAT)
      REAL(8)   ,INTENT(OUT):: C(3,NAT,3,NAT)
      INTEGER(4)            :: IAT
      REAL(8)               :: DLENG,TMASS
      REAL(8)               :: X,Y,Z
!     ******************************************************************
      DLENG=SQRT(X_**2+Y_**2+Z_**2)
      X=X_/DLENG
      Y=Y_/DLENG
      Z=Z_/DLENG
      TMASS=0.D0
      DO IAT=1,NAT
        IF(TMEMBER(IAT)) TMASS=TMASS+RMASS(IAT)
      ENDDO
      CALL CONSTRAINTS_INIT(3*NAT,A,B,C) 
      DO IAT=1,NAT
        IF(TMEMBER(IAT)) THEN
          B(1,IAT)=RMASS(IAT)/TMASS*X
          B(2,IAT)=RMASS(IAT)/TMASS*Y
          B(3,IAT)=RMASS(IAT)/TMASS*Z
        END IF
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_ORIENTATIONVECTOR(NAT,TMEMBER,PHI,R0,RMASS,VEC)
!     ******************************************************************
!     **                                                              **
!     **  TRANSFORMS THE CONSTRAINT FOR FIXING THE ORIENTATION        **
!     **    PHI*[SUM_I:M_I*XREF_I CROSS X_I]                          **
!     **    WITH  X_I=R_I-[SUM_J:M_JR_J]/[SUM_J:M_J]                  **
!     **                                                              **
!     ** INTO A VECTOR FOR A LINEAR CONSTRAINT                        **
!     **    [SUM_I:VEC(I)R(I)]=0                                      **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: PHI(3)
      LOGICAL(4),INTENT(IN) :: TMEMBER(NAT)
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)
      REAL(8)   ,INTENT(IN) :: RMASS(NAT)
      REAL(8)   ,INTENT(OUT):: VEC(3,NAT)
      INTEGER(4)            :: IAT
      REAL(8)               :: PHILEN,PHIX,PHIY,PHIZ
      REAL(8)               :: TOTM
      REAL(8)               :: COG(3)
      REAL(8)               :: DX,DY,DZ
!     ******************************************************************
!
!     =================================================================
!     == NORMALIZE ANGLE VECTOR                                      ==
!     =================================================================
      PHILEN=SQRT(PHI(1)**2+PHI(2)**2+PHI(3)**2)
      PHIX=PHI(1)/PHILEN
      PHIY=PHI(2)/PHILEN
      PHIZ=PHI(3)/PHILEN
!     
!     =================================================================
!     == CALCULATE CENTER OF GRAVITY =================================
!     =================================================================
      TOTM=0.D0
      COG(1)=0.D0
      COG(2)=0.D0
      COG(3)=0.D0
      DO IAT=1,NAT
        IF(TMEMBER(IAT)) THEN
          TOTM=TOTM+RMASS(IAT)
          COG(1)=COG(1)+RMASS(IAT)*R0(1,IAT)
          COG(2)=COG(2)+RMASS(IAT)*R0(2,IAT)
          COG(3)=COG(3)+RMASS(IAT)*R0(3,IAT)
        END IF
      ENDDO
      COG(1)=COG(1)/TOTM
      COG(2)=COG(2)/TOTM
      COG(3)=COG(3)/TOTM
!     
!     =================================================================
!     == CALCULATE CENTER OF GRAVITY =================================
!     =================================================================
      VEC(:,:)=0.D0
      DO IAT=1,NAT
        IF(TMEMBER(IAT)) THEN
          DX=R0(1,IAT)-COG(1)
          DY=R0(2,IAT)-COG(2)
          DZ=R0(3,IAT)-COG(3)
          VEC(1,IAT)=RMASS(IAT)*(PHIY*DZ-PHIZ*DY)
          VEC(2,IAT)=RMASS(IAT)*(PHIZ*DX-PHIX*DZ)
          VEC(3,IAT)=RMASS(IAT)*(PHIX*DY-PHIY*DX)
        END IF
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CONSTRAINTS_ROT$VALUE &
     &           (NAT,X_,Y_,Z_,TMEMBER,R0,RM,DELTA,RMASS,VALUE_)
!     **************************************************************************
!     **  ENFORCE ANGULAR MOMENTUM TO BE ZERO                                 **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: X_,Y_,Z_
      LOGICAL(4),INTENT(IN) :: TMEMBER(NAT)
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)
      REAL(8)   ,INTENT(IN) :: RM(3,NAT)
      REAL(8)   ,INTENT(IN) :: DELTA
      REAL(8)   ,INTENT(IN) :: RMASS(NAT)
      REAL(8)   ,INTENT(OUT):: VALUE_
!     **************************************************************************
!     == CHANGED BY P. BLOECHL JULY 13, 2007. THIS ROUTINE SHALL ENFORCE THE 
!     == ANGULAR MOMENTUM TO BE ZERO AND NOT TO A CONSTANT VALUE. THEREFORE 
!     == THE VALUE SHALL BE EQUAL TO ZERO AND NOT EQUAL TO THE ANGULAR MOMENTUM.
      VALUE_=0.D0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CONSTRAINTS_ROT$SET(NAT,X_,Y_,Z_,TMEMBER,R0,DELTA,RMASS,A,B,C)
!     **************************************************************************
!     **  SET A AND B FOR ANGULAR MOMENTUM CONSTRAINT                         **
!     **                                                                      **
!     **  THE ANGULAR MOMENTUM IS SET TO ZERO IN THE REFERENCE FRAME OF THE   **
!     **  CENTER OF GRAVITY. THUS ROTATION AND TRANSLATIONS CAN BE ENFORCES   **
!     **  INDEPENDENTLY                                                       **
!     **                                                                      **
!     **  PMARGL:                                                             **
!     **  CAUTION!!! THIS ROUTINE CREATES A SINGULARITY IN                    **
!     **  CONSTRAINTS_APPLY IF THE DIRECTIONAL VECTOR X,Y,Z POINTS            **
!     **  IN THE CIRECTION OF R0-COG, WHICH HAPPENS QUITE FREQUENTLY,         **
!     **  ESPECIALLY AT THE ONSET OF A SIMULATION.                            **
!     **  THERE IS PROBABLY A WAY AROUND THAT IN CONSTRAINTS_APPLY            **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: X_,Y_,Z_
      LOGICAL(4),INTENT(IN) :: TMEMBER(NAT)
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)
      REAL(8)   ,INTENT(IN) :: DELTA
      REAL(8)   ,INTENT(IN) :: RMASS(NAT)
      REAL(8)   ,INTENT(OUT):: A
      REAL(8)   ,INTENT(OUT):: B(3,NAT)
      REAL(8)   ,INTENT(OUT):: C(3,NAT,3,NAT)
      INTEGER(4)            :: IAT
      REAL(8)               :: XLEN
      REAL(8)               :: X,Y,Z ! NORMALIZED ROTATION AXIS
      REAL(8)               :: TOTM  ! TOTAL MASS 
      REAL(8)               :: COG(3)! CENTER OF GRAVITY
      REAL(8)               :: U(3)
      REAL(8)               :: DX,DY,DZ
!     **************************************************************************
      A=0.D0
      B(:,:)=0.D0
      C(:,:,:,:)=0.D0
!
!     == DETERMINE NORMALIZED ROTATION AXIS ====================================
      XLEN=SQRT(X_**2+Y_**2+Z_**2)
      X=X_/XLEN
      Y=Y_/XLEN
      Z=Z_/XLEN
!
!     == CALCULATE CENTER OF GRAVITY COG AND TOTAL MASS TOTM ===================
      TOTM=0.D0
      COG(:)=0.D0
      DO IAT=1,NAT
        IF(TMEMBER(IAT)) THEN
          TOTM=TOTM+RMASS(IAT)
          COG(:)=COG(:)+RMASS(IAT)*R0(:,IAT)
        END IF
      ENDDO
      COG(:)=COG(:)/TOTM
!     
!     ===  M*(AXIS X (R-COG))
      CALL CONSTRAINTS_INIT(3*NAT,A,B,C)
      U(:)=0.D0
      DO IAT=1,NAT
        IF(TMEMBER(IAT)) THEN
          DX=R0(1,IAT)-COG(1)
          DY=R0(2,IAT)-COG(2)
          DZ=R0(3,IAT)-COG(3)
          B(1,IAT)=RMASS(IAT)*(Y*DZ-Z*DY)
          B(2,IAT)=RMASS(IAT)*(Z*DX-X*DZ)
          B(3,IAT)=RMASS(IAT)*(X*DY-Y*DX)
          U(:)=U(:)+B(:,IAT)
        END IF
      ENDDO
!
!     === SUBTRACT TRANSLATIONAL MOTION =========================================
      U(:)=U(:)/TOTM
      DO IAT=1,NAT
        IF(TMEMBER(IAT)) THEN
          B(:,IAT)=(B(:,IAT)-U(:)*RMASS(IAT))/DELTA
        END IF
      ENDDO
        
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_FIXATOMVALUE(NAT,IAT,R,VALUE_)
!     ******************************************************************
!     **  CONSTRAINTS_FIXATOMVALUE                                    **
!     **  CALCULATE VALUE OF CONSTRAINT FROM THE REFERENCE STRUCTURE  **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: IAT
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      REAL(8)   ,INTENT(OUT):: VALUE_(3)
      INTEGER(4)            :: I
!     ******************************************************************
      DO I=1,3
        VALUE_(I)=R(I,IAT)
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_SETFIXATOM(NAT,IAT,A,B,C)
!     ******************************************************************
!     **  CONSTRAINTS_SETFIXATOM                                      **
!     **  APPROXIMATES CONSTRAINT FUNCTION (WITHOUT VALUE)            **
!     **  BY AN QUADRATIC POLYNOM IN A+B*R+0.5*R*C*R                  **
!     **  THE CONSTRAINT IS THE POLYNOMIAL MINUS THE VALUE FROM THE   **
!     **  REFERENCE STRUCTURE                                         **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: IAT
      REAL(8)   ,INTENT(OUT):: A(3)
      REAL(8)   ,INTENT(OUT):: B(3,NAT,3)
      REAL(8)   ,INTENT(OUT):: C(3,NAT,3,NAT,3)
      INTEGER(4)            :: IC
!     ******************************************************************
      DO IC=1,3
        CALL CONSTRAINTS_INIT(3*NAT,A(IC),B(1,1,IC),C(1,1,1,1,IC))
        B(IC,IAT,IC)=1.D0
      ENDDO
      RETURN
      END
!
! PMARGL START
!
!     ...................................................................
      SUBROUTINE CONSTRAINTS_COGSEPVALUE(NAT,GROUP1,GROUP2,R0,VALUE_)
!     ******************************************************************
!     **  COG DISTANCE                                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NAT
      CHARACTER(*),INTENT(IN) :: GROUP1  ! NAME OF FIRST GROUP
      CHARACTER(*),INTENT(IN) :: GROUP2  ! NAME OF SECOND GROUP
      REAL(8)     ,INTENT(IN) :: R0(3,NAT)  ! ATOMIC POSITIONS
      REAL(8)     ,INTENT(OUT):: VALUE_     ! DISTANCE OF COG
      LOGICAL(4)              :: TMEMBER1(NAT)
      LOGICAL(4)              :: TMEMBER2(NAT)
      INTEGER(4)              :: I,IAT
      REAL(8)                 :: RMASS(NAT)
      REAL(8)                 :: SUMMASS1
      REAL(8)                 :: SUMMASS2
      REAL(8)                 :: COG1(3)
      REAL(8)                 :: COG2(3)
      REAL(8)                 :: SVAR
!     ******************************************************************
!
!     ==================================================================
!     == DETERMINE GROUPS AND ATOM-MASSES
!     ==================================================================
      CALL ATOMLIST$GETR8A('MASS',0,NAT,RMASS)
      CALL GROUPLIST$MEMBERS(GROUP1,NAT,TMEMBER1)
      CALL GROUPLIST$MEMBERS(GROUP2,NAT,TMEMBER2)
!
!     ==================================================================
!     == CALCULATE CENTER OF GRAVITY FOR GROUPS
!     ==================================================================
      SUMMASS1=0.D0
      SUMMASS2=0.D0
      DO I=1,3
        COG1(I)=0.D0
        COG2(I)=0.D0
      ENDDO
      DO IAT=1,NAT
        IF (TMEMBER1(IAT)) THEN
!          IF (TMEMBER2(IAT)) THEN
!            CALL ERROR$MSG('ATOM CAN ONLY BE MEMBER OF ONE GROUP')
!            CALL ERROR$STOP('CONSTRAINTS_COGSEPVALUE')
!          ENDIF
          SUMMASS1=SUMMASS1+RMASS(IAT)
          DO I=1,3
            COG1(I)=COG1(I)+R0(I,IAT)*RMASS(IAT)
          ENDDO
        ELSE IF (TMEMBER2(IAT)) THEN
          SUMMASS2=SUMMASS2+RMASS(IAT)
          DO I=1,3
            COG2(I)=COG2(I)+R0(I,IAT)*RMASS(IAT)
          ENDDO
        ENDIF
      ENDDO
      DO I=1,3
        COG1(I)=COG1(I)/SUMMASS1
        COG2(I)=COG2(I)/SUMMASS2
      ENDDO
!
!     ==================================================================
!     == CALCULATE DISTANCE BETWEEN CENTERS OF GRAVITY
!     ==================================================================
      SVAR=0.D0
      DO I=1,3
        SVAR=SVAR+(COG1(I)-COG2(I))**2.D0
      ENDDO
      VALUE_=SQRT(SVAR)
      PRINT*,'COG GROUP1: ',TMEMBER1
      PRINT*,'COG GROUP2: ',TMEMBER2
      PRINT*,'COG  VALUE: ',VALUE_
!
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_SETCOGSEP(NAT,GROUP1,GROUP2,A,B,C,R0)
!     ******************************************************************
!     **  SET UP CONSTRAINT IN COORDINATE SPACE                       **
!     **  CONSTRAINT IS THE DISTANCE BETWEEN THE CENTER OF GRAVITIES  **
!     **  OF TWO NAMED GROUPS                                         **
!     **                                                              **
!     **  THE CONSTRAINT IS EXPRESSED AS A+B*R+0.5*R*C*R=0            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NAT
      CHARACTER(*),INTENT(IN) :: GROUP1  ! NAME OF FIRST GROUP
      CHARACTER(*),INTENT(IN) :: GROUP2  ! NAME OF SECOND GROUP
      REAL(8)     ,INTENT(OUT):: A        
      REAL(8)     ,INTENT(OUT):: B(3,NAT)
      REAL(8)     ,INTENT(OUT):: C(3,NAT,3,NAT)
      REAL(8)     ,INTENT(IN) :: R0(3,NAT)
      REAL(8)                 :: RMASS(NAT)
      LOGICAL(4)              :: TMEMBER1(NAT)
      LOGICAL(4)              :: TMEMBER2(NAT)
      REAL(8)                 :: SUMMASS1,SUMMASS2
      REAL(8)                 :: COG1(3)   ! COG OF FIRST GROUP
      REAL(8)                 :: COG2(3)   ! COG OF SECOND GROUP
      REAL(8)                 :: X(3)      ! COG1-COG2
      REAL(8)                 :: DADX(3)   
      REAL(8)                 :: D2ADX2(3,3)   
      INTEGER(4)              :: I,J,IAT,IAT1,IAT2
!     ******************************************************************
      CALL ATOMLIST$GETR8A('MASS',0,NAT,RMASS)
      CALL GROUPLIST$MEMBERS(GROUP1,NAT,TMEMBER1)
      CALL GROUPLIST$MEMBERS(GROUP2,NAT,TMEMBER2)
      CALL CONSTRAINTS_INIT(3*NAT,A,B,C)
!
!     ==================================================================
!     ==  CALCULATE CENTERS OF GRAVITY                                ==
!     ==================================================================
      SUMMASS1=0.D0
      SUMMASS2=0.D0
      COG1(:)=0.D0
      COG2(:)=0.D0
      DO IAT=1,NAT
        IF(TMEMBER1(IAT)) THEN
          SUMMASS1=SUMMASS1+RMASS(IAT)
          COG1(:)=COG1(:)+RMASS(IAT)*R0(:,IAT)
        ENDIF
        IF(TMEMBER2(IAT)) THEN
          SUMMASS2=SUMMASS2+RMASS(IAT)
          COG2(:)=COG2(:)+RMASS(IAT)*R0(:,IAT)
        ENDIF
      ENDDO
      COG1(:)=COG1(:)/SUMMASS1
      COG2(:)=COG2(:)/SUMMASS2
!
!     =============================================================
!     ==  SEPARATION BETWEEN THE CENTERS OF GRAVITY              ==
!     =============================================================
      X(:)=COG1(:)-COG2(:)
!
!     =============================================================
!     ==  CALCULATE A(X) AND ITS DERIVATIVES                     ==
!     =============================================================
      A=SQRT(DOT_PRODUCT(X(:),X(:)))
      DADX(:)=X(:)/A
      D2ADX2(:,:)=0.D0
      DO I=1,3
        DO J=1,3
          D2ADX2(I,J)=-X(I)*X(J)/A**3
        ENDDO
        D2ADX2(I,I)=D2ADX2(I,I)+1.D0/A
      ENDDO  
!
!     =============================================================
!     ==  TRANSFORM FROM X TO R                                  ==
!     =============================================================
      B(:,:)=0.D0
      C(:,:,:,:)=0.D0
      DO IAT1=1,NAT
        IF(TMEMBER1(IAT1)) THEN
          B(:,IAT1)=B(:,IAT1)+DADX(:)*RMASS(IAT1)/SUMMASS1
        END IF
        IF(TMEMBER2(IAT1)) THEN
          B(:,IAT1)=B(:,IAT1)-DADX(:)*RMASS(IAT1)/SUMMASS2
        END IF
        DO IAT2=1,NAT
!         == BOTH ATOMS IN FIRST GROUP
          IF(TMEMBER1(IAT1).AND.TMEMBER1(IAT2)) THEN
            C(:,IAT1,:,IAT2)=C(:,IAT1,:,IAT2) &
       &               +D2ADX2(:,:)*RMASS(IAT1)*RMASS(IAT2)/SUMMASS1**2
          END IF
!         == BOTH ATOMS IN SECOND GROUP
          IF(TMEMBER2(IAT1).AND.TMEMBER2(IAT2)) THEN
            C(:,IAT1,:,IAT2)=C(:,IAT1,:,IAT2) &
       &               +D2ADX2(:,:)*RMASS(IAT1)*RMASS(IAT2)/SUMMASS2**2
          END IF
!         == IAT1 IN FIRST AND IAT2 IN SECOND GROUP
          IF(TMEMBER1(IAT1).AND.TMEMBER2(IAT2)) THEN
            C(:,IAT1,:,IAT2)=C(:,IAT1,:,IAT2) &
       &         -D2ADX2(:,:)*RMASS(IAT1)*RMASS(IAT2)/(SUMMASS1*SUMMASS2)
          END IF
!         == IAT1 IN SECOND AND IAT2 IN FIRST GROUP
          IF(TMEMBER2(IAT1).AND.TMEMBER1(IAT2)) THEN
            C(:,IAT1,:,IAT2)=C(:,IAT1,:,IAT2) &
       &       -D2ADX2(:,:)*RMASS(IAT1)*RMASS(IAT2)/(SUMMASS1*SUMMASS2)
          END IF 
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  TRANSFORM FROM DERIVATIVES TO POLYNOMIAL COEFFICIENTS       ==
!     ==================================================================
      CALL CONSTRAINTS_DERIVTOPOLYNOM(3*NAT,R0,A,B,C)
!
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_BONDANGLEVALUE(NAT,IAT1,IAT2,IAT3,R0,VALUE_)
!     ******************************************************************
!     **
!     **  CONSTRAINS THE ANGLE ENCLOSED BETWEEN ATOMS IAT1--IAT2--IAT3
!     **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: IAT1
      INTEGER(4),INTENT(IN) :: IAT2
      INTEGER(4),INTENT(IN) :: IAT3
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)
      REAL(8)   ,INTENT(OUT):: VALUE_
      REAL(8)               :: U(3),V(3),UV,UU,VV
      INTEGER(4)            :: I
!     ******************************************************************
      DO I=1,3
        U(I)=R0(I,IAT1)-R0(I,IAT2)
        V(I)=R0(I,IAT3)-R0(I,IAT2)
      ENDDO
      UV=DOT_PRODUCT(U,V)
      UU=DOT_PRODUCT(U,U)
      VV=DOT_PRODUCT(V,V)
      VALUE_=DACOS(UV/SQRT(UU*VV))
      PRINT*,'BONDANGLE IN REFERENCE',VALUE_
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_SETBONDANGLE(NAT,IAT1,IAT2,IAT3,A,B,C,R0)
!     ******************************************************************
!     **
!     **  CONSTRAINS THE ANGLE ENCLOSED BETWEEN ATOMS IAT1--IAT2--IAT3
!     **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: IAT1
      INTEGER(4),INTENT(IN) :: IAT2
      INTEGER(4),INTENT(IN) :: IAT3
      REAL(8)   ,INTENT(OUT):: A
      REAL(8)   ,INTENT(OUT):: B(3,NAT)
      REAL(8)   ,INTENT(OUT):: C(3,NAT,3,NAT)
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)
      REAL(8)               :: X(9),DI(9),DIDJ(9,9),B1(9)
      REAL(8)               :: PREFAC
      REAL(8)               :: PHI
      INTEGER(4)            :: I,J,ICHK
!     ******************************************************************
!
!     ==================================================================
!     ==  SET UP CONSTRAINT IN COORDINATE SPACE                       ==
!     ==================================================================
!
!     ==================================================================
!     == CALCULATE VALUE AND FIRST TWO DERIVATIVES OF                 ==
!     == PHI=<X1-X2|X3-X2>/SQRT(<X1-X2|X1-X2><X3-X2|X3-X2>)           ==
!     ==================================================================
! --- TRANSFER R0 TO X ARRAY TO MAKE HANDLING AND PASSING EASIER
      DO I=1,3
       X(I)  =R0(I,IAT1)
       X(I+3)=R0(I,IAT2)
       X(I+6)=R0(I,IAT3)
      ENDDO
!
! --- CALCULATE PHI,D PHI/D X_I, D^2 PHI/D X_I D X_J
      CALL CONSTRAINTS_BONDANGLEDERIVS(X,PHI,DI,DIDJ,ICHK)
!PRINT*,'PHI FROM BONDANGLEDERIVS',PHI 
!PRINT*,'DI ',DI 
!DO I=1,9
!  DO J=1,9
!    IF(ABS(DIDJ(I,J)).GT.1.D+5) THEN
!      PRINT*,'WARNING! ',I,J,DIDJ(I,J)
!      STOP
!    END IF
!  ENDDO
!ENDDO
!
!
! --- INITIALIZE CONSTRAINTS EQUATION ARRAYS
      CALL CONSTRAINTS_INIT(3*NAT,A,B,C)
!
!     ==================================================================
!     ==  TRANSFORM DERIVATIVES TO POLYNOMIAL COEFFICIENTS            ==
!     ==================================================================
      PRINT*,'CURRENT BOND ANGLE FROM DERIVS:',PHI
      A=PHI
! --- + -\SUM_1^9 (D PHI/D X_I) *X_I0
      DO I=1,9
       A=A-DI(I)*X(I)
      ENDDO
! --- + 1/PREFAC \SUM_{I<=J} (D^2 PHI / DX_I DX_J)* X_I,0 X_J,0      
      DO J=1,9
        DO I=1,J
          IF (I .EQ. J) THEN
            PREFAC=0.5D0
          ELSE
           PREFAC=1.D0
          ENDIF
          A=A+PREFAC*DIDJ(I,J)*X(I)*X(J)
        ENDDO
      ENDDO
! -------------------------------------------------------------
! --- CALCULATE B
! -------------------------------------------------------------
! --- -\SUM_J=1^9 D^2 PHI/DI,J * X(J)
      DO I=1,9
        B1(I)=0.D0
        DO J=1,9
          B1(I)=B1(I)-DIDJ(I,J)*X(J)
        ENDDO
      ENDDO
! --- B(I)=B(I)+D PHI/DX(I)
      DO I=1,9
        B1(I)=B1(I)+DI(I)
      ENDDO
! -------------------------------------------------------------
! --- CALCULATE C
! -------------------------------------------------------------
! --- C IS SIMPLY THE MATRIX OF SECOND DERIVATIVES
      DO I=1,3
        B(I,IAT1)=B1(I)
        B(I,IAT2)=B1(I+3)
        B(I,IAT3)=B1(I+6)
        DO J=1,3
          C(I,IAT1,J,IAT1)=DIDJ(I,J)
          C(I,IAT2,J,IAT2)=DIDJ(I+3,J+3)
          C(I,IAT3,J,IAT3)=DIDJ(I+6,J+6)
          C(I,IAT1,J,IAT2)=DIDJ(I,J+3)
          C(I,IAT1,J,IAT3)=DIDJ(I,J+6)
          C(I,IAT2,J,IAT1)=DIDJ(I+3,J)
          C(I,IAT3,J,IAT1)=DIDJ(I+6,J)
          C(I,IAT2,J,IAT3)=DIDJ(I+3,J+6)
          C(I,IAT3,J,IAT2)=DIDJ(I+6,J+3)
        ENDDO
      ENDDO
!
!     CALL FILEHANDLER$UNIT('PROT',NFILO)
!     WRITE(NFILO,FMT='("BONDANGLE CONSTRAINT",F10.5)')SVAR
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_BONDANGLEDERIVS(X,PHI,DI,DIDJ,ICHK)
!     -------------------------------------------------------------
!     --- CALCULATE ANGLE AND ANALYTICAL DERIVATIVES             --
!     --- 9DIM VECTOR X IS HANDED OVER, PHI,DI,DIDJ ARE RETURNED -- 
!     -------------------------------------------------------------
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: X(9)
      REAL(8)   ,INTENT(OUT):: PHI
      REAL(8)   ,INTENT(OUT):: DI(9)
      REAL(8)   ,INTENT(OUT):: DIDJ(9,9)
      INTEGER(4),INTENT(OUT):: ICHK
      REAL(8)               :: A    ! U*U
      REAL(8)               :: B    ! V*V
      REAL(8)               :: C    ! U*V
      REAL(8)               :: D    ! 1/SQRT(A*B)
      REAL(8)               :: DADX(9)
      REAL(8)               :: DBDX(9)
      REAL(8)               :: DCDX(9)
      REAL(8)               :: DDDX(9)
      REAL(8)               :: DJDX(9)
      REAL(8)               :: DVDX(9)
      REAL(8)               :: DGDX(9)    ! D(A*B)/DX
      REAL(8)               :: D2ADX2(9,9)
      REAL(8)               :: D2BDX2(9,9)
      REAL(8)               :: D2CDX2(9,9)
      REAL(8)               :: D2DDX2(9,9)
      REAL(8)               :: AB,AB2N,V,PREFAC
      REAL(8)               :: COSPHI,SINPHI
      REAL(8)               :: DIVDX,AI
      INTEGER(4)            :: I,J
!     ******************************************************************
      ICHK=0
!
!     ----------------------------------------------------------------
!     --- R1=X(1:3)   R2=X(4:6)  R3=X(7:9)
!     --- U=R1-R2 V=R3-R2
!     --- A=U*U  B=V*V  C=U*V
!     ----------------------------------------------------------------
      A=(X(1)-X(4))**2 + (X(2)-X(5))**2 + (X(3)-X(6))**2
      B=(X(7)-X(4))**2 + (X(8)-X(5))**2 + (X(9)-X(6))**2
      C=(X(1)-X(4))*(X(7)-X(4)) &
     & +(X(2)-X(5))*(X(8)-X(5)) &
     & +(X(3)-X(6))*(X(9)-X(6))
!     ---  CALCULATE DADX
      DADX(1)= 2.D0*(X(1)-X(4))
      DADX(2)= 2.D0*(X(2)-X(5))
      DADX(3)= 2.D0*(X(3)-X(6))
      DADX(4)=-2.D0*(X(1)-X(4))
      DADX(5)=-2.D0*(X(2)-X(5))
      DADX(6)=-2.D0*(X(3)-X(6))
      DADX(7)= 0.D0
      DADX(8)= 0.D0
      DADX(9)= 0.D0
!     ---  CALCULATE DBDX
      DBDX(1)= 0.D0
      DBDX(2)= 0.D0
      DBDX(3)= 0.D0
      DBDX(4)=-2.D0*(X(7)-X(4))
      DBDX(5)=-2.D0*(X(8)-X(5))
      DBDX(6)=-2.D0*(X(9)-X(6))
      DBDX(7)= 2.D0*(X(7)-X(4))
      DBDX(8)= 2.D0*(X(8)-X(5))
      DBDX(9)= 2.D0*(X(9)-X(6))
! --- CALCULATE DCDX
      DCDX(1)= X(7)-X(4)
      DCDX(2)= X(8)-X(5)
      DCDX(3)= X(9)-X(6)
      DCDX(4)=-X(7)+2.D0*X(4)-X(1)
      DCDX(5)=-X(8)+2.D0*X(5)-X(2)
      DCDX(6)=-X(9)+2.D0*X(6)-X(3)
      DCDX(7)= X(1)-X(4)
      DCDX(8)= X(2)-X(5)
      DCDX(9)= X(3)-X(6)
! ----------------------------------------------------------------
! --- CALCULATE SECOND DERIVATIVES
! ----------------------------------------------------------------
! --- FILL UP WITH ZERO, SINCE MATRICES ARE SPARSE
      D2ADX2(:,:)=0.D0
      D2BDX2(:,:)=0.D0
      D2CDX2(:,:)=0.D0
! --- CALCULATE SECOND DERIVATIVES OF A
      D2ADX2(1,4)=-2.D0
      D2ADX2(1,1)=+2.D0
      D2ADX2(2,5)=-2.D0
      D2ADX2(2,2)=+2.D0
      D2ADX2(3,6)=-2.D0
      D2ADX2(3,3)=+2.D0
      D2ADX2(4,4)=+2.D0
      D2ADX2(5,5)=+2.D0
      D2ADX2(6,6)=+2.D0
! --- CALCULATE SECOND DERIVATIVES OF B
      D2BDX2(4,4)=+2.D0
      D2BDX2(4,7)=-2.D0
      D2BDX2(5,5)=+2.D0
      D2BDX2(5,8)=-2.D0
      D2BDX2(6,6)=+2.D0
      D2BDX2(6,9)=-2.D0
      D2BDX2(7,7)=+2.D0
      D2BDX2(8,8)=+2.D0
      D2BDX2(9,9)=+2.D0
! --- CALCULATE SECOND DERIVATIVES OF C
      D2CDX2(1,4)=-1.D0
      D2CDX2(1,7)=+1.D0
      D2CDX2(2,5)=-1.D0
      D2CDX2(2,8)=+1.D0
      D2CDX2(3,6)=-1.D0
      D2CDX2(3,9)=+1.D0
      D2CDX2(4,4)=+2.D0
      D2CDX2(4,7)=-1.D0
      D2CDX2(5,5)=+2.D0
      D2CDX2(5,8)=-1.D0
      D2CDX2(6,6)=+2.D0
      D2CDX2(6,9)=-1.D0
! --- MIRROR THROUGH DIAGONAL
      DO J=1,8
       DO I=J+1,9
        D2ADX2(I,J)=D2ADX2(J,I)
        D2BDX2(I,J)=D2BDX2(J,I)
        D2CDX2(I,J)=D2CDX2(J,I)
       ENDDO
      ENDDO
!
!     ==================================================================
!     ==  CALCULATE PHI ETC  ===========================================
!     ==================================================================
      AB=A*B
      D=SQRT(AB)
      COSPHI=C/D       !=COS(PHI)
      PHI=DACOS(COSPHI)
!     --- CALCULATE DERIVATIVE OF ARCCOS
      SINPHI=SQRT(1.D0-COSPHI**2)
      V=-1.D0/SINPHI                 !=-SIN(PHI)=D(ACOS(PHI))/D(COSPHI)
! --- CALCULATE FIRST DERIVATIVES OF G=A*B
      DO I=1,9
        DGDX(I)=DADX(I)*B+A*DBDX(I)
      ENDDO
! --- CALCULATE DDDX
      AB2N=0.5D0/D
      DO I=1,9
        DDDX(I)=AB2N*DGDX(I)
      ENDDO
! --- CALCULATE FIRST DERIVATIVES OF J
      PREFAC=1.D0/AB
      DO I=1,9
        DJDX(I)=(DCDX(I)*D-C*DDDX(I))*PREFAC
      ENDDO
! --- CALCULATE FIRST DERIVATIVE
      DO I=1,9
        DI(I)=DJDX(I)*V
      ENDDO
! --- CALCULATE SECOND DERIVATIVES OF D
      PREFAC=-0.5D0*AB2N/AB
      DO I=1,9
        DO J=1,9
          D2DDX2(I,J)=PREFAC*(DGDX(J))*(DGDX(I)) &
     &      + AB2N*(D2ADX2(I,J)*B+DADX(I)*DBDX(J)+DADX(J)*DBDX(I)&
     &      + A*D2BDX2(I,J))
!        D2DDX2(I,J)=-0.25*(AB)**(-1.5D0)*DGDX(J)*DGDX(I) &
!     &   + AB2N*(D2ADX2(I,J)*B+DADX(I)*DBDX(J)+DADX(J)*DBDX(I) &
!     &   + A*D2BDX2(I,J))
        ENDDO
      ENDDO
!
! --- CALCULATE FIRST DERIVATIVES OF V
      PREFAC=COSPHI/SINPHI**3
      DO I=1,9
        DVDX(I)=-DJDX(I)*PREFAC
      ENDDO
!
! --- NOW FINALLY CALCULATE SECOND DERIVATIVES OF PHI
      PREFAC=1.D0/AB
      DO I=1,9
        AI=DJDX(I)*AB
        DO J=1,9
          DIVDX=(D2CDX2(I,J)*D+DCDX(I)*DDDX(J)-DCDX(J)*DDDX(I) &
     &         -C*D2DDX2(I,J))*V+ AI*DVDX(J)
          DIDJ(I,J)=(DIVDX-DI(I)*DGDX(J))*PREFAC
        ENDDO
      ENDDO
      RETURN
      END
!
!     .................................................................
      SUBROUTINE CONSTRAINTS_TORSIONVALUE(NAT,IAT1,IAT2,IAT3,IAT4,R0,VALUE_)
!     ******************************************************************
!     **
!     **  CONSTRAINS THE TORSION OF IAT1-IAT2---IAT3-IAT4
!     **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: IAT1
      INTEGER(4),INTENT(IN) :: IAT2
      INTEGER(4),INTENT(IN) :: IAT3
      INTEGER(4),INTENT(IN) :: IAT4
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)
      REAL(8)   ,INTENT(OUT):: VALUE_
      REAL(8)               :: T(3),U(3),R(3),W(3),P(3),Q(3) 
      REAL(8)               :: CX,AX,BX,ABYB,D,COVERD
!     ******************************************************************
      T(:)=R0(:,IAT1)-R0(:,IAT2)
      U(:)=R0(:,IAT3)-R0(:,IAT2)
      R(:)=R0(:,IAT4)-R0(:,IAT3)
      W(:)=R0(:,IAT2)-R0(:,IAT3)
      CALL VPRD(T,U,P)                     ! P=CROSS-PRODUCT(T,U)
      CALL VPRD(W,R,Q)
      CX=DOT_PRODUCT(P,Q)
      AX=DOT_PRODUCT(P,P)
      BX=DOT_PRODUCT(Q,Q)
      ABYB=AX*BX
      D=SQRT(ABYB)
      COVERD=CX/D
      VALUE_=DACOS(COVERD)
      PRINT*,'TORSION IN REFERENCE',VALUE_
!
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_SETTORSION(NAT,IAT1,IAT2,IAT3,IAT4,A,B,C,R0)
!     ******************************************************************
!     **
!     **  CONSTRAINS THE TORSION OF IAT1-IAT2---IAT3-IAT4
!     **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: IAT1
      INTEGER(4),INTENT(IN) :: IAT2
      INTEGER(4),INTENT(IN) :: IAT3
      INTEGER(4),INTENT(IN) :: IAT4
      REAL(8)   ,INTENT(IN) :: R0(3,NAT)
      REAL(8)   ,INTENT(OUT):: A
      REAL(8)   ,INTENT(OUT):: B(3,NAT)
      REAL(8)   ,INTENT(OUT):: C(3,NAT,3,NAT)
      REAL(8)               :: X(12),DI(12),DIDJ(12,12),B1(12)
      REAL(8)               :: PHI
      REAL(8)               :: PREFAC
      INTEGER(4)            :: I,J,ICHK
!     ******************************************************************
!
! --- TRANSFER R0 TO X ARRAY TO MAKE HANDLING AND PASSING EASIER
      DO I=1,3
        X(I)=R0(I,IAT1)
        X(I+3)=R0(I,IAT2)
        X(I+6)=R0(I,IAT3)
        X(I+9)=R0(I,IAT4)
      ENDDO
!
! --- CALCULATE PHI,D PHI/D X_I, D^2 PHI/D X_I D X_J
      CALL CONSTRAINTS_TORSIONDERIVS(X,PHI,DI,DIDJ,ICHK)
!
! --- INITIALIZE CONSTRAINTS EQUATION ARRAYS
      CALL CONSTRAINTS_INIT(3*NAT,A,B,C)
! -------------------------------------------------------------
! --- CALCULATE A
! -------------------------------------------------------------
      PRINT*,'CURRENT TORSION FROM DERIVS:',PHI
      A=PHI
! --- + -\SUM_1^9 (D PHI/D X_I) *X_I0
      DO I=1,12
       A=A-DI(I)*X(I)
      ENDDO
! --- + 1/PREFAC \SUM_{I<=J} (D^2 PHI / DX_I DX_J)* X_I,0 X_J,0      
      DO J=1,12
       DO I=1,J
        IF (I .EQ. J) THEN
         PREFAC=0.5D0
        ELSE
         PREFAC=1.D0
        ENDIF
        A=A+PREFAC*DIDJ(I,J)*X(I)*X(J)
       ENDDO
      ENDDO
! -------------------------------------------------------------
! --- CALCULATE B
! -------------------------------------------------------------
! --- -\SUM_J=1^9 D^2 PHI/DI,J * X(J)
      DO I=1,12
       B1(I)=0.D0
       DO J=1,12
        B1(I)=B1(I)-DIDJ(I,J)*X(J)
       ENDDO
      ENDDO
! --- B(I)=B(I)+D PHI/DX(I)
      DO I=1,12
       B1(I)=B1(I)+DI(I)
      ENDDO
! -------------------------------------------------------------
! --- CALCULATE C
! -------------------------------------------------------------
! --- C IS SIMPLY THE MATRIX OF SECOND DERIVATIVES
      DO I=1,3
       B(I,IAT1)=B1(I)
       B(I,IAT2)=B1(I+3)
       B(I,IAT3)=B1(I+6)
       B(I,IAT4)=B1(I+9)
       DO J=1,3
        C(I,IAT1,J,IAT1)=DIDJ(I,J)
        C(I,IAT2,J,IAT2)=DIDJ(I+3,J+3)
        C(I,IAT3,J,IAT3)=DIDJ(I+6,J+6)
        C(I,IAT4,J,IAT4)=DIDJ(I+9,J+9)
        C(I,IAT1,J,IAT2)=DIDJ(I,J+3)
        C(I,IAT1,J,IAT3)=DIDJ(I,J+6)
        C(I,IAT1,J,IAT4)=DIDJ(I,J+9)
        C(I,IAT2,J,IAT1)=DIDJ(I+3,J)
        C(I,IAT3,J,IAT1)=DIDJ(I+6,J)
        C(I,IAT4,J,IAT1)=DIDJ(I+9,J)
        C(I,IAT2,J,IAT3)=DIDJ(I+3,J+6)
        C(I,IAT2,J,IAT4)=DIDJ(I+3,J+9)
        C(I,IAT3,J,IAT2)=DIDJ(I+6,J+3)
        C(I,IAT4,J,IAT2)=DIDJ(I+9,J+3)
        C(I,IAT3,J,IAT4)=DIDJ(I+6,J+9)
        C(I,IAT4,J,IAT3)=DIDJ(I+9,J+6)
       ENDDO
      ENDDO
!
!     CALL FILEHANDLER$UNIT('PROT',NFILO)
!     WRITE(NFILO,FMT='("BONDANGLE CONSTRAINT",F10.5)')SVAR
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_TORSIONDERIVS(X,TORS,DI,DIDJ,ICHK)
!     -------------------------------------------------------------
!     --- CALCULATE TORSION AND ANALYTICAL DERIVATIVES           --
!     --- 12DIM VECTOR X IS HANDED OVER, PHI,DI,DIDJ ARE RETURNED-- 
!     -------------------------------------------------------------
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: X(12)
      REAL(8)   ,INTENT(OUT):: TORS
      REAL(8)   ,INTENT(OUT):: DI(12)
      REAL(8)   ,INTENT(OUT):: DIDJ(12,12)
      INTEGER(4),INTENT(OUT):: ICHK
      REAL(8) :: P(3),Q(3),T(3),U(3),R(3),W(3)
      REAL(8) :: DADX(12),DBDX(12),DCDX(12)
      REAL(8) :: DDDX(12),DJDX(12),DVDX(12),DGDX(12),DPDX(3,12),DQDX(3,12)
      REAL(8) :: D2QDX2(3,12,12),D2PDX2(3,12,12)
      REAL(8) :: D2BDX2(12,12),D2ADX2(12,12),D2CDX2(12,12),D2DDX2(12,12)
      REAL(8) :: D,COVERD,ABYB2N,ABYB,C,A,B,AI,DIVDX,PREFAC,V
      INTEGER(4) :: I,J,K
!     ******************************************************************
      ICHK=0
      T(:)=X(1:3)  -X(4:6)
      U(:)=X(7:9)  -X(4:6)
      R(:)=X(10:12)-X(7:9)
      W(:)=X(4:6)  -X(7:9)
      CALL VPRD(T,U,P)
      CALL VPRD(W,R,Q)
      C=DOT_PRODUCT(P,Q)
      A=DOT_PRODUCT(P,P)
      B=DOT_PRODUCT(Q,Q)
!
      ABYB=A*B
      D=ABYB**0.5D0
      COVERD=C/D
      ABYB2N=0.5D0/D
!
      TORS=DACOS(COVERD)
!
! ----------------------------------------------------------------
! --- CALCULATE FIRST DERIVATIVES
! ----------------------------------------------------------------
! --- CALCULATE DERIVATIVE OF ARCCOS
      IF (ABS(ABS(COVERD)-1.D0) < 1.D-9) THEN
        PRINT*,'WARNING: DERIVATIVE OF TORSION AT ZERO AND PI IS UNDEFINED!'
        V=0.D0
      ELSE
        V=-1.D0/SQRT(1.D0-(COVERD)**2.D0)
      ENDIF
! --- CALCULATE DPDX
! --- ZERO OUT SPARSE MATRIX
      DO J=1,3
       DO I=1,12
        DPDX(J,I)=0.D0
        DQDX(J,I)=0.D0
       ENDDO
      ENDDO
      DPDX(1,2)=X(9)-X(6)
      DPDX(1,3)=X(5)-X(8)
      DPDX(1,5)=X(3)-X(9)
      DPDX(1,6)=X(8)-X(2)
      DPDX(1,8)=X(6)-X(3)
      DPDX(1,9)=X(2)-X(5)
      DPDX(2,1)=X(6)-X(9)
      DPDX(2,3)=X(7)-X(4)
      DPDX(2,4)=X(9)-X(3)
      DPDX(2,6)=X(1)-X(7)
      DPDX(2,7)=X(3)-X(6)
      DPDX(2,9)=X(4)-X(1)
      DPDX(3,1)=X(8)-X(5)
      DPDX(3,2)=X(4)-X(7)
      DPDX(3,4)=X(2)-X(8)
      DPDX(3,5)=X(7)-X(1)
      DPDX(3,7)=X(5)-X(2)
      DPDX(3,8)=X(1)-X(4)
! --- CALCULATE DQDX
      DQDX(1,5)=X(12)-X(9)
      DQDX(1,6)=X(8)-X(11)
      DQDX(1,8)=X(6)-X(12)
      DQDX(1,9)=X(11)-X(5)
      DQDX(1,11)=X(9)-X(6)
      DQDX(1,12)=X(5)-X(8)
      DQDX(2,4)=X(9)-X(12)
      DQDX(2,6)=X(10)-X(7)
      DQDX(2,7)=X(12)-X(6)
      DQDX(2,9)=X(4)-X(10)
      DQDX(2,10)=X(6)-X(9)
      DQDX(2,12)=X(7)-X(4)
      DQDX(3,4)=X(11)-X(8)
      DQDX(3,5)=X(7)-X(10)
      DQDX(3,7)=X(5)-X(11)
      DQDX(3,8)=X(10)-X(4)
      DQDX(3,10)=X(8)-X(5)
      DQDX(3,11)=X(4)-X(7)
! ---  CALCULATE DADX
      DO I=1,12
       DADX(I)=0.D0
       DBDX(I)=0.D0
       DCDX(I)=0.D0
       DO J=1,3
        DADX(I)=DADX(I)+2.D0*P(J)*DPDX(J,I)
        DBDX(I)=DBDX(I)+2.D0*Q(J)*DQDX(J,I)
        DCDX(I)=DCDX(I)+DPDX(J,I)*Q(J)+P(J)*DQDX(J,I)
       ENDDO
      ENDDO
! --- CALCULATE FIRST DERIVATIVES OF G
      DO I=1,12
       DGDX(I)=DADX(I)*B+A*DBDX(I)
      ENDDO
! --- CALCULATE DDDX
      DO I=1,12
       DDDX(I)=ABYB2N*DGDX(I)
      ENDDO
! --- CALCULATE FIRST DERIVATIVES OF J
      PREFAC=1.D0/ABYB
      DO I=1,12
       DJDX(I)=(DCDX(I)*D-C*DDDX(I))*PREFAC
      ENDDO
! --- CALCULATE FIRST DERIVATIVE
      DO I=1,12
       DI(I)=DJDX(I)*V
      ENDDO
! ----------------------------------------------------------------
! --- CALCULATE SECOND DERIVATIVES
! ----------------------------------------------------------------
! --- FILL UP WITH ZERO, SINCE MATRICES ARE SPARSE
      DO J=1,3
       DO I=1,12
        DO K=1,12
         D2PDX2(J,I,K)=0.D0
         D2QDX2(J,I,K)=0.D0
        ENDDO
       ENDDO
      ENDDO
! --- CALCULATE D2PDX2 
      D2PDX2(1,2,9)=1.D0
      D2PDX2(1,2,6)=-1.D0
      D2PDX2(1,3,8)=-1.D0
      D2PDX2(1,3,5)=1.D0
      D2PDX2(1,5,9)=-1.D0
      D2PDX2(1,6,8)=1.D0
!
      D2PDX2(2,1,6)=1.D0
      D2PDX2(2,1,9)=-1.D0
      D2PDX2(2,3,4)=-1.D0
      D2PDX2(2,3,7)=1.0
      D2PDX2(2,4,9)=1.D0
      D2PDX2(2,6,7)=-1.D0
!
      D2PDX2(3,1,5)=-1.D0
      D2PDX2(3,1,8)=1.D0
      D2PDX2(3,2,4)=1.D0
      D2PDX2(3,2,7)=-1.D0
      D2PDX2(3,4,8)=-1.D0
      D2PDX2(3,5,7)=1.D0
!
! --- CALCULATE D2QDX2
      D2QDX2(1,5,9)=-1.D0
      D2QDX2(1,5,12)=1.D0
      D2QDX2(1,6,8)=1.D0
      D2QDX2(1,6,11)=-1.D0
      D2QDX2(1,8,12)=-1.D0
      D2QDX2(1,9,11)=1.D0
!
      D2QDX2(2,4,12)=-1.D0
      D2QDX2(2,4,9)=1.D0
      D2QDX2(2,6,10)=1.D0
      D2QDX2(2,6,7)=-1.D0
      D2QDX2(2,7,12)=1.D0
      D2QDX2(2,9,10)=-1.D0
!
      D2QDX2(3,4,11)=1.D0
      D2QDX2(3,4,8)=-1.D0
      D2QDX2(3,5,10)=-1.D0
      D2QDX2(3,5,7)=1.D0
      D2QDX2(3,7,11)=-1.D0
      D2QDX2(3,8,10)=1.D0
!
      DO J=1,3
       DO I=1,11
        DO K=I+1,12
         D2QDX2(J,K,I)=D2QDX2(J,I,K)
         D2PDX2(J,K,I)=D2PDX2(J,I,K)
        ENDDO
       ENDDO
      ENDDO
!
      DO I=1,12
       DO J=1,12
        D2ADX2(I,J)=0.D0
        D2BDX2(I,J)=0.D0
        D2CDX2(I,J)=0.D0
       ENDDO
      ENDDO
! --- CALCULATE SECOND DERIVATIVES OF A,B,C
      DO J=1,3
       DO I=1,12
        DO K=1,12
         D2ADX2(I,K)=D2ADX2(I,K) &
     &    + 2.D0*(DPDX(J,K)*DPDX(J,I)+P(J)*D2PDX2(J,I,K))
         D2BDX2(I,K)=D2BDX2(I,K) &
     &    + 2.D0*(DQDX(J,K)*DQDX(J,I)+Q(J)*D2QDX2(J,I,K))
         D2CDX2(I,K)=D2CDX2(I,K) &
     &    + (D2PDX2(J,I,K)*Q(J)+DPDX(J,I)*DQDX(J,K) &
     &    +  DPDX(J,K)*DQDX(J,I)+P(J)*D2QDX2(J,I,K))
        ENDDO
       ENDDO
      ENDDO
! --- MIRROR THROUGH DIAGONAL
      DO J=1,11
       DO I=J+1,12
        D2ADX2(I,J)=D2ADX2(J,I)
        D2BDX2(I,J)=D2BDX2(J,I)
        D2CDX2(I,J)=D2CDX2(J,I)
       ENDDO
      ENDDO
! --- CALCULATE SECOND DERIVATIVES OF D
      PREFAC=-0.5D0*ABYB2N/ABYB
      DO I=1,12
       DO J=1,12
        D2DDX2(I,J)=PREFAC*(DGDX(J))*(DGDX(I)) &
     &    + ABYB2N*(D2ADX2(I,J)*B+DADX(I)*DBDX(J)+DADX(J)*DBDX(I) &
     &    + A*D2BDX2(I,J))
!        D2DDX2(I,J)=-0.25*(ABYB)**(-1.5D0)*DGDX(J)*DGDX(I) &
!     &  + ABYB2N*(D2ADX2(I,J)*B+DADX(I)*DBDX(J)+DADX(J)*DBDX(I) &
!     &  + A*D2BDX2(I,J))
       ENDDO
      ENDDO
! --- CALCULATE FIRST DERIVATIVES OF V
      PREFAC=COVERD*(1.D0-COVERD**2.D0)**(-1.5D0)
      DO I=1,12
       DVDX(I)=-DJDX(I)*PREFAC
      ENDDO
! --- NOW FINALLY CALCULATE SECOND DERIVATIVES OF PHI
      PREFAC=1.D0/ABYB
      DO I=1,12
       AI=DJDX(I)*ABYB
       DO J=1,12
        DIVDX=(D2CDX2(I,J)*D+DCDX(I)*DDDX(J)-DCDX(J)*DDDX(I) &
     & - C*D2DDX2(I,J))*V+ AI*DVDX(J)
        DIDJ(I,J)=(DIVDX-DI(I)*DGDX(J))*PREFAC
       ENDDO
      ENDDO
!
!
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE VPRD(X,Y,Z)
!     ******************************************************************
!     **  CROSS PRODUCT (VECTOR PRODUCT OF TWO VECTORS)               **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: X(3)
      REAL(8),INTENT(IN) :: Y(3)
      REAL(8),INTENT(OUT):: Z(3)
!     ******************************************************************
      Z(1)=X(2)*Y(3)-X(3)*Y(2)
      Z(2)=X(3)*Y(1)-X(1)*Y(3)
      Z(3)=X(1)*Y(2)-X(2)*Y(1)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_INIT(NX,A,B,C)
!     ******************************************************************
!     **                                                              **
!     **  SET  CONSTRAINT TO ZERO                                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NX
      REAL(8)   ,INTENT(OUT):: A
      REAL(8)   ,INTENT(OUT):: B(NX)
      REAL(8)   ,INTENT(OUT):: C(NX*NX)
      INTEGER(4)            :: I
!     ******************************************************************
      A=0.D0
      DO I=1,NX
        B(I)=0.D0
      ENDDO
      DO I=1,NX*NX
        C(I)=0.D0
      ENDDO
      RETURN
      END       
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_RIGID$BASIS(NAT,TMEMBER,R,IAT1_,IAT2_,IAT3_)
!     ******************************************************************
!     **  CONSTRAINTS_RIGID$BASIS                                     **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      LOGICAL(4),INTENT(IN) :: TMEMBER(NAT)
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      INTEGER(4),INTENT(OUT):: IAT1_
      INTEGER(4),INTENT(OUT):: IAT2_
      INTEGER(4),INTENT(OUT):: IAT3_
      REAL(8)               :: DETMAX,DET
      INTEGER(4)            :: IAT1,IAT2,IAT3
      REAL(8)               :: X21,Y21,Z21
      REAL(8)               :: X31,Y31,Z31
!     ******************************************************************
      DETMAX=0.D0
      IAT1_=0
      IAT2_=0
      IAT3_=0
      DO IAT1=1,NAT
        IF(TMEMBER(IAT1)) THEN
          DO IAT2=IAT1+1,NAT
            IF(TMEMBER(IAT2)) THEN
              X21=R(1,IAT2)-R(1,IAT1)
              Y21=R(2,IAT2)-R(2,IAT1) 
              Z21=R(3,IAT2)-R(3,IAT1)
              DO IAT3=IAT2+1,NAT
                IF(TMEMBER(IAT3)) THEN
                  X31=R(1,IAT3)-R(1,IAT1)
                  Y31=R(2,IAT3)-R(2,IAT1) 
                  Z31=R(3,IAT3)-R(3,IAT1)
                  DET=(Y21*Z31-Z21*Y31)**2 &
     &               +(Z21*X31-X21*Z31)**2 &
     &               +(X21*Y31-Y21*X31)**2
                  IF(DET.GT.DETMAX) THEN
                    IAT1_=IAT1
                    IAT2_=IAT2
                    IAT3_=IAT3
                    DETMAX=DET
                  END IF
                END IF
              ENDDO
            END IF
          ENDDO
        END IF
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_RIGID$VALUE(NAT,TMEMBER,R,IAT1_,IAT2_,IAT3_ &
     &                          ,NC,VALUE_)
!     ******************************************************************
!     **  CONSTRAINTS_RIGID$VALUE                                     **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      LOGICAL(4),INTENT(IN) :: TMEMBER(NAT)
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      INTEGER(4),INTENT(IN) :: IAT1_
      INTEGER(4),INTENT(IN) :: IAT2_
      INTEGER(4),INTENT(IN) :: IAT3_
      INTEGER(4),INTENT(IN) :: NC
      REAL(8)   ,INTENT(OUT):: VALUE_(NC)
      REAL(8)               :: X21,Y21,Z21
      REAL(8)               :: X31,Y31,Z31
      REAL(8)               :: X32,Y32,Z32
      INTEGER(4)            :: IC,IAT
      REAL(8)               :: XI1,YI1,ZI1
!     ******************************************************************
      X21=R(1,IAT2_)-R(1,IAT1_)
      Y21=R(2,IAT2_)-R(2,IAT1_) 
      Z21=R(3,IAT2_)-R(3,IAT1_)
      X31=R(1,IAT3_)-R(1,IAT1_)
      Y31=R(2,IAT3_)-R(2,IAT1_) 
      Z31=R(3,IAT3_)-R(3,IAT1_)
      X32=R(1,IAT3_)-R(1,IAT2_)
      Y32=R(2,IAT3_)-R(2,IAT2_) 
      Z32=R(3,IAT3_)-R(3,IAT2_)
      IC=0
      IC=IC+1
      VALUE_(IC)=X21*X21+Y21*Y21+Z21*Z21
      IC=IC+1
      VALUE_(IC)=X31*X31+Y31*Y31+Z31*Z31
      IC=IC+1
      VALUE_(IC)=X32*X32+Y32*Y32+Z32*Z32
      DO IAT=1,NAT
        IF(TMEMBER(IAT)) THEN
          IF(IAT.NE.IAT1_.AND.IAT.NE.IAT2_.AND.IAT.NE.IAT3_) THEN
            IF(IC+3.GT.NC) THEN
              CALL ERROR$MSG('NC SMALLER THAN EXPECTED')
              CALL ERROR$I4VAL('NC',NC)
              CALL ERROR$I4VAL('IC',IC)
              CALL ERROR$STOP('CONSTRAINTS_RIGID$VALUE')
            END IF
            XI1=R(1,IAT)-R(1,IAT1_)
            YI1=R(2,IAT)-R(2,IAT1_) 
            ZI1=R(3,IAT)-R(3,IAT1_)
            IC=IC+1
            VALUE_(IC)=XI1*X21+YI1*Y21+ZI1*Z21
            IC=IC+1
            VALUE_(IC)=XI1*X31+YI1*Y31+ZI1*Z31
            IC=IC+1
            VALUE_(IC)=XI1*(Y21*Z31-Z21*Y31) &
     &                +YI1*(Z21*X31-X21*Z31) &
     &                +ZI1*(X21*Y31-Y21*X31)
          END IF
        END IF
      ENDDO
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE CONSTRAINTS_RIGID$SET(NAT,TMEMBER,R,IAT1_,IAT2_,IAT3_ &
     &                       ,NC,A,B,C)
!     ******************************************************************
!     **  CONSTRAINTS_RIGID$SET                                       **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      LOGICAL(4),INTENT(IN) :: TMEMBER(NAT)
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      INTEGER(4),INTENT(IN) :: IAT1_
      INTEGER(4),INTENT(IN) :: IAT2_
      INTEGER(4),INTENT(IN) :: IAT3_
      INTEGER(4),INTENT(IN) :: NC
      REAL(8)   ,INTENT(OUT):: A(NC)
      REAL(8)   ,INTENT(OUT):: B(3,NAT,NC)
      REAL(8)   ,INTENT(OUT):: C(3,NAT,3,NAT,NC)
      INTEGER(4)            :: IAT1,IAT2,IAT3
      REAL(8)               :: X21,Y21,Z21
      REAL(8)               :: X31,Y31,Z31
!      REAL(8)               :: X32,Y32,Z32
      INTEGER(4)            :: IC,IAT,I,J
      REAL(8)               :: XI1,YI1,ZI1
!     ******************************************************************
      IAT1=IAT1_
      IAT2=IAT2_
      IAT3=IAT3_
      X21=R(1,IAT2_)-R(1,IAT1_)
      Y21=R(2,IAT2_)-R(2,IAT1_) 
      Z21=R(3,IAT2_)-R(3,IAT1_)
      X31=R(1,IAT3_)-R(1,IAT1_)
      Y31=R(2,IAT3_)-R(2,IAT1_) 
      Z31=R(3,IAT3_)-R(3,IAT1_)
!      X32=R(1,IAT3_)-R(1,IAT2_)
!      Y32=R(2,IAT3_)-R(2,IAT2_) 
!      Z32=R(3,IAT3_)-R(3,IAT2_)
!     
      IC=0
!     ================================================================
!     == BONDLENGTH CONSTRAINTS                                     ==
!     ================================================================
!     __(R2-R1)*(R2-R1)_______________________________________________
      IC=IC+1
      CALL CONSTRAINTS_INIT(3*NAT,A(IC),B(1,1,IC),C(1,1,1,1,IC))
      DO I=1,3
        C(I,IAT2,I,IAT2,IC)= 2.D0
        C(I,IAT2,I,IAT1,IC)=-2.D0
        C(I,IAT1,I,IAT2,IC)=-2.D0
        C(I,IAT1,I,IAT1,IC)= 2.D0
      ENDDO
!     __(R3-R1)*(R3-R1)_______________________________________________
      IC=IC+1
      CALL CONSTRAINTS_INIT(3*NAT,A(IC),B(1,1,IC),C(1,1,1,1,IC))
      DO I=1,3
        C(I,IAT3,I,IAT3,IC)= 2.D0
        C(I,IAT3,I,IAT1,IC)=-2.D0
        C(I,IAT1,I,IAT3,IC)=-2.D0
        C(I,IAT1,I,IAT1,IC)= 2.D0
      ENDDO
!     __(R3-R2)*(R3-R2)_______________________________________________
      IC=IC+1
      CALL CONSTRAINTS_INIT(3*NAT,A(IC),B(1,1,IC),C(1,1,1,1,IC))
      DO I=1,3
        C(I,IAT3,I,IAT3,IC)= 2.D0
        C(I,IAT3,I,IAT2,IC)=-2.D0
        C(I,IAT2,I,IAT3,IC)=-2.D0
        C(I,IAT2,I,IAT2,IC)= 2.D0
      ENDDO
!     
!     ================================================================
!     ==                                                            ==
!     ================================================================
      DO IAT=1,NAT
        IF(TMEMBER(IAT)) THEN
          IF(IAT.NE.IAT1.AND.IAT.NE.IAT2.AND.IAT.NE.IAT3) THEN
            XI1=R(1,IAT)-R(1,IAT1)
            YI1=R(2,IAT)-R(2,IAT1) 
            ZI1=R(3,IAT)-R(3,IAT1)
!           __(RI-R1)*(R2-R1)_________________________________________
            IC=IC+1
            CALL CONSTRAINTS_INIT(3*NAT,A(IC),B(1,1,IC),C(1,1,1,1,IC))
            DO I=1,3
              C(I,IAT,I,IAT2,IC) =+1.D0
              C(I,IAT2,I,IAT,IC) =+1.D0
              C(I,IAT,I,IAT1,IC) =-1.D0
              C(I,IAT1,I,IAT,IC) =-1.D0
              C(I,IAT2,I,IAT1,IC)=-1.D0
              C(I,IAT1,I,IAT2,IC)=-1.D0
              C(I,IAT1,I,IAT1,IC)=+2.D0
            ENDDO
!           __(RI-R1)*(R3-R1)_________________________________________
            IC=IC+1
            CALL CONSTRAINTS_INIT(3*NAT,A(IC),B(1,1,IC),C(1,1,1,1,IC))
            DO I=1,3
              C(I,IAT,I,IAT3,IC) =+1.D0
              C(I,IAT3,I,IAT,IC) =+1.D0
              C(I,IAT,I,IAT1,IC) =-1.D0
              C(I,IAT1,I,IAT,IC) =-1.D0
              C(I,IAT3,I,IAT1,IC)=-1.D0
              C(I,IAT1,I,IAT3,IC)=-1.D0
              C(I,IAT1,I,IAT1,IC)=+2.D0
            ENDDO
!           __DET[(RI-R1),(R3-R1),(R2-R1)]____________________________
            IC=IC+1
            CALL CONSTRAINTS_INIT(3*NAT,A(IC),B(1,1,IC),C(1,1,1,1,IC))
            A(IC)=XI1*(Y21*Z31-Z21*Y31) &
     &           +YI1*(Z21*X31-X21*Z31) &
     &           +ZI1*(X21*Y31-Y21*X31) 
!     
            B(1,IAT,IC) =Y21*Z31-Z21*Y31
            B(2,IAT,IC) =Z21*X31-X21*Z31
            B(3,IAT,IC) =X21*Y31-Y21*X31
            B(1,IAT2,IC)=Y31*ZI1-Z31*YI1
            B(2,IAT2,IC)=Z31*XI1-X31*ZI1
            B(3,IAT2,IC)=X31*YI1-Y31*XI1
            B(1,IAT3,IC)=YI1*Z21-ZI1*Y21
            B(2,IAT3,IC)=ZI1*X21-XI1*Z21
            B(3,IAT3,IC)=XI1*Y21-YI1*X21
            B(1,IAT1,IC)=-B(1,IAT,IC)-B(1,IAT2,IC)-B(1,IAT3,IC)
            B(2,IAT1,IC)=-B(2,IAT,IC)-B(2,IAT2,IC)-B(2,IAT3,IC)
            B(3,IAT1,IC)=-B(3,IAT,IC)-B(3,IAT2,IC)-B(3,IAT3,IC)
!     
            C(1,IAT,2,IAT2,IC)=+Z31
            C(1,IAT,3,IAT2,IC)=-Y31
            C(2,IAT,3,IAT2,IC)=+X31
            C(1,IAT,2,IAT3,IC)=-Z21
            C(1,IAT,3,IAT3,IC)=+Y21
            C(2,IAT,3,IAT3,IC)=-X21
            C(1,IAT2,2,IAT3,IC)=+ZI1
            C(1,IAT2,3,IAT3,IC)=-YI1
            C(2,IAT2,3,IAT3,IC)=+XI1
            DO I=1,3
              DO J=I+1,3
                C(J,IAT,I,IAT2,IC) =-C(I,IAT,J,IAT2,IC)
                C(J,IAT,I,IAT3,IC) =-C(I,IAT,J,IAT3,IC)
                C(J,IAT2,I,IAT3,IC)=-C(I,IAT2,J,IAT3,IC)
              ENDDO
            ENDDO
            DO I=1,3
              DO J=1,3
                C(I,IAT2,J,IAT,IC) = C(J,IAT,I,IAT2,IC)
                C(I,IAT3,J,IAT,IC) = C(J,IAT,I,IAT3,IC)
                C(I,IAT3,J,IAT2,IC)= C(J,IAT2,I,IAT3,IC)
                C(I,IAT,J,IAT1,IC) =-C(I,IAT,J,IAT2,IC) &
     &                              -C(I,IAT,J,IAT3,IC)
                C(I,IAT2,J,IAT1,IC)=-C(I,IAT2,J,IAT,IC) &
     &                              -C(I,IAT2,J,IAT3,IC)
                C(I,IAT3,J,IAT1,IC)=-C(I,IAT3,J,IAT,IC) &
     &                              -C(I,IAT3,J,IAT2,IC)
                C(J,IAT1,I,IAT,IC)= C(I,IAT,J,IAT1,IC)
                C(J,IAT1,I,IAT2,IC)=C(I,IAT2,J,IAT1,IC)
                C(J,IAT1,I,IAT3,IC)=C(I,IAT3,J,IAT1,IC)
                C(I,IAT1,J,IAT1,IC)=-C(I,IAT1,J,IAT,IC) &
     &                              -C(I,IAT1,J,IAT2,IC) &
     &                              -C(I,IAT1,J,IAT3,IC)
              ENDDO
            ENDDO
            CALL CONSTRAINTS_DERIVTOPOLYNOM(3*NAT,R &
     &                               ,A(IC),B(1,1,IC),C(1,1,1,1,IC))
!     
          END IF
        END IF
      ENDDO
      RETURN
!
      END
!     
!     ..................................................................
      SUBROUTINE CONSTRAINTS_DERIVTOPOLYNOM(N,R,A,B,C)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N
      REAL(8)   ,INTENT(IN)   :: R(N)
      REAL(8)   ,INTENT(INOUT):: A
      REAL(8)   ,INTENT(INOUT):: B(N)
      REAL(8)   ,INTENT(IN)   :: C(N,N)
      INTEGER(4)              :: I,J
!     ******************************************************************
      DO I=1,N
        A=A-B(I)*R(I)
        DO J=1,N
          A=A+0.5D0*R(I)*C(I,J)*R(J)
          B(I)=B(I)-C(I,J)*R(J)
        ENDDO
      ENDDO
      RETURN
      END
!     
!     ..................................................................
      SUBROUTINE CONSTRAINTS_TESTSET(N,R,A,B,C,TOL)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: R(N)
      REAL(8)   ,INTENT(IN) :: A
      REAL(8)   ,INTENT(IN) :: B(N)
      REAL(8)   ,INTENT(IN) :: C(N,N)
      REAL(8)   ,INTENT(IN) :: TOL
      REAL(8)               :: SVAR
      INTEGER(4)            :: I,J
!     ******************************************************************
      SVAR=A
      DO I=1,N
        SVAR=SVAR+B(I)*R(I)
        DO J=1,N
          SVAR=SVAR+0.5D0*R(I)*C(I,J)*R(J)
        ENDDO
      ENDDO
      IF(ABS(SVAR).GT.TOL) THEN
        WRITE(*,FMT='("TEST",3E15.5)')SVAR
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTRAINTS_APPLY(NX,X0,XP,XMASS,DELT,NC,A,B,C,RLAM)
!     ******************************************************************
!     **                                                              **
!     **   A+B*X+0.5*X*C*X=0                                          **
!     **                                                              **
!     **   CALCULATES VARIATIONS OF THE RESULT WITH CHANGES OF A(IC)  ** 
!     **   FOR A NUMBER (NDER) OF SPECIFIED CONSTRINTS                **
!     **                                                              **
!     ******************************************************************
      USE MPE_MODULE
      IMPLICIT NONE
      LOGICAL(4)               :: TSWITCH=.TRUE.
      REAL(8)   ,PARAMETER     :: TOL=1.D-9
      LOGICAL(4),PARAMETER     :: TPR=.FALSE.
      LOGICAL(4),PARAMETER     :: TTEST=.FALSE.
      INTEGER(4),PARAMETER     :: ITERX=100
      INTEGER(4),INTENT(IN)    :: NX
      INTEGER(4),INTENT(IN)    :: NC
      REAL(8)   ,INTENT(IN)    :: X0(NX)
      REAL(8)   ,INTENT(IN)    :: DELT
      REAL(8)   ,INTENT(IN)    :: XMASS(NX)
      REAL(8)   ,INTENT(IN)    :: A(NC)
      REAL(8)   ,INTENT(IN)    :: B(NX,NC)
      REAL(8)   ,INTENT(IN)    :: C(NX,NX,NC)
      REAL(8)   ,INTENT(INOUT) :: XP(NX)
      REAL(8)   ,INTENT(OUT)   :: RLAM(NC)
      REAL(8)                  :: G(NX,NC)
      REAL(8)                  :: SCALE(NC)
      REAL(8)                  :: AA(NC)
      REAL(8)                  :: BB(NC,NC)
      REAL(8)                  :: CC(NC,NC,NC)
      REAL(8)                  :: VEC(NC)
      REAL(8)                  :: RMAT(NC,NC)
      REAL(8)   ,ALLOCATABLE   :: WORK2D(:,:)
      REAL(8)   ,ALLOCATABLE   :: WORK3D(:,:,:)
      LOGICAL(4)               :: TCONV
      LOGICAL(4)               :: TLINDEP
      INTEGER(4)               :: THISTASK,NTASKS
      INTEGER(4)               :: IC,IC1,IC2,IC3,I,J,ICOUNT,ITER
      REAL(8)                  :: SVAR
!     ******************************************************************
!
!     ==================================================================
!     == PRINT                                                        ==
!     ==================================================================
      IF(TPR) THEN
        PRINT*,'INFORMATION FROM CONSTRAINTS_APPLY'
        PRINT*,'=================================='
        PRINT*,'NUMBER OF CONSTRAINTS        ',NC
        PRINT*,'NUMBER OF DEGREES OF FREEDOM ',NX
        PRINT*,'X0'
        WRITE(*,FMT='(6E14.3)')(X0(I),I=1,NX)
        PRINT*,'XP'
        WRITE(*,FMT='(6E14.3)')(XP(I),I=1,NX)
        PRINT*,'A'
        WRITE(*,FMT='(6E14.3)')(A(IC),IC=1,NC)
        PRINT*,'B'
        DO IC=1,NC
          WRITE(*,FMT='(6E14.3)')(B(I,IC),I=1,NX)
        ENDDO
        PRINT*,'C'
        DO IC=1,NC
          PRINT*,'CONSTRAINT ',IC
          DO I=1,NX
            DO J=1,NX
              IF(ABS(C(I,J,IC)).GT.1.D-6) &
     &           WRITE(*,FMT='("I,J",2I5," C ",F10.5)')I,J,C(I,J,IC)
            ENDDO
          ENDDO
        ENDDO
        PRINT*,'XMASS'
        WRITE(*,FMT='(6E14.3)')(XMASS(I),I=1,NX)
      END IF
!
!     ==================================================================
!     == FORCES OF CONSTRAINT (ALREADY DIVIDED BY THE MASS)           ==
!     ==================================================================
      DO IC=1,NC
        SVAR=0.D0
        DO I=1,NX
          G(I,IC)=B(I,IC)
          DO J=1,NX
            G(I,IC)=G(I,IC)+0.5D0*(C(I,J,IC)+C(J,I,IC))*X0(J)
          ENDDO
          G(I,IC)=G(I,IC)*DELT**2/XMASS(I)
          SVAR=SVAR+G(I,IC)**2
        ENDDO
        SCALE(IC)=1.D0/SQRT(SVAR)
        DO I=1,NX
          G(I,IC)=G(I,IC)*SCALE(IC)
        ENDDO
      ENDDO 
!
      IF(TPR) THEN
        PRINT*,'G'
        DO IC=1,NC
          WRITE(*,FMT='(6E14.3)')(G(I,IC),I=1,NX)
        ENDDO
        WRITE(*,FMT='("SCALE",6E14.3)')(SCALE(IC),IC=1,NC)
      END IF
!
!     ==================================================================
!     == CHECK LINEAR DEPENDENCY OF CONSTRAINT FORCES                 ==
!     == IF LINEAR DEPENDENCY IS DETECTED, SINGULAR VALUE             ==
!     == DECOMPOSITION IS USED TO SOLVE FOR LAGRANGE MULTIPLIERS      ==
!     ==================================================================
      TLINDEP=.FALSE.
      DO IC1=1,NC
        DO IC2=IC1+1,NC
          SVAR=0.D0
          DO I=1,NX
            SVAR=SVAR+G(I,IC1)*G(I,IC2)
          ENDDO
          IF(ABS(SVAR).GT.0.99D0) THEN
            TLINDEP=.TRUE.
            GOTO 90
          END IF
        ENDDO
      ENDDO
 90   CONTINUE
      IF(TLINDEP) PRINT*,'WARNING! CONSTRAINTS NEARLY LINEAR DEPENDENT'
!
!     ==================================================================
!     == EQUATION FOR LAMBDA                                          ==
!     ==  AA(N)+BB(N,I)*LAMBDA(I)+0.5*LAMBDA(I)*C(I,N,J)*LAMBDA(J)=0  ==
!     ==================================================================
IF(TSWITCH) THEN
!
!     ==================================================================
!     == EQUATION FOR LAMBDA                                          ==
!     ==  AA(N)+BB(N,I)*LAMBDA(I)+0.5*LAMBDA(I)*CC(I,N,J)*LAMBDA(J)=0 ==
!     ==  AA(N)=A(N)+XP(I)*B(I,N)+0.5*XP(I)*C(I,J,N)*XP(J)
!     ==  BB(N1,N2)=B(I,N1)*G(I,N2)+XP(I)*C(I,J,N1)*G(J,N2)
!     ==  CC(N2,N1,N3)=G(I,N2)*C(I,J,N1)*G(J,N3)
!     ==================================================================
      ALLOCATE(WORK3D(NX,NC,NC)) 
!     ==   WORK3D(J,N1,N2)=SUM_I C(I,J,N1)*G(I,N2)
      CALL LIB$SCALARPRODUCTR8(.FALSE.,NX,NX*NC,C,NC,G,WORK3D)
!     ==   CC(N3,N1,N2)=SUM_I G(I,N3)*WORK3D(I,N1,N2)
      CALL LIB$SCALARPRODUCTR8(.FALSE.,NX,NC,G,NC*NC,WORK3D,CC)
!    
!     ==   BB(N1,N2)=SUM_I XP(I)*WORK3D(I,N1,N2)
      CALL LIB$SCALARPRODUCTR8(.FALSE.,NX,1,XP,NC*NC,WORK3D,BB)
      DEALLOCATE(WORK3D)
!
      ALLOCATE(WORK2D(NC,NC))
!     ===  WORK2D(N1,N2)=SUM_I B(I,N1)*G(I,N2)
      CALL LIB$SCALARPRODUCTR8(.FALSE.,NX,NC,B,NC,G,WORK2D)
!     ===  BB(N1,N2)=BB(N1,N2)+WORK2D(N1,N2)
      BB(:,:)=BB(:,:)+WORK2D(:,:)
      DEALLOCATE(WORK2D)
!
      ALLOCATE(WORK2D(NX,NC))
!     == WORK3D(J,N1)=B(J,N1)+SUM_I XP(I)*C(I,J,N1)
      CALL LIB$SCALARPRODUCTR8(.FALSE.,NX,1,XP,NX*NC,C,WORK2D)
      WORK2D(:,:)=0.5D0*WORK2D(:,:)+B(:,:)
!     == AA=A+SUM_I XP(I)*WORK3D(I,N)
      CALL LIB$SCALARPRODUCTR8(.FALSE.,NX,1,XP,NC,WORK2D,AA)
      AA(:)=A(:)+AA(:)
      DEALLOCATE(WORK2D)
ELSE
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      ICOUNT=0
      DO IC=1,NC
        ICOUNT=ICOUNT+1
!       __ SELECTION FOR PARALLEL PROCESSING____________________________
        IF(MOD(ICOUNT-1,NTASKS).EQ.THISTASK-1) THEN
          SVAR=A(IC)
          DO I=1,NX
            SVAR=SVAR+B(I,IC)*XP(I)
            DO J=1,NX
              SVAR=SVAR+0.5D0*XP(I)*C(I,J,IC)*XP(J)
            ENDDO
          ENDDO
          AA(IC)=SVAR
!       __ELSE SELECTION FOR PARALLEL PROCESSING________________________
        ELSE
          AA(IC)=0.D0
!        __SELECTION OF PARALLELPROCESSING FINISHED______________________
        END IF
      ENDDO           
      CALL MPE$COMBINE('MONOMER','+',AA)
!
      ICOUNT=0
      DO IC1=1,NC
        DO IC2=1,NC
          ICOUNT=ICOUNT+1
!         __ SELECTION FOR PARALLEL PROCESSING__________________________
          IF(MOD(ICOUNT-1,NTASKS).EQ.THISTASK-1) THEN
            SVAR=0.D0
            DO I=1,NX
              SVAR=SVAR+B(I,IC1)*G(I,IC2)
              DO J=1,NX
                SVAR=SVAR+XP(I)*0.5D0*(C(I,J,IC1)+C(J,I,IC1))*G(J,IC2)
              ENDDO
            ENDDO  
            BB(IC1,IC2)=SVAR
!         __ELSE SELECTION FOR PARALLEL PROCESSING______________________
          ELSE
            BB(IC1,IC2)=0.D0
!         __SELECTION OF PARALLELPROCESSING FINISHED____________________
          END IF
        ENDDO
      ENDDO
      CALL MPE$COMBINE('MONOMER','+',BB)
! 
      ICOUNT=0
      DO IC1=1,NC
        DO IC2=1,NC
          DO IC3=1,NC
            ICOUNT=ICOUNT+1
!           __ SELECTION FOR PARALLEL PROCESSING________________________
            IF(MOD(ICOUNT-1,NTASKS).EQ.THISTASK-1) THEN
              SVAR=0.D0
              DO I=1,NX
                DO J=1,NX
                  SVAR=SVAR+G(I,IC1)*0.5D0*(C(I,J,IC2)+C(J,I,IC2))*G(J,IC3)
                ENDDO
              ENDDO
              CC(IC1,IC2,IC3)=SVAR
!           __ELSE SELECTION FOR PARALLEL PROCESSING____________________
            ELSE
              CC(IC1,IC2,IC3)=0.D0
!           __SELECTION OF PARALLELPROCESSING FINISHED__________________
            END IF
          ENDDO
        ENDDO
      ENDDO         
      CALL MPE$COMBINE('MONOMER','+',CC)
END IF
!
      IF(TPR) THEN
        PRINT*,'AA '
        WRITE(*,FMT='(7E10.2)')(AA(IC),IC=1,NC)
        PRINT*,'BB '
        DO IC1=1,NC
          WRITE(*,FMT='(7E10.2)')(BB(IC1,IC2),IC2=1,NC)
        ENDDO
        PRINT*,'CC '
        DO IC1=1,NC
        DO IC2=1,NC
          WRITE(*,FMT='(7E10.2)')(CC(IC1,IC2,IC3),IC3=1,NC)
        ENDDO
        ENDDO
      END IF
!
!     ==================================================================
!     == SET UP ITERATION FOR LAGRANGE MULTIPLIERS                    ==
!     ==================================================================
!
      DO IC1=1,NC
        DO IC2=1,NC
          RMAT(IC1,IC2)=BB(IC1,IC2)
        ENDDO
      ENDDO 
!
!      IF(NC.GT.1.AND.(.NOT.TLINDEP)) THEN
!        CALL DGEF(RMAT,NC,NC,IPVT)
!      END IF
!
!     ==================================================================
!     == ITERATE TO OBTAIN LAGRANGE MULTIPLIERS                       ==
!     ==================================================================
      DO IC=1,NC
        RLAM(IC)=0.D0
      ENDDO
!
      DO ITER=1,ITERX
!
        DO IC1=1,NC
          SVAR=AA(IC1)
          DO IC2=1,NC
            SVAR=SVAR+BB(IC1,IC2)*RLAM(IC2)
            DO IC3=1,NC
              SVAR=SVAR+0.5D0*RLAM(IC2)*CC(IC2,IC1,IC3)*RLAM(IC3)
            ENDDO
          ENDDO 
          VEC(IC1)=SVAR
        ENDDO
!
        IF(TPR) THEN
          WRITE(*,FMT='(72("="))')
          WRITE(*,FMT='("VEC IN ",10E9.2)')(VEC(IC1),IC1=1,NC)
        END IF
!
!       ==  TEST CONVERGENCE
        SVAR=0.D0
        DO IC=1,NC
          SVAR=MAX(SVAR,ABS(VEC(IC)))
        ENDDO
        TCONV=(SVAR.LT.TOL)
!       PRINT*,'DEVIATION IN CONSTRAINT',ABS(VEC)
        PRINT*,'DEVIATION IN CONSTRAINT',SVAR,ITER
        IF(TPR) PRINT*,'DEVIATION IN CONSTRAINT',SVAR
        IF(.NOT.(SVAR.LE.0.OR.SVAR.GT.0)) THEN
          CALL ERROR$MSG('DEVIATION VECTOR IS NOT A NUMERICAL QUANTITY')
          CALL ERROR$R8VAL('DEVIATION',SVAR)
          CALL ERROR$STOP('CONSTRAINTS_APPLY')
        END IF
        IF(TCONV) GOTO 100
!
!       == SOLVE VEC+RMAT*DRLAM=0 FOR DRLAM ============================
!       == DGESM('N',A,NC,NC,IPVT,B,NC,1) SOLVES A*X=B ; B=X ===========
        IF(NC.GT.1) THEN
!         __ CHECK LINEAR DEPENDENCE OF CONSTRAINTS_____________________
          IF(TLINDEP) THEN
            CALL LIB$MATRIXSOLVER8(NC,NC,1,RMAT,VEC,VEC)
          ELSE
            CALL LIB$MATRIXSOLVER8(NC,NC,1,RMAT,VEC,VEC)
          END IF
        ELSE
          VEC(1)=VEC(1)/RMAT(1,1)
        END IF
!
        IF(TPR) THEN
          WRITE(*,FMT='("VEC OUT",10E9.2)')(VEC(IC1),IC1=1,NC)
        END IF
 
!       == UPDATE RLAM =================================================
!       == VEC CONTAINS -D(RLAM)! ======================================
        DO IC=1,NC
          RLAM(IC)=RLAM(IC)-VEC(IC)
        ENDDO 
        IF(TPR) THEN
          WRITE(*,FMT='("RLAM   ",10E9.2)')(RLAM(IC1),IC1=1,NC)
        END IF
      ENDDO
      CALL ERROR$MSG('LOOP NOT CONVERGED')
      CALL ERROR$I4VAL('ITER',ITER)
      CALL ERROR$STOP('CONSTRAINTS_APPLY')
 100  CONTINUE         

!
!     ================================================================== 
!     ==  ADJUST XP FOR CONSTRAINTS                                   ==
!     ================================================================== 
!     == XP=XBAR+G*RLAM
      DO IC=1,NC
        DO I=1,NX
          XP(I)=XP(I)+G(I,IC)*RLAM(IC)
        ENDDO
      ENDDO
      IF(TPR) THEN
        PRINT*,'XP'
        WRITE(*,FMT='(6E14.3)')(XP(I),I=1,NX)
      END IF
!
!     ================================================================== 
!     ==  RESCALE LAGRANGE PARAMETERS                                 ==
!     ================================================================== 
      DO IC=1,NC
        RLAM(IC)=RLAM(IC)*SCALE(IC)
      ENDDO
!
!     ================================================================== 
!     ==  TEST CONSTRAINTS                                            ==
!     ================================================================== 
      IF(TTEST) THEN
        DO IC=1,NC
          SVAR=A(IC)
          DO I=1,NX
            SVAR=SVAR+XP(I)*B(I,IC)
            DO J=1,NX
              SVAR=SVAR+0.5D0*C(I,J,IC)*XP(I)*XP(J)
            ENDDO
          ENDDO
          PRINT*,'CONSTRAINT TEST ',IC,SVAR
          IF(TPR) THEN
            PRINT*,'CONSTRAINT TEST ',IC,SVAR
          END IF
          IF(ABS(SVAR).GT.TOL) THEN
            CALL ERROR$MSG('CONSTRAINT TEST FAILED')
            CALL ERROR$I4VAL('IC (CONSTRAINT NO.)',IC)
            CALL ERROR$R8VAL('ERROR',SVAR)
            CALL ERROR$STOP('CONSTRAINTS_APPLY')
          END IF
        ENDDO
      ENDIF
!
      RETURN
      END
!
!     ..................................................................
      MODULE SYMMETRIZE_MODULE
      LOGICAL(4) :: TOFF=.TRUE.
      INTEGER(4) :: NSYM
      INTEGER(4) :: ISYMDEF
      INTEGER(4) :: NAT
      REAL(8)   ,ALLOCATABLE :: ORIGIN(:,:) !(3,NSYM)
      REAL(8)   ,ALLOCATABLE :: RMAT(:,:,:) !(3,3,NSYM)
      INTEGER(4),ALLOCATABLE :: MAP(:,:)    !(NAT,NSYM)
      END MODULE SYMMETRIZE_MODULE
!
!     ..................................................................
      SUBROUTINE SYMMETRIZE$INITIALIZE(NSYM_,NAT_)
!     ******************************************************************
!     ******************************************************************
      USE SYMMETRIZE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSYM_
      INTEGER(4),INTENT(IN) :: NAT_
!     ******************************************************************
      TOFF=.TRUE.
      NSYM=NSYM_
      NAT=NAT_
      ALLOCATE(ORIGIN(3,NSYM))
      ALLOCATE(RMAT(3,3,NSYM))
      ALLOCATE(MAP(NAT,NSYM))
      ISYMDEF=0
      RETURN
      END 
!
!     ..................................................................
      SUBROUTINE SYMMETRIZE$SETSYMMETRY(RMAT_,ORIGIN_,NAT_,MAP_)
!     ******************************************************************
!     ******************************************************************
      USE SYMMETRIZE_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RMAT_(3,3)
      REAL(8)   ,INTENT(IN) :: ORIGIN_(3)
      INTEGER(4),INTENT(IN) :: NAT_
      INTEGER(4),INTENT(IN) :: MAP_(NAT_)
      INTEGER(4)            :: IAT,I,J
      INTEGER(4)            :: NFILO
!     ******************************************************************
      IF(.NOT.TOFF) THEN
        CALL ERROR$STOP('SYMMETRIZE$SETSYMMETRY')
      END IF
      IF(NAT_.NE.NAT) THEN
        CALL ERROR$STOP('SYMMETRIZE$SETSYMMETRY')
      END IF
      ISYMDEF=ISYMDEF+1   
      DO I=1,3
        ORIGIN(I,ISYMDEF)=ORIGIN_(I)
        DO J=1,3
          RMAT(I,J,ISYMDEF)=RMAT_(I,J)
        ENDDO
      ENDDO
      DO IAT=1,NAT
        MAP(IAT,ISYMDEF)=MAP_(IAT)
      ENDDO
      TOFF=(ISYMDEF.EQ.NSYM)
      IF(TOFF) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,FMT='(''NUMBER OF SYMMETRY OPERATIONS:'',I5)')NSYM
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SYMMETRIZE$FORCE(NAT_,F_)
!     ******************************************************************
!     ******************************************************************
      USE SYMMETRIZE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NAT_
      REAL(8)   ,INTENT(INOUT):: F_(3,NAT_)
      INTEGER(4)              :: IAT,IAT1,ISYM,I,J
      REAL(8)                 :: SVAR
      REAL(8)                 :: XNEW(3,NAT)
!     ******************************************************************
      IF(TOFF) RETURN
      DO IAT=1,NAT
        DO I=1,3
          XNEW(I,IAT)=0.D0
        ENDDO
      ENDDO
      DO ISYM=1,NSYM
        DO IAT=1,NAT
          DO I=1,3
            SVAR=0.D0
            DO J=1,3
              SVAR=SVAR+RMAT(I,J,ISYM)*F_(J,IAT)
            ENDDO
            IAT1=MAP(IAT,ISYM)
            XNEW(I,IAT1)=XNEW(I,IAT1)+SVAR/DBLE(NSYM)
          ENDDO
        ENDDO
      ENDDO
      DO IAT=1,NAT
        DO I=1,3
          F_(I,IAT)=XNEW(I,IAT)
        ENDDO
      ENDDO
      RETURN
      END
!       
!
!     ..................................................................
      SUBROUTINE SYMMETRIZE$POSITION(NAT_,R_)
!     ******************************************************************
!     ******************************************************************
      USE SYMMETRIZE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NAT_
      REAL(8)   ,INTENT(INOUT):: R_(3,NAT_)
      INTEGER(4)              :: IAT,IAT1,ISYM,I,J
      REAL(8)                 :: SVAR
      REAL(8)                 :: XNEW(3,NAT)
!     ******************************************************************
      IF(TOFF) RETURN
      IF(NAT_.NE.NAT) THEN
        CALL ERROR$STOP('SYMMETRIZE$POSITION')
      END IF
      DO IAT=1,NAT
        DO I=1,3
          XNEW(I,IAT)=0.D0
        ENDDO
      ENDDO
      DO ISYM=1,NSYM
        DO IAT=1,NAT
          DO I=1,3
            SVAR=0.D0
            DO J=1,3
              SVAR=SVAR+RMAT(I,J,ISYM)*(R_(J,IAT)-ORIGIN(J,ISYM))
            ENDDO
            IAT1=MAP(IAT,ISYM)
            SVAR=SVAR+ORIGIN(I,ISYM)
            XNEW(I,IAT1)=XNEW(I,IAT1)+SVAR
          ENDDO
        ENDDO
      ENDDO
      DO IAT=1,NAT
        DO I=1,3
          R_(I,IAT)=XNEW(I,IAT)/DBLE(NSYM)
        ENDDO
      ENDDO
      RETURN
      END

