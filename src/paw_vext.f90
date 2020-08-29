! this module does not make sense and shall be removed....

MODULE VEXT_MODULE
!***********************************************************************
!**                                                                   **
!**                                                                   **
!**                                                                   **
!**                                                                   **
!**                                                                   **
!**                                                                   **
!***********************************************************************
TYPE VPAIR_TYPE
  INTEGER(4)               :: IAT1
  INTEGER(4)               :: IAT2
  REAL(8)                  :: E0
  TYPE(VPAIR_TYPE),POINTER :: NEXT
END TYPE VPAIR_TYPE
!
LOGICAL(4)               :: TON=.FALSE.
TYPE(VPAIR_TYPE),POINTER :: FIRSTVPAIR
CONTAINS
!      .................................................................
       SUBROUTINE VEXT_INITIALIZE
       IMPLICIT NONE
!      *****************************************************************
       TON=.TRUE.
       NULLIFY(FIRSTVPAIR)
       RETURN
       END SUBROUTINE VEXT_INITIALIZE
END MODULE VEXT_MODULE
!
!      .................................................................
       SUBROUTINE VEXT$SETUNBIND(IAT1,IAT2,E0)
!      *****************************************************************
!      **                                                             **
!      **                                                             **
!      **                                                             **
!      **                                                             **
!      **                                                             **
!      *****************************************************************
       USE VEXT_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)    :: IAT1
       INTEGER(4),INTENT(IN)    :: IAT2
       REAL(8)   ,INTENT(IN)    :: E0         
       TYPE(VPAIR_TYPE),POINTER :: XXX
!      *****************************************************************
       IF(.NOT.TON)CALL VEXT_INITIALIZE
       IF(.NOT.ASSOCIATED(FIRSTVPAIR)) THEN
         ALLOCATE(FIRSTVPAIR)
         XXX=>FIRSTVPAIR
       ELSE
         XXX=>FIRSTVPAIR
         DO WHILE(ASSOCIATED(XXX%NEXT))
           XXX=>XXX%NEXT
         ENDDO
         ALLOCATE(XXX%NEXT)
         XXX=>XXX%NEXT 
       END IF
       XXX%IAT1=IAT1
       XXX%IAT2=IAT2
       XXX%E0=E0
       NULLIFY(XXX%NEXT)
       RETURN
       END
!
!      .................................................................
       SUBROUTINE VEXT$APPLY
!      *****************************************************************
!      **                                                             **
!      **                                                             **
!      **                                                             **
!      **                                                             **
!      **                                                             **
!      *****************************************************************
       USE VEXT_MODULE
       IMPLICIT NONE
       INTEGER(4)                  :: NAT
       REAL(8)                     :: RBAS(3,3)
       REAL(8)         ,ALLOCATABLE:: R0(:,:)
       REAL(8)                     :: ENERGY
       REAL(8)         ,ALLOCATABLE:: FORCE(:,:)
       REAL(8)         ,ALLOCATABLE:: FORCE1(:,:)
!      *****************************************************************
       IF(.NOT.TON) RETURN
       CALL ATOMLIST$NATOM(NAT)
       ALLOCATE(R0(3,NAT))
       ALLOCATE(FORCE(3,NAT))
       ALLOCATE(FORCE1(3,NAT))
       CALL CELL$GETR8A('T(0)',9,RBAS)
       CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
       CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,FORCE)
       CALL VEXT_UNBIND(RBAS,NAT,R0,ENERGY,FORCE1)
       FORCE(:,:)=FORCE(:,:)+FORCE1(:,:)
       CALL ENERGYLIST$SET('EXTERNAL POTENTIAL',ENERGY)
       CALL ENERGYLIST$ADD('TOTAL ENERGY',ENERGY)
       CALL ATOMLIST$SETR8A('FORCE',0,3*NAT,FORCE)
       DEALLOCATE(R0)
       DEALLOCATE(FORCE)
       DEALLOCATE(FORCE1)
       RETURN
       END
!
!      .................................................................
       SUBROUTINE VEXT_UNBIND(RBAS,NAT,R0,ENERGY,FORCE)
!      *****************************************************************
!      **                                                             **
!      **  APPLIES A PAIR POTENTIAL BETWEEN SELECTED ATOM PAIRS       **
!      **  THAT VANISHES BEYOND THE SUM OF VANDER WAALS RADII,        **
!      **  AND INCREASES QUADRATICALLY FOR SHORTER DISTANCES.         **
!      **  AT THE  COVALENT BODN DISTANCE THE MAGNITUDE OF THE        **
!      **  POTENTIAL IS GIVEN.                                        **
!      **                                                             **
!      *****************************************************************
       USE VEXT_MODULE
       USE PERIODICTABLE_MODULE
       IMPLICIT NONE
       INTEGER(4)      ,INTENT(IN) :: NAT
       REAL(8)         ,INTENT(IN) :: RBAS(3,3)
       REAL(8)         ,INTENT(IN) :: R0(3,NAT)
       REAL(8)         ,INTENT(OUT):: ENERGY
       REAL(8)         ,INTENT(OUT):: FORCE(3,NAT)
       TYPE(VPAIR_TYPE),POINTER    :: VPAIR
       INTEGER(4)                  :: IAT1,IAT2
       REAL(8)                     :: E0
       REAL(8)                     :: DR(3)
       REAL(8)                     :: DIS
       REAL(8)                     :: aez
       REAL(8)                     :: DCOV,DVDW
       REAL(8)                     :: SVAR
       INTEGER(4)                  :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
       REAL(8)                     :: RVDW(NAT)
       REAL(8)                     :: RCOV(NAT)
       INTEGER(4)                  :: IAT,IT1,IT2,IT3,ISP
!      *****************************************************************
       ENERGY=0.D0
       FORCE(:,:)=0.D0
!
!      =================================================================
!      ==  COLLECT DATA                                               ==
!      =================================================================
       DO IAT=1,NAT
         CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
         CALL SETUP$ISELECT(ISP)
         CALL SETUP$GETR8('AEZ',AEZ)
         CALL SETUP$unSELECT()
         CALL PERIODICTABLE$GET(AEZ,'R(VDW)',RVDW(IAT))
         CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV(IAT))
       ENDDO
!
!      =================================================================
!      ==  DESTABILIZE COVALENT BONDS                                 ==
!      =================================================================
       IF(ASSOCIATED(FIRSTVPAIR)) THEN
         VPAIR=>FIRSTVPAIR
         DO 
           IAT1=VPAIR%IAT1
           IAT2=VPAIR%IAT2
           E0=VPAIR%E0
!
!          == CALCULATE VAN DER WAALS DISTANCE =========================
           DVDW=RVDW(IAT1)+RVDW(IAT2)
           DCOV=RCOV(IAT1)+RCOV(IAT2)
!
!          == FIND LATTICE DISPLACEMENTS ==============================
           DR(:)=R0(:,IAT2)-R0(:,IAT2)
           CALL BOXSPH(RBAS,-DR(1),-DR(2),-DR(3),DVDW &
       &               ,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
!
!          == CALCULATE DISTANCE =======================================
           DO IT1=MIN1,MAX1
             DO IT2=MIN2,MAX2
               DO IT3=MIN3,MAX3
                 DR(:)=R0(:,IAT2)-R0(:,IAT1) &
         &            +RBAS(:,1)*DBLE(IT1)+RBAS(:,2)*DBLE(IT2)+RBAS(:,3)*DBLE(IT3)
                 DIS=SQRT(DOT_PRODUCT(DR(:),DR(:)))
                 IF(DIS.LT.1.D-6) CYCLE
                 IF(DIS.LT.DVDW) THEN
                   SVAR=E0/(DCOV-DVDW)**2
                   ENERGY=ENERGY+SVAR*(DIS-DVDW)**2
                   SVAR=SVAR*2.D0*(DIS-DVDW)
                   FORCE(:,IAT1)=FORCE(:,IAT1)+SVAR*DR(:)/DIS
                   FORCE(:,IAT2)=FORCE(:,IAT2)-SVAR*DR(:)/DIS
                 END IF
               ENDDO
             ENDDO
           ENDDO
!
           IF(.NOT.ASSOCIATED(VPAIR%NEXT))EXIT
           VPAIR=>VPAIR%NEXT
         ENDDO
       ENDIF    
       RETURN
       END
