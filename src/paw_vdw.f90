!***************************************************************************
!**  van der waals interaction described by an aditional pair potential   **
!**                                                                       **
!**  Stefan Grimme, J. Comput.Chem. 25, 1463 (2004)                       **
!**                                                                       **
!**  (see S. Grimme, JCP132, 154104 (2010) for a more recent version)     **
!**                                                                       **
!***************************************************************************
MODULE VDW_MODULE
LOGICAL              :: TON=.FALSE.
LOGICAL              :: TINI=.FALSE.
LOGICAL              :: TPERIODIC=.TRUE.
REAL(8)   ,PARAMETER :: ALPHA=23
REAL(8)              :: S6=1.D0
CHARACTER(16)        :: FUNCTIONAL='NONE'
REAL(8)              :: C6FAC
REAL(8)              :: R0FAC
END MODULE VDW_MODULE
!
!       ...................................................................
        SUBROUTINE VDW$SETL4(ID,VAL)
        USE VDW_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID
        LOGICAL(4)  ,INTENT(IN) :: VAL
!       ********************************************************************
        IF(ID.EQ.'ON') THEN
          TON=VAL
        ELSE IF(ID.EQ.'PERIODIC') THEN
          TPERIODIC=VAL
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL(TRIM(ID),ID)
          CALL ERROR$STOP('VDW$SETL4')
        END IF
        RETURN
        END
!
!       ...................................................................
        SUBROUTINE VDW$GETL4(ID,VAL)
        USE VDW_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID
        LOGICAL(4)  ,INTENT(out) :: VAL
!       ********************************************************************
        IF(ID.EQ.'ON') THEN
          val=TON
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL(TRIM(ID),ID)
          CALL ERROR$STOP('VDW$GETL4')
        END IF
        RETURN
        END
!
!       ...................................................................
        SUBROUTINE VDW$SETCH(ID,VAL)
        USE VDW_MODULE
        IMPLICIT NONE
        CHARACTER(*),INTENT(IN) :: ID
        CHARACTER(*),INTENT(IN) :: VAL
!       ********************************************************************
        IF(ID.EQ.'FUNCTIONAL') THEN
          FUNCTIONAL=VAL
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL(TRIM(ID),ID)
          CALL ERROR$STOP('VDW$SETCH')
        END IF
        RETURN
        END
!
!       ...................................................................
        SUBROUTINE VDW_INITIALIZE()
        USE VDW_MODULE
        USE CONSTANTS_MODULE
        IMPLICIT NONE
        real(8)        :: svar
!       ********************************************************************
        IF(TINI) RETURN
!
        IF(FUNCTIONAL.EQ.'PBE') THEN
          S6=0.7D0
        ELSE IF(FUNCTIONAL.EQ.'BLYP') THEN
          S6=1.4D0
         ELSE IF(FUNCTIONAL.EQ.'BP86') THEN
          S6=1.3D0
        ELSE
          CALL ERROR$MSG('FUNCTIONAL NOT RECOGNIZED')
          CALL ERROR$CHVAL(TRIM(FUNCTIONAL),FUNCTIONAL)
          CALL ERROR$STOP('VDW_INUTIALIZE')
        END IF
!
        C6FAC=1.D0
        CALL CONSTANTS$GET('KJ/MOL',SVAR)
        C6FAC=C6FAC*SVAR*1.d-3
        CALL CONSTANTS$GET('METER',SVAR)
        C6FAC=C6FAC*SVAR**6
        CALL CONSTANTS$GET('NANO',SVAR)
        C6FAC=C6FAC*SVAR**6
!
        R0FAC=1.D0
        CALL CONSTANTS$GET('METER',SVAR)
        R0FAC=R0FAC*SVAR
        CALL CONSTANTS$GET('PICO',SVAR)
        R0FAC=R0FAC*SVAR

        TINI=.TRUE.
        RETURN
        END
!
!       ...................................................................
        SUBROUTINE VDW$INTERFACE(NAT,R,E,F)
        USE VDW_MODULE
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: NAT
        REAL(8)   ,INTENT(IN) :: R(3,NAT)
        REAL(8)   ,INTENT(OUT):: E
        REAL(8)   ,INTENT(OUT):: F(3,NAT)
        INTEGER(4)            :: IAT1,IAT2,IAT
        REAL(8)               :: DR(3),DR0(3)
        REAL(8)               :: R0(NAT)
        REAL(8)               :: C6(NAT)
        REAL(8)               :: R0PAIR
        REAL(8)               :: C6PAIR
        REAL(8)               :: E1
        REAL(8)               :: F1(3)
        REAL(8)               :: AEZ
        INTEGER(4)            :: IZ
        INTEGER(4)            :: ISP
        INTEGER(4)            :: It1min,it1max
        INTEGER(4)            :: It2min,it2max
        INTEGER(4)            :: It3min,it3max
        INTEGER(4)            :: It1,it2,it3
        real(8)               :: rbas(3,3)
!       ********************************************************************
        E=0.D0
        F(:,:)=0.D0
        IF(.NOT.TON) RETURN
        CALL VDW_INITIALIZE()
!
!       ====================================================================
!       ==  NOW CALCULATE TOTAL ENERGY                                    ==
!       ====================================================================
        CALL CELL$GETR8A('T(0)',9,RBAS)
        DO IAT=1,NAT
          CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
          CALL SETUP$AEZ(ISP,AEZ)
          IZ=NINT(AEZ)
!         __ GRIMME TABLE 1
          IF(IZ.EQ.1) THEN
            C6(IAT)=0.16D0
            R0(IAT)=111.D0
          ELSE IF(IZ.EQ.6) THEN
            C6(IAT)=1.65D0
            R0(IAT)=161.D0
          ELSE IF(IZ.EQ.7) THEN
            C6(IAT)=1.11D0
            R0(IAT)=155.D0
          ELSE IF(IZ.EQ.8) THEN
            C6(IAT)=0.70D0
            R0(IAT)=149.D0
          ELSE IF(IZ.EQ.9) THEN
            C6(IAT)=0.57D0
            R0(IAT)=143.D0
          ELSE IF(IZ.EQ.19) THEN
            C6(IAT)=0.45D0
            R0(IAT)=138.D0
          ELSE
            C6(IAT)=-1.D0
            R0(IAT)=0.D0
          END IF
        ENDDO
        C6(:)=C6(:)*C6FAC
        R0(:)=R0(:)*R0FAC
!
!       ====================================================================
!       ==  NOW CALCULATE TOTAL ENERGY                                    ==
!       ====================================================================
        IF(TPERIODIC) THEN
          it1min=-1
          it1max=1
          it2min=-1
          it2max=1
          it3min=-1
          it3max=1
        else
          it1min=0
          it1max=0
          it2min=0
          it2max=0
          it3min=0
          it3max=0
        END IF
        DO IAT1=1,NAT
          IF(C6(IAT1).LE.0.d0) CYCLE
          DO IAT2=IAT1,NAT
            IF(C6(IAT2).LE.0.d0) CYCLE
!
            DR0(:)=R(:,IAT2)-R(:,IAT1)
            R0PAIR=R0(IAT1)+R0(IAT2)
!           __ EQ.4 IN GRIMME
            C6PAIR=2.D0/(1.D0/C6(IAT1)+1.D0/C6(IAT2))
            do it1=it1min,it1max
              do it2=it2min,it2max
                do it3=it3min,it3max
                  if(iat1.eq.iat2) then
                    if(it1.eq.0.and.it2.eq.0.and.it3.eq.0) cycle
                  end if
                  dr(:)=dr0(:)+rbas(:,1)*real(it1,kind=8) &
      &                       +rbas(:,2)*real(it2,kind=8) &
      &                       +rbas(:,3)*real(it3,kind=8) 
                  CALL VDW_PAIR(DR,S6,C6PAIR,R0PAIR,ALPHA,E1,F1)
                  E=E+E1
                  F(:,IAT2)=F(:,IAT2)+F1(:)
                  F(:,IAT1)=F(:,IAT1)-F1(:)
                enddo
              enddo
            ENDDO
          ENDDO
        ENDDO
        RETURN
        END

!       ...................................................................
        SUBROUTINE VDW_PAIR(R,S6,C6,R0,ALPHA,E,F)
        IMPLICIT NONE
        real(8)   ,INTENT(IN) :: R(3)
        REAL(8)   ,INTENT(IN) :: S6
        REAL(8)   ,INTENT(IN) :: C6
        REAL(8)   ,INTENT(IN) :: R0
        REAL(8)   ,INTENT(IN) :: ALPHA
        REAL(8)   ,INTENT(OUT):: E
        REAL(8)   ,INTENT(OUT):: F(3)
        REAL(8)               :: D ! ACTUAL DISTANCE
        REAL(8)               :: FDMP
        REAL(8)               :: DFDMPBYFDMP
        REAL(8)               :: FAC !-1/D*(DE/DD)
!       ********************************************************************
        D=sqrt(SUM(R(:)**2))
!       == EQ.3 IN GRIMME
        FDMP=1.D0/(1.D0+EXP(-ALPHA*(D/R0-1.d0)))
        DFDMPBYFDMP=-FDMP*(-ALPHA/R0*EXP(-ALPHA*(D/R0-1.d0)))
!       == Eq.2 IN GRIMME
        E=-S6*C6/D**6*FDMP
        FAC=-(-6.D0*E/D+E*DFDMPBYFDMP)/D
        F(:)=FAC*R(:)
        RETURN
        END


     
