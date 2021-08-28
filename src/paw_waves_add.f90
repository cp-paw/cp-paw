!
!     ..................................................................
      SUBROUTINE WAVES$1CFORCECORR(SIGN)
!     ******************************************************************
!     **  CALCULATE 1-CENTER CONTRIBUTION TO THE FORCE                **
!     ******************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      REAL(8)     ,INTENT(IN):: SIGN 
      REAL(8)                :: OCC(NB,NKPT,NSPIN)        
      REAL(8)   ,ALLOCATABLE :: FORCE(:,:)   !(3,NAT)     
      REAL(8)   ,ALLOCATABLE :: FORCET(:,:)  !(3,NAT)     
      INTEGER(4),ALLOCATABLE :: ISPECIES(:)  !(NAT)       
      INTEGER(4),ALLOCATABLE :: LMNX(:)      !(NSP)       
      INTEGER(4)             :: ISP,I,IKPT,ISPIN,IAT
      INTEGER(4)             :: NTASKS,THISTASK
!     ******************************************************************
CALL ERROR$MSG('ROUTINE IS NOT FINISHED')
CALL ERROR$MSG('WATCH FOR MPE CALL TO "NONE"')
CALL ERROR$STOP('WAVES$1CFORCECORR')
                             CALL TRACE$PUSH('WAVES$1CFORCE')
      CALL SETUP$GETI4('NSP',NSP) !FORMER CALL SETUP$NSPECIES(NSP)
      ALLOCATE(LMNX(NSP))      
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LMNX',LMNX(ISP)) !FORMERCALL SETUP$LMNX(ISP,LMNX(ISP))
      ENDDO
      CALL SETUP$ISELECT(0)
      ALLOCATE(ISPECIES(NAT))  
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
!     == GET OCCUPATIONS =============================================
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
!     
!     ==================================================================
!     ==  CALCULATE FORCES DUE TO 1CENTER TERMS                       ==
!     ==================================================================
      ALLOCATE(FORCET(3,NAT)) 
      DO IAT=1,NAT
        DO I = 1,3
          FORCET(I,IAT)=0.D0
        ENDDO
      ENDDO
      CALL MPE$QUERY('NONE',NTASKS,THISTASK)
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IAT=THISTASK,NAT,NTASKS
            ISP=ISPECIES(IAT)              
            CALL WAVES_1CFORCECORR(IAT,NAT,NB,LMNXX &
     &                ,NB,LMNX(ISP),WKPT,OCC(:,IKPT,ISPIN) &
     &                ,RLAM0(:,:,IKPT,ISPIN) &
     &                ,DFNL(:,:,:,:,IKPT,ISPIN),HNL(:,:,:,IKPT,ISPIN) &
     &                ,ONL(:,:,:,IKPT,ISPIN),FORCET(:,IAT))
          ENDDO
        ENDDO
      ENDDO
      CALL MPE$COMBINE('NONE','+',FORCET)
!     
!     ==================================================================
!     ==  ADD THESE FORCES TO ATOMLIST                                ==
!     ==================================================================
      ALLOCATE(FORCE(3,NAT))   
      CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,FORCE)
      DO IAT=1,NAT
        DO I=1,3
          FORCE(I,IAT)=FORCE(I,IAT)+SIGN*FORCET(I,IAT)
        ENDDO
      ENDDO
      CALL ATOMLIST$SETR8A('FORCE',0,3*NAT,FORCE)
!     
!     ==================================================================
!     ==  PRINT FOR TEST                                              ==
!     ==================================================================
      IF(TPR) THEN
        PRINT*,'FORCE AFTER WAVES_1CFORCE' 
        DO IAT=1,NAT
          WRITE(*,FMT='(3F20.10)')(FORCE(I,IAT),I=1,3)
        ENDDO
      END IF
!     
!     ==================================================================
!     ==  CLOSE ROUTINE                                               ==
!     ==================================================================
      DEALLOCATE(FORCE)    
      DEALLOCATE(FORCET)   
      DEALLOCATE(ISPECIES) 
      DEALLOCATE(LMNX)     
                           CALL TRACE$POP
      RETURN
      END
!
!     .....................................................INVOV........
      SUBROUTINE WAVES_1CFORCECORR(IAT,NAT,NX,LMNXX &
     &                     ,NB,LMNX,WKPT,F,H,DFNL,HNL,ONL,FION)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE FORCE DUE TO THE POSITION DEPENDENCE         **
!     **  OF THE OVERLAP OPERATOR                                     **
!     **                                                              **
!     *******************************************P.E. BLOECHL, (1992)***
      IMPLICIT NONE
      INTEGER(4),PARAMETER  :: IPR=0
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: NX
      INTEGER(4),INTENT(IN) :: LMNXX
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: LMNX
      REAL(8)   ,INTENT(IN) :: WKPT
      REAL(8)   ,INTENT(IN) :: F(NB)
      REAL(8)   ,INTENT(IN) :: H(NX,NX)   ! ORTHO-LAGRANGE MULTIPLIER
      REAL(8)   ,INTENT(IN) :: DFNL(NAT,NB,LMNXX,3)
      REAL(8)   ,INTENT(IN) :: HNL(NAT,NB,LMNXX)  
      REAL(8)   ,INTENT(IN) :: ONL(NAT,NB,LMNXX)
      REAL(8)   ,INTENT(INOUT):: FION(3)
      INTEGER(4)            :: IDIR,IB,LMN,IB1,IB2,I
      REAL(8)               :: FAC,SVAR
      REAL(8)               :: F1,F2
!     ******************************************************************
!
!     ==================================================================
!     ==  PRINTOUT FOR TEST                                           ==
!     ==================================================================
      IF(IPR.EQ.1) THEN
        PRINT*,'OVFOR START'
        WRITE(*,FMT='("FORCE ON ATOM ",I3,":",3F15.5)') &
     &               IAT,(FION(I),I=1,3)
      ENDIF
!
!     ==================================================================
!     ==  2*TR[ <PSI|P>DO<GRAD*P|PSI> * LAMBDA ]                      ==
!     ==================================================================
      DO LMN=1,LMNX
        DO IB1=1,NB
          F1=F(IB1)
          SVAR=0.D0
          DO IB2=1,NB
            F2=F(IB2)
            SVAR=SVAR + H(IB1,IB2)* ONL(IAT,IB2,LMN) &
     &                * 2.D0*F1*F2/(F1+F2+1.D-12)
          ENDDO
          SVAR=SVAR*2.D0*WKPT
          FION(1)=FION(1)+DFNL(IAT,IB1,LMN,1)*SVAR
          FION(2)=FION(2)+DFNL(IAT,IB1,LMN,2)*SVAR
          FION(3)=FION(3)+DFNL(IAT,IB1,LMN,3)*SVAR
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  PRINTOUT FOR TEST                                           ==
!     ==================================================================
      IF(IPR.EQ.1) THEN
        PRINT*,'AFTER OVFOR'
        WRITE(*,FMT='("FORCE ON ATOM ",I3,":",3F15.5)') &
     &               IAT,(FION(I),I=1,3)
      ENDIF 
      RETURN
      END
      

