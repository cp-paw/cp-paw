!      .................................................................
       SUBROUTINE STPPARMS_RESOLVE(LL_STRC_,LL_PARMS)
!      *****************************************************************
!      **  
!      *****************************************************************
       USE LINKEDLIST_MODULE
       IMPLICIT NONE
       TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
       TYPE(LL_TYPE)            :: LL_STRC 
       TYPE(LL_TYPE)            :: LL_PARMS
       INTEGER(4)               :: NSP
       INTEGER(4)               :: ISP
       LOGICAL(4)               :: TCHK
       LOGICAL(4)               :: TFOUND
       CHARACTER(32)            :: IDSTRC   ! SHORT NAME FOR THIS ATOM
       CHARACTER(32)            :: IDPARM   ! SHORT NAME FOR THIS ATOM
       CHARACTER(256)           :: STPROOT  ! ROOT NAME FOR SETUPS
       INTEGER(4)               :: NSETUP,ISETUP
       REAL(8)                  :: R8VAL
       CHARACTER(256)           :: CHVAL
       REAL(8)       ,ALLOCATABLE :: R8ARRAY(:)
!      *****************************************************************
       LL_STRC=LL_STRC_
       CALL LINKEDLIST$SELECT(LL_STRC,'~')
       CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
       CALL LINKEDLIST$NLISTS(LL_STRC,'SPECIES',NSP)
       DO ISP=1,NSP
         CALL LINKEDLIST$SELECT(LL_STRC,'SPECIES',ISP)
!
!        ===============================================================
!        ==  GET PARAMETER SHORT NAME IF PRESENT                      ==
!        ===============================================================
         CALL LINKEDLIST$EXISTD(LL_STRC,'PARMID',TCHK)
         IF(.NOT.TCHK) THEN
           CALL LINKEDLIST$SELECT(LL_STRC,'..')
           CYCLE
         END IF
         CALL LINKEDLIST$GET(LL_STRC,'PARMID',1,IDSTRC)
!
!        ===============================================================
!        ==  SEARCH FOR ENTRY IN PARAMETER FILE                       ==
!        ===============================================================
         CALL LINKEDLIST$SELECT(LL_PARMS,'~')
         CALL LINKEDLIST$SELECT(LL_PARMS,'STPPARMS')
         CALL LINKEDLIST$GET(LL_PARMS,'ROOT',STPROOT)
         CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NSETUP)
         TFOUND=.FALSE.
         DO ISETUP=1,NSETUP
           CALL LINKEDLIST$SELECT(LL_PARMS,'ATOM',ISETUP)
           CALL LINKEDLIST$GET(LL_PARMS,'PARMID',1,IDPARM)
           IF(IDSTRC.EQ.IDPARM) THEN
             TFOUND=.TRUE.
!
!            ===========================================================
!            == TRANSFER INFORMATION                                  ==
!            ===========================================================
!
!            == TRANSFER ATOMIC NUMBER =================================
             CALL LINKEDLIST$EXISTD(LL_STRC,'Z',1,TCHK)
             IF(.NOT.TCHK) THEN
               CALL LINKEDLIST$GET(LL_PARMS,'Z',1,R8VAL)      
               CALL LINKEDLIST$GET(LL_STRC,'Z',1,R8VAL)      
             END IF
!
!            == TRANSFER #(VALENCE ELECTRONS) ===========================
             CALL LINKEDLIST$GET(LL_PARMS,'ZV',1,R8VAL)      
             CALL LINKEDLIST$GET(LL_STRC,'ZV',1,R8VAL)      
!
!            == TRANSFER SETUP FILE =====================================
             CALL LINKEDLIST$EXISTD(LL_STRC,'FILE',1,TCHK)
             IF(.NOT.TCHK) THEN
               CALL ERROR$MSG('FILE MUST NOT BE PRESENT')
               CALL ERROR$STOP('STPPARMS_RESOLVE')
             END IF
             CALL LINKEDLIST$GET(LL_PARMS,'FILE',1,CHVAL)      
             CHVAL=TRIM(STPROOT)//TRIM(CHVAL)
             CALL LINKEDLIST$GET(LL_STRC,'FILE',1,TRIM(CHVAL))
! 
!            == TRANSFER NPRO ============================================
             CALL LINKEDLIST$EXISTD(LL_STRC,'NPRO',1,TCHK)
             IF(.NOT.TCHK) THEN
               CALL LINKEDLIST$SIZE(LL_PARMS,'NPRO',1,LENG)      
               ALLOCATE(I4ARRAY(LENG))
               CALL LINKEDLIST$GET(LL_PARMS,'NPRO',1,I4ARRAY)      
               CALL LINKEDLIST$SET(LL_STRC,'NPRO',1,I4ARRAY)
               DEALLOCATE(I4ARRAY)
             END IF
! 
!            == TRANSFER LRHOX ============================================
             CALL LINKEDLIST$EXISTD(LL_STRC,'LRHOX',1,TCHK)
             IF(.NOT.TCHK) THEN
               CALL LINKEDLIST$GET(LL_PARMS,'LRHOX',1,I4VAL)      
               CALL LINKEDLIST$SET(LL_STRC,'LRHOX',1,I4VAL)
             END IF
           END IF
           CALL LINKEDLIST$SELECT(LL_STRC,'..')
         ENDDO         
         IF(.NOT.TFOUND) THEN
           CALL ERROR$MSG('SETUP ENTRY IN PARAMETER FILE NOT FOUND')
           CALL ERROR$STOP('STPPARMS_RESOLVE')
         END IF
         CALL LINKEDLIST$SELECT(LL_STRC,'..')
       ENDDO
       RETURN
       END


