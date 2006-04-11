!***********************************************************************
!***********************************************************************
!**                                                                   **
!**  INTERFACES TO SCIENTIFIC LIBRARY ROUTINES                        **
!**  USES THE ESSL LIBRARY                                            **
!**                                                                   **
!***********************************************************************
!***********************************************************************
! CPPVAR_FFTW      USE FFTW FOR FOURIRT TRANSFORMS
! CPPVAR_FFT_ESSL      USE ESSL FOR FOURIER TRANSFORMS
! CPPVAR_FFTPACK      USE EXPLICIT FFT
!
! CPPVAR_BLAS_ATLAS    USE ATLAS BLAS
! CPPVAR_BLAS_ESSL     USE ESSL
! CPPVAR_BLAS_EXPLICIT USE EXPLICIT BLAS
!
! CPPVAR_LANGEXT_XLF   LANGUAGE EXTENSIONS XLF
!
! CPPVAR_SUPPORT_XLF
! CPPVAR_U77
! CPPVAR_SUPPORT_DEC
!
! CPPVAR_ERF_EXPLICIT  USE EXPLICIT ERF AND ERFC
! CPPVAR_ERF_IMPLICIT  USE IMPLICIT ERF AND ERFC
!
! CPPVAR_USAGE_EXIST  C-ROUTINE GETRUSAGE AVAILABLE
!
!=============== ENVIRONMENT ABSOFT COMPILER =========================
!    IT ASSUMES THE ATLAS LAPACK AND BLAS ROUTINES TO BE LINKED
!    IT ASSUMES THE ABSOFT SUPPORT LIBRARY LIBU77.A  TO BE LINKED
!    IT ASSUMES THE FAST FOURIER TRANSFORM LIBRARY LIBFFTW.A  TO BE LINKED
!#IF DEFINED(CPPVARIABLE_ABS)
!#DEFINE CPPVAR_FFTW
!#DEFINE CPPVAR_BLAS_ATLAS
!#DEFINE CPPVAR_U77
!=============== ENVIRONMENT IBM AIX ===================================
!    IT ASSUMES THE ESSL LIBRARY TO BE LINKED
!    IT USES XLF SPECIFIC SUPPORT ROUTINES AND LANGUAGE EXTENSIONS
!#ELIF DEFINED(CPPVARIABLE_XLF)
!#DEFINE CPPVAR_FFT_ESSL
!#DEFINE CPPVAR_BLAS_ESSL 
!#DEFINE CPPVAR_SUPPORT_XLF
!#DEFINE CPPVAR_LANGEXT_XLF  
!#DEFINE CPPVAR_USAGE_EXIST
!#UNDEFINE EXLICITERF
!=============== ENVIRONMENT DEC ALPHA =================================
!#ELIF DEFINED(CPPVARIABLE_DEC)
!#DEFINE CPPVAR_FFTW
!#DEFINE CPPVAR_BLAS_ATLAS
!=============== FREE ENVIRONMENT ===== =================================
! MAKES NO ASSUMPTIONS ABOUT THE ENVIRONMENT
!#ELSE
!#DEFINE CPPVAR_FFTPACK
!#DEFINE CPPVAR_BLAS_EXPLICIT
!#ENDIF 
!#IF DEFINED(IBMLICENSE)
!
!     .................................................................
      SUBROUTINE LIB$GETUSAGE(ID,VALUE)
!     *****************************************************************
!     **                                                             **
!     **  PROVIDES INFORMATION  ON THE USAGE OF SYSTEM RESOURCES     **
!     **                                                             **
!     **  USES STANDARD C-LIBRARY ROUTINE GETRUSAGE                  **
!     **                                                             **
!     **                                                             **
!     **                                                             **
!     *****************************************************************
      IMPLICIT NONE
      TYPE USG_TYPE 
        SEQUENCE
        INTEGER(4) :: UTIME(2)   ! USER TIME    (SECONDS,MICROSECONDS)
        INTEGER(4) :: STIME(2)   ! SYSTEM TIME  (SECONDS,MICROSECONDS)
        INTEGER(4) :: MAXRSS     ! MAX SIZE [KBYTE] 
        INTEGER(4) :: IXRSS      ! INTEGRAL SHARED MEMORY SIZE [KBYTE*SEC]         
        INTEGER(4) :: IDRSS      ! INTEGRAL UNSHARED DATA [KBYTE*SEC]    
        INTEGER(4) :: ISRSS      ! INTEGRAL UNSHARED STACK [KBYTE*SEC]    
        INTEGER(4) :: MINFLT     ! #(PAGE RECLAIMS)
        INTEGER(4) :: MAJFLT     ! #(PAGE FAULTS)
        INTEGER(4) :: NSWAP      ! #(SWAPPING PROCESS OUT OF MAIN)            
        INTEGER(4) :: INBLOCK    ! #(FILE INPUTS)                             
        INTEGER(4) :: OUBLOCK    ! #(FILE OUTPUTS)                            
        INTEGER(4) :: MSGSND     ! #(IPC MESSAGES SEND)                       
        INTEGER(4) :: MSGRCV     ! #(IPC MESSAGES RECEIVED)                   
        INTEGER(4) :: NSIGNALS   ! #(SIGNALS SENT)                            
        INTEGER(4) :: NVCSW      ! #(VOLUNTARY CONTEXT SWITCHES)
        INTEGER(4) :: NIVCSW     ! #(INVOLUNTARY CONTEXT SWITCHES)
      END TYPE USG_TYPE
      CHARACTER(*),INTENT(IN)  :: ID
      REAL(8)     ,INTENT(OUT) :: VALUE
      TYPE(USG_TYPE)           :: USG
      REAL(8)                  :: KBYTE=2.D0**10
      INTEGER(4)               :: GETRUSAGE
      INTEGER(4)               :: RC
      REAL(8)                  :: CPUTIME
#IF DEFINED(CPPVAR_USAGE_EXIST)
      EXTERNAL GETRUSAGE
#ENDIF
!     ******************************************************************
      USG%UTIME=(/0,0/)
      USG%STIME=(/0,0/)
      USG%MAXRSS=0     
      USG%IXRSS=0      
      USG%IDRSS=0      
      USG%ISRSS=0      
      USG%MINFLT=0     
      USG%MAJFLT=0     
      USG%NSWAP=0      
      USG%INBLOCK=0    
      USG%OUBLOCK=0    
      USG%MSGSND=0     
      USG%MSGRCV=0     
      USG%NSIGNALS=0   
      USG%NVCSW=0      
      USG%NIVCSW=0     
!     ==================================================================
!     ==================================================================
#IF DEFINED(CPPVAR_USAGE_EXIST)
      RC=GETRUSAGE(%VAL(0),USG)    ! C-ROUTINE
      IF(RC.NE.0) THEN
        CALL ERROR$MSG('ERROR CALLING MYGETRUSAGE')
        CALL ERROR$STOP('MEMORY$GET')
      END IF
#ENDIF
!
      CPUTIME=(REAL(USG%UTIME(1)+USG%STIME(1),KIND=8) &
     &        +REAL(USG%UTIME(2)+USG%STIME(2),KIND=8))*1.D-6
      CPUTIME=MAX(CPUTIME,1.D-6)
      IF(ID.EQ.'MAXMEM') THEN
        VALUE=REAL(USG%MAXRSS,KIND=8)*KBYTE
      ELSE IF(ID.EQ.'USRTIME') THEN
        VALUE=(REAL(USG%UTIME(1),KIND=8)+REAL(USG%UTIME(2),KIND=8))*1.D-6
      ELSE IF(ID.EQ.'SYSTIME') THEN
        VALUE=(REAL(USG%STIME(1),KIND=8)+REAL(USG%STIME(2),KIND=8))*1.D-6
      ELSE IF(ID.EQ.'CPUTIME') THEN
        VALUE=CPUTIME
      ELSE IF(ID.EQ.'PAGEFAULTRATE') THEN
        VALUE=REAL(USG%MAJFLT,KIND=8)/CPUTIME
      ELSE IF(ID.EQ.'SWAPRATE') THEN
        VALUE=REAL(USG%NSWAP,KIND=8)/CPUTIME
      ELSE IF(ID.EQ.'SWITCHRATE') THEN
        VALUE=REAL(USG%NVCSW+USG%NIVCSW,KIND=8)/CPUTIME
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('MEMORY$GET')
      END IF
      RETURN
      END
#IF DEFINED(CPPVAR_LANGEXT_XLF)
!     ..................................................................
      SUBROUTINE LIB$ERFR8(X,Y)
!     ******************************************************************
!     **  ERROR FUNCTION                                              **
!     **    Y=2/SQRT(PI)  INT_0^X DZ EXP(-Z**2)                       **
!     **    Y=(INFINITY)=1                                            **
!     ******************************************************************
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: Y
!     ******************************************************************
      Y=DERF(X) 
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$ERFCR8(X,Y)
!     ****************************************************************** 
!     **  COMPLEMENTARY ERROR FUNCTION                                **
!     **  Y=1-ERF(X)                                                  **
!     ******************************************************************
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: Y
!     ******************************************************************
      Y=DERFC(X)
      RETURN
      END
#ELSE
!
!     ..................................................................
      SUBROUTINE LIB$ERFR8(X,Y)
!     ******************************************************************
!     **  COPYRIGHT(C) 1996 TAKUYA OOURA                              **
!     **  (EMAIL: OOURA@MMM.T.U-TOKYO.AC.JP).                         **
!     **  YOU MAY USE, COPY, MODIFY THIS CODE FOR ANY PURPOSE AND     **
!     **  WITHOUT FEE. YOU MAY DISTRIBUTE THIS ORIGINAL PACKAGE.        **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: Y
      REAL(8)            :: W
      REAL(8)            :: T
      INTEGER(4)         :: K,I
      REAL(8)            :: A(0:64)
      REAL(8)            :: B(0:64)
!     ******************************************************************
      DATA (A(I), I = 0, 12) / &
     &    0.00000000005958930743D0, -0.00000000113739022964D0, &
     &    0.00000001466005199839D0, -0.00000016350354461960D0, &
     &    0.00000164610044809620D0, -0.00001492559551950604D0, &
     &    0.00012055331122299265D0, -0.00085483269811296660D0, &
     &    0.00522397762482322257D0, -0.02686617064507733420D0, &
     &    0.11283791670954881569D0, -0.37612638903183748117D0, &
     &    1.12837916709551257377D0 / 
      DATA (A(I), I = 13, 25) /                                &
     &    0.00000000002372510631D0, -0.00000000045493253732D0, &
     &    0.00000000590362766598D0, -0.00000006642090827576D0, &
     &    0.00000067595634268133D0, -0.00000621188515924000D0, &
     &    0.00005103883009709690D0, -0.00037015410692956173D0, &
     &    0.00233307631218880978D0, -0.01254988477182192210D0, &
     &    0.05657061146827041994D0, -0.21379664776456006580D0, &
     &    0.84270079294971486929D0 / 
      DATA (A(I), I = 26, 38) /                                &
     &    0.00000000000949905026D0, -0.00000000018310229805D0, &
     &    0.00000000239463074000D0, -0.00000002721444369609D0, &
     &    0.00000028045522331686D0, -0.00000261830022482897D0, &
     &    0.00002195455056768781D0, -0.00016358986921372656D0, &
     &    0.00107052153564110318D0, -0.00608284718113590151D0, &
     &    0.02986978465246258244D0, -0.13055593046562267625D0, &
     &    0.67493323603965504676D0 / 
      DATA (A(I), I = 39, 51) /                                &
     &    0.00000000000382722073D0, -0.00000000007421598602D0, &
     &    0.00000000097930574080D0, -0.00000001126008898854D0, &
     &    0.00000011775134830784D0, -0.00000111992758382650D0, &
     &    0.00000962023443095201D0, -0.00007404402135070773D0, &
     &    0.00050689993654144881D0, -0.00307553051439272889D0, &
     &    0.01668977892553165586D0, -0.08548534594781312114D0, &
     &    0.56909076642393639985D0 / 
      DATA (A(I), I = 52, 64) /                                &
     &    0.00000000000155296588D0, -0.00000000003032205868D0, &
     &    0.00000000040424830707D0, -0.00000000471135111493D0, &
     &    0.00000005011915876293D0, -0.00000048722516178974D0, &
     &    0.00000430683284629395D0, -0.00003445026145385764D0, &
     &    0.00024879276133931664D0, -0.00162940941748079288D0, &
     &    0.00988786373932350462D0, -0.05962426839442303805D0, &
     &    0.49766113250947636708D0 / 
      DATA (B(I), I = 0, 12) /                                 &
     &   -0.00000000029734388465D0,  0.00000000269776334046D0, &
     &   -0.00000000640788827665D0, -0.00000001667820132100D0, & 
     &   -0.00000021854388148686D0,  0.00000266246030457984D0, &
     &    0.00001612722157047886D0, -0.00025616361025506629D0, &
     &    0.00015380842432375365D0,  0.00815533022524927908D0, &
     &   -0.01402283663896319337D0, -0.19746892495383021487D0, & 
     &    0.71511720328842845913D0 / 
      DATA (B(I), I = 13, 25) /                                 &
     &   -0.00000000001951073787D0, -0.00000000032302692214D0,  &
     &    0.00000000522461866919D0,  0.00000000342940918551D0,  &
     &   -0.00000035772874310272D0,  0.00000019999935792654D0,  &
     &    0.00002687044575042908D0, -0.00011843240273775776D0,  &
     &   -0.00080991728956032271D0,  0.00661062970502241174D0,  &
     &    0.00909530922354827295D0, -0.20160072778491013140D0,  &
     &    0.51169696718727644908D0 / 
      DATA (B(I), I = 26, 38) /                                 &
     &    0.00000000003147682272D0, -0.00000000048465972408D0,  &
     &    0.00000000063675740242D0,  0.00000003377623323271D0,  &
     &   -0.00000015451139637086D0, -0.00000203340624738438D0,  &
     &    0.00001947204525295057D0,  0.00002854147231653228D0,  &
     &   -0.00101565063152200272D0,  0.00271187003520095655D0,  &
     &    0.02328095035422810727D0, -0.16725021123116877197D0,  &
     &    0.32490054966649436974D0 / 
      DATA (B(I), I = 39, 51) /                                &
     &    0.00000000002319363370D0, -0.00000000006303206648D0, & 
     &   -0.00000000264888267434D0,  0.00000002050708040581D0, & 
     &    0.00000011371857327578D0, -0.00000211211337219663D0, & 
     &    0.00000368797328322935D0,  0.00009823686253424796D0, & 
     &   -0.00065860243990455368D0, -0.00075285814895230877D0, & 
     &    0.02585434424202960464D0, -0.11637092784486193258D0, & 
     &    0.18267336775296612024D0 / 
      DATA (B(I), I = 52, 64) /                                &
     &   -0.00000000000367789363D0,  0.00000000020876046746D0, & 
     &   -0.00000000193319027226D0, -0.00000000435953392472D0, & 
     &    0.00000018006992266137D0, -0.00000078441223763969D0, & 
     &   -0.00000675407647949153D0,  0.00008428418334440096D0, & 
     &   -0.00017604388937031815D0, -0.00239729611435071610D0, & 
     &    0.02064129023876022970D0, -0.06905562880005864105D0, & 
     &    0.09084526782065478489D0 / 
      W = ABS(X)
      IF (W .LT. 2.2D0) THEN
          T = W * W
          K = INT(T)
          T = T - REAL(K,KIND=8)
          K = K * 13
          Y = (((((((((((( A(K    )*T + A(K+ 1)) * T + &
     &        A(K+ 2))*T + A(K+ 3))*T + A(K+ 4)) * T + &
     &        A(K+ 5))*T + A(K+ 6))*T + A(K+ 7)) * T + &
     &        A(K+ 8))*T + A(K+ 9))*T + A(K+10)) * T + &
     &        A(K+11))*T + A(K+12))*W
      ELSE IF (W .LT. 6.9D0) THEN
          K = INT(W)
          T = W - REAL(K,KIND=8)
          K = 13 * (K - 2)
          Y = (((((((((((B(K) * T + B(K + 1)) * T + &
     &        B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T + &
     &        B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T + &
     &        B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T + &
     &        B(K + 11)) * T + B(K + 12)
          Y = Y * Y
          Y = Y * Y
          Y = Y * Y
          Y = 1 - Y * Y
      ELSE
          Y = 1
      END IF
      IF(X.LT.0) Y=-Y
      RETURN
      END
!     .......................................................................
      SUBROUTINE LIB$ERFCR8(X,Y)
!     **  COPYRIGHT(C) 1996 TAKUYA OOURA                              **
!     **  (EMAIL: OOURA@MMM.T.U-TOKYO.AC.JP).                         **
!     **  YOU MAY USE, COPY, MODIFY THIS CODE FOR ANY PURPOSE AND     **
!     **  WITHOUT FEE. YOU MAY DISTRIBUTE THIS ORIGINAL PACKAGE.        **
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: Y
      REAL(8)     ,PARAMETER :: PA  = 3.97886080735226000D+00
      REAL(8)     ,PARAMETER :: P0  = 2.75374741597376782D-01
      REAL(8)     ,PARAMETER :: P1  = 4.90165080585318424D-01 
      REAL(8)     ,PARAMETER :: P2  = 7.74368199119538609D-01 
      REAL(8)     ,PARAMETER :: P3  = 1.07925515155856677D+00 
      REAL(8)     ,PARAMETER :: P4  = 1.31314653831023098D+00 
      REAL(8)     ,PARAMETER :: P5  = 1.37040217682338167D+00 
      REAL(8)     ,PARAMETER :: P6  = 1.18902982909273333D+00 
      REAL(8)     ,PARAMETER :: P7  = 8.05276408752910567D-01 
      REAL(8)     ,PARAMETER :: P8  = 3.57524274449531043D-01 
      REAL(8)     ,PARAMETER :: P9  = 1.66207924969367356D-02 
      REAL(8)     ,PARAMETER :: P10 =-1.19463959964325415D-01 
      REAL(8)     ,PARAMETER :: P11 =-8.38864557023001992D-02
      REAL(8)     ,PARAMETER :: P12 = 2.49367200053503304D-03 
      REAL(8)     ,PARAMETER :: P13 = 3.90976845588484035D-02 
      REAL(8)     ,PARAMETER :: P14 = 1.61315329733252248D-02 
      REAL(8)     ,PARAMETER :: P15 =-1.33823644533460069D-02 
      REAL(8)     ,PARAMETER :: P16 =-1.27223813782122755D-02 
      REAL(8)     ,PARAMETER :: P17 = 3.83335126264887303D-03 
      REAL(8)     ,PARAMETER :: P18 = 7.73672528313526668D-03 
      REAL(8)     ,PARAMETER :: P19 =-8.70779635317295828D-04 
      REAL(8)     ,PARAMETER :: P20 =-3.96385097360513500D-03 
      REAL(8)     ,PARAMETER :: P21 = 1.19314022838340944D-04 
      REAL(8)     ,PARAMETER :: P22 = 1.27109764952614092D-03
      REAL(8)                :: T,U
!     **********************************************************
      T = PA / (PA + ABS(X))
      U = T - 0.5D0
      Y = (((((((((P22 * U + P21) * U + P20) * U + &
     &    P19) * U + P18) * U + P17) * U + P16) * U + &
     &    P15) * U + P14) * U + P13) * U + P12 
      Y = ((((((((((((Y * U + P11) * U + P10) * U + &
     &    P9) * U + P8) * U + P7) * U + P6) * U + P5) * U + &
     &    P4) * U + P3) * U + P2) * U + P1) * U + P0) * T * &
     &    EXP(-X * X)
      IF (X .LT. 0) Y = 2 - Y
      RETURN
      END
#ENDIF
!
!.......................................................................
MODULE RANDOM_MODULE
  INTEGER(8),SAVE        :: SEED=1
  INTEGER(8),PARAMETER   :: RANFAC1=22222221
  INTEGER(8),PARAMETER   :: RANFAC2=2**24
! CHOICE EXPLICIT CPPVARIABLE_XLF RNG
! INTEGER(8),SAVE        :: SEED=1_8
! INTEGER(8),PARAMETER   :: RANFAC1=44485709377909_8
! INTEGER(8),PARAMETER   :: RANFAC2=2_8**48
END MODULE RANDOM_MODULE
!
!     ..................................................................
      SUBROUTINE LIB$RANDOM(X)
!     ****************************************************************** 
!     **  RETURNS A RANDOM NUMBER                                     **
!     ******************************************************************
      USE RANDOM_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: X
!     ******************************************************************
      SEED=MODULO(RANFAC1*SEED,RANFAC2)
      X=REAL(SEED,KIND=8)/REAL(RANFAC2,KIND=8)
!     == CHOICE EXPLICIT CPPVARIABLE_XLF RNG
!     CALL RANDOM_NUMBER(Y)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$RANDOMSEED
!     ****************************************************************** 
!     **  CHOOSE A SEED  FOR THE NADOM NUMBER GENERATOR               **
!     ******************************************************************
      USE RANDOM_MODULE
      IMPLICIT NONE
!     ******************************************************************
      SEED=1
!     == CHOICE EXPLICIT CPPVARIABLE_XLF RNG
!     SEED=1_8
!     == CHOICE ORIGINAL CPPVARIABLE_XLF RNG
!     CALL RANDOM_SEED(GENERATOR=2)
!     CALL RANDOM_SEED
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$INVERTR8(N,A,AINV)
!     ******************************************************************
!     **                                                              **
!     **  INVERTS THE REAL, SQUARE MATRIX A                           **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: DGEICD                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: A(N,N)
      REAL(8)   ,INTENT(OUT):: AINV(N,N)
      INTEGER(4)            :: NAUX
      REAL(8)               :: AUX(100*N)
      INTEGER(4)            :: I,J
#IF DEFINED(CPPVAR_BLAS_ESSL)
      REAL(8)               :: RCOND
      REAL(8)               :: DET(2)
#ELSE 
      INTEGER(4)            :: IPIV(N)
      INTEGER(4)            :: INFO
#ENDIF
!     ******************************************************************
      NAUX=100*N
      AINV(1:N,1:N)=A(1:N,1:N)

#IF DEFINED(CPPVAR_BLAS_ESSL)
      CALL DGEICD(AINV,N,N,0,RCOND,DET,AUX,NAUX) !ESSL
#ELSE 
      CALL DGETRF(N,N,AINV,N,IPIV,INFO) !LAPACK
      CALL DGETRI(N,AINV,N,IPIV,AUX,NAUX,INFO) !LAPACK
#ENDIF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$DIAGR8(N,H,E,U)
!     ******************************************************************
!     **                                                              **
!     **  DIAGONALIZES THE REAL, SQUARE MATRIX H AFTER SYMMETRIZATION **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                   **
!     **                                                              **
!     **         U(K,I)*H(K,L)*U(L,J)=DELTA(I,J)*E(I)                 **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL DSPEV  MATRIX DIAGONALIZATION P727                   **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
#IF DEFINED(CPPVAR_ESSL)
      INTERFACE
        SUBROUTINE EINFO(ICODE,INF1,INF2) !ESSL ERROR HANDLING ROUTINE
        INTEGER                       :: ICODE
        INTEGER ,INTENT(OUT),OPTIONAL :: INF1
        INTEGER ,INTENT(OUT),OPTIONAL :: INF2
        END SUBROUTINE EINFO
      END INTERFACE
#ENDIF
      LOGICAL(4) ,PARAMETER :: TESSLERR=.FALSE.
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      REAL(8)   ,INTENT(OUT):: U(N,N)
      REAL(8)               :: WORK1((N*(N+1))/2)
#IF DEFINED(CPPVAR_BLAS_ESSL)
      REAL(8)               :: WORK2(2*N)
#ELSE
      REAL(8)               :: WORK2(3*N)
#ENDIF
      INTEGER(4)            :: K,I,J
      CHARACTER(8)          :: SAV2101
      INTEGER(4)            :: I1,I2
      INTEGER(4)            :: INFO
      INTEGER               :: INF1
      INTEGER               :: INF2
!     ******************************************************************
!
!     ==================================================================
!     ====  STORE IN LOWER PACKED STORAGE MODE FOR DIAGONALIZATION    ==
!     ==================================================================
      K=0
      DO J=1,N
        DO I=J,N
          K=K+1
          WORK1(K)=0.5D0*(H(I,J)+H(J,I))
        ENDDO
      ENDDO
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
#IF DEFINED(CPPVAR_BLAS_ESSL)
      IF(TESSLERR) THEN
        CALL EINFO(0,INF1,INF2)
        CALL ERRSAV(2101,SAV2101)
        CALL ERRSET(2101,255,0,0,0,2101)
!       CALL DSPEV(1,WORK1,E,U,N,N,WORK2,2*N,400) !->ESSL
        CALL DSPEV(1,WORK1,E,U,N,N,WORK2,2*N) !->ESSL
      ELSE
        CALL DSPEV(1,WORK1,E,U,N,N,WORK2,2*N) !->ESSL
      END IF
#ELSE
      CALL DSPEV('V','L',N,WORK1,E,U,N,WORK2,INFO) !->LAPACK
      IF(INFO.NE.0) THEN
        CALL ERROR$MSG('DIAGONALIZATION NOT CONVERGED')
        CALL ERROR$STOP('DIAG')
      END IF
#ENDIF
      RETURN
!
!     ==================================================================
!     == ESSL ERROR HANDLING                                          ==
!     ==================================================================
#IF DEFINED(CPPVAR_BLAS_ESSL)
 400  CONTINUE
      CALL ERROR$MSG('DIAGONALIZATION NOT CONVERGED')
      IF(TESSLERR) THEN
        CALL EINFO(2101,I1,I2)
        CALL ERROR$I4VAL('EIGENVALUE FAILED TO CONVERGE',I1)
        CALL ERROR$I4VAL('NUMBER OF ITERATIONS',I2)
      END IF
      DO I=1,N
        PRINT*,'H ',I,H(1:3,I)
      ENDDO
      CALL ERROR$STOP('LIB$DIAGR8')
      STOP
#ENDIF
      END
!
!     ..................................................................
      SUBROUTINE LIB$DIAGC8(N,H,E,U)
!     ******************************************************************
!     **                                                              **
!     **  DIAGONALIZES THE HERMITEAN, SQUARE MATRIX H                 **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                   **
!     **                                                              **
!     **      CONJG(U(K,I))*H(K,L)*U(L,J)=DELTA(I,J)*E(I)             **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: ZHPEV :  MATRIX DIAGONALIZATION  P727               **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTERFACE
        SUBROUTINE EINFO(ICODE,INF1,INF2) !ESSL ERROR HANDLING ROUTINE
        INTEGER                       :: ICODE
        INTEGER ,INTENT(OUT),OPTIONAL :: INF1
        INTEGER ,INTENT(OUT),OPTIONAL :: INF2
        END SUBROUTINE EINFO
      END INTERFACE
      LOGICAL(4) ,PARAMETER :: TESSLERR=.FALSE.
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      COMPLEX(8),INTENT(OUT):: U(N,N)
      COMPLEX(8)            :: WORK1((N*(N+1))/2)
#IF DEFINED(CPPVAR_BLAS_ESSL)
      REAL(8)               :: RWORK(4*N)
#ELSE
      COMPLEX(8)            :: CWORK(2*N)
      REAL(8)               :: RWORK(3*N)
#ENDIF
      INTEGER(4)            :: K,I,J
      CHARACTER(8)          :: SAV2101
      INTEGER(4)            :: I1,I2
      LOGICAL(4),PARAMETER  :: TTEST=.FALSE.
      COMPLEX(8)            :: CSVAR
      LOGICAL(4)            :: TCHK
      COMPLEX(8),ALLOCATABLE:: WORK3(:,:)
      INTEGER(4)            :: INFO
!CLEMENS
      INTEGER               :: INF1
      INTEGER               :: INF2
!     ******************************************************************

!     ==================================================================
!     ====  ESSL ERROR HANDLING                                       ==
!     ==================================================================
!
!     ==================================================================
!     ====  STORE IN LOWER PACKED STORAGE MODE FOR DIAGONALIZATION    ==
!     ==================================================================
      K=0
      DO J=1,N
        DO I=J,N
          K=K+1
          WORK1(K)=0.5D0*(H(I,J)+CONJG(H(J,I)))
        ENDDO
      ENDDO
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
#IF DEFINED(CPPVAR_BLAS_ESSL)
      IF(TESSLERR) THEN
        CALL EINFO(0,INF1,INF2)
        CALL ERRSAV(2101,SAV2101)
        CALL ERRSET(2101,255,0,0,0,2101)
      END IF
      CALL ZHPEV(1,WORK1,E,U,N,N,RWORK,4*N)  !ESSL
#ELSE
      CALL ZHPEV('V','L',N,WORK1,E,U,N,CWORK,RWORK,INFO) !LAPACK
      IF(INFO.NE.0) THEN
        CALL ERROR$MSG('DIAGONALIZATION NOT CONVERGED')
        CALL ERROR$STOP('DIAG')
      END IF
#ENDIF
!
!     ==================================================================
!     ====  OPTIONAL TEST                                             ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(WORK3(N,N)) 
        TCHK=.TRUE.
        DO I=1,N
          DO J=1,N
            CSVAR=(0.D0,0.D0)
            DO K=1,N
              CSVAR=CSVAR+U(I,K)*E(K)*CONJG(U(J,K))
            ENDDO
            IF(ABS(CSVAR-H(I,J)).GT.1.D-8) THEN
              WRITE(*,FMT='(2I5,5F10.5)')I,J,CSVAR,H(I,J),ABS(CSVAR-H(I,J))
              TCHK=.FALSE.
            END IF
          ENDDO
        ENDDO
        DEALLOCATE(WORK3)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$STOP('LIB$DIAGC8')
        END IF
      END IF
      RETURN
!
!     ==================================================================
!     == ESSL ERROR HANDLING                                          ==
!     ==================================================================
#IF DEFINED(CPPVAR_BLAS_ESSL)
 400  CONTINUE
      CALL ERROR$MSG('DIAGONALIZATION NOT CONVERGED')
      IF(TESSLERR) THEN
        CALL EINFO(2101,I1,I2)
        CALL ERROR$I4VAL('EIGENVALUE FAILED TO CONVERGE',I1)
        CALL ERROR$I4VAL('NUMBER OF ITERATIONS',I2)
      END IF
      DO I=1,N
        PRINT*,'H ',I,H(1:3,I)
      ENDDO
      CALL ERROR$STOP('LIB$DIAGC8')
      STOP
#ENDIF
      END
!
!     ......................................................................
      subroutine lib$generaleigenvaluer8(n,a,b,e,vec)
!     == solves the generalized, real non-symmetric eigenvalue problem   **
!     ** not tested!!!!!!                                                **
      implicit none
      integer(4),intent(in) :: n
      real(8)   ,intent(in) :: a(n,n)
      real(8)   ,intent(in) :: b(n,n)
      real(8)   ,intent(out):: e(n)
      real(8)   ,intent(out) :: vec(n,n)
      integer                :: lwork     
      real(8)                :: work(16*n)
      real(8)                :: alphar(n)
      real(8)                :: alphai(n)
      real(8)                :: beta(n)
      integer(4)             :: info
      real(8)                :: a1(n,n)
      real(8)                :: b1(n,n)
!     *********************************************************************
      lwork=16*n
      a1=a   !dggev overwrites input!!
      b1=b
      call dggev('N','V',n,a1,n,b1,n,alphar,alphai,beta &
     &           ,vec,n,vec,n,work,lwork,info)
      if(abs(maxval(alphai(:))).gt.0.d0)print*,'maxval(alphai(:))',maxval(alphai(:))
      e(:)=alphar(:)/beta(:)
      return
      end
!
!     ..................................................................
      SUBROUTINE LIB$FFTC8(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     **  USE FFTW AS STANDARD                                        **
!     **  USE FFTESSL IF ESSL IS INSTALLED                            **
!     **  USE FFTPACK AS BACKUP IF C-ROUTINES CANNOT BE LINKED OR     **
!     **      FFTW IS NOT AVAILABLE                                   **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
!     *******************************************************************
#IF DEFINED(CPPVAR_FFT_ESSL)
      CALL LIB_FFTESSL(DIR,LEN,NFFT,X,Y)                  
#ELIF DEFINED(CPPVAR_FFT_FFTW)
      CALL LIB_FFTW(DIR,LEN,NFFT,X,Y)
#ELIF DEFINED(CPPVAR_FFT_PACK)
      CALL LIB_FFTPACK(DIR,LEN,NFFT,X,Y)
#ELIF DEFINED(CPPVAR_FFT_CXML)
      CALL LIB_FFTCXML(DIR,LEN,NFFT,X,Y)
#ELSE
      CALL ERROR$MSG('NO FFT PACKAGE SELECTED DURING COMPILATION')
      CALL ERROR$STOP('LIB$FFTC8')
#ENDIF

      RETURN
      END
!
#IF DEFINED(CPPVAR_FFT_CXML)
!
!  ATTENTION !!!
!  DO NOT USE !!
!  CLEMENS FOERST
!     ..................................................................
      SUBROUTINE LIB_FFTCXML(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INCLUDE 'DXMLDEF.FOR'
      RECORD /DXML_Z_FFT_STRUCTURE/  :: FFT_STRUCT
      CHARACTER(1)            :: DIR1
      INTEGER(4)              :: STATUS 
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      CHARACTER(4),SAVE       :: DIRSAVE=''
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      REAL(8)     ,SAVE       :: SCALE
      INTEGER(4)              :: I
!     ******************************************************************
!
!      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE.OR.NFFT.NE.NFFTSAVE) THEN

         IF(DIR.EQ.'GTOR') THEN
            DIR1='B'
         ELSE IF(DIR.EQ.'RTOG') THEN
            DIR1='F'
         ELSE 
            CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
            CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
            CALL ERROR$CHVAL('DIR',TRIM(DIR))
            CALL ERROR$STOP('FFT')
         END IF

!         DIRSAVE=DIR
!         LENSAVE=LEN
!         NFFTSAVE=NFFT
!      END IF
      STATUS=ZFFT_INIT(LEN,FFT_STRUCT,.TRUE.)
      DO I=1,NFFT
         STATUS=ZFFT_APPLY('C','C',DIR1,X(:,I),Y(:,I),FFT_STRUCT,1)
      END DO
      STATUS=ZFFT_EXIT_GRP(FFT_STRUCT)
      RETURN
      END
#ENDIF
!
#IF DEFINED(CPPVAR_FFT_FFTW)
!     ..................................................................
      SUBROUTINE LIB_FFTW(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INCLUDE 'FFTW_F77.I'
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      CHARACTER(4),SAVE       :: DIRSAVE=''
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      INTEGER(4)  ,SAVE       :: ISIGN
      REAL(8)     ,SAVE       :: SCALE
      COMPLEX(8)              :: XDUMMY(LEN,NFFT)
      INTEGER(4)              :: I
      INTEGER(4),SAVE         :: NP=0
      INTEGER(4),PARAMETER    :: NPX=10 ! #(DIFFERENT FFT PLANS)
#IF DEFINED(CPPVAR_64BIT)
      INTEGER(8),SAVE         :: PLANS(NPX,2),PLAN=-1
#ELSE
      INTEGER(4),SAVE         :: PLANS(NPX,2),PLAN=-1
#ENDIF
      LOGICAL                 :: DEF
!     ******************************************************************
      XDUMMY=X
!
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE) THEN
        IF (DIR.EQ.'GTOR') THEN
          ISIGN=1
        ELSE IF (DIR.EQ.'RTOG') THEN
          ISIGN=-1
        ELSE
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('1D-FFTW')
        END IF
!
!       == FIND PLAN IN THE LIST
        DEF=.FALSE.
        DO I=1,NP
          IF((LEN*ISIGN).EQ.PLANS(I,1)) THEN
            DEF=.TRUE.
            PLAN=PLANS(I,2)
            EXIT
          END IF
        END DO
!
!       == CREATE NEW PLAN IF NOT IN THE LIST ==========================
        IF(.NOT.DEF) THEN
          WRITE(*,*) 'FFTW CREATE PLAN FOR: ', ISIGN,LEN,NP
          NP=NP+1
          IF(NP.GE.NPX) NP=NPX ! ALLOW ONLY NPX PLANS
          CALL FFTW_F77_CREATE_PLAN(PLANS(NP,2),LEN,ISIGN,FFTW_MEASURE)
          PLANS(NP,1)=ISIGN*LEN
          PLAN=PLANS(NP,2)
        END IF
        LENSAVE=LEN
        DIRSAVE=DIR
        SCALE=1.D0/REAL(LEN,KIND=8)
      END IF
!
!     ==================================================================
!     ==  NOW PERFORM FFT                                             ==
!     ==================================================================
      CALL FFTW_F77(PLAN,NFFT,XDUMMY(1,1),1,LEN,Y(1,1),1,LEN)
      IF (DIR.EQ.'RTOG') THEN
        Y(:,:)=Y(:,:)*SCALE
      END IF
      RETURN
      END
#ENDIF
!
!DCFT:  1-D FFT P765
#IF DEFINED(CPPVAR_FFT_ESSL)
!     ..................................................................
      SUBROUTINE LIB_FFTESSL(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      INTEGER(4)  ,PARAMETER  :: NAUX1=20000
      INTEGER(4)  ,PARAMETER  :: NAUX2=20000
      REAL(8)     ,SAVE       :: AUX2(NAUX2)
      REAL(8)     ,SAVE       :: AUX1(NAUX1)
      CHARACTER(4),SAVE       :: DIRSAVE=''
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      INTEGER(4)  ,SAVE       :: ISIGN
      REAL(8)     ,SAVE       :: SCALE
      INTEGER(4)              :: IFFT,I
!     ******************************************************************
!  
!     ==================================================================
!     == INITIALIZATION PHASE                                         ==
!     ==================================================================
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE.OR.NFFT.NE.NFFTSAVE) THEN
        IF(LEN.GT.8192) THEN
          CALL ERROR$MSG('FFT TOO LONG')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('FFT')
        END IF
        IF(DIR.EQ.'GTOR') THEN
          ISIGN=-1
          SCALE=1.D0
        ELSE IF(DIR.EQ.'RTOG') THEN
          ISIGN=1
          SCALE=1.D0/REAL(LEN,KIND=8)
        ELSE 
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('FFT')
        END IF
        CALL DCFT(1,X,1,LEN,Y,1,LEN,LEN,NFFT,ISIGN,SCALE,AUX1,NAUX1,AUX2,NAUX2)
        DIRSAVE=DIR
        LENSAVE=LEN
        NFFTSAVE=NFFT
      END IF
!  
!     ==================================================================
!     == FOURIER TRANSFORM                                            ==
!     ==================================================================
      CALL DCFT(0,X,1,LEN,Y,1,LEN,LEN,NFFT,ISIGN,SCALE,AUX1,NAUX1,AUX2,NAUX2)
      RETURN
      END
#ENDIF
!
!DCFT:  1-D FFT P765
#IF DEFINED(CPPVAR_FFT_PACK)
!     ..................................................................
      SUBROUTINE LIB_FFTPACK(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      INTEGER(4)  ,PARAMETER  :: NAUX2=2000
      REAL(8)     ,SAVE       :: AUX2(NAUX2)
      INTEGER(4)  ,SAVE       :: IAUX(30)
      REAL(8)                 :: SEQUENCE(2*LEN)
      CHARACTER(4),SAVE       :: DIRSAVE=''
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      INTEGER(4)  ,SAVE       :: ISIGN
      REAL(8)     ,SAVE       :: SCALE
      INTEGER(4)              :: IFFT,I
!     ******************************************************************
!  
!     ==================================================================
!     == INITIALIZATION PHASE                                         ==
!     ==================================================================
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE) THEN
        IF(DIR.EQ.'GTOR') THEN
          ISIGN=-1
          SCALE=1.D0
        ELSE IF(DIR.EQ.'RTOG') THEN
          ISIGN=1
          SCALE=1.D0/REAL(LEN,KIND=8)
        ELSE 
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('LIB_FFTPACK')
        END IF
        DIRSAVE=DIR
        LENSAVE=LEN
        IF(LEN.EQ.1) RETURN
        IF(NAUX2.LT.2*LEN) THEN
          CALL ERROR$MSG('AUXILIARY ARRAY TOO SMALL: INCREASE NAUX2')
          CALL ERROR$I4VAL('NAUX2',NAUX2)
          CALL ERROR$I4VAL('2*LEN',2*LEN)
          CALL ERROR$STOP('LIB_FFTPACK')
        END IF
        CALL CFFTI(LEN,AUX2,IAUX)
      END IF
!  
!     ==================================================================
!     == FOURIER TRANSFORM                                            ==
!     ==================================================================
      IF(LEN.EQ.1) RETURN
      IF(ISIGN.EQ.1) THEN
        DO IFFT=1,NFFT
          DO I=1,LEN
            SEQUENCE(2*I-1)= REAL(X(I,IFFT),KIND=8)
            SEQUENCE(2*I  )=AIMAG(X(I,IFFT))
          ENDDO
          CALL CFFTF(LEN,SEQUENCE,AUX2,IAUX)
          DO I=1,LEN
            Y(I,IFFT)=CMPLX(SEQUENCE(2*I-1),SEQUENCE(2*I),KIND=8)*SCALE
          ENDDO
        ENDDO
      ELSE
        DO IFFT=1,NFFT
          DO I=1,LEN
            SEQUENCE(2*I-1)= REAL(X(I,IFFT),KIND=8)
            SEQUENCE(2*I  )=AIMAG(X(I,IFFT))
          ENDDO
          CALL CFFTB(LEN,SEQUENCE,AUX2,IAUX)
          DO I=1,LEN
            Y(I,IFFT)=CMPLX(SEQUENCE(2*I-1),SEQUENCE(2*I),KIND=8)*SCALE
          ENDDO
        ENDDO
      END IF
      RETURN
      END
#ENDIF
!
!     ...................................................FESSL..........
      SUBROUTINE LIB$FFTADJUSTGRD(NR)
!     ******************************************************************
!     **  THIS ROUTINE RETURNS THE ALLOWED FOURIER TRANSFORM LENGTH   ** 
!     **  THAT IS EQUAL OR LARGER THAN THE LENGTH SUPPLIED, BUT       ** 
!     **  BUT OTHERWISE AS SMALL AS POSSIBLE.                         **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT):: NR
      INTEGER(4),PARAMETER  :: MAXI=300
      LOGICAL(4),SAVE       :: TINIT=.TRUE.
      INTEGER(4),SAVE       :: COUNT
      INTEGER(4),SAVE       :: IFR(MAXI)
      INTEGER(4)            :: H,I,J,K,M
      INTEGER(4)            :: ISVAR
      REAL(8)               :: SVAR
!     ******************************************************************
      IF (TINIT) THEN
        TINIT=.FALSE.
#IF DEFINED(CPPVAR_FFT_ESSL)
!       == ALLOWED LENGTH FOR THE ESSL FFT =============================
        COUNT=0
        OUTER: DO H=1,25
          DO I=0,2
            DO J=0,1
              DO K=0,1
                DO M=0,1
                  IF(COUNT.GE.MAXI) EXIT OUTER
                  SVAR = 2.D0**H * 3.D0**I * 5.D0**J * 7.D0**K *11.D0**M
                  IF(SVAR.GT.37748736.D0) CYCLE
                  COUNT=COUNT+1
                  IFR(COUNT)=NINT(SVAR)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO OUTER
#ELIF DEFINED(CPPVAR_FFT_FFTW)
        COUNT=0
        OUTER: DO H=1,25
          DO I=0,2
            DO J=0,1
              DO K=0,1
                DO M=0,1
                  IF(COUNT.GE.MAXI) EXIT OUTER
                  SVAR = 2.D0**H * 3.D0**I * 5.D0**J * 7.D0**K *11.D0**M
                  IF(SVAR.GT.37748736.D0) CYCLE
                  COUNT=COUNT+1
                  IFR(COUNT)=NINT(SVAR)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO OUTER
#ELSE
!       == THE GENERAL LENGTH HAS BEEN REMOVED BECAUSE PASSF AND PASSB
!       == DO NOT WORK
        COUNT=0
        OUTER: DO H=1,25
          DO I=0,4
            DO J=0,2
              IF(COUNT.GE.MAXI) EXIT OUTER
              SVAR = 2.D0**H * 3.D0**I * 5.D0**J 
              IF(SVAR.GT.37748736.D0) CYCLE
              COUNT=COUNT+1
              IFR(COUNT)=NINT(SVAR)
            ENDDO
          ENDDO
        ENDDO OUTER
#ENDIF
        DO I=1,COUNT
          DO J=I+1,COUNT
            IF(IFR(I).GT.IFR(J)) THEN
              ISVAR=IFR(I)
              IFR(I)=IFR(J)
              IFR(J)=ISVAR
            END IF
          ENDDO
        ENDDO
      ENDIF
!
      DO I=2,COUNT
        IF(IFR(I).GE.NR) THEN
          NR=IFR(I)
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('REQUESTED GRIDPOINT OUT OF RANGE')
      CALL ERROR$I4VAL('NR',NR)
      CALL ERROR$STOP('PLANEWAVE$ADJUSTFFTGRD')
      STOP
      END
!
!     ..................................................................
      SUBROUTINE LIB$3DFFTC8(DIR,N1,N2,N3,X,Y)
!     ******************************************************************
!     **  3-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **    USES THE 3D FFTW ROUTINES                                 **
!     **                                        CLEMENS FOERST, 2001  **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
!     ******************************************************************
#IF DEFINED(CPPVAR_FFT_ESSL)
      CALL LIB_3DFFT_ESSL(DIR,N1,N2,N3,X,Y)
#ELIF DEFINED(CPPVAR_FFT_FFTW)
      CALL LIB_3DFFTW(DIR,N1,N2,N3,X,Y)
#ELIF DEFINED(CPPVAR_FFT_PACK)
      CALL LIB_3DFFTPACK(DIR,N1,N2,N3,X,Y)
#ELSE
      CALL ERROR$MSG('NO FFT PACKAGE SELECTED DURING COMPILATION')
      CALL ERROR$STOP('LIB$3DFFTC8')
#ENDIF
      RETURN
      END
!
#IF DEFINED(CPPVAR_FFT_FFTW)
!     ..................................................................
      SUBROUTINE LIB_3DFFTW(DIR,N1,N2,N3,X,Y)
!     ******************************************************************
!     **  3-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **    USES THE 3D FFTW ROUTINES                                 **
!     **                                        CLEMENS FOERST, 2001  **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: DIM(3)
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      INTEGER(4)              :: PLAN
      INTEGER(4)  ,SAVE       :: NP=0
      INTEGER(4),PARAMETER    :: NPX=10
      INTEGER(4)  ,SAVE       :: PLANS(NPX,4)
      REAL(8)     ,SAVE       :: SCALE
      INTEGER(4)  ,SAVE       :: DIMSAVE(3)=0
      CHARACTER(4),SAVE       :: DIRSAVE=''
      LOGICAL                 :: DEF
      INTEGER(4)              :: I
      INTEGER(4)  ,SAVE       :: ISIGN
      INCLUDE 'FFTW_F77.I'
!     ******************************************************************
      DIM(1)=N1
      DIM(2)=N2
      DIM(3)=N3
      IF(DIM(1).NE.DIMSAVE(1).OR.DIM(2).NE.DIMSAVE(2).OR. &
     &  DIM(3).NE.DIMSAVE(3).OR.DIR.NE.DIRSAVE) THEN
        IF (DIR.EQ.'GTOR') THEN
          ISIGN=1
        ELSE IF (DIR.EQ.'RTOG') THEN
          ISIGN=-1
        ELSE
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('3D-FFTW')
        END IF
!
!       == FIND PLAN IN THE LIST
        DEF=.FALSE.
        DO I=1,NP
          IF((DIM(1)*ISIGN).EQ.PLANS(I,1).AND.(DIM(2)*ISIGN).EQ.PLANS(I,2)&
     &                     .AND.(DIM(3)*ISIGN).EQ.PLANS(I,3)) THEN
            DEF=.TRUE.
            PLAN=PLANS(I,4)
            EXIT
          END IF
        END DO
!
!       == CREATE NEW PLAN IF NOT IN THE LIST ==========================
        IF(.NOT.DEF) THEN
          WRITE(*,*) '3D-FFTW CREATE PLAN FOR: ', ISIGN,DIM
          NP=NP+1
          IF(NP.GE.NPX) NP=NPX ! ALLOW ONLY NPX PLANS
          CALL FFTWND_F77_CREATE_PLAN(PLAN,3,DIM,ISIGN,FFTW_ESTIMATE)
          PLANS(NP,1:3)=ISIGN*DIM
          PLANS(NP,4)=PLAN
        END IF
        DIMSAVE=DIM
        DIRSAVE=DIR
        IF(DIR.EQ.'RTOG') SCALE=1.D0/REAL(N1*N2*N3,KIND=8)
      END IF
!
!     ==================================================================
!     ==  NOW PERFORM FFT                                             ==
!     ==================================================================
      CALL FFTWND_F77_ONE(PLAN,X,Y)
      IF (DIR.EQ.'RTOG') Y=Y*SCALE
      RETURN
      END
!
#ELIF DEFINED(CPPVAR_FFT_ESSL)
!     ..................................................................
      SUBROUTINE LIB_3DFFT_ESSL(DIR,N1,N2,N3,X,Y)
!     ******************************************************************
!     **  3-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **    USES THE 3D FFT ROUTINE OF ESSL DCFT3                     **
!     **                                        PETER BLOECHL, 2001   **
!     ******************************************************************
!     **    NOT TESTED                                                **
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: DIM(3)
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      INTEGER(4)              :: NAUX       
      REAL(8)   ,ALLOCATABLE  :: AUX(:)
      INTEGER(4)              :: ISIGN
      REAL(8)                 :: SCALE
      INTEGER(4)              :: S,PSI,LAMBDA
!     ******************************************************************
      NAUX=60000
      IF(N1.GT.2048) NAUX=NAUX+4.56*N1
      IF(N3.LT.252) THEN
        IF(N2.GE.252) THEN
          S=MIN(64,N1)
          LAMBDA=(2*N2+256)*(S+4.56)
          NAUX=NAUX+LAMBDA
        END IF
      ELSE
        IF(N2.GE.252) THEN
          S=MIN(64,N1*N2)
          PSI=(2*N3+256)*(S+4.56)
          NAUX=NAUX+PSI
        ELSE
          S=MIN(64,N1*N2)
          PSI=(2*N3+256)*(S+4.56)
          S=MIN(64,N1)
          LAMBDA=(2*N2+256)*(S+4.56)
          NAUX=NAUX+MAX(PSI,LAMBDA)
        END IF
      END IF
      IF (DIR.EQ.'GTOR') THEN
        ISIGN=1
        SCALE=1.D0
      ELSE IF (DIR.EQ.'RTOG') THEN
        ISIGN=-1
        SCALE=1.D0/REAL(N1*N2*N3,KIND=8)
      ELSE
        CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
        CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
        CALL ERROR$CHVAL('DIR',TRIM(DIR))
        CALL ERROR$STOP('3D-FFTW')
      END IF
      ALLOCATE(AUX(NAUX))
      CALL DCFT3(X,N1,N1*N2,Y,N1,N1*N2,N1,N2,N3,ISIGN,SCALE,AUX,NAUX)
      DEALLOCATE(AUX)
      RETURN
      END
!
#ELIF DEFINED(CPPVAR_FFT_PACK)
!     ..................................................................
      SUBROUTINE LIB_3DFFTPACK(DIR,N1,N2,N3,X,Y)
!     ******************************************************************
!     **  3-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **    USES THE 3D FFT ROUTINE OF ESSL DCFT3                     **
!     **                                        PETER BLOECHL, 2001   **
!     ******************************************************************
!     **    NOT TESTED                                                **
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: DIM(3)
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      COMPLEX(8)              :: WORK1(N1*N2*N3),WORK2(N1*N2*N3)
      INTEGER(4)              :: I,J,K,IND
!     ******************************************************************
      CALL LIB_FFTPACK(DIR,N1,N2*N3,X,WORK2)
      IND=0
      DO K=1,N3
        DO J=1,N2
          DO I=1,N1
            IND=IND+1
            WORK1(J+N2*(K-1+N1*(I-1)))=WORK2(IND)
          ENDDO
        ENDDO
      ENDDO
      CALL LIB_FFTPACK(DIR,N2,N1*N3,WORK1,WORK2)
      IND=0
      DO I=1,N1
        DO K=1,N3
          DO J=1,N2
            IND=IND+1
            WORK1(K+N3*(I-1+N1*(J-1)))=WORK2(IND)
          ENDDO
        ENDDO
      ENDDO
      CALL LIB_FFTPACK(DIR,N3,N1*N2,WORK1,WORK2)
      IND=0
      DO J=1,N2
        DO I=1,N1
          DO K=1,N3
            IND=IND+1
            Y(I,J,K)=WORK2(IND)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
#ENDIF
!
!DGEMUL: MATRIX MULTIPLICATION P441
!     ..................................................................
      SUBROUTINE LIB$MATMULR8(N,M,L,A,B,C)
!     ******************************************************************
!     **  MATRIX MULTPLICATION   A*B=C                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(IN) :: B(M,L)
      REAL(8)   ,INTENT(OUT):: C(N,L)
      CHARACTER(8),PARAMETER:: LIB='ESSL'
      REAL(8)     ,PARAMETER:: ONE=1.D0
      REAL(8)     ,PARAMETER:: ZERO=0.D0
!     ******************************************************************
#IF DEFINED(CPPVAR_BLAS_ESSL)
      CALL DGEMUL(A,N,'N',B,M,'N',C,N,N,M,L)
!#ELIF DEFINED(CPPVAR_BLAS_ATLAS)
!      CALL DGEMM('N','N',N,L,M,ONE,A,N,B,M,ZERO,C,N)
#ELSE
      CALL DGEMM('N','N',N,L,M,ONE,A,N,B,M,ZERO,C,N)
#ENDIF
      RETURN
      END
!
!ZGEMUL: MATRIX MULTIPLICATION P441
!     . ................................................................
      SUBROUTINE LIB$MATMULC8(N,M,L,A,B,C)
!     ******************************************************************
!     **  MATRIX MULTPLICATION   A*B=C                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: L
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(IN) :: B(M,L)
      COMPLEX(8),INTENT(OUT):: C(N,L)
      CHARACTER(8),PARAMETER:: LIB='ESSL'
      COMPLEX(8)  ,PARAMETER:: ONE=(1.D0,0.D0)
      COMPLEX(8)  ,PARAMETER:: ZERO=(0.D0,0.D0)
!     ******************************************************************
#IF DEFINED(CPPVAR_BLAS_ESSL)
      CALL ZGEMUL(A,N,'N',B,M,'N',C,N,N,M,L)
!#ELIF DEFINED(CPPVAR_BLAS_ATLAS)
!      CALL ZGEMM('N','N',N,L,M,ONE,A,N,B,M,ZERO,C,N)
#ELSE
      C(:,:)=(0.D0,0.D0)
      CALL ZGEMM('N','N',N,L,M,ONE,A,N,B,M,ZERO,C,N)
#ENDIF
      RETURN
      END
!
!ZGEMM(TRANSA,TRANSB,L,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC):
!C=BETA*C+ALPHA*A*B P456
!     ..................................................................
      SUBROUTINE LIB$ADDPRODUCTC8(TID,N,M,L,A,B,C)
!     ******************************************************************
!     **  ADD THE PRODUCT OF TWO MATRICES                             **
!     **  C=C+MATMUL(A,B)                                             **
!     **  FOR TID=.TRUE., A AND C MAY USE IDENTICAL STORAGE           **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)   :: TID
      INTEGER(4),INTENT(IN)   :: N
      INTEGER(4),INTENT(IN)   :: M
      INTEGER(4),INTENT(IN)   :: L
      COMPLEX(8),INTENT(IN)   :: A(N,M)
      COMPLEX(8),INTENT(IN)   :: B(M,L)
      COMPLEX(8),INTENT(INOUT):: C(N,L)
      COMPLEX(8)              :: SUM
      COMPLEX(8),ALLOCATABLE  :: WORK(:,:)
      INTEGER(4)              :: I,J,K
!     ******************************************************************
      IF(TID) THEN
        ALLOCATE(WORK(N,M))
        WORK=A
!       ==  C=C+MATMUL(WORK,B) 
        CALL ZGEMM('N','N',N,L,M,(1.D0,0.D0),WORK,N,B,M,(1.D0,0.D0),C,N)
        DEALLOCATE(WORK)
      ELSE
!       ==  C=C+MATMUL(A,B) 
        CALL ZGEMM('N','N',N,L,M,(1.D0,0.D0),A,N,B,M,(1.D0,0.D0),C,N)
      END IF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$DYADSUMR8(LEN1,LEN2,N,PSI1,PSI2,OPERATOR)
!     ******************************************************************
!     **   OPERATOR = SUM_I |PSI1(I)><PSI2(I)|                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LEN1
      INTEGER(4),INTENT(IN) :: LEN2
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: PSI1(LEN1,N)
      REAL(8)   ,INTENT(IN) :: PSI2(LEN2,N)
      REAL(8)   ,INTENT(OUT):: OPERATOR(LEN1,LEN2)
      INTEGER(4)            :: I,J,K
      REAL(8)               :: SUM
!     ******************************************************************
#IF DEFINED(CPPVAR_BLAS_ESSL)
      CALL DGEMUL(PSI1,LEN1,'N',PSI2,LEN2,'T',OPERATOR,LEN1,LEN1,N,LEN2)
#ELSE 
!     == OPERATOR=MATMUL(PSI1,TRANSPOSE(PSI2))
      CALL DGEMM('N','T',LEN1,LEN2,N,1.D0,PSI1,LEN1,PSI2,LEN2,0.D0 &
     &          ,OPERATOR,LEN1)
#ENDIF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$DYADSUMC8(LEN1,LEN2,N,PSI1,PSI2,OPERATOR)
!     ******************************************************************
!     **   OPERATOR = SUM_I |PSI1(I)><PSI2(I)|                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LEN1
      INTEGER(4),INTENT(IN) :: LEN2
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: PSI1(LEN1,N)
      COMPLEX(8),INTENT(IN) :: PSI2(LEN2,N)
      COMPLEX(8),INTENT(OUT):: OPERATOR(LEN1,LEN2)
      INTEGER(4)            :: I,J,K
      COMPLEX(8)            :: SUM
!     ******************************************************************
#IF DEFINED(CPPVAR_BLAS_ESSL)
      CALL ZGEMUL(PSI1,LEN1,'N',PSI2,LEN2,'C',OPERATOR,LEN1,LEN1,N,LEN2)
#ELSE 
!     == OPERATOR=MATMUL(PSI1,TRANSPOSE(PSI2))
      OPERATOR(:,:)=(0.D0,0.D0)
      CALL ZGEMM('N','C',LEN1,LEN2,N,(1.D0,0.D0) &
     &          ,PSI1(:,:),LEN1,PSI2(:,:),LEN2,(0.D0,0.D0),OPERATOR,LEN1)
#ENDIF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$SCALARPRODUCTR8(TID,LEN,N1,PSI1,N2,PSI2,OVERLAP)
!     ******************************************************************
!     ** PERFORMS THE SCALAR PRODUCTS OF TWO ARRAYS OF COMPLEX        **
!     ** STATE VECTORS                                                **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TID
      INTEGER(4),INTENT(IN) :: LEN
      INTEGER(4),INTENT(IN) :: N1
      INTEGER(4),INTENT(IN) :: N2
      REAL(8)   ,INTENT(IN) :: PSI1(LEN,N1)
      REAL(8)   ,INTENT(IN) :: PSI2(LEN,N2)
      REAL(8)   ,INTENT(OUT):: OVERLAP(N1,N2)
      INTEGER(4)            :: I,J,K
      REAL(8)               :: SUM
!     ******************************************************************
      IF(TID.AND.N1.NE.N2) THEN
        CALL ERROR$MSG('PSI2 AND PSI1 DIFFER FOR TID=.TRUE.')
        CALL ERROR$STOP('LIB$SCALARPRODUCTR8')
      END IF
      IF(TID) THEN
!       ==  OVERLAP(I,J) = 0.D0*OVERLAP+1.D0*SUM_K:PSI1(K,I)*PSI1(K,J) =
        CALL DSYRK('U','T',N1,LEN,1.D0,PSI1,LEN,0.D0,OVERLAP,N1)
        DO I=1,N1
          DO J=I+1,N2
            OVERLAP(J,I)=OVERLAP(I,J)
          ENDDO
        ENDDO
      ELSE
#IF DEFINED(CPPVAR_BLAS_ESSL)
        CALL DGEMUL(PSI1,LEN,'T',PSI2,LEN,'N',OVERLAP,N1,N1,LEN,N2)
#ELSE 
         CALL DGEMM('T','N',N1,N2,LEN,1.D0,PSI1(:,:),LEN,PSI2(:,:),LEN &
     &             ,0.D0,OVERLAP,N1)
#ENDIF
      END IF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$SCALARPRODUCTC8(TID,LEN,N1,PSI1,N2,PSI2,OVERLAP)
!     ******************************************************************
!     ** PERFORMS THE SCALAR PRODUCTS OF TWO ARRAYS OF COMPLEX        **
!     ** STATE VECTORS                                                **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TID
      INTEGER(4),INTENT(IN) :: LEN
      INTEGER(4),INTENT(IN) :: N1
      INTEGER(4),INTENT(IN) :: N2
      COMPLEX(8),INTENT(IN) :: PSI1(LEN,N1)
      COMPLEX(8),INTENT(IN) :: PSI2(LEN,N2)
      COMPLEX(8),INTENT(OUT):: OVERLAP(N1,N2)
      INTEGER(4)            :: I,J,K
      COMPLEX(8)            :: SUM
!     ******************************************************************
      IF(TID.AND.N1.NE.N2) THEN
        CALL ERROR$MSG('PSI2 AND PSI1 DIFFER FOR TID=.TRUE.')
        CALL ERROR$STOP('LIB$SCALARPRODUCTC8')
      END IF
      OVERLAP(:,:)=(0.D0,0.D0)
      IF(TID) THEN
!       == ATTENTION: SCALAR FACTORS ARE SUPPOSED TO BE REAL AS THEY ARE
        CALL ZHERK('U','C',N1,LEN,1.D0,PSI1,LEN,0.D0,OVERLAP,N1)
        DO I=1,N1
          DO J=I+1,N2
            OVERLAP(J,I)=CONJG(OVERLAP(I,J))
          ENDDO
        ENDDO
      ELSE
#IF DEFINED(CPPVAR_BLAS_ESSL)
        CALL ZGEMUL(PSI1,LEN,'C',PSI2,LEN,'N',OVERLAP,N1,N1,LEN,N2)
#ELSE 
        CALL ZGEMM('C','N',N1,N2,LEN,(1.D0,0.D0),PSI1,LEN,PSI2,LEN &
     &            ,(0.D0,0.D0),OVERLAP,N1)
#ENDIF
      END IF
      RETURN
      END
!
!ZAXPY(N,ALPHA,X,INCX,Y,INCY)  Y=Y+ALPHA*X P267
!     ..................................................................
      SUBROUTINE LIB$VECTORADDC8(N,X,FAC,Y)
!     ******************************************************************
!     **  ADD TWO VECTORS                                             **
!     **  X=X+FAC*Y  WHERE X,Y ARE VECTORS AND FAC IS A SCALAR        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N
      COMPLEX(8),INTENT(INOUT) :: X(N)
      COMPLEX(8),INTENT(IN)    :: FAC
      COMPLEX(8),INTENT(IN)    :: Y(N)
      INTEGER(4)               :: I
      CHARACTER(8),PARAMETER   :: LIB='ESSL'
!     ******************************************************************
      CALL ZAXPY(N,FAC,X,1,Y,1)
      RETURN
      END
!
!DAXPY(N,ALPHA,X,INCX,Y,INCY)  Y=Y+ALPHA*X P267
!DVES(N,X,INCX,Y,INCY,X,INCZ)  Z=X-Y       P313
!     ..................................................................
      SUBROUTINE LIB$VECTORADDR8(N,X,FAC,Y)
!     ******************************************************************
!     **  ADD TWO VECTORS                                             **
!     **  X=X+FAC*Y  WHERE X,Y ARE VECTORS AND FAC IS A SCALAR        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N
      REAL(8)   ,INTENT(INOUT) :: X(N)
      REAL(8)   ,INTENT(IN)    :: FAC
      REAL(8)   ,INTENT(IN)    :: Y(N)
      INTEGER(4)               :: I
      CHARACTER(8),PARAMETER   :: LIB='ESSL'
!     ******************************************************************
      IF(FAC.EQ.-1.D0) THEN
#IF DEFINED(CPPVAR_BLAS_ESSL)
        CALL DVES(N,X,1,Y,1,X,1)
#ELSE
        CALL DAXPY(N,FAC,X,1,Y,1)
#ENDIF
      ELSE
        CALL DAXPY(N,FAC,X,1,Y,1)
      END IF
      RETURN
      END
!
!DGEF(A,LDA,N,IPVT)            MATRIX FACTORIZATION P507
!DGES(A,LDA,N,IPVT,BX)         BX=A^{-1}*BX (USES IPVT FROM DGEF)
!DGESVF                        SINGULAR VALUE DECOMPOSITION P696
!DGESVS  LEAST SQUARES SOLUTION USING SINGULAR VALUE DECOMPOSITION P703
!     ..................................................................
      SUBROUTINE LIB$MATRIXSOLVEnew(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  SOLVES THE LINEAR EQUATION SYSTEM AX=B                      **
!     **  WHERE A IS A(N,M)                                           **
!     **  IF A IS NOT SQUARE, THE EQUATION IS SOLVED IN A LEAST       **
!     **  SQUARE SENSE                                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(OUT):: X(M,NEQ)
      REAL(8)   ,INTENT(IN) :: B(N,NEQ)
      real(8)               :: a1(n,m)
      real(8)               :: b1(n,neq)
      integer               :: info
      real(8)               :: rcond=1.d-12
      integer               :: ldwork
      integer               :: irank
      real(8)               :: sing(n)
      real(8)   ,allocatable:: work(:)
      integer               :: n1,m1,neq1
!     ******************************************************************
      ldwork=3*min(M,N)+max(2*min(M,N),max(M,N),Neq)
      allocate(work(ldwork))
      n1=n
      m1=m
      neq1=neq
      a1=a 
      b1=b
      call dgelss(n1,m1,neq1,a1,n,b1,n1,sing,rcond,irank,work,ldwork,info)
      x=b1
      deallocate(work)
      if(info.gt.0) then
        call error$msg('ith argument has illegal value')
        call error$i4val('i',info)
        call error$stop('LIB$MATRIXSOLVEnew')
      else if(info.lt.0) then
        call error$msg('singlar value decomposition not converged')
        call error$msg('i off-diagonal elements of intermediate')
        call error$msg('bidiagonal form did not converge to zero')
        call error$i4val('i',info)
        call error$stop('LIB$MATRIXSOLVEnew')
      end if
      return
      end
!
!     ..................................................................
      SUBROUTINE LIB$MATRIXSOLVE(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  SOLVES THE LINEAR EQUATION SYSTEM AX=B                      **
!     **  WHERE A IS A(N,M)                                           **
!     **  IF A IS NOT SQUARE, THE EQUATION IS SOLVED IN A LEAST       **
!     **  SQUARE SENSE                                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(OUT):: X(M,NEQ)
      REAL(8)   ,INTENT(IN) :: B(N,NEQ)
      REAL(8)   ,ALLOCATABLE:: AFACT(:,:)
      REAL(8)   ,ALLOCATABLE:: BFACT(:,:)
      INTEGER(4),ALLOCATABLE:: IPVT(:)
      REAL(8)   ,ALLOCATABLE:: AUX(:)
      REAL(8)   ,ALLOCATABLE:: S(:)
      INTEGER(4)            :: NAUX
      INTEGER(4)            :: NM
      INTEGER(4)            :: IRANK
      REAL(8)               :: TAU=1.D-6
      INTEGER(4)            :: INFO
      INTEGER(4)            :: I
!     ******************************************************************
      IF(N.EQ.M) THEN  !MATRIX FACTORIZATION
        ALLOCATE(AFACT(N,N))
        ALLOCATE(IPVT(N)) 
        AFACT=A
#IF DEFINED(CPPVAR_BLAS_ESSL)
        CALL DGEF(AFACT,N,N,IPVT)
        X=B
        CALL DGES(AFACT,N,N,IPVT,X,0)
#ELSE 
        CALL DGEFA(AFACT,N,N,IPVT,INFO)
        X=B
        CALL DGESL(AFACT,N,N,IPVT,X,0)
#ENDIF
        DEALLOCATE(AFACT)
        DEALLOCATE(IPVT)
      ELSE   !SINGULAR VALUE DECOMPOSITION
        NM=MAX(1,MAX(N,M))
#IF DEFINED(CPPVAR_BLAS_ESSL)
        ALLOCATE(AFACT(NM,M))
        ALLOCATE(BFACT(N,NEQ))
        AFACT(1:N,:)=A
        AFACT(N+1:,:)=0.D0
        BFACT=B
        NAUX=2*M+MAX(NM,NEQ)
        ALLOCATE(S(M))    !SINGULAR VALUES
        ALLOCATE(AUX(NAUX))
        CALL DGESVF(2,AFACT,N,BFACT,N,NEQ,S,N,M,AUX,NAUX) !->ESSL
        DEALLOCATE(AUX)
        CALL DGESVS(AFACT,NM,BFACT,N,NEQ,S,X,M,N,M,TAU) !->ESSL
        DEALLOCATE(S)
        DEALLOCATE(BFACT)
        DEALLOCATE(AFACT)
#ELSE
        ALLOCATE(AFACT(N,M))
        AFACT(:,:)=A(:,:)
        ALLOCATE(BFACT(NM,NEQ))
        BFACT(1:N,:)=B(:,:)
        ALLOCATE(S(M))    !SINGULAR VALUES
        NAUX=3*NM+MAX(2*MIN(M,N),NM,NEQ)
        ALLOCATE(AUX(NAUX))
        CALL DGELSS(N,M,NEQ,AFACT,NM,BFACT,N,S,TAU,IRANK,AUX,NAUX,INFO) !->LAPACK
        IF(INFO.NE.0) THEN
          IF(INFO.LT.0) THEN
            CALL ERROR$MSG('INCORRECT ARGUMENT FOR SINGULAR VALUE DECOMPOSITION')
            CALL ERROR$I4VAL('WRONG ARGUMENT NR.',-I)
          ELSE IF(INFO.GT.0) THEN
            CALL ERROR$MSG('SINGULAR VALUE DECOMPOSITION FAILED TO CONVERGE')
            CALL ERROR$I4VAL('#(SING. VALUES WITHOUT CONVERGENCE)',I)
          END IF
          CALL ERROR$STOP('LIB$MATRIXSOLVE')
        END IF
        DEALLOCATE(AUX)
        X(:,:)=BFACT(1:M,:)
        DEALLOCATE(S)
        DEALLOCATE(BFACT)
        DEALLOCATE(AFACT)
#ENDIF
      END IF
      RETURN
      END
!
!ISORT  SORTS ELEMENTS IN A SEQUENCE P904
!     .................................................................
      SUBROUTINE LIB$SORTI4(N,X)
!     ******************************************************************
!     **  SORTS THE ARRAY X IN ASCENDING ORDER                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N
      INTEGER(4),INTENT(INOUT):: X(N)
      INTEGER(4)              :: IHELP
      INTEGER(4)              :: I,J
!     ******************************************************************
#IF DEFINED(CPPVAR_BLAS_ESSL)
      CALL ISORT(X,1,N)
#ELSE
      DO I=1,N-1
        DO J=I+1,N
          IF(X(I).GT.X(J)) THEN
            IHELP=X(I)
            X(I) =X(J)
            X(J) =IHELP
          ENDIF
        ENDDO
      ENDDO
#ENDIF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$GETHOSTNAME(HOSTNAME)
!     *********************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                **
!     *********************************************************************
      CHARACTER(*),INTENT(OUT)  :: HOSTNAME
      INTEGER(4)                :: RC
!     *********************************************************************
#IF DEFINED(CPPVAR_SUPPORT_XLF)
      RC=HOSTNM_(HOSTNAME)    ! XLF SUPPORT LIBRARY
      IF(RC.NE.0)HOSTNAME='UNKNOWN'
#ELSE
      HOSTNAME='UNKNOWN'
!     HOSTNM_=GETHOSTNAME(HOSTNAME,%VAL(LEN(HOSTNAME)))
#ENDIF
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE LIB$FLUSHFILE(N)
!     *********************************************************************
!     ** FLUSHES THE BUFFER FOR THE FILE CONNECTED TO FORTRAN UNIT N     **
!     *********************************************************************
      INTEGER(4),INTENT(IN) :: N
!     *********************************************************************
#IF DEFINED(CPPVAR_SUPPORT_XLF)
      CALL FLUSH_(N)  ! XLF USPPORT LIBRARY
#ELIF DEFINED(CPPVAR_U77)
      CALL FLUSH(N) ! FROM ABSOFT SUPPORT LIBRARY (UNDERSCORE)
#ENDIF
      RETURN
      END 
!
#IFNDEF CPPVAR_BLAS_ESSL
!
!     ..................................................................
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
!     ******************************************************************
!     ** DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN          **
!     ** ELIMINATION.                                                 **
!     **                                                               **
!     ** DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED       **
!     ** DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.     **
!     ** (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA).               **
!     **                                                              **
!     **  ON RETURN, A IS AN UPPER TRIANGULAR MATRIX                  **
!     **  AND THE MULTIPLIERS, WHICH WERE USED TO OBTAIN IT.          **
!     **  THE FACTORIZATION CAN BE WRITTEN  A = L*U,                  **
!     **  WHERE L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER         **
!     **  TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.            **
!     **                                                              **
!     **  INFO                                                        **
!     **     = 0  NORMAL VALUE.                                       **
!     **     = K  IF U(K,K).EQ.0.0 , THIS IS NOT AN ERROR CONDITION   **
!     **          FOR THIS SUBROUTINE, BUT IT DOES INDICATE           **
!     **          THAT DGESL OR DGEDI WILL DIVIDE BY ZERO IF CALLED.  **
!     **          USE RCOND IN DGECO FOR A RELIABLE INDICATION        **
!     **          OF SINGULARITY.                                     **
!     **                                                              **
!     **  SUBROUTINES AND FUNCTIONS CALLED:                           **
!     **     BLAS ROUTINES: DAXPY,DSCAL,IDAMAX                        **
!     **                                                              **
!     **  LINPACK. THIS VERSION DATED 08/14/78 .                      **
!     **  CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.**
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N        ! ORDER OF MATRIX A
      INTEGER(4),INTENT(IN)   :: LDA      ! LEADING DIM. OF MATRIX A
      REAL(8)   ,INTENT(INOUT):: A(LDA,N) ! MATRIX TO BE FACTORED
      INTEGER(4),INTENT(OUT)  :: IPVT(N)  ! PIVOT INDICES
      INTEGER(4),INTENT(OUT)  :: INFO     ! RETURN CODE
      REAL(8)                 :: T
      INTEGER(4)              :: IDAMAX,J,K,KP1,L,NM1
!     ******************************************************************
!     ==================================================================
!     ==  GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                  ==
!     ==================================================================
      INFO=0
      NM1 =N-1
      IF(NM1.LT.1) THEN
        IPVT(N)=N
         IF(A(N,N).EQ.0.0D0) INFO=N
        RETURN
      END IF
!
      DO K=1,NM1
        KP1=K+1
!       ==  FIND L = PIVOT INDEX ======================================
        L=IDAMAX(N-K+1,A(K,K),1) + K-1  !IDAMAX<-BLAS LIBRARY
        IPVT(K)=L
!       == ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED ======
        IF(A(L,K).NE.0.0D0) THEN
!         ==  INTERCHANGE IF NECESSARY ===============================
          IF(L.NE.K) THEN
            T=A(L,K)
            A(L,K)=A(K,K)
            A(K,K)=T
          END IF
!         == COMPUTE MULTIPLIERS ======================================
          T=-1.0D0/A(K,K)
          CALL DSCAL(N-K,T,A(K+1,K),1)
!         ==  ROW ELIMINATION WITH COLUMN INDEXING  ===================
          DO J=KP1,N
            T=A(L,J)
            IF(L.NE.K) THEN
              A(L,J)=A(K,J)
              A(K,J)=T
            END IF
            CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
          ENDDO
        ELSE
          INFO=K
        END IF
      ENDDO
      IPVT(N)=N
      IF(A(N,N).EQ.0.0D0) INFO=N
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
!     ******************************************************************
!     ** DGESL SOLVES THE DOUBLE PRECISION SYSTEM                     **
!     **       A*X=B  OR  TRANS(A)*X=B                                **
!     **  USING THE FACTORS COMPUTED BY DGECO OR DGEFA.               **
!     **                                                              **
!     **  JOB=0       TO SOLVE  A*X = B ,                             **
!     **     =NONZERO TO SOLVE  TRANS(A)*X = B,                       **
!     **                        ,WHERE TRANS(A)  IS THE TRANSPOSE.    **
!     **                                                              **
!     **  ERROR CONDITION: A DIVISION BY ZERO WILL OCCUR              **
!     **    IF THE INPUT FACTOR  CONTAINS A ZERO ON THE DIAGONAL.     **
!     **    TECHNICALLY THIS INDICATES SINGULARITY BUT IT IS OFTEN    **
!     **    CAUSED BY IMPROPER ARGUMENTS OR IMPROPER SETTING OF LDA.  **
!     **    IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY **
!     **    AND IF DGECO HAS SET RCOND>0. OR DGEFA HAS SET INFO=0.    **
!     **                                                              **
!     **  TO COMPUTE  INVERSE(A)*C WHERE C IS A MATRIX WITH P COLUMNS **
!     **    CALL DGECO(A,LDA,N,IPVT,RCOND,Z)                          **
!     **    IF (RCOND IS TOO SMALL) THEN                              **
!     **      PRINT ERROR MESSAGE AND STOP...                         **
!     **    END IF                                                    **
!     **    DO J=1,P                                                  **
!     **      CALL DGESL(A,LDA,N,IPVT,C(1,J),0)                       **
!     **    ENDDO                                                     **
!     **                                                              **
!     **   SUBROUTINES AND FUNCTIONS CALLED                           **
!     **      BLAS ROUTINES: DAXPY,DDOT                               **
!     **                                                              **
!     **  LINPACK. THIS VERSION DATED 08/14/78 .                      **
!     **  CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.**
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N        ! ORDER OF A
      INTEGER(4),INTENT(IN)   :: LDA      ! LEADING DIMENSION OF A
      REAL(8)   ,INTENT(IN)   :: A(LDA,N) ! OUTPUT OF DGECO OR DGEFA
      REAL(8)   ,INTENT(INOUT):: B(N)     ! INPUT:  RIGHT HAND SIDE B
                                          ! OUTPUT: SOLUTION VECTOR X
      INTEGER(4),INTENT(IN)   :: IPVT(N)  ! PIVOT VECTOR FROM DGECO OR DEGFA
      INTEGER(4),INTENT(IN)   :: JOB      ! A*X=B OR TRANS(A)X=B
      REAL(8)                 :: DDOT
      REAL(8)                 :: T
      INTEGER(4)              :: K,KB,L,NM1
!     ******************************************************************
      NM1=N-1
!     ==================================================================
!     == JOB=0: SOLVE A*X=B                                           ==
!     ==================================================================
      IF(JOB.EQ.0) THEN
!       == FIRST SOLVE  L*Y = B ========================================
        IF(NM1.GE.1) THEN
          DO K=1,NM1
            L = IPVT(K)
            T = B(L)
            IF(L.NE.K) THEN
              B(L) = B(K)
              B(K) = T
            END IF
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
          ENDDO
        END IF
!       ==  NOW SOLVE  U*X = Y ==========================================
        DO KB=1,N
          K =N+1-KB
          B(K)=B(K)/A(K,K)
          T =-B(K)
          CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
        ENDDO
!
!     ==================================================================
!     == JOB = NONZERO, SOLVE  TRANS(A) * X = B                       ==
!     ==================================================================
     ELSE
        
!       == FIRST SOLVE  TRANS(U)*Y=B ===================================
        DO K = 1, N
          T = DDOT(K-1,A(1,K),1,B(1),1)
          B(K) = (B(K) - T)/A(K,K)
        ENDDO
!       == NOW SOLVE TRANS(L)*X = Y =====================================
        IF (NM1.GE.1) THEN
          DO KB=1,NM1
            K = N-KB
            B(K)=B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF(L.NE.K) THEN
              T = B(L)
              B(L) = B(K)
              B(K) = T
            END IF
          ENDDO
        END IF
      END IF
      RETURN
      END
#ENDIF
#IF DEFINED(CPPVAR_FFT_PACK)
!
!     ..................................................................
      SUBROUTINE CFFTB(N,C,WA,IFAC)
!     ******************************************************************
!     ** CFFTB COMPUTES THE BACKWARD COMPLEX DISCRETE FOURIER         **
!     ** TRANSFORM (THE FOURIER SYNTHESIS).                           **
!     ** EQUIVALENTLY, CFFTB COMPUTES A COMPLEX PERIODIC SEQUENCE     **
!     ** FROM ITS FOURIER COEFFICIENTS.                               **
!     **                                                              **
!     **  THE TRANSFORM IS NOT NORMALIZED. TO OBTAIN A NORMALIZED     **
!     **  TRANSFORM THE OUTPUT OF CFFTF! MUST BE DIVIDED BY N.        **
!     **  OTHERWISE A CALL OF CFFTF FOLLOWED BY A CALL OF CFFTB       **
!     **  WILL MULTIPLY THE SEQUENCE BY N.                            **
!     **                                                              **
!     **  THE ARRAYS WA AND IFAC WHICH ARE USED BY SUBROUTINE CFFTF   **
!     **  MUST BE INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE).   **
!     **  THE ARRAYS DEPEND ON THE LENGTH N OF THE TRANSFORM          **
!     **  THE SAME ARRAYS CAN BE USED FOR FORWARD AN BACK TRANSFORM   **
!     **  THE WORK ARRAYS ARE NOT CHANGED BY THE FOURIER TANSFORM     **
!     **  ROUTINES CFFTF AND CFFTB. AND CAN BE KEPT FOR SUBSEQUENT    **
!     **  TRANSFORMS WITH THE SAME LENGTH                             **
!     **                                                              **
!     **  THE METHOD IS MORE EFFICIENT IF THE LENGTH OF THE TRANSFORM **
!     **  IS A PRODUCT OF SMALL PRIME NUMBERS                         **
!     **                                                              **
!     **  THE COMPLEX SEQUENCE C OF LENGTH N IS TREATED INTERNALLY    **
!     **  AS A REAL ARRAY OF LENGTH 2*N                               **
!     **                                                              **
!     **  RESULT                                                      **
!     **     FOR K,J=1,...,N AND I=SQRT(-1)                           **
!     **     C(J)=SUM_{K=1}^N C(K)*EXP(+I*(J-1)*(K-1)*2*PI/N)         **
!     **                                                              **
!     **  THE ARRAYS WA,IFAC CONTAIN INITIALIZATION CALCULATIONS      **
!     **  WHICH MUST NOT BE DESTROYED BETWEEN CALLS                   **
!     **  OF SUBROUTINE CFFTF OR CFFTB                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N
      REAL(8)   ,INTENT(INOUT):: C(2*N)
      REAL(8)   ,INTENT(IN)   :: WA(2*N)
      INTEGER(4),INTENT(IN)   :: IFAC(30)
      REAL(8)                 :: CH(2*N)
      INTEGER(4)              :: NF,NA,L1,IW,L2,IDO,IDOT,IDL1,N2
      INTEGER(4)              :: IP
      INTEGER(4)              :: IX2,IX3,IX4
      INTEGER(4)              :: K1
      INTEGER(4)              :: I
      INTEGER(4)              :: NAC
!     ******************************************************************
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1=1,NF
        IP = IFAC(K1+2)
        L2 = IP*L1
        IDO = N/L2
        IDOT = IDO+IDO
        IDL1 = IDOT*L1
        IF(IP.EQ.4) THEN
          IX2 = IW+IDOT
          IX3 = IX2+IDOT
          IF(NA.EQ.0) THEN
            CALL PASSB4(IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
          ELSE
            CALL PASSB4(IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
          END IF
          NA = 1-NA
        ELSE IF(IP.EQ.2) THEN
          IF(NA.EQ.0) THEN
            CALL PASSB2(IDOT,L1,C,CH,WA(IW))
          ELSE
            CALL PASSB2(IDOT,L1,CH,C,WA(IW))
          END IF
          NA = 1-NA
        ELSE IF(IP.EQ.3) THEN
          IX2 = IW+IDOT
          IF(NA.EQ.0) THEN
            CALL PASSB3(IDOT,L1,C,CH,WA(IW),WA(IX2))
          ELSE
            CALL PASSB3(IDOT,L1,CH,C,WA(IW),WA(IX2))
          END IF
          NA = 1-NA
        ELSE IF(IP.EQ.5) THEN
          IX2 = IW+IDOT
          IX3 = IX2+IDOT
          IX4 = IX3+IDOT
          IF(NA.EQ.0) THEN
            CALL PASSB5(IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
          ELSE
            CALL PASSB5(IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
          END IF
          NA = 1-NA
        ELSE
          IF(NA.EQ.0) THEN
            CALL PASSB(NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
          ELSE
            CALL PASSB(NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
          END IF
          IF(NAC.NE.0) NA = 1-NA
        END IF
        L1 = L2
        IW = IW+(IP-1)*IDOT
      ENDDO
!
      IF(NA.NE.0) THEN
        N2 = N+N
        DO I=1,N2
          C(I) = CH(I)
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CFFTF(N,C,WA,IFAC)
!     ******************************************************************
!     **  CFFTF COMPUTES THE FORWARD COMPLEX DISCRETE FOURIER         **
!     **  TRANSFORM (THE FOURIER ANALYSIS).                           **
!     **  EQUIVALENTLY, CFFTF COMPUTES THE FOURIER COEFFICIENTS       **
!     **  OF A COMPLEX PERIODIC SEQUENCE.                             **
!     **                                                              **
!     **  THE TRANSFORM IS NOT NORMALIZED. TO OBTAIN A NORMALIZED     **
!     **  TRANSFORM THE OUTPUT MUST BE DIVIDED BY N.                  **
!     **  OTHERWISE A CALL OF CFFTF FOLLOWED BY A CALL OF CFFTB       **
!     **  WILL MULTIPLY THE SEQUENCE BY N.                            **
!     **                                                              **
!     **  THE ARRAYS WA AND IFAC WHICH ARE USED BY SUBROUTINE CFFTF   **
!     **  MUST BE INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE).   **
!     **  THE ARRAYS DEPEND ON THE LENGTH N OF THE TRANSFORM          **
!     **  THE SAME ARRAYS CAN BE USED FOR FORWARD AN BACK TRANSFORM   **
!     **  THE WORK ARRAYS ARE NOT CHANGED BY THE FOURIER TANSFORM     **
!     **  ROUTINES CFFTF AND CFFTB. AND CAN BE KEPT FOR SUBSEQUENT    **
!     **  TRANSFORMS WITH THE SAME LENGTH                             **
!     **                                                              **
!     **  THE METHOD IS MORE EFFICIENT IF THE LENGTH OF THE TRANSFORM **
!     **  IS A PRODUCT OF SMALL PRIME NUMBERS                         **
!     **                                                              **
!     **  THE COMPLEX SEQUENCE C OF LENGTH N IS TREATED INTERNALLY    **
!     **  AS A REAL ARRAY OF LENGTH 2*N                               **
!     **                                                              **
!     **  RESULT                                                      **
!     **     FOR K,J=1,...,N AND I=SQRT(-1)                           **
!     **     C(J)=SUM_{K=1}^N C(K)*EXP(-I*(J-1)*(K-1)*2*PI/N)         **
!     **                                                              **
!     **  THE ARRAYS WA,IFAC CONTAIN INITIALIZATION CALCULATIONS      **
!     **  WHICH MUST NOT BE DESTROYED BETWEEN CALLS                   **
!     **  OF SUBROUTINE CFFTF OR CFFTB                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N        ! LENGTH OF FFT
      REAL(8)   ,INTENT(INOUT) :: C(2*N)   ! SEQUENCE TO BE TRANSFORMED
      REAL(8)   ,INTENT(IN)    :: WA(2*N)  ! WORK ARRAY
      INTEGER(4),INTENT(IN)    :: IFAC(30) ! WORK ARRAY
      REAL(8)                  :: CH(2*N)  ! WORK ARRAY
      INTEGER(4)               :: NF,NA,L1,IW
      INTEGER(4)               :: IP,L2,IDO,IDOT,IDL1,IX2,IX3,IX4
      INTEGER(4)               :: N2
      INTEGER(4)               :: I,K1
      INTEGER(4)               :: NAC
!     ******************************************************************
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1=1,NF
        IP = IFAC(K1+2)
        L2 = IP*L1
        IDO = N/L2
        IDOT = IDO+IDO
        IDL1 = IDOT*L1
        IF(IP.EQ.4) THEN
          IX2 = IW+IDOT
          IX3 = IX2+IDOT
          IF(NA.EQ.0) THEN
            CALL PASSF4(IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
          ELSE
            CALL PASSF4(IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
          END IF
          NA = 1-NA
        ELSE IF(IP.EQ.2) THEN
          IF(NA.EQ.0) THEN
            CALL PASSF2(IDOT,L1,C,CH,WA(IW))
          ELSE
            CALL PASSF2(IDOT,L1,CH,C,WA(IW))
          END IF
          NA = 1-NA
        ELSE IF(IP.EQ.3) THEN
          IX2 = IW+IDOT
          IF(NA.EQ.0) THEN
            CALL PASSF3(IDOT,L1,C,CH,WA(IW),WA(IX2))
          ELSE
            CALL PASSF3(IDOT,L1,CH,C,WA(IW),WA(IX2))
          END IF
          NA = 1-NA
        ELSE IF(IP.EQ.5) THEN
          IX2 = IW+IDOT
          IX3 = IX2+IDOT
          IX4 = IX3+IDOT
          IF(NA.EQ.0) THEN
            CALL PASSF5(IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
          ELSE
            CALL PASSF5(IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
          ENDIF
          NA = 1-NA
        ELSE
          IF(NA.NE.0) THEN
            CALL PASSF(NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
          ELSE
            CALL PASSF(NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
          END IF
          IF (NAC .NE. 0) NA = 1-NA
        END IF
        L1 = L2
        IW = IW+(IP-1)*IDOT
      ENDDO
      IF(NA.NE.0) THEN
        N2 = N+N
        DO I=1,N2
          C(I) = CH(I)
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CFFTI(N,WA,IFAC)
!     ******************************************************************
!     **  CFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN BOTH     **
!     **  CFFTF AND CFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH **
!     **  A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED    **
!     **  AND STORED IN WA AND IFAC.                                  **
!     **                                                              **
!     **  THE ARRAYS WA,IFAC CONTAIN INITIALIZATION CALCULATIONS      **
!     **  WHICH MUST NOT BE DESTROYED BETWEEN CALLS                   **
!     **  OF SUBROUTINE CFFTF OR CFFTB                                **
!     **                                                              **
!     ** IFAC CONTAINS (1) THE LENGTH OF THE TRANSFORM,               **
!     **               (2) THE NUMBER OF FACTORS AND                  **
!     **           AND (3,NF+2) NF FACTORS                            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N        ! LENGTH OF THE FOURIER TRANSFORM
      INTEGER(4),INTENT(OUT) :: IFAC(30) ! PRIME FACTORIZATION
      REAL(8)   ,INTENT(OUT) :: WA(2*N)  ! TRIGONOMETRIC FUNCTIONS
      INTEGER(4),PARAMETER   :: NTRYH(4)=(/3,4,2,5/)
      INTEGER(4)             :: NL,NF,NQ,NR
      INTEGER(4)             :: NTRY,IB,L1,IP,LD,L2,IDO,IDOT,IPM,I1
      INTEGER(4)             :: I,J,II,K1
      REAL(8)                :: TPI,ARGH,FI,ARGLD,ARG
      REAL(8)                :: PI
!     ******************************************************************
!
!     ==================================================================
!     == PRIME FACTORIZATION                                          ==
!     ==================================================================
      NL = N    ! REMAINDER OF FACTORED LENGTH OF TRANSFORM
      NF = 0    !#(FACTORS)
      J=0
 101  CONTINUE
!     __ SELECT A NEW PRIME NUMBER NTRY FROM THE ARRAY 3,4,2,5,7,9,11,...
      J=J+1
      IF(J.LE.4) THEN
        NTRY=NTRYH(J)
      ELSE
        NTRY=NTRY+2
      END IF
104   CONTINUE
!     __  CHECK IF PRIME NUMBER IS A FACTOR_____________________________
      NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF(NR.NE.0) GOTO 101 ! IF NTRY IS NOT A FACTOR PICK A LARGER PRIME
!     __NEW FACTOR FOUND, ADD TO ARRAY IFAC_____________________________
      NF = NF+1
      IF(NF.GT.13) THEN
        CALL ERROR$MSG('ARRAY SIZE FOR FACTORIZATION TOO SMALL')
        CALL ERROR$MSG('INCREASE ARRAY SIZE FOR IFAC IN CFFTI,CFFTB,CFFTF')
        CALL ERROR$STOP('CFFTI')
      END IF
      IFAC(NF+2) = NTRY
      NL = NQ
!     __PLACE ANY PRIME FACTOR 2 UPFRONT________________________________
      IF(NTRY.EQ.2.AND.NF.NE.1) THEN  
        DO I=2,NF
          IB = NF-I+2
          IFAC(IB+2) = IFAC(IB+1)
        ENDDO
        IFAC(3) = 2
      END IF
      IF(NL.NE.1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
!
!     ==================================================================
!     == TABULATE TRIGONOMETRIC FUNCTIONS                             ==
!     ==================================================================
!     PI=3.14159265358979323846264338327950288419716939937510582097494D0
      PI=4.D0*DATAN(1.D0)
      TPI = 2.0D0*PI
      ARGH = TPI/FLOAT(N)
      I = 2
      L1 = 1
      DO K1=1,NF
        IP = IFAC(K1+2)
        LD = 0
        L2 = L1*IP
        IDO = N/L2
        IDOT = IDO+IDO+2
        IPM = IP-1
        DO J=1,IPM
          I1 = I
          WA(I-1) = 1.0D0
          WA(I)   = 0.0D0
          LD = LD+L1
          FI = 0.0D0
          ARGLD = FLOAT(LD)*ARGH
          DO II=4,IDOT,2
            I = I+2
            FI = FI+1.0D0
            ARG = FI*ARGLD
            WA(I-1)=COS(ARG)
            WA(I)  =SIN(ARG)
          ENDDO
          IF(IP.GT.5) THEN
            WA(I1-1) = WA(I-1)
            WA(I1)   = WA(I)
          END IF
        ENDDO
        L1 = L2
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSB(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT)  :: NAC
      INTEGER(4),INTENT(IN)   :: IDO
      INTEGER(4),INTENT(IN)   :: IP
      INTEGER(4),INTENT(IN)   :: L1
      INTEGER(4),INTENT(IN)   :: IDL1
      REAL(8)   ,INTENT(OUT)  :: CH(IDO,L1,IP)
      REAL(8)   ,INTENT(IN)   :: CC(IDO,IP,L1)
      REAL(8)   ,INTENT(OUT)  :: C1(IDO,L1,IP)
      REAL(8)   ,INTENT(OUT)  :: C2(IDL1,IP)
      REAL(8)   ,INTENT(INOUT):: CH2(IDL1,IP)
      REAL(8)   ,INTENT(IN)   :: WA(*) 
!MAX(1+(IP-1)*IDO,2+((IP-1)/2)*IDO)??
      INTEGER(4)              :: IDOT,NT,IPP2,IPPH,IDP,IDL,INC,IDLJ,IDIJ,IDJ
      INTEGER(4)              :: JC,LC
      INTEGER(4)              :: I,J,K,L,IK
      REAL(8)                 :: WAR,WAI
!     ******************************************************************
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
!
      IF(IDO.GE.L1) THEN
        DO J=2,IPPH
          JC = IPP2-J
          DO K=1,L1
            DO I=1,IDO
              CH(I,K,J)  = CC(I,J,K)+CC(I,JC,K)
              CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
            ENDDO
          ENDDO
        ENDDO
        DO K=1,L1
          DO I=1,IDO
            CH(I,K,1) = CC(I,1,K)
          ENDDO
        ENDDO
      ELSE
        DO J=2,IPPH
          JC = IPP2-J
          DO I=1,IDO
            DO K=1,L1
              CH(I,K,J)  = CC(I,J,K)+CC(I,JC,K)
              CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
            ENDDO
          ENDDO
        ENDDO
        DO I=1,IDO
          DO K=1,L1
            CH(I,K,1) = CC(I,1,K)
          ENDDO
        ENDDO
      ENDIF
!
      IDL = 2-IDO
      INC = 0
      DO L=2,IPPH
        LC = IPP2-L
        IDL = IDL+IDO
        DO IK=1,IDL1
          C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
          C2(IK,LC)=           WA(IDL)  *CH2(IK,IP)
        ENDDO
        IDLJ = IDL
        INC = INC+IDO
        DO J=3,IPPH
          JC = IPP2-J
          IDLJ = IDLJ+INC
          IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
          WAR = WA(IDLJ-1)
          WAI = WA(IDLJ)
          DO IK=1,IDL1
            C2(IK,L)  = C2(IK,L) +WAR*CH2(IK,J)
            C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
          ENDDO
        ENDDO
      ENDDO
!
      DO J=2,IPPH
        DO IK=1,IDL1
          CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
        ENDDO
      ENDDO
!
      DO J=2,IPPH
        JC = IPP2-J
        DO IK=2,IDL1,2
          CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
          CH2(IK-1,JC)= C2(IK-1,J)+C2(IK,JC)
          CH2(IK,J)   = C2(IK,J)  +C2(IK-1,JC)
          CH2(IK,JC)  = C2(IK,J)  -C2(IK-1,JC)
        ENDDO
      ENDDO
!
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO IK=1,IDL1
        C2(IK,1) = CH2(IK,1)
      ENDDO
      DO J=2,IP
        DO K=1,L1
          C1(1,K,J) = CH(1,K,J)
          C1(2,K,J) = CH(2,K,J)
        ENDDO
      ENDDO
!
      IF(IDOT.LE.L1) THEN
        IDIJ = 0
        DO J=2,IP
         IDIJ = IDIJ+2
         DO I=4,IDO,2
           IDIJ = IDIJ+2
           DO K=1,L1
             C1(I-1,K,J)=WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
             C1(I,K,J)  =WA(IDIJ-1)*CH(I,K,J)  +WA(IDIJ)*CH(I-1,K,J)
           ENDDO
         ENDDO
       ENDDO
      ELSE
        IDJ = 2-IDO
        DO J=2,IP
          IDJ = IDJ+IDO
          DO K=1,L1
            IDIJ = IDJ
            DO I=4,IDO,2
              IDIJ = IDIJ+2
              C1(I-1,K,J)= WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
              C1(I,K,J)  = WA(IDIJ-1)*CH(I,K,J)  +WA(IDIJ)*CH(I-1,K,J)
            ENDDO
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSB2(IDO,L1,CC,CH,WA1)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,2,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,2)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      INTEGER(4)            :: K,I
      REAL(8)               :: TR2,TI2
!     ******************************************************************
      IF(IDO.LE.2) THEN
        DO K=1,L1
          CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
          CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
          CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
          CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSB3 (IDO,L1,CC,CH,WA1,WA2)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,3,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,3)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      REAL(8)   ,INTENT(IN) :: WA2(IDO)
      REAL(8)               :: TR2,TI2
      REAL(8)               :: CR2,CI2
      REAL(8)               :: CR3,CI3
      REAL(8)               :: DR2,DI2
      REAL(8)               :: DR3,DI3
      REAL(8)   ,PARAMETER  :: TAUR=-0.5D0
      REAL(8)   ,PARAMETER  :: TAUI=0.866025403784439D0
      INTEGER(4)            :: K,I
!     ******************************************************************
      IF(IDO.EQ.2) THEN
        DO K=1,L1
          TR2 = CC(1,2,K)+CC(1,3,K)
          CR2 = CC(1,1,K)+TAUR*TR2
          CH(1,K,1) = CC(1,1,K)+TR2
          TI2 = CC(2,2,K)+CC(2,3,K)
          CI2 = CC(2,1,K)+TAUR*TI2
          CH(2,K,1) = CC(2,1,K)+TI2
          CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
          CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
          CH(1,K,2) = CR2-CI3
          CH(1,K,3) = CR2+CI3
          CH(2,K,2) = CI2+CR3
          CH(2,K,3) = CI2-CR3
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSB4(IDO,L1,CC,CH,WA1,WA2,WA3)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,4,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,4)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      REAL(8)   ,INTENT(IN) :: WA2(IDO)
      REAL(8)   ,INTENT(IN) :: WA3(IDO)
      REAL(8)               :: TR1,TI1
      REAL(8)               :: TR2,TI2
      REAL(8)               :: TR3,TI3
      REAL(8)               :: TR4,TI4
      REAL(8)               :: CR2,CI2
      REAL(8)               :: CR3,CI3
      REAL(8)               :: CR4,CI4
      INTEGER(4)            :: I,K
!     ******************************************************************
      IF(IDO.EQ.2) THEN
        DO K=1,L1
          TI1 = CC(2,1,K)-CC(2,3,K)
          TI2 = CC(2,1,K)+CC(2,3,K)
          TR4 = CC(2,4,K)-CC(2,2,K)
          TI3 = CC(2,2,K)+CC(2,4,K)
          TR1 = CC(1,1,K)-CC(1,3,K)
          TR2 = CC(1,1,K)+CC(1,3,K)
          TI4 = CC(1,2,K)-CC(1,4,K)
          TR3 = CC(1,2,K)+CC(1,4,K)
          CH(1,K,1) = TR2+TR3
          CH(1,K,3) = TR2-TR3
          CH(2,K,1) = TI2+TI3
          CH(2,K,3) = TI2-TI3
          CH(1,K,2) = TR1+TR4
          CH(1,K,4) = TR1-TR4
          CH(2,K,2) = TI1+TI4
          CH(2,K,4) = TI1-TI4
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSB5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,5,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,5)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      REAL(8)   ,INTENT(IN) :: WA2(IDO)
      REAL(8)   ,INTENT(IN) :: WA3(IDO)
      REAL(8)   ,INTENT(IN) :: WA4(IDO)
      REAL(8)   ,PARAMETER  :: TR11= 0.309016994374947D0
      REAL(8)   ,PARAMETER  :: TI11= 0.951056516295154D0
      REAL(8)   ,PARAMETER  :: TR12=-0.809016994374947D0
      REAL(8)   ,PARAMETER  :: TI12= 0.587785252292473D0
      INTEGER(4)            :: I,K
      REAL(8)               :: TR2,TI2
      REAL(8)               :: TR3,TI3
      REAL(8)               :: TR4,TI4
      REAL(8)               :: TR5,TI5
      REAL(8)               :: CR2,CI2
      REAL(8)               :: CR3,CI3
      REAL(8)               :: CR4,CI4
      REAL(8)               :: CR5,CI5
      REAL(8)               :: DR2,DI2
      REAL(8)               :: DR3,DI3
      REAL(8)               :: DR4,DI4
      REAL(8)               :: DR5,DI5
!     ******************************************************************
      IF(IDO.EQ.2) THEN
      DO K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1)   = CC(I,1,K)  +TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)  +TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)  +TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2)   = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3)   = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4)   = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5)   = WA4(I-1)*DI5+WA4(I)*DR5
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSF(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT)  :: NAC
      INTEGER(4),INTENT(IN)   :: IDO
      INTEGER(4),INTENT(IN)   :: IP
      INTEGER(4),INTENT(IN)   :: L1
      INTEGER(4),INTENT(IN)   :: IDL1
      REAL(8)   ,INTENT(IN)   :: CC(IDO,IP,L1)
      REAL(8)   ,INTENT(OUT)  :: C1(IDO,L1,IP)
      REAL(8)   ,INTENT(OUT)  :: C2(IDL1,IP)
      REAL(8)   ,INTENT(OUT)  :: CH(IDO,L1,IP)
      REAL(8)   ,INTENT(INOUT):: CH2(IDL1,IP)
      REAL(8)   ,INTENT(IN)   :: WA(*)
      INTEGER(4)              :: IDOT,NT,IPP2,IPPH,IDP,INC
      INTEGER(4)              :: IDJ,IDIJ,IDL,IDLJ
      INTEGER(4)              :: JC,LC
      INTEGER(4)              :: I,J,K,L,IK
      REAL(8)                 :: WAR,WAI
!     ******************************************************************
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
!
      IF(IDO.GE.L1) THEN
        DO J=2,IPPH
          JC = IPP2-J
          DO K=1,L1
            DO I=1,IDO
              CH(I,K,J)  = CC(I,J,K)+CC(I,JC,K)
              CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
            ENDDO
          ENDDO
        ENDDO
        DO K=1,L1
          DO I=1,IDO
            CH(I,K,1) = CC(I,1,K)
          ENDDO
        ENDDO
      ELSE
        DO J=2,IPPH
          JC = IPP2-J
          DO I=1,IDO
            DO K=1,L1
              CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
              CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
            ENDDO
          ENDDO
        ENDDO
        DO I=1,IDO
          DO K=1,L1
            CH(I,K,1) = CC(I,1,K)
          ENDDO
        ENDDO
      END IF
!
      IDL = 2-IDO
      INC = 0
      DO L=2,IPPH
        LC = IPP2-L
        IDL = IDL+IDO
        DO IK=1,IDL1
          C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
          C2(IK,LC)=          -WA(IDL)  *CH2(IK,IP)
        ENDDO
        IDLJ = IDL
        INC = INC+IDO
        DO J=3,IPPH
          JC = IPP2-J
          IDLJ = IDLJ+INC
          IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
          WAR = WA(IDLJ-1)
          WAI = WA(IDLJ)
          DO IK=1,IDL1
            C2(IK,L)  = C2(IK,L) +WAR*CH2(IK,J)
            C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
          ENDDO
        ENDDO
      ENDDO
!
      DO J=2,IPPH
        DO IK=1,IDL1
          CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
        ENDDO
      ENDDO
!
      DO J=2,IPPH
        JC = IPP2-J
        DO IK=2,IDL1,2
          CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
          CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
          CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
          CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
        ENDDO
      ENDDO
!
      NAC = 1
      IF(IDO.EQ.2) RETURN
      NAC = 0
      DO IK=1,IDL1
        C2(IK,1) = CH2(IK,1)
      ENDDO
      DO J=2,IP
        DO K=1,L1
          C1(1,K,J) = CH(1,K,J)
          C1(2,K,J) = CH(2,K,J)
        ENDDO
      ENDDO
!
      IF(IDOT.LE.L1) THEN
        IDIJ = 0
        DO J=2,IP
          IDIJ = IDIJ+2
          DO I=4,IDO,2
            IDIJ = IDIJ+2
            DO K=1,L1
              C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
              C1(I,K,J)   = WA(IDIJ-1)*CH(I,K,J)  -WA(IDIJ)*CH(I-1,K,J)
            ENDDO
          ENDDO
        ENDDO
        RETURN
      END IF
!
      IDJ = 2-IDO
      DO J=2,IP
        IDJ = IDJ+IDO
        DO K=1,L1
          IDIJ = IDJ
          DO I=4,IDO,2
            IDIJ = IDIJ+2
            C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
            C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSF2(IDO,L1,CC,CH,WA1)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,2,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,2)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      INTEGER(4)            :: K,I
      REAL(8)               :: TR2,TI2
!     ******************************************************************
      IF(IDO.LE.2) THEN
        DO K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSF3(IDO,L1,CC,CH,WA1,WA2)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,3,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,3)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      REAL(8)   ,INTENT(IN) :: WA2(IDO)
      INTEGER(4)            :: K,I
      REAL(8)               :: TR2,TI2
      REAL(8)               :: CR2,CI2
      REAL(8)               :: CR3,CI3
      REAL(8)               :: DR2,DI2
      REAL(8)               :: DR3,DI3
      REAL(8)   ,PARAMETER  :: TAUR=-0.5D0
      REAL(8)   ,PARAMETER  :: TAUI=-0.866025403784439D0
!     ******************************************************************
      IF(IDO.EQ.2) THEN
        DO  K=1,L1
          TR2 = CC(1,2,K)+CC(1,3,K)
          CR2 = CC(1,1,K)+TAUR*TR2
          CH(1,K,1) = CC(1,1,K)+TR2
          TI2 = CC(2,2,K)+CC(2,3,K)
          CI2 = CC(2,1,K)+TAUR*TI2
          CH(2,K,1) = CC(2,1,K)+TI2
          CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
          CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
          CH(1,K,2) = CR2-CI3
          CH(1,K,3) = CR2+CI3
          CH(2,K,2) = CI2+CR3
          CH(2,K,3) = CI2-CR3
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSF4(IDO,L1,CC,CH,WA1,WA2,WA3)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,4,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,4)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      REAL(8)   ,INTENT(IN) :: WA2(IDO)
      REAL(8)   ,INTENT(IN) :: WA3(IDO)
      INTEGER(4)            :: K,I
      REAL(8)               :: TR1,TI1
      REAL(8)               :: TR2,TI2
      REAL(8)               :: TR3,TI3
      REAL(8)               :: TR4,TI4
      REAL(8)               :: CR2,CI2
      REAL(8)               :: CR3,CI3
      REAL(8)               :: CR4,CI4
!     ******************************************************************
      IF(IDO.EQ.2) THEN
        DO K=1,L1
          TI1 = CC(2,1,K)-CC(2,3,K)
          TI2 = CC(2,1,K)+CC(2,3,K)
          TR4 = CC(2,2,K)-CC(2,4,K)
          TI3 = CC(2,2,K)+CC(2,4,K)
          TR1 = CC(1,1,K)-CC(1,3,K)
          TR2 = CC(1,1,K)+CC(1,3,K)
          TI4 = CC(1,4,K)-CC(1,2,K)
          TR3 = CC(1,2,K)+CC(1,4,K)
          CH(1,K,1) = TR2+TR3
          CH(1,K,3) = TR2-TR3
          CH(2,K,1) = TI2+TI3
          CH(2,K,3) = TI2-TI3
          CH(1,K,2) = TR1+TR4
          CH(1,K,4) = TR1-TR4
          CH(2,K,2) = TI1+TI4
          CH(2,K,4) = TI1-TI4
        ENDDO
        RETURN
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,5,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,5)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      REAL(8)   ,INTENT(IN) :: WA2(IDO)
      REAL(8)   ,INTENT(IN) :: WA3(IDO)
      REAL(8)   ,INTENT(IN) :: WA4(IDO)
      INTEGER(4)            :: K,I
      REAL(8)               :: TR1,TI1
      REAL(8)               :: TR2,TI2
      REAL(8)               :: TR3,TI3
      REAL(8)               :: TR4,TI4
      REAL(8)               :: TR5,TI5
      REAL(8)               :: CR2,CI2
      REAL(8)               :: CR3,CI3
      REAL(8)               :: CR4,CI4
      REAL(8)               :: CR5,CI5
      REAL(8)               :: DR2,DI2
      REAL(8)               :: DR3,DI3
      REAL(8)               :: DR4,DI4
      REAL(8)               :: DR5,DI5
      REAL(8)   ,PARAMETER  :: TR11=+0.309016994374947D0
      REAL(8)   ,PARAMETER  :: TI11=-0.951056516295154D0
      REAL(8)   ,PARAMETER  :: TR12=-0.809016994374947D0
      REAL(8)   ,PARAMETER  :: TI12=-0.587785252292473D0
!     ******************************************************************
      IF(IDO.NE.2) THEN
        DO K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
          ENDDO
        ENDDO
      END IF
      RETURN
      END
#ENDIF
!#ENDIF(IBMLICENSE)










