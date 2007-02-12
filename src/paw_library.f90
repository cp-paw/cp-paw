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
!     .....................................................................
      SUBROUTINE LIB$GENERALEIGENVALUER8(N,H,S,E,VEC)
!     **                                                                 **
!     ** SOLVES THE GENERALIZED, REAL NON-SYMMETRIC EIGENVALUE PROBLEM   **
!     **      [H(:,:)-E(I)*S(:,:)]*VEC(:,I)=0                            **
!     ** NOT TESTED!!!!!!                                                **
!     **                                                                 **
!     ** REMARK: H AND S MUST BE SYMMETRIC                              **
!     **         S MUST BE POSITIVE DEFINITE                             **
!     **         EIGENVECTORS ARE ORTHONORMAL IN THE SENSE               **
!     **             MATMUL(TRANSPOSE(VEC),MATMUL(S,VEC))=IDENTITY       **
!     **                                                                 **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N
      REAL(8)   ,INTENT(IN)  :: H(N,N)    ! HAMILTON MATRIX
      REAL(8)   ,INTENT(IN)  :: S(N,N)    ! OVERLAP MATRIX
      REAL(8)   ,INTENT(OUT) :: E(N)      ! EIGENVALUES
      REAL(8)   ,INTENT(OUT) :: VEC(N,N)  ! EIGENVECTORS
      INTEGER                :: LWORK     ! SIZE OF WORK ARRAY
      REAL(8)                :: WORK(16*N)! WORK ARRAY FOR LAPACK ROUTINE
      INTEGER                :: INFO      ! ERROR CODE OF LAPACK ROUTINE
      REAL(8)                :: B(N,N)    ! COPY OF OVERLAP MATRIX
      LOGICAL  ,PARAMETER    :: TTEST=.TRUE. ! IF TRUE TEST RESULT
      REAL(8)                :: DEV       ! DEVIATION
      INTEGER                :: I 
!     *********************************************************************
!
!     ========================================================================
!     == TAKE CARE OF TRIVIAL CASES                                         ==
!     ========================================================================
      IF(N.EQ.1) THEN
        E(1)=H(1,1)/S(1,1)
        VEC(1,1)=1.D0
        RETURN
      ELSE IF(N.EQ.0) THEN
        RETURN
      END IF
!
!     ========================================================================
!     == TEST IF INPUT MATRICES ARE SYMMETRIC                               ==
!     ========================================================================
      IF(TTEST) THEN
        DEV=SUM(ABS(H-TRANSPOSE(H)))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('HAMILTON MATRIX NOT SYMMETRIC')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
        END IF
        DEV=SUM(ABS(S-TRANSPOSE(S)))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('OVERLAP MATRIX NOT SYMMETRIC')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
        END IF
      END IF
!
!     ========================================================================
!     == CALL LAPACK ROUTINE                                                ==
!     ========================================================================
      LWORK=16*N     
      VEC=H          ! LAPACK ROUTINE OVERWRITES HAMILTONIAN WITH EIGENVECTORS
      B=S            ! LAPACK ROUTINE OVERWRITES INPUT
      CALL DSYGV(1,'V','U',N,VEC,N,B,N,E,WORK,LWORK,INFO)
      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('ITH ARGUMENT OF DSYGV HAS ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$MSG('THE QZ ITERATION FAILED')
        CALL ERROR$I4VAL('INFO',INFO)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
      END IF
!
!     ========================================================================
!     == TEST RESULT OF THE ROUTINE                                         ==
!     ========================================================================
      IF(TTEST) THEN
        B=MATMUL(TRANSPOSE(VEC),MATMUL(S,VEC))
        DO I=1,N
          B(I,I)=B(I,I)-1.D0
          DEV=SUM(ABS(MATMUL(H-E(I)*S,VEC(:,I))))
          IF(DEV.GT.1.D-7) THEN
            CALL ERROR$MSG('EIGENVALUE PROBLEM FAILED')
            CALL ERROR$R8VAL('DEV',DEV)
            CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
          END IF
        ENDDO
!       == THIS IS NOT FULFILLED FOR TSYMMETRIC=.FALSE.
        DEV=SUM(ABS(B))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('EIGENSTATES NOT ORTHONORMAL')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
        END IF
      END IF
      RETURN
      END SUBROUTINE LIB$GENERALEIGENVALUER8
!
!     ......................................................................
      SUBROUTINE LIB$GENERALEIGENVALUEC8(N,H,S,E,VEC)
!     **                                                                 **
!     ** SOLVES THE GENERALIZED, REAL NON-SYMMETRIC EIGENVALUE PROBLEM   **
!     **      [H(:,:)-E(I)*S(:,:)]*VEC(:,I)=0                            **
!     **                                                                 **
!     ** REMARK: H AND S MUST BE HERMITEANC                              **
!     **         S MUST BE POSITIVE DEFINITE                             **
!     **         EIGENVECTORS ARE ORTHONORMAL IN THE SENSE               **
!     **             MATMUL(TRANSPOSE(VEC),MATMUL(S,VEC))=IDENTITY       **
!     **                                                                 **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)    ! HAMITON MATRIX
      COMPLEX(8),INTENT(IN) :: S(N,N)    ! OVERLAP MATRIX
      REAL(8)   ,INTENT(OUT):: E(N)      ! EIGENVALUES
      COMPLEX(8),INTENT(OUT):: VEC(N,N)  ! EIGENVECTORS
      INTEGER               :: LDWORK
      COMPLEX(8)            :: WORK(N*N)
      COMPLEX(8)            :: S1(N,N)
      REAL(8)               :: RWORK(3*N-2)
      INTEGER               :: INFO
      LOGICAL   ,PARAMETER  :: TTEST=.TRUE.
      REAL(8)               :: DEV
      INTEGER               :: I
!     *********************************************************************
!
!     ========================================================================
!     == TEST IF INPUT MATRICES ARE SYMMETRIC                               ==
!     ========================================================================
      IF(TTEST) THEN
        DEV=SUM(ABS(H-TRANSPOSE(CONJG(H))))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('HAMILTON MATRIX NOT HERMITEAN')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
        END IF
        DEV=SUM(ABS(S-TRANSPOSE(CONJG(S))))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('OVERLAP MATRIX NOT HERMITEAN')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
        END IF
      END IF
!
!     ========================================================================
!     == CALL LAPACK ROUTINE                                                ==
!     ========================================================================
      LDWORK=N*N
      VEC=H     ! LAPACK ROUTINE OVERWRITES HAMILTONIAN WITH EIGENVECTORS
      S1=S
      CALL ZHEGV(1,'V','U',N,VEC,N,S1,N,E,WORK,LDWORK,RWORK,INFO)
      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('ITH ARGUMENT OF ZHGEV HAS ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$MSG('FAILED')
        CALL ERROR$I4VAL('INFO',INFO)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
      END IF
!
!     ========================================================================
!     == TEST RESULT OF THE ROUTINE                                         ==
!     ========================================================================
      IF(TTEST) THEN
        S1=MATMUL(TRANSPOSE(CONJG(VEC)),MATMUL(S,VEC))
        DEV=0.D0
        DO I=1,N
          S1(I,I)=S1(I,I)-(1.D0,0.D0)
          DEV=MAX(DEV,MAXVAL(ABS(MATMUL(H-E(I)*S,VEC(:,I)))))
        ENDDO
        IF(DEV.GT.1.D-6) THEN
          CALL ERROR$MSG('GENERAL EIGENVALUE TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
        END IF
        DEV=SUM(ABS(S1))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('EIGENSTATES NOT ORTHONORMAL')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
        END IF
      END IF
      RETURN
      END
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
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      CHARACTER(4),SAVE       :: DIRSAVE=''
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      INTEGER     ,SAVE       :: ISIGN
      REAL(8)     ,SAVE       :: SCALE
      COMPLEX(8)              :: XDUMMY(LEN,NFFT)
      INTEGER(4)              :: I
      INTEGER(4),SAVE         :: NP=0
      INTEGER(4),PARAMETER    :: NPX=10 ! #(DIFFERENT FFT PLANS)
      INTEGER(8),SAVE         :: PLANS(NPX,2),PLAN=-1
      LOGICAL                 :: DEF
!     INCLUDE 'FFTW_F77.I'
!     ***********  fftw_f77.i *******************************************
!     This file contains PARAMETER statements for various constants
!     that can be passed to FFTW routines.  You should include
!     this file in any FORTRAN program that calls the fftw_f77
!     routines (either directly or with an #include statement
!     if you use the C preprocessor).
      integer,parameter :: FFTW_FORWARD=-1 ! sign in the exponent of the forward ft
      integer,parameter :: FFTW_BACKWARD=1 ! sign in the exponent of the backward ft
      integer,parameter :: FFTW_REAL_TO_COMPLEX=-1
      integer,parameter :: FFTW_COMPLEX_TO_REAL=1
      integer,parameter :: FFTW_ESTIMATE=0
      integer,parameter :: FFTW_MEASURE=1
      integer,parameter :: FFTW_OUT_OF_PLACE=0
      integer,parameter :: FFTW_IN_PLACE=8
      integer,parameter :: FFTW_USE_WISDOM=16
      integer,parameter :: FFTW_THREADSAFE=128
!     Constants for the MPI wrappers:
      integer,parameter :: FFTW_TRANSPOSED_ORDER=1
      integer,parameter :: FFTW_NORMAL_ORDER=0
      integer,parameter :: FFTW_SCRAMBLED_INPUT=8192
      integer,parameter :: FFTW_SCRAMBLED_OUTPUT=16384
!     ******************************************************************
!
!     ==================================================================
!     ==  initialize fft                                              ==
!     ==================================================================
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
!         CALL FFTW_F77_CREATE_PLAN(PLAN,LEN,fftw_forward,FFTW_MEASURE)
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
      XDUMMY=X
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
      SUBROUTINE LIB$MATRIXSOLVENEW(N,M,NEQ,A,X,B)
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
      REAL(8)               :: A1(N,M)
      REAL(8)               :: B1(N,NEQ)
      INTEGER               :: INFO
      REAL(8)               :: RCOND=-1.D0
      INTEGER               :: LDWORK
      INTEGER               :: IRANK
      REAL(8)               :: SING(N)
      REAL(8)   ,ALLOCATABLE:: WORK(:)
      INTEGER               :: N1,M1,NEQ1
!     ******************************************************************
      IF(N.EQ.1.AND.M.EQ.1) THEN
        X(1,:)=B(1,:)/A(1,1)
        RETURN
      END IF
      LDWORK=3*MIN(M,N)+MAX(2*MIN(M,N),MAX(M,N),NEQ)
!     -- USE 3*M+3*N+NEQ
      ALLOCATE(WORK(LDWORK))
      N1=N
      M1=M
      NEQ1=NEQ
      A1=A 
      B1=B
      CALL DGELSS(N1,M1,NEQ1,A1,N1,B1,N1,SING,RCOND,IRANK,WORK,LDWORK,INFO)
      X=B1
      DEALLOCATE(WORK)
      IF(INFO.GT.0) THEN
        CALL ERROR$MSG('ITH ARGUMENT HAS ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',INFO)
        CALL ERROR$STOP('LIB$MATRIXSOLVENEW')
      ELSE IF(INFO.LT.0) THEN
        CALL ERROR$MSG('SINGLAR VALUE DECOMPOSITION NOT CONVERGED')
        CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF INTERMEDIATE')
        CALL ERROR$MSG('BIDIAGONAL FORM DID NOT CONVERGE TO ZERO')
        CALL ERROR$I4VAL('I',INFO)
        CALL ERROR$STOP('LIB$MATRIXSOLVENEW')
      END IF
      RETURN
      END
!     ..................................................................
      SUBROUTINE LIB$MATRIXSOLVENEWC8(N,M,NEQ,A,X,B)
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
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(OUT):: X(M,NEQ)
      COMPLEX(8),INTENT(IN) :: B(N,NEQ)
      COMPLEX(8)            :: A1(N,M)
      COMPLEX(8)            :: B1(N,NEQ)
      INTEGER               :: INFO
      REAL(8)               :: RCOND=-1.D0
      INTEGER               :: LDWORK
      INTEGER               :: IRANK
      INTEGER               :: IPIVOT(N)
      REAL(8)               :: SING(N)
      REAL(8)   ,ALLOCATABLE:: WORK(:)
      INTEGER               :: N1,M1,NEQ1
!     ******************************************************************
      IF(N.EQ.1.AND.M.EQ.1) THEN
        X(1,:)=B(1,:)/A(1,1)
        RETURN
      END IF
      LDWORK=3*MIN(M,N)+MAX(2*MIN(M,N),MAX(M,N),NEQ)
!     -- USE 3*M+3*N+NEQ
      ALLOCATE(WORK(LDWORK))
      N1=N
      M1=M
      NEQ1=NEQ
      A1=A 
      B1=B
      IF(N1.EQ.M1) THEN
        CALL ZGESV(N1,NEQ1,A1,N1,IPIVOT,B1,N1,INFO)
        X=B1
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('ITH ARGUMENT HAS ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB$MATRIXSOLVENEWC8')
        ELSE IF(INFO.GT.0) THEN
          CALL ERROR$MSG('PROBLEM IS SINGULAR. NO SOLUTION CAN BE COMPUTED')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB$MATRIXSOLVENEW')
        END IF
      ELSE
        CALL ERROR$MSG('N.NEQ.M NOT IMPLEMENTED')
        CALL ERROR$STOP('LIB$MATRIXSOLVENEWC8')
!        CALL ZGELSD(N,M,NEQ,A,N,B,N,S,RCOND,RANK,WORK,LWORK,RWORK,IWORK,INFO)
      END IF
      DEALLOCATE(WORK)
      RETURN
      END
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
      LOGICAL(4)            :: TTEST=.FALSE.
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
        CALL DGESV(N,NEQ,AFACT,N,IPVT,B,N,INFO)
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT TO DGESV HAS ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB$MATRIXSOLVE')
        ELSE IF(INFO.GT.0) THEN
          CALL ERROR$MSG('PROBLEM FOR DGESV IS SINGULAR')
          CALL ERROR$MSG('U(I,I)=0, WHERE A=P*L*U')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB$MATRIXSOLVE')
        END IF          
!!$        == OLD STUFF SHALL BE REMOVED
!!$        CALL DGEFA(AFACT,N,N,IPVT,INFO)
!!$        X=B
!!$        CALL DGESL(AFACT,N,N,IPVT,X,0)
#ENDIF
        DEALLOCATE(AFACT)
        DEALLOCATE(IPVT)
      ELSE   !SINGULAR VALUE DECOMPOSITION
#IF DEFINED(CPPVAR_BLAS_ESSL)
!       ===========================================================
!       == SINGULAR VALUE DECOMPOSION WITH ESSL                  ==
!       ===========================================================
        NM=MAX(1,MAX(N,M))
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
!       ===========================================================
!       == SINGULAR VALUE DECOMPOSION WITH LAPACK DRIVER ROUTINE ==
!       ===========================================================
        NM=MAX(1,MAX(N,M))
        ALLOCATE(AFACT(N,M))
        AFACT(:,:)=A(:,:)
        ALLOCATE(BFACT(NM,NEQ))
        BFACT(1:N,:)=B(:,:)
        ALLOCATE(S(M))    !SINGULAR VALUES
        NAUX=3*NM+MAX(2*MIN(M,N),NM,NEQ)
        ALLOCATE(AUX(NAUX))
!       == CALL LAPACK
        CALL DGELSS(N,M,NEQ,AFACT,NM,BFACT,N,S,TAU,IRANK,AUX,NAUX,INFO) 
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT TO DGELSY HAS ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB$MATRIXSOLVE')
        ELSE IF(INFO.LT.0) THEN
          CALL ERROR$MSG('THE ALGORITHM FOR COMPUTING THE')
          CALL ERROR$MSG('SINGULAR VALUE DECOMPOSITION FAILED TO CONVERGE')
          CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
          CALL ERROR$MSG('BIDIAGONAL FORM DID NOT CONVERGE TO ZERO.')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB$MATRIXSOLVE')
        END IF
        DEALLOCATE(AUX)
        X(:,:)=BFACT(1:M,:)
        DEALLOCATE(S)
        DEALLOCATE(BFACT)
        DEALLOCATE(AFACT)
#ENDIF
      END IF
!     =================================================================
!     ==  TEST                                                       ==
!     =================================================================
      IF(TTEST) THEN
        ALLOCATE(AUX(1))
        AUX(:)=MAXVAL(ABS(MATMUL(A,X)-B))
        WRITE(*,*)'MAX ERROR OF LIB$MATRIXSOLVE ',AUX(1)
        DEALLOCATE(AUX)
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
