      SUBROUTINE INTSTP(LQUAD,N, F,FTOL,GTOL, STP0,
     ;                  X,G,S, ALPHA,NFEV, INFO)
!
!     Interpolates a step length
!     http://trilinos.sandia.gov/packages/docs/r5.0/packages/nox/doc/html/classNOX_1_1LineSearch_1_1Polynomial.html
      implicit none
      REAL*8, INTENT(INOUT) :: X(N)
      REAL*8, INTENT(IN) :: G(N), S(N), FTOL, GTOL, F, STP0
      INTEGER*4, INTENT(IN) :: N
      LOGICAL*4, INTENT(IN) :: LQUAD
      INTEGER*4, INTENT(INOUT) :: NFEV, INFO
      REAL*8, INTENT(OUT) :: ALPHA
      REAL*8, ALLOCATABLE :: WORK(:)
      REAL*8 LMAT(2,2), P(2), PHI0, PHI1, PHI2, DPHI0, STP1, STP2, 
     ;       A, B, FTEST1, FTEST2 
      REAL*8 DDOT 
      SAVE  
!
!----------------------------------------------------------------------!
!
      print *, lquad,n, f,ftol,gtol, stp0
      IF (INFO ==-1) GOTO 45 
      IF (.NOT.ALLOCATED(WORK)) ALLOCATE(WORK(N))  
      CALL DCOPY(N,X,1,WORK,1) !save initial search direction
      PHI0 = F
      DPHI0 = DDOT(N,G,1,S,1) 
      IF (DPHI0.GE.0.D0) THEN
         WRITE(*,*) 'intstp: This is not a descent direction!'
         INFO = 7
      ENDIF
      NFEV = NFEV + 1 
      ALPHA = STP0
      INFO =-1 
      X(1:N) = WORK(1:N) + ALPHA*S(1:N)
      RETURN   !calculate new objective function 
   45 CONTINUE !second 
!
!.... quadratic model
      IF (LQUAD) THEN
         IF (NFEV.EQ.1) THEN !interpolate step
            NFEV = NFEV + 1
            PHI1 = F
            STP1 = ALPHA
            ALPHA =-PHI0*STP1**2/(2.D0*(PHI1 - PHI0 - DPHI0*STP0))
            INFO =-1
            X(1:N) = WORK(1:N) + ALPHA*S(1:N)
            RETURN
         ELSE !convergence test
            INFO = 1
            FTEST1 = PHI0 + ALPHA*FTOL*DPHI0 
            FTEST2 = GTOL*DABS(DPHI0)
            IF (F.LE.FTEST1) !sufficient decrease
     ;      WRITE(*,*) 'intstp: Step has sufficient decrease'
            IF (DABS(DDOT(N,G,1,F,1)).LE.FTEST2) !curvature
     ;      WRITE(*,*) 'intstp: Step satisfies curvature condition'
            IF (F.GT.PHI0) THEN
               INFO = 3
               WRITE(*,*) 'intstp: Warning this is a bad step'
               IF (PHI1 < PHI0) THEN
                  INFO = 1
                  WRITE(*,*) 'intstp: Overriding to stp0'
                  ALPHA = STP0
               ENDIF
            ENDIF
            X(1:N) = WORK(1:N) + ALPHA*S(1:N)
         ENDIF 
         IF (ALLOCATED(WORK)) DEALLOCATE(WORK)
         RETURN
      ENDIF
      IF (NFEV.EQ.1) THEN 
         NFEV = NFEV + 1
         PHI1 = F
         STP1 = STP0/2.D0
         ALPHA = STP1 
         INFO =-1
         X(1:N) = WORK(1:N) + ALPHA*S(1:N)
         RETURN
      ELSEIF (NFEV.EQ.2) THEN
         NFEV = NFEV + 1
         PHI2 = F 
         INFO = 1 
!
!....... cubic model
         LMAT(1,1) = (1.D0/STP1)**2
         LMAT(1,2) =-(1.D0/STP0)**2
         LMAT(2,1) =-STP0*(1.D0/STP1)**2
         LMAT(2,2) = STP1*(1.D0/STP0)**2
         P(1) = PHI1 - PHI0 - DPHI0*STP1 
         P(2) = PHI2 - PHI0 - DPHI0*STP2 
         A = 1.D0/(STP1 - STP2)*(LMAT(1,1)*P(1) + LMAT(1,2)*P(2))
         B = 1.D0/(STP1 - STP2)*(LMAT(2,1)*P(1) + LMAT(2,2)*P(2))
         ALPHA =-B + DSQRT(B**2 - 3.D0*A*DPHI0)/(3.D0*A)
         X(1:N) = WORK(1:N) + ALPHA*S(1:N)
         RETURN
      ELSE !convergence test 
         FTEST1 = PHI0 + ALPHA*FTOL*DPHI0
         FTEST2 = GTOL*DABS(DPHI0)
         IF (F.LE.FTEST1) !sufficient decrease
     ;   WRITE(*,*) 'intstp: Cubic step has sufficient decrease'
         IF (DABS(DDOT(N,G,1,F,1)).LE.FTEST2) !curvature
     ;   WRITE(*,*) 'intstp: Cubic step satisfies curvature condition'
         IF (F > PHI0) THEN
            WRITE(*,*) 'intstp: This is a bad step'
            IF (PHI1 < PHI0 .OR. PHI2 < PHI0) THEN
               INFO = 1
               IF (PHI1 < PHI2) THEN
                  WRITE(*,*) 'intstp: Overriding to stp 1'
                  ALPHA = STP1
               ELSE
                  WRITE(*,*) 'intstp: Overriding to stp2'
                  ALPHA = STP0
               ENDIF 
            ENDIF
            X(1:N) = WORK(1:N) + ALPHA*S(1:N)
         ELSE
            INFO = 3
         ENDIF
      ENDIF 
      IF (ALLOCATED(WORK)) DEALLOCATE(WORK)
      RETURN
      END 
C    ------------------------------------------------------------------
C
C     **************************
C     LINE SEARCH ROUTINE MCSRCH
C     **************************
C
      SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,GTOL,XTOL,STPMIN,STPMAX,
     ;                  MAXFEV,INFO,NFEV,WA,LP)
      INTEGER N,MAXFEV,INFO,NFEV
      DOUBLE PRECISION F,STP,FTOL,GTOL,XTOL,STPMIN,STPMAX
      DOUBLE PRECISION X(N),G(N),S(N),WA(N)
!     COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
      SAVE
C
C                     SUBROUTINE MCSRCH
C                
C     A slight modification of the subroutine CSRCH of More' and Thuente.
C     The changes are to allow reverse communication, and do not affect
C     the performance of the routine. 
C
C     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
C     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
C
C     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
C     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
C     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
C     MINIMIZER OF THE MODIFIED FUNCTION
C
C          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
C
C     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
C     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
C     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
C     CONTAINS A MINIMIZER OF F(X+STP*S).
C
C     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
C     THE SUFFICIENT DECREASE CONDITION
C
C           F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),
C
C     AND THE CURVATURE CONDITION
C
C           ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).
C
C     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
C     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
C     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
C     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
C     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
C     SATISFIES THE SUFFICIENT DECREASE CONDITION.
C
C     THE SUBROUTINE STATEMENT IS
C
C        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA)
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF VARIABLES.
C
C       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
C         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
C         X + STP*S.
C
C       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
C         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
C
C       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
C         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
C         OF F AT X + STP*S.
C
C       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
C         SEARCH DIRECTION.
C
C       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
C         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
C         STP CONTAINS THE FINAL ESTIMATE.
C
C       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
C         communication implementation GTOL is defined in a COMMON
C         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
C         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
C         SATISFIED.
C
C       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
C         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C         IS AT MOST XTOL.
C
C       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
C         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
C         communication implementatin they are defined in a COMMON
C         statement).
C
C       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
C         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
C         MAXFEV BY THE END OF AN ITERATION.
C
C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
C
C         INFO = 0  IMPROPER INPUT PARAMETERS.
C
C         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
C
C         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
C                   DIRECTIONAL DERIVATIVE CONDITION HOLD.
C
C         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C                   IS AT MOST XTOL.
C
C         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
C
C         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
C
C         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
C
C         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
C                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
C                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
C                   TOLERANCES MAY BE TOO SMALL.
C
C       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
C         CALLS TO FCN.
C
C       WA IS A WORK ARRAY OF LENGTH N.
C
C     SUBPROGRAMS CALLED
C
C       MCSTEP
C
C       FORTRAN-SUPPLIED...ABS,MAX,MIN
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
C     JORGE J. MORE', DAVID J. THUENTE
C
C     **********
      INTEGER INFOC,J
      LOGICAL BRACKT,STAGE1
      DOUBLE PRECISION DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,
     *       FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY,
     *       STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO
      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,4.0D0,0.0D0/
      IF(INFO.EQ.-1) GO TO 45
      INFOC = 1
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (N .LE. 0) THEN
         WRITE(*,*) 'MCSRCH ERROR: N <= 0'
         RETURN
      ENDIF
      IF (STP .LE. ZERO) THEN
         WRITE(*,*) 'MCSRCH ERROR: STP <= 0'
         RETURN
      ENDIF
      IF (FTOL .LT. ZERO) THEN
         WRITE(*,*) 'MCSRCH ERROR: FTOL < 0'
         RETURN
      ENDIF 
      IF (GTOL .LT. ZERO) THEN
         WRITE(*,*) 'MCSRCH ERROR: GTOL < 0'
         RETURN
      ENDIF
      IF (XTOL .LT. ZERO) THEN
         WRITE(*,*) 'MCSRCH ERROR: XTOL < 0'
         RETURN
      ENDIF
      IF (STPMIN .LT. ZERO) THEN
         WRITE(*,*) 'MCSRCH ERROR: STPMIN < 0'
         RETURN
      ENDIF
      IF (STPMAX .LT. STPMIN) THEN
         WRITE(*,*) 'MCSRCH STPMAX < STPMIN'
         RETURN
      ENDIF
      IF (MAXFEV .LE. 0) THEN
         WRITE(*,*) 'MCSRCH MAXFEV <= 0'
         RETURN
      ENDIF
C
C     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
C     AND CHECK THAT S IS A DESCENT DIRECTION.
C
      DGINIT = ZERO
      DO 10 J = 1, N
         DGINIT = DGINIT + G(J)*S(J)
   10    CONTINUE
      IF (DGINIT .GE. ZERO) then
         write(LP,15)
   15    FORMAT(/'  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION')
         RETURN
         ENDIF
C
C     INITIALIZE LOCAL VARIABLES.
C
      BRACKT = .FALSE.
      STAGE1 = .TRUE.
      NFEV = 0
      FINIT = F
      DGTEST = FTOL*DGINIT
      WIDTH = STPMAX - STPMIN
      WIDTH1 = WIDTH/P5
      DO 20 J = 1, N
         WA(J) = X(J)
   20    CONTINUE
C
C     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
C     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
C     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
C     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
C     THE INTERVAL OF UNCERTAINTY.
C     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
C     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
C
      STX = ZERO
      FX = FINIT
      DGX = DGINIT
      STY = ZERO
      FY = FINIT
      DGY = DGINIT
C
C     START OF ITERATION.
C
   30 CONTINUE
C
C        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
C        TO THE PRESENT INTERVAL OF UNCERTAINTY.
C
         IF (BRACKT) THEN
            STMIN = MIN(STX,STY)
            STMAX = MAX(STX,STY)
         ELSE
            STMIN = STX
            STMAX = STP + XTRAPF*(STP - STX)
            END IF
C
C        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
C
         STP = MAX(STP,STPMIN)
         STP = MIN(STP,STPMAX)
C
C        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
C        STP BE THE LOWEST POINT OBTAINED SO FAR.
C
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     *      .OR. NFEV .GE. MAXFEV-1 .OR. INFOC .EQ. 0
     *      .OR. (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)) STP = STX
C
C        EVALUATE THE FUNCTION AND GRADIENT AT STP
C        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
C        We return to main program to obtain F and G.
C
         DO 40 J = 1, N
            X(J) = WA(J) + STP*S(J)
   40       CONTINUE
         INFO=-1
         RETURN
C
   45    INFO=0
         NFEV = NFEV + 1
         DG = ZERO
         DO 50 J = 1, N
            DG = DG + G(J)*S(J)
   50       CONTINUE
         FTEST1 = FINIT + STP*DGTEST
C
C        TEST FOR CONVERGENCE.
C
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     *      .OR. INFOC .EQ. 0) INFO = 6
         IF (STP .EQ. STPMAX .AND.
     *       F .LE. FTEST1 .AND. DG .LE. DGTEST) INFO = 5
         IF (STP .EQ. STPMIN .AND.
     *       (F .GT. FTEST1 .OR. DG .GE. DGTEST)) INFO = 4
         IF (NFEV .GE. MAXFEV) INFO = 3
         IF (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX) INFO = 2
         IF (F .LE. FTEST1 .AND. ABS(DG) .LE. GTOL*(-DGINIT)) INFO = 1
C
C        CHECK FOR TERMINATION.
C
         IF (INFO .NE. 0) RETURN
C
C        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
C        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
C
         IF (STAGE1 .AND. F .LE. FTEST1 .AND.
     *       DG .GE. MIN(FTOL,GTOL)*DGINIT) STAGE1 = .FALSE.
C
C        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
C        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
C        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
C        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
C        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
C
         IF (STAGE1 .AND. F .LE. FX .AND. F .GT. FTEST1) THEN
C
C           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
C
            FM = F - STP*DGTEST
            FXM = FX - STX*DGTEST
            FYM = FY - STY*DGTEST
            DGM = DG - DGTEST
            DGXM = DGX - DGTEST
            DGYM = DGY - DGTEST
C
C           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
C           AND TO COMPUTE THE NEW STEP.
C
            CALL MCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,
     *                 BRACKT,STMIN,STMAX,INFOC)
C
C           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
C
            FX = FXM + STX*DGTEST
            FY = FYM + STY*DGTEST
            DGX = DGXM + DGTEST
            DGY = DGYM + DGTEST
         ELSE
C
C           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
C           AND TO COMPUTE THE NEW STEP.
C
            CALL MCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG,
     *                 BRACKT,STMIN,STMAX,INFOC)
            END IF
C
C        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
C        INTERVAL OF UNCERTAINTY.
C
         IF (BRACKT) THEN
            IF (ABS(STY-STX) .GE. P66*WIDTH1)
     *         STP = STX + P5*(STY - STX)
            WIDTH1 = WIDTH
            WIDTH = ABS(STY-STX)
            END IF
C
C        END OF ITERATION.
C
         GO TO 30
C
C     LAST LINE OF SUBROUTINE MCSRCH.
C
      END
      SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
     *                 STPMIN,STPMAX,INFO)
      INTEGER INFO
      DOUBLE PRECISION STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX
      LOGICAL BRACKT,BOUND
C
C     SUBROUTINE MCSTEP
C
C     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
C     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
C     A MINIMIZER OF THE FUNCTION.
C
C     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
C     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
C     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
C     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
C     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
C     WITH ENDPOINTS STX AND STY.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
C                        STPMIN,STPMAX,INFO)
C
C     WHERE
C
C       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
C         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
C         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
C         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
C         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
C
C       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
C         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
C         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
C         UPDATED APPROPRIATELY.
C
C       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,
C         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
C         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
C         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
C
C       BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER
C         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
C         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
C         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
C
C       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
C         AND UPPER BOUNDS FOR THE STEP.
C
C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
C         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
C         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
C         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
C     JORGE J. MORE', DAVID J. THUENTE
C
      DOUBLE PRECISION GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
      INFO = 0
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF ((BRACKT .AND. (STP .LE. MIN(STX,STY) .OR.
     *    STP .GE. MAX(STX,STY))) .OR.
     *    DX*(STP-STX) .GE. 0.D0 .OR. STPMAX .LT. STPMIN) RETURN
C
C     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
C
      SGND = DP*(DX/ABS(DX))
C
C     FIRST CASE. A HIGHER FUNCTION VALUE.
C     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
C     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
C     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
C
      IF (FP .GT. FX) THEN
         INFO = 1
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .LT. STX) GAMMA = -GAMMA
         P = (GAMMA - DX) + THETA
         Q = ((GAMMA - DX) + GAMMA) + DP
         R = P/Q
         STPC = STX + R*(STP - STX)
         STPQ = STX + ((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP - STX)
         IF (ABS(STPC-STX) .LT. ABS(STPQ-STX)) THEN
            STPF = STPC
         ELSE
           STPF = STPC + (STPQ - STPC)/2
           END IF
         BRACKT = .TRUE.
C
C     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
C     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
C     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
C     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
C
      ELSE IF (SGND .LT. 0.D0) THEN
         INFO = 2
         BOUND = .FALSE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = ((GAMMA - DP) + GAMMA) + DX
         R = P/Q
         STPC = STP + R*(STX - STP)
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF (ABS(STPC-STP) .GT. ABS(STPQ-STP)) THEN
            STPF = STPC
         ELSE
            STPF = STPQ
            END IF
         BRACKT = .TRUE.
C
C     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
C     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
C     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
C     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
C     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
C     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
C     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
C     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
C
      ELSE IF (ABS(DP) .LT. ABS(DX)) THEN
         INFO = 3
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
C
C        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
C        TO INFINITY IN THE DIRECTION OF THE STEP.
C
         GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = (GAMMA + (DX - DP)) + GAMMA
         R = P/Q
         IF (R .LT. 0.0 .AND. GAMMA .NE. 0.0) THEN
            STPC = STP + R*(STX - STP)
         ELSE IF (STP .GT. STX) THEN
            STPC = STPMAX
         ELSE
            STPC = STPMIN
            END IF
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF (BRACKT) THEN
            IF (ABS(STP-STPC) .LT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
               END IF
         ELSE
            IF (ABS(STP-STPC) .GT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
               END IF
            END IF
C
C     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
C     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
C     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
C     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
C
      ELSE
         INFO = 4
         BOUND = .FALSE.
         IF (BRACKT) THEN
            THETA = 3*(FP - FY)/(STY - STP) + DY + DP
            S = MAX(ABS(THETA),ABS(DY),ABS(DP))
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))
            IF (STP .GT. STY) GAMMA = -GAMMA
            P = (GAMMA - DP) + THETA
            Q = ((GAMMA - DP) + GAMMA) + DY
            R = P/Q
            STPC = STP + R*(STY - STP)
            STPF = STPC
         ELSE IF (STP .GT. STX) THEN
            STPF = STPMAX
         ELSE
            STPF = STPMIN
            END IF
         END IF
C
C     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
C     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
C
      IF (FP .GT. FX) THEN
         STY = STP
         FY = FP
         DY = DP
      ELSE
         IF (SGND .LT. 0.0) THEN
            STY = STX
            FY = FX
            DY = DX
            END IF
         STX = STP
         FX = FP
         DX = DP
         END IF
C
C     COMPUTE THE NEW STEP AND SAFEGUARD IT.
C
      STPF = MIN(STPMAX,STPF)
      STPF = MAX(STPMIN,STPF)
      STP = STPF
      IF (BRACKT .AND. BOUND) THEN
         IF (STY .GT. STX) THEN
            STP = MIN(STX+0.66D0*(STY-STX),STP)
         ELSE
            STP = MAX(STX+0.66D0*(STY-STX),STP)
            END IF
         END IF
      RETURN
C
C     LAST LINE OF SUBROUTINE MCSTEP.
C
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE MCSRCH_ERROR(INFO,MAXFEV,MYDEST,IERR)
!
!     Checks for and writes errors from mcsrch
      IMPLICIT NONE
      INTEGER*4, INTENT(IN) :: INFO, MAXFEV
      INTEGER*4, INTENT(OUT) :: MYDEST, IERR 
!
!----------------------------------------------------------------------!
!
      IERR = 0
      MYDEST = 2000
      IF (INFO == 0) THEN !error messages
         WRITE(*,*) 'mcsrch_error: Improper input parameters'
         IERR = 1  
         !GOTO 500
         MYDEST = 500
      ENDIF
      IF (INFO == 2) THEN 
         WRITE(*,*) 'mcsrch_error: Relative width of uncertainty too sma
     ;ll'
         IERR = 1
         !GOTO 500
         MYDEST = 500
      ENDIF
      IF (INFO == 3) THEN 
         WRITE(*,*) 'mcsrch_error: More than MAXFEV function evaluations
     ; required', MAXFEV
         IERR = 1
         !GOTO 500
         MYDEST = 500
      ENDIF  
      IF (INFO == 4) THEN 
         WRITE(*,*) 'mcsrch_error: The step size is too small'
         IERR = 1
         !GOTO 500
         MYDEST = 500
       ENDIF   
       IF (INFO == 5) THEN 
          WRITE(*,*) 'mcsrch_error: The step size is too large'
          IERR = 1
          !GOTO 500
          MYDEST = 500
       ENDIF
       IF (INFO == 6) THEN 
          WRITE(*,*) 'mcsrch_error: Rounding errors prevent further prog
     ;ress'
          IERR = 1
          !GOTO 500
          MYDEST = 500
       ENDIF
!
!..... line search iteration
       IF (INFO ==-1) MYDEST = 3000
       RETURN
       END




