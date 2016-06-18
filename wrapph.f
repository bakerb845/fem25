      REAL*4 FUNCTION WRAPPH(IDEG,ANGLE)
!
!     Wraps the phase to be between [-pi,pi)
!     http://stackoverflow.com/questions/4633177/
!     c-how-to-wrap-a-float-to-the-interval-pi-pi
!
!     INPUT     MEANING
!     -----     -------
!     ANGLE     angle to wrap
!     IDEG      =1 input is is in degrees
!
!.... variable declarations
      IMPLICIT NONE
      REAL*4, INTENT(IN) :: ANGLE
      INTEGER*4, INTENT(IN) :: IDEG 
      REAL*4 FANG, PI, TWOPI, PI180, XSIGN 
      PARAMETER(PI = 3.141592653589793E0)
      PARAMETER(TWOPI = 6.283185307179586E0) 
      PARAMETER(PI180 = 0.017453292519943295E0) 
!
!----------------------------------------------------------------------!
!
      IF (IDEG.EQ.1) THEN !degrees -> radians
         FANG = ANGLE*PI180 
      ELSE !already in radians
         FANG = ANGLE
      ENDIF 
      IF (ABS(FANG).LE.PI) THEN
         WRAPPH = FANG 
         RETURN
      ELSE 
         IF (FANG + PI.LT.0.E0) THEN
            XSIGN =-1.E0
         ELSE
            XSIGN = 1.E0
         ENDIF 
         WRAPPH = MOD(ABS(FANG + PI),TWOPI) - PI   
         WRAPPH = XSIGN*WRAPPH
      ENDIF
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*8 FUNCTION DWRAPPH(IDEG,ANGLE)
!
!     Wraps the phase to be between [-pi,pi)
!     http://stackoverflow.com/questions/4633177/
!     c-how-to-wrap-a-float-to-the-interval-pi-pi
!
!     INPUT     MEANING
!     -----     -------
!     ANGLE     angle to wrap
!     IDEG      =1 input is is in degrees
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: ANGLE
      INTEGER*4, INTENT(IN) :: IDEG 
      REAL*8 FANG, PI, TWOPI, PI180, XSIGN 
      PARAMETER(PI = 3.141592653589793D0)
      PARAMETER(TWOPI = 6.283185307179586D0) 
      PARAMETER(PI180 = 0.017453292519943295D0) 
!
!----------------------------------------------------------------------!
!
      IF (IDEG.EQ.1) THEN !degrees -> radians
         FANG = ANGLE*PI180 
      ELSE !already in radians
         FANG = ANGLE
      ENDIF 
      IF (DABS(FANG).LE.PI) THEN
         DWRAPPH = FANG 
         RETURN
      ELSE 
         IF (FANG + PI.LT.0.D0) THEN
            XSIGN =-1.D0
         ELSE
            XSIGN = 1.D0
         ENDIF 
         DWRAPPH = MOD(DABS(FANG + PI),TWOPI) - PI   
         DWRAPPH = XSIGN*DWRAPPH
      ENDIF
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DPUNWR(X,LX,IRD)
! 
!     SIMPLE PHASE UNWRAPPING ROUTINE, USEFUL WHERE PHASE WAS COMPUTED 
!       USING ATAN2 FUNCTION AND IN SIMILAR SITUATIONS WHERE PHASE 
!       IS WRAPPED INTO THE RANGE FROM -PI TO PI RADIANS
!     X(0:LX) = SEQUENCE OF PHASE ANGLES IN RADIANS OR DEGREES, ALL
!               ASSUMED TO BE OUTPUTS OF ATAN2 IN RANGE (-PI,PI) RAD.
!     ITRD = 1 TO INDICATE X IS IN RADIANS, OR 2 TO INDICATE DEGREES
!     THE ROUTINE INSERS A CORRECTION OF 2*PI OR 360 DGREES WHEREVER 
!       THE PHASE JUMPS MORE THAN PI RADIANS
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8, INTENT(INOUT) :: X(0:LX)
      INTEGER*4, INTENT(IN) :: LX, IRD
!
!----------------------------------------------------------------------!
!
      ANGL = 180.D0
      IF (IRD.EQ.1) ANGL = 4.D0*DATAN(1.D0)
      COR = 0.D0 
      DO 2 K=1,LX
         DX = X(K) - (X(K-1) - COR) 
         IF (DABS(DX).LE.ANGL) GOTO 1
         COR = COR - DSIGN(2.D0*ANGL,DX)
    1    X(K) = X(K) + COR
    2 CONTINUE
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SPUNWR(X,LX,IRD)
! 
!     SIMPLE PHASE UNWRAPPING ROUTINE, USEFUL WHERE PHASE WAS COMPUTED 
!       USING ATAN2 FUNCTION AND IN SIMILAR SITUATIONS WHERE PHASE 
!       IS WRAPPED INTO THE RANGE FROM -PI TO PI RADIANS
!     X(0:LX) = SEQUENCE OF PHASE ANGLES IN RADIANS OR DEGREES, ALL
!               ASSUMED TO BE OUTPUTS OF ATAN2 IN RANGE (-PI,PI) RAD.
!     ITRD = 1 TO INDICATE X IS IN RADIANS, OR 2 TO INDICATE DEGREES
!     THE ROUTINE INSERS A CORRECTION OF 2*PI OR 360 DGREES WHEREVER 
!       THE PHASE JUMPS MORE THAN PI RADIANS
      REAL*4, INTENT(INOUT) :: X(0:LX)
      INTEGER*4, INTENT(IN) :: LX, IRD 
!
!----------------------------------------------------------------------!
!
      ANGL = 180.0
      IF (IRD.EQ.1) ANGL = 4.0*ATAN(1.0)
      COR = 0.0 
      DO 2 K=1,LX
         DX = X(K) - (X(K-1) - COR) 
         IF (ABS(DX).LE.ANGL) GOTO 1
         COR = COR - SIGN(2.0*ANGL,DX)
    1    X(K) = X(K) + COR 
    2 CONTINUE
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*8 FUNCTION DPHASE(Q)
!     Calculates the phase angle (double precision)
      COMPLEX*16, INTENT(IN) :: Q
      DPHASE = DATAN2(DIMAG(Q),DREAL(Q))
      RETURN
      END
      REAL*4 FUNCTION SPHASE(Q)
!     Calculates the phase angle (single precision)
      COMPLEX*8, INTENT(IN) :: Q
      SPHASE = ATAN2(IMAG(Q),REAL(Q))
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      COMPLEX*16 FUNCTION ZPHM2CM(XMAG,PHASE) 
!     Converts phase, magnitude to a complex number (double precision)
      REAL*8, INTENT(IN) :: XMAG, PHASE
      ZPHM2CM = DCMPLX(XMAG,0.D0)*CDEXP(DCMPLX(0.D0,PHASE)) 
      RETURN
      END
      COMPLEX*8 FUNCTION CPHM2CM(XMAG,PHASE)
!     Converts phase, magnitude to a complex number (single precision)
      REAL*4, INTENT(IN) :: XMAG, PHASE
      CPHM2CM = CMPLX(XMAG,0.D0)*CEXP(CMPLX(0.0,PHASE))
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CROTATE(THETA,Q1,Q2, RQ1,RQ2)
!     If theta > 0 this is a counter-clockwise rotation
      COMPLEX*8, INTENT(IN) :: Q1, Q2 
      REAL*4, INTENT(IN) :: THETA 
      COMPLEX*8, INTENT(OUT) :: RQ1, RQ2
      COMPLEX*8 COST, SINT
      REAL*4 PI180 
      PARAMETER(PI180 = 0.017453292519943295)
      COST = CMPLX(COS(THETA*PI180),0.0)
      SINT = CMPLX(SIN(THETA*PI180),0.0) 
      RQ1 = Q1*COST - Q2*SINT
      RQ2 = Q1*SINT + Q2*COST
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE ZROTATE(THETA,Q1,Q2, RQ1,RQ2)
!     If theta > 0 this is a counter-clockwise rotation
      COMPLEX*16, INTENT(IN) :: Q1, Q2  
      REAL*8, INTENT(IN) :: THETA 
      COMPLEX*16, INTENT(OUT) :: RQ1, RQ2
      COMPLEX*16 COST, SINT
      REAL*8 PI180
      PARAMETER(PI180 = 0.017453292519943295D0)
      COST = DCMPLX(DCOS(THETA*PI180),0.D0)
      SINT = DCMPLX(DSIN(THETA*PI180),0.D0) 
      RQ1 = Q1*COST - Q2*SINT
      RQ2 = Q1*SINT + Q2*COST
      RETURN
      END
 
