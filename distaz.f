      subroutine bjdaz2(zlat0,alon0,zlat,alon,delt,az1,az2,lattype)
!
!     Calculate geocentric distance, azimuth, and back azimuth from
!     a reference point to another point on Earth's surface
!
!       bruce julian     29 dec 1976  
!       modified by d gubbins 23 july 1977  
!       further modified by s roecker 1981  
!       This version allows inputs of lat or colat
! 
!     Note: 1) All angles are in radians
!           2) North latitude and east longitude are positive
!           3) Azimuth clockwise from north
!           4) If one point is at north or south pole, azimuth
!              from that point will be the limit of the azimuth as 
!              the pole is approached along the merdidian whose 
!              longitude is given 
!
!     input      meaning
!     -----      ------- 
!     alon0      longtitude of reference point (radians)
!     alon       longitude of point (radians)
!     lattype    =1 -> zlat0 = geocentric latitude of reference point
!                      zlat  = geocentric latitude of point
!                =2 -> zlat0 = geocentric colatitude of reference pt
!                      zlat  = geocentric colatitude of reference pt
!      zlat0     latitude or co-latitude of reference point (radians)
!      zlat      latitude or co-latitude of point (radians)
!
!      output    meaning
!      ------    ------- 
!      az1       azimuth of point from the reference point (radians)
!      az2       azimuth of reference point from point (radians)
!      delt      epicentral distance 
!  
!.... variable declarations
      implicit double precision (a-h,o-z)
      real*8, intent(in) :: zlat0,zlat,alon0,alon
      integer*4, intent(in) :: lattype
      real*8, intent(out) :: az1,az2,delt
!.... local variables
      data pii/3.1415926535897931d0/
      data twopi/6.2831853071795862d0/
!
!----------------------------------------------------------------------!
!
      if (lattype.eq.1) then
         alat0=pii*0.5d0-dble(zlat0)
         alat=pii*0.5d0-dble(zlat)
      else
         alat0 = dble(zlat0)
         alat = dble(zlat)
      end if
      st0 = dsin(alat0)
      ct0 = dcos(alat0)
      ct1 = dcos(alat)
      s0c1 = st0*ct1
      st1 = dsin(alat)
      s1c0 = st1*ct0
      dlon = dble(alon) - dble(alon0)
      sdlon = dsin(dlon)
      cdlon = dcos(dlon)
      cdelt = st0*st1*cdlon + ct0*ct1
      b = s0c1 - s1c0*cdlon
      a = st1*sdlon
      sdelt = dsqrt(b*b+a*a)
      ddelt = datan2(sdelt, cdelt)
      aze = 0.d0
      if (sdelt.ne.0.d0) aze = datan2(a,b)
!.... calculate back azimuth  
      a = -sdlon*st0
      b = s1c0 - s0c1*cdlon
      azs =pii
      if (sdelt.ne.0.d0) azs = datan2(a,b)
!.... make  0 < azimuth < twopi  
      if(aze .lt. 0.d0) aze = aze + twopi
      if(azs .lt. 0.d0) azs = azs + twopi
      delt = ddelt
      az1 = aze
      az2 = azs
      return
      end
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*8 FUNCTION GLAT(HLAT)
!
!     Convert geographic latitude to geocentric latitude
!
!     INPUT      MEANING
!     -----      ------- 
!     HLAT       geographic latitude in radians (north positive)
!
!     RESULT     MEANING
!     ------     ------- 
!     GLAT       geocentric latitude in radians (north positive)
!
!.... variable declarations
      IMPLICIT REAL*8 (A-H, O-Z)
      DATA HALFPI /1.5707963267948966D0/,
     ;     POLEFAC/0.010632D0/, ELFAC/0.993277D0/
!----------------------------------------------------------------------!
      IF (HALFPI - DABS(HLAT).GE.0.05D0) THEN
         GLAT = DATAN(ELFAC*DSIN(HLAT)/DCOS(HLAT))
      ELSE !special formula near pole
         GLAT = HLAT/ELFAC - DSIGN(POLEFAC, HLAT)
      ENDIF
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*8 FUNCTION GLATINV(HLAT)
!
!     Convert geocentric latitude back to geographic latitude
!
!     INPUT      MEANING
!     -----      ------- 
!     HLAT       geographic latitude in radians (north postive)
!
!     RESULT     MEANING
!     ------     -------
!     GLATINV    geographic latitude in radisn (north postivie)
!
!.... variable declarations
      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*8, INTENT(IN) :: HLAT
      DATA HALFPI /1.5707963267948966D0/,
     ;     POLEFAC/0.010632D0/, ELFAC/0.993277D0/
!----------------------------------------------------------------------!
      IF (HALFPI - DABS(HLAT).GE.0.05D0) THEN
         GLATINV = DATAN(DSIN(HLAT)/DCOS(HLAT)/ELFAC)
      ELSE !special formula near pole
         GLATINV = (HLAT + DSIGN(POLEFAC, HLAT))*ELFAC
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE VINCENTY(LVERB, XLAT0,XLON0,XLAT1,XLON1,  
     ;                    DSTDEG,DIST,AZ,BAZ,IERR) 
! 
!     Calculates the distance between two points using Vincenty's 
!     formula.  For more see: 
!       http://en.wikipedia.org/wiki/Vincenty's_formulae
!       http://www.movable-type.co.uk/scripts/latlong-vincenty.html
!       http://www.icsm.gov.au/gda/gdatm/gdav2.3.pdf, pp. 15
!
!     If speed is an issue we can use Newton's method and 
!     Squire-Trapp for high precision numerical differentiation
!
!     - B. Baker January 2013
!
!     INPUT      MEANING
!     -----      ------- 
!     LVERB      verbosity level
!     XLAT0      latitude of source (degrees) 
!     XLAT1      latitude of receiver (degrees)
!     XLON0      longitude of source (degrees)
!     XLON1      longitude of receiver (degrees)
!
!     OUTPUT     MEANING
!     ------     -------
!     AZ         azimuth (degrees)
!     BAZ        back azimuth (degrees) 
!     DIST       distance between points (km)
!     DSTDEG     distance between points (degrees)
!     IERR       error flag; failure to converge
!
!.... variable declarations
      REAL*8, INTENT(IN)  :: XLAT0,XLAT1, XLON0,XLON1 
      LOGICAL*4, INTENT(IN) :: LVERB
      REAL*8, INTENT(OUT) :: DSTDEG, DIST, AZ, BAZ 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8 F, A, B, RLAM1, RLAM2, PHI1, PHI2, U1, U2, 
     ;       COSU1, COSU2, SINU1, SINU2, SINL, COSL, SINS, COSS, 
     ;       COS2S, COS2S2, SINA, COS2A, C, DS, SIGMA, 
     ;       A1, B1, A2, DL, RLAM, RLAM0 
      REAL*8 A84, F84, PI180, TWOPI, TOL, EPS, PI, HPI 
      LOGICAL*4 LEQ 
      PARAMETER(A84 = 6378137.D0) !WGS-84
      PARAMETER(F84 = 0.0033528106647474805D0) !1/298.257223563; WGS-84
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180 
      PARAMETER(TWOPI = 6.2831853071795862D0) !2*pi
      PARAMETER(PI = 3.1415926535897931D0) !pi
      PARAMETER(HPI = 1.5707963267948966D0) !pi/2
      PARAMETER(TOL = 1.D-12) !~0.006 milli-meters accurate
      PARAMETER(EPS = 2.22D-15) 
      PARAMETER(MAXIT = 1000) 
!
!----------------------------------------------------------------------!
!
!.... check the points arent the same
      IERR = 0
      IF (DABS(XLON1 - XLON0).LT.EPS .AND.
     ;    DABS(XLAT1 - XLAT0).LT.EPS) THEN
         IF (LVERB) WRITE(*,*) 'vincenty: Warning points are same!'
         DIST = 0.D0
         AZ = 0.D0
         BAZ = 0.D0
         RETURN
      ENDIF
!.... initialize projection
      F = F84 !flattening
      A = A84 !radius at equator
!.... initialize constants
      B = (1.D0 - F)*A !radius at the poles
      RLAM1 = XLON0*PI180 !longitude of source
      RLAM2 = XLON1*PI180 !longitude of receiver 
      PHI1  = XLAT0*PI180 !latitude of source
      PHI2  = XLAT1*PI180 !latitude of receiver
!.... correct for errors at poles by adjusting to 0.6mm (vdist.m)
      IF (DABS(HPI - DABS(PHI1)).LT.1.D-10) 
     ;PHI1 = DSIGN(PHI1,1.D0)*(HPI - 1.D-10) 
      IF (DABS(HPI - DABS(PHI2)).LT.1.D-10) 
     ;PHI2 = DSIGN(PHI2,1.D0)*(HPI - 1.D-10) 
!.... correction at equator 
      LEQ = .FALSE. !points aren't on equator
      IF (DABS(PHI1).LT.TOL .AND. DABS(PHI2).LT.TOL) LEQ = .TRUE.
      U1 = DATAN2((1.D0 - F)*DTAN(PHI1),1.D0) !reduced latitude
      U2 = DATAN2((1.D0 - F)*DTAN(PHI2),1.D0) !reduced latitude
      COSU1 = DCOS(U1)
      COSU2 = DCOS(U2)
      SINU1 = DSIN(U1)
      SINU2 = DSIN(U2)
      DL = RLAM2 - RLAM1 !difference in latitude is first guess 
      !IF (DL > PI) DL = TWOPI - DL
!
!.... iterative loop
      RLAM = DL !initial guess
      RLAM0 = RLAM
      DO 100 K=1,MAXIT
         SINL = DSIN(RLAM)
         COSL = DCOS(RLAM)
         SINS = DSQRT( (COSU2*SINL)**2
     ;               + (COSU1*SINU2 - SINU1*COSU2*COSL)**2)
         COSS = SINU1*SINU2 + COSU1*COSU2*COSL
         SIGMA = DATAN2(SINS,COSS)
         SINA = COSU1*COSU2*SINL/SINS
         COS2A = 1.D0 - SINA**2
         IF (.NOT.LEQ) THEN
            COS2S = COSS - 2.D0*SINU1*SINU2/COS2A
         ELSE
            COS2S =-1.D0
         ENDIF
         COS2S2 = COS2S**2
         C = F/16.D0*COS2A*(4.D0 + F*(4.D0 - 3.D0*COS2A))
         RLAM = DL + (1.D0 - C)*F*SINA
     ;        *(SIGMA + C*SINS*(COS2S + C*COSS*(-1.D0 + 2.D0*COS2S2)))
         IF (DABS(RLAM - RLAM0).LT.TOL) GOTO 105
         IF (DABS(SINS).LT.EPS) THEN
            WRITE(*,*) 'vincenty: Points are diametrically opposed'
            RLAM = PI  
            GOTO 105 
         ENDIF
         RLAM0 = RLAM
  100 CONTINUE
      WRITE(*,*) 'vincenty: Failed to converge!'
      IERR = 1
      DIST = 0.D0
      AZ = 0.D0
      BAZ = 0.D0
      RETURN
  105 CONTINUE
!
!.... evaluate constants corresponding to lambda 
      U2 = COS2A*( (A**2 - B**2)/B**2 )
      A1 = 1.D0
     ;   + U2/16384.D0*(4096.D0 + U2*(-768.D0 + U2*(320.D0 -175.D0*U2)))
      B1 = U2/1024.D0*(256.D0 + U2*(-128.D0 + U2*(74.D0 - 47.D0*U2)))
      DS = B1*SINS*(COS2S + 0.25D0*B1*(COSS*(-1.D0 + 2.D0*COS2S**2)
     ;                    - 1.D0/6.D0*B1*COS2S*(-3.D0 + 4.D0*SINS**2)
     ;                     *(-3.D0 + 4.D0*COS2S**2)))
!.... distance, azimuth, back-azimuth
      DSTDEG = (SIGMA - DS)/PI180 !radians -> degrees
      DIST = B*A1*(SIGMA - DS)*1.D-3 !m -> km
      A1 = DATAN2(COSU2*SINL, COSU1*SINU2 - SINU1*COSU2*COSL)
      A2 = DATAN2(COSU1*SINL,-SINU1*COSU2 + COSU1*SINU2*COSL)
      AZ = A1/PI180 !azimuth point 1 -> point 2
!.... a2/pi180 is angle point 2 -> point 1, want it facing other way
!.... so subtract 180, but this is a back azimuth so add 360
      BAZ = 180.D0 + A2/PI180 !360 - 180 + angle point 2 - > point 1 

      IF (AZ .LT.0.D0  ) AZ  = AZ  + 360.D0
      IF (AZ .GT.360.D0) AZ  = AZ  - 360.D0
      IF (BAZ.LT.0.D0  ) BAZ = BAZ + 360.D0
      IF (BAZ.GT.360.D0) BAZ = BAZ - 360.D0

      IF (LVERB) THEN
         WRITE(*,905) DIST, AZ, BAZ
  905    FORMAT(' vincenty: Distance in km between points:',F12.3,/, 
     ;          '           Azimuth in degrees:',F8.3,/, 
     ;          '           Back-azimuth degrees:',F8.3,/)
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE VINDIRECT(XLAT0,XLON0,S,AZ, XLAT1,XLON1,A21,IERR)
!
!     Vincenty's direct algorithm:
!     Given a starting latitude and longitude, distance, and azimuth
!     computes the new latitude and longitude along that azimtuh 
!     for the given length.  for more see:
!     http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html
!     http://en.wikipedia.org/wiki/Vincenty's_formulae#Direct_Problem
!     vreckon.m
!
!     INPUT      MEANING
!     -----      ------- 
!     AZ         azimuth to travel along (degrees)
!     S          distance to travel (m) 
!     XLAT0      initial latitude 
!     XLON0      initial longitude 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     A21        back-azimuth (degrees)
!     IERR       error flag, failed to converge
!     XLAT1      destination latitude
!     XLON1      destination longitude
!   
!.... variable declarations 
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: XLAT0,XLON0,S,AZ
      REAL*8, INTENT(OUT) :: XLAT1,XLON1,A21
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8 F,A,B,C,A1,B1, PHI0,RLAM0,PHI1,RLAM1,AZ1,SINA1,COSA1, 
     ;       SIGMA,SIGMAP,SIGMA1,U2,DS, 
     ;       SINA,SIN2A,COS2A,TANU1,COSU1,SINU1,
     ;       SINS,COSS,COS2S,COS2S2, XNUM,XDEN,ARG,DL   
      INTEGER*4 K
      REAL*8 A84, B84, F84, PI180, TWOPI, TOL, PI, EPS 
      INTEGER*4 MAXIT
      PARAMETER(A84 = 6378137.D0) !WGS-84
      PARAMETER(B84 = 6356752.3142D0) !WGS-84 radius at pole
      PARAMETER(F84 = 0.0033528106647474805D0) !1/298.257223563; WGS-84
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180 
      PARAMETER(PI = 3.1415926535897931D0) !pi 
      PARAMETER(TWOPI = 6.2831853071795862D0) !2*pi
      PARAMETER(TOL = 1.D-12) !~0.006 milli-meters accurate
      PARAMETER(EPS = 2.22D-15)
      PARAMETER(MAXIT = 1000)

!
!----------------------------------------------------------------------!
!
!.... early return at same point
      IERR = 0
      IF (DABS(S).LT.1.11D-7) THEN
         XLAT1 = XLAT0
         XLON1 = XLON0
         A21 = 0.D0
         RETURN
      ENDIF
!.... initialize projection
      IERR = 0 
      F = F84 !flattening
      A = A84 !radius at equator
!.... initialize constants
      B = (1.D0 - F)*A !radius at the poles
      PHI0  = XLAT0*PI180 !latitude of source
      RLAM0 = XLON0*PI180 !longitude of source
      AZ1   = AZ*PI180    !azimuth of source to receiver
      SINA1 = DSIN(AZ1)
      COSA1 = DCOS(AZ1)
!.... initialize
      TANU1 = (1.D0 - F)*DTAN(PHI0)
      COSU1 = 1.D0/(DSQRT(1.D0 + TANU1**2))
      SINU1 = 1.D0*TANU1*COSU1
      SIGMA1 = DATAN2(TANU1,COSA1)

      SINA = COSU1*SINA1
      SIN2A = SINA**2
      COS2A = 1.D0 - SIN2A
      U2 = COS2A*(A**2 - B**2)/(B**2)
      A1 = 1.D0
     ;   + U2/16384.D0*(4096.D0 + U2*(-768.D0 + U2*(320.D0 -175.D0*U2)))
      B1 = U2/1024.D0*(256.D0 + U2*(-128.D0 + U2*(74.D0 - 47.D0*U2)))
      SIGMA = S/(B*A1)
      SIGMAP = TWOPI
      DO 1 K=1,MAXIT
         COS2S = DCOS(2.D0*SIGMA1 + SIGMA)
         SINS = DSIN(SIGMA)
         COSS = DCOS(SIGMA)
         COS2S2 = COS2S**2
         DS = B1*SINS*(COS2S + 0.25D0*B1*(COSS*(-1.D0 + 2.D0*COS2S2)
     ;                       - B1/6.D0*COS2S*(-3.D0 + 4.D0*SINS**2)
     ;                        *(-3.D0 + 4.D0*COS2S2)))
         SIGMAP = SIGMA
         SIGMA = S/(B*A1) + DS
         IF (DABS(SIGMAP - SIGMA).LT.TOL) GOTO 100 !convergence?
    1 CONTINUE
      WRITE(*,*) 'vindirect: Failed to converge!'
      IERR = 1
      XLAT1 = 0.D0
      XLON1 = 0.D0
      A21 = 0.D0
      RETURN
  100 CONTINUE
!
!.... finish
      XNUM = SINU1*COSS + COSU1*SINS*COSA1
      ARG = SIN2A + (SINU1*SINS - COSU1*COSS*COSA1)**2
      XDEN = (1.D0 - F)*DSQRT(ARG)
      PHI1  = DATAN2(XNUM,XDEN)
      XNUM = SINS*SINA1
      XDEN = COSU1*COSS - SINU1*SINS*COSA1
      RLAM1 = DATAN2(XNUM,XDEN)
      C = F/16.D0*COS2A*(4.D0 + F*(4.D0 - 3.D0*COS2A))
      DL = RLAM1 - (1.D0 - C)*F*SINA
     ;    *(SIGMA + C*SINS*(COS2S + C*COSS*(-1.D0 + 2.D0*COS2S2)))
      XLAT1 = PHI1/PI180
      XLON1 = (RLAM0 + DL)/PI180
      !XLON1 = DMOD(XLON1,360.D0) ![0,360] convention 
      IF (XLON1.LT.  0.D0) XLON1 = XLON1 + 360.D0
      IF (XLON1.GT.360.D0) XLON1 = XLON1 - 360.D0

      A21 = DATAN2(SINA,-SINU1*SINS + COSU1*COSS*COSA1)
      A21 = 180.D0 + A21/PI180 !direction reversal
      IF (A21.LT.  0.D0) A21 = A21 + 360.D0
      IF (A21.GT.360.D0) A21 = A21 - 360.D0
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE XPTS2DIST3(NXPTS, CSIDE,SLAT,SLON, XMOD0,XMOD1, 
     ;                      XLATMIN,XLONMIN, XLATMAX,XLONMAX, XLOCS, 
     ;                      XDISTS,AZIMS,IERR) 
!
!     This routine calculates the latitude and longitude of points 
!     along the azimuth from xmod0 to xmod1 or xmod1 to xmod0 
!     depending on which way the source is approaching.  Additionally
!     we calculate the azimuths to each point. 
!
!     INPUT      MEANING
!     -----      -------
!     SLAT       source latitude (degrees)
!     SLON       source longitude (degrees)
!     XLATMAX    latitude at xmod1 (degrees)
!     XLATMIN    latitude at xmod0 (degrees)
!     XLOCS      x locations to evaluate distance and azimtuth (m)
!     XLONMAX    longitude of xmod1 (degrees)
!     XLONMIN    longitude at xmod0 (degrees)
!     XMOD0      left Bielak/absorbing boundary x position (m)
!     XMOD1      right Bielak/absorbing boundary x position (m)
!
!     OUTPUT     MEANING
!     ------     ------- 
!     AZIMS      azimuths from source to model points (deg)
!     IERR       error flag
!     XDISTS     distances from source to model points (km)
! 
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: CSIDE 
      REAL*8, INTENT(IN) :: XLOCS(NXPTS), SLAT,SLON, XMOD0,XMOD1, 
     ;                      XLATMIN,XLONMIN, XLATMAX,XLONMAX 
      INTEGER*4, INTENT(IN) :: NXPTS
      REAL*8, INTENT(OUT) :: XDISTS(NXPTS), AZIMS(NXPTS)
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      REAL*8 DIST,AZ,BAZ,AZPTH, XLAT0,XLON0, XLAT1,XLON1, S, 
     ;       DSTDEG
      INTEGER*4 IX 
!
!----------------------------------------------------------------------!
!
      IERR = 0 
      IF (CSIDE.EQ.'L') THEN
         CALL VINCENTY(.FALSE., XLATMIN,XLONMIN,XLATMAX,XLONMAX, 
     ;                 DSTDEG,DIST,AZPTH,BAZ,IERR)  
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'x2daz: Error calling vincent1'
            RETURN
         ENDIF
         XLAT0 = XLATMIN
         XLON0 = XLONMIN
      ELSE
         CALL VINCENTY(.FALSE., XLATMAX,XLONMAX,XLATMIN,XLONMIN, 
     ;                 DSTDEG,DIST,AZPTH,BAZ,IERR) 
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'x2daz: Error calling vincent2'
            RETURN
         ENDIF
         XLAT0 = XLATMAX
         XLON0 = XLONMAX
      ENDIF
!
!.... calculate distance using bjdaz2 convention
      DO 1 IX=1,NXPTS
         IF (CSIDE.EQ.'L') THEN
            S = DABS(XLOCS(IX) - XMOD0)
         ELSE
            S = DABS(XMOD1 - XLOCS(IX))
         ENDIF
         CALL VINDIRECT(XLAT0,XLON0,S,AZPTH, XLAT1,XLON1,BAZ,IERR)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'x2daz: Error calling vindirect3'
            XDISTS(1:NXPTS) = 0.D0
            AZIMS(1:NXPTS) = 0.D0
            RETURN
         ENDIF
!
!....... now calculate the distance and azimuth from the source
         CALL VINCENTY(.FALSE., SLAT,SLON, XLAT1,XLON1, 
     ;                 DSTDEG,DIST,AZ,BAZ,IERR)  
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'x2daz: Error calling vincenty4'
            XDISTS(1:NXPTS) = 0.D0
            AZIMS(1:NXPTS) = 0.D0
         ENDIF
         XDISTS(IX) = DIST
         AZIMS(IX)  = AZ 
         !print *, xdists(ix), azims(ix) 
c        RLON = XLON1*PI180
c        RLAT = HPI - GLAT(XLAT1*PI180)
c        CALL BJDAZ2(SLAT0,SLON0, RLAT,RLON, DELT1,AZ1,AZ2,0)
!        XDISTS(IX) = DELT1/PI180*KMPERDEG
!        AZIMS(IX) = AZ
    1 CONTINUE
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE XPTS2DIST2(NXPTS, CSIDE,SLAT,SLON, XMOD0,XMOD1, 
     ;                      XLATMIN,XLONMIN, XLATMAX,XLONMAX, XLOCS,
     ;                      XDIST,IERR) 
!
!     This routine calculates the latitude and longitude away from the 
!     appropriate Bielak/absorbing side, then uses BJDAZ2 to calculate 
!     the great circle distance.  We rely on BJDAZ2 since it is used 
!     in Steve's code and I want things to be consistent.  Fortunately
!     we aren't working at the poles so we can only expect to be off 
!     a few kilometers.  - B. Baker January 2013
!
!     INPUT      MEANING
!     -----      -------
!     CSIDE      using left or right side of model
!     NXPTS      number of x points to convert
!     XLATMAX    latitude (deg) of right bielak/abs bdry
!     XLATMIN    latitude (deg) of left bielak/abs bdry
!     XLONMAX    longitude (deg) of right bielak/abs bdry
!     XLONMIN    longitude (deg) of left bielak/abs bdry
!     SLAT       source latitude
!     SLON       source longitude
!     XLOCS      x locations (m) in bielak layer to calc gc distances 
!     XMOD0      x bielak/absorbing boundary left
!     XMOD1      x bielak/absorbing boundary right
!
!     OUTPUT     MEANING
!     ------     -------
!     IERR       an errror occurred in Vincenty direct calculation
!     XDIST      gc distance (degrees) from model side to bielak node
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: CSIDE 
      REAL*8, INTENT(IN) :: XLOCS(NXPTS), SLAT,SLON, XLATMIN,XLONMIN, 
     ;                      XLATMAX,XLONMAX, XMOD0, XMOD1 
      INTEGER*4, INTENT(IN) :: NXPTS
      REAL*8, INTENT(OUT) :: XDIST(NXPTS)  
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8 SLAT0, SLON0, AZ, BAZ, XLAT0, XLON0, XLAT1, XLON1,  
     ;       RLAT, RLON, S, AZ1, AZ2, DELT1, DELT, DIST, DSTDEG
      REAL*8 PI180, HPI, KMPERDEG
      INTEGER*4 IX
      PARAMETER(PI180 = 0.017453292519943295D0)
      PARAMETER(KMPERDEG = 111.069365447154D0)
      PARAMETER(HPI = 1.5707963267948966D0)
      REAL*8 GLAT 
!
!----------------------------------------------------------------------!
!
!.... get the azimuth between the points 
      IERR = 0
      SLAT0 = HPI - GLAT(SLAT*PI180)
      SLON0 = SLON*PI180
      IF (CSIDE.EQ.'L') THEN
         CALL VINCENTY(.FALSE., XLATMIN,XLONMIN,XLATMAX,XLONMAX, 
     ;                 DSTDEG,DIST,AZ,BAZ,IERR)  
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'xpts2dist2: Error calling vincent1'
            RETURN
         ENDIF
         XLAT0 = XLATMIN
         XLON0 = XLONMIN
      ELSE
         CALL VINCENTY(.FALSE., XLATMAX,XLONMAX,XLATMIN,XLONMIN, 
     ;                 DSTDEG,DIST,AZ,BAZ,IERR) 
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'xpts2dist2: Error calling vincent2'
            RETURN
         ENDIF
         XLAT0 = XLATMAX
         XLON0 = XLONMAX
      ENDIF
      RLAT = HPI - GLAT(XLAT0*PI180)
      RLON = XLON0*PI180
      CALL BJDAZ2(SLAT0,SLON0, RLAT,RLON, DELT,AZ1,AZ2,0) 
!
!.... calculate distance using bjdaz2 convention
      DO 1 IX=1,NXPTS
         IF (CSIDE.EQ.'L') THEN
            S = DABS(XLOCS(IX) - XMOD0)
         ELSE
            S = DABS(XMOD1 - XLOCS(IX))
         ENDIF
         CALL VINDIRECT(XLAT0,XLON0,S,AZ, XLAT1,XLON1,BAZ,IERR)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'xpts2dist2: Error calling vindirect2'
            XDIST(1:NXPTS) = 0.D0
            RETURN
         ENDIF
         RLON = XLON1*PI180
         RLAT = HPI - GLAT(XLAT1*PI180)
         CALL BJDAZ2(SLAT0,SLON0, RLAT,RLON, DELT1,AZ1,AZ2,0) 
         XDIST(IX) = DELT1/PI180*KMPERDEG 
    1 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE XPTS2DIST(NXPTS,LZOFF,POFF, SLAT,SLON, XLATMIN,XLONMIN,
     ;                     XLATMAX,XLONMAX,XMOD0,XMOD1, XLOCS, XDIST)
!
!     Converts x locations to their distance from the source.  
!     This works, because as it happens, surface waves travel 
!     along the surface and can be converted to kilometers.  
!  
!     Points are linearly interpolated with dtform
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: XLOCS(NXPTS), SLAT,SLON, XLATMIN,XLONMIN,
     ;                      XLATMAX,XLONMAX, XMOD0,XMOD1, POFF
      INTEGER*4, INTENT(IN) :: NXPTS
      LOGICAL*4, INTENT(IN) :: LZOFF
      REAL*8, INTENT(OUT) :: XDIST(NXPTS)
!.... local variables
      REAL*8 SLAT0,SLON0,GRLAT,GRLON, DELT1,DELT2,AZ1,AZ2,
     ;       XMOD0KM,XMOD1KM,XKM
      INTEGER*4 IX
      REAL*8 PI180, HPI, KMPERDEG
      PARAMETER(PI180 = 0.017453292519943295D0)
      PARAMETER(KMPERDEG = 111.069365447154D0)
      PARAMETER(HPI = 1.5707963267948966D0)
      REAL*8 GLAT, DTFORM
!
!----------------------------------------------------------------------!
!
      SLAT0 = HPI - GLAT(SLAT*PI180)
      SLON0 = SLON*PI180
      GRLAT = HPI - GLAT(XLATMIN*PI180)
      GRLON = XLONMIN*PI180
      CALL BJDAZ2(SLAT0,SLON0, GRLAT,GRLON, DELT1,AZ1,AZ2,0)
      GRLAT = HPI - GLAT(XLATMAX*PI180)
      GRLON = XLONMAX*PI180
      CALL BJDAZ2(SLAT0,SLON0, GRLAT,GRLON, DELT2,AZ1,AZ2,0)
!
!.... convert radians -> degrees -> kilometers
      DELT1 = DELT1/PI180*KMPERDEG
      DELT2 = DELT2/PI180*KMPERDEG
!.... convert grid from meters -> kilomters
      IF (LZOFF) THEN
         XMOD0KM = (XMOD0 - POFF)*1.D-3 
         XMOD1KM = (XMOD1 - POFF)*1.D-3
      ELSE
         XMOD0KM = XMOD0*1.D-3
         XMOD1KM = XMOD1*1.D-3
      ENDIF
!
!.... stick onto grid
      DO 1 IX=1,NXPTS
         XKM = XLOCS(IX)*1.D-3
         XDIST(IX) = DTFORM(XMOD0KM,XMOD1KM, DELT1,DELT2, XKM)
    1 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*8 FUNCTION DTFORM(A,B, C,D,X)
! 
!     Converts numbers on the interval x in [A,B] -> xi in [c,d]
!
!     From Hughes we have [A,B] -> [-1,1] is given by solving 
!     c_1, c_2 in 
!       -1 = c_1 + A*c_2
!        1 = c_1 + B*c_2 
!     so more generally we can do 
!        c = c_1 + A*c_2 
!        b = c_2 + B*c_2  
! 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8, INTENT(IN) :: A,B, C,D,X
      IF (B.EQ.A) THEN
         WRITE(*,*) 'dtform: Error A = B, determinant undefined'
         DTFORM = 1.D0
         RETURN
      ENDIF
      DET = 1.D0/(B - A)
      C1 = DET*(B*C - A*D)
      C2 = DET*(-C + D)
      XI = C1 + X*C2
      DTFORM = XI
      RETURN
      END


