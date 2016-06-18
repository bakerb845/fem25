!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDREC_HD(PROJNM, NREC,IERR)
!
!     Reads number of frequencies in .freqs file
!
!     INPUT      MEANING
!     -----      ------- 
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     NREC       number of receivers 
!
!.... variable declarations 
      CHARACTER(80), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(OUT) :: NREC,IERR
      CHARACTER(80) FLNAME
      LOGICAL*4 LEX 
      PARAMETER(IUNIT = 44) 
!
!----------------------------------------------------------------------!
!
      IERR = 0 
      NREC = 0
      FLNAME(1:80) = ' ' 
      FLNAME = TRIM(PROJNM)//'.recv'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rdrec_hd: Error cannot locate recv file'
         IERR =-1 
         RETURN
      ENDIF
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME))
      READ(IUNIT,*,END=60) !header
      READ(IUNIT,*,END=60) NREC 
      CLOSE(IUNIT)
      RETURN
   60 CONTINUE
      IERR = 2
      WRITE(*,*) 'rdrec_hd: Error premature end of recv file'
      CLOSE(IUNIT) 
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDREC(PROJNM, NREC, LFSURF,XREC,YREC,ZREC, IERR)
!
!     Reads number of frequencies in .freqs file
!
!     INPUT      MEANING
!     -----      ------- 
!     NREC       number of receivers
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     LFSURF     True -> force receiver to be at free surface
!     XREC       x reciever location
!     YREC       y receiver location
!     ZREC       z receiver location
!
!.... variable declarations 
      CHARACTER(80), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(IN) :: NREC
      REAL*8, INTENT(OUT) :: XREC(NREC), YREC(NREC), ZREC(NREC)
      LOGICAL*4, INTENT(OUT) :: LFSURF(NREC)
      INTEGER*4, INTENT(OUT) :: IERR
      CHARACTER(80) FLNAME
      PARAMETER(IUNIT = 44)
!
!----------------------------------------------------------------------!
!
!.... initialize
      IERR = 0 
      FLNAME(1:80) = ' '
      FLNAME = TRIM(PROJNM)//'.recv'
      FLNAME = ADJUSTL(FLNAME)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME))
      READ(IUNIT,*,END=60) !header
      READ(IUNIT,*,END=60) !NREC
      READ(IUNIT,*,END=60) !header 2
      DO 1 IREC=1,NREC
         READ(IUNIT,*,END=60) LFSURF(IREC), 
     ;                        XREC(IREC),YREC(IREC),ZREC(IREC)
    1 CONTINUE
      CLOSE(IUNIT)
      RETURN
   60 CONTINUE
      IERR = 1 
      WRITE(*,*) 'rdrec: Premature end of file'
      CLOSE(IUNIT)
      RETURN 
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDSRC_EQ(PROJNM, XMLAT0,XMLON0,XMLAT1,XMLON1,
     ;                    AZMOD, SRC,IERR)
!
!     Reads source list for body waves and corrects back-azimuth.  
!     Now handles memory allocation within routine. 
!
!     If the model is striking east-west, (i.e., azmod = 0),
!     then  X is positive North, and Y is positive East.  Hence if the 
!     back azimuth (baz) is 180 degrees the incoming wave will be 
!     advancing in the +X (pure north) direction. In general, the 
!     slowness vector (px, py, pz) will be:
!        (sin(ai)cos(baz-pi), sin(ai)sin(baz-pi), cos(ai))/v
!     The azimuth of the model is defined by the direction of the +X 
!     axis relative to north.  Hence for a non-zero azmod we modify 
!     above form for (px,py,pz) to be
!        (sin(ai)cos(baz-pi-azmod), sin(ai)sin(baz-pi-azmod), cos(ai))/v
!     However, bazn has been corrected when read in 

!
!     INPUT      MEANING
!     -----      ------- 
!     AZMOD      model azimuth (degrees)
!     PROJNM     project name
!     XMLAT0     latitude of left model point (degrees)
!     XMLAT1     latitude of right model point degrees)
!     XMLON0     longitude of left point (degrees)
!     XMLON1     longitude of right model (degrees)
!
!     OUTPUT     MEANING
!     ------     ------- 
!     BAZN       corrected back-azimuth of source; see above
!     CSIDE      'L' source coming from left, 'R' source from right
!     DIP        dip of fault +down from horizontal (degrees) [0,90]
!     IERR       error flag
!     MODE       for surface waves, the mode to model
!     NSRC       number of sources
!     SLAT       source latitude (degrees)
!     SLON       source longitude (degrees)
!     RAKE       angle between slip and strike (degrees) [-180,180]
!     SMAG       magnitude in dyne-cm (1E-7 newton-meters) 
!     STRIKE     strike of fault +clockwise from north (degrees) [0,360]
!
!.... variable declarations 
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: PROJNM
      REAL*8, INTENT(IN) :: AZMOD, XMLON0, XMLAT0, XMLAT1, XMLON1 
      INCLUDE 'fwd_struc.h'
      TYPE (SRC_INFO) SRC
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      CHARACTER(80) FLNAME
      REAL*8 BAZ, TLAT, TLON, AZ1, DELT, PI180, HPI, 
     ;       DSTDEG, KMPERDEG, AZ0, BAZ0, DELT0, DELT1, BAZ1
      REAL*8 GLAT
      REAL*8 MW2DCM 
      INTEGER*4 ISRC, IUNIT 
      PARAMETER(IUNIT = 44) 
      PARAMETER(HPI = 1.5707963267948966D0)     !pi/2
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
      PARAMETER(KMPERDEG = 111.069365447154D0)
!       
!----------------------------------------------------------------------!
!   
      IERR = 0
      FLNAME(1:80) = ' '
      FLNAME = TRIM(PROJNM)//'.src'
      FLNAME = ADJUSTL(FLNAME)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME))
      READ(IUNIT,*,END=60) !header
      READ(IUNIT,*,END=60) src%NSRC
      READ(IUNIT,*,END=60) !header 2
!
!.... set space
      ALLOCATE(src%SRCTYP(SRC%NSRC))  !character descriptor of source
      ALLOCATE(src%CSIDE(SRC%NSRC))   !'L' -> source from left, 
                                      !'R' -> source from right
      ALLOCATE(src%AOI(SRC%NSRC))     !angle of incidence 
      ALLOCATE(src%BAZN(SRC%NSRC))    !azimuth (corrected)
      ALLOCATE(src%SLAT(SRC%NSRC))    !source latitude(deg)
      ALLOCATE(src%SLON(SRC%NSRC))    !source longitude (deg)
      ALLOCATE(src%SDEP(SRC%NSRC))    !source depth (km)
      ALLOCATE(src%STRIKE(SRC%NSRC))  !strike (deg)
      ALLOCATE(src%DIP(SRC%NSRC))     !dip (deg)
      ALLOCATE(src%RAKE(SRC%NSRC))    !rake (deg)
      ALLOCATE(src%SMAG(SRC%NSRC))    !magnitude (dyne-cm)
      ALLOCATE(src%MODE(SRC%NSRC))    !mode to model 

      src%NSRC_SRF = 0 
      src%NSRC_BDY = 0
      TLON = XMLON0*PI180
      TLAT = HPI - GLAT(XMLAT0*PI180) 
      WRITE(*,*)
      DO 1 ISRC=1,src%NSRC
         READ(IUNIT,*,END=60) src%SRCTYP(ISRC), 
     ;                        src%SLAT(ISRC),   src%SLON(ISRC),
     ;                        src%SDEP(ISRC),   src%AOI(ISRC),
     ;                        src%STRIKE(ISRC), src%DIP(ISRC), 
     ;                        src%RAKE(ISRC),   src%SMAG(ISRC), 
     ;                        src%MODE(ISRC)
!....... capitalize
         IF (src%SRCTYP(ISRC)(1:1).EQ.'p') src%SRCTYP(ISRC)(1:1) = 'P'
         IF (src%SRCTYP(ISRC)(1:1).EQ.'s') src%SRCTYP(ISRC)(1:1) = 'S'
         IF (src%SRCTYP(ISRC)(2:2).EQ.'p') src%SRCTYP(ISRC)(2:2) = 'P'
         IF (src%SRCTYP(ISRC)(2:2).EQ.'s') src%SRCTYP(ISRC)(2:2) = 'S'
         IF (src%SRCTYP(ISRC)(2:2).EQ.'r') src%SRCTYP(ISRC)(2:2) = 'R'
         IF (src%SRCTYP(ISRC)(2:2).EQ.'l') src%SRCTYP(ISRC)(2:2) = 'L'
         IF (src%SRCTYP(ISRC)(2:2).EQ.'v') src%SRCTYP(ISRC)(2:2) = 'V'
         IF (src%SRCTYP(ISRC)(2:2).EQ.'b') src%SRCTYP(ISRC)(2:2) = 'B'
!
!....... calculate the azimuth for the wave traveling to the closest 
!....... point in the model 
         CALL VINCENTY(.FALSE.,src%SLAT(ISRC),src%SLON(ISRC), 
     ;                 XMLAT0,XMLON0, DSTDEG,DELT0,AZ0,BAZ0,IERR)
         IF (IERR /= 0) RETURN
         CALL VINCENTY(.FALSE.,src%SLAT(ISRC),src%SLON(ISRC), 
     ;                 XMLAT1,XMLON1, DSTDEG,DELT1,AZ1,BAZ1,IERR)
         IF (IERR /= 0) RETURN
         IF (DELT0.LT.DELT1) THEN !left model point closer than right 
            src%CSIDE(ISRC) = 'L'
            CALL VINCENTY(.FALSE.,src%SLAT(ISRC),src%SLON(ISRC), 
     ;                    XMLAT0,XMLON0, DSTDEG,DELT,AZ1,BAZ,IERR)
            IF (IERR /= 0) RETURN
         ELSE !right model point closer than left
            src%CSIDE(ISRC) = 'R' 
            CALL VINCENTY(.FALSE., src%SLAT(ISRC),src%SLON(ISRC), 
     ;                    XMLAT1,XMLON1, DSTDEG,DELT,AZ1,BAZ,IERR)
            IF (IERR /= 0) RETURN
         ENDIF
         DELT = DSTDEG 
!         IF (IERR.NE.0) THEN !copied from steve 
!            RLON = SRC%SLON(ISRC)*PI180 !the source is the "receiver"
!            RLAT = HPI - GLAT(SRC%SLAT(ISRC)*PI180) 
!!           CALL BJDAZ2(TLAT,TLON,RLAT,RLON,DELT,AZ1,AZ2,0) 
!            CALL BJDAZ2(RLAT,RLON,TLAT,RLON,DELT,AZ1,AZ2,0) 
!            AZ1 = AZ1/PI180
!            BAZ = AZ2/PI180
!            DELT = DELT/PI180
!         ENDIF
          src%BAZN(ISRC) = BAZ - 180.D0 - AZMOD !corrected azimuth
          IF (src%BAZN(ISRC).LT.0.D0) 
     ;    src%BAZN(ISRC) = src%BAZN(ISRC) + 360.D0
!
!........ round body waves to nearest integeer
          IF (src%SRCTYP(ISRC)(1:1) == 'P')  
     ;    src%BAZN(ISRC) = DFLOAT( INT(src%BAZN(ISRC) + 0.5) ) 
!
!....... LEFT -> RIGHT is direction of source
!        IF (src%CSIDE(ISRC) == 'L') THEN !model points LEFT to RIGHT
!           src%BAZN(ISRC) = AZMOD - AZ1 !corrected azimuth
!           IF (src%BAZN(ISRC).LT.0.D0) 
!    ;      src%BAZN(ISRC) = src%BAZN(ISRC) + 360.D0
!           print *, src%bazn(isrc)
!        ELSE !model points RIGHT to LEFT; flip AZMOD
!           src%BAZN(ISRC) = AZ1 - (AZMOD + 180.D0) 
!           print *, src%bazn(isrc),azmod+180.d0
!           IF (src%BAZN(ISRC) < 0.D0) THEN !go to upper half plane
!              src%BAZN(ISRC) = 180.D0 + DABS(src%BAZN(ISRC))
!           ELSE !to lower half plane
!              src%BAZN(ISRC) = 180.D0 - DABS(src%BAZN(ISRC))
!           ENDIF
!           print *, src%bazn(isrc), azmod+180.d0
!        ENDIF
c        IF (src%BAZN(ISRC).LT.0.D0) 
c    ;   src%BAZN(ISRC) = src%BAZN(ISRC) + 360.D0
c        IF (src%BAZN(ISRC).GT.360.D0) 
c    ;   src%BAZN(ISRC) = src%BAZN(ISRC) - 360.D0 
!        IF (DCOS(SRC%BAZN(ISRC)*PI180).LT.0.D0) THEN 
!           SRC%CSIDE(ISRC) = 'R'
         IF (src%CSIDE(ISRC) == 'L') THEN
            WRITE(*,805) DELT,BAZ,src%BAZN(ISRC) 
         ELSE
            WRITE(*,806) DELT,BAZ,src%BAZN(ISRC)
         ENDIF
  805    FORMAT(' rdsrc_eq: Wave approaching from model left',/, 
     ;          '           Distance to model:'     ,F12.4,' degrees',/,
     ;          '           Event to model azimuth:',F12.4,' degrees',/,
     ;          '           Corrected azimuth:'     ,F12.4,' degrees',/)
  806    FORMAT(' rdsrc_eq: Wave approaching from model right',/, 
     ;          '           Distance to model:'     ,F12.4,' degrees',/,
     ;          '           Event to model azimuth:',F12.4,' degrees',/,
     ;          '           Corrected azimuth:'     ,F12.4,' degrees',/)
         IF (src%STRIKE(ISRC).LT.0.D0) THEN
            WRITE(*,*) 'rdsrc_eq: Warning setting strike to 0'
            src%STRIKE(ISRC) = 0.D0
         ENDIF
         IF (src%STRIKE(ISRC).GT.360.D0) THEN
            WRITE(*,*) 'rdsrc_eq: Warning setting strike to 360'
            SRC%STRIKE(ISRC) = 360.D0
         ENDIF 
         IF (src%DIP(ISRC).LT.0.D0) THEN
            WRITE(*,*) 'rdsrc_eq: Warning setting dip to 0'
            src%DIP(ISRC) = 0.D0
         ENDIF
         IF (src%DIP(ISRC).GT.90.D0) THEN
            WRITE(*,*) 'rdsrc_eq: Warning setting dip to 90'
            src%DIP(ISRC) = 90.D0
         ENDIF
         IF (src%RAKE(ISRC).LT.-180.D0) THEN
            WRITE(*,*) 'rdsrc_eq: Warning setting rake to -180'
            SRC%RAKE(ISRC) =-180.D0
         ENDIF
         IF (src%RAKE(ISRC).GT. 180.D0) THEN
            WRITE(*,*) 'rdsrc_eq: Warning setting rake to +180'
            src%RAKE(ISRC) = 180.D0
         ENDIF
         IF (src%SMAG(ISRC).LE.0.D0) THEN
            WRITE(*,*) 'rdsrc_eq: Warning magnitude <= 0'
            WRITE(*,*) 'rdsrc_eq: Setting to Mw 5.5'
         ENDIF
         IF (src%SMAG(ISRC).LT.10.D0) THEN
            WRITE(*,*) 'rdsrc_eq: Converting to seismic moment'
            src%SMAG(ISRC) = MW2DCM(SRC%SMAG(ISRC)) 
         ENDIF
         IF (src%SRCTYP(ISRC)(1:1) == 'S') THEN
            src%NSRC_SRF = src%NSRC_SRF + 1
            src%AOI(ISRC) = 90.D0
            IF (src%MODE(ISRC).GT.0) THEN
               WRITE(*,*) 'rdsrc_eq: Warning cant do mode > 0'
               WRITE(*,*) 'rdsrc_eq: Resetting to fundamental mode'
               src%MODE(ISRC) = 0 
            ENDIF
            IF (src%SRCTYP(ISRC)(2:2) /= 'B' .AND. !both 
     ;          src%SRCTYP(ISRC)(2:2) /= 'R' .AND. !rayleigh
     ;          src%SRCTYP(ISRC)(2:2) /= 'V' .AND. !vertical
     ;          src%SRCTYP(ISRC)(2:2) /= 'L') THEN !love
               IERR = 1
               GOTO 60
               WRITE(*,*) 'rdsrc_eq: Invalid surface wave type!'
            ENDIF
         ELSE 
            src%NSRC_BDY = src%NSRC_BDY + 1
            src%MODE(ISRC) = 0
            IF (src%AOI(ISRC).LT.0.D0) THEN
               WRITE(*,*) 'rdsrc_eq: Warning setting aoi = 0.'
               src%AOI(ISRC) = 0.D0
            ENDIF
            IF (src%AOI(ISRC).GT.90.D0) THEN
               WRITE(*,*) 'rdsrc_eq: Warning setting aoi = 90.'
               src%AOI(ISRC) = 90.D0
            ENDIF
         ENDIF
    1 CONTINUE
      IF (src%NSRC_SRF > 0) WRITE(*,905) src%NSRC_SRF
      IF (src%NSRC_BDY > 0) WRITE(*,906) src%NSRC_BDY 
  905 FORMAT(' rdsrc_eq: Number of surface wave sources:',I4)
  906 FORMAT(' rdsrc_eq: Number of body wave sources:',I4)  
      CLOSE(IUNIT)
      RETURN
   60 CONTINUE 
      CLOSE(IUNIT) 
      IERR = 2
      WRITE(*,*) 'rdsrc_eq: Premature end of src file'
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*8 FUNCTION MW2DCM(MW) 
!
!     Converts Mw scale to dyne-centimeters where 
!     Mw = 2/3(log_10 M0 - 9.1) 
!
!     INPUT      MEANING
!     -----      ------- 
!     MW         moment magnitude 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     MW2DCM     seismic Moment (M0) in dyne-centimeters
!
!.... variable declarations
      REAL*8, INTENT(IN) :: MW 
      REAL*8 ARG, NM
!
!----------------------------------------------------------------------!
!
      ARG = 3.D0/2.D0*MW + 9.1D0
      NM = 10.D0**ARG             !M0 in Newton-Meters 
      MW2DCM = NM*10.D0**7        !convert to dyne-centimeters
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDFREQ(PROJNM, FRQ,IERR)
!
!     Reads number of frequencies in .freqs file.  We now set space
!     for FREQ in the routine. - B Baker March 2013 
!
!     INPUT      MEANING
!     -----      ------- 
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     CFTYPE     'B' body wave, 'S' surface wave
!     FREQ       frequencies to model (Hz)
!     IERR       error flag
!     NFREQ      number of frequencies
!
!.... variable declarations 
      CHARACTER(80), INTENT(IN) :: PROJNM
      INCLUDE 'fwd_struc.h'
      TYPE (FRQ_INFO) FRQ 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(80) FLNAME
      CHARACTER(1) CFTYPE
      REAL*8 FREQ
      INTEGER*4 ITYPE, IFREQ, JFREQ
      PARAMETER(IUNIT = 44)
!       
!----------------------------------------------------------------------!
!   
      IERR = 0
      FLNAME(1:80) = ' '
      FLNAME = TRIM(PROJNM)//'.freqs'
      FLNAME = ADJUSTL(FLNAME)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME))
      READ(IUNIT,*,END=60) !header
      READ(IUNIT,*,END=60) frq%NFREQ
      ALLOCATE(frq%FREQ(frq%NFREQ))
      ALLOCATE(frq%CFTYPE(frq%NFREQ)) 
      frq%NFREQ_SRF = 0
      frq%NFREQ_BDY = 0 
      JFREQ = 0
      DO 1 ITYPE=1,2
         IF (ITYPE.EQ.2) THEN
            READ(IUNIT,*,END=60) !header
            READ(IUNIT,*,END=60) !nfreq
         ENDIF
         DO 2 IFREQ=1,frq%NFREQ
            READ(IUNIT,*,END=60) CFTYPE,FREQ 
            IF (ITYPE.EQ.1 .AND. CFTYPE.EQ.'S') THEN
               JFREQ = JFREQ + 1
               frq%NFREQ_SRF = frq%NFREQ_SRF + 1
               frq%CFTYPE(JFREQ) = CFTYPE 
               frq%FREQ(JFREQ) = FREQ
            ELSEIF (ITYPE.EQ.2 .AND. CFTYPE == 'B') THEN
               JFREQ = JFREQ + 1
               frq%NFREQ_BDY = frq%NFREQ_BDY + 1
               frq%CFTYPE(JFREQ) = CFTYPE
               frq%FREQ(JFREQ) = FREQ
            ELSE
               IF (CFTYPE.NE.'B' .AND. CFTYPE.NE.'S') THEN
                  IERR = 1 
                  WRITE(*,*) 'rdfreq: Cant determine frequency type:',
     ;                       CFTYPE 
                  GOTO 10 
               ENDIF
            ENDIF !end check on type
    2    CONTINUE !loop on frequencies
         REWIND(IUNIT) 
    1 CONTINUE !loop on types
   10 CONTINUE
      IF (frq%NFREQ_BDY + frq%NFREQ_SRF.NE.frq%NFREQ) THEN
         WRITE(*,*) 'rdfreq: Error reading frequency file1'
         IERR = 1
      ENDIF
      IF (IERR.NE.0) THEN
         IF (frq%NFREQ_SRF.GT.0) THEN 
            WRITE(*,*) 'rdfreq: Number of surface wave frequencies:',
     ;                 frq%NFREQ_SRF
         ENDIF
         IF (frq%NFREQ_BDY.GT.0) THEN
            WRITE(*,*) 'rdfreq: Number of body wave frequencies:',
     ;                 frq%NFREQ_BDY
         ENDIF
         WRITE(*,*) 'rdfreq: Number of frequencies:',frq%NFREQ
      ENDIF
      CLOSE(IUNIT)
      RETURN
   60 CONTINUE
      IERR = 2
      WRITE(*,*) 'rdfreq: Premature end of frequency file'
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDFREQI_BLHD(PROJNM,NBLOCKS,IERR) 
!
!     Reads the inversion frequency block table headers
!
!     INPUT      MEANING
!     -----      ------- 
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     NBLOCKS    number of blocks in inversion
!
!.... variable declarations 
      CHARACTER(80), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(OUT) :: NBLOCKS,IERR
      CHARACTER(80) FLNAME
      LOGICAL*4 LEX 
      PARAMETER(IUNIT = 44) 
!       
!----------------------------------------------------------------------!
!   
      IERR = 0 
      FLNAME(1:80) = ' ' 
      FLNAME = TRIM(PROJNM)//'.freqinv'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX) 
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rdfreqi_blhd: Error cannot locate freqinv file'
         IERR = 1 
         RETURN
      ENDIF
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME)) 
      READ(IUNIT,*,END=60) !header
      READ(IUNIT,*,END=60) NBLOCKS
      CLOSE(IUNIT)
      RETURN
   60 CONTINUE
      IERR = 2 
      WRITE(*,*) 'rdfreqi_blhd: Premature end of frequency file'
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDFREQI_FHD(PROJNM,IBLOCK, NOMINV,IERR) 
!
!     Reads the inversion frequency frequency block header 
!
!     INPUT      MEANING
!     -----      ------- 
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     NOMINV     number of frequencies to invert for in this block
!
!.... variable declarations 
      CHARACTER(80), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(OUT) :: NOMINV,IERR
      CHARACTER(80) FLNAME
      LOGICAL*4 LEX 
      PARAMETER(IUNIT = 44) 
!       
!----------------------------------------------------------------------!
!   
      IERR = 0 
      FLNAME(1:80) = ' ' 
      FLNAME = TRIM(PROJNM)//'.freqinv'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX) 
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rdfreqi_fhd: Error cannot locate freqinv file'
         IERR = 1 
         RETURN
      ENDIF
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME)) 
      READ(IUNIT,*,END=60) !header
      READ(IUNIT,*,END=60) NBLOCKS
      DO 1 JBLOCK=1,NBLOCKS
         READ(IUNIT,*,END=60) !header
         READ(IUNIT,*,END=60) NOMINV 
         IF (JBLOCK.EQ.IBLOCK) THEN
            CLOSE(IUNIT)
            RETURN 
         ENDIF
         DO 2 IFREQ=1,NOMINV
            READ(IUNIT,*,END=60)
    2    CONTINUE 
    1 CONTINUE  
      WRITE(*,*) 'rdfreqi_fhd: Could not find frequency block!'
      IERR = 2 
      CLOSE(IUNIT)
   60 CONTINUE
      IERR = 3
      WRITE(*,*) 'rdfreqi_fhd: Premature end of frequency file'
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RD_JFREQ_INV_HD(PROJNM, NBLOCKS,IERR)
!
!     Reads the joint frequency inversion list header
!
!     INPUT      MEANING
!     -----      ------- 
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     NBLOCKS    number of inversion blocks
!
!.... variable declarations
      CHARACTER(80), INTENT(IN) :: PROJNM 
      INTEGER*4, INTENT(OUT) :: NBLOCKS, IERR 
      CHARACTER(80) FLNAME
      LOGICAL*4 LEX 
      PARAMETER(IUNIT = 55) 
!
!----------------------------------------------------------------------!
!
      IERR = 0 
      FLNAME(1:80) = ' ' 
      FLNAME = TRIM(ADJUSTL(PROJNM))//'.jfinv'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX) 
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rd_jfreq_inv_hd: Error .jfinv file does not exist!'
         IERR = 1 
         RETURN
      ENDIF
!
!.... open and read
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME))
      READ(IUNIT,*,END=60) !header
      READ(IUNIT,*,END=60) NBLOCKS
      CLOSE(IUNIT)
      RETURN
   60 CLOSE(IUNIT) 
      WRITE(*,*) 'rd_jfreq_inv_hd: Premature end of file!' 
      IERR = 2 
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RD_JFREQ_INV(PROJNM,IBLOCK,LWNDO_SRF,LWNDO_BDY,  
     ;                        FRQ,IERR) 
!
!     Reads the joint frequency inversion list
!
!     INPUT      MEANING
!     -----      -------
!     IBLOCK     block number to read
!     PROJNM     project name
!     
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     CFREQI     'B' body wave inversion frequency, 'S' surface wave
!     FREQINV    inversion frequency list
!     LINVF      True -> frequency is an inversion frequency
!     LWNDO_BDY  True -> windowing body waves
!     LWNDO_SRF  True -> windowing surface waves
!     NFINV      number of inversion frequencies 
!
!.... variable declarations
      INCLUDE 'fwd_struc.h'
      TYPE (FRQ_INFO) FRQ 
      CHARACTER(*), INTENT(IN) :: PROJNM 
      INTEGER*4, INTENT(IN) :: IBLOCK 
      LOGICAL*4, INTENT(INOUT) :: LWNDO_SRF, LWNDO_BDY 
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CFTYPE(:)
      REAL*8, ALLOCATABLE :: FREQ(:) 
      LOGICAL*4, ALLOCATABLE :: LINVB(:) 
      CHARACTER(80) FLNAME
      CHARACTER(1) CF_IN
      REAL*8 F 
      LOGICAL*4 LEX, LINVF, LSAVE_SRF, LSAVE_BDY, LINIT, LRETRY
      PARAMETER(IUNIT = 55) 
!
!----------------------------------------------------------------------!
!
      IERR = 0 
      FLNAME(1:80) = ' ' 
      FLNAME = TRIM(ADJUSTL(PROJNM))//'.jfinv'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX) 
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rd_jfreq_inv: Error .jfinv file does not exist!'
         IERR = 1 
         RETURN
      ENDIF
!
!.... open and read surface waves then body waves
      ITRY = 0
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME))
  100 CONTINUE !retry
      frq%NFREQ_BDY = 0
      frq%NFREQ_SRF = 0
      frq%NFREQ_INV = 0
      frq%NFREQ_SRF_INV = 0
      frq%NFREQ_BDY_INV = 0
      IFREQ = 0
      LINIT = .FALSE.
      DO 1 ITYPE=1,2
         READ(IUNIT,*,END=60) !header
         READ(IUNIT,*,END=60) NBLOCKS
         DO 2 JBLOCK=1,NBLOCKS
            READ(IUNIT,*,END=60) !header
            READ(IUNIT,*,END=60) NFB
            DO 3 IFB=1,NFB
               IF (IBLOCK.NE.JBLOCK) THEN
                  READ(IUNIT,*,END=60) !CFTYPE,LINVF,F
               ELSE
                  IF (.NOT.LINIT) THEN
                     ALLOCATE(CFTYPE(NFB))
                     ALLOCATE(LINVB(NFB))
                     ALLOCATE(FREQ(NFB)) 
                     LINIT = .TRUE.
                  ENDIF
                  LSAVE_BDY = .FALSE.
                  LSAVE_SRF = .FALSE.
                  READ(IUNIT,*,END=60) CF_IN,LINVF,F
                  IF (CF_IN.NE.'B' .AND. CF_IN.NE.'S') THEN
                     WRITE(*,*) 'rd_jfreq_inv: Cant determine wave type'
                     IERR = 3
                     GOTO 60
                  ENDIF
                  IF (CF_IN.EQ.'S' .AND. ITYPE.EQ.1) THEN
                     IF (LWNDO_SRF) THEN
                        LSAVE_SRF = .TRUE.
                     ELSE
                        IF (LINVF) LSAVE_SRF = .TRUE.
                     ENDIF
                     IF (LSAVE_SRF) THEN
                        IFREQ = IFREQ + 1
                        frq%NFREQ_SRF = frq%NFREQ_SRF + 1
                        FREQ(IFREQ) = F
                        CFTYPE(IFREQ) = CF_IN 
                        LINVB(IFREQ) = LINVF 
                        IF (LINVB(IFREQ)) THEN
                           frq%NFREQ_SRF_INV = frq%NFREQ_SRF_INV + 1
                           frq%NFREQ_INV = frq%NFREQ_INV + 1
                        ENDIF !end check on inversion frequency
                     ENDIF !end check on saving
                  ENDIF !end check on surface wave
                  IF (CF_IN.EQ.'B' .AND. ITYPE.EQ.2) THEN
                     IF (LWNDO_BDY) THEN
                        LSAVE_BDY = .TRUE.
                     ELSE
                        IF (LINVF) LSAVE_BDY = .TRUE.
                     ENDIF
                     IF (LSAVE_BDY) THEN
                        IFREQ = IFREQ + 1
                        frq%NFREQ_BDY = frq%NFREQ_BDY + 1
                        FREQ(IFREQ) = F 
                        CFTYPE(IFREQ) = CF_IN
                        LINVB(IFREQ) = LINVF 
                        print *, linvf,cftype(ifreq)
                        IF (LINVB(IFREQ)) THEN 
                           frq%NFREQ_BDY_INV = frq%NFREQ_BDY_INV + 1
                           frq%NFREQ_INV = frq%NFREQ_INV + 1
                        ENDIF !end check on inversion frequency
                     ENDIF !end check on saving 
                  ENDIF !end check on body waves
               ENDIF
    3       CONTINUE !loop on frequencies
            IF (IBLOCK.EQ.JBLOCK) GOTO 20
    2    CONTINUE !loop on frequency blocks
         WRITE(*,*) 'rd_jfreq_inv: Couldnt locaate freq block!'
         IERR = 1
         CLOSE(IUNIT)
         RETURN
   20    CONTINUE  
         REWIND(IUNIT)
    1 CONTINUE
!
!.... copy
      frq%NFREQ = frq%NFREQ_SRF + frq%NFREQ_BDY
      IF (frq%NFREQ == 0) THEN
         WRITE(*,*) 'rd_jfreq_inv: No inversion frequencies!'
         IERR = 1
         RETURN
      ENDIF
      ALLOCATE(frq%FREQ(frq%NFREQ))
      ALLOCATE(frq%LINVF(frq%NFREQ))
      ALLOCATE(frq%CFTYPE(frq%NFREQ))
      DO 4 IFREQ=1,frq%NFREQ
         frq%FREQ(IFREQ) = FREQ(IFREQ)
         frq%LINVF(IFREQ) = LINVB(IFREQ)
         frq%CFTYPE(IFREQ) = CFTYPE(IFREQ)
    4 CONTINUE 
      DEALLOCATE(FREQ)
      DEALLOCATE(LINVB)
      DEALLOCATE(CFTYPE)
      IF (frq%NFREQ_INV.NE.frq%NFREQ_BDY_INV + frq%NFREQ_SRF_INV) THEN 
         WRITE(*,*) 'rd_jfreq_inv: Inversion frequency mismatch!'
         IERR = 1
      ENDIF
!
!.... there is a switch to override having no inversion frequencies
      LRETRY = .FALSE.
      IF (frq%NFREQ_SRF_INV == 0 .AND. ITRY.LT.1) THEN
         WRITE(*,*) 'rd_jfreq_inv: Assuming you dont want surface waves'
         frq%NFREQ_SRF = 0
         IF (LWNDO_SRF) LWNDO_SRF = .FALSE.
         LRETRY = .TRUE.
         DEALLOCATE(frq%FREQ)
         DEALLOCATE(frq%LINVF) 
         DEALLOCATE(frq%CFTYPE)
      ENDIF
      IF (frq%NFREQ_BDY_INV == 0 .AND. ITRY.LT.1) THEN
         WRITE(*,*) 'rd_jfreq_inv: Assuming you dont want body waves'
         frq%NFREQ_BDY = 0
         IF (LWNDO_BDY) LWNDO_BDY = .FALSE.
         LRETRY = .TRUE.
         DEALLOCATE(frq%FREQ)
         DEALLOCATE(frq%LINVF)
         DEALLOCATE(frq%CFTYPE)
      ENDIF
!
!.... to get sizes right, just repeat the activity
      IF (LRETRY.AND.ITRY.LT.1) THEN
         ITRY = ITRY + 1
         REWIND(IUNIT) 
         GOTO 100
      ENDIF
!
!.... error handling
      IF (frq%NFREQ_INV == 0) THEN
         WRITE(*,*) 'rd_jfreq_inv: Error no inversion frequencies!'
         IERR = 1
         RETURN
      ENDIF
      IF (frq%NFREQ_SRF + frq%NFREQ_BDY == 0) THEN
         WRITE(*,*) 'rd_jfreq_inv: Error no modeling frequencies!'
         IERR = 1
         RETURN
      ENDIF 

!
!.... write some stuff
      WRITE(*,*)
     ;'rd_jfreq_inv: Number of surface wave modeling frequencies:',
     ;frq%NFREQ_SRF
      WRITE(*,*)
     ;'rd_jfreq_inv: Number of surface wave inversion frequencies:',
     ;frq%NFREQ_SRF_INV
      WRITE(*,*)
     ;'rd_jfreq_inv: Number of body wave modeling frequencies:',
     ;frq%NFREQ_BDY
      WRITE(*,*) 
     ;'rd_jfreq_inv: Number of body wave inversion frequencies:',
     ;frq%NFREQ_BDY_INV
      WRITE(*,*) 'rd_jfreq_inv: Number of modeling frequencies:',
     ;frq%NFREQ
      WRITE(*,*) 
     ;'rd_jfreq_inv: Number of inversion frequcnies:',
     ;frq%NFREQ_INV
      CLOSE(IUNIT)
      RETURN
   60 WRITE(*,*) 'rd_jfreq_inv: Error reading .jfinv file'
      IERR = 2
      CLOSE(IUNIT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !

      SUBROUTINE LINV_SB(PROJNM, IBLOCK, NSRC_SRF,NSRC_BDY, 
     ;                   LSURF,LBODY,IERR)
!
!     Determines if this is a surface wave and/or body wave inversion
      IMPLICIT NONE
      CHARACTER(80), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(IN) :: NSRC_SRF, NSRC_BDY, IBLOCK
      LOGICAL*4, INTENT(OUT) :: LSURF, LBODY
      INTEGER*4, INTENT(OUT) :: IERR
!... .local variables
      CHARACTER(80) FLNAME
      CHARACTER(1) CTYPE
      INTEGER*4 NBLOCKS, NFREQ, IUNIT, JBLOCK, IFREQ
      LOGICAL*4 LINVF
      PARAMETER(IUNIT = 65)
!
!----------------------------------------------------------------------!
!
!.... loop on frequency blocks
      IERR = 0
      LSURF = .FALSE.
      LBODY = .FALSE.
      FLNAME(1:80) = ' '  
      FLNAME = TRIM(ADJUSTL(PROJNM))//'.jfinv'
      FLNAME = ADJUSTL(FLNAME)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME))
      READ(IUNIT,*,END=60)
      READ(IUNIT,*,END=60) NBLOCKS
      DO 1 JBLOCK=1,IBLOCK !NBLOCKS
         READ(IUNIT,*,END=60)
         READ(IUNIT,*,END=60) NFREQ
         DO 2 IFREQ=1,NFREQ
            READ(IUNIT,*,END=60) CTYPE, LINVF 
            IF (IBLOCK.EQ.JBLOCK) THEN
               IF (CTYPE == 'S') THEN
                  IF (NSRC_SRF.GT.0 .AND. LINVF) LSURF = .TRUE.
               ELSEIF (CTYPE == 'B') THEN
                  IF (NSRC_BDY.GT.0 .AND. LINVF) LBODY = .TRUE. 
               ELSE
                  IERR = 1
                  WRITE(*,*) 'linv_sb: Invalid frequency type'
                  GOTO 50
               ENDIF
               IF (LSURF .AND. LBODY) GOTO 50
            ENDIF
    2    CONTINUE
    1 CONTINUE
   50 CONTINUE
      CLOSE(IUNIT)
      RETURN
   60 CONTINUE
      WRITE(*,*) 'linv_sb: Premature end of file'
      IERR = 1
      CLOSE(IUNIT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDFREQI(PROJNM,IBLOCK, FRQ,IERR)
!
!     Reads number of frequencies in .freqinv file
!
!     INPUT      MEANING
!     -----      ------- 
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     FREQ       frequencies to invert at (Hz)
!     IERR       error flag
!     NFREQ      number of frequencies to invert at
!
!.... variable declarations 
      CHARACTER(80), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(IN) :: IBLOCK 
      INCLUDE 'fwd_struc.h'
      TYPE (FRQ_INFO) FRQ 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(80) FLNAME
      PARAMETER(IUNIT = 44)
!       
!----------------------------------------------------------------------!
!   
      IERR = 0
      FLNAME(1:80) = ' '
      FLNAME = TRIM(PROJNM)//'.freqinv'
      FLNAME = ADJUSTL(FLNAME)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME))
      READ(IUNIT,*,END=60) !header
      READ(IUNIT,*,END=60) NBLOCKS 
      DO 1 JBLOCK=1,NBLOCKS
         READ(IUNIT,*,END=60) !header
         READ(IUNIT,*,END=60) NOMINV
         IF (IBLOCK.EQ.JBLOCK) THEN
            frq%NFREQ = NOMINV
            frq%NFREQ_INV = NOMINV
            ALLOCATE(frq%FREQ(frq%NFREQ)) 
         ENDIF
         DO 2 IFREQ=1,NOMINV
            IF (IBLOCK.EQ.JBLOCK) THEN
               READ(IUNIT,*,END=60) frq%FREQ(IFREQ)
            ELSE
               READ(IUNIT,*,END=60) 
            ENDIF
    2    CONTINUE
         IF (IBLOCK.EQ.JBLOCK) THEN
            CLOSE(IUNIT)
            RETURN
         ENDIF
    1 CONTINUE
      WRITE(*,*) 'rdfreqi: Error could not find block number' 
      IERR = 2
      CLOSE(IUNIT)
      RETURN
   60 CONTINUE
      IERR = 3
      WRITE(*,*) 'rdfreqi: Premature end of frequency file'
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDMESHBK(PROJNM, NELEME,NABS, MSH, IERR) 
!
!     Reads the Bielak mesh.  I now handle memory allocation through
!     pointers and structures to make modifications easy.  Quite 
!     simply, if you have a parameter you think should belong to 
!     the mesh design update 'fwd_struc.h' MESH_INFO and this routine 
!     - B. Baker March 2013
!
!     INPUT      MEANING
!     -----      ------- 
!     PROJNM     project name
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     AZMOD      model azimuth  (deg) 
!     CDOMAIN    holds elements domain: 
!                  'A' -> absorbing PML domain
!                  'E' -> Bielak domain
!                  'I' -> interior domain of interest
!     CNNP       holds the anchor node descriptor
!                  'AA' -> Purely absorbing domain node
!                  'BA' -> Bielak/absorbing domain interface
!                  'BI' -> Bielak/interior domain interface
!                  'EE' -> Purely Bielak domain node 
!                  'FB' -> Fixed bottom; 0 displacement
!                  'FL' -> Fixed left hand side; 0 displacement
!                  'FR' -> Fixed right hand side; 0 displacement 
!                  'II' -> Purely interior domain node 
!     DENS       density at anchor nodes (kg/m**3)
!     ECOEFF     elastic stiffness coefficients (Pa)
!     IENG       anchor node pointer on element
!     LISISO     True -> simulation isotropic 
!     NABS       number of elements in absorbing boundary
!     NELEM      number of elements in mesh
!     NELEME     number of elements in Bielak layer
!     NNPG       number of anchor nodes
!     QP         P quality factor
!     QS         S quality factor
!     XD         anchor node depth into PML in x (m)
!     XLATMAX    latitude of Bielak/Absorbing right boundary (degrees)
!     XLATMIN    latitude of Bielak/Absorbing left boundary (degrees)
!     XLONMAX    longitude of Bielak/Absorbing right boundary (degrees) 
!     XLONMIN    longitude of Bielak/Absorbing left boundary (degrees)
!     XLOCS      x locations of anchor nodes (m)
!     XWIDTH     x width of PML (m)
!     ZD         anchor node depth into PML in z (m)
!     ZLOCS      z locations of anchor nodes (m)
!     ZWIDTH     z width of PML (m)
! 
!.... variable declarations
      INCLUDE 'fwd_struc.h'
      CHARACTER(80), INTENT(IN) :: PROJNM
      TYPE (MESH_INFO) MSH 
      INTEGER*4, INTENT(OUT) :: NELEME,NABS, IERR 
!.... local variables
      CHARACTER(80) MESHFL, NODESFL
      REAL*8, ALLOCATABLE :: VARS(:)
      REAL*8 FAZMOD 
      LOGICAL*4 LEX
      PARAMETER(IUNIT = 44)
!
!----------------------------------------------------------------------!
!
!.... read the nodes file 
      IERR = 0 
      NODESFL(1:80) = ' ' 
      NODESFL = TRIM(PROJNM)//'.nodes'
      NODESFL = ADJUSTL(NODESFL) 
      INQUIRE(FILE=TRIM(NODESFL),EXIST=LEX) 
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rdmeshbk: Nodes files does notexist1',TRIM(NODESFL)
         IERR = 1
         RETURN
      ENDIF
      OPEN(FILE=TRIM(NODESFL),UNIT=IUNIT)
      READ(IUNIT,*,END=60) MSH%NELEM
!.... set space
      ALLOCATE(MSH%CDOMAIN(MSH%NELEM))       !domain descriptor
      ALLOCATE(MSH%IENG(NGNOD,MSH%NELEM))    !anchor node IEN pointer

      NELEME = 0 
      NABS = 0 
      DO 1 IELEM=1,MSH%NELEM
         READ(IUNIT,*,END=60) msh%CDOMAIN(IELEM), msh%IENG(1:4,IELEM)
         IF (msh%CDOMAIN(IELEM).EQ.'E') NELEME = NELEME + 1
         IF (msh%CDOMAIN(IELEM).EQ.'A') NABS = NABS + 1
    1 CONTINUE
      CLOSE(IUNIT)
!
!.... read the mesh
      MESHFL(1:80) = ' '
      MESHFL = TRIM(PROJNM)//'.mesh'
      MESHFL = ADJUSTL(MESHFL)
      INQUIRE(FILE=TRIM(MESHFL),EXIST=LEX) 
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rdmeshbk: Mesh files does not exist!',TRIM(MESHFL)
         IERR = 1
         RETURN 
      ENDIF
      OPEN(FILE=TRIM(MESHFL),UNIT=IUNIT)
      READ(IUNIT,*,END=61) msh%NNPG, INPTYP, msh%NORD, msh%IITYPE, IDISP
      READ(IUNIT,*,END=61) msh%XLATMIN,msh%XLONMIN, 
     ;                     msh%XLATMAX,msh%XLONMAX 
      NCOEFF = 2
      msh%LISISO = .TRUE.
      IF (INPTYP.GT.2) THEN
         NCOEFF = 5
         msh%LISISO = .FALSE.
      ENDIF
      IF (NCOEFF.GT.2 .AND. IDISP.GT.0) THEN
         WRITE(*,*) 'rdmeshbk: Anisotropic dispersion turned off'
         IDISP = 0
      ENDIF
      IF (IDISP.EQ.0) msh%FREQ0 = 0.D0 !turn off
      IF (IDISP.EQ.0) THEN
         msh%LDISP = .FALSE.
         msh%FREQ0 = 0.D0    !turn off
      ELSE
         msh%LDISP = .TRUE.  !dispersion is on in assembly and 1D solns
      ENDIF
      IF (msh%LDISP) WRITE(*,*) 'rdmeshbk: Simulation is dispersive'
      msh%NCOEFF = NCOEFF
!.... set space
      msh%NLXI =  msh%NORD + 1
      msh%NLETA = msh%NLXI !mesh is conforming
      msh%NEN   = msh%NLXI*msh%NLETA !number of element nodes
      ALLOCATE(msh%CNNPG(msh%NNPG))          !anchor node descriptor
      ALLOCATE(msh%ECOEFF(msh%NNPG,NCOEFF))  !elastic coefficitions (Pa)
      ALLOCATE(msh%DENS(msh%NNPG))           !density (kg/m**3)
      ALLOCATE(msh%XLOCS(msh%NNPG))          !x locations (m)
      ALLOCATE(msh%ZLOCS(msh%NNPG))          !z locations (m)
      ALLOCATE(msh%XD(msh%NNPG))             !x depths into PML (m)
      ALLOCATE(msh%ZD(msh%NNPG))             !z depths into PML (m)
      ALLOCATE(msh%QP(msh%NNPG))             !P quality factor
      ALLOCATE(msh%QS(msh%NNPG))             !S quality factor
      ALLOCATE(VARS(NCOEFF))
      msh%XWIDTH = 0.D0
      msh%ZWIDTH = 0.D0
      DO 2 INPG=1,MSH%NNPG
         IF (IDISP.GT.0) THEN
            READ(IUNIT,*,END=61) msh%CNNPG(INPG), 
     ;                           msh%XLOCS(INPG), msh%ZLOCS(INPG), 
     ;                           msh%XD(INPG),    msh%ZD(INPG), 
     ;                           VARS(1:NCOEFF),  msh%DENS(INPG), 
     ;                           msh%QP(INPG),    msh%QS(INPG)  
         ELSE
            READ(IUNIT,*,END=61) msh%CNNPG(INPG), 
     ;                           msh%XLOCS(INPG), msh%ZLOCS(INPG), 
     ;                           msh%XD(INPG),    msh%ZD(INPG), 
     ;                           VARS(1:NCOEFF),  msh%DENS(INPG)
            msh%QP(INPG) = 9999.D0 
            msh%QS(INPG) = 9999.D0 !like infinite quality 
         ENDIF
         IF (INPTYP.EQ.1) THEN !(p,s) -> lambda, mu
            msh%ECOEFF(INPG,2) = VARS(2)**2*msh%DENS(INPG) 
            msH%ECOEFF(INPG,1) = VARS(1)**2*msh%DENS(INPG) 
     ;                         - 2.D0*msh%ECOEFF(INPG,2)
         ELSE 
            msh%ECOEFF(INPG,1:NCOEFF) = VARS(1:NCOEFF)
         ENDIF
         msh%XWIDTH = DMAX1(MSH%XWIDTH,MSH%XD(INPG))
         msh%ZWIDTH = DMAX1(MSH%ZWIDTH,MSH%ZD(INPG)) 
    2 CONTINUE
      CLOSE(IUNIT)
      DEALLOCATE(VARS)
!
!.... calculate model azimuth
      MSH%AZMOD = FAZMOD(msh%XLATMIN,msh%XLONMIN, 
     ;                   msh%XLATMAX,msh%XLONMAX)
      RETURN
!
!.... error handling
   60 CONTINUE
      WRITE(*,*) 'rdmeshbk: Premature end of nodesfl' 
      IERR = 11
      CLOSE(IUNIT)
      RETURN
   61 CONTINUE 
      WRITE(*,*) 'rdmeshbk: Premature end of meshfl'
      IERR = 12
      CLOSE(IUNIT) 
      IF (ALLOCATED(VARS)) DEALLOCATE(VARS)  
      RETURN 
      END
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*8 FUNCTION FAZMOD(XLATMIN,XLONMIN, XLATMAX,XLONMAX)
!
!     Calcaultes the model azimuth from start and end Bielak 
!     Absorbing boundary locations defined by lat,lon
!
!     INPUT      MEANING
!     -----      ------- 
!     XLATMAX    right Bielak/absorbing latitude (degrees)
!     XLATMIN    left Bielak/absorbing latitude (degrees)
!     XLONMIN    right Bielak/absorbing longitude (degrees) 
!     XLONMIN    left Bielak/absorbing longitude (degrees)
!
!     OUTPUT     MEANING
!     ------     ------- 
!     FAZMOD     model azimuth (degrees)
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: XLATMIN, XLONMIN, XLATMAX, XLONMAX
!.... local variables
      REAL*8 TLON,TLAT, RLON,RLAT, DELT,AZ1,AZ2, HPI, PI180, GLAT, 
     ;       DSTDEG, KMPERDEG
      INTEGER*4 IERR
      PARAMETER(HPI = 1.5707963267948966D0)     !pi/2
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
      PARAMETER(KMPERDEG = 111.069365447154D0)
!
!----------------------------------------------------------------------!
! 
!.... first try for a high accuracy calculation
      CALL VINCENTY(.FALSE., XLATMIN,XLONMIN,XLATMAX,XLONMAX,  
     ;              DSTDEG,DELT,AZ1,AZ2,IERR) 
      FAZMOD = AZ1 !take the azimuth from point to point 2
      IF (IERR.NE.0) THEN
         TLON = XLONMIN*PI180   !left is source
         TLAT = HPI - GLAT(XLATMIN*PI180) 
         RLON = XLONMAX*PI180   !right is receiver
         RLAT = HPI - GLAT(XLATMAX*PI180)
         CALL BJDAZ2(TLAT,TLON,RLAT,RLON, DELT,AZ1,AZ2,0)
         FAZMOD = AZ1/PI180 !azimuth
         WRITE(*,905)  DELT/PI180*KMPERDEG, FAZMOD
      ELSE
         WRITE(*,905) DELT,AZ1 
      ENDIF
  905 FORMAT(/,' fazmod: Model distance: ',F9.3,' (km)',/, 
     ;         '         Model azimuth:  ',F9.3,' (degrees)',/)
      RETURN
      END 
