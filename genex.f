      SUBROUTINE GENEX(TMPDIR,LNSRF, 
     ;                 SRCTYP, MDIM, NL1D_LT,NL1D_RT, 
     ;                 NREC,CSIDE,
     ;                 ISRC,MYID,LNEW, MODE,  
     ;                 FREQ0,FREQ,AOI,BAZN, XMOD0,XMOD1, 
     ;                 XMLAT0,XMLON0, XMLAT1,XMLON1,
     ;                 SLAT,SLON,SDEP,SMAG, 
     ;                 STRIKE,DIP,RAKE,
     ;                 XREC,YREC, 
     ;                 VP1D_LT,VS1D_LT,RH1D_LT,Z1D_LT,
     ;                 VP1D_RT,VS1D_RT,RH1D_RT,Z1D_RT,
     ;                 VPD_RLLT,VSD_RLLT,ROD_RLLT,HDD_RLLT,          
     ;                 VPD_LVLT,VSD_LVLT,ROD_LVLT,HDD_LVLT,        
     ;                 VPD_RLRT,VSD_RLRT,ROD_RLRT,HDD_RLRT,      
     ;                 VPD_LVRT,VSD_LVRT,ROD_LVRT,HDD_LVRT,
     ;                 QP1D_LT,QP1D_RT, QS1D_LT,QS1D_RT, 
     ;                 STF, PERIOD,CCRAY,CCLOV, 
     ;                 EXACT,IERR)
!
!     Calculates 1D response at receivers
!
!     INPUT      MEANING
!     -----      ------- 
!     AOI        angle of incidence degrees
!     BAZN       corrected back azimuth
!     CSIDE      model side to take as 1D base model
!     FREQ       frequency of interest (Hz)
!     FREQ0      reference frequency
!     MDIM       leading dimension for IDOFSE
!     MODE       if modeling surface waves the mode
!     NDOF       number of degrees of freedom
!     NL1D_LT    number of points in left 1D model
!     NL1D_RT    number of points in right 1D model
!     NREC       number of receivers
!     RH1D_LT    density 1d model on left
!     RH1D_RT    density 1d model on right
!     SRCTYP     source type
!     STF        source time function to convolve
!     VP1D_LT    vp 1d model on left
!     VP1D_RT    vp 1d model on right
!     VS1D_LT    vs 1d model on left
!     VS1D_RT    vs 1d model on right
!     XMOD0      left Bielak/Internal boundary position in x
!     XMOD1      right Bielak/Internal boundary position in x
!     XREC       x receiver positions
!     YREC       y receiver positions
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     EXACT      1D response at receivers for this frequency  
! 
!.... variable declarations
      implicit none
      CHARACTER(*), INTENT(IN) :: TMPDIR
      CHARACTER(2), INTENT(IN) :: SRCTYP
      CHARACTER(1), INTENT(IN) :: CSIDE
      COMPLEX*8, INTENT(IN) :: STF
      REAL*8, INTENT(IN) :: VP1D_LT(NL1D_LT),VS1D_LT(NL1D_LT),
     ;                      VP1D_RT(NL1D_RT),VS1D_RT(NL1D_RT), 
     ;                      RH1D_LT(NL1D_LT),Z1D_LT(NL1D_LT), 
     ;                      RH1D_RT(NL1D_RT),Z1D_RT(NL1D_RT), 
     ;                      QP1D_RT(NL1D_RT),QP1D_LT(NL1D_LT), 
     ;                      QS1D_RT(NL1D_RT),QS1D_LT(NL1D_LT), 
     ;                      VPD_RLLT(NL1D_LT),VSD_RLLT(NL1D_LT), 
     ;                      ROD_RLLT(NL1D_LT),HDD_RLLT(NL1D_LT),  
     ;                      VPD_LVLT(NL1D_LT),VSD_LVLT(NL1D_LT),
     ;                      ROD_LVLT(NL1D_LT),HDD_LVLT(NL1D_LT), 
     ;                      VPD_RLRT(NL1D_RT),VSD_RLRT(NL1D_RT), 
     ;                      ROD_RLRT(NL1D_RT),HDD_RLRT(NL1D_RT), 
     ;                      VPD_LVRT(NL1D_RT),VSD_LVRT(NL1D_RT), 
     ;                      ROD_LVRT(NL1D_RT),HDD_LVRT(NL1D_RT),  
     ;                      XREC(NREC), YREC(NREC), 
     ;                      FREQ0,FREQ,AOI,BAZN, XMOD0,XMOD1, 
     ;                      XMLAT0,XMLON0, XMLAT1,XMLON1, 
     ;                      SLAT,SLON,SDEP,SMAG, STRIKE,DIP,RAKE
      INTEGER*4, INTENT(IN) :: MDIM, NL1D_LT,NL1D_RT, NREC, MYID, ISRC, 
     ;                         MODE 
      LOGICAL*4, INTENT(IN) :: LNEW, LNSRF
      COMPLEX*8, INTENT(OUT) :: EXACT(MDIM,NREC)
      REAL*8, INTENT(OUT) :: PERIOD(20),CCRAY(20),CCLOV(20)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(80) DISPFL, PHASEFL, RAYFL, LOVFL, PHASE_RAY, PHASE_LOV 
      CHARACTER(5) CID 
      COMPLEX*16, ALLOCATABLE ::  UGRN1D(:), VGRN1D(:), WGRN1D(:)
      REAL*8, ALLOCATABLE :: VP1D(:), VS1D(:), RH1D(:), Z1D(:)
      REAL*8, ALLOCATABLE :: VPD_RL(:), VSD_RL(:), ROD_RL(:), HDD_RL(:),
     ;                       VPD_LV(:), VSD_LV(:), ROD_LV(:), HDD_LV(:),
     ;                       QA(:), QB(:), XDIST(:), AZIMS(:) 
      COMPLEX*16 CCBAZ, CSBAZ, DPH, U, V, W, STF16 
      REAL*8 PERRAYW(20),PERLOVW(20),CCRAYW(20),CCLOVW(20), 
     ;       VBASE, POFF, DTT, OMEGA, PX, CBAZ, SBAZ, XOFF, YOFF,
     ;       TOL, TWOPI, PI180, VPBASE,VSBASE,RHBASE, PY,  
     ;       QPP, QSS, C1, C2, TREF, SSEC, ZKM, T0, ARG  
      INTEGER*4 IEVTIME(4), NSRCIN, NL, NBRAN1, NBRAN2, MDMIN, MDMAX,
     ;          IEFL   
      INTEGER*4 NF,JCOM,IREC,NXPTS,NBRAN2P,IZ,IORD,IORD1 
      LOGICAL*4 LPRATT,LBIN,LFLIP,LDEL,LEX 
      PARAMETER(TOL = 1.11D-7) 
      PARAMETER(TWOPI = 6.2831853071795862D0)
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
      PARAMETER(QPP = 600.D0, QSS = 300.D0) !Qp and Qs quality factors 
      PARAMETER(LFLIP = .TRUE.)
      PARAMETER(LPRATT = .FALSE.) !don't use pratt's convention on F.T.
      PARAMETER(LBIN = .TRUE.) !use binary files or not
      PARAMETER(NBRAN1 = 0) !min mode for earthsr
      PARAMETER(MDMIN = 0) !min mode for srgram
      PARAMETER(NSRCIN = 1) !only one source
      PARAMETER(C1  = 0.D0  , C2 = 0.D0) !min and max phase velocities
      PARAMETER(IEFL = 0)     ! = 1 apply earth flattening correction,
                              ! = 0 do not, pry use this b/c mesh is 
                              ! already flat
      PARAMETER(T0   = 0.D0)  !start time is 0 in srgramf
      !PARAMETER(TREF = 10.D0) !ref period for material disp correction
      DATA IEVTIME/2008,1,0,0/ !year, julian day, hour minute
      DATA SSEC/0.D0/
!
!----------------------------------------------------------------------!
!
!
!.... set phase shift offset based on direction of propagation
      IERR = 0
      CID(1:5) = ' '
      WRITE(CID,'(I5)') MYID 
      CID = ADJUSTL(CID)
      VPBASE = 0.D0
      VSBASE = 0.D0
      RHBASE = 0.D0
      IF (CSIDE.EQ.'L') THEN !wave moving in +x, left is zero time 
         POFF = XMOD0
         NL = NL1D_LT
         IF (SRCTYP(1:1).EQ.'P') THEN !body wave
            ALLOCATE(VP1D(NL))
            ALLOCATE(VS1D(NL))
            ALLOCATE(RH1D(NL))
            ALLOCATE(Z1D (NL))
            ALLOCATE(QA(NL))
            ALLOCATE(QB(NL))
            VP1D(1:NL) = VP1D_LT(1:NL)
            VS1D(1:NL) = VS1D_LT(1:NL)
            RH1D(1:NL) = RH1D_LT(1:NL)
            Z1D (1:NL) =  Z1D_LT(1:NL)
            QA(1:NL)   = QP1D_LT(1:NL) !QPP 
            QB(1:NL)   = QS1D_LT(1:NL) !QSS
            VPBASE = VP1D(NL)
            VSBASE = VS1D(NL)
            RHBASE = RH1D(NL)
         ELSE !surface wave 
            ALLOCATE(VPD_RL(NL))
            ALLOCATE(VSD_RL(NL))
            ALLOCATE(ROD_RL(NL))
            ALLOCATE(HDD_RL(NL)) 
            ALLOCATE(VPD_LV(NL))
            ALLOCATE(VSD_LV(NL))
            ALLOCATE(ROD_LV(NL))
            ALLOCATE(HDD_LV(NL))
            ALLOCATE(QA(NL))
            ALLOCATE(QB(NL))
            ALLOCATE(Z1D(NL))
            VPD_RL(1:NL) = VPD_RLLT(1:NL)
            VSD_RL(1:NL) = VSD_RLLT(1:NL)
            ROD_RL(1:NL) = ROD_RLLT(1:NL)
            HDD_RL(1:NL) = HDD_RLLT(1:NL)
            VPD_LV(1:NL) = VPD_LVLT(1:NL)
            VSD_LV(1:NL) = VSD_LVLT(1:NL)
            ROD_LV(1:NL) = ROD_LVLT(1:NL)
            HDD_LV(1:NL) = HDD_LVLT(1:NL)
            Z1D(1:NL)    =   Z1D_LT(1:NL)
            QA(1:NL)     =  QP1D_LT(1:NL) !QPP 
            QB(1:NL)     =  QS1D_LT(1:NL) !QSS 
         ENDIF
      ELSE !wave moving in -x, right is zero time
         POFF = XMOD1  
         NL = NL1D_RT
         IF (SRCTYP(1:1).EQ.'P') THEN !body wave
            ALLOCATE(VP1D(NL))
            ALLOCATE(VS1D(NL))
            ALLOCATE(RH1D(NL))
            ALLOCATE(Z1D (NL))
            ALLOCATE(QA(NL))
            ALLOCATE(QB(NL))
            VP1D(1:NL) = VP1D_RT(1:NL)
            VS1D(1:NL) = VS1D_RT(1:NL)
            RH1D(1:NL) = RH1D_RT(1:NL)
            Z1D (1:NL) =  Z1D_RT(1:NL)
            QA(1:NL)   = QP1D_RT(1:NL)  !QPP
            QB(1:NL)   = QS1D_RT(1:NL)  !QSS
            VPBASE = VP1D(NL)
            VSBASE = VS1D(NL)
            RHBASE = RH1D(NL)
         ELSE
            ALLOCATE(VPD_RL(NL))
            ALLOCATE(VSD_RL(NL))
            ALLOCATE(ROD_RL(NL))
            ALLOCATE(HDD_RL(NL))
            ALLOCATE(VPD_LV(NL))
            ALLOCATE(VSD_LV(NL))
            ALLOCATE(ROD_LV(NL))
            ALLOCATE(HDD_LV(NL))
            ALLOCATE(QA(NL))
            ALLOCATE(QB(NL))
            ALLOCATE(Z1D(NL))
            VPD_RL(1:NL) = VPD_RLRT(1:NL)
            VSD_RL(1:NL) = VSD_RLRT(1:NL)
            ROD_RL(1:NL) = ROD_RLRT(1:NL)
            HDD_RL(1:NL) = HDD_RLRT(1:NL)
            VPD_LV(1:NL) = VPD_LVRT(1:NL)
            VSD_LV(1:NL) = VSD_LVRT(1:NL)
            ROD_LV(1:NL) = ROD_LVRT(1:NL)
            HDD_LV(1:NL) = HDD_LVRT(1:NL)
            Z1D(1:NL)    =   Z1D_RT(1:NL)
            QA(1:NL)     =  QP1D_RT(1:NL)  !QPP
            QB(1:NL)     =  QS1D_RT(1:NL)  !QSS
         ENDIF
      ENDIF
      CBAZ = DCOS(BAZN*PI180) !> 0 moving left to right
      SBAZ = DSIN(BAZN*PI180)
      IF (DABS(CBAZ).LT.1.11D-15) CBAZ = 0.D0
      IF (DABS(SBAZ).LT.1.11D-15) SBAZ = 0.D0
      CCBAZ = DCMPLX(CBAZ,0.D0)
      CSBAZ = DCMPLX(SBAZ,0.D0)
      OMEGA = TWOPI*FREQ
!
!.... body waves
      STF16 = DCMPLX(STF) 
      IF (SRCTYP(1:1).EQ.'P') THEN
         ALLOCATE(UGRN1D(1))
         ALLOCATE(WGRN1D(1)) 
         IF (FREQ0.GT.0.D0) THEN
            CALL HASKATTN(NL,1,SRCTYP,LFLIP, FREQ,FREQ0,AOI,
     ;                    Z1D,VP1D,VS1D,RH1D,QA,QB, Z1D(1), 
     ;                    UGRN1D,WGRN1D, IERR)
            IF (IERR.NE.0) THEN
               WRITE(*,*) 'genex: Error in haskattn'
               EXACT(1:3,1:NREC) = CMPLX(0.0,0.0)
            ENDIF
         ELSE 
            CALL HASKGRN(NL,1,SRCTYP,LFLIP, FREQ,AOI, Z1D,VP1D,
     ;                   VS1D,RH1D,Z1D(1), UGRN1D,WGRN1D, IERR)
            IF (IERR.NE.0) THEN
               WRITE(*,*) 'genex: Error in haskgrn'
               EXACT(1:3,1:NREC) = CMPLX(0.,0.)
               GOTO 50
            ENDIF
         ENDIF
         VBASE = VPBASE
         OMEGA = TWOPI*FREQ 
         PX = DSIN(AOI*PI180)*DCOS(BAZN*PI180)/VBASE
         PY = DSIN(AOI*PI180)*DSIN(BAZN*PI180)/VBASE  
!
!....... loop on receivers and shift 
         DO 3 IREC=1,NREC
            XOFF = XREC(IREC) - POFF
            YOFF = YREC(IREC) - 0.D0
            DTT = XOFF*PX
            DPH = CDEXP(DCMPLX(0.D0,-OMEGA*DTT))
     ;           *CDEXP(DCMPLX(0.D0,+OMEGA*PY*YREC(IREC))) 
            EXACT(1,IREC) = CMPLX(UGRN1D(1)*DPH*CCBAZ)*STF
            EXACT(2,IREC) = CMPLX(UGRN1D(1)*DPH*CSBAZ)*STF
            EXACT(3,IREC) =-CMPLX(WGRN1D(1)*DPH)*STF
    3    CONTINUE
      ELSE
!
!....... set mode info 
         IF (FREQ0.GT.0.D0) THEN
            TREF = 1.D0/FREQ0
         ELSE
            TREF = 10.D0
         ENDIF 
         NBRAN2 = MODE 
         MDMAX = MODE  
!
!....... null out phase velocity info
         PERIOD(1:20) = 0.D0
         CCRAY(1:20)  = 0.D0
         CCLOV(1:20)  = 0.D0
!
!....... set phase vels, source depths, receiver depths for srgramf
         NF = 1
         NBRAN2P = NBRAN2
         JCOM = 1 !rayleigh
         LDEL = .FALSE. !i need those phase velocities, dont delete 
         IF (LFLIP) THEN 
            ZKM = MAXVAL(Z1D) - Z1D(1)
         ELSE 
            ZKM = Z1D(1)
         ENDIF
         ZKM = ZKM*1.D-3 
         DISPFL(1:80) = ' '
         PHASEFL(1:80) = ' '
         RAYFL(1:80) = ' '
         DISPFL = TRIM(ADJUSTL(TMPDIR))//'surf_ray_'//TRIM(CID)//'.out'
         RAYFL = TRIM(ADJUSTL(TMPDIR))//'ray-'//TRIM(CID)
         IF (CSIDE == 'L') THEN
            PHASEFL = TRIM(ADJUSTL(TMPDIR))//'surf_rayL_'//
     ;                TRIM(CID)//'.out'
         ELSE
            PHASEFL = TRIM(ADJUSTL(TMPDIR))//'surf_rayR_'// 
     ;                TRIM(CID)//'.out'
         ENDIF
         PHASE_RAY(1:80) = ' '
         DISPFL = ADJUSTL(DISPFL)
         PHASEFL = ADJUSTL(PHASEFL)
         RAYFL = ADJUSTL(RAYFL) 
         PHASE_RAY = PHASEFL
         IF (LNEW) THEN
            INQUIRE(FILE=TRIM(PHASEFL),EXIST=LEX)
            IF (LEX) THEN
               OPEN(UNIT=66,FILE=TRIM(PHASEFL),STATUS='OLD')
               CLOSE(66,STATUS='DELETE')
            ENDIF
         ENDIF
         CALL EARTHSR(NL,JCOM, NSRCIN,NF,IEFL,NBRAN1,NBRAN2P,C1,C2,     
     ;            TREF,ZKM,FREQ,SDEP, HDD_RL,VPD_RL,VSD_RL,ROD_RL,QA,QB,
     ;            LBIN,LDEL,DISPFL,PHASEFL,RAYFL, 
     ;            PERRAYW,PERLOVW,CCRAYW,CCLOVW, IERR)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'genex: Error in earthsr1'
            GOTO 500
         ENDIF
         !IF (LFLIP .AND. IL == NNP1D .OR. .NOT.LFLIP .AND. IL == 1) THEN
            DO 10 IORD=0,NBRAN2P
               IORD1 = IORD + 1
               IF (IORD.EQ.0) PERIOD(1) = PERRAYW(IORD1)
               CCRAY (IORD1) = CCRAYW (IORD1)
   10       CONTINUE  
         !ENDIF
!
!....... repeat for love waves
         NBRAN2P = NBRAN2 !max mode, can be altered  
         JCOM = 2 !Love
         DISPFL(1:80) = ' ' 
         LOVFL(1:80) = ' '
         PHASEFL(1:80) = ' '
         DISPFL = TRIM(ADJUSTL(TMPDIR))//'surf_lov_'//TRIM(CID)//'.out'
         LOVFL = TRIM(ADJUSTL(TMPDIR))//'lov-'//TRIM(CID)
         IF (CSIDE == 'L') THEN
            PHASEFL = TRIM(ADJUSTL(TMPDIR))//'surf_lovL_'//
     ;                TRIM(CID)//'.out'
         ELSE
            PHASEFL = TRIM(ADJUSTL(TMPDIR))//'surf_lovR_'//
     ;                TRIM(CID)//'.out'
         ENDIF
         PHASE_LOV(1:80) = ' '
         DISPFL = ADJUSTL(DISPFL)
         PHASEFL = ADJUSTL(PHASEFL)
         LOVFL = ADJUSTL(LOVFL)
         PHASE_LOV = PHASEFL
!        IF (LNEW) THEN 
             INQUIRE(FILE=TRIM(PHASEFL),EXIST=LEX)
             IF (LEX) THEN 
                OPEN(UNIT=66,FILE=TRIM(PHASEFL),STATUS='OLD')
                CLOSE(66,STATUS='DELETE')
             ENDIF
!        ENDIF
         CALL EARTHSR(NL,JCOM, NSRCIN,NF,IEFL,NBRAN1,NBRAN2P,C1,C2,    
     ;            TREF,ZKM,FREQ,SDEP, HDD_RL,VPD_RL,VSD_RL,ROD_RL,QA,QB,
     ;            LBIN,LDEL,DISPFL,PHASEFL,LOVFL, 
     ;            PERRAYW,PERLOVW,CCRAYW,CCLOVW, IERR)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'genex: Error in earthsr2'
            GOTO 500 
         ENDIF
         !IF (LFLIP .AND. IL == NNP1D .OR. .NOT.LFLIP .AND. IL == 1) THEN
            DO 11 IORD=0,NBRAN2P
               IORD1 = IORD + 1
               !PERIOD(IORD1) = PERLOVW(IORD1)
               CCLOV (IORD1) = CCLOVW (IORD1)
   11       CONTINUE
         !ENDIF

!
!....... calculate offsets to model
         NXPTS = NREC 
         ALLOCATE(XDIST(NXPTS)) 
         ALLOCATE(AZIMS(NXPTS)) 
         CALL XPTS2DIST3(NXPTS, CSIDE,SLAT,SLON, XMOD0,XMOD1,
     ;                   XMLAT0,XMLON0, XMLAT1,XMLON1, XREC,
     ;                   XDIST,AZIMS,IERR)  
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'genex: Error calculating distances to points!'
            RETURN
         ENDIF 
         ALLOCATE(UGRN1D(NREC))
         ALLOCATE(VGRN1D(NREC))
         ALLOCATE(WGRN1D(NREC)) 
!
!....... calculate the surface wave greens functions
         NF = 1 
         IZ = 1
         CALL SRGRAMF(NF,LPRATT,LBIN, IZ,NXPTS, IEVTIME(1:4),SSEC,  
     ;                TMPDIR,RAYFL,LOVFL, MDMIN,MDMAX, .TRUE.,      
     ;                SLAT,SLON,SDEP, 
     ;                STRIKE,DIP,RAKE,SMAG,T0,   
     ;                ISRC,XMLAT0,XMLON0,XMLAT1,XMLON1,     
     ;                MYID,XDIST,AZIMS,FREQ,  
     ;                NXPTS,UGRN1D,VGRN1D,WGRN1D, IERR)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'genex: Error in srgramf'
            GOTO 500
         ENDIF
         DEALLOCATE(XDIST) 
         DEALLOCATE(AZIMS) 
!
!
!        IF (LNEW) THEN
            INQUIRE(FILE=TRIM(ADJUSTL(PHASE_RAY)),EXIST=LEX)
            IF (LEX) THEN
               OPEN(UNIT=66,FILE=TRIM(ADJUSTL(PHASE_RAY)),STATUS='OLD')
               CLOSE(66,STATUS='DELETE')
            ENDIF
            INQUIRE(FILE=TRIM(PHASE_LOV),EXIST=LEX)
            IF (LEX) THEN
               OPEN(UNIT=66,FILE=TRIM(ADJUSTL(PHASE_LOV)),STATUS='OLD')
               CLOSE(66,STATUS='DELETE')
            ENDIF
!        ENDIF
         IF (LNSRF) THEN
            DO 18 IREC=1,NREC
               IF (CDABS(UGRN1D(IREC)).GT.0.D0) 
     ;         UGRN1D(IREC) = UGRN1D(IREC)/CDABS(UGRN1D(IREC))
               IF (CDABS(VGRN1D(IREC)).GT.0.D0) 
     ;         VGRN1D(IREC) = VGRN1D(IREC)/CDABS(VGRN1D(IREC))
               IF (CDABS(WGRN1D(IREC)).GT.0.D0)
     ;         WGRN1D(IREC) = WGRN1D(IREC)/CDABS(WGRN1D(IREC))
   18       CONTINUE 
         ENDIF
!
!....... convolve source time function
         CALL ZSCAL(NREC,STF16,UGRN1D,1)
         CALL ZSCAL(NREC,STF16,VGRN1D,1)
         CALL ZSCAL(NREC,STF16,WGRN1D,1)
!
!....... read the py table
         EXACT(1:3,1:NREC) = CMPLX(0.0,0.0) 
         DO 20 IREC=1,NREC
!
!.......... modify greens functions
            IF (SRCTYP(2:2) == 'R') THEN !rayleigh only
               VGRN1D(IREC) = DCMPLX(0.D0,0.D0)
            ELSEIF (SRCTYP(2:2) == 'L') THEN !love only 
               UGRN1D(IREC) = DCMPLX(0.D0,0.D0)
               WGRN1D(IREC) = DCMPLX(0.D0,0.D0)
            ELSEIF (SRCTYP(2:2) == 'V') THEN !vertical only
               UGRN1D(IREC) = DCMPLX(0.D0,0.D0) !vertical only
               VGRN1D(IREC) = DCMPLX(0.D0,0.D0) !vertical only
            ELSE
               IF (SRCTYP(2:2) /= 'B') THEN
                  WRITE(*,*) 'grns_srf: Unknown source type!'
                  IERR = 1
                  GOTO 500
               ENDIF
            ENDIF
            PX = 0.D0
            PY = 0.D0
            IF (SRCTYP(2:2) == 'L') THEN
               IF (CCLOV(1).GT.0.D0) THEN
                  PX = DCOS(BAZN*PI180)/(CCLOV(1)*1.D3)
                  PY = DSIN(BAZN*PI180)/(CCLOV(1)*1.D3)
               ENDIF
            ELSE 
               IF (CCRAY(1).GT.0.D0) THEN
                  PX = DCOS(BAZN*PI180)/(CCRAY(1)*1.D3)
                  PY = DSIN(BAZN*PI180)/(CCRAY(1)*1.D3)
               ENDIF
            ENDIF 
            XOFF = XREC(IREC) - POFF
            YOFF = YREC(IREC) - 0.D0 
            ARG =-OMEGA*(XOFF*PX + YOFF*PY) !-omega*p.x
            DPH = CDEXP(DCMPLX(0.D0,ARG)) 
            U = UGRN1D(IREC)*CCBAZ - VGRN1D(IREC)*CSBAZ
            V = UGRN1D(IREC)*CSBAZ + VGRN1D(IREC)*CCBAZ
            W = WGRN1D(IREC) !already point up 
            EXACT(1,IREC) = CMPLX(U)
            EXACT(2,IREC) = CMPLX(V)
            EXACT(3,IREC) = CMPLX(W)
   20    CONTINUE !Loop on receivers
      ENDIF 
   50 CONTINUE
  500 CONTINUE
      IF (ALLOCATED(VP1D))   DEALLOCATE(VP1D)
      IF (ALLOCATED(VS1D))   DEALLOCATE(VS1D)
      IF (ALLOCATED(RH1D))   DEALLOCATE(RH1D)
      IF (ALLOCATED(Z1D))    DEALLOCATE(Z1D ) 
      IF (ALLOCATED(VPD_RL)) DEALLOCATE(VPD_RL)
      IF (ALLOCATED(VSD_RL)) DEALLOCATE(VSD_RL)
      IF (ALLOCATED(ROD_RL)) DEALLOCATE(ROD_RL)
      IF (ALLOCATED(HDD_RL)) DEALLOCATE(HDD_RL)
      IF (ALLOCATED(VPD_LV)) DEALLOCATE(VPD_LV)
      IF (ALLOCATED(VSD_LV)) DEALLOCATE(VSD_LV)
      IF (ALLOCATED(ROD_LV)) DEALLOCATE(ROD_LV)
      IF (ALLOCATED(HDD_LV)) DEALLOCATE(HDD_LV)
      IF (ALLOCATED(QA))     DEALLOCATE(QA)
      IF (ALLOCATED(QB))     DEALLOCATE(QB) 
      IF (ALLOCATED(UGRN1D)) DEALLOCATE(UGRN1D)
      IF (ALLOCATED(VGRN1D)) DEALLOCATE(VGRN1D)
      IF (ALLOCATED(WGRN1D)) DEALLOCATE(WGRN1D) 
      RETURN
      END
