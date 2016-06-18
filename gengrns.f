      SUBROUTINE GENGRNS(PROJNM,SRCTYP, MDIM,NDOF,NNPE, 
     ;                   NL1D_LT,NL1D_RT, CSIDE, NDIM,ISRC, NFREQ, 
     ;                   FREQ,AOI,BAZN, XBLKL,XBLKR, XMOD0,XMOD1, 
     ;                   STF, IDOFSE, 
     ;                   VP1D_LT,VS1D_LT,RH1D_LT,Z1D_LT, 
     ;                   VP1D_RT,VS1D_RT,RH1D_RT,Z1D_RT, XLOCSE,ZLOCSE, 
     ;                   UE,IERR)
!
!     Calculates the Green's functions at nodal points in the 'E' 
!     domain.  
!
!     Also note that the convention is Z (or w) positive up in finite 
!     elements.  If the model is striking east-west, (i.e., azmod = 0),
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
!     AOI        angle of incidence degrees
!     BAZN       corrected back azimuth, degrees (see above) 
!     CSIDE      model side to take as 1D base model
!     IDOFSE     DOF pointer for nodes elements in Bielak domain
!     ISRC       source number
!     FREQ       frequency of interest (Hz)
!     MDIM       leading dimension for IDOFSE
!     NDIM       number of spatial dimensions
!     NDOF       number of degrees of freedom
!     NL1D_LT    number of points in left 1D model
!     NL1D_RT    number of points in right 1D model
!     NNPE       number of nodal points in Bielak domain
!     PROJNM     project name
!     RH1D_LT    density 1d model on left
!     RH1D_RT    density 1d model on right
!     SRCTYP     source type
!     STF        source time function to convolve
!     VP1D_LT    vp 1d model on left
!     VP1D_RT    vp 1d model on right
!     VS1D_LT    vs 1d model on left
!     VS1D_RT    vs 1d model on right
!     XLOCSE     x locations of points in Bielak domain
!     XMOD0      left Bielak/Internal boundary position in x
!     XMOD1      right Bielak/Internal boundary position in x
!     ZLOCSE     z locations of points in Bielak domain
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     UE         Greens functions at nodal points in Bielak domain 
! 
!.... variable declarations
      CHARACTER(80), INTENT(IN) :: PROJNM
      CHARACTER(2), INTENT(IN) :: SRCTYP
      CHARACTER(1), INTENT(IN) :: CSIDE
      COMPLEX*16, INTENT(IN) :: STF
      REAL*8, INTENT(IN) :: XLOCSE(NNPE), ZLOCSE(NNPE),
     ;                      VP1D_LT(NL1D_LT),VS1D_LT(NL1D_LT),
     ;                      VP1D_RT(NL1D_RT),VS1D_RT(NL1D_RT), 
     ;                      RH1D_LT(NL1D_LT),Z1D_LT(NL1D_LT), 
     ;                      RH1D_RT(NL1D_RT),Z1D_RT(NL1D_RT), 
     ;                      FREQ,AOI,BAZN, XBLKL,XBLKR, XMOD0,XMOD1
      INTEGER*4, INTENT(IN) :: IDOFSE(MDIM,*), MDIM,NDOF,NNPE,
     ;                         NL1D_LT,NL1D_RT, NDIM, NFREQ, ISRC
      COMPLEX*8, INTENT(OUT) :: UE(NDOF)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CLP(:) 
      COMPLEX*16, ALLOCATABLE :: USPEC(:), VSPEC(:), WSPEC(:),
     ;                           UGRN1D(:), WGRN1D(:)
      REAL*8, ALLOCATABLE :: VP1D(:), VS1D(:), RH1D(:), Z1D(:),
     ;                       ZPTS1D(:)
      REAL*8, ALLOCATABLE :: ZPTS_SHIFT(:), Z1D_SHIFT(:)
      INTEGER*4, ALLOCATABLE :: INPINI(:), INPSAV(:) 
      COMPLEX*16 CCBAZ, CSBAZ, U,V,W, XC, YC
      REAL*8 VBASE, POFF, OMEGA, PX, PY, CBAZ, SBAZ, XOFF, YOFF, ARG, 
     ;       TOL, TWOPI, PI180, VPBASE,VSBASE,RHBASE, ZMAX, ZBASE, ZTOP
      INTEGER*4 NL, NNPERD,ISRCRD,NFREQRD, NNP1D, INP1D,INPE,ISIDE,I, 
     ;          IDOF
      LOGICAL*4 LCHECK, LFLIP
      PARAMETER(TOL = 1.11D-7) 
      PARAMETER(TWOPI = 6.2831853071795862D0)
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
      PARAMETER(LFLIP = .TRUE.)
!
!----------------------------------------------------------------------!
!
!.... body wave section
      IERR = 0
      UE(1:NDOF) = CMPLX(0.D0,0.D0)
!
!.... set phase shift offset based on direction of propagation
      IF (CSIDE.EQ.'L') THEN !wave moving in +x, left is zero time 
         POFF = XMOD0
      ELSE !wave moving in -x, right is zero time
         POFF = XMOD1  
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
      IF (SRCTYP(1:1).EQ.'P' .OR. SRCTYP(1:1).EQ.'p') THEN
!
!....... force base layers to match on sides 
         IF (CSIDE.EQ.'L') THEN !wave moving +x
            VPBASE = VP1D_LT(NL1D_LT)
            VSBASE = VS1D_LT(NL1D_LT)
            RHBASE = RH1D_LT(NL1D_LT)
            ZBASE  = Z1D_LT (NL1D_LT)
            ZTOP   = Z1D_LT (1)
         ELSE !wave moving -x
            VPBASE = VP1D_RT(NL1D_RT)
            VSBASE = VS1D_RT(NL1D_RT)
            RHBASE = RH1D_RT(NL1D_RT)
            ZBASE  = Z1D_RT (NL1D_RT)
            ZTOP   = Z1D_RT (1)
         ENDIF
!
!....... classify the points
         ALLOCATE(CLP(NNPE)) 
         DO 1 INPE=1,NNPE
            IF (CSIDE.EQ.'L') THEN
               IF (XLOCSE(INPE).GE.XBLKR - TOL) THEN 
                  CLP(INPE) = 'R' 
               ELSE
                  CLP(INPE) = 'L' 
               ENDIF
            ELSE
               IF (XLOCSE(INPE).LE.XBLKL + TOL) THEN
                  CLP(INPE) = 'L' 
               ELSE
                  CLP(INPE) = 'R' 
               ENDIF 
            ENDIF
    1    CONTINUE 
!
!....... outer loop on sides
         ALLOCATE(INPINI(NNPE)) 
         INPINI(1:NNPE) = 0
         ALLOCATE(INPSAV(NNPE))
         ALLOCATE(ZPTS1D(NNPE))
         DO 2 ISIDE=1,2
            IF (ISIDE.EQ.1) THEN !set model to left 1D model
               NL = NL1D_LT
               ALLOCATE(VP1D(NL)) 
               ALLOCATE(VS1D(NL)) 
               ALLOCATE(RH1D(NL))
               ALLOCATE(Z1D (NL)) 
               VP1D(1:NL) = VP1D_LT(1:NL)
               VS1D(1:NL) = VS1D_LT(1:NL)
               RH1D(1:NL) = RH1D_LT(1:NL) 
               Z1D (1:NL) = Z1D_LT (1:NL) 
            ELSE !set model to right 1D model
               NL = NL1D_RT 
               ALLOCATE(VP1D(NL))
               ALLOCATE(VS1D(NL))
               ALLOCATE(RH1D(NL))
               ALLOCATE(Z1D (NL)) 
               VP1D(1:NL) = VP1D_RT(1:NL)
               VS1D(1:NL) = VS1D_RT(1:NL)
               RH1D(1:NL) = RH1D_RT(1:NL)
               Z1D (1:NL) = Z1D_RT (1:NL)
            ENDIF
            VP1D(NL) = VPBASE !force base to be consistent
            VS1D(NL) = VSBASE 
            RH1D(NL) = RHBASE  
!
!.......... set the Greens fns at unique points
            ZMAX =-HUGE(1.D0)
            ZPTS1D(1:NNPE) = 0.D0
            INP1D = 0
            DO 3 INPE=1,NNPE
               LCHECK = .FALSE.
               IF (ISIDE.EQ.1 .AND. CLP(INPE).EQ.'L') LCHECK = .TRUE. 
               IF (ISIDE.EQ.2 .AND. CLP(INPE).EQ.'R') LCHECK = .TRUE.
               IF (LCHECK) THEN 
                  IF (INPE.LT.NNPE) THEN
                     IF (ZLOCSE(INPE).NE.ZLOCSE(INPE+1) .OR. 
     ;                   CLP(INPE)   .NE.CLP   (INPE+1)) THEN
                        INP1D = INP1D + 1
                        ZPTS1D(INP1D) = ZLOCSE(INPE)
                        INPSAV(INP1D) = INPE
                     ENDIF
                  ELSE
                     IF (ZLOCSE(NNPE).NE.ZLOCSE(NNPE-1) .OR.
     ;                   CLP(NNPE)   .NE.CLP   (NNPE-1)) THEN
                        INP1D = INP1D + 1
                        ZPTS1D(INP1D) = ZLOCSE(NNPE)
                        INPSAV(INP1D) = NNPE
                     ENDIF
                  ENDIF
                  ZMAX = DMAX1(ZMAX,ZLOCSE(INPE))
               ENDIF
    3       CONTINUE
!
!.......... i can miss the last point
            IF (DABS(ZPTS1D(INP1D) - ZMAX).GT.TOL .AND. 
     ;          ZMAX.NE.-HUGE(1.D0)) THEN
               INP1D = INP1D + 1
               ZPTS1D(INP1D) = ZMAX
            ENDIF
 
!
!.......... set the Greens functions at unique points
            NNP1D = INP1D
            ALLOCATE(UGRN1D(NNP1D)) 
            ALLOCATE(WGRN1D(NNP1D)) 
!
!.......... P body wave
            IF (SRCTYP(2:2).EQ.'P' .OR. SRCTYP(2:2).EQ.'p') THEN
               VBASE = VP1D(NL)
               ALLOCATE(Z1D_SHIFT(NL))
               ALLOCATE(ZPTS_SHIFT(NNP1D))
               do i=1,nl
                  Z1D_SHIFT(i) = Z1D(i) + zbase 
               enddo
               do i=1,nnp1d
                  ZPTS_SHIFT(i) = ZPTS_SHIFT(i)  + zbase
               enddo
               CALL HASKGRN(NL,NNP1D,SRCTYP,LFLIP, FREQ,AOI, Z1D,VP1D,
     ;                  VS1D,RH1D,ZPTS1D(1:nnp1d), UGRN1D,WGRN1D, IERR)
!              CALL HASKGRN(NL,NNP1D,SRCTYP,.false.,FREQ,AOI, Z1D_SHIFT,
!    ;             VP1D,VS1D,RH1D,ZPTS_SHIFT,      UGRN1D,WGRN1D, IERR)
               DEALLOCATE(Z1D_SHIFT)
               DEALLOCATE(ZPTS_SHIFT)
               IF (IERR.NE.0) THEN
                  WRITE(*,*) 'gengrns: Error calling haskgrn'
                  GOTO 50
               ENDIF
               
            ELSE
               VBASE = VS1D(NL)
               IERR = 1
               WRITE(*,*) 'gengrns: Not programmed yet!'
            ENDIF
            PX = DSIN(AOI*PI180)*DCOS(BAZN*PI180)/VBASE
            PY = DSIN(AOI*PI180)*DSIN(BAZN*PI180)/VBASE 
!
!.......... expand the greens functions
            INP1D = 1
            DO 11 INPE=1,NNPE
               LCHECK = .FALSE.
               IF (ISIDE.EQ.1 .AND. CLP(INPE).EQ.'L') LCHECK = .TRUE. 
               IF (ISIDE.EQ.2 .AND. CLP(INPE).EQ.'R') LCHECK = .TRUE.
               IF (LCHECK) THEN 
!
!................ try to reuse values
c                 IF (DABS(ZLOCSE(INPE) - ZPTS1D(INP1D)).GT.1.D-8) THEN
c                    IF (INP1D.LT.NNP1D) THEN
c                       IF (DABS(ZLOCSE(INPE) - ZPTS1D(INP1D+1)).LT.
c    ;                      1.D-8) THEN
c                          INP1D = INP1D + 1
c                       ENDIF
c                    ELSE
                        INP1D = LOCATE8(NNP1D,TOL,ZLOCSE(INPE),ZPTS1D)
c                    ENDIF
c                 ENDIF
                  IF (INPINI(INPE).EQ.1) THEN
                     WRITE(*,*) 'gengrns: Warning point initialized'
                  ENDIF
                  INPINI(INPE) = 1
                  XOFF = XLOCSE(INPE) - POFF
                  YOFF = 0.D0 !we say that  y is in the plane
                  !this is a wave number vector where omega*p = k
                  !haskell has given us the response at depths  
                  !so we don't need z
                  ARG =-OMEGA*(XOFF*PX + YOFF*PY) !-omega*p.x
                  U = UGRN1D(INP1D)*CCBAZ !u component, pure radial
                  V = UGRN1D(INP1D)*CSBAZ !v component, pure transverse
!                 u = ugrn1d(inp1d)
!                 v = dcmplx(0.d0,0.d0)
                  W = WGRN1D(INP1D)       !w component, pure vertical
                  U = U*CDEXP(DCMPLX(0.D0,ARG))*STF!phase shift/convolve
                  V = V*CDEXP(DCMPLX(0.D0,ARG))*STF!phase shift/convolve
                  W = W*CDEXP(DCMPLX(0.D0,ARG))*STF!phase shift/convolve 
                  DO 12 I=1,NDIM
                     IDOF = IDOFSE(I,INPE)
                     IF (IDOF.GT.0) THEN
                      !if (i.eq.1) write(45,*) xlocse(inpe),zlocse(inpe)
                        IF (I.EQ.1) THEN
                           IF (SRCTYP(2:2).EQ.'P' .OR.
     ;                         SRCTYP(2:2).EQ.'p' .OR.
     ;                         SRCTYP(2:2).EQ.'S' .OR.
     ;                         SRCTYP(2:2).EQ.'s') THEN!u=ugrn*cos(bazn)
                              UE(IDOF) = CMPLX(U)
                           ELSE !sh; u =-ugrn*sin(bazn)
                              UE(IDOF) =-CMPLX(V)
                           ENDIF
                        ELSEIF (I.EQ.2) THEN
                           IF (SRCTYP(2:2).EQ.'P' .OR.
     ;                         SRCTYP(2:2).EQ.'p' .OR.
     ;                         SRCTYP(2:2).EQ.'S' .OR.
     ;                         SRCTYP(2:2).EQ.'s') THEN!v=ugrn*sin(bazn)
                              UE(IDOF) = CMPLX(V)
                           ELSE !sh; v = ugrn*cos(bazn)
                              UE(IDOF) = CMPLX(U)
                           ENDIF
                        ELSE
                           UE(IDOF) =-CMPLX(W)!point up now
                        ENDIF !end check on component
                     ENDIF !ned check on DOF
   12             CONTINUE !Loop on spatial dimensions
               ENDIF !end check on appropriate side
   11       CONTINUE !loop on nodal points
!.......... clean for next pass
            DEALLOCATE(UGRN1D)
            DEALLOCATE(WGRN1D) 
            DEALLOCATE(VP1D)
            DEALLOCATE(VS1D)
            DEALLOCATE(RH1D)
            DEALLOCATE(Z1D) 
    2     CONTINUE !loop on sides
!
!........ check if I got them all
          DO 15 INPE=1,NNPE
             IF (INPINI(INPE).EQ.0) THEN
                WRITE(*,*) 'gengrns: Error initializing point',INPE
                IERR = 1
             ENDIF
   15     CONTINUE 
          DEALLOCATE(INPINI)
          DEALLOCATE(INPSAV)
          DEALLOCATE(CLP)

!.... surface wave section
      ELSEIF (SRCTYP(1:1).EQ.'S' .OR. SRCTYP(1:1).EQ.'s') THEN
         ALLOCATE(USPEC(NNPE))
         ALLOCATE(VSPEC(NNPE))
         ALLOCATE(WSPEC(NNPE))
         CALL READSF_HD(PROJNM,ISRC, NNPERD,ISRCRD,NFREQRD)
         IF (NFREQRD.NE.NFREQ)
     ;   WRITE(*,*) 'gengrns: Warning nfreq != nfreq_rd'
         IF (ISRCRD .NE.ISRC )
     ;   WRITE(*,*) 'gengrns: Warning isrc != isrc_rd'
         IF (NNPERD .NE.NNPE)
     ;   WRITE(*,*) 'gengrns: Warning nnpe != nnpe_rd'
         CALL READSF(PROJNM, NNPE,NFREQ, ISRC, FREQ,
     ;               USPEC,VSPEC,WSPEC, IERR)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'gengrns: Error reading Greens functions'
            GOTO 50
         ENDIF
!
!....... tag the dofs
         DO 41 INPE=1,NNPE
!
!.......... modify greens functions
            IF (SRCTYP(2:2).EQ.'R'.OR.SRCTYP(2:2).EQ.'r') THEN
               VSPEC(INPE) = DCMPLX(0.D0,0.D0) !rayleigh only
            ELSEIF (SRCTYP(2:2).EQ.'L'.OR.SRCTYP(2:2).EQ.'l') THEN
               USPEC(INPE) = DCMPLX(0.D0,0.D0) !love only
               WSPEC(INPE) = DCMPLX(0.D0,0.D0)
            ELSEIF (SRCTYP(2:2).EQ.'V'.OR.SRCTYP(2:2).EQ.'v') THEN
               USPEC(INPE) = DCMPLX(0.D0,0.D0) !vertical onl
               VSPEC(INPE) = DCMPLX(0.D0,0.D0) !vertical only
            ENDIF
!
!.......... rotate into plane
            XC = USPEC(INPE)*CCBAZ - VSPEC(INPE)*CSBAZ
            YC = USPEC(INPE)*CSBAZ + VSPEC(INPE)*CCBAZ
            DO 42 I=1,NDIM
               IDOF = IDOFSE(I,INPE)
               IF (IDOF.GT.0) THEN 
                  IF (I.EQ.1) THEN
                     UE(IDOF) = CMPLX(XC*STF) !convolve STF
                  ELSEIF (I.EQ.2) THEN 
                     UE(IDOF) = CMPLX(YC*STF)
                  ELSE
                     UE(IDOF) = CMPLX(WSPEC(INPE)*STF)
                  ENDIF
               ENDIF !end check on dof 
   42       CONTINUE !Loop on spatial dimensions
   41    CONTINUE !loop on nodal points
         DEALLOCATE(USPEC)
         DEALLOCATE(VSPEC)
         DEALLOCATE(WSPEC)
      ELSE
         IERR = 1
         WRITE(*,*) 'gengrns: Error cannot determine source type'
      ENDIF
   50 CONTINUE !break ahead for an error
      RETURN
      END
