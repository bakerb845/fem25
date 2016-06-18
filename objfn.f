!----------------------------------------------------------------------!
!     A set of routines for calculating the objective function and     !
!     residuals to backpropagate                                       !
!----------------------------------------------------------------------!
      REAL*4 FUNCTION COBJ4(MDIM,MFREQ,MREC, NFREQ,NREC,NSRC, NDIM, 
     ;                      LVERB,NORM,IRESTP, EPSS,THRESHS,DTMAX, 
     ;                      FREQ, SWGHT,OBS,EST)
! 
!     Calculates the objective function 
! 
!     INPUT      MEANING
!     -----      ------- 
!     EPSS       thresholding criteria, ex: 0.2(mean|d_{obs}_i|) 
!     EST        estimate wavefield (should be rotated by now!)
!     IRESTP     =1 -> Phase only 
!                =2 -> Amplitude only 
!                =3 -> phase and amplitude (default)
!                False -> data is stored complex 
!     LVERB      True -> verbose, False -> quiet
!     MFREQ      max number of 
!     MREC       max number of receivers
!     NFREQ      number of frequencies
!     NORM       =1 ->1 norm
!                =2 ->2 norm (default)
!                =3 ->Huber norm
!                =4 ->Combination 
!     NREC       number of receivers 
!     OBS        observations in observation coordinates
!     SWGHT      observation weights
!     THRESH     threshholding criteria for l1/l2 norm 
!  
!.... variable declarations 
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,MFREQ,MREC,*), 
     ;                         EST(MDIM,MFREQ,MREC,*)
      REAL*8, INTENT(IN) :: FREQ(NFREQ)
      REAL*4, INTENT(IN) :: SWGHT(MDIM,MFREQ,MREC,*), EPSS,THRESHS,
     ;                      DTMAX
      INTEGER*4, INTENT(IN) :: MDIM,MFREQ,MREC,NFREQ,NREC,NSRC, 
     ;                         NDIM,NORM,IRESTP
      LOGICAL*4 LVERB
!.... local variables 
      COMPLEX*8 DIF,S
      REAL*4 THRESH,EPS,RESMN, EPS2,C,ARG,COBJF, ARMS, PRMS, 
     ;       DPHI, DAMP, PI180I, SPHASE 
      INTEGER*4 IFREQ,IREC,ISRC,I, NOBS
      COMPLEX*8 CRESID
      REAL*4 CRESMN
      PARAMETER(PI180I = 57.2957802) !180/pi
      REAL*4, PARAMETER :: TWOPI = 6.283185307179586
      REAL*4, PARAMETER :: PI = 3.1415926535897931
! 
!----------------------------------------------------------------------!
! 
!.... loop on sources and receivers and sum objective function
      EPS = 0.0
      THRESH = 0.0
      IF (NORM.GT.2) THEN
         RESMN = CRESMN(MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC, 
     ;                  IRESTP, OBS,EST) 
         EPS = EPSS*RESMN
         THRESH = THRESHS*RESMN 
      ENDIF
      EPS2 = EPS**2
      NOBS = 0
      COBJ4 = 0.0
      ARMS = 0.0
      PRMS = 0.0
      DO 1 IFREQ=1,NFREQ
         IF (LVERB) COBJF = 0.0
         DO 2 IREC=1,NREC
            DO 3 ISRC=1,NSRC
               DO 4 I=1,NDIM
                  IF (CABS(OBS(I,IFREQ,IREC,ISRC)).EQ.0.0) GOTO 400 
                  NOBS = NOBS + 1
                  DIF = CRESID(IRESTP,OBS(I,IFREQ,IREC,ISRC),
     ;                                EST(I,IFREQ,IREC,ISRC))
! 
!................ calculate residual 
                  S = CMPLX(SWGHT(I,IFREQ,IREC,ISRC),0.0)
                  IF (NORM.EQ.1) THEN !L1 
                     C = CABS(S*DIF)
                  ELSEIF (NORM.EQ.3) THEN !Huber
                     ARG = CABS(S*DIF)
                     IF (ARG.LE.EPS) THEN
                        C = ARG**2/(2.0*EPS)
                     ELSE
                        C = ARG - EPS/2.0
                     ENDIF
                  ELSEIF (NORM.EQ.4) THEN !L1/L2 criterion
                     ARG = CABS(DIF)
                     IF (ARG.LE.THRESH) THEN
                        C = CABS(S*DIF)**2/(2.0*EPS2)
                     ELSE
                        C = CABS(S*DIF)/EPS
                     ENDIF
                  ELSE !default to L2
                     DAMP = CABS(OBS(I,IFREQ,IREC,ISRC)) 
     ;                    - CABS(EST(I,IFREQ,IREC,ISRC)) 
                     DPHI = SPHASE(OBS(I,IFREQ,IREC,ISRC))*PI180I 
     ;                    - SPHASE(EST(I,IFREQ,IREC,ISRC))*PI180I
                     C = CABS(S)**2*REAL(CONJG(DIF)*DIF) !CABS(DIF)**2
                     ARMS = ARMS + CABS(S)**2*DAMP**2
                     PRMS = PRMS + CABS(S)**2*DPHI**2
                     IF (DTMAX > 0.0) THEN
                        DPHI = SPHASE(OBS(I,IFREQ,IREC,ISRC))
     ;                       - SPHASE(EST(I,IFREQ,IREC,ISRC))
                        IF (DPHI <-PI) DPHI = DPHI + TWOPI
                        IF (DPHI >+PI) DPHI = DPHI - TWOPI 
                        IF (ABS(DPHI/(TWOPI*FREQ(IFREQ))) > DTMAX) THEN
                          write(*,*) 'bad',dtmax,freq(ifreq),
     ;                               dphi/(twopi*freq(ifreq))
                          write(*,*) i,ifreq,irec,isrc,
     ;obs(i,ifreq,irec,isrc),est(i,ifreq,irec,isrc)
                          IF (IRESTP.EQ.1) THEN
                             C = CABS(S)**2*4.0
                          ELSE
                             C = CABS(S)**2*(DAMP)**2
                          ENDIF
                          DPHI = DPHI*PI180I
                          PRMS = PRMS - CABS(S)**2*DPHI**2
                          PRMS = PRMS + CABS(S)**2*90.0**2
                        ENDIF
                     ENDIF
                  ENDIF
                  COBJ4 = COBJ4 + C
                  IF (LVERB) COBJF = COBJF + C
  400             CONTINUE !nothing to do
    4          CONTINUE !loop on components
    3       CONTINUE !loop on sources
    2    CONTINUE !loop on receivers
         IF (LVERB) WRITE(*,900) IFREQ,COBJF
  900    FORMAT(' cobj4: Cumulative residual for frequency: ',I4,' is:',
     ;          E18.8)
    1 CONTINUE
! 
!.... very least write out objective function
      IF (NORM.EQ.1) THEN
         WRITE(*,901) COBJ4
  901    FORMAT(' cobj4: L1 Objective function is: ',E18.8)
      ELSEIF (NORM.EQ.3) THEN
         WRITE(*,902) COBJ4
  902    FORMAT(' cobj4: Huber Objective function is: ',E18.8)
      ELSEIF (NORM.EQ.4) THEN
         WRITE(*,903) COBJ4
  903    FORMAT(' cobj4: L1/L2 Objective function is: ',E18.8)
      ELSE
         IF (NOBS.EQ.0) THEN
            WRITE(*,*) 'cobj4: Warning no observations!'
         ELSE
            WRITE(*,*) 'cobj4: Number of observations:',NOBS
         ENDIF 
         IF (IRESTP.EQ.2 .OR. IRESTP.EQ.3) THEN
            IF (NOBS.GT.0) THEN
               ARMS = SQRT(ARMS/FLOAT(NOBS))
               WRITE(*,905) ARMS
  905          FORMAT(' cobj4: Amplitude RMS:',E18.8)
            ENDIF
         ENDIF
         IF (IRESTP.EQ.1 .OR. IRESTP.EQ.3) THEN 
            IF (NOBS.GT.0) THEN
               PRMS = SQRT(PRMS/FLOAT(NOBS))
               WRITE(*,906) PRMS 
  906          FORMAT(' cobj4: Phase RMS:',E18.8,' degrees')
            ENDIF
         ENDIF
         COBJ4 = COBJ4/2.0 !technically it is 1/2 |delta d|**2
         WRITE(*,907) COBJ4
  907    FORMAT(' cobj4: L2 Objective function is: ',E18.8)
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*4 FUNCTION COBJ4_COVD(MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC, 
     ;                           IRESTP,IBPHASE, LCOVD,LDWGHT, 
     ;                           DTMAX,AZMOD,DX,FREQ, WGHTS,OBS,EST, 
     ;                           IERR)
!
!     Calculates the objective function 
!       E = 1/2 (d - u)* Sigma^{-1} C_D^{-1} Sigma^{-1} (d - u)
!
!
!     INPUT     MEANING
!     -----     -------
!     AZMOD     model azimuth
!     DX        receiver distances
!     DTMAX     max phase residual (s)
!     EST       estimates
!     FREQ      frequencies (Hz)
!     IBPHASE   residual backpropagation correction (pry 0)
!     IRESTP    residual type, use 2 only for now
!     LCOVD     False -> set data covariance matrix to identity
!     LDWGHT    True -> use distances when calculating data covar
!     MDIM      leading dimension
!     MFREQ     leading dimension
!     MREC      leading dimension
!     NDIM      number of components (3)
!     NFREQ     number of frequencies
!     NREC      number of receivers
!     NSRC      number of sources
!     OBS       observations
!     WGHTS     weights
! 
!     OUTPUT    MEANING
!     ------    -------
!     COBJ      objective function
!     IERR      error flag
!   
!.... variable declarations
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,MFREQ,MREC,*),
     ;                         EST(MDIM,MFREQ,MREC,*)
      REAL*8, INTENT(IN) :: FREQ(NFREQ), DX(NREC-1), AZMOD
      REAL*4, INTENT(IN) :: WGHTS(MDIM,MFREQ,MREC,*), DTMAX
      INTEGER*4, INTENT(IN) :: MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC,
     ;                         IRESTP,IBPHASE 
      LOGICAL*4, INTENT(IN) :: LCOVD, LDWGHT
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables 
      REAL*8, ALLOCATABLE :: COVD8(:,:)
      REAL*4, ALLOCATABLE :: COVD(:,:), W1(:,:), RESIDV(:), WORK(:)
      INTEGER*4, ALLOCATABLE :: IPERM_REC(:)
      REAL*4 HALF, OBJL, DT, DPHI, DAMP, PHASE_RMS, AMP_RMS, TWOPI, 
     ;       PI180I, SPHASE, SDOT
      INTEGER*4 NOBS, NOBS2, M, LDC, IFREQ, ISRC, IREC, JREC, 
     ;          IROW, JCOL, LROW, KCOL, I, J
      PARAMETER(HALF = 0.50)
      PARAMETER(PI180I = 57.2957802) !180/pi
      PARAMETER(TWOPI = 6.283185307179586)
!
!----------------------------------------------------------------------!
!
!.... set space and calculate laplacian
      COBJ4_COVD = HUGE(1.0)
      NOBS2 = 2*NDIM*NREC
      M = NOBS2
      ALLOCATE(IPERM_REC(NDIM*NREC))
      ALLOCATE(RESIDV(NOBS2))
      ALLOCATE(WORK(NOBS2))
      ALLOCATE(W1(NDIM,NREC))
      LDC = M
      ALLOCATE(COVD8(LDC,M))
      IF (.NOT.LCOVD) THEN
!        CALL LAPLACE_1D(LDC,NREC,DX, COVD8,IERR)
!        IF (IERR.NE.0) THEN
!           IERR = 0
!           WRITE(*,*) 'cobj4_covd: Error setting Laplacian!'
!           WRITE(*,*) 'cobj4_covd: Returning to identity matrix'
!           COVD8(1:M,1:M) = 0.D0
!           DO 1 I=1,M
!              COVD8(I,I) = 1.D0
!   1      CONTINUE
!        ENDIF
!     ELSE
         COVD8(1:M,1:M) = 0.D0
         DO 1 I=1,M
            COVD8(I,I) = 1.D0
    1    CONTINUE 
      ENDIF
      ALLOCATE(COVD(LDC,M))
!     DO 3 I=1,M
!        DO 4 J=1,M
!           COVD(I,J) = SNGL(COVD8(I,J)) 
!   4    CONTINUE
!   3 CONTINUE
!
!.... calculate objective function
      NOBS = 0
      OBJL = 0.0
      PHASE_RMS = 0.0
      AMP_RMS = 0.0
      DO 11 IFREQ=1,NFREQ
         DO 12 ISRC=1,NSRC
            RESIDV(1:NOBS2) = 0.0
            IF (MAXVAL(CABS(OBS(1:NDIM,IFREQ,1:NREC,ISRC))).GT.0.0) THEN
!
!............. generate a permutation list [w obs, v obs, u obs; null obs]
               CALL PERM_OBS(NDIM,NDIM,NREC, 
     ;                       OBS(1:NDIM,IFREQ,1:NREC,ISRC), 
     ;                       IPERM_REC,IERR)
               IF (IERR /= 0) THEN 
                  WRITE(*,*) 'cobj_covd8: Error calling perm_obs!'
                  RETURN
               ENDIF
               CALL LAPLACE_1DV2(LDC,NDIM,NREC,LDWGHT, DX,  
     ;                           OBS(1:NDIM,IFREQ,1:NREC,ISRC), 
     ;                           COVD8,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'cobj_covd8: Error calling laplace_1d!'
                  DO 9 I=1,M
                     COVD8(I,I) = 1.D0
    9             CONTINUE
               ENDIF

!
!............. use kronecker delta on weights, will use covar matrix
               DO 112 IREC=1,NREC
                  DO 113 I=1,NDIM
                     IF (WGHTS(I,IFREQ,IREC,ISRC).GT.0.0) THEN
                        W1(I,IREC) = 1.0
                     ELSE
                        W1(I,IREC) = 0.0
                     ENDIF   
  113             CONTINUE
  112          CONTINUE
               CALL OE2RESV(NDIM, NOBS2, NDIM,NREC, DTMAX,  
     ;                      IRESTP,IBPHASE, AZMOD,FREQ(IFREQ),
     ;                      IPERM_REC, W1,EST(1:NDIM,IFREQ,1:NREC,ISRC),
     ;                         OBS(1:NDIM,IFREQ,1:NREC,ISRC),
     ;                      WORK,IERR) 
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'cobj4_covd: Error calling oe2resv!'
                  RETURN
               ENDIF
!
!............. scale columns; post multiplication by diagonal matrix
               I = 1
               IREC = 0
               DO 121 IROW=1,M
                  IREC = IREC + 1
                  IF (IREC.GT.NREC) THEN
                     I = I + 1
                     IREC = 1
                  ENDIF
                  IF (IROW.EQ.NDIM*NREC + 1) I = 1
                  IF (IROW > NDIM*NREC) THEN
                     LROW = NDIM*NREC + IPERM_REC( (IREC - 1)*NDIM + I) 
                  ELSE
                     LROW =             IPERM_REC( (IREC - 1)*NDIM + I)
                  ENDIF
                  J = 1
                  JREC = 0
                  DO 122 JCOL=1,M
                     JREC = JREC + 1
                     IF (JREC.GT.NREC) THEN
                        J = J + 1
                        JREC = 1
                     ENDIF
                     IF (JCOL.EQ.NDIM*NREC + 1) J = 1
                     IF (JCOL > NDIM*NREC) THEN
                        KCOL = NDIM*NREC + IPERM_REC((JREC - 1)*NDIM+J)
                     ELSE
                        KCOL =             IPERM_REC((JREC - 1)*NDIM+J)
                     ENDIF
                     COVD(LROW,KCOL) = SQRT(WGHTS(I,IFREQ,IREC,ISRC))
     ;                                *SNGL(COVD8(LROW,KCOL))
     ;                                *SQRT(WGHTS(J,IFREQ,JREC,ISRC))
  122             CONTINUE
  121          CONTINUE
!
!............. apply resid = Sigma^{-1} C_D^{-1} Sigma^{-1} (d - Au)
               CALL SGEMV('N',M,M, 1.0, COVD,LDC, WORK,1, 0.0,RESIDV,1)
               OBJL = OBJL + HALF*SDOT(NOBS2,RESIDV,1,RESIDV,1)
               DO 13 IREC=1,NREC
                  DO 14 I=1,NDIM
                     IF (CABS(OBS(I,IFREQ,IREC,ISRC)).GT.0.0) THEN 
                        NOBS = NOBS + 1 
                        DPHI = SPHASE(OBS(I,IFREQ,IREC,ISRC))
     ;                       - SPHASE(EST(I,IFREQ,IREC,ISRC))
                        DT = DPHI/(TWOPI*SNGL(FREQ(IFREQ)))
                        PHASE_RMS = PHASE_RMS 
     ;                            + (WGHTS(I,IFREQ,IREC,ISRC)*DPHI)**2
                        DAMP = CABS(OBS(I,IFREQ,IREC,ISRC)) 
     ;                       - CABS(EST(I,IFREQ,IREC,ISRC))
                        AMP_RMS = AMP_RMS 
     ;                          + (WGHTS(I,IFREQ,IREC,ISRC)*DAMP)**2
                     ENDIF
   14             CONTINUE
   13          CONTINUE
            ENDIF !end check on active frequency/source pair
   12    CONTINUE
   11 CONTINUE 
      COBJ4_COVD = OBJL
      IF (NOBS.EQ.0) THEN
         WRITE(*,*) 'cobj4_covd: Error no observations!'
         IERR = 1
         RETURN
      ENDIF
      WRITE(*,805) NOBS
      WRITE(*,900) COBJ4_COVD 
      WRITE(*,905) SQRT(PHASE_RMS/PI180I/FLOAT(NOBS))
      WRITE(*,910) SQRT(AMP_RMS/FLOAT(NOBS))
  805 FORMAT(' cobj4_covd: Number of observations:',I6)
  900 FORMAT(' cobj4_covd: L2 objective function is:',E18.8)
  905 FORMAT(' cobj4_covd: Phase RMS:',E18.8,' degrees')
  910 FORMAT(' cobj4_covd: Amplitude RMS:',E18.8) 
      DEALLOCATE(IPERM_REC)
      DEALLOCATE(COVD)
      DEALLOCATE(WORK)
      DEALLOCATE(RESIDV)
      DEALLOCATE(W1)
      DEALLOCATE(COVD8)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE UNWRAP_OEST(MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC, 
     ;                       OBS_IN,EST_IN, OBS,EST) 
!
!     Unpacks the estimates/observations in (amplitude,phase) form 
!     for a residual calculation 
!
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: EST_IN(MDIM,MFREQ,MREC,*), 
     ;                         OBS_IN(MDIM,MFREQ,MREC,*) 
      INTEGER*4, INTENT(IN) :: MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC
      COMPLEX*8, INTENT(OUT) :: EST(MDIM,MFREQ,MREC,*),
     ;                          OBS(MDIM,MFREQ,MREC,*) 

      INTEGER*4 IFREQ, IREC, ISRC, I 
!
!----------------------------------------------------------------------!
!
      DO 1 IFREQ=1,NFREQ
         DO 2 IREC=1,NREC
            DO 3 ISRC=1,NSRC
               DO 4 I=1,NDIM
                  EST(I,IFREQ,IREC,ISRC) = EST_IN(I,IFREQ,IREC,ISRC)
                  OBS(I,IFREQ,IREC,ISRC) = OBS_IN(I,IFREQ,IREC,ISRC)
    4          CONTINUE 
    3      CONTINUE !loop on sources 
    2    CONTINUE 
    1 CONTINUE 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*4 FUNCTION CRESMN(MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC, 
     ;                       IRESTP, OBS,EST)
!
!     Calculates average magnitude magnitude of residuals
!
!     INPUT      MEANING
!     -----      ------- 
!     EST        estimate wavefield (should be rotated by now!)
!     IRESTP     =1 -> Phase only 
!                =2 -> Amplitude only 
!                =3 -> phase and amplitude (default)
!     MFREQ      max number of 
!     MREC       max number of receivers
!     NFREQ      number of frequencies
!     NREC       number of receivers 
!     OBS        observations in observation coordinates
!     THRESH     threshholding criteria for l1/l2 norm 
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,MFREQ,MREC,*), 
     ;                         EST(MDIM,MFREQ,MREC,*)
      INTEGER*4, INTENT(IN) :: MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC,
     ;                         IRESTP
!.... local variables
      COMPLEX*8 DIF, CRESID
      INTEGER*4 IFREQ,IREC,ISRC,I, NOBS
!
!----------------------------------------------------------------------!
!
      NOBS = 0
      CRESMN = 0.0
      DO 1 IFREQ=1,NFREQ
         DO 2 IREC=1,NREC
            DO 3 ISRC=1,NSRC
               DO 4 I=1,NDIM
                  IF (CABS(OBS(I,IFREQ,IREC,ISRC)).EQ.0.0) GOTO 400
                  DIF = CRESID(IRESTP,OBS(I,IFREQ,IREC,ISRC),
     ;                                EST(I,IFREQ,IREC,ISRC))
                  NOBS = NOBS + 1
                  CRESMN = CRESMN + CABS(DIF)
  400             CONTINUE 
    4          CONTINUE
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
      IF (NOBS.GT.0) THEN
         CRESMN = CRESMN/FLOAT(NOBS)  
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE ROTEST(MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC, 
     ;                  AZMOD,EST) 
!
!     This routine rotates our estimate wavefield (u,v,w) to 
!     observation coordinates (North,East,Longitudnal).  Note, we 
!     already have z pointing up as a consequence of FEM.  Rotation 
!     conventions are positive counter-clockwise: 
!     http://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions
!
!     INPUT      MEANING
!     -----      ------- 
!     EST        estimtes in (u,v,w) coordinates
!     MDIM       leading dimension for est
!     MFREQ      leading dimension for est
!     MREC       leading dimension for est
!     NDIM       number of components in solution
!     NFREQ      number of frequencies 
!     NREC       number of receivers
!     NSRC       number of sources 
! 
!     OUTPUT     MEANING 
!     ------     ------- 
!     EST        estimates rotated into observation (N,E,Z) coords
!
!.... variable declarations
      COMPLEX*8, INTENT(INOUT) :: EST(MDIM,MFREQ,MREC,*)
      REAL*8, INTENT(IN) :: AZMOD
      INTEGER*4, INTENT(IN) :: MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC
      COMPLEX*8 COSAZ,SINAZ, U,V,W, NEST,EEST,ZEST
      REAL*8 PI180
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
!
!----------------------------------------------------------------------!
!
      COSAZ = CMPLX(DCOS(AZMOD*PI180),0.D0)
      SINAZ = CMPLX(DSIN(AZMOD*PI180),0.D0)
      DO 1 IFREQ=1,NFREQ
         DO 2 ISRC=1,NSRC
            DO 3 IREC=1,NREC
               IF (NDIM > 2) THEN
                  U = EST(1,IFREQ,IREC,ISRC)
                  V = EST(2,IFREQ,IREC,ISRC)
                  W = EST(3,IFREQ,IREC,ISRC)
               ELSE !could use in a 2D analysis
                  U = EST(1,IFREQ ,IREC,ISRC)
                  V = CMPLX(0.,0.)
                  W = EST(2,IFREQ ,IREC,ISRC)
               ENDIF
!              NEST = U*COSAZ - V*SINAZ
!              EEST = U*SINAZ + V*COSAZ
               CALL CROTATE(SNGL(AZMOD),U,V, NEST,EEST) 
               ZEST = W !finite elements already in correct frame
               IF (NDIM > 2) THEN
                  EST(1,IFREQ,IREC,ISRC) = NEST
                  EST(2,IFREQ,IREC,ISRC) = EEST
                  EST(3,IFREQ,IREC,ISRC) = ZEST
               ELSE !could use in a 2d anlaysis
                  EST(1,IFREQ,IREC,ISRC) = NEST
                  EST(2,IFREQ,IREC,ISRC) = ZEST
               ENDIF
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE LOBSLST(MDIM, NREC,NDIM, OBS, NOBS,LOBS) 
!
!     Tag if an observation is active for the Jacobian calculation
!
!     INPUT      MEANING
!     -----      ------- 
!     MDIM       leading dimension 
!     NDIM       number of components 
!     NREC       number of receivers
!     
!
!     OUTPUT     MEANING
!     ------     ------- 
!     LOBS       True -> observation exists at [component,receiver]
!     NOBS       number of observations
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,*)
      INTEGER*4 MDIM, NDIM, NREC 
      INTEGER*4, INTENT(OUT) :: NOBS
      LOGICAL*4, INTENT(OUT) :: LOBS(MDIM,*) 
!.... local variables
      INTEGER*4 IREC, I
!
!----------------------------------------------------------------------!
!
!.... loop on obsrevations list 
      NOBS = 0
      LOBS(1:MDIM,1:NREC) = .FALSE.
      DO 1 IREC=1,NREC
         DO 2 I=1,NDIM 
            IF (CABS(OBS(I,IREC)).GT.0.0) THEN
               NOBS = NOBS + 1
               LOBS(I,IREC) = .TRUE. 
            ENDIF
    2    CONTINUE
    1 CONTINUE
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      INTEGER*4 FUNCTION NOBSLST(MDIM, NREC,NDIM, OBS) 
!
!     Tag if an observation is active for the Jacobian calculation
!
!     INPUT      MEANING
!     -----      ------- 
!     MDIM       leading dimension 
!     NDIM       number of components 
!     NREC       number of receivers
!     
!
!     OUTPUT     MEANING
!     ------     ------- 
!     NOBSLST    number of observations
!
!.... variable declarations
      IMPLICIT NONE 
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,*)
      INTEGER*4 MDIM, NDIM, NREC 
!.... local variables
      INTEGER*4 NOBS, IREC, I
!
!----------------------------------------------------------------------!
!
!.... loop on obsrevations list 
      NOBS = 0
      DO 1 IREC=1,NREC
         DO 2 I=1,NDIM 
            IF (CABS(OBS(I,IREC)).GT.0.0) THEN 
               NOBS = NOBS + 1
            ENDIF
    2    CONTINUE
    1 CONTINUE
      NOBSLST = NOBS
      RETURN
      END 
!                                                                      1
!======================================================================!
!                                                                      !
      SUBROUTINE BPRESID_COVD4(MDIM,LDC, NDIM,NREC,IRESTP, IPERM_REC, 
     ;                         DTMAX,FREQ, COVD8, SWGHT, OBS,EST, RESID)
!
!     Data backpropagation with a model covariance matrix
!
!     INPUT      MEANING
!     -----      -------
!     COVD8      data covariance matrix
!     DTMAX      max residual in backpropagation
!     EST        estimates
!     FREQ       current frequency (Hz)
!     IRESTP     residual type (1) phase, (2) amp, (3) both
!     LDC        leading dimension
!     MDIM       leading dimension
!     NDIM       number of components
!     NREC       number of receivers
!     OBS        observations
!     SWGHT      receiver weights
! 
!     OUTPUT     MEANING
!     ------     -------
!     RESID      residuals 
!
!.... variable declarations
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,*), EST(MDIM,*)
      REAL*8, INTENT(IN) :: FREQ
      REAL*4, INTENT(IN) :: SWGHT(MDIM,*), DTMAX
      INTEGER*4, INTENT(IN) :: IPERM_REC(NDIM*NREC), MDIM, LDC, NDIM, 
     ;                         NREC, IRESTP
      COMPLEX*8, INTENT(OUT) :: RESID(MDIM,*)
      REAL*8, INTENT(IN) :: COVD8(LDC,*)
!.... local variables
      REAL*4, ALLOCATABLE :: COVD(:,:), WORK(:), RESIDV(:)
      COMPLEX*8 DIF, CRESID
      REAL*4 DPHI, SPHASE, TWOPI, PI
      INTEGER*4 M, IREC, I, ILOC, JLOC, IROW, JCOL, LROW, KCOL
      PARAMETER(TWOPI = 6.283185307179586)
      PARAMETER(PI = 3.1415926535897931)
!
!----------------------------------------------------------------------!
!
!.... set space and calculate residuals
      M = 2*NREC*NDIM
      ALLOCATE(COVD(M,M))
      ALLOCATE(RESIDV(M))
      ALLOCATE(WORK(M))
      WORK(1:M) = 0.0
      DO 1 IREC=1,NREC
         DO 2 I=1,NDIM
            IF (CABS(OBS(I,IREC)).GT.0.0) THEN
               DIF = CRESID(IRESTP,OBS(I,IREC),EST(I,IREC))
               ILOC =             IPERM_REC( (IREC - 1)*NDIM + I )
               JLOC = NREC*NDIM + IPERM_REC( (IREC - 1)*NDIM + I )
               WORK(ILOC) = REAL(DIF)
               WORK(JLOC) = IMAG(DIF)
               IF (DTMAX.GT.0.0) THEN
                  DPHI = SPHASE(OBS(I,IREC)) - SPHASE(EST(I,IREC))
                  IF (DPHI <-PI) DPHI = DPHI + TWOPI
                  IF (DPHI >+PI) DPHI = DPHI - TWOPI
                  DPHI = DPHI/(TWOPI*SNGL(FREQ)) !convert to seconds
                  IF (ABS(DPHI).GT.DTMAX) THEN
                     WORK(ILOC) = 0.0
                     WORK(JLOC) = 0.0
                  ENDIF
               ENDIF
            ENDIF
    2    CONTINUE
    1 CONTINUE
!
!.... scale columns; post multiplication by diagonal matrix
      I = 1
      IREC = 0
      DO 3 IROW=1,M
         IREC = IREC + 1
         IF (IREC.GT.NREC) THEN
            IREC = 1
            I = I + 1
         ENDIF
         IF (IROW == NDIM*NREC + 1) I = 1
         IF (IROW > NDIM*NREC) THEN
            LROW = NDIM*NREC + IPERM_REC( (IREC - 1)*NDIM + I)
         ELSE
            LROW =             IPERM_REC( (IREC - 1)*NDIM + I)
         ENDIF
         J = 1
         JREC = 0
         DO 4 JCOL=1,M
            JREC = JREC + 1
            IF (JREC.GT.NREC) THEN
               JREC = 1
               J = J + 1
            ENDIF
            IF (JCOL == NDIM*NREC + 1) J = 1
            IF (JCOL > NDIM*NREC) THEN
               KCOL = NDIM*NREC + IPERM_REC( (JREC - 1)*NDIM + J)
            ELSE
               KCOL =             IPERM_REC( (JREC - 1)*NDIM + J) 
            ENDIF
            COVD(LROW,KCOL) = SWGHT(I,IREC)*SNGL(COVD8(LROW,KCOL)) 
     ;                       *SWGHT(J,JREC)
    4    CONTINUE
    3 CONTINUE
!
!.... calculate Sigma^{-1} C_D^{-1} Sigma^{-1} (d - Au)
      CALL SGEMV('N',M,M, 1.0, COVD,M, WORK,1, 0.0,RESIDV,1)
!
!.... put result into RESID
      DO 7 IREC=1,NREC
         DO 8 I=1,NDIM
            ILOC =             IPERM_REC( (IREC - 1)*NDIM + I )
            JLOC = NREC*NDIM + IPERM_REC( (IREC - 1)*NDIM + I )
            RESID(I,IREC) = CMPLX(RESIDV(ILOC),RESIDV(JLOC))
    8    CONTINUE
    7 CONTINUE
      DEALLOCATE(COVD)
      DEALLOCATE(WORK)
      DEALLOCATE(RESIDV)
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BPRESID4(MDIM,NREC,NDIM, NORM,IRESTP, EPSS,
     ;                    FREQ,DTMAX, SWGHT, OBS,EST, RESID)
! 
!     Calculates residuals to backpropagate in choice of norm for 
!     sources and receivers.  If the data is weighted we do that in 
!     the backpropagation routine.  Source: 
!        Which data residual norm for robust elastic frequency-domain 
!        full waveform inversion - Brossier, Operto, and Virieux  
! 
!     INPUT      MEANING
!     -----      ------- 
!     EST        estimate responses for current source/frequency
!     EPSS       scales thresholding criteria (mean|d_{obs}_i|)
!     IRESTP     residual type
!                1->phase
!                2->amplitude 
!                3->phase and amplitude (default) 
!     MDIM       leading dimension for obs and est
!     NDIM       number of components in observations/estimates etc.
!     NORM       =1 ->1 norm
!                =2 ->2 norm (default)
!                =3 ->Huber norm
!                =4 ->Combination
!     NREC       number of receivers
!     OBS        observations for current source/frequency
!     SWGHT      data weights for current source/frequency 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     RESID      residuals to backpropagate for this source 
!
!.... variable declarations
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,*), EST(MDIM,*)
      REAL*8, INTENT(IN) :: FREQ
      REAL*4, INTENT(IN) :: SWGHT(MDIM,*), EPSS, DTMAX
      INTEGER*4, INTENT(IN) :: MDIM,NREC,NDIM, NORM,IRESTP
      COMPLEX*8, INTENT(OUT) :: RESID(MDIM,*) 
!.... local variables
      COMPLEX*8 DIF 
      REAL*4 EPS,EPS2,DENON,WEIGHT2,SPHASE,TWOPI,DPHI
      COMPLEX*8 CRESID,WEIGHT  
      PARAMETER(TWOPI = 6.283185307179586)
      PARAMETER(PI = 3.1415926535897931)
! 
!----------------------------------------------------------------------!
! 
!.... 
      EPS = 0.0
      if (norm > 2) print *, 'bpresid4 not done yet'
!.... calculate residuals 
      EPS2 = EPSS*EPS**2
      DO 1 IREC=1,NREC
         DO 2 I=1,NDIM
            RESID(I,IREC) = CMPLX(0.0,0.0)
            IF (CABS(OBS(I,IREC)).EQ.0.0) GOTO 200 !no observations
            DIF = CRESID(IRESTP,OBS(I,IREC),EST(I,IREC))
            WEIGHT = CMPLX(SWGHT(I,IREC),0.0)
            IF (NORM.EQ.1) THEN !L1 norm
               RESID(I,IREC) = CONJG(WEIGHT)*DIF
     ;                        /CMPLX(CABS(DIF),0.0)
            ELSEIF (NORM.EQ.3) THEN !Huber norm
               IF (CABS(WEIGHT*DIF).LE.EPS) THEN
                  WEIGHT2 = CABS(WEIGHT)**2/EPS
                  RESID(I,IREC) = CMPLX(WEIGHT2,0.0)*DIF
               ELSE
                  RESID(I,IREC) = CONJG(WEIGHT)*DIF
     ;                           /CMPLX(CABS(DIF),0.0)
               ENDIF
            ELSEIF (NORM.EQ.4) THEN !Combination norm 
               WEIGHT2 = CABS(WEIGHT)**2
               DENON = EPS2*SQRT(1.0 + CABS(WEIGHT*DIF)**2/EPS2)
               RESID(I,IREC) = CMPLX(WEIGHT2/DENON,0.0)*DIF
            ELSE !default to L2 
               WEIGHT2 = CABS(WEIGHT)**2
               RESID(I,IREC) = CMPLX(WEIGHT2,0.0)*DIF
               IF (DTMAX > 0.0) THEN
                  DPHI = SPHASE(OBS(I,IREC)) - SPHASE(EST(I,IREC))
                  IF (DPHI <-PI) DPHI = DPHI + TWOPI
                  IF (DPHI >+PI) DPHI = DPHI - TWOPI
                  IF (ABS(DPHI/(TWOPI*FREQ)) > DTMAX) THEN
                     RESID(I,IREC) = CMPLX(0.0,0.0)
                  ENDIF
               ENDIF
            ENDIF
  200       CONTINUE !no observations
    2    CONTINUE
    1 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE FLLRHSGN(MDIM,MFREQ,MREC, NOBS, 
     ;                    NDIM,NFREQ,NREC,NSRC, 
     ;                    NORM,IRESTP, EPSS, SWGHT,
     ;                    OBS,EST, RHS)
! 
!     Fills the RHS for the Gauss-Newton search direction solve
! 
!     INPUT      MEANING
!     -----      ------- 
!     EST        estimate responses (in appropriate N,E,Z frame) 
!     EPSS       scales thresholding criteria (mean|d_{obs}_i|)
!     IRESTP     residual type
!                1->phase
!                2->amplitude 
!                3->phase and amplitude (default) 
!     MDIM       leading dimension for obs and est
!     NDIM       number of components in observations/estimates etc.
!     NOBS       number of obsrevations
!     NORM       =1 ->1 norm
!                =2 ->2 norm (default)
!                =3 ->Huber norm
!                =4 ->Combination
!     NREC       number of receivers
!     OBS        observations 
!     SWGHT      data weights for current source/frequency 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     RHS        RHS for Gauss-Newton or migration
!
!.... variable declarations
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,MFREQ,MREC,*), 
     ;                         EST(MDIM,MFREQ,MREC,*)
      REAL*4, INTENT(IN) :: SWGHT(MDIM,MFREQ,MREC,*), EPSS
      INTEGER*4, INTENT(IN) :: MDIM,NREC,NDIM, NORM,IRESTP
      REAL*4, INTENT(OUT) :: RHS(NOBS)
!.... local variables
      COMPLEX*8 DIF, RESID
      REAL*4 EPS,EPS2,DENON,WEIGHT2
      COMPLEX*8 CRESID,WEIGHT
! 
!----------------------------------------------------------------------!
! 
!.... 
      EPS = 0.0
      if (norm > 2) print *, 'fllrhsgn: not done yet'
!.... calculate if residual exists
      EPS2 = EPSS*EPS**2
      JOBS = 0 
      ILOC = 0 
      DO 1 IFREQ=1,NFREQ
         DO 2 ISRC=1,NSRC 
            JOBS = 0 
            DO 3 IREC=1,NREC
               DO 4 I=1,NDIM
                  DIF = CRESID(IRESTP,OBS(I,IFREQ,IREC,ISRC),
     ;                                        EST(I,IFREQ,IREC,ISRC))
                  IF (CABS(DIF).GT.0.0) JOBS = JOBS + 1 
    4          CONTINUE  
    3       CONTINUE 
            IF (JOBS == 0) GOTO 200 !nothing to do here
            DO 5 IREC=1,NREC
               DO 6 I=1,NDIM 
                  IF (CABS(OBS(I,IFREQ,IREC,ISRC)).EQ.0.0) GOTO 600
                  DIF = CRESID(IRESTP,OBS(I,IFREQ,IREC,ISRC),
     ;                                        EST(I,IFREQ,IREC,ISRC)) 
                  ILOC = ILOC + 1
                  JLOC = ILOC + JOBS
                  WEIGHT = CMPLX(SWGHT(I,IFREQ,IREC,ISRC),0.0)
                  IF (NORM.EQ.1) THEN !L1 norm
                     RESID = CONJG(WEIGHT)*DIF
     ;                      /CMPLX(CABS(DIF),0.0)
                  ELSEIF (NORM.EQ.3) THEN !Huber norm
                     IF (CABS(WEIGHT*DIF).LE.EPS) THEN
                        WEIGHT2 = CABS(WEIGHT)**2/EPS
                        RESID = CMPLX(WEIGHT2,0.0)*DIF
                     ELSE
                        RESID = CONJG(WEIGHT)*DIF
     ;                         /CMPLX(CABS(DIF),0.0)
                     ENDIF
                  ELSEIF (NORM.EQ.4) THEN !Combination norm 
                     WEIGHT2 = CABS(WEIGHT)**2
                     DENON =EPS2*SQRT(1.0 + CABS(WEIGHT*DIF)**2/EPS2)
                     RESID = CMPLX(WEIGHT2/DENON,0.0)*DIF
                  ELSE !default to L2 
                    WEIGHT2 = CABS(WEIGHT)**2
                    RESID = CMPLX(WEIGHT2,0.0)*DIF
                  ENDIF
                  RHS(ILOC) = REAL(RESID)
                  RHS(JLOC) = IMAG(RESID) 
  600             CONTINUE !nothing to do 
    6          CONTINUE !Loop on components
    5       CONTINUE !loop on receivers
            ILOC = JLOC !advance ahead for next block Re, Im block
  200       CONTINUE !dead source
    2    CONTINUE !loop on sources
    1 CONTINUE !loop on frequencies
      print *, iloc-jobs,jloc,nobs
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      COMPLEX*8 FUNCTION CRESID(IRESTP,OBS,EST)
! 
!     Calculates the difference between an observed and estimated value 
! 
!     INPUT      MEANING
!     -----      -------
!     EST        estimated value  
!     IRESTP     residual type: 
!                = 1 -> Phase only 
!                = 2 -> Amplitude only 
!                Otherwise phase and amplitude
!     OBS        observed value 
! 
!     RESULT     MEANING
!     ------     ------- 
!     CRESID     difference between observed and estimated 
! 
!.... variable declarations 
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: OBS,EST
      INTEGER*4, INTENT(IN) :: IRESTP
!.... local variables
      REAL*4 ROBS,ZOBS,REST,ZEST, ADATA,PDATA,AEST,PEST!, T1,T2,PHI,DELTA
      !REAL*4 WRAPPH
! 
!----------------------------------------------------------------------!
! 
!.... easy cases
      IF (OBS.EQ.CMPLX(0.0,0.0)) THEN !no observations set to zero
         CRESID = CMPLX(0.0,0.0)
         RETURN
      ENDIF
! 
!.... calculate residual
      IF (IRESTP.EQ.1 .OR. IRESTP.EQ.2) THEN
         ROBS = REAL(OBS)
         ZOBS = IMAG(OBS)
         REST = REAL(EST)
         ZEST = IMAG(EST)
         IF (IRESTP.EQ.2) THEN !amplitude only
            ADATA = SQRT(ROBS**2 + ZOBS**2)
            PDATA = 0.0
            AEST  = SQRT(REST**2 + ZEST**2)
            PEST  = 0.0
         ELSE !default should probably be phase 
            ADATA = 1.0
            PDATA = ATAN2(ZOBS,ROBS)
            AEST  = 1.0
            PEST  = ATAN2(ZEST,REST)
            !A e^phi_d - A e^phi_e = 2Ai sin(phi_d - phi_e)
         ENDIF
         CRESID = CMPLX(ADATA,0.0)*CEXP(CMPLX(0.0,PDATA))
     ;          - CMPLX(AEST ,0.0)*CEXP(CMPLX(0.0,PEST ))
       ELSE !default is everything, for now
         CRESID = OBS - EST
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      INTEGER*4 FUNCTION ICOBS(MDIM,MFREQ,MREC,NFREQ,NREC,NSRC,NDIM,OBS)
! 
!     Counts number of observations for this frequency
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,MFREQ,MREC,*)
      INTEGER*4, INTENT(IN) :: MDIM,MFREQ,MREC, NFREQ,NREC,NSRC,NDIM
      COMPLEX*8 CZERO
      PARAMETER(CZERO = CMPLX(0.0,0.0))
! 
!----------------------------------------------------------------------!
!
      ICOBS = 0
      DO 1 IFREQ=1,NFREQ
         DO 2 IREC=1,NREC
            DO 3 ISRC=1,NSRC
               DO 4 I=1,NDIM
                  IF (OBS(I,IFREQ,IREC,ISRC).NE.CZERO) ICOBS = ICOBS + 1
    4          CONTINUE
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
      ICOBS = 2*ICOBS !remember each observation has a real and complex
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RESID_TABLE(MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC, 
     ;                       PROJNM, IBLOCK,K, AZMOD, CSIDE, 
     ;                       FREQ,XREC, EST,OBS) 
!
!     Generates residual tables.  The tables are designed if the wave 
!     approaches from the  left, the closest frequency is zero offset
!     and if the wave is from the right, the right receiver is zero
!     offset.  This is because I expect my residuals to increase
!     away from the  
!
!     INPUT      MEANING
!     -----      -------
!     AZMOD      model azimuth (degrees)
!     CSIDE      'R' wave approaching from right, 'L' wave from left
!     EST        estimates (u,v,w) frame
!     IBLOCK     inversion block number
!     K          loop iteration 
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MREC       leading dimension
!     NDIM       number of components 
!     NFREQ      number of frequencies
!     NREC       number of recievers
!     NSRC       number of sources
!     OBS        observations (N,E,Z) frame  
!     PROJNM     project name
!     XREC       x receiver locations
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(80), INTENT(IN) :: PROJNM 
      CHARACTER(1), INTENT(IN) :: CSIDE(NSRC)  
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,MFREQ,MREC,*),
     ;                         EST(MDIM,MFREQ,MREC,*) 
      REAL*8, INTENT(IN) :: FREQ(NFREQ), XREC(NREC), AZMOD 
      INTEGER*4, INTENT(IN) :: MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC,
     ;                         IBLOCK, K   
!.... local variables
      REAL*8, ALLOCATABLE :: XWORK(:)
      INTEGER*4, ALLOCATABLE :: IPERM(:) 
      CHARACTER(80) FILENM
      CHARACTER(5) CBLOCK, CK
      COMPLEX*8 ESTR(3), OBST(3), DIF, QU, QV, QW
      REAL*8 XRMIN, XRMAX, PI180
      REAL*4 OPHASE(3), EPHASE(3), PRESID(3), 
     ;       OAMP(3), EAMP(3), ARESID(3), PI180I, SPHASE
      LOGICAL*4 LISDIR, LEX
      INTEGER*4 ISRC, IFREQ, IREC, I, JREC, JRB, JRE, JDIR,
     ;          IUNITP1, IUNITA1, IUNITP2, IUNITA2, IUNITP3, IUNITA3
      PARAMETER(IUNITP1 = 54, IUNITP2 = 55, IUNITP3 = 56, 
     ;          IUNITA1 = 64, IUNITA2 = 65, IUNITA3 = 66)
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180 
      PARAMETER(PI180I = 57.2957802) !180/pi
!
!----------------------------------------------------------------------!
!
!.... set directory and stuff for file names 
      LEX = LISDIR('./resid_table')
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./resid_table')
      CBLOCK(1:5) = ' '
      CK(1:5) = ' '
      WRITE(CBLOCK,'(I5)') IBLOCK
      WRITE(CK    ,'(I5)') K
      CBLOCK = ADJUSTL(CBLOCK)
      CK = ADJUSTL(CK)
!....... set file name
      FILENM(1:80) = ' '
      FILENM = './resid_table/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'N.pdat'
      FILENM = ADJUSTL(FILENM)
      OPEN(FILE=TRIM(FILENM),UNIT=IUNITP1)
      FILENM(1:80) = ' '
      FILENM = './resid_table/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'E.pdat'
      FILENM = ADJUSTL(FILENM)
      OPEN(FILE=TRIM(FILENM),UNIT=IUNITP2)
      FILENM(1:80) = ' '
      FILENM = './resid_table/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'Z.pdat'
      FILENM = ADJUSTL(FILENM)
      OPEN(FILE=TRIM(FILENM),UNIT=IUNITP3)

      FILENM(1:80) = ' '
      FILENM = './resid_table/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'N.adat'
      FILENM = ADJUSTL(FILENM)
      OPEN(FILE=TRIM(FILENM),UNIT=IUNITA1)
      FILENM(1:80) = ' '
      FILENM = './resid_table/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'E.adat'
      FILENM = ADJUSTL(FILENM)
      OPEN(FILE=TRIM(FILENM),UNIT=IUNITA2)
      FILENM(1:80) = ' '
      FILENM = './resid_table/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'Z.adat'
      FILENM = ADJUSTL(FILENM)
      OPEN(FILE=TRIM(FILENM),UNIT=IUNITA3)


!
!.... sort 
      XRMIN = MINVAL(XREC)
      XRMAX = MAXVAL(XREC) 
      ALLOCATE(XWORK(NREC))
      ALLOCATE(IPERM(NREC))
      XRMIN = XREC(1)
      XRMAX = XREC(1)
      DO 1 IREC=1,NREC
         XWORK(IREC) = XREC(IREC)
         IPERM(IREC) = IREC
         XRMIN = DMIN1(XREC(IREC),XRMIN)
         XRMAX = DMAX1(XREC(IREC),XRMAX)
    1 CONTINUE
      CALL SHELLR8I2(NREC, XWORK,IPERM)
      DEALLOCATE(XWORK)
!
!.... each source has a residual table (frequency,position,residual)
      DO 2 ISRC=1,NSRC

!....... determine loop direction on receivers
         IF (CSIDE(ISRC).EQ.'R') THEN
            JRB = NREC
            JRE = 1
            JDIR =-1
         ELSE
            JRB = 1
            JRE = NREC
            JDIR = 1
         ENDIF
         DO 3 IFREQ=1,NFREQ
            DO 4 JREC=1,NREC !JRB,JRE,JDIR
               !IREC = IPERM(JREC)
               IREC = JREC
               QU = EST(1,IFREQ,IREC,ISRC)
               QV = EST(2,IFREQ,IREC,ISRC)
               QW = EST(3,IFREQ,IREC,ISRC) 
               CALL CROTATE(SNGL(AZMOD),QU,QV, ESTR(1),ESTR(2))
               ESTR(3) = QW
               OBST(1) = OBS(1,IFREQ,IREC,ISRC)
               OBST(2) = OBS(2,IFREQ,IREC,ISRC)
               OBST(3) = OBS(3,IFREQ,IREC,ISRC) 
               EAMP(1) = CABS(ESTR(1))
               EAMP(2) = CABS(ESTR(2))
               EAMP(3) = CABS(ESTR(3))
               EPHASE(1) = SPHASE(ESTR(1))
               EPHASE(2) = SPHASE(ESTR(2))
               EPHASE(3) = SPHASE(ESTR(3)) 
!
!............. loop on components
               DO 5 I=1,NDIM
                  IF (CABS(OBS(I,IFREQ,IREC,ISRC)).GT.0.0) THEN
                     DIF = OBST(I) - ESTR(I) 
                     PRESID(I) = ATAN2(IMAG(DIF),REAL(DIF))
                     ARESID(I) = CABS(DIF)
                  ELSE
                     PRESID(I) = 0.0
                     ARESID(I) = 0.0
                     OPHASE(I) = 0.0
                     OAMP(I) = 0.0 
                  ENDIF
    5          CONTINUE !loop on components
               WRITE(IUNITP1,905) CSIDE(ISRC), 
     ;                            SNGL(FREQ(IFREQ)),XREC(IREC),
     ;                            IREC,ISRC,
     ;                            OPHASE(1),EPHASE(1),PRESID(1)
               WRITE(IUNITP2,905) CSIDE(ISRC), 
     ;                            SNGL(FREQ(IFREQ)),XREC(IREC),
     ;                            IREC,ISRC,
     ;                            OPHASE(2),EPHASE(2),PRESID(2)
               WRITE(IUNITP3,905) CSIDE(ISRC), 
     ;                            SNGL(FREQ(IFREQ)),XREC(IREC),
     ;                            IREC,ISRC,
     ;                            OPHASE(3),EPHASE(3),PRESID(3)
  905          FORMAT(A1,' ',F12.5,F12.2,2I5,3E14.5) 
               WRITE(IUNITA1,906) CSIDE(ISRC),
     ;                            SNGL(FREQ(IFREQ)),XREC(IREC),
     ;                            IREC,ISRC,
     ;                            OAMP(1),EAMP(1),ARESID(1)
               WRITE(IUNITA2,906) CSIDE(ISRC), 
     ;                            SNGL(FREQ(IFREQ)),XREC(IREC), 
     ;                            IREC,ISRC,
     ;                            OAMP(2),EAMP(2),ARESID(2)
               WRITE(IUNITA3,906) CSIDE(ISRC),
     ;                            SNGL(FREQ(IFREQ)),XREC(IREC),
     ;                            IREC,ISRC,
     ;                            OAMP(3),EAMP(3),ARESID(3)
  906          FORMAT(A1,' ',F12.5,F12.2,2I5,3E16.5) 
    4       CONTINUE !loop on receivers
            WRITE(IUNITP1,*) !leave a blank line
            WRITE(IUNITP2,*) 
            WRITE(IUNITP3,*) !leave a blank line
            WRITE(IUNITA1,*)
            WRITE(IUNITA2,*) !leave a blank line
            WRITE(IUNITA3,*)
    3    CONTINUE !loop on frequencies 
    2 CONTINUE !loop on sources
      DEALLOCATE(IPERM)
      CLOSE(IUNITP1)
      CLOSE(IUNITP2)
      CLOSE(IUNITP3)
      CLOSE(IUNITA1)
      CLOSE(IUNITA2)
      CLOSE(IUNITA3) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RTABLE_RMS(MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC,
     ;                       PROJNM, IBLOCK,K,NORM, AZMOD, CSIDE, 
     ;                       FREQ,XREC,YREC, SWGHT,EST,OBS)
!
!     Plots the RMS misfits for phase and amplitude at frequency
!     The tables are designed to be in order from model left to right 
!
!     INPUT      MEANING
!     -----      -------
!     AZMOD      model azimuth (degrees)
!     EST        estimates (u,v,w) frame
!     IBLOCK     inversion block number
!     K          loop iteration 
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MREC       leading dimension
!     NDIM       number of components 
!     NFREQ      number of frequencies
!     NORM       (1) L1 norm, otherwise (2) norm
!     NREC       number of recievers
!     NSRC       number of sources
!     OBS        observations (N,E,Z) frame  
!     PROJNM     project name
!     XREC       x receiver locations
!     YREC       y receiver locations
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(80), INTENT(IN) :: PROJNM 
      CHARACTER(1), INTENT(IN) :: CSIDE(NSRC)
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,MFREQ,MREC,*),
     ;                         EST(MDIM,MFREQ,MREC,*) 
      REAL*8, INTENT(IN) :: FREQ(NFREQ), XREC(NREC), YREC(NREC), AZMOD 
      REAL*4, INTENT(IN) :: SWGHT(MDIM,MFREQ,MREC,*) 
      INTEGER*4, INTENT(IN) :: MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC,
     ;                         NORM, IBLOCK, K 
!.... local variables
      REAL*8, ALLOCATABLE :: XWORK(:)
      REAL*4, ALLOCATABLE :: RMSPTAB(:,:), RMSATAB(:,:),  RMSTAB(:,:)
      INTEGER*4, ALLOCATABLE :: NOBSTAB(:,:), IPERM(:)
      CHARACTER(80) FILENM
      CHARACTER(5) CBLOCK, CK
      COMPLEX*8 ESTR(3), DIFA, DIFP, DIF, CRESID, COSAZ, SINAZ, S
      REAL*8 XRMIN, XRMAX, PI180
      REAL*4 XOFF, PI180I, RESIDP, RESIDA, RESID
      LOGICAL*4 LISDIR, LEX
      INTEGER*4 NOBS, IOBS, ISRC, IFREQ, IREC, I, JREC, 
     ;          KREC, IUNITPL,IUNITAL,IUNITL, IUNITPR,IUNITAR,IUNITR,
     ;          IUNITP,IUNITA,IUNITB
      PARAMETER(IUNITPL = 54, IUNITAL = 55, IUNITL = 56,  
     ;          IUNITPR = 57, IUNITAR = 58, IUNITR = 59, 
     ;          IUNITP  = 60, IUNITA  = 61, IUNITB = 62)
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180 
      PARAMETER(PI180I = 57.2957802) !180/pi
!
!----------------------------------------------------------------------!
!
!.... set directory and stuff for file names 
      LEX = LISDIR('./resid_table')
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./resid_table')
      CBLOCK(1:5) = ' '
      CK(1:5) = ' '
      WRITE(CBLOCK,'(I5)') IBLOCK
      WRITE(CK    ,'(I5)') K
      CBLOCK = ADJUSTL(CBLOCK)
      CK = ADJUSTL(CK)
!.... waves from left side
      FILENM(1:80) = ' ' 
      FILENM = './resid_table/'//TRIM(ADJUSTL(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'_RMSphaseL.dat'
      FILENM = ADJUSTL(FILENM)
      OPEN(UNIT=IUNITPL,FILE=TRIM(FILENM))
      FILENM(1:80) = ' ' 
      FILENM = './resid_table/'//TRIM(ADJUSTL(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'_RMSampL.dat'
      FILENM = ADJUSTL(FILENM)
      OPEN(UNIT=IUNITAL,FILE=TRIM(FILENM))
      FILENM(1:80) = ' ' 
      FILENM = './resid_table/'//TRIM(ADJUSTL(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'_RMSbothL.dat'
      FILENM = ADJUSTL(FILENM)
      OPEN(UNIT=IUNITL,FILE=TRIM(FILENM))
!.... waves from right
      FILENM(1:80) = ' ' 
      FILENM = './resid_table/'//TRIM(ADJUSTL(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'_RMSphaseR.dat'
      FILENM = ADJUSTL(FILENM)
      OPEN(UNIT=IUNITPR,FILE=TRIM(FILENM))
      FILENM(1:80) = ' '
      FILENM = './resid_table/'//TRIM(ADJUSTL(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'_RMSampR.dat'
      FILENM = ADJUSTL(FILENM)
      OPEN(UNIT=IUNITAR,FILE=TRIM(FILENM))
      FILENM(1:80) = ' '
      FILENM = './resid_table/'//TRIM(ADJUSTL(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'_RMSbothR.dat'
      FILENM = ADJUSTL(FILENM)
      OPEN(UNIT=IUNITR,FILE=TRIM(FILENM))
!.... both 
      FILENM(1:80) = ' ' 
      FILENM = './resid_table/'//TRIM(ADJUSTL(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'_RMSphase.dat'
      FILENM = ADJUSTL(FILENM)
      OPEN(UNIT=IUNITP,FILE=TRIM(FILENM))
      FILENM(1:80) = ' ' 
      FILENM = './resid_table/'//TRIM(ADJUSTL(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'_RMSamp.dat'
      FILENM = ADJUSTL(FILENM)
      OPEN(UNIT=IUNITA,FILE=TRIM(FILENM))
      FILENM(1:80) = ' ' 
      FILENM = './resid_table/'//TRIM(ADJUSTL(PROJNM))//'-'//
     ;         TRIM(CBLOCK)//'-'//TRIM(CK)//'_RMSboth.dat'
      FILENM = ADJUSTL(FILENM)
      OPEN(UNIT=IUNITB,FILE=TRIM(FILENM))
!
!.... rotation angle
      COSAZ = CMPLX(SNGL(DCOS(AZMOD*PI180)),0.0)
      SINAZ = CMPLX(SNGL(DSIN(AZMOD*PI180)),0.0)
!
!.... sort on distance from left most receiver 
      XRMIN = MINVAL(XREC) 
      XRMAX = MAXVAL(XREC) 
      ALLOCATE(XWORK(NREC))
      ALLOCATE(IPERM(NREC))
      DO 1 IREC=1,NREC
         XWORK(IREC) = DSQRT((XREC(IREC) - XRMIN)**2 + YREC(IREC)**2)
         IPERM(IREC) = IREC
    1 CONTINUE
      CALL SHELLR8I2(NREC, XWORK,IPERM)
      DEALLOCATE(XWORK)
!
!.... set space
      ALLOCATE(RMSPTAB(3*NREC,NDIM)) !phase RMS difference table
      ALLOCATE(RMSATAB(3*NREC,NDIM)) !amplitude RMS difference table
      ALLOCATE(RMSTAB(3*NREC,NDIM))  !phase and amp RMS difference table
      ALLOCATE(NOBSTAB(3*NREC,NDIM)) !counter for number of observations 
!
!.... loop on receivers
      DO 2 IFREQ=1,NFREQ
         DO 3 JREC=1,NREC
            IREC = IPERM(JREC) 
            DO 4 I=1,NDIM
               KREC = JREC !left most receiver is "close"
               RMSPTAB(KREC,I) = 0.0
               RMSATAB(KREC,I) = 0.0
               RMSTAB (KREC,I) = 0.0
               NOBSTAB(KREC,I) = 0
               KREC = 2*NREC + 1 - JREC !left most receiver is "far" 
               RMSPTAB(KREC,I) = 0.0 
               RMSATAB(KREC,I) = 0.0 
               RMSTAB (KREC,I) = 0.0 
               NOBSTAB(KREC,I) = 0 
               KREC = 2*NREC + JREC !just keep them all
               RMSPTAB(KREC,I) = 0.0 
               RMSATAB(KREC,I) = 0.0 
               RMSTAB (KREC,I) = 0.0 
               NOBSTAB(KREC,I) = 0 
    4       CONTINUE
            DO 5 ISRC=1,NSRC
               ESTR(1) = EST(1,IFREQ,IREC,ISRC)*COSAZ
     ;                 - EST(2,IFREQ,IREC,ISRC)*SINAZ
               ESTR(2) = EST(1,IFREQ,IREC,ISRC)*SINAZ
     ;                 + EST(2,IFREQ,IREC,ISRC)*COSAZ
               ESTR(3) = EST(3,IFREQ,IREC,ISRC)
               DO 6 I=1,NDIM
                  DIFP = CRESID(1,OBS(I,IFREQ,IREC,ISRC),
     ;                            ESTR(I))
                  DIFA = CRESID(2,OBS(I,IFREQ,IREC,ISRC),
     ;                            ESTR(I)) 
                  DIF  = CRESID(3,OBS(I,IFREQ,IREC,ISRC),
     ;                            ESTR(I)) 
                  DIFP = CMPLX(ATAN2(IMAG(DIFP),REAL(DIFP))*PI180I,0.0)
                  S = CMPLX(SWGHT(I,IFREQ,IREC,ISRC),0.0)
                  IF (CABS(DIF).GT.0.0) THEN
                     IF (NORM.EQ.1) THEN
                        RESIDP = CABS(S*DIFP)
                        RESIDA = CABS(S*DIFA)
                        RESID  = CABS(S*DIF) 
                     ELSE
                        RESIDP = CABS(S)**2*CABS(DIFP)**2
                        RESIDA = CABS(S)**2*CABS(DIFA)**2
                        RESID  = CABS(S)**2*CABS(DIF)**2
                     ENDIF
                     IOBS = 1
                  ELSE 
                     RESIDP = 0.0
                     RESIDA = 0.0
                     RESID = 0.0
                     IOBS = 0 
                  ENDIF
                  IF (CSIDE(ISRC).EQ.'L') THEN !left is close
                     KREC = JREC
                  ELSE !left is far
                     KREC = 2*NREC + 1 - JREC
                  ENDIF
                  RMSPTAB(KREC,I) = RMSPTAB(KREC,I) + RESIDP**2
                  RMSATAB(KREC,I) = RMSATAB(KREC,I) + RESIDA**2
                  RMSTAB (KREC,I) = RMSTAB (KREC,I) + RESID**2 
                  NOBSTAB(KREC,I) = NOBSTAB(KREC,I) + IOBS 
                  KREC = 2*NREC + JREC 
                  RMSPTAB(KREC,I) = RMSPTAB(KREC,I) + RESIDP**2
                  RMSATAB(KREC,I) = RMSATAB(KREC,I) + RESIDA**2
                  RMSTAB (KREC,I) = RMSTAB (KREC,I) + RESID**2
                  NOBSTAB(KREC,I) = NOBSTAB(KREC,I) + IOBS 
    6          CONTINUE !Loop on componenets
    5       CONTINUE !loop on sources
    3    CONTINUE !loop on receivers
!
!....... now calculate and dump RMS 
         DO 7 JREC=1,NREC
            DO 8 I=1,NDIM
               KREC = JREC
               NOBS = NOBSTAB(KREC,I)
               IF (NOBS.GT.0) THEN  
                  RMSPTAB(KREC,I) = SQRT(RMSPTAB(KREC,I)/FLOAT(NOBS))
                  RMSATAB(KREC,I) = SQRT(RMSATAB(KREC,I)/FLOAT(NOBS))
                  RMSTAB (KREC,I) = SQRT(RMSTAB (KREC,I)/FLOAT(NOBS))
               ENDIF 
               KREC = 2*NREC + 1 - JREC
               NOBS = NOBSTAB(KREC,I)
               IF (NOBS.GT.0) THEN 
                  RMSPTAB(KREC,I) = SQRT(RMSPTAB(KREC,I)/FLOAT(NOBS))
                  RMSATAB(KREC,I) = SQRT(RMSATAB(KREC,I)/FLOAT(NOBS))
                  RMSTAB (KREC,I) = SQRT(RMSTAB (KREC,I)/FLOAT(NOBS))
               ENDIF 
               KREC = 2*NREC + JREC
               NOBS = NOBSTAB(KREC,I)
               IF (NOBS.GT.0) THEN
                  RMSPTAB(KREC,I) = SQRT(RMSPTAB(KREC,I)/FLOAT(NOBS))
                  RMSATAB(KREC,I) = SQRT(RMSATAB(KREC,I)/FLOAT(NOBS))
                  RMSTAB (KREC,I) = SQRT(RMSTAB (KREC,I)/FLOAT(NOBS))
               ENDIF
    8       CONTINUE !loop on components 
            IREC = IPERM(JREC) 
            KREC = JREC
            XOFF = SNGL(DSQRT((XREC(IREC) - XRMIN)**2 + YREC(IREC)**2))
            WRITE(IUNITPL,*) XOFF,FREQ(IFREQ),IREC, RMSPTAB(KREC,1:NDIM)
            WRITE(IUNITAL,*) XOFF,FREQ(IFREQ),IREC, RMSATAB(KREC,1:NDIM)
            WRITE(IUNITL ,*) XOFF,FREQ(IFREQ),IREC, RMSTAB (KREC,1:NDIM)
            KREC = 2*NREC + 1 - JREC
            XOFF = SNGL(DSQRT((XRMAX - XREC(IREC))**2 + YREC(IREC)**2)) 
            WRITE(IUNITPR,*) XOFF,SNGL(FREQ(IFREQ)),IREC, 
     ;                       RMSPTAB(KREC,1:NDIM)
            WRITE(IUNITAR,*) XOFF,SNGL(FREQ(IFREQ)),IREC,
     ;                       RMSATAB(KREC,1:NDIM)
            WRITE(IUNITR ,*) XOFF,SNGL(FREQ(IFREQ)),IREC, 
     ;                       RMSTAB (KREC,1:NDIM)
            KREC = 2*NREC + JREC 
            XOFF = SNGL(DSQRT((XREC(IREC) - XRMIN)**2 + YREC(IREC)**2))
            WRITE(IUNITP,*) XOFF,SNGL(FREQ(IFREQ)),IREC,  
     ;                      RMSPTAB(KREC,1:NDIM)
            WRITE(IUNITA,*) XOFF,SNGL(FREQ(IFREQ)),IREC,
     ;                      RMSATAB(KREC,1:NDIM)
            WRITE(IUNITB,*) XOFF,SNGL(FREQ(IFREQ)),IREC, 
     ;                      RMSTAB(KREC,1:NDIM) 
    7    CONTINUE !loop on receivers 
         WRITE(IUNITPL,*) !leave a blank line 
         WRITE(IUNITAL,*)
         WRITE(IUNITL ,*)
         WRITE(IUNITPR,*)
         WRITE(IUNITAR,*)
         WRITE(IUNITR ,*)
         WRITE(IUNITP,*)
         WRITE(IUNITA,*)
         WRITE(IUNITB,*)
    2 CONTINUE !loop on frequencies
!
!.... close files
      CLOSE(IUNITPL)
      CLOSE(IUNITAL)
      CLOSE(IUNITL)
      CLOSE(IUNITPR)
      CLOSE(IUNITAR)
      CLOSE(IUNITR)  
      CLOSE(IUNITP)
      CLOSE(IUNITA)
      CLOSE(IUNITB)
      DEALLOCATE(RMSPTAB)
      DEALLOCATE(RMSATAB)  
      DEALLOCATE(RMSTAB)
      DEALLOCATE(NOBSTAB)       
      DEALLOCATE(IPERM) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CONVEST_STF(MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC,
     ;                       SOURCE, EST) 
!
!     Convolves estimates with a new source time function 
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: SOURCE(MFREQ,nsrc) 
      INTEGER*4, INTENT(IN) :: MDIM, MFREQ, MREC, NDIM, NFREQ, NREC, 
     ;                         NSRC
      COMPLEX*8, INTENT(INOUT) :: EST(MDIM,MFREQ,MREC,*)
!.... local variables
      INTEGER*4 IFREQ, IREC, ISRC, I 
!
!----------------------------------------------------------------------!
!
!.... multiply each receiver source pari by new STF
      DO 1 IFREQ=1,NFREQ
         DO 2 IREC=1,NREC
            DO 3 ISRC=1,NSRC
               DO 4 I=1,NDIM
                  EST(I,IFREQ,IREC,ISRC) = EST(I,IFREQ,IREC,ISRC)
     ;                                    *SOURCE(IFREQ,ISRC)
    4          CONTINUE !loop on components
    3       CONTINUE !loop on sources
    2    CONTINUE !loop on receivers
    1 CONTINUE !loop on frequencies
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CONVEST_RRF(MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC,
     ;                       RECV, EST) 
!
!     Convolves estimates with new receiver response functions 
      IMPLICIT NONE 
      COMPLEX*8, INTENT(IN) :: RECV(MDIM,MFREQ,nrec) 
      INTEGER*4, INTENT(IN) :: MDIM, MFREQ, MREC, NDIM, NFREQ, NREC, 
     ;                         NSRC
      COMPLEX*8, INTENT(INOUT) :: EST(MDIM,MFREQ,MREC,*)
!.... local variables
      INTEGER*4 IFREQ, IREC, ISRC, I 
!
!----------------------------------------------------------------------!
!
!.... multiply each receiver source pari by new RRF 
      DO 1 IFREQ=1,NFREQ
         DO 2 IREC=1,NREC
            DO 3 ISRC=1,NSRC
               DO 4 I=1,NDIM
                  EST(I,IFREQ,IREC,ISRC) = EST(I,IFREQ,IREC,ISRC)
     ;                                    *RECV(I,IFREQ,IREC)
    4          CONTINUE !loop on components
    3       CONTINUE !loop on sources
    2    CONTINUE !loop on receivers
    1 CONTINUE !loop on frequencies
      RETURN
      END  
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DUMP_RESID(PROJNM,MDIM,MFREQ,MREC, 
     ;                      NDIM,NFREQ,NREC,NSRC, 
     ;                      LSURF,LVERT,IRESTP,K, 
     ;                      AZMOD,LINVF,CFTYPE,SRCTYP,FREQ, 
     ;                      WGHTS,OBS,EST) 
!
!     Writes the residuals to file for viewing
!
!     INPUT      MEANING
!     -----      ------- 
!     AZMOD      model azimuth
!     CFTYPE     character descriptor for frequency type
!     EST        estimates (u,v,w)
!     FREQ       frequency list
!     IRESTP     residual type (1) phase (2) amplitude (3) both
!     K          iteration number
!     LINVF      True -> frequency is an inversion frequency
!     LSURF      True -> surface wave residuals
!     LVERT      True -> vertical component only
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MREC       leading dimension
!     NDIM       number of components
!     NFREQ      number of frequencies
!     NREC       number of receivers
!     NSRC       number of sources
!     OBS        observations (N,E,Z)
!     PROJNM     project name
!     SRCTYP     character descriptor for source type
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: PROJNM
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC)
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,MFREQ,MREC,NSRC),  
     ;                         EST(MDIM,MFREQ,MREC,NSRC)
      REAL*8, INTENT(IN) :: FREQ(NFREQ), AZMOD
      REAL*4, INTENT(IN) :: WGHTS(MDIM,MFREQ,MREC,NSRC)
      INTEGER*4, INTENT(IN) :: MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC, 
     ;                         IRESTP, K
      LOGICAL*4, INTENT(IN) :: LINVF(NFREQ), LSURF, LVERT 
!.... local variables
      CHARACTER(80) OUTFL
      CHARACTER(5) CSRC, CK
      COMPLEX*8 U, V ,W, CN, CE, CZ, COBS, CEST
      REAL*8 TWOPI, PI180
      REAL*4 OMEGA, PHASEO, PHASEE, AMPO, AMPE, DA, DT, SPHASE, THETA,
     ;       WT
      INTEGER*4 ISRC, IFREQ, IREC, I, IUNIT1, IUNIT2, IUNIT3
      LOGICAL*4 LISDIR
      PARAMETER(TWOPI = 6.283185307179586D0)
      PARAMETER(PI180 = 0.017453292519943295D0)
      PARAMETER(IUNIT1 = 66)
      PARAMETER(IUNIT2 = 67)
      PARAMETER(IUNIT3 = 68)
!
!----------------------------------------------------------------------!
!
!.... initialize
      THETA = SNGL(AZMOD*PI180)
      IF (.NOT.LISDIR('./dump_resid')) CALL SYSTEM('mkdir ./dump_resid')
      CK(1:5) = ' '
      WRITE(CK,'(I5)') K
      CK = ADJUSTL(TRIM(CK))
!
!.... each source gets a residual file
      DO 1 ISRC=1,NSRC
         OUTFL(1:80) = ' '
         CSRC(1:5) = ' '
         WRITE(CSRC,'(I5)') ISRC
         CSRC = ADJUSTL(TRIM(CSRC))
         IF (     LSURF.AND.SRCTYP(ISRC)(1:1).NE.'S') GOTO 110
         IF (.NOT.LSURF.AND.SRCTYP(ISRC)(1:1).NE.'B') GOTO 110
         IF (LVERT) THEN
            IF (LSURF) THEN
               OUTFL = './dump_resid/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;                 TRIM(ADJUSTL(CK))//'-'//TRIM(ADJUSTL(CSRC))//
     ;                 '_srf.dat'
            ELSE
               OUTFL = './dump_resid/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;                 TRIM(ADJUSTL(CK))//'-'//TRIM(ADJUSTL(CSRC))//
     ;                 '_bdy.dat'
            ENDIF
            OPEN(UNIT=IUNIT3,FILE=TRIM(ADJUSTL(OUTFL)))
            DO 2 IFREQ=1,NFREQ
               IF (     LSURF.AND.CFTYPE(IFREQ).EQ.'S' .OR. 
     ;             .NOT.LSURF.AND.CFTYPE(IFREQ).EQ.'B') THEN
                  IF (LINVF(IFREQ)) THEN
                     OMEGA = SNGL(TWOPI*FREQ(IFREQ))
                     DO 3 IREC=1,NREC
                        AMPO = CABS(OBS(3,IFREQ,IREC,ISRC))
                        AMPE = CABS(EST(3,IFREQ,IREC,ISRC))
                        PHASEO = SPHASE(OBS(3,IFREQ,IREC,ISRC))
                        PHASEE = SPHASE(EST(3,IFREQ,IREC,ISRC))
                        WT = WGHTS(3,IFREQ,IREC,ISRC)
                        DA = AMPO - AMPE
                        DT = (PHASEO - PHASEE)/OMEGA
                        IF (IRESTP.EQ.1) DA = 0.0
                        IF (IRESTP.EQ.2) DT = 0.0
                        IF (CABS(OBS(3,IFREQ,IREC,ISRC)).EQ.0.0) THEN
                           DA = 0.0
                           DT = 0.0
                           WT = 0.0
                        ENDIF
                        WRITE(IUNIT3,950) ISRC,IREC,
     ;                                    SNGL(FREQ(IFREQ)),WT,DA,DT
    3                CONTINUE
                     WRITE(IUNIT3,*)
                  ENDIF
               ENDIF
    2       CONTINUE
            CLOSE(IUNIT3)
         ELSE
            IF (LSURF) THEN
               OUTFL = './dump_resid/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;                 TRIM(ADJUSTL(CK))//'-'//TRIM(ADJUSTL(CSRC))//
     ;                '-1_srf.dat'
            ELSE
               OUTFL = './dump_resid/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;                 TRIM(ADJUSTL(CK))//'-'//TRIM(ADJUSTL(CSRC))//
     ;                 '-1_bdy.dat'
            ENDIF
            OPEN(UNIT=IUNIT1,FILE=TRIM(ADJUSTL(OUTFL)))
            IF (LSURF) THEN
               OUTFL = './dump_resid/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;                 TRIM(ADJUSTL(CK))//'-'//TRIM(ADJUSTL(CSRC))//
     ;                '-2_srf.dat'
            ELSE
               OUTFL = './dump_resid/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;                 TRIM(ADJUSTL(CK))//'-'//TRIM(ADJUSTL(CSRC))//
     ;                 '-2_bdy.dat'
            ENDIF
            OPEN(UNIT=IUNIT2,FILE=TRIM(ADJUSTL(OUTFL)))
            IF (LSURF) THEN
               OUTFL = './dump_resid/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;                 TRIM(ADJUSTL(CK))//'-'//TRIM(ADJUSTL(CSRC))//
     ;                '-3_srf.dat'
            ELSE
               OUTFL = './dump_resid/'//ADJUSTL(TRIM(PROJNM))//'-'//
     ;                 TRIM(ADJUSTL(CK))//'-'//TRIM(ADJUSTL(CSRC))//
     ;                 '-3_bdy.dat'
            ENDIF
            OPEN(UNIT=IUNIT3,FILE=TRIM(ADJUSTL(OUTFL)))
            DO 4 I=1,NDIM
               DO 5 IFREQ=1,NFREQ
                  IF (     LSURF.AND.CFTYPE(IFREQ).EQ.'S' .OR. 
     ;                .NOT.LSURF.AND.CFTYPE(IFREQ).EQ.'B') THEN 
                     IF (LINVF(IFREQ)) THEN
                        OMEGA = SNGL(TWOPI*FREQ(IFREQ))
                        DO 6 IREC=1,NREC
                           U = CMPLX(EST(1,IFREQ,IREC,ISRC))
                           V = CMPLX(EST(2,IFREQ,IREC,ISRC))
                           W = CMPLX(EST(3,IFREQ,IREC,ISRC))
                           CALL CROTATE(THETA,U,V, CN,CE)
                           CZ = W
                           IF (I.EQ.1) THEN
                              COBS = OBS(I,IFREQ,IREC,ISRC)
                              CEST = CN
                           ELSEIF (I.EQ.2) THEN
                              COBS = OBS(I,IFREQ,IREC,ISRC)
                              CEST = CE
                           ELSE
                              COBS = OBS(I,IFREQ,IREC,ISRC)
                              CEST = CZ
                           ENDIF
                           AMPO = CABS(COBS)
                           AMPE = CABS(CEST)
                           PHASEO = SPHASE(COBS)
                           PHASEE = SPHASE(CEST)
                           WT = WGHTS(I,IFREQ,IREC,ISRC)
                           DA = AMPO - AMPE
                           DT = (PHASEO - PHASEE)/OMEGA
                           IF (IRESTP == 1) DA = 0.0
                           IF (IRESTP == 2) DT = 0.0
                           IF (IRESTP == 1) DA = 0.0
                           IF (IRESTP == 2) DT = 0.0
                           IF (CABS(COBS).EQ.0.0) THEN
                              DA = 0.0
                              DT = 0.0
                              WT = 0.0
                           ENDIF 
                           IF (I.EQ.1)
     ;                     WRITE(IUNIT1,950) ISRC,IREC,
     ;                                       SNGL(FREQ(IFREQ)),WT,DA,DT
                           IF (I.EQ.2)
     ;                     WRITE(IUNIT2,950) ISRC,IREC,
     ;                                       SNGL(FREQ(IFREQ)),WT,DA,DT
                           IF (I.EQ.3)
     ;                     WRITE(IUNIT3,950) ISRC,IREC,
     ;                                       SNGL(FREQ(IFREQ)),WT,DA,DT
    6                   CONTINUE !loop on receivers
                        IF (I.EQ.1) WRITE(IUNIT1,*)
                        IF (I.EQ.2) WRITE(IUNIT2,*)
                        IF (I.EQ.3) WRITE(IUNIT3,*)
                     ENDIF !end check on inversion frequency
                  ENDIF !end check on frequency/type match
    5          CONTINUE !Loop on frequencies
    4       CONTINUE
            CLOSE(IUNIT1)
            CLOSE(IUNIT2)
            CLOSE(IUNIT3)
         ENDIF
  110    CONTINUE !source mismatch
    1 CONTINUE
  950 FORMAT(I4,I4,F12.5,F12.5,E16.4,F12.5)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE NULL_OBSFD(MDIM,MFREQ,MREC, 
     ;                      NDIM,NFREQ,NREC,NSRC, NCUT, LVERB,
     ;                      LCOVD_SRF,LCOVD_BDY, CFTYPE,SRCTYP, 
     ;                      WGHTS,OBS)
!
!     If we are finite differencing the phaes observations we need a minimum 
!     number of observations.  If we do not meet this minimum we kill the frequency
!
!     CFTYPE     frequency type; body or surface wave
!     LCOVD_BDY  True -> using body wave data covar matrix
!     LCOVD_SRF  True -> using surface wave data covar matrix
!     LVERB      contrls verbosity
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MREC       leading dimension
!     NCUT       minimum number of active receivers for a (source,freq)
!     NDIM       number of components
!     NFREQ      number of frequencies
!     NREC       number of receivers
!     NSRC       number of sources
!     SRCTYP     source type 
!     
!     INPUT      MEANING
!     OBS        updated observations
!     WGHTS      updated weights
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC)
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      INTEGER*4, INTENT(IN) :: MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC, 
     ;                         NCUT
      LOGICAL*4, INTENT(IN) :: LVERB, LCOVD_SRF, LCOVD_BDY
      COMPLEX*8, INTENT(INOUT) :: OBS(MDIM,MFREQ,MREC,*)
      REAL*4, INTENT(INOUT) :: WGHTS(MDIM,MFREQ,MREC,*)
!.... local variables
      INTEGER*4 NOBS, IFREQ,ISRC,IREC,I
      COMPLEX*8, PARAMETER :: CZERO = CMPLX(0.0,0.0)
!
!----------------------------------------------------------------------!
!
!.... early return
      IF (.NOT.LCOVD_SRF .AND. .NOT.LCOVD_BDY) RETURN
      IF (LVERB .AND. LCOVD_SRF) THEN 
         WRITE(*,*) 'null_obsfd: Surface wave threshold is',
     ;   NCUT,'observations'
      ENDIF
      IF (LVERB .AND. LCOVD_BDY) THEN
         WRITE(*,*) 'null_obsfd: Body wave threshold is',
     ;   NCUT,'observations'
      ENDIF
!
!.... loop on frequencies then sources
      DO 1 IFREQ=1,NFREQ
         DO 2 ISRC=1,NSRC
!
!.......... check we are interested in potentially muting
            IF (CFTYPE(IFREQ).EQ.'S' .AND. SRCTYP(ISRC)(1:1).EQ.'S' 
     ;          .AND. LCOVD_SRF .OR. 
     ;          CFTYPE(IFREQ).EQ.'B' .AND. SRCTYP(ISRC)(1:1).EQ.'P' 
     ;          .AND. LCOVD_BDY) THEN
               DO 3 I=1,NDIM 
!
!................ make sure the receivers have at least ncut obs
                  NOBS = 0
                  DO 4 IREC=1,NREC
                     IF (CABS(OBS(I,IFREQ,IREC,ISRC)) > 0.0) 
     ;               NOBS = NOBS + 1
    4             CONTINUE
                  IF (NOBS.LT.NCUT .AND. NOBS.GT.0) THEN
                     IF (LVERB) THEN 
                        WRITE(*,*) 'null_obsfd: Nulling frequency source
     ; pair:',IFREQ,ISRC
                     ENDIF
                     OBS(I,IFREQ,1:NREC,ISRC) = CZERO
                     WGHTS(I,IFREQ,1:NREC,ISRC) = 0.0
                  ENDIF
    3          CONTINUE
            ENDIF
    2    CONTINUE
    1 CONTINUE
      RETURN
      END
 
