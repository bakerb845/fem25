      SUBROUTINE RECUPD(MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC, 
     ;                  LUNWRAP,IRESTP, AZMOD, EST,OBS, RECV)
!
!     Updates the receiver statics.  Here we update the receiver statics
!     on all components.  This is a full Newton step and does not 
!     require iteration since the objective function is from an 
!     inner product,e L2 norm, and amounts to finding the minima of 
!     a parabaloid. - B. Baker February 2013 
!
!     INPUT      MEANING
!     -----      ------- 
!     AZMOD      model azimuth (degrees)
!     EST        estimates at receivers (already rotated!)
!     IRESTP     =1 -> Phase only 
!                =2 -> Amplitude only 
!                =3 -> phase and amplitude (default)
!     LUNWRAP    True -> data is unwrapped and stored (amp,phase)
!                False -> data is stored as a complex number 
!     MDIM       leading dimension for estimates/observations
!     MFREQ      leading dimension
!     NDIM       number of components for receivers
!     NFREQ      number of frequencies 
!     NREC       number of receivers
!     NSRC       number of sources
!     OBS        observations at receivers 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     RECV       updated receiver statics (u,v,w)
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,MFREQ,MREC,*), 
     ;                         EST(MDIM,MFREQ,MREC,*)  
      REAL*8, INTENT(IN) :: AZMOD 
      INTEGER*4, INTENT(IN) :: MFREQ, MDIM, MREC, NFREQ, NSRC, NREC, 
     ;                         NDIM, IRESTP
      LOGICAL*4, INTENT(IN) :: LUNWRAP
      COMPLEX*8, INTENT(OUT) :: RECV(MDIM,MFREQ,*)
!.... local variables 
      COMPLEX*16, ALLOCATABLE :: ESTR(:,:) 
      COMPLEX*16 ZNUM, ZZERO, U, V, W, D
      COMPLEX*8 RRFNEZ(3), COSAZ, SINAZ, CZERO
      REAL*8 RDEN, ROBS, QOBS, REST, QEST, ADATA, PDATA, AEST, PEST,
     ;       OAMP, OPHS, EAMP, EPHS, DPH, PI180, UMAG, VMAG, WMAG, 
     ;       UPHASE, VPHASE, WPHASE  
!     REAL*8 DWRAPPH
      INTEGER*4 IFREQ, ISRC, IREC, I
      LOGICAL*4 LDEAD 
      PARAMETER(ZZERO = DCMPLX(0.D0,0.D0), CZERO = CMPLX(0.0,0.0))
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
!
!----------------------------------------------------------------------!
!
!.... initialize rotation angles
      ALLOCATE(ESTR(NDIM,NSRC))
      COSAZ = CMPLX(DCOS(AZMOD*PI180),0.D0)
      SINAZ = CMPLX(DSIN(AZMOD*PI180),0.D0)
!.... loop on frequencies
      DO 100 IFREQ=1,NFREQ 
!
!....... loop on receivers
         DO 200 IREC=1,NREC 
!
!.......... rotate estaimtes
            DO 201 ISRC=1,NSRC
!
!............. rotate estimates counter-clockwise
               IF (LUNWRAP) THEN
                  UMAG = DBLE(REAL(EST(1,IFREQ,IREC,ISRC)))
                  VMAG = DBLE(REAL(EST(2,IFREQ,IREC,ISRC)))
                  WMAG = DBLE(REAL(EST(3,IFREQ,IREC,ISRC)))
                  UPHASE = DBLE(IMAG(EST(1,IFREQ,IREC,ISRC)))
                  VPHASE = DBLE(IMAG(EST(2,IFREQ,IREC,ISRC)))
                  WPHASE = DBLE(IMAG(EST(3,IFREQ,IREC,ISRC)))
                  U = DCMPLX(UMAG,0.D0)*CDEXP(DCMPLX(0.D0,UPHASE))
                  V = DCMPLX(VMAG,0.D0)*CDEXP(DCMPLX(0.D0,VPHASE))
                  W = DCMPLX(WMAG,0.D0)*CDEXP(DCMPLX(0.D0,WPHASE))
               ELSE
                  U = DCMPLX(EST(1,IFREQ,IREC,ISRC))
                  V = DCMPLX(EST(2,IFREQ,IREC,ISRC))
                  W = DCMPLX(EST(3,IFREQ,IREC,ISRC))
               ENDIF
               ESTR(1,ISRC) = U*COSAZ - V*SINAZ
               ESTR(2,ISRC) = U*SINAZ + V*COSAZ
               ESTR(3,ISRC) = W 
  201       CONTINUE !Loop on sources
!
!.......... loop on components
            LDEAD = .FALSE.
            DO 300 I=1,NDIM 
!
!............. loop on sources 
               ZNUM = ZZERO
               RDEN = 0.D0
               DO 400 ISRC=1,NSRC 
                  IF (OBS(I,IFREQ,IREC,ISRC).EQ.CZERO) THEN
                     GOTO 450 !nothing here
                  ENDIF
!
!................ estimates are unwrapped
                  IF (LUNWRAP) THEN
                     OPHS = DBLE(IMAG(OBS(I,IFREQ,IREC,ISRC)))
                     EPHS = DBLE(IMAG(EST(I,IFREQ,IREC,ISRC))) 
                     IF (IRESTP.EQ.1) THEN 
                        OAMP = 1.D0
                        EAMP = 1.D0
                     ELSE
                        OAMP = DBLE(REAL(OBS(I,IFREQ,IREC,ISRC)))
                        EAMP = CDABS(ESTR(I,ISRC))
                     ENDIF
                     !DPH = DWRAPPH(0,EPHS - OPHS)
                     DPH =-EPHS + OPHS !u*w
                     ZNUM = ZNUM + DCMPLX(EAMP*OAMP,0.D0)
     ;                            *CDEXP(DCMPLX(0.D0,DPH))
                     RDEN = RDEN + EAMP**2
                  ELSE !estimates are complex
                     ROBS = DBLE(REAL(OBS(I,IFREQ,IREC,ISRC)))
                     QOBS = DBLE(IMAG(OBS(I,IFREQ,IREC,ISRC)))
                     REST = DREAL(ESTR(I,ISRC))
                     QEST = DIMAG(ESTR(I,ISRC))
                     IF (IRESTP.EQ.1 .OR. IRESTP.EQ.2) THEN
                        IF (IRESTP.EQ.2) THEN !amplitude only
                           ADATA = DSQRT(ROBS**2 + QOBS**2)
                           PDATA = 0.D0
                           AEST  = DSQRT(REST**2 + QEST**2)
                           PEST  = 0.D0
                        ELSE !phase is default here
                           ADATA = 1.D0
                           PDATA = DATAN2(QOBS,ROBS)
                           AEST  = 1.D0
                           PEST  = DATAN2(QEST,REST)
                        ENDIF
                        D = DCMPLX(ADATA,0.D0)*CDEXP(DCMPLX(0.D0,PDATA))
                        U = DCMPLX(AEST ,0.D0)*CDEXP(DCMPLX(0.D0,PEST ))
                     ELSE !phase and amplitude 
                        D = DCMPLX(OBS(I,IFREQ,IREC,ISRC))
                        U = DCMPLX(ESTR(I,ISRC))
                     ENDIF
                     ZNUM = ZNUM + DCONJG(U)*D
                     RDEN = RDEN + DREAL(DCONJG(U)*U) !|u|^2 = conj(u)*u
                  ENDIF !end check on residual packing
  450             CONTINUE !break ahead, no data 
  400          CONTINUE !loop on sources 
               IF (RDEN.GT.0.D0) THEN
                  RRFNEZ(I) = CMPLX(ZNUM/DCMPLX(RDEN,0.D0))
               ELSE
                  LDEAD = .TRUE.
                  WRITE(*,*) 'recupd: Receiver',IREC,'component',I,
     ;                       'is dead'
                  WRITE(*,*) '        Setting response to zero'
                  RRFNEZ(I) = CMPLX(0.0,0.0)
               ENDIF  
  300       CONTINUE !loop on components
!
!.......... put back into (u,v,w)
            RECV(1,IFREQ,IREC) = RRFNEZ(1)*COSAZ + RRFNEZ(2)*SINAZ
            RECV(2,IFREQ,IREC) =-RRFNEZ(1)*SINAZ + RRFNEZ(2)*COSAZ
            RECV(3,IFREQ,IREC) = RRFNEZ(3)
  200    CONTINUE !loop on receivers 
  100 CONTINUE !loop on frequencies
      DEALLOCATE(ESTR) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE WTRECST(PROJNM, MDIM,MFREQ, NDIM,NFREQ,NREC, 
     ;                   IBLOCK,ITER, AZMOD, FREQ,RECV, IERR) 
!
!     Writes the receiver static corrections for this inversion blcok
!     and iteration
!
!     INPUT      MEANING
!     -----      ------- 
!     AZMOD      model azimtuth
!     FREQ       frequency list
!     IBLOCK     block number in inversion
!     ITER       iteration in inversion
!     MDIM       leading dimension for RECV
!     MFREQ      leading dimension for RECV
!     NDIM       number of components on receiver
!     NFREQ      number of frequencies
!     NREC       number of receivers 
!     PROJNM     project name
!     RECV       receiver statics in (u,v,w) frame 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag 
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: PROJNM
      COMPLEX*8, INTENT(IN) :: RECV(MDIM,MFREQ,*) 
      REAL*8, INTENT(IN) :: FREQ(NFREQ), AZMOD 
      INTEGER*4, INTENT(IN) :: MDIM, MFREQ, NDIM, NFREQ, NREC, IBLOCK, 
     ;                         ITER
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(80) FLNAME 
      CHARACTER(5) CITER
      CHARACTER(3) CBLOCK
      CHARACTER(4), ALLOCATABLE :: CDAT(:)  
      COMPLEX*8 N, E, Z
      REAL*8 PI180 
      REAL*4 R, Q
      INTEGER*4 NBYTES, LWORK, MYEND, INDX, IFREQ, IREC, I, J, IUNIT
      LOGICAL*4 LEX, LISDIR, LSWAP
      CHARACTER(4) PACKI4, PACKR4
      INTEGER*4 ENDIAN 
      PARAMETER(IUNIT = 68) 
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
!
!----------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0 
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE. 
!.... allocate space and pack header
      NBYTES = 12 + 4*NFREQ + 8*NDIM*NFREQ*NREC 
      LWORK = NBYTES/4
      ALLOCATE(CDAT(LWORK))
      CDAT(1) = PACKI4(LSWAP,NFREQ)
      CDAT(2) = PACKI4(LSWAP,NREC)
      CDAT(3) = PACKI4(LSWAP,NDIM) 
      INDX = 4
      DO 1 IFREQ=1,NFREQ
         CDAT(INDX) = PACKR4(LSWAP,SNGL(FREQ(IFREQ)))
         INDX = INDX + 1
         DO 2 IREC=1,NREC
            !rotate (u,v) counterclockise to (N,E)
            CALL CROTATE(SNGL(AZMOD),
     ;                   RECV(1,IFREQ,IREC),RECV(2,IFREQ,IREC),N,E)
            Z = RECV(3,IFREQ,IREC) 
            DO 3 I=1,NDIM
               IF (I.EQ.1) THEN
                  R = REAL(N)
                  Q = IMAG(N) 
               ELSEIF (I.EQ.2) THEN
                  R = REAL(E)
                  Q = IMAG(E)
               ELSE 
                  R = REAL(Z)
                  Q = IMAG(Z) 
               ENDIF
               CDAT(INDX) = PACKR4(LSWAP,R)
               INDX = INDX + 1
               CDAT(INDX) = PACKR4(LSWAP,Q)
               INDX = INDX + 1
    3       CONTINUE !loop on components
    2    CONTINUE !loop on sources
    1 CONTINUE !loop on frequencies
      INDX = INDX - 1
      IF (INDX*4.NE.NBYTES) THEN
         WRITE(*,*) 'ststf: Error this file is the wrong size!'
         IERR = 1
      ENDIF
!
!.... file handling
      LEX = LISDIR('./recrsp')
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./recrsp')
!.... set file name and write
      FLNAME(1:80) = ' ' 
      CITER(1:5) = ' ' 
      CBLOCK(1:3) = ' ' 
      WRITE(CITER,'(I5)') ITER
      WRITE(CBLOCK,'(I3)') IBLOCK
      CITER = ADJUSTL(CITER)
      CBLOCK = ADJUSTL(CBLOCK) 
      FLNAME = './recrsp/'//TRIM(PROJNM)//'-'//TRIM(CBLOCK)//'-'//
     ;         TRIM(CITER)//'.wrec'
      FLNAME = ADJUSTL(FLNAME)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED', 
     ;     STATUS='REPLACE',ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      WRITE(IUNIT,REC=1,IOSTAT=IERR) (CDAT(J),J=1,LWORK) 
      CLOSE(IUNIT)
      DEALLOCATE(CDAT) 
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE MODREC(MDIM,MFREQ, NDIM,NFREQ,NREC, IRESTP, RECV) 
!
!     Modifies the receiver response functions if inverting for phase 
!     or amplitude only.  This is a superfluous task.
!
!     INPUT      MEANING
!     -----      -------
!     IRESTP     =1 -> Phase only 
!                =2 -> Amplitude only 
!                =3 -> phase and amplitude (default)
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     NDIM       number of components 
!     NFREQ      number of frequencies
!     NREC       number of receivers
!
!     OUTPUT     MEANING
!     ------     -------
!     RECV       modified receiver response function
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(INOUT) :: RECV(MDIM,MFREQ,*) 
      INTEGER*4, INTENT(IN) :: MDIM, MFREQ, NFREQ, NREC, NDIM, IRESTP 
!.... local variables
      COMPLEX*8 CPHM2CM 
      REAL*4 SPHASE, PHASE, RMAG 
      INTEGER*4 IFREQ, IREC, I 
!
!----------------------------------------------------------------------!
!
      IF (IRESTP.NE.1 .AND .IRESTP.NE.2) RETURN !Nothing to do
      DO 1 IFREQ=1,NFREQ
         DO 2 IREC=1,NREC 
            DO 3 I=1,NDIM
               RMAG  =   CABS(RECV(I,IFREQ,IREC))  
               PHASE = SPHASE(RECV(I,IFREQ,IREC))  
               IF (RMAG.GT.0.0) THEN
                  IF (IRESTP.EQ.1) RMAG  = 1.0 !phase only inversion
                  IF (IRESTP.EQ.2) PHASE = 0.0 !amplitude only inversion
                  RECV(I,IFREQ,IREC) = CPHM2CM(RMAG,PHASE)
               ELSE
                  RECV(I,IFREQ,IREC) = CPHM2CM(0.0 ,0.0)
               ENDIF
   3        CONTINUE !loop on components
    2    CONTINUE
    1 CONTINUE  
      RETURN
      END
