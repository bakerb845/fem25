!
!     This is a testing utility.  It generates the 1D analytic solutions at the 
!     receiver DOFs for the mesh.  -  B. Baker February 2013
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO)  MSH  !mesh, model parameters
      TYPE (MOD1D_INFO) M1D  !1D models
      TYPE (SRC_INFO)   SRC  !source information
      TYPE (RECV_INFO)  RCV  !receiver info
      TYPE (FRQ_INFO)   FRQ  !frequency list
      COMPLEX*8, ALLOCATABLE :: EXACT(:,:,:,:), EST_LOC(:,:,:,:) 
      REAL*8, ALLOCATABLE :: XREC(:), ZREC(:), FREQ_LOC(:) 
      REAL*4, ALLOCATABLE :: WIGGLE(:,:,:,:), PRTAB(:,:), PLTAB(:,:), FREQ4(:), PERIOD(:)
      LOGICAL*4, ALLOCATABLE :: LFSURF(:)
      CHARACTER(80) PROJNM, FILENM, EXFL, TMPDIR 
      REAL*8 PW(20), CRW(20), CLW(20)
      COMPLEX*16 CARG
      COMPLEX*8 CONE
      REAL*8 VBASE, PI180, DT8, START8, OMEGA, TWOPI
      REAL*4 DT4_SRF, DT4_BDY, START4_SRF, START4_BDY
      INTEGER*4 ISRC,I,IREC,IFREQ,IDIR,IANS,NSAMP_SRF,NSAMP_BDY,NSRCIN, &
                JFREQ, JSRC, IERR
      PARAMETER(TWOPI = 6.2831853071795862D0)
      PARAMETER(PI180 = 0.017453292519943295D0)
      PARAMETER(CONE = CMPLX(1.0,0.0))
      REAL*8 CVBASE
!.... MPI
      INCLUDE 'mpif.h'
      REAL*8 TSSIM,TESIM
      INTEGER*4 MASTER,MYID,MPIERR,  NPROCS,NPGROUPS,NPARTS
      PARAMETER(MASTER = 0)
!.... mesh variables
      LOGICAL*4 LNEW, LRECEX
      INTEGER*4 NABS,NELEME
!.... miscellaneous
      LOGICAL*4 LNSRF
      LOGICAL*4 LFILES !unused
!
!----------------------------------------------------------------------------------------!
!
!.... initialize MPI 
      CALL MPI_INIT(MPIERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, MPIERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, MPIERR)
!----------------------------------------------------------------------------------------!
!                    Read information on mesh and model                                  !
!----------------------------------------------------------------------------------------!
!
!.... have head node read model information
      TSSIM = MPI_WTIME() 
      IERR = 0
      IF (MYID == MASTER) THEN
         WRITE(*,9407) 
 9407    FORMAT(' --------------------------------------------------------------',/, &
                ' -   xexact1d: A utility for calculating analytic solutions   -',/, &
                ' -             at receivers                                   -',/, &
                ' --------------------------------------------------------------',/) 
         PROJNM(1:80) = ' '
         WRITE(*,*) 'xexact1d: Enter project name:'
         READ(*,'(A)') PROJNM
         PROJNM = ADJUSTL(PROJNM)         
         NPGROUPS = 1
!
!....... get variables from the .ini files
         CALL READ_FWD_INI(PROJNM,.TRUE., TMPDIR,LFILES,LNSRF, MSH, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xbielak25: Error I cannot locate your spec file!'
            GOTO 500 
         ENDIF
!        msh%AZTOL = 10.D0
!        msh%AOITOL = 5.D0
!        msh%FREQ0 = 20.D0
         NPARTS = NPROCS/NPGROUPS
         WRITE(*,*) 'xexact25: Reading mesh...'
         CALL RDMESHBK(PROJNM, NELEME,NABS, MSH, IERR) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xexact25: Error reading mesh files!'
            GOTO 500 
         ENDIF 

         WRITE(*,*) 'xbielak25: Setting interpolation points...'
         ALLOCATE(msh%XIPTS (msh%NLXI))
         ALLOCATE(msh%ETAPTS(msh%NLETA)) 
         CALL GENPTS(msh%IITYPE, msh%NLXI,msh%NLETA, msh%XIPTS,msh%ETAPTS) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xexact25: Error reading mesh files!'
            GOTO 500 
         ENDIF 

         CALL PLOT_ELMOD_VTK(PROJNM,NGNOD,msh%NNPG,msh%NNPG, msh%NELEM, msh%LISISO, &
                             msh%IENG,msh%XLOCS,msh%ZLOCS, msh%DENS,msh%ECOEFF)
!
!....... read the frequency list
         CALL RDFREQ(PROJNM, FRQ,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xexact25: Error reading frequency list'
            GOTO 500
         ENDIF
!
!....... read the source list
         CALL RDSRC_EQ(PROJNM, msh%XLATMIN,msh%XLONMIN,msh%XLATMAX,msh%XLONMAX, &
                       msh%AZMOD, SRC,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xexact25: Error calling rdsrc_eq'
            GOTO 500
         ENDIF
!
!....... check for a receiver file
         CALL RDREC_HD(PROJNM, rcv%NREC,IERR)
         IF (IERR /= 0) THEN
            IF (IERR > 0) THEN
               WRITE(*,*) 'xexact25: Error reading recv header'
               GOTO 500
            ELSE
               WRITE(*,*) 'xexact25: Error there is no receiver file!'
               GOTO 500
            ENDIF
         ELSE !read the receiver file
            ALLOCATE(LFSURF  (rcv%NREC))
            ALLOCATE(    XREC(rcv%NREC))
            ALLOCATE(rcv%YREC(rcv%NREC))
            ALLOCATE(    ZREC(rcv%NREC))
            CALL RDREC(PROJNM, rcv%NREC, LFSURF,XREC,rcv%YREC,ZREC, IERR)
         ENDIF
!
!....... get the body wave freuqency info
         DT4_SRF = 0.0
         START4_SRF = 0.0
         IF (frq%NFREQ_SRF > 0 .AND. src%NSRC_SRF > 0) THEN
            FILENM(1:80) = ' '
            FILENM = TRIM(PROJNM)//'_srf'
            FILENM = ADJUSTL(FILENM)
            CALL RDSRCT_HD(FILENM, NSRCIN,NSAMP_SRF,DT8,START8,IERR) 
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xexact25: Could not read .srcts file'
               WRITE(*,*) 'xexact25: Enter (1) for a time domain solution'
               READ *, IANS
               IF (IANS == 1) THEN 
                  WRITE(*,*) 'xexact25: Enter the number of samples'
                  READ *, NSAMP_SRF
                  WRITE(*,*) 'xexact25: Enter the sampling period in seconds'
                  READ *, DT4_SRF
               ELSE
                  NSAMP_SRF = 0 
                  DT4_SRF = 0.0
               ENDIF
            ELSE
               DT4_SRF = REAL(DT8)
               START4_SRF = REAL(START8) 
            ENDIF
         ENDIF
!
!....... get the body wave freuqency info
         DT4_BDY = 0.0 
         START4_BDY = 0.0 
         IF (frq%NFREQ_BDY > 0 .AND. src%NSRC_BDY > 0) THEN
            FILENM(1:80) = ' '
            FILENM = TRIM(PROJNM)//'_bdy'
            FILENM = ADJUSTL(FILENM)
            CALL RDSRCT_HD(FILENM, NSRCIN,NSAMP_BDY,DT8,START8,IERR) 
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xexact25: Could not read .srcts file'
               WRITE(*,*) 'xexact25: Enter (1) for a time domain solution'
               READ *, IANS
               IF (IANS == 1) THEN 
                  WRITE(*,*) 'xexact25: Enter the number of samples'
                  READ *, NSAMP_BDY
                  WRITE(*,*) 'xexact25: Enter the sampling period in seconds'
                  READ *, DT4_BDY
               ELSE
                  NSAMP_BDY = 0 
                  DT4_BDY = 0.0 
               ENDIF
            ELSE
               DT4_BDY = REAL(DT8)
               START4_BDY = REAL(START8) 
            ENDIF
         ENDIF
!
!....... check for a source file
         ALLOCATE(src%SOURCE(frq%NFREQ,src%NSRC))
         CALL SRCSUB_SB(PROJNM,frq%NFREQ, frq%NFREQ,src%NSRC,                   &   
                        frq%NFREQ_SRF,frq%NFREQ_BDY,src%NSRC_SRF,src%NSRC_BDY,  &
                        frq%CFTYPE,src%SRCTYP, frq%FREQ, src%SOURCE,IERR)
         IF (IERR /= 0) THEN
            IF (IERR < 0) THEN
               WRITE(*,*) 'xexact25: Warning Greens functions are unity'
               WRITE(*,*) '          You can convolve an appropriate STF later'
               IERR = 0
            ELSE
               WRITE(*,*) 'xexact25: Error calling srcsub'
               GOTO 500 
            ENDIF
         ENDIF
!
!....... read the frequency response list
         ALLOCATE(rcv%RECV(NDIM,frq%NFREQ,rcv%NREC))
         CALL RDREC_RESP(PROJNM, NDIM,frq%NFREQ,NDIM,frq%NFREQ,rcv%NREC, msh%AZMOD, &
                         frq%FREQ, LRECEX,rcv%RECV)
         IF (.NOT.LRECEX) THEN
            WRITE(*,*) 'xexact25: Receiver response functions are set to unity'
            IERR = 0
         ENDIF
!
!....... generate the 1D models
         IF (MYID == MASTER) WRITE(*,*) 'xexact25: Generating 1D models...'
         CALL FILL1D(MSH,M1D) 
      ENDIF 
!----------------------------------------------------------------------------------------!
!                        Graph generation                                                !
!----------------------------------------------------------------------------------------!
!
!.... generate a graph
      IF (MYID == MASTER) THEN
         WRITE(*,*) 'xexact25: Generating graph...'
         CALL GEN_GRAPH25(.TRUE.,NPARTS, LFSURF,XREC,ZREC, MSH,RCV,IERR) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xexact25: Error generating graph!'
            GOTO 500
         ENDIF
         WRITE(*,*)
         WRITE(*,9408) msh%NORD,msh%NELEM,NABS,NELEME, &
                       msh%NNPG,msh%NNPE,msh%NDOF,msh%NZERO  
 9408    FORMAT(' xexact25: Polynomial order:'                    ,I4 ,/,        &
                '           Number of elements in mesh:'          ,I10,/,         &
                '           Number of absorbing elements:'        ,I8 ,/,         &
                '           Number of Bielak elements:'           ,I8 ,/,         &
                '           Number of anchor nodes in mesh:'      ,I10,/,         & 
                '           Number of nodes in Bielak boundary:'  ,I10,/,         &
                '           Number of degrees of freedom:'        ,I14,/,         &
                '           Number of non-zeros in global matrix:',I16,/)
!
!....... now i can delete a bunch of stuff
         DEALLOCATE(msh%CNNPG)
         DEALLOCATE(msh%CDOMAIN)
         DEALLOCATE(msh%CNP) 
         DEALLOCATE(msh%ECOEFF)
         DEALLOCATE(msh%DENS)
         DEALLOCATE(msh%XD)
         DEALLOCATE(msh%ZD)
         DEALLOCATE(msh%XLOCS)
         DEALLOCATE(msh%ZLOCS)
         DEALLOCATE(msh%XLOCSE)  
         DEALLOCATE(msh%ZLOCSE)
         DEALLOCATE(msh%XIPTS)
         DEALLOCATE(msh%ETAPTS)
         DEALLOCATE(msh%LM)
         DEALLOCATE(msh%IENG)
         DEALLOCATE(msh%IDOFSE)
         DEALLOCATE(msh%IRPTR)
         DEALLOCATE(msh%JCPTR)
         DEALLOCATE(msh%PART)
         DEALLOCATE(LFSURF)
      ENDIF
      IF (MYID == MASTER) WRITE(*,*) 'xexact25: Broadcasting 1D models...'
      CALL BCAST_MOD1D_INFO(MYID,MPI_COMM_WORLD,MASTER, MSH,M1D)
      IF (MYID == MASTER) WRITE(*,*) 'xexact25: Broadcasting frequency information...'
      CALL BCAST_FRQ_INFO(MYID,MPI_COMM_WORLD,MASTER, FRQ)
      IF (MYID == MASTER) WRITE(*,*) 'xexact25: Broadcasting source details...'
      CALL BCAST_SRC_INFO(MYID,MPI_COMM_WORLD,MASTER, SRC)
      IF (MYID == MASTER) WRITE(*,*) 'xexact25: Broadcasting receiver response...'
      CALL BCAST_RCV_INFO(MYID,MPI_COMM_WORLD,MASTER, RCV)
      IF (.NOT.ALLOCATED(XREC)) ALLOCATE(XREC(rcv%NREC))
      IF (.NOT.ALLOCATED(ZREC)) ALLOCATE(ZREC(rcv%NREC))
      CALL MPI_BCAST(XREC,rcv%NREC,MPI_DOUBLE_PRECISION, MASTER,MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(ZREC,rcv%NREC,MPI_DOUBLE_PRECISION, MASTER,MPI_COMM_WORLD,MPIERR)
      CALL BCAST_RCV_RESP(MYID,MPI_COMM_WORLD,MASTER, .TRUE.,frq%NFREQ, RCV)
!
!.... generate greens functions 
      IF (MYID == MASTER) THEN
!
!....... set space and null out exactx solution 
         ALLOCATE(EXACT(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
         DO 50 IFREQ=1,frq%NFREQ
            DO 51 ISRC=1,src%NSRC
               DO 52 I=1,NDIM
                  EXACT(I,IFREQ,1:rcv%NREC,ISRC) = CMPLX(0.0,0.0)
    52         CONTINUE  
    51      CONTINUE
    50   CONTINUE 
!
!....... loop on frequencies
         WRITE(*,*) 'xexact1d: Generating exact solution'
         ALLOCATE(src%PYTAB(frq%NFREQ,src%NSRC))
         LNEW = .TRUE.
         DO 100 IFREQ=1,frq%NFREQ
            DO 101 ISRC=1,src%NSRC
               IF (frq%CFTYPE(IFREQ) == 'S' .AND. src%SRCTYP(ISRC)(1:1) == 'P') GOTO 110
               IF (frq%CFTYPE(IFREQ) == 'B' .AND. src%SRCTYP(ISRC)(1:1) == 'S') GOTO 110
               CALL GENEX(TMPDIR,LNSRF, &
                          src%SRCTYP(ISRC), NDIM, m1d%NL1D_LT,m1d%NL1D_RT, &
                          rcv%NREC, src%CSIDE(ISRC),  &
                          ISRC,MYID, LNEW, src%MODE(ISRC), &
                          msh%FREQ0,frq%FREQ(IFREQ),src%AOI(ISRC),src%BAZN(ISRC), msh%XMOD0,msh%XMOD1, &
                          msh%XLATMIN,msh%XLONMIN, msh%XLATMAX,msh%XLONMAX, &
                          src%SLAT(ISRC),src%SLON(ISRC),src%SDEP(ISRC),src%SMAG(ISRC), &
                          src%STRIKE(ISRC),src%DIP(ISRC),src%RAKE(ISRC), &
                          XREC,rcv%YREC, &
                          m1d%VP1D_LT, m1d%VS1D_LT, m1d%RH1D_LT,m1d%Z1D_LT,   &   
                          m1d%VP1D_RT, m1d%VS1D_RT, m1d%RH1D_RT,m1d%Z1D_RT,   &   
                          m1d%VPD_RLLT,m1d%VSD_RLLT,m1d%ROD_RLLT,m1d%HDD_RLLT,   &   
                          m1d%VPD_LVLT,m1d%VSD_LVLT,m1d%ROD_LVLT,m1d%HDD_LVLT,  &   
                          m1d%VPD_RLRT,m1d%VSD_RLRT,m1d%ROD_RLRT,m1d%HDD_RLRT, &   
                          m1d%VPD_LVRT,m1d%VSD_LVRT,m1d%ROD_LVRT,m1d%HDD_LVRT, &
                          m1d%QP1D_LT, m1d%QP1D_RT, m1d%QS1D_LT, m1d%QS1D_RT,  &
                          CONE, PW,CRW,CLW,  &
                          EXACT(1:NDIM,IFREQ,1:rcv%NREC,ISRC),IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xexact1d: Error calling genex!'
                  CALL MPI_ABORT(MPI_COMM_WORLD,30,MPIERR)
               ENDIF
               LNEW = .FALSE.
               IF (src%SRCTYP(ISRC)(1:1) == 'S') THEN
                  IF (src%SRCTYP(ISRC)(1:1) == 'L') THEN
                     IF (CLW(1) == 0.D0) THEN
                        PRINT *, 'ERR1'
                        CALL MPI_ABORT(MPI_COMM_WORLD,30,MPIERR)
                     ENDIF
                     src%PYTAB(IFREQ,ISRC) = DSIN(src%BAZN(ISRC)*PI180)/(CLW(1)*1.D3)
                  ELSE
                     IF (CRW(1) == 0.D0) THEN
                        PRINT *, 'ERR2' 
                        CALL MPI_ABORT(MPI_COMM_WORLD,30,MPIERR)
                     ENDIF
                     src%PYTAB(IFREQ,ISRC) = DSIN(src%BAZN(ISRC)*PI180)/(CRW(1)*1.D3)
                  ENDIF
               ELSE
                  VBASE = CVBASE(src%SRCTYP(ISRC),src%CSIDE(ISRC), m1d%NL1D_LT,m1d%NL1D_RT,   &
                                 m1d%VP1D_LT,m1d%VS1D_LT, m1d%VP1D_RT,m1d%VS1D_RT)  
                  src%PYTAB(IFREQ,ISRC) = DSIN(src%BAZN(ISRC)*PI180)*DSIN(src%AOI(ISRC)*PI180)/VBASE 
               ENDIF
  110          CONTINUE
  101       CONTINUE  !loop on sources 
  100    CONTINUE  !loop on frequencies 
!
!....... extract the py table
         ALLOCATE(PERIOD(frq%NFREQ))
         ALLOCATE(PRTAB(frq%NFREQ,src%NSRC))
         ALLOCATE(PLTAB(frq%NFREQ,src%NSRC))
!        VFAST = MAXVAL(VP1D_RT) 
!        print *, NFREQ, NFREQ,NSRC,NL1D_LT,NL1D_RT, NPROCS
!        print *, mode
!        print *, bazn
!        print *, aoi
!        print *, vfast 
!        CALL CONCAT_PYTAB(NFREQ, NFREQ,NSRC,NL1D_LT,NL1D_RT, NPROCS, &
!                          CSIDE,SRCTYP, MODE, VP1D_LT,VS1D_LT,VP1D_RT,VS1D_RT, &
!                          BAZN,AOI, VFAST,PERIOD,PRTAB,PLTAB,IERR)
!        ALLOCATE(PYTAB(NFREQ,NSRC)) 
         DO 200 IFREQ=1,frq%NFREQ
            DO 201 ISRC=1,src%NSRC
!              IF (SRCTYP(ISRC)(1:1) == 'P') THEN
!                 PYTAB(IFREQ,ISRC) = PRTAB(IFREQ,ISRC)
!              ELSE 
!                 IF (SRCTYP(ISRC)(2:2) == 'B' .OR. &
!                     SRCTYP(ISRC)(2:2) == 'V' .OR. &
!                     SRCTYP(ISRC)(2:2) == 'R') THEN 
!                    PYTAB(IFREQ,ISRC) = PRTAB(IFREQ,ISRC)
!                    IF (SRCTYP(ISRC)(2:2) == 'B') THEN 
!                       WRITE(*,*) 'gen_grns: Dont know what to do here, ask steve'
!                       WRITE(*,*) '          Defaulting to Rayleigh pytab'
!                    ENDIF 
!                 ELSE 
!                    PYTAB(IFREQ,ISRC) = PLTAB(IFREQ,ISRC)
!                 ENDIF
!              ENDIF
!
!............. phase shift and convolve instrument response and stf
               !print *, ifreq,isrc, src%source(ifreq,isrc)
               DO 202 IREC=1,rcv%NREC
                  OMEGA = TWOPI*frq%FREQ(IFREQ)
                  CARG = CDEXP(DCMPLX(0.D0,+OMEGA*src%PYTAB(IFREQ,ISRC)*rcv%YREC(IREC)))
                  DO 203 I=1,NDIM
                     EXACT(I,IFREQ,IREC,ISRC) = EXACT(I,IFREQ,IREC,ISRC)*CMPLX(CARG)
                     EXACT(I,IFREQ,IREC,ISRC) = EXACT(I,IFREQ,IREC,ISRC)      &
                                               *rcv%RECV(I,IFREQ,IREC)        &
                                               *src%SOURCE(IFREQ,ISRC)
  203             CONTINUE
  202          CONTINUE !loop on receivers
!              if (imag(exact(1,ifreq,1,1)) > 0.0) then
!                 write(*,844) exact(1,ifreq,1,1)
!              else
!                 write(*,843) real(exact(1,ifreq,1,1)),abs(imag(exact(1,ifreq,1,1)))
!              endif
!844           format(' ',f10.7,' +',f10.7,'i')
!843           format(' ',f10.7,' -',f10.7,'i')
  201       CONTINUE !loop on sources
  200    CONTINUE !lop o nfrequencies
!        DEALLOCATE(PERIOD)
!        DEALLOCATE(PRTAB) 
!        DEALLOCATE(PLTAB) 
!        WRITE(*,*) 'xexact1d: Writing exact solution'
!        EXFL(1:80) = ' '
!        EXFL = TRIM(ADJUSTL(PROJNM))//'_exact'
!        EXFL = ADJUSTL(EXFL)
!        CALL WTEEST25(EXFL, frq%NFREQ,rcv%NREC,frq%NFREQ, rcv%NREC,src%NSRC,            &
!                      frq%FREQ, EXACT(1,:,:,:),EXACT(2,:,:,:),EXACT(3,:,:,:), IERR)
         IF (frq%NFREQ_SRF > 0 .AND. src%NSRC_SRF > 0) THEN
            WRITE(*,*) 'xexact1d: Inverse transforming surface waves...'
            ALLOCATE(EST_LOC(NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF))
            EST_LOC(:,:,:,:) = CMPLX(0.0,0.0) 
            ALLOCATE(FREQ_LOC(frq%NFREQ_SRF))
            JFREQ = 0
            DO IFREQ=1,frq%NFREQ
               IF (frq%CFTYPE(IFREQ) == 'S') THEN 
                  JFREQ = JFREQ + 1
                  FREQ_LOC(JFREQ) = SNGL(frq%FREQ(IFREQ))
                  JSRC = 0
                  DO ISRC=1,src%NSRC
                     IF (src%SRCTYP(ISRC)(1:1) == 'S') THEN
                        JSRC = JSRC + 1
                        DO IREC=1,rcv%NREC
                           DO I=1,NDIM
                              EST_LOC(I,JFREQ,IREC,JSRC) = EXACT(I,IFREQ,IREC,ISRC) 
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            IDIR = 1 !inverse F.T.
            EXFL(1:80) = ' '
            EXFL = TRIM(ADJUSTL(PROJNM))//'_exact_srf'
            EXFL = ADJUSTL(EXFL)
            CALL WTEEST25(EXFL, frq%NFREQ_SRF,rcv%NREC,frq%NFREQ_SRF,  &   
                          rcv%NREC,src%NSRC_SRF, FREQ_LOC,   &   
                          EST_LOC(1,:,:,:),EST_LOC(2,:,:,:),EST_LOC(3,:,:,:), IERR)
            IF (NSAMP_SRF > 0) THEN
               ALLOCATE(FREQ4(frq%NFREQ_SRF))
               ALLOCATE(WIGGLE(NDIM,NSAMP_SRF,rcv%NREC,src%NSRC_SRF))
               FREQ4(1:frq%NFREQ_SRF) = SNGL(FREQ_LOC(1:frq%NFREQ_SRF))
               CALL DFTWIG4(NDIM,NSAMP_SRF,frq%NFREQ_SRF,rcv%NREC, &
                       NSAMP_SRF,NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF,IDIR,DT4_SRF,&
                       FREQ4,EST_LOC, WIGGLE) 
               print *, minval(wiggle), maxval(wiggle)
               CALL WTSEISM(EXFL, NDIM,NSAMP_SRF,rcv%NREC, &
                                  NDIM,NSAMP_SRF,rcv%NREC,src%NSRC_SRF,    &
                            DT4_SRF,START4_SRF, WIGGLE, IERR)
               IF (IERR /= 0) WRITE(*,*) 'xwiggle25: Error writing .pest file'
               DEALLOCATE(WIGGLE)
               DEALLOCATE(FREQ4)
            ENDIF
            DEALLOCATE(EST_LOC)
            DEALLOCATE(FREQ_LOC)
         ENDIF
!
!....... repeat for body waves
         IF (frq%NFREQ_BDY > 0 .AND. src%NSRC_BDY > 0) THEN
            WRITE(*,*) 'xexact1d: Inverse transforming body waves...'
            ALLOCATE(EST_LOC(NDIM,frq%NFREQ_BDY,rcv%NREC,src%NSRC_BDY))
            EST_LOC(:,:,:,:) = CMPLX(0.0,0.0)
            ALLOCATE(FREQ_LOC(frq%NFREQ_BDY))
            JFREQ = 0 
            DO IFREQ=1,frq%NFREQ
               IF (frq%CFTYPE(IFREQ) == 'B') THEN 
                  JFREQ = JFREQ + 1 
                  FREQ_LOC(JFREQ) = SNGL(frq%FREQ(IFREQ))
                  JSRC = 0 
                  DO ISRC=1,src%NSRC
                     IF (src%SRCTYP(ISRC)(1:1) == 'P') THEN
                        JSRC = JSRC + 1 
                        DO IREC=1,rcv%NREC
                           DO I=1,NDIM
                              EST_LOC(I,JFREQ,IREC,JSRC) = EXACT(I,IFREQ,IREC,ISRC) 
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            IDIR = 1 !inverse F.T.
            EXFL(1:80) = ' ' 
            EXFL = TRIM(ADJUSTL(PROJNM))//'_exact_bdy'
            EXFL = ADJUSTL(EXFL)
            CALL WTEEST25(EXFL, frq%NFREQ_BDY,rcv%NREC,frq%NFREQ_BDY,  &   
                          rcv%NREC,src%NSRC_BDY, FREQ_LOC,   &   
                          EST_LOC(1,:,:,:),EST_LOC(2,:,:,:),EST_LOC(3,:,:,:), IERR)
            IF (NSAMP_BDY > 0) THEN
               ALLOCATE(FREQ4(frq%NFREQ_BDY))
               ALLOCATE(WIGGLE(NDIM,NSAMP_BDY,rcv%NREC,src%NSRC_BDY))
               FREQ4(1:frq%NFREQ_BDY) = SNGL(FREQ_LOC(1:frq%NFREQ_BDY))
               CALL DFTWIG4(NDIM,NSAMP_BDY,frq%NFREQ_BDY,rcv%NREC, &
                       NSAMP_BDY,NDIM,frq%NFREQ_BDY,rcv%NREC,src%NSRC_BDY,IDIR,DT4_BDY,&
                       FREQ4,EST_LOC, WIGGLE) 
               CALL WTSEISM(EXFL, NDIM,NSAMP_BDY,rcv%NREC,  &
                                  NDIM,NSAMP_BDY,rcv%NREC,src%NSRC_BDY,    &   
                            DT4_BDY,START4_BDY, WIGGLE, IERR)
               IF (IERR /= 0) WRITE(*,*) 'xwiggle25: Error writing .pest file'
               DEALLOCATE(WIGGLE)
               DEALLOCATE(FREQ4)
            ENDIF
            DEALLOCATE(EST_LOC)
            DEALLOCATE(FREQ_LOC)
         ENDIF
      ENDIF
!
!.... clean receiver info
      IF (ALLOCATED(XREC)) DEALLOCATE(XREC)
      IF (ASSOCIATED(rcv%YREC)) DEALLOCATE(rcv%YREC) 
      IF (ALLOCATED(ZREC)) DEALLOCATE(ZREC)
!
!.... clean the 1D models
      IF (ASSOCIATED(rcv%RECV))     DEALLOCATE(rcv%RECV)
      IF (ASSOCIATED(m1d%VP1D_LT))  DEALLOCATE(m1d%VP1D_LT)
      IF (ASSOCIATED(m1d%VS1D_LT))  DEALLOCATE(m1d%VS1D_LT)
      IF (ASSOCIATED(m1d%RH1D_LT))  DEALLOCATE(m1d%RH1D_LT)
      IF (ASSOCIATED(m1d% Z1D_LT))  DEALLOCATE(m1d% Z1D_LT)

      IF (ASSOCIATED(m1d%VP1D_RT))  DEALLOCATE(m1d%VP1D_RT)
      IF (ASSOCIATED(m1d%VS1D_RT))  DEALLOCATE(m1d%VS1D_RT)
      IF (ASSOCIATED(m1d%RH1D_RT))  DEALLOCATE(m1d%RH1D_RT)
      IF (ASSOCIATED(m1d% Z1D_RT))  DEALLOCATE(m1d% Z1D_RT)

      IF (ASSOCIATED(m1d%VPD_RLLT)) DEALLOCATE(m1d%VPD_RLLT)
      IF (ASSOCIATED(m1d%VSD_RLLT)) DEALLOCATE(m1d%VSD_RLLT)
      IF (ASSOCIATED(m1d%ROD_RLLT)) DEALLOCATE(m1d%ROD_RLLT)
      IF (ASSOCIATED(m1d%HDD_RLLT)) DEALLOCATE(m1d%HDD_RLLT)

      IF (ASSOCIATED(m1d%VPD_RLRT)) DEALLOCATE(m1d%VPD_RLRT)
      IF (ASSOCIATED(m1d%VSD_RLRT)) DEALLOCATE(m1d%VSD_RLRT)
      IF (ASSOCIATED(m1d%ROD_RLRT)) DEALLOCATE(m1d%ROD_RLRT)
      IF (ASSOCIATED(m1d%HDD_RLRT)) DEALLOCATE(m1d%HDD_RLRT)

      IF (ASSOCIATED(m1d%VPD_LVLT)) DEALLOCATE(m1d%VPD_LVLT)
      IF (ASSOCIATED(m1d%VSD_LVLT)) DEALLOCATE(m1d%VSD_LVLT)
      IF (ASSOCIATED(m1d%ROD_LVLT)) DEALLOCATE(m1d%ROD_LVLT)
      IF (ASSOCIATED(m1d%HDD_LVLT)) DEALLOCATE(m1d%HDD_LVLT)

      IF (ASSOCIATED(m1d%VPD_LVRT)) DEALLOCATE(m1d%VPD_LVRT)
      IF (ASSOCIATED(m1d%VSD_LVRT)) DEALLOCATE(m1d%VSD_LVRT)
      IF (ASSOCIATED(m1d%ROD_LVRT)) DEALLOCATE(m1d%ROD_LVRT)
      IF (ASSOCIATED(m1d%HDD_LVRT)) DEALLOCATE(m1d%HDD_LVRT)
!.... clean source
      IF (ASSOCIATED(src%SOURCE)) DEALLOCATE(src%SOURCE)
      IF (ASSOCIATED(src%SRCTYP)) DEALLOCATE(src%SRCTYP)
      IF (ASSOCIATED(src%BAZN))   DEALLOCATE(src%BAZN)
      IF (ASSOCIATED(src%AOI))    DEALLOCATE(src%AOI)
      IF (ASSOCIATED(src%SLAT))   DEALLOCATE(src%SLAT)
      IF (ASSOCIATED(src%SLON))   DEALLOCATE(src%SLON)
      IF (ASSOCIATED(src%SDEP))   DEALLOCATE(src%SDEP)
      IF (ASSOCIATED(src%STRIKE)) DEALLOCATE(src%STRIKE)
      IF (ASSOCIATED(src%DIP))    DEALLOCATE(src%DIP)
      IF (ASSOCIATED(src%RAKE))   DEALLOCATE(src%RAKE)
      IF (ASSOCIATED(src%SMAG))   DEALLOCATE(src%SMAG)
      IF (ASSOCIATED(src%MODE))   DEALLOCATE(src%MODE)
      IF (ASSOCIATED(src%CSIDE))  DEALLOCATE(src%CSIDE)
      IF (ASSOCIATED(src%PYTAB))  DEALLOCATE(src%PYTAB) 
!.... clean frequencies
      IF (ASSOCIATED(frq%FREQ))   DEALLOCATE(frq%FREQ) 
!----------------------------------------------------------------------------------------!
!                 Error handling and memory cleaning                                     !
!----------------------------------------------------------------------------------------!
  500 CONTINUE !break ahead for errors
      IF (IERR /= 0) THEN
         WRITE(*,9406) MYID
 9406    FORMAT(' -------------------------------------------------------------',/, &
                ' -               An Error was Detected on Process ',I6,'     -',/, &
                ' -                  I will now abort the program             -',/, &
                ' -------------------------------------------------------------',/)
         CALL MPI_ABORT(MPI_COMM_WORLD,30,MPIERR)  
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
!
!.... all done 
      TESIM = MPI_WTIME()
      IF (MYID == MASTER) THEN
         WRITE(*,9405) (TESIM - TSSIM)/60.D0
 9405    FORMAT(' xexact25: Simulation time:',F8.2,' minutes')
      ENDIF 
      CALL MPI_FINALIZE(MPIERR)
      STOP
      END
