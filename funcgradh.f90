      SUBROUTINE FUNCGRADH25(MASTER,MYID,MYNID,IPGROUP,NPGROUPS,             &
                             MGCOMM,MYSLV_COMM,MYHD_COMM, LPGNEW,            &
                             VFAST, INV,MSH,SRC,RCV,FRQ, MID, IERR) 
!
!     Calculates the function, gradient, and pre-conditioner for the Hessian
!
!     INPUT      MEANING
!     -----      ------- 
!     AZMOD      model azimuth
!     IRESTP     residual type, 1 -> phase only, 2 -> amplitude 3 -> both (default) 
!     IPGROUP    process group, 0,1,2,...ngroups - 1
!     LFUNC      True -> Calculate an objective function
!     LGRAD      True -> Calculate a gradient
!     LPGNEW     True -> Calculate a diagonal Hessian gradient preconditioner
!     LGNEWT     True -> Gauss-Newton approximation, save J on first freq. group
!     MASTER     master process ID for communicators
!     MGCOMM     global communicator for all process
!     MID        MUMPS data structure
!     MYHD_COMM  communicator for head processes 
!     MYID       original process ID
!     MYNID      frequency group process ID
!     MYSLV_COMM communicator for head process and its slaves
!     NORM       residual norm - 1 -> L1, 2 -> default L2 
!     VFAST      a max velocity in 1D models for calculating py tolerances
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     FOBJ       objective function, through common
!     GRAD       gradient, through common
!     HESS       diagonal hessian, through common
!
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'cmumps_struc.h'
      INCLUDE 'fwd_struc.h'
      TYPE (CMUMPS_STRUC) MID  !MUMPS 
      TYPE (INV_INFO)     INV  !inversion information 
      TYPE (MESH_INFO)    MSH  !mesh, model parameters
      TYPE (SRC_INFO)     SRC  !source information
      TYPE (RECV_INFO)    RCV  !receiver information
      TYPE (FRQ_INFO)     FRQ  !frequency informatio 
      REAL*8, INTENT(IN) :: VFAST
      INTEGER*4, INTENT(IN) :: MASTER,MYID,MYNID, MGCOMM,MYSLV_COMM, MYHD_COMM, &
                               NPGROUPS,IPGROUP 
      LOGICAL*4 LPGNEW
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*8, ALLOCATABLE :: WAVE(:), UE(:), FMAT_DIST(:),  &
                                BDHESS(:) 
      REAL*4, ALLOCATABLE :: GRADL(:), HESSB(:), WORK(:)   
      LOGICAL*4, ALLOCATABLE :: LOBS(:,:) 
      COMPLEX*8 CZERO 
      REAL*8 PI180, PY
      REAL*4 EPSS, THRESHS, WMIN, WMAX
      INTEGER*4 IFREQ, IFREQL, ISG, JPGROUP, JSRC, ISRC,  &
                IBPHASE, NSGMAX, NOBS, MPIERR 
      REAL*4 COBJ4
      LOGICAL*4 LEXIST, LISDIR, LWORK 
      PARAMETER(CZERO = CMPLX(0.0,0.0))
      PARAMETER(IBPHASE = 0)   !we shouldn't have to apply geometric spreading correction
      PARAMETER(EPSS = 0.2)    !scales thresholding criteria in Huber norm, untested
      PARAMETER(THRESHS = 0.2) !scales thresholding criteria for L1/L2 norm, untested
      PARAMETER(PI180 = 0.017453292519943295D0) 

      complex*8, allocatable :: est(:,:,:,:), obs(:,:,:,:)
      complex*8 cphm2cm
      real*4 umag, vmag, wmag, uph, vph, wph, nmag, emag, zmag, nph, eph, zph
      integer*4 irec, i 
! 
!----------------------------------------------------------------------------------------!
! 
!.... nothing may be required, in which case, early return
      IERR = 0
      CALL MPI_BCAST(inv%LFUNC ,1,MPI_LOGICAL,MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(inv%LGRAD ,1,MPI_LOGICAL,MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(inv%LGNEWT,1,MPI_LOGICAL,MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(LPGNEW,1,MPI_LOGICAL,MASTER,MGCOMM,MPIERR)
      IF (.NOT.inv%LFUNC .AND. .NOT.inv%LGRAD .AND. &
          .NOT.LPGNEW.AND. .NOT.inv%LGNEWT) THEN
         IF (MYID == MASTER) &
         WRITE(*,*) 'funcgradh: No function, gradient, or Hessian required, returning...'
         RETURN
      ELSE
         IF (MYID == MASTER) THEN
            WRITE(*,*) ''
            IF (inv%LFUNC) WRITE(*,*) 'funcgradh: Will calculate function'
            IF (inv%LGRAD) WRITE(*,*) 'funcgradh: Will calculate gradient'
            IF (LPGNEW)    WRITE(*,*) 'funcgradh: Will calculate preconditioner'
            IF (inv%LGNEWT)WRITE(*,*) 'funcgradh: Will calculate and save Jacobians'
         ENDIF
      ENDIF
!
!.... set null the estimate wavefield for reduction 
      IF (MYNID == MASTER) THEN
         IF (MYID == MASTER) WRITE(*,*) 'funcgradh: Nulling estimates...'
         inv%EST(1:NDIM,1:frq%NFREQ,1:rcv%NREC,1:src%NSRC) = CMPLX(0.,0.)
      ENDIF
   
!
!.... set space for gradient
      IF (inv%LGRAD) THEN
         IF (MYNID == MASTER) THEN
            ALLOCATE(GRADL(inv%NA35)) 
            GRADL(1:inv%NA35) = 0.0 
            IF (MYID == MASTER) THEN
               inv%GRAD(1:inv%NA35) = 0.0
            ENDIF
         ELSE
            ALLOCATE(GRADL(1)) 
            GRADL(1) = 0.0 
         ENDIF
      ENDIF
!
!.... set space for hessian preconditioner
      IF (LPGNEW) THEN
         IF (MYNID == MASTER) THEN
            ALLOCATE(HESSB(inv%NHSIZE)) 
            HESSB(1:inv%NHSIZE) = 0.0 
         ELSE
            ALLOCATE(HESSB(1))
            HESSB(1) = 0.0 
         ENDIF
      ELSE
         ALLOCATE(HESSB(1))
         HESSB(1) = 0.0 
      ENDIF
!
!.... set space for local Jacobian
      IF (inv%LGNEWT) THEN
         IF (MYID == MASTER) THEN
            LEXIST = LISDIR('./scratch')
            IF (LEXIST) CALL SYSTEM(TRIM('rm -r ./scratch'))
            CALL SYSTEM(TRIM('mkdir ./scratch'))
         ENDIF
         IF (.NOT.inv%LFUNC) THEN
            WRITE(*,*) 'funcgradh: Error need residuals for a Gauss-Newton step!'
            IERR = 1 
            RETURN
         ENDIF
      ENDIF
      CALL MPI_BARRIER(MGCOMM,MPIERR)
!
!.... need to tabulate which rows of Jacobian are non-zero (no observations)
      IF (LPGNEW .OR. inv%LGNEWT) ALLOCATE(LOBS(NDIM,rcv%NREC))
!
!.... loop on frequencies in group
      ALLOCATE(UE(msh%NDOF)) 
      ALLOCATE(WAVE(msh%NDOF))  
      IF (inv%LGRAD .OR. LPGNEW .OR. inv%LGNEWT) ALLOCATE(FMAT_DIST(inv%NZ_FDIST))
      DO 1000 IFREQL=1,frq%NFREQ
!
!....... get frequency
         IFREQ = (IFREQL - 1)*NPGROUPS + IPGROUP + 1  
         IF ( (IFREQL - 1)*NPGROUPS + 1 > frq%NFREQ) GOTO 1005 !we are done
!
!....... set and send source group info
         IF (MYNID == MASTER) THEN
            IF (IFREQ <= frq%NFREQ) THEN
               CALL FILL_SRCPRM(IFREQ, VFAST,msh%AZTOL,msh%AOITOL, SRC) 
            ELSE
               src%NSG = 0 
            ENDIF
            CALL MPI_ALLREDUCE(src%NSG,NSGMAX,1,MPI_INTEGER, MPI_MAX,MYHD_COMM,MPIERR)
         ENDIF
         CALL MPI_BCAST(NSGMAX,1,MPI_INTEGER,MASTER, MYSLV_COMM,MPIERR)
         IF (IFREQ <= frq%NFREQ) CALL BCAST_SRC_PTRS(MYNID,MYSLV_COMM,MASTER, SRC)    
!
!....... loop on source groups
         DO 2000 ISG=1,NSGMAX
            DO 201 JPGROUP=0,NPGROUPS-1
               IF (IPGROUP == JPGROUP .AND. MYNID == MASTER .AND. &
                   IFREQ <= frq%NFREQ .AND. ISG <= src%NSG) THEN
                  IF (MYID == MASTER) WRITE(*,*) ''
                  WRITE(*,9409) IPGROUP+1,ISG,src%PYAVG(ISG),frq%FREQ(IFREQ)
 9409             FORMAT(' ---------------------------------------------------------',/, &
                    ' - funcgradh: Group:',I4,'                                 -'   ,/, &
                    ' -            Propagating source group ',I3,'               -'  ,/, &
                    ' -            With average slowness in y',E12.4,'     - '       ,/, &
                    ' -            At frequency',F12.5,' Hz                -'        ,/, &
                    ' ---------------------------------------------------------',/)
               ENDIF !end check on printing
               CALL MPI_BARRIER(MGCOMM,MPIERR) 
  201       CONTINUE
            IF (ISG > src%NSG) GOTO 2001 
            IF (IFREQ > frq%NFREQ) GOTO 2001 !this group is done with things to do
!
!.......... assemble and factor
            IF (MYID == MASTER) WRITE(*,*) 'funcgradh: Assembling impedance matrix...'
            CALL ASMBLE_DRIVER(msh%NZLOC,frq%CFTYPE(IFREQ),frq%FREQ(IFREQ),  &
                               src%PYAVG(ISG), MSH, mid%A_LOC,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'funcgradh: Error calling asmble_driver!'
               GOTO 500
            ENDIF
            IF (MYID == MASTER) WRITE(*,*) 'funcgradh: Factoring impedance matrix...'
            mid%JOB = 2
            CALL CMUMPS(MID) 
            IF (mid%INFO(1) < 0) THEN
               WRITE(*,*) 'funcgradh: An error occurred in the factorization!'
               IERR = 1
               GOTO 500
            ENDIF
!
!.......... loop on sources in source group
            DO 3000 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1)-1
               ISRC = src%ISRCPRM(JSRC) !extract original source ID 
!
!............. load the 1D solution
               IF (MYNID == MASTER) THEN
                  IF (MYID == MASTER) WRITE(*,9410) ISRC, src%AOI(ISRC), &
                                                    src%BAZN(ISRC), src%PYTAB(IFREQ,ISRC)
 9410             FORMAT(/,' funcgradh: Calculating ue field for source ',I3,/,          &
                           '            With angle of incidence ',F8.3,' degrees',/,     &
                           '            And adjusted azimuth ',F8.3,' degrees',/,        &
                           '            And y slowness ',G10.3,' s/m',/)
!
!................ this is where we convolve the source time function 
                  CALL LOAD_GRNS25(NDIM,msh%NDOF,msh%NNPE,NDIM, src%SRCTYP(ISRC),ISRC,   &
                                   frq%FREQ(IFREQ), src%SOURCE(IFREQ,ISRC),msh%IDOFSE,   &
                                   UE,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'funcgradh: Error calling load_grns25!'
                     GOTO 500 
                  ENDIF
                  IF (MYID == MASTER) &
                  WRITE(*,*) 'funcgradh: Calculating effective forces...'
               ENDIF !end check on mynid
!
!............. calculate the force distribution
               CALL MPI_BCAST(UE,msh%NDOF,MPI_COMPLEX, MASTER,MYSLV_COMM,MPIERR)
               CALL ASMB_PEFF(MASTER,MYSLV_COMM, msh%NDOF,msh%NDOFL,msh%NZLOC, msh%CNP,  &
                              msh%MYDOFS,msh%IRPTR_LOC,msh%JCPTR_LOC, mid%A_LOC,UE, WAVE)
               IF (MYNID == MASTER) CALL CCOPY(msh%NDOF,WAVE,1,mid%RHS,1) 
               mid%ICNTL(9) = 1 !always solve Ax = b
               mid%JOB = 3 !solve phase
               IF (MYID == MASTER) WRITE(*,*) 'funcgradh: Solving Ax=b...'
               CALL CMUMPS(MID)
               IF (mid%INFO(1) < 0) THEN 
                  WRITE(*,*) 'funcgradh: An error occurred in the solution!',MYID
                  IERR = 1  
                  RETURN
               ENDIF
!
!............. head processes extract response
               IF (MYNID == MASTER) THEN
                  CALL ADDBLK(msh%NDOF,msh%CNP,UE, mid%RHS) !add background field in
                  IF (rcv%NREC > 0) THEN
                     IF (MYID == MASTER) WRITE(*,*) 'funcgradh: Extracting responses...'
                     IF (inv%LUNWRAP) THEN
                        CALL EXRESP_MP_2(NDIM,rcv%NREC,msh%NDOF,                         &
                                      1,rcv%NREC,NDIM, frq%FREQ(IFREQ),                  &
                                      src%PYTAB(IFREQ,ISRC),rcv%YREC,mid%RHS, rcv%MRDOF, &
                                      rcv%RECV(1:NDIM,IFREQ,1:rcv%NREC),                 &
                                      inv%EST(1:NDIM,IFREQ,1:rcv%NREC,ISRC))
                     ELSE
                        CALL EXRESP25D(NDIM,rcv%NREC,msh%NDOF,                           &
                                     1,rcv%NREC,NDIM, frq%FREQ(IFREQ),                   &
                                     src%PYTAB(IFREQ,ISRC),rcv%YREC, MID%RHS, rcv%MRDOF, &
                                     rcv%RECV(1:NDIM,IFREQ,1:rcv%NREC),                  &
                                     inv%EST(1:NDIM,IFREQ,1:rcv%NREC,ISRC))
                     ENDIF
                  ENDIF
                     !CALL PLOT_ELCRESP_VTK(NGNOD,NDIM,msh%NEN,NDIM,msh%NNPG,msh%NDOF,   &   
                     !                    msh%NLXI,msh%NLETA,'good', msh%NELEM,3,ISRC,   &   
                     !                    3,frq%FREQ(IFREQ), msh%LM,msh%IENG,            &   
                     !                    msh%XLOCS,msh%ZLOCS, mid%RHS) 

               ENDIF !end check on myid
!
!............. gradient and hessian need a distributed F matrix 
               IF (inv%LGRAD .OR. LPGNEW .OR. inv%LGNEWT) THEN
                  IF (MYID == MASTER) &
                  WRITE(*,*) 'funcgradh: Calculating distributed F matrix...'
                  IF (MYNID == MASTER) CALL CCOPY(msh%NDOF,mid%RHS,1,WAVE,1) 
                  CALL MPI_BCAST(WAVE,msh%NDOF,MPI_COMPLEX, MASTER,MYSLV_COMM,MPIERR)
                  PY = src%PYTAB(IFREQ,ISRC)
                  CALL CFMAT_DIST(inv%NZ_FDIST,msh%NDOF, frq%FREQ(IFREQ),PY, &
                                  WAVE, MSH,INV, FMAT_DIST,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'funcgradh: Error calling cfmat_dist!'
                     GOTO 500
                  ENDIF
               ENDIF
!
!............. gradient computation: adj[dS/dm_1 u, ..., dS/dm_m u]
               IF (inv%LGRAD) THEN
                  IF (MYID == MASTER) THEN
                     WRITE(*,*) 'funcgradh: Calculating gradient...'
                     WRITE(*,*)
                  ENDIF
                  CALL ADJGRAD(MASTER,MYID,MYNID,MYSLV_COMM, inv%NZ_FDIST,msh%NDOF,   &
                               .TRUE.,.FALSE.,IFREQ,ISRC, msh%AZMOD,frq%FREQ(IFREQ),  &
                               EPSS,0.0, FMAT_DIST, RCV,MSH,INV,MID, GRADL,IERR) 
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'funcgradh: Error calling adjgrad!'
                     GOTO 500
                  ENDIF 
                  IF (MYID == MASTER) WRITE(*,*) 
               ENDIF
!
!............. calculate pre-conditioner and/or save local jacobians
               IF (LPGNEW .OR. inv%LGNEWT) THEN
!
!................ see which observations are active [rows of J are non-zero]
                  IF (MYNID == MASTER) THEN 
                     !write(myid+49,*), inv%obs(1:ndim,ifreq,1:rcv%nrec,isrc)
                     CALL LOBSLST(NDIM, rcv%NREC,NDIM,                            &
                                  inv%OBS(1:NDIM,IFREQ,1:rcv%NREC,ISRC), NOBS,LOBS)
                  ENDIF
                  CALL MPI_BCAST(NOBS,1,MPI_INTEGER, MASTER,MYSLV_COMM,MPIERR)
                  CALL BCAST_LOBS(MYSLV_COMM,MASTER, NDIM,NDIM,rcv%NREC, LOBS) 
                  LWORK = .TRUE.
                  IF (NOBS == 0) LWORK = .FALSE. !no observations here
                  IF (MYID == MASTER) WRITE(*,*) 'funcgradh: Calling pchadj...'
                  IF (LWORK) THEN
                     CALL PCHADJ(MYID,MYNID,MASTER,MYSLV_COMM, inv%NZ_FDIST,              &
                                 .FALSE.,.FALSE.,LPGNEW,IFREQ,ISRC,frq%FREQ(IFREQ), LOBS, &
                                 FMAT_DIST, RCV,INV,MID, HESSB,IERR) 
                     IF (IERR /= 0) THEN
                        WRITE(*,*) 'funcgradh: Error calling pchadj!'
                        GOTO 500
                     ENDIF
                  ELSE
                     WRITE(*,9501) IPGROUP+1,frq%FREQ(IFREQ),ISRC
 9501                FORMAT(' funcgradh: pachadj notice',/, &
                            '            Group:',I5,/, &
                            '            Skipping frequency:',F12.5,' and source',I5,/)
                  ENDIF
               ENDIF
 3000       CONTINUE !loop on sources in group 
 2001       CONTINUE !break ahead for loop on source group
 2000    CONTINUE !loop on sources in source group
         IF (src%NSG > 0 .AND. IFREQ <= frq%NFREQ) THEN
            IF (ASSOCIATED(src%ISGPTR))   DEALLOCATE(src%ISGPTR)
            IF (ASSOCIATED(src%ISRCPRM))  DEALLOCATE(src%ISRCPRM)
            IF (ASSOCIATED(src%PYAVG))    DEALLOCATE(src%PYAVG)
         ENDIF
         CALL MPI_BARRIER(MGCOMM,MPIERR)  
 1000 CONTINUE !loop on frequencies in inversion 
 1005 CONTINUE !we are out of frequencies
  500 CONTINUE !break ahead for errors
      IF (IERR /= 0) THEN 
         WRITE(*,*) 'funcgradh: An error was detected on process:',MYID
         RETURN
      ENDIF
!
!.... reduce estimates
      IF (MYNID == MASTER) THEN
!
!....... count up LGNEWT residuals
         !IF (MYID == MASTER .AND. inv%LGNEWT) THEN
         !   inv%NOBS = ICOBS(NDIM,frq%NFREQ,rcv%NREC,frq%NFREQ,rcv%NREC,src%NSRC,NDIM, &
         !                    inv%OBS)
         !   ALLOCATE(inv%RHS(inv%NOBS))
         !   inv%RHS(1:inv%NOBS) = 0.0 
         !   WRITE(*,*) 'funcgradh: Number of observations:',inv%NOBS
         !ENDIF
!
!....... calculate objective function
         IF (inv%LFUNC) THEN 
            IF (MYID == MASTER) WRITE(*,*) 'funcgradh: Reducing estimates...'
            CALL REDEST4(MYID,MASTER,MYHD_COMM,NPGROUPS,IPGROUP,           &   
                         NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC,1, inv%EST)
            IF (MYID == MASTER) THEN 
               ALLOCATE(EST(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
               ALLOCATE(OBS(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
               DO IFREQ=1,frq%NFREQ
                  DO IREC=1,rcv%NREC
                     DO ISRC=1,src%NSRC
                        IF (inv%LUNWRAP) THEN
                           UMAG = REAL(inv%EST(1,IFREQ,IREC,ISRC))
                           VMAG = REAL(inv%EST(2,IFREQ,IREC,ISRC))
                           WMAG = REAL(inv%EST(3,IFREQ,IREC,ISRC))
                           UPH  = IMAG(inv%EST(1,IFREQ,IREC,ISRC))
                           VPH  = IMAG(inv%EST(2,IFREQ,IREC,ISRC))
                           WPH  = IMAG(inv%EST(3,IFREQ,IREC,ISRC))
                           EST(1,IFREQ,IREC,ISRC) = CPHM2CM(UMAG,UPH)
                           EST(2,IFREQ,IREC,ISRC) = CPHM2CM(VMAG,VPH)
                           EST(3,IFREQ,IREC,ISRC) = CPHM2CM(WMAG,WPH)

                           NMAG = REAL(inv%OBS(1,IFREQ,IREC,ISRC))
                           EMAG = REAL(inv%OBS(2,IFREQ,IREC,ISRC))
                           ZMAG = REAL(inv%OBS(3,IFREQ,IREC,ISRC))
                           NPH  = IMAG(inv%OBS(1,IFREQ,IREC,ISRC))
                           EPH  = IMAG(inv%OBS(2,IFREQ,IREC,ISRC))
                           ZPH  = IMAG(inv%OBS(3,IFREQ,IREC,ISRC))
                           OBS(1,IFREQ,IREC,ISRC) = CPHM2CM(NMAG,NPH)
                           OBS(2,IFREQ,IREC,ISRC) = CPHM2CM(EMAG,EPH)
                           OBS(3,IFREQ,IREC,ISRC) = CPHM2CM(ZMAG,ZPH)
                        ELSE
                           DO I=1,NDIM
                              EST(I,IFREQ,IREC,ISRC) = inv%EST(I,IFREQ,IREC,ISRC)
                              OBS(I,IFREQ,IREC,ISRC) = inv%OBS(I,IFREQ,IREC,ISRC)
                           ENDDO
                        ENDIF 
                     ENDDO
                  ENDDO
               ENDDO  
               WRITE(*,*) 'funcgradh: Rotating estimates to (N,E,Z) frame...'
               CALL ROTEST(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC,   &
                           msh%AZMOD,EST) 
               WRITE(*,*) 'funcgradh: Calculating objective function...'
               inv%FOBJ = COBJ4(NDIM,frq%NFREQ,rcv%NREC, frq%NFREQ,rcv%NREC,src%NSRC, NDIM,  &
                                .FALSE.,inv%NORM,inv%IRESTP, EPSS,THRESHS,0.0,               &
                                frq%FREQ, inv%WGHTS,OBS,EST)
               !IF (inv%LGNEWT) THEN
               !   WRITE(*,*) 'funcgradh: Filling RHS least squares vector...'
               !   CALL FLLRHSGN(NDIM,frq%NFREQ,rcv%NREC, inv%NOBS, &
               !                 NDIM,frq%NFREQ,rcv%NREC,src%NSRC,  &
               !                 inv%NORM,inv%IRESTP, EPSS, inv%WGHTS, &
               !                 OBS,EST, inv%RHS)
               !ENDIF
               DEALLOCATE(EST)
               DEALLOCATE(OBS)
            ENDIF !end check on myid
         ENDIF !end check on function evaluation
!
!....... gradient objective fn reduction
         IF (inv%LGRAD) THEN 
            IF (MYID == MASTER) THEN
               ALLOCATE(WORK(inv%NA35))
            ELSE
               ALLOCATE(WORK(1))
            ENDIF
            IF (MYID == MASTER) WRITE(*,*) 'funcgradh: Reducing gradient...'
            CALL MPI_REDUCE(GRADL,WORK,inv%NA35,MPI_REAL,MPI_SUM, &
                            MASTER,MYHD_COMM,MPIERR)
            IF (MYID == MASTER) THEN
               CALL SCOPY(inv%NA35,WORK,1,inv%GRAD,1)
               WMIN = MINVAL(inv%WMASK)
               WMAX = MAXVAL(inv%WMASK)
               IF (WMIN /= 1.0 .AND. WMAX /= 1.0) THEN
                  WRITE(*,*) 'funcgradh: Multiplying gradient mask...'
                  CALL WTGRAD(inv%NA35,inv%WMASK, inv%GRAD)
               ENDIF
               IF (inv%NCASC > 0 .and. .not.inv%lpgrad) THEN
                  WRITE(*,*) 'funcgradh: Smoothing gradient...'
                  CALL AVGRAD(NGNOD,msh%NNPG, inv%NA35,msh%NNPG,     &
                              NGNOD,inv%NCON,inv%NVINV,inv%NCASC,    &
                              inv%MASKG,inv%MCONN,msh%IENG, inv%GRAD)  
               ENDIF
            ENDIF
            DEALLOCATE(WORK)
         ENDIF
!
!....... hessian calculation
         IF (LPGNEW) THEN
            IF (MYID == MASTER) THEN
               IF (inv%NVINV == 1) THEN
                  WRITE(*,*) 'funcgradh: Reducing diagonal Hessian...'
               ELSE
                  WRITE(*,*) 'funcgradh: Reducing block diagonal Hessian...'
               ENDIF
               ALLOCATE(BDHESS(inv%NHSIZE))
            ELSE
               ALLOCATE(BDHESS(1))
            ENDIF 
!           CALL MPI_REDUCE(CHESSB,BDHESS,inv%NHSIZE,MPI_COMPLEX,MPI_SUM, &
!                           MASTER,MYHD_COMM,MPIERR)
            CALL MPI_REDUCE(HESSB,BDHESS,inv%NHSIZE,MPI_REAL,MPI_SUM, &
                            MASTER,MYHD_COMM,MPIERR)
!
!.......... regularize the pre-conditioner 
            IF (MYID == MASTER) THEN
               WRITE(*,*) 'funcgradh: Regularizing gradient pre-conditioner...'
               CALL REGPRECON(inv%NHSIZE,inv%NA35, inv%NNPINV,inv%NVINV, inv%PTHRESH, &
                              BDHESS, inv%GRADPC) 
            ENDIF !end check on myid
            DEALLOCATE(BDHESS)
         ENDIF
      ENDIF !end check on mynid
      CALL MPI_BARRIER(MGCOMM,MPIERR) 
      IF (MYID == MASTER) WRITE(*,*) 'funcgradh: Freeing memory...'
!
!.... free memory
      IF (ALLOCATED(UE))        DEALLOCATE(UE) 
      IF (ALLOCATED(WAVE))      DEALLOCATE(WAVE) 
      IF (ALLOCATED(GRADL))     DEALLOCATE(GRADL)
      IF (ALLOCATED(HESSB))     DEALLOCATE(HESSB)
      IF (ALLOCATED(FMAT_DIST)) DEALLOCATE(FMAT_DIST)
      IF (ALLOCATED(LOBS))      DEALLOCATE(LOBS)
      RETURN
      END 
