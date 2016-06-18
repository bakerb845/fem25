!
!     This routine makes use of the Harwell Subroutine Library's L-BFGS 
!     inversion scheme to invert for model parameters.  Nominally, 
!     L-BFGS needs to be able to evaluate a function, gradient, and if desired
!     gradient preconditioner.  This is handled in funcgradh.f90 
!     - B. Baker March 2013
!
!.... includes 
      INCLUDE 'mpif.h'
      INCLUDE 'cmumps_struc.h'
      INCLUDE 'fwd_struc.h'
      TYPE (CMUMPS_STRUC) MID
      TYPE (SRC_INFO) SRC  
      TYPE (INV_INFO) INV 
      TYPE (RECV_INFO) RCV 
      TYPE (MESH_INFO) MSH 
      TYPE (FRQ_INFO) FRQ
      TYPE (MOD1D_INFO) M1D
!.... project name
      CHARACTER(80) PROJNM
!.... receiver information
      REAL*8, ALLOCATABLE :: XREC(:), ZREC(:) 
      LOGICAL*4, ALLOCATABLE :: LFSURF(:) 
      LOGICAL*4 LINREC
!.... funcgradh stuff
      LOGICAL*4 LPGNEW
      LOGICAL*4 LSRCEX,LRECEX 
!.... MPI stuff
      REAL*8 TSSIM, TESIM
      INTEGER*4 MASTER,MYID,MYNID, MYSLV_COMM,MYHD_COMM,IPGROUP,MPIERR, &
                NPGROUPS, NPROCS, NPARTS, MYDEST   
      PARAMETER(MASTER = 0)
!.... miscellaneous stuff
      CHARACTER(80) ESTFL, TMPDIR
      CHARACTER(5) CK, CBLOCK
      INTEGER*4 NSRC_HD, NFREQ_HD, NREC_HD, NDIM_HD, NELEME, NABS, IOMINV, &
                IBLOCK, NBLOCKS, K, NFN35, NH35, ITYPE, ICALPHA, NBLOCK,   &
                ISTOP_PT, IERR  
      LOGICAL*4 LCGRNS, LFILES, LKBDRY, LNSRF
      PARAMETER(LFILES = .FALSE.) 
!.... va35s parameters
      REAL*4, ALLOCATABLE :: DIAG(:)    !contains diagonal H_k^0 values
      REAL*4, ALLOCATABLE :: GRAD(:)    !gradient
      REAL*4, ALLOCATABLE :: XMOD(:)    !model to invert for
      REAL*4, ALLOCATABLE :: HESS(:)    !for plotting terms in the Hessian
      REAL*4, ALLOCATABLE :: SEARCH(:)  !search direction from L-BFGS
      REAL*8, ALLOCATABLE :: X8(:)      !model 
      REAL*8, ALLOCATABLE :: X0(:)      !old model
      REAL*8, ALLOCATABLE :: G8(:)      !gradient
      REAL*8, ALLOCATABLE :: G0(:)      !old gradient direction
      REAL*8, ALLOCATABLE :: D8(:)      !diagonal inverse hessian
      REAL*8, ALLOCATABLE :: P8(:)      !search direction 
      REAL*8, ALLOCATABLE :: Y(:)       !past differenced gradient directions
      REAL*8, ALLOCATABLE :: S(:)       !past search directions
      REAL*8, ALLOCATABLE :: WA(:)      !workspace for line search
      INTEGER*4, ALLOCATABLE :: IPIV(:) !holds pivots for double preicison hessian diag
      REAL*8 EPS8, GTOL8, XTOL8, FTOL8, STPMIN, STPMAX, STP, DS0, F8 
      INTEGER*4 INFO, MLBFGS, NWORK35, IFLAG  
!.... functions
      INTEGER*4 ICELEML
!
!----------------------------------------------------------------------------------------!
!
!.... initialize MPI
      CALL MPI_INIT(MPIERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, MPIERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, MPIERR)
      TSSIM = MPI_WTIME()
!
!.... these are specific to each program
      inv%LMIGR  = .FALSE. !this is not a migration program
      inv%LGNEWT = .FALSE. !do not perform a Gauss-Newton step
!
!.... head node reads model information
      IERR = 0
      IF (MYID == MASTER) THEN
         PROJNM(1:80) = ' ' 
         WRITE(*,8000) 
 8000    FORMAT(' ---------------------------------------------------------',/, &
                ' -  xlbfgs_hsl: A massively parallel 2.5D unstructured   -',/, &
                ' -              elastic inversion algorithm utilizing    -',/, & 
                ' -              the HSL L-BFGS inversion algorithm       -',/, &
                ' ---------------------------------------------------------',/)
         WRITE(*,*) 'xlbfgs_hsl: Enter project name:' 
         READ(*,'(A)') PROJNM
         PROJNM = ADJUSTL(PROJNM)
         WRITE(*,*) 'xlbfgs_hsl: Enter number of process groups:' 
         READ *, NPGROUPS

         mlbfgs = 3  !memory for L-BFGS correction 
         msh%freq0 = 20.d0
         msh%aztol = 10.d0
         msh%aoitol = 5.d0
         inv%cinvtype = 'PP' !one at a time is enough 
         inv%pthresh =  5.  !percent of singular values below which are inflated
         eps8 = 1.d-5
         gtol8 = 9.d-1  
         xtol8 = EPSILON(1.D0)*100.D0
         ftol8 = 1.d-3 
         stpmin = 1.d-20
         stpmax = 1.d+20
         maxfev = 6 
         lcgrns = .false. !true calculate new greens fns, false read old ones
         inv%maxit = 5  !max number of iterations
         inv%norm = 2  !1 -> L1 norm, 2 -> L2 norm (default)
         inv%lunwrap = .TRUE. !try working with unwrapped phase in data
         inv%irestp = 1  !1 -> phase only, 2 -> amplitude only, 3 -> both (default) 
         inv%ibphase = 0 !no geometric spreading correction in backpropagation 
         inv%lpgrad = .TRUE. !calculate a pseudo hessian for gradient pre-conditioning
         inv%ncasc = 0
         inv%imodsrc = 3 !1 -> phase only, 2 -> amplitude only, 3 -> both (default)
         inv%imodrec = 3 !1 -> phase only, 2 -> amplitude only, 3 -> both (default)
         linrec = .true. !include receivers in inversion points?
         lkbdry = .true. !include nodes on interior/bielak bdry
         ds0  = 55.d0    !the gradient (or preconditioned gradient) are poorly scaled, 
                          !this is the initial rescale factor (m/s).   
         istop_pt = 0 !stop point; n/a here
!....... initial warnings on NPGROUPS
         IF (MOD(NPROCS,NPGROUPS) /= 0) THEN 
            WRITE(*,*) 'xlbfgs_hsl: Error I cant divide nprocs by npgroups evenly!'
            IERR = 1 
            GOTO 500 
         ENDIF
         NPARTS = NPROCS/NPGROUPS
!
!....... check the inversion scheme
         IF (inv%CINVTYPE == 'PP' .OR. inv%CINVTYPE == 'pp' .OR.      &
             inv%CINVTYPE == 'SS' .OR. inv%CINVTYPE == 'ss' .OR.      &
             inv%CINVTYPE == 'PS' .OR. inv%CINVTYPE == 'ps')  THEN
            inv%NVINV = 2 
            IF (inv%CINVTYPE == 'PP' .OR. inv%CINVTYPE == 'pp' .OR.  &
                inv%CINVTYPE == 'SS' .OR. inv%CINVTYPE == 'ss') THEN
               inv%NVINV = 1
            ENDIF
         ELSE
            WRITE(*,*) 'xlbfgs_hsl: I do not know what to invert for:',inv%CINVTYPE
            IERR = 1 
            GOTO 500 
         ENDIF 
!
!....... read the mesh
         WRITE(*,*) 'xlbfgs_hsl: Reading mesh...'
         CALL RDMESHBK(PROJNM, NELEME,NABS, MSH, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xlbfgs_hsl: Error reading mesh!'
            IERR = 1 
            GOTO 500 
         ENDIF
         ALLOCATE(msh%XIPTS (msh%NLXI))
         ALLOCATE(msh%ETAPTS(msh%NLETA)) 
         CALL GENPTS(msh%IITYPE, msh%NLXI,msh%NLETA, msh%XIPTS,msh%ETAPTS) 
!
!....... initial check on observation file 
         WRITE(*,*) 'xlbfgs_hsl: Checking headers on observation file...'
         CALL RDTOBS_HD(PROJNM, NDIM_HD,NFREQ_HD,NREC_HD,NSRC_HD, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xlbfgs_hsl: There was an error with your observation file'
            GOTO 500
         ENDIF 
         IF (NDIM_HD /= NDIM) THEN
            WRITE(*,*) 'xlbfgs_hsl: Error component number mismatch!',NDIM_HD,NDIM
            IERR = 1
            GOTO 500
         ENDIF
         IF (NFREQ_HD <= 0) THEN
            WRITE(*,*) 'xlbfgs_hsl: Error observations have 0 frequencies!',NFREQ_HD
            IERR = 1
            GOTO 500
         ENDIF
!
!....... get the number of blocks
         CALL RDFREQI_BLHD(PROJNM,NBLOCKS,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xlbfgs_hsl: Error cannot locate the frequency block file'
            GOTO 500
         ENDIF
!
!....... read the source list
         CALL RDSRC_EQ(PROJNM, msh%XLATMIN,msh%XLONMIN,msh%XLATMAX,msh%XLONMAX, &
                       msh%AZMOD, SRC,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xlbfgs_hsl: Error calling rdsrc_eq!'
            GOTO 500 
         ENDIF
         IF (NSRC_HD /= src%NSRC) THEN
            WRITE(*,*) 'xlbfgs_hsl: Error source number mismatch!',NSRC_HD,src%NSRC
            IERR = 1
            GOTO 500
         ENDIF
         inv%LBODY = .FALSE.
         inv%LSURF = .FALSE.
         IF (src%SRCTYP(1)(1:1) == 'P') THEN
            inv%LBODY = .TRUE.
         ELSE
            inv%LSURF = .TRUE.
         ENDIF
!
!....... check for a receiver file
         CALL RDREC_HD(PROJNM, rcv%NREC,IERR)
         IF (IERR /= 0) THEN
            IF (IERR > 0) THEN
               WRITE(*,*) 'xlbfgs_hsl: Error reading recv header'
               GOTO 500 
            ELSE
               WRITE(*,*) 'xlbfgs_hsl: There is no receiver file!  Cannot make estimates!'
               IERR = 1 
               GOTO 500
            ENDIF 
         ENDIF 
         IF (rcv%NREC <= 0) THEN
            WRITE(*,*) 'xlbfgs_hsl: There are no receivers!  Aborting!'
            IERR = 1
            GOTO 500
         ENDIF
         IF (NREC_HD /= rcv%NREC) THEN
            WRITE(*,*) 'xlbfgs_hsl: Error recevier number mismatch!',NREC_HD,rcv%NREC  
            IERR = 1
            GOTO 500
         ENDIF
         ALLOCATE(  LFSURF(rcv%NREC))
         ALLOCATE(    XREC(rcv%NREC))
         ALLOCATE(rcv%YREC(rcv%NREC))
         ALLOCATE(    ZREC(rcv%NREC))
         CALL RDREC(PROJNM, rcv%NREC, LFSURF,XREC,rcv%YREC,ZREC, IERR)
         ALLOCATE(inv%DX(rcv%NREC-1))
         inv%DX(1:rcv%NREC-1) = 0.D0 
         IF (inv%LCOVD_SRF .OR. inv%LCOVD_BDY) THEN 
            CALL REC_SPACE1D(rcv%NREC,XREC, inv%DX,IERR)
            IF (IERR /= 0) THEN 
               WRITE(*,*) 'xlbfgs_hsl: Error detected in rec_space1d'
               IF (inv%LCOVD_SRF) WRITE(*,*) 'xlbfgs_hsl: Setting lcovd_srf to false'
               IF (inv%LCOVD_BDY) WRITE(*,*) 'xlbfgs_hsl: Setting lcovd_bdy to false'
            ENDIF
         ENDIF
 
         WRITE(*,*) 'xlbfgs_hsl: Splitting process groups...'
      ENDIF
      CALL MKMPI_GROUPS(MASTER,MYID,NPGROUPS,NPROCS,MPI_COMM_WORLD, &
                        IPGROUP,MYNID,MYHD_COMM,MYSLV_COMM)
!----------------------------------------------------------------------------------------!
!     MUMPS initialization phase and graph reordering utilities.  The mesh should not    !
!     have to change.  If it does the user should re-mesh then re-run the inversion      !
!     software.                                                                          !
!----------------------------------------------------------------------------------------!
      IF (MYID == MASTER) WRITE(*,*) 'xlbfgs_hsl: Initializing MUMPS...'
      mid%COMM = MYSLV_COMM !set communicator
      mid%SYM = 0 !unsymmetric
      mid%JOB =-1 !initialize
      mid%PAR = 1 !host host working 
      CALL CMUMPS(MID)
      mid%ICNTL(3) = 0 !suppress output
      mid%ICNTL(4) = 1 !only error messages
!
!.... generate a graph
      IF (MYID == MASTER) THEN
         WRITE(*,*) 'xlbfgs_hsl: Generating graph...'
         CALL GEN_GRAPH25(.TRUE.,NPARTS, LFSURF,XREC,ZREC, MSH,RCV,IERR) 
         IF (inv%LUNWRAP) THEN 
            WRITE(*,*) 'xsrcrec25: Generating free surface information...'
            CALL GEN_GIDOFFS(MSH,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xsrcrec25: Error calling gen_gidoffs!'
               GOTO 500
            ENDIF
         ENDIF
         mid%N  = msh%NDOF
         mid%NZ = msh%NZERO
         CALL PLOT_MESHP_VTK(PROJNM,NGNOD,NDIM,msh%NEN, msh%NNPG,msh%NELEM,msh%NDOF,     &
                             NDIM,msh%NEN, msh%PART,msh%IENG,msh%LM, msh%XLOCS,msh%ZLOCS)
         WRITE(*,*)
         WRITE(*,9408) msh%NORD,msh%NELEM,NABS,NELEME, &
                       msh%NNPG,msh%NNPE,msh%NDOF,msh%NZERO  
 9408    FORMAT(' xlbfgs_hsl: Polynomial order:'                    ,I4 ,/,        &   
                '             Number of elements in mesh:'          ,I10,/,         &   
                '             Number of absorbing elements:'        ,I8 ,/,         &   
                '             Number of Bielak elements:'           ,I8 ,/,         &   
                '             Number of anchor nodes in mesh:'      ,I10,/,         &   
                '             Number of nodes in Bielak boundary:'  ,I10,/,         &   
                '             Number of degrees of freedom:'        ,I14,/,         &   
                '             Number of non-zeros in global matrix:',I16,/)
         IF (ALLOCATED(LFSURF)) DEALLOCATE(LFSURF) 
         IF (ALLOCATED(ZREC)) DEALLOCATE(ZREC)
!
!....... create the pointers associated w/ the inverse problem
         WRITE(*,*) 'xlbfgs_hsl: Generating mask for inverse problem...'
         CALL GENGMASK(PROJNM, NDIM,msh%NEN,NGNOD, msh%NDOF,msh%NNPG,msh%NELEM,     &
                       msh%NLXI,msh%NLETA,LINREC,LKBDRY,                            &
                       msh%CDOMAIN,msh%CNNPG,msh%PART,                              &
                       msh%LM,msh%IENG, RCV,INV, IERR) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xlbfgs_hsl: Error generating gradient mask!'
            GOTO 500
         ENDIF
         CALL ELEM_WTS(MSH,INV,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xlbfgs_hsl: Error calling ELEM_WTS'
            GOTO 500
         ENDIF
         inv%NA35 = inv%NNPINV*inv%NVINV
         inv%MCJLOC = rcv%NREC*NDIM
         ALLOCATE(DIAG(inv%NA35))
         DIAG(1:inv%NA35) = 1.0
         ALLOCATE(inv%GRAD(inv%NA35))
         ALLOCATE(XMOD(inv%NA35)) 
         ALLOCATE(SEARCH(inv%NA35)) !search direction from L-BFGS 
         ALLOCATE(GRAD(inv%NA35)) !gradient for L-BFGS
         GRAD(1:inv%NA35) = 0.0
         !NWORK35 = inv%NA35*(2*MLBFGS + 1) + 2*MLBFGS
         NWORK35 = MLBFGS*inv%NA35 
         ALLOCATE(Y(NWORK35))   !saved differenced gradient directions
         ALLOCATE(S(NWORK35))   !saved search directions
         ALLOCATE(WA(inv%NA35)) !workspace for line search
         ALLOCATE(X0(inv%NA35)) !old model 
         ALLOCATE(G0(inv%NA35)) !old gradient direction
         ALLOCATE(X8(inv%NA35)) !model
         ALLOCATE(G8(inv%NA35)) !gradient
         ALLOCATE(P8(inv%NA35)) !search direction
!
!....... also set up the gradient preconditioner 
         inv%NHSIZE = 0 
         IF (inv%LPGRAD) THEN
            NBLOCK = inv%NVINV**2            !each submatrix is nvinv x nvinv
            inv%NHSIZE = inv%NNPINV*NBLOCK   !total number of non-zeros
            ALLOCATE(inv%GRADPC(inv%NHSIZE)) !gradient pre-conditioner
            ALLOCATE(inv%IPIVH(inv%NA35))    !pivots in LU factorization (could use 
                                             !cholesky, but small terms worry me) 
            NHSIZE = inv%NHSIZE 
         ELSE
            NHSIZE = inv%NA35 
         ENDIF
         ALLOCATE(D8(NHSIZE))
         ALLOCATE(IPIV(inv%NA35)) 
      ENDIF
!----------------------------------------------------------------------------------------!
!                          Broadcast Inverse Problem Parameters                          !
!----------------------------------------------------------------------------------------!
      IF (MYID == MASTER) WRITE(*,*) 'xlbfgs_hsl: Broadcasting 2D model...' 
      CALL MPI_BCAST(NPARTS ,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
      
      CALL BCAST_MESH_INFO(MYID,MPI_COMM_WORLD,MASTER, &
                           LNSRF, PROJNM,TMPDIR, MSH)
      IF (MYNID == MASTER) THEN
         IF (MYID == MASTER) WRITE(*,*) 'xlbfgs_hsl: Broadcasting receiver locations...'
         CALL BCAST_RCV_INFO(MYID,MYHD_COMM,MASTER, RCV)
      ENDIF
      IF (MYID == MASTER) WRITE(*,*) 'xlbfgs_hsl: Broadcasting inversion parameters...'
      CALL MPI_BCAST(rcv%NREC,1,MPI_INTEGER, MASTER,MYSLV_COMM,MPIERR) !for jacobian
      CALL BCAST_INV_PARMS(MYID,MPI_COMM_WORLD,MASTER, msh%NNPG,msh%NELEM,NBLOCKS, &
                           rcv%NREC, ISTOP_PT, INV) 
      IF (MYNID == MASTER .AND. inv%LUNWRAP) THEN
         IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Broadcasting free surface DOFS...'
         CALL BCAST_FS_INFO(MYID,MYHD_COMM,MASTER, MSH)
      ENDIF
!
!.... figure out the local gradient graph
      IF (MYID == MASTER) WRITE(*,*) 'xlbfgs_hsl: Generating gradient graph...'
      CALL GRADPTRS(MYNID,MASTER,MYSLV_COMM, NDIM,msh%NEN, msh%NEN,msh%NNPG, msh%LM, &
                    INV,IERR)
      IF (IERR /= 0) THEN
         IF (MYNID == MASTER) WRITE(*,*) 'xlbfgs_hsl: An error occurred in gradptrs!'
         GOTO 500
      ENDIF 
!----------------------------------------------------------------------------------------!
!                 Graph distribution.  Note we do not reorder the matrix as MUMPS        !
!                 calls METIS_NodeWND which seems to be a bit more clever than           !
!                 METIS_NodeND when the matrix is distributed                            !
!----------------------------------------------------------------------------------------!
!
!.... have masters save matrix sizes 
      IF (MYNID == MASTER) THEN
         mid%N = msh%NDOF
         mid%NZ = msh%NZERO
         mid%ICNTL(18) = 3 !matrix distributed by user 
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
      IF (MYID == MASTER) WRITE(*,*) 'xlbfgs_hsl: Generating local graphs...'
      CALL ICNDZ_LOC(MYNID,1,msh%NDOF, msh%PART,msh%IRPTR, msh%NDOFL,msh%NZLOC)
      ALLOCATE(msh%MYDOFS(msh%NDOFL))
      ALLOCATE(msh%IRPTR_LOC(msh%NDOFL+1))
      ALLOCATE(msh%JCPTR_LOC(msh%NZLOC))
      CALL PART2CRSL(msh%NZERO,msh%NZLOC, msh%NDOF,msh%NDOFL, MYNID,1,             &   
                     msh%IRPTR,msh%JCPTR,msh%PART, msh%MYDOFS,msh%IRPTR_LOC,msh%JCPTR_LOC)
      msh%NELEML = ICELEML(MYNID, NDIM,msh%NEN,msh%NDOF, msh%NELEM, 1,             &   
                           msh%NEN,NDIM, msh%PART,msh%LM)
      ALLOCATE(msh%MYELEM(msh%NELEML))
      CALL GENELEML(MYNID, NDIM,msh%NEN,msh%NDOF, msh%NELEM,msh%NELEML, 1,         &   
                    msh%NEN,NDIM, msh%PART,msh%LM, msh%MYELEM)
      IF (ASSOCIATED(msh%IRPTR)) DEALLOCATE(msh%IRPTR)
      IF (ASSOCIATED(msh%JCPTR)) DEALLOCATE(msh%JCPTR)
      IF (ASSOCIATED(msh%PART))  DEALLOCATE(msh%PART)
      MID%NZ_LOC = msh%NZLOC
      ALLOCATE(mid%IRN_LOC(mid%NZ_LOC))
      ALLOCATE(mid%JCN_LOC(mid%NZ_LOC)) 
      CALL CRS2COOLOC(msh%NZLOC,msh%NDOFL, msh%MYDOFS,msh%IRPTR_LOC,msh%JCPTR_LOC,       &   
                      mid%IRN_LOC,mid%JCN_LOC)
      mid%JOB = 1 
      CALL CMUMPS(MID)
      ALLOCATE(mid%A_LOC(mid%NZ_LOC))
      IF (MYNID == MASTER) THEN
         mid%LRHS = mid%N
         ALLOCATE(mid%RHS(mid%N)) 
      ENDIF
      F0 = 0.0 
!----------------------------------------------------------------------------------------!
!     This is the loop on frequency blocks.  To help linearize the inversion we adopt    !
!     a multiscale strategy.  The option is certainly here, but I don't recommend        !
!     using it and would instead recommend shell scripts because once an artifcact is    !
!     introduced it is not coming out.                                                   ! 
!----------------------------------------------------------------------------------------!
      DO 1000 IBLOCK=1,NBLOCKS
!
!....... read the frequencies, source weights, observations in block
         IF (MYID == MASTER) THEN
            WRITE(*,*) 
            WRITE(*,*) 'xlbfgs_hsl: Filling 1D models...'
            CALL FILL1D(MSH,M1D)
         ENDIF
         CALL MPI_BCAST(msh%XMOD0  ,1,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,MPIERR) 
         CALL MPI_BCAST(msh%XMOD1  ,1,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,MPIERR) 
         CALL MPI_BCAST(msh%XBLKL  ,1,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(msh%XBLKR  ,1,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,MPIERR)
         IF (MYID == MASTER) THEN
            WRITE(*,*) 'xlbfgs_hsl: Reading frequency block...'
            CALL RDFREQI(PROJNM,IBLOCK, FRQ,IERR) 
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xlbfgs_hsl: Cannot read inversion frequency block!'
               GOTO 500
            ENDIF   
            WRITE(*,*) 'xlbfgs_hsl: Reading observation file...'
            ALLOCATE(inv%  EST(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
            ALLOCATE(inv%  OBS(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
            ALLOCATE(inv%WGHTS(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
            CALL RDTOBS25(PROJNM, NDIM,frq%NFREQ,rcv%NREC,frq%NFREQ, NDIM,rcv%NREC,    &
                          src%NSRC, inv%LUNWRAP, frq%FREQ, inv%WGHTS,inv%OBS, IERR)
            IF (IERR /= 0) THEN
               DEALLOCATE(inv%OBS)
               DEALLOCATE(inv%WGHTS) 
               WRITE(*,*) 'xlbfgs_hsl: Error reading observation file1'
               GOTO 500
            ENDIF
            WRITE(*,*) 'xlbfgs_hsl: Checking for srcinv file...'
            ALLOCATE(src%SOURCE(frq%NFREQ,src%NSRC))
            CALL RDSRCINV(frq%NFREQ,PROJNM, frq%NFREQ,src%NSRC, frq%FREQ, &
                          LSRCEX,src%SOURCE)
            IF (.NOT.LSRCEX) THEN
               WRITE(*,*) 'xlbfgs_hsl: No srcinv file detected'
               IF (inv%IMODSRC == 0) THEN
                  WRITE(*,*) 'xlbfgs_hsl: Overriding and estimating source' 
                  WRITE(*,*) 'xlbfgs_hsl: Overriding to invert for source phase and amplitude' 
                  inv%IMODSRC = 3
               ENDIF
            ELSE
               WRITE(*,*) 'xlbfgs_hsl: Source inversion file read'
            ENDIF
            ALLOCATE(rcv%RECV(NDIM,frq%NFREQ,rcv%NREC)) 
            CALL RDREC_RESP(PROJNM, NDIM,frq%NFREQ, NDIM,frq%NFREQ,rcv%NREC, msh%AZMOD, &
                            frq%FREQ, LRECEX,rcv%RECV)
            IF (.NOT.LRECEX) THEN
               WRITE(*,*) 'xlbfgs_hsl: Receiver response functions are set to unity'
               IF (inv%IMODREC == 0) THEN
                  WRITE(*,*) 'xlbfgs_hsl: Overriding and estimating receiver response'
                  WRITE(*,*) 'xlbfgs_hsl: Overriding to invert for receiver phase and amplitude'
                  inv%IMODREC = 3
               ENDIF
               IERR = 0
            ENDIF
            WRITE(*,8100) IBLOCK 
            DO 1101 IOMINV=1,frq%NFREQ 
               WRITE(*,8101) IOMINV,frq%FREQ(IOMINV) 
 1101       CONTINUE !loop on frequency blocks 
            WRITE(*,8102) 
 8100       FORMAT(/,' -------------------------------------------------',/, &
                     ' -   Frequencies in Block:',I3,'                    -',/, &
                     ' -                                               -')
 8101       FORMAT(  ' -   Frequency Number:',I3,1X,F12.5,' (Hz)      -')
 8102       FORMAT(  ' -                                               -',/, &
                     ' -------------------------------------------------',/)
         ENDIF
!
!....... broadcast frequency block and observations
         IF (MYID == MASTER) WRITE(*,*) 'xlbfgs_hsl: Broadcasting 1D models...'
         CALL BCAST_MOD1D_INFO(MYID,MPI_COMM_WORLD,MASTER, MSH,M1D)
         IF (MYID == MASTER) WRITE(*,*) 'xlbfgs_hsl: Broadcasting source details...'
         CALL BCAST_SRC_INFO(MYID,MPI_COMM_WORLD,MASTER, SRC)
         IF (MYID == MASTER) WRITE(*,*) 'xlbfgs_hsl: Broadcasting frequency information...'
         CALL BCAST_FRQ_INFO(MYID,MPI_COMM_WORLD,MASTER, FRQ)
         IF (MYNID == MASTER) THEN
            IF (MYID == MASTER) &
            WRITE(*,*) 'xlbfgs_hsl: Broadcasting observations...'
            CALL BCAST_OBS_INFO(MYID,MYHD_COMM,MASTER, frq%NFREQ,rcv%NREC,src%NSRC, &
                                INV) 
            IF (MYID /= MASTER) ALLOCATE(inv%EST(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
            WRITE(*,*) 'xlbfgs_hsl: Broadcasting receiver response functions...'
            CALL BCAST_RCV_RESP(MYID,MYHD_COMM,MASTER, .TRUE.,frq%NFREQ,RCV)
            IF (MYID == MASTER) &
            WRITE(*,*) 'xlbfgs_hsl: Broadcasting source time function...'
            CALL BCAST_SRC_STF(MYID,MYHD_COMM,MASTER,  .TRUE.,frq%NFREQ, SRC)  
         ENDIF
!----------------------------------------------------------------------------------------!
!                            Calculate 1D models solutions                               !
!----------------------------------------------------------------------------------------!
         IF (MYID == MASTER) WRITE(*,*) 'xlbfgs_hsl: Calculating 1D solutions...'
         CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
         ALLOCATE(src%PYTAB(frq%NFREQ,src%NSRC))
         src%PYTAB(1:frq%NFREQ,1:src%NSRC) = 0.D0
         CALL GEN_GRNS(MYID,MASTER,MPI_COMM_WORLD, TMPDIR,LCGRNS,LNSRF, &
                       SRC,M1D,FRQ,MSH, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xlbfgs_hsl: Error in gen_grns!'
            GOTO 500 
         ENDIF
!
!....... free source info
         IF (IBLOCK == NBLOCKS) THEN
            IF (MYNID /= MASTER) THEN
               IF (ASSOCIATED(src%SRCTYP))  DEALLOCATE(src%SRCTYP)
               IF (ASSOCIATED(src%BAZN))    DEALLOCATE(src%BAZN)
               IF (ASSOCIATED(src%AOI))     DEALLOCATE(src%AOI)
            ENDIF
            IF (ASSOCIATED(src%SLAT))    DEALLOCATE(src%SLAT)
            IF (ASSOCIATED(src%SLON))    DEALLOCATE(src%SLON)
            IF (ASSOCIATED(src%SDEP))    DEALLOCATE(src%SDEP)
            IF (ASSOCIATED(src%STRIKE))  DEALLOCATE(src%STRIKE)
            IF (ASSOCIATED(src%DIP))     DEALLOCATE(src%DIP)
            IF (ASSOCIATED(src%RAKE))    DEALLOCATE(src%RAKE)
            IF (ASSOCIATED(src%SMAG))    DEALLOCATE(src%SMAG)
            IF (ASSOCIATED(src%MODE))    DEALLOCATE(src%MODE)
         ENDIF
!
!....... this is the initialization step 
         inv%LFUNC = .TRUE. !Most certainly need a function evaluation
         inv%LGRAD = .TRUE. !VA35 takes a gradient first go around
         IF (MYID == MASTER) THEN
            IF (inv%LPGRAD) THEN
               WRITE(*,*) 'xlbfgs_hsl: Initializing objective function, gradient,'
               WRITE(*,*) '            and preconditioner...'
               LPGNEW = .TRUE.
            ELSE
               WRITE(*,*) 'xlbfgs_hsl: Initializing objective function and gradient...'
               LPGNEW = .FALSE.
            ENDIF
         ENDIF
         CALL FUNCGRADH25(MASTER,MYID,MYNID,IPGROUP,NPGROUPS,               &
                          MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM, LPGNEW,      &
                          m1d%VFAST, INV,MSH,SRC,RCV,FRQ, MID, IERR) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xlbfgs_hsl: An error occurred while calling funcgradh!'
            GOTO 500
         ENDIF
         IF (MYID == MASTER) THEN
            DIAG(1:inv%NA35) = 1.0
            CALL GENXMOD(msh%NNPG,inv%NA35,msh%NNPG, inv%CINVTYPE,   &
                         inv%MASKG,msh%DENS,msh%ECOEFF, XMOD,IERR)
            IF (inv%LGRAD) THEN
               CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,   &
                                    inv%CINVTYPE, 1,IBLOCK,0,0,            &
                                    inv%MASKG,msh%IENG,msh%XLOCS,msh%ZLOCS, inv%GRAD)
            ENDIF
            IF (LPGNEW) THEN
               WRITE(*,*) 'xlbfgs_hsl: Factoring gradient preconditioner...'
               D8(1:inv%NHSIZE) = DBLE(inv%GRADPC(1:inv%NHSIZE)) 
               CALL DINVHESS(NHSIZE,inv%NA35, inv%NNPINV,inv%NVINV, D8,IPIV,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xlbfgs_hsl: An error occurred in INVHESS'
                  GOTO 500
               ENDIF
               !IF (inv%NVINV == 1) THEN
               !   WRITE(*,*) 'xlbfgs_hsl: Averaging block diagonal...'
               !   print *, NGNOD,msh%NNPG, inv%NA35,msh%NNPG, NGNOD,inv%NCON
               !   print *, inv%NVINV,inv%NCASC
               !   CALL AVGRAD8(NGNOD,msh%NNPG, inv%NA35,msh%NNPG, NGNOD,inv%NCON, &
               !                inv%NVINV,inv%NCASC, &
               !                inv%MASKG,inv%MCONN,msh%IENG, D8) 
               !ENDIF

               ALLOCATE(HESS(inv%NA35))
               IF (inv%CINVTYPE == 'PP' .OR. inv%CINVTYPE == 'pp' .OR. &
                   inv%CINVTYPE == 'SS' .OR. inv%CINVTYPE == 'ss') THEN
                  WRITE(*,*) 'xlbfgs_hsl: Copying diagonal...'
                  CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,          &
                                       inv%CINVTYPE, 2,IBLOCK,0,0,            &   
                                       inv%MASKG,msh%IENG,msh%XLOCS,msh%ZLOCS, inv%GRADPC) 
               ELSEIF (inv%CINVTYPE == 'PS' .OR. inv%CINVTYPE == 'ps') THEN 
                  K = 0 
                  WRITE(*,*) 'xlbfgs_hsl: Ben you should redo this, refactoring...'
                  CALL INVHESS(inv%NHSIZE,inv%NA35, inv%NNPINV,inv%NVINV,  &
                               inv%GRADPC,inv%IPIVH,IERR)
                  DO ITYPE=1,2
                     CALL GETBDIAG_HESS(inv%NHSIZE,inv%NNPINV,inv%NVINV,inv%NA35, ITYPE, &
                                        inv%IPIVH,inv%GRADPC, HESS)   
                     I1 = (ITYPE -1)*inv%NA35 + 1
                     I2 = I1 + inv%NA35 - 1
                     WRITE(*,*) 'xlbfgs_hsl: Copying diagonal...'
                     DIAG(I1:I2) = 1./HESS(1:inv%NA35) 
                     IF (ITYPE == 1) THEN
                        CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,          &
                                             'PS', 2,IBLOCK,0,0,                &
                                             inv%MASKG,msh%IENG,msh%XLOCS,msh%ZLOCS, HESS) 
                     ELSE 
                        CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,          &   
                                             'SP', 2,IBLOCK,0,0,                &   
                                             inv%MASKG,msh%IENG,msh%XLOCS,msh%ZLOCS, HESS) 
                     ENDIF
                  ENDDO !loop on types  
               ELSE
                  WRITE(*,*) 'xlbfgs_hsl: Anisotropy not yet done!'
               ENDIF
               DEALLOCATE(HESS) 
               F8 = DBLE(inv%FOBJ) 
               X8(1:inv%NA35) = DBLE(XMOD(1:inv%NA35)) 
               G8(1:inv%NA35) = DBLE(inv%GRAD(1:inv%NA35))
               CALL DCOPY(inv%NA35,X8,1,X0,1)
               CALL DCOPY(inv%NA35,G8,1,G0,1)
            ENDIF
            IF (inv%LGRAD) THEN
               WRITE(*,*) 'xlbfgs_hsl: Copying gradient...'
               CALL SCOPY(inv%NA35,inv%GRAD,1,GRAD,1)
            ENDIF
!
!.......... rescale the pre-conditioned diagonal
            print *, maxval(grad), minval(grad)  
!
!.......... flash initial objective function
            F0 = inv%FOBJ
            WRITE(*,8145) IBLOCK, inv%FOBJ 
 8145       FORMAT(' ++++++++++++++++++++++++++++++++++++++++++++++++',/, &
                   ' +  xlbfgs_hsl: Beginning iterative inversion   +',/, &
                   ' +              Inversion block:',I4,'            +',/, &
                   ' +              Objective function:',E12.4,' +',/, & 
                   ' ++++++++++++++++++++++++++++++++++++++++++++++++',/)

         ENDIF !end check on myid
         NH35  = 1
         NFN35 = 1
         inv%LFUNC = .FALSE. 
         inv%LGRAD = .FALSE.
         LPGNEW = .FALSE. 
         CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR) 
!----------------------------------------------------------------------------------------!
!                         This is the iterative inversion loop                           !
!----------------------------------------------------------------------------------------!
         NFUN = 1 !initialization is 1 function evaluation 
         IFLAG = 0
         ICALPHA = 0 
         DO 2000 K=1,inv%MAXIT 
!
!.......... calculate obj. fn, gradient or hessian?
            IF (LPGNEW) THEN
               CALL FUNCGRADH25(MASTER,MYID,MYNID,IPGROUP,NPGROUPS,                &   
                                MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM, LPGNEW,       &   
                                m1d%VFAST, INV,MSH,SRC,RCV,FRQ, MID, IERR) 
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xlbfgs_hsl: An error occurred calling FUNCGRADH25'
                  WRITE(*,*) 'xlbfgs_hsl: The error was on process:',MYID
                  GOTO 500
               ENDIF
               LPGNEW = .FALSE.
               IF (MYID ==  MASTER) THEN
                  WRITE(*,*) 'xlbfgs_hsl: Factoring gradient preconditioner...'
                  D8(1:NHSIZE) = DBLE(inv%GRADPC(1:inv%NHSIZE)) 
                  CALL DINVHESS(NHSIZE,inv%NA35, inv%NNPINV,inv%NVINV, D8,IPIV,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xlbfgs_hsl: An error occurred in INVHESS'
                     GOTO 500 
                  ENDIF
               ENDIF
            ELSE
               IF (MYID == MASTER .AND. K > 1) &
               WRITE(*,*) 'xlbfgs_hsl: Recyling old diagonal Hessian'
            ENDIF !end check on function/gradient or hessian calculation
!
!.......... call L-BFGS
            IF (MYID == MASTER) THEN
               WRITE(*,*) 'xlbfgs_hsl: Calculating search direction...'
               CALL LBFGS_UPD(NWORK35,NHSIZE,inv%NA35, K,MLBFGS,1,  &
                              inv%NNPINV,inv%NVINV, &
                              IPIV, X8,X0, G8,G0, D8, P8, S,Y, IERR) 
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xlbfgs_hsl: Error calling LBFGS_UPD!'
                  GOTO 500
               ENDIF
               WRITE(*,*) 'xlbfgs_hsl: Averaging search direction...'
               CALL AVGRAD8(NGNOD,msh%NNPG, inv%NA35,msh%NNPG, NGNOD,inv%NCON, &
                            inv%NVINV,inv%NCASC, &
                            inv%MASKG,inv%MCONN,msh%IENG, P8) 
               SEARCH(1:inv%NA35) = SNGL(P8(1:inv%NA35)) 
               CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,          &   
                                    inv%CINVTYPE, 3,IBLOCK,K,0,                &
                                    inv%MASKG,msh%IENG,msh%XLOCS,msh%ZLOCS, SEARCH)
            ENDIF
!----------------------------------------------------------------------------------------!
!                                This is the linesearch loop                             !
!----------------------------------------------------------------------------------------!
!
!.......... begin iterative inversion
            ICALPHA = 0 
            NFEV = 0 
            IF (MYID == MASTER) THEN
               print *, 's8 minimax', minval(p8),maxval(p8)
               !IF (K == 1) THEN
                  STP = 1.D0/MAX(1.D0,MAXVAL(ABS(P8)))*DS0
               !ELSE !try superconvergence?
               !   STP = 1.D0
               !ENDIF
               print *, 'stp0',stp
            ENDIF 
  3000      CONTINUE !back for another
               IF (inv%LFUNC .OR. inv%LGRAD) THEN
                  CALL FUNCGRADH25(MASTER,MYID,MYNID,IPGROUP,NPGROUPS,                &
                                   MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM, LPGNEW,       &
                                   m1d%VFAST, INV,MSH,SRC,RCV,FRQ, MID, IERR) 
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xlbfgs_hsl: An error occurred calling FUNCGRADH25'
                     WRITE(*,*) 'xlbfgs_hsl: The error was on process:',MYID
                     GOTO 500 
                  ENDIF
                  IF (MYID == MASTER) THEN
                     F8 = DBLE(inv%FOBJ)
                     G8(1:inv%NA35) = DBLE(inv%GRAD(1:inv%NA35)) 
                     NFUN = NFUN + 1
                  ENDIF
                  inv%LFUNC = .FALSE.
                  inv%LGRAD = .FALSE. 
               ENDIF
               IF (MYID == MASTER) THEN
                  CALL MCSRCH(inv%NA35,X8,F8,G8,P8, STP,       &
                              FTOL8,GTOL8,XTOL8,STPMIN,STPMAX, &
                              MAXFEV,INFO,NFEV,WA,LP)  
                  XMOD(1:inv%NA35) = SNGL(X8(1:inv%NA35)) !model update
                  MYDEST =-1 !set to an error
                  IF (INFO == 0) THEN !error messages
                     WRITE(*,*) 'Improper input parameters'
                     IERR = 1 
                     GOTO 500
                  ENDIF
                  IF (INFO == 2) THEN
                     WRITE(*,*) 'xlbfgs_hsl: Relative width of uncertainty too small'
                     IERR = 1
                     GOTO 500
                  ENDIF
                  IF (INFO == 3) THEN
                     WRITE(*,*) 'xlbfgs_hsl: More than MAXFEV function evaluations required',MAXFEV
                     IERR = 1
                     GOTO 500
                  ENDIF
                  IF (INFO == 4) THEN
                     WRITE(*,*) 'xlbfgs_hsl: The step size is too small'
                     IERR = 1
                     GOTO 500
                  ENDIF
                  IF (INFO == 5) THEN
                     WRITE(*,*) 'xlbfgs_hsl: The step size is too large'
                     IERR = 1
                     GOTO 500
                  ENDIF
                  IF (INFO == 6) THEN
                     WRITE(*,*) 'xlbfgs_hsl: Rounding errors prevent further progress'
                     IERR = 1
                     GOTO 500  
                  ENDIF
!
!................ convergence test
                  IF (INFO == 1) THEN
                     ICALPHA = 0 
                     MYDEST = 2000 !default, continue with iterative loop 
                     WRITE(*,8151) K,NFUN,STP,inv%FOBJ,100. - inv%FOBJ/F0*100.0
 8151                FORMAT(' +++++++++++++++++++++++++++++++++++++++++++++++++',/, &
                            ' + xlbfgs_hsl: Iteration:',I3,' finished            +',/,&
                            ' +             Function evaluations:',I3,'          +',/,&
                            ' +             Step length:',E12.4,'          +',/, &
                            ' +             Objective function:',E12.4,'   +',/,&
                            ' +             Percent reduction:',F8.3,'        +',/,& 
                            ' +++++++++++++++++++++++++++++++++++++++++++++++++',/)
                     IF (inv%FOBJ/F0*100. < 20.0) THEN
                        WRITE(*,*) 'xlbfgs_hsl: Convergence achieved'
                        MYDEST = 1005
                     ENDIF 
                     !IF (inv%LPGRAD) LPGNEW = .TRUE.  
                     LPGNEW = .FALSE. !recycle preconditioner
                  ENDIF 
!
!................ this is the line search
                  IF (INFO ==-1) THEN 
                     MYDEST = 3000 
                     inv%LFUNC = .TRUE.
                     inv%LGRAD = .TRUE.
                     LPGNEW = .FALSE.
                     ICALPHA = ICALPHA + 1 !this is another step calculation
                     WRITE(*,8150) ICALPHA,NFUN,STP,inv%FOBJ !STPOUT,inv%FOBJ
 8150                FORMAT(' ++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
                            ' +  xlbfgs_hsl: Step length itereration:',I3,'       +',/,&
                            ' +              Function evaluations:',I3,'          +',/,&
                            ' +              Step length:',E12.4,'          +',/,&
                            ' +              Objective function:',E12.4,'   +',/,&
                            ' ++++++++++++++++++++++++++++++++++++++++++++++++++',/)
                  ENDIF
!
!................ slipped through the cracks? impossible, would be a memory error
                  IF (MYDEST ==-1) THEN
                     WRITE(*,*) 'xlbfgs_hsl: Error could not determine MYDEST',MYDEST
                     IERR = 1
                     GOTO 500
                  ENDIF
!
!................ update the models
                  WRITE(*,*) 'xlbfgs_hsl: Updating models...'
                  print *, lpgnew, inv%lgrad,inv%lfunc
                  CALL UPDMOD(PROJNM,LFILES, IBLOCK,K,ICALPHA, XMOD, INV,MSH, IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xlbfgs_hsl: Error updating model!'
                     GOTO 500
                  ENDIF       
                  CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,   &   
                                       inv%CINVTYPE, 1,IBLOCK,K,ICALPHA,            &
                                       inv%MASKG,msh%IENG,msh%XLOCS,msh%ZLOCS, inv%GRAD)

               ENDIF !end check on MYID
               CALL MPI_BCAST(inv%LFUNC,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
               CALL MPI_BCAST(inv%LGRAD,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
               CALL MPI_BCAST(   LPGNEW,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
               CALL MPI_BCAST(MYDEST,   1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
               CALL BCAST_MOD_UPD(MPI_COMM_WORLD,MASTER, MSH)
               IF (MYDEST == 3000) GOTO 3000 
!----------------------------------------------------------------------------------------!
!                             This concludes the line search                             !
!----------------------------------------------------------------------------------------!
            IF (MYNID == MASTER) THEN 
               IF (MYID == MASTER) THEN !update the L-BFGS memory
                  CALL LBFGS_UPD(NWORK35,NHSIZE,inv%NA35, K,MLBFGS,2, &
                                 inv%NNPINV,inv%NVINV, &
                                 IPIV, X8,X0, G8,G0, D8, P8, S,Y, IERR) 
               ENDIF 
               IF (MYID == MASTER .AND. inv%IMODSRC > 0) THEN
                  !WRITE(*,*) 'xlbfgs_hsl: Updating source time functions...'
                  !CALL SRCUPD(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                  !            inv%IMODSRC, msh%AZMOD, inv%EST,inv%OBS, src%SOURCE)
               ENDIF
               IF (MYID == MASTER .AND. inv%IMODREC > 0) THEN
!                 WRITE(*,*) 'xlbfgs_hsl: Updating receiver response functions...'
!                 CALL RECUPD(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
!                             inv%IMODREC, msh%AZMOD, inv%EST,inv%OBS, rcv%RECV)
               ENDIF 
!
!............. write the residual tables
               IF (MYID == MASTER) THEN
                  WRITE(*,*) 'xlbfgs_hsl: Writing residual tables...'
                  CALL RESID_TABLE(NDIM,frq%NFREQ,rcv%NREC, &
                                   NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                                   PROJNM, IBLOCK,K, msh%AZMOD, src%CSIDE, &
                                   inv%LUNWRAP, frq%FREQ,XREC,rcv%YREC, inv%EST,inv%OBS)
!                 WRITE(*,*) 'xlbfgs_hsl: Writing residual RMS tables...'
!                 CALL RTABLE_RMS(NDIM,frq%NFREQ,rcv%NREC,  &
!                                 NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
!                                 PROJNM, IBLOCK,K,2, msh%AZMOD, src%CSIDE, &
!                                 inv%LUNWRAP, frq%FREQ,XREC,rcv%YREC, inv%WGHTS,inv%EST ,inv%OBS)
               ENDIF

!              WRITE(*,*) 'xlbfgs_hsl: Broadcasting receiver response functions...'
!              CALL BCAST_RCV_RESP(MYID,MYHD_COMM,MASTER, .FALSE.,frq%NFREQ,RCV)
!              IF (MYID == MASTER) &
!              WRITE(*,*) 'xlbfgs_hsl: Broadcasting source time function...'
!              CALL BCAST_SRC_STF(MYID,MYHD_COMM,MASTER,  .FALSE.,frq%NFREQ,SRC)
            ENDIF !end check on process ID 
            IF (MYDEST == 1005) GOTO 1005 !we've finished
            IF (MYDEST == 2005) GOTO 2005 !we've hit the iteraiton limit
!----------------------------------------------------------------------------------------!
!                       The iteration has updated, broadcast source                      !
!----------------------------------------------------------------------------------------!
            IF (MYID == MASTER) THEN
               WRITE(*,*) 'xlbfgs_hsl: Saving old model and gradient...'
               CALL DCOPY(inv%NA35,X8,1,X0,1) 
               CALL DCOPY(inv%NA35,G8,1,G0,1)  
            ENDIF
!
!.......... lastly save our estimates
            IF (MYID == MASTER) THEN
               WRITE(*,*) 'xlbfgs_hsl: Writing estimates...'
               ESTFL(1:80) = ' ' 
               CK(1:5) = ' '
               CBLOCK(1:5) = ' '
               WRITE(CK,'(I5)') K 
               WRITE(CBLOCK,'(I5)') IBLOCK
               CK = ADJUSTL(CK) 
               CBLOCK = ADJUSTL(CBLOCK) 
               ESTFL = TRIM(ADJUSTL(PROJNM))//'_est-'//TRIM(CBLOCK)//'-'//TRIM(CK)
               ESTFL = ADJUSTL(ESTFL)
               CALL WTEEST25(ESTFL, frq%NFREQ,rcv%NREC,frq%NFREQ, rcv%NREC,src%NSRC,     &
                             frq%FREQ,inv%EST(1,:,:,:),inv%EST(2,:,:,:),inv%EST(3,:,:,:),&
                             IERR)
            ENDIF
 2000    CONTINUE  !loop on iterations
 2005    CONTINUE !We have reached the max number of iterations
         IF (MYID == MASTER .AND. K == inv%MAXIT) &
         WRITE(*,*) 'xlbfgs_hsl: Iteration limit reached!'
 1005    CONTINUE !this block as converged
!
!....... and while we have them, get the latest STF and RRF estimates
         IF (MYID == MASTER .AND. inv%IMODSRC > 0) THEN
            WRITE(*,*) 'xlbfgs_hsl: Updating source time functions...'
            CALL SRCUPD(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                        inv%LUNWRAP,inv%IMODSRC, msh%AZMOD, inv%EST,inv%OBS, src%SOURCE)
            WRITE(*,*) 'xlbfgs_hsl: Writing new source time functions...'
            CALL WTSTF(PROJNM, frq%NFREQ, frq%NFREQ,src%NSRC, IBLOCK,K, &
                       frq%FREQ,src%SOURCE, IERR)
            WRITE(*,*) 'xlbfgs_hsl: Convolving new STFs...'
            CALL CONVEST_STF(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                             src%SOURCE, inv%EST)
         ENDIF
         IF (MYID == MASTER .AND. inv%IMODREC > 0) THEN
            WRITE(*,*) 'xlbfgs_hsl: Updating receiver response functions...'
            CALL RECUPD(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                        inv%LUNWRAP,inv%IMODREC, msh%AZMOD, inv%EST,inv%OBS, rcv%RECV)
            WRITE(*,*) 'xlbfgs_hsl: Writing new RRFs...'
            CALL WTRECST(PROJNM, NDIM,frq%NFREQ, NDIM,frq%NFREQ,rcv%NREC,  &
                         IBLOCK,K, msh%AZMOD, frq%FREQ,rcv%RECV, IERR)

            WRITE(*,*) 'xlbfgs_hsl: Convolving new RRFs...'
            CALL CONVEST_RRF(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                             rcv%RECV, inv%EST)
         ENDIF 
!
!....... write the final estimates
         IF (MYID == MASTER) THEN
            WRITE(*,*) 'xlbfgs_hsl: Writing estimates...'
            ESTFL(1:80) = ' '
            CBLOCK(1:5) = ' '
            WRITE(CBLOCK,'(I5)') IBLOCK
            CBLOCK = ADJUSTL(CBLOCK)
            ESTFL = TRIM(ADJUSTL(PROJNM))//'_est-'//TRIM(CBLOCK)
            ESTFL = ADJUSTL(ESTFL)
            CALL WTEEST25(ESTFL, frq%NFREQ,rcv%NREC,frq%NFREQ, rcv%NREC,src%NSRC,     &
                          frq%FREQ,inv%EST(1,:,:,:),inv%EST(2,:,:,:),inv%EST(3,:,:,:),&
                          IERR)
         ENDIF
!
!....... free space for next block 
         IF (MYNID == MASTER) THEN
            IF (ASSOCIATED(inv%WGHTS))  DEALLOCATE(inv%WGHTS)
            IF (ASSOCIATED(inv%OBS))    DEALLOCATE(inv%OBS) 
            IF (ASSOCIATED(inv%EST))    DEALLOCATE(inv%EST)  
            IF (ASSOCIATED(src%SOURCE)) DEALLOCATE(src%SOURCE)
         ENDIF
         IF (ASSOCIATED(frq%FREQ))   DEALLOCATE(frq%FREQ) 
 1000 CONTINUE !this is the loop on frequency blocks
!----------------------------------------------------------------------------------------!
!                    This concludes the outer block on frequency groups                  !
!----------------------------------------------------------------------------------------! 
!
!.... break ahead for errors
  500 CONTINUE
      IF (IERR /= 0) THEN
         WRITE(*,8050) MYID
 8050    FORMAT(' -------------------------------------------------------------',/, &
                ' -               An Error was Detected on Process ',I6,'     -',/, &
                ' -                  I will now abort the program             -',/, &
                ' -------------------------------------------------------------',/)
         CALL MPI_ABORT(MPI_COMM_WORLD,30,MPIERR)
      ENDIF
!
!.... free memory
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR) 
      IF (MYID == MASTER) WRITE(*,*) 'xlbfgs_hsl: Freeing memory...'
!.... model
      IF (ASSOCIATED(msh%ECOEFF))    DEALLOCATE(msh%ECOEFF)
      IF (ASSOCIATED(msh%DENS))      DEALLOCATE(msh%DENS)
      IF (ASSOCIATED(msh%XD))        DEALLOCATE(msh%XD)
      IF (ASSOCIATED(msh%ZD))        DEALLOCATE(msh%ZD)
      IF (ASSOCIATED(msh%XLOCS))     DEALLOCATE(msh%XLOCS)
      IF (ASSOCIATED(msh%ZLOCS))     DEALLOCATE(msh%ZLOCS)
      IF (ASSOCIATED(msh%XLOCSE))    DEALLOCATE(msh%XLOCSE)
      IF (ASSOCIATED(msh%ZLOCSE))    DEALLOCATE(msh%ZLOCSE)
      IF (ASSOCIATED(msh%XIPTS))     DEALLOCATE(msh%XIPTS)
      IF (ASSOCIATED(msh%ETAPTS))    DEALLOCATE(msh%ETAPTS) 
      IF (ASSOCIATED(msh%CNNPG))     DEALLOCATE(msh%CNNPG)
      IF (ASSOCIATED(msh%CDOMAIN))   DEALLOCATE(msh%CDOMAIN) 
      IF (ASSOCIATED(msh%CNP))       DEALLOCATE(msh%CNP) 
!.... graph
      IF (ASSOCIATED(msh%LM))        DEALLOCATE(msh%LM)
      IF (ASSOCIATED(msh%IENG))      DEALLOCATE(msh%IENG)
      IF (ASSOCIATED(msh%IDOFSE))    DEALLOCATE(msh%IDOFSE)
      IF (ASSOCIATED(msh%MYDOFS))    DEALLOCATE(msh%MYDOFS)
      IF (ASSOCIATED(msh%MYELEM))    DEALLOCATE(msh%MYELEM)    
      IF (ASSOCIATED(msh%IRPTR_LOC)) DEALLOCATE(msh%IRPTR_LOC)
      IF (ASSOCIATED(msh%JCPTR_LOC)) DEALLOCATE(msh%JCPTR_LOC)
!.... source
      IF (ASSOCIATED(src%PYTAB))   DEALLOCATE(src%PYTAB) 
      IF (MYNID == MASTER) THEN
         IF (ASSOCIATED(src%SOURCE))  DEALLOCATE(src%SOURCE)
         IF (ASSOCIATED(src%SRCTYP))  DEALLOCATE(src%SRCTYP)
         IF (ASSOCIATED(src%BAZN))    DEALLOCATE(src%BAZN)
         IF (ASSOCIATED(src%AOI))     DEALLOCATE(src%AOI)
      ENDIF
!
!.... inversion stuff
      IF (MYNID == MASTER) THEN
         IF (MYID == MASTER) THEN
            DEALLOCATE(XMOD) 
            IF (inv%LPGRAD) THEN
               DEALLOCATE(inv%GRADPC) 
               DEALLOCATE(inv%IPIVH) 
            ENDIF
            DEALLOCATE(DIAG)
            DEALLOCATE(SEARCH) 
            IF (ALLOCATED(GRAD))   DEALLOCATE(GRAD) 
            IF (ALLOCATED(X0))     DEALLOCATE(X0)
            IF (ALLOCATED(G0))     DEALLOCATE(G0)
            IF (ALLOCATED(X8))     DEALLOCATE(X8)
            IF (ALLOCATED(G8))     DEALLOCATE(G8)  
            IF (ALLOCATED(D8))     DEALLOCATE(D8) 
            IF (ALLOCATED(P8))     DEALLOCATE(P8)
            IF (ALLOCATED(Y))      DEALLOCATE(Y)
            IF (ALLOCATED(S))      DEALLOCATE(S) 
            IF (ALLOCATED(WA))     DEALLOCATE(WA) 
            IF (ALLOCATED(XREC))   DEALLOCATE(XREC)
         ENDIF
      ENDIF 
      IF (ASSOCIATED(inv%MASKG))      DEALLOCATE(inv%MASKG)  
      IF (ASSOCIATED(inv%IGPART))     DEALLOCATE(inv%IGPART)
      IF (ASSOCIATED(inv%MCONN))      DEALLOCATE(inv%MCONN) 
      IF (ASSOCIATED(inv%WMASK))      DEALLOCATE(inv%WMASK) 
      IF (ASSOCIATED(inv%MYGRAD))     DEALLOCATE(inv%MYGRAD)
      IF (ASSOCIATED(inv%ICSC_FDIST)) DEALLOCATE(inv%ICSC_FDIST)
      IF (ASSOCIATED(inv%JCSC_FDIST)) DEALLOCATE(inv%JCSC_FDIST) 

      IF (ASSOCIATED(mid%IRN_LOC)) DEALLOCATE(mid%IRN_LOC)
      IF (ASSOCIATED(mid%JCN_LOC)) DEALLOCATE(mid%JCN_LOC)
      IF (ASSOCIATED(mid%A_LOC))   DEALLOCATE(mid%A_LOC) 
      IF (MYNID == MASTER) THEN
         IF (ASSOCIATED(mid%RHS)) DEALLOCATE(MID%RHS) 
      ENDIF
      MID%JOB =-2
      CALL CMUMPS(MID)
!
!.... all done 
      TESIM = MPI_WTIME()
      IF (MYID == MASTER) THEN
         WRITE(*,9405) (TESIM - TSSIM)/3600.D0
 9405    FORMAT(' xlbfgs_hsl: Inversion time:',F8.2,' hours')
      ENDIF
      CALL MPI_FINALIZE(MPIERR)
      STOP
      END 
