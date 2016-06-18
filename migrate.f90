!
!     This routine drives the data migration.  To migrate data we 
!
!       (1) Read the data 
!       (2) Calculate the forward problem for Greens functions at each source 
!           and frequency and offload the Jacobians to disk.  
!       (3) Solve for a source time function for each source
!       (4) Solve adj(J) J m =-adj(J) d where J, for each frequency and source 
!           Here we can solve adj(J) J a few different ways.  
!           a) We can assume it is identity, i.e., the migrated image is just 
!              m =-adj(J) d   
!              but we need to remember to divide the data d by the 
!              magnitude of the source time function estimate
!           b) We can make adj(j) J block diagonal and solve 
!              diag(adj(J) J) m =-adj(J) d 
!              If we are just migrating P or S independtly the block diagonal 
!              matrix will reduce to a diagonal.  
!           c) We can solve adj(J) J m =-adj(J) d a la conjugate gradient method
!           d) We can solve adj(J) J m =-adj(J) d with the preconditioned 
!              conjugate gradient method where our preconditioner is the block 
!              diagonal matrix (see b)  
!           e) We can also solve the overdetermined problem 
!               {J_1} 
!               { . } m = {d} 
!               { . } 
!               { . } 
!               {J_n}
!              via LSQR where n = number of frequencies * number of sources.
!           **Note that at each frequency and source J must be multiplied by the 
!           source time function estimate.
!       (5) The images are then stacked over the frequency band to produce the 
!           final migrated image
!
!     Ben Baker February 2013
!
!.... includes 
      INCLUDE 'source.inc'
      INCLUDE 'mesh_inv.inc'
      INCLUDE 'mesh.inc'
      INCLUDE 'mpif.h'
      INCLUDE 'cmumps_struc.h'
      TYPE (CMUMPS_STRUC) MID
!.... project name
      CHARACTER(80) PROJNM
!.... frequency block
      REAL*8, POINTER :: FREQINV(:) 
      INTEGER*4 NFPASS 
!.... receiver information
      REAL*8, ALLOCATABLE :: XREC(:), ZREC(:) 
      LOGICAL*4, ALLOCATABLE :: LFSURF(:) 
      LOGICAL*4 LINREC
!.... funcgradh stuff
      LOGICAL*4 LFUNC,LGRAD,LHESS,LPHESS
!.... source info and inversion 
      REAL*8 AZMOD,AZTOL,AOITOL
      LOGICAL*4 LSRCEX,LMODSRC,LSRCINV
!.... MPI stuff
      REAL*8 TSSIM, TESIM
      PARAMETER(MASTER = 0)
!.... va35s parameters
      COMPLEX*8, ALLOCATABLE :: CPGRAD(:)  !complex valued preconditioned gradient
      COMPLEX*8, ALLOCATABLE :: CSEARCH(:) !comlex valued search direction
      REAL*4, ALLOCATABLE :: SEARCH(:)     !search direction
      REAL*4, ALLOCATABLE :: PGRAD(:)      !pre-conditioned gradient 
      REAL*4, ALLOCATABLE :: XMIGR(:)      !migration model 
      REAL*4, ALLOCATABLE :: HESS(:)       !for plotting terms in the Hessian
      REAL*4 EPS4   
      logical*4 lgnewt
!.... this will be inescapable 
      COMMON /IFCOMM_BLOCK/ FREQINV,NFPASS 
      LOGICAL*4 LMIGR/.TRUE./ !we are migrating data so add STF and Jacobian multiplies
!.... parameters and local variables
      REAL*8 PI180 
      COMPLEX*8 QN, QE, QZ, SINAZ, COSAZ  
      PARAMETER(PI180 = 0.017453292519943295D0)
!
!----------------------------------------------------------------------------------------!
!
!.... initialize MPI
      CALL MPI_INIT(MPIERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, MPIERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, MPIERR)
      TSSIM = MPI_WTIME()
!
!.... head node reads model information
      IERR = 0
      IF (MYID == MASTER) THEN
         PROJNM(1:80) = ' ' 
 8000    FORMAT(' ---------------------------------------------------------',/, &
                ' -  xmigr25: A massively parallel 2.5D unstructured   -',/, &
                ' -           elastic migration algorithm utilizing    -',/, & 
                ' ---------------------------------------------------------',/)
         projnm = 'simple'
         projnm = 'tester'
         npgroups = 1
         WRITE(*,*) 'xmigr25: Enter project name:' 
         !READ(*,'(A)') PROJNM
         PROJNM = ADJUSTL(PROJNM)
         WRITE(*,*) 'xmigr25: Enter number of process groups:' 
         !READ *, NPGROUPS
         !npgroups = 1 
         IBLOCK = 1
         azmod = 0.d0 !+ west to east
         aztol = 5.d0
         aoitol = 5.d0
         cinvtype = 'PS' !test them both
         pthresh = 2.0  !percent of singular values below which are inflated
         eps4 = 1.e-5
         maxit = 40  !max number of iterations
         norm = 2    !1 -> L1 norm, 2 -> L2 norm (default)
         irestp = 3  !1 -> phase only, 2 -> amplitude only, 3 -> both (default)  
         lphess = .true. !calculate a pseudo hessian for gradient pre-conditioning
         lgnewt = .true.
         lmodsrc = .true.
         imodsrc = 3 !1 -> phase only, 2 -> amplitude only, 3 -> both (default)
         linrec = .true. !include receivers in inversion points?
!....... initial warnings on NPGROUPS
         IF (MOD(NPROCS,NPGROUPS) /= 0) THEN 
            WRITE(*,*) 'xmigr25: Error I cant divide nprocs by npgroups evenly!'
            IERR = 1 
            GOTO 500 
         ENDIF
         NPARTS = NPROCS/NPGROUPS
!
!....... check the inversion scheme
         IF (CINVTYPE == 'PP' .OR. CINVTYPE == 'pp' .OR.      &
             CINVTYPE == 'SS' .OR. CINVTYPE == 'ss' .OR.      &
             CINVTYPE == 'PS' .OR. CINVTYPE == 'ps')  THEN
            NVINV = 2 
            IF (CINVTYPE == 'PP' .OR. CINVTYPE == 'pp' .OR.  &
                CINVTYPE == 'SS' .OR. CINVTYPE == 'ss') THEN
               NVINV = 1
            ENDIF
         ELSE
            WRITE(*,*) 'xmigr25: I do not know what to invert for:',CINVTYPE
            IERR = 1 
            GOTO 500 
         ENDIF 
!
!....... read the mesh
         WRITE(*,*) 'xmigr25: Reading mesh...'
         CALL RDMESH_BKHD(PROJNM, LISISO,NELEM,NNPG,NORD,IITYPE,   &
                          XLATMIN,XLONMIN, XLATMAX,XLONMAX, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xmigr25: Error reading mesh headers'
            IERR = 1 
            GOTO 500 
         ENDIF
         NLXI = NORD + 1 
         NLETA = NLXI  !mesh must be conforming
         NEN = NLXI*NLETA
         ALLOCATE(XIPTS(NLXI))
         ALLOCATE(ETAPTS(NLETA)) 
         CALL GENPTS(IITYPE, NLXI,NLETA, XIPTS,ETAPTS) 
         ALLOCATE(CDOMAIN(NELEM))
         ALLOCATE(CNNPG(NNPG))
         ALLOCATE(IENG(NGNOD,NELEM))
         ALLOCATE(XLOCS(NNPG))
         ALLOCATE(ZLOCS(NNPG))
         ALLOCATE(XD(NNPG)) 
         ALLOCATE(ZD(NNPG)) 
         ALLOCATE(DENS(NNPG))
         IF (LISISO) THEN
            ALLOCATE(ECOEFF(NNPG,2))
         ELSE
            ALLOCATE(ECOEFF(NNPG,5)) 
         ENDIF
         CALL RDMESHBK(NGNOD,NNPG, PROJNM, NNPG,NELEM, NELEME,NABS,        &   
                       CDOMAIN,CNNPG, IENG, XLOCS,ZLOCS,XD,ZD,DENS,        &   
                       ECOEFF,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xbielak25: Error reading mesh files!'
            GOTO 500 
         ENDIF 
         XWIDTH = MAXVAL(XD)
         ZWIDTH = MAXVAL(ZD)
!
!....... initial check on observation file 
         WRITE(*,*) 'xmigr25: Checking headers on observation file...'
         CALL RDTOBS_HD(PROJNM, NDIMOUTT,NFREQ,NREC,NSRC, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xmigr25: There was an error with your observation file'
            GOTO 500
         ENDIF 
!
!....... get the number of blocks
         CALL RDFREQI_BLHD(PROJNM,NBLOCKS,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xmigr25: Error cannot locate the frequency block file'
            GOTO 500
         ENDIF
!
!....... read the source list
         CALL RDSRC_HD(PROJNM, NSRC,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xmigr25: Error reading src header'
            GOTO 500 
         ENDIF
         ALLOCATE(SRCTYP(NSRC))  !character descriptor of source
         ALLOCATE(CSIDE(NSRC))   !'L' -> source from left, 'R' source from right
         ALLOCATE(AOI(NSRC))     !angle of incidence 
         ALLOCATE(BAZN(NSRC))    !back azimuth (correct)
         ALLOCATE(SLAT(NSRC))    !source latitude
         ALLOCATE(SLON(NSRC))    !source longitude
         ALLOCATE(SDEP(NSRC))    !source depth
         ALLOCATE(STRIKE(NSRC))  !strike
         ALLOCATE(DIP(NSRC))     !dip
         ALLOCATE(RAKE(NSRC))    !rake 
         ALLOCATE(SMAG(NSRC))    !magnitude
         ALLOCATE(MODE(NSRC))    !mode to model
         AZMOD = FAZMOD(XLATMIN,XLONMIN, XLATMAX,XLONMAX)
         CALL RDSRC_EQ(PROJNM,NSRC,XLATMIN,XLONMIN,XLATMAX,XLONMAX, &
                       AZMOD, CSIDE,SRCTYP,                         &
                       AOI,BAZN,SLAT,SLON,SDEP,STRIKE,DIP,RAKE,     &
                       SMAG,MODE, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,* 0'xmigr25: Error calling rdsrc_eq!'
            GOTO 500 
         ENDIF
!
!....... check for a receiver file
         CALL RDREC_HD(PROJNM, NREC,IERR)
         IF (IERR /= 0) THEN
            IF (IERR > 0) THEN
               WRITE(*,*) 'xmigr25: Error reading recv header'
               GOTO 500 
            ELSE
               WRITE(*,*) 'xmigr25: There is no receiver file!  Cannot migrate data!'
               IERR = 1 
               GOTO 500
            ENDIF 
         ENDIF
         ALLOCATE(LFSURF(NREC))
         ALLOCATE(XREC(NREC))
         ALLOCATE(YREC(NREC))
         ALLOCATE(ZREC(NREC))
         CALL RDREC(PROJNM, NREC, LFSURF,XREC,YREC,ZREC, IERR)
!
!....... fill the 1D models
         CALL FILL1D()
         WRITE(*,*) 'xmigr25: Splitting process groups...'
      ENDIF
      CALL MKMPI_GROUPS(MASTER,MYID,NPGROUPS,NPROCS,MPI_COMM_WORLD, &
                        IPGROUP,MYNID,MYHD_COMM,MYSLV_COMM)
      IF (MYID == MASTER) WRITE(*,*) 'xmigr25: Initializing MUMPS...'
      MID%COMM = MYSLV_COMM !set communicator
      MID%SYM = 0 !unsymmetric
      MID%JOB =-1 !initialize
      MID%PAR = 1 !host host working 
      CALL CMUMPS(MID)
      MID%ICNTL(3) = 0 !suppress output
      MID%ICNTL(4) = 1 !only error messages
!
!.... generate a graph
      IF (MYID == MASTER) THEN
         WRITE(*,*) 'xmigr25: Generating graph...'
         CALL GEN_GRAPH25(.TRUE.,NELEME,NPARTS, LFSURF,XREC,ZREC) 
         MID%N = NDOF
         MID%NZ = NZERO
         CALL PLOT_MESHP_VTK(PROJNM,NGNOD,NDIM,NEN, NNPG,NELEM,NDOF,   &   
                             NDIM,NEN, PART,IENG,LM, XLOCS,ZLOCS) 
         WRITE(*,*)
         WRITE(*,9408) NORD,NELEM,NABS,NELEME, NNPG,NNPE,NDOF,NZERO  
 9408    FORMAT(' xmigrate25: Polynomial order:'                    ,I4 ,/,        &   
                '             Number of elements in mesh:'          ,I10,/,         &   
                '             Number of absorbing elements:'        ,I8 ,/,         &   
                '             Number of Bielak elements:'           ,I8 ,/,         &   
                '             Number of anchor nodes in mesh:'      ,I10,/,         &   
                '             Number of nodes in Bielak boundary:'  ,I10,/,         &   
                '             Number of degrees of freedom:'        ,I14,/,         &   
                '             Number of non-zeros in global matrix:',I16,/)
         IF (ALLOCATED(LFSURF)) DEALLOCATE(LFSURF) 
         IF (ALLOCATED(XREC)) DEALLOCATE(XREC)
         IF (ALLOCATED(ZREC)) DEALLOCATE(ZREC)
!
!....... create the pointers associated w/ the inverse problem
         WRITE(*,*) 'xmigrate25: Generating mask for inverse problem...'
         CALL GENGMASK(PROJNM, NDIM,NEN,NGNOD, NDOF,NNPG,NELEM,NGNOD,     &
                       NLXI,NLETA,NDIM,LINREC,LKBDRY,  CDOMAIN,CNNPG,PART,LM,IENG, IERR) 
         NA35 = NNPINV*NVINV
         ALLOCATE(GRAD(NA35))
         ALLOCATE(CGRAD(NA35))
         ALLOCATE(XMIGR(NA35)) 
!
!....... also set up the Hessian
         NHSIZE = 0
         IF (LPHESS) THEN
            NBLOCK = NVINV**2        !each submatrix is nvinv x nvinv
            NHSIZE = NNPINV*NBLOCK   !total number of non-zeros
            ALLOCATE(PCHESS(NHSIZE)) !Hessian pre-conditioner
            ALLOCATE(IPIVH(NA35))    !pivots in LU factorization (could use cholesky, 
                                     !but small terms worry me) 
         ENDIF
      ENDIF
!----------------------------------------------------------------------------------------!
!                          Broadcast Inverse Problem Parameters                          !
!----------------------------------------------------------------------------------------!
!
!.... broadcast section 1
      CALL BCAST_INV_P1(MPI_COMM_WORLD,MASTER, PROJNM,CINVTYPE, IITYPE,   &
                        LISISO,LPHESS,LGNEWT,                             &
                        NPARTS,NNPINV,NVINV,NA35,NHSIZE,           &
                        MAXIT,NSRC,NSG,NREC, NORD,NLXI,NLETA,NEN,               &
                        NNPG,NNPE,NELEM,NBLOCKS, NL1D_LT,NL1D_RT,         &
                        NDOF,NZERO,NCON)
      IF (MYID /= MASTER) THEN
         ALLOCATE(XIPTS(NLXI))
         ALLOCATE(ETAPTS(NLETA))
         IF (LISISO) THEN
            NCOEFF = 2   
         ELSE
            NCOEFF = 5 
         ENDIF
         ALLOCATE(ECOEFF(NNPG,NCOEFF))
         ALLOCATE(DENS(NNPG))
         ALLOCATE(XD(NNPG))
         ALLOCATE(ZD(NNPG))
         ALLOCATE(XLOCS(NNPG))
         ALLOCATE(ZLOCS(NNPG))
         ALLOCATE(XLOCSE(NNPE))
         ALLOCATE(ZLOCSE(NNPE))
         ALLOCATE(WMASK(NNPINV)) 

         ALLOCATE(LM(NDIM,NEN,NELEM))
         ALLOCATE(IENG(NGNOD,NELEM))
         ALLOCATE(IDOFSE(NDIM,NNPE)) 
         ALLOCATE(MASKG(NNPG)) 
         ALLOCATE(IGPART(NNPINV)) 
         ALLOCATE(MCONN(NNPG,NCON))

         ALLOCATE(CDOMAIN(NELEM))
         ALLOCATE(CNP(NDOF))
         ALLOCATE(CNNPG(NNPG))

         ALLOCATE(VP1D_LT(NL1D_LT))
         ALLOCATE(VS1D_LT(NL1D_LT))
         ALLOCATE(RH1D_LT(NL1D_LT))
         ALLOCATE( Z1D_LT(NL1D_LT))
         ALLOCATE(VPD_RLLT(NL1D_LT))
         ALLOCATE(VSD_RLLT(NL1D_LT))
         ALLOCATE(ROD_RLLT(NL1D_LT))
         ALLOCATE(HDD_RLLT(NL1D_LT))
         ALLOCATE(VPD_LVLT(NL1D_LT))
         ALLOCATE(VSD_LVLT(NL1D_LT))
         ALLOCATE(ROD_LVLT(NL1D_LT))
         ALLOCATE(HDD_LVLT(NL1D_LT))

         ALLOCATE(VP1D_RT(NL1D_RT))
         ALLOCATE(VS1D_RT(NL1D_RT))
         ALLOCATE(RH1D_RT(NL1D_RT))
         ALLOCATE( Z1D_RT(NL1D_RT))
         ALLOCATE(VPD_RLRT(NL1D_RT))
         ALLOCATE(VSD_RLRT(NL1D_RT))
         ALLOCATE(ROD_RLRT(NL1D_RT))
         ALLOCATE(HDD_RLRT(NL1D_RT))
         ALLOCATE(VPD_LVRT(NL1D_RT))
         ALLOCATE(VSD_LVRT(NL1D_RT))
         ALLOCATE(ROD_LVRT(NL1D_RT))
         ALLOCATE(HDD_LVRT(NL1D_RT))

         ALLOCATE(CSIDE(NSRC))
         ALLOCATE(MRDOF(NDIM,NREC))
         ALLOCATE(YREC(NREC))
 
         ALLOCATE(SRCTYP(NSRC))
         ALLOCATE(BAZN(NSRC))
         ALLOCATE(AOI(NSRC)) 
      ELSE
         WRITE(*,*) 'xmigrate25: Broadcasting model and assembly matrices...'
      ENDIF
      CALL BCAST_INV_MD(MPI_COMM_WORLD,MASTER, .FALSE.,LISISO,             &
                        NDIM,NEN,NGNOD,NNPG,                      &
                        NDIM,NEN,NGNOD,NNPG, NNPE,NELEM,NDOF, NLXI,NLETA,  &
                        NNPINV,NCON, AZMOD,XWIDTH,ZWIDTH,                  &
                        CDOMAIN,CNP, CNNPG, IGPART,MASKG,     &
                        IDOFSE,IENG,LM,   &
                        MCONN, WMASK, XIPTS,ETAPTS, XLOCS,ZLOCS,XD,ZD,     &
                        XLOCSE,ZLOCSE, DENS,ECOEFF) 
      IF (MYNID == MASTER) THEN
         WRITE(*,*) 'xmigrate25: Broadcasting receiver DOFs and y locations...'
         CALL MPI_BCAST(LFILES,1,MPI_LOGICAL, MASTER,MYHD_COMM,MPIERR)
         CALL BCAST_REC(MYHD_COMM,MASTER, NREC,NREC,NDIM, MRDOF,YREC)
         CALL MPI_BCAST(  NORM ,1,MPI_INTEGER, MASTER,MYHD_COMM,MPIERR)
         CALL MPI_BCAST(IRESTP ,1,MPI_INTEGER, MASTER,MYHD_COMM,MPIERR)
      ENDIF
!
!.... figure out the local gradient graph
      IF (MYID == MASTER) WRITE(*,*) 'xmigrate25: Generating gradient graph...'
      CALL GRADPTRS(MYNID,MASTER,MYSLV_COMM, NDIM,NEN, NDIM,NEN,NNPG, LM, IERR)  
      IF (IERR /= 0) THEN
         IF (MYNID == MASTER) WRITE(*,*) 'xmigrate25: An error occurred in gradptrs!'
         GOTO 500
      ENDIF 
!----------------------------------------------------------------------------------------!
!                 Graph distribution.  Note we do not reorder the matrix as MUMPS        !
!                 calls METIS_NodeWND which seems to be a bit more clever than           !
!                 METIS_NodeND when the matrix is distributed                            !
!----------------------------------------------------------------------------------------!
      IF (MYID /= MASTER) THEN
         ALLOCATE(IRPTR(NDOF+1))
         ALLOCATE(JCPTR(NZERO))  
         ALLOCATE(PART(NDOF))
      ELSE
         WRITE(*,*) 'xmigrate25: Broadcasting mesh partition and global assembly...'
      ENDIF
      CALL BCAST_GL_GRAPH(MPI_COMM_WORLD,MASTER, NDOF,NZERO,             &   
                          PART,IRPTR,JCPTR) 
!.... have masters save matrix sizes 
      IF (MYNID == MASTER) THEN
         MID%N = NDOF
         MID%NZ = NZERO
         MID%ICNTL(18) = 3 !matrix distributed by user 
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
      IF (MYID == MASTER) WRITE(*,*) 'xmigrate25: Generating local graphs...'
      CALL ICNDZ_LOC(MYNID,1,NDOF, PART,IRPTR, NDOFL,NZLOC)
      ALLOCATE(MYDOFS(NDOFL))
      ALLOCATE(IRPTR_LOC(NDOFL+1))
      ALLOCATE(JCPTR_LOC(NZLOC))
      CALL PART2CRSL(NZERO,NZLOC, NDOF,NDOFL, MYNID,1,             &   
                     IRPTR,JCPTR,PART, MYDOFS,IRPTR_LOC,JCPTR_LOC) 
      NELEML = ICELEML(MYNID, NDIM,NEN,NDOF, NELEM, 1,             &   
                       NEN,NDIM, PART,LM)
      ALLOCATE(MYELEM(NELEML))
      CALL GENELEML(MYNID, NDIM,NEN,NDOF, NELEM,NELEML, 1,         &   
                    NEN,NDIM, PART,LM, MYELEM)
      DEALLOCATE(IRPTR)
      DEALLOCATE(JCPTR)
      DEALLOCATE(PART)
      MID%NZ_LOC = NZLOC
      ALLOCATE(MID%IRN_LOC(MID%NZ_LOC))
      ALLOCATE(MID%JCN_LOC(MID%NZ_LOC)) 
      CALL CRS2COOLOC(NZLOC,NDOFL, MYDOFS,IRPTR_LOC,JCPTR_LOC,       &   
                      MID%IRN_LOC,MID%JCN_LOC)
      MID%JOB = 1 
      CALL CMUMPS(MID)
      ALLOCATE(MID%A_LOC(MID%NZ_LOC))
      IF (MYNID == MASTER) THEN
         MID%LRHS = MID%N
         ALLOCATE(MID%RHS(MID%N)) 
      ENDIF
!----------------------------------------------------------------------------------------!
!     The first step is to run the forward problem with just the Greens fns so that we   !
!     make estimates at the receiver locations                                           !
!----------------------------------------------------------------------------------------!
!
!.... for the desired block pull out the frequency list 
      IF (MYID == MASTER) THEN  
         WRITE(*,*) 
         WRITE(*,*) 'xmigrate25: Reading frequency block:',IBLOCK
         CALL RDFREQI_FHD(PROJNM,IBLOCK, NOMINV,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xmigrate25: Error rading frequency block header!'
            GOTO 500
         ENDIF
         ALLOCATE(FREQINV(NOMINV))
         CALL RDFREQI(PROJNM,IBLOCK,NOMINV, FREQINV,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xmigrate25: Error reading frequency block!'
            GOTO 500
         ENDIF 
         ALLOCATE(  EST(NDIM,NOMINV,NREC,NSRC)) !response estimates (no STF) 
         ALLOCATE(  OBS(NDIM,NOMINV,NREC,NSRC)) !observations
         ALLOCATE(WGHTS(NDIM,NOMINV,NREC,NSRC)) !data weights
         CALL RDTOBS25(PROJNM, NDIM,NOMINV,NREC,NOMINV, NDIM,NREC,    &
                       NSRC, FREQINV, WGHTS,OBS, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xmigrate25: Error reading observation file!'
            GOTO 500
         ENDIF
         ALLOCATE(SOURCE(NOMINV,NSRC)) 
         SOURCE(1:NOMINV,1:NSRC) = CMPLX(1.D0,0.D0) !will extract the medium response 
      ENDIF 
!
!.... hardwire constants for FUNCGRADH 
      LFUNC = .FALSE.  !don't need a function evaluation 
      LGRAD = .FALSE.  !I'll set the RHS externally to control residuals 
      LPHESS = .TRUE.  !want a pre-conditioner on hand
      LGNEWT = .TRUE.  !need Jacobians for solving LSQR problem 
!
!.... now run the forward problem, and save the  
      CALL FUNCGRADH(MASTER,MYID,MYNID,IPGROUP,NPGROUPS,            &
                     MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM,           &
                     LFUNC,LGRAD,LPHESS,LGNEWT, PROJNM,NORM,IRESTP, &
                     AZMOD, MID, IERR) 
      IF (MYID == MASTER) THEN
         WRITE(*,*) 'xmigrate25: Error calling funcgradh!'
         GOTO 500 
      ENDIF
!
!.... fill data vector
      IF (MYID == MASTER) THEN
         INDX = 0
         NWORK = 3*NOMINV*NSRC*NREC 
         WRITE(*,*) 'xmigrate25: Filling RHS vector...'
         COSAZ = CMPLX(DCOS(AZMOD*PI180),0.D0)
         SINAZ = CMPLX(DSIN(AZMOD*PI180),0.D0)
         DO IFREQ=1,NOMINV
            CALL FILL_SRCPRM(IFREQ,AZTOL,AOITOL)
            DO ISG=1,NSG 
               DO JSRC=1ISGPTR(ISG),ISGPTR(ISG+1)-1
                  ISRC = ISRCPRM(JSRC)
                  DO IREC=1,NREC
                     QN = OBS(1,IFREQ,IREC,ISRC)
                     QE = OBS(2,IFREQ,IREC,ISRC)
                     QZ = OBS(3,IFREQ,IREC,ISRC)
                     P(INDX+1) = QN*COSAZ + QE*SINAZ
                     P(INDX+2) =-QN*SINAZ + QE*COSAZ
                     P(INDX+3) = QZ
                     IF (SRCTYP(1:1) == 'S') THEN
                        IF (SRCTYP(2:2) == 'L') THEN !rayleigh
                           P(INDX+2) = CMPLX(0.0,0.0) 
                        ELSEIF (SRCTYP(2:2) == 'R') THEN !love
                           P(INDX+1) = CMPLX(0.0,0.0)
                           P(INDX+3) = CMPLX(0.0,0.0)
                        ELSEIF (SRCTYP(2:2) == 'V') THEN !vertical   
                           P(INDX+1) = CMPLX(0.0,0.0)
                           P(INDX+2) = CMPLX(0.0,0.0)
                        ENDIF
                     ENDIF
                     INDX = INDX = 3
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
!.... set RHS
      CALL JAC25_DRIVER(MYID,MYNID,MASTER,MYHD_COMM,MYSLV_COMM, &
                        .TRUE.,IPGROUP,NPGROUPS, 3,P,?,IERR)
!
!.... master process estimates the source time function
      IF (MYID == MASTER) THEN
         CALL SRCUPD(NDIM,NOMINV,NREC, NDIM,NOMINV,NREC,NSRC, IRESTP, &
                     EST,OBS, SOURCE) 
      ENDIF
!
!.... invert the diagonal hessian
      IF (MYID == MASTER) THEN
         CALL INVHESS(NHSIZE,NA35, NNPINV,NVINV, PCHESS,IPIVH,IERR)
         IF (IERR /= 0) THEN 
            WRITE(*,*) 'xmigrate25: An error occurred in invhess!'
            GOTO 500
         ENDIF  
      ENDIF 
!
!.... broadcast STF to all process
      IF (.NOT.ALLOCATED(SOURCE)) ALLOCATE(SOURCE(NOMINV,NSRC)) 
      DO 100 IFREQ=1,NFREQ
         CALL MPI_BCAST(SOURCE(IFREQ,1:NSRC),NSRC,MPI_DOUBLE_COMPLEX,  &
                        MASTER,MPI_COMM_WORLD,MPIERR)
  100 CONTINUE  
!
!.... solve gauss-newton problem
      CALL CGITER_CF(MYID,MYNID,MASTER,MYHD_COMM,MYSLV_COMM,               &
                     MPI_COMM_WORLD,IPGROUP,NPGROUPS, NNPINV,NVINV,NHSIZE, &
                     LMIGR, IPIVH,PCHESS,NA35,.FALSE., CGRAD,XMIGR,IERR) 
!----------------------------------------------------------------------------------------!
!     This is the loop on frequency blocks.  To help linearize the inversion we adopt    !
!     a multiscale strategy.  The option is certainly here, but I don't recommend        !
!     using it and would instead recommend shell scripts because once an artifcact is    !
!     introduced it is not coming out.                                                   ! 
!----------------------------------------------------------------------------------------!
!
!.... free memory
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR) 
      IF (MYID == MASTER) WRITE(*,*) 'xlbfgs_hsl: Freeing memory...'
      DEALLOCATE(ECOEFF) 
      DEALLOCATE(DENS)
      DEALLOCATE(XD)
      DEALLOCATE(ZD)
      DEALLOCATE(XLOCS)
      DEALLOCATE(ZLOCS)
      DEALLOCATE(XLOCSE)
      DEALLOCATE(ZLOCSE)
      DEALLOCATE(XIPTS)
      DEALLOCATE(ETAPTS)

      DEALLOCATE(LM)
      DEALLOCATE(IENG)
      DEALLOCATE(IDOFSE)
      DEALLOCATE(MYDOFS)
      DEALLOCATE(MYELEM)
      DEALLOCATE(IRPTR_LOC)
      DEALLOCATE(JCPTR_LOC)

      DEALLOCATE(CDOMAIN)
      DEALLOCATE(CNP)
      DEALLOCATE(CNNPG)
      IF (MYNID == MASTER) THEN
         DEALLOCATE(VP1D_LT)
         DEALLOCATE(VS1D_LT)
         DEALLOCATE(RH1D_LT)
         DEALLOCATE( Z1D_LT)

         DEALLOCATE(VP1D_RT)
         DEALLOCATE(VS1D_RT)
         DEALLOCATE(RH1D_RT)
         DEALLOCATE( Z1D_RT)
         DEALLOCATE(CSIDE)
         DEALLOCATE(MRDOF)
         DEALLOCATE(YREC)

         DEALLOCATE(SRCTYP)
         DEALLOCATE(BAZN) 
         DEALLOCATE(AOI)
         IF (MYID == MASTER) THEN
            DEALLOCATE(XMIGR) 
            DEALLOCATE(GRAD)
            DEALLOCATE(CGRAD) 
            IF (LPHESS) THEN
               DEALLOCATE(PCHESS) 
               DEALLOCATE(IPIVH) 
            ENDIF
            IF (ALLOCATED(CSEARCH)) DEALLOCATE(CSEARCH)
            IF (ALLOCATED(SEARCH))  DEALLOCATE(SEARCH)
            IF (ALLOCATED(PGRAD))   DEALLOCATE(PGRAD) 
            IF (ALLOCATED(CPGRAD))  DEALLOCATE(CPGRAD) 
         ENDIF
      ENDIF 
      DEALLOCATE(IGPART)
      DEALLOCATE(MASKG)  
      DEALLOCATE(ICSC_FDIST)
      DEALLOCATE(JCSC_FDIST) 
      DEALLOCATE(MYGRAD) 

      DEALLOCATE(MID%IRN_LOC)
      DEALLOCATE(MID%JCN_LOC)
      DEALLOCATE(MID%A_LOC) 
      IF (MYNID == MASTER) DEALLOCATE(MID%RHS) 
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
