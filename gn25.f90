!
!     This routine makes use of the Gauss-Newton step and line search of More and 
!     Thuente.  Nominally, MCSRCH needs to be able to calculate a function and
!     gradient.  Addtionally, we may want a gradient pre-conditioner.  This is handled 
!     by funcgradh.f90. - March 2013 (Happy St. Patty's Day) 
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
      TYPE (WIN_INFO) WIN
      TYPE (MOD1D_INFO) M1D
!.... project name
      CHARACTER(80) PROJNM
!.... receiver information
      REAL*8, ALLOCATABLE :: XREC(:), ZREC(:) 
      LOGICAL*4, ALLOCATABLE :: LFSURF(:) 
      LOGICAL*4 LINREC
      LOGICAL*4 LPGNEW
!.... MPI stuff
      REAL*8 TSSIM, TESIM
      INTEGER*4 STAT(MPI_STATUS_SIZE), MASTER,MYID,MYNID, &
                MYSLV_COMM,MYHD_COMM,IPGROUP,MPIERR, &
                NPPGRP, NPGROUPS, NPROCS, NPARTS, MYDEST, MYTAG, &
                IDEST, IPROCS  
      PARAMETER(MASTER = 0)
!.... BLACS stuff
      REAL*8, ALLOCATABLE :: HLOC(:,:) !local Re{ adj(J) J } matrix
      !REAL*8, ALLOCATABLE :: HFAC(:,:) !local factored Re{ adj(J) J } matrix
      REAL*8, ALLOCATABLE :: B(:,:)    !local RHS for Bv = g
      REAL*8, ALLOCATABLE :: DIAG(:)   !diagonal of Re{ adj(J) J} for debugging
      REAL*8, ALLOCATABLE :: COVM(:,:) 
      REAL*8 RLAM          !for padding diagonal
      INTEGER*4 DESCH(9)   !BLACS descriptor
      INTEGER*4 DESCB(9)   !Vector descriptor for RHS
      INTEGER*4 ICTXT      !BLACS context 
      INTEGER*4 NPROW      !number of process rows; nprow*npcol = nprocs
      INTEGER*4 NPCOL      !number of process columns; nprow*npcol = nprocs
      INTEGER*4 MYROW      !processes' row in BLACS grid
      INTEGER*4 MYCOL      !processes' column in BLACS grid        
      INTEGER*4 MBJ, NBJ   !row and column block sizes for jacobian
      INTEGER*4 MBJT,NBJT  !row and column block sizes for adjoint jacobian
      INTEGER*4 MBH ,NBH   !row and column block sizes for Hessian
      INTEGER*4 MH, NH     !number of local rows/columns to be held in HLOC
      INTEGER*4 NUMROC     !calculates NUmber of Rows Or Columns
      INTEGER*4 M, N       !number of rows and columns of Jacobian
      INTEGER*4 NBHEAD     !ID on MPI_COMM_WORLD of (0,0) process on BLACS grid
      INTEGER*4 LDB, LDH   !leading dimensions for HLOC and B
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &    
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 ) 
!.... miscellaneous stuff
      CHARACTER(80) ESTFL, TMPDIR
      CHARACTER(5) CBLOCK 
      COMPLEX*8, ALLOCATABLE :: OBS(:,:,:,:), EST(:,:,:,:) 
      INTEGER*4 NSRC_HD, NFREQ_HD, NREC_HD, NDIM_HD, NELEME, NABS, IOMINV, &
                IBLOCK, NBLOCKS, K, NFUN, NFEV, ICALPHA, ISRC, ISTOP_PT, IERR  
      LOGICAL*4 LUPDSRC, LUPDREC, LCGRNS, LSRCEX, LRECEX, LBLACS, LINIT, LSTF, &
                LFILES, LPRIOR, LKBDRY, LNSRF
      REAL*4 COBJ4, EPSS, THRESHS, F0
      PARAMETER(EPSS = 0.2, THRESHS = 0.2) 
!.... line search variables  
      REAL*4, ALLOCATABLE :: XMOD(:)    !model to invert for
      REAL*4, ALLOCATABLE :: GRAD(:)    !gradient 
      REAL*4, ALLOCATABLE :: SEARCH(:)  !search direction, solution of H p =-g
      REAL*4, ALLOCATABLE :: HESS(:)    !for plotting hessian
      REAL*8, ALLOCATABLE :: SEARCH8(:) !double precision search direction
      REAL*8, ALLOCATABLE :: DX(:)      !x_n - x_prior
      REAL*8, ALLOCATABLE :: XPRIOR(:)  !a prior moddel
      REAL*8, ALLOCATABLE :: X8(:)      !double precision model for MCSRC
      REAL*8, ALLOCATABLE :: P8(:)      !double precision search direction for MCSRCH
      REAL*8, ALLOCATABLE :: G8(:)      !double precision gradient for MCSRCH 
      REAL*8, ALLOCATABLE :: WA(:)      !workspace for MCSRCH 
      REAL*8 STPMIN, STPMAX, F8, EPS8, GTOL8, XTOL8, FTOL8, STP, DS0 
      INTEGER*4 MAXFEV, LP, INFO 
      PARAMETER(LP = 6) 
!.... functions
      INTEGER*4 ICELEML
      real*8, allocatable :: grad8(:), p1(:), p2(:)
      parameter(nrec_list = 3) 
      integer*4 irec_list(nrec_list)
      data irec_list/3,15,30/
      REAL*8 DELTA, SIGMAX, SIGMAZ, VPVAR
      integer*4 descc(9), kernel
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
      inv%LGNEWT = .TRUE.  !definiately perform a Gauss-Newton step
      inv%LPGRAD = .FALSE. !no preconditioner
      LBLACS = .FALSE.     !BLACS grid not yet initiated
!
!.... head node reads model information
      IERR = 0
      IF (MYID == MASTER) THEN
         PROJNM(1:80) = ' ' 
         WRITE(*,8000) 
 8000    FORMAT(' -----------------------------------------------------',/, &
                ' -  xgn25: A massively parallel 2.5D unstructured    -',/, &
                ' -         Gauss Newton elastic inversion algorithm  -',/, &
                ' -----------------------------------------------------',/)
         WRITE(*,*) 'xgn25: Enter project name:' 
         READ(*,'(A)') PROJNM
         PROJNM = ADJUSTL(PROJNM)
         WRITE(*,*) 'xgn25: Enter number of process groups:' 
         READ *, NPGROUPS
!
!....... get variables from the .ini files
         CALL READ_FWD_INI(PROJNM,.TRUE., TMPDIR,LFILES,LNSRF, MSH, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xbielak25: Error I cannot locate your spec file!'
            GOTO 500 
         ENDIF
!        msh%freq0 = 20.d0 
!        msh%aztol = 10.d0
!        msh%aoitol = 5.d0
         inv%cinvtype = 'PP'
!....... line search parameters
         eps8 = 1.d-5
         gtol8 = 9.d-1  
         xtol8 = EPSILON(1.D0)*100.D0
         ftol8 = 1.d-3 
         stpmin = 1.d-20
         stpmax = 1.d+20
         delta = dsqrt( 6000.d0*5000.d0 ) !2500 -> ( 50 m/s )^2, 5000~# of inversion nodes
         sigmax = 2.50d3  !characteristic distance in Gaussian kernel
         sigmaz = sigmax
         vpvar = 100.d0  !Vp model variance
         kernel = 1      !gaussian kernel
         maxfev = 6  
         lcgrns = .false. !true calculate new greens fns, false read old ones
         inv%maxit = 5  !max number of iterations
         inv%norm = 2  !1 -> L1 norm, 2 -> L2 norm (default)
         inv%lunwrap = .false. !try working with unwrapped phase in data
         inv%irestp = 3  !1 -> phase only, 2 -> amplitude only, 3 -> both (default) 
         inv%ibphase = 0 !no geometric spreading correction in backpropagation 
         inv%lpgrad = .FALSE. !calculate a pseudo hessian for gradient pre-conditioning
         inv%ncasc = 0
         inv%imodsrc = 3 !1 -> phase only, 2 -> amplitude only, 3 -> both (default)
         inv%imodrec = 0 !1 -> phase only, 2 -> amplitude only, 3 -> both (default)
         linrec = .true. !include receivers in inversion points?
         lkbdry = .true. !include nodes on interior/bielak bdry
         ds0 = 50.d0     !initial step will allow this for a max step, m/s
         lprior = .false. !no prior models
         istop_pt = 0  !stop point n/a here
!....... initial warnings on NPGROUPS
         IF (MOD(NPROCS,NPGROUPS) /= 0) THEN 
            WRITE(*,*) 'xgn25: Error I cant divide nprocs by npgroups evenly!'
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
            WRITE(*,*) 'xgn25: I do not know what to invert for:',inv%CINVTYPE
            IERR = 1 
            GOTO 500 
         ENDIF 
!
!....... read the mesh
         WRITE(*,*) 'xgn25: Reading mesh...'
         CALL RDMESHBK(PROJNM, NELEME,NABS, MSH, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xgn25: Error reading mesh!'
            IERR = 1 
            GOTO 500 
         ENDIF
         ALLOCATE(msh%XIPTS (msh%NLXI))
         ALLOCATE(msh%ETAPTS(msh%NLETA)) 
         CALL GENPTS(msh%IITYPE, msh%NLXI,msh%NLETA, msh%XIPTS,msh%ETAPTS) 
!
!....... initial check on observation file 
         WRITE(*,*) 'xgn25: Checking headers on observation file...'
         CALL RDTOBS_HD(PROJNM, NDIM_HD,NFREQ_HD,NREC_HD,NSRC_HD, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xgn25: There was an error with your observation file'
            GOTO 500
         ENDIF 
         IF (NDIM_HD /= NDIM) THEN
            WRITE(*,*) 'xgn25: Error component number mismatch!',NDIM_HD,NDIM
            IERR = 1
            GOTO 500
         ENDIF
         IF (NFREQ_HD <= 0) THEN
            WRITE(*,*) 'xgn25: Error observations have 0 frequencies!',NFREQ_HD
            IERR = 1
            GOTO 500
         ENDIF
!
!....... get the number of blocks
         CALL RD_JFREQ_INV_HD(PROJNM, NBLOCKS,IERR)
         IF (IERR /= 0) THEN 
            WRITE(*,*) 'xjoint: Error cannot locate the frequency block file'
            GOTO 500
         ENDIF
         !CALL RDFREQI_BLHD(PROJNM,NBLOCKS,IERR)
         !IF (IERR /= 0) THEN
         !   WRITE(*,*) 'xgn25: Error cannot locate the frequency block file'
         !   GOTO 500
         !ENDIF
!
!....... read the source list
         CALL RDSRC_EQ(PROJNM, msh%XLATMIN,msh%XLONMIN,msh%XLATMAX,msh%XLONMAX, &
                       msh%AZMOD, SRC,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xgn25: Error calling rdsrc_eq!'
            GOTO 500 
         ENDIF
         IF (NSRC_HD /= src%NSRC) THEN
            WRITE(*,*) 'xgn25: Error source number mismatch!',NSRC_HD,src%NSRC
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
               WRITE(*,*) 'xgn25: Error reading recv header'
               GOTO 500 
            ELSE
               WRITE(*,*) 'xgn25: There is no receiver file!  Cannot make estimates!'
               IERR = 1 
               GOTO 500
            ENDIF 
         ENDIF 
         IF (rcv%NREC <= 0) THEN
            WRITE(*,*) 'xgn25: There are no receivers!  Aborting!'
            IERR = 1
            GOTO 500
         ENDIF
         IF (NREC_HD /= rcv%NREC) THEN
            WRITE(*,*) 'xgn25: Error recevier number mismatch!',NREC_HD,rcv%NREC  
            IERR = 1
            GOTO 500
         ENDIF
         ALLOCATE(  LFSURF(rcv%NREC))
         ALLOCATE(    XREC(rcv%NREC))
         ALLOCATE(rcv%YREC(rcv%NREC))
         ALLOCATE(    ZREC(rcv%NREC))
         CALL RDREC(PROJNM, rcv%NREC, LFSURF,XREC,rcv%YREC,ZREC, IERR)
         WRITE(*,*) 'xgn25: Splitting process groups...' 
         ALLOCATE(inv%DX(rcv%NREC-1))
         inv%DX(1:rcv%NREC-1) = 0.D0 
         IF (inv%LCOVD_SRF .OR. inv%LCOVD_BDY) THEN 
            CALL REC_SPACE1D(rcv%NREC,XREC, inv%DX,IERR)
            IF (IERR /= 0) THEN 
               WRITE(*,*) 'xgn25: Error detected in rec_space1d'
               IF (inv%LCOVD_SRF) WRITE(*,*) 'xgn25: Setting lcovd_srf to false'
               IF (inv%LCOVD_BDY) WRITE(*,*) 'xgn25: Setting lcovd_bdy to false'
            ENDIF
         ENDIF
      ENDIF
      CALL MKMPI_GROUPS(MASTER,MYID,NPGROUPS,NPROCS,MPI_COMM_WORLD, &
                        IPGROUP,MYNID,MYHD_COMM,MYSLV_COMM)
!----------------------------------------------------------------------------------------!
!     MUMPS initialization phase and graph reordering utilities.  The mesh should not    !
!     have to change.  If it does the user should re-mesh then re-run the inversion      !
!     software.                                                                          !
!----------------------------------------------------------------------------------------!
      IF (MYID == MASTER) WRITE(*,*) 'xgn25: Initializing MUMPS...'
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
         WRITE(*,*) 'xgn25: Generating graph...'
         CALL GEN_GRAPH25(.TRUE.,NPARTS, LFSURF,XREC,ZREC, MSH,RCV,IERR) 
         IF (inv%LUNWRAP) THEN 
            WRITE(*,*) 'xgn25: Generating free surface information...'
            CALL GEN_GIDOFFS(MSH,IERR)
            IF (IERR /= 0) THEN 
               WRITE(*,*) 'xgn25: Error calling gen_gidoffs!'
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
 9408    FORMAT(' xgn25: Polynomial order:'                    ,I4 ,/,        &   
                '        Number of elements in mesh:'          ,I10,/,         &   
                '        Number of absorbing elements:'        ,I8 ,/,         &   
                '        Number of Bielak elements:'           ,I8 ,/,         &   
                '        Number of anchor nodes in mesh:'      ,I10,/,         &   
                '        Number of nodes in Bielak boundary:'  ,I10,/,         &   
                '        Number of degrees of freedom:'        ,I14,/,         &   
                '        Number of non-zeros in global matrix:',I16,/)
         IF (ALLOCATED(LFSURF)) DEALLOCATE(LFSURF) 
         IF (ALLOCATED(ZREC)) DEALLOCATE(ZREC)
!
!....... create the pointers associated w/ the inverse problem
         WRITE(*,*) 'xgn25: Generating mask for inverse problem...'
         CALL GENGMASK(PROJNM, NDIM,msh%NEN,NGNOD, msh%NDOF,msh%NNPG,msh%NELEM,     &
                       msh%NLXI,msh%NLETA,LINREC,LKBDRY,                            &
                       msh%CDOMAIN,msh%CNNPG,msh%PART,                              &
                       msh%LM,msh%IENG, RCV,INV, IERR) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xgn25: Error generating gradient mask!'
            GOTO 500
         ENDIF
         CALL ELEM_WTS(MSH,INV,IERR)  
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xgn25: Error calling ELEM_WTS'
            GOTO 500
         ENDIF
         inv%NA35 = inv%NNPINV*inv%NVINV
         inv%MCJLOC = rcv%NREC*NDIM
         ALLOCATE(inv%GRAD(inv%NA35))
         ALLOCATE(XMOD(inv%NA35)) 
         ALLOCATE(SEARCH(inv%NA35)) !search direction, solution of Gauss-Newton system
         ALLOCATE(WA(inv%NA35))
         ALLOCATE(X8(inv%NA35)) 
         ALLOCATE(XPRIOR(inv%NA35)) 
         ALLOCATE(G8(inv%NA35))
         ALLOCATE(P8(inv%NA35))
         ALLOCATE(GRAD(inv%NA35)) 
         ALLOCATE(HESS(inv%NA35)) 
      ENDIF
!----------------------------------------------------------------------------------------!
!                          Broadcast Inverse Problem Parameters                          !
!----------------------------------------------------------------------------------------!
      IF (MYID == MASTER) WRITE(*,*) 'xgn25: Broadcasting 2D model...' 
      CALL MPI_BCAST(NPARTS ,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
      
      CALL BCAST_MESH_INFO(MYID,MPI_COMM_WORLD,MASTER, &
                           LNSRF, PROJNM,TMPDIR, MSH)
      IF (MYNID == MASTER) THEN
         IF (MYID == MASTER) WRITE(*,*) 'xgn25: Broadcasting receiver locations...'
         CALL BCAST_RCV_INFO(MYID,MYHD_COMM,MASTER, RCV)
      ENDIF
      IF (MYID == MASTER) WRITE(*,*) 'xgn25: Broadcasting inversion parameters...'
      CALL MPI_BCAST(rcv%NREC,1,MPI_INTEGER, MASTER,MYSLV_COMM,MPIERR) !for jacobian
      CALL BCAST_INV_PARMS(MYID,MPI_COMM_WORLD,MASTER, msh%NNPG,msh%NELEM,NBLOCKS, &
                           rcv%NREC, ISTOP_PT, INV) 
      IF (MYNID == MASTER .AND. inv%LUNWRAP) THEN 
         IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Broadcasting free surface DOFS...'
         CALL BCAST_FS_INFO(MYID,MYHD_COMM,MASTER, MSH) 
      ENDIF
!
!.... figure out the local gradient graph
      IF (MYID == MASTER) WRITE(*,*) 'xgn25: Generating gradient graph...'
      CALL GRADPTRS(MYNID,MASTER,MYSLV_COMM, NDIM,msh%NEN, msh%NEN,msh%NNPG, msh%LM, &
                    INV,IERR)
      IF (IERR /= 0) THEN
         IF (MYNID == MASTER) WRITE(*,*) 'xgn25: An error occurred in gradptrs!'
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
      IF (MYID == MASTER) WRITE(*,*) 'xgn25: Generating local graphs...'
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
      mid%NZ_LOC = msh%NZLOC
      ALLOCATE(mid%IRN_LOC(mid%NZ_LOC))
      ALLOCATE(mid%JCN_LOC(mid%NZ_LOC)) 
      CALL CRS2COOLOC(msh%NZLOC,msh%NDOFL, msh%MYDOFS,msh%IRPTR_LOC,msh%JCPTR_LOC,       &   
                      mid%IRN_LOC,mid%JCN_LOC)
      mid%JOB = 1 
      CALL CMUMPS(MID)
      ALLOCATE(mid%A_LOC(mid%NZ_LOC))
      IF (MYNID == MASTER) THEN
         mid%LRHS = mid%N
         ALLOCATE(mid%RHS(MID%N)) 
      ENDIF
!----------------------------------------------------------------------------------------!
!                   Initialize Scalapack/BLACS for the dense factorization               !
!----------------------------------------------------------------------------------------! 
      IF (MYID == MASTER) THEN
         N = inv%NA35       !adj(J) J will be na35 x na35 
         M = rcv%NREC*NDIM  !assume local Jacobians will be nobs x na35
         M = 2*M            !this holds the real and imaginary parts
         CALL GEN_PGRID(N,N,NPROCS, NPROW,NPCOL,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xanewton: Error calling gen_pgrid!'
            GOTO 500 
         ENDIF
         ! block sizes for Hessian
         !IF (N/NPROW > 64) THEN
         !   MBH = 64
         !ELSE
            MBH = MIN(32,N/NPROW)
         !ENDIF
         !IF (N/NPCOL > 64) THEN
         !   NBH = 64
         !ELSE
            NBH = MIN(32,N/NPCOL)
         !ENDIF
         ! block sizes for Jacobian 
         !IF (M/NPROW > 64) THEN
         !   MBJ = 64
         !ELSE
            MBJ = MIN(32,M/NPROW)
         !ENDIF
         NBJ  = NBH !figured this out with Hmat block size
         MBJT = NBJ !transpose block sizes 
         NBJT = MBJ 
      ENDIF
      CALL MPI_BCAST(M,    1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(N,    1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(NPROW,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(NPCOL,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(MBH  ,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(NBH  ,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(MBJ  ,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(NBJ  ,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(MBJT ,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(NBJT ,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
!
!.... initialize BLACS and set context
      LBLACS = .TRUE.
      CALL SL_INIT(ICTXT, NPROW, NPCOL) !make a process grid in row major ordering   
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL) !set BLACS grid
      IF (MYROW /=-1 .AND. MYCOL /=-1) THEN !i'm i the process grid
         MH = NUMROC(N, MBH,MYROW, 0,NPROW) !number of rows in local H matrix
         NH = NUMROC(N, NBH,MYCOL, 0,NPCOL) !number of columns in local H matrix
         CALL DESCINIT(DESCH,N,N, MBH,NBH, 0,0, ICTXT, MAX0(1,MH),INFO)
         CALL DESCINIT(DESCB,N,1, MBH,  1, 0,0, ICTXT, MAX0(1,MH),INFO)
         LDH = DESCH(LLD_)
         LDB = DESCB(LLD_)
         ALLOCATE(HLOC(LDH,MAX(NH,1)))
         ALLOCATE(B(LDB,1))
         B(1:LDB,1) = 0.D0               !null out, but will be filled later 
      ENDIF
!
!.... let MPI_COMM_WORLD know who BLACS (0,0) grid is
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         NBHEAD = MYID 
         DO 30 IPROCS=1,NPROCS
            IDEST = IPROCS - 1
            IF (IDEST /= MYID) & 
            CALL MPI_SEND(NBHEAD,1,MPI_INTEGER, IDEST,NBHEAD, MPI_COMM_WORLD,MPIERR)
   30    CONTINUE 
      ELSE
         CALL MPI_RECV(NBHEAD,1,MPI_INTEGER, MPI_ANY_TAG,MPI_ANY_SOURCE, &
                       MPI_COMM_WORLD,STAT,MPIERR) 
      ENDIF
      IF (MYID == NBHEAD) THEN
         IF (.NOT.ALLOCATED(GRAD))    ALLOCATE(GRAD(inv%NA35))  
         IF (.NOT.ALLOCATED(SEARCH))  ALLOCATE(SEARCH(inv%NA35)) 
         IF (.NOT.ALLOCATED(SEARCH8)) ALLOCATE(SEARCH8(inv%NA35))
         IF (.NOT.ALLOCATED(HESS))    ALLOCATE(HESS(inv%NA35)) 
         IF (.NOT.ALLOCATED(XMOD))    ALLOCATE(XMOD(inv%NA35)) 
         IF (.NOT.ALLOCATED(X8))      ALLOCATE(X8(inv%NA35))
         IF (.NOT.ALLOCATED(XPRIOR))  ALLOCATE(XPRIOR(inv%NA35)) 
      ELSE
         IF (.NOT.ALLOCATED(GRAD))    ALLOCATE(GRAD(1))
         IF (.NOT.ALLOCATED(SEARCH8)) ALLOCATE(SEARCH8(1)) 
         IF (.NOT.ALLOCATED(SEARCH))  ALLOCATE(SEARCH(1)) 
         IF (.NOT.ALLOCATED(XMOD))    ALLOCATE(XMOD(1)) 
         IF (.NOT.ALLOCATED(X8))      ALLOCATE(X8(1)) 
         IF (.NOT.ALLOCATED(XPRIOR))  ALLOCATE(XPRIOR(1)) 
      ENDIF
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
            WRITE(*,*) 'xgn25: Filling 1D models...'
            CALL FILL1D(MSH,M1D) 
         ENDIF
         CALL MPI_BCAST(msh%XMOD0  ,1,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,MPIERR) 
         CALL MPI_BCAST(msh%XMOD1  ,1,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,MPIERR) 
         CALL MPI_BCAST(msh%XBLKL  ,1,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(msh%XBLKR  ,1,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,MPIERR)
         IF (MYID == MASTER) THEN
            win%DT_SRF = 0.0
            win%DT_BDY = 0.0
            win%START_SRF = 0.0
            win%START_BDY = 0.0
            win%NSAMP_SRF = 0
            win%NSAMP_BDY = 0
            win%NPCT_SRF = 0
            win%NPCT_BDY = 0
            win%LWNDO_BDY = .FALSE.
            win%LWNDO_SRF = .FALSE.
            WRITE(*,*) 'xgn25: Reading frequency block...'
            !CALL RDFREQI(PROJNM,IBLOCK, FRQ,IERR) 
            CALL RD_JFREQ_INV(PROJNM,IBLOCK,win%LWNDO_SRF,win%LWNDO_BDY,  &
                              FRQ,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xgn25: Cannot read inversion frequency block!'
               GOTO 500
            ENDIF   
            WRITE(*,*) 'xgn25: Reading observation file...'
            ALLOCATE(frq%CFTYPE(frq%NFREQ)) 
            IF (src%SRCTYP(1)(1:1) == 'S') THEN
               frq%CFTYPE(1:frq%NFREQ) = 'S'
            ELSE
               frq%CFTYPE(1:frq%NFREQ) = 'B'
            ENDIF
            ALLOCATE(inv%  EST(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
            ALLOCATE(inv%  OBS(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
            ALLOCATE(inv%WGHTS(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
            CALL RDTOBS25(PROJNM, NDIM,frq%NFREQ,rcv%NREC,frq%NFREQ, NDIM,rcv%NREC,    &
                          src%NSRC, inv%LUNWRAP, frq%FREQ, inv%WGHTS,inv%OBS, IERR)
            IF (IERR /= 0) THEN
               DEALLOCATE(inv%OBS)
               DEALLOCATE(inv%WGHTS) 
               WRITE(*,*) 'xgn25: Error reading observation file1'
               GOTO 500
            ENDIF
            WRITE(*,*) 'xgn25: Checking for srcinv file...'
            LUPDSRC = .FALSE. 
            ALLOCATE(src%SOURCE(frq%NFREQ,src%NSRC))
            CALL RDSRCINV(frq%NFREQ,PROJNM, frq%NFREQ,src%NSRC, frq%FREQ, &
                          LSRCEX,src%SOURCE)
            IF (.NOT.LSRCEX) THEN
               LUPDSRC = .TRUE. !you will never succeed without diong this
               WRITE(*,*) 'xgn25: No srcinv file detected'
               IF (inv%IMODSRC == 0) THEN 
                  LUPDSRC = .FALSE. !fine, you asked for it
                  WRITE(*,*) 'xgn25: Warning you are proceeding without an STF estimate'
               ENDIF
            ELSE
               WRITE(*,*) 'xgn25: Source inversion file read'
            ENDIF
            IF (inv%IRESTP == 1) THEN 
               WRITE(*,*) 'xgn25: Normalizing source scale to unity...'
            ELSEIF (inv%IRESTP == 2) THEN 
               WRITE(*,*) 'xgn25: Setting source phase to 0...'
            ENDIF 
            CALL MODSRC(frq%NFREQ, frq%NFREQ,src%NSRC,inv%IRESTP, src%SOURCE) 
            LUPDREC = .FALSE.
            ALLOCATE(rcv%RECV(NDIM,frq%NFREQ,rcv%NREC)) 
            CALL RDREC_RESP(PROJNM, NDIM,frq%NFREQ, NDIM,frq%NFREQ,rcv%NREC, msh%AZMOD, &
                            frq%FREQ, LRECEX,rcv%RECV)
            IF (.NOT.LRECEX) THEN
               LUPDREC = .TRUE.
               IF (inv%IMODREC == 0) THEN
                  LUPDREC = .FALSE. !not as much of a sin as forgetting the STF
                  WRITE(*,*) 'xgn25: Note you are proceeding without an RRF estimate'
               ENDIF 
               IERR = 0 
            ENDIF
            IF (inv%IRESTP == 1) THEN
               WRITE(*,*) 'xgn25: Normalziing receiver responses to unity...' 
            ELSEIF (inv%IRESTP == 2) THEN
               WRITE(*,*) 'xgn25: Setting receiver response phase to 0...'
            ENDIF
            CALL MODREC(NDIM,frq%NFREQ, NDIM,frq%NFREQ,rcv%NREC, inv%IRESTP, rcv%RECV)

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
         IF (MYID == MASTER) WRITE(*,*) 'xgn25: Broadcasting 1D models...'
         CALL BCAST_MOD1D_INFO(MYID,MPI_COMM_WORLD,MASTER, MSH,M1D)
         IF (MYID == MASTER) WRITE(*,*) 'xgn25: Broadcasting source details...'
         CALL BCAST_SRC_INFO(MYID,MPI_COMM_WORLD,MASTER, SRC)
         IF (MYID == MASTER) WRITE(*,*) 'xgn25: Broadcasting frequency information...'
         CALL BCAST_FRQ_INFO(MYID,MPI_COMM_WORLD,MASTER, FRQ)
         IF (MYNID == MASTER) THEN
            IF (MYID == MASTER) &
            WRITE(*,*) 'xgn25: Broadcasting observations...'
            CALL BCAST_OBS_INFO(MYID,MYHD_COMM,MASTER, frq%NFREQ,rcv%NREC,src%NSRC, &
                                INV) 
            IF (MYID /= MASTER) ALLOCATE(inv%EST(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
            WRITE(*,*) 'xgn25: Broadcasting receiver response functions...'
            CALL BCAST_RCV_RESP(MYID,MYHD_COMM,MASTER, .TRUE.,frq%NFREQ,RCV)
            IF (MYID == MASTER) &
            WRITE(*,*) 'xgn25: Broadcasting source time function...'
            CALL BCAST_SRC_STF(MYID,MYHD_COMM,MASTER,  .TRUE.,frq%NFREQ,SRC)  
         ENDIF
!----------------------------------------------------------------------------------------!
!                            Calculate 1D models solutions                               !
!----------------------------------------------------------------------------------------!
         IF (MYID == MASTER) WRITE(*,*) 'xgn25: Calculating 1D solutions...'
         CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
         ALLOCATE(src%PYTAB(frq%NFREQ,src%NSRC))
         src%PYTAB(1:frq%NFREQ,1:src%NSRC) = 0.D0
         CALL GEN_GRNS(MYID,MASTER,MPI_COMM_WORLD, TMPDIR,LCGRNS,LNSRF, &
                       SRC,M1D,FRQ,MSH, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xgn25: Error in gen_grns!'
            GOTO 500
         ENDIF
!....... source
         IF (IBLOCK == NBLOCKS) THEN
            IF (MYNID /= MASTER) THEN
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
         CALL MPI_BCAST(LUPDSRC,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(LUPDREC,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
         inv%LFUNC = .TRUE.   !Most certainly need a function evaluation
         inv%LGRAD = .TRUE.   !MCSRCH will take a gradient first go around
         inv%LGNEWT = .TRUE.  !Set Jacobians for Gauss-Newton step
         LPGNEW = .FALSE.     !dont have funcgradh set a preconditioner
         IF (LUPDREC .OR. LUPDSRC) inv%LGRAD = .FALSE.
         IF (MYID == MASTER) THEN
            IF (inv%LGRAD) THEN
               WRITE(*,*) 'xgn25: Initializing objective function, gradient, and ',  &
                          'Jacobians...'
            ELSE
               WRITE(*,*) 'xgn25: Initializing objective function and Jacobians...'
            ENDIF
         ENDIF
         F0 = 0.0 

         call mpi_bcast(kernel,1,mpi_integer         , master,mpi_comm_world,mpierr)
         call mpi_bcast(sigmax,1,mpi_double_precision, master,mpi_comm_world,mpierr)
         call mpi_bcast(sigmaz,1,mpi_double_precision, master,mpi_comm_world,mpierr)
         call mpi_bcast(vpvar, 1,mpi_double_precision, master,mpi_comm_world,mpierr)
         if (myrow >= 0 .and. mycol >= 0) then
            ldc = ldh
            allocate(covm(ldc,MAX(NH,1)))
            descc(1:9) = desch(1:9)
            print *, kernel,sigma,vpvar
            call covmodel_dist(ldc, nprow,npcol,ictxt,  &
                               kernel,sigmax,sigmaz,vpvar, descc, &
                               inv,msh, covm,ierr) 
            print *, minval(covm),maxval(covm) 
         endif  
!----------------------------------------------------------------------------------------!
!                               Begin the iterative loop                                 !
!----------------------------------------------------------------------------------------! 
         NFUN = 1 
         LINIT = .TRUE.  !initialization
         LSTF  = .FALSE. !assume Jacobians do not STFs
         DO 2000 K=1,inv%MAXIT

!           if (.not.allocated(p1)) then
!              if (myid == master) then
!                 search(:) = 1.0 
!                 allocate(p1(inv%na35))
!                 call scopy(inv%na35,search,1,p1,1)
!              else
!                 allocate(p1(1))
!              endif
!           endif
!           if (.not.allocated(p2)) then
!              if (myid == master) then
!                 allocate(p2(inv%na35))
!              else
!                 allocate(p2(1))
!              endif
!           endif

!           CALL QPMFREE(MASTER,MYID,MYNID,IPGROUP,NPGROUPS,    &
!                        MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM,   &
!                        m1d%VFAST,P1, INV,MSH,SRC,RCV,FRQ, MID, P2,IERR)

!           goto 500


            CALL FUNCGRADH25(MASTER,MYID,MYNID,IPGROUP,NPGROUPS,               &
                             MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM, LPGNEW,      &
                             m1d%VFAST, INV,MSH,SRC,RCV,FRQ, MID, IERR) 
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xgn25: An error occurred while calling funcgradh!'
               GOTO 500
            ENDIF
!
!.......... may need to get an STF, RRF, and/or gradient
            IF (MYID == MASTER .AND. LINIT) THEN
               IF (LUPDSRC) THEN !update source 
                  WRITE(*,*) 'xgn25: Updating STFs...'
                  CALL SRCUPD(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                              inv%LUNWRAP,inv%IMODSRC, msh%AZMOD,inv%EST,inv%OBS,        &
                              src%SOURCE)
                  IF (inv%IRESTP == 1) THEN
                     WRITE(*,*) 'xgn25: Normalizing source scale to unity...'
                  ELSEIF (inv%IRESTP == 2) THEN
                     WRITE(*,*) 'xgn25: Setting source phase to 0...'
                  ENDIF 
                  CALL MODSRC(frq%NFREQ, frq%NFREQ,src%NSRC,inv%IRESTP, src%SOURCE)
                  WRITE(*,*) 'xgn25: Updating estimates...'
                  CALL CONVEST_STF(NDIM,frq%NFREQ,rcv%NREC,           &
                                   NDIM,frq%NFREQ,rcv%NREC,src%NSRC,  &
                                   src%SOURCE, inv%EST)
                  WRITE(*,*) 'xgn25: Writing STFs...'
                  CALL WTSTF(PROJNM, frq%NFREQ, frq%NFREQ,src%NSRC, IBLOCK,0, &
                             frq%FREQ,src%SOURCE, IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xgn25: Error writing source time functions!'
                     GOTO 500
                  ENDIF
                  LSTF = .TRUE. !need to convolve STFs into Jacobians
               ENDIF
               IF (LUPDREC) THEN !update receiver
                  WRITE(*,*) 'xgn25: Updating RRFs...'
                  CALL RECUPD(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                              inv%LUNWRAP,inv%IMODREC, msh%AZMOD, inv%EST,inv%OBS,       &
                              rcv%RECV)
                  IF (inv%IRESTP == 1) THEN 
                     WRITE(*,*) 'xgn25: Normalizing receiver responses to unity...' 
                  ELSEIF (inv%IRESTP == 2) THEN 
                     WRITE(*,*) 'xgn25: Setting receiver response phase to 0...'
                  ENDIF
                  CALL MODREC(NDIM,frq%NFREQ, NDIM,frq%NFREQ,rcv%NREC,inv%IRESTP,rcv%RECV)
                  WRITE(*,*) 'xgn25: Updating estimates...'
                  CALL CONVEST_RRF(NDIM,frq%NFREQ,rcv%NREC,          &
                                   NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                                   rcv%RECV, inv%EST)
                  WRITE(*,*) 'xgn25: Writing RRFs...'
                  CALL WTRECST(PROJNM, NDIM,frq%NFREQ, NDIM,frq%NFREQ,rcv%NREC,  &
                               IBLOCK,0, msh%AZMOD, frq%FREQ,rcv%RECV, IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xgn25: Error writing receiver estimates!'
                     GOTO 500
                  ENDIF
               ENDIF
               IF (.NOT.inv%LGRAD) THEN !need a gradient
                  WRITE(*,*) 'xgn25: Generating initial gradient...'
                  NPPGRP = NPROCS/NPGROUPS 
                  CALL JLOC2GRAD(NDIM,frq%NFREQ,rcv%NREC,                             &
                                 NDIM,frq%NFREQ,rcv%NREC,inv%NA35, src%NSRC,          &
                                 inv%IRESTP,inv%IBPHASE,NPPGRP,inv%LUNWRAP,LSTF,      &
                                 msh%AZMOD, frq%FREQ, src%SOURCE, inv%EST,inv%OBS,    &
                                 inv%GRAD,IERR)
                  IF (IERR /= 0) THEN 
                     WRITE(*,*) 'xgn25: Error generating gradient!'
                     GOTO 500 
                  ENDIF 
                  WRITE(*,*) 'funcgradh: Smoothing gradient...'
                  CALL AVGRAD(NGNOD,msh%NNPG, inv%NA35,msh%NNPG, &
                              NGNOD,inv%NCON,inv%NVINV,inv%NCASC, &
                              inv%MASKG,inv%MCONN,msh%IENG, inv%GRAD)  
                  CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,               &
                                       inv%CINVTYPE, 1,IBLOCK,1,0,                     &
                                       msh%CDOMAIN,inv%MASKG,msh%IENG,                 &
                                       msh%XLOCS,msh%ZLOCS, inv%GRAD)

               ENDIF
               IF (LUPDSRC .OR. LUPDREC) THEN !need a new objective function
                  ALLOCATE(EST(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
                  ALLOCATE(OBS(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
                  CALL UNWRAP_OEST(NDIM,frq%NFREQ,rcv%NREC,  &
                                   NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                                   inv%OBS,inv%EST, OBS,EST) 
                  WRITE(*,*) 'xgn25: Rotating estimates to (N,E,Z) frame...'
                  CALL ROTEST(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC,   &
                              msh%AZMOD,EST) 
                  WRITE(*,*) 'xgn25: Calculating objective function...'
                  inv%FOBJ = COBJ4(NDIM,frq%NFREQ,rcv%NREC,                              &
                                   frq%NFREQ,rcv%NREC,src%NSRC, NDIM,                    &
                                   .FALSE.,inv%NORM,inv%IRESTP, EPSS,THRESHS,0.0,        &
                                   frq%FREQ, inv%WGHTS,OBS,EST)
                  DEALLOCATE(OBS)
                  DEALLOCATE(EST) 
               ENDIF
               nppgrp = nprocs/npgroups
               call plot_jacob_drive(nrec_list,rcv%nrec, nppgrp,lstf, &
                                     msh,src,frq,inv, irec_list, ierr) 
!
!............. intialize model and gradeint
               CALL GENXMOD(msh%NNPG,inv%NA35,msh%NNPG, inv%CINVTYPE,   &
                            inv%MASKG,msh%DENS,msh%ECOEFF, XMOD,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xgn25: Error setting xmod'
                  GOTO 500
               ENDIF 
!
!............. flash initial objective function
               IF (LINIT) THEN
                  WRITE(*,8145) IBLOCK, inv%FOBJ 
                  LINIT = .FALSE.
 8145             FORMAT(' +++++++++++++++++++++++++++++++++++++++++++',/, &
                         ' +  xgn25: Beginning iterative inversion   +',/, &
                         ' +         Inversion block:',I4,'            +',/, &
                         ' +         Objective function:',E12.4,' +',/, & 
                         ' +++++++++++++++++++++++++++++++++++++++++++',/)
               ENDIF
!
!............. set the initial model, gradient, and objective function
               CALL SCOPY(inv%NA35,inv%GRAD,1,GRAD,1)
               X8(1:inv%NA35) = DBLE(XMOD(1:inv%NA35))
               G8(1:inv%NA35) = DBLE(GRAD(1:inv%NA35))
               IF (LPRIOR) THEN
                  print *, 'xgn25 this isnt done no prior models'
                  ierr = 1
                  goto 500 
               ELSE
                  XPRIOR(1:inv%NA35) = DBLE(XMOD(1:inv%NA35))
               ENDIF
               F0 = inv%FOBJ
            ENDIF !end check on MYID
!
!.......... get new masters the source/receiver corrections
            IF (MYNID == MASTER) THEN
               IF (LUPDREC) THEN
                  IF (MYID == MASTER) &
                  WRITE(*,*) 'xgn25: Broadcasting receiver response functions...'
                  CALL BCAST_RCV_RESP(MYID,MYHD_COMM,MASTER, .FALSE.,frq%NFREQ,RCV)
               ENDIF
               IF (LUPDSRC) THEN
                  IF (MYID == MASTER) &
                  WRITE(*,*) 'xgn25: Broadcasting source time function...'
                  CALL BCAST_SRC_STF(MYID,MYHD_COMM,MASTER,  .FALSE.,frq%NFREQ,SRC) 
               ENDIF
            ENDIF
            IF (LUPDSRC) LUPDSRC = .FALSE.
            IF (LUPDREC) LUPDREC = .FALSE.
            CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
!----------------------------------------------------------------------------------------!
!                          Begin the Scalapack/BLACS factorization                       !
!----------------------------------------------------------------------------------------!
! 
!.......... have the head process send GRAD to the (0,0) grid 
            IF (MYID == MASTER) THEN
               WRITE(*,*) 'xgn25: Sending gradient to (0,0) grid process...'
               IF (MYID /= NBHEAD) THEN
                  IDEST = NBHEAD 
                  MYTAG = MYID 
                  CALL MPI_SEND(LSTF,1,MPI_LOGICAL, IDEST,MYTAG,      &
                                MPI_COMM_WORLD,MPIERR)
                  CALL MPI_SEND(GRAD,inv%NA35,MPI_REAL, IDEST,MYTAG,  &
                                MPI_COMM_WORLD,MPIERR)
                  CALL MPI_SEND(XMOD,inv%NA35,MPI_REAL, IDEST,MYTAG,  &
                                MPI_COMM_WORLD,MPIERR) 
                  CALL MPI_SEND(XPRIOR,inv%NA35,MPI_DOUBLE_PRECISION, IDEST,MYTAG, &
                                MPI_COMM_WORLD,MPIERR)
                  IF (LSTF) THEN
                     DO ISRC=1,src%NSRC
                        CALL MPI_SEND(src%SOURCE(:,ISRC),frq%NFREQ,MPI_COMPLEX,  &
                                      IDEST,MYTAG, MPI_COMM_WORLD,MPIERR)
                     ENDDO
                  ENDIF
               ENDIF 
            ELSE 
               IF (MYROW == 0 .AND. MYCOL == 0) THEN
                  CALL MPI_RECV(LSTF,1,MPI_LOGICAL, IDEST,MYTAG,      &
                                MPI_COMM_WORLD,MPIERR)
                  CALL MPI_RECV(GRAD,inv%NA35,MPI_REAL, MPI_ANY_TAG,MPI_ANY_SOURCE, &
                                MPI_COMM_WORLD,STAT,MPIERR)  
                  CALL MPI_RECV(XMOD,inv%NA35,MPI_REAL, MPI_ANY_TAG,MPI_ANY_SOURCE, &
                                MPI_COMM_WORLD,STAT,MPIERR)
                  CALL MPI_RECV(XPRIOR,inv%NA35,MPI_DOUBLE_PRECISION, MPI_ANY_TAG,  &
                                MPI_ANY_SOURCE, MPI_COMM_WORLD,STAT,MPIERR) 
                  IF (LSTF) THEN
                     IF (MYNID /= MASTER) ALLOCATE(src%SOURCE(frq%NFREQ,src%NSRC))
                     DO ISRC=1,src%NSRC
                        CALL MPI_RECV(src%SOURCE(:,ISRC),frq%NFREQ,MPI_COMPLEX, &
                                      MPI_ANY_TAG,MPI_ANY_SOURCE, &
                                      MPI_COMM_WORLD,STAT,MPIERR)
                     ENDDO
                  ENDIF 
               ENDIF
            ENDIF
            IF (MYROW ==-1 .OR. MYCOL ==-1) GOTO 2100
            IF (MYROW == 0 .AND. MYCOL == 0) THEN 
               ALLOCATE(DX(inv%NA35)) 
               DX(1:inv%NA35) = DBLE( XMOD(1:inv%NA35) - XPRIOR(1:inv%NA35) )
            ELSE
               ALLOCATE(DX(1)) 
            ENDIF
!
!.......... calculate the size of the local H matrices
            IF (MYROW == 0 .AND. MYCOL == 0) THEN
               WRITE(*,*) 'xgn25: Generating hloc...'
               NPPGRP = NPROCS/NPGROUPS
            ENDIF
            CALL GENHLOC(ICTXT,NPROW,NPCOL, MBJ,NBJ,MBJT,NBJT, DESCH, &
                         LDH,M,N, LSTF, NPPGRP,FRQ,RCV,INV,SRC, HLOC, IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xgn25: Error generating local H matrix!'
               GOTO 500 
            ENDIF
!
!.......... get the diagonal 
            ALLOCATE(DIAG(N)) 
            DIAG(1:N) = 0.D0 
            CALL GET_DIAG(ICTXT,NPROW,NPCOL, LDH,N, .TRUE.,DESCH,HLOC, DIAG,IERR) 
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xgn25: Error calling get_diag!' 
               GOTO 500 
            ENDIF
            IF (MYROW == 0 .AND. MYCOL == 0) THEN
               HESS(1:N) = SNGL(DIAG(1:N))
               ALLOCATE(GRAD8(N)) 
               DO I=1,N
                  GRAD8(I) = DBLE(GRAD(I))
               ENDDO
            ELSE
               ALLOCATE(GRAD8(1)) 
            ENDIF
            CALL TARANTOLA_394(ICTXT, LDH,LDC,NPROW,NPCOL, DESCH,DESCC, GRAD8,DX, &
                               HLOC,COVM, SEARCH8,IERR)  
            if (ierr /= 0) goto 500
            DEALLOCATE(DX) 
!
!.......... calculate the trust region (B + lambda I) p_1 =-g 
!           CALL NWALG_442(ICTXT,NPROW,NPCOL, K, LDH,DESCH,DELTA,HLOC,  GRAD8,  &
!                          RLAM,SEARCH8, IERR)  
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xng25: Error calling nwalg_442!'
               GOTO 500
            ENDIF
            IF (MYROW == 0 .AND. MYCOL == 0) SEARCH(1:N) = SNGL(SEARCH8(1:N))  
            DEALLOCATE(GRAD8)
 2100       CONTINUE !break ahead not in process grid
            CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
!
!.......... send solution and diagonal back to master
            IF (MYROW == 0 .AND. MYCOL == 0) THEN
               MYDEST = MASTER
               MYTAG  = MYID 
               IF (MYID /= NBHEAD) THEN 
                  CALL MPI_SEND(SEARCH,inv%NA35,MPI_REAL, IDEST,MYTAG, &
                                MPI_COMM_WORLD,MPIERR)
                  CALL MPI_SEND(HESS  ,inv%NA35,MPI_REAL, IDEST,MYTAG, &
                                MPI_COMM_WORLD,MPIERR)
               ENDIF
            ELSE
               IF (MYID == MASTER) THEN 
                  CALL MPI_RECV(SEARCH,inv%NA35,MPI_REAL, MPI_ANY_TAG,MPI_ANY_SOURCE,  &
                                MPI_COMM_WORLD,STAT,MPIERR)
                  CALL MPI_RECV(HESS  ,inv%NA35,MPI_REAL, MPI_ANY_TAG,MPI_ANY_SOURCE,  &
                                MPI_COMM_WORLD,STAT,MPIERR) 
               ENDIF
            ENDIF
!----------------------------------------------------------------------------------------!
!                           End the Scalapack/BLACS factorization                        !
!----------------------------------------------------------------------------------------!
            IF (MYID == MASTER) THEN !set model and plot
               P8(1:inv%NA35) = DBLE(SEARCH(1:inv%NA35))
               CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,               &
                                    inv%CINVTYPE, 1,IBLOCK,K,0,                     &
                                    msh%CDOMAIN,inv%MASKG,msh%IENG,                 &
                                    msh%XLOCS,msh%ZLOCS, inv%GRAD)
               CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,               &
                                    inv%CINVTYPE, 2,IBLOCK,K,0,                     &
                                    msh%CDOMAIN,inv%MASKG,msh%IENG,                 &
                                    msh%XLOCS,msh%ZLOCS, HESS) 
               CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,               &
                                    inv%CINVTYPE, 3,IBLOCK,K,0,                     &
                                    msh%CDOMAIN,inv%MASKG,msh%IENG,                 &
                                    msh%XLOCS,msh%ZLOCS, SEARCH)
            ENDIF !end check on myid
            inv%LFUNC  = .FALSE. !toggle off function evaluation
            inv%LGRAD  = .FALSE. !toggle off gradient evaluation
            inv%LGNEWT = .FALSE. !toggle off Jacobian evaluation
            LPGNEW = .FALSE.     !this should never be on 
            CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR) 
            goto 500
!----------------------------------------------------------------------------------------!
!                          Approximate higher order scattering                           !
!----------------------------------------------------------------------------------------!
            if (.not.allocated(p1)) then
               if (myid == master) then
                  allocate(p1(inv%na35))
               else
                  allocate(p1(1))
               endif
            endif
            if (.not.allocated(p2)) then 
               if (myid == master) then
                  allocate(p2(inv%na35))
               else
                  allocate(p2(1))
               endif
            endif
            if (myid == master) call scopy(inv%na35,search,1,p1,1) 
!           CALL QPMFREE(MASTER,MYID,MYNID,IPGROUP,NPGROUPS,    &
!                        MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM,   &
!                        m1d%VFAST,P1, INV,MSH,SRC,RCV,FRQ, MID, P2,IERR)
            IF (IERR /= 0) THEN
               IF (MYROW >= 0 .AND. MYCOL >= 0) THEN
                  diag(1:n) = 1.d0 
                  CALL ADD_DIAG(ICTXT,NPROW,NPCOL, LDH,N, DESCH, RLAM,DIAG, HLOC)  
                  print *, 'factoring matrix'
                  if (myrow == 0 .and. mycol == 0) diag(:) = dble(p2(:))
                  CALL SET_SOL8(ICTXT,NPROW,NPCOL, M, DESCB, DIAG, B) !set B to -g
                  CALL PDPOTRF('L',N, HLOC,1,1, DESCH, INFO)
                  print *, 'solving'
                  CALL PDTRSM('L', 'L', 'N', 'N', N, 1, & 
                              1.D0,HLOC, 1,1, DESCH, B, 1, 1, DESCB) !solve R^T q = p
                  CALL GET_SOL8(ICTXT,NPROW,NPCOL, N, DESCB, B, DIAG) !send q back to (0,0)
               ENDIF
            ENDIF
            IF (MYID == MASTER) THEN
               CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,               &    
                                    inv%CINVTYPE, 5,IBLOCK,K,0,                     &    
                                    msh%CDOMAIN,inv%MASKG,msh%IENG,                 &
                                    msh%XLOCS,msh%ZLOCS, P2)
               do i=1,inv%na35
                  search(i) =-inv%grad(i)/hess(i)
               enddo
               CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,               &    
                                    inv%CINVTYPE, 3,IBLOCK,K+1,0,                     &    
                                    msh%CDOMAIN,inv%MASKG,msh%IENG,                 &
                                    msh%XLOCS,msh%ZLOCS, SEARCH)
               search(:) = sngl(diag(:))
               call plot_shgrad_vtk(projnm,ngnod,msh%nnpg, msh%nelem, &
                                    inv%cinvtype, 5,iblock,k+2,0, &
                                    msh%cdomain,inv%maskg,msh%ieng, &
                                    msh%xlocs,msh%zlocs, search)
            ENDIF
            goto 500
!----------------------------------------------------------------------------------------!
!                              This is the line search loop                              !
!----------------------------------------------------------------------------------------!
            ICALPHA = 0 
            NFEV = 0 
            IF (MYID == MASTER) THEN
               STP = 1.D0/MAX(1.D0,MAXVAL(ABS(P8)))*DS0
               print *, 'stp0',stp
            ENDIF 
 3000       CONTINUE !loop on line search
               IF (inv%LFUNC .OR. inv%LGRAD) THEN 
                  CALL FUNCGRADH25(MASTER,MYID,MYNID,IPGROUP,NPGROUPS,                &
                                   MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM, LPGNEW,       &
                                   m1d%VFAST, INV,MSH,SRC,RCV,FRQ, MID, IERR) 
                  IF (IERR /= 0) THEN 
                     WRITE(*,*) 'xgn25: An error occurred calling FUNCGRADH25'
                     WRITE(*,*) 'xgn25: The error was on process:',MYID
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
                     WRITE(*,*) 'xgn25: Improper input parameters'
                     IERR = 1  
                     GOTO 500
                  ENDIF
                  IF (INFO == 2) THEN 
                     WRITE(*,*) 'xgn25: Relative width of uncertainty too small'
                     IERR = 1
                     GOTO 500
                  ENDIF
                  IF (INFO == 3) THEN 
                     WRITE(*,*) 'xgn25: More than MAXFEV function evaluations required', &
                                MAXFEV
                     IERR = 1
                     GOTO 500
                  ENDIF   
                  IF (INFO == 4) THEN
                     WRITE(*,*) 'xgn25: The step size is too small'
                     IERR = 1
                     GOTO 500
                  ENDIF   
                  IF (INFO == 5) THEN
                     WRITE(*,*) 'xgn25: The step size is too large'
                     IERR = 1
                     GOTO 500
                  ENDIF
                  IF (INFO == 6) THEN
                     WRITE(*,*) 'xgn25: Rounding errors prevent further progress'
                     IERR = 1
                     GOTO 500
                  ENDIF
!
!................ convergence test
                  IF (INFO == 1) THEN
                     ICALPHA = 0
                     MYDEST = 2000 !default, continue with iterative loop 
                     WRITE(*,8151) K,NFUN,STP,inv%FOBJ,100. - inv%FOBJ/F0*100.0
 8151                FORMAT(' ++++++++++++++++++++++++++++++++++++++++++++',/, &
                            ' + xgn25: Iteration:',I3,' finished            +',/,&
                            ' +        Function evaluations:',I3,'          +',/,&
                            ' +        Step length:',E12.4,'          +',/, &
                            ' +        Objective function:',E12.4,'   +',/,&
                            ' +        Percent reduction:',F8.3,'        +',/,&
                            ' +++++++++++++++++++++++++++++++++++++++++++',/)
                     IF (inv%FOBJ/F0*100. < 20.0) THEN
                        WRITE(*,*) 'xgn25: Convergence achieved'
                        MYDEST = 1005
                     ENDIF
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
 8150                FORMAT(' ++++++++++++++++++++++++++++++++++++++++++++',/, &
                            ' +  xgn25: Step length itereration:',I3,'       +',/,&
                            ' +         Function evaluations:',I3,'          +',/,&
                            ' +         Step length:',E12.4,'          +',/,&
                            ' +         Objective function:',E12.4,'   +',/,&
                            ' +++++++++++++++++++++++++++++++++++++++++++',/)
                  ENDIF
!
!................ slipped through the cracks? impossible, would be a memory error
                  IF (MYDEST ==-1) THEN
                     WRITE(*,*) 'xgn25: Error could not determine MYDEST',MYDEST
                     IERR = 1
                     GOTO 500
                  ENDIF
!
!................ update the models
                  WRITE(*,*) 'xgn25: Updating models...'
                  print *, lpgnew, inv%lgrad,inv%lfunc
                  CALL UPDMOD(PROJNM,LFILES, IBLOCK,K,ICALPHA, XMOD, INV,MSH, IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xgn25: Error updating model!'
                     GOTO 500
                  ENDIF
                  CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,                 &
                                       inv%CINVTYPE, 1,IBLOCK,K,ICALPHA,                 &
                                       msh%CDOMAIN,inv%MASKG,msh%IENG,                   &
                                       msh%XLOCS,msh%ZLOCS, inv%GRAD)

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
            IF (MYDEST == 1005) GOTO 1005 !we've finished
            IF (MYDEST == 2005) GOTO 2005 !we've hit the iteraiton limit
!----------------------------------------------------------------------------------------!
!                       The iteration has updated, broadcast source                      !
!----------------------------------------------------------------------------------------!
 2000    CONTINUE !iterative loop
 2005    CONTINUE !we have reached the max number of iterations
         IF (MYID == MASTER) &
         WRITE(*,*) 'xgn25: Iteration limit reached!'
 1005    CONTINUE !this block as converged
!
!....... and while we have them, get the latest STF and RRF estimates
         IF (MYID == MASTER .AND. inv%IMODSRC > 0) THEN
            WRITE(*,*) 'xgn25: Updating source time functions...'
            CALL SRCUPD(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                        inv%LUNWRAP,inv%IMODSRC, msh%AZMOD, inv%EST,inv%OBS, src%SOURCE)
            WRITE(*,*) 'xgn25: Convolving new STFs...'
            CALL CONVEST_STF(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                             src%SOURCE, inv%EST)
            WRITE(*,*) 'xgn25: Writing new source time functions...'
            CALL WTSTF(PROJNM, frq%NFREQ, frq%NFREQ,src%NSRC, IBLOCK,K, &
                       frq%FREQ,src%SOURCE, IERR)
         ENDIF
         IF (MYID == MASTER .AND. inv%IMODREC > 0) THEN
            WRITE(*,*) 'xgn25: Updating receiver response functions...'
            CALL RECUPD(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                        inv%LUNWRAP,inv%IMODREC, msh%AZMOD, inv%EST,inv%OBS, rcv%RECV)
            WRITE(*,*) 'xgn25: Convolving new RRFs...'
            CALL CONVEST_RRF(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
                             rcv%RECV, inv%EST)
            WRITE(*,*) 'xgn25: Writing new RRFs...'
            CALL WTRECST(PROJNM, NDIM,frq%NFREQ, NDIM,frq%NFREQ,rcv%NREC,  &
                         IBLOCK,K, msh%AZMOD, frq%FREQ,rcv%RECV, IERR)
         ENDIF
!
!....... write the final estimates
         IF (MYID == MASTER) THEN
            WRITE(*,*) 'xgn25: Writing estimates...'
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
         IF (ASSOCIATED(frq%CFTYPE)) DEALLOCATE(frq%CFTYPE) 
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
      IF (MYID == MASTER) WRITE(*,*) 'xgn25: Freeing memory...'
!.... model
      IF (ASSOCIATED(msh%ECOEFF))    DEALLOCATE(msh%ECOEFF)
      IF (ASSOCIATED(msh%DENS))      DEALLOCATE(msh%DENS)
      IF (ASSOCIATED(msh%QP))        DEALLOCATE(msh%QP)
      IF (ASSOCIATED(msh%QS))        DEALLOCATE(msh%QS) 
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
      IF (ASSOCIATED(src%SRCTYP))  DEALLOCATE(src%SRCTYP)
      IF (ASSOCIATED(src%PYTAB))   DEALLOCATE(src%PYTAB) 
      IF (MYNID == MASTER) THEN
         IF (ASSOCIATED(src%SOURCE))  DEALLOCATE(src%SOURCE)
         IF (ASSOCIATED(src%SRCTYP))  DEALLOCATE(src%SRCTYP)
         IF (ASSOCIATED(src%BAZN))    DEALLOCATE(src%BAZN)
         IF (ASSOCIATED(src%AOI))     DEALLOCATE(src%AOI)
      ENDIF

      IF (MYNID == MASTER) THEN
         IF (MYID == MASTER) THEN
            IF (ALLOCATED(XMOD))   DEALLOCATE(XMOD) 
            IF (ALLOCATED(X8))     DEALLOCATE(X8)
            IF (ALLOCATED(G8))     DEALLOCATE(G8)
            IF (ALLOCATED(P8))     DEALLOCATE(P8)
            IF (ALLOCATED(WA))     DEALLOCATE(WA)
            IF (ALLOCATED(XREC))   DEALLOCATE(XREC)
         ENDIF
      ENDIF 
      IF (ALLOCATED(SEARCH)) DEALLOCATE(SEARCH)

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
!
!.... BLACS/ScaLapack stuff
      IF (ALLOCATED(DIAG)) DEALLOCATE(DIAG) 
      IF (ALLOCATED(HESS)) DEALLOCATE(HESS)
      IF (ALLOCATED(HLOC)) DEALLOCATE(HLOC)
      IF (ALLOCATED(B))    DEALLOCATE(B) 
      mid%JOB =-2
      CALL CMUMPS(MID)
!
!.... release the BLACS grid and destroy BLACS data structures
      IF (LBLACS) THEN 
         IF (MYROW >= 0 .AND. MYCOL >= 0) CALL BLACS_GRIDEXIT(ICTXT)
         CALL BLACS_EXIT(1) !still more MPI calls so wait 
      ENDIF
!
!.... all done 
      TESIM = MPI_WTIME()
      IF (MYID == MASTER) THEN
         WRITE(*,9405) (TESIM - TSSIM)/3600.D0
 9405    FORMAT(' xgn25: Inversion time:',F8.2,' hours')
      ENDIF
      CALL MPI_FINALIZE(MPIERR)
      STOP
      END 
