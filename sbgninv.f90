
!
!     Surface/body wave joint inversion using the Gauss-Newton approximation
      IMPLICIT NONE
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
      TYPE (WIN_INFO) WIN
!.... project name
      CHARACTER(80) PROJNM, PROJNM_SRF, PROJNM_BDY, TMPDIR
!.... receiver information
      REAL*8, ALLOCATABLE :: XREC(:), ZREC(:)
      LOGICAL*4, ALLOCATABLE :: LFSURF(:) 
!.... MPI stuff
      REAL*8 TSSIM, TESIM
      INTEGER*4 STAT(MPI_STATUS_SIZE), MASTER,MYID,MYNID, &
                MYSLV_COMM,MYHD_COMM,IPGROUP,MPIERR, &
                NPPGRP, NPGROUPS, NPROCS, NPARTS, MYDEST, MYTAG, &
                IDEST, IPROCS
      PARAMETER(MASTER = 0)
!.... model gradient search direction stuff
      REAL*8, ALLOCATABLE :: XMOD8(:)      !double precision model
      REAL*8, ALLOCATABLE :: XMOD_SRF8(:)  !double precision surface wave model 
      REAL*8, ALLOCATABLE :: XMOD_BDY8(:)  !double precision body wave model
      REAL*8, ALLOCATABLE :: GRAD8(:)      !double precision surface wave + body wave gradient
      REAL*8, ALLOCATABLE :: GRAD_SRF8(:)  !double precision surface wave gradient
      REAL*8, ALLOCATABLE :: GRAD_BDY8(:)  !double precision body wave gradient 
      REAL*8, ALLOCATABLE :: SRCH8(:)      !double precision search direction
      REAL*8, ALLOCATABLE :: SRCH_SRF8(:)  !double precision surface wave search direction
      REAL*8, ALLOCATABLE :: SRCH_BDY8(:)  !double precision body wave search direction
      REAL*8, ALLOCATABLE :: DX(:)         !(x - x_prior)  
      REAL*8, ALLOCATABLE :: XSAVE(:)      !saves model
      REAL*4, ALLOCATABLE :: SEARCH(:)     !search direction for file IO 
!.... possibly quantify higher order scattering
      REAL*4, ALLOCATABLE :: P1(:)         !search direction; body + surface waves
      REAL*4, ALLOCATABLE :: PSSRF(:)      !HOT hessian search direction for surface waves
      REAL*4, ALLOCATABLE :: PSBDY(:)      !HOT hessian search direction for body waves
 
!.... BLACS stuff
      REAL*8, ALLOCATABLE :: HLOC_SRF(:,:) !local Re{ adj(J) J } matrix, surface waves
      REAL*8, ALLOCATABLE :: HLOC_BDY(:,:) !local Re{ adj(J) J } matrix, body waves
      REAL*8, ALLOCATABLE :: COVM_SRF(:,:) !Surface wave covariance matrix
      REAL*8, ALLOCATABLE :: COVM_BDY(:,:) !Body wave covariance matrix
      !REAL*8 RLAM          !for padding diagonal
      REAL*8 VPVAR_SRF     !variance for Vp surface waves
      REAL*8 VPVAR_BDY     !variance for Vp body waves
      REAL*8 SIGMA_SRFX    !surf wave characteristic distance in Gaussian kernel in x
      REAL*8 SIGMA_SRFZ    !surf wave characteristic distance in Guussian kernel in z
      REAL*8 SIGMA_BDYX    !body wave characteristic distance in Gaussian kernel in x
      REAL*8 SIGMA_BDYZ    !body wave characteristic distance in Gaussian kernel in x
      INTEGER*4 KERNEL     !
      INTEGER*4 DESCC(9)   !BLACS descriptor for covariance matrices
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
      INTEGER*4 LDB, LDH, &!leading dimensions for B, HLOC, and covariance
                LDC      
      LOGICAL*4 LBLACS     !True -> BLACS has been initialized
      LOGICAL*4 LPRIOR     !True -> using a priod model in updates
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &    
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 ) 
!.... miscellaneous
      LOGICAL*4, ALLOCATABLE :: LFSAVE(:) !saves original copy of inversion freqs
      COMPLEX*8 CZERO
      REAL*8 XAVG, GET_STP0, VPMAX_INV, VPMIN_INV, FMIN
      REAL*4 F0_BDY, F0_SRF, F0, ALPHA
      INTEGER*4 NSRC_BDY_HD, NSRC_SRF_HD, NFREQ_SRF_HD, NFREQ_BDY_HD,  &
                NREC_HD, NDIM_HD, NELEME, NABS, IOMINV, &
                IBLOCK, NBLOCKS, K, NFUN, NFUN_SRF, NFUN_BDY, NFEV, NFEV_SRF, NFEV_BDY,  &
                ICALPHA, ICALPHA_SRF, ICALPHA_BDY, &
                ISRC, ICGRNS, ISTOP_PT, IERR  
      LOGICAL*4 LINREC, LFILES, LSRCEX, LRECEX, LUPDSRC, LUPDREC, LFUNC,  &
                LGRAD, LINIT, LSTF, LCGRNS, LEXS, LEXB, LKBDRY, LINV_BDY, &
                LINV_SRF, LNSRF
      PARAMETER(CZERO = CMPLX(0.0,0.0))
!.... line search
      REAL*8 STPMIN, STPMAX, FSRF8, FOBJ8, FBDY8, GTOL8, XTOL8, FTOL8, STP, STP_SRF, &
             STP_BDY, DS0, DP0, STPMIN_SRF, STPMIN_BDY, STP0
      REAL*8, ALLOCATABLE :: WA(:) !workspace for line search
      INTEGER*4 MAXFEV, LP, INFO 
      PARAMETER(LP = 6) 
!.... functions
      INTEGER*4 ICELEML
      logical*4 ljsrch
      !real*8 delta
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
!.... head node information
      IERR = 0
      IF (MYID == MASTER) THEN 
         PROJNM(1:80) = ' '  
         PROJNM_SRF(1:80) = ' '
         PROJNM_BDY(1:80) = ' '
         WRITE(*,8000) 
 8000    FORMAT(' ------------------------------------------------------',/, &
                ' -  xjoint: A massively parallel 2.5D unstructured    -',/, &
                ' -          grid joint inversion method for elastic   -',/, &
                ' -          parameters using the Gauss-Newton method  -',/, &
                ' ------------------------------------------------------',/)
         WRITE(*,*) 'xjoint: Enter project name:' 
         READ(*,'(A)') PROJNM
         PROJNM = ADJUSTL(PROJNM)
         PROJNM_SRF = TRIM(PROJNM)//'_surf'
         PROJNM_BDY = TRIM(PROJNM)//'_body'
         PROJNM_SRF = ADJUSTL(PROJNM_SRF)
         PROJNM_BDY = ADJUSTL(PROJNM_BDY)
         WRITE(*,*) 'xjoint: Enter number of process groups:' 
         READ *, NPGROUPS
         WRITE(*,*) 'xbielak25: Enter 1 to load the Greens functions from disk'
         WRITE(*,*) '           Any other number will calculate Greens functions'
         READ *, ICGRNS
         IF (ICGRNS == 1) THEN
            LCGRNS = .FALSE.
         ELSE
            LCGRNS = .TRUE.
         ENDIF
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
         CALL READ_INV_INI(PROJNM,.TRUE., &
                           LPRIOR,LJSRCH,LINREC,LKBDRY, &
                           ISTOP_PT,MAXFEV,KERNEL,                        &   
                           GTOL8,XTOL8,FTOL8, STPMIN,STPMAX, DS0,DP0,     &   
                           SIGMA_BDYX,SIGMA_BDYZ, SIGMA_SRFX,SIGMA_SRFZ,  &
                           VPVAR_BDY,VPVAR_SRF, VPMAX_INV,VPMIN_INV,      &
                           INV,IERR) 
         !delta = dsqrt( 2000.d0*50.d0**2 ) !relic
!....... lazy error
         IF (inv%CINVTYPE == 'PS' .OR. inv%CINVTYPE == 'ps') then
            WRITE(*,*) 'xjoint: Cannot perform simulataneous P/S inversion yet!'
            IERR = 1
            GOTO 500
         ENDIF
!....... initial warnings on NPGROUPS
         IF (MOD(NPROCS,NPGROUPS) /= 0) THEN
            WRITE(*,*) 'xjoint: Error I cant divide nprocs by npgroups evenly!'
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
            WRITE(*,*) 'xjoint: I do not know what to invert for:',inv%CINVTYPE
            IERR = 1
            GOTO 500
         ENDIF
!
!....... read the mesh
         WRITE(*,*) 'xjoint: Reading mesh...'
         CALL RDMESHBK(PROJNM, NELEME,NABS, MSH, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xjoint: Error reading mesh!'
            IERR = 1
            GOTO 500
         ENDIF
         ALLOCATE(msh%XIPTS (msh%NLXI))
         ALLOCATE(msh%ETAPTS(msh%NLETA))
         CALL GENPTS(msh%IITYPE, msh%NLXI,msh%NLETA, msh%XIPTS,msh%ETAPTS)
!
!....... initial check on observation file 
         WRITE(*,*) 'xjoint: Checking headers on observation file...'
         !CALL RDTOBS_HD(PROJNM, NDIM_HD,NFREQ_HD,NREC_HD,NSRC_HD, IERR)
         CALL RDTOBS_SB_HD(PROJNM, NDIM_HD,NFREQ_SRF_HD,NFREQ_BDY_HD,  & 
                           NREC_HD,NSRC_SRF_HD,NSRC_BDY_HD, IERR)
         IF (IERR /= 0) THEN 
            WRITE(*,*) 'xjoint: There was an error with your observation file'
            GOTO 500  
         ENDIF 
         IF (NDIM_HD /= NDIM) THEN 
            WRITE(*,*) 'xjoint: Error component number mismatch!',NDIM_HD,NDIM
            IERR = 1  
            GOTO 500  
         ENDIF
         IF (NFREQ_SRF_HD + NFREQ_BDY_HD <= 0) THEN 
            WRITE(*,*) 'xjoint: Error observations have 0 frequencies!', &
                       NFREQ_SRF_HD, NFREQ_BDY_HD 
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
!
!....... read the source list
         CALL RDSRC_EQ(PROJNM, msh%XLATMIN,msh%XLONMIN,msh%XLATMAX,msh%XLONMAX, &
                       msh%AZMOD, SRC,IERR)
         IF (IERR /= 0) THEN 
            WRITE(*,*) 'xjoint: Error calling rdsrc_eq!'
            GOTO 500  
         ENDIF
         IF (NSRC_SRF_HD /= src%NSRC_SRF) THEN 
            WRITE(*,*) 'xjoint: Surface wave source number mismatch!', &
                       NSRC_SRF_HD,src%NSRC_SRF
            IERR = 1
            GOTO 500
         ENDIF
         IF (NSRC_BDY_HD /= src%NSRC_BDY) THEN 
            WRITE(*,*) 'xjoint: Body wave source number mismatch!', &
                       NSRC_BDY_HD,src%NSRC_BDY
            GOTO 500
         ENDIF
         inv%LSURF = .FALSE.
         inv%LBODY = .FALSE.
         CALL LINV_SB(PROJNM, 1, src%NSRC_SRF,src%NSRC_BDY,  &
                      inv%LSURF,inv%LBODY,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xjoint: Error calling linvb_sb!'
            GOTO 500
         ENDIF
         IF (.NOT.inv%LSURF .AND. .NOT.inv%LBODY) THEN
            WRITE(*,*) 'xjoint: We arent inverting anything!'
            IERR = 1
            GOTO 500
         ENDIF
         IF (inv%LSURF .AND. inv%LBODY) THEN
            WRITE(*,*) 'xjoint: This is a joint inversion' 
         ELSE
            IF (inv%LSURF) THEN
               WRITE(*,*) 'xjoint: This is a surface wave inverison' 
            ELSE
               WRITE(*,*) 'xjoint: This is a body wave inversion'
            ENDIF
         ENDIF
!
!....... check for a receiver file
         CALL RDREC_HD(PROJNM, rcv%NREC,IERR)
         IF (IERR /= 0) THEN
            IF (IERR > 0) THEN
               WRITE(*,*) 'xjoint: Error reading recv header'
               GOTO 500
            ELSE
               WRITE(*,*) 'xjoint: There is no receiver file!  Cannot make estimates!'
               IERR = 1
               GOTO 500
            ENDIF
         ENDIF
         IF (rcv%NREC <= 0) THEN
            WRITE(*,*) 'xjoint: There are no receivers!  Aborting!'
            IERR = 1
            GOTO 500
         ENDIF
         IF (NREC_HD /= rcv%NREC) THEN
            WRITE(*,*) 'xjoint: Error recevier number mismatch!',NREC_HD,rcv%NREC
            IERR = 1
            GOTO 500
         ENDIF
         ALLOCATE(  LFSURF(rcv%NREC))
         ALLOCATE(    XREC(rcv%NREC))
         ALLOCATE(rcv%YREC(rcv%NREC))
         ALLOCATE(    ZREC(rcv%NREC))
         CALL RDREC(PROJNM, rcv%NREC, LFSURF,XREC,rcv%YREC,ZREC, IERR)
         WRITE(*,*) 'xjoint: Splitting process groups...' 
         ALLOCATE(inv%DX(rcv%NREC-1))
         inv%DX(1:rcv%NREC-1) = 0.D0
         IF (inv%LCOVD_SRF .OR. inv%LCOVD_BDY) THEN
            CALL REC_SPACE1D(rcv%NREC,XREC, inv%DX,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xjoint: Error detected in rec_space1d'
               IF (inv%LCOVD_SRF) WRITE(*,*) 'xjoint: Setting lcovd_srf to false'
               IF (inv%LCOVD_BDY) WRITE(*,*) 'xjoint: Setting lcovd_bdy to false'
            ENDIF
         ENDIF
      ENDIF
      CALL MKMPI_GROUPS(MASTER,MYID,NPGROUPS,NPROCS,MPI_COMM_WORLD, &
                        IPGROUP,MYNID,MYHD_COMM,MYSLV_COMM)
      NPPGRP = NPROCS/NPGROUPS
!----------------------------------------------------------------------------------------!
!     MUMPS initialization phase and graph reordering utilities.  The mesh should not    !
!     have to change.  If it does the user should re-mesh then re-run the inversion      !
!     software.                                                                          !
!----------------------------------------------------------------------------------------!
      IF (MYID == MASTER) WRITE(*,*) 'xjoint: Initializing MUMPS...'
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
         WRITE(*,*) 'xjoint: Generating graph...'
         CALL GEN_GRAPH25(.TRUE.,NPARTS, LFSURF,XREC,ZREC, MSH,RCV,IERR) 
         IF (inv%LUNWRAP) THEN 
            WRITE(*,*) 'xjoint: Generating free surface information...'
            CALL GEN_GIDOFFS(MSH,IERR)
            IF (IERR /= 0) THEN 
               WRITE(*,*) 'xjoint: Error calling gen_gidoffs!'
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
 9408    FORMAT(' xjoint: Polynomial order:'                    ,I4 ,/,        &    
                '         Number of elements in mesh:'          ,I10,/,         &    
                '         Number of absorbing elements:'        ,I8 ,/,         &    
                '         Number of Bielak elements:'           ,I8 ,/,         &    
                '         Number of anchor nodes in mesh:'      ,I10,/,         &    
                '         Number of nodes in Bielak boundary:'  ,I10,/,         &    
                '         Number of degrees of freedom:'        ,I14,/,         &    
                '         Number of non-zeros in global matrix:',I16,/)
         IF (ALLOCATED(LFSURF)) DEALLOCATE(LFSURF) 
         IF (ALLOCATED(ZREC)) DEALLOCATE(ZREC)
!
!....... create the pointers associated w/ the inverse problem
         WRITE(*,*) 'xjoint: Generating mask for inverse problem...'
         CALL GENGMASK(PROJNM, NDIM,msh%NEN,NGNOD, msh%NDOF,msh%NNPG,msh%NELEM,     &    
                       msh%NLXI,msh%NLETA,LINREC,LKBDRY,                            &
                       msh%CDOMAIN,msh%CNNPG,msh%PART,  &
                       msh%LM,msh%IENG, RCV,INV, IERR) 
         IF (IERR /= 0) THEN 
            WRITE(*,*) 'xjoint: Error generating gradient mask!'
            GOTO 500
         ENDIF
         CALL ELEM_WTS(MSH,INV,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xjoint: Error calling ELEM_WTS!'
            GOTO 500
         ENDIF
         inv%NA35 = inv%NNPINV*inv%NVINV
         inv%MCJLOC = rcv%NREC*NDIM
         ALLOCATE(DX(inv%NA35))
         DX(1:inv%NA35) = 0.D0
         ALLOCATE(XSAVE(inv%NA35))
         XSAVE(:) = 0.D0
         IF (inv%LSURF) THEN
            ALLOCATE(XMOD_SRF8(inv%NA35))
            ALLOCATE(GRAD_SRF8(inv%NA35))
            ALLOCATE(SRCH_SRF8(inv%NA35))
            ALLOCATE(inv%GRAD_SRF(inv%NA35))
            XMOD_SRF8(:) = 0.D0
            GRAD_SRF8(:) = 0.D0
            SRCH_SRF8(:) = 0.D0
            inv%GRAD_SRF(:) = 0.0
            CALL GENXMOD8(msh%NNPG,inv%NA35,msh%NNPG, inv%CINVTYPE,    &
                         inv%MASKG,msh%DENS,msh%ECOEFF, XMOD_SRF8,IERR)
         ENDIF
         IF (inv%LBODY) THEN
            ALLOCATE(XMOD_BDY8(inv%NA35))
            ALLOCATE(GRAD_BDY8(inv%NA35))
            ALLOCATE(SRCH_BDY8(inv%NA35))
            ALLOCATE(inv%GRAD_BDY(inv%NA35))
            XMOD_BDY8(:) = 0.D0
            GRAD_BDY8(:) = 0.D0
            SRCH_BDY8(:) = 0.D0
            inv%GRAD_BDY(:) = 0.0
            CALL GENXMOD8(msh%NNPG,inv%NA35,msh%NNPG, inv%CINVTYPE,    &
                          inv%MASKG,msh%DENS,msh%ECOEFF, XMOD_BDY8,IERR)
         ENDIF
         IF (inv%LSURF .AND. inv%LBODY) THEN
            ALLOCATE(GRAD8(inv%NA35))
            ALLOCATE(XMOD8(inv%NA35))
            ALLOCATE(SRCH8(inv%NA35))
            ALLOCATE(inv%GRAD(inv%NA35))
            GRAD8(1:inv%NA35) = 0.D0
            XMOD8(1:inv%NA35) = 0.D0
            SRCH8(1:inv%NA35) = 0.D0
            inv%GRAD(:) = 0.0
       
            CALL GENXMOD8(msh%NNPG,inv%NA35,msh%NNPG, inv%CINVTYPE,    &
                          inv%MASKG,msh%DENS,msh%ECOEFF, XMOD8,IERR)
         ENDIF
         ALLOCATE(WA(inv%NA35)) !workspace for line search 
!        ALLOCATE(XMOD(inv%NA35))
!        ALLOCATE(SEARCH(inv%NA35)) !search direction, solution of Gauss-Newton system
!        ALLOCATE(WA(inv%NA35))
!        ALLOCATE(X8(inv%NA35))
!        ALLOCATE(XPRIOR(inv%NA35))
!        ALLOCATE(G8(inv%NA35))
!        ALLOCATE(P8(inv%NA35))
!        ALLOCATE(GRAD(inv%NA35))
!        ALLOCATE(HESS(inv%NA35))
      ELSE
!        ALLOCATE(inv%GRAD(1)) 
!        ALLOCATE(inv%GRAD_SRF(1))
!        ALLOCATE(inv%GRAD_BDY(1)) 
      ENDIF
!----------------------------------------------------------------------------------------!
!                          Broadcast Inverse Problem Parameters                          !
!----------------------------------------------------------------------------------------!
      CALL MPI_BCAST(LJSRCH ,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
      IF (MYID == MASTER) WRITE(*,*) 'xjoint: Broadcasting 2D model...' 
      CALL MPI_BCAST(NPARTS ,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
     
      CALL BCAST_MESH_INFO(MYID,MPI_COMM_WORLD,MASTER, &
                           LNSRF, PROJNM,TMPDIR, MSH) 
      !IF (MYNID == MASTER) THEN 
      IF (MYID == MASTER) WRITE(*,*) 'xjoint: Broadcasting receiver locations...'
      CALL BCAST_RCV_INFO(MYID,MPI_COMM_WORLD,MASTER, RCV) 
      !ENDIF
      IF (MYID == MASTER) WRITE(*,*) 'xjoint: Broadcasting inversion parameters...'
      CALL MPI_BCAST(rcv%NREC,1,MPI_INTEGER, MASTER,MYSLV_COMM,MPIERR) !for jacobian
      CALL BCAST_INV_PARMS(MYID,MPI_COMM_WORLD,MASTER, msh%NNPG,msh%NELEM,NBLOCKS, &
                           rcv%NREC, ISTOP_PT, INV) 
      IF (MYNID == MASTER .AND. inv%LUNWRAP) THEN 
         IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Broadcasting free surface DOFS...'
         CALL BCAST_FS_INFO(MYID,MYHD_COMM,MASTER, MSH) 
      ENDIF
!
!.... figure out the local gradient graph
      IF (MYID == MASTER) WRITE(*,*) 'xjoint: Generating gradient graph...'
      CALL GRADPTRS(MYNID,MASTER,MYSLV_COMM, NDIM,msh%NEN, msh%NEN,msh%NNPG, msh%LM, &
                    INV,IERR)
      IF (IERR /= 0) THEN 
         IF (MYNID == MASTER) WRITE(*,*) 'xjoint: An error occurred in gradptrs!'
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
      IF (MYID == MASTER) WRITE(*,*) 'xjoint: Generating local graphs...'
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
         IF (NPROCS > 3) THEN
            CALL GEN_PGRID_SQR(N,N,NPROCS, NPROW,NPCOL,IERR)
         ELSE
           CALL GEN_PGRID(N,N,NPROCS, NPROW,NPCOL,IERR)
         ENDIF 
         IF (IERR /= 0) THEN 
            WRITE(*,*) 'xjoint: Error calling gen_pgrid!'
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
         LDC = LDH
         DESCC(1:9) = DESCH(1:9)
         IF (inv%LSURF) THEN 
            ALLOCATE(HLOC_SRF(LDH,MAX(NH,1)))
            HLOC_SRF(:,:) = 0.D0
            ALLOCATE(COVM_SRF(LDC,MAX(NH,1)))
            COVM_SRF(:,:) = 0.D0
         ENDIF
         IF (inv%LBODY) THEN
            ALLOCATE(HLOC_BDY(LDH,MAX(NH,1)))
            HLOC_BDY(:,:) = 0.D0
            ALLOCATE(COVM_BDY(LDC,MAX(NH,1)))
            COVM_BDY(:,:) = 0.D0
         ENDIF
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
!        IF (.NOT.ALLOCATED(GRAD_BDY))  ALLOCATE(GRAD_SRF(inv%NA35))
!        IF (.NOT.ALLOCATED(GRAD_BDY))  ALLOCATE(GRAD_BDY(inv%NA35))
!        IF (.NOT.ALLOCATED(SRCH_SRF))  ALLOCATE(SRCH_SRF(inv%NA35))
!        IF (.NOT.ALLOCATED(SRCH_BDY))  ALLOCATE(SRCH_BDY(inv%NA35))
!        IF (.NOT.ALLOCATED(SEARCH8)) ALLOCATE(SEARCH8(inv%NA35))
!        IF (.NOT.ALLOCATED(HESS))    ALLOCATE(HESS(inv%NA35))
!        IF (.NOT.ALLOCATED(XMOD))    ALLOCATE(XMOD(inv%NA35))
!        IF (.NOT.ALLOCATED(X8))      ALLOCATE(X8(inv%NA35))
!        IF (.NOT.ALLOCATED(XPRIOR))  ALLOCATE(XPRIOR(inv%NA35))
      ELSE
!        IF (.NOT.ALLOCATED(GRAD))    ALLOCATE(GRAD(1))
!        IF (.NOT.ALLOCATED(SEARCH8)) ALLOCATE(SEARCH8(1))
!        IF (.NOT.ALLOCATED(SEARCH))  ALLOCATE(SEARCH(1))
!        IF (.NOT.ALLOCATED(XMOD))    ALLOCATE(XMOD(1))
!        IF (.NOT.ALLOCATED(X8))      ALLOCATE(X8(1))
!        IF (.NOT.ALLOCATED(XPRIOR))  ALLOCATE(XPRIOR(1))
      ENDIF
      IF (ISTOP_PT == 1) THEN !stop after initializaiton
         IF (MYID == MASTER) WRITE(*,*) 'xjoint: Terminating at stop_pt 1'
         GOTO 500 
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
            WRITE(*,*) 'xjoint: Filling 1D models..'
            CALL FILL1D(MSH,M1D) 
         ENDIF
         CALL MPI_BCAST(msh%XMOD0  ,1,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,MPIERR) 
         CALL MPI_BCAST(msh%XMOD1  ,1,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,MPIERR) 
         CALL MPI_BCAST(msh%XBLKL  ,1,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(msh%XBLKR  ,1,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,MPIERR)
         IF (MYID == MASTER) THEN
            WRITE(*,*) 'xjoint: Reading windows...'
            ALLOCATE(win%TWIN(rcv%NREC,src%NSRC,2)) 
            CALL RDWIN(PROJNM, rcv%NREC,src%NSRC, src%SRCTYP,  WIN,IERR) 
            IF (IERR /= 0) THEN
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xjoint: Window file does not exist'
               ELSE
                  WRITE(*,*) 'xjoint: There was an error reading the window file'
                  GOTO 500
               ENDIF
            ENDIF
            WRITE(*,*)
            WRITE(*,*) 'xjoint: Reading frequency block...'
            CALL RD_JFREQ_INV(PROJNM,IBLOCK,win%LWNDO_SRF,win%LWNDO_BDY,  &
                              FRQ,IERR)
            IF (IERR /= 0) THEN 
               WRITE(*,*) 'xjoint: Error calling rd_jfreq_inv!'
               GOTO 500
            ENDIF
            ALLOCATE(LFSAVE(frq%NFREQ)) 
!           IF (win%LWNDO_SRF) THEN
!              NSPACE = frq%NFREQ_SRF
!           ELSE
!              NSPACE = frq%NFREQ_SRF_INV
!           ENDIF 
!           IF (win%LWNDO_BDY) THEN
!              NSPACE = NSPACE + frq%NFREQ_BDY 
!           ELSE
!              NSPACE = NSPACE + frq%NFREQ_BDY_INV
!           ENDIF
!           NFREQ_OBS = NSPACE
            WRITE(*,*) 'xjoint: Reading observation file...'
            print *, inv%dx 
            ALLOCATE(inv%  EST(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
            ALLOCATE(inv%  OBS(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
            ALLOCATE(inv%WGHTS(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
            inv%OBS(:,:,:,:) = CZERO
            inv%EST(:,:,:,:) = CZERO
            inv%WGHTS(:,:,:,:) = 0.0
            CALL RDTOBS_SB25(PROJNM, RCV,SRC,FRQ,INV, IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xjoint: Error reading observation file!'
               GOTO 500
            ENDIF
            CALL NULL_OBSFD(NDIM,frq%NFREQ,rcv%NREC,                                   &
                            NDIM,frq%NFREQ,rcv%NREC,src%NSRC, inv%NREC_MIN, .TRUE.,    &
                            inv%LCOVD_SRF,inv%LCOVD_BDY, frq%CFTYPE,src%SRCTYP,        &
                            inv%WGHTS,inv%OBS)
            WRITE(*,*) 'xjoint: Checking for srcinv file...'
            LUPDSRC = .FALSE. 
            ALLOCATE(src%SOURCE(frq%NFREQ,src%NSRC))
            CALL RDSRCINV_SB(PROJNM, frq%NFREQ,frq%NFREQ,src%NSRC,                   &
                             frq%NFREQ_SRF,frq%NFREQ_BDY, src%NSRC_SRF,src%NSRC_BDY, &
                             frq%CFTYPE,src%SRCTYP, frq%FREQ,  LEXS,LEXB,src%SOURCE,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xjoint: Error calling rdsrcinv_sb!'
               GOTO 500
            ENDIF
            LSRCEX = .TRUE.
            IF (inv%LSURF .AND. inv%LBODY) THEN
               IF (.NOT.LEXS .OR. .NOT.LEXB) LSRCEX = .FALSE.
            ELSE
               IF (inv%LSURF) THEN
                  IF (.NOT.LEXS) LSRCEX = .FALSE.
               ELSE
                  IF (.NOT.LEXB) LSRCEX = .FALSE.
               ENDIF
            ENDIF
            IF (.NOT.LSRCEX) THEN
               WRITE(*,*) 'xjoint: Attempting to read time domain sources...'
               CALL SRCSUB_SB(PROJNM,frq%NFREQ, frq%NFREQ,src%NSRC,                   &
                              frq%NFREQ_SRF,frq%NFREQ_BDY,src%NSRC_SRF,src%NSRC_BDY,  &
                              frq%CFTYPE,src%SRCTYP, frq%FREQ, src%SOURCE,IERR) 
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xjoint: Error calling srcsub_sb!'
                  LUPDSRC = .TRUE. !you will never succeed without diong this
                  IF (inv%IMODSRC == 0) THEN 
                     LUPDSRC = .FALSE. !fine, you asked for it
                     WRITE(*,*) 'xjoint: Warning youre proceeding without an STF estimate'
                  ENDIF
                  IERR = 0
               ENDIF
            ELSE 
               WRITE(*,*) 'xjoint: Source inversion file read'
            ENDIF
!           IF (inv%IMODSRC > 0) THEN
!              LGRNS = .TRUE.
!              DO IFREQ=1,frq%NFREQ 
!                 DO ISRC=1,src%NSRC
!                    IF (src%SOURCE(IFREQ,ISRC) /= CMPLX(1.0,0.0)) LGRNS = .FALSE.
!                 ENDDO
!              ENDDO
!           ENDIF
      !     IF (inv%IRESTP == 1) THEN 
      !        WRITE(*,*) 'xjoint: Normalizing source scale to unity...'
      !     ELSEIF (inv%IRESTP == 2) THEN 
      !        WRITE(*,*) 'xjoint: Setting source phase to 0...'
      !     ENDIF 
      !     CALL MODSRC(frq%NFREQ, frq%NFREQ,src%NSRC,inv%IRESTP, src%SOURCE) 
            WRITE(*,*) 'xjoint: Checking for wrec file'
            LUPDREC = .FALSE.
            ALLOCATE(rcv%RECV(NDIM,frq%NFREQ,rcv%NREC)) 
            CALL RDREC_RESP(PROJNM, NDIM,frq%NFREQ, NDIM,frq%NFREQ,rcv%NREC, msh%AZMOD, &
                            frq%FREQ, LRECEX,rcv%RECV)
            IF (.NOT.LRECEX) THEN 
               LUPDREC = .TRUE.
               IF (inv%IMODREC == 0) THEN 
                  LUPDREC = .FALSE. !not as much of a sin as forgetting the STF
                  WRITE(*,*) 'xjoint: Note you are proceeding without an RRF estimate'
               ENDIF 
               IERR = 0  
            ELSE
               WRITE(*,*) 'xjoint: Receiver inversion file read'
            ENDIF
      !     IF (inv%IRESTP == 1) THEN
      !        WRITE(*,*) 'xjoint: Normalizing receiver responses to unity...'
      !     ELSEIF (inv%IRESTP == 2) THEN
      !        WRITE(*,*) 'xjoint: Setting receiver response phase to 0...'
      !     ENDIF
      !     CALL MODREC(NDIM,frq%NFREQ, NDIM,frq%NFREQ,rcv%NREC, inv%IRESTP, rcv%RECV)

            WRITE(*,8100) IBLOCK 
            DO 1101 IOMINV=1,frq%NFREQ 
               WRITE(*,8101) IOMINV,frq%CFTYPE(IOMINV),frq%LINVF(IOMINV),frq%FREQ(IOMINV) 
 1101       CONTINUE !loop on frequency blocks 
            WRITE(*,8102) 
 8100       FORMAT(/,' -------------------------------------------------',/, &
                     ' -   Frequencies in Block:',I4,'                   -',/, &
                     ' -                                               -')
 8101       FORMAT(  ' -   Frequency Number:',I3,1X,A1,1X,L1,1X,F12.5,' (Hz)  -')
 8102       FORMAT(  ' -                                               -',/, &
                     ' -------------------------------------------------',/)
         ENDIF
!
!....... broadcast frequency block and observations
         CALL MPI_BCAST(rcv%NREC ,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
!        CALL MPI_BCAST(NFREQ_OBS,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
         IF (MYID == MASTER) WRITE(*,*) 'xjoint: Broadcasting 1D models...'
         CALL BCAST_MOD1D_INFO(MYID,MPI_COMM_WORLD,MASTER, MSH,M1D)
         IF (MYID == MASTER) WRITE(*,*) 'xjoint: Broadcasting source details...'
         CALL BCAST_SRC_INFO(MYID,MPI_COMM_WORLD,MASTER, SRC)
         IF (MYID == MASTER) &
         WRITE(*,*) 'xjoint: Broadcasting inversion frequency information...'
         CALL BCAST_FRQ_JOINT_INFO(MYID,MPI_COMM_WORLD,MASTER, FRQ)
         IF (MYNID == MASTER) THEN
            IF (MYID == MASTER) WRITE(*,*) 'xjoint: Broadcasting window information...'
            CALL BCAST_WIN(MYID,MYHD_COMM,MASTER,rcv%NREC,src%NSRC,rcv%NREC,src%NSRC, WIN)
         ENDIF
         IF (MYID == MASTER) &
         WRITE(*,*) 'xjoint: Broadcasting observations...'
         CALL BCAST_OBS_INFO(MYID,MPI_COMM_WORLD,MASTER, frq%NFREQ,rcv%NREC,src%NSRC, &
                             INV)
         IF (MYNID == MASTER) THEN
            IF (MYID /= MASTER) ALLOCATE(inv%EST(NDIM, frq%NFREQ,rcv%NREC,src%NSRC))
            IF (MYID == MASTER) &
            WRITE(*,*) 'xjoint: Broadcasting receiver response functions...'
            CALL BCAST_RCV_RESP(MYID,MYHD_COMM,MASTER, .TRUE.,frq%NFREQ,RCV)
            IF (MYID == MASTER) &
            WRITE(*,*) 'xjoint: Broadcasting source time function...'
            CALL BCAST_SRC_STF(MYID,MYHD_COMM,MASTER,  .TRUE.,frq%NFREQ,SRC)
         ENDIF
!----------------------------------------------------------------------------------------!
!                            Calculate 1D models solutions                               !
!----------------------------------------------------------------------------------------!
         IF (MYID == MASTER) WRITE(*,*) 'xjoint: Calculating 1D solutions...'
         CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
         ALLOCATE(src%PYTAB(frq%NFREQ,src%NSRC))
         src%PYTAB(1:frq%NFREQ,1:src%NSRC) = 0.D0
         CALL GEN_GRNS_SB(MYID,MASTER,MPI_COMM_WORLD, TMPDIR,LCGRNS,LNSRF, &
                          SRC,M1D,FRQ,MSH, IERR)
         IF (IERR /= 0) THEN 
            WRITE(*,*) 'xjoint: Error in gen_grns!'
            GOTO 500
         ENDIF
!....... free source
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
!....... generate the covariance matrix
         CALL MPI_BCAST(KERNEL    ,1,MPI_INTEGER         ,MASTER,MPI_COMM_WORLD,MPIERR)
         !CALL MPI_BCAST(DELTA     ,1,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(SIGMA_SRFX,1,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(SIGMA_SRFZ,1,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(VPVAR_SRF ,1,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(SIGMA_BDYX,1,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(SIGMA_BDYZ,1,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(VPVAR_BDY ,1,MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,MPIERR)
         IF (MYROW >= 0 .AND. MYCOL >= 0) THEN
            IF (inv%LSURF) THEN
               IF (MYROW == 0 .AND. MYCOL == 0) &
               WRITE(*,*) 'xjoint: Generating covariance matrix for surface waves...'
               CALL COVMODEL_DIST(LDC, NPROW,NPCOL,ICTXT,                       &
                                  KERNEL,SIGMA_SRFX,SIGMA_SRFZ,VPVAR_SRF,DESCC, &
                                  INV,MSH, COVM_SRF,IERR)
               IF (IERR /= 0) GOTO 500
               print *, minval(covm_srf),maxval(covm_srf)
            ENDIF
            IF (inv%LBODY) THEN
               IF (MYROW == 0 .AND. MYCOL == 0) &
               WRITE(*,*) 'xjoint: Generating covariance matrix for body waves...'
               CALL COVMODEL_DIST(LDC, NPROW,NPCOL,ICTXT,                        &
                                  KERNEL,SIGMA_BDYX,SIGMA_BDYZ,VPVAR_BDY,DESCC, &
                                  INV,MSH, COVM_BDY,IERR)
               IF (IERR /= 0) GOTO 500
               print *, minval(covm_bdy),maxval(covm_bdy) 
            ENDIF
         ENDIF  
!----------------------------------------------------------------------------------------!
!                               Begin the iterative loop                                 !
!----------------------------------------------------------------------------------------! 
!
!....... initializaiton
         IERR = 0
         NFUN = 1
         NFUN_SRF = 1
         NFUN_BDY = 1
         LINIT = .TRUE.  !initialization
         CALL MPI_BCAST(LSRCEX ,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(LRECEX ,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(LUPDREC,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
         CALL MPI_BCAST(LUPDSRC,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
         inv%FOBJ_SRF = 0.0
         inv%FOBJ_BDY = 0.0
         inv%FOBJ = 0.0
         F0_BDY = 0.0
         F0_SRF = 0.0
         F0 = 0.0
         DO 2000 K=1,inv%MAXIT
            ICALPHA = 0
            LFUNC = .TRUE.      !need an objective function calculation
            LGRAD = .TRUE.      !need a gradient
            inv%LGNEWT = .TRUE. !inv%LGNEWT = .TRUE. !want Jacobians
            CALL FGH_WINDOW(MASTER,MYID,MYNID,IPGROUP,NPGROUPS, NPPGRP,           &
                            MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM,                  &
                            LUPDREC,LUPDSRC,LSRCEX,LRECEX, LFUNC,LGRAD,           &
                            IBLOCK,K,ICALPHA,                                     &
                            MID,INV,MSH,SRC,RCV,FRQ,M1D,WIN,  IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xjoint: Error calling FGH_WINDOW:',MYID
               GOTO 500
            ENDIF
            IF (MYID == MASTER) THEN
               LSTF = .FALSE.
               IF (LUPDREC) LUPDREC = .FALSE.
               IF (LUPDSRC) THEN
                  LSTF = .TRUE.
                  LUPDSRC = .FALSE. 
               ENDIF
               LSRCEX = .TRUE.
               LRECEX = .TRUE. 
               F0_BDY = inv%FOBJ_BDY
               F0_SRF = inv%FOBJ_SRF
               F0 = inv%FOBJ !F0_BDY + F0_SRF
               if (ljsrch) then !dumps files for jacobian regularzation search program
                  if (inv%lsurf) then
                     call dump_xzinv_dat(projnm,ngnod,msh%nnpg, msh%nelem, inv%na35, &
                                         inv%nvinv,ngnod, .true.,                    &
                                         msh%cdomain,inv%maskg,msh%ieng,             &
                                         msh%xlocs,msh%zlocs, inv%wmask,inv%grad_srf)
                  endif
                  if (inv%lbody) then
                     call dump_xzinv_dat(projnm,ngnod,msh%nnpg, msh%nelem, inv%na35, &
                                         inv%nvinv,ngnod,.false.,                    &
                                         msh%cdomain,inv%maskg,msh%ieng,             &
                                         msh%xlocs,msh%zlocs, inv%wmask,inv%grad_bdy)
                  endif
                  call dump_oe_est(projnm, ndim,frq%nfreq,rcv%nrec, &
                                   ndim,frq%nfreq,rcv%nrec,src%nsrc, &
                                   frq%nfreq_srf_inv,frq%nfreq_bdy_inv, &
                                   src%nsrc_srf,src%nsrc_bdy, inv%irestp, &
                                   msh%azmod, src%pytab, frq%cftype,src%srctyp, &
                                   frq%freq, inv%wghts,inv%obs,inv%est)
               endif
            ENDIF
!
!.......... get gradient onto master 
            IF (MYID == MASTER) THEN
               IF (inv%LSURF) &
               GRAD_SRF8(1:inv%NA35) = DBLE(inv%GRAD_SRF(1:inv%NA35)) !copy gradient
               IF (inv%LBODY) &
               GRAD_BDY8(1:inv%NA35) = DBLE(inv%GRAD_BDY(1:inv%NA35)) !copy gradient
               IF (MYID /= NBHEAD) THEN
                  IDEST = NBHEAD
                  MYTAG = MYID
                  CALL MPI_SEND(LSTF,1,MPI_LOGICAL, IDEST,MYTAG,      &
                                MPI_COMM_WORLD,MPIERR)
                  IF (inv%LSURF) THEN
                     CALL MPI_SEND(GRAD_SRF8,inv%NA35,MPI_DOUBLE_PRECISION, IDEST,MYTAG, &
                                   MPI_COMM_WORLD,MPIERR)
                  ENDIF
                  IF (inv%LBODY) THEN
                     CALL MPI_SEND(GRAD_BDY8,inv%NA35,MPI_DOUBLE_PRECISION, IDEST,MYTAG, &
                                   MPI_COMM_WORLD,MPIERR)
                  ENDIF
                  CALL MPI_SEND(DX,inv%NA35,MPI_DOUBLE_PRECISION, IDEST,MYTAG,  &
                                MPI_COMM_WORLD,MPIERR)
                  IF (LSTF) THEN
                     DO ISRC=1,src%NSRC
                        CALL MPI_SEND(src%SOURCE(:,ISRC),frq%NFREQ,MPI_COMPLEX, &
                                      IDEST,MYTAG, MPI_COMM_WORLD,MPIERR)
                     ENDDO
                  ENDIF
               ENDIF
            ELSE !not master
               IF (MYROW == 0 .AND. MYCOL == 0) THEN
                  IF (.NOT.ALLOCATED(DX)) ALLOCATE(DX(inv%NA35))
                  IF (inv%LSURF) THEN 
                     IF (.NOT.ALLOCATED(GRAD_SRF8)) ALLOCATE(GRAD_SRF8(inv%NA35)) 
                     IF (.NOT.ALLOCATED(SRCH_SRF8)) ALLOCATE(SRCH_SRF8(inv%NA35))
                  ENDIF
                  IF (inv%LBODY) THEN 
                     IF (.NOT.ALLOCATED(GRAD_BDY8)) ALLOCATE(GRAD_BDY8(inv%NA35))
                     IF (.NOT.ALLOCATED(SRCH_BDY8)) ALLOCATE(SRCH_BDY8(inv%NA35)) 
                  ENDIF
                  CALL MPI_RECV(LSTF,1,MPI_LOGICAL, IDEST,MYTAG,      &
                                MPI_COMM_WORLD,STAT,MPIERR)
                  IF (inv%LSURF) THEN
                     CALL MPI_RECV(GRAD_SRF8,inv%NA35,MPI_DOUBLE_PRECISION,   &
                                   MPI_ANY_TAG,MPI_ANY_SOURCE,            &
                                   MPI_COMM_WORLD,STAT,MPIERR)
                  ENDIF
                  IF (inv%LBODY) THEN
                     CALL MPI_RECV(GRAD_BDY8,inv%NA35,MPI_DOUBLE_PRECISION,   &
                                   MPI_ANY_TAG,MPI_ANY_SOURCE,                &
                                   MPI_COMM_WORLD,STAT,MPIERR)
                  ENDIF
                  CALL MPI_RECV(DX,inv%NA35,MPI_DOUBLE_PRECISION,   &
                                MPI_ANY_TAG,MPI_ANY_SOURCE,         &
                                MPI_COMM_WORLD,STAT,MPIERR)
                  IF (LSTF) THEN
                     IF (MYNID /= MASTER) ALLOCATE(src%SOURCE(frq%NFREQ,src%NSRC))
                     DO ISRC=1,src%NSRC
                        CALL MPI_RECV(src%SOURCE(:,ISRC),frq%NFREQ,MPI_COMPLEX,  &
                                      MPI_ANY_TAG,MPI_ANY_SOURCE,                &
                                      MPI_COMM_WORLD,STAT,MPIERR)
                     ENDDO
                  ENDIF !end check on new STF
               ENDIF !end check on process grid locaiton
               IF (MYROW >= 0 .AND. MYCOL >= 0) THEN
                  IF (inv%LSURF) THEN
                     IF (.NOT.ALLOCATED(GRAD_SRF8)) ALLOCATE(GRAD_SRF8(1))
                     IF (.NOT.ALLOCATED(SRCH_SRF8)) ALLOCATE(SRCH_SRF8(1))
                  ENDIF
                  IF (inv%LBODY) THEN
                     IF (.NOT.ALLOCATED(GRAD_BDY8)) ALLOCATE(GRAD_BDY8(1))
                     IF (.NOT.ALLOCATED(SRCH_BDY8)) ALLOCATE(SRCH_BDY8(1))
                  ENDIF 
                  IF (.NOT.ALLOCATED(DX)) ALLOCATE(DX(1))
               ENDIF
            ENDIF !end check on myid
            CALL MPI_BCAST(LSTF,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR) !block
!
!.......... generate Gauss-Newton matrix
            IF (MYROW >= 0 .AND. MYCOL >= 0) THEN
               NPPGRP = NPROCS/NPGROUPS
               IF (inv%LSURF) THEN
                  IF (MYROW == 0 .AND. MYCOL == 0) &
                  WRITE(*,*) 'xjoint: Generating surface wave approximate Hessian...'
                  CALL GENHLOC_SB(ICTXT,NPROW,NPCOL, MBJ,NBJ,MBJT,NBJT, DESCH,         &
                                  LDH,M,N, LSTF,.TRUE., NPPGRP,FRQ,RCV,SRC,INV,        &
                                  HLOC_SRF,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xjoint: Error generating local H matrix!'
                     GOTO 500
                  ENDIF
                  IF (MYROW == 0 .AND. MYCOL == 0) THEN 
                     WRITE(*,*) 'xjoint: Solving (B + inv(C))dp =-g_s + inv(C) (x - x0)'
                     SRCH_SRF8(1:inv%NA35) = 0.D0 
                  ENDIF
                  !CALL DIST_LRANK(LDH, NPROW,NPCOL,ICTXT, .true.,K, 10.d0, DESCH,  &
                  !                Hloc_srf, IERR) 
                  !if (ierr /= 0) then
                  !   write(*,*) 'an error occurred in dist_lrank',myid,ierr
                  !   goto 500
                  !endif
                  !CALL TARANTOLA_394_TR(ICTXT,NPROW,NPCOL, K,LDH,DESCH, &
                  !                      msh%NNPG,inv%NVINV,inv%MASKG, msh%XLOCS,msh%ZLOCS, &
                  !                      kernel,vpvar_srf, sigma_srfx,sigma_srfz, &
                  !                      delta,hloc_srf, grad_srf8, rlam,srch_srf8, ierr) 
                  !call NWALG_442(ICTXT,NPROW,NPCOL, K, LDH,DESCH,DELTA,HLOC_srf, Grad_srf8,  &
                  !               RLAM,Srch_srf8, IERR) 

 
                  !call PICARD_DIST(ICTXT,NPROW,NPCOL, LDH, DESCH,HLOC_srf,Grad_srf8, IERR)
                  CALL TARANTOLA_394(ICTXT, LDH,LDC,NPROW,NPCOL,DESCH,DESCC,GRAD_SRF8,DX,&
                                     HLOC_SRF,COVM_SRF, SRCH_SRF8,IERR)
                  IF (IERR /= 0) THEN 
                     WRITE(*,*) 'xjoint: Error calling tarantola_394 1',MYID
                     GOTO 500
                  ENDIF
                  call est_sigma(ictxt,nprow,npcol, ldh,ldc,.true.,0, desch,descc, & 
                                 hloc_srf,covm_srf, ierr)
               ENDIF !end check on surface wave inversion
               IF (inv%LBODY) THEN
                  IF (MYROW == 0 .AND. MYCOL == 0) &
                  WRITE(*,*) 'xjoint: Generating body wave approximate Hessian...'
                  CALL GENHLOC_SB(ICTXT,NPROW,NPCOL, MBJ,NBJ,MBJT,NBJT, DESCH,         &
                                  LDH,M,N, LSTF,.FALSE., NPPGRP,FRQ,RCV,SRC,INV,       &
                                  HLOC_BDY, IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xjoint: Error generating local H matrix!'
                     GOTO 500
                  ENDIF
                  IF (MYROW == 0 .AND. MYCOL == 0) THEN 
                     WRITE(*,*) 'xjoint: Solving (B + inv(C))dp =-g_b + inv(C) (x - x0)'
                     SRCH_BDY8(1:inv%NA35) = 0.D0
                  ENDIF
                  CALL TARANTOLA_394(ICTXT, LDH,LDC,NPROW,NPCOL,DESCH,DESCC,GRAD_BDY8,DX,&
                                     HLOC_BDY,COVM_BDY, SRCH_BDY8,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xjoint: Error calling tarantola_394 2',MYID
                     GOTO 500
                  ENDIF
                  call est_sigma(ictxt,nprow,npcol, ldh,ldc,.false.,0,desch,descc, &
                                 hloc_bdy,covm_bdy, ierr)
               ENDIF !end check on body wave inversion
               IF (inv%LBODY .AND. inv%LSURF) THEN
                  IF (MYROW == 0 .AND. MYCOL == 0) &
                  WRITE(*,*) 'xjoint: Generating total Hessian...'
                  HLOC_SRF(:,:) = HLOC_SRF(:,:) + inv%SCLBDY*HLOC_BDY(:,:)
                  IF (MYROW == 0 .AND. MYCOL == 0) THEN
                     IF (.NOT.ALLOCATED(GRAD8)) THEN
                        ALLOCATE(GRAD8(inv%NA35))
                        ALLOCATE(SRCH8(inv%NA35))
                     ENDIF
                     GRAD8(1:inv%NA35) = GRAD_SRF8(1:inv%NA35) &
                                       + inv%SCLBDY*GRAD_BDY8(1:inv%NA35)
                     SRCH8(1:inv%NA35) = 0.D0
                  ELSE
                     IF (.NOT.ALLOCATED(GRAD8)) THEN
                        ALLOCATE(GRAD8(1))
                        ALLOCATE(SRCH8(1))
                     ENDIF
                  ENDIF
                  CALL TARANTOLA_394(ICTXT, LDH,LDC,NPROW,NPCOL,DESCH,DESCC,GRAD8,DX, &
                                     HLOC_SRF,COVM_SRF, SRCH8,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xjoint: Error calling tarantola_394 3',MYID
                     GOTO 500
                  ENDIF 
               ENDIF
            ENDIF
!
!.......... get search directions back onto head node and clean space
            IF (MYID == MASTER) THEN !head node recieving
               IF (MYID /= NBHEAD) THEN
                  IF (inv%LSURF) THEN
                     CALL MPI_RECV(SRCH_SRF8,inv%NA35,MPI_DOUBLE_PRECISION,  &
                                   MPI_ANY_TAG,MPI_ANY_SOURCE,               &
                                   MPI_COMM_WORLD,STAT,MPIERR)
                     !CALL WTSRCH8(inv%NA35,inv%WMASK, SRCH_SRF8)
                  ENDIF
                  IF (inv%LBODY) THEN
                     CALL MPI_RECV(SRCH_BDY8,inv%NA35,MPI_DOUBLE_PRECISION,  &
                                   MPI_ANY_TAG,MPI_ANY_SOURCE,               &
                                   MPI_COMM_WORLD,STAT,MPIERR)
                     !CALL WTSRCH8(inv%NA35,inv%WMASK, SRCH_BDY8)
                  ENDIF
                  IF (inv%LSURF .AND. inv%LBODY) THEN
                     CALL MPI_RECV(SRCH8,inv%NA35,MPI_DOUBLE_PRECISION, &
                                   MPI_ANY_TAG,MPI_ANY_SOURCE,          &
                                   MPI_COMM_WORLD,STAT,MPIERR)
                  ENDIF 
               ENDIF
            ELSE 
               IF (MYID == NBHEAD) THEN !send back to master
                  IF (inv%LSURF) THEN
                     CALL MPI_SEND(SRCH_SRF8,inv%NA35,MPI_DOUBLE_PRECISION, IDEST,MYTAG, &
                                   MPI_COMM_WORLD,MPIERR)
                  ENDIF
                  IF (inv%LBODY) THEN
                     CALL MPI_SEND(SRCH_BDY8,inv%NA35,MPI_DOUBLE_PRECISION, IDEST,MYTAG, &
                                   MPI_COMM_WORLD,MPIERR)
                  ENDIF
                  IF (inv%LSURF .AND. inv%LBODY) THEN
                     CALL MPI_SEND(SRCH8,inv%NA35,MPI_DOUBLE_PRECISION, IDEST,MYTAG, &
                                   MPI_COMM_WORLD,MPIERR)
                  ENDIF
               ENDIF
               IF (ALLOCATED(SRCH_SRF8)) DEALLOCATE(SRCH_SRF8)
               IF (ALLOCATED(SRCH_BDY8)) DEALLOCATE(SRCH_BDY8)
               IF (ALLOCATED(GRAD_SRF8)) DEALLOCATE(GRAD_SRF8)
               IF (ALLOCATED(GRAD_BDY8)) DEALLOCATE(GRAD_BDY8)
               IF (ALLOCATED(GRAD8))     DEALLOCATE(GRAD8)
               IF (ALLOCATED(SRCH8))     DEALLOCATE(SRCH8)
               IF (MYID /= MASTER) THEN !master holds onto the apriori model
                  IF (ALLOCATED(DX)) DEALLOCATE(DX)
               ENDIF
               IF (MYNID /= MASTER .AND. MYID == NBHEAD) DEALLOCATE(src%SOURCE) 
            ENDIF !end check on myid
!
!.......... dump gradient
            IF (MYID == MASTER) THEN
               ALLOCATE(SEARCH(inv%NA35)) 
               IF (inv%LSURF) THEN
                  SEARCH(1:inv%NA35) = SNGL(SRCH_SRF8(1:inv%NA35)) 
                  CALL PLOT_SHGRAD_VTK(PROJNM_SRF,NGNOD,msh%NNPG, msh%NELEM, &
                                      inv%CINVTYPE, 1,IBLOCK,K,ICALPHA, &
                                      msh%CDOMAIN,inv%MASKG,msh%IENG,   &
                                      msh%XLOCS,msh%ZLOCS,inv%GRAD_SRF)
                  CALL PLOT_SHGRAD_VTK(PROJNM_SRF,NGNOD,msh%NNPG, msh%NELEM,  &
                                      inv%CINVTYPE, 3,IBLOCK,K,ICALPHA,  &
                                      msh%CDOMAIN,inv%MASKG,msh%IENG,    &
                                      msh%XLOCS,msh%ZLOCS, SEARCH)

               ENDIF
               IF (inv%LBODY) THEN
                  SEARCH(1:inv%NA35) = SNGL(SRCH_BDY8(1:inv%NA35))
                  CALL PLOT_SHGRAD_VTK(PROJNM_BDY,NGNOD,msh%NNPG, msh%NELEM, &
                                      inv%CINVTYPE, 1,IBLOCK,K,ICALPHA, &
                                      msh%CDOMAIN,inv%MASKG,msh%IENG,   &
                                      msh%XLOCS,msh%ZLOCS,inv%GRAD_BDY)
                  CALL PLOT_SHGRAD_VTK(PROJNM_BDY,NGNOD,msh%NNPG, msh%NELEM,  &
                                      inv%CINVTYPE, 3,IBLOCK,K,ICALPHA,  &
                                      msh%CDOMAIN,inv%MASKG,msh%IENG,    &
                                      msh%XLOCS,msh%ZLOCS, SEARCH)

               ENDIF
               IF (inv%LSURF .AND. inv%LBODY) THEN
                  CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM, &
                                       inv%CINVTYPE, 1,IBLOCK,K,ICALPHA, &
                                       msh%CDOMAIN,inv%MASKG,msh%IENG,   &
                                       msh%XLOCS,msh%ZLOCS,inv%GRAD) 
                  SEARCH(1:inv%NA35) = SNGL(SRCH8(1:inv%NA35))
                  CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM, &
                                       inv%CINVTYPE, 3,IBLOCK,K,ICALPHA, &
                                       msh%CDOMAIN,inv%MASKG,msh%IENG,   &
                                       msh%XLOCS,msh%ZLOCS, SEARCH)
               ENDIF
               DEALLOCATE(SEARCH) 
            ENDIF 
            INFO = 1 !assume everything is a-okay with line search
            ICALPHA = 0 
            NFEV = 0
            STPMIN = 0.D0
            ICALPHA_SRF = 0
            ICALPHA_BDY = 0
            NFEV_SRF = 0
            NFEV_BDY = 0
            STPMIN_SRF = 0.D0
            STPMIN_BDY = 0.D0
            ALPHA =-1.D0 !set to fail
            IF (ISTOP_PT == 2) THEN
               IF (MYID == MASTER) WRITE(*,*) 'xjoint: Stopping at stop point 2!'
               GOTO 500
            ENDIF
            IF (inv%LSURF .AND. inv%LBODY) THEN
               IF (MYID == MASTER) THEN
                  IF (MAXVAL(DABS(SRCH8)).EQ.0.D0) THEN
                     WRITE(*,*) 'xjoint: Undefined search direction!'
                     IERR = 1
                     GOTO 500
                  ENDIF
                  STP = GET_STP0(inv%NA35,DP0, XMOD8,SRCH8)
                  STPMIN = 5.D0/MAXVAL(DABS(SRCH8),MASK=DABS(SRCH8) > 0.D0)
                  IF (STPMIN > STP) THEN
                     WRITE(*,*) 'xjoint: Youve made a mistake in stpmin'
                     WRITE(*,*) 'xjoint: Defaulting stmin to 1.E-20'
                     STPMIN = 1.D-20
                  ENDIF
                  STP0 = STP
               ENDIF
               FMIN = DBLE(inv%FOBJ)
               GOTO 3020
            ENDIF
            IF (MYID == MASTER) THEN 
               LINV_BDY = inv%LBODY !save
               LINV_SRF = inv%LSURF !save
               IF (inv%LSURF .AND. inv%LBODY) inv%LBODY = .FALSE.  !turn off to start
               STP_SRF = 0.D0
               STP_BDY = 0.D0
               IF (inv%LSURF) THEN
                  !XAVG = SUM(DABS(SRCH_SRF8))/DFLOAT(inv%NA35)
                  !IF (XAVG == 0.D0) THEN
                  !   WRITE(*,*) 'xjoint: Undefined surface wave search direction!'
                  !   IERR = 1
                  !   GOTO 500
                  !ENDIF
                  !STP_SRF = 1.D0/XAVG*DS0 
                  !IF (MAXVAL(DABS(SRCH_SRF8)) == 0.D0) THEN
                  !   WRITE(*,*) 'xjoint: Division by zero in srch_srf8!'
                  !   IERR = 1
                  !   GOTO 500
                  !ENDIF
                  STP_SRF = GET_STP0(inv%NA35,DP0, XMOD_SRF8,SRCH_SRF8)
                  STPMIN_SRF = 5.D0/MAXVAL(DABS(SRCH_SRF8),MASK=DABS(SRCH_SRF8) > 0.D0) 
                  IF (STPMIN_SRF > STP_SRF) THEN
                     WRITE(*,*) 'xjoint: Youve made a mistake in stpmin'
                     WRITE(*,*) 'xjoint: Defaulting stpmin to 1.E-20'
                     STPMIN_SRF = 1.D-20
                  ENDIF
                  IF (MAXVAL(DABS(SRCH_SRF8)).EQ.0.D0) THEN
                     WRITE(*,*) 'xjoint: Undefined surface wave search direction!'
                     IERR = 1
                     GOTO 500
                  ENDIF
                  !STP_SRF = 1.D0/MAXVAL(DABS(SRCH_SRF8))*DS0
               ENDIF
               IF (inv%LBODY) THEN
                  XAVG = SUM(DABS(SRCH_BDY8))/DFLOAT(inv%NA35)
                  IF (XAVG == 0.D0) THEN
                     WRITE(*,*) 'xjoint: Undefined body wave search direction!'
                     IERR = 1
                     GOTO 500
                  ENDIF
                  STP_BDY = 1.D0/XAVG*DS0
                  STPMIN_BDY = 5.D0/MAXVAL(DABS(SRCH_BDY8),MASK=DABS(SRCH_BDY8) > 0.D0)
                  IF (STPMIN_BDY > STP_BDY) THEN
                     WRITE(*,*) 'xjoint: Youve made a mistake in stpmin'
                     WRITE(*,*) 'xjoint: Defaulting stpmin_bdy to 1.E-20'
                     STPMIN_BDY = 1.D-20
                  ENDIF
                  IF (MAXVAL(DABS(SRCH_BDY8)).EQ.0.D0) THEN
                     WRITE(*,*) 'xjoint: Undefined body wave search direction!'
                     IERR = 1
                     GOTO 500
                  ENDIF
                  !STP_BDY = 1.D0/MAXVAL(DABS(SRCH_BDY8))*DS0 
               ENDIF
               CALL MASK_LIFREQ(frq%NFREQ,.TRUE.,1, frq%CFTYPE, frq%LINVF,LFSAVE) 
               print *, 'stp0',stp_srf,stp_bdy
               FSRF8 = 0.D0
               FBDY8 = 0.D0
               STPMIN = 1.D-20
               IF (inv%LSURF) THEN
                  FSRF8 = DBLE(inv%FOBJ_SRF)
                  STPMIN = STPMIN_SRF
                  FMIN = FSRF8
               ELSE
                  FBDY8 = DBLE(inv%FOBJ_BDY)
                  STPMIN = STPMIN_BDY
                  FMIN = FBDY8
               ENDIF
               IF (inv%LBODY) FBDY8 = DBLE(inv%FOBJ_BDY)
            ENDIF
            CALL MPI_BCAST(inv%LBODY,        1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
            CALL MPI_BCAST(inv%LSURF,        1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
            CALL MPI_BCAST(frq%LINVF,frq%NFREQ,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
!----------------------------------------------------------------------------------------!
!                         This is the surface wave line search loop                      !
!----------------------------------------------------------------------------------------!
!
!............. this is the surface wave section
  3000         CONTINUE !loop for another surface wave iteration
               IF (inv%LSURF) THEN
                  CALL MPI_BCAST(inv%LFUNC,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
                  CALL MPI_BCAST(inv%LGRAD,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
                  IF (inv%LFUNC .OR. inv%LGRAD) THEN
                     LUPDREC = .FALSE.   !no receiver update 
                     LUPDSRC = .FALSE.   !no source update
                     LFUNC = .TRUE.      !need an objective function calculation
                     LGRAD = .TRUE.      !need a gradient
                     inv%LGNEWT = .FALSE.!dont want Jacobians; unless windowing
                     CALL FGH_WINDOW(MASTER,MYID,MYNID,IPGROUP,NPGROUPS, NPPGRP,         &
                                     MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM,                &
                                     LUPDREC,LUPDSRC,LSRCEX,LRECEX, LFUNC,LGRAD,         &
                                     IBLOCK,K,ICALPHA_SRF,                               &
                                     MID,INV,MSH,SRC,RCV,FRQ,M1D,WIN,  IERR)
                     IF (IERR /= 0) THEN
                        WRITE(*,*) 'xjoint: An error occurred in FGH_WINDOW srf2'
                        GOTO 500
                     ENDIF
                     IF (MYID == MASTER) THEN
                        FSRF8 = DBLE(inv%FOBJ_SRF)
                        GRAD_SRF8(1:inv%NA35) = DBLE(inv%GRAD_SRF(1:inv%NA35))
                        NFUN_SRF = NFUN_SRF + 1
                        NFUN = NFUN + 1
                        CALL PLOT_SHGRAD_VTK(PROJNM_SRF,NGNOD,msh%NNPG, msh%NELEM, &
                                             inv%CINVTYPE,1,IBLOCK,K,ICALPHA_SRF,  &
                                             msh%CDOMAIN,inv%MASKG,msh%IENG,       &
                                             msh%XLOCS,msh%ZLOCS,inv%GRAD_SRF)
                        IF (FSRF8 < FMIN) THEN
                           ALPHA = STP_SRF 
                           FMIN = FSRF8 
                           CALL DCOPY(inv%NA35,XMOD_SRF8,1,XSAVE,1) 
                        ENDIF
                     ENDIF
                     inv%LFUNC = .FALSE. 
                     inv%LGRAD = .FALSE.
                  ENDIF

                  IF (MYID == MASTER) THEN
                     print *, minval(xmod_srf8), maxval(xmod_srf8)
                     print *, fsrf8, minval(grad_srf8), maxval(grad_srf8)
                     CALL MCSRCH(inv%NA35,XMOD_SRF8,FSRF8,GRAD_SRF8,SRCH_SRF8, STP_SRF,  &
                                 FTOL8,GTOL8,XTOL8,STPMIN,STPMAX,                        &
                                 MAXFEV,INFO,NFEV_SRF,WA,LP)
                     print *, minval(xmod_srf8), maxval(xmod_srf8), info
                     MYDEST =-1
                     CALL MCSRCH_ERROR(INFO,MAXFEV,MYDEST,IERR)
                     IF (INFO == 3) THEN !max function evaluations
                        PRINT *, ALPHA
                        IF (ALPHA /=-1.D0) THEN
                           WRITE(*,*) 'xjoint: Overriding MCSRCH and taking min obj fn'
                           STP_SRF = ALPHA
                           inv%FOBJ_SRF = SNGL(FMIN)
                           CALL DCOPY(inv%NA35,XSAVE,1,XMOD_SRF8,1)
                           INFO = 1
                           IERR = 0 
                        ENDIF
                     ENDIF 
                     print *, 'mydest =',mydest
                     IF (IERR /= 0) GOTO 500 
!
!................... convergence test
                     IF (INFO == 1) THEN
                        ICALPHA_SRF = 0
                        MYDEST = 2000 !default, continue with iterative loop 
                        WRITE(*,8151) K,NFUN,NFUN_SRF,STP_SRF,inv%FOBJ_SRF, &
                                      100. - inv%FOBJ_SRF/F0_SRF*100.0
                        IF (inv%FOBJ_SRF/F0_SRF*100. < 20.0) THEN
                           WRITE(*,*) 'xjoint: Convergence achieved'
                           MYDEST = 1005
                        ENDIF
                     ENDIF 
!
!................... this is the line search
                     IF (INFO ==-1) THEN
                        MYDEST = 3000
                        inv%LFUNC = .TRUE.
                        inv%LGRAD = .TRUE.
                        ICALPHA_SRF = ICALPHA_SRF + 1 !this is another step calculation
                        WRITE(*,8160) ICALPHA_SRF,NFUN_SRF,STP_SRF,inv%FOBJ_SRF 
                     ENDIF 
                     WRITE(*,*) 'xjoint: Updating surface wave model...'
                     CALL UPDMOD8(PROJNM_SRF,LFILES,inv%NA35, IBLOCK,K,ICALPHA_SRF,   &
                                  VPMIN_INV,VPMAX_INV, XMOD_SRF8, INV,MSH,IERR)
                     IF (IERR /= 0) THEN
                        WRITE(*,*) 'xjoint: Error calling updmod!'
                        IERR = 1 
                        GOTO 500
                     ENDIF
                  ENDIF !end check on surface wave inversion
                  CALL MPI_BCAST(MYDEST,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
                  CALL BCAST_MOD_UPD(MPI_COMM_WORLD,MASTER, MSH)
                  IF (MYDEST == 3000) GOTO 3000
!
!................ converged write residuals
                  IF (MYID == MASTER) THEN
                     CALL DUMP_RESID(PROJNM,NDIM,frq%NFREQ,rcv%NREC,                  &
                               NDIM,frq%NFREQ,rcv%NREC,src%NSRC,                      &
                               .TRUE.,.TRUE.,inv%IRESTP,K,                            & 
                               msh%AZMOD,frq%LINVF,frq%CFTYPE,src%SRCTYP,frq%FREQ,    &
                               inv%WGHTS,inv%OBS,inv%EST) 
                  ENDIF
               ENDIF !end check on surface wave inversion
!----------------------------------------------------------------------------------------!
!                                 Surface wave line search done                          !
!----------------------------------------------------------------------------------------!
!
!.......... continue forward onto body waves
            IF (MYID == MASTER) THEN
               inv%LBODY = LINV_BDY
               inv%LSURF = LINV_SRF
               IF (inv%LSURF.AND.inv%LBODY) inv%LSURF = .FALSE. !turn off surface waves
               IF (inv%LBODY) THEN
                  STPMIN = STPMIN_BDY
                  FMIN = FBDY8 
                  IF (inv%LSURF) THEN
                     WRITE(*,*) 'xjoint: Restoring original model...'
                     CALL UPDMOD8(PROJNM_BDY,LFILES,inv%NA35, IBLOCK,K,ICALPHA_BDY, &   
                                  VPMIN_INV,VPMAX_INV, XMOD_BDY8, INV,MSH,IERR)
                  ENDIF  
               ENDIF 
               CALL MASK_LIFREQ(frq%NFREQ,.TRUE. ,2, frq%CFTYPE, frq%LINVF,LFSAVE)
               CALL MASK_LIFREQ(frq%NFREQ,.FALSE.,1, frq%CFTYPE, frq%LINVF,LFSAVE)
            ENDIF 
            CALL MPI_BCAST(inv%LBODY,        1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
            CALL MPI_BCAST(inv%LSURF,        1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
            CALL MPI_BCAST(frq%LINVF,frq%NFREQ,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
            IF (inv%LSURF .AND. inv%LBODY) &
            CALL BCAST_MOD_UPD(MPI_COMM_WORLD,MASTER, MSH)
!----------------------------------------------------------------------------------------!
!                          This is the body wave line search loop                        !
!----------------------------------------------------------------------------------------!
!
!............. this is the body wave section
 3010          CONTINUE  
               IF (inv%LBODY) THEN
                  CALL MPI_BCAST(inv%LFUNC,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
                  CALL MPI_BCAST(inv%LGRAD,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
                  IF (inv%LFUNC .OR. inv%LGRAD) THEN 
                     LUPDREC = .FALSE.   !no receiver update 
                     LUPDSRC = .FALSE.   !no source update
                     LFUNC = .TRUE.      !need an objective function calculation
                     LGRAD = .TRUE.      !need a gradient
                     inv%LGNEWT = .FALSE.!dont want Jacobians; unless windowing
                     CALL FGH_WINDOW(MASTER,MYID,MYNID,IPGROUP,NPGROUPS, NPPGRP,         &
                                     MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM,                &
                                     LUPDREC,LUPDSRC,LSRCEX,LRECEX, LFUNC,LGRAD,         &
                                     IBLOCK,K,ICALPHA_BDY,                               &
                                     MID,INV,MSH,SRC,RCV,FRQ,M1D,WIN,  IERR)
                     IF (IERR /= 0) THEN 
                        WRITE(*,*) 'xjoint: An error occurred in FGH_WINDOW srf2'
                        GOTO 500
                     ENDIF
                     IF (MYID == MASTER) THEN 
                        FBDY8 = DBLE(inv%FOBJ_BDY)
                        GRAD_BDY8(1:inv%NA35) = DBLE(inv%GRAD_BDY(1:inv%NA35))
                        NFUN_BDY = NFUN_BDY + 1
                        NFUN = NFUN + 1
                        CALL PLOT_SHGRAD_VTK(PROJNM_BDY,NGNOD,msh%NNPG, msh%NELEM, &
                                             inv%CINVTYPE,1,IBLOCK,K,ICALPHA_BDY,  &
                                             msh%CDOMAIN,inv%MASKG,msh%IENG,       &
                                             msh%XLOCS,msh%ZLOCS,inv%GRAD_BDY)
                        IF (FBDY8 < FMIN) THEN
                           ALPHA = STP_BDY
                           FMIN = FBDY8
                           CALL DCOPY(inv%NA35,XMOD_BDY8,1,XSAVE,1)
                        ENDIF
                     ENDIF
                     inv%LFUNC = .FALSE. 
                     inv%LGRAD = .FALSE.
                  ENDIF

                  IF (MYID == MASTER) THEN
                     print *, minval(xmod_bdy8), maxval(xmod_bdy8)
                     CALL MCSRCH(inv%NA35,XMOD_BDY8,FBDY8,GRAD_BDY8,SRCH_BDY8, STP_BDY, &
                                 FTOL8,GTOL8,XTOL8,STPMIN,STPMAX,                        &
                                 MAXFEV,INFO,NFEV_BDY,WA,LP)
                     print *, minval(xmod_bdy8), maxval(xmod_bdy8), info
                     MYDEST =-1
                     CALL MCSRCH_ERROR(INFO,MAXFEV,MYDEST,IERR)
                     print *, 'mydest =',mydest
                     IF (INFO == 3) THEN !max function evaluations
                        PRINT *, ALPHA
                        IF (ALPHA /=-1.D0) THEN
                           WRITE(*,*) 'xjoint: Overriding MCSRCH and taking min obj fn'
                           STP_BDY = ALPHA
                           inv%FOBJ_BDY = SNGL(FMIN)
                           CALL DCOPY(inv%NA35,XSAVE,1,XMOD_BDY8,1)
                           INFO = 1
                           IERR = 0
                        ENDIF
                     ENDIF
                     IF (IERR /= 0) GOTO 500
!
!................... convergence test
                     IF (INFO == 1) THEN
                        ICALPHA_BDY = 0
                        MYDEST = 2000 !default, continue with iterative loop 
                        WRITE(*,8152) K,NFUN,NFUN_BDY,STP_BDY,inv%FOBJ_BDY, &
                                      100. - inv%FOBJ_BDY/F0_BDY*100.0
                        IF (inv%FOBJ_BDY/F0_BDY*100. < 20.0) THEN
                           WRITE(*,*) 'xjoint: Convergence achieved'
                           MYDEST = 1005
                        ENDIF
                     ENDIF
!
!................... this is the line search
                     IF (INFO ==-1) THEN
                        MYDEST = 3010
                        inv%LFUNC = .TRUE.
                        inv%LGRAD = .TRUE.
                        ICALPHA_BDY = ICALPHA_BDY + 1 !this is another step calculation
                        WRITE(*,8161) ICALPHA_BDY,NFUN_BDY,STP_BDY,inv%FOBJ_BDY
                     ENDIF
                     WRITE(*,*) 'xjoint: Updating body wave model...'
                     CALL UPDMOD8(PROJNM_BDY,LFILES,inv%NA35, IBLOCK,K,ICALPHA_BDY,   &
                                  VPMIN_INV,VPMAX_INV, XMOD_BDY8, INV,MSH,IERR)
                     IF (IERR /= 0) THEN
                        WRITE(*,*) 'xjoint: Error calling updmod!'
                        IERR = 1
                        GOTO 500
                     ENDIF
                  ENDIF !end check on surface wave inversion
                  CALL MPI_BCAST(MYDEST,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
                  CALL BCAST_MOD_UPD(MPI_COMM_WORLD,MASTER, MSH)
                  IF (MYDEST == 3010) GOTO 3010
!
!................ converged write residuals
                  IF (MYID == MASTER) THEN 
                     CALL DUMP_RESID(PROJNM,NDIM,frq%NFREQ,rcv%NREC,                &
                               NDIM,frq%NFREQ,rcv%NREC,src%NSRC,                    &
                               .FALSE.,.FALSE.,inv%IRESTP,K,                        &
                               msh%AZMOD,frq%LINVF,frq%CFTYPE,src%SRCTYP,frq%FREQ,  &
                               inv%WGHTS,inv%OBS,inv%EST) 
                  ENDIF
               ENDIF !end check on body wave inversion
!----------------------------------------------------------------------------------------!
!                                    Body wave line search                               !
!----------------------------------------------------------------------------------------!
            IF (MYID == MASTER) THEN
               inv%LBODY = LINV_BDY
               inv%LSURF = LINV_SRF
               CALL MASK_LIFREQ(frq%NFREQ,.FALSE.,2, frq%CFTYPE, frq%LINVF,LFSAVE)
            ENDIF
            CALL MPI_BCAST(inv%LBODY,        1,MPI_LOGICAL,MASTER,MPI_COMM_WORLD,MPIERR)
            CALL MPI_BCAST(inv%LSURF,        1,MPI_LOGICAL,MASTER,MPI_COMM_WORLD,MPIERR)
            CALL MPI_BCAST(frq%LINVF,frq%NFREQ,MPI_LOGICAL,MASTER,MPI_COMM_WORLD,MPIERR)
!----------------------------------------------------------------------------------------!
!                                    This is the joint loop                              !
!----------------------------------------------------------------------------------------!
!
!............. this is the joint inversion section
 3020          CONTINUE  
               IF (inv%LSURF .AND. inv%LBODY) THEN 
                  CALL MPI_BCAST(inv%LFUNC,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
                  CALL MPI_BCAST(inv%LGRAD,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
                  IF (inv%LFUNC .OR. inv%LGRAD) THEN 
                     LUPDREC = .FALSE.   !no receiver update 
                     LUPDSRC = .FALSE.   !no source update
                     LFUNC = .TRUE.      !need an objective function calculation
                     LGRAD = .TRUE.      !need a gradient
                     inv%LGNEWT = .FALSE.!dont want Jacobians; unless windowing
                     CALL FGH_WINDOW(MASTER,MYID,MYNID,IPGROUP,NPGROUPS, NPPGRP,         &    
                                     MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM,                &    
                                     LUPDREC,LUPDSRC,LSRCEX,LRECEX, LFUNC,LGRAD,         &    
                                     IBLOCK,K,ICALPHA,                                   &
                                     MID,INV,MSH,SRC,RCV,FRQ,M1D,WIN,  IERR)
                     IF (IERR /= 0) THEN 
                        WRITE(*,*) 'xjoint: An error occurred in FGH_WINDOW joint'
                        GOTO 500
                     ENDIF
                     IF (MYID == MASTER) THEN
                        FOBJ8 = DBLE(inv%FOBJ) 
                        print *, 'fobj8:',fobj8
                        GRAD8(1:inv%NA35) = DBLE(inv%GRAD(1:inv%NA35)) 
                        NFUN = NFUN + 1
                        CALL PLOT_SHGRAD_VTK(PROJNM,NGNOD,msh%NNPG, msh%NELEM,     &
                                             inv%CINVTYPE,1,IBLOCK,K,ICALPHA,      &
                                             msh%CDOMAIN,inv%MASKG,msh%IENG,       &
                                             msh%XLOCS,msh%ZLOCS,inv%GRAD)
                        IF (FOBJ8 < FMIN) THEN
                          print *, 'fobj8',fobj8
                           ALPHA = STP
                           FMIN = FOBJ8
                           CALL DCOPY(inv%NA35,XMOD8,1,XSAVE,1)
                        ENDIF
                     ENDIF
                     inv%LFUNC = .FALSE.
                     inv%LGRAD = .FALSE.
                  ENDIF
                  IF (MYID == MASTER) THEN
                     print *, minval(xmod8), maxval(xmod8)
                     CALL MCSRCH(inv%NA35,XMOD8,FOBJ8,GRAD8,SRCH8, STP,  &
                                 FTOL8,GTOL8,XTOL8,STPMIN,STPMAX,       &
                                 MAXFEV,INFO,NFEV,WA,LP)
                     !CALL INTSTP(.TRUE.,inv%NA35, FOBJ8,FTOL8,GTOL8, STP0, &
                     !            XMOD8,GRAD8,SRCH8,  STP,NFEV, INFO)
                     print *, minval(xmod8), maxval(xmod8), info
                     MYDEST =-1
                     CALL MCSRCH_ERROR(INFO,MAXFEV,MYDEST,IERR)
                     print *, 'mydest =',mydest,info
                     IF (INFO == 3) THEN !max function evaluations
                        PRINT *, ALPHA
                        IF (ALPHA /=-1.D0) THEN
                           WRITE(*,*) 'xjoint: Overriding MCSRCH and taking min obj fn'
                           STP = ALPHA
                           inv%FOBJ = SNGL(FMIN)
                           CALL DCOPY(inv%NA35,XSAVE,1,XMOD8,1)
                           INFO = 1
                           IERR = 0
                        ENDIF
                     ENDIF
                     IF (IERR /= 0) GOTO 500
!
!................... convergence test
                     IF (INFO == 1) THEN
                        ICALPHA = 0
                        MYDEST = 2000 !default, continue with iterative loop 
                        WRITE(*,8153) K,NFUN,NFUN,STP,inv%FOBJ, &
                                      100. - inv%FOBJ/F0*100.0
                        IF (inv%FOBJ/F0*100. < 20.0) THEN
                           WRITE(*,*) 'xjoint: Convergence achieved'
                           MYDEST = 1005
                        ENDIF
                     ENDIF
!
!................... this is the line search
                     IF (INFO ==-1) THEN
                        MYDEST = 3020
                        inv%LFUNC = .TRUE.
                        inv%LGRAD = .TRUE.
                        ICALPHA = ICALPHA + 1 !this is another step calculation
                        WRITE(*,8162) ICALPHA,NFUN,STP,inv%FOBJ,inv%FOBJ_SRF,inv%FOBJ_BDY
                     ENDIF
                     print *, 'mydest =',mydest
                     WRITE(*,*) 'xjoint: Updating joint model...'
                     CALL UPDMOD8(PROJNM,LFILES,inv%NA35, IBLOCK,K,ICALPHA,   &
                                  VPMIN_INV,VPMAX_INV, XMOD8, INV,MSH,IERR)
                     IF (IERR /= 0) THEN
                        WRITE(*,*) 'xjoint: Error calling updmod!'
                        IERR = 1
                        GOTO 500
                     ENDIF
                  ENDIF !end check on surface wave inversion
                  CALL MPI_BCAST(MYDEST,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR)
                  CALL BCAST_MOD_UPD(MPI_COMM_WORLD,MASTER, MSH)
                  IF (MYDEST == 3020) GOTO 3020
               ENDIF !end check on surface + body wave inversion
!----------------------------------------------------------------------------------------!
!                      This is the end of the joint inversion line search                !
!----------------------------------------------------------------------------------------!
!
!.......... flash initial objective function
            IF (MYID == MASTER) THEN
               IF (LINIT) THEN 
                  WRITE(*,8145) IBLOCK, inv%FOBJ_SRF, inv%FOBJ_BDY, inv%FOBJ 
                  LINIT = .FALSE.
 8145             FORMAT(' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
                         ' + xjoint: Beginning iterative inversion                 +',/, &
                         ' +         Inversion block:',I6,'                        +',/, &
                         ' +         Surface wave objective function:',E12.4,'  +',/, &
                         ' +         Body wave objective function:',E12.4,'     +',/,&
                         ' +         Objective function:',E12.4,'               +',/, & 
                         ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/)
               ENDIF
            ENDIF
!----------------------------------------------------------------------------------------!
!                   May want to approximate higher order scattering                      !
!----------------------------------------------------------------------------------------!
            IF (inv%LAHESS_HOT) THEN
               IF (MYID == MASTER) THEN
                  xmod_bdy8(:) = xmod_bdy8(:) - stp_bdy*srch_bdy8(:) !restore model
                  ALLOCATE(P1(inv%NA35))
                  P1(1:inv%NA35) = 0.0
                  IF (inv%LSURF) THEN
                     ALLOCATE(PSSRF(inv%NA35))
                     P1(1:inv%NA35) = SNGL(STP_SRF*SRCH_SRF8(1:inv%NA35))
                  ELSE
                     ALLOCATE(PSSRF(1))
                  ENDIF
                  IF (inv%LBODY) THEN
                     ALLOCATE(PSBDY(inv%NA35))
                     P1(1:inv%NA35) = P1(1:inv%NA35) &
                                    + SNGL(STP_BDY*SRCH_BDY8(1:inv%NA35))
                  ELSE
                     ALLOCATE(PSBDY(1)) 
                  ENDIF
                  print *, minval(p1), maxval(p1)
               ELSE
                  ALLOCATE(P1(1))
                  ALLOCATE(PSSRF(1))
                  ALLOCATE(PSBDY(1))
                  P1(1) = 0.0
                  PSSRF(1) = 0.0
                  PSBDY(1) = 0.0
               ENDIF
               if (myid == master) then
                  call updmod8('restore',LFILES,inv%NA35, iblock,K,0,        &
                               vpmin_inv,vpmax_inv, xmod_bdy8, inv,msh,ierr)
               endif
               call bcast_mod_upd(mpi_comm_world,master, msh)

               CALL QPMFREE(MASTER,MYID,MYNID,IPGROUP,NPGROUPS,NPPGRP,   &
                            MPI_COMM_WORLD,MYSLV_COMM,MYHD_COMM, PROJNM, &
                            P1, WIN,INV,MSH,SRC,RCV,FRQ,M1D, MID,        &
                            PSSRF,PSBDY, IERR)
               deallocate(p1)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xjoint: Error calling qpmfree!'
                  GOTO 500
               ENDIF
               if (myid == master) then
                  idest = nbhead
                  mytag = myid
                  if (myrow /= 0 .or. mycol /= 0) then 
                     if (inv%LSURF) then
                        if (.not.allocated(grad_srf8)) allocate(grad_srf8(inv%na35))
                        grad_srf8(1:inv%NA35) = dble(pssrf(1:inv%NA35))
                        call mpi_send(grad_srf8,inv%NA35,mpi_double_precision, idest,mytag, &
                                      mpi_comm_world,mpierr)
                     endif
                     if (inv%LBODY) then
                        if (.not.allocated(grad_bdy8)) allocate(grad_bdy8(inv%NA35))
                        grad_bdy8(1:inv%NA35) = dble(psbdy(1:inv%NA35))
                        call mpi_send(grad_bdy8,inv%NA35,mpi_double_precision, idest,mytag, &
                                      mpi_comm_world,mpierr)
                     endif
                  else
                     if (inv%LSURF) then
                        if (.not.allocated(grad_srf8)) allocate(grad_srf8(inv%na35))
                        grad_srf8(1:inv%NA35) = dble(pssrf(1:inv%NA35))
                     endif
                     if (inv%LBODY) then
                        if (.not.allocated(grad_bdy8)) allocate(grad_bdy8(inv%na35))
                        grad_bdy8(1:inv%NA35) = dble(psbdy(1:inv%NA35))
                     endif
                  endif 
               else 
                  if (myid == nbhead) then 
                     if (inv%LSURF) then
                        if (.not.allocated(grad_srf8)) allocate(grad_srf8(inv%NA35))
                        call mpi_recv(grad_srf8,inv%NA35,mpi_double_precision,   &
                                      mpi_any_tag,mpi_any_source,                &
                                      mpi_comm_world,stat,mpierr)
                     endif
                     if (inv%LBODY) THEN
                        if (.not.allocated(grad_bdy8)) allocate(grad_bdy8(inv%NA35))
                        call mpi_recv(grad_bdy8,inv%NA35,mpi_double_precision,   &
                                      mpi_any_tag,mpi_any_source,                &
                                      mpi_comm_world,stat,mpierr)
                     endif
                  endif
               endif
!
!............. solve p2 =-inv(B)Qp
               if (myrow >= 0 .and. mycol >= 0) then 
                  if (myrow == 0 .and. mycol == 0) then
                     if (inv%lsurf) then
                        if (.not.allocated(srch_srf8)) allocate(srch_srf8(inv%NA35))
                     endif
                     if (inv%lbody) then 
                        if (.not.allocated(srch_bdy8)) allocate(srch_bdy8(inv%NA35))
                     endif
                     if (.not.allocated(dx)) allocate(dx(inv%NA35))
                     dx(1:inv%NA35) = 0.D0
                  else
                     if (inv%lsurf) then
                        if (.not.allocated(srch_srf8)) allocate(srch_srf8(1))
                        if (.not.allocated(grad_srf8)) allocate(grad_srf8(1))
                     endif
                     if (inv%lbody) then
                        if (.not.allocated(srch_bdy8)) allocate(srch_bdy8(1))
                        if (.not.allocated(grad_bdy8)) allocate(grad_bdy8(1))
                     endif
                     if (.not.allocated(dx)) allocate(dx(1))
                  endif
                  if (inv%LSURF) then
                     call tarantola_394(ictxt,ldh,ldc,nprow,npcol,desch,descc,grad_srf8,dx,&
                                        hloc_srf,covm_srf, srch_srf8,ierr)
                     if (ierr /= 0) then
                        write(*,*) 'xjoint: ERror calling tarantola_394!'
                        goto 500
                     endif
                  endif
                  if (inv%LBODY) then
                     call tarantola_394(ictxt,ldh,ldc,nprow,npcol,desch,descc,grad_bdy8,dx,&
                                        hloc_bdy,covm_bdy, srch_bdy8,ierr)
                     if (ierr /= 0) then
                        write(*,*) 'xjoint: ERror calling tarantola_394'
                        goto 500
                     endif
                  endif
               endif
               call mpi_barrier(mpi_comm_world,mpierr)
!
!............. get the updates onto head node 
               if (myid == master) then
                  if (myid /= nbhead) then
                     if (inv%LSURF) then
                        call mpi_recv(srch_srf8,inv%NA35,mpi_double_precision,  &
                                      mpi_any_tag,mpi_any_source,               &
                                      mpi_comm_world,stat,mpierr)
                     endif
                     if (inv%LBODY) then
                        call mpi_recv(srch_bdy8,inv%NA35,mpi_double_precision,  &
                                      mpi_any_tag,mpi_any_source,               &
                                      mpi_comm_world,stat,mpierr)
                     endif
                  endif
               else
                  if (myid == nbhead) then !send back to master
                     if (inv%LSURF) then
                        call mpi_send(srch_srf8,inv%NA35,mpi_double_precision, idest,mytag, &
                                      mpi_comm_world,mpierr)
                     endif
                     if (inv%LBODY) then
                        call mpi_send(srch_bdy8,inv%NA35,mpi_double_precision, idest,mytag, &
                                      mpi_comm_world,mpierr)
                     endif
                  endif
               endif
!
!............. copy onto search directions and plot 
               IF (MYID == MASTER) THEN
                  IF (inv%LSURF) THEN
                     PSSRF(1:inv%NA35) = SNGL(SRCH_SRF8(1:inv%NA35))
                     CALL PLOT_SHGRAD_VTK(PROJNM_SRF,NGNOD,msh%NNPG, msh%NELEM,      &
                                          inv%CINVTYPE, 5,IBLOCK,K,0,                &
                                          msh%CDOMAIN,inv%MASKG,msh%IENG,            &    
                                          msh%XLOCS,msh%ZLOCS, PSSRF)
                  ENDIF
                  IF (inv%LBODY) THEN
                     PSBDY(1:inv%NA35) = sngl(SRCH_BDY8(1:inv%NA35))
                     CALL PLOT_SHGRAD_VTK(PROJNM_BDY,NGNOD,msh%NNPG, msh%NELEM,      &
                                          inv%CINVTYPE, 5,IBLOCK,K,0,                &
                                          msh%CDOMAIN,inv%MASKG,msh%IENG,            &
                                          msh%XLOCS,msh%ZLOCS, PSBDY)
                  ENDIF
               ENDIF
               IF (ALLOCATED(SRCH_SRF8)) DEALLOCATE(SRCH_SRF8)
               IF (ALLOCATED(SRCH_BDY8)) DEALLOCATE(SRCH_BDY8)
               IF (ALLOCATED(GRAD_SRF8)) DEALLOCATE(GRAD_SRF8)
               IF (ALLOCATED(GRAD_BDY8)) DEALLOCATE(GRAD_BDY8)
               IF (ALLOCATED(DX))        DEALLOCATE(DX)
               IF (ALLOCATED(PSSRF)) DEALLOCATE(PSSRF)
               IF (ALLOCATED(PSBDY)) DEALLOCATE(PSBDY)
            ENDIF

 2000    CONTINUE !end iterative loop
!----------------------------------------------------------------------------------------!
!                          This concludes the iterative loop                             !
!----------------------------------------------------------------------------------------!
!
!....... free space for next block 
         IF (MYID == MASTER) THEN
            IF (ALLOCATED(LFSAVE)) DEALLOCATE(LFSAVE)
         ENDIF
         IF (ASSOCIATED(inv%WGHTS))  DEALLOCATE(inv%WGHTS)
         IF (ASSOCIATED(inv%OBS))    DEALLOCATE(inv%OBS)
         IF (MYNID == MASTER) THEN
            IF (ASSOCIATED(inv%EST))    DEALLOCATE(inv%EST)
            IF (ASSOCIATED(src%SOURCE)) DEALLOCATE(src%SOURCE)
         ENDIF
         IF (ASSOCIATED(frq%FREQ))   DEALLOCATE(frq%FREQ)
         IF (ASSOCIATED(frq%CFTYPE)) DEALLOCATE(frq%CFTYPE)
 1000 CONTINUE !loop on frequency blocks
!----------------------------------------------------------------------------------------!
!                   This concludes the outer block on frequency groups                   !
!----------------------------------------------------------------------------------------!
!
!.... break ahead for errors
  500 CONTINUE 
      IF (IERR /= 0) THEN
         WRITE(*,*) 'xjoint: There was an error on process:',MYID
         CALL MPI_ABORT(MPI_COMM_WORLD,30,MPIERR)
      ENDIF 
      IF (MYID == MASTER) WRITE(*,*) 'xjoint: Freeing space...'
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
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
      IF (ASSOCIATED(src%SRCTYP))   DEALLOCATE(src%SRCTYP)

      IF (ALLOCATED(GRAD_SRF8))     DEALLOCATE(GRAD_SRF8)
      IF (ALLOCATED(GRAD_BDY8))     DEALLOCATE(GRAD_BDY8)
      IF (ALLOCATED(SRCH_SRF8))     DEALLOCATE(SRCH_SRF8)
      IF (ALLOCATED(SRCH_BDY8))     DEALLOCATE(SRCH_BDY8) 
      IF (ASSOCIATED(inv%GRAD_SRF)) DEALLOCATE(inv%GRAD_SRF)
      IF (ASSOCIATED(inv%GRAD_BDY)) DEALLOCATE(inv%GRAD_BDY) 

      IF (ASSOCIATED(inv%MASKG))      DEALLOCATE(inv%MASKG)  
      IF (ASSOCIATED(inv%IGPART))     DEALLOCATE(inv%IGPART)
      IF (ASSOCIATED(inv%MCONN))      DEALLOCATE(inv%MCONN) 
      IF (ASSOCIATED(inv%WMASK))      DEALLOCATE(inv%WMASK) 
      IF (ASSOCIATED(inv%DX))         DEALLOCATE(inv%DX)
      IF (ASSOCIATED(inv%MYGRAD))     DEALLOCATE(inv%MYGRAD)
      IF (ASSOCIATED(inv%ICSC_FDIST)) DEALLOCATE(inv%ICSC_FDIST)
      IF (ASSOCIATED(inv%JCSC_FDIST)) DEALLOCATE(inv%JCSC_FDIST) 

      IF (ASSOCIATED(mid%IRN_LOC)) DEALLOCATE(mid%IRN_LOC)
      IF (ASSOCIATED(mid%JCN_LOC)) DEALLOCATE(mid%JCN_LOC)
      IF (ASSOCIATED(mid%A_LOC))   DEALLOCATE(mid%A_LOC) 
      IF (MYNID == MASTER) THEN
         IF (ASSOCIATED(mid%RHS)) DEALLOCATE(mid%RHS) 
      ENDIF

!
!.... finish MUMPS
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
 9405    FORMAT(' xjoint: Inversion time:',F8.2,' hours')
      ENDIF
      CALL MPI_FINALIZE(MPIERR)
!
!.... some format statements
 8160 FORMAT(' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
             ' +  xjoint: Surface wave step length itereration:',I3,'     +',/,&
             ' +          Surface wave function evaluations:',I3,'        +',/,&
             ' +          Surface wave step length:',E12.4,'        +',/,&
             ' +          Objective function:',E12.4,'              +',/,&
             ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/)
 8161 FORMAT(' ++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
             ' +  xjoint: Body wave step length itereration:',I3,'     +',/,&
             ' +          Body wave function evaluations:',I3,'        +',/,&
             ' +          Body wave step length:',E12.4,'        +',/,&
             ' +          Body wave objective function:',E12.4,' +',/,&
             ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++',/)
 8162 FORMAT(' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
             ' +  xjoint: Step length iteration:',I3,'                     +',/, &
             ' +          Function evaluations:',I3,'                      +',/, &
             ' +          Step length:',E12.4,'                      +',/, &
             ' +          Objective function:',E12.4,'               +',/, &
             ' +          Surface wave objective function:',E12.4,'  +',/, &
             ' +          Body wave objective function:',E12.4,'     +',/, &
             ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/)

 8151 FORMAT(' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
             ' + xjoint: Surface wave iteration:',I3,' finished          +',/,&
             ' +         Function evaluations:',I3,'                     +',/,&
             ' +         Surface wave function evaluations:',I3,'        + '/,&
             ' +         Surface wave step length:',E12.4,'        +',/, &
             ' +         Surface wave objective function:',E12.4,' +',/,&
             ' +         Percent reduction:',F8.3,'                   +',/,&
             ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/)

 8152 FORMAT(' +++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
             ' + xjoint: Body wave iteration:',I3,' finished          +',/,&
             ' +         Function evaluations:',I3,'                  +',/,&
             ' +         Body wave function evaluations:',I3,'        +',/,&
             ' +         Body wave step length:',E12.4,'        +',/, &
             ' +         Body wave objective function:',E12.4,' +',/,&
             ' +         Percent reduction:',F8.3,'                +',/,&
             ' +++++++++++++++++++++++++++++++++++++++++++++++++++++',/)

 8153 FORMAT(' +++++++++++++++++++++++++++++++++++++++++++++++++++++',/, &
             ' + xjoint: Joint iteration:',I3,' finished              +',/,&
             ' +         Function evaluations:',I3,'                  +',/,&
             ' +         Joint function evaluations:',I3,'            +',/,&
             ' +         Joint step length:',E12.4,'            +',/, &
             ' +         Joint objective function:',E12.4,'     +',/,&
             ' +         Percent reduction:',F8.3,'                +',/,&
             ' +++++++++++++++++++++++++++++++++++++++++++++++++++++',/)


      STOP
      END
