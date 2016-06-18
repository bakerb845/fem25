!
!     Updates the source time function and/or receiver response functions.    
!     - B. Baker March Madness 2013
      IMPLICIT NONE
!.... MUMPS
      INCLUDE 'cmumps_struc.h'
      TYPE (CMUMPS_STRUC) MID
!.... MPI
      INCLUDE 'mpif.h'
      REAL*8 TSSIM,TESIM
      INTEGER*4 MASTER,MYID,MYNID, MYSLV_COMM,MYHD_COMM,IPGROUP,JPGROUP,MPIERR,  &
                NPGROUPS,NPROCS,NPARTS
      PARAMETER(MASTER = 0)
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO)  MSH  !mesh, model parameters 
      TYPE (MOD1D_INFO) M1D  !1D models
      TYPE (SRC_INFO)   SRC  !source information
      TYPE (RECV_INFO)  RCV  !receiver information
      TYPE (FRQ_INFO)   FRQ  !frequency information 
      TYPE (WIN_INFO)   WIN  !window information
      TYPE (INV_INFO)   INV  !inversion information
!.... local varaiables
      CHARACTER(80) PROJNM, FILENM, TMPDIR
      COMPLEX*8, ALLOCATABLE :: EST_WRK(:,:,:,:), EST_ROT(:,:,:,:), &
                                OBS_WRK(:,:,:,:), SRCWRK(:,:), &
                                RECWRK(:,:,:), UE(:), WAVE(:) 
      REAL*8, ALLOCATABLE :: XREC(:), ZREC(:), FREQ_LOC(:) 
      REAL*4, ALLOCATABLE :: WGHT_WRK(:,:,:,:) 
      INTEGER*4, ALLOCATABLE :: NSG_LSTB(:), NSG_LST(:), ISRCLST(:) 
      LOGICAL*4, ALLOCATABLE :: LFSURF(:)
      COMPLEX*8 CONE, CZERO
      REAL*8 PYAVG
      REAL*4 FOBJ, FOBJ_SRF, FOBJ_BDY, COBJ4
      INTEGER*4 STAT(MPI_STATUS_SIZE), NABS, NELEME, NSGMAX, IFREQ, JFREQ,           &  
                IFREQL, ISRC, JSRC, KSRC, ISG, JOB, IBLOCK, NBLOCKS, NDIM_HD,        &
                NFREQ_SRF_HD, NFREQ_BDY_HD, NREC_HD, NSRC_SRF_HD, NSRC_BDY_HD,       &
                IMODREC, IMODSRC, IERR, MYDEST, &
                MYTAG, MYSRC, NSG_GRP
      LOGICAL*4 LFILES, LSRCEX, LRECEX, LCGRNS, LNSRF
      PARAMETER(CONE = CMPLX(1.0,0.0), CZERO = CMPLX(0.0,0.0))
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
!----------------------------------------------------------------------------------------!
!                    Read information on mesh and model                                  !
!----------------------------------------------------------------------------------------!
!
!.... have head node read model information
      TSSIM = MPI_WTIME() 
      IERR = 0 
      LFILES = .FALSE. 
      LCGRNS = .false. !.TRUE.  !always calculate new greens fns beginning STF/RRF udate
      IF (MYID == MASTER) THEN
         WRITE(*,9407) 
 9407    FORMAT(' ---------------------------------------------------------------',/, &
                ' -   xsrcrec25: A massively parallel program for estimating    -',/, &
                ' -              source and response functions for plane and    -',/, &
                ' -              surface waves in locally hetergeneous media    -',/, &
                ' ---------------------------------------------------------------',/) 
         PROJNM(1:80) = ' ' 
         WRITE(*,*) 'xsrcrec25: Enter project name:'
         READ(*,'(A)') PROJNM

         PROJNM = ADJUSTL(PROJNM)    
         WRITE(*,*) 'xsrcrec25: Enter number of process groups:'
         READ *, NPGROUPS
         msh%freq0 = 20.d0
         msh%aztol = 10.d0
         msh%aoitol = 5.d0
         inv%LUNWRAP = .FALSE. !use an unwrapped phase?
!....... get forward modeling parameters
         CALL READ_FWD_INI(PROJNM,.TRUE., TMPDIR,LFILES,LNSRF, MSH, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xbielak25: Error I cannot locate your spec file!'
            GOTO 500 
         ENDIF
!....... initial warnings on NPGROUPS
         IF (MOD(NPROCS,NPGROUPS) /= 0) THEN
            WRITE(*,*) 'xsrcrec25: Error I cant divide nprocs by npgroups evenly!'
            IERR = 1
            GOTO 500
         ENDIF
!....... also read the job
         WRITE(*,*) 'xsrcrec25: Enter job type'
         WRITE(*,*) '           (1) is for a source time function estimate (default)'
         WRITE(*,*) '           (2) is for a receiver response estimate'
         WRITE(*,*) '           (3) is for both' 
         READ *, JOB
         IF (JOB < 0 .OR. JOB > 3) THEN
            WRITE(*,*) 'xsrcrec25: Overriding to job = 1'
            JOB = 1
         ENDIF
         WRITE(*,*) 'xsrcrec25: Enter inversion block number'
         READ *, IBLOCK 
         IMODSRC = 3 !1 -> phase only, 2 -> amplitude only, 3 -> both (default)
         IMODREC = 3 !1 -> phase only, 2 -> amplitude only, 3 -> both (default) 
         IF (JOB == 1 .OR. JOB == 3) THEN
            WRITE(*,*) 'xsrcrec25: Enter source inversion type'
            WRITE(*,*) '           (1) Phase only'
            WRITE(*,*) '           (2) Amplitude only'
            WRITE(*,*) '           (3) Phase and amplitude (default)'
            READ *, IMODSRC
         ENDIF
         IF (IMODSRC < 0 .OR. IMODSRC > 3) THEN
            WRITE(*,*) 'xsrcrec25: Overriding IMODSRC to 3' 
            IMODSRC = 3
         ENDIF 
         IF (JOB == 2 .OR. JOB == 3) THEN
            WRITE(*,*) 'xsrcrec25: Enter receiver inversion type'
            WRITE(*,*) '           (1) Phase only'
            WRITE(*,*) '           (2) Amplitude only'
            WRITE(*,*) '           (3) Phase and amplitude (default)'
            READ *, IMODREC 
         ENDIF
         IF (IMODREC < 0 .OR. IMODREC > 3) THEN
            WRITE(*,*) 'xsrcrec25: Overriding IMODREC to 3'
            IMODREC = 3 
         ENDIF
!....... read the mesh
         NPARTS = NPROCS/NPGROUPS
         WRITE(*,*) 'xsrcrec25: Reading mesh...'
         CALL RDMESHBK(PROJNM, NELEME,NABS, MSH, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xsrcrec25: Error reading mesh files!'
            GOTO 500
         ENDIF

         WRITE(*,*) 'xsrcrec25: Setting interpolation points...'
         ALLOCATE(msh%XIPTS (msh%NLXI))
         ALLOCATE(msh%ETAPTS(msh%NLETA))
         print *, msh%iitype
         CALL GENPTS(msh%IITYPE, msh%NLXI,msh%NLETA, msh%XIPTS,msh%ETAPTS)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xsrcrec25: Error reading mesh files!'
            GOTO 500
         ENDIF
!
!....... initial check on observation file 
         WRITE(*,*) 'xsrcrec25: Checking headers on observation file...'
         !CALL RDTOBS_HD(PROJNM, NDIM_HD,NFREQ_HD,NREC_HD,NSRC_HD, IERR)
         CALL RDTOBS_SB_HD(PROJNM, NDIM_HD,NFREQ_SRF_HD,NFREQ_BDY_HD,  & 
                           NREC_HD,NSRC_SRF_HD,NSRC_BDY_HD, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xsrcrec25: There was an error with your observation file'
            GOTO 500 
         ENDIF 
         IF (NDIM_HD /= NDIM) THEN
            WRITE(*,*) 'xsrcrec25: Error component number mismatch!',NDIM_HD,NDIM
            IERR = 1 
            GOTO 500 
         ENDIF
         IF (NFREQ_SRF_HD + NFREQ_BDY_HD <= 0) THEN
            WRITE(*,*) 'xsrcrec25: Error observations have 0 frequencies!', &
                       NFREQ_SRF_HD, NFREQ_BDY_HD 
            IERR = 1 
            GOTO 500 
         ENDIF
!
!....... check desired blocks fits in the number of blocks 
         !CALL RDFREQI_BLHD(PROJNM,NBLOCKS,IERR)
         CALL RD_JFREQ_INV_HD(PROJNM,NBLOCKS,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xsrcrec25: Error cannot locate the frequency block file'
            GOTO 500
         ENDIF
         IF (IBLOCK > NBLOCKS) THEN
            WRITE(*,*) 'xsrcrec25: Error IBLOCK > NBLOCKS!'
            IERR = 1
            GOTO 500
         ENDIF 
!
!....... read the source list
         CALL RDSRC_EQ(PROJNM, msh%XLATMIN,msh%XLONMIN,msh%XLATMAX,msh%XLONMAX, &
                       msh%AZMOD, SRC,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xsrcrec25: Error calling rdsrc_eq!'
            GOTO 500
         ENDIF
         IF (NSRC_SRF_HD /= src%NSRC_SRF) THEN
            WRITE(*,*) 'xsrcrec25: Surface wave source number mismatch!', &
                       NSRC_SRF_HD,src%NSRC_SRF
            IERR = 1
            GOTO 500
         ENDIF
         IF (NSRC_BDY_HD /= src%NSRC_BDY) THEN
            WRITE(*,*) 'xsrcrec25: Body wave source number mismatch!', &
                       NSRC_BDY_HD,src%NSRC_BDY
            GOTO 500
         ENDIF
!
!....... check for a receiver file
         CALL RDREC_HD(PROJNM, rcv%NREC,IERR)
         IF (IERR /= 0) THEN
            IF (IERR > 0) THEN
               WRITE(*,*) 'xsrcrec25: Error reading recv header'
               GOTO 500
            ELSE
               WRITE(*,*) 'xsrcrec25: There is no receiver file!  Cannot make estimates!'
               IERR = 1
               GOTO 500
            ENDIF
         ENDIF
         IF (rcv%NREC <= 0) THEN
            WRITE(*,*) 'xsrcrec25: There are no receivers!  Aborting!'
            GOTO 500
         ENDIF
         IF (NREC_HD /= rcv%NREC) THEN
            WRITE(*,*) 'xsrcrec25: Error recevier number mismatch!',NREC_HD,rcv%NREC
            IERR = 1
            GOTO 500
         ENDIF
         ALLOCATE(  LFSURF(rcv%NREC))
         ALLOCATE(    XREC(rcv%NREC))
         ALLOCATE(rcv%YREC(rcv%NREC))
         ALLOCATE(    ZREC(rcv%NREC))
         CALL RDREC(PROJNM, rcv%NREC, LFSURF,XREC,rcv%YREC,ZREC, IERR)
!
!....... check if there are windows
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
!
!....... read the frequency block
         WRITE(*,*) 'xsrcrec25: Reading frequency block...'
         !CALL RDFREQI(PROJNM,IBLOCK, FRQ,IERR) 
         CALL RD_JFREQ_INV(PROJNM,IBLOCK,win%LWNDO_SRF,win%LWNDO_BDY,  &
                           FRQ,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xsrcrec25: Cannot read inversion frequency block!'
            GOTO 500
         ENDIF
!
!....... get the observations
         WRITE(*,*) 'xsrcrec25: Reading observation file...'
         ALLOCATE(  inv%OBS(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
         ALLOCATE(inv%WGHTS(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
         inv%OBS(1:NDIM,1:frq%NFREQ,1:rcv%NREC,1:src%NSRC) = CZERO
         inv%WGHTS(1:NDIM,1:frq%NFREQ,1:rcv%NREC,1:src%NSRC) = 0.0
!        CALL RDTOBS25(PROJNM, NDIM,frq%NFREQ,rcv%NREC,frq%NFREQ, NDIM,rcv%NREC,    &   
!                      src%NSRC, LUNWRAP, frq%FREQ, WGHTS,OBS, IERR)
         CALL RDTOBS_SB25(PROJNM, RCV,SRC,FRQ,INV, IERR)
         IF (IERR /= 0) THEN
            DEALLOCATE(inv%OBS)
            DEALLOCATE(inv%WGHTS) 
            WRITE(*,*) 'xsrcrec25: Error reading observation file1'
            GOTO 500 
         ENDIF
!
!....... check for a receiver response file
         ALLOCATE(rcv%RECV(NDIM,frq%NFREQ,rcv%NREC)) 
         IF (JOB == 1) THEN !okay for RRFs 
            CALL RDREC_RESP(PROJNM, NDIM,frq%NFREQ, NDIM,frq%NFREQ,rcv%NREC, msh%AZMOD, &
                            frq%FREQ, LRECEX,rcv%RECV)
            IF (.NOT.LRECEX) THEN
               WRITE(*,*) 'xlbfgs_hsl: Receiver response functions are set to unity'
               IERR = 0
            ENDIF
         ELSE
            rcv%RECV(1:NDIM,1:frq%NFREQ,1:rcv%NREC) = CONE 
         ENDIF
         ALLOCATE(src%SOURCE(frq%NFREQ,src%NSRC))
         IF (JOB == 2) THEN !okay for STFs 
            CALL RDSRCINV(frq%NFREQ,PROJNM, frq%NFREQ,src%NSRC, frq%FREQ, &
                          LSRCEX,src%SOURCE)
            IF (.NOT.LSRCEX) THEN
               WRITE(*,*) 'xsrcrec25: STFs were set to unity'
            ENDIF
         ELSE 
            WRITE(*,*) 'xsrcrec25: Setting STFs to unity...'
            src%SOURCE(1:frq%NFREQ,1:src%NSRC) = CONE 
         ENDIF 
!
!....... fill the 1D models
         CALL FILL1D(MSH,M1D)
         WRITE(*,*) 'xsrcrec25: Splitting process groups...'
      ENDIF
      CALL MKMPI_GROUPS(MASTER,MYID,NPGROUPS,NPROCS,MPI_COMM_WORLD, &
                        IPGROUP,MYNID,MYHD_COMM,MYSLV_COMM)
!----------------------------------------------------------------------------------------!
!     MUMPS initialization phase and graph reordering utilities.                         !
!----------------------------------------------------------------------------------------!
      IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Initializing MUMPS...'
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
         WRITE(*,*) 'xsrcrec25: Generating graph...'
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
         WRITE(*,*)
         WRITE(*,9408) msh%NORD,msh%NELEM,NABS,NELEME, &
                       msh%NNPG,msh%NNPE,msh%NDOF,msh%NZERO
 9408    FORMAT(' xsrcrec25: Polynomial order:'                    ,I4 ,/,        &
                '            Number of elements in mesh:'          ,I10,/,         &
                '            Number of absorbing elements:'        ,I8 ,/,         &
                '            Number of Bielak elements:'           ,I8 ,/,         &
                '            Number of anchor nodes in mesh:'      ,I10,/,         &
                '            Number of nodes in Bielak boundary:'  ,I10,/,         &
                '            Number of degrees of freedom:'        ,I14,/,         &
                '            Number of non-zeros in global matrix:',I16,/)
!
!....... free some unwatend variables 
         IF (ALLOCATED(LFSURF)) DEALLOCATE(LFSURF) 
         IF (ALLOCATED(ZREC))   DEALLOCATE(ZREC)
         mid%ICNTL(18) = 3 !i will tell MUMPS how to split the mesh
      ENDIF
!----------------------------------------------------------------------------------------!
!                          Broadcast Forward Problem Parameters                          !
!----------------------------------------------------------------------------------------!
      IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Broadcasting 2D models...'
      CALL MPI_BCAST(NPARTS,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR) 
      CALL BCAST_MESH_INFO(MYID,MPI_COMM_WORLD,MASTER, &
                          LNSRF, PROJNM,TMPDIR, MSH) 
      IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Broadcasting 1D models...'
      CALL BCAST_MOD1D_INFO(MYID,MPI_COMM_WORLD,MASTER, MSH,M1D) 
      IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Broadcasting source details...' 
      CALL BCAST_SRC_INFO(MYID,MPI_COMM_WORLD,MASTER, SRC)
      IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Broadcasting frequency information...'
      CALL BCAST_FRQ_JOINT_INFO(MYID,MPI_COMM_WORLD,MASTER, FRQ)
      CALL MPI_BCAST(inv%LUNWRAP,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR) 
      IF (MYNID == MASTER) THEN
         IF (MYID == MASTER) &
         WRITE(*,*) 'xsrcrec25: Broadcasting receiver response functions...'
         CALL BCAST_RCV_INFO(MYID,MYHD_COMM,MASTER, RCV)
         CALL BCAST_RCV_RESP(MYID,MYHD_COMM,MASTER, .TRUE.,frq%NFREQ,RCV)
         IF (MYID == MASTER) &
         WRITE(*,*) 'xsrcrec25: Broadcasting source time function...'
         CALL BCAST_SRC_STF(MYID,MYHD_COMM,MASTER,  .TRUE.,frq%NFREQ, SRC)
         IF (inv%LUNWRAP) THEN
            IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Broadcasting free surface DOFS...'
            CALL BCAST_FS_INFO(MYID,MYHD_COMM,MASTER, MSH)  
         ENDIF 
      ENDIF
      CALL MPI_BCAST(win%LWNDO_SRF,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(win%LWNDO_BDY,1,MPI_LOGICAL, MASTER,MPI_COMM_WORLD,MPIERR)

!----------------------------------------------------------------------------------------!
!                            Calculate 1D models solutions                               !
!----------------------------------------------------------------------------------------!
      IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Calculating 1D solutions...'
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
      ALLOCATE(src%PYTAB(frq%NFREQ,src%NSRC))
      src%PYTAB(1:frq%NFREQ,1:src%NSRC) = 0.D0
      !CALL GEN_GRNS(MYID,MASTER,MPI_COMM_WORLD, TMPDIR,LCGRNS,LNSRF, &
      !              SRC,M1D,FRQ,MSH, IERR)
      CALL GEN_GRNS_SB(MYID,MASTER,MPI_COMM_WORLD, TMPDIR,LCGRNS,LNSRF, &
                       SRC,M1D,FRQ,MSH, IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'xsrcrec25: Error calling gen_grns_sb!'
         GOTO 500 
      ENDIF
!.... clean some space 
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
!.... count number of local non-zeros
      IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Generating local graphs...'
      CALL ICNDZ_LOC(MYNID,1,msh%NDOF, msh%PART,msh%IRPTR, msh%NDOFL,msh%NZLOC)
!.... generate distributed CSR structure
      ALLOCATE(msh%MYDOFS(msh%NDOFL))
      ALLOCATE(msh%IRPTR_LOC(msh%NDOFL+1))
      ALLOCATE(msh%JCPTR_LOC(msh%NZLOC))
      CALL PART2CRSL(msh%NZERO,msh%NZLOC, msh%NDOF,msh%NDOFL, MYNID,1,             &   
                     msh%IRPTR,msh%JCPTR,msh%PART, msh%MYDOFS,msh%IRPTR_LOC,msh%JCPTR_LOC)
!.... counter number of local elements
      msh%NELEML = ICELEML(MYNID, NDIM,msh%NEN,msh%NDOF, msh%NELEM, 1,             &   
                           msh%NEN,NDIM, msh%PART,msh%LM)
!.... have process ave its elements 
      ALLOCATE(msh%MYELEM(msh%NELEML))
      CALL GENELEML(MYNID, NDIM,msh%NEN,msh%NDOF, msh%NELEM,msh%NELEML, 1,         &   
                    msh%NEN,NDIM, msh%PART,MSH%LM, msh%MYELEM)
      IF (ASSOCIATED(msh%IRPTR)) DEALLOCATE(msh%IRPTR)
      IF (ASSOCIATED(msh%JCPTR)) DEALLOCATE(msh%JCPTR)
      IF (ASSOCIATED(msh%PART))  DEALLOCATE(msh%PART)
!.... create the distributed COO storage for MUMPS
      mid%NZ_LOC = msh%NZLOC
      ALLOCATE(mid%IRN_LOC(mid%NZ_LOC))
      ALLOCATE(mid%JCN_LOC(mid%NZ_LOC)) 
      CALL CRS2COOLOC(msh%NZLOC,msh%NDOFL, msh%MYDOFS,msh%IRPTR_LOC,msh%JCPTR_LOC,      &
                      mid%IRN_LOC,mid%JCN_LOC)
!.... reorder the matrix
      IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Reordering matrix...'
      mid%JOB = 1 
      CALL CMUMPS(MID) 
      ALLOCATE(mid%A_LOC(mid%NZ_LOC))
      IF (MYNID == MASTER) THEN
         mid%LRHS = mid%N
         ALLOCATE(mid%RHS(mid%N))
         ALLOCATE(inv%EST(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
         IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Nulling estimates...'
         inv%EST(1:NDIM,1:frq%NFREQ,1:rcv%NREC,1:src%NSRC) = CZERO
      ENDIF
      IF (MYID == MASTER) WRITE(*,*)
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
!
!.... set space
      ALLOCATE(UE(msh%NDOF))
      IF (MYNID == MASTER) THEN
         ALLOCATE(WAVE(msh%NDOF))
         ALLOCATE(NSG_LSTB(NPGROUPS))
         ALLOCATE(NSG_LST(NPGROUPS))
         ALLOCATE(ISRCLST(src%NSRC)) 
      ELSE
         ALLOCATE(WAVE(1))
      ENDIF
!----------------------------------------------------------------------------------------!
!                 Loop on frequencies                                                    !
!----------------------------------------------------------------------------------------!
  !goto 8
      DO 1000 IFREQL=1,frq%NFREQ
!
!....... get frequency
         IFREQ = (IFREQL - 1)*NPGROUPS + IPGROUP + 1
         IF ( (IFREQL - 1)*NPGROUPS + 1 > frq%NFREQ) GOTO 1000 !done
!
!....... set and send source group info
         src%NSG = 0
         IF (MYNID == MASTER) THEN
            IF (IFREQ <= frq%NFREQ) THEN
               IF (frq%CFTYPE(IFREQ) == 'S') THEN
                  IF (win%LWNDO_SRF .OR. frq%LINVF(IFREQ)) THEN
                     CALL FILL_SRCPRM_SB(IFREQ, .TRUE., m1d%VFAST_SRF,msh%AZTOL,  &
                                         msh%AOITOL, SRC)
                  ENDIF
               ELSEIF (frq%CFTYPE(IFREQ) == 'B') THEN 
                  IF (win%LWNDO_BDY .OR. frq%LINVF(IFREQ)) THEN
                     CALL FILL_SRCPRM_SB(IFREQ,.FALSE., m1d%VFAST_BDY,msh%AZTOL,  &
                                         msh%AOITOL, SRC)
                  ENDIF
               ELSE
                  WRITE(*,*) 'xsrcrec25: Cant determine frequency type!'
                  IERR = 1
                  GOTO 500
               ENDIF
            ELSE
               src%NSG = 0
            ENDIF
            NSG_LSTB(1:NPGROUPS) = 0 
            NSG_LSTB(IPGROUP+1) = src%NSG
            CALL MPI_ALLREDUCE(src%NSG ,NSGMAX,        1,MPI_INTEGER,MPI_MAX, &
                               MYHD_COMM,MPIERR)
            CALL MPI_ALLREDUCE(NSG_LSTB,NSG_LST,NPGROUPS,MPI_INTEGER,MPI_SUM, &
                               MYHD_COMM,MPIERR)
         ENDIF
         CALL MPI_BCAST(NSGMAX ,1,MPI_INTEGER,MASTER, MYSLV_COMM,MPIERR)
         IF (IFREQ <= frq%NFREQ) CALL BCAST_SRC_PTRS(MYNID,MYSLV_COMM,MASTER, SRC)
!
!....... loop on source groups
         DO 2000 ISG=1,NSGMAX
            IF (MYID == MASTER) THEN
               MYSRC = 0
               PYAVG = 0.D0
               DO 201 JPGROUP=0,NPGROUPS-1
                  JFREQ = (IFREQL - 1)*NPGROUPS + JPGROUP + 1
                  NSG_GRP = NSG_LST(JPGROUP+1)
                  IF (JFREQ <= frq%NFREQ .AND. NSG_GRP > 0 .AND. ISG <= NSG_GRP) THEN 
                     IF (frq%CFTYPE(JFREQ) == 'S') THEN 
                        IF (frq%LINVF(JFREQ) .OR. win%LWNDO_SRF) THEN 
                           IF (JPGROUP > 0) THEN
                              MYSRC = JPGROUP*NPROCS/NPGROUPS
                              CALL MPI_RECV(PYAVG,1,MPI_DOUBLE_PRECISION,MYSRC,     &
                                            MPI_ANY_TAG, MPI_COMM_WORLD,STAT,MPIERR)
                           ELSE
                              PYAVG = src%PYAVG(ISG)
                           ENDIF
                           WRITE(*,9410) JPGROUP+1,ISG,PYAVG,frq%FREQ(JFREQ)
                        ENDIF
                     ELSE 
                        IF (frq%LINVF(JFREQ) .OR. win%LWNDO_BDY) THEN
                           IF (JPGROUP > 0) THEN
                              MYSRC = JPGROUP*NPROCS/NPGROUPS
                              CALL MPI_RECV(PYAVG,1,MPI_DOUBLE_PRECISION,MYSRC,     &
                                            MPI_ANY_TAG, MPI_COMM_WORLD,STAT,MPIERR)
                           ELSE
                              PYAVG = src%PYAVG(ISG)
                           ENDIF
                           WRITE(*,9411) JPGROUP+1,ISG,PYAVG,frq%FREQ(JFREQ) 
                        ENDIF
                     ENDIF !end check on frequency type 
                  ENDIF !end check on ranges
  201          CONTINUE !loop on process groups
            ELSE !slave
               IF (MYNID == MASTER) THEN !process group master
                  MYDEST = MASTER
                  MYTAG  = MYID
                  IF (IFREQ <= frq%NFREQ .AND. src%NSG > 0 .AND. ISG <= src%NSG) THEN
                     IF (frq%CFTYPE(IFREQ) == 'S') THEN
                        IF (frq%LINVF(IFREQ) .OR. win%LWNDO_SRF) THEN
                           CALL MPI_SEND(src%PYAVG(ISG),1,MPI_DOUBLE_PRECISION,MASTER, &
                                         MYTAG, MPI_COMM_WORLD,MPIERR)
                        ENDIF
                     ELSE !body wave
                        IF (frq%LINVF(IFREQ) .OR. win%LWNDO_BDY) THEN
                           CALL MPI_SEND(src%PYAVG(ISG),1,MPI_DOUBLE_PRECISION,MASTER, &
                                         MYTAG, MPI_COMM_WORLD,MPIERR)

                        ENDIF
                     ENDIF !end check on frequency type 
                  ENDIF !end check on range
               ENDIF !end check on mynid
            ENDIF
!
!.......... assemble impedance matrix for average py value in group 
            IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Assembling impedance matrix...'
            IF (ISG <= src%NSG .AND. IFREQ <= frq%NFREQ) THEN
               CALL ASMBLE_DRIVER(msh%nzloc,frq%CFTYPE(IFREQ),frq%FREQ(IFREQ), &
                                  src%PYAVG(ISG), MSH, mid%A_LOC,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xsrcrec25: An error occurred calling asmble_driver',MYID
                  GOTO 500
               ENDIF
            ENDIF
            IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Factoring impedance matrix...'
            IF (ISG <= src%NSG .AND. IFREQ <= frq%NFREQ) THEN
               mid%JOB = 2 !factorization
               CALL CMUMPS(MID)
               IF (mid%INFO(1) < 0) THEN
                  WRITE(*,*) 'xsrcrec25: An error occurred in the factorization!',MYID
                  IERR = 1
                  GOTO 500
               ENDIF
            ENDIF
!
!.......... have master tell us which sources are being used
            IF (MYID == MASTER) THEN
               MYSRC = 0
               ISRCLST(1:src%NSRC) = 0
               DO 290 JPGROUP=0,NPGROUPS-1
                  JFREQ = (IFREQL - 1)*NPGROUPS + JPGROUP + 1
                  NSG_GRP = NSG_LST(JPGROUP+1)
                  IF (JFREQ <= frq%NFREQ .AND. NSG_GRP > 0 .AND. ISG <= NSG_GRP) THEN
                     IF (frq%CFTYPE(JFREQ) == 'S') THEN
                        IF (frq%LINVF(JFREQ) .OR. win%LWNDO_SRF) THEN
                           IF (JPGROUP > 0) THEN !receiving
                              MYSRC = JPGROUP*NPROCS/NPGROUPS
                              CALL MPI_RECV(ISRCLST,src%NSRC,MPI_INTEGER,               &
                                            MYSRC,MPI_ANY_TAG,MPI_COMM_WORLD,STAT,MPIERR)
                           ELSE !mine
                              KSRC = 0
                              DO 291 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1) - 1
                                 KSRC = KSRC + 1
                                 ISRCLST(KSRC) = src%ISRCPRM(JSRC)
  291                         CONTINUE
                           ENDIF
                           DO 292 KSRC=1,src%NSRC
                              ISRC = ISRCLST(KSRC)
                              IF (ISRC == 0) GOTO 282
                              WRITE(*,9412) ISRC, src%AOI(ISRC), src%BAZN(ISRC), &
                                            src%PYTAB(JFREQ,ISRC)
  292                      CONTINUE 
  282                      CONTINUE !break ahead
                        ENDIF
                     ELSE !body wave frequency
                        IF (frq%LINVF(JFREQ) .OR. win%LWNDO_BDY) THEN
                           IF (JPGROUP > 0) THEN
                              MYSRC = JPGROUP*NPROCS/NPGROUPS
                              CALL MPI_RECV(ISRCLST,src%NSRC,MPI_INTEGER,                &
                                            MYSRC,MPI_ANY_TAG,MPI_COMM_WORLD,STAT,MPIERR)
                           ELSE !mine
                              KSRC = 0
                              DO 293 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1) - 1
                                 KSRC = KSRC + 1
                                 ISRCLST(KSRC) = src%ISRCPRM(JSRC)
  293                         CONTINUE
                           ENDIF !end check on source
                           DO 294 KSRC=1,src%NSRC
                              ISRC = ISRCLST(KSRC)
                              IF (ISRC == 0) GOTO 284
                              WRITE(*,9413) ISRC, src%AOI(ISRC), src%BAZN(ISRC), &
                                            src%PYTAB(JFREQ,ISRC)
  294                      CONTINUE 
  284                      CONTINUE !break ahead 
                        ENDIF !end check on working
                     ENDIF !end check on frequency type
                  ENDIF !end check on bounds
  290          CONTINUE !loop on process groups
            ELSE !slave
               IF (MYNID == MASTER) THEN
                  MYDEST = MASTER
                  MYTAG = MYID
                  ISRCLST(1:src%NSRC) = 0
                  IF (IFREQ <= frq%NFREQ .AND. src%NSG > 0 .AND. ISG <= src%NSG) THEN
                     IF (frq%CFTYPE(IFREQ) == 'S') THEN !surface wave
                        IF (frq%LINVF(IFREQ) .OR. win%LWNDO_SRF) THEN
                           KSRC = 0
                           DO 297 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1)-1
                              KSRC = KSRC + 1
                              ISRCLST(KSRC) = src%ISRCPRM(JSRC)
  297                      CONTINUE
                           CALL MPI_SEND(ISRCLST,src%NSRC,MPI_INTEGER, MASTER,MYTAG, &
                                         MPI_COMM_WORLD,MPIERR)

                        ENDIF !end check on working
                     ELSE !body wave
                        IF (frq%LINVF(IFREQ) .OR. win%LWNDO_BDY) THEN
                           KSRC = 0
                           DO 298 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1)-1
                              KSRC = KSRC + 1
                              ISRCLST(KSRC) = src%ISRCPRM(JSRC)
  298                      CONTINUE !loop on sources
                           CALL MPI_SEND(ISRCLST,src%NSRC,MPI_INTEGER, MASTER,MYTAG, &
                                         MPI_COMM_WORLD,MPIERR)
                        ENDIF !end check on working
                     ENDIF !end check on frequency type
                  ENDIF !end check on bounds
               ENDIF !end check on mynid
            ENDIF !end check on myid
            IF (ISG > src%NSG)     GOTO 2001 !no solves, leave
            IF (IFREQ > frq%NFREQ) GOTO 2001 !no solves, leave
!----------------------------------------------------------------------------------------!
!                                 Loop on sources in source group                        ! 
!----------------------------------------------------------------------------------------!
!
!.......... loop on sources in source group
            DO 3000 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1)-1
               ISRC = src%ISRCPRM(JSRC) !extract original source ID
!
!............. generate 1D solution
               IF (MYNID == MASTER) THEN
                  CALL LOAD_GRNS25(NDIM,msh%NDOF,msh%NNPE,NDIM, src%SRCTYP(ISRC),ISRC,   &
                                   frq%FREQ(IFREQ), src%SOURCE(IFREQ,ISRC),msh%IDOFSE,   &
                                   UE,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xsrcrec25: Error in load_grns25!'
                     GOTO 500
                  ENDIF
               ENDIF  !end check on master
               CALL MPI_BCAST(UE,msh%NDOF,MPI_COMPLEX, MASTER,MYSLV_COMM,MPIERR)
               CALL ASMB_PEFF(MASTER,MYSLV_COMM, msh%NDOF,msh%NDOFL,msh%NZLOC, msh%CNP, &
                              msh%MYDOFS,msh%IRPTR_LOC,msh%JCPTR_LOC, mid%A_LOC,UE, WAVE)
               IF (MYNID == MASTER) mid%RHS(1:msh%NDOF) = WAVE(1:msh%NDOF)
               mid%JOB = 3 !solve phase
               CALL CMUMPS(MID)
               IF (mid%INFO(1) < 0) THEN
                  WRITE(*,*) 'xsrcrec25: An error occurred in the solve phase!',MYID
                  IERR = 1
                  GOTO 500
               ENDIF
!
!............. head processes extract solution
               IF (MYNID == MASTER) THEN
                  CALL ADDBLK(msh%NDOF,msh%CNP,UE, mid%RHS) !add background field in
                  IF (inv%LUNWRAP) THEN
                     CALL EXRESP_MP_2(NDIM,rcv%NREC,msh%NDOF,  &
                                    1,rcv%NREC,NDIM, frq%FREQ(IFREQ),  &
                                    src%PYTAB(IFREQ,ISRC),rcv%YREC, mid%RHS, rcv%MRDOF,  &
                                    rcv%RECV(1:NDIM,IFREQ,1:rcv%NREC),                   &
                                    inv%EST(1:NDIM,IFREQ,1:rcv%NREC,ISRC))

                     IF (IERR /= 0) THEN
                        WRITE(*,*) 'xsrcrec25: Error unwrapping phase!'
                        GOTO 500
                     ENDIF
                  ELSE
                     CALL EXRESP25D(NDIM,rcv%NREC,msh%NDOF,  &
                                    1,rcv%NREC,NDIM, frq%FREQ(IFREQ),  &
                                    src%PYTAB(IFREQ,ISRC),rcv%YREC, mid%RHS, rcv%MRDOF,  &
                                    rcv%RECV(1:NDIM,IFREQ,1:rcv%NREC),                   &
                                    inv%EST(1:NDIM,IFREQ,1:rcv%NREC,ISRC))
                  ENDIF
               ENDIF !end check on master
 3000       CONTINUE !loop on sources
 2001       CONTINUE !break ahead for loop on source group
 2000    CONTINUE !loop on source groups
!
!....... clean source info
         IF (src%NSG > 0 .AND. IFREQ <= frq%NFREQ) THEN
            IF (ASSOCIATED(src%ISRCPRM)) DEALLOCATE(src%ISRCPRM)
            IF (ASSOCIATED(src%ISGPTR))  DEALLOCATE(src%ISGPTR)
            IF (ASSOCIATED(src%PYAVG))   DEALLOCATE(src%PYAVG)
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
 1000 CONTINUE !Loop on frequencies 
 8 continue
!----------------------------------------------------------------------------------------!
!                     Estimate reduction and STF/RRF updates                             !
!----------------------------------------------------------------------------------------!
!
!.... reduce estimates
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
      IF (MYNID == MASTER) THEN
         IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Reducing estimates...'
         CALL REDEST4(MYID,MASTER,MYHD_COMM,NPGROUPS,IPGROUP,                            &
                      NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC,1,inv%EST)
         IF (MYID == MASTER) THEN
            FOBJ_SRF = 0.0
            FOBJ_BDY = 0.0
            !reality check on objective function
            IF (frq%NFREQ_SRF > 0 .AND. src%NSRC_SRF > 0) THEN
               ALLOCATE( EST_WRK(NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF))
               ALLOCATE( EST_ROT(NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF)) 
               ALLOCATE( OBS_WRK(NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF))
               ALLOCATE(FREQ_LOC(frq%NFREQ_SRF))
               ALLOCATE(WGHT_WRK(NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF))
               CALL GET_OE_SB(.TRUE.,  NDIM,frq%NFREQ_SRF,rcv%NREC, INV,FRQ,SRC,RCV,  &
                              WGHT_WRK, EST_WRK,OBS_WRK, IERR)
               IF (IERR /= 0) THEN 
                  WRITE(*,*) 'xsrcrec25: Error calling get_oe_sb 1' 
                  GOTO 500
               ENDIF
               CALL GET_FREQ_MOD(frq%NFREQ, frq%NFREQ_SRF, .TRUE., frq%CFTYPE,frq%FREQ, &
                                 FREQ_LOC,IERR) 
               IF (IERR /= 0) THEN 
                  WRITE(*,*) 'xsrcrec25: Error calling get_freq_mod' 
                  GOTO 500
               ENDIF
               FOBJ_SRF = COBJ4(NDIM,frq%NFREQ_SRF,rcv%NREC,                &
                                frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF, NDIM,  & 
                               .FALSE.,2,3, .FALSE., 0.2,0.2,0.0,           &
                                FREQ_LOC, WGHT_WRK,OBS_WRK,EST_WRK)
               WRITE(*,*) 'xsrcrec25: Surface wave objective function:',FOBJ_SRF
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xsrcrec25: Error calling get_oe_sb 1'
                  GOTO 500
               ENDIF
               IF (JOB == 1 .OR. JOB == 3) THEN
                  WRITE(*,*) 'xsrcrec25: Solving new surface wave STFs...'
                  ALLOCATE(SRCWRK(frq%NFREQ_SRF,src%NSRC_SRF)) 
                  SRCWRK(:,:) = CZERO
                  CALL GET_SRC_MOD(frq%NFREQ,frq%NFREQ_SRF, frq%NFREQ,src%NSRC, &
                                   frq%NFREQ_SRF,src%NSRC_SRF,                  &
                                   .TRUE.,frq%CFTYPE, src%SRCTYP, src%SOURCE,   &
                                   SRCWRK,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xsrcrec25: Error calling get_src_mod 1'
                     GOTO 500
                  ENDIF
                  CALL SRCUPD(NDIM,frq%NFREQ_SRF,rcv%NREC,                               &
                              NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF,                  &
                              inv%LUNWRAP,inv%IMODSRC, msh%AZMOD, EST_WRK,OBS_WRK, SRCWRK)
               
                  CALL MODSRC(frq%NFREQ_SRF, frq%NFREQ_SRF,src%NSRC_SRF,inv%IRESTP,  &
                              SRCWRK)

                  CALL CONVEST_STF(NDIM,frq%NFREQ_SRF,rcv%NREC,               &
                                   NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF,  &
                                   SRCWRK, EST_WRK)
                  FILENM(1:80) = ' '
                  FILENM = TRIM(ADJUSTL(PROJNM))//'_srf'
                  CALL WTSTF(FILENM, frq%NFREQ_SRF, frq%NFREQ_SRF,src%NSRC_SRF, IBLOCK,0,&
                             FREQ_LOC,SRCWRK, IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xsrcrec25: Error writing STF 1'
                     GOTO 500
                  ENDIF
                  CALL SET_SRC_MOD(frq%NFREQ,frq%NFREQ_SRF, frq%NFREQ,src%NSRC, &
                                   frq%NFREQ_SRF,src%NSRC_SRF,                  &
                                   .TRUE.,frq%CFTYPE, src%SRCTYP, SRCWRK,       &
                                   src%SOURCE,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xsrcrec25: Error calling set_src_inv 1' 
                     GOTO 500
                  ENDIF
                  CALL ROTEST(NDIM,frq%NFREQ_SRF,rcv%NREC,                 &
                              NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF,    &
                              msh%AZMOD,EST_ROT)
                  FOBJ_SRF = COBJ4(NDIM,frq%NFREQ_SRF,rcv%NREC,                &
                                   frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF, NDIM,  &
                                   .FALSE.,2,3, 0.2,0.2,0.0,                   &
                                   FREQ_LOC, WGHT_WRK,OBS_WRK,EST_WRK)
                  WRITE(*,*) 'xsrec25: Updated STF surface wave objective function', &
                             FOBJ_SRF 
               ENDIF
               IF (JOB == 2 .OR. JOB == 3) THEN
                  ALLOCATE(RECWRK(NDIM,frq%NFREQ_SRF,rcv%NREC))
                  RECWRK(:,:,:) = CZERO
                  WRITE(*,*) 'xsrcrec25: Solving new RRFs...'
                  CALL GET_RCV_MOD(NDIM,frq%NFREQ,frq%NFREQ_SRF,               &
                                   frq%NFREQ,NDIM,frq%NFREQ_SRF, rcv%NREC,     &
                                   .TRUE.,frq%CFTYPE, rcv%RECV, RECWRK, IERR)
                  CALL RECUPD(NDIM,frq%NFREQ_SRF,rcv%NREC,                               &
                              NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF,                  &
                              inv%LUNWRAP,inv%IMODREC, msh%AZMOD, EST_WRK,OBS_WRK, RECWRK)
                  IF (inv%IRESTP == 1) THEN 
                     WRITE(*,*) 'xsrcrec25: Setting receiver amplitude to unity...'
                  ELSEIF (inv%IRESTP == 2) THEN 
                     WRITE(*,*) 'xsrcrec25: Setting receiver phase to zero...'
                  ENDIF
                  CALL MODREC(NDIM,frq%NFREQ_SRF, NDIM,frq%NFREQ_SRF,rcv%NREC,  &
                              inv%IRESTP, RECWRK)
                  CALL CONVEST_RRF(NDIM,frq%NFREQ_SRF,rcv%NREC,              &    
                                   NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF, &
                                   RECWRK, EST_WRK)
                  WRITE(*,*) 'xsrcrec25: Writing new RRFs...'
                  FILENM(1:80) = ' '
                  FILENM = TRIM(ADJUSTL(FILENM))//'_srf'
                  CALL WTRECST(FILENM, NDIM,frq%NFREQ_SRF, NDIM,frq%NFREQ_SRF,rcv%NREC,  &
                               IBLOCK,0, msh%AZMOD, FREQ_LOC,RECWRK, IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xsrcrec25: Error calling wtrecst 2'
                     GOTO 500
                  ENDIF
                  CALL SET_RCV_MOD(NDIM,frq%NFREQ,frq%NFREQ_SRF,               &
                                   frq%NFREQ,NDIM,frq%NFREQ_SRF, rcv%NREC,     &
                                   .TRUE.,frq%CFTYPE, RECWRK, rcv%RECV, IERR)
                  IF (IERR /= 0) THEN 
                     WRITE(*,*) 'xsrcrec25: Error calling set_rcv_mod 1!'
                     GOTO 500
                  ENDIF
                  WRITE(*,*) 'xsrcrec25: Updating objective function...'
                  CALL ROTEST(NDIM,frq%NFREQ_SRF,rcv%NREC,                 &
                              NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF,    &
                              msh%AZMOD,EST_ROT)
                  FOBJ_SRF = COBJ4(NDIM,frq%NFREQ_SRF,rcv%NREC,                &
                                   frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF, NDIM,  &
                                   .FALSE.,2,3, 0.2,0.2,0.0,                   &
                                   FREQ_LOC, WGHT_WRK,OBS_WRK,EST_WRK)
                  WRITE(*,*) 'xsrec25: Updated RRF surface wave objective function', &
                             FOBJ_SRF
               ENDIF 
               IF (ALLOCATED( EST_WRK)) DEALLOCATE(EST_WRK)
               IF (ALLOCATED( EST_ROT)) DEALLOCATE(EST_ROT) 
               IF (ALLOCATED( OBS_WRK)) DEALLOCATE(OBS_WRK)
               IF (ALLOCATED(FREQ_LOC)) DEALLOCATE(FREQ_LOC) 
               IF (ALLOCATED(WGHT_WRK)) DEALLOCATE(WGHT_WRK)
               IF (ALLOCATED(SRCWRK))   DEALLOCATE(SRCWRK)
               IF (ALLOCATED(RECWRK))   DEALLOCATE(RECWRK)
            ENDIF
!----------------------------------------------------------------------------------------!
!                                 Repeat for body waves                                  !
!----------------------------------------------------------------------------------------!
            IF (frq%NFREQ_BDY > 0 .AND. src%NSRC_BDY > 0) THEN
               ALLOCATE( EST_WRK(NDIM,frq%NFREQ_BDY,rcv%NREC,src%NSRC_BDY))
               ALLOCATE( EST_ROT(NDIM,frq%NFREQ_BDY,rcv%NREC,src%NSRC_BDY))
               ALLOCATE( OBS_WRK(NDIM,frq%NFREQ_BDY,rcv%NREC,src%NSRC_BDY))
               ALLOCATE(FREQ_LOC(frq%NFREQ_BDY)) 
               ALLOCATE(WGHT_WRK(NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF))
               CALL GET_OE_SB(.FALSE.,  NDIM,frq%NFREQ_BDY,rcv%NREC, INV,FRQ,SRC,RCV,  &
                              WGHT_WRK,EST_WRK,OBS_WRK, IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xsrcrec25: Error calling get_oe_sb 2' 
                  GOTO 500
               ENDIF
               CALL GET_FREQ_MOD(frq%NFREQ, frq%NFREQ_BDY, .FALSE., frq%CFTYPE,frq%FREQ, &
                                 FREQ_LOC,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xsrcrec25: Error calling get_freq_mod'
                  GOTO 500
               ENDIF
               IF (JOB == 1 .OR. JOB == 3) THEN
                  WRITE(*,*) 'xsrcrec25: Solving new surface wave STFs...'
                  ALLOCATE(SRCWRK(frq%NFREQ_BDY,src%NSRC_BDY))
                  SRCWRK(:,:) = CZERO
                  CALL GET_SRC_MOD(frq%NFREQ,frq%NFREQ_BDY, frq%NFREQ,src%NSRC,  &
                                   frq%NFREQ_BDY,src%NSRC_BDY,                   &
                                   .FALSE.,frq%CFTYPE, src%SRCTYP, src%SOURCE,   &
                                   SRCWRK,IERR)

                  CALL SRCUPD(NDIM,frq%NFREQ_BDY,rcv%NREC,  &
                              NDIM,frq%NFREQ_BDY,rcv%NREC,src%NSRC_BDY,  &
                              inv%LUNWRAP,inv%IMODSRC, msh%AZMOD, EST_WRK,OBS_WRK, SRCWRK)

                  CALL MODSRC(frq%NFREQ_BDY, frq%NFREQ_BDY,src%NSRC_BDY,inv%IRESTP,  &
                              SRCWRK)

                  CALL CONVEST_STF(NDIM,frq%NFREQ_BDY,rcv%NREC,               &
                                   NDIM,frq%NFREQ_BDY,rcv%NREC,src%NSRC_BDY,  &
                                   SRCWRK, EST_WRK)

                  CALL SET_SRC_MOD(frq%NFREQ,frq%NFREQ_BDY, frq%NFREQ,src%NSRC, &
                                   frq%NFREQ_BDY,src%NSRC_BDY,                  &
                                   .FALSE.,frq%CFTYPE, src%SRCTYP, SRCWRK,      &
                                   src%SOURCE,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xsrcrec25: Error calling set_src_inv 2'
                     GOTO 500
                  ENDIF
                  FILENM(1:80) = ' '
                  FILENM = TRIM(ADJUSTL(PROJNM))//'_bdy'
                  CALL WTSTF(FILENM, frq%NFREQ_BDY, frq%NFREQ_BDY,src%NSRC_BDY, IBLOCK,0,&
                             FREQ_LOC,SRCWRK, IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xsrcrec25: Error writing STF 1'
                     GOTO 500
                  ENDIF
               ENDIF
               IF (JOB == 2 .OR. JOB == 3) THEN
                  ALLOCATE(RECWRK(NDIM,frq%NFREQ_BDY,rcv%NREC))
                  RECWRK(:,:,:) = CZERO 
                  WRITE(*,*) 'upd_srcrec_wndo: Solving new RRFs...'
                  CALL GET_RCV_MOD(NDIM,frq%NFREQ,frq%NFREQ_BDY,               &
                                   frq%NFREQ,NDIM,frq%NFREQ_BDY, rcv%NREC,     &
                                   .FALSE.,frq%CFTYPE, rcv%RECV, RECWRK, IERR)
                  CALL RECUPD(NDIM,frq%NFREQ_BDY,rcv%NREC,  &
                              NDIM,frq%NFREQ_BDY,rcv%NREC,src%NSRC_BDY, &
                              inv%LUNWRAP,inv%IMODREC, msh%AZMOD, EST_WRK,OBS_WRK, RECWRK)
                  IF (inv%IRESTP == 1) THEN
                     WRITE(*,*) 'xsrcrec25: Setting receiver amplitude to unity...'
                  ELSEIF (inv%IRESTP == 2) THEN
                     WRITE(*,*) 'xsrcrec25: Setting receiver phase to zero...'
                  ENDIF
                  CALL MODREC(NDIM,frq%NFREQ_BDY, NDIM,frq%NFREQ_BDY,rcv%NREC,  &
                              inv%IRESTP, RECWRK)
                  CALL CONVEST_RRF(NDIM,frq%NFREQ_BDY,rcv%NREC,              &
                                   NDIM,frq%NFREQ_BDY,rcv%NREC,src%NSRC_BDY, &
                                   RECWRK, EST_WRK)
                  WRITE(*,*) 'xsrcrec25: Writing new RRFs...'
                  FILENM(1:80) = ' '
                  FILENM = TRIM(ADJUSTL(FILENM))//'_bdy'
                  CALL WTRECST(FILENM, NDIM,frq%NFREQ_BDY, NDIM,frq%NFREQ_BDY,rcv%NREC,  &
                               IBLOCK,0, msh%AZMOD, FREQ_LOC,RECWRK, IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xsrcrec25: Error calling wtrecst 2'
                     GOTO 500
                  ENDIF
                  CALL SET_RCV_MOD(NDIM,frq%NFREQ,frq%NFREQ_BDY,               &
                                   frq%NFREQ,NDIM,frq%NFREQ_BDY, rcv%NREC,     &
                                   .FALSE.,frq%CFTYPE, RECWRK, rcv%RECV, IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xsrcrec25: Error calling set_rcv_mod 2!'
                     GOTO 500
                  ENDIF
               ENDIF
               IF (ALLOCATED( EST_WRK)) DEALLOCATE(EST_WRK)
               IF (ALLOCATED( EST_ROT)) DEALLOCATE(EST_ROT) 
               IF (ALLOCATED( OBS_WRK)) DEALLOCATE(OBS_WRK)
               IF (ALLOCATED(FREQ_LOC)) DEALLOCATE(FREQ_LOC) 
               IF (ALLOCATED(WGHT_WRK)) DEALLOCATE(WGHT_WRK)
               IF (ALLOCATED(SRCWRK))   DEALLOCATE(SRCWRK)
               IF (ALLOCATED(RECWRK))   DEALLOCATE(RECWRK)
            ENDIF !end check on body waves
            FOBJ = FOBJ_SRF + FOBJ_BDY
            WRITE(*,*) 'xsrcrec: Initial objective function will be:',FOBJ
         ENDIF !end check on myid
      ENDIF !end check on mynid
!
!.......... write the residual tables
!           WRITE(*,*) 'xsrcrec25: Writing residual tables...'
!           CALL RESID_TABLE(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
!                            PROJNM, IBLOCK,0, msh%AZMOD, src%CSIDE, &
!                            frq%FREQ,XREC, inv%EST,inv%OBS) 
            !WRITE(*,*) 'xsrcrec25: Writing residual RMS tables...'
            !CALL RTABLE_RMS(NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC, &
            !                PROJNM, IBLOCK,0,2, msh%AZMOD, src%CSIDE, &
            !                frq%FREQ,XREC,rcv%YREC, WGHTS,EST ,OBS)  
            !DEALLOCATE(OBS)
!        ENDIF
!     ENDIF
!----------------------------------------------------------------------------------------!
!                 Error handling and memory cleaning                                     !
!----------------------------------------------------------------------------------------!
  500 CONTINUE !break ahead for errors or out of frequencies
      IF (IERR /= 0) THEN 
         WRITE(*,9406) MYID 
 9406    FORMAT(' -------------------------------------------------------------',/, &
                ' -               An Error was Detected on Process ',I6,'     -',/, &
                ' -                  I will now abort the program             -',/, &
                ' -------------------------------------------------------------',/)
         CALL MPI_ABORT(MPI_COMM_WORLD,30,MPIERR)
      ENDIF
!
!.... free memory
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
      IF (MYID == MASTER) WRITE(*,*) 'xsrcrec25: Freeing memory...'
      IF (MYNID == MASTER) THEN
         DEALLOCATE(NSG_LSTB)
         DEALLOCATE(NSG_LST)
         DEALLOCATE(ISRCLST) 
      ENDIF
!.... background field
      IF (ALLOCATED(UE)) DEALLOCATE(UE)
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
!.... receiver
      IF (MYNID == MASTER) THEN
         IF (ASSOCIATED(rcv%YREC))  DEALLOCATE(rcv%YREC)
         IF (ASSOCIATED(rcv%MRDOF)) DEALLOCATE(rcv%MRDOF)
         IF (ASSOCIATED(rcv%RECV))  DEALLOCATE(rcv%RECV)
      ENDIF
!.... mumps
      IF (ASSOCIATED(mid%A_LOC))   DEALLOCATE(MID%A_LOC)
      IF (ASSOCIATED(mid%IRN_LOC)) DEALLOCATE(MID%IRN_LOC)
      IF (ASSOCIATED(mid%JCN_LOC)) DEALLOCATE(MID%JCN_LOC)
      IF (ALLOCATED(WAVE))         DEALLOCATE(WAVE)
      IF (MYNID == MASTER) THEN
         IF (ASSOCIATED(mid%RHS))   DEALLOCATE(mid%RHS)
         IF (ASSOCIATED(inv%EST))   DEALLOCATE(inv%EST)
      ENDIF
      IF (MYID == MASTER) THEN
         IF (ALLOCATED(XREC))     DEALLOCATE(XREC) 
         IF (ASSOCIATED(inv%OBS)) DEALLOCATE(inv%OBS) 
      ENDIF
!.... free surface
      IF (MYNID == MASTER .AND. inv%LUNWRAP) DEALLOCATE(msh%IDOF_FS) 
      MID%JOB =-2
      CALL CMUMPS(MID)
!
!.... all done 
      TESIM = MPI_WTIME()
      IF (MYID == MASTER) THEN
         WRITE(*,9405) (TESIM - TSSIM)/60.D0
 9405    FORMAT(' xsrcrec25: Inversion time:',F8.2,' minutes')
      ENDIF
      CALL MPI_FINALIZE(MPIERR)
!----------------------------------------------------------------------------------------!
!                                     Format statements                                  !
!----------------------------------------------------------------------------------------!
 9410 FORMAT(' ---------------------------------------------------------'   ,/, &
             ' -  xsrcrec25: Group:',I4,'                                -' ,/, &
             ' -             Propagating surface wave source group ',I3,' -',/, &
             ' -             With average slowness in y',E12.4,'    -'      ,/, &
             ' -             At frequency',F12.5,' Hz               -'      ,/, &
             ' ---------------------------------------------------------',/)
 9411 FORMAT(' ---------------------------------------------------------'   ,/, &
             ' -  xsrcrec25: Group:',I4,'                                -' ,/, &
             ' -             Propagating body wave source group',I3,'     -',/, &
             ' -             With average slowness in y',E12.4,'    -'      ,/, &
             ' -             At frequency',F12.5,' Hz               -'      ,/, &
             ' ---------------------------------------------------------',/) 
 9412 FORMAT(/,' xsrcrec25: Will process surface wave source ',I3,/,           &   
               '            With angle of incidence ',F8.3,' degrees',/,  &
               '            Corrected azimuth ',F8.3,' degrees',/,    &   
               '            And y slowness ',G10.3,' s/m',/)
 9413 FORMAT(/,' xsrcrec25: Will process body wave source ',I3,/,           &   
               '            With angle of incidence ',F8.3,' degrees',/,  &
               '            Corrected azimuth ',F8.3,' degrees',/,    &   
               '            And y slowness ',G10.3,' s/m',/)

      STOP
      END
