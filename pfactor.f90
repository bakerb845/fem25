!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE PFACTOR(MYID,MGCOMM,NPGROUPS,MASTER, PROJNM,IBLOCK,KITER, &
                         LPAD, INV,FRQ,SRC,RCV, P, IERR) 
!   
!     This is a direct method for solving the normal equations.  This routine:
!       (1) Sets a BLACS process grid.  
!       (2) Extracts the 2D block cyclic matrices for the Jacobian (frequency,source)
!       (3) Multiplies and adds adj(J) J
!       (4) Factors and solve H p = g
!     This is double precision since adj(J) J is likely  poorly conditioned
!     - B. Baker March 2013
!
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO) INV
      TYPE (FRQ_INFO) FRQ
      TYPE (SRC_INFO) SRC
      TYPE (RECV_INFO) RCV
      CHARACTER(80), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(IN) :: MYID, NPGROUPS, MGCOMM, MASTER, IBLOCK,KITER 
      LOGICAL*4, INTENT(IN) :: LPAD 
      REAL*4, INTENT(OUT) :: P(*) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(80) FLNAME
      CHARACTER(12) CFREQ
      CHARACTER(5) CID,CSRC 
      REAL*8, ALLOCATABLE :: JACL(:,:) 
      REAL*4, ALLOCATABLE :: JACSYS(:,:), RESID4(:), WORKR4(:)
      REAL*4 WORKQ4 
      INTEGER*4, ALLOCATABLE :: NOBSV(:) 
      INTEGER*4 DESCJ(9)  !Matrix descriptor of block cyclic Jacobian 
      INTEGER*4 DESCB(9)  !Vector descriptor for RHS 
      INTEGER*4 ICTXT !BLACS needs a context, like a handle

      real*8, allocatable :: h8(:,:), p8(:)
       real*4, allocatable :: h(:,:)
       integer*4, allocatable :: ipiv(:) 
       real*8 xpad
      INTEGER*4 STAT(MPI_STATUS_SIZE), MYSEND_ID, IUNIT, NCJ,NRJ, N,M,K,I,J,  &
                NPROW,NPCOL,NPROW0,NPCOL0, MBL,NBL, &
                MBJ,NBJ, IDEST, INPGL, IFREQ,IOBS,IREC,ISRC, IPPGRP,  &
                NPPGRP, MYROW,MYCOL, NBSAVE, MJL, NJL, KJT, KJ,  INFO, NPROCS,  &
                NOBS_IN, NDIM_IN, NREC_IN, NNPGL_IN, NWORK, NWORKL, NRHS, &
                NDIV, LDA, LDAL, LDB, NOBS, NOBSLMX, LWORK, ISAVE, MPIERR 
      integer*4 i1,i2,j1,j2,iter 
      LOGICAL*4 LAPACK 
      INTEGER*4 NUMROC
      PARAMETER(NRHS = 1) !Only one rhs here
!
!----------------------------------------------------------------------------------------!
!
!.... get size 
      CALL MPI_COMM_SIZE(MGCOMM, NPROCS, MPIERR) !get communicator size 
!
!.... may need to switch to pure LAPACK
      LAPACK = .FALSE.
      IF (NPROCS == 1) THEN 
         WRITE(*,*) 'pfactor: Only one process available!'
         WRITE(*,*) 'pfactor: Switching to LAPACK...'
         LAPACK = .TRUE.
      ENDIF
      write(*,*) 'override'
      lapack = .true. 
!
!.... i need to classify the problem
      NPPGRP = NPROCS/NPGROUPS !number of process per group, evenly divisible 
      IF (MYID == MASTER) THEN  
         CALL ICNOBS_JMAT(frq%NFREQ,src%NSRC, MYID, frq%FREQ, NOBSLMX,NOBS,IERR)
         WRITE(*,*) 'pfactor: Number of observations:',NOBS 
         WRITE(*,*) 'pfactor: Number of inversion variables',inv%NA35
      ENDIF 
      IF (MYID == MASTER .AND. LAPACK) THEN
         WRITE(*,*) 'pfactor: Assembling global Jacobian...',inv%NA35
         N = inv%NA35 !Jacobian columns
         M = NOBS     !Jacobian rows 
         IF (M /= inv%NOBS) THEN
            WRITE(*,*) 'pfactor: Error observation count mismatch!',m,inv%nobs
            IERR = 1
            RETURN
         ENDIF
         ALLOCATE(JACSYS(M,N),STAT=IERR) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'pfactor: Error setting space for Jacobian for LAPACK!'
            RETURN
         ENDIF 
         ALLOCATE(NOBSV(frq%NFREQ*src%NSRC))
         CALL ASMBLE_JACSYS(M,NOBSLMX,inv%NA35,frq%NFREQ, NDIM,src%NSRC,rcv%NREC,  &
                            NPPGRP, frq%FREQ, NOBSV,JACSYS,IERR) 
!        DEALLOCATE(NOBSV) 
         LDB = MAX(M,N)
         ALLOCATE(RESID4(LDB)) 
         RESID4(1:M) = inv%RHS(1:M) 
         LWORK =-1
!        CALL SGELS('N', M, N, 1, JACSYS, NOBS, RESID4, LDB, WORKQ4,LWORK,  &
!                   INFO)
!        IF (INFO /= 0) THEN
!           WRITE(*,*) 'pfactor: Error calling querying SGELS!',INFO
!           IERR = 1
!           RETURN
!        ENDIF
!        LWORK = INT(WORKQ4)
         WRITE(*,*) 'pfactor: Workspace size:',LWORK
!        ALLOCATE(WORKR4(LWORK))
         WRITE(*,*) 'pfactor: Calling SGELS...'
!        CALL SGELS('N', M, N, 1, JACSYS, NOBS, RESID4, LDB, WORKR4,LWORK, &
!                   INFO)  
         IF (IERR /= 0) THEN
            WRITE(*,*) 'pfactor: Error calling SGELS!',INFO
            IERR = 1
            RETURN
         ENDIF
         P(1:N) = RESID4(1:N) !copy inversion back 
!        DEALLOCATE(WORKR4) 

         p(1:n) = inv%grad(1:n) 
         print *, n,inv%na35
         allocate(ipiv(n))
         allocate(h(n,n))
         h(:,:) = 0.0
         print *, 'mamtul'
         !h = matmul(transpose(jacsys),jacsys)
         call sgemm('T','N', n,n,m, 1.0,jacsys,m, jacsys,m, 0.0,h,n)
         !pad the diagonal
         write(*,*) 'sorting...'
         open(unit=44,file='diag_hess.dat')
         allocate(inv%gradpc(inv%na35))
         xpad = 1.d0
         do i=1,n
            inv%gradpc(i) = h(i,i)
            write(44,*) i,inv%gradpc(i)
            xpad = xpad*inv%gradpc(i)**(1.d0/dfloat(n)) !try geometric mean
         enddo
         close(44)
         !call rshell1(inv%na35,inv%gradpc)
         write(*,*) 'padding diagonal',xpad, minval(inv%gradpc),maxval(inv%gradpc)
         xpad = sum(inv%gradpc)/dfloat(n)
         write(*,*) 'trying trace',xpad
         j = int(inv%na35*.1)
         do i=1,n
           !if (h(i,i) < inv%gradpc(j)) h(i,i) = inv%gradpc(j) 
           h(i,i) = h(i,i) !+ xpad
         enddo
          
         print *, 'dump mat'
         call plot_aphess_vtk(projnm, n,n, iblock,kiter, h)
!        open(unit=44,file='hess_dump.dat')
!        do i=1,n
!           do j=1,n
!              write(44,*) i,n+1-j,h(i,j)
!           enddo
!           write(44,*)
!        enddo
!        close(44) 
         print *, 'call sgetrf'
         i1 = 1
         iter = 0 
!        do ifreq=1,frq%nfreq
!           do isrc=1,src%nsrc
!              iter = iter + 1
!              i2 = i1 + nobsv(iter) - 1
!              j1 = i2 + 1
!              j2 = j2 + nobsv(iter) - 1 
!              h = h + matmul(transpose(jacsys(i1:i2,:)),jacsys(i1:i2,:)) 
!              h = h + matmul(transpose(jacsys(j1:j2,:)),jacsys(j1:j2,:))
!              i1 = j2 + 1
!           enddo
!        enddo 
         allocate(h8(n,n))
         allocate(p8(n))
         h8(:,:) = dble(h(:,:))
         call dgetrf(n,n,h8,n,ipiv,info)  
         print *, 'back',info
         if (info == 0) then
            print *, 'calling solve'
            p8(1:n) = dble(inv%grad(1:n))
            call dgetrs('N',n,1,h8,n,ipiv, p8,n,info)
            p(1:n) =-sngl(p8(1:n))
         endif 

         DEALLOCATE(RESID4)
         DEALLOCATE(JACSYS)
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
      IF (LAPACK) RETURN 
      CALL MPI_BCAST(    NOBS,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR) 
      CALL MPI_BCAST(inv%NA35,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR) 
      NCJ = inv%NA35 !Jacobian columns
      NRJ = NOBS     !number of observations comprises rows of Jacobian
      M = NRJ        !number of rows of J, nobs observations 
      N = NCJ        !number of columns of J, m model parameters
!
!.... broadcast the RHS to all processes
!     ALLOCATE(BW(inv%NA35)) 
!     IF (MYID == MASTER) BW(1:inv%NA35) = DBLE(inv%GRAD(1:inv%NA35))
!     CALL MPI_BCAST(BW,inv%NA35,MPI_DOUBLE_PRECISION, MASTER,MGCOMM,MPIERR)
!
!.... have head process figure out process row/column size 
      IERR = 0
      IF (MYID == MASTER) THEN
!
!....... i'm capping this at 36 processes which ~size of a receiver array
         WRITE(*,*) 'pfactor: Estimating process grid sizes...'
         NPROW = 1 
         NPCOL = 1
         IF (NPROCS > 36) THEN
            NPROW = 6
            NPCOL = 6
            GOTO 500
         ENDIF
!
!....... try to get squarish block sizes
         NPCOL0 = 1 
         NPROW = 1 
         DO 1 I=1,8
            NPCOL = NPROW 
            IF (NPROW == NRJ .OR. NPCOL == NCJ) GOTO 500 !this would be bad
            DO 2 J=1,2
               IF (NPCOL*NPROW < NPROCS) THEN !save
                  NPCOL0 = NPCOL
                  NPROW0 = NPROW
               ELSEIF (NPCOL*NPROW > NPROCS) THEN 
                  NPROW = NPROW0 
                  NPCOL = NPCOL0
                  GOTO 500 
               ELSE !exactly equal, awesome
                  GOTO 500 
               ENDIF
               NPCOL = NPCOL + 1 
    2       CONTINUE 
            NPROW = NPROW + 1 
    1    CONTINUE 
  500    CONTINUE !break ahead
         IF (NPROW*NPCOL > NPROCS) THEN
            WRITE(*,*) 'Serious error nprow*npcol > nprocs!'
            IERR = 1
            RETURN
         ENDIF
!
!....... tall matrices would want Pr > Pc 
         IF (NOBS > NCJ) THEN
            ISAVE = NPROW
            NPROW = NPCOL  
            NPCOL = ISAVE 
         ENDIF 
!
!....... try to get a block size around 50 - 100
         IF (NRJ/NPROW > 50) THEN
            MBJ = MIN(50,NRJ/NPROW)
         ELSE
            MBJ = MIN(MAX(50,NRJ/NPROW),NRJ)
         ENDIF 
         IF (NCJ/NPCOL > 50) THEN
            NBJ = MIN(50,NCJ/NPCOL)
         ELSE
            NBJ = MIN(MAX(50,NCJ/NPCOL),NCJ)
         ENDIF
      ENDIF
      CALL MPI_BCAST(NPROW,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(NPCOL,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR)  
      CALL MPI_BCAST(NBJ  ,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(MBJ  ,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR)
      IF (MYID == MASTER) WRITE(*,9005) NPROW,NPCOL,MBJ,NBJ
 9005 FORMAT(' pfactor: Number of process rows:',I4,   /, &
             '          Number of process columns:',I4,/, &
             '          Jacobian row block size:',I6,/,   &
             '          Jacobian column block size:',I6,/)
!
!.... initialize the grid
      CALL SL_INIT(ICTXT, NPROW, NPCOL) !make a process grid in row major ordering
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL) !set BLACS grid
      IF (MYROW ==-1 .OR. MYCOL ==-1) GOTO 1000 !didn't make the team, maybe next time
!.... calculate the NUMber of Rows Or Columns (NUMROC) 
      MJL  = NUMROC(M,   MBJ,  MYROW, 0, NPROW) !number of local rows of J 
      NJL  = NUMROC(N,   NBJ,  MYCOL, 0, NPCOL) !number of local columms of J
      MBL  = NUMROC(N,   NBJ,  MYROW, 0, NPROW) !number of local rows of B
      NBL  = NUMROC(NRHS,  1,  MYCOL, 0, NPCOL) !number of local columns of RHS matrix B
      NBSAVE = NBL
!     print *, 'rows,cols of jt',myid,mjt,kjt !mjt,kjt,kj,nj
      !print *, 'kjt',mjt,kjt, mbjt,nbjt
!.... create the descriptors for each matrix, adj(J), J, and adj(J)J
      INFO = 0
      CALL DESCINIT(DESCJ ,M,N, MJL ,NJL,  0,0, ICTXT, MAX0(1,MJL),INFO)
      IF (INFO < 0) WRITE(*,*) 'pfactor: DESCINIT 2 illegal argument:',INFO
!     CALL DESCINIT(DESCB ,N,1, NJ , 1,  0,0, ICTXT, MAX0(1,MBL),INFO)   
      !print *, descj
      print *, myrow,mycol,mjl,njl
!
!.... set space
      ALLOCATE(JACL(DESCJ(9),MAX0(NJL,1)))    !distributed Jacobian 
      JACL(1:DESCJ(9),1:MAX0(NJL,1)) = 0.D0 
!
!.... have (0,0) grid read the jacobian 
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         WRITE(*,*) 'pfactor: Reading Jacobian...'
      ENDIF
!     ALLOCATE(AJACR(DESCJT(9),MAX(KJT,1))) !distributed real transpose jacobian
!     ALLOCATE( JACR(DESCJ(9) ,MAX(NJ ,1))) !distributed real jacobian
!     ALLOCATE(AJACI(DESCJT(9),MAX(KJT,1))) !distributed imaginary transpose jacobian
!     ALLOCATE( JACI(DESCJ(9) ,MAX(NJ ,1))) !distributed imaginary jacobian
!     ALLOCATE(HMAT(DESCH(9) ,MAX(NJ ,1))) !complex version of Hmat
!     ALLOCATE(   B(DESCB(9) ,MAX(NB ,1))) !distributed right hand side matrix
!     HMAT(:,:) = 0.D0 !this needs to be zero'd out b/c it is stacked
!     AJACR(:,:) = 0.D0 
!     JACR(:,:) = 0.D0 
!     AJACI(:,:) = 0.D0
!     JACI(:,:) = 0.D0
!
!.... set reading space
!     ALLOCATE(LOBS(NDIM,rcv%NREC))
!     ALLOCATE(RJAC(NRJ,NCJ))
!     ALLOCATE(CJAC(NRJ,NCJ))
!     RJAC(:,:) = 0.0
!     CJAC(:,:) = 0.0
!     ALLOCATE(ZJAC(NRJ,NCJ))
!     ALLOCATE(MYDEST(NCJ))  
!     MYDEST(:) = 0
!     XPAD = 0.D0 !regularization estimate
!     NDIV = 0
!
!.... now calculate the space required to hold each matrix
!     DO 101 IFREQ=1,frq%NFREQ
!        CFREQ(1:12) = ' '
!        WRITE(CFREQ,'(F12.5)') frq%FREQ(IFREQ)
!        CFREQ = ADJUSTL(CFREQ)
!        DO 102 ISRC=1,src%NSRC
!           CSRC(1:5) = ' '
!           WRITE(CSRC,'(I5)') ISRC 
!           CSRC = ADJUSTL(CSRC)
!
!.......... loop on process groups to fill entire Jacobian 
!           ZJAC(:,:) = ZZERO
!           DO 103 IPPGRP=1,NPPGRP
!              IUNIT = 10 + MYID
!              WRITE(CID,'(I5)') IPPGRP - 1
!              CID = ADJUSTL(CID) 
!              FLNAME(1:80) = ' '
!              FLNAME = './scratch/jac_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'//TRIM(CSRC)//'.dat'
!              FLNAME = ADJUSTL(FLNAME)
!              OPEN(FILE=TRIM(FLNAME),UNIT=IUNIT,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=IERR)
!              CALL LOAD_JAC(NDIM,NRJ, NDIM,rcv%NREC, IUNIT,.TRUE.,.TRUE., MYID,  &
!                            NOBS_IN,NNPGL_IN,LOBS,RJAC,CJAC,MYDEST, IERR) 
!              IF (IERR /= 0) THEN
!                 WRITE(*,*) 'pfactor: Error loading Jacobian!'
!                 RETURN
!              ENDIF 
!
!............. copy the jacobian
!              IOBS = 0
!              K = 0 
!              DO 104 IREC=1,rcv%NREC
!                 DO 105 I=1,NDIM
!                    K = K + 1 
!                    DO 106 INPGL=1,NNPGL_IN
!                       IDEST = MYDEST(INPGL)
!                       IF (IPPGRP == 1) THEN
!                          CALL SSCAL(inv%NA35,0.0,GJACR(K,1:inv%NA35),1) 
!                          CALL SSCAL(inv%NA35,0.0,GJACI(K,1:inv%NA35),1) 
!                       ENDIF 
!                       IF (LOBS(I,IREC)) THEN 
!                          IOBS = IOBS + 1
!                          NDIV = NDIV + 1
!                          GJACR(K,IDEST) = RJAC(IOBS,INPGL)
!                          GJACI(K,IDEST) = CJAC(IOBS,INPGL)
!                          XPAD = XPAD + DJACR(K,IDEST)**2 + DJACI(K,IDEST)**2 
!                       ENDIF
! 106                CONTINUE !loop on 
! 105             CONTINUE !loop on components
! 104          CONTINUE !loop on receivers
! 103       CONTINUE !loop on processes
!
!.......... fill Re(J), Re(J^T), then Im(J), Im(J^T) 
!           IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'pfactor: Filling real Jacobians...'
!           LDATL = DESCJT(9) !leading dimension 
!           LDAL  = DESCJ(9)  !leading dimension
!           LDA  = NRJ       !number of rows in Jacobian 
!           CALL FILL_DMATS_2DC(LDATL,LDAL,LDA, MYROW,MYCOL, NPROW,NPCOL,  &
!                               DESCJT, DJACR, AJACR,JACR,IERR)
!           IF (IERR /= 0) THEN
!              WRITE(*,*) 'pfactor: An error occurred on process:',MYID
!              RETURN
!           ENDIF
!           IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) &
!           'pfactor: Filling imaginary Jacobians...'
!           CALL FILL_DMATS_2DC(LDALT,LDAL,LDA, MYROW,MYCOL, NPROW,NPCOL,  &
!                               DESCJT, DJACI, AJACI,JACI,IERR)
!           IF (IERR /= 0) THEN 
!              WRITE(*,*) 'pfactor: An error occurred on process:',MYID
!              RETURN
!           ENDIF

!           IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'pfactor: Filling search direction...'
!           NWORKL = DESCB(9) !leading dimension
!
!.......... multiply adj(J) J and update H = adj(J)*J + H
!           IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'pfactor: Multiplying adj(J)J...'
!           CALL PDGEMM('N', 'N', M, N, K, 1.D0,       &   
!                       AJACR, 1, 1, DESCJT, JACR, 1, 1, DESCJ,   &   
!                       1.D0, HMAT, 1, 1, DESCH)
!           CALL PDGEMM('N', 'N', M, N, K, 1.D0,       &
!                       AJACI, 1, 1, DESCJT, JACI, 1, 1, DESCJ,   &
!                       1.D0, HMAT, 1, 1, DESCH)
!           IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) ''
!           print *, myid,minval(dabs(hmat)),maxval(dabs(hmat))
! 102    CONTINUE
! 101 CONTINUE 
!     XPAD = DSQRT(XPAD)/DFLOAT(NDIV) 
 
!.... free some space
!     DEALLOCATE(RJAC) 
!     DEALLOCATE(CJAC) 
!     DEALLOCATE(ZJAC) 
!     DEALLOCATE(AJAC) 
!     DEALLOCATE( JAC) 
!     DEALLOCATE(LOBS)
!     DEALLOCATE(MYDEST)
!
!.... fill the RHS 
!     IF (NBSAVE > 0) THEN 
!        CALL FILL_VEC_2DC(NWORKL,N, MYROW, NPROW,  DESCB, BW, &
!                          B,IERR) 
!        IF (IERR /= 0) THEN
!           WRITE(*,*) 'pfactor: Error setting RHS on process:',MYID
!           RETURN
!        ENDIF
!     ELSE
!        B(:,:) = 0.D0
!     ENDIF
!     ALLOCATE(WORK(DESCJT(3))) 
!     CALL PSLAPRNT(M, K, AJAC, 1, 1, DESCJT, 0, 0,  &
!                    'A', 6, WORK)
!     IF (ALLOCATED(WORK)) DEALLOCATE(WORK) 
!     ALLOCATE(WORK(DESCB(3))) 
!     CALL PDLAPRNT(N, 1, B, 1, 1, DESCB, 0, 0,  &
!                    'X', 6, WORK)
!     IF (ALLOCATED(WORK)) DEALLOCATE(WORK) 

!     ALLOCATE(WORK(DESCH(3)))
!     CALL PDLAPRNT(M, N, HMAT, 1, 1, DESCH, 0, 0,  &
!                   'A', 6, WORK)
!     IF (ALLOCATED(WORK)) DEALLOCATE(WORK) 
!     ALLOCATE(DIAG(NCJ)) 
!     DIAG(:) = 0.D0
!     LDAL = DESCH(9)  
!     CALL GET_DIAG_2DC(LDAL, MYROW,MYCOL, NPROW,NPCOL, DESCH, HMAT,DIAG,IERR)
!     XPAD = HUGE(1.D0) 
!     DO I=1,NCJ
!        IF (DIAG(I) > 0.D0) XPAD = DMIN1(DIAG(I),XPAD)
!        if (diag(i) < 0.d0) write(*,*) 'error'
!     ENDDO
!     print *, 'xpad1',xpad 
!     XPAD = (MAXVAL(DIAG)*0.25D0 + XPAD*0.75D0)/2.D0
!     DEALLOCATE(DIAG) 
!     print *, 'xpad2',xpad
!
!.... pad matrix for this example
!     IF (LPAD) THEN
!        LDAL = DESCH(9)
!        CALL PAD_MAT_2DC(LDAL, MYROW,MYCOL, NPROW,NPCOL, DESCH, XPAD,HMAT,IERR)
!        IF (IERR /= 0) THEN
!           WRITE(*,*) 'pfactor: Error padding matrix!',MYID,MYROW,MYCOL
!           RETURN
!        ENDIF
!     ENDIF 
!
!.... factor
!     NWORK = NUMROC(M,MBJT, MYROW, 0, NPROW) + DESCB(5) 
!     ALLOCATE(IPIV(NWORK)) 
!     IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'pfactor: Solving Ax=b...'
!     NRHS = 1
!     CALL PDGESV(M, NRHS, HMAT, 1, 1, DESCH, IPIV, B, 1, 1, DESCB,  &
!                 INFO)   
!     IF (INFO /= 0) THEN
!        IF (INFO < 0) THEN
!           WRITE(*,*) 'pfactor: Invalid input on process:',MYID,INFO
!           IERR = 1
!        ELSE
!           WRITE(*,*) 'pfactor: The matrix is singular at index:',MYID,INFO 
!           IERR =-1
!        ENDIF
!     ENDIF
!
!.... retrieve solution
!     IF (MYROW == 0 .AND. MYCOL == 0) THEN
!        WRITE(*,*) 'pfactor: Collecting solution onto (0,0) process...'
!        ALLOCATE(BGLOB(N))
!     ELSE
!        ALLOCATE(BGLOB(1))
!     ENDIF
!     LDB = DESCB(9)
!     CALL FCH_SOL_2DC(ICTXT,LDB,DESCB, NPROW,NPCOL, MYROW,MYCOL, B,  &
!                      BGLOB,IERR)
!     IF (MYROW == 0 .AND. MYCOL == 0) then
!        do i=1,n
!           write(*,*) i,bglob(i)
!        enddo
!     endif
!     DEALLOCATE(BGLOB)

!     ALLOCATE(WORK(DESCB(3))) 
!     CALL PDLAPRNT(N, 1, B, 1, 1, DESCB, 0, 0,  &
!                    'B', 6, WORK)
!     IF (ALLOCATED(WORK)) DEALLOCATE(WORK) 
      CALL BLACS_GRIDEXIT(ICTXT) !free the BLACS context 
 1000 CONTINUE !nothing for this process to do 
      CALL MPI_BARRIER(MGCOMM,MPIERR)
!
!.... make sure master process gets solution
!     MYSEND_ID =-1 
!     IF (MYROW == 0 .AND. MYCOL == 0) MYSEND_ID = MYID
!     CALL MPI_BCAST(MYSEND_ID,1,MPI_INTEGER, MYID,MGCOMM,MPIERR)  
!     IF (MYSEND_ID /= MASTER) THEN
!        IF (MYID == MYSEND_ID) THEN
!           CALL MPI_SEND(BGLOB,inv%NA35,MPI_DOUBLE_PRECISION, MASTER,MYSEND_ID, &
!                         MGCOMM,MPIERR) 
!        ELSE
!           IF (MYID == MASTER) THEN
!              WRITE(*,*) 'pfactor: Receiving search direction from process:',MYSEND_ID
!              ALLOCATE(BGLOB(inv%NA35)) 
!              CALL MPI_RECV(BGLOB,inv%NA35,MPI_DOUBLE_PRECISION, MYSEND_ID, &
!                            MPI_ANY_TAG,MGCOMM, STAT,MPIERR) 
!           ENDIF
!        ENDIF 
!     ENDIF
!     IF (MYID == MASTER) THEN
!        WRITE(*,*) 'pfactor: Filling search direction...' 
!        P(1:inv%NA35) =-SNGL(BGLOB(1:inv%NA35)) !negative for minimization 
!     ENDIF
      CALL BLACS_EXIT(1) !i will continue with after BLACS closes (otherwise 0)
      CALL MPI_BARRIER(MGCOMM,MPIERR) !block until everyone else is here
      RETURN 
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE DIST_MAT2D( ) 
!
!     Distributes the global Jacobian to local block cyclic structure

      RETURN 
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE FILL_DMATS_2DC(LDATL,LDAL,LDA, MYROW,MYCOL, NPROW,NPCOL,  &
                                DESCA, A, ALOCT,ALOC,IERR) 
!
!     Takes a full dense matrix then extracts the corresponding Scalapack local matrix 
!     for A and transpose(A)
!
!     INPUT      MEANING
!     -----      -------
!     A          matrix to extract
!     DESCA      BLACS descriptor of Aloc
!     LDA        leading dimension of A
!     LDAL       leading dimension of 2d block cyclic matrix ALOC
!     MYCOL      process column in BLACS process grid 
!     MYROW      process row in BLACS process grid
!     NPCOL      number of process rows in BLACS grid
!     NPROW      number of process columsn in BLACS grid
!
!     OUTPUT     MEANING
!     ------     ------- 
!     ALOC       processes' 2d block cyclic part of A 
!     ALOCT      processes' 2d block cyclic part of A^T
!     IERR       error flag
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: A(LDA,*)
      INTEGER*4, INTENT(IN) :: DESCA(9), LDATL, LDAL, LDA, MYROW, MYCOL, NPROW, NPCOL   
      REAL*8, INTENT(OUT) :: ALOCT(LDATL,*), ALOC(LDAL,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 M, N, IC, IOFFC, I, IR, IOFFR, J, IK, JK, IAROW, IACOL,  &
                IRNUM, ICNUM, MB, NB, MEND, NEND,      &
                MRROW, MRCOL, MOFF, NOFF, IROW_GLOB, JCOL_GLOB 
      INTEGER*4 NUMROC, ICEIL
!
!----------------------------------------------------------------------------------------!
!
!.... initialize
      IERR = 0
      M  = DESCA(3)    !number of rows of matrix
      N  = DESCA(4)    !number of columns of matrix
      MB = DESCA(5)    !row block size
      NB = DESCA(6)    !column block size
      IRNUM = NUMROC(M,MB,MYROW,0,NPROW)
      ICNUM = NUMROC(N,NB,MYCOL,0,NPCOL) 
      IAROW = DESCA(7) !process row over which first row of A distributed
      IACOL = DESCA(8) !process col over which first col of A distributed
      IF (MB <= 0) THEN
         WRITE(*,*) 'fill_mat_2dc: Serious error 1, mb <= 0',MB 
         IERR = 1 
         RETURN
      ENDIF
      IF (NB <= 0) THEN
         WRITE(*,*) 'fill_mat_2dc: Serious error 2, nb <= 0',NB
         IERR = 2
         RETURN
      ENDIF
      MOFF = NPROW/MB !#process rows/block size
      NOFF = NPCOL/NB !#process columns/block size
      MEND = ICEIL(IRNUM, MB) + MOFF
      NEND = ICEIL(ICNUM, NB) + NOFF
      MRROW = MOD( NPROW+MYROW-IAROW, NPROW )
      MRCOL = MOD( NPCOL+MYCOL-IACOL, NPCOL )
!
!.... fill local matrices from global matrix
      JK = 1
      DO 290 IC = NOFF+1, NEND !loop on column blocks
         IOFFC = ((IC-1)*NPCOL+MRCOL)*NB
         DO 280 I = 1, NB !loop on columns in block 
            IF(JK .GT. ICNUM) GO TO 300 !done with rows
            IK = 1 
            DO 260 IR = MOFF+1, MEND !loop on row blocks
               IOFFR = ((IR-1)*NPROW+MRROW)*MB
               DO 250 J = 1, MB !Loop on rows in block
                  IF( IK .GT. IRNUM ) GO TO 270 !out of local rows
                  IROW_GLOB = IOFFR + J !row of matrix
                  JCOL_GLOB = IOFFC + I 
                  ALOCT(IK,JK) = DBLE(A(JCOL_GLOB,IROW_GLOB))
                  ALOC(IK,JK)  = DBLE(A(IROW_GLOB,JCOL_GLOB))
                  IK = IK + 1 
  250          CONTINUE !loop on rows in block
  260       CONTINUE !loop on row blocks
  270       CONTINUE !out of local rows
            JK = JK + 1 
  280    CONTINUE !loop on rows
  290 CONTINUE !loop on column blocks 
  300 CONTINUE !done
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE FCH_SOL_2DC(ICTXT,LDB,DESCB, NPROW,NPCOL, MYROW,MYCOL, B,  &
                             BGLOB,IERR) 
!
!     Fetch the solution from the 2D blocky cyclic structure.  For now only works on 
!     vectors
!
!     INPUT      MEANING
!     -----      ------- 
!     B          2d block cylclic vector
!     DESCB      B BLACS descriptor
!     ICTXT      BLACS grid context
!     LDB        leading dimension for B vector
!     MYCOL      processes' column
!     MYROW      processes' row 
!     NPCOL      number of process columns
!     NPROW      number of process rows
! 
!     OUTPUT     MEANING
!     ------     -------
!     BGLOB      global vector on (0,0) process
!     IERR       error flag
!
!.... variable declarations
      IMPLICIT NONE 
      REAL*8, INTENT(IN) :: B(LDB,1)
      INTEGER*4, INTENT(IN) :: DESCB(9), ICTXT,LDB, NPROW,NPCOL, MYROW,MYCOL
      REAL*8, INTENT(OUT) :: BGLOB(*)  
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables 
      REAL*8, ALLOCATABLE :: BUFF(:) 
      INTEGER*4 M,MB, IAROW,MOFF, IPROW,IPCOL, IRNUM, MEND,MRROW, IK, IR, J, &
                IOFFR, IROW_GLOB, MREC, IWORK 
      INTEGER*4 NUMROC, ICEIL 
!
!----------------------------------------------------------------------------------------!
!
!.... initialize
      IERR = 0 
      M  = DESCB(3)    !M_
      MB = DESCB(5)    !row block size
      IF (MB <= 0) THEN
         WRITE(*,*) 'fill_vec_2dc: Serious error 1, mb <= 0',MB
         IERR = 1
         RETURN
      ENDIF
      IAROW = DESCB(7) !process row over which first row of A distributed
      MOFF = NPROW/MB  !#process rows/block size
      ALLOCATE(BUFF(M))
      BUFF(1:M) = 0.D0
!
!.... master process does fetching
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
!
!....... loop on process grid
         DO 100 IPROW=0,NPROW-1
            DO 101 IPCOL=0,NPCOL-1
               IRNUM = NUMROC(M,MB,IPROW,0,NPROW)
               MEND = ICEIL(IRNUM, MB) + MOFF
               MRROW = MOD( NPROW+IPROW-IAROW, NPROW )
               IWORK = NUMROC(1,1,IPCOL,0,NPCOL)
               IF (IPCOL == 0 .AND. IPROW == 0) THEN !just getting it from myself
                  MREC = DESCB(9) !leading dimensions 
                  BUFF(1:MREC) = B(1:MREC,1)
               ELSE !receiving from other processes
                  CALL IGERV2D(ICTXT,1,1,IWORK,1,IPROW,IPCOL)
                  IF (IWORK > 0) THEN
                     CALL IGERV2D(ICTXT,   1,1,MREC,   1,IPROW,IPCOL)
                     CALL DGERV2D(ICTXT,MREC,1,BUFF,MREC,IPROW,IPCOL)
                  ENDIF
               ENDIF
!
!............. check if there is something to fill, then fill global vector 
               IF (IWORK > 0) THEN
                  IK = 1
                  DO 260 IR = MOFF+1, MEND !loop on row blocks
                     IOFFR = ((IR-1)*NPROW+MRROW)*MB
                     DO 250 J = 1, MB !Loop on rows in block
                        IF ( IK .GT. IRNUM ) GO TO 270 !out of local rows
                        IROW_GLOB = IOFFR + J !row of matrix
                        BGLOB(IROW_GLOB) = BUFF(IK)
                        IK = IK + 1
  250                   CONTINUE !loop on rows in block
  260                CONTINUE !loop on row blocks
  270             CONTINUE !out of local rows
               ENDIF !end check if there is anything to fill
  101       CONTINUE !loop on columns in BLACS grid
  100    CONTINUE !loop on rows in BLACS grid
      ELSE !other process just send info
         IWORK = NUMROC(1,1,MYCOL,0,NPCOL)
         CALL IGESD2D(ICTXT,1,1,IWORK,1,0,0)
         IF (IWORK > 0) THEN
            MREC = DESCB(9) !leading dimension  
            CALL IGESD2D(ICTXT,   1,1,       MREC,    1,0,0)
            CALL DGESD2D(ICTXT,MREC,1,B(1:MREC,1), MREC,0,0)
         ENDIF
      ENDIF
      DEALLOCATE(BUFF)
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE PAD_MAT_2DC(LDAL, MYROW,MYCOL, NPROW,NPCOL, DESCA, XLAM,ALOC,IERR) 
!
!     Pads the main diagonal of the matrix A in 2d block cycle storage
!
!     INPUT      MEANING
!     -----      -------
!     ALOC       matrix to pad
!     DESCA      BLACS descriptor of Aloc
!     LDAL       leading dimension of 2d block cyclic matrix ALOC
!     MYCOL      process column in BLACS process grid 
!     MYROW      process row in BLACS process grid
!     NPCOL      number of process rows in BLACS grid
!     NPROW      number of process columsn in BLACS grid
!     XLAM       padding factor
!
!     OUTPUT     MEANING
!     ------     ------- 
!     ALOC       ALOC with main diagonal padded by xlam
!     IERR       error flag
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: ALOC(LDAL,*)
      REAL*8, INTENT(IN) :: XLAM
      INTEGER*4, INTENT(IN) :: DESCA(9), LDAL, MYROW, MYCOL, NPROW, NPCOL   
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 M, N, IC, IOFFC, I, IR, IOFFR, J, IK, JK, IAROW, IACOL,  &
                IRNUM, ICNUM, MB, NB, MEND, NEND,      &   
                MRROW, MRCOL, MOFF, NOFF, IROW_GLOB, JCOL_GLOB 
      INTEGER*4 NUMROC, ICEIL
!
!----------------------------------------------------------------------------------------!
!
!.... initialize
      IERR = 0
      M  = DESCA(3)    !number of rows of matrix
      N  = DESCA(4)    !number of columns of matrix
      MB = DESCA(5)    !row block size
      NB = DESCA(6)    !column block size
      IRNUM = NUMROC(M,MB,MYROW,0,NPROW)
      ICNUM = NUMROC(N,NB,MYCOL,0,NPCOL)
      IAROW = DESCA(7) !process row over which first row of A distributed
      IACOL = DESCA(8) !process col over which first col of A distributed
      IF (MB <= 0) THEN
         WRITE(*,*) 'pad_mat_2dc: Serious error 1, mb <= 0',MB
         IERR = 1
         RETURN
      ENDIF
      IF (NB <= 0) THEN
         WRITE(*,*) 'pad_mat_2dc: Serious error 2, nb <= 0',NB
         IERR = 2
         RETURN
      ENDIF
      MOFF = NPROW/MB !#process rows/block size
      NOFF = NPCOL/NB !#process columns/block size
      MEND = ICEIL(IRNUM, MB) + MOFF
      NEND = ICEIL(ICNUM, NB) + NOFF
      MRROW = MOD( NPROW+MYROW-IAROW, NPROW )
      MRCOL = MOD( NPCOL+MYCOL-IACOL, NPCOL )
!
!.... fill local matrices from global matrix
      JK = 1
      DO 290 IC = NOFF+1, NEND !loop on column blocks
         IOFFC = ((IC-1)*NPCOL+MRCOL)*NB
         DO 280 I = 1, NB !loop on columns in block 
            IF(JK .GT. ICNUM) GO TO 300 !done with rows
            IK = 1
            DO 260 IR = MOFF+1, MEND !loop on row blocks
               IOFFR = ((IR-1)*NPROW+MRROW)*MB
               DO 250 J = 1, MB !Loop on rows in block
                  IF( IK .GT. IRNUM ) GO TO 270 !out of local rows
                  IROW_GLOB = IOFFR + J !row of matrix
                  JCOL_GLOB = IOFFC + I
                  IF (JCOL_GLOB == IROW_GLOB) ALOC(IK,JK) = ALOC(IK,JK) + XLAM 
                  IK = IK + 1
  250          CONTINUE !loop on rows in block
  260       CONTINUE !loop on row blocks
  270       CONTINUE !out of local rows
            JK = JK + 1
  280    CONTINUE !loop on rows
  290 CONTINUE !loop on column blocks 
  300 CONTINUE !done
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE GET_DIAG_2DC(LDAL, MYROW,MYCOL, NPROW,NPCOL, DESCA, ALOC,DIAG,IERR)
!
!     Pads the main diagonal of the matrix A in 2d block cycle storage
!
!     INPUT      MEANING
!     -----      -------
!     ALOC       matrix to get diagonal of 
!     DESCA      BLACS descriptor of Aloc
!     LDAL       leading dimension of 2d block cyclic matrix ALOC
!     MYCOL      process column in BLACS process grid 
!     MYROW      process row in BLACS process grid
!     NPCOL      number of process rows in BLACS grid
!     NPROW      number of process columsn in BLACS grid
!
!     OUTPUT     MEANING
!     ------     ------- 
!     DIAG       diagonal of matrix A 
!     IERR       error flag
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: ALOC(LDAL,*)
      INTEGER*4, INTENT(IN) :: DESCA(9), LDAL, MYROW, MYCOL, NPROW, NPCOL
      REAL*8, INTENT(OUT) :: DIAG(*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 M, N, IC, IOFFC, I, IR, IOFFR, J, IK, JK, IAROW, IACOL,  &
                IRNUM, ICNUM, MB, NB, MEND, NEND,      &
                MRROW, MRCOL, MOFF, NOFF, IROW_GLOB, JCOL_GLOB
      INTEGER*4 NUMROC, ICEIL
!
!----------------------------------------------------------------------------------------!
!
!.... initialize
      IERR = 0
      M  = DESCA(3)    !number of rows of matrix
      N  = DESCA(4)    !number of columns of matrix
      MB = DESCA(5)    !row block size
      NB = DESCA(6)    !column block size
      IRNUM = NUMROC(M,MB,MYROW,0,NPROW)
      ICNUM = NUMROC(N,NB,MYCOL,0,NPCOL)
      IAROW = DESCA(7) !process row over which first row of A distributed
      IACOL = DESCA(8) !process col over which first col of A distributed
      IF (MB <= 0) THEN
         WRITE(*,*) 'get_diag_2dc: Serious error 1, mb <= 0',MB
         IERR = 1
         RETURN
      ENDIF
      IF (NB <= 0) THEN
         WRITE(*,*) 'get_diag_2dc: Serious error 2, nb <= 0',NB
         IERR = 2
         RETURN
      ENDIF
      MOFF = NPROW/MB !#process rows/block size
      NOFF = NPCOL/NB !#process columns/block size
      MEND = ICEIL(IRNUM, MB) + MOFF
      NEND = ICEIL(ICNUM, NB) + NOFF
      MRROW = MOD( NPROW+MYROW-IAROW, NPROW )
      MRCOL = MOD( NPCOL+MYCOL-IACOL, NPCOL )
      DIAG(1:M) = 0.D0 
!
!.... fill local matrices from global matrix
      JK = 1
      DO 290 IC = NOFF+1, NEND !loop on column blocks
         IOFFC = ((IC-1)*NPCOL+MRCOL)*NB
         DO 280 I = 1, NB !loop on columns in block 
            IF(JK .GT. ICNUM) GO TO 300 !done with rows
            IK = 1
            DO 260 IR = MOFF+1, MEND !loop on row blocks
               IOFFR = ((IR-1)*NPROW+MRROW)*MB
               DO 250 J = 1, MB !Loop on rows in block
                  IF( IK .GT. IRNUM ) GO TO 270 !out of local rows
                  IROW_GLOB = IOFFR + J !row of matrix
                  JCOL_GLOB = IOFFC + I
                  DIAG(IROW_GLOB) = ALOC(IK,JK) 
                  IK = IK + 1
  250          CONTINUE !loop on rows in block
  260       CONTINUE !loop on row blocks
  270       CONTINUE !out of local rows
            JK = JK + 1
  280    CONTINUE !loop on rows
  290 CONTINUE !loop on column blocks 
  300 CONTINUE !done
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE FILL_VEC_2DC(NWORKL,NWORK, MYROW, NPROW,  DESCB, B, &
                              BLOC,IERR) 
!
!     Takes a full vector then extracts the corresponding Scalapack local vector.
!
!     INPUT      MEANING
!     -----      ------- 
!     B          global vector B
!     DESCB      BLACS descriptor of vector B
!     MYROW      process row in BLACS grid
!     NPROW      number of process rows
!     NWORK      length of B
!     NWORKL     size of BLOC 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     BLOC       block cyclic vector B
! 
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: B(NWORK)
      INTEGER*4, INTENT(IN) :: DESCB(9), NWORK, NWORKL, MYROW, NPROW
      REAL*8, INTENT(OUT) :: BLOC(NWORKL)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 M, MB, IR, IRNUM, IAROW, MOFF, MEND, MRROW, IK, IOFFR, IROW_GLOB, J
      INTEGER*4 NUMROC, ICEIL
!
!----------------------------------------------------------------------------------------!
!
!.... initialize
      IERR = 0
      M  = DESCB(3)    !number of rows of matrix
      MB = DESCB(5)    !row block size
      IRNUM = NUMROC(M,MB,MYROW,0,NPROW)
      IAROW = DESCB(7) !process row over which first row of A distributed
      IF (MB <= 0) THEN
         WRITE(*,*) 'fill_vec_2dc: Serious error 1, mb <= 0',MB
         IERR = 1
         RETURN
      ENDIF
      MOFF = NPROW/MB !#process rows/block size
      MEND = ICEIL(IRNUM, MB) + MOFF
      MRROW = MOD( NPROW+MYROW-IAROW, NPROW )
!
!.... fill local vector 
      IK = 1
      DO 260 IR = MOFF+1, MEND !loop on row blocks
         IOFFR = ((IR-1)*NPROW+MRROW)*MB
         DO 250 J = 1, MB !Loop on rows in block
            IF ( IK .GT. IRNUM ) GO TO 270 !out of local rows
            IROW_GLOB = IOFFR + J !row of matrix
            BLOC(IK) = B(IROW_GLOB)
            IK = IK + 1
  250    CONTINUE !loop on rows in block
  260 CONTINUE !loop on row blocks
  270 CONTINUE !out of local rows
      RETURN
      END

      SUBROUTINE ASMBLE_JACSYS(LDJS,NOBSLMX,NA35,NFREQ, NDIM,NSRC,NREC, &
                               NPPGRP, FREQ, NOBSV,JACSYS,IERR) 
!
!     Assembles the global Jacobian onto the head process. From here we can distribute
!     to other processes on grid or factor. 
!
!     INPUT      MEANING
!     -----      ------- 

      IMPLICIT NONE  
      REAL*8, INTENT(IN) :: FREQ(NFREQ)
      INTEGER*4, INTENT(IN) :: LDJS, NOBSLMX, NA35, NFREQ, NDIM, NSRC, NREC, NPPGRP  
      REAL*4, INTENT(OUT) :: JACSYS(LDJS,na35) 
      INTEGER*4, INTENT(OUT) :: NOBSV(NFREQ*NSRC) 
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      CHARACTER(80) FLNAME
      CHARACTER(12) CFREQ
      CHARACTER(5) CID, CSRC 
      REAL*4, ALLOCATABLE :: CJACL(:,:), RJACL(:,:) 
      INTEGER*4, ALLOCATABLE :: MYDEST(:)
      INTEGER*4 LDJAC, IUNIT, IFREQ, IREC, ISRC, IPPGRP, IBEG, IEND, LDOBS, &
                NOBS_IN,NNPGL_IN, K, NSUM
      LOGICAL*4 LOBS 
      PARAMETER(IUNIT = 30) 

      LDJAC = NOBSLMX/2 !load real and imaginary
      ALLOCATE(RJACL(LDJAC,NA35)) 
      ALLOCATE(CJACL(LDJAC,NA35)) 
      ALLOCATE(MYDEST(NA35)) 
      RJACL(:,:) = 0.0
      CJACL(:,:) = 0.0
      MYDEST(:) = 0 
      LDOBS = 1  
      IBEG = 1
      IEND = 0
      K = 0
      NSUM = 0 
      DO 1 IFREQ=1,NFREQ
         CFREQ(1:12) = ' '
         WRITE(CFREQ,'(F12.5)') FREQ(IFREQ)
         CFREQ = ADJUSTL(CFREQ)
         DO 2 ISRC=1,NSRC
            K = K + 1
            CSRC(1:5) = ' '
            WRITE(CSRC,'(I5)') ISRC
            CSRC = ADJUSTL(CSRC)
            NOBSV(K) = 0 
            DO 3 IPPGRP=1,NPPGRP
               FLNAME(1:80) = ' '
               CID(1:5) = ' '
               WRITE(CID,'(I5)') IPPGRP - 1
               CID = ADJUSTL(CID)
               FLNAME = './scratch/jac_'//TRIM(CID)//'-'// &
                        TRIM(CFREQ)//'-'//TRIM(CSRC)//'.dat'
               FLNAME = ADJUSTL(FLNAME)
               WRITE(*,*) 'asmble_jacsys: Reading:',TRIM(FLNAME)
               OPEN(FILE=TRIM(FLNAME),UNIT=IUNIT,STATUS='OLD', &
                    FORM='UNFORMATTED',IOSTAT=IERR)
               CALL LOAD_JAC(LDOBS,LDJAC, NDIM,NREC, IUNIT,.FALSE.,.TRUE., 0, & 
                             NOBS_IN,NNPGL_IN,LOBS,RJACL,CJACL,MYDEST, IERR) 
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'asmble_jacsys: Error loading Jacobian!'
                  RETURN
               ENDIF
               IF (IPPGRP == 1) NOBSV(K) = NOBS_IN
               CALL FILL_CJACSYS(LDJS,LDJAC, NNPGL_IN,NOBS_IN,IBEG, MYDEST, RJACL,CJACL, &
                                 IEND,JACSYS)
    3       CONTINUE
            IBEG = IEND + 1 !next source 
            NSUM = NSUM + NOBSV(K)*2 
    2    CONTINUE
    1 CONTINUE 
      IF (NSUM /= IEND) THEN
         print *, nsum,iend   
         WRITE(*,*) 'asmble_jacsys: Warning may be an observation mismatch' 
      ENDIF
      DEALLOCATE(RJACL)
      DEALLOCATE(CJACL)
      DEALLOCATE(MYDEST)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE LOAD_JAC(LDOBS,LDJAC, NDIM,NREC, IUNIT,LOBSIN,LDSTIN, MYID,  &
                          NOBS_IN,NNPGL_IN,LOBS,RJAC,CJAC,MYDEST, IERR) 
!
!     Loads the Jacobians, and optionally whether the observation exists at a 
!     receiver/component pair, and the destination in the global inversion problem
!
!.... variable declarations
      INTEGER*4, INTENT(IN) :: IUNIT 
      LOGICAL*4, INTENT(IN) :: LOBSIN, LDSTIN 
      REAL*4, INTENT(OUT) :: RJAC(LDJAC,*), CJAC(LDJAC,*)
      LOGICAL*4, INTENT(OUT) :: LOBS(LDOBS,*)
      INTEGER*4, INTENT(OUT) :: MYDEST(*), IERR, NOBS_IN, NNPGL_IN 
!.... local variables
      LOGICAL*4 LOBST
      INTEGER*4 NREC_IN, NDIM_IN, I, IOBS, INPGL  
!
!----------------------------------------------------------------------------------------!
! 
      IERR = 0 
      READ(IUNIT,IOSTAT=IERR) NOBS_IN,NNPGL_IN, NDIM_IN, NREC_IN
      IF (IERR /= 0) THEN
         WRITE(*,*) 'load_jac: Error reading header',MYID
         CLOSE(IUNIT)
         RETURN
      ENDIF
      IF (NDIM_IN /= NDIM) THEN
         WRITE(*,*) 'load_jac: Serious error, ndim_in /= ndim',MYID,NDIM_IN,NDIM
         IERR = 1
         CLOSE(IUNIT)
         RETURN
      ENDIF
      IF (NREC_IN /= NREC) THEN
         WRITE(*,*) 'load_jac: Serious error, nrec_in /= nrec',MYID,NREC_IN,NREC
         IERR = 1
         CLOSE(IUNIT)
         RETURN
      ENDIF
      IF (LOBSIN) THEN
         READ(IUNIT,IOSTAT=IERR) ((LOBS(I,IREC),I=1,NDIM),IREC=1,NREC_IN)
      ELSE
         READ(IUNIT,IOSTAT=IERR) ((LOBST,I=1,NDIM),IREC=1,NREC_IN)
      ENDIF
      IF (IERR /= 0) THEN
         WRITE(*,*) 'load_jac: Error reading LOBS',MYID
         CLOSE(IUNIT)
         RETURN
      ENDIF
      READ(IUNIT,IOSTAT=IERR)((RJAC(IOBS,INPGL),IOBS=1,NOBS_IN),INPGL=1,NNPGL_IN)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'load_jac: Error reading real Jacobian!',MYID
         CLOSE(IUNIT)
         RETURN
      ENDIF
      READ(IUNIT,IOSTAT=IERR)((CJAC(IOBS,INPGL),IOBS=1,NOBS_IN),INPGL=1,NNPGL_IN)
      IF (IERR /= 0) THEN 
         WRITE(*,*) 'load_jac: Error reading complex Jacobian!',MYID
         CLOSE(IUNIT)
         RETURN
      ENDIF
      IF (LDSTIN) THEN
         READ(IUNIT,IOSTAT=IERR)(MYDEST(INPGL),INPGL=1,NNPGL_IN)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'load_jac: Error reading inversion destinations!',MYID
            CLOSE(IUNIT)
            RETURN
         ENDIF
      ENDIF
      CLOSE(IUNIT) !done with file 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE ICNOBS_JMAT(NFREQ,NSRC, MYID, FREQ, NOBSLMX,NOBS,IERR) 
!
!     INPUT     MEANING
!     -----     ------- 
!     FREQ      frequency list  
!     MYID      process ID (prevents MPI from grabbing same file and unit number)
!     NFREQ     number of frequencies in inversion block
!     NSRC      number of sources to invert
!
!     OUTPUT    MEANING
!     ------    ------- 
!     IERR      error flag
!     NOBS      number of observations
!     NOBSLMX   max number of observations for all local Jacobians
!     
!.... variable declarations
      IMPLICIT NONE 
      REAL*8, INTENT(IN) :: FREQ(NFREQ) 
      INTEGER*4, INTENT(IN) :: NFREQ, NSRC, MYID
      INTEGER*4, INTENT(OUT) :: NOBSLMX, NOBS, IERR
!.... variable declarations
      CHARACTER(80) FLNAME 
      CHARACTER(12) CFREQ
      CHARACTER(5) CID, CSRC 
      INTEGER*4 NOBS_IN, IFREQ, ISRC, IPPGRP, IUNIT
!
!----------------------------------------------------------------------------------------!
!
      IERR = 0
      NOBS = 0 
      NOBSLMX = 0
      DO 1 IFREQ=1,NFREQ
         CFREQ(1:12) = ' ' 
         WRITE(CFREQ,'(F12.5)') FREQ(IFREQ)
         CFREQ = ADJUSTL(CFREQ) 
         DO 2 ISRC=1,NSRC
            CSRC(1:5) = ' '
            WRITE(CSRC,'(I5)') ISRC
            CSRC = ADJUSTL(CSRC)
            IPPGRP = 1 !only need one group 
            IUNIT = 10 + MYID 
            FLNAME(1:80) = ' '
            CID(1:5) = ' ' 
            WRITE(CID,'(I5)') IPPGRP - 1
            CID = ADJUSTL(CID)
            FLNAME = './scratch/jac_'//TRIM(CID)//'-'// &
                     TRIM(CFREQ)//'-'//TRIM(CSRC)//'.dat'
            FLNAME = ADJUSTL(FLNAME)
            OPEN(FILE=TRIM(FLNAME),UNIT=IUNIT,STATUS='OLD', &
                 FORM='UNFORMATTED',IOSTAT=IERR)
            READ(IUNIT,IOSTAT=IERR) NOBS_IN
            IF (IERR /= 0) THEN
               WRITE(*,*) 'icnobs_jmat: Error reading:',TRIM(FLNAME)
               CLOSE(IUNIT)
               RETURN
            ENDIF
            CLOSE(IUNIT)
            NOBS = NOBS + NOBS_IN
            NOBSLMX = MAX0(NOBSLMX,NOBS_IN)
    2    CONTINUE
    1 CONTINUE 
      NOBS = NOBS*2 !have complex equations as well
      NOBSLMX = NOBSLMX*2 !have complex equations as well
      RETURN
      END  
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE FILL_CJACSYS(LDJS,LDJL, NNPGL_IN,NOBS_IN,IBEG, MYDEST, RJACL,CJACL,  &
                              IEND,JACSYS)
!
!     Fills in the Jacobian which is ordered 
!         [ Re{J}_1 ] 
!         [ Im{J}_1 ] 
!         [ Re{J}_2 ]
!         [ Im{J}_2 ]
!         [    .    ]
!         [    .    ]
!         [    .    ]
!
!     INPUT      MEANING
!     -----      ------- 
!     CJACL      complex local Jacobian read from disk
!     IBEG       index to begin filling global Jacobian 
!     LDJL       leading dimension of local Jacobian from disk 
!     LDJS       leading dimension for Jacobian system 
!     MYDEST     maps local Jacobian inversion point to global inversion point
!     NNPGL_IN   number of local inversion variables in this Jacobian 
!     NOBS_IN    number of obsrevations in local jacobian
!     RJACL      real local Jacobian read from disk 
!
!     OUTPUT     MEANING
!     ------     -------
!     IEND       index Jacobian fill finished on 
!     JACSYS     updated Jacobian 
!
!.... variable declarations 
      IMPLICIT NONE
      REAL*4, INTENT(INOUT) :: JACSYS(LDJS,*)
      REAL*4, INTENT(IN) :: RJACL(LDJL,NNPGL_IN), CJACL(LDJL,NNPGL_IN)
      INTEGER*4, INTENT(IN) :: MYDEST(NNPGL_IN), LDJS,LDJL, NNPGL_IN, NOBS_IN, IBEG
      INTEGER*4, INTENT(OUT) :: IEND 
      INTEGER*4 ILOC, JLOC, IDEST, IOBS, INPGL 
!
!----------------------------------------------------------------------------------------!
!
      ILOC = IBEG - 1 
      JLOC = ILOC + NOBS_IN 
      DO 1 IOBS=1,NOBS_IN
         ILOC = ILOC + 1
         JLOC = JLOC + 1
         DO 2 INPGL=1,NNPGL_IN
            IDEST = MYDEST(INPGL)
            JACSYS(ILOC,IDEST) = RJACL(IOBS,INPGL)
            JACSYS(JLOC,IDEST) = CJACL(IOBS,INPGL)
    2    CONTINUE
    1 CONTINUE   
      IEND = JLOC 
      RETURN
      END
