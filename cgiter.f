      SUBROUTINE CGITER_CF(MYID,MYNID,MASTER,MYHD_COMM,MYSLV_COMM,
     ;                     MGCOMM,IPGROUP,NPGROUPS, LGUESS,  
     ;                     FRQ,SRC,INV, X,IERR)   
!
!     Applies the reverse communication complex conjugate gradient 
!     method from CERFACS
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO)  INV
      TYPE (FRQ_INFO)  FRQ
      TYPE (SRC_INFO)  SRC
      INTEGER*4, INTENT(IN) :: MYID, MYNID, MASTER, MYHD_COMM,
     ;                         MYSLV_COMM, MGCOMM,IPGROUP,NPGROUPS
      LOGICAL*4, INTENT(IN) :: LGUESS
      REAL*4, INTENT(INOUT) :: X(*) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: WORK(:)
      REAL*4, ALLOCATABLE :: TX(:), TZ(:)
      REAL*8 CNTL(3), RINFO(3), DDOT 
      INTEGER*4 IRC(4), ICNTL(7), INFO(3), ICOLX, ICOLY, ICOLZ, IERR1, 
     ;          LWORK, NLOC, MPIERR 
      INTEGER*4 MATVEC, PRECONLEFT, DOTPROD, REVCOM, MAXIT
      PARAMETER(MATVEC = 1, PRECONLEFT = 2, DOTPROD = 3)
      !PARAMETER(MAXIT = 500) 
!  
!----------------------------------------------------------------------!
!
!.... initialize 
      IERR = 0 
      ICOLX = 1
      ICOLY = 1
      ICOLZ = 2
      CALL MPI_BARRIER(MGCOMM,MPIERR)
      IF (MYID.EQ.MASTER) THEN
         N = inv%NA35 
         MAXIT = N
         CALL INIT_DCG(ICNTL,CNTL)
         INFO(1:3) = 0
         icntl(1) = 6
         icntl(2) = 6
         icntl(3) = 6
         ICNTL(4) = 1 !left pre-conditioning
         IF (LGUESS) THEN
            ICNTL(5) = 1 !initial guess
         ELSE
            ICNTL(5) = 0 !no initial guess
         ENDIF
         ICNTL(6) = MAXIT !max number of iterations
         CNTL(1) = 1.19209290E-05 !stopping criteria, 100*macheps
         NLOC = N !CG is handled seqeuntially
         IF (ICNTL(7).EQ.0) THEN
            LWORK = 6*N + 1 
         ELSE
            LWORK = 6*N + 1 + 2*(MAXIT + 1)
         ENDIF
         ICOLX = 1 
         ICOLY = 1 + N
         ALLOCATE(TX(N))
         ALLOCATE(TZ(N))
         ALLOCATE(WORK(LWORK))
         IF (LGUESS) THEN
            WORK(1:N) = DBLE(X(1:N)) 
         ELSE
            WORK(1:N) = 0.D0
         ENDIF
         !CALL SCOPY(N,B,1,WORK(ICOLY),1) !intialize RHS
         !CALL SCOPY(N,inv%GRAD(1:N),1,WORK(ICOLY),1) !initialize RHS
         WORK(N+1:2*N) = DBLE(inv%GRAD(1:N))
      ELSE
         ALLOCATE(WORK(2)) 
         WORK(1:2) = 0.D0
      ENDIF

   10 CONTINUE 
      IF (MYID.EQ.MASTER) THEN
         CALL DRIVE_DCG(N,NLOC,LWORK,WORK,
     ;                  IRC,ICNTL,CNTL,INFO,RINFO)  
         IF (INFO(1).EQ.-1) THEN
            WRITE(*,*) 'cgiter_cf: Error n < 1'
         ELSEIF (INFO(1).EQ.-2) THEN
            WRITE(*,*) 'cgiter_cf: Error LWORK is too small!'
         ELSEIF (INFO(1).EQ.-3) THEN
            WRITE(*,*) 'cgiter_cf: Error iteration limit reached!'
         ELSEIF (INFO(1).EQ.-4) THEN
            WRITE(*,*) 'cgiter_cf: Error preconditioner not set!'
         ENDIF
         REVCOM = IRC(1) 
         ICOLX = IRC(2)
         ICOLY = IRC(3)
         ICOLZ = IRC(4) 
      ENDIF
      CALL MPI_BCAST(REVCOM,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR) 
      !print *, revcom,ICOLX,ICOLZ
      IF (REVCOM.EQ.MATVEC) THEN !matrix vector multiply
         CALL GNMFREE(MYID,MYNID,MASTER,MYHD_COMM,MYSLV_COMM, 
     ;                IPGROUP,NPGROUPS, FRQ,SRC,INV, WORK(ICOLX),
     ;                WORK(ICOLZ),IERR1) 
         IF (IERR1.NE.0) THEN 
            WRITE(*,*) 'cgiter_cf: Error occured on',MYID,' in gnmfree'
            RETURN
         ENDIF
         CALL MPI_ALLREDUCE(IERR1,IERR,1,MPI_INTEGER,MPI_MAX, 
     ;                      MGCOMM,MPIERR)
         IF (IERR.EQ.0) GOTO 10
      ELSEIF (REVCOM.EQ.PRECONLEFT) THEN  !left preconditioning
         IF (MYID.EQ.MASTER) THEN
            IF (.NOT.inv%LPGRAD) THEN
               CALL DCOPY(N,WORK(ICOLX),1,WORK(ICOLZ),1)
            ELSE
               TX(1:N) = WORK(ICOLX:ICOLX+N-1)
               CALL PREGRAD(inv%NHSIZE,inv%NA35, inv%NNPINV,inv%NVINV,  
     ;                      inv%IPIVH,inv%GRADPC,TX, TZ,IERR) 
               WORK(ICOLZ:ICOLZ+N-1) = TZ(1:N)
               IF (IERR.NE.0) 
     ;         WRITE(*,*) 'cgiter_cf: Error applying preconditioner'
            ENDIF
         ENDIF
         CALL MPI_BCAST(IERR,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR)
         IF (IERR.EQ.0) GOTO 10 
      ELSEIF (REVCOM.EQ.DOTPROD) THEN !dot product
         IF (MYID.EQ.MASTER) 
     ;   WORK(ICOLZ) = DDOT(N,WORK(ICOLX),1,WORK(ICOLY),1)
         GOTO 10
      ENDIF

      IF (MYID.EQ.MASTER) THEN
         IF (ICNTL(7).EQ.1) THEN
            WRITE(*,*) 'cgiter: Estimation of smallest eigenvalue:',
     ;                 RINFO(2)
            WRITE(*,*) 'cgiter: Estimation of largest eigenvalue:', 
     ;                 RINFO(3)
         ENDIF
         !CALL SCOPY(N,WORK(ICOLX),1,X,1)
         X(1:N) = SNGL(WORK(ICOLX:ICOLX+N-1)) 
         IF (ALLOCATED(TX)) DEALLOCATE(TX)
         IF (ALLOCATED(TZ)) DEALLOCATE(TZ) 
      ENDIF
      DEALLOCATE(WORK) 
      RETURN
      END 
 
      SUBROUTINE CGITER_GN(MYID,MYNID,MASTER,MYHD_COMM,MYSLV_COMM, 
     ;                     MGCOMM,IPGROUP,NPGROUPS, LGUESS,  
     ;                     FRQ,SRC,INV, XIN,IERR) 
!
!     Applies the pre-conditioned conjugate gradient method to the 
!     iterative problem.  We initialize x to 0 on the recommendation 
!     of Nocedal and Wright pg. 139 and Demmel page 317. For more see
!     Netlib templates and A Set of Conjugate Gradient Routines for 
!     Real and Complex Arithmetics - Cerfacs 
!
!     INPUT      MEANING
!     -----      ------- 
!     B          gradient, right hand side of Ax = b
!     MASTER     master process ID
!     MYCOMM     MPI communicatr 
!     MYID       process ID on MGCOMM 
!     N          order of problem 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     X          search direction, solution of Ax = b.  But only real
!
!.... variable declarations
      implicit none
      INCLUDE 'fwd_struc.h'
      INCLUDE 'mpif.h'
      TYPE (INV_INFO)  INV
      TYPE (FRQ_INFO)  FRQ
      TYPE (SRC_INFO)  SRC
      INTEGER*4, INTENT(IN) :: MYID, MYNID, MASTER, MYHD_COMM,
     ;                         MYSLV_COMM, MGCOMM,IPGROUP,NPGROUPS
      LOGICAL*4, INTENT(IN) :: LGUESS 
      REAL*4, INTENT(INOUT) :: XIN(*) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: R(:), P(:), Q(:), Z(:), B(:), X(:)
      REAL*8 ZERO, ONE, RHO, RHO1, BETA, ALPHA, ALPHAD
      REAL*8 BNRM2, ERROR, DNRM2, TOL, GOLD, DDOT
      INTEGER*4 I, K, MYDEST, MPIERR, MAXIT, MMUL, N 
      PARAMETER(TOL  = 1.19209290E-05) !mach eps*100
      PARAMETER(GOLD = 1.6180339887498949) !(1 + sqrt(5))/2 
      PARAMETER(MAXIT = 4000) 
      PARAMETER(ZERO = 0.D0, ONE = 1.D0)
!
!----------------------------------------------------------------------!
!
!.... initialize
      IERR = 0
      RHO  = ZERO
      RHO1 = ZERO
      IF (MYID.EQ.MASTER) THEN
         N = inv%NA35
         ALLOCATE(B(N)) 
         ALLOCATE(X(N)) 
         !CALL SCOPY(N,inv%GRAD,1,B,1) 
         B(1:N) = inv%GRAD(1:N) 
         BNRM2 = DNRM2(N,B,1) 
         IF (BNRM2.EQ.0.D0) THEN
            WRITE(*,*) 'cgiter: The right hand side is zero!'
            CALL DSCAL(N,ZERO,X,1)
            IERR =-1
            GOTO 1500 
         ENDIF 
         ALLOCATE(R(N))  !residual
         ALLOCATE(P(N))  !search direction
         ALLOCATE(Q(N))  !result A p_k
         ALLOCATE(Z(N))  !pre-conditioned solution of inv(~H) r 
         IF (.NOT.LGUESS) THEN !zero guess, 139 Nocedal and Wright
            X(1:N) = DBLE(XIN(1:N)) 
            CALL DSCAL(N,ZERO,X,1) 
            CALL DCOPY(N,B,1,R,1)   !Residual r = b - A x = b
            MMUL = 0
         ELSE 
            MMUL = 1
         ENDIF !otherwise X
         CALL DCOPY(N,B,1,P,1) 
      ELSE
         ALLOCATE(P(1))
         ALLOCATE(Q(1))
      ENDIF
      CALL MPI_BCAST(IERR,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR)
      IF (IERR.NE.0) GOTO 2000
      CALL MPI_BCAST(MMUL,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR)
!
!.... if we have an initial guess, need a matrix vector multiply
      IF (MMUL.EQ.1) THEN
         CALL GNMFREE(MYID,MYNID,MASTER,MYHD_COMM,MYSLV_COMM,
     ;                IPGROUP,NPGROUPS, FRQ,SRC,INV, X,
     ;                Q,IERR)
         IF (IERR.NE.0) GOTO 2000 
         IF (MYID.EQ.MASTER) R(1:N) = B(1:N) - Q(1:N)
      ENDIF
      IF (MYID.EQ.MASTER) THEN 
         ERROR = 1.0
         print *, 'beginning with:',minval(dabs(r)),maxval(dabs(r))
         print *, 'beginning with:',minval(dabs(b(1:n))),
     ; maxval(abs(b(1:n)))
      ENDIF
 1500 CONTINUE !break ahead for an error 
      CALL MPI_BCAST(IERR,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR)
      IF (IERR.NE.0) GOTO 2000 
!
!.... iterative loop in preconditioned CG
      DO 1000 K=1,MAXIT
         IF (MYID.EQ.MASTER) THEN
            IF (inv%LPGRAD) THEN !solve M z_k = r_k
               CALL PREGRAD(inv%NHSIZE,inv%NA35, inv%NNPINV,inv%NVINV,
     ;                      inv%IPIVH,inv%GRADPC,R,  Z,IERR)
               IF (IERR.NE.0) THEN
                  WRITE(*,*) 'cgiter: An preconditioning error occurred'
                  MYDEST = 2000
                  GOTO 1200 !break ahead to broadcast error 
               ENDIF
            ELSE !preconditioner identity
               Z(1:N) = R(1:N) !CALL DCOPY(N,R,1,Z,1) 
            ENDIF
            RHO = DOT_PRODUCT(Z,R) !DDOT(N,Z,1,R,1)
            IF (K.EQ.1) THEN
               P(1:N) = Z(1:N)
            ELSE 
               BETA = RHO/RHO1 
               P(1:N) = Z(1:N) + BETA*P(1:N) 
c              CALL DAXPY(N, BETA,Z,1,P,1)
c              DO 1 I=1,N !update search
c                 P(I) = Z(I) + BETA*P(I)
c   1          CONTINUE
            ENDIF
         ENDIF 
 1200    CONTINUE
         CALL MPI_BCAST(MYDEST,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR)
         IF (MYDEST.EQ.2000) GOTO 2000
         CALL GNMFREE(MYID,MYNID,MASTER,MYHD_COMM,MYSLV_COMM, 
     ;                IPGROUP,NPGROUPS, FRQ,SRC,INV, P,
     ;                Q,IERR)  !q = A p_k
         IF (IERR /= 0) THEN
            WRITE(*,*) 'cgiter: Error calling gnmfree on process!',MYID
            GOTO 1000
         ENDIF
         MYDEST = 1000
         IF (MYID.EQ.MASTER) THEN
            ALPHAD = DOT_PRODUCT(P,Q) !ALPHAD = DDOT(N,P,1,Q,1) 
            ALPHA = RHO/ALPHAD
            X(1:N) = X(1:N) + ALPHA*P(1:N) !CALL DAXPY(N, ALPHA,P,1,X,1) !x = x + alpha p; approx
            R(1:N) = R(1:N) - ALPHA*Q(1:N) !CALL DAXPY(N,-ALPHA,Q,1,R,1) !r = r - alpha q; residual
            ERROR = DSQRT(DOT_PRODUCT(R,R))/BNRM2 !DNRM2(N,R,1)/BNRM2
            print *, 'error:',error, dsqrt(dot_product(r,r))
            RHO1 = RHO
            IF (ERROR.LT.TOL) MYDEST = 2000
         ENDIF 
         CALL MPI_BCAST(MYDEST,1,MPI_INTEGER, MASTER,MGCOMM,MPIERR)
         IF (MYDEST.EQ.2000) GOTO 2000 
 1000 CONTINUE !loop on max iterations
      IF (MYID.EQ.MASTER) THEN
         WRITE(*,*) 'cgiter: Max iteration limit reached!'
         !IERR = 1 
      ENDIF
 2000 CONTINUE !we've converged (or converged to an error :-/)
      IF (MYID.EQ.MASTER .AND. IERR.EQ.0) THEN
         WRITE(*,*) 'cgiter: Converged after iteration',K
      ENDIF
      IF (IERR.EQ.-1) IERR = 0
      IF (MYID.EQ.MASTER) THEN
         XIN(1:N) = SNGL(X(1:N))
         IF (ALLOCATED(R)) DEALLOCATE(R)
         IF (ALLOCATED(Z)) DEALLOCATE(Z)
         IF (ALLOCATED(X)) DEALLOCATE(X)
      ENDIF
      IF (ALLOCATED(P)) DEALLOCATE(P)
      IF (ALLOCATED(Q)) DEALLOCATE(Q)
      RETURN
      END

