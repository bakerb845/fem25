      MODULE EQN439_DAT
         REAL*8, ALLOCATABLE :: Z(:,:)  !eigenvectors
         REAL*8, ALLOCATABLE :: W(:)    !eigenvalues 
         REAL*8, ALLOCATABLE :: GLOC(:) !vector g in matrix vector multiply
         REAL*8, ALLOCATABLE :: YLOC(:) !result y = Zg 
         REAL*8, ALLOCATABLE :: PLOC(:) !result p = Z^T L Z g
         REAL*8 DELTA_PASS              !trust region size
         INTEGER*4 DESCZ(9)             !descriptor for eigenvector matrix
         INTEGER*4 DESCG(9)             !descriptor for vector g
         INTEGER*4 ICTXT_PASS           !BLACS context  
         INTEGER*4 NPROW_PASS           !number of rows in process grid
         INTEGER*4 NPCOL_PASS           !number of columns in process grid
         INTEGER*4 MYROW                !my process row 
         INTEGER*4 MYCOL                !my process column
         INTEGER*4 N                    !size of z 
         LOGICAL*4 LDONE                !done working
         LOGICAL*4 LVERB                !(0,0) rank prints updates in optimization
      END MODULE EQN439_DAT
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE EST_SIGMA(ICTXT,NPROW,NPCOL, LDH,LDC,LSURF,K, DESCH,DESCC, &
                           HLOC,COVM, IERR) 
!
!     Estimates the spatial correlation length by matching the eigenspectrum of
!     C^{-1} and adj(J) J in H

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: HLOC(LDH,*), COVM(LDC,*) 
      INTEGER*4, INTENT(IN) :: DESCH(9), DESCC(9), ICTXT, LDH, LDC, NPROW, NPCOL, K
      LOGICAL*4, INTENT(IN) :: LSURF
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: A(:,:), Z(:,:), W(:), WORK(:)
      CHARACTER(80) FLNAME
      CHARACTER(5) CK
      REAL*8 WORK8 
      INTEGER*4 DESCZ(9), N, LDA, LWORK, LIWORK, MYROW, MYCOL, I
!.... BLACS stuff
      !INTEGER*4 DESCA(9)   !A matrix descriptor
      INTEGER*4 NUMROC     !calculates NUmber of Rows Or Columns
      INTEGER*4 INFO       !error flag from ScaLapack 
      INTEGER*4 M          !number of rows and columns of Jacobian
      INTEGER*4 MH         !number of local rows in h
      INTEGER*4 NH         !number of local columns in h
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
!
!----------------------------------------------------------------------------------------!
!
!.... set A 
      IERR = 0
      LDA = DESCH(LLD_)
      M = DESCH(M_)
      N = DESCH(N_)
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
      MH = NUMROC(M, DESCH(MB_),MYROW, 0,NPROW) !number of rows in local H matrix
      NH = NUMROC(N, DESCH(NB_),MYCOL, 0,NPCOL) !number of columns in local H matrix
      ALLOCATE(A(LDA,MAX(NH,1)),STAT=IERR)
      A(:,:) = 0.D0
      IF (IERR /= 0) THEN
         WRITE(*,*) 'est_sigma: Error setting space for A on grid',MYROW,MYCOL
         RETURN
      ENDIF
      ALLOCATE(Z(1,1)) !LDA,MAX(NH,1)),STAT=IERR)
      Z(:,:) = 0.D0
      IF (IERR /= 0) THEN
         WRITE(*,*) 'est_sigma: Error setting space for Z on grid',MYROW,MYCOL
         RETURN
      ENDIF
      DESCZ(1:9) = DESCH(1:9)
      CALL COPY_HMAT_UL('U',ICTXT,NPROW,NPCOL, LDH,N,DESCH, HLOC,A)

      LWORK =-1  !space inquiry 
      LIWORK =-1 !space inquiry
      ALLOCATE(W(N)) !eigenvalues
      CALL PDSYEV('N','U',N,A,1,1,DESCH, W,Z,1,1, &
                   DESCZ, WORK8,LWORK, INFO)
      IF (INFO /= 0) THEN
          WRITE(*,*) 'est_sigma: Error initializing pdsyev',INFO
          IERR = 1
          RETURN
      ENDIF
      LWORK = INT(WORK8)
      ALLOCATE(WORK(LWORK))
      CALL PDSYEV('N','U',N,A,1,1,DESCH, W,Z,1,1, &
                   DESCZ, WORK,LWORK, INFO)
      IF (INFO /= 0) THEN
         WRITE(*,*) 'est_sigma: Error calling pdsyev',INFO
         IERR = 1
         RETURN
      ENDIF
!     IF (ALLOCATED(WORK))  DEALLOCATE(WORK)
!     IF (ALLOCATED(Z))     DEALLOCATE(Z) 
!
!.... save the eigenspectrum, b/c it is informative
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         FLNAME(1:80) = ' '
         CK(1:5) = ' '
         WRITE(CK,'(I5)') K
         CK = ADJUSTL(CK)
         IF (LSURF) THEN
            FLNAME = 'eigen_spc_srf-'//TRIM(CK)//'.txt'
         ELSE
            FLNAME = 'eigen_spc_bdy-'//TRIM(CK)//'.txt'
         ENDIF
         FLNAME = ADJUSTL(FLNAME)
         OPEN(UNIT=45,FILE=TRIM(FLNAME))
         DO 1 I=1,N
            WRITE(45,*) I, W(I)
    1    CONTINUE
         CLOSE(45)
      ENDIF
!
!.... get the eigenspectrum of C
      DEALLOCATE(A)
      LDA = DESCC(LLD_) 
      MH = NUMROC(M, DESCC(MB_),MYROW, 0,NPROW) !number of rows in local H matrix
      NH = NUMROC(N, DESCC(NB_),MYCOL, 0,NPCOL) !number of columns in local H matrix
      ALLOCATE(A(LDA,MAX(NH,1)),STAT=IERR)
      DESCZ(1:9) = DESCC(1:9)
      A(:,:) = 0.D0
      CALL COPY_HMAT_UL('U',ICTXT,NPROW,NPCOL, LDC,N,DESCC, COVM,A)
      CALL PDSYEV('N','U',N,A,1,1,DESCC, W,Z,1,1, &
                   DESCZ, WORK,LWORK, INFO)
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         FLNAME(1:80) = ' '
         IF (LSURF) THEN
             FLNAME = 'covar_spc_srf.txt'
         ELSE
             FLNAME = 'covar_spc_bdy.txt'
         ENDIF 
         FLNAME = ADJUSTL(FLNAME)
         OPEN(UNIT=45,FILE=TRIM(FLNAME))
         DO 2 I=1,N
            WRITE(45,*) I,W(I)
    2    CONTINUE
         CLOSE(45)
      ENDIF
 
      IF (ALLOCATED(WORK)) DEALLOCATE(WORK)
      IF (ALLOCATED(Z))    DEALLOCATE(Z)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE TARANTOLA_394_TR(ICTXT,NPROW,NPCOL,  K,LDH,DESCH,                   &
                                  NNPG,NVINV,MASKG, XLOCS,ZLOCS,                     &
                                  KERNEL,VELVAR, SIGMAX,SIGMAZ,                      &
                                  DELTA,HLOC, G,  RLAM,SEARCH, IERR) 
!
!     Trust region implementation of Tarantola's (3.94) Gauss-Newton equation: 
!        A = G^T inv(C_D) G + inv(C_M) = B + inv(C_M)
!        y = G^T inv(C_D) r + inv(C_M) (m - m0) = g + inv(C_M) (m - m0)
!     This can be rewritten as 
!        A = C_M G^T inv(C_D) G + lambda I = C_M B + lambda I
!        y = C_M G^T inv(_D) r + (m - m0) = C_M g + (m - m0)
!
!     Now we identify a lambda so that we are in the trust region
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: HLOC(LDH,*), XLOCS(NNPG), ZLOCS(NNPG), G(*), VELVAR,  &
                            SIGMAX, SIGMAZ, DELTA
      INTEGER*4, INTENT(IN) :: MASKG(NNPG), DESCH(9), ICTXT, NPROW, NPCOL, LDH, K, &
                               NNPG, NVINV, KERNEL 
      REAL*8, INTENT(OUT) :: SEARCH(*),RLAM
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: A(:,:), COVM(:,:), GLOC(:), YLOC(:), GWORK(:)
      REAL*8 D1, D2, COV, SMALL
      REAL*8, PARAMETER :: ONE  = 1.D0
      REAL*8, PARAMETER :: ZERO = 0.D0
      REAL*8, PARAMETER :: SQRT2PI = 2.5066282746310002D0 !sqrt(2*pi)
      INTEGER*4 LDA, INPINV, IA35, JNPINV, JA35, INPG, JNPG, IIA,JJA, IAROW,IACOL, &
                IVINV, JVINV
      INTEGER*4 DESCA(9)   !A matrix descriptor
      INTEGER*4 DESCC(9)   !model covariance matrix descriptor 
      INTEGER*4 DESCG(9)   !local gradient descriptor
      INTEGER*4 M,N        !H [m x n]
      INTEGER*4 MH,NH      !Number of local rows/columsn in H
      INTEGER*4 NUMROC     !calculates NUmber of Rows Or Columns
      INTEGER*4 INFO       !error flag from ScaLapack
      INTEGER*4 MYROW      !my process row
      INTEGER*4 MYCOL      !my process column
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &    
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )  
      logical*4, parameter :: ldebug = .true.
      real*8 z, work8
      real*8, allocatable :: w(:), work(:) 
      integer*4 lwork
!
!----------------------------------------------------------------------------------------!
!
!.... get sizes and recall our blacs information
      IERR = 0
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
      M = DESCH(M_)
      N = DESCH(N_) 
      IF (M /= N) THEN 
         WRITE(*,*) 'tarantola_394_tr: Error m /= n!'
         IERR = 1
         RETURN
      ENDIF
      LDA = DESCH(LLD_) 
      MH = NUMROC(M, DESCH(MB_),MYROW, 0,NPROW) !number of rows in local H matrix
      NH = NUMROC(N, DESCH(NB_),MYCOL, 0,NPCOL) !number of columns in local H matrix
      CALL ICOPY(9,DESCH,1,DESCA,1)
      CALL ICOPY(9,DESCH,1,DESCC,1)
!
!.... set the model covariance matrix
      IF (MYROW == 0 .AND. MYCOL == 0) &
      WRITE(*,*) 'tarantola_394_tr: Setting model covariance matrix...'
      ALLOCATE(COVM(LDA,MAX(1,NH)),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'tarantola_394_tr: Error couldnt set space for covm'
         RETURN
      ENDIF
      SMALL = EPSILON(1.D0)
      COVM(:,:) =-1.D0
      ALLOCATE(A(LDA,MAX(1,NH)),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'tarantola_394_tr: Error couldnt set space for A'
         RETURN 
      ENDIF
      DO 1 INPG=1,NNPG
         INPINV = MASKG(INPG)
         IF (INPINV > 0) THEN !inversion variable exists
            DO 2 IVINV=1,NVINV
               IA35 = (INPINV - 1)*NVINV + IVINV
!
!............. work on upper diagonal, copy to lower diagonal
               DO 3 JNPG=INPG,NNPG
                  COV =-1.D0 !error
                  D1 = DSQRT( (XLOCS(INPG) - XLOCS(JNPG))**2 )
                  D2 = DSQRT( (ZLOCS(INPG) - ZLOCS(JNPG))**2 )
                  IF (KERNEL.EQ.2) THEN !exponetial 
                     COV = VELVAR**2/SQRT2PI*DEXP(-D1/SIGMAX)*DEXP(-D2/SIGMAZ)
                  ELSE
                     COV = VELVAR**2*DEXP( -0.5D0*(D1/SIGMAX)**2 ) &
                                    *DEXP( -0.5D0*(D2/SIGMAZ)**2 )  !tarantola pg 112
                  ENDIF
                  !IF (COV < SMALL) COV = 0.D0
                  JNPINV = MASKG(JNPG)
                  IF (JNPINV > 0) THEN !inversion variable exists 
                     DO 4 JVINV=1,NVINV
                        JA35 = (JNPINV - 1)*NVINV + JVINV
                        CALL INFOG2L(IA35,JA35,DESCC, NPROW,NPCOL, MYROW,MYCOL,  &
                                     IIA,JJA, IAROW,IACOL)
                        IF (MYROW == IAROW .AND. MYCOL == IACOL) &
                        COVM(IIA,JJA) = COV
                        CALL INFOG2L(JA35,IA35,DESCC, NPROW,NPCOL, MYROW,MYCOL,  &
                                     IIA,JJA, IAROW,IACOL)
                        IF (MYROW == IAROW .AND. MYCOL == IACOL) &
                        COVM(IIA,JJA) = COV
    4                CONTINUE !loop on inversion variables, j
                  ENDIF !end check on inversion variables 
    3          CONTINUE !loop on anchor nodes, j 
    2       CONTINUE !loop on inversion variables, i 
         ENDIF !end check on inversion variable
    1 CONTINUE !loop on anchor nodes, i
!
!.... matrix multiply A = C_M H 
      IF (MYROW == 0 .AND. MYCOL == 0) &
      WRITE(*,*) 'tarantola_394_tr: Multiplying A = C_M H...'
      CALL PDGEMM('N','N', N,N,N, ONE,COVM,1,1,DESCC, HLOC,1,1,DESCH, &
                  ZERO, A,1,1,DESCA)
      if (ldebug) then
         lwork =-1  !space inquiry 
         allocate(w(n))
         call pdsyev ('N','U',N, A,1,1,DESCA,W,Z,1,1, &
                      descc, work8,lwork, info) !IWORK4,LIWORK, INFO) 
         lwork  = int(work8)
         allocate(work(max(lwork,1)))  
         call pdsyev ('N','U',N, A,1,1,DESCA,W,Z,1,1, &
                      descc, work ,lwork, info)
         if (myrow == 0 .and. mycol == 0) then
            do ia35=1,n
               write(24,*) ia35,w(ia35)
            enddo
            close(24)
         endif
      endif 
    
       
!
!.... matrix vector multiply C_M g 
      IF (MYROW == 0 .AND. MYCOL == 0) &
      WRITE(*,*) 'tarantola_394_tr: Multiplying C_M g...'
      CALL DESCINIT(DESCG,M,1, DESCA(MB_),1, 0,0, ICTXT, MAX(1,MH),INFO)
      ALLOCATE(GLOC(MAX(1,MH)))
      ALLOCATE(YLOC(MAX(1,MH)))
      CALL SET_SOL8(ICTXT,NPROW,NPCOL, M, DESCG, G, GLOC) !g -> local g
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         ALLOCATE(GWORK(N))
      ELSE
         ALLOCATE(GWORK(1))
         GWORK(1) = 0.D0
      ENDIF
      CALL PDGEMV('N',N,N,1.D0,COVM,1,1, DESCC, GLOC,1,1,DESCG,1, 0.D0,YLOC,1,1,DESCG,1) 
      CALL GET_SOL8(ICTXT,NPROW,NPCOL, M, DESCG, YLOC, GWORK) !C_M g = gwork onto (0,0)
      DEALLOCATE(COVM)
      DEALLOCATE(GLOC)
      DEALLOCATE(YLOC)
!
!.... solve the trust region problem
      IF (MYROW == 0 .AND. MYCOL == 0) &
      WRITE(*,*) 'tarantola_394_tr: Beginning trust region problem...'
      CALL NWALG_442(ICTXT,NPROW,NPCOL, K, LDA,DESCA,DELTA,A, GWORK, &
                     RLAM,SEARCH, IERR)
      DEALLOCATE(A)
      RETURN
      END
!                                                                                        !
!========================================================================================! 
!                                                                                        !
      SUBROUTINE NWALG_442(ICTXT,NPROW,NPCOL, K, LDH,DESCH,DELTA,HLOC, G,  &
                           RLAM,SEARCH, IERR) 
!
!     This is the 1D root finding algorithm, equation 4.42, in Nocedal and Wright 
      USE EQN439_DAT 
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: HLOC(LDH,*), G(*), DELTA
      INTEGER*4, INTENT(IN) :: DESCH(9), ICTXT, NPROW, NPCOL, LDH, K 
      REAL*8, INTENT(OUT) :: SEARCH(*), RLAM 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: A(:,:), WORK(:)
      INTEGER*4, ALLOCATABLE :: IWORK(:) 
      CHARACTER(80) FLNAME
      CHARACTER(5) CK 
      REAL*8 WORK8 
      INTEGER*4, PARAMETER :: N1=1  !optimization in lambda which is one variable
      INTEGER*4, PARAMETER :: LWA=8 !hybrdy workspace size (n*(3*n+13))/2; n = 1
      REAL*8 X(N1), FVEC(N1), WA(LWA), TOL
      INTEGER*4 LDA, LIWORK, LWORK, IWORK4, I, IIA,JJA, IAROW,IACOL
!.... BLACS stuff
      INTEGER*4 DESCA(9)   !A matrix descriptor
      INTEGER*4 NUMROC     !calculates NUmber of Rows Or Columns
      INTEGER*4 INFO       !error flag from ScaLapack 
      INTEGER*4 M          !number of rows and columns of Jacobian
      INTEGER*4 MH         !number of local rows in h
      INTEGER*4 NH         !number of local columns in h
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &    
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )  
      EXTERNAL EQN439
!
!----------------------------------------------------------------------------------------!
!
!.... recall our process grid and sizes 
      IERR = 0
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
      M = DESCH(M_)
      N = DESCH(N_) 
      IF (M /= N) THEN
         WRITE(*,*) 'nwalg_442: Error m /= n!'
         IERR = 1
         RETURN
      ENDIF
      LDA = DESCH(LLD_) 
      MH = NUMROC(M, DESCH(MB_),MYROW, 0,NPROW) !number of rows in local H matrix
      NH = NUMROC(N, DESCH(NB_),MYCOL, 0,NPCOL) !number of columns in local H matrix
      CALL ICOPY(9,DESCH,1,DESCA,1)
      CALL ICOPY(9,DESCH,1,DESCZ,1)
!
!.... set matrix to factor 
      ALLOCATE(A(LDA,MAX(NH,1)),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'nwalg_442: Error setting space for A on grid',MYROW,MYCOL
         RETURN
      ENDIF 
      CALL COPY_HMAT(ICTXT,NPROW,NPCOL, LDH,N,DESCH, HLOC,A)
!
!.... set space for eigenvalue decomposition
      IWORK4 = 0
      ALLOCATE(W(N)) 
      ALLOCATE(Z(LDA,MAX(NH,1)),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'nwalg_442: Error setting space for Z on grid',MYROW,MYCOL
         RETURN
      ENDIF
      LWORK =-1  !space inquiry 
      LIWORK =-1 !space inquiry
      CALL PDSYEV ('V','U',N, A,1,1,DESCA,W,Z,1,1, &
                   DESCZ, WORK8,LWORK, INFO) !IWORK4,LIWORK, INFO) 
      IF (INFO /= 0) THEN
         WRITE(*,*) 'nwalg_442: Error in space estimate:',INFO, MYROW,MYCOL
         IERR = 1
         RETURN
      ENDIF
      LWORK  = INT(WORK8)
      ALLOCATE(WORK(MAX(LWORK,1)),STAT=IERR)  
      IF (IERR /= 0) THEN
         WRITE(*,*) 'nwalg_442: Error setting space for work', MYROW,MYCOL
         RETURN
      ENDIF
      !LIWORK = IWORK4
      !ALLOCATE(IWORK(MAX(LIWORK,1)),STAT=IERR)
      !IF (IERR /= 0) THEN
      !   WRITE(*,*) 'nwalg_442: Error setting space for iwork',MYROW,MYCOL
      !   RETURN
      !ENDIF
!
!.... calculate the eigenvalues and eigenvectors
      IF (MYROW == 0 .AND. MYCOL == 0) &
      WRITE(*,*) 'nwalg_442: Calculating eigendecomposition...'
      CALL PDSYEV ('V','U',N, A,1,1,DESCA,W,Z,1,1, &
                   DESCZ, WORK,LWORK, INFO) !IWORK,LIWORK, INFO) 
      IF (INFO /= 0) THEN
         WRITE(*,*) 'nwalg_442: Error in pdyevd:',INFO, MYROW,MYCOL
         IERR = 1
         RETURN
      ENDIF
      DEALLOCATE(WORK)
      IF (ALLOCATED(IWORK)) DEALLOCATE(IWORK)
      DEALLOCATE(A)
!
!.... save the eigenspectrum, b/c it is informative
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         FLNAME(1:80) = ' '
         CK(1:5) = ' '
         WRITE(CK,'(I5)') K
         CK = ADJUSTL(CK)
         FLNAME = 'eigen_spectrum-'//TRIM(CK)//'.txt'
         FLNAME = ADJUSTL(FLNAME)
         OPEN(UNIT=45,FILE=TRIM(FLNAME))
         DO 1 I=1,N 
            WRITE(45,*) I, W(I)
    1    CONTINUE  
         CLOSE(45) 
         WRITE(*,*) 'nwalg_442: Setting the RHS...'
      ENDIF
!
!.... initialize gradient in root finding
      ICTXT_PASS = ICTXT
      NPROW_PASS = NPROW
      NPCOL_PASS = NPCOL 
      CALL DESCINIT(DESCG,M,1, DESCH(MB_),1, 0,0, ICTXT, MAX(1,MH),INFO)
      ALLOCATE(GLOC(MAX(1,MH)))
      ALLOCATE(YLOC(MAX(1,MH)))
      ALLOCATE(PLOC(MAX(1,MH)))
      GLOC(:) = 0.D0
      YLOC(:) = 0.D0 
      PLOC(:) = 0.D0 
      CALL SET_SOL8(ICTXT,NPROW,NPCOL, M, DESCG, G, GLOC)
!
!.... begin the root finding
      LDONE = .FALSE.
      LVERB = .TRUE.   !track optimization
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         INFO = 1
         DELTA_PASS = DELTA
         print *, delta_pass
         X(1) = SUM(W)/DFLOAT(N) !guess lambda as the average of the spectrum
        
         WRITE(*,*) 'nwalg_442: Beginning lambda optimization with rlam:',X(1) 
         TOL = 1.D-8 !we dont need to do a very good job
         FVEC(1) = 0.D0 
         CALL HYBRD1(EQN439,N1,X,FVEC,TOL, INFO,WA,LWA)   
         LDONE = .TRUE.
         RLAM = X(1) 
         IF (INFO == 0) THEN
            WRITE(*,*) 'nwalg_442: Error Improper inputs for hybrd1'
            IERR = 1
         ELSEIF (INFO == 2) THEN
            WRITE(*,*) 'nwalg_442: Error number of calls to fcn reached or exceed limit'
            IERR = 1
         ELSEIF (INFO == 3) THEN
            WRITE(*,*) 'nwalg_442: Warning tol is too small to improve solution'
         ELSEIF (INFO == 4) THEN 
            WRITE(*,*) 'nwalg_442: Error iteration is not making good progress' 
            IERR = 1
         ELSE
            WRITE(*,*) 'nwalg_442: Converged to solution:',RLAM
            WRITE(*,*) '           With error:',FVEC(1) 
         ENDIF
         IF (IERR /= 0) THEN
            WRITE(*,*) 'nwalg_442: Continuing with lambda =',RLAM
         ENDIF
         CALL EQN439(N1,X,FVEC,IERR) !tell remaining processes we are done 
      ELSE !just go to the subroutine
         CALL EQN439(N1,X,FVEC,IERR) 
      ENDIF
!
!.... report errors
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         CALL IGEBS2D(ICTXT,'A',' ',1,1,IERR,1)
      ELSE
         CALL IGEBR2D(ICTXT,'A',' ',1,1,IERR,1,0,0)
      ENDIF
!
!.... invert in 3 steps: 1) y = Q^T g; 2) y <- y/(Gamma + rlam I); g =-Qy 
      IF (IERR == 0) THEN
         CALL PDGEMV('T',N,N, 1.D0,Z,1,1, DESCZ, GLOC,1,1,DESCG,1, 0.D0,YLOC,1,1,DESCG,1)
         DO 2 I=1,N 
            CALL INFOG2L(I,1,DESCG, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                         IAROW,IACOL)
            IF (MYROW == IAROW .AND. MYCOL == IACOL) YLOC(IIA) = YLOC(IIA)/(W(I) + RLAM)
    2    CONTINUE 
         CALL PDGEMV('N',N,N,-1.D0,Z,1,1, DESCZ, YLOC,1,1,DESCG,1, 0.D0,GLOC,1,1,DESCG,1)
         CALL GET_SOL8(ICTXT,NPROW,NPCOL, N, DESCG, GLOC, SEARCH) !send p back to (0,0)
      ENDIF
!
!.... clean
      DEALLOCATE(GLOC)
      DEALLOCATE(YLOC)
      DEALLOCATE(PLOC) 
      DEALLOCATE(Z)
      DEALLOCATE(W)  
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE EQN439(N1,X,FVEC,IFLAG) 
!
!     Applies equation 4.39 of Nocedal and Wright
!
!     INPUT      MEANING
!     -----      ------- 
!     N1         length of vectors; 1 
!     X          current estimate of root
!
!     OUTPUT     MEANING
!     ------     ------- 
!     FVEC       result 
!     IFLAG      error flag
!
!.... variable declarations
      USE EQN439_DAT
      IMPLICIT NONE 
      REAL*8, INTENT(IN) :: X(N1)
      INTEGER*4, INTENT(IN) :: N1 
      REAL*8, INTENT(OUT) :: FVEC(N1) 
      INTEGER*4, INTENT(OUT) :: IFLAG 
!.... local variables
      REAL*8 RLAM, DELTA, XSUM
      INTEGER*4 NPROW, NPCOL, ICTXT, IAROW, IACOL, IIA, JJA, I, JOB  
!
!----------------------------------------------------------------------------------------!
!
!.... intialize
      IFLAG = 0 
      ICTXT = ICTXT_PASS
      NPROW = NPROW_PASS
      NPCOL = NPCOL_PASS 
      IF (.NOT.LDONE) THEN
         JOB = 1
      ELSE
         JOB = 0
      ENDIF
      IF (MYROW == 0 .AND. MYCOL == 0) RLAM = X(1)
  100 CONTINUE !return
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         DELTA = DELTA_PASS 
         CALL IGEBS2D(ICTXT,'A',' ',1,1, JOB,  1)
         CALL DGEBS2D(ICTXT,'A',' ',1,1, RLAM, 1)
         CALL DGEBS2D(ICTXT,'A',' ',1,1, DELTA,1)
      ELSE !block until the (0,0) location tells us what to do 
         CALL IGEBR2D(ICTXT,'A',' ',1,1, JOB,  1,0,0)
         CALL DGEBR2D(ICTXT,'A',' ',1,1, RLAM, 1,0,0)
         CALL DGEBR2D(ICTXT,'A',' ',1,1, DELTA,1,0,0)
      ENDIF
      IF (JOB == 0) RETURN !we are done
!
!.... this is the only potential problem
      IF (RLAM < 0.D0) THEN
         WRITE(*,*) 'eqn439: Error rlam is negative1',RLAM
         IFLAG = 1
         RETURN
      ENDIF
!
!.... perform transpose matrix vector multiply this will hold q_j g
      CALL PDGEMV('T',N,N,1.D0,Z,1,1, DESCZ, GLOC,1,1,DESCG,1, 0.D0,YLOC,1,1,DESCG,1)
      DO 1 I=1,N
         CALL INFOG2L(I,1,DESCG, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                      IAROW,IACOL)
         IF (MYROW == IAROW .AND. MYCOL == IACOL) YLOC(IIA) = YLOC(IIA)/(W(I) + RLAM)
    1 CONTINUE 
      CALL PDGEMV('N',N,N,1.D0,Z,1,1, DESCZ, YLOC,1,1,DESCG,1, 0.D0,PLOC,1,1,DESCG,1)
      XSUM = 0.D0
      CALL PDNRM2(N, XSUM,PLOC,1,1,DESCG,1) !calculate 2 norm of soln vector 
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         CALL DGEBS2D(ICTXT,'A',' ',1,1, XSUM,1)
      ELSE
         CALL DGEBR2D(ICTXT,'A',' ',1,1, XSUM,1,0,0)
      ENDIF
      FVEC(1) = XSUM - DELTA 
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         IF (LVERB) WRITE(*,*) 'eqn439: rlam,|p|_2,error:', &
         SNGL(RLAM),SNGL(XSUM),SNGL(FVEC(1))
         RETURN
      ELSE
         GOTO 100 
      ENDIF
      RETURN
      END  
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE NWALG41(ICTXT,NPROW,NPCOL, LDH,LDB,NORM,DESCH,DESCB,  &
                         DELTA,GRAD,DIAG,HLOC, B, SEARCH, IERR) 
!
!     Implements Agorithm 4.1 of Nocedal and Wright for normalized solution of 
!      
!
      IMPLICIT NONE
      INCLUDE 'mpif.h'

      REAL*8, INTENT(IN) :: HLOC(LDH,*), GRAD(*), DIAG(*), DELTA 
      REAL*8, INTENT(INOUT) :: B(LDB) 
      INTEGER*4, INTENT(IN) :: DESCH(9), DESCB(9), ICTXT, NPROW, NPCOL, LDH, LDB,  &
                               NORM 
      REAL*8, INTENT(OUT) :: SEARCH(*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: HFAC(:,:)  !Upper triangular matrix to factor
      REAL*8, ALLOCATABLE :: PGLOB(:)   !solution of R^T R p =-g
      REAL*8, ALLOCATABLE :: QGLOB(:)   !solution of R^T q = p
      REAL*8, ALLOCATABLE :: DIAG2(:)   !diag R^T R squared
      REAL*8 RLAM                       !damping factor to solve for
      REAL*8 PMAG, QMAG, DNRM2, RLAM0  
      INTEGER*4 K
!.... BLACS stuff
      INTEGER*4 MYROW      !processes' row in BLACS grid
      INTEGER*4 MYCOL      !processes' column in BLACS grid        
      INTEGER*4 NUMROC     !calculates NUmber of Rows Or Columns
      INTEGER*4 INFO       !error flag from ScaLapack 
      INTEGER*4 M, N       !number of rows and columns of Jacobian
      INTEGER*4 NH         !number of local columsn in h
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 ) 
!
!----------------------------------------------------------------------------------------!
!
!.... recall our process grid 
      IERR = 0
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!
!.... get the diagonal and broadcast it to every other process in group
      M = DESCH(M_)
      N = DESCH(N_) 
      NH = NUMROC(N, DESCH(NB_),MYCOL, 0,NPCOL) !number of columns in local H matrix
      print *, m,n,nh
      ALLOCATE(HFAC(LDH,MAX(NH,1)))
!
!.... initially just set rlam to the average of the diagonal 
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         CALL DCOPY(N,GRAD,1,SEARCH,1)
         RLAM = SUM(DIAG(1:N))/DFLOAT(N)  
         ALLOCATE(PGLOB(N))
         ALLOCATE(QGLOB(N))
      ELSE
         ALLOCATE(PGLOB(1)) 
         ALLOCATE(QGLOB(1))
      ENDIF
      ALLOCATE(DIAG2(N))
      DO 1 K=1,N
         DIAG2(K) = DIAG(K)!**2 
    1 CONTINUE
!
!.... this is the iterative solution with safeguard 
      RLAM0 = 0.D0
      DO 2000 K=1,50 !MAXIT
!
!....... broadcast lambda
         IF (MYROW == 0 .AND. MYCOL == 0) THEN 
            WRITE(*,*) 'nwalg41: Initializing iteration:',K,RLAM
            RLAM0 = RLAM
            CALL DGEBS2D(ICTXT,'A',' ',1,1,RLAM,1)
         ELSE
            CALL DGEBR2D(ICTXT,'A',' ',1,1,RLAM,1, 0,0)
         ENDIF
         CALL COPY_HMAT(ICTXT,NPROW,NPCOL, LDH,N,DESCH, HLOC,HFAC)
         CALL ADD_DIAG(ICTXT,NPROW,NPCOL, LDH,N, DESCH, RLAM,DIAG2, HFAC)
         IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'nwalg41: Factoring...'
         CALL PDPOTRF('L',N, HFAC,1,1, DESCH, INFO)
         IF (INFO /= 0 .AND. RLAM > 0.D0) THEN
            WRITE(*,*) 'nwalg41: Error the matrix is not SPD!',INFO,RLAM
            GOTO 2500 
         ENDIF  
         CALL SET_SOL8(ICTXT,NPROW,NPCOL, M, DESCB, GRAD, B) !set B to -g
         IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'nwalg41: Solving R^T R p =-g'
         CALL PDPOTRS('L',N, 1, HFAC,1,1, DESCH, B, 1,1, DESCB,INFO) !solve R^T R p =-g
         IF (INFO /= 0) THEN
            WRITE(*,*) 'nwalg41: There was an error in the solve phase!'
            IERR = 1
            GOTO 2500
         ENDIF
         CALL GET_SOL8(ICTXT,NPROW,NPCOL, N, DESCB, B, PGLOB) !send p back to (0,0)
         IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'nwalg41: Solving R^T q = p'
         CALL PDTRSM('L', 'L', 'N', 'N', N, 1, &
                     1.D0,HLOC, 1,1, DESCH, B, 1, 1, DESCB) !solve R^T q = p
         CALL GET_SOL8(ICTXT,NPROW,NPCOL, N, DESCB, B, QGLOB) !send q back to (0,0)
         IF (MYROW == 0 .AND. MYCOL == 0) THEN
            IF (NORM == 1) THEN
               PMAG = SUM(DABS(PGLOB))
               QMAG = SUM(DABS(QGLOB)) 
            ELSEIF (NORM == 2) THEN
               PMAG = DNRM2(N,PGLOB,1)  
               QMAG = DNRM2(N,QGLOB,1)
            ELSE !infinity norm
               PMAG = MAXVAL(DABS(PGLOB))
               QMAG = MAXVAL(DABS(QGLOB))    
            ENDIF
            RLAM = RLAM + (PMAG/QMAG)**2*(PMAG - DELTA)/DELTA 
            print *, 'dif',pmag,qmag,rlam - rlam0, rlam
            IF (RLAM <= 0.D0) RLAM = RLAM0/2.D0
  
         ENDIF
 2000 CONTINUE
      IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'nwalg41: Warning hit iteration limit!'
      IERR =-1
 2500 CONTINUE !break ahead for errors
      IF (MYROW == 0 .AND. MYCOL == 0) CALL DCOPY(N,PGLOB,1,SEARCH,1) 
!
!.... free the process grid
      DEALLOCATE(HFAC) 
      DEALLOCATE(PGLOB)
      DEALLOCATE(QGLOB)
      DEALLOCATE(DIAG2) 
      RETURN 
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE TARANTOLA_394(ICTXT, LDH,LDC,NPROW,NPCOL, DESCH,DESCC, GRAD,XMOD, &
                               HLOC,COVM, P,IERR ) 
!
!     Solves the quasi-Newton algorithm (3.93) of Tarantola with Hloc the distributed 
!     matrix Re{ adj(J) J } and Covm the inverse of the covariance matrix
!
!.... varaible declarations
      implicit none
      REAL*8, INTENT(IN) :: HLOC(LDH,*),  COVM(LDC,*), GRAD(*), XMOD(*) 
      INTEGER*4, INTENT(IN) :: DESCH(9), DESCC(9), ICTXT,LDH,LDC, NPROW,NPCOL  
      REAL*8, INTENT(OUT) :: P(*)
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      REAL*8, ALLOCATABLE :: ALOC(:,:), BLOC(:), XLOC(:)
      !INTEGER*4, ALLOCATABLE :: IPIV(:) 
      INTEGER*4 DESCA(9), DESCB(9), NUMROC, INFO, I,J,IIA,JJA,IAROW,IACOL,MYROW,MYCOL,  &
                M,N,MB,NB,MH,NH,IRSRC,ICSRC!, NWORK !calculate NUmber of Rows Or Columns
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
!
!----------------------------------------------------------------------------------------!
!
!.... initial check
      IERR = 0
      DO 1 I=1,9
         IF (DESCC(I) /= DESCH(I)) THEN
            WRITE(*,*) 'tarantola_394: Error cant add matrices!'
            IERR = 1
         ENDIF
         DESCA(I) = DESCH(I)
    1 CONTINUE 
!
!.... recall the BLACS grid
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL) 
      M  = DESCH(M_) 
      N  = DESCH(N_) 
      IF (M /= N) THEN
         WRITE(*,*) 'tarantola_394: Matrices are not square!'
         IERR = 1
         RETURN
      ENDIF
      IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'tarantola_394: Setting matrix...'
      MB = DESCH(MB_)
      NB = DESCH(NB_)
      IRSRC = DESCH(RSRC_)
      ICSRC = DESCH(CSRC_)
      MH = NUMROC(M,MB,MYROW,IRSRC,NPROW) 
      NH = NUMROC(N,NB,MYCOL,ICSRC,NPCOL) 
      CALL DESCINIT(DESCB,N,1, MB, 1, IRSRC,ICSRC, ICTXT, MAX(1,MH),INFO)  
      ALLOCATE(ALOC(LDH,NH))
      ALLOCATE(BLOC(DESCB(LLD_)))
      ALLOCATE(XLOC(DESCB(LLD_)))
      BLOC(:) = 0.D0 
      ALOC(:,:) = 0.D0 
      DO 2 I=1,N
         DO 3 J=1,N 
            CALL INFOG2L(I,J,DESCA, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                         IAROW,IACOL)
            IF (MYROW == IAROW .AND. MYCOL == IACOL) & 
            ALOC(IIA,JJA) = HLOC(IIA,JJA) + COVM(IIA,JJA)
    3    CONTINUE
    2 CONTINUE 
      print *, minval(aloc),maxval(aloc)
      IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'tarantola_394: Setting RHS...'
      CALL SET_SOL8(ICTXT,NPROW,NPCOL, N, DESCB, XMOD,XLOC) !x <- (m - m_{prior}) 
      CALL SET_SOL8(ICTXT,NPROW,NPCOL, N, DESCB, GRAD,BLOC) !local gradient
      !(adj(J)J + inv(CM)dp =-grad - inv(Cm)*(m - m_{prior}) =-Ax - b for PBLAS
      CALL PDGEMV('N',N,N,-1.D0,COVM,1,1,DESCC, XLOC,1,1,DESCB,1,-1.D0,BLOC,1,1,DESCB,1) 
!
!.... cholesky factorization
      CALL PDPOTRF('U',N,ALOC,1,1,DESCA,INFO)
      IF (INFO /= 0) THEN
         WRITE(*,*) 'tarantola_394: Error in pdpotrf!',INFO
         IERR = 1
         RETURN
      ENDIF
      CALL PDPOTRS('U',N,1, ALOC,1,1,DESCA, BLOC,1,1, DESCB, INFO)
      IF (INFO /= 0) THEN
         WRITE(*,*) 'tarantola_394: Error in pdpotrs!',INFO
         IERR = 1
         RETURN
      ENDIF
!
!.... LU factorization 
!     IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'tarantola_394: Factoring matrix...'
!     NWORK = MAX( NUMROC(M,MB,MYROW,IRSRC,NPROW) + MB, 1 )
!     ALLOCATE(IPIV(NWORK))
!     CALL PDGETRF(N,N, ALOC,1,1, DESCA, IPIV,INFO)
!     IF (INFO /= 0) THEN
!        WRITE(*,*) 'tarantola_394: Error getting LU factor!'
!        IERR = 1
!        RETURN
!     ENDIF
!
!.... solve Ax = b  and fetch solution onto (0,0) 
!     IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'tarantola_394: Solution phase...'
!     CALL PDGETRS('N',N,1, ALOC,1,1,DESCA, IPIV,BLOC, 1,1, DESCB, INFO)  
!     IF (INFO /= 0) THEN
!        WRITE(*,*) 'tarantola_394: Error calling pdgetrs!'
!        IERR = 1
!        RETURN
!     ENDIF
      CALL GET_SOL8(ICTXT,NPROW,NPCOL, N, DESCB, BLOC, P)
!
!.... clean space
      IF (ALLOCATED(ALOC)) DEALLOCATE(ALOC)
      IF (ALLOCATED(BLOC)) DEALLOCATE(BLOC)
      IF (ALLOCATED(XLOC)) DEALLOCATE(XLOC) 
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE COVMODEL_DIST(LDC, NPROW,NPCOL,ICTXT,                  &
                               KERNEL,SIGMAX,SIGMAZ,VPVAR, DESCC,       &
                               INV,MSH, COVM,IERR) 
!
!     Calculates the inverse of the distributed model covariance matrix.  This 
!     routine inverts covm via Cholesky which will check covm is positive
!     definite.  COVM must be positive definite so that 
! 
!      [ J_1^+ J_2^+ ... J_m^+ C_m^{-1/2}][     J_1    ] dm =-g 
!                                         [     J_2    ] 
!                                         [      .     ]
!                                         [      .     ]
!                                         [      .     ]
!                                         [     J_m    ] 
!                                         [ C_m^{-1/2} ]
!     makes sense 
!     
!
!     INPUT      MEANING
!     -----      ------- 
!     DESCC      BLACS descriptor for covariance matrix
!     ICTXT      BLACS context
!     INV        inversion variables
!     KERNEL     1 -> Gaussian (default), 2 -> exponential
!     LDC        leading dimension for C
!     MSH        mesh variables
!     NPCOL      number of process columns
!     NPROW      number of process rows 
!     SIGMAX     standard deviation distance in x
!     SIGMAZ     standard deviation distance in z 
!     VPVAR      Vp velocity variance (m/s)
!
!     OUTPUT     MEANING
!     ------     -------
!     COVM       covariance matrix 
!     IERR       error flag; will occur if Covm is not positive definite
!
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO)  INV
      TYPE (MESH_INFO) MSH
      REAL*8, INTENT(IN) :: SIGMAX, SIGMAZ, VPVAR
      INTEGER*4, INTENT(IN) :: DESCC(9), LDC, NPROW, NPCOL, KERNEL, ICTXT
      REAL*8, INTENT(OUT) :: COVM(LDC,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      !REAL*8, ALLOCATABLE :: WORK(:)
      !INTEGER*4, ALLOCATABLE :: IPIV(:), IWORK(:) 
      REAL*8 D1, D2, COV, VARI2, VARI !, WORK8
      INTEGER*4 INPINV,IVINV, JNPINV,JVINV, IA35, JA35,  &
                IIA, JJA, IAROW, IACOL, MYROW, MYCOL, INFO, N, M, MB, IRSRC, INPG, &
                JNPG, JAROW, JACOL, IIAT, JJAT !,LWORK,LIWORK,IWORK4,NWORK
      LOGICAL*4 LCOL
      REAL*8 SQRT2PI, SMALL 
      PARAMETER(SQRT2PI = 2.5066282746310002D0) !sqrt(2*pi)
      !INTEGER*4 NUMROC !calculate NUmber of Rows Or Columns
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL*8, PARAMETER :: ONE = 1.D0
      !real*8, allocatable :: test(:,:), w(:)
!
!----------------------------------------------------------------------------------------!
!
!.... intial check
      IERR = 0
      IF (inv%NVINV > 1) & 
      WRITE(*,*) 'covmodel_dist: Warning nvinv > 1, this isnt tested!'
      if (inv%nvinv > 1) print *, 'covmodel_dist: need new covariances!'
!
!.... recall the BLACS grid
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
      SMALL = EPSILON(1.D0) 
      IF (MYROW == 0 .AND. MYCOL == 0) & 
      WRITE(*,*) 'covmodel_dist: Setting model covariance matrix...'
      !if (myrow == 0 .and. mycol == 0) allocate(test(inv%na35,inv%na35))
      LCOL = .TRUE. !i want to use the cholesky decomposition
  100 CONTINUE
!
!.... set the model covariance
      DO 1 INPG=1,msh%NNPG 
         INPINV = inv%MASKG(INPG) 
         IF (INPINV > 0) THEN !inversion variable exists
            DO 2 IVINV=1,inv%NVINV 
               IA35 = (INPINV - 1)*inv%NVINV + IVINV
!
!............. work on upper diagonal
               DO 3 JNPG=INPG,msh%NNPG
!
!................ master defines distance and covariance
                  COV =-1.D0 !error
                  IF (MYROW == 0 .AND. MYCOL == 0 .OR. LCOL) THEN
                     D1 = DSQRT( (msh%XLOCS(INPG) - msh%XLOCS(JNPG))**2 )
                     D2 = DSQRT( (msh%ZLOCS(INPG) - msh%ZLOCS(JNPG))**2 )
                     IF (KERNEL.EQ.2) THEN !exponetial 
                        COV = 1.D0**2/SQRT2PI*DEXP(-D1/SIGMAX)*DEXP(-D2/SIGMAZ)
                     ELSE 
                        COV = 1.D0**2*DEXP( -0.5D0*(D1/SIGMAX)**2 ) &
                                     *DEXP( -0.5D0*(D2/SIGMAZ)**2 )  !tarantola pg 112
                     ENDIF
                     !IF (COV < SMALL) COV = 0.D0 
                  ENDIF
                  JNPINV = inv%MASKG(JNPG)
                  IF (JNPINV > 0) THEN !inversion variable exists 
                     DO 4 JVINV=1,inv%NVINV
!
!...................... upper half
                        JA35 = (JNPINV - 1)*inv%NVINV + JVINV
                        CALL INFOG2L(IA35,JA35,DESCC, NPROW,NPCOL, MYROW,MYCOL, &
                                     IIA,JJA, IAROW,IACOL) 
                        IF (MYROW == IAROW .AND. MYCOL == IACOL) &
                        COVM(IIA,JJA) = COV*inv%WMASK(IA35)*inv%WMASK(JA35)
!
!...................... lower half
                        CALL INFOG2L(JA35,IA35,DESCC, NPROW,NPCOL, MYROW,MYCOL,  &
                                     IIA,JJA, IAROW,IACOL)
                        IF (MYROW == IAROW .AND. MYCOL == IACOL) &
                        COVM(IIA,JJA) = COV*inv%WMASK(IA35)*inv%WMASK(JA35)
    4                CONTINUE !loop on inversion variables
                  ENDIF !end check if inversion variable exists
    3          CONTINUE !loop on anchor nodes 2
    2       CONTINUE ! loop on inversion variables 
         ENDIF !end check if inversion variable exists
    1 CONTINUE !Loop on anchor nodes 1
!
!.... invert the positive definite matrix 
      IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'covmodel_dist: Factoring matrix...'
      M = DESCC(M_)
      N = DESCC(N_) 
      MB = DESCC(MB_) 
      IRSRC = DESCC(RSRC_)
      IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'covmodel_dist: Factoring...'
      CALL PDPOTRF('U',N,COVM,1,1,DESCC, INFO) 
      IF (INFO /= 0) THEN
         IF (MYROW == 0 .AND. MYCOL == 0) &
         WRITE(*,*) 'covmodel_dist: pdpotrf error!',INFO
         IERR = 1
         RETURN
      ENDIF
      IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'covmodel_dist: Inverting...'
      CALL PDPOTRI('U',N,COVM,1,1,DESCC,INFO)
      IF (INFO /= 0) THEN
         IF (MYROW == 0 .AND. MYCOL == 0) &
         WRITE(*,*) 'covmodel_dist: pdpotri error!',INFO
         IERR = 1
         RETURN
      ENDIF
!
!.... communicate inverse to lower half of matrix
      IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'covmodel_dist filling lower half...'
      DO 11 IA35=1,inv%NA35
         DO 12 JA35=IA35+1,inv%NA35
            CALL INFOG2L(IA35,JA35,DESCC, NPROW,NPCOL, MYROW,MYCOL, &
                         IIA ,JJA,  IAROW,IACOL)
            CALL INFOG2L(JA35,IA35,DESCC, NPROW,NPCOL, MYROW,MYCOL, &
                         IIAT,JJAT, JAROW,JACOL)
            IF (MYROW == IAROW .AND. MYCOL == IACOL) THEN
               IF (MYROW == JAROW .AND. MYCOL == JACOL) THEN
                  COVM(IIAT,JJAT) = COVM(IIA,JJA)
               ELSE !send to (iarow,iacol)
                  CALL DGESD2D(ICTXT,1,1,COVM(IIA,JJA),1, JAROW,JACOL)
               ENDIF
            ELSE
               IF (MYROW == JAROW .AND. MYCOL == JACOL) &
               CALL DGERV2D(ICTXT,1,1,COVM(IIAT,JJAT),1, IAROW,IACOL)
            ENDIF
   12    CONTINUE
   11 CONTINUE 
!
!.... include the variance or scale factor
      VARI2 = (1.D0/VPVAR)**2
      VARI  = 1.D0/VPVAR
      DO 13 IA35=1,inv%NA35
         DO 14 JA35=1,inv%NA35
            CALL INFOG2L(IA35,JA35,DESCC, NPROW,NPCOL, MYROW,MYCOL, &
                         IIA ,JJA, IAROW,IACOL)
            IF (MYROW == IAROW .AND. MYCOL == IACOL) COVM(IIA,JJA) = VARI2*COVM(IIA,JJA)
   14    CONTINUE
   13 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE DIST_LRANK(LDH, NPROW,NPCOL,ICTXT, LSURF,K, THRESH, DESCH,  &
                            HMAT, IERR) 
!
!     Solves Re{adj(J) J} dm =-g via the truncated singular value decomposition.  
!     This may be useful for surface waves.  By theorem 3.3 of Demmel if A is 
!     symmetric with svd A = U S V^T then it also has an eigendecomposition 
!     A = U L U^T.  The signgular values are abs(L) and the right singular
!     vectors are sign(L_i) u_i.  Of course, we ignore any singular vectors 
!     with singular values less than 0 so the last point is moot. 

      IMPLICIT NONE 
      REAL*8, INTENT(INOUT) :: HMAT(LDH,*)
      REAL*8, INTENT(IN) :: THRESH
      INTEGER*4, INTENT(IN) :: DESCH(9), LDH, NPROW, NPCOL, ICTXT, K
      LOGICAL*4, INTENT(IN) :: LSURF
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      REAL*8, ALLOCATABLE :: A(:,:), Z(:,:), W(:), S(:), WORK(:)
      CHARACTER(80) FLNAME
      CHARACTER(5) CK 
      REAL*8 WORK8, CUT, XFACT, ONE, ZERO
      INTEGER*4 DESCA(9), DESCZ(9), LWORK, MH, NH, MA, NA, MB, NB, &
                I, J, I1, J1, IAROW, IACOL, JAROW, JACOL, IIA, JJA, IIAD, JJAD, LDA,   &
                MYROW, MYCOL, INFO, N, M, NKEEP, MZ, NZ, LDZ
      INTEGER*4 NUMROC !calculate NUmber of Rows Or Columns
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &    
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      PARAMETER(ONE = 1.D0, ZERO = 0.D0)
!
!----------------------------------------------------------------------------------------!
!
!.... figure out cutoff
      IERR = 0
      CUT = THRESH
      IF (THRESH.LE.0.D0) CUT = EPSILON(1.D0)
!
!.... set space and copy upper triangle
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
      M = DESCH(M_)
      N = DESCH(N_)
      LDA = DESCH(LLD_)
      MH = NUMROC(M, DESCH(MB_),MYROW, 0,NPROW) !number of rows in local H matrix
      NH = NUMROC(N, DESCH(NB_),MYCOL, 0,NPCOL) !number of columns in local H matrix
      ALLOCATE(A(LDA,MAX(1,NH)),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'dist_lrank: Error setting space for A'
         RETURN
      ENDIF
      ALLOCATE(Z(LDA,MAX(1,NH)),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'dist_lrank: Error setting space for Z'
         RETURN
      ENDIF
      ALLOCATE(W(N))
      A(:,:) = 0.D0
      DESCA(1:9) = DESCH(1:9)
      DESCZ(1:9) = DESCH(1:9)
      DO 1 I=1,N
         DO 2 J=I,N
            CALL INFOG2L(I,J ,DESCZ, NPROW,NPCOL, MYROW,MYCOL, IIA, JJA,  &
                         IAROW,IACOL)
            IF (MYROW == IAROW .AND. MYCOL == IACOL) &
            A(IIA,JJA) = HMAT(IIA,JJA)
    2    CONTINUE
    1 CONTINUE
!
!.... calculate the eigenvalue decomposition
      IF (MYROW == 0 .AND. MYCOL == 0) &
      WRITE(*,*) 'dist_lrank: Calculating eigendecomposition...'
      LWORK =-1
      CALL PDSYEV('V','U',N,A,1,1,DESCA, W,  &
                  Z,1,1, DESCZ, WORK8,LWORK, INFO)
      LWORK = INT(WORK8)
      ALLOCATE(WORK(LWORK),STAT=IERR)
      IF (IERR /= 0) THEN
         IF (MYROW == 0 .AND. MYCOL == 0) &
         WRITE(*,*) 'dist_lrank: Error allocating space for work',IERR
         RETURN
      ENDIF
      CALL PDSYEV('V','U',N,A,1,1,DESCA, W,  &
                  Z,1,1, DESCZ, WORK ,LWORK, INFO)
      IF (INFO /= 0) THEN
         IF (MYROW == 0 .AND. MYCOL == 0) &
         WRITE(*,*) 'dist_lrank: Error calling pdysev',INFO
         IERR = 1
         RETURN
      ENDIF
      DEALLOCATE(WORK) !done with workspace 
!
!.... count the singular values and reset A
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         FLNAME(1:80) = ' '
         CK(1:5) = ' '
         WRITE(CK,'(I5)') K
         CK = ADJUSTL(CK)
         IF (LSURF) THEN
            FLNAME = 'eigen_spc_srf-'//TRIM(CK)//'.txt'
         ELSE
            FLNAME = 'eigen_spc_bdy-'//TRIM(CK)//'.txt'
         ENDIF
         FLNAME = ADJUSTL(FLNAME)
         OPEN(UNIT=64,FILE=TRIM(FLNAME))
      ENDIF
      ALLOCATE(S(N))
      S(1:N) = 0.D0
      I1 = 0
      DO 21 I=N,1,-1
         IF (MYROW == 0 .AND. MYCOL == 0) WRITE(64,*) W(I)
         IF (W(I).GT.CUT) THEN
            I1 = I1 + 1
            S(I1) = DABS(W(I))
         ENDIF
   21 CONTINUE
      IF (MYROW == 0 .AND. MYCOL == 0) CLOSE(64)
      NKEEP = I1
!
!.... reset vectors
      DEALLOCATE(A)
      MB = DESCH(MB_)
      NB = MIN(DESCH(NB_),NKEEP) !may be way less columsn now
      MA = NUMROC(N,    MB,MYROW,0,NPROW) !unchange; local rows in A
      NA = NUMROC(NKEEP,NB,MYCOL,0,NPCOL) !changed (possibly); local columns in A
      CALL DESCINIT(DESCA,N,NKEEP, MB,NB, 0,0, ICTXT, MAX0(1,MA),INFO)
      IF (INFO /= 0) THEN
         WRITE(*,*)' dist_lrank: descinit 1 mistake'
         IERR = 1
         RETURN
      ENDIF
      LDA = DESCA(LLD_)
      ALLOCATE(A(LDA,MAX(1,NA)),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'dist_lrank: Error setting space for A2'
         RETURN
      ENDIF
      A(:,:) = 0.D0
      write(*,*) 'nkeep=',nkeep
!
!.... recall the eigenvalues are in ascending order, so work backwards
      IF (MYROW == 0 .AND. MYCOL == 0) &
      WRITE(*,*) 'dist_lrank: Reducing U...'
      J1 = 0
      DO 23 J=N,N-NKEEP+1,-1 !loop on columns
         IF (W(J).GT.CUT) THEN
            J1 = J1 + 1
            DO 24 I=1,N !save this row
               CALL INFOG2L(I,J ,DESCZ, NPROW,NPCOL, MYROW,MYCOL, IIA, JJA,  &
                            IAROW,IACOL)
               CALL INFOG2L(I,J1,DESCA, NPROW,NPCOL, MYROW,MYCOL, IIAD,JJAD, &
                            JAROW,JACOL)
               IF (MYROW == IAROW .AND. MYCOL == IACOL) THEN
                  IF (MYROW == JAROW .AND. MYCOL == JACOL) THEN
                     A(IIAD,JJAD) = Z(IIA,JJA)
                  ELSE
                     CALL DGESD2D(ICTXT,1,1,Z(IIA,JJA),1, JAROW,JACOL)
                  ENDIF
               ELSE
                  IF (MYROW == JAROW .AND. MYCOL == JACOL) & !receive from iarow,iacol
                  CALL DGERV2D(ICTXT,1,1,A(IIAD,JJAD),1, IAROW,IACOL)
               ENDIF

   24       CONTINUE
         ENDIF
         CALL BLACS_BARRIER(ICTXT,'A')
   23 CONTINUE
!
!.... reset Z to carry truncated V
      DEALLOCATE(Z) 
      MZ = NUMROC(NKEEP,MB,MYROW,0,NPROW) !number of local rows in Z 
      NZ = NUMROC(N    ,NB,MYCOL,0,NPCOL) !number of local columsn in Z
      CALL DESCINIT(DESCZ,NKEEP,N, MB,NB, 0,0, ICTXT, MAX0(1,MZ),INFO)
      IF (INFO /= 0) THEN 
         WRITE(*,*)' dist_lrank: descinit 2 mistake'
         IERR = 1
         RETURN
      ENDIF
      LDZ = DESCZ(LLD_)
      ALLOCATE(Z(LDZ,MAX(1,NZ)),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'dist_lrank: Error setting space for Z2'
         RETURN
      ENDIF
      Z(:,:) = 0.D0
      print *, 'zset',nkeep

!
!.... multiply U S and copy U -> V whilst fixing signs
      J1 = N + 1
      DO 26 J=1,NKEEP
         J1 = J1 - 1 
         IF (W(J1) < 0.D0) THEN
            XFACT =-1.D0
         ELSE
            XFACT = 1.D0
         ENDIF
         DO 27 I=1,N !loop on column and scale by j'th singular value
            CALL INFOG2L(I,J,DESCA, NPROW,NPCOL, MYROW,MYCOL, IIA, JJA,  &
                         IAROW,IACOL)
            CALL INFOG2L(J,I,DESCZ, NPROW,NPCOL, MYROW,MYCOL, IIAD,JJAD, &
                         JAROW,JACOL)  
            IF (MYROW == IAROW .AND. MYCOL == IACOL) THEN
               IF (MYROW == JAROW .AND. MYCOL == JACOL) THEN
                  Z(IIAD,JJAD) = A(IIA,JJA)
               ELSE
                  CALL DGESD2D(ICTXT,1,1,A(IIA,JJA),1, JAROW,JACOL)
               ENDIF
            ELSE
               IF (MYROW == JAROW .AND. MYCOL == JACOL) & !receive from iarow,iacol
               CALL DGERV2D(ICTXT,1,1,Z(IIAD,JJAD),1, IAROW,IACOL)
            ENDIF
            !each column is scaled by the singular value
            IF (MYROW == IAROW .AND. MYCOL == IACOL) A(IIA ,JJA ) = A(IIA,JJA)*S(J) 
            !each column of V (row of V^T) is multiplied by -sign(eigenvalue) 
            IF (MYROW == JAROW .AND. MYCOL == JACOL) Z(IIAD,JJAD) = XFACT*Z(IIAD,JJAD)
   27    CONTINUE
         CALL BLACS_BARRIER(ICTXT,'A')
   26 CONTINUE
      print *, 'z ifinished',myrow,mycol
!
!.... multiply USV^T -> Hmat
      IF (MYROW == 0 .AND. MYCOL == 0) &
      WRITE(*,*) 'dist_lrank: Multiplying H = U_p S_p V_p^T...'
      CALL BLACS_BARRIER(ICTXT,'A')
      CALL PDGEMM('N','N', N,N,NKEEP, ONE,A,1,1,DESCA, Z,1,1,DESCZ, & 
                  ZERO, HMAT,1,1,DESCH)
      DEALLOCATE(A)
      DEALLOCATE(Z)
      DEALLOCATE(S)
      DEALLOCATE(W)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE TSVD_DIST(LDH, NPROW,NPCOL,ICTXT, LSURF,K, THRESH, DESCH,G,  &
                           HMAT, SEARCH,IERR) 
!
!     Solves Re{adj(J) J} dm =-g via the truncated singular value decomposition.  
!     This may be useful for surface waves.  By theorem 3.3 of Demmel if A is 
!     symmetric with svd A = U S V^T then it also has an eigendecomposition 
!     A = U L U^T.  The signgular values are abs(L) and the right singular
!     vectors are sign(L_i) u_i.  Of course, we ignore any singular vectors 
!     with singular values less than 0 so the last point is moot. 

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: HMAT(LDH,*), G(*), THRESH
      INTEGER*4, INTENT(IN) :: DESCH(9), LDH, NPROW, NPCOL, ICTXT, K
      LOGICAL*4, INTENT(IN) :: LSURF
      REAL*8, INTENT(OUT) :: SEARCH(*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: A(:,:), Z(:,:), W(:), S(:), WORK(:), GLOC(:), YLOC(:), Y(:) 
      CHARACTER(80) FLNAME
      CHARACTER(5) CK 
      REAL*8 WORK8, CUT
      INTEGER*4 DESCA(9), DESCZ(9), DESCG(9), DESCY(9), LWORK, MH, NH, MA, NA, MB, NB, &
                I, J, I1, J1, IAROW, IACOL, JAROW, JACOL, IIA, JJA, IIAD, JJAD, LDA,  &
                MYROW, MYCOL, INFO, N, M, NKEEP, NY, NG
      LOGICAL*4 LDEBUG
      INTEGER*4 NUMROC !calculate NUmber of Rows Or Columns
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      PARAMETER(LDEBUG = .FALSE.)

!
!----------------------------------------------------------------------------------------!
!
!.... figure out cutoff
      IERR = 0
      CUT = THRESH
      IF (THRESH.LE.0.D0) CUT = EPSILON(1.D0)
!
!.... set space and copy upper triangle
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
      M = DESCH(M_)
      N = DESCH(N_)
      LDA = DESCH(LLD_)
      MH = NUMROC(M, DESCH(MB_),MYROW, 0,NPROW) !number of rows in local H matrix
      NH = NUMROC(N, DESCH(NB_),MYCOL, 0,NPCOL) !number of columns in local H matrix
      IF (LDEBUG) GOTO 1001 !debugging
      ALLOCATE(A(LDA,MAX(1,NH)),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'tsvd_dist: Error setting space for A'
         RETURN
      ENDIF
      ALLOCATE(Z(LDA,MAX(1,NH)),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'tsvd_dist: Error setting space for Z'
         RETURN
      ENDIF
      ALLOCATE(W(N))
      A(:,:) = 0.D0
      DESCA(1:9) = DESCH(1:9)
      DESCZ(1:9) = DESCH(1:9)
      DO 1 I=1,N
         DO 2 J=I,N
            CALL INFOG2L(I,J ,DESCZ, NPROW,NPCOL, MYROW,MYCOL, IIA, JJA,  &
                         IAROW,IACOL)
            IF (MYROW == IAROW .AND. MYCOL == IACOL) &
            A(IIA,JJA) = HMAT(IIA,JJA)
    2    CONTINUE
    1 CONTINUE
!
!.... calculate the eigenvalue decomposition
      IF (MYROW == 0 .AND. MYCOL == 0) &
      WRITE(*,*) 'tsvd_dist: Calculating eigendecomposition...'
      LWORK =-1 
      CALL PDSYEV('V','U',N,A,1,1,DESCA, W,  &
                  Z,1,1, DESCZ, WORK8,LWORK, INFO)
      LWORK = INT(WORK8)
      ALLOCATE(WORK(LWORK),STAT=IERR)
      IF (IERR /= 0) THEN
         IF (MYROW == 0 .AND. MYCOL == 0) &
         WRITE(*,*) 'tsvd_dist: Error allocating space for work',IERR
         RETURN
      ENDIF
      CALL PDSYEV('V','U',N,A,1,1,DESCA, W,  &
                  Z,1,1, DESCZ, WORK ,LWORK, INFO)
      IF (INFO /= 0) THEN
         IF (MYROW == 0 .AND. MYCOL == 0) &
         WRITE(*,*) 'tsvd_dist: Error calling pdysev',INFO
         IERR = 1
         RETURN
      ENDIF
      DEALLOCATE(WORK) !done with workspace 
!
!.... count the singular values and reset A
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         FLNAME(1:80) = ' '
         CK(1:5) = ' '
         WRITE(CK,'(I5)') K
         CK = ADJUSTL(CK)
         IF (LSURF) THEN
            FLNAME = 'eigen_spc_srf-'//TRIM(CK)//'.txt'
         ELSE
            FLNAME = 'eigen_spc_bdy-'//TRIM(CK)//'.txt'
         ENDIF
         FLNAME = ADJUSTL(FLNAME)
         OPEN(UNIT=64,FILE=TRIM(FLNAME))
      ENDIF
      ALLOCATE(S(N))
      S(1:N) = 0.D0
      I1 = 0
      DO 21 I=N,1,-1
         IF (MYROW == 0 .AND. MYCOL == 0) WRITE(64,*) W(I)
         IF (W(I).GT.CUT) THEN
            I1 = I1 + 1
            S(I1) = DABS(W(I))
         ENDIF
   21 CONTINUE
      IF (MYROW == 0 .AND. MYCOL == 0) CLOSE(64)
      NKEEP = I1
!
!.... reset vectors
      DEALLOCATE(A) 
      MB = DESCH(MB_)
      NB = MIN(DESCH(NB_),NKEEP) !may be way less columsn now
      MA = NUMROC(N,    MB,MYROW,0,NPROW) !unchange; local rows in A
      NA = NUMROC(NKEEP,NB,MYCOL,0,NPCOL) !changed (possibly); local columns in A
      CALL DESCINIT(DESCA,N,NKEEP, MB,NB, 0,0, ICTXT, MAX0(1,MA),INFO)
      IF (INFO /= 0) THEN
         WRITE(*,*)' tsvd_dist: descinit 1 mistake'
         IERR = 1
         RETURN
      ENDIF
      NY = NUMROC(NKEEP,NB,MYROW,0,NPROW) !number of local rows in y 
      CALL DESCINIT(DESCY,NKEEP,1, NB, 1, 0,0, ICTXT, MAX0(1,NY),INFO)
      IF (INFO /= 0) THEN
         WRITE(*,*) 'tsvd_dist: descinit 2 mistake'
         IERR = 1
         RETURN
      ENDIF 
      NG = NUMROC(N,MB,MYROW,0,NPROW) !number of local rows in gradient 
      CALL DESCINIT(DESCG,N    ,1, MB, 1, 0,0, ICTXT, MAX0(1,NG),INFO)
      IF (INFO /= 0) THEN
         WRITE(*,*) 'tsvd_dist: descinit 3 mistake'
         IERR = 1
         RETURN
      ENDIF
!  MH = NUMROC(N, MBH,MYROW, 0,NPROW) !number of rows in local H
!  NH = NUMROC(N, NBH,MYCOL, 0,NPCOL) !number of columns in local H
!  CALL DESCINIT(DESCH,N,N, MBH,NBH, 0,0, ICTXT, MAX0(1,MH),INFO)
!  CALL DESCINIT(DESCB,N,1, MBH,  1, 0,0, ICTXT, MAX0(1,MH),INFO)

      LDA = DESCA(LLD_)
      ALLOCATE(A(LDA,MAX(1,NA)),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'tsvd_dist: Error setting space for A2'
         RETURN
      ENDIF
      ALLOCATE(YLOC(MAX(1,NY)))
      ALLOCATE(GLOC(MAX(1,NG)))
      CALL SET_SOL8(ICTXT,NPROW,NPCOL, N, DESCG, G, GLOC)
      A(:,:) = 0.D0
      write(*,*) 'nkeep=',nkeep
!
!.... recall the eigenvalues are in ascending order, so work backwards
      IF (MYROW == 0 .AND. MYCOL == 0) &
      WRITE(*,*) 'tsvd_dist: Reducing U...'
      J1 = 0
      DO 23 J=N,N-NKEEP+1,-1 !loop on columns
         IF (W(J).GT.CUT) THEN
            J1 = J1 + 1
            DO 24 I=1,N !save this row
               CALL INFOG2L(I,J ,DESCZ, NPROW,NPCOL, MYROW,MYCOL, IIA, JJA,  &
                            IAROW,IACOL)
               CALL INFOG2L(I,J1,DESCA, NPROW,NPCOL, MYROW,MYCOL, IIAD,JJAD, &
                            JAROW,JACOL)
               IF (MYROW == IAROW .AND. MYCOL == IACOL) THEN
                  IF (MYROW == JAROW .AND. MYCOL == JACOL) THEN
                     A(IIAD,JJAD) = Z(IIA,JJA)
                  ELSE
                     CALL DGESD2D(ICTXT,1,1,Z(IIA,JJA),1, JAROW,JACOL)
                  ENDIF
               ELSE 
                  IF (MYROW == JAROW .AND. MYCOL == JACOL) & !receive from iarow,iacol
                  CALL DGERV2D(ICTXT,1,1,A(IIAD,JJAD),1, IAROW,IACOL)
               ENDIF
   24       CONTINUE
         ENDIF
         CALL BLACS_BARRIER(ICTXT,'A') 
   23 CONTINUE
!
!.... multiply U^T g; [n_keep x n] x [n x 1] = [n_keep x 1]
      CALL PDGEMV('T',N,NKEEP, 1.D0,A,1,1,DESCA, GLOC,1,1,DESCG,1, 0.D0, &
                  YLOC,1,1,DESCY,1)
!
!.... invert
      DO 31 I=1,NKEEP
         CALL INFOG2L(I,1,DESCY, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                      IAROW,IACOL)
         IF (MYROW == IAROW .AND. MYCOL == IACOL) &
         YLOC(IIA) = YLOC(IIA)/S(I)
   31 CONTINUE
!
!.... get the signs right on V
      DO 32 J=1,N !nkeep !loop on columns; n in case we want V^T V res matrix
         IF (W(J) < 0.D0) THEN
            DO 33 I=1,N !negate row
               CALL INFOG2L(I,J,DESCA, NPROW,NPCOL, MYROW,MYCOL, IIA, JJA, &
                            IAROW,IACOL)
               IF (MYROW == IAROW .AND. MYCOL == IACOL) A(IIA,JJA) =-A(IIA,JJA)
   33       CONTINUE
         ENDIF
   32 CONTINUE
!
!.... multiply by transpose(U) [n_keep x n]; negative for minimizing direction
      CALL PDGEMV('N',N,NKEEP,-1.D0,A,1,1,DESCA, YLOC,1,1,DESCY,1, 0.D0, &
                  GLOC,1,1,DESCG,1)
      CALL GET_SOL8(ICTXT,NPROW,NPCOL, N, DESCG, GLOC, SEARCH) !send p back to (0,0)
!
!.... clean
      DEALLOCATE(A)
      DEALLOCATE(Z)
      DEALLOCATE(YLOC)
      DEALLOCATE(GLOC)
      DEALLOCATE(S)
      DEALLOCATE(W)
!
!.... debug 
 1001 CONTINUE
      IF (LDEBUG) THEN
         IF (MYROW == 0 .AND. MYCOL == 0) THEN
            ALLOCATE(A(N,N),STAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'tsvd_dist: Debug failed to allocate A'
               RETURN
            ENDIF
            ALLOCATE(W(N))
            LWORK =-1
            CALL DSYEV('V','U',N, A,N, W,WORK8,LWORK, INFO)
            LWORK = INT(WORK8)  
            ALLOCATE(WORK(LWORK),STAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'tsvd_dist: Debug failed to allocated work'
               RETURN
            ENDIF
            WRITE(*,*) 'tsvd_dist: Getting A...' 
            A(:,:) = 0.D0
         ENDIF !end check on (0,0) 
         DO 101 I=1,N
            DO 102 J=I,N
               CALL INFOG2L(I,J,DESCH, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                            IAROW,IACOL)
               IF (MYROW == 0 .AND. MYCOL == 0) THEN
                  IF (MYROW == IAROW .AND. MYCOL == IACOL) THEN
                     A(I,J) = HMAT(IIA,JJA)
                  ELSE
                     CALL DGERV2D(ICTXT,1,1,A(I,J),1, IAROW,IACOL)
                  ENDIF
               ELSE
                  IF (MYROW == IAROW .AND. MYCOL == IACOL) &
                  CALL DGESD2D(ICTXT,1,1,HMAT(IIA,JJA),1, 0,0)
               ENDIF 
  102       CONTINUE
            CALL BLACS_BARRIER(ICTXT,'A')
  101    CONTINUE  
         IF (MYROW == 0 .AND. MYCOL == 0) THEN
            WRITE(*,*) 'tsvd_dist: Debug performing eigendecomposition...' 
            CALL DSYEV('V','U',N, A,N, W,WORK,LWORK, INFO) 
            IF (INFO /= 0) THEN
               WRITE(*,*) 'tsvd_dist: Eigendecomposition failed!'
               RETURN
            ENDIF
            nkeep = 0
            do i=1,n
               if (w(i) > cut) nkeep = nkeep + 1
            enddo
            print *, nkeep
!
!.......... multiply U^T g
            ALLOCATE(Y(NKEEP))
            I1 = N - NKEEP + 1
            !CALL DGEMV('T',N,NKEEP, 1.D0,A(1:N,I1:N),N, G,1, 0.D0,Y,1)
            Y = MATMUL(TRANSPOSE(A(1:N,I1:N)),G(1:N))
!
!.......... scale
            I1 = N + 1
            DO 103 I=NKEEP,1,-1
               I1 = I1 - 1
               Y(I) = Y(I)/W(I1)
  103       CONTINUE
!
!.......... fix V
            J1 = N + 1 
            DO 104 J=1,N !loop on columns
               J1 = J1 - 1
               write(34,*) w(j)
               IF (W(J1) < 0.D0) THEN
                  DO 105 I=1,N !negate row of U
                     A(I,J) =-A(I,J)
  105             CONTINUE
               ENDIF
  104       CONTINUE
!
!.......... multiply Vy =-search
            I1 = N - NKEEP + 1
            !CALL DGEMV('N',N,NKEEP,-1.D0,A(1:N,I1:N),N, Y,1, 0.D0,SEARCH,1)
            SEARCH(1:N) =-MATMUL(A(1:N,I1:N),Y)
            DEALLOCATE(Y)
            DEALLOCATE(A) 
         ENDIF !end check on myrow/mycol
      ENDIF !end check on debugging
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE PLOT_JACOB_DRIVE(NREC_LIST,NREC,NPPGRP,LSTF,  &
                                  MSH,SRC,FRQ,INV, IREC_LIST, IERR)

      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO) MSH 
      TYPE (SRC_INFO)  SRC
      TYPE (FRQ_INFO)  FRQ
      TYPE (INV_INFO)  INV 
      INTEGER*4, INTENT(IN) :: IREC_LIST(NREC_LIST), NREC_LIST, NPPGRP, NREC
      LOGICAL*4, INTENT(IN) :: LSTF
      INTEGER*4, INTENT(OUT) :: IERR  

      REAL*4, ALLOCATABLE :: JMATL(:,:), RJAC(:,:), QJAC(:,:) 
      INTEGER*4, ALLOCATABLE :: IPERM_REC(:)
      COMPLEX*8 CRES
      CHARACTER(80) FLNAME 
      CHARACTER(12) CFREQ
      CHARACTER(5)  CSRC, CID 
      INTEGER*4 M,N,LDJAC, IFREQ,ISRC,IPPGRP, JREC,IREC, IA35, IROW,JROW, I 
!
!.... set space
      IERR = 0
      M = 2*NDIM*NREC   !number of observations (real,complex) 
      N = inv%NA35      !number of inversion points
      LDJAC = M     !leading dimension 
      ALLOCATE(JMATL(LDJAC,N))
      ALLOCATE(RJAC(NDIM,N))
      ALLOCATE(QJAC(NDIM,N))
      ALLOCATE(IPERM_REC(NDIM*NREC))
!
!.... loop on frequencies
      DO 1 IFREQ=1,frq%NFREQ
         WRITE(*,*) 'plot_jacob_drive: Processing frequency:',frq%FREQ(IFREQ)
         CFREQ(1:12) = ' '
         WRITE(CFREQ,'(F12.5)') frq%FREQ(IFREQ)
         CFREQ = ADJUSTL(CFREQ)
!
!....... loop on sources
         DO 2 ISRC=1,src%NSRC
            CSRC(1:5) = ' '
            WRITE(CSRC,'(I5)') ISRC
            CSRC = ADJUSTL(CSRC)
            IF (MAXVAL(CABS(inv%OBS(1:NDIM,IFREQ,1:NREC,ISRC))) == 0.0) GOTO 20
            CALL PERM_OBS(NDIM,NDIM,NREC, inv%OBS(1:NDIM,IFREQ,1:NREC,ISRC), &
                          IPERM_REC,IERR)
            IF (IERR /= 0) THEN 
               WRITE(*,*) 'muljac_resid: Error calling perm_obs!'
               RETURN
            ENDIF
!
!.......... load the jacobians
            DO 3 IPPGRP=1,NPPGRP
               CID(1:5) = ' '
               WRITE(CID,'(I5)') IPPGRP - 1
               CID = ADJUSTL(CID)
               FLNAME(1:80) = ' '
               FLNAME = './scratch/jac_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'// &
                       TRIM(CSRC)//'.dat'
               FLNAME = ADJUSTL(FLNAME)
               CALL LOAD_FJAC(FLNAME,LDJAC, NDIM,NREC, IPERM_REC, JMATL,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'plot_jacob_drive: Error unable to load Jacobian', &
                             frq%FREQ(IFREQ),ISRC
                  RETURN
               ENDIF
    3       CONTINUE
!
!.......... extract and write
            DO 4 JREC=1,NREC_LIST
               IREC = IREC_LIST(JREC)
               DO 5 I=1,NDIM
                  IROW = (IREC - 1)*NDIM + I
                  JROW = IROW + NDIM*NREC 
                  DO 6 IA35=1,inv%NA35 
                     IF (LSTF) THEN
                        CRES = src%SOURCE(IFREQ,ISRC) &
                              *CMPLX(JMATL(IROW,IA35),JMATL(JROW,IA35))
                        RJAC(I,IA35) = REAL(CRES)
                        QJAC(I,IA35) = IMAG(CRES) 
                     ELSE
                        RJAC(I,IA35) = JMATL(IROW,IA35) 
                        QJAC(I,IA35) = JMATL(JROW,IA35) 
                     ENDIF 
    6             CONTINUE
    5          CONTINUE
!
!............. write
               CALL PLOT_JACOB(NDIM,NGNOD, NDIM,msh%NNPG,msh%NELEM,IREC,ISRC,  &
                               msh%CDOMAIN,inv%MASKG,msh%IENG, RJAC,QJAC,      &
                               frq%FREQ(IFREQ), msh%XLOCS,msh%ZLOCS) 
    4       CONTINUE
   20       CONTINUE
    2   CONTINUE !loop on sources
    1 CONTINUE !Loop on frequencies
      DEALLOCATE(RJAC)
      DEALLOCATE(QJAC)
      DEALLOCATE(JMATL) 
      DEALLOCATE(IPERM_REC)
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE JLOC2GRAD(MDIM,MFREQ,MREC,                                              &
                           NDIM,NFREQ,NREC,NA35, NSRC,                                   &
                           IRESTP,IBPHASE,NPPGRP,LUNWRAP,LSTF,                           &
                           AZMOD,FREQ, SOURCE, EST,OBS,                                  &
                           GRAD,IERR) 
!
!     From the local Jacobians generates a gradients by matrix multiplication with 
!     the residuals
!
!     INPUT      MEANING
!     -----      ------- 
!     AZMOD      model azimuth (degrees)
!     EST        frequency domain estaimtes (U,V,W) format
!     IBPHASE    geometric spreading correction for residuals
!     IRESTP     residual type (1 phase only, 2 amplitude only, 3 phase and amplitude)
!     FREQ       frequency list
!     LSTF       True -> convolve source time function
!     LUNWRAP    True -> estimates,observations in amplitude phase format
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MREC       leading dimension
!     NA35       number of inversion variables
!     NDIM       number of components
!     NFREQ      number of frequencies 
!     NPPGRP     number of processes per process group
!     NREC       number of receivers
!     NSRC       number of sources
!     OBS        frequency domain observations (N,E,Z) format
!     SOURCE     source time function
!
!     OUTPUT     MEANING
!     ------     ------- 
!     GRAD       gradient 
!     IERR       error flag
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: EST(MDIM,MFREQ,MREC,*), OBS(MDIM,MFREQ,MREC,*),  &
                               SOURCE(MFREQ,*) 
      REAL*8, INTENT(IN) :: FREQ(NFREQ), AZMOD  
      INTEGER*4, INTENT(IN) :: MDIM, MFREQ, MREC, NDIM, NFREQ, NREC, NA35, NSRC, NPPGRP, &
                                     IRESTP, IBPHASE
      LOGICAL*4, INTENT(IN) :: LUNWRAP, LSTF
      REAL*4, INTENT(OUT) :: GRAD(NA35) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables 
      CHARACTER(80) FLNAME 
      CHARACTER(12) CFREQ 
      CHARACTER(5) CID, CSRC  
      REAL*4, ALLOCATABLE :: JMATL(:,:), RESIDV(:)
      INTEGER*4, ALLOCATABLE :: IPERM_REC(:)
      COMPLEX*8 CRESID,HDIFFER,CPHM2CM, Q1,Q2,Q3, RN,RE, OBSN,OBSE,OBSZ, QN,QE,QZ, &
                HFACT  
      REAL*4 UMAG,VMAG,WMAG, NMAG,EMAG,ZMAG, UPH,VPH,WPH, NPH,EPH,ZPH, ALPHA  
      INTEGER*4 M, N, LDJAC, ML, IFREQ, IREC, ISRC, I, I1, I2, IPPGRP, INDX, JNDX
      LOGICAL*4 LOBS, LZERO 
!
!----------------------------------------------------------------------------------------! 
!
!.... initialize
      IERR = 0
      M = 2*NDIM*NREC   !number of observations (real,complex) 
      N = NA35          !number of inversion points
      LDJAC = M     !leading dimension 
      ALLOCATE(JMATL(LDJAC,N)) 
      ALLOCATE(RESIDV(M))   !real then complex block 
      ALLOCATE(IPERM_REC(NDIM*NREC))
      GRAD(1:NA35) = 0.0 !0 out b/c we will be summing 
      DO IREC=1,NREC
         DO I=1,NDIM
            IPERM_REC( (IREC - 1)*NDIM + I) = (IREC - 1)*NDIM + I
         ENDDO
      ENDDO
!
!.... loop on frequencies
      LZERO = .TRUE.
      DO 1 IFREQ=1,NFREQ
         WRITE(*,*) 'jloc2grad: Processing frequency:',FREQ(IFREQ)
         CFREQ(1:12) = ' '
         WRITE(CFREQ,'(F12.5)') FREQ(IFREQ)
         CFREQ = ADJUSTL(CFREQ)
         HFACT = HDIFFER(IBPHASE,FREQ(IFREQ)) 
!
!....... loop on sources
         DO 2 ISRC=1,NSRC
!
!.......... zero source time function
            IF (SOURCE(IFREQ,ISRC) == CMPLX(0.0,0.0)) THEN
               WRITE(*,*) 'jloc2grad: Warning STF is zero, skipping...'
               GOTO 200
            ENDIF
!
!.......... set source character name
            CSRC(1:5) = ' '
            WRITE(CSRC,'(I5)') ISRC
            CSRC = ADJUSTL(CSRC)
!
!.......... check we have observations for this source
            LOBS = .FALSE.
            DO 3 IREC=1,NREC 
               DO 4 I=1,NDIM
                  IF (CABS(OBS(I,IFREQ,IREC,ISRC)) > 0.0) THEN
                     LOBS = .TRUE.
                     GOTO 30
                  ENDIF  
    4          CONTINUE
    3       CONTINUE
            WRITE(*,*) 'jloc2grad: No observations for frequency source pair:', &
                        FREQ(IFREQ),ISRC
            GOTO 200
   30       CONTINUE !observation exists, Jacobian should exist
!
!.......... load the jacobians
            DO 5 IPPGRP=1,NPPGRP
               CID(1:5) = ' '
               WRITE(CID,'(I5)') IPPGRP - 1
               CID = ADJUSTL(CID)
               FLNAME(1:80) = ' '
               FLNAME = './scratch/jac_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'// &
                       TRIM(CSRC)//'.dat'
               FLNAME = ADJUSTL(FLNAME)
               CALL LOAD_FJAC(FLNAME,LDJAC, NDIM,NREC, IPERM_REC, JMATL,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'jloc2grad: Error unable to load Jacobian',FREQ(IFREQ),ISRC
                  RETURN 
               ENDIF
    5       CONTINUE  
!
!.......... fill in the residual vectors
            INDX = 0
            JNDX = NREC*NDIM
            DO 6 IREC=1,NREC
               IF (LUNWRAP) THEN
                  UMAG = REAL(EST(1,IFREQ,IREC,ISRC))
                  VMAG = REAL(EST(2,IFREQ,IREC,ISRC)) 
                  WMAG = REAL(EST(3,IFREQ,IREC,ISRC)) 
                  UPH  = IMAG(EST(1,IFREQ,IREC,ISRC))
                  VPH  = IMAG(EST(2,IFREQ,IREC,ISRC))
                  WPH  = IMAG(EST(3,IFREQ,IREC,ISRC)) 
                  Q1 = CPHM2CM(UMAG,UPH)
                  Q2 = CPHM2CM(VMAG,VPH)
                  QZ = CPHM2CM(WMAG,WPH)  

                  NMAG = REAL(OBS(1,IFREQ,IREC,ISRC))
                  EMAG = REAL(OBS(2,IFREQ,IREC,ISRC))
                  ZMAG = REAL(OBS(3,IFREQ,IREC,ISRC))
                  NPH  = IMAG(OBS(1,IFREQ,IREC,ISRC))
                  EPH  = IMAG(OBS(2,IFREQ,IREC,ISRC))
                  ZPH  = IMAG(OBS(3,IFREQ,IREC,ISRC))
                  OBSN = CPHM2CM(NMAG,NPH) 
                  OBSE = CPHM2CM(EMAG,EPH)
                  OBSZ = CPHM2CM(ZMAG,ZPH)
               ELSE
                  Q1 = EST(1,IFREQ,IREC,ISRC) 
                  Q2 = EST(2,IFREQ,IREC,ISRC)
                  QZ = EST(3,IFREQ,IREC,ISRC)
                  OBSN = EST(1,IFREQ,IREC,ISRC)
                  OBSE = EST(2,IFREQ,IREC,ISRC)
                  OBSZ = EST(3,IFREQ,IREC,ISRC)
               ENDIF
!
!............. rotate (u,v,w) -> (N,E,Z)
               CALL CROTATE( SNGL(AZMOD),Q1,Q2, QN,QE)
!
!............. calculate residuals
               RN = CRESID(IRESTP,OBSN,QN)
               RE = CRESID(IRESTP,OBSE,QE)
               Q3 = CRESID(IRESTP,OBSZ,QZ)
!
!............. rotate residuals to (u,v,w) frame for backpropagation 
               CALL CROTATE(-SNGL(AZMOD),RN,RE, Q1,Q2) 
               Q1 = Q1*HFACT 
               Q2 = Q2*HFACT 
               Q3 = Q3*HFACT

               INDX = INDX + 1 
               JNDX = JNDX + 1 
               RESIDV(INDX) = REAL(Q1)
               RESIDV(JNDX) = IMAG(Q1)
               INDX = INDX + 1
               JNDX = JNDX + 1 
               RESIDV(INDX) = REAL(Q2)
               RESIDV(JNDX) = IMAG(Q2)
               INDX = INDX + 1
               JNDX = JNDX + 1
               RESIDV(INDX) = REAL(Q3)
               RESIDV(JNDX) = IMAG(Q3) 
    6       CONTINUE !loop on receivers
!
!.......... convolve STF
            IF (LSTF) CALL SCALE_JLOC(LDJAC,NREC*NDIM,NA35, SOURCE(IFREQ,ISRC),JMATL)
!
!.......... now update grad = grad - Re[ adj(J) delta d ] =-Re(J^T) Re(d) - Im(J^T) Im(d)
            I1 = 1
            I2 = NDIM*NREC 
            ML = I2 - I1 + 1
            ALPHA =-1.0 
            CALL SGEMV('T',ML,N, ALPHA,JMATL(I1:I2,:),ML, RESIDV(I1:I2),1,1.0,GRAD,1)
            !GRAD(:) = GRAD(:) + MATMUL(TRANSPOSE(JMATL(I1:I2,:)),RESIDV(I1:I2))
            I1 = NDIM*NREC + 1
            I2 = 2*NDIM*NREC 
            ML = I2 - I1 + 1 
            CALL SGEMV('T',ML,N, ALPHA,JMATL(I1:I2,:),ML, RESIDV(I1:I2),1,1.0,GRAD,1)
            !GRAD(:) = GRAD(:) - MATMUL(TRANSPOSE(JMATL(I1:I2,:)),RESIDV(I1:I2)) 
            LZERO = .FALSE.
  200       CONTINUE !break ahead for no observations
    2    CONTINUE !loop on sources
    1 CONTINUE !loop on frequencies
!
!.... make sure gradient isnt zero 
      IF (LZERO) THEN
         WRITE(*,*) 'jloc2grad: Error gradient is 0!'
         IERR = 1
      ENDIF 
!
!.... free space
      DEALLOCATE(RESIDV)
      DEALLOCATE(JMATL) 
      DEALLOCATE(IPERM_REC)
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE SETRHS_DIST_428(LDB,M, NPROW,NPCOL,ICTXT, DESCB, RHS_GLOB, RHS) 
!  
!     Sets the distributed RHS 
!
!     INPUT      MEANING
!     -----      ------- 
!     DESCB      BLACS descriptor for RHS 
!     ICTXT      BLACS context
!     LDB        leading dimension of B; workspace size
!     M          number of rows
!     NPCOL      number of process columns
!     NPROW      number of process rows 
!     RHS_GLOB   holds the global RHS on (0,0) to distribute
!
!     OUTPUT     MEANING
!     ------     ------- 
!     RHS        local RHS 
!
!.... variable declarations
      IMPLICIT NONE
      REAL*4, INTENT(IN) :: RHS_GLOB(*)
      INTEGER*4, INTENT(IN) :: DESCB(9), LDB, M, NPROW, NPCOL, ICTXT
      REAL*8, INTENT(OUT) :: RHS(LDB)
!.... local variables
      REAL*4 VAL
      INTEGER*4 I, MYROW, MYCOL, IIA, JJA, IAROW, IACOL 
!
!----------------------------------------------------------------------------------------!
!
!.... get my rank in the BLACS grid
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!
!.... loop on rows
      DO 1 I=1,M
         CALL INFOG2L(I,1,DESCB, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                      IAROW,IACOL) 
         IF (MYROW == 0 .AND. MYCOL == 0) THEN
            IF (IAROW == 0 .AND. IACOL == 0) THEN !belongs to me
               RHS(IIA) = DBLE(RHS_GLOB(I))
            ELSE !send
               CALL SGESD2D(ICTXT,1,1,RHS_GLOB(I),1, IAROW,IACOL) !send to (iarow,iacol) 
            ENDIF
         ELSE
            IF (MYROW == IAROW .AND. MYCOL == IACOL) THEN
               CALL SGERV2D(ICTXT,1,1,VAL,1, 0,0) !receiver from (0,0)
               RHS(IIA) = DBLE(VAL)
            ENDIF
         ENDIF 
    1 CONTINUE
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE SCALE_JLOC(LDJAC,NOBS,NA35, CSCLR,JMATL) 
!
!     Scales the Jacobian by a complex number 
!
!     INPUT      MEANING
!     -----      -------
!     CSCLR      complex scalar to multiply Jacobian by
!     JMATL      initial Jacobian[real,imag derivatives, inversion vars]
!     LDJAC      leading dimension for Jacobian
!     NA35       number of inversion variables
!     NOBS       number of observations 
!    
!     OUTPUT     MEANING
!     ------     -------
!     JMATL      Jacobian scaled by CSCLR
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: CSCLR  
      INTEGER*4, INTENT(IN) :: LDJAC, NOBS, NA35 
      REAL*4, INTENT(INOUT) :: JMATL(LDJAC,*) 
!.... local variables
      COMPLEX*8 CRES
      INTEGER*4 JOBS, IOBS, IA35
!
!----------------------------------------------------------------------------------------!
!
      IF (CSCLR == CMPLX(1.0,0.0)) RETURN !would be just multiplying by one
      DO 1 IOBS=1,NOBS
         JOBS = NOBS + IOBS !offset to imaginary part  
         DO 2 IA35=1,NA35
            CRES = CSCLR*CMPLX(JMATL(IOBS,IA35),JMATL(JOBS,IA35)) 
            JMATL(IOBS,IA35) = REAL(CRES) 
            JMATL(JOBS,IA35) = IMAG(CRES) 
    2    CONTINUE !loop on inversion variables
    1 CONTINUE 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GENHLOC_SB(ICTXT,NPROW,NPCOL, MBJ,NBJ,MBJT,NBJT, DESCH,       &
                            LDH,M,N, LSTF,LSURF, NPPGRP,FRQ,RCV,SRC,INV,       &
                            HLOC, IERR)
!
!     Generates the local approximate Hessian from the saved Jacobians
!
!     INPUT      MEANING
!     -----      -------
!     ICTXT      BLACS context 
!     LDH        leading dimension for local Hessian matrix
!     LSTF       True -> convolve source time function
!     LSURF      True -> generate surface wave approximate hessian 
!     M          number of rows in Jacobian
!     MBJ        block size for jacobian
!     MBJT       block size for rows of transpose jacobian
!     NBJ        block size for jacobian
!     NBJT       block size for columns of transpose jacobian
!     N          number of columns in Jacobian 
!     NPCOL      number of process columns
!     NPROW      number of process rows
!
!     OUTPUT     MEANING
!     ------     -------
!     HLOC       local Hessian matrix
!     IERR       error flag 
! 
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO)  INV
      TYPE (FRQ_INFO)  FRQ
      TYPE (RECV_INFO) RCV
      TYPE (SRC_INFO)  SRC
      INTEGER*4, INTENT(IN) :: DESCH(9), ICTXT, NPROW, NPCOL, MBJ,NBJ,MBJT,NBJT, &
                               LDH, M, N, NPPGRP
      LOGICAL*4, INTENT(IN) :: LSTF, LSURF
      REAL*8, INTENT(INOUT) :: HLOC(LDH,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: FREQ_LOC(:)
      REAL*4, ALLOCATABLE :: WGT(:,:,:,:)
      INTEGER*4, ALLOCATABLE :: IPERM_REC(:)
      CHARACTER(80) FLNAME
      CHARACTER(12) CFREQ
      CHARACTER(5) CID, CSRC
      COMPLEX*8 STF
      REAL*8 VAL, DMAG2, ONE
      INTEGER*4 LDC, IFREQ, ISRC, IPPGRP, I, J, IROW,JCOL, IREC, IIA,JJA,IAROW,IACOL, &
                LDJAC, MM, NM, KM, M2, NFREQ, NSRC, MCL, NCL, ISKIP, JFREQ,JSRC, JREC, &
                LROW, KCOL
      LOGICAL*4 LGN, LBODY 
      PARAMETER(ONE = 1.D0)
      REAL*8, PARAMETER :: ZERO = 0.D0
!.... BLACS variables
      INTEGER*4 DESCJT(9)  !BLACS descriptor for Jacobian^T
      INTEGER*4 DESCJ(9)   !BLACS descriptor for Jacobian 
      INTEGER*4 DESCC(9)   !BLACS descriptor for data covariance matrix
      COMPLEX*8, ALLOCATABLE :: OBS(:,:,:,:)
      REAL*8, ALLOCATABLE :: COVD(:,:), COVD_DIST(:,:), JLOC_CD(:,:) 
      REAL*8, ALLOCATABLE :: JTLOC(:,:) !will hold transpose(Jacobian) 
      REAL*8, ALLOCATABLE :: JLOC(:,:)  !will holds Jacobian
      REAL*4, ALLOCATABLE :: JMATL(:,:) !single precision Jacobian for reading 
      INTEGER*4 MYROW, MYCOL, MJT, MJ, NJT, NJ, NH, INFO, NUMROC
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
!
!----------------------------------------------------------------------------------------!
!
!.... quick return
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
      IF (LSURF .AND. frq%NFREQ_SRF_INV == 0) RETURN
      IF (LBODY .AND. frq%NFREQ_BDY_INV == 0) RETURN
!
!.... check which processes are in the grid
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
      IF (MYROW == -1 .OR. MYCOL == -1) GOTO 1000 !nothing for this process to do 
!
!.... set the frequency list and weights
      IF (LSURF) THEN 
         NFREQ = frq%NFREQ_SRF_INV
         NSRC  = src%NSRC_SRF
      ELSE 
         NFREQ = frq%NFREQ_BDY_INV
         NSRC  = src%NSRC_BDY
      ENDIF
      ALLOCATE(FREQ_LOC(NFREQ)) 
      CALL GET_FREQ_INV(frq%NFREQ, NFREQ, LSURF, frq%CFTYPE,frq%LINVF,frq%FREQ,   &
                        FREQ_LOC,IERR)
      ALLOCATE(OBS(NDIM,NFREQ,rcv%NREC,NSRC))
      ALLOCATE(WGT(NDIM,NFREQ,rcv%NREC,NSRC))
      WGT(:,:,:,:) = 0.0
      JFREQ = 0
      DO 30 IFREQ=1,frq%NFREQ
         IF (     LSURF .AND. frq%CFTYPE(IFREQ) == 'S' .OR.   &
             .NOT.LSURF .AND. frq%CFTYPE(IFREQ) == 'B') THEN
            IF (frq%LINVF(IFREQ)) THEN
               JFREQ = JFREQ + 1
               JSRC = 0
               DO 31 ISRC=1,src%NSRC
                  IF (     LSURF .AND. src%SRCTYP(ISRC)(1:1) == 'S' .OR. &
                      .NOT.LSURF .AND. src%SRCTYP(ISRC)(1:1) == 'P') THEN
                     JSRC = JSRC + 1
                     DO 32 IREC=1,rcv%NREC
                        DO 33 I=1,NDIM
                           OBS(I,JFREQ,IREC,JSRC) = inv%OBS(I,IFREQ,IREC,ISRC) 
                           WGT(I,JFREQ,IREC,JSRC) = inv%WGHTS(I,IFREQ,IREC,ISRC)
   33                   CONTINUE
   32                CONTINUE
                  ENDIF
   31          CONTINUE
            ENDIF
         ENDIF
   30 CONTINUE
!
!.... set sizes and space
      MJT = NUMROC(N, MBJT, MYROW, 0, NPROW) !block rows in J^T  
      MJ  = NUMROC(M, MBJ , MYROW, 0, NPROW) !block rows in J 
      NJT = NUMROC(M, NBJT, MYCOL, 0, NPCOL) !block columns in J^T
      NJ  = NUMROC(N, NBJ , MYCOL, 0, NPCOL) !block columsn in J
      MCL = NUMROC(M, MBJ , MYROW, 0, NPROW) !block rows in C_D^{-1}
      NCL = NUMROC(M, NBJ , MYCOL, 0, NPCOL) !block columns in C_D^{-1}
      CALL DESCINIT(DESCJT,N,M, MBJT,NBJT,0,0, ICTXT, MAX0(1,MJT),INFO)
      IF (INFO /= 0) THEN
         WRITE(*,*) 'genhloc_sb: An error occurred creating DESCAT',MYROW,MYCOL
         IERR = 1
         RETURN
      ENDIF
      CALL DESCINIT(DESCJ ,M,N, MBJ ,NBJ ,0,0, ICTXT, MAX0(1,MJ ),INFO)
      IF (INFO /= 0) THEN
         WRITE(*,*) 'genhloc_sb: An error occurred creating DESCA',MYROW,MYCOL
         IERR = 1
         RETURN
      ENDIF
      CALL DESCINIT(DESCC, M,M, MBJ, NBJ ,0,0, ICTXT, MAX0(1,MCL),INFO)
      IF (INFO /= 0) THEN 
         WRITE(*,*) 'genhloc_sb: Error setting descinit descc'
         IERR = 1
         RETURN 
      ENDIF
      LDC = DESCC(LLD_)
      ALLOCATE(COVD_DIST(LDC,MAX0(1,NCL)))
      ALLOCATE(COVD(LDC,MAX0(1,NCL)))
      ALLOCATE(IPERM_REC(NDIM*rcv%NREC))
      COVD_DIST(1:LDC,1:NCL) = 0.D0
      COVD(1:LDC,1:NCL) = 0.D0 
      IF (     LSURF .AND. .NOT.inv%LCOVD_SRF .OR. &
          .NOT.LSURF .AND. .NOT.inv%LCOVD_BDY) THEN 
!        CALL LAPLACE_DIST1D(ICTXT,NPROW,NPCOL, LDC,rcv%NREC,MCL,NCL,  &
!                            DESCC, inv%DX, COVD_DIST,IERR)
!        IF (IERR /= 0) THEN
!           IERR = 0
!           WRITE(*,*) 'genhloc_sb: Error in laplace_dist1d! Setting to idenity'
!           COVD_DIST(1:LDC,1:NCL) = 0.D0
!           DO I=1,M 
!              CALL INFOG2L(I,I,DESCC, NPROW,NPCOL, MYROW,MYCOL, &
!                           IIA,JJA, IAROW,IACOL)
!              IF (MYROW == IAROW == MYROW .AND. MYCOL == IACOL) COVD_DIST(IIA,JJA) = 1.D0
!           ENDDO 
!        ENDIF
!     ELSE
         COVD_DIST(1:LDC,1:NCL) = 0.D0
         DO I=1,M
            CALL INFOG2L(I,I,DESCC, NPROW,NPCOL, MYROW,MYCOL, &
                         IIA,JJA, IAROW,IACOL)
            IF (MYROW == IAROW .AND. MYCOL == IACOL) COVD_DIST(IIA,JJA) = 1.D0
         ENDDO
      ENDIF

      MM = N
      NM = N
      KM = M
      M2 = M/2 !number of obsrevations 
      ALLOCATE(JTLOC  (DESCJT(LLD_),MAX0(1,NJT)))
      ALLOCATE(JLOC   (DESCJ (LLD_),MAX0(1,NJ )))
      ALLOCATE(JLOC_CD(DESCJ (LLD_),MAX0(1,NJ )))
      LDJAC = M
      ALLOCATE(JMATL(LDJAC,N))
      NH = NUMROC(N, DESCH(NB_),MYCOL, 0,NPCOL) !number of columns in local H matrix
      HLOC(1:LDH,1:MAX(1,NH)) = 0.D0 
      JTLOC(:,:) = 0.D0
      JLOC (:,:) = 0.D0
      JLOC_CD(:,:) = 0.D0
      JMATL(:,:) = 0.0
!
!.... loop on frequencies
      LGN = .FALSE.
      DO 100 IFREQ=1,NFREQ
         IF (MYROW == 0 .AND. MYCOL == 0) THEN
            WRITE(*,*) 'genhloc_sb: Summing frequency:',SNGL(FREQ_LOC(IFREQ))
            CFREQ(1:12) = ' '
            WRITE(CFREQ,'(F12.5)') FREQ_LOC(IFREQ)
            CFREQ = ADJUSTL(CFREQ)
         ENDIF
!
!....... loop on sources
         JSRC = 0
         DO 200 ISRC=1,src%NSRC
            IF (LSURF .AND. src%SRCTYP(ISRC)(1:1) == 'P') GOTO 205 
            IF (LBODY .AND. src%SRCTYP(ISRC)(1:1) == 'S') GOTO 205
            JSRC = JSRC + 1
            IF (MYROW == 0 .AND. MYCOL == 0) THEN
               CSRC(1:5) = ' '
               WRITE(CSRC,'(I5)') ISRC
               CSRC = ADJUSTL(CSRC)
            ENDIF
!
!.......... set peruation vector
            IF (MAXVAL(CABS(OBS(1:NDIM,IFREQ,1:rcv%NREC,JSRC))) == 0.0) GOTO 250
!
!.......... generate a permutation list [w obs, v obs, u obs; null obs]
            CALL PERM_OBS(NDIM,NDIM,rcv%NREC, OBS(1:NDIM,IFREQ,1:rcv%NREC,JSRC), &
                          IPERM_REC,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'genhloc_sb: Error calling perm_obs!'
               RETURN
            ENDIF
            IF (     LSURF .AND. inv%LCOVD_SRF .OR. &
                .NOT.LSURF .AND. inv%LCOVD_BDY) THEN
               CALL LAPLACE_DIST1DV2(ICTXT,NPROW,NPCOL, &
                                     LDC,NDIM,rcv%NREC,MCL,NCL,DESCC, inv%LDWGHT,inv%DX, &
                                     OBS(1:NDIM,IFREQ,1:rcv%NREC,JSRC), COVD_DIST,IERR)
               IF (IERR /= 0) THEN
                  COVD_DIST(1:LDC,1:NCL) = 0.D0
                  DO 201 I=1,M
                     CALL INFOG2L(I,I,DESCC, NPROW,NPCOL, MYROW,MYCOL, &
                                  IIA,JJA, IAROW,IACOL)
                     IF (MYROW == IAROW .AND. MYCOL == IACOL) &
                     COVD_DIST(IIA,JJA) = 1.D0
  201             CONTINUE
               ENDIF
            ENDIF
!
!.......... (0,0) process reads and broadcats
            ISKIP = 0
            IF (MYROW == 0 .AND. MYCOL == 0) THEN
!
!............. check STF is non-zero
               IF (LSTF) THEN
                  IF (src%SOURCE(IFREQ,ISRC) == CMPLX(0.0,0.0)) THEN
                     WRITE(*,*) 'genhloc_sb: STF is zero, skipping!'
                     ISKIP = 1
                     GOTO 350
                  ENDIF
               ENDIF
!
!............. load the Jacobian for this source frequency pair
               JMATL(1:LDJAC,1:N) = 0.0
               DO 300 IPPGRP=1,NPPGRP
                  CID(1:5) = ' '
                  WRITE(CID,'(I5)') IPPGRP - 1
                  CID = ADJUSTL(CID)
                  FLNAME(1:80) = ' '
                  IF (LSURF) THEN
                     FLNAME = './scratch/jac_surf_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'// &
                              TRIM(CSRC)//'.dat'
                  ELSE
                     FLNAME = './scratch/jac_body_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'// &
                              TRIM(CSRC)//'.dat'
                  ENDIF
                  FLNAME = ADJUSTL(FLNAME)
                  !WRITE(*,*) 'genhloc_sb: Reading:',TRIM(FLNAME)
                  CALL LOAD_FJAC(FLNAME,LDJAC, NDIM,rcv%NREC, IPERM_REC, JMATL,IERR)
                  IF (IERR < 0) THEN
                     WRITE(*,*)'genhloc_sb: No local Jacobian for frequency source pair',&
                                SNGL(frq%FREQ(IFREQ)),ISRC, TRIM(FLNAME)
                     ISKIP = 1
                     IERR = 0
                  ENDIF
                  IF (IERR > 0) THEN
                     WRITE(*,*) 'genhloc_sb: Error reading file:'
                     RETURN
                  ENDIF
  300          CONTINUE !loop on process groups
!
!............. adj( s(w) J) s(w) J = |s|^2 adj(J) J
               IF (LSTF) THEN
                  STF = src%SOURCE(IFREQ,ISRC)
                  DMAG2 = DBLE( REAL(CONJG(STF)*STF) )
               ELSE
                  DMAG2 = 1.D0
               ENDIF
  350          CONTINUE !break ahead for no source time function 
               CALL IGEBS2D(ICTXT,'A','I', 1,1, ISKIP,1)
               IF (ISKIP /= 1) THEN
                  CALL DGEBS2D(ICTXT,'A','I', 1,1, DMAG2,1)
                  CALL SGEBS2D(ICTXT,'A','I', M,N, JMATL,LDJAC)
               ENDIF
            ELSE
               CALL IGEBR2D(ICTXT,'A','I', 1,1, ISKIP,1    ,0,0)
               IF (ISKIP /= 1) THEN
                  CALL DGEBR2D(ICTXT,'A','I', 1,1, DMAG2,1    ,0,0)
                  CALL SGEBR2D(ICTXT,'A','I', M,N, JMATL,LDJAC,0,0)
               ENDIF
            ENDIF
            IF (ISKIP == 1) GOTO 250 !break ahead, no jacobian
!
!.......... put local jacobians onto local matrices
            DO 301 I=1,M !loop on rows of Jacobian
               DO 302 J=1,N !loop on columns of Jacobian
                  VAL = DBLE(JMATL(I,J))
                  CALL PDELSET(JTLOC,J,I, DESCJT,VAL)
                  CALL PDELSET(JLOC ,I,J, DESCJ ,VAL)
  302          CONTINUE !loop on local columns (inversion points)
  301       CONTINUE !loop on rows of Jacobian (observations)
!
!.......... set covariance matrix by multiplying by diagonal weights 
            IREC = 0
            I = 1
            DO 303 IROW=1,M
               IREC = IREC + 1
               IF (IREC > rcv%NREC) THEN
                  IREC = 1
                  I = I + 1
               ENDIF
               IF (IROW == NDIM*rcv%NREC + 1) I = 1
               IF (IROW > NDIM*rcv%NREC) THEN
                  LROW = NDIM*rcv%NREC + IPERM_REC( (IREC - 1)*NDIM + I) 
               ELSE
                  LROW = IPERM_REC( (IREC - 1)*NDIM + I)
               ENDIF
               JREC = 0
               J = 1
               DO 304 JCOL=1,M
                  JREC = JREC + 1
                  IF (JREC > rcv%NREC) THEN
                     JREC = 1
                     J = J + 1
                  ENDIF
                  IF (JCOL == NDIM*rcv%NREC + 1) J = 1
                  IF (JCOL > NDIM*rcv%NREC) THEN
                     KCOL = NDIM*rcv%NREC + IPERM_REC( (JREC - 1)*NDIM + J)
                  ELSE
                     KCOL = IPERM_REC( (JREC - 1)*NDIM + J)
                  ENDIF
                  CALL INFOG2L(LROW,KCOL,DESCC, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                               IAROW,IACOL)
                  IF (IAROW == MYROW .AND. IACOL == MYCOL) THEN
                     COVD(IIA,JJA) = DBLE(SQRT(WGT(I,IFREQ,IREC,JSRC))) &
                                    *COVD_DIST(IIA,JJA) & 
                                    *DBLE(SQRT(WGT(J,IFREQ,JREC,JSRC)))
                  ENDIF
  304           CONTINUE
  303       CONTINUE
!
!.......... matrix multiply H = Re[ adj(J) C_D^{-1} J ]
            MM = M !rows of data covariance matrix  C [n_obs x n_obs]
            NM = N !number of columns of matrix B [n_obs x na35]
            KM = M !number of of columns of matrix C [n_obs x n_obs] 
            CALL PDGEMM('N','N', MM,NM,KM, ONE,COVD,1,1,DESCC, &
                        JLOC,1,1,DESCJ, ZERO, JLOC_CD,1,1,DESCJ)

!
!.......... matrix mutliply multiply H = Re[ adj(J) C_D^{-1} J ]
            MM = N !rows of adj(J) [na35 x n_{obs}] 
            NM = N !columns of matrix C_D^{-1} J [n_{obs} x na35]
            KM = M !columns of matrix adj(J) [na35 x n_{obs}]
            CALL PDGEMM('N','N', MM,NM,KM, DMAG2,JTLOC,1,1,DESCJT, &
                        JLOC_CD,1,1,DESCJ, ONE, HLOC,1,1,DESCH)
            LGN = .TRUE.
  250       CONTINUE !break ahead, no local jacobian
  205       CONTINUE !source / body wave surface wave mismatch
  200    CONTINUE !loop on sources
  100 CONTINUE !loop on frequencies
!
!.... check to make sure Re{ adj(J)J } isnt the zero matrix
      IF (.NOT.LGN) THEN
         WRITE(*,*) 'genhloc_sb: adj(J) J is the 0 matrix!',MYROW,MYCOL
         IERR = 1
      ENDIF
!
!.... free space
      DEALLOCATE(OBS)
      DEALLOCATE(WGT)
      DEALLOCATE(JTLOC)
      DEALLOCATE(IPERM_REC)
      DEALLOCATE(JLOC)
      DEALLOCATE(JLOC_CD)
      DEALLOCATE(JMATL)
      DEALLOCATE(FREQ_LOC)
 1000 CONTINUE !nothing for this process to do (would be strange)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GENHLOC(ICTXT,NPROW,NPCOL, MBJ,NBJ,MBJT,NBJT, DESCH,  &
                         LDH,M,N, LSTF, NPPGRP,FRQ,RCV,INV,SRC, HLOC, IERR)
!
!     Generates the local approximate Hessian from the saved Jacobians
!
!     INPUT      MEANING
!     -----      -------
!     ICTXT      BLACS context 
!     LDH        leading dimension for local Hessian matrix
!     LSTF       True -> convolve source time function
!     M          number of rows in Jacobian
!     MBJ        block size for jacobian
!     MBJT       block size for rows of transpose jacobian
!     NBJ        block size for jacobian
!     NBJT       block size for columns of transpose jacobian
!     N          number of columns in Jacobian 
!     NPCOL      number of process columns
!     NPROW      number of process rows
!
!     OUTPUT     MEANING
!     ------     -------
!     HLOC       local Hessian matrix
!     IERR       error flag 
! 
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (FRQ_INFO)  FRQ
      TYPE (RECV_INFO) RCV
      TYPE (SRC_INFO)  SRC
      TYPE (INV_INFO)  INV
      INTEGER*4, INTENT(IN) :: DESCH(9), ICTXT, NPROW, NPCOL, MBJ,NBJ,MBJT,NBJT, & 
                               LDH, M, N, NPPGRP
      LOGICAL*4, INTENT(IN) :: LSTF
      REAL*8, INTENT(INOUT) :: HLOC(LDH,*) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4, ALLOCATABLE :: IPERM_REC(:)
      CHARACTER(80) FLNAME 
      CHARACTER(12) CFREQ
      CHARACTER(5) CID, CSRC  
      COMPLEX*8 STF
      REAL*8 VAL, DMAG2, ONE
      INTEGER*4 IFREQ, ISRC, IPPGRP, I, J, LDJAC, MM, NM, KM, M2, ISKIP
      LOGICAL*4 LGN  
      PARAMETER(ONE = 1.D0) 
!.... BLACS variables
      INTEGER*4 DESCJT(9)  !BLACS descriptor for Jacobian^T
      INTEGER*4 DESCJ(9)   !BLACS descriptor for Jacobian 
      REAL*8, ALLOCATABLE :: JTLOC(:,:) !will hold transpose(Jacobian) 
      REAL*8, ALLOCATABLE :: JLOC(:,:)  !will holds Jacobian
      REAL*4, ALLOCATABLE :: JMATL(:,:) !single precision Jacobian for reading 
      INTEGER*4 MYROW, MYCOL, MJT, MJ, NJT, NJ, NH, INFO, NUMROC 
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &   
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &   
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 ) 
!.... debug parameters
      real*8, allocatable :: hwork(:,:)
      logical*4 ldebug 
      parameter(ldebug = .false.)
!
!----------------------------------------------------------------------------------------!
!
!.... check which processes are in the grid
      IERR = 0
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
      IF (MYROW == -1 .OR. MYCOL == -1) GOTO 1000 !nothing for this process to do 
!
!.... set sizes and space
      MJT = NUMROC(N, MBJT, MYROW, 0, NPROW) !block rows in J^T  
      MJ  = NUMROC(M, MBJ , MYROW, 0, NPROW) !block rows in J 
      NJT = NUMROC(M, NBJT, MYCOL, 0, NPCOL) !block columns in J^T
      NJ  = NUMROC(N, NBJ , MYCOL, 0, NPCOL) !block columsn in J
      CALL DESCINIT(DESCJT,N,M, MBJT,NBJT,0,0, ICTXT, MAX0(1,MJT),INFO)
      IF (INFO /= 0) THEN
         WRITE(*,*) 'genhloc: An error occurred creating DESCAT',MYROW,MYCOL
         IERR = 1
         RETURN
      ENDIF 
      CALL DESCINIT(DESCJ ,M,N, MBJ ,NBJ ,0,0, ICTXT, MAX0(1,MJ ),INFO)
      IF (INFO /= 0) THEN
         WRITE(*,*) 'genhloc: An error occurred creating DESCA',MYROW,MYCOL
         IERR = 1
         RETURN
      ENDIF 
      MM = N
      NM = N
      KM = M
      M2 = M/2 !number of obsrevations 
      ALLOCATE(JTLOC(DESCJT(LLD_),MAX0(1,NJT)))
      ALLOCATE(JLOC (DESCJ (LLD_),MAX0(1,NJ )))
      ALLOCATE(IPERM_REC(NDIM*rcv%NREC))
      LDJAC = M 
      ALLOCATE(JMATL(LDJAC,N)) 
      JTLOC(:,:) = 0.D0
      JLOC (:,:) = 0.D0 
      JMATL(:,:) = 0.0
!
!.... extract the local Jacobians 
      IF (MYROW == 0 .AND. MYCOL == 0) WRITE(*,*) 'genhloc: Assembling Hessian...'
      if (myrow == 0 .and. mycol == 0 .and. ldebug) then
         write(*,*) 'genhloc: Debug initializing hwork'
         print *, 'n=',n
         allocate(hwork(n,n))
         hwork(:,:) = 0.d0 
         print *, 't1 passed'
      endif
      NH = NUMROC(N, DESCH(NB_),MYCOL, 0,NPCOL)
      HLOC(1:DESCH(LLD_),MAX0(1,NH)) = 0.D0  !initialize to 0 b/c we add it
!
!.... loop on frequencies
      LGN = .FALSE.
      DO 100 IFREQ=1,frq%NFREQ
         IF (MYROW == 0 .AND. MYCOL == 0) THEN
            WRITE(*,*) 'genhloc: Summing frequency:',SNGL(frq%FREQ(IFREQ))
            CFREQ(1:12) = ' '
            WRITE(CFREQ,'(F12.5)') frq%FREQ(IFREQ) 
            CFREQ = ADJUSTL(CFREQ) 
         ENDIF
!
!....... loop on sources
         DO 200 ISRC=1,src%NSRC 
            IF (MYROW == 0 .AND. MYCOL == 0) THEN
               CSRC(1:5) = ' '
               WRITE(CSRC,'(I5)') ISRC
               CSRC = ADJUSTL(CSRC) 
            ENDIF
!
!.......... set peruation vector
            IF (MAXVAL(CABS(inv%OBS(1:NDIM,IFREQ,1:rcv%NREC,ISRC))) == 0.0) GOTO 250
!
!.......... generate a permutation list [w obs, v obs, u obs; null obs]
            CALL PERM_OBS(NDIM,NDIM,rcv%NREC, inv%OBS(1:NDIM,IFREQ,1:rcv%NREC,ISRC), &
                          IPERM_REC,IERR)
            IF (IERR /= 0) THEN 
               WRITE(*,*) 'genhloc_sb: Error calling perm_obs!'
               RETURN
            ENDIF
!
!.......... (0,0) process reads and broadcats
            ISKIP = 0
            IF (MYROW == 0 .AND. MYCOL == 0) THEN
!
!............. check STF is non-zero
               IF (src%SOURCE(IFREQ,ISRC) == CMPLX(0.0,0.0)) THEN
                  WRITE(*,*) 'genhloc: STF is zero, skipping!'
                  ISKIP = 1
                  GOTO 350
               ENDIF
!
!............. load the Jacobian for this source frequency pair
               JMATL(1:LDJAC,1:N) = 0.0 
               DO 300 IPPGRP=1,NPPGRP
                  CID(1:5) = ' ' 
                  WRITE(CID,'(I5)') IPPGRP - 1  
                  CID = ADJUSTL(CID) 
                  FLNAME(1:80) = ' '
                  FLNAME = './scratch/jac_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'// &
                          TRIM(CSRC)//'.dat'
                  FLNAME = ADJUSTL(FLNAME)
                  CALL LOAD_FJAC(FLNAME,LDJAC, NDIM,rcv%NREC, IPERM_REC, JMATL,IERR)
                  IF (IERR < 0) THEN
                     WRITE(*,*) 'genhloc: No local Jacobian for frequency source pair',  &
                                 frq%FREQ(IFREQ),ISRC 
                     ISKIP = 1
                     IERR = 0
                  ENDIF  
                  IF (IERR > 0) THEN 
                     WRITE(*,*) 'genhloc: Error reading file:'
                     RETURN
                  ENDIF
  300          CONTINUE 
!
!............. possibly rescale by STF  
!              IF (LSTF) CALL SCALE_JLOC(LDJAC,M2,N, src%SOURCE(IFREQ,ISRC),JMATL) 
!
!............. adj( s(w) J) s(w) J = |s|^2 adj(J) J
               IF (LSTF) THEN
                  STF = src%SOURCE(IFREQ,ISRC)
                  DMAG2 = DBLE( REAL(CONJG(STF)*STF) )
               ELSE
                  DMAG2 = 1.D0
               ENDIF
  350          CONTINUE !break ahead for no source time function 
               CALL IGEBS2D(ICTXT,'A','I', 1,1, ISKIP,1) 
               IF (ISKIP /= 1) THEN
                  CALL DGEBS2D(ICTXT,'A','I', 1,1, DMAG2,1)
                  CALL SGEBS2D(ICTXT,'A','I', M,N, JMATL,LDJAC)
               ENDIF
               if (ldebug) then 
                  write(*,*) 'genhloc: Debug multiplying jmatl jmatl'
                  do i=1,n 
                     do j=1,n
                        hwork(i,j) = hwork(i,j) + dmag2*dot_product(jmatl(:,i),jmatl(:,j))
                     enddo
                  enddo
               endif
            ELSE
               CALL IGEBR2D(ICTXT,'A','I', 1,1, ISKIP,1    ,0,0)
               IF (ISKIP /= 1) THEN
                  CALL DGEBR2D(ICTXT,'A','I', 1,1, DMAG2,1    ,0,0) 
                  CALL SGEBR2D(ICTXT,'A','I', M,N, JMATL,LDJAC,0,0)
               ENDIF
            ENDIF
            IF (ISKIP == 1) GOTO 250 !break ahead, no jacobian
!
!.......... put local jacobians onto local matrices
            DO 301 I=1,M !loop on rows of Jacobian
               DO 302 J=1,N !loop on columns of Jacobian
                  VAL = DBLE(JMATL(I,J)) 
                  CALL PDELSET(JTLOC,J,I, DESCJT,VAL)
                  CALL PDELSET(JLOC ,I,J, DESCJ ,VAL) 
  302          CONTINUE !loop on local columns (inversion points)
  301       CONTINUE !loop on rows of Jacobian (observations)
!
!.......... matrix multiply 
            CALL PDGEMM('N','N', MM,NM,KM, DMAG2,JTLOC,1,1,DESCJT, JLOC,1,1,DESCJ, &
                        ONE, HLOC,1,1,DESCH) 
            LGN = .TRUE.
  250       CONTINUE !break ahead, no local jacobian
  200    CONTINUE !loop on sources
  100 CONTINUE !loop on frequencies
!
!.... check to make sure Re{ adj(J)J } isnt the zero matrix
      IF (.NOT.LGN) THEN
         WRITE(*,*) 'genhloc: adj(J) J is the 0 matrix!',MYROW,MYCOL
         IERR = 1
      ENDIF 
!
!.... debugging on diagonal extraction 
      if (myrow == 0 .and. mycol == 0 .and. ldebug) then
         open(unit=58,file='diag_centralized.txt')
         do i=1,n 
            write(58,*) i,hwork(i,i)
            !hwork(i,i) = 2.0*hwork(i,i) 
         enddo
         write(*,*) 'checking for spd 1'
         call dpotrf('U',n,hwork,n,info)
         print *, 'info=',info
         deallocate(hwork)
         close(58)
      endif
!
!.... free space
      DEALLOCATE(IPERM_REC)
      DEALLOCATE(JTLOC)
      DEALLOCATE(JLOC) 
      DEALLOCATE(JMATL)
 1000 CONTINUE !nothing for this process to do (would be strange)
      RETURN  
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE COPY_HMAT(ICTXT,NPROW,NPCOL, LDH,M,DESCH, HLOC,HFAC) 
      REAL*8, INTENT(IN) :: HLOC(LDH,*)
      INTEGER*4, INTENT(IN) :: DESCH(9), ICTXT, NPROW, NPCOL, LDH, M 
      REAL*8, INTENT(OUT) :: HFAC(LDH,*)
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL) 
      DO 301 I=1,M !loop on rows of hessian 
         DO 302 J=1,M !loop on columns of hessian 
            CALL INFOG2L(I,J,DESCH, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                         IAROW,IACOL)
            IF (MYROW == IAROW .AND. MYCOL == IACOL) HFAC(IIA,JJA) = HLOC(IIA,JJA)
  302    CONTINUE !loop on local columns (inversion points)
  301 CONTINUE !loop on rows of Jacobian (observations)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE COPY_HMAT_UL(JOBZ,ICTXT,NPROW,NPCOL, LDH,M,DESCH, HLOC,HFAC)
      CHARACTER(1), INTENT(IN) :: JOBZ
      REAL*8, INTENT(IN) :: HLOC(LDH,*)
      INTEGER*4, INTENT(IN) :: DESCH(9), ICTXT, NPROW, NPCOL, LDH, M
      REAL*8, INTENT(OUT) :: HFAC(LDH,*)
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
      IF (JOBZ == 'L') THEN
         DO 301 I=1,M !loop on rows of hessian 
            DO 302 J=1,I !loop on columns of hessian 
               CALL INFOG2L(I,J,DESCH, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                            IAROW,IACOL)
               IF (MYROW == IAROW .AND. MYCOL == IACOL) HFAC(IIA,JJA) = HLOC(IIA,JJA)
  302       CONTINUE !loop on local columns (inversion points)
  301    CONTINUE !loop on rows of Jacobian (observations)
      ELSEIF (JOBZ == 'U') THEN
         DO 303 I=1,M !loop on rows of hessian 
            DO 304 J=I,M !loop on columns of hessian 
               CALL INFOG2L(I,J,DESCH, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                            IAROW,IACOL)
               IF (MYROW == IAROW .AND. MYCOL == IACOL) HFAC(IIA,JJA) = HLOC(IIA,JJA)
  304       CONTINUE !loop on local columns (inversion points)
  303    CONTINUE !loop on rows of Jacobian (observations)
      ELSE 
         WRITE(*,*) 'copy_hmat_ul: I dont know what you want'
         CALL COPY_HMAT(ICTXT,NPROW,NPCOL, LDH,M,DESCH, HLOC,HFAC)
      ENDIF
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GET_SOL_824(ICTXT,NPROW,NPCOL, M, DESCB, BLOC, RHS4)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: BLOC(*)
      INTEGER*4, INTENT(IN) :: DESCB(9), ICTXT, NPROW, NPCOL, M 
      REAL*4, INTENT(OUT) :: RHS4(*)
!.... local variables
      REAL*8 VAL
      INTEGER*4 I, IIA, JJA, MYROW, MYCOL, IAROW, IACOL 
!
!----------------------------------------------------------------------------------------!
! 
!.... get my rank in the BLACS grid
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL) 
!
!.... loop on rows 
      DO 1 I=1,M
         CALL INFOG2L(I,1,DESCB, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                      IAROW,IACOL) 
         IF (IAROW == 0 .AND. IACOL == 0 .AND.  &
             MYROW == 0 .AND. MYCOL == 0) THEN !on (0,0) process
            RHS4(I) = SNGL(BLOC(IIA))
         ELSE
            IF (MYROW == 0 .AND. MYCOL == 0) THEN !(0,0) receives from (myrow,mycol)
               CALL DGERV2D(ICTXT,1,1,VAL,1, IAROW,IACOL)
               RHS4(I) = SNGL(VAL) 
            ELSE !block diagonal on process (myrow,mycol)
               IF (MYROW == IAROW .AND. MYCOL == IACOL) &
               CALL DGESD2D(ICTXT,1,1,BLOC(IIA),1, 0,0) 
            ENDIF
         ENDIF
    1 CONTINUE !loop on rows
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GET_SOL8(ICTXT,NPROW,NPCOL, M, DESCB, BLOC, RHS)
      IMPLICIT NONE 
      REAL*8, INTENT(IN) :: BLOC(*)
      INTEGER*4, INTENT(IN) :: DESCB(9), ICTXT, NPROW, NPCOL, M 
      REAL*8, INTENT(OUT) :: RHS(*)
!.... local variables
      INTEGER*4 I, IIA, JJA, MYROW, MYCOL, IAROW, IACOL 
!
!----------------------------------------------------------------------------------------!
! 
!.... get my rank in the BLACS grid
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL) 
!
!.... loop on rows 
      DO 1 I=1,M
         CALL INFOG2L(I,1,DESCB, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                      IAROW,IACOL) 
         IF (IAROW == 0 .AND. IACOL == 0 .AND.  &
             MYROW == 0 .AND. MYCOL == 0) THEN !on (0,0) process
            RHS(I) = BLOC(IIA)
         ELSE 
            IF (MYROW == 0 .AND. MYCOL == 0) THEN !(0,0) receives from (myrow,mycol)
               CALL DGERV2D(ICTXT,1,1,RHS(I),1, IAROW,IACOL)
            ELSE !block diagonal on process (myrow,mycol)
               IF (MYROW == IAROW .AND. MYCOL == IACOL) &
               CALL DGESD2D(ICTXT,1,1,BLOC(IIA),1, 0,0) 
            ENDIF
         ENDIF
    1 CONTINUE !loop on rows
      RETURN
      END  
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE SET_SOL8(ICTXT,NPROW,NPCOL, M, DESCB, B, BLOC)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: B(*) 
      INTEGER*4, INTENT(IN) :: DESCB(9), ICTXT, NPROW, NPCOL, M
      REAL*8, INTENT(OUT) :: BLOC(*)
!.... local variables
      INTEGER*4 I, IIA, JJA, MYROW, MYCOL, IAROW, IACOL
!
!----------------------------------------------------------------------------------------!
! 
!.... get my rank in the BLACS grid
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!
!.... loop on rows 
      DO 1 I=1,M
         CALL INFOG2L(I,1,DESCB, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                      IAROW,IACOL)
         IF (MYROW == 0 .AND. MYCOL == 0) THEN 
            IF (IAROW == 0 .AND. IACOL == 0) THEN !mine
               BLOC(IIA) = B(I)
            ELSE !(0,0) sends to (iarow,iacol) 
               CALL DGESD2D(ICTXT,1,1,B(I),1, IAROW,IACOL)
            ENDIF 
         ELSE !(myrow,mycol) receivers diagonal
            IF (MYROW == IAROW .AND. MYCOL == IACOL) &
            CALL DGERV2D(ICTXT,1,1,BLOC(IIA),1, 0,0)
         ENDIF
    1 CONTINUE !loop on rows
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE DSCALE_HLOC(ICTXT,NPROW,NPCOL,LINV, LDH,M, DESCH,DIAG, HLOC) 
!
!     Scales the matrix hloc by the diagonal matrix D s.t. B <- D B D
!
!.... variable declarations
      IMPLICIT NONE 
      REAL*8, INTENT(INOUT) :: HLOC(LDH,*)
      REAL*8, INTENT(IN) :: DIAG(M) 
      INTEGER*4, INTENT(IN) :: DESCH(9), LDH, M, NPROW, NPCOL, ICTXT 
      LOGICAL*4 LINV 
!.... local variables
      REAL*8 RMIN, D1, D2
      INTEGER*4 MYROW, MYCOL, I,J, IIA,JJA, IAROW, IACOL 
!----------------------------------------------------------------------------------------!
!
!.... get my rank in the BLACS grid
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!
!.... warning 
      IF (LINV) THEN
         IF (MINVAL(DIAG) == 0.D0) THEN
            WRITE(*,*) 'dscale_hloc: Warning minval(diag) == 0!'
            RMIN = HUGE(1.D0) 
            DO 1 I=1,M
               IF (DABS(DIAG(I)) > 0.D0) RMIN = DMIN1(DABS(DIAG(I)),RMIN)
    1       CONTINUE
            WRITE(*,*) 'dscale_hloc: rmin is ',RMIN
            RMIN = 1.D0/RMIN
         ENDIF
       ENDIF
              
!
!.... scale by diagonal matrices
      DO 3 I=1,M
         DO 4 J=1,M
            CALL INFOG2L(I,J,DESCH, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                         IAROW,IACOL)
            IF (IAROW == MYROW .AND. IACOL == MYCOL) THEN 
               IF (LINV) THEN
                  IF (DIAG(I) > 0.D0) THEN
                     D1 = 1.D0/DIAG(I)
                  ELSE
                     D1 = RMIN
                  ENDIF
                  IF (DIAG(J) > 0.D0) THEN
                     D2 = 1.D0/DIAG(J)
                  ELSE
                     D2 = RMIN 
                  ENDIF
                  HLOC(IIA,JJA) = D1*HLOC(IIA,JJA)*D2
               ELSE
                  HLOC(IIA,JJA) = DIAG(I)*HLOC(IIA,JJA)*DIAG(J)
               ENDIF
            ENDIF
    4    CONTINUE
    3 CONTINUE
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GET_DIAG(ICTXT,NPROW,NPCOL, LDH,M, LBCAST, DESCH,HLOC, DIAG,IERR)
!
!     Gets the diagonal of the global Hessian from the local hessian
!
!     INPUT      MEANING
!     -----      ------- 
!     DESCH      BLACS descriptor for matrix H
!     ICTXT      BLACS grid context 
!     LBCAST     True -> send diagonal to all other processes in grid
!     LDH        leading dimension of HLOC (probably DESCH(9))
!     M          number of rows/columns in matrix
!     NPCOL      number of process columns
!     NPROW      number of process rows
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     DIAG       block diagonal (safeguarded so all elements > 0'
!     IERR       error flag, diagonal < 0.  error in matrix multiplication
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: HLOC(LDH,*)
      INTEGER*4, INTENT(IN) :: DESCH(9), LDH, M, NPROW, NPCOL, ICTXT 
      LOGICAL*4, INTENT(IN) :: LBCAST 
      REAL*8, INTENT(OUT) :: DIAG(*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8 DMIN
      INTEGER*4 MYROW, MYCOL, I, IIA, JJA, IAROW, IACOL 
      LOGICAL*4 LDEBUG
      PARAMETER(LDEBUG = .FALSE.)
!----------------------------------------------------------------------------------------!
!
!.... get my rank in the BLACS grid
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL) 
      IF (MYROW == 0 .AND. MYCOL == 0) DMIN = HUGE(1.D0)
!
!.... debugging
      IF (MYROW == 0 .AND. MYCOL == 0 .AND. LDEBUG) OPEN(UNIT=57,FILE='diag_dist.txt')
!
!.... loop on diagonal
      DO 1 I=1,M
         CALL INFOG2L(I,I,DESCH, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                      IAROW,IACOL) 
         IF (IAROW == 0 .AND. IACOL == 0 .AND.  &
             MYROW == 0 .AND. MYCOL == 0) THEN !on (0,0) process
            DIAG(I) = HLOC(IIA,JJA)
         ELSE
            IF (MYROW == 0 .AND. MYCOL == 0) THEN !(0,0) receives from (myrow,mycol)
               CALL DGERV2D(ICTXT,1,1,DIAG(I),1, IAROW,IACOL)
            ELSE !block diagonal on process (myrow,mycol)
               IF (MYROW == IAROW .AND. MYCOL == IACOL) &
               CALL DGESD2D(ICTXT,1,1,HLOC(IIA,JJA),1, 0,0)
            ENDIF
         ENDIF
         IF (MYROW == 0 .AND. MYCOL == 0 .AND. LDEBUG) WRITE(57,*) I,DIAG(I) 
         IF (MYROW == 0 .AND. MYCOL == 0) THEN
            IF (DIAG(I) > 0.D0) DMIN = DMIN1(DMIN,DIAG(I))
         ENDIF
    1 CONTINUE
      IF (LBCAST) THEN
         IERR = 0
         IF (MYROW == 0 .AND. MYCOL == 0) THEN
            DO 2 I=1,M 
               IF (DIAG(I) <= 0.D0) THEN
                  WRITE(*,*) 'get_diag: Serious error diag <= 0'
                  IERR = 1
               ENDIF
               IF (DIAG(I) == 0.D0) DIAG(I) = DMIN
    2       CONTINUE  
            CALL DGEBS2D(ICTXT,'A','I',M,1,DIAG,M)
         ELSE 
            CALL DGEBR2D(ICTXT,'A','I',M,1,DIAG,M, 0,0) 
         ENDIF
      ENDIF
      IF (MYROW == 0 .AND. MYCOL == 0 .AND. LDEBUG) CLOSE(57)
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE ADD_DIAG(ICTXT,NPROW,NPCOL, LDH,M, DESCH, RLAM,DIAG, HLOC) 
!
!     Adds rlam*diagonal to local Hessian hloc
!
!     INPUT      MEANING
!     -----      ------- 
!     DIAG       diagonal
!     DESCH      BLACS descriptor for matrix H
!     ICTXT      BLACS grid context 
!     LDH        leading dimension of HLOC (probably DESCH(9))
!     M          number of rows/columns in matrix
!     NPCOL      number of process columns
!     NPROW      number of process rows
!     RLAM       amount to add to diagonal
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     DIAG       block diagonal on (0,0) process 
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: HLOC(LDH,*)
      REAL*8, INTENT(IN) :: DIAG(M), RLAM 
      INTEGER*4, INTENT(IN) :: DESCH(9), LDH, M, NPROW, NPCOL, ICTXT 
!.... local variables
      INTEGER*4 MYROW, MYCOL, I, IIA, JJA, IAROW, IACOL 
!
!----------------------------------------------------------------------------------------!
!
!.... get my rank in the BLACS grid
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!
!.... add the scaled diagonal back in 
      DO 1 I=1,M
         CALL INFOG2L(I,I,DESCH, NPROW,NPCOL, MYROW,MYCOL, IIA,JJA, &
                      IAROW,IACOL)
         IF (MYROW == IAROW .AND. MYCOL == IACOL) &
         HLOC(IIA,JJA) = HLOC(IIA,JJA) + RLAM*DIAG(I)
    1 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE LOAD_FJAC(FLNAME,LDJAC, NDIM,NREC, IPERM_REC, JMATL,IERR) 
!
!     Loads the full Jacobians, and optionally whether the observation exists at a 
!     receiver/component pair, and the destination in the global inversion problem
!
!.... variable declarations
      CHARACTER(*), INTENT(IN) :: FLNAME
      INTEGER*4, INTENT(IN) :: IPERM_REC(NDIM*NREC), NDIM, NREC
      REAL*4 JMATL(LDJAC,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*4, ALLOCATABLE :: RJAC(:,:), CJAC(:,:) 
      INTEGER*4, ALLOCATABLE :: MYDEST(:)
      LOGICAL*4, ALLOCATABLE :: LOBS(:,:) 
      INTEGER*4 NREC_IN, NDIM_IN, I, J, K, L, IREC, IOBS, INPGL, IUNIT  
      LOGICAL*4 LEX 
      PARAMETER(IUNIT = 65) 
!
!----------------------------------------------------------------------------------------!
!
!.... check the file exists
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'load_fjac: No local Jacobian file found'
         IERR =-1
         RETURN
      ENDIF 
!.... open file then read
      IERR = 0  
      OPEN(FILE=TRIM(FLNAME),UNIT=IUNIT,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'load_fjac: Error opening file:',TRIM(FLNAME)
         RETURN
      ENDIF
      READ(IUNIT,IOSTAT=IERR) NOBS_IN,NNPGL_IN, NDIM_IN, NREC_IN
      IF (IERR /= 0) THEN 
         WRITE(*,*) 'load_fjac: Error reading header'
         CLOSE(IUNIT)
         RETURN
      ENDIF
      ALLOCATE(RJAC(NOBS_IN,NNPGL_IN))
      ALLOCATE(CJAC(NOBS_IN,NNPGL_IN)) 
      ALLOCATE(MYDEST(NNPGL_IN))
      ALLOCATE(LOBS(NDIM,NREC_IN)) 
      IF (NDIM_IN /= NDIM) THEN 
         WRITE(*,*) 'load_fjac: Serious error, ndim_in /= ndim',NDIM_IN,NDIM
         IERR = 1
         CLOSE(IUNIT)
         RETURN
      ENDIF
      IF (NREC_IN /= NREC) THEN 
         WRITE(*,*) 'load_fjac: Serious error, nrec_in /= nrec',NREC_IN,NREC
         IERR = 1
         CLOSE(IUNIT)
         RETURN
      ENDIF
      READ(IUNIT,IOSTAT=IERR) ((LOBS(I,IREC),I=1,NDIM),IREC=1,NREC_IN)
      IF (IERR /= 0) THEN 
         WRITE(*,*) 'load_fjac: Error reading LOBS'
         CLOSE(IUNIT)
         RETURN
      ENDIF
      READ(IUNIT,IOSTAT=IERR)((RJAC(IOBS,INPGL),IOBS=1,NOBS_IN),INPGL=1,NNPGL_IN)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'load_fjac: Error reading real Jacobian!'
         CLOSE(IUNIT)
         RETURN
      ENDIF
      READ(IUNIT,IOSTAT=IERR)((CJAC(IOBS,INPGL),IOBS=1,NOBS_IN),INPGL=1,NNPGL_IN)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'load_fjac: Error reading complex Jacobian!'
         CLOSE(IUNIT)
         RETURN
      ENDIF
      MYDEST(1:NNPGL_IN) = 0 
      READ(IUNIT,IOSTAT=IERR)(MYDEST(INPGL),INPGL=1,NNPGL_IN)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'load_fjac: Error reading inversion destinations!'
         CLOSE(IUNIT)
         RETURN
      ENDIF
      CLOSE(IUNIT) !done with file 
!
!.... put into jacobian
      I = 0 
      K = NDIM*NREC
      L = 0
      DO 1 IREC=1,NREC_IN
         DO 2 J=1,NDIM
            INDX = (IREC - 1)*NDIM + J
            I = IPERM_REC(INDX)
            K = NDIM*NREC + I
            L = L + 1
            DO 3 INPGL=1,NNPGL_IN
               INDX = MYDEST(INPGL)
               IF (LOBS(J,IREC)) THEN
                  JMATL(I,INDX) = RJAC(L,INPGL)
                  JMATL(K,INDX) = CJAC(L,INPGL)
               ELSE
                  JMATL(I,INDX) = 0.0
                  JMATL(K,INDX) = 0.0
               ENDIF
    3       CONTINUE 
    2    CONTINUE 
    1 CONTINUE 
      DEALLOCATE(CJAC)
      DEALLOCATE(RJAC) 
      DEALLOCATE(MYDEST) 
      DEALLOCATE(LOBS) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE PICARD_DIST(ICTXT,NPROW,NPCOL, LDH, DESCH,HLOC,G, IERR) 
!
!     Generates a Picard plot.  I'm out of ideas for surface waves. 
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: HLOC(LDH,*), G(*)
      INTEGER*4, INTENT(IN) :: DESCH(9), ICTXT, NPROW, NPCOL, LDH 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: A(:,:), Z(:,:), WORK(:), GLOC(:), YLOC(:), SOL(:), W(:)
      REAL*8 WORK8
      INTEGER*4 DESCA(9), DESCZ(9), DESCG(9), I , MYROW, MYCOL, LDA, LWORK
      INTEGER*4 NUMROC     !calculates NUmber of Rows Or Columns
      INTEGER*4 INFO       !error flag from ScaLapack 
      INTEGER*4 N,M        !number of rows and columns
      INTEGER*4 MH         !number of local rows in h
      INTEGER*4 NH         !number of local columns in h
      INTEGER*4 BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER(BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,    &    
                CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,   &
                RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
!
!----------------------------------------------------------------------------------------!
!
!.... set A 
      IERR = 0
      LDA = DESCH(LLD_)
      M = DESCH(M_)
      N = DESCH(N_)
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
      MH = NUMROC(M, DESCH(MB_),MYROW, 0,NPROW) !number of rows in local H matrix
      NH = NUMROC(N, DESCH(NB_),MYCOL, 0,NPCOL) !number of columns in local H matrix
      DESCA(1:9) = DESCH(1:9)
      ALLOCATE(A(LDA,MAX(NH,1)),STAT=IERR)
      A(:,:) = 0.D0 
      IF (IERR /= 0) THEN 
         WRITE(*,*) 'picard_dist: Error setting space for A on grid',MYROW,MYCOL
         RETURN
      ENDIF
      A(:,:) = 0.D0
      ALLOCATE(Z(LDA,MAX(NH,1)),STAT=IERR)
      IF (IERR /= 0) THEN 
         WRITE(*,*) 'picard_dist: Error setting space for Z on grid',MYROW,MYCOL
         RETURN
      ENDIF
      Z(:,:) = 0.D0
      DESCZ(1:9) = DESCH(1:9)
      CALL COPY_HMAT_UL('U',ICTXT,NPROW,NPCOL, LDH,N,DESCH, HLOC,A)
!
!.... calculate the eigenvalue decomposition
      IF (MYROW == 0 .AND. MYCOL == 0) &
      WRITE(*,*) 'picard_dist: Calculating eigendecomposition...'
      LWORK =-1
      ALLOCATE(W(N))
      CALL PDSYEV('V','U',N,A,1,1,DESCA, W,  &
                  Z,1,1, DESCZ, WORK8,LWORK, INFO)
      LWORK = INT(WORK8)
      ALLOCATE(WORK(LWORK),STAT=IERR)
      IF (IERR /= 0) THEN
         IF (MYROW == 0 .AND. MYCOL == 0) &
         WRITE(*,*) 'picard_dist: Error allocating space for work',IERR
         RETURN
      ENDIF
      CALL PDSYEV('V','U',N,A,1,1,DESCA, W,  &
                  Z,1,1, DESCZ, WORK ,LWORK, INFO)
      IF (INFO /= 0) THEN
         IF (MYROW == 0 .AND. MYCOL == 0) &
         WRITE(*,*) 'picard_dist: Error calling pdysev',INFO
         IERR = 1
         RETURN
      ENDIF
      DEALLOCATE(WORK) !done with workspace 

!
!.... dot a row of U with g; all rows is matrix vecotr multiplication  
      CALL DESCINIT(DESCG,M,1, DESCA(MB_),1, 0,0, ICTXT, MAX(1,MH),INFO)
      ALLOCATE(GLOC(MAX(1,MH)))
      ALLOCATE(YLOC(MAX(1,MH)))
      CALL SET_SOL8(ICTXT,NPROW,NPCOL, M, DESCG, G, GLOC) !g -> local g
      IF (MYROW == 0 .AND. MYCOL == 0) THEN
         ALLOCATE(SOL(N))
      ELSE
         ALLOCATE(SOL(1))
         SOL(1) = 0.D0
      ENDIF
      CALL PDGEMV('T',N,N,1.D0,Z,1,1, DESCZ, GLOC,1,1,DESCG,1, 0.D0,YLOC,1,1,DESCG,1)
      CALL GET_SOL8(ICTXT,NPROW,NPCOL, N, DESCG, YLOC, SOL) !send g back to (0,0)
      DEALLOCATE(Z)
      DEALLOCATE(A) 
      DEALLOCATE(GLOC)
      DEALLOCATE(YLOC)
!
!.... output picard plot
      IF (MYROW == 0 .AND. MYCOL == 0) THEN 
         OPEN(UNIT=24,FILE='picard.txt')
         DO 100 I=M,1,-1 !loop on rows; eigendecomposition is backwards
            WRITE(24,900) M + 1 - I, W(I), ABS(SOL(I)), ABS(SOL(I))/W(I)
  900       FORMAT(I6,G12.4,G12.4,G12.4)
  100    CONTINUE
         CLOSE(24) 
      ENDIF
      DEALLOCATE(SOL)
      DEALLOCATE(W)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GEN_PGRID(M,N,NPROCS, NPROW,NPCOL,IERR)
!
!     Generates a square-ish process grid nprocs ~= nprow*npcol.  
!     At most we can have 8x8 grids.  This may change
!
!     INPUT      MEANING
!     -----      ------- 
!     M          number of rows in matrix
!     N          number of columsn in matrix
!     NPROCS     number of processes available 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     NPCOL      number of process columns
!     NPROW      number of process rows
!    
!.... variable declarations
      IMPLICIT NONE
      INTEGER*4, INTENT(IN) :: M,N,NPROCS
      INTEGER*4, INTENT(OUT) :: NPROW, NPCOL, IERR
      INTEGER*4 NPROW0, NPCOL0, NRJ, NCJ, I,J
!
!----------------------------------------------------------------------------------------!
! 
      IERR = 0
      NPROW = 1
      NPCOL = 1
      IF (NPROCS >= 64) THEN
         WRITE(*,*) 'gen_pgrid: Warning overriding to 6x6 grid'
         NPROW = 8
         NPCOL = 8
         GOTO 500
      ENDIF
      NRJ = M
      NCJ = N
!
!.... try to get squarish block sizes 
      NPROW0 = 0 
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
    2    CONTINUE
         NPROW = NPROW + 1
    1 CONTINUE
  500 CONTINUE !break ahead
      IF (NPROW*NPCOL > NPROCS) THEN
         WRITE(*,*) 'gen_pgrid: Serious error nprow*npcol > nprocs!'
         IERR = 1
         RETURN
      ENDIF
      IF (NPROW == 0 .OR. NPCOL == 0) THEN
         WRITE(*,*) 'gen_pgrid: Serious error process grid is size 0!'
         IERR = 2
         RETURN
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GEN_PGRID_SQR(M,N,NPROCS, NPROW,NPCOL,IERR)
!
!     Generates a square process grid nprocs ~= nprow*npcol.  
!     At most we can have 8x8 grids.  This may change
!
!     INPUT      MEANING
!     -----      ------- 
!     M          number of rows in matrix
!     N          number of columsn in matrix
!     NPROCS     number of processes available 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     NPCOL      number of process columns
!     NPROW      number of process rows
!    
!.... variable declarations
      IMPLICIT NONE 
      INTEGER*4, INTENT(IN) :: M,N,NPROCS
      INTEGER*4, INTENT(OUT) :: NPROW, NPCOL, IERR 
      INTEGER*4 NPROW0, NPCOL0, NRJ, NCJ, I,J
!
!----------------------------------------------------------------------------------------!
! 
      IERR = 0
      NPROW = 1
      NPCOL = 1
      IF (NPROCS >= 64) THEN 
         WRITE(*,*) 'gen_pgrid: Warning overriding to 6x6 grid'
         NPROW = 8
         NPCOL = 8
         GOTO 500
      ENDIF
      NRJ = M
      NCJ = N
!
!.... try to get square block sizes 
      NPROW0 = 0
      NPCOL0 = 1
      NPROW = 1
      DO 1 I=1,8
         NPCOL = NPROW
         IF (NPROW == NRJ .OR. NPCOL == NCJ) GOTO 500 !this would be bad
         DO 2 J=1,1
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
    2    CONTINUE
         NPROW = NPROW + 1
    1 CONTINUE
  500 CONTINUE !break ahead
      IF (NPROW*NPCOL > NPROCS) THEN
         WRITE(*,*) 'gen_pgrid: Serious error nprow*npcol > nprocs!'
         IERR = 1
         RETURN
      ENDIF
      IF (NPROW == 0 .OR. NPCOL == 0) THEN
         WRITE(*,*) 'gen_pgrid: Serious error process grid is size 0!'
         IERR = 2
         RETURN
      ENDIF
      RETURN
      END

