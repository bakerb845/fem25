      SUBROUTINE INVHESS(NHSIZE,NA35, NNPINV,NVINV, GRADPC,IPIVH,IERR) 
!
!     Inverts the Hessian precoditioner.  Note that hess is stored in 
!     CRS format but since each block is a submatrix we don't bother 
!     making the other pointers 
!
!     INPUT      MEANING
!     -----      ------- 
!     GRADPC     non-inverted gradient pre-conditioner
!     NA35       size of problem = NNPINV*NVINV 
!     NNPINV     number of nodal points in inversion
!     NHSIZE     number of non-zeros in block-diagonal Hessian
!     NVINV      number of variables to invert for
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag 
!     IPIVH      pivots on factorization of block diagonal Hessian
!     GRADPC     LU factored gradient preconditioner
!  
!.... variable declarations
      REAL*4, INTENT(INOUT) :: GRADPC(NHSIZE)  
      INTEGER*4, INTENT(IN) :: NHSIZE, NA35, NNPINV, NVINV 
      INTEGER*4, INTENT(OUT) :: IPIVH(NA35), IERR  
!.... local variables
      REAL*4, ALLOCATABLE :: BDIAG(:,:) 
      PARAMETER(EPS = 1.19209290E-07)
!
!----------------------------------------------------------------------!
!
!.... this matrix is diagonal
      IERR = 0 
      WRITE(*,*) 'invhess: Inverting Hessian preconditioner...'
      IF (NVINV.EQ.1) THEN
         DO 1 IA35=1,NA35
            IF (GRADPC(IA35).EQ.0.0) THEN
               WRITE(*,*) 'invhess: Warning division by zero!'
               GRADPC(IA35) = EPS
            ELSEIF (GRADPC(IA35).LT.0.0) THEN
               WRITE(*,*) 'invhess: Warning Hessian is less than 0!'
               GRADPC(IA35) = EPS
            ENDIF
            GRADPC(IA35) = 1.0/GRADPC(IA35) 
            IPIVH(IA35) = 1 
    1    CONTINUE
      ELSE
         ALLOCATE(BDIAG(NVINV,NVINV))  
         INDX = 0 
         JNDX = 0
         DO 2 INPINV=1,NNPINV !loop on inversion point nodes 
            DO 3 IVINV=1,NVINV !extract a block 
               DO 4 JVINV=1,NVINV 
                  INDX = INDX + 1 
                  BDIAG(IVINV,JVINV) = GRADPC(INDX)
    4          CONTINUE !loop on columns
    3       CONTINUE !loop on rows of submatrix 
            I1 = (INPINV - 1)*NVINV + 1 
            I2 = INPINV*NVINV
            CALL SGETRF(NVINV,NVINV,BDIAG,NVINV,IPIVH(I1:I2),INFO) 
            IF (INFO.LT.0) THEN 
               WRITE(*,*) 'invhess: Error in cgetrf input:',INFO
               IERR = 1
               GOTO 750
            ENDIF
            IF (INFO.GT.0) THEN
               WRITE(*,*) 'invhess: Hessian block is singular!',INFO
               IERR = 2
               GOTO 750
            ENDIF
!
!.......... copy the inverted hessian back
            DO 5 IVINV=1,NVINV
               DO 6 JVINV=1,NVINV
                  JNDX = JNDX + 1
                  GRADPC(JNDX) = BDIAG(IVINV,JVINV)
    6          CONTINUE
    5       CONTINUE  
    2    CONTINUE !loop on anchor nodes in inversion
         DEALLOCATE(BDIAG) 
      ENDIF
  750 CONTINUE !break ahead for an error
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DID_BDIAG(NHSIZE,NA35, NNPINV,NVINV, IPIV,D)
!     Sets block diagonal to unity
      IMPLICIT NONE
      INTEGER*4, INTENT(IN) :: NNPINV, NVINV, NHSIZE, NA35
      REAL*8, INTENT(OUT) :: D(NHSIZE)
      INTEGER*4, INTENT(OUT) :: IPIV(NA35)
!.... local variables
      INTEGER*4 IHSIZE, IA35, IVINV, JVINV, INPINV 
!
!----------------------------------------------------------------------!
!
      IHSIZE = 0 
      IA35 = 0 
      DO 1 INPINV=1,NNPINV 
         DO 2 IVINV=1,NVINV
            IA35 = IA35 + 1
            IPIV(IA35) = 1
            DO 3 JVINV=1,NVINV
               IHSIZE = IHSIZE + 1 
               IF (IVINV.EQ.JVINV) THEN
                  D(IHSIZE) = 1.D0
               ELSE
                  D(IHSIZE) = 0.D0
               ENDIF
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DINVHESS(NHSIZE,NA35, NNPINV,NVINV, GRADPC,IPIVH,IERR) 
!
!     Inverts the Hessian precoditioner.  Note that hess is stored in 
!     CRS format but since each block is a submatrix we don't bother 
!     making the other pointers.  Double precision 
!
!     INPUT      MEANING
!     -----      ------- 
!     GRADPC     non-inverted gradient pre-conditioner
!     NA35       size of problem = NNPINV*NVINV 
!     NNPINV     number of nodal points in inversion
!     NHSIZE     number of non-zeros in block-diagonal Hessian
!     NVINV      number of variables to invert for
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag 
!     IPIVH      pivots on factorization of block diagonal Hessian
!     GRADPC     LU factored gradient preconditioner
!  
!.... variable declarations
      REAL*8, INTENT(INOUT) :: GRADPC(NHSIZE)  
      INTEGER*4, INTENT(IN) :: NHSIZE, NA35, NNPINV, NVINV 
      INTEGER*4, INTENT(OUT) :: IPIVH(NA35), IERR  
!.... local variables
      REAL*8, ALLOCATABLE :: BDIAG(:,:) 
      REAL*8 EPS 
      PARAMETER(EPS = 4.440892098500626D-016)
!
!----------------------------------------------------------------------!
!
!.... this matrix is diagonal
      IERR = 0
      WRITE(*,*) 'invhess: Inverting Hessian preconditioner...'
      IF (NVINV.EQ.1) THEN
         DO 1 IA35=1,NA35
            IF (GRADPC(IA35).EQ.0.D0) THEN
               WRITE(*,*) 'dinvhess: Warning division by zero!'
               GRADPC(IA35) = EPS
            ELSEIF (GRADPC(IA35).LT.0.D0) THEN
               WRITE(*,*) 'dinvhess: Warning Hessian is less than 0!'
               GRADPC(IA35) = EPS
            ENDIF
            GRADPC(IA35) = 1.D0/GRADPC(IA35)
            IPIVH(IA35) = 1
    1    CONTINUE
      ELSE
         ALLOCATE(BDIAG(NVINV,NVINV))
         INDX = 0
         JNDX = 0
         DO 2 INPINV=1,NNPINV !loop on inversion point nodes 
            DO 3 IVINV=1,NVINV !extract a block 
               DO 4 JVINV=1,NVINV
                  INDX = INDX + 1
                  BDIAG(IVINV,JVINV) = GRADPC(INDX)
    4          CONTINUE !loop on columns
    3       CONTINUE !loop on rows of submatrix 
            I1 = (INPINV - 1)*NVINV + 1
            I2 = INPINV*NVINV
            CALL DGETRF(NVINV,NVINV,BDIAG,NVINV,IPIVH(I1:I2),INFO)
            IF (INFO.LT.0) THEN
               WRITE(*,*) 'dinvhess: Error in cgetrf input:',INFO
               IERR = 1
               GOTO 750
            ENDIF
            IF (INFO.GT.0) THEN
               WRITE(*,*) 'dinvhess: Hessian block is singular!',INFO
               IERR = 2
               GOTO 750
            ENDIF
!
!.......... copy the inverted hessian back
            DO 5 IVINV=1,NVINV
               DO 6 JVINV=1,NVINV
                  JNDX = JNDX + 1
                  GRADPC(JNDX) = BDIAG(IVINV,JVINV)
    6          CONTINUE
    5       CONTINUE
    2    CONTINUE !loop on anchor nodes in inversion
         DEALLOCATE(BDIAG)
      ENDIF
  750 CONTINUE !break ahead for an error
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE PREGRAD(NHSIZE,NA35, NNPINV,NVINV, IPIVH,GRADPC,GRAD,
     ;                   PGRAD,IERR)
!
!     Applies the inverse of the block diagonal Hessian to the gradient 
!
!     INPUT      MEANING
!     -----      ------- 
!     IPIVH      holds the pivots on the pre-conditioner 
!     GRAD       gradient (not-preconditioner)
!     GRADPC     block diagonal LU factored gradient preconditioner
!     NA35       number of points in inversion = NVINV*NNPINV
!     NHSIZE     number of elements in GRADPC 
!     NNPINV     number of nodal points in inversion
!     NVINV      number of variables in inversion
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     PGRAD      pre-conditioned gradient
!
!.... variable declarations
      REAL*4, INTENT(IN) :: GRADPC(NHSIZE), GRAD(NA35)  
      INTEGER*4, INTENT(IN) :: IPIVH(NA35), NHSIZE, NA35, NNPINV, NVINV 
      REAL*4, INTENT(OUT) :: PGRAD(NA35) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*4, ALLOCATABLE :: BDIAG(:,:)
!
!----------------------------------------------------------------------!
!
!.... this matrix is diagonal
      IERR = 0 
      IF (NVINV.EQ.1) THEN 
         DO 1 IA35=1,NA35
            PGRAD(IA35) = GRADPC(IA35)*GRAD(IA35)
    1    CONTINUE 
      ELSE !this matrix is block diagonal
         CALL SCOPY(NA35,GRAD,1,PGRAD,1) 
         ALLOCATE(BDIAG(NVINV,NVINV))  
         INDX = 0 
         DO 2 INPINV=1,NNPINV !loop on inversion point nodes 
            DO 3 IVINV=1,NVINV !extract a block 
               DO 4 JVINV=1,NVINV 
                  INDX = INDX + 1 
                  BDIAG(IVINV,JVINV) = GRADPC(INDX)
    4          CONTINUE !loop on columns
    3       CONTINUE !loop on rows of submatrix 
            I1 = (INPINV - 1)*NVINV + 1
            I2 = INPINV*NVINV 
            CALL SGETRS('N',NVINV,1,BDIAG,NVINV,IPIVH(I1:I2), 
     ;                  PGRAD(I1:I2),NVINV,INFO)
            IF (INFO.LT.0) THEN
               WRITE(*,*) 'pregrad: Error in argument:',INFO
               IERR = 1 
               RETURN
            ENDIF
            I1 = (INPINV - 1)*NVINV + 1
            I2 = INPINV*NVINV 
    2    CONTINUE 
         DEALLOCATE(BDIAG)
      ENDIF !end check on type of matrix
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DPREGRAD(NHSIZE,NA35, NNPINV,NVINV, IPIVH,GRADPC,GRAD,
     ;                    PGRAD,IERR)
!
!     Applies the inverse of the block diagonal Hessian to the gradient
!     This is the double precision version.   
!
!     INPUT      MEANING
!     -----      ------- 
!     IPIVH      holds the pivots on the pre-conditioner 
!     GRAD       gradient (not-preconditioner)
!     GRADPC     block diagonal LU factored gradient preconditioner
!     NA35       number of points in inversion = NVINV*NNPINV
!     NHSIZE     number of elements in GRADPC 
!     NNPINV     number of nodal points in inversion
!     NVINV      number of variables in inversion
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     PGRAD      pre-conditioned gradient
!
!.... variable declarations
      REAL*8, INTENT(IN) :: GRADPC(NHSIZE), GRAD(NA35)  
      INTEGER*4, INTENT(IN) :: IPIVH(NA35), NHSIZE, NA35, NNPINV, NVINV 
      REAL*8, INTENT(OUT) :: PGRAD(NA35) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: BDIAG(:,:)
!
!----------------------------------------------------------------------!
!
!.... this matrix is diagonal
      IERR = 0 
      IF (NVINV.EQ.1) THEN 
         DO 1 IA35=1,NA35
            PGRAD(IA35) = GRADPC(IA35)*GRAD(IA35)
    1    CONTINUE 
      ELSE !this matrix is block diagonal
         CALL DCOPY(NA35,GRAD,1,PGRAD,1)
         ALLOCATE(BDIAG(NVINV,NVINV))
         INDX = 0
         DO 2 INPINV=1,NNPINV !loop on inversion point nodes 
            DO 3 IVINV=1,NVINV !extract a block 
               DO 4 JVINV=1,NVINV
                  INDX = INDX + 1
                  BDIAG(IVINV,JVINV) = GRADPC(INDX)
    4          CONTINUE !loop on columns
    3       CONTINUE !loop on rows of submatrix 
            I1 = (INPINV - 1)*NVINV + 1
            I2 = INPINV*NVINV
            CALL DGETRS('N',NVINV,1,BDIAG,NVINV,IPIVH(I1:I2),
     ;                  PGRAD(I1:I2),NVINV,INFO)
            IF (INFO.LT.0) THEN
               WRITE(*,*) 'dpregrad: Error in argument:',INFO
               IERR = 1
               RETURN
            ENDIF
            I1 = (INPINV - 1)*NVINV + 1
            I2 = INPINV*NVINV
    2    CONTINUE
         DEALLOCATE(BDIAG)
      ENDIF !end check on type of matrix
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE REGPRECON(NHSIZE,NA35, NNPINV,NVINV, PTHRESH,
     ;                     BDHESS, GRADPC)
!
!     Regularizes the pre-conditioner by thresholding the smallest 
!     singular values to the bottom thresh.  On input we have the 
!     complex valued block diagonal adj(J)J.  So we reguarlize this 
!     and output the real since the inversion is defined by 
!       Re{ adj(J) J } p =-Re{ adj(J) delta d} 
!
!     INPUT      MEANING
!     -----      ------- 
!     BDHESS     block diagonal hessian
!     NA35       number of variables in inversion = NNPINV*NVINV
!     NHSIZE     number of non-zeros in hessian preconditioner 
!     NNPINV     number of nodal points in inversion
!     NVINV      number of variables to invert for each nodal point
!     PTHRESH    (0,100) if a singular value falls below the thresh 
!                percentage it will be inflated
!
!     OUTPUT     MEANING
!     ------     ------- 
!     GRADPC     gradient preconditioner
!
!.... variable declarations
      IMPLICIT NONE
      REAL*4, INTENT(IN) :: BDHESS(NHSIZE), PTHRESH 
      INTEGER*4, INTENT(IN) :: NHSIZE, NA35, NNPINV, NVINV
      REAL*4, INTENT(OUT) :: GRADPC(NHSIZE)
!.... local variables
      REAL*4, ALLOCATABLE :: BDIAG(:,:), EYE(:,:), WMAT(:,:), 
     ;                       U(:,:), VT(:,:), WORK(:), SORT(:)  
      INTEGER*4, ALLOCATABLE :: IWORK(:) 
      REAL*4 THRESH, CUT
      INTEGER*4 LWORK, ICUT, INPINV, IVINV, JVINV, INDX, JNDX, 
     ;          IA35, I1, I2, INFO  
!  
!----------------------------------------------------------------------!
! 
!.... quick check on thresholding
      THRESH = PTHRESH/100.0
      IF (THRESH.LE.0.0) THEN
         WRITE(*,*) 'regprecon: Threshholding <= 0 percent!  Skipping!'
         RETURN
      ENDIF
      IF (THRESH.GE.1.0) THEN 
         WRITE(*,*) 'regprecon: Thresholding >= 100 percent! Skipping!'
         RETURN
      ENDIF
!
!.... determine if problem is diagonal or block diagonal
      IF (NVINV.EQ.1) THEN
         ALLOCATE(SORT(NA35))
         CALL SCOPY(NA35,BDHESS,1,SORT,1)
         CALL RSHELL1(NA35,SORT)
         ICUT = IFIX(THRESH*FLOAT(NA35) + 0.5) 
         CUT = SORT(ICUT) 
         WRITE(*,*) 'regprecon: Inflating small singular values by:',CUT
         DO 1 IA35=1,NA35
            GRADPC(IA35) = BDHESS(IA35) + CUT 
   1     CONTINUE 
         DEALLOCATE(SORT) 
      ELSE
         ALLOCATE(BDIAG(NVINV,NVINV)) 
         ALLOCATE(SORT(NA35))
         LWORK = 3*NVINV + 4*NVINV*NVINV + 4*NVINV
         ALLOCATE(WORK(LWORK)) 
         ALLOCATE(IWORK(8*NVINV)) !8*min(M,N)
         ALLOCATE(U(NVINV,NVINV))
         ALLOCATE(VT(NVINV,NVINV))
         INDX = 0 
         DO 2 INPINV=1,NNPINV !loop on inversion point nodes 
            DO 3 IVINV=1,NVINV !extract a block 
               DO 4 JVINV=1,NVINV 
                  INDX = INDX + 1 
                  BDIAG(IVINV,JVINV) = BDHESS(INDX)
    4          CONTINUE !loop on columns
    3       CONTINUE !loop on rows of submatrix 
            I1 = (INPINV - 1)*NVINV + 1
            I2 = INPINV*NVINV
            CALL SGESDD('N',NVINV,NVINV,BDIAG,NVINV, SORT(I1:I2), 
     ;                  U,NVINV,VT,NVINV, WORK,LWORK,IWORK,INFO)
            IF (INFO.LT.0) THEN
               WRITE(*,*) 'regprecon: i-th argument illegal1:',INFO
            ELSEIF (INFO.GT.0) THEN
               WRITE(*,*) 'regprecon: SBDSDC didnt converge'
            ENDIF 
    2    CONTINUE
         CALL RSHELL1(NA35,SORT)
         ICUT = IFIX(THRESH*FLOAT(NA35) + 0.5) 
         CUT = SORT(ICUT)
         WRITE(*,*) 'regprecon: Inflating small singular values by:',CUT
         INDX = 0
         JNDX = 0
         ALLOCATE(EYE(NVINV,NVINV))
         ALLOCATE(WMAT(NVINV,NVINV))
         EYE(1:NVINV,1:NVINV) = 0.0 !CMPLX(0.0,0.0)
         DO 5 INPINV=1,NNPINV !loop on inversion point nodes 
            DO 6 IVINV=1,NVINV !extract a block 
               DO 7 JVINV=1,NVINV
                  INDX = INDX + 1
                  BDIAG(IVINV,JVINV) = BDHESS(INDX)
    7          CONTINUE !loop on columns
    6       CONTINUE !loop on rows of submatrix 
            I1 = (INPINV - 1)*NVINV + 1
            I2 = INPINV*NVINV
            CALL SGESDD('A',NVINV,NVINV,BDIAG,NVINV, SORT(1:NVINV), 
     ;                  U,NVINV,VT,NVINV, WORK,LWORK,IWORK,INFO)
            IF (INFO.LT.0) THEN
               WRITE(*,*) 'regprecon: i-th argument illegal2:',INFO
            ELSEIF (INFO.GT.0) THEN
               WRITE(*,*) 'regprecon: CBDSDC didnt converge'
            ENDIF
            IF (INFO.NE.0) THEN
               WRITE(*,*) 'regprecon: Skipping...'
               GOTO 50
            ENDIF
            DO 8 IVINV=1,NVINV
               EYE(IVINV,IVINV) = SORT(IVINV) + CUT 
    8       CONTINUE  
            WMAT = MATMUL(EYE,VT)
            BDIAG = MATMUL(U,WMAT)
            DO 11 IVINV=1,NVINV
               DO 12 JVINV=1,NVINV
                  JNDX = JNDX + 1
                  GRADPC(JNDX) = BDIAG(IVINV,JVINV) 
   12          CONTINUE
   11       CONTINUE 
   50       CONTINUE
    5    CONTINUE
         DEALLOCATE(BDIAG)
         DEALLOCATE(SORT)
         DEALLOCATE(WORK)
         DEALLOCATE(IWORK)
         DEALLOCATE(EYE)
         DEALLOCATE(WMAT)
         DEALLOCATE(U)
         DEALLOCATE(VT)
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RSCLVEC(N, C,D, XMIN,XMAX, X,IERR)
!
!     Rescales a vector to into a specified range [c,d] while 
!     preserving relative distance from end points
!
!     INPUT      MEANING
!     -----      ------- 
!     C          lower bound to rescale vector 
!     D          upper bound to rescale vector
!     N          number of points in gradient
!     X          vector to rescale
!
!     OUTPUT     MEANING
!     ------     ------- 
!     X          vector rescaled to [xmin,xmax]
!     XMAX       upper value to rescale to  max(abs(x))
!     XMIN       lower value to rescale to -max(abs(x)) 
! 
!.... variable declarations
      IMPLICIT NONE
      REAL*4, INTENT(INOUT) :: X(N), XMIN, XMAX
      REAL*4, INTENT(IN) :: C, D
      INTEGER*4, INTENT(IN) :: N 
      INTEGER*4, INTENT(OUT) :: IERR
      REAL*4 STFORM
      INTEGER*4 I 
!
!----------------------------------------------------------------------!
!
!.... want to preserve relative distance from endpoints 
      IF (XMIN.EQ.0.0 .AND. XMAX.EQ.0.0) THEN
         XMAX = 0.0 
         DO 1 I=1,N
            XMAX = AMAX1(ABS(X(I)),XMAX)  
    1    CONTINUE 
!
!....... this would be bad, determinant undefined
         IERR = 0
         IF (XMAX.EQ.0.0) THEN
            WRITE(*,*) 'rcsclvec: Error x is all zeros!'
            IERR = 1
            RETURN 
         ENDIF 
         XMIN =-XMAX 
      ENDIF
!
!.... now rescale
      DO 2 I=1,N
         X(I) = STFORM(XMIN,XMAX, C,D,X(I))   
    2 CONTINUE 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*4 FUNCTION STFORM(A,B, C,D,X)
!     
!     transforms x in x in [a,b] to interval [c,d] 
!
!     From Hughes we have [A,B] -> [-1,1] is given by solving 
!     c_1, c_2 in 
!       -1 = c_1 + A*c_2
!        1 = c_1 + B*c_2 
!     so more generally we can do 
!        c = c_1 + A*c_2 
!        b = c_2 + B*c_2  
! 
      REAL*4, INTENT(IN) :: A,B, C,D,X
      IF (B == A) THEN 
         WRITE(*,*) 'stform; Error A = B, determinant undefined'
         STFORM = 1.0  
         RETURN
      ENDIF
      DET = 1.0/(B - A) 
      C1 = DET*(B*C - A*D) 
      C2 = DET*(-C + D) 
      XI = C1 + X*C2 
      STFORM = XI 
      RETURN
      END  
