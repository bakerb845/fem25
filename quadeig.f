      SUBROUTINE QPEEIG2(LDA,LDX,N, A0,A1,A2, EIGS,X, IERR ) 
!
!     Fortran implementation of polynomial eigenvalue solver:
!     http://mind.cog.jhu.edu/courses/680/octave/Installers/
!     Octave/Octave.OSX10.6/Applications/MATLAB_R2009b.app/
!     toolbox/matlab/matfun/polyeig.m
!  
!     For (A_0 + lambda A_1 + lambda^2 A_2)x = 0
!
!     Step 1: Build the [n*p x n*p] matrices
!
!     A = [A0  0]   B = [-A1  -A2]
!         [ 0  I]     = [  I    0]
!
      COMPLEX*16, INTENT(IN) :: A0(LDA,*), A1(LDA,*), A2(LDA,*)
      COMPLEX*16, INTENT(OUT) :: X(LDX,*), EIGS(*)
      INTEGER*4, INTENT(IN) :: LDA, LDX, N 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*16, ALLOCATABLE :: A(:,:), B(:,:), VL(:,:), VR(:,:), 
     ;            R(:,:), EYES(:,:), ALPHA(:), BETA(:), 
     ;            WORK(:), R1(:) 
      REAL*8, ALLOCATABLE :: RWORK(:)  
      LOGICAL*4, ALLOCATABLE :: LKEEP(:)
      CHARACTER(1) JOBVL, JOBVR 
      COMPLEX*16 ZR2, CZERO, CONE, ZDOTC 
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0))
      PARAMETER(CONE  = DCMPLX(1.D0,0.D0)) 
!
!.... build the two [n*p x n*p] matrices
      IERR = 0
      NN = (IP + 1)*N 
      ALLOCATE(A(NN,NN))
      ALLOCATE(B(NN,NN)) 
     
      I1 = 1  
      DO 1 I=1,NN
         I1 = I1 + 1
         IF (I1.GT.N) I1 = 1 
         JTRIP = 1 
         IP = 0 
         DO 2 J=1,NN
            J1 = J1 + 1
            IF (J1.GT.N) THEN
               IP = IP + 1 
               J1 = 1 
            ENDIF
            IF (I.LE.N .AND. J.LE.N) THEN
               A(I,J) = A0(I1,J1) 
            ELSE
               A(I,J) = CZERO 
            ENDIF 
            B(I,J) = CZERO
            IF (I.LE.N) THEN 
               IF (J.LE.N) THEN !fill matrices across top
                  IF (IP.EQ.0) THEN
                     B(I,J) =-A1(I1,J1) 
                  ELSEIF(IP.EQ.1) THEN
                     B(I,J) =-A2(I1,J1)
                  ELSE
                     WRITE(*,*) 'qpeeig2: Error ip > 2!' 
                     IERR = 1
                     RETURN
                  ENDIF 
               ENDIF 
            ELSE
               IF (I.GE.(IP+1)*N .AND. I.LE.(IP+2)*N) THEN
                  IF (I1.EQ.J1) B(I,J) = CONE
               ENDIF
            ENDIF
    2    CONTINUE 
    1 CONTINUE  
!
!.... solve the right generalized eigenvalue problem
      JOBVL = 'N'
      JOBVR = 'V' 
      LDVL = 1
      LDVR = NN
      ALLOCATE(VL(1,NN))
      ALLOCATE(VR(LDVR,N))
      ALLOCATE(ALPHA(NN)) 
      ALLOCATE(BETA(NN)) 
      LWORK = MAX(1,2*NN) + 10 !add some extra
      ALLOCATE(WORK(LWORK)) 
      ALLOCATE(RWORK(8*NN))  
      CALL ZGGEV(JOBVL,JOBVF,NN,A,NN,B,NN, ALPHA,BETA, 
     ;           VL, LDVL, VR, LDVR, WORK,LWORK, RWORK, INFO)  
      DEALLOCATE(WORK)
      DEALLOCATE(RWORK)
      DEALLOCATE(VL)
      IF (INFO.LT.0) THEN
         WRITE(*,*) 'qpeeig2: Error illegal argument:',INFO
         IERR = 1
         RETURN
      ELSE
         IF (INFO.GE.1 .AND.INFO.LE.NN) THEN
            WRITE(*,*) 'qpeeig2: No eigenvectors, eigenvalues up to',
     ;      INFO
            IERR = 1
            RETURN
         ELSE 
            IF (INFO.EQ.NN+1) THEN
               WRITE(*,*) 'qpeeig2: zggev failed in dhgeqz'
            ELSE
               WRITE(*,*) 'qpeeig2: zggev failed in dtgevc'
            ENDIF
            IERR = 1
            RETURN
         ENDIF
      ENDIF
      ALLOCATE(LKEEP(NN)) 
      DO 10 I=1,NN
         IF (CDABS(BETA(I)).GT.1.11D-16) THEN 
            EIGS(I) = ALPHA(I)/BETA(I)
            LKEEP(I) = .TRUE.
         ELSE
            EIGS(I) = CZERO
            LKEEP(I) = .FALSE.
         ENDIF
   10 CONTINUE  
      DEALLOCATE(ALPHA)
      DEALLOCATE(BETA) 
!
!.... initialize space and identity matrix for next step 
      ALLOCATE(R(N,N)) 
      ALLOCATE(R1(N)) 
      ALLOCATE(EYES(N,N)) 
      DO 15 I=1,N
         DO 16 J=1,N
            IF (I.EQ.J) THEN
               EYES(I,J) = CONE
            ELSE
               EYES(I,J) = CZERO
            ENDIF
   16   CONTINUE
   15 CONTINUE
!
!.... extract eigvec/eigval pair of big eigenvect X w/ smallest resid 
      DO 20 J=1,NN
         DO 21 I=1,N
            CALL ZCOPY(N,A2(I,1:N),1,R(I,1:N),1)
   21    CONTINUE !R = varagin(p+1)
         IF (LKEEP(J)) THEN
            DO 22 K=2,1,-1
               IF (K.EQ.2) THEN !
                  CALL ZGEMM('N','N',N,N,N,CONE,A1,LDA,EYES,N, CONE,R,N)
               ELSE
                  CALL ZGEMM('N','N',N,N,N,CONE,A0,LDA,EYES,N, CONE,R,N)
               ENDIF 
   22       CONTINUE
         ENDIF
         CALL ZGEMV('N',N,N,CONE,R,N,VR(1:NN,J),1, CZERO,R1,1) 
         RESN = 0.D0
         RESD = 0.D0
         DO 23 I=1,N 
            RESN = RESN + CDABS(R1(I)) 
            RESD = RESD + CDABS(X(I,J))
   23    CONTINUE
         IND = 1
         RES = RESN/RESD
         R2 = DREAL(ZDOTC(N,V,1,V,1))  
         ZR2 = DCMPLX(R2,0.D0)
         DO 24 I=1,N
            X(I,J) = V(I)/ZR2
   24    CONTINUE 
   20 CONTINUE  
      DEALLOCATE(VL)  
      RETURN
      END
