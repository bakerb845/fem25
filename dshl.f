      SUBROUTINE DSHL(MEN,MINTX,MINTZ, IITYPE,NINTX,NINTZ, NLXI,NLETA,
     ;                XIIN,ETAIN, XIGLL,ETAGLL,SHL)
! 
!     Calculate N_a(xi,eta), N_a,xi(xi,eta), N_a,eta,(xi,eta) 
! 
!     Based on page 149 of Hughes, so we'll use the 
! 
!     N_a    (xi,eta)     SHL(3,I,IX,IZ) 
!     N_a,xi (xi,eta)     SHL(1,I,IX,IZ)     NROWSH X NEN X NLX X NLZ 
!     N_a,eta(xi,eta)     SHL(2,I,IX,IZ)        
! 
!     Given the standard element: 
!
!     -   -   -   -   - 
!     4   12  11  10  3
! 
!     -   -   -   -   -
!     13  18  19  20  9
!
!     -   -   -   -   -
!     14  15  16  17  8
!
!     -   -   -   -   - 
!     1   5   6   7   2
! 
!     Then we can define 
! 
!     N_1  = L_1(xi) L_1(eta)  
!     N_2  = L_5(xi) L_1(eta) 
!     N_3  = L_5(xi) L_4(eta) 
!     N_4  = L_1(xi) L_4(eta) 
!     N_5  = L_2(xi) L_1(eta) 
!      . 
!      . 
!      . 
!     N_20 = L_4(xi) L_3(eta) 
! 
!     However, to reduce the amount of computations we will be 
!     methodical and order the element
!
!     -   -   -   -   - 
!     16  17  18  19  20
! 
!     -   -   -   -   -
!     11  12  13  14  15  
!
!     -   -   -   -   -
!     6   7   8   9   10 
!
!     -   -   -   -   - 
!     1   2   3   4   5
! 
!     N_1  = L_1(xi) L_1(eta)  
!     N_2  = L_2(xi) L_1(eta) 
!     N_3  = L_3(xi) L_1(eta) 
!     N_4  = L_4(xi) L_1(eta) 
!     N_5  = L_5(xi) L_1(eta) 
!     N_6  = L_1(xi) L_2(eta)
!      . 
!      . 
!      . 
!     N_20 = L_5(xi) L_4(eta) 
!
!     INPUT      MEANING 
!     -----      ------- 
!     IITYPE     = 1 gauss lobatto legendre interpolation points
!                = 2 chebyshev gauss lobatto interpolation points
!                = 3 custom points 
!     NLETA      number of GLL points in eta
!     NLXI       number of GLL points in xi 
!     MEN        max number of element nodes; NLX*NLZ
!     MINTX      max number of integration points in x 
!     MINTZ      max number of integration points in z
!     NEN        number of element nodes
!     NINTX      number of integration points in x 
!     NINTZ      number of integration points in z  
!
!     OUTPUT     MEANING 
!     ------     ------- 
!     ETAGLL     integration abscissas/weights in eta
!     NEN        number of element nodes = NLX*NLZ 
!     NINT1      number of integration points = NLX*NLZ 
!     SHL        lagrange shape functions, see above  
!     XIGLL      integration abscissas/weights in xi
! 
!     VARIABLE   MEANING 
!     --------   ------- 
!     COFXI      coefficients of lag. polys in xi 
!     COFETA     coefficients of lag. polys in eta 
! 
!.... variable declarations 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XIIN(NLXI), ETAIN(NLETA) 
      REAL*8 XIGLL(MINTX,2), ETAGLL(MINTZ,2), SHL(3,MEN,MINTX,*)
!.... local variables 
      ALLOCATABLE XIL(:), ETAL(:), COFXI(:), COFETA(:), 
     ;            S(:), WTWORK(:), LAGR(:), XIPTS(:),
     ;            ETAPTS(:) 
      REAL*8 ENDPTS(2), ALPHA, BETA, LAGR, XIPTS, ETAPTS, 
     ;       XIL, ETAL, COFXI, COFETA, S, SUMD, EPS
      REAL*8 SHAPEF
      PARAMETER(ALPHA = 0.D0, BETA = 0.D0) 
      PARAMETER(EPS = 2.11D-16)
      DATA ENDPTS /-1.D0, 1.D0/ 
      PARAMETER(IKIND = 5, KPTS = 2) 
! 
!----------------------------------------------------------------------!
! 
!.... generate the interpolation (GLL or CGL) points for the element
      ALLOCATE(XIPTS (NLXI))
      ALLOCATE(ETAPTS(NLETA))
      ALLOCATE(COFXI (NLXI))
      ALLOCATE(COFETA(NLETA))
      IF (NLXI.GT.NLETA) THEN
         ALLOCATE(WTWORK(NLXI))
         ALLOCATE(S(NLXI))
         ALLOCATE(LAGR(NLXI))
      ELSE
         ALLOCATE(WTWORK(NLETA)) 
         ALLOCATE(S(NLETA))
         ALLOCATE(LAGR(NLETA)) 
      ENDIF 
      IF (IITYPE.EQ.1) THEN !Gauss Legendre Lobatto points
         CALL GAUSSQ(IKIND,NLXI, ALPHA,BETA,KPTS,ENDPTS, 
     ;               XIPTS, WTWORK)
         CALL GAUSSQ(IKIND,NLETA,ALPHA,BETA,KPTS,ENDPTS,
     ;               ETAPTS,WTWORK)
      ELSEIF(IITYPE.EQ.2) THEN !Chebyshev Gauss Lobatto points 
         CALL CBLPTS(NLXI , XIPTS)
         CALL CBLPTS(NLETA, ETAPTS)  
      ELSE
         XIPTS(1:NLXI) = XIIN(1:NLXI)
         ETAPTS(1:NLETA) = ETAIN(1:NLETA)  
      ENDIF 
      IF (MOD(NLXI, 2).NE.0)  XIPTS( (NLXI -1)/2 + 1) = 0.D0 
      IF (MOD(NLETA,2).NE.0) ETAPTS( (NLETA-1)/2 + 1) = 0.D0 
! 
!.... set the integration points and weights on GLL points 
      ALLOCATE(XIL (NINTX)) !derivative of xi shape fns at int pts
      ALLOCATE(ETAL(NINTZ)) !derivative of eta shape fns at int pts 
      CALL GAUSSQ(IKIND,NINTX,ALPHA,BETA,KPTS,ENDPTS, 
     ;            XIGLL(:,1), XIGLL(:,2))
      CALL GAUSSQ(IKIND,NINTZ,ALPHA,BETA,KPTS,ENDPTS, 
     ;            ETAGLL(:,1),ETAGLL(:,2))
      IF (MOD(NINTX,2).NE.0) XIGLL( (NINTX - 1)/2 + 1,1) = 0.D0
      IF (MOD(NINTZ,2).NE.0)ETAGLL( (NINTZ - 1)/2 + 1,1) = 0.D0
! 
!.... loop over all the elements nodes in eta 
      IA = 0 !element node counter
      DO 1 IETA=1,NLETA !loop in z direction 
! 
!....... set the interpolation points for eta lagrange polynomial 
         DO 3 K=1,NLETA
            IF (K.EQ.IETA) THEN
               LAGR(K) = 1.D0
            ELSE
               LAGR(K) = 0.D0
            ENDIF
    3    CONTINUE
         CALL POLCOE(ETAPTS,LAGR,NLETA, COFETA, S)
! 
!....... differentiate w.r.t. eta and evaluate at integration points
         DO 4 INTZ=1,NINTZ
            SUMD = 0.D0
            DO 5 K=NLETA,2,-1
               DK = DFLOAT(K - 1)
               SUMD = SUMD + DK*COFETA(K)*ETAGLL(INTZ,1)**(K - 2)
    5       CONTINUE
            IF (DABS(SUMD).LT.EPS) SUMD = 0.D0
            ETAL(INTZ) = SUMD
    4    CONTINUE 
! 
!....... loop over xi points and repeat
         DO 2 IXI=1,NLXI
            IA = IA + 1
            DO 6 K=1,NLXI
               IF (K.EQ.IXI) THEN 
                  LAGR(K) = 1.D0
               ELSE
                  LAGR(K) = 0.D0
               ENDIF
    6       CONTINUE
            CALL POLCOE(XIPTS,LAGR,NLXI, COFXI, S)
            DO 7 INTX=1,NINTX
               SUMD = 0.D0
               DO 8 K=NLXI,2,-1
                  DK = DFLOAT(K - 1)
                  SUMD = SUMD + DK*COFXI(K)*XIGLL(INTX,1)**(K - 2) 
    8          CONTINUE
               IF (DABS(SUMD).LT.EPS) SUMD = 0.D0
               XIL(INTX) = SUMD
    7       CONTINUE
! 
!.......... now evaluate shape functions at integration points
            L = 0 
            DO 10 INTZ=1,NINTZ
               DO 11 INTX=1,NINTX
! 
!................ my mass terms will have kronecker delta
                  L = L + 1
                  IF (NINTX.EQ.NLXI .AND. NINTZ.EQ.NLETA) THEN
                     IF (L.EQ.IA) THEN
                        SHL(3,IA,INTX,INTZ) = 1.D0
                     ELSE
                        SHL(3,IA,INTX,INTZ) = 0.D0
                     ENDIF 
                  ELSE 
                     SNA = SHAPEF(NLETA,ETAGLL(INTZ,1),COFETA)
                     SNB = SHAPEF(NLXI , XIGLL(INTX,1),COFXI )
                     SHL(3,IA,INTX,INTZ) = SNA*SNB !N_a*N_b 
                     IF (DABS(SHL(3,IA,INTX,INTZ)).LT.EPS) THEN
                        SHL(3,IA,INTX,INTZ) = 0.D0 
                     ENDIF 
                  ENDIF
                  SNA = SHAPEF(NLETA,ETAGLL(INTZ,1),COFETA)
                  SHL(1,IA,INTX,INTZ) = SNA*XIL (INTX) !N_a,xi *N_a(eta)
                  SNA = SHAPEF(NLXI ,XIGLL (INTX,1), COFXI)
                  SHL(2,IA,INTX,INTZ) = SNA*ETAL(INTZ) !N_a,eta*N_a(xi)
                  IF (DABS(SHL(1,IA,INTX,INTZ)).LT.EPS) 
     ;               SHL(1,IA,INTX,INTZ) = 0.D0
                  IF (DABS(SHL(2,IA,INTX,INTZ)).LT.EPS) 
     ;               SHL(2,IA,INTX,INTZ) = 0.D0
   11          CONTINUE
   10       CONTINUE
    2    CONTINUE !end loop on eta
    1 CONTINUE !end loop on xi
! 
!.... clean up
c     IF (IITYPE.NE.3) THEN
c        XIIN(1:NLXI) = XIPTS(1:NLXI)
c        ETAIN(1:NLETA) = ETAPTS(1:NLETA) 
c     ENDIF   
      DEALLOCATE(XIPTS)
      DEALLOCATE(ETAPTS)
      DEALLOCATE(COFXI )
      DEALLOCATE(COFETA)
      DEALLOCATE(WTWORK)
      DEALLOCATE(S)
      DEALLOCATE(XIL)
      DEALLOCATE(ETAL) 
      DEALLOCATE(LAGR)
      RETURN 
      END 
!                                                                      !
!======================================================================!
!                                                                      ! 
      FUNCTION SHAPEF(N,XI,COFS) 
! 
!     evaluates the shape function 
!     cof(n)*xi**n + cof(n-1)*xi**(n-1) + ... + cof(1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 COFS(N)
      PARAMETER(EPS = 1.11D-15)
      SHAPEF = 0.D0
      DO 1 I=N,2,-1
         SHAPEF = SHAPEF + COFS(I)*XI**(I - 1)
    1 CONTINUE
      SHAPEF = SHAPEF + COFS(1)
      IF (DABS(SHAPEF).LT.EPS) SHAPEF = 0.D0
      RETURN  
      END
!                                                                      ! 
!======================================================================!
!                                                                      !

      SUBROUTINE POLCOE(X, Y, N, COF, S) 
! 
!     this method is of order n^2 however it is less stable than 
!     polcof which is of order n^3 because of the vandermonde matrix
! 
!     given the arrays x(1:n) and y(1:n) containing a tabulated fn. 
!     y_i=f(x_i), routine returns the coefficents cof(1:n) s.t. 
!     y_i = sum_j cof_j x_i^{j-1}.  algorithm from numerical recipes 
!     pg 115  
! 
!     modification: we don't need to tabulate y here b/c its either 0 
!     or 1 
! 
!     varible declarations 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 COF(N), X(N), Y(N), S(N)  
      INTEGER*4 I,J,K 
      REAL*8 B, FF, PHI
! 
!----------------------------------------------------------------------!
! 
!.... zero initialization 
      DO 1 I=1,N 
         S(I)=0.D0  
         COF(I)=0.D0 
    1 CONTINUE 
      S(N)=-X(1) 
!
!.... coeffs s_i of the master polynomial p(x) are found by recurrence 
      DO 2 I=2,N 
         DO 3 J=N+1-I,N-1 
            S(J)=S(J)-X(I)*S(J+1) 
    3    CONTINUE 
         S(N)=S(N)-X(I) 
    2 CONTINUE 
      DO 4 J=1,N 
         PHI=DFLOAT(N) 
         DO 5 K=N-1,1,-1 
            PHI=DFLOAT(K)*S(K+1)+X(J)*PHI 
    5    CONTINUE 
         FF=Y(J)/PHI 
         B=1.D0 
         DO 6 K=N,1,-1 
            COF(K)=COF(K)+B*FF 
            B=S(K)+X(J)*B 
    6    CONTINUE 
    4 CONTINUE 
      RETURN 
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CBLPTS(NPTS, CPTS) 
!     Generates the Chebyshev-Gauss-Lobatto grid points 
!     http://en.wikipedia.org/wiki/Chebyshev_polynomials 
      INTEGER*4, INTENT(IN) :: NPTS
      REAL*8, INTENT(OUT) :: CPTS(NPTS) 
      REAL*8 PI, EPS, ARG
      PARAMETER(PI = 3.1415926535897931D0) 
      PARAMETER(EPS = 1.11D-16)
      DO 1 I=0,NPTS-1
         K = I + 1 
         ARG = DFLOAT(I)/DFLOAT(NPTS - 1)  
         CPTS(K) =-DCOS(ARG*PI) 
         IF (DABS(CPTS(K)).LT.EPS) CPTS(K) = 0.D0
    1 CONTINUE
      RETURN 
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE EQLPTS(NPTS, EPTS)
!     generates the npt equally spaced points [-1,1]
      INTEGER*4, INTENT(IN) :: NPTS
      REAL*8, INTENT(OUT) :: EPTS(NPTS)
      REAL*8 EPS
      PARAMETER(EPS = 1.11D-16) 
      DO 1 I=1,NPTS
         EPTS(I) = -1.D0 + DFLOAT(I - 1)*2.D0/DFLOAT(NPTS -1)
         IF (DABS(EPTS(I)).LT.EPS) EPTS(I) = 0.D0
    1 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GENPTS(IITYPE, NLXI,NLETA, XIPTS,ETAPTS) 
!     
!     INPUT      MEANING
!     -----      ------- 
!     IITYPE     interpolation type (1) GLL default (2) Cheby (3) equal
!     NLETA      number of lagrange interpolant points in eta
!     NLXI       number of lagrange interpolant points in xi 
!  
!     OUTPUT     MEANING
!     ------     ------- 
!     ETAPTS     lagrange interpolant points in eta
!     XIPTS      lagrange interpolant points in xi 
!  
!.... variable declarations
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8, INTENT(OUT) :: XIPTS(NLXI),ETAPTS(NLETA)
      INTEGER*4, INTENT(IN) :: IITYPE, NLXI,NLETA
!.... local variables
      ALLOCATABLE WTWORK(:)
      REAL*8 ENDPTS(2),WTWORK,ALPHA,BETA,EPS
      PARAMETER(ALPHA = 0.D0, BETA = 0.D0) 
      PARAMETER(EPS = 1.11D-16)
      DATA ENDPTS /-1.D0, 1.D0/ 
      PARAMETER(IKIND = 5, KPTS = 2)  
! 
!----------------------------------------------------------------------!
! 
!.... generate the interpolation (GLL or CGL) points for the element
      IF (IITYPE.EQ.3) THEN
         CALL EQLPTS(NLXI ,XIPTS )
         CALL EQLPTS(NLETA,ETAPTS)
      ELSEIF(IITYPE.EQ.2) THEN !Chebyshev Gauss Lobatto points 
         CALL CBLPTS(NLXI , XIPTS)
         CALL CBLPTS(NLETA, ETAPTS)  
      ELSE !default GLL nodes
         IF (IITYPE.NE.1) WRITE(*,*) 'genpts: Defaulting to GLL nodes'
         ALLOCATE(WTWORK(MAX0(NLXI,NLETA))) 
         CALL GAUSSQ(IKIND,NLXI, ALPHA,BETA,KPTS,ENDPTS, 
     ;               XIPTS, WTWORK)
         CALL GAUSSQ(IKIND,NLETA,ALPHA,BETA,KPTS,ENDPTS, 
     ;               ETAPTS,WTWORK)
         DEALLOCATE(WTWORK) 
      ENDIF 
      IF (MOD(NLXI, 2).NE.0)  XIPTS( (NLXI -1)/2 + 1) = 0.D0 
      IF (MOD(NLETA,2).NE.0) ETAPTS( (NLETA-1)/2 + 1) = 0.D0 
      RETURN
      END
