      SUBROUTINE FOBJWRK(N,LFUNC,LGRAD,LHESS,F,X,G,H)

      REAL*4 X(N), G(N), H(N), F
      LOGICAL*4, INTENT(IN) :: LFUNC, LGRAD, LHESS
      IF (LFUNC) THEN
         F = 0.0
         DO 1 J=1,N,2
            T1 = 1.0 - X(J)
            T2 = 1.E1*(X(J+1) - X(J)**2)
            F = F + T1**2 + T2**2
    1    CONTINUE
      ENDIF
      IF (LGRAD) THEN
         DO 2 J=1,N,2
            T1 = 1.0 - X(J)
            T2 = 1.E1*(X(J+1) - X(J)**2)
            G(J+1) = 2.E1*T2
            G(J) =-2.0*(X(J)*G(J+1) + T1) 
    2    CONTINUE
      ENDIF
      IF (LHESS) THEN
         DO 3 J=1,N,2
            H(J+1) = 200.0 
            H(J) =-400.*(X(J+1) - 3.*X(J)**2) + 2.*X(J)
    3    CONTINUE
         h(1:n) = 1.0/h(1:n)
      ENDIF
      RETURN
      END
