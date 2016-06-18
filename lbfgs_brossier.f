      SUBROUTINE LBFGS( ) 

      DO 1 I=1,N
         WORK(I) = DIAG(I)*G(I)
    1 CONTINUE 
!
!.... shift memory
      IF (ITER > M) THEN
         DO 2 MEM=M,2,-1
            I1 = (MEM - 1)*N + 1 
            I2 = I1 + N - 1
            J1 = (MEM - 2)*N + 1
            J2 = J1 + N - 1
            CALL SCOPY(N,Y(I1:I2),1,Y(J1:J2),1)
            CALL SCOPY(N,S(I1:I2),1,S(J1:J2),1)
    2    CONTINUE
      ENDIF 
      I1 = MIN((ITER - 1)*N + 1,(M - 1)*N)  
      I2 = I1 + N - 1
      IF (ITER.GT.1) THEN 
         S(I1:I2) = X(1:I2) - XOLD(I1:I2)  
         Y(I1:I2) = D(1:I2) - DOLD(I1:I2)   
      ENDIF
!
!....
      SUBROUTINE TWO_LOOP( ) 
      ALLOCATE(WORK(N)) 
      CALL SCOPY(N,D,1,WORK,1) !q <- d
      JBEG = ITER
      JEND = MIN(1,ITER-M) 
      I = MIN(M,ITER) + 1  
      DO 1 J=JBEG,JEND-M,-1
         I = I - 1
         I1 = (I - 1)*N + 1
         I2 = I1 + N - 1
         RHO(I) = SDOT(N,Y(I1:I2),1,S(I1:I2),1) 
         ALPHA(I) = RHO*SDOT(N,S(I1:I2),1,WORK,1)
         CALL SAXPY(N,-ALPHA,Y(I1:I2),1,WORK,1) !q = q - alpha y 
    1 CONTINUE 
      IF (ITER.GT.1) THEN
         I1 = (M - 1)*N + 1
         I2 = I1 + N - 1
         GAMMA1 = SDOT(N,S(I1:I2),1,Y(I1:I2),1)
     ;           /SDOT(N,Y(I1:I2),1,Y(I1;i2),1)
      ENDIF
      I = 0
      CALL SCAL(N,GAMMA1,WORK,1) !r = gamma q
      DO 2 J=JBEG,JEND,+1
         I = I + 1
         I1 = (I - 1)*N + 1
         I2 = I1 + N - 1
         BETA = RHO(I)*SDOT(N,Y(I1:I2),1,WORK,1)
         CALL SAXPY(N,(ALPHA(I) - BETA),S(I1:I2),1,WORK,1) !r=r+(a-b)*s
    2 CONTINUE 
      CALL SSCAL(N,-1.0,WORK,1) !p =-r
      CALL SCOPY(N,WORK,1,P,1)  !finish with search direction
      DEALLOCATE(WORK) 
      RETURN
      END
