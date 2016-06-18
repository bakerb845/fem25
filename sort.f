      SUBROUTINE ISHELL1(N, A)
!
!     sorts integer vector A.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
      INTEGER*4, INTENT(INOUT) :: A(N) 
      INTEGER*4, INTENT(IN) :: N
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             ITEMP  = A(I) 
             J = I  
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.ITEMP) THEN 
                A(J) = A(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = ITEMP
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END  
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE ISHELL2(N, A,B) 
!
!     sorts integer vector A.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
      INTEGER*4 A(N),B(N) 
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             ITEMP  = A(I)
             ITEMP2 = B(I)
             J = I
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.ITEMP) THEN
                A(J) = A(J-INC)
                B(J) = B(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = ITEMP
             B(J) = ITEMP2
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SHELLID4(N, A,B,C,D,E)
!
!     sorts integer vectors.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
!     INPUT        MEANING
!     -----        -------
!     A            unsorted list to sort on 
!     N            length of vectors
!
!.... formal variables  
      REAL*8 A(N), B(N) 
      INTEGER*4 C(N), D(N), E(N)
      REAL*8 TEMP, TEMP2 
      INTEGER*4 TEMP3, TEMP4, TEMP5
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             TEMP  = A(I)
             TEMP2 = B(I)
             TEMP3 = C(I)
             TEMP4 = D(I)
             TEMP5 = E(I)
             J = I
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.TEMP) THEN
                A(J) = A(J-INC)
                B(J) = B(J-INC)
                C(J) = C(J-INC)
                D(J) = D(J-INC)
                E(J) = E(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = TEMP
             B(J) = TEMP2
             C(J) = TEMP3
             D(J) = TEMP4
             E(J) = TEMP5
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RSHELL1(N, A)
!
!     sorts integer vector A.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
      INTEGER*4, INTENT(IN) :: N
      REAL*4, INTENT(INOUT) :: A(N) 
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             TEMP  = A(I) 
             J = I  
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.TEMP) THEN 
                A(J) = A(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = TEMP
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END  
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CSHELL1(N, A)
!
!     sorts integer vector A.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
      INTEGER*4, INTENT(IN) :: N 
      COMPLEX*8, INTENT(INOUT) :: A(N) 
      COMPLEX*8 TEMP
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             TEMP  = A(I) 
             J = I  
!
!........... inner loop of straight insertion
    3        IF (CABS(A(J-INC)).GT.CABS(TEMP)) THEN 
                A(J) = A(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = TEMP 
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END  
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SHELL4(N, A,B,C,D)
!
!     sorts integer vectors.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
!     INPUT        MEANING
!     -----        -------
!     A            unsorted list to sort on 
!     N            length of vectors
!
!     OUTPUT       MEANING 
!     ------       ------- 
!     A            sorted list 
!     B            sorted with A 
!     C            sorted with A
!     D            sorted with A
!
!.... formal variables  
      REAL*8 A(N), B(N), C(N), D(N)
      REAL*8 TEMP, TEMP2, TEMP3, TEMP4
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             TEMP  = A(I)
             TEMP2 = B(I)
             TEMP3 = C(I)
             TEMP4 = D(I)
             J = I
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.TEMP) THEN
                A(J) = A(J-INC)
                B(J) = B(J-INC)
                C(J) = C(J-INC)
                D(J) = D(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = TEMP
             B(J) = TEMP2
             C(J) = TEMP3
             D(J) = TEMP4
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SHELL5(N, A,B,C,D,E)
!
!     sorts integer vectors.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
!     INPUT        MEANING
!     -----        -------
!     A            unsorted list to sort on 
!     N            length of vectors
!
!     OUTPUT       MEANING 
!     ------       ------- 
!     A            sorted list 
!     B            sorted with A 
!     C            sorted with A
!     D            sorted with A
!     E            sorted with A
!
!.... formal variables  
      REAL*8 A(N), B(N), C(N), D(N), E(N)
      REAL*8 TEMP, TEMP2, TEMP3, TEMP4, TEMP5
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             TEMP  = A(I)
             TEMP2 = B(I)
             TEMP3 = C(I)
             TEMP4 = D(I)
             TEMP5 = E(I)
             J = I 
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.TEMP) THEN
                A(J) = A(J-INC)
                B(J) = B(J-INC)
                C(J) = C(J-INC)
                D(J) = D(J-INC)
                E(J) = E(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = TEMP
             B(J) = TEMP2
             C(J) = TEMP3
             D(J) = TEMP4
             E(J) = TEMP5
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SHELL6(N, A,B,C,D,E,F)
!
!     sorts integer vectors.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
!     INPUT        MEANING
!     -----        -------
!     A            unsorted list to sort on 
!     N            length of vectors
!
!     OUTPUT       MEANING 
!     ------       ------- 
!     A            sorted list 
!     B            sorted with A 
!     C            sorted with A
!     D            sorted with A
!     E            sorted with A
!     F            sorted with A
!
!.... formal variables  
      REAL*8 A(N), B(N), C(N), D(N), E(N), F(N)
      REAL*8 TEMP, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             TEMP  = A(I)
             TEMP2 = B(I)
             TEMP3 = C(I)
             TEMP4 = D(I)
             TEMP5 = E(I)
             TEMP6 = F(I)
             J = I
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.TEMP) THEN
                A(J) = A(J-INC)
                B(J) = B(J-INC)
                C(J) = C(J-INC)
                D(J) = D(J-INC)
                E(J) = E(J-INC)
                F(J) = F(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = TEMP
             B(J) = TEMP2
             C(J) = TEMP3
             D(J) = TEMP4
             E(J) = TEMP5
             F(J) = TEMP6
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SHELL9(N, A,B,C,D,E,F,G,H,O)
!
!     sorts integer vectors.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
!     INPUT        MEANING
!     -----        -------
!     A            unsorted list to sort on 
!     N            length of vectors
!
!     OUTPUT       MEANING 
!     ------       ------- 
!     A            sorted list 
!     B            sorted with A 
!     C            sorted with A
!     D            sorted with A
!     E            sorted with A
!
!.... formal variables  
      REAL*8 A(N), B(N), C(N), D(N), E(N), F(N), G(N), H(N), O(N)
      REAL*8 TEMP, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6, TEMP7, TEMP8, 
     ;       TEMP9 
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             TEMP  = A(I)
             TEMP2 = B(I)
             TEMP3 = C(I)
             TEMP4 = D(I)
             TEMP5 = E(I)
             TEMP6 = F(I)
             TEMP7 = G(I)
             TEMP8 = H(I)
             TEMP9 = O(I)
             J = I
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.TEMP) THEN
                A(J) = A(J-INC)
                B(J) = B(J-INC)
                C(J) = C(J-INC)
                D(J) = D(J-INC)
                E(J) = E(J-INC)
                F(J) = F(J-INC)
                G(J) = G(J-INC)
                H(J) = H(J-INC)
                O(J) = O(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = TEMP
             B(J) = TEMP2
             C(J) = TEMP3
             D(J) = TEMP4
             E(J) = TEMP5
             F(J) = TEMP6
             G(J) = TEMP7
             H(J) = TEMP8
             O(J) = TEMP9
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE ISHELL4(N, A,B,C,D)
!
!     sorts integer vectors.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
!     INPUT        MEANING
!     -----        -------
!     A            unsorted list to sort on 
!     N            length of vectors
!
!     OUTPUT       MEANING 
!     ------       ------- 
!     A            sorted list 
!     B            sorted with A 
!     C            sorted with A
!     D            sorted with A
!
!.... formal variables  
      INTEGER*4 A(N), B(N), C(N), D(N)
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             ITEMP  = A(I)
             ITEMP2 = B(I)
             ITEMP3 = C(I)
             ITEMP4 = D(I)
             J = I
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.ITEMP) THEN
                A(J) = A(J-INC)
                B(J) = B(J-INC)
                C(J) = C(J-INC)
                D(J) = D(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = ITEMP
             B(J) = ITEMP2
             C(J) = ITEMP3
             D(J) = ITEMP4
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE ISHELL5(N, A,B,C,D,E)
!
!     sorts integer vectors.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
!     INPUT        MEANING
!     -----        -------
!     A            unsorted list to sort on 
!     N            length of vectors
!
!     OUTPUT       MEANING 
!     ------       ------- 
!     A            sorted list 
!     B            sorted with A 
!     C            sorted with A
!     D            sorted with A
!     E            sorted with A
!
!.... formal variables  
      INTEGER*4 A(N), B(N), C(N), D(N), E(N)
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             ITEMP  = A(I)
             ITEMP2 = B(I)
             ITEMP3 = C(I)
             ITEMP4 = D(I)
             ITEMP5 = E(I)
             J = I 
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.ITEMP) THEN
                A(J) = A(J-INC)
                B(J) = B(J-INC)
                C(J) = C(J-INC)
                D(J) = D(J-INC)
                E(J) = E(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = ITEMP
             B(J) = ITEMP2
             C(J) = ITEMP3
             D(J) = ITEMP4
             E(J) = ITEMP5
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SHELL1(N, A) 
!
!     sorts integer vector A.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
      INTEGER*4, INTENT(IN) ::N 
      REAL*8, INTENT(INOUT) :: A(N)
      REAL*8 TEMP
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             TEMP  = A(I) 
             J = I   
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.TEMP) THEN 
                A(J) = A(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = TEMP
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SHELL2(N, A,B) 
!
!     sorts integer vector A.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
      INTEGER*4, INTENT(IN) :: N
      REAL*8, INTENT(INOUT) :: A(N),B(N)
      REAL*8 TEMP,TEMP2 
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             TEMP  = A(I) 
             TEMP2 = B(I) 
             J = I   
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.TEMP) THEN 
                A(J) = A(J-INC)
                B(J) = B(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = TEMP
             B(J) = TEMP2 
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SHELLR8I2(N, A,B) 
!
!     sorts integer vector A.  program from ch 8.2, shells method
!     for quicksorting, numerical recipes
!
      INTEGER*4, INTENT(IN) :: N
      REAL*8, INTENT(INOUT) :: A(N)
      INTEGER*4, INTENT(INOUT) :: B(N)
      REAL*8 TEMP
      INTEGER*4 TEMP2
!
!----------------------------------------------------------------------!
!
!.... determine starting increment
      INC=1
    1 INC=3*INC+1
      IF(INC.LE.N) GOTO 1
!
!.... loop over the partial sorts
    2 CONTINUE
          INC=INC/3
!........ outer loop of straight insertion
          DO 11 I=INC+1,N
             TEMP  = A(I) 
             TEMP2 = B(I) 
             J = I   
!
!........... inner loop of straight insertion
    3        IF (A(J-INC).GT.TEMP) THEN 
                A(J) = A(J-INC)
                B(J) = B(J-INC)
                J=J-INC
                IF (J.LE.INC) GOTO 4
             GOTO 3
             ENDIF
    4        A(J) = TEMP
             B(J) = TEMP2 
   11     CONTINUE
      IF (INC.GT.1) GOTO 2
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DHPSRT1(N,X)
!
!     Purpose:
!     This subroutine sorts the  array X using HeapSort.  It rearranges
!     the elements of Y and Z at the same time.
!
!     Parameters:
!     N = the length of the arrays (input).
!     X = the primary array to be used in sorting (input/output).
!
!     Noel M. Nachtigal
!     October 4, 1990
!
!**********************************************************************
!
      INTEGER N
      DOUBLE PRECISION X(N)
!
!     Local variables.
!
      INTEGER I, J, K, L
      DOUBLE PRECISION TMPX
!
      IF (N.LE.1) RETURN
!
      L = N / 2 + 1
      K = N
 10   IF (L.GT.1) THEN
         L = L - 1
         TMPX = X(L)
      ELSE
         TMPX = X(K)
         X(K) = X(1)
         K = K - 1
         IF (K.LE.1) THEN
            X(1) = TMPX
            RETURN
         END IF
      END IF
      I = L
      J = L + L
 20   IF (J.LE.K) THEN
         IF (J.LT.K) THEN
            IF (X(J).LT.X(J+1)) J = J + 1
         END IF
         IF (TMPX.LT.X(J)) THEN
            X(I) = X(J)
            I    = J
            J    = J + J
         ELSE
            J = K + 1
         END IF
         GO TO 20
      END IF
      X(I) = TMPX
      GO TO 10
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DHPSRT2(N,X,Y)
!
!     Purpose:
!     This subroutine sorts the  array X using HeapSort.  It rearranges
!     the elements of Y and Z at the same time.
!
!     Parameters:
!     N = the length of the arrays (input).
!     X = the primary array to be used in sorting (input/output).
!     Y = another array, sorted in the same order as X (input/output).
!
!     Noel M. Nachtigal
!     October 4, 1990
!
!**********************************************************************
!
      INTEGER N
      DOUBLE PRECISION X(N), Y(N)
!
!     Local variables.
!
      INTEGER I, J, K, L
      DOUBLE PRECISION TMPX, TMPY
!
      IF (N.LE.1) RETURN
!
      L = N / 2 + 1
      K = N
 10   IF (L.GT.1) THEN
         L = L - 1
         TMPX = X(L)
         TMPY = Y(L)
      ELSE
         TMPX = X(K)
         TMPY = Y(K)
         X(K) = X(1)
         Y(K) = Y(1)
         K = K - 1
         IF (K.LE.1) THEN
            X(1) = TMPX
            Y(1) = TMPY
            RETURN
         END IF
      END IF
      I = L
      J = L + L
 20   IF (J.LE.K) THEN
         IF (J.LT.K) THEN
            IF (X(J).LT.X(J+1)) J = J + 1
         END IF
         IF (TMPX.LT.X(J)) THEN
            X(I) = X(J)
            Y(I) = Y(J)
            I    = J
            J    = J + J
         ELSE
            J = K + 1
         END IF
         GO TO 20
      END IF
      X(I) = TMPX
      Y(I) = TMPY
      GO TO 10
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DHPSRT3(N,X,Y,Z)
!
!     Purpose:
!     This subroutine sorts the  array X using HeapSort.  It rearranges
!     the elements of Y and Z at the same time.
!
!     Parameters:
!     N = the length of the arrays (input).
!     X = the primary array to be used in sorting (input/output).
!     Y = another array, sorted in the same order as X (input/output).
!     Z = another array, sorted in the same order as X (input/output).
!     Noel M. Nachtigal
!     October 4, 1990
!
!**********************************************************************
!
      INTEGER N
      DOUBLE PRECISION X(N), Y(N), Z(N)
!
!     Local variables.
!
      INTEGER I, J, K, L
      DOUBLE PRECISION TMPX, TMPY, TMPZ
!
      IF (N.LE.1) RETURN
!
      L = N / 2 + 1 
      K = N 
 10   IF (L.GT.1) THEN
         L = L - 1 
         TMPX = X(L)
         TMPY = Y(L)
         TMPZ = Z(L)
      ELSE
         TMPX = X(K)
         TMPY = Y(K)
         TMPZ = Z(K)
         X(K) = X(1)
         Y(K) = Y(1)
         Z(K) = Z(1)
         K = K - 1 
         IF (K.LE.1) THEN
            X(1) = TMPX
            Y(1) = TMPY
            Z(1) = TMPZ
            RETURN
         END IF
      END IF
      I = L
      J = L + L
 20   IF (J.LE.K) THEN 
         IF (J.LT.K) THEN 
            IF (X(J).LT.X(J+1)) J = J + 1
         END IF
         IF (TMPX.LT.X(J)) THEN 
            X(I) = X(J) 
            Y(I) = Y(J) 
            Z(I) = Z(J)
            I    = J  
            J    = J + J
         ELSE 
            J = K + 1
         END IF
         GO TO 20
      END IF
      X(I) = TMPX 
      Y(I) = TMPY 
      Z(I) = TMPZ
      GO TO 10
      END  
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE IHPSRT3(N,IX,IY,IZ)
!
!     Purpose:
!     This subroutine sorts the  array X using HeapSort.  It rearranges
!     the elements of Y and Z at the same time.
!
!     Parameters:
!     N = the length of the arrays (input).
!     X = the primary array to be used in sorting (input/output).
!     Y = another array, sorted in the same order as X (input/output).
!
!     Noel M. Nachtigal
!     October 4, 1990
!
!**********************************************************************
!
      IMPLICIT NONE
      INTEGER N
      INTEGER*4 IX(N), IY(N), IZ(N)
!
!     Local variables.
!
      INTEGER*4 I, J, K, L
      INTEGER*4 ITMPX, ITMPY, ITMPZ
!
      IF (N.LE.1) RETURN
!
      L = N / 2 + 1 
      K = N 
 10   IF (L.GT.1) THEN
         L = L - 1 
         ITMPX = IX(L)
         ITMPY = IY(L)
         ITMPZ = IZ(L)
      ELSE
         ITMPX = IX(K)
         ITMPY = IY(K)
         ITMPZ = IZ(K)
         IX(K) = IX(1)
         IY(K) = IY(1)
         IZ(K) = IZ(1)
         K = K - 1 
         IF (K.LE.1) THEN
            IX(1) = ITMPX
            IY(1) = ITMPY
            IZ(1) = ITMPZ
            RETURN
         END IF
      END IF
      I = L
      J = L + L
 20   IF (J.LE.K) THEN
         IF (J.LT.K) THEN
            IF (IX(J).LT.IX(J+1)) J = J + 1
         END IF
         IF (ITMPX.LT.IX(J)) THEN
            IX(I) = IX(J)
            IY(I) = IY(J)
            IZ(I) = IZ(J)
            I    = J
            J    = J + J
         ELSE
            J = K + 1
         END IF
         GO TO 20
      END IF
      IX(I) = ITMPX
      IY(I) = ITMPY
      IZ(I) = ITMPZ
      GO TO 10
      END
!                                                                      !
!======================================================================!
!                                                                      !
      INTEGER*4 FUNCTION IBSECT(N,IX,IVEC) 
! 
!     Searches an integer array IVEC for index IX via bisection
!     Numerical Recipes pg 111.  
      INTEGER*4 IVEC(N) 
      IBSECT =-1 !flag an error
! 
!.... get the ends out of the way first
      IF (IX.EQ.IVEC(1)) THEN
         IBSECT = 1 
         RETURN
      ENDIF
      IF (IX.EQ.IVEC(N)) THEN
         IBSECT = N !n.r. has n-1, not appropriate here 
         RETURN
      ENDIF 
! 
!.... begin bisection search
      JL = 0 
      JU = N + 1 
    1 IF (JU - JL.GT.1) THEN 
         JM = (JU + JL)/2 !Midoint
         IF ( (IVEC(N).GE.IVEC(1)).EQV.(IX.GE.IVEC(JM)) ) THEN 
            JL = JM  
         ELSE
            JU = JM
         ENDIF
         GOTO 1 
      ENDIF 
      IF (IX.EQ.IVEC(JL)) THEN
         IBSECT = JL
      ELSE
         PRINT *, 'ibsect: Error could not find index',IX
         IBSECT =-1 
      ENDIF
      RETURN 
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      INTEGER*4 FUNCTION IBSECT8(N,JOB,TOL,X,VEC)
! 
!     Searches an integer array IVEC for index IX via bisection
!     Numerical Recipes pg 111.  
      REAL*8 VEC(N),X,TOL
      IBSECT8 =-1 !flag an error
! 
!.... get the ends out of the way first
      !IF (X.EQ.VEC(1)) THEN
      IF (DABS(X - VEC(1)).LT.TOL) THEN
         IBSECT8 = 1 
         RETURN !definitely done
      ENDIF
      !IF (X.EQ.VEC(N)) THEN
      IF (DABS(X - VEC(N)).LT.TOL) THEN
         IBSECT8 = N !n.r. has n-1, not appropriate here 
         GOTO 50 !may still need to count backwards
      ENDIF 
! 
!.... begin bisection search
      JL = 0 
      JU = N + 1 
    1 IF (JU - JL.GT.1) THEN 
         JM = (JU + JL)/2 !Midoint
         IF ( (VEC(N).GE.VEC(1)).EQV.(X.GE.VEC(JM)) ) THEN
            JL = JM
         ELSE
            JU = JM
         ENDIF
         GOTO 1 
      ENDIF 
      IF (JL.LE.0) THEN
         IBSECT8 =-1 
         RETURN
      ENDIF
!
!.... bracketing is approximate
      IF (DABS(VEC(JL) - X).LT.TOL) THEN
         IBSECT8 = JL
      ELSEIF (DABS(VEC(JM) - X).LT.TOL) THEN
         IBSECT8 = JM
      ELSEIF (DABS(VEC(JU) - X).LT.TOL) THEN
         IBSECT8 = JU
      ELSE
         IBSECT8 =-1 
      ENDIF
! 
!.... possibly count backwards if list has repeats
   50 CONTINUE
      IF (JOB.EQ.1) THEN
         DO 40 ILOC=IBSECT8,2,-1 
            IF (DABS(VEC(ILOC) - VEC(ILOC-1)).GT.TOL) THEN
               IBSECT8 = ILOC
               RETURN
            ENDIF 
  40     CONTINUE 
         PRINT *, 'ibsect8: shouldnt be here'
      ENDIF
      RETURN 
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      INTEGER*4 FUNCTION IBSECT2(N,IX,IVEC)
! 
!     Searches an integer array IVEC for index IX via bisection
!     Numerical Recipes pg 111.  
      INTEGER*4 IVEC(N)
      IBSECT2 =-1 !flag an error
! 
!.... get the ends out of the way first
      IF (IX.EQ.IVEC(1)) THEN
         IBSECT2 = 1 
         RETURN
      ENDIF
      IF (IX.EQ.IVEC(N)) THEN
         IBSECT2 = N !n.r. has n-1, not appropriate here 
         RETURN
      ENDIF
! 
!.... begin bisection search
      JL = 0 
      JU = N + 1 
    1 IF (JU - JL.GT.1) THEN
         JM = (JU + JL)/2 !Midoint
         IF ( (IVEC(N).GE.IVEC(1)).EQV.(IX.GE.IVEC(JM)) ) THEN
            JL = JM
         ELSE
            JU = JM
         ENDIF
         GOTO 1
      ENDIF
      IF (JL.LE.0) THEN
         IBSECT2 =-1 
         RETURN
      ENDIF
      IF (IX.EQ.IVEC(JL)) THEN
         IBSECT2 = JL
      ELSE
         IBSECT2 =-1 
      ENDIF
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      INTEGER*4 FUNCTION LOCATE8(N,THRESH, X,VEC) 
! 
!     Searches an array VEC for index X via bisection.
!     Numerical Recipes pg 111.  
      REAL*8 VEC(N),THRESH,X
      INTEGER*4 N
      INTEGER*4 JM,JL,JU
      LOCATE8 =-1 !flag an error
! 
!.... get the ends out of the way first
      !IF (X.EQ.IVEC(1)) THEN
      IF (DABS(VEC(1)-X).LT.THRESH) THEN
         LOCATE8 = 1 
         RETURN
      ENDIF
      !IF (IX.EQ.IVEC(N)) THEN
      IF (DABS(VEC(N)-X).LT.THRESH) THEN
         LOCATE8 = N !n.r. has n-1, not appropriate here 
         RETURN
      ENDIF 
! 
!.... begin bisection search
      JL = 0 
      JU = N + 1 
    1 IF (JU - JL.GT.1) THEN 
         JM = (JU + JL)/2 !Midoint
         IF ( (VEC(N).GE.VEC(1)).EQV.(X.GE.VEC(JM)) ) THEN 
            JL = JM  
         ELSE
            JU = JM
         ENDIF
         GOTO 1 
      ENDIF 
      !IF (X.EQ.VEC(JL)) THEN
      IF (DABS(VEC(JL)-X).LT.THRESH) THEN
         LOCATE8 = JL
      ELSEIF (DABS(VEC(JU)-X).LT.THRESH) THEN
         LOCATE8 = JU
      ELSE
         PRINT *, 'locate8: Error could not find index',X
         LOCATE8 =-1 
      ENDIF
      RETURN 
      END 

