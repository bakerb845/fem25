C=====================================================================C 
C                   a few string manipulation utilities               C
C=====================================================================C

      SUBROUTINE NULLS(LEN,ALINE)
C
C     nulls out a string
      CHARACTER*1 ALINE(LEN)
      DO 1 I=1,LEN
         ALINE(I) = ' '
    1 CONTINUE
      RETURN
      END
C
C=====================================================================C
C

      SUBROUTINE DBLANK(INLEN, LINE, LEN)
C
C     removes blanks from a string, code modified from steves original
C
C     INPUT      MEANING
C     -----      -------
C     INLEN      input length of string line
C
C     I/O        MEANING
C     ---        -------
C     LINE       string to remove blanks from
C
C     OUTPUT     MEANING
C     ------     -------
C     LEN        length of line with removed blanks
C
C     variable declarations
      CHARACTER*1 LINE(*)
C
C------------------------------------------------------------------C
C
C.... look for first character (possible quick return)
      DO 1 I=1,INLEN
         J=INLEN+1-I
         IF (LINE(J).NE.' ') GOTO 60
    1 CONTINUE
      LEN=0
      RETURN
C.... find first blank
   60 DO 2 I=1,J
C....... if there's a blank shift everything in ahead back 1
         IF (LINE(I).EQ.' ') THEN
            DO 3 K=I,J+1
               LINE(K)=LINE(K+1)
    3       CONTINUE
         ENDIF
    2 CONTINUE
C.... look for the first character
      DO 4 I=1,INLEN
         J=INLEN+1-I
         IF (LINE(J).NE.' ') GOTO 61
    4 CONTINUE
C.... find next blank and remove it
   61 LEN=J
      DO 5 I=1,J
         IF (LINE(I).EQ.' ') GOTO 60
    5 CONTINUE
      RETURN
      END

C
C====================================================================C
C
      SUBROUTINE COPYS(LEN,STRING,ALINE)
C
C     copies string of len into aline
      CHARACTER*1 STRING(LEN), ALINE(LEN)
      DO 1 I=1,LEN
         ALINE(I) = STRING(I)
    1 CONTINUE
      RETURN
      END

C
C====================================================================C
C

      SUBROUTINE DTB(STRING,LS, N,RC)
C
C     converts the right justified numeral in STRING to an integer
C     in N (Decimal to Binary).  code from mike kupferschmid
C
C     INPUT    MEANING
C     -----    -------
C     LS       number of characters in the string
C     STRING   the string of numerals to be converted
C
C     OUTPUT   MEANING
C     ------   -------
C     N        the integer returned
C     RC       return code; 0 => ok, 1 => error
C
C.... formal parameters
      CHARACTER*1 STRING(LS)
      INTEGER*4 RC
C.... local variables
      INTEGER*4 SIGN
      CHARACTER*1 NUMERL(10)/'0','1','2','3','4',
     ;                       '5','6','7','8','9'/
C
C--------------------------------------------------------------------C
C
C.... if the string is empty, something is wrong
      RC = 1
      N = 0
      IF (LS.LE.0) RETURN
C.... skip leading blanks
      DO 1 K=1,LS
         KFDGT = K
         IF (STRING(K).NE.' ') GOTO 2
    1 CONTINUE
C.... the string is all blanks
      RETURN
C.... determine the sign and the index of the first digit
    2 IF (STRING(KFDGT).EQ.'-') THEN
         SIGN = -1
         KFDGT = KFDGT + 1
         GOTO 3
      ENDIF
      SIGN = +1
      IF (STRING(KFDGT).EQ.'+') KFDGT = KFDGT + 1
    3 IF (KFDGT.GT.LS) RETURN
C.... examine the digits and construct the value of the number
      DO 4 K=KFDGT,LS
C....... add in the positional value of the digit
         DO 5 L=1,10
            IF (STRING(K).NE.NUMERL(L)) GOTO 5
            N = N + (L - 1)*10**(LS-K)
            GOTO 4
    5    CONTINUE
C....... the character is not a digit
         N = 0
         RETURN
    4 CONTINUE
C.... give the value the correct sign
      N = SIGN*N
      RC = 0
      RETURN
      END
C
C=====================================================================C
C
      SUBROUTINE CENTRA(ZULU,IMIN,IMAX)
C 
C     returns positions of the first and last non-blank characters of 
C     zulu in imin and imax 
C
      CHARACTER *80 ZULU
      DO 1 I=1,80
         IF (ZULU(I:I).NE.' ') GOTO 11
    1 CONTINUE
   11 IMIN=I
      DO 2 I=1,80
         IF (ZULU(81-I:81-I).NE.' ') GOTO 22
    2 CONTINUE
   22 IMAX=81-I
      RETURN
      END 
C 
C======================================================================C
C 
      INTEGER FUNCTION ITRMLN(STRING)
c
c This function returns the length of a character string
c with all trailing blanks removed. It also appends a null
C after the last character.
c It reads the character string backwards until a 
c character is encountered that is neither a null nor a space.
c
      CHARACTER*(*) STRING
      INTEGER L,I 
c
c check the length of the character string.
c
      L = LEN(STRING)
c
      IF (L.LE.0) THEN
         ITRMLN=0
         GOTO 30
      END IF
C
      I = L 
10    IF (STRING(I:I) .NE. ' '.AND.STRING(I:I).NE.CHAR(0)) GO TO 20
        I = I - 1 
        IF (I.GT.0) THEN
           GO TO 10
        ELSE
          GO TO 20
        END IF
C
20    CONTINUE
      IF (I.GT.0) THEN
         ITRMLN = I 
      ELSE
         ITRMLN = 0 
      END IF
30    CONTINUE
      IF (ITRMLN.LT.L) STRING(ITRMLN+1:ITRMLN+1) = CHAR(0)
      RETURN
      END
C 
C======================================================================C
C 
      FUNCTION NOBLNK(STRING)
C
C Deletes all spaces:
C
C      IMPLICIT NONE
      INTEGER       I1, I2, LTRM, I
      INTEGER       ITRMLN
      CHARACTER*(*) STRING
      CHARACTER*(*) NOBLNK
      NOBLNK = STRING
      LTRM = ITRMLN(NOBLNK)
      DO 8500 I = 1, LTRM
         IF (NOBLNK(I:I).NE.' ') GOTO 8505
 8500 CONTINUE
 8505 NOBLNK(1:) = NOBLNK(I:)
 9000 CONTINUE
         I1 = INDEX(NOBLNK,' ')
         LTRM = ITRMLN(NOBLNK)
         IF (I1.EQ.0.OR.I1.GT.LTRM) GOTO 9050
         I2 = I1 + 1
 9010    IF (NOBLNK(I2:I2).EQ.' ') THEN
            I2 = I2 + 1
            GOTO 9010
         END IF
         I2 = MIN(I2,LTRM)
         NOBLNK = NOBLNK(1:I1-1)//NOBLNK(I2:LTRM)
         GOTO 9000
 9050 CONTINUE
      END
C 
C======================================================================C
C 
      INTEGER FUNCTION CHAR_LEN(STRING)

C     this function returns the length of a character string
C     with all trailing blanks removed. it also appends a null 
c     after the last character.
C     it reads the character string backwards until a 
C     character is encountered that is neither a null nor a space.

      CHARACTER*(*) STRING
      INTEGER I, ILOOP

      DO 1 I = LEN(STRING), 1, -1
         ILOOP = I
         IF (STRING(I:I) .NE. ' ' .AND. 
     ;       STRING(I:I) .NE. CHAR(0) ) GOTO 10
    1 CONTINUE 

      STRING(1:1) = CHAR(0)     
      CHAR_LEN = 0  
      RETURN

   10 CONTINUE
      IF (ILOOP .LT. LEN(STRING)) THEN 
         STRING(ILOOP+1:ILOOP+1) = CHAR(0)
      ENDIF
      CHAR_LEN = ILOOP
      RETURN
      END

