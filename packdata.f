!----------------------------------------------------------------------!
!                functions for packing/unpacking data                  !
!----------------------------------------------------------------------!
      CHARACTER(4) FUNCTION PACKI4(LSWAP,IN4)
      IMPLICIT NONE
      INTEGER*4, INTENT(IN) :: IN4 
      LOGICAL*4, INTENT(IN) :: LSWAP
      CHARACTER(4) C4
      INTEGER*4 I,I4
      EQUIVALENCE(C4,I4)
      I4 = IN4 
      IF (LSWAP) THEN
         DO 1 I=1,4
            PACKI4(I:I) = C4(5-I:5-I)
    1    CONTINUE
      ELSE
         DO 2 I=1,4
            PACKI4(I:I) = C4(I:I)
    2    CONTINUE
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      CHARACTER(4) FUNCTION PACKR4(LSWAP,RIN4)
      IMPLICIT NONE
      REAL*4, INTENT(IN) :: RIN4
      LOGICAL*4, INTENT(IN) :: LSWAP
      CHARACTER(4) C4
      REAL*4 R4
      INTEGER*4 I
      EQUIVALENCE(C4,R4)
      R4 = RIN4
      IF (LSWAP) THEN
         DO 1 I=1,4
            PACKR4(I:I) = C4(5-I:5-I)
    1    CONTINUE
      ELSE
         DO 2 I=1,4
            PACKR4(I:I) = C4(I:I)
    2    CONTINUE
      ENDIF
      RETURN    
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      CHARACTER(8) FUNCTION PACKR8(LSWAP,RIN8)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: RIN8
      LOGICAL*4, INTENT(IN) :: LSWAP
      CHARACTER(8) C8
      REAL*8 R8
      INTEGER*4 I
      EQUIVALENCE(C8,R8)
      R8 = RIN8
      IF (LSWAP) THEN
         DO 1 I=1,8
            PACKR8(I:I) = C8(9-I:9-I)
    1    CONTINUE
      ELSE
         DO 2 I=1,8
            PACKR8(I:I) = C8(I:I)
    2    CONTINUE
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      INTEGER*4 FUNCTION UNPACKI4(LSWAP,CIN4)
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: CIN4(4)
      LOGICAL*4, INTENT(IN) :: LSWAP
      CHARACTER(4) C4
      INTEGER*4 I,I4
      EQUIVALENCE(C4,I4)
      IF (LSWAP) THEN
         DO 1 I=1,4
            C4(I:I) = CIN4(5-I)
    1    CONTINUE
      ELSE
         DO 2 I=1,4
            C4(I:I) = CIN4(I)
    2    CONTINUE
      ENDIF
      UNPACKI4 = I4
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*4 FUNCTION UNPACKR4(LSWAP,CIN4)
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: CIN4(4)
      LOGICAL*4, INTENT(IN) :: LSWAP
      REAL*4 R4
      INTEGER*4 I
      CHARACTER(4) C4
      EQUIVALENCE(C4,R4)
      IF (LSWAP) THEN
         DO 1 I=1,4
            C4(I:I) = CIN4(5-I)
    1    CONTINUE
      ELSE
         DO 2 I=1,4
            C4(I:I) = CIN4(I)
    2    CONTINUE
      ENDIF
      UNPACKR4 = R4
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*8 FUNCTION UNPACKR8(LSWAP,CIN8)
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: CIN8(8)
      LOGICAL*4, INTENT(IN) :: LSWAP
      CHARACTER(8) C8
      REAL*8 R8
      INTEGER*4 I
      EQUIVALENCE(C8,R8)
      IF (LSWAP) THEN
         DO 1 I=1,8
            C8(I:I) = CIN8(9-I)
    1    CONTINUE
      ELSE
         DO 2 I=1,8
            C8(I:I) = CIN8(I)
    2    CONTINUE
      ENDIF
      UNPACKR8 = R8
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      FUNCTION ENDIAN()
! This file: http://ftp.aset.psu.edu/pub/ger/fortran/hdk/endian.f90
!
! by Code Tuning co-guide, 1998 Lahey Fortran Users' Conference
! Check what endian this program is running on. 
! 
!     RESULTS    MEANING 
!     -------    ------- 
!     endian     = 0 -> little endian 
!                = 1 -> big endian 
!                =-1 -> mixed endian 
! 
!----------------------------------------------------------------------!
! 
      INTEGER*4 ENDIAN, IEND
      INTEGER*4 I, ASCII0, ASCII1, ASCII2, ASCII3
      PARAMETER(ASCII0 = 48, ASCII1 = 49, ASCII2 = 50, ASCII3 = 51)
      COMMON // I
      I = ASCII0 + ASCII1*256 + ASCII2*(256**2) + ASCII3*(256**3)
      CALL SUBEND(IEND)
      ENDIAN = IEND
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SUBEND(IEND)
      INTEGER*4 IEND
      CHARACTER(4) I
      COMMON // I
      IF (I == '0123') THEN !little endian
        IEND = 0
        RETURN
      ENDIF
      IF (I == '3210') THEN !big endian
        IEND = 1
        RETURN
      ENDIF
      IEND = -1 !mixed
      RETURN
      END

