      LOGICAL*4 FUNCTION LISDIR(DIRNM)
!     Based on the INTEL or GNU calling sequence, determines
!     if a directory name, DIRNM, exists
      CHARACTER(*), INTENT(IN) :: DIRNM
      LOGICAL*4 LEX
#ifdef INTEL
      INQUIRE(DIRECTORY=TRIM(ADJUSTL(DIRNM)),EXIST=LEX)
#else
      INQUIRE(FILE=TRIM(ADJUSTL(DIRNM)),EXIST=LEX) 
#endif
      LISDIR = LEX
      RETURN
      END
