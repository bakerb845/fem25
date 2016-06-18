      logical*4 lex
      integer*4 b
      real*4, allocatable :: a(:,:) 
      integer*4, allocatable :: ipiv(:)

      n = 4000
      CALL INIT_RANDOM_SEED() !create a random seed
      write(*,*) 'initializing matrix...' 
      allocate(a(n,n))
      allocate(ipiv(n)) 
      do 1 i=1,n
         do 2 j=i,n
            call random_number(t)
            if (i.eq.j) then
               a(i,j) = t + 2.0 
            else
               a(i,j) = T 
               a(j,i) = 0.0 
            endif
    2    continue
    1 continue
      write(*,*) 'factoring matrix...'
      print *, n,n
      call sgetrf(n,n,a,n,ipiv,info)
      print *, 'info',info
      STOP
      END

      SUBROUTINE INIT_RANDOM_SEED()
      INTEGER*4 I, N, CLOCK
      INTEGER*4, ALLOCATABLE :: SEED(:)
      CALL RANDOM_SEED(SIZE = N)
      ALLOCATE(SEED(N))
      CALL SYSTEM_CLOCK(COUNT=CLOCK)
      SEED = CLOCK + 37 * (/ (I - 1, I = 1, N) /)
      CALL RANDOM_SEED(PUT = SEED)
      DEALLOCATE(SEED)
      RETURN
      END 
