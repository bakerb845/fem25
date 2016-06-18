      SUBROUTINE CSF2D(MDIM,NGNOD,XI,ETA, SF, DSF, IERR)
! 
!     the standard element  
!
!     1 ^   4 - - - 7 - - - 3
!       |   |       |       |
!       |   |       |       |
!     0 xi  8 - - - 9 - - - 6
!       |   |       |       |
!       |   |       |       |
!    -1     1 - - - 5 - - - 2
!
!            ----- eta ---->
!          -1       0       1
! 
!     INPUT    MEANING 
!     -----    ------- 
!     ETA      eta location 
!     MDIM     max number of dimensions 
!     NGNOD    number of global nodes on element 
!     XI       xi location 
!   
!     OUTPUT   MEANING 
!     ------   ------- 
!     DSF      derivative of shape fns at xi / eta 
!     IERR     error flag 
!     SF       shape functions at xi / eta   
! 
!.... variable declarations 
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: XI, ETA 
      INTEGER*4, INTENT(IN) :: MDIM, NGNOD 
      REAL*8, INTENT(OUT) :: SF(NGNOD), DSF(MDIM,*) 
      INTEGER*4, INTENT(OUT) :: IERR 
      INTEGER*4 I 
      REAL*8 SUM1, SUM2, SUM3, EPS 
      PARAMETER(EPS = 2.2204D-16)
! 
!----------------------------------------------------------------------!
!     
!.... 4 noded element 
      IERR = 0 
      IF (NGNOD.EQ.4) THEN 
!
!....... hat shape functions 
         SF(1) = .25D0*(1.D0 - XI)*(1.D0 - ETA)
         SF(2) = .25D0*(1.D0 + XI)*(1.D0 - ETA)
         SF(3) = .25D0*(1.D0 + XI)*(1.D0 + ETA)
         SF(4) = .25D0*(1.D0 - XI)*(1.D0 + ETA)
         IF (DABS(SF(1)).LT.EPS) SF(1) = 0.D0
         IF (DABS(SF(2)).LT.EPS) SF(2) = 0.D0
         IF (DABS(SF(3)).LT.EPS) SF(3) = 0.D0
         IF (DABS(SF(4)).LT.EPS) SF(4) = 0.D0
!
!....... their derivatives w.r.t. xi 
         DSF(1,1) = .25D0*(-1.D0)*(1.D0 - ETA)
         DSF(1,2) = .25D0*( 1.D0)*(1.D0 - ETA)
         DSF(1,3) = .25D0*( 1.D0)*(1.D0 + ETA)
         DSF(1,4) = .25D0*(-1.D0)*(1.D0 + ETA)
! 
!....... their derivatives w.r.t. eta 
         DSF(2,1) = .25D0*(1.D0 - XI)*(-1.D0)
         DSF(2,2) = .25D0*(1.D0 + XI)*(-1.D0)
         DSF(2,3) = .25D0*(1.D0 + XI)*( 1.D0)
         DSF(2,4) = .25D0*(1.D0 - XI)*( 1.D0)
! 
!.... 9 noded element 
      ELSEIF(NGNOD.EQ.9) THEN
! 
!....... bubble functions 
         SF(1) = .25D0*XI*(XI - 1.D0)*ETA*(ETA - 1.D0)
         SF(2) = .25D0*XI*(XI + 1.D0)*ETA*(ETA - 1.D0)
         SF(3) = .25D0*XI*(XI + 1.D0)*ETA*(ETA + 1.D0)
         SF(4) = .25D0*XI*(XI - 1.D0)*ETA*(ETA + 1.D0)
         SF(5) = .5D0*(1.D0 - XI**2)*ETA*(ETA - 1.D0)
         SF(7) = .5D0*(1.D0 - XI**2)*ETA*(ETA + 1.D0)
         SF(6) = .5D0*XI*(XI + 1.D0)*(1.D0 - ETA**2)
         SF(8) = .5D0*XI*(XI - 1.D0)*(1.D0 - ETA**2)
         SF(9) =  (1.D0 - XI**2)*(1.D0 - ETA**2)
! 
!....... their derivatives w.r.t xi 
         DSF(1,1) = .5D0*(XI - .5D0)*ETA*(ETA - 1.D0)
         DSF(1,2) = .5D0*(XI + .5D0)*ETA*(ETA - 1.D0)
         DSF(1,3) = .5D0*(XI + .5D0)*ETA*(ETA + 1.D0)
         DSF(1,4) = .5D0*(XI - .5D0)*ETA*(ETA + 1.D0)
         DSF(1,5) = -XI*ETA*(ETA - 1.D0)
         DSF(1,7) = -XI*ETA*(ETA + 1.D0)
         DSF(1,6) =  (XI + .5D0) * (1.D0 - ETA**2)
         DSF(1,8) =  (XI - .5D0) * (1.D0 - ETA**2)
         DSF(1,9) = -2.D0*XI*(1.D0 - ETA**2)
! 
!....... their derivatives w.r.t. eta 
         DSF(2,1) = .5D0*XI*(XI - 1.D0)*(ETA - .5D0)
         DSF(2,2) = .5D0*XI*(XI + 1.D0)*(ETA - .5D0)
         DSF(2,3) = .5D0*XI*(XI + 1.D0)*(ETA + .5D0)
         DSF(2,4) = .5D0*XI*(XI - 1.D0)*(ETA + .5D0)
         DSF(2,5) =  (1.D0 - XI**2)*(ETA - .5D0)
         DSF(2,7) =  (1.D0 - XI**2)*(ETA + .5D0)
         DSF(2,6) =  XI*(XI + 1.D0)*(-ETA)
         DSF(2,8) =  XI*(XI - 1.D0)*(-ETA)
         DSF(2,9) = -2.D0*ETA*(1.D0 - XI**2)
      ELSE !no serendipity nodes yet and pry ever
         WRITE(*,*) 'csf2d: only 4 or 9 nodes supported'
         IERR = 1
         RETURN
      ENDIF
! 
!.... i like this simple check in komatisch and tromp 
      SUM1 = 0.D0
      SUM2 = 0.D0
      SUM3 = 0.D0
      DO 1 I=1,NGNOD
         IF (DABS(SF(I)).LT.EPS) SF(I) = 0.D0
         IF (DABS(DSF(1,I)).LT.EPS) DSF(1,I) = 0.D0
         IF (DABS(DSF(2,I)).LT.EPS) DSF(2,I) = 0.D0
         SUM1 = SUM1 + SF(I)
         SUM2 = SUM2 + DSF(1,I)
         SUM3 = SUM3 + DSF(2,I)
    1 CONTINUE
      IF (DABS(SUM2).GT.EPS) THEN
         WRITE(*,*) 'csf2d: Error in xi shape functions'
         IERR = 1
         RETURN
      ENDIF
      IF (DABS(SUM3).GT.EPS) THEN
         WRITE(*,*) 'csf2d: Error in eta shape functions'
         IERR = 1
         RETURN
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CSF2DND(NGNOD,XI,ETA, SF, IERR)
! 
!     Calculates 2D shape functions for interpolation 
!
!     1 ^   4 - - - 7 - - - 3
!       |   |       |       |
!       |   |       |       |
!     0 xi  8 - - - 9 - - - 6
!       |   |       |       |
!       |   |       |       |
!    -1     1 - - - 5 - - - 2
!
!            ----- eta ---->
!          -1       0       1
! 
!     INPUT    MEANING 
!     -----    ------- 
!     ETA      eta location 
!     NGNOD    number of global nodes on element 
!     XI       xi location 
!   
!     OUTPUT   MEANING 
!     ------   ------- 
!     IERR     error flag 
!     SF       shape functions at xi / eta   
! 
!.... variable declarations 
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: XI, ETA
      INTEGER*4, INTENT(IN) :: NGNOD
      REAL*8, INTENT(OUT) :: SF(NGNOD)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 I 
      REAL*8 EPS 
      PARAMETER(EPS = 2.2204D-16)
! 
!----------------------------------------------------------------------!
!     
!.... 4 noded element 
      IERR = 0 
      IF (NGNOD.EQ.4) THEN !4 node element; hat shape fns 
         SF(1) = .25D0*(1.D0 - XI)*(1.D0 - ETA)
         SF(2) = .25D0*(1.D0 + XI)*(1.D0 - ETA)
         SF(3) = .25D0*(1.D0 + XI)*(1.D0 + ETA)
         SF(4) = .25D0*(1.D0 - XI)*(1.D0 + ETA)
      ELSEIF(NGNOD.EQ.9) THEN !9 node element; bubble shape fns
         SF(1) = .25D0*XI*(XI - 1.D0)*ETA*(ETA - 1.D0)
         SF(2) = .25D0*XI*(XI + 1.D0)*ETA*(ETA - 1.D0)
         SF(3) = .25D0*XI*(XI + 1.D0)*ETA*(ETA + 1.D0)
         SF(4) = .25D0*XI*(XI - 1.D0)*ETA*(ETA + 1.D0)
         SF(5) = .5D0*(1.D0 - XI**2)*ETA*(ETA - 1.D0)
         SF(7) = .5D0*(1.D0 - XI**2)*ETA*(ETA + 1.D0)
         SF(6) = .5D0*XI*(XI + 1.D0)*(1.D0 - ETA**2)
         SF(8) = .5D0*XI*(XI - 1.D0)*(1.D0 - ETA**2)
         SF(9) =  (1.D0 - XI**2)*(1.D0 - ETA**2)
      ELSE !no serenditity nodes as of now
         WRITE(*,*) 'csf2dnd: only 4 or 9 nodes supported'
         IERR = 1
         RETURN
      ENDIF
      DO 1 I=1,NGNOD
         IF (DABS(SF(I)).LT.EPS) SF(I) = 0.D0
    1 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CSFT2D(MDIM,NGNOD, R,S, SF,DSF)
! 
!     Triangle shape functions from Hughes pg 167-168
!       
!       2         (1,3,1) x
!       |\                |\
!       | \               | \
!       5  4      (1,2,1) x  x (2,2,1)
!       |   \             |   \
!       |    \            |    \
!       3 -6- 1   (1,1,3) x -x- x (3,1,1)
!                         (2,1,2)
! 
!     INPUT      MEANING
!     -----      ------- 
!     MDIM       leading dimension for DSF
!     NGNOD      number of anchor nodes (3, or 6, 3 default)
!     R          r local coordinate 
!     S          s local coordinate 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     DSF        derivatives of shape functions in r,s at (r,s)
!     SF         shape functions evaluated at (r,s) 
! 
!.... variable declarations
      INTEGER*4, INTENT(IN) :: MDIM,NGNOD 
      REAL*8, INTENT(IN) :: R, S
      REAL*8, INTENT(OUT) :: DSF(MDIM,NGNOD), SF(NGNOD)
      REAL*8 T, DTDS, DTDR, SUM1, SUM2, SUM3, EPS
      PARAMETER(EPS = 2.2204D-16) 
! 
!----------------------------------------------------------------------1
! 
!.... quadratic
      IF (NGNOD.EQ.6) THEN
         T = 1.D0 - R - S 
         SF(1) = R*(2.D0*R - 1.D0)
         SF(2) = S*(2.D0*S - 1.D0)
         SF(3) = T*(2.D0*T - 1.D0)
         SF(4) = 4.D0*R*S
         SF(5) = 4.D0*S*T
         SF(6) = 4.D0*R*T
         DTDR =-1.D0
         DTDS =-1.D0
         DSF(1,1) = 4.D0*R - 1.D0
         DSF(1,2) = 0.D0
         DSF(1,3) = (4.D0*T - 1.D0)*DTDR
         DSF(1,4) = 4.D0*S 
         DSF(1,5) = 4.D0*S*DTDR
         DSF(1,6) = 4.D0*(T + R*DTDR)
         DSF(2,1) = 0.D0
         DSF(2,2) = 4.D0*S - 1.D0 
         DSF(2,3) = (4.D0*T - 1.D0)*DTDS
         DSF(2,4) = 4.D0*R 
         DSF(2,5) = 4.D0*(T + S*DTDS)
         DSF(2,6) = 4.D0*R*DTDS
      ELSE !linear default
         SF(1) = R 
         SF(2) = S 
         SF(3) = 1.D0 - R - S 
         DSF(1,1) = 1.D0
         DSF(1,2) = 0.D0  
         DSF(1,3) =-1.D0
         DSF(2,1) = 0.D0
         DSF(2,2) = 1.D0
         DSF(2,3) =-1.D0
      ENDIF
! 
!.... clean up
      SUM1 = 0.D0
      SUM2 = 0.D0
      SUM3 = 0.D0
      DO 1 I=1,NGNOD
         IF (DABS(SF(I)).LT.EPS) SF(I) = 0.D0
         IF (DABS(DSF(1,I)).LT.EPS) DSF(1,I) = 0.D0
         IF (DABS(DSF(2,I)).LT.EPS) DSF(2,I) = 0.D0
         SUM1 = SUM1 + SF(I)
         SUM2 = SUM2 + DSF(1,I)
         SUM3 = SUM3 + DSF(2,I)
    1 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CSFT2DND(NGNOD, R,S, SF)
! 
!     Triangle shape functions from Hughes pg 167-168
!       
!       2         (1,3,1) x
!       |\                |\
!       | \               | \
!       5  4      (1,2,1) x  x (2,2,1)
!       |   \             |   \
!       |    \            |    \
!       3 -6- 1   (1,1,3) x -x- x (3,1,1)
!                         (2,1,2)
! 
!     INPUT      MEANING
!     -----      ------- 
!     NGNOD      number of anchor nodes (3, or 6, 3 default)
!     R          r local coordinate 
!     S          s local coordinate 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     SF         shape functions evaluated at (r,s) 
! 
!.... variable declarations
      REAL*8, INTENT(IN) :: R, S
      INTEGER*4, INTENT(IN) :: NGNOD
      REAL*8, INTENT(OUT) :: SF(NGNOD) 
      REAL*8 T, EPS 
      PARAMETER(EPS = 2.2204D-16) 
! 
!----------------------------------------------------------------------1
! 
!.... quadratic
      IF (NGNOD.EQ.6) THEN
         T = 1.D0 - R - S
         SF(1) = R*(2.D0*R - 1.D0)
         SF(2) = S*(2.D0*S - 1.D0)
         SF(3) = T*(2.D0*T - 1.D0)
         SF(4) = 4.D0*R*S
         SF(5) = 4.D0*S*T
         SF(6) = 4.D0*R*T
      ELSE !linear default
         SF(1) = R
         SF(2) = S
         SF(3) = 1.D0 - R - S
      ENDIF
! 
!.... clean up
      DO 1 I=1,NGNOD
         IF (DABS(SF(I)).LT.EPS) SF(I) = 0.D0
    1 CONTINUE
      RETURN
      END

