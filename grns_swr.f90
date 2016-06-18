
      SUBROUTINE GRNS_SWR(SRCTYP, MDIM,NDOF,NNPE, LSAVE,                    &
                          NL1D_LT,NL1D_RT, CSIDE, NDIM, ISRC,               &
                          FREQ,FREQ0,AOI,BAZN,                              &
                          XBLKL,XBLKR,ZBASE_INT,                   &
                          XMOD0,XMOD1, STF, IDOFSE,                         &
                          VP1D_LT,VS1D_LT,RH1D_LT,  &
                          QP1D_LT,QS1D_LT,Z1D_LT,   &
                          VP1D_RT,VS1D_RT,RH1D_RT,  &
                          QP1D_RT,QS1D_RT,Z1D_RT,   &
                          XLOCSE,ZLOCSE, UE,IERR)

      IMPLICIT NONE 
      CHARACTER(2), INTENT(IN) :: SRCTYP
      CHARACTER(1), INTENT(IN) :: CSIDE
      COMPLEX*16, INTENT(IN) :: STF
      REAL*8, INTENT(IN) :: XLOCSE(NNPE), ZLOCSE(NNPE),&    
              VP1D_LT(NL1D_LT), VS1D_LT(NL1D_LT), RH1D_LT(NL1D_LT),Z1D_LT(NL1D_LT),     &    
              VP1D_RT(NL1D_RT), VS1D_RT(NL1D_RT), RH1D_RT(NL1D_RT),Z1D_RT(NL1D_RT),     &    
              QP1D_RT(NL1D_RT), QP1D_LT(NL1D_LT), QS1D_RT(NL1D_RT),QS1D_LT(NL1D_LT),    &    
              FREQ,FREQ0,AOI,BAZN, XBLKL,XBLKR, XMOD0,XMOD1, ZBASE_INT
      INTEGER*4, INTENT(IN) :: IDOFSE(MDIM,*), MDIM,NDOF,NNPE,    &    
                               NL1D_LT,NL1D_RT, NDIM, ISRC 
      LOGICAL*4, INTENT(IN) :: LSAVE 
      COMPLEX*8, INTENT(OUT) :: UE(*)
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CLP(:)
      COMPLEX*16, ALLOCATABLE :: UGRN1D(:), WGRN1D(:), UGRN1D_WHL(:), WGRN1D_WHL(:)
      COMPLEX*8, ALLOCATABLE :: UGRNS(:), VGRNS(:), WGRNS(:)
      REAL*8, ALLOCATABLE :: ZPTS1D(:)
      LOGICAL*4, ALLOCATABLE :: LNPINI(:)
      COMPLEX*16 CCAZ, CSAZ, U,V,W
      REAL*8 Z1D_WHL(2), RH1D_WHL(2), VP1D_WHL(2), VS1D_WHL(2), &
             POFF, OMEGA, PX, PY, CAZ, SAZ, XOFF, YOFF, ARG,    &
             TOL, TWOPI, PI180, VPBASE,VSBASE,     &
             SIGNPX, SIGNPY
      INTEGER*4 NNP1D, INP1D,INPE
      INTEGER*4 LOCATE8
      LOGICAL*4 LFLIP, LDISP
      PARAMETER(TOL = 1.11D-7)
      PARAMETER(TWOPI = 6.2831853071795862D0)
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
      PARAMETER(LFLIP = .TRUE.)
      PARAMETER(YOFF = 0.D0) !we say that  y is in the plane
!
!----------------------------------------------------------------------------------------!
!
!.... input errors
      IERR = 0
      IF (SRCTYP(2:2) == 'H') THEN
         WRITE(*,*) 'grns_swr: Error srctyp:',SRCTYP,' not yet programmed'
         IERR = 1
         RETURN
      ENDIF
!.... initialize 
      ALLOCATE(UGRNS(NNPE))
      ALLOCATE(VGRNS(NNPE))
      ALLOCATE(WGRNS(NNPE))
      UGRNS(1:NNPE) = CMPLX(0.,0.)
      VGRNS(1:NNPE) = CMPLX(0.,0.)
      WGRNS(1:NNPE) = CMPLX(0.,0.)
      ALLOCATE(LNPINI(NNPE))
      LNPINI(1:NNPE) = .FALSE.
      ALLOCATE(ZPTS1D(NNPE+1))
      ALLOCATE(CLP(NNPE))
!
!.... set phase shift offset based on direction of propagation
      IF (CSIDE == 'L') THEN !wave moving in +x, left is zero time 
         POFF = XMOD0
         SIGNPX = 1.D0
      ELSE !wave moving in -x, right is zero time
         POFF = XMOD1
         SIGNPX =-1.D0
      ENDIF
      SIGNPY = 1.D0
      CAZ = SIGNPX*DABS(DCOS(BAZN*PI180)) !>0 moving left to right
      SAZ = DSIN(BAZN*PI180)
      IF (DABS(CAZ) < 1.11D-15) CAZ = 0.D0
      IF (DABS(SAZ) < 1.11D-15) SAZ = 0.D0
      CCAZ = DCMPLX(CAZ,0.D0)
      CSAZ = DCMPLX(SAZ,0.D0)
      OMEGA = TWOPI*FREQ
!
!.... extract baesment velocity
      IF (CSIDE == 'L') THEN
         VPBASE = VP1D_LT(NL1D_LT)
         VSBASE = VS1D_LT(NL1D_LT)
      ELSE
         VPBASE = VP1D_RT(NL1D_RT)
         VSBASE = VS1D_RT(NL1D_RT)
      ENDIF
!
!.... calculate py
      IF (SRCTYP(2:2) == 'P') THEN
         PX = DSIN(AOI*PI180)*CAZ/VPBASE
         PY = DSIN(AOI*PI180)*SAZ/VPBASE
      ELSE
         PX = DSIN(AOI*PI180)*CAZ/VSBASE
         PY = DSIN(AOI*PI180)*SAZ/VSBASE
      ENDIF

      LDISP = .TRUE. !assume dispersion
      IF (FREQ0 == 0.D0) LDISP = .FALSE.
!
!.... generate the 1D solutions for the appropriate model
      IF (CSIDE == 'L') THEN
         CLP(1:NNPE) = 'L'
         CALL SET_GRNS1D(NNPE,1,CLP,TOL,ZLOCSE, NNP1D,ZPTS1D)
         ALLOCATE(UGRN1D(NNP1D))
         ALLOCATE(WGRN1D(NNP1D))
         IF (LDISP) THEN
            CALL HASKATTN(NL1D_LT,NNP1D,SRCTYP,LFLIP, FREQ,FREQ0,AOI,      &
                          Z1D_LT,VP1D_LT,VS1D_LT,RH1D_LT,QP1D_LT,QS1D_LT,  &
                          ZPTS1D(1:NNP1D), UGRN1D,WGRN1D, IERR)
         ELSE
            CALL HASKGRN(NL1D_LT,NNP1D,SRCTYP,LFLIP, FREQ,AOI, Z1D_LT,VP1D_LT, &
                         VS1D_LT,RH1D_LT,ZPTS1D(1:NNP1D), UGRN1D,WGRN1D, IERR)
         ENDIF
         IF (IERR /= 0) THEN
            WRITE(*,*) 'grns_bdy: Error calling haskgrn 1'
            RETURN
         ENDIF
      ELSE
         CLP(1:NNPE) = 'R'
         CALL SET_GRNS1D(NNPE,2,CLP,TOL,ZLOCSE, NNP1D,ZPTS1D)
         ALLOCATE(UGRN1D(NNP1D))
         ALLOCATE(WGRN1D(NNP1D))
         IF (LDISP) THEN
            CALL HASKATTN(NL1D_RT,NNP1D,SRCTYP,LFLIP, FREQ,FREQ0,AOI,      &
                          Z1D_RT,VP1D_RT,VS1D_RT,RH1D_RT,QP1D_RT,QS1D_RT,  &
                          ZPTS1D(1:NNP1D), UGRN1D,WGRN1D, IERR)
         ELSE
            CALL HASKGRN(NL1D_RT,NNP1D,SRCTYP,LFLIP, FREQ,AOI, Z1D_RT,VP1D_RT, &
                         VS1D_RT,RH1D_RT,ZPTS1D(1:NNP1D), UGRN1D,WGRN1D, IERR)
         ENDIF
         IF (IERR /= 0) THEN
            WRITE(*,*) 'grns_bdy: Error calling haskgrn 1'
            RETURN
         ENDIF
      ENDIF
      DEALLOCATE(CLP)
!
!.... now generate the solution in the whole space
      ALLOCATE(UGRN1D_WHL(NNP1D))
      ALLOCATE(WGRN1D_WHL(NNP1D))
      UGRN1D_WHL(:) = DCMPLX(0.D0,0.D0)
      WGRN1D_WHL(:) = DCMPLX(0.D0,0.D0) 
      IF (CSIDE == 'L') THEN
         Z1D_WHL(1) = Z1D_LT(1) 
         Z1D_WHL(2) = Z1D_LT(NL1D_LT)
         VP1D_WHL(1:2) = VP1D_LT(NL1D_LT)
         VS1D_WHL(1:2) = VS1D_LT(NL1D_LT)
         RH1D_WHL(1:2) = RH1D_LT(NL1D_LT)
      ELSE
         Z1D_WHL(1) = Z1D_RT(1) 
         Z1D_WHL(2) = Z1D_RT(NL1D_RT)
         VP1D_WHL(1:2) = VP1D_RT(NL1D_RT)
         VS1D_WHL(1:2) = VS1D_RT(NL1D_RT)
         RH1D_WHL(1:2) = RH1D_RT(NL1D_RT)
      ENDIF
      CALL HASKINF(2,NNP1D,SRCTYP,LFLIP, FREQ,AOI, Z1D_WHL,VP1D_WHL,  &
                   VS1D_WHL,RH1D_WHL,ZPTS1D(1:NNP1D), UGRN1D_WHL,WGRN1D_WHL, IERR) 
      IF (IERR /= 0) THEN
         WRITE(*,*) 'grns_swr: Error calling haskinf!'
         RETURN
      ENDIF
!
!.... set points
      DO 100 INPE=1,NNPE
         INP1D = LOCATE8(NNP1D,TOL,ZLOCSE(INPE),ZPTS1D)
         IF (INP1D < 1) THEN
            WRITE(*,*) 'grns_swr: Could not locate point!'
            IERR = 1
            GOTO 55
         ENDIF
         LNPINI(INPE) = .TRUE.
!
!....... calculate offset and phase shifter
         XOFF = XLOCSE(INPE) - POFF
         ARG =-OMEGA*(XOFF*PX + YOFF*PY) !-omega*p.x
!
!....... may want to use points in middle
         IF (XLOCSE(INPE) + 1.D-6 > XBLKL .AND. XLOCSE(INPE) - 1.D-6 < XBLKR) THEN
            U = UGRN1D_WHL(INP1D)*CCAZ !u component, pure radial
            V = UGRN1D_WHL(INP1D)*CSAZ !v component, pure transverse
            W = WGRN1D_WHL(INP1D)      !w component, pure vertical
         ELSE !Otherwise use 1D layered model with free surface on sides
            U = UGRN1D(INP1D)*CCAZ !u component, pure radial
            V = UGRN1D(INP1D)*CSAZ !v component, pure transverse
            W = WGRN1D(INP1D)      !w component, pure vertical
         ENDIF
!
!....... convolve STF and shift
         U = U*CDEXP(DCMPLX(0.D0,ARG))*STF!phase shift/convolve
         V = V*CDEXP(DCMPLX(0.D0,ARG))*STF!phase shift/convolve
         W = W*CDEXP(DCMPLX(0.D0,ARG))*STF!phase shift/convolve 
         IF (SRCTYP(2:2) == 'P' .OR. SRCTYP(2:2) == 'S') THEN !incoming P/SV
            UGRNS(INPE) = CMPLX(U)
            VGRNS(INPE) = CMPLX(V)
         ELSE !incoming SH polarization
            UGRNS(INPE) =-CMPLX(V)
            VGRNS(INPE) = CMPLX(U)
         ENDIF
         WGRNS(INPE) =-CMPLX(W)!point up now 
  100 CONTINUE
   55 CONTINUE
!
!.... fill in the DOFs  or save ?
      IF (.NOT.LSAVE) THEN 
         CALL FILL_GRNS(MDIM,NDOF,NNPE, NDIM,LNPINI,IDOFSE, UGRNS,VGRNS,WGRNS, UE,IERR)
         IF (IERR /= 0) THEN 
            WRITE(*,*) 'grns_swr: Serious warning in fill_grns'
            IERR = 0
         ENDIF
      ELSE 
         DO 45 INPE=1,NNPE
            IF (.NOT.LNPINI(INPE)) & 
            WRITE(*,*) 'grns_swr: Serious warning in grns_bdy, unitialized point',INPE
   45    CONTINUE 
         CALL SAVE_GRNS(NDIM, NNPE,NDIM,1,ISRC,.FALSE., FREQ,IDOFSE,UGRNS,VGRNS,WGRNS) 
      ENDIF
      DEALLOCATE(UGRNS)
      DEALLOCATE(WGRNS) 
      DEALLOCATE(UGRN1D)
      DEALLOCATE(WGRN1D)
      DEALLOCATE(ZPTS1D)
      DEALLOCATE(UGRN1D_WHL)
      DEALLOCATE(WGRN1D_WHL)
      RETURN
      END  
