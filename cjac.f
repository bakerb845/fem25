      SUBROUTINE CJMK_TVI(MEN,MINTX,MINTZ,MNPG, LISISO,
     ;                    NINTX,NINTZ,NEN,NGNOD,NNPG, 
     ;                    IENG,XIGLL,ETAGLL,
     ;                    XLOCS,ZLOCS,DENS,ECOEFF,SHL, 
     ;                    DET,RHO,DMAT,SHG,IERR) 
!
!     Calculates the Jacobian and interpolated densities for 
!     integration of the stiffness and mass matrices for transverse
!     vertical isotropy 
!
!     INPUT      INPUT  
!     -----      ----- 
!     DENS       densities at the global nodes
!     ECOEFF     elastic coefficients at anchor nodes, see DMAT
!     ETAGLL     holds integration points/weights in eta
!     IENG       global anchor node pointer
!     LISISO     True -> simulation isotrpic, False -> anisotropic
!     MEN        max element nodes
!     MINTX      max integration points in x 
!     MINTZ      max integration points in z 
!     MNPG       max number of global anchor nodes
!     NEN        number of element nodes 
!     NGNOD      number of global element nodes (anchor nodes, 4 or 9)
!     NINTX      number of integration points in x 
!     NINTZ      number of integration points in z 
!     NNPG       number of nodal points in global mesh
!     SHL        lagrange shape function array 
!     XIGLL      holds integration points/weights in xi 
!     XLOCS      physical x locations of anchor nodes 
!     ZLOCS      physical z locations of anchor nodes
!  
!     OUTPUT     MEANING 
!     ------     ------- 
!     DET        the jacobian (determinant) at integration points 
!     DMAT       elastic coefficient matrix at points  
!                if isiso  DMAT(intx,intz,1) = lambda 
!                          DMAT(intx,intz,2) = mu
!                otherwise DMAT(intx,intz,1) = C_11 = C_22 = c_{1111} 
!                          DMAT(intx,intz,2) = C_13 = C_23 = c_{1133} 
!                          DMAT(intx,intz,3) = C_33        = c_{3333} 
!                          DMAT(intx,intz,4) = C_44 = C_55 = c_{2323}
!                          DMAT(intx,intz,5) = C_66 = c_{6666}
!                Recall that C_12 = C_11 - 2 C_66
!     IERR       = 1 then jacobian is undefined or negative 
!     RHO        densities at integration points
!     SHG        shape functions at integration points 
! 
!.... variable declarations 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8, INTENT(IN) :: SHL(3,MEN,MINTX,*), ECOEFF(MNPG,*), 
     ;                      XLOCS(NNPG), ZLOCS(NNPG), DENS(NNPG), 
     ;                      XIGLL(NINTX), ETAGLL(NINTZ) 
      INTEGER*4, INTENT(IN) :: IENG(NGNOD), MEN,MINTX,MINTZ,MNPG, 
     ;                         NINTX,NINTZ,NEN,NGNOD,NNPG
      LOGICAL*4, INTENT(IN) :: LISISO
      REAL*8, INTENT(OUT) :: SHG(3,MEN,MINTX,*), DMAT(MINTX,MINTZ,*), 
     ;                       DET(MINTX,*), RHO(MINTX,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables 
      REAL*8 DSF(2,9), SF(9), X, Z, DXDXI, DXDETA, DZDXI, DZDETA, RJAC,
     ;       XI, ETA, EPS
      PARAMETER(EPS = 2.22D-16)
! 
!----------------------------------------------------------------------!
!
!.... check number of coefficients
      IERR = 0
      IF (LISISO) THEN
         NCOEFF = 2
      ELSE
         NCOEFF = 5
      ENDIF
!
!.... loop on integration points in x 
      DO 1 INTX=1,NINTX
         XI = XIGLL(INTX) 
         DO 2 INTZ=1,NINTZ
! 
!.......... create the 2D shape functions at integration points 
            ETA = ETAGLL(INTZ)
            CALL CSF2D(2,NGNOD,XI,ETA, SF,DSF,IERR)
            IF (IERR.NE.0) THEN
               WRITE(*,*) 'cjmk_tvi: Error calling csf2d'
               RETURN
            ENDIF
! 
!.......... find the 2D shape fns and their derivatives 
            X = 0.D0
            Z = 0.D0
            RHO(INTX,INTZ) = 0.D0
            DO 3 INDX=1,NCOEFF 
               DMAT(INTX,INTZ,INDX) = 0.D0
    3       CONTINUE
            DXDXI = 0.D0
            DXDETA = 0.D0
            DZDXI = 0.D0
            DZDETA = 0.D0
            DO 4 IA=1,NGNOD !loop over anchor nodes  
! 
!............. interpolate x, z, rho, and elastic coefficients 
               ILOC = IENG(IA)
               X = X + SF(IA)*XLOCS(ILOC)
               Z = Z + SF(IA)*ZLOCS(ILOC)
               DO 5 INDX=1,NCOEFF 
                  DMAT(INTX,INTZ,INDX) = DMAT(INTX,INTZ,INDX)
     ;                                 + SF(IA)*ECOEFF(ILOC,INDX)
    5          CONTINUE
               RHO(INTX,INTZ) = RHO(INTX,INTZ) + SF(IA)*DENS(ILOC)
               DXDXI  = DXDXI  + DSF(1,IA)*XLOCS(ILOC)
               DXDETA = DXDETA + DSF(2,IA)*XLOCS(ILOC)
               DZDXI  = DZDXI  + DSF(1,IA)*ZLOCS(ILOC)
               DZDETA = DZDETA + DSF(2,IA)*ZLOCS(ILOC)
    4       CONTINUE
! 
!.......... check on jacobian 
            RJAC = DXDXI*DZDETA - DXDETA*DZDXI
            DET(INTX,INTZ) = RJAC
            IF (DABS(RJAC).LE.EPS .OR. RJAC.LT.0.D0) THEN
               WRITE(*,*) 'cjkm_tvi: Jacobian undefined!'
               WRITE(*,*) 'cjkm_tvi: Jacobian is ', RJAC
               IF (RJAC.LT.0.D0)
     ;         WRITE(*,*)'cjkm_tvi: Jacobian can not be negative'
               IERR = 1
               RETURN
            ENDIF
! 
!.......... loop over element nodes and fill out global shape fns 
            DO 6 IAE=1,NEN
               DRNXI  = SHL(1,IAE,INTX,INTZ) !dN_a/d_xi
               DRNETA = SHL(2,IAE,INTX,INTZ) !dN_a/d_eta
               SHG(3,IAE,INTX,INTZ) = SHL(3,IAE,INTX,INTZ)
               SHG(1,IAE,INTX,INTZ) = (DRNXI*DZDETA - DRNETA*DZDXI)/RJAC
               SHG(2,IAE,INTX,INTZ) =-(DRNXI*DXDETA - DRNETA*DXDXI)/RJAC
    6       CONTINUE
    2    CONTINUE !end loop on integration points in z 
    1 CONTINUE !end loop on integration points in x
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CJKM_QUAL(MEN,MINTX,MNPG, LISISO,
     ;                     NINTX,NINTZ,NEN,NGNOD,NNPG,
     ;                     IENG,XIGLL,ETAGLL,
     ;                     XLOCS,ZLOCS,DENS,QP,QS,ECOEFF,SHL,
     ;                     DET,VP,VS,QPP,QSS,SHG,IERR) 
!
!     Calculates the Jacobian, Quality factors, and velocities at 
!     at integration pionts
! 
!     INPUT      INPUT  
!     -----      ----- 
!     DENS       densities at the global nodes (kg/m**3)
!     ECOEFF     elastic coefficients at anchor nodes, (lambda,mu); Pa
!     ETAGLL     holds integration points/weights in eta
!     IENG       global anchor node pointer
!     LISISO     True -> simulation isotrpic, False -> error
!     MEN        max element nodes
!     MINTX      max integration points in x 
!     MINTZ      max integration points in z 
!     MNPG       max number of global anchor nodes
!     NEN        number of element nodes 
!     NGNOD      number of global element nodes (anchor nodes, 4 or 9)
!     NINTX      number of integration points in x 
!     NINTZ      number of integration points in z 
!     NNPG       number of nodal points in global mesh
!     QP         Vp quality factor at anchor nodes
!     QS         Vs quality factor at anchor nodes 
!     SHL        lagrange shape function array 
!     XIGLL      holds integration points/weights in xi 
!     XLOCS      physical x locations of anchor nodes 
!     ZLOCS      physical z locations of anchor nodes
!
!  
!     OUTPUT     MEANING 
!     ------     ------- 
!     DET        the jacobian (determinant) at integration points 
!     IERR       = 1 then jacobian is undefined or negative
!     QPP        Vp quality factor at integration points
!     QSS        Vs quality factor at integration points 
!     SHG        shape functions at integration points 
!     VP         compressional velocities at integration points (m/s)
!     VS         shear velocities at integration points (m/s)
!
!.... variable declarations 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8, INTENT(IN) :: SHL(3,MEN,MINTX,*), ECOEFF(MNPG,*), 
     ;                      XLOCS(NNPG), ZLOCS(NNPG), DENS(NNPG),
     ;                      QP(NNPG), QS(NNPG), 
     ;                      XIGLL(NINTX), ETAGLL(NINTZ)
      INTEGER*4, INTENT(IN) :: IENG(NGNOD), MEN,MINTX,MNPG,
     ;                         NINTX,NINTZ,NEN,NGNOD,NNPG
      LOGICAL*4, INTENT(IN) :: LISISO
      REAL*8, INTENT(OUT) :: SHG(3,MEN,MINTX,*), VP(MINTX,*), 
     ;                       VS(MINTX,*), QPP(MINTX,*), QSS(MINTX,*), 
     ;                       DET(MINTX,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables 
      REAL*8 DSF(2,9), SF(9), X, Z, DXDXI, DXDETA, DZDXI, DZDETA, RJAC,
     ;       XI, ETA, RLAM, RMU, VPP, VSS, EPS
      PARAMETER(EPS = 2.22D-16)
!
!----------------------------------------------------------------------!
!
      IERR = 0 
      IF (.NOT.LISISO) THEN
         WRITE(*,*) 'cjkm_qual: Error anisotropy not an option!' 
         IERR = 1
      ENDIF
!
!.... loop on integration points in x 
      DO 1 INTX=1,NINTX
         XI = XIGLL(INTX) 
         DO 2 INTZ=1,NINTZ
! 
!.......... create the 2D shape functions at integration points 
            ETA = ETAGLL(INTZ)
            CALL CSF2D(2,NGNOD,XI,ETA, SF,DSF,IERR)
            IF (IERR.NE.0) THEN
               WRITE(*,*) 'cjmk_qual: Error calling csf2d'
               RETURN
            ENDIF
! 
!.......... find the 2D shape fns and their derivatives 
            X = 0.D0
            Z = 0.D0
            VP(INTX,INTZ) = 0.D0
            VS(INTX,INTZ) = 0.D0
            QPP(INTX,INTZ) = 0.D0
            QSS(INTX,INTZ) = 0.D0
            DXDXI = 0.D0
            DXDETA = 0.D0
            DZDXI = 0.D0
            DZDETA = 0.D0
            DO 4 IA=1,NGNOD !loop over anchor nodes   
! 
!............. interpolate x, z, rho, and elastic coefficients 
               ILOC = IENG(IA) 
               X = X + SF(IA)*XLOCS(ILOC)
               Z = Z + SF(IA)*ZLOCS(ILOC)
               RLAM = ECOEFF(ILOC,1)
               RMU  = ECOEFF(ILOC,2)
               VPP = DSQRT( (RLAM + 2.D0*RMU)/DENS(ILOC) )
               VSS = DSQRT(RMU/DENS(ILOC))
               VP(INTX,INTZ)  = VP(INTX,INTZ) + SF(IA)*VPP 
               VS(INTX,INTZ)  = VS(INTX,INTZ) + SF(IA)*VSS
               QPP(INTX,INTZ) = QPP(INTX,INTZ) + SF(IA)*QP(ILOC)
               QSS(INTX,INTZ) = QSS(INTX,INTZ) + SF(IA)*QS(ILOC)
               DXDXI  = DXDXI  + DSF(1,IA)*XLOCS(ILOC)
               DXDETA = DXDETA + DSF(2,IA)*XLOCS(ILOC)
               DZDXI  = DZDXI  + DSF(1,IA)*ZLOCS(ILOC)
               DZDETA = DZDETA + DSF(2,IA)*ZLOCS(ILOC)
    4       CONTINUE
! 
!.......... check on jacobian 
            RJAC = DXDXI*DZDETA - DXDETA*DZDXI
            DET(INTX,INTZ) = RJAC
            IF (DABS(RJAC).LE.EPS .OR. RJAC.LT.0.D0) THEN
               WRITE(*,*) 'cjkm_qual: Jacobian undefined!'
               WRITE(*,*) 'cjkm_qual: Jacobian is ', RJAC
               IF (RJAC.LT.0.D0)
     ;         WRITE(*,*)'cjkm_qual: Jacobian can not be negative'
               IERR = 1
               RETURN
            ENDIF
! 
!.......... loop over element nodes and fill out global shape fns 
            DO 6 IAE=1,NEN
               DRNXI  = SHL(1,IAE,INTX,INTZ) !dN_a/d_xi
               DRNETA = SHL(2,IAE,INTX,INTZ) !dN_a/d_eta
               SHG(3,IAE,INTX,INTZ) = SHL(3,IAE,INTX,INTZ)
               SHG(1,IAE,INTX,INTZ) = (DRNXI*DZDETA - DRNETA*DZDXI)/RJAC
               SHG(2,IAE,INTX,INTZ) =-(DRNXI*DXDETA - DRNETA*DXDXI)/RJAC
    6       CONTINUE
    2    CONTINUE !end loop on integration points in z 
    1 CONTINUE !end loop on integration points in x
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CJMK_VPVS(MEN,MINTX,
     ;                     NINTX,NINTZ,NEN,NGNOD,NNPG,
     ;                     IENG,XIGLL,ETAGLL,
     ;                     XLOCS,ZLOCS,VPANC,VSANC,SHL,
     ;                     DET,VP,VS,SHG,IERR)
!
!     This is for debugging the inverse problem.  Calculates the vp 
!     and vs at integration nodes along with the Jacobian 
!     - B. Baker April 2013 
!
!     INPUT     INPUG 
!     -----     -----
!     ETAGLL    eta integration abscissas
!     IENG      IEN pointer for anchor nodes
!     MEN       leading dimension
!     MINTX     leading dimension
!     NINTX     number of integration points in xi
!     NINTZ     number of integration points in eta
!     NEN       number of element nodes
!     NGNOD     number of anchor nodes on an element
!     NNPG      number of anchor nodes in mesh
!     SHL       local shape fns and derivatives at int. pts 
!     VPANC     p velocity at anchor nodes
!     VSANC     s velocity at anchor nodes
!     XIGLL     xi integration absicssas
!     XLOCS     x anchor node locations
!     ZLOCS     z anchor node locations
!
!     OUTPUT    MEANING
!     ------    -------
!     DET       element determinant
!     IERR      error flag
!     SHG       global shape fns and derivatives at int. pts
!     VP        p velocity at integration points
!     VS        s velocity at integration points  
!
!.... variable declarations 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8, INTENT(IN) :: SHL(3,MEN,MINTX,*),VPANC(NNPG),VSANC(NNPG),
     ;                      XLOCS(NNPG), ZLOCS(NNPG), 
     ;                      XIGLL(NINTX), ETAGLL(NINTZ)
      INTEGER*4, INTENT(IN) :: IENG(NGNOD), MEN,MINTX,
     ;                         NINTX,NINTZ,NEN,NGNOD,NNPG
      REAL*8, INTENT(OUT) :: SHG(3,MEN,MINTX,*), VP(MINTX,*),
     ;                       VS(MINTX,*), DET(MINTX,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables 
      REAL*8 DSF(2,9), SF(9), X, Z, DXDXI, DXDETA, DZDXI, DZDETA, RJAC,
     ;       XI, ETA, EPS
      PARAMETER(EPS = 2.22D-16)
!
!----------------------------------------------------------------------!
!
!.... loop on integration points in x 
      IERR = 0
      DO 1 INTX=1,NINTX
         XI = XIGLL(INTX) 
         DO 2 INTZ=1,NINTZ
! 
!.......... create the 2D shape functions at integration points 
            ETA = ETAGLL(INTZ)
            CALL CSF2D(2,NGNOD,XI,ETA, SF,DSF,IERR)
            IF (IERR.NE.0) THEN
               WRITE(*,*) 'cjmk_vpvs: Error calling csf2d'
               RETURN
            ENDIF
! 
!.......... find the 2D shape fns and their derivatives 
            X = 0.D0
            Z = 0.D0
            VP(INTX,INTZ) = 0.D0
            VS(INTX,INTZ) = 0.D0
            DXDXI = 0.D0
            DXDETA = 0.D0
            DZDXI = 0.D0
            DZDETA = 0.D0
            DO 4 IA=1,NGNOD !loop over anchor nodes   
! 
!............. interpolate x, z, rho, and elastic coefficients 
               ILOC = IENG(IA)
               X = X + SF(IA)*XLOCS(ILOC)
               Z = Z + SF(IA)*ZLOCS(ILOC)
               VP(INTX,INTZ)  = VP(INTX,INTZ) + SF(IA)*VPANC(ILOC)
               VS(INTX,INTZ)  = VS(INTX,INTZ) + SF(IA)*VSANC(ILOC)
               DXDXI  = DXDXI  + DSF(1,IA)*XLOCS(ILOC)
               DXDETA = DXDETA + DSF(2,IA)*XLOCS(ILOC)
               DZDXI  = DZDXI  + DSF(1,IA)*ZLOCS(ILOC)
               DZDETA = DZDETA + DSF(2,IA)*ZLOCS(ILOC)
    4       CONTINUE
! 
!.......... check on jacobian 
            RJAC = DXDXI*DZDETA - DXDETA*DZDXI
            DET(INTX,INTZ) = RJAC
            IF (DABS(RJAC).LE.EPS .OR. RJAC.LT.0.D0) THEN
               WRITE(*,*) 'cjkm_vpvs: Jacobian undefined!'
               WRITE(*,*) 'cjkm_vpvs: Jacobian is ', RJAC
               IF (RJAC.LT.0.D0)
     ;         WRITE(*,*)'cjkm_vpvs: Jacobian can not be negative'
               IERR = 1
               RETURN
            ENDIF
! 
!.......... loop over element nodes and fill out global shape fns 
            DO 6 IAE=1,NEN
               DRNXI  = SHL(1,IAE,INTX,INTZ) !dN_a/d_xi
               DRNETA = SHL(2,IAE,INTX,INTZ) !dN_a/d_eta
               SHG(3,IAE,INTX,INTZ) = SHL(3,IAE,INTX,INTZ)
               SHG(1,IAE,INTX,INTZ) = (DRNXI*DZDETA - DRNETA*DZDXI)/RJAC
               SHG(2,IAE,INTX,INTZ) =-(DRNXI*DXDETA - DRNETA*DXDXI)/RJAC
    6       CONTINUE
    2    CONTINUE !end loop on integration points in z 
    1 CONTINUE !end loop on integration points in x
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CJACCQ(MEN,MINTX,MINTZ,MNPG, LISISO, 
     ;                  LCONV,NINTX,NINTZ,NEN,NGNOD,NNPG, 
     ;                  OMEGA,XWIDTH,ZWIDTH, 
     ;                  RCOEFF,DMAX_IN,RKAPPA,FCENT, 
     ;                  IENG, XIGLL,ETAGLL, 
     ;                  XD,ZD,XLOCS,ZLOCS, DENS,ECOEFF,
     ;                  SHL, DET,RHO,DMAT,SHG, DAMPX,DAMPZ, IERR)
! 
!     Calculates the global shape functions, derivatives, and jacobian
!     at each integration point.  This is for use with the element 
!     damping and element mass 
! 
!     INPUT      MEANING 
!     -----      ------- 
!     DENS       density at each grid point 
!     DSF        workspace array; holds derivative of shape fns
!     ECOEFF     elastic coefficients 
!     ETAGLL     holds integration points/weights in eta
!     IENG       points to global anchor node numbers
!     LCONV      True -> use convolution C-PML, False -> old pml
!     LISISO     True -> isotropic simulation; False -> anisotropic
!     MEN        max number of element nodes 
!     MINTX      max number of integration points in x 
!     MINTZ      max number of integration points in z 
!     NGNOD      number of anchor nodes 
!     NINTX      number of integration points in x 
!     NINTZ      number of integration points in z  
!     SHL        holds shape functions evaluated at integration points
!     SF         workspace array; holds shape functions 
!     XIGLL      holds integration points/weights in xi 
!     XD         x distance into PML 
!     XLOCS      x anchor nodes to build Jacobian on 
!     XWIDTH     x width of PML
!     ZD         z distance into PML 
!     ZLOCS      z locations of anchor nodes
!     ZWIDTH     z width of PML
! 
!     OUTPUT     MEANING 
!     ------     ------- 
!     DAMPX      damping function in x at integration points
!     DAMPZ      damping function in z at integration points
!     DET        jacobian at integration points 
!     DMAT       elastic coefficient matrix at integration points 
!     DMAT       elastic coefficient matrix at points  
!                if isiso  DMAT(L,1) = lambda 
!                          DMAT(L,2) = mu
!                otherwise DMAT(intx,intz,1) = C_11 = C_22 = c_{1111} 
!                          DMAT(intx,intz,2) = C_13 = C_23 = c_{1133} 
!                          DMAT(intx,intz,3) = C_33        = c_{3333} 
!                          DMAT(intx,intz,4) = C_44 = C_55 = c_{2323}
!                          DMAT(intx,intz,5) = C_66 = c_{6666}
!     DAMPX      damping parameter in x direction 
!     DAMPZ      damping parameter in z direction 
!     IERR       = 1 then jacobian is undefined 
!     SHG        global shape functions 
! 
!.... variable declarations 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8, INTENT(IN) :: SHL(3,MEN,MINTX,*), ECOEFF(MNPG,*),
     ;                      XIGLL(NINTX), ETAGLL(NINTZ), 
     ;                      XLOCS(NNPG), ZLOCS(NNPG), XD(NNPG), 
     ;                      ZD(NNPG), DENS(NNPG), XWIDTH,ZWIDTH,OMEGA, 
     ;                      RCOEFF, RKAPPA, DMAX_IN, FCENT
      INTEGER*4, INTENT(IN) :: IENG(NGNOD),MEN,MINTX,MINTZ,NGNOD,NNPG,
     ;                         MNPG,NINTX,NINTZ,NEN
      LOGICAL*4, INTENT(IN) :: LCONV,LISISO
      COMPLEX*16, INTENT(OUT) :: DAMPX(MINTX,*), DAMPZ(MINTX,*)
      REAL*8, INTENT(OUT) :: DMAT(MINTX,MINTZ,*), SHG(3,MEN,MINTX,*),
     ;                       DET(MINTX,*), RHO(MINTX,*) 
      INTEGER*4, INTENT(OUT) :: IERR
      COMPLEX*16 CONE
      REAL*8 DSF(2,9), SF(9)
      REAL*8 ALPHA2,EPS,X,Z,XDEP,ZDEP,RJAC,XI,ETA,VEL,
     ;       DXDXI,DXDETA,DZDXI,DZDETA,DRNXI,DRNETA
      INTEGER*4 NCOEFF,INTX,INTZ,IAE,IA,I,ILOC
      PARAMETER(EPS = 2.22D-16)
      PARAMETER(CONE = DCMPLX(1.D0,0.D0))
      COMPLEX*16 ZDAMPFN
! 
!----------------------------------------------------------------------!
!
!.... possible warning 
      IERR = 0
      IF (.NOT.LISISO) THEN
         IF (.NOT.LCONV) THEN
            WRITE(*,*) 'cjaccq: Warning using quasi-P velocity'
         ENDIF
         NCOEFF = 5
      ELSE
         NCOEFF = 2
      ENDIF
! 
!.... loop on integration points
      DO 1 INTX=1,NINTX
         XI = XIGLL(INTX)
         DO 2 INTZ=1,NINTZ
! 
!.......... create the 2D shape functions at integration points 
            ETA = ETAGLL(INTZ)
            CALL CSF2D(2,NGNOD,XI,ETA, SF,DSF,IERR)
            IF (IERR.NE.0) THEN
               WRITE(*,*) 'cjaccq: Error calling csf2d'
               RETURN
            ENDIF
! 
!.......... find the 2D shape fns and their derivatives 
            X = 0.D0
            Z = 0.D0
            XDEP = 0.D0
            ZDEP = 0.D0
            VEL = 0.D0
            RHO(INTX,INTZ) = 0.D0
            DO 3 I=1,NCOEFF
               DMAT(INTX,INTZ,I) = 0.D0
    3       CONTINUE
            DXDXI = 0.D0
            DXDETA = 0.D0
            DZDXI = 0.D0
            DZDETA = 0.D0
            DO 4 IA=1,NGNOD !loop over anchor nodes  
! 
!............. use shape fns to interpolate x, z, and elastic coeffs
               ILOC = IENG(IA)
               X = X + SF(IA)*XLOCS(ILOC)
               Z = Z + SF(IA)*ZLOCS(ILOC)
               DO 5 I=1,NCOEFF
                  DMAT(INTX,INTZ,I) = DMAT(INTX,INTZ,I)
     ;                              + SF(IA)*ECOEFF(ILOC,I)
    5          CONTINUE
               RHO(INTX,INTZ) = RHO(INTX,INTZ) + SF(IA)*DENS(ILOC)
! 
!............. interpolate PML damping depth
               XDEP = XDEP + SF(IA)*XD(ILOC)
               ZDEP = ZDEP + SF(IA)*ZD(ILOC)
               ALPHA2 = (ECOEFF(ILOC,1) + 2.D0*ECOEFF(ILOC,2))
     ;                 /DENS(ILOC)
               VEL = VEL + SF(IA)*DSQRT(ALPHA2)
! 
!............. and derivatives
               DXDXI  = DXDXI  + DSF(1,IA)*XLOCS(ILOC)
               DXDETA = DXDETA + DSF(2,IA)*XLOCS(ILOC)
               DZDXI  = DZDXI  + DSF(1,IA)*ZLOCS(ILOC)
               DZDETA = DZDETA + DSF(2,IA)*ZLOCS(ILOC)
    4       CONTINUE !loop on anchor nodes
! 
!.......... evaluate damping profile at location 
            IF (XDEP.GT.EPS) THEN
               DAMPX(INTX,INTZ) = ZDAMPFN(XDEP,VEL,XWIDTH, DMAX_IN, 
     ;                                    RCOEFF,RKAPPA,FCENT,OMEGA)
            ELSE
               DAMPX(INTX,INTZ) = CONE 
            ENDIF
            IF (ZDEP.GT.EPS) THEN
               DAMPZ(INTX,INTZ) = ZDAMPFN(ZDEP,VEL,ZWIDTH, DMAX_IN,
     ;                                    RCOEFF,RKAPPA,FCENT,OMEGA)
            ELSE
               DAMPZ(INTX,INTZ) = CONE 
            ENDIF
! 
!.......... check on jacobian 
            RJAC = DXDXI*DZDETA - DXDETA*DZDXI
            DET(INTX,INTZ) = RJAC
            IF (DABS(RJAC).LE.EPS .OR. RJAC.LT.0.D0) THEN
               WRITE(*,*) 'cjaccq: Jacobian undefined!'
               WRITE(*,*) 'cjaccq: Jacobian is ', RJAC
               IF (RJAC.LT.0.D0)
     ;         WRITE(*,*)'cjaccq: Jacobian can not be negative'
               IERR = 1
               RETURN
            ENDIF
! 
!.......... loop over element nodes and fill out global shape fns 
            DO 6 IAE=1,NEN
               DRNXI  = SHL(1,IAE,INTX,INTZ) !dN_a/d_xi
               DRNETA = SHL(2,IAE,INTX,INTZ) !dN_a/d_eta
               SHG(3,IAE,INTX,INTZ) = SHL(3,IAE,INTX,INTZ)
               SHG(1,IAE,INTX,INTZ) = (DRNXI*DZDETA - DRNETA*DZDXI)/RJAC
               SHG(2,IAE,INTX,INTZ) =-(DRNXI*DXDETA - DRNETA*DXDXI)/RJAC
    6       CONTINUE
    2    CONTINUE !end loop on integration points in z 
    1 CONTINUE !end loop on integration points in x
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CJCM_QUAL(MEN,MINTX,MNPG, LISISO,
     ;                     NINTX,NINTZ,NEN,NGNOD,NNPG,
     ;                     IENG,XIGLL,ETAGLL,
     ;                     OMEGA,XWIDTH,ZWIDTH, 
     ;                     RCOEFF,DMAX_IN,RKAPPA,FCENT, XD,ZD,  
     ;                     XLOCS,ZLOCS,DENS,QP,QS,ECOEFF,SHL,
     ;                     DET,VP,VS,QPP,QSS,SHG, DAMPX,DAMPZ, IERR)
!
!     Calculates the Jacobian, Quality factors, velocities, and
!     PML damping at integration points 
! 
!     INPUT      INPUT  
!     -----      ----- 
!     DENS       densities at the global nodes (kg/m**3)
!     ECOEFF     elastic coefficients at anchor nodes, (lambda,mu); Pa
!     ETAGLL     holds integration points/weights in eta
!     IENG       global anchor node pointer
!     LISISO     True -> simulation isotrpic, False -> error
!     MEN        max element nodes
!     MINTX      max integration points in x 
!     MINTZ      max integration points in z 
!     MNPG       max number of global anchor nodes
!     NEN        number of element nodes 
!     NGNOD      number of global element nodes (anchor nodes, 4 or 9)
!     NINTX      number of integration points in x 
!     NINTZ      number of integration points in z 
!     NNPG       number of nodal points in global mesh
!     OMEGA      angular frequency
!     QP         Vp quality factor at anchor nodes
!     QS         Vs quality factor at anchor nodes 
!     SHL        lagrange shape function array 
!     XIGLL      holds integration points/weights in xi 
!     XD         x distance into PML
!     XLOCS      physical x locations of anchor nodes 
!     XWIDTH     x width of PML
!     ZD         z distance into PML
!     ZLOCS      physical z locations of anchor nodes
!     ZWIDTH     z width of PML 
!
!  
!     OUTPUT     MEANING 
!     ------     ------- 
!     DAMPX      PML damping profile in x at integration points
!     DAMPZ      PML damping profile in z at integration points 
!     DET        the jacobian (determinant) at integration points 
!     IERR       = 1 then jacobian is undefined or negative
!     QPP        Vp quality factor at integration points
!     QSS        Vs quality factor at integration points 
!     SHG        shape functions at integration points 
!     VP         compressional velocities at integration points (m/s)
!     VS         shear velocities at integration points (m/s)
!
!.... variable declarations 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8, INTENT(IN) :: SHL(3,MEN,MINTX,*), ECOEFF(MNPG,*),
     ;                      XLOCS(NNPG), ZLOCS(NNPG), DENS(NNPG),
     ;                      QP(NNPG), QS(NNPG), XD(NNPG), ZD(NNPG), 
     ;                      XIGLL(NINTX), ETAGLL(NINTZ), OMEGA, 
     ;                      XWIDTH, ZWIDTH, RCOEFF,DMAX_IN,RKAPPA,FCENT
      INTEGER*4, INTENT(IN) :: IENG(NGNOD), MEN,MINTX,MNPG,
     ;                         NINTX,NINTZ,NEN,NGNOD,NNPG
      LOGICAL*4, INTENT(IN) :: LISISO
      COMPLEX*16, INTENT(OUT) :: DAMPX(MINTX,*), DAMPZ(MINTX,*)
      REAL*8, INTENT(OUT) :: SHG(3,MEN,MINTX,*), VP(MINTX,*),
     ;                       VS(MINTX,*), QPP(MINTX,*), QSS(MINTX,*),
     ;                       DET(MINTX,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables 
      COMPLEX*16 CONE, ZDAMPFN 
      REAL*8 DSF(2,9), SF(9), X, Z, DXDXI, DXDETA, DZDXI, DZDETA, RJAC,
     ;       XI, ETA, RLAM, RMU, VPP, VSS, VEL, EPS,  
     ;       XDEP, ZDEP, DRNXI, DRNETA 
      INTEGER*4 INTX, INTZ, IAE, IA, ILOC 
      PARAMETER(CONE = DCMPLX(1.D0,0.D0))
      PARAMETER(EPS = 2.22D-16)
!
!----------------------------------------------------------------------!
!
      IERR = 0
      IF (.NOT.LISISO) THEN
         WRITE(*,*) 'cjkm_qual: Error anisotropy not an option!'
         IERR = 1
      ENDIF
!
!.... loop on integration points in x 
      DO 1 INTX=1,NINTX
         XI = XIGLL(INTX)
         DO 2 INTZ=1,NINTZ
! 
!.......... create the 2D shape functions at integration points 
            ETA = ETAGLL(INTZ)
            CALL CSF2D(2,NGNOD,XI,ETA, SF,DSF,IERR)
            IF (IERR.NE.0) THEN
               WRITE(*,*) 'cjmk_tvi: Error calling csf2d'
               RETURN
            ENDIF
! 
!.......... find the 2D shape fns and their derivatives 
            X = 0.D0
            Z = 0.D0
            XDEP = 0.D0
            ZDEP = 0.D0
            VP(INTX,INTZ) = 0.D0
            VS(INTX,INTZ) = 0.D0
            QPP(INTX,INTZ) = 0.D0
            QSS(INTX,INTZ) = 0.D0
            DXDXI = 0.D0
            DXDETA = 0.D0
            DZDXI = 0.D0
            DZDETA = 0.D0
            DO 4 IA=1,NGNOD !loop over anchor nodes   
! 
!............. interpolate x, z, rho, and elastic coefficients 
               ILOC = IENG(IA)
               X = X + SF(IA)*XLOCS(ILOC)
               Z = Z + SF(IA)*ZLOCS(ILOC)
               RLAM = ECOEFF(ILOC,1)
               RMU  = ECOEFF(ILOC,2)
               VPP = DSQRT( (RLAM + 2.D0*RMU)/DENS(ILOC) )
               VSS = DSQRT(RMU/DENS(ILOC))
               VP(INTX,INTZ)  = VP(INTX,INTZ) + SF(IA)*VPP
               VS(INTX,INTZ)  = VS(INTX,INTZ) + SF(IA)*VSS
               QPP(INTX,INTZ) = QPP(INTX,INTZ) + SF(IA)*QP(ILOC)
               QSS(INTX,INTZ) = QSS(INTX,INTZ) + SF(IA)*QS(ILOC)
               XDEP = XDEP + SF(IA)*XD(ILOC)
               ZDEP = ZDEP + SF(IA)*ZD(ILOC)
               DXDXI  = DXDXI  + DSF(1,IA)*XLOCS(ILOC)
               DXDETA = DXDETA + DSF(2,IA)*XLOCS(ILOC)
               DZDXI  = DZDXI  + DSF(1,IA)*ZLOCS(ILOC)
               DZDETA = DZDETA + DSF(2,IA)*ZLOCS(ILOC)
    4       CONTINUE
! 
!.......... check on jacobian 
            RJAC = DXDXI*DZDETA - DXDETA*DZDXI
            DET(INTX,INTZ) = RJAC
            IF (DABS(RJAC).LE.EPS .OR. RJAC.LT.0.D0) THEN
               WRITE(*,*) 'cjcm_qual: Jacobian undefined!'
               WRITE(*,*) 'cjcm_qual: Jacobian is ', RJAC
               IF (RJAC.LT.0.D0)
     ;         WRITE(*,*)'cjcm_qual: Jacobian can not be negative'
               IERR = 1
               RETURN
            ENDIF
! 
!.......... evaluate damping profile at location 
            VEL = VP(INTX,INTZ) 
            IF (XDEP.GT.EPS) THEN
               DAMPX(INTX,INTZ) = ZDAMPFN(XDEP,VEL,XWIDTH, DMAX_IN,
     ;                                    RCOEFF,RKAPPA,FCENT,OMEGA)
            ELSE
               DAMPX(INTX,INTZ) = CONE
            ENDIF
            IF (ZDEP.GT.EPS) THEN
               DAMPZ(INTX,INTZ) = ZDAMPFN(ZDEP,VEL,ZWIDTH, DMAX_IN,
     ;                                    RCOEFF,RKAPPA,FCENT,OMEGA)
            ELSE
               DAMPZ(INTX,INTZ) = CONE
            ENDIF
! 
!.......... loop over element nodes and fill out global shape fns 
            DO 6 IAE=1,NEN
               DRNXI  = SHL(1,IAE,INTX,INTZ) !dN_a/d_xi
               DRNETA = SHL(2,IAE,INTX,INTZ) !dN_a/d_eta
               SHG(3,IAE,INTX,INTZ) = SHL(3,IAE,INTX,INTZ)
               SHG(1,IAE,INTX,INTZ) = (DRNXI*DZDETA - DRNETA*DZDXI)/RJAC
               SHG(2,IAE,INTX,INTZ) =-(DRNXI*DXDETA - DRNETA*DXDXI)/RJAC
    6       CONTINUE
    2    CONTINUE !end loop on integration points in z 
    1 CONTINUE !end loop on integration points in x
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !

      SUBROUTINE FDAMP1(RN,ALPHA,DELTA, RKAP,RGAM)
      REAL*8 RN,ALPHA,DELTA, RKAP,RGAM, R
      PARAMETER(R = 0.001D0)
      RKAP = 1.D0
      RGAM = RN**2*ALPHA*3.D0/(2.D0*DELTA**3)*DLOG10(1.D0/R)
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE FDAMP3(RN,VP,DELTA,OMEGA, RKAP,RGAM)
      REAL*8 RN,VP,DELTA,OMEGA,FCENT, RKAP,RGAM
      REAL*8 RKAPPA,DMAX,D,ALPHA
c     real*8 d1,d2,beta
      REAL*8 RCOEFF, PI , TWOPI
      PARAMETER(RCOEFF = 1.D-3, PI = 3.14159265358979323846D0, 
     ;          TWOPI = 6.2831853071795862D0, FCENT = 0.D0)
      PARAMETER(N1 = 3, N2 = 0, N3 = 1) 
c     real*8 amax,bmax,kmax,beta
c     REAL*8 RKAPPA /1.00D0/ !could be 1.1 if necessary

      DMAX =-3.D0*VP/dsqrt(3.d0)*DLOG(RCOEFF)/(2.D0*DELTA) !log10(1/r) =-log10(r)
      D = DMAX*(RN/DELTA)**3
      ALPHA = PI*fcent*(1.D0 - RN/DELTA)!DELTA - RN) 
      rkappa = 1.d0 !+ (2.5d0 - 1.d0)*(rn/delta)**2
      RKAP = RKAPPA + ALPHA/(ALPHA**2 + OMEGA**2)
      RGAM =-OMEGA*D/(ALPHA**2 + OMEGA**2)
c     rkap = omega**2/(omega**2 + d**2)
c     rgam =-omega*d/(omega**2 + d**2)

!     RKAP = 1.D0
!     RGAM = RN**2*VP*3.D0/(2.D0*DELTA**3)*DLOG10(1.D0/Rcoeff)


c     amax = twopi*fcent
c     bmax =-dfloat(1 + n1 + n2)*vp*dlog10(RCOEFF)/(2.D0*delta)
c     kmax =-dfloat(1 + n1)*10.d0*delta*dlog10(RCOEFF)/(2.d0*delta)
c     rkappa = 1.d0 + kmax*(rn/delta)**n1
c     beta = bmax*(rn/delta)**(n1 + n2)
c     alpha = amax*( (delta - rn)/delta )**n3
c     rkap = rkappa + alpha*beta/(alpha**2 + omega**2)
c     rgam =-omega*beta/(alpha**2 + omega**2)
c     return

c     d1 = 1.d0 + .5d0*dsin(rn*pi/(2.d0*delta))
c     d2 = 11.d0*dsin(rn*pi/(2.d0*delta))
c     rkap = omega**2*d1/(omega**2*d1**2 + d2**2)
c     rgam =-omega   *d2/(omega**2*d1**2 + d2**2)
      RETURN
      END 

      REAL*8 FUNCTION DIN(N,ALPHA,DELTA,R)
      REAL*8 N,ALPHA,DELTA,R
!     DIN = N**2*ALPHA*3.D0/(2.D0*DELTA**3)*DLOG(1.D0/R) 
      DIN = N**2*ALPHA*3.D0/(2.D0*DELTA**3)*DLOG10(1.D0/R)
      RETURN
      END

      COMPLEX*16 FUNCTION ZDAMPFN(RN,VP,DELTA, DMAX_IN,
     ;                            RCOEFF,RKAPPA,FCENT,OMEGA)
!
!     Pick a damping function and implement it, mix and match as you
!     please
!
!     INPUT      MEANING
!     -----      ------- 
!     DELTA      PML width
!     DMAX_IN    an override for max damping
!     FCENT      center frequency
!     OMEGA      angular frequency 
!     RN         distance into PML  
!     RCOEFF     theoretical refleciton coefficient
!     RKAPPA     seems to work with damping
!     VP         Vp velocity 
!
!     OUTPUT     MEANING
!     ------     -------
!     ZDAMPFN    damping function evaluated at RN 
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: RN, VP, DELTA, OMEGA, RKAPPA, FCENT, 
     ;                      RCOEFF, DMAX_IN 
      REAL*8 DMAX, D, ALPHA
!.... local variables
      REAL*8 PI , TWOPI
      INTEGER*4 N1, N2, N3
      PARAMETER(PI     = 3.1415926535897931D0, 
     ;          TWOPI  = 6.2831853071795862D0) 
      PARAMETER(N1 = 3, N2 = 0, N3 = 1)  
!
!----------------------------------------------------------------------1
! 
      IF (DMAX_IN > 0.0) THEN
         DMAX = DMAX_IN 
      ELSE
         DMAX =-3.D0*VP*DLOG(RCOEFF)/(2.D0*DELTA)
      ENDIF
      D = DMAX*(RN/DELTA)**2
      ALPHA = PI*FCENT*(1.D0 - RN/DELTA)
      ZDAMPFN = DCMPLX(RKAPPA,0.D0) 
     ;        + DCMPLX(DABS(D),0.D0)/DCMPLX(ALPHA,OMEGA)
      ZDAMPFN = DCMPLX(1.D0,0.D0)/ZDAMPFN
      RETURN
c     R = 0.D0
c     Z = 0.D0
c     IF (IOPT.EQ.2) THEN !brossier
c        D1 = 1.D0 + 0.5D0*DCOS(RN*PI/(2.D0*DELTA)) !0.5d0
c        D2 = 11.D0*DCOS(RN*PI/(2.D0*DELTA)) !11.d0
c        !R = OMEGA**2*D1/(OMEGA**2*D1**2 + D2**2)
c        !Z =-OMEGA   *D2/(OMEGA**2*D1**2 + D2**2)
c        zdampfn = dcmplx(1.d0,0.d0)/dcmplx(d1,-d2/omega)
c        return
c     ELSEIF (IOPT.EQ.2) THEN !brossier modified
c        D1 = 1.D0 + 0.5D0*DCOS(RN*PI/(2.D0*DELTA))
c        D2 = 11.D0*DSIN(RN*PI/(2.D0*DELTA)) 
c        R = OMEGA**2*D1/(OMEGA**2*D1**2 + D2**2)
c        Z =-OMEGA   *D2/(OMEGA**2*D1**2 + D2**2)
c     ELSEIF (IOPT.EQ.3) THEN !more general komatisch, dont remember
c        AMAX = TWOPI*FCENT
c        BMAX =-DFLOAT(1 + N1 + N2)*VP*DLOG10(RCOEFF)/(2.D0*DELTA)
c        KMAX =-DFLOAT(1 + N1)*10.D0*DELTA*DLOG10(RCOEFF)/(2.D0*DELTA)
c        RKAPPA = 1.D0 + KMAX*(RN/DELTA)**N1
c        BETA = BMAX*(RN/DELTA)**(N1 + N2)
c        ALPHA = AMAX*( (DELTA - RN)/DELTA )**N3
c        R = RKAPPA + ALPHA*BETA/(ALPHA**2 + OMEGA**2)
c        Z =-OMEGA*BETA/(ALPHA**2 + OMEGA**2)
c     ELSEIF (IOPT.EQ.4) THEN
c        R = 1.D0 + 0.d0*(RN/DELTA)
c        Z =-20.d0*(RN/DELTA) / (OMEGA*DELTA/(VP/DSQRT(3.D0)))
c     ELSEIF (IOPT.EQ.5) THEN
c        PMLF = 3.D0*DLOG(1.D0/RCOEFF)/(2.D0*DELTA**3) 
c        ALPHA = TWOPI*FCENT*(1.D0 - RN/DELTA)
c        RKAPPA = 1.D0 + (KAPPA - 1.D0)*(RN/DELTA)**2
c        D = RN**2*VP*PMLF 
c        ZDAMPFN = DCMPLX(ALPHA,-OMEGA) 
c    ;            /DCMPLX(RKAPPA*ALPHA + D,-RKAPPA*OMEGA)
c        zdampfn = dconjg(zdampfn) 
c        !zdampfn = dcmplx(1.d0,0.d0)/zdampfn
c        RETURN
c     ELSE !komatisch, no convolution  
         IF (DMAX_IN > 0.0) THEN
            DMAX = DMAX_IN 
         ELSE
            DMAX =-3.D0*VP*DLOG(RCOEFF)/(2.D0*DELTA)
         ENDIF
c        dmax = 2.d0
         D = DMAX*(RN/DELTA)**2
         ALPHA = PI*FCENT*(1.D0 - RN/DELTA)
c        RKAPPA = 1.5D0
c        R = RKAPPA + D*ALPHA/(ALPHA**2 + OMEGA**2)
c        Z =-OMEGA*D/(ALPHA**2 + OMEGA**2)
c        !zdampfn = dcmplx(1.d0,0.d0)/dcmplx(1.d0,-dabs(d)/omega)
c        zdampfn = dcmplx(1.d0,0.d0)/dcmplx(rkappa+alpha,-dabs(d)/omega)
         ZDAMPFN = DCMPLX(RKAPPA,0.D0) 
     ;           + DCMPLX(DABS(D),0.D0)/DCMPLX(ALPHA,OMEGA)
         ZDAMPFN = DCMPLX(1.D0,0.D0)/ZDAMPFN
         RETURN
c     ENDIF
c     ZDAMPFN = DCMPLX(R,-Z)
      RETURN
      END
