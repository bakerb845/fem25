      SUBROUTINE DJACAB(MEN,MINTX,MINTZ,MNPG, LISISO,
     ;                  NINTX,NINTZ,NEN,NGNOD,NNPG, INPG, IENG,
     ;                  XIGLL,ETAGLL,XLOCS,ZLOCS,DENS,ECOEFF,SHL, 
     ;                  DET,DAB,SHG,IERR) 
!
!     Jacobian, global shape functions, and material derivatives in 
!     alpha and beta at integration points. The idea is that we 
!     collocate on the inversion nodal point INPG and calculate its 
!     value at an integration abscissa with the hat shape functions.  
!     - B. Baker Dec 2012
!
!     INPUT      INPUT  
!     -----      ----- 
!     DENS       densities at the global nodes
!     ECOEFF     elastic coefficients at anchor nodes 
!     ETAGLL     holds integration points/weights in eta
!     IENG       global anchor node pointer
!     INPG       nodal point to collocate derivative on
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
!
!     OUTPUT     MEANING
!     ------     ------- 
!     DAB        INPG alpha and beta derivative at integration points 
!     DET        jacobian at integration points 
!     IERR       error flag
!     SHG        global shape fns and derivatives at integration points
!
!.... variable declarations 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8, INTENT(IN) :: SHL(3,MEN,MINTX,*), ECOEFF(MNPG,*), 
     ;                      XLOCS(NNPG), ZLOCS(NNPG), DENS(NNPG), 
     ;                      XIGLL(NINTX), ETAGLL(NINTZ) 
      INTEGER*4, INTENT(IN) :: IENG(NGNOD), MEN,MINTX,MINTZ,MNPG, 
     ;                         NINTX,NINTZ,NEN,NGNOD,NNPG, INPG
      LOGICAL*4, INTENT(IN) :: LISISO
      REAL*8, INTENT(OUT) :: SHG(3,MEN,MINTX,*), DAB(MINTX,MINTZ,*), 
     ;                       DET(MINTX,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: DSF(:,:), SF(:)
      REAL*8 XI,ETA,DXDXI,DXDETA,DZDXI,DZDETA, TWO  
      PARAMETER(EPS = 2.22D-16, TWO = 2.D0)
!
!----------------------------------------------------------------------!
!
!.... initialize and check
      IERR = 0 
      IF (.NOT.LISISO) THEN
         WRITE(*,*)  'djacab: You cant call djacab without isotropy'
         IERR = 1
         RETURN
      ENDIF
      ALLOCATE(SF(NGNOD))
      ALLOCATE(DSF(2,NGNOD))
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
               WRITE(*,*) 'djacab: Error calling csf2d'
               RETURN
            ENDIF
! 
!.......... find the 2D shape fns and their derivatives 
            X = 0.D0
            Z = 0.D0
            DXDXI = 0.D0
            DXDETA = 0.D0
            DZDXI = 0.D0
            DZDETA = 0.D0
            DAB(INTX,INTZ,1) = 0.D0
            DAB(INTX,INTZ,2) = 0.D0
            DO 4 IA=1,NGNOD !loop over anchor nodes 
! 
!............. interpolate x, z, rho, and elastic coefficients 
               ILOC = IENG(IA)
               IF (ILOC.EQ.INPG) THEN
                  ALPHA = DSQRT( (ECOEFF(ILOC,1) + 2.D0*ECOEFF(ILOC,2))
     ;                          /DENS(ILOC) )
                  BETA  = DSQRT( ECOEFF(ILOC,2)/DENS(ILOC) ) 
               ELSE 
                  ALPHA = 0.D0
                  BETA  = 0.D0 
               ENDIF 
               DAB(INTX,INTZ,1) = DAB(INTX,INTZ,1) + TWO*ALPHA*SF(IA)
               DAB(INTX,INTZ,2) = DAB(INTX,INTZ,2) + TWO*BETA *SF(IA)
               X = X + SF(IA)*XLOCS(ILOC)
               Z = Z + SF(IA)*ZLOCS(ILOC)
               DXDXI  = DXDXI  + DSF(1,IA)*XLOCS(ILOC)
               DXDETA = DXDETA + DSF(2,IA)*XLOCS(ILOC)
               DZDXI  = DZDXI  + DSF(1,IA)*ZLOCS(ILOC)
               DZDETA = DZDETA + DSF(2,IA)*ZLOCS(ILOC)
    4       CONTINUE
!           IF (DABS(XLOCS(INPG) - X).GT.1.D-7 .OR.  
!    ;          DABS(ZLOCS(INPG) - Z).GT.1.D-7) DAB(INTX,INTZ,1:2)=0.D0
!
!.......... check on jacobian
            RJAC = DXDXI*DZDETA - DXDETA*DZDXI
            DET(INTX,INTZ) = RJAC
            IF (DABS(RJAC).LE.EPS .OR. RJAC.LT.0.D0) THEN
               WRITE(*,*) 'djacab: Jacobian undefined!'
               WRITE(*,*) 'djacab: Jacobian is ', RJAC
               IF (RJAC.LT.0.D0) 
     ;         WRITE(*,*) 'djacab: Jacobian can not be negative'
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
      DEALLOCATE(SF)
      DEALLOCATE(DSF)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DJACAB2(MEN,MINTX,MINTZ,MNPG, LISISO,
     ;                   NINTX,NINTZ,NEN,NGNOD,NNPG, INPG, IENG,
     ;                   XIGLL,ETAGLL,XLOCS,ZLOCS, 
     ;                   DENS,ECOEFF,SHL,
     ;                   DET, DAB,DAB2,SHG,IERR) 
!
!     Jacobian, global shape functions, and material first and 2nd 
!     derivatives in alpha and beta at integration points. The idea is 
!     that we collocate on the inversion nodal point INPG then calculate 
!     its value at an integration abscissa with the hat shape functions.  
!     - B. Baker April 2013
!
!     INPUT      INPUT  
!     -----      ----- 
!     DENS       density at anchor nodes
!     ECOEFF     lambda, mu at anchor nodes
!     ETAGLL     holds integration points/weights in eta
!     IENG       global anchor node pointer
!     INPG       nodal point to collocate derivative on
!     LISISO     True -> simulation isotrpic, False -> anisotropic
!     MEN        max element nodes
!     MINTX      max integration points in x 
!     MINTZ      max integration points in z 
!     NEN        number of element nodes 
!     NGNOD      number of global element nodes (anchor nodes, 4 or 9)
!     NINTX      number of integration points in x 
!     NINTZ      number of integration points in z 
!     NNPG       number of nodal points in global mesh
!     SHL        lagrange shape function array 
!     XIGLL      holds integration points/weights in xi 
!     XLOCS      physical x locations of anchor nodes 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     DAB        1st derivatives of INPG alpha and beta derivative 
!                at integration points
!     DAB2       2nd derivatives of INPG alpha and beta derivative 
!                at integration points 
!     DET        jacobian at integration points 
!     IERR       error flag
!     SHG        global shape fns and derivatives at integration points
!
!.... variable declarations 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8, INTENT(IN) :: SHL(3,MEN,MINTX,*), ECOEFF(MNPG,*), 
     ;                      DENS(NNPG), XLOCS(NNPG),ZLOCS(NNPG),
     ;                      XIGLL(NINTX), ETAGLL(NINTZ)
      INTEGER*4, INTENT(IN) :: IENG(NGNOD), MEN,MINTX,MINTZ,MNPG,
     ;                         NINTX,NINTZ,NEN,NGNOD,NNPG, INPG
      LOGICAL*4, INTENT(IN) :: LISISO
      REAL*8, INTENT(OUT) :: SHG(3,MEN,MINTX,*), DAB(MINTX,MINTZ,*), 
     ;                       DAB2(MINTX,MINTZ,*), DET(MINTX,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: DSF(:,:), SF(:)
      REAL*8 XI,ETA,DXDXI,DXDETA,DZDXI,DZDETA, DRNXI,DRNETA, 
     ;       ALPHA,BETA, DALPHA2,DBETA2, RJAC, EPS, TWO
      INTEGER*4 INTX, INTZ, ILOC, IAE, IA
      PARAMETER(EPS = 2.22D-16, TWO = 2.D0)
!
!----------------------------------------------------------------------!
!
!.... initialize and check
      IERR = 0
      IF (.NOT.LISISO) THEN
         WRITE(*,*)  'djacab2: You cant call djacab2 without isotropy'
         IERR = 1
         RETURN
      ENDIF
      ALLOCATE(SF(NGNOD))
      ALLOCATE(DSF(2,NGNOD))
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
               WRITE(*,*) 'djacab2: Error calling csf2d'
               RETURN
            ENDIF
! 
!.......... find the 2D shape fns and their derivatives 
            DXDXI = 0.D0
            DXDETA = 0.D0
            DZDXI = 0.D0
            DZDETA = 0.D0
            DAB(INTX,INTZ,1) = 0.D0
            DAB(INTX,INTZ,2) = 0.D0
            DAB2(INTX,INTZ,1) = 0.D0
            DAB2(INTX,INTZ,2) = 0.D0
            DO 4 IA=1,NGNOD !loop over anchor nodes 
! 
!............. interpolate x, z, rho, and elastic coefficients 
               ILOC = IENG(IA)
               IF (ILOC.EQ.INPG) THEN
                  DALPHA2 = 2.D0
                  DBETA2  = 2.D0
                  ALPHA = DSQRT( (ECOEFF(ILOC,1) + 2.D0*ECOEFF(ILOC,2))
     ;                          /DENS(ILOC) )
                  BETA  = DSQRT( ECOEFF(ILOC,2)/DENS(ILOC) )
               ELSE
                  DALPHA2 = 0.D0
                  DBETA2  = 0.D0
                  ALPHA   = 0.D0
                  BETA    = 0.D0
               ENDIF
               DAB (INTX,INTZ,1) = DAB (INTX,INTZ,1) + TWO*ALPHA*SF(IA)
               DAB (INTX,INTZ,2) = DAB (INTX,INTZ,2) + TWO*BETA *SF(IA)
               DAB2(INTX,INTZ,1) = DAB2(INTX,INTZ,1) + DALPHA2*SF(IA)
               DAB2(INTX,INTZ,2) = DAB2(INTX,INTZ,2) + DBETA2 *SF(IA)
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
               WRITE(*,*) 'djacab2: Jacobian undefined!'
               WRITE(*,*) 'djacab2: Jacobian is ', RJAC
               IF (RJAC.LT.0.D0)
     ;         WRITE(*,*) 'djacab2: Jacobian can not be negative'
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
      DEALLOCATE(SF)
      DEALLOCATE(DSF)
      RETURN
      END

