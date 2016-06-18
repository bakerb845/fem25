      SUBROUTINE KM25ISOQ(MEN,MEE,MINTX, NEN,NINTX,NINTZ, 
     ;                    OMEGA,PY, XIWTS,ETAWTS, VLAM,VMU,
     ;                    RHO,DET,SHG, ME,KE)
! 
!     Calculates all the element stiffness and mass matrices for 
!     the 2.5D elastic wave problem.  This is for isotropic 
!     elements. 
! 
!     INPUT      MEANING 
!     -----      ------- 
!     DET        holds the jacobian at the integration points
!     ETAWTS     holds integration weights in eta
!     MEE        max number of element equations
!     MEN        max number of element nodes
!     NEN        number of element nodes
!     NINTX      number of integration points in x 
!     NINTZ      number of integration points in z  
!     OMEGA      angular frequency
!     RHO        density at integration points
!     SHG        global shape fns at int. pts for each node
!     VLAM       holds the lambda lame parameter at integration points
!     VMU        holds the mu lame parameter at integration points
!     XIWTS      holds integration weights in xi
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     KE         element stiffness matrix
!     ME         element mass matrix 
!  
!.... variable declarations 
      REAL*8 SHG(3,MEN,MINTX,*), VLAM(MINTX,*), VMU(MINTX,*), 
     ;       RHO(MINTX,*), DET(MINTX,*), XIWTS(NINTX), ETAWTS(NINTZ),
     ;       OMEGA, PY
      COMPLEX*16 KE(MEE,*)
      REAL*8 ME(MEE,*)
!.... local variables
      COMPLEX*16 CZERO
      REAL*8 WPY, PY2, W2, WEIGHT, LAMBDA, MU, RHOW, 
     ;       PHIB, PHIBX, PHIBZ, PHIA, PHIAX, PHIAZ, 
     ;       KE11, KE12, KE13, KE21, KE22, KE23, KE31, KE32, KE33,
     ;       COEFFM1, COEFFM2, COEFFM3, PHIBR1,PHIBR2,PHIBR3 
      REAL*8 C11,C12,C13, C22,C23, C33, C44,C55,C66
      REAL*8 KE11T1,KE11T2, KE12T1,KE12T2, KE13T1,KE13T2,
     ;       KE21T1,KE21T2, KE22T1,KE22T2, KE23T1,KE23T2, 
     ;       KE31T1,KE31T2, KE32T1,KE32T2, KE33T1,KE33T2
      LOGICAL*4 LHERM
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0))
      PARAMETER(NDIM = 3) 
      PARAMETER(LHERM = .FALSE.) 
! 
!----------------------------------------------------------------------!
! 
!.... null out element matrices
      NEE = NEN*NDIM
      DO 1 IPE=1,NEE
         DO 2 IQE=IPE,NEE
            ME(IPE,IQE) = 0.D0
            KE(IPE,IQE) = CZERO
            ME(IQE,IPE) = 0.D0
            KE(IQE,IPE) = CZERO
    2    CONTINUE
    1 CONTINUE 
! 
!.... loop on integration points
      WPY = OMEGA*PY
      PY2 = PY**2
      W2 = OMEGA**2
      DO 3 INTX=1,NINTX
         DO 4 INTZ=1,NINTZ
            WEIGHT = XIWTS(INTX)*ETAWTS(INTZ)*DET(INTX,INTZ)
            LAMBDA = VLAM(INTX,INTZ)*WEIGHT
            MU = VMU(INTX,INTZ)*WEIGHT

            C11 = LAMBDA + 2.D0*MU
            C12 = LAMBDA
            C13 = LAMBDA

            C22 = LAMBDA + 2.D0*MU
            C23 = LAMBDA

            C33 = LAMBDA + 2.D0*MU

            C44 = MU
            C55 = C44
            C66 = C44 

            RHOW = RHO(INTX,INTZ)*WEIGHT
            COEFFM1 =-W2*(RHOW - PY2*C66)
            COEFFM2 =-W2*(RHOW - PY2*C22)
            COEFFM3 =-W2*(RHOW - PY2*C44)
! 
!.......... loop on element nodes
            DO 5 IBE=1,NEN
               PHIB  = SHG(3,IBE,INTX,INTZ)
               PHIBX = SHG(1,IBE,INTX,INTZ)
               PHIBZ = SHG(2,IBE,INTX,INTZ)
               PHIBR1 = PHIB*COEFFM1
               PHIBR2 = PHIB*COEFFM2
               PHIBR3 = PHIB*COEFFM3

               KE11T1 = C11*PHIBX
               KE11T2 = C55*PHIBZ
               KE13T1 = C13*PHIBZ
               KE13T2 = C55*PHIBX
               KE22T1 = C66*PHIBX
               KE22T2 = C44*PHIBZ
               KE31T1 = C13*PHIBX 
               KE31T2 = C55*PHIBZ
               KE33T1 = C33*PHIBZ 
               KE33T2 = C55*PHIBX

               KE12T1 = WPY*C12*PHIB
               KE12T2 = WPY*C66*PHIBX

               KE21T1 = WPY*C12*PHIBX 
               KE21T2 = WPY*C66*PHIB
               KE23T1 = WPY*C23*PHIBZ
               KE23T2 = WPY*C44*PHIB

               KE32T1 = WPY*C23*PHIB 
               KE32T2 = WPY*C44*PHIBZ

               IF (LHERM) THEN
                  IAE1 = IBE
               ELSE
                  IAE1 = 1
               ENDIF 
               DO 6 IAE=IAE1,NEN 
                  PHIA  = SHG(3,IAE,INTX,INTZ)
                  PHIAX = SHG(1,IAE,INTX,INTZ)
                  PHIAZ = SHG(2,IAE,INTX,INTZ)
                  !indices
                  IPE = NDIM*(IAE - 1) + 1
                  IQE = NDIM*(IBE - 1) + 1 
                  !mass matrix
                  ME(IPE  ,IQE  ) = ME(IPE  ,IQE  ) + PHIBR1*PHIA
                  ME(IPE+1,IQE+1) = ME(IPE+1,IQE+1) + PHIBR2*PHIA
                  ME(IPE+2,IQE+2) = ME(IPE+2,IQE+2) + PHIBR3*PHIA
                  !stiffness, real
                  KE11 = PHIAX*KE11T1 + PHIAZ*KE11T2
                  KE13 = PHIAX*KE13T1 + PHIAZ*KE13T2
                  KE22 = PHIAX*KE22T1 + PHIAZ*KE22T2
                  KE31 = PHIAZ*KE31T1 + PHIAX*KE31T2
                  KE33 = PHIAZ*KE33T1 + PHIAX*KE33T2
                  !stiffness, imaginary 
                  KE12 = PHIAX*KE12T1 - PHIA *KE12T2 
                  KE21 = PHIA *KE21T1 - PHIAX*KE21T2
                  KE23 = PHIA *KE23T1 - PHIAZ*KE23T2 
                  KE32 = PHIAZ*KE32T1 - PHIA *KE32T2

                  KE(IPE  ,IQE  ) = KE(IPE  ,IQE  ) + DCMPLX(KE11,0.D0)
                  KE(IPE  ,IQE+1) = KE(IPE  ,IQE+1) - DCMPLX(0.D0,KE12)
                  KE(IPE  ,IQE+2) = KE(IPE  ,IQE+2) + DCMPLX(KE13,0.D0)
                  KE(IPE+1,IQE  ) = KE(IPE+1,IQE  ) + DCMPLX(0.D0,KE21)
                  KE(IPE+1,IQE+1) = KE(IPE+1,IQE+1) + DCMPLX(KE22,0.D0)
                  KE(IPE+1,IQE+2) = KE(IPE+1,IQE+2) + DCMPLX(0.D0,KE23)
                  KE(IPE+2,IQE  ) = KE(IPE+2,IQE  ) + DCMPLX(KE31,0.D0)
                  KE(IPE+2,IQE+1) = KE(IPE+2,IQE+1) - DCMPLX(0.D0,KE32)
                  KE(IPE+2,IQE+2) = KE(IPE+2,IQE+2) + DCMPLX(KE33,0.D0)
    6          CONTINUE !loop on nodes in iae
    5       CONTINUE !loop on nodes in ibe
    4    CONTINUE !loop on integration points in eta 
    3 CONTINUE !loop on integration points in xi
! 
!.... fill in lower half
      IF (LHERM) THEN
         DO 10 IPE=1,NEE
            DO 11 IQE=IPE-1,1,-1
               ME(IPE,IQE) = ME(IQE,IPE)
               KE(IPE,IQE) = DCONJG(KE(IQE,IPE))
   11       CONTINUE
   10    CONTINUE
      ELSE
!        DO 20 IPE=1,NEE
!           DO 21 IQE=IPE-1,1,-1
!              XDIF = CDABS( KE(IPE,IQE) - DCONJG(KE(IQE,IPE)) )
!    ;               /CDABS(KE(IPE,IQE))  
!              IF (XDIF.GT.1.D-10)
!    ;         WRITE(*,*) 'lherm not right',ke(ipe,iqe),ke(iqe,ipe)
!  21       CONTINUE
!  20    CONTINUE 
c        do ipe=1,nee
c           do iqe=ipe+1,nee
c              rtemp = me(ipe,iqe) 
c              ztemp = ke(ipe,iqe)
c              me(ipe,iqe) = me(iqe,ipe)
c              me(iqe,ipe) = rtemp
c              ke(ipe,iqe) = ke(iqe,ipe)
c              ke(iqe,ipe) = ztemp
c           enddo
c       enddo
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CM25ISOQ(MEN,MEE,MINTX, NEN,NINTX,NINTZ, 
     ;                    OMEGA,PY, XIWTS,ETAWTS, VLAM,VMU,
     ;                    RHO,DET, SHG, GAMMAX,GAMMAZ, ME,CE)
! 
!     Calculates all the element damping and mass matrices for 
!     the 2.5D elastic wave problem.  This is for isotropic elements
!     B. Baker - Demember 2012 
!
!     This is now templated for anisotropy and a bug was fixed, 
!     i w py was set as w py.  This almost caused me to hulk smash 
!     something. B. Baker - January 2013 
! 
!     INPUT      MEANING 
!     -----      ------- 
!     DET        holds the jacobian at the integration points
!     ETAWTS     holds integration weights in eta
!     GAMMAX     1/gx damping function at integration points
!     GAMMAZ     1/gz damping function at integration points
!     MEE        max number of element equations
!     MEN        max number of element nodes
!     NEN        number of element nodes
!     MINTX      leading dimension 
!     NINTX      number of integration points in x 
!     NINTZ      number of integration points in z  
!     OMEGA      angular frequency
!     PY         aparrant slowness in y (s/m)
!     RHO        density at integration points
!     SHG        global shape fns at int. pts for each node
!     VLAM       holds the lambda lame parameter at integration points
!     VMU        holds the mu lame parameter at integration points
!     XIWTS      holds integration weights in xi
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     KE         element stiffness matrix
!     ME         element mass matrix 
!  
!.... variable declarations 
      IMPLICIT NONE
      COMPLEX*16, INTENT(IN) :: GAMMAX(MINTX,*), GAMMAZ(MINTX,*)
      REAL*8, INTENT(IN) :: SHG(3,MEN,MINTX,*), VLAM(MINTX,*), 
     ;        VMU(MINTX,*), RHO(MINTX,*), DET(MINTX,*), XIWTS(NINTX), 
     ;        ETAWTS(NINTZ), OMEGA, PY
      INTEGER*4, INTENT(IN) :: MEN,MEE,MINTX, NEN,NINTX,NINTZ 
      COMPLEX*16 CE(MEE,*)
      REAL*8 ME(MEE,*)
!.... local variables
      COMPLEX*16 PHIBX, PHIBZ, PHIAX, PHIAZ,
     ;           CE11, CE12, CE13, CE21, CE22, CE23, CE31, CE32, CE33,
     ;           CZERO, CONE, GX, GZ, GX2, GXZ, GZ2, CPHIB, CPHIA 
!    ;           ZPHIBR1, ZPHIBR2, ZPHIBR3
      COMPLEX*16 C11,C12,C13, C22,C23,C33, C44,C55,C66
      COMPLEX*16 KE11T1,KE11T2,KE12T1,KE12T2,KE13T1,KE13T2,
     ;           KE21T1,KE21T2,KE22T1,KE22T2,KE23T1,KE23T2,
     ;           KE31T1,KE31T2,KE32T1,KE32T2,KE33T1,KE33T2, CWPY !, 
c    ;           ME11, ME22, ME33
      REAL*8 PY2, W2, WEIGHT, LAMBDA, MU, RHOW, COEFFM1, COEFFM2, 
     ;       COEFFM3, PHIB, PHIA, PHIBR1,PHIBR2,PHIBR3, 
     ;       RC11, RC12, RC13, RC22, RC23, RC33, RC44, RC55, RC66, 
     ;       WPY  
      INTEGER*4 NDIM, NEE, IAE, IBE, INTX, INTZ, IPE, IQE 
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0), 
     ;          CONE  = DCMPLX(1.D0,0.D0))
      PARAMETER(NDIM = 3)  
! 
!----------------------------------------------------------------------!
! 
!.... null out element matrices
      NEE = NEN*NDIM
      DO 1 IPE=1,NEE
         DO 2 IQE=1,NEE
            ME(IPE,IQE) = 0.D0
            CE(IPE,IQE) = CZERO 
    2    CONTINUE
    1 CONTINUE 
! 
!.... loop on integration points 
      WPY = OMEGA*PY 
      CWPY = DCMPLX(0.D0,WPY)
      PY2 = PY**2
      W2 = OMEGA**2
      DO 3 INTX=1,NINTX 
         DO 4 INTZ=1,NINTZ 
            WEIGHT = XIWTS(INTX)*ETAWTS(INTZ)*DET(INTX,INTZ)
            GX = GAMMAX(INTX,INTZ) !already 1/gx
            GZ = GAMMAZ(INTX,INTZ) !already 1/gz 
            GX2 = GX**2
            GXZ = GX*GZ 
            GZ2 = GZ**2 

            LAMBDA = VLAM(INTX,INTZ)*WEIGHT 
            MU = VMU(INTX,INTZ)*WEIGHT 

            RC11 = LAMBDA + 2.D0*MU
            RC12 = LAMBDA
            RC13 = RC12
   
            RC22 = RC11
            RC23 = RC13
 
            RC33 = RC11
  
            RC44 = MU
            RC55 = RC44
            RC66 = RC44

            C11 = DCMPLX(RC11,0.D0)
            C12 = DCMPLX(RC12,0.D0)
            C13 = DCMPLX(RC13,0.D0)
            C22 = DCMPLX(RC22,0.D0)
            C23 = DCMPLX(RC23,0.D0)
            C33 = DCMPLX(RC33,0.D0)
            C44 = DCMPLX(RC44,0.D0)
            C55 = DCMPLX(RC55,0.D0)
            C66 = DCMPLX(RC66,0.D0)

            RHOW = RHO(INTX,INTZ)*WEIGHT
            COEFFM1 =-W2*(RHOW - PY2*RC66)
            COEFFM2 =-W2*(RHOW - PY2*RC22)
            COEFFM3 =-W2*(RHOW - PY2*RC44)
            DO 5 IBE=1,NEN
               PHIB  = SHG(3,IBE,INTX,INTZ)
               PHIBX = DCMPLX(SHG(1,IBE,INTX,INTZ),0.D0)*GX
               PHIBZ = DCMPLX(SHG(2,IBE,INTX,INTZ),0.D0)*GZ
               CPHIB = DCMPLX(PHIB,0.D0)
               PHIBR1 = COEFFM1*PHIB 
               PHIBR2 = COEFFM2*PHIB
               PHIBR3 = COEFFM3*PHIB

               KE11T1 = C11*PHIBX
               KE11T2 = C55*PHIBZ
               KE13T1 = C13*PHIBZ
               KE13T2 = C55*PHIBX
               KE22T1 = C66*PHIBX
               KE22T2 = C44*PHIBZ
               KE31T1 = C13*PHIBX
               KE31T2 = C55*PHIBZ
               KE33T1 = C33*PHIBZ
               KE33T2 = C55*PHIBX
!
!............. only damp the real part
c              PHIBX = DCMPLX(SHG(1,IBE,INTX,INTZ),0.D0)
c              PHIBZ = DCMPLX(SHG(2,IBE,INTX,INTZ),0.D0)
               KE12T1 = CWPY*C12*CPHIB
               KE12T2 = CWPY*C66*PHIBX
               KE21T1 = CWPY*C12*PHIBX
               KE21T2 = CWPY*C66*CPHIB
               KE23T1 = CWPY*C23*PHIBZ
               KE23T2 = CWPY*C44*CPHIB
               KE32T1 = CWPY*C23*CPHIB
               KE32T2 = CWPY*C44*PHIBZ

               DO 6 IAE=1,NEN 
                  PHIA  = SHG(3,IAE,INTX,INTZ)
                  PHIAX = DCMPLX(SHG(1,IAE,INTX,INTZ),0.D0)*GX
                  PHIAZ = DCMPLX(SHG(2,IAE,INTX,INTZ),0.D0)*GZ
                  CPHIA = DCMPLX(PHIA,0.D0) 
                  !indices
                  IPE = NDIM*(IAE - 1) + 1
                  IQE = NDIM*(IBE - 1) + 1
                  !mass matrix

                  ME(IPE  ,IQE  ) = ME(IPE  ,IQE  ) + PHIA*PHIBR1
                  ME(IPE+1,IQE+1) = ME(IPE+1,IQE+1) + PHIA*PHIBR2
                  ME(IPE+2,IQE+2) = ME(IPE+2,IQE+2) + PHIA*PHIBR3

                  CE11 = PHIAX*KE11T1 + PHIAZ*KE11T2 !+ ME11
                  CE13 = PHIAX*KE13T1 + PHIAZ*KE13T2
                  CE22 = PHIAX*KE22T1 + PHIAZ*KE22T2 !+ ME22
                  CE31 = PHIAZ*KE31T1 + PHIAX*KE31T2
                  CE33 = PHIAZ*KE33T1 + PHIAX*KE33T2 !+ ME33
                  !stiffness, imaginary 
c                 PHIAX = DCMPLX(SHG(1,IAE,INTX,INTZ),0.D0)
c                 PHIAZ = DCMPLX(SHG(2,IAE,INTX,INTZ),0.D0)
                  CE12 = PHIAX*KE12T1 - CPHIA*KE12T2
                  CE21 = CPHIA*KE21T1 - PHIAX*KE21T2
                  CE23 = CPHIA*KE23T1 - PHIAZ*KE23T2
                  CE32 = PHIAZ*KE32T1 - CPHIA*KE32T2

                  !stiffness
                  CE(IPE  ,IQE  ) = CE(IPE  ,IQE  ) + CE11 
                  CE(IPE  ,IQE+1) = CE(IPE  ,IQE+1) - CE12
                  CE(IPE  ,IQE+2) = CE(IPE  ,IQE+2) + CE13

                  CE(IPE+1,IQE  ) = CE(IPE+1,IQE  ) + CE21
                  CE(IPE+1,IQE+1) = CE(IPE+1,IQE+1) + CE22 
                  CE(IPE+1,IQE+2) = CE(IPE+1,IQE+2) + CE23

                  CE(IPE+2,IQE  ) = CE(IPE+2,IQE  ) + CE31
                  CE(IPE+2,IQE+1) = CE(IPE+2,IQE+1) - CE32
                  CE(IPE+2,IQE+2) = CE(IPE+2,IQE+2) + CE33 
    6          CONTINUE !loop on nodes in iae
    5       CONTINUE !loop on nodes in ibe
    4    CONTINUE !loop on integration points in eta 
    3 CONTINUE !loop on integration points in xi
! 
!.... fill in lower half
c     DO 10 IPE=1,NEE
c        DO 11 IQE=IPE-1,1,-1
c           ME(IPE,IQE) = ME(IQE,IPE)
c           CE(IPE,IQE) = CE(IQE,IPE)
c           ME(IPE,IQE) = ME(IQE,IPE)
c           CE(IPE,IQE) = dconjg(CE(IQE,IPE))
c  11    CONTINUE
c  10 CONTINUE

c        do ipe=1,nee
c           do iqe=ipe+1,nee
c              rtemp = me(ipe,iqe) 
c              ztemp = ce(ipe,iqe)
c              me(ipe,iqe) = me(iqe,ipe)
c              me(iqe,ipe) = rtemp
c              ce(ipe,iqe) = ce(iqe,ipe)
c              ce(iqe,ipe) = ztemp
c           enddo
c       enddo

c     do ipe=1,nee
c        do iqe=1,nee
c           me(ipe,iqe) = dreal(mew(ipe,iqe))
c           ce(ipe,iqe) = ce(ipe,iqe) + dcmplx(0.d0,dimag(mew(ipe,iqe)))
c        enddo
c     enddo

      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CM25ANISO(MEN,MEE,MINTX, NEN,NINTX,NINTZ,
     ;                    OMEGA,PY, XIWTS,ETAWTS, VLAM,VMU,
     ;                    RHO,DET, SHG, GAMMAX,gammay,GAMMAZ, ME,CE)

      complex*16, intent(in) :: gammay(mintx,*)
      COMPLEX*16, INTENT(IN) :: GAMMAX(MINTX,*), GAMMAZ(MINTX,*)
      REAL*8, INTENT(IN) :: SHG(3,MEN,MINTX,*), VLAM(MINTX,*),
     ;        VMU(MINTX,*), RHO(MINTX,*), DET(MINTX,*), XIWTS(NINTX),
     ;        ETAWTS(NINTZ), OMEGA, PY
      INTEGER*4, INTENT(IN) :: MEN,MEE,MINTX, NEN,NINTX,NINTZ
      COMPLEX*16 CE(MEE,*)
      REAL*8 ME(MEE,*)
!.... local variables
      COMPLEX*16 D(9,9), BA(3,9), BB(9,3), S(9,3), K(3,3), 
     ;           GX, GY, GZ, PHIBX, PHIBZ, PHIAX, PHIAZ, WPY 
      REAL*8 RLAM, RMU, RHOW, WEIGHT, COEFFM1, COEFFM2, COEFFM3, 
     ;       RM11, RM22, RM33, PHIA, PHIB, PY2, W2
      COMPLEX*16 CFOUR, CTWO, CONE, CZERO
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0))
      PARAMETER(CONE  = DCMPLX(1.D0,0.D0)) 
      PARAMETER(CTWO  = DCMPLX(2.D0,0.D0))
      PARAMETER(CFOUR = DCMPLX(4.D0,0.D0)) 
      PARAMETER(NDIM = 3)
!
!----------------------------------------------------------------------!
!
      NEE = NDIM*NEN
      DO 1 IPE=1,NEE
         DO 2 IQE=1,NEE
            ME(IPE,IQE) = 0.D0
            CE(IPE,IQE) = CZERO
    2    CONTINUE
    1 CONTINUE 
      D(1:9,1:9) = CZERO
      BA(1:3,1:9) = CZERO
      BB(1:9,1:3) = CZERO 
      S(1:9,1:3) = CZERO
!
!.... loop on integration points
      WPY = DCMPLX(0.D0,OMEGA*PY)
      PY2  = PY**2
      W2 = OMEGA**2

      DO 3 INTX=1,NINTX
         DO 4 INTZ=1,NINTZ 
            GX = GAMMAX(INTX,INTZ)
            GY = GAMMAY(INTX,INTZ)
            GZ = GAMMAZ(INTX,INTZ) 
            WEIGHT = XIWTS(INTX)*ETAWTS(INTZ)*DET(INTX,INTZ)
            RLAM = VLAM(INTX,INTZ)*WEIGHT
            RMU  = VMU (INTX,INTZ)*WEIGHT
            RHOW = RHO(INTX,INTZ)*WEIGHT
 
            D(1,1) = DCMPLX( (RLAM + 2.D0*RMU),0.D0)*GY*GZ/GX
            D(1,2) = DCMPLX(RLAM,0.D0)*GZ
            D(1,3) = DCMPLX(RLAM,0.D0)*GY

            D(2,1) = D(1,2)
            D(2,2) = DCMPLX( (RLAM + 2.D0*RMU),0.D0)*GX*GZ/GY
            D(2,3) = DCMPLX(RLAM,0.D0)*GX 
           
            D(3,1) = D(1,3)
            D(3,2) = D(2,3)
            D(3,3) = DCMPLX( (RLAM + 2.D0*RMU),0.D0)*GX*GY/GZ

            D(4,4) = DCMPLX(RMU,0.D0) 
     ;              *( GZ/CTWO + GY*GZ/(CFOUR*GX) + GX*GZ/(CFOUR*GY))
            D(4,7) = DCMPLX(RMU/4.D0,0.D0)*(GY*GZ/GX - GX*GZ/GY)
            D(7,4) = D(4,7)

            D(5,5) = DCMPLX(RMU,0.D0) 
     ;              *( GY/CTWO + GY*GZ/(CFOUR*GX) + GX*GY/(CFOUR*GZ))
            D(5,8) = DCMPLX(RMU/4.D0,0.D0)*(GY*GZ/GX - GX*GY/GZ) 
            D(8,5) = D(5,8)

            D(6,6) = DCMPLX(RMU,0.D0)
     ;              *( GX/CTWO + GX*GZ/(CFOUR*GY) + GX*GY/(CFOUR*GZ))
            D(6,9) = DCMPLX(RMU/4.D0,0.D0)*(GY*GZ/GX - GX*GY/GZ)
            D(9,6) = D(6,9) 

            D(7,7) = DCMPLX(RMU,0.D0)
     ;              *(-GZ/CTWO + GY*GZ/(CFOUR*GX) + GX*GZ/(CFOUR*GY))
            D(8,8) = DCMPLX(RMU,0.D0)
     ;              *(-GY/CTWO + GY*GZ/(CFOUR*GX) + GX*GY/(CFOUR*GZ))
            D(9,9) = DCMPLX(RMU,0.D0)
     ;              *(-GX/CTWO + GX*GZ/(CFOUR*GY) + GX*GY/(CFOUR*GZ))

            COEFFM1 =-W2*RHOW
            COEFFM2 =-W2*RHOW
            COEFFM3 =-W2*RHOW 

            DO 5 IBE=1,NEN
               PHIB  = SHG(3,IBE,INTX,INTZ)
               PHIBX = DCMPLX(SHG(1,IBE,INTX,INTZ),0.D0)
               PHIBZ = DCMPLX(SHG(2,IBE,INTX,INTZ),0.D0)

               BB(1,1) = PHIBX
               BB(2,2) =-WPY*DCMPLX(PHIB,0.D0)
               BB(3,3) = PHIBZ
               BB(4,1) =-WPY*DCMPLX(PHIB,0.D0)
               BB(4,2) = PHIBX
               BB(5,1) = PHIBZ
               BB(5,3) = PHIBX
               BB(6,2) = PHIBZ 
               BB(6,3) =-WPY*DCMPLX(PHIB,0.D0)
               BB(7,1) =-WPY*DCMPLX(PHIB,0.D0)
               BB(7,2) =-PHIBX
               BB(8,1) = PHIBZ
               BB(8,3) =-PHIBX
               BB(9,2) = PHIBZ
               BB(9,3) = WPY*DCMPLX(PHIB,0.D0) 

               S = MATMUL(D,BB)

               RM11 = COEFFM1*PHIB
               RM22 = COEFFM2*PHIB
               RM33 = COEFFM3*PHIB 

               DO 6 IAE=1,NEN

                  PHIA  = SHG(3,IAE,INTX,INTZ)
                  PHIAX = DCMPLX(SHG(1,IAE,INTX,INTZ),0.D0)
                  PHIAZ = DCMPLX(SHG(2,IAE,INTX,INTZ),0.D0)

                  BA(1,1) = PHIAX
                  BA(2,2) =-WPY*DCMPLX(PHIA,0.D0)
                  BA(3,3) = PHIAZ
                  BA(1,4) =-WPY*DCMPLX(PHIA,0.D0)
                  BA(2,4) = PHIAX
                  BA(1,5) = PHIAZ
                  BA(3,5) = PHIAX
                  BA(2,6) = PHIAZ
                  BA(3,6) =-WPY*DCMPLX(PHIA,0.D0)
                  BA(1,7) =-WPY*DCMPLX(PHIA,0.D0)
                  BA(2,7) =-PHIAX
                  BA(1,8) = PHIAZ
                  BA(3,8) =-PHIAX
                  BA(2,9) = PHIAZ
                  BA(3,9) = WPY*DCMPLX(PHIA,0.D0)

                  K = MATMUL(BA,S)

                  IPE = NDIM*(IAE - 1) + 1
                  IQE = NDIM*(IBE - 1) + 1

                  ME(IPE  ,IQE  ) = ME(IPE  ,IQE  ) + RM11*PHIA 
                  ME(IPE+1,IQE+1) = ME(IPE+1,IQE+1) + RM22*PHIA
                  ME(IPE+2,IQE+2) = ME(IPE+2,IQE+2) + RM33*PHIA 

                  CE(IPE  ,IQE  ) = CE(IPE  ,IQE  ) + K(1,1)
                  CE(IPE  ,IQE+1) = CE(IPE  ,IQE+1) + K(1,2)
                  CE(IPE  ,IQE+2) = cE(IPE  ,IQE+2) + K(1,3) 

                  CE(IPE+1,IQE  ) = CE(IPE+1,IQE  ) + K(2,1)
                  CE(IPE+1,IQE+1) = CE(IPE+1,IQE+1) + K(2,2)
                  CE(IPE+1,IQE+2) = CE(IPE+1,IQE+2) + K(2,3) 

                  CE(IPE+2,IQE  ) = CE(IPE+2,IQE  ) + K(3,1)
                  CE(IPE+2,IQE+1) = CE(IPE+2,IQE+1) + K(3,2)
                  CE(IPE+2,IQE+2) = CE(IPE+2,IQE+2) + K(3,3) 
    6          CONTINUE
    5       CONTINUE
    4    CONTINUE
    3 CONTINUE 
      RETURN 
      END
