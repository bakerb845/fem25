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
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0))
      PARAMETER(NDIM = 3) 
! 
!----------------------------------------------------------------------!
! 
!.... null out element matrices
      NEE = NEN*NDIM
      DO 1 IPE=1,NEE
         DO 2 IQE=1,NEE !IPE,NEE
            ME(IPE,IQE) = 0.D0
            KE(IPE,IQE) = CZERO
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

               DO 6 IAE=1,NEN !IBE
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
!     DO 10 IPE=1,NEE
!        DO 11 IQE=IPE-1,1,-1
!           ME(IPE,IQE) = ME(IQE,IPE)
!           KE(IPE,IQE) = KE(IQE,IPE)
!  11    CONTINUE
!  10 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CM25ISOQ(MEN,MEE,MINTX, NEN,NINTX,NINTZ, 
     ;                    OMEGA,PY, XIWTS,ETAWTS, VLAM,VMU,
     ;                    RHO,DET, SHG, DAMPX,DAMPZ, ME,CE)
! 
!     Calculates all the element damping and mass matrices for 
!     the 2.5D elastic wave problem.  This is for isotropic 
!     elements. 
! 
!     INPUT      MEANING 
!     -----      ------- 
!     DET        holds the jacobian at the integration points
!     ETAWTS     holds integration weights in eta
!     MEE        max number of element equations
!     MEN        max number of element nodes
!     MINTX      leading dimension 
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
      COMPLEX*16 DAMPX(MINTX,*), DAMPZ(MINTX,*) 
      REAL*8 SHG(3,MEN,MINTX,*), VLAM(MINTX,*), VMU(MINTX,*), 
     ;       RHO(MINTX,*), DET(MINTX,*), XIWTS(NINTX), ETAWTS(NINTZ),
     ;       OMEGA, PY
      COMPLEX*16 CE(MEE,*)
      REAL*8 ME(MEE,*)
!.... local variables
      COMPLEX*16 PHIBX, PHIBZ, PHIAX, PHIAZ,
     ;           CE11, CE12, CE13, CE21, CE22, CE23, CE31, CE32, CE33,
     ;           CLAM, CMU, CZERO, CWPY, EX, EZ
      COMPLEX*16 C11,C12,C13, C22,C23,C33, C44,C55,C66
      COMPLEX*16 KE11T1,KE11T2,KE12T1,KE12T2,KE13T1,KE13T2,
     ;           KE21T1,KE21T2,KE22T1,KE22T2,KE23T1,KE23T2,
     ;           KE31T1,KE31T2,KE32T1,KE32T2,KE33T1,KE33T2
      REAL*8 PY2, W2, WEIGHT, LAMBDA, MU, RHOW, COEFFM1, COEFFM2, 
     ;       COEFFM3, PHIB, PHIA, WPY, PHIBR1,PHIBR2,PHIBR3
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0))
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
      !CWPY = DCMPLX(WPY,0.D0)
      PY2 = PY**2
      W2 = OMEGA**2
      DO 3 INTX=1,NINTX 
         DO 4 INTZ=1,NINTZ 
            WEIGHT = XIWTS(INTX)*ETAWTS(INTZ)*DET(INTX,INTZ)
!           EX = DCMPLX(DAMPX(INTX,INTZ,1),DAMPX(INTX,INTZ,2)) !1/ex
!           EZ = DCMPLX(DAMPZ(INTX,INTZ,1),DAMPZ(INTX,INTZ,2)) !1/ez
            EX = DAMPX(INTX,INTZ) !already 1/ex
            EZ = DAMPZ(INTX,INTZ) !already 1/ez
            LAMBDA = VLAM(INTX,INTZ)*WEIGHT 
            MU = VMU(INTX,INTZ)*WEIGHT 
            CLAM = DCMPLX(LAMBDA,0.D0) 
            CMU = DCMPLX(MU,0.D0)

            C11 = CLAM + CMU + CMU !DCMPLX(LAMBDA + 2.D0*MU,0.D0)
            C12 = CLAM !DCMPLX(LAMBDA,0.D0)
            C13 = CLAM !DCMPLX(LAMBDA,0.D0)

            C22 = CLAM + CMU + CMU !DCMPLX(LAMBDA + 2.D0*MU,0.D0)
            C23 = CLAM !DCMPLX(LAMBDA,0.D0)

            C33 = CLAM + CMU + CMU !DCMPLX(LAMBDA + 2.D0*MU,0.D0)

            C44 = CMU !DCMPLX(MU,0.D0)
            C55 = C44
            C66 = C44 !DCMPLX(MU,0.D0)

            RHOW = RHO(INTX,INTZ)*WEIGHT
            COEFFM1 =-W2*(RHOW - PY2*MU)
            COEFFM2 =-W2*(RHOW - PY2*(LAMBDA + 2.D0*MU))
            COEFFM3 =-W2*(RHOW - PY2*MU)
            DO 5 IBE=1,NEN
               PHIB  = SHG(3,IBE,INTX,INTZ)
               PHIBX = DCMPLX(SHG(1,IBE,INTX,INTZ),0.D0)*EX
               PHIBZ = DCMPLX(SHG(2,IBE,INTX,INTZ),0.D0)*EZ
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

               KE12T1 = CWPY*C12*DCMPLX(PHIB,0.D0)
               KE12T2 = CWPY*C66*PHIBX

               KE21T1 = CWPY*C12*PHIBX 
               KE21T2 = CWPY*C66*DCMPLX(PHIB,0.D0)
               KE23T1 = CWPY*C23*PHIBZ
               KE23T2 = CWPY*C44*DCMPLX(PHIB,0.D0)

               KE32T1 = CWPY*C23*DCMPLX(PHIB,0.D0)
               KE32T2 = CWPY*C44*PHIBZ

               DO 6 IAE=1,NEN 
                  PHIA  = SHG(3,IAE,INTX,INTZ)
                  PHIAX = DCMPLX(SHG(1,IAE,INTX,INTZ),0.D0)*EX
                  PHIAZ = DCMPLX(SHG(2,IAE,INTX,INTZ),0.D0)*EZ
                  !indices
                  IPE = NDIM*(IAE - 1) + 1
                  IQE = NDIM*(IBE - 1) + 1
                  !mass matrix
                  ME(IPE  ,IQE  ) = ME(IPE  ,IQE  ) + PHIA*PHIBR1
                  ME(IPE+1,IQE+1) = ME(IPE+1,IQE+1) + PHIA*PHIBR2
                  ME(IPE+2,IQE+2) = ME(IPE+2,IQE+2) + PHIA*PHIBR3

                  CE11 = PHIAX*KE11T1 + PHIAZ*KE11T2
                  CE13 = PHIAX*KE13T1 + PHIAZ*KE13T2
                  CE22 = PHIAX*KE22T1 + PHIAZ*KE22T2
                  CE31 = PHIAZ*KE31T1 + PHIAX*KE31T2
                  CE33 = PHIAZ*KE33T1 + PHIAX*KE33T2
                  !stiffness, imaginary 
                  CE12 = PHIAX*KE12T1 - DCMPLX(PHIA,0.D0)*KE12T2 
                  CE21 = DCMPLX(PHIA,0.D0)*KE21T1 - PHIAX*KE21T2
                  CE23 = DCMPLX(PHIA,0.D0)*KE23T1 - PHIAZ*KE23T2 
                  CE32 = PHIAZ*KE32T1 - DCMPLX(PHIA,0.D0)*KE32T2

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
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE KM_VPVS(MEN,MEE,MINTX, NEN,NINTX,NINTZ,
     ;                   FREQ,PY, XIWTS,ETAWTS, VPSQR,VSSQR,
     ;                   DET,SHG, ME,KE)
!
!     2.5D isotropic elastic element assembly with P and S wave 
!     velocities.  Useful for inversion  - B. Baker May 2013
!
!     INPUT      MEANING
!     -----      ------- 
!     DET        jacobian at integration points
!     ETAWTS     eta integration weights
!     FREQ       current frequency (Hz)
!     FREQ0      reference frequency for damping 
!     MEE        leading dimension
!     MEN        leading dimension
!     MINTX      leading dimension
!     NEN        number of element nodes
!     NINTX      number of integration points in xi
!     NINTZ      number of integration points in eta 
!     PY         apparent slowness in y (s/m)
!     SHG        global shape fns and derivatives at int. pts
!     XIWTS      xi integration weights
!     VPSQR      Vp^2 velocities at integration points (m/s) 
!     VSSQR      Vs^2  velocities at integration points (m/s)
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     KE         damped element stiffness with complex element mass term
!     ME         real element mass
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: SHG(3,MEN,MINTX,*),VPSQR(MINTX,*), 
     ;                      VSSQR(MINTX,*),
     ;                      DET(MINTX,*), XIWTS(NINTX),ETAWTS(NINTZ), 
     ;                      FREQ,PY
      INTEGER*4, INTENT(IN) :: MEN, MEE, MINTX, NEN,NINTX,NINTZ
      COMPLEX*16, INTENT(OUT) :: KE(MEE,*)
      REAL*8, INTENT(OUT) :: ME(MEE,*)
!.... local variables
      COMPLEX*16 ZZERO 
      REAL*8 TWOPI, OMEGA, WPY, PY2, W2, WEIGHT, RVP2,RVS2,
     ;       RMCOEFF1,RMCOEFF2,RMCOEFF3,
     ;       RK11T1,RK11T2, RK13T1,RK13T2,RK13T3,
     ;       RK22T1,RK22T2, RK31T1,RK31T2,RK31T3,
     ;       RK33T1,RK33T2, RK12T1,RK12T2,RK12T3,
     ;       RK21T1,RK21T2,RK21T3, RK23T1,RK23T2,RK23T3,
     ;       RK32T1,RK32T2,RK32T3, RM11T1,RM22T1,RM33T1,
     ;       RKE11, RKE12, RKE13, RKE21, RKE22, RKE23,
     ;       RKE31, RKE32, RKE33, RME11, RME22, RME33,
     ;       PHIA,PHIB, PHIAX,PHIBX,PHIAZ,PHIBZ
      INTEGER*4 NEE,NDIM, IPE,IQE, INTX,INTZ, IAE,IBE
      PARAMETER(ZZERO = DCMPLX(0.D0,0.D0))
      PARAMETER(TWOPI = 6.2831853071795862D0)
      PARAMETER(NDIM = 3)
! 
!----------------------------------------------------------------------!
! 
!.... null out element matrices
      NEE = NEN*NDIM
      DO 1 IPE=1,NEE
         DO 2 IQE=1,NEE
            ME(IPE,IQE) = 0.D0
            KE(IPE,IQE) = ZZERO
    2    CONTINUE
    1 CONTINUE
! 
!.... loop on integration points
      OMEGA = TWOPI*FREQ
      WPY = OMEGA*PY
      PY2 = PY**2
      W2 = OMEGA**2
      DO 3 INTX=1,NINTX
         DO 4 INTZ=1,NINTZ
            WEIGHT = XIWTS(INTX)*ETAWTS(INTZ)*DET(INTX,INTZ)
            RVP2 = VPSQR(INTX,INTZ)*WEIGHT
            RVS2 = VSSQR(INTX,INTZ)*WEIGHT
 
            RMCOEFF1 =-W2*WEIGHT + W2*PY2*RVS2
            RMCOEFF2 =-W2*WEIGHT + W2*PY2*RVP2
            RMCOEFF3 = RMCOEFF1
            DO 5 IBE=1,NEN
               PHIB  = SHG(3,IBE,INTX,INTZ)
               PHIBX = SHG(1,IBE,INTX,INTZ)
               PHIBZ = SHG(2,IBE,INTX,INTZ)

               RM11T1 = RMCOEFF1*PHIB
               RM22T1 = RMCOEFF2*PHIB
               RM33T1 = RM11T1

               RK11T1 = RVP2*PHIBX
               RK11T2 = RVS2*PHIBZ

               RK13T1 = RVP2*PHIBZ
               RK13T2 = RVS2*PHIBX
               RK13T3 =-2.D0*RVS2*PHIBZ

               RK22T1 = RVS2*PHIBX
               RK22T2 = RVS2*PHIBZ

               RK31T1 = RVP2*PHIBX
               RK31T2 = RVS2*PHIBZ
               RK31T3 =-2.D0*RVS2*PHIBX

               RK33T1 = RVP2*PHIBZ
               RK33T2 = RVS2*PHIBX

               RK12T1 = WPY*RVP2*PHIB
               RK12T2 =-2.D0*WPY*RVS2*PHIB
               RK12T3 =     -WPY*RVS2*PHIBX

               RK21T1 = WPY*RVP2*PHIBX
               RK21T2 =-2.D0*WPY*RVS2*PHIBX
               RK21T3 =     -WPY*RVS2*PHIB

               RK23T1 = WPY*RVP2*PHIBZ
               RK23T2 =-2.D0*WPY*RVS2*PHIBZ
               RK23T3 =     -WPY*RVS2*PHIB

               RK32T1 = WPY*RVP2*PHIB
               RK32T2 =-2.D0*WPY*RVS2*PHIB
               RK32T3 =     -WPY*RVS2*PHIBZ

               DO 6 IAE=1,NEN
                  PHIA  = SHG(3,IAE,INTX,INTZ)
                  PHIAX = SHG(1,IAE,INTX,INTZ)
                  PHIAZ = SHG(2,IAE,INTX,INTZ)

                  RME11 = RM11T1*PHIA
                  RME22 = RM22T1*PHIA
                  RME33 = RME11

                  RKE11 = RK11T1*PHIAX + RK11T2*PHIAZ
                  RKE13 = RK13T1*PHIAX + RK13T2*PHIAZ + RK13T3*PHIAX
                  RKE22 = RK22T1*PHIAX + RK22T2*PHIAZ
                  RKE31 = RK31T1*PHIAZ + RK31T2*PHIAX + RK31T3*PHIAZ
                  RKE33 = RK33T1*PHIAZ + RK33T2*PHIAX

                  RKE12 = RK12T1*PHIAX + RK12T2*PHIAX + RK12T3*PHIA
                  RKE21 = RK21T1*PHIA  + RK21T2*PHIA  + RK21T3*PHIAX

                  RKE23 = RK23T1*PHIA  + RK23T2*PHIA  + RK23T3*PHIAZ
                  RKE32 = RK32T1*PHIAZ + RK32T2*PHIAZ + RK32T3*PHIA

                  !indices
                  IPE = NDIM*(IAE - 1) + 1
                  IQE = NDIM*(IBE - 1) + 1
                  !mass 
                  ME(IPE  ,IQE  ) = ME(IPE  ,IQE  ) + RME11
                  ME(IPE+1,IQE+1) = ME(IPE+1,IQE+1) + RME22
                  ME(IPE+2,IQE+2) = ME(IPE+2,IQE+2) + RME33
                  !main diagonal stiffenss
                  KE(IPE  ,IQE  ) = KE(IPE  ,IQE  ) + DCMPLX(RKE11,0.D0)
                  KE(IPE+1,IQE+1) = KE(IPE+1,IQE+1) + DCMPLX(RKE22,0.D0)
                  KE(IPE+2,IQE+2) = KE(IPE+2,IQE+2) + DCMPLX(RKE33,0.D0)
                  !real off diagonal 
                  KE(IPE  ,IQE+2) = KE(IPE  ,IQE+2) + DCMPLX(RKE13,0.D0)
                  KE(IPE+2,IQE  ) = KE(IPE+2,IQE  ) + DCMPLX(RKE31,0.D0)
                  !complex off diagonal
                  KE(IPE  ,IQE+1) = KE(IPE  ,IQE+1) - DCMPLX(RKE12,0.D0)
                  KE(IPE+1,IQE  ) = KE(IPE+1,IQE  ) + DCMPLX(RKE21,0.D0)
                  KE(IPE+1,IQE+2) = KE(IPE+1,IQE+2) + DCMPLX(RKE23,0.D0)
                  KE(IPE+2,IQE+1) = KE(IPE+2,IQE+1) - DCMPLX(RKE32,0.D0)
    6          CONTINUE !loop on iae
    5       CONTINUE !Loop on ibe
    4    CONTINUE !loop on eta integration
    3 CONTINUE !loop on xi integration
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE KMDAMP(MEN,MEE,MINTX, NEN,NINTX,NINTZ,
     ;                  FREQ,FREQ0,PY, XIWTS,ETAWTS, VP,VS,
     ;                  QPP,QSS, DET,SHG, ME,KE)
!
!     2.5D isotropic elastic element assembly with viscous damping 
!     - B. Baker April 2013
!
!     INPUT      MEANING
!     -----      ------- 
!     DET        jacobian at integration points
!     ETAWTS     eta integration weights
!     FREQ       current frequency (Hz)
!     FREQ0      reference frequency for damping 
!     MEE        leading dimension
!     MEN        leading dimension
!     MINTX      leading dimension
!     NEN        number of element nodes
!     NINTX      number of integration points in xi
!     NINTZ      number of integration points in eta 
!     PY         apparent slowness in y (s/m)
!     QPP        Vp quality factor at integration points
!     QSS        Vs quality factor at integration points 
!     SHG        global shape fns and derivatives at int. pts
!     XIWTS      xi integration weights
!     VP         Vp velocities at integration points (m/s)
!     VS         Vs velocities at integration points (m/s)
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     KE         damped element stiffness with complex element mass term
!     ME         real element mass
!     
!.... variable declarations
      REAL*8, INTENT(IN) :: SHG(3,MEN,MINTX,*),VP(MINTX,*),VS(MINTX,*),
     ;                      DET(MINTX,*),QPP(MINTX,*),QSS(MINTX,*), 
     ;                      XIWTS(NINTX),ETAWTS(NINTZ), FREQ,FREQ0,PY 
      INTEGER*4, INTENT(IN) :: MEN, MEE, MINTX, NEN,NINTX,NINTZ
      COMPLEX*16, INTENT(OUT) :: KE(MEE,*) 
      REAL*8, INTENT(OUT) :: ME(MEE,*) 
!.... local variables
      COMPLEX*16 VCMPLX, ZVP2, ZVS2, ZZERO
      REAL*8 TWOPI, OMEGA, WPY, PY2, W2, WEIGHT, RVP2,QVP2,RVS2,QVS2, 
     ;       RMCOEFF1,QMCOEFF1, RMCOEFF2,QMCOEFF2, RMCOEFF3,QMCOEFF3, 
     ;       RK11T1,QK11T1,RK11T2,QK11T2, 
     ;       RK13T1,QK13T1,RK13T2,QK13T2,RK13T3,QK13T3, 
     ;       RK22T1,QK22T1,RK22T2,QK22T2, 
     ;       RK31T1,QK31T1,RK31T2,QK31T2,RK31T3,QK31T3,
     ;       RK33T1,QK33T1,RK33T2,QK33T2,   
     ;       RK12T1,QK12T1,RK12T2,QK12T2,RK12T3,QK12T3, 
     ;       RK21T1,QK21T1,RK21T2,QK21T2,RK21T3,QK21T3, 
     ;       RK23T1,QK23T1,RK23T2,QK23T2,RK23T3,QK23T3, 
     ;       RK32T1,QK32T1,RK32T2,QK32T2,RK32T3,QK32T3,
     ;       RM11T1,QM11T1,RM22T1,QM22T1,RM33T1,QM33T1, 
     ;       RKE11,QKE11, RKE12,QKE12, RKE13,QKE13, 
     ;       RKE21,QKE21, RKE22,QKE22, RKE23,QKE23, 
     ;       RKE31,QKE31, RKE32,QKE32, RKE33,QKE33,  
     ;       RME11,QME11, RME22,QME22, RME33,QME33, 
     ;       PHIA,PHIB, PHIAX,PHIBX,PHIAZ,PHIBZ   
      INTEGER*4 NEE,NDIM, IPE,IQE, INTX,INTZ, IAE,IBE
      PARAMETER(ZZERO = DCMPLX(0.D0,0.D0)) 
      PARAMETER(TWOPI = 6.2831853071795862D0)
      PARAMETER(NDIM = 3)
! 
!----------------------------------------------------------------------!
! 
!.... null out element matrices
      NEE = NEN*NDIM
      DO 1 IPE=1,NEE
         DO 2 IQE=1,NEE
            ME(IPE,IQE) = 0.D0
            KE(IPE,IQE) = ZZERO
    2    CONTINUE
    1 CONTINUE
! 
!.... loop on integration points
      OMEGA = TWOPI*FREQ
      WPY = OMEGA*PY
      PY2 = PY**2
      W2 = OMEGA**2
      DO 3 INTX=1,NINTX
         DO 4 INTZ=1,NINTZ
            WEIGHT = XIWTS(INTX)*ETAWTS(INTZ)*DET(INTX,INTZ)
            ZVP2 = VCMPLX(FREQ,FREQ0,VP(INTX,INTZ),QPP(INTX,INTZ))
            ZVS2 = VCMPLX(FREQ,FREQ0,VS(INTX,INTZ),QSS(INTX,INTZ))
            ZVP2 = ZVP2**2
            ZVS2 = ZVS2**2
            
            RVP2 = DREAL(ZVP2)*WEIGHT
            QVP2 = DIMAG(ZVP2)*WEIGHT
            RVS2 = DREAL(ZVS2)*WEIGHT
            QVS2 = DIMAG(ZVS2)*WEIGHT
            
            RMCOEFF1 =-W2*WEIGHT + W2*PY2*RVS2
            QMCOEFF1 =             W2*PY2*QVS2
            RMCOEFF2 =-W2*WEIGHT + W2*PY2*RVP2
            QMCOEFF2 =             W2*PY2*QVP2
            RMCOEFF3 = RMCOEFF1
            QMCOEFF3 = QMCOEFF1

            DO 5 IBE=1,NEN

               PHIB  = SHG(3,IBE,INTX,INTZ)
               PHIBX = SHG(1,IBE,INTX,INTZ)
               PHIBZ = SHG(2,IBE,INTX,INTZ)

               RM11T1 = RMCOEFF1*PHIB
               QM11T1 = QMCOEFF1*PHIB 
               RM22T1 = RMCOEFF2*PHIB
               QM22T1 = QMCOEFF2*PHIB
               RM33T1 = RM11T1
               QM33T1 = QM11T1

               RK11T1 = RVP2*PHIBX
               QK11T1 = QVP2*PHIBX
               RK11T2 = RVS2*PHIBZ
               QK11T2 = QVS2*PHIBZ

               RK13T1 = RVP2*PHIBZ
               QK13T1 = QVP2*PHIBZ
               RK13T2 = RVS2*PHIBX
               QK13T2 = QVS2*PHIBX
               RK13T3 =-2.D0*RVS2*PHIBZ
               QK13T3 =-2.D0*QVS2*PHIBZ

               RK22T1 = RVS2*PHIBX
               QK22T1 = QVS2*PHIBX
               RK22T2 = RVS2*PHIBZ
               QK22T2 = QVS2*PHIBZ

               RK31T1 = RVP2*PHIBX
               QK31T1 = QVP2*PHIBX
               RK31T2 = RVS2*PHIBZ
               QK31T2 = QVS2*PHIBZ
               RK31T3 =-2.D0*RVS2*PHIBX
               QK31T3 =-2.D0*QVS2*PHIBX

               RK33T1 = RVP2*PHIBZ
               QK33T1 = QVP2*PHIBZ
               RK33T2 = RVS2*PHIBX
               QK33T2 = QVS2*PHIBX

               RK12T1 = WPY*RVP2*PHIB
               QK12T1 = WPY*QVP2*PHIB
               RK12T2 =-2.D0*WPY*RVS2*PHIB
               QK12T2 =-2.D0*WPY*QVS2*PHIB
               RK12T3 =     -WPY*RVS2*PHIBX
               QK12T3 =     -WPY*QVS2*PHIBX

               RK21T1 = WPY*RVP2*PHIBX
               QK21T1 = WPY*QVP2*PHIBX
               RK21T2 =-2.D0*WPY*RVS2*PHIBX 
               QK21T2 =-2.D0*WPY*QVS2*PHIBX
               RK21T3 =     -WPY*RVS2*PHIB
               QK21T3 =     -WPY*QVS2*PHIB

               RK23T1 = WPY*RVP2*PHIBZ
               QK23T1 = WPY*QVP2*PHIBZ
               RK23T2 =-2.D0*WPY*RVS2*PHIBZ
               QK23T2 =-2.D0*WPY*QVS2*PHIBZ
               RK23T3 =     -WPY*RVS2*PHIB
               QK23T3 =     -WPY*QVS2*PHIB

               RK32T1 = WPY*RVP2*PHIB 
               QK32T1 = WPY*QVP2*PHIB
               RK32T2 =-2.D0*WPY*RVS2*PHIB 
               QK32T2 =-2.D0*WPY*QVS2*PHIB
               RK32T3 =     -WPY*RVS2*PHIBZ
               QK32T3 =     -WPY*QVS2*PHIBZ

               DO 6 IAE=1,NEN
                  PHIA  = SHG(3,IAE,INTX,INTZ)
                  PHIAX = SHG(1,IAE,INTX,INTZ)
                  PHIAZ = SHG(2,IAE,INTX,INTZ)

                  RME11 = RM11T1*PHIA
                  QME11 = QM11T1*PHIA
                  RME22 = RM22T1*PHIA
                  QME22 = QM22T1*PHIA 
                  RME33 = RME11
                  QME33 = QME11 

                  RKE11 = RK11T1*PHIAX + RK11T2*PHIAZ
                  QKE11 = QK11T1*PHIAX + QK11T2*PHIAZ

                  RKE13 = RK13T1*PHIAX + RK13T2*PHIAZ + RK13T3*PHIAX
                  QKE13 = QK13T1*PHIAX + QK13T2*PHIAZ + QK13T3*PHIAX

                  RKE22 = RK22T1*PHIAX + RK22T2*PHIAZ
                  QKE22 = QK22T1*PHIAX + QK22T2*PHIAZ

                  RKE31 = RK31T1*PHIAZ + RK31T2*PHIAX + RK31T3*PHIAZ
                  QKE31 = QK31T1*PHIAZ + QK31T2*PHIAX + QK31T3*PHIAZ

                  RKE33 = RK33T1*PHIAZ + RK33T2*PHIAX 
                  QKE33 = QK33T1*PHIAZ + QK33T2*PHIAX 

                  RKE12 = RK12T1*PHIAX + RK12T2*PHIAX + RK12T3*PHIA
                  QKE12 = QK12T1*PHIAX + QK12T2*PHIAX + QK12T3*PHIA
                  RKE21 = RK21T1*PHIA  + RK21T2*PHIA  + RK21T3*PHIAX
                  QKE21 = QK21T1*PHIA  + QK21T2*PHIA  + QK21T3*PHIAX 
  
                  RKE23 = RK23T1*PHIA  + RK23T2*PHIA  + RK23T3*PHIAZ
                  QKE23 = QK23T1*PHIA  + QK23T2*PHIA  + QK23T3*PHIAZ
                  RKE32 = RK32T1*PHIAZ + RK32T2*PHIAZ + RK32T3*PHIA
                  QKE32 = QK32T1*PHIAZ + QK32T2*PHIAZ + QK32T3*PHIA 

                  !indices
                  IPE = NDIM*(IAE - 1) + 1 
                  IQE = NDIM*(IBE - 1) + 1 
                  !mass 
                  ME(IPE  ,IQE  ) = ME(IPE  ,IQE  ) + RME11
                  ME(IPE+1,IQE+1) = ME(IPE+1,IQE+1) + RME22
                  ME(IPE+2,IQE+2) = ME(IPE+2,IQE+2) + RME33 
                  !add the complex part of the mass matrix in 
                  KE(IPE  ,IQE  ) = KE(IPE  ,IQE  )
     ;                            + DCMPLX(RKE11,QKE11+QME11)
                  KE(IPE+1,IQE+1) = KE(IPE+1,IQE+1) 
     ;                            + DCMPLX(RKE22,QKE22+QME22)
                  KE(IPE+2,IQE+2) = KE(IPE+2,IQE+2) 
     ;                            + DCMPLX(RKE33,QKE33+QME33)
                  !off diagonal
                  KE(IPE  ,IQE+2) = KE(IPE  ,IQE+2) 
     ;                            + DCMPLX(RKE13,QKE13)
                  KE(IPE+2,IQE  ) = KE(IPE+2,IQE  ) 
     ;                            + DCMPLX(RKE31,QKE31)
                  !multiply by -i hence the switch
                  KE(IPE  ,IQE+1) = KE(IPE  ,IQE+1) 
     ;                            - DCMPLX(-QKE12,RKE12)
                  KE(IPE+1,IQE  ) = KE(IPE+1,IQE  ) 
     ;                            + DCMPLX(-QKE21,RKE21)
                  KE(IPE+1,IQE+2) = KE(IPE+1,IQE+2) 
     ;                            + DCMPLX(-QKE23,RKE23)
                  KE(IPE+2,IQE+1) = KE(IPE+2,IQE+1) 
     ;                            - DCMPLX(-QKE32,RKE32)

                   
    6          CONTINUE
    5       CONTINUE  
    4    CONTINUE !Loop on integration points in z 
    3 CONTINUE !Loop on points in x 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CMDAMP(MEN,MEE,MINTX, NEN,NINTX,NINTZ,
     ;                  FREQ,FREQ0,PY, XIWTS,ETAWTS, VP,VS,
     ;                  QPP,QSS, DET,SHG, DAMPX,DAMPZ, ME,CE)
!
!     INPUT      MEANING
!     -----      ------- 
!     DAMPX      PML stretching factor in x at int. pts
!     DAMPZ      PML stretching factor in z at int. pts 
!     DET        jacobian at integration points
!     ETAWTS     eta integration weights
!     FREQ       current frequency (Hz)
!     FREQ0      reference frequency for damping 
!     MEE        leading dimension
!     MEN        leading dimension
!     MINTX      leading dimension
!     NEN        number of element nodes
!     NINTX      number of integration points in xi
!     NINTZ      number of integration points in eta 
!     PY         apparent slowness in y (s/m)
!     QPP        Vp quality factor at integration points
!     QSS        Vs quality factor at integration points 
!     SHG        global shape fns and derivatives at int. pts
!     XIWTS      xi integration weights
!     VP         Vp velocities at integration points (m/s)
!     VS         Vs velocities at integration points (m/s)
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     KE         damped element stiffness with complex element mass term
!     ME         real element mass
!
!.... variable declarations
      COMPLEX*16, INTENT(IN) :: DAMPX(MINTX,*), DAMPZ(MINTX,*)  
      REAL*8, INTENT(IN) :: SHG(3,MEN,MINTX,*),VP(MINTX,*),VS(MINTX,*),
     ;                      DET(MINTX,*),QPP(MINTX,*),QSS(MINTX,*),
     ;                      XIWTS(NINTX),ETAWTS(NINTZ), FREQ,FREQ0,PY
      INTEGER*4, INTENT(IN) :: MEN, MEE, MINTX, NEN,NINTX,NINTZ
      COMPLEX*16, INTENT(OUT) :: CE(MEE,*)
      REAL*8, INTENT(OUT) :: ME(MEE,*)
!.... local variables
      COMPLEX*16 ZTWO, ZZERO, VCMPLX, EX, EZ, ZVP2, ZVS2, 
     ;           CE11T1,CE11T2, CE22T1,CE22T2, CE33T1,CE33T2, 
     ;           CE13T1,CE13T2,CE13T3, CE31T1,CE31T2,CE31T3, 
     ;           CE12T1,CE12T2,CE12T3, CE21T1,CE21T2,CE21T3, 
     ;           CE23T1,CE23T2,CE23T3, CE32T1,CE32T2,CE32T3, 
     ;           CE11,CE12,CE13, CE21,CE22,CE23, CE31,CE32,CE33, 
     ;           PHIAX,PHIBX, PHIAZ,PHIBZ, CWPY   
      REAL*8 OMEGA, TWOPI, WPY, PY2, W2, RVP2, QVP2, RVS2, QVS2, WEIGHT,
     ;       PHIA, PHIB, RMCOEFF1, RMCOEFF2, RMCOEFF3, 
     ;       QMCOEFF1, QMCOEFF2, QMCOEFF3, 
     ;       RM11T1, QM11T1, RM22T1, QM22T1, RM33T1, QM33T1, 
     ;       RME11,QME11, RME22,QME22, RME33,QME33   
      INTEGER*4 NEE, NDIM, INTX, INTZ, IAE, IBE, IPE, IQE
      PARAMETER(ZTWO = DCMPLX(2.D0,0.D0), ZZERO = DCMPLX(0.D0,0.D0))  
      PARAMETER(TWOPI = 6.2831853071795862D0)
      PARAMETER(NDIM = 3) 
! 
!----------------------------------------------------------------------!
! 
!.... null out element matrices
      NEE = NEN*NDIM
      DO 1 IPE=1,NEE
         DO 2 IQE=1,NEE
            ME(IPE,IQE) = 0.D0
            CE(IPE,IQE) = ZZERO
    2    CONTINUE
    1 CONTINUE
! 
!.... loop on integration points
      OMEGA = TWOPI*FREQ
      WPY = OMEGA*PY
      CWPY = DCMPLX(0.D0,WPY)
      !CWPY = DCMPLX(WPY,0.D0)
      PY2 = PY**2
      W2 = OMEGA**2
      DO 3 INTX=1,NINTX
         DO 4 INTZ=1,NINTZ
            WEIGHT = XIWTS(INTX)*ETAWTS(INTZ)*DET(INTX,INTZ)
            EX = DAMPX(INTX,INTZ) !already 1/ex
            EZ = DAMPZ(INTX,INTZ) !already 1/ez
            ZVP2 = VCMPLX(FREQ,FREQ0,VP(INTX,INTZ),QPP(INTX,INTZ))
            ZVS2 = VCMPLX(FREQ,FREQ0,VS(INTX,INTZ),QSS(INTX,INTZ))
            ZVP2 = ZVP2**2*DCMPLX(WEIGHT,0.D0)
            ZVS2 = ZVS2**2*DCMPLX(WEIGHT,0.D0)

            RVP2 = DREAL(ZVP2)
            QVP2 = DIMAG(ZVP2)
            RVS2 = DREAL(ZVS2)
            QVS2 = DIMAG(ZVS2)

            RMCOEFF1 =-W2*WEIGHT + W2*PY2*RVS2
            QMCOEFF1 =             W2*PY2*QVS2
            RMCOEFF2 =-W2*WEIGHT + W2*PY2*RVP2
            QMCOEFF2 =             W2*PY2*QVP2
            RMCOEFF3 = RMCOEFF1
            QMCOEFF3 = QMCOEFF1

            DO 5 IBE=1,NEN

               PHIB  = SHG(3,IBE,INTX,INTZ)
               PHIBX = SHG(1,IBE,INTX,INTZ)*EX
               PHIBZ = SHG(2,IBE,INTX,INTZ)*EZ

               RM11T1 = RMCOEFF1*PHIB
               QM11T1 = QMCOEFF1*PHIB
               RM22T1 = RMCOEFF2*PHIB
               QM22T1 = QMCOEFF2*PHIB
               RM33T1 = RM11T1
               QM33T1 = QM11T1

               CE11T1 = ZVP2*PHIBX
               CE11T2 = ZVS2*PHIBZ

               CE13T1 = ZVP2*PHIBZ
               CE13T2 = ZVS2*PHIBX
               CE13T3 =-ZTWO*ZVS2*PHIBZ

               CE22T1 = ZVS2*PHIBX 
               CE22T2 = ZVS2*PHIBZ

               CE31T1 = ZVP2*PHIBX
               CE31T2 = ZVS2*PHIBZ
               CE31T3 =-ZTWO*ZVS2*PHIBX

               CE33T1 = ZVP2*PHIBZ 
               CE33T2 = ZVS2*PHIBX

               CE12T1 =      CWPY*ZVP2*DCMPLX(PHIB,0.D0) 
               CE12T2 =-ZTWO*CWPY*ZVS2*DCMPLX(PHIB,0.D0)
               CE12T3 =     -CWPY*ZVS2*PHIBX
 
               CE21T1 =      CWPY*ZVP2*PHIBX
               CE21T2 =-ZTWO*CWPY*ZVS2*PHIBX 
               CE21T3 =     -CWPY*ZVS2*DCMPLX(PHIB,0.D0)

               CE23T1 =      CWPY*ZVP2*PHIBZ
               CE23T2 =-ZTWO*CWPY*ZVS2*PHIBZ
               CE23T3 =     -CWPY*ZVS2*DCMPLX(PHIB,0.D0) 

               CE32T1 =      CWPY*ZVP2*DCMPLX(PHIB,0.D0) 
               CE32T2 =-ZTWO*CWPY*ZVS2*DCMPLX(PHIB,0.D0) 
               CE32T3 =     -CWPY*ZVS2*PHIBZ

               DO 6 IAE=1,NEN
                  PHIA  = SHG(3,IAE,INTX,INTZ)
                  PHIAX = SHG(1,IAE,INTX,INTZ)*EX
                  PHIAZ = SHG(2,IAE,INTX,INTZ)*EZ

                  RME11 = RM11T1*PHIA
                  QME11 = QM11T1*PHIA
                  RME22 = RM22T1*PHIA
                  QME22 = QM22T1*PHIA 
                  RME33 = RME11
                  QME33 = QME11 

                  CE11 = CE11T1*PHIAX + CE11T2*PHIAZ
                  CE13 = CE13T1*PHIAX + CE13T2*PHIAZ + CE13T3*PHIAX
                  CE22 = CE22T1*PHIAX + CE22T2*PHIAZ 
                  CE31 = CE31T1*PHIAZ + CE31T2*PHIAX + CE31T3*PHIAZ
                  CE33 = CE33T1*PHIAZ + CE33T2*PHIAX 

                  CE12 = CE12T1*PHIAX 
     ;                 + CE12T2*PHIAX
     ;                 + CE12T3*DCMPLX(PHIA,0.D0)

                  CE21 = CE21T1*DCMPLX(PHIA,0.D0)
     ;                 + CE21T2*DCMPLX(PHIA,0.D0)
     ;                 + CE21T3*PHIAX

                  CE23 = CE23T1*DCMPLX(PHIA,0.D0)
     ;                 + CE23T2*DCMPLX(PHIA,0.D0)
     ;                 + CE23T3*PHIAZ

                  CE32 = CE32T1*PHIAZ
     ;                 + CE32T2*PHIAZ
     ;                 + CE32T3*DCMPLX(PHIA,0.D0)

                  !indices
                  IPE = NDIM*(IAE - 1) + 1 
                  IQE = NDIM*(IBE - 1) + 1 
                  !mass
                  ME(IPE  ,IQE  ) = ME(IPE  ,IQE  ) + RME11
                  ME(IPE+1,IQE+1) = ME(IPE+1,IQE+1) + RME22
                  ME(IPE+2,IQE+2) = ME(IPE+2,IQE+2) + RME33
                  !stiffness 
                  CE(IPE  ,IQE  ) = CE(IPE  ,IQE  ) + CE11
     ;                            + DCMPLX(0.D0,QME11)
                  CE(IPE  ,IQE+1) = CE(IPE  ,IQE+1) - CE12
                  CE(IPE  ,IQE+2) = CE(IPE  ,IQE+2) + CE13

                  CE(IPE+1,IQE  ) = CE(IPE+1,IQE  ) + CE21
                  CE(IPE+1,IQE+1) = CE(IPE+1,IQE+1) + CE22
     ;                            + DCMPLX(0.D0,QME22)
                  CE(IPE+1,IQE+2) = CE(IPE+1,IQE+2) + CE23

                  CE(IPE+2,IQE  ) = CE(IPE+2,IQE  ) + CE31
                  CE(IPE+2,IQE+1) = CE(IPE+2,IQE+1) - CE32
                  CE(IPE+2,IQE+2) = CE(IPE+2,IQE+2) + CE33
     ;                            + DCMPLX(0.D0,QME33) 

    6          CONTINUE !loop on element nodes, iae
    5       CONTINUE !loop on element nodes, ibe 
    4   CONTINUE !loop on eta integration points 
    3 CONTINUE !loop on xi integration points
      RETURN
      END
