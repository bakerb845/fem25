      SUBROUTINE DABMKISOB25(MEN,MEE,MINTX,MINTZ, NEN,NINTX,NINTZ, 
     ;                       LISP,LISS, OMEGA,PY, 
     ;                       XIWTS,ETAWTS, 
     ;                       DET,DAB,SHG, DMDA,DMDB,DKDA,DKDB)  
!
!     Derivative of the stiffness matrix with respect to alpha and 
!     beta 
!
!     INPUT      MEANING
!     -----      ------- 
!     DAB        derivatives w.r.t. to alpha and beta at int. pts
!     DET        jacobian 
!     ETAWTS     eta integration weights
!     LISP       True -> calculate P derivative element matrices
!     LISS       True -. calculate S derivative element matrices
!     MEE        leading dimension for element matrices 
!     MEN        leading dimension for SHG 
!     MINTX      leading dimension DET, DAB, SHG  
!     MINTZ      leading dimension for DAB
!     NEN        number of element nodes
!     NINTX      number of integration points in xi 
!     NINTZ      number of integration points in eta
!     OMEGA      angular frequency 
!     PY         apparaent slowness in y 
!     SHG        global shape fns and derivatives at int. pts
!     XIWTS      xi integration weights
!  
!     OUTPUT     MEANING
!     ------     ------- 
!     DKDA       derivative of the stiffness matrix w.r.t. alpha
!     DKDB       derivative of the stiffness matrix w.r.t. beta
!     DMDA       derivative of the mass matrix w.r.t. alpha
!     DMDB       derivative of the mass matrix w.r.t. beta
!      
!.... variable declarations
      REAL*8, INTENT(IN) :: SHG(3,MEN,MINTX,*), DAB(MINTX,MINTZ,*), 
     ;                      DET(MINTX,*), XIWTS(NINTX), ETAWTS(NINTZ),
     ;                      OMEGA,PY
      INTEGER*4, INTENT(IN) :: MEN,MEE,MINTX,MINTZ, NEN,NINTX,NINTZ  
      LOGICAL*4, INTENT(IN) :: LISP, LISS 
      COMPLEX*16, INTENT(OUT) :: DKDA(MEE,*), DKDB(MEE,*) 
      REAL*8, INTENT(OUT) :: DMDA(MEE,*), DMDB(MEE,*)
!.... local variables
      COMPLEX*16 DKE11DB,DKE12DB,DKE13DB, DKE21DB,DKE22DB,DKE23DB,
     ;           DKE31DB,DKE32DB,DKE33DB, DKE11DA,DKE12DA,DKE13DA, 
     ;           DKE21DA,DKE22DA,DKE23DA, DKE31DA,DKE32DA,DKE33DA, CZERO
      REAL*8 WPY,WPY2, PHIB,PHIBX,PHIBZ, PHIA,PHIAX,PHIAZ, 
     ;       DALPHA,DBETA, 
     ;       DMASSDA,DMASSDB, DMDAT1,DMDBT1, WEIGHT,
     ;       DK11DA,DK12DA,DK13DA, DK21DA,DK22DA,DK23DA, 
     ;       DK31DA,DK32DA,DK33DA,
     ;       DK11DBT1,DK12DBT1,DK12DBT2,DK13DBT1,DK13DBT2, 
     ;       DK21DBT1,DK21DBT2,DK22DBT1,DK22DBT2,DK23DBT1,DK23DBT2,
     ;       DK31DBT1,DK31DBT2,DK32DBT1,DK32DBT2,DK33DBT1
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0))
      PARAMETER(NDIM = 3) 
! 
!----------------------------------------------------------------------!
! 
!.... null out element matrices
      NEE = NEN*NDIM
      DO 1 IPE=1,NEE
         DO 2 IQE=1,NEE
            DMDA(IPE,IQE) = 0.D0
            DMDB(IPE,IQE) = 0.D0 
            DKDA(IPE,IQE) = CZERO
            DKDB(IPE,IQE) = CZERO
    2    CONTINUE
    1 CONTINUE 
! 
!.... loop on integration points
      WPY = OMEGA*PY
      WPY2 = (OMEGA*PY)**2
      DO 3 INTX=1,NINTX
         DO 4 INTZ=1,NINTZ
            WEIGHT = XIWTS(INTX)*ETAWTS(INTZ)*DET(INTX,INTZ)
            DALPHA = DAB(INTX,INTZ,1)*WEIGHT
            DBETA  = DAB(INTX,INTZ,2)*WEIGHT
            DMASSDA = WPY2*DALPHA
            DMASSDB = WPY2*DBETA
!
!.......... loop on element nodes
            DO 5 IBE=1,NEN
               PHIB  = SHG(3,IBE,INTX,INTZ)
               PHIBX = SHG(1,IBE,INTX,INTZ)
               PHIBZ = SHG(2,IBE,INTX,INTZ)

               IF (LISP) THEN
                  DMDAT1 = DMASSDA*PHIB

                  DK11DA = DALPHA*PHIBX
                  DK13DA = DALPHA*PHIBZ
                  DK22DA = 0.D0
                  DK31DA = DALPHA*PHIBX
                  DK33DA = DALPHA*PHIBZ

                  DK12DA = WPY*DALPHA*PHIB
                  DK21DA = WPY*DALPHA*PHIBX
                  DK23DA = WPY*DALPHA*PHIBZ
                  DK32DA = WPY*DALPHA*PHIB
               ENDIF
               IF (LISS) THEN
                  DMDBT1 = DMASSDB*PHIB

                  DK11DBT1 =      DBETA*PHIBZ
                  DK13DBT1 =      DBETA*PHIBX 
                  DK13DBT2 =-2.D0*DBETA*PHIBZ
                  DK22DBT1 =      DBETA*PHIBX
                  DK22DBT2 =      DBETA*PHIBZ
                  DK31DBT1 =      DBETA*PHIBZ
                  DK31DBT2 =-2.D0*DBETA*PHIBX
                  DK33DBT1 =      DBETA*PHIBX

                  DK12DBT1 = 2.D0*WPY*DBETA*PHIB
                  DK12DBT2 =      WPY*DBETA*PHIBX
                  DK21DBT1 = 2.D0*WPY*DBETA*PHIBX
                  DK21DBT2 =      WPY*DBETA*PHIB
                  DK23DBT1 = 2.D0*WPY*DBETA*PHIBZ
                  DK23DBT2 =      WPY*DBETA*PHIB
                  DK32DBT1 = 2.D0*WPY*DBETA*PHIB
                  DK32DBT2 =      WPY*DBETA*PHIBZ
               ENDIF

               DO 6 IAE=1,NEN !IBE
                  PHIA  = SHG(3,IAE,INTX,INTZ)
                  PHIAX = SHG(1,IAE,INTX,INTZ)
                  PHIAZ = SHG(2,IAE,INTX,INTZ)
                  !indices
                  IPE = NDIM*(IAE - 1) + 1
                  IQE = NDIM*(IBE - 1) + 1
                  !dmass/dalpha
                  IF (LISP) THEN
                     DMDA(IPE  ,IQE  ) = 0.D0 
                     DMDA(IPE+1,IQE+1) = DMDA(IPE+1,IQE+1) + DMDAT1*PHIA
                     DMDA(IPE+2,IQE+2) = 0.D0
                  ENDIF
                  !dmass/dbeta
                  IF (LISS) THEN
                     DMDB(IPE  ,IQE  ) = DMDB(IPE,IQE) + DMDBT1*PHIA 
                     DMDB(IPE+1,IQE+1) = 0.D0 
                     DMDB(IPE+2,IQE+2) = DMDB(IPE,IQE) 
                  ENDIF

                  IF (LISP) THEN
                     !dstiffness/dalpha, real
                     DKE11DA = DCMPLX(PHIAX*DK11DA,0.D0) 
                     DKE13DA = DCMPLX(PHIAX*DK13DA,0.D0) 
                     DKE22DA = CZERO
                     DKE31DA = DCMPLX(PHIAZ*DK31DA,0.D0)
                     DKE33DA = DCMPLX(PHIAZ*DK33DA,0.D0)
                     !dstiffness/dalpha, imaginary
                     DKE12DA = DCMPLX(0.D0,PHIAX*DK12DA)
                     DKE21DA = DCMPLX(0.D0,PHIA *DK21DA)
                     DKE23DA = DCMPLX(0.D0,PHIA *DK23DA)
                     DKE32DA = DCMPLX(0.D0,PHIAZ*DK32DA)

                     DKDA(IPE  ,IQE  ) = DKDA(IPE  ,IQE  ) + DKE11DA
                     DKDA(IPE  ,IQE+1) = DKDA(IPE  ,IQE+1) - DKE12DA
                     DKDA(IPE  ,IQE+2) = DKDA(IPE  ,IQE+2) + DKE13DA
                     DKDA(IPE+1,IQE  ) = DKDA(IPE+1,IQE  ) + DKE21DA
                     DKDA(IPE+1,IQE+1) = DKDA(IPE+1,IQE+1) + DKE22DA
                     DKDA(IPE+1,IQE+2) = DKDA(IPE+1,IQE+2) + DKE23DA
                     DKDA(IPE+2,IQE  ) = DKDA(IPE+2,IQE  ) + DKE31DA
                     DKDA(IPE+2,IQE+1) = DKDA(IPE+2,IQE+1) - DKE32DA
                     DKDA(IPE+2,IQE+2) = DKDA(IPE+2,IQE+2) + DKE33DA
                  ENDIF

                  IF (LISS) THEN
                     !dstiffness/dalpha, real
                     DKE11DB = DCMPLX(PHIAZ*DK11DBT1,0.D0) 
                     DKE13DB = DCMPLX(PHIAZ*DK13DBT1 
     ;                              + PHIAX*DK13DBT2,0.D0)
                     DKE22DB = DCMPLX(PHIAX*DK22DBT1 
     ;                              + PHIAZ*DK22DBT2,0.D0)
                     DKE31DB = DCMPLX(PHIAX*DK31DBT1 
     ;                              + PHIAZ*DK31DBT2,0.D0)
                     DKE33DB = DCMPLX(PHIAX*DK33DBT1,0.D0)
                     !dstiffness/dalpha, imaginary
                     DKE12DB = DCMPLX(0.D0,PHIAX*DK12DBT1 
     ;                                   + PHIA *DK12DBT2)
                     DKE21DB = DCMPLX(0.D0,PHIA *DK21DBT1 
     ;                                   + PHIAX*DK21DBT2)
                     DKE23DB = DCMPLX(0.D0,PHIA *DK23DBT1 
     ;                                   + PHIAZ*DK23DBT2)
                     DKE32DB = DCMPLX(0.D0,PHIAZ*DK32DBT1 
     ;                                   + PHIA *DK32DBT2)
  
                     DKDB(IPE  ,IQE  ) = DKDB(IPE  ,IQE  ) + DKE11DB
                     DKDB(IPE  ,IQE+1) = DKDB(IPE  ,IQE+1) + DKE12DB
                     DKDB(IPE  ,IQE+2) = DKDB(IPE  ,IQE+2) + DKE13DB
                     DKDB(IPE+1,IQE  ) = DKDB(IPE+1,IQE  ) - DKE21DB
                     DKDB(IPE+1,IQE+1) = DKDB(IPE+1,IQE+1) + DKE22DB
                     DKDB(IPE+1,IQE+2) = DKDB(IPE+1,IQE+2) - DKE23DB
                     DKDB(IPE+2,IQE  ) = DKDB(IPE+2,IQE  ) + DKE31DB
                     DKDB(IPE+2,IQE+1) = DKDB(IPE+2,IQE+1) + DKE32DB
                     DKDB(IPE+2,IQE+2) = DKDB(IPE+2,IQE+2) + DKE33DB
                  ENDIF
    6          CONTINUE
    5       CONTINUE !loop on element nodes in b
    4    CONTINUE !loop on integration points in z
    3 CONTINUE !loop on integration poitns in x
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !

      SUBROUTINE DABM2KISOB25(MEN,MEE,MINTX,MINTZ, NEN,NINTX,NINTZ,
     ;                        LISP,LISS, OMEGA,PY,
     ;                        XIWTS,ETAWTS,
     ;                        DET,DAB,DAB2,SHG, DMDA,DMDB,DKDA,DKDB, 
     ;                        DMDA2,DMDB2,DKDA2,DKDB2)

      REAL*8, INTENT(IN) :: SHG(3,MEN,MINTX,*), DAB(MINTX,MINTZ,*), 
     ;                      DAB2(MINTX,MINTZ,*),  DET(MINTX,*), 
     ;                      XIWTS(NINTX), ETAWTS(NINTZ), OMEGA,PY
      INTEGER*4, INTENT(IN) :: MEN,MEE,MINTX,MINTZ, NEN,NINTX,NINTZ  
      LOGICAL*4, INTENT(IN) :: LISP, LISS 
      COMPLEX*16, INTENT(OUT) :: DKDA(MEE,*), DKDB(MEE,*), DKDA2(MEE,*),
     ;                           DKDB2(MEE,*)  
      REAL*8, INTENT(OUT) :: DMDA(MEE,*), DMDB(MEE,*), DMDA2(MEE,*), 
     ;                       DMDB2(MEE,*)
!.... local variables
      COMPLEX*16 DKE11DB,DKE12DB,DKE13DB, DKE21DB,DKE22DB,DKE23DB,
     ;           DKE31DB,DKE32DB,DKE33DB, DKE11DA,DKE12DA,DKE13DA, 
     ;           DKE21DA,DKE22DA,DKE23DA, DKE31DA,DKE32DA,DKE33DA, 
     ;           DKE11DB2,DKE12DB2,DKE13DB2, DKE21DB2,DKE22DB2,DKE23DB2,
     ;           DKE31DB2,DKE32DB2,DKE33DB2, DKE11DA2,DKE12DA2,DKE13DA2,
     ;           DKE21DA2,DKE22DA2,DKE23DA2, DKE31DA2,DKE32DA2,DKE33DA2,
     ;           CZERO
      REAL*8 WPY,WPY2, PHIB,PHIBX,PHIBZ, PHIA,PHIAX,PHIAZ, 
     ;       DALPHA,DBETA,DALPHA2,DBETA2,  
     ;       DMASSDA,DMASSDB, DMASSDA2,DMASSDB2, 
     ;       DMDAT1(2),DMDBT1(2), WEIGHT,
     ;       DK11DA(2),DK12DA(2),DK13DA(2), 
     ;       DK21DA(2),DK22DA(2),DK23DA(2), 
     ;       DK31DA(2),DK32DA(2),DK33DA(2),
     ;       DK11DBT1(2),DK12DBT1(2),DK12DBT2(2),
     ;       DK13DBT1(2),DK13DBT2(2),
     ;       DK21DBT1(2),DK21DBT2(2),DK22DBT1(2),
     ;       DK22DBT2(2),DK23DBT1(2),DK23DBT2(2),
     ;       DK31DBT1(2),DK31DBT2(2),DK32DBT1(2),
     ;       DK32DBT2(2),DK33DBT1(2)
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0))
      PARAMETER(NDIM = 3)  
! 
!----------------------------------------------------------------------!
! 
!.... null out element matrices
      NEE = NEN*NDIM
      DO 1 IPE=1,NEE
         DO 2 IQE=1,NEE
            DMDA(IPE,IQE) = 0.D0
            DMDB(IPE,IQE) = 0.D0 
            DMDA2(IPE,IQE) = 0.D0
            DMDB2(IPE,IQE) = 0.D0
            DKDA(IPE,IQE) = CZERO
            DKDB(IPE,IQE) = CZERO
            DKDA2(IPE,IQE) = CZERO
            DKDB2(IPE,IQE) = CZERO
    2    CONTINUE
    1 CONTINUE 
! 
!.... loop on integration points
      WPY = OMEGA*PY
      WPY2 = (OMEGA*PY)**2
      DO 3 INTX=1,NINTX
         DO 4 INTZ=1,NINTZ
            WEIGHT = XIWTS(INTX)*ETAWTS(INTZ)*DET(INTX,INTZ)
            DALPHA = DAB(INTX,INTZ,1)*WEIGHT
            DBETA  = DAB(INTX,INTZ,2)*WEIGHT
            DALPHA2 = DAB2(INTX,INTZ,1)*WEIGHT
            DBETA2  = DAB2(INTX,INTZ,2)*WEIGHT 
            DMASSDA = WPY2*DALPHA
            DMASSDB = WPY2*DBETA
            DMASSDA2 = WPY2*DALPHA2
            DMASSDB2 = WPY2*DBETA2 
!
!.......... loop on element nodes
            DO 5 IBE=1,NEN
               PHIB  = SHG(3,IBE,INTX,INTZ)
               PHIBX = SHG(1,IBE,INTX,INTZ)
               PHIBZ = SHG(2,IBE,INTX,INTZ)

               IF (LISP) THEN
                  DMDAT1(1) = DMASSDA*PHIB
                  !real stiffness; 1st derivatives
                  DK11DA(1) = DALPHA*PHIBX
                  DK13DA(1) = DALPHA*PHIBZ
                  DK22DA(1) = 0.D0
                  DK31DA(1) = DALPHA*PHIBX
                  DK33DA(1) = DALPHA*PHIBZ
                  !imaginary stiffness; 1st derivatives
                  DK12DA(1) = WPY*DALPHA*PHIB
                  DK21DA(1) = WPY*DALPHA*PHIBX
                  DK23DA(1) = WPY*DALPHA*PHIBZ
                  DK32DA(1) = WPY*DALPHA*PHIB
                  !mass; 2nd derivatives
                  DMDAT1(2) = DMASSDA2*PHIB
                  !real stiffness; 2nd derivatives
                  DK11DA(2) = DALPHA2*PHIBX
                  DK13DA(2) = DALPHA2*PHIBZ
                  DK22DA(2) = 0.D0
                  DK31DA(2) = DALPHA2*PHIBX
                  DK33DA(2) = DALPHA2*PHIBZ
                  !imaginary stiffness; 2nd derivatives
                  DK12DA(2) = WPY*DALPHA2*PHIB
                  DK21DA(2) = WPY*DALPHA2*PHIBX
                  DK23DA(2) = WPY*DALPHA2*PHIBZ
                  DK32DA(2) = WPY*DALPHA2*PHIB
               ENDIF
               IF (LISS) THEN
                  DMDBT1(1) = DMASSDB*PHIB
                  !real stiffness; 1st derivatives
                  DK11DBT1(1) =      DBETA*PHIBZ
                  DK13DBT1(1) =      DBETA*PHIBX
                  DK13DBT2(1) =-2.D0*DBETA*PHIBZ
                  DK22DBT1(1) =      DBETA*PHIBX
                  DK22DBT2(1) =      DBETA*PHIBZ
                  DK31DBT1(1) =      DBETA*PHIBZ
                  DK31DBT2(1) =-2.D0*DBETA*PHIBX
                  DK33DBT1(1) =      DBETA*PHIBX
                  !imaginary stiffness; 1st derivatives
                  DK12DBT1(1) = 2.D0*WPY*DBETA*PHIB
                  DK12DBT2(1) =      WPY*DBETA*PHIBX
                  DK21DBT1(1) = 2.D0*WPY*DBETA*PHIBX
                  DK21DBT2(1) =      WPY*DBETA*PHIB
                  DK23DBT1(1) = 2.D0*WPY*DBETA*PHIBZ
                  DK23DBT2(1) =      WPY*DBETA*PHIB
                  DK32DBT1(1) = 2.D0*WPY*DBETA*PHIB
                  DK32DBT2(1) =      WPY*DBETA*PHIBZ
                  !mass; 2nd derivatives
                  DMDBT1(2) = DMASSDB2*PHIB
                  !real stiffness; 1st derivatives
                  DK11DBT1(2) =      DBETA2*PHIBZ
                  DK13DBT1(2) =      DBETA2*PHIBX
                  DK13DBT2(2) =-2.D0*DBETA2*PHIBZ
                  DK22DBT1(2) =      DBETA2*PHIBX
                  DK22DBT2(2) =      DBETA2*PHIBZ
                  DK31DBT1(2) =      DBETA2*PHIBZ
                  DK31DBT2(2) =-2.D0*DBETA2*PHIBX
                  DK33DBT1(2) =      DBETA2*PHIBX
                  !imaginary stiffness; 1st derivatives
                  DK12DBT1(2) = 2.D0*WPY*DBETA2*PHIB
                  DK12DBT2(2) =      WPY*DBETA2*PHIBX
                  DK21DBT1(2) = 2.D0*WPY*DBETA2*PHIBX
                  DK21DBT2(2) =      WPY*DBETA2*PHIB
                  DK23DBT1(2) = 2.D0*WPY*DBETA2*PHIBZ
                  DK23DBT2(2) =      WPY*DBETA2*PHIB
                  DK32DBT1(2) = 2.D0*WPY*DBETA2*PHIB
                  DK32DBT2(2) =      WPY*DBETA2*PHIBZ
               ENDIF

               DO 6 IAE=1,NEN !IBE
                  PHIA  = SHG(3,IAE,INTX,INTZ)
                  PHIAX = SHG(1,IAE,INTX,INTZ)
                  PHIAZ = SHG(2,IAE,INTX,INTZ)
                  !indices
                  IPE = NDIM*(IAE - 1) + 1
                  IQE = NDIM*(IBE - 1) + 1
                  !dmass/dalpha
                  IF (LISP) THEN
                     DMDA(IPE  ,IQE  ) = 0.D0
                     DMDA(IPE+1,IQE+1) = DMDA(IPE+1,IQE+1) 
     ;                                 + DMDAT1(1)*PHIA
                     DMDA(IPE+2,IQE+2) = 0.D0

                     DMDA2(IPE  ,IQE  ) = 0.D0
                     DMDA2(IPE+1,IQE+1) = DMDA2(IPE+1,IQE+1)
     ;                                  + DMDAT1(2)*PHIA
                     DMDA2(IPE+2,IQE+2) = 0.D0
                  ENDIF
                  !dmass/dbeta
                  IF (LISS) THEN
                     DMDB(IPE  ,IQE  ) = DMDB(IPE,IQE) + DMDBT1(1)*PHIA
                     DMDB(IPE+1,IQE+1) = 0.D0
                     DMDB(IPE+2,IQE+2) = DMDB(IPE,IQE)

                     DMDB2(IPE  ,IQE  ) = DMDB2(IPE,IQE) 
     ;                                  + DMDBT1(2)*PHIA
                     DMDB2(IPE+1,IQE+1) = 0.D0
                     DMDB2(IPE+2,IQE+2) = DMDB2(IPE,IQE)
                  ENDIF

                  IF (LISP) THEN
                     !dstiffness/dalpha, real
                     DKE11DA = DCMPLX(PHIAX*DK11DA(1),0.D0)
                     DKE13DA = DCMPLX(PHIAX*DK13DA(1),0.D0)
                     DKE22DA = CZERO
                     DKE31DA = DCMPLX(PHIAZ*DK31DA(1),0.D0)
                     DKE33DA = DCMPLX(PHIAZ*DK33DA(1),0.D0)
                     !dstiffness/dalpha, imaginary
                     DKE12DA = DCMPLX(0.D0,PHIAX*DK12DA(1))
                     DKE21DA = DCMPLX(0.D0,PHIA *DK21DA(1))
                     DKE23DA = DCMPLX(0.D0,PHIA *DK23DA(1))
                     DKE32DA = DCMPLX(0.D0,PHIAZ*DK32DA(1))

                     DKDA(IPE  ,IQE  ) = DKDA(IPE  ,IQE  ) + DKE11DA
                     DKDA(IPE  ,IQE+1) = DKDA(IPE  ,IQE+1) - DKE12DA
                     DKDA(IPE  ,IQE+2) = DKDA(IPE  ,IQE+2) + DKE13DA
                     DKDA(IPE+1,IQE  ) = DKDA(IPE+1,IQE  ) + DKE21DA
                     DKDA(IPE+1,IQE+1) = DKDA(IPE+1,IQE+1) + DKE22DA
                     DKDA(IPE+1,IQE+2) = DKDA(IPE+1,IQE+2) + DKE23DA
                     DKDA(IPE+2,IQE  ) = DKDA(IPE+2,IQE  ) + DKE31DA
                     DKDA(IPE+2,IQE+1) = DKDA(IPE+2,IQE+1) - DKE32DA
                     DKDA(IPE+2,IQE+2) = DKDA(IPE+2,IQE+2) + DKE33DA

                     !dstiffness/dalpha, real
                     DKE11DA2 = DCMPLX(PHIAX*DK11DA(2),0.D0)
                     DKE13DA2 = DCMPLX(PHIAX*DK13DA(2),0.D0)
                     DKE22DA2 = CZERO
                     DKE31DA2 = DCMPLX(PHIAZ*DK31DA(2),0.D0)
                     DKE33DA2 = DCMPLX(PHIAZ*DK33DA(2),0.D0)
                     !dstiffness/dalpha, imaginary
                     DKE12DA2 = DCMPLX(0.D0,PHIAX*DK12DA(2))
                     DKE21DA2 = DCMPLX(0.D0,PHIA *DK21DA(2))
                     DKE23DA2 = DCMPLX(0.D0,PHIA *DK23DA(2))
                     DKE32DA2 = DCMPLX(0.D0,PHIAZ*DK32DA(2))

                     DKDA2(IPE  ,IQE  ) = DKDA2(IPE  ,IQE  ) + DKE11DA2
                     DKDA2(IPE  ,IQE+1) = DKDA2(IPE  ,IQE+1) - DKE12DA2
                     DKDA2(IPE  ,IQE+2) = DKDA2(IPE  ,IQE+2) + DKE13DA2
                     DKDA2(IPE+1,IQE  ) = DKDA2(IPE+1,IQE  ) + DKE21DA2
                     DKDA2(IPE+1,IQE+1) = DKDA2(IPE+1,IQE+1) + DKE22DA2
                     DKDA2(IPE+1,IQE+2) = DKDA2(IPE+1,IQE+2) + DKE23DA2
                     DKDA2(IPE+2,IQE  ) = DKDA2(IPE+2,IQE  ) + DKE31DA2
                     DKDA2(IPE+2,IQE+1) = DKDA2(IPE+2,IQE+1) - DKE32DA2
                     DKDA2(IPE+2,IQE+2) = DKDA2(IPE+2,IQE+2) + DKE33DA2
                  ENDIF
                  IF (LISS) THEN
                     !dstiffness/dalpha, real
                     DKE11DB = DCMPLX(PHIAZ*DK11DBT1(1),0.D0)
                     DKE13DB = DCMPLX(PHIAZ*DK13DBT1(1)
     ;                              + PHIAX*DK13DBT2(1),0.D0)
                     DKE22DB = DCMPLX(PHIAX*DK22DBT1(1)
     ;                              + PHIAZ*DK22DBT2(1),0.D0)
                     DKE31DB = DCMPLX(PHIAX*DK31DBT1(1)
     ;                              + PHIAZ*DK31DBT2(1),0.D0)
                     DKE33DB = DCMPLX(PHIAX*DK33DBT1(1),0.D0)
                     !dstiffness/dalpha, imaginary
                     DKE12DB = DCMPLX(0.D0,PHIAX*DK12DBT1(1)
     ;                                   + PHIA *DK12DBT2(1))
                     DKE21DB = DCMPLX(0.D0,PHIA *DK21DBT1(1)
     ;                                   + PHIAX*DK21DBT2(1))
                     DKE23DB = DCMPLX(0.D0,PHIA *DK23DBT1(1)
     ;                                   + PHIAZ*DK23DBT2(1))
                     DKE32DB = DCMPLX(0.D0,PHIAZ*DK32DBT1(1)
     ;                                   + PHIA *DK32DBT2(1))

                     DKDB(IPE  ,IQE  ) = DKDB(IPE  ,IQE  ) + DKE11DB
                     DKDB(IPE  ,IQE+1) = DKDB(IPE  ,IQE+1) + DKE12DB
                     DKDB(IPE  ,IQE+2) = DKDB(IPE  ,IQE+2) + DKE13DB
                     DKDB(IPE+1,IQE  ) = DKDB(IPE+1,IQE  ) - DKE21DB
                     DKDB(IPE+1,IQE+1) = DKDB(IPE+1,IQE+1) + DKE22DB
                     DKDB(IPE+1,IQE+2) = DKDB(IPE+1,IQE+2) - DKE23DB
                     DKDB(IPE+2,IQE  ) = DKDB(IPE+2,IQE  ) + DKE31DB
                     DKDB(IPE+2,IQE+1) = DKDB(IPE+2,IQE+1) + DKE32DB
                     DKDB(IPE+2,IQE+2) = DKDB(IPE+2,IQE+2) + DKE33DB

                     !dstiffness/dalpha^2, real
                     DKE11DB2 = DCMPLX(PHIAZ*DK11DBT1(2),0.D0)
                     DKE13DB2 = DCMPLX(PHIAZ*DK13DBT1(2)
     ;                               + PHIAX*DK13DBT2(2),0.D0)
                     DKE22DB2 = DCMPLX(PHIAX*DK22DBT1(2)
     ;                               + PHIAZ*DK22DBT2(2),0.D0)
                     DKE31DB2 = DCMPLX(PHIAX*DK31DBT1(2)
     ;                               + PHIAZ*DK31DBT2(2),0.D0)
                     DKE33DB2 = DCMPLX(PHIAX*DK33DBT1(2),0.D0)
                     !dstiffness/dalpha^2, imaginary
                     DKE12DB2 = DCMPLX(0.D0,PHIAX*DK12DBT1(2)
     ;                                    + PHIA *DK12DBT2(2))
                     DKE21DB2 = DCMPLX(0.D0,PHIA *DK21DBT1(2)
     ;                                    + PHIAX*DK21DBT2(2))
                     DKE23DB2 = DCMPLX(0.D0,PHIA *DK23DBT1(2)
     ;                                    + PHIAZ*DK23DBT2(2))
                     DKE32DB2 = DCMPLX(0.D0,PHIAZ*DK32DBT1(2)
     ;                                    + PHIA *DK32DBT2(2))

                     DKDB2(IPE  ,IQE  ) = DKDB2(IPE  ,IQE  ) + DKE11DB2
                     DKDB2(IPE  ,IQE+1) = DKDB2(IPE  ,IQE+1) + DKE12DB2
                     DKDB2(IPE  ,IQE+2) = DKDB2(IPE  ,IQE+2) + DKE13DB2
                     DKDB2(IPE+1,IQE  ) = DKDB2(IPE+1,IQE  ) - DKE21DB2
                     DKDB2(IPE+1,IQE+1) = DKDB2(IPE+1,IQE+1) + DKE22DB2
                     DKDB2(IPE+1,IQE+2) = DKDB2(IPE+1,IQE+2) - DKE23DB2
                     DKDB2(IPE+2,IQE  ) = DKDB2(IPE+2,IQE  ) + DKE31DB2
                     DKDB2(IPE+2,IQE+1) = DKDB2(IPE+2,IQE+1) + DKE32DB2
                     DKDB2(IPE+2,IQE+2) = DKDB2(IPE+2,IQE+2) + DKE33DB2
                  ENDIF
    6          CONTINUE
    5       CONTINUE !loop on element nodes in b
    4    CONTINUE !loop on integration points in z
    3 CONTINUE !loop on integration poitns in x
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DKM_VPVS(MEN,MEE,MINTX, NEN,NINTX,NINTZ,
     ;                    OMEGA,PY, XIWTS,ETAWTS,
     ;                    DVP,DVS, DET,SHG, DME,DKE)
!
!     2.5D isotropic elastic element assembly with P and S wave 
!     velocities.  Useful for inversion  - B. Baker May 2013
!
!     INPUT      MEANING
!     -----      ------- 
!     DET        jacobian at integration points
!     DVP        derivatives w.r.t. p velocity 
!     DVS        derivatives w.r.t. p velocity/vpvs^2
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
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     DKE        derivative element stiffness 
!     DME        derivative element mass
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: SHG(3,MEN,MINTX,*),DVP(MINTX,*), 
     ;                      DVS(MINTX,*),
     ;                      DET(MINTX,*), XIWTS(NINTX),ETAWTS(NINTZ), 
     ;                      OMEGA,PY
      INTEGER*4, INTENT(IN) :: MEN, MEE, MINTX, NEN,NINTX,NINTZ
      COMPLEX*16, INTENT(OUT) :: DKE(MEE,*)
      REAL*8, INTENT(OUT) :: DME(MEE,*)
!.... local variables
      COMPLEX*16 ZZERO 
      REAL*8 WPY, PY2, W2, WEIGHT, DRVP,DRVS,
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
      PARAMETER(NDIM = 3)
! 
!----------------------------------------------------------------------!
! 
!.... null out element matrices
      NEE = NEN*NDIM
      DO 1 IPE=1,NEE
         DO 2 IQE=1,NEE
            DME(IPE,IQE) = 0.D0
            DKE(IPE,IQE) = ZZERO
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
            DRVP = DVP(INTX,INTZ)*WEIGHT
            DRVS = DVS(INTX,INTZ)*WEIGHT

            RMCOEFF1 = W2*PY2*DRVS
            RMCOEFF2 = W2*PY2*DRVP
            RMCOEFF3 = RMCOEFF1
            DO 5 IBE=1,NEN
               PHIB  = SHG(3,IBE,INTX,INTZ)
               PHIBX = SHG(1,IBE,INTX,INTZ)
               PHIBZ = SHG(2,IBE,INTX,INTZ)

               RM11T1 = RMCOEFF1*PHIB
               RM22T1 = RMCOEFF2*PHIB
               RM33T1 = RM11T1

               RK11T1 = DRVP*PHIBX
               RK11T2 = DRVS*PHIBZ

               RK13T1 = DRVP*PHIBZ
               RK13T2 = DRVS*PHIBX
               RK13T3 =-2.D0*DRVS*PHIBZ

               RK22T1 = DRVS*PHIBX
               RK22T2 = DRVS*PHIBZ

               RK31T1 = DRVP*PHIBX
               RK31T2 = DRVS*PHIBZ
               RK31T3 =-2.D0*DRVS*PHIBX

               RK33T1 = DRVP*PHIBZ
               RK33T2 = DRVS*PHIBX

               RK12T1 = WPY*DRVP*PHIB
               RK12T2 =-2.D0*WPY*DRVS*PHIB
               RK12T3 =     -WPY*DRVS*PHIBX

               RK21T1 = WPY*DRVP*PHIBX
               RK21T2 =-2.D0*WPY*DRVS*PHIBX
               RK21T3 =     -WPY*DRVS*PHIB

               RK23T1 = WPY*DRVP*PHIBZ
               RK23T2 =-2.D0*WPY*DRVS*PHIBZ
               RK23T3 =     -WPY*DRVS*PHIB

               RK32T1 = WPY*DRVP*PHIB
               RK32T2 =-2.D0*WPY*DRVS*PHIB
               RK32T3 =     -WPY*DRVS*PHIBZ
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
                  DME(IPE  ,IQE  ) = DME(IPE  ,IQE  ) + RME11
                  DME(IPE+1,IQE+1) = DME(IPE+1,IQE+1) + RME22
                  DME(IPE+2,IQE+2) = DME(IPE+2,IQE+2) + RME33
                  !main diagonal stiffenss
                  DKE(IPE  ,IQE  )=DKE(IPE  ,IQE  ) + DCMPLX(RKE11,0.D0)
                  DKE(IPE+1,IQE+1)=DKE(IPE+1,IQE+1) + DCMPLX(RKE22,0.D0)
                  DKE(IPE+2,IQE+2)=DKE(IPE+2,IQE+2) + DCMPLX(RKE33,0.D0)
                  !real off diagonal 
                  DKE(IPE  ,IQE+2)=DKE(IPE  ,IQE+2) + DCMPLX(RKE13,0.D0)
                  DKE(IPE+2,IQE  )=DKE(IPE+2,IQE  ) + DCMPLX(RKE31,0.D0)
                  !complex off diagonal
                  DKE(IPE  ,IQE+1)=DKE(IPE  ,IQE+1) - DCMPLX(RKE12,0.D0)
                  DKE(IPE+1,IQE  )=DKE(IPE+1,IQE  ) + DCMPLX(RKE21,0.D0)
                  DKE(IPE+1,IQE+2)=DKE(IPE+1,IQE+2) + DCMPLX(RKE23,0.D0)
                  DKE(IPE+2,IQE+1)=DKE(IPE+2,IQE+1) - DCMPLX(RKE32,0.D0)
    6          CONTINUE !loop on iae
    5       CONTINUE !Loop on ibe
    4    CONTINUE !loop on eta integration
    3 CONTINUE !loop on xi integration
      RETURN
      END

