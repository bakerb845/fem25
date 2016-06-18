      SUBROUTINE VFORCE25(MNZROW,MEN,MINTX,MINTZ,MNPG,MDIM,MGNOD,
     ;                           NEN,NINTX,NINTZ,NNPG,     NGNOD,
     ;                    NNPGL,NCON,NVINV, NZ_FDIST,NDOF,NELEM, 
     ;                    CINVTYPE,LISISO, INPG,INPGL, 
     ;                    JCSC_FDIST,ICSC_FDIST, IENG,MCONN,LM,
     ;                    OMEGA,PY,VPVS, XLOCS,ZLOCS,DENS,ECOEFF,
     ;                    XIGLL,ETAGLL,SHL, ELEM_WTS, WAVE,
     ;                    FORCEV, IERR)
!
!     Calculates the virtural forces where the i'th virtural force 
!     which is nominally of length ndof is given by dS/dm_i u.  However,
!     since S is sparse consequently dS/dm_i is sparse. dS/dm_i is 
!     then stored in CSC format 
!
!     INPUT      MEANING
!     -----      ------- 
! 
!     OUTPUT     MEANING 
!     ------     -------
!     FORCEV     columns of virtural force matrix corresponding INPG 
!     IERR       error flag
!
!.... variable declarations
      CHARACTER(2), INTENT(IN) :: CINVTYPE 
      COMPLEX*8, INTENT(IN) :: WAVE(NDOF) 
      REAL*8, INTENT(IN) :: SHL(3,MEN,MINTX,*), ECOEFF(MNPG,*), 
     ;                      ETAGLL(MINTX,*), XIGLL(MINTZ,*), 
     ;                      XLOCS(NNPG), ZLOCS(NNPG), DENS(NNPG), 
     ;                      ELEM_WTS(NELEM), OMEGA, PY, VPVS
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), MCONN(MNPG,*), 
     ;           IENG(MGNOD,*),JCSC_FDIST(NNPGL+1),ICSC_FDIST(NZ_FDIST),
     ;           MNZROW,MEN,MINTX,MINTZ,MNPG,MDIM,MGNOD, 
     ;                  NEN,NINTX,NINTZ,NNPG,     NGNOD,
     ;           NNPGL,NCON,NVINV,  NZ_FDIST,NDOF, INPG, INPGL  
      LOGICAL*4, INTENT(IN) :: LISISO
      COMPLEX*8, INTENT(OUT) :: FORCEV(MNZROW,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*16, ALLOCATABLE :: DKDA(:,:), DKDB(:,:)
      COMPLEX*8, ALLOCATABLE :: CARG(:)
      REAL*8, ALLOCATABLE :: SHG(:,:,:,:), DAB(:,:,:), DET(:,:),
     ;                       DMDA(:,:), DMDB(:,:)
      COMPLEX*8 CZERO, CWGHT
      INTEGER*4 NDIM,NEE,NNZROW, I1,I2,IDOF,JDOF, IAE,I,IPE, IBE,J,IQE,
     ;          IVINV, ICON, IELEM, INZROW  
      INTEGER*4 IBSECT 
      LOGICAL*4 LISP, LISS, LDEBUG
      PARAMETER(CZERO = CMPLX(0.0,0.0))
      PARAMETER(NDIM = 3)
      PARAMETER(LDEBUG = .FALSE.) 
      COMPLEX*16, ALLOCATABLE :: KE1(:,:), KE(:,:), DKDP_FD(:,:) 
      REAL*8, ALLOCATABLE :: ME1(:,:), ME(:,:), DMDP_FD(:,:), 
     ;                       VP(:,:), VS(:,:), QPP(:,:), QSS(:,:), 
     ;                       DVP(:,:), DVS(:,:), VPANC(:), VSANC(:) 
      REAL*8 FREQ, FREQ0, TWOPI, VPSAVE, VSSAVE, SMALL, VPVS2 
      INTEGER*4 JNPG
      PARAMETER(TWOPI = 6.2831853071795862D0)
      PARAMETER(SMALL = 1.D-7) 
!     PARAMETER(VPVS2 = 3.D0) !assume a poisson solid Vs = Vp/sqrt(3)
!
!----------------------------------------------------------------------!
!
!.... first error check anisotropy
      IERR = 0 
      VPVS2 = VPVS**2
      IF (CINVTYPE.NE.'PP' .AND. CINVTYPE.NE.'pp' .AND. 
     ;    CINVTYPE.NE.'SS' .AND. CINVTYPE.NE.'ss' .AND. 
     ;    CINVTYPE.NE.'PS' .AND. CINVTYPE.NE.'ps') THEN
         WRITE(*,*) 'vforce25: Error anisotropy not programmed!'
         IERR = 1 
         RETURN
      ENDIF
      LISP = .FALSE.
      LISS = .FALSE.
      IF (CINVTYPE.EQ.'PP' .OR. CINVTYPE.EQ.'pp') LISP = .TRUE.
      IF (CINVTYPE.EQ.'SS' .OR. CINVTYPE.EQ.'ss') LISS = .TRUE.
      IF (CINVTYPE.EQ.'PS' .OR. CINVTYPE.EQ.'ps') THEN
         LISP = .TRUE.
         LISS = .TRUE.
      ENDIF
      IF (LISP .AND. LISS) THEN
         WRITE(*,*) 'vforce25: Error S and P derivatives not done yet'
         IERR = 1
         RETURN
      ENDIF
!     IF (LISS) THEN
!        WRITE(*,*) 'vforce25: Error S derivatives not yet programmed!'
!        IERR = 1
!        RETURN
!     ENDIF
!
!.... set space
      NEE = NDIM*NEN
      ALLOCATE(SHG(3,NEN,NINTX,NINTZ))
      IF (NVINV.LE.2) THEN
         ALLOCATE(DAB(NINTX,NINTZ,2))
      ELSE
         WRITE(*,*) 'vforce25: Really dont know what do w/ anisotropy' 
         ALLOCATE(DAB(NINTX,NINTZ,5)) 
      ENDIF 
      ALLOCATE(DET(NINTX,NINTZ))
      ALLOCATE(DMDA(NEE,NEE))
      ALLOCATE(DMDB(NEE,NEE))
      ALLOCATE(DKDA(NEE,NEE))
      ALLOCATE(DKDB(NEE,NEE))
      ALLOCATE(CARG(NVINV))
      ALLOCATE(DVP(NINTX,NINTZ))
      ALLOCATE(DVS(NINTX,NINTZ)) 
      DMDA(:,:) = 0.D0
      DMDB(:,:) = 0.D0
      DKDA(:,:) = DCMPLX(0.D0,0.D0)
      DKDB(:,:) = DCMPLX(0.D0,0.D0)
!
!.... debugging
      FREQ = OMEGA/TWOPI
      IF (LDEBUG) THEN
         ALLOCATE(VPANC(NNPG))
         ALLOCATE(VSANC(NNPG)) 
         ALLOCATE(VP(NINTX,NINTZ))
         ALLOCATE(VS(NINTX,NINTZ)) 
         ALLOCATE(ME1(NEE,NEE))
         ALLOCATE(KE1(NEE,NEE)) 
         ALLOCATE(ME (NEE,NEE))
         ALLOCATE(KE (NEE,NEE)) 
         ALLOCATE(DMDP_FD(NEE,NEE))
         ALLOCATE(DKDP_FD(NEE,NEE)) 
         ALLOCATE(QPP(NINTX,NINTZ))
         ALLOCATE(QSS(NINTX,NINTZ)) 
         DO JNPG=1,NNPG
            VPANC(JNPG) = (ECOEFF(JNPG,1) + 2.D0*ECOEFF(JNPG,2))
     ;                   /DENS(JNPG)
            VSANC(JNPG) = ECOEFF(JNPG,2)/DENS(JNPG)
            VPANC(JNPG) = DSQRT(VPANC(JNPG))
            VSANC(JNPG) = DSQRT(VSANC(JNPG)) 
         ENDDO 
         QPP(:,:) = 9999.D0 !turn off damping
         QSS(:,:) = 9999.D0 !turn off damping
         FREQ = OMEGA/TWOPI
         FREQ0 = 0.D0 !turn off damping
      ENDIF
!
!.... loop on element connectivity
      JCON = 0 
      DO 200 ICON=1,NCON
         IELEM = MCONN(INPG,ICON) 
         IF (IELEM.LE.0) GOTO 250 !done with connections
         JCON = JCON + 1 
!
!....... integrate element
         IF (LISISO) THEN !P and S gradients 
!
!.......... evaluate model derivatives at integration points
            CALL DJACAB(NEN,NINTX,NINTZ,NNPG, LISISO,              
     ;                  NINTX,NINTZ,NEN,NGNOD,NNPG, INPG,   
     ;                  IENG(1:NGNOD,IELEM),   
     ;                  XIGLL(1:NINTX,1),ETAGLL(1:NINTZ,1), 
     ;                  XLOCS,ZLOCS,DENS,ECOEFF(1:NNPG,1:2),
     ;                  SHL(1:3,1:NEN,1:NINTX,1:NINTZ),
     ;                  DET,DAB,SHG,IERR) 
            IF (IERR /= 0) THEN
               WRITE(*,*) 'vforce25: Error calling djacab'
               GOTO 730 
            ENDIF
            IF (LISP .AND. .NOT.LISS) THEN
               DO INTX=1,NINTX
                  DO INTZ=1,NINTZ
                     DVP(INTX,INTZ) = DAB(INTX,INTZ,1) 
                     DVS(INTX,INTZ) = DAB(INTX,INTZ,1)/VPVS2
                  ENDDO
               ENDDO   
            ENDIF
            IF (LISS .AND. .NOT.LISP) THEN
               DO INTX=1,NINTX
                  DO INTZ=1,NINTZ
                     DVP(INTX,INTZ) = DAB(INTX,INTZ,2)*VPVS2 
                     DVS(INTX,INTZ) = DAB(INTX,INTZ,2)
                  ENDDO
               ENDDO   
            ENDIF
            CALL DKM_VPVS(NEN,NEE,NINTX, NEN,NINTX,NINTZ, 
     ;                    OMEGA,PY, XIGLL(1:NINTX,2),ETAGLL(1:NINTZ,2),
     ;                    DVP,DVS, DET,SHG, DMDA,DKDA)  
!
!.......... numerically integrate element matrices
!           CALL DABMKISOB25(NEN,NEE,NINTX,NINTZ, NEN,NINTX,NINTZ,   
!    ;                       LISP,LISS,OMEGA,PY,  
!    ;                       XIGLL(1:NINTX,2),ETAGLL(1:NINTZ,2),
!    ;                       DET,DAB,SHG, DMDA,DMDB,DKDA,DKDB)  
         ELSE  
            WRITE(*,*) 'vforce25: No anisotropy yet!'
            IERR = 1 
            GOTO 730 
         ENDIF 
!
!....... debug
         IF (LDEBUG) THEN
            IF (LISP) THEN 
               CALL CJMK_VPVS(NEN,NINTX,
     ;                       NINTX,NINTZ,NEN,NGNOD,NNPG,
     ;                       IENG(1:NGNOD,IELEM),
     ;                       XIGLL(1:NINTX,1),ETAGLL(1:NINTZ,1),
     ;                       XLOCS,ZLOCS,VPANC,VSANC,SHL,
     ;                       DET,VP,VS,SHG,IERR)
               CALL KMDAMP(NEN,NEE,NINTX, NEN,NINTX,NINTZ, 
     ;                     FREQ,FREQ0,PY,  
     ;                     XIGLL(1:NINTX,2),ETAGLL(1:NINTZ,2), VP,VS,
     ;                     QPP,QSS, DET,SHG, ME,KE)
               VPSAVE = VPANC(INPG) 
               VPANC(INPG) = VPANC(INPG) + SMALL
               CALL CJMK_VPVS(NEN,NINTX,
     ;                       NINTX,NINTZ,NEN,NGNOD,NNPG, 
     ;                       IENG(1:NGNOD,IELEM),
     ;                       XIGLL(1:NINTX,1),ETAGLL(1:NINTZ,1),
     ;                       XLOCS,ZLOCS,VPANC,VSANC,SHL, 
     ;                       DET,VP,VS,SHG,IERR) 
               CALL KMDAMP(NEN,NEE,NINTX, NEN,NINTX,NINTZ, 
     ;                     FREQ,FREQ0,PY,  
     ;                     XIGLL(1:NINTX,2),ETAGLL(1:NINTZ,2), VP,VS, 
     ;                     QPP,QSS, DET,SHG, ME1,KE1) 
               VPANC(INPG) = VPSAVE 
               CALL FD_EMATS_25(NEE,NEN, SMALL,ME,ME1, KE,KE1,
     ;                          DMDP_FD,DKDP_FD)
               DO IPE=1,NEE
                  DO IQE=1,NEE
                     PRINT *, DMDP_FD(IPE,IQE), DMDA(IPE,IQE)
c                    PRINT *, CMPLX(DKDP_FD(IPE,IQE)),
c    ;                        CMPLX(DKDA(IPE,IQE)) 
c                    PRINT *, SNGL(DMDP_FD(IPE,IQE) - DMDA(IPE,IQE)),
c    ;                        CMPLX(DKDP_FD(IPE,IQE) - DKDA(IPE,IQE))
                  ENDDO
               ENDDO
              PAUSE
            ENDIF
            IF (LISS) THEN
               CALL CJMK_VPVS(NEN,NINTX,
     ;                       NINTX,NINTZ,NEN,NGNOD,NNPG,
     ;                       IENG(1:NGNOD,IELEM),
     ;                       XIGLL(1:NINTX,1),ETAGLL(1:NINTZ,1),
     ;                       XLOCS,ZLOCS,VPANC,VSANC,SHL,
     ;                       DET,VP,VS,SHG,IERR)
               CALL KMDAMP(NEN,NEE,NINTX, NEN,NINTX,NINTZ, 
     ;                     FREQ,FREQ0,PY,  
     ;                     XIGLL(1:NINTX,2),ETAGLL(1:NINTZ,2), VP,VS,
     ;                     QPP,QSS, DET,SHG, ME,KE)
               VSSAVE = VSANC(INPG)
               VSANC(INPG) = VSANC(INPG) + SMALL
               CALL CJMK_VPVS(NEN,NINTX,
     ;                       NINTX,NINTZ,NEN,NGNOD,NNPG,
     ;                       IENG(1:NGNOD,IELEM),
     ;                       XIGLL(1:NINTX,1),ETAGLL(1:NINTZ,1),
     ;                       XLOCS,ZLOCS,VPANC,VSANC,SHL,
     ;                       DET,VP,VS,SHG,IERR)
               CALL KMDAMP(NEN,NEE,NINTX, NEN,NINTX,NINTZ,
     ;                     FREQ,FREQ0,PY,
     ;                     XIGLL(1:NINTX,2),ETAGLL(1:NINTZ,2), VP,VS, 
     ;                     QPP,QSS, DET,SHG, ME1,KE1) 
               VSANC(INPG) = VSSAVE
               CALL FD_EMATS_25(NEE,NEN, SMALL,ME,ME1, KE,KE1,
     ;                          DMDP_FD,DKDP_FD)
               WRITE(*,*) 
               DO IPE=1,NEE
                  DO IQE=1,NEE
                     PRINT *, DMDP_FD(IPE,IQE),DMDB(IPE,IQE)
c                    PRINT *, CMPLX(DKDP_FD(IPE,IQE)),
c    ;                        CMPLX(DKDB(IPE,IQE))
c                    PRINT *, SNGL(DMDP_FD(IPE,IQE) - DMDB(IPE,IQE)),
c    ;                        CMPLX(DKDP_FD(IPE,IQE) - DKDPB(IPE,IQE))
                  ENDDO
               ENDDO
               PAUSE
            ENDIF
         ENDIF !end check on debugging
!
!....... null out these columns of adj(F) 
         I1 = JCSC_FDIST(INPGL)
         I2 = JCSC_FDIST(INPGL+1) - 1
         NNZROW = I2 - I1 + 1
         DO 150 IVINV=1,NVINV
            FORCEV(1:NNZROW,IVINV) = CZERO
  150    CONTINUE

         CWGHT = CMPLX(SNGL(ELEM_WTS(IELEM)),0.0)
!
!....... loop on rows
         IPE = 0
         DO 201 IAE=1,NEN
            DO 202 I=1,NDIM
               IPE = IPE + 1
               IDOF = LM(I,IAE,IELEM)
               IF (IDOF.LE.0) GOTO 220 !not a dof
               !get the non-zero row index
               INZROW = IBSECT(NNZROW,IDOF,ICSC_FDIST(I1:I2)) 
               IF (INZROW.LE.0 .OR. INZROW.GT.NNZROW) THEN
                  WRITE(*,*) 'vforce: Cant find anchor node!'
                  IERR = 1
                  GOTO 730
               ENDIF
!
!............. loop on columns of S to perform row vector multiply
               IQE = 0
               DO 203 IBE=1,NEN
                  DO 204 J=1,NDIM
                     IQE = IQE + 1
                     JDOF = LM(J,IBE,IELEM)
                     IF (JDOF.LE.0) GOTO 240 !not a dof 
                     CARG(1:NVINV) = CZERO 
                     IF (LISISO) THEN
                        IVINV = 0
                        IF (LISP) THEN !dS/dalpha
                           IVINV = IVINV + 1
                           CARG(IVINV)= CMPLX(DCMPLX(DMDA(IPE,IQE),0.D0)
     ;                                             + DKDA(IPE,IQE))
                        ENDIF
                        IF (LISS) THEN !dS/dbeta
                           IVINV = IVINV + 1
                           CARG(IVINV)= CMPLX(DCMPLX(DMDA(IPE,IQE),0.D0)
     ;                                             + DKDA(IPE,IQE))
                        ENDIF
                     ELSE
                        WRITE(*,*) 'adjfv_dist25_2: No anisotropy!'
                        IERR = 1
                        RETURN
                     ENDIF
!
!................... multiply dS/dm_i u and save to virtual force f_i   
                     DO 205 IVINV=1,NVINV
                        FORCEV(INZROW,IVINV) = FORCEV(INZROW,IVINV) 
     ;                                       + CARG(IVINV)*WAVE(JDOF)
     ;                                       *CWGHT
  205                CONTINUE !loop on components
  240                CONTINUE !not a DOF
  204             CONTINUE !loop on components
  203          CONTINUE !loop on element nodes 
  220          CONTINUE !not a DOF
  202       CONTINUE !loop on components
  201    CONTINUE !loop on element nodes
  200 CONTINUE !loop on connectivity 
  250 CONTINUE !done with connections
  730 CONTINUE !break ahead for an error
!     IF (JCON > 0) FORCEV(:,1:NVINV) = FORCEV(:,1:NVINV)
!    ;                                 /CMPLX(FLOAT(JCON),0.0)
!
!.... free space
      DEALLOCATE(SHG)
      DEALLOCATE(DAB)
      DEALLOCATE(DET)
      DEALLOCATE(DMDA)
      DEALLOCATE(DMDB)
      DEALLOCATE(DKDA)
      DEALLOCATE(DKDB)
      DEALLOCATE(CARG)
      IF (LDEBUG) THEN
         DEALLOCATE(VPANC)
         DEALLOCATE(VSANC) 
         DEALLOCATE(VP)
         DEALLOCATE(VS) 
         DEALLOCATE(ME1)
         DEALLOCATE(KE1) 
         DEALLOCATE(ME )
         DEALLOCATE(KE ) 
         DEALLOCATE(DMDP_FD)
         DEALLOCATE(DKDP_FD)
         DEALLOCATE(QPP)
         DEALLOCATE(QSS)
      ENDIF
      DEALLOCATE(DVP) 
      DEALLOCATE(DVS) 
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE VFORCE25_2(MNZROW,MDIM,MINTX,MINTZ,MGNOD,MEN,MNPG,  
     ;                     NNPG,NDOF, NINTX,NINTZ,NGNOD,NEN,NELEMG, 
     ;                     NVINV, CINVTYPE,LISISO,INPG, MYGELEM,IENG,LM,
     ;                     OMEGA,PY, XLOCS,ZLOCS,DENS,ECOEFF, 
     ;                     XIGLL,ETAGLL,SHL, WAVE, 
     ;                     NNZROW,MYROWS,FORCEV, IERR)
!
!     Calculates the virtual forces where the i'th virtual force
!     which is nominally of length ndof is given by dS/dm_i u.  
!     However, since S is sparse consequently dS/dm_i is sparse
!
!     INTPUT     MEANING
!     ------     ------- 
!     CINVTYPE   'PP' -> alpha velocity inversion
!                'SS  -> beta velocity inversion
!                'PS' -> alpha and beta velocity inversion
!     DENS       density at anchor noes
!     ECOEFF     elastic coefficients at anchor nodes
!     ETAGLL     eta GLL abscissas and integration weights
!     IENG       element to global anchor node pointer
!     INPG       anchor node to which we are generating a virtual force
!     LISISO     True -> isotropic simulation
!     LM         element to global DOF pointer
!     MDIM       leading dimension for LM
!     MEN        leading dimension for LM and SHL
!     MGNOD      leading dimension for IENG
!     MINTX      leading dimension for SHL and XIGLL
!     MINTZ      leading dimension for ETAGLL
!     MNPG       leading dimension for ECOEFF
!     MYGELEM    elements attached to INPG for gradient calculation
!     MNZROW     leading dimension for FORCEV
!     NDOF       number of degrees of freedom
!     NELEMG     number of elements in MYGELEMG
!     NEN        number of element nodes
!     NGNOD      number of anchor nodes
!     NINTX      number of integration points in xi 
!     NINTZ      number of integration points in eta
!     NNPG       number of anchor nodes in mesh
!     NVINV      number of variables in inverse problem
!     OMEGA      angular frequency
!     PY         apparaent slowness in y 
!     SHL        element lagrange shape fns
!     WAVE       wavefield from solving Su = f 
!     XIGLL      xi GLL abscissas and integration weights
!     XLOCS      x locations of anchor nodes 
!     ZLOCS      z locations of anchor nodes 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag 
!     FORCEV     virtual forces for each model derivative, sparse
!     MYROWS     global location of non-zero rows in f_i 
!     NNZROW     number of non-zeros in row
! 
!.... variable declarations
      !implicit none
      CHARACTER(2), INTENT(IN) :: CINVTYPE
      COMPLEX*8, INTENT(IN) :: WAVE(NDOF)
      REAL*8, INTENT(IN) :: SHL(3,MEN,MINTX,*), ECOEFF(MNPG,*),
     ;                      XIGLL(MINTX,*), ETAGLL(MINTZ,*),
     ;                      XLOCS(NNPG), ZLOCS(NNPG), DENS(NNPG),
     ;                      OMEGA, PY  
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), IENG(MGNOD,*), 
     ;                         MYGELEM(NELEMG), 
     ;                         MNZROW,MDIM,MINTX,MINTZ,MGNOD,MEN,MNPG,
     ;                         NNPG,NDOF, NINTX,NINTZ,NGNOD,NEN,NELEMG,
     ;                         NVINV, INPG 
      LOGICAL*4, INTENT(IN) :: LISISO 
      COMPLEX*8, INTENT(OUT) :: FORCEV(MNZROW,*)
      INTEGER*4, INTENT(OUT) :: MYROWS(MNZROW), NNZROW, IERR  
!.... local variables
      COMPLEX*16, ALLOCATABLE :: DKDA(:,:), DKDB(:,:)
      COMPLEX*8, ALLOCATABLE :: CARG(:) 
      REAL*8, ALLOCATABLE :: SHG(:,:,:,:), DAB(:,:,:), DET(:,:), 
     ;                       DMDA(:,:), DMDB(:,:)
      COMPLEX*8 CZERO 
      INTEGER*4 NEE,NDIM, IELEMG,IELEM, IAE,I,IBE,J, IDOF,JDOF,INZROW, 
     ;          JNZROW,IPE,IQE,IVINV  
      LOGICAL*4 LISP, LISS 
      INTEGER*4 IBSECT
      PARAMETER(CZERO = CMPLX(0.,0.)) 
      PARAMETER(NDIM = 3) 
!
!----------------------------------------------------------------------!
!  
!.... loop on elements in gradient attached to anchor node
      IERR = 0 
      INZROW = 0 
      MYROWS(1:MNZROW) = 0 
      DO 1 IELEMG=1,NELEMG !elements attached to anchor node 
         IELEM = MYGELEM(IELEMG)
         DO 2 IAE=1,NEN !element nodes
            DO 3 I=1,NDIM !spatial components
               IDOF = LM(I,IAE,IELEM) !extract DOF
               IF (IDOF.GT.0) THEN !dof
                  DO 4 JNZROW=1,INZROW !prevent repeats
                     IF (MYROWS(JNZROW).EQ.IDOF) GOTO 40 !repeat
    4             CONTINUE 
                  INZROW = INZROW + 1
                  IF (INZROW.GT.MNZROW) THEN
                     WRITE(*,*) 'vforce: Error myrows is too small' 
                     IERR = 1
                     RETURN
                  ENDIF
                  MYROWS(INZROW) = IDOF 
   40             CONTINUE !repeat, break ahead
               ENDIF
    3       CONTINUE !Loop on components 
    2    CONTINUE !loop on element nodes
    1 CONTINUE !Loop on elements
      NNZROW = INZROW 
      CALL ISHELL1(NNZROW,MYROWS) !sort so we can find it
!
!.... set space in derivative matrices and element shape fns 
      NEE = NEN*NDIM
      IF (LISISO) THEN 
         LISP = .FALSE.
         LISS = .FALSE.
         IF (CINVTYPE == 'PP' .OR. CINVTYPE == 'pp') THEN 
            LISP = .TRUE.
         ELSEIF (CINVTYPE == 'SS' .OR. CINVTYPE == 'ss') THEN 
            LISS = .TRUE.
         ELSE
            LISP = .TRUE.
            LISS = .TRUE.
         ENDIF
         ALLOCATE(DKDA(NEE,NEE))
         ALLOCATE(DKDB(NEE,NEE)) 
         ALLOCATE(DMDA(NEE,NEE)) 
         ALLOCATE(DMDB(NEE,NEE))
         ALLOCATE(DAB(MINTX,MINTZ,2)) 
      ELSE
         WRITE(*,*) 'vforce: Anisotropy not yet programmed!'
         IERR = 1
         RETURN
      ENDIF
      ALLOCATE(SHG(3,MEN,MINTX,MINTZ))
      ALLOCATE(DET(MINTX,MINTZ)) 
      ALLOCATE(CARG(NVINV))
  
!
!.... null out force vector 
      DO 5 INZROW=1,NNZROW 
         FORCEV(INZROW,1:NVINV) = CZERO  
    5 CONTINUE 
!
!.... loop on partitions elements in gradient
      DO 10 IELEMG=1,NELEMG  
         IELEM = MYGELEM(IELEMG)
         IF (LISISO) THEN !P and S gradients 
!
!.......... evaluate model derivatives at integraiton points
            CALL DJACAB(MEN,MINTX,MINTZ,MNPG, LISISO,        
     ;                  NINTX,NINTZ,NEN,NGNOD,NNPG, INPG, 
     ;                  IENG(1:NGNOD,IELEM), 
     ;                  XIGLL(1:NINTX,1),ETAGLL(1:NINTZ,1),
     ;                  XLOCS,ZLOCS,DENS,ECOEFF,SHL, 
     ;                  DET,DAB,SHG,IERR) 
            IF (IERR.NE.0) THEN
               WRITE(*,*) 'vforce: Error calling djacab'
               RETURN
            ENDIF
!
!.......... numerically integrate element matrices
            CALL DABMKISOB25(MEN,NEE,MINTX,MINTZ, NEN,NINTX,NINTZ, 
     ;                       LISP,LISS,OMEGA,PY, 
     ;                       XIGLL(1:NINTX,2),ETAGLL(1:NINTZ,2),
     ;                       DET,DAB,SHG, DMDA,DMDB,DKDA,DKDB)  
         ELSE  
            WRITE(*,*) 'vforce: No anistropy yet!'
            IERR = 1
            RETURN 
         ENDIF
!
!....... matrix vector multiply of dS/dmi u
         IPE = 0
         DO 11 IAE=1,NEN !loop on element nodes
            DO 12 I=1,NDIM !loop on components
               IPE = IPE + 1
               IDOF = LM(I,IAE,IELEM) !extra row
               IF (IDOF.LE.0) GOTO 120 !not a dof
               INZROW = IBSECT(NNZROW,IDOF,MYROWS) !non-zero row index 
               IF (INZROW.LE.0 .OR. INZROW.GT.NDOF) THEN
                  WRITE(*,*) 'vforce25: Error calling IBSECT'
                  IERR = 1
               ENDIF 
               IF (MYROWS(INZROW).NE.IDOF) THEN
                  WRITE(*,*) 'vforce25: Error locating DOF!'
                  IERR = 1
               ENDIF
              
!
!............. loop on columns associated with row
               IQE = 0
               DO 13 IBE=1,NEN  
                  DO 14 J=1,NDIM
                     IQE = IQE + 1
                     JDOF = LM(J,IBE,IELEM)
                     IF (JDOF.LE.0) GOTO 140
                     CARG(1:NVINV) = CMPLX(0.0,0.0)
                     IVINV = 0
                     IF (LISISO) THEN
                        IF (LISP) THEN !dS/dalpha
                           IVINV = IVINV + 1
                           CARG(IVINV) =CMPLX(DCMPLX(DMDA(IPE,IQE),0.D0)
     ;                                             + DKDA(IPE,IQE))
                        ENDIF
                        IF (LISS) THEN !dS/dbeta
                           IVINV = IVINV + 1
                           CARG(IVINV) =CMPLX(DCMPLX(DMDB(IPE,IQE),0.D0)
     ;                                             + DKDB(IPE,IQE))
                        ENDIF 
                     ELSE
                        WRITE(*,*) 'vforce25: Error no ansitropy!'
                        IERR = 1
                        RETURN
                     ENDIF
                     DO 15 IVINV=1,NVINV
                        FORCEV(INZROW,IVINV) = FORCEV(INZROW,IVINV) 
     ;                                       + CARG(IVINV)*WAVE(JDOF)
   15                CONTINUE !loop on variables in inversion
  140                CONTINUE !break ahead, not a DOF 
   14             CONTINUE !loop on components
   13          CONTINUE !loop element nodes

  120          CONTINUE !break ahed, not a DOF
   12       CONTINUE !loop on components
   11    CONTINUE !loop on element nodes
   10 CONTINUE !loop on elements
!
!.... free memory 
      IF (ALLOCATED(DAB))  DEALLOCATE(DAB)
      IF (ALLOCATED(DMDA)) DEALLOCATE(DMDA)
      IF (ALLOCATED(DMDB)) DEALLOCATE(DMDB)
      IF (ALLOCATED(DKDA)) DEALLOCATE(DKDA)
      IF (ALLOCATED(DKDB)) DEALLOCATE(DKDB) 
      IF (ALLOCATED(CARG)) DEALLOCATE(CARG) 
      IF (ALLOCATED(SHG))  DEALLOCATE(SHG)
      IF (ALLOCATED(DET))  DEALLOCATE(DET) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE FD_EMATS_25(MEE,NEN, SMALL,ME,ME1, KE,KE1, 
     ;                       DMDP_FD,DKDP_FD)
!
!     Finite differences the element matrices

      IMPLICIT NONE
      COMPLEX*16, INTENT(IN) :: KE1(MEE,*), KE(MEE,*)  
      REAL*8, INTENT(IN) :: ME1(MEE,*), ME(MEE,*), SMALL
      INTEGER*4, INTENT(IN) :: MEE,NEN 
      COMPLEX*16, INTENT(OUT) :: DKDP_FD(MEE,*)
      REAL*8, INTENT(OUT) :: DMDP_FD(MEE,*) 
!.... local variables 
      REAL*8 KDIF11,KDIF12,KDIF13, KDIF21,KDIF22,KDIF23, 
     ;       KDIF31,KDIF32,KDIF33, MDIF11,MDIF22,MDIF33
      INTEGER*4 NDIM, IAE, IBE, IPE, IQE 
      PARAMETER(NDIM = 3) 
!
!----------------------------------------------------------------------!
!
      DO 1 IAE=1,NEN
         DO 2 IBE=1,NEN
            IPE = (IAE - 1)*NDIM + 1
            IQE = (IBE - 1)*NDIM + 1
            !mass
            MDIF11 = ME1(IPE  ,IQE  ) - ME(IPE  ,IQE  )
            MDIF22 = ME1(IPE+1,IQE+1) - ME(IPE+1,IQE+1)
            MDIF33 = ME1(IPE+2,IQE+2) - ME(IPE+2,IQE+2)
            DMDP_FD(IPE  ,IQE  ) = MDIF11/SMALL
            DMDP_FD(IPE  ,IQE+1) = 0.D0
            DMDP_FD(IPE  ,IQE+2) = 0.D0
            DMDP_FD(IPE+1,IQE  ) = 0.D0
            DMDP_FD(IPE+1,IQE+1) = MDIF22/SMALL
            DMDP_FD(IPE+1,IQE+2) = 0.D0
            DMDP_FD(IPE+2,IQE  ) = 0.D0
            DMDP_FD(IPE+2,IQE+1) = 0.D0
            DMDP_FD(IPE+2,IQE+2) = MDIF33/SMALL
            !stiffness
            KDIF11 = DREAL(KE1(IPE  ,IQE  ) - KE(IPE  ,IQE  ))
            KDIF12 = DIMAG(KE1(IPE  ,IQE+1) - KE(IPE  ,IQE+1))
            KDIF13 = DREAL(KE1(IPE  ,IQE+2) - KE(IPE  ,IQE+2))
            KDIF21 = DIMAG(KE1(IPE+1,IQE  ) - KE(IPE+1,IQE  ))
            KDIF22 = DREAL(KE1(IPE+1,IQE+1) - KE(IPE+1,IQE+1))
            KDIF23 = DIMAG(KE1(IPE+1,IQE+2) - KE(IPE+1,IQE+2))
            KDIF31 = DREAL(KE1(IPE+2,IQE  ) - KE(IPE+2,IQE  ))
            KDIF32 = DIMAG(KE1(IPE+2,IQE+1) - KE(IPE+2,IQE+1))
            KDIF33 = DREAL(KE1(IPE+2,IQE+2) - KE(IPE+2,IQE+2))
            DKDP_FD(IPE  ,IQE  ) = DCMPLX(KDIF11/SMALL,0.D0)
            DKDP_FD(IPE  ,IQE+1) = DCMPLX(0.D0,KDIF12/SMALL)
            DKDP_FD(IPE  ,IQE+2) = DCMPLX(KDIF13/SMALL,0.D0)
            DKDP_FD(IPE+1,IQE  ) = DCMPLX(0.D0,KDIF21/SMALL)
            DKDP_FD(IPE+1,IQE+1) = DCMPLX(KDIF22/SMALL,0.D0)
            DKDP_FD(IPE+1,IQE+2) = DCMPLX(0.D0,KDIF23/SMALL)
            DKDP_FD(IPE+2,IQE  ) = DCMPLX(KDIF31/SMALL,0.D0)
            DKDP_FD(IPE+2,IQE+1) = DCMPLX(0.D0,KDIF32/SMALL)
            DKDP_FD(IPE+2,IQE+2) = DCMPLX(KDIF33/SMALL,0.D0)
    2    CONTINUE !Loop on element nodes 
    1 CONTINUE !Loop on element nodes
      RETURN
      END

