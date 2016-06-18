!----------------------------------------------------------------------!
!     Updated utilities for handling the Bielak conditions on          !
!     unstructured grids.  It turns out we can be a bit more clever    !
!     with how we generate a mesh and can really simplify things       !
!     for the solver                                                   !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_1D(CSIDE,MYID, NL,POFF,Z1D,VP1D,VS1D,RH1D)
!
!     Writes the 1D model for a specific side
      CHARACTER(1), INTENT(IN) :: CSIDE
      REAL*8, INTENT(IN) :: Z1D(NL), VP1D(NL), VS1D(NL), RH1D(NL),POFF 
      INTEGER*4, INTENT(IN) :: NL
      CHARACTER(80) FLNAME
      CHARACTER(4) CMYID
!
!----------------------------------------------------------------------!
!
      FLNAME(1:80) = ' ' 
      CMYID(1:4) = ' ' 
      WRITE(CMYID,'(I4)') MYID
      CMYID = ADJUSTL(CMYID)
      IF (CSIDE.EQ.'L' .OR. CSIDE.EQ.'l') THEN
         FLNAME = 'left_model-' //TRIM(CMYID)//'.txt'
      ELSE
         FLNAME = 'right_model-'//TRIM(CMYID)//'.txt'
      ENDIF
      FLNAME = ADJUSTL(FLNAME) 
      IUNIT = 44 + MYID
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME))
      WRITE(IUNIT,*) NL
      WRITE(IUNIT,*) POFF
      DO 1 I=1,NL
         WRITE(IUNIT,*) Z1D(I),VP1D(I),VS1D(I),RH1D(I)
    1 CONTINUE
      CLOSE(IUNIT) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE READ1D_HD(CSIDE,MYID, NL,IERR)
!     Reads the header off the appropriate model side and returns 
!     the number of layers
      CHARACTER(1), INTENT(IN) :: CSIDE
      INTEGER*4, INTENT(IN) :: MYID 
      INTEGER*4, INTENT(OUT) :: NL,IERR
      CHARACTER(80) FLNAME
      CHARACTER(4) CMYID
      LOGICAL*4 LEX
!
!----------------------------------------------------------------------!
!
      IERR = 0 
      FLNAME(1:80) = ' '
      CMYID(1:4) = ' ' 
      WRITE(CMYID,'(I4)') MYID
      CMYID = ADJUSTL(CMYID)
      IF (CSIDE.EQ.'L' .OR. CSIDE.EQ.'l') THEN
         FLNAME = 'left_model-' //TRIM(CMYID)//'.txt'
      ELSE
         FLNAME = 'right_model-'//TRIM(CMYID)//'.txt'
      ENDIF
      FLNAME = ADJUSTL(FLNAME) 
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'read1d_hd: 1D model file does not exist!'
         IERR = 1 
         RETURN
      ENDIF
      IUNIT = 44 + MYID
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME))
      READ(IUNIT,*,END=60) NL
      CLOSE(IUNIT) 
      RETURN
   60 CONTINUE
      IERR = 1
      WRITE(*,*) 'read1d_hd: Premature end of 1D input file'
      CLOSE(IUNIT) 
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE READ1D(CSIDE,MYID, NL,POFF,Z1D,VP1D,VS1D,RH1D,IERR)
!     Reads the header off the appropriate model side and returns 
!     the number of layers
      CHARACTER(1), INTENT(IN) :: CSIDE
      INTEGER*4, INTENT(IN) :: MYID,NL
      REAL*8, INTENT(OUT) :: Z1D(NL),VP1D(NL),VS1D(NL),RH1D(NL),POFF  
      INTEGER*4, INTENT(OUT) :: IERR
      CHARACTER(80) FLNAME
      CHARACTER(4) CMYID
!
!----------------------------------------------------------------------!
!
      IERR = 0 
      FLNAME(1:80) = ' ' 
      CMYID(1:4) = ' ' 
      WRITE(CMYID,'(I4)') MYID
      CMYID = ADJUSTL(CMYID)
      IF (CSIDE.EQ.'L' .OR. CSIDE.EQ.'l') THEN
         FLNAME = 'left_model-' //TRIM(CMYID)//'.txt'
      ELSE
         FLNAME = 'right_model-'//TRIM(CMYID)//'.txt'
      ENDIF
      FLNAME = ADJUSTL(FLNAME) 
      IUNIT = 44 + MYID
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME))
      READ(IUNIT,*,END=60) !NL
      READ(IUNIT,*,END=60) POFF
      DO 1 I=1,NL
         READ(IUNIT,*,END=60) Z1D(I),VP1D(I),VS1D(I),RH1D(I)
    1 CONTINUE
      CLOSE(IUNIT) 
      RETURN
   60 CONTINUE
      IERR = 1 
      WRITE(*,*) 'read1d: Premature end of 1D input file'
      CLOSE(IUNIT) 
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GXZPTS_BLK2(MDIM,MEN,MGNOD, NDOF,NNPG,NELEM,NLXI,NLETA,
     ;                       NWORK, NGNOD,NDIM, CDOMAIN,CNP, IENG,LM, 
     ;                       XIPTS,ETAPTS, XLOCS,ZLOCS,  
     ;                       NNPE,XLOCSE,ZLOCSE,IERR)

      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM), CNP(NDOF)
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG), XIPTS(NLXI),
     ;                      ETAPTS(NLETA)
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), IENG(MGNOD,*), 
     ;           MDIM,MEN,MGNOD, NDOF,NNPG,NELEM,NLXI,NLETA,NWORK, 
     ;           NGNOD,NDIM 
      REAL*8, INTENT(OUT) :: XLOCSE(NWORK), ZLOCSE(NWORK)
      INTEGER*4, INTENT(OUT) :: NNPE,IERR
!.... local variables
      REAL*8, ALLOCATABLE :: SF(:)
      REAL*8 X,Z,XI,ETA, TOL
      LOGICAL*4 LWORK
      PARAMETER(TOL = 1.D-7)  
!
!----------------------------------------------------------------------!
!
!.... initialize 
      IERR = 0
      INPE = 0
      ALLOCATE(SF(NGNOD))
!
!.... loop on domain
      DO 1 IELEM=1,NELEM
!
!....... element belongs to bielak domain
         IF (CDOMAIN(IELEM).EQ.'E') THEN
            DO 2 ILETA=1,NLETA
               ETA = ETAPTS(ILETA)
               DO 3 ILXI=1,NLXI 
                  IAE = (ILETA - 1)*NLXI + ILXI
                  XI = XIPTS(ILXI)
                  LWORK = .FALSE.
                  IF (MINVAL(LM(1:NDIM,IAE,IELEM)).LE.0) GOTO 30
                  DO 4 I=1,NDIM
                     IDOF = LM(I,IAE,IELEM)
                     IF (IDOF.GT.0) THEN
                        IF (CNP(IDOF).EQ.'E' .OR. 
     ;                      CNP(IDOF).EQ.'B') THEN
                           LWORK = .TRUE.
                           GOTO 40
                        ENDIF
                     ENDIF
    4             CONTINUE 
   40             CONTINUE
                  IF (LWORK) THEN
                     X = 0.D0 
                     Z = 0.D0
                     CALL CSF2DND(NGNOD,XI,ETA, SF,IERR)
                     IF (IERR.NE.0) THEN
                        WRITE(*,*) 'gxzpts_blk2: Error in csf2dnd'
                        GOTO 730
                     ENDIF
                     DO 5 IA=1,NGNOD
                        ILOC = IENG(IA,IELEM)
                        X = X + SF(IA)*XLOCS(ILOC)
                        Z = Z + SF(IA)*ZLOCS(ILOC) 
    5                CONTINUE
!
!................... check if the node exists
                     DO 6 JNPE=1,INPE
                        IF (DABS(XLOCSE(JNPE) - X).LT.TOL .AND.
     ;                      DABS(ZLOCSE(JNPE) - Z).LT.TOL) GOTO 60
    6                CONTINUE
                     INPE = INPE + 1
                     XLOCSE(INPE) = X
                     ZLOCSE(INPE) = Z
   60                CONTINUE !break ahead, elements
                  ENDIF !end check on working
   30             CONTINUE !break ahead for new node
    3          CONTINUE
    2       CONTINUE
         ENDIF
    1 CONTINUE
  730 CONTINUE
      DEALLOCATE(SF) 
      IF (IERR.GT.0) RETURN
!.... sort on z 
      WRITE(*,*) 'gxzpts_blk2: Sorting z locations...'
      NNPE = INPE
      CALL DHPSRT2(NNPE,ZLOCSE(1:NNPE),XLOCSE(1:NNPE))
!.... sort on x 
      WRITE(*,*) 'gxzpts_blk2: Sorting x locations...'
      I2 = 0
      DO 21 JNPE=1,NNPE
         I1 = I2 + 1
         DO 22 KNPE=I1,NNPE
            IF (KNPE.LT.NNPE) THEN
               IF (DABS(ZLOCSE(KNPE) - ZLOCSE(KNPE+1)).GT.TOL) THEN
                  I2 = KNPE
                  NVAR = I2 - I1 + 1
                  IF (NVAR.GT.1) CALL SHELL1(NVAR, XLOCSE(I1:I2))
                  GOTO 110
               ENDIF
            ELSE
               ! here the final point, nnpe is obviously sorted
               IF (DABS(ZLOCSE(NNPE) - ZLOCSE(NNPE-1)).GT.TOL) THEN
                  I2 = NNPE - 1
                  NVAR = I2 - I1 + 1
                  IF (NVAR.GT.1) CALL SHELL1(NVAR, XLOCSE(I1:I2))
               ELSE !here the last point is included in the sort
                  I2 = NNPE
                  NVAR = I2 - I1 + 1
                  IF (NVAR.GT.1) CALL SHELL1(NVAR, XLOCSE(I1:I2))
               ENDIF
               GOTO 100
            ENDIF
   22    CONTINUE !loop on points in bielak domain
  110    CONTINUE !Out of search loop
   21 CONTINUE !loop on points in bielak domain
  100 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CDOF2XZE2(MDIM,MEN,MGNOD, NDOF,NNPG,NELEM,NLXI,NLETA,
     ;                     NNPE, NGNOD,NDIM, CDOMAIN,CNP, IENG,LM, 
     ;                     XIPTS,ETAPTS, XLOCS,ZLOCS, XLOCSE,ZLOCSE, 
     ;                     IDOFSE,IERR)

      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM), CNP(NDOF)
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG), XLOCSE(NNPE),
     ;                      ZLOCSE(NNPE), XIPTS(NLXI), ETAPTS(NLETA)
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), IENG(MGNOD,*), 
     ;           MDIM,MEN,MGNOD, NDOF,NNPG,NELEM,NLXI,NLETA,NNPE, 
     ;           NGNOD,NDIM 
      INTEGER*4, INTENT(OUT) :: IDOFSE(MDIM,*), IERR
!.... local variables
      REAL*8, ALLOCATABLE :: SF(:)
      REAL*8 X,Z,XI,ETA, TOL 
      LOGICAL*4 LWORK
      PARAMETER(TOL = 1.D-7)  
      PARAMETER(LNULL =-5) 
!
!----------------------------------------------------------------------!
!
!.... initialize 
      IERR = 0 
      INPE = 0 
      ALLOCATE(SF(NGNOD))
      IDOFSE(1:NDIM,1:NNPE) = LNULL
!
!.... loop on domain
      DO 1 IELEM=1,NELEM
!
!....... element belongs to bielak domain
         IF (CDOMAIN(IELEM).EQ.'E') THEN
            DO 2 ILETA=1,NLETA
               ETA = ETAPTS(ILETA)
               DO 3 ILXI=1,NLXI
                  IAE = (ILETA - 1)*NLXI + ILXI
                  XI = XIPTS(ILXI)
                  LWORK = .FALSE.
                  IF (MINVAL(LM(1:NDIM,IAE,IELEM)).LE.0) GOTO 30
                  DO 4 I=1,NDIM
                     IDOF = LM(I,IAE,IELEM)
                     IF (IDOF.GT.0) THEN
                        IF (CNP(IDOF).EQ.'E' .OR.
     ;                      CNP(IDOF).EQ.'B') THEN
                           LWORK = .TRUE.
                           GOTO 40
                        ENDIF
                     ENDIF
    4             CONTINUE
   40             CONTINUE
                  IF (LWORK) THEN
                     X = 0.D0
                     Z = 0.D0
                     CALL CSF2DND(NGNOD,XI,ETA, SF,IERR)
                     IF (IERR.NE.0) THEN
                        WRITE(*,*) 'gxzpts_blk2: Error in csf2dnd'
                        GOTO 730
                     ENDIF
                     DO 5 IA=1,NGNOD
                        ILOC = IENG(IA,IELEM)
                        X = X + SF(IA)*XLOCS(ILOC)
                        Z = Z + SF(IA)*ZLOCS(ILOC)
    5                CONTINUE
!
!................... check if the node exists
                     JBEG = IBSECT8(NNPE,1,TOL,Z,ZLOCSE)
                     IF (JBEG.GT.0) THEN
                        DO 6 INPE=JBEG,NNPE
                           IF (DABS(XLOCSE(INPE) - X).LT.TOL .AND.
     ;                         DABS(ZLOCSE(INPE) - Z).LT.TOL) THEN 
                              DO 7 I=1,NDIM
                                 IDOFSE(I,INPE) = LM(I,IAE,IELEM)
    7                         CONTINUE
                              GOTO 60
                           ENDIF
    6                   CONTINUE
                     ENDIF 
                     WRITE(*,*) 'cdof2xze2: Warning brute force!'
                     DO 8 INPE=1,NNPE
                        IF (DABS(XLOCSE(INPE) - X).LT.TOL .AND.
     ;                      DABS(ZLOCSE(INPE) - Z).LT.TOL) THEN
                           DO 9 I=1,NDIM
                              IDOFSE(I,INPE) = LM(I,IAE,IELEM)
    9                      CONTINUE
                           GOTO 60
                        ENDIF 
    8                CONTINUE
                     WRITE(*,*) 'cdof2xze2: Cant locate node!'
                     IERR = 1
                     GOTO 730
   60                CONTINUE !break ahead, elements
                  ENDIF !end check on working
   30             CONTINUE !break ahead for new node
    3          CONTINUE
    2       CONTINUE
         ENDIF
    1 CONTINUE
  730 CONTINUE
      DEALLOCATE(SF)
      IF (IERR.NE.0) RETURN
!
!.... quality control
      DO 21 INPE=1,NNPE
         IF (MINVAL(IDOFSE(1:NDIM,INPE)).EQ.LNULL) THEN
            WRITE(*,*) 'cdof2xze2: Failed to initialize inpe:',INPE
            IERR = 1
         ENDIF
   21 CONTINUE
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GXZPTS_BLK(MGNOD,NELEM,NNPG,NWORK,NLXI,NLETA, NGNOD,
     ;                  NELEME, CNNPG,CDOMAIN,IENG, XIPTS,ETAPTS, 
     ;                  XLOCS,ZLOCS, NNPE,XLOCSE,ZLOCSE,IERR)
!
!     Generates the x and z points in the bielak domain
      IMPLICIT REAL*8 (A-H,O-Z) 
      CHARACTER(2), INTENT(IN) :: CNNPG(NNPG) 
      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM)
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG), XIPTS(NLXI), 
     ;        ETAPTS(NLETA) 
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MGNOD, NELEM, NNPG, NWORK,
     ;           NLXI, NLETA, NGNOD
      REAL*8, INTENT(OUT) :: XLOCSE(NWORK), ZLOCSE(NWORK)
      INTEGER*4, INTENT(OUT) :: NNPE, IERR
!.... local variables
      REAL*8, ALLOCATABLE :: XS(:), ZS(:), SF(:) 
      REAL*8 V1(3), V2(3) 
      PARAMETER(TOL = 1.D-10) 
!
!----------------------------------------------------------------------!
! 
!.... initialize 
      IERR = 0 
      INPE = 0
      INPS = 0
      NTEMP = 4*MAX0(NLXI,NLETA)*NELEME
      ALLOCATE(SF(NGNOD))
      ALLOCATE(XS(NTEMP))
      ALLOCATE(ZS(NTEMP)) 
!
!.... calculate (x,z) locations of all points within Bielak domain 
      DO 1 IELEM=1,NELEM
         IF (CDOMAIN(IELEM).EQ.'E') THEN !check we are bielak domain
            IA1 = 0 
            IA2 = 0
            KEEP = 0 
            DO 2 IA=1,NGNOD
               ILOC = IENG(IA,IELEM)
               IF (CNNPG(ILOC).EQ.'BA') THEN 
                  KEEP = KEEP + 1
                  IA1 = ILOC
                  DO 3 JA=1,NGNOD
                     JLOC = IENG(JA,IELEM)
                     IF (CNNPG(JLOC).EQ.'FS' .OR. 
     ;                   CNNPG(JLOC).EQ.'BA') THEN
                        IF (JLOC.NE.ILOC) THEN
                           KEEP = KEEP + 1
                           IA2 = JLOC
                           GOTO 45
                        ENDIF !end check on exists 
                     ENDIF !check on 'BA' 'FS' condition
    3             CONTINUE !second loop on anchor noes
                  GOTO 45
               ENDIF !first check on 'BA' condition
    2       CONTINUE !loop on anchor nodes 
   45       CONTINUE !break ahead
            V1(1:3) = 0.D0
            IF (KEEP.EQ.2) THEN
               V1(1) = XLOCS(IA2) - XLOCS(IA1)
               V1(2) = 0.D0
               V1(3) = ZLOCS(IA2) - ZLOCS(IA1)
            ENDIF
!
!.......... loop on element nodes 
            DO 6 ILETA=1,NLETA 
               ETA = ETAPTS(ILETA)
               DO 7 ILXI=1,NLXI
                  XI = XIPTS(ILXI) 
                  CALL CSF2DND(NGNOD,XI,ETA, SF, IERR)
                  IF (IERR.NE.0) THEN
                     WRITE(*,*) 'gxzpts_blk: Error calling csf2dnd'
                     IERR = 1
                     GOTO 70
                  ENDIF
                  X = 0.D0
                  Z = 0.D0
                  DO 8 IA=1,NGNOD
                     ILOC = IENG(IA,IELEM)
                     X = X + SF(IA)*XLOCS(ILOC)
                     Z = Z + SF(IA)*ZLOCS(ILOC)
    8             CONTINUE !loop on anchor nodes
!
!................ just one anchor node attached to boundary
                  IF (KEEP.EQ.1) THEN
                     IF (DABS(X - XLOCS(IA1)).LT.TOL .AND. 
     ;                   DABS(Z - ZLOCS(IA1)).LT.TOL) GOTO 400 
                  ELSEIF (KEEP.EQ.2) THEN !along line?
                     V2(1) = X - XLOCS(IA1) !line cant be parallel
                     V2(2) = 0.D0
                     V2(3) = Z - ZLOCS(IA1)
                     DET = V1(1)*V2(3) - V1(3)*V2(1)
                     IF (DABS(DET).LT.TOL) GOTO 400 
                  ENDIF
!
!................ don't locaitons on element sides
                  IF (ILXI .EQ.1 .OR. ILXI .EQ.NLXI   .OR. 
     ;                ILETA.EQ.1 .OR. ILETA.EQ.NLETA) THEN
                     IF (INPS.EQ.0) THEN
                        INPS = INPS + 1
                        XS(INPS) = X
                        ZS(INPS) = Z 
                        INPE = INPE + 1
                        XLOCSE(INPE) = X
                        ZLOCSE(INPE) = Z
                     ELSE
                        DO 9 JNPS=1,INPS
                           IF (DABS(XS(JNPS) - X).LT.TOL .AND. 
     ;                         DABS(ZS(JNPS) - Z).LT.TOL) GOTO 50 
    9                   CONTINUE
                        INPS = INPS + 1
                        XS(INPS) = X
                        ZS(INPS) = Z
                        INPE = INPE + 1
                        IF (INPE.GT.NWORK) THEN
                           WRITE(*,*) 'gxzpts_blk: nwork too small!' 
                           IERR = 2
                           GOTO 70 
                        ENDIF 
                        XLOCSE(INPE) = X
                        ZLOCSE(INPE) = Z
   50                   CONTINUE !break ahead; point exists
                     ENDIF !end check on inps
                  ELSE !element interior, always new
                     INPE = INPE + 1
                     IF (INPE.GT.NWORK) THEN
                        WRITE(*,*) 'gxzpts_blk: nwork too small!'
                        IERR = 2
                        GOTO 70
                     ENDIF
                     XLOCSE(INPE) = X
                     ZLOCSE(INPE) = Z
                  ENDIF !end check on interior/boundary
  400             CONTINUE !break ahead for new xi point
    7          CONTINUE !loop on xi  
    6       CONTINUE
         ENDIF
    1 CONTINUE !loop on elements
   70 CONTINUE !break ahead for errros
      DEALLOCATE(SF)
      DEALLOCATE(XS)
      DEALLOCATE(ZS)
      IF (IERR.GT.0) RETURN 
!.... sort on z 
      WRITE(*,*) 'gxzpts_blk: Sorting z locations...'
      NNPE = INPE
      CALL DHPSRT2(NNPE,ZLOCSE(1:NNPE),XLOCSE(1:NNPE))
!.... sort on x 
      WRITE(*,*) 'gxzpts_blk: Sorting x locations...'
      I2 = 0 
      DO 21 JNPE=1,NNPE 
         I1 = I2 + 1 
         DO 22 KNPE=I1,NNPE
            IF (KNPE.LT.NNPE) THEN 
               IF (DABS(ZLOCSE(KNPE) - ZLOCSE(KNPE+1)).GT.TOL) THEN
                  I2 = KNPE  
                  NVAR = I2 - I1 + 1
                  IF (NVAR.GT.1) CALL SHELL1(NVAR, XLOCSE(I1:I2)) 
                  GOTO 110
               ENDIF
            ELSE
               ! here the final point, nnpe is obviously sorted
               IF (DABS(ZLOCSE(NNPE) - ZLOCSE(NNPE-1)).GT.TOL) THEN
                  I2 = NNPE - 1 
                  NVAR = I2 - I1 + 1 
                  IF (NVAR.GT.1) CALL SHELL1(NVAR, XLOCSE(I1:I2))
               ELSE !here the last point is included in the sort
                  I2 = NNPE 
                  NVAR = I2 - I1 + 1 
                  IF (NVAR.GT.1) CALL SHELL1(NVAR, XLOCSE(I1:I2)) 
               ENDIF 
               GOTO 100
            ENDIF
   22    CONTINUE !loop on points in bielak domain
  110    CONTINUE !Out of search loop
   21 CONTINUE !loop on points in bielak domain
  100 CONTINUE
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !

      SUBROUTINE CDOF2XZE(MDIM,MEN,MGNOD, NELEM,NNPG,NLXI,NLETA,NNPE, 
     ;                    NGNOD,NDIM, CNNPG,CDOMAIN, IENG,LM, 
     ;                    XIPTS,ETAPTS, XLOCS,ZLOCS,XLOCSE,ZLOCSE, 
     ;                    IDOFSE, IERR)
 
! 
!     Connects the DOF numbers to (xlocse,zlocse) locations read from 
!     disk in directory GRNS 
!
!.... variable declarations
      IMPLICIT REAL*8 (A-H,O-Z) 
      CHARACTER(2), INTENT(IN) :: CNNPG(NNPG)  
      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM)
      REAL*8, INTENT(IN) :: XLOCS(NNPG),ZLOCS(NNPG), 
     ;        XLOCSE(NNPE),ZLOCSE(NNPE), XIPTS(NLXI),ETAPTS(NLETA)  
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), IENG(MGNOD,*), 
     ;           MDIM,MEN,MGNOD, NELEM,NNPG,NLXI,NLETA,NNPE, 
     ;           NGNOD,NDIM 
      INTEGER*4, INTENT(OUT) :: IDOFSE(MDIM,*), IERR 
!.... local variables 
      REAL*8, ALLOCATABLE :: SF(:)
      REAL*8 V1(3), V2(3) 
      REAL*8 TOL/1.D-10/
      PARAMETER(LNULL =-5) 
      INTEGER*4 IBSECT8
! 
!----------------------------------------------------------------------!
!
!.... null out IDOFSE
      IERR = 0
      ALLOCATE(SF(NGNOD))
      DO 1 INPE=1,NNPE
         IDOFSE(1:NDIM,INPE) = LNULL 
    1 CONTINUE
!
!.... loop on points
      DO 2 IELEM=1,NELEM
         IF (CDOMAIN(IELEM).EQ.'E') THEN
            IA1 = 0 
            IA2 = 0
            KEEP = 0
            DO 4 IA=1,NGNOD
               ILOC = IENG(IA,IELEM)
               IF (CNNPG(ILOC).EQ.'BA') THEN
                  KEEP = KEEP + 1
                  IA1 = ILOC
                  DO 5 JA=1,NGNOD
                     JLOC = IENG(JA,IELEM)
                     IF (CNNPG(JLOC).EQ.'FS' .OR.
     ;                   CNNPG(JLOC).EQ.'BA') THEN
                        IF (JLOC.NE.ILOC) THEN
                           KEEP = KEEP + 1
                           IA2 = JLOC
                           GOTO 45
                        ENDIF !end check on exists 
                     ENDIF !check on 'BA' 'FS' condition
    5             CONTINUE !second loop on anchor noes
                  GOTO 45
               ENDIF !first check on 'BA' condition
    4       CONTINUE !loop on anchor nodes 
   45       CONTINUE !break ahead
!
!.......... set determinant
            V1(1:3) = 0.D0
            IF (KEEP.EQ.2) THEN
               V1(1) = XLOCS(IA2) - XLOCS(IA1)
               V1(2) = 0.D0
               V1(3) = ZLOCS(IA2) - ZLOCS(IA1)
            ENDIF
!
!.......... loop on element nodes 
            DO 6 ILETA=1,NLETA
               ETA = ETAPTS(ILETA)
               DO 7 ILXI=1,NLXI
                  IAE = (ILETA - 1)*NLXI + ILXI
                  XI = XIPTS(ILXI)
                  CALL CSF2DND(NGNOD,XI,ETA, SF, IERR)
                  IF (IERR.NE.0) THEN
                     WRITE(*,*) 'gxzpts_blk: Error calling csf2dnd'
                     IERR = 1
                     GOTO 70 
                  ENDIF
                  X = 0.D0
                  Z = 0.D0
                  DO 8 IA=1,NGNOD
                     ILOC = IENG(IA,IELEM)
                     X = X + SF(IA)*XLOCS(ILOC)
                     Z = Z + SF(IA)*ZLOCS(ILOC)
    8             CONTINUE !loop on anchor nodes
!
!................ just one anchor node attached to boundary
                  IF (KEEP.EQ.1) THEN
                     IF (DABS(X - XLOCS(IA1)).LT.TOL .AND.
     ;                   DABS(Z - ZLOCS(IA1)).LT.TOL) GOTO 400
                  ELSEIF (KEEP.EQ.2) THEN !along line?
                     V2(1) = X - XLOCS(IA1) !line cant be parallel
                     V2(2) = 0.D0
                     V2(3) = Z - ZLOCS(IA1)
                     DET = V1(1)*V2(3) - V1(3)*V2(1)
                     IF (DABS(DET).LT.TOL) GOTO 400
                  ENDIF
!
!................ go looking for (z,x) 
                  JNPE = IBSECT8(NNPE,1,TOL,Z,ZLOCSE)  
                  DO 9 INPE=JNPE,NNPE
                     IF (DABS(ZLOCSE(INPE) - Z).LT.TOL .AND. 
     ;                   DABS(XLOCSE(INPE) - X).LT.TOL) THEN 
                        GOTO 90 
                     ENDIF
    9             CONTINUE !Loop on nodal points
                  WRITE(*,*) 'cdof2xze: Warning brute forcing...'
                  print *, jnpe,nnpe
                  DO 10 INPE=1,NNPE
                     IF (DABS(ZLOCSE(INPE) - Z).LT.TOL .AND.
     ;                   DABS(XLOCSE(INPE) - X).LT.TOL) THEN
                        GOTO 90
                     ENDIF 
   10             CONTINUE !Loop on brute force
                  WRITE(*,*) 'cdof2xze: Error cant locate point'
                  IERR = 1
                  GOTO 70 
   90             CONTINUE !found location
                  DO 11 I=1,NDIM
                     IDOFSE(I,INPE) = LM(I,IAE,IELEM) 
   11             CONTINUE
  400             CONTINUE !break ahead; not in right domain
    7          CONTINUE !loop on xi
    6       CONTINUE !loop on eta
         ENDIF !end check on domain 
    2 CONTINUE !Loop on elements
   70 CONTINUE !error break
      DEALLOCATE(SF) 
      IF (IERR.NE.0) RETURN
      return
!
!.... fidelity check
      DO 21 INPE=1,NNPE
         DO 22 I=1,NDIM
            IF (IDOFSE(I,INPE).EQ.LNULL) THEN
               WRITE(*,*) 'codf2xze: Could not initialize point:',INPE
               IERR = 1
            ENDIF
   22    CONTINUE
   21 CONTINUE 
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE G1DBLK(MNPG,MGNOD, NNPG,NWORK,NELEM,NLXI, NGNOD, CSIDE,
     ;                  CDOMAIN,CNNPG, IENG, XIPTS,XLOCS,ZLOCS,DENS,
     ;                  QP,QS,ECOEFF, 
     ;                  NL1D, XBLK,POFF, Z1D,VP1D,VS1D,RH1D,QP1D,QS1D)
!
!     Generates the 1D model corresponding to the side of the Bielak
!     interior domain 
! 
!     INPUT      MEANING
!     -----      ------- 
!     CDOMAIN    holds the element's domain
!     CNNPG      holds the anchor nodes domain
!     CSIDE      side of Bielak boundary we are interested in 
!     DENS       density at anchor nodes
!     ECOEFF     elastic coefficients at anchor nodes
!     IENG       global IEN matrix 
!     MGNOD      leading dimension for anchor node
!     MNPG       leading dimension for ECOEFF 
!     NELEM      number of elements in model
!     NGNOD      number of anchor nodes on element
!     NLXI       number of lagrange interpolant points in xi/eta
!     NNPG       number of anchor nodes in mesh
!     NWORK      max size for 1D models
!     QP         P quality factor
!     QS         S quality factor
!     XIPTS      xi interpolation points
!     XLOCS      x locations of anchor nodes
!     ZLOCS      z locations of anchor nodes
!
!     OUTPUT     MEANING 
!     ------     ------- 
!     NL1D       number of layers in 1D model
!     POFF       offset distance to use in haskel functions
!     RH1D       1D density model
!     VP1D       1D vp velocity model
!     VS1D       1D vs velocity model
!     XBLK       x location of the x bielak/internal boundary
!     QP1D       1D P quality factor
!     QS1D       1D S quality factor
!     Z1D        1D interface depths
!
!.... variable declarations
      CHARACTER(2), INTENT(IN) :: CNNPG(NNPG)
      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM), CSIDE
      REAL*8, INTENT(IN) :: ECOEFF(MNPG,*), DENS(NNPG), XLOCS(NNPG), 
     ;                      QP(NNPG), QS(NNPG), ZLOCS(NNPG), XIPTS(NLXI)
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MNPG,MGNOD, NNPG, 
     ;                         NWORK,NELEM,NLXI, NGNOD 
      REAL*8, INTENT(OUT) :: Z1D(NWORK), VP1D(NWORK), VS1D(NWORK), 
     ;                       RH1D(NWORK), QP1D(NWORK), QS1D(NWORK), 
     ;                       XBLK, POFF  
      INTEGER*4, INTENT(OUT) :: NL1D
!.... local variables
      REAL*8, ALLOCATABLE :: ZWORK(:), VPWRK(:), VSWRK(:), RHWRK(:), 
     ;                       QPWRK(:), QSWRK(:)
      REAL*8 DSF1D(2), SF1D(2), XTARG, Z,RLAM,RMU,RHO,XI, 
     ;       ALPHA,BETA,QPP,QSS,VPBASE,VSBASE,DNBASE,ZMIN,TOL
      INTEGER*4 IANCH2(2) 
      PARAMETER(TOL = 1.D-10) 
!
!----------------------------------------------------------------------!
! 
!.... generate points around boundary
      NIBLK = ICNIBLK(NNPG,CNNPG) !number of anchor nodes on Bielak bdry
      NELBLK = NIBLK !- 1 if top isn't connected so this is a worst case
      VPBASE = 0.D0
      VSBASE = 0.D0
      DNBASE = 0.D0
      IF (CSIDE.EQ.'L'.OR.CSIDE.EQ.'l') THEN
         XTARG = HUGE(1.D0) 
         POFF = HUGE(1.D0)
         DO 1 INPG=1,NNPG
            IF (CNNPG(INPG).EQ.'BI') XTARG = DMIN1(XTARG,XLOCS(INPG))
            IF (CNNPG(INPG).EQ.'BA') POFF  = DMIN1(POFF ,XLOCS(INPG))
    1    CONTINUE
         XBLK = XTARG 
      ELSE
         XTARG =-HUGE(1.D0)
         POFF =-HUGE(1.D0)
         DO 3 INPG=1,NNPG
            IF (CNNPG(INPG).EQ.'BI') XTARG = DMAX1(XTARG,XLOCS(INPG))
            IF (CNNPG(INPG).EQ.'BA') POFF  = DMAX1(POFF ,XLOCS(INPG))
    3    CONTINUE
         XBLK = XTARG
      ENDIF
      ZMIN = HUGE(1.D0)
      DO 4 INPG=1,NNPG
         IF (CNNPG(INPG).EQ.'BA') ZMIN = DMIN1(ZMIN,ZLOCS(INPG)) 
    4 CONTINUE 
      VPBASE = 0.D0
      VSBASE = 0.D0
      DNBASE = 0.D0
      NSPACE = (NELBLK - 1)*(NLXI - 1) + NLXI + 1 !+1 for base add on
      ALLOCATE(ZWORK(NSPACE))
      ALLOCATE(VPWRK(NSPACE))
      ALLOCATE(VSWRK(NSPACE))
      ALLOCATE(RHWRK(NSPACE)) 
      ALLOCATE(QPWRK(NSPACE))
      ALLOCATE(QSWRK(NSPACE)) 
!
!.... loop on bielak elements
      INPS = 0 
      DO 11 IELEM=1,NELEM 
         IF (CDOMAIN(IELEM).EQ.'E') THEN !bielak domain
            IANCH = 0
            IANCH2(1:2) = 0 
            DO 12 IA=1,NGNOD
               ILOC = IENG(IA,IELEM)
               IF (DABS(XLOCS(ILOC) - XTARG).LT.TOL) THEN 
                  IANCH = IANCH + 1
                  IANCH2(IANCH) = ILOC
               ENDIF
   12       CONTINUE 
!
!.......... along two anchor points interpolate velocity and position
            IF (IANCH.EQ.2) THEN
               DO 13 ILXI=1,NLXI
                  XI = XIPTS(ILXI)
                  CALL CSF1D(XI, SF1D,DSF1D)
                  Z = 0.D0
                  RLAM = 0.D0
                  RMU  = 0.D0
                  RHO  = 0.D0
                  QPP  = 0.D0
                  QSS  = 0.D0 
                  DO 14 IA1D=1,2
                     ILOC1D = IANCH2(IA1D) 
                     Z = Z + SF1D(IA1D)*ZLOCS(ILOC1D)
                     RLAM = RLAM + SF1D(IA1D)*ECOEFF(ILOC1D,1)
                     RMU  = RMU + SF1D(IA1D)*ECOEFF(ILOC1D,2)
                     RHO  = RHO + SF1D(IA1D)*DENS(ILOC1D)
                     QPP  = QPP + SF1D(IA1D)*QP(ILOC1D)
                     QSS  = QSS + SF1D(IA1D)*QS(ILOC1D)
   14             CONTINUE
                  ALPHA = DSQRT( (RLAM + 2.D0*RMU)/RHO )
                  BETA = DSQRT(RMU/RHO)
                  IF (INPS.EQ.0) THEN
                     INPS = INPS + 1 
                     ZWORK(INPS) = Z 
                     VPWRK(INPS) = ALPHA
                     VSWRK(INPS) = BETA
                     RHWRK(INPS) = RHO 
                     QPWRK(INPS) = QPP
                     QSWRK(INPS) = QSS
                  ELSE 
                     DO 15 JNPS=1,INPS
                        IF (DABS(Z - ZWORK(JNPS)).LT.TOL) GOTO 50
   15                CONTINUE
                     INPS = INPS + 1 
                     ZWORK(INPS) = Z 
                     VPWRK(INPS) = ALPHA
                     VSWRK(INPS) = BETA
                     RHWRK(INPS) = RHO 
                     QPWRK(INPS) = QPP
                     QSWRK(INPS) = QSS
   50                CONTINUE !break ahead, point exist 
                  ENDIF 
   13          CONTINUE !loop on xi points
            ENDIF !end check on side
         ENDIF !end check on bielak domain 
   11 CONTINUE 
!
!.... sort 
      NNPS = INPS
      CALL SHELL6(NNPS, ZWORK,VPWRK,VSWRK,RHWRK,QPWRK,QSWRK)
!
!.... initalize base 
      IL1D = 1
      Z1D (1) = ZMIN !MINVAL(ZLOCS)
      VP1D(1) = VPWRK(1) !tecnically should be a halfspace
      VS1D(1) = VSWRK(1) !below so we have no added waves 
      RH1D(1) = RHWRK(1) !coming in later
      QP1D(1) = QPWRK(1) 
      QS1D(1) = QSWRK(1) 
!.... check that you aren't dressed like the guy in front of you
      DO 21 INPS=1,NNPS-1
         IF (DABS( VPWRK(INPS+1) - VPWRK(INPS)).GT.TOL .OR.
     ;       DABS( VSWRK(INPS+1) - VSWRK(INPS)).GT.TOL .OR.
     ;       DABS( RHWRK(INPS+1) - RHWRK(INPS)).GT.TOL .OR. 
     ;       DABS( QPWRK(INPS+1) - QPWRK(INPS)).GT.TOL .OR.
     ;       DABS( QSWRK(INPS+1) - QSWRK(INPS)).GT.TOL) THEN
            IL1D = IL1D + 1
            Z1D (IL1D) = ZWORK(INPS)
            VP1D(IL1D) = VPWRK(INPS)
            VS1D(IL1D) = VSWRK(INPS)
            RH1D(IL1D) = RHWRK(INPS)
            QP1D(IL1D) = QPWRK(INPS) 
            QS1D(IL1D) = QSWRK(INPS) 
         ENDIF
   21 CONTINUE
!.... cap off a good night
      NL1D = IL1D + 1
      Z1D (NL1D) = ZWORK(NNPS)
      VP1D(NL1D) = VPWRK(NNPS) 
      VS1D(NL1D) = VSWRK(NNPS) 
      RH1D(NL1D) = RHWRK(NNPS) 
      QP1D(NL1D) = QPWRK(NNPS)
      QS1D(NL1D) = QSWRK(NNPS)
!.... clean
      DEALLOCATE(ZWORK)
      DEALLOCATE(VPWRK)
      DEALLOCATE(VSWRK)
      DEALLOCATE(RHWRK)
      DEALLOCATE(QPWRK)
      DEALLOCATE(QSWRK)
!
!.... sort it out starting at base and working up
      CALL SHELL6(NL1D, Z1D,VP1D,VS1D,RH1D, QP1D,QS1D)
!
!.... swap for printing and easy reading by seismologists
      J = NL1D + 1
      DO 22 I=1,NL1D/2
         J = J - 1
         CALL SWAPR8(Z1D (I),Z1D (J))
         CALL SWAPR8(VP1D(I),VP1D(J))
         CALL SWAPR8(VS1D(I),VS1D(J))
         CALL SWAPR8(RH1D(I),RH1D(J))
         CALL SWAPR8(QP1D(I),QP1D(J))
         CALL SWAPR8(QS1D(I),QS1D(J))
   22 CONTINUE
!
!.... write out a summary
      WRITE(*,*)  
      IF (CSIDE.EQ.'L' .OR. CSIDE.EQ.'l') THEN
         WRITE(*,*) 'g1dblk: Left Bielak 1D model summary:'
      ELSE
         WRITE(*,*) 'g1dblk: Right Bielak 1D model summary:'
      ENDIF
      WRITE(*,*) 
      WRITE(*,900)
      WRITE(*,901)
      DO 23 I=1,NL1D
         WRITE(*,902) I,Z1D(I),VP1D(I),VS1D(I),RH1D(I),QP1D(I),QS1D(I)
   23 CONTINUE
      WRITE(*,*) 
  900 FORMAT('     Interface      Z         Vp        Vs      Density   
     ; Qp       Qs')
  901 FORMAT('     ---------  ---------   -------   -------   -------  
     ;------  ------')
  902 FORMAT(6X,I5,1F14.2,3F10.2,2F8.2)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CNMESH(MDIM,MEN, NDOF,NELEM, NDIM,NEN, CDOMAIN,LM, 
     ;                  CNP,IERR)
!
!     Assigns a value to each degree of freedom if it is in the 'I'
!     interior domain, 'E' exterior, 'A' artificial boundary, or 'B'
!     for boundary of 'I' and 'E'
!
!     INPUT      MEANING
!     -----      ------- 
!     CDOMAIN    elmenet domain
!     LM         location matrix
!     MDIM       leading dimension for LM
!     MEN        leading dimension for LM
!     NDIM       number of spatial dimensions of solution
!     NDOF       nubmer of degrees of freedom
!     
!     OUTPUT     MEANING
!     ------     ------- 
!     CNP        each DOFs position in mesh, see above
!     IERR       error flag
!
!.... variable declarations
      CHARACTER*1, INTENT(IN) :: CDOMAIN(NELEM)
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), MDIM,MEN, NDIM,NELEM 
      CHARACTER*1, INTENT(OUT) :: CNP(NDOF)
      INTEGER*4, INTENT(OUT) :: IERR
! 
!----------------------------------------------------------------------!
!
!.... first tag the absorbing conditions
      IERR = 0 
      IDOF = 0 !compiler not happy how i programmed this
      CNP(1:NDOF) = ' ' 
      DO 1 IELEM=1,NELEM
         IF (CDOMAIN(IELEM).EQ.'A') THEN
            DO 2 IAE=1,NEN 
               DO 3 I=1,NDIM
                  IDOF = LM(I,IAE,IELEM)
                  IF (IDOF.GT.0) CNP(IDOF) = 'A' 
    3          CONTINUE
    2       CONTINUE
         ENDIF !end check domain is in absorbing layer
    1 CONTINUE
!
!.... now tag the extended elements
      DO 4 IELEM=1,NELEM
         IF (CDOMAIN(IELEM).EQ.'E') THEN
            DO 5 IAE=1,NEN
               DO 6 I=1,NDIM
                  IDOF = LM(I,IAE,IELEM)
                  IF (IDOF.GT.0) THEN
                     IF (CNP(IDOF).EQ.' ') CNP(IDOF) = 'E' 
                  ENDIF
    6          CONTINUE
    5       CONTINUE
          ENDIF !end check domain is in exterior
    4 CONTINUE
!
!.... now tag the elements in the interior
      DO 7 IELEM=1,NELEM
         IF (CDOMAIN(IELEM).EQ.'I') THEN
            DO 8 IAE=1,NEN
               DO 9 I=1,NDIM
                  IDOF = LM(I,IAE,IELEM)
                  IF (IDOF.GT.0) THEN
                     IF (CNP(IDOF).EQ.' ') THEN
                        CNP(IDOF) = 'I'
                     ELSE
                        IF (CNP(IDOF).EQ.'E') THEN
                           CNP(IDOF) = 'B'
                        ELSE !element exists, override? 
                           IF (CNP(IDOF).EQ.'A') THEN
                              WRITE(*,*) 'cnmesh: Error, what???'
                              IERR = 1
                           ENDIF !should not be an absorbing element
                        ENDIF !end check if dof exists in exterior
                     ENDIF !end check if dof exists 
                  ENDIF
    9          CONTINUE !loop on spatial dimensions
    8       CONTINUE !loop on element nodes
         ENDIF !end that element is domain interior
    7 CONTINUE !loop on elements
!
!.... fidelity check, make sure all DOFs initialized
      DO 10 IDOF=1,NDOF
         IF (CNP(IDOF).EQ.' ') THEN
            WRITE(*,*) 'cnmesh: Error DOF',IDOF,'not initialized'
            IERR = 2
         ENDIF
   10 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      INTEGER*4 FUNCTION ICNIBLK(NNPG,CNNPG)
!     Counts number of anchor nodes at the Bielak/Interior boundary
      CHARACTER*2, INTENT(IN) :: CNNPG(NNPG)
      INTEGER*4, INTENT(IN) :: NNPG
      ICNIBLK = 0
      DO 1 INPG=1,NNPG
         IF (CNNPG(INPG).EQ.'BI') ICNIBLK = ICNIBLK + 1
    1 CONTINUE 
      RETURN
      END  
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GET_XZ_IABLK(NNPG,NIBLK, CNNPG, XLOCS,ZLOCS, XBLK,ZBLK)
!
!     Gets the anchor noes on the Bielak interior boundary
!
!     INPUT      MEANING
!     -----      ------- 
!     CNNPG      holds anchor nodes domain
!     NIBLK      number of anchor nodes on bielak interior bdry 
!     NNPG       number of anchor nodes in mesh
!     XLOCS      x locations of anchor nodes
!     ZLOCS      z locations of anchor nodes
!
!     OUTPUT     MEANING
!     ------     ------- 
!     XBLK       x anchor nodes on Bielak/interior boundary
!     ZBLK       z anchor nodes on Bielak/interior boundary 
!
!.... variable declaratinos
      CHARACTER*2, INTENT(IN) :: CNNPG(NNPG)
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG)
      INTEGER*4, INTENT(IN) :: NNPG, NIBLK
      REAL*8, INTENT(OUT) :: XBLK(NIBLK), ZBLK(NIBLK)
!
!----------------------------------------------------------------------!
!
      IABLK = 0 
      DO 1 INPG=1,NNPG
         IF (CNNPG(INPG).EQ.'BI') THEN
            IABLK = IABLK + 1
            XBLK(IABLK) = XLOCS(INPG)
            ZBLK(IABLK) = ZLOCS(INPG)
         ENDIF
    1 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CSF1D(XI, SF,DSF)
!
!     Evaluates linear shape fn and derivative at xi 
      REAL*8, INTENT(IN) :: XI
      REAL*8, INTENT(OUT) :: SF(2), DSF(2)
      SF(1) = 0.5D0*(1.D0 - XI) 
      SF(2) = 0.5D0*(1.D0 + XI) 
      DSF(1) =-0.5D0
      DSF(2) = 0.5D0
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SWAPR8(X,Y) !swaps x and y
      REAL*8, INTENT(INOUT) :: X, Y
      REAL*8 TEMP
      TEMP = X 
      X = Y 
      Y = TEMP
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SWAPI4(X,Y) !swaps x and y
      INTEGER*4, INTENT(INOUT) :: X, Y
      INTEGER*4 TEMP
      TEMP = X 
      X = Y 
      Y = TEMP
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE ADDBLK(NDOF,CNP,U0, SOL) 
!     Adds the Bielak conditions back into the total solution
      CHARACTER*1, INTENT(IN) :: CNP(NDOF)
      COMPLEX*8, INTENT(IN) :: U0(NDOF)
      INTEGER*4, INTENT(IN) :: NDOF
      COMPLEX*8, INTENT(INOUT) :: SOL(NDOF)
      DO 1 IDOF=1,NDOF
         IF (CNP(IDOF).EQ.'E') SOL(IDOF) = SOL(IDOF) + U0(IDOF)
    1 CONTINUE
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*8 FUNCTION DZBASE(NNPG,JOB,CNNPG,ZLOCS)
!
!     Calculates the z model base at the bielak interior boundary 
!     (job = 1) or z model base at the bielak absorbing boundary 
!     (default).  
!
!     INPUT       MEANING
!     -----       -------
!     CNNPG       character descirptor of anchor nodes
!     JOB         = 1 base of interior domain
!                 otherwise base of bielak domain
!     NNPG        number of anchor nodes
!     ZLOCS       z locations of anchor nodes
!
!     OUTPUT      MEANING
!     ------      ------- 
!     DZBASE      base of model at desired point
!
!.... variable declarations
      IMPLICIT NONE 
      CHARACTER(2) CNNPG(NNPG) 
      REAL*8, INTENT(IN) :: ZLOCS(NNPG) 
      INTEGER*4, INTENT(IN) :: NNPG, JOB 
!.... local variables
      REAL*8 ZBASE
      INTEGER*4 INPG 
!
!----------------------------------------------------------------------!
!
      ZBASE = HUGE(1.D0)
      DO 1 INPG=1,NNPG
         IF (JOB.EQ.1) THEN !bielak interior boudnary
            IF (CNNPG(INPG) == 'BI') THEN !bielak interior z bdry 
               ZBASE = DMIN1(ZBASE,ZLOCS(INPG))
            ENDIF
         ELSE
            IF (CNNPG(INPG) == 'BA') THEN !bielak absorbing z bdry
               ZBASE = DMIN1(ZBASE,ZLOCS(INPG))
            ENDIF
         ENDIF
    1 CONTINUE
      DZBASE = ZBASE
      RETURN
      END
