      SUBROUTINE SRCLIST(NSRC,NL1D_LT,NL1D_RT,  AZTOL,AOITOL, 
     ;                SRCTYP, BAZN,AOI, VP1D_LT,VS1D_LT,VP1D_RT,VS1D_RT,
     ;                NSG, CSIDE,ISGPTR,ISRCPRM, PYAVG,IERR) 
!
!     Here we generate a source list to minimize the number of 
!     fractorizations needed.  Effectively we say that  similar 
!     py values will result in similar factorizations.  This 
!     spares the user, aka me, from thinking.  Also, for simplicity 
!     just take which side of the model things are coming from
!
!     INPUT      MEANING
!
!     OUTPUT     MEANING
!
!.... variable declarations  
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC) 
      REAL*8, INTENT(IN) :: VP1D_RT(NL1D_RT),VS1D_RT(NL1D_RT), 
     ;                      VP1D_LT(NL1D_LT),VS1D_LT(NL1D_LT), 
     ;                      BAZN(NSRC), AOI(NSRC), AZTOL, AOITOL 
      INTEGER*4, INTENT(IN) :: NSRC,NL1D_LT,NL1D_RT 
      CHARACTER(1), INTENT(OUT) :: CSIDE(NSRC) 
      REAL*8, INTENT(OUT) :: PYAVG(NSRC)
      INTEGER*4, INTENT(OUT) :: ISGPTR(NSRC+1), ISRCPRM(NSRC), NSG
!.... local variables
      REAL*8, ALLOCATABLE :: PY(:) 
      REAL*8 PYT, PY0, PTOL, VBASE, VFAST, PI180
      PARAMETER(PI180 = 0.017453292519943295D0)
      REAL*8 CVBASE
      LOGICAL*4 LABS
      PARAMETER(LABS = .TRUE.)
!
!----------------------------------------------------------------------!
!
!.... calculate the py for every source 
      IERR = 0
      VFAST = 0.D0
      ALLOCATE(PY(NSRC))
      J1 = 1
      JSRC = 0 
      DO 1 ISGN=1,2
         DO 2 ISRC=1,NSRC
            IF (DCOS(BAZN(ISRC)*PI180).LT.0.D0) THEN
               CSIDE(ISRC) = 'R'
            ELSE
               CSIDE(ISRC) = 'L'
            ENDIF
            VBASE = CVBASE(SRCTYP(ISRC),CSIDE(ISRC), NL1D_LT,NL1D_RT, 
     ;                     VP1D_LT,VS1D_LT, VP1D_RT,VS1D_RT) 
!           IF (CSIDE(ISRC).EQ.'L') THEN
!              IF (SRCTYP(ISRC)(2:2).EQ.'P' .OR. 
!    ;             SRCTYP(ISRC)(2:2).EQ.'p') THEN
!                 VBASE = VP1D_LT(NL1D_LT)
!              ELSE
!                 VBASE = VS1D_LT(NL1D_LT)
!              ENDIF
!           ELSE
!              IF (SRCTYP(ISRC)(2:2).EQ.'P' .OR.
!    ;             SRCTYP(ISRC)(2:2).EQ.'p') THEN
!                 VBASE = VP1D_RT(NL1D_RT)
!              ELSE
!                 VBASE = VS1D_RT(NL1D_RT)
!              ENDIF
!           ENDIF
            VFAST = DMAX1(VFAST,VBASE)
            PYT = DSIN(BAZN(ISRC)*PI180)*DSIN(AOI(ISRC)*PI180)/VBASE
            if (labs) pyt = dabs(pyt) !only need to come from one direction
            IF (DABS(PYT).LT.1.D-10) PYT = 0.D0
            IF (ISGN.EQ.1.AND.PYT.GE.0.D0) THEN
               JSRC = JSRC + 1
               PY(JSRC) = PYT
               ISRCPRM(JSRC) = ISRC
            ENDIF
            IF (ISGN.EQ.2.AND.PYT.LT.0.D0) THEN
               JSRC = JSRC + 1
               PY(JSRC) = PYT
               ISRCPRM(JSRC) = ISRC
            ENDIF
    2    CONTINUE 
         IF (ISGN.EQ.1) THEN
            NV = JSRC - J1 + 1
            IF (NV.GT.0) THEN
               CALL SHELLR8I2(NV,PY(J1:JSRC),ISRCPRM(J1:NSRC))
               J1 = JSRC + 1
            ELSE
               J1 = 1
            ENDIF
         ELSE
            NV = NSRC - J1 + 1
            IF (NV.GT.0) THEN
               CALL SHELLR8I2(NV,PY(J1:NSRC),ISRCPRM(J1:NSRC))
               K1 = NSRC + 1
               DO 3 I1=J1,(J1+NSRC)/2
                  K1 = K1 - 1
                  CALL SWAPR8(PY(I1),PY(K1))
                  CALL SWAPI4(ISRCPRM(I1),ISRCPRM(K1))
    3          CONTINUE
            ENDIF
          ENDIF 
    1 CONTINUE !loop on sources
!
!.... now calculate a py tolerance, vfast makes this small 
      PTOL = DABS(DSIN(AOITOL*PI180)*DSIN(AZTOL*PI180)/VFAST)  
!
!.... generate a pointer list to loop through sources
      ISGPTR(1:NSRC+1) = 0
      NSG = 1
      ISGPTR(1) = 1
      IF (NSRC.EQ.1) THEN
         ISGPTR(2) = 2
      ELSE
         PY0 = DMAX1(PY(1),1.D-10) !won't switch signs at 0
         DO 11 JSRC=1,NSRC-1
            IF (DABS(PY(JSRC+1) - PY0).GT.PTOL .OR.
     ;          PY(JSRC+1)*PY0.LT.0.D0) THEN
               ISGPTR(NSG+1) = JSRC + 1
               NSG = NSG + 1 
               PY0 = PY(JSRC+1)
            ENDIF
  11     CONTINUE
         IF (ISGPTR(NSG).EQ.NSRC) THEN
            NSG = MIN0(NSRC,NSG) 
            ISGPTR(NSG+1) = NSRC + 1
         ELSE
            NSG = MIN0(NSRC,NSG) 
            ISGPTR(NSG+1) = NSRC + 1 
         ENDIF
      ENDIF 
!
!.... create the average py list
      DO 21 ISG=1,NSG
         PYAVG(ISG) = 0.D0
         NHIT = 0
         DO 22 JSRC=ISGPTR(ISG),ISGPTR(ISG+1)-1 
            NHIT = NHIT + 1
            PYAVG(ISG) = PYAVG(ISG) + PY(JSRC) 
   22    CONTINUE
         IF (NHIT.EQ.0) THEN
            WRITE(*,*) 'srclist: You have an empty source group'
            IERR = 1
            RETURN
         ENDIF
         PYAVG(ISG) = PYAVG(ISG)/DFLOAT(NHIT)
         IF (LABS) PYAVG(ISG) =-PYAVG(ISG)
   21 CONTINUE 
      DEALLOCATE(PY)
!
!.... quality control, generate a source list permulation list
      WRITE(*,*) 'srclist: Permuation list:'
      DO 31 ISG=1,NSG
         DO 32 JSRC=ISGPTR(ISG),ISGPTR(ISG+1)-1
            ISRC = ISRCPRM(JSRC)
            VBASE = CVBASE(SRCTYP(ISRC),CSIDE(ISRC), NL1D_LT,NL1D_RT,
     ;                     VP1D_LT,VS1D_LT, VP1D_RT,VS1D_RT) 
            PYT = DSIN(BAZN(ISRC)*PI180)*DSIN(AOI(ISRC)*PI180)/VBASE
            IF (DABS(PYT).LT.1.11D-16) PYT = 0.D0
            IF (LABS) PYT = DABS(PYT) 
            WRITE(*,900) ISRC,JSRC,PYT,ISG
   32    CONTINUE
   31 CONTINUE
  900 FORMAT('          Source',I3,' is now source',I3,
     ;       ' with py',E12.4,' and in group',I3)
      WRITE(*,*) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SRCLIST2(NSRC,NL1D_LT,NL1D_RT,  AZTOL,AOITOL, 
     ;                SRCTYP, BAZN,AOI, VP1D_LT,VS1D_LT,VP1D_RT,VS1D_RT,
     ;                NSG, CSIDE,ISGPTR,ISRCPRM, PYAVG,IERR) 
!
!     Here we generate a source list to minimize the number of 
!     factorizations.  Effectively we do not care which direction the 
!     wave is coming from, but only into which py bin it falls.  Now 
!     we must be mindful of course at the end of processing to consider 
!     the wave directionality for the appropriate station correction
!     - B Baker December 2012
!
!.... variable declarations  
      IMPLICIT REAL*8 (A-H,O-Z) 
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC)
      REAL*8, INTENT(IN) :: VP1D_RT(NL1D_RT),VS1D_RT(NL1D_RT),
     ;                      VP1D_LT(NL1D_LT),VS1D_LT(NL1D_LT),
     ;                      BAZN(NSRC), AOI(NSRC), AZTOL, AOITOL
      INTEGER*4, INTENT(IN) :: NSRC,NL1D_LT,NL1D_RT
      CHARACTER(1), INTENT(OUT) :: CSIDE(NSRC)
      REAL*8, INTENT(OUT) :: PYAVG(NSRC)
      INTEGER*4, INTENT(OUT) :: ISGPTR(NSRC+1),ISRCPRM(NSRC),NSG 
!.... local variables
      REAL*8, ALLOCATABLE :: PY(:)
      INTEGER*4, ALLOCATABLE :: MYGROUP(:) 
      REAL*8 TOL
      PARAMETER(PI180 = 0.017453292519943295D0)
!
!----------------------------------------------------------------------!
!
!.... this would be pointless 
      IF (NSRC.LE.0) THEN
         WRITE(*,*) 'srclist: You do not have any sources!'
         IERR = 1
         RETURN
      ENDIF
!
!.... loop on the sources, calculates pys, and vbase for ptol 
      IERR = 0
      TOL = EPSILON(1.D0)*20.D0
      VFAST = 0.D0
      ALLOCATE(PY(NSRC))
      DO 1 ISRC=1,NSRC
         IF (DCOS(BAZN(ISRC)*PI180).LT.0.D0) THEN
            CSIDE(ISRC) = 'R' 
         ELSE
            CSIDE(ISRC) = 'L' 
         ENDIF
         VBASE = CVBASE(SRCTYP(ISRC),CSIDE(ISRC), NL1D_LT,NL1D_RT, 
     ;                  VP1D_LT,VS1D_LT, VP1D_RT,VS1D_RT) 
         PY(ISRC) = DSIN(BAZN(ISRC)*PI180)*DSIN(AOI(ISRC)*PI180)/VBASE
         PY(ISRC) = DABS(PY(ISRC)) !dont care for source groups
         IF (DABS(PY(ISRC)).LT.TOL) PY(ISRC) = 0.D0
         ISRCPRM(ISRC) = ISRC 
         VFAST = DMAX1(VFAST,VBASE)
    1 CONTINUE
!
!.... define py tolerance and sort
      PTOL = DABS(DSIN(AOITOL*PI180)*DSIN(AZTOL*PI180)/VFAST)
      CALL SHELLR8I2(NSRC,PY,ISRCPRM)
!
!.... i'm going to do this like an idiot, but it needs to be right
      PY0 = PY(1) 
      ALLOCATE(MYGROUP(NSRC+1))
      MYGROUP(1:NSRC+1) = 0
      ISG = 1
      MYGROUP(1) = 1
      DO 2 ISRC=2,NSRC
         IF (PY(ISRC).GT.PY0+PTOL) THEN
            PY0 = PY(ISRC)
            ISG = ISG + 1
         ENDIF
         MYGROUP(ISRC) = ISG
    2 CONTINUE
      IF (MINVAL(MYGROUP(1:NSRC)).EQ.0) THEN
         WRITE(*,*) 'srclist: Error generating source list!'
         IERR = 1
         GOTO 750
      ENDIF 
!
!.... calculate the source groups
      ISGPTR(1:NSRC+1) = 0
      NSG = ISG 
      ISGPTR(1) = 1
      DO 3 ISRC=1,NSRC
         IF (MYGROUP(ISRC).NE.MYGROUP(ISRC+1)) THEN
            ISG = MYGROUP(ISRC)
            ISGPTR(ISG+1) = ISRC + 1
         ENDIF
    3 CONTINUE 
!
!.... calculate the averages
      WRITE(*,*)
      PYAVG(1:NSRC) = 0.D0 
      DO 4 ISG=1,NSG
         J1 = ISGPTR(ISG)
         J2 = ISGPTR(ISG+1) - 1
         DEN = DFLOAT(J2 - J1 + 1)
         PYAVG(ISG) = 0.D0
         DO 5 JSRC=J1,J2
            PYAVG(ISG) = PYAVG(ISG) + PY(JSRC) 
    5    CONTINUE
         PYAVG(ISG) = PYAVG(ISG)/DEN
         WRITE(*,905) ISG,PYAVG(ISG)
  905    FORMAT(' srclist: Source group:',I3,' average py:',E12.4)
    4 CONTINUE 
!
!.... tabulate my results
      DO 6 ISG=1,NSG
         J1 = ISGPTR(ISG)
         J2 = ISGPTR(ISG+1) - 1
         DO 7 JSRC=J1,J2
            ISRC = ISRCPRM(JSRC)
            WRITE(*,910) JSRC,ISRC
    7    CONTINUE 
  910    FORMAT(' srclist: Source:',I3,' is now source:',I3) 
    6 CONTINUE
  750 CONTINUE !break ahead for errors
      DEALLOCATE(PY) 
      DEALLOCATE(MYGROUP)
      WRITE(*,*) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SRCLIST3(NSRC,LVERB, VFAST,AZTOL,AOITOL, PYTAB, 
     ;                    NSG,ISGPTR,ISRCPRM,PYAVG)
!
!     Here we generate a source list to minimize the number of 
!     factorizations.  Effectively we do not care which direction the 
!     wave is coming from, but only into which py bin it falls.  Now 
!     we must be mindful of course at the end of processing to consider 
!     the wave directionality for the appropriate station correction. 
!     The difference between this and SRCLST2 is that the surface wave
!     calculation requires us to find the py values a priori so we 
!     simply tabulate and make use of them here.
!     - B Baker January 2013
!
!     INPUT      MEANING
!     -----      ------- 
!     AOITOL     angle of incidence tolerance (degrees)
!     AZTOL      azimuth tolerance (degrees)
!     LVERB      toggles verbosity level
!     NSRC       number of sources
!     PYTAB      apparent slowness (s/m) for these sources 
!     VFAST      fastest velocity in 1D models to be strict w/ ptol
!  
!     OUTPUT     MEANING
!     ------     -------
!     ISGPTR     source group ointer
!     ISRCPRM    source permuation list
!     NSG        number of source groups
!     PYAVG      average py in group
! 
!.... variable declarations
      REAL*8, INTENT(IN) :: PYTAB(NSRC), VFAST,  AZTOL,AOITOL 
      INTEGER*4, INTENT(IN) :: NSRC
      LOGICAL*4, INTENT(IN) :: LVERB
      REAL*8, INTENT(OUT) :: PYAVG(NSRC) 
      INTEGER*4, INTENT(OUT) :: ISGPTR(NSRC+1), ISRCPRM(NSRC), NSG 
!.... local variables
      REAL*8, ALLOCATABLE :: PY(:) 
      INTEGER*4, ALLOCATABLE :: MYGROUP(:)
      REAL*8 PY0, PTOL, DEN, PI180
      PARAMETER(PI180 = 0.017453292519943295D0)  
!
!----------------------------------------------------------------------!
!
!.... copy the absolute value of py to work
      ALLOCATE(PY(NSRC))
      DO 1 ISRC=1,NSRC
         PY(ISRC) = DABS(PYTAB(ISRC)) 
         ISRCPRM(ISRC) = ISRC 
    1 CONTINUE 
!
!.... define py tolerance and sort
      PTOL = DABS(DSIN(AOITOL*PI180)*DSIN(AZTOL*PI180)/VFAST)
      CALL SHELLR8I2(NSRC,PY,ISRCPRM)
!
!.... i'm going to do this like an idiot, but it needs to be right
      PY0 = PY(1) 
      ALLOCATE(MYGROUP(NSRC+1))
      MYGROUP(1:NSRC+1) = 0 
      ISG = 1 
      MYGROUP(1) = 1 
      DO 2 ISRC=2,NSRC
         IF (PY(ISRC).GT.PY0+PTOL) THEN
            PY0 = PY(ISRC)
            ISG = ISG + 1 
         ENDIF
         MYGROUP(ISRC) = ISG 
    2 CONTINUE
      IF (MINVAL(MYGROUP(1:NSRC)).EQ.0) THEN
         WRITE(*,*) 'srclist3: Error generating source list!'
         IERR = 1 
         GOTO 750 
      ENDIF 
!
!.... calculate the source groups
      ISGPTR(1:NSRC+1) = 0
      NSG = ISG
      ISGPTR(1) = 1
      DO 3 ISRC=1,NSRC
         IF (MYGROUP(ISRC).NE.MYGROUP(ISRC+1)) THEN
            ISG = MYGROUP(ISRC)
            ISGPTR(ISG+1) = ISRC + 1
         ENDIF
    3 CONTINUE
!
!.... calculate the averages
      IF (LVERB) WRITE(*,*)
      PYAVG(1:NSRC) = 0.D0
      DO 4 ISG=1,NSG
         J1 = ISGPTR(ISG)
         J2 = ISGPTR(ISG+1) - 1
         DEN = DFLOAT(J2 - J1 + 1)
         PYAVG(ISG) = 0.D0
         DO 5 JSRC=J1,J2
            PYAVG(ISG) = PYAVG(ISG) + PY(JSRC)
    5    CONTINUE
         PYAVG(ISG) = PYAVG(ISG)/DEN
         IF (LVERB) WRITE(*,905) ISG,PYAVG(ISG)
  905    FORMAT(' srclist: Source group:',I3,' average py:',E12.4)
    4 CONTINUE
!
!.... tabulate my results
      DO 6 ISG=1,NSG
         J1 = ISGPTR(ISG)
         J2 = ISGPTR(ISG+1) - 1
         DO 7 JSRC=J1,J2
            ISRC = ISRCPRM(JSRC)
            IF (LVERB) WRITE(*,910) JSRC,ISRC
    7    CONTINUE
  910    FORMAT(' srclist: Source:',I3,' is now source:',I3)
    6 CONTINUE
  750 CONTINUE !break ahead for errors
      DEALLOCATE(PY)
      DEALLOCATE(MYGROUP)
      IF (LVERB) WRITE(*,*)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SRCLIST4(NSRC,LVERB, VFAST,AZTOL,AOITOL, PYTAB,
     ;                    NSG,ISGPTR,ISRCPRM,PYAVG)
!
!     Here we generate a source list to minimize the number of 
!     factorizations.  The surface waves have frequency dependent py 
!     values which have already been tabulated.   
!     - B Baker January 2013
!
!     INPUT      MEANING
!     -----      ------- 
!     AOITOL     angle of incidence tolerance (degrees)
!     AZTOL      azimuth tolerance (degrees)
!     LVERB      toggles verbosity level
!     NSRC       number of sources
!     PYTAB      apparent slowness (s/m) for these sources 
!     VFAST      fastest velocity in 1D models to be strict w/ ptol
!  
!     OUTPUT     MEANING
!     ------     -------
!     ISGPTR     source group ointer
!     ISRCPRM    source permuation list
!     NSG        number of source groups
!     PYAVG      average py in group
! 
!.... variable declarations
      REAL*8, INTENT(IN) :: PYTAB(NSRC), VFAST,  AZTOL,AOITOL
      INTEGER*4, INTENT(IN) :: NSRC
      LOGICAL*4, INTENT(IN) :: LVERB
      REAL*8, INTENT(OUT) :: PYAVG(NSRC)
      INTEGER*4, INTENT(OUT) :: ISGPTR(NSRC+1), ISRCPRM(NSRC), NSG
!.... local variables
      REAL*8, ALLOCATABLE :: PY(:)
      INTEGER*4, ALLOCATABLE :: MYGROUP(:)
      REAL*8 PY0, PTOL, DEN, PI180
      LOGICAL*4 LNEG
      PARAMETER(PI180 = 0.017453292519943295D0)
!
!----------------------------------------------------------------------!
!
!.... copy the py values to workspace, and define tolerance
      PTOL = DABS(DSIN(AOITOL*PI180)*DSIN(AZTOL*PI180)/VFAST)
      ALLOCATE(PY(NSRC)) 
      DO 1 ISRC=1,NSRC
         ISRCPRM(ISRC) = ISRC
         PY(ISRC) = PYTAB(ISRC)  
    1 CONTINUE 
      CALL SHELLR8I2(NSRC,PY,ISRCPRM) !sort
!
!.... we now have two cases, all positive/negative or mixed
      ALLOCATE(MYGROUP(NSRC+1)) 
      MYGROUP(1:NSRC+1) = 0 !extra 0 indicates we are out
      ISG = 1 
      IF (MAXVAL(PY).LE.0.D0 .OR.  !all positive 
     ;    MINVAL(PY).GE.0.D0) THEN !all negative
         PY0 = PY(1) 
         DO 2 ISRC=1,NSRC 
            IF (PY(ISRC).GT.PY0+PTOL) THEN
               PY0 = PY(ISRC)
               ISG = ISG + 1
            ENDIF
            MYGROUP(ISRC) = ISG
    2    CONTINUE 
      ELSE !mixed  
         PY0 = PY(1) 
         LNEG = .TRUE. !start negative
         DO 3 ISRC=1,NSRC 
            IF (LNEG) THEN 
               IF (PY(ISRC).GT.PY0+PTOL .OR. PY(ISRC).GT.0.D0) THEN
                  PY0 = PY(ISRC) 
                  ISG = ISG + 1
                  IF (PY(ISRC).GT.0.D0) LNEG = .FALSE.
               ENDIF 
            ELSE
               IF (PY(ISRC).GT.PY0+PTOL) THEN
                  PY0 = PY(ISRC)
                  ISG = ISG + 1
               ENDIF
            ENDIF
            MYGROUP(ISRC) = ISG 
    3    CONTINUE  
      ENDIF
!
!.... quality check
      IF (MINVAL(MYGROUP(1:NSRC)).EQ.0) THEN
         WRITE(*,*) 'srclist4: Error generating source list!'
         IERR = 1
         GOTO 750
      ENDIF
!
!.... each source now belongs to a group, make CRS pointer
      ISGPTR(1:NSRC+1) = 0
      NSG = ISG
      ISGPTR(1) = 1
      DO 4 ISRC=1,NSRC
         IF (MYGROUP(ISRC).NE.MYGROUP(ISRC+1)) THEN
            ISG = MYGROUP(ISRC)
            ISGPTR(ISG+1) = ISRC + 1
         ENDIF
    4 CONTINUE
!
!.... calculate the average py values in each group
      IF (LVERB) WRITE(*,*)
      PYAVG(1:NSRC) = 0.D0
      DO 5 ISG=1,NSG
         J1 = ISGPTR(ISG)
         J2 = ISGPTR(ISG+1) - 1
         DEN = DFLOAT(J2 - J1 + 1)
         PYAVG(ISG) = 0.D0
         DO 6 JSRC=J1,J2
            PYAVG(ISG) = PYAVG(ISG) + PY(JSRC)
    6    CONTINUE
         PYAVG(ISG) = PYAVG(ISG)/DEN
         IF (LVERB) WRITE(*,905) ISG,PYAVG(ISG)
  905    FORMAT(' srclist: Source group:',I3,' average py:',E12.4)
    5 CONTINUE
!
!.... show the permuatation list 
      IF (LVERB) THEN
         DO 7 ISG=1,NSG
            J1 = ISGPTR(ISG)
            J2 = ISGPTR(ISG+1) - 1
            DO 8 JSRC=J1,J2
               ISRC = ISRCPRM(JSRC)
               IF (LVERB) WRITE(*,910) JSRC,ISRC
    8       CONTINUE
  910       FORMAT(' srclist: Source:',I3,' is now source:',I3)
    7    CONTINUE
      ENDIF
  750 CONTINUE !break ahead for errors
      IF (ALLOCATED(PY)) DEALLOCATE(PY)
      IF (ALLOCATED(MYGROUP)) DEALLOCATE(MYGROUP)
      IF (LVERB) WRITE(*,*)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*8 FUNCTION CVBASE(SRCTYP,CSIDE, NL1D_LT,NL1D_RT, 
     ;                       VP1D_LT,VS1D_LT, VP1D_RT,VS1D_RT)
!
!     Calculates velocity at base of model.  Useful for when  
!     calculating py
      CHARACTER(2), INTENT(IN) :: SRCTYP
      CHARACTER(1), INTENT(IN) :: CSIDE 
      REAL*8, INTENT(IN) :: VP1D_LT(NL1D_LT), VS1D_LT(NL1D_LT), 
     ;                      VP1D_RT(NL1D_RT), VS1D_RT(NL1D_RT)
      CVBASE = 0.D0
      IF (CSIDE.EQ.'L') THEN
         IF (SRCTYP(2:2).EQ.'P' .OR. 
     ;       SRCTYP(2:2).EQ.'p') THEN
            CVBASE = VP1D_LT(NL1D_LT)
         ELSE
            CVBASE = VS1D_LT(NL1D_LT)
         ENDIF
      ELSE
         IF (SRCTYP(2:2).EQ.'P' .OR.
     ;       SRCTYP(2:2).EQ.'p') THEN
            CVBASE = VP1D_RT(NL1D_RT)
         ELSE
            CVBASE = VS1D_RT(NL1D_RT)
         ENDIF
      ENDIF
      RETURN
      END 
