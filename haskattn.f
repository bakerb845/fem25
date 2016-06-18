c----------------------------------------------------------------------
c
c       Compute the Green's functions for a 1D medium using TH (Thomson-Haskell).  
c
c       This version allows for attenuating medium by using complex wavespeeds.
c
c       It also normalizes the stress using the shear modulus of the half space
c       in order to balance the magnitudes of the entries in the 'a' propagaion
c       matrix.
c
c       This routine is the computative "engine" of the T-H method.
c       It impliments the concepts in Haskell (1953; H53) and Haskell
c       (1962; H62) to compute the Green's functions in layered media
c       do to an external plane wave source of unit amplitude
c
c       Arguments:
c       df      The frequency 
c       isrc    source number used to determine source type (P or S)
c       IERR    error flag
c
c       This version for 1 frequency = df
c
c       This version allows the addition of a fluid layer (vs = 0)
c
c       Author:  S. Roecker 4/09
c
        subroutine haskattn(nl,nzpts,srctyp,lflip, df,freq0,ain, 
     ;                      z1d,vp1d,vs1d,rh1d,qa1d,qb1d, depth, 
     ;                      ugrn1f,wgrn1f, ierr)

        implicit none
        character(2), intent(in) :: srctyp
        real*8, intent(in) :: z1d(nl),vp1d(nl),vs1d(nl),rh1d(nl),
     ;                        qa1d(nl), qb1d(nl), depth(nzpts),ain,df,
     ;                        freq0
        integer*4, intent(in) :: nl,nzpts
        logical*4, intent(in) :: lflip
        complex*16, intent(out) :: ugrn1f(nzpts),wgrn1f(nzpts)
        integer*4, intent(out) :: ierr


        complex*16, allocatable :: r1(:), r2(:), aag(:), aag1(:), 
     ;                             xh(:), xh1(:), 
     ;                             s01(:),c01(:), s02(:), c02(:) 
        real*8, allocatable :: h1d(:), h(:), tt(:) 

        complex*16 alph, beta

        real*8 vsvpnl
        real*8 coa

        character(1) ctype

        complex*16  x,y,ag,ag1,tm,
     &          pd,qd,pm,qm,cp,cq
        complex*16  c11,c12,c21,c22,
     &          c31,c32,c41,c42
        complex*16  b11,b12,b21,b22,
     &          b31,b32,b41,b42
        complex*16  d11,d12,d21,d22,
     &          d31,d32,d41,d42
        complex*16  a11,a12,a13,a14,a21,a22,a23,a24,
     &          a31,a32,a33,a34,a41,a42,a43,a44

        complex*16 e11, e13, e22, e24, e33, e31, e42, e44,
     ;             z01, z02, eye, up, wp, cova, covb

        real*8 cv, pic, hc, dept, dloc, ztop, pi180, pi
        parameter(pi180 = 0.017453292519943295d0, 
     ;            pi = 3.1415926535897932384626434d0)
        parameter(eye = dcmplx(0.d0, 1.d0))

        integer*4 i, m, iz
        integer*4 ibov, ibovo, ifin, idir, izb, ize
        complex*16 vcmplx 

!
!...... check angle of incidence
        if (ain.ge.90.d0 .or. ain.lt.0.d0) then
           write(*,*) 'haskattn: Angle of incidence invalid:',ain
           ugrn1f(1:nzpts) = dcmplx(0.d0,0.d0)
           wgrn1f(1:nzpts) = dcmplx(0.d0,0.d0)
           ierr = 1
           return
        endif
c----incident wave type
        ierr = 0
        ctype = 'P'
        if (srctyp.eq.'PS' .or. srctyp.eq.'ps' .or. 
     ;      srctyp.eq.'pS' .or. srctyp.eq.'Ps') ctype = 'S'


c----horizontal wavespeed
c  NB:  it seems to me that if we do the dispersion correction we use the apparent
c  wavespeed from the original vp1d/vs1d and ain(isrc), so to be consistent we should
c  use vp1d and vs1d here.  Otherwise we should either redefine ain, or define cv elsewhere.
        if (ctype.eq.'P') then
          cv = vp1d(nl)/dsin(ain*pi180)
        else
          cv = vs1d(nl)/dsin(ain*pi180)
        endif

        pic = 2.d0*pi*df/cv
c       pic = df/cv

c       write(*,*) ' cv, pic = ', cv, pic
! 
!...... may need to flip coordinate system in z
        allocate(h(nl-1))
        h(1:nl-1) = 0.d0
        allocate(h1d(nl))
        h1d(1:nl) = 0.d0
        ztop = z1d(1)
        do i=1,nl
           if (lflip) then
              h1d(i) = ztop - z1d(i)
           else
              h1d(i) = z1d(i)
           endif
        enddo
c----thickness and complex wavespeeds array
        do i = 1, nl-1
          h(i) = h1d(i+1) - h1d(i)
        enddo
!
!...... set space
        allocate(r1(nl+1))
        allocate(r2(nl+1)) 
        allocate(aag(nl+1)) 
        allocate(aag1(nl+1)) 
        allocate(xh(nl+1)) 
        allocate(xh1(nl+1)) 
        allocate(s01(nl+1)) 
        allocate(c01(nl+1)) 
        allocate(s02(nl+1)) 
        allocate(c02(nl+1))
        allocate(tt(nl+1))
c----Loop over the model, setting up layer specific values
        do m = 1,nl

          alph = vcmplx(df,freq0,vp1d(m),qa1d(m)) !cmplx(vp1d(m), vp1di(m))
          beta = vcmplx(df,freq0,vs1d(m),qb1d(m)) !cmplx(vs1d(m), vs1di(m))

c--- c/alpha and c/beta
          cova = dcmplx(cv,0.d0)/alph
c---set up a bogus value for covb in case of vs = 0
          if (cdabs(beta).ne.0.d0) then
            covb = dcmplx(cv,0.d0)/beta
          else
            covb = dcmplx(1.d0,0.d0)
          endif 

c----r1 and r2 are r(alpha) and r(beta)
c----aag is gamma, aag1 is gamma-1
          x = cova*cova - dcmplx(1.d0,0.d0)
          r1(m) = cdsqrt(x)
          x = covb*covb - dcmplx(1.d0,0.d0)
          r2(m) = cdsqrt(x)
          aag(m)  = dcmplx(2.d0,0.d0)/(covb*covb)
          aag1(m) = aag(m) - dcmplx(1.d0,0.d0)

          if(m.lt.nl) then
            tt(m) = (rh1d(m)/rh1d(nl))*(cv/vs1d(nl))*(cv/vs1d(nl))
            xh(m) = aag(m)*aag(m)
            xh1(m) = aag1(m)*aag1(m)
            hc = pic*h(m)

c----z01, z02 are Pm and Qm in H53
            z01 = dcmplx(hc,0.d0)*r1(m)
            z02 = dcmplx(hc,0.d0)*r2(m)

            x = cdexp(eye*z01)
            c01(m) = dcmplx(0.5d0,0.d0)*(x + dcmplx(1.d0,0.d0)/x)
            s01(m) = dcmplx(0.5d0,0.d0)*(x - dcmplx(1.d0,0.d0)/x)

            x = cdexp(eye*z02)
            c02(m) = dcmplx(0.5d0,0.d0)*(x + dcmplx(1.d0,0.d0)/x)
            s02(m) = dcmplx(0.5d0,0.d0)*(x - dcmplx(1.d0,0.d0)/x)

c           write(*,*) m,  c01(m), s01(m), e2(m), c02(m), s02(m)
c            read(*,*) pause
          endif
        enddo

c---Now propagate to the surface to compute the surface displacements

c---The following should be equivalent to the inverse E matrix for
c       the half space (see Eqn 2.16 of H53). 

        vsvpnl = vs1d(nl)/vp1d(nl)
        vsvpnl = vsvpnl*vsvpnl
        coa = cv/vp1d(nl)

        e11 = -dcmplx(2.d0*vsvpnl,0.d0)
        e13 =  dcmplx(vsvpnl,0.d0)
        e22 =  dcmplx(coa*coa,0.d0)*aag1(nl)/r1(nl)
        e24 =  dcmplx(vsvpnl,0.d0)/r1(nl)
        e31 =  aag1(nl)/(aag(nl)*r2(nl))
        e33 = -dcmplx(1.d0,0.d0)/(dcmplx(2.d0,0.d0)*r2(nl))
        e42 =  dcmplx(1.d0,0.d0)
        e44 =  dcmplx(0.5d0,0.d0)

c       write(*,*) ' nl, x, e11, e33 = ', nl, x, e11, e33
c----The 4 x 2 b matrix will hold the multiplication
c    of all the a matrices.  This initialization allows
c    the first two columns of the a matrix to be saved
c    the first time around.

        b11 = dcmplx(1.d0,0.d0)
        b12 = dcmplx(0.d0,0.d0)
        b21 = dcmplx(0.d0,0.d0)
        b22 = dcmplx(1.d0,0.d0)
        b31 = dcmplx(0.d0,0.d0)
        b32 = dcmplx(0.d0,0.d0)
        b41 = dcmplx(0.d0,0.d0)
        b42 = dcmplx(0.d0,0.d0)

c----Loop over all layers
        do 7 m = 1,nl-1
          pd = s01(m)/r1(m)
          qd = s02(m)/r2(m)
          pm = s01(m)*r1(m)
          qm = s02(m)*r2(m)
c
c       c1, cp  = cos(Pm)
c       c2, cq  = cos(Qm)
c       s1      = sin(Pm)
c       s2      = sin(Qm)
c  The following lines are a clever way of
c  upgrading the sin and cosine for the next
c  frequency.  Essentally we are evaluating:
c       cos(f + df) = cos(f)cos(df) - sin(f)sin(df)
c       sin(f + df) = sin(f)cos(df) + cos(f)sin(df)
c
          cp = c01(m)
          cq = c02(m)
          tm = dcmplx(tt(m),0.d0)
          ag = aag(m)
          ag1 = aag1(m)
          x = xh(m)
          y = xh1(m)
c
c  These "a" terms those shown on page 21 of H53
c
          if (vs1d(m).ne.0.d0) then
            a11 = ag*cp - ag1*cq
            a12 = ag1*pd + ag*qm
            a13 = (cq - cp)/tm
            a14 = (pd + qm)/tm
            a21 = -ag*pm - ag1*qd
            a22 = ag*cq - ag1*cp
            a23 = (pm + qd)/tm
            a24 = a13
            a31 = tm*ag*ag1*(cp - cq)
            a32 = tm*(y*pd + x*qm)
            a33 = a22
            a34 = a12
            a41 = tm*(x*pm + y*qd)
            a42 = a31
            a43 = a21
            a44 = a11
          else
c---special forms for fluid layer from eqn 6.3 of H53.  Note there
c   seems to be a sign error on a12 in his paper, however.
            a11 = dcmplx(0.d0,0.d0)
            a12 = -pd
            a13 = -cp/tm
            a14 = dcmplx(0.d0,0.d0)
            a21 = dcmplx(0.d0,0.d0)
            a22 = cp
            a23 = pm/tm
            a24 = dcmplx(0.d0,0.d0)
            a31 = dcmplx(0.d0,0.d0)
            a32 = tm*pd
            a33 = cp
            a34 = dcmplx(0.d0,0.d0)
            a41 = dcmplx(0.d0,0.d0)
            a42 = dcmplx(0.d0,0.d0)
            a43 = dcmplx(0.d0,0.d0)
            a44 = dcmplx(0.d0,0.d0)
          endif
c
c       D is the accumulation of the multiplication of the A
c       matrices as we descend through the layers.  Note that
c       because we start at the surface, the stresses are 0 and
c       so the vector we multiply is (u, w, 0, 0).  This means 
c       that for the entire problem we need save only the first
c       two columns of the multiplication.
c
c
          d11 =  a11*b11 + a12*b21 + a13*b31 + a14*b41
          d12 =  a11*b12 + a12*b22 + a13*b32 + a14*b42
          d21 =  a21*b11 + a22*b21 + a23*b31 + a24*b41
          d22 =  a21*b12 + a22*b22 + a23*b32 + a24*b42
          d31 =  a31*b11 + a32*b21 + a33*b31 + a34*b41
          d32 =  a31*b12 + a32*b22 + a33*b32 + a34*b42
          d41 =  a41*b11 + a42*b21 + a43*b31 + a44*b41
          d42 =  a41*b12 + a42*b22 + a43*b32 + a44*b42
          b11 = d11
          b12 = d12
          b21 = d21
          b22 = d22
          b31 = d31
          b32 = d32
          b41 = d41
          b42 = d42

7       continue
c
c       Now multiply A by inverse E of the last layer.
c       C is then the same as J in H53.

        c11 = e11*b11 + e13*b31
        c12 = e11*b12 + e13*b32
        c21 = e22*b21 + e24*b41
        c22 = e22*b22 + e24*b42
        c31 = e31*b11 + e33*b31
        c32 = e31*b12 + e33*b32
        c41 = e42*b21 + e44*b41
        c42 = e42*b22 + e44*b42

c  ag id the denominator
        ag = (c11 - c21)*(c32-c42) - (c12 - c22)*(c31 - c41)

c---t for P waves
        if (ctype.eq.'P') then
          tm = dcmplx(2.d0*cv,0.d0)/(dcmplx(vp1d(nl),0.d0)*ag)
c---t for S waves
        else
          tm = dcmplx(cv,0.d0)/(dcmplx(vs1d(nl),0.d0)*ag)
        endif
c
c       wr, wi are the real and imaginary components of vertical
c               velocity
c       ur, ui are the real and imaginary components of horizontal
c               velocity
c       See equations 5 and 6 of H62.  Note that D = |D| e(i theta), so
c       1/D = e(-i theta)/|D| = (x - iy)/|D|^2 = D*/(DD*)
c
c
c       rs1 and rs2 are the real and imaginary parts of the
c       velocity ratio u/w for P and w/u for S.  Note that
c       u/w = uw*/(ww*) = (ur*wr + ui*wi + i(ui*wr - ur*wi))/ww*
c       w/u = wu*/(uu*) = (wr*ur + wi*ui + i(wi*ur - wr*ui))/uu*
c
c---P terms: The "W" terms are reversed in sign because a postive
c       "L" component will will be in the (-Z, +R) direction as defined by
c       H53.  Note that this does NOT change the sense of Z - it is still 
c       positive down.  
c
c       The idea is that the impulse is (sind, cosd) which would be a pulse in
c       the +X, +Z direction, but we want the response to (+x, -z) or (sind, -cosd)
c       so we revers the sign on the Z component.
c
        if (ctype.eq.'P') then
c--this from equations 10 and 11 of H62, modified as discussed at the top of page 4752
          wp = tm*(c31 - c41)
          up = tm*(c42 - c32)
        else
c---S terms: these are the Haskell conventions on W and U.  A source in the +H
c   direction would be (+Z, +R) in the Haskell sense, so we retain the signs
c   on W and U.  Equations 16 and 17 in H62.
          wp = tm*(c21 - c11)
          up = tm*(c12 - c22)
        endif

c---these are the surface displacements
c       write(*,*) ' wp, up = ', wp, up
c       read(*,*) pause

c---now propagate backwards, looping over the grid points
        ibovo  = 1
        ifin = 0

        ugrn1f(1:nzpts) = dcmplx(0.d0,0.d0)
        wgrn1f(1:nzpts) = dcmplx(0.d0,0.d0)

        if (lflip) then
           izb = nzpts 
           ize = 1
           idir =-1
        else
           izb = 1
           ize = nzpts
           idir = 1
        endif
!
!...... need to be careful, distributed models may not start at top
        if (dabs(depth(izb)- z1d(1)).lt.1.11d-7) then
           ugrn1f(izb) = up
           wgrn1f(izb) = wp
           if (lflip) izb = izb + idir
        endif

c---copy over response at the surface

        do iz =izb,ize,idir ! 2, nz
          !dept = dz*(iz-1) + zorig
          if (lflip) then
             dept = ztop - depth(iz)
          else
             dept = depth(iz)
          endif
          
c----Loop over the model, setting up layer specific values
c----dep(1) should the top of the model (e.g., zero) so start at dep(2)
          do i = 2,nl
            if (dept.lt.h1d(i)) go to 4
          enddo
          i = nl
4         ibov = i - 1
          dloc = dept - h1d(i-1)
          if (dloc.lt.0.d0) then
            WRITE(*,*) 'haskattn: Error in haskgrn dloc is < 0! ',dloc
            ierr = 1
            return
          endif
c         write(*,*) ' iz, dept, dloc, ibov = ', iz, dept, dloc, ibov

c----see if we have entered a new layer.  If so, finish off the previous one
          if (ibov.ne.ibovo) then
            ibovo = ibov
            if (ibov.gt.1) then
              ifin = ifin + 1
              hc = pic*h(ifin)
c----z01, z02 are Pm and Qm in H53
              z01 = dcmplx(hc,0.d0)*r1(ifin)
              z02 = dcmplx(hc,0.d0)*r2(ifin)

              x = cdexp(eye*z01)
              c01(ifin) = dcmplx(0.5d0,0.d0)*(x + dcmplx(1.d0,0.d0)/x)
              s01(ifin) = dcmplx(0.5d0,0.d0)*(x - dcmplx(1.d0,0.d0)/x)

              x = cdexp(eye*z02)
              c02(ifin) = dcmplx(0.5d0,0.d0)*(x + dcmplx(1.d0,0.d0)/x)
              s02(ifin) = dcmplx(0.5d0,0.d0)*(x - dcmplx(1.d0,0.d0)/x)
            end if
          end if

c----current layer terms
          hc = pic*dloc
c----z01, z02 are Pm and Qm in H53
          z01 = dcmplx(hc,0.d0)*r1(ibov)
          z02 = dcmplx(hc,0.d0)*r2(ibov)
          x = cdexp(eye*z01)
          c01(ibov) = dcmplx(0.5d0,0.d0)*(x + dcmplx(1.d0,0.d0)/x)
          s01(ibov) = dcmplx(0.5d0,0.d0)*(x - dcmplx(1.d0,0.d0)/x)

          x = cdexp(eye*z02)
          c02(ibov) = dcmplx(0.5d0,0.d0)*(x + dcmplx(1.d0,0.d0)/x)
          s02(ibov) = dcmplx(0.5d0,0.d0)*(x - dcmplx(1.d0,0.d0)/x)

c---The following should be equivalent to the inverse E matrix for
c         the half space (see Eqn 2.16 of H53).

          if (ibov.eq.nl) then
            vsvpnl = vs1d(nl)/vp1d(nl)
            vsvpnl = vsvpnl*vsvpnl
            coa = cv/vp1d(nl)
            e11 = -dcmplx(2.d0*vsvpnl,0.d0)
            e13 =  dcmplx(vsvpnl,0.d0)
            e22 =  dcmplx(coa*coa,0.d0)*aag1(nl)/r1(nl)
            e24 =  dcmplx(vsvpnl,0.d0)/r1(nl)
            e31 =  aag1(nl)/(aag(nl)*r2(nl))
            e33 = -dcmplx(1.d0,0.d0)/(dcmplx(2.d0,0.d0)*r2(nl))
            e42 =  dcmplx(1.d0,0.d0)
            e44 =  dcmplx(0.5d0,0.d0)
          endif


c----The 4 x 2 b matrix will hold the multiplication
c    of all the a matrices.  This initialization allows
c    the first two columns of the a matrix to be saved
c    the first time around.

c---this part you only need do when you enter a new layer.  Save
c   results up to ifin then then do one ore iteration for the current layer.

          b11 = dcmplx(1.d0,0.d0)
          b12 = dcmplx(0.d0,0.d0)
          b21 = dcmplx(0.d0,0.d0)
          b22 = dcmplx(1.d0,0.d0)
          b31 = dcmplx(0.d0,0.d0)
          b32 = dcmplx(0.d0,0.d0)
          b41 = dcmplx(0.d0,0.d0)
          b42 = dcmplx(0.d0,0.d0)

c----Loop over all layers
          do  m = 1,ibov
            pd = s01(m)/r1(m)
            qd = s02(m)/r2(m)
            pm = s01(m)*r1(m)
            qm = s02(m)*r2(m)
c
c       c1, cp  = cos(Pm)
c       c2, cq  = cos(Qm)
c       s1      = sin(Pm)
c       s2      = sin(Qm)
c  The following lines are a clever way of
c  upgrading the sin and cosine for the next
c  frequency.  Essentally we are evaluating:
c       cos(f + df) = cos(f)cos(df) - sin(f)sin(df)
c       sin(f + df) = sin(f)cos(df) + cos(f)sin(df)
c
            cp = c01(m)
            cq = c02(m)
            tm = dcmplx(tt(m),0.d0)
            ag = aag(m)
            ag1 = aag1(m)
            x = xh(m)
            y = xh1(m)
c
c  These "a" terms those shown on page 21 of H53
c
          if (vs1d(m).ne.0.d0) then
            a11 = ag*cp - ag1*cq
            a12 = ag1*pd + ag*qm
            a13 = (cq - cp)/tm
            a14 = (pd + qm)/tm
            a21 = -ag*pm - ag1*qd
            a22 = ag*cq - ag1*cp
            a23 = (pm + qd)/tm
            a24 = a13
            a31 = tm*ag*ag1*(cp - cq)
            a32 = tm*(y*pd + x*qm)
            a33 = a22
            a34 = a12
            a41 = tm*(x*pm + y*qd)
            a42 = a31
            a43 = a21
            a44 = a11
          else
c---special forms for fluid layer from eqn 6.3 of H53.  Note there
c   seems to be a sign error on a12 in his paper, however.
            a11 = dcmplx(0.d0,0.d0)
            a12 = -pd
            a13 = -cp/tm
            a14 = dcmplx(0.d0,0.d0)
            a21 = dcmplx(0.d0,0.d0)
            a22 = cp
            a23 = pm/tm
            a24 = dcmplx(0.d0,0.d0)
            a31 = dcmplx(0.d0,0.d0)
            a32 = tm*pd
            a33 = cp
            a34 = dcmplx(0.d0,0.d0)
            a41 = dcmplx(0.d0,0.d0)
            a42 = dcmplx(0.d0,0.d0)
            a43 = dcmplx(0.d0,0.d0)
            a44 = dcmplx(0.d0,0.d0)
          endif

c       D is the accumulation of the multiplication of the A
c       matrices as we descend through the layers.  Note that
c       because we start at the surface, the stresses are 0 and
c       so the vector we multiply is (u, w, 0, 0).  This means 
c       that for the entire problem we need save only the first
c       two columns of the multiplication.
c
c
            d11 =  a11*b11 + a12*b21 + a13*b31 + a14*b41
            d12 =  a11*b12 + a12*b22 + a13*b32 + a14*b42
            d21 =  a21*b11 + a22*b21 + a23*b31 + a24*b41
            d22 = +a21*b12 + a22*b22 + a23*b32 + a24*b42
            d31 =  a31*b11 + a32*b21 + a33*b31 + a34*b41
            d32 =  a31*b12 + a32*b22 + a33*b32 + a34*b42
            d41 =  a41*b11 + a42*b21 + a43*b31 + a44*b41
            d42 =  a41*b12 + a42*b22 + a43*b32 + a44*b42
            b11 = d11
            b12 = d12
            b21 = d21
            b22 = d22
            b31 = d31
            b32 = d32
            b41 = d41
            b42 = d42

          enddo

          !ugrn1f(iz) = b11*ugrn1f(1) + b12*wgrn1f(1)
          !wgrn1f(iz) = b21*ugrn1f(1) + b22*wgrn1f(1)
          ugrn1f(iz) = b11*up + b12*wp
          wgrn1f(iz) = b21*up + b22*wp

        enddo !end loop over nodes
!
!...... free space 
        deallocate(r1)
        deallocate(r2)
        deallocate(aag)
        deallocate(aag1)
        deallocate(xh) 
        deallocate(xh1) 
        deallocate(s01)
        deallocate(c01)
        deallocate(s02)
        deallocate(c02)
        deallocate(tt)
        return
        end

      COMPLEX*16 FUNCTION VCMPLX(FREQ,FREQ0,V,Q) 
!
!     Calculates the complex velocity given a reference frequency, 
!     velocity, and the reciprocal of the quality factor.  Based on 
!     SWRs dispersion.f and Aki and Richards pg 175 formula 5.94. 
!
!     INPUT     MEANING
!     -----     -------
!     FREQ      current frequency (Hz)
!     FREQ0     reference frequency (Hz)
!     V         velocity Vp or Vs (m/s)
!     Q         dispersion factor (Qp or Qs)
!
!     OUTPUT    MEANING
!     ------    -------  
!     VCMPLX    complex velocity
!
!.... variable declarations
      REAL*8, INTENT(IN) :: FREQ, FREQ0, V, Q
      REAL*8 FAC1, FACT, VRNEW, VINEW
      REAL*8 PI 
      PARAMETER(PI = 3.1415926535897932384626434D0) 
!
!----------------------------------------------------------------------!
!
!.... checks
      IF (FREQ0.EQ.0.D0 .OR. DABS(Q - 9999.D0).LT.1.D-5) THEN
         VCMPLX = DCMPLX(V,0.D0)
         RETURN
      ENDIF 
      IF (Q.EQ.0.D0) THEN
         VCMPLX = DCMPLX(V,0.D0)
         WRITE(*,*) 'vcmplx: 0 quality factor!, returning real velocity'
         RETURN
      ENDIF
      IF (FREQ.GT.FREQ0) 
     ;WRITE(*,*) 'vcmplx: Warning freq > freq0, velocity increase!'
      FAC1 = DLOG(FREQ/FREQ0) 
      FACT = 1.D0 + 1.D0/(PI*Q)*FAC1
      IF (FACT.LT.0.1D0) FACT = 0.1D0 !a 90 percent reduction in vel
      VRNEW = V*FACT     !Re{ Aki and Richards 5.94}
      VINEW =-V/(2.D0*Q) !Im{ Aki And Richards 5.94}
      VCMPLX = DCMPLX(VRNEW,VINEW) 
      RETURN
      END 
