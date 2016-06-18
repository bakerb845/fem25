
        subroutine haskinf(nl,nzpts,srctyp,lflip, df,ain, z1d,vp1d,
     ;                     vs1d,rh1d,depth, ugrn1f,wgrn1f, ierr)

c       Arguments:
c       df      The frequency
c       isrc    The source number
c
c       Note: this version tested for P in multiple layers but NOT for  S in any case.
c

        implicit none
        character(2), intent(in) :: srctyp
        real*8, intent(in) :: z1d(nl),vp1d(nl),vs1d(nl),rh1d(nl),
     ;                        depth(nzpts),ain,df
        integer*4, intent(in) :: nl,nzpts
        logical*4, intent(in) :: lflip
        complex*16, intent(out) :: ugrn1f(nzpts),wgrn1f(nzpts)
        integer*4, intent(out) :: ierr
!...... local variables
        complex*16, allocatable :: ra(:), rb(:), cosp(:), sinp(:), 
     ;                             cosq(:), sinq(:)  
        real*8, allocatable :: g(:), gm1(:), h(:), h1d(:)
        complex*16 a(4,4), b(4,4), ce(4,4), d(4,4), e(4,4)
        complex*16 u, w, u1, w1, sig, taw, pm, qm,
     ;             bij, cij, dij, den, du1, wu1,
     ;             j11, j12, j21, j22, j31, j32, j41, j42, ram, rbm, 
     ;             cospm, sinpm, cosqm, sinqm 
        real*8 trans(4) 
        real*8 x, gm, gmm1, tm, pic, hc, pcsq, cv, csq, 
     ;         cova, covb, dloc, amsq, bmsq, dept, rhom, ztop,  
     ;         pi, pi180  
        integer*4 i, j, k, m, iz, iz1, izb, ize, idir, ibov,  
     ;            ibovo, ifin
        character(1) ctype
        parameter(pi180 = 0.017453292519943295d0,
     ;            pi = 3.1415926535897932384626434d0)

c----transform array for 1st interface (u, w, sig, taw)
c----We reverse the sign on w because we go from +z to -z.   The signs
c    on the stresses are difference because the continuous traction condition depends
c    on the the sign of the normal to the interace, which also reverses. Hence, sig(zz) 
c    has a double negative, while sig(xz) = tau has a single negative.
        ierr = 0
        trans(1) =  1.d0
        trans(2) = -1.d0
        trans(3) =  1.d0
        trans(4) = -1.d0

c----incident wave type
        ctype = 'P'
        if (srctyp.eq.'PS' .or. srctyp.eq.'ps') ctype = 'S'


c----apparent p wavespeed (1/px)
        cv = vp1d(nl)/dsin(ain*pi180)
        csq = cv*cv

c---apparent wavenumber (kx)
        pic = 2.d0*pi*df/cv
c       pic = df/cv


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
! 
!...... thickness array where layer 1 is the model basement
        do i = 1, nl-1
          h(i) = h1d(i+1) - h1d(i)
        enddo
!
!...... set space
        allocate(ra(nl))
        allocate(rb(nl)) 
        allocate(g(nl))
        allocate(cosp(nl))
        allocate(sinp(nl))
        allocate(cosq(nl))
        allocate(sinq(nl))
        allocate(gm1(nl))
c---Loop over the model, setting up layer specific values
        do m = 1, nl

c--- c/alpha and c/beta
          cova = cv/vp1d(m)
c---set up a bogus value for covb in case of vs = 0
          if (vs1d(m).ne.0.d0) then
            covb = cv/vs1d(m)
          else
            covb = 1.d0
          endif

c----r1 and r2 are r(alpha) and r(beta)
          x = cova*cova - 1.d0
          if (x.gt.0.d0) then
            ra(m) = dsqrt(x)
          else
            ra(m) = -dcmplx(0.d0, 1.d0)*dsqrt(-x)
          endif
          x = covb*covb - 1.d0
          if (x.gt.0.d0) then
            rb(m) = dsqrt(x)
          else
            rb(m) = -dcmplx(0.d0, 1.d0)*dsqrt(-x)
          endif
          g(m)   = 2.d0/(covb*covb)
          gm1(m) = g(m) - 1.d0

          if(m.lt.nl) then
            hc = pic*h(m)
            pm = hc*ra(m)
            qm = hc*rb(m)

c---NB; if pm or qm is complex, cos turns into cosh so we do not have to 
c   to anything special
            cosp(m) = cdcos(pm)
            sinp(m) = cdsin(pm)
            cosq(m) = cdcos(qm)
            sinq(m) = cdsin(qm)
          endif
        enddo

        amsq = vp1d(nl)*vp1d(nl)
        bmsq = vs1d(nl)*vs1d(nl)
        ram  = ra(nl)
        rbm  = rb(nl)
        gm   = g(nl)
        gmm1 = gm - 1.d0
        rhom = rh1d(nl)
        pcsq = rhom*csq

c---- Inverse E matrix for last layer
        e(1,1) = -2.d0*(bmsq/amsq)
        e(1,2) =  0.d0
        e(1,3) =  1.d0/(rhom*amsq)
        e(1,4) =  0.d0
        e(2,1) =  0.d0
        e(2,2) =  csq*gmm1/(amsq*ram)
        e(2,3) =  0.d0
        e(2,4) =  1.d0/(rhom*amsq*ram)
        e(3,1) =  gmm1/(gm*rbm)
        e(3,2) =  0.d0
        e(3,3) =  -1.d0/(pcsq*gm*rbm)
        e(3,4) =  0.d0
        e(4,1) =  0.d0
        e(4,2) =  1.d0
        e(4,3) =  0.d0
        e(4,4) =  1.d0/(pcsq*gm)

c---initialize product matrix b
        do i = 1, 4
          do j = 1, 4
            b(j,i) = dcmplx(0.d0, 0.d0)
          enddo
        enddo
        do i = 1, 4
          b(i,i) = dcmplx(1.d0, 0.d0)
        enddo

c---Loop over all layers below the first interface to form a(n-1)a(n-2)...a(3)a(2)
        do m = 2, nl-1
          gm   = g(m)
          gmm1 = g(m) - 1.d0
          cospm = cosp(m)
          sinpm = sinp(m)
          cosqm = cosq(m)
          sinqm = sinq(m)
          ram = ra(m)
          rbm = rb(m)
          pcsq = rh1d(m)*csq

c-- a matrix
          a(1,1) =  gm*cospm - gmm1*cosqm
          a(1,2) =  dcmplx(0.d0,1.d0)*(gmm1*sinpm/ram + gm*rbm*sinqm)
          a(1,3) = -(cospm - cosqm)/pcsq
          a(1,4) =  dcmplx(0.d0,1.d0)*(sinpm/ram + rbm*sinqm)/pcsq
          a(2,1) = -dcmplx(0.d0,1.d0)*(gm*ram*sinpm + gmm1*sinqm/rbm)
          a(2,2) = -gmm1*cospm + gm*cosqm
          a(2,3) =  dcmplx(0.d0,1.d0)*(ram*sinpm + sinqm/rbm)/pcsq
          a(2,4) =  a(1,3)
          a(3,1) =  pcsq*gm*gmm1*(cospm - cosqm)
          a(3,2) =  dcmplx(0.d0,1.d0)*pcsq*(gmm1*gmm1*sinpm/ram 
     ;           + gm*gm*rbm*sinqm)
          a(3,3) =  a(2,2)
          a(3,4) =  a(1,2)
          a(4,1) =  dcmplx(0.d0,1.d0)*pcsq*(gm*gm*ram*sinpm 
     ;           + gmm1*gmm1*sinqm/rbm)
          a(4,2) =  a(3,1)
          a(4,3) =  a(2,1)
          a(4,4) =  a(1,1)

          do i = 1, 4
            do j = 1, 4
              dij = 0.0
              do k = 1, 4
                dij = dij + a(i,k)*b(k,j)
              enddo
              d(i,j) = dij
            enddo
          enddo
          do i = 1, 4
            do j = 1, 4
              b(i,j) = d(i,j)
            enddo
          enddo

        enddo

c---- Now premultiply by inverse E of the last layer
        do i = 1, 4
          do j = 1, 4
            cij = dcmplx(0.d0,0.d0)
            do k = 1, 4
              cij = cij + e(i,k)*b(k,j)
            enddo
            ce(i,j) = cij
          enddo
        enddo

        amsq = vp1d(1)*vp1d(1)
        bmsq = vs1d(1)*vs1d(1)
        ram  =  ra(1)
        rbm  =  rb(1)
        gm   =  g(1)
        gmm1 =  gm - 1.d0
        rhom =  rh1d(1)

c---- E matrix for 1st layer
        e(1,1) = -amsq/csq
        e(1,2) =  dcmplx(0.d0,0.d0)
        e(1,3) = -gm*rbm
        e(1,4) =  dcmplx(0.d0,0.d0)
        e(2,1) =  dcmplx(0.d0,0.d0)
        e(2,2) = -amsq*ram/csq
        e(2,3) =  dcmplx(0.d0,0.d0)
        e(2,4) =  gm
        e(3,1) = -rhom*amsq*gmm1
        e(3,2) =  dcmplx(0.d0,0.d0)
c---error in Haskell
        e(3,3) = -rhom*csq*gm*gm*rbm
        e(3,4) =  dcmplx(0.d0,0.d0)
        e(4,1) =  dcmplx(0.d0,0.d0)
        e(4,2) =  rhom*amsq*gm*ram
        e(4,3) =  dcmplx(0.d0,0.d0)
        e(4,4) = -rhom*csq*gm*gmm1

c---transform coordinates
        do i = 1, 4
          do j = 1, 4
            e(i,j) = e(i,j)*trans(i)
          enddo
        enddo

c---post multiply by E of first layer; the result is the "K" matrix in the notes
        do i = 1, 4
          do j = 1, 4
            bij = 0.0
            do k = 1, 4
              bij = bij + ce(i,k)*e(k,j)
            enddo
            b(i,j) = bij
          enddo
        enddo

c---now form the "J" terms

        j11 = b(1,1) + b(1,2)
        j21 = b(2,1) + b(2,2)
        j31 = b(3,1) + b(3,2)
        j41 = b(4,1) + b(4,2)
        j12 = b(1,3) + b(1,4)
        j22 = b(2,3) + b(2,4)
        j32 = b(3,3) + b(3,4)
        j42 = b(4,3) + b(4,4)

c----denominator

        den = (j11 - j21)*(j32 - j42) - (j12 - j22)*(j31 - j41)

c--- these are the potenials in layer 1.   Note that the factor of "2" in the notes is
c    included in the "tm" term (need to check for S waves)

        du1 = (j32 - j42)/den
        wu1 = (j41 - j31)/den

c---factor to scale to displacement
c---t for P waves
        if (ctype.eq.'P') then
          tm = 2.d0*cv/vp1d(nl)
c---t for S waves
        else
          tm = cv/vs1d(nl)
        endif

c---now propagate backwards, looping over the grid points
c
c       Note that in the first layer, we multiply the potentials in layer 1 by the appropriate
c       "D" to recover velocities at any particular depth.  When we get to the first interface,
c       we can recover them in the usual way.
c
        amsq = vp1d(1)*vp1d(1)
        bmsq = vs1d(1)*vs1d(1)
        ram  = ra(1)
        rbm  = rb(1)
        gm   = g(1)
        gmm1 = gm - 1.d0
        rhom = rh1d(1)

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
        ugrn1f(1:nzpts) = dcmplx(0.d0,0.d0)
        wgrn1f(1:nzpts) = dcmplx(0.d0,0.d0)
!       if (dabs(depth(izb)- z1d(1)).lt.1.11d-7) then
!          ugrn1f(izb) = quresp
!          wgrn1f(izb) = qwresp
!          if (lflip) izb = izb + idir
!       endif

        do iz = izb,ize,idir !1, nz
          !dept = dz*(iz-1) + zorig
          if (lflip) then
             dept = ztop - depth(iz)
          else
             dept = depth(iz)
          endif
          if (dept.gt.h1d(2)) go to 44
          dloc = h1d(2) - dept
          if (dloc.lt.0.d0) then
            write(*,*) 'haskinf: Error dloc is < 0! '
            ierr = 1
            return 
          endif
          hc = pic*dloc
c----Pm and Qm in H53
          pm = hc*ra(1)
          qm = hc*rb(1)
          cospm = cdcos(pm)
          sinpm = cdsin(pm)
          cosqm = cdcos(qm)
          sinqm = cdsin(qm)
c---- D matrix for this depth
          d(1,1) = -(amsq/csq)*cospm
          d(1,2) =  dcmplx(0.d0, 1.d0)*(amsq/csq)*sinpm
          d(1,3) = -gm*rbm*cosqm
          d(1,4) =  dcmplx(0.d0, 1.d0)*gm*rbm*sinqm
          d(2,1) =  dcmplx(0.d0, 1.d0)*(amsq/csq)*ram*sinpm
          d(2,2) = -(amsq/csq)*ram*cospm
          d(2,3) = -dcmplx(0.d0, 1.d0)*gm*sinqm
          d(2,4) =  gm*cosqm
          d(3,1) = -rhom*amsq*gmm1*cospm 
          d(3,2) =  dcmplx(0.d0, 1.d0)*rhom*amsq*gmm1*sinpm
          d(3,3) = -rhom*csq*gm*gm*rbm*cospm 
          d(3,4) =  dcmplx(0.d0, 1.d0)*rhom*csq*gm*gm*rbm*sinqm
          d(4,1) = -dcmplx(0.d0, 1.d0)*rhom*amsq*gm*ram*sinpm
          d(4,2) =  rhom*amsq*gm*ram*cospm 
          d(4,3) =  dcmplx(0.d0, 1.d0)*rhom*csq*gm*gmm1*sinqm
          d(4,4) = -rhom*csq*gm*gmm1*cosqm 
c---velocities at this grid point
          u = du1*(d(1,1) + d(1,2)) + wu1*(d(1,3) + d(1,4))
          w = du1*(d(2,1) + d(2,2)) + wu1*(d(2,3) + d(2,4))
c---Note we retain the sign on w because it should be reversed anyway (it gets negated below as well). The
c   idea is that a positive longitudinal wave will have a (-Z, +R) displacement.
c   The sign reversal on u is in accordance with my notes on how to convert potential to velocity.
          ugrn1f(iz) = -u*tm
          wgrn1f(iz) =  w*tm 
c         write (*,*) ' iz, u, w, tm = ', iz, u, w, tm
c         read(*,*) pause
        end do

c---finished with first layer; compute the response at interface 1
44      iz1 = iz
c---- E matrix ; note this is just D above at h = 0

        e(1,1) = -amsq/csq
        e(1,2) =  dcmplx(0.d0,0.d0)
        e(1,3) = -gm*rbm
        e(1,4) =  dcmplx(0.d0,0.d0)
        e(2,1) =  dcmplx(0.d0,0.d0)
        e(2,2) = -amsq*ram/csq
        e(2,3) =  dcmplx(0.d0,0.d0)
        e(2,4) =  gm
        e(3,1) = -rhom*amsq*gmm1
c       e(3,1) = -rhom*amsq*gmm1/csq
        e(3,2) =  dcmplx(0.d0,0.d0)
c---Note there is a mistake in Haskell 1953, eqn 2.12 for this term
        e(3,3) = -rhom*csq*gm*gm*rbm
c       e(3,3) = -rhom*gm*gm*rbm
        e(3,4) =  dcmplx(0.d0,0.d0)
        e(4,1) =  dcmplx(0.d0,0.d0)
        e(4,2) =  rhom*amsq*gm*ram
c       e(4,2) =  rhom*amsq*gm*ram/csq
        e(4,3) =  dcmplx(0.d0,0.d0)
        e(4,4) = -rhom*csq*gm*gmm1
c       e(4,4) = -rhom*gm*gmm1

        u1  = du1*(e(1,1) + e(1,2)) + wu1*(e(1,3) + e(1,4))
        w1  = du1*(e(2,1) + e(2,2)) + wu1*(e(2,3) + e(2,4))
        sig = du1*(e(3,1) + e(3,2)) + wu1*(e(3,3) + e(3,4))
        taw = du1*(e(4,1) + e(4,2)) + wu1*(e(4,3) + e(4,4))

c---this part applies the trans array to the displacements
        w1 = -w1
        taw = -taw

        ibovo  = 1
        ifin = 0
        ibovo  = 2
        ifin = 1

        do iz =iz1,ize,idir ! iz1, nz
          !dept = dz*(iz-1) + zorig
          if (lflip) then
             dept = ztop - depth(iz)
          else
             dept = depth(iz)
          endif

c----Loop over the model, setting up layer specific values
c
c       Note that in the first layer, we multiply the potentials in layer 1 by the appropriate
c       "D" to recover velocities at any particular depth.  When we get to the first interface,
c       we can recover them in the usual way.
c
c----dep(1) should the top of the model (e.g., zero) so start at dep(2)
          do i = 2,nl
            if (dept.lt.h1d(i)) go to 4
          enddo
          i = nl
4         ibov = i - 1
          dloc = dept - h1d(i-1)
          if (dloc.lt.0.d0) then
            write(*,*) 'haskinf Error: dloc is < 0! '
            ierr = 1 
            return 
          endif
c         write(*,*) ' iz, dept, dloc, ibov, ibovo = ', iz, dept, dloc, ibov, ibovo
c         read(*,*) pause

c----see if we have entered a new layer.  If so, finish off the previous one
          if (ibov.ne.ibovo) then
            ibovo = ibov
            if (ibov.gt.1) then
              ifin = ifin + 1
              hc = pic*h(ifin)
c----Pm and Qm in H53
              pm = hc*ra(ifin)
              qm = hc*rb(ifin)
              cosp(ifin) = cdcos(pm)
              sinp(ifin) = cdsin(pm)
              cosq(ifin) = cdcos(qm)
              sinq(ifin) = cdsin(qm)
            end if
          end if

c----current layer terms
          hc = pic*dloc
c----z01, z02 are Pm and Qm in H53
          pm = hc*ra(ibov)
          qm = hc*rb(ibov)
c----Case of r(alpha) being imaginary
          cosp(ibov) = cdcos(pm)
          sinp(ibov) = cdsin(pm)
          cosq(ibov) = cdcos(qm)
          sinq(ibov) = cdsin(qm)

          if (ibov.eq.nl) then
            amsq = vp1d(nl)*vp1d(nl)
            bmsq = vs1d(nl)*vs1d(nl)
            rbm  = rb(nl)
            ram  = ra(nl)
            rbm  = rb(nl)
            gm   = g(nl)
            gmm1 = gm - 1
            rhom = rh1d(nl)
            pcsq = rhom*csq

c---- Inverse E matrix
            e(1,1) = -2.d0*(bmsq/amsq)
            e(1,2) =  dcmplx(0.d0,0.d0)
            e(1,3) =  1.d0/(rhom*amsq)
            e(1,4) =  dcmplx(0.d0,0.d0)
            e(2,1) =  dcmplx(0.d0,0.d0)
            e(2,2) =  csq*gmm1/(amsq*ram)
            e(2,3) =  dcmplx(0.d0,0.d0)
            e(2,4) =  1.d0/(rhom*amsq*ram)
            e(3,1) =  gmm1/(gm*rbm)
            e(3,2) =  dcmplx(0.d0,0.d0)
            e(3,3) =  -1.d0/(pcsq*gm*rbm)
            e(3,4) =  dcmplx(0.d0,0.d0)
            e(4,1) =  dcmplx(0.d0,0.d0)
            e(4,2) =  dcmplx(1.d0,0.d0)
            e(4,3) =  dcmplx(0.d0,0.d0)
            e(4,4) =  1.d0/(pcsq*gm)
          endif

c---initialize product matrix b
          do i = 1, 4
            do j = 1, 4
              b(j,i) = dcmplx(0.d0,0.d0)
            enddo
          enddo
          do i = 1, 4
            b(i,i) = dcmplx(1.d0,0.d0)
          enddo

c---Loop over all layers below the first interface
c         do m = 1, ibov
          do m = 2, ibov
            gm   = g(m)
            gmm1 = g(m) - 1.d0
            cospm = cosp(m)
            sinpm = sinp(m)
            cosqm = cosq(m)
            sinqm = sinq(m)
            ram = ra(m)
            rbm = rb(m)
            pcsq = rh1d(m)*csq
            rhom = rh1d(m)

c-- a matrix
            a(1,1) =  gm*cospm - gmm1*cosqm
            a(1,2) =  dcmplx(0.d0,1.d0)*(gmm1*sinpm/ram + gm*rbm*sinqm)
            a(1,3) = -(cospm - cosqm)/pcsq
c           a(1,3) = -(cospm - cosqm)/rhom
            a(1,4) =  dcmplx(0.d0,1.d0)*(sinpm/ram + rbm*sinqm)/pcsq
c           a(1,4) =  cmplx(0.,1.)*(sinpm/ram + rbm*sinqm)/rhom
            a(2,1) = -dcmplx(0.d0,1.d0)*(gm*ram*sinpm + gmm1*sinqm/rbm)
            a(2,2) = -gmm1*cospm + gm*cosqm
            a(2,3) =  dcmplx(0.d0,1.d0)*(ram*sinpm + sinqm/rbm)/pcsq
c           a(2,3) =  cmplx(0.,1.)*(ram*sinpm + sinqm/rbm)/rhom
            a(2,4) =  a(1,3)
            a(3,1) =  pcsq*gm*gmm1*(cospm - cosqm)
            a(3,2) =  dcmplx(0.d0,1.d0)*pcsq*(gmm1*gmm1*sinpm/ram 
     ;             + gm*gm*rbm*sinqm)
            a(3,3) =  a(2,2)
            a(3,4) =  a(1,2)
            a(4,1) =  dcmplx(0.d0,1.d0)*pcsq*(gm*gm*ram*sinpm 
     ;             + gmm1*gmm1*sinqm/rbm)
            a(4,2) =  a(3,1)
            a(4,3) =  a(2,1)
            a(4,4) =  a(1,1)

            do i = 1, 4
              do j = 1, 4
                dij = dcmplx(0.d0,0.d0)
                do k = 1, 4
                dij = dij + a(i,k)*b(k,j)
                enddo
                d(i,j) = dij
              enddo
            enddo
            do i = 1, 4
              do j = 1, 4
                b(i,j) = d(i,j)
              enddo
            enddo

          enddo

          u = u1*b(1,1) + w1*b(1,2) + sig*b(1,3) + taw*b(1,4)
          w = u1*b(2,1) + w1*b(2,2) + sig*b(2,3) + taw*b(2,4)
        
c---The sign on w is reversed because a positive longitudinal wave will have a (-Z, +R) displacement.
c   As above, the sign reversal on u is in accordance with my notes on how to convert potential to velocity.
          ugrn1f(iz) = -u*tm
          wgrn1f(iz) = -w*tm

        enddo
c---end loop over nodes
        if (allocated(h))    deallocate(h)
        if (allocated(h1d))  deallocate(h1d)
        if (allocated(ra))   deallocate(ra)
        if (allocated(rb))   deallocate(rb) 
        if (allocated(g))    deallocate(g)
        if (allocated(cosp)) deallocate(cosp)
        if (allocated(sinp)) deallocate(sinp)
        if (allocated(cosq)) deallocate(cosq)
        if (allocated(sinq)) deallocate(sinq)
        if (allocated(gm1))  deallocate(gm1) 
        return
        end



