**********************************************************************
      subroutine bcond_lin(a,c, za1,za0,zc1,zc0, ga,gc)
      real*8 a, c, za1, za0, zc1, zc0, ga, gc
      real*8 q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
c  --------------------------------------------------------------
c  Bourdary conditions (all parameter in common /bc/)
c         should be specfied in the subroutine
c  [a,c] : two coordinates of boundary points
c         za1 u'(a) + za0 u(a) = ga
c         zc1 u'(c) + zc0 u(c) = gc
c  G(x,t) = gr(x) gl(t) / wron  or  gr(t) gl(x) / wron,
c      where wron(zl,zr) = zl(x) zrp(x) - zlp(x) zr(x) = constant
c  ---------------------------
c  How to choose Green's functions:
c
c  If za0 * zc0 * (c-a) + ( zc1*za0 - za1*zc0 ) > 1e-4 then use
c    # linear green's functions (q0=0) and linear inhomogeneous term
c    gl = gl1*x + gl0
c    gr = gr1*x + gr0
c    ui = ui1*x + ui0
c
c  otherwise (for Neumann type problems)
c    # cosh,sinh for green's functions (q0=-1) and quadratic u_i(x)
c    gl = gl1 * dcosh(x) - gl0 * dsinh(x)
c    gr = gr1 * dcosh(x) - gr0 * dsinh(x)
c
c    ui = ui2*x**2 + ui1*x + ui0 where
c    ( 2 a za1 + a^2 za0 , za1 + a za0 , za0 ) (ui2)   (ga)
c    (                                       ) (ui1) = (  )
c    ( 2 c za1 + c^2 zc0 , zc1 + c zc0 , zc0 ) (ui0)   (gc)
c
c    - It solves (ui2,ui0) first to get (ui0) and then
c         solves (ui2,ui1) with (ga - za0 ui0, gc - zc0 ui0)
c    - This will be performed only when | det(ui1,ui0) | < 1d-4
c         i.e. (za1 + a za0, zc1 + c zc0) // (za0, zc0)
c      So det(ui2,ui1)=0 implies that
c         (2a za1 + a^2 za0, 2c za1 + c^2 zc0 ) // (za0, zc0)
c         and the rank of the linear system is ONE.
c      Therefore it is solvable only when (ga, gc) is in the range space.
c
c  TO DO -------------------------------
c  1. BETTER WAY TO FIND TO LINEARLY INDEPENDENT GREEN'S FUNCTION
c  2. Using a SVD routine to find a best possible quadratic inhomogeneous term
**********************************************************************
      real*8 det

      wron = za0 * zc0 * (c-a) + ( zc1*za0 - za1*zc0 )
      if ( abs(wron) .ge. 1d-4 ) then
         q0 = 0d0
         gl1 = za0
         gr1 = zc0
         gl0 = -a * za0 - za1
         gr0 = -c * zc0 - zc1

         ui2 = 0d0
         ui1 = ( gc*za0 - ga*zc0 ) / wron
         ui0 = ( (c*zc0+zc1)*ga - (a*za0+za1)*gc ) / wron
         return
      end if

      wron = (za0*zc0-za1*zc1)*dsinh(c-a) + (zc1*za0-za1*zc0)*dcosh(c-a)
      if ( abs(wron) .ge. 1d-4 ) then
         q0 = -1d0
         gl1 = za1*dcosh(a) + za0*dsinh(a)
         gl0 = za1*dsinh(a) + za0*dcosh(a)
         gr1 = zc1*dcosh(c) + zc0*dsinh(c)
         gr0 = zc1*dsinh(c) + zc0*dcosh(c)

         det = za0*zc0*(a*a-c*c) + ( 2d0*a*za1*zc0 - 2d0*c*zc1*za0 )
         if ( abs(det) .lt. 1d-4 ) then
            ui0 = 0d0
         else
            ui0 = ((a*a*za0+2d0*a*za1)*gc-(c*c*zc0+2d0*c*zc1)*ga) / det
         end if
         det = a*c*(a-c)*za0*zc0 + 2d0*(a-c)*za1*zc1
     &       + c*(2d0*a-c)*za1*zc0 + a*(a-2d0*c)*za0*zc1
         ui2 = (  (zc1+c*zc0)*(ga-ui0*za0)
     &          - (za1+a*za0)*(gc-ui0*zc0) ) / det
         ui1 = (  (2d0*a*za1+a*a*za0)*(gc-ui0*zc0)
     &          - (2d0*c*zc1+c*c*zc0)*(ga-ui0*za0) ) / det
         return
      end if

      print*, 'Given boundary conditions are Neumann Type', wron
      print*, 'This mixed boundary condition make some troubles'
      print*, 'to form two linearly independent Green''s functions'
      print*, 'using q0 = 0 or -1. Please make a new choice'
      stop

      end

**********************************************************************
c   gl(x), gr(x), gl'(x), gr'(x) :
c   Green's functions and it's derivatives
c   will be used for function evaluation (int_density.f)
**********************************************************************
      real*8 function glp_x(x)
      real*8 x, q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      if ( q0 .eq. 0 ) then
         glp_x = gl1
      else
         glp_x = gl1 * dsinh(x) - gl0 * dcosh(x)
      end if
      end

      real*8 function grp_x(x)
      real*8 x, q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      if ( q0 .eq. 0 ) then
         grp_x = gr1
      else
         grp_x = gr1 * dsinh(x) - gr0 * dcosh(x)
      end if
      end

      real*8 function gl_x(x)
      real*8 x, q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      if ( q0 .eq. 0 ) then
         gl_x = gl1*x + gl0
      else
         gl_x = gl1 * dcosh(x) - gl0 * dsinh(x)
      end if
      end

      real*8 function gr_x(x)
      real*8 x, q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      if ( q0 .eq. 0 ) then
         gr_x = gr1*x + gr0
      else
         gr_x = gr1 * dcosh(x) - gr0 * dsinh(x)
      end if
      end

**************************************************************************
      subroutine int_eqn(j,nnd,a,c,xnd,pt,qt,ft,gl,gr,phil,phir,eta)
      integer j, nnd, l
      real*8 x, a, c, xnd(*), pt(*), qt(*), ft(*)
      real*8 gl(*), gr(*), phil(*), phir(*), eta(*)
      real*8 q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
c  -----------------------------------------------------------------------
c  This routine defines a local integral equation from ODE
c  int_eqn: with inhomogeneous term, int_eig: without inhomogeneous term
c
c  all arguments except j and nnd are defind between 1+nnd*(j-1) and nnd*j
c  xnd, pt, qt, ft is a mandatory input
c  gl, gr, eta, phil, phir are mandatory outputs
c  -----------------------------------------------------------------------
c    ur(x) * wron = vl(x) = gl1 x + gl0 or gl1 * dcosh(x) - gl0 * dsinh(x)
c    ul(x) * wron = vr(x) = gr1 x + gr0 or gr1 * dcosh(x) - gr0 * dsinh(x)
**************************************************************************
      if ( q0 .eq. 0 ) then
         do l = 1 + nnd*(j-1), nnd*j
            x = xnd(l)
            gl(l) = gl1*x + gl0
            gr(l) = gr1*x + gr0
            phil(l) = ( pt(l) * gr1 + qt(l) * gr(l) ) / wron
            phir(l) = ( pt(l) * gl1 + qt(l) * gl(l) ) / wron
            eta(l) = ft(l) - ( pt(l)*ui1 + qt(l)*(ui1*x+ui0) )
         end do
      else
         do l = 1 + nnd*(j-1), nnd*j
            x = xnd(l)
            gl(l) = gl1 * dcosh(x) - gl0 * dsinh(x)
            gr(l) = gr1 * dcosh(x) - gr0 * dsinh(x)
            phil(l) = ( pt(l) * ( gl1*dsinh(x) - gl0*dcosh(x) )
     &                  + (qt(l)-q0) * gr(l) ) / wron
            phir(l) = ( pt(l) * ( gr1*dsinh(x) - gr0*dcosh(x) )
     &                  + (qt(l)-q0) * gl(l) ) / wron
            eta(l) = ft(l) - ( ui2 + pt(l)*(2d0*ui2*x+ui1)
     &                             + qt(l)*(ui2*x**2+ui1*x+ui0) )
         end do
      end if

      end

**************************************************************************
      subroutine int_eig(j,nnd,a,c,xnd,pt,qt,wcopy,gl,gr,psil,psir)
      integer j, nnd, l
      real*8 x, a, c, xnd(*), pt(*), qt(*)
      real*8 wcopy, gl(*), gr(*), psil(*), psir(*)
      real*8 q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
c  -----------------------------------------------------------------------
c  This routine defines a local integral equation from ODE
c  int_eqn: with inhomogeneous term, int_eig: without inhomogeneous term
c
c  all arguments except j and nnd are defind between 1+nnd*(j-1) and nnd*j
c  xnd,pt, qt, ft is a mandatory input
c  wcopy = a copy of wron for shifted eigenfunction iteration.
c  gl, gr, eta, phil, phir are mandatory outputs
c  -----------------------------------------------------------------------
c    ur(x) * wron = vl(x) = gl1 x + gl0 or gl1 * dcosh(x) - gl0 * dsinh(x)
c    ul(x) * wron = vr(x) = gr1 x + gr0 or gr1 * dcosh(x) - gr0 * dsinh(x)
**************************************************************************
      wcopy = wron
      if ( q0 .eq. 0 ) then
         do l = 1 + nnd*(j-1), nnd*j
            x = xnd(l)
            gl(l) = gl1*x + gl0
            gr(l) = gr1*x + gr0
            psil(l) = ( pt(l) * gr1 + qt(l) * gr(l) ) / wron
            psir(l) = ( pt(l) * gl1 + qt(l) * gl(l) ) / wron
         end do
      else
         do l = 1 + nnd*(j-1), nnd*j
            x = xnd(l)
            gl(l) = gl1 * dcosh(x) - gl0 * dsinh(x)
            gr(l) = gr1 * dcosh(x) - gr0 * dsinh(x)
            psil(l) = ( pt(l) * ( gl1*dsinh(x) - gl0*dcosh(x) )
     &                  + (qt(l)-q0) * gr(l) ) / wron
            psil(l) = ( pt(l) * ( gr1*dsinh(x) - gr0*dcosh(x) )
     &                  + (qt(l)-q0) * gl(l) ) / wron
         end do
      end if

      end

**********************************************************************
      subroutine initmesh(meshtype, br, nsub, leafid, a, c, b, blength)
      integer meshtype, nsub, nmxsub
      parameter (nmxsub=1024)
      integer leafid(nmxsub)
      real*8 br, a, c, b(nmxsub), blength(nmxsub)
c -----------------------------------------
c  Mesh type:
c    0 - Equal space
c    1,2 - Polynomial order 1 and 2
c    -1(br) - centered, br = b(center+1)/b(center)
c    -2(br) - 0 divides the intervals (nsub=even) : br = b(center+1)/b(center)
**********************************************************************
      integer kp, j
      real*8 bsum

      do kp = 1, nsub
         j = leafid(kp)
         if ( meshtype .eq. 0 ) then
            b(j) = ( dble(kp-1)*c + dble(nsub-kp+1)* a ) / dble(nsub)
            blength(j) = (c-a) / dble(nsub)
         else if (meshtype .eq. 1 ) then
            b(j) = dble(kp*(kp-1)) / (nsub*(nsub+1))
            blength(j) = 2d0 * kp / (nsub*(nsub+1))
         else if (meshtype .eq. 2 ) then
            b(j) = dble((kp-1)*kp*(2*kp-1)) / (nsub*(nsub+1)*(2*nsub+1))
            blength(j) = 6d0 * kp**2 / (nsub*(nsub+1)*(2*nsub+1))
         else if ( meshtype.eq.-1 .and. nsub/2*2.eq.nsub ) then
            bsum = 2d0 * br * (br**nsub-1d0) / (br**2-1d0)
            blength(j) = (c-a)/bsum * br**iabs(2*kp-nsub-1)
            if ( kp .le. nsub/2 ) then
               b(j) = (c+a)/2d0 - (c-a)/bsum
     &               * (br**(nsub+3-2*kp)-br)/(br**2-1d0)
            else
               b(j) = (c+a)/2d0 + (c-a)/bsum
     &               * (br**(2*kp-nsub-1)-br)/(br**2-1d0)
            end if
         else if ( meshtype.eq.-1 .and. nsub/2*2.ne.nsub ) then
            bsum = 2d0 * (br**(nsub+1)-1d0) / (br**2-1d0) - 1d0
            blength(j) = (c-a)/bsum * br**iabs(2*kp-nsub-1)
            if ( kp .le. (nsub+1)/2 ) then
               b(j) = (c+a)/2d0 - (c-a)/bsum
     &               * ( (br**(nsub+3-2*kp)-1d0)/(br**2-1d0) -.5d0)
            else
               b(j) = (c+a)/2d0 + (c-a)/bsum
     &               * ( (br**(2*kp-nsub-1)-1d0)/(br**2-1d0) -.5d0)
            end if
         else if ( meshtype.eq.-2 .and. nsub/2*2.eq.nsub ) then
            bsum = br * (br**nsub-1d0) / (br**2-1d0)
            if ( kp .le. nsub/2 ) then
               blength(j) = -a/bsum * br**iabs(2*kp-nsub-1)
               b(j) = a/bsum * (br**(nsub+3-2*kp)-br)/(br**2-1d0)
            else
               blength(j) = c/bsum * br**iabs(2*kp-nsub-1)
               b(j) = c/bsum * (br**(2*kp-nsub-1)-br)/(br**2-1d0)
            end if
         else
            print*,'ERROR: Unknown Mesh Initialization' 
         end if
      end do

      end
