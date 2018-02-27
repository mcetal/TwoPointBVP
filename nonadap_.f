************************************************************************
      program nonadap
      integer nmxnd, nmxsub, nmxnn, nsub, nnd
      parameter (nmxnd=24,nmxsub=1024)
      parameter (nmxnn=nmxnd*nmxsub)
      integer leafid(nmxsub), leaflist(nmxsub), parentlist(nmxsub)
      integer leafid(nmxsub), leaflist(nmxsub)
      real*8 a, c, blength(nmxsub), b(nmxsub)
      real*8 xnd(nmxnn), gl(nmxnn),gr(nmxnn)
      real*8 pt(nmxnn),qt(nmxnn),ft(nmxnn)
      real*8 phil(nmxnn),phir(nmxnn),eta(nmxnn)
      real*8 ucom(nmxnn),upcom(nmxnn),sigma(nmxnn)
      real*8 ugvn(nmxnn), discond, mergcond
c----------------------------------------------------------------------
c this program implements the "fast adaptive numerical method for
c two-point boundary value problems" by June-Yub Lee and Leslie Greengard
c   NON-ADAPTIVE VERSION....
c 
c Read solve.f/solve_{sigma|eval} first to understand the
c non-adaptive ode solver by L. Greengard and V. Rokhlin in CPAM, 1991
c You will find the definitions of the above variables in solve.f
************************************************************************
      integer ieval, j, kp, verbose 

c--------------------------------------------------
c parameters set-up ; read HowToRun for details
c--------------------------------------------------
      a = 0d0
      c = 5d0
      call bcond_lin(a,c, 1d0,0d0,0d0, 1d0,0d0,1d0)
      ! boundary conditions (a,c,za1,za0,zc1,zc0,ga,gc)

      nnd = 16                 ! 16th order method
      call chsetup(nnd)
      ieval = 3                ! evaluate both u and u' at Chebyshev nodes

c ------------------------------------
c  generate equal length initial mesh
c ------------------------------------
      nsub = 20                ! Start with ten intervals
      call gentree(nsub, leafid, leaflist, parentlist)

      do kp = 1, nsub
         j = leafid(kp)
         b(j) = ( dble(kp-1)*c + dble(nsub-kp+1)* a ) / dble(nsub)
         blength(j) = (c-a) / dble(nsub)
      end do

c -------------------------------
c  ODE SOLVER
c -------------------------------
         do kp = 1, nsub
            j = leafid(kp)
            call chnodes(nnd,j,b,blength,xnd)
            call tabfns(j, nnd, xnd, pt, qt, ft, ugvn)
            call int_eqn(j,nnd,a,c,xnd,pt,qt,ft,gl,gr,phil,phir,eta)
         end do
         call solve_sigma(nnd, nsub, leaflist, leafid,
     &                    nsub, leaflist, leafid, parentlist,
     &                    blength,gl,gr,phil,phir,eta,
     &                    mergcond, discond, sigma)

         call solve_eval(ieval, nnd, nsub, leaflist, leafid,
     &                    blength,xnd,gl,gr, sigma, upcom, ucom)

c ------------------------------------
c  post-processing (Optional)
c ------------------------------------

      verbose = 5              ! If Ugvn is not defined
      verbose = 7              ! If Ugvn is defined
      if ( mod(verbose/2,2) .eq. 1 ) then
         call prerr2(nsub,nnd,blength,ucom,ugvn,leafid)
      end if

      call writefn(nsub,nnd,ieval,xnd,ucom,ugvn,upcom,sigma,leafid)

      stop
      end

**************************************************************************
      subroutine tabfns(j, nnd,xnd,pt,qt,ft,ugvn)
      integer j, nnd, l
      real*8 xnd(*), pt(*), qt(*), ft(*), ugvn(*)
      real*8 x, pi, eps1
      parameter (pi=3.141592653589793238d0)
c  --------------------------------------------------------------
c    input : j, nnd, xnd
c    output : pt(x), qt(x), ft(x), (optional) ugvn(x)
**************************************************************************
      eps1 = 1d0
      do l = 1 + nnd*(j-1), nnd*j
         x = xnd(l)
         pt(l) = 1d0 / x
         qt(l) = 1d0 - eps1**2 / x**2
         ft(l) = 0d0
         ugvn(l) = 1d0     ! J_eps1(x)
      end do
      end

