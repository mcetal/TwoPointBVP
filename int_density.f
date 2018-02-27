*************************************************************************
      subroutine evalsigma(nnd,sigma,located,leaflist,leafid,
     &           lambda1,lambda2,phil,phir,eta)
      integer nmxsub, nnd
      parameter (nmxsub=1024)
      integer located, leaflist(nmxsub), leafid(nmxsub)
      real*8 lambda1(2*nmxsub), lambda2(2*nmxsub)
      real*8 sigma(*), phil(*), phir(*), eta(*)
c -------------------------------------------------
c   Step 2.C of the algorithm by Lee and Greengard, SISC 1995
c -------------------------------------------------
c     sigma = (def) eta + lambda1*phi_l + lambda2*phi_r
*************************************************************************
      integer i, kp, k, j

      do kp = 1, located
           k = leaflist(kp)
           j = leafid(kp)
           do i = (j-1)*nnd+1, j*nnd
              sigma(i) = eta(i)+lambda1(k)*phil(i)+lambda2(k)*phir(i)
           end do
      end do
      end

*************************************************************************
      subroutine defint_intv(nsub,leaflist, deffl,deffr,
     &  lambda1,lambda2,alpha11,alpha12,alpha21,alpha22,delta_l,delta_r)
      integer nmxsub, nsub
      parameter (nmxsub=1024)
      integer leaflist(nmxsub)
      real*8 deffl(nmxsub+1), deffr(nmxsub+1)
      real*8 alpha11(2*nmxsub),alpha12(2*nmxsub)
      real*8 alpha21(2*nmxsub),alpha22(2*nmxsub)
      real*8 delta_l(2*nmxsub),delta_r(2*nmxsub)
      real*8 lambda1(2*nmxsub),lambda2(2*nmxsub)
c -------------------------------------------------
c   Step 3.A of the algorithm by Lee and Greengard, SISC 1995
c -------------------------------------------------
c  deffl(j) = \int_A^B(j-1) { vl sigma }
c  deffr(j) = \int_B(j)^C { vr sigma }
c      where sigma = eta + lambda1*phi_l + lambda2*phi_r
*************************************************************************
      integer k, kp

        deffl(1) = 0d0
        do kp = 1, nsub
           k = leaflist(kp)
           deffl(kp+1) = deffl(kp) +  delta_l(k)
     &            + lambda1(k)*alpha11(k) + lambda2(k)*alpha12(k)
        end do

        deffr(nsub+1) = 0d0
        do kp = nsub, 1, -1
           k = leaflist(kp)
           deffr(kp) = deffr(kp+1) + delta_r(k)
     &            + lambda1(k)*alpha21(k) + lambda2(k)*alpha22(k)
        end do
      end

*************************************************************************
      subroutine int_cheby(ieval,nnd,located,leafid,
     &           blength,xeval,gl,gr,sigma,deffl,deffr,ucom,upcom)
      integer ieval, nmxnd, nmxsub, nnd
      parameter (nmxnd=24,nmxsub=1024)
      integer located, leafid(nmxsub)
      real*8 blength(nmxsub), xeval(*)
      real*8 glp_x, grp_x
      real*8 gl(*), gr(*), ucom(*), upcom(*), sigma(*)
      real*8 deffl(nmxsub+1), deffr(nmxsub+1)
      real*8 spdef, spbx, spxb
      common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
      real*8 q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
c------------------------------------------------------------------------
c   Step 3.B of the algorithm by Lee and Greengard, SISC 1995
c------------------------------------------------------------------------
c    Given Sigma (a density funtion) and G0 (through ul,ur,vl,vr)
c        Evaluate phi on all of chebyshev nodal points (nnd*located)
c    Phi(x) = ul(x) \int_a^x{vl sigma} + ur(x) \int_x^b{vr sigma}
c    Sol(x) = Ui(x) + Uh(x) = ui2 * x**2 + ui1 * x + ui0 + Phi(x)
c
c  G(x,t) = gr(x) gl(t) / wron  or  gr(t) gl(x) / wron,
c      where wron(gl,gr) = gl(x) grp(x) - glp(x) gr(x) = constant
c------------------------------------------------------------------------
c   input variables: nmxsub, nnd, located, leafid
c                   sigma, gl, gr, xnd, blength, deffl, deffr
c   (optional) output variables : ucom, upcom
*************************************************************************
      integer i, j, k, kp, ip, lp
      real*8 scale, sum_r, sum_l

      do kp = 1, located
           j = leafid(kp)

c ----------------------------------------------------------------
c  ucom(x) = \int G0(x,t) sigma(t) dt
c         = ul(x) * int_A^x {vl*sigma} + ur(x) int_x^C {vr*sigma}
c  upcom(x) = \int G1(x,t) sigma(t) dt
c         = ulp(x) * int_A^x {vl*sigma} + urp(x) int_x^C {vr*sigma}
c ----------------------------------------------------------------
          scale = 0.5d0 * blength(j)
          do i = 1, nnd
             ip = i + nnd * (j-1)
             sum_l = 0d0
             sum_r = 0d0
             do k = 1, nnd
                lp = k + nnd * (j-1)
                sum_l = sum_l + spbx(i,k)*gl(lp)* sigma(lp)
                sum_r = sum_r + spxb(i,k)*gr(lp)* sigma(lp)
             end do
             if ( ieval/2 - ieval/4*2 .eq. 1 ) then
                upcom(ip) = ( grp_x(xeval(ip))*(deffl(kp)+scale*sum_l)
     &                      + glp_x(xeval(ip))*(scale*sum_r+deffr(kp+1))
     &                      ) / wron  + 2d0*ui2*xeval(ip) + ui1
             end if
             if ( ieval - ieval/2*2 .eq. 1 ) then
                ucom(ip) = ( gr(ip)*(deffl(kp)+scale*sum_l)
     &                     + gl(ip)*(scale*sum_r+deffr(kp+1)) ) / wron
     &                   + ui2*xeval(ip)**2 + ui1*xeval(ip) + ui0
             end if
	  end do
      end do
      end

*************************************************************************
      subroutine int_bnodes(ieval,located,leafid,
     &           blength,xeval,gl,gr,sigma,deffl,deffr,ucom,upcom)
      integer ieval, nmxsub
      parameter (nmxsub=1024)
      integer located, leafid(nmxsub)
      real*8 blength(nmxsub), xeval(*)
      real*8 glp_x, grp_x, gl_x, gr_x
      real*8 gl(*), gr(*), ucom(*), upcom(*), sigma(*)
      real*8 deffl(nmxsub+1), deffr(nmxsub+1)
      real*8 q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
c------------------------------------------------------------------------
c   Step 3.B of the algorithm by Lee and Greengard, SISC 1995
c------------------------------------------------------------------------
c    Given Sigma (a density funtion) and G0 (through ul,ur,vl,vr)
c        Evaluate phi on each nodal point b(1)=a, ... , b(i), ... , c
c    Phi(x) = ul(x) \int_a^x{vl sigma} + ur(x) \int_x^b{vr sigma}
c    Sol(x) = Ui(x) + Uh(x) = ui2 * x**2 + ui1 * x + ui0 + Phi(x)
c
c  G(x,t) = gr(x) gl(t) / wron  or  gr(t) gl(x) / wron,
c      where wron(gl,gr) = gl(x) grp(x) - glp(x) gr(x) = constant
c------------------------------------------------------------------------
c   input variables: nmxsub, nnd, located, leafid
c                   sigma, gl, gr, xnd, blength, deffl, deffr
c   (optional) output variables : ucom, upcom
*************************************************************************
      integer j, kp
      real*8 x

      do kp = 1, located+1
         if ( kp .ne. located+1) then
            j = leafid(kp)
            x = xeval(j)
         else
            j = leafid(located)
            x = xeval(j) + blength(j)
         end if

c ----------------------------------------------------------------
c  ucom(x) = \int G0(x,t) sigma(t) dt
c         = ul(x) * int_A^x {vl*sigma} + ur(x) int_x^C {vr*sigma}
c  upcom(x) = \int G1(x,t) sigma(t) dt
c         = ulp(x) * int_A^x {vl*sigma} + urp(x) int_x^C {vr*sigma}
c ----------------------------------------------------------------
         if ( ieval/2 - ieval/4*2 .eq. 1 ) then
            upcom(kp) = ( grp_x(x)*deffl(kp) + glp_x(x)*deffr(kp) )
     &                  / wron + 2d0*ui2*x + ui1
         end if
         if ( ieval - ieval/2*2 .eq. 1 ) then
            ucom(kp) = ( gr_x(x)*deffl(kp) + gl_x(x)*deffr(kp) )
     &                 / wron + ui2*x**2 + ui1*x + ui0
         end if
      end do
      end
*************************************************************************
      subroutine int_bnodes1(ieval,located,leafid,
     &           blength,xeval,gl,gr,sigma,deffl,deffr,ucom,upcom)
c
c  DOES END POINTS ONLY!
c
      integer ieval, nmxsub
      parameter (nmxsub=1024)
      integer located, leafid(nmxsub)
      real*8 blength(nmxsub), xeval(*)
      real*8 glp_x, grp_x, gl_x, gr_x
      real*8 gl(*), gr(*), ucom(*), upcom(*), sigma(*)
      real*8 deffl(nmxsub+1), deffr(nmxsub+1)
      real*8 q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
c------------------------------------------------------------------------
c   Step 3.B of the algorithm by Lee and Greengard, SISC 1995
c------------------------------------------------------------------------
c    Given Sigma (a density funtion) and G0 (through ul,ur,vl,vr)
c        Evaluate phi on each nodal point b(1)=a, ... , b(i), ... , c
c    Phi(x) = ul(x) \int_a^x{vl sigma} + ur(x) \int_x^b{vr sigma}
c    Sol(x) = Ui(x) + Uh(x) = ui2 * x**2 + ui1 * x + ui0 + Phi(x)
c
c  G(x,t) = gr(x) gl(t) / wron  or  gr(t) gl(x) / wron,
c      where wron(gl,gr) = gl(x) grp(x) - glp(x) gr(x) = constant
c------------------------------------------------------------------------
c   input variables: nmxsub, nnd, located, leafid
c                   sigma, gl, gr, xnd, blength, deffl, deffr
c   (optional) output variables : ucom, upcom
*************************************************************************
      integer j, kp
      real*8 x

      do kp = 1, located+1,located
         if ( kp .ne. located+1) then
            j = leafid(kp)
            x = xeval(j)
         else
ccc            kp = located+1
            j = leafid(located)
            x = xeval(j) + blength(j)
         end if

c ----------------------------------------------------------------
c  ucom(x) = \int G0(x,t) sigma(t) dt
c         = ul(x) * int_A^x {vl*sigma} + ur(x) int_x^C {vr*sigma}
c  upcom(x) = \int G1(x,t) sigma(t) dt
c         = ulp(x) * int_A^x {vl*sigma} + urp(x) int_x^C {vr*sigma}
c ----------------------------------------------------------------
         if ( ieval/2 - ieval/4*2 .eq. 1 ) then
            upcom(kp) = ( grp_x(x)*deffl(kp) + glp_x(x)*deffr(kp) )
     &                  / wron + 2d0*ui2*x + ui1
         end if
         if ( ieval - ieval/2*2 .eq. 1 ) then
            ucom(kp) = ( gr_x(x)*deffl(kp) + gl_x(x)*deffr(kp) )
     &                 / wron + ui2*x**2 + ui1*x + ui0
         end if
      end do
      end

