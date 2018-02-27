************************************************************************
      subroutine solve_sigma(nnd, updated, updatelist, updateid,
     &                    located, leaflist, leafid, parentlist,
     &                    blength,gl,gr,phil,phir,eta,
     &                    mergcond, discond, sigma)
      integer nmxsub, nnd, updated, located
      parameter (nmxsub=1024)
      integer updatelist(nmxsub), updateid(nmxsub)
      integer leafid(nmxsub), leaflist(nmxsub), parentlist(nmxsub)
      real*8 blength(nmxsub), discond, mergcond
      real*8 gl(*), gr(*), phil(*), phir(*), eta(*), sigma(*)
c ----------------------------------------------------------------
c
c    Solving an integeral equaion with unknown SIGMA
c    To get solution U(*), use solve_eval
c
c input:
c    nnd : order of solver
c    updated, updatelist, updateid for local integral equation solver
c    located, leaflist, leafid : list of subintervals
c        may be linked with the variables for updated subintervals
c    parentlist for hierarchical solver of global integral equation
c    blength(j), j = leafid(kp), kp = 1..located
c    gl, gr, phil, phir, eta at xnd(l+nnd*(j-1)), l = 1..nnd
c
c output :
c    discond, mergcond : Pivoting Not implement
c    sigma : the solution of the problem at xnd(l+nnd*(j-1)), l = 1..nnd
c
c common variables :
c    Used in discrete.f & initialized in chsetup(nnd) in chtools.f
c        common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
c    Used in recur.f & initialized in make*, free*, gentree in tree.f
c        common /hierarchy/ tree, node, leaf
c ----------------------------------------------------------------
      real*8 alpha11(2*nmxsub),alpha12(2*nmxsub)
      real*8 alpha21(2*nmxsub),alpha22(2*nmxsub)
      real*8 delta_l(2*nmxsub),delta_r(2*nmxsub)
      real*8 lambda1(2*nmxsub),lambda2(2*nmxsub)
      common /inner/ alpha11, alpha12, alpha21, alpha22,
     &                delta_l, delta_r, lambda1, lambda2
************************************************************************

      integer kp, k, j
      real*8 localcond(nmxsub)

c ---------------------------------------------
c  Solve local integral equations for eta, phi
c ---------------------------------------------
         do kp = 1, updated
            k = updatelist(kp)
            j = updateid(kp)
            call discret(nnd,j,blength, gl,gr,localcond,eta,phil,phir)
            call mkcoef(nnd,k,j,blength,gl,gr,eta,phil,phir,
     &              alpha11,alpha21,alpha12,alpha22,delta_l,delta_r)
         end do

c --------------------------------
c  Condtion numer of local solver
c --------------------------------
         discond = 1d0
         do kp = 1, located
            discond = dmin1( discond, localcond(leafid(kp)) )
         end do

c -------------------------------
c  Solve global integral equation
c -------------------------------
         call upward_sweep(located,parentlist,mergcond,
     &           alpha11,alpha12,alpha21,alpha22,delta_l,delta_r)
         call down_sweep(located,parentlist,alpha11,alpha12,
     &           alpha21,alpha22, delta_l,delta_r, lambda1,lambda2)
         call evalsigma(nnd,sigma,located,leaflist,leafid,
     &           lambda1,lambda2,phil,phir,eta)
      end

************************************************************************
      subroutine solve_eval(ieval, nnd, located, leaflist, leafid,
     &                    blength,xeval,gl,gr, sigma, upcom, ucom)
      integer nmxsub
      parameter (nmxsub=1024)
      integer ieval, nnd, located, leafid(nmxsub), leaflist(nmxsub)
      real*8 blength(nmxsub), xeval(*)
      real*8 gl(*), gr(*), ucom(*), upcom(*), sigma(*)
c ----------------------------------------------------------------
c
c    Using Sigma and G0,G1 (through gl,gr),
c    Evaluate U(*) / U'(*) at specified point(s)
c 
c input:
c -> ieval : solution to evaluate (applying Green's function with sigma)
c        (xeval:bnodes=1,chnodes=0)(upcom:yes=1,no=0)(ucom:yes=1,no=0)_2
c        Examples) 3 = both u and u' at chebyshev nodes
c        Examples) 4 + XY_2 = u/u' at b nodes
c    located, leaflist, leafid : list of subintervals
c        may be linked with the variables for updated subintervals
c    blength(j), j = leafid(kp), kp = 1..located
c -> xeval : nodal points where solution will be evaluated
c        ieval=0 then won't be referenced
c        btest(ieval,3) = 0 then xnd(l+nnd*(j-1)) = chebyshev nodes
c        btest(ieval,3) = 1 then b(j), j = leafid(1..located) and c
c        cf. btest(i,j) = mod(i/2**(j-1),2)
c    gl, gr at xnd(l+nnd*(j-1)), l = 1..nnd
c    sigma : the solution of the problem at xnd(l+nnd*(j-1)), l = 1..nnd
c
c output :
c -> upcom, ucom : solution at xeval()
c        will be referenced only when requested by ieval
c
c common variables :
c    Used in int_density.f & initialized in chsetup(nnd) in chtools.f
c        common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
c    Used in int_density.f & initialized in set_examples(a,b) in examples.f
c        common /bc/ q0, wron, zl1, zr1, zl0, zr0, ui2, ui1, ui0
c ----------------------------------------------------------------
      real*8 alpha11(2*nmxsub),alpha12(2*nmxsub)
      real*8 alpha21(2*nmxsub),alpha22(2*nmxsub)
      real*8 delta_l(2*nmxsub),delta_r(2*nmxsub)
      real*8 lambda1(2*nmxsub),lambda2(2*nmxsub)
      common /inner/ alpha11, alpha12, alpha21, alpha22,
     &                delta_l, delta_r, lambda1, lambda2
************************************************************************
      real*8 deffl(nmxsub+1), deffr(nmxsub+1)

c -------------------------------
c   Get the solution of ode
c -------------------------------
         if ( ieval .ne. 0 ) then
            call defint_intv(located, leaflist,
     &           deffl,deffr,lambda1,lambda2,
     &           alpha11,alpha12,alpha21,alpha22,delta_l,delta_r)
         end if

         if ( ieval/4 - ieval/8*2 .eq. 0 ) then
            call int_cheby(ieval,nnd,located,leafid,
     &           blength,xeval,gl,gr,sigma,deffl,deffr,ucom,upcom)
         else if ( ieval/4 - ieval/8*2 .eq. 1 ) then
ccc            call int_bnodes(ieval,located,leafid,
ccc     &           blength,xeval,gl,gr,sigma,deffl,deffr,ucom,upcom)
c
c  EVALUATE AT END POINTS ONLY!!
c
            call int_bnodes1(ieval,located,leafid,
     &           blength,xeval,gl,gr,sigma,deffl,deffr,ucom,upcom)
         endif

	end

