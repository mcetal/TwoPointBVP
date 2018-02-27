**************************************************************************
      subroutine discret(nnd,j,blength,v_l,v_r,localcond,eta,phil,phir)
      integer nmxnd, nmxsub, nnd, j
      parameter (nmxnd=24,nmxsub=1024)
      real*8 blength(nmxsub), localcond(nmxsub)
      real*8 v_l(*), v_r(*), phil(*), phir(*), eta(*)
      real*8 spdef, spbx, spxb
      common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
c------------------------------------------------------------
c this subroutine implements the
c Step 1 of the adaptive ode solver by Lee and Greengard, 1995 SISC 
c------------------------------------------------------------
c for each subinterval at the finest level, it sets up a nnd*nnd
c matrix that is derived from discretizing the integral equation 
c solution. then linpack subroutines are used to solve three linear
c systems with the matrix.
c------------------------------------------------------------
c input  : nmxnd = number of maximum nodes per subinterval
c          nnd = number of nodes specified by user
c          v_l, v_r, blength
c in/out : eta = solution with the right hand side f
c          phil = solution with the right hand side psil
c          phir = solution with the right hand side phir
c output : localcond = inverse condition number of local integral operator
**************************************************************************
      real*8 pcc(nmxnd,nmxnd), z(nmxnd)
      integer i, ipvt(nmxnd), k, kp, ip

        do i = 1, nnd
           ip = i + nnd*(j-1)
           do k = 1, nnd
              kp = k + nnd*(j-1)
              pcc(k,i) = blength(j)/2d0 * 
     &         ( phil(kp)*spbx(k,i)*v_l(ip)
     &         + phir(kp)*spxb(k,i)*v_r(ip) )
           end do
           pcc(i,i) = 1d0 + pcc(i,i)
        end do

c ----------------------------------------------------------
c sets up the right hand side vectors for f, psil, psir
c    solves P_i(eta_i) = f_i
c    solves P_i(phil_i) = psil_i
c    solves P_i(phir_i) = phir_i
c and copy it to the right place
c ----------------------------------------------------------
c       call dgefa(pcc,nmxnd,nnd,ipvt,info)
        call dgeco(pcc,nmxnd,nnd,ipvt,localcond(j),z)
        call dgesl(pcc,nmxnd,nnd,ipvt,eta(1+nnd*(j-1)),0)
        call dgesl(pcc,nmxnd,nnd,ipvt,phil(1+nnd*(j-1)),0)
        call dgesl(pcc,nmxnd,nnd,ipvt,phir(1+nnd*(j-1)),0)
      end

****************************************************************************
      subroutine mkcoef(nnd,k,j,blength,v_l,v_r,eta,phil,phir,
     &            alpha11,alpha21,alpha12,alpha22,delta_l,delta_r)
      integer nmxnd, nmxsub, nnd, k, j
      parameter (nmxnd=24,nmxsub=1024)
      real*8 blength(nmxsub)
      real*8 v_r(*), v_l(*), phil(*),phir(*),eta(*)
      real*8 alpha11(2*nmxsub),alpha12(2*nmxsub)
      real*8 alpha21(2*nmxsub),alpha22(2*nmxsub)
      real*8 delta_l(2*nmxsub),delta_r(2*nmxsub)
      real*8 spdef, spbx, spxb
      common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
c------------------------------------------------------------
c this subroutine implements the
c Step 1 of the adaptive ode solver by Lee and Greengard, 1995 SISC 
c------------------------------------------------------------
c it uses chebychev quadrature method to evaluate inner products on
c each subinterval
c------------------------------------------------------------
c input  : nmxnd = number of maximum nodes per subinterval
c          nnd = number of nodes specified by user
c          nsub = number of subintervals specified by user
c          k = node number in hierarchy tree
c          j = node id to update
c          blength, v_l, v_r
c          eta = stores solution with f as the right hand side 
c          phil = stores solution with psil as the right hand side 
c          phir= stores solution with psir as the right hand side 
c output : alpha11 = inner product of V_l and phi_l
c          alpha12 = inner product of V_l and phi_r
c          alpha21 = inner product of V_r and phi_l
c          alpha22 = inner product of V_r and phi_r
c          delta_l = inner product of V_l and eta
c          delta_r = inner product of V_r and eta
****************************************************************************
      integer i, l
      real*8 sc

      sc = blength(j) / 2d0
      alpha11(k) = 0d0
      alpha12(k) = 0d0
      alpha21(k) = 0d0
      alpha22(k) = 0d0
      delta_l(k) = 0d0
      delta_r(k) = 0d0

      do i = 1, nnd
         l = i + nnd*(j-1)
         alpha11(k) = alpha11(k) + v_l(l)* phil(l)*spdef(i)*sc
         alpha12(k) = alpha12(k) + v_l(l)* phir(l)*spdef(i)*sc
         alpha21(k) = alpha21(k) + v_r(l)* phil(l)*spdef(i)*sc
         alpha22(k) = alpha22(k) + v_r(l)* phir(l)*spdef(i)*sc
         delta_l(k) = delta_l(k) + v_l(l)* eta(l) *spdef(i)*sc
         delta_r(k) = delta_r(k) + v_r(l)* eta(l) *spdef(i)*sc
      end do
      end
