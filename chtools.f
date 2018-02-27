***********************************************************************
      subroutine chsetup(nnd)
      integer nnd, nmxnd
      parameter (nmxnd=24)
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      common /geofn/ theta, ns, c, chnd01
      real*8 cftm, citm, spdef, spbx, spxb
      common /chmat/ cftm(nmxnd,nmxnd), citm(nmxnd,nmxnd)
      common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
c------------------------------------------------------------
c input  :
c    nnd = number of nodes specified by a user
c output :
c    theta = Chebychev nodes on [0,pi] were stored in decreasing order
c            so that corresponding nodes in [-1,1] are in increasing order
c    chnd01 = Chebyshev nodes between [0,1]
c    cftm : ch(1:n) = cftm(1:n,1:n) * fn(1:n)
c    citm : fn(1:n) = citm(1:n,1:n)' * fn(1:n)
c    SPBX sigma  = Cn^{-1} SPINT_BXn Cn = \INT_{-1}^{x} sigma(t) dt
c    SPXB sigma  = Cn^{-1} SPINT_XBn Cn = \INT_{x}^{1} sigma(t) dt
***********************************************************************
        integer i,k
        real*8 work(nmxnd), defint

        do i = 1, nnd
           theta(i) = (nnd-i+.5d0)/nnd * pi
           ns(i) = dsin(theta(i)) * dsqrt(2d0)
           c(i) = dcos(theta(i))
           chnd01(i) = ( dcos(theta(i)) + 1d0 ) / 2d0
        end do

c------------------------------------------------------------
c computes terms in cosine inverse transform
c and stores them into matrix citm in the transposed format
c to facilitate data access during matrix multiplications.
c------------------------------------------------------------
        do  k = 1,nnd
           do i = 1, nnd
              citm(i,k) = dcos((i-1)*theta(k))
           end do
        end do

c------------------------------------------------------------
c computes terms in cosine forward transform with scaling factors
c and stores them into matrix cftm.
c    cftm(i,k) = 2d0/dble(nnd) * dcos((i-1)*theta(k))
c------------------------------------------------------------
        do k = 1,nnd
           cftm(1,k) = 1d0/dble(nnd)
           do i = 2, nnd
              cftm(i,k) = 2d0/dble(nnd)*citm(i,k)
            end do
        end do

c------------------------------------------------------------
c Spectral integration matrix
c   ctfm(:,i) = Cn En where En = ( 0, 0, ... , 1=i-th, ... , 0 )
c   SPBX Sigma  = Cn^{-1} SPINT_BXn Cn = \INT_{-1}^{x} Sigma(t) dt
c   SPXB Sigma  = Cn^{-1} SPINT_XBn Cn = \INT_{x}^{1} Sigma(t) dt
c------------------------------------------------------------
        do i = 1, nnd
           call chindef_bx(work,cftm(1,i),nnd)
           call chbtrans(spbx(1,i),work,nnd)
           call chindef_xb(work,cftm(1,i),nnd)
           call chbtrans(spxb(1,i),work,nnd)
           spdef(i) = defint(cftm(1,i),nnd)
        end do
      end
       
***********************************************************************
      subroutine chnodes(nnd,j,b,blength,xnd)
      integer nmxnd, nmxsub, nnd, j
      parameter (nmxnd=24,nmxsub=1024)
      real*8 xnd(nnd), b(nmxsub),blength(nmxsub)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      common /geofn/ theta, ns, c, chnd01
c --------------------------------------------------------------------
c return the scaled chebychev nodes on a subinterval [b(j),b(j)+blength(j)]
c input : b(j), blength(j)
c common : chnd01(nnd)
c output : xnd(nnd)
***********************************************************************
      integer i, l
        do i = 1, nnd
           l = i + nnd*(j-1)
           xnd(l) = b(j) + blength(j) * chnd01(i)
        end do
      end

*************************************************************************
      subroutine chftrans(ch, fn, nnd)
      integer nmxnd, nnd, i, j
      parameter (nmxnd=24)
      real*8 ch(nnd), fn(nnd)
      real*8 cftm, citm
      common /chmat/ cftm(nmxnd,nmxnd), citm(nmxnd,nmxnd)
c------------------------------------------------------------------------
c  Given function values at the chebyshev nodes (physical space)
c  calculate the corresponding chebyshev coefficients (fourier space).
*************************************************************************
        do i = 1, nnd
           ch(i) = 0d0
           do j = 1, nnd
              ch(i) = ch(i) + cftm(i,j)*fn(j) 
           end do
        end do
      end

*********************************************************************
      subroutine chbtrans(fn, ch, nnd)
      integer nmxnd, nnd, i, j
      parameter (nmxnd=24)
      real*8 ch(nnd), fn(nnd)
      real*8 cftm, citm
      common /chmat/ cftm(nmxnd,nmxnd), citm(nmxnd,nmxnd)
c------------------------------------------------------------------------
c  from the chebyshev coefficients to function values at the chebyshev nodes
*************************************************************************
        do i = 1, nnd
           fn(i) = 0d0
           do j = 1, nnd
              fn(i) = fn(i) + citm(j,i)*ch(j)
           end do
        end do
      end

*********************************************************************
      function defint(coef,nnd)
      integer  nnd
      real*8 defint,coef(0:nnd-1)
*********************************************************************
      integer i
        defint =  2d0 * coef(0)
        do i = 2, nnd-1, 2
           defint = defint - 2d0 * coef(i) / dble(i*i-1)
        end do
      end
	
*********************************************************************
      subroutine chindef_bx(chintfl,coeffl,nnd)
      integer i,nnd
      real*8 coeffl(0:nnd-1), chintfl(0:nnd-1)
c  ------------------------------------------------------------------
c this routine calculates the chebyshev coefficients of
c the left (normal) indefinite integrals
*********************************************************************

        chintfl(1) = coeffl(0) - coeffl(2)/2d0
        do i = 2,nnd-2
           chintfl(i) = (coeffl(i-1)-coeffl(i+1)) / dble(2*i)
        end do
        chintfl(nnd-1) = coeffl(nnd-2) / dble(2*nnd-2)

        chintfl(0) = (-1)**(nnd+1) * coeffl(nnd-1) / dble(2*nnd)
        do i = nnd-1, 1, -1
           chintfl(0) = chintfl(0) - (-1)**i * chintfl(i)
        end do
      end

*********************************************************************
      subroutine chindef_xb(chintfr,coeffr,nnd)
      integer i,nnd
      real*8 coeffr(0:nnd-1),chintfr(0:nnd-1)
c  ------------------------------------------------------------------
c this routine calculates the chebyshev coefficients of
c the right (backward) indefinite integrals
*********************************************************************

        chintfr(1) = coeffr(2)/2d0 - coeffr(0)
        do i = 2, nnd-2
           chintfr(i) = (coeffr(i+1)-coeffr(i-1)) / dble(2*i)
        end do
        chintfr(nnd-1) = -coeffr(nnd-2) / dble(2*nnd-2)

        chintfr(0) = coeffr(nnd-1) / dble(2*nnd)
        do i = 1, nnd-1
           chintfr(0) = chintfr(0) -  chintfr(i)
        end do
      end

