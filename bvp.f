c---------------
      subroutine BVP (nsub,nnd,a,za1,za0,ga,c,zc1,zc0,gc,xnd,pt,qt,ft,
     *                leafid,ucom,upcom,uprimea,uprimec,w,iw)
c---------------
c  Solves the BVP using June's code for the given arrays PT,QT,FT
c  at the points X = XND:  
c      u'' + pt(x) u'(x) + q(x) u(x) = f(x)
c  The coefficients ZA1, ZA0 of the boundary condition at X = A 
c  are given by 
c     za1 u'(a) + za0 u(a) = GA.
c  The coefficients ZC1 and ZC0 of the boundary condition at X = C are
c  given by 
c     zc1 u'(c) + zc0 u(c) = GC.
c
c  The solution U and U' are returned at the chebychev nodes in 
c  UCOM and UPCOM, as is the value of U' at x = A in UPRIMEA and
c  U' at x = C in UPRIMEC.
c
      implicit real*8 (a-h,o-z)
      parameter (nmxsub=1024)
      dimension xnd(nsub*nnd),pt(nsub*nnd),qt(nsub*nnd),ft(nsub*nnd),
     *          ucom(nsub*nnd),upcom(nsub*nnd)
      dimension w(*)
      integer leafid(nmxsub),iw(*)
c
c  carve up workspace
c
         npnts = nsub*nnd
         ib = 1
         iblength = ib + nmxsub
         ibupcom = iblength + nmxsub
         igl = ibupcom + npnts
         igr = igl + npnts
         iphil = igr + npnts
         iphir = iphil + npnts
         ieta = iphir + npnts
         isigma = ieta + npnts
         ileaflist = 1
         iparentlist = ileaflist + nmxsub
c
ccc         a = 0d0
ccc         c = 1d0
ccc         ga = 0.d0
ccc         gc = dexp(1.d0) +  0.25d0
ccc         ga = 0.d0
ccc         gc = 0.d0
         call bcond_lin(a, c, za1, za0, zc1, zc0, ga, gc)
            ! boundary conditions (a,c,za1,za0,zc1,zc0,ga,gc)
c
c -------------------------------
c  ODE SOLVER
c -------------------------------
         do kp = 1, nsub
            j = leafid(kp)
            call int_eqn(j,nnd,a,c,xnd,pt,qt,ft,w(igl),w(igr),w(iphil),
     *                   w(iphir),w(ieta))
         end do
         call solve_sigma(nnd, nsub, iw(ileaflist), leafid, nsub, 
     &                    iw(ileaflist), leafid, iw(iparentlist),
     &                    w(iblength),w(igl),w(igr),w(iphil),w(iphir),
     &                    w(ieta), mergcond, discond, w(isigma))
c
c   get u, u' at chebychev nodes
c
         ieval = 3
         call solve_eval(ieval, nnd, nsub, iw(ileaflist), leafid,
     &                   w(iblength),xnd,w(igl),w(igr), w(isigma),
     &                   upcom, ucom)
c
c   get u' at end point (as well as all other b nodes)
c
         ieval = 6
         call solve_eval(ieval, nnd, nsub, iw(ileaflist), leafid,
     &                    w(iblength),w(ib),w(igl),w(igr), w(isigma),
     &                    w(ibupcom), bucom)
         uprimea = w(ibupcom)
         uprimec = w(ibupcom+nsub)
 1000    format (f8.3,1x,e12.6,1x,e12.6,5x,e12.6,1x,e12.6)
c
ccc         call CHECK (nsub,nnd,leafid,xnd,ucom,upcom,upc)
c
      return
      end
