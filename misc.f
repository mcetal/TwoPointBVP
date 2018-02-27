**********************************************************************
      subroutine writefn(nsub,nnd,ieval,xnd,u,ugvn,up,sigma,leafid)
      integer nsub, nnd, ieval, nmxsub
      parameter (nmxsub=1024)
      integer leafid(nmxsub)
      real*8 xnd(*), u(*), ugvn(*), up(*), sigma(*)
**********************************************************************
      integer kp, j, l
         open(file='exout.plt',unit=23,status='unknown')

c        write (23,'(a,i4,a9,i3)') '# nsub =',nsub,'nnd =',nnd
c        write (23,'(a,4a15)') '# x','Sol','Real','Err','Sigma'
         do kp = 1, nsub
            j = leafid(kp)
            do l = 1 + nnd*(j-1), nnd*j
               if ( ieval/2 - ieval/4*2 .eq. 1 ) then
                  write (23,'(f10.5,f15.9,f15.9,e15.6,e15.7)')
     &              xnd(l),u(l), ugvn(l), up(l), sigma(l)
               else
                  write (23,'(f10.5,f15.9,f15.9,e15.6,e15.7)')
     &              xnd(l),u(l), ugvn(l), u(l)-ugvn(l), sigma(l)
               end if
            end do
         end do

         close(unit=23, status='keep')
      end

**********************************************************************
      subroutine writemesh(countup, located,leafid,b,c)
      integer nmxsub
      parameter (nmxsub=1024)
      integer countup, located, leafid(nmxsub), inited, i
      real*8 b(nmxsub), c
      data inited /0/
**********************************************************************
      if ( inited .eq. 0 ) then
         open(file='exout.tree',unit=25,status='unknown')
         inited = 1
      end if

      do i = 1, located
         write (25,'(f8.4,i3)')  b(leafid(i)), countup
      end do
      write (25,'(f8.4,i3)')  c, countup

      end

**********************************************************************
      subroutine writestatus(countup, l2err, errest, trun,
     &                divided, updated, located, discond, mergcond)
      integer countup, divided, updated, located, inited
      real*8 l2err, errest, trun, discond, mergcond
      data inited /0/
**********************************************************************
      if ( inited .eq. 0 ) then
         open(file='exout.conv',unit=27,status='unknown')
c        write(27,'(a2,2a4,a4,3a11,2a15)'),
c    &        '#C','Div','Del','Loc',
c    &        'P_i','Merge','Trun', 'Conv Err','E: Real,Conv'
         inited = 1
      end if

      write(27,'(i2,2i4,i4,3e11.4,2e15.8)')
     &        countup,divided,max0(updated-2*divided,0),located,
     &        discond, mergcond, trun, errest, l2err
      end

**********************************************************************
      subroutine prerr2(nsub,nnd,blength,u,ugvn,leafid)
      integer nsub,nnd, nmxnd, nmxsub
      parameter (nmxnd=24,nmxsub=1024)
      integer leafid(nmxsub)
      real*8 dsqrt, blength(nmxsub), u(*), ugvn(*)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      common /geofn/ theta, ns, c, chnd01
c -------------------------------------------------
c  print Rel_2_Error = |u-ugvn|_2 / |ugvn|_2
**********************************************************************
      integer kp, j, l, lp
      real*8 err2, l2norm
         err2 = 0d0
         l2norm = 0d0
         do kp = 1, nsub
            j = leafid(kp)
            do lp = 1, nnd
                l = lp + nnd*(j-1)
                err2 = err2 + (ugvn(l)-u(l))**2 * blength(j) * ns(lp)
                l2norm = l2norm + ugvn(l)**2 * blength(j)
            end do
         end do
         print'(a,i4,a5,i2,a5,i5,a,e15.8)',
     &       'm=',nsub,'p=',nnd,'n=',nsub*nnd,
     &       ',   Relative L2 Error =',dsqrt(err2/l2norm)
      end

**********************************************************************
      function err2(nsub,nnd,blength,u,ugvn,leafid)
      integer nmxnd, nmxsub, nsub, nnd
      parameter (nmxnd=24,nmxsub=1024)
      integer leafid(nmxsub)
      real*8 err2, dsqrt, blength(nmxsub), u(*), ugvn(*)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      common /geofn/ theta, ns, c, chnd01
c -------------------------------------------------
c  print Rel_2_Error = |u-ugvn|_2 / |ugvn|_2
**********************************************************************
      integer kp, j, l, lp
      real*8 err, l2norm
         err = 0d0
         l2norm = 0d0
         do kp = 1, nsub
            j = leafid(kp)
            do lp = 1, nnd
                l = lp + nnd*(j-1)
                err = err + (ugvn(l)-u(l))**2 * blength(j) * ns(lp)
                l2norm = l2norm + ugvn(l)**2 * blength(j)
            end do
         end do
         err2 = dsqrt(err/l2norm)
      end

************************************************************************
      subroutine two2one(nnd,combined,first,second)
      integer nmxnd, nnd
      real*8 combined(nnd), first(nnd), second(nnd)
      parameter (nmxnd=24)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      common /geofn/ theta, ns, c, chnd01
c  -------------------------------------------------------------------
c   when two functions are defined on adjacent intervals with same length,
c   evaluate functional values on chevshev nodes of merged interval 
************************************************************************
      real*8 work(nmxnd), ti, sum
      integer i,k

c ----------------------
c f(x<0) = fl(2x+1)
c ----------------------
      call chftrans(work, first, nnd)
      do i = 1, nnd/2
         ti = dacos(2d0*c(i)+1d0)
         sum = 0d0
         do k = 1, nnd
            sum = sum + dcos( (k-1) * ti ) * work(k)
         end do
         combined(i) = sum
      end do
          
c -------------------------------------------------------
c f(x_mid = boundary between two intevals) = (fl(1)+fr(-1)) / 2d0
c f(x=0) = (fl(1)+fr(-1)) / 2d0 when nnd = odd
c And chebychev transfromation for second interval
c -------------------------------------------------------
      if ( (nnd+1)/2*2 .eq. nnd+1 ) then
         sum = 0d0
         do k = 1, nnd
            sum = sum + work(k)
         end do
      end if

      call chftrans(work, second, nnd)

      if ( (nnd+1)/2*2 .eq. nnd+1 ) then
         do k = 1, nnd*3/4
            sum = sum + (-1)**(k-1) * work(k)
         end do
         combined((nnd+1)/2) = sum / 2d0
      end if

c ----------------------
c f(x>0) = fr(2x-1)
c ----------------------
      do i = (nnd+1)/2+1, nnd
         ti = dacos(2d0*c(i)-1d0)
         sum = 0d0
         do k = 1, nnd*3/4
            sum = sum + dcos( (k-1) * ti ) * work(k)
         end do
         combined(i) = sum
      end do

      end
         
**********************************************************************
      function updateerr(located,leafid,updated,updateid,updatedone,
     &        blength,nnd,eval,seval)
      integer nmxnd, nmxsub, nnd
      parameter (nmxnd=24,nmxsub=1024)
      integer located, leafid(nmxsub)
      integer updated, updateid(nmxsub), updatedone(nmxsub)
      real*8 updateerr, eval(*), seval(*)
      real*8 combined(nmxnd), blength(nmxsub)
c -------------------------------------------------
c  updateerr = |eval-seval|_2 / |(eval+seval)/2|_2
c  Unfortuately, if there were updates in mesh then
c  eval(l) and seval(l) represent different nodal points
c  So use two2one subroutine for updated intervals
**********************************************************************
      integer lp, up, j, jl, jr
      real*8 bj, num, den

      lp = 1
      up = 1
      num = 0d0
      den = 0d0

      do while (lp .le. located)
         j = leafid(lp)
         jl = (j-1)*nnd + 1
c wasn't updated
         if ( up.gt.updated .or.
     &        j.ne.updateid(up) ) then
            bj = blength(j)
            call l2sum(num, den, seval(jl), eval(jl), nnd, bj)
            lp = lp + 1
c divided
         else if ( leafid(lp).eq.updateid(up) .and.
     &             updatedone(up).le.0 ) then
            bj = blength(j) * 2d0
            jr = (leafid(lp+1)-1)*nnd +1
            call two2one(nnd,combined, eval(jl),eval(jr))
            call l2sum(num, den, seval(jl), combined, nnd, bj)
            up = up + 2
            lp = lp + 2
         else
c merged
            bj = blength(j)
            jr = (updatedone(up)-1)*nnd +1
            call two2one(nnd,combined, seval(jl), seval(jr))
            call l2sum(num, den, combined, eval(jl), nnd, bj)
            up = up + 1
            lp = lp + 1
         end if
      end do

      updateerr = dsqrt(num/den)
      end

**********************************************************************
      subroutine l2sum(num, den, s1, s2, nnd, blength)
      integer nnd, l, nmxnd
      parameter (nmxnd=24)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      common /geofn/ theta, ns, c, chnd01
      real*8 s1(nnd), s2(nnd), num, den, blength
c -------------------------------------
c  Called by updateerr
**********************************************************************
      do l = 1, nnd
         num = num + (s1(l)-s2(l))**2 * blength * ns(l)
         den = den + (s1(l)+s2(l))**2 / 4d0 * blength 
      end do
      end

**********************************************************************
      subroutine savenodeval(located, leafid, nnd, val, sval)
      integer located, nmxsub, nnd
      parameter (nmxsub=1024)
      integer leafid(nmxsub)
      real*8 val(*), sval(*)
**********************************************************************
      integer kp, j, l

      do kp = 1, located
         j = leafid(kp)
         do l = 1 + nnd*(j-1), nnd*j
            sval(l) = val(l)
         end do
      end do
      end

