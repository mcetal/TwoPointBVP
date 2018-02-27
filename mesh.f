**********************************************************************
      function monfn(nnd, j, sigma)
      integer nnd, j, nmxnd
      parameter (nmxnd=24)
      real*8 monfn, sigma(*), ch(nmxnd)
c------------------------------------------------------------------------
c   A monitor function (phi) is defined as
c   phi(x) ^ (nnd-2)
c   = (d/dx)^(nnd-1) \int^x sigma(t) dt
c   = (d/dx)&(nnd-1) { ch_(nnd-1) T_(nnd-1) + ... }
c             where ch_(nnd-1) = sigma_(nnd-2) / (2*nnd-2)
c   = ch_(nnd-1) * (nnd-1)! * 2^(nnd-1)
c   \propto sigma_(nnd-2)
c ------------------------------------------------------------
c  TO DO :
c  1. find better monitor function
**********************************************************************
        call chftrans(ch, sigma(1+nnd*(j-1)), nnd)
        monfn = dabs(ch(nnd-1)) + dabs(ch(nnd)-ch(nnd-2))
*       print*,j, dabs(ch(nnd-1)),dabs(ch(nnd)-ch(nnd-2))
      end


**********************************************************************
      subroutine doublemesh(b,blength,
     &       divided, updated,updateid,updatelist,updatedone,
     &       located, leafid,leaflist)
      integer nmxsub
      parameter (nmxsub=1024)
      real*8 b(nmxsub),blength(nmxsub)
      integer divided,updated,located, leafid(nmxsub),leaflist(nmxsub)
      integer updateid(nmxsub), updatelist(nmxsub), updatedone(nmxsub)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
c------------------------------------------------------------------------
c Divide all intervals halves
c ---------------------------
c   in/out : b, blength
c   input :  located, leafid, leaflist
c   output : divided, updated,updateid,updatelist,updatedone
c   common /hierarchy/ updated though make/free children
**********************************************************************
      integer j, k, kp

c ----------------------
c  Divide
c ----------------------
      updated = 0
      divided = 0
      if ( 2*located .gt. nmxsub ) then
         print*, 'EXHAUST all subintervals, can''t double mesh'
         return
      end if
      do kp = 1, located
         k = leaflist(kp)
         j = leafid(kp)
             call makechildren(k)
             updateid(updated+1) = j
             updateid(updated+2) = tree(3,tree(3,k))
             updatelist(updated+1) = tree(2,k)
             updatelist(updated+2) = tree(3,k)
             updatedone(updated+1) = 0
             updatedone(updated+2) = 0
             blength(tree(3,tree(3,k))) = blength(j) / 2d0
             b(tree(3,tree(3,k))) = b(j) + blength(j) / 2d0
             blength(j) = blength(j) / 2d0
             updated = updated + 2
             divided = divided + 1
      end do
             
      end

**********************************************************************
      subroutine topmesh(tolratio, sumtrun, nnd, b,blength,sigma,
     &       divided, updated,updateid,updatelist,updatedone,
     &       located, leafid,leaflist)
      real*8 sigma(*)
      integer nmxsub, nnd
      parameter (nmxsub=1024)
      real*8 tolratio,maxtol, mintol, maxtrun, sumtrun
      real*8 b(nmxsub), blength(nmxsub), monitor(nmxsub), monfn
      integer located, leafid(nmxsub), leaflist(nmxsub)
      integer divided, updated
      integer updateid(nmxsub), updatelist(nmxsub), updatedone(nmxsub)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
c------------------------------------------------------------------------
c Divide if Mi >= max(Mmax*tolratio)
c Merge if Mi+Mj <= min(M_div) / 2**nnd
c ---------------------------
c   input : nnd, tolratio, sigma
c   output : sumtrun
c   in/out : b, blength
c   input : nnd, located, leafid, leaflist
c   output : divided, updated,updateid,updatelist,updatedone
c   common /hierarchy/ updated though make/free children
c ---------------------------
c  TO DO :
c  1. Can you make an automatic mesh refinement routine without
c     user specified constant TOLRATIO ?
**********************************************************************
      integer j, k, kp, prnted

      prnted = 0

c ----------------------
c  Scan Monitor function
c ----------------------
      maxtrun = 0d0
      sumtrun = 0d0
      do kp = 1, located
         j = leafid(kp)
         monitor(j) = monfn(nnd, j, sigma)
         maxtrun = dmax1(maxtrun, monitor(j))
         sumtrun = sumtrun + blength(j)**(nnd-1)*monitor(j)
*        print'(f12.5,g10.3,i6,f12.7)',b(j),blength(j),j,monitor(j)
      end do

      updated = 0
      divided = 0
      mintol = 1d38
      maxtol = maxtrun * tolratio

c ----------------------
c  Divide
c ----------------------
      do kp = 1, located
         k = leaflist(kp)
         j = leafid(kp)
         if ( monitor(j).ge.maxtol ) then
            mintol = dmin1(mintol, monitor(j)/2d0**nnd)
            if ( prnted.eq.0 .and. leaf(0).gt.nmxsub ) then
               print*,'EXHAUST all subintervals, can''t divide it'
               prnted = 1
            else
               call makechildren(k)
               updateid(updated+1) = j
               updateid(updated+2) = tree(3,tree(3,k))
               updatelist(updated+1) = tree(2,k)
               updatelist(updated+2) = tree(3,k)
               updatedone(updated+1) = 0
               updatedone(updated+2) = 0
               blength(tree(3,tree(3,k))) = blength(j) / 2d0
               b(tree(3,tree(3,k))) = b(j) + blength(j) / 2d0
               blength(j) = blength(j) / 2d0
               updated = updated + 2
               divided = divided + 1
            end if
         end if
      end do

c ----------------------
c  Merge
c ----------------------
      do kp = 1, located
         k = leaflist(kp)
         j = leafid(kp)
         if ( tree(3,tree(1,k)).eq.k .and.
     &      tree(2,tree(1,k)).eq.leaflist(kp-1) .and.
     &      monitor(j)+monitor(leafid(kp-1)) .lt. mintol ) then
            call freechildren(tree(1,k))
            updateid(updated+1) = leafid(kp-1)
            updatelist(updated+1) = tree(1,k)
            updatedone(updated+1) = j
            blength(leafid(kp-1)) = blength(leafid(kp-1)) + blength(j)
            updated = updated + 1
         end if
      end do
             
      end

