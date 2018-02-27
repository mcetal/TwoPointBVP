*****************************************************************
      subroutine upward_sweep(nsub,parentlist,mergcond,
     1        alpha11,alpha12,alpha21,alpha22,delta_l,delta_r)
      integer nmxsub, nsub
      parameter (nmxsub=1024)
      integer parentlist(nmxsub)
      real*8 mergcond
      real*8 alpha11(2*nmxsub),alpha12(2*nmxsub)
      real*8 alpha21(2*nmxsub),alpha22(2*nmxsub)
      real*8 delta_l(2*nmxsub),delta_r(2*nmxsub)
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
      common /hierarchy/ tree, node, leaf
c  -------------------------------------------------------------
c     Merge Condition Number = Max (1-A21b*A12a) in absolute value 
c  -------------------------------------------------------------
c    Step 2.A in the paper by Lee and Greengard, SISC 1995
c       Given alpha, delta in the finest level (nlevel)
c       recusively calculate alpha and delta in all levels
*****************************************************************
      integer k, pc, pa, pb
      real*8 deter

      mergcond = 0d0
c    ------------------------------------------------
c              pc = parentlist(k)
c              alpha_Coarse = alpha(pc)
c              alpha_finer_A = alpha(tree(2),pc)
c              alpha_finer_B = alpha(tree(3),pc)
c              See Eq. 67 & 61
c    ------------------------------------------------

      do k = nsub-1, 1, -1
         pc = parentlist(k)
         pa = tree(2,pc)
         pb = tree(3,pc)

         deter = 1d0 - alpha21(pb) * alpha12(pa)
         mergcond = dmax1(mergcond, dabs(deter))

         alpha11(pc) = alpha11(pb) + (1d0-alpha11(pb)) *
     &         (alpha11(pa)-alpha12(pa)*alpha21(pb)) / deter
         alpha21(pc) = alpha21(pa) +
     &         alpha21(pb)*(1d0-alpha22(pa))*(1d0-alpha11(pa))/deter
         alpha12(pc) = alpha12(pb) +
     &         alpha12(pa)*(1d0-alpha22(pb))*(1d0-alpha11(pb))/deter
         alpha22(pc) = alpha22(pa) + (1d0-alpha22(pa)) *
     &         (alpha22(pb)-alpha12(pa)*alpha21(pb)) / deter

         delta_l(pc) = (1d0-alpha11(pb)) / deter * delta_l(pa) +
     &         delta_l(pb) +
     &         (alpha11(pb)-1d0)*alpha12(pa)/deter * delta_r(pb)
         delta_r(pc) = (1d0-alpha22(pa)) / deter * delta_r(pb) +
     &         delta_r(pa) +
     &         (alpha22(pa)-1d0)*alpha21(pb)/deter * delta_l(pa)

      end do
      end

*****************************************************************
	subroutine down_sweep(nsub,parentlist,alpha11,alpha12,
     1     alpha21,alpha22, delta_l,delta_r, lambda1,lambda2)
      integer nmxsub, nsub
      parameter (nmxsub=1024)
      real*8 alpha11(2*nmxsub),alpha12(2*nmxsub)
      real*8 alpha21(2*nmxsub),alpha22(2*nmxsub)
      real*8 delta_l(2*nmxsub),delta_r(2*nmxsub)
      real*8 lambda1(2*nmxsub),lambda2(2*nmxsub)
      integer parentlist(nmxsub)
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
      common /hierarchy/ tree, node, leaf
c  -------------------------------------------------------------
c    Step 2.B in the paper by Lee and Greengard, SISC 1995
c       Given alpha, delta in all levels
c       recusively calculate lambda in all levels
*****************************************************************
	integer k, pc, pa, pb
	real*8 deter

c    ---------------------
c      initial data
c      lambda^{0th-level}_{1,2} = 0
c    ---------------------
        lambda1(1+1) = 0d0
        lambda2(1+1) = 0d0

c    ------------------------------------------------
c              pc = parentlist(k)
c              lambda_Coarse = lambda(pc)
c              lambda_finer_A = lambda(tree(2),pc)
c              lambda_finer_B = lambda(tree(3),pc)
c              See Eq. 73 & 75
c    ------------------------------------------------
        do k = 1, nsub-1
           pc = parentlist(k)
           pa = tree(2,pc)
           pb = tree(3,pc)
           deter = 1d0 - alpha21(pb) * alpha12(pa)

           lambda1(pa) = lambda1(pc)
           lambda2(pa) = -(1d0-alpha11(pa))*alpha21(pb)/deter *
     &         lambda1(pc) + (1d0-alpha22(pb))/deter*lambda2(pc) +
     &         alpha21(pb)/deter*(delta_l(pa)-alpha12(pa)*delta_r(pb))
     &             - delta_r(pb)
           lambda1(pb) = +(1d0-alpha11(pa))/deter*lambda1(pc) -
     &         (1d0-alpha22(pb))*alpha12(pa)/deter*lambda2(pc) +
     &         alpha12(pa)/deter*(delta_r(pb)-alpha21(pb)*delta_l(pa))
     &             - delta_l(pa)
           lambda2(pb) = lambda2(pc)
        end do
      end

