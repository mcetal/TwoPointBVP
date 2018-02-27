**********************************************************************
      block data nodetree
      integer nmxsub
      parameter (nmxsub=1024)
c ---------------------------------------------------------
c  node(1) : next look up pointer, number of used nodes+2 (initial=2)
c  node(2) : the whole interval
c  ordinary node : 2 <= j <= 2*nmxsub
c      node(j) = 0   : never used, allocate j-th node
c      non-zero(ge2) : next available node number
c ---------------------------------------------------------
c  grand-parent: node number = 2
c     tree(1,1) = 0, tree(2,1) = 0, tree(3,1) = 0
c     tree(1,2) = no parent(1), tree(2,2) = 3, tree(3,2) = 4
c  a node with children: node number = j
c     tree(1,j) = parent of j-th node
c     tree(2,j) = left child of j-th node
c     tree(3,j) = right child of j-th node
c  a leaf, node without child: node number = j
c     tree(1,j) = parent of j-th node
c     tree(2,j) = 0
c     tree(3,j) = leaf number
c ---------------------------------------------------------
c  leaf(0) : next look up pointer, number of used leaves
c  leaf( leaf(0) ) : current avaible leaf if not zero
**********************************************************************
      integer tree(3,2*nmxsub),node(2*nmxsub), leaf(0:nmxsub)
      common /hierarchy/ tree, node, leaf
      data node(1),node(2),node(3) /3,0,0/
      data leaf(0),leaf(1),leaf(2) /2,0,0/
      data tree(1,2),tree(2,2),tree(3,2) /1,0,1/
      end

**********************************************************************
      subroutine gentree(n, updateid, updatelist, parentlist)
      integer nmxsub, n, j
      parameter (nmxsub=1024)
      integer updateid(nmxsub), updatelist(nmxsub), parentlist(nmxsub)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
c ----------------------------------------------
c   Example:n=1    Example:n=4     Example:n=7
c   2              2 - 3 - 5       2 - 3 - 5 - 9
c                    |     6         |   |     10
c                    - 4 - 7         |   - 6 - 11
c                          8         |         12
c                                    - 4 - 7 - 13
c                                        |     14
c                                        - 8
**********************************************************************
        do j = 2, n
           tree(1,j) = (j+1) / 2
           tree(2,j) = 2*j-1
           tree(3,j) = 2*j
        end do

        do j = n+1, 2*n
           tree(1,j) = (j+1) / 2 
           tree(2,j) = 0
           tree(3,j) = j - n
        end do

        node(1) = 2*n + 1
        leaf(0) = n + 1
        call listtree(2, n, updateid, updatelist, parentlist)
      end

**********************************************************************
      subroutine makechildren(j)
      integer nmxsub, j
      parameter (nmxsub=1024)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub),node(2*nmxsub), leaf(0:nmxsub)
c ----------------------------------------------
c  A LEAF becomes a parent from geting two children
c  one leaf comes from leaf stack and the other from parent
**********************************************************************
      integer leftchild, rightchild, newleaf

        leftchild = node(node(1))
           if ( leftchild .le. 0 ) leftchild = node(1) 
        rightchild = node(node(1)+1)
           if ( rightchild .le. 0 ) rightchild = node(1) + 1
        newleaf = leaf(leaf(0))
           if ( newleaf .le. 0 ) newleaf = leaf(0)

        node(1) = node(1) + 2
        leaf(0) = leaf(0) + 1

        tree(2,j) = leftchild
        tree(1,leftchild) = j
        tree(2,leftchild) = 0
        tree(3,leftchild) = tree(3,j)

        tree(3,j) = rightchild
        tree(1,rightchild) = j
        tree(2,rightchild) = 0
        tree(3,rightchild) = newleaf
      end

**********************************************************************
      subroutine makeoffspring(allocate, n)
      integer nmxsub, n, nextnode, allocate
      parameter (nmxsub=1024)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub),node(2*nmxsub), leaf(0:nmxsub)
c ----------------------------------------------
c  if there is a node with children kill it
c  use depth first search to find a type 3 node
c  deletion order is inverse order of search
**********************************************************************
      integer j, jold, next, kind, allocated

        allocated = 0

        do while ( allocated .lt. allocate)
           jold = tree(1,n)
           j = n
           do while ( allocated .lt. allocate .and.
     &                j.ne.tree(1,n) .and. j.ne.tree(3,tree(1,n)) )
              next = nextnode(kind, j, jold)
              if ( abs(kind) .eq. 4 ) then
                 call makechildren(j)
                 allocated = allocated + 1
              end if
              jold = j
              j = next
           end do
        end do
      end

**********************************************************************
      subroutine freechildren(j)
      integer nmxsub, j
      parameter (nmxsub=1024)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub),node(2*nmxsub), leaf(0:nmxsub)
c ----------------------------------------------
c  parent release two children which SHOULD BE leaves
c  push 2 children nodes goes back to the node stack
c  left child's leaf goes back to the leaf stack
**********************************************************************
        leaf(0) = leaf(0) - 1
        leaf(leaf(0)) = tree(3,tree(3,j))

        node( node(1)-1 ) = tree(3,j)
        node( node(1)-2 ) = tree(2,j)
        node(1) = node(1) - 2

        tree(3,j) = tree(3,tree(2,j))
        tree(2,j) = 0
      end

**********************************************************************
      subroutine freeoffspring(released, n)
      integer nmxsub, n, nextnode, released
      parameter (nmxsub=1024)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub),node(2*nmxsub), leaf(0:nmxsub)
c ----------------------------------------------
c  if there is a node with children kill it
c  use depth first search to find a type 3 node
c  deletion order is inverse order of search
**********************************************************************
      integer j, jold, next, kind

        released = 0
        jold = tree(1,n)
        j = n

        do while ( j.ne.tree(1,n) .and. j.ne.tree(3,tree(1,n)) )
           next = nextnode(kind, j, jold)
           if ( abs(kind) .eq. 3 ) then
              call freechildren(j)
              released = released + 1
           end if
           jold = j
           j = next
        end do
      end

**********************************************************************
      function nextnode(kind,j,jold)
      integer nmxsub, nextnode, kind, j, jold
      parameter (nmxsub=1024)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
c -----------------------------------------------------------------
c  finding next node in "depth first search algorithm"
c  brife description of the algorithm:
c     1. parent with children            (downward pass)
c     3. parent visited both children    (upward pass)
c     4. a leaf                          (finest level pass)
c.   -3. If you returned root node, Stop
c.   -4. If you have only root node (a leaf), Stop
c -----------------------------------------------------------------
c   Example:n=1
c   2               jold  : 1*
c                   j     : 2
c                   next  : 1*
c                   kind  : -4
c   Example:n=4
c   2 - 3 - 5       jold  : 1* 2  3  5  6  3  4  7  8  4
c     |     6       j     : 2  3  5  6  3  4  7  8  4  2
c     - 4 - 7       next  : 3  5  6  3  4  7  8  4  2  1*
c           8       kind  : 1  1  4  4  3  1  4  4  3  -3
**********************************************************************
c ----- A leaf so return to parent, TYPE 4 or -4
      if ( tree(2,j).eq.0 ) then
         if ( j .eq. 2 ) then
            kind = -4
            nextnode = 1
         else if ( j .eq. tree(2,tree(1,j)) ) then
            kind = 4
            nextnode = tree(3,tree(1,j))
         else
            kind = 4
            nextnode = tree(1,j)
         end if
c ----- Move to parent from right child, TYPE 3 or -3
      else if ( jold.eq.tree(3,j) ) then
         if ( j.eq.2 ) then
            kind = -3
            nextnode = 1
         else if ( j .eq. tree(2,tree(1,j)) ) then
            kind = 3
            nextnode = tree(3,tree(1,j))
         else
            kind = 3
            nextnode = tree(1,j)
         end if
c ----- Move to left child, TYPE 1
      else
         kind = 1
         nextnode = tree(2,j)
      end if
      end

**********************************************************************
      subroutine listtree(n, located, leafid, leaflist, parentlist)
      integer nmxsub, n, located, nextnode
      parameter (nmxsub=1024)
      integer leafid(nmxsub), leaflist(nmxsub), parentlist(nmxsub)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
**********************************************************************
      integer j, jold, next, kind, parentlocate

        located = 0
        parentlocate = 0

        jold = tree(1,n)
        j = n
        do while ( j.ne.tree(1,n) .and. j.ne.tree(3,tree(1,n)) )
           next = nextnode(kind, j, jold)

           if ( kind .eq. 1 ) then
              parentlocate = parentlocate + 1
              parentlist( parentlocate ) = j
           else if ( abs(kind) .eq. 3 ) then
           else
              located = located + 1
              leaflist(located) = j
              leafid(located) = tree(3,j)
           end if

           jold = j
           j = next
        end do

        if ( located-1 .ne. parentlocate ) then
           print*,'Internal Error: located-1 .ne. parentlocate'
        end if
      end

**********************************************************************
      subroutine printtree(n)
      integer nmxsub, n, nextnode
      parameter (nmxsub=1024)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
**********************************************************************
      integer k, j, jold, depth_pr, depth, next, kind

        print'(i3,a,a,i3)',node(1)-2,' nodes used,   ',
     &            'next node is',node(node(1))
        print'(i3,a,a,i3)',leaf(0)-1,' leaves used,   ',
     &            'next leaf is',leaf(leaf(0))

        jold = tree(1,n)
        j = n
        depth = 0
        depth_pr = 0

        do while ( j.ne.tree(1,n) .and. j.ne.tree(3,tree(1,n)) )
           next = nextnode(kind, j, jold)

           if ( kind .eq. 1 ) then
              do k = depth_pr, depth-1
                 print'(a4,$)',':'
              end do
              depth = depth + 1
              depth_pr = depth
              print'(i4,$)',j
           else if ( abs(kind) .eq. 3 ) then
              depth = depth - 1
           else
              do k = depth_pr, depth-1
                 print'(a4,$)',':'
              end do
              print'(i4,a,i3,a)',j,'   (',tree(3,j),' )'
              depth_pr = 0
           end if

           jold = j
           j = next
        end do
      end

