      module injectorData
      real*8 pmaxb, pminb
      real*8 pmaxh, pminh
      integer, allocatable, dimension(:):: seednodes
      real*8, allocatable, dimension(:,:) ::activenodes
      logical injstart
      end module
      
      subroutine identifyInjectorPatch(iBC, BC, mark,x,injectorDirect
     &      ,markpatch)
      use injectorData
      include "common.h"
      include "mpif.h"
      real*8 BC(nshg,ndofBC), x(numnp,nsd)
      integer iBC(nshg)
      logical mark(nshg)
      logical markpatch(nshg)
      integer injectorDirect

      real*8 velmag
      integer node1, node2
      real*8 maxb(1), minb(1)
      real*8 maxh(1), minh(1)
      real*8 dist
      integer count
      
      mark = .false.
      markpatch = .false.
      do i=1,nshg
         velmag = sqrt(BC(i,3)**2.0 + BC(i,4)**2.0 + BC(i,5)**2.0)
         if(iBC(i) .eq. 56 .and. velmag .gt. 0.0)then
!     if iBC == 56 means that velocity is specfied on this particular node
!     in simmodeler or chef's Simplified_SPJ_file or temporal BC (this remains to be checked)
!     Thus this ensures that the node belongs to inlet patch
            mark(i) = .true.
            markpatch(i) = .true.
         endif         
      enddo

      
!     Determine the direction of inlet patch normal
      node1 = 0
      node2 = 0
      if(any(mark(:)))then
         do i=1,nshg
            if(mark(i) .eq. .true. .and. node1 .ne. 0)then
               node2 = i
            elseif(mark(i) .eq. .true. .and. node1 .eq. 0)then
               node1 = i
            elseif(node1 .ne. 0 .and. node2 .ne. 0)then
               exit
            endif
         enddo
         if(node1 .ne. 0 .and. node2 .ne. 0)then
            if(abs(x(node1,1)-x(node2,1)) .lt. 1e-8)then
               injectorDirect = 1
            elseif(abs(x(node1,2)-x(node2,2)) .lt. 1e-8)then
               injectorDirect = 2
            else
               injectorDirect = 3
            endif
         endif
      endif

      if(myrank .ne. master .and. injectorDirect .ne. 0)then
         call MPI_Send(injectorDirect,1,MPI_INTEGER,master,myrank
     &        ,MPI_COMM_WORLD,ierr)
      elseif(myrank .eq. master)then
         call MPI_Recv(injectorDirect,1,MPI_INTEGER,MPI_ANY_SOURCE,
     &        MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS,ierr)
      endif

      call MPI_Bcast(injectorDirect,1,MPI_INTEGER,
     &                    master,MPI_COMM_WORLD,ierr)


      
!Following section to determine the spatial limits of injector patch
      if(injectorDirect .eq. 1)then
         maxb = maxval(x(:,2))
         minb = minval(x(:,2))
         maxh = maxval(x(:,3))
         minh = minval(x(:,3))
      elseif(injectorDirect .eq. 2)then
         maxb = maxval(x(:,1))
         minb = minval(x(:,1))
         maxh = maxval(x(:,3))
         minh = minval(x(:,3))
      elseif(injectorDirect .eq. 3)then
         maxb = maxval(x(:,1))
         minb = minval(x(:,1))
         maxh = maxval(x(:,2))
         minh = minval(x(:,2))
      endif

      call MPI_ALLREDUCE(maxb(1), pmaxb, 1,
     &     MPI_DOUBLE_PRECISION, MPI_MAX,
     &     MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(minb(1), pminb, 1,
     &     MPI_DOUBLE_PRECISION, MPI_MIN,
     &     MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(maxh(1), pmaxh, 1,
     &     MPI_DOUBLE_PRECISION, MPI_MAX,
     &     MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(minh(1), pminh, 1,
     &     MPI_DOUBLE_PRECISION, MPI_MIN,
     &     MPI_COMM_WORLD, ierr)
      


!     Unmark the points for which bubble will be too close to edge
      dist = 0.0
      count = 0
      k=0
      if(any(mark(:)))then
         do i=1,nshg
            if(mark(i) .eq. .true.)then
               if(injectorDirect .eq. 1)then
                  dist = min(abs(pmaxb-x(i,2)),abs(pminb-x(i,2)),
     &                 abs(pmaxh-x(i,3)),abs(pminh-x(i,3)))
               elseif(injectorDirect .eq. 2)then
                  dist = min(abs(pmaxb-x(i,1)),abs(pminb-x(i,1)),
     &                 abs(pmaxh-x(i,3)),abs(pminh-x(i,3)))
               elseif(injectorDirect .eq. 3)then
                  dist = min(abs(pmaxb-x(i,1)),abs(pminb-x(i,1)),
     &                 abs(pmaxh-x(i,2)),abs(pminh-x(i,2)))
               endif

               if(dist .lt. 1.25*rinj)then
                  mark(i) = .false.
               endif
            endif
         enddo

         do i=1,nshg
            if(mark(i) .eq. .true.) count=count+1
         enddo

         if(count .gt. 0)then
            allocate(seednodes(count))
            do i=1,nshg
               if(mark(i) .eq. .true.)then
                  k=k+1
                  seednodes(k) = i
               endif
            enddo
         endif
      endif

      !if(any(mark(:)))write(*,*)"proc is marked ",myrank
      

      
      end

      subroutine addnode()
!The new node is added at the top of the list
!Previous ones are pushed down
      use injectorData
      include "common.h"
      include "mpif.h"

      real*8, allocatable, dimension(:,:) ::tempnodes

      integer newsize

      newsize = size(activenodes,1) + 1
      if(newsize .gt. 1)then
         allocate(tempnodes(newsize,6))
         tempnodes(2:newsize,:) = activenodes(1:newsize-1,:)
         deallocate(activenodes)
         allocate(activenodes(newsize,6))
         activenodes = tempnodes
         deallocate(tempnodes)  
      else
         allocate(activenodes(1,6))
      endif
      activenodes(1,:) = 0.0
      end

      subroutine deletenode(inode)
      use injectordata
      include "common.h"
      include "mpif.h"

      integer inode

      integer isize, k
      real*8, allocatable, dimension(:,:) ::tempnodes
      
      isize = size(activenodes,1)
      k=1
      if(isize .eq. 1)then
         deallocate(activenodes)
      else
         allocate(tempnodes(isize-1,6))
         do i=1,isize
            if(i .ne. inode)then
               tempnodes(k,:) = activenodes(i,:)
               k=k+1
            endif
         enddo
         deallocate(activenodes)
         allocate(activenodes(isize-1,6))
         activenodes = tempnodes
         deallocate(tempnodes)
      endif
      end


      
      subroutine injector(y, iBC, BC, x, mark,injectorDirect,markpatch,
     &      injectortime1, injectortime2)

!The injector patch will most likely be the patch
!where velocity is specified in simmodeler
!(or in chef - have to check this)
!(or if time varying boundary conditions are specified)
!This subroutine works on this assumption i.e., seed points are nodes where iBC = 56
      
      use injectorData
      include "common.h"
      include "mpif.h"
      real*8 y(nshg,ndof), BC(nshg,ndofBC),
     &     x(numnp,nsd)
      integer iBC(nshg)
      integer injectorDirect
      integer injectortime1, injectortime2

      real*8 velmag, temp
      logical mark(nshg),markpatch(nshg)
      integer chosen, itemp
      real*8 tperiodinj         !input variable
      real*8 position(3)
      real*8 velocity(3)
      integer ilast
      real*8 dist, dist1
!     integer, allocatable, dimension(:):: procseednum
      integer procseednum(numpe)
      real*8, allocatable, dimension(:) ::temphi
      
      if(time .gt. 1.0e-12) injstart = .true.
      
      chosen = 0
      injectortime1 = floor(time/tperiodinj)
      if(injstart .eq. .false. .or. injectortime1.ne.injectortime2)then
         injstart = .true.
         injectortime2 = injectortime1
         call addnode()
         ilast = size(activenodes,1)
!Choose a random seed point from this processor, if any
!Ensure that the chosen point is sufficiently far from prevoius seed points
         if(any(mark(:)))then
            if(ilast .gt. 1)then
               do i=1,10 !10 iterations should be enough to choose a valid seed
                  call random_number(temp)
                  temp = 1.0 + (size(seednodes) - 1)*temp
                  chosen = floor(temp)
                  dist = 100.0
                  do j=1,ilast-1
                     if(injectorDirect .eq. 1)then
                        dist1 = sqrt(
     &              (x(seednodes(chosen),2)-activenodes(j,2))**2.0 +
     &              (x(seednodes(chosen),3)-activenodes(j,3))**2.0)
                     elseif(injectorDirect .eq. 2)then
                        dist1 = sqrt(
     &              (x(seednodes(chosen),1)-activenodes(j,1))**2.0 +
     &              (x(seednodes(chosen),3)-activenodes(j,3))**2.0)
                     elseif(injectorDirect .eq. 3)then
                        dist1 = sqrt(
     &              (x(seednodes(chosen),2)-activenodes(j,2))**2.0 +
     &              (x(seednodes(chosen),1)-activenodes(j,1))**2.0)
                     endif
                     if(dist1 .lt. dist)dist = dist1
                  enddo
                  if(dist .gt. 2.5*rinj)exit
               enddo
            else
               call random_number(temp)
               temp = 1.0 + (size(seednodes) - 1)*temp
               chosen = floor(temp)
               !write(*,*)" proc chose ",myrank,chosen
            endif

            
         endif

!Assemble list of all chosen points on master (may include 0)
!         if(myrank .eq. master)then
!            allocate(procseednum(numpe))
!         endif
         !if(myrank .eq. master)write(*,*)"reached here1",master
         call MPI_Gather(chosen,1,MPI_INTEGER,procseednum,1,
     &        MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

         !if(myrank .eq. master)write(*,*)"reached here2",master
!Let the master choose a processor from whoch to select the next seed
         if(myrank .eq. master)then
            call random_number(temp)
            temp = 1.0 + (numpe - 1)*temp
            itemp = floor(temp)
            do i=1,numpe
               if(procseednum(itemp) .ne. 0)then
                  exit
               else
                  if(itemp .eq. numpe) itemp=0
                  itemp = itemp +1
               endif
            enddo
            itemp = itemp-1
!            deallocate(procseednum)
         endif
         !if(myrank .eq. master)write(*,*)"reached here3",master
!The master now broadcasts the elected processor to all
         call MPI_Bcast(itemp,1,MPI_INTEGER,
     &        master,MPI_COMM_WORLD,ierr)
         !if(myrank .eq. master)write(*,*)"reached here4",master,itemp
!The elected processor sends info about the seed to all
         if(myrank .eq. itemp)then
           activenodes(1,1:3)=x(seednodes(chosen),1:3) !position of seed
           activenodes(1,4:6)=BC(seednodes(chosen),3:5) !velocity of seed
         endif
         call MPI_Bcast(activenodes(1,:),6,
     &        MPI_DOUBLE_PRECISION,itemp,MPI_COMM_WORLD,ierr)
         !if(myrank .eq. master)write(*,*)"reached here5",master
!     Adjust position of seed such that bubble lies just outside patch (on all procs)
!Correction needed here. It is assumed that inlet patches are:
!left for x-direction
!bottom for y-direction
!back for z-direction
         if(injectorDirect .eq. 1)then
            activenodes(1,1) = activenodes(1,1)-rinj
         elseif(injectorDirect .eq. 2)then
            activenodes(1,2) = activenodes(1,2)-rinj
         elseif(injectorDirect .eq. 3)then
            activenodes(1,3) = activenodes(1,3)-rinj
         endif

         if(myrank .eq. master)then
            write(*,*)"Seed point selected from proc ",itemp
            write(*,*)activenodes(1,1), activenodes(1,2),
     &           activenodes(1,3)
            write(*,*)activenodes(1,4), activenodes(1,5),
     &            activenodes(1,6)
         endif
      endif

!     Modify LS field on inlet patch
      ilast = size(activenodes,1)
      if(ilast .gt. 0)then
         if(any(markpatch(:)))then
            allocate(temphi(ilast))            
            do i=1,nshg
               if(markpatch(i) .eq. .true.)then
                  do j=1,ilast
                     position(1:3) = activenodes(j,1:3) + time
     &                 *activenodes(j,4:6)
                     temphi(1) = sqrt((x(i,1)-position(1))**2.0 +
     &                    (x(i,2)-position(2))**2.0 +
     &                    (x(i,3)-position(3))**2.0)-rinj
                  enddo
                  !Will work as long as seed points do not overlap
                  y(i,6) = minval(temphi)
               endif
            enddo
            deallocate(temphi)
         endif
      endif

!Find and remove the seeds for which bubble is completely inside
      ilast = size(activenodes,1)
      if(ilast .gt. 0)then
         do i=1,ilast
            position(1:3) = activenodes(i,1:3) + time
     &           *activenodes(i,4:6)
            if(injectDirect .eq. 1)then
               if(abs(position(1)-activenodes(i,1)) .gt. 2.0*rinj)then
                  call deletenode(i)
               endif
            elseif(injectDirect .eq. 2)then
               if(abs(position(2)-activenodes(i,2)) .gt. 2.0*rinj)then
                  call deletenode(i)
               endif
            elseif(injectDirect .eq. 3)then
               if(abs(position(3)-activenodes(i,3)) .gt. 2.0*rinj)then
                  call deletenode(i)
               endif
            endif
         enddo
      endif
      
      
 100  FORMAT(1I10,5X,E10.5,5X,E10.5,5X,E10.5,5X,E10.5)
      end
