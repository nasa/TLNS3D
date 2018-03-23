c------------------------------------------------------------------------------
c The TLNS3D code was developed to solve Reynolds-averaged Navier-Stokes
c Equations to simulate turbulent, viscous flows over three-dimensional
c configurations. A general multiblock grid approach is used to model
c complex configurations.  A multi-stage Runge-Kutta pseudo-time stepping
c scheme is coupled with residual smoothing and multigrid acceleration
c techniques to form an efficient algorithm for solving transonic viscous
c flows over aerodynamic configurations of practical interest.
c
c The TLNS3D framework is licensed under the Apache License, Version 2.0
c (the "License"); you may not use this application except in compliance
c with the License. You may obtain a copy of the License at
c http://www.apache.org/licenses/LICENSE-2.0. 

c Unless required by applicable law or agreed to in writing, software
c distributed under the License is distributed on an "AS IS" BASIS,
c WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
c See the License for the specific language governing permissions and
c limitations under the License.
c----------------------------------------------------------------------------------------------
c
      subroutine disdr(imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,
     .                 imp2,jmp2,kmp2,m1cc,m3np,mbloc,nbloc,
     .                 mx1dwk,mxdstf,mxdsti,mxsizc,mxsiz3n,
     .                 imap,nseg,msegt,msegn,mgrlev,nlev,
     .                 x,wk3d,xs,smin,wk3d5,iwrk,nsurf,nodes   )
c
c     object : To serve as a driver routine for computing
c              distance to closest solid surface
c             (note that search region involves information
c              across all the blocks for all solid segments)
c
c     initial coding : by v.n.vatsa of NASA Langley (March, 1993 )
c
c     modified from disdr routine for distributed computing by
c     March 1995. Added xs array to argument
c     to provide space for surface points.
c

#if defined(BUILD_PVM) || defined(BUILD_MPI)
#ifdef BUILD_PVM
#     include "fpvm3.h"
#  ifdef TLN_DBL
#    define RTYPE REAL8
#  else
#    define RTYPE REAL4
#  endif
#else
      include "mpif.h"
#  ifdef TLN_DBL
#    define RTYPE MPI_DOUBLE_PRECISION
#  else
#    define RTYPE MPI_REAL
#  endif
      dimension istat(MPI_STATUS_SIZE)
#endif
#     include "tags.h"
#endif
      character*128  errmsg
      common/dstcom/ errmsg,master,myrank,mycomm,ndlist,nnodes
      dimension      nodes(1)
c
      dimension     im  (mgrlev,mbloc), jm  (mgrlev,mbloc),
     .              km  (mgrlev,mbloc), imp1(mgrlev,mbloc),
     .              jmp1(mgrlev,mbloc), kmp1(mgrlev,mbloc),
     .              imp2(mgrlev,mbloc), jmp2(mgrlev,mbloc),
     .              kmp2(mgrlev,mbloc)
c
      dimension     imap(msegt,msegn,mbloc,mgrlev), nseg(mbloc)
c
      dimension     m1cc(mgrlev,mbloc), m3np(mgrlev,mbloc)
c
c
      dimension x(mxsiz3n),wk3d(mxsizc),smin(mxsizc)
      dimension xs(mx1dwk*mx1dwk*3)
      dimension wk3d5(mxdstf), iwrk(mxdsti)
      dimension tim(2)
      integer vlist,bbdef,surf
      common/tresv/ timrv
c
c--------------------  mapping function description ---------------------
c
c      imap  : mapping function containing topological information
c      msegt : maximum no. of types of operations/boundary conditons
c              
c      msegn : maximum no. of segments permitted on all faces of a block
c              
c
c      imap(1 , ) : specifies boundary/operation type
c      imap(2 , ) : specifies face number (1-6)
c      imap(3 , ) : specifies beginning of direction 1
c      imap(4 , ) : specifies end       of direction 1
c      imap(5 , ) : specifies beginning of direction 2
c      imap(6 , ) : specifies end       of direction 2
c
c--------------------  boundary/operation type description ---------------
c
c      nbctype    = imap(1 , )
c                 = 0  corresponds to an inner cut that maps to same block
c                 = 1  corresponds to a cut that maps to another block
c                 = 2  corresponds to a slip (inviscid) wall
c                 = 3  corresponds to a noslip (viscous) wall
c
c***************************************************************************
c
      tim(2) = 0.0
      iss = 0
      igrid   = 1
c
c     initialize the distance function
c
      do 10 i=1,mxsizc
      smin(i) = 1.e+20
   10 continue
c
c	print *,'mxdstf,mxdsti,nsurf: ',mxdstf,mxdsti,nsurf
      call initi(mxdsti)
      call initf(mxdstf)
      surf = ifalloc(3*nsurf)
      ntri = iialloc(nsurf)
      iptri = iialloc(8*nsurf)
      isurf=0
c	print *,'loop 1000'
c
c---  begin outer loop on the blocks  ----------------------
c
      do 1000 ibloc = 1,nbloc
c     
         ns = nseg(ibloc)
c     
c---- begin outer loop on segments  -----------------------
c     
         do 100 iseg = 1,ns
c     
            nbctype  =  imap(1 ,iseg ,ibloc,igrid)
c     
            if (nbctype.eq.3.or.nbctype.eq.102) then
               nface  =  imap(2 ,iseg ,ibloc,igrid)
               n1beg  =  imap(3 ,iseg ,ibloc,igrid)
               n1end  =  imap(4 ,iseg ,ibloc,igrid)
               n2beg  =  imap(5 ,iseg ,ibloc,igrid)
               n2end  =  imap(6 ,iseg ,ibloc,igrid)
c     
               n1pts  =  n1end - n1beg + 1
               n2pts  =  n2end - n2beg + 1
               npts   =  n1pts * n2pts * 3
c     
c------on node where ibloc is local
c     
# if defined BUILD_MPI
               if (nodes(ibloc)-1.eq.myrank) then
# else
               if (nodes(ibloc).eq.myrank) then
# endif
c     
c     get surface points into work array
c     
c	print *,'getpts: ',ibloc,iseg
                  call getpts (imn,jmn,kmn,
     .               im(igrid,ibloc),jm(igrid,ibloc),
     .               km(igrid,ibloc),imp1(igrid,ibloc),
     .               jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .               x(m3np(igrid,ibloc)),mx1dwk,nface,
     .               n1beg,n1end,n2beg,n2end,xs)
c     
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c     
c     send surface points to other nodes
c     
                  do inode = 1,nnodes
# if defined BUILD_MPI
                     if (nodes(ndlist+inode)-1.ne.myrank) then
# else
                     if (nodes(ndlist+inode).ne.myrank) then
# endif
#ifdef BUILD_PVM
                        call PVMFpsend (nodes(ndlist+inode),TAG_DIST,
     .                        xs,npts,RTYPE,ierr)
#else 
                        call MPI_Send (xs,npts,RTYPE,
     .                       nodes(ndlist+inode)-1,TAG_DIST,
     .                       mycomm,ierr)
#endif
                     endif
                  enddo
c     
c------on other nodes, receive surface points
c     
               else
#ifdef BUILD_PVM
                  call PVMFprecv (nodes(ibloc),TAG_DIST,
     .                    xs,npts,RTYPE,
     .                    itid,itag,ilen,ierr)
#else 
                  call MPI_Recv (xs,npts,RTYPE,
     .                   nodes(ibloc)-1,TAG_DIST,
     .                   mycomm,istat,ierr)
#endif
#endif
               endif
c     
c	print *,'collect_surf: ',ibloc,iseg
               call collect_surf (xs,
     .            n1beg ,n1end ,n2beg , n2end, 
     .            nsurf, wk3d5(surf), isurf,
     .            iwrk(ntri), iwrk(iptri) )
c	print *,'done: ',ibloc,iseg
c     
            endif
c     
c---- end loop on segment
 100     continue
c---- end loop on blocks
 1000 continue
c	print *,'end loop 1000'
c
c to help locate identical points on boundaries, we sort with
c respect to x coordinate.  Better to use a stable sort on
c x,y,z and ntri consecutively, but this is good enough
c for time being.
c
      call sort_x(nsurf,wk3d5(surf),iwrk(ntri),iwrk(iptri),
     .     wk3d5,iwrk)
c
c Now run Mac Ice routines for making boxes that surround
c groups of viscous surface points
c
c      time1=second()
      time1=0.0
      minbox=sqrt(real(nsurf))
      minbox=max(minbox,50)
      nbb = 3*nsurf/minbox
      bbdef = ifalloc(6*nbb)
      ipv = iialloc(2*nbb)
      vlist = iialloc(nsurf)
c	print *,'call makebb'
      call makebb(minbox,10,
     .     nsurf,wk3d5(surf),wk3d5(surf+nsurf),wk3d5(surf+2*nsurf),
     .     nbb, ntotv,
     .     wk3d5(bbdef),iwrk(ipv),iwrk(vlist),iwrk)
c      time2=second()
      time2=1.0
c      print *,"time to set up boxes=",time2-time1
c
c------ on all nodes - begin inner loop on the blocks
c
c	print *,'loop 1010'
      do 1010 iibloc = 1,nbloc
c
# if defined BUILD_MPI
         if (nodes(iibloc)-1.eq.myrank)
# else
         if (nodes(iibloc).eq.myrank)
# endif
     .      call calc_dist(imn,jmn,kmn,
     .        im  (igrid,iibloc),jm  (igrid,iibloc),
     .        km  (igrid,iibloc),
     .        imp1(igrid,iibloc),jmp1(igrid,iibloc),
     .        kmp1(igrid,iibloc),
     .        imp2(igrid,iibloc),jmp2(igrid,iibloc),
     .        kmp2(igrid,iibloc),
     .        x(m3np(igrid,iibloc)),smin(m1cc(igrid,iibloc)),
     .        nsurf, wk3d5(surf), nbb, wk3d5(bbdef), iwrk(ipv),
     .        iwrk(vlist),
     .        iwrk(ntri), iwrk(iptri) ,wk3d5,iwrk)
 1010 continue
c	print *,'end loop 1010'
c
c      computation of distance to node points is now complete
c      on the finest grid
c
      call ifree(nsurf)
      call ifree(2*nbb)
      call ffree(6*nbb)
      call ffree(3*nsurf)
      call ifree(8*nsurf)
      call ifree(nsurf)
c	print *,'done with free'
c
c     computation of distance to node points is now complete
c     on the finest grid
c
c     generate distance to nodes on coarse grids 
c
      do 3120 ibloc = 1,nbloc
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) then
# else
      if (nodes(ibloc).eq.myrank) then
# endif
      do 3105 igrid = 1,nlev-1
c
      call distcg ( imn,jmn,kmn,
     . im  (igrid  ,ibloc), jm  (igrid  ,ibloc), km  (igrid  ,ibloc),
     . imp1(igrid  ,ibloc), jmp1(igrid  ,ibloc), kmp1(igrid  ,ibloc),
     . imp2(igrid  ,ibloc), jmp2(igrid  ,ibloc), kmp2(igrid  ,ibloc),
     . im  (igrid+1,ibloc), jm  (igrid+1,ibloc), km  (igrid+1,ibloc),
     . imp1(igrid+1,ibloc), jmp1(igrid+1,ibloc), kmp1(igrid+1,ibloc),
     . imp2(igrid+1,ibloc), jmp2(igrid+1,ibloc), kmp2(igrid+1,ibloc),
     . smin(m1cc(igrid,ibloc)), smin(m1cc(igrid+1,ibloc))           )
c
 3105 continue
      endif
 3120 continue
c
c     generate distance to cell-centers 
c
      do 3150 ibloc = 1,nbloc
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) then
# else
      if (nodes(ibloc).eq.myrank) then
# endif
      do 3110 igrid = 1,nlev
c
      call distcc ( imn,jmn,kmn,
     . im  (igrid  ,ibloc), jm  (igrid  ,ibloc), km  (igrid  ,ibloc),
     . imp1(igrid  ,ibloc), jmp1(igrid  ,ibloc), kmp1(igrid  ,ibloc),
     . imp2(igrid  ,ibloc), jmp2(igrid  ,ibloc), kmp2(igrid  ,ibloc),
     . smin(m1cc(igrid,ibloc)), wk3d                                 )
c
 3110 continue
      endif
 3150 continue
c
      return
      end
