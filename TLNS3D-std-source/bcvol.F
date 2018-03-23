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
      subroutine bcvol(imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,
     .                 imp2,jmp2,kmp2,m1cc,
     .                 mgrlev,igrid,isoln,mbloc,nbloc,
     .                 mxsizc,mx1dwk,
     .                 imap,msegt,msegn,nseg,
     .                 vol,wk2da,
     .                 mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .                 mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .                 lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .                 kpitmb1,kpitmb2,frc,
     .                 m1pch1,m1pch2,lswpchb,lswpche,
     .                 ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes )
c-----------------------------------------------------------------------
c
c     object : To set ghost-cell values for cell volumes
c              based on topological mapping funcions, 'imap'
c
c     initial coding : by v.n.vatsa (March, 1990)
c     latest mods    : To accommodate patched inter-faces (Vatsa: Nov. 1994)
c                      nbctyp < 0 for patched boundaries
c
c     modified for distributed computing : April, 1995
c
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
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      data iitmbeg /1/
c
      dimension     im  (mgrlev,mbloc), jm  (mgrlev,mbloc),
     .              km  (mgrlev,mbloc), imp1(mgrlev,mbloc),
     .              jmp1(mgrlev,mbloc), kmp1(mgrlev,mbloc),
     .              imp2(mgrlev,mbloc), jmp2(mgrlev,mbloc),
     .              kmp2(mgrlev,mbloc),
     .              imap(msegt,msegn,mbloc), nseg(mbloc)
c
      dimension     m1cc(mgrlev,mbloc)
c
      dimension     vol(mxsizc)
c
c     array declaration for temporary work-space
c     this needs to be at least size (mx1dwk,mx1dwk,2)
c
      dimension   wk2da(1)
c
c
c     patcher related information/arrays
c
      common/facetr/ ifacetr(6),ipatchg
c
      dimension   ntpchcb(mbloc,mgrlev), npchcbf(6,mbloc,mgrlev),
     .            m1pch1 (mbloc,mgrlev), m1pch2 (mbloc,mgrlev),
     .          ipatchc(mxtpchc), jpatchc(mxtpchc),
     .          kpatchc(mxtpchc)
c
c
      dimension nswpchb  (mbloc,mgrlev),   npchitm(mxpchs,mgrlev),
     .          lswpchb  (mbloc,mgrlev),   
     .          lswpche  (mbloc,mgrlev),
     .          lspchb1 (mxpchs,mgrlev),   lspchf1(mxpchs,mgrlev),
     .          lspchb2 (mxpchs,mgrlev),   lspchf2(mxpchs,mgrlev),
     .          ipitmb1(mpchitm),   ipitmb2(mpchitm),
     .          jpitmb1(mpchitm),   jpitmb2(mpchitm),
     .          kpitmb1(mpchitm),   kpitmb2(mpchitm),
     .          frc    (mpchitm),
     .          ipitmbs(mpchitm),   jpitmbs(mpchitm),
     .          kpitmbs(mpchitm),   iitmsa (mxpchs,mgrlev)
                            
c
c
c
c---  lpchs is cumulative value of surface segments with patched b.c.
c---  on current grid level
c---  litmbeg is starting no. in cumulative patched items at "lpchs" patch
c---  lpchcb  is starting (global) location for patched cells on a block
c
      lpchs    = 0
c     litmbeg  = 1
      lpchcb   = 1
c
c     nghost is number of ghost cells layers
c     2 on finest grid, 1 on coarse grids
c
      nghost   = 1
      if (isoln.eq.igrid) nghost = 2
c
c-------  begin outer loop on the blocks for interface boundaries --------
c
      do 1000 ibloc = 1,nbloc
c
       ns       =  nseg(ibloc)
c
c----------  begin outer loop on the segments  ----------------------------
c
       do 100 iseg = 1,ns
c
       nbctype  =  imap(1 ,iseg ,ibloc)
c
       if (nbctype.eq.1 .or. nbctype.eq.0) then
c
         nface   =  imap(2 ,iseg ,ibloc)
         n1beg   =  imap(3 ,iseg ,ibloc)
         n1end   =  imap(4 ,iseg ,ibloc)
         n2beg   =  imap(5 ,iseg ,ibloc)
         n2end   =  imap(6 ,iseg ,ibloc)
         nblocs  =  imap(7 ,iseg ,ibloc)
         nfaces  =  imap(8 ,iseg ,ibloc)
         n1begs  =  imap(9 ,iseg ,ibloc)
         n1ends  =  imap(10,iseg ,ibloc)
         n2begs  =  imap(11,iseg ,ibloc)
         n2ends  =  imap(12,iseg ,ibloc)
c
         n1cnt   =  iabs (n1ends - n1begs) + 2
         n2cnt   =  iabs (n2ends - n2begs) + 2
         if (nfaces.lt.0) then
c
c          source and target directions do not match
c
           ncnt  =  n1cnt
           n1cnt =  n2cnt
           n2cnt =  ncnt
         endif
         ncnt    =  n1cnt * n2cnt * nghost
c
c----------------  block-interface/inner cut    ---------------------
c
c        get ghost cell volumes from source block
c
# if defined BUILD_MPI
         if (nodes(nblocs)-1.eq.myrank) then
# else
         if (nodes(nblocs).eq.myrank) then
# endif
c
           call bccutget (imn,jmn,kmn,
     .     im  (igrid,nblocs), jm  (igrid,nblocs), km  (igrid,nblocs),
     .     imp2(igrid,nblocs), jmp2(igrid,nblocs), kmp2(igrid,nblocs),
     .     vol(m1cc(igrid,nblocs)),
     .     nfaces,n1begs,n1ends,n2begs,n2ends,
     .     nghost,wk2da                                 )
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c          if target is not local, send ghost cell volumes to node
c
# if defined BUILD_MPI
           if (nodes(ibloc)-1.ne.myrank) then
# else
           if (nodes(ibloc).ne.myrank) then
# endif
#ifdef BUILD_PVM
             call PVMFpsend (nodes(ibloc),TAG_VOL,
     .                       wk2da,ncnt,RTYPE,ierr)
#else
             call MPI_Send (wk2da,ncnt,RTYPE,
     .                      nodes(ibloc)-1,TAG_VOL,mycomm,ierr)
#endif
           endif
#endif
c
         endif
c
c        update ghost cell volumes on target block
c
# if defined BUILD_MPI
         if (nodes(ibloc)-1.eq.myrank) then
# else
         if (nodes(ibloc).eq.myrank) then
# endif
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c          receive ghost cell volumes from node if not already local
c
# if defined BUILD_MPI
           if (nodes(nblocs)-1.ne.myrank) then
# else
           if (nodes(nblocs).ne.myrank) then
# endif
#ifdef BUILD_PVM
             call PVMFprecv (nodes(nblocs),TAG_VOL,
     .                       wk2da,ncnt,RTYPE,
     .                       itid,itag,ilen,ierr)
#else
             call MPI_Recv (wk2da,ncnt,RTYPE,
     .                      nodes(nblocs)-1,TAG_VOL,mycomm,istat,ierr)
#endif
           endif
#endif
c
           call bccutset (imn,jmn,kmn,
     .     im  (igrid,ibloc),jm  (igrid,ibloc),km  (igrid,ibloc),
     .     imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .     vol(m1cc(igrid,ibloc )),
     .     nface,n1beg,n1end,n2beg,n2end,
     .     nghost,wk2da,n1cnt,n2cnt                     )
         endif
c
       endif
c
c----      end loop on segments
  100  continue
c
c---------- initialize variables on patched boundaries
c
       if (ipatchg.eq.0) go to 101
       if (ntpchcb(ibloc,igrid).le.0) go to 101
c
# if defined BUILD_MPI
       if (nodes(ibloc)-1.eq.myrank) then
# else
       if (nodes(ibloc).eq.myrank) then
# endif
c
          if( (m1pch1(ibloc,igrid)+ntpchcb(ibloc,igrid)).gt.mxtpchc) 
     .    then
             write (iwrit,'(2x,"dimension conflict for mxtpchc "/)')
             write (iwrit,'(2x,"mxtpchc m1pch1 ntpchc igrid ibloc"/)')
             write (iwrit,'(2x,5i7)') mxtpchc,m1pch1(ibloc,igrid),
     .       ntpchcb(ibloc,igrid),igrid,ibloc
             write (iwrit,'(2x,"stop in bcvol  sending ipatchc"/)')
c
             call ERREXIT (nodes)

          endif
         call initpgr (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .    vol(m1cc(igrid,ibloc)),
     .    npchcbf(1,ibloc,igrid),       ipatchc(m1pch1(ibloc,igrid)),
     .    jpatchc(m1pch1(ibloc,igrid)), kpatchc(m1pch1(ibloc,igrid)),
     .    igrid, isoln)
c     else
c       do iface=1,6
c         if (npchcbf(iface,ibloc,igrid).gt.0)
c    .      lpchcb = lpchcb + npchcbf(iface,ibloc,igrid)
c       enddo
      endif
c
c------------- patched interface ---------------------------------------
c
c
      litmbeg = m1pch2(ibloc,igrid)
      do 120 lpchs=lswpchb(ibloc,igrid)+1,lswpche(ibloc,igrid)
cvn   iitmbeg = litmbeg
ccvn  iitmsa(lpchs,igrid) = iitmbeg
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
#else
c*/IN-CORE/
      iitmbeg = litmbeg
#endif
c

c
c      convert face numbers to tlns3d's convention
c      note: ibloc1 and ibloc are equal (this was already checked)
c
       ibloc1 = lspchb1(lpchs,igrid)
       iface1 = ifacetr(lspchf1 (lpchs,igrid))
       ibloc2 = lspchb2 (lpchs,igrid)
       iface2 = ifacetr(lspchf2 (lpchs,igrid))
       litems = npchitm (lpchs,igrid)
       ncnt   = litems * nghost
c
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c      send the source indices info. to node containing ibloc2
c
# if defined BUILD_MPI
         if (nodes(ibloc)-1.eq.myrank.and.nodes(ibloc2)-1.ne.myrank)
# else
         if (nodes(ibloc).eq.myrank.and.nodes(ibloc2).ne.myrank)
# endif
     .     then
#ifdef BUILD_PVM
c
      write (iwrit,'("PVM not supported for patched grids")')
      call ERREXIT (nodes)

#else
c         call MPI_Send (litmbeg,1,MPI_INTEGER,
c    .                   nodes(ibloc2)-1,TAG_FLOW,
c    .                   mycomm,ierr)
c
          if( (litmbeg+litems).gt.mpchitm) then
             write (iwrit,'(2x,"dimension conflict for mpchitm "/)')
             write (iwrit,'(2x,"mpchitm litmbeg litems igrid ibloc"/)')
             write (iwrit,'(2x,5i7)')mpchitm,litmbeg,litems,igrid,ibloc
             write (iwrit,'(2x,"stop in  bcvol before sending itms"/)')
c
             call ERREXIT (nodes)

          endif
c
          call MPI_Send (ipitmb2(litmbeg),litems,MPI_INTEGER,
     .                   nodes(ibloc2)-1,TAG_FLOW,
     .                   mycomm,ierr)
          call MPI_Send (jpitmb2(litmbeg),litems,MPI_INTEGER,
     .                   nodes(ibloc2)-1,TAG_FLOW,
     .                   mycomm,ierr)
          call MPI_Send (kpitmb2(litmbeg),litems,MPI_INTEGER,
     .                   nodes(ibloc2)-1,TAG_FLOW,
     .                   mycomm,ierr)
#endif
         endif
#endif
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c      put the source indices info. at the node containing ibloc2
c
# if defined BUILD_MPI
         if (nodes(ibloc2)-1.eq.myrank) then
# else
         if (nodes(ibloc2).eq.myrank) then
# endif
cvn          iitmbeg = 1
# if defined BUILD_MPI
             if (nodes(ibloc)-1.eq.myrank) then
# else
             if (nodes(ibloc).eq.myrank) then
# endif
c               iitmbeg = litmbeg
c               iibeg = iitmbeg
c               iiend = iitmbeg + litems - 1
c               do ii = iibeg,iiend
c                 ipitmbs(ii) = ipitmb2(ii)
c                 jpitmbs(ii) = jpitmb2(ii)
c                 kpitmbs(ii) = kpitmb2(ii)
c               enddo
             else
#ifdef BUILD_PVM
c
      write (iwrit,'("PVM not supported for patched grids")')
      call ERREXIT (nodes)

#else
c         call MPI_Recv (iitmbeg,1,MPI_INTEGER,
c    .                   nodes(ibloc1)-1,TAG_FLOW,
c    .                   mycomm,istat,ierr)
cvn (3-14-97)
cvn       iitmbeg = 1
          call MPI_Recv (ipitmbs(iitmbeg),litems,MPI_INTEGER,
     .                   nodes(ibloc1)-1,TAG_FLOW,
     .                   mycomm,istat,ierr)
          call MPI_Recv (jpitmbs(iitmbeg),litems,MPI_INTEGER,
     .                   nodes(ibloc1)-1,TAG_FLOW,
     .                   mycomm,istat,ierr)
          call MPI_Recv (kpitmbs(iitmbeg),litems,MPI_INTEGER,
     .                   nodes(ibloc1)-1,TAG_FLOW,
     .                   mycomm,istat,ierr)
#endif
             endif
         endif
#endif

c
# if defined BUILD_MPI
       if (nodes(ibloc2)-1.eq.myrank) then
# else
       if (nodes(ibloc2).eq.myrank) then
# endif
c
cvn (10-22-97)
# if defined BUILD_MPI
         if (nodes(ibloc)-1.eq.myrank) then
# else
         if (nodes(ibloc).eq.myrank) then
# endif
         call bcpchget (imn,jmn,kmn,
     .     im  (igrid,ibloc2), jm  (igrid,ibloc2), km  (igrid,ibloc2),
     .     imp2(igrid,ibloc2), jmp2(igrid,ibloc2), kmp2(igrid,ibloc2),
     .     vol(m1cc(igrid,ibloc2)),
     .     iface2,ipitmb2(litmbeg),jpitmb2(litmbeg),kpitmb2(litmbeg),
     .     litems,nghost,wk2da)
c
         else
c
         call bcpchget (imn,jmn,kmn,
     .     im  (igrid,ibloc2), jm  (igrid,ibloc2), km  (igrid,ibloc2),
     .     imp2(igrid,ibloc2), jmp2(igrid,ibloc2), kmp2(igrid,ibloc2),
     .     vol(m1cc(igrid,ibloc2)),
     .     iface2,ipitmbs(iitmbeg),jpitmbs(iitmbeg),kpitmbs(iitmbeg),
     .     litems,nghost,wk2da)
c
         endif

#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
# if defined BUILD_MPI
         if (nodes(ibloc)-1.ne.myrank) then
# else
         if (nodes(ibloc).ne.myrank) then
# endif
c (vatsa: 3-18-97)
           iitmbeg = iitmbeg+litems
           iitmsa(lpchs,igrid) = iitmbeg
c
          if( iitmbeg.gt.mpchitm) then
             write (iwrit,'(2x,"dimension conflict for mpchitm "/)')
             write (iwrit,'(2x,"mpchitm litmbeg litems igrid iblc2"/)')
             write (iwrit,'(2x,5i7)')mpchitm,iitmbeg,litems,igrid,ibloc
             write (iwrit,'(2x,"stop in  bcvol before sending itms"/)')
c
             call ERREXIT (nodes)

          endif
c
c
c
#ifdef BUILD_PVM
           call PVMFpsend (nodes(ibloc),TAG_VOL,
     .                     wk2da,ncnt,RTYPE,ierr)
#else
           call MPI_Send (wk2da,ncnt,RTYPE,
     .                    nodes(ibloc)-1,TAG_VOL,mycomm,ierr)
#endif
         endif
#endif
c
       endif
c
# if defined BUILD_MPI
       if (nodes(ibloc)-1.eq.myrank) then
# else
       if (nodes(ibloc).eq.myrank) then
# endif
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
# if defined BUILD_MPI
         if (nodes(ibloc2)-1.ne.myrank) then
# else
         if (nodes(ibloc2).ne.myrank) then
# endif
#ifdef BUILD_PVM
           call PVMFprecv (nodes(ibloc2),TAG_VOL,
     .                     wk2da,ncnt,RTYPE,
     .                     itid,itag,ilen,ierr)
#else
           call MPI_Recv (wk2da,ncnt,RTYPE,
     .                    nodes(ibloc2)-1,TAG_VOL,mycomm,istat,ierr)
#endif
         endif
#endif
c
         call bcpchset (imn,jmn,kmn,
     .     im  (igrid,ibloc),jm  (igrid,ibloc),km  (igrid,ibloc),
     .     imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .     vol(m1cc(igrid,ibloc)),
     .     iface1,ipitmb1(litmbeg),jpitmb1(litmbeg),kpitmb1(litmbeg),
     .     litems,nghost,frc(litmbeg),wk2da)
c      endif
c
       litmbeg = litmbeg +litems
       endif
 120  continue
c
 101  continue
 1000 continue
c
      do 1002 ibloc=1,nbloc
c---  fill in edges (corners) of block boundaries with extrapolation b.c
c
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) then
# else
      if (nodes(ibloc).eq.myrank) then
# endif
c         write (iwrit,'(2x,"finish bcvol for igrid ibloc "/)')
c         write (iwrit,'(2x,2i7)') igrid,ibloc
cvn       write (iwrit,'(2x,"finish bcvol for",
cvn  .                   2x,4i5)') igrid,ibloc,nodes(ibloc),myrank

          call bcedggr (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .    vol(m1cc(igrid,ibloc)),
     .    igrid, isoln              )
      endif
c
c----      end loop on blocks
 1002 continue
c
c
c-------  begin outer loop on the blocks for non-interface boundaries  ----
c
      do 2000 ibloc = 1,nbloc
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) then
# else
      if (nodes(ibloc).eq.myrank) then
# endif
c
        ns       =  nseg(ibloc)
c
c----------  begin outer loop on the segments  ----------------------------
c
        do 200 iseg = 1,ns
c
        nbctype  =  imap(1 ,iseg ,ibloc)
c
        if (nbctype .gt. 1) then
c
          nface    =  imap(2 ,iseg ,ibloc)
          n1beg    =  imap(3 ,iseg ,ibloc)
          n1end    =  imap(4 ,iseg ,ibloc)
          n2beg    =  imap(5 ,iseg ,ibloc)
          n2end    =  imap(6 ,iseg ,ibloc)
c
          if (nbctype .eq. 4) nsym = imap(8 ,iseg ,ibloc)
c
c--------------- extrapolation condition    ---------------------------
c---  same treatment is used at this time for cell-volumes
c---  for symmetry, walls and far-field  type boundary conditions
c---  ghost cell value is set equal to its interior neighbor
c
          call bcextgr (imn,jmn,kmn,
     .     im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .     imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .     vol(m1cc(igrid,ibloc )),
     .     nface ,n1beg ,n1end ,n2beg ,n2end , isoln   )
c
         endif
c
c----   end loop on segments
  200   continue
c
c----  end loop on blocks
       endif
 2000  continue
c
c
       return
       end
