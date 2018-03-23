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
      subroutine bcturb(imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2,m1cc,m5cc,
     .                  mgrlev,igrid,isoln,mbloc,nbloc,
     .                  mxsizc,mxsiz5c,mx1dwk,imap,msegt,msegn,nseg,
     .                  wk2d,eomu,turv1,turv2,w,p,smin,
     .                  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,
     .                  uplwfk2,tnuwfi1,tnuwfi2,tnuwfj1,tnuwfj2,tnuwfk1,
     .                  tnuwfk2,m1is,m1js,m1ks,
     .              mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .              mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .              lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .              kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .              ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes )
c----------------------------------------------------------------------
c
c     object : To set ghost-cell values for eddy-viscosity
c              based on topological mapping funcions, 'imap'
c
c     initial coding : by v.n.vatsa (May, 1991)
c     latest mods    : To accommodate patched inter-faces (Vatsa: Nov. 1994)
c                      nbctyp < 0 for patched boundaries
c
c     modified for distributed computing : April, 1995
c
#undef ASYNC
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
      common/trbfs/ anutinf,zkinf,ominf
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c
      dimension     im  (mgrlev,mbloc), jm  (mgrlev,mbloc),
     .              km  (mgrlev,mbloc), imp1(mgrlev,mbloc),
     .              jmp1(mgrlev,mbloc), kmp1(mgrlev,mbloc),
     .              imp2(mgrlev,mbloc), jmp2(mgrlev,mbloc),
     .              kmp2(mgrlev,mbloc),
     .              imap(msegt,msegn,mbloc), nseg(mbloc)
c
      dimension     m1cc(mgrlev,mbloc), m5cc(mgrlev,mbloc)
c
      dimension     w(mxsiz5c),   p(mxsizc) , smin(mxsizc)
c
      dimension     eomu(mxsizc), turv1(mxsizc), turv2(mxsizc)
c
c     wall function related arrays
c
      dimension     m1is(mgrlev,mbloc), m1js(mgrlev,mbloc),
     .              m1ks(mgrlev,mbloc)

      dimension     uplwfi1(1),uplwfi2(1),uplwfj1(1),uplwfj2(1),
     .              uplwfk1(1),uplwfk2(1)
c
      dimension     tnuwfi1(1),tnuwfi2(1),tnuwfj1(1),tnuwfj2(1),
     .              tnuwfk1(1),tnuwfk2(1)
c
c     patcher related arrays
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
      dimension   wk2d(mx1dwk*mx1dwk*15)
c
#include "async.dc"
#ifdef ASYNC
c dana recv buffer
      integer idanareq(ibufct)
      integer danasa(MPI_STATUS_SIZE,ibufct)
      dimension danarecv(ibufct*ibufsz)
      integer isegct,isegnum
      common /async/ idanareq,danasa,danarecv,isegct,isegnum
#include "bcturb.Fasync"
#else
#include "bcturb.Fsync"
#endif
c
c-------  begin outer loop on the blocks for non-interface boundaries  ----
c
cvn       write (iwrit,'(2x,"enter bcturb 2000",
cvn  .                   2x,4i5)') igrid,myrank
      do 2000 ibloc = 1,nbloc
# ifdef BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) then
# else
      if (nodes(ibloc).eq.myrank) then
# endif
cvn       write (iwrit,'(2x,"enter bcturb 2000",
cvn  .                   2x,4i5)') igrid,ibloc,nodes(ibloc),myrank
c
       ns       =  nseg(ibloc)
c
c----------  begin outer loop on the segments  ----------------------------
c
       do 200 iseg = 1,ns
c
       nbctype  =  imap(1 ,iseg ,ibloc)
c
       if(nbctype .gt. 1) then
c
           nface    =  imap(2 ,iseg ,ibloc)
           n1beg    =  imap(3 ,iseg ,ibloc)
           n1end    =  imap(4 ,iseg ,ibloc)
           n2beg    =  imap(5 ,iseg ,ibloc)
           n2end    =  imap(6 ,iseg ,ibloc)
c
c---------------  extrapolation condition        -----------------------
c---  same treatment is used at this time for eddy-viscosity
c---  for symmetry, and walls type boundary conditions
c---  ghost cell values are set equal to its interior neighbor
c
cvn       if(nbctype.ne.6.and.nbctype.ne.3.and.nbctype.ne.10)
          if(nbctype.ne.6.and.nbctype.ne.3.and.nbctype.ne.8.
     .    and.nbctype.ne.9.and.nbctype.ne.10.and.nbctype.ne.102)
     .    call bcextev (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .    eomu(m1cc(igrid,ibloc )), 
     .    turv1(m1cc(igrid,ibloc )), turv2(m1cc(igrid,ibloc)),
     .    nface ,n1beg ,n1end ,n2beg ,n2end , isoln            )
c
          if(nbctype.eq.3)
     .    call bcvwev  (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc )),   smin(m1cc(igrid,ibloc)),
     .    turv1(m1cc(igrid,ibloc )), turv2(m1cc(igrid,ibloc)),
     .    nface ,n1beg ,n1end ,n2beg ,n2end ,isoln,mx1dwk,wk2d )
c
          if(nbctype.eq.102)
     .    call bcvwevwf  (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc )),   smin(m1cc(igrid,ibloc)),
     .    turv1(m1cc(igrid,ibloc )), turv2(m1cc(igrid,ibloc)),
     .    uplwfi1 (m1is(igrid,ibloc)), uplwfi2 (m1is(igrid,ibloc)),
     .    uplwfj1 (m1js(igrid,ibloc)), uplwfj2 (m1js(igrid,ibloc)),
     .    uplwfk1 (m1ks(igrid,ibloc)), uplwfk2 (m1ks(igrid,ibloc)),
     .    tnuwfi1 (m1is(igrid,ibloc)), tnuwfi2 (m1is(igrid,ibloc)),
     .    tnuwfj1 (m1js(igrid,ibloc)), tnuwfj2 (m1js(igrid,ibloc)),
     .    tnuwfk1 (m1ks(igrid,ibloc)), tnuwfk2 (m1ks(igrid,ibloc)),
     .    nface ,n1beg ,n1end ,n2beg ,n2end ,isoln,mx1dwk,wk2d )
c
cvn       if(nbctype.eq.6.or.nbctype.eq.10)
          if(nbctype.eq.6.or.nbctype.ge.8)
     .    call bcfarev (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .    eomu(m1cc(igrid,ibloc )), 
     .    turv1(m1cc(igrid,ibloc )), turv2(m1cc(igrid,ibloc)),
     .    nface ,n1beg ,n1end ,n2beg ,n2end , isoln            )
c
       endif
c
c----      end loop on segments
  200  continue
c
c----      end loop on blocks
       endif
 2000  continue
c
c---      fill in edges (corners) of block boundaries with extrapolation b.c
c
c     do 2002 ibloc=1,nbloc

c     if (nodes(ibloc)-1.eq.myrank) then


c         call bcedgev (imn,jmn,kmn,
c    .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
c    .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
c    .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
c    .    eomu(m1cc(igrid,ibloc)),
c    .    turv1(m1cc(igrid,ibloc)), turv2(m1cc(igrid,ibloc)),
c    .    igrid, isoln              )
c
c     endif
c
c----      end loop on blocks
 2002 continue
c
c
      return
      end
