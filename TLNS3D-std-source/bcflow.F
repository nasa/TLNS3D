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
      subroutine bcflow(imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2,m1np,m2np,m3np,m4np,m5np,
     .                  m1cc,m2cc,m3cc,m4cc,m5cc,m5cg,
     .                  m1is,m1js,m1ks,m3is,m3js,m3ks,
     .                  mgrlev,igrid,isoln,mbloc,nbloc,
     .                  mxsizc,mxsiz5c,mxsizn,mxsiz3n,mx1dwk,
     .                  imap,msegt,msegn,nseg,w,p,eomu,turv1,turv2,vol,
     .                  x,six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                  ri1,ri2,rj1,rj2,rk1,rk2,
     .                  tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .                  uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
     .                  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .                  wk2d,const1,const2,const3,const4,const5,
     .                  mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .                  mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .                  lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .               kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .           ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,iunsteady,totptim,
     .           smin )
c----------------------------------------------------------------------
c
c     object : To set ghost-cell values for flow quantities
c              based on topological mapping funcions, 'imap'
c
c     initial coding : by v.n.vatsa (March, 1990)
c     latest mods    : To accommodate patched inter-faces (Vatsa: Sept. 1994)
c                      nbctyp < 0 for patched boundaries
c
c     modified for distributed computing : April, 1995
c
c     modified May,95  to remove eomu and p transmission since
c     eomu is updated for turbulent flows in bcturb, and p can be computed
c     from w
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
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      common/facetr/ ifacetr(6),ipatchg
c
      dimension     im  (mgrlev,mbloc), jm  (mgrlev,mbloc),
     .              km  (mgrlev,mbloc), imp1(mgrlev,mbloc),
     .              jmp1(mgrlev,mbloc), kmp1(mgrlev,mbloc),
     .              imp2(mgrlev,mbloc), jmp2(mgrlev,mbloc),
     .              kmp2(mgrlev,mbloc),
     .              imap(msegt,msegn,mbloc), nseg(mbloc)
c
      dimension     m1np(mgrlev,mbloc), m2np(mgrlev,mbloc),
     .              m3np(mgrlev,mbloc), m4np(mgrlev,mbloc),
     .              m5np(mgrlev,mbloc),
     .              m1cc(mgrlev,mbloc), m2cc(mgrlev,mbloc),
     .              m3cc(mgrlev,mbloc), m4cc(mgrlev,mbloc),
     .              m5cc(mgrlev,mbloc), m5cg(mgrlev,mbloc),
     .              m1is(mgrlev,mbloc), m1js(mgrlev,mbloc),
     .              m1ks(mgrlev,mbloc),
     .              m3is(mgrlev,mbloc), m3js(mgrlev,mbloc),
     .              m3ks(mgrlev,mbloc)
c
      dimension     w (mxsiz5c), p(mxsizc ), smin(mxsizc ),
     .              vol(mxsizc), x(mxsiz3n), eomu(mxsizc ),
     .                                      turv1(mxsizc ),
     .                                      turv2(mxsizc ),
     .              six(mxsizn),siy(mxsizn),siz(mxsizn),
     .              sjx(mxsizn),sjy(mxsizn),sjz(mxsizn),
     .              skx(mxsizn),sky(mxsizn),skz(mxsizn)
c
      dimension     ivisb(3,mbloc),      const1(msegn,mbloc),
     .              const2(msegn,mbloc),const3(msegn,mbloc),
     .              const4(msegn,mbloc),const5(msegn,mbloc)
c
      dimension     ri1(1),ri2(1),rj1(1),rj2(1),rk1(1),rk2(1)
c
      dimension     tauwfi1(1),tauwfi2(1),tauwfj1(1),tauwfj2(1),
     .              tauwfk1(1),tauwfk2(1)
      dimension     uplwfi1(1),uplwfi2(1),uplwfj1(1),uplwfj2(1),
     .              uplwfk1(1),uplwfk2(1)
      dimension     uswfi1(1),uswfi2(1),uswfj1(1),uswfj2(1),
     .              uswfk1(1),uswfk2(1)
c
c     array declaration for temporary work-space
c     this needs to be at least size (mx1dwk,mx1dwk,14)
c     corresponding to 5 variables for w, p and eomu times
c     2 ghost cells on finest grid.

      dimension   wk2d(mx1dwk*mx1dwk*15)
c
c     patcher related arrays
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
c     integer * 8 loc, locw

C ***********
      save jsegnummxs,   jsegnummxr,   jdanaszmxs,   jdanaszmxr
      data jsegnummxs/0/,jsegnummxr/0/,jdanaszmxs/0/,jdanaszmxr/0/
C ***********

#include "async.dc"
#ifdef ASYNC
c dana recv buffer
      integer idanareq(ibufct)
      integer danasa(MPI_STATUS_SIZE,ibufct)
      dimension danarecv(ibufct*ibufsz)
      integer isegct,isegnum
      common /async/ idanareq,danasa,danarecv,isegct,isegnum
#include "bcflow.Fasync"
#else
#include "bcflow.Fsync"
#endif
c
c-------  begin outer loop on the blocks for non-interface boundaries  ----
c
      do 2000 ibloc = 1,nbloc
# ifdef BUILD_MPI
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
        if (nbctype.gt.1) then
c
          nface    =  imap(2 ,iseg ,ibloc)
          n1beg    =  imap(3 ,iseg ,ibloc)
          n1end    =  imap(4 ,iseg ,ibloc)
          n2beg    =  imap(5 ,iseg ,ibloc)
          n2end    =  imap(6 ,iseg ,ibloc)
          extyp    =  imap(7 ,iseg ,ibloc)
c
          if (nbctype .eq. 4) nsym  = imap(8 ,iseg ,ibloc)
c
c---------------  slip (inviscid wall )  condition     -------------------
        if(nbctype .eq. 2) 
     .    call bcinwfl (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .     six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .     siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .     sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .     skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .     skz(m1np(igrid,ibloc)),
     .     ri1(m3is(igrid,ibloc)), ri2(m3is(igrid,ibloc)),
     .     rj1(m3js(igrid,ibloc)), rj2(m3js(igrid,ibloc)),
     .     rk1(m3ks(igrid,ibloc)), rk2(m3ks(igrid,ibloc)),
     .     mx1dwk ,wk2d,
     .     nface ,n1beg ,n1end ,n2beg ,n2end , igrid, isoln,  extyp )
c
c---------------  slip (inviscid wall )  with wall-function     --------------
        if(nbctype .eq. 102) 
     .    call bcinwfun (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    smin(m1cc(igrid,ibloc)),    x(m3np(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .     six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .     siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .     sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .     skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .     skz(m1np(igrid,ibloc)), vol(m1cc(igrid,ibloc)),
     .     ri1(m3is(igrid,ibloc)), ri2(m3is(igrid,ibloc)),
     .     rj1(m3js(igrid,ibloc)), rj2(m3js(igrid,ibloc)),
     .     rk1(m3ks(igrid,ibloc)), rk2(m3ks(igrid,ibloc)),
     .     tauwfi1(m1is(igrid,ibloc)), tauwfi2(m1is(igrid,ibloc)),
     .     tauwfj1(m1js(igrid,ibloc)), tauwfj2(m1js(igrid,ibloc)),
     .     tauwfk1(m1ks(igrid,ibloc)), tauwfk2(m1ks(igrid,ibloc)),
     .     uswfi1 (m1is(igrid,ibloc)), uswfi2 (m1is(igrid,ibloc)),
     .     uswfj1 (m1js(igrid,ibloc)), uswfj2 (m1js(igrid,ibloc)),
     .     uswfk1 (m1ks(igrid,ibloc)), uswfk2 (m1ks(igrid,ibloc)),
     .     uplwfi1(m1is(igrid,ibloc)), uplwfi2(m1is(igrid,ibloc)),
     .     uplwfj1(m1js(igrid,ibloc)), uplwfj2(m1js(igrid,ibloc)),
     .     uplwfk1(m1ks(igrid,ibloc)), uplwfk2(m1ks(igrid,ibloc)),
     .     mx1dwk ,wk2d,
     .     nface ,n1beg ,n1end ,n2beg ,n2end , igrid, isoln,
     .     extyp, const1(iseg,ibloc)  )
c
c---------------  no-slip (viscous wall )  condition     ----------------
        if(nbctype .eq. 3) then
cvn      if(const5(iseg,ibloc).le.0.) then
         if(const5(iseg,ibloc).le.1000.) then
          call bcvwfl (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .                            turv2(m1cc(igrid,ibloc)),
     .     six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .     siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .     sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .     skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .     skz(m1np(igrid,ibloc)),
     .     ri1(m3is(igrid,ibloc)), ri2(m3is(igrid,ibloc)),
     .     rj1(m3js(igrid,ibloc)), rj2(m3js(igrid,ibloc)),
     .     rk1(m3ks(igrid,ibloc)), rk2(m3ks(igrid,ibloc)),
     .     mx1dwk ,wk2d,
     .     nface ,n1beg ,n1end ,n2beg ,n2end , igrid, isoln,  extyp,
     .     const1(iseg,ibloc),const2 (iseg,ibloc),
     .     const3(iseg,ibloc),const4(iseg,ibloc),const5(iseg,ibloc),
     .     totptim, iunsteady, nodes,ibloc,myrank          )
         else
c         radiation b.c. (note: this b.c.is implemented for jmn only)
          call bcvwfl2 (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .     six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .     siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .     sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .     skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .     skz(m1np(igrid,ibloc)),
     .     ri1(m3is(igrid,ibloc)), ri2(m3is(igrid,ibloc)),
     .     rj1(m3js(igrid,ibloc)), rj2(m3js(igrid,ibloc)),
     .     rk1(m3ks(igrid,ibloc)), rk2(m3ks(igrid,ibloc)),
     .     mx1dwk ,wk2d,
     .     nface ,n1beg ,n1end ,n2beg ,n2end , igrid, isoln,  extyp,
     .     const1(iseg,ibloc),const2 (iseg,ibloc),
     .     const3(iseg,ibloc),const4(iseg,ibloc),
     .     const5(iseg,ibloc),
     .     wk2d(1+mx1dwk**2),x(m3np(igrid,ibloc)), nodes )
          endif
c
        endif
c
c
c---------------  symmetry condition        ----------------------------
        if(nbctype .eq. 4) 
     .    call bcsymfl (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .     ri1(m3is(igrid,ibloc)), ri2(m3is(igrid,ibloc)),
     .     rj1(m3js(igrid,ibloc)), rj2(m3js(igrid,ibloc)),
     .     rk1(m3ks(igrid,ibloc)), rk2(m3ks(igrid,ibloc)),
     .    nface,n1beg,n1end,n2beg,n2end,igrid,isoln,nsym,extyp       )
c
c------------  downstream/extrapolation condition      -----------------
        if(nbctype .eq. 5) 
     .    call bcextfl (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .     six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .     siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .     sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .     skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .     skz(m1np(igrid,ibloc)),
     .       x(m3np(igrid,ibloc)),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .    nface ,n1beg ,n1end ,n2beg ,n2end, igrid, isoln,  extyp,
     .    const1(iseg,ibloc) )
c
c-----------  far-field (1-d characteristics)  condition  ---------------
        if(nbctype .eq. 6) then
c         locw = loc (w(m5cc(igrid,ibloc) + 3405))
          call bcfarfl (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .     ri1(m3is(igrid,ibloc)), ri2(m3is(igrid,ibloc)),
     .     rj1(m3js(igrid,ibloc)), rj2(m3js(igrid,ibloc)),
     .     rk1(m3ks(igrid,ibloc)), rk2(m3ks(igrid,ibloc)),
     .     nface ,n1beg ,n1end ,n2beg ,n2end, igrid, isoln,  extyp ,
     .     iunsteady )
        endif
c
c------------   boundary condition for singular lines      -----------------
c       if(nbctype .eq. 7) 
c    .    call bcsinfl (imn,jmn,kmn,
        if(nbctype .eq. 7) 
     .    call bcextfl (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .     six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .     siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .     sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .     skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .     skz(m1np(igrid,ibloc)),
     .       x(m3np(igrid,ibloc)),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .    nface ,n1beg ,n1end ,n2beg ,n2end, igrid, isoln,  extyp,
     .    const1(iseg,ibloc) )
c
c----------  fan-upstream (str. tube contraction ratio) condition      -------
        if(nbctype .eq. 8) 
     .    call bcfnufl (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .    nface ,n1beg ,n1end ,n2beg ,n2end, igrid, isoln,  extyp,
     .     const1(iseg,ibloc) )
c
c---------------  engine/fan (mass-flow) condition      -------------------
        if(nbctype .eq. 9) 
     .    call bcengmfl (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .     six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .     siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .     sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .     skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .     skz(m1np(igrid,ibloc)),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .       x(m3np(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .    nface ,n1beg ,n1end ,n2beg ,n2end, igrid, isoln,  extyp,
     .     const1(iseg,ibloc),const2 (iseg,ibloc),nodes )
c
c------------  jet-nozzle inflow/outflow condition      -----------------
        if(nbctype .eq. 10) 
     .    call bcjetfl (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .     ri1(m3is(igrid,ibloc)), ri2(m3is(igrid,ibloc)),
     .     rj1(m3js(igrid,ibloc)), rj2(m3js(igrid,ibloc)),
     .     rk1(m3ks(igrid,ibloc)), rk2(m3ks(igrid,ibloc)),
     .     nface ,n1beg ,n1end ,n2beg ,n2end, igrid, isoln,  extyp,
     .     const1(iseg,ibloc),const2(iseg,ibloc),const3(iseg,ibloc),
     .     const4(iseg,ibloc),const5(iseg,ibloc) )
c
c
c-----  flow-streams with total pressure, temp. &  inflow angles specified -----
c       
        if(nbctype .eq. 11) 
     .    call bcpttfla (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .       x(m3np(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .     nface ,n1beg ,n1end ,n2beg ,n2end, igrid, isoln,  extyp,
     .     const1(iseg,ibloc),const2(iseg,ibloc),const3(iseg,ibloc),
     .     const4(iseg,ibloc),const5(iseg,ibloc))
c
c
c----------  flow-streams with total pressue & temp. specification  ----
c           ( inflow stream is aligned with grid-lines )
        if(nbctype .eq. 12)
     .    call bcpttflb (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .       x(m3np(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .     ri1(m3is(igrid,ibloc)), ri2(m3is(igrid,ibloc)),
     .     rj1(m3js(igrid,ibloc)), rj2(m3js(igrid,ibloc)),
     .     rk1(m3ks(igrid,ibloc)), rk2(m3ks(igrid,ibloc)),
     .     nface ,n1beg ,n1end ,n2beg ,n2end, igrid, isoln,  extyp,
     .     const1(iseg,ibloc),const2(iseg,ibloc),const3(iseg,ibloc))
c
c
c-----------  Moving Diaphragm condition   --------------------
        if(nbctype .eq. 13)
     .    call bcdiaph (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .     six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .     siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .     sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .     skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .     skz(m1np(igrid,ibloc)),
     .     ri1(m3is(igrid,ibloc)), ri2(m3is(igrid,ibloc)),
     .     rj1(m3js(igrid,ibloc)), rj2(m3js(igrid,ibloc)),
     .     rk1(m3ks(igrid,ibloc)), rk2(m3ks(igrid,ibloc)),
     .     mx1dwk ,wk2d,
     .     nface ,n1beg ,n1end ,n2beg ,n2end , igrid, isoln,  extyp,
     .     const1(iseg,ibloc),const4(iseg,ibloc),const5(iseg,ibloc),
     .     totptim,x(m3np(igrid,ibloc))        )
c
c
      end if
c
c-----      end loop on segments
  200  continue
c
c-----      end loop on blocks
      endif
 2000 continue
      do 2002 ibloc=1,nbloc
c
c---      fill in edges (corners) of block boundaries with extrapolation b.c
c

c     if (nodes(ibloc)-1.eq.myrank) then

c         call bcedgfl (imn,jmn,kmn,
c    .    im  (igrid,ibloc),jm  (igrid,ibloc),km  (igrid,ibloc),
c    .    imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
c    .    imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
c    .       w(m5cc(igrid,ibloc)),    p(m1cc(igrid,ibloc)),
c    .    eomu(m1cc(igrid,ibloc)),
c    .    igrid, isoln              )
c
c     endif
c-----      end loop on blocks
 2002 continue

c
       return
       end
