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
      subroutine turbdr(imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2,m1np,m2np,m3np,m4np,m5np,
     .                  m1cc,m2cc,m3cc,m4cc,m5cc,m5cg,
     .                  m1is,m1js,m1ks,m3is,m3js,m3ks,
     .                  mgrlev,igrid,isoln,mbloc,nbloc,mxsizc,mxsiz4c,
     .                  mxsiz5c,mxsizn,mxsiz3n,mx1dwk,mx3dwk,
     .                  imap,msegt,msegn,nseg,w,p,eomu,turv1,turv2,
     .                  vol,x,six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                  ri1,ri2,rj1,rj2,rk1,rk2,
     .                  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .                  tnuwfi1,tnuwfi2,tnuwfj1,tnuwfj2,tnuwfk1,tnuwfk2,
     .                  wk1d, wk2d, wk3d, fv, dtl, dti, dtj, dtk,
     .                  dtvi,dtvj,dtvk,mxszfg,m1fg,ivisb,mstrt,smin,
     .                  iturbb,itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,
     .             mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .             mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .             lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .             kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .             ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,mtorder,
     .             ntorder,iturv1,iturv2,m1fgt,tv1old,tv2old,dtphy,
     .             mebdf,ibdf2opt, tv1bar, tv2bar,
     .             m1fgtrk,m5fgtrk,irkstg,nrkstg,mrkstg,
     .             rkphysa,rkphyse,dtv1ork,dtv2ork,itime )
c----------------------------------------------------------------------
c
c     object : To serve as a driver routine for evaluating
c              eddy-viscosity using a selected turbulence
c              model over a selected range on the fine mesh
c              based on topological mapping funcions, 'imap'
c
c     initial coding : by v.n.vatsa (April, 1991)
c
c     modified for distributed computing : April, 1995
c

      character*128  errmsg
      common/dstcom/ errmsg,master,myrank,mycomm,ndlist,nnodes
      dimension      nodes(1)
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension     im  (mgrlev,mbloc), jm  (mgrlev,mbloc),
     .              km  (mgrlev,mbloc), imp1(mgrlev,mbloc),
     .              jmp1(mgrlev,mbloc), kmp1(mgrlev,mbloc),
     .              imp2(mgrlev,mbloc), jmp2(mgrlev,mbloc),
     .              kmp2(mgrlev,mbloc),
     .              imap(msegt,msegn,mbloc), nseg(mbloc)
c
      dimension     ivisb(3,mbloc)
c
      dimension     itrb1(mgrlev,mbloc), itrb2(mgrlev,mbloc),
     .              jtrb1(mgrlev,mbloc), jtrb2(mgrlev,mbloc),
     .              ktrb1(mgrlev,mbloc), ktrb2(mgrlev,mbloc), 
     .              iturbb(mbloc)
c
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
      dimension     w (mxsiz5c), p(mxsizc ),
     .              vol(mxsizc), x(mxsiz3n), eomu(mxsizc),
     .              smin(mxsizc), turv1(mxsizc),turv2(mxsizc),
     .              six(mxsizn),siy(mxsizn),siz(mxsizn),
     .              sjx(mxsizn),sjy(mxsizn),sjz(mxsizn),
     .              skx(mxsizn),sky(mxsizn),skz(mxsizn)
c
c     array allocation for time dependent quantities
c
      dimension     tv1old(mxszfg*ntorder*iturv1+1),
     .              tv2old(mxszfg*ntorder*iturv2+1)
      dimension     tv1bar(mxszfg*(ntorder-1)*iturv1+1),
     .              tv2bar(mxszfg*(ntorder-1)*iturv2+1)
c
      dimension     m1fgt(mtorder,mbloc)
c
      dimension     m1fg(mbloc)
c
c#ifdef  RKPHY
      dimension     m1fgtrk(nrkstg ,mbloc), m5fgtrk(nrkstg ,mbloc)
c
      dimension     rkphysa(mrkstg,mrkstg),rkphyse(mrkstg)
c
      dimension     dtv1ork (mxszfg*nrkstg*iturv1+1),
     .              dtv2ork (mxszfg*nrkstg*iturv2+1)
c#endif
c
      dimension     dtvi(mxszfg), dtvj(mxszfg), dtvk(mxszfg)
      dimension     dti(mxszfg), dtj(mxszfg), dtk(mxszfg)
     .             ,dtl(1)
c
      dimension     fv(1)
c
      dimension     ri1(1),ri2(1),rj1(1),rj2(1),rk1(1),rk2(1)
c
      dimension     uplwfi1(1),uplwfi2(1),uplwfj1(1),uplwfj2(1),
     .              uplwfk1(1),uplwfk2(1)
      dimension     tnuwfi1(1),tnuwfi2(1),tnuwfj1(1),tnuwfj2(1),
     .              tnuwfk1(1),tnuwfk2(1)
c
c     array declaration for temporary work-space
c
c     1-d arrays
c
      dimension   wk1d(mx1dwk,27)
c
c     2-d arrays
c
      dimension   wk2d(mx1dwk*mx1dwk,15)
c
c     3-d arrays
c
      dimension   wk3d(mx3dwk,5)
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
      data init/0/
c
c     In this subroutine and in subsequent calls to vort and turbulent
c     routines the 3-d work arrays represent the following variables
c
c     wk3d(1) = up
c     wk3d(2) = vp
c     wk3d(3) = wp
c     wk3d(4) = ts
c     wk3d(5) = amu
c
c***********************************************************************
c--------------------  mapping function description ---------------------
c
c      imap  : mapping function containing topological information
c      msegt : maximum no. of types of operations/boundary conditons
c              (currently set to 20)
c      msegn : maximum no. of segments permitted on all faces of a block
c              (currently set to 20)
c
c      imap(1 , ) : specifies boundary/operation type
c      imap(2 , ) : specifies face number (1-6)
c      imap(3 , ) : specifies beginning of direction 1
c      imap(4 , ) : specifies end       of direction 1
c      imap(5 , ) : specifies beginning of direction 2
c      imap(6 , ) : specifies end       of direction 2
c
c      imap(8,  ) : if the boundary type is symmetry
c                   takes the value 1,2 or 3 to indicate symmetry 
c                   about x, y or z = constant planes, respectively
c 
c  if the boundary type is a cut/interface
c  additional information described below is required 
c
c      imap(7 , ) : specifies block no. of source segment
c     |imap(8 , )|: specifies face no. of source segment
c      imap(8 , ) > 0 corresponds to direction 1 of source segment
c                      matching with direction 1 of target segment
c                      and same  for direction 2
c      imap(8 , ) < 0 corresponds to direction 2 of source segment
c                      matching with direction 1 of target segment
c                      and vice-versa
c
c      imap(9 , ) : specifies beginning of direction 1 of source segment
c      imap(10, ) : specifies end       of direction 1 of source segment
c      imap(11, ) : specifies beginning of direction 2 of source segment
c      imap(12, ) : specifies end       of direction 2 of source segment
c
c      turbulence/transitional related information
c    
c      imap(13, ) : specifies if there is turbulent flow on this segment
c                   1  corresponds to turbulent flow
c                   0  corresponds to laminar   flow
c      imap(14, ) : begining index in direction 1 for turbulent flow
c      imap(15, ) : ending   index in direction 1 for turbulent flow
c      imap(16, ) : begining index in direction 2 for turbulent flow
c      imap(17, ) : ending   index in direction 2 for turbulent flow
c      imap(18, ) : begining index in normal direction for fmax 
c      imap(19, ) : ending   index in normal direction for fmax 
c      imap(20, ) : ending   index in normal direction for turb. flow 
c
c--------------------  boundary/operation type description ---------------
c
c      nbctype    = imap(1 , )
c                 = 0  corresponds to an inner cut that maps to same block
c                 = 1  corresponds to a cut that maps to another block
c                 = 2  corresponds to a slip (inviscid) wall
c                 = 3  corresponds to a noslip (viscous) wall
c                 = 4  symmetry condition (imap(8) tells about which plane)
c                 = 5  downstream/extrapolation condition
c                 = 6  far-field condition (Riemann invariants)
c
c***************************************************************************
c
c----   begin subiteration on field-equations ----
c
      nsubit = 2
c     nsubit = 8
      nsubit = 1
      if(init.le.5.and.mstrt.eq.0) nsubit = 5
      if(iturb.lt.2) nsubit = 1
c
c************************************************************************
c     mods to correct subiteration in field equations (vatsa, July 95)
      nsubitr = nsubit
c     (pass nsubitr instead of isubit to turbsa and turbsst as argument)
c************************************************************************
c
c     do 2000 isubit = 1,nsubit
c
c----------  begin outer loop on the blocks  ----------------------------
c
       do 1000 ibloc = 1,nbloc
# if defined BUILD_MPI
       if (nodes(ibloc)-1.eq.myrank) then
# else
       if (nodes(ibloc).eq.myrank) then
# endif
c
       n1cc = m1cc(igrid,ibloc)
       n5cc = m5cc(igrid,ibloc)
       n1np = m1np(igrid,ibloc)
       n3np = m3np(igrid,ibloc)
c
c      compute the vorticity distribution on the entire block
c      unless it is flagged to be a laminar block
c      in that case fill-in some required arrays 
c      and then return
c
       call  vort (imn,jmn,kmn,
     . im  (igrid,ibloc),jm  (igrid,ibloc),km  (igrid,ibloc),
     . imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     . imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     . w(n5cc),p(n1cc),eomu(n1cc),vol(n1cc),x(n3np),
     . six(n1np),siy(n1np),siz(n1np),sjx(n1np),sjy(n1np),
     . sjz(n1np),skx(n1np),sky(n1np),skz(n1np), 
     . wk1d(1,1), wk1d(1,2), wk1d(1,3), wk1d(1,4), wk1d(1,5),
     . wk1d(1,6), wk1d(1,7), wk1d(1,8), wk1d(1,9), wk1d(1,10),
     . wk1d(1,11),wk1d(1,12),wk1d(1,13),wk1d(1,14),wk1d(1,15),
     . wk1d(1,16),wk1d(1,17),wk1d(1,18),wk1d(1,19),wk1d(1,20),
     . wk1d(1,21),wk1d(1,22),wk1d(1,23),wk1d(1,24),wk1d(1,25),
     . wk1d(1,26),wk1d(1,27),wk2d(1,1),wk2d(1,2),wk2d(1,3),
     . wk2d(1,4),wk2d(1,5),wk2d(1,6),wk2d(1,7),wk2d(1,8),wk2d(1,9),
     . wk3d(1,1),wk3d(1,2),wk3d(1,3),wk3d(1,4),wk3d(1,5),fv,
     . ivisb(1,ibloc) ,iturbb(ibloc)                           )
c
       if(iturbb(ibloc).ne.1) go to 1000
       if(iturb.le.0) go to 1000
       ns       =  nseg(ibloc)
c
c----------  begin outer loop on the segments  ----------------------------
c
c      bypass the call to turbulence routine on blocks where
c      viscous fluxes are not evaluated
c
       if(ivisb(1,ibloc).eq.0.and.ivisb(2,ibloc).eq.0.and.
     .    ivisb(3,ibloc).eq.0)  go to 1000
c
c      iss    is a counter on segments with viscous surfaces
       iss      = 0
       do 100 iseg = 1,ns
c
       nbctype  =  imap(1 ,iseg ,ibloc)
       nface    =  imap(2 ,iseg ,ibloc)
       n1beg    =  imap(3 ,iseg ,ibloc)
       n1end    =  imap(4 ,iseg ,ibloc)
       n2beg    =  imap(5 ,iseg ,ibloc)
       n2end    =  imap(6 ,iseg ,ibloc)
c
       nturb    =  imap(13,iseg ,ibloc)
       n1begt   =  imap(14,iseg ,ibloc)
       n1endt   =  imap(15,iseg ,ibloc)
       n2begt   =  imap(16,iseg ,ibloc)
       n2endt   =  imap(17,iseg ,ibloc)
c
       nnstrt   =  imap(18,iseg ,ibloc)
       nnstop   =  imap(19,iseg ,ibloc)
       nnendt   =  imap(20,iseg ,ibloc)
c
c
c
       if(iturb.eq.1) 
     .    call turbmb (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .    w(n5cc), p(n1cc), x(n3np), eomu(n1cc),
     .    ri1(m3is(igrid,ibloc)), ri2(m3is(igrid,ibloc)),
     .    rj1(m3js(igrid,ibloc)), rj2(m3js(igrid,ibloc)),
     .    rk1(m3ks(igrid,ibloc)), rk2(m3ks(igrid,ibloc)),
     .    wk1d(1,1),wk1d(1,2),wk1d(1,3),wk1d(1,4),wk1d(1,5),
     .    wk3d(1,1),wk3d(1,2),wk3d(1,3),wk3d(1,4),wk3d(1,5),
     .    mx1dwk, nface ,n1beg ,n1end ,n2beg , n2end , 
     .    nbctype,nturb ,n1begt,n1endt,n2begt, n2endt ,
     .    nnstrt ,nnstop,nnendt, fv  )
c
c
c----      end loop on segments
  100  continue
c
       if ( iturb.eq.2 )
     .  call turbsa(imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .    w(n5cc),p(n1cc),eomu(n1cc),turv1(n1cc),vol(n1cc),x(n3np),
     .    six(n1np),siy(n1np),siz(n1np),sjx(n1np),sjy(n1np),
     .    sjz(n1np),skx(n1np),sky(n1np),skz(n1np), 
     .    wk2d(1,1), wk2d(1,2), wk2d(1,3), wk2d(1,4), wk2d(1,5),
     .    wk2d(1,6), wk2d(1,7), wk2d(1,8), wk2d(1,9), wk2d(1,10),
     .    wk2d(1,11),wk2d(1,12),wk2d(1,13),wk2d(1,14),wk2d(1,15), fv,
     .    dtl,dti(m1fg(ibloc)),dtj(m1fg(ibloc)),dtk(m1fg(ibloc)),
     .    dtvi(m1fg(ibloc)),dtvj(m1fg(ibloc)),dtvk(m1fg(ibloc)),
     .    mx1dwk,wk3d(1,1),wk3d(1,2),wk3d(1,3),wk3d(1,4),wk3d(1,5),
     .    wk1d(1,1),wk1d(1,2),mstrt,smin(n1cc),ibloc,nsubitr,
     .    itrb1(igrid,ibloc),itrb2(igrid,ibloc),jtrb1(igrid,ibloc),
     .    jtrb2(igrid,ibloc),ktrb1(igrid,ibloc),ktrb2(igrid,ibloc),
     .    ntorder,tv1old(m1fgt(1,ibloc)),dtphy ,
     .    mebdf,ibdf2opt,tv1bar(m1fgtrk(1,ibloc)),
     .    irkstg,nrkstg,mrkstg,
     .    rkphysa,rkphyse,dtv1ork(m1fgtrk(1,ibloc)),itime )
c
       if ( iturb.eq.3 )
     .  call turbsst(imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .    w(n5cc),p(n1cc),eomu(n1cc),turv1(n1cc),turv2(n1cc), 
     .    vol(n1cc),x(n3np),six(n1np),siy(n1np),siz(n1np),sjx(n1np),
     .    sjy(n1np),sjz(n1np),skx(n1np),sky(n1np),skz(n1np), 
     .    wk2d(1,1), wk2d(1,2), wk2d(1,3), wk2d(1,4), wk2d(1,5),
     .    wk2d(1,6), wk2d(1,7), wk2d(1,8), wk2d(1,9), wk2d(1,10),
     .    wk2d(1,11),wk2d(1,12),wk2d(1,13),wk2d(1,14),wk2d(1,15), fv,
     .    dtl,dti(m1fg(ibloc)),dtj(m1fg(ibloc)),dtk(m1fg(ibloc)),
     .    dtvi(m1fg(ibloc)),dtvj(m1fg(ibloc)),dtvk(m1fg(ibloc)),
     .    mx1dwk,wk3d(1,1),wk3d(1,2),wk3d(1,3),wk3d(1,4),
     .    wk1d(1,1),wk1d(1,2),mstrt,smin(n1cc),ibloc,nsubitr,
     .    itrb1(igrid,ibloc),itrb2(igrid,ibloc),jtrb1(igrid,ibloc),
     .    jtrb2(igrid,ibloc),ktrb1(igrid,ibloc),ktrb2(igrid,ibloc),
     .    ntorder,tv1old(m1fgt(1,ibloc)),tv2old(m1fgt(1,ibloc)),dtphy,
     .    irkstg,nrkstg,mrkstg,
     .    rkphysa,rkphyse,dtv1ork(m1fgtrk(1,ibloc)),
     .                    dtv2ork(m1fgtrk(1,ibloc)) )
c----------------------------------------------------------------------
      endif
c----      end loop on blocks
 1000 continue
c
c
c----------------------------------------------------------------
c     fill-in the boundary information on eddy-viscosity
c     computed in the subroutine 'turb'
c
       if(iturb.le.0) go to 2000

      call bcturb (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,
     .  imp2,jmp2,kmp2,m1cc,m5cc,mgrlev,igrid,isoln,mbloc,nbloc,
     .  mxsizc,mxsiz5c,mx1dwk,imap,msegt,msegn,nseg,
     .                wk2d(1,1),eomu,turv1,turv2,w,p,smin,
     .  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .  tnuwfi1,tnuwfi2,tnuwfj1,tnuwfj2,tnuwfk1,tnuwfk2,m1is,m1js,
     .  m1ks,mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .  mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .  lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .  kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .  ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes )
c
 2000 continue
      init = init + 1
c     end of subiteration loop
c----------------------------------------------------------------
c
c
       return
       end
