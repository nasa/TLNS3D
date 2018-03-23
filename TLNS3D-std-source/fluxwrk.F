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
      subroutine fluxwrk(imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,
     .                 imp2,jmp2,kmp2,m1np,m2np,m3np,m4np,m5np,
     .                 m1cc,m2cc,m3cc,m4cc,m5cc,m1fg,m4fg,m5fg,m5cg,
     .                 m1is,m1js,m1ks,
     .                 m3is,m3js,m3ks,igrid,isoln,nbloc,mres,
     .                 mx1dwk,mx3dwk,mxsizc,mxsiz3c,mxsiz4c,mxsiz5c,
     .                 mxszfg,mxsz4fg,mxsz5fg,mxsz5cg,mxsizn,mxsiz3n,
     .                 imap,nseg,msegt,msegn,mbloc,mgrlev,
     .                 w,p,ws,wr,delw,fw,fv,eomu,turv1,turv2,vol,x,fqs,
     .                 six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                 ri1,ri2,rj1,rj2,rk1,rk2,
     .                 tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .                 uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
     .                 uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .                 tnuwfi1,tnuwfi2,tnuwfj1,tnuwfj2,tnuwfk1,tnuwfk2,
     .                 dtl,dti,dtj,dtk,dtmin,ratioij,ratiojk,eprs,
     .                 wk1d,wk2d,wk1d5,
     .                 rtrmsg,nitr,iforfn,ilast,ncyc,ibcup,ilev,nres,
     .                 ivisb,const1,const2,const3,const4,const5,mstrt,
     .                 smin,iturbb,itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,
     .                 mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .                 mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .                 lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .              kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .              ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,
     .              dtvi,dtvj,dtvk,fbeta2,fbetav2,wold,tv1old,tv2old,
     .              m1fgt,m5fgt,iturv1,iturv2,mtorder,ntorder,iunsteady,
     .              dtphy,totptim,m1fgtrk,m5fgtrk,irkstg,nrkstg,mrkstg,
     .              rkphysa,rkphyse,delwork,dtv1ork,dtv2ork )
c#ifdef  RKPHY
c# else
c    .              dtphy,totptim )
c# endif
c
c     object : To evaluate fluxes of viscous flow equations
c              for use in physical R-K and MEBDF time-stepping scheme
c
c     initial coding : by v.n.vatsa (Jan. 2001)
c     latest mods: Nov. 2003
c

      character*128  errmsg
      common/dstcom/ errmsg,master,myrank,mycomm,ndlist,nnodes
      dimension      nodes(1)
c
      common/conv/  rtmaxg,hmaxg,rtrms2, hrms2, totpts,
     .              irtmxg,jrtmxg,krtmxg,mrtbloc,
     .               ihmxg, jhmxg, khmxg, mhbloc
c
      common/rkerr/ erkmaxl(5),erkmaxg(5),erkrms2(5),erkrmsg(5),
     .              totptse
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      common/rk/    cfl,c(6),qfil(6),beta(6),vt,hm,mstage
c
      common/rkdis/ vis0,vis2,vis4,zeta,vepsn,vepsl,enteps,icau
c
      common/rkrsm/ smoopi,smoopj,smoopk,smoopic,smoopjc,smoopkc
c
c-------   preconditioning related quantities      -----------------
c
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubcon2,iprec
      common/unsprec/ npcons,nexp,iprecg,ijac,icord,inav
c
c
c
      dimension     im  (mgrlev,mbloc), jm  (mgrlev,mbloc),
     .              km  (mgrlev,mbloc), imp1(mgrlev,mbloc),
     .              jmp1(mgrlev,mbloc), kmp1(mgrlev,mbloc),
     .              imp2(mgrlev,mbloc), jmp2(mgrlev,mbloc),
     .              kmp2(mgrlev,mbloc)
c
      dimension     imap(msegt,msegn,mbloc,mgrlev), nseg(mbloc)
c
      dimension     ivisb(3,mbloc),      const1(msegn,mbloc),
     .              const2(msegn,mbloc),const3(msegn,mbloc),
     .              const4(msegn,mbloc),const5(msegn,mbloc)
c
      dimension     itrb1(mgrlev,mbloc), itrb2(mgrlev,mbloc),
     .              jtrb1(mgrlev,mbloc), jtrb2(mgrlev,mbloc),
     .              ktrb1(mgrlev,mbloc), ktrb2(mgrlev,mbloc), 
     .              iturbb(mbloc)
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
      dimension     m1fg(mbloc),m4fg(mbloc),m5fg(mbloc)
c
      dimension      w(mxsiz5c),vol(mxsizc ),p(mxsizc ),eomu(mxsizc),
     .                      smin(mxsizc), turv1(mxsizc),turv2(mxsizc)
c
c     array allocation for time dependent quantities
c
      dimension     wold (mxsz5fg*ntorder+1),
     .              tv1old(mxszfg*ntorder*iturv1+1),
     .              tv2old(mxszfg*ntorder*iturv2+1)
c
      dimension     m1fgt(mtorder,mbloc), m5fgt(mtorder,mbloc)
c
c#ifdef  RKPHY
      dimension     m1fgtrk(nrkstg ,mbloc), m5fgtrk(nrkstg ,mbloc)
c
      dimension     rkphysa(mrkstg,mrkstg),rkphyse(mrkstg)
c
      dimension     delwork (mxsz5fg*nrkstg+1),
     .              dtv1ork (mxszfg*nrkstg*iturv1+1),
     .              dtv2ork (mxszfg*nrkstg*iturv2+1)
c#endif
c
      dimension      x(mxsiz3n),six(mxsizn ),siy(mxsizn ),siz(mxsizn ),
     .                          sjx(mxsizn ),sjy(mxsizn ),sjz(mxsizn ),
     .                          skx(mxsizn ),sky(mxsizn ),skz(mxsizn )
c
      dimension      wr(mxsz5cg)
c
      dimension      ws(mxsz5fg), fw(mxsz5fg),  fv(mxsz4fg),
     .              dti(mxszfg ),dtj(mxszfg ), dtk(mxszfg )
cprec
      dimension    dtvi(mxsz5fg),   dtvj(mxsz5fg),  dtvk(mxsz4fg),
     .             fbeta2(mxszfg ), fbetav2(mxszfg )
c
      dimension     delw(mx3dwk,5)
c
      dimension     dtl(1),ratioij(1),ratiojk(1),eprs(1)
c
      dimension     fqs(mx1dwk,mx1dwk,6)
c
      dimension     ri1(1),ri2(1),rj1(1),rj2(1),rk1(1),rk2(1)
      dimension     tauwfi1(1),tauwfi2(1),tauwfj1(1),tauwfj2(1),
     .              tauwfk1(1),tauwfk2(1)
      dimension     uswfi1(1),uswfi2(1),uswfj1(1),uswfj2(1),
     .              uswfk1(1),uswfk2(1)
      dimension     uplwfi1(1),uplwfi2(1),uplwfj1(1),uplwfj2(1),
     .              uplwfk1(1),uplwfk2(1)
      dimension     tnuwfi1(1),tnuwfi2(1),tnuwfj1(1),tnuwfj2(1),
     .              tnuwfk1(1),tnuwfk2(1)
c
      dimension     rtrmsg(mres)
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
c     one-dimensional work-arrays  with 5-elements
c
      dimension   wk1d5(5*mx1dwk,5)
c
c     patching related arrays
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
c**************************************************************** 
c
c--------- if ksmoop.gt.0 smooth residuals at all stages
c          if ksmoop.lt.0 smooth residuals at alternate stages
c
      ksmoop    = 0
      if (smoopi.ne.0..or.smoopj.ne.0..or.smoopk.ne.0.) ksmoop = 1
      if (smoopi.lt.0..or.smoopj.lt.0..or.smoopk.lt.0.) ksmoop = -1
      ismoop    = iabs(ksmoop)
      if (mod(mstage,2).eq.0) ismoop = ksmoop
c
c--------------------------------------------------------------
c     compute time step on all blocks of the current  grid
c
c     initialize dtmin here instead of relying on the block
c     number to be 1 in step since on most nodes the block
c     number will never be 1
c
      dtmin = 1.e+10
      do 500 ibloc = 1,nbloc
c
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank)
# else
      if (nodes(ibloc).eq.myrank)
# endif
     .call step (imn,jmn,kmn,
     .   im  (igrid,ibloc),jm  (igrid,ibloc),km  (igrid,ibloc),
     .   imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .   imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),ibloc,
     .   w(m5cc(igrid,ibloc)),   p(m1cc(igrid,ibloc)),
     . vol(m1cc(igrid,ibloc)),   x(m3np(igrid,ibloc)),
     . six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     . siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     . sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     . skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     . skz(m1np(igrid,ibloc)),eomu(m1cc(igrid,ibloc)),dti(m1fg(ibloc)),
     . dtj(m1fg(ibloc)), dtk(m1fg(ibloc)), dtmin,ivisb(1,ibloc),
     . wk1d(1,1),wk1d(1,2),mx1dwk,
     . ntorder, dtphy,
     . dtvi(m1fg(ibloc)), dtvj(m1fg(ibloc)), dtvk(m1fg(ibloc)),
     . fbeta2(m1fg(ibloc)), fbetav2(m1fg(ibloc)) )
c
  500 continue
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c     if using global minimum time step, need to get minimum
c     across all the nodes
c
      if (vt.ne.0. .and. nnodes.gt.1)
     .  call getdtmin (dtmin,nodes)
#endif
c
c 
c---------- compute the fluxes for turbulence model on all blocks
c           via the driver routine fluxtrk
c
       call fluxtrk (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     . m1np,m2np,m3np,m4np,m5np,m1cc,m2cc,m3cc,m4cc,m5cc,m5cg,
     . m1is,m1js,m1ks,m3is,m3js,m3ks,mgrlev,igrid,isoln,mbloc,nbloc,
     . mxsizc,mxsiz4c,mxsiz5c,mxsizn,mxsiz3n,mx1dwk,mx3dwk,
     . imap(1,1,1,igrid),msegt,msegn,nseg,w,p,eomu,turv1,turv2,vol,x,
     . six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     . ri1,ri2,rj1,rj2,rk1,rk2,
     . uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     . tnuwfi1,tnuwfi2,tnuwfj1,tnuwfj2,tnuwfk1,tnuwfk2,
     . wk1d,wk2d,delw,
     . fv, dtl, dti, dtj, dtk, dtvi,dtvj,dtvk,
     . mxszfg, m1fg, ivisb,mstrt,smin,
     .             iturbb,itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,
     . mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     . mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     . lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     . kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     . ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,mtorder,
     . ntorder,iturv1,iturv2,m1fgt,tv1old,tv2old,dtphy ,
     . m1fgtrk,m5fgtrk,irkstg,nrkstg,mrkstg,
     . rkphysa,rkphyse,dtv1ork,dtv2ork )
c
c
      nstage  =  1
      jbloc   =  0
c
c
      do 520 ibloc = 1,nbloc
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) then
# else
      if (nodes(ibloc).eq.myrank) then
# endif
        jbloc = jbloc + 1
c
c     set up fqs array, which determines whether a given boundary
c     is an inviscid solid surface or not
c
        do 512 ia=1,6
        do 512 n2=1,mx1dwk
        do 512 n1=1,mx1dwk
        fqs(n1,n2,ia) = 1.0
  512   continue
c
        do 516 ns=1,nseg(ibloc)
        nbctype = imap(1,ns,ibloc,igrid)
c       if(nbctype.ge.2.and.nbctype.le.3) then
        if(nbctype.eq.2.or.nbctype.eq.102) then
          nface = imap(2,ns,ibloc,igrid)
          nbeg1 = imap(3,ns,ibloc,igrid)
          nend1 = imap(4,ns,ibloc,igrid)
          nbeg2 = imap(5,ns,ibloc,igrid)
          nend2 = imap(6,ns,ibloc,igrid)
c
          if(nface.le.2) then
             if( nbeg1 .ne. jmn            ) nbeg1 = nbeg1 +1
             if( nend1 .eq. jm(igrid,ibloc)) nend1 = nend1 +1
             if( nbeg2 .ne. kmn            ) nbeg2 = nbeg2 +1
             if( nend2 .eq. km(igrid,ibloc)) nend2 = nend2 +1
          else if(nface.ge.5) then
             if( nbeg1 .ne. imn            ) nbeg1 = nbeg1 +1
             if( nend1 .eq. im(igrid,ibloc)) nend1 = nend1 +1
             if( nbeg2 .ne. jmn            ) nbeg2 = nbeg2 +1
             if( nend2 .eq. jm(igrid,ibloc)) nend2 = nend2 +1
          else
             if( nbeg1 .ne. kmn            ) nbeg1 = nbeg1 +1
             if( nend1 .eq. km(igrid,ibloc)) nend1 = nend1 +1
             if( nbeg2 .ne. imn            ) nbeg2 = nbeg2 +1
             if( nend2 .eq. im(igrid,ibloc)) nend2 = nend2 +1
          end if
c
          do 514 n2=nbeg2,nend2
          do 514 n1=nbeg1,nend1
          fqs(n1,n2,nface) = 0.0
  514     continue
        end if
  516   continue
c
c      set up time step-related terms for this block
c
      call blkstp (imn,jmn,kmn,
     .     im  (igrid,ibloc),jm  (igrid,ibloc),km  (igrid,ibloc),
     .     imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .     imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     . vol(m1cc(igrid,ibloc)), dti(m1fg(ibloc)), dtj(m1fg(ibloc)),
     .         dtk(m1fg(ibloc)), dtvi(m1fg(ibloc)), dtvj(m1fg(ibloc)),
     .         dtvk(m1fg(ibloc)), dtl,ratioij,ratiojk,eprs  ,
#ifdef  RKPHY
cvn  .       ntorder, dtphy ,rkphysa(irkstg,irkstg) )
     .       ntorder, dtphy ,rkphysa(2,2) )
# else
     .       ntorder, dtphy )
# endif
c
c      compute dissipative and viscous fluxes
c
c
c         compute dissipative flux (fw) on the solution mesh
c
          call filtf (ncyc,nstage,imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,
     .         imp2,jmp2,kmp2,m1cc,mgrlev,igrid,mbloc,ibloc,mxsizc,
     .         mx1dwk,imap(1,1,ibloc,igrid),msegt,msegn,nseg(ibloc),p,
     .         im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .         imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .         imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),   p(m1cc(igrid,ibloc)),
     .     vol(m1cc(igrid,ibloc)),
     .     six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .     siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .     sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .     skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .     skz(m1np(igrid,ibloc)),       dti(m1fg(ibloc)),
     .           dtj(m1fg(ibloc)),       dtk(m1fg(ibloc)), 
     .           dtl, ratioij, ratiojk,
     .      wk1d5(1,1),wk1d5(1,2),wk1d5(1,3),wk1d5(1,4),wk1d5(1,5),
     .      wk1d(1,1),wk1d(1,2),wk1d(1,3),
     .      wk2d(1,1),wk2d(1,2),wk2d(1,3),wk2d(1,4),delw,
     .      fv(m4fg(ibloc)), fw(m5fg(ibloc)) ,
     .      dtvi(m1fg(ibloc)), dtvj(m1fg(ibloc)), dtvk(m1fg(ibloc)),
     .      fbetav2(m1fg(ibloc)) )

c
c         compute viscous flux (fv)
c
          call navier (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .    imap(1,1,ibloc,igrid), msegt, msegn, mbloc, nseg(ibloc),
     .      w(m5cc(igrid,ibloc)),   p(m1cc(igrid,ibloc)),
     .    vol(m1cc(igrid,ibloc)),   x(m3np(igrid,ibloc)),
     .    six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .    siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .    sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .    skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .    skz(m1np(igrid,ibloc)),eomu(m1cc(igrid,ibloc)),
     .    tauwfi1(m1is(igrid,ibloc)), tauwfi2(m1is(igrid,ibloc)),
     .    tauwfj1(m1js(igrid,ibloc)), tauwfj2(m1js(igrid,ibloc)),
     .    tauwfk1(m1ks(igrid,ibloc)), tauwfk2(m1ks(igrid,ibloc)),
     .    wk1d(1,1),wk1d(1,2),wk1d(1,3),wk1d(1,4),wk1d(1,5),
     .    wk2d(1,1),wk2d(1,2),wk2d(1,3),wk2d(1,4),wk2d(1,5),
     .    wk2d(1, 6),wk2d(1, 7),wk2d(1, 8),wk2d(1, 9),
     .    wk2d(1,10),wk2d(1,11),wk2d(1,12),wk2d(1,13), mx1dwk,
     .    delw(1,1),delw(1,2),delw(1,3),delw(1,4),delw(1,5),
     .     fv(m4fg(ibloc)), ivisb(1,ibloc) )
c
c
c     evaluate convective fluxes without dt scaling
c
      call deltaw (ncyc,nstage,igrid,imn,jmn,kmn,
     .  im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .  imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .  imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc), mx1dwk,
     .     w(m5cc(igrid,ibloc)),   p(m1cc(igrid,ibloc)),
     .   six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .   siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .   sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .   skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .   skz(m1np(igrid,ibloc)),fqs,
     .   wk1d(1, 1),wk1d(1, 2),wk1d(1, 3),wk1d(1, 4),wk1d(1, 5),
     .   wk1d(1, 6),wk1d(1, 7),wk1d(1, 8),wk1d(1, 9),wk1d(1,10),
     .   wk1d(1,11),wk1d(1,12),wk1d(1,13),wk1d(1,14),wk1d(1,15),
     .   wk2d(1,1),wk2d(1,2),wk2d(1,3),wk2d(1,4),wk2d(1,5), delw, 0  )
c
c     add physical dissipative (fv) flux to convective flux (delw)
c
         call adfluxv (imn,jmn,kmn,
     .   im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .   imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .   imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .   delw, fv(m4fg(ibloc)) )
c
c     add artificial dissipative (fw) flux to convective flux (delw)
c        ( conservative )
c
      if (iprec.ne.3)
     .   call adfluxd (nstage,imn,jmn,kmn,
     .   im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .   imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .   imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .   delw, fw(m5fg(ibloc)) )
c
c     put delw into delwork array
c
         call equw1w (imn,jmn,kmn,
     .         imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .         imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .         delw, delwork(m5fgtrk(irkstg,ibloc)),5 )
c
c    if at the last stage compute error estimate for
c    temporal error physical RK
c
#ifdef  RKPHY
c
       if (irkstg.eq.nrkstg)
     . call errphrk(imn,jmn,kmn,
     .   im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .   imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .              delw,dtl,vol(m1cc(igrid,ibloc)),dtphy,
cvn  .              mx1dwk,wk1d(1,1),ibloc,nodes(ibloc),
     .              mx1dwk,wk1d(1,1),jbloc,nodes(ibloc),
     .              delwork(m5fgtrk(1,ibloc)),
     .              rkphysa,rkphyse,mrkstg,nrkstg,irkstg)
c
#endif
c
c     end the local block test on do 520
      endif
  520 continue
c
      return
      end
