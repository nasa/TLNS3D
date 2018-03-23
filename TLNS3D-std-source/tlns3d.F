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
      subroutine tlns3d (mbloc,mgrlev,mxid,mxjd,mxkd,mx1d,mxsizc,mxsizn,
     .                   mxszcg,mxszfg,mxs2dc,mxsurf,mxpchs,mxtpchc,
     .                   mpchitm,mxpchc,mxpitm,mx3dwk,mx2dwk,
     .                   nbufct,nbufsz,mres,msegt,msegn,mgrleg,
     .                   mxidp3,mxjdp3,mxkdp3,mx1dwk,mxijp3,
     .                   mx1dwk5,mx3dwk3,mx3dwk4,mx3dwk5,mxsiz2c,
     .                   mxsiz3c,mxsiz4c,mxsiz5c,mxsz5cg,mxsz4fg,
     .                   mxsz5fg,mxsiz2n,mxsiz3n,mxsiz4n,mxsiz5n,
     .                   mxs2d3c,mxdstf,mxdsti,mxinpln,mxint1,
     .                   mxint2,mxintln ,mtorder,ntorder,iunsteady,
     .                   ngroup,iturv1,iturv2,mrkstg,nrkstg,mebdf,
     .                   ibdf2opt)

c 
c     analysis of transonic flow over 3-d configurations
c     by solution of unsteady thin-layer Navier-Stokes equations
c     finite volume method with central difference scheme
c     using full-multigrid acceleration
c     multiple blocks (structured) of grids accomodated
c     initial design and coding : by v.n.vatsa of nasa langley, feb.1990
c     
c     distributed version initial coding by b.w.wedan, 1995
c     asynchrnous version developed by d.p.hammond & v.n.vatsa, 1997
c
c     mods for yaw, new engine and multiple jet-stream
c     conditions : by v.n.vatsa, Jan. 1998
c     yaw angle: based on Youn Oh's (Hughes) work
c
c     Dynamic memory version: Sept. 1999
c     (follows the work of Hamid, Massey in PAB3D & Biedron in CFL3D)
c
c
c     time-dependent version mods : May 2000 by v.n.vatsa
c
c     w(i,j,k,1)  = density 
c     w(i,j,k,2)  = momentum in x direction 
c     w(i,j,k,3)  = momentum in y direction 
c     w(i,j,k,4)  = momentum in z direction 
c     w(i,j,k,5)  = total energy
c     p(i,j,k)    = pressure
c
c     The distributed version implementation consists of a master
c     or root process (which has rank 0) and any number of node
c     processes from 1 up to the total number of blocks. These
c     processes are number consectively starting from 1. The rank
c     of a process is given by the variable myrank, and the total
c     number of node processes by nnodes. These are both defined
c     in the common block dstcom.
c
c     The master process is responsible for all the input/output,
c     sending the input data to the nodes, and collecting the
c     output data from the nodes. It performs only limited
c     calculations. The parts of the code associated with the master
c     process are those identified by myrank.eq.master (0)
c
c     The node processes, on the otherhand, do all the calculations
c     but no input or output. The code is structured such that any
c     number of blocks may run on any node. This allows many small
c     blocks to be collected together with computations done by
c     a single process. In order to identify which blocks are
c     associated with which process, the array nodes(mblock) is used.
c     The value of nodes(ibloc) is the rank of the process for which
c     the block, ibloc, is local. Thus, the statement:
c         if (nodes(ibloc).eq.myrank)
c     will be true only if the block is local.
c
c     The above test is used by the program in all the loops on the
c     blocks to determine which blocks are to be calculated by
c     which processes.
c

c****************************************************************
c     parameters  for order of time-accurate calculations
c
c     ntorder  = 0 steady state calculation
c     ntorder >= 1 for 1st order
c                2 for 2nd order
c                3 for 3rd order
c                4 for 4th order
c                5 for 5th order
c                6 for 6th order
c
c     mtorder  = max(1,ntorder)  : used for dimensioning of 2-d pointers
c****************************************************************
c****************************************************************
c
c     parameters  for multigrid on multiple blocks 
c
c     (The code is  now self-sizing, and does not need parms.h file)
c      All the required parameters are computed in main.F and
c      passed to tlns3d routine through the argument list
c
c****************************************************************
c     maximum dimensions in terms of node-points for any one block
c
c     note: for 2-d option mxkd should be set equal to 2
c
c---- warning :   mx1d should be set = max(mxid,mxjd,mxkd) -----
c
c     parameter  (mxidp3= mxid+3, mxjdp3= mxjd+3, mxkdp3= mxkd+3)
c     parameter  (mxijp3  = mxidp3*mxjdp3       )
c     parameter  (mx1dwk  = mx1d+3  )
c     parameter  (mx3dwk  = mxijkp3 ,mx1dwk5 = mx1dwk*5)
c     parameter  (mx3dwk3 = mx3dwk*3,mx3dwk4 = mx3dwk*4)
c     parameter  (mx3dwk5 = mx3dwk*5)
c
c****************************************************************
c     maximum size allowed for 3-d arrays (with 1-5 subelements)
c     on all blocks for all grid levels used in multi-grid
c     ( this is a cumulative number and includes ghost cells)
c     mxsizc is for quantities stored at cell-centers
c     mxsizn is for quantities stored at nodes
c
c     parameter  (mxsiz2c = mxsizc*2, mxsiz3c = mxsizc*3)
c     parameter  (mxsiz4c = mxsizc*4, mxsiz5c = mxsizc*5)
c     parameter  (mxsz5cg = mxszcg*5)
c     parameter  (mxsz4fg = mxszfg*4, mxsz5fg = mxszfg*5)
c     parameter  (mxsiz2n = mxsizn*2, mxsiz3n = mxsizn*3)
c     parameter  (mxsiz4n = mxsizn*4, mxsiz5n = mxsizn*5)
c     parameter  (mxs2d3c = mxs2dc*3)
c
c     parameters for new distance function
c
c     parameter  (mxdstf = mxsz4fg + 14*mxsurf)
c     parameter  (mxdsti = mxszfg  + 12*mxsurf)
c
c****************************************************************
c     control parameters for multiblock/multigrid
c----------------------------------------------------------------
c
      common/conv/  rtmaxg,hmaxg,rtrms2, hrms2, totpts,
     .              irtmxg,jrtmxg,krtmxg,mrtbloc,
     .               ihmxg, jhmxg, khmxg, mhbloc
c
      common/rkerr/ erkmaxl(5),erkmaxg(5),erkrms2(5),erkrmsg(5),
     .              totptse
c
      dimension     im  (mgrlev,mbloc), jm  (mgrlev,mbloc),
     .              km  (mgrlev,mbloc), imp1(mgrlev,mbloc),
     .              jmp1(mgrlev,mbloc), kmp1(mgrlev,mbloc),
     .              imp2(mgrlev,mbloc), jmp2(mgrlev,mbloc),
     .              kmp2(mgrlev,mbloc),
     .              idum(mbloc), jdum(mbloc), kdum(mbloc)
c
      pointer      (ip_im  ,im  ), (ip_jm  ,jm  ),
     .             (ip_km  ,km  ), (ip_imp1,imp1),
     .             (ip_jmp1,jmp1), (ip_kmp1,kmp1),
     .             (ip_imp2,imp2), (ip_jmp2,jmp2),
     .             (ip_kmp2,kmp2), (ip_idum,idum),
     .             (ip_jdum,jdum), (ip_kdum,kdum) 
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
     .              m3ks(mgrlev,mbloc),
     .              m1fg(mbloc),m4fg(mbloc),m5fg(mbloc),
     .              m1fgt  (mtorder,mbloc), m5fgt  (mtorder,mbloc)
c#ifdef  RKPHY
      dimension     m1fgtrk(nrkstg ,mbloc), m5fgtrk(nrkstg ,mbloc)
      pointer      (ip_m1fgtrk,m1fgtrk), (ip_m5fgtrk,m5fgtrk)
c#endif
c
      pointer      (ip_m1np,m1np), (ip_m2np,m2np),
     .             (ip_m3np,m3np), (ip_m4np,m4np),
     .             (ip_m5np,m5np),
     .             (ip_m1cc,m1cc), (ip_m2cc,m2cc),
     .             (ip_m3cc,m3cc), (ip_m4cc,m4cc),
     .             (ip_m5cc,m5cc), (ip_m5cg,m5cg),
     .             (ip_m1is,m1is), (ip_m1js,m1js),
     .             (ip_m1ks,m1ks), 
     .             (ip_m3is,m3is), (ip_m3js,m3js),
     .             (ip_m3ks,m3ks), (ip_m1fg,m1fg),
     .             (ip_m4fg,m4fg), (ip_m5fg,m5fg),
     .             (ip_m1fgt  ,m1fgt  ), (ip_m5fgt  ,m5fgt  ) 
c
c
      dimension     itrb1(mgrlev,mbloc), itrb2(mgrlev,mbloc),
     .              jtrb1(mgrlev,mbloc), jtrb2(mgrlev,mbloc),
     .              ktrb1(mgrlev,mbloc), ktrb2(mgrlev,mbloc), 
     .              iturbb(mbloc)
c
      pointer      (ip_itrb1,itrb1), (ip_itrb2,itrb2),
     .             (ip_jtrb1,jtrb1), (ip_jtrb2,jtrb2),
     .             (ip_ktrb1,ktrb1), (ip_ktrb2,ktrb2),
     .             (ip_iturbb,iturbb)
c
c
      dimension     icoars(mgrleg,mgrlev), ifine(mgrleg,mgrlev),
     .              ibegr (mgrleg,mgrlev),
     .              mcyc(mgrlev) ,isoln(mgrlev), nrleg(mgrlev),
     .              gitr(mgrlev) ,gitp (mgrlev),mbcupd(mgrlev),
     .              mcnvout(mgrlev), nsteps(mgrlev)
c
      pointer      (ip_icoars,icoars), (ip_ifine  ,ifine ),
     .             (ip_ibegr ,ibegr ), (ip_mcyc   ,mcyc  ),
     .             (ip_isoln ,isoln ), (ip_nrleg  ,nrleg ),
     .             (ip_gitr  ,gitr  ), (ip_gitp   ,gitp  ),
     .             (ip_mbcupd,mbcupd), (ip_mcnvout,mcnvout),
     .             (ip_nsteps,nsteps)
c
c
      dimension     imap(msegt,msegn,mbloc,mgrlev), nseg(mbloc)
c
      pointer      (ip_imap ,imap ), (ip_nseg ,nseg )
c
c-------   convergence related quantities      -----------------
c     (last 5 lines are added for yaw/roll capability)
c
      dimension     wcount(mres),  rtrmsg(mres), 
     .              hrmsg (mres),  nsupg (mres),
     .              resout(mres),  supout(mres),
     .              clg   (mres),  cmxg  (mres),
     .              cmyg  (mres),  cmzg  (mres),
     .              cdtg  (mres),  cdvg  (mres),
     .              czg   (mres),  cz    (mbloc),
     .              cl    (mbloc), cmx   (mbloc),
     .              cmy   (mbloc), cmz   (mbloc),
     .              cdt   (mbloc), cdv   (mbloc),
     .              cx    (mbloc), cxg   (mres), cng(mres),
     .              cy    (mbloc), cyg   (mres),
     .              cxv   (mbloc), cxvg  (mres),
     .              cyv   (mbloc), cyvg  (mres),
     .              czv   (mbloc), czvg  (mres)
c
      pointer      (ip_wcount,wcount), (ip_rtrmsg,rtrmsg),
     .             (ip_hrmsg ,hrmsg ), (ip_nsupg ,nsupg ),
     .             (ip_resout,resout), (ip_supout,supout),
     .             (ip_clg   ,clg   ), (ip_cmxg  ,cmxg  ),
     .             (ip_cmyg  ,cmyg  ), (ip_cmzg  ,cmzg  ),
     .             (ip_cdtg  ,cdtg  ), (ip_cdvg  ,cdvg  ),
     .             (ip_czg   ,czg   ), (ip_cz    ,cz    ),
     .             (ip_cl    ,cl    ), (ip_cmx   ,cmx   ),
     .             (ip_cmy   ,cmy   ), (ip_cmz   ,cmz   ),
     .             (ip_cdt   ,cdt   ), (ip_cdv   ,cdv   ),
     .             (ip_cx    ,cx    ), (ip_cxg   ,cxg   ),
     .             (ip_cng   ,cng   ),
     .             (ip_cy    ,cy    ), (ip_cyg   ,cyg   ),
     .             (ip_cxv   ,cxv   ), (ip_cxvg  ,cxvg  ),
     .             (ip_cyv   ,cyv   ), (ip_cyvg  ,cyvg  ),
     .             (ip_czv   ,czv   ), (ip_czvg  ,czvg  )

c
c-------   patching related quantities      -----------------
c
      common/facetr/ ifacetr(6),ipatchg
c
c
      dimension   ntpchcb(mbloc,mgrlev), npchcbf(6,mbloc,mgrlev),
     .            m1pch1 (mbloc,mgrlev), m1pch2 (mbloc,mgrlev),
     .            ipatchc(mxtpchc), jpatchc(mxtpchc),
     .            kpatchc(mxtpchc)
c
      pointer  (ip_ntpchcb,ntpchcb), (ip_npchcbf,npchcbf),
     .         (ip_m1pch1 ,m1pch1 ), (ip_m1pch2 ,m1pch2 ),
     .         (ip_ipatchc,ipatchc), (ip_jpatchc,jpatchc),
     .         (ip_kpatchc,kpatchc)
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
     .          frc    (mpchitm)
     .         ,iitmsa (mxpchs,mgrlev)
c
      pointer  (ip_nswpchb,nswpchb), (ip_npchitm,npchitm),
     .         (ip_lswpchb,lswpchb), (ip_lswpche,lswpche),
     .         (ip_lspchb1,lspchb1), (ip_lspchf1,lspchf1),
     .         (ip_lspchb2,lspchb2), (ip_lspchf2,lspchf2),
     .         (ip_ipitmb1,ipitmb1), (ip_ipitmb2,ipitmb2),
     .         (ip_jpitmb1,jpitmb1), (ip_jpitmb2,jpitmb2),
     .         (ip_kpitmb1,kpitmb1), (ip_kpitmb2,kpitmb2),
     .         (ip_frc    ,frc    ), (ip_iitmsa ,iitmsa )
c
c     mods for low-memory patching arrays
#if defined(BUILD_PVM) || defined(BUILD_MPI)

      dimension ipitmbs(mpchitm),jpitmbs(mpchitm),kpitmbs(mpchitm)
c
      pointer  (ip_ipitmbs,ipitmbs), (ip_jpitmbs,jpitmbs),
     .         (ip_kpitmbs,kpitmbs)
#else

      dimension ipitmbs(mpchitm),jpitmbs(mpchitm),kpitmbs(mpchitm)
c
      pointer  (ip_ipitmbs,ipitmbs), (ip_jpitmbs,jpitmbs),
     .         (ip_kpitmbs,kpitmbs)

#endif

c
c-------   preconditioning related quantities      -----------------
c     Define the Namelist "auxinput" for additional inputs
c
      NAMELIST/auxinput/ diffac,rsvis,omega,betvis,betvisv,
     .                   anin,aninb,bninb,
     .                   npcons,nexp,ijac,icord,inav,iprbet
c
c        diffac is factor for viscous addition to time step
c        rsvis used in blkstp for viscous correction to residual smoothing
c        omega used in blkstp for residual smoothing
c        betvis=0 no viscous contribution to beta^2
c        anin used for physical time step contribution to time step
c        npcons used in advflow to choose update option
c        nexp used in advflow for explicit treatment of physical time step
c             only used for iprec=2
c        ijac=1 use point Jacobi preconditioning
c        icord=1 use Cord Rossow's "implicit smoothing"
c        icord=-1 use Cord Rossow's "implicit smoothing" with kstage
c        inav=1 call Navier on coarse grids
c        iprbet=1 calls betaout
c
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubcon2,iprec
      common/timek/  ct(9)
      common/vis   / diffac,rsvis,omega
      common/unscom1/ anin,aninb,bninb,betvis,betvisv
      common/unsprec/ npcons,nexp,iprecg,ijac,icord,inav
c
c-------   turbulence related quantities      -----------------
c
      common/trbfs/ anutinf,zkinf,ominf
c
cold  common/radbc/ phirad
c
c***********************************************************************
c
c     Instead of using the variable names twotref,cmassw,
c     cporous, pchambr and epsrad, which are meaningful only
c     for viscous walls, we are renaming these 5 variables
c     that are read from map file as const1.....const5
c     with the following meanings (Jan. 1998)
c
c     For viscous walls (bctyp=3) :
c         const1   = twotref
c         const2   = cmassw
c         const3   = cporous
c         const4   = pchambr
c         const5   = epsrad
c
c     For extrapolation condition (bctyp=5) :
c
c         const1   = pback (back pressure)
c        (if const1 is zero, pressure is also extrapolated)
c
c     For fan-upstream engine boundary condition (bctyp=8):
c
c         const1   = epsfu (stream-tube contraction ratio)
c
c     For engine/fan (massflow) boundary condition (bctyp=9):
c
c         const1   = epseng (engine/fan massflow to free-stream ratio)
c         const2   = relax (relaxation factor >= 0.05 to 1.0)
c
c     For jet/nozzle boundary condition (bctyp=10) :
c
c         const1   = ptjet (jet total pressure/ref. static pressure)
c         const2   = ttjet (jet total temp./ref. static temp.)
c         const3   = rmjet (nominal jet Mach number)
c         const4   = alfaj (jet inclination angle in x-y plane)
c         const5   = betaj (jet inclination angle in x-z plane)
c        (alfaj, betaj are usually set to zero)
c
c     For flow-streams with different total pressure & temp. (bctyp=11)
c         and specified inflow angles
c
c         const1   = ptstr (stream total pressure/ref. static pressure)
c         const2   = ttstr (stream total temp./ref. static temp.)
c         const3   = relax (relaxation factor : default = 1.0)
c         const4   = alfstr (stream inclination angle in x-y plane)
c         const5   = betstr (stream inclination angle in x-z plane)
c        (alfstr, betstr are usually set to zero)
c
c     For flow-streams with different total pressure & temp. (bctyp=12)
c         and inflow angles aligned with grid lines
c
c         const1   = ptstr (stream total pressure/ref. static pressure)
c         const2   = ttstr (stream total temp./ref. static temp.)
c         const3   = relax (relaxation factor : default = 1.0)
c
      dimension     ivisb(3,mbloc),      const1(msegn,mbloc),
     .              const2(msegn,mbloc),const3(msegn,mbloc),
     .              const4(msegn,mbloc),const5(msegn,mbloc)
c
      pointer      (ip_ivisb ,ivisb ), (ip_const1,const1),
     .             (ip_const2,const2), (ip_const3,const3),
     .             (ip_const4,const4), (ip_const5,const5)
c****************************************************************************
c
      dimension     w (mxsiz5c),p(mxsizc ),vol(mxsizc ),eomu(mxsizc),
     .              smin(mxsizc*iturv1+1),turv1(mxsizc*iturv1+1),
     .              turv2(mxsizc*iturv2+1)
c
      pointer      (ip_w    ,w    ), (ip_p    ,p    ),
     .             (ip_vol  ,vol  ), (ip_eomu ,eomu ),
     .             (ip_smin ,smin ), (ip_turv1,turv1),
     .             (ip_turv2,turv2)
c
c     array allocation for time dependent quantities
c
      dimension     wold (mxsz5fg*ntorder+1),
     .              tv1old(mxszfg*ntorder*iturv1+1),
     .              tv2old(mxszfg*ntorder*iturv2+1)
      pointer      (ip_wold ,wold ), (ip_tv1old  ,tv1old  ),
     .                               (ip_tv2old  ,tv2old  )
c
c     array allocation for mebdf related quantities
c
      dimension     wbar (mxsz5fg*(ntorder-1)*mebdf+1),
     .              tv1bar(mxszfg*(ntorder-1)*mebdf*iturv1+1),
     .              tv2bar(mxszfg*(ntorder-1)*mebdf*iturv2+1)
      pointer      (ip_wbar ,wbar ), (ip_tv1bar  ,tv1bar  ),
     .                               (ip_tv2bar  ,tv2bar  )
c
c     array allocation for time averaged quantities
c
      dimension     wavg (mxsz5fg),pavg(mxszfg),eomuavg(mxszfg)
      pointer      (ip_wavg ,wavg ), (ip_pavg ,pavg ),
     .             (ip_eomuavg ,eomuavg)
c
c#ifdef  RKPHY
      dimension     delwork (mxsz5fg*nrkstg+1),
     .              dtv1ork (mxszfg*nrkstg*iturv1+1),
     .              dtv2ork (mxszfg*nrkstg*iturv2+1)
c
      dimension     rkphysa(mrkstg,mrkstg),rkphyse(mrkstg)
c
      pointer      (ip_delwork ,delwork ), (ip_dtv1ork  ,dtv1ork  ),
     .                                     (ip_dtv2ork  ,dtv2ork  )
c
      pointer      (ip_rkphysa, rkphysa ), (ip_rkphyse  ,rkphyse  )
c
c#endif
c
c     array allocation for time dependent output quantities
c
      dimension    ibeggrp(ngroup),iendgrp(ngroup),iskpgrp(ngroup),
     .             jbeggrp(ngroup),jendgrp(ngroup),jskpgrp(ngroup),
     .             kbeggrp(ngroup),kendgrp(ngroup),kskpgrp(ngroup),
     .             iblkgrp(ngroup)
c
      pointer      (ip_ibeggrp,ibeggrp), (ip_iendgrp,iendgrp),
     .             (ip_iskpgrp,iskpgrp), (ip_jbeggrp,jbeggrp),
     .             (ip_jendgrp,jendgrp), (ip_jskpgrp,jskpgrp),
     .             (ip_kbeggrp,kbeggrp), (ip_kendgrp,kendgrp),
     .             (ip_kskpgrp,kskpgrp), (ip_iblkgrp,iblkgrp) 
c
# ifdef BUILD_MPI
      dimension buff(mxsiz3n)
      dimension buffw(5*mxszfg),buffp(mxszfg),buffe(mxszfg),
     . buff1(mxszfg),buff2(mxszfg),buffs(mxszfg)

      integer tbuf1,tbuf2,tbuf3,tbuf4,tbuf5,tbuf6,tbuf7,tbuf8,tbuf9
      dimension tbuf1(mxpchc),tbuf2(mxpchc),tbuf3(mxpchc)
      dimension tbuf4(mxpitm),tbuf5(mxpitm),tbuf6(mxpitm),
     .  tbuf7(mxpitm),tbuf8(mxpitm),tbuf9(mxpitm),tbuf10(mxpitm)
c
      pointer      (ip_buff  ,buff  ), (ip_buffw ,buffw ),
     .             (ip_buffp ,buffp ), (ip_buffe ,buffe ),
     .             (ip_buff1 ,buff1 ), (ip_buff2 ,buff2 ),
     .             (ip_buffs ,buffs ), (ip_tbuf1 ,tbuf1 ),
     .             (ip_tbuf2 ,tbuf2 ), (ip_tbuf3 ,tbuf3 ),
     .             (ip_tbuf4 ,tbuf4 ), (ip_tbuf5 ,tbuf5 ),
     .             (ip_tbuf6 ,tbuf6 ), (ip_tbuf7 ,tbuf7 ),
     .             (ip_tbuf8 ,tbuf8 ), (ip_tbuf9 ,tbuf9 ),
     .             (ip_tbuf10,tbuf10)

      dimension buffwo(mx3dwk5*ntorder+1),
     .          bufftv1(mx3dwk*ntorder*iturv1+1),
     .          bufftv2(mx3dwk*ntorder*iturv2+1)
      pointer      (ip_buffwo ,buffwo ), (ip_bufftv1 ,bufftv1),
     .                                   (ip_bufftv2 ,bufftv2)
# endif
      dimension     x(mxsiz3n),six(mxsizn),siy(mxsizn),siz(mxsizn),
     .                         sjx(mxsizn),sjy(mxsizn),sjz(mxsizn),
     .                         skx(mxsizn),sky(mxsizn),skz(mxsizn)
c
      pointer      (ip_x   ,x   ), (ip_six ,six ),
     .             (ip_siy ,siy ), (ip_siz ,siz ),
     .             (ip_sjx ,sjx ), (ip_sjy ,sjy ),
     .             (ip_sjz ,sjz ), (ip_skx ,skx ),
     .             (ip_sky ,sky ), (ip_skz ,skz )
c
c
      dimension     w1(mxsz5cg),wr(mxsz5cg)
c
      pointer      (ip_w1  ,w1  ), (ip_wr  ,wr  ) 
c
      dimension     ws(mxsz5fg), fw(mxsz5fg),  fv(mxsz4fg),
     .              dti(mxszfg ),dtj(mxszfg ), dtk(mxszfg )
c
      pointer      (ip_ws  ,ws  ), (ip_fw  ,fw  ), 
     .             (ip_fv  ,fv  ), (ip_dti ,dti ),
     .             (ip_dtj ,dtj ), (ip_dtk ,dtk ) 
c
cprec
      dimension     dtvi(mxsz5fg),   dtvj(mxsz5fg),  dtvk(mxsz4fg),
     .              fbeta2(mxszfg ), fbetav2(mxszfg )
c
      pointer      (ip_dtvi   ,dtvi   ), (ip_dtvj  ,dtvj  ), 
     .             (ip_dtvk   ,dtvk   ), (ip_fbeta2 ,fbeta2 ),
     .                                   (ip_fbetav2 ,fbetav2 )
c
      dimension     delw(mx3dwk,5),   dtl(mx3dwk ),
     .              eprs(mx3dwk3),ratioij(mx3dwk ),ratiojk(mx3dwk )
c
      pointer      (ip_delw   ,delw   ), (ip_dtl     ,dtl     ), 
     .             (ip_eprs   ,eprs   ), (ip_ratioij ,ratioij ),
     .             (ip_ratiojk,ratiojk)
c
      dimension     ri1(mxs2d3c),ri2(mxs2d3c),rj1(mxs2d3c),
     .              rj2(mxs2d3c),rk1(mxs2d3c),rk2(mxs2d3c) 
c
      pointer      (ip_ri1    ,ri1    ), (ip_ri2   ,ri2   ), 
     .             (ip_rj1    ,rj1    ), (ip_rj2   ,rj2   ),
     .             (ip_rk1    ,rk1    ), (ip_rk2   ,rk2   )
c
c     work arrays for wall-functions
c
      dimension     tauwfi1(mxs2dc), uswfi1(mxs2dc),  uplwfi1(mxs2dc),
     .              tauwfi2(mxs2dc), uswfi2(mxs2dc),  uplwfi2(mxs2dc),
     .              tauwfj1(mxs2dc), uswfj1(mxs2dc),  uplwfj1(mxs2dc),
     .              tauwfj2(mxs2dc), uswfj2(mxs2dc),  uplwfj2(mxs2dc),
     .              tauwfk1(mxs2dc), uswfk1(mxs2dc),  uplwfk1(mxs2dc),
     .              tauwfk2(mxs2dc), uswfk2(mxs2dc),  uplwfk2(mxs2dc), 
     .              tnuwfi1(mxs2dc), tnuwfi2(mxs2dc), tnuwfj1(mxs2dc),
     .              tnuwfj2(mxs2dc), tnuwfk1(mxs2dc), tnuwfk2(mxs2dc)
c
      pointer      (ip_tauwfi1, tauwfi1), (ip_uswfi1, uswfi1),
     .             (ip_tauwfi2, tauwfi2), (ip_uswfi2, uswfi2),
     .             (ip_tauwfj1, tauwfj1), (ip_uswfj1, uswfj1),
     .             (ip_tauwfj2, tauwfj2), (ip_uswfj2, uswfj2),
     .             (ip_tauwfk1, tauwfk1), (ip_uswfk1, uswfk1),
     .             (ip_tauwfk2, tauwfk2), (ip_uswfk2, uswfk2) 
c
      pointer      (ip_uplwfi1, uplwfi1), (ip_uplwfi2, uplwfi2),
     .             (ip_uplwfj1, uplwfj1), (ip_uplwfj2, uplwfj2),
     .             (ip_uplwfk1, uplwfk1), (ip_uplwfk2, uplwfk2) 
c
      pointer      (ip_tnuwfi1, tnuwfi1), (ip_tnuwfi2, tnuwfi2),
     .             (ip_tnuwfj1, tnuwfj1), (ip_tnuwfj2, tnuwfj2),
     .             (ip_tnuwfk1, tnuwfk1), (ip_tnuwfk2, tnuwfk2) 
c
c     work arrays for distance function
c
      dimension     fwrk(mxdstf), iwrk(mxdsti)
c
      pointer      (ip_fwrk   ,fwrk   ), (ip_iwrk  ,iwrk  )
c
c****************************************************************
c     dimension statements for temporary 1-, 2- and 3-D work-arrays
c
c     needed large contiguous data block of size (mx2dwk,15)
c     for data buffering for sends and receives. Thus concatenated
c     all the 2-D buffers into 1 space and equivalenced the
c     original variables to this data block. Also did the same
c     for the 1-D and 3-D buffers.
c
      dimension   wk1d(mx1dwk,27)
c
      pointer    (ip_wk1d  ,wk1d  )
c
c
c
c---------- allocating square 2-d space --------------------------
c
      dimension   wk2d(mx2dwk,15)
c
      pointer    (ip_wk2d  ,wk2d  )
c
c
c
c     space allocation for 5-element, one-dimensional work-arrays 
c
      dimension   wk1d5(mx1dwk5,5)
c
      pointer    (ip_wk1d5  ,wk1d5  )
c
c
c
      dimension   fqs(mx1dwk,mx1dwk,6)
c
      pointer    (ip_fqs  ,fqs  )
c
c dana 060895: changed time to always be SP; dimension to real*4
      real*4 tim(3,3),tsum(3)
      dimension twk1da (4)
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c---- distributed computing related ----------------------------
c
#  ifdef BUILD_PVM
#     include "fpvm3.h"
#     ifdef TLN_DBL
#       define RTYPE REAL8
#     else
#       define RTYPE REAL4
#     endif
#  else
      include "mpif.h"
#     ifdef TLN_DBL
#       define RTYPE MPI_DOUBLE_PRECISION
#     else
#       define RTYPE MPI_REAL
#     endif
      dimension istat(MPI_STATUS_SIZE)
#  endif
#     include "tags.h"
c
c     parameter (mxint1 = 10*mbloc+2)
c     parameter (mxint2 = 6*mbloc*mgrlev)
c     parameter (mxintln =
c    .  ((mxint1/mxint2)*mxint1 + (mxint2/mxint1)*mxint2) /
c    .  ((mxint1/mxint2)        + (mxint2/mxint1)       ))
c
      dimension intval(mxintln),nodeid(mbloc)
c
      pointer  (ip_intval,intval), (ip_nodeid,nodeid)
#else
c
c     need these defined for in-core version
c
      integer TAG_ERROR,TAG_OK
      parameter (TAG_ERROR=0,TAG_OK=1)
#endif
c
cold  parameter (mxinpln = 15+mgrlev)
c
      dimension fpval(8,mxinpln)
c
      pointer  (ip_fpval,fpval)
c
c***************************************************************
c     common statements for global variables
c
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c
      common/rk/    cfl,c(6),qfil(6),beta(6),vt,hm,mstage
c
      common/ma/    amachg
c
      common/rkdis/ vis0,vis2,vis4,zeta,vepsn,vepsl,enteps,icau
c
      common/rkrsm/ smoopi,smoopj,smoopk,smoopic,smoopjc,smoopkc
c
c---------------------------------------------------------------------
c     Common block for distributed computing.
c     The variables in the common block are:
c
c         master  - the ID of the master process
c         myrank  - the ID for the local process
c         mycomm  - communicator (only used by MPI)
c         ndlist  - offset in nodes array for list
c                   of process ID's in order (set to mbloc)
c         nnodes  - number of node processes (including the master)
c         errmsg  - error message string
c         nodes   - the array of block-to-node mappings and the
c                   list of node ID's. The first nbloc entries
c                   contain the node ID's for each block. The
c                   entries (ndlist+1,...ndlist+nnodes) contain
c                   the list of node ID's
c
cNOTE:The variable nodes() must be the last variable in the common
c     block dstcom, since in the subroutines in which it is used
c     the size is unknown and declared as nodes(1). Only in the main
c     routine is the size given as nodes(mbloc*2). This means if you
c     want to add something to the common block, add it BEFORE nodes.
c
c     Sept. 1999 : In the dynamic version, nodes is passed
c                  through argument list instead of common/dstcom/
c
      character*128  errmsg
      common/dstcom/ errmsg,master,myrank,mycomm,ndlist,nnodes
      dimension      nodes(mbloc*2)
c
c***************************************************************
c     data statements to initialize selected variables
c
      data          imn/2/, jmn/2/, kmn/2/
c****************************************************************
c 
      character*80 title,grdfil,rstfil,pltfil,savfil,mapfil,
     .             pchfil1,pchfil2,rhtall,rhtphy,
     .             gpltfil,qpltfil,datfil,pntfil,avgfil
c
c---- unit assignments
c
c     igrdfil = unit for grid input
c     irstfil = unit for restart file input
c     isavfil = unit for solution output
c     iread   = unit for input data
c     iwrit   = unit for print output
c     ipltfil = unit for plot output (used for RK errors)
c     imapfil = unit for map file input
c
c     additional files for time-dependent version
c     irhtall = unit fot complete residual history at all time/pseudo times
c     irhtphy = unit for residual history at physical time steps
c     igpfil  = unit for storing plot3d grid for time-dependent solution
c               (all node points at selected physical time-steps)
c     iqpfil  = unit for storing plot3d solution for time-dependent solution
c               (all node points at selected physical time-steps)
c     idatfil = unit for storing x,y,z,w(1-5),p
c               (selected node points at all physical time-steps)
c     ipntfil = unit for specifying at which ibloc,i,j,k the info. is
c               to be written on idatfil (ngroup set of values)
c
c     iavgfil = unit for writing average flow quantities
c
      ndlist     = mbloc
c
c***************** assign storage for arrays through pointers ***********
c
      icumsiz    =  0
      iname      = 20

      call umalloc (ip_im,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_jm,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_km,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_imp1,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_jmp1,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_kmp1,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_imp2,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_jmp2,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_kmp2,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_idum,mbloc,1,iname,icumsiz)
      call umalloc (ip_jdum,mbloc,1,iname,icumsiz)
      call umalloc (ip_kdum,mbloc,1,iname,icumsiz)
c
      call umalloc (ip_m1np,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m2np,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m3np,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m4np,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m5np,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m1cc,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m2cc,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m3cc,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m4cc,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m5cc,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m5cg,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m1is,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m1js,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m1ks,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m3is,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m3js,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m3ks,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m1fg,mbloc,1,iname,icumsiz)
      call umalloc (ip_m4fg,mbloc,1,iname,icumsiz)
      call umalloc (ip_m5fg,mbloc,1,iname,icumsiz)
      call umalloc (ip_m1fgt,mtorder*mbloc,1,iname,icumsiz)
      call umalloc (ip_m5fgt,mtorder*mbloc,1,iname,icumsiz)
c
c#ifdef  RKPHY
      call umalloc (ip_m1fgtrk,nrkstg*mbloc,1,iname,icumsiz)
      call umalloc (ip_m5fgtrk,nrkstg*mbloc,1,iname,icumsiz)
c#endif
c
      call umalloc (ip_itrb1,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_itrb2,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_jtrb1,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_jtrb2,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_ktrb1,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_ktrb2,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_iturbb,mbloc,1,iname,icumsiz)
c
      call umalloc (ip_icoars,mgrlev*mgrleg,1,iname,icumsiz)
      call umalloc (ip_ifine ,mgrlev*mgrleg,1,iname,icumsiz)
      call umalloc (ip_ibegr ,mgrlev*mgrleg,1,iname,icumsiz)
      call umalloc (ip_mcyc  ,mgrlev,1,iname,icumsiz)
      call umalloc (ip_isoln ,mgrlev,1,iname,icumsiz)
      call umalloc (ip_nrleg ,mgrlev,1,iname,icumsiz)
      call umalloc (ip_gitr  ,mgrlev,0,iname,icumsiz)
      call umalloc (ip_gitp  ,mgrlev,0,iname,icumsiz)
      call umalloc (ip_mbcupd,mgrlev,1,iname,icumsiz)
      call umalloc (ip_mcnvout,mgrlev,1,iname,icumsiz)
      call umalloc (ip_nsteps,mgrlev,1,iname,icumsiz)
c
      call umalloc (ip_ibeggrp,ngroup,1,iname,icumsiz)
      call umalloc (ip_iendgrp,ngroup,1,iname,icumsiz)
      call umalloc (ip_iskpgrp,ngroup,1,iname,icumsiz)
      call umalloc (ip_jbeggrp,ngroup,1,iname,icumsiz)
      call umalloc (ip_jendgrp,ngroup,1,iname,icumsiz)
      call umalloc (ip_jskpgrp,ngroup,1,iname,icumsiz)
      call umalloc (ip_kbeggrp,ngroup,1,iname,icumsiz)
      call umalloc (ip_kendgrp,ngroup,1,iname,icumsiz)
      call umalloc (ip_kskpgrp,ngroup,1,iname,icumsiz)
      call umalloc (ip_iblkgrp,ngroup,1,iname,icumsiz)
c
      call umalloc (ip_imap,msegt*msegn*mbloc*mgrlev,1,iname,icumsiz)
      call umalloc (ip_nseg,mbloc,1,iname,icumsiz)
c
      call umalloc (ip_wcount,mres ,0,iname,icumsiz)
      call umalloc (ip_rtrmsg,mres ,0,iname,icumsiz)
      call umalloc (ip_hrmsg ,mres ,0,iname,icumsiz)
      call umalloc (ip_nsupg ,mres ,1,iname,icumsiz)
      call umalloc (ip_resout,mres ,0,iname,icumsiz)
      call umalloc (ip_supout,mres ,0,iname,icumsiz)
      call umalloc (ip_clg   ,mres ,0,iname,icumsiz)
      call umalloc (ip_cmxg  ,mres ,0,iname,icumsiz)
      call umalloc (ip_cmyg  ,mres ,0,iname,icumsiz)
      call umalloc (ip_cmzg  ,mres ,0,iname,icumsiz)
      call umalloc (ip_cdtg  ,mres ,0,iname,icumsiz)
      call umalloc (ip_cdvg  ,mres ,0,iname,icumsiz)
      call umalloc (ip_cxg   ,mres ,0,iname,icumsiz)
      call umalloc (ip_cyg   ,mres ,0,iname,icumsiz)
      call umalloc (ip_czg   ,mres ,0,iname,icumsiz)
      call umalloc (ip_cng   ,mres ,0,iname,icumsiz)
      call umalloc (ip_cxvg  ,mres ,0,iname,icumsiz)
      call umalloc (ip_cyvg  ,mres ,0,iname,icumsiz)
      call umalloc (ip_czvg  ,mres ,0,iname,icumsiz)
      call umalloc (ip_cmx   ,mbloc,0,iname,icumsiz)
      call umalloc (ip_cmy   ,mbloc,0,iname,icumsiz)
      call umalloc (ip_cmz   ,mbloc,0,iname,icumsiz)
      call umalloc (ip_cdt   ,mbloc,0,iname,icumsiz)
      call umalloc (ip_cdv   ,mbloc,0,iname,icumsiz)
      call umalloc (ip_cx    ,mbloc,0,iname,icumsiz)
      call umalloc (ip_cy    ,mbloc,0,iname,icumsiz)
      call umalloc (ip_cz    ,mbloc,0,iname,icumsiz)
      call umalloc (ip_cl    ,mbloc,0,iname,icumsiz)
      call umalloc (ip_cxv   ,mbloc,0,iname,icumsiz)
      call umalloc (ip_cyv   ,mbloc,0,iname,icumsiz)
      call umalloc (ip_czv   ,mbloc,0,iname,icumsiz)
c
      call umalloc (ip_ntpchcb,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m1pch1 ,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_m1pch2 ,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_npchcbf,mgrlev*mbloc*6,1,iname,icumsiz)
      call umalloc (ip_ipatchc ,mxtpchc,1,iname,icumsiz)
      call umalloc (ip_jpatchc ,mxtpchc,1,iname,icumsiz)
      call umalloc (ip_kpatchc ,mxtpchc,1,iname,icumsiz)
c
      call umalloc (ip_nswpchb,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_lswpchb,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_lswpche,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (ip_npchitm,mgrlev*mxpchs,1,iname,icumsiz)
      call umalloc (ip_lspchb1,mgrlev*mxpchs,1,iname,icumsiz)
      call umalloc (ip_lspchb2,mgrlev*mxpchs,1,iname,icumsiz)
      call umalloc (ip_lspchf1,mgrlev*mxpchs,1,iname,icumsiz)
      call umalloc (ip_lspchf2,mgrlev*mxpchs,1,iname,icumsiz)
      call umalloc (ip_iitmsa ,mgrlev*mxpchs,1,iname,icumsiz)
      call umalloc (ip_ipitmb1,mpchitm,1,iname,icumsiz)
      call umalloc (ip_ipitmb2,mpchitm,1,iname,icumsiz)
      call umalloc (ip_jpitmb1,mpchitm,1,iname,icumsiz)
      call umalloc (ip_jpitmb2,mpchitm,1,iname,icumsiz)
      call umalloc (ip_kpitmb1,mpchitm,1,iname,icumsiz)
      call umalloc (ip_kpitmb2,mpchitm,1,iname,icumsiz)
      call umalloc (ip_ipitmbs,mpchitm,1,iname,icumsiz)
      call umalloc (ip_jpitmbs,mpchitm,1,iname,icumsiz)
      call umalloc (ip_kpitmbs,mpchitm,1,iname,icumsiz)
      call umalloc (ip_frc    ,mpchitm,0,iname,icumsiz)
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)

      call umalloc (ip_intval,mxintln,1,iname,icumsiz)
      call umalloc (ip_nodeid,mbloc  ,1,iname,icumsiz)

      call umalloc (ip_buff ,mxsiz3n,0,iname,icumsiz)
      call umalloc (ip_buffw,mxszfg*5,0,iname,icumsiz)
      call umalloc (ip_buffp,mxszfg ,0,iname,icumsiz)
      call umalloc (ip_buffe,mxszfg ,0,iname,icumsiz)
      call umalloc (ip_buff1,mxszfg ,0,iname,icumsiz)
      call umalloc (ip_buff2,mxszfg ,0,iname,icumsiz)
      call umalloc (ip_buffs,mxszfg ,0,iname,icumsiz)

      call umalloc (ip_tbuf1 ,mxpchc ,1,iname,icumsiz)
      call umalloc (ip_tbuf2 ,mxpchc ,1,iname,icumsiz)
      call umalloc (ip_tbuf3 ,mxpchc ,1,iname,icumsiz)
      call umalloc (ip_tbuf4 ,mxpitm ,1,iname,icumsiz)
      call umalloc (ip_tbuf5 ,mxpitm ,1,iname,icumsiz)
      call umalloc (ip_tbuf6 ,mxpitm ,1,iname,icumsiz)
      call umalloc (ip_tbuf7 ,mxpitm ,1,iname,icumsiz)
      call umalloc (ip_tbuf8 ,mxpitm ,1,iname,icumsiz)
      call umalloc (ip_tbuf9 ,mxpitm ,1,iname,icumsiz)
      call umalloc (ip_tbuf10,mxpitm ,0,iname,icumsiz)

      call umalloc (ip_buffwo,mx3dwk5*ntorder+1 ,0,iname, icumsiz)
      call umalloc (ip_bufftv1,mx3dwk*ntorder*iturv1+1 ,0,iname,
     .                                                    icumsiz)
      call umalloc (ip_bufftv2,mx3dwk*ntorder*iturv2+1 ,0,iname,
     .                                                    icumsiz)

#endif
c
      call umalloc (ip_ivisb,mbloc*3,1,iname,icumsiz)
      call umalloc (ip_const1,mbloc*msegn,0,iname,icumsiz)
      call umalloc (ip_const2,mbloc*msegn,0,iname,icumsiz)
      call umalloc (ip_const3,mbloc*msegn,0,iname,icumsiz)
      call umalloc (ip_const4,mbloc*msegn,0,iname,icumsiz)
      call umalloc (ip_const5,mbloc*msegn,0,iname,icumsiz)
c
      call umalloc (ip_w    ,mxsiz5c,0,iname,icumsiz)
      call umalloc (ip_p    ,mxsizc ,0,iname,icumsiz)
      call umalloc (ip_vol  ,mxsizc ,0,iname,icumsiz)
      call umalloc (ip_eomu ,mxsizc ,0,iname,icumsiz)
      call umalloc (ip_smin ,mxsizc*iturv1+1 ,0,iname,icumsiz)
      call umalloc (ip_turv1,mxsizc*iturv1+1 ,0,iname,icumsiz)
      call umalloc (ip_turv2,mxsizc*iturv2+1 ,0,iname,icumsiz)
c
      call umalloc (ip_wold   ,mxsz5fg*ntorder+1 ,0,iname,icumsiz)
      call umalloc (ip_tv1old,mxszfg*ntorder*iturv1+1 ,0,iname,icumsiz)
      call umalloc (ip_tv2old,mxszfg*ntorder*iturv2+1 ,0,iname,icumsiz)
c
      n_mebdf = mebdf*(ntorder-1)
      call umalloc (ip_wbar  ,mxsz5fg*n_mebdf+1 ,0,iname,icumsiz)
      call umalloc (ip_tv1bar,mxszfg*n_mebdf*iturv1+1 ,0,iname,icumsiz)
      call umalloc (ip_tv2bar,mxszfg*n_mebdf*iturv2+1 ,0,iname,icumsiz)
c
      call umalloc (ip_wavg   ,mxsz5fg ,0,iname,icumsiz)
      call umalloc (ip_pavg   ,mxszfg  ,0,iname,icumsiz)
      call umalloc (ip_eomuavg,mxszfg  ,0,iname,icumsiz)
c
c#ifdef  RKPHY
      call umalloc (ip_delwork,mxsz5fg*nrkstg+1 ,0,iname,icumsiz)
      call umalloc (ip_dtv1ork,mxszfg*nrkstg*iturv1+1 ,0,iname,icumsiz)
      call umalloc (ip_dtv2ork,mxszfg*nrkstg*iturv2+1 ,0,iname,icumsiz)
      call umalloc (ip_rkphysa,mrkstg*mrkstg ,0,iname,icumsiz)
      call umalloc (ip_rkphyse,mrkstg        ,0,iname,icumsiz)
c#endif
c
      call umalloc (ip_x    ,mxsiz3n ,0,iname,icumsiz)
      call umalloc (ip_six  ,mxsizn  ,0,iname,icumsiz)
      call umalloc (ip_siy  ,mxsizn  ,0,iname,icumsiz)
      call umalloc (ip_siz  ,mxsizn  ,0,iname,icumsiz)
      call umalloc (ip_sjx  ,mxsizn  ,0,iname,icumsiz)
      call umalloc (ip_sjy  ,mxsizn  ,0,iname,icumsiz)
      call umalloc (ip_sjz  ,mxsizn  ,0,iname,icumsiz)
      call umalloc (ip_skx  ,mxsizn  ,0,iname,icumsiz)
      call umalloc (ip_sky  ,mxsizn  ,0,iname,icumsiz)
      call umalloc (ip_skz  ,mxsizn  ,0,iname,icumsiz)
c
      call umalloc (ip_w1   ,mxsz5cg ,0,iname,icumsiz)
      call umalloc (ip_wr   ,mxsz5cg ,0,iname,icumsiz)
      call umalloc (ip_ws   ,mxsz5fg ,0,iname,icumsiz)
      call umalloc (ip_fw   ,mxsz5fg ,0,iname,icumsiz)
      call umalloc (ip_fv   ,mxsz4fg ,0,iname,icumsiz)
      call umalloc (ip_dti  ,mxszfg  ,0,iname,icumsiz)
      call umalloc (ip_dtj  ,mxszfg  ,0,iname,icumsiz)
      call umalloc (ip_dtk  ,mxszfg  ,0,iname,icumsiz)
      call umalloc (ip_dtvi ,mxsz5fg ,0,iname,icumsiz)
      call umalloc (ip_dtvj ,mxsz5fg ,0,iname,icumsiz)
      call umalloc (ip_dtvk ,mxsz5fg ,0,iname,icumsiz)
      call umalloc (ip_fbeta2 ,mxszfg ,0,iname,icumsiz)
      call umalloc (ip_fbetav2,mxszfg ,0,iname,icumsiz)
c
      call umalloc (ip_delw  ,mx3dwk*5 ,0,iname,icumsiz)
      call umalloc (ip_dtl   ,mx3dwk ,0,iname,icumsiz)
      call umalloc (ip_eprs  ,mx3dwk3,0,iname,icumsiz)
      call umalloc (ip_ratioij ,mx3dwk,0,iname,icumsiz)
      call umalloc (ip_ratiojk ,mx3dwk,0,iname,icumsiz)
c
      call umalloc (ip_ri1   ,mxs2d3c,0,iname,icumsiz)
      call umalloc (ip_ri2   ,mxs2d3c,0,iname,icumsiz)
      call umalloc (ip_rj1   ,mxs2d3c,0,iname,icumsiz)
      call umalloc (ip_rj2   ,mxs2d3c,0,iname,icumsiz)
      call umalloc (ip_rk1   ,mxs2d3c,0,iname,icumsiz)
      call umalloc (ip_rk2   ,mxs2d3c,0,iname,icumsiz)
c
      call umalloc (ip_tauwfi1 ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_tauwfi2 ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_tauwfj1 ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_tauwfj2 ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_tauwfk1 ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_tauwfk2 ,mxs2dc ,0,iname,icumsiz)
c
      call umalloc (ip_uswfi1  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_uswfi2  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_uswfj1  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_uswfj2  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_uswfk1  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_uswfk2  ,mxs2dc ,0,iname,icumsiz)
c
      call umalloc (ip_uplwfi1  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_uplwfi2  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_uplwfj1  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_uplwfj2  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_uplwfk1  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_uplwfk2  ,mxs2dc ,0,iname,icumsiz)
c
      call umalloc (ip_tnuwfi1  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_tnuwfi2  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_tnuwfj1  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_tnuwfj2  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_tnuwfk1  ,mxs2dc ,0,iname,icumsiz)
      call umalloc (ip_tnuwfk2  ,mxs2dc ,0,iname,icumsiz)
c
      call umalloc (ip_fwrk  ,mxdstf ,0,iname,icumsiz)
      call umalloc (ip_iwrk  ,mxdsti ,1,iname,icumsiz)
c
      call umalloc (ip_wk1d  ,mx1dwk*27 ,0,iname,icumsiz)
      call umalloc (ip_wk2d  ,mx2dwk*15 ,0,iname,icumsiz)
      call umalloc (ip_wk1d5 ,mx1dwk5*5 ,0,iname,icumsiz)
      call umalloc (ip_fqs   ,mx1dwk*mx1dwk*6   ,0,iname,icumsiz)
      call umalloc (ip_fpval ,mxinpln*8   ,0,iname,icumsiz)
c
      if (myrank.eq.master) then
       memtotal  = icumsiz/1.e+06
       write (iwrit,'("memory (MB) needed in tlns3d =",i8)') memtotal
      endif
c
      totptim = 0.
      totptse = 0
      do n=1,5
        erkmaxl(n) = 0.
        erkmaxg(n) = 0.
        erkrms2(n) = 0.
        erkrmsg(n) = 0.
      enddo
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
#  ifdef BUILD_PVM
c
c---- PVM initialization
c
      write (iwrit,'("PVM not supported")')
      stop
#  else
c
c---- MPI initialization
c
c     call xgetpid (mypid)
c     write (iwrit, 6105) myrank, mypid
c6105  format ("rank ", i2, ": mypid = ", i8)
c     call flush (iwrit)
c     call sleep (300)
c
c     kill any extra processes
c
      if (nnodes.gt.mbloc) then
        do inode=mbloc+1,nnodes
          call MPI_Send (ndum,0,MPI_INTEGER,
cjm     .                   inode,TAG_ERROR,
     .                   inode-1,TAG_ERROR,
     .                   mycomm,ierr)
        enddo
        nnodes = nbloc
      endif
c
      do inode=1,nnodes
        nodeid(inode)       = inode
        nodes(ndlist+inode) = inode
      enddo
#  endif
#else
c
c---- IN-CORE initialization
c
      do ibloc=1,mbloc
        nodes(ibloc)  = 0
      enddo
      nodes(ndlist+1) = 0
#endif
c
c---- initialize for timings
c
c     timing array modifed to track user and system time
c     tim(1,1)  = total user time
c     tim(2,1)  = total system time
c     tim(3,1)  = total wall clock time
c     tim(1,2)  = user time since last call to cputim()
c     tim(2,2)  = system time since last call to cputim()
c     tim(3,2)  = wall clock time since last call to cputim()
c     tim(1,3),tim(2,3) and tim(3,3) used for intermediate results
c
      do j=1,3
      do i=1,3
        tim(i,j) = 0.
      enddo
      enddo
      call cputim (tim)
      rad       = 45./atan(1.)
c
c         coefficients of leading term for BDF
c         index is ntorder+1
c
      ct(1) = 0.
      ct(2) = 1.
      ct(3) = 1.5
      ct(4) = 11./6.
      ct(5) = 25./12.
      ct(6) = 137./60.
      ct(7) = 49./20
c            tolerance to stop MG cycles
      resend    = 1.e-17
c
c     mod for bdf2opt
c
c     if(ibdf2opt.eq.1) ct(4) = 1.69
      if(ibdf2opt.eq.1) ct(4) = 1.5d0 + 0.58d0/3.0d0
c
c####################################################################
c master process (rank 0)
c     read input and send to nodes
c####################################################################
c
      if (myrank.eq.master) then
c
      rewind iread
c
c     read the standard TLNS3D input. The data is read into
c     an array for simplicity since it is to be sent to all
c     the processes. The array is set up for 8 input values
c     per line.
c
      do j=1,mxinpln
        do i=1,8
          fpval(i,j) = 0.
        enddo
      enddo
c
      read (iread,'(a80)') title
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') grdfmt,fstrt,fsave,fplot,fcplpr,fnout
c                            ,favg
      read (iread,'(8f10.6)') (fpval(n,1),n=1,7)
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') flev,fgrids,fmgtyp,ftorder,dtphy,
c                             frkcase,ftolrk
      read (iread,'(8f10.6)') (fpval(n,2),n=1,7)
c
c     nlev   = number of fmg levels used
c     ngrid  = no. of grids to be used in multigrid
c              for finest fmg level
c
      nlev      = fpval(1,2)
      if (nlev .le. 0) nlev = 1
      ngrid     = fpval(2,2)
      if (nlev .gt. mgrlev) then
        write (iwrit,'(" fmg levels = ",i3,2x,"specified for this",
     .       1x,"run exceeds the maximum dimension in code = ",1x,i3)' )
     .           nlev,mgrlev
        call ERREXIT (nodes)
      endif
      if (ngrid .gt. mgrlev) then
        write (iwrit,'(" no. grids = ",i3,2x,"specified for this",
     .       1x,"run exceeds the maximum dimension in code = ",1x,i3)' )
     .           ngrid,mgrlev
        call ERREXIT (nodes)
      endif
c
c     currently number grids limited to 8 just to make input easier
c
      if (ngrid .gt. 8) then
        write (iwrit,'(" no. grids = ",i3," greater than 8")') ngrid
        call ERREXIT (nodes)
      endif
c
      read (iread,'(1x)')
      do ilev=1,nlev
c       read (iread,'(8f10.6)') fcyc,fbcupd,fcnvout,fsteps,ftol
        read (iread,'(8f10.6)') (fpval(n,2+ilev),n=1,5)
      enddo
c
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') (gitr(igrid),igrid=1,ngrid)
      read (iread,'(8f10.6)') (fpval(n,3+nlev),n=1,ngrid)
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') (gitp(igrid),igrid=1,ngrid)
      read (iread,'(8f10.6)') (fpval(n,4+nlev),n=1,ngrid)
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') cflf,hmf,vis2,vis4,zeta
      read (iread,'(8f10.6)') (fpval(n,5+nlev),n=1,5)
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') c(1),c(2),c(3),c(4),c(5),c(6)
      read (iread,'(8f10.6)') (fpval(n,6+nlev),n=1,6)
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') qfil(1),qfil(2),qfil(3),
c    .                        qfil(4),qfil(5),qfil(6)
      read (iread,'(8f10.6)') (fpval(n,7+nlev),n=1,6)
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') beta(1),beta(2),beta(3),
c    .                        beta(4),beta(5),beta(6)
      read (iread,'(8f10.6)') (fpval(n,8+nlev),n=1,6)
c
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') smoopi,smoopj,smoopk,vepsn,vepsl,cau,
c    .       enteps
      read (iread,'(8f10.6)') (fpval(n,9+nlev),n=1,7)
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') cflc,hmc,vis0,smoopic,smoopjc,smoopkc
      read (iread,'(8f10.6)') (fpval(n,10+nlev),n=1,6)
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') rm,reyl,al,yaw,roll,fi2dfl
      read (iread,'(8f10.6)') (fpval(n,11+nlev),n=1,6)
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') gamma,pr,prt,suthc,tref,fiturb,fpatchg
      read (iread,'(8f10.6)') (fpval(n,12+nlev),n=1,7)
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') xref,yref,zref,sref
      read (iread,'(8f10.6)') (fpval(n,13+nlev),n=1,4)
c     read (iread,'(1x)')
c     read (iread,'(8f10.6)') feng,epsf,ptj,ttj,fjet,ptjopfs,ttjott0,
c    .                         rmjet
cprec
c     read preconditioning related parameters now
c
c     iprec = 0; No preconditioning, dissipation use original variables
c             1; preconditioning, dissipation in p,u,v,w,T variables
c             2; preconditioning, dissipation in conservation variables
c             3; Non-conservative preconditioning
c     iprecg controls preconditioning in iteration procedure
c                (not artificial viscosity)
c
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') fiprec,cfls,ubmin,upk,ulprec,
c    .                         deltpre,ubcon
c---  ubmin is modified by free-stream velocity later on ----
      read (iread,'(8f10.6)') (fpval(n,14+nlev),n=1,7)
c
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') diffac,rsvis,omega,betvis,betvisv,
c    .                        npcons,nexp
      read (iread,'(8f10.6)') (fpval(n,15+nlev),n=1,7)
c
      read (iread,'(1x)')
c     read (iread,'(8f10.6)') anin,aninb,bninb,ijac,icord,inav,iprbet
      read (iread,'(8f10.6)') (fpval(n,16+nlev),n=1,7)
c
c---- read the file names
c
      read  (iread,'(a80)') grdfil
      call fixstr (grdfil,80)
      read  (iread,'(a80)') mapfil
      call fixstr (mapfil,80)
      read  (iread,'(a80)') rstfil
      call fixstr (rstfil,80)
      read  (iread,'(a80)') pltfil
      call fixstr (pltfil,80)
      read  (iread,'(a80)') savfil
      call fixstr (savfil,80)
      read  (iread,'(a80)') pchfil1
      call fixstr (pchfil1,80)
      read  (iread,'(a80)') pchfil2
      call fixstr (pchfil2,80)
      print *,'files=',grdfil,mapfil,rstfil,pltfil,
     .  savfil,pchfil1,pchfil2
c
c     read Namelist file with auxiliary data input
c
c     iauxfil = 55
c     open (iauxfil,form='formatted')
c     rewind iauxfil
c     read (iauxfil,auxinput)
c     write(*, auxinput)
c
c     The parameters in this namelist are now read through standard input file
c
      diffac  = fpval(1,15+nlev)
      rsvis   = fpval(2,15+nlev)
      omega   = fpval(3,15+nlev)
      betvis  = fpval(4,15+nlev)
      betvisv = fpval(5,15+nlev)
      fnpcons = fpval(6,15+nlev)
      fnexp   = fpval(7,15+nlev)
c
      anin    = fpval(1,16+nlev)
      aninb   = fpval(2,16+nlev)
      bninb   = fpval(3,16+nlev)
      fijac   = fpval(4,16+nlev)
      ficord  = fpval(5,16+nlev)
      finav   = fpval(6,16+nlev)
      fiprbet = fpval(7,16+nlev)
c     re-interpret some input variables as integers
      npcons  = fnpcons
      nexp    = fnexp
      ijac    = fijac
      icord   = ficord
      inav    = finav
      iprbet  = fiprbet
c
c
c        iprec=0  no preconditioning - conservation variables
c        iprec=1  AV and update based on conservation variables
c        iprec=2  AV based on primitive variables
c                 update depends on npcons
c
c        iprecg<0   no preconditioning in update stage only in AV
c
c        nexp relevant only for dual time stepping
c        nexp=1 explicit treatment of physical time dependent RHS
c        nexp=0 implicit treatment of physical time dependent RHS
c        nexp is used only in advflow with iprec=2
c
c        if npcons=0  no implicit correction term of time accurate derivative
c        npcons=0 always explicit
c        if iprec=1 and npcons=1,2  update conservation variables
c        iprec=1 npcons=1 always implicit
c        if iprec=2 and npcons=1  update conservation variables using
c           linearization from (p,u,v,w,T) variables
c           explicit/implicit depends on nexp
c        if iprec=2 and npcons=2 update primitive variables
c        iprec=2 npcons=2 always implicit
c
c
c           anin  is factor for inclusion of dtphy for time step
c           aninb is factor for inclusion of dtphy for beta
c           bninb is factor for inclusion of dtphy for betav
      if (ntorder.eq.0) then
        anin     = 0.0
        aninb    = 0.0
        bninb    = 0.0
      endif
c     write(*, auxinput)
c
#ifdef BUILD_MPI
          call MPI_Bcast (npcons,1,MPI_INTEGER,master,mycomm,ierr)
          call MPI_Bcast (nexp,1,MPI_INTEGER,master,mycomm,ierr)
          call MPI_Bcast (iprbet,1,MPI_INTEGER,master,mycomm,ierr)
          call MPI_Bcast (ijac,1,MPI_INTEGER,master,mycomm,ierr)
          call MPI_Bcast (icord,1,MPI_INTEGER,master,mycomm,ierr)
          call MPI_Bcast (inav,1,MPI_INTEGER,master,mycomm,ierr)
          call MPI_Bcast (diffac,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (rsvis,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (omega,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (anin,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (aninb,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (bninb,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (betvis,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (betvisv,1,RTYPE,master,
     .                    mycomm,ierr)
#endif
c
c
c---- get number blocks and block sizes from grid file
c
      len     = istlen (grdfil,80)
      igrdfmt = fpval(1,1)
      i2dfl   = fpval(6,11+nlev)
c       i2dfl = 0 for 3d case
c       i2dfl = 1 for 2d case
      if (igrdfmt.eq.1) then
        open (igrdfil,file=grdfil(1:len),form='formatted')
        rewind igrdfil
        read (igrdfil,*,end=6000) nbloc
        if(i2dfl.eq.1) then
          read (igrdfil,*) (im(1,ibloc),jm(1,ibloc),ibloc=1,nbloc)
        else
          read (igrdfil,*) (im(1,ibloc),jm(1,ibloc),
     .                      km(1,ibloc),ibloc=1,nbloc)
        endif
      else
        open (igrdfil,file=grdfil(1:len),form='unformatted') 
        rewind igrdfil
        read (igrdfil,end=6000) nbloc
        if(i2dfl.eq.1) then
          read (igrdfil) (im(1,ibloc),jm(1,ibloc),ibloc=1,nbloc)
        else
          read (igrdfil) (im(1,ibloc),jm(1,ibloc),
     .                    km(1,ibloc),ibloc=1,nbloc)
        endif
      endif
      go to 6001
 6000 write (iwrit,'("grid file not found or empty")')
      call ERREXIT (nodes)
c
c     check the dimensions regarding  # of blocks
c
 6001 if (nbloc.gt.mbloc) then
       write (iwrit,*)
       write (iwrit,'(2x,"the no. of blocks input for this run = ",
     .                i2)' ) nbloc
       write (iwrit,'(2x,"exceeds the max block size, which is = ",
     .                i2)' ) mbloc
       call ERREXIT (nodes)
      endif
c
c     set k dimension for 2-D case
c
      if (i2dfl.eq.1) then
        do ibloc=1,nbloc
          km(1,ibloc) = 2
        enddo
      endif
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c---- set up block to node process mapping
c
      call mapblks (mgrlev,mbloc,nbloc,im,jm,km,
     .              nnodes,nodes,intval,65,iwrit,ierr)
      if (ierr.ne.0) call ERREXIT (nodes)
c
c     build the list of nodes that are being used
c     and kill any unused nodes
c
      do inode=1,nnodes
        intval(inode) = 0
      enddo
      do ibloc=1,nbloc
        inode         = nodes(ibloc)
        nodes(ibloc)  = nodeid(inode)
        intval(inode) = intval(inode) + 1
      enddo
      icnt = 0
      do inode=1,nnodes
        if (intval(inode).eq.0) then
#  ifdef BUILD_PVM
          call PVMFpsend (nodeid(inode),TAG_ERROR,
     .                    ndum,1,INTEGER4,ierr)
#  else
cjm
         if (nodeid(inode)-1.ne.master) then
          call MPI_Send (ndum,1,MPI_INTEGER,
cjm     .                   nodeid(inode),TAG_ERROR,
     .                   nodeid(inode)-1,TAG_ERROR,
     .                   mycomm,ierr)
cjm
         end if
#  endif
        else
          icnt               = icnt + 1
          nodes(ndlist+icnt) = nodeid(inode)
        endif
      enddo
      nnodes = icnt
c
c---- send input data to nodes and wait for OK response
c     the last 8 values in the fpval array are used to
c     send master process parameter values
c
      lines = 15 +nlev
      fpval(1,lines) = nbloc
      fpval(2,lines) = mbloc
      fpval(3,lines) = mgrlev
      fpval(4,lines) = msegt
      fpval(5,lines) = msegn
      fpval(6,lines) = mxtpchc
      fpval(7,lines) = mxpchs
      fpval(8,lines) = mpchitm

      do inode=1,nnodes
#  ifdef BUILD_PVM
        call PVMFpsend (nodes(ndlist+inode),TAG_INPUT,
     .                  fpval,lines*8,RTYPE,ierr)
        call PVMFprecv (nodes(ndlist+inode),-1,
     .                  errmsg,128,BYTE1,
     .                  itid,itag,ilen,ierr)
#  else
cjm
       if (nodeid(inode)-1.ne.master) then
        call MPI_Send (fpval,lines*8,RTYPE,
cjm     .                 nodes(ndlist+inode),TAG_INPUT,
     .                 nodes(ndlist+inode)-1,TAG_INPUT,
     .                 mycomm,ierr)
        call MPI_Recv (errmsg,128,MPI_CHARACTER,
cjm     .                 nodes(ndlist+inode),MPI_ANY_TAG,
     .                 nodes(ndlist+inode)-1,MPI_ANY_TAG,
     .                 mycomm,istat,ierr)
        itag = istat(MPI_TAG)
c
       end if
#  endif
cjm
      if (nodeid(inode)-1.ne.master) then
        if (itag.eq.TAG_ERROR) then
          len = istlen (errmsg,128)
          write (iwrit,'(/," error from node",i4,": ",a)')
     .      inode,errmsg(1:len)
          call ERREXIT (nodes)
        endif
cjm
      end if
      enddo
c
c#######################################################################
cjm
c node processes (for PVM)
c all processes (for MPI)
c     receive input from master process
c#######################################################################
c
#  ifdef BUILD_PVM
      else
#  else
      end if
#  endif
c
#  ifdef BUILD_PVM
      call PVMFprecv (master,-1,fpval,mxinpln*8,RTYPE,
     .               itid,itag,ilen,ierr)
#  else
      if (myrank.ne.master) then
      call MPI_Recv (fpval,mxinpln*8,RTYPE,master,
     .               MPI_ANY_TAG,mycomm,istat,ierr)
      itag = istat(MPI_TAG)
      end if
      if (myrank.ne.master) then
          call MPI_Bcast (npcons,1,MPI_INTEGER,master,mycomm,ierr)
          call MPI_Bcast (nexp,1,MPI_INTEGER,master,mycomm,ierr)
          call MPI_Bcast (iprbet,1,MPI_INTEGER,master,mycomm,ierr)
          call MPI_Bcast (ijac,1,MPI_INTEGER,master,mycomm,ierr)
          call MPI_Bcast (icord,1,MPI_INTEGER,master,mycomm,ierr)
          call MPI_Bcast (inav,1,MPI_INTEGER,master,mycomm,ierr)
          call MPI_Bcast (diffac,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (rsvis,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (omega,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (anin,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (aninb,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (bninb,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (betvis,1,RTYPE,master,
     .                    mycomm,ierr)
          call MPI_Bcast (betvisv,1,RTYPE,master,
     .                    mycomm,ierr)
      end if
#  endif
cjm
      if (myrank.ne.master) then
      if (itag.ne.TAG_INPUT)
     .  call ERREXIT (nodes)
cjm
      end if
      nlev      = fpval(1,2)
      if (nlev .le. 0) nlev = 1
      ngrid     = fpval(2,2)
      ipatchg   = fpval(7,12+nlev)
c
c---- check parameter values in case different than master
c     the following requirements could be relaxed if the data
c     was packed before being sent, rather than sent intact
c
      lines   = 15 +nlev
      ibloc   = fpval(2,lines) + .1
      igrlev  = fpval(3,lines) + .1
      isegt   = fpval(4,lines) + .1
      isegn   = fpval(5,lines) + .1
      ixtpchc = fpval(6,lines) + .1
      ixpchs  = fpval(7,lines) + .1
      ipchitm = fpval(8,lines) + .1
      itag    = TAG_OK
c
      if (nlev.gt.mgrlev .or. ngrid.gt.mgrlev) then
        itag = TAG_ERROR
        write (errmsg,'("fmg levels or no. grids exceeds the ",
     .    "maximum dimension in code")')
      endif
c
c     need msegt,msegn and mbloc the same for imap array
c
      if (isegt.ne.msegt .or. isegn.ne.msegn .or. ibloc.ne.mbloc) then
        itag = TAG_ERROR
        write (errmsg,'("msegt, msegn or mbloc parameter ",
     .    "not equal to the master process values")')
      endif
c
c     need mxtpchc,mxpchs and mpchitm the same for patched grid data
c
      if (ipatchg.gt.0 .and. (ixtpchc.ne.mxtpchc .or.
     .    ixpchs.ne.mxpchs .or. ipchitm.ne.mpchitm)) then
        itag = TAG_ERROR
        write (errmsg,'("mxtpchc, mxpchs or mpchitm parameter ",
     .    "not equal to the master process values")')
      endif
c
c     send error message or OK tag to master
c
#  ifdef BUILD_PVM
      call PVMFpsend (master,itag,errmsg,128,BYTE1,ierr)
#  else
cjm
      if (myrank.ne.master) then
      call MPI_Send (errmsg,128,MPI_CHARACTER,master,
     .               itag,mycomm,ierr)
cjm
      end if
#  endif
c
      if (itag.eq.TAG_ERROR) then
#  ifdef BUILD_PVM
        call PVMFprecv (master,TAG_ERROR,ndum,1,INTEGER4,
     .                  itid,itag,ilen,ierr)
#  else
cjm
      if (myrank.ne.master) then
        call MPI_Recv (ndum,1,MPI_INTEGER,master,
     .                 TAG_ERROR,mycomm,istat,ierr)
cjm
      end if
#  endif
cjm
      if (myrank.ne.master) then
        call ERREXIT (nodes)
cjm
      end if
      endif
c
c     if fnout > 0, then open a file for output on the node.
c     Since we need to account for cross-mounted file systems, we name
c     the file: node.myrank
c
      if (fpval(6,1).gt.0.) then
        if (myrank.gt.99) then
          len = 8
          write (title,'("node.",i3)') myrank
        else if (myrank.gt.9) then
          len = 7
          write (title,'("node.",i2)') myrank
        else
          len = 6
          write (title,'("node.",i1)') myrank
        endif
        do n=len+1,80
          title(n:n) = ' '
        enddo
        open (iwrit,file=title(1:len))
      else
       if (myrank.ne.master) then
        iwrit = 0
       end if
      endif
c
#endif
c
#  ifdef BUILD_MPI
#  else
      endif
#  endif
c
c#######################################################################
c all processes
c     variable assignment and data initialization
c     This makes sure all processes are using the same data
c#######################################################################
c
      grdfmt  = fpval(1,1)
      fstrt   = fpval(2,1)
      fsave   = fpval(3,1)
      fplot   = fpval(4,1)
      fcplpr  = fpval(5,1)
      fnout   = fpval(6,1)
      favg    = fpval(7,1)
c
      flev    = fpval(1,2)
      fgrids  = fpval(2,2)
      fmgtyp  = fpval(3,2)
      ftorder = fpval(4,2)
      dtphy   = fpval(5,2)
      frkcase = fpval(6,2)
      irkcase = frkcase
      ftolrk  = fpval(7,2)
c
c     initialise and read coefficients for implicit
c     RK multi stage scheme for physical time stepping
c
#ifdef  RKPHY
      if(irkcase.gt.0) then
         do jj=1,mrkstg
            rkphyse(jj)=0.
            do ii=1,mrkstg
               rkphysa(ii,jj) = 0.
            enddo
         enddo
c
         call rkphfil (irkcase,mrkstg,rkphysa,rkphyse,sorder)
      endif
#endif
c
      do ilev=1,nlev
        mcyc(ilev)    = fpval(1,2+ilev)
        mbcupd(ilev)  = fpval(2,2+ilev)
        mcnvout(ilev) = fpval(3,2+ilev)
        nsteps(ilev)  = fpval(4,2+ilev)
        if(nsteps(ilev).le.0) nsteps(ilev) = 1
        ftol          = fpval(5,2+ilev)
      enddo
c
      do igrid=1,ngrid
        gitr(igrid) = fpval(igrid,3+nlev)
        gitp(igrid) = fpval(igrid,4+nlev)
      enddo
c
      cflf    = fpval(1,5+nlev)
      hmf     = fpval(2,5+nlev)
      vis2    = fpval(3,5+nlev)
      vis4    = fpval(4,5+nlev)
      zeta    = fpval(5,5+nlev)
c
      do n=1,6
        c(n)    = fpval(n,6+nlev)
        qfil(n) = fpval(n,7+nlev)
        beta(n) = fpval(n,8+nlev)
      enddo
c
      smoopi  = fpval(1,9+nlev)
      smoopj  = fpval(2,9+nlev)
      smoopk  = fpval(3,9+nlev)
      vepsn   = fpval(4,9+nlev)
      vepsl   = fpval(5,9+nlev)
      cau     = fpval(6,9+nlev)
      enteps  = fpval(7,9+nlev)
c
      cflc    = fpval(1,10+nlev)
      hmc     = fpval(2,10+nlev)
      vis0    = fpval(3,10+nlev)
      smoopic = fpval(4,10+nlev)
      smoopjc = fpval(5,10+nlev)
      smoopkc = fpval(6,10+nlev)
c
      rm      = fpval(1,11+nlev)
      reyl    = fpval(2,11+nlev)
      al      = fpval(3,11+nlev)
      yaw     = fpval(4,11+nlev)
      roll    = fpval(5,11+nlev)
      fi2dfl  = fpval(6,11+nlev)
      amachg  = rm*rm
c
      gamma   = fpval(1,12+nlev)
      pr      = fpval(2,12+nlev)
      prt     = fpval(3,12+nlev)
      suthc   = fpval(4,12+nlev)
      tref    = fpval(5,12+nlev)
      fiturb  = fpval(6,12+nlev)
      fpatchg = fpval(7,12+nlev)
c
      xref    = fpval(1,13+nlev)
      yref    = fpval(2,13+nlev)
      zref    = fpval(3,13+nlev)
      sref    = fpval(4,13+nlev)
c
cprec
c
      fiprec  = fpval(1,14+nlev)
      cfls    = fpval(2,14+nlev)
      ubmin   = fpval(3,14+nlev)
      upk     = fpval(4,14+nlev)
      ulprec  = fpval(5,14+nlev)
      deltpre = fpval(6,14+nlev)
      ubcon   = fpval(7,14+nlev)
c
c---------------------------------------------------------------------
c     echo input to output
c---------------------------------------------------------------------
c
      if (iwrit.gt.0) then
c
      write (iwrit,'(/,
     .  " Program : TLNS3D - Multiblock Version")')
#ifdef BUILD_PVM
      write (iwrit,'(
     .  " Distributed Computing using PVM message passing")')
#endif
#ifdef BUILD_MPI
      write (iwrit,'(
     .  " Distributed Computing using MPI message passing")')
#endif
#ifdef TLN_DBL
      write (iwrit,'(" Double Precision")')
#endif

      write (iwrit,'(
     .  " Multistage Time-Stepping scheme using",
     .  " Multigrid acceleration technique for solving",/,
     .  " Thin-Layer Navier-Stokes Equations with",
     .  " Baldwin-Lomax and one and two-equation Turb. Models",/)')
c
      len = istlen (title,80)
      write (iwrit,'(1x,a)') title(1:len)
c
      write (iwrit,'(/," Echo the input file to flow code",/)')
c
      write (iwrit,'(5x,"  grdfmt      mstrt       msave       mplot  ",
     .                "    mcplpr      nodeout     favg      nrkstg" )')
      write (iwrit,'(1x,7f12.5,i7)') grdfmt,fstrt,fsave,fplot,fcplpr,
     .                               fnout,favg,nrkstg
      write (iwrit,'(5x,"   nlev       ngrids      fmgtyp",3x,
     .                 "   ntorder    dtphy       frkcase     ftolrk")')
      write (iwrit,'(1x,7f12.5)') flev,fgrids,fmgtyp,ftorder,dtphy,
     .                            frkcase,ftolrk
c
      write (iwrit,'(5x,"   ilev        ncyc       mbcupd  ",
     .                  "   mcnvout     nsteps     ftol  " )')
      do 12 ilev=1,nlev
c
c       fcyc    : no. of multigrid cycles at this level in fmg cycle
c       fbcupd  = 0 update b.c. on every stage of Runge-Kutta
c             .gt.0 update b.c. only on last stage of Runge-Kutta
c       fcnvout = 0 donot write residual convergence on plot file
c             .gt.0       write residual convergence on plot file
c
      fcyc    = mcyc(ilev)
      fbcupd  = mbcupd(ilev)
      fcnvout = mcnvout(ilev)
      fsteps  = nsteps(ilev)
      write (iwrit,'(6x,i5,2x,6f12.5)') ilev,fcyc,fbcupd,fcnvout,
     .                                  fsteps,ftol
  12  continue
c
      write (iwrit,'(5x," nitr(1)     nitr(2)     nitr(3)    ",
     .                  " nitr(4)     nitr(5)     nitr(6)"    )')
      write (iwrit,'(1x,7f12.5)') (gitr(igrid),igrid=1,ngrid)
      write (iwrit,'(5x," nitp(1)     nitp(2)     nitp(3)    ",
     .                  " nitp(4)     nitp(5)     nitp(6)    ")')
      write (iwrit,'(1x,7f12.5)') (gitp(igrid),igrid=1,ngrid)
      write (iwrit,'(5x,"  cflf         hmf        vis2      ",
     .                  "  vis4        zeta    ")')
      write (iwrit,'(1x,7f12.5)') cflf,hmf,vis2,vis4,zeta
      write (iwrit,'(5x,"   c(1)        c(2)        c(3)     ",
     .                  "   c(4)        c(5)        c(6)"    )')
      write (iwrit,'(1x,7f12.5)') c(1),c(2),c(3),c(4),c(5),c(6)
      write (iwrit,'(5x," qfil(1)     qfil(2)     qfil(3)    ",
     .                  " qfil(4)     qfil(5)     qfil(6)"    )')
      write (iwrit,'(1x,7f12.5)') qfil(1),qfil(2),qfil(3),
     .                            qfil(4),qfil(5),qfil(6)
      write (iwrit,'(5x," beta(1)     beta(2)     beta(3)    ",
     .                  " beta(4)     beta(5)     beta(6)"    )')
      write (iwrit,'(1x,7f12.5)') beta(1),beta(2),beta(3),
     .                            beta(4),beta(5),beta(6)
c
      write (iwrit,'(5x,"  smoopi      smoopj      smoopk  ",
     .                  "    vepsn       vepsl        cau        enteps"
     .                  )')
      write (iwrit,'(1x,7f12.5)') smoopi,smoopj,smoopk,vepsn,vepsl,cau,
     .       enteps
      write (iwrit,'(5x,"   cflc         hmc        vis0     ",
     .                  " smoopic     smoopjc     smoopkc"  )')
      write (iwrit,'(1x,7f12.5)') cflc,hmc,vis0,smoopic,smoopjc,smoopkc
      write (iwrit,'(5x,"    rm         reyl       alpha    ",
     .                  "    yaw         roll      i2dfl")')
      write (iwrit,'(1x,6f12.5)') rm,reyl,al,yaw,roll,fi2dfl
      write (iwrit,'(5x,"   gamma       pr          prt   ",
     .             "    suthc        tref       iturb       ipatchg")')
      write (iwrit,'(1x,7f12.5)') gamma,pr,prt,suthc,tref,fiturb,fpatchg
      write (iwrit,'(5x,"   xref        yref        zref    ",
     .                  "   sref"   )')
      write (iwrit,'(1x,7f12.5)') xref,yref,zref,sref
c
cprec
      write (iwrit,'(5x,"  iprec        cfls       ubmin     ",
     .                  "   upk       ulprec      deltpre      ubcon")')
      write (iwrit,'(1x,7f12.5)') fiprec,cfls,ubmin,upk,ulprec,
     .                            deltpre,ubcon
c
      write (iwrit,'(5x,"  diffac       rsvis      omega     ",
     .                  "  betvis     betvisv     npcons      nexp")')
      write (iwrit,'(1x,7f12.5)') diffac,rsvis,omega,betvis,betvisv,
     .                            fnpcons,fnexp
c
      write (iwrit,'(5x,"  anin         aninb       bninb    ",
     .                  " ijac        icord       inav        iprbet")')
      write (iwrit,'(1x,7f12.5)') anin,aninb,bninb,fijac,ficord,
     .                            finav,fiprbet
c
      endif
c
c---- determine the no. of runge-kutta stages from input
c
      mstage    = 1
      do 14 k=1,6
      if (c(k).ne.0.) mstage = k
   14 continue
      fstage    = mstage
c
c---- set up flags
c
c     igrdfmt: grid input format
c         igrdfmt.eq. 1  read formatted grid
c         igrdfmt.ne. 1  read unformatted grid
c     mstrt  : restart option
c         mstrt  .le. 0  start from scratch
c         mstrt  .ge. 1  start from previous solution
c     (for unsteday flows, mstrt=2 starts from steady-state solution)
c     msave  : save option for future restart
c         msave  > 0 frequency of updating restart file
c         msave  < 0 save restart only at end of run
c     mplot  > 0 save plot file
c     mcplpr > 0 produce line-plot of convergence history
c     mgtype : used for selecting v or w cycles
c         mgtype = 0  donot use multigrid
c                  1  use v-cycle for multigrid
c                  2  use w-cycle for multigrid
c     icau
c     i2dfl
c     ipatchg
c     iturb
c     ieng
c     ijet
c     iforce
c
      igrdfmt   = grdfmt
      mstrt     = fstrt 
      msave     = fsave 
      mavg      = favg 
      mplot     = fplot 
      mcplpr    = fcplpr
      mgtype    = fmgtyp
      icau      = cau
      i2dfl     = fi2dfl
      ipatchg   = fpatchg
      iturb     = fiturb
      iforce    = 1
c
cprec
      iprec     = fiprec

c     if iprecg is negative do not precondition in precong
      iprecg    = iprec
      iprec     = abs(iprec)
c
c---- initialize convergence output parameters
c
      nres      = 0
      ncyc      = 0
      totw      = 0.
c
c---- initialize global quantities
c
      scal      = 1.
      if (sref.le.0.0) sref = 1.0
c
cyho    add yaw and roll (beta and phi) capability

c       Two angles, either alpha and beta or alpha-total
c       and phi (roll), are sufficient to specify an attitude.
c       Specifying all three angles have a redundency.
c       Ordinarily one should enter alpha total and phi (roll)
c       or alpha and beta (yaw), but not all three.

c       'yaw' will be treated as beta, not alpha-y (or alpha-z here)
c       that is, beta is the rotation angle around vertical axis
c       (lift axis) of the wind axes system.
c       positive beta means turning nose to the left of pilot
c       that is, positive beta will result in positive w0 in a
c       typical plot3d grid sign convention (positive x-axis
c       points aft of the body from the front). (see Etkin, p11)

c       'roll' is the rotation angle around x axis of the body axes
c       system,
c       phi (roll) is positive anti-clock-wise measured from the
c       maneuvering plane to the body positive y axis, from
c       front looking aft.  When x positive axis is pointing aft of the
c       body, phi (roll) can be considered as an angle measured from
c       the positive y axis to the manevering plane with right handed
c       sign convention.
c
c       alpha, beta  and alpha-total, phi have following relationship;
c
c       alpha = atan (v0/u0) = atan (tan(alphat)*cos(phi))
c       beta  = asin (w0/V ) = asin (sin(alphat)*sin(phi))
c       where V = sqrt (u0**2 + v0**2 + w0**2)
c


      alpha     = al/rad
      yaw       = yaw/rad
      roll      = roll/rad
      ca        = cos(alpha)
      sa        = sin(alpha)
cyho    ca and sa will be redefined later
      vt        = 1.
      if (cflf.lt.0.) vt = 0.
      cflf      = abs(cflf)
      if (cflc.le.0.) cflc = cflf 
      cfl       = cflf
      gm1       = gamma - 1.
      dgm1      = 1./gm1
      gogm1     = gamma/gm1                                             
      rho0      = 1.
      p0        = 1.
      c0        = sqrt(gamma*p0/rho0)
      ei0       = p0/(gm1*rho0)
      rmq2      = rm*rm/2.                                              
      ptot0     = p0*(1.+gm1*rmq2)**gogm1                               
      rey       = reyl*scal*1.e+06
c
c **********  Jan, 1998   *******
c   put the complete form for pitch,roll and yaw angles
c
cyho
c     u0        = rm*c0*ca
c     v0        = rm*c0*sa
c     w0        = 0.
c
      if(alpha.eq.0.0.and.yaw.eq.0.0) then
        phi = roll
        alphas = alpha
        betas  = yaw
      else
        if(roll.ne.0.0) then
c
          phi    = roll
          alphat = alpha
          if(yaw.ne.0.0) then
c       convert this mixed up situation to alpha total & roll
            tanalpha =  sqrt(tan(alpha)**2+tan(yaw)**2
     .                  +(tan(alpha)*tan(yaw))**2 )
            alphat   = atan(tanalpha)
            phi = ((1.0-sign(1.,tan(alpha)))*0.5)*180./rad+
     .          sign(1.,tan(alpha))*asin(tan(yaw)*
     .          sqrt(1.+tan(alpha)**2)/tanalpha)
            phi      = phi+roll
          end if
          alphas  = atan(tan(alphat)*cos(phi))
          betas   = asin(sin(phi)*sin(alphat))
        else
c       specified roll is zero
          alphas        = alpha
          betas         = yaw
            tanalpha =  sqrt(tan(alpha)**2+tan(yaw)**2
     .                  +(tan(alpha)*tan(yaw))**2 )
            alphat   = atan(tanalpha)
            phi = ((1.0-sign(1.,tan(alpha)))*0.5)*180./rad+
     .          sign(1.,tan(alpha))*asin(tan(yaw)*
     .          sqrt(1.+tan(alpha)**2)/tanalpha)
        end if
      end if
c
c       (see NASA SP-3070)
c       direction cosines for unit velocity vector
c
c       li      = cos(alphat)
c               = cos(betas)*cos(alphas)
c       vertical
c       lj      = sin(alphat)*cos(phi)
c         or    = cos(alphat)*tan(alphas)
c         or    = cos(betas)*sin(alphas)
c       lateral
c       lk      = sin(alphat)*sin(phi)
c         or    = sin(betas)
c
        u0      = rm*c0*cos(betas)*cos(alphas)
        v0      = rm*c0*cos(betas)*sin(alphas)
        w0      = rm*c0*sin(betas)

c  redefine ca and sa for lift and drag calculation later

        ca      = cos(alphas)
        sa      = sin(alphas)


      h0        = gamma*ei0 + .5*(u0*u0 + v0*v0 + w0*w0)
c
cprec
c
c     scale ubmin, a preconditioning parameter by free-stream velocity
c
      ubmin2    = (ubmin*rm)**2
      ubcon2    = (ubcon*rm)**2
      ammax     = amachg
c
c     following moved here from data statement
c
      anutinf   = 1.341946
      zkinf     = 9.e-09
      ominf     = 1.e-06
c
c---- conditions for jet-nozzle (upstream)
c
c******************************************************
c
c     using the input quantities available at this point regarding
c     the fmg levels and the number of grids to be used in multigrid,
c     and the type of mg cycle, generate the information about
c     the mg-legs, beginning, ending, and solution grid-levels for
c     restriction and prolongation parts of the multigrid cycle
c
c
      call setseq ( mstrt,mcyc,mgrlev,ngrid,nlev,mgtype,mgrleg,
     .              nrleg,isoln,icoars,ifine,ibegr )
c 
c######################################################################
c master process
c     read and check input/output files and
c     send grid sizes to nodes
c######################################################################
c
      if (myrank.eq.master) then
c          
c---- grid file
c
      len = istlen (grdfil,80)
      if (igrdfmt.eq.1) then
        write (iwrit,'(/," formatted grid file   : ",a)')
     .                   grdfil(1:len)
      else
        write (iwrit,'(/," unformatted grid file : ",a)')
     .                   grdfil(1:len)
      endif
c
c---- map file
c
      len = istlen (mapfil,80)
      open (imapfil,file=mapfil(1:len),form='formatted')
      write (iwrit,'(" topological map file  : ",a)')
     .               mapfil(1:len)
      rewind imapfil
      read (imapfil,'(1x)',end=6002)
      read (imapfil,*,end=6002) ibloc
      if (ibloc.ne.nbloc) then
         write(iwrit,'(2x,"no. of blocks on map file ",i3,
     .     "do not match the no. of blocks on grid file",i3)')
     .     ibloc,nbloc
         call ERREXIT (nodes)
      endif
      go to 6003
 6002 write (iwrit,'("map file not found or empty")')
      call ERREXIT (nodes)
 6003 rewind imapfil
c
c---- restart file
c
      if (mstrt.ge.1) then
        len = istlen (rstfil,80)
        open (irstfil,file=rstfil(1:len),form='unformatted') 
        write (iwrit,'(" solution restart file : ",a)')
     .                 rstfil(1:len)
        rewind irstfil
        read (irstfil,end=6004) nres,jturb,j2dfl,nnorder,totptim
#ifdef BUILD_MPI
          call MPI_Bcast (nnorder,1,MPI_INTEGER,master,mycomm,ierr)
#endif
        if(mstrt.eq.2) totptim = 0.
c
        if (nres .gt. mres) then
          write (iwrit,'(" the no. of residuals for this run = ",
     .                   i2)' ) nres
          write (iwrit,'(" exceed the max residual dimension = ",
     .                i2)' ) mres
          call ERREXIT (nodes)
        endif
c
        if( nnorder .lt. ntorder. and. mstrt.lt.2 ) then
          write (iwrit,*)
          write (iwrit,'(2x,"order of time derivative specified (",
     .                   i2," )")' ) ntorder
          write (iwrit,'(2x,"is greater than the restart file (",
     .                i2," ) can handle.")' ) nnorder
c         call ERREXIT (nodes)
        endif
c
c
cvndifturb     if (jturb.ne.iturb) then
        if (jturb.ne.iturb.and.jturb.lt.iturb) then
          write (iwrit,'(" iturb from restart    = ",i2)') jturb
          write (iwrit,'(" not equal to this run = ",i2)') iturb
c       broadcast jturb to all processors
c
#ifdef BUILD_MPI
          call MPI_Bcast (jturb,1,MPI_INTEGER,master,mycomm,ierr)
#endif
c         call ERREXIT (nodes)
        endif
c
c
        read (irstfil,end=6004) (fpval(n,1),n=1,8)
        read (irstfil,end=6004) (fpval(n,2),n=1,5)
        read (irstfil,end=6004)
     .    (wcount(ires),rtrmsg(ires),hrmsg(ires),
     .    clg(ires),cmzg(ires),cdtg(ires),cdvg(ires),nsupg(ires),
     .    ires=1,nres  )
c
        totw = wcount(nres)
c
        read (irstfil,end=6004) nnbloc
        if (nnbloc .ne. nbloc) then
          write(iwrit,'(2x,"no. of blocks",i3," on restart file",
     .       " do not match no. of blocks",i3," on grid file")')
     .         nnbloc,nbloc
          call ERREXIT (nodes)
        endif
c
c       read in the dimensions from restart file and
c       check to see that dimensions on restart file for all
c       blocks are compatible with the dimensions on grid file
c
        read (irstfil,end=6004)
     .    (idum(ibloc),jdum(ibloc),kdum(ibloc),ibloc=1,nbloc)
c
c       do ibloc = 1,nbloc
c         if (idum(ibloc).ne.im(1,ibloc)+2.or.
c    .        jdum(ibloc).ne.jm(1,ibloc)+2.or.
c    .        kdum(ibloc).ne.km(1,ibloc)+2) then
c           write (iwrit,'(2x,"the dimensions",3i5,1x,"on restart",
c    .       1x,"file are incompatible with grid file dimensions",
c    .       1x,3i5,1x,"for the block = ",i3)') idum(ibloc),
c    .       jdum(ibloc),kdum(ibloc),im(1,ibloc)+2,
c    .       jm(1,ibloc)+2,km(1,ibloc)+2,ibloc
c           call ERREXIT (nodes)
c         endif
c       enddo
        go to 6005
 6004   write (iwrit,'("restart file not found or empty")')
        call ERREXIT (nodes)
c
c       can't do FMG with a restart
c
 6005   if (mstrt.gt.1 .and. nlev.gt.1) then
          write (iwrit,'(" Can Not Restart and use Full Multigrid")')
          write (iwrit,'(" Program Terminated")')
          call ERREXIT (nodes)
        endif
c
      endif
c
c---- plot output file (for RK error estimates)
c
cvn   if (mplot.gt.0) then
        len = istlen (pltfil,80)
        open (ipltfil,file=pltfil(1:len),form='formatted')
        write (iwrit,'(" plot output file      : ",a)')
     .                 pltfil(1:len)
        rewind ipltfil
cvn   endif
c
c---- save file
c
      if (msave.ne.0) then
        len = istlen (savfil,80)
        open (isavfil,file=savfil(1:len),form='unformatted')
        write (iwrit,'(" solution save file    : ",a)')
     .                 savfil(1:len)
        rewind isavfil
      endif
c
c     read and store the pointer arrays for time-dependent output data
c
      if(ntorder.ge.1) then
        rewind ipntfil
        do igrp=1,ngroup
          read(ipntfil,*) iblkgrp(igrp),
     .         ibeggrp(igrp),iendgrp(igrp),iskpgrp(igrp),
     .         jbeggrp(igrp),jendgrp(igrp),jskpgrp(igrp),
     .         kbeggrp(igrp),kendgrp(igrp),kskpgrp(igrp)
        enddo
c
c       move the pointers from the external designation to the
c       internal designation by adding 1 to their values
c
      write(iwrit,'("special output for time accurate calculations at:",
     .              /," ibloc  ibeg  iend iskip  jbeg  jend  jskip  ",
     .                        "kbeg  kend kskip")')
        do igrp=1,ngroup
          ibeggrp(igrp) = ibeggrp(igrp) + imn - 1
          iendgrp(igrp) = iendgrp(igrp) + imn - 1
          jbeggrp(igrp) = jbeggrp(igrp) + jmn - 1
          jendgrp(igrp) = jendgrp(igrp) + jmn - 1
          if(i2dfl.eq.1) then
            kbeggrp(igrp) = kmn + 1
            kendgrp(igrp) = kmn + 1
          else
            kbeggrp(igrp) = kbeggrp(igrp) + kmn - 1
            kendgrp(igrp) = kendgrp(igrp) + kmn - 1
          end if
          write(iwrit,'(10i6)') iblkgrp(igrp),
     .            ibeggrp(igrp),iendgrp(igrp),iskpgrp(igrp),
     .            jbeggrp(igrp),jendgrp(igrp),jskpgrp(igrp),
     .            kbeggrp(igrp),kendgrp(igrp),kskpgrp(igrp)
        enddo
      endif
c
c
c---- patched grid files
c
      if (ipatchg.gt.0) then
        len = istlen (pchfil1,80)
        open (ipfil1,file=pchfil1(1:len),form='formatted')
        write (iwrit,'(" grid patching file 1  : ",a)')
     .                 pchfil1(1:len)
        rewind ipfil1
        read (ipfil1,*,end=6006) ngridp
        if (ngridp.gt.mgrlev .and. ngridp.lt.ngrid) then
          write (iwrit,'(2x,"grid level  conflicts dimensions ")')
          write (iwrit,'(2x,"ngridp ngrid  mgrlev ")')
          write (iwrit,'(2x,3i7)') ngridp,ngrid,mgrlev
          call ERREXIT (nodes)
        endif
c
        len = istlen (pchfil2,80)
        open (ipfil2,file=pchfil2(1:len),form='formatted')
        write (iwrit,'(" grid patching file 2  : ",a)')
     .                 pchfil2(1:len)
        rewind ipfil2
        read (ipfil2,*,end=6006) ngridpp
        if (ngridpp.gt.mgrlev .and. ngridpp.lt.ngrid) then
          write (iwrit,'(2x,"grid level  conflicts dimensions ")')
          write (iwrit,'(2x,"ngridpp ngrid  mgrlev ")')
          write (iwrit,'(2x,3i7)') ngridpp,ngrid,mgrlev
          call ERREXIT (nodes)
        endif

        go to 6007
 6006   write (iwrit,'("patch file not found or empty")')
        call ERREXIT (nodes)
 6007   rewind ipfil1
        rewind ipfil2
      endif
c
c---- write out restart data
c
      if (mstrt.eq.0) then
        write (iwrit,'(/," Start solution from scratch")')
      else
        write (iwrit,'(/," restart solution:",/,
     .   "         rm          reyl       alpha       gamma",
     .    "       tref         scal       sref")')
        write (iwrit,'(1x,7f12.5)') fpval(2,1),fpval(6,1),
     .    fpval(3,1)*rad,fpval(1,1),fpval(8,1),fpval(1,2),
     .    fpval(5,2)
      endif
c
      write (iwrit,'(//)')
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c---- write grid block assignments
c
      write (iwrit,'(" grid block assignments",/,
     .               "  node  blocks    points")')
      do inode=1,nnodes
        nblk = 0
        npts = 0
        do ibloc=1,nbloc
          if (nodes(ibloc).eq.nodes(ndlist+inode)) then
            nblk = nblk +1
            npts = npts +im(1,ibloc)*jm(1,ibloc)*km(1,ibloc)
          endif
        enddo
        write (iwrit,'(i5,i7,i12)') inode,nblk,npts
      enddo
#endif
c
c---- read imap array
c
      do ibloc=1,nbloc
c
c     adjust the dimensions of the blocks to account for
c     starting point being other than 1 on each face
c
        im(1,ibloc) = im(1,ibloc) + imn - 1
        jm(1,ibloc) = jm(1,ibloc) + jmn - 1
        km(1,ibloc) = km(1,ibloc) + kmn - 1
c
c       need these on master for start and savsol2
c
        imp2(1,ibloc) = im(1,ibloc) +2
        jmp2(1,ibloc) = jm(1,ibloc) +2
        kmp2(1,ibloc) = km(1,ibloc) +2
      enddo
c
      igrid = 1
      call imapin (imn,jmn,kmn,im,jm,km,
     .             mgrlev,ngrid,igrid,mbloc,nbloc,imap,msegt,
     .             msegn,nseg,mx1dwk,const1,const2,const3,const4,
     .             const5,wk2d(1,1),wk2d(1,2),wk2d(1,3),
     .             wk2d(1,4),wk2d(1,5),wk2d(1,6),ivisb,
     .             iturbb,itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,nodes)
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c---- send block-to-process mapping and grid sizes to nodes
c     where dimensions are checked. Wait for OK response.
c     If a message with the tag TAG_ERROR is received, then
c     the dimensions were exceeded on the node process, and
c     we exit.
c
      intval(1) = nbloc
      intval(2) = nnodes
      icnt      = 2
      do ibloc=1,nbloc
        intval(icnt+1) = nodes(ibloc)
        intval(icnt+2) = im(1,ibloc)
        intval(icnt+3) = jm(1,ibloc)
        intval(icnt+4) = km(1,ibloc)
        intval(icnt+5) = nseg(ibloc)
        intval(icnt+6) = iturbb(ibloc)
        intval(icnt+7) = ivisb(1,ibloc)
        intval(icnt+8) = ivisb(2,ibloc)
        intval(icnt+9) = ivisb(3,ibloc)
        icnt           = icnt + 9
      enddo
      do inode=1,nnodes
        icnt         = icnt + 1
        intval(icnt) = nodes(ndlist+inode)
      enddo
c
      do inode=1,nnodes
#  ifdef BUILD_PVM
        call PVMFpsend (nodes(ndlist+inode),TAG_DIMS,
     .                  intval,icnt,INTEGER4,ierr)
        call PVMFprecv (nodes(ndlist+inode),-1,
     .                  errmsg,128,BYTE1,
     .                  itid,itag,ilen,ierr)
#  else
      if (nodeid(inode)-1.ne.master) then
        call MPI_Send (intval,icnt,MPI_INTEGER,
cjm     .                 nodes(ndlist+inode),TAG_DIMS,
     .                 nodes(ndlist+inode)-1,TAG_DIMS,
     .                 mycomm,ierr)
c
c       send the array pointers for time-dependent output to nodes
c
        if (ntorder.ge.1) then
          call MPI_Send (iblkgrp,ngroup,MPI_INTEGER,
     .    nodes(ndlist+inode)-1,TAG_GRP,mycomm,ierr)

          call MPI_Send (ibeggrp,ngroup,MPI_INTEGER,
     .    nodes(ndlist+inode)-1,TAG_GRP+1,mycomm,ierr)
          call MPI_Send (iendgrp,ngroup,MPI_INTEGER,
     .    nodes(ndlist+inode)-1,TAG_GRP+2,mycomm,ierr)
          call MPI_Send (iskpgrp,ngroup,MPI_INTEGER,
     .    nodes(ndlist+inode)-1,TAG_GRP+3,mycomm,ierr)

          call MPI_Send (jbeggrp,ngroup,MPI_INTEGER,
     .    nodes(ndlist+inode)-1,TAG_GRP+4,mycomm,ierr)
          call MPI_Send (jendgrp,ngroup,MPI_INTEGER,
     .    nodes(ndlist+inode)-1,TAG_GRP+5,mycomm,ierr)
          call MPI_Send (jskpgrp,ngroup,MPI_INTEGER,
     .    nodes(ndlist+inode)-1,TAG_GRP+6,mycomm,ierr)

          call MPI_Send (kbeggrp,ngroup,MPI_INTEGER,
     .    nodes(ndlist+inode)-1,TAG_GRP+7,mycomm,ierr)
          call MPI_Send (kendgrp,ngroup,MPI_INTEGER,
     .    nodes(ndlist+inode)-1,TAG_GRP+8,mycomm,ierr)
          call MPI_Send (kskpgrp,ngroup,MPI_INTEGER,
     .    nodes(ndlist+inode)-1,TAG_GRP+9,mycomm,ierr)
        endif

        call MPI_Recv (errmsg,128,MPI_CHARACTER,
cjm     .                 nodes(ndlist+inode),MPI_ANY_TAG,
     .                 nodes(ndlist+inode)-1,MPI_ANY_TAG,
     .                 mycomm,istat,ierr)
        itag = istat(MPI_TAG)
cjm
      end if
#  endif
cjm
      if (inode.ne.1) then
        if (itag.eq.TAG_ERROR) then
          len = istlen (errmsg,128)
          write (iwrit,'(/," error from node",i4,": ",a)')
     .      inode,errmsg(1:len)
          call ERREXIT (nodes)
        endif
cjm
      end if
      enddo
c
c######################################################################
c node processes
c     receive block sizes from root, set up pointers
c     and check dimensions
c######################################################################
c
#  ifdef BUILD_PVM
      else
#  else
      end if
#  endif
c
#  ifdef BUILD_PVM
      call PVMFprecv (master,-1,intval,mxintln,INTEGER4,
     .                itid,itag,ilen,ierr)
#  else
      if (myrank.ne.master) then
      call MPI_Recv (intval,mxintln,MPI_INTEGER,master,
     .               MPI_ANY_TAG,mycomm,istat,ierr)
      itag = istat(MPI_TAG)
      end if
#  endif
cjm
      if (myrank.ne.master) then
      if (itag.ne.TAG_DIMS)
     .  call ERREXIT (nodes)
cjm
      end if
c
c       receive the array pointers for time-dependent output from master
      if (myrank.ne.master) then
        if (ntorder.ge.1) then
          call MPI_Recv (iblkgrp,ngroup,MPI_INTEGER,
     .    master,TAG_GRP,mycomm,istat,ierr)

          call MPI_Recv (ibeggrp,ngroup,MPI_INTEGER,
     .    master,TAG_GRP+1,mycomm,istat,ierr)
          call MPI_Recv (iendgrp,ngroup,MPI_INTEGER,
     .    master,TAG_GRP+2,mycomm,istat,ierr)
          call MPI_Recv (iskpgrp,ngroup,MPI_INTEGER,
     .    master,TAG_GRP+3,mycomm,istat,ierr)

          call MPI_Recv (jbeggrp,ngroup,MPI_INTEGER,
     .    master,TAG_GRP+4,mycomm,istat,ierr)
          call MPI_Recv (jendgrp,ngroup,MPI_INTEGER,
     .    master,TAG_GRP+5,mycomm,istat,ierr)
          call MPI_Recv (jskpgrp,ngroup,MPI_INTEGER,
     .    master,TAG_GRP+6,mycomm,istat,ierr)

          call MPI_Recv (kbeggrp,ngroup,MPI_INTEGER,
     .    master,TAG_GRP+7,mycomm,istat,ierr)
          call MPI_Recv (kendgrp,ngroup,MPI_INTEGER,
     .    master,TAG_GRP+8,mycomm,istat,ierr)
          call MPI_Recv (kskpgrp,ngroup,MPI_INTEGER,
     .    master,TAG_GRP+9,mycomm,istat,ierr)

        endif
      endif

      nbloc  = intval(1)
      nnodes = intval(2)
      icnt   = 2
      do ibloc=1,nbloc
        nodes(ibloc)   = intval(icnt+1)
        im(1,ibloc)    = intval(icnt+2)
        jm(1,ibloc)    = intval(icnt+3)
        km(1,ibloc)    = intval(icnt+4)
        nseg(ibloc)    = intval(icnt+5)
        iturbb(ibloc)  = intval(icnt+6)
        ivisb(1,ibloc) = intval(icnt+7)
        ivisb(2,ibloc) = intval(icnt+8)
        ivisb(3,ibloc) = intval(icnt+9)
        icnt           = icnt + 9
      enddo
      do inode=1,nnodes
        icnt                = icnt + 1
        nodes(ndlist+inode) = intval(icnt)
      enddo
#endif
c
c     In the next do loop generate the pointers for any blocks
c     which are on this node. The pointers for blocks on other
c     nodes are simply set to 1, since they are not used by
c     this process. The variable jbloc is used as the block
c     counter on a given node.
c
      jbloc = 0
      do 100 ibloc = 1,nbloc
c
c     initialize the pointers for the finest grid 
c
# ifdef BUILD_MPI
      if (nodes(ibloc)-1.ne.myrank .or. jbloc.eq.0) then
# else
      if (nodes(ibloc).ne.myrank .or. jbloc.eq.0) then
# endif
        m1np(1,ibloc) = 1
        m2np(1,ibloc) = 1
        m3np(1,ibloc) = 1
        m4np(1,ibloc) = 1
        m5np(1,ibloc) = 1
c
        m1cc(1,ibloc) = 1
        m2cc(1,ibloc) = 1
        m3cc(1,ibloc) = 1
        m4cc(1,ibloc) = 1
        m5cc(1,ibloc) = 1
        m5cg(1,ibloc) = 1
        m1fg(  ibloc) = 1
        m4fg(  ibloc) = 1
        m5fg(  ibloc) = 1
        m1fgt(1,ibloc)= 1
        m5fgt(1,ibloc)= 1
c#ifdef  RKPHY
        m1fgtrk(1,ibloc)= 1
        m5fgtrk(1,ibloc)= 1
c#endif
        do iord=2,mtorder
         m1fgt(iord,ibloc) = m1fgt(iord-1,ibloc) +    (im(1,ibloc)+2)
     .                            * (jm(1,ibloc)+2) * (km(1,ibloc)+2)
         m5fgt(iord,ibloc) = m5fgt(iord-1,ibloc) +    (im(1,ibloc)+2)
     .                            * (jm(1,ibloc)+2) * (km(1,ibloc)+2) *5
        enddo
c
c#ifdef  RKPHY
      if(nrkstg.gt.1) then
        do istg=2,nrkstg
         m1fgtrk(istg,ibloc) = m1fgtrk(istg-1,ibloc) +   (im(1,ibloc)+2)
     .                            * (jm(1,ibloc)+2) * (km(1,ibloc)+2)
         m5fgtrk(istg,ibloc) = m5fgtrk(istg-1,ibloc) +   (im(1,ibloc)+2)
     .                            * (jm(1,ibloc)+2) * (km(1,ibloc)+2) *5
        enddo
       endif
c#endif
c
        m1is(1,ibloc) = 1
        m1js(1,ibloc) = 1
        m1ks(1,ibloc) = 1
        m3is(1,ibloc) = 1
        m3js(1,ibloc) = 1
        m3ks(1,ibloc) = 1
c
      else
        ltempnp = imp1(ngrid,jbloc)*jmp1(ngrid,jbloc)*
     .            kmp1(ngrid,jbloc)
        ltempcc = imp2(ngrid,jbloc)*jmp2(ngrid,jbloc)*
     .            kmp2(ngrid,jbloc)
c
        m1np(1,ibloc) = m1np(ngrid,jbloc) + ltempnp
        m2np(1,ibloc) = m2np(ngrid,jbloc) + ltempnp*2
        m3np(1,ibloc) = m3np(ngrid,jbloc) + ltempnp*3
        m4np(1,ibloc) = m4np(ngrid,jbloc) + ltempnp*4
        m5np(1,ibloc) = m5np(ngrid,jbloc) + ltempnp*5
c
        m1cc(1,ibloc) = m1cc(ngrid,jbloc) + ltempcc
        m2cc(1,ibloc) = m2cc(ngrid,jbloc) + ltempcc*2
        m3cc(1,ibloc) = m3cc(ngrid,jbloc) + ltempcc*3
        m4cc(1,ibloc) = m4cc(ngrid,jbloc) + ltempcc*4
        m5cc(1,ibloc) = m5cc(ngrid,jbloc) + ltempcc*5
        m5cg(1,ibloc) = m5cg(ngrid,jbloc) + ltempcc*5
        m1fg(  ibloc) = m1fg(      jbloc) + imp2(1,jbloc)
     .                    * jmp2(1,jbloc) * kmp2(1,jbloc)
        m4fg(  ibloc) = m4fg(      jbloc) + imp2(1,jbloc)
     .                    * jmp2(1,jbloc) * kmp2(1,jbloc) * 4
        m5fg(  ibloc) = m5fg(      jbloc) + imp2(1,jbloc)
     .                    * jmp2(1,jbloc) * kmp2(1,jbloc) * 5
        m1fgt(1   ,ibloc) = m1fgt(mtorder,jbloc) + imp2(1,jbloc)
     .                    * jmp2(1,jbloc) * kmp2(1,jbloc)
        m5fgt(1   ,ibloc) = m5fgt(mtorder,jbloc) + imp2(1,jbloc)
     .                    * jmp2(1,jbloc) * kmp2(1,jbloc)*5
c
c#ifdef  RKPHY
        m1fgtrk(1   ,ibloc) = m1fgtrk(nrkstg,jbloc) + imp2(1,jbloc)
     .                    * jmp2(1,jbloc) * kmp2(1,jbloc)
        m5fgtrk(1   ,ibloc) = m5fgtrk(nrkstg,jbloc) + imp2(1,jbloc)
     .                    * jmp2(1,jbloc) * kmp2(1,jbloc)*5
c#endif
c
        do iord=2,mtorder
         m1fgt(iord,ibloc) = m1fgt(iord-1,ibloc) +    (im(1,ibloc)+2)
     .                            * (jm(1,ibloc)+2) * (km(1,ibloc)+2)
         m5fgt(iord,ibloc) = m5fgt(iord-1,ibloc) +    (im(1,ibloc)+2)
     .                            * (jm(1,ibloc)+2) * (km(1,ibloc)+2) *5
        enddo
c
c#ifdef  RKPHY
      if(nrkstg.gt.1) then
        do istg=2,nrkstg
         m1fgtrk(istg,ibloc) = m1fgtrk(istg-1,ibloc) +   (im(1,ibloc)+2)
     .                            * (jm(1,ibloc)+2) * (km(1,ibloc)+2)
         m5fgtrk(istg,ibloc) = m5fgtrk(istg-1,ibloc) +   (im(1,ibloc)+2)
     .                            * (jm(1,ibloc)+2) * (km(1,ibloc)+2) *5
        enddo
      endif
c#endif
c
c       pointers for 2-d (surface) arrays
c       for i,j and k planes
c
        m1is(1,ibloc) = m1is(ngrid,jbloc) +
     .                  jmp1(ngrid,jbloc) * kmp1(ngrid,jbloc) 
        m1js(1,ibloc) = m1js(ngrid,jbloc) +
     .                  imp1(ngrid,jbloc) * kmp1(ngrid,jbloc) 
        m1ks(1,ibloc) = m1ks(ngrid,jbloc) +
     .                  imp1(ngrid,jbloc) * jmp1(ngrid,jbloc) 
        m3is(1,ibloc) = m3is(ngrid,jbloc) +
     .                  jmp1(ngrid,jbloc) * kmp1(ngrid,jbloc) *3
        m3js(1,ibloc) = m3js(ngrid,jbloc) +
     .                  imp1(ngrid,jbloc) * kmp1(ngrid,jbloc) *3
        m3ks(1,ibloc) = m3ks(ngrid,jbloc) +
     .                  imp1(ngrid,jbloc) * jmp1(ngrid,jbloc) *3
c
      endif
c
c     set the pointers for remaining grids for ibloc
c     through the subroutine "pointer"
c
      call pointer (ngrid,
     .     m1np(1,ibloc), m2np(1,ibloc), m3np(1,ibloc),
     .     m4np(1,ibloc), m5np(1,ibloc), 
     .     m1cc(1,ibloc), m2cc(1,ibloc), m3cc(1,ibloc),
     .     m4cc(1,ibloc), m5cc(1,ibloc), m5cg(1,ibloc),
     .     m1is(1,ibloc), m1js(1,ibloc), m1ks(1,ibloc),
     .     m3is(1,ibloc), m3js(1,ibloc), m3ks(1,ibloc),
     .     im  (1,ibloc), jm  (1,ibloc), km  (1,ibloc),
     .     imp1(1,ibloc), jmp1(1,ibloc), kmp1(1,ibloc),
     .     imp2(1,ibloc), jmp2(1,ibloc), kmp2(1,ibloc) )
c
#ifdef BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) jbloc = ibloc
#else
      if (nodes(ibloc).eq.myrank) jbloc = ibloc
#endif
  100 continue
c
c     check to see if size of the node centered, cell centered,
c           coarse grid, and surface grid arrays are within
c           the allocated size
c
      m1npmx = m1np(ngrid,jbloc)   +imp1(ngrid,jbloc)*jmp1(ngrid,jbloc)
     .                                               *kmp1(ngrid,jbloc)
      m1ccmx = m1cc(ngrid,jbloc)   +imp2(ngrid,jbloc)*jmp2(ngrid,jbloc)
     .                                               *kmp2(ngrid,jbloc)
      m5cgmx = m5cg(ngrid,jbloc) +5*imp2(ngrid,jbloc)*jmp2(ngrid,jbloc)
     .                                               *kmp2(ngrid,jbloc)
      m3ismx = m3is(ngrid,jbloc) +3*jmp1(ngrid,jbloc)*kmp1(ngrid,jbloc)
      m3jsmx = m3js(ngrid,jbloc) +3*kmp1(ngrid,jbloc)*imp1(ngrid,jbloc)
      m3ksmx = m3ks(ngrid,jbloc) +3*imp1(ngrid,jbloc)*jmp1(ngrid,jbloc)
c
      if (m1npmx .gt. mxsizn) then
        itag = TAG_ERROR
        write (errmsg,'("the storage required for the node centered "
     .    "variables is larger than the allocated size")')
      else if (m1ccmx .gt. mxsizc) then
        itag = TAG_ERROR
        write (errmsg,'("the storage required for the cell centered "
     .     "variables is larger than the allocated size")')
      else if (m5cgmx .gt. 5*mxsizc) then
        itag = TAG_ERROR
        write (errmsg,'("the storage required for the coarse grid "
     .    "variables is larger than the allocated size")')
      else if (m3ismx .gt. mxs2d3c) then
        itag = TAG_ERROR
        write (errmsg,'("the storage required for the i=imin face "
     .    "variables is larger than the allocated size")')
      else if (m3jsmx .gt. mxs2d3c) then
        itag = TAG_ERROR
        write (errmsg,'("the storage required for the j=jmin face "
     .    "variables is larger than the allocated size")')
      else if (m3ksmx .gt. mxs2d3c) then
        itag = TAG_ERROR
        write (errmsg,'("the storage required for the k=kmin face "
     .    "variables is larger than the allocated size")')
      else
        itag = TAG_OK
        write (errmsg,'("storage OK for node",i5)') myrank
      endif
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c     send error message or OK tag to master
c
#  ifdef BUILD_PVM
      call PVMFpsend (master,itag,errmsg,128,BYTE1,ierr)
#  else
cjm
      if (myrank.ne.master) then
      call MPI_Send (errmsg,128,MPI_CHARACTER,master,
     .               itag,mycomm,ierr)
cjm
      end if
#  endif
c
      if (iwrit.gt.0) then
        len = istlen (errmsg,128)
        write (iwrit,'(//," ",a)') errmsg(1:len)
      endif
#else
c
c     IN-CORE version - stop if error
c
      if (itag.eq.TAG_ERROR) then
        len = istlen (errmsg,128)
        write (iwrit,'(" ",a)') errmsg(1:len)
        stop
      endif
#endif
c
#  ifdef BUILD_MPI
#  else
      endif
#  endif
c
c######################################################################
c all processes
c     print grid block dimensions
c######################################################################
c
      if (iwrit.gt.0) then
c
      write (iwrit,'(/," grid block dimensions",/,
     .  " block  imax   jmax   kmax     total   node   nodeid")')
      do ibloc=1,nbloc
# ifdef BUILD_MPI
        if (myrank.eq.master .or. nodes(ibloc)-1.eq.myrank) then
# else
        if (myrank.eq.master .or. nodes(ibloc).eq.myrank) then
# endif
          do i=1,nnodes
            if (nodes(ibloc).eq.nodes(ndlist+i)) then
              inode = i
              go to 1001
            endif
          enddo
 1001     continue
          ipts = im(1,ibloc) - imn + 1
          jpts = jm(1,ibloc) - jmn + 1
          kpts = km(1,ibloc) - kmn + 1
          write (iwrit,'(1x,i4,3i7,i10,i6,i10)') ibloc,
     .      ipts,jpts,kpts,ipts*jpts*kpts,inode,nodes(ibloc)
        endif
      enddo
c
      endif
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c 
c######################################################################
c master process
c     dimensions are OK, so send imap data to all nodes
c######################################################################
c
      if (myrank.eq.master) then
c
c     just in case:
      if (mxsiz5c.lt.5*msegn*nbloc) then
        write(iwrit,'("STRANGE: w array size < 5*msegn*nbloc")')
        call ERREXIT (nodes)
      endif
c
      icnt = 0
      jcnt = 0
      do ibloc=1,nbloc
        do igrid=1,ngrid
          intval(icnt+1) = itrb1(igrid,ibloc)
          intval(icnt+2) = itrb2(igrid,ibloc)
          intval(icnt+3) = jtrb1(igrid,ibloc)
          intval(icnt+4) = jtrb2(igrid,ibloc)
          intval(icnt+5) = ktrb1(igrid,ibloc)
          intval(icnt+6) = ktrb2(igrid,ibloc)
          icnt           = icnt + 6
        enddo
        do iseg=1,msegn
          w(jcnt+1)      = const1(iseg,ibloc)
          w(jcnt+2)      = const2(iseg,ibloc)
          w(jcnt+3)      = const3(iseg,ibloc)
          w(jcnt+4)      = const4(iseg,ibloc)
          w(jcnt+5)      = const5(iseg,ibloc)
          jcnt           = jcnt + 5
        enddo
      enddo
      icnt = 6*ngrid*nbloc
      jcnt = 5*msegn*nbloc
      ncnt = msegt*msegn*mbloc*ngrid
c
#  ifdef BUILD_PVM
      call PVMFinitsend (PvmDataInPlace,ierr)
      call PVMFpack (INTEGER4,imap,ncnt,1,ierr)
      call PVMFmcast (nnodes,nodes(ndlist+1),TAG_IMAP,ierr)
c
      call PVMFinitsend (PvmDataInPlace,ierr)
      call PVMFpack (INTEGER4,intval,icnt,1,ierr)
      call PVMFpack (RTYPE,w,jcnt,1,ierr)
      call PVMFmcast (nnodes,nodes(ndlist+1),TAG_IMAP,ierr)
#  else
cjm      do inode=1,nnodes
      do inode=2,nnodes
        call MPI_Send (imap,ncnt,MPI_INTEGER,
cjm     .                 nodes(ndlist+inode),TAG_IMAP,
     .                 nodes(ndlist+inode)-1,TAG_IMAP,
     .                 mycomm,ierr)
        call MPI_Send (intval,icnt,MPI_INTEGER,
cjm     .                 nodes(ndlist+inode),TAG_IMAP,
     .                 nodes(ndlist+inode)-1,TAG_IMAP,
     .                 mycomm,ierr)
        call MPI_Send (w,jcnt,RTYPE,
cjm     .                 nodes(ndlist+inode),TAG_IMAP,
     .                 nodes(ndlist+inode)-1,TAG_IMAP,
     .                 mycomm,ierr)
      enddo
#  endif
c
c######################################################################
c node processes
c     receive the imap data from root
c######################################################################
c
#  ifdef BUILD_PVM
      else
#  else
      end if
#  endif
c
      icnt = msegt*msegn*mbloc*ngrid
#  ifdef BUILD_PVM
      call PVMFprecv (master,-1,imap,icnt,INTEGER4,
     .                itid,itag,ilen,ierr)
#  else
cjm
      if (myrank.ne.master) then
      call MPI_Recv (imap,icnt,MPI_INTEGER,master,
     .               MPI_ANY_TAG,mycomm,istat,ierr)
      itag = istat(MPI_TAG)
      end if
#  endif
cjm
      if (myrank.ne.master) then
      if (itag.ne.TAG_IMAP)
     .  call ERREXIT (nodes)
cjm
      end if
c
#  ifdef BUILD_PVM
      call PVMFrecv (master,TAG_IMAP,ierr)
      call PVMFunpack (INTEGER4,intval,6*ngrid*nbloc,1,ierr)
      call PVMFunpack (RTYPE,w,5*msegn*nbloc,1,ierr)
#  else
cjm
      if (myrank.ne.master) then
      call MPI_Recv (intval,6*ngrid*nbloc,MPI_INTEGER,master,
     .               TAG_IMAP,mycomm,istat,ierr)
      call MPI_Recv (w,5*msegn*nbloc,RTYPE,master,
     .               TAG_IMAP,mycomm,istat,ierr)
      end if
#  endif
c
      icnt = 0
      jcnt = 0
      do ibloc=1,nbloc
        do igrid=1,ngrid
          itrb1(igrid,ibloc)  = intval(icnt+1)
          itrb2(igrid,ibloc)  = intval(icnt+2)
          jtrb1(igrid,ibloc)  = intval(icnt+3)
          jtrb2(igrid,ibloc)  = intval(icnt+4)
          ktrb1(igrid,ibloc)  = intval(icnt+5)
          ktrb2(igrid,ibloc)  = intval(icnt+6)
          icnt                = icnt + 6
        enddo
        do iseg=1,msegn
          const1(iseg,ibloc)  = w(jcnt+1)
          const2(iseg,ibloc)  = w(jcnt+2)
          const3(iseg,ibloc)  = w(jcnt+3)
          const4(iseg,ibloc)  = w(jcnt+4)
          const5(iseg,ibloc)  = w(jcnt+5)
          jcnt                = jcnt + 5
        enddo
      enddo
c
#  ifdef BUILD_PVM
      endif
#  endif
#endif
c
c######################################################################
c all processes
c     print turbulent zones
c######################################################################
c
      if (iwrit.gt.0) then
c
        if (iturb.eq.1) then
          write(iwrit,'(/," using Baldwin-Lomax turbulence model")')
        else if (iturb.eq.2) then
          write(iwrit,'(/," using Spalart turbulence model")')
        else if (iturb.eq.3) then
          write(iwrit,'(/," using K-Omega turbulence model")')
        else
          write(iwrit,'(/," turbulence not computed")')
        endif
        write(iwrit,'(" turbulent zones:",/,
     .    "              Index Range          ",
     .    "            Dimensions ",/,
     .    "  blk   i1   i2   j1   j2   k1   k2",
     .    "        im   jm   km turb")')
c
        do ibloc=1,nbloc
# ifdef BUILD_MPI
          if (myrank.eq.master .or. nodes(ibloc)-1.eq.myrank)
# else
          if (myrank.eq.master .or. nodes(ibloc).eq.myrank)
# endif
     .      write(iwrit,'(1x,i4,6i5,5x,4i5)') ibloc,
     .        itrb1(1,ibloc),itrb2(1,ibloc),
     .        jtrb1(1,ibloc),jtrb2(1,ibloc),
     .        ktrb1(1,ibloc),ktrb2(1,ibloc),
     .        im(1,ibloc),jm(1,ibloc),km(1,ibloc),iturbb(ibloc)
        enddo
c
      endif
c
c
c----------------------------------------------------------------
c
c     initial setup for patched interfaces
c     read the patcher information from ipfil1 and ipfil2 and
c     rearrange it in a form suitable for the tlns3d flow code
c
      if (ipatchg.le.0) go to 106
c
c######################################################################
c master process
c     read and check dimensions for patched grid related items
c     send this information to nodes
c######################################################################
c
      if (myrank.eq.master) then
c
c       read the patch files to create the pointers
c
        call setpchi (imn,jmn,kmn,mgrlev,mbloc,ngrid,nbloc,
     .         m1pch1,m1pch2,lswpchb,lswpche,
     .         mxtpchc,ntpchcb,npchcbf,
     .         mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .         lspchb2,lspchf2,nodes )
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
#  ifdef BUILD_PVM
c
      write (iwrit,'("PVM not supported for patched grids")')
      call ERREXIT (nodes)
#  else
c
cjm      do inode=1,nnodes
      do inode=2,nnodes
          call MPI_Send (ntpchcb,mbloc*ngrid,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   nodes(ndlist+inode)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (npchcbf,6*mbloc*ngrid,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   nodes(ndlist+inode)-1,TAG_PATCH,
     .                   mycomm,ierr)
c
          call MPI_Send (lswpchb,mbloc*ngrid,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   nodes(ndlist+inode)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (lswpche,mbloc*ngrid,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   nodes(ndlist+inode)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (nswpchb,mbloc*ngrid,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   nodes(ndlist+inode)-1,TAG_PATCH,
     .                   mycomm,ierr)
c
          call MPI_Send (m1pch1,mbloc*ngrid,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   nodes(ndlist+inode)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (m1pch2,mbloc*ngrid,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   nodes(ndlist+inode)-1,TAG_PATCH,
     .                   mycomm,ierr)
c
          call MPI_Send (lspchb1,mxpchs*ngrid,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   nodes(ndlist+inode)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (lspchb2,mxpchs*ngrid,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   nodes(ndlist+inode)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (lspchf1,mxpchs*ngrid,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   nodes(ndlist+inode)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (lspchf2,mxpchs*ngrid,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   nodes(ndlist+inode)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (npchitm,mxpchs*ngrid,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   nodes(ndlist+inode)-1,TAG_PATCH,
     .                   mycomm,ierr)
c
      enddo
c
#  endif
c
#  ifdef BUILD_PVM
      else
#  else
      end if
#  endif
c######################################################################
c node processes
c     receive dimensions for patched grid related items
c     for the compute nodes
c######################################################################
c
c
#  ifdef BUILD_PVM
c
      write (iwrit,'("PVM not supported for patched grids")')
      call ERREXIT (nodes)
#  else
c
cjm
      if (myrank.ne.master) then
        call MPI_Recv (ntpchcb,mbloc*ngrid,MPI_INTEGER,master,
     .                 MPI_ANY_TAG,mycomm,istat,ierr)
        itag = istat(MPI_TAG)
        if (itag.ne.TAG_PATCH)
     .    call ERREXIT (nodes)
        call MPI_Recv (npchcbf,6*mbloc*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
c
        call MPI_Recv (lswpchb,mbloc*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
        call MPI_Recv (lswpche,mbloc*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
        call MPI_Recv (nswpchb,mbloc*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
c
        call MPI_Recv (m1pch1,mbloc*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
        call MPI_Recv (m1pch2,mbloc*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
c
        call MPI_Recv (lspchb1,mxpchs*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
        call MPI_Recv (lspchb2,mxpchs*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
        call MPI_Recv (lspchf1,mxpchs*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
        call MPI_Recv (lspchf2,mxpchs*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
        call MPI_Recv (npchitm,mxpchs*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
cjm
      end if
c
#  endif
c
#endif
c
c     All nodes execute the next do loop (do 103) to
c     generate pointers for patching related arrays
c     for any blocks that reside on this node. The pointers
c     for blocks on other nodes are simply set to 1 since they
c     are not used by this process. The variable jbloc is used
c     as the block counter on a given node
c
      jbloc = 0
      do 103 ibloc = 1,nbloc
c       
c       initialize the pointer for the finest grid
# ifdef BUILD_MPI
        if (nodes(ibloc)-1.ne.myrank .or. jbloc.eq.0) then
# else
        if (nodes(ibloc).ne.myrank .or. jbloc.eq.0) then
#endif
          m1pch1(ibloc,1) = 1
          m1pch2(ibloc,1) = 1
        else
          ncnt1      = 0
          ncnt2      = 0
          if(ntpchcb(jbloc,ngrid).gt.0) then
cjun97    if(ntpchcb(jbloc,ngrid).gt.0) then
c           ncnt1    = 0
            do iface = 1,6
              ncnt1  = ncnt1 + npchcbf(iface,jbloc,ngrid)
            enddo
c   
c           ncnt2    = 0
            do lpchs = lswpchb(jbloc,ngrid)+1,lswpche(jbloc,ngrid)
              ncnt2  = ncnt2 + npchitm(lpchs,ngrid)
            enddo
          endif
cjun97    endif
c   
            m1pch1(ibloc,1) = m1pch1(jbloc,ngrid) + ncnt1
            m1pch2(ibloc,1) = m1pch2(jbloc,ngrid) + ncnt2
c
        endif
c   
c       set the pointers for remaining grids in multigrid
c       for ibloc
c
        do igrid=2,ngrid
c
          ncnt1      = 0
          ncnt2      = 0
#if defined BUILD_MPI
          if(myrank.eq.nodes(ibloc)-1) then
#else
          if(myrank.eq.nodes(ibloc)) then
#endif
          if(ntpchcb(ibloc,igrid).gt.0) then
cjun97    if(ntpchcb(ibloc,igrid).gt.0) then
c           ncnt1    = 0
            do iface = 1,6
              ncnt1  = ncnt1 + npchcbf(iface,ibloc,igrid-1)
            enddo
c
c           ncnt2   = 0
            do lpchs=lswpchb(ibloc,igrid-1)+1,lswpche(ibloc,igrid-1)
              ncnt2 = ncnt2 + npchitm(lpchs,igrid-1)
            enddo
cvn       endif
cvn       endif
c   
            m1pch1(ibloc,igrid) = m1pch1(ibloc,igrid-1) + ncnt1
            m1pch2(ibloc,igrid) = m1pch2(ibloc,igrid-1) + ncnt2
cnov97    endif
          else
            m1pch1(ibloc,igrid) = m1pch1(ibloc,igrid-1)
            m1pch2(ibloc,igrid) = m1pch2(ibloc,igrid-1)
          endif
cjun97    endif
          endif
        enddo
c
#ifdef BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) jbloc = ibloc
#else
      if (nodes(ibloc).eq.myrank) jbloc = ibloc
#endif
c     if (nodes(ibloc).eq.myrank) then
c         jbloc  = ibloc
c         do igrid=1,ngrid
c         call MPI_Send (m1pch1(ibloc,igrid),1,MPI_INTEGER,
c    .                   master,TAG_PATCH,
c    .                   mycomm,ierr)
c         call MPI_Send (m1pch2(ibloc,igrid),1,MPI_INTEGER,
c    .                   master,TAG_PATCH,
c    .                   mycomm,ierr)
c         enddo
c     endif
  103 continue
c
c
#if defined BUILD_MPI
#else
      endif
#endif
c
c------------------------------------------------------
c     this completes the computation of pointers for patching
c     on all processors. Next we read detailed patching info.
c------------------------------------------------------
c
c
c     print *,'before mpi_barrier'
c     call MPI_Barrier(MPI_COMM_WORLD,ierr1)
c     print *,'after  mpi_barrier'
c
c
c######################################################################
c  start a new master process to
c     read and check dimensions for patched grid related items
c     send this information to nodes
c######################################################################
c
      if (myrank.eq.master) then
c
c      do ibloc=1,nbloc
c      do igrid=1,ngrid
c      call MPI_Recv (m1pch1(ibloc,igrid),1,MPI_INTEGER,nodes(ibloc),
c    .                 TAG_PATCH,mycomm,istat,ierr)
c      call MPI_Recv (m1pch2(ibloc,igrid),1,MPI_INTEGER,nodes(ibloc),
c    .                 TAG_PATCH,mycomm,istat,ierr)
c      enddo
c      enddo
c
       rewind ipfil1
       rewind ipfil2

c      read the information about which cells are patched
c      from the file ipfil1. In addition,
c
c      read the information from ipfil2 about source cells which
c      contribute to cells identified in setpch1 from the file ipfil1
c
       read (ipfil1,*) iidum
       read (ipfil2,*) iidum
         do igrid = 1,ngrid
           read (ipfil1,*) iidum
           read (ipfil2,*) iidum
           read (ipfil2,*) iidum
           do 104 ibloc=1,nbloc
c
cvn        if(ntpchcb(ibloc,igrid).le.0) go to 104
c
c          ntpchcb(ibloc,igrid) values of ipatchc, jpatchc and kpatchc
c          are read in the subroutine setpch1
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)

#  ifdef BUILD_PVM
             call setpch1 (imn,jmn,kmn,mgrlev,mbloc,ngrid,nbloc,
     .            ibloc,igrid,ntpchcb,npchcbf,
     .            ipatchc,jpatchc,kpatchc,nodes)
#  else
         if (nodes(ibloc)-1.ne.myrank) then
             call setpch1 (imn,jmn,kmn,mgrlev,mbloc,ngrid,nbloc,
     .            ibloc,igrid,ntpchcb,npchcbf,
     .            tbuf1,tbuf2,tbuf3,nodes)
         else
             call setpch1 (imn,jmn,kmn,mgrlev,mbloc,ngrid,nbloc,
     .            ibloc,igrid,ntpchcb,npchcbf,
     .            ipatchc(m1pch1(ibloc,igrid)),
     .            jpatchc(m1pch1(ibloc,igrid)),
     .            kpatchc(m1pch1(ibloc,igrid)),nodes )
         end if
#  endif
#else
             call setpch1 (imn,jmn,kmn,mgrlev,mbloc,ngrid,nbloc,
     .            ibloc,igrid,ntpchcb,npchcbf,
     .            ipatchc(m1pch1(ibloc,igrid)),
     .            jpatchc(m1pch1(ibloc,igrid)),
     .            kpatchc(m1pch1(ibloc,igrid)),nodes )

#endif
               ncnt1 = 0
               do iface = 1,6
               ncnt1 = ncnt1 + npchcbf(iface,ibloc,igrid)
               enddo
c
c           npchitm(lpchs,igrid) values of ipitmb1,ipitmb2,
c           jpitmb1,jpitmb2,kpitmb1,kpitmb1,kpitmb2,frc
c           are read in by setpch2 for the patched surfaces
c           on current block
c      note: lswpchb,lspche loop is now inside subroutine setpch2
c
c                 litmbeg = m1pch2(ibloc,igrid)
#if defined(BUILD_PVM) || defined(BUILD_MPI)
#  ifdef BUILD_PVM
                  litmbeg = 1
                  call setpch2 (imn,jmn,kmn,mgrlev,mbloc,ngrid,nbloc,
     .            ibloc,igrid,lswpchb,lswpche,
     .            mxpchs,mpchitm,nswpchb,npchitm,
     .            lspchb1,lspchf1,lspchb2,lspchf2,
     .            ipitmb1(litmbeg),ipitmb2(litmbeg),
     .            jpitmb1(litmbeg),jpitmb2(litmbeg),
     .            kpitmb1(litmbeg),kpitmb2(litmbeg),frc(litmbeg),nodes )
#  else
         if (nodes(ibloc)-1.ne.myrank) then
                  call setpch2 (imn,jmn,kmn,mgrlev,mbloc,ngrid,nbloc,
     .            ibloc,igrid,lswpchb,lswpche,
     .            mxpchs,mpchitm,nswpchb,npchitm,
     .            lspchb1,lspchf1,lspchb2,lspchf2,
     ,            tbuf4,tbuf5,tbuf6,tbuf7,tbuf8,tbuf9,tbuf10,nodes)
         else
                  litmbeg = m1pch2(ibloc,igrid)
                  call setpch2 (imn,jmn,kmn,mgrlev,mbloc,ngrid,nbloc,
     .            ibloc,igrid,lswpchb,lswpche,
     .            mxpchs,mpchitm,nswpchb,npchitm,
     .            lspchb1,lspchf1,lspchb2,lspchf2,
     .            ipitmb1(litmbeg),ipitmb2(litmbeg),
     .            jpitmb1(litmbeg),jpitmb2(litmbeg),
     .            kpitmb1(litmbeg),kpitmb2(litmbeg),frc(litmbeg),nodes )
         end if
#  endif
#else
                  litmbeg = m1pch2(ibloc,igrid)
                  call setpch2 (imn,jmn,kmn,mgrlev,mbloc,ngrid,nbloc,
     .            ibloc,igrid,lswpchb,lswpche,
     .            mxpchs,mpchitm,nswpchb,npchitm,
     .            lspchb1,lspchf1,lspchb2,lspchf2,
     .            ipitmb1(litmbeg),ipitmb2(litmbeg),
     .            jpitmb1(litmbeg),jpitmb2(litmbeg),
     .            kpitmb1(litmbeg),kpitmb2(litmbeg),frc(litmbeg),nodes )
#endif
c           
               ncnt2 = 0
               do lpchs=lswpchb(ibloc,igrid)+1,lswpche(ibloc,igrid)
               ncnt2 = ncnt2 + npchitm(lpchs,igrid)
               enddo

c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
#  ifdef BUILD_PVM
c
       write (iwrit,'("PVM not supported for patched grids")')
       call ERREXIT (nodes)
#  else
c
c         n1b  = m1pch1(ibloc,igrid)
c         n2b  = m1pch2(ibloc,igrid)
          n1b  = 1
          n2b  = 1
          if (nodes(ibloc)-1.ne.master) then
          call MPI_Send (tbuf1,ncnt1,MPI_INTEGER,
     .                   nodes(ibloc)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (tbuf2,ncnt1,MPI_INTEGER,
     .                   nodes(ibloc)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (tbuf3,ncnt1,MPI_INTEGER,
     .                   nodes(ibloc)-1,TAG_PATCH,
     .                   mycomm,ierr)
c
          call MPI_Send (tbuf4,ncnt2,MPI_INTEGER,
     .                   nodes(ibloc)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (tbuf5,ncnt2,MPI_INTEGER,
     .                   nodes(ibloc)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (tbuf6,ncnt2,MPI_INTEGER,
     .                   nodes(ibloc)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (tbuf7,ncnt2,MPI_INTEGER,
     .                   nodes(ibloc)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (tbuf8,ncnt2,MPI_INTEGER,
     .                   nodes(ibloc)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (tbuf9,ncnt2,MPI_INTEGER,
     .                   nodes(ibloc)-1,TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (tbuf10,ncnt2,RTYPE,
     .                   nodes(ibloc)-1,TAG_PATCH,
     .                   mycomm,ierr)
          endif
c
#  endif
#endif
      
 104        continue
c           enddo

         enddo
c
c
c
c######################################################################
c node processes
c     receive detailed patching info. from master
c######################################################################
c
#ifdef BUILD_MPI
      end if
#else
      else
#endif


c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
      do igrid = 1,ngrid
         do ibloc=1,nbloc
         if (nodes(ibloc)-1.eq.myrank) then
c
           ncnt1 = 0
           do iface = 1,6
           ncnt1 = ncnt1 + npchcbf(iface,ibloc,igrid)
           enddo
c           
           ncnt2 = 0
           do lpchs=lswpchb(ibloc,igrid)+1,lswpche(ibloc,igrid)
           ncnt2 = ncnt2 + npchitm(lpchs,igrid)
           enddo
c
#  ifdef BUILD_PVM
c
       write (iwrit,'("PVM not supported for patched grids")')
       call ERREXIT (nodes)
#  else
c
c
        if (nodes(ibloc)-1.ne.master) then
        call MPI_Recv (ipatchc(m1pch1(ibloc,igrid)),ncnt1,
     .                 MPI_INTEGER,master,TAG_PATCH,mycomm,
     .                 istat,ierr )
        call MPI_Recv (jpatchc(m1pch1(ibloc,igrid)),ncnt1,
     .                 MPI_INTEGER,master,TAG_PATCH,mycomm,
     .                 istat,ierr )
        call MPI_Recv (kpatchc(m1pch1(ibloc,igrid)),ncnt1,
     .                 MPI_INTEGER,master,TAG_PATCH,mycomm,
     .                 istat,ierr )
c
        call MPI_Recv (ipitmb1(m1pch2(ibloc,igrid)),ncnt2,
     .                 MPI_INTEGER,master,TAG_PATCH,mycomm,
     .                 istat,ierr )
        call MPI_Recv (ipitmb2(m1pch2(ibloc,igrid)),ncnt2,
     .                 MPI_INTEGER,master,TAG_PATCH,mycomm,
     .                 istat,ierr )
        call MPI_Recv (jpitmb1(m1pch2(ibloc,igrid)),ncnt2,
     .                 MPI_INTEGER,master,TAG_PATCH,mycomm,
     .                 istat,ierr )
        call MPI_Recv (jpitmb2(m1pch2(ibloc,igrid)),ncnt2,
     .                 MPI_INTEGER,master,TAG_PATCH,mycomm,
     .                 istat,ierr )
        call MPI_Recv (kpitmb1(m1pch2(ibloc,igrid)),ncnt2,
     .                 MPI_INTEGER,master,TAG_PATCH,mycomm,
     .                 istat,ierr )
        call MPI_Recv (kpitmb2(m1pch2(ibloc,igrid)),ncnt2,
     .                 MPI_INTEGER,master,TAG_PATCH,mycomm,
     .                 istat,ierr )
        call MPI_Recv (frc(m1pch2(ibloc,igrid)),ncnt2,
     .                 RTYPE,master,TAG_PATCH,mycomm,
     .                 istat,ierr )
        end if
c
#  endif
         endif
          enddo
c
       enddo

#  endif
c
#ifdef BUILD_MPI
#else
       endif
#endif
 106   continue
c-----------------------------------------------------------------
c      This completes all processes for patching
c-----------------------------------------------------------------
c
c
c
c######################################################################
c master process
c     read and send grid points to nodes
c######################################################################
c
      if (myrank.eq.master) then
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c     initialize x array

      do nn=1,mxsiz3n
        x(nn) = 0.
        buff(nn) = 0.
      enddo
c
      do ibloc=1,nbloc
        ipts = im(1,ibloc) - imn + 1
        jpts = jm(1,ibloc) - jmn + 1
        kpts = km(1,ibloc) - kmn + 1
        npts = ipts * jpts * kpts * 3
c
        call gridin (igrdfmt,1,1,1,ipts,jpts,kpts,ipts,jpts,kpts,
cjm
     .   buff)
c
#  ifdef BUILD_PVM
        call PVMFpsend (nodes(ibloc),TAG_GRID,buff,npts,RTYPE,ierr)
#  else
cjm        call MPI_Send (x,npts,RTYPE,nodes(ibloc),
      if (myrank.ne.nodes(ibloc)-1) then
        call MPI_Send (buff,npts,RTYPE,nodes(ibloc)-1,
     .                 TAG_GRID,mycomm,ierr)
      else
      do n=0,2
        nn = n * kpts
        ni = n * kmp1(1,ibloc)
        do k=0,kpts-1
          kn = (k + nn) * jpts
          ki = (k + ni + kmn - 1) * jmp1(1,ibloc)
          do j=0,jpts-1
            in = (kn + j) * ipts
            ii = (ki + j + jmn - 1) * imp1(1,ibloc) +
     .           imn - 1 + m3np(1,ibloc) - 1
            do i=1,ipts
              x(ii+i) = buff(in+i)
            enddo
          enddo
        enddo
      enddo
      end if
#  endif
      enddo
c
c######################################################################
c node processes
c     receive grids assigned to this process from host
c     use w array as work space and then
c     fill in the x array with appropiate offsets
c######################################################################
c
#  ifdef BUILD_PVM
      else
#  else
      end if
#  endif
c
      if (myrank.ne.master) then
      do 112 ibloc=1,nbloc
      if (nodes(ibloc)-1.eq.myrank) then
      ipts = im(1,ibloc) - imn + 1
      jpts = jm(1,ibloc) - jmn + 1
      kpts = km(1,ibloc) - kmn + 1
      npts = ipts * jpts * kpts * 3
c
#  ifdef BUILD_PVM
      call PVMFprecv (master,-1,w,npts,RTYPE,
     .                itid,itag,ilen,ierr)
#  else
      call MPI_Recv (w,npts,RTYPE,master,
     .               MPI_ANY_TAG,mycomm,istat,ierr)
      itag = istat(MPI_TAG)
#  endif
      if (itag.ne.TAG_GRID)
     .  call ERREXIT (nodes)
c
      do n=0,2
        nn = n * kpts
        ni = n * kmp1(1,ibloc)
        do k=0,kpts-1
          kn = (k + nn) * jpts
          ki = (k + ni + kmn - 1) * jmp1(1,ibloc)
          do j=0,jpts-1
            in = (kn + j) * ipts
            ii = (ki + j + jmn - 1) * imp1(1,ibloc) +
     .           imn - 1 + m3np(1,ibloc) - 1
            do i=1,ipts
              x(ii+i) = w(in+i)
            enddo
          enddo
        enddo
      enddo
c
      endif
  112 continue
      end if
c
#else	/* IN-CORE version */
c
      do 112 ibloc=1,nbloc
      call gridin (igrdfmt,imn,jmn,kmn,
     .     im  (1,ibloc),jm  (1,ibloc),km  (1,ibloc),
     .     imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
     .     x(m3np(1,ibloc)) )
  112 continue
c
#endif
c
      do 111 ibloc=1,nbloc
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) then
# else
      if (nodes(ibloc).eq.myrank) then
# endif
c
c     generate coarse grids for ibloc
c
      do 105 igrid = 1,ngrid-1
c
      call gridc ( imn,jmn,kmn,
     . im  (igrid  ,ibloc), jm  (igrid  ,ibloc), km  (igrid  ,ibloc),
     . imp1(igrid  ,ibloc), jmp1(igrid  ,ibloc), kmp1(igrid  ,ibloc),
     . im  (igrid+1,ibloc), jm  (igrid+1,ibloc), km  (igrid+1,ibloc),
     . imp1(igrid+1,ibloc), jmp1(igrid+1,ibloc), kmp1(igrid+1,ibloc),
     . x(m3np(igrid,ibloc)), x(m3np(igrid+1,ibloc))          )
c
 105  continue
c
c     compute the volumes and metrics on all grid-levels
c     for ibloc
c
c     the following sets all volumes on all grid levels for
c     the given block to -1 for later detection of volumes
c     which do not get computed
c
      i1 = m1cc(1,ibloc)
      i2 = m1cc(ngrid,ibloc)+imp2(ngrid,ibloc)*
     .     jmp2(ngrid,ibloc)*kmp2(ngrid,ibloc)-1
      do i=i1,i2
        vol(i) = -1.0
      enddo
c
      do 110 igrid = 1,ngrid
c
c     if(igrid.eq.1)
        call volume (imn,jmn,kmn,
     .  im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .  imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .  imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .  vol(m1cc(igrid,ibloc)), x(m3np(igrid,ibloc)),igrid,ibloc)
c
c     if(igrid.gt.1)
c    . call restrvol (igrid,imn,jmn,kmn,
c    .  im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
c    .  imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
c    .  imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
c    .  im  (igrid  ,ibloc), jm  (igrid  ,ibloc), km  (igrid  ,ibloc),
c    .  imp1(igrid  ,ibloc), jmp1(igrid  ,ibloc), kmp1(igrid  ,ibloc),
c    .  imp2(igrid  ,ibloc), jmp2(igrid  ,ibloc), kmp2(igrid  ,ibloc),
c    .    vol(m1cc(igrid-1,ibloc)),  vol(m1cc(igrid,ibloc)) )
c
      call metric (imn,jmn,kmn,
     .  im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .  imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .  x  (m3np(igrid,ibloc)),
     .  six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .  siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .  sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .  skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .  skz(m1np(igrid,ibloc)), 
     .  ri1(m3is(igrid,ibloc)), ri2(m3is(igrid,ibloc)),
     .  rj1(m3js(igrid,ibloc)), rj2(m3js(igrid,ibloc)),
     .  rk1(m3ks(igrid,ibloc)), rk2(m3ks(igrid,ibloc)),ibloc,igrid,
     .  nodes )
c
  110 continue
c
      endif
  111 continue
c
c----------------------------------------------------------------
c
      ilev = nlev
      do igrid=1,ngrid
         do ipchs = 1,mxpchs
           iitmsa(ipchs,igrid) = 1
         enddo
      enddo
c
      do 115 igrid = 1,ngrid
cvn     write (iwrit,'("ilev,igrid,isoln in bcvol",3i3)')
cvn  .  ilev,igrid,isoln(igrid)
c
c     fill-in ghost cell values for cell volumes based on
c     topological information in mapping function 'imap'
c
      call bcvol (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,imp2,jmp2,kmp2, 
     .  m1cc,mgrlev,igrid,igrid,mbloc,nbloc,mxsizc,mx1dwk,
     .  imap(1,1,1,igrid),msegt,msegn,nseg,vol,wk2d,
     .  mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .  mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .  lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .  kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .  ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes)
cvn   if(igrid.lt.ngrid) iitmsa(igrid+1) = iitmsa(igrid)
c
  115 continue
c
c     evaluate the distance to the closest solid wall for the
c     entire field
c
c	print *,'iturb,mstrt: ',iturb,mstrt
c	print *,'calling disdr'
      if (iturb.ge.2 .and. mstrt.lt.1) then
        do iisrf=1,mxdstf
           fwrk(iisrf) = 0.
        enddo
c
        call disdr (imn, jmn, kmn, im,  jm,  km,
     .              imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .              m1cc,m3np,mbloc,nbloc,
     .              mx1dwk,mxdstf,mxdsti,mxsizc,mxsiz3n,
     .              imap,nseg,msegt,msegn,mgrlev,nlev,
     .              x,dtl,wk2d,smin,fwrk,iwrk,mxsurf,nodes   )
      endif
c	print *,'disdr done'
c
c     do 116 igrid=1,ngrid
c     do 116 ibloc=1,nbloc
# if defined BUILD_MPI
c     if (nodes(iblock)-1.eq.myrank) then
# else
c     if (nodes(iblock).eq.myrank) then
# endif
c     call volchk (imn,jmn,kmn,
c    .  im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
c    .  imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
c    .  imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
c    .  vol(m1cc(igrid,ibloc)),x(m3np(igrid,ibloc)),igrid,ibloc)
c     endif
c 116 continue
c     stop
c
#  ifdef BUILD_MPI
#  else
      endif
#  endif
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
#  ifdef BUILD_PVM

#  else
#  endif
#endif
c---------------------------------------------------------------------
c
c     at this point grid-coordinates,cell volumes metrics and
c     normals for all grid levels are available on all blocks
c
c**********************************************************
c     initialize flow quantities and convergence histories
c     mstrt .eq. 0 , free-stream
c     mstrt .ne. 0 , from restart file
c
c######################################################################
c master process
c     if this is a restart, read the restart data and send to
c     the appropriate node process
c######################################################################
c
c
c     the following initializes all flow quantities on all grids
c
      do ibloc=1,nbloc
# if defined BUILD_MPI
        if (nodes(ibloc)-1.eq.myrank) then
#  else
        if (nodes(ibloc).eq.myrank) then
#  endif
          i1 = m5cc(1,ibloc)
          i2 = m5cc(ngrid,ibloc)+imp2(ngrid,ibloc)*
     .         jmp2(ngrid,ibloc)*kmp2(ngrid,ibloc)*5-1
          do i=i1,i2
            w(i) = 0.0
          enddo
          i1 = m1cc(1,ibloc)
          i2 = m1cc(ngrid,ibloc)+imp2(ngrid,ibloc)*
     .         jmp2(ngrid,ibloc)*kmp2(ngrid,ibloc)-1
          do i=i1,i2
            p(i)     = 0.0
            eomu(i)  = 0.0
            if(iturb.ge.2) turv1(i) = 0.0
            if(iturb.gt.2) turv2(i) = 0.0
          enddo
        endif
      enddo
c
      if (myrank.eq.master) then
c	print *,'all grid stuff done'
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
      if (mstrt.ne.0) then
        igrid = isoln(1)
        do ibloc = 1,nbloc
c
#  ifdef BUILD_PVM
          call start ( mstrt,imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w, p, eomu, turv1, turv2, smin,
     .    wold,tv1old,tv2old,ntorder,iunsteady,jturb,nnorder )
#  else

        if (myrank.eq.nodes(ibloc)-1) then
          call start ( mstrt,imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .    w(m5cc(igrid,ibloc)), p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)), turv1(m1cc(igrid,ibloc)),
     .    turv2(m1cc(igrid,ibloc)),smin(m1cc(igrid,ibloc)) ,
     .    wold(m5fgt(1,ibloc)),tv1old(m1fgt(1,ibloc)),
     .    tv2old(m1fgt(1,ibloc)),ntorder,iunsteady,jturb,nnorder )
        else
          call start ( mstrt,imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .    buffw, buffp, buffe, buff1, buff2, buffs,
     .    buffwo,bufftv1,bufftv2,ntorder,iunsteady,jturb,nnorder )
        end if
 
#  endif
c
          ncnt = imp2(igrid,ibloc) * jmp2(igrid,ibloc) *
     .           kmp2(igrid,ibloc)
#  ifdef BUILD_PVM
          call PVMFpsend (nodes(ibloc),TAG_RSTRT,
     .                    w,ncnt*5,RTYPE,ierr)
          call PVMFinitsend (PvmDataInPlace,ierr)
          call PVMFpack (RTYPE,p,ncnt,1,ierr)
          call PVMFpack (RTYPE,eomu,ncnt,1,ierr)
          if (iturb.eq.2 .or. iturb.eq.3) then
            call PVMFpack (RTYPE,smin,ncnt,1,ierr)
            call PVMFpack (RTYPE,turv1,ncnt,1,ierr)
            if (iturb.eq.3)
     .        call PVMFpack (RTYPE,turv2,ncnt,1,ierr)
          endif
          call PVMFsend (nodes(ibloc),TAG_RSTRT,ierr)
#  else
        if (nodes(ibloc)-1.ne.myrank) then
          call MPI_Send (buffw,ncnt*5,RTYPE,
     .                   nodes(ibloc)-1,TAG_RSTRT,mycomm,ierr)
          call MPI_Send (buffp,ncnt,RTYPE,
     .                   nodes(ibloc)-1,TAG_RSTRT,mycomm,ierr)
          call MPI_Send (buffe,ncnt,RTYPE,
     .                   nodes(ibloc)-1,TAG_RSTRT,mycomm,ierr)
          if (iturb.eq.2 .or. iturb.eq.3) then
            call MPI_Send (buffs,ncnt,RTYPE,
     .                     nodes(ibloc)-1,TAG_RSTRT,mycomm,ierr)
            call MPI_Send (buff1,ncnt,RTYPE,
     .                     nodes(ibloc)-1,TAG_RSTRT,mycomm,ierr)
            if (iturb.eq.3)
     .        call MPI_Send (buff2,ncnt,RTYPE,
     .                       nodes(ibloc)-1,TAG_RSTRT,mycomm,ierr)
          end if
c
c         logic for time-dependent quantities (old time values)
c
          if(ntorder.gt.0) then
            call MPI_Send (totptim,1,RTYPE,
     .                     nodes(ibloc)-1,TAG_RSTRT,mycomm,ierr)
            call MPI_Send (buffwo,ncnt*5*ntorder,RTYPE,
     .                     nodes(ibloc)-1,TAG_RSTRT,mycomm,ierr)
            if (iturb.eq.2 .or. iturb.eq.3) then
              call MPI_Send (bufftv1,ncnt*ntorder,RTYPE,
     .                       nodes(ibloc)-1,TAG_RSTRT,mycomm,ierr)
              if (iturb.eq.3)
     .        call MPI_Send (bufftv2,ncnt*ntorder,RTYPE,
     .                       nodes(ibloc)-1,TAG_RSTRT,mycomm,ierr)
            endif
          endif
        end if
#  endif
        enddo
c
      endif
c
c######################################################################
c node processes
c     initialize flowfield and receive restart data from root
c######################################################################
c
#ifdef BUILD_PVM
      else
#endif
#ifdef BUILD_MPI
      end if
#endif
#endif
      igrid = isoln(1)
ccvn    write (iwrit,'("ilev,igrid,isoln in init",3i3)')
ccvn .  ilev,igrid,isoln(igrid)
      do 120 ibloc = 1,nbloc
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) then
# else
      if (nodes(ibloc).eq.myrank) then
# endif
c
c       start from free-stream conditions
c
        if (mstrt.eq.0) then
cvn       totptim = 0.
          call init (imn,jmn,kmn,
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)), p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)), smin(m1cc(igrid,ibloc)),
     .    x(m3np(igrid,ibloc)),
     .    turv1(m1cc(igrid,ibloc)),turv2(m1cc(igrid,ibloc)),ibloc ,
     .    wold(m5fgt(1,ibloc)),tv1old(m1fgt(1,ibloc)),
     .    tv2old(m1fgt(1,ibloc)),ntorder,iunsteady,igrid,dtphy)
c
c       restart from previous solution
c
        else
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c       receive restart data from master process
c
          ncnt = imp2(igrid,ibloc) * jmp2(igrid,ibloc) *
     .           kmp2(igrid,ibloc)
#  ifdef BUILD_PVM
          call PVMFprecv (master,-1,w(m5cc(igrid,ibloc)),ncnt*5,
     .                    RTYPE,itid,itag,ilen,ierr)
#  else
c--------1---------2---------3---------4---------5---------6---------7-2
        if (myrank.ne.master) then
          call MPI_Recv (w(m5cc(igrid,ibloc)),ncnt*5,
     .                   RTYPE,master,MPI_ANY_TAG,
     .                   mycomm,istat,ierr)
          itag = istat(MPI_TAG)
        end if
#  endif
        if (myrank.ne.master) then
          if (itag.ne.TAG_RSTRT)
     .      call ERREXIT (nodes)
        end if
#  ifdef BUILD_PVM
          call PVMFrecv (master,TAG_RSTRT,ierr)
          call PVMFunpack (RTYPE,p(m1cc(igrid,ibloc)),
     .                     ncnt,1,ierr)
          call PVMFunpack (RTYPE,eomu(m1cc(igrid,ibloc)),
     .                     ncnt,1,ierr)
          if (iturb.eq.2 .or. iturb.eq.3) then
            call PVMFunpack (RTYPE,smin(m1cc(igrid,ibloc)),
     .                       ncnt,1,ierr)
            call PVMFunpack (RTYPE,turv1(m1cc(igrid,ibloc)),
     .                       ncnt,1,ierr)
            if (iturb.eq.3)
     .        call PVMFunpack (RTYPE,turv2(m1cc(igrid,ibloc)),
     .                         ncnt,1,ierr)
          endif
#  else
        if (myrank.ne.master) then
          call MPI_Recv (p(m1cc(igrid,ibloc)),ncnt,
     .                   RTYPE,master,TAG_RSTRT,
     .                   mycomm,istat,ierr)
          call MPI_Recv (eomu(m1cc(igrid,ibloc)),ncnt,
     .                   RTYPE,master,TAG_RSTRT,
     .                   mycomm,istat,ierr)
          if (iturb.eq.2 .or. iturb.eq.3) then
            call MPI_Recv (smin(m1cc(igrid,ibloc)),ncnt,
     .                     RTYPE,master,TAG_RSTRT,
     .                     mycomm,istat,ierr)
            call MPI_Recv (turv1(m1cc(igrid,ibloc)),ncnt,
     .                     RTYPE,master,TAG_RSTRT,
     .                     mycomm,istat,ierr)
            if (iturb.eq.3)
     .        call MPI_Recv (turv2(m1cc(igrid,ibloc)),ncnt,
     .                       RTYPE,master,TAG_RSTRT,
     .                       mycomm,istat,ierr)
          end if
c
c         logic for time-dependent quantities (old time values)
c
          if (ntorder.gt.0) then
            call MPI_Recv (totptim,1,
     .                     RTYPE,master,TAG_RSTRT,
     .                     mycomm,istat,ierr)
            call MPI_Recv (wold(m5fgt(1,ibloc)),ncnt*5*ntorder,
     .                     RTYPE,master,TAG_RSTRT,
     .                     mycomm,istat,ierr)
            if (iturb.eq.2 .or. iturb.eq.3) then
              call MPI_Recv (tv1old(m1fgt(1,ibloc)),ncnt*ntorder,
     .                       RTYPE,master,TAG_RSTRT,
     .                       mycomm,istat,ierr)
              if (iturb.eq.3)
     .        call MPI_Recv (tv2old(m1fgt(1,ibloc)),ncnt*ntorder,
     .                       RTYPE,master,TAG_RSTRT,
     .                       mycomm,istat,ierr)
             endif
           endif
        endif
#  endif
#else	/* IN-CORE version */
c
          call start ( mstrt,imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .    w(m5cc(igrid,ibloc)), p(m1cc(igrid,ibloc)),
     .    eomu(m1cc(igrid,ibloc)),turv1(m1cc(igrid,ibloc)),
     .    turv2(m1cc(igrid,ibloc)),smin(m1cc(igrid,ibloc)),
     .    wold(m5fgt(1,ibloc)),tv1old(m1fgt(1,ibloc)),
     .    tv2old(m1fgt(1,ibloc)),ntorder,iunsteady,jturb,nnorder )
c
#endif
        endif
c
      endif
 120  continue
c
cvn       write (iwrit,'(" ready to call bcflow")')

      call bcflow (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  m1np,m2np,m3np,m4np,m5np,m1cc,m2cc,m3cc,m4cc,m5cc,m5cg,
     .  m1is,m1js,m1ks,m3is,m3js,m3ks,mgrlev,igrid,igrid,mbloc,nbloc,
     .  mxsizc,mxsiz5c,mxsizn,mxsiz3n,mx1dwk,
     .  imap(1,1,1,igrid),msegt,msegn,nseg,w,p,eomu,turv1,turv2,vol,x,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .  ri1,ri2,rj1,rj2,rk1,rk2,
     .  tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .  uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
     .  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .  wk2d,const1,const2,const3,const4,const5,
     .  mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .  mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .  lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .  kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .  ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,iunsteady,totptim,smin)

ccvn      write (iwrit,'(" finished bcflow")')
c
c     fill-in the boundary information on eddy-viscosity
c     and turbulent quantities (for field equations)
c
cvn       write (iwrit,'(" ready to call bcturb")')

      call bcturb (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,
     .  imp2,jmp2,kmp2,m1cc,m5cc,mgrlev,igrid,igrid,mbloc,nbloc,
     .  mxsizc,mxsiz5c,mx1dwk,imap(1,1,1,igrid),msegt,msegn,nseg,
     .               wk2d(1,1),eomu,turv1,turv2,w,p,smin,
     .  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .  tnuwfi1,tnuwfi2,tnuwfj1,tnuwfj2,tnuwfk1,tnuwfk2,m1is,m1js,m1ks,
     .  mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .  mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .  lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .  kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche ,
     .  ipitmbs,jpitmbs,kpitmbs ,iitmsa,nodes)
c
ccvn       write (iwrit,'(" finished bcturb",2x,i5)') igrid
#ifdef BUILD_MPI
#else
      endif
#endif
c
c
c######################################################################
c all processes
c     initialization of flow quantities is now complete
c######################################################################
c
c     write out the initial solution to the plotting files
c
        if(mplot .gt. 0. and. iunsteady.gt.0) then
           call savtim(imn,jmn,kmn,
     .                 im,jm,km,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .                 mxsizc,mxsiz5c,mxsiz3n,mgrlev,mbloc,msegt,msegn,
     .                 m5cc,m1cc,m3np,imap,nbloc,nseg,
# if defined BUILD_MPI
     .                 mxszfg,buffw,buffp,buff,
# endif
     .                 x,w,p,delw,alpha,totptim,1,nodes)
           call savdata(ibeggrp,iendgrp,iskpgrp,
     .                  jbeggrp,jendgrp,jskpgrp,
     .                  kbeggrp,kendgrp,kskpgrp,
     .                  iblkgrp,ngroup,
     .                  imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2,
     .                  mxsizc,mxsiz5c,mxsiz3n,mgrlev,mbloc,msegt,msegn,
     .                  m5cc,m1cc,m3np,imap,nbloc,nseg,
# if defined BUILD_MPI
     .                  mxszfg,buffw,buffp,buff,
# endif
     .                  x,w,p,delw,alpha,totptim,1,nodes)
        end if
c
c The master process (myrank = 0) executes only the outer-most
c 2 loops so that it can synchronize with the nodes, receive
c the convergence results from the nodes and print the results
c
c The node processes execute all the loops, however calculations
c are only performed for the blocks which are local to the node
c process (nodes(ibloc) = myrank)
c
c jbloc has been introduced for some of the loops over the blocks
c below. This variable counts the blocks that are actually on a
c given node and is used in place of ibloc as an argument to the
c routines:
c    restrct, restrev, prolong, and chkconv
c It is probably unneccessary in the first 3 routines, since it
c is only used to change the value of iftf, which is no longer
c used. It is needed by chkconv, however, since initialization
c is done in that routine when ibloc = 1.
c
c     write heading for start of calculations
c
      if (iwrit.gt.0) then
c
      write (iwrit,'(//,
     .  " Steady-state Solution on",i4," blocks")') nbloc
      if (nlev.gt.1) then
        write (iwrit,'(" Full Multigrid using",i3," grid levels")')
     .    nlev
        if (mgtype.eq.1) then
          write (iwrit,'(" and V-Cycles on each grid level")')
        else
          write (iwrit,'(" and W-Cycles on each grid level")')
        endif
      else
        if (mgtype.eq.1) then
          write (iwrit,'(" Multigrid with V-Cycles")')
        else if (mgtype.eq.2) then
          write (iwrit,'(" Multigrid with W-Cycles")')
        endif
      endif
c
      write (iwrit,'(/,4x,"    rm          reyl        alpha   ",
     .                   "    yaw         roll    ")')
      write (iwrit,'(1x,7f12.5)') rm,reyl,al,yaw,roll
c
      endif
c
c     initialize timings for grids
c
      do n=1,3
        tsum(n) = 0.
      enddo
      call cputim (tim)
      wtot = totw
c
c******************  begin outer fmg loop  ******************
c
      do 1000 ilev = 1,nlev
c
c     write headings for multi-grid solution on this mesh
c
      if (iwrit.gt.0) then
c
      write (iwrit,'(//)')
      if (nlev.gt.1)
     .  write (iwrit,'(" Multigrid level  =",i4)') ilev
      write (iwrit,'(" Number of Cycles =",i4)') mcyc(ilev)
      write (iwrit,'(" Solution Grid    =",i4)') isoln(ilev)
      write (iwrit,'(/," ibegr  icoars  ifine")')
      do irleg=1,nrleg(ilev)
        write (iwrit,'(2x,i2,6x,i2,5x,i2)')
     .    ibegr(irleg,ilev),icoars(irleg,ilev),ifine(irleg,ilev)
      enddo
      write (iwrit,'(/," cycle  work  ",
     .              " max dr/dt i   j  k  n avg dr/dt",
     .              " max enthpy i   j   k  n max_rate",
     .              " avg_rate  cl       cdt         cdv        cmz")')
c    .                " max enthpy  i   j   k   n avg enthpy",
c    .                " n sup     cl      cdt      cdv      cmz")')
      if (mstrt.ne.0 .and. myrank.eq.master)
     .  write (iwrit,651) nres,wcount(nres),rtrmsg(nres),hrmsg(nres),
     .    nsupg(nres),clg(nres),cdtg(nres),cdvg(nres),cmzg(nres)
c
      endif
c
  650 format(i4,f8.2,e10.3,2i4,2i3,2e10.3,2i4,2i3,2f6.3,
     .       4e12.4)
c 650 format(1x,i4,1x,f10.3,e10.3,4i4,2e10.3,4i4,e10.3,5f9.5)
c 650 format(1x,i4,f10.3,e10.3,4i4,2e10.3,4i4,e10.3,i7,4f9.5)
  651 format(1x,i4,f10.3,5x,"restart",7x,e10.3,7x,"restart",7x,
     .       e10.3,i7,4e12.4)
c
c######################################################################
c node processes
c     before proceeding to the next fmg level, use interpolation
c     to get starting solution on the fine mesh
c######################################################################
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
cjm      if (myrank.ne.master) then
#endif
      if (ilev.gt.1) then
c
        igrid = isoln(ilev)
        do 130 ibloc = 1,nbloc
c
# if defined BUILD_MPI
        if (nodes(ibloc)-1.eq.myrank)
# else
        if (nodes(ibloc).eq.myrank)
# endif
     .  call refine ( imn, jmn, kmn,
     .  imp2(igrid,  ibloc),jmp2(igrid,  ibloc),kmp2(igrid,  ibloc),
     .  imp2(igrid+1,ibloc),jmp2(igrid+1,ibloc),kmp2(igrid+1,ibloc),
     .               w(m5cc(igrid,ibloc)), w(m5cc(igrid+1,ibloc)),
     .               p(m1cc(igrid,ibloc)), p(m1cc(igrid+1,ibloc)),
     .         eomu(m1cc(igrid,ibloc)), eomu(m1cc(igrid+1,ibloc)),
     .        turv1(m1cc(igrid,ibloc)),turv1(m1cc(igrid+1,ibloc)),
     .        turv2(m1cc(igrid,ibloc)),turv2(m1cc(igrid+1,ibloc)))
  130   continue
ccvn     write (iwrit,'(" finish refine:ig",2i5)') igrid,isoln(ilev)
c
c----------------------------------------------------------------
c   apply the boundary conditions on conserved 
c   variables 'w'  at this grid level
c
      call bcflow (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  m1np,m2np,m3np,m4np,m5np,m1cc,m2cc,m3cc,m4cc,m5cc,m5cg,
     .  m1is,m1js,m1ks,m3is,m3js,m3ks,mgrlev,igrid,igrid,mbloc,nbloc,
     .  mxsizc,mxsiz5c,mxsizn,mxsiz3n,mx1dwk,
     .  imap(1,1,1,igrid),msegt,msegn,nseg,w,p,eomu,turv1,turv2,vol,x,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .  ri1,ri2,rj1,rj2,rk1,rk2,
     .  tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .  uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
     .  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .  wk2d,const1,const2,const3,const4,const5,
     .  mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .  mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .  lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .  kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .  ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,iunsteady,totptim,smin)
ccvn       write (iwrit,'(" finished bcflow2")')
c
c     fill-in the boundary information on eddy-viscosity
c     and turbulent quantities (for field equations)
c
ccvn       write (iwrit,'(" ready to call bcturb",i5)') isoln(ilev)

      call bcturb (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,
     .  imp2,jmp2,kmp2,m1cc,m5cc,mgrlev,igrid,isoln(ilev),mbloc,nbloc,
     .  mxsizc,mxsiz5c,mx1dwk,imap(1,1,1,igrid),msegt,msegn,nseg,
     .               wk2d(1,1),eomu,turv1,turv2,w,p,smin,
     .  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .  tnuwfi1,tnuwfi2,tnuwfj1,tnuwfj2,tnuwfk1,tnuwfk2,m1is,m1js,m1ks,
     .  mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .  mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .  lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .  kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche ,
     .  ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes )

      endif
c
c
c######################################################################
c all processes
c######################################################################
c
      work  = 0.
      ires0 = nres +1
c
c
c**************  begin physical time stepping loop *************
c
c
      do 1005 itime = 1,nsteps(ilev)
c
c     increment the physical time (stored in totptim)
c
c     if(myrank.eq.master) totptim = totptim + dtphy
      if(ntorder .ge. 1) then
c
          totptim = totptim + dtphy
        if(myrank.eq.master) then
cvn       totptim = totptim + dtphy
          write (irhtall  ,'("#  cycle   work   ",
     .                " max dr/dt  i   j   k   n avg dr/dt",
     .                " max enthpy  i  j    k n max_rate",
     .                " avg_rate   cl      cdt      cdv       cmz")')
          if (itime .eq. 1 .and. ilev .eq. 1)
     .    write (irhtphy ,'("# time step",3x,"time",13x,"rumsey t",
     .        7x,"cl",15x,"cdt",14x,"cdv",14x,"cmz")')
       endif
c#ifdef RKPHY
      if(irkcase.gt.0.or.mebdf.gt.0) then
        if(myrank.eq.master) then
          if (itime .eq. 1 .and. ilev .eq. 1)
     .    write (ipltfil ,'("# time step",3x,"time",13x,"rumsey t",
     .        7x,"erkmx1",12x,"erkmx2",11x,"erkmx3",11x,"erkmx5",
     .        9x,"erkrms1",11x,"erkrms2",10x,"erkrms3",10x,"erkrms5")')
        endif
      endif
c#endif
c
c       code for accelerated time-accurate calculations
c
c
        do 1002 ibloc=1,nbloc
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) then
# else
      if (nodes(ibloc).eq.myrank) then
# endif
          igrid=1
cvn       if (itime .gt. 1) then
cvn       call equw2w (imn,jmn,kmn,
cvn  .         imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
cvn  .         imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
cvn  .         wold(m5fgt(iorder-1,ibloc)),
cvn  .         w(m5cc(igrid,ibloc)), ws(m5fg(ibloc)) ,5)
cvn       else
          call equw1w (imn,jmn,kmn,
     .         imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .         imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .         w(m5cc(igrid,ibloc)), ws(m5fg(ibloc)) ,5)
cvn       endif
c
          if(ntorder.gt.1) then
            do iorder=ntorder,2,-1
            call equw1w (imn,jmn,kmn,
     .           imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
     .           imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
     .           wold(m5fgt(iorder-1,ibloc)),
     .           wold(m5fgt(iorder  ,ibloc)) ,5)
            if(iturb.ge.2)
     .      call equw1w (imn,jmn,kmn,
     .           imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
     .           imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
     .           tv1old(m1fgt(iorder-1,ibloc)),
     .           tv1old(m1fgt(iorder  ,ibloc)) ,1)
            if(iturb.eq.3)
     .      call equw1w (imn,jmn,kmn,
     .           imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
     .           imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
     .           tv2old(m1fgt(iorder-1,ibloc)),
     .           tv2old(m1fgt(iorder  ,ibloc)) ,1)
            enddo
          endif
c
         if(ntorder.ge.1) then
          call equw1w (imn,jmn,kmn,
     .         imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
     .         imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
     .         w(m5cc(1,ibloc)),  wold(m5fgt(1,ibloc)) ,5)
          if(iturb.ge.2) 
     .    call equw1w (imn,jmn,kmn,
     .         imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
     .         imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
     .         turv1(m1cc(1,ibloc)),tv1old(m1fgt(1,ibloc)) ,1)
          if(iturb.eq.3) 
     .    call equw1w (imn,jmn,kmn,
     .         imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
     .         imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
     .         turv2(m1cc(1,ibloc)),tv2old(m1fgt(1,ibloc)) ,1)
#ifdef RKPHY
cvn      if(irkcase.ge.1.and.itime.gt.1)
cvn  .   call equw1w (imn,jmn,kmn,
cvn  .         imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
cvn  .         imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
cvn  .         delw, delwork(m5fgtrk(1,ibloc)),5 )
#endif
        endif
c
c       end the test for local block  for do 1002
c
        endif
 1002 continue
c
c
        if (myrank.eq.master)
     .  ttrum = totptim * sqrt(gamma)
c       write (iwrit,'("0",10x,"time step = ",i5,5x,"total time = ",
c    .               1pe12.5," rumsey time = ",e12.5)')
c    .       itime,totptim,ttrum
c
c     end the test on ntorder to fill-in time-dependent arrays
c     started before the do 1002
      endif
c
      if (ntorder.gt.0 .and. itime.gt.1) then
c
c     scale ubmin2, a preconditioning parameter by maximum velocity
c            amachg is maximum Mach number over field
c
        ubmin2    = ubmin*ubmin*amachg
        ubcon2    = ubcon*ubcon*amachg
        ammax     = max(ammax,amachg)
        if (iwrit.gt.0)
     .    write (iwrit,'("time =",i6," maxMach = ",e10.3,
     .      " global maxMach = ",e10.3," ubmin2 = ",e10.3,
     .      " ubcon2 = ",e10.3)')
     .    itime,sqrt(amachg),sqrt(ammax),ubmin2,ubcon2
      endif
c
c
#ifdef  RKPHY
c    Initialize flux (delwork) for physical RK-scheme
c
cvn   if (irkcase.ge.1 .and. itime.eq.1) then
      if (irkcase.ge.1) then
      igrid = 1
      cfl   = cflf
      irkstg= 1

      call fluxwrk (imn, jmn, kmn, im,  jm,  km,
     .            imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .            m1np,m2np,m3np,m4np,m5np, m1cc,m2cc,m3cc,m4cc,m5cc,
     .            m1fg,m4fg,m5fg,m5cg,m1is,m1js,m1ks,m3is,m3js,m3ks,
     .            igrid,isoln(ilev),nbloc,mres,
     .            mx1dwk,mx3dwk,mxsizc,mxsiz3c,mxsiz4c,mxsiz5c,
     .            mxszfg,mxsz4fg,mxsz5fg,mxsz5cg,mxsizn,mxsiz3n,
     .            imap,nseg,msegt,msegn,mbloc,mgrlev,
     .            w,p,ws,wr,delw,fw,fv,eomu,turv1,turv2,vol,x,fqs,
     .            six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .            ri1,ri2,rj1,rj2,rk1,rk2,
     .            tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .            dtl,dti,dtj,dtk,dtmin,ratioij,ratiojk,eprs,
     .            wk1d,wk2d,wk1d5,rtrmsg,
     .            nitr,iforfn,ilast,ncyc,ibcup,ilev,nres,
     .            ivisb,const1,const2,const3,const4,const5,mstrt,
     .            smin,iturbb,itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,
     .            mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .            mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .            lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .            kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .            ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,
     .            dtvi,dtvj,dtvk,fbeta2,wold,tv1old,tv2old,
     .            m1fgt,m5fgt,iturv1,iturv2,mtorder,ntorder,iunsteady,
     .            dtphy,m1fgtrk,m5fgtrk,irkstg,nrkstg,mrkstg,
     .            rkphysa,rkphyse,delwork,dtv1ork,dtv2ork )
c
c     put the computed flux into delwork array
c
      do ibloc=1,nbloc
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) then
# else
      if (nodes(ibloc).eq.myrank) then
# endif

cvn      call equw1w (imn,jmn,kmn,
cvn  .         imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
cvn  .         imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
cvn  .         delw, delwork(m5fgtrk(1,ibloc)),5 )
c 
c      end of local block if-test
       endif
      enddo
c
c     end initialization of physical R-K first stage
      endif
#endif

c
#ifdef  RKPHY
ch************* begin physical RK-loop
c     This do loop is outside of do 1010 loop on subiterations
c     in pseudo-time and ends right after "1010 continue"
c
      do 2004 irkstg = 2, nrkstg
# else
      do 2004 irkstg = 1, nrkstg
# endif
c
      if(myrank.eq.master)
     .  write(iwrit,'("# irkstg = ",i5)') irkstg
c
#ifdef  RKPHY
c     utilize a better guess on "w" from delwold arrays from predictrk
c
      do ibloc=1,nbloc
# if defined BUILD_MPI
        if (nodes(ibloc)-1.eq.myrank) then
# else
        if (nodes(ibloc).eq.myrank) then
# endif
cvn       if(irkcase.eq.5.and.itime.gt.1) then
c     deactivate next call
cvn       if(irkcase.eq.5.and.itime.gt.1.and.irkstg.le.2) then
cvn       call equw2w (imn,jmn,kmn,
cvn  .         im  (1,ibloc),jm  (1,ibloc),km  (1,ibloc),
cvn  .         imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
cvn  .         imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
cvn  .         delwork(m5fgtrk(irkstg,ibloc)),
cvn  .         w(m5cc(1,ibloc)), 5)
cvn       endif
cvn     endif
       endif
      enddo
c
c     update the b.c.'s on new estimates of "w"
c
c     igrid = 1
c     call bcflow (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
c    .  m1np,m2np,m3np,m4np,m5np,m1cc,m2cc,m3cc,m4cc,m5cc,m5cg,
ccvn .  m5cg,m3is,m3js,m3ks,mgrlev,igrid,isoln(1),mbloc,nbloc,
c    .  m1is,m1js,m1ks,m3is,m3js,m3ks,mgrlev,igrid,igrid,mbloc,nbloc,
c    .  mxsizc,mxsiz5c,mxsizn,mxsiz3n,mx1dwk,
c    .  imap(1,1,1,igrid),msegt,msegn,nseg,w,p,eomu,turv1,turv2,vol,x,
c    .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
c    .  ri1,ri2,rj1,rj2,rk1,rk2,
c    .  tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
c    .  uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
c    .  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
c    .  wk2d,const1,const2,const3,const4,const5,
c    .  mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
c    .  mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
c    .  lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
c    .  kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
c    .  ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,iunsteady,totptim,smin)

c     end of b.c. setup for RKPHY with new estimates on stage values
#endif
c
      if(mebdf.gt.0.and.irkstg.eq.1.and.itime.gt.1) then
c
c     utilize a better guess on"w" from MEBDF stage II
c
        do ibloc = 1,nbloc
c
# if defined BUILD_MPI
          if (nodes(ibloc)-1.eq.myrank) then
# else
          if (nodes(ibloc).eq.myrank) then
# endif
            igrid  = 1

            if(irkstg.eq.1)
     .      call equw1w (imn,jmn,kmn,
     .            imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .            imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .            wbar(m5fgtrk(irkstg+1,ibloc)),w(m5cc(igrid,ibloc)),5 )
c         update pressure
            if(irkstg.eq.1)
     .      call press (imn,jmn,kmn,
     .        im  (igrid,ibloc),jm  (igrid,ibloc),km  (igrid,ibloc),
     .        imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .        imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .        w(m5cc(igrid,ibloc)),               p(m1cc(igrid,ibloc)) )
c
            if(irkstg.eq.1.and.iturb.ge.2)
     .      call equw1w (imn,jmn,kmn,
     .            imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .            imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .      tv1bar(m1fgtrk(irkstg+1,ibloc)),turv1(m1cc(igrid,ibloc)),1 )
c
            if(irkstg.eq.1.and.iturb.eq.3)
     .      call equw1w (imn,jmn,kmn,
     .            imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .            imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .      tv2bar(m1fgtrk(irkstg+1,ibloc)),turv2(m1cc(igrid,ibloc)),1 )
c
          endif
        enddo
c       end of do loop on blocks for updating of "w"
c
c       update the b.c.'s on new estimates of "w"
c
        call bcflow (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  m1np,m2np,m3np,m4np,m5np,m1cc,m2cc,m3cc,m4cc,m5cc,m5cg,
ccvn .  m5cg,m3is,m3js,m3ks,mgrlev,igrid,isoln(1),mbloc,nbloc,
     .  m1is,m1js,m1ks,m3is,m3js,m3ks,mgrlev,igrid,igrid,mbloc,nbloc,
     .  mxsizc,mxsiz5c,mxsizn,mxsiz3n,mx1dwk,
     .  imap(1,1,1,igrid),msegt,msegn,nseg,w,p,eomu,turv1,turv2,vol,x,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .  ri1,ri2,rj1,rj2,rk1,rk2,
     .  tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .  uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
     .  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .  wk2d,const1,const2,const3,const4,const5,
     .  mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .  mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .  lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .  kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .  ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,iunsteady,totptim,smin)

c       end of b.c. setup for new estimates on "w" for MEBDF4
c
      endif
c
c**************  begin multigrid cycle iteration loop *************
c
      do 1010 icyc = 1,mcyc(ilev)
c     call flush (iwrit)
      ncyc  = ncyc +1
      if (nres.lt.mres) nres  = nres +1
      iftf  = 1
      ibcup = mbcupd(ilev)
c
c######################################################################
c node processes
c######################################################################
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
cjm      if (myrank.ne.master) then
#endif
c
c     the do loop 1010 is the outer most loop for
c     multigrid for a given fmg level
c     this loop wraps around the do loops on
c     mg legs, restriction and prolongation operators
c     in this do loop, specified (from input) number of
c     complete mg cycles are performed for a given level
c     of fmg cycle
c
c--------------  loop through the legs of mg cycle
c
      do 1020 irleg = 1,nrleg(ilev)
c
c     a complete multigrid leg consist of a restriction
c     and a prolongation leg
c     the starting and ending grid levels for these legs
c     were computed in "setseq"
c
c--------------  begin restriction loop    ------------------
c
      do 500 igrid = ibegr(irleg,ilev),icoars(irleg,ilev)
c
c   use cflf if on finest grid level; else use cflc
c
      if (igrid.eq.isoln(ilev)) then
        cfl = cflf
        hm  = hmf
      else
        cfl = cflc
        hm  = hmc
      endif
c     if (igrid.gt.1) cfl = cflc
      if (icord.ne.0 .and. icyc.le.9) cfl = min(cfl,40.)
c
c
c   skip the forcing function computation for the
c   finest mesh (on the given fmg level), since
c   the forcing function
c
c    h
c   f  = 0 ,  for the finest mesh
c
c
        if (igrid .eq. ibegr(irleg,ilev) ) go to 200
c
c              h   h
c   compute   L  (u  ) , for coarser grids
c   which consists of convective, dissipative and physical fluxes
c
        jbloc = 0
        do 210 ibloc = 1,nbloc
# if defined BUILD_MPI
        if (nodes(ibloc)-1.eq.myrank) then
# else
        if (nodes(ibloc).eq.myrank) then
# endif
        jbloc = jbloc + 1
c
c     set up fqs array, which determines whether a given boundary
c     is a solid surface or not
c
        do 202 ia=1,6
        do 202 n2=1,mx1dwk
        do 202 n1=1,mx1dwk
        fqs(n1,n2,ia) = 1.0
  202   continue
c
        do 206 ns=1,nseg(ibloc)
        nbctype = imap(1,ns,ibloc,igrid-1)
c       if(nbctype.ge.2.and.nbctype.le.3) then
        if(nbctype.eq.2.or.nbctype.eq.102) then
          nface = imap(2,ns,ibloc,igrid-1)
          nbeg1 = imap(3,ns,ibloc,igrid-1) + 1
          nend1 = imap(4,ns,ibloc,igrid-1)
          nbeg2 = imap(5,ns,ibloc,igrid-1) + 1
          nend2 = imap(6,ns,ibloc,igrid-1)
          do 204 n2=nbeg2,nend2
          do 204 n1=nbeg1,nend1
          fqs(n1,n2,nface) = 0.0
  204     continue
        end if
  206   continue
c
c      set up time step-related terms for this block
c
      call blkstp (imn,jmn,kmn,
     .     im  (igrid-1,ibloc),jm  (igrid-1,ibloc),km  (igrid-1,ibloc),
     .     imp1(igrid-1,ibloc),jmp1(igrid-1,ibloc),kmp1(igrid-1,ibloc),
     .     imp2(igrid-1,ibloc),jmp2(igrid-1,ibloc),kmp2(igrid-1,ibloc),
     . vol(m1cc(igrid-1,ibloc)),   dti(m1fg(ibloc)),   dtj(m1fg(ibloc)),
     .         dtk(m1fg(ibloc)), dtvi(m1fg(ibloc)), dtvj(m1fg(ibloc)),
     .         dtvk(m1fg(ibloc)), dtl,ratioij,ratiojk,eprs  ,
#ifdef  RKPHY
     .       ntorder, dtphy ,rkphysa(irkstg,irkstg) )
# else
     .       ntorder, dtphy )
# endif
c
c
c   evaluate artificial and physical dissipation terms on
c   a grid 1-level finer than 'igrid'
c
c
        if ( (igrid - isoln(ilev)) .eq. 1 ) then
c
ccc      if(icyc.eq.1)
          call filtf (ncyc,1,imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,
     .         imp2,jmp2,kmp2,m1cc,mgrlev,igrid-1,mbloc,ibloc,
     .         mxsizc,mx1dwk,imap(1,1,ibloc,igrid-1),
     .         msegt,msegn,nseg(ibloc),p,
     .    im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .    imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .    imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .       w(m5cc(igrid-1,ibloc)),   p(m1cc(igrid-1,ibloc)),
     .     vol(m1cc(igrid-1,ibloc)),
     .     six(m1np(igrid-1,ibloc)), siy(m1np(igrid-1,ibloc)),
     .     siz(m1np(igrid-1,ibloc)), sjx(m1np(igrid-1,ibloc)),
     .     sjy(m1np(igrid-1,ibloc)), sjz(m1np(igrid-1,ibloc)),
     .     skx(m1np(igrid-1,ibloc)), sky(m1np(igrid-1,ibloc)),
     .     skz(m1np(igrid-1,ibloc)),         dti(m1fg(ibloc)), 
     .             dtj(m1fg(ibloc)),         dtk(m1fg(ibloc)), 
     .             dtl,ratioij, ratiojk,
     .      wk1d5(1,1),wk1d5(1,2),wk1d5(1,3),wk1d5(1,4),wk1d5(1,5),
     .      wk1d(1,1),wk1d(1,2),wk1d(1,3),
     .      wk2d(1,1), wk2d(1,2), wk2d(1,3), wk2d(1,4), delw,
     .      fv(m4fg(ibloc)),  fw(m5fg(ibloc)) ,
     .      dtvi(m1fg(ibloc)), dtvj(m1fg(ibloc)), dtvk(m1fg(ibloc)),
     .      fbetav2(m1fg(ibloc)) )
c
c
        call navier (imn,jmn,kmn,
     .    im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .    imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .    imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .    imap(1,1,ibloc,igrid-1), msegt, msegn, mbloc, nseg(ibloc),
     .      w(m5cc(igrid-1,ibloc)),   p(m1cc(igrid-1,ibloc)),
     .    vol(m1cc(igrid-1,ibloc)),   x(m3np(igrid-1,ibloc)),
     .    six(m1np(igrid-1,ibloc)), siy(m1np(igrid-1,ibloc)),
     .    siz(m1np(igrid-1,ibloc)), sjx(m1np(igrid-1,ibloc)),
     .    sjy(m1np(igrid-1,ibloc)), sjz(m1np(igrid-1,ibloc)),
     .    skx(m1np(igrid-1,ibloc)), sky(m1np(igrid-1,ibloc)),
     .    skz(m1np(igrid-1,ibloc)),eomu(m1cc(igrid-1,ibloc)),
     .    tauwfi1(m1is(igrid-1,ibloc)), tauwfi2(m1is(igrid-1,ibloc)),
     .    tauwfj1(m1js(igrid-1,ibloc)), tauwfj2(m1js(igrid-1,ibloc)),
     .    tauwfk1(m1ks(igrid-1,ibloc)), tauwfk2(m1ks(igrid-1,ibloc)),
     .    wk1d(1,1),wk1d(1,2),wk1d(1,3),wk1d(1,4),wk1d(1,5),
     .    wk2d(1,1),wk2d(1,2),wk2d(1,3),wk2d(1,4),wk2d(1,5),
     .    wk2d(1, 6),wk2d(1, 7),wk2d(1, 8),wk2d(1, 9),
     .    wk2d(1,10),wk2d(1,11),wk2d(1,12),wk2d(1,13), mx1dwk,
     .    delw(1,1),delw(1,2),delw(1,3),delw(1,4),delw(1,5),
     .    fv(m4fg(ibloc)),ivisb(1,ibloc) )
c
        else
c
         call filtc (1,imn,jmn,kmn,
     .    im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .    imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .    imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .      w(m5cc(igrid-1,ibloc)),   p(m1cc(igrid-1,ibloc)),
     .    vol(m1cc(igrid-1,ibloc)),         dti(m1fg(ibloc)),
     .            dtj(m1fg(ibloc)),         dtk(m1fg(ibloc)),
     .    dtl, ratioij, ratiojk,
     .    mx1dwk,wk1d(1,1),wk1d(1,2),wk1d(1,3),
     .    wk1d(1,4),wk1d(1,5),wk1d(1,6),wk1d(1,7),
     .    fw(m5fg(ibloc)) ,
     .    dtvi(m1fg(ibloc)), dtvj(m1fg(ibloc)), dtvk(m1fg(ibloc)),
     .    fbetav2(m1fg(ibloc)) )
c   
        if (inav.eq.1) call navier (imn,jmn,kmn,
     .    im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .    imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .    imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .    imap(1,1,ibloc,igrid-1), msegt, msegn, mbloc, nseg(ibloc),
     .      w(m5cc(igrid-1,ibloc)),   p(m1cc(igrid-1,ibloc)),
     .    vol(m1cc(igrid-1,ibloc)),   x(m3np(igrid-1,ibloc)),
     .    six(m1np(igrid-1,ibloc)), siy(m1np(igrid-1,ibloc)),
     .    siz(m1np(igrid-1,ibloc)), sjx(m1np(igrid-1,ibloc)),
     .    sjy(m1np(igrid-1,ibloc)), sjz(m1np(igrid-1,ibloc)),
     .    skx(m1np(igrid-1,ibloc)), sky(m1np(igrid-1,ibloc)),
     .    skz(m1np(igrid-1,ibloc)),eomu(m1cc(igrid-1,ibloc)),
     .    tauwfi1(m1is(igrid-1,ibloc)), tauwfi2(m1is(igrid-1,ibloc)),
     .    tauwfj1(m1js(igrid-1,ibloc)), tauwfj2(m1js(igrid-1,ibloc)),
     .    tauwfk1(m1ks(igrid-1,ibloc)), tauwfk2(m1ks(igrid-1,ibloc)),
     .    wk1d(1,1),wk1d(1,2),wk1d(1,3),wk1d(1,4),wk1d(1,5),
     .    wk2d(1,1),wk2d(1,2),wk2d(1,3),wk2d(1,4),wk2d(1,5),
     .    wk2d(1, 6),wk2d(1, 7),wk2d(1, 8),wk2d(1, 9),
     .    wk2d(1,10),wk2d(1,11),wk2d(1,12),wk2d(1,13), mx1dwk,
     .    delw(1,1),delw(1,2),delw(1,3),delw(1,4),delw(1,5),
     .    fv(m4fg(ibloc)),ivisb(1,ibloc) )
c
        endif
c
c     evaluate convective fluxes without dt scaling on
c     a grid 1-level finer than 'igrid'
c
        call deltaw (ncyc,1,igrid-1,imn,jmn,kmn,
     .    im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .    imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .    imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .       mx1dwk,
     .       w(m5cc(igrid-1,ibloc)),   p(m1cc(igrid-1,ibloc)),
     .     six(m1np(igrid-1,ibloc)), siy(m1np(igrid-1,ibloc)),
     .     siz(m1np(igrid-1,ibloc)), sjx(m1np(igrid-1,ibloc)),
     .     sjy(m1np(igrid-1,ibloc)), sjz(m1np(igrid-1,ibloc)),
     .     skx(m1np(igrid-1,ibloc)), sky(m1np(igrid-1,ibloc)),
     .     skz(m1np(igrid-1,ibloc)), fqs,
     .     wk1d(1,1),wk1d(1,2),wk1d(1,3),wk1d(1,4),
     .     wk1d(1,5),wk1d(1,6),wk1d(1,7),wk1d(1,8),
     .     wk1d(1,9),wk1d(1,10),wk1d(1,11),wk1d(1,12),
     .     wk1d(1,13),wk1d(1,14),wk1d(1,15),
     .     wk2d(1,1),wk2d(1,2),wk2d(1,3),wk2d(1,4),wk2d(1,5), delw, 1)
c
cprec
c        change the sequence to add viscous fluxes first for preconditioning
c       (mods by vatsa : Feb. 1998)
c
c   add viscous flux (fv) to the convective flux (delw)
c
        if ( (igrid - isoln(ilev)) .eq. 1 )
     .    call adfluxv (imn,jmn,kmn,
     .    im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .    imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .    imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .   delw, fv(m4fg(ibloc)) )
c
c
c   add dissipative flux (fw) to convective flux (delw)
c
      if (iprec.ne.3)
     .    call adfluxd (1,imn,jmn,kmn,
     .    im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .    imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .    imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .       delw,                  fw(m5fg(ibloc)) )
c
c     the time-derivative terms are moved here because of rescaling
c     of fluxes in R-K physical time-stepping scheme
c
c     add physical time derivative term depending on whether you
c     are on the coarse or fine grid
c
      if (ntorder.gt.0) then
      if ( (igrid - isoln(ilev)) .eq. 1 ) then
         call addeltf (imn,jmn,kmn,
     .    im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .    imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .    imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .    ntorder,delw,
     .    vol(m1cc(igrid-1,ibloc)), w(m5cc(igrid-1,ibloc)),
     .    wold(m5fgt(1,ibloc)),dtphy,mebdf,ibdf2opt,icyc,
     .    rkphysa,mrkstg,nrkstg,irkstg,delwork(m5fgtrk(1,ibloc)),
     .    wbar(m5fgtrk(1,ibloc)),itime )
c#ifdef  RKPHY
c# else
c    .    wold(m5fgt(1,ibloc)),dtphy )
c# endif
      else
         call addeltc (imn,jmn,kmn,
     .    im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .    imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .    imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .    ntorder,delw,
     .    vol(m1cc(igrid-1,ibloc)), w(m5cc(igrid-1,ibloc)),dtphy,
     .    rkphysa,mrkstg,nrkstg,irkstg )
c#ifdef  RKPHY
c# else
c    .    vol(m1cc(igrid-1,ibloc)), w(m5cc(igrid-1,ibloc)),dtphy )
c# endif
      end if ! igrid
      end if ! ntorder
c
c       add dissipative flux (fw) to convective flux (delw)
c       for non-conservative preconditioning
c
      if (iprec.eq.3)
     .    call adfluxd (1,imn,jmn,kmn,
     .    im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .    imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .    imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .       delw,                  fw(m5fg(ibloc)) )
c

c   add the forcing function to the residual on meshes
c   coarser than the solution mesh. i.e. evaluate
c
c                 h   h  h
c               (f - L  u )
c
        if ( (igrid-1) .ne. isoln(ilev) )
     .  call adforfn (imn,jmn,kmn,
     .   im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .   imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .   imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .      delw,          wr(m5cg(igrid-1,ibloc)))
c
c
c                 2h  h        2h  h    h h
c   calculate    I  (u )  and I  (f  - L u )
c                 h            h
c
c   save the fine-grid values restricted to coarse-grids
c   so that they can be subtracted from the iterated
c   coarse-grid values during prolongotation (FAS)
c
c
       call restrct (iftf,igrid,imn,jmn,kmn,
     .  im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .  imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .  imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .  im  (igrid  ,ibloc), jm  (igrid  ,ibloc), km  (igrid  ,ibloc),
     .  imp1(igrid  ,ibloc), jmp1(igrid  ,ibloc), kmp1(igrid  ,ibloc),
     .  imp2(igrid  ,ibloc), jmp2(igrid  ,ibloc), kmp2(igrid  ,ibloc),
     .    w(m5cc(igrid-1,ibloc)),  w(m5cc(igrid,ibloc)),
     .    p(m1cc(igrid-1,ibloc)),  p(m1cc(igrid,ibloc)),
     .         delw,              wr(m5cg(igrid,ibloc)),
     .  vol(m1cc(igrid-1,ibloc)),jbloc                 )
c
       call restrev (iftf,igrid,imn,jmn,kmn,
     .  im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .  imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .  imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .  im  (igrid  ,ibloc), jm  (igrid  ,ibloc), km  (igrid  ,ibloc),
     .  imp1(igrid  ,ibloc), jmp1(igrid  ,ibloc), kmp1(igrid  ,ibloc),
     .  imp2(igrid  ,ibloc), jmp2(igrid  ,ibloc), kmp2(igrid  ,ibloc),
     .    eomu(m1cc(igrid-1,ibloc)),  eomu(m1cc(igrid,ibloc)),
     .  vol(m1cc(igrid-1,ibloc)),jbloc                 )
c
c
c     convert the residuals to "p,u,v,w,T" variables (iprec.ge.2)
c                     and precondition
c
      if (iprec.ge.2)
     .   call transct (imn,jmn,kmn,
     .   im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
     .   imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
     .   imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
     .   w(m5cc(igrid-1,ibloc)),   p(m1cc(igrid-1,ibloc)), delw )
ccc   if (ijac.eq.1 .and. iprec.gt.0) iprecg = - iprec
ccc   if (iprec.ge.1)
ccc  .   call precong (imn,jmn,kmn,
ccc  .   im  (igrid-1,ibloc), jm  (igrid-1,ibloc), km  (igrid-1,ibloc),
ccc  .   imp1(igrid-1,ibloc), jmp1(igrid-1,ibloc), kmp1(igrid-1,ibloc),
ccc  .   imp2(igrid-1,ibloc), jmp2(igrid-1,ibloc), kmp2(igrid-1,ibloc),
ccc  .   six(m1np(igrid-1,ibloc)), siy(m1np(igrid-1,ibloc)),
ccc  .   siz(m1np(igrid-1,ibloc)), sjx(m1np(igrid-1,ibloc)),
ccc  .   sjy(m1np(igrid-1,ibloc)), sjz(m1np(igrid-1,ibloc)),
ccc  .   skx(m1np(igrid-1,ibloc)), sky(m1np(igrid-1,ibloc)),
ccc  .   skz(m1np(igrid-1,ibloc)),
ccc  .   w(m5cc(igrid-1,ibloc)),   p(m1cc(igrid-1,ibloc)), delw ,
ccc  .      fbeta2(m1fg(ibloc)) )
c
c
c---- end test for local block
c
      endif
c
  210 continue
c
c----------------------------------------------------------------
c   apply the boundary conditions on conserved 
c   variables 'w' (and eomu, turv1) at this grid level
c
      call bcflow (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  m1np,m2np,m3np,m4np,m5np,m1cc,m2cc,m3cc,m4cc,m5cc,m5cg,
     .  m1is,m1js,m1ks,m3is,m3js,m3ks,mgrlev,igrid,igrid,mbloc,nbloc,
     .  mxsizc,mxsiz5c,mxsizn,mxsiz3n,mx1dwk,
     .  imap(1,1,1,igrid),msegt,msegn,nseg,w,p,eomu,turv1,turv2,vol,x,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .  ri1,ri2,rj1,rj2,rk1,rk2,
     .  tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .  uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
     .  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .  wk2d,const1,const2,const3,const4,const5,
     .  mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .  mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .  lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .  kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .  ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,iunsteady,totptim,smin)
c
c-----------------------------------------------------------------
c   save the initial guess for the new restricted w values (w1)
c
      do 215 ibloc = 1,nbloc
c
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank)
# else
      if (nodes(ibloc).eq.myrank)
# endif
     .call equw1w (imn,jmn,kmn,
     .            imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .            imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .            w(m5cc(igrid,ibloc)), w1(m5cg(igrid,ibloc)),5 )
c
c
  215 continue
c
c----------------------------------------------------------------
c   end the logic for evaluating forcing function
c   on meshes coarser than the solution mesh
c
  200 continue
c
c------------ begin the subiteration in restriction loop ---------
c
      mitr = gitr(igrid)
      if(igrid .eq. isoln(ilev)) mitr = 1
c
      do 510 nitr = 1,mitr
c
c----------------------------------------------------------------
c
c     use the Runge-Kutta solver for advancing the solution in time
c
      iforfn = 1
      ilast  = 0
      if(igrid.eq.ibegr(irleg,ilev)) iforfn = 0
c
      call solve (imn, jmn, kmn, im,  jm,  km,
     .            imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .            m1np,m2np,m3np,m4np,m5np, m1cc,m2cc,m3cc,m4cc,m5cc,
     .            m1fg,m4fg,m5fg,m5cg,m1is,m1js,m1ks,m3is,m3js,m3ks,
     .            igrid,isoln(ilev),nbloc,mres,
     .            mx1dwk,mx3dwk,mxsizc,mxsiz3c,mxsiz4c,mxsiz5c,
     .            mxszfg,mxsz4fg,mxsz5fg,mxsz5cg,mxsizn,mxsiz3n,
     .            imap,nseg,msegt,msegn,mbloc,mgrlev,
     .            w,p,ws,wr,delw,fw,fv,eomu,turv1,turv2,vol,x,fqs,
     .            six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .            ri1,ri2,rj1,rj2,rk1,rk2,
     .            tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .            uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
     .            uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .            tnuwfi1,tnuwfi2,tnuwfj1,tnuwfj2,tnuwfk1,tnuwfk2,
     .            dtl,dti,dtj,dtk,dtmin,ratioij,ratiojk,eprs,
     .            wk1d,wk2d,wk1d5,rtrmsg,
     .            nitr,iforfn,ilast,ncyc,ibcup,ilev,nres,
     .            ivisb,const1,const2,const3,const4,const5,mstrt,
     .            smin,iturbb,itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,
     .            mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .            mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .            lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .            kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .            ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,
     .            dtvi,dtvj,dtvk,fbeta2,fbetav2,wold,tv1old,tv2old,
     .            m1fgt,m5fgt,iturv1,iturv2,mtorder,ntorder,iunsteady,
     .            dtphy,totptim,m1fgtrk,m5fgtrk,irkstg,nrkstg,mrkstg,
     .            rkphysa,rkphyse,delwork,dtv1ork,dtv2ork,
     .            mebdf,ibdf2opt,icyc,itime,wbar,tv1bar,tv2bar)
c#ifdef  RKPHY
c# else
c     .            dtphy,totptim )
c# endif
c
      fgrd = 1./(2.**real((igrid-1)))
      work  = work +fgrd*fgrd*fgrd
c
      if (nitr.eq.1 .and. igrid.eq.isoln(ilev)) then
c
         wcount(nres) = totw + work
c
c        convergence parameters over all blocks
c
         jbloc = 0
         do 270 ibloc = 1,nbloc
# if defined BUILD_MPI
         if (nodes(ibloc)-1.eq.myrank) then
# else
         if (nodes(ibloc).eq.myrank) then
# endif
         jbloc = jbloc + 1
         call ckhconv (imn,jmn,kmn,
     .     im  (igrid,ibloc),   jm  (igrid,ibloc),km  (igrid,ibloc),
     .     imp2(igrid,ibloc),   jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .   w(m5cc(igrid,ibloc)),p(m1cc(igrid,ibloc)),
     .   nsupg(nres),hrmsg(nres),mx1dwk,wk1d(1,1),jbloc,nodes(ibloc))
         ammax     = max(ammax,amachg)
         endif
  270    continue
c
c        compute force/moment coefficients 
c 
         call force (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .   m1cc,m5cc,m1np,m3np,m1is,m1js,m1ks,m3is,m3js,m3ks,
     .   mgrlev,igrid,mbloc,nbloc,
     .   mxsizc,mxsiz5c,mxsizn,mxsiz3n,mx1dwk,
     .   imap(1,1,1,igrid),msegt,msegn,nseg,w,p,vol,x,
     .   six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .   ri1,ri2,rj1,rj2,rk1,rk2,
     .   tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .   wk2d(1,1),wk2d(1,2),wk2d(1,3),wk2d(1,4),
     .   scal,xref,yref,zref,sref,
     .   cxg(nres),cyg(nres),czg(nres),
     .   cxvg(nres),cyvg(nres),czvg(nres),
     .   cmxg(nres),cmyg(nres),cmzg(nres),
     .   cx,cy,cz,cxv,cyv,czv,cmx,cmy,cmz, ivisb, iforce,nodes )
c
c       changed above argument list
c       dropped clg(nres),cdtg(nres),cdvg(nres) and in that place,
c       added cxvg(nres),cyvg(nres),czvg(nres),cxg(nres),cyg(nres)
c       also dropped cl, cdt and cdv and added there
c       cxv,cyv,czv,cx,cy
c

      endif
c
c------ end sub-iteration loop -----------------------------
  510 continue
c
c------ end restriction loop   -----------------------------
  500 continue
c
c-----------------------------------------------------------
c       begin prolongation loop  starting with the
c       coarsest mesh and ending at the finest mesh
c-----------------------------------------------------------
c
      if(mgtype.eq.0) go to 601
      do 600 igrid = icoars(irleg,ilev)-1, ifine(irleg,ilev), -1
c
c
c                 h      h      h   2h    2h  h
c     calculate  w    = w    + I  (w   - I  (w )   )
c                 new    old    2h        h     old
c
c
      jbloc = 0
      do 310 ibloc=1,nbloc
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) then
# else
      if (nodes(ibloc).eq.myrank) then
# endif
        jbloc = jbloc +1
c
        call prolong (iftf,igrid,imn,jmn,kmn,
     .  im  (igrid  ,ibloc), jm  (igrid  ,ibloc), km  (igrid  ,ibloc),
     .  imp1(igrid  ,ibloc), jmp1(igrid  ,ibloc), kmp1(igrid  ,ibloc),
     .  imp2(igrid  ,ibloc), jmp2(igrid  ,ibloc), kmp2(igrid  ,ibloc),
     .  im  (igrid+1,ibloc), jm  (igrid+1,ibloc), km  (igrid+1,ibloc),
     .  imp1(igrid+1,ibloc), jmp1(igrid+1,ibloc), kmp1(igrid+1,ibloc),
     .  imp2(igrid+1,ibloc), jmp2(igrid+1,ibloc), kmp2(igrid+1,ibloc),
     .    w(m5cc(igrid,ibloc)), w(m5cc(igrid+1,ibloc)),
     .       delw, wk2d(1,1),wk2d(1,2),
     .    p(m1cc(igrid,ibloc)),w1(m5cg(igrid+1,ibloc)),jbloc )
c
      endif
  310 continue
c
c----------------------------------------------------------------
c   apply the boundary conditions on conserved 
c   variables 'w' (and eomu, turv1) at this grid level
c
      call bcflow (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  m1np,m2np,m3np,m4np,m5np,m1cc,m2cc,m3cc,m4cc,m5cc,m5cg,
     .  m1is,m1js,m1ks,m3is,m3js,m3ks,mgrlev,igrid,igrid,mbloc,nbloc,
     .  mxsizc,mxsiz5c,mxsizn,mxsiz3n,mx1dwk,
     .  imap(1,1,1,igrid),msegt,msegn,nseg,w,p,eomu,turv1,turv2,vol,x,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .  ri1,ri2,rj1,rj2,rk1,rk2,
     .  tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .  uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
     .  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .  wk2d,const1,const2,const3,const4,const5,
     .  mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .  mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .  lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .  kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .  ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,iunsteady,totptim,smin)
c
c----------------------------------------------------------------
c     begin subiteration loop on prolongation
c----------------------------------------------------------------
c
      mitp = gitp(igrid)
c
      if (mitp .eq. 0) go to 600
      if (igrid .eq. isoln(ilev)) mitp = 1
c
      do 610 nitr=1,mitp
c
c     perform fixed no. of iterations in the prolongation loop
c     of the multigrid cycle based on input parameters
c
c----------------------------------------------------------------
c
c     use the Runge-Kutta solver for advancing the solution in time
c
      iforfn = 0
      ilast  = 1
c
      call solve (imn, jmn, kmn, im,  jm,  km,
     .            imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .            m1np,m2np,m3np,m4np,m5np, m1cc,m2cc,m3cc,m4cc,m5cc,
     .            m1fg,m4fg,m5fg,m5cg,m1is,m1js,m1ks,m3is,m3js,m3ks,
     .            igrid,isoln(ilev),nbloc,mres,
     .            mx1dwk,mx3dwk,mxsizc,mxsiz3c,mxsiz4c,mxsiz5c,
     .            mxszfg,mxsz4fg,mxsz5fg,mxsz5cg,mxsizn,mxsiz3n,
     .            imap,nseg,msegt,msegn,mbloc,mgrlev,
     .            w,p,ws,wr,delw,fw,fv,eomu,turv1,turv2,vol,x,fqs,
     .            six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .            ri1,ri2,rj1,rj2,rk1,rk2,
     .            tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .            uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
     .            uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .            tnuwfi1,tnuwfi2,tnuwfj1,tnuwfj2,tnuwfk1,tnuwfk2,
     .            dtl,dti,dtj,dtk,dtmin,ratioij,ratiojk,eprs,
     .            wk1d,wk2d,wk1d5,rtrmsg,
     .            nitr,iforfn,ilast,ncyc,ibcup,ilev,nres,
     .            ivisb,const1,const2,const3,const4,const5,mstrt,
     .            smin,iturbb,itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,
     .            mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .            mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .            lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .            kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .            ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,
     .            dtvi,dtvj,dtvk,fbeta2,fbetav2,wold,tv1old,tv2old,
     .            m1fgt,m5fgt,iturv1,iturv2,mtorder,ntorder,iunsteady,
     .            dtphy,totptim,m1fgtrk,m5fgtrk,irkstg,nrkstg,mrkstg,
     .            rkphysa,rkphyse,delwork,dtv1ork,dtv2ork,
     .            mebdf,ibdf2opt,icyc,itime,wbar,tv1bar,tv2bar)
c#ifdef  RKPHY
c# else
c     .            dtphy,totptim )
c# endif
c
      fgrd = 1./(2.**real((igrid-1)))
      work  = work +fgrd*fgrd*fgrd
c
c
c------ end sub-iteration loop -----------------------------
  610 continue
c
c------ end prolongation loop  -----------------------------
  600 continue
  601 continue
c
c------ end the loop on legs of the mg cycle
 1020 continue
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c------ end the node process execution
cjm      endif
c
c######################################################################
c all processes
c     collect convergence results and print and save intermediate
c     solution if msave > 0
c######################################################################
c
      call getconv (mbloc,nbloc,wcount(nres),rtrmsg(nres),hrmsg(nres),
     .              nsupg(nres),cxg(nres),cyg(nres),czg(nres),
     .              cxvg(nres),cyvg(nres),czvg(nres),
     .              cmxg(nres),cmyg(nres),cmzg(nres),
     .              cx,cy,cz,cxv,cyv,czv,cmx,cmy,cmz,totw,wk2d,nodes)
c
      call MPI_Bcast (amachg,1,RTYPE,master,
     .                mycomm,ierr)
#endif
c
      if (myrank.eq.master) then

cyho    new variables (for global forces)

c       cng     normal force coefficient. represents force vector
c               component on the maneuvering plane and normal to the x-axis
c               maneuvering plane is defined as a plane which include x-(body)
c               axis and velociy (total) vector.
c       vnx     absolute value of total lateral velocity vector normal to u0
c       vnx     = sqrt(v0**2 + w0**2)
c       if(vnx.gt.1.e-3) then
c         cosp  = v0/vnx
c         sinp  = w0/vnx
c       else
c         cosp  = 1.0
c         sinp  = 0.0
c       endif
c  coordinate transformation
        cosp = cos(phi)
        sinp = sin(phi)
        cng (nres) = cyg(nres)*cosp + czg(nres)*sinp
        clg (nres) = cyg(nres)*ca   - cxg(nres)*sa
        cdtg(nres) = (cyg(nres)*sa  + cxg(nres)*ca)*cos(betas)
     .                  + czg(nres)*sin(betas)
        cdvg(nres) = (cyvg(nres)*sa  + cxvg(nres)*ca)*cos(betas)
     .                  + czvg(nres)*sin(betas)
c
        if (icyc.eq.1) rtrmsc0 = rtrmsg(nres)
        if (icyc.le.2) then
            ratemax = 1.0
            raterms = 1.0
            rtmaxgp = rtmaxg
            rtrmsgp = rtrmsg(nres)
            rtmaxg0 = rtmaxgp
            rtrmsg0 = rtrmsgp
           
        else
cvn         raterms = rtrmsg(nres)/rtrmsgp
cvn         ratemax = rtmaxg/rtmaxgp
            raterms = (rtrmsg(nres)/rtrmsg0)**(1./real(icyc  -1))
            ratemax = (rtmaxg/rtmaxg0)**(1./real(icyc  -1))
            rtmaxgp = rtmaxg
            rtrmsgp = rtrmsg(nres)
        end if
c
      if (iwrit.gt.0)
     .    write (iwrit,650) icyc,wcount(nres),
     .    rtmaxg,irtmxg,jrtmxg,krtmxg,mrtbloc,rtrmsg(nres),
     .    hmaxg, ihmxg, jhmxg, khmxg, mhbloc,  ratemax,
     .    raterms,clg(nres),cdtg(nres),cdvg(nres),cmzg(nres)
c    .    hmaxg, ihmxg, jhmxg, khmxg, mhbloc,  hrmsg(nres),
c    .    nsupg(nres),clg(nres),cdtg(nres),cdvg(nres),cmzg(nres)

          write (irhtall,650) icyc,wcount(nres),
     .    rtmaxg,irtmxg,jrtmxg,krtmxg,mrtbloc,rtrmsg(nres),
     .    hmaxg, ihmxg, jhmxg, khmxg, mhbloc,  ratemax,
     .    raterms,clg(nres),cdtg(nres),cdvg(nres),cmzg(nres)
c    .    hmaxg, ihmxg, jhmxg, khmxg, mhbloc,  hrmsg(nres),
c    .    nsupg(nres),clg(nres),cdtg(nres),cdvg(nres),cmzg(nres)
cjm
      call flush(iwrit)
      endif
c
c------ dump out flow field for future restart
c
      if (msave.gt.0) then
c
      if (mod(icyc,msave).eq.0.and.irkstg.eq.nrkstg)
     .  call savsol2(imn,jmn,kmn,imp2,jmp2,kmp2,m1cc,m5cc,mxsizc,
     .               mxsiz5c,mgrlev,mbloc,nbloc,mres,nres,
     .               w,p,eomu,turv1,turv2,smin,wcount,rtrmsg,
     .               hrmsg,clg,cmzg,cdtg,cdvg,nsupg,
     .               mxszfg,wold,tv1old,tv2old,m5fgt,m1fgt,mtorder,
     .               ntorder,iturv1,iturv2,iunsteady,totptim,
# if defined BUILD_MPI
     .               scal,xref,yref,zref,sref,
     .               buffw,buffp,buffe,buff1,buff2,buffs,
     .               mx3dwk,buffwo,bufftv1,bufftv2,nodes )
# else
     .               scal,xref,yref,zref,sref,nodes)
# endif
c
      endif
c
c     collect fluxes of gov. equations for RKPHY and MEBDF schemes
c
cvn   if(irkcase.gt.0.or.(mebdf.gt.0.and.irkstg.lt.nrkstg)) then
      if(irkcase.gt.0.or.(mebdf.gt.0.and.irkstg.lt.nrkstg).or.
     .                   (mebdf.eq.0.and.irkstg.eq.1) ) then
        igrid = 1
        cfl   = cflf

c       call fluxwrk (imn, jmn, kmn, im,  jm,  km,
c    .            imp1,jmp1,kmp1,imp2,jmp2,kmp2,
c    .            m1np,m2np,m3np,m4np,m5np, m1cc,m2cc,m3cc,m4cc,m5cc,
c    .            m1fg,m4fg,m5fg,m5cg, m1is,m1js,m1ks,m3is,m3js,m3ks,
c    .            igrid,isoln(ilev),nbloc,mres,
c    .            mx1dwk,mx3dwk,mxsizc,mxsiz3c,mxsiz4c,mxsiz5c,
c    .            mxszfg,mxsz4fg,mxsz5fg,mxsz5cg,mxsizn,mxsiz3n,
c    .            imap,nseg,msegt,msegn,mbloc,mgrlev,
c    .            w,p,ws,wr,delw,fw,fv,eomu,turv1,turv2,vol,x,fqs,
c    .            six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
c    .            ri1,ri2,rj1,rj2,rk1,rk2,
c    .            tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
c    .            uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
c    .            uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
c    .            tnuwfi1,tnuwfi2,tnuwfj1,tnuwfj2,tnuwfk1,tnuwfk2,
c    .            dtl,dti,dtj,dtk,dtmin,ratioij,ratiojk,eprs,
c    .            wk1d,wk2d,wk1d5,rtrmsg,
c    .            nitr,iforfn,ilast,ncyc,ibcup,ilev,nres,
c    .            ivisb,const1,const2,const3,const4,const5,mstrt,
c    .            smin,iturbb,itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,
c    .            mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
c    .            mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
c    .            lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
c    .            kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
c    .            ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,
c    .            dtvi,dtvj,dtvk,fbeta2,fbetav2,wold,tv1old,tv2old,
c    .            m1fgt,m5fgt,iturv1,iturv2,mtorder,ntorder,iunsteady,
c    .            dtphy,totptim,m1fgtrk,m5fgtrk,irkstg,nrkstg,mrkstg,
c    .            rkphysa,rkphyse,delwork,dtv1ork,dtv2ork )
c
      endif
c
#ifdef RKPHY
      if (rtmaxg.lt.ftolrk.and.icyc.gt.2) then
         go to 2003
      else if (icyc.eq.mcyc(ilev) .and. ntorder.ne.0) then
         if(myrank.eq.master)
     .   write(*, '(2x,"no convergence in multigrid cycle of this
     .       physical time step")')
cvn      stop
      endif
# endif
c     end of logic for RKPHY exit criterion
c
c------ end the loop on mg cycles (for a fixed fmg level)
c     call flush (iwrit)
c
        igrid = 1
        cfl   = cflf
        jbloc   =  0

        do ibloc = 1,nbloc
# if defined BUILD_MPI
          if (nodes(ibloc)-1.eq.myrank) then
# else
          if (nodes(ibloc).eq.myrank) then
# endif
          jbloc = jbloc + 1
c

          if(mebdf.gt.0.and.irkstg.eq.3)
     .    call errmebdf(imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .              delw,dtl,vol(m1cc(igrid,ibloc)),dtphy,
     .              mx1dwk,wk1d(1,1),jbloc,nodes(ibloc),
     .              w(m5cc(igrid,ibloc)),
     .              wbar(m5fgtrk(1,ibloc)),ntorder )

          if(mebdf.eq.0.and.ibdf2opt.eq.0.and.ntorder.eq.3)
     .    call errbdf3(imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .              delw,dtl,vol(m1cc(igrid,ibloc)),dtphy,
     .              mx1dwk,wk1d(1,1),jbloc,nodes(ibloc),
     .              w(m5cc(igrid,ibloc)),
     .              wold(m5fgt(1,ibloc)),
     .              delwork(m5fgtrk(1,ibloc)),ntorder )

          if(ibdf2opt.eq.1.and.ntorder.eq.3)
     .    call errbdf2opt(imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .              delw,dtl,vol(m1cc(igrid,ibloc)),dtphy,
     .              mx1dwk,wk1d(1,1),jbloc,nodes(ibloc),
     .              w(m5cc(igrid,ibloc)),
     .              wold(m5fgt(1,ibloc)),
     .              delwork(m5fgtrk(1,ibloc)),ntorder )

          if(ibdf2opt.eq.0.and.ntorder.eq.2)
     .    call errbdf2(imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .              delw,dtl,vol(m1cc(igrid,ibloc)),dtphy,
     .              mx1dwk,wk1d(1,1),jbloc,nodes(ibloc),
     .              w(m5cc(igrid,ibloc)),
     .              wold(m5fgt(1,ibloc)),
     .              delwork(m5fgtrk(1,ibloc)),ntorder )

          endif
c         end the local block test
        enddo
c     end of error estimation logic for MEBDF and BDF3
c
c     start logic for MEBDF exit criterion
c
      if(myrank.eq.master) then

      if(mebdf.gt.0.and.itime.gt.1) then
        if(icyc.le.3) go to 1010
         errest = ftol*erkrmsg(1)
         if(rtrmsg(nres).lt.errest) then
           write(*, '(2x,"ftol,erkrmsg,rtrmsg", 3e12.5)')
     .                    ftol,erkrmsg(1),rtrmsg(nres)
            go to 2003
         else if (icyc.eq.mcyc(ilev).and.myrank.eq.master) then
           write(*, '(2x,"no convergence in multigrid cycle of this
     .     physical time step for MEBDF")')
           write(*, '(2x,"ftol,erkrmsg,rtrmsg", 3e12.5)')
     .                    ftol,erkrmsg(1),rtrmsg(nres)
         endif
      endif
c
c     end logic for MEBDF exit criterion
c
c     start logic for BDF3 exit criterion
c
      if(mebdf.eq.0.and.(ntorder.eq.3.or.ntorder.eq.2).
     .   and.itime.gt.1) then
        if(icyc.le.3) go to 1010
         errest = ftol*erkrmsg(1)
         if(rtrmsg(nres).lt.errest) then
           write(*, '(2x,"ftol,erkrmsg,rtrmsg", 3e12.5)')
     .                    ftol,erkrmsg(1),rtrmsg(nres)
           go to 2003
         else if (icyc.eq.mcyc(ilev).and.myrank.eq.master) then
           write(*, '(2x,"no convergence in multigrid cycle of this
     .     physical time step for BDF3/BDF2OPT")')
           write(*, '(2x,"ftol,erkrmsg,rtrmsg", 3e12.5)')
     .                    ftol,erkrmsg(1),rtrmsg(nres)
         endif
      endif
      endif
c
c     end logic for BDF3/BDF2OPT exit criterion
c
      resrat    = rtrmsg(nres)/rtrmsc0
      if (resrat .le. resend) go to 2003
 1010 continue
c
 2003 continue
c
      if(mebdf.gt.0) then
c
c     Store the latest stage values for better estimates in MEBDF scheme
c
        do ibloc = 1,nbloc
c
# if defined BUILD_MPI
          if (nodes(ibloc)-1.eq.myrank) then
# else
          if (nodes(ibloc).eq.myrank) then
# endif
            igrid  = 1

            if(irkstg.le.2)
     .      call equw1w (imn,jmn,kmn,
     .            imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .            imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .            w(m5cc(igrid,ibloc)), wbar(m5fgtrk(irkstg,ibloc)),5 )
c
            if(irkstg.le.2.and.iturb.ge.2)
     .      call equw1w (imn,jmn,kmn,
     .            imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .            imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .      turv1(m1cc(igrid,ibloc)), tv1bar(m1fgtrk(irkstg,ibloc)),1 )
c
            if(irkstg.le.2.and.iturb.eq.3)
     .      call equw1w (imn,jmn,kmn,
     .            imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .            imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .      turv2(m1cc(igrid,ibloc)), tv2bar(m1fgtrk(irkstg,ibloc)),1 )
c
c
c           Use a better guess on initial solution for MEBDF step III
c
            if(irkstg.eq.2)
     .      call equw1w (imn,jmn,kmn,
     .            imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .            imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .            wbar(m5fgtrk(irkstg-1,ibloc)),w(m5cc(igrid,ibloc)),5 )
c         update pressure
            if(irkstg.eq.2)
     .      call press (imn,jmn,kmn,
     .        im  (igrid,ibloc),jm  (igrid,ibloc),km  (igrid,ibloc),
     .        imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .        imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .        w(m5cc(igrid,ibloc)),               p(m1cc(igrid,ibloc)) )
c
            if(irkstg.eq.2.and.iturb.ge.2)
     .      call equw1w (imn,jmn,kmn,
     .            imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .            imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .      tv1bar(m1fgtrk(irkstg-1,ibloc)),turv1(m1cc(igrid,ibloc)),1 )
c
            if(irkstg.eq.2.and.iturb.eq.3)
     .      call equw1w (imn,jmn,kmn,
     .            imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
     .            imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .      tv2bar(m1fgtrk(irkstg-1,ibloc)),turv2(m1cc(igrid,ibloc)),1 )
c
          endif
        enddo
c
c       update the b.c.'s on new estimates of "w"
c
        if(irkstg.eq.2)
     .  call bcflow (imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  m1np,m2np,m3np,m4np,m5np,m1cc,m2cc,m3cc,m4cc,m5cc,m5cg,
ccvn .  m5cg,m3is,m3js,m3ks,mgrlev,igrid,isoln(1),mbloc,nbloc,
     .  m1is,m1js,m1ks,m3is,m3js,m3ks,mgrlev,igrid,igrid,mbloc,nbloc,
     .  mxsizc,mxsiz5c,mxsizn,mxsiz3n,mx1dwk,
     .  imap(1,1,1,igrid),msegt,msegn,nseg,w,p,eomu,turv1,turv2,vol,x,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .  ri1,ri2,rj1,rj2,rk1,rk2,
     .  tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .  uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
     .  uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .  wk2d,const1,const2,const3,const4,const5,
     .  mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .  mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .  lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .  kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .  ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,iunsteady,totptim,smin)

c       end of b.c. setup for new estimates on "w" for MEBDF4
c
      endif
c     end physical RK-stage (or) MEBDF stage loop
 2004 continue
#ifdef RKPHY
c     use the info. stored in delwold array to estimate
c     the starting conditions for next step
c
      do ibloc=1,nbloc
# if defined BUILD_MPI
      if (nodes(ibloc)-1.eq.myrank) then
# else
      if (nodes(ibloc).eq.myrank) then
# endif
      igrid = 1
cvn   call predictrk(imn,jmn,kmn,
cvn  .            im  (igrid,ibloc),jm  (igrid,ibloc),km  (igrid,ibloc),
cvn  .            imp1(igrid,ibloc),jmp1(igrid,ibloc),kmp1(igrid,ibloc),
cvn  .            imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
cvn  .            irkcase,nrkstg,dtphy,dtphy,rkphysa,
cvn  .            vol(m1cc(igrid,ibloc)),
cvn  .            w(m5cc(igrid,ibloc)),delwork(m5fgtrk(igrid,ibloc)) )
      endif
      enddo
#endif
c
c     average flow quantities for unsteady flow
c
cvnv  if(iunsteady.gt.0.and.mstrt.gt.0.and.mavg.gt.0) then
      if(iunsteady.gt.0.and.mavg.gt.0) then
        do ibloc=1,nbloc
# if defined BUILD_MPI
        if (nodes(ibloc)-1.eq.myrank) then
# else
        if (nodes(ibloc).eq.myrank) then
# endif
        igrid = 1
        call flowavg(imn,jmn,kmn,
     .            imp2(igrid,ibloc),jmp2(igrid,ibloc),kmp2(igrid,ibloc),
     .            itime,p(m1cc(igrid,ibloc)),pavg(m1fg(ibloc)),
     .            w(m5cc(igrid,ibloc)),wavg(m5fg(ibloc)) ,
     .            eomu(m1cc(igrid,ibloc)),eomuavg(m1fg(ibloc)) )
        endif
        enddo
      endif
c
c     write out the current solution to the plotting files
c
        if(mplot .gt. 0. and. iunsteady.gt.0) then
           if(mod(itime,mplot).eq.0)
     .     call savtim(imn,jmn,kmn,
     .                 im,jm,km,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .                 mxsizc,mxsiz5c,mxsiz3n,mgrlev,mbloc,msegt,msegn,
     .                 m5cc,m1cc,m3np,imap,nbloc,nseg,
# if defined BUILD_MPI
     .                 mxszfg,buffw,buffp,buff,
# endif
     .                 x,w,p,delw,alpha,totptim,0,nodes)
           call savdata(ibeggrp,iendgrp,iskpgrp,
     .                  jbeggrp,jendgrp,jskpgrp,
     .                  kbeggrp,kendgrp,kskpgrp,
     .                  iblkgrp,ngroup,
     .                  imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2,
     .                  mxsizc,mxsiz5c,mxsiz3n,mgrlev,mbloc,msegt,msegn,
     .                  m5cc,m1cc,m3np,imap,nbloc,nseg,
# if defined BUILD_MPI
     .                  mxszfg,buffw,buffp,buff,
# endif
     .                  x,w,p,delw,alpha,totptim,0,nodes)
        end if
c
c------ dump out the average flow field for unsteady flows
c
      if (mavg.gt.0. and.iunsteady.gt.0.and.itime.gt.1) then
c
      if (mod(itime,mavg).eq.0.and.irkstg.eq.nrkstg)
     .  call savsolavg(imn,jmn,kmn,imp2,jmp2,kmp2,m1fg,m5fg,
     .               mxszfg,mxsz5fg,mgrlev,mbloc,nbloc,
     .               wavg,pavg,eomuavg,wcount(nres),rtrmsg(nres),
     .               hrmsg(nres),clavg,cmzavg,cdtavg,cdvavg,nsupg(nres),
     .               iunsteady,totptim,
# if defined BUILD_MPI
     .               scal,xref,yref,zref,sref,
     .               buffw,buffp,buffe,nodes )
# else
     .               scal,xref,yref,zref,sref,nodes)
# endif
c
      endif
c
c
        if(myrank.eq.master) then
          ttrum = totptim * sqrt(gamma) * rm
          cdp   = cdtg(nres) - cdvg(nres)
          if(iunsteady.gt.0.and.mavg.gt.0) then
            call forcavg(itime,clg(nres),cdtg(nres),cdvg(nres),
     .                   cmzg(nres),clavg,cdtavg,cdvavg,cmzavg )
            write (irhtphy ,'(i5,2e11.4,8e13.5)')
     .           itime,totptim,ttrum,clg(nres),cdtg(nres),
     .           cdvg(nres),cmzg(nres),clavg,cdtavg,
     .           erkmaxg(1),erkrmsg(1)
          else
            write (irhtphy ,'(i5,2e11.4,6e17.9)')
     .             itime,totptim,ttrum,clg(nres),cdtg(nres),
     .             cdvg(nres),cmzg(nres),
     .             erkmaxg(1),erkrmsg(1)
c    .             erkmaxg(2),erkrmsg(2)
c           write (irhtphy ,'(i5,2e11.4,8e13.5)')
c    .             itime,totptim,ttrum,clg(nres),cdtg(nres),
c    .             cdvg(nres),cmzg(nres),
c    .             erkmaxg(1),erkrmsg(1),
c    .             erkmaxg(2),erkrmsg(2)
          endif
        endif
c#ifdef RKPHY
c     if(irkcase.gt.0.or.mebdf.gt.0) then
c       if(myrank.eq.master) then
c         write (ipltfil ,'(i6,10e12.5)')
c    .         itime,totptim,ttrum,
c    .         erkmaxg(1),erkmaxg(2),erkmaxg(3),erkmaxg(5),
c    .         erkrmsg(1),erkrmsg(2),erkrmsg(3),erkrmsg(5) 
c       endif
c     endif
c#endif
c
c------ end the loop on physical time loop
 1005  continue
c
c######################################################################
c node processes
c     do one final iteration on the finest grid level to 
c        give some additional smoothing to the solution.
c######################################################################
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
cjm      if (myrank.ne.master) then
#endif
c
      igrid = isoln(ilev)
      cfl   = cflf
      nitr  = 1
c
c----------------------------------------------------------------
c     use the Runge-Kutta solver for advancing the solution in time
c
      iforfn = 1
      ilast  = 1
c
#ifdef  RKPHY
# else
cvn   if(irkcase.eq.0)
c fix to avoid call to solve here in time-accurate mode
cvn   if(irkcase.eq.0.and.iunsteady.eq.0 )
      if(irkcase.eq.0)
     .call solve (imn, jmn, kmn, im,  jm,  km,
     .            imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .            m1np,m2np,m3np,m4np,m5np, m1cc,m2cc,m3cc,m4cc,m5cc,
     .            m1fg,m4fg,m5fg,m5cg,m1is,m1js,m1ks,m3is,m3js,m3ks,
     .            igrid,isoln(ilev),nbloc,mres,
     .            mx1dwk,mx3dwk,mxsizc,mxsiz3c,mxsiz4c,mxsiz5c,
     .            mxszfg,mxsz4fg,mxsz5fg,mxsz5cg,mxsizn,mxsiz3n,
     .            imap,nseg,msegt,msegn,mbloc,mgrlev,
     .            w,p,ws,wr,delw,fw,fv,eomu,turv1,turv2,vol,x,fqs,
     .            six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .            ri1,ri2,rj1,rj2,rk1,rk2,
     .            tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .            uswfi1,uswfi2,uswfj1,uswfj2,uswfk1,uswfk2,
     .            uplwfi1,uplwfi2,uplwfj1,uplwfj2,uplwfk1,uplwfk2,
     .            tnuwfi1,tnuwfi2,tnuwfj1,tnuwfj2,tnuwfk1,tnuwfk2,
     .            dtl,dti,dtj,dtk,dtmin,ratioij,ratiojk,eprs,
     .            wk1d,wk2d,wk1d5,rtrmsg,
     .            nitr,iforfn,ilast,ncyc,ibcup,ilev,nres,
     .            ivisb,const1,const2,const3,const4,const5,mstrt,
     .            smin,iturbb,itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,
     .            mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .            mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .            lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .            kpitmb1,kpitmb2,frc,m1pch1,m1pch2,lswpchb,lswpche,
     .            ipitmbs,jpitmbs,kpitmbs,iitmsa,nodes,
     .            dtvi,dtvj,dtvk,fbeta2,fbetav2,wold,tv1old,tv2old,
     .            m1fgt,m5fgt,iturv1,iturv2,mtorder,ntorder,iunsteady,
     .            dtphy,totptim,m1fgtrk,m5fgtrk,irkstg,nrkstg,mrkstg,
     .            rkphysa,rkphyse,delwork,dtv1ork,dtv2ork,
     .            mebdf,ibdf2opt,icyc,itime,wbar,tv1bar,tv2bar)
# endif
c
      fgrd = 1./(2.**real((igrid-1)))
      work  = work +fgrd*fgrd*fgrd
c
      call cputim (tim)
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c---- send timings and work on this grid
c
      do n=1,3
        twk1da(n) = tim(n,2)
      enddo
      twk1da(4) = work
#  ifdef BUILD_PVM
      call PVMFpsend (master,TAG_TIME,twk1da,4,RTYPE,ierr)
#  else
      if (myrank.ne.master) then
      call MPI_Send (twk1da,4,RTYPE,master,
     .               TAG_TIME,mycomm,ierr)
      end if
#  endif
c
c######################################################################
c master process
c     collect timings from the nodes
c######################################################################
c
      if (myrank.eq.master) then
c
      do n=1,3
        tim(n,2) = 0.
      enddo
      work = 0.
c
      do inode=1,nnodes
#  ifdef BUILD_PVM
        call PVMFprecv (nodes(ndlist+inode),TAG_TIME,
     .                  twk1da,4,RTYPE,
     .                  itid,itag,ilen,ierr)
#  else
      if (nodes(ndlist+inode)-1.ne.myrank) then
        call MPI_Recv (twk1da,4,RTYPE,
     .                 nodes(ndlist+inode)-1,TAG_TIME,
     .                 mycomm,istat,ierr)
      end if
#  endif
        do n=1,3
          tim(n,2) = tim(n,2) + twk1da(n) / real(nnodes)
        enddo
        work = work + twk1da(4) / real(nnodes)
      enddo
c
      endif
#endif
c
c######################################################################
c all processes
c     print results on this grid
c######################################################################
c
      do n=1,3
        tsum(n) = tsum(n) + tim(n,2)
      enddo
      totw = totw + work
c
      if (iwrit.gt.0) then
c
      icyc      = max(mcyc(ilev),2)
      rtrms0    = rtrmsg(ires0)
c dph rate1     = (rtrmsg(nres)/rtrms0)**   (1./(work  -1.))
      rate1     = (rtrmsg(nres)/rtrms0)**abs(1./(work-1.))
      rate2     = (rtrmsg(nres)/rtrms0)**(1./real(icyc  -1))
      write (iwrit,'(/,"       time        work    time/cyc",
     .                "    work/cyc       rtrms 1     rtrms 2",
     .                "   reductn/work reductn/cyc")')
      write (iwrit,'(1x,f10.2,3f12.4,2x,2e12.4,2f12.4)')
     .  tim(1,2),work,tim(1,2)/real(mcyc(ilev)),
     .  work/real(mcyc(ilev)),rtrms0,rtrmsg(nres),rate1,rate2
c dph check 2f12.4 format for rate1 (causes ********)
c
      endif
c
c------ end the loop on fmg levels
c
 1000 continue
c
c******************************************************************
c
c     totals for all grids
c
      if (iwrit.gt.0) then
c
      if (nlev.gt.1) then
        work      = totw - wtot
        icyc      = max(ncyc,2)
        rtrms0    = rtrmsg(1)
c dph   rate1     = (rtrmsg(nres)/rtrms0)**   (1./(work-1.))
        rate1     = (rtrmsg(nres)/rtrms0)**abs(1./(work-1.))
        rate2     = (rtrmsg(nres)/rtrms0)**(1./real(icyc  -1))
        write (iwrit,'(//," convergence for all grid levels")')
        write (iwrit,'(/,"       time        work    time/cyc",
     .                  "    work/cyc       rtrms 1     rtrms 2",
     .                  "   reductn/work reductn/cyc")')
        write (iwrit,'(1x,f10.2,3f12.4,2x,2e12.4,2f12.4)')
     .    tsum(1),work,tsum(1)/real(ncyc),
     .    work/real(ncyc),rtrms0,rtrmsg(nres),rate1,rate2
      endif
c
      if (mcplpr.gt.0)
     .  call rplot (iwrit,nres,nsupg,rtrmsg,wcount,resout,supout)
c 
      write (iwrit,'(//,
     .  " Force and moment results",//,
     .  3x," mach number ang attack   ref area  ",
     .  " ref length   x mom ref   y mom ref ")')
      write (iwrit,'(1x,7f12.5)') rm,al,sref,scal,xref,yref 
c
      cd        = cdtg(nres) +cdvg(nres)
      vld1      = 0.
      if (abs(cdtg(nres)).gt.1.e-6) vld1 = clg(nres)/cdtg(nres)
      vld       = 0.
      if (abs(cd).gt.1.e-6) vld = clg(nres)/cd
      write (iwrit,'(3x,"       cl      cd form   cd friction",
     .                  "       cd     l/d form       l/d    ")')
      write (iwrit,'(1x,7f12.5)') clg(nres),cdtg(nres),
     .  cdvg(nres),cd,vld1,vld
      write (iwrit,'(3x,"  cm pitch     cm roll      cm yaw      cn ")')
      write (iwrit,'(1x,7f12.5)') cmzg(nres),cmxg(nres),cmyg(nres)
     .                                                 ,cng (nres)
c
      endif
c 
c--------- dump out flowfield at the end of run
c
      if (msave.ne.0)
     .  call savsol2(imn,jmn,kmn,imp2,jmp2,kmp2,m1cc,m5cc,mxsizc,
     .               mxsiz5c,mgrlev,mbloc,nbloc,mres,nres,
     .               w,p,eomu,turv1,turv2,smin,wcount,rtrmsg,
     .               hrmsg,clg,cmzg,cdtg,cdvg,nsupg,
     .               mxszfg,wold,tv1old,tv2old,m5fgt,m1fgt,mtorder,
     .               ntorder,iturv1,iturv2,iunsteady,totptim,
# if defined BUILD_MPI
     .               scal,xref,yref,zref,sref,
     .               buffw,buffp,buffe,buff1,buff2,buffs,
     .               mx3dwk,buffwo,bufftv1,bufftv2,nodes )
# else
     .               scal,xref,yref,zref,sref,nodes)
# endif
c
c------ dump out the average flow field for unsteady flows
c
      if (mavg.gt.0) then
c
        call savsolavg(imn,jmn,kmn,imp2,jmp2,kmp2,m1fg,m5fg,
     .               mxszfg,mxsz5fg,mgrlev,mbloc,nbloc,
     .               wavg,pavg,eomuavg,wcount(nres),rtrmsg(nres),
     .               hrmsg(nres),clavg,cmzavg,cdtavg,cdvavg,nsupg(nres),
     .               iunsteady,totptim,
# if defined BUILD_MPI
     .               scal,xref,yref,zref,sref,
     .               buffw,buffp,buffe,nodes )
# else
     .               scal,xref,yref,zref,sref,nodes)
# endif
c
      endif
c
c
c---- save plot file
c
      if (mplot.gt.0 .and. myrank.eq.master)
     .  call savplt (imn,jmn,kmn,im,jm,km,
     .               mgrlev,mbloc,nbloc,mres,nres,
     .               wcount,rtrmsg,hrmsg,nsupg,
     .               clg,cmxg,cmyg,cmzg,cdtg,cdvg,czg,
     .               cl ,cmx ,cmy ,cmz ,cdt ,cdv ,cz )
c 
c---- output timings for this run
c
      call cputim (tim)
c
      if (myrank.eq.master) then
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
        write (iwrit,'(//,
     .    "                      time in seconds",/,
     .    "  node      user    system     total    wall clock")')
        write (iwrit,'(i6,3f10.2,f12.2)')
     .    0,tim(1,1),tim(2,1),tim(1,1)+tim(2,1),tim(3,1)
        do inode=2,nnodes
#  ifdef BUILD_PVM
          call PVMFprecv (nodes(ndlist+inode),TAG_TIME,
     .                    tsum,3,REAL4,
     .                    itid,itag,ilen,ierr)
#  else
      if (nodes(ndlist+inode)-1.ne.myrank) then
          call MPI_Recv (tsum,3,MPI_REAL,
     .                   nodes(ndlist+inode)-1,TAG_TIME,
     .                   mycomm,istat,ierr)
      else
        do ntim=1,3
        tsum(ntim) = tim(ntim,1)
        end do
      end if
#  endif
          write (iwrit,'(i6,3f10.2,f12.2)')
     .      inode-1,tsum(1),tsum(2),tsum(1)+tsum(2),tsum(3)
        enddo
c
cjm      else
      end if
c
#  ifdef BUILD_PVM
        call PVMFpsend (master,TAG_TIME,tim,3,REAL4,ierr)
#  else
      if (myrank.ne.master) then
        call MPI_Send (tim,3,MPI_REAL,master,
     .                 TAG_TIME,mycomm,ierr)
      end if
#  endif
#endif
        if (iwrit.gt.0) then
          write (iwrit,'(//,
     .      " timings in seconds :      user    system     total",
     .      "    wall clock")')
          write (iwrit,'(" time on all grids  :",3f10.2,f12.2)')
     .      tsum(1),tsum(2),tsum(1)+tsum(2),tsum(3)
          write (iwrit,'(" time for execution :",3f10.2,f12.2)')
     .      tim(1,1),tim(2,1),tim(1,1)+tim(2,1),tim(3,1)
        endif
c
# if defined BUILD_MPI
# else
      endif
# endif
c
      nhrs  = tim(1,1) / 3600
      nsecs = tim(1,1) - nhrs * 3600
      nmins = nsecs / 60
      nsecs = nsecs - nmins * 60
      if (iwrit.gt.0)
     .  write (iwrit,'(/," total run time = ",i4," hours ",
     .    i4," minutes ",i4," seconds")') nhrs,nmins,nsecs
c
c     stop
      end
