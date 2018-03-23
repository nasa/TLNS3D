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
      subroutine filtf(ncyc,nstage,imn,jmn,kmn,im,jm,km,
     .                 imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .                 m1cc,mgrlev,igrid,mbloc,ibloc,
     .                 mxsizc,mx1dwk,imap,msegt,msegn,nseg,p,
     .                 imx,jmx,kmx,imxp1,jmxp1,kmxp1,
     .                 imxp2,jmxp2,kmxp2,wn,pn,vol,
     .                 six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                 dti,dtj,dtk,dtl,ratioij,ratiojk,
     .                 dw,ew,fs,df2,df4,dp,ep,phi,
     .                 wk2da,wk2db,wk2dc,wk2dd,wk3d,fv,fw,
     .                 dtvi,dtvj,dtvk,fbeta2 )
c 
c 
c           controller routine for
c           artificial dissipations (for fine mesh of mg)
c           blended second and fourth diferrences
c
c           initial coding: by m.d.sanetrik (Sept, 1990)
c           latest mods: by v.n.vatsa for making it compatible 
c                        with preconditioned version (Feb. 1998)
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
c
      dimension     im  (mgrlev,mbloc), jm  (mgrlev,mbloc),
     .              km  (mgrlev,mbloc), imp1(mgrlev,mbloc),
     .              jmp1(mgrlev,mbloc), kmp1(mgrlev,mbloc),
     .              imp2(mgrlev,mbloc), jmp2(mgrlev,mbloc),
     .              kmp2(mgrlev,mbloc),
     .              imap(msegt,msegn)
c
      dimension     m1cc(mgrlev,mbloc)
c
      dimension     p(mxsizc )
c
      dimension    wn(imxp2,jmxp2,kmxp2,5),   pn(imxp2,jmxp2,kmxp2),
     .             fw(imxp2,jmxp2,kmxp2,5),   fv(imxp2,jmxp2,kmxp2),
     .            vol(imxp2,jmxp2,kmxp2),
     .            dtl(imxp2,jmxp2,kmxp2),    dti(imxp2,jmxp2,kmxp2),
     .            dtj(imxp2,jmxp2,kmxp2),    dtk(imxp2,jmxp2,kmxp2),
     .        ratioij(imxp2,jmxp2,kmxp2),ratiojk(imxp2,jmxp2,kmxp2)
c
      dimension   dtvi(imxp2,jmxp2,kmxp2),  dtvj(imxp2,jmxp2,kmxp2),
     .            dtvk(imxp2,jmxp2,kmxp2),
     .            fbeta2(imxp2,jmxp2,kmxp2)
c
      dimension   six(imxp1,jmxp1,kmxp1),siy(imxp1,jmxp1,kmxp1),
     .            siz(imxp1,jmxp1,kmxp1),sjx(imxp1,jmxp1,kmxp1),
     .            sjy(imxp1,jmxp1,kmxp1),sjz(imxp1,jmxp1,kmxp1),
     .            skx(imxp1,jmxp1,kmxp1),sky(imxp1,jmxp1,kmxp1),
     .            skz(imxp1,jmxp1,kmxp1)
c
c
c     array declaration for temporary work-space
c
c     1-d, 2-d & 3-d arrays
c
      dimension   dw(mx1dwk,5),ew(mx1dwk,5),fs(mx1dwk,5),df2(mx1dwk,5),
     .           df4(mx1dwk,5),dp(mx1dwk)  ,ep(mx1dwk)  ,phi(mx1dwk)
c
      dimension  wk2da(mx1dwk,mx1dwk),wk2db(mx1dwk,mx1dwk),
     .           wk2dc(mx1dwk,mx1dwk),wk2dd(mx1dwk,mx1dwk)
      dimension  wk3d(imxp2,jmxp2,kmxp2,5)
c
c
      call getdifi(imn,jmn,kmn,im,jm,km,
     .             imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .             m1cc,mgrlev,igrid,mbloc,ibloc,mxsizc,mx1dwk,
     .             imap,msegt,msegn,nseg,p,
     .             imx,jmx,kmx,imxp1,jmxp1,kmxp1,imxp2,jmxp2,kmxp2,
     .             wn,pn,fv,
     .             wk2da,wk2db,wk2dc,wk2dd,wk3d)
c
      call filtfx(ncyc,nstage,imn,jmn,kmn,imx,jmx,kmx,
     .            imxp1,jmxp1,kmxp1,imxp2,jmxp2,kmxp2,wn,pn,vol,
     .            six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .            dtl,dti,dtj,dtk,ratioij,ratiojk,
     .            dw,ew,fs,df2,df4,dp,ep,phi,fw,fv,
     .            mx1dwk,wk3d,dtvi,dtvj,dtvk,fbeta2)
c
      call getdifj(nstage,imn,jmn,kmn,im,jm,km,
     .             imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .             m1cc,mgrlev,igrid,mbloc,ibloc,mxsizc,mx1dwk,
     .             imap,msegt,msegn,nseg,p,
     .             imx,jmx,kmx,imxp1,jmxp1,kmxp1,imxp2,jmxp2,kmxp2,
     .             wn,pn,fv,
     .             wk2da,wk2db,wk2dc,wk2dd,wk3d  )
c
      call filtfy(ncyc,nstage,imn,jmn,kmn,imx,jmx,kmx,
     .            imxp1,jmxp1,kmxp1,imxp2,jmxp2,kmxp2,wn,pn,vol,
     .            six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .            dtl,dti,dtj,dtk,ratioij,ratiojk,
     .            dw,ew,fs,df2,df4,dp,ep,phi,fw,fv,
     .            mx1dwk,wk3d,dtvi,dtvj,dtvk,fbeta2)
c
      if(i2dfl.eq.1) return
c
      call getdifk(ncyc,nstage,imn,jmn,kmn,im,jm,km,
     .             imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .             m1cc,mgrlev,igrid,mbloc,ibloc,mxsizc,mx1dwk,
     .             imap,msegt,msegn,nseg,p,
     .             imx,jmx,kmx,imxp1,jmxp1,kmxp1,imxp2,jmxp2,kmxp2,
     .             wn,pn,fv,
     .             wk2da,wk2db,wk2dc,wk2dd,wk3d    )
c
      call filtfz(ncyc,nstage,imn,jmn,kmn,imx,jmx,kmx,
     .            imxp1,jmxp1,kmxp1,imxp2,jmxp2,kmxp2,wn,pn,vol,
     .            six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .            dtl,dti,dtj,dtk,ratioij,ratiojk,
     .            dw,ew,fs,df2,df4,dp,ep,phi,fw,fv,
     .            mx1dwk,wk3d,dtvi,dtvj,dtvk,fbeta2)
c
c
      return
      end
