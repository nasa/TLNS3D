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
      subroutine jacobi (imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .  eomu,ivisb,ws,
     .  w,p,delw,vol,ft,dtphy,ntorder,
     .  dtvi,dtvj,dtvk,fbeta2)
c
c        point Jacobi preconditioning
c
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
      common/rk/    cfl,c(6),qfil(6),beta(6),vt,hm,mstage
      common/rkdis/ vis0,vis2,vis4,zeta,vepsn,vepsl,enteps,icau
      common/timek/  ct(9)

      common/vis   / diffac,rsvis,omega
      common/unscom1/ anin,aninb,bninb,betvis,betvisv
c
c-------   preconditioning related quantities      -----------------
c
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubcon2,iprec
c
      dimension    w(imp2,jmp2,kmp2,5),   ws(imp2,jmp2,kmp2,5),
     .             delw(imp2,jmp2,kmp2,5), p(imp2,jmp2,kmp2),
     .             dtvi(imp2,jmp2,kmp2),
     .             dtvj(imp2,jmp2,kmp2),dtvk(imp2,jmp2,kmp2),
     .             vol(imp2,jmp2,kmp2), eomu(imp2,jmp2,kmp2)

      dimension    fbeta2(imp2,jmp2,kmp2)
c
      dimension    six(imp1,jmp1,kmp1),
     .             siy(imp1,jmp1,kmp1), siz(imp1,jmp1,kmp1),
     .             sjx(imp1,jmp1,kmp1), sjy(imp1,jmp1,kmp1),
     .             sjz(imp1,jmp1,kmp1), skx(imp1,jmp1,kmp1),
     .             sky(imp1,jmp1,kmp1), skz(imp1,jmp1,kmp1)
c
      dimension    ivisb(3)

      if (i2dfl.eq.1) then

        if (iprec.eq.0) call jac2dp0(
     .  imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,ivisb,ws,
     .  w,p,delw,vol,ft,dtphy,ntorder,dtvi,dtvj)

        if (iprec.eq.1) call jac2dp1(
     .  imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,ivisb,ws,
     .  w,p,delw,vol,ft,dtphy,ntorder,dtvi,dtvj,fbeta2)

        if (iprec.eq.2) call jac2dp2(
     .  imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,ivisb,ws,
     .  w,p,delw,vol,ft,dtphy,ntorder,dtvi,dtvj,fbeta2)

      endif

      if (i2dfl.eq.0) then

        if (iprec.eq.0) call jac3dp0(
     .  imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,ivisb,ws,
     .  w,p,delw,vol,ft,dtphy,ntorder,dtvi,dtvj,dtvk)

        if (iprec.eq.1) call jac3dp1(
     .  imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,ivisb,ws,
     .  w,p,delw,vol,ft,dtphy,ntorder,dtvi,dtvj,dtvk,fbeta2)

        if (iprec.eq.2) call jac3dp2(
     .  imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,ivisb,ws,
     .  w,p,delw,vol,ft,dtphy,ntorder,dtvi,dtvj,dtvk,fbeta2)

      endif

      return
      end
