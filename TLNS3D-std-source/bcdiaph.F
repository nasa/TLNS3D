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
      subroutine  bcdiaph(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                   imp2,jmp2,kmp2,w,p,eomu,turv1,
     .                   six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                   ri1,ri2,rj1,rj2,rk1,rk2,mx1dwk,dp,
     .                   nface,n1beg,n1end,n2beg,n2end,
     .             igrid,isoln,extyp,twotref,amplit,tperiod,totptim,x)
c=======================================================================
c 
c     object : to compute the ghost cell values for flow variables
c              at moving diaphragm walls
c     Note: This routine was created for a very specific geometric
c           configuration - it is not a general purpose routine
c
c     initial coding : by v.n.vatsa, nasa langley ( October, 2000)
c     last mods: by v.n.vatsa, nasa langley (  Sept. 2001)
c
c     w(i,j,k,1)  = density 
c     w(i,j,k,2)  = momentum in x direction 
c     w(i,j,k,3)  = momentum in y direction 
c     w(i,j,k,4)  = momentum in z direction 
c     w(i,j,k,5)  = total energy
c     p(i,j,k)    = pressure
c     eomu(i,j,k) = eddy-viscosity
c
c
c
      dimension      w(imp2,jmp2,kmp2,5),
     .               p(imp2,jmp2,kmp2)  ,  eomu(imp2,jmp2,kmp2),
     .                                    turv1(imp2,jmp2,kmp2)
      dimension      x(imp1,jmp1,kmp1,3)
c
      dimension      six(imp1,jmp1,kmp1), siy(imp1,jmp1,kmp1),
     .               siz(imp1,jmp1,kmp1), sjx(imp1,jmp1,kmp1),
     .               sjy(imp1,jmp1,kmp1), sjz(imp1,jmp1,kmp1),
     .               skx(imp1,jmp1,kmp1), sky(imp1,jmp1,kmp1),
     .               skz(imp1,jmp1,kmp1)
c
      dimension      ri1(jmp1,kmp1,3)   , ri2(jmp1,kmp1,3),
     .               rj1(imp1,kmp1,3)   , rj2(imp1,kmp1,3),
     .               rk1(imp1,jmp1,3)   , rk2(imp1,jmp1,3)
c
      dimension      dp(mx1dwk,mx1dwk)
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c
c
c------------------------------------------------------------------
c     extyp  =   0  zero order extrapolation for p at wall boundary
c                1  linear     extrapolation for p at wall boundary
c               -1  normal momentum for pressure at wall boundary
c
c------------------------------------------------------------------
c
c
c     extyp = 0.0
c
c      check the permissibility of the 'extyp' 
c
      if(extyp .gt.1.01 .or. extyp .lt. -1.01 )   then
         write (iwrit,'(2x," extyp = ",f8.4, " not permissible")')
         stop
      endif
c
c      define a switch for adiabatic or specified wall temperature
c      condition
c
         iadb = 0
         if (twotref.le.0.0) iadb = 1
c
c      branch-cut for different type of boundary faces
c
      go to (100,200,300,400,500,600) nface
c
c
  100 continue
c--------------------- i = imn  face ----------------------------
c
c     if( n1beg .ne. jmn) n1beg = n1beg +1
c     if( n1end .eq. jmx) n1end = n1end +1
c     if( n2beg .ne. kmn) n2beg = n2beg +1
c     if( n2end .eq. kmx) n2end = n2end +1
      n1beg = n1beg +1
      n2beg = n2beg +1
c
      i    =    imn
c
c     set the parameters for nozzle/test case
c
      pi      = acos(-1.)
c     uinfdim = rm*cinfdim
      freqnd  = 2.*pi/tperiod
      uumax   = freqnd*amplit*sin(freqnd*totptim)
cvn   uumax   = uumax*rm
      
      do 170 k= n2beg, n2end
      delty        =    x(i,n1end,k,2)-x(i,n1beg,k,2)
      do 170 j= n1beg, n1end
      eomu(i,j,k)  =    eomu(i+1,j,k)
ccc   turv1(i,j,k) =    turv1(i+1,j,k)
      yloc         =    .5*(x(i,j,k,2)+x(i,j-1,k,2))-x(i,n1beg,k,2)
      yrat         =    abs(yloc/delty)
      ugn          =    uumax*(sin(pi*yrat)**2 )
cvn   ugn          =    uumax* sin(pi*yrat)
      rhoi         =    1./w(i+1,j,k,1)
      ugc          =    2.*ugn - w(i+1,j,k,2)*rhoi
cvn   vgc          =    0.
cvn   wgc          =    0.
      vgc          =   -w(i+1,j,k,3)*rhoi
      wgc          =   -w(i+1,j,k,4)*rhoi
      w(i,j,k,1)   =    w(i+1,j,k,1)
      rhoi         =    1./w(i,j,k,1)
      tsgc         =    twotref
      p(i,j,k)     =    w(i,j,k,1)*tsgc
      w(i,j,k,2)   =    w(i,j,k,1)*ugc
      w(i,j,k,3)   =    w(i,j,k,1)*vgc
      w(i,j,k,4)   =    w(i,j,k,1)*wgc
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  170 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         i   =  imn - 1
c
         do 180 j= n1beg, n1end
         do 180 k= n2beg, n2end
         w(i,j,k,1)  = 3.*w(i+1,j,k,1) -3.*w(i+2,j,k,1) +w(i+3,j,k,1)
         w(i,j,k,2)  = 3.*w(i+1,j,k,2) -3.*w(i+2,j,k,2) +w(i+3,j,k,2)
         w(i,j,k,3)  = 3.*w(i+1,j,k,3) -3.*w(i+2,j,k,3) +w(i+3,j,k,3)
         w(i,j,k,4)  = 3.*w(i+1,j,k,4) -3.*w(i+2,j,k,4) +w(i+3,j,k,4)
         w(i,j,k,5)  = 3.*w(i+1,j,k,5) -3.*w(i+2,j,k,5) +w(i+3,j,k,5)
         p(i,j,k)    = 3.*p(i+1,j,k)   -3.*p(i+2,j,k)   +p(i+3,j,k)
         eomu(i,j,k) = eomu(i+1,j,k)
ccc      turv1(i,j,k)=turv1(i+1,j,k)
  180    continue
      endif
      return
c
c
  200 continue
c--------------------- i = imx  face ----------------------------
c
c     if( n1beg .ne. jmn) n1beg = n1beg +1
c     if( n1end .eq. jmx) n1end = n1end +1
c     if( n2beg .ne. kmn) n2beg = n2beg +1
c     if( n2end .eq. kmx) n2end = n2end +1
c
      n1beg = n1beg +1
      n2beg = n2beg +1

      i   =   imx + 1
c
c     set the parameters specific to the nozzle/test case
c     (distances are in meters, speeds in m/sec)
c
      pi      = acos(-1.)
      cinfdim = 340.
      dispmx  = 0.00058
      dnozzle = 0.000508
      hertz   = 500.
      uinfdim = rm*cinfdim
      amplit  = 2.*pi*hertz*dispmx/uinfdim
      freqnd  = 2.*pi*hertz*dnozzle/uinfdim
      uumax   = amplit*sin(freqnd*totptim)*rm*sqrt(gamma)
      
      do 270 k= n2beg, n2end
      delty        =    x(i-1,n1end,k,2)-x(i-1,n1beg,k,2)
      do 270 j= n1beg, n1end
      eomu(i,j,k)  =    eomu(i-1,j,k)
ccc   turv1(i,j,k) =    turv1(i-1,j,k)
      yloc         =    .5*(x(i-1,j,k,2)+x(i-1,j-1,k,2))
     .             -        x(i-1,n1beg,k,2)
      yrat         =    abs(yloc/delty)
      ugn          =   -uumax*(sin(pi*yrat)**2 )
      ugc          =    2.*ugn - w(i-1,j,k,2)/w(i-1,j,k,1)
      vgc          =    0.
      wgc          =    0.
      w(i,j,k,1)   =    w(i-1,j,k,1)
      rhoi         =    1./w(i,j,k,1)
      tsgc         =    twotref
      p(i,j,k)     =    w(i,j,k,1)*tsgc
      w(i,j,k,2)   =    w(i,j,k,1)*ugc
      w(i,j,k,3)   =    w(i,j,k,1)*vgc
      w(i,j,k,4)   =    w(i,j,k,1)*wgc
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  270 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         i   =  imx + 2
c
         do 280 j= n1beg, n1end
         do 280 k= n2beg, n2end
         w(i,j,k,1)  = 3.*w(i-1,j,k,1) -3.*w(i-2,j,k,1) +w(i-3,j,k,1)
         w(i,j,k,2)  = 3.*w(i-1,j,k,2) -3.*w(i-2,j,k,2) +w(i-3,j,k,2)
         w(i,j,k,3)  = 3.*w(i-1,j,k,3) -3.*w(i-2,j,k,3) +w(i-3,j,k,3)
         w(i,j,k,4)  = 3.*w(i-1,j,k,4) -3.*w(i-2,j,k,4) +w(i-3,j,k,4)
         w(i,j,k,5)  = 3.*w(i-1,j,k,5) -3.*w(i-2,j,k,5) +w(i-3,j,k,5)
         p(i,j,k)    = 3.*p(i-1,j,k)   -3.*p(i-2,j,k)   +p(i-3,j,k)
         eomu(i,j,k) = eomu(i-1,j,k)
ccc      turv1(i,j,k)=turv1(i-1,j,k)
  280    continue
      endif
      return
c
c
  300 continue
c--------------------- j = jmn  face ----------------------------
c
      if( n1beg .ne. kmn) n1beg = n1beg +1
      if( n1end .eq. kmx) n1end = n1end +1
      if( n2beg .ne. imn) n2beg = n2beg +1
      if( n2end .eq. imx) n2end = n2end +1
c
      j   =   jmn
c
      if (extyp .lt. 0.) go to 350
      if (extyp .gt. 0.) go to 340
c
c--------     zero order extrapolation for surface pressure
      do 330 k= n1beg, n1end
      do 330 i= n2beg, n2end
      dp(i,k)      = 0.
 330  continue
      go to 369
c
c--------     linear extrapolation for surface pressure
 340  continue
      do 345 k= n1beg, n1end
      do 345 i= n2beg, n2end
      dp(i,k)      = p(i,j+2,k) - p(i,j+1,k)
 345  continue
      go to 369
 350  continue
c
c--------     normal momentum equation for surface pressure
c
      do 355 k= n1beg, n1end
      kp      = k + 1
      km      = k - 1
      do 355 i= n2beg, n2end
c
      sixa      = (six(i,j+1,k)  +six(i-1,j+1,k))
      siya      = (siy(i,j+1,k)  +siy(i-1,j+1,k))
      siza      = (siz(i,j+1,k)  +siz(i-1,j+1,k))
      sjxa      = 2.*sjx(i,j,k)
      sjya      = 2.*sjy(i,j,k)
      sjza      = 2.*sjz(i,j,k)
      skxa      = (skx(i,j+1,k)  +skx(i,j+1,k-1))
      skya      = (sky(i,j+1,k)  +sky(i,j+1,k-1))
      skza      = (skz(i,j+1,k)  +skz(i,j+1,k-1))
      rxi       = .5*(rj1(i+1,k,1)   -rj1(i-1,k,1)) 
      ryi       = .5*(rj1(i+1,k,2)   -rj1(i-1,k,2)) 
      rzi       = .5*(rj1(i+1,k,3)   -rj1(i-1,k,3)) 
      pi        = .5*(p(i+1,j+1,k)     -p(i-1,j+1,k))
      rxk       = .5*(rj1(i,kp,1)    -rj1(i,km,1))
      ryk       = .5*(rj1(i,kp,2)    -rj1(i,km,2))
      rzk       = .5*(rj1(i,kp,3)    -rj1(i,km,3))
      pk        = .5*(p(i,j+1,kp)      -p(i,j+1,km))
      si        = rj1(i,k,1)*sixa  +rj1(i,k,2)*siya  +rj1(i,k,3)*siza 
      sj        = rj1(i,k,1)*sjxa  +rj1(i,k,2)*sjya  +rj1(i,k,3)*sjza 
      sk        = rj1(i,k,1)*skxa  +rj1(i,k,2)*skya  +rj1(i,k,3)*skza 
      qi        = w(i,j+1,k,2)*sixa  +w(i,j+1,k,3)*siya
     .          + w(i,j+1,k,4)*siza
      qk        = w(i,j+1,k,2)*skxa  +w(i,j+1,k,3)*skya 
     .          + w(i,j+1,k,4)*skza
      dp(i,k)   = (qi*(w(i,j+1,k,2)*rxi  +w(i,j+1,k,3)*ryi
     .                +w(i,j+1,k,4)*rzi) /w(i,j+1,k,1) 
     .            +qk*(w(i,j+1,k,2)*rxk  +w(i,j+1,k,3)*ryk
     .                +w(i,j+1,k,4)*rzk) /w(i,j+1,k,1) 
     .            -si*pi  -sk*pk)/sj
  355 continue
      if (extyp.eq.-1.) go to 369
c
c--------- improved extrapolation using quadratic fit 
c          and normal pressure gradient at the surface
c
      do 360 k= n1beg, n1end
      do 360 i= n2beg,n2end
      dp(i,k)   = .75*dp(i,k) + .25*(p(i,j+2,k) - p(i,j+1,k))
  360 continue
c
c
  369 continue
      do 370 k= n1beg, n1end
      do 370 i= n2beg, n2end
ccc   w(i,j,k,1)   =    w(i,j+1,k,1)
ccc   w(i,j,k,2)   =   -w(i,j+1,k,2)
ccc   w(i,j,k,3)   =   -w(i,j+1,k,3)
ccc   w(i,j,k,4)   =   -w(i,j+1,k,4)
      eomu(i,j,k)  =    eomu(i,j+1,k)
ccc   turv1(i,j,k) =   turv1(i,j+1,k)
      p(i,j,k)     =    dim(p(i,j+1,k),dp(i,k))
      rhoi         =    1./w(i,j+1,k,1)
      ugc          =   -w(i,j+1,k,2)*rhoi
      vgc          =   -w(i,j+1,k,3)*rhoi
      wgc          =   -w(i,j+1,k,4)*rhoi
      temp         =    p(i,j+1,k)*rhoi
      tsgc         =     max((2.*twotref - temp),0.8*twotref)*(1-iadb)
     .             +    iadb*temp
      p(i,j,k)     =     max(p(i,j,k),0.001)
      w(i,j,k,1)   =    p(i,j,k)/tsgc
      w(i,j,k,2)   =    w(i,j,k,1)*ugc
      w(i,j,k,3)   =    w(i,j,k,1)*vgc
      w(i,j,k,4)   =    w(i,j,k,1)*wgc
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  370 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         j   =  jmn - 1
c
         do 380 k= n1beg, n1end
         do 380 i= n2beg, n2end
         w(i,j,k,1)  = 3.*w(i,j+1,k,1) -3.*w(i,j+2,k,1) +w(i,j+3,k,1)
         w(i,j,k,2)  = 3.*w(i,j+1,k,2) -3.*w(i,j+2,k,2) +w(i,j+3,k,2)
         w(i,j,k,3)  = 3.*w(i,j+1,k,3) -3.*w(i,j+2,k,3) +w(i,j+3,k,3)
         w(i,j,k,4)  = 3.*w(i,j+1,k,4) -3.*w(i,j+2,k,4) +w(i,j+3,k,4)
         w(i,j,k,5)  = 3.*w(i,j+1,k,5) -3.*w(i,j+2,k,5) +w(i,j+3,k,5)
         p(i,j,k)    = 3.*p(i,j+1,k)   -3.*p(i,j+2,k)   +p(i,j+3,k)
         p(i,j,k)    =     max(p(i,j,k),0.001)
         eomu(i,j,k) = eomu(i,j+1,k)
  380    continue
c
c        do 380 k= n1beg, n1end
c        do 380 i= n2beg, n2end
c        w(i,j,k,1)   =    w(i,j+3,k,1)
c        w(i,j,k,2)   =   -w(i,j+3,k,2)
c        w(i,j,k,3)   =   -w(i,j+3,k,3)
c        w(i,j,k,4)   =   -w(i,j+3,k,4)
c        w(i,j,k,5)   =    w(i,j+3,k,5)
c        eomu(i,j,k)  =    eomu(i,j+3,k)
c        turv1(i,j,k) =   turv1(i,j+3,k)
c        p(i,j,k)     =    dim(p(i,j+1,k),dp(i,k))
c        vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
c    .                     w(i,j,k,4)**2
c        w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
c 380    continue
      endif
      return
c
c
  400 continue
c--------------------- j = jmx  face ----------------------------
c
      if( n1beg .ne. kmn) n1beg = n1beg +1
      if( n1end .eq. kmx) n1end = n1end +1
      if( n2beg .ne. imn) n2beg = n2beg +1
      if( n2end .eq. imx) n2end = n2end +1
c
      j   =   jmx + 1
c
      if (extyp .lt. 0.) go to 450
      if (extyp .gt. 0.) go to 440
c
c--------     zero order extrapolation for surface pressure
      do 430 k= n1beg, n1end
      do 430 i= n2beg, n2end
      dp(i,k)      = 0.
 430  continue
      go to 469
c
c--------     linear extrapolation for surface pressure
 440  continue
      do 445 k= n1beg, n1end
      do 445 i= n2beg, n2end
      dp(i,k)      = p(i,j-2,k) - p(i,j-1,k)
 445  continue
      go to 469
 450  continue
c
c--------     normal momentum equation for surface pressure
c
      do 455 k= n1beg, n1end
      kp      = k + 1
      km      = k - 1
      do 455 i= n2beg, n2end
c
      sixa      = (six(i,j-1,k)  +six(i-1,j-1,k))
      siya      = (siy(i,j-1,k)  +siy(i-1,j-1,k))
      siza      = (siz(i,j-1,k)  +siz(i-1,j-1,k))
      sjxa      = 2.*sjx(i,j,k)
      sjya      = 2.*sjy(i,j,k)
      sjza      = 2.*sjz(i,j,k)
      skxa      = (skx(i,j-1,k)  +skx(i,j-1,k-1))
      skya      = (sky(i,j-1,k)  +sky(i,j-1,k-1))
      skza      = (skz(i,j-1,k)  +skz(i,j-1,k-1))
      rxi       = .5*(rj2(i+1,k,1)   -rj2(i-1,k,1)) 
      ryi       = .5*(rj2(i+1,k,2)   -rj2(i-1,k,2)) 
      rzi       = .5*(rj2(i+1,k,3)   -rj2(i-1,k,3)) 
      pi        = .5*(p(i+1,j-1,k)     -p(i-1,j-1,k))
      rxk       = .5*(rj2(i,kp,1)    -rj2(i,km,1))
      ryk       = .5*(rj2(i,kp,2)    -rj2(i,km,2))
      rzk       = .5*(rj2(i,kp,3)    -rj2(i,km,3))
      pk        = .5*(p(i,j-1,kp)      -p(i,j-1,km))
      si        = rj2(i,k,1)*sixa  +rj2(i,k,2)*siya  +rj2(i,k,3)*siza 
      sj        = rj2(i,k,1)*sjxa  +rj2(i,k,2)*sjya  +rj2(i,k,3)*sjza 
      sk        = rj2(i,k,1)*skxa  +rj2(i,k,2)*skya  +rj2(i,k,3)*skza 
      qi        = w(i,j-1,k,2)*sixa  +w(i,j-1,k,3)*siya
     .          + w(i,j-1,k,4)*siza
      qk        = w(i,j-1,k,2)*skxa  +w(i,j-1,k,3)*skya
     .          + w(i,j-1,k,4)*skza
      dp(i,k)   = (qi*(w(i,j-1,k,2)*rxi  +w(i,j-1,k,3)*ryi
     .                +w(i,j-1,k,4)*rzi) /w(i,j-1,k,1) 
     .            +qk*(w(i,j-1,k,2)*rxk  +w(i,j-1,k,3)*ryk
     .                +w(i,j-1,k,4)*rzk) /w(i,j-1,k,1) 
     .            -si*pi  -sk*pk)/sj
  455 continue
      if (extyp.eq.-1.) go to 469
c
c--------- improved extrapolation using quadratic fit 
c          and normal pressure gradient at the surface
c
      do 460 k= n1beg, n1end
      do 460 i= n2beg,n2end
      dp(i,k)   = .75*dp(i,k) + .25*(p(i,j-2,k) - p(i,j-1,k))
  460 continue
c
c
  469 continue
      do 470 k= n1beg, n1end
      do 470 i= n2beg, n2end
ccc   w(i,j,k,1)   =    w(i,j-1,k,1)
ccc   w(i,j,k,2)   =   -w(i,j-1,k,2)
ccc   w(i,j,k,3)   =   -w(i,j-1,k,3)
ccc   w(i,j,k,4)   =   -w(i,j-1,k,4)
      eomu(i,j,k)  =    eomu(i,j-1,k)
ccc   turv1(i,j,k) =   turv1(i,j-1,k)
      p(i,j,k)     =    dim(p(i,j-1,k),dp(i,k))
      rhoi         =    1./w(i,j-1,k,1)
      ugc          =   -w(i,j-1,k,2)*rhoi
      vgc          =   -w(i,j-1,k,3)*rhoi
      wgc          =   -w(i,j-1,k,4)*rhoi
      temp         =    p(i,j-1,k)*rhoi
      tsgc         =     max((2.*twotref - temp),0.8*twotref)*(1-iadb)
     .             +    iadb*temp
      p(i,j,k)     =     max(p(i,j,k),0.001)
      w(i,j,k,1)   =    p(i,j,k)/tsgc
      w(i,j,k,2)   =    w(i,j,k,1)*ugc
      w(i,j,k,3)   =    w(i,j,k,1)*vgc
      w(i,j,k,4)   =    w(i,j,k,1)*wgc
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  470 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         j   =  jmx + 2
c
         do 480 k= n1beg, n1end
         do 480 i= n2beg, n2end
         w(i,j,k,1)  = 3.*w(i,j-1,k,1) -3.*w(i,j-2,k,1) +w(i,j-3,k,1)
         w(i,j,k,2)  = 3.*w(i,j-1,k,2) -3.*w(i,j-2,k,2) +w(i,j-3,k,2)
         w(i,j,k,3)  = 3.*w(i,j-1,k,3) -3.*w(i,j-2,k,3) +w(i,j-3,k,3)
         w(i,j,k,4)  = 3.*w(i,j-1,k,4) -3.*w(i,j-2,k,4) +w(i,j-3,k,4)
         w(i,j,k,5)  = 3.*w(i,j-1,k,5) -3.*w(i,j-2,k,5) +w(i,j-3,k,5)
         p(i,j,k)    = 3.*p(i,j-1,k)   -3.*p(i,j-2,k)   +p(i,j-3,k)
         eomu(i,j,k) = eomu(i,j-1,k)
ccc      turv1(i,j,k)=turv1(i,j-1,k)
  480    continue
      endif
      return
c
c
  500 continue
c--------------------- k = kmn  face ----------------------------
c
      if( n1beg .ne. imn) n1beg = n1beg +1
      if( n1end .eq. imx) n1end = n1end +1
      if( n2beg .ne. jmn) n2beg = n2beg +1
      if( n2end .eq. jmx) n2end = n2end +1
c
      k   =   kmn
c
      if (extyp .lt. 0.) go to 550
      if (extyp .gt. 0.) go to 540
c
c--------     zero order extrapolation for surface pressure
      do 530 j= n2beg, n2end
      do 530 i= n1beg, n1end
      dp(i,j)      = 0.
 530  continue
      go to 569
c
c--------     linear extrapolation for surface pressure
 540  continue
      do 545 j= n2beg, n2end
      do 545 i= n1beg, n1end
      dp(i,j)      = p(i,j,k+2) - p(i,j,k+1)
 545  continue
      go to 569
c
 550  continue
c
c--------- extrapolation using normal momentum equation at surface
c
      do 555 j= n2beg, n2end
      jp        = j + 1
      jm        = j - 1
      do 555 i= n1beg, n1end
      ip        = i + 1
      im        = i - 1
      sixa      = (six(i,j,k+1) + six(i-1,j,k+1))
      siya      = (siy(i,j,k+1) + siy(i-1,j,k+1))
      siza      = (siz(i,j,k+1) + siz(i-1,j,k+1))
      sjxa      = (sjx(i,j,k+1) + sjx(i,j-1,k+1))
      sjya      = (sjy(i,j,k+1) + sjy(i,j-1,k+1))
      sjza      = (sjz(i,j,k+1) + sjz(i,j-1,k+1))
      skxa      = 2.*skx(i,j,k)
      skya      = 2.*sky(i,j,k)
      skza      = 2.*skz(i,j,k)
      rxi       = .5*(rk1(ip,j,1) - rk1(im,j,1))
      ryi       = .5*(rk1(ip,j,2) - rk1(im,j,2))
      rzi       = .5*(rk1(ip,j,3) - rk1(im,j,3))
      pi        = .5*(p(ip,j,k+1)   - p(im,j,k+1))
      rxj       = .5*(rk1(i,jp,1) - rk1(i,jm,1))
      ryj       = .5*(rk1(i,jp,2) - rk1(i,jm,2))
      rzj       = .5*(rk1(i,jp,3) - rk1(i,jm,3))
      pj        = .5*(p(i,jp,k+1)   - p(i,jm,k+1))
      si        = rk1(i,j,1)*sixa + rk1(i,j,2)*siya + rk1(i,j,3)*siza 
      sj        = rk1(i,j,1)*sjxa + rk1(i,j,2)*sjya + rk1(i,j,3)*sjza 
      sk        = rk1(i,j,1)*skxa + rk1(i,j,2)*skya + rk1(i,j,3)*skza 
      qi        = w(i,j,k+1,2)*sixa + w(i,j,k+1,3)*siya 
     .          + w(i,j,k+1,4)*siza
      qj        = w(i,j,k+1,2)*sjxa + w(i,j,k+1,3)*sjya
     .          + w(i,j,k+1,4)*sjza
      dp(i,j)   = (qi*(w(i,j,k+1,2)*rxi + w(i,j,k+1,3)*ryi
     .                +w(i,j,k+1,4)*rzi)/ w(i,j,k+1,1) 
     .            +qj*(w(i,j,k+1,2)*rxj + w(i,j,k+1,3)*ryj
     .                +w(i,j,k+1,4)*rzj)/ w(i,j,k+1,1) 
     .            -si*pi - sj*pj)/sk
  555 continue
      if (extyp .eq. -1.) go to 569
c
c--------- improved extrapolation using quadratic fit 
c          and normal pressure gradient at the surface
c
      do 560 j= n2beg,n2end
      do 560 i= n1beg,n1end
      dp(i,j)   = .75*dp(i,j) + .25*(p(i,j,k+2) - p(i,j,k+1))
  560 continue
c
  569 continue
c
      do 570 i= n1beg, n1end
      do 570 j= n2beg, n2end
ccc   w(i,j,k,1)   =    w(i,j,k+1,1)
ccc   w(i,j,k,2)   =   -w(i,j,k+1,2)
ccc   w(i,j,k,3)   =   -w(i,j,k+1,3)
ccc   w(i,j,k,4)   =   -w(i,j,k+1,4)
      eomu(i,j,k)  =    eomu(i,j,k+1)
ccc   turv1(i,j,k) =   turv1(i,j,k+1)
      p(i,j,k)     =    dim(p(i,j,k+1),dp(i,j))
      rhoi         =    1./w(i,j,k+1,1)
      ugc          =   -w(i,j,k+1,2)*rhoi
      vgc          =   -w(i,j,k+1,3)*rhoi
      wgc          =   -w(i,j,k+1,4)*rhoi
      temp         =    p(i,j,k+1)*rhoi
      tsgc         =     max((2.*twotref - temp),0.8*twotref)*(1-iadb)
     .             +    iadb*temp
      p(i,j,k)     =     max(p(i,j,k),0.001)
      w(i,j,k,1)   =    p(i,j,k)/tsgc
      w(i,j,k,2)   =    w(i,j,k,1)*ugc
      w(i,j,k,3)   =    w(i,j,k,1)*vgc
      w(i,j,k,4)   =    w(i,j,k,1)*wgc
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  570 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         k   =  kmn - 1
c
         do 580 i= n1beg, n1end
         do 580 j= n2beg, n2end
         w(i,j,k,1)  = 3.*w(i,j,k+1,1) -3.*w(i,j,k+2,1) +w(i,j,k+3,1)
         w(i,j,k,2)  = 3.*w(i,j,k+1,2) -3.*w(i,j,k+2,2) +w(i,j,k+3,2)
         w(i,j,k,3)  = 3.*w(i,j,k+1,3) -3.*w(i,j,k+2,3) +w(i,j,k+3,3)
         w(i,j,k,4)  = 3.*w(i,j,k+1,4) -3.*w(i,j,k+2,4) +w(i,j,k+3,4)
         w(i,j,k,5)  = 3.*w(i,j,k+1,5) -3.*w(i,j,k+2,5) +w(i,j,k+3,5)
         p(i,j,k)    = 3.*p(i,j,k+1)   -3.*p(i,j,k+2)   +p(i,j,k+3)
         eomu(i,j,k) = eomu(i,j,k+1)
ccc      turv1(i,j,k)=turv1(i,j,k+1)
  580    continue
      endif
      return
c
c
  600 continue
c--------------------- k = kmx  face ----------------------------
c
      if( n1beg .ne. imn) n1beg = n1beg +1
      if( n1end .eq. imx) n1end = n1end +1
      if( n2beg .ne. jmn) n2beg = n2beg +1
      if( n2end .eq. jmx) n2end = n2end +1
c
      k   =   kmx + 1
c
      if (extyp .lt. 0.) go to 650
      if (extyp .gt. 0.) go to 640
c
c--------     zero order extrapolation for surface pressure
      do 630 j= n2beg, n2end
      do 630 i= n1beg, n1end
      dp(i,j)      = 0.
 630  continue
      go to 669
c
c--------     linear extrapolation for surface pressure
 640  continue
      do 645 j= n2beg, n2end
      do 645 i= n1beg, n1end
      dp(i,j)      = p(i,j,k-2) - p(i,j,k-1)
 645  continue
      go to 669
c
 650  continue
c
c--------- extrapolation using normal momentum equation at surface
c
      do 655 j= n2beg, n2end
      jp        = j + 1
      jm        = j - 1
      do 655 i= n1beg, n1end
      ip        = i + 1
      im        = i - 1
      sixa      = (six(i,j,k-1) + six(i-1,j,k-1))
      siya      = (siy(i,j,k-1) + siy(i-1,j,k-1))
      siza      = (siz(i,j,k-1) + siz(i-1,j,k-1))
      sjxa      = (sjx(i,j,k-1) + sjx(i,j-1,k-1))
      sjya      = (sjy(i,j,k-1) + sjy(i,j-1,k-1))
      sjza      = (sjz(i,j,k-1) + sjz(i,j-1,k-1))
      skxa      = 2.*skx(i,j,k)
      skya      = 2.*sky(i,j,k)
      skza      = 2.*skz(i,j,k)
      rxi       = .5*(rk2(ip,j,1) - rk2(im,j,1))
      ryi       = .5*(rk2(ip,j,2) - rk2(im,j,2))
      rzi       = .5*(rk2(ip,j,3) - rk2(im,j,3))
      pi        = .5*(p(ip,j,k-1)   - p(im,j,k-1))
      rxj       = .5*(rk2(i,jp,1) - rk2(i,jm,1))
      ryj       = .5*(rk2(i,jp,2) - rk2(i,jm,2))
      rzj       = .5*(rk2(i,jp,3) - rk2(i,jm,3))
      pj        = .5*(p(i,jp,k-1)   - p(i,jm,k-1))
      si        = rk2(i,j,1)*sixa + rk2(i,j,2)*siya + rk2(i,j,3)*siza 
      sj        = rk2(i,j,1)*sjxa + rk2(i,j,2)*sjya + rk2(i,j,3)*sjza 
      sk        = rk2(i,j,1)*skxa + rk2(i,j,2)*skya + rk2(i,j,3)*skza 
      qi        = w(i,j,k-1,2)*sixa + w(i,j,k-1,3)*siya 
     .          + w(i,j,k-1,4)*siza
      qj        = w(i,j,k-1,2)*sjxa + w(i,j,k-1,3)*sjya
     .          + w(i,j,k-1,4)*sjza
      dp(i,j)   = (qi*(w(i,j,k-1,2)*rxi + w(i,j,k-1,3)*ryi
     .                +w(i,j,k-1,4)*rzi)/ w(i,j,k-1,1) 
     .            +qj*(w(i,j,k-1,2)*rxj + w(i,j,k-1,3)*ryj
     .                +w(i,j,k-1,4)*rzj)/ w(i,j,k-1,1) 
     .            -si*pi - sj*pj)/sk
  655 continue
      if (extyp .eq. -1.) go to 669
c
c--------- improved extrapolation using quadratic fit 
c          and normal pressure gradient at the surface
c
      do 660 j= n2beg,n2end
      do 660 i= n1beg,n1end
      dp(i,j)   = .75*dp(i,j) + .25*(p(i,j,k-2) - p(i,j,k-1))
  660 continue
c
c
  669 continue
c
      do 670 i= n1beg, n1end
      do 670 j= n2beg, n2end
      w(i,j,k,1)   =    w(i,j,k-1,1)
      w(i,j,k,2)   =   -w(i,j,k-1,2)
      w(i,j,k,3)   =   -w(i,j,k-1,3)
      w(i,j,k,4)   =   -w(i,j,k-1,4)
c     w(i,j,k,5)   =    w(i,j,k-1,5)
      eomu(i,j,k)  =    eomu(i,j,k-1)
ccc   turv1(i,j,k) =   turv1(i,j,k-1)
      p(i,j,k)     =    dim(p(i,j,k-1),dp(i,j))
      rhoi         =    1./w(i,j,k-1,1)
      ugc          =   -w(i,j,k-1,2)*rhoi
      vgc          =   -w(i,j,k-1,3)*rhoi
      wgc          =   -w(i,j,k-1,4)*rhoi
      temp         =    p(i,j,k-1)*rhoi
      tsgc         =     max((2.*twotref - temp),0.8*twotref)*(1-iadb)
     .             +    iadb*temp
      p(i,j,k)     =     max(p(i,j,k),0.001)
      w(i,j,k,1)   =    p(i,j,k)/tsgc
      w(i,j,k,2)   =    w(i,j,k,1)*ugc
      w(i,j,k,3)   =    w(i,j,k,1)*vgc
      w(i,j,k,4)   =    w(i,j,k,1)*wgc
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  670 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         k   =  kmx + 2
c
         do 680 i= n1beg, n1end
         do 680 j= n2beg, n2end
         w(i,j,k,1)  = 3.*w(i,j,k-1,1) -3.*w(i,j,k-2,1) +w(i,j,k-3,1)
         w(i,j,k,2)  = 3.*w(i,j,k-1,2) -3.*w(i,j,k-2,2) +w(i,j,k-3,2)
         w(i,j,k,3)  = 3.*w(i,j,k-1,3) -3.*w(i,j,k-2,3) +w(i,j,k-3,3)
         w(i,j,k,4)  = 3.*w(i,j,k-1,4) -3.*w(i,j,k-2,4) +w(i,j,k-3,4)
         w(i,j,k,5)  = 3.*w(i,j,k-1,5) -3.*w(i,j,k-2,5) +w(i,j,k-3,5)
         p(i,j,k)    = 3.*p(i,j,k-1)   -3.*p(i,j,k-2)   +p(i,j,k-3)
         eomu(i,j,k) = eomu(i,j,k-1)
ccc      turv1(i,j,k)=turv1(i,j,k-1)
  680    continue
      endif
c
c
      return
      end
