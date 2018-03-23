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
      subroutine  bcjetfl(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                    imp2,jmp2,kmp2,w,p,eomu,turv1,
     .                    ri1,ri2,rj1,rj2,rk1,rk2,
     .                    nface,n1beg,n1end,n2beg,n2end,
     .                    igrid,isoln,extyp,
     .                    ptjopfs,ttjotfs,rmjet,alfaj,betaj )
c=======================================================================
c 
c     object : to compute the ghost cell values for flow variables
c              at nozzle/jet inflow/outflow boundary 
c              presently Riemann invariants based on 1-d
c              characteristics are used
c
c     initial coding : by v.n.vatsa, nasa langley ( Sept., 1993)
c     revision       : 10-15-93, to enforce pt at inflow
c
c     latest version : Jan. 1998, to allow multiple jet/nozzle
c     streams. The parameters ptjopfs,ttjotfs,rmjet,alfaj,betaj
c     are read in mapfile as const1,..... const5
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
c
c
      dimension      ri1(jmp1,kmp1,3)   , ri2(jmp1,kmp1,3),
     .               rj1(imp1,kmp1,3)   , rj2(imp1,kmp1,3),
     .               rk1(imp1,jmp1,3)   , rk2(imp1,jmp1,3)
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
      common/trbfs/ anutinf,zkinf,ominf
c
c
c     ptjopfs    : ratio of total pressure of jet and 
c                  free-stream (reference) static pressure (input)
c     ttjotfs    : ratio of total temperatures in jet and 
c                  free-stream (reference) static temperature (input)
c   
c     rmjet      : inflow jet Mach number
c
c     alfaj,betaj: jet inclination angles at upstream (input)
c     
c
c------------------------------------------------------------------
c
c     compute free stream conditions upstream of jet
c
      ptotjet=  p0*ptjopfs
      if (ptjopfs.eq.0.0) ptotjet = p0*(1.+gm1*rm*rm*.5)**gogm1
      t0     =  p0/rho0
      pjet   =  ptotjet/((1.+.5*gm1*rmjet*rmjet)**gogm1)
      tjet   =  t0*ttjotfs/(1.+.5*gm1*rmjet*rmjet)
      rhojet =  pjet/tjet
      cjet   =  sqrt(gamma*pjet/rhojet)
      eijet  =  pjet/(gm1*rhojet)
c
c
      ujet   =  rmjet*cjet*cos(alfaj)*cos(betaj)
      vjet   =  rmjet*cjet*sin(alfaj)*cos(betaj)
      wjet   =  rmjet*cjet*sin(betaj)
      hjet   =  gamma*eijet + .5*(ujet**2+vjet**2+wjet**2)
      sjet  =    rhojet**gamma/pjet
c--------------------------------------------------------------------
c
      evinf     =  0.0
c
      if(iturb.eq.2) then
         chi    =  anutinf
         fv1    =  chi**3/(chi**3 + 7.1**3)
         evinf  =  fv1*anutinf
      endif
c
      if(iturb.eq.3) then
         evinf  =  zkinf/ominf
      endif
c
c      branch-cut for different type of boundary faces
c
      go to (100,200,300,400,500,600) nface
c
c
  100 continue
c--------------------- i = imn  face ----------------------------
c
      if( n1beg .ne. jmn) n1beg = n1beg +1
      if( n1end .eq. jmx) n1end = n1end +1
      if( n2beg .ne. kmn) n2beg = n2beg +1
      if( n2end .eq. kmx) n2end = n2end +1
c
      i    =    imn
c
c
      if(rmjet.lt.1.0) then
         do 110 k= n2beg, n2end
cdir$ ivdep
         do 110 j= n1beg, n1end
         uf        = w(i+1,j,k,2)/w(i+1,j,k,1) 
         vf        = w(i+1,j,k,3)/w(i+1,j,k,1) 
         wf        = w(i+1,j,k,4)/w(i+1,j,k,1) 
         qnf       = uf*ri1(j,k,1)  +vf*ri1(j,k,2)  +wf*ri1(j,k,3)
         qn0       = ujet*ri1(j,k,1)+vjet*ri1(j,k,2)+wjet*ri1(j,k,3)
         c         = sqrt(gamma*p(i+1,j,k)/w(i+1,j,k,1)) 
         er        = qnf - 2.*dgm1*c
         fr        = qn0 + 2.*dgm1*cjet
         c         = .25*(er - fr)*gm1
         qn        = .50*(er + fr)
c 
c--------- inflow boundary
c 
         if (qn.gt.0.) then
           uf        = ujet + (qn - qn0)*ri1(j,k,1)
           vf        = vjet + (qn - qn0)*ri1(j,k,2)
           wf        = wjet + (qn - qn0)*ri1(j,k,3)
           s         = sjet
c 
c--------- outflow boundary 
c 
         else
           uf        = uf + (qn - qnf)*ri1(j,k,1)
           vf        = vf + (qn - qnf)*ri1(j,k,2)
           wf        = wf + (qn - qnf)*ri1(j,k,3)
           s         = w(i+1,j,k,1)**gamma/p(i+1,j,k)
         endif 
c 
         cc          = c*c/gamma
ccc      w(i,j,k,1)  = (s*cc)**dgm1
         vtots       =  uf**2 + vf**2 + wf**2
         rmloc       = sqrt(vtots/c**2)
         p(i,j,k)    = ptotjet/((1.+.5*gm1*rmloc*rmloc)**gogm1)
ccc      tloc        = t0*ttjotfs/(1.+.5*gm1*rmloc*rmloc)
ccc      w(i,j,k,1)  = p(i,j,k)/tloc
         w(i,j,k,1)  = p(i,j,k)/cc
         w(i,j,k,2)  = w(i,j,k,1)*uf
         w(i,j,k,3)  = w(i,j,k,1)*vf
         w(i,j,k,4)  = w(i,j,k,1)*wf
ccc      p(i,j,k)    = w(i,j,k,1)*cc
c
c        energy is overspecified
c     w(i,j,k,5)  = -p(i,j,k) + h0*w(i,j,k,1)
c
         vtots       =     w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                     w(i,j,k,4)**2
         w(i,j,k,5)  =     p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
         eomu(i,j,k) =     evinf
  110    continue
      else
c        supersonic flow condition
         do 115 k= n2beg, n2end
         do 115 j= n1beg, n1end
         w(i,j,k,1)  = rhojet
         w(i,j,k,2)  = rhojet*ujet
         w(i,j,k,3)  = rhojet*vjet
         w(i,j,k,4)  = rhojet*wjet
         w(i,j,k,5)  = rhojet*hjet -pjet
         p(i,j,k)    = pjet
         eomu(i,j,k) = evinf
  115    continue
      end if
c
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         i   =  imn - 1
c
         do 120 k= n2beg, n2end
         do 120 j= n1beg, n1end
         w(i,j,k,1)  = 2.*w(i+1,j,k,1) -w(i+2,j,k,1)
         w(i,j,k,2)  = 2.*w(i+1,j,k,2) -w(i+2,j,k,2)
         w(i,j,k,3)  = 2.*w(i+1,j,k,3) -w(i+2,j,k,3)
         w(i,j,k,4)  = 2.*w(i+1,j,k,4) -w(i+2,j,k,4)
         w(i,j,k,5)  = 2.*w(i+1,j,k,5) -w(i+2,j,k,5)
         p(i,j,k)    = 2.*p(i+1,j,k)   -p(i+2,j,k)
         eomu(i,j,k) = evinf
  120    continue
      endif
      return
c
c
  200 continue
c--------------------- i = imx  face ----------------------------
c
      if( n1beg .ne. jmn) n1beg = n1beg +1
      if( n1end .eq. jmx) n1end = n1end +1
      if( n2beg .ne. kmn) n2beg = n2beg +1
      if( n2end .eq. kmx) n2end = n2end +1
c
      i   =   imx + 1
c
      if(rmjet.lt.1.0) then
      do 210 k= n2beg, n2end
cdir$ ivdep
      do 210 j= n1beg, n1end
      uf        = w(i-1,j,k,2)/w(i-1,j,k,1) 
      vf        = w(i-1,j,k,3)/w(i-1,j,k,1) 
      wf        = w(i-1,j,k,4)/w(i-1,j,k,1) 
      qnf       = uf*ri2(j,k,1)  +vf*ri2(j,k,2)  +wf*ri2(j,k,3) 
      qn0       = ujet*ri2(j,k,1)+vjet*ri2(j,k,2)+wjet*ri2(j,k,3) 
      c         = sqrt(gamma*p(i-1,j,k)/w(i-1,j,k,1)) 
      er        = qnf + 2.*dgm1*c
      fr        = qn0 - 2.*dgm1*cjet
      c         = .25*(er - fr)*gm1
      qn        = .50*(er + fr)
c 
c--------- inflow boundary
c 
      if (qn.le.0.) then
        uf        = ujet + (qn - qn0)*ri2(j,k,1)
        vf        = vjet + (qn - qn0)*ri2(j,k,2)
        wf        = wjet + (qn - qn0)*ri2(j,k,3)
        s         = sjet
c 
c--------- outflow boundary 
c 
      else
        uf        = uf + (qn - qnf)*ri2(j,k,1)
        vf        = vf + (qn - qnf)*ri2(j,k,2)
        wf        = wf + (qn - qnf)*ri2(j,k,3)
        s         = w(i-1,j,k,1)**gamma/p(i-1,j,k)
      endif 
c 
      cc          = c*c/gamma
ccc      w(i,j,k,1)  = (s*cc)**dgm1
         vtots       =  uf**2 + vf**2 + wf**2
         rmloc       = sqrt(vtots/c**2)
         p(i,j,k)    = ptotjet/((1.+.5*gm1*rmloc*rmloc)**gogm1)
ccc      tloc        = t0*ttjotfs/(1.+.5*gm1*rmloc*rmloc)
ccc      w(i,j,k,1)  = p(i,j,k)/tloc
         w(i,j,k,1)  = p(i,j,k)/cc
         w(i,j,k,2)  = w(i,j,k,1)*uf
         w(i,j,k,3)  = w(i,j,k,1)*vf
         w(i,j,k,4)  = w(i,j,k,1)*wf
ccc      p(i,j,k)    = w(i,j,k,1)*cc
c
c        energy is overspecified
c     w(i,j,k,5) = -p(i,j,k) + h0*w(i,j,k,1)
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  210 continue
      else
c        supersonic flow condition
         do 215 k= n2beg, n2end
         do 215 j= n1beg, n1end
         w(i,j,k,1)  = rhojet
         w(i,j,k,2)  = rhojet*ujet
         w(i,j,k,3)  = rhojet*vjet
         w(i,j,k,4)  = rhojet*wjet
         w(i,j,k,5)  = rhojet*hjet -pjet
         p(i,j,k)    = pjet
         eomu(i,j,k) = evinf
  215    continue
      end if
c
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         i   =  imx + 2
c
         do 220 k= n2beg, n2end
         do 220 j= n1beg, n1end
         w(i,j,k,1)  = 2.*w(i-1,j,k,1) -w(i-2,j,k,1)
         w(i,j,k,2)  = 2.*w(i-1,j,k,2) -w(i-2,j,k,2)
         w(i,j,k,3)  = 2.*w(i-1,j,k,3) -w(i-2,j,k,3)
         w(i,j,k,4)  = 2.*w(i-1,j,k,4) -w(i-2,j,k,4)
         w(i,j,k,5)  = 2.*w(i-1,j,k,5) -w(i-2,j,k,5)
         p(i,j,k)    = 2.*p(i-1,j,k)   -p(i-2,j,k)
         eomu(i,j,k) = evinf
  220    continue
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
      if(rmjet.lt.1.0) then
      do 310 k= n1beg, n1end
cdir$ ivdep
      do 310 i= n2beg, n2end
      uf        = w(i,j+1,k,2)/w(i,j+1,k,1) 
      vf        = w(i,j+1,k,3)/w(i,j+1,k,1) 
      wf        = w(i,j+1,k,4)/w(i,j+1,k,1) 
      qnf       = uf*rj1(i,k,1)  +vf*rj1(i,k,2)  +wf*rj1(i,k,3) 
      qn0       = ujet*rj1(i,k,1)+vjet*rj1(i,k,2)+wjet*rj1(i,k,3) 
      c         = sqrt(gamma*p(i,j+1,k)/w(i,j+1,k,1)) 
      er        = qnf - 2.*dgm1*c
      fr        = qn0 + 2.*dgm1*cjet
      c         = .25*(er - fr)*gm1
      qn        = .50*(er + fr)
c 
c--------- inflow boundary
c 
      if (qn.gt.0.) then
        uf        = ujet + (qn - qn0)*rj1(i,k,1)
        vf        = vjet + (qn - qn0)*rj1(i,k,2)
        wf        = wjet + (qn - qn0)*rj1(i,k,3)
        s         = sjet
c 
c--------- outflow boundary 
c 
      else
        uf        = uf + (qn - qnf)*rj1(i,k,1)
        vf        = vf + (qn - qnf)*rj1(i,k,2)
        wf        = wf + (qn - qnf)*rj1(i,k,3)
        s         = w(i,j+1,k,1)**gamma/p(i,j+1,k)
      endif 
c 
      cc          = c**2/gamma
ccc      w(i,j,k,1)  = (s*cc)**dgm1
         vtots       =  uf**2 + vf**2 + wf**2
         rmloc       = sqrt(vtots/c**2)
         p(i,j,k)    = ptotjet/((1.+.5*gm1*rmloc*rmloc)**gogm1)
ccc      tloc        = t0*ttjotfs/(1.+.5*gm1*rmloc*rmloc)
ccc      w(i,j,k,1)  = p(i,j,k)/tloc
         w(i,j,k,1)  = p(i,j,k)/cc
         w(i,j,k,2)  = w(i,j,k,1)*uf
         w(i,j,k,3)  = w(i,j,k,1)*vf
         w(i,j,k,4)  = w(i,j,k,1)*wf
ccc      p(i,j,k)    = w(i,j,k,1)*cc
c
c        energy is overspecified
c     w(i,j,k,5) = -p(i,j,k) + h0*w(i,j,k,1)
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  310 continue
c
      else
c        supersonic flow condition
         do 315 k= n1beg, n1end
         do 315 i= n2beg, n2end
         w(i,j,k,1)  = rhojet
         w(i,j,k,2)  = rhojet*ujet
         w(i,j,k,3)  = rhojet*vjet
         w(i,j,k,4)  = rhojet*wjet
         w(i,j,k,5)  = rhojet*hjet -pjet
         p(i,j,k)    = pjet
         eomu(i,j,k) = evinf
  315    continue
      end if
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         j   =  jmn - 1
c
         do 320 k= n1beg, n1end
         do 320 i= n2beg, n2end
         w(i,j,k,1)  = 2.*w(i,j+1,k,1) -w(i,j+2,k,1)
         w(i,j,k,2)  = 2.*w(i,j+1,k,2) -w(i,j+2,k,2)
         w(i,j,k,3)  = 2.*w(i,j+1,k,3) -w(i,j+2,k,3)
         w(i,j,k,4)  = 2.*w(i,j+1,k,4) -w(i,j+2,k,4)
         w(i,j,k,5)  = 2.*w(i,j+1,k,5) -w(i,j+2,k,5)
         p(i,j,k)    = 2.*p(i,j+1,k)   -p(i,j+2,k)
         eomu(i,j,k) = evinf
  320    continue
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
      if(rmjet.lt.1.0) then
      do 410 k= n1beg, n1end
cdir$ ivdep
      do 410 i= n2beg, n2end
      uf        = w(i,j-1,k,2)/w(i,j-1,k,1) 
      vf        = w(i,j-1,k,3)/w(i,j-1,k,1) 
      wf        = w(i,j-1,k,4)/w(i,j-1,k,1) 
      qnf       = uf*rj2(i,k,1)  +vf*rj2(i,k,2)  +wf*rj2(i,k,3) 
      qn0       = u0*rj2(i,k,1)  +v0*rj2(i,k,2)  +w0*rj2(i,k,3) 
      c         = sqrt(gamma*p(i,j-1,k)/w(i,j-1,k,1)) 
      er        = qnf + 2.*dgm1*c
      fr        = qn0 - 2.*dgm1*cjet
      c         = .25*(er - fr)*gm1
      qn        = .50*(er + fr)
c 
c--------- inflow boundary
c 
      if (qn.le.0.) then
        uf        = ujet + (qn - qn0)*rj2(i,k,1)
        vf        = vjet + (qn - qn0)*rj2(i,k,2)
        wf        = wjet + (qn - qn0)*rj2(i,k,3)
        s         = sjet
c 
c--------- outflow boundary 
c 
      else
        uf        = uf + (qn - qnf)*rj2(i,k,1)
        vf        = vf + (qn - qnf)*rj2(i,k,2)
        wf        = wf + (qn - qnf)*rj2(i,k,3)
        s         = w(i,j-1,k,1)**gamma/p(i,j-1,k)
      endif 
c 
      cc          = c*c/gamma
ccc      w(i,j,k,1)  = (s*cc)**dgm1
         vtots       =  uf**2 + vf**2 + wf**2
         rmloc       = sqrt(vtots/c**2)
         p(i,j,k)    = ptotjet/((1.+.5*gm1*rmloc*rmloc)**gogm1)
ccc      tloc        = t0*ttjotfs/(1.+.5*gm1*rmloc*rmloc)
ccc      w(i,j,k,1)  = p(i,j,k)/tloc
         w(i,j,k,1)  = p(i,j,k)/cc
         w(i,j,k,2)  = w(i,j,k,1)*uf
         w(i,j,k,3)  = w(i,j,k,1)*vf
         w(i,j,k,4)  = w(i,j,k,1)*wf
ccc      p(i,j,k)    = w(i,j,k,1)*cc
c
c        energy is overspecified
c     w(i,j,k,5) = -p(i,j,k) + h0*w(i,j,k,1)
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  410 continue
c
      else
c        supersonic flow condition
         do 415 k= n1beg, n1end
         do 415 i= n2beg, n2end
         w(i,j,k,1)  = rhojet
         w(i,j,k,2)  = rhojet*ujet
         w(i,j,k,3)  = rhojet*vjet
         w(i,j,k,4)  = rhojet*wjet
         w(i,j,k,5)  = rhojet*hjet -pjet
         p(i,j,k)    = pjet
         eomu(i,j,k) = evinf
  415    continue
      end if
c
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         j   =  jmx + 2
c
         do 420 k= n1beg, n1end
         do 420 i= n2beg, n2end
         w(i,j,k,1)  = 2.*w(i,j-1,k,1) -w(i,j-2,k,1)
         w(i,j,k,2)  = 2.*w(i,j-1,k,2) -w(i,j-2,k,2)
         w(i,j,k,3)  = 2.*w(i,j-1,k,3) -w(i,j-2,k,3)
         w(i,j,k,4)  = 2.*w(i,j-1,k,4) -w(i,j-2,k,4)
         w(i,j,k,5)  = 2.*w(i,j-1,k,5) -w(i,j-2,k,5)
         p(i,j,k)    = 2.*p(i,j-1,k)   -p(i,j-2,k)
         eomu(i,j,k) = evinf
  420    continue
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
      if(rmjet.lt.1.0) then
      do 510 j= n2beg, n2end
cdir$ ivdep
      do 510 i= n1beg, n1end
      uf        = w(i,j,k+1,2)/w(i,j,k+1,1) 
      vf        = w(i,j,k+1,3)/w(i,j,k+1,1) 
      wf        = w(i,j,k+1,4)/w(i,j,k+1,1) 
      qnf       = uf*rk1(i,j,1)  +vf*rk1(i,j,2)  +wf*rk1(i,j,3) 
      qn0       = ujet*rk1(i,j,1)+vjet*rk1(i,j,2)+wjet*rk1(i,j,3) 
      c         = sqrt(gamma*p(i,j,k+1)/w(i,j,k+1,1)) 
      er        = qnf - 2.*dgm1*c
      fr        = qn0 + 2.*dgm1*cjet
      c         = .25*(er - fr)*gm1
      qn        = .50*(er + fr)
c 
c--------- inflow boundary
c 
      if (qn.gt.0.) then
        uf        = ujet + (qn - qn0)*rk1(i,j,1)
        vf        = vjet + (qn - qn0)*rk1(i,j,2)
        wf        = wjet + (qn - qn0)*rk1(i,j,3)
        s         = sjet
c 
c--------- outflow boundary 
c 
      else
        uf        = uf + (qn - qnf)*rk1(i,j,1)
        vf        = vf + (qn - qnf)*rk1(i,j,2)
        wf        = wf + (qn - qnf)*rk1(i,j,3)
        s         = w(i,j,k+1,1)**gamma/p(i,j,k+1)
      endif 
c 
      cc          = c**2/gamma
ccc      w(i,j,k,1)  = (s*cc)**dgm1
         vtots       =  uf**2 + vf**2 + wf**2
         rmloc       = sqrt(vtots/c**2)
         p(i,j,k)    = ptotjet/((1.+.5*gm1*rmloc*rmloc)**gogm1)
ccc      tloc        = t0*ttjotfs/(1.+.5*gm1*rmloc*rmloc)
ccc      w(i,j,k,1)  = p(i,j,k)/tloc
         w(i,j,k,1)  = p(i,j,k)/cc
         w(i,j,k,2)  = w(i,j,k,1)*uf
         w(i,j,k,3)  = w(i,j,k,1)*vf
         w(i,j,k,4)  = w(i,j,k,1)*wf
ccc      p(i,j,k)    = w(i,j,k,1)*cc
c
c        energy is overspecified
c     w(i,j,k,5) = -p(i,j,k) + h0*w(i,j,k,1)
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  510 continue
c
      else
c        supersonic flow condition
         do 515 j= n2beg, n2end
         do 515 i= n1beg, n1end
         w(i,j,k,1)  = rhojet
         w(i,j,k,2)  = rhojet*ujet
         w(i,j,k,3)  = rhojet*vjet
         w(i,j,k,4)  = rhojet*wjet
         w(i,j,k,5)  = rhojet*hjet -pjet
         p(i,j,k)    = pjet
         eomu(i,j,k) = evinf
  515    continue
      end if
c
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         k   =  kmn - 1
c
         do 520 j= n2beg, n2end
         do 520 i= n1beg, n1end
         w(i,j,k,1)  = 2.*w(i,j,k+1,1) -w(i,j,k+2,1)
         w(i,j,k,2)  = 2.*w(i,j,k+1,2) -w(i,j,k+2,2)
         w(i,j,k,3)  = 2.*w(i,j,k+1,3) -w(i,j,k+2,3)
         w(i,j,k,4)  = 2.*w(i,j,k+1,4) -w(i,j,k+2,4)
         w(i,j,k,5)  = 2.*w(i,j,k+1,5) -w(i,j,k+2,5)
         p(i,j,k)    = 2.*p(i,j,k+1)   -p(i,j,k+2)
         eomu(i,j,k) = evinf
  520    continue
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
      if(rmjet.lt.1.0) then
      do 610 j= n2beg, n2end
cdir$ ivdep
      do 610 i= n1beg, n1end
      uf        = w(i,j,k-1,2)/w(i,j,k-1,1) 
      vf        = w(i,j,k-1,3)/w(i,j,k-1,1) 
      wf        = w(i,j,k-1,4)/w(i,j,k-1,1) 
      qnf       = uf*rk2(i,j,1)  +vf*rk2(i,j,2)  +wf*rk2(i,j,3) 
      qn0       = u0*rk2(i,j,1)  +v0*rk2(i,j,2)  +w0*rk2(i,j,3) 
      c         = sqrt(gamma*p(i,j,k-1)/w(i,j,k-1,1)) 
      er        = qnf + 2.*dgm1*c
      fr        = qn0 - 2.*dgm1*cjet
      c         = .25*(er - fr)*gm1
      qn        = .50*(er + fr)
c 
c--------- inflow boundary
c 
      if (qn.le.0.) then
        uf        = ujet + (qn - qn0)*rk2(i,j,1)
        vf        = vjet + (qn - qn0)*rk2(i,j,2)
        wf        = wjet + (qn - qn0)*rk2(i,j,3)
        s         = sjet
c 
c--------- outflow boundary 
c 
      else
        uf        = uf + (qn - qnf)*rk2(i,j,1)
        vf        = vf + (qn - qnf)*rk2(i,j,2)
        wf        = wf + (qn - qnf)*rk2(i,j,3)
        s         = w(i,j,k-1,1)**gamma/p(i,j,k-1)
      endif 
c 
      cc          = c**2/gamma
ccc      w(i,j,k,1)  = (s*cc)**dgm1
         vtots       =  uf**2 + vf**2 + wf**2
         rmloc       = sqrt(vtots/c**2)
         p(i,j,k)    = ptotjet/((1.+.5*gm1*rmloc*rmloc)**gogm1)
ccc      tloc        = t0*ttjotfs/(1.+.5*gm1*rmloc*rmloc)
ccc      w(i,j,k,1)  = p(i,j,k)/tloc
         w(i,j,k,1)  = p(i,j,k)/cc
         w(i,j,k,2)  = w(i,j,k,1)*uf
         w(i,j,k,3)  = w(i,j,k,1)*vf
         w(i,j,k,4)  = w(i,j,k,1)*wf
ccc      p(i,j,k)    = w(i,j,k,1)*cc
c
c        energy is overspecified
c     w(i,j,k,5) = -p(i,j,k) + h0*w(i,j,k,1)
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  610 continue
c
      else
c        supersonic flow condition
         do 615 j= n2beg, n2end
         do 615 i= n1beg, n1end
         w(i,j,k,1)  = rhojet
         w(i,j,k,2)  = rhojet*ujet
         w(i,j,k,3)  = rhojet*vjet
         w(i,j,k,4)  = rhojet*wjet
         w(i,j,k,5)  = rhojet*hjet -pjet
         p(i,j,k)    = pjet
         eomu(i,j,k) = evinf
  615    continue
      end if
c
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         k   =  kmx + 2
c
         do 620 j= n2beg, n2end
         do 620 i= n1beg, n1end
         w(i,j,k,1)  = 2.*w(i,j,k-1,1) -w(i,j,k-2,1)
         w(i,j,k,2)  = 2.*w(i,j,k-1,2) -w(i,j,k-2,2)
         w(i,j,k,3)  = 2.*w(i,j,k-1,3) -w(i,j,k-2,3)
         w(i,j,k,4)  = 2.*w(i,j,k-1,4) -w(i,j,k-2,4)
         w(i,j,k,5)  = 2.*w(i,j,k-1,5) -w(i,j,k-2,5)
         p(i,j,k)    = 2.*p(i,j,k-1)   -p(i,j,k-2)
         eomu(i,j,k) = evinf
  620    continue
      endif
c
c
      return
      end
