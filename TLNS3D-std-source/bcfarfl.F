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
      subroutine  bcfarfl(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                    imp2,jmp2,kmp2,w,p,eomu,turv1,
     .                    ri1,ri2,rj1,rj2,rk1,rk2,
     .                    nface,n1beg,n1end,n2beg,n2end,
     .                    igrid,isoln,extyp,iunsteady                )
c=======================================================================
c 
c     object : to compute the ghost cell values for flow variables
c              at far-field boundaries
c              presently Riemann invariants based on 1-d
c              characteristics are used
c
c     initial coding : by v.n.vatsa, nasa langley ( March, 1990)
c     latest mods    : to use appropriate characteristic variables
c                      for preconditioning (Feb. 1998)
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
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubcon2,iprec
c
      common/trbfs/ anutinf,zkinf,ominf
      common/unsprec/ npcons,nexp,iprecg,ijac,icord,inav
c
c
      if (iunsteady.eq.0) then
         iiprec = iprec
      else
c        turn low speed preconditiong off for unsteady flow
         iiprec = 0
      endif
c
c------------------------------------------------------------------
c
      s0    =    rho0**gamma/p0
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
c----- mods by vatsa (10-3-94) to compute only interior cells
c
cvnv   n1beg = n1beg + 1
cvnv   n2beg = n2beg + 1
c-----
c      branch-cut for different type of boundary faces
c
      go to (100,200,300,400,500,600) nface
c
c
  100 continue
c--------------------- i = imn  face ----------------------------
c
c----- mods by vatsa (10-3-94) to compute only interior cells
      if( n1beg .ne. jmn) n1beg = n1beg +1
      if( n1end .eq. jmx) n1end = n1end +1
      if( n2beg .ne. kmn) n2beg = n2beg +1
      if( n2end .eq. kmx) n2end = n2end +1
c
      i    =    imn
c
c
      if(rm.lt.1.0) then
         do 110 k= n2beg, n2end
         if(iiprec.le.0) then
c-----   subsonic free stream (no preconditioning)
cdir$ ivdep
         do 105 j= n1beg, n1end
         uf        = w(i+1,j,k,2)/w(i+1,j,k,1) 
         vf        = w(i+1,j,k,3)/w(i+1,j,k,1) 
         wf        = w(i+1,j,k,4)/w(i+1,j,k,1) 
         qnf       = uf*ri1(j,k,1)  +vf*ri1(j,k,2)  +wf*ri1(j,k,3)
         qn0       = u0*ri1(j,k,1)  +v0*ri1(j,k,2)  +w0*ri1(j,k,3)
         c         = sqrt(gamma*p(i+1,j,k)/w(i+1,j,k,1)) 
         er        = qnf - 2.*dgm1*c
         fr        = qn0 + 2.*dgm1*c0
         c         = .25*(er - fr)*gm1
         qn        = .50*(er + fr)
c
c--------- inflow boundary (subsonic)
c 
         if (qn.gt.0.) then
           uf        = u0 + (qn - qn0)*ri1(j,k,1)
           vf        = v0 + (qn - qn0)*ri1(j,k,2)
           wf        = w0 + (qn - qn0)*ri1(j,k,3)
           s         = s0
c 
c--------- outflow boundary (subsonic)
c 
         else
           uf        = uf + (qn - qnf)*ri1(j,k,1)
           vf        = vf + (qn - qnf)*ri1(j,k,2)
           wf        = wf + (qn - qnf)*ri1(j,k,3)
           s         = w(i+1,j,k,1)**gamma/p(i+1,j,k)
         endif 
c 
         cc          = c*c/gamma
         w(i,j,k,1)  = (s*cc)**dgm1
         w(i,j,k,2)  = w(i,j,k,1)*uf
         w(i,j,k,3)  = w(i,j,k,1)*vf
         w(i,j,k,4)  = w(i,j,k,1)*wf
         p(i,j,k)    = w(i,j,k,1)*cc
c
c        specify energy 
c     w(i,j,k,5)  = -p(i,j,k) + h0*w(i,j,k,1)
c
         vtots       =     w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                     w(i,j,k,4)**2
         w(i,j,k,5)  =     p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
         eomu(i,j,k) =     evinf
  105    continue
c
      else
c
c-----   subsonic free stream (with preconditioning)
c
cdir$ ivdep
      do 108 j= n1beg, n1end
      uf        = w(i+1,j,k,2)/w(i+1,j,k,1)
      vf        = w(i+1,j,k,3)/w(i+1,j,k,1)
      wf        = w(i+1,j,k,4)/w(i+1,j,k,1)
      qnf       = uf*ri1(j,k,1) + vf*ri1(j,k,2) + wf*ri1(j,k,3)
c
c--------- inflow boundary (subsonic)
c
      if (qnf.gt.0.0) then
        uf        = u0
        vf        = v0
        wf        = w0
        Tf        = p0/rho0
        p(i,j,k)  = p(i+1,j,k)
        w(i,j,k,1)  = p(i,j,k)/Tf
c
c--------- outflow boundary (subsonic)
c
      else
        uf        = w(i+1,j,k,2)/w(i+1,j,k,1)
        vf        = w(i+1,j,k,3)/w(i+1,j,k,1)
        wf        = w(i+1,j,k,4)/w(i+1,j,k,1)
        Tf        = p(i+1,j,k)  /w(i+1,j,k,1)
        p(i,j,k)  = p0
        w(i,j,k,1)  = p(i,j,k)/Tf
      endif
c
      w(i,j,k,2)  = w(i,j,k,1)*uf
      w(i,j,k,3)  = w(i,j,k,1)*vf
      w(i,j,k,4)  = w(i,j,k,1)*wf
c
c       compute energy
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k)*dgm1 + .5*vtots/w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  108 continue
c     end of subsonic inflow/outflow condition
      endif
  110 continue
      else
c
c-----   supersonic free stream
c
         do 115 k= n2beg, n2end
         do 115 j= n1beg, n1end
         uf        = w(i+1,j,k,2)/w(i+1,j,k,1) 
         vf        = w(i+1,j,k,3)/w(i+1,j,k,1) 
         wf        = w(i+1,j,k,4)/w(i+1,j,k,1) 
         qnf       = uf*ri1(j,k,1)  +vf*ri1(j,k,2)  +wf*ri1(j,k,3)
         qn0       = u0*ri1(j,k,1)  +v0*ri1(j,k,2)  +w0*ri1(j,k,3)
         c         = sqrt(gamma*p(i+1,j,k)/w(i+1,j,k,1)) 
         er        = qnf - 2.*dgm1*c
         fr        = qn0 + 2.*dgm1*c0
         c         = .25*(er - fr)*gm1
         qn        = .50*(er + fr)
c 
c--------- inflow boundary (supersonic)
c 
       if (qn.gt.0.) then
         w(i,j,k,1)  = rho0
         w(i,j,k,2)  = rho0*u0
         w(i,j,k,3)  = rho0*v0
         w(i,j,k,4)  = rho0*w0
         w(i,j,k,5)  = rho0*h0 -p0
         p(i,j,k)    = p0
         eomu(i,j,k) = evinf
       else
c
c------- outflow boundary (supersonic)
c 
         w(i,j,k,1)  = w(i+1,j,k,1)
         w(i,j,k,2)  = w(i+1,j,k,2)
         w(i,j,k,3)  = w(i+1,j,k,3)
         w(i,j,k,4)  = w(i+1,j,k,4)
         w(i,j,k,5)  = w(i+1,j,k,5)
         p(i,j,k)    = p(i+1,j,k)
         eomu(i,j,k) = eomu(i+1,j,k)
c
       endif
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
c----- mods by vatsa (10-3-94) to compute only interior cells
      if( n1beg .ne. jmn) n1beg = n1beg +1
      if( n1end .eq. jmx) n1end = n1end +1
      if( n2beg .ne. kmn) n2beg = n2beg +1
      if( n2end .eq. kmx) n2end = n2end +1
c
      i   =   imx + 1
c
      if(rm.lt.1.0) then
c------ subsonic free stream
      do 210 k= n2beg, n2end
      if(iiprec.le.0) then
c-----   subsonic free stream (no preconditioning)
cdir$ ivdep
      do 205 j= n1beg, n1end
      uf        = w(i-1,j,k,2)/w(i-1,j,k,1) 
      vf        = w(i-1,j,k,3)/w(i-1,j,k,1) 
      wf        = w(i-1,j,k,4)/w(i-1,j,k,1) 
      qnf       = uf*ri2(j,k,1)  +vf*ri2(j,k,2)  +wf*ri2(j,k,3) 
      qn0       = u0*ri2(j,k,1)  +v0*ri2(j,k,2)  +w0*ri2(j,k,3) 
      c         = sqrt(gamma*p(i-1,j,k)/w(i-1,j,k,1)) 
      er        = qnf + 2.*dgm1*c
      fr        = qn0 - 2.*dgm1*c0
      c         = .25*(er - fr)*gm1
      qn        = .50*(er + fr)
c 
c--------- inflow boundary (subsonic)
c 
      if (qn.le.0.) then
        uf        = u0 + (qn - qn0)*ri2(j,k,1)
        vf        = v0 + (qn - qn0)*ri2(j,k,2)
        wf        = w0 + (qn - qn0)*ri2(j,k,3)
        s         = s0
c 
c--------- outflow boundary (subsonic)
c 
      else
        uf        = uf + (qn - qnf)*ri2(j,k,1)
        vf        = vf + (qn - qnf)*ri2(j,k,2)
        wf        = wf + (qn - qnf)*ri2(j,k,3)
        s         = w(i-1,j,k,1)**gamma/p(i-1,j,k)
      endif 
c 
      cc          = c*c/gamma
      w(i,j,k,1)  = (s*cc)**dgm1
      w(i,j,k,2)  = w(i,j,k,1)*uf
      w(i,j,k,3)  = w(i,j,k,1)*vf
      w(i,j,k,4)  = w(i,j,k,1)*wf
      p(i,j,k)    = w(i,j,k,1)*cc
c
c        specify energy 
c     w(i,j,k,5) = -p(i,j,k) + h0*w(i,j,k,1)
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  205 continue
c
      else
c
c-----   subsonic free stream (with preconditioning)
c
cdir$ ivdep
      do 208 j= n1beg, n1end
      uf        = w(i-1,j,k,2)/w(i-1,j,k,1)
      vf        = w(i-1,j,k,3)/w(i-1,j,k,1)
      wf        = w(i-1,j,k,4)/w(i-1,j,k,1)
      qnf       = uf*ri2(j,k,1) + vf*ri2(j,k,2) + wf*ri2(j,k,3)
c
c--------- inflow boundary (subsonic)
c
      if (qnf.le.0.0) then
        uf        = u0
        vf        = v0
        wf        = w0
        Tf        = p0/rho0
        p(i,j,k)  = p(i-1,j,k)
        w(i,j,k,1)  = p(i,j,k)/Tf
c
c--------- outflow boundary (subsonic)
c
      else
        uf        = w(i-1,j,k,2)/w(i-1,j,k,1)
        vf        = w(i-1,j,k,3)/w(i-1,j,k,1)
        wf        = w(i-1,j,k,4)/w(i-1,j,k,1)
        Tf        = p(i-1,j,k)  /w(i-1,j,k,1)
        p(i,j,k)  = p0
        w(i,j,k,1)  = p(i,j,k)/Tf
      endif
c
      w(i,j,k,2)  = w(i,j,k,1)*uf
      w(i,j,k,3)  = w(i,j,k,1)*vf
      w(i,j,k,4)  = w(i,j,k,1)*wf
c
c       compute energy
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k)*dgm1 + .5*vtots/w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  208 continue
c     end of subsonic inflow/outflow condition
      endif
  210 continue
      else
c        supersonic flow condition
         do 215 k= n2beg, n2end
         do 215 j= n1beg, n1end
         uf        = w(i-1,j,k,2)/w(i-1,j,k,1) 
         vf        = w(i-1,j,k,3)/w(i-1,j,k,1) 
         wf        = w(i-1,j,k,4)/w(i-1,j,k,1) 
         qnf       = uf*ri2(j,k,1)  +vf*ri2(j,k,2)  +wf*ri2(j,k,3) 
         qn0       = u0*ri2(j,k,1)  +v0*ri2(j,k,2)  +w0*ri2(j,k,3) 
         c         = sqrt(gamma*p(i-1,j,k)/w(i-1,j,k,1)) 
         er        = qnf + 2.*dgm1*c
         fr        = qn0 - 2.*dgm1*c0
         c         = .25*(er - fr)*gm1
         qn        = .50*(er + fr)
c 
c--------- inflow boundary (supersonic)
c 
       if (qn.gt.0.) then
         w(i,j,k,1)  = rho0
         w(i,j,k,2)  = rho0*u0
         w(i,j,k,3)  = rho0*v0
         w(i,j,k,4)  = rho0*w0
         w(i,j,k,5)  = rho0*h0 -p0
         p(i,j,k)    = p0
         eomu(i,j,k) = evinf
       else
c
c------- outflow boundary (supersonic)
c 
         w(i,j,k,1)  = w(i-1,j,k,1)
         w(i,j,k,2)  = w(i-1,j,k,2)
         w(i,j,k,3)  = w(i-1,j,k,3)
         w(i,j,k,4)  = w(i-1,j,k,4)
         w(i,j,k,5)  = w(i-1,j,k,5)
         p(i,j,k)    = p(i-1,j,k)
         eomu(i,j,k) = eomu(i-1,j,k)
c
       endif
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
c----- mods by vatsa (10-3-94) to compute only interior cells
      if( n1beg .ne. kmn) n1beg = n1beg +1
      if( n1end .eq. kmx) n1end = n1end +1
      if( n2beg .ne. imn) n2beg = n2beg +1
      if( n2end .eq. imx) n2end = n2end +1
c
      j   =   jmn
c
      if(rm.lt.1.0) then
      do 310 k= n1beg, n1end
      if(iiprec.le.0) then
c-----   subsonic free stream (no preconditioning)
cdir$ ivdep
      do 305 i= n2beg, n2end
      uf        = w(i,j+1,k,2)/w(i,j+1,k,1) 
      vf        = w(i,j+1,k,3)/w(i,j+1,k,1) 
      wf        = w(i,j+1,k,4)/w(i,j+1,k,1) 
      qnf       = uf*rj1(i,k,1)  +vf*rj1(i,k,2)  +wf*rj1(i,k,3) 
      qn0       = u0*rj1(i,k,1)  +v0*rj1(i,k,2)  +w0*rj1(i,k,3) 
      c         = sqrt(gamma*p(i,j+1,k)/w(i,j+1,k,1)) 
      er        = qnf - 2.*dgm1*c
      fr        = qn0 + 2.*dgm1*c0
      c         = .25*(er - fr)*gm1
      qn        = .50*(er + fr)
c 
c--------- inflow boundary (subsonic)
c 
      if (qn.gt.0.) then
        uf        = u0 + (qn - qn0)*rj1(i,k,1)
        vf        = v0 + (qn - qn0)*rj1(i,k,2)
        wf        = w0 + (qn - qn0)*rj1(i,k,3)
        s         = s0
c 
c--------- outflow boundary (subsonic) 
c 
      else
        uf        = uf + (qn - qnf)*rj1(i,k,1)
        vf        = vf + (qn - qnf)*rj1(i,k,2)
        wf        = wf + (qn - qnf)*rj1(i,k,3)
        s         = w(i,j+1,k,1)**gamma/p(i,j+1,k)
      endif 
c 
      cc          = c**2/gamma
      w(i,j,k,1)  = (s*cc)**dgm1
      w(i,j,k,2)  = w(i,j,k,1)*uf
      w(i,j,k,3)  = w(i,j,k,1)*vf
      w(i,j,k,4)  = w(i,j,k,1)*wf
      p(i,j,k)    = w(i,j,k,1)*cc
c
c        compute energy
c     w(i,j,k,5) = -p(i,j,k) + h0*w(i,j,k,1)
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  305 continue
c
      else
c
c-----   subsonic free stream (with preconditioning)
c
cdir$ ivdep
      do 308 i= n2beg, n2end
      uf        = w(i,j+1,k,2)/w(i,j+1,k,1)
      vf        = w(i,j+1,k,3)/w(i,j+1,k,1)
      wf        = w(i,j+1,k,4)/w(i,j+1,k,1)
      qnf       = uf*rj1(i,k,1) + vf*rj1(i,k,2) + wf*rj1(i,k,3)
c
c--------- inflow boundary (subsonic)
c
      if (qnf.gt.0.0) then
        uf        = u0
        vf        = v0
        wf        = w0
        Tf        = p0/rho0
        p(i,j,k)  = p(i,j+1,k)
        w(i,j,k,1)  = p(i,j,k)/Tf
c
c--------- outflow boundary (subsonic)
c
      else
        uf        = w(i,j+1,k,2)/w(i,j+1,k,1)
        vf        = w(i,j+1,k,3)/w(i,j+1,k,1)
        wf        = w(i,j+1,k,4)/w(i,j+1,k,1)
        Tf        = p(i,j+1,k)  /w(i,j+1,k,1)
        p(i,j,k)  = p0
        w(i,j,k,1)  = p(i,j,k)/Tf
      endif
c
      w(i,j,k,2)  = w(i,j,k,1)*uf
      w(i,j,k,3)  = w(i,j,k,1)*vf
      w(i,j,k,4)  = w(i,j,k,1)*wf
c
c       compute energy
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k)*dgm1 + .5*vtots/w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  308 continue
c     end of subsonic inflow/outflow condition
      endif
  310 continue
      else
c        supersonic flow condition
         do 315 k= n1beg, n1end
         do 315 i= n2beg, n2end
         uf        = w(i,j+1,k,2)/w(i,j+1,k,1) 
         vf        = w(i,j+1,k,3)/w(i,j+1,k,1) 
         wf        = w(i,j+1,k,4)/w(i,j+1,k,1) 
         qnf       = uf*rj1(i,k,1)  +vf*rj1(i,k,2)  +wf*rj1(i,k,3) 
         qn0       = u0*rj1(i,k,1)  +v0*rj1(i,k,2)  +w0*rj1(i,k,3) 
         c         = sqrt(gamma*p(i,j+1,k)/w(i,j+1,k,1)) 
         er        = qnf - 2.*dgm1*c
         fr        = qn0 + 2.*dgm1*c0
         c         = .25*(er - fr)*gm1
         qn        = .50*(er + fr)
c 
c--------- inflow boundary (supersonic)
c 
      if (qn.gt.0.) then
         w(i,j,k,1)  = rho0
         w(i,j,k,2)  = rho0*u0
         w(i,j,k,3)  = rho0*v0
         w(i,j,k,4)  = rho0*w0
         w(i,j,k,5)  = rho0*h0 -p0
         p(i,j,k)    = p0
         eomu(i,j,k) = evinf
       else
c
c------- outflow boundary (supersonic)
c 
         w(i,j,k,1)  = w(i,j+1,k,1)
         w(i,j,k,2)  = w(i,j+1,k,2)
         w(i,j,k,3)  = w(i,j+1,k,3)
         w(i,j,k,4)  = w(i,j+1,k,4)
         w(i,j,k,5)  = w(i,j+1,k,5)
         p(i,j,k)    = p(i,j+1,k)
         eomu(i,j,k) = eomu(i,j+1,k)
c
       endif
  315    continue
      end if
c
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
c----- mods by vatsa (10-3-94) to compute only interior cells
      if( n1beg .ne. kmn) n1beg = n1beg +1
      if( n1end .eq. kmx) n1end = n1end +1
      if( n2beg .ne. imn) n2beg = n2beg +1
      if( n2end .eq. imx) n2end = n2end +1
c
      j   =   jmx + 1
c
      if(rm.lt.1.0) then
c-----   subsonic free stream (no preconditioning)
      do 410 k= n1beg, n1end
      if(iiprec.le.0) then
cdir$ ivdep
      do 405 i= n2beg, n2end
      uf        = w(i,j-1,k,2)/w(i,j-1,k,1) 
      vf        = w(i,j-1,k,3)/w(i,j-1,k,1) 
      wf        = w(i,j-1,k,4)/w(i,j-1,k,1) 
      qnf       = uf*rj2(i,k,1)  +vf*rj2(i,k,2)  +wf*rj2(i,k,3) 
      qn0       = u0*rj2(i,k,1)  +v0*rj2(i,k,2)  +w0*rj2(i,k,3) 
      c         = sqrt(gamma*p(i,j-1,k)/w(i,j-1,k,1)) 
      er        = qnf + 2.*dgm1*c
      fr        = qn0 - 2.*dgm1*c0
      c         = .25*(er - fr)*gm1
      qn        = .50*(er + fr)
c 
c--------- inflow boundary (subsonic)
c 
      if (qn.le.0.) then
        uf        = u0 + (qn - qn0)*rj2(i,k,1)
        vf        = v0 + (qn - qn0)*rj2(i,k,2)
        wf        = w0 + (qn - qn0)*rj2(i,k,3)
        s         = s0
c 
c--------- outflow boundary (subsonic) 
c 
      else
        uf        = uf + (qn - qnf)*rj2(i,k,1)
        vf        = vf + (qn - qnf)*rj2(i,k,2)
        wf        = wf + (qn - qnf)*rj2(i,k,3)
        s         = w(i,j-1,k,1)**gamma/p(i,j-1,k)
      endif 
c 
      cc          = c*c/gamma
      w(i,j,k,1)  = (s*cc)**dgm1
      w(i,j,k,2)  = w(i,j,k,1)*uf
      w(i,j,k,3)  = w(i,j,k,1)*vf
      w(i,j,k,4)  = w(i,j,k,1)*wf
      p(i,j,k)    = w(i,j,k,1)*cc
c
c        energy is overspecified
c     w(i,j,k,5) = -p(i,j,k) + h0*w(i,j,k,1)
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  405 continue
c
      else
c
c-----   subsonic free stream (with preconditioning)
c
cdir$ ivdep
      do 408 i= n2beg, n2end
      uf        = w(i,j-1,k,2)/w(i,j-1,k,1) 
      vf        = w(i,j-1,k,3)/w(i,j-1,k,1) 
      wf        = w(i,j-1,k,4)/w(i,j-1,k,1) 
      qnf       = uf*rj2(i,k,1) + vf*rj2(i,k,2) + wf*rj2(i,k,3) 
c 
c--------- inflow boundary (subsonic)
c 
      if (qnf.le.0.0) then
        uf        = u0
        vf        = v0
        wf        = w0
        Tf        = p0/rho0
        p(i,j,k)  = p(i,j-1,k)
        w(i,j,k,1)  = p(i,j,k)/Tf
c 
c--------- outflow boundary (subsonic) 
c 
      else
        uf        = w(i,j-1,k,2)/w(i,j-1,k,1)
        vf        = w(i,j-1,k,3)/w(i,j-1,k,1)
        wf        = w(i,j-1,k,4)/w(i,j-1,k,1)
        Tf        = p(i,j-1,k)  /w(i,j-1,k,1)
        p(i,j,k)  = p0
        w(i,j,k,1)  = p(i,j,k)/Tf
      endif 
c 
      w(i,j,k,2)  = w(i,j,k,1)*uf
      w(i,j,k,3)  = w(i,j,k,1)*vf
      w(i,j,k,4)  = w(i,j,k,1)*wf
c
c       compute energy
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k)*dgm1 + .5*vtots/w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  408 continue
c     end of subsonic inflow/outflow condition
      endif
  410 continue
      else
c        supersonic flow condition
         do 415 k= n1beg, n1end
         do 415 i= n2beg, n2end
         uf        = w(i,j-1,k,2)/w(i,j-1,k,1) 
         vf        = w(i,j-1,k,3)/w(i,j-1,k,1) 
         wf        = w(i,j-1,k,4)/w(i,j-1,k,1) 
         qnf       = uf*rj2(i,k,1)  +vf*rj2(i,k,2)  +wf*rj2(i,k,3) 
         qn0       = u0*rj2(i,k,1)  +v0*rj2(i,k,2)  +w0*rj2(i,k,3) 
         c         = sqrt(gamma*p(i,j-1,k)/w(i,j-1,k,1)) 
         er        = qnf + 2.*dgm1*c
         fr        = qn0 - 2.*dgm1*c0
         c         = .25*(er - fr)*gm1
         qn        = .50*(er + fr)
c 
c--------- inflow boundary (supersonic)
c 
      if (qn.le.0.) then
         w(i,j,k,1)  = rho0
         w(i,j,k,2)  = rho0*u0
         w(i,j,k,3)  = rho0*v0
         w(i,j,k,4)  = rho0*w0
         w(i,j,k,5)  = rho0*h0 -p0
         p(i,j,k)    = p0
         eomu(i,j,k) = evinf
       else
c
c------- outflow boundary (supersonic)
c 
         w(i,j,k,1)  = w(i,j-1,k,1)
         w(i,j,k,2)  = w(i,j-1,k,2)
         w(i,j,k,3)  = w(i,j-1,k,3)
         w(i,j,k,4)  = w(i,j-1,k,4)
         w(i,j,k,5)  = w(i,j-1,k,5)
         p(i,j,k)    = p(i,j-1,k)
         eomu(i,j,k) = eomu(i,j-1,k)
c
       endif
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
c----- mods by vatsa (10-3-94) to compute only interior cells
      if( n1beg .ne. imn) n1beg = n1beg +1
      if( n1end .eq. imx) n1end = n1end +1
      if( n2beg .ne. jmn) n2beg = n2beg +1
      if( n2end .eq. jmx) n2end = n2end +1
c
      k   =   kmn
c
      if(rm.lt.1.0) then
      do 510 j= n2beg, n2end
      if(iiprec.le.0) then
c-----   subsonic free stream (no preconditioning)
cdir$ ivdep
      do 505 i= n1beg, n1end
      uf        = w(i,j,k+1,2)/w(i,j,k+1,1) 
      vf        = w(i,j,k+1,3)/w(i,j,k+1,1) 
      wf        = w(i,j,k+1,4)/w(i,j,k+1,1) 
      qnf       = uf*rk1(i,j,1)  +vf*rk1(i,j,2)  +wf*rk1(i,j,3) 
      qn0       = u0*rk1(i,j,1)  +v0*rk1(i,j,2)  +w0*rk1(i,j,3) 
      c         = sqrt(gamma*p(i,j,k+1)/w(i,j,k+1,1)) 
      er        = qnf - 2.*dgm1*c
      fr        = qn0 + 2.*dgm1*c0
      c         = .25*(er - fr)*gm1
      qn        = .50*(er + fr)
c 
c--------- inflow boundary (subsonic)
c 
      if (qn.gt.0.) then
        uf        = u0 + (qn - qn0)*rk1(i,j,1)
        vf        = v0 + (qn - qn0)*rk1(i,j,2)
        wf        = w0 + (qn - qn0)*rk1(i,j,3)
        s         = s0
c 
c--------- outflow boundary (subsonic)
c 
      else
        uf        = uf + (qn - qnf)*rk1(i,j,1)
        vf        = vf + (qn - qnf)*rk1(i,j,2)
        wf        = wf + (qn - qnf)*rk1(i,j,3)
        s         = w(i,j,k+1,1)**gamma/p(i,j,k+1)
      endif 
c 
      cc          = c**2/gamma
      w(i,j,k,1)  = (s*cc)**dgm1
      w(i,j,k,2)  = w(i,j,k,1)*uf
      w(i,j,k,3)  = w(i,j,k,1)*vf
      w(i,j,k,4)  = w(i,j,k,1)*wf
      p(i,j,k)    = w(i,j,k,1)*cc
c
c        specify energy 
c     w(i,j,k,5) = -p(i,j,k) + h0*w(i,j,k,1)
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  505 continue
c
      else
c
c-----   subsonic free stream (with preconditioning)
c
cdir$ ivdep
      do 508 i= n1beg, n1end
      uf        = w(i,j,k+1,2)/w(i,j,k+1,1)
      vf        = w(i,j,k+1,3)/w(i,j,k+1,1)
      wf        = w(i,j,k+1,4)/w(i,j,k+1,1)
      qnf       = uf*rk1(i,j,1) + vf*rk1(i,j,2) + wf*rk1(i,j,3)
c
c--------- inflow boundary (subsonic)
c
      if (qnf.gt.0.0) then
        uf        = u0
        vf        = v0
        wf        = w0
        Tf        = p0/rho0
        p(i,j,k)  = p(i,j,k+1)
        w(i,j,k,1)  = p(i,j,k)/Tf
c
c--------- outflow boundary (subsonic)
c
      else
        uf        = w(i,j,k+1,2)/w(i,j,k+1,1)
        vf        = w(i,j,k+1,3)/w(i,j,k+1,1)
        wf        = w(i,j,k+1,4)/w(i,j,k+1,1)
        Tf        = p(i,j,k+1)  /w(i,j,k+1,1)
        p(i,j,k)  = p0
        w(i,j,k,1)  = p(i,j,k)/Tf
      endif
c
      w(i,j,k,2)  = w(i,j,k,1)*uf
      w(i,j,k,3)  = w(i,j,k,1)*vf
      w(i,j,k,4)  = w(i,j,k,1)*wf
c
c       compute energy
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k)*dgm1 + .5*vtots/w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  508 continue
c     end of subsonic inflow/outflow condition
      endif
  510 continue
      else
c        supersonic flow condition
         do 515 j= n2beg, n2end
         do 515 i= n1beg, n1end
         uf        = w(i,j,k+1,2)/w(i,j,k+1,1) 
         vf        = w(i,j,k+1,3)/w(i,j,k+1,1) 
         wf        = w(i,j,k+1,4)/w(i,j,k+1,1) 
         qnf       = uf*rk1(i,j,1)  +vf*rk1(i,j,2)  +wf*rk1(i,j,3) 
         qn0       = u0*rk1(i,j,1)  +v0*rk1(i,j,2)  +w0*rk1(i,j,3) 
         c         = sqrt(gamma*p(i,j,k+1)/w(i,j,k+1,1)) 
         er        = qnf - 2.*dgm1*c
         fr        = qn0 + 2.*dgm1*c0
         c         = .25*(er - fr)*gm1
         qn        = .50*(er + fr)
c 
c--------- inflow boundary (supersonic)
c 
      if (qn.gt.0.) then
         w(i,j,k,1)  = rho0
         w(i,j,k,2)  = rho0*u0
         w(i,j,k,3)  = rho0*v0
         w(i,j,k,4)  = rho0*w0
         w(i,j,k,5)  = rho0*h0 -p0
         p(i,j,k)    = p0
         eomu(i,j,k) = evinf

      else
c
c------- outflow boundary (supersonic)
c 
         w(i,j,k,1)  = w(i,j,k+1,1)
         w(i,j,k,2)  = w(i,j,k+1,2)
         w(i,j,k,3)  = w(i,j,k+1,3)
         w(i,j,k,4)  = w(i,j,k+1,4)
         w(i,j,k,5)  = w(i,j,k+1,5)
         p(i,j,k)    = p(i,j,k+1)
         eomu(i,j,k) = eomu(i,j,k+1)
c
       endif
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
c----- mods by vatsa (10-3-94) to compute only interior cells
      if( n1beg .ne. imn) n1beg = n1beg +1
      if( n1end .eq. imx) n1end = n1end +1
      if( n2beg .ne. jmn) n2beg = n2beg +1
      if( n2end .eq. jmx) n2end = n2end +1
c
      k   =   kmx + 1
c
      if(rm.lt.1.0) then
      do 610 j= n2beg, n2end
      if(iiprec.le.0) then
c-----   subsonic free stream (no preconditioning)
cdir$ ivdep
      do 605 i= n1beg, n1end
      uf        = w(i,j,k-1,2)/w(i,j,k-1,1) 
      vf        = w(i,j,k-1,3)/w(i,j,k-1,1) 
      wf        = w(i,j,k-1,4)/w(i,j,k-1,1) 
      qnf       = uf*rk2(i,j,1)  +vf*rk2(i,j,2)  +wf*rk2(i,j,3) 
      qn0       = u0*rk2(i,j,1)  +v0*rk2(i,j,2)  +w0*rk2(i,j,3) 
      c         = sqrt(gamma*p(i,j,k-1)/w(i,j,k-1,1)) 
      er        = qnf + 2.*dgm1*c
      fr        = qn0 - 2.*dgm1*c0
      c         = .25*(er - fr)*gm1
      qn        = .50*(er + fr)
c 
c--------- inflow boundary (subsonic)
c 
      if (qn.le.0.) then
        uf        = u0 + (qn - qn0)*rk2(i,j,1)
        vf        = v0 + (qn - qn0)*rk2(i,j,2)
        wf        = w0 + (qn - qn0)*rk2(i,j,3)
        s         = s0
c 
c--------- outflow boundary (subsonic)
c 
      else
        uf        = uf + (qn - qnf)*rk2(i,j,1)
        vf        = vf + (qn - qnf)*rk2(i,j,2)
        wf        = wf + (qn - qnf)*rk2(i,j,3)
        s         = w(i,j,k-1,1)**gamma/p(i,j,k-1)
      endif 
c 
      cc          = c**2/gamma
      w(i,j,k,1)  = (s*cc)**dgm1
      w(i,j,k,2)  = w(i,j,k,1)*uf
      w(i,j,k,3)  = w(i,j,k,1)*vf
      w(i,j,k,4)  = w(i,j,k,1)*wf
      p(i,j,k)    = w(i,j,k,1)*cc
c
c        specify energy
c     w(i,j,k,5) = -p(i,j,k) + h0*w(i,j,k,1)
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  605 continue
      else
c
c-----   subsonic free stream (with preconditioning)
c
cdir$ ivdep
      do 608 i= n2beg, n2end
      uf        = w(i,j,k-1,2)/w(i,j,k-1,1)
      vf        = w(i,j,k-1,3)/w(i,j,k-1,1)
      wf        = w(i,j,k-1,4)/w(i,j,k-1,1)
      qnf       = uf*rk2(i,j,1) + vf*rk2(i,j,2) + wf*rk2(i,j,3)
c
c--------- inflow boundary (subsonic)
c
      if (qnf.le.0.0) then
        uf        = u0
        vf        = v0
        wf        = w0
        Tf        = p0/rho0
        p(i,j,k)  = p(i,j,k-1)
        w(i,j,k,1)  = p(i,j,k)/Tf
c
c--------- outflow boundary (subsonic)
c
      else
        uf        = w(i,j,k-1,2)/w(i,j,k-1,1)
        vf        = w(i,j,k-1,3)/w(i,j,k-1,1)
        wf        = w(i,j,k-1,4)/w(i,j,k-1,1)
        Tf        = p(i,j,k-1)  /w(i,j,k-1,1)
        p(i,j,k)  = p0
        w(i,j,k,1)  = p(i,j,k)/Tf
      endif
c
      w(i,j,k,2)  = w(i,j,k,1)*uf
      w(i,j,k,3)  = w(i,j,k,1)*vf
      w(i,j,k,4)  = w(i,j,k,1)*wf
c
c       compute energy
c
      vtots        =      w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                    w(i,j,k,4)**2
      w(i,j,k,5)   =      p(i,j,k)*dgm1 + .5*vtots/w(i,j,k,1)
      eomu(i,j,k)  =      evinf
  608 continue
c     end of subsonic inflow/outflow condition
      endif
  610 continue
      else
c        supersonic flow condition
         do 615 j= n2beg, n2end
         do 615 i= n1beg, n1end
         uf        = w(i,j,k-1,2)/w(i,j,k-1,1) 
         vf        = w(i,j,k-1,3)/w(i,j,k-1,1) 
         wf        = w(i,j,k-1,4)/w(i,j,k-1,1) 
         qnf       = uf*rk2(i,j,1)  +vf*rk2(i,j,2)  +wf*rk2(i,j,3) 
         qn0       = u0*rk2(i,j,1)  +v0*rk2(i,j,2)  +w0*rk2(i,j,3) 
         c         = sqrt(gamma*p(i,j,k-1)/w(i,j,k-1,1)) 
         er        = qnf + 2.*dgm1*c
         fr        = qn0 - 2.*dgm1*c0
         c         = .25*(er - fr)*gm1
         qn        = .50*(er + fr)
c 
c--------- inflow boundary (supersonic)
c 
      if (qn.le.0.) then
         w(i,j,k,1)  = rho0
         w(i,j,k,2)  = rho0*u0
         w(i,j,k,3)  = rho0*v0
         w(i,j,k,4)  = rho0*w0
         w(i,j,k,5)  = rho0*h0 -p0
         p(i,j,k)    = p0
         eomu(i,j,k) = evinf

      else
c
c------- outflow boundary (supersonic)
c 
         w(i,j,k,1)  = w(i,j,k-1,1)
         w(i,j,k,2)  = w(i,j,k-1,2)
         w(i,j,k,3)  = w(i,j,k-1,3)
         w(i,j,k,4)  = w(i,j,k-1,4)
         w(i,j,k,5)  = w(i,j,k-1,5)
         p(i,j,k)    = p(i,j,k-1)
         eomu(i,j,k) = eomu(i,j,k-1)
c
       endif
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
