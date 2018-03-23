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
      subroutine  bcpttfl(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2,
     .                  w,p,x,eomu,turv1,nface,n1beg,n1end,
     .                  n2beg,n2end,igrid,isoln,extyp,
     .                  ptstr,ttstr,relax ,alfastr,betastr)
c=======================================================================
c
c     object : to compute the ghost cell values for flow variables
c              by imposing specified total pressure and temperature
c
c     initial coding : by v.n.vatsa, nasa langley ( Jan. 1998)
c
c     patterned after Ted Reyhner's work to handle flows into the
c     computational field given the total temperature and total pressure
c
c     it assumes that static pressure (in computational field) is correct
c     and adjust the other variables to give the proper total conditions
c
c     ptstr & ttstr can differ for each stream (multiple stream capability)
c     relax = underelaxation factor (not used currently, i.e. relax=1.)
c    (ptstr,ttstr,relax are read as const1,const2 and const3 from mapfile)
c    
c
c     coding may not be correct for j and k, but it was installed anyway
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
c------------------------------------------------------------------
c
c     extyp  =   0  zero order extrapolation
c                1  linear     extrapolation  (not yet coded)
c
c------------------------------------------------------------------
c
c     extyp = 0.0
c
c      check the permissibility of the 'extyp'
c
cBWW moved this test to imapin.f
cBWW  if(extyp .gt.1.01 .or. extyp .lt. -0.01 )   then
cBWW     write (iwrit,'(2x," extyp = ",f8.4, " not permissible")')
cBWW     stop
cBWW  endif
c
c     set constants so that they are not constantly recalculated
      gogm1i=1.0/gogm1
      gm2i=2.0/gm1
      fact1=relax
      fact2=1.0-relax
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
c     branch-cut for different type of boundary faces
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
      i   =   imn
c
c
        do 110 k= n2beg, n2end
        do 110 j= n1beg, n1end
c
        pstr    = min(p(i+1,j,k),0.999*ptstr)
        rmstrsq = gm2i*((ptstr/pstr)**gogm1i - 1.)
        tstr    = ttstr/(1.+.5*rmstrsq*gm1)
        rmstr   = sqrt(rmstrsq)
        rhostr  = pstr/tstr
        cstr    = sqrt(gamma*tstr)
        vtotstr = rmstr*cstr
c
        w(i,j,k,1)  = rhostr*fact1 + fact2*w(i,j,k,1)
        rq          = w(i,j,k,1)*vtotstr
        w(i,j,k,2)  = rq*cos(alfastr)*cos(betastr)*fact1 +
     .                fact2*w(i,j,k,2)
        w(i,j,k,3)  = rq*sin(alfastr)*cos(betastr)*fact1 +
     .                fact2*w(i,j,k,3)
        w(i,j,k,4)  = rq*sin(betastr)*fact1 +
     .                fact2*w(i,j,k,4)
        p(i,j,k)    = pstr
        w(i,j,k,5)  =(p(i,j,k) *dgm1 + .5*(vtotstr**2)*w(i,j,k,1))*
     .                fact1 + fact2*w(i,j,k,5)
        eomu(i,j,k) = evinf
  110   continue
c
c
        if( igrid .eq. isoln) then
c
c          second layer of ghost cells for fine grids
c
           i   =  imn - 1
c
           do 120 k= n2beg, n2end
           do 120 j= n1beg, n1end
           w(i,j,k,1)   =    w(i+1,j,k,1)
           w(i,j,k,2)   =    w(i+1,j,k,2)
           w(i,j,k,3)   =    w(i+1,j,k,3)
           w(i,j,k,4)   =    w(i+1,j,k,4)
           w(i,j,k,5)   =    w(i+1,j,k,5)
           p(i,j,k)     =    p(i+1,j,k)
           eomu(i,j,k)  =    eomu(i+1,j,k)
ccc        turv1(i,j,k) =   turv1(i+1,j,k)
  120      continue
        endif
c
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
        do 210 k= n2beg, n2end
        do 210 j= n1beg, n1end
c
        pstr    = min(p(i-1,j,k),0.999*ptstr)
        rmstrsq = gm2i*((ptstr/pstr)**gogm1i - 1.)
        tstr    = ttstr/(1.+.5*rmstrsq*gm1)
        rmstr   = sqrt(rmstrsq)
        rhostr  = pstr/tstr
        cstr    = sqrt(gamma*tstr)
        vtotstr = rmstr*cstr
c
        w(i,j,k,1)  = rhostr*fact1 + fact2*w(i,j,k,1)
        rq          = w(i,j,k,1)*vtotstr
        w(i,j,k,2)  = rq*cos(alfastr)*cos(betastr)*fact1 +
     .                fact2*w(i,j,k,2)
        w(i,j,k,3)  = rq*sin(alfastr)*cos(betastr)*fact1 +
     .                fact2*w(i,j,k,3)
        w(i,j,k,4)  = rq*sin(betastr)*fact1 +
     .                fact2*w(i,j,k,4)
        p(i,j,k)    = pstr
        w(i,j,k,5)  =(p(i,j,k) *dgm1 + .5*(vtotstr**2)*w(i,j,k,1))*
     .                fact1 + fact2*w(i,j,k,5)
        eomu(i,j,k) = evinf
  210   continue
c
c
        if( igrid .eq. isoln) then
c
c          second layer of ghost cells for fine grids
c
           i   =  imx + 2
c
           do 220 k= n2beg, n2end
           do 220 j= n1beg, n1end
           w(i,j,k,1)   =    w(i-1,j,k,1)
           w(i,j,k,2)   =    w(i-1,j,k,2)
           w(i,j,k,3)   =    w(i-1,j,k,3)
           w(i,j,k,4)   =    w(i-1,j,k,4)
           w(i,j,k,5)   =    w(i-1,j,k,5)
           p(i,j,k)     =    p(i-1,j,k)
           eomu(i,j,k)  =    eomu(i-1,j,k)
ccc        turv1(i,j,k) =   turv1(i-1,j,k)
  220      continue
        endif
c
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
c
c
        do 310 i= n2beg, n2end
        do 310 k= n1beg, n1end
c
        pstr    = min(p(i,j+1,k),0.999*ptstr)
        rmstrsq = gm2i*((ptstr/pstr)**gogm1i - 1.)
        tstr    = ttstr/(1.+.5*rmstrsq*gm1)
        rmstr   = sqrt(rmstrsq)
        rhostr  = pstr/tstr
        cstr    = sqrt(gamma*tstr)
        vtotstr = rmstr*cstr
c
        w(i,j,k,1)  = rhostr*fact1 + fact2*w(i,j,k,1)
        rq          = w(i,j,k,1)*vtotstr
        w(i,j,k,2)  = rq*cos(alfastr)*cos(betastr)*fact1 +
     .                fact2*w(i,j,k,2)
        w(i,j,k,3)  = rq*sin(alfastr)*cos(betastr)*fact1 +
     .                fact2*w(i,j,k,3)
        w(i,j,k,4)  = rq*sin(betastr)*fact1 +
     .                fact2*w(i,j,k,4)
        p(i,j,k)    = pstr
        w(i,j,k,5)  =(p(i,j,k) *dgm1 + .5*(vtotstr**2)*w(i,j,k,1))*
     .                fact1 + fact2*w(i,j,k,5)
        eomu(i,j,k) = evinf
  310   continue
c
        if( igrid .eq. isoln) then
c
c          second layer of ghost cells for fine grids
c
           j   =  jmn - 1
c
           do 320 k= n1beg, n1end
           do 320 i= n2beg, n2end
           w(i,j,k,1)   =    w(i,j+1,k,1)
           w(i,j,k,2)   =    w(i,j+1,k,2)
           w(i,j,k,3)   =    w(i,j+1,k,3)
           w(i,j,k,4)   =    w(i,j+1,k,4)
           w(i,j,k,5)   =    w(i,j+1,k,5)
           p(i,j,k)     =    p(i,j+1,k)
           eomu(i,j,k)  =    eomu(i,j+1,k)
ccc        turv1(i,j,k) =   turv1(i,j+1,k)
  320      continue
        endif
c
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
c
        do 410 k= n1beg, n1end
        do 410 i= n2beg, n2end
c
        pstr    = min(p(i,j-1,k),0.999*ptstr)
        rmstrsq = gm2i*((ptstr/pstr)**gogm1i - 1.)
        tstr    = ttstr/(1.+.5*rmstrsq*gm1)
        rmstr   = sqrt(rmstrsq)
        rhostr  = pstr/tstr
        cstr    = sqrt(gamma*tstr)
        vtotstr = rmstr*cstr
c
        w(i,j,k,1)  = rhostr*fact1 + fact2*w(i,j,k,1)
        rq          = w(i,j,k,1)*vtotstr
        w(i,j,k,2)  = rq*cos(alfastr)*cos(betastr)*fact1 +
     .                fact2*w(i,j,k,2)
        w(i,j,k,3)  = rq*sin(alfastr)*cos(betastr)*fact1 +
     .                fact2*w(i,j,k,3)
        w(i,j,k,4)  = rq*sin(betastr)*fact1 +
     .                fact2*w(i,j,k,4)
        p(i,j,k)    = pstr
        w(i,j,k,5)  =(p(i,j,k) *dgm1 + .5*(vtotstr**2)*w(i,j,k,1))*
     .                fact1 + fact2*w(i,j,k,5)
        eomu(i,j,k) = evinf
  410   continue

c
        if( igrid .eq. isoln) then
c
c          second layer of ghost cells for fine grids
c
           j   =  jmx + 2
c
           do 420 k= n1beg, n1end
           do 420 i= n2beg, n2end
           w(i,j,k,1)   =    w(i,j-1,k,1)
           w(i,j,k,2)   =    w(i,j-1,k,2)
           w(i,j,k,3)   =    w(i,j-1,k,3)
           w(i,j,k,4)   =    w(i,j-1,k,4)
           w(i,j,k,5)   =    w(i,j-1,k,5)
           p(i,j,k)     =    p(i,j-1,k)
           eomu(i,j,k)  =    eomu(i,j-1,k)
ccc        turv1(i,j,k) =   turv1(i,j-1,k)
  420      continue
        endif
c
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
c
        do 510 j= n2beg, n2end
        do 510 i= n1beg, n1end
c
        pstr    = min(p(i,j,k+1),0.999*ptstr)
        rmstrsq = gm2i*((ptstr/pstr)**gogm1i - 1.)
        tstr    = ttstr/(1.+.5*rmstrsq*gm1)
        rmstr   = sqrt(rmstrsq)
        rhostr  = pstr/tstr
        cstr    = sqrt(gamma*tstr)
        vtotstr = rmstr*cstr
c
        w(i,j,k,1)  = rhostr*fact1 + fact2*w(i,j,k,1)
        rq          = w(i,j,k,1)*vtotstr
        w(i,j,k,2)  = rq*cos(alfastr)*cos(betastr)*fact1 +
     .                fact2*w(i,j,k,2)
        w(i,j,k,3)  = rq*sin(alfastr)*cos(betastr)*fact1 +
     .                fact2*w(i,j,k,3)
        w(i,j,k,4)  = rq*sin(betastr)*fact1 +
     .                fact2*w(i,j,k,4)
        p(i,j,k)    = pstr
        w(i,j,k,5)  =(p(i,j,k) *dgm1 + .5*(vtotstr**2)*w(i,j,k,1))*
     .                fact1 + fact2*w(i,j,k,5)
        eomu(i,j,k) = evinf
  510   continue
c
        if( igrid .eq. isoln) then
c
c          second layer of ghost cells for fine grids
c
           k   =  kmn - 1
c
           do 520 j= n2beg, n2end
           do 520 i= n1beg, n1end
           w(i,j,k,1)   =    w(i,j,k+1,1)
           w(i,j,k,2)   =    w(i,j,k+1,2)
           w(i,j,k,3)   =    w(i,j,k+1,3)
           w(i,j,k,4)   =    w(i,j,k+1,4)
           w(i,j,k,5)   =    w(i,j,k+1,5)
           p(i,j,k)     =    p(i,j,k+1)
           eomu(i,j,k)  =    eomu(i,j,k+1)
ccc        turv1(i,j,k) =   turv1(i,j,k+1)
  520      continue
        endif
c
      return
c
c
  600 continue
c--------------------- k = kmx  face ----------------------------
c
      if( n1beg .ne. jmn) n1beg = n1beg +1
      if( n1end .eq. jmx) n1end = n1end +1
      if( n2beg .ne. kmn) n2beg = n2beg +1
      if( n2end .eq. kmx) n2end = n2end +1
c
      k   =   kmx + 1
c
c
        do 610 j= n2beg, n2end
        do 610 i= n1beg, n1end
c
        pstr    = min(p(i,j,k-1),0.999*ptstr)
        rmstrsq = gm2i*((ptstr/pstr)**gogm1i - 1.)
        tstr    = ttstr/(1.+.5*rmstrsq*gm1)
        rmstr   = sqrt(rmstrsq)
        rhostr  = pstr/tstr
        cstr    = sqrt(gamma*tstr)
        vtotstr = rmstr*cstr
c
        w(i,j,k,1)  = rhostr*fact1 + fact2*w(i,j,k,1)
        rq          = w(i,j,k,1)*vtotstr
        w(i,j,k,2)  = rq*cos(alfastr)*cos(betastr)*fact1 +
     .                fact2*w(i,j,k,2)
        w(i,j,k,3)  = rq*sin(alfastr)*cos(betastr)*fact1 +
     .                fact2*w(i,j,k,3)
        w(i,j,k,4)  = rq*sin(betastr)*fact1 +
     .                fact2*w(i,j,k,4)
        p(i,j,k)    = pstr
        w(i,j,k,5)  =(p(i,j,k) *dgm1 + .5*(vtotstr**2)*w(i,j,k,1))*
     .                fact1 + fact2*w(i,j,k,5)
        eomu(i,j,k) = evinf
  610   continue
c
        if( igrid .eq. isoln) then
c
c          second layer of ghost cells for fine grids
c
c
           k   =  kmx + 2
c
           do 620 j= n2beg, n2end
           do 620 i= n1beg, n1end
           w(i,j,k,1)   =    w(i,j,k-1,1)
           w(i,j,k,2)   =    w(i,j,k-1,2)
           w(i,j,k,3)   =    w(i,j,k-1,3)
           w(i,j,k,4)   =    w(i,j,k-1,4)
           w(i,j,k,5)   =    w(i,j,k-1,5)
           p(i,j,k)     =    p(i,j,k-1)
           eomu(i,j,k)  =    eomu(i,j,k-1)
ccc        turv1(i,j,k) =   turv1(i,j,k-1)
  620      continue
        endif
c
c
      return
      end
