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
      subroutine advflow(imn, jmn, kmn, imx, jmx, kmx,
     .                   imp1,jmp1,kmp1,imp2,jmp2,kmp2,ws,delw,w,p,
     .                   dtl,ft,gt,dtphy,ntorder ,fbeta2)
c
c     object : Advance the flow variables in time
c              during the Runge-Kutta stages
c              w (i,j,k,.) is current  and new RK level
c              ws(i,j,k,.) is original RK level (k=0)
c              delw is residual multiplied by alpha_k delta-tau P
c
c     initial coding : v.n.vatsa of nasa langley, Feb. 1990
c            mod     : June 1995, to accommodate preconditioning
c     latest mods    : Aug  2000, to accommodate time accuracy
c
c
c-------   preconditioning related quantities      -----------------
c
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubcon2,iprec
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
      common/timek/  ct(9)
      common/unsprec/ npcons,nexp,iprecg,ijac,icord,inav
c
      dimension delw(imp2,jmp2,kmp2,5), ws(imp2,jmp2,kmp2,5),
     .             w(imp2,jmp2,kmp2,5),  p(imp2,jmp2,kmp2),
     .             dtl(imp2,jmp2,kmp2)
c
      dimension    fbeta2(imp2,jmp2,kmp2)
c
c        nexp relevant only for dual time stepping
c        nexp=1 explicit treatment of physical time dependent RHS
c        nexp=0 implicit treatment of physical time dependent RHS
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
      delmx         = 0.250
      tol           = 1.e-12
      delta         = deltpre
      if (iprecg .lt.0) delta=0.
      cp            = gamma*dgm1
      dcp           = gm1/gamma
 
      if (iprec.le.1) go to 100
c
c            iprec=2
c            (p,u,v,w,T) variables
c
      do 10 k=kmn+1,kmx
      do 10 j=jmn+1,jmx
      do 10 i=imn+1,imx
      rho           = ws(i,j,k,1)
      rhoi          = 1./rho
      uc            = ws(i,j,k,2)*rhoi
      vc            = ws(i,j,k,3)*rhoi
      wc            = ws(i,j,k,4)*rhoi
      qq            = .5*(uc**2 + vc**2 + wc**2)
      ps            = gm1*(ws(i,j,k,5) - rho*qq)
      ts            = ps*rhoi
c
c     make sure p and T do not change by more than
c     delmx% of its previous value anywhere in the field
c
      dp            = delw(i,j,k,1)
      dratio        = abs(dp/ps) + tol
      check         = dratio - delmx
      sfac          = delmx/dratio
      sfac1         = cvmgm(1.,sfac,check)
c
      dt            = delw(i,j,k,5)
      dratio        = abs(dt/ts) + tol
      check         = dratio - delmx
      sfac          = delmx/dratio
      sfac2         = cvmgm(1.,sfac,check)
c
      tmp1          = min(sfac1,sfac2)
      delw(i,j,k,1) = delw(i,j,k,1)*tmp1
      delw(i,j,k,2) = delw(i,j,k,2)*tmp1
      delw(i,j,k,3) = delw(i,j,k,3)*tmp1
      delw(i,j,k,4) = delw(i,j,k,4)*tmp1
      delw(i,j,k,5) = delw(i,j,k,5)*tmp1
   10 continue
c
      if (npcons.eq.0) then
c
c       no modification of time accurate terms
c       physical time derivative treated explicitly
c       this requires modification to artificial time step in blkstp
c
      do 20 k=kmn+1,kmx
      do 20 j=jmn+1,jmx
      do 20 i=imn+1,imx
        rho         = ws(i,j,k,1)
        rhoi        = 1./rho
        uc          = ws(i,j,k,2)*rhoi
        vc          = ws(i,j,k,3)*rhoi
        wc          = ws(i,j,k,4)*rhoi
        qq          = .5*(uc**2 + vc**2 + wc**2)
        ps          = gm1*(ws(i,j,k,5) - rho*qq)
        ts          = ps*rhoi
c
        pnew        = ps - delw(i,j,k,1)
        unew        = uc - delw(i,j,k,2)
        vnew        = vc - delw(i,j,k,3)
        wnew        = wc - delw(i,j,k,4)
        tnew        = ts - delw(i,j,k,5)
c
        w(i,j,k,1)  = pnew/tnew
        w(i,j,k,2)  = unew*w(i,j,k,1)
        w(i,j,k,3)  = vnew*w(i,j,k,1)
        w(i,j,k,4)  = wnew*w(i,j,k,1)
        w(i,j,k,5)  = pnew*dgm1 
     .              + .5*w(i,j,k,1)*(unew**2 + vnew**2 + wnew**2)
        p(i,j,k)    = pnew
   20 continue
      endif ! npcons=0
c
      if (npcons.eq.1) then
c         npcons=1   convert to conservation variables
      if (ntorder.eq.0 .or. nexp.eq.1) then
c         explicit treatment of physical time derivative
c                 no correction term
      do 30 k=kmn+1,kmx
      do 30 j=jmn+1,jmx
      do 30 i=imn+1,imx
        rho         = ws(i,j,k,1)
        rhoi        = 1./rho
        uc          = ws(i,j,k,2)*rhoi
        vc          = ws(i,j,k,3)*rhoi
        wc          = ws(i,j,k,4)*rhoi
        qq          = .5*(uc**2 + vc**2 + wc**2)
        ps          = gm1*(ws(i,j,k,5) - rho*qq)
        ts          = ps*rhoi
c
c       use linearization to directly calculate conservation variables
c
        drho        = rho*(delw(i,j,k,1)/ps - delw(i,j,k,5)/ts)
        drhou       = uc*drho + rho*delw(i,j,k,2)
        drhov       = vc*drho + rho*delw(i,j,k,3)
        drhow       = wc*drho + rho*delw(i,j,k,4)
        dE          = delw(i,j,k,1)*dgm1 + qq*drho
     .              + ws(i,j,k,2)*delw(i,j,k,2)
     .              + ws(i,j,k,3)*delw(i,j,k,3)
     .              + ws(i,j,k,4)*delw(i,j,k,4)
c
        w(i,j,k,1)  = ws(i,j,k,1) - drho
        w(i,j,k,2)  = ws(i,j,k,2) - drhou
        w(i,j,k,3)  = ws(i,j,k,3) - drhov
        w(i,j,k,4)  = ws(i,j,k,4) - drhow
        w(i,j,k,5)  = ws(i,j,k,5) - dE
        qq          = (w(i,j,k,2)**2 + w(i,j,k,3)**2
     .              +  w(i,j,k,4)**2 )/w(i,j,k,1)
        p(i,j,k)    = gm1*(w(i,j,k,5) - .5*qq )
        p(i,j,k)    = max(p(i,j,k),0.001)
   30 continue
      else
c               nexp=0
c         convert to conservation variables
c         implicit treatment of physical time derivative
c         ntorder.gt.0 .and. nexp.eq.0
      do 40 k=kmn+1,kmx
      do 40 j=jmn+1,jmx
      do 40 i=imn+1,imx
        rhoc        = ws(i,j,k,1)
        rhoi        = 1./rhoc
        uc          = ws(i,j,k,2)*rhoi
        vc          = ws(i,j,k,3)*rhoi
        wc          = ws(i,j,k,4)*rhoi
        qq          = .5*(uc**2 + vc**2 + wc**2)
        ps          = gm1*(ws(i,j,k,5) - rhoc*qq)
        ts          = ps*rhoi
        c2          = gamma*ts
        hc          = cp*ts + qq
        bc          = fbeta2(i,j,k)
        if (iprecg.le.0) bc = 1.
c          precondition w on RHS by P_C
        y1          = gm1*(qq*w(i,j,k,1) - uc*w(i,j,k,2)
     .              - vc*w(i,j,k,3) - wc*w(i,j,k,4)
     .              + w(i,j,k,5) )/c2
c
        dw1         = (bc-1.)*y1 + delta*(y1-w(i,j,k,1))
        w1          = w(i,j,k,1) + dw1
        w2          = w(i,j,k,2) + uc*dw1
        w3          = w(i,j,k,3) + vc*dw1
        w4          = w(i,j,k,4) + wc*dw1
        w5          = w(i,j,k,5) + hc*dw1
c
c       use linearization to calculate conservation residuals
c              from (p,u,v,w,T) residuals
c
        drho        = rhoc*(delw(i,j,k,1)/ps - delw(i,j,k,5)/ts)
        drhou       = uc*drho + rho*delw(i,j,k,2)
        drhov       = vc*drho + rho*delw(i,j,k,3)
        drhow       = wc*drho + rho*delw(i,j,k,4)
        dE          = delw(i,j,k,1)*dgm1 + qq*drho
     .              + ws(i,j,k,2)*delw(i,j,k,2)
     .              + ws(i,j,k,3)*delw(i,j,k,3)
     .              + ws(i,j,k,4)*delw(i,j,k,4)
c
c           similar to loop 150
c
        dt          = (ft*dtl(i,j,k) + gt)*ct(ntorder+1)/dtphy
        w1          = ws(i,j,k,1) - drho  + w1*dt
        w2          = ws(i,j,k,2) - drhou + w2*dt
        w3          = ws(i,j,k,3) - drhov + w3*dt
        w4          = ws(i,j,k,4) - drhow + w4*dt
        w5          = ws(i,j,k,5) - dE    + w5*dt
c
c          precondition w on LHS by (I + dt*P_C)
c          i.e. multiply RHS by (I + dt*P_C)^(-1)
c
c            "correct formula"
        dw1         = ((1.-bc)*dt*y1 - delta*dt*(y1-w1))/(1.+bc*dt)
        ddt         = 1./(1.+dt)
        w(i,j,k,1)  = (w1 +    dw1)*ddt
        w(i,j,k,2)  = (w2 + uc*dw1)*ddt
        w(i,j,k,3)  = (w3 + vc*dw1)*ddt
        w(i,j,k,4)  = (w4 + wc*dw1)*ddt
        w(i,j,k,5)  = (w5 + hc*dw1)*ddt
        qq          = (w(i,j,k,2)**2 + w(i,j,k,3)**2
     .              +  w(i,j,k,4)**2 )/w(i,j,k,1)
        p(i,j,k)    =  gm1*(w(i,j,k,5) - .5*qq)
        p(i,j,k)    = max(p(i,j,k),0.001)
   40 continue
      endif ! nexp
      endif ! npcons=1
c
      if (npcons.eq.2) then
         if (nexp.eq.1) then
c        explicit treatment in primitive variables
      do 50 k=kmn+1,kmx
      do 50 j=jmn+1,jmx
      do 50 i=imn+1,imx
        rhoc        = ws(i,j,k,1)
        rhoi        = 1./rhoc
        uc          = ws(i,j,k,2)*rhoi
        vc          = ws(i,j,k,3)*rhoi
        wc          = ws(i,j,k,4)*rhoi
        qq          = .5*(uc**2 + vc**2 + wc**2)
        ps          = gm1*(ws(i,j,k,5) - rhoc*qq)
        ts          = ps*rhoi
c
        pnew        = ps - delw(i,j,k,1)
        unew        = uc - delw(i,j,k,2)
        vnew        = vc - delw(i,j,k,3)
        wnew        = wc - delw(i,j,k,4)
        tnew        = ts - delw(i,j,k,5)
c
        w(i,j,k,1)  = pnew/tnew
        w(i,j,k,2)  = unew*w(i,j,k,1)
        w(i,j,k,3)  = vnew*w(i,j,k,1)
        w(i,j,k,4)  = wnew*w(i,j,k,1)
        w(i,j,k,5)  = pnew*dgm1 
     .              + .5*w(i,j,k,1)*(unew**2 + vnew**2 + wnew**2)
        p(i,j,k)    = pnew
   50 continue
      else
c        nexp=0  implicit treatment in primitive variables
      do 90 k=kmn+1,kmx
      do 90 j=jmn+1,jmx
      do 90 i=imn+1,imx
        rhoc        = ws(i,j,k,1)
        rhoi        = 1./rhoc
        uc          = ws(i,j,k,2)*rhoi
        vc          = ws(i,j,k,3)*rhoi
        wc          = ws(i,j,k,4)*rhoi
        qq          = .5*(uc**2 + vc**2 + wc**2)
        ps          = gm1*(ws(i,j,k,5) - rhoc*qq)
        ts          = ps*rhoi
c
        rhos        = 1./w(i,j,k,1)
        uu          = w(i,j,k,2)*rhos
        vv          = w(i,j,k,3)*rhos
        ww          = w(i,j,k,4)*rhos
        pp          = p(i,j,k)
        tt          = pp*rhos
c
        bc          = fbeta2(i,j,k)
        if (iprecg.le.0) bc = 1.
c
c          precondition Q on RHS by P_Q
c
        w1          = bc*pp
        w2          = uu
        w3          = vv
        w4          = ww
        w5          = tt + dcp*(bc-1.)*rhoi*pp
c
        dt          = (ft*dtl(i,j,k) + gt)*ct(ntorder+1)/dtphy
        pnew        = ps - delw(i,j,k,1) + w1*dt
        unew        = uc - delw(i,j,k,2) + w2*dt
        vnew        = vc - delw(i,j,k,3) + w3*dt
        wnew        = wc - delw(i,j,k,4) + w4*dt
        tnew        = ts - delw(i,j,k,5) + w5*dt
c
c          precondition w on LHS by (I + dt*P_Q)
c          i.e. multiply RHS by (I + dt*P_Q)^(-1)
c
        ddt         = 1./(1.+dt)
        e           = dt*(1.-bc)/(1.+dt*bc)
        unew        = unew*ddt
        vnew        = vnew*ddt
        wnew        = wnew*ddt
        tnew        = (tnew + dcp*rhoi*e*pnew)*ddt
        pnew        = (1. + e)*pnew*ddt
c
        w(i,j,k,1)  = pnew/tnew
        w(i,j,k,2)  = unew*w(i,j,k,1)
        w(i,j,k,3)  = vnew*w(i,j,k,1)
        w(i,j,k,4)  = wnew*w(i,j,k,1)
        w(i,j,k,5)  = pnew*dgm1 
     .              + .5*w(i,j,k,1)*(unew**2 + vnew**2 + wnew**2)
        p(i,j,k)    = pnew
   90 continue
      endif ! nexp
      endif ! npcons=2
  
      return
c
  100 continue
c
c        (iprec.le.1) conservation variables
c
c     make sure rho does not change by more than
c     delmx% of its previous value anywhere in the field
c
      do 110 k =kmn+1,kmx
      do 110 j =jmn+1,jmx
      do 110 i =imn+1,imx
      rho           = ws(i,j,k,1)
      rhoi          = 1./rho
      ut            = ws(i,j,k,2)*rhoi
      vt            = ws(i,j,k,3)*rhoi
      wt            = ws(i,j,k,4)*rhoi
      drho          = delw(i,j,k,1)
c
      dratio        = abs(drho*rhoi) + tol
      check         = dratio - delmx
      sfac          = delmx/dratio
      sfac1         = cvmgm(1.,sfac,check)
c
      dratio        = abs(delw(i,j,k,5)/ws(i,j,k,5)) + tol
      check         = dratio - delmx
      sfac          = delmx/dratio
      sfac2         = cvmgm(1.,sfac,check)
      sfac1         = min(sfac1,sfac2)
c
      qq            = 0.5*(ut*ut +vt*vt +wt*wt)
      ps            = gm1*(ws(i,j,k,5) - rho*qq)
      dp            = gm1*(qq*delw(i,j,k,1) - ut*delw(i,j,k,2)
     .              - vt*delw(i,j,k,3) - wt*delw(i,j,k,4)
     .              + delw(i,j,k,5))
c
      dratio        = abs(dp/ps) + tol
      check         = dratio - delmx
      sfac          = delmx/dratio
      sfac2         = cvmgm(1.,sfac,check)
c
      tmp1          = min(sfac1,sfac2)
      delw(i,j,k,1) = delw(i,j,k,1)*tmp1
      delw(i,j,k,2) = delw(i,j,k,2)*tmp1
      delw(i,j,k,3) = delw(i,j,k,3)*tmp1
      delw(i,j,k,4) = delw(i,j,k,4)*tmp1
      delw(i,j,k,5) = delw(i,j,k,5)*tmp1
  110 continue
c
      if (ntorder.eq.0 .or. npcons.eq.0) then
c
c       no modification of time accurate terms
c       physical time derivative treated explicitly
c       this requires modification to artificial time step in blkstp
c
        do k  = kmn+1,kmx
        do j  = jmn+1,jmx
        do i  = imn+1,imx
          w(i,j,k,1) = ws(i,j,k,1) - delw(i,j,k,1)
          w(i,j,k,2) = ws(i,j,k,2) - delw(i,j,k,2)
          w(i,j,k,3) = ws(i,j,k,3) - delw(i,j,k,3)
          w(i,j,k,4) = ws(i,j,k,4) - delw(i,j,k,4)
          w(i,j,k,5) = ws(i,j,k,5) - delw(i,j,k,5)
          qq         = (w(i,j,k,2)**2 + w(i,j,k,3)**2 
     .               +  w(i,j,k,4)**2 )/w(i,j,k,1)
          p(i,j,k)   =  gm1*(w(i,j,k,5) - .5*qq)
          p(i,j,k)   = max(p(i,j,k),0.001)
        enddo
        enddo
        enddo
      else
c
c        correction to make physical time derivative implicit
c
      if (iprec.eq.0 .or. iprecg.lt.0) then
c                  no preconditioning
c         implicit treatment of physical time step
        do k  = kmn+1,kmx
        do j  = jmn+1,jmx
        do i  = imn+1,imx
        dt          = (ft*dtl(i,j,k) + gt)*ct(ntorder+1)/dtphy
        ddt         = 1./(1.+dt)
        w(i,j,k,1)  = (ws(i,j,k,1) - delw(i,j,k,1)
     .              +  w(i,j,k,1)*dt)*ddt
        w(i,j,k,2)  = (ws(i,j,k,2) - delw(i,j,k,2)
     .              +  w(i,j,k,2)*dt)*ddt
        w(i,j,k,3)  = (ws(i,j,k,3) - delw(i,j,k,3)
     .              +  w(i,j,k,3)*dt)*ddt
        w(i,j,k,4)  = (ws(i,j,k,4) - delw(i,j,k,4)
     .              +  w(i,j,k,4)*dt)*ddt
        w(i,j,k,5)  = (ws(i,j,k,5) - delw(i,j,k,5)
     .              +  w(i,j,k,5)*dt)*ddt
        qq          = (w(i,j,k,2)**2 +w(i,j,k,3)**2 +
     .                 w(i,j,k,4)**2 )/w(i,j,k,1)
        p(i,j,k)    =  gm1*(w(i,j,k,5) - .5*qq)
        p(i,j,k)    = max(p(i,j,k),0.001)
        enddo
        enddo
        enddo
      else
c              iprec =1
c       time derivative includes preconditioning
c              assumes alpha  = 0
      do 150 k =kmn+1,kmx
      do 150 j =jmn+1,jmx
      do 150 i =imn+1,imx
        rhoi        = 1./ws(i,j,k,1)
        uc          = ws(i,j,k,2)*rhoi
        vc          = ws(i,j,k,3)*rhoi
        wc          = ws(i,j,k,4)*rhoi
        qq          = 0.5*(uc**2 + vc**2 + wc**2)
        ts          = p(i,j,k)*rhoi
        c2          = gamma*ts
        hc          = cp*ts + qq
        bc          = fbeta2(i,j,k)
c       if (iprecg.le.0) bc = 1.
c          precondition w on RHS by P_C
        y1          = gm1*(qq*w(i,j,k,1) - uc*w(i,j,k,2)
     .              - vc*w(i,j,k,3) - wc*w(i,j,k,4)
     .              + w(i,j,k,5) )/c2
c
        dw1         = (bc-1.)*y1 + delta*(y1-w(i,j,k,1))
        w1          = w(i,j,k,1) + dw1
        w2          = w(i,j,k,2) + uc*dw1
        w3          = w(i,j,k,3) + vc*dw1
        w4          = w(i,j,k,4) + wc*dw1
        w5          = w(i,j,k,5) + hc*dw1
c
        dt          = (ft*dtl(i,j,k) + gt)*ct(ntorder+1)/dtphy
        w1          = ws(i,j,k,1) - delw(i,j,k,1) + w1*dt
        w2          = ws(i,j,k,2) - delw(i,j,k,2) + w2*dt
        w3          = ws(i,j,k,3) - delw(i,j,k,3) + w3*dt
        w4          = ws(i,j,k,4) - delw(i,j,k,4) + w4*dt
        w5          = ws(i,j,k,5) - delw(i,j,k,5) + w5*dt
c
c          precondition w on LHS by (I + dt*P_C)
c          i.e. multiply RHS by (I + dt*P_C)^(-1)
c
c       dw1         = ((1.-bc)*dt*y1 - delta*dt*(y1-w1)/(1.+bc*dt))
c       dw1         = ((1.-bc)*dt*y1 - delta*dt*(y1-w1))/(1.+bc*dt)
c                correct formula   ----- doesn't seem to work
CCC     dw1         = (1.-bc)*dt*y1/(1.+bc*dt)
c                kludge
        dw1         = (1.-bc)*dt*y1
c       bct         = bc*dt
        ddt         = 1./(1.+dt)
        w(i,j,k,1)  = (w1 +    dw1)*ddt
        w(i,j,k,2)  = (w2 + uc*dw1)*ddt
        w(i,j,k,3)  = (w3 + vc*dw1)*ddt
        w(i,j,k,4)  = (w4 + wc*dw1)*ddt
        w(i,j,k,5)  = (w5 + hc*dw1)*ddt
        qq          = (w(i,j,k,2)**2 + w(i,j,k,3)**2
     .              +  w(i,j,k,4)**2 )/w(i,j,k,1)
        p(i,j,k)    = gm1*(w(i,j,k,5) - .5*qq)
        p(i,j,k)    = max(p(i,j,k),0.001)
  150 continue
      end if ! ntorder
      end if ! iprec
c
      return
      end
