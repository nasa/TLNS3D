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
      subroutine precong (imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .  w,p,delw,fbeta2)
c
c     convert fluxes from conserved variables to "p,u,v,w,t" variables
c     this is accomplished by multiplying fluxes with preconditioning
c     matrix gamma4
c     
c
c     initial coding : v.n.vatsa of nasa langley, June, 1995
c
c
c-------   preconditioning related quantities      -----------------
c
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubcon2,iprec
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb

      common/unsprec/ npcons,nexp,iprecg,ijac,icord,inav
c
      dimension      w(imp2,jmp2,kmp2,5), delw(imp2,jmp2,kmp2,5),
     .               p(imp2,jmp2,kmp2), six(imp1,jmp1,kmp1),
     .             siy(imp1,jmp1,kmp1), siz(imp1,jmp1,kmp1),
     .             sjx(imp1,jmp1,kmp1), sjy(imp1,jmp1,kmp1),
     .             sjz(imp1,jmp1,kmp1), skx(imp1,jmp1,kmp1),
     .             sky(imp1,jmp1,kmp1), skz(imp1,jmp1,kmp1)
c
      dimension    fbeta2(imp2,jmp2,kmp2)
c
      if (iprecg.le.0 .and. iprec.eq.1)  return
      cp         = gamma/gm1
      dcp        = gm1/gamma
      delta      = deltpre
      if (iprecg.le.0) delta=0.
c
c
      if (iprec.eq.1) go to 100
      if (delta.eq.1.0) go to 30
      if (delta.gt.0.0) go to 15
c            iprec=2     primitive variables
c     call transct (imn,jmn,kmn,imx,jmx,kmx,
c    .              imp1,jmp1,kmp1,imp2,jmp2,kmp2,w,p,delw )

c                delta=0
      do 10 k = kmn+1,kmx
      do 10 j = jmn+1,jmx
      do 10 i = imn+1,imx
        rhoi      = 1./w(i,j,k,1)
        uc        = w(i,j,k,2)*rhoi
        vc        = w(i,j,k,3)*rhoi
        wc        = w(i,j,k,4)*rhoi
        q2        = 0.5*(uc**2 + vc**2 + wc**2)
        ts        = p(i,j,k)*rhoi
        c2        = gamma*ts
c   
c       transfrom from conservative to (p,u,v,w,T) variables
c
        dw1       = gm1*(q2*delw(i,j,k,1) -uc*delw(i,j,k,2)
     .            -      vc*delw(i,j,k,3) -wc*delw(i,j,k,4)
     .            +         delw(i,j,k,5) )
        dw2       = (delw(i,j,k,2) -uc*delw(i,j,k,1))*rhoi
        dw3       = (delw(i,j,k,3) -vc*delw(i,j,k,1))*rhoi
        dw4       = (delw(i,j,k,4) -wc*delw(i,j,k,1))*rhoi
        dw5       = (dw1 -ts*delw(i,j,k,1))*rhoi
c
c        evaluate alpha and beta for preconditioning
c
        bc        = fbeta2(i,j,k)
c       bc        = max(fbeta2(i,j,k),
c      .                fbeta2(i+1,j,k),fbeta2(i-1,j,k),
c      .                fbeta2(i,j+1,k),fbeta2(i,j-1,k),
c      .                fbeta2(i,j,k+1),fbeta2(i,j,k-1))
        alp       = max(ulprec*(1.-bc),0.)
c
        alpr      = (alp*rhoi/c2)*dw1
        fac       = dcp
c       fac       = 0.
        delw(i,j,k,1) = bc*dw1
        delw(i,j,k,2) = dw2 - alpr*uc
        delw(i,j,k,3) = dw3 - alpr*vc
        delw(i,j,k,4) = dw4 - alpr*wc
        delw(i,j,k,5) = dw5 + fac*rhoi*(delw(i,j,k,1) - dw1)
   10 continue
      return
   15 do 20 k = kmn+1,kmx
      do 20 j = jmn+1,jmx
      do 20 i = imn+1,imx
        rhoi      = 1./w(i,j,k,1)
        uc        = w(i,j,k,2)*rhoi
        vc        = w(i,j,k,3)*rhoi
        wc        = w(i,j,k,4)*rhoi
        q2        = 0.5*(uc**2 + vc**2 + wc**2)
        ts        = p(i,j,k)*rhoi
        c2        = gamma*ts
c   
c       transfrom from conservative to (p,u,v,w,T) variables
c
        dw1       = gm1*(q2*delw(i,j,k,1) -uc*delw(i,j,k,2)
     .            -      vc*delw(i,j,k,3) -wc*delw(i,j,k,4)
     .            +         delw(i,j,k,5) )
        dw2       = (delw(i,j,k,2) -uc*delw(i,j,k,1))*rhoi
        dw3       = (delw(i,j,k,3) -vc*delw(i,j,k,1))*rhoi
        dw4       = (delw(i,j,k,4) -wc*delw(i,j,k,1))*rhoi
        dw5       = (dw1 -ts*delw(i,j,k,1))*rhoi
c
c        evaluate alpha and beta for preconditioning
c
        bc        = fbeta2(i,j,k)
c       bc        = max(fbeta2(i,j,k),
c      .                fbeta2(i+1,j,k),fbeta2(i-1,j,k),
c      .                fbeta2(i,j+1,k),fbeta2(i,j-1,k),
c      .                fbeta2(i,j,k+1),fbeta2(i,j,k-1))
        alp       = max(ulprec*(1.-bc),0.)
c
c       am        = q2/c2
        g         = 1. +gm1*delta
        z         = g*dw1 - gamma*w(i,j,k,1)*delta*dw5
        alpr      = (alp*rhoi/c2)*z
        delw(i,j,k,1) = bc*z
        delw(i,j,k,2) = dw2 - alpr*uc
        delw(i,j,k,3) = dw3 - alpr*vc
        delw(i,j,k,4) = dw4 - alpr*wc
        delw(i,j,k,5) = dw5 + (delw(i,j,k,1) - dw1)*rhoi*dcp
   20 continue
      return
c
   30 continue
c                    delta=1
      do 50 k = kmn+1,kmx
      do 50 j = jmn+1,jmx
      do 50 i = imn+1,imx
      rhoi      = 1./w(i,j,k,1)
      uc        = w(i,j,k,2)*rhoi
      vc        = w(i,j,k,3)*rhoi
      wc        = w(i,j,k,4)*rhoi
      us        = uc**2 + vc**2 + wc**2
      ts        = p(i,j,k)  *rhoi
      hc        = cp*ts
      c2        = gamma*ts
c
c      evaluate alpha and beta for preconditioning
c
      bc        = fbeta2(i,j,k)
      if (iprecg.le.0) bc = 1.
      beta2     = bc*c2
c     beta2     = fbeta2(i,j,k)
c     bc        = beta2/c2
      alp       = max(ulprec*(1.-bc),0.)
      alp1      = (1.+alp)*delw(i,j,k,1)
c
      dw1       = beta2*delw(i,j,k,1)
      dw2       = (delw(i,j,k,2) - alp1*uc)*rhoi
      dw3       = (delw(i,j,k,3) - alp1*vc)*rhoi
      dw4       = (delw(i,j,k,4) - alp1*wc)*rhoi
      tp        = dcp*ts/p(i,j,k)
      dw5       =((beta2+.5*us-hc)*delw(i,j,k,1)-uc*delw(i,j,k,2)
     .          - vc*delw(i,j,k,3) - wc*delw(i,j,k,4)
     .          + delw(i,j,k,5) )*tp
c
      delw(i,j,k,1) = dw1
      delw(i,j,k,2) = dw2
      delw(i,j,k,3) = dw3
      delw(i,j,k,4) = dw4
      delw(i,j,k,5) = dw5
  50  continue
      return
c
c                   iprec=1
c      Preconditioning matrix P in conservative variables
c
  100 continue
      do 110 k = kmn+1,kmx
      do 110 j = jmn+1,jmx
      do 110 i = imn+1,imx
        rhoi      = 1./w(i,j,k,1)
        uc        = w(i,j,k,2)*rhoi
        vc        = w(i,j,k,3)*rhoi
        wc        = w(i,j,k,4)*rhoi
        q2        = 0.5*(uc**2 + vc**2 + wc**2)
        ts        = p(i,j,k)*rhoi
        c2        = gamma*ts
c
        hc        = cp*ts + q2
        bc        = fbeta2(i,j,k)
        alp       = max(ulprec*(1.-bc),0.)
c
        z         = gm1*(q2*delw(i,j,k,1) - uc*delw(i,j,k,2)
     .            - vc*delw(i,j,k,3) - wc*delw(i,j,k,4)
     .            + delw(i,j,k,5) )/c2
        dw1       = (1. - bc)*z - delta*(z - delw(i,j,k,1))
        dw2       = (1. + alp - bc)*z
        dw3       = (hc*dw1 + alp*q2)*z
        delw(i,j,k,1) = delw(i,j,k,1) - dw1
        delw(i,j,k,2) = delw(i,j,k,2) - uc*dw2
        delw(i,j,k,3) = delw(i,j,k,3) - vc*dw2
        delw(i,j,k,4) = delw(i,j,k,4) - wc*dw2
        delw(i,j,k,5) = delw(i,j,k,5) - dw3
  110 continue
c
      return
      end
