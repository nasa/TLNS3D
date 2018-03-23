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
      subroutine  filtfz(ncyc,nstage,imn,jmn,kmn,imx,jmx,kmx,
     .                   imp1,jmp1,kmp1,imp2,jmp2,kmp2,w,p,vol,
     .                   six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                   dtl,dti,dtj,dtk,ratioij,ratiojk,
     .                   dw,ew,fs,df2,df4,dp,ep,phi,fw,fv,
     .                   mx1dwk,wk3d,
     .                   dtvi,dtvj,dtvk,fbeta2)
c=======================================================================
c 
c           artificial dissipations (for fine mesh of mg)
c           blended second and fourth diferrences
c
c     latest mods    : by v.n.vatsa for preconditioning (Feb. 1998)
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
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubfil2,iprec
c
      common/vis   / diffac,rsvis,omega
c
      dimension    fw(imp2,jmp2,kmp2,5),
     .             fv(imp2,jmp2,kmp2  ),
     .              w(imp2,jmp2,kmp2,5), 
     .              p(imp2,jmp2,kmp2),vol(imp2,jmp2,kmp2),
     .            dtl(imp2,jmp2,kmp2),dti(imp2,jmp2,kmp2),
     .            dtj(imp2,jmp2,kmp2),dtk(imp2,jmp2,kmp2),
     .            ratioij(imp2,jmp2,kmp2),ratiojk(imp2,jmp2,kmp2),
     .       dtvi(imp2,jmp2,kmp2),dtvj(imp2,jmp2,kmp2),
     .       dtvk(imp2,jmp2,kmp2),
     .     fbeta2(imp2,jmp2,kmp2)
c
      dimension   six(imp1,jmp1,kmp1),siy(imp1,jmp1,kmp1),
     .            siz(imp1,jmp1,kmp1),sjx(imp1,jmp1,kmp1),
     .            sjy(imp1,jmp1,kmp1),sjz(imp1,jmp1,kmp1),
     .            skx(imp1,jmp1,kmp1),sky(imp1,jmp1,kmp1),
     .            skz(imp1,jmp1,kmp1)
c
c     space allocation for working arrays
c
      dimension  dw(mx1dwk,5),ew(mx1dwk,5),fs(mx1dwk,5),df2(mx1dwk,5),
     .          df4(mx1dwk,5),dp(mx1dwk  ),ep(mx1dwk  ),phi(mx1dwk  )
      dimension   wk3d(imp2,jmp2,kmp2,5)
c
c
      tol       = 1.e-14
      fis2      = abs(vis2) 
      fis4      = vis4/32.
      bet       = beta(nstage)
      cp        = gamma*dgm1
      dcp       = 1./cp
      delta     = deltpre
c     delta     = 0.
c
      fill      = 1.
      if (qfil(1).lt.0.) fill = 0.
c
c--------- dissipation in k direction 
c
Cmic$ Doall Autoscope
      do 299 j=jmn+1,jmx
      do 299 i=imn+1,imx
c
      do k=kmn-1,kmp2
        dp(k)     = fv(i,j,k)
      enddo
c
c         variables depending on iprec (see subroutine getdifk)
c         iprec=0,1  wk3d is first difference of conservation variables
c         iprec=2,3  wk3d is first difference of (p,u,v,w,T) variables
c
      do k=kmn-1,kmp1
        dw(k,1)   = wk3d(i,j,k,1)
        dw(k,2)   = wk3d(i,j,k,2)
        dw(k,3)   = wk3d(i,j,k,3)
        dw(k,4)   = wk3d(i,j,k,4)
        dw(k,5)   = wk3d(i,j,k,5)
      enddo
c
      do k=kmn,kmx
        ep(k)     = max(dp(k-1),dp(k),dp(k+1),dp(k+2))
      enddo
c
      do k=kmn,kmp1
        reigjk    = ratiojk(i,j,k)                  
        reigik    = ratioij(i,j,k)*ratiojk(i,j,k)
        phi(k)    =  1. + reigik + reigjk
      enddo
c 
      do k=kmn,kmx
        ew(k,1)   = dw(k+1,1)-2.*dw(k,1)+dw(k-1,1)
        ew(k,2)   = dw(k+1,2)-2.*dw(k,2)+dw(k-1,2)
        ew(k,3)   = dw(k+1,3)-2.*dw(k,3)+dw(k-1,3)
        ew(k,4)   = dw(k+1,4)-2.*dw(k,4)+dw(k-1,4)
        ew(k,5)   = dw(k+1,5)-2.*dw(k,5)+dw(k-1,5)
      enddo
c
      if (vepsn .ge. 1.0) then
c                  SCALAR VISCOSITY
        do k=kmn,kmx
        fila      = 0.5*((dtk(i,j,k+1)+diffac*dtvk(i,j,k+1))*phi(k+1)
     .                 + (dtk(i,j,k  )+diffac*dtvk(i,j,k  ))*phi(k  ))
        dis2      = fis2*ep(k)
        dis4      = dim(fis4,dis2)*fila
        dis2      = dis2*fila
        fs(k,1)   = dis2*dw(k,1) - dis4*ew(k,1)
        fs(k,2)   = dis2*dw(k,2) - dis4*ew(k,2)
        fs(k,3)   = dis2*dw(k,3) - dis4*ew(k,3)
        fs(k,4)   = dis2*dw(k,4) - dis4*ew(k,4)
        fs(k,5)   = dis2*dw(k,5) - dis4*ew(k,5)
        enddo
c
c                 MATRIX VISCOSITY
c
      else if (iprec.le.0) then
c                 no preconditioning
c                 standard (E) matrix viscosity
        do 230 k=kmn,kmx
        rhoa    = 0.5*(w(i,j,k+1,1)+w(i,j,k,1))
        rhoi    = 1./rhoa
        uu      = 0.5*(w(i,j,k+1,2)/w(i,j,k+1,1)+w(i,j,k,2)/w(i,j,k,1))
        vv      = 0.5*(w(i,j,k+1,3)/w(i,j,k+1,1)+w(i,j,k,3)/w(i,j,k,1))
        ww      = 0.5*(w(i,j,k+1,4)/w(i,j,k+1,1)+w(i,j,k,4)/w(i,j,k,1))
        vels    = 0.5*(uu*uu+vv*vv+ww*ww)
        pa      = 0.5*(p(i,j,k+1) + p(i,j,k))
        c2      = gamma*pa*rhoi
        he      = c2*dgm1+vels
c              Roe average
c       srhol   = sqrt(w(i,j,k,1))
c       srhor   = sqrt(w(i,j,k+1,1))
c       anorm   = 1./(srhol+srhor)
c       rhoa    = srhol*srhor
c       rhoi    = 1./rhoa
c       uu      = (srhol*w(i,j,k,2)/w(i,j,k,1) +
c    .             srhor*w(i,j,k+1,2)/w(i,j,k+1,1))*anorm
c       vv      = (srhol*w(i,j,k,3)/w(i,j,k,1) +
c    .             srhor*w(i,j,k+1,3)/w(i,j,k+1,1))*anorm
c       ww      = (srhol*w(i,j,k,4)/w(i,j,k,1) +
c    .             shror*w(i,j,k+1,4)/w(i,j,k+1,1))*anorm
c       vels    = .5*(uu*uu+vv+vv+ww*ww)
c       pa      = (srhol*p(i,j,k) + srhor*p(i,j,k+1))*anorm
c
        sx      = skx(i,j,k)
        sy      = sky(i,j,k)
        sz      = skz(i,j,k)
c
c       fix for collapsing faces
c
        stot    = sx*sx+sy*sy+sz*sz
        stotnz  = max(tol,stot)
        srat    = stot/stotnz
        sqtrat  = sqrt(srat)
c
        qs      = sx*uu+sy*vv+sz*ww
        cc      = sqrt(c2)
        cs      = sqrt(stotnz*c2)
c       cs      = sqrt(stot*c2)
c
c                  eigenvalues
c                dtlp is lambda_+
c                dtlm is lambda_-
c
        dtlp    = qs + cs*sqtrat
        dtlm    = qs - cs*sqtrat
c       dtlp    = qs + cs
c       dtlm    = qs - cs
c
        specn   = vepsn*(abs(qs)+cs*sqtrat)
        specl   = vepsl*(abs(qs)+cs*sqtrat)
c       specn   = vepsn*(abs(qs)+cs)
c       specl   = vepsl*(abs(qs)+cs)
c
        e1      = max(abs(dtlp),specn)
        e2      = max(abs(dtlm),specn)
        e3      = max(abs(qs),specl)
c       e1      = max(abs(dtlp),specn) + diffac*dtvk(i,j,k)
c       e2      = max(abs(dtlm),specn) + diffac*dtvk(i,j,k)
c       e3      = max(abs(qs),specl)   + diffac*dtvk(i,j,k)
c
c                 JST switch
c
        phia    = 0.5*(phi(k+1) + phi(k))
        dis2    = fis2*ep(k)
        dis4    = dim(fis4,dis2)*phia
        dis2    = dis2*phia
        q1      = dis2*dw(k,1) - dis4*ew(k,1)
        q2      = dis2*dw(k,2) - dis4*ew(k,2)
        q3      = dis2*dw(k,3) - dis4*ew(k,3)
        q4      = dis2*dw(k,4) - dis4*ew(k,4)
        q5      = dis2*dw(k,5) - dis4*ew(k,5)
c
c                abs(C)
c
        s1      = 0.5*(e1+e2)
        s2      = 0.5*(e1-e2)
        s3      = e3
        ak1     = (s1-s3)*gm1/c2
c
c       fix for collapsing area
c
c       ak2     = s2/cs
c       abw     = c2*dgm1/(sx*sx+sy*sy+sz*sz)
        ak2     = s2/cs*sqtrat
        abw     = c2*dgm1/stotnz*srat
c
        aa      = sx*abw
        bb      = sy*abw
        cc      = sz*abw
        wc      = qs*abw
c
        vd      = vels*q1-uu*q2-vv*q3-ww*q4+q5
        qd      = qs*q1-sx*q2-sy*q3-sz*q4
        fs(k,1) = ak1*vd-ak2*qd+s3*q1
        fs(k,2) = ak1*(uu*vd-aa*qd)+ak2*(gm1*sx*vd-uu*qd)+s3*q2
        fs(k,3) = ak1*(vv*vd-bb*qd)+ak2*(gm1*sy*vd-vv*qd)+s3*q3
        fs(k,4) = ak1*(ww*vd-cc*qd)+ak2*(gm1*sz*vd-ww*qd)+s3*q4
        fs(k,5) = ak1*(he*vd-wc*qd)+ak2*(gm1*qs*vd-he*qd)+s3*q5
  230   continue
      elseif (iprec.eq.1) then
c           matrix + preconditioning conservation variables
c           assumes delta = 0
        delta     = deltpre
        do 240 k=kmn,kmx
        rhoa    = 0.5*(w(i,j,k+1,1)+w(i,j,k,1))
        rhoi    = 1./rhoa
        uu      = 0.5*(w(i,j,k+1,2)/w(i,j,k+1,1)+w(i,j,k,2)/w(i,j,k,1))
        vv      = 0.5*(w(i,j,k+1,3)/w(i,j,k+1,1)+w(i,j,k,3)/w(i,j,k,1))
        ww      = 0.5*(w(i,j,k+1,4)/w(i,j,k+1,1)+w(i,j,k,4)/w(i,j,k,1))
        vels    = 0.5*(uu*uu+vv*vv+ww*ww)
        pa      = 0.5*(p(i,j,k+1) + p(i,j,k))
c              Roe average
c       srhol   = sqrt(w(i,j,k,1))
c       srhor   = sqrt(w(i,j,k+1,1))
c       anorm   = 1./(srhol+srhor)
c       rhoa    = srhol*srhor
c       rhoi    = 1./rhoa
c       uu      = (srhol*w(i,j,k,2)/w(i,j,k,1) +
c    .             srhor*w(i,j,k+1,2)/w(i,j,k+1,1))*anorm
c       vv      = (srhol*w(i,j,k,3)/w(i,j,k,1) +
c    .             srhor*w(i,j,k+1,3)/w(i,j,k+1,1))*anorm
c       ww      = (srhol*w(i,j,k,4)/w(i,j,k,1) +
c    .             shror*w(i,j,k+1,4)/w(i,j,k+1,1))*anorm
c       vels    = .5*(uu*uu+vv+vv+ww*ww)
c       pa      = (srhol*p(i,j,k) + srhor*p(i,j,k+1)*anorm
c
        c2      = gamma*pa*rhoi
        cc      = sqrt(c2)
c
        sx      = skx(i,j,k)
        sy      = sky(i,j,k)
        sz      = skz(i,j,k)
c
c       fix for collapsing faces
c
        stot    = sx*sx+sy*sy+sz*sz
        stotnz  = max(tol,stot)
        sqtot   = sqrt(stotnz)
        dsqtot  = 1./sqtot
c       srat    = stot/stotnz
c       sqtrat  = sqrt(srat)
c
        sxh     = sx*dsqtot
        syh     = sy*dsqtot
        szh     = sz*dsqtot
        qs      = sx*uu+sy*vv+sz*ww
        qsh     = sxh*uu+syh*vv+szh*ww
c
c       bc      = 0.5*(fbeta2(i,j,k) + fbeta2(i,j,k+1))
        bc      = max (fbeta2(i,j,k) , fbeta2(i,j,k+1))
        beta2   = bc*c2
        alp     = max(ulprec*(1.-bc),0.)
        cs1     = 0.5*(1. - alp + bc)
        csq     = cs1*qs
        csqh    = cs1*qsh
        cs      = sqrt(csq**2  + stotnz*beta2 - qs*qs*bc)
        csh     = sqrt(csqh**2 +        beta2 - qsh*qsh*bc)
c
c                  eigenvalues
c                dtlp is lambda_+
c                dtlm is lambda_-
c
c       dtlp    = csq  + cs*sqtrat
c       dtlph   = csqh + csh*sqtrat
c       dtlm    = csq  - cs*sqtrat
c       dtlmh   = csqh - csh*sqtrat
        dtlp    = csq  + cs
        dtlph   = csqh + csh
        dtlm    = csq  - cs
        dtlmh   = csqh - csh
c
c       specn   = vepsn*(abs(csq)+cs*sqtrat)
c       specl   = vepsl*(abs(csq)+cs*sqtrat)
        specn   = vepsn*(abs(csq)+cs)
        specl   = vepsl*(abs(csq)+cs)
        e1      = max(abs(dtlp),specn)
        e2      = max(abs(dtlm),specn)
        e3      = max(abs(qs),specl)
c       e1      = max(abs(dtlp),specn) + diffac*dtvk(i,j,k)
c       e2      = max(abs(dtlm),specn) + diffac*dtvk(i,j,k)
c       e3      = max(abs(qs),specl)   + diffac*dtvk(i,j,k)
c
c                 JST switch
c
        phia    = 0.5*(phi(k+1) + phi(k))
        dis2    = fis2*ep(k)
        dis4    = dim(fis4,dis2)*phia
        dis2    = dis2*phia
        qq1     = dis2*dw(k,1) - dis4*ew(k,1)
        qq2     = dis2*dw(k,2) - dis4*ew(k,2)
        qq3     = dis2*dw(k,3) - dis4*ew(k,3)
        qq4     = dis2*dw(k,4) - dis4*ew(k,4)
        qq5     = dis2*dw(k,5) - dis4*ew(k,5)
c
c         multiply by Jacobian of w_0 with respect to w_1
c
        z       = gm1*(vels*qq1 - uu*qq2 - vv*qq3 -ww*qq4 + qq5)
        q1      = z
        q2      = (qq2 - uu*qq1)*rhoi
        q3      = (qq3 - vv*qq1)*rhoi
        q4      = (qq4 - ww*qq1)*rhoi
        q5      = z - c2*qq1
c
c                    abs(C)
c
c         R_{ij} are the elements of R^{-1}
c         the matrix whose colmuns are the eigenvectors
c
        R21     = (beta2*sxh-alp*uu*dtlph)/(dtlph-qsh)
        R31     = (beta2*syh-alp*vv*dtlph)/(dtlph-qsh)
        R41     = (beta2*szh-alp*ww*dtlph)/(dtlph-qsh)
        R22     = (beta2*sxh-alp*uu*dtlmh)/(dtlmh-qsh)
        R32     = (beta2*syh-alp*vv*dtlmh)/(dtlmh-qsh)
        R42     = (beta2*szh-alp*ww*dtlmh)/(dtlmh-qsh)
c
c          S_{ij} are the elements of R
c
        det     = 1./(sxh*(R21-R22) + syh*(R31-R32) + szh*(R41-R42))
        yy      = (sxh*q2 + syh*q3 + szh*q4)*det
        rb      = rhoa*beta2
        S11     = - (sxh*R22 + syh*R32 + szh*R42)*det/rb
c       S21     = 1./rb - S11
c
        z1      = S11*q1 + yy
c       z2      = S21*q1 - yy
        z2      = q1/rb - z1
        e13     = (e1-e3)*z1
        e23     = (e2-e3)*z2
c
c       multiply by diagonal eigenvalues
c       then multiply by R^{-1} and finally by inverse Jacobian
c
        fs1     = rb*(e13 + e23)    + e3*q1
        fs2     = R21*e13 + R22*e23 + e3*q2
        fs3     = R31*e13 + R32*e23 + e3*q3
        fs4     = R41*e13 + R42*e23 + e3*q4
        fs5     = e3*q5
c
c         multiply by Jacobian of w_1 with respect to w_0
c
        z       = (fs1 - fs5)/c2
        fs(k,1) = z
        fs(k,2) = rhoa*fs2 + uu*z
        fs(k,3) = rhoa*fs3 + vv*z
        fs(k,4) = rhoa*fs4 + ww*z
        fs(k,5) = vels*z + fs1*dgm1
     .          + rhoa*(uu*fs2 + vv*fs3 + ww*fs4)
c 240   continue
c            for conservative formulation multiply by P^{-1}
c     do 250 k=kmn,kmx
c       rhoi      = 1./w(i,j,k,1)
c       uu        = w(i,j,k,2)*rhoi
c       vv        = w(i,j,k,3)*rhoi
c       ww        = w(i,j,k,4)*rhoi
c       vels      = 0.5*(uc**2 + vv**2 + ww**2)
c       ts        = p(i,j,k)*rhoi
c       c2        = gamma*ts
        hc        = cp*pa*rhoi + vels
c       bc        = fbeta2(i,j,k)
        z         = gm1*(vels*fs(k,1) - uu*fs(k,2)
     .            - vv*fs(k,3) - ww*fs(k,4)
     .            + fs(k,5) )/c2
        dw1       = ((1. - bc)*z - delta*(z - fs(k,1)))/bc
        fs(k,1)   = fs(k,1) + dw1
        fs(k,2)   = fs(k,2) + uu*dw1
        fs(k,3)   = fs(k,3) + vv*dw1
        fs(k,4)   = fs(k,4) + ww*dw1
        fs(k,5)   = fs(k,5) + hc*dw1
c 250 continue
  240 continue
c
      elseif (iprec.ge.2) then
c           matrix + preconditioning (p,u,v,w,T) variables
c           assumes delta = 0
c
        do 260 k=kmn,kmx
        rhoa    = 0.5*(w(i,j,k+1,1)+w(i,j,k,1))
        rhoi    = 1./rhoa
        uu      = 0.5*(w(i,j,k+1,2)/w(i,j,k+1,1)+w(i,j,k,2)/w(i,j,k,1))
        vv      = 0.5*(w(i,j,k+1,3)/w(i,j,k+1,1)+w(i,j,k,3)/w(i,j,k,1))
        ww      = 0.5*(w(i,j,k+1,4)/w(i,j,k+1,1)+w(i,j,k,4)/w(i,j,k,1))
        vels    = 0.5*(uu*uu+vv*vv+ww*ww)
        pa      = 0.5*(p(i,j,k+1) + p(i,j,k))
c              Roe average
c       srhol   = sqrt(w(i,j,k,1))
c       srhor   = sqrt(w(i,j,k+1,1))
c       anorm   = 1./(srhol+srhor)
c       rho     = srhol*srhor
c       rhoi    = 1./rho
c       uu      = (srhol*w(i,j,k,2)/w(i,j,k,1) +
c    .             srhor*w(i,j,k+1,2)/w(i,j,k+1,1))*anorm
c       vv      = (srhol*w(i,j,k,3)/w(i,j,k,1) +
c    .             srhor*w(i,j,k+1,3)/w(i,j,k+1,1))*anorm
c       ww      = (srhol*w(i,j,k,4)/w(i,j,k,1) +
c    .             shror*w(i,j,k+1,4)/w(i,j,k+1,1))*anorm
c       vels    = .5*(uu*uu+vv+vv+ww*ww)
c       pa      = (srhol*p(i,j,k) + srhor*p(i,j,k+1)*anorm
c
        c2      = gamma*pa*rhoi
        cc      = sqrt(c2)
c
        sx      = skx(i,j,k)
        sy      = sky(i,j,k)
        sz      = skz(i,j,k)
c
c       fix for collapsing faces
c
        stot    = sx*sx+sy*sy+sz*sz
        stotnz  = max(tol,stot)
        sqtot   = sqrt(stotnz)
        dsqtot  = 1./sqtot
c       srat    = stot/stotnz
c       sqtrat  = sqrt(srat)
c
        sxh     = sx*dsqtot
        syh     = sy*dsqtot
        szh     = sz*dsqtot
        qs      = sx*uu+sy*vv+sz*ww
        qsh     = sxh*uu+syh*vv+szh*ww
c
c       bc      = 0.5*(fbeta2(i,j,k) + fbeta2(i,j,k+1))
        bc      = max (fbeta2(i,j,k) , fbeta2(i,j,k+1))
        beta2   = bc*c2
        alp     = max(ulprec*(1.-bc),0.)
        cs1     = 0.5*(1. - alp + bc)
        csq     = cs1*qs
        csqh    = cs1*qsh
        cs      = sqrt(csq**2  + stotnz*beta2 - qs*qs*bc)
        csh     = sqrt(csqh**2 +        beta2 - qsh*qsh*bc)
c
c                  eigenvalues
c                dtlp is lambda_+
c                dtlm is lambda_-
c
c       dtlp    = csq  + cs*sqtrat
c       dtlph   = csqh + csh*sqtrat
c       dtlm    = csq  - cs*sqtrat
c       dtlmh   = csqh - csh*sqtrat
        dtlp    = csq  + cs
        dtlph   = csqh + csh
        dtlm    = csq  - cs
        dtlmh   = csqh - csh
c
c       specn   = vepsn*(abs(csq)+cs*sqtrat)
c       specl   = vepsl*(abs(csq)+cs*sqtrat)
        specn   = vepsn*(abs(csq)+cs)
        specl   = vepsl*(abs(csq)+cs)
        e1      = max(abs(dtlp),specn)
        e2      = max(abs(dtlm),specn)
        e3      = max(abs(qs),specl)
c       e1      = max(abs(dtlp),specn) + diffac*dtvk(i,j,k)
c       e2      = max(abs(dtlm),specn) + diffac*dtvk(i,j,k)
c       e3      = max(abs(qs),specl)   + diffac*dtvk(i,j,k)
c
c                 JST switch
c
        phia    = 0.5*(phi(k+1) + phi(k))
        dis2    = fis2*ep(k)
        dis4    = dim(fis4,dis2)*phia
        dis2    = dis2*phia
c
        q1      = dis2*dw(k,1) - dis4*ew(k,1)
        q2      = dis2*dw(k,2) - dis4*ew(k,2)
        q3      = dis2*dw(k,3) - dis4*ew(k,3)
        q4      = dis2*dw(k,4) - dis4*ew(k,4)
c         multiply by Jacobian of w_0 with respect to w_1
        qq5     = dis2*dw(k,5) - dis4*ew(k,5)
        q5      = gamma*qq5 - gm1*rhoi*q1
c
c                abs(C)
c
c         R_{ij} are the elements of R^{-1}
c         the matrix whose colmuns are the eigenvectors
c
        R21     = (beta2*sxh-alp*uu*dtlph)/(dtlph-qsh)
        R31     = (beta2*syh-alp*vv*dtlph)/(dtlph-qsh)
        R41     = (beta2*szh-alp*ww*dtlph)/(dtlph-qsh)
        R22     = (beta2*sxh-alp*uu*dtlmh)/(dtlmh-qsh)
        R32     = (beta2*syh-alp*vv*dtlmh)/(dtlmh-qsh)
        R42     = (beta2*szh-alp*ww*dtlmh)/(dtlmh-qsh)
c
c          S_{ij} are the elements of R
c
        det     = 1./(sxh*(R21-R22) + syh*(R31-R32) + szh*(R41-R42))
        yy      = (sxh*q2 + syh*q3 + szh*q4)*det
        rb      = rhoa*beta2
        S11     = - (sxh*R22 + syh*R32 + szh*R42)*det/rb
c       S21     = 1./rb - S11
c
        z1      = S11*q1 + yy
c       z2      = S21*q1 - yy
        z2      = q1/rb - z1
        e13     = (e1-e3)*z1
        e23     = (e2-e3)*z2
c
c         multiply by Jacobian of w_1 with respect to w_0
c
        fs(k,1) = rb*(e1*z1+e2*z2)
        fs(k,2) = R21*e13 + R22*e23 + e3*q2
        fs(k,3) = R31*e13 + R32*e23 + e3*q3
        fs(k,4) = R41*e13 + R42*e23 + e3*q4
        fs5     =                     e3*q5
        fs(k,5) = dcp*rhoi*fs(k,1) + fs5/gamma
  260 continue
c
      if (iprec.eq.2) then
c         for conservative formulation multiply by Gamma^{-1}
c
      if (delta.eq.0.) then
      do 270 k=kmn,kmx
        rhoa    = 0.5*(w(i,j,k,1)+w(i,j,k+1,1))
        rhoi    = 1./rhoa
        rhou    = 0.5*(w(i,j,k,2) + w(i,j,k+1,2))
        rhov    = 0.5*(w(i,j,k,3) + w(i,j,k+1,3))
        rhow    = 0.5*(w(i,j,k,4) + w(i,j,k+1,4))
        ts      = 0.5*(p(i,j,k)   + p(i,j,k+1))*rhoi
        hqs     = 0.5*(rhou*rhou+rhov*rhov+rhow*rhow)*rhoi*rhoi
        c2      = gamma*ts
        dc2     = 1./c2
        hh      = c2*dgm1 + hqs
c
c               evaluate alpha and beta for preconditioning
c
c       bc      = 0.5*(fbeta2(i,j,k) + fbeta2(i,j,k+1))
        bc      = max (fbeta2(i,j,k) , fbeta2(i,j,k+1))
        alp     = max(ulprec*(1.-bc),0.)
        dbeta2  = 1./bc
        theta   = (gm1 + dbeta2)*dc2
        rtq5    = rhoa*fs(k,5)/ts
        aq1     = alp*fs(k,1)*dbeta2*dc2
c
        dw1     = theta*fs(k,1) - rtq5
        dwpt    = rhoi*(dw1 + aq1)
        dw2     = rhoa*fs(k,2) + rhou*dwpt
        dw3     = rhoa*fs(k,3) + rhov*dwpt
        dw4     = rhoa*fs(k,4) + rhow*dwpt
        dw5     = (theta*hh-1.)*fs(k,1)
     .          + rhou*fs(k,2) + rhov*fs(k,3) + rhow*fs(k,4)
     .          + hqs*(2.*aq1 - rtq5)
c
        fs(k,1) = dw1
        fs(k,2) = dw2
        fs(k,3) = dw3
        fs(k,4) = dw4
        fs(k,5) = dw5
  270 continue
      else
c             delta=1
      do 280 k=kmn,kmx
        rhoa    = 0.5*(w(i,j,k,1)+w(i,j,k+1,1))
        rhoi    = 1./rhoa
        uu      = 0.5*(w(i,j,k,2) + w(i,j,k+1,2))*rhoi
        vv      = 0.5*(w(i,j,k,3) + w(i,j,k+1,3))*rhoi
        ww      = 0.5*(w(i,j,k,4) + w(i,j,k+1,4))*rhoi
        ts      = 0.5*(p(i,j,k)   + p(i,j,k+1))*rhoi
        hh      = 0.5*(w(i,j,k,5) + w(i,j,k+1,5))*rhoi + ts
        c2      = gamma*ts
        us      = uu*uu + vv*vv + ww*ww
c
c               evaluate alpha and beta for preconditioning
c
c       bc      = 0.5*(fbeta2(i,j,k) + fbeta2(i,j,k+1))
        bc      = max (fbeta2(i,j,k) , fbeta2(i,j,k+1))
        beta2   = bc*c2
        alp     = max(ulprec*(1.-bc),0.)
        dbeta2  = 1./beta2
        alp1    = (1.+alp)*dbeta2*fs(k,1)
c
        dw1     = dbeta2*fs(k,1)
        dw2     = rhoa*fs(k,2) + alp1*uu
        dw3     = rhoa*fs(k,3) + alp1*vv
        dw4     = rhoa*fs(k,4) + alp1*ww
        g41     = dbeta2*(hh + alp*us) -1.
        dw5     = g41*fs(k,1)
     .          + rhoa*(uu*fs(k,2) + vv*fs(k,3) + ww*fs(k,4)
     .          + cp*fs(k,5))
        fs(k,1) = dw1
        fs(k,2) = dw2
        fs(k,3) = dw3
        fs(k,4) = dw4
        fs(k,5) = dw5
  280 continue
      endif ! delta
      endif ! iprec=2
      endif ! iprec .ge. 2
c
      do k=kmn+1,kmx
        fw(i,j,k,1) = fw(i,j,k,1)  + bet*(fs(k-1,1)-fs(k,1))
        fw(i,j,k,2) = fw(i,j,k,2)  + bet*(fs(k-1,2)-fs(k,2))
        fw(i,j,k,3) = fw(i,j,k,3)  + bet*(fs(k-1,3)-fs(k,3))
        fw(i,j,k,4) = fw(i,j,k,4)  + bet*(fs(k-1,4)-fs(k,4))
        fw(i,j,k,5) = fw(i,j,k,5)  + bet*(fs(k-1,5)-fs(k,5))
      enddo
  299 continue
      return
      end
