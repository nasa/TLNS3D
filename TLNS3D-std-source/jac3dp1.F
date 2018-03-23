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
      subroutine jac3dp1 (imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .  ivisb,ws,
     .  w,p,delw,vol,ft,dtphy,ntorder,dtvi,dtvj,dtvk,fbeta2)
c
c        combined Jacobi and low speed preconditioning with iprec=1
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
      dimension    w(imp2,jmp2,kmp2,5),    ws(imp2,jmp2,kmp2,5),
     .             delw(imp2,jmp2,kmp2,5),  p(imp2,jmp2,kmp2),
     .             dtvi(imp2,jmp2,kmp2), dtvj(imp2,jmp2,kmp2),
     .             dtvk(imp2,jmp2,kmp2),  vol(imp2,jmp2,kmp2)
c
      dimension    six(imp1,jmp1,kmp1),
     .             siy(imp1,jmp1,kmp1), siz(imp1,jmp1,kmp1),
     .             sjx(imp1,jmp1,kmp1), sjy(imp1,jmp1,kmp1),
     .             sjz(imp1,jmp1,kmp1), skx(imp1,jmp1,kmp1),
     .             sky(imp1,jmp1,kmp1), skz(imp1,jmp1,kmp1)
c
      dimension    fbeta2(imp2,jmp2,kmp2)
      dimension    ivisb(3)
      dimension    a(4,4)
c
      vepsnn    = 0.3
      vepsll    = 0.3
      difac     = diffac
      tol       = 1.e-14
      dt        = 0.
      if (ntorder.gt.0) dt = ct(ntorder+1)/dtphy
c
      vismultx = 0.
      vismulty = 0.
      vismultz = 0.
c
      if (ivisb(1).gt.0.) vismultx = 1.
      if (ivisb(2).gt.0.) vismulty = 1.
      if (ivisb(3).gt.0.) vismultz = 1.
c
      do 100 k=kmn+1,kmx
      do 100 j=jmn+1,jmx
      do 100 i=imn+1,imx
c
c         transformation to Q = (dp/rho*beta*c , u,v,w, dp - c^2 d rho)
c         transformation from conservative variables
c
        rhoa     = w(i,j,k,1)
        rhoi     = 1./rhoa
        ua       = w(i,j,k,2)*rhoi
        va       = w(i,j,k,3)*rhoi
        wa       = w(i,j,k,4)*rhoi
        qq       = 0.5*(ua*ua + va*va + wa*wa)
        beta2    = fbeta2(i,j,k)
        c2       = gamma*p(i,j,k)*rhoi
        rhoc     = rhoa*sqrt(beta2*c2)
c
        dp       = gm1*(qq*delw(i,j,k,1) - ua*delw(i,j,k,2)
     .           - va*delw(i,j,k,3) - wa*delw(i,j,k,4) + delw(i,j,k,5))
        t1       = dp/rhoc
        t2       = (delw(i,j,k,2) - ua*delw(i,j,k,1))*rhoi
        t3       = (delw(i,j,k,3) - va*delw(i,j,k,1))*rhoi
        t4       = (delw(i,j,k,4) - wa*delw(i,j,k,1))*rhoi
        t5       = dp - c2*delw(i,j,k,1)
c
c          coefficients of Jacobi matrix in "x" direction
c
        sx       = 0.5*(six(i,j,k) + six(i-1,j,k))
        sy       = 0.5*(siy(i,j,k) + siy(i-1,j,k))
        sz       = 0.5*(siz(i,j,k) + siz(i-1,j,k))
c
c       fix for collapsing faces
c
        stot     = sx*sx + sy*sy + sz*sz
        stotnz   = max(tol,stot)
        sqtot    = sqrt(stotnz)
        dsqtot   = 1./sqtot
        cs2      = stotnz*c2
        cs       = sqrt(cs2)
        betac2   = beta2*cs2
        betac    = sqrt(betac2)
c
        sxh      = sx*dsqtot
        syh      = sy*dsqtot
        szh      = sz*dsqtot
        qs       = sx*ua+sy*va+sz*wa
c
c                  eigenvalues
c                dtlp is lambda_+
c                dtlm is lambda_-
c
        sq       = sqrt(((beta2-1.)*qs)**2 + 4.*betac2)
        bqs      = (beta2 + 1.)*qs
        dtlp     = 0.5*(bqs + sq)
        dtlm     = 0.5*(bqs - sq)

c                 absolute values
        qmax     = 0.5*(abs(bqs) + sq)
        specn    = vepsnn*qmax
        specl    = vepsll*qmax
        adtlp    = max(abs(dtlp),specn)
        adtlm    = max(abs(dtlm),specn)
        aqs      = max(abs(qs),specl)
c
        denom    = 1./ (dtlp - dtlm)
        R1       = ((dtlp-qs)*adtlp - (dtlm-qs)*adtlm)*denom
        R1h      = ((dtlp-qs)*adtlm - (dtlm-qs)*adtlp)*denom
        R2       = betac*(adtlp-adtlm)*denom
        R1u      = R1h - aqs
c
c          only half are calculated by symmetry
c          fifth row&column is zero except for diagonal
c                so it is not calculated
c
c      diagonal comes from identity plus time dependent terms
c
        dtvis    = difac*(vismultx*dtvi(i,j,k) + vismulty*dtvj(i,j,k)
     .           +        vismultz*dtvk(i,j,k))
c
        dtl      = anin*dt*vol(i,j,k)
        dtlv     = dtvis + dtl
c
        a11      = beta2*dtl + dtvis + R1
        a12      =             sxh*R2
        a13      =             syh*R2
        a14      =             szh*R2
        a22      = dtlv + aqs + sxh*sxh*R1u
        a23      =             sxh*syh*R1u
        a24      =             sxh*szh*R1u
        a33      = dtlv + aqs + syh*syh*R1u
        a34      =             syh*szh*R1u
        a44      = dtlv + aqs + szh*szh*R1u
        a55      = dtlv + aqs
c
c                get Jacobi coefficients in "y" direction
c
        sx       = 0.5*(sjx(i,j,k) + sjx(i,j-1,k))
        sy       = 0.5*(sjy(i,j,k) + sjy(i,j-1,k))
        sz       = 0.5*(sjz(i,j,k) + sjz(i,j-1,k))
c
c       fix for collapsing faces
c
        stot     = sx*sx + sy*sy + sz*sz
        stotnz   = max(tol,stot)
        sqtot    = sqrt(stotnz)
        dsqtot   = 1./sqtot
        cs2      = stotnz*c2
        cs       = sqrt(cs2)
        betac2   = beta2*cs2
        betac    = sqrt(betac2)
c
        sxh      = sx*dsqtot
        syh      = sy*dsqtot
        szh      = sz*dsqtot
        qs       = sx*ua+sy*va+sz*wa
c
c                  eigenvalues
c                dtlp is lambda_+
c                dtlm is lambda_-
c
        sq       = sqrt(((beta2-1.)*qs)**2 + 4.*betac2)
        bqs      = (beta2 + 1.)*qs
        dtlp     = 0.5*(bqs + sq)
        dtlm     = 0.5*(bqs - sq)

c                 absolute values
        qmax     = 0.5*(abs(bqs) + sq)
        specn    = vepsnn*qmax
        specl    = vepsll*qmax
        adtlp    = max(abs(dtlp),specn)
        adtlm    = max(abs(dtlm),specn)
        aqs      = max(abs(qs),specl)
c
        denom    = 1./ (dtlp - dtlm)
        R1       = ((dtlp-qs)*adtlp - (dtlm-qs)*adtlm)*denom
        R1h      = ((dtlp-qs)*adtlm - (dtlm-qs)*adtlp)*denom
        R2       = betac*(adtlp-adtlm)*denom
        R1u      = R1h - aqs
c
c          only half are calculated by symmetry
c          fifth row&column is zero except for diagonal
c                so it is not calculated
c
        a11      = a11 + R1
        a12      = a12       + sxh*R2
        a13      = a13       + syh*R2
        a14      = a14       + szh*R2
        a22      = a22 + aqs + sxh*sxh*R1u
        a23      = a23       + sxh*syh*R1u
        a24      = a24       + sxh*szh*R1u
        a33      = a33 + aqs + syh*syh*R1u
        a34      = a34       + syh*szh*R1u
        a44      = a44 + aqs + szh*szh*R1u
        a55      = a55 + aqs
c
c                get Jacobi coefficients in "z" direction
c
        sx       = 0.5*(skx(i,j,k) + skx(i,j,k-1))
        sy       = 0.5*(sky(i,j,k) + sky(i,j,k-1))
        sz       = 0.5*(skz(i,j,k) + skz(i,j,k-1))
c
c       fix for collapsing faces
c
        stot     = sx*sx + sy*sy + sz*sz
        stotnz   = max(tol,stot)
        sqtot    = sqrt(stotnz)
        dsqtot   = 1./sqtot
        cs2      = stotnz*c2
        cs       = sqrt(cs2)
        betac2   = beta2*cs2
        betac    = sqrt(betac2)
c
        sxh      = sx*dsqtot
        syh      = sy*dsqtot
        szh      = sz*dsqtot
        qs       = sx*ua+sy*va+sz*wa
c
c                  eigenvalues
c                dtlp is lambda_+
c                dtlm is lambda_-
c
        sq       = sqrt(((beta2-1.)*qs)**2 + 4.*betac2)
        bqs      = (beta2 + 1.)*qs
        dtlp     = 0.5*(bqs + sq)
        dtlm     = 0.5*(bqs - sq)

c                 absolute values
        qmax     = 0.5*(abs(bqs) + sq)
        specn    = vepsnn*qmax
        specl    = vepsll*qmax
        adtlp    = max(abs(dtlp),specn)
        adtlm    = max(abs(dtlm),specn)
        aqs      = max(abs(qs),specl)
c
        denom    = 1./ (dtlp - dtlm)
        R1       = ((dtlp-qs)*adtlp - (dtlm-qs)*adtlm)*denom
        R1h      = ((dtlp-qs)*adtlm - (dtlm-qs)*adtlp)*denom
        R2       = betac*(adtlp-adtlm)*denom
        R1u      = R1h - aqs
c
c          only half are calculated by symmetry
c          fifth row&column is zero except for diagonal
c                so it is not calculated
c
        a11      = a11 + R1
        a12      = a12       + sxh*R2
        a13      = a13       + syh*R2
        a14      = a14       + szh*R2
        a22      = a22 + aqs + sxh*sxh*R1u
        a23      = a23       + sxh*syh*R1u
        a24      = a24       + sxh*szh*R1u
        a33      = a33 + aqs + syh*syh*R1u
        a34      = a34       + syh*szh*R1u
        a44      = a44 + aqs + szh*szh*R1u
        a55      = a55 + aqs
c
c           Multiply RHS by P
c
       t1  = beta2*t1
c
c           solve C*x = T
c
c                 get L and U factors
c                 store in same place
c
clu     a(1,1)   = a11
clu     a(1,2)   = a12
clu     a(1,3)   = a13
clu     a(1,4)   = a14
clu     a(2,1)   = a12/a11
clu     a(2,2)   = a22 - a(2,1)*a12
clu     a(2,3)   = a23 - a(2,1)*a13
clu     a(2,4)   = a24 - a(2,1)*a13
clu     a(3,1)   = a13/a11
clu     a(3,2)   = (a23-a(3,1)*a12)/a(2,2)
clu     a(3,3)   = a33 - a(3,1)*a13 - a(3,2)*a(2,3)
c
c                 get U = U^t factors
c                 store in place
c
        a(1,1)   = sqrt(a11)
        a(1,2)   = a12/a(1,1)
        a(1,3)   = a13/a(1,1)
        a(1,4)   = a14/a(1,1)
        a(2,2)   = sqrt(a22 - a(1,2)*a(1,2))
        a(2,3)   = (a23 - a(1,2)*a(1,3))/a(2,2)
        a(2,4)   = (a24 - a(1,2)*a(1,4))/a(2,2)
        a(3,3)   = sqrt(a33 - a(1,3)*a(1,3) - a(2,3)*a(2,3))
        a(3,4)   = (a34 - a(1,3)*a(1,4) - a(2,3)*a(2,4))/a(3,3)
        a(4,4)   = a44 - a(1,4)*a(1,4) - a(2,4)*a(2,4) - a(3,4)*a(3,4)
c
c                    forward substitution
c
clu     t2       = t2 - a(2,1)*t1
clu     t3       = (t3 - a(3,1)*t1 - a(3,2)*t2)/a(3,3)
clu     t4       = (t4 - a(4,1)*t1 - a(4,2)*t2 - a(4,3)*t3)/a(4,4)
  
        t1       = t1/a(1,1)
        t2       = (t2 - a(1,2)*t1)/a(2,2)
        t3       = (t3 - a(1,3)*t1 - a(2,3)*t2)/a(3,3)
        t4       = (t4 - a(1,4)*t1 - a(2,4)*t2 - a(3,4)*t3)/a(4,4)
c 
c                    backward substitution
c
        t3       = (t3 - a(3,4)*t4)/a(3,3)
        t2       = (t2 - a(2,3)*t3 - a(2,4)*t4)/a(2,2)
        t1       = (t1 - a(1,2)*t2 - a(1,3)*t3 - a(1,4)*t4)/a(1,1)

        t5       = t5/a55
c
c
c          transform RHS from Q= (dp/rhoc, du,dv, dS) variables
c                     to w=conservative variables
c                  i.e. multiply delw by dw/dQ
c
        bt1              = (rhoc*t1 - t5)/c2
        t5               = qq*bt1 + rhoc*t1*dgm1
     .                   + w(i,j,k,2)*t2
     .                   + w(i,j,k,3)*t3

        t1               = bt1
        t2               = rhoa*t2 + ua*bt1
        t3               = rhoa*t3 + va*bt1
        t4               = rhoa*t4 + wa*bt1

c       tp               = gm1*(qq*t1 - ua*t2 - va*t3)/gamma
c       t5               = t5 - tp
c
        delw(i,j,k,1)    = ft*t1
        delw(i,j,k,2)    = ft*t2
        delw(i,j,k,3)    = ft*t3
        delw(i,j,k,4)    = ft*t4
        delw(i,j,k,5)    = ft*t5
  100 continue
      return
      end
