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
      subroutine  step(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                 imp2,jmp2,kmp2,ibloc,w,p,vol,x,
     .                 six,siy,siz,sjx,sjy,sjz,skx,sky,skz,eomu,
     .                 dti,dtj,dtk,dtmin,ivisb,amu,etemp,mx1dwk,
     .                 ntorder, dtphy,
     .                 dtvi,dtvj,dtvk,fbeta2,fbetav2 )
c=======================================================================
c
c     object :  compute local time step and eigen-values for
c               dissipation and resudual smoothing routines
c
c     initial coding : by V.N.Vatsa of NASA Langley (Feb. 1990)
c     latest mods    : July 1995, to accommodate preconditioning
c
c
      dimension w(imp2,jmp2,kmp2,5),p(imp2,jmp2,kmp2),
     .        vol(imp2,jmp2,kmp2),dti(imp2,jmp2,kmp2),
     .        dtj(imp2,jmp2,kmp2),dtk(imp2,jmp2,kmp2),
     .       eomu(imp2,jmp2,kmp2),dtvi(imp2,jmp2,kmp2),
     .       dtvj(imp2,jmp2,kmp2),dtvk(imp2,jmp2,kmp2),
     .     fbeta2(imp2,jmp2,kmp2),fbetav2(imp2,jmp2,kmp2)
c
      dimension                  x(imp1,jmp1,kmp1,3),
     .     six(imp1,jmp1,kmp1),siy(imp1,jmp1,kmp1),siz(imp1,jmp1,kmp1),
     .     sjx(imp1,jmp1,kmp1),sjy(imp1,jmp1,kmp1),sjz(imp1,jmp1,kmp1),
     .     skx(imp1,jmp1,kmp1),sky(imp1,jmp1,kmp1),skz(imp1,jmp1,kmp1)
c
      dimension ivisb(3), amu(mx1dwk)
c
c************* warning: this may not be standard on some compilers  *******
cvn   dimension emax (imp1,jmp1), etemp(imp1)
      dimension emax (imp1,jmp1), etemp(mx1dwk)
c********************************
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c-------   preconditioning related quantities      -----------------
c
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubfil2,iprec
c
      common/rk/    cfl,c(6),qfil(6),beta(6),vt,hm,mstage
c
      common/rkdis/ vis0,vis2,vis4,zeta,vepsn,vepsl,enteps,icau
c
      common/rkrsm/ smoopi,smoopj,smoopk,smoopic,smoopjc,smoopkc
c
      common/timek/  ct(9)
c
      common/vis   / diffac,rsvis,omega
      common/unscom1/ anin,aninb,bninb,betvis,betvisv
      common/unsprec/ npcons,nexp,iprecg,ijac,icord,inav
c
      difac    = diffac
      ai2      = real(i2dfl)
      ai2m     = real(1-i2dfl)
      ivis     = 0
      vismultx = 0.
      vismulty = 0.
      vismultz = 0.
c              ivis=1 if it is viscous in some direction
      if (ivisb(1).gt.0.)  then
        vismultx = 1.
        ivis     = 1
      endif
      if (ivisb(2).gt.0.)  then
        vismulty = 1.
        ivis     = 1
      endif
      if (ivisb(3).gt.0.)  then
        vismultz = 1.
        ivis     = 1
      endif
c
      do 25 k=kmn+1,kmx
c******************************************************************
c     determine maximum value of eddy viscosity for diffusion limit
c
      do 5 i=imn+1,imx
      do j=jmn,jmx+1
        etemp(j) = eomu(i,j,k)
      enddo
      do j=jmn+1,jmx
        emax(i,j)  = max(etemp(j-1),etemp(j),etemp(j+1)) + 1.
      enddo
    5 continue
c*******************************************************************
      do 25 j=jmn+1,jmx
c 
c--------- compute molecular viscosity
c
      do i=imn+1,imx
        ts        =  p(i,j,k)/w(i,j,k,1) 
        amu(i)    =  ts
      enddo
      if(suthc.gt.0.) then
c
        visc1     = suthc/tref
        do i=imn+1,imx
          ts      =  p(i,j,k)/w(i,j,k,1) 
          amu(i)  = ((1.+visc1)/(ts+visc1))*ts**1.5
	enddo
      endif
c
      if (iprec.eq.0) then
c               no preconditioning
c
      do 10 i=imn+1,imx
        rhoi    = 1./w(i,j,k,1)
        c2      = gamma*p(i,j,k)*rhoi
c
c---------   step limit for differences in i direction
c
        sx      = 0.5*(six(i,j,k) + six(i-1,j,k))
        sy      = 0.5*(siy(i,j,k) + siy(i-1,j,k))
        sz      = 0.5*(siz(i,j,k) + siz(i-1,j,k))
        qs      = (w(i,j,k,2)*sx + w(i,j,k,3)*sy + w(i,j,k,4)*sz)*rhoi
        anux    = sx**2 + sy**2 + sz**2
        xix     = abs(sx)
        xiy     = abs(sy)
        xiz     = abs(sz)
        cs      = c2*anux
        dtli    = abs(qs) + sqrt(cs) 
c
c---------   step limit for differences in j direction
c
        sx      = 0.5*(sjx(i,j,k) + sjx(i,j-1,k))
        sy      = 0.5*(sjy(i,j,k) + sjy(i,j-1,k))
        sz      = 0.5*(sjz(i,j,k) + sjz(i,j-1,k))
        qs      = (w(i,j,k,2)*sx + w(i,j,k,3)*sy + w(i,j,k,4)*sz)*rhoi
        anuy    = sx**2 + sy**2 + sz**2
        etax    = abs(sx)
        etay    = abs(sy)
        etaz    = abs(sz)
        cs      = c2*anuy
        dtlj    = abs(qs) + sqrt(cs)  
c
c---------   step limit for differences in k direction
c
        sx      = 0.5*(skx(i,j,k) + skx(i,j,k-1))
        sy      = 0.5*(sky(i,j,k) + sky(i,j,k-1))
        sz      = 0.5*(skz(i,j,k) + skz(i,j,k-1))
        qs      = (w(i,j,k,2)*sx + w(i,j,k,3)*sy + w(i,j,k,4)*sz)*rhoi
        anuz    = sx**2 + sy**2 + sz**2
        zetax   = abs(sx)
        zetay   = abs(sy)
        zetaz   = abs(sz)
        cs      = c2*anuz
        dtlk    = abs(qs) + sqrt(cs) 
c
c--------    diffusion time step limit
c
c       anu     = sqrt(gamma)*rm*gamma/
c       anu     = sqrt(gamma)*rm*gamma*emax(i,j)/
c       anu     = sqrt(gamma)*rm*gamma*(1.+eomu(i,j,k))/
        anu     = sqrt(gamma)*rm*gamma*emax(i,j)/
     .                (w(i,j,k,1)*pr*rey)
        anu     = anu*amu(i)
c
        dtlvi   = vismultx*anu*(anux + pr/(3.*gamma)*( abs(xix*  xiy) +
     .                     abs(xiy*  xiz  ) +   abs(xiz*  xix) ))
        dtlvj   = vismulty*anu*(anuy + pr/(3.*gamma)*( abs(etax*etay) +
     .                     abs(etay*etaz )  +  abs(etaz*etax) ))
        dtlvk   = vismultz*anu*(anuz + pr/(3.*gamma)*( abs(zetax*zetay)+
     .                     abs(zetay*zetaz) +  abs(zetaz*zetax)))
c
c
        dtvi(i,j,k) = dtlvi/vol(i,j,k)
        dtvj(i,j,k) = dtlvj/vol(i,j,k)
        dtvk(i,j,k) = dtlvk/vol(i,j,k)
c
c       dti(i,j,k) = dtli + difac*dtlvi/vol(i,j,k)
c       dtj(i,j,k) = dtlj + difac*dtlvj/vol(i,j,k)
c       dtk(i,j,k) = dtlk + difac*dtlvk/vol(i,j,k)
        dti(i,j,k) = dtli
        dtj(i,j,k) = dtlj
        dtk(i,j,k) = dtlk
c
        fbeta2(i,j,k)  = 1.
        fbetav2(i,j,k) = 1.
   10 continue
c
        else

c
c          time step for preconditioning
c
c        if iprec>0 and iprecg<0 then no preconditioning in iteration
c               but preconditioning needed for artificial viscosity
c               iprec.gt.0    and all iprecg
      do 20 i=imn+1,imx
        rhoi    = 1./w(i,j,k,1)
        c2      = gamma*p(i,j,k)*rhoi
c
        sx      = 0.5*(six(i,j,k) + six(i-1,j,k))
        sy      = 0.5*(siy(i,j,k) + siy(i-1,j,k))
        sz      = 0.5*(siz(i,j,k) + siz(i-1,j,k))
        qsx     = (w(i,j,k,2)*sx + w(i,j,k,3)*sy + w(i,j,k,4)*sz)*rhoi
        anux    = sx**2 + sy**2 + sz**2
        xix     = abs(sx)
        xiy     = abs(sy)
        xiz     = abs(sz)
c
        sx      = 0.5*(sjx(i,j,k) + sjx(i,j-1,k))
        sy      = 0.5*(sjy(i,j,k) + sjy(i,j-1,k))
        sz      = 0.5*(sjz(i,j,k) + sjz(i,j-1,k))
        qsy     = (w(i,j,k,2)*sx + w(i,j,k,3)*sy + w(i,j,k,4)*sz)*rhoi
        anuy    = sx**2 + sy**2 + sz**2
        etax    = abs(sx)
        etay    = abs(sy)
        etaz    = abs(sz)
c
        sx      = 0.5*(skx(i,j,k) + skx(i,j,k-1))
        sy      = 0.5*(sky(i,j,k) + sky(i,j,k-1))
        sz      = 0.5*(skz(i,j,k) + skz(i,j,k-1))
        qsz     = (w(i,j,k,2)*sx + w(i,j,k,3)*sy + w(i,j,k,4)*sz)*rhoi
        anuz    = sx**2 + sy**2 + sz**2
        zetax   = abs(sx)
        zetay   = abs(sy)
        zetaz   = abs(sz)
c
c--------- step limit for differences in i direction
c
        uc      = w(i,j,k,2)*rhoi
        vc      = w(i,j,k,3)*rhoi
        wc      = w(i,j,k,4)*rhoi
        us0     = uc**2 + vc**2 + wc**2
c
c       evaluate alpha and beta for preconditioning
c       store bc in the array fbeta2 for use in precong routine
c
c       vs      = us0
        dvol    = anux + anuy + ai2m*anuz
c       vsold   = (qsx**2 + qsy**2 + ai2m*qsz**2)/dvol
        vsnew   = qsx*qsx/anux + qsy*qsy/anuy + ai2m*qsz*qsz/anuz
c       vs      = vsold
        vs      = vsnew
        eig     = ct(ntorder+1)*vol(i,j,k) /(dvol*sqrt(c2)*dtphy)
        eig     = sqrt(eig)
        am      = vs/c2
c          bc is beta^2/c^2 ~ M^2
c       bc      = max(upk*am,ubmin2)
        bc      = upk*am + aninb*eig + ubmin2
        bcv     = upk*am + bninb*eig + ubfil2
        bc      = min(bc,1.)
        if (iprecg .le. 0) bc=1.
        bcv     = min(bcv,1.)
c       bcv     = bc
        fbeta2(i,j,k)  = bc
        fbetav2(i,j,k) = bcv
c
        beta2   = bc*c2
        alp     = max(ulprec*(1.-bc),0.)
        cs1     = 0.5*(1. - alp + bc)
c
c---------   step limit for differences in i direction
c
        csq     = cs1*qsx
        cs      = sqrt(csq**2 + anux*beta2 - qsx*qsx*bc)
        dtli    = abs(csq) + cs
c
c---------   step limit for differences in j direction
c
        csq     = cs1*qsy
        cs      = sqrt(csq**2 + anuy*beta2 - qsy*qsy*bc)
        dtlj    = abs(csq) + cs
c
c--------- step limit for differences in k direction
c
        csq     = cs1*qsz
        cs      = sqrt(csq**2 + anuz*beta2 - qsz*qsz*bc)
        dtlk    = abs(csq) + cs
c
        dti(i,j,k) = dtli
        dtj(i,j,k) = dtlj
        dtk(i,j,k) = dtlk
c
       if (ivis.gt.0) then
c
c--------   diffusion time step limit
c
c       anu     = sqrt(gamma)*rm*gamma/
cvnv    anu     = sqrt(gamma)*rm*gamma*(1.+eomu(i,j,k))/
        anu     = sqrt(gamma)*rm*gamma*emax(i,j)/
     .                 (w(i,j,k,1)*pr*rey)
        anu     = anu*amu(i)
c
        dtlvi   = vismultx*anu*(anux + pr/(3.*gamma)*( abs(xix*  xiy) +
     .                     abs(xiy*  xiz  ) +   abs(xiz*  xix) ))
        dtlvj   = vismulty*anu*(anuy + pr/(3.*gamma)*( abs(etax*etay) +
     .                     abs(etay*etaz )  +  abs(etaz*etax) ))
        dtlvk   = vismultz*anu*(anuz + pr/(3.*gamma)*( abs(zetax*zetay)+
     .                     abs(zetay*zetaz) +  abs(zetaz*zetax)))
c
c                 viscous contribution to beta
c
        dtlvis     = max(dtlvi,dtlvj,dtlvk)/vol(i,j,k)
        dtlinv     = max(dtli,dtlj,dtlk)
c       dtlinv     = max(abs(qsx),abs(qsy),abs(qsz))
c         reyc is inverse of local cell Reynolds number
        reyc       = dtlvis/dtlinv
c
        bc         = (1. + 3.0*betvis*sqrt(difac*reyc))*fbeta2(i,j,k)
        bcv        = (1. + 3.0*betvisv*sqrt(difac*reyc))*fbetav2(i,j,k)
        bc      = min(bc,1.)
        bcv     = min(bcv,1.)
c
        fbeta2(i,j,k)  = bc
        fbetav2(i,j,k) = bcv
c
c--------   diffusion time step limit
c
        dtvi(i,j,k) = dtlvi/vol(i,j,k)
        dtvj(i,j,k) = dtlvj/vol(i,j,k)
        dtvk(i,j,k) = dtlvk/vol(i,j,k)
c
       else
        dtvi(i,j,k) = 0.
        dtvj(i,j,k) = 0.
        dtvk(i,j,k) = 0.
       endif
   20 continue
      endif
c
   25 continue
c 
c 
c     set values along boundaries by extrapolation
c
      imxp = imx + 1
      jmxp = jmx + 1
      kmxp = kmx + 1
      imnp = imn + 1
      jmnp = jmn + 1
      kmnp = kmn + 1
c
      do k=kmnp,kmx
      do j=jmnp,jmx
        dti(imn ,j,k) = dti(imnp,j,k)
        dtj(imn ,j,k) = dtj(imnp,j,k)
        dtk(imn ,j,k) = dtk(imnp,j,k)
        dti(imxp,j,k) = dti(imx ,j,k)
        dtj(imxp,j,k) = dtj(imx ,j,k)
        dtk(imxp,j,k) = dtk(imx ,j,k)
c
        dtvi(imn ,j,k) = dtvi(imnp,j,k)
        dtvj(imn ,j,k) = dtvj(imnp,j,k)
        dtvk(imn ,j,k) = dtvk(imnp,j,k)
        fbeta2(imn ,j,k)  = fbeta2(imnp,j,k)
        fbetav2(imn ,j,k)  = fbetav2(imnp,j,k)
        dtvi(imxp,j,k) = dtvi(imx ,j,k)
        dtvj(imxp,j,k) = dtvj(imx ,j,k)
        dtvk(imxp,j,k) = dtvk(imx ,j,k)
        fbeta2(imxp,j,k)  = fbeta2(imx ,j,k)
        fbetav2(imxp,j,k)  = fbetav2(imx ,j,k)
      enddo
      enddo
c
      do k=kmnp,kmx
      do i=imn ,imxp
        dti(i,jmn ,k) = dti(i,jmnp,k)
        dtj(i,jmn ,k) = dtj(i,jmnp,k)
        dtk(i,jmn ,k) = dtk(i,jmnp,k)
        dti(i,jmxp,k) = dti(i,jmx ,k)
        dtj(i,jmxp,k) = dtj(i,jmx ,k)
        dtk(i,jmxp,k) = dtk(i,jmx ,k)
c
        dtvi(i,jmn ,k) = dtvi(i,jmnp,k)
        dtvj(i,jmn ,k) = dtvj(i,jmnp,k)
        dtvk(i,jmn ,k) = dtvk(i,jmnp,k)
        fbeta2(i,jmn ,k)  = fbeta2(i,jmnp,k)
        fbetav2(i,jmn ,k)  = fbetav2(i,jmnp,k)
        dtvi(i,jmxp,k) = dtvi(i,jmx ,k)
        dtvj(i,jmxp,k) = dtvj(i,jmx ,k)
        dtvk(i,jmxp,k) = dtvk(i,jmx ,k)
        fbeta2(i,jmxp,k)  = fbeta2(i,jmx ,k)
        fbetav2(i,jmxp,k)  = fbetav2(i,jmx ,k)
      enddo
      enddo
c
      do j=jmn,jmxp
      do i=imn,imxp
        dti(i,j,kmn ) = dti(i,j,kmnp)
        dtj(i,j,kmn ) = dtj(i,j,kmnp)
        dtk(i,j,kmn ) = dtk(i,j,kmnp)
        dti(i,j,kmxp) = dti(i,j,kmx )
        dtj(i,j,kmxp) = dtj(i,j,kmx )
        dtk(i,j,kmxp) = dtk(i,j,kmx )
c
        dtvi(i,j,kmn ) = dtvi(i,j,kmnp)
        dtvj(i,j,kmn ) = dtvj(i,j,kmnp)
        dtvk(i,j,kmn ) = dtvk(i,j,kmnp)
        fbeta2(i,j,kmn )  = fbeta2(i,j,kmnp)
        fbetav2(i,j,kmn )  = fbetav2(i,j,kmnp)
        dtvi(i,j,kmxp) = dtvi(i,j,kmx )
        dtvj(i,j,kmxp) = dtvj(i,j,kmx )
        dtvk(i,j,kmxp) = dtvk(i,j,kmx )
        fbeta2(i,j,kmxp)  = fbeta2(i,j,kmx )
        fbetav2(i,j,kmxp)  = fbetav2(i,j,kmx )
      enddo
      enddo
c
c 
      if(ibloc .eq. 1) dtmin     = 1.e+10
      if (vt.eq.0.) return
      do k=kmn+1,kmx
      do j=jmn+1,jmx
      do i=imn+1,imx
c       dtli   =  1./dti(i,j,k)
c       dtlj   =  1./dtj(i,j,k)
c       dtlk   =  1./dtk(i,j,k)
c
c       dtl    = vol(i,j,k)*dtli*dtlj*dtlk/
c    .             (dtli*dtlj  +dtlj*dtlk  +dtlk*dtli)
c       dtl    = vol(i,j,k)/(dti(i,j,k)+dtj(i,j,k)+dtk(i,j,k))
c
c       dtl    = ai2m*dtl + ai2*vol(i,j,k)*dtli*dtlj/(dtli+dtlj)
c       dtl    = ai2m*dtl + ai2*vol(i,j,k)/(dti(i,j,k)+dtj(i,j,k))
c
        dtli   =  dti(i,j,k)
        dtlj   =  dtj(i,j,k)
        dtlk   =  dtk(i,j,k)
        dtl    =  vol(i,j,k)/(dtli + dtlj + ai2m*dtlk)
c
        dtmin  = min(dtmin,dtl)
      enddo
      enddo
      enddo
c
      return
      end
