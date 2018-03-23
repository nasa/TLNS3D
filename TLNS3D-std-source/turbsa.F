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
      subroutine turbsa(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2,w,p,eomu,turv1,vol,x,
     .                  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                  bx,cx,dx,by,cy,dy,bz,cz,dz,
     .                  fx,fy,fz,wk2dx,wk2dy,wk2dz,
     .                  fv, dtl, dti, dtj, dtk,dtvi,dtvj,dtvk, 
     .                  mx1dwk,turv1s, anu, amu, rhs, rhsimp,
     .                  wk1da,wk1db,mstrt,smin ,ibloc,iter1,
     .                  itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,
     .                  ntorder,tv1old,dtphy ,
     .                  mebdf,ibdf2opt,tv1bar,
     .                  irkstg,nrkstg,mrkstg,
     .                  rkphysa,rkphyse,dtv1ork,itime)
c#ifdef  RKPHY
c# else
c    .                  ntorder,tv1old,dtphy )
c# endif
c----------------------------------------------------------------------
c 
c 
c       object : solve the 1-equation Spalart's turbulence model
c       initial coding : by v.n.vatsa of nasa langley (november -1992)
c 
c 
c 
c

      common/fld/  gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,hh0,c0,u0,v0,w0,
     .             ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubcon2,iprec
      common/trbfs/ anutinf,zkinf,ominf
      common/timek/  ct(9)
      common/vis   / diffac,rsvis,omega
      common/unscom1/ anin,aninb,bninb,betvis,betvisv
      common/unsprec/ npcons,nexp,iprecg,ijac,icord,inav
c
      dimension 
     .    w(imp2,jmp2,kmp2,5),p(imp2,jmp2,kmp2),vol(imp2,jmp2,kmp2),
     .    x(imp1,jmp1,kmp1,3),fv(imp2,jmp2,kmp2,4),
     . turv1(imp2,jmp2,kmp2),turv1s(imp2,jmp2,kmp2),
     .  anu(imp2,jmp2,kmp2),amu(imp2,jmp2,kmp2),eomu(imp2,jmp2,kmp2),
     .  six(imp1,jmp1,kmp1),siy(imp1,jmp1,kmp1),siz(imp1,jmp1,kmp1),
     .  sjx(imp1,jmp1,kmp1),sjy(imp1,jmp1,kmp1),sjz(imp1,jmp1,kmp1),
     .  skx(imp1,jmp1,kmp1),sky(imp1,jmp1,kmp1),skz(imp1,jmp1,kmp1),
     .  dti(imp2,jmp2,kmp2),dtj(imp2,jmp2,kmp2),dtk(imp2,jmp2,kmp2),
     .  dtvi(imp2,jmp2,kmp2),dtvj(imp2,jmp2,kmp2),dtvk(imp2,jmp2,kmp2),
     .  dtl(imp2,jmp2,kmp2),smin(imp2,jmp2,kmp2),
     .  rhs(imp2,jmp2,kmp2),rhsimp(imp2,jmp2,kmp2)
c
      dimension bx(mx1dwk,mx1dwk),by(mx1dwk,mx1dwk),bz(mx1dwk,mx1dwk),
     .          cx(mx1dwk,mx1dwk),cy(mx1dwk,mx1dwk),cz(mx1dwk,mx1dwk),
     .          dx(mx1dwk,mx1dwk),dy(mx1dwk,mx1dwk),dz(mx1dwk,mx1dwk),
     .          fx(mx1dwk,mx1dwk),fy(mx1dwk,mx1dwk),fz(mx1dwk,mx1dwk),
     .          wk2dx(mx1dwk,mx1dwk),wk2dy(mx1dwk,mx1dwk),
     .          wk2dz(mx1dwk,mx1dwk),wk1da(mx1dwk),wk1db(mx1dwk)
c
c     array allocation for time dependent quantities
c
      dimension     tv1old(imp2,jmp2,kmp2,ntorder)
      dimension     tv1bar(imp2,jmp2,kmp2,ntorder-1)
c
c#ifdef  RKPHY
c
      dimension     rkphysa(mrkstg,mrkstg),rkphyse(mrkstg)
c
      dimension     dtv1ork (imp2,jmp2,kmp2,nrkstg)
c#endif
c
      data init/0/, cfltrb/25./
cvn   data init/0/, cfltrb/10./
c
      tol   = 1.e-07
c
c     fv(i,j,k,1) : used for vorticity
c     fv(i,j,k,3) : used for ?
c     fv(i,j,k,4) : used for ?
c
c
c     temporary for debugging purposes only
      if (init.eq.-1) then
        idim = imx - (imn-1)
        jdim = jmx - (jmn-1)
        kdim = kmx - (kmn-1)
        open(71,file="/scr/vatsa/dis.q",form="unformatted")
        rewind 71
        idum = 1
        write (71) idum
        write (71) idim,jdim,kdim
        dum = 0.5
        write (71) dum,dum,dum,dum
        write (71) (((fv(i,j,k,1),i=imn,imx),j=jmn,jmx),k=kmn,kmx),
     .             (((fv(i,j,k,2),i=imn,imx),j=jmn,jmx),k=kmn,kmx),
     .             ((( x(i,j,k,1),i=imn,imx),j=jmn,jmx),k=kmn,kmx),
     .             ((( x(i,j,k,2),i=imn,imx),j=jmn,jmx),k=kmn,kmx),
     .             ((( x(i,j,k,3),i=imn,imx),j=jmn,jmx),k=kmn,kmx)
        stop 5000
      endif
c
c
c Set up constants
      akarman=.41
      cb1=.1355
      sigma=2./3.
      cb2=.622
      cw1=cb1/akarman**2+(1.+cb2)/sigma
      cw2=.30
      cw3=2.0
      cw4=(1.+cb2)*akarman**2/(sigma*cb1)
      cv1=7.1
      ct1=1.0
      ct2=2.0
      ct3=1.2
      ct4=0.5
cold  ct3=1.1
cold  ct4=2.0
cnew  ct3=1.2
cnew  ct4=0.5
c  Set up some other needed parameters
      sfac = sqrt(gamma)*rm/(sigma*rey)
      sfacb= sqrt(gamma)*rm/rey
c
c Set number of subiterations to solve turbulence field eqn per iteration
c (usually, 1 is sufficient... but if residual diverges then may need more)
c NOTE:  Extra subiterations are expensive! (approx 3.e-4 sec/gridpoint
c per extra subiteration on CRAY YMP)
c
c     iter  = iter1
c     nsubit=2
c     if(init.le.5.and.mstrt.eq.0) nsubit = 30
ccc   nsubit = 1
c
c***********************************************************************
c     mods to correct the subiteration error (vatsa, July 95)
      nsubit = iter1
c***********************************************************************
c
c     compute the time step
c
      ai2  = real(i2dfl)
      ai2m = real(1-i2dfl)
      difac    = diffac
c     if (iprecg.ge.1) difac=8
      if (iprec.ge.1) difac=8
      do k=kmn,kmp1
      do j=jmn,jmp1
      do i=imn,imp1
c
c
c     dtli     =  1./dti(i,j,k)
c     dtlj     =  1./dtj(i,j,k)
c     dtlk     =  1./dtk(i,j,k)
c
c     dtl(i,j,k)  = vol(i,j,k)*dtli*dtlj*dtlk/
c    .             (dtli*dtlj  +dtlj*dtlk  +dtlk*dtli)
c
c     dtl(i,j,k)  = (1.-i2dfl)*dtl(i,j,k) + i2dfl*vol(i,j,k)*
c    .                                      dtli*dtlj/(dtli+dtlj)
        dtli       =  dti(i,j,k) + difac*dtvi(i,j,k)
        dtlj       =  dtj(i,j,k) + difac*dtvj(i,j,k)
        dtlk       =  dtk(i,j,k) + difac*dtvk(i,j,k)
        eigt       =  anin*ct(ntorder+1)*vol(i,j,k) / dtphy
        dtl(i,j,k) =  vol(i,j,k)/(eigt + dtli + dtlj + ai2m*dtlk)
      enddo
      enddo
      enddo
c
c     compute molecular viscosities
c
      visc1     = suthc/tref
c
      do 20 k=kmn,kmp1
      do 20 j=jmn,jmp1
      do i=imn,imp1
        ts         =  p(i,j,k)/w(i,j,k,1) 
        amu(i,j,k) = ts
      enddo
c
      if(suthc.gt.0.) then
c
        do i=imn,imp1
          ts        =  p(i,j,k)/w(i,j,k,1) 
          amu(i,j,k)=((1.+visc1)/(ts+visc1))*ts**1.5
        enddo
      endif
c
      do 20 i=imn,imp1
      anu(i,j,k) = amu(i,j,k)/w(i,j,k,1)
   20 continue
c
c     compute eddy-viscosity from turv1 
c  
      do k=kmn+1,kmx
        do j=jmn+1,jmx
          do i=imn+1,imx
c
          chi     = turv1(i,j,k)/anu(i,j,k)
          fv1     = chi**3/(chi**3 + cv1**3)
          eomu(i,j,k)  = fv1*turv1(i,j,k)*w(i,j,k,1)/amu(i,j,k)
          enddo
        enddo
      enddo
c
c
c     start outer iteration loop for the 1-equation model
c
      do 1000 iter = 1,nsubit
c
c     Evaluate  Explicit part of F-xi-xi viscous terms
c
      do 110 k=kmn+1,kmx
      do 110 j=jmn+1,jmx
c
      do 115 i=imn+1,imx
c
        xp    =  six(i,j,k)/(.5*(vol(i,j,k)+vol(i+1,j,k)))
        yp    =  siy(i,j,k)/(.5*(vol(i,j,k)+vol(i+1,j,k)))
        zp    =  siz(i,j,k)/(.5*(vol(i,j,k)+vol(i+1,j,k)))
        xm    =  six(i-1,j,k)/(.5*(vol(i,j,k)+vol(i-1,j,k)))
        ym    =  siy(i-1,j,k)/(.5*(vol(i,j,k)+vol(i-1,j,k)))
        zm    =  siz(i-1,j,k)/(.5*(vol(i,j,k)+vol(i-1,j,k)))
        xa    = 0.5*(six(i,j,k)+six(i-1,j,k))/vol(i,j,k)
        ya    = 0.5*(siy(i,j,k)+siy(i-1,j,k))/vol(i,j,k)
        za    = 0.5*(siz(i,j,k)+siz(i-1,j,k))/vol(i,j,k)
c
c       ttp = xix(i)*xix(i+1/2)+xiy(i)*xiy(i+1/2)+xiz(i)*xiz(i+1/2)
c       ttm = xix(i-1)*xix(i+1/2)+xiy(i-1)*xiy(i+1/2)+xiz(i-1)*xiz(i+1/2)
c
        ttp   = xp*xa + yp*ya + zp*za
        ttm   = xm*xa + ym*ya + zm*za
        cnud  = -cb2*turv1(i,j,k)*sfac
        cap   = ttp*cnud
        cam   = ttm*cnud
        anutp = 0.5*(turv1(i,j,k) + turv1(i+1,j,k))
        anutm = 0.5*(turv1(i,j,k) + turv1(i-1,j,k))
        anup  = 0.5*(anu  (i,j,k) + anu  (i+1,j,k))
        anum  = 0.5*(anu  (i,j,k) + anu  (i-1,j,k))
        cdp   = ( anup+(1.+cb2)*anutp ) * ttp *sfac
        cdm   = ( anum+(1.+cb2)*anutm ) * ttm *sfac
        bzz   =  -max(cdm+cam, 0.)
        czz   =   max(cdp+cap, 0.) + max(cdm+cam, 0.)
        dzz   =  -max(cdp+cap, 0.)
        rhs(i,j,k) = -bzz*turv1(i-1,j,k) -czz*turv1(i,j,k)
     .               -dzz*turv1(i+1,j,k)
c
c     Evaluate  Explicit part of F-xi advective terms
c
        uu    = (xa*w(i,j,k,2) + ya*w(i,j,k,3) + za*w(i,j,k,4))/
     .              w(i,j,k,1)
        sgnu  = sign(1., uu)
        app   = 0.5*(1. + sgnu)
        apm   = 0.5*(1. - sgnu)
        rhs(i,j,k) = rhs(i,j,k)
     .             - uu*(app*(turv1(i  ,j,k) - turv1(i-1,j,k))+
     .                   apm*(turv1(i+1,j,k) - turv1(i  ,j,k)) )
  115   continue
  110 continue
c
c     Evaluate  Explicit part of F-eta-eta viscous terms
c
      do 210 k=kmn+1,kmx
      do 210 j=jmn+1,jmx
c
      do 215 i=imn+1,imx
c
c*****************************************************************
c       changes in diffusion terms suggested by Clay (2-22-93)
c         do 211 jj=jmn+1,jmx
c         wk1da(jj) = turv1(i,jj,k)
c211      continue
c         jturvmx   = ismax(jmx-jmn-1,wk1da(jmn+1),1) + jmn 
c         turvmx    = wk1da(jturvmx)
c         do 212 jj=jmn+1,jturvmx
c         wk1db(jj) = turv1(i,jj,k)
c212      continue
c         do 213 jj=jturvmx+1,jmx
c         wk1db(jj) = turvmx
c213      continue
c         wk1db(jmn) = turv1(i,jmn,k)
c         wk1db(jmx+1) = wk1db(jmx)
c*****************************************************************
c
        xp    =  sjx(i,j,k)/(.5*(vol(i,j,k)+vol(i,j+1,k)))
        yp    =  sjy(i,j,k)/(.5*(vol(i,j,k)+vol(i,j+1,k)))
        zp    =  sjz(i,j,k)/(.5*(vol(i,j,k)+vol(i,j+1,k)))
        xm    =  sjx(i,j-1,k)/(.5*(vol(i,j,k)+vol(i,j-1,k)))
        ym    =  sjy(i,j-1,k)/(.5*(vol(i,j,k)+vol(i,j-1,k)))
        zm    =  sjz(i,j-1,k)/(.5*(vol(i,j,k)+vol(i,j-1,k)))
        xa    = 0.5*(sjx(i,j,k)+sjx(i,j-1,k))/vol(i,j,k)
        ya    = 0.5*(sjy(i,j,k)+sjy(i,j-1,k))/vol(i,j,k)
        za    = 0.5*(sjz(i,j,k)+sjz(i,j-1,k))/vol(i,j,k)
c
c  ttp = etax(j)*etax(j+1/2)+etay(j)*etay(j+1/2)+etaz(j)*etaz(j+1/2)
c  ttm = etax(j-1)*etax(j+1/2)+etay(j-1)*etay(j+1/2)+etaz(j-1)*xiz(j+1/2)
c
        ttp   = xp*xa + yp*ya + zp*za
        ttm   = xm*xa + ym*ya + zm*za
        cnud  = -cb2*turv1(i,j,k)*sfac
clay    cnud  = -cb2*wk1db(j)*sfac
        cap   = ttp*cnud
        cam   = ttm*cnud
        anutp = 0.5*(turv1(i,j,k) + turv1(i,j+1,k))
        anutm = 0.5*(turv1(i,j,k) + turv1(i,j-1,k))
c*****************************************************************
c       changes in diffusion terms suggested by Clay (2-22-93)
c       anutp = 0.5*(wk1db(j) + wk1db(j+1))
c       anutm = 0.5*(wk1db(j) + wk1db(j-1))
c*****************************************************************
        anup  = 0.5*(anu  (i,j,k) + anu  (i,j+1,k))
        anum  = 0.5*(anu  (i,j,k) + anu  (i,j-1,k))
        cdp   = ( anup+(1.+cb2)*anutp ) * ttp *sfac
        cdm   = ( anum+(1.+cb2)*anutm ) * ttm *sfac
        bzz   =  -max(cdm+cam, 0.)
        czz   =   max(cdp+cap, 0.) + max(cdm+cam, 0.)
        dzz   =  -max(cdp+cap, 0.)
        rhs(i,j,k) =  rhs(i,j,k)
     .               -bzz*turv1(i,j-1,k) -czz*turv1(i,j,k)
     .               -dzz*turv1(i,j+1,k)
c
c     Evaluate  Explicit part of F-eta advective terms
c
        uu    = (xa*w(i,j,k,2) + ya*w(i,j,k,3) + za*w(i,j,k,4))/
     .              w(i,j,k,1)
        sgnu  = sign(1., uu)
        app   = 0.5*(1. + sgnu)
        apm   = 0.5*(1. - sgnu)
        rhs(i,j,k) = rhs(i,j,k)
     .             - uu*(app*(turv1(i,j,k  ) - turv1(i,j-1,k))+
     .                   apm*(turv1(i,j+1,k) - turv1(i,j,k  )) )
  215   continue
  210 continue
c
c    Evaluate  Explicit part of F-zeta-zeta viscous terms
c
      if(i2dfl.eq.1) go to 311
c
      do 310 k=kmn+1,kmx
      do 310 j=jmn+1,jmx
c
      do 315 i=imn+1,imx
c
        xp    =  skx(i,j,k)/(.5*(vol(i,j,k)+vol(i,j,k+1)))
        yp    =  sky(i,j,k)/(.5*(vol(i,j,k)+vol(i,j,k+1)))
        zp    =  skz(i,j,k)/(.5*(vol(i,j,k)+vol(i,j,k+1)))
        xm    =  skx(i,j,k-1)/(.5*(vol(i,j,k)+vol(i,j,k-1)))
        ym    =  sky(i,j,k-1)/(.5*(vol(i,j,k)+vol(i,j,k-1)))
        zm    =  skz(i,j,k-1)/(.5*(vol(i,j,k)+vol(i,j,k-1)))
        xa    = 0.5*(skx(i,j,k)+skx(i,j,k-1))/vol(i,j,k)
        ya    = 0.5*(sky(i,j,k)+sky(i,j,k-1))/vol(i,j,k)
        za    = 0.5*(skz(i,j,k)+skz(i,j,k-1))/vol(i,j,k)
c
c  ttp = zetax(k)*zetax(k+1/2)+zetay(k)*zetay(k+1/2)+zetaz(k)*zetaz(k+1/2)
c  ttm = zetax(k-1)*zetax(k+1/2)+zetay(k-1)*zetay(k+1/2)+zetaz(k-1)*xiz(k+1/2)
c
        ttp   = xp*xa + yp*ya + zp*za
        ttm   = xm*xa + ym*ya + zm*za
        cnud  = -cb2*turv1(i,j,k)*sfac
        cap   = ttp*cnud
        cam   = ttm*cnud
        anutp = 0.5*(turv1(i,j,k) + turv1(i,j,k+1))
        anutm = 0.5*(turv1(i,j,k) + turv1(i,j,k-1))
        anup  = 0.5*(anu  (i,j,k) + anu  (i,j,k+1))
        anum  = 0.5*(anu  (i,j,k) + anu  (i,j,k-1))
        cdp   = ( anup+(1.+cb2)*anutp ) * ttp *sfac
        cdm   = ( anum+(1.+cb2)*anutm ) * ttm *sfac
        bzz   =  -max(cdm+cam, 0.)
        czz   =   max(cdp+cap, 0.) + max(cdm+cam, 0.)
        dzz   =  -max(cdp+cap, 0.)
        rhs(i,j,k) =  rhs(i,j,k)
     .               -bzz*turv1(i,j,k-1) -czz*turv1(i,j,k)
     .               -dzz*turv1(i,j,k+1)
c
c     Evaluate  Explicit part of F-zeta advective terms
c
        uu    = (xa*w(i,j,k,2) + ya*w(i,j,k,3) + za*w(i,j,k,4))/
     .              w(i,j,k,1)
        sgnu  = sign(1., uu)
        app   = 0.5*(1. + sgnu)
        apm   = 0.5*(1. - sgnu)
        rhs(i,j,k) = rhs(i,j,k)
     .             - uu*(app*(turv1(i,j,k  ) - turv1(i,j,k-1))+
     .                   apm*(turv1(i,j,k+1) - turv1(i,j,k  )) )
  315   continue
  310 continue
c
  311 continue
c
c     Evaluate source terms
c
      do 350 k=kmn+1,kmx
      do 350 i=imn+1,imx
c
      cutoff  = 1.
c
c
c     transition specified via itrb1...ktrb2 (9-29-93)
c
      if(i.lt.itrb1.or.i.gt.itrb2) cutoff = 0.
      if(k.lt.ktrb1.or.k.gt.ktrb2) cutoff = 0.
c
      do 352 j=jmn+1,jmx
      wk1da(j) = 1.
      if(j.lt.jtrb1.or.j.gt.jtrb2) wk1da(j) = 0.
 352  continue
c
      do 355 j=jmn+1,jmx
c
        ss    = fv(i,j,k,1)
        ssmin = smin(i,j,k)
c
        chi   = turv1(i,j,k)/anu(i,j,k)
        fv1   = chi**3/(chi**3 + cv1**3)
        fv2   = 1. - chi/(1.+chi*fv1)
        sst   = ss + turv1(i,j,k)*fv2*sfacb/((akarman*ssmin)**2)
cnew
        sst   = max (sst,tol)
        rr    = turv1(i,j,k)*sfacb/(sst*(akarman*ssmin)**2)
        rr    = min (rr,10.)
        gg    = rr + cw2*(rr**6 - rr)
cvn     gg    = max(gg,tol)
c
cvn     fw    = ( (gg**(-6)+cw3**(-6))/(1.+cw3**(-6)) )**(-1./6.)
c from fun3d (may 02)
        fw    = gg*((1. + cw3**6)/(gg**6 + cw3**6))**(1./6.)
        ft2   = ct3*exp(-ct4*chi**2)
        term1 = cb1 * (1.-ft2) * ss
        term2 = cb1 * ((1.-ft2)*fv2 + ft2)/akarman**2 - cw1*fw
        dist2i= sfacb/(ssmin*ssmin + 1.e-20)
        tt    = term1*turv1(i,j,k) + term2*turv1(i,j,k)**2*dist2i
c    
c       update the rhs
c
        rhs   (i,j,k) = rhs(i,j,k) + tt*cutoff*wk1da(j)
c
c       save part of rhs to be added to implicit LHS
c
cvb     rhsimp(i,j,k) =2.*term2*turv1(i,j,k)*dist2i*cutoff*wk1da(j)
        rhsimp(i,j,k) =2.*term2*turv1(i,j,k)*dist2i
c
c       compute source terms arising from chain rule on dissipation
c      (a-la Rumsey, Dec. 97)
c      coded by Vatsa (Dec. 2000)
c
       dfv1  = (fv1 - fv1*fv1)*3./turv1(i,j,k)
       dfv2  = (fv2 - 1.)/turv1(i,j,k) + ((1. - fv2)**2)*
     .         (fv1/turv1(i,j,k) + dfv1)
       dft2  =-(2.*ct4*turv1(i,j,k)/(anu(i,j,k)**2))*ft2
       drr   =  rr/turv1(i,j,k) -rr*rr*(fv2/turv1(i,j,k) + dfv2)
       dgg   = (1.- cw2 + 6.*cw2*(rr**5))*drr
       gg    =  max(gg,tol*10.)
c      dfw   = (gg**(-7)/(1. + (cw3**(-6))))*(fw**7)*dgg
cvn
c mods corresponding to fun3d def. of fw
       dfw   = ((1. + cw3**6)/(gg**6 + cw3**6))**(1./6.) -
     .         ( (1. + cw3**6)**(1./6.)/( (gg**6 + cw3**6)**(7./6.)) )
     .         *gg**6
       dfw   = dfw*dgg
c      modified expression for implicit LHS
       rhsimp(i,j,k) = rhsimp(i,j,k) +
     .                 dist2i*(turv1(i,j,k)**2)*(cb1/(akarman**2)*
     .                (dfv2 - ft2*dfv2 - fv2*dft2 + dft2 ) - cw1*dfw)
       rhsimp(i,j,k) = rhsimp(i,j,k)*cutoff*wk1da(j)
  355   continue
  350 continue
c
c     Evaluate  Implicit part of F-eta-eta viscous terms
c
      do 510 k=kmn+1,kmx
      do 525 i=imn+1,imx
c
c*****************************************************************
c       changes in diffusion terms suggested by Clay (2-22-93)
c         do 511 jj=jmn+1,jmx
c         wk1da(jj) = turv1(i,jj,k)
c511      continue
c         jturvmx   = ismax(jmx-jmn-1,wk1da(jmn+1),1) + jmn 
c         turvmx    = wk1da(jturvmx)
c         do 512 jj=jmn+1,jturvmx
c         wk1db(jj) = turv1(i,jj,k)
c512      continue
c         do 513 jj=jturvmx+1,jmx
c         wk1db(jj) = turvmx
c513      continue
c         wk1db(jmn) = turv1(i,jmn,k)
c         wk1db(jmx+1) = wk1db(jmx)
c*****************************************************************
      do 515 j=jmn+1,jmx
c
        xp    =  sjx(i,j,k)/(.5*(vol(i,j,k)+vol(i,j+1,k)))
        yp    =  sjy(i,j,k)/(.5*(vol(i,j,k)+vol(i,j+1,k)))
        zp    =  sjz(i,j,k)/(.5*(vol(i,j,k)+vol(i,j+1,k)))
        xm    =  sjx(i,j-1,k)/(.5*(vol(i,j,k)+vol(i,j-1,k)))
        ym    =  sjy(i,j-1,k)/(.5*(vol(i,j,k)+vol(i,j-1,k)))
        zm    =  sjz(i,j-1,k)/(.5*(vol(i,j,k)+vol(i,j-1,k)))
        xa    = 0.5*(sjx(i,j,k)+sjx(i,j-1,k))/vol(i,j,k)
        ya    = 0.5*(sjy(i,j,k)+sjy(i,j-1,k))/vol(i,j,k)
        za    = 0.5*(sjz(i,j,k)+sjz(i,j-1,k))/vol(i,j,k)
c
c  ttp = etax(j)*etax(j+1/2)+etay(j)*etay(j+1/2)+etaz(j)*etaz(j+1/2)
c  ttm = etax(j-1)*etax(j+1/2)+etay(j-1)*etay(j+1/2)+etaz(j-1)*xiz(j+1/2)
c
        ttp   = xp*xa + yp*ya + zp*za
        ttm   = xm*xa + ym*ya + zm*za
        cnud  = -cb2*turv1(i,j,k)*sfac
clay    cnud  = -cb2*wk1db(j)*sfac
        cap   = ttp*cnud
        cam   = ttm*cnud
        anutp = 0.5*(turv1(i,j,k) + turv1(i,j+1,k))
        anutm = 0.5*(turv1(i,j,k) + turv1(i,j-1,k))
c*****************************************************************
c       changes in diffusion terms suggested by Clay (2-22-93)
c       anutp = 0.5*(wk1db(j) + wk1db(j+1))
c       anutm = 0.5*(wk1db(j) + wk1db(j-1))
c*****************************************************************
        anup  = 0.5*(anu  (i,j,k) + anu  (i,j+1,k))
        anum  = 0.5*(anu  (i,j,k) + anu  (i,j-1,k))
        cdp   = ( anup+(1.+cb2)*anutp ) * ttp *sfac
        cdm   = ( anum+(1.+cb2)*anutm ) * ttm *sfac
        by(j,i)   =  -max(cdm+cam, 0.)
        cy(j,i)   =   max(cdp+cap, 0.) + max(cdm+cam, 0.)
        dy(j,i)   =  -max(cdp+cap, 0.)
c
c     Evaluate  Implicit part of F-eta advective terms
c
        uu    = (xa*w(i,j,k,2) + ya*w(i,j,k,3) + za*w(i,j,k,4))/
     .              w(i,j,k,1)
        sgnu  = sign(1., uu)
        app   = 0.5*(1. + sgnu)
        apm   = 0.5*(1. - sgnu)
        by(j,i)   = by(j,i) - uu*app
        cy(j,i)   = cy(j,i) + uu*(app-apm)
        dy(j,i)   = dy(j,i) + uu*apm
c
c       add part of the source term to the diagonal in LHS
c
        cy(j,i)   = cy(j,i)
     .            - .5*( 1. - sign(1.,rhsimp(i,j,k)) )*rhsimp(i,j,k)
  515   continue
c
      if(mebdf.eq.1) go to 511
      if(ntorder.eq.0) then
        do 520 j=jmn+1,jmx
cvnv    fact  = 10.*vol(i,j,k)/dtl(i,j,k)
        fact  = cfltrb*dtl(i,j,k)
        fact  = min(fact,100.)
        by(j,i)   = by(j,i)*fact
        cy(j,i)   = cy(j,i)*fact + 1.
        dy(j,i)   = dy(j,i)*fact
        fy(j,i)   = rhs(i,j,k)*fact
  520   continue
        else if(ntorder .eq. 1) then
#ifdef  RKPHY
        if(nrkstg.gt.2) then
c       save the residual without time-terms for later use
          if(irkstg.lt.nrkstg) then
            do j=jmn+1,jmx
             dtv1ork(i,j,k,irkstg) = rhs(i,j,k)
            enddo
          endif
c
          do 519 j=jmn+1,jmx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          rhsold= 0.
          if(irkstg.gt.1) then
            do lstg = 1,irkstg-1
            rhsold = rhsold + rkphysa(irkstg,lstg)*dtv1ork(i,j,k,lstg)
            enddo
          endif
          by(j,i) = by(j,i)*fact*rkphysa(irkstg,irkstg)
          cy(j,i) = cy(j,i)*fact*rkphysa(irkstg,irkstg) + 1.
          dy(j,i) = dy(j,i)*fact*rkphysa(irkstg,irkstg)
          fy(j,i) = rhs(i,j,k)*fact*rkphysa(irkstg,irkstg)
     .            -  fact/dtphy*(turv1(i,j,k)-tv1old(i,j,k,1))
     .            +  fact*rhsold
  519     continue
        endif
# else
          do 521 j=jmn+1,jmx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          by(j,i) = by(j,i)*fact
          cy(j,i) = cy(j,i)*fact + 1.
          dy(j,i) = dy(j,i)*fact
          fy(j,i) = rhs(i,j,k)*fact
     .              -fact/dtphy*(turv1(i,j,k)-tv1old(i,j,k,1))
  521     continue
# endif
        else if(ntorder .eq. 2) then
cvn     temp fix for variable dt
         if(itime.ge.2) then
           h0 = dtphy
           h1 = dtphy
         else
cvn     temp fix for variable dt
           h0 = dtphy
           h1 = dtphy
c          h1 = 0.015625d0
         endif
c
      A4  =  (H1+2*H0)/(H0*(H1+H0))
      A3  = -(H1+H0)/(H0*H1)
      A2  =   H0/(H1*(H1+H0))
      a4dt=   a4*dtphy
c
          do 522 j=jmn+1,jmx
cvn       fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k)*1.5)
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+a4dt*cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          by(j,i) = by(j,i)*fact
          cy(j,i) = cy(j,i)*fact + 1.
          dy(j,i) = dy(j,i)*fact
c         fy(j,i) = rhs(i,j,k)*fact
c    .             -fact/dtphy*(1.5*( turv1(i,j,k)  -tv1old(i,j,k,1))
c    .                         -0.5*(tv1old(i,j,k,1)-tv1old(i,j,k,2)))

          tdir    = a4*turv1(i,j,k) + a3*tv1old(i,j,k,1) +
     .              a2*tv1old(i,j,k,2)

          fy(j,i) = rhs(i,j,k)*fact -fact*tdir
  522     continue
        else if(ntorder .eq. 3) then
         if(itime.ge.3) then
           dt0 = dtphy
           dt1 = dtphy
           dt2 = dtphy
           a4 = (dt1*dt2+2*dt0*dt2+dt1**2+4*dt0*dt1+3*dt0**2)/
     .          (dt0*(dt1+dt0)*(dt2+dt1+dt0))
           a3 = -(dt1*dt2+dt0*dt2+dt1**2+2*dt0*dt1+dt0**2)/
     .           (dt0*dt1*(dt2+dt1))
           a2 = (dt0*dt2+dt0*dt1+dt0**2)/(dt1*(dt1+dt0)*dt2)
           a1 = -(dt0*dt1+dt0**2)/(dt2*(dt2+dt1)*(dt2+dt1+dt0))
c  
           if(ibdf2opt.eq.1) then
             bdfoptc = -0.58d0/3.0d0
             a4  = ( 1.5d0 - bdfoptc)/dtphy
             a3  = (-2.0d0 + bdfoptc*3.0d0)/dtphy
             a2  = ( 0.5d0 - bdfoptc*3.0d0)/dtphy
             a1  = bdfoptc/dtphy
           endif
         else
cvn        temp fix for variable dt
           dt0 = dtphy
           dt1 = dtphy
           dt2 = dtphy
c          dt1 = 0.015625d0
c          dt2 = 0.015625d0
           if(itime.eq.2) dt1 = dtphy
c
           a4 = (dt1*dt2+2*dt0*dt2+dt1**2+4*dt0*dt1+3*dt0**2)/
     .          (dt0*(dt1+dt0)*(dt2+dt1+dt0))
           a3 = -(dt1*dt2+dt0*dt2+dt1**2+2*dt0*dt1+dt0**2)/
     .           (dt0*dt1*(dt2+dt1))
           a2 = (dt0*dt2+dt0*dt1+dt0**2)/(dt1*(dt1+dt0)*dt2)
           a1 = -(dt0*dt1+dt0**2)/(dt2*(dt2+dt1)*(dt2+dt1+dt0))
c
         endif
c
c     scale a4 with dtphy
c
          a4dt = a4*dtphy
          do 523 j=jmn+1,jmx
c         fact  = 6.*cfltrb*dtl(i,j,k)*dtphy/( 6.*dtphy +
c    .                                        11.*cfltrb*dtl(i,j,k) )
          fact=cfltrb*dtl(i,j,k)*dtphy/(dtphy+a4dt*cfltrb*dtl(i,j,k))
c
          fact  = min(fact,100.)
          by(j,i) = by(j,i)*fact
          cy(j,i) = cy(j,i)*fact + 1.
          dy(j,i) = dy(j,i)*fact
          fy(j,i) = rhs(i,j,k)*fact
     .             -fact*(a4*turv1 (i,j,k)    + a3*tv1old(i,j,k,1) +
     .                    a2*tv1old(i,j,k,2)  + a1*tv1old(i,j,k,3) )
c    .             -fact/(6.*dtphy)*
c    .                   (11.*turv1 (i,j,k)   -18.*tv1old(i,j,k,1) +
c    .                     9.*tv1old(i,j,k,2) - 2.*tv1old(i,j,k,3) )
  523     continue
c
        end if
c
  511   continue
c
      if(mebdf.ne.1) go to 525
c
c     MEBDF4 scheme
c
      if(itime.ge.3) then
        H0 = dtphy
        H1 = dtphy
        H2 = dtphy
      else
c     temp fix for variable time
        H0 = dtphy
        H1 = dtphy
        H2 = dtphy
c       H1 = 0.015625d0
c       H2 = 0.015625d0
        if(itime.eq.2) H1 = dtphy
      endif
        
c
      if(irkstg.eq.1) then
c
c     MEBDF4 : stage 1
c
c     Save the fluxes for later stages before adding time terms
c
          do j = jmn+1,jmx
            dtv1ork(i,j,k,irkstg) = rhs(i,j,k)
          enddo
c
      A4 = (H1*H2+2*H0*H2+H1**2+4*H0*H1+3*H0**2)/(H0*(H1+H0)*(H2+H1+H0))
      A3 = -(H1*H2+H0*H2+H1**2+2*H0*H1+H0**2)/(H0*H1*(H2+H1))
      A2 = (H0*H2+H0*H1+H0**2)/(H1*(H1+H0)*H2)
      A1 = -(H0*H1+H0**2)/(H2*(H2+H1)*(H2+H1+H0))

          do j = jmn+1,jmx
            tdir = A4    *  turv1(i,j,k)
     .           + A3    * tv1old(i,j,k,1)
     .           + A2    * tv1old(i,j,k,2)
     .           + A1    * tv1old(i,j,k,3) 
            rhs(i,j,k) = rhs(i,j,k) - tdir
          enddo
c
      else if(irkstg.eq.2) then
c     MEBDF4 : stage 2
c
c     Save the fluxes for later stages before adding time terms
c
          do j = jmn+1,jmx
            dtv1ork(i,j,k,irkstg) = rhs(i,j,k)
          enddo
c
      HF = dtphy
      A4 = (H0*H1+2*HF*H1+H0**2+4*HF*H0+3*HF**2)/(HF*(H0+HF)*(H1+H0+HF))
      A3 = -(H0*H1+HF*H1+H0**2+2*HF*H0+HF**2)/(HF*H0*(H1+H0))
      A2 = (HF*H1+HF*H0+HF**2)/(H0*(H0+HF)*H1)
      A1 = -(HF*H0+HF**2)/(H1*(H1+H0)*(H1+H0+HF))

          do j = jmn+1,jmx
            tdir = A4    *  turv1(i,j,k)
     .           + A3    * tv1bar(i,j,k,1)
     .           + A2    * tv1old(i,j,k,1)
     .           + A1    * tv1old(i,j,k,2) 
            rhs(i,j,k) = rhs(i,j,k) - tdir
          enddo
c
      else if(irkstg.eq.3) then
c     MEBDF4 : stage 3
c
        RbetaK = (H1*H2+2*H0*H2+H1**2+4*H0*H1+3*H0**2)/(H0*(H1+H0)*
     &           (H2+H1+H0))
      HF  = dtphy
      t1  = 0.
      t2  = t1 + h2
      t3  = t2 + h1
      t4  = t3 + h0
      t5  = t4 + hf
      DENOMIN = (t5-t4)*(12*t4**2*t5**2-8*t3*t4*t5**2-8*t2*t4*t5**2+4*t2
     1   *t3*t5**2-8*t3*t4**2*t5-8*t2*t4**2*t5+6*t3**2*t4*t5+10*t2*t3*t4
     2   *t5+6*t2**2*t4*t5-3*t2*t3**2*t5-3*t2**2*t3*t5+4*t2*t3*t4**2-3*t
     3   2*t3**2*t4-3*t2**2*t3*t4+2*t2**2*t3**2)
      betaKH=(t4*(t4-t2)*(t4-t3)*(4*t5**3-3*t4*t5**2-3*t3*t5**2-3*t2*t5*
     1   *2+2*t3*t4*t5+2*t2*t4*t5+2*t2*t3*t5-t2*t3*t4) )/ denomin
      betaKP1H=(-t4**2*(t4-t2)**2*(t4-t3)**2 )/ denomin
      AH4 = 1
      AH3 =(-t4**2*(t4-t2)**2*(t5-t4)*(4*t5**2-2*t4*t5-3*t2*t5+t2*t4)/(t
     1   3*(t3-t2)) )/ denomin
      AH2 =(t4**2*(t4-t3)**2*(t5-t4)*(4*t5**2-2*t4*t5-3*t3*t5+t3*t4)/(t2
     1   *(t3-t2)) )/ denomin
      AH1 =(-(t4-t2)**2*(t4-t3)**2*(t5-t4)*(4*t5**2-2*t4*t5-3*t3*t5-3*t2
     1   *t5+t3*t4+t2*t4+2*t2*t3)/(t2*t3) )/ denomin

      A4 = aH4*Rbetak
      A3 = aH3*Rbetak
      A2 = aH2*Rbetak
      A1 = aH1*Rbetak
      B1  = betaKH*RbetaK - 1
      B2  = betaKp1H*RbetaK
c
          do j = jmn+1,jmx
            tdir = A4    *  turv1(i,j,k)
     .           + A3    * tv1old(i,j,k,1)
     .           + A2    * tv1old(i,j,k,2)
     .           + A1    * tv1old(i,j,k,3) 
     .           - B1 * dtv1ork(i,j,k,1) - B2 * dtv1ork(i,j,k,2)
cvn  .           + B1 * dtv1ork(i,j,k,1) + B2 * dtv1ork(i,j,k,2)
            rhs(i,j,k) = rhs(i,j,k) - tdir
          enddo
c
      endif
c
c       scale a4 with dtphy
c
        a4dt = a4*dtphy
        do j = jmn+1,jmx
          fact=cfltrb*dtl(i,j,k)*dtphy/(dtphy+a4dt*cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          by(j,i) = by(j,i)*fact
          cy(j,i) = cy(j,i)*fact + 1.
          dy(j,i) = dy(j,i)*fact
          fy(j,i) = rhs(i,j,k)*fact
        enddo
c
c     end of logic for MEBDF4
c
c
  525 continue
c
c     solve the scalar tridiagonal equations in j-direction
c
      call tridag(mx1dwk,mx1dwk,jmn+1,jmx,imn+1,imx,wk2dy,by,cy,dy,fy)
c
c     update the r.h.s. for A-F scheme
c
      do 530 i=imn+1,imx
      do 530 j=jmn+1,jmx
      rhs(i,j,k)   = fy(j,i)
  530 continue
c
  510 continue
c
c
c     Evaluate  Implicit part of F-xi-xi viscous terms
c
      do 410 k=kmn+1,kmx
      do 425 j=jmn+1,jmx
c
      do 415 i=imn+1,imx
c
        xp    =  six(i,j,k)/(.5*(vol(i,j,k)+vol(i+1,j,k)))
        yp    =  siy(i,j,k)/(.5*(vol(i,j,k)+vol(i+1,j,k)))
        zp    =  siz(i,j,k)/(.5*(vol(i,j,k)+vol(i+1,j,k)))
        xm    =  six(i-1,j,k)/(.5*(vol(i,j,k)+vol(i-1,j,k)))
        ym    =  siy(i-1,j,k)/(.5*(vol(i,j,k)+vol(i-1,j,k)))
        zm    =  siz(i-1,j,k)/(.5*(vol(i,j,k)+vol(i-1,j,k)))
        xa    = 0.5*(six(i,j,k)+six(i-1,j,k))/vol(i,j,k)
        ya    = 0.5*(siy(i,j,k)+siy(i-1,j,k))/vol(i,j,k)
        za    = 0.5*(siz(i,j,k)+siz(i-1,j,k))/vol(i,j,k)
c
c       ttp = xix(i)*xix(i+1/2)+xiy(i)*xiy(i+1/2)+xiz(i)*xiz(i+1/2)
c       ttm = xix(i-1)*xix(i+1/2)+xiy(i-1)*xiy(i+1/2)+xiz(i-1)*xiz(i+1/2)
c
        ttp   = xp*xa + yp*ya + zp*za
        ttm   = xm*xa + ym*ya + zm*za
        cnud  = -cb2*turv1(i,j,k)*sfac
        cap   = ttp*cnud
        cam   = ttm*cnud
        anutp = 0.5*(turv1(i,j,k) + turv1(i+1,j,k))
        anutm = 0.5*(turv1(i,j,k) + turv1(i-1,j,k))
        anup  = 0.5*(anu  (i,j,k) + anu  (i+1,j,k))
        anum  = 0.5*(anu  (i,j,k) + anu  (i-1,j,k))
        cdp   = ( anup+(1.+cb2)*anutp ) * ttp *sfac
        cdm   = ( anum+(1.+cb2)*anutm ) * ttm *sfac
        bx(i,j)   =  -max(cdm+cam, 0.)
        cx(i,j)   =   max(cdp+cap, 0.) + max(cdm+cam, 0.)
        dx(i,j)   =  -max(cdp+cap, 0.)
c
c     Evaluate  Implicit part of F-xi advective terms
c
        uu    = (xa*w(i,j,k,2) + ya*w(i,j,k,3) + za*w(i,j,k,4))/
     .              w(i,j,k,1)
        sgnu  = sign(1., uu)
        app   = 0.5*(1. + sgnu)
        apm   = 0.5*(1. - sgnu)
        bx(i,j)   = bx(i,j) - uu*app
        cx(i,j)   = cx(i,j) + uu*(app-apm)
        dx(i,j)   = dx(i,j) + uu*apm
  415   continue
c
      if(ntorder .eq. 0) then
        do 420 i=imn+1,imx
cvnv    fact  = 10.*vol(i,j,k)/dtl(i,j,k)
        fact  = cfltrb*dtl(i,j,k)
        fact  = min(fact,100.)
        bx(i,j)   = bx(i,j)*fact
        cx(i,j)   = cx(i,j)*fact + 1.
        dx(i,j)   = dx(i,j)*fact
        fx(i,j)   = rhs(i,j,k)
  420   continue
        else if(ntorder .eq. 1) then
#ifdef  RKPHY
        if(nrkstg.gt.2) then
          do 419 i=imn+1,imx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          bx(i,j) = bx(i,j)*fact*rkphysa(irkstg,irkstg)
          cx(i,j) = cx(i,j)*fact*rkphysa(irkstg,irkstg) + 1.
          dx(i,j) = dx(i,j)*fact*rkphysa(irkstg,irkstg)
          fx(i,j) = rhs(i,j,k)
  419     continue
        endif
# else
          do 421 i=imn+1,imx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          bx(i,j) = bx(i,j)*fact
          cx(i,j) = cx(i,j)*fact + 1.
          dx(i,j) = dx(i,j)*fact
          fx(i,j) = rhs(i,j,k)
  421     continue
# endif
        else if(ntorder .eq. 2) then
         if(itime.ge.2) then
           h0 = dtphy
           h1 = dtphy
         else
cvn     temp fix for variable dt
           h0 = dtphy
           h1 = dtphy
c          h1 = 0.015625d0
         endif
         A4  =  (H1+2*H0)/(H0*(H1+H0))
         a4dt=   a4*dtphy
          do 422 i=imn+1,imx
cvn       fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k)*1.5)
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+a4dt*cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          bx(i,j) = bx(i,j)*fact
          cx(i,j) = cx(i,j)*fact + 1.
          dx(i,j) = dx(i,j)*fact
          fx(i,j) = rhs(i,j,k)
  422     continue
        else if(ntorder .eq. 3) then
         if(itime.ge.3) then
           dt0 = dtphy
           dt1 = dtphy
           dt2 = dtphy
           a4 = (dt1*dt2+2*dt0*dt2+dt1**2+4*dt0*dt1+3*dt0**2)/
     .          (dt0*(dt1+dt0)*(dt2+dt1+dt0))
c  
           if(ibdf2opt.eq.1) then
             bdfoptc = -0.58d0/3.0d0
             a4  = ( 1.5d0 - bdfoptc)/dtphy
           endif
         else
cvn        temp fix for variable dt
           dt0 = dtphy
           dt1 = dtphy
           dt2 = dtphy
c          dt1 = 0.015625d0
c          dt2 = 0.015625d0
           if(itime.eq.2) dt1 = dtphy
c
           a4 = (dt1*dt2+2*dt0*dt2+dt1**2+4*dt0*dt1+3*dt0**2)/
     .          (dt0*(dt1+dt0)*(dt2+dt1+dt0))
c
         endif
c
c     scale a's with dtphy
c
          a4dt = a4*dtphy
          do 423 i=imn+1,imx
c         fact  = 6.*cfltrb*dtl(i,j,k)*dtphy/( 6.*dtphy +
c    .                                        11.*cfltrb*dtl(i,j,k) )
          fact=cfltrb*dtl(i,j,k)*dtphy/(dtphy+a4dt*cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          bx(i,j) = bx(i,j)*fact
          cx(i,j) = cx(i,j)*fact + 1.
          dx(i,j) = dx(i,j)*fact
          fx(i,j) = rhs(i,j,k)
  423     continue
c
        end if
c
  425 continue
c
c     solve the scalar tridiagonal equations in i-direction
c
      call tridag(mx1dwk,mx1dwk,imn+1,imx,jmn+1,jmx,wk2dx,bx,cx,dx,fx)
c
c     update the r.h.s. for A-F scheme
c
      do 430 j=jmn+1,jmx
      do 430 i=imn+1,imx
      rhs(i,j,k)   = fx(i,j)
  430 continue
c
  410 continue
c
c     Evaluate  Implicit part of F-zeta-zeta viscous terms
c
      if(i2dfl.eq.1) go to 611
c
      do 610 i=imn+1,imx
      do 625 j=jmn+1,jmx
c
      do 615 k=kmn+1,kmx
c
        xp    =  skx(i,j,k)/(.5*(vol(i,j,k)+vol(i,j,k+1)))
        yp    =  sky(i,j,k)/(.5*(vol(i,j,k)+vol(i,j,k+1)))
        zp    =  skz(i,j,k)/(.5*(vol(i,j,k)+vol(i,j,k+1)))
        xm    =  skx(i,j,k-1)/(.5*(vol(i,j,k)+vol(i,j,k-1)))
        ym    =  sky(i,j,k-1)/(.5*(vol(i,j,k)+vol(i,j,k-1)))
        zm    =  skz(i,j,k-1)/(.5*(vol(i,j,k)+vol(i,j,k-1)))
        xa    = 0.5*(skx(i,j,k)+skx(i,j,k-1))/vol(i,j,k)
        ya    = 0.5*(sky(i,j,k)+sky(i,j,k-1))/vol(i,j,k)
        za    = 0.5*(skz(i,j,k)+skz(i,j,k-1))/vol(i,j,k)
c
c  ttp = zetax(k)*zetax(k+1/2)+zetay(k)*zetay(k+1/2)+zetaz(k)*zetaz(k+1/2)
c  ttm = zetax(k-1)*zetax(k+1/2)+zetay(k-1)*zetay(k+1/2)+zetaz(k-1)*xiz(k+1/2)
c
        ttp   = xp*xa + yp*ya + zp*za
        ttm   = xm*xa + ym*ya + zm*za
        cnud  = -cb2*turv1(i,j,k)*sfac
        cap   = ttp*cnud
        cam   = ttm*cnud
        anutp = 0.5*(turv1(i,j,k) + turv1(i,j,k+1))
        anutm = 0.5*(turv1(i,j,k) + turv1(i,j,k-1))
        anup  = 0.5*(anu  (i,j,k) + anu  (i,j,k+1))
        anum  = 0.5*(anu  (i,j,k) + anu  (i,j,k-1))
        cdp   = ( anup+(1.+cb2)*anutp ) * ttp *sfac
        cdm   = ( anum+(1.+cb2)*anutm ) * ttm *sfac
        bz(k,j)   =  -max(cdm+cam, 0.)
        cz(k,j)   =   max(cdp+cap, 0.) + max(cdm+cam, 0.)
        dz(k,j)   =  -max(cdp+cap, 0.)
c
c     Evaluate  Implicit part of F-zeta advective terms
c
        uu    = (xa*w(i,j,k,2) + ya*w(i,j,k,3) + za*w(i,j,k,4))/
     .              w(i,j,k,1)
        sgnu  = sign(1., uu)
        app   = 0.5*(1. + sgnu)
        apm   = 0.5*(1. - sgnu)
        bz(k,j)   = bz(k,j) - uu*app
        cz(k,j)   = cz(k,j) + uu*(app-apm)
        dz(k,j)   = dz(k,j) + uu*apm
c
  615   continue
c
      if(ntorder .eq. 0) then
        do 620 k=kmn+1,kmx
cvnv    fact  = 10.*vol(i,j,k)/dtl(i,j,k)
        fact  = cfltrb*dtl(i,j,k)
        fact  = min(fact,100.)
        bz(k,j)   = bz(k,j)*fact
        cz(k,j)   = cz(k,j)*fact + 1.
        dz(k,j)   = dz(k,j)*fact
        fz(k,j)   = rhs(i,j,k)
  620   continue
        else if(ntorder .eq. 1) then
#ifdef  RKPHY
        if(nrkstg.gt.2) then
          do 619 k=kmn+1,kmx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          bz(k,j) = bz(k,j)*fact*rkphysa(irkstg,irkstg)
          cz(k,j) = cz(k,j)*fact*rkphysa(irkstg,irkstg) + 1.
          dz(k,j) = dz(k,j)*fact*rkphysa(irkstg,irkstg)
          fz(k,j) = rhs(i,j,k)
  619     continue
        endif
# else
          do 621 k=kmn+1,kmx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          bz(k,j) = bz(k,j)*fact
          cz(k,j) = cz(k,j)*fact + 1.
          dz(k,j) = dz(k,j)*fact
          fz(k,j) = rhs(i,j,k)
  621     continue
# endif
        else if(ntorder .eq. 2) then
         if(itime.ge.2) then
           h0 = dtphy
           h1 = dtphy
         else
cvn     temp fix for variable dt
           h0 = dtphy
           h1 = dtphy
c          h1 = 0.015625d0
         endif
         A4  =  (H1+2*H0)/(H0*(H1+H0))
         a4dt=   a4*dtphy
          do 622 k=kmn+1,kmx
cvn       fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k)*1.5)
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+a4dt*cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          bz(k,j) = bz(k,j)*fact
          cz(k,j) = cz(k,j)*fact + 1.
          dz(k,j) = dz(k,j)*fact
          fz(k,j) = rhs(i,j,k)
  622     continue
        else if(ntorder .eq. 3) then
         if(itime.ge.3) then
           dt0 = dtphy
           dt1 = dtphy
           dt2 = dtphy
           a4 = (dt1*dt2+2*dt0*dt2+dt1**2+4*dt0*dt1+3*dt0**2)/
     .          (dt0*(dt1+dt0)*(dt2+dt1+dt0))
c  
           if(ibdf2opt.eq.1) then
             bdfoptc = -0.58d0/3.0d0
             a4  = ( 1.5d0 - bdfoptc)/dtphy
           endif
         else
cvn        temp fix for variable dt
           dt0 = dtphy
           dt1 = dtphy
           dt2 = dtphy
c          dt1 = 0.015625d0
c          dt2 = 0.015625d0
           if(itime.eq.2) dt1 = dtphy
c
           a4 = (dt1*dt2+2*dt0*dt2+dt1**2+4*dt0*dt1+3*dt0**2)/
     .          (dt0*(dt1+dt0)*(dt2+dt1+dt0))
c
         endif
c
c     scale a4 with dtphy
c
          a4dt = a4*dtphy
          do 623 k=kmn+1,kmx
c         fact  = 6.*cfltrb*dtl(i,j,k)*dtphy/( 6.*dtphy +
c    .                                        11.*cfltrb*dtl(i,j,k) )
          fact=cfltrb*dtl(i,j,k)*dtphy/(dtphy+a4dt*cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          bz(k,j) = bz(k,j)*fact
          cz(k,j) = cz(k,j)*fact + 1.
          dz(k,j) = dz(k,j)*fact
          fz(k,j) = rhs(i,j,k)
  623     continue
c
        end if
c
  625 continue
c
c     solve the scalar tridiagonal equations in k-direction
c
      call tridag(mx1dwk,mx1dwk,kmn+1,kmx,jmn+1,jmx,wk2dz,bz,cz,dz,fz)
c
c     update the r.h.s. for A-F scheme
c
      do 630 j=jmn+1,jmx
      do 630 k=kmn+1,kmx
      rhs(i,j,k)   = fz(k,j)
  630 continue
c
  610 continue
c
  611 continue
c
c
c     smooth corrections through residual smoothing
c
cvn   call psmoot (imn,  jmn,  kmn,  imx,  jmx,kmx,
cvn  .             imp1,jmp1,kmp1,imp2,jmp2,kmp2,
cvn  .             wk2dx,wk2dy,rhs,0.3,0.3,0.3)
c
c     update turv1 ( spalart's "nu~")
c
      sumn     = 0.
      negn     = 0
      do 700 k=kmn+1,kmx
      do 700 j=jmn+1,jmx
      do 700 i=imn+1,imx
      temp     = rhs(i,j,k) + turv1(i,j,k)
      turv1(i,j,k) = max(temp,1.e-12)
      rhs(i,j,k)   = max(rhs(i,j,k),0.)
      incrneg      = (1 - ifix (sign(1.,temp)) )/2
      negn         = negn + incrneg
      sumn         = sumn + (rhs(i,j,k))**2
 700  continue
c
c     write the residual history for this equation  on unit 92
c
cBWW  idebug       = 1
      idebug       = 1
      if(idebug.eq.1) then
        sumn       = sqrt(sumn)/real((kmx-kmn)*(jmx-jmn)*(imx-imn))
        write (91,' (" iteration,ibloc, log res, negn",
     .            /,2i5,e15.5,i5)') iter,ibloc,sumn, negn
c    .            /,i5,e15.5,i5)') iter,ibloc, alog10(sumn), negn
      endif
c
 1000 continue
c
c     compute eddy-viscosity from turv1 
c  
      do 800 k=kmn+1,kmx
      do 800 j=jmn+1,jmx
      do 800 i=imn+1,imx
c
      chi     = turv1(i,j,k)/anu(i,j,k)
      fv1     = chi**3/(chi**3 + cv1**3)
      eomu(i,j,k)  = fv1*turv1(i,j,k)*w(i,j,k,1)/amu(i,j,k)
 800  continue
c
      idebug = 0
      if(idebug.gt.0) then
      if(init.gt.5.and.init.le.10) then
        do 810 k=kmn+1,kmx,4
        do 810 i=imn+1,imx,8
        do 810 j=jmn+1,jmx
        write(91,901) k,i,j,amu(i,j,k),eomu(i,j,k),turv1(i,j,k),
     .                fv(i,j,k,1),fv(i,j,k,2),x(i,j,k,1),
     .                x(i,j,k,2),x(i,j,k,3)
 810    continue
 901    format(2x,'k,i,j,amu,eomu,turv1,vor,smin',3i5,8e12.3)
      endif
      endif
      init = init + 1
      return
      end
