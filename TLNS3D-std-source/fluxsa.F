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
      subroutine fluxsa(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2,w,p,eomu,turv1,vol,x,
     .                  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                  bx,cx,dx,by,cy,dy,bz,cz,dz,
     .                  fx,fy,fz,wk2dx,wk2dy,wk2dz,
     .                  fv, dtl, dti, dtj, dtk,dtvi,dtvj,dtvk, 
     .                  mx1dwk,turv1s, anu, amu, rhs, rhsimp,
     .                  wk1da,wk1db,mstrt,smin ,ibloc,iter1,
     .                  itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,
     .                  ntorder,tv1old,dtphy ,
     .                  irkstg,nrkstg,mrkstg,
     .                  rkphysa,rkphyse,dtv1ork)
c#ifdef  RKPHY
c# else
c    .                  ntorder,tv1old,dtphy )
c# endif
c----------------------------------------------------------------------
c 
c 
c       object : compute fluxes for the 1-equation Spalart's
c                turbulence model and save it in the dtv1ork array
c       initial coding : by v.n.vatsa of nasa langley (Freb. 2001)
c                        last update (Nov. 2003)
c 
c 
c

      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
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
      if(irkstg.lt.1) return
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
c
c     start outer iteration loop for the 1-equation model
c
cvn   do 1000 iter = 1,nsubit
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
        gg    = max(gg,tol)
c
c       fw    = ( (gg**(-6)+cw3**(-6))/(1.+cw3**(-6)) )**(-1./6.)
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
c#ifdef  RKPHY
c       save the rhs in dtv1ork array
c
        dtv1ork(i,j,k,irkstg) = rhs(i,j,k)
c
c#endif
  355   continue
  350 continue
c
      return
      end
