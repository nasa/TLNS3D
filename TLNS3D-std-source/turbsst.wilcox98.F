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
      subroutine turbsst(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2,w,p,eomu,turv1,turv2,vol,x,
     .                  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                  bx,cx,dx,by,cy,dy,bz,cz,dz,
     .                  fx,fy,fz,wk2dx,wk2dy,wk2dz,
     .                  fv, dtl, dti, dtj, dtk,dtvi,dtvj,dtvk,mx1dwk,
     .                  anu, amu, damp1, damp2,
     .                  wk1da,wk1db,mstrt,smin ,ibloc,iter1,
     .                  itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,
#ifdef  RKPHY
     .                  ntorder,tv1old,tv2old,dtphy ,
     .                  irkstg,nrkstg,mrkstg,
     .                  rkphysa,rkphyse,dtv1ork,dtv2ork )
# else
     .                  ntorder,tv1old,tv2old,dtphy )
# endif
c----------------------------------------------------------------------
c 
c 
c       object : solve the 2-equation Menter's turbulence model
c       initial coding : by v.n.vatsa of nasa langley (April, 1993)
c       modified       : (Sep. 1993) doesnot use Menter's suggestions
c                        about c-d term in damping
c                        However uses new limit on Pk
c                        (ala AIAA Paper 93-2906)
c                        rest of the coefficients are from original paper
c
c      latest mods: make it compatible with CFL3D and Menter's suggestions
c                   on cross-diffusion terms (June 2002)
c 
c       turv1 is k
c       turv2 is omega
c 
c 
c

      common/trbfs/ anutinf,zkinf,ominf
      common/timek/  ct(9)
      common/vis   / diffac
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension 
     .    w(imp2,jmp2,kmp2,5),p(imp2,jmp2,kmp2),vol(imp2,jmp2,kmp2),
     .    x(imp1,jmp1,kmp1,3),fv(imp2,jmp2,kmp2,4),
     . turv1(imp2,jmp2,kmp2),turv2(imp2,jmp2,kmp2),
     .  anu(imp2,jmp2,kmp2),amu(imp2,jmp2,kmp2),eomu(imp2,jmp2,kmp2),
     .  six(imp1,jmp1,kmp1),siy(imp1,jmp1,kmp1),siz(imp1,jmp1,kmp1),
     .  sjx(imp1,jmp1,kmp1),sjy(imp1,jmp1,kmp1),sjz(imp1,jmp1,kmp1),
     .  skx(imp1,jmp1,kmp1),sky(imp1,jmp1,kmp1),skz(imp1,jmp1,kmp1),
     .  dti(imp2,jmp2,kmp2),dtj(imp2,jmp2,kmp2),dtk(imp2,jmp2,kmp2),
     .  dtvi(imp2,jmp2,kmp2),dtvj(imp2,jmp2,kmp2),dtvk(imp2,jmp2,kmp2),
     .  dtl(imp2,jmp2,kmp2),smin(imp2,jmp2,kmp2),
     .  damp1(imp2,jmp2,kmp2),damp2(imp2,jmp2,kmp2)
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
      dimension     tv1old(imp2,jmp2,kmp2,ntorder),
     .              tv2old(imp2,jmp2,kmp2,ntorder)
c
#ifdef  RKPHY
c
      dimension     rkphysa(mrkstg,mrkstg),rkphyse(mrkstg)
c
      dimension     dtv1ork (imp2,jmp2,kmp2,nrkstg)
      dimension     dtv2ork (imp2,jmp2,kmp2,nrkstg)
#endif
c
c     data init/0/, cfltrb/50./
c     data init/0/, cfltrb/25./
c     data init/0/, cfltrb/10./
      data init/0/, cfltrb/5./
cwilcox   data init/0/, cfltrb/5./
      iwilcox = 1
c     iwilcox = 0
c
      if (iprec.gt.1) then
        ipy    = 1
        ipn    = 0
      else
        ipy    = 0
        ipn    = 1
      endif
      anin     = 9.*real(ipy)
c        does not account for physical time derivative in Delta tau
      anin     = 0.
c
c     fv(i,j,k,1) : used for vorticity
c     fv(i,j,k,2) : used to store turbulent-vioscosity locally
c     fv(i,j,k,3) : used to store k values during iterative process
c     fv(i,j,k,4) : used to store omega values during iterative process
c
c     damp1(i,j,k): used to store cross-diffusion term of Menter
c     damp2(i,j,k): used to store blending function F1
c
c
c     Set up universal constants
c
      vk      = 0.41
cold  a1      = 0.3
      a1      = 0.31
c
c     set up constants for set 1
c
      beta1   = 0.075
      betstr1 = 0.09
      sigk1   = 0.85
      if(iwilcox.eq.1) sigk1 = 0.5
c   ***rumseyWilcox98
      if(iwilcox.eq.1) beta1   = 0.072
c   *********
ccold sigo1   = 0.65
cnew  sigo1   = 0.50
      sigo1   = 0.65
      if(iwilcox.eq.1) sigo1 = 0.5
      alp1    = beta1/betstr1 - sigo1*vk*vk/sqrt(betstr1)
c
c     set up constants for set 2
c
      beta2   = 0.0828
      betstr2 = 0.09
      sigk2   = 1.00
      sigo2   = 0.856
      alp2    = beta2/betstr2 - sigo2*vk*vk/sqrt(betstr2)
c
c     set up some other convenient groupings
c
      sfac = sqrt(gamma)*rm/rey
c
c Set number of subiterations to solve turbulence field eqn per iteration
c (usually, 1 is sufficient... but if residual diverges then may need more)
c NOTE:  Extra subiterations are expensive! (approx 3.e-4 sec/gridpoint
c per extra subiteration on CRAY YMP)
c
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
      ai2m     = real(1-i2dfl)

      tfac = 1.

      do 10 k=kmn,kmp1
      do 10 j=jmn,jmp1
      do 10 i=imn,imp1
c
c
c
c     dtl(i,j,k)  = vol(i,j,k)*dtli*dtlj*dtlk/
c    .             (dtli*dtlj  +dtlj*dtlk  +dtlk*dtli)
c
c     dtl(i,j,k)  = (1.-i2dfl)*dtl(i,j,k) + i2dfl*vol(i,j,k)*
c    .                                      dtli*dtlj/(dtli+dtlj)
c
c     dtl(i,j,k)  = vol(i,j,k)/(dtli + dtlj + ai2m*dtlk)
        dtli       =  dti(i,j,k) + tfac*diffac*dtvi(i,j,k)
        dtlj       =  dtj(i,j,k) + tfac*diffac*dtvj(i,j,k)
        dtlk       =  dtk(i,j,k) + tfac*diffac*dtvk(i,j,k)
        eigt       =  anin*ct(ntorder+1)*vol(i,j,k) / dtphy
        dtl(i,j,k) =  vol(i,j,k)/(eigt + dtli + dtlj + ai2m*dtlk)
  10  continue
c
c     compute molecular viscosities and put initial values
c     of k and omega in fv(..3) and fv(..4) arrays
c     put turbulent viscosity in fv(..2)
c
c
      visc1     = suthc/tref
c
      do 20 k=kmn,kmp1
      do 20 j=jmn,jmp1
      do 15 i=imn,imp1
      fv(i,j,k,3) = turv1(i,j,k)
      fv(i,j,k,4) = turv2(i,j,k)
      ts        =  p(i,j,k)/w(i,j,k,1) 
      amu(i,j,k)= ts
   15 continue
c
      if(suthc.gt.0.) then
c
        do 16 i=imn,imp1
        ts        =  p(i,j,k)/w(i,j,k,1) 
        amu(i,j,k)=((1.+visc1)/(ts+visc1))*ts**1.5
   16   continue
      endif
c
      do 20 i=imn,imp1
      fv(i,j,k,2) = amu(i,j,k)*eomu(i,j,k)
      anu(i,j,k)  = amu(i,j,k)/w(i,j,k,1)
   20 continue
c
c
c     start outer iteration loop for the 2-equation model
c
      do 1000 iter = 1,nsubit
c
c     compute turbulent viscosity : fv(i,j,k,2) at interior cells
c     using current values of k and omega
c
c     go to 106
c     if(iter.gt.1) go to 106
      do 105 k=kmn+1,kmx
      do 105 j=jmn+1,jmx
      tran    = 1.0
c     if(ibloc.eq.4. and.(i-imn).le.4) tran = 0.
c     if(ibloc.eq.1) tran = 0.
      if(iwilcox.eq.1) then
        do i=imn+1,imx
          fv(i,j,k,2)  = w(i,j,k,1)*turv1(i,j,k)/turv2(i,j,k)
          fv(i,j,k,2)  = min(fv(i,j,k,2),100000.)
          enddo
      else
c
        do i=imn+1,imx
c
c     if(init.eq.0) write(6,901) i,j,k,turv1(i,j,k),turv2(i,j,k),
c    .                           smin(i,j,k),sfac
c901  format(2x,'t1,t2,s',3i4,4e12.5)
      denom1  = 2.*sqrt(turv1(i,j,k))/(.09*turv2(i,j,k)*smin(i,j,k))
     .            *sfac
      denom2  = 500.*amu(i,j,k)*sfac*sfac/(w(i,j,k,1)*turv2(i,j,k)*
     .                                     smin(i,j,k)*smin(i,j,k))
      arg2    = max(denom1,denom2)
      f2      = tanh(arg2*arg2)
      denom   = max(a1*turv2(i,j,k),sfac*fv(i,j,k,1)*f2)
      fv(i,j,k,2) = tran*a1*w(i,j,k,1)*turv1(i,j,k)/denom
      fv(i,j,k,2) = min(fv(i,j,k,2),100000.)
        enddo
      endif
 105  continue
 106  continue
c
c   ***rumseyWilcox98
c     compute damp1 = CD = cross-derivative terms for SST or Wilcox98 model
c   *********
c
      do k=kmn+1,kmx
        do j=jmn+1,jmx
          do i=imn+1,imx
c
c     xi-derivatives
c
          xa    = 0.5*(six(i,j,k)+six(i-1,j,k))/vol(i,j,k)
          ya    = 0.5*(siy(i,j,k)+siy(i-1,j,k))/vol(i,j,k)
          za    = 0.5*(siz(i,j,k)+siz(i-1,j,k))/vol(i,j,k)
          tt    = xa*xa + ya*ya + za*za
          caa   = 2.0*sigo2*tt*sfac/turv2(i,j,k)
          damp1(i,j,k) = 0.25*caa*(turv1(i+1,j,k)-turv1(i-1,j,k))*
     .                            (turv2(i+1,j,k)-turv2(i-1,j,k))
c
c     eta derivatives
c
          xa    = 0.5*(sjx(i,j,k)+sjx(i,j-1,k))/vol(i,j,k)
          ya    = 0.5*(sjy(i,j,k)+sjy(i,j-1,k))/vol(i,j,k)
          za    = 0.5*(sjz(i,j,k)+sjz(i,j-1,k))/vol(i,j,k)
          tt    = xa*xa + ya*ya + za*za
          caa   = 2.0*sigo2*tt*sfac/turv2(i,j,k)
          damp1(i,j,k) = damp1(i,j,k) +
     .                   0.25*caa*(turv1(i,j+1,k)-turv1(i,j-1,k))*
     .                            (turv2(i,j+1,k)-turv2(i,j-1,k))
c
c     zeta derivatives
c
          xa    = 0.5*(skx(i,j,k)+skx(i,j,k-1))/vol(i,j,k)
          ya    = 0.5*(sky(i,j,k)+sky(i,j,k-1))/vol(i,j,k)
          za    = 0.5*(skz(i,j,k)+skz(i,j,k-1))/vol(i,j,k)
          tt    = xa*xa + ya*ya + za*za
          caa   = 2.0*sigo2*tt*sfac/turv2(i,j,k)
          damp1(i,j,k) = damp1(i,j,k) +
     .                   0.25*caa*(turv1(i,j,k+1)-turv1(i,j,k-1))*
     .                            (turv2(i,j,k+1)-turv2(i,j,k-1))*
     .                            (1.-i2dfl)
c
          enddo
        enddo
      enddo
c   ***rumseyWilcox98
c     Wilcox cross derivative term includes the above, but is different:
      if (iwilcox.eq.1) then
      do k=kmn+1,kmx
        do j=jmn+1,jmx
          do i=imn+1,imx
            temp=damp1(i,j,k)*sfac/(2.*sigo2*
     +           turv2(i,j,k)**2)
            temp=max(temp,0.)
            damp1(j,k,i)=(1.+680.*temp*temp)/(1.+400.*temp*temp)
          enddo
        enddo
      enddo
      end if
c   *********
c
c     evaluate damping terms f1
c
c
      do 50 k=kmn+1,kmx
      do 50 j=jmn+1,jmx
      if(iwilcox.eq.1) then
         do i=imn+1,imx
         damp2(i,j,k)    =  1.0
         enddo
      else
        do i=imn+1,imx
c
c       evaluate F1 and store it in damp2
c
        arg1  =  sqrt(turv1(i,j,k))/(0.09*turv2(i,j,k)*
     .           smin(i,j,k)) * sfac
        arg2  =  500.*amu(i,j,k)*sfac*sfac/(w(i,j,k,1)*turv2(i,j,k)*
     .                                     smin(i,j,k)*smin(i,j,k))
        arga  =  max(arg1,arg2)
        temp  =  max(damp1(i,j,k)/sfac,1.e-20)
        argb  =  4.*sigo2*turv1(i,j,k)/(temp*smin(i,j,k)*smin(i,j,k))
        arg   =  min(arga,argb)
        damp2(i,j,k)    =  tanh(arg*arg*arg*arg)
        enddo
      endif

  50    continue
c
c     Evaluate  Explicit part of F-xi-xi viscous terms
c
      do 110 k=kmn+1,kmx
      do 110 j=jmn+1,jmx
c
      do 115 i=imn+1,imx
c
        iu    =  min(i+1,imx)
        il    =  max(i-1,imn+1)
        dfacep=  0.5*(damp2(i,j,k)+damp2(iu,j,k) )
        dfacem=  0.5*(damp2(i,j,k)+damp2(il,j,k) )
        sigkp =  dfacep*sigk1 + (1.-dfacep)*sigk2
        sigkm =  dfacem*sigk1 + (1.-dfacem)*sigk2
        sigop =  dfacep*sigo1 + (1.-dfacep)*sigo2
        sigom =  dfacem*sigo1 + (1.-dfacem)*sigo2
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
        anutp = 0.5*(fv(i,j,k,2) + fv(i+1,j,k,2))
        anutm = 0.5*(fv(i,j,k,2) + fv(i-1,j,k,2))
        anup  = 0.5*(amu  (i,j,k) + amu  (i+1,j,k))
        anum  = 0.5*(amu  (i,j,k) + amu  (i-1,j,k))
        cdp   = ( anup+sigop*anutp ) * ttp *sfac/w(i,j,k,1)
        cdm   = ( anum+sigom*anutm ) * ttm *sfac/w(i,j,k,1)
c   extra terms in sst model
c       zkp   = 0.5*(turv1(i,j,k) + turv1(i+1,j,k))
c       zkm   = 0.5*(turv1(i,j,k) + turv1(i-1,j,k))
c       cap   = 2.*(1.-f1)*sigo2*zkp*ttp*sfac/turv2(i,j,k)
c       cam   = 2.*(1.-f1)*sigo2*zkm*ttm*sfac/turv2(i,j,k)
c
c       bxxo  =  -cdm - cam
c       cxxo  =   cdp + cap + cdm + cam
c       dxxo  =  -cdp - cap
        bxxo  =  -cdm
        cxxo  =   cdp + cdm
        dxxo  =  -cdp 
c
        cdp   = ( anup+sigkp*anutp ) * ttp *sfac/w(i,j,k,1)
        cdm   = ( anum+sigkm*anutm ) * ttm *sfac/w(i,j,k,1)
        bxxk  =  -cdm
        cxxk  =   cdp + cdm
        dxxk  =  -cdp
        fv(i,j,k,3)  = -bxxk*turv1(i-1,j,k) -cxxk*turv1(i,j,k)
     .                 -dxxk*turv1(i+1,j,k)
        fv(i,j,k,4)  = -bxxo*turv2(i-1,j,k) -cxxo*turv2(i,j,k)
     .                 -dxxo*turv2(i+1,j,k)
c
c     Evaluate  Explicit part of F-xi advective terms
c
        uu    = (xa*w(i,j,k,2) + ya*w(i,j,k,3) + za*w(i,j,k,4))/
     .              w(i,j,k,1)
        sgnu  = sign(1., uu)
        app   = 0.5*(1. + sgnu)
        apm   = 0.5*(1. - sgnu)
        fv(i,j,k,3)  = fv(i,j,k,3)
     .               - uu*(app*(turv1(i  ,j,k) - turv1(i-1,j,k))+
     .                     apm*(turv1(i+1,j,k) - turv1(i  ,j,k)) )
        fv(i,j,k,4)  = fv(i,j,k,4)
     .               - uu*(app*(turv2(i  ,j,k) - turv2(i-1,j,k))+
     .                     apm*(turv2(i+1,j,k) - turv2(i  ,j,k)) )
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
        ju    =  min(j+1,jmx)
        jl    =  max(j-1,jmn+1)
        dfacep=  0.5*(damp2(i,j,k)+damp2(i,ju,k) )
        dfacem=  0.5*(damp2(i,j,k)+damp2(i,jl,k) )
        sigkp =  dfacep*sigk1 + (1.-dfacep)*sigk2
        sigkm =  dfacem*sigk1 + (1.-dfacem)*sigk2
        sigop =  dfacep*sigo1 + (1.-dfacep)*sigo2
        sigom =  dfacem*sigo1 + (1.-dfacem)*sigo2
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
        anutp = 0.5*(fv(i,j,k,2) + fv(i,j+1,k,2))
        anutm = 0.5*(fv(i,j,k,2) + fv(i,j-1,k,2))
        anup  = 0.5*(amu  (i,j,k) + amu  (i,j+1,k))
        anum  = 0.5*(amu  (i,j,k) + amu  (i,j-1,k))
        cdp   = ( anup+sigop*anutp ) * ttp *sfac/w(i,j,k,1)
        cdm   = ( anum+sigom*anutm ) * ttm *sfac/w(i,j,k,1)
c   extra terms in sst model
c       zkp   = 0.5*(turv1(i,j,k) + turv1(i,j+1,k))
c       zkm   = 0.5*(turv1(i,j,k) + turv1(i,j-1,k))
c       cap   = 2.*(1.-f1)*sigo2*zkp*ttp*sfac/turv2(i,j,k)
c       cam   = 2.*(1.-f1)*sigo2*zkm*ttm*sfac/turv2(i,j,k)
c
c       byyo  =  -cdm - cam
c       cyyo  =   cdp + cap + cdm + cam
c       dyyo  =  -cdp - cap
        byyo  =  -cdm
        cyyo  =   cdp + cdm 
        dyyo  =  -cdp 
c
        cdp   = ( anup+sigkp*anutp ) * ttp *sfac/w(i,j,k,1)
        cdm   = ( anum+sigkm*anutm ) * ttm *sfac/w(i,j,k,1)
        byyk  =  -cdm
        cyyk  =   cdp + cdm
        dyyk  =  -cdp
        fv(i,j,k,3)  =  fv(i,j,k,3)
     .                 -byyk*turv1(i,j-1,k) -cyyk*turv1(i,j,k)
     .                 -dyyk*turv1(i,j+1,k)
        fv(i,j,k,4)  =  fv(i,j,k,4)
     .                 -byyo*turv2(i,j-1,k) -cyyo*turv2(i,j,k)
     .                 -dyyo*turv2(i,j+1,k)
c
c     Evaluate  Explicit part of F-eta advective terms
c
        uu    = (xa*w(i,j,k,2) + ya*w(i,j,k,3) + za*w(i,j,k,4))/
     .              w(i,j,k,1)
        sgnu  = sign(1., uu)
        app   = 0.5*(1. + sgnu)
        apm   = 0.5*(1. - sgnu)
        fv(i,j,k,3)  = fv(i,j,k,3)
     .               - uu*(app*(turv1(i  ,j,k) - turv1(i,j-1,k))+
     .                     apm*(turv1(i,j+1,k) - turv1(i  ,j,k)) )
        fv(i,j,k,4)  = fv(i,j,k,4)
     .               - uu*(app*(turv2(i  ,j,k) - turv2(i,j-1,k))+
     .                     apm*(turv2(i,j+1,k) - turv2(i  ,j,k)) )
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
        ku    =  min(k+1,kmx)
        kl    =  max(k-1,kmn+1)
        dfacep=  0.5*(damp2(i,j,k)+damp2(i,j,ku) )
        dfacem=  0.5*(damp2(i,j,k)+damp2(i,j,kl) )
        sigkp =  dfacep*sigk1 + (1.-dfacep)*sigk2
        sigkm =  dfacem*sigk1 + (1.-dfacem)*sigk2
        sigop =  dfacep*sigo1 + (1.-dfacep)*sigo2
        sigom =  dfacem*sigo1 + (1.-dfacem)*sigo2
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
c
        anutp = 0.5*(fv(i,j,k,2) + fv(i,j,k+1,2))
        anutm = 0.5*(fv(i,j,k,2) + fv(i,j,k-1,2))
        anup  = 0.5*(amu  (i,j,k) + amu  (i,j,k+1))
        anum  = 0.5*(amu  (i,j,k) + amu  (i,j,k-1))
        cdp   = ( anup+sigop*anutp ) * ttp *sfac/w(i,j,k,1)
        cdm   = ( anum+sigom*anutm ) * ttm *sfac/w(i,j,k,1)
c   extra terms in sst model
c       zkp   = 0.5*(turv1(i,j,k) + turv1(i,j,k+1))
c       zkm   = 0.5*(turv1(i,j,k) + turv1(i,j,k-1))
c       cap   = 2.*(1.-f1)*sigo2*zkp*ttp*sfac/turv2(i,j,k)
c       cam   = 2.*(1.-f1)*sigo2*zkm*ttm*sfac/turv2(i,j,k)
c
c       bzzo  =  -cdm - cam
c       czzo  =   cdp + cap + cdm + cam
c       dzzo  =  -cdp - cap
        bzzo  =  -cdm 
        czzo  =   cdp + cdm
        dzzo  =  -cdp 
c
        cdp   = ( anup+sigkp*anutp ) * ttp *sfac/w(i,j,k,1)
        cdm   = ( anum+sigkm*anutm ) * ttm *sfac/w(i,j,k,1)
        bzzk  =  -cdm
        czzk  =   cdp + cdm
        dzzk  =  -cdp
        fv(i,j,k,3)  =  fv(i,j,k,3)
     .                 -bzzk*turv1(i,j,k-1) -czzk*turv1(i,j,k)
     .                 -dzzk*turv1(i,j,k+1)
        fv(i,j,k,4)  =  fv(i,j,k,4)
     .                 -bzzo*turv2(i,j,k-1) -czzo*turv2(i,j,k)
     .                 -dzzo*turv2(i,j,k+1)
c
c     Evaluate  Explicit part of F-zeta advective terms
c
        uu    = (xa*w(i,j,k,2) + ya*w(i,j,k,3) + za*w(i,j,k,4))/
     .              w(i,j,k,1)
        sgnu  = sign(1., uu)
        app   = 0.5*(1. + sgnu)
        apm   = 0.5*(1. - sgnu)
        fv(i,j,k,3) = fv(i,j,k,3)
     .              - uu*(app*(turv1(i  ,j,k) - turv1(i,j,k-1))+
     .                    apm*(turv1(i,j,k+1) - turv1(i  ,j,k)) )
        fv(i,j,k,4) = fv(i,j,k,4)
     .              - uu*(app*(turv2(i  ,j,k) - turv2(i,j,k-1))+ 
     .                    apm*(turv2(i,j,k+1) - turv2(i  ,j,k)) )
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
c     transition via Barth
c     presently set to fully turbulent
c
      cutoff  = 1.
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
c   ***rumseyWilcox98
      if (iwilcox.eq.1) then
c   Wilcox98 source terms
      do 354 j=jmn+1,jmx
c
        ss     = fv(i,j,k,1)
c
c       add source terms to r.h.s
c
        pk          = fv(i,j,k,2)*sfac*ss*ss/w(i,j,k,1)
        dk          = damp1(i,j,k)*betstr1*turv1(i,j,k)
     .              * turv2(i,j,k)/sfac
        pk          = min(pk,(20.*dk))
        fv(i,j,k,3) = fv(i,j,k,3)
     .              + cutoff*wk1da(j)*pk - dk
        fv(i,j,k,4) = fv(i,j,k,4)
     .              + cutoff*wk1da(j)*alp1*sfac*ss*ss
     .              - beta1 *turv2(i,j,k)*turv2(i,j,k)/sfac
  354   continue
      else
c   *********
c   SST source terms
      do 355 j=jmn+1,jmx
c
        ss     = fv(i,j,k,1)
c
        f1     =  damp2(i,j,k)
        betax  =  f1*beta1   + (1.-f1)*beta2
        betstr =  f1*betstr1 + (1.-f1)*betstr2
        alp    =  f1*alp1    + (1.-f1)*alp2
c
c       add source terms to r.h.s
c
        pk          = fv(i,j,k,2)*sfac*ss*ss/w(i,j,k,1)
        dk          = betstr*turv1(i,j,k)*turv2(i,j,k)/sfac
        pk          = min(pk,(20.*dk))
        fv(i,j,k,3) = fv(i,j,k,3)   
     .              + cutoff*wk1da(j)*pk - dk
        fv(i,j,k,4) = fv(i,j,k,4)   
     .              + cutoff*wk1da(j)*alp*sfac*ss*ss
     .              - betax *turv2(i,j,k)*turv2(i,j,k)/sfac
     .              + (1.-f1)*damp1(i,j,k)
  355   continue
c   ***rumseyWilcox98
      end if
c   *********
c
  350 continue
c        write(6,'(2x,"turbulent zone in i,j,k-direction",6i4,
c    .                 2x,"dimension of block are",3i5,2f10.3)')
c    .                 itrb1,itrb2,
c    .                 jtrb1,jtrb2,
c    .                 ktrb1,ktrb2,
c    .                 imx,jmx,kmx,cutoff,wk1da(jmn+1)
c
c
c     Evaluate  Implicit part of F-eta-eta viscous terms
c
      do 510 k=kmn+1,kmx
      do 525 i=imn+1,imx
c
      do 515 j=jmn+1,jmx
c
        ju    =  min(j+1,jmx)
        jl    =  max(j-1,jmn+1)
        dfacep=  0.5*(damp2(i,j,k)+damp2(i,ju,k) )
        dfacem=  0.5*(damp2(i,j,k)+damp2(i,jl,k) )
        sigkp =  dfacep*sigk1 + (1.-dfacep)*sigk2
        sigkm =  dfacem*sigk1 + (1.-dfacem)*sigk2
        sigop =  dfacep*sigo1 + (1.-dfacep)*sigo2
        sigom =  dfacem*sigo1 + (1.-dfacem)*sigo2
        f1    =  damp2(i,j,k)
        sigk  =  f1*sigk1 + (1.-f1)*sigk2
        sigo  =  f1*sigo1 + (1.-f1)*sigo2
        betax  =  f1*beta1   + (1.-f1)*beta2
        betstr =  f1*betstr1 + (1.-f1)*betstr2
        alp    =  f1*alp1    + (1.-f1)*alp2
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
        anutp = 0.5*(fv(i,j,k,2) + fv(i,j+1,k,2))
        anutm = 0.5*(fv(i,j,k,2) + fv(i,j-1,k,2))
        anup  = 0.5*(amu  (i,j,k) + amu  (i,j+1,k))
        anum  = 0.5*(amu  (i,j,k) + amu  (i,j-1,k))
        cdp   = ( anup+sigop*anutp ) * ttp *sfac/w(i,j,k,1)
        cdm   = ( anum+sigom*anutm ) * ttm *sfac/w(i,j,k,1)
c   extra terms in sst model
c       zkp   = 0.5*(turv1(i,j,k) + turv1(i,j+1,k))
c       zkm   = 0.5*(turv1(i,j,k) + turv1(i,j-1,k))
c       cap   = 2.*(1.-f1)*sigo2*zkp*ttp*sfac/turv2(i,j,k)
c       cam   = 2.*(1.-f1)*sigo2*zkm*ttm*sfac/turv2(i,j,k)
c
c  coefficients for omega-equation
        bz(j,i)  =  -cdm 
        cz(j,i)  =   cdp + cdm 
        dz(j,i)  =  -cdp 
c  coefficients for k-equation
        cdp   = ( anup+sigkp*anutp ) * ttp *sfac/w(i,j,k,1)
        cdm   = ( anum+sigkm*anutm ) * ttm *sfac/w(i,j,k,1)
        by(j,i)  =  -cdm
        cy(j,i)  =   cdp + cdm
        dy(j,i)  =  -cdp
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
        bz(j,i)   = bz(j,i) - uu*app
        cz(j,i)   = cz(j,i) + uu*(app-apm)
        dz(j,i)   = dz(j,i) + uu*apm
c
c       add part of the source term to the diagonal in LHS
c
c       if(iwilcox.eq.1) then
c         cy(j,i)   = cy(j,i)
c    .              + betstr1*turv2(i,j,k)/sfac
c         cz(j,i)   = cz(j,i)
c    .              + 2.*beta1*turv2(i,j,k)/sfac
c       else

c   ***rumseyWilcox98
        if(iwilcox.eq.1) then
          cy(j,i)   = cy(j,i)
     .              + damp1(i,j,k)*betstr*turv2(i,j,k)/sfac
        else
c   *********
          cy(j,i)   = cy(j,i)
     .              + betstr*turv2(i,j,k)/sfac
c   ***rumseyWilcox98
        end if
c   *********
          cz(j,i)   = cz(j,i)
     .              + 2.*betax*turv2(i,j,k)/sfac
        if(iwilcox.eq.0)  
     .    cz(j,i)   = cz(j,i)
     .              + abs(damp1(i,j,k))/turv2(i,j,k)
c       endif
  515   continue
c
      if(ntorder.eq.0) then
        do 520 j=jmn+1,jmx
        fact  = cfltrb*dtl(i,j,k)
        fact  = min(fact,100.)
        by(j,i)   = by(j,i)*fact
        cy(j,i)   = cy(j,i)*fact + 1.
        dy(j,i)   = dy(j,i)*fact
        fy(j,i)   = fv(i,j,k,3)*fact
c
        bz(j,i)   = bz(j,i)*fact
        cz(j,i)   = cz(j,i)*fact + 1.
        dz(j,i)   = dz(j,i)*fact
        fz(j,i)   = fv(i,j,k,4)*fact
c
  520   continue
        else if(ntorder .eq. 1) then
#ifdef  RKPHY
        if(nrkstg.gt.2) then
          do 519 j=jmn+1,jmx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          rhsold1= 0.
          if(irkstg.gt.1) then
            do lstg = 1,irkstg-1
            rhsold1 = rhsold1 + rkphysa(irkstg,lstg)*dtv1ork(i,j,k,lstg)
            enddo
          endif
          by(j,i) = by(j,i)*fact*rkphysa(irkstg,irkstg)
          cy(j,i) = cy(j,i)*fact*rkphysa(irkstg,irkstg) + 1.
          dy(j,i) = dy(j,i)*fact*rkphysa(irkstg,irkstg)
          fy(j,i) = fv(i,j,k,3)*fact*rkphysa(irkstg,irkstg)
     .            -  fact/dtphy*(turv1(i,j,k)-tv1old(i,j,k,1))
     .            +  fact*rhsold1
c
          rhsold2= 0.
          if(irkstg.gt.1) then
            do lstg = 1,irkstg-1
            rhsold2 = rhsold2 + rkphysa(irkstg,lstg)*dtv2ork(i,j,k,lstg)
            enddo
          endif
          bz(j,i) = bz(j,i)*fact*rkphysa(irkstg,irkstg)
          cz(j,i) = cz(j,i)*fact*rkphysa(irkstg,irkstg) + 1.
          dz(j,i) = dz(j,i)*fact*rkphysa(irkstg,irkstg)
          fz(j,i) = fv(i,j,k,4)*fact*rkphysa(irkstg,irkstg)
     .            -  fact/dtphy*(turv2(i,j,k)-tv2old(i,j,k,1))
     .            +  fact*rhsold2
  519     continue
        endif
# else
         do 521 j=jmn+1,jmx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          by(j,i) = by(j,i)*fact
          cy(j,i) = cy(j,i)*fact + 1.
          dy(j,i) = dy(j,i)*fact
          fy(j,i) = fv(i,j,k,3)*fact
     .              -fact/dtphy*(turv1(i,j,k)-tv1old(i,j,k,1))
c
          bz(j,i) = bz(j,i)*fact
          cz(j,i) = cz(j,i)*fact + 1.
          dz(j,i) = dz(j,i)*fact
          fz(j,i) = fv(i,j,k,4)*fact
     .              -fact/dtphy*(turv2(i,j,k)-tv2old(i,j,k,1))
  521     continue
# endif
        else if(ntorder .eq. 2) then
          do 522 j=jmn+1,jmx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k)*1.5)
          fact  = min(fact,100.)
          by(j,i) = by(j,i)*fact
          cy(j,i) = cy(j,i)*fact + 1.
          dy(j,i) = dy(j,i)*fact
          fy(j,i) = fv(i,j,k,3)*fact
     .             -fact/dtphy*(1.5*( turv1(i,j,k)  -tv1old(i,j,k,1))
     .                         -0.5*(tv1old(i,j,k,1)-tv1old(i,j,k,2)))
c
          bz(j,i) = bz(j,i)*fact
          cz(j,i) = cz(j,i)*fact + 1.
          dz(j,i) = dz(j,i)*fact
          fz(j,i) = fv(i,j,k,4)*fact
     .             -fact/dtphy*(1.5*( turv2(i,j,k)  -tv2old(i,j,k,1))
     .                         -0.5*(tv2old(i,j,k,1)-tv2old(i,j,k,2)))
  522     continue
        end if
c
  525 continue
c
c     solve the scalar tridiagonal equations in j-direction
c
      call tridag(mx1dwk,mx1dwk,jmn+1,jmx,imn+1,imx,wk2dy,by,cy,dy,fy)
      call tridag(mx1dwk,mx1dwk,jmn+1,jmx,imn+1,imx,wk2dz,bz,cz,dz,fz)
c
c     update the r.h.s. for A-F scheme
c
      do 530 i=imn+1,imx
      do 530 j=jmn+1,jmx
      fv(i,j,k,3)   = fy(j,i)
      fv(i,j,k,4)   = fz(j,i)
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
        iu    =  min(i+1,imx)
        il    =  max(i-1,imn+1)
        dfacep=  0.5*(damp2(i,j,k)+damp2(iu,j,k) )
        dfacem=  0.5*(damp2(i,j,k)+damp2(il,j,k) )
        sigkp =  dfacep*sigk1 + (1.-dfacep)*sigk2
        sigkm =  dfacem*sigk1 + (1.-dfacem)*sigk2
        sigop =  dfacep*sigo1 + (1.-dfacep)*sigo2
        sigom =  dfacem*sigo1 + (1.-dfacem)*sigo2
        f1    =  damp2(i,j,k)
        sigk  =  f1*sigk1 + (1.-f1)*sigk2
        sigo  =  f1*sigo1 + (1.-f1)*sigo2
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
        anutp = 0.5*(fv(i,j,k,2) + fv(i+1,j,k,2))
        anutm = 0.5*(fv(i,j,k,2) + fv(i-1,j,k,2))
        anup  = 0.5*(amu  (i,j,k) + amu  (i+1,j,k))
        anum  = 0.5*(amu  (i,j,k) + amu  (i-1,j,k))
        cdp   = ( anup+sigop*anutp ) * ttp *sfac/w(i,j,k,1)
        cdm   = ( anum+sigom*anutm ) * ttm *sfac/w(i,j,k,1)
c   extra terms in sst model
c       zkp   = 0.5*(turv1(i,j,k) + turv1(i+1,j,k))
c       zkm   = 0.5*(turv1(i,j,k) + turv1(i-1,j,k))
c       cap   = 2.*(1.-f1)*sigo2*zkp*ttp*sfac/turv2(i,j,k)
c       cam   = 2.*(1.-f1)*sigo2*zkm*ttm*sfac/turv2(i,j,k)
c  coefficients for omega-equation
        bz(i,j)  =  -cdm 
        cz(i,j)  =   cdp + cdm
        dz(i,j)  =  -cdp
c  coefficients for k-equation
        cdp   = ( anup+sigkp*anutp ) * ttp *sfac/w(i,j,k,1)
        cdm   = ( anum+sigkm*anutm ) * ttm *sfac/w(i,j,k,1)
        by(i,j)  =  -cdm
        cy(i,j)  =   cdp + cdm
        dy(i,j)  =  -cdp
c
c     Evaluate  Implicit part of F-eta advective terms
c
        uu    = (xa*w(i,j,k,2) + ya*w(i,j,k,3) + za*w(i,j,k,4))/
     .              w(i,j,k,1)
        sgnu  = sign(1., uu)
        app   = 0.5*(1. + sgnu)
        apm   = 0.5*(1. - sgnu)
        by(i,j)   = by(i,j) - uu*app
        cy(i,j)   = cy(i,j) + uu*(app-apm)
        dy(i,j)   = dy(i,j) + uu*apm
c
        bz(i,j)   = bz(i,j) - uu*app
        cz(i,j)   = cz(i,j) + uu*(app-apm)
        dz(i,j)   = dz(i,j) + uu*apm
c
  415   continue
c
      if(ntorder .eq. 0) then
        do 420 i=imn+1,imx
        fact  = cfltrb*dtl(i,j,k)
        fact  = min(fact,100.)
        by(i,j)   = by(i,j)*fact
        cy(i,j)   = cy(i,j)*fact + 1.
        dy(i,j)   = dy(i,j)*fact
        fy(i,j)   = fv(i,j,k,3)
c
        bz(i,j)   = bz(i,j)*fact
        cz(i,j)   = cz(i,j)*fact + 1.
        dz(i,j)   = dz(i,j)*fact
        fz(i,j)   = fv(i,j,k,4)
  420   continue
        else if(ntorder .eq. 1) then
#ifdef  RKPHY
        if(nrkstg.gt.2) then
          do 419 i=imn+1,imx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          by(i,j) = by(i,j)*fact*rkphysa(irkstg,irkstg)
          cy(i,j) = cy(i,j)*fact*rkphysa(irkstg,irkstg) + 1.
          dy(i,j) = dy(i,j)*fact*rkphysa(irkstg,irkstg)
          fy(i,j) = fv(i,j,k,3)
c
          bz(i,j) = bz(i,j)*fact*rkphysa(irkstg,irkstg)
          cz(i,j) = cz(i,j)*fact*rkphysa(irkstg,irkstg) + 1.
          dz(i,j) = dz(i,j)*fact*rkphysa(irkstg,irkstg)
          fz(i,j) = fv(i,j,k,4)
  419     continue
        endif
# else
          do 421 i=imn+1,imx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          by(i,j) = by(i,j)*fact
          cy(i,j) = cy(i,j)*fact + 1.
          dy(i,j) = dy(i,j)*fact
          fy(i,j) = fv(i,j,k,3)
c
          bz(i,j) = bz(i,j)*fact
          cz(i,j) = cz(i,j)*fact + 1.
          dz(i,j) = dz(i,j)*fact
          fz(i,j) = fv(i,j,k,4)
  421     continue
# endif
        else if(ntorder .eq. 2) then
          do 422 i=imn+1,imx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k)*1.5)
          fact  = min(fact,100.)
          by(i,j) = by(i,j)*fact
          cy(i,j) = cy(i,j)*fact + 1.
          dy(i,j) = dy(i,j)*fact
          fy(i,j) = fv(i,j,k,3)
c
          bz(i,j) = bz(i,j)*fact
          cz(i,j) = cz(i,j)*fact + 1.
          dz(i,j) = dz(i,j)*fact
          fz(i,j) = fv(i,j,k,4)
  422     continue
        end if
c
  425 continue
c
c     solve the scalar tridiagonal equations in i-direction
c
      call tridag(mx1dwk,mx1dwk,imn+1,imx,jmn+1,jmx,wk2dy,by,cy,dy,fy)
      call tridag(mx1dwk,mx1dwk,imn+1,imx,jmn+1,jmx,wk2dz,bz,cz,dz,fz)
c
c     update the r.h.s. for A-F scheme
c
      do 430 j=jmn+1,jmx
      do 430 i=imn+1,imx
      fv(i,j,k,3)   = fy(i,j)
      fv(i,j,k,4)   = fz(i,j)
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
        ku    =  min(k+1,kmx)
        kl    =  max(k-1,kmn+1)
        dfacep=  0.5*(damp2(i,j,k)+damp2(i,j,ku) )
        dfacem=  0.5*(damp2(i,j,k)+damp2(i,j,kl) )
        sigkp =  dfacep*sigk1 + (1.-dfacep)*sigk2
        sigkm =  dfacem*sigk1 + (1.-dfacem)*sigk2
        sigop =  dfacep*sigo1 + (1.-dfacep)*sigo2
        sigom =  dfacem*sigo1 + (1.-dfacem)*sigo2
        f1    =  damp2(i,j,k)
        sigk  =  f1*sigk1 + (1.-f1)*sigk2
        sigo  =  f1*sigo1 + (1.-f1)*sigo2
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
        anutp = 0.5*(fv(i,j,k,2) + fv(i,j,k+1,2))
        anutm = 0.5*(fv(i,j,k,2) + fv(i,j,k-1,2))
        anup  = 0.5*(amu  (i,j,k) + amu  (i,j,k+1))
        anum  = 0.5*(amu  (i,j,k) + amu  (i,j,k-1))
        cdp   = ( anup+sigop*anutp ) * ttp *sfac/w(i,j,k,1)
        cdm   = ( anum+sigom*anutm ) * ttm *sfac/w(i,j,k,1)
c   extra terms in sst model
c       zkp   = 0.5*(turv1(i,j,k) + turv1(i,j,k+1))
c       zkm   = 0.5*(turv1(i,j,k) + turv1(i,j,k-1))
c       cap   = 2.*(1.-f1)*sigo2*zkp*ttp*sfac/turv2(i,j,k)
c       cam   = 2.*(1.-f1)*sigo2*zkm*ttm*sfac/turv2(i,j,k)
c
c  coefficients for omega-equation
        bz(k,j)  =  -cdm 
        cz(k,j)  =   cdp + cdm 
        dz(k,j)  =  -cdp 
c  coefficients for k-equation
        cdp   = ( anup+sigkp*anutp ) * ttp *sfac/w(i,j,k,1)
        cdm   = ( anum+sigkm*anutm ) * ttm *sfac/w(i,j,k,1)
        by(k,j)  =  -cdm
        cy(k,j)  =   cdp + cdm
        dy(k,j)  =  -cdp
c
c     Evaluate  Implicit part of F-eta advective terms
c
        uu    = (xa*w(i,j,k,2) + ya*w(i,j,k,3) + za*w(i,j,k,4))/
     .              w(i,j,k,1)
        sgnu  = sign(1., uu)
        app   = 0.5*(1. + sgnu)
        apm   = 0.5*(1. - sgnu)
        by(k,j)   = by(k,j) - uu*app
        cy(k,j)   = cy(k,j) + uu*(app-apm)
        dy(k,j)   = dy(k,j) + uu*apm
c
        bz(k,j)   = bz(k,j) - uu*app
        cz(k,j)   = cz(k,j) + uu*(app-apm)
        dz(k,j)   = dz(k,j) + uu*apm
c
c
  615   continue
c
      if(ntorder .eq. 0) then
        do 620 k=kmn+1,kmx
        fact  = cfltrb*dtl(i,j,k)
        fact  = min(fact,100.)
        by(k,j)   = by(k,j)*fact
        cy(k,j)   = cy(k,j)*fact + 1.
        dy(k,j)   = dy(k,j)*fact
        fy(k,j)   = fv(i,j,k,3)
c
        bz(k,j)   = bz(k,j)*fact
        cz(k,j)   = cz(k,j)*fact + 1.
        dz(k,j)   = dz(k,j)*fact
        fz(k,j)   = fv(i,j,k,4)
  620   continue
        else if(ntorder .eq. 1) then
#ifdef  RKPHY
        if(nrkstg.gt.2) then
          do 619 k=kmn+1,kmx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          by(k,j) = by(k,j)*fact*rkphysa(irkstg,irkstg)
          cy(k,j) = cy(k,j)*fact*rkphysa(irkstg,irkstg) + 1.
          dy(k,j) = dy(k,j)*fact*rkphysa(irkstg,irkstg)
          fy(k,j) = fv(i,j,k,3)
c
          bz(k,j) = bz(k,j)*fact*rkphysa(irkstg,irkstg)
          cz(k,j) = cz(k,j)*fact*rkphysa(irkstg,irkstg) + 1.
          dz(k,j) = dz(k,j)*fact*rkphysa(irkstg,irkstg)
          fz(k,j) = fv(i,j,k,4)
  619     continue
        endif
# else
          do 621 k=kmn+1,kmx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k))
          fact  = min(fact,100.)
          by(k,j) = by(k,j)*fact
          cy(k,j) = cy(k,j)*fact + 1.
          dy(k,j) = dy(k,j)*fact
          fy(k,j) = fv(i,j,k,3)
c
          bz(k,j) = bz(k,j)*fact
          cz(k,j) = cz(k,j)*fact + 1.
          dz(k,j) = dz(k,j)*fact
          fz(k,j) = fv(i,j,k,4)
  621     continue
# endif
        else if(ntorder .eq. 2) then
          do 622 k=kmn+1,kmx
          fact  = cfltrb*dtl(i,j,k)*dtphy/(dtphy+cfltrb*dtl(i,j,k)*1.5)
          fact  = min(fact,100.)
          by(k,j) = by(k,j)*fact
          cy(k,j) = cy(k,j)*fact + 1.
          dy(k,j) = dy(k,j)*fact
          fy(k,j) = fv(i,j,k,3)
c
          bz(k,j) = bz(k,j)*fact
          cz(k,j) = cz(k,j)*fact + 1.
          dz(k,j) = dz(k,j)*fact
          fz(k,j) = fv(i,j,k,4)
  622     continue
        end if
c
  625 continue
c
c     solve the scalar tridiagonal equations in k-direction
c
      call tridag(mx1dwk,mx1dwk,kmn+1,kmx,jmn+1,jmx,wk2dy,by,cy,dy,fy)
      call tridag(mx1dwk,mx1dwk,kmn+1,kmx,jmn+1,jmx,wk2dz,bz,cz,dz,fz)
c
c     update the r.h.s. for A-F scheme
c
      do 630 j=jmn+1,jmx
      do 630 k=kmn+1,kmx
      fv(i,j,k,3)   = fy(k,j)
      fv(i,j,k,4)   = fz(k,j)
  630 continue
c
  610 continue
c
  611 continue
c
c
c     smooth corrections through residual smoothing
c
ccc   call psmoot (imn,  jmn,  kmn,  imx,  jmx,kmx,
ccc  .             imp1,jmp1,kmp1,imp2,jmp2,kmp2,
ccc  .             wk2dx,wk2dy,rhs,0.3,0.3,0.3)
c
c     do k=1,kmp2
c       do j=1,jmp2
c         do i=1,imp2
c           damp2(i,j,k) = max(fv(i,j,k,3),0.)
c         enddo
c       enddo
c     enddo
c
c     call psmoot (imn,  jmn,  kmn,  imx,  jmx,kmx,
c    .             imp1,jmp1,kmp1,imp2,jmp2,kmp2,
c    .             wk2dx,wk2dy,damp2,0.1,0.1,0.1)
c
c     do k=1,kmp2
c       do j=1,jmp2
c         do i=1,imp2
c           fv(i,j,k,3) = damp2(i,j,k)
c         enddo
c       enddo
c     enddo
c
c
c     do k=1,kmp2
c       do j=1,jmp2
c         do i=1,imp2
c           damp2(i,j,k) = max(fv(i,j,k,4),0.)
c         enddo
c       enddo
c     enddo
c
c     call psmoot (imn,  jmn,  kmn,  imx,  jmx,kmx,
c    .             imp1,jmp1,kmp1,imp2,jmp2,kmp2,
c    .             wk2dx,wk2dy,damp2,0.1,0.1,0.1)
c
c     do k=1,kmp2
c       do j=1,jmp2
c         do i=1,imp2
c           fv(i,j,k,4) = damp2(i,j,k)
c         enddo
c       enddo
c     enddo
c
c     update turv1  and turv2 (k and omega)
c
      sumnk     = 0.
      sumno     = 0.
      negnk     = 0
      negno     = 0
      do 700 k=kmn+1,kmx
      do 700 j=jmn+1,jmx
      do 700 i=imn+1,imx
      temp1        = fv(i,j,k,3) + turv1(i,j,k)
      turv1(i,j,k) = max(temp1,1.e-14)
      fv(i,j,k,3)  = max(fv(i,j,k,3),0.)
      incnegk      = (1 - ifix (sign(1.,temp1)) )/2
      negnk        = negnk + incnegk
      sumnk        = sumnk + (fv(i,j,k,3))**2
c
      temp2        = fv(i,j,k,4) + turv2(i,j,k)
      turv2(i,j,k) = max(temp2,1.e-14)
      fv(i,j,k,4)  = max(fv(i,j,k,4),0.)
      incnego      = (1 - ifix (sign(1.,temp2)) )/2
      negno        = negno + incnego
      sumno        = sumno + (fv(i,j,k,4))**2
 700  continue
c
c     call psmoot (imn,  jmn,  kmn,  imx,  jmx,kmx,
c    .             imp1,jmp1,kmp1,imp2,jmp2,kmp2,
c    .             wk2dx,wk2dy,turv1,0.1,0.1,0.1)
c
c     call psmoot (imn,  jmn,  kmn,  imx,  jmx,kmx,
c    .             imp1,jmp1,kmp1,imp2,jmp2,kmp2,
c    .             wk2dx,wk2dy,turv2,0.1,0.1,0.1)
c
c     write the residual history for these equations on unit 91
c
cBWW  idebug       = 1
      idebug       = 0
      idebug       = 1
      if(idebug.eq.1) then
        sumnk      = sqrt(sumnk)/real((kmx-kmn)*(jmx-jmn)*(imx-imn))
        sumno      = sqrt(sumno)/real((kmx-kmn)*(jmx-jmn)*(imx-imn))
        write (92,' (" iteraion,ibloc, log res1,res2, negnk,negno",
     .        /,2i5,2e15.5,4i7)') iter ,ibloc,sumnk,sumno,negnk,negno
c    .                           ,imx,jmx
c    .            /,i5,e15.5,i5)') iter,ibloc, alog10(sumn), negn
      endif
c
 1000 continue
c
c     update  eddy-viscosity
c  
      do 800 k=kmn+1,kmx
      do 800 j=jmn+1,jmx
      tran    = 1.0
c     if(ibloc.eq.4. and.(i-imn).le.4) tran = 0.
c     if(ibloc.eq.1) tran = 0.
      if(iwilcox.eq.1) then
        do i=imn+1,imx
          fv(i,j,k,2)  = w(i,j,k,1)*turv1(i,j,k)/turv2(i,j,k)
          fv(i,j,k,2)  = min (fv(i,j,k,2),100000.)
          enddo
      else
c
        do 810 i=imn+1,imx
c
      denom1  = 2.*sqrt(turv1(i,j,k))/(.09*turv2(i,j,k)*smin(i,j,k))
     .            *sfac
      denom2  = 500.*amu(i,j,k)*sfac*sfac/(w(i,j,k,1)*turv2(i,j,k)*
     .                                     smin(i,j,k)*smin(i,j,k))
      arg2    = max(denom1,denom2)
      f2      = tanh(arg2*arg2)
      denom   = max(a1*turv2(i,j,k),sfac*fv(i,j,k,1)*f2)
      fv(i,j,k,2)  = tran*a1*w(i,j,k,1)*turv1(i,j,k)/denom
      fv(i,j,k,2)  = min(fv(i,j,k,2),100000.)
c
c     eomu(i,j,k)  = fv(i,j,k,2)/amu(i,j,k)
 810    continue
      endif
c
        relax = 0.5
        relax = 1.0
        do i=imn+1,imx
          eomu(i,j,k) = relax*fv(i,j,k,2)/amu(i,j,k) +
     .                 (1.-relax)*eomu(i,j,k)
        enddo
 800  continue
c     call psmoot (imn,  jmn,  kmn,  imx,  jmx,kmx,
c    .             imp1,jmp1,kmp1,imp2,jmp2,kmp2,
c    .             wk2dx,wk2dy,eomu,0.1,0.1,0.1)
c
c     init = init + 1
      return
      end
