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
c=======================================================================
      subroutine psmcord(igrid,nstage,imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,
     .                   kmp1,imp2,jmp2,kmp2,ivisb,w,p,vol,x,
     .                   six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                   eomu,dti,dtj,dtk,dtl,dw)
c
c
c******  implicit smoothing by approximate jacobian *********************
c
      
c
      common/rk/    cfl,c(6),qfil(6),beta(6),vt,hm,mstage
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,prl,prt,rey,suthc,tref,i2dfl,iturb
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubcon2,iprec
      common/unsprec/ npcons,nexp,iprecg,ijac,icord,inav
c
      dimension   w(imp2,jmp2,kmp2,5), p(imp2,jmp2,kmp2),
     .          vol(imp2,jmp2,kmp2),eomu(imp2,jmp2,kmp2)
      dimension dtl(imp2,jmp2,kmp2), dti(imp2,jmp2,kmp2),
     .          dtj(imp2,jmp2,kmp2), dtk(imp2,jmp2,kmp2)
      dimension dw (imp2,jmp2,kmp2,5)

      dimension ivisb(3)
      dimension                  x(imp1,jmp1,kmp1,3),
     .     six(imp1,jmp1,kmp1),siy(imp1,jmp1,kmp1),siz(imp1,jmp1,kmp1),
     .     sjx(imp1,jmp1,kmp1),sjy(imp1,jmp1,kmp1),sjz(imp1,jmp1,kmp1),
     .     skx(imp1,jmp1,kmp1),sky(imp1,jmp1,kmp1),skz(imp1,jmp1,kmp1)
c
c           temporary storage
      dimension amu(imp2,jmp2,kmp2)
      dimension svi(imp2,jmp2,kmp2,4),svj(imp2,jmp2,kmp2,4),
     .          svk(imp2,jmp2,kmp2,4)
      dimension r00(imp2,jmp2,kmp2)  ,p00(imp2,jmp2,kmp2),
     .          u00(imp2,jmp2,kmp2)  ,v00(imp2,jmp2,kmp2),
     .          w00(imp2,jmp2,kmp2)
      dimension rs (imp2,jmp2,kmp2,5), sc(imp2,jmp2,kmp2,3)
      dimension dwn(imp2,jmp2,kmp2,5)
c
      dimension xx(5,5),rhs(5)
c
      iiprec  = iprec
      if (abs(icord).eq.2) iiprec  = 1
cvn   if (abs(icord).eq.2) iiprec  = 0
      if (icord.gt.0) then
         mrk = nstage
      else
         mrk = 1
      endif
c     mrk   = nstage
c     mrk   = 1
      ivis  = 0
      if (ivisb(1).gt.0. .or. ivisb(2).gt.0 .or. ivisb(3).gt.0) ivis = 1 
c
      ai2m     = real(1-i2dfl)
      difac    = 1.
c     difac    = 2.
      sgrmrei  = sqrt(gamma)*rm/rey
c
      imode    = ivis                 ! Navier Stokes
c
c  next two values are Cord's values
      epsx     = 0.4
c     if (igrid.eq.1) epsx = 0.5
cvn   if (igrid.eq.1) epsx = 0.45
      epsy     = 1.0
      epsim    = epsx
      epse0    = 0.15
c                c-mesh
      imesh    = 1
cSTD  niter    = 06
      niter    = 6
      tol      = 1.e-14
      do k=kmn,kmp1
      do j=jmn,jmp1
      do i=imn,imp1
c
c  creating new residuals
c    inserted 7/15/05
c
       dwn(i,j,k,1)= dw(i,j,k,1)
       dwn(i,j,k,2)= dw(i,j,k,2)
       dwn(i,j,k,3)= dw(i,j,k,3)
       dwn(i,j,k,4)= dw(i,j,k,4)
       dwn(i,j,k,5)= dw(i,j,k,5)
      enddo
      enddo
      enddo
c
c------  first step: compute approximate jacobian entries----------------
      if (mrk .eq. 1) then
c                                 first stage only
c
      do k=kmn,kmp1
      do j=jmn,jmp1
      do i=imn,imp1
        r00(i,j,k) = w(i,j,k,1)
        u00(i,j,k) = w(i,j,k,2)/w(i,j,k,1)
        v00(i,j,k) = w(i,j,k,3)/w(i,j,k,1)
        w00(i,j,k) = w(i,j,k,4)/w(i,j,k,1)
        p00(i,j,k) = p(i,j,k)
c
c  creating new residuals
c    inserted 7/15/05
c
c      dwn(i,j,k,1)= dw(i,j,k,1)
c      dwn(i,j,k,2)= dw(i,j,k,2)
c      dwn(i,j,k,3)= dw(i,j,k,3)
c      dwn(i,j,k,4)= dw(i,j,k,4)
c      dwn(i,j,k,5)= dw(i,j,k,5)
      enddo
      enddo
      enddo
c
      if (suthc.gt.0.) then
        visc1     = suthc/tref
        do k=kmn,kmp1
        do j=jmn,jmp1
        do i=imn,imp1
          ts         =  p(i,j,k)/w(i,j,k,1)
          amu(i,j,k) = ((1.+visc1)/(ts+visc1))*ts*sqrt(ts)
        enddo
        enddo
        enddo
      else
        do k=kmn,kmp1
        do j=jmn,jmp1
        do i=imn,imp1
          ts         =  p(i,j,k)/w(i,j,k,1)
          amu(i,j,k) =  ts
        enddo
        enddo
        enddo
      endif
      if (iturb.gt.0) then
        do k=kmn,kmp1
        do j=jmn,jmp1
        do i=imn,imp1
          amu(i,j,k) = amu(i,j,k)*(1.+eomu(i,j,k))
        enddo
        enddo
        enddo
      endif
c

c
      if (imode .ge. 1) then
c
c-----------navier-stokes
        if (iiprec .eq. 0) then
c  next 2 values are standard values used by Cord
cSTD     hfl1 = 0.3
cSTD     hfl2 = 0.1
c        hfl1 = 0.2
         hfl2 = 0.05
         hfl1 = 0.15
ccs
cvn      hfl1 = 0.2
cvn      hfl2 = 0.1
        else
c  Note: rm needs to be defined !!!!
cSTD     hfl1 = 0.30*rm
cSTD     hfl2 = 0.10*rm
c        hfl1 = 0.15*rm
c        hfl2 = 0.05*rm
         hfl1 = 0.20*rm
         hfl2 = 0.10*rm
c
c        hfl1 = max(hfl1,.01)
c        hfl2 = max(hfl2,.0025)
c         hfl1 = 0.07
c         hfl2 = 0.01
        end if
c
      else if (imode.eq.0) then
c-----------euler
        if (iiprec.eq.0) then
         hfl1 = 0.20
c        hfl1 = 0.15
         hfl2 = 0.02
        else
         hfl1 = 0.10*rm
         hfl2 = 0.01*rm
        end if
c  finish of imode end if (NS)
      end if
c
c      **********************************************
c      ******* coefficients in i-direction   ********
c      **********************************************
      do 100 k=kmn+1,kmx
      do 100 j=jmn+1,jmx
      do 100 i=imn+1,imp1  
c
        sxa    = six(i-1,j,k)
        sya    = siy(i-1,j,k)
        sza    = siz(i-1,j,k)
        sar    = sqrt(sxa**2+sya**2+sza**2)
        svol   = 0.5*(vol(i-1,j,k)+vol(i,j,k))
c
c       left-right variables
        rl     = r00(i-1,j,k)
        pl     = p00(i-1,j,k)
        ul     = u00(i-1,j,k)
        vl     = v00(i-1,j,k)
        wl     = w00(i-1,j,k)
c
        rr     = r00(i,j,k)
        pr     = p00(i,j,k)
        ur     = u00(i,j,k)
        vr     = v00(i,j,k)
        wr     = w00(i,j,k)
c
c----  arithmetic averaging
        rhoa   = 0.5*(rl+rr)
        us     = 0.5*(ul+ur)
        vs     = 0.5*(vl+vr)
        ws     = 0.5*(wl+wr)
        ps     = 0.5*(pl+pr)
        qq     = us*us + vs*vs + ws*ws
        qs     = (sxa*us + sya*vs + sza*ws)/sar
        sigm   = sign(1.0,qs)
        qa     = abs(qs)
        ccs    = ps/rhoa*gamma
        cs     = sqrt(ccs)
c
c--------------Mach  numbers
        fms    = qa/cs
        fm0    = min(fms,1.0)
c
        epse   = max( 0.0, min(epse0, 1.+epse0-fms) )
        fme    = min( 1.-epse, fm0 )
c
c        f1m0  = 1.-fm0
        f1m0   = 1.-fme
        fm0    = fm0*sigm
c
c--------------prevent  zero dissipation at qs=0
c
        dti1   =      sqrt(dtj(i  ,j,k)/dti(i  ,j,k))
     .         + ai2m*sqrt(dtk(i  ,j,k)/dti(i  ,j,k))
        dtim   =      sqrt(dtj(i-1,j,k)/dti(i-1,j,k))
     .         + ai2m*sqrt(dtk(i-1,j,k)/dti(i-1,j,k))
        fasp   = max(0.5*(dti1+dtim)    ,1.e-30)
        dfak   = max(fasp,1.)
        dcut   = min(fasp,1.)
c  hfl1 = 0.3, hfl2 = 0.1 (Cord's values)
        dcut   = max(hfl1*dcut,hfl2)
        fma    = qa/cs
        fmd    = max(fma,dcut)
        sqa    = fmd*cs*dfak * epsy
c--------------------------------------------------------
        cprec  = cs
c              for preconditioning only
        if (iiprec.ge.1) then
c                      ubmin2    = (ubmin*rm)**2
          qref2 = min( max( upk*qq, ubmin2*ccs), ccs)
c         alp   = (1.-qref2/ccs)/2.0
          alp   = 0.
          cprec = sqrt(qs*qs * alp*alp + qref2)
        endif
c
c---------viscous  jacobian (primitive)
c
        scale  = sar*sar/(svol*rhoa) * epsy
c  next line inserted 7/17/05 to account for nondimensionalization
        scale  = difac*scale*sgrmrei
        if (ivisb(1) .eq. 0) scale = 0.0
        rmue   = 0.5*(amu (i-1,j,k) + amu (i,j,k))
        rmuem  = rmue*scale
c
c--------------store  relevant face quantities
c
        svi(i,j,k,1) = f1m0
        svi(i,j,k,2) = fm0
        svi(i,j,k,3) = sqa
        svi(i,j,k,4) = rmuem
        sc (i,j,k,1) = cprec
c
  100 continue
c
c      **********************************************
c      ******* coefficients in j-direction   ********
c      **********************************************
      do 200 k=kmn+1,kmx
      do 200 j=jmn+1,jmp1
      do 200 i=imn+1,imx
c
        sxa    = sjx(i,j-1,k)
        sya    = sjy(i,j-1,k)
        sza    = sjz(i,j-1,k)
        sar    = sqrt(sxa**2+sya**2+sza**2)
        svol   = 0.5*(vol(i,j-1,k)+vol(i,j,k))
c
c
c       left-right variables
        rl     = r00(i,j-1,k)
        pl     = p00(i,j-1,k)
        ul     = u00(i,j-1,k)
        vl     = v00(i,j-1,k)
        wl     = w00(i,j-1,k)
c
        rr     = r00(i,j,k)
        pr     = p00(i,j,k)
        ur     = u00(i,j,k)
        vr     = v00(i,j,k)
        wr     = w00(i,j,k)
c
c----  arithmetic averaging
        rhoa   = 0.5*(rl+rr)
        us     = 0.5*(ul+ur)
        vs     = 0.5*(vl+vr)
        ws     = 0.5*(wl+wr)
        ps     = 0.5*(pl+pr)
        qq     = us*us + vs*vs + ws*ws
        qs     = (sxa*us + sya*vs + sza*ws)/sar
        sigm   = sign(1.0,qs)
        qa     = abs(qs)
        ccs    = ps/rhoa*gamma
        cs     = sqrt(ccs)
c
c--------------Mach  numbers
        fms    = qa/cs
        fm0    = min(fms,1.0)
c
        epse   = max( 0.0, min(epse0, 1.+epse0-fms) )
        fme    = min( 1.-epse, fm0 )
c
c        f1m0  = 1.-fm0
        f1m0   = 1.-fme
        fm0    = fm0*sigm
c
c--------------prevent  zero dissipation at qs=0
c
        dtj1   =      sqrt(dti(i,j  ,k)/dtj(i,j  ,k))
     .         + ai2m*sqrt(dtk(i,j  ,k)/dtj(i,j  ,k))
        dtjm   =      sqrt(dti(i,j-1,k)/dtj(i,j-1,k))
     .         + ai2m*sqrt(dtk(i,j-1,k)/dtj(i,j-1,k))
        fasp   = max(0.5*(dtj1+dtjm)    ,1.e-30)
        dfak   = max(fasp,1.)
c******************************************************   ??
        dcut   = min(fasp,1.)
        dcut   = max(hfl1*dcut,hfl2)
        fma    = qa/cs
        fmd    = max(fma,dcut)
        sqa    = fmd*cs*dfak * epsy
c--------------------------------------------------------
        cprec  = cs
c               for preconditioning only
        if (iiprec.ge.1) then
          qref2 = min( max( upk*qq, ubmin2*ccs), ccs)
c         alp   = (1.-qref2/ccs)/2.0
          alp   = 0.
          cprec = sqrt(qs*qs * alp*alp + qref2)
        endif
c
c--------------------------------------------------------
c
c---------viscous  jacobian (primitive)
c
        scale  = sar*sar/(svol*rhoa) * epsy
c  next line inserted 7/17/05 to account for nondimensionalization
        scale  = difac*scale*sgrmrei
        if (ivisb(2) .eq. 0) scale = 0.0
        rmue   = 0.5*(amu (i,j-1,k) + amu (i,j,k))
        rmuem  = rmue*scale
c
c--------------store  relevant face quantities
c
        svj(i,j,k,1) = f1m0
        svj(i,j,k,2) = fm0
        svj(i,j,k,3) = sqa
        svj(i,j,k,4) = rmuem
        sc (i,j,k,2) = cprec
c
  200 continue
c
c      **********************************************
c      ******* coefficients in k-direction   ********
c      **********************************************
      if (i2dfl.eq.1) go to 301
      do 300 k=kmn+1,kmp1
      do 300 j=jmn+1,jmx
      do 300 i=imn+1,imx
c
        sxa    = skx(i,j,k-1)
        sya    = sky(i,j,k-1)
        ska    = skz(i,j,k-1)
        sar    = sqrt(sxa**2+sya**2+sza**2)
        svol   = 0.5*(vol(i,j,k-1)+vol(i,j,k))
c
c
c       left-right variables
        rl     = r00(i,j,k-1)
        pl     = p00(i,j,k-1)
        ul     = u00(i,j,k-1)
        vl     = v00(i,j,k-1)
        wl     = w00(i,j,k-1)
c
        rr     = r00(i,j,k)
        pr     = p00(i,j,k)
        ur     = u00(i,j,k)
        vr     = v00(i,j,k)
        wr     = w00(i,j,k)
c
c---- arithmetic averaging
        rhoa   = 0.5*(rl+rr)
        us     = 0.5*(ul+ur)
        vs     = 0.5*(vl+vr)
        ws     = 0.5*(wl+wr)
        qq     = us*us + vs*vs + ws*ws
        ps     = 0.5*(pl+pr)
        qs     = (sxa*us + sya*vs + sza*ws)/sar
        sigm   = sign(1.0,qs)
        qa     = abs(qs)
        ccs    = ps/rhoa*gamma
        cs     = sqrt(ccs)
c
c--------------Mach  numbers
        fms    = qa/cs
        fm0    = min(fms,1.0)
c
        epse   = max( 0.0, min(epse0, 1.+epse0-fms) )
        fme    = min( 1.-epse, fm0 )
c
c        f1m0  = 1.-fm0
        f1m0   = 1.-fme
        fm0    = fm0*sigm
c
c--------------prevent  zero dissipation at qs=0
c
        dtk1   = 1. + sqrt(dti(i,j,k  )/dtk(i,j,k  ))
     .         +      sqrt(dtj(i,j,k  )/dtk(i,j,k  ))
        dtkm   = 1. + sqrt(dti(i,j,k-1)/dtk(i,j,k-1))
     .         +      sqrt(dtj(i,j,k-1)/dtk(i,j,k-1))
        fasp   = max(0.5*(dtk1+dtkm)-1.0,1.e-30)
        dfak   = max(fasp,1.)
c******************************************************   ??
        dcut   = min(fasp,1.)
        dcut   = max(hfl1*dcut,hfl2)
        fma    = qa/cs
        fmd    = max(fma,dcut)
        sqa    = fmd*cs*dfak * epsy
c
c-------------------------------------------------------
        cprec  = cs
c             for preconditioning only
        if (iiprec.ge.1) then
          qref2 = min( max( upk*qq, ubmin2*ccs), ccs)
          alp   = (1.-qref2/ccs)/2.0
          alp   = 0.
          cprec = sqrt(qs*qs * alp*alp + qref2)
        endif
c
c--------------------------------------------------------
c
c---------viscous  jacobian (primitive)
c
        scale  = sar*sar/(svol*rhoa) * epsy
c  next line inserted 7/17/05 to account for nondimensionalization
        scale  = difac*scale* sgrmrei
        if (ivisb(3) .eq. 0) scale = 0.0
        rmue   = 0.5*(amu (i,j,k-1) + amu (i,j,k))
        rmuem  = rmue*scale
c
c--------------store  relevant face quantities
c
        svk(i,j,k,1) = f1m0
        svk(i,j,k,2) = fm0
        svk(i,j,k,3) = sqa
        svk(i,j,k,4) = rmuem
        sc (i,j,k,3) = cprec
c
  300 continue
  301 continue
c          finish mrk=1
      end if 
c            end of "frozen" material at first stage
c
c
c------  second step: set up and solve linear implicit system----
c
      if (abs(iprec).le.1 ) then
        call res_cons_to_prim(imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,w,dwn)
      else
        call res_puvt_to_prim(imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,w,p,dwn)
      endif
c
c      **********************************************
c      ******* preset smoothed residuals ************
c      **********************************************
      call bcresi(imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,dwn)
c
      if (iiprec.eq.0) then
cvn   if (iiprec.lt.0) then
        do k=kmn,kmp1 
        do j=jmn,jmp1 
        do i=imn,imp1 
          rs(i,j,k,1) = dwn(i,j,k,1)
          rs(i,j,k,2) = dwn(i,j,k,2)
          rs(i,j,k,3) = dwn(i,j,k,3)
          rs(i,j,k,4) = dwn(i,j,k,4)
          rs(i,j,k,5) = dwn(i,j,k,5)
        end do
        end do
        end do
      else
        do k=kmn,kmp1 
        do j=jmn,jmp1 
        do i=imn,imp1 
          rs(i,j,k,1) = 0.
          rs(i,j,k,2) = 0.
          rs(i,j,k,3) = 0.
          rs(i,j,k,4) = 0.
          rs(i,j,k,5) = 0.
        end do
        end do
        end do
      endif

c
c      **********************************************
c      ******* solution by gauss-seidel iteration**********
c      **********************************************
c      goto 3000
c
      idir =  -1
c     icount = 0
c     ilimit = 0
c
c      write(6,*) ' '
c      write(6,*) ' starting iteration '
c
      do 1000 iter = 1,niter
c
c     icount = icount+1
c
c     if(icount.gt.ilimit) then
c       idir = -idir
c       icount = 0
c     endif

      idir = -idir

c     write(6,*) ' iter = ', iter, ' idir = ', idir
c
      call bcresj(imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,rs)
c
      chr     = 0.0
      chu     = 0.0
      chv     = 0.0
      chw     = 0.0
      che     = 0.0
      chr0    = 1.0
      chu0    = 1.0
      chv0    = 1.0
      chw0    = 1.0
      che0    = 1.0
c
c
      if (idir .eq. 1) then
        ia = imn+1
        ja = jmn+1
        ka = kmn+1
        ib = imx
        jb = jmx
        kb = kmx
      else
        ia = imx
        ja = jmx
        ka = kmx
        ib = imn+1
        jb = jmn+1
        kb = kmn+1
      end if
 
      do 2000 k=ka,kb,idir
      do 2000 j=ja,jb,idir
      do 2000 i=ia,ib,idir
c
c       if (idir .eq. 1) then
c         j      = no
c         i      = ni
c       elseif (idir .eq. -1) then
c         j      = nj_e - no+2
c         i      = ni_e - ni+2
c       else
c         write(6,*) 'error in gauss-seidel with directions'
c         stop
c       endif
c
        do m=1,5
        do n=1,5
          xx(m,n) = 0.
        end do
        end do
c
        dtv2   = 0.5*cfl*dtl(i,j,k)/vol(i,j,k) * epsim
c
c.........left  matrix
        a1  = rs(i-1,j,k,1)
        a2  = rs(i-1,j,k,2)
        a3  = rs(i-1,j,k,3)
        a4  = rs(i-1,j,k,4)
        a5  = rs(i-1,j,k,5)
c
        sigf   = -1.0
c
        f1m0   = svi(i,j,k,1)
        fm0    = svi(i,j,k,2) * sigf
        qsa    = svi(i,j,k,3)
        rmuem  = svi(i,j,k,4)
        cp     = sc (i,j,k,1)
        sxa    = six(i-1,j,k)
        sya    = siy(i-1,j,k)
        sza    = siz(i-1,j,k)
c
        sar    = sqrt(sxa**2+sya**2+sza**2)
        vnx    = sxa/sar * sigf
        vny    = sya/sar * sigf
        vnz    = sza/sar * sigf
c
c       left-right variables
        rl     = r00(i-1,j,k)
        pl     = p00(i-1,j,k)
        ul     = u00(i-1,j,k)
        vl     = v00(i-1,j,k)
        wl     = w00(i-1,j,k)
c
        rr     = r00(i,j,k)
        pr     = p00(i,j,k)
        ur     = u00(i,j,k)
        vr     = v00(i,j,k)
        wr     = w00(i,j,k)
c
        call rhs_face(vnx, vny, vnz, sar, dtv2,
     .                f1m0, fm0, qsa, rmuem, cp,
     .                rl, pl, ul, vl, wl,
     .                rr, pr, ur, vr, wr,
     .                a1, a2, a3, a4, a5,
     .                rl1, rl2, rl3, rl4, rl5,
     .                xx)

c
c.........right  matrix
        a1  = rs(i+1,j,k,1)
        a2  = rs(i+1,j,k,2)
        a3  = rs(i+1,j,k,3)
        a4  = rs(i+1,j,k,4)
        a5  = rs(i+1,j,k,5)
c
        sigf   =  1.0
c
        f1m0   = svi(i+1,j,k,1)
        fm0    = svi(i+1,j,k,2) * sigf
        qsa    = svi(i+1,j,k,3)
        rmuem  = svi(i+1,j,k,4)
        cp     = sc (i+1,j,k,1)
        sxa    = six(i,j,k)
        sya    = siy(i,j,k)
        sza    = siz(i,j,k)
c
        sar    = sqrt(sxa**2+sya**2+sza**2)
        vnx    = sxa/sar * sigf
        vny    = sya/sar * sigf
        vnz    = sza/sar * sigf
c
c       left-right variables
        rl     = r00(i,j,k)
        pl     = p00(i,j,k)
        ul     = u00(i,j,k)
        vl     = v00(i,j,k)
        wl     = w00(i,j,k)
c
        rr     = r00(i+1,j,k)
        pr     = p00(i+1,j,k)
        ur     = u00(i+1,j,k)
        vr     = v00(i+1,j,k)
        wr     = w00(i+1,j,k)
c
c
        call rhs_face(vnx, vny, vnz, sar, dtv2,
     .                f1m0, fm0, qsa, rmuem, cp,
     .                rl, pl, ul, vl, wl,
     .                rr, pr, ur, vr, wr,
     .                a1, a2, a3, a4, a5,
     .                rr1, rr2, rr3, rr4, rr5,
     .                xx)
c
c.........y bottom  matrix
c
        a1  = rs(i,j-1,k,1)
        a2  = rs(i,j-1,k,2)
        a3  = rs(i,j-1,k,3)
        a4  = rs(i,j-1,k,4)
        a5  = rs(i,j-1,k,5)
c
        sigf   = -1.0
c
        f1m0   = svj(i,j,k,1)
        fm0    = svj(i,j,k,2) * sigf
        qsa    = svj(i,j,k,3)
        rmuem  = svj(i,j,k,4)
        cp     = sc (i,j,k,2)
        sxa    = sjx(i,j-1,k)
        sya    = sjy(i,j-1,k)
        sza    = sjz(i,j-1,k)
c
        sar    = sqrt(sxa**2+sya**2+sza**2)
        vnx    = sxa/sar * sigf
        vny    = sya/sar * sigf
        vnz    = sza/sar * sigf
c
c       left-right variables
          rl      = r00(i,j-1,k)
          pl      = p00(i,j-1,k)
          ul      = u00(i,j-1,k)
          vl      = v00(i,j-1,k)
          wl      = w00(i,j-1,k)
c
          rr      = r00(i,j,k)
          pr      = p00(i,j,k)
          ur      = u00(i,j,k)
          vr      = v00(i,j,k)
          wr      = w00(i,j,k)
c
        call rhs_face(vnx, vny, vnz, sar, dtv2,
     .                f1m0, fm0, qsa, rmuem, cp,
     .                rl, pl, ul, vl, wl,
     .                rr, pr, ur, vr, wr,
     .                a1, a2, a3, a4, a5,
     .                rb1, rb2, rb3, rb4, rb5,
     .                xx)
c
c.........y top  matrix
        a1  = rs(i,j+1,k,1)
        a2  = rs(i,j+1,k,2)
        a3  = rs(i,j+1,k,3)
        a4  = rs(i,j+1,k,4)
        a5  = rs(i,j+1,k,5)
c
        sigf   =  1.0
c
        f1m0   = svj(i,j+1,k,1)
        fm0    = svj(i,j+1,k,2) * sigf
        qsa    = svj(i,j+1,k,3)
        rmuem  = svj(i,j+1,k,4)
        cp     = sc (i,j+1,k,2)
        sxa    = sjx(i,j,k)
        sya    = sjy(i,j,k)
        sza    = sjz(i,j,k)
c
        sar    = sqrt(sxa**2+sya**2+sza**2)
        vnx    = sxa/sar * sigf
        vny    = sya/sar * sigf
        vnz    = sza/sar * sigf
c
c       left-right variables
          rl      = r00(i,j,k)
          pl      = p00(i,j,k)
          ul      = u00(i,j,k)
          vl      = v00(i,j,k)
          wl      = w00(i,j,k)
c
          rr      = r00(i,j+1,k)
          pr      = p00(i,j+1,k)
          ur      = u00(i,j+1,k)
          vr      = v00(i,j+1,k)
          wr      = w00(i,j+1,k)
c
c
        call rhs_face(vnx, vny, vnz, sar, dtv2,
     .                f1m0, fm0, qsa, rmuem, cp,
     .                rl, pl, ul, vl, wl,
     .                rr, pr, ur, vr, wr,
     .                a1, a2, a3, a4, a5,
     .                rt1, rt2, rt3, rt4, rt5,
     .                xx)

      if (i2dfl.ne.1) then
c       three dimensional case
c
c.........k-bottom  matrix
        a1  = rs(i,j,k-1,1)
        a2  = rs(i,j,k-1,2)
        a3  = rs(i,j,k-1,3)
        a4  = rs(i,j,k-1,4)
        a5  = rs(i,j,k-1,5)
c
        sigf   = -1.0
c
        f1m0   = svk(i,j,k,1)
        fm0    = svk(i,j,k,2) * sigf
        qsa    = svk(i,j,k,3)
        rmuem  = svk(i,j,k,4)
        cp     = sc (i,j,k,3)
        sxa    = skx(i,j,k-1)
        sya    = sky(i,j,k-1)
        sza    = skz(i,j,k-1)
c
        sar    = sqrt(sxa**2+sya**2+sza**2)
        vnx    = sxa/sar * sigf
        vny    = sya/sar * sigf
        vnz    = sza/sar * sigf
c
c       left-right variables
          rl      = r00(i,j,k-1)
          pl      = p00(i,j,k-1)
          ul      = u00(i,j,k-1)
          vl      = v00(i,j,k-1)
          wl      = w00(i,j,k-1)
c
          rr      = r00(i,j,k)
          pr      = p00(i,j,k)
          ur      = u00(i,j,k)
          vr      = v00(i,j,k)
          wr      = w00(i,j,k)
c
        call rhs_face(vnx, vny, vnz, sar, dtv2,
     .                f1m0, fm0, qsa, rmuem, cp,
     .                rl, pl, ul, vl, wl,
     .                rr, pr, ur, vr, wr,
     .                a1, a2, a3, a4, a5,
     .                ry1, ry2, ry3, ry4, ry5,
     .                xx)
c
c.........k-top  matrix
        a1  = rs(i,j,k+1,1)
        a2  = rs(i,j,k+1,2)
        a3  = rs(i,j,k+1,3)
        a4  = rs(i,j,k+1,4)
        a5  = rs(i,j,k+1,5)
c
        sigf   =  1.0
c
        f1m0   = svk(i,j,k+1,1)
        fm0    = svk(i,j,k+1,2) * sigf
        qsa    = svk(i,j,k+1,3)
        rmuem  = svk(i,j,k+1,4)
        cp     = sc (i,j,k+1,3)
        sxa    = skx(i,j,k)
        sya    = sky(i,j,k)
        sza    = skz(i,j,k)
c
        sar    = sqrt(sxa**2+sya**2+sza**2)
        vnx    = sxa/sar * sigf
        vny    = sya/sar * sigf
        vnz    = sza/sar * sigf
c
c       left-right variables
          rl      = r00(i,j,k)
          pl      = p00(i,j,k)
          ul      = u00(i,j,k)
          vl      = v00(i,j,k)
          wl      = w00(i,j,k)
c
          rr      = r00(i,j,k+1)
          pr      = p00(i,j,k+1)
          ur      = u00(i,j,k+1)
          vr      = v00(i,j,k+1)
          wr      = w00(i,j,k+1)
c
        call rhs_face(vnx, vny, vnz, sar, dtv2,
     .                f1m0, fm0, qsa, rmuem, cp,
     .                rl, pl, ul, vl, wl,
     .                rr, pr, ur, vr, wr,
     .                a1, a2, a3, a4, a5,
     .                rz1, rz2, rz3, rz4, rz5,
     .                xx)
      endif
c
c************assemble  complete right-hand-side
      if (i2dfl.ne.1) then
        rhs(1) =      dwn(i,j,k,1) - rl1 - rr1 - rb1 - rt1 - ry1 - rz1
        rhs(2) =      dwn(i,j,k,2) - rl2 - rr2 - rb2 - rt2 - ry2 - rz2
        rhs(3) =      dwn(i,j,k,3) - rl3 - rr3 - rb3 - rt3 - ry3 - rz3
        rhs(4) =      dwn(i,j,k,4) - rl4 - rr4 - rb4 - rt4 - ry4 - rz4
        rhs(5) =      dwn(i,j,k,5) - rl5 - rr5 - rb5 - rt5 - ry5 - rz5
      else
c       two dimensional case
        rhs(1) =      dwn(i,j,k,1) - rl1 - rr1 - rb1 - rt1
        rhs(2) =      dwn(i,j,k,2) - rl2 - rr2 - rb2 - rt2
        rhs(3) =      dwn(i,j,k,3) - rl3 - rr3 - rb3 - rt3
        rhs(4) =      dwn(i,j,k,4) - rl4 - rr4 - rb4 - rt4
        rhs(5) =      dwn(i,j,k,5) - rl5 - rr5 - rb5 - rt5

      endif
c
c--------*----------*---------*---------*---------*---------*---------*---------*
c
        xx(1,1)   = 1. + xx(1,1)
        xx(2,2)   = 1. + xx(2,2)
        xx(3,3)   = 1. + xx(3,3)
        xx(4,4)   = 1. + xx(4,4)
        xx(5,5)   = 1. + xx(5,5)
c
        difr      = rs(i,j,k,1)
        difu      = rs(i,j,k,2)
        difv      = rs(i,j,k,3)
        difw      = rs(i,j,k,4)
        dife      = rs(i,j,k,5)
c
        rhs1      = rhs(1)
        rhs2      = rhs(2)
        rhs3      = rhs(3)
        rhs4      = rhs(4)
        rhs5      = rhs(5)

        if (i2dfl.eq.1) call solve4_rs(rhs,xx)
        if (i2dfl.ne.1) call solve5_rs(rhs,xx)
c
        rs(i,j,k,1) = rhs(1)
        rs(i,j,k,2) = rhs(2)
        rs(i,j,k,3) = rhs(3)
        rs(i,j,k,4) = rhs(4)
        rs(i,j,k,5) = rhs(5)
c
        difr   = abs(abs(difr)-abs(rs(i,j,k,1)))
        difu   = abs(abs(difu)-abs(rs(i,j,k,2)))
        difv   = abs(abs(difv)-abs(rs(i,j,k,3)))
        difw   = abs(abs(difw)-abs(rs(i,j,k,4)))
        dife   = abs(abs(dife)-abs(rs(i,j,k,5)))
c
cR      chr    = max(difr,chr,1.e-12)
cR      chu    = max(difu,chu,1.e-12)
cR      chv    = max(difv,chv,1.e-12)
cR      chw    = max(difw,chv,1.e-12)
cR      che    = max(dife,che,1.e-12)
        chr    = max(difr,chr*chr0,1.e-12)
        chu    = max(difu,chu*chu0,1.e-12)
        chv    = max(difv,chv*chv0,1.e-12)
        chw    = max(difw,chw*chw0,1.e-12)
        che    = max(dife,che*che0,1.e-12)
c
c       chr    = max(difr,chr)
c       chu    = max(difu,chu)
c       chv    = max(difv,chv)
c       chw    = max(difw,chw)
c       che    = max(dife,che)
c
 2000 continue
c
      if (iter .eq. 1) then
        chr0  = chr
        chu0  = chu
        chv0  = chv
        chw0  = chw
        che0  = che
      end if
c
      chr    = chr/chr0
      chu    = chu/chu0
      chv    = chv/chv0
      chw    = chw/chw0
      che    = che/che0
      check  = max(chr,chu,chv,chw,che)
c     check  = chr
c     write (530,931) iter,check,tol
c931  format (5x,'iter,check,tol = ',i5,2e15.6)
      if (check .le. tol) go to 3000
c
 1000 continue                ! end GS
c
 3000 continue
c
      do k=kmn+1,kmx
      do j=jmn+1,jmx
      do i=imn+1,imx
        dwn(i,j,k,1) = rs(i,j,k,1)
        dwn(i,j,k,2) = rs(i,j,k,2)
        dwn(i,j,k,3) = rs(i,j,k,3)
        dwn(i,j,k,4) = rs(i,j,k,4)
        dwn(i,j,k,5) = rs(i,j,k,5)
      end do
      end do
      end do
c
      if (abs(iprec).le.1) then
        call res_prim_to_cons(imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,w,dw,dwn)
      else
        call res_prim_to_puvt(imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,w,p,dw,dwn)
      endif

      return
      end
c
      subroutine  rhs_face(vnx, vny, vnz, sar, dtv2,
     .                     f1m0, fm0, qsa, rmuem, cp,
     .                     rl, pl, ul, vl, wl,
     .                     rr, pr, ur, vr, wr,
     .                     a1, a2, a3, a4, a5,
     .                     rhs_f1, rhs_f2, rhs_f3, rhs_f4, rhs_f5,
     .                     xx)
c
c******  computation of right-hand-side face-coefficients for sgs********
c
c     implicit double precision (a-h,o-z)
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,prl,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension     xx(5,5)
c
c
c----  arithmetic averaging
         rhoa  = 0.5*(rl+rr)
         us    = 0.5*(ul+ur)
         vs    = 0.5*(vl+vr)
         ws    = 0.5*(wl+wr)
         ps    = 0.5*(pl+pr)
         qs    = vnx*us + vny*vs + vnz*ws
         ccs   = ps/rhoa*gamma
         cs    = sqrt(ccs)
         hs    = ccs/gm1
         rc    = 1./cp
         ccp   = cp*cp
         rmuem = rmuem             *dtv2
c        grm   = prl/gogm1*gamma * rmuem
c  CHECK !!!!!!!
c        grm   = gamma/prl * rmuem
         grm   = gamma/(prl*gm1) * rmuem
c
         sqs   = sar*(qs-qsa)
         sc0   = sar*(1.-fm0)
         sc1   = sar*rc*f1m0
c
         dqs   = vnx*a3 + vny*a4 + vnz*a5
c
         sc0   = sc0               *dtv2
         sc0rr = sc0/rhoa
         sc0p  = sc0rr*a2
c
         sc0r  = sc0*rhoa
         sc0gp = sc0*gamma*ps
         sc0q  = sc0*dqs
c
         sc1   = sc1               *dtv2
         sc1gh = sc1*gm1*hs
c        sc1gh = sc1*ccs
         sc1p  = sc1*a2
         sc1ghp= sc1gh*a2
         sc1v  = sc1*ccp+rmuem/3.
         sc1vx = sc1v*vnx
         sc1vy = sc1v*vny
         sc1vz = sc1v*vnz
c
         sqs   = sqs               *dtv2
         sqv   = sqs-rmuem
         sqp   = sqs-grm

c
         rhs_f1   = sqs*a1 - sc1p     + sc0r *dqs
         rhs_f2   = sqp*a2 - sc1ghp   + sc0gp*dqs
         rhs_f3   = sqv*a3 + vnx*sc0p - sc1vx*dqs
         rhs_f4   = sqv*a4 + vny*sc0p - sc1vy*dqs
         rhs_f5   = sqv*a5 + vnz*sc0p - sc1vz*dqs
c
c
         sqs   = sar*(-qs-qsa)
         sc0   = sar*(1.+fm0)
c
         sc0   = sc0               *dtv2
         sc0rr = sc0/rhoa
c
         sc0r  = sc0*rhoa
         sc0gp = sc0*gamma*ps
c
         sqs   = sqs               *dtv2
         sqv   = sqs-rmuem
         sqp   = sqs-grm
c

      xx(1,1)   = -sqs                         + xx(1,1)
      xx(1,2)   =  sc1                         + xx(1,2)
      xx(1,3)   =  vnx*sc0r                    + xx(1,3)
      xx(1,4)   =  vny*sc0r                    + xx(1,4)
      xx(1,5)   =  vnz*sc0r                    + xx(1,5)
c
      xx(2,1)   =  0.0
      xx(2,2)   = -sqp + sc1gh                 + xx(2,2)
      xx(2,3)   =                    vnx*sc0gp + xx(2,3)
      xx(2,4)   =                    vny*sc0gp + xx(2,4)
      xx(2,5)   =                    vnz*sc0gp + xx(2,5)
c
c
      xx(3,1)   =  0.0
      xx(3,2)   =        vnx*sc0rr             + xx(3,2)
      xx(3,3)   = -sqv             + vnx*sc1vx + xx(3,3)
      xx(3,4)   =                    vny*sc1vx + xx(3,4)
      xx(3,5)   =                    vnz*sc1vx + xx(3,5)
c
      xx(4,1)   =  0.0
      xx(4,2)   =        vny*sc0rr             + xx(4,2)
      xx(4,3)   =                    vnx*sc1vy + xx(4,3)
      xx(4,4)   = -sqv             + vny*sc1vy + xx(4,4)
      xx(4,5)   =                    vnz*sc1vy + xx(4,5)

      xx(5,1)   =  0.0
      xx(5,2)   =        vnz*sc0rr             + xx(5,2)
      xx(5,3)   =                    vnx*sc1vz + xx(5,3)
      xx(5,4)   =                    vny*sc1vz + xx(5,4)
      xx(5,5)   = -sqv             + vnz*sc1vz + xx(5,5)
c
      return
      end
c
c
c
      subroutine  bcresi(imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,dwn)
c
c*****   set boundary values to zero and implement cut bc's
c*****    for cell centered variables
c
c
      dimension dwn(imp2,jmp2,kmp2,5)
c
      do 10 k= kmn,kmp1
      do 10 i= imn,imp1
        dwn(i,jmn,k,1)  = 0.0
        dwn(i,jmn,k,2)  = 0.0
        dwn(i,jmn,k,3)  = 0.0
        dwn(i,jmn,k,4)  = 0.0
        dwn(i,jmn,k,5)  = 0.0
c
        dwn(i,jmp1,k,1) = 0.0
        dwn(i,jmp1,k,2) = 0.0
        dwn(i,jmp1,k,3) = 0.0
        dwn(i,jmp1,k,4) = 0.0
        dwn(i,jmp1,k,5) = 0.0
   10 continue
c
      do 20 k= kmn,kmp1
      do 20 j= jmn,jmp1
        dwn(imn,j,k,1)  = 0.0
        dwn(imn,j,k,2)  = 0.0
        dwn(imn,j,k,3)  = 0.0
        dwn(imn,j,k,4)  = 0.0
        dwn(imn,j,k,5)  = 0.0
c
        dwn(imp1,j,k,1) = 0.0
        dwn(imp1,j,k,2) = 0.0
        dwn(imp1,j,k,3) = 0.0
        dwn(imp1,j,k,4) = 0.0
        dwn(imp1,j,k,5) = 0.0
   20 continue
c
      do 30 j= jmn,jmp1
      do 30 i= imn,imp1
        dwn(i,j,kmn,1)  = 0.0
        dwn(i,j,kmn,2)  = 0.0
        dwn(i,j,kmn,3)  = 0.0
        dwn(i,j,kmn,4)  = 0.0
        dwn(i,j,kmn,5)  = 0.0
c
        dwn(i,j,kmp1,1) = 0.0
        dwn(i,j,kmp1,2) = 0.0
        dwn(i,j,kmp1,3) = 0.0
        dwn(i,j,kmp1,4) = 0.0
        dwn(i,j,kmp1,5) = 0.0
   30 continue
c
      return
      end
c
c
c
      subroutine  bcresj(imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,rs)
c
c*****   set boundary values to zero and implement cut bc's
c*****  for face values of residuals
c
c
      dimension   rs(imp2,jmp2,kmp2,5)
c
      do 10 k=kmn,kmp1
      do 10 i=imn,imp1
        rs(i,jmn,k,1)  = 0.0
        rs(i,jmn,k,2)  = 0.0
        rs(i,jmn,k,3)  = 0.0
        rs(i,jmn,k,4)  = 0.0
        rs(i,jmn,k,5)  = 0.0
c
        rs(i,jmp1,k,1) = 0.0
        rs(i,jmp1,k,2) = 0.0
        rs(i,jmp1,k,3) = 0.0
        rs(i,jmp1,k,4) = 0.0
        rs(i,jmp1,k,5) = 0.0
   10 continue
c
      do 20 k=kmn,kmp1
      do 20 j=jmn,jmp1
        rs(imn,j,k,1)  = 0.0
        rs(imn,j,k,2)  = 0.0
        rs(imn,j,k,3)  = 0.0
        rs(imn,j,k,4)  = 0.0
        rs(imn,j,k,5)  = 0.0
c
        rs(imp1,j,k,1) = 0.0
        rs(imp1,j,k,2) = 0.0
        rs(imp1,j,k,3) = 0.0
        rs(imp1,j,k,4) = 0.0
        rs(imp1,j,k,5) = 0.0
   20 continue
c
      do 30 j=jmn,jmp1
      do 30 i=imn,imp1
        rs(i,j,kmn,1)  = 0.0
        rs(i,j,kmn,2)  = 0.0
        rs(i,j,kmn,3)  = 0.0
        rs(i,j,kmn,4)  = 0.0
        rs(i,j,kmn,5)  = 0.0
c
        rs(i,j,kmp1,1) = 0.0
        rs(i,j,kmp1,2) = 0.0
        rs(i,j,kmp1,3) = 0.0
        rs(i,j,kmp1,4) = 0.0
        rs(i,j,kmp1,5) = 0.0
   30 continue
c
      return
      end
c
c   ********** change conservative to primitive variables
c
      subroutine  res_cons_to_prim(imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,w,dwn)
c
c     conservative residuals dw changed to nonconservative 
c     residuals
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,prl,prt,rey,suthc,tref,i2dfl,iturb
      dimension w(imp2,jmp2,kmp2,5),dwn(imp2,jmp2,kmp2,5)
c
c
      do 900 k=kmn+1,kmx
      do 900 j=jmn+1,jmx
      do 900 i=imn+1,imx
        res_r  = dwn(i,j,k,1)
        res_ru = dwn(i,j,k,2)
        res_rv = dwn(i,j,k,3)
        res_rw = dwn(i,j,k,4)
        res_re = dwn(i,j,k,5)
c
        rr     = 1./w(i,j,k,1)
        u      = w(i,j,k,2) * rr
        v      = w(i,j,k,3) * rr
        ww     = w(i,j,k,4) * rr
        qq2    = 0.5 * (u*u + v*v + ww*ww)
        res_u  = rr*(res_ru - u *res_r)
        res_v  = rr*(res_rv - v *res_r)
        res_w  = rr*(res_rw - ww*res_r)
        res_p  = gm1*(qq2*res_r - u*res_ru - v*res_rv - ww*res_rw
     .                + res_re )
c
        dwn(i,j,k,1) = res_r
        dwn(i,j,k,2) = res_p
        dwn(i,j,k,3) = res_u
        dwn(i,j,k,4) = res_v
        dwn(i,j,k,5) = res_w
c
  900 continue
c
       return
       end
c
c
c   ********** change primitive to conservative variables
      subroutine res_prim_to_cons(imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,w,dw,dwn)
c
c     nonconservative residuals dw changed to conservative 
c     residuals
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,prl,prt,rey,suthc,tref,i2dfl,iturb
      dimension  w(imp2,jmp2,kmp2,5)
      dimension dw(imp2,jmp2,kmp2,5),dwn(imp2,jmp2,kmp2,5)
c
      do 900 k=kmn+1,kmx
      do 900 j=jmn+1,jmx
      do 900 i=imn+1,imx
        res_r  = dwn(i,j,k,1)
        res_p  = dwn(i,j,k,2)
        res_u  = dwn(i,j,k,3)
        res_v  = dwn(i,j,k,4)
        res_w  = dwn(i,j,k,5)
c
        rho    = w(i,j,k,1)
        rr     = 1./rho
        u      = w(i,j,k,2) * rr
        v      = w(i,j,k,3) * rr
        ww     = w(i,j,k,4) * rr
        qq2    = 0.5 * ( u*u + v*v +ww*ww)
c
        res_ru = rho*res_u + u*res_r
        res_rv = rho*res_v + v*res_r
        res_rw = rho*res_w + ww*res_r
        res_re = res_p/gm1 - qq2*res_r + u*res_ru + v*res_rv + ww*res_rw
c
        dw(i,j,k,1) = res_r
        dw(i,j,k,2) = res_ru
        dw(i,j,k,3) = res_rv
        dw(i,j,k,4) = res_rw
        dw(i,j,k,5) = res_re
c
  900 continue
c
       return
       end
c
c
      subroutine  res_puvt_to_prim(imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,w,p,dwn)
c
c
c   ********** change (p,u,v,w,T) to primitive residuals
c
      dimension w(imp2,jmp2,kmp2,5),dwn(imp2,jmp2,kmp2,5)
      dimension p(imp2,jmp2,kmp2)
c
c
      do 900 k=kmn+1,kmx
      do 900 j=jmn+1,jmx
      do 900 i=imn+1,imx
        res_p = dwn(i,j,k,1)
        res_u = dwn(i,j,k,2)
        res_v = dwn(i,j,k,3)
        res_w = dwn(i,j,k,4)
        res_T = dwn(i,j,k,5)
        
c
        rho   = w(i,j,k,1)
        overT = rho/p(i,j,k)
        res_r = (res_p - rho*res_T)*overT
c
        dwn(i,j,k,1) = res_r
        dwn(i,j,k,2) = res_p
        dwn(i,j,k,3) = res_u
        dwn(i,j,k,4) = res_v
        dwn(i,j,k,5) = res_w
c
  900 continue
c
       return
       end
c
c

      subroutine res_prim_to_puvt(imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,w,p,dw,dwn)
c
c   ********** change primitive to (p,u,v,w,T) residuals
c
      dimension  w(imp2,jmp2,kmp2,5),  p(imp2,jmp2,kmp2)
      dimension dw(imp2,jmp2,kmp2,5),dwn(imp2,jmp2,kmp2,5)
c
      do 900 k=kmn+1,kmx
      do 900 j=jmn+1,jmx
      do 900 i=imn+1,imx
        res_r  = dwn(i,j,k,1)
        res_p  = dwn(i,j,k,2)
        res_u  = dwn(i,j,k,3)
        res_v  = dwn(i,j,k,4)
        res_w  = dwn(i,j,k,5)
c
        rr     = 1./w(i,j,k,1)
        T      = p(i,j,k)*rr
        res_T  = (res_p - T*res_r)*rr
c
        dw(i,j,k,1) = res_p
        dw(i,j,k,2) = res_u
        dw(i,j,k,3) = res_v
        dw(i,j,k,4) = res_w
        dw(i,j,k,5) = res_T
c
  900 continue
c
       return
       end
c
       subroutine solve4_rs(f,b)
c
c          two dimensional
c
c	solve 4x4 system b*x = f
c	replaces x in place of f
c
cc      assumes b(2,1)=b(3,1)=b(4,1) = 0
c           leaves out #5 position in 5x5 matrix
c
        dimension f(5), b(5,5)

        b(1,1) = 1./b(1,1)
        b(1,2) = b(1,2)*b(1,1)
        b(1,3) = b(1,3)*b(1,1)
        b(1,4) = b(1,4)*b(1,1)
c
        b(2,2) = 1./(b(2,2)-b(2,1)*b(1,2))
        b(2,3) = (b(2,3)-b(2,1)*b(1,3))*b(2,2)
        b(2,4) = (b(2,4)-b(2,1)*b(1,4))*b(2,2)
        b(3,2) = b(3,2)-b(3,1)*b(1,2)
        b(3,3) = 1./(b(3,3)-b(3,1)*b(1,3)-b(3,2)*b(2,3))
        b(3,4) = (b(3,4)-b(3,1)*b(1,4)-b(3,2)*b(2,4))*b(3,3)
        b(4,2) = b(4,2)-b(4,1)*b(1,2)
        b(4,3) = b(4,3)-b(4,1)*b(1,3)-b(4,2)*b(2,3)
        b(4,4) = 1./(b(4,4)-b(4,1)*b(1,4)-b(4,2)*b(2,4)-b(4,3)*b(3,4))
c
c		fwd substitution
c
 	d1     = f(1)*b(1,1)
        d2     = (f(2)-b(2,1)*d1)*b(2,2)
        d3     = (f(3)-b(3,1)*d1-b(3,2)*d2)*b(3,3)
        d4     = (f(4)-b(4,1)*d1-b(4,2)*d2-b(4,3)*d3)*b(4,4)
c
c		bwd substitution
c
        f(5)   = 0.
 	f(4)   = d4
        f(3)   = d3-b(3,4)*f(4)
        f(2)   = d2-b(2,3)*f(3)-b(2,4)*f(4)
        f(1)   = d1-b(1,2)*f(2)-b(1,3)*f(3)-b(1,4)*f(4)
c
        return
c
        end

c
      subroutine solve5_rs(f,b)
c
c          three dimensional
c      assumes b(2,1)=b(3,1)=b(4,1) = b(5,1) = 0
c
      dimension f(5), b(5,5)

      b(1,1)     = 1./b(1,1)
      b(1,2)     = b(1,2)*b(1,1)
      b(1,3)     = b(1,3)*b(1,1)
      b(1,4)     = b(1,4)*b(1,1)
      b(1,5)     = b(1,5)*b(1,1)
c
      b(2,2)     = 1./(b(2,2) -b(2,1)*b(1,2))
      b(2,3)     = (b(2,3)-b(2,1)*b(1,3))*b(2,2)
      b(2,4)     = (b(2,4)-b(2,1)*b(1,4))*b(2,2)
      b(2,5)     = (b(2,5)-b(2,1)*b(1,5))*b(2,2)

      b(3,2)     = b(3,2)-b(3,1)*b(1,2)
      b(3,3)     = 1./(b(3,3)-b(3,1)*b(1,3)-b(3,2)*b(2,3))
      b(3,4)     = (b(3,4)-b(3,1)*b(1,4)-b(3,2)*b(2,4))*b(3,3)
      b(3,5)     = (b(3,5)-b(3,1)*b(1,5)-b(3,2)*b(2,5))*b(3,3)

      b(4,2)     = b(4,2)-b(4,1)*b(1,2)
      b(4,3)     = b(4,3)-b(4,1)*b(1,3)-b(4,2)*b(2,3)
      b(4,4)     = 1./(b(4,4)-b(4,1)*b(1,4)-b(4,2)*b(2,4)
     .              -b(4,3)*b(3,4))
      b(4,5)     = (b(4,5)-b(4,1)*b(1,5)-b(4,2)*b(2,5)
     .                    -b(4,3)*b(3,5))*b(4,4)

      b(5,2)     = b(5,2)-b(5,1)*b(1,2)
      b(5,3)     = b(5,3)-b(5,1)*b(1,3)-b(5,2)*b(2,3)
      b(5,4)     = b(5,4)-b(5,1)*b(1,4)-b(5,2)*b(2,4)-b(5,3)*b(3,4)
      b(5,5)     = 1./(b(5,5)-b(5,1)*b(1,5)-b(5,2)*b(2,5)
     .                       -b(5,3)*b(3,5)-b(5,4)*b(4,5))
c
c            fwd substitution
c
      d1         = f(1)*b(1,1)
      d2         = (f(2)-b(2,1)*d1)*b(2,2)
      d3         = (f(3)-b(3,1)*d1-b(3,2)*d2)*b(3,3)
      d4         = (f(4)-b(4,1)*d1-b(4,2)*d2-b(4,3)*d3)*b(4,4)
      d5         = (f(5)-b(5,1)*d1-b(5,2)*d2-b(5,3)*d3-b(5,4)*d4)*b(5,5)
c
c            bwd substitution
c
      f(5)       = d5
      f(4)       = d4-b(4,5)*f(5)
      f(3)       = d3-b(3,4)*f(4)-b(3,5)*f(5)
      f(2)       = d2-b(2,3)*f(3)-b(2,4)*f(4)-b(2,5)*f(5)
      f(1)       = d1-b(1,2)*f(2)-b(1,3)*f(3)-b(1,4)*f(4)-b(1,5)*f(5)
c
      return
c
      end
