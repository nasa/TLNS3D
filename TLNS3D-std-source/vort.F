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
      subroutine  vort (imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2,w,p,eomu,vol,x,
     .                  six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                  duxi,dvxi,dwxi,duyi,dvyi,dwyi,duzi,dvzi,dwzi,
     .                  duxj,dvxj,dwxj,duyj,dvyj,dwyj,duzj,dvzj,dwzj,
     .                  duxm,dvxm,dwxm,duym,dvym,dwym,duzm,dvzm,dwzm,
     .                  duxk,dvxk,dwxk,duyk,dvyk,dwyk,duzk,dvzk,dwzk,
     .                  up,vp,wp,ts,amu,fv,ivisb,iturbb             )
c=======================================================================
c 
c 
c       object : compute 3-D vorticity distribution for
c                turbulence model; also initialize eddy-viscosity
c       initial coding : by v.n.vatsa of nasa langley (april -1991)
c 
c 
c 

      dimension  ivisb(3)
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension 
     .    w(imp2,jmp2,kmp2,5),p(imp2,jmp2,kmp2),vol(imp2,jmp2,kmp2),
     .    x(imp1,jmp1,kmp1,3),fv(imp2,jmp2,kmp2,4),
     .   up(imp1,jmp1,kmp1), vp(imp1,jmp1,kmp1), wp(imp1,jmp1,kmp1),
     .   ts(imp1,jmp1,kmp1),amu(imp1,jmp1,kmp1),eomu(imp2,jmp2,kmp2),
     .  six(imp1,jmp1,kmp1),siy(imp1,jmp1,kmp1),siz(imp1,jmp1,kmp1),
     .  sjx(imp1,jmp1,kmp1),sjy(imp1,jmp1,kmp1),sjz(imp1,jmp1,kmp1),
     .  skx(imp1,jmp1,kmp1),sky(imp1,jmp1,kmp1),skz(imp1,jmp1,kmp1),
     .  dvxi(imp1),dwxi(imp1),duyi(imp1),dwyi(imp1),duzi(imp1),
     .  dvxj(imp1),dwxj(imp1),duyj(imp1),dwyj(imp1),duzj(imp1),
     .  dvxm(imp1),dwxm(imp1),duym(imp1),dwym(imp1),duzm(imp1),
     .  dvzi(imp1),dvzj(imp1),dvzm(imp1),
     .  duxi(imp1),dvyi(imp1),dwzi(imp1),
     .  duxj(imp1),dvyj(imp1),dwzj(imp1),
     .  duxm(imp1),dvym(imp1),dwzm(imp1),
     .  dvxk(imp1,jmp1),dwxk(imp1,jmp1),duyk(imp1,jmp1),
     .  dwyk(imp1,jmp1),duzk(imp1,jmp1),dvzk(imp1,jmp1),
     .  duxk(imp1,jmp1),dvyk(imp1,jmp1),dwzk(imp1,jmp1) 
c
c
c
      imnp = imn + 1
      jmnp = jmn + 1
      kmnp = kmn + 1
c 
c--------- compute velocities and molecular viscosity
c
      do 10 k=kmn,kmp1
      do 10 j=jmn,jmp1
      do 10 i=imn,imp1
      up(i,j,k) =  w(i,j,k,2)/w(i,j,k,1)
      vp(i,j,k) =  w(i,j,k,3)/w(i,j,k,1)
      wp(i,j,k) =  w(i,j,k,4)/w(i,j,k,1)
      ts(i,j,k) =  p(i,j,k)/w(i,j,k,1) 
      amu(i,j,k)= ts(i,j,k) 
      fv(i,j,k,1) = 0.
c
c     initialize fv(i,j,k,2) for smin in spalart model
c
      fv(i,j,k,2) = 1.e+20
   10 continue
c
      if(suthc.gt.0.) then
c
        visc1     = suthc/tref
        do 6 k=kmn,kmp1
        do 6 j=jmn,jmp1
        do 6 i=imn,imp1
        amu(i,j,k)=((1.+visc1)/(ts(i,j,k)+visc1))*ts(i,j,k)**1.5
    6   continue
      endif
c
c     initialize eddy-viscosity for Baldwin-Lomax model
c
      if (iturb.le.1) then
        do 15 k=kmn-1,kmp2
        do 15 j=jmn-1,jmp2
        do 15 i=imn-1,imp2
        eomu(i,j,k) = 0.
        fv(i,j,k,1) = 0.
        fv(i,j,k,3) = 0.
        fv(i,j,k,4) = 0.
   15   continue
      endif
c 
c********* compute vorticity for later use in turbulence routine
c 
      if(iturbb.ne.1) return
c     
      if(ivisb(1).eq.0.and.ivisb(2).eq.0.and.ivisb(3).eq.0) return
c
      do 20 j=jmnp,jmx
      do 20 i=imnp,imx
      rvol      = 2./(vol(i,j,kmnp)  +vol(i,j,kmn))
      duxk(i,j) = skx(i,j,kmn)*(up(i,j,kmnp) - up(i,j,kmn))*rvol
      dvxk(i,j) = skx(i,j,kmn)*(vp(i,j,kmnp) - vp(i,j,kmn))*rvol
      dwxk(i,j) = skx(i,j,kmn)*(wp(i,j,kmnp) - wp(i,j,kmn))*rvol
      duyk(i,j) = sky(i,j,kmn)*(up(i,j,kmnp) - up(i,j,kmn))*rvol
      dvyk(i,j) = sky(i,j,kmn)*(vp(i,j,kmnp) - vp(i,j,kmn))*rvol
      dwyk(i,j) = sky(i,j,kmn)*(wp(i,j,kmnp) - wp(i,j,kmn))*rvol
      duzk(i,j) = skz(i,j,kmn)*(up(i,j,kmnp) - up(i,j,kmn))*rvol
      dvzk(i,j) = skz(i,j,kmn)*(vp(i,j,kmnp) - vp(i,j,kmn))*rvol
      dwzk(i,j) = skz(i,j,kmn)*(wp(i,j,kmnp) - wp(i,j,kmn))*rvol
   20 continue
c
      do 60 k=kmnp,kmx
      do 30 i=imnp,imx
      rvol      = 2./(vol(i,jmnp,k) + vol(i,jmn,k))
      duxj(i)   = sjx(i,jmn,k)*(up(i,jmnp,k) - up(i,jmn,k))*rvol
      dvxj(i)   = sjx(i,jmn,k)*(vp(i,jmnp,k) - vp(i,jmn,k))*rvol
      dwxj(i)   = sjx(i,jmn,k)*(wp(i,jmnp,k) - wp(i,jmn,k))*rvol
      duyj(i)   = sjy(i,jmn,k)*(up(i,jmnp,k) - up(i,jmn,k))*rvol
      dvyj(i)   = sjy(i,jmn,k)*(vp(i,jmnp,k) - vp(i,jmn,k))*rvol
      dwyj(i)   = sjy(i,jmn,k)*(wp(i,jmnp,k) - wp(i,jmn,k))*rvol
      duzj(i)   = sjz(i,jmn,k)*(up(i,jmnp,k) - up(i,jmn,k))*rvol
      dvzj(i)   = sjz(i,jmn,k)*(vp(i,jmnp,k) - vp(i,jmn,k))*rvol
      dwzj(i)   = sjz(i,jmn,k)*(wp(i,jmnp,k) - wp(i,jmn,k))*rvol
   30 continue
c
      do 60 j=jmnp,jmx
      do 40 i=imn ,imx
      rvol      = 2./(vol(i+1,j,k) + vol(i,j,k))
      duxi(i)   = six(i,j,k)*(up(i+1,j,k) - up(i,j,k))*rvol
      dvxi(i)   = six(i,j,k)*(vp(i+1,j,k) - vp(i,j,k))*rvol
      dwxi(i)   = six(i,j,k)*(wp(i+1,j,k) - wp(i,j,k))*rvol
      duyi(i)   = siy(i,j,k)*(up(i+1,j,k) - up(i,j,k))*rvol
      dvyi(i)   = siy(i,j,k)*(vp(i+1,j,k) - vp(i,j,k))*rvol
      dwyi(i)   = siy(i,j,k)*(wp(i+1,j,k) - wp(i,j,k))*rvol
      duzi(i)   = siz(i,j,k)*(up(i+1,j,k) - up(i,j,k))*rvol
      dvzi(i)   = siz(i,j,k)*(vp(i+1,j,k) - vp(i,j,k))*rvol
      dwzi(i)   = siz(i,j,k)*(wp(i+1,j,k) - wp(i,j,k))*rvol
   40 continue
c
      do 50 i=imnp,imx
      duxm(i)   = duxi(i) + duxi(i-1) + duxj(i) + duxk(i,j)
      dvxm(i)   = dvxi(i) + dvxi(i-1) + dvxj(i) + dvxk(i,j)
      dwxm(i)   = dwxi(i) + dwxi(i-1) + dwxj(i) + dwxk(i,j)
      duym(i)   = duyi(i) + duyi(i-1) + duyj(i) + duyk(i,j)
      dvym(i)   = dvyi(i) + dvyi(i-1) + dvyj(i) + dvyk(i,j)
      dwym(i)   = dwyi(i) + dwyi(i-1) + dwyj(i) + dwyk(i,j)
      duzm(i)   = duzi(i) + duzi(i-1) + duzj(i) + duzk(i,j)
      dvzm(i)   = dvzi(i) + dvzi(i-1) + dvzj(i) + dvzk(i,j)
      dwzm(i)   = dwzi(i) + dwzi(i-1) + dwzj(i) + dwzk(i,j)
   50 continue
c
      do 60 i=imnp,imx
      rvol      = 2./(vol(i,j+1,k) + vol(i,j,k))
      duxj(i)   = sjx(i,j,k)*(up(i,j+1,k) - up(i,j,k))*rvol
      dvxj(i)   = sjx(i,j,k)*(vp(i,j+1,k) - vp(i,j,k))*rvol
      dwxj(i)   = sjx(i,j,k)*(wp(i,j+1,k) - wp(i,j,k))*rvol
      duyj(i)   = sjy(i,j,k)*(up(i,j+1,k) - up(i,j,k))*rvol
      dvyj(i)   = sjy(i,j,k)*(vp(i,j+1,k) - vp(i,j,k))*rvol
      dwyj(i)   = sjy(i,j,k)*(wp(i,j+1,k) - wp(i,j,k))*rvol
      duzj(i)   = sjz(i,j,k)*(up(i,j+1,k) - up(i,j,k))*rvol
      dvzj(i)   = sjz(i,j,k)*(vp(i,j+1,k) - vp(i,j,k))*rvol
      dwzj(i)   = sjz(i,j,k)*(wp(i,j+1,k) - wp(i,j,k))*rvol
c 
      rvol      = 2./(vol(i,j,k+1) + vol(i,j,k))
      duxk(i,j) = skx(i,j,k)*(up(i,j,k+1) - up(i,j,k))*rvol
      dvxk(i,j) = skx(i,j,k)*(vp(i,j,k+1) - vp(i,j,k))*rvol
      dwxk(i,j) = skx(i,j,k)*(wp(i,j,k+1) - wp(i,j,k))*rvol
      duyk(i,j) = sky(i,j,k)*(up(i,j,k+1) - up(i,j,k))*rvol
      dvyk(i,j) = sky(i,j,k)*(vp(i,j,k+1) - vp(i,j,k))*rvol
      dwyk(i,j) = sky(i,j,k)*(wp(i,j,k+1) - wp(i,j,k))*rvol
      duzk(i,j) = skz(i,j,k)*(up(i,j,k+1) - up(i,j,k))*rvol
      dvzk(i,j) = skz(i,j,k)*(vp(i,j,k+1) - vp(i,j,k))*rvol
      dwzk(i,j) = skz(i,j,k)*(wp(i,j,k+1) - wp(i,j,k))*rvol
c
      dudx      = duxm(i) + duxj(i) + duxk(i,j)
      dvdx      = dvxm(i) + dvxj(i) + dvxk(i,j)
      dwdx      = dwxm(i) + dwxj(i) + dwxk(i,j)
      dudy      = duym(i) + duyj(i) + duyk(i,j)
      dvdy      = dvym(i) + dvyj(i) + dvyk(i,j)
      dwdy      = dwym(i) + dwyj(i) + dwyk(i,j)
      dudz      = duzm(i) + duzj(i) + duzk(i,j)
      dvdz      = dvzm(i) + dvzj(i) + dvzk(i,j)
      dwdz      = dwzm(i) + dwzj(i) + dwzk(i,j)
c
c    the derivatives here are twice teir actual values
c
c    vorticiy
c
      vortl  = .5*sqrt((dudy-dvdx)*(dudy-dvdx)+
     .                      (dvdz-dwdy)*(dvdz-dwdy)+
     .                      (dwdx-dudz)*(dwdx-dudz))
c
c    strain tensor coded by vatsa (2-16-93)
c
      strain = 2.*(dudx*dudx + dvdy*dvdy + dwdz*dwdz)
     .       +    (dudy + dvdx)**2 + (dudz + dwdx)**2
     .       +    (dvdz + dwdy)**2
      strain = .5*sqrt(strain)
c
c
ccc   fv(i,j,k,1) = min(vortl,strain)
ccc   fv(i,j,k,1) = vortl
      fv(i,j,k,1) = strain
   60 continue
c 
c
      return
      end 
