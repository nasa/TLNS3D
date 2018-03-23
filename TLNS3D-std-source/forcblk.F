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
      subroutine  forcblk(imn, jmn, kmn, imx, jmx, kmx,
     .                    imp1,jmp1,kmp1,imp2,jmp2,kmp2 ,
     .                    w,p,vol,x,
     .                    six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                  tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .                    ri1,ri2,rj1,rj2,rk1,rk2,
     .                    mx1dwk, cp, xp, yp, zp,
     .                    nface ,n1beg ,n1end ,n2beg ,n2end ,
     .                    scal, xref, yref, zref,
     .                   cx,cy,cz,cxv,cyv,czv,cmx,cmy,cmz,nbctype,ivisb)
c=======================================================================
c 
c     object : to compute force/moment coefficients on specified
c              face of a given block
c
c     initial coding : by v.n.vatsa of nasa langley (April 1990)
c 
c     latest revision: by v.n.vatsa (Jan. 1998)
c     (To account for forces due to yaw angle: after Oh's work)
c
c
      dimension      w(imp2,jmp2,kmp2,5), vol(imp2,jmp2,kmp2),
     .               x(imp1,jmp1,kmp1,3),   p(imp2,jmp2,kmp2)
c
      dimension      six(imp1,jmp1,kmp1), siy(imp1,jmp1,kmp1),
     .               siz(imp1,jmp1,kmp1), sjx(imp1,jmp1,kmp1),
     .               sjy(imp1,jmp1,kmp1), sjz(imp1,jmp1,kmp1),
     .               skx(imp1,jmp1,kmp1), sky(imp1,jmp1,kmp1),
     .               skz(imp1,jmp1,kmp1)
c
      dimension      ri1(jmp1,kmp1,3)   , ri2(jmp1,kmp1,3),
     .               rj1(imp1,kmp1,3)   , rj2(imp1,kmp1,3),
     .               rk1(imp1,jmp1,3)   , rk2(imp1,jmp1,3)
c
      dimension    tauwfi1(jmp1,kmp1)  , tauwfi2(jmp1,kmp1),
     .             tauwfj1(imp1,kmp1)  , tauwfj2(imp1,kmp1),
     .             tauwfk1(imp1,jmp1)  , tauwfk2(imp1,jmp1)
c
c
      dimension      cp(mx1dwk,mx1dwk) , xp(mx1dwk,mx1dwk),
     .               yp(mx1dwk,mx1dwk) , zp(mx1dwk,mx1dwk)
c
      dimension      ivisb(3)
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c    
      sgrmi                = 1./(rm*sqrt(gamma))
      gmm                  = gamma*rm*rm
      constv               = 0.
      if(ivisb(1).gt.0.or.ivisb(2).gt.0.or.ivisb(3).gt.0)
     . constv = 4.*sgrmi/rey
      isuthl  = 0
      if(suthc.gt.0.) isuthl = 1
      visc1     = suthc/tref
c
c
c     branch-cut for different type of boundary face
c
      go to (100,200,300,400,500,600) nface
c
  100 continue
c--------------------- i = imn face ------------------------------------
c
      i    =  imn
c
c
      do 110 k = n2beg, n2end
      do 110 j = n1beg, n1end
      cp(j,k)   = (p(i,j,k) + p(i+1,j,k) - 2.*p0)/gmm
  110 continue
c
c
      do 150 k = n2beg+1, n2end
c
      do 130 j = n1beg+1, n1end
      xp(j,k)   = .25*scal*(x(i,j  ,k,1) + x(i,j  ,k-1,1) 
     .                     +x(i,j-1,k,1) + x(i,j-1,k-1,1))
      yp(j,k)   = .25*scal*(x(i,j  ,k,2) + x(i,j  ,k-1,2) 
     .                     +x(i,j-1,k,2) + x(i,j-1,k-1,2))
      zp(j,k)   = .25*scal*(x(i,j  ,k,3) + x(i,j  ,k-1,3) 
     .                     +x(i,j-1,k,3) + x(i,j-1,k-1,3))
  130 continue
c
        dcxv      = 0.
        dcyv      = 0.
        dczv      = 0.
c
        do 140 j = n1beg+1, n1end
        dcx       = -cp(j,k)*six(i,j,k) 
        dcy       = -cp(j,k)*siy(i,j,k) 
        dcz       = -cp(j,k)*siz(i,j,k) 
c
c     add viscous contribution to body forces
c
        if(nbctype .eq. 3) then
          tstat     =  p(i+1,j,k)/w(i+1,j,k,1)
c
c         use Linear or Sutherland's law based on input
c
c         amuw      =  tstat
          amuw    = (1.-isuthl)*tstat + 
     .              isuthl*((1.+visc1)/(tstat+visc1))*tstat**1.5
c
          qn        =  w(i+1,j,k,2)*ri1(j,k,1)
     .              +  w(i+1,j,k,3)*ri1(j,k,2)
     .              +  w(i+1,j,k,4)*ri1(j,k,3)
          sisq  =  six(i,j,k)**2+siy(i,j,k)**2+siz(i,j,k)**2
          sipsq =  six(i+1,j,k)**2+siy(i+1,j,k)**2+siz(i+1,j,k)**2
          si    =  sqrt(sisq)
          siav  =  .5* ( si + sqrt(sipsq) )
          const =  constv*amuw*si*siav/(vol(i+1,j,k)*w(i+1,j,k,1))
c
c     a multiplicative factor of 2 was used in constv since vol/siav
c     gives the distance to node center ,which is twice of 
c     distance to cell center
c     the other 2 comes from definition of force coefficients
c     and definition of reference area
c
          dcxv      =  const*(w(i+1,j,k,2) - qn*ri1(j,k,1))
          dcyv      =  const*(w(i+1,j,k,3) - qn*ri1(j,k,2))
          dczv      =  const*(w(i+1,j,k,4) - qn*ri1(j,k,3))
        end if
c
        cxv       = cxv + dcxv 
        cyv       = cyv + dcyv
        czv       = czv + dczv
        dcx       = dcx + dcxv
        dcy       = dcy + dcyv
        dcz       = dcz + dczv
        cx        = cx  + dcx 
        cy        = cy  + dcy 
        cz        = cz  + dcz 
        cmx       = cmx + dcy*(zp(j,k) -zref)  -dcz*(yp(j,k) -yref) 
        cmy       = cmy + dcz*(xp(j,k) -xref)  -dcx*(zp(j,k) -zref)
        cmz       = cmz + dcx*(yp(j,k) -yref)  -dcy*(xp(j,k) -xref) 
  140 continue
c
  150 continue
      return
c
  200 continue
c--------------------- i = imx face ------------------------------------
c
      i    =  imx
c
c
      do 210 k = n2beg, n2end
      do 210 j = n1beg, n1end
      cp(j,k)   = (p(i,j,k) + p(i+1,j,k) - 2.*p0)/gmm
  210 continue
c
c
      do 250 k = n2beg+1, n2end
c
      do 230 j = n1beg+1, n1end
      xp(j,k)   = .25*scal*(x(i,j  ,k,1) + x(i,j  ,k-1,1) 
     .                     +x(i,j-1,k,1) + x(i,j-1,k-1,1))
      yp(j,k)   = .25*scal*(x(i,j  ,k,2) + x(i,j  ,k-1,2) 
     .                     +x(i,j-1,k,2) + x(i,j-1,k-1,2))
      zp(j,k)   = .25*scal*(x(i,j  ,k,3) + x(i,j  ,k-1,3) 
     .                     +x(i,j-1,k,3) + x(i,j-1,k-1,3))
  230 continue
c
        dcxv      = 0.
        dcyv      = 0.
        dczv      = 0.
c
        do 240 j = n1beg+1, n1end
        dcx       =  cp(j,k)*six(i,j,k) 
        dcy       =  cp(j,k)*siy(i,j,k) 
        dcz       =  cp(j,k)*siz(i,j,k) 
c
c     add viscous contribution to body forces
c
        if(nbctype .eq. 3) then
          tstat     =  p(i,j,k)/w(i,j,k,1)
c
c
c         amuw      =  tstat
          amuw    = (1.-isuthl)*tstat + 
     .              isuthl*((1.+visc1)/(tstat+visc1))*tstat**1.5
c
          qn        =  w(i,j,k,2)*ri2(j,k,1)
     .              +  w(i,j,k,3)*ri2(j,k,2)
     .              +  w(i,j,k,4)*ri2(j,k,3)
          sisq  =  six(i,j,k)**2+siy(i,j,k)**2+siz(i,j,k)**2
          simsq =  six(i-1,j,k)**2+siy(i-1,j,k)**2+siz(i-1,j,k)**2
          si    =  sqrt(sisq)
          siav  =  .5* ( si + sqrt(simsq) )
          const =  constv*amuw*si*siav/(vol(i,j,k)*w(i,j,k,1))
c
          dcxv      =  const*(w(i,j,k,2) - qn*ri2(j,k,1))
          dcyv      =  const*(w(i,j,k,3) - qn*ri2(j,k,2))
          dczv      =  const*(w(i,j,k,4) - qn*ri2(j,k,3))
        end if
c
        cxv       = cxv + dcxv 
        cyv       = cyv + dcyv
        czv       = czv + dczv
        dcx       = dcx + dcxv
        dcy       = dcy + dcyv
        dcz       = dcz + dczv
        cx        = cx  + dcx 
        cy        = cy  + dcy 
        cz        = cz  + dcz 
        cmx       = cmx + dcy*(zp(j,k) -zref)  -dcz*(yp(j,k) -yref) 
        cmy       = cmy + dcz*(xp(j,k) -xref)  -dcx*(zp(j,k) -zref)
        cmz       = cmz + dcx*(yp(j,k) -yref)  -dcy*(xp(j,k) -xref) 
  240 continue
c
  250 continue
      return
c
  300 continue
c--------------------- j = jmn face ------------------------------------
c
      j    =  jmn
c
c
      do 310 k = n1beg, n1end
      do 310 i = n2beg, n2end
      cp(i,k)   = (p(i,j,k) + p(i,j+1,k) - 2.*p0)/gmm
  310 continue
c
c
      do 350 k = n1beg+1, n1end
c
      do 330 i = n2beg+1, n2end
      xp(i,k)   = .25*scal*(x(i  ,j,k,1) + x(i  ,j,k-1,1) 
     .                     +x(i-1,j,k,1) + x(i-1,j,k-1,1))
      yp(i,k)   = .25*scal*(x(i  ,j,k,2) + x(i  ,j,k-1,2) 
     .                     +x(i-1,j,k,2) + x(i-1,j,k-1,2))
      zp(i,k)   = .25*scal*(x(i  ,j,k,3) + x(i  ,j,k-1,3) 
     .                     +x(i-1,j,k,3) + x(i-1,j,k-1,3))
  330 continue
c
        dcxv      = 0.
        dcyv      = 0.
        dczv      = 0.
c
        do 340 i = n2beg+1, n2end
        dcx       = -cp(i,k)*sjx(i,j,k) 
        dcy       = -cp(i,k)*sjy(i,j,k) 
        dcz       = -cp(i,k)*sjz(i,j,k) 
c
c     add viscous contribution to body forces
c
        if(nbctype .eq. 3) then
          tstat     =  p(i,j+1,k)/w(i,j+1,k,1)
c
c
c         amuw      =  tstat
          amuw    = (1.-isuthl)*tstat + 
     .              isuthl*((1.+visc1)/(tstat+visc1))*tstat**1.5
c
          qn        =  w(i,j+1,k,2)*rj1(i,k,1)
     .              +  w(i,j+1,k,3)*rj1(i,k,2)
     .              +  w(i,j+1,k,4)*rj1(i,k,3)
          sjsq  =  sjx(i,j,k)**2+sjy(i,j,k)**2+sjz(i,j,k)**2
          sjpsq =  sjx(i,j+1,k)**2+sjy(i,j+1,k)**2+sjz(i,j+1,k)**2
          sj    =  sqrt(sjsq)
          sjav  =  .5* ( sj + sqrt(sjpsq) )
          const =  constv*amuw*sj*sjav/(vol(i,j+1,k)*w(i,j+1,k,1))
c
          dcxv      =  const*(w(i,j+1,k,2) - qn*rj1(i,k,1))
          dcyv      =  const*(w(i,j+1,k,3) - qn*rj1(i,k,2))
          dczv      =  const*(w(i,j+1,k,4) - qn*rj1(i,k,3))
        end if
c
c     add viscous contribution to body forces for wall-function
c
        if(nbctype .eq. 102) then
c
          uslip   = w(i,j,k,2)/w(i,j,k,1)
          vslip   = w(i,j,k,3)/w(i,j,k,1)
          wslip   = w(i,j,k,4)/w(i,j,k,1)
          vtslip  = sqrt(uslip**2 + vslip**2 + wslip**2)
          vtslip  = max (1.e-10,vtslip)
c
          sjsq  =  sjx(i,j,k)**2+sjy(i,j,k)**2+sjz(i,j,k)**2
          sj    =  sqrt(sjsq)
          const =  0.5*constv*sj*tauwfj1(i,k)/vtslip
c
          dcxv      =  const*uslip
          dcyv      =  const*vslip
          dczv      =  const*wslip
        end if
c
c
        cxv       = cxv + dcxv 
        cyv       = cyv + dcyv
        czv       = czv + dczv
        dcx       = dcx + dcxv
        dcy       = dcy + dcyv
        dcz       = dcz + dczv
        cx        = cx  + dcx 
        cy        = cy  + dcy 
        cz        = cz  + dcz 
        cmx       = cmx + dcy*(zp(i,k) -zref)  -dcz*(yp(i,k) -yref) 
        cmy       = cmy + dcz*(xp(i,k) -xref)  -dcx*(zp(i,k) -zref)
        cmz       = cmz + dcx*(yp(i,k) -yref)  -dcy*(xp(i,k) -xref) 
  340 continue
  350 continue
      return
c
  400 continue
c--------------------- j = jmx face ------------------------------------
c
      j    =  jmx
c
c
      do 410 k = n1beg, n1end
      do 410 i = n2beg, n2end
      cp(i,k)   = (p(i,j,k) + p(i,j+1,k) - 2.*p0)/gmm
  410 continue
c
c
      do 450 k = n1beg+1, n1end
c
      do 430 i = n2beg+1, n2end
      xp(i,k)   = .25*scal*(x(i  ,j,k,1) + x(i  ,j,k-1,1) 
     .                     +x(i-1,j,k,1) + x(i-1,j,k-1,1))
      yp(i,k)   = .25*scal*(x(i  ,j,k,2) + x(i  ,j,k-1,2) 
     .                     +x(i-1,j,k,2) + x(i-1,j,k-1,2))
      zp(i,k)   = .25*scal*(x(i  ,j,k,3) + x(i  ,j,k-1,3) 
     .                     +x(i-1,j,k,3) + x(i-1,j,k-1,3))
  430 continue
c
        dcxv      = 0.
        dcyv      = 0.
        dczv      = 0.
c
        do 440 i = n2beg+1, n2end
        dcx       =  cp(i,k)*sjx(i,j,k) 
        dcy       =  cp(i,k)*sjy(i,j,k) 
        dcz       =  cp(i,k)*sjz(i,j,k) 
c
c     add viscous contribution to body forces
c
        if(nbctype .eq. 3) then
          tstat     =  p(i,j,k)/w(i,j,k,1)
c
c
c         amuw      =  tstat
          amuw    = (1.-isuthl)*tstat + 
     .              isuthl*((1.+visc1)/(tstat+visc1))*tstat**1.5
c
          qn        =  w(i,j,k,2)*rj2(i,k,1)
     .              +  w(i,j,k,3)*rj2(i,k,2)
     .              +  w(i,j,k,4)*rj2(i,k,3)
          sjsq  =  sjx(i,j,k)**2+sjy(i,j,k)**2+sjz(i,j,k)**2
          sjmsq =  sjx(i,j-1,k)**2+sjy(i,j-1,k)**2+sjz(i,j-1,k)**2
          sj    =  sqrt(sjsq)
          sjav  =  .5* ( sj + sqrt(sjmsq) )
          const =  constv*amuw*sj*sjav/(vol(i,j,k)*w(i,j,k,1))
c
          dcxv      =  const*(w(i,j,k,2) - qn*rj2(i,k,1))
          dcyv      =  const*(w(i,j,k,3) - qn*rj2(i,k,2))
          dczv      =  const*(w(i,j,k,4) - qn*rj2(i,k,3))
        end if
c
        cxv       = cxv + dcxv 
        cyv       = cyv + dcyv
        czv       = czv + dczv
        dcx       = dcx + dcxv
        dcy       = dcy + dcyv
        dcz       = dcz + dczv
        cx        = cx  + dcx 
        cy        = cy  + dcy 
        cz        = cz  + dcz 
        cmx       = cmx + dcy*(zp(i,k) -zref)  -dcz*(yp(i,k) -yref) 
        cmy       = cmy + dcz*(xp(i,k) -xref)  -dcx*(zp(i,k) -zref)
        cmz       = cmz + dcx*(yp(i,k) -yref)  -dcy*(xp(i,k) -xref) 
  440 continue
  450 continue
      return
c
  500 continue
c--------------------- k = kmn face ------------------------------------
c
      k    =  kmn
c
c
      do 510 j = n2beg, n2end
      do 510 i = n1beg, n1end
      cp(i,j)   = (p(i,j,k) + p(i,j,k+1) - 2.*p0)/gmm
  510 continue
c
c
      do 550 j = n2beg+1, n2end
c
      do 530 i = n1beg+1, n1end
      xp(i,j)   = .25*scal*(x(i  ,j,k,1) + x(i  ,j-1,k,1) 
     .                     +x(i-1,j,k,1) + x(i-1,j-1,k,1))
      yp(i,j)   = .25*scal*(x(i  ,j,k,2) + x(i  ,j-1,k,2) 
     .                     +x(i-1,j,k,2) + x(i-1,j-1,k,2))
      zp(i,j)   = .25*scal*(x(i  ,j,k,3) + x(i  ,j-1,k,3) 
     .                     +x(i-1,j,k,3) + x(i-1,j-1,k,3))
  530 continue
c
        dcxv      = 0.
        dcyv      = 0.
        dczv      = 0.
c
        do 540 i = n1beg+1, n1end
        dcx       = -cp(i,j)*skx(i,j,k) 
        dcy       = -cp(i,j)*sky(i,j,k) 
        dcz       = -cp(i,j)*skz(i,j,k) 
c
c     add viscous contribution to body forces
c
        if(nbctype .eq. 3) then
          tstat     =  p(i,j,k+1)/w(i,j,k+1,1)
c
c
c         amuw      =  tstat
          amuw    = (1.-isuthl)*tstat + 
     .              isuthl*((1.+visc1)/(tstat+visc1))*tstat**1.5
c
          qn        =  w(i,j,k+1,2)*rk1(i,j,1)
     .              +  w(i,j,k+1,3)*rk1(i,j,2)
     .              +  w(i,j,k+1,4)*rk1(i,j,3)
          sksq  =  skx(i,j,k)**2+sky(i,j,k)**2+skz(i,j,k)**2
          skpsq =  skx(i,j,k+1)**2+sky(i,j,k+1)**2+skz(i,j,k+1)**2
          sk    =  sqrt(sksq)
          skav  =  .5* ( sk + sqrt(skpsq) )
          const =  constv*amuw*sk*skav/(vol(i,j,k+1)*w(i,j,k+1,1))
c
          dcxv      =  const*(w(i,j,k+1,2) - qn*rk1(i,j,1))
          dcyv      =  const*(w(i,j,k+1,3) - qn*rk1(i,j,2))
          dczv      =  const*(w(i,j,k+1,4) - qn*rk1(i,j,3))
        end if
c
        cxv       = cxv + dcxv 
        cyv       = cyv + dcyv
        czv       = czv + dczv
        dcx       = dcx + dcxv
        dcy       = dcy + dcyv
        dcz       = dcz + dczv
        cx        = cx  + dcx 
        cy        = cy  + dcy 
        cz        = cz  + dcz 
        cmx       = cmx + dcy*(zp(i,j) -zref)  -dcz*(yp(i,j) -yref) 
        cmy       = cmy + dcz*(xp(i,j) -xref)  -dcx*(zp(i,j) -zref)
        cmz       = cmz + dcx*(yp(i,j) -yref)  -dcy*(xp(i,j) -xref) 
  540 continue
  550 continue
      return
c
  600 continue
c--------------------- k = kmx face ------------------------------------
c
      k    =  kmx
c
c
      do 610 j = n2beg, n2end
      do 610 i = n1beg, n1end
      cp(i,j)   = (p(i,j,k) + p(i,j,k+1) - 2.*p0)/gmm
  610 continue
c
c
      do 650 j = n2beg+1, n2end
c
      do 630 i = n1beg+1, n1end
      xp(i,j)   = .25*scal*(x(i  ,j,k,1) + x(i  ,j-1,k,1) 
     .                     +x(i-1,j,k,1) + x(i-1,j-1,k,1))
      yp(i,j)   = .25*scal*(x(i  ,j,k,2) + x(i  ,j-1,k,2) 
     .                     +x(i-1,j,k,2) + x(i-1,j-1,k,2))
      zp(i,j)   = .25*scal*(x(i  ,j,k,3) + x(i  ,j-1,k,3) 
     .                     +x(i-1,j,k,3) + x(i-1,j-1,k,3))
  630 continue
c
        dcxv      = 0.
        dcyv      = 0.
        dczv      = 0.
c
        do 640 i = n1beg+1, n1end
        dcx       =  cp(i,j)*skx(i,j,k) 
        dcy       =  cp(i,j)*sky(i,j,k) 
        dcz       =  cp(i,j)*skz(i,j,k) 
c
c     add viscous contribution to body forces
c
        if(nbctype .eq. 3) then
          tstat     =  p(i,j,k)/w(i,j,k,1)
c
c
c         amuw      =  tstat
          amuw    = (1.-isuthl)*tstat + 
     .              isuthl*((1.+visc1)/(tstat+visc1))*tstat**1.5
c
          qn        =  w(i,j,k,2)*rk2(i,j,1)
     .              +  w(i,j,k,3)*rk2(i,j,2)
     .              +  w(i,j,k,4)*rk2(i,j,3)
          sksq  =  skx(i,j,k)**2+sky(i,j,k)**2+skz(i,j,k)**2
          skmsq =  skx(i,j,k-1)**2+sky(i,j,k-1)**2+skz(i,j,k-1)**2
          sk    =  sqrt(sksq)
          skav  =  .5* ( sk + sqrt(skmsq) )
          const =  constv*amuw*sk*skav/(vol(i,j,k)*w(i,j,k,1))
c
          dcxv      =  const*(w(i,j,k,2) - qn*rk2(i,j,1))
          dcyv      =  const*(w(i,j,k,3) - qn*rk2(i,j,2))
          dczv      =  const*(w(i,j,k,4) - qn*rk2(i,j,3))
        end if
c
        cxv       = cxv + dcxv 
        cyv       = cyv + dcyv
        czv       = czv + dczv
        dcx       = dcx + dcxv
        dcy       = dcy + dcyv
        dcz       = dcz + dczv
        cx        = cx  + dcx 
        cy        = cy  + dcy 
        cz        = cz  + dcz 
        cmx       = cmx + dcy*(zp(i,j) -zref)  -dcz*(yp(i,j) -yref) 
        cmy       = cmy + dcz*(xp(i,j) -xref)  -dcx*(zp(i,j) -zref)
        cmz       = cmz + dcx*(yp(i,j) -yref)  -dcy*(xp(i,j) -xref) 
  640 continue
  650 continue
      return
c
      end
