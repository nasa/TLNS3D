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
      subroutine  turbmb (imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2,w,p,x,eomu,
     .                  ri1,ri2,rj1,rj2,rk1,rk2,
     .                  xn,yn,zn,fbl,eomui,up,vp,wp,ts,amu,
     .                  mx1dwk, nface ,n1beg ,n1end ,n2beg ,n2end , 
     .                  nbctype,nturb ,n1begt,n1endt,n2begt, n2endt,
     .                  nnstrt ,nnstop,nnendt, fv                  )
c=======================================================================
c
c     original coding: by v.n.vatsa (1991)
c    (Degani-Schiff correction can be activated by setting ids=1)
c     latest mods: by m.d.sanetrik to get identical b.l.edge on min/max faces
c      b.l. edge on min/max faces   (1997)
c

      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      data          aplus/26./,ccp/1.6/,ckleb/.3/,cwk/1./,vk/.4/,
     .              clauser/.0168/
c 
c
      dimension     w(imp2,jmp2,kmp2,5),    p(imp2,jmp2,kmp2),
     .              x(imp1,jmp1,kmp1,3), eomu(imp2,jmp2,kmp2),
     .             fv(imp2,jmp2,kmp2,4),
     .             up(imp1,jmp1,kmp1), vp(imp1,jmp1,kmp1),
     .             wp(imp1,jmp1,kmp1),amu(imp1,jmp1,kmp1),
     .             ts(imp1,jmp1,kmp1),
     .              xn(mx1dwk),   yn(mx1dwk), zn(mx1dwk),
     .             fbl(mx1dwk),eomui(mx1dwk)
c
      dimension     ri1(jmp1,kmp1,3)   ,  ri2(jmp1,kmp1,3),
     .              rj1(imp1,kmp1,3)   ,  rj2(imp1,kmp1,3),
     .              rk1(imp1,jmp1,3)   ,  rk2(imp1,jmp1,3)
c
c     fv(i,j,k,1) is used for storing the 3-d vorticity distribution
c     fv(i,j,k,2) is used for storing normal distance from wall 
c     fv(i,j,k,3) is used for storing sum of 1/fv(i,j,k,2)**2
c     fv(i,j,k,4) is used for storing sum of eomu/fv(i,j,k,2)**2
c 
c     ids       = 1
      ids       = 0
c     ids       = 1
c
      sgrmi     = 1./(rm*sqrt(gamma)) 
      ckout     = rey*sgrmi*clauser*ccp 
c
c     branch-cut for different type of boundary faces
c
      if (nturb.le.0) return 
      go to (100,200,300,400,500,600) nface
c
  100 continue
c--------------------------  i = imn  face   ---------------------------
c
      nmn    = imn
      nmx    = imx
c
c     make sure search/computational range is within the grid
c
      if(nnstrt .lt. (nmn+1) .or. nnstrt .gt. nmx) nnstrt = nmn+1
      if(nnstop .gt. nmx .or. nnstop .lt. (nmn+1)) nnstop = nmx
      if(nnendt .gt. nmx .or. nnendt .lt. (nmn+1)) nnendt = nmx
c
c     set fqs = 0, if imn is a solid wall boundary
c     fqs is used for damping of eddy-viscosity near walls
c
      fqs    = 1.0
      if(nbctype.eq.3) fqs = 0.
c
      do 110 j=n1begt+1,n1endt
      do 120 k=n2begt+1,n2endt
c 
      fv(nmn,j,k,1) = 0.
      fv(nmn,j,k,2) = 0.
c 
      do 122 i=nmn,nmx
      fbl(i)        = fv(i,j,k,1)
      eomu(i,j,k)   = 0.
  122 continue
c
c--------- evaluate normal distance to i=imn wall
c 
      do 125 i=nmn,nmx
      xn(i)       = .25*(x(i,j,k  ,1) + x(i,j-1,k  ,1)
     .                  +x(i,j,k-1,1) + x(i,j-1,k-1,1))
      yn(i)       = .25*(x(i,j,k  ,2) + x(i,j-1,k  ,2)
     .                  +x(i,j,k-1,2) + x(i,j-1,k-1,2))
      zn(i)       = .25*(x(i,j,k  ,3) + x(i,j-1,k  ,3)
     .                  +x(i,j,k-1,3) + x(i,j-1,k-1,3))
  125 continue
      xnw         = xn(nmn)
      ynw         = yn(nmn)
      znw         = zn(nmn)
      do 126 i=nmn+1,nmx
      xnc         = .5*(xn(i) + xn(i-1)) 
      ync         = .5*(yn(i) + yn(i-1)) 
      znc         = .5*(zn(i) + zn(i-1)) 
      fv(i,j,k,2) = abs((xnc-xnw)*ri1(j,k,1) +(ync-ynw)*ri1(j,k,2)
     .                 +(znc-znw)*ri1(j,k,3))
      fv(i,j,k,2) = max  (fv(i,j,k,2),fv(i-1,j,k,2))
  126 continue
c 
c--------- use linear extrapolation to obtain wall values
c 
      vortw     = 1.5*fv(nmn+1,j,k,1)  -.5*fv(nmn+2,j,k,1) 
      rhow      = 1.5*w(nmn+1,j,k,1)   -.5*w(nmn+2,j,k,1) 
      amuw      = 1.5*amu(nmn+1,j,k)   -.5*amu(nmn+2,j,k) 
      tauw      = amuw*vortw
c
      nvortm  = ismax(nnstop-2,fbl(nmn+1),1) + nmn
      vortm   = fbl(nvortm)
      yplusoy   = abs(rhow*vortm*rey*sgrmi/amuw)
      yplusoy   = sqrt(yplusoy) 
      utmax     = 0.
      fbl(nmn)    = 0.
c
      do 130 i=nmn+1,nmx
      yplus     = yplusoy*fv(i,j,k,2) 
      yoa       = (1.-fqs)*yplus/aplus  +fqs*1000.
      yoa       = min  (1000.,yoa)
      damp      = 1.-exp(-yoa)
      fbl(i)    = fv(i,j,k,1)*fv(i,j,k,2)*damp
      amixl     = vk*fv(i,j,k,2)*damp 
      eomui(i)  = rey*sgrmi*w(i,j,k,1)*fv(i,j,k,1)*amixl*amixl/
     .            amu(i,j,k)
      utot      = up(i,j,k)**2  +vp(i,j,k)**2  +wp(i,j,k)**2
      utmax     = max  (utmax,utot)
      eomui(i)  = max  (0.,eomui(i))
  130 continue
c
      fmax      = 0.
      ymax      = 0.
      fmax1     = 0.
      nfmax     = nnstrt
      do 135 i=nnstrt,nnstop
      if (fbl(i).gt.fmax) then
         fmax      = fbl(i)
         ymax      = fv(i,j,k,2)
         fmax1     = .9*fmax
         nfmax     = i
c
c     pickup first peak in fbl - due to degani/schiff
c
      else if (ids .eq. 1 .and. fbl(i) .lt. fmax1) then
         go to 131
      endif
  135 continue
c
  131 continue
      if (nfmax.lt.nnstop) then
      dfm     =   fbl(nfmax) - fbl(nfmax-1)
      dfp     =   fbl(nfmax) - fbl(nfmax+1)
      if(fbl(nfmax-1).lt.fbl(nfmax+1)) then
        ymax  = fv(nfmax,j,k,2) + .5*(fv(nfmax+1,j,k,2)-
     .                               fv(nfmax,j,k,2))*(1.-dfp/dfm)
      else
        ymax  = fv(nfmax,j,k,2) - .5*(fv(nfmax,j,k,2)  -
     .                               fv(nfmax-1,j,k,2))*(1.-dfm/dfp)
      end if
      end if
      utmax       = sqrt(utmax) 
      utmin       = sqrt(up(nmn+1,j,k)**2  +vp(nmn+1,j,k)**2 
     .                                     +wp(nmn+1,j,k)**2 )
      udif        = utmax -utmin
      fwake       = cwk*udif*udif*ymax/fmax 
      fwake       = min  (ymax*fmax,fwake)
      do 140 i=nmn+1,nnendt
      fkleb       = 1./(1.+5.5*(ckleb*fv(i,j,k,2)/ymax)**6 )
      eomuo       = ckout*w(i,j,k,1)/amu(i,j,k)*fkleb*fwake 
      eomu(i,j,k) = max  (0.,eomuo)
  140 continue
c
      do 150 i=nmn+1,nnendt
      nedge       = i
      if (eomui(i).ge.eomu(i,j,k)) go to 160
  150 continue
  160 do 170 i=nmn+1,nedge-1
  170 eomu(i,j,k) = eomui(i)
c 
c--------- set ghost cell values of eddy-viscosity
c 
cc    eomu(nmn,j,k)= eomu(nmn+1,j,k)
cc    if (eomu(nmn,j,k).le.0.) eomu(nmn,j,k) = 0.
c
c     add the contribution of this segment to the summation in
c     fv(i,j,k,3) and fv(i,j,k,4) and update the value of
c     eomu including the contribution from the current segment
c
      do 175 i=nmn+1,nnendt
      fv(i,j,k,3) = fv(i,j,k,3) + 1./(fv(i,j,k,2)**2)
      fv(i,j,k,4) = fv(i,j,k,4) + eomu(i,j,k)/(fv(i,j,k,2)**2)
  175 continue
      do 180 i=nmn+1,nnendt
      eomu(i,j,k) = fv(i,j,k,4)/fv(i,j,k,3)
  180 continue
c
  120 continue
c
c--------- end of outer loop on k
c
  110 continue
c
c--------- end of outer loop on j
c
c
c--------- end eddy-viscosity computations for imn face    -------------
c
      return
c 
  200 continue
c--------------------------  i = imx  face   ---------------------------
c
      nmn    = imn
      nmx    = imx
c
c     make sure search/computational range is within the grid
c
      if(nnstrt.lt.(nmn+1).or.nnstrt.gt.nmx) nnstrt = nmn+1
      if(nnstop.gt.nmx.or.nnstop.lt.(nmn+1)) nnstop = nmx
      if(nnendt.gt.nmx.or.nnendt.lt.(nmn+1)) nnendt = nmx
c
c     reverse the search/computational count for imx face
c
      nnstrt = nmx - nnstrt + nmn + 1
      nnstop = nmx - nnstop + nmn + 1
      nnendt = nmx - nnendt + nmn + 1
c
c     set fqs = 0, if imx is a solid wall boundary
c     fqs is used for damping of eddy-viscosity near walls
c
      fqs    = 1.0
      if(nbctype.eq.3) fqs = 0.
c
      do 210 j=n1begt+1,n1endt
      do 220 k=n2begt+1,n2endt
c 
      fv(nmx+1,j,k,1) = 0.
      fv(nmx+1,j,k,2) = 0.
c 
      do 222 i=nmn,nmx
      fbl(i)        = fv(i,j,k,1)
      eomu(i,j,k)   = 0.
  222 continue
c
c--------- evaluate normal distance to i=imx wall
c 
      do 225 i=nmn,nmx
      xn(i)       = .25*(x(i,j,k  ,1) + x(i,j-1,k  ,1)
     .                  +x(i,j,k-1,1) + x(i,j-1,k-1,1))
      yn(i)       = .25*(x(i,j,k  ,2) + x(i,j-1,k  ,2)
     .                  +x(i,j,k-1,2) + x(i,j-1,k-1,2))
      zn(i)       = .25*(x(i,j,k  ,3) + x(i,j-1,k  ,3)
     .                  +x(i,j,k-1,3) + x(i,j-1,k-1,3))
  225 continue
      xnw         = xn(nmx)
      ynw         = yn(nmx)
      znw         = zn(nmx)
      do 226 i=nmx,nmn+1,-1
      xnc         = .5*(xn(i) + xn(i-1)) 
      ync         = .5*(yn(i) + yn(i-1)) 
      znc         = .5*(zn(i) + zn(i-1)) 
      fv(i,j,k,2) = abs((xnc-xnw)*ri2(j,k,1) +(ync-ynw)*ri2(j,k,2)
     .                 +(znc-znw)*ri2(j,k,3))
      fv(i,j,k,2) = max  (fv(i,j,k,2),fv(i+1,j,k,2))
  226 continue
c 
c--------- use linear extrapolation to obtain wall values
c 
      vortw     = 1.5*fv(nmx,j,k,1)  -.5*fv(nmx-1,j,k,1) 
      rhow      = 1.5*w(nmx,j,k,1)   -.5*w(nmx-1,j,k,1) 
      amuw      = 1.5*amu(nmx,j,k)   -.5*amu(nmx-1,j,k) 
      tauw      = amuw*vortw
c
      nvortm  = ismax(nmx-nnstop+1,fbl(nnstop),1) + nnstop -1
      vortm   = fbl(nvortm)
      yplusoy   = abs(rhow*vortm*rey*sgrmi/amuw)
      yplusoy   = sqrt(yplusoy) 
      utmax     = 0.
      fbl(nmx+1)    = 0.
c
      do 230 i=nmn+1,nmx
      yplus     = yplusoy*fv(i,j,k,2) 
      yoa       = (1.-fqs)*yplus/aplus  +fqs*1000.
      yoa       = min  (1000.,yoa)
      damp      = 1.-exp(-yoa)
      fbl(i)    = fv(i,j,k,1)*fv(i,j,k,2)*damp
      amixl     = vk*fv(i,j,k,2)*damp 
      eomui(i)  = rey*sgrmi*w(i,j,k,1)*fv(i,j,k,1)*amixl*amixl/
     .            amu(i,j,k)
      utot      = up(i,j,k)**2  +vp(i,j,k)**2  +wp(i,j,k)**2
      utmax     = max  (utmax,utot)
      eomui(i)  = max  (0.,eomui(i))
  230 continue
c
      fmax      = 0.
      ymax      = 0.
      fmax1     = 0.
      nfmax     = nnstrt
      do 235 i=nnstrt,nnstop,-1
      if (fbl(i).gt.fmax) then
         fmax      = fbl(i)
         ymax      = fv(i,j,k,2)
         fmax1     = .9*fmax
         nfmax     = i
c
c     pickup first peak in fbl - due to degani/schiff
c
      else if (ids .eq. 1 .and. fbl(i) .lt. fmax1) then
         go to 231
      endif
  235 continue
c
  231 continue
      if (nfmax.gt.nnstop) then
      dfm     =   fbl(nfmax) - fbl(nfmax-1)
      dfp     =   fbl(nfmax) - fbl(nfmax+1)
      if(fbl(nfmax-1).lt.fbl(nfmax+1)) then
        ymax  = fv(nfmax,j,k,2) + .5*(fv(nfmax+1,j,k,2)-
     .                               fv(nfmax,j,k,2))*(1.-dfp/dfm)
      else
        ymax  = fv(nfmax,j,k,2) - .5*(fv(nfmax,j,k,2)  -
     .                               fv(nfmax-1,j,k,2))*(1.-dfm/dfp)
      end if
      end if
      utmax       = sqrt(utmax) 
      utmin       = sqrt(up(nmx,j,k)**2  +vp(nmx,j,k)**2 
     .                                   +wp(nmx,j,k)**2 )
      udif        = utmax -utmin
      fwake       = cwk*udif*udif*ymax/fmax 
      fwake       = min  (ymax*fmax,fwake)
      do 240 i=nmx,nnendt,-1
      fkleb       = 1./(1.+5.5*(ckleb*fv(i,j,k,2)/ymax)**6 )
      eomuo       = ckout*w(i,j,k,1)/amu(i,j,k)*fkleb*fwake 
      eomu(i,j,k) = max  (0.,eomuo)
  240 continue
c
      do 250 i=nmx,nnendt,-1
      nedge       = i
      if (eomui(i).ge.eomu(i,j,k)) go to 260
  250 continue
  260 do 270 i=nmx,nedge+1,-1
  270 eomu(i,j,k) = eomui(i)
c 
c--------- set ghost cell values of eddy-viscosity
c 
cc    eomu(nmx+1,j,k)= eomu(nmx,j,k)
cc    if (eomu(nmx+1,j,k).le.0.) eomu(nmx+1,j,k) = 0.
c
c     add the contribution of this segment to the summation in
c     fv(i,j,k,3) and fv(i,j,k,4) and update the value of
c     eomu including the contribution from the current segment
c
      do 275 i=nmx,nnendt,-1
      fv(i,j,k,3) = fv(i,j,k,3) + 1./(fv(i,j,k,2)**2)
      fv(i,j,k,4) = fv(i,j,k,4) + eomu(i,j,k)/(fv(i,j,k,2)**2)
  275 continue
      do 280 i=nmx,nnendt,-1
      eomu(i,j,k) = fv(i,j,k,4)/fv(i,j,k,3)
  280 continue
c
  220 continue
c
c--------- end of outer loop on k
c
  210 continue
c
c--------- end of outer loop on j
c
c
c--------- end eddy-viscosity computations for imx face    -------------
c
      return
c 
  300 continue
c--------------------------  j = jmn  face   ---------------------------
c
c 
      nmn    = jmn
      nmx    = jmx
c
c     make sure search/computational range is within the grid
c
      if(nnstrt .lt. (nmn+1) .or. nnstrt .gt. nmx) nnstrt = nmn+1
      if(nnstop .gt. nmx .or. nnstop .lt. (nmn+1)) nnstop = nmx
      if(nnendt .gt. nmx .or. nnendt .lt. (nmn+1)) nnendt = nmx
c
c     set fqs = 0, if jmn is a solid wall boundary
c     fqs is used for damping of eddy-viscosity near walls
c
      fqs    = 1.0
      if(nbctype.eq.3) fqs = 0.
c
      do 310 k=n1begt+1,n1endt
      do 320 i=n2begt+1,n2endt
c 
      fv(i,nmn,k,1) = 0.
      fv(i,nmn,k,2) = 0.
c 
      do 322 j=nmn,nmx
      fbl(j)        = fv(i,j,k,1)
      eomu(i,j,k)   = 0.
  322 continue
c
c--------- evaluate normal distance to j=jmn wall
c 
      do 325 j=nmn,nmx
      xn(j)       = .25*(x(i,j,k  ,1) + x(i-1,j,k  ,1)
     .                  +x(i,j,k-1,1) + x(i-1,j,k-1,1))
      yn(j)       = .25*(x(i,j,k  ,2) + x(i-1,j,k  ,2)
     .                  +x(i,j,k-1,2) + x(i-1,j,k-1,2))
      zn(j)       = .25*(x(i,j,k  ,3) + x(i-1,j,k  ,3)
     .                  +x(i,j,k-1,3) + x(i-1,j,k-1,3))
  325 continue
      xnw         = xn(nmn)
      ynw         = yn(nmn)
      znw         = zn(nmn)
      do 326 j=nmn+1,nmx
      xnc         = .5*(xn(j) + xn(j-1)) 
      ync         = .5*(yn(j) + yn(j-1)) 
      znc         = .5*(zn(j) + zn(j-1)) 
c     fv(i,j,k,2) =     (xnc-xnw)*rj1(i,k,1)      +(ync-ynw)*rj1(i,k,2) 
c    .                 +(znc-znw)*rj1(i,k,3) 
      fv(i,j,k,2) = abs((xnc-xnw)*rj1(i,k,1)      +(ync-ynw)*rj1(i,k,2)
     .             +    (znc-znw)*rj1(i,k,3))
c     fv(i,j,k,2) = abs((xnc-xnw)*rj1(i,k,1)) +abs((ync-ynw)*rj1(i,k,2))
c    .             +abs((znc-znw)*rj1(i,k,3))
c     fv(i,j,k,2) = sqrt((xnc-xnw)**2 +(ync-ynw)**2 +(znc-znw)**2)
      fv(i,j,k,2) = max  (fv(i,j,k,2),fv(i,j-1,k,2))
  326 continue
c 
c--------- use linear extrapolation to obtain wall values
c 
      vortw     = 1.5*fv(i,nmn+1,k,1)  -.5*fv(i,nmn+2,k,1) 
      rhow      = 1.5*w(i,nmn+1,k,1)   -.5*w(i,nmn+2,k,1) 
      amuw      = 1.5*amu(i,nmn+1,k)   -.5*amu(i,nmn+2,k) 
      tauw      = amuw*vortw
c
      nvortm  = ismax(nnstop-2,fbl(nmn+1),1) + nmn
      vortm   = fbl(nvortm)
      yplusoy   = abs(rhow*vortm*rey*sgrmi/amuw)
      yplusoy   = sqrt(yplusoy) 
      utmax     = 0.
      fbl(nmn)    = 0.
c
      do 330 j=nmn+1,nmx
      yplus     = yplusoy*fv(i,j,k,2) 
      yoa       = (1.-fqs)*yplus/aplus  +fqs*1000.
      yoa       = min  (1000.,yoa)
      damp      = 1.-exp(-yoa)
      fbl(j)    = fv(i,j,k,1)*fv(i,j,k,2)*damp
      amixl     = vk*fv(i,j,k,2)*damp 
      eomui(j)  = rey*sgrmi*w(i,j,k,1)*fv(i,j,k,1)*amixl*amixl/
     .            amu(i,j,k)
      utot      = up(i,j,k)**2  +vp(i,j,k)**2  +wp(i,j,k)**2
      utmax     = max  (utmax,utot)
      eomui(j)  = max  (0.,eomui(j))
  330 continue
c
      fmax      = 0.
      ymax      = 0.
      fmax1     = 0.
      nfmax     = nnstrt
      do 335 j=nnstrt,nnstop
      if (fbl(j).gt.fmax) then
         fmax   = fbl(j)
         ymax   = fv(i,j,k,2)
         fmax1  = .9*fmax
         nfmax  = j
c
c     pickup first peak in fbl - due to degani/schiff
c
      else if (ids .eq. 1 .and. fbl(j) .lt. fmax1) then
         go to 331
      endif
  335 continue
c
  331 continue
      if (nfmax.lt.nnstop) then
      dfm     =   fbl(nfmax) - fbl(nfmax-1)
      dfp     =   fbl(nfmax) - fbl(nfmax+1)
      if(fbl(nfmax-1).lt.fbl(nfmax+1)) then
        ymax  = fv(i,nfmax,k,2) + .5*(fv(i,nfmax+1,k,2)-
     .                               fv(i,nfmax,k,2))*(1.-dfp/dfm)
      else
        ymax  = fv(i,nfmax,k,2) - .5*(fv(i,nfmax,k,2)  -
     .                               fv(i,nfmax-1,k,2))*(1.-dfm/dfp)
      end if
      end if
      utmax       = sqrt(utmax) 
      utmin       = sqrt(up(i,nmn+1,k)**2  +vp(i,nmn+1,k)**2 
     .                                     +wp(i,nmn+1,k)**2 )
      udif        = utmax -utmin
      fwake       = cwk*udif*udif*ymax/fmax 
      fwake       = min  (ymax*fmax,fwake)
      do 340 j=nmn+1,nnendt
      fkleb       = 1./(1.+5.5*(ckleb*fv(i,j,k,2)/ymax)**6 )
      eomuo       = ckout*w(i,j,k,1)/amu(i,j,k)*fkleb*fwake 
      eomu(i,j,k) = max  (0.,eomuo)
  340 continue
c
      do 350 j=nmn+1,nnendt
      nedge       = j
      if (eomui(j).ge.eomu(i,j,k)) go to 360
  350 continue
  360 do 370 j=nmn+1,nedge-1
  370 eomu(i,j,k) = eomui(j)
c 
c--------- set ghost cell values of eddy-viscosity
c 
cc    eomu(i,nmn,k)= eomu(i,nmn+1,k)
cc    if (eomu(i,nmn,k).le.0.) eomu(i,nmn,k) = 0.
c
c     add the contribution of this segment to the summation in
c     fv(i,j,k,3) and fv(i,j,k,4) and update the value of
c     eomu including the contribution from the current segment
c
      do 375 j=nmn+1,nnendt
      fv(i,j,k,3) = fv(i,j,k,3) + 1./(fv(i,j,k,2)**2)
      fv(i,j,k,4) = fv(i,j,k,4) + eomu(i,j,k)/(fv(i,j,k,2)**2)
  375 continue
      do 380 j=nmn+1,nnendt
      eomu(i,j,k) = fv(i,j,k,4)/fv(i,j,k,3)
  380 continue
c
  320 continue
c
c--------- end of outer loop on i
c
c
  310 continue
c
c--------- end of outer loop on k
c
c--------- end eddy-viscosity computations for jmn face    -------------
c
      return
c
  400 continue
c--------------------------  j = jmx  face   ---------------------------
c
      nmn    = jmn
      nmx    = jmx
c
c     make sure search/computational range is within the grid
c
      if(nnstrt.lt.(nmn+1).or.nnstrt.gt.nmx) nnstrt = nmn+1
      if(nnstop.gt.nmx.or.nnstop.lt.(nmn+1)) nnstop = nmx
      if(nnendt.gt.nmx.or.nnendt.lt.(nmn+1)) nnendt = nmx
c
c     reverse the search/computational count for jmx face
c
      nnstrt = nmx - nnstrt + nmn + 1
      nnstop = nmx - nnstop + nmn + 1
      nnendt = nmx - nnendt + nmn + 1
c
c     set fqs = 0, if jmx is a solid wall boundary
c     fqs is used for damping of eddy-viscosity near walls
c
      fqs    = 1.0
      if(nbctype.eq.3) fqs = 0.
c
      do 410 k=n1begt+1,n1endt
      do 420 i=n2begt+1,n2endt
c 
      fv(i,nmx+1,k,1) = 0.
      fv(i,nmx+1,k,2) = 0.
c 
      do 422 j=nmn,nmx
      fbl(j)        = fv(i,j,k,1)
      eomu(i,j,k)   = 0.
  422 continue
c
c--------- evaluate normal distance to j=jmx wall
c 
      do 425 j=nmn,nmx
      xn(j)       = .25*(x(i,j,k  ,1) + x(i-1,j,k  ,1)
     .                  +x(i,j,k-1,1) + x(i-1,j,k-1,1))
      yn(j)       = .25*(x(i,j,k  ,2) + x(i-1,j,k  ,2)
     .                  +x(i,j,k-1,2) + x(i-1,j,k-1,2))
      zn(j)       = .25*(x(i,j,k  ,3) + x(i-1,j,k  ,3)
     .                  +x(i,j,k-1,3) + x(i-1,j,k-1,3))
  425 continue
      xnw         = xn(nmx)
      ynw         = yn(nmx)
      znw         = zn(nmx)
      do 426 j=nmx,nmn+1,-1
      xnc         = .5*(xn(j) + xn(j-1)) 
      ync         = .5*(yn(j) + yn(j-1)) 
      znc         = .5*(zn(j) + zn(j-1)) 
      fv(i,j,k,2) = abs((xnc-xnw)*rj2(i,k,1) +(ync-ynw)*rj2(i,k,2)
     .                 +(znc-znw)*rj2(i,k,3))
      fv(i,j,k,2) = max  (fv(i,j,k,2),fv(i,j+1,k,2))
  426 continue
c 
c--------- use linear extrapolation to obtain wall values
c 
      vortw     = 1.5*fv(i,nmx,k,1)  -.5*fv(i,nmx-1,k,1) 
      rhow      = 1.5*w(i,nmx,k,1)   -.5*w(i,nmx-1,k,1) 
      amuw      = 1.5*amu(i,nmx,k)   -.5*amu(i,nmx-1,k) 
      tauw      = amuw*vortw
c
      nvortm  = ismax(nmx-nnstop+1,fbl(nnstop),1) + nnstop -1
      vortm   = fbl(nvortm)
      yplusoy   = abs(rhow*vortm*rey*sgrmi/amuw)
      yplusoy   = sqrt(yplusoy) 
      utmax     = 0.
      fbl(nmx+1)    = 0.
c
      do 430 j=nmn+1,nmx
      yplus     = yplusoy*fv(i,j,k,2) 
      yoa       = (1.-fqs)*yplus/aplus  +fqs*1000.
      yoa       = min  (1000.,yoa)
      damp      = 1.-exp(-yoa)
      fbl(j)    = fv(i,j,k,1)*fv(i,j,k,2)*damp
      amixl     = vk*fv(i,j,k,2)*damp 
      eomui(j)  = rey*sgrmi*w(i,j,k,1)*fv(i,j,k,1)*amixl*amixl/
     .            amu(i,j,k)
      utot      = up(i,j,k)**2  +vp(i,j,k)**2  +wp(i,j,k)**2
      utmax     = max  (utmax,utot)
      eomui(j)  = max  (0.,eomui(j))
  430 continue
c
      fmax      = 0.
      ymax      = 0.
      fmax1     = 0.
      nfmax     = nnstrt
      do 435 j=nnstrt,nnstop,-1
      if (fbl(j).gt.fmax) then
         fmax      = fbl(j)
         ymax      = fv(i,j,k,2)
         fmax1     = .9*fmax
         nfmax     = j
c
c     pickup first peak in fbl - due to degani/schiff
c
      else if (ids .eq. 1 .and. fbl(j) .lt. fmax1) then
         go to 431
      endif
  435 continue
c
  431 continue
      if (nfmax.gt.nnstop) then
      dfm     =   fbl(nfmax) - fbl(nfmax-1)
      dfp     =   fbl(nfmax) - fbl(nfmax+1)
      if(fbl(nfmax-1).lt.fbl(nfmax+1)) then
        ymax  = fv(i,nfmax,k,2) + .5*(fv(i,nfmax+1,k,2)-
     .                               fv(i,nfmax,k,2))*(1.-dfp/dfm)
      else
        ymax  = fv(i,nfmax,k,2) - .5*(fv(i,nfmax,k,2)  -
     .                               fv(i,nfmax-1,k,2))*(1.-dfm/dfp)
      end if
      end if
      utmax       = sqrt(utmax) 
      utmin       = sqrt(up(i,nmx,k)**2  +vp(i,nmx,k)**2 
     .                                   +wp(i,nmx,k)**2 )
      udif        = utmax -utmin
      fwake       = cwk*udif*udif*ymax/fmax 
      fwake       = min  (ymax*fmax,fwake)
      do 440 j=nmx,nnendt,-1
      fkleb       = 1./(1.+5.5*(ckleb*fv(i,j,k,2)/ymax)**6 )
      eomuo       = ckout*w(i,j,k,1)/amu(i,j,k)*fkleb*fwake 
      eomu(i,j,k) = max  (0.,eomuo)
  440 continue
c
      do 450 j=nmx,nnendt,-1
      nedge       = j
      if (eomui(j).ge.eomu(i,j,k)) go to 460
  450 continue
  460 do 470 j=nmx,nedge+1,-1
  470 eomu(i,j,k) = eomui(j)
c 
c--------- set ghost cell values of eddy-viscosity
c 
cc    eomu(i,nmx+1,k)= eomu(i,nmx,k)
cc    if (eomu(i,nmx+1,k).le.0.) eomu(i,nmx+1,k) = 0.
c
c     add the contribution of this segment to the summation in
c     fv(i,j,k,3) and fv(i,j,k,4) and update the value of
c     eomu including the contribution from the current segment
c
      do 475 j=nmx,nnendt,-1
      fv(i,j,k,3) = fv(i,j,k,3) + 1./(fv(i,j,k,2)**2)
      fv(i,j,k,4) = fv(i,j,k,4) + eomu(i,j,k)/(fv(i,j,k,2)**2)
  475 continue
      do 480 j=nmx,nnendt,-1
      eomu(i,j,k) = fv(i,j,k,4)/fv(i,j,k,3)
  480 continue
c
  420 continue
c
c--------- end of outer loop on i
c
  410 continue
c
c--------- end of outer loop on k
c
c
c--------- end eddy-viscosity computations for jmx face    -------------
c 
      return
c
  500 continue
c--------------------------  k = kmn  face   ---------------------------
c
c     2-d flows
c
      if(i2dfl.eq.1) return
c
      nmn    = kmn
      nmx    = kmx
c
c     make sure search/computational range is within the grid
c
      if(nnstrt .lt. (nmn+1) .or. nnstrt .gt. nmx) nnstrt = nmn+1
      if(nnstop .gt. nmx .or. nnstop .lt. (nmn+1)) nnstop = nmx
      if(nnendt .gt. nmx .or. nnendt .lt. (nmn+1)) nnendt = nmx
c
c     set fqs = 0, if kmn is a solid wall boundary
c     fqs is used for damping of eddy-viscosity near walls
c
      fqs    = 1.0
      if(nbctype.eq.3) fqs = 0.
c
      do 510 i=n1begt+1,n1endt
      do 520 j=n2begt+1,n2endt
c 
      fv(i,j,nmn,1) = 0.
      fv(i,j,nmn,2) = 0.
c 
      do 522 k=nmn,nmx
      fbl(k)        = fv(i,j,k,1)
      eomu(i,j,k)   = 0.
  522 continue
c
c--------- evaluate normal distance to k=kmn wall
c 
      do 525 k=nmn,nmx
      xn(k)       = .25*(x(i,j,k  ,1) + x(i-1,j,k  ,1)
     .                  +x(i,j-1,k,1) + x(i-1,j-1,k,1))
      yn(k)       = .25*(x(i,j,k  ,2) + x(i-1,j,k  ,2)
     .                  +x(i,j-1,k,2) + x(i-1,j-1,k,2))
      zn(k)       = .25*(x(i,j,k  ,3) + x(i-1,j,k  ,3)
     .                  +x(i,j-1,k,3) + x(i-1,j-1,k,3))
  525 continue
      xnw         = xn(nmn)
      ynw         = yn(nmn)
      znw         = zn(nmn)
      do 526 k=nmn+1,nmx
      xnc         = .5*(xn(k) + xn(k-1)) 
      ync         = .5*(yn(k) + yn(k-1)) 
      znc         = .5*(zn(k) + zn(k-1)) 
      fv(i,j,k,2) = abs((xnc-xnw)*rk1(i,j,1) +(ync-ynw)*rk1(i,j,2)
     .                 +(znc-znw)*rk1(i,j,3))
      fv(i,j,k,2) = max  (fv(i,j,k,2),fv(i,j,k-1,2))
  526 continue
c 
c--------- use linear extrapolation to obtain wall values
c 
      vortw     = 1.5*fv(i,j,nmn+1,1)  -.5*fv(i,j,nmn+2,1) 
      rhow      = 1.5*w(i,j,nmn+1,1)   -.5*w(i,j,nmn+2,1) 
      amuw      = 1.5*amu(i,j,nmn+1)   -.5*amu(i,j,nmn+2) 
      tauw      = amuw*vortw
c
      nvortm  = ismax(nnstop-2,fbl(nmn+1),1) + nmn
      vortm   = fbl(nvortm)
      yplusoy   = abs(rhow*vortm*rey*sgrmi/amuw)
      yplusoy   = sqrt(yplusoy) 
      utmax     = 0.
      fbl(nmn)    = 0.
c
      do 530 k=nmn+1,nmx
      yplus     = yplusoy*fv(i,j,k,2) 
      yoa       = (1.-fqs)*yplus/aplus  +fqs*1000.
      yoa       = min  (1000.,yoa)
      damp      = 1.-exp(-yoa)
      fbl(k)    = fv(i,j,k,1)*fv(i,j,k,2)*damp
      amixl     = vk*fv(i,j,k,2)*damp 
      eomui(k)  = rey*sgrmi*w(i,j,k,1)*fv(i,j,k,1)*amixl*amixl/
     .            amu(i,j,k)
      utot      = up(i,j,k)**2  +vp(i,j,k)**2  +wp(i,j,k)**2
      utmax     = max  (utmax,utot)
      eomui(k)  = max  (0.,eomui(k))
  530 continue
c
      fmax      = 0.
      ymax      = 0.
      fmax1     = 0.
      nfmax     = nnstrt
      do 535 k=nnstrt,nnstop
      if (fbl(k).gt.fmax) then
         fmax      = fbl(k)
         ymax      = fv(i,j,k,2)
         fmax1     = .9*fmax
         nfmax     = k
c
c     pickup first peak in fbl - due to degani/schiff
c
      else if (ids .eq. 1 .and. fbl(k) .lt. fmax1) then
         go to 531
      endif
  535 continue
c
  531 continue
      if (nfmax.lt.nnstop) then
      dfm     =   fbl(nfmax) - fbl(nfmax-1)
      dfp     =   fbl(nfmax) - fbl(nfmax+1)
      if(fbl(nfmax-1).lt.fbl(nfmax+1)) then
        ymax  = fv(i,j,nfmax,2) + .5*(fv(i,j,nfmax+1,2)-
     .                               fv(i,j,nfmax,2))*(1.-dfp/dfm)
      else
        ymax  = fv(i,j,nfmax,2) - .5*(fv(i,j,nfmax,2)  -
     .                               fv(i,j,nfmax-1,2))*(1.-dfm/dfp)
      end if
      end if
      utmax       = sqrt(utmax) 
      utmin       = sqrt(up(i,j,nmn+1)**2  +vp(i,j,nmn+1)**2 
     .                                     +wp(i,j,nmn+1)**2 )
      udif        = utmax -utmin
      fwake       = cwk*udif*udif*ymax/fmax 
      fwake       = min  (ymax*fmax,fwake)
      do 540 k=nmn+1,nnendt
      fkleb       = 1./(1.+5.5*(ckleb*fv(i,j,k,2)/ymax)**6 )
      eomuo       = ckout*w(i,j,k,1)/amu(i,j,k)*fkleb*fwake 
      eomu(i,j,k) = max  (0.,eomuo)
  540 continue
c
      do 550 k=nmn+1,nnendt
      nedge       = k
      if (eomui(k).ge.eomu(i,j,k)) go to 560
  550 continue
  560 do 570 k=nmn+1,nedge-1
  570 eomu(i,j,k) = eomui(k)
c 
c--------- set ghost cell values of eddy-viscosity
c 
cc    eomu(i,j,nmn)= eomu(i,j,nmn+1)
cc    if (eomu(i,j,nmn).le.0.) eomu(i,j,nmn) = 0.
c
c     add the contribution of this segment to the summation in
c     fv(i,j,k,3) and fv(i,j,k,4) and update the value of
c     eomu including the contribution from the current segment
c
      do 575 k=nmn+1,nnendt
      fv(i,j,k,3) = fv(i,j,k,3) + 1./(fv(i,j,k,2)**2)
      fv(i,j,k,4) = fv(i,j,k,4) + eomu(i,j,k)/(fv(i,j,k,2)**2)
  575 continue
      do 580 k=nmn+1,nnendt
      eomu(i,j,k) = fv(i,j,k,4)/fv(i,j,k,3)
  580 continue
c
  520 continue
c
c--------- end of outer loop on j
c
c
  510 continue
c
c--------- end of outer loop on i
c
c
c--------- end eddy-viscosity computations for kmn face    -------------
c
      return
c
  600 continue
c--------------------------  k = kmx  face   ---------------------------
c
c     2-d flows
c
      if(i2dfl.eq.1) return
c
      nmn    = kmn
      nmx    = kmx
c
c     make sure search/computational range is within the grid
c
      if(nnstrt.lt.(nmn+1).or.nnstrt.gt.nmx) nnstrt = nmn+1
      if(nnstop.gt.nmx.or.nnstop.lt.(nmn+1)) nnstop = nmx
      if(nnendt.gt.nmx.or.nnendt.lt.(nmn+1)) nnendt = nmx
c
c     reverse the search/computational count for kmx face
c
      nnstrt = nmx - nnstrt + nmn + 1
      nnstop = nmx - nnstop + nmn + 1
      nnendt = nmx - nnendt + nmn + 1
c
c     set fqs = 0, if kmx is a solid wall boundary
c     fqs is used for damping of eddy-viscosity near walls
c
      fqs    = 1.0
      if(nbctype.eq.3) fqs = 0.
c
      do 610 i=n1begt+1,n1endt
      do 620 j=n2begt+1,n2endt
c 
      fv(i,j,nmx+1,1) = 0.
      fv(i,j,nmx+1,2) = 0.
c 
      do 622 k=nmn,nmx
      fbl(k)        = fv(i,j,k,1)
      eomu(i,j,k)   = 0.
  622 continue
c
c--------- evaluate normal distance to k=kmx wall
c 
      do 625 k=nmn,nmx
      xn(k)       = .25*(x(i,j,k  ,1) + x(i-1,j,k  ,1)
     .                  +x(i,j-1,k,1) + x(i-1,j-1,k,1))
      yn(k)       = .25*(x(i,j,k  ,2) + x(i-1,j,k  ,2)
     .                  +x(i,j-1,k,2) + x(i-1,j-1,k,2))
      zn(k)       = .25*(x(i,j,k  ,3) + x(i-1,j,k  ,3)
     .                  +x(i,j-1,k,3) + x(i-1,j-1,k,3))
  625 continue
      xnw         = xn(nmx)
      ynw         = yn(nmx)
      znw         = zn(nmx)
      do 626 k=nmx,nmn+1,-1
      xnc         = .5*(xn(k) + xn(k-1)) 
      ync         = .5*(yn(k) + yn(k-1)) 
      znc         = .5*(zn(k) + zn(k-1)) 
      fv(i,j,k,2) = abs((xnc-xnw)*rk2(i,j,1) +(ync-ynw)*rk2(i,j,2)
     .                 +(znc-znw)*rk2(i,j,3))
      fv(i,j,k,2) = max  (fv(i,j,k,2),fv(i,j,k+1,2))
  626 continue
c 
c--------- use linear extrapolation to obtain wall values
c 
      vortw     = 1.5*fv(i,j,nmx,1)  -.5*fv(i,j,nmx-1,1) 
      rhow      = 1.5*w(i,j,nmx,1)   -.5*w(i,j,nmx-1,1) 
      amuw      = 1.5*amu(i,j,nmx)   -.5*amu(i,j,nmx-1) 
      tauw      = amuw*vortw
c
      nvortm  = ismax(nmx-nnstop+1,fbl(nnstop),1) + nnstop -1
      vortm   = fbl(nvortm)
      yplusoy   = abs(rhow*vortm*rey*sgrmi/amuw)
      yplusoy   = sqrt(yplusoy) 
      utmax     = 0.
      fbl(nmx+1)    = 0.
c
      do 630 k=nmn+1,nmx
      yplus     = yplusoy*fv(i,j,k,2) 
      yoa       = (1.-fqs)*yplus/aplus  +fqs*1000.
      yoa       = min  (1000.,yoa)
      damp      = 1.-exp(-yoa)
      fbl(k)    = fv(i,j,k,1)*fv(i,j,k,2)*damp
      amixl     = vk*fv(i,j,k,2)*damp 
      eomui(k)  = rey*sgrmi*w(i,j,k,1)*fv(i,j,k,1)*amixl*amixl/
     .            amu(i,j,k)
      utot      = up(i,j,k)**2  +vp(i,j,k)**2  +wp(i,j,k)**2
      utmax     = max  (utmax,utot)
      eomui(k)  = max  (0.,eomui(k))
  630 continue
c
      fmax      = 0.
      ymax      = 0.
      fmax1     = 0.
      nfmax     = nnstrt
      do 635 k=nnstrt,nnstop,-1
      if (fbl(k).gt.fmax) then
         fmax      = fbl(k)
         ymax      = fv(i,j,k,2)
         fmax1     = .9*fmax
         nfmax     = k
c
c     pickup first peak in fbl - due to degani/schiff
c
      else if (ids .eq. 1 .and. fbl(k) .lt. fmax1) then
         go to 631
      endif
  635 continue
c
  631 continue
      if (nfmax.gt.nnstop) then
      dfm     =   fbl(nfmax) - fbl(nfmax-1)
      dfp     =   fbl(nfmax) - fbl(nfmax+1)
      if(fbl(nfmax-1).lt.fbl(nfmax+1)) then
        ymax  = fv(i,j,nfmax,2) + .5*(fv(i,j,nfmax+1,2)-
     .                               fv(i,j,nfmax,2))*(1.-dfp/dfm)
      else
        ymax  = fv(i,j,nfmax,2) - .5*(fv(i,j,nfmax,2)  -
     .                               fv(i,j,nfmax-1,2))*(1.-dfm/dfp)
      end if
      end if
      utmax       = sqrt(utmax) 
      utmin       = sqrt(up(i,j,nmx)**2  +vp(i,j,nmx)**2 
     .                                   +wp(i,j,nmx)**2 )
      udif        = utmax -utmin
      fwake       = cwk*udif*udif*ymax/fmax 
      fwake       = min  (ymax*fmax,fwake)
      do 640 k=nmx,nnendt,-1
      fkleb       = 1./(1.+5.5*(ckleb*fv(i,j,k,2)/ymax)**6 )
      eomuo       = ckout*w(i,j,k,1)/amu(i,j,k)*fkleb*fwake 
      eomu(i,j,k) = max  (0.,eomuo)
  640 continue
c
      do 650 k=nmx,nnendt,-1
      nedge       = k
      if (eomui(k).ge.eomu(i,j,k)) go to 660
  650 continue
  660 do 670 k=nmx,nedge+1,-1
  670 eomu(i,j,k) = eomui(k)
c 
c--------- set ghost cell values of eddy-viscosity
c 
cc    eomu(i,j,nmx+1)= eomu(i,j,nmx)
cc    if (eomu(i,j,nmx+1).le.0.) eomu(i,j,nmx+1) = 0.
c
c     add the contribution of this segment to the summation in
c     fv(i,j,k,3) and fv(i,j,k,4) and update the value of
c     eomu including the contribution from the current segment
c
      do 675 k=nmx,nnendt,-1
      fv(i,j,k,3) = fv(i,j,k,3) + 1./(fv(i,j,k,2)**2)
      fv(i,j,k,4) = fv(i,j,k,4) + eomu(i,j,k)/(fv(i,j,k,2)**2)
  675 continue
      do 680 k=nmx,nnendt,-1
      eomu(i,j,k) = fv(i,j,k,4)/fv(i,j,k,3)
  680 continue
c
  620 continue
c
c--------- end of outer loop on j
c
  610 continue
c
c--------- end of outer loop on i
c
c
c--------- end eddy-viscosity computations for kmx face    -------------
c
      return
      end 
