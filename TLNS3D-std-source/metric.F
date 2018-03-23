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
      subroutine  metric (imn,jmn,kmn,
     .                    imx,jmx,kmx,imp1,jmp1,kmp1,x,
     .                    six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                    ri1,ri2,rj1,rj2,rk1,rk2,ibloc,igrid,nodes)
c=======================================================================
c
c             calculates metrics and surface normals
c
c             six,siy .. skz  etc. are actual surface areas
c             and contain no multiplicative factors
c 

      dimension
     .  six(imp1,jmp1,kmp1), siy(imp1,jmp1,kmp1), siz(imp1,jmp1,kmp1),
     .  sjx(imp1,jmp1,kmp1), sjy(imp1,jmp1,kmp1), sjz(imp1,jmp1,kmp1),
     .  skx(imp1,jmp1,kmp1), sky(imp1,jmp1,kmp1), skz(imp1,jmp1,kmp1),
     .  ri1(jmp1,kmp1,3)   , ri2(jmp1,kmp1,3)   , rj1(imp1,kmp1,3)   ,
     .  rj2(imp1,kmp1,3)   , rk1(imp1,jmp1,3)   , rk2(imp1,jmp1,3)   ,
     .   x(imp1,jmp1,kmp1,3)
c
      dimension    nodes(1)
c
      data init/0/
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      tol  = 1.e-17
cvn   tolsqt = sqrt(tol)
      tol  = -tol
      sf   = 0.01
      ineg = 0
      ineg = 0
c     write(6,900) ibloc,igrid,imx,jmx,kmx,init
c900  format(2x,"ibloc,igrid,imx,jmx,kmx,init",6i5)
      init = init + 1
c
c--------- metric terms in i-direction
c
      do 10 k=kmn+1,kmx
      do 10 j=jmn+1,jmx
      do 11 i=imn,imx
      six(i,j,k) = (x(i,j,k-1,2) - x(i,j-1,k  ,2)) 
     .            *(x(i,j,k  ,3) - x(i,j-1,k-1,3)) 
     .            -(x(i,j,k-1,3) - x(i,j-1,k  ,3)) 
     .            *(x(i,j,k  ,2) - x(i,j-1,k-1,2)) 
      siy(i,j,k) = (x(i,j,k-1,3) - x(i,j-1,k  ,3)) 
     .            *(x(i,j,k  ,1) - x(i,j-1,k-1,1)) 
     .            -(x(i,j,k-1,1) - x(i,j-1,k  ,1)) 
     .            *(x(i,j,k  ,3) - x(i,j-1,k-1,3)) 
      siz(i,j,k) = (x(i,j,k-1,1) - x(i,j-1,k  ,1)) 
     .            *(x(i,j,k  ,2) - x(i,j-1,k-1,2)) 
     .            -(x(i,j,k-1,2) - x(i,j-1,k  ,2)) 
     .            *(x(i,j,k  ,1) - x(i,j-1,k-1,1)) 
c
      six(i,j,k) = .5*six(i,j,k)
      siy(i,j,k) = .5*siy(i,j,k)
      siz(i,j,k) = .5*siz(i,j,k)
c
   11 continue
      do 12 i=imn,imx
      stot = six(i,j,k)**2 + siy(i,j,k)**2 + siz(i,j,k)**2
      if(stot.le.tol.and.i.lt.imx) then
         six(i,j,k) = six(i+1,j,k)*sf
         siy(i,j,k) = siy(i+1,j,k)*sf
         siz(i,j,k) = siz(i+1,j,k)*sf
       stotn = six(i,j,k)**2 + siy(i,j,k)**2 + siz(i,j,k)**2
       if(stotn.le.tol*sf*sf) then
       ineg = ineg + 1
       write (6,901) i,j,k,init,stotn
c      if(ineg.ge.20) stop 901
       if(ineg.ge.20) call ERREXIT (nodes)
 901   format(2x,"i,j,k,stoti",4i5,2x,e15.6)
       endif

      endif
      if(stot.le.tol.and.i.eq.imx) then
         six(i,j,k) = six(i-1,j,k)*sf
         siy(i,j,k) = siy(i-1,j,k)*sf
         siz(i,j,k) = siz(i-1,j,k)*sf
       if(stotn.le.tol*sf*sf) then
       ineg = ineg + 1
       write (6,901) i,j,k,init,stotn
       if(ineg.ge.20) call ERREXIT (nodes)
       endif
      endif
   12 continue
c
   10 continue
c
c--------- metric terms in j-direction
c
      do 20 k=kmn+1,kmx
      do 20 i=imn+1,imx
      do 21 j=jmn,jmx
      sjx(i,j,k) = (x(i  ,j,k-1,2) - x(i-1,j,k,2)) 
     .            *(x(i-1,j,k-1,3) - x(i  ,j,k,3)) 
     .            -(x(i  ,j,k-1,3) - x(i-1,j,k,3)) 
     .            *(x(i-1,j,k-1,2) - x(i  ,j,k,2)) 
      sjy(i,j,k) = (x(i  ,j,k-1,3) - x(i-1,j,k,3)) 
     .            *(x(i-1,j,k-1,1) - x(i  ,j,k,1)) 
     .            -(x(i  ,j,k-1,1) - x(i-1,j,k,1)) 
     .            *(x(i-1,j,k-1,3) - x(i  ,j,k,3)) 
      sjz(i,j,k) = (x(i  ,j,k-1,1) - x(i-1,j,k,1)) 
     .            *(x(i-1,j,k-1,2) - x(i  ,j,k,2)) 
     .            -(x(i  ,j,k-1,2) - x(i-1,j,k,2)) 
     .            *(x(i-1,j,k-1,1) - x(i  ,j,k,1)) 
c
      sjx(i,j,k) = .5*sjx(i,j,k)
      sjy(i,j,k) = .5*sjy(i,j,k)
      sjz(i,j,k) = .5*sjz(i,j,k)
c
   21 continue
      do 22 j=jmn,jmx
      stot = sjx(i,j,k)**2 + sjy(i,j,k)**2 + sjz(i,j,k)**2
      if(stot.le.tol.and.j.lt.jmx) then
         sjx(i,j,k) = sjx(i,j+1,k)*sf
         sjy(i,j,k) = sjy(i,j+1,k)*sf
         sjz(i,j,k) = sjz(i,j+1,k)*sf
       stotn = sjx(i,j,k)**2 + sjy(i,j,k)**2 + sjz(i,j,k)**2
       if(stotn.le.tol*sf*sf) then
       ineg = ineg + 1
       write (6,902) i,j,k,stotn
       if(ineg.ge.20) call ERREXIT (nodes)
 902   format(2x,"i,j,k,stotj",3i5,2x,e15.6)
       endif
      endif
      if(stot.le.tol.and.j.eq.jmx) then
         sjx(i,j,k) = sjx(i,j-1,k)*sf
         sjy(i,j,k) = sjy(i,j-1,k)*sf
         sjz(i,j,k) = sjz(i,j-1,k)*sf
       stotn = sjx(i,j,k)**2 + sjy(i,j,k)**2 + sjz(i,j,k)**2
       if(stotn.le.tol*sf*sf) then
       ineg = ineg + 1
       write (6,902) i,j,k,stotn
       if(ineg.ge.20) call ERREXIT (nodes)
       endif
      endif
   22 continue
c
   20 continue
c
c--------- metric terms in k-direction
c
      do 30 j=jmn+1,jmx
      do 30 i=imn+1,imx
      do 31 k=kmn,kmx
      skx(i,j,k) = (x(i  ,j,k,2) - x(i-1,j-1,k,2)) 
     .            *(x(i-1,j,k,3) - x(i  ,j-1,k,3)) 
     .            -(x(i  ,j,k,3) - x(i-1,j-1,k,3)) 
     .            *(x(i-1,j,k,2) - x(i  ,j-1,k,2)) 
      sky(i,j,k) = (x(i  ,j,k,3) - x(i-1,j-1,k,3)) 
     .            *(x(i-1,j,k,1) - x(i  ,j-1,k,1)) 
     .            -(x(i  ,j,k,1) - x(i-1,j-1,k,1)) 
     .            *(x(i-1,j,k,3) - x(i  ,j-1,k,3)) 
      skz(i,j,k) = (x(i  ,j,k,1) - x(i-1,j-1,k,1)) 
     .            *(x(i-1,j,k,2) - x(i  ,j-1,k,2)) 
     .            -(x(i  ,j,k,2) - x(i-1,j-1,k,2)) 
     .            *(x(i-1,j,k,1) - x(i  ,j-1,k,1)) 
c
      skx(i,j,k) = .5*skx(i,j,k)
      sky(i,j,k) = .5*sky(i,j,k)
      skz(i,j,k) = .5*skz(i,j,k)
c
   31 continue
      do 32 k=kmn,kmx
      stot = skx(i,j,k)**2 + sky(i,j,k)**2 + skz(i,j,k)**2
      if(stot.le.tol.and.k.lt.kmx) then
         skx(i,j,k) = skx(i,j,k+1)*sf
         sky(i,j,k) = sky(i,j,k+1)*sf
         skz(i,j,k) = skz(i,j,k+1)*sf
       stotn = skx(i,j,k)**2 + sky(i,j,k)**2 + skz(i,j,k)**2
       if(stotn.le.tol*sf*sf) then
       ineg = ineg + 1
       write (6,903) i,j,k,stotn
       if(ineg.ge.20) call ERREXIT (nodes)
 903   format(2x,"i,j,k,stotk",3i5,2x,e15.6)
       endif
      endif
      if(stot.le.tol.and.k.eq.kmx) then
         skx(i,j,k) = skx(i,j,k-1)*sf
         sky(i,j,k) = sky(i,j,k-1)*sf
         skz(i,j,k) = skz(i,j,k-1)*sf
       stotn = skx(i,j,k)**2 + sky(i,j,k)**2 + skz(i,j,k)**2
       if(stotn.le.tol*sf*sf) then
       ineg = ineg + 1
       write (6,903) i,j,k,stotn
       if(ineg.ge.20) call ERREXIT (nodes)
       endif
      endif
   32 continue
c
   30 continue
c
c--------- surface normals at i=imn and i=imx
c
      do 40 k=kmn+1,kmx
      do 40 j=jmn+1,jmx
      if(six(imn,j,k).eq.0.0 .and. siy(imn,j,k).eq.0.0 .and.
     &                             siz(imn,j,k).eq.0.0) then
         aden = 1.0
      else
         aden = sqrt(six(imn,j,k)**2+siy(imn,j,k)**2+siz(imn,j,k)**2)
      end if
      a     =  1./aden
      ri1(j,k,1)  = a*six(imn,j,k)
      ri1(j,k,2)  = a*siy(imn,j,k)
      ri1(j,k,3)  = a*siz(imn,j,k)
c
      if(six(imx,j,k).eq.0.0 .and. siy(imx,j,k).eq.0.0 .and.
     &                             siz(imx,j,k).eq.0.0) then
         aden = 1.0
      else
         aden = sqrt(six(imx,j,k)**2+siy(imx,j,k)**2+siz(imx,j,k)**2)
      end if
      a     =  1./aden
      ri2(j,k,1)  = a*six(imx,j,k)
      ri2(j,k,2)  = a*siy(imx,j,k)
      ri2(j,k,3)  = a*siz(imx,j,k)
   40 continue
c
c------- fill-in ghost cell values using linear extrapolation
c
      do 42 k=kmn+1,kmx
      if(jmx.gt.jmn+1) then
c
      ri1(jmn  ,k,1) = 2.*ri1(jmn+1,k,1) - ri1(jmn+2,k,1)
      ri1(jmn  ,k,2) = 2.*ri1(jmn+1,k,2) - ri1(jmn+2,k,2)
      ri1(jmn  ,k,3) = 2.*ri1(jmn+1,k,3) - ri1(jmn+2,k,3)
      ri2(jmn  ,k,1) = 2.*ri2(jmn+1,k,1) - ri2(jmn+2,k,1)
      ri2(jmn  ,k,2) = 2.*ri2(jmn+1,k,2) - ri2(jmn+2,k,2)
      ri2(jmn  ,k,3) = 2.*ri2(jmn+1,k,3) - ri2(jmn+2,k,3)
c
      ri1(jmx+1,k,1) = 2.*ri1(jmx  ,k,1) - ri1(jmx-1,k,1)
      ri1(jmx+1,k,2) = 2.*ri1(jmx  ,k,2) - ri1(jmx-1,k,2)
      ri1(jmx+1,k,3) = 2.*ri1(jmx  ,k,3) - ri1(jmx-1,k,3)
      ri2(jmx+1,k,1) = 2.*ri2(jmx  ,k,1) - ri2(jmx-1,k,1)
      ri2(jmx+1,k,2) = 2.*ri2(jmx  ,k,2) - ri2(jmx-1,k,2)
      ri2(jmx+1,k,3) = 2.*ri2(jmx  ,k,3) - ri2(jmx-1,k,3)
c
      else
c
      ri1(jmx+1,k,1) = ri1(jmx  ,k,1)
      ri1(jmx+1,k,2) = ri1(jmx  ,k,2)
      ri1(jmx+1,k,3) = ri1(jmx  ,k,3)
      ri2(jmx+1,k,1) = ri2(jmx  ,k,1)
      ri2(jmx+1,k,2) = ri2(jmx  ,k,2)
      ri2(jmx+1,k,3) = ri2(jmx  ,k,3)
c
      ri1(jmn  ,k,1) = ri1(jmn+1,k,1)
      ri1(jmn  ,k,2) = ri1(jmn+1,k,2)
      ri1(jmn  ,k,3) = ri1(jmn+1,k,3)
      ri2(jmn  ,k,1) = ri2(jmn+1,k,1)
      ri2(jmn  ,k,2) = ri2(jmn+1,k,2)
      ri2(jmn  ,k,3) = ri2(jmn+1,k,3)
c
      endif
   42 continue
c 
         if(i2dfl.ne.1.and.kmx.gt.kmn+1) then
      do 44 j=jmn,jmx+1
c
      ri1(j,kmn  ,1) = 2.*ri1(j,kmn+1,1) - ri1(j,kmn+2,1)
      ri1(j,kmn  ,2) = 2.*ri1(j,kmn+1,2) - ri1(j,kmn+2,2)
      ri1(j,kmn  ,3) = 2.*ri1(j,kmn+1,3) - ri1(j,kmn+2,3)
      ri2(j,kmn  ,1) = 2.*ri2(j,kmn+1,1) - ri2(j,kmn+2,1)
      ri2(j,kmn  ,2) = 2.*ri2(j,kmn+1,2) - ri2(j,kmn+2,2)
      ri2(j,kmn  ,3) = 2.*ri2(j,kmn+1,3) - ri2(j,kmn+2,3)
c      
      ri1(j,kmx+1,1) = 2.*ri1(j,kmx  ,1) - ri1(j,kmx-1,1)
      ri1(j,kmx+1,2) = 2.*ri1(j,kmx  ,2) - ri1(j,kmx-1,2)
      ri1(j,kmx+1,3) = 2.*ri1(j,kmx  ,3) - ri1(j,kmx-1,3)
      ri2(j,kmx+1,1) = 2.*ri2(j,kmx  ,1) - ri2(j,kmx-1,1)
      ri2(j,kmx+1,2) = 2.*ri2(j,kmx  ,2) - ri2(j,kmx-1,2)
      ri2(j,kmx+1,3) = 2.*ri2(j,kmx  ,3) - ri2(j,kmx-1,3)
  44  continue
c
         else
      do 45 j=jmn,jmx+1
c
      ri1(j,kmn  ,1) = ri1(j,kmn+1,1)
      ri1(j,kmn  ,2) = ri1(j,kmn+1,2)
      ri1(j,kmn  ,3) = ri1(j,kmn+1,3)
      ri2(j,kmn  ,1) = ri2(j,kmn+1,1)
      ri2(j,kmn  ,2) = ri2(j,kmn+1,2)
      ri2(j,kmn  ,3) = ri2(j,kmn+1,3)
c      
      ri1(j,kmx+1,1) = ri1(j,kmx  ,1) 
      ri1(j,kmx+1,2) = ri1(j,kmx  ,2) 
      ri1(j,kmx+1,3) = ri1(j,kmx  ,3) 
      ri2(j,kmx+1,1) = ri2(j,kmx  ,1) 
      ri2(j,kmx+1,2) = ri2(j,kmx  ,2) 
      ri2(j,kmx+1,3) = ri2(j,kmx  ,3) 
   45 continue
         endif
c
c--------- surface normals at j=jmn and j=jmx
c
      do 50 k=kmn+1,kmx
      do 50 i=imn+1,imx
      if(sjx(i,jmn,k).eq.0.0 .and. sjy(i,jmn,k).eq.0.0 .and.
     &                             sjz(i,jmn,k).eq.0.0) then
         aden = 1.0
      else
         aden = sqrt(sjx(i,jmn,k)**2+sjy(i,jmn,k)**2+sjz(i,jmn,k)**2)
      end if
      a     =  1./aden
      rj1(i,k,1)  = a*sjx(i,jmn,k)
      rj1(i,k,2)  = a*sjy(i,jmn,k)
      rj1(i,k,3)  = a*sjz(i,jmn,k)
c
      if(sjx(i,jmx,k).eq.0.0 .and. sjy(i,jmx,k).eq.0.0 .and.
     &                             sjz(i,jmx,k).eq.0.0) then
         aden = 1.0
      else
         aden = sqrt(sjx(i,jmx,k)**2+sjy(i,jmx,k)**2+sjz(i,jmx,k)**2)
      end if
      a     =  1./aden
      rj2(i,k,1)  = a*sjx(i,jmx,k)
      rj2(i,k,2)  = a*sjy(i,jmx,k)
      rj2(i,k,3)  = a*sjz(i,jmx,k)
   50 continue
c
c------- fill-in ghost cell values using linear extrapolation
c
      do 52 k=kmn+1,kmx
      if (imx.gt.imn+1) then
c
      rj1(imn  ,k,1) = 2.*rj1(imn+1,k,1) - rj1(imn+2,k,1)
      rj1(imn  ,k,2) = 2.*rj1(imn+1,k,2) - rj1(imn+2,k,2)
      rj1(imn  ,k,3) = 2.*rj1(imn+1,k,3) - rj1(imn+2,k,3)
      rj2(imn  ,k,1) = 2.*rj2(imn+1,k,1) - rj2(imn+2,k,1)
      rj2(imn  ,k,2) = 2.*rj2(imn+1,k,2) - rj2(imn+2,k,2)
      rj2(imn  ,k,3) = 2.*rj2(imn+1,k,3) - rj2(imn+2,k,3)
c
      rj1(imx+1,k,1) = 2.*rj1(imx  ,k,1) - rj1(imx-1,k,1)
      rj1(imx+1,k,2) = 2.*rj1(imx  ,k,2) - rj1(imx-1,k,2)
      rj1(imx+1,k,3) = 2.*rj1(imx  ,k,3) - rj1(imx-1,k,3)
      rj2(imx+1,k,1) = 2.*rj2(imx  ,k,1) - rj2(imx-1,k,1)
      rj2(imx+1,k,2) = 2.*rj2(imx  ,k,2) - rj2(imx-1,k,2)
      rj2(imx+1,k,3) = 2.*rj2(imx  ,k,3) - rj2(imx-1,k,3)
c
      else
c
      rj1(imx+1,k,1) = rj1(imx  ,k,1)
      rj1(imx+1,k,2) = rj1(imx  ,k,2)
      rj1(imx+1,k,3) = rj1(imx  ,k,3)
      rj2(imx+1,k,1) = rj2(imx  ,k,1)
      rj2(imx+1,k,2) = rj2(imx  ,k,2)
      rj2(imx+1,k,3) = rj2(imx  ,k,3)
c
      rj1(imn  ,k,1) = rj1(imn+1,k,1)
      rj1(imn  ,k,2) = rj1(imn+1,k,2)
      rj1(imn  ,k,3) = rj1(imn+1,k,3)
      rj2(imn  ,k,1) = rj2(imn+1,k,1)
      rj2(imn  ,k,2) = rj2(imn+1,k,2)
      rj2(imn  ,k,3) = rj2(imn+1,k,3)
c
      endif
c
   52 continue
c 
         if(i2dfl.ne.1.and.kmx.gt.kmn+1) then
      do 54 i=imn,imx+1
c
      rj1(i,kmn  ,1) = 2.*rj1(i,kmn+1,1) - rj1(i,kmn+2,1)
      rj1(i,kmn  ,2) = 2.*rj1(i,kmn+1,2) - rj1(i,kmn+2,2)
      rj1(i,kmn  ,3) = 2.*rj1(i,kmn+1,3) - rj1(i,kmn+2,3)
      rj2(i,kmn  ,1) = 2.*rj2(i,kmn+1,1) - rj2(i,kmn+2,1)
      rj2(i,kmn  ,2) = 2.*rj2(i,kmn+1,2) - rj2(i,kmn+2,2)
      rj2(i,kmn  ,3) = 2.*rj2(i,kmn+1,3) - rj2(i,kmn+2,3)
c
      rj1(i,kmx+1,1) = 2.*rj1(i,kmx  ,1) - rj1(i,kmx-1,1)
      rj1(i,kmx+1,2) = 2.*rj1(i,kmx  ,2) - rj1(i,kmx-1,2)
      rj1(i,kmx+1,3) = 2.*rj1(i,kmx  ,3) - rj1(i,kmx-1,3)
      rj2(i,kmx+1,1) = 2.*rj2(i,kmx  ,1) - rj2(i,kmx-1,1)
      rj2(i,kmx+1,2) = 2.*rj2(i,kmx  ,2) - rj2(i,kmx-1,2)
      rj2(i,kmx+1,3) = 2.*rj2(i,kmx  ,3) - rj2(i,kmx-1,3)
   54 continue
c
         else
      do 55 i=imn,imx+1
c
      rj1(i,kmn  ,1) = rj1(i,kmn+1,1) 
      rj1(i,kmn  ,2) = rj1(i,kmn+1,2) 
      rj1(i,kmn  ,3) = rj1(i,kmn+1,3)
      rj2(i,kmn  ,1) = rj2(i,kmn+1,1)
      rj2(i,kmn  ,2) = rj2(i,kmn+1,2)
      rj2(i,kmn  ,3) = rj2(i,kmn+1,3)
c
      rj1(i,kmx+1,1) = 2.*rj1(i,kmx  ,1) 
      rj1(i,kmx+1,2) = 2.*rj1(i,kmx  ,2) 
      rj1(i,kmx+1,3) = 2.*rj1(i,kmx  ,3) 
      rj2(i,kmx+1,1) = 2.*rj2(i,kmx  ,1) 
      rj2(i,kmx+1,2) = 2.*rj2(i,kmx  ,2) 
      rj2(i,kmx+1,3) = 2.*rj2(i,kmx  ,3) 
   55 continue
         endif
c
c
c--------- surface normal at k =kmn and k=kmx
c
      do 60 j=jmn+1,jmx
      do 60 i=imn+1,imx
      if(skx(i,j,kmn).eq.0.0 .and. sky(i,j,kmn).eq.0.0 .and.
     &                             skz(i,j,kmn).eq.0.0) then
         aden = 1.0
      else
         aden = sqrt(skx(i,j,kmn)**2+sky(i,j,kmn)**2+skz(i,j,kmn)**2)
      end if
      a     =  1./aden
      rk1(i,j,1)  = a*skx(i,j,kmn)
      rk1(i,j,2)  = a*sky(i,j,kmn)
      rk1(i,j,3)  = a*skz(i,j,kmn)
c
      if(skx(i,j,kmx).eq.0.0 .and. sky(i,j,kmx).eq.0.0 .and.
     &                             skz(i,j,kmx).eq.0.0) then
         aden = 1.0
      else
         aden = sqrt(skx(i,j,kmx)**2+sky(i,j,kmx)**2+skz(i,j,kmx)**2)
      end if
      a     =  1./aden
      rk2(i,j,1)  = a*skx(i,j,kmx)
      rk2(i,j,2)  = a*sky(i,j,kmx)
      rk2(i,j,3)  = a*skz(i,j,kmx)
   60 continue
c
c------- fill-in ghost cell values using linear extrapolation
c
      do 62 j=jmn+1,jmx
      if(imx.gt.imn+1) then
c
      rk1(imn  ,j,1) = 2.*rk1(imn+1,j,1) - rk1(imn+2,j,1)
      rk1(imn  ,j,2) = 2.*rk1(imn+1,j,2) - rk1(imn+2,j,2)
      rk1(imn  ,j,3) = 2.*rk1(imn+1,j,3) - rk1(imn+2,j,3)
      rk2(imn  ,j,1) = 2.*rk2(imn+1,j,1) - rk2(imn+2,j,1)
      rk2(imn  ,j,2) = 2.*rk2(imn+1,j,2) - rk2(imn+2,j,2)
      rk2(imn  ,j,3) = 2.*rk2(imn+1,j,3) - rk2(imn+2,j,3)
c
      rk1(imx+1,j,1) = 2.*rk1(imx  ,j,1) - rk1(imx-1,j,1)
      rk1(imx+1,j,2) = 2.*rk1(imx  ,j,2) - rk1(imx-1,j,2)
      rk1(imx+1,j,3) = 2.*rk1(imx  ,j,3) - rk1(imx-1,j,3)
      rk2(imx+1,j,1) = 2.*rk2(imx  ,j,1) - rk2(imx-1,j,1)
      rk2(imx+1,j,2) = 2.*rk2(imx  ,j,2) - rk2(imx-1,j,2)
      rk2(imx+1,j,3) = 2.*rk2(imx  ,j,3) - rk2(imx-1,j,3)
c
      else
c
      rk1(imx+1,j,1) = rk1(imx  ,j,1)
      rk1(imx+1,j,2) = rk1(imx  ,j,2)
      rk1(imx+1,j,3) = rk1(imx  ,j,3)
      rk2(imx+1,j,1) = rk2(imx  ,j,1)
      rk2(imx+1,j,2) = rk2(imx  ,j,2)
      rk2(imx+1,j,3) = rk2(imx  ,j,3)
c
      rk1(imn  ,j,1) = rk1(imn+1,j,1)
      rk1(imn  ,j,2) = rk1(imn+1,j,2)
      rk1(imn  ,j,3) = rk1(imn+1,j,3)
      rk2(imn  ,j,1) = rk2(imn+1,j,1)
      rk2(imn  ,j,2) = rk2(imn+1,j,2)
      rk2(imn  ,j,3) = rk2(imn+1,j,3)
c
      endif
c
   62 continue
c 
      do 64 i=imn,imx+1
      if(jmx.gt.jmn+1) then
c
      rk1(i,jmn  ,1) = 2.*rk1(i,jmn+1,1) - rk1(i,jmn+2,1)
      rk1(i,jmn  ,2) = 2.*rk1(i,jmn+1,2) - rk1(i,jmn+2,2)
      rk1(i,jmn  ,3) = 2.*rk1(i,jmn+1,3) - rk1(i,jmn+2,3)
      rk2(i,jmn  ,1) = 2.*rk2(i,jmn+1,1) - rk2(i,jmn+2,1)
      rk2(i,jmn  ,2) = 2.*rk2(i,jmn+1,2) - rk2(i,jmn+2,2)
      rk2(i,jmn  ,3) = 2.*rk2(i,jmn+1,3) - rk2(i,jmn+2,3)
c      
      rk1(i,jmx+1,1) = 2.*rk1(i,jmx  ,1) - rk1(i,jmx-1,1)
      rk1(i,jmx+1,2) = 2.*rk1(i,jmx  ,2) - rk1(i,jmx-1,2)
      rk1(i,jmx+1,3) = 2.*rk1(i,jmx  ,3) - rk1(i,jmx-1,3)
      rk2(i,jmx+1,1) = 2.*rk2(i,jmx  ,1) - rk2(i,jmx-1,1)
      rk2(i,jmx+1,2) = 2.*rk2(i,jmx  ,2) - rk2(i,jmx-1,2)
      rk2(i,jmx+1,3) = 2.*rk2(i,jmx  ,3) - rk2(i,jmx-1,3)
c
      else
c
      rk1(i,jmx+1,1) = rk1(i,jmx  ,1)
      rk1(i,jmx+1,2) = rk1(i,jmx  ,2)
      rk1(i,jmx+1,3) = rk1(i,jmx  ,3)
      rk2(i,jmx+1,1) = rk2(i,jmx  ,1)
      rk2(i,jmx+1,2) = rk2(i,jmx  ,2)
      rk2(i,jmx+1,3) = rk2(i,jmx  ,3)
c
      rk1(i,jmn  ,1) = rk1(i,jmn+1,1)
      rk1(i,jmn  ,2) = rk1(i,jmn+1,2)
      rk1(i,jmn  ,3) = rk1(i,jmn+1,3)
      rk2(i,jmn  ,1) = rk2(i,jmn+1,1)
      rk2(i,jmn  ,2) = rk2(i,jmn+1,2)
      rk2(i,jmn  ,3) = rk2(i,jmn+1,3)
c
      endif
c
   64 continue
c
c
      return
      end
