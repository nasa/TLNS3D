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
      subroutine  deltaw(ncyc,nstage,igrid,imn,jmn,kmn,imx,jmx,kmx,
     .                   imp1,jmp1,kmp1,imp2,jmp2,kmp2,mx1dwk,w,p,
     .                   six,siy,siz,sjx,sjy,sjz,skx,sky,skz,fqs,
     .                   fs1,fs2,fs3,fs4,fs5,gs1,gs2,gs3,gs4,gs5,
     .                   dw1,dw2,dw3,dw4,dw5,hs1,hs2,hs3,hs4,hs5,delw,
     .                   iflag)
c=======================================================================
c 
c     compute convective fluxes witout dt scaling
c
c     w(i,j,k,1)  = density 
c     w(i,j,k,2)  = momentum in x direction 
c     w(i,j,k,3)  = momentum in y direction 
c     w(i,j,k,4)  = momentum in z direction 
c     w(i,j,k,5)  = total energy
c     p(i,j,k)    = pressure
c
c
c

      dimension      w(imp2,jmp2,kmp2,5), delw(imp2,jmp2,kmp2,5),
     .               p(imp2,jmp2,kmp2), six(imp1,jmp1,kmp1),
     .             siy(imp1,jmp1,kmp1), siz(imp1,jmp1,kmp1),
     .             sjx(imp1,jmp1,kmp1), sjy(imp1,jmp1,kmp1),
     .             sjz(imp1,jmp1,kmp1), skx(imp1,jmp1,kmp1),
     .             sky(imp1,jmp1,kmp1), skz(imp1,jmp1,kmp1)
c
c     dimension of temporary work arrays
c
      dimension     hs1(imp1,jmp1),hs2(imp1,jmp1),hs3(imp1,jmp1),
     .              hs4(imp1,jmp1),hs5(imp1,jmp1)
c
      dimension     fs1(imp1),fs2(imp1),fs3(imp1),fs4(imp1),fs5(imp1),
     .              gs1(imp1),gs2(imp1),gs3(imp1),gs4(imp1),gs5(imp1),
     .              dw1(imp1),dw2(imp1),dw3(imp1),dw4(imp1),dw5(imp1)
c
c     fqs = 0. =>  solid surface boundary
c         = 1. =>  any boundary execpt solid surface
c
      dimension     fqs(mx1dwk,mx1dwk,6)
c
c---------  initialize the fluxes 
c
      do 10 ia=1,5
      do 10 k=1,kmp2
      do 10 j=1,jmp2
      do 10 i=1,imp2
c
      delw(i,j,k,ia) = 0.
  10  continue
c
c
c--------- calculate flux at k=kmn plane
c
      do 20 j=jmn+1,jmx
      do 20 i=imn+1,imx
      ra        = 0.5*(w(i,j,kmn+1,1) + w(i,j,kmn,1))
      rua       = 0.5*(w(i,j,kmn+1,2) + w(i,j,kmn,2))
      rva       = 0.5*(w(i,j,kmn+1,3) + w(i,j,kmn,3))
      rwa       = 0.5*(w(i,j,kmn+1,4) + w(i,j,kmn,4))
      rea       = 0.5*(w(i,j,kmn+1,5) + w(i,j,kmn,5))
      pa        = 0.5*(p(i,j,kmn+1)   + p(i,j,kmn))
      qs        = fqs(i,j,5)*(rua*skx(i,j,kmn) + rva*sky(i,j,kmn)
     .                       +rwa*skz(i,j,kmn))/ra
      hs1(i,j)  = qs*ra
      hs2(i,j)  = qs*rua + pa*skx(i,j,kmn)
      hs3(i,j)  = qs*rva + pa*sky(i,j,kmn)
      hs4(i,j)  = qs*rwa + pa*skz(i,j,kmn)
      hs5(i,j)  = qs*(rea + pa)
   20 continue
c
c********* outer loop on k
c
      do 50 k=kmn+1,kmx
c
c--------- calculate flux at j=jmn plane
c
      do 30 i=imn+1,imx
      ra        = 0.5*(w(i,jmn+1,k,1) + w(i,jmn,k,1))
      rua       = 0.5*(w(i,jmn+1,k,2) + w(i,jmn,k,2))
      rva       = 0.5*(w(i,jmn+1,k,3) + w(i,jmn,k,3))
      rwa       = 0.5*(w(i,jmn+1,k,4) + w(i,jmn,k,4))
      rea       = 0.5*(w(i,jmn+1,k,5) + w(i,jmn,k,5))
      pa        = 0.5*(p(i,jmn+1,k)   + p(i,jmn,k))
      qs        = fqs(k,i,3)*(rua*sjx(i,jmn,k) + rva*sjy(i,jmn,k)
     .                       +rwa*sjz(i,jmn,k))/ra
      gs1(i)    = qs*ra
      gs2(i)    = qs*rua + pa*sjx(i,jmn,k)
      gs3(i)    = qs*rva + pa*sjy(i,jmn,k)
      gs4(i)    = qs*rwa + pa*sjz(i,jmn,k)
      gs5(i)    = qs*(rea + pa)
   30 continue
c
c********* outer loop on j
c
      do 50 j=jmn+1,jmx
c
c--------- calculate flux at i=imn face
c
      ra        = 0.5*(w(imn+1,j,k,1) + w(imn,j,k,1))
      rua       = 0.5*(w(imn+1,j,k,2) + w(imn,j,k,2))
      rva       = 0.5*(w(imn+1,j,k,3) + w(imn,j,k,3))
      rwa       = 0.5*(w(imn+1,j,k,4) + w(imn,j,k,4))
      rea       = 0.5*(w(imn+1,j,k,5) + w(imn,j,k,5))
      pa        = 0.5*(p(imn+1,j,k)   + p(imn,j,k))
      qs        = fqs(j,k,1)*(rua*six(imn,j,k) + rva*siy(imn,j,k)
     .                       +rwa*siz(imn,j,k))/ra
      fs1(imn)    = qs*ra
      fs2(imn)    = qs*rua + pa*six(imn,j,k)
      fs3(imn)    = qs*rva + pa*siy(imn,j,k)
      fs4(imn)    = qs*rwa + pa*siz(imn,j,k)
      fs5(imn)    = qs*(rea + pa)
c
c--------- flux in i direction
c
      do 40 i=imn+1,imx
      ra        = 0.5*(w(i+1,j,k,1) + w(i,j,k,1))
      rua       = 0.5*(w(i+1,j,k,2) + w(i,j,k,2))
      rva       = 0.5*(w(i+1,j,k,3) + w(i,j,k,3))
      rwa       = 0.5*(w(i+1,j,k,4) + w(i,j,k,4))
      rea       = 0.5*(w(i+1,j,k,5) + w(i,j,k,5))
      pa        = 0.5*(p(i+1,j,k)   + p(i,j,k))
      qs        = (rua*six(i,j,k) + rva*siy(i,j,k)
     .            +rwa*siz(i,j,k))/ra
      fimxb     = real(i/imx)
      qs        = (1. + (fqs(j,k,2) -1.)*fimxb)*qs
      fs1(i)    = qs*ra
      fs2(i)    = qs*rua + pa*six(i,j,k)
      fs3(i)    = qs*rva + pa*siy(i,j,k)
      fs4(i)    = qs*rwa + pa*siz(i,j,k)
      fs5(i)    = qs*(rea + pa)
   40 continue
c
c--------- accumulate flux in i direction 
c          and incoming fluxes in j and k directions
c
      do 42 i=imn+1,imx
      dw1(i)    = fs1(i) - fs1(i-1) - gs1(i) - hs1(i,j)
      dw2(i)    = fs2(i) - fs2(i-1) - gs2(i) - hs2(i,j)
      dw3(i)    = fs3(i) - fs3(i-1) - gs3(i) - hs3(i,j)
      dw4(i)    = fs4(i) - fs4(i-1) - gs4(i) - hs4(i,j)
      dw5(i)    = fs5(i) - fs5(i-1) - gs5(i) - hs5(i,j)
   42 continue
      do 44 i=imn+1,imx
c
c--------- flux in j direction
c
      ra        = 0.5*(w(i,j+1,k,1) + w(i,j,k,1))
      rua       = 0.5*(w(i,j+1,k,2) + w(i,j,k,2))
      rva       = 0.5*(w(i,j+1,k,3) + w(i,j,k,3))
      rwa       = 0.5*(w(i,j+1,k,4) + w(i,j,k,4))
      rea       = 0.5*(w(i,j+1,k,5) + w(i,j,k,5))
      pa        = 0.5*(p(i,j+1,k)   + p(i,j,k))
      qs        = (rua*sjx(i,j,k) + rva*sjy(i,j,k)
     .            +rwa*sjz(i,j,k))/ra
      fjmxb     = real(j/jmx)
      qs        = (1. + (fqs(k,i,4) -1.)*fjmxb)*qs
      gs1(i)    = qs*ra
      gs2(i)    = qs*rua + pa*sjx(i,j,k)
      gs3(i)    = qs*rva + pa*sjy(i,j,k)
      gs4(i)    = qs*rwa + pa*sjz(i,j,k)
      gs5(i)    = qs*(rea + pa)
c
c
c--------- flux in k direction
c
      ra        = 0.5*(w(i,j,k+1,1) + w(i,j,k,1))
      rua       = 0.5*(w(i,j,k+1,2) + w(i,j,k,2))
      rva       = 0.5*(w(i,j,k+1,3) + w(i,j,k,3))
      rwa       = 0.5*(w(i,j,k+1,4) + w(i,j,k,4))
      rea       = 0.5*(w(i,j,k+1,5) + w(i,j,k,5))
      pa        = 0.5*(p(i,j,k+1)   + p(i,j,k))
      qs        = (rua*skx(i,j,k) + rva*sky(i,j,k)
     .            +rwa*skz(i,j,k))/ra
      fkmxb     = real(k/kmx)
      qs        = (1. + (fqs(i,j,6) -1.)*fkmxb)*qs
      hs1(i,j)  = qs*ra
      hs2(i,j)  = qs*rua + pa*skx(i,j,k)
      hs3(i,j)  = qs*rva + pa*sky(i,j,k)
      hs4(i,j)  = qs*rwa + pa*skz(i,j,k)
      hs5(i,j)  = qs*(rea + pa)
   44 continue
c
c********* calculate residuals
c
      do 50 i=imn+1,imx
      delw(i,j,k,1) = dw1(i) + gs1(i) + hs1(i,j)
      delw(i,j,k,2) = dw2(i) + gs2(i) + hs2(i,j)
      delw(i,j,k,3) = dw3(i) + gs3(i) + hs3(i,j)
      delw(i,j,k,4) = dw4(i) + gs4(i) + hs4(i,j)
      delw(i,j,k,5) = dw5(i) + gs5(i) + hs5(i,j)
c
   50 continue
c
c 
      return
      end
