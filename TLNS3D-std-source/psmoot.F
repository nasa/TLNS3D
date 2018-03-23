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
      subroutine psmoot (imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                   imp2,jmp2,kmp2,wk2da,wk2db,dw,
     .                   smoopi,smoopj,smoopk)
c=======================================================================
c 
c           implicit residual smoothing with constant coefficients
c
c

      dimension dw(imp2,jmp2,kmp2,1),
     .          wk2da(imx,jmx), wk2db(imx,kmx)
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c
c
      imnp = imn+1
      jmnp = jmn+1
      kmnp = kmn+1
c
c--------- smoothing in i direction 
c
      if(smoopi.le.0.00001) return
      a           = smoopi
      do 40 k=kmnp,kmx
      do 10 j=jmnp,jmx
      t           = 1./(1.+a+a)
c     neumann b.c.
      t           = 1./(1.+a)
      wk2da(imnp,j)  = t*a
      dw(imnp,j,k,1) = t*dw(imnp,j,k,1) 
   10 continue
      do 20 i=imn+2,imx
cdir$ ivdep
      do 20 j=jmnp ,jmx
      t           = 1./(1.+a+a -a*wk2da(i-1,j))
      wk2da(i,j)  = t*a
      dw(i,j,k,1) = t*(dw(i,j,k,1) + a*dw(i-1,j,k,1)) 
c
c     neumann b.c.
      if(i.eq.imx) then
         dw(i,j,k,1) = dw(i,j,k,1)/(1.-wk2da(i,j))
      endif
   20 continue
      do 40 i=imx-1,imnp,-1
      do 40 j=jmnp ,jmx
      dw(i,j,k,1) = dw(i,j,k,1) + wk2da(i,j)*dw(i+1,j,k,1)
   40 continue
c
c--------- smoothing in j direction 
c
      a           = smoopj
      do 70 k=kmnp,kmx
      do 50 i=imnp,imx
      t           = 1./(1.+a+a)
c     neumann b.c.
      t           = 1./(1.+a)
      wk2da(i,jmnp)  = t*a
      dw(i,jmnp,k,1) = t*dw(i,jmnp,k,1) 
   50 continue
      do 60 j=jmn+2,jmx
cdir$ ivdep
      do 60 i=imnp ,imx
      t           = 1./(1.+a+a -a*wk2da(i,j-1))
cmds  t           = 1./(1.+a+a)
      wk2da(i,j)  = t*a
      dw(i,j,k,1) = t*(dw(i,j,k,1) + a*dw(i,j-1,k,1)) 
c
c     neumann b.c.
      if(j.eq.jmx) then
         dw(i,j,k,1) = dw(i,j,k,1)/(1.-wk2da(i,j))
      endif
   60 continue
      do 70 j=jmx-1,jmnp,-1
cdir$ ivdep
      do 70 i=imnp ,imx
      dw(i,j,k,1) = dw(i,j,k,1) + wk2da(i,j)*dw(i,j+1,k,1)
   70 continue
c
c--------- smoothing in k direction 
c
      if(i2dfl.eq.1) return
c
      a           = smoopk
      do 100 j=jmnp,jmx
      do 80  i=imnp,imx
      t           = 1./(1.+a+a)
c
c     neumann b.c.
      t           = 1./(1.+a)
      wk2db(i,kmnp)  = t*a
      dw(i,j,kmnp,1) = t*dw(i,j,kmnp,1) 
   80 continue
      do 90 k=kmn+2,kmx
cdir$ ivdep
      do 90 i=imnp ,imx
      t           = 1./(1.+a+a-a*wk2db(i,k-1))
      wk2db(i,k)  = t*a
      dw(i,j,k,1) = t*(dw(i,j,k,1) + a*dw(i,j,k-1,1)) 
c
c     neumann b.c.
      if(k.eq.kmx) then
         dw(i,j,k,1) = dw(i,j,k,1)/(1.-wk2db(i,k))
      endif
   90 continue
      do 100 k=kmx-1,kmnp,-1
cdir$ ivdep
      do 100 i=imnp ,imx
      dw(i,j,k,1) = dw(i,j,k,1) + wk2db(i,k)*dw(i,j,k+1,1)
  100 continue
c
c
      return
      end
