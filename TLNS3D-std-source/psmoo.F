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
      subroutine psmoo(ncyc,nstage,imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,
     .                 kmp1,imp2,jmp2,kmp2,d,eprs,dw)
c=======================================================================
c 
c              implicit residual smoothing
c
c 
c

      dimension eprs(imp2,jmp2,kmp2,3), dw(imp2,jmp2,kmp2,5),
     .          d(imx,jmx)
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      imnp  = imn +1
      jmnp  = jmn +1
      kmnp  = kmn +1
      imnp2 = imn +2
      jmnp2 = jmn +2
      kmnp2 = kmn +2
c
c--------- smoothing in i direction 
c
c     call flowmark('loop 40'l)
      do 40 k=kmnp,kmx
c     call flowmark('loop 10'l)
      do 10 j=jmnp,jmx
      a              = eprs(imnp,j,k,1)
      t              = 1./(1.+a+a)
      d(imnp,j)      = t*a
      dw(imnp,j,k,1) = t*dw(imnp,j,k,1) 
      dw(imnp,j,k,2) = t*dw(imnp,j,k,2) 
      dw(imnp,j,k,3) = t*dw(imnp,j,k,3) 
      dw(imnp,j,k,4) = t*dw(imnp,j,k,4) 
      dw(imnp,j,k,5) = t*dw(imnp,j,k,5) 
   10 continue
c     call flowmark(0)
c
c     call flowmark('loop 20'l)
      do 20 i=imnp2,imx
cdir$ ivdep
      do 20 j=jmnp ,jmx
      a           = eprs(i,j,k,1)
      t           = 1./(1.+a+a -a*d(i-1,j))
      d(i,j)      = t*a
      dw(i,j,k,1) = t*(dw(i,j,k,1) + a*dw(i-1,j,k,1)) 
      dw(i,j,k,2) = t*(dw(i,j,k,2) + a*dw(i-1,j,k,2)) 
      dw(i,j,k,3) = t*(dw(i,j,k,3) + a*dw(i-1,j,k,3)) 
      dw(i,j,k,4) = t*(dw(i,j,k,4) + a*dw(i-1,j,k,4)) 
      dw(i,j,k,5) = t*(dw(i,j,k,5) + a*dw(i-1,j,k,5)) 
   20 continue
c     call flowmark(0)
c
      do 40 i=imx-1,imnp,-1
      do 40 j=jmnp ,jmx
      dw(i,j,k,1) = dw(i,j,k,1) + d(i,j)*dw(i+1,j,k,1)
      dw(i,j,k,2) = dw(i,j,k,2) + d(i,j)*dw(i+1,j,k,2)
      dw(i,j,k,3) = dw(i,j,k,3) + d(i,j)*dw(i+1,j,k,3)
      dw(i,j,k,4) = dw(i,j,k,4) + d(i,j)*dw(i+1,j,k,4)
      dw(i,j,k,5) = dw(i,j,k,5) + d(i,j)*dw(i+1,j,k,5)
   40 continue
c     call flowmark(0)
c
c--------- smoothing in j direction 
c
c     call flowmark('loop 70'l)
      do 70 k=kmnp,kmx
      do 50 i=imnp,imx
      a              = eprs(i,jmnp,k,2)
      t              = 1./(1.+a+a)
      d(i,jmnp)      = t*a
      dw(i,jmnp,k,1) = t*dw(i,jmnp,k,1) 
      dw(i,jmnp,k,2) = t*dw(i,jmnp,k,2) 
      dw(i,jmnp,k,3) = t*dw(i,jmnp,k,3) 
      dw(i,jmnp,k,4) = t*dw(i,jmnp,k,4) 
      dw(i,jmnp,k,5) = t*dw(i,jmnp,k,5) 
   50 continue
c
      do 60 j=jmnp2,jmx
cdir$ ivdep
      do 60 i=imnp ,imx
      a           = eprs(i,j,k,2)
      t           = 1./(1.+a+a -a*d(i,j-1))
      d(i,j)      = t*a
      dw(i,j,k,1) = t*(dw(i,j,k,1) + a*dw(i,j-1,k,1)) 
      dw(i,j,k,2) = t*(dw(i,j,k,2) + a*dw(i,j-1,k,2)) 
      dw(i,j,k,3) = t*(dw(i,j,k,3) + a*dw(i,j-1,k,3)) 
      dw(i,j,k,4) = t*(dw(i,j,k,4) + a*dw(i,j-1,k,4)) 
      dw(i,j,k,5) = t*(dw(i,j,k,5) + a*dw(i,j-1,k,5)) 
   60 continue
c
      do 70 j=jmx-1,jmnp,-1
cdir$ ivdep
      do 70 i=imnp ,imx
      dw(i,j,k,1) = dw(i,j,k,1) + d(i,j)*dw(i,j+1,k,1)
      dw(i,j,k,2) = dw(i,j,k,2) + d(i,j)*dw(i,j+1,k,2)
      dw(i,j,k,3) = dw(i,j,k,3) + d(i,j)*dw(i,j+1,k,3)
      dw(i,j,k,4) = dw(i,j,k,4) + d(i,j)*dw(i,j+1,k,4)
      dw(i,j,k,5) = dw(i,j,k,5) + d(i,j)*dw(i,j+1,k,5)
   70 continue
c     call flowmark(0)
c
c--------- smoothing in k direction 
c
c     call flowmark('loop100'l)
c
      if(i2dfl.eq.1) return
c
      do 100 j=jmnp,jmx
      do 80  i=imnp,imx
      a              = eprs(i,j,kmnp,3)
      t              = 1./(1.+a+a)
c
c     neumann condition
      t              = 1./(1.+a)
      d(i,kmnp)      = t*a
      dw(i,j,kmnp,1) = t*dw(i,j,kmnp,1) 
      dw(i,j,kmnp,2) = t*dw(i,j,kmnp,2) 
      dw(i,j,kmnp,3) = t*dw(i,j,kmnp,3) 
      dw(i,j,kmnp,4) = t*dw(i,j,kmnp,4) 
      dw(i,j,kmnp,5) = t*dw(i,j,kmnp,5) 
   80 continue
c
      do 90 k=kmnp2,kmx
cdir$ ivdep
      do 90 i=imnp ,imx
      a           = eprs(i,j,k,3)
      t           = 1./(1.+a+a-a*d(i,k-1))
      d(i,k)      = t*a
      dw(i,j,k,1) = t*(dw(i,j,k,1) + a*dw(i,j,k-1,1)) 
      dw(i,j,k,2) = t*(dw(i,j,k,2) + a*dw(i,j,k-1,2)) 
      dw(i,j,k,3) = t*(dw(i,j,k,3) + a*dw(i,j,k-1,3)) 
      dw(i,j,k,4) = t*(dw(i,j,k,4) + a*dw(i,j,k-1,4)) 
      dw(i,j,k,5) = t*(dw(i,j,k,5) + a*dw(i,j,k-1,5)) 
c
c     neumann b.c.
      if(k.eq.kmx) then
         dw(i,j,k,1) = dw(i,j,k,1)/(1.-d(i,k))
         dw(i,j,k,2) = dw(i,j,k,2)/(1.-d(i,k))
         dw(i,j,k,3) = dw(i,j,k,3)/(1.-d(i,k))
         dw(i,j,k,4) = dw(i,j,k,4)/(1.-d(i,k))
         dw(i,j,k,5) = dw(i,j,k,5)/(1.-d(i,k))
      endif
   90 continue
c
      do 100 k=kmx-1,kmnp,-1
cdir$ ivdep
      do 100 i=imnp ,imx
      dw(i,j,k,1) = dw(i,j,k,1) + d(i,k)*dw(i,j,k+1,1)
      dw(i,j,k,2) = dw(i,j,k,2) + d(i,k)*dw(i,j,k+1,2)
      dw(i,j,k,3) = dw(i,j,k,3) + d(i,k)*dw(i,j,k+1,3)
      dw(i,j,k,4) = dw(i,j,k,4) + d(i,k)*dw(i,j,k+1,4)
      dw(i,j,k,5) = dw(i,j,k,5) + d(i,k)*dw(i,j,k+1,5)
  100 continue
c     call flowmark(0)
c
c
      return
      end
