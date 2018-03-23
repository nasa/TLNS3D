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
      subroutine hdamp(imn, jmn, kmn, imx, jmx, kmx,
     .                 imp1,jmp1,kmp1,imp2,jmp2,kmp2,hm,w,p )
c
c     object : use enthalpy-damping technique to enhance
c              numerical stability
c
c     initial coding : v.n.vatsa of nasa langley, Feb. 1990
c
c

      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension    w(imp2,jmp2,kmp2,5),  p(imp2,jmp2,kmp2)
c
c
      do 10 k = kmn+1,kmx
      do 10 j = jmn+1,jmx
      do 10 i = imn+1,imx
      f          = 1./(1. +hm*((w(i,j,k,5) +p(i,j,k))/w(i,j,k,1) -h0))
      w(i,j,k,1) = f*w(i,j,k,1)
      w(i,j,k,2) = f*w(i,j,k,2)
      w(i,j,k,3) = f*w(i,j,k,3)
      w(i,j,k,4) = f*w(i,j,k,4)
      w(i,j,k,5) =  (w(i,j,k,5)  - hm*p(i,j,k))/(1. + hm)
     .           +h0*w(i,j,k,1)*(1. - 1./(f*(1.+hm)))
  10  continue
c
c     evaluate the pressure distribution
c
      do 20 k = kmn+1,kmx
      do 20 j = jmn+1,jmx
      do 20 i = imn+1,imx
c
      qq          = (w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .               w(i,j,k,4)**2 )/w(i,j,k,1)
      p(i,j,k)    =  gm1*(w(i,j,k,5) - .5*qq )
  20  continue
c
c
      return
      end
