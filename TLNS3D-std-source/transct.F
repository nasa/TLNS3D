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
      subroutine transct (imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,w,p,delw )
c
c     transform from conserved variables to "p,u,v,w,T" variables
c     
c
c     initial coding : v.n.vatsa of nasa langley, July, 1995
c
c
c-------   preconditioning related quantities      -----------------
c
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubfil2,iprec
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension      w(imp2,jmp2,kmp2,5), delw(imp2,jmp2,kmp2,5),
     .               p(imp2,jmp2,kmp2)
c
c
      do 10 k = kmn+1,kmx
      do 10 j = jmn+1,jmx
      do 10 i = imn+1,imx
        rhoi          = 1./w(i,j,k,1)
        ua            = w(i,j,k,2)*rhoi
        va            = w(i,j,k,3)*rhoi
        wa            = w(i,j,k,4)*rhoi
        q2            = 0.5*(ua**2 + va**2 + wa**2)
        ts            = p(i,j,k)*rhoi
c
c
        dw1           = gm1*(q2*delw(i,j,k,1) - ua*delw(i,j,k,2)
     .                -       va*delw(i,j,k,3) -wa*delw(i,j,k,4)
     .                +          delw(i,j,k,5) )
        delw(i,j,k,2) =(delw(i,j,k,2) - ua*delw(i,j,k,1))*rhoi
        delw(i,j,k,3) =(delw(i,j,k,3) - va*delw(i,j,k,1))*rhoi
        delw(i,j,k,4) =(delw(i,j,k,4) - wa*delw(i,j,k,1))*rhoi
        dw5           =(dw1 - ts*delw(i,j,k,1))*rhoi
c
        delw(i,j,k,1) = dw1
        delw(i,j,k,5) = dw5
  10  continue
c
c
      return
      end
