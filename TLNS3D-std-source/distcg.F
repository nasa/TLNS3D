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
      subroutine distcg(imn,jmn,kmn,imf,jmf,kmf,imp1f,jmp1f,kmp1f,
     .               imp2f,jmp2f,kmp2f, imc,jmc,kmc,
     .               imp1c,jmp1c,kmp1c,imp2c,jmp2c,kmp2c,sminf,sminc)
c
c     object : given the distance from wall on a given size mesh,
c              pick corresponding distances on one level coarser mesh 
c              sminf and sminc are the distances on fine & coarse grids
c
c     initial coding : by v.n.vatsa of NASA Langley (March, 1993)
c

      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension sminf(imp2f,jmp2f,kmp2f),sminc(imp2c,jmp2c,kmp2c)
c
c
c     accomodate two cell faces at each boundary
c
      if(i2dfl.eq.1) go to 100
c
        do 10 kf=kmn,kmf,2
      kc= (kf+2)/2
        do 10 jf=jmn,jmf,2
      jc= (jf+2)/2
        do 10 if=imn,imf,2
      ic= (if+2)/2
c
      sminc(ic,jc,kc)=sminf(if,jf,kf)
c
   10 continue
c
      return
c
  100 continue
c
c     coarse grid for 2-d flows
c
        do 20 kc=kmn,kmf
        do 20 jf=jmn,jmf,2
      jc= (jf+2)/2
        do 20 if=imn,imf,2
      ic= (if+2)/2
c
      sminc(ic,jc,kc)=sminf(if,jf,kmn)
c
   20 continue
      return
      end
