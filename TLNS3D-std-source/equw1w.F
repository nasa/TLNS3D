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
      subroutine equw1w(imn,jmn,kmn,imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2,w,w1,nvar)
c
c     object : equate 'w1' to 'w'
c
c     initial coding : m.d.sanetrik of as&m, dec. 1990
c

      dimension w(imp2,jmp2,kmp2,nvar),w1(imp2,jmp2,kmp2,nvar)
c
c
      do 10 nv=1,nvar
      do 10 k = kmn,kmp1
      do 10 j = jmn,jmp1
      do 10 i = imn,imp1
      w1(i,j,k,nv) = w(i,j,k,nv)
  10  continue
c
c
      return
      end
