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
      subroutine equate (imn,jmn,kmn,imx,jmx,kmx,
     .                   imp2,jmp2,kmp2,w,ws)
c
c     object : equate 'ws' to 'w'
c
c     initial coding : v.n.vatsa of nasa langley, Feb. 1990
c

      dimension w(imp2,jmp2,kmp2,5),ws(imp2,jmp2,kmp2,5)
c
c
      do 10 nv=1,5
      do 10 k = kmn+1,kmx
      do 10 j = jmn+1,jmx
      do 10 i = imn+1,imx
      ws(i,j,k,nv) = w(i,j,k,nv)
  10  continue
c
c
      return
      end
