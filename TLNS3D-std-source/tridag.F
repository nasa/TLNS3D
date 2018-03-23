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
       subroutine tridag(idim,jdim,imin,imax,jmin,jmax,x,a,b,c,f)
c
c      object : to solve the tridiagonal equations (scalar)
c               equations are implicit in i- direction
c
c      initial cooding : by v.n.vatsa (nov. 1992)
c
c      a,b and c are the lower,middle and upper diagonal coefficients
c      f is the r.h.s when this routine is called
c      (f is replaced by the solution before exiting this routine)
c      x is working array
c

       dimension x(idim,jdim),a(idim,jdim),b(idim,jdim),
     .           c(idim,jdim),f(idim,jdim)
c
       do 10 j=jmin,jmax
       x(imin,j)   =  c(imin,j)/b(imin,j)
       f(imin,j)   =  f(imin,j)/b(imin,j)
  10   continue
c
       imnp1     = imin + 1
       do 100 ii = imnp1,imax
c
         do 110 j  = jmin,jmax
         d         = 1./(b(ii,j) - a(ii,j)*x(ii-1,j))
         x(ii,j)   = c(ii,j)*d
         f(ii,j)   =(f(ii,j) - a(ii,j)*f(ii-1,j))*d
  110    continue
c
  100  continue
c
       itotal = imin + imax
       do 200 ii = imnp1,imax
         i    = itotal - ii
         do 210 j  = jmin,jmax
         f(i,j)    = f(i,j) - x(i,j)*f(i+1,j)
  210    continue
c
  200  continue
c
       return
       end
