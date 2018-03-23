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
      subroutine segment ( d, x, y, z, x1, y1, z1 )
c  compute distance from point (x,y,z) to segment [(0,0,0);(x1,y1,z1)],
c  knowing that we are closer to the (0,0,0) end.
c

c  find out if projection falls outside segment. If so, the minimum distance
c  is to the end point, and we already have taken that into account.
      if ( x*x1 + y*y1 + z*z1 .lt. 0. ) return
c
c  now we know the projection falls within the segment; compute the distance
c  to the projection, by a cross product.
      anorm = sqrt ( x1**2 + y1**2 + z1**2 )
      if ( anorm .eq. 0. ) return
      xx = y*z1 - z*y1
      yy = z*x1 - x*z1
      zz = x*y1 - y*x1
      d = min ( d, sqrt ( xx**2 + yy**2 + zz**2 ) / anorm )
c
      return
c
      end
