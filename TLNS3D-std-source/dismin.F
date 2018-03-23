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
      subroutine dismin (imn,jmn,kmn,imx,jmx,kmx,
     .                   imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .                   x,smin,ipts,jpts,xs,mxsizc,d2)
c=======================================================================
c
c
c    object : compute distance at node points from the closest
c             solid surface
c
c    initial coding : by v.n.vatsa of NASA Langley (March, 1993)
c
c
c    modified by: P. Spalart and W. Wilkinson  (September, 1993)
c    to reflect the true distance to the wall, even with skewed grids
c    This mod is kludgy in terms of the work array d2.
c    This can be fixed by somebody that knows the code structure and available
c    work space.
c
c    modified for distributed computing in March 1995
c    d2 is passed as wk3d from caller which is dimension mxsizc which
c    is the same size as smin. xs is now a 2-D array of x,y,z triplets
c    containing the surface points.
c

      dimension x(imp1,jmp1,kmp1,3),smin(imp2,jmp2,kmp2)
      dimension xs(3,ipts,jpts),d2(mxsizc)
c
c  loop on field point
c
      do 120 k = kmn,kmx
      do 120 j = jmn,jmx
      do 120 i = imn,imx
      xnc       = x(i,j,k,1)
      ync       = x(i,j,k,2)
      znc       = x(i,j,k,3)
c
c  compute distance squared to all wall points
c  put it into a 1D array for long vector and good use of work space
c
      do 110 jw = 1,jpts
        do 110 iw = 1,ipts
 110      d2(iw+ipts*(jw-1)) = (xnc-xs(1,iw,jw))**2
     .                        +(ync-xs(2,iw,jw))**2
     .                        +(znc-xs(3,iw,jw))**2
c
c  find closest wall grid point within this face.
c
      ijw = ismin ( ipts * jpts, d2, 1 )
      dd = sqrt ( d2(ijw) )
      iw = mod ( ijw, ipts )
      jw = ijw / ipts + 1
      if (iw.eq.0) then
        iw = ipts
        jw = jw - 1
      endif
c
c  find distance to the four triangles that have that point as an apex
c
      if (iw.lt.ipts .and. jw.lt.jpts)
     .  call triang ( dd, xnc, ync, znc,
     .    xs(1,iw,jw),   xs(2,iw,jw),   xs(3,iw,jw),
     .    xs(1,iw+1,jw), xs(2,iw+1,jw), xs(3,iw+1,jw),
     .    xs(1,iw,jw+1), xs(2,iw,jw+1), xs(3,iw,jw+1) )
      if (iw.lt.ipts .and. jw.gt.1)
     .  call triang ( dd, xnc, ync, znc,
     .    xs(1,iw,jw),   xs(2,iw,jw),   xs(3,iw,jw),
     .    xs(1,iw,jw-1), xs(2,iw,jw-1), xs(3,iw,jw-1),
     .    xs(1,iw+1,jw), xs(2,iw+1,jw), xs(3,iw+1,jw) )
      if (iw.gt.1 .and. jw.lt.jpts)
     .  call triang ( dd, xnc, ync, znc,
     .    xs(1,iw,jw),   xs(2,iw,jw),   xs(3,iw,jw),
     .    xs(1,iw,jw+1), xs(2,iw,jw+1), xs(3,iw,jw+1),
     .    xs(1,iw-1,jw), xs(2,iw-1,jw), xs(3,iw-1,jw) )
      if (iw.gt.1 .and. jw.gt.1)
     .  call triang ( dd, xnc, ync, znc,
     .    xs(1,iw,jw),   xs(2,iw,jw),   xs(3,iw,jw),
     .    xs(1,iw-1,jw), xs(2,iw-1,jw), xs(3,iw-1,jw),
     .    xs(1,iw,jw-1), xs(2,iw,jw-1), xs(3,iw,jw-1) )
c
c  put the result into the distance array
c
      smin(i,j,k) = min ( smin(i,j,k), dd )
  120 continue
c
      return
      end
