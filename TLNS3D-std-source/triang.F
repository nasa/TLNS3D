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
      subroutine triang(dist,pp1,pp2,pp3,
     . aa1,aa2,aa3,bb1,bb2,bb3,cc1,cc2,cc3)
      a1=bb1-aa1
      a2=bb2-aa2
      a3=bb3-aa3
      b1=cc1-aa1
      b2=cc2-aa2
      b3=cc3-aa3
      p1=pp1-aa1
      p2=pp2-aa2
      p3=pp3-aa3
      daa=a1*a1+a2*a2+a3*a3
      dab=a1*b1+a2*b2+a3*b3
      dbb=b1*b1+b2*b2+b3*b3
      den=dab**2-daa*dbb
      if (den .eq. 0) go to 100
      dap=a1*p1+a2*p2+a3*p3
      dbp=b1*p1+b2*p2+b3*p3
      s=(dab*dbp-dbb*dap)/den
      t=(dab*dap-daa*dbp)/den
      if (s .lt. 0 .or. t .lt. 0 .or. t+s .gt. 1) go to 100
      r1=p1-s*a1-t*b1
      r2=p2-s*a2-t*b2
      r3=p3-s*a3-t*b3
      dist=min(dist,sqrt(r1**2+r2**2+r3**2))
      return
 100  continue
      dsq=dist*dist
      if (daa .eq. 0) go to 200
      dap=a1*p1+a2*p2+a3*p3
      t=dap/daa
      if (t .lt. 0 .or. t .gt. 1) go to 200
      r1=p1-t*a1
      r2=p2-t*a2
      r3=p3-t*a3
      dsq=min(dsq,r1**2+r2**2+r3**2)
 200  continue
      if (dbb .eq. 0) go to 300
      dbp=b1*p1+b2*p2+b3*p3
      t=dbp/dbb
      if (t .lt. 0. .or. t .gt. 1) go to 300
      r1=p1-t*b1
      r2=p2-t*b2
      r3=p3-t*b3
      dsq=min(dsq,r1**2+r2**2+r3**2)
 300  continue
      p1=pp1-bb1
      p2=pp2-bb2
      p3=pp3-bb3
      a1=cc1-bb1
      a2=cc2-bb2
      a3=cc3-bb3
      daa=a1*a1+a2*a2+a3*a3
      if (daa .eq. 0) go to 400
      dap=a1*p1+a2*p2+a3*p3
      t=dap/daa
      if (t .lt. 0 .or. t .gt. 1) go to 400
      r1=p1-t*a1
      r2=p2-t*a2
      r3=p3-t*a3
      dsq=min(dsq,r1**2+r2**2+r3**2)
400   continue
      if (dsq .lt. dist*dist) dist=sqrt(dsq)
      return
      end
