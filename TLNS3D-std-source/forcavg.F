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
      subroutine forcavg(itime,clg,cdtg,cdvg,cmzg,
     .                   clavg,cdtavg,cdvavg,cmzavg )
c
c     object : average force coefficients for unsteady flows
c
c     initial coding : v.n.vatsa (June 2001)
c

c
      if(itime.eq.1) then
        clavg  = clg
        cdtavg = cdtg
        cdvavg = cdvg
        cmzavg = cmzg
c
      else
c
        wf1  = real(itime-1)/real(itime)
        wf2  = 1. - wf1
        clavg  = wf1*clavg  + wf2*clg
        cdtavg = wf1*cdtavg + wf2*cdtg
        cdvavg = wf1*cdvavg + wf2*cdvg
        cmzavg = wf1*cmzavg + wf2*cmzg
c
      endif
c
      return
      end
