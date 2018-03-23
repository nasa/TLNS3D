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
      subroutine flowavg(imn,jmn,kmn,imp2,jmp2,kmp2,
     .                   itime,p,pavg,w,wavg,eomu,eomuavg)
c
c     object : average flow variables for unsteady flows
c
c     initial coding : v.n.vatsa (June 2001)
c

      dimension w(imp2,jmp2,kmp2,5),wavg(imp2,jmp2,kmp2,5)
      dimension p(imp2,jmp2,kmp2)  ,pavg(imp2,jmp2,kmp2)
      dimension eomu(imp2,jmp2,kmp2)  ,eomuavg(imp2,jmp2,kmp2)
c
c
      if(itime.eq.1) then
        do k = 1,kmp2
          do j = 1,jmp2
            do i = 1,imp2
            wavg(i,j,k,1)  = w(i,j,k,1)
            wavg(i,j,k,2)  = w(i,j,k,2)
            wavg(i,j,k,3)  = w(i,j,k,3)
            wavg(i,j,k,4)  = w(i,j,k,4)
            wavg(i,j,k,5)  = w(i,j,k,5)
            pavg(i,j,k  )  = p(i,j,k  )
            eomuavg(i,j,k) = eomu(i,j,k)
            enddo
          enddo
        enddo
c
      else
c
        wf1  = real(itime-1)/real(itime)
        wf2  = 1. - wf1
        do k = 1,kmp2
          do j = 1,jmp2
            do i = 1,imp2
            wavg(i,j,k,1)  = wf1*wavg(i,j,k,1) + wf2*w(i,j,k,1)
            wavg(i,j,k,2)  = wf1*wavg(i,j,k,2) + wf2*w(i,j,k,2)
            wavg(i,j,k,3)  = wf1*wavg(i,j,k,3) + wf2*w(i,j,k,3)
            wavg(i,j,k,4)  = wf1*wavg(i,j,k,4) + wf2*w(i,j,k,4)
            wavg(i,j,k,5)  = wf1*wavg(i,j,k,5) + wf2*w(i,j,k,5)
            pavg(i,j,k  )  = wf1*pavg(i,j,k  ) + wf2*p(i,j,k  )
            eomuavg(i,j,k) = wf1*eomuavg(i,j,k) + wf2*eomu(i,j,k)
            enddo
          enddo
        enddo
c
      endif
c
      return
      end
