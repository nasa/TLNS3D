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
      subroutine collect_surf (xs,
     .                  n1beg ,n1end ,n2beg ,n2end ,
     .                  nsurf, surf, isurf,
     .                  ntri, iptri )
c============================================================================

      dimension     xs(3,*),surf(nsurf,3)
      dimension     ntri(nsurf),iptri(nsurf,8)
c
c Collect surface points.
c Now add collection of local triangles.
c Important side effect is that isurf is incremented properly
c
      isurf_beg=isurf
      ipts = n1end - n1beg + 1
c
      do k=n2beg,n2end
         do j=n1beg,n1end
            isurf=isurf+1
            i = 1 + (j-n1beg) + ipts*(k-n2beg)
            surf(isurf,1)=xs(1,i)
            surf(isurf,2)=xs(2,i)
            surf(isurf,3)=xs(3,i)
         end do
      end do
c
      isurf=isurf_beg
      n1inc=1
      n2inc=n1end-n1beg+1
      do j=n2beg,n2end
      do i=n1beg,n1end
        isurf=isurf+1
        itri=0
        do jj=max(j-1,n2beg),min(j+1,n2end)
        do ii=max(i-1,n1beg),min(i+1,n1end)
          if (ii .ne. i .and. jj .ne. j) then
          itri=itri+1
          iptri(isurf,2*itri-1) = isurf+(ii-i)*n1inc
          iptri(isurf,2*itri)   = isurf+(jj-j)*n2inc
          end if
        end do
        end do
        ntri(isurf)=itri
      end do
      end do
      return
      end
