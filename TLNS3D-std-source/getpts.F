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
c=====================================================================
      subroutine getpts (imn,jmn,kmn,imx,jmx,kmx,
     .                   imp1,jmp1,kmp1,x,mx1dwk,
     .                   nface,n1beg,n1end,n2beg,n2end,xs)
c=====================================================================
c
c     this subroutine collects the solid surface points from a
c     given segment on a block into the array xs.
c     The points are put in the array as triplets (x,y,z) so
c     that the addressing is xs(3,npts).
c

      dimension x(imp1,jmp1,kmp1,3),xs(3,mx1dwk*mx1dwk)
c
      ipts    = n1end - n1beg + 1
      jpts    = n2end - n2beg + 1
c
c---- branch based on block face number
c
      go to (100,200,300,400,500,600) nface
c
  100 iw = imn
      go to 201
  200 iw = imx
  201 do kw=n2beg,n2end
        do jw=n1beg,n1end
          i = 1 +jw-n1beg +ipts*(kw-n2beg)
          do n=1,3
            xs(n,i) = x(iw,jw,kw,n)
          enddo
        enddo
      enddo
      return
c
  300 jw = jmn
      go to 401
  400 jw = jmx
  401 do iw=n2beg,n2end
        do kw=n1beg,n1end
          i = 1 +kw-n1beg +ipts*(iw-n2beg)
          do n=1,3
            xs(n,i) = x(iw,jw,kw,n)
          enddo
        enddo
      enddo
      return
c
  500 kw = kmn
      go to 601
  600 kw = kmx
  601 do jw=n2beg,n2end
        do iw=n1beg,n1end
          i = 1 +iw-n1beg +ipts*(jw-n2beg)
          do n=1,3
            xs(n,i) = x(iw,jw,kw,n)
          enddo
        enddo
      enddo
      return
      end
