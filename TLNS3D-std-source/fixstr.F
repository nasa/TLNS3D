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
      subroutine fixstr (str, len)
c
      character*1 str(len),quote,blank
      data quote/''''/, blank/' '/
c
      do n=1,len
        if (str(n).ne.blank) go to 10
      enddo
      return
   10 if (str(n).eq.quote) n = n + 1
      if (n.gt.1) then
        j = 1
        do i=n,len
          str(j) = str(i)
          j = j + 1
        enddo
      endif
      do n=1,len
        if (str(n).eq.quote .or. str(n).eq.blank) go to 20
      enddo
      return
   20 do i=n,len
        str(i) = blank
      enddo
      return
      end
