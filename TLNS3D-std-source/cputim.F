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
       subroutine  cputim(tim)
c=======================================================================
c
c dana 060895: changed time to always be SP; dimension to real(4)

       real*4 tim, tm, y
       dimension tim(3,3),tm(3),ia(3)
c
c     etime returns elapsed time as:
c       tm(1) = user time
c       tm(2) = system time
c     itime returns wall clock time as:
c       ia(1) = hour (0-23)
c       ia(2) = minute (0-59)
c       ia(3) = second (0-59)
c
       y = etime (tm)
       call itime (ia)
       tm(3) = ia(3) +60*(ia(2) +60*ia(1))
c
c     first time here
c
       tm(1) = max(tm(1), 1.e-7)
       if (tim(1,3).lt.1.e-8) then
         do n=1,3
           tim(n,3) = tm(n)
         enddo
       endif
c
c     check if wall clock passed 24:00 (midnight)
c     since the last call to this routine
c
       if (tm(3).lt.tim(3,3)) tm(3) = tm(3) + 24*3600
c
c     get totals
c
       do n=1,3
         tim(n,2) = tm(n) - tim(n,3)
         tim(n,1) = tim(n,1) + tim(n,2)
         tim(n,3) = tm(n)
       enddo
       return
       end
