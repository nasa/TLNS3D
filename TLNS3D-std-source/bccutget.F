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
      subroutine bccutget(imn,jmn,kmn,imxs,jmxs,kmxs,
     .                    imds,jmds,kmds,src,
     .                    nfaces,n1begs,n1ends,n2begs,n2ends,
     .                    nghost,wrk                         )
c=======================================================================
c 
c     object : to get ghost cell values and fill in work array
c              at block interfaces and interior cuts
c
c     initial coding : April, 1995
c

      dimension src(imds,jmds,kmds)
      dimension wrk(1)
c
c
c***********************************************************************
c     extract the information about ghost cells
c     adjacent to the block/cut interface boundary
c***********************************************************************
c
c     find the increment for directions 1 and 2 of source segment
c
      if (n1begs .le. n1ends) then
        n1dirs =  1
        n1strt = n1begs
        n1stop = n1ends +1
      else
        n1dirs = -1
        n1strt = n1begs +1
        n1stop = n1ends
      endif
c
      if (n2begs .le. n2ends) then
        n2dirs =  1
        n2strt = n2begs
        n2stop = n2ends +1
      else
        n2dirs = -1
        n2strt = n2begs +1
        n2stop = n2ends
      endif
c
c     index counter
c
      n = 0
c
c     store the required variables from source segment into
c     temporary working arrays
c     before doing so, test for the face type of the source segment
c     and set the prper indices for the face in question
c
      if (abs(nfaces).le.2) then
c
c---------------------- i = constant boundary  ------------------
c
        if (abs(nfaces).eq.1) then
          i    = imn + 1
          ninc = 1
        else
          i    = imxs
          ninc = -1
        endif
c
        do 19 ighost = 1, nghost
c
c        direction 1 of source maps to direction 1 of target
c
         if (nfaces.gt.0) then
           do 10 k = n2strt, n2stop, n2dirs
           do 10 j = n1strt, n1stop, n1dirs
           n      = n + 1
           wrk(n) = src(i,j,k)
  10       continue
c
c        direction 1 of source maps to direction 2 of target
c
         else
           do 15 j = n1strt, n1stop, n1dirs
           do 15 k = n2strt, n2stop, n2dirs
           n      = n + 1
           wrk(n) = src(i,j,k)
  15       continue
         endif
c
         i      = i + ninc
  19     continue
c
      else if (abs(nfaces).le.4) then
c
c---------------------- j = constant boundary  ------------------
c
        if (abs(nfaces).eq.3) then
          j    = jmn + 1
          ninc = 1
        else
          j    = jmxs
          ninc = -1
        endif
c
        do 29 ighost = 1, nghost
c
c        direction 1 of source maps to direction 1 of target
c
         if (nfaces.gt.0) then
           do 20 i = n2strt, n2stop, n2dirs
           do 20 k = n1strt, n1stop, n1dirs
           n      = n + 1
           wrk(n) = src(i,j,k)
  20       continue
c
c        direction 1 of source maps to direction 2 of target
c
         else
           do 25 k = n1strt, n1stop, n1dirs
           do 25 i = n2strt, n2stop, n2dirs
           n      = n + 1
           wrk(n) = src(i,j,k)
  25       continue
         endif
c
         j      = j + ninc
  29     continue
c
      else
c
c---------------------- k = constant boundary  ------------------
c
        if (abs(nfaces).eq.5) then
          k    = kmn + 1
          ninc = 1
        else
          k    = kmxs
          ninc = -1
        endif
c
        do 39 ighost = 1, nghost
c
c        direction 1 of source maps to direction 1 of target
c
         if (nfaces.gt.0) then
           do 30 j = n2strt, n2stop, n2dirs
           do 30 i = n1strt, n1stop, n1dirs
           n      = n + 1
           wrk(n) = src(i,j,k)
  30       continue
c
c        direction 1 of source maps to direction 2 of target
c
         else
           do 35 i = n1strt, n1stop, n1dirs
           do 35 j = n2strt, n2stop, n2dirs
           n      = n + 1
           wrk(n) = src(i,j,k)
  35       continue
         endif
c
         k      = k + ninc
  39     continue
c
      endif
c
      return
      end
