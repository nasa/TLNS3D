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
      subroutine bcpchget(imn,jmn,kmn,imxs,jmxs,kmxs,
     .                    imds,jmds,kmds,src,
     .                    nfaces,ndir1s,ndir2s,ndir3s,
     .                    nitems,nghost,wrk)
c=======================================================================
c 
c     object : to get ghost cell values and fill in work array
c              for patched block interfaces
c
c     initial coding :  May, 1995
c

      dimension src(imds,jmds,kmds)
      dimension ndir1s(1),ndir2s(1),ndir3s(1)
      dimension wrk(1)
c
c
c***********************************************************************
c     extract the information about ghost cells
c     adjacent to the patched interface boundary
c***********************************************************************
c
      n       = 0
      nitmbeg = 1
      nitmend = nitmbeg +nitems -1
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
        do 15 ighost=1,nghost
          do 10 nn=nitmbeg,nitmend
            j      = ndir2s(nn)
            k      = ndir3s(nn)
            n      = n +1
            wrk(n) = src(i,j,k)
  10      continue
          i      = i +ninc
  15    continue
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
        do 25 ighost=1,nghost
          do 20 nn=nitmbeg,nitmend
            k      = ndir3s(nn)
            i      = ndir1s(nn)
            n      = n +1
            wrk(n) = src(i,j,k)
  20      continue
          j      = j +ninc
  25    continue
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
        do 35 ighost=1,nghost
          do 30 nn=nitmbeg,nitmend
            i      = ndir1s(nn)
            j      = ndir2s(nn)
            n      = n +1
            wrk(n) = src(i,j,k)
  30      continue
          k      = k +ninc
  35    continue
c
      endif
c
      return
      end
