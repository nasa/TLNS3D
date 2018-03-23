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
      subroutine bcpchset(imn,jmn,kmn,imx,jmx,kmx,
     .                    imd,jmd,kmd,dst,
     .                    nface,ndir1,ndir2,ndir3,
     .                    nitems,nghost,frc,wrk)
c=======================================================================
c 
c     object : to set the ghost cell values from the work array
c              for pathed block interfaces
c
c     initial coding : April, 1995
c

      dimension dst(imd,jmd,kmd)
      dimension ndir1(1),ndir2(1),ndir3(1)
      dimension frc(1),wrk(1)
c
      n       = 0
      nitmbeg = 1
      nitmend = nitmbeg +nitems -1
c
c**********************************************************************
c     transfering the ghost cell values to target segment
c**********************************************************************
c
c     test for the face type of the target segment and
c     set the prper indices for the face in question
c
      if (nface.le.2) then
c
c---------------------- i = constant boundary  ------------------
c
        if (nface .eq. 1) then
          i    = imn
          ninc = -1
        else
          i    = imx + 1
          ninc = 1
        endif
c
        do 41 ighost=1,nghost
          do 40 nn=nitmbeg,nitmend
            j          = ndir2(nn)
            k          = ndir3(nn)
            n          = n +1
            dst(i,j,k) = dst(i,j,k) + frc(nn)*wrk(n)
  40      continue
          i          = i + ninc
  41    continue
c 
      else if (nface .le. 4) then
c
c---------------------- j = constant boundary  ------------------
c
        if (nface .eq. 3) then
          j    = jmn
          ninc = -1
        else
          j    = jmx + 1
          ninc = 1
        endif
c
        do 51 ighost=1,nghost
          do 50 nn=nitmbeg,nitmend
            k          = ndir3(nn)
            i          = ndir1(nn)
            n          = n +1
            dst(i,j,k) = dst(i,j,k) + frc(nn)*wrk(n)
  50      continue
          j          = j + ninc
  51    continue
c 
      else
c
c---------------------- k = constant boundary  ------------------
c
        if (nface .eq. 5) then
          k    = kmn
          ninc = -1
        else
          k    = kmx + 1
          ninc = 1
        endif
c
        do 61 ighost=1,nghost
          do 60 nn=nitmbeg,nitmend
            i          = ndir1(nn)
            j          = ndir2(nn)
            n          = n +1
            dst(i,j,k) = dst(i,j,k) + frc(nn)*wrk(n)
  60      continue
          k          = k + ninc
  61    continue
c
      endif
c
      return
      end
