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
      subroutine bccutset(imn,jmn,kmn,imx,jmx,kmx,
     .                    imd,jmd,kmd,dst,
     .                    nface,m1beg,m1end,m2beg,m2end,
     .                    nghost,wrk,n1cnt,n2cnt        )
c=======================================================================
c
c     object : to set the ghost cell values from the work array
c              at block interfaces and interior cuts
c
c     initial coding : April, 1995
c

      dimension dst(imd,jmd,kmd)
      dimension wrk(n1cnt,n2cnt,nghost)
c
      n1beg = m1beg+1
      n1end = m1end
      n2beg = m2beg+1
      n2end = m2end
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
c       if (n1beg .ne. jmn) n1beg = n1beg +1
c       if (n1end .eq. jmx) n1end = n1end +1
c       if (n2beg .ne. kmn) n2beg = n2beg +1
c       if (n2end .eq. kmx) n2end = n2end +1
c
        if (nface .eq. 1) then
          i    = imn
          ninc = -1
        else
          i    = imx + 1
          ninc = 1
        endif
c
        n1b = 0
        if (n1beg .ne. jmn) n1b = 1
        n2b = 0
        if (n2beg .ne. kmn) n2b = 1
c
        do 41 ighost = 1, nghost
        n2         = n2b
        do 40 k = n2beg, n2end
        n2         = n2 + 1
        n1         = n1b
        do 40 j = n1beg, n1end
        n1         = n1 + 1
        dst(i,j,k) = wrk(n1,n2,ighost)
  40    continue
c
        i          = i + ninc
  41    continue
c
      else if (nface .le. 4) then
c
c---------------------- j = constant boundary  ------------------
c
c       if (n1beg .ne. kmn) n1beg = n1beg +1
c       if (n1end .eq. kmx) n1end = n1end +1
c       if (n2beg .ne. imn) n2beg = n2beg +1
c       if (n2end .eq. imx) n2end = n2end +1
c
        if (nface .eq. 3) then
          j    = jmn
          ninc = -1
        else
          j    = jmx + 1
          ninc = 1
        endif
c
        n1b = 0
        if (n1beg .ne. kmn) n1b = 1
        n2b = 0
        if (n2beg .ne. imn) n2b = 1
c
        do 51 ighost = 1, nghost
        n2         = n2b
        do 50 i = n2beg, n2end
        n2         = n2 + 1
        n1         = n1b
        do 50 k = n1beg, n1end
        n1         = n1 + 1
        dst(i,j,k) = wrk(n1,n2,ighost)
  50    continue
c
        j          = j + ninc
  51    continue
c
      else
c
c---------------------- k = constant boundary  ------------------
c
c       if (n1beg .ne. imn) n1beg = n1beg +1
c       if (n1end .eq. imx) n1end = n1end +1
c       if (n2beg .ne. jmn) n2beg = n2beg +1
c       if (n2end .eq. jmx) n2end = n2end +1
c
        if (nface .eq. 5) then
          k    = kmn
          ninc = -1
        else
          k    = kmx + 1
          ninc = 1
        endif
c
        n1b = 0
        if (n1beg .ne. imn) n1b = 1
        n2b = 0
        if (n2beg .ne. jmn) n2b = 1
c
        do 61 ighost = 1, nghost
        n2         = n2b
        do 60 j = n2beg, n2end
        n2         = n2 + 1
        n1         = n1b
        do 60 i = n1beg, n1end
        n1         = n1 + 1
        dst(i,j,k) = wrk(n1,n2,ighost)
  60    continue
c
        k          = k + ninc
  61    continue
c
      endif
c
      return
      end
