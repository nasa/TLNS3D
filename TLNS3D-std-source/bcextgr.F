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
      subroutine  bcextgr(imn,jmn,kmn,
     .                    imx ,jmx ,kmx ,imp2 ,jmp2 ,kmp2 , vol,
     .                    nface ,n1beg ,n1end ,n2beg ,n2end ,isoln )
c=======================================================================
c 
c     object : to compute the ghost cell values for cell-volumes
c              by simple extrpolation
c
c     initial coding : by v.n.vatsa, nasa langley ( March, 1990)
c
c
c

      dimension      vol(imp2,jmp2,kmp2)
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c
c
c      branch-cut for different type of boundary faces
c
      go to (100,200,300,400,500,600) nface
c
c
  100 continue
c--------------------- i = imn  face ----------------------------
c
      if( n1beg .ne. jmn) n1beg = n1beg +1
      if( n1end .eq. jmx) n1end = n1end +1
      if( n2beg .ne. kmn) n2beg = n2beg +1
      if( n2end .eq. kmx) n2end = n2end +1
c
      i   =   imn
c
      do 110 k= n2beg, n2end
      do 110 j= n1beg, n1end
      vol(i,j,k)   =    vol(i+1,j,k)
  110 continue
c
      return
c
c
  200 continue
c--------------------- i = imx  face ----------------------------
c
      if( n1beg .ne. jmn) n1beg = n1beg +1
      if( n1end .eq. jmx) n1end = n1end +1
      if( n2beg .ne. kmn) n2beg = n2beg +1
      if( n2end .eq. kmx) n2end = n2end +1
c
      i   =   imx + 1
c
      do 210 k= n2beg, n2end
      do 210 j= n1beg, n1end
      vol(i,j,k)   =    vol(i-1,j,k)
  210 continue
c
      return
c
c
  300 continue
c--------------------- j = jmn  face ----------------------------
c
      if( n1beg .ne. kmn) n1beg = n1beg +1
      if( n1end .eq. kmx) n1end = n1end +1
      if( n2beg .ne. imn) n2beg = n2beg +1
      if( n2end .eq. imx) n2end = n2end +1
c
      j   =   jmn
c
      do 310 k= n1beg, n1end
      do 310 i= n2beg, n2end
      vol(i,j,k)   =    vol(i,j+1,k)
  310 continue
c
      return
c
c
  400 continue
c--------------------- j = jmx  face ----------------------------
c
      if( n1beg .ne. kmn) n1beg = n1beg +1
      if( n1end .eq. kmx) n1end = n1end +1
      if( n2beg .ne. imn) n2beg = n2beg +1
      if( n2end .eq. imx) n2end = n2end +1
c
      j   =   jmx + 1
c
      do 410 k= n1beg, n1end
      do 410 i= n2beg, n2end
      vol(i,j,k)   =    vol(i,j-1,k)
  410 continue
c
      return
c
c
  500 continue
c--------------------- k = kmn  face ----------------------------
c
      if( n1beg .ne. imn) n1beg = n1beg +1
      if( n1end .eq. imx) n1end = n1end +1
      if( n2beg .ne. jmn) n2beg = n2beg +1
      if( n2end .eq. jmx) n2end = n2end +1
c
      k   =   kmn
c
      do 510 j= n2beg, n2end
      do 515 i= n1beg, n1end
      vol(i,j,k)   =    vol(i,j,k+1)
  515 continue
      if(i2dfl.eq.1) then
        do 520 i= n1beg, n1end
        vol(i,j,k-1)   =    vol(i,j,k+1)
  520   continue
      endif
c
  510 continue
c
      return
c
c
  600 continue
c--------------------- k = kmx  face ----------------------------
c
      if( n1beg .ne. imn) n1beg = n1beg +1
      if( n1end .eq. imx) n1end = n1end +1
      if( n2beg .ne. jmn) n2beg = n2beg +1
      if( n2end .eq. jmx) n2end = n2end +1
c
      k   =   kmx + 1
c
      do 610 j= n2beg, n2end
      do 615 i= n1beg, n1end
      vol(i,j,k)   =    vol(i,j,k-1)
  615 continue
      if(i2dfl.eq.1) then
        do 620 i= n1beg, n1end
        vol(i,j,k+1)   =    vol(i,j,k-1)
  620   continue
      endif
c
  610 continue
c
c
c
      return
      end
