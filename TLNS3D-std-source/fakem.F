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
      subroutine initf(isize)
      common /alloc/ iptr,imax,ifptr,ifmax
c
      ifptr = 1
      ifmax = isize
c      print *,'ifmax = ',ifmax
c
      return
      end

      subroutine initi(isize)
      common /alloc/ iptr,imax,ifptr,ifmax
c
      iptr = 1
      imax = isize
c      print *,'imax = ',imax
c
      return
      end

      subroutine ifree(isize)
      common /alloc/ iptr,imax,ifptr,ifmax
c
      if (isize.ge.iptr) then
	 print *,'freeing too much: ',isize,iptr
	 stop
      endif
      iptr = iptr - isize
c      print *,'ifree: iptr = ',iptr
c
      return
      end

      subroutine ffree(isize)
      common /alloc/ iptr,imax,ifptr,ifmax
c
      if (isize.ge.ifptr) then
	 print *,'freeing too much: ',isize,ifptr
	 stop
      endif
      ifptr = ifptr - isize
c      print *,'ffree: ifptr = ',ifptr
c
      return
      end

      integer function iialloc(isize)
      common /alloc/ iptr,imax,ifptr,ifmax
c
      iialloc = iptr
      iptr = iptr + isize
      if (iptr.gt.(imax+1)) then
	 print *,'iialloc failed: ',isize,(iptr-isize)
	 stop
      endif
c      print *,'iialloc: iptr0,size,iptr = ',iialloc,
c     .      isize,iptr
c
      return
      end

      integer function ifalloc(isize)
      common /alloc/ iptr,imax,ifptr,ifmax
c
      ifalloc = ifptr
      ifptr = ifptr + isize
      if (ifptr.gt.(ifmax+1)) then
	 print *,'ifalloc failed: ',isize,(ifptr-isize)
	 stop
      endif
c      print *,'ifalloc: ifptr0,size,ifptr = ',ifalloc,
c     .      isize,ifptr
c
      return
      end

      integer function isrcheq(n,x,incx,target)
      dimension x(*)
c
      j = 1
      isrcheq = 0
      if (n.le.0) return
      if (incx.lt.0) j = 1 - (n-1)*incx
      do 100 i=1,n
        if (x(j).eq.target) goto 200
	j = j + incx
 100  continue
 200  isrcheq = i
c
      return
      end

