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
      subroutine sort_x(nsurf,surf,ntri,iptr,wk3d5,iwrk)
      dimension surf(nsurf,3),ntri(nsurf),iptr(nsurf,8)
      dimension wk3d5(*), iwrk(*)
      integer temp
c      pointer (p_temp,temp(nsurf)),(p_itemp,itemp(nsurf))
c      pointer (p_iperm,iperm(nsurf)),(p_inv_iperm,inv_iperm(nsurf))
c      p_iperm=my_malloc(nsurf)
c      p_temp=my_malloc(nsurf)
      iperm = iialloc(nsurf)
      temp = ifalloc(nsurf)
      call heap_sort(nsurf,surf,iwrk(iperm))
      do j=1,3
         call move_real(nsurf,surf(1,j),iwrk(iperm),wk3d5(temp))
      end do
c      call my_free(p_temp,nsurf)
      call ffree(nsurf)
c      p_itemp=my_malloc(nsurf)
      itemp = iialloc(nsurf)
      call move_integer(nsurf,ntri,iwrk(iperm),iwrk(itemp))
      do j=1,8
         call move_integer(nsurf,iptr(1,j),iwrk(iperm),iwrk(itemp))
      end do
c      call my_free(p_itemp,nsurf)
      call ifree(nsurf)
c      p_inv_iperm=my_malloc(nsurf)
      inv_iperm = iialloc(nsurf)
      do i=1,nsurf
	 ii = iwrk(iperm+i-1)
         iwrk(inv_iperm+ii-1)=i
      end do
c      call my_free(p_iperm,nsurf)
      do i=1,nsurf
         numtri=ntri(i)
         do j=1,2*numtri
	    ii = iptr(i,j)
            iptr(i,j)=iwrk(inv_iperm+ii-1)
         end do
      end do
c      call my_free(p_inv_iperm,nsurf)
      call ifree(nsurf)
      call ifree(nsurf)
      return
      end
