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
c=======================================================================
      subroutine ERREXIT (nodes)
c=======================================================================

      dimension      nodes(1)

#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c     kill all processes and exit
c
#ifdef BUILD_PVM
#     include "fpvm3.h"
#else
      include "mpif.h"
#endif
#     include "tags.h"
c
      character*128  errmsg
      common/dstcom/ errmsg,master,myrank,mycomm,ndlist,nnodes
c
      if (myrank.eq.master) then
        do inode=1,nnodes
#ifdef BUILD_PVM
          call PVMFpsend (nodes(ndlist+inode),TAG_ERROR,
     .                    ndum,1,INTEGER4,ierr)
#else
          if (nodes(ndlist+inode)-1.ne.master) then
          call MPI_Send (ndum,1,MPI_INTEGER,
cjm     .                   nodes(ndlist+inode),TAG_ERROR,
     .                   nodes(ndlist+inode)-1,TAG_ERROR,
     .                   mycomm,ierr)
          end if
#endif
        enddo
      endif
c
#ifdef BUILD_PVM
      call PVMFexit (ierr)
#else
      call MPI_Finalize (ierr)
#endif
#endif
c
      stop
      end
