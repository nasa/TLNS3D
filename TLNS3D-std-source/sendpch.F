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
c====================================================================
      subroutine sendpch (mgrlev,mbloc,ngrid,nbloc,
     .         mxtpchc,ntpchcb,npchcbf,ipatchc,jpatchc,kpatchc,
     .         mxpchs,mpchitm,nswpchb,npchitm,lspchb1,lspchf1,
     .         lspchb2,lspchf2,ipitmb1,ipitmb2,jpitmb1,jpitmb2,
     .         kpitmb1,kpitmb2,frc,m1pch1,lpchitm,lswpchb,lswpche ,
     .         nodes )
c====================================================================
c
c     called by PVM or MPI to send/receive patched grid data
c

#if defined(BUILD_PVM) || defined(BUILD_MPI)
#  ifdef BUILD_PVM
#     include "fpvm3.h"
#  ifdef TLN_DBL
#    define RTYPE REAL8
#  else
#    define RTYPE REAL4
#  endif
#  else
      include "mpif.h"
#  ifdef TLN_DBL
#    define RTYPE MPI_DOUBLE_PRECISION
#  else
#    define RTYPE MPI_REAL
#  endif
      dimension istat(MPI_STATUS_SIZE)
#  endif
#     include "tags.h"
c
      character*128  errmsg
      common/dstcom/ errmsg,master,myrank,mycomm,ndlist,nnodes
      dimension      nodes(1)
c
      dimension   ntpchcb(mbloc,mgrlev), npchcbf(6,mbloc,mgrlev),
     .            m1pch1 (mbloc,mgrlev),
     .          ipatchc(mxtpchc), jpatchc(mxtpchc),
     .          kpatchc(mxtpchc)
c    .          ipatchc(mxtpchc,mgrlev), jpatchc(mxtpchc,mgrlev),
c    .          kpatchc(mxtpchc,mgrlev)
c
      dimension nswpchb  (mbloc,mgrlev),   npchitm(mxpchs,mgrlev),
     .          lswpchb  (mbloc,mgrlev),   lpchitm(mxpchs,mgrlev),
     .          lswpche  (mbloc,mgrlev),
     .          lspchb1 (mxpchs,mgrlev),   lspchf1(mxpchs,mgrlev),
     .          lspchb2 (mxpchs,mgrlev),   lspchf2(mxpchs,mgrlev),
     .          ipitmb1(mpchitm),   ipitmb2(mpchitm),
     .          jpitmb1(mpchitm),   jpitmb2(mpchitm),
     .          kpitmb1(mpchitm),   kpitmb2(mpchitm),
     .          frc    (mpchitm)
c    .          ipitmb1(mpchitm,mgrlev),   ipitmb2(mpchitm,mgrlev),
c    .          jpitmb1(mpchitm,mgrlev),   jpitmb2(mpchitm,mgrlev),
c    .          kpitmb1(mpchitm,mgrlev),   kpitmb2(mpchitm,mgrlev),
c    .          frc   (mpchitm,mgrlev)
c
      if (myrank.eq.master) then
c
#  ifdef BUILD_PVM
        call PVMFinitsend (PvmDataInPlace,ierr)
        call PVMFpack (INTEGER4,ntpchcb,mbloc*ngrid,1,ierr)
        call PVMFmcast (nnodes,nodes(ndlist+1),TAG_PATCH,ierr)
c
        call PVMFinitsend (PvmDataInPlace,ierr)
        call PVMFpack (INTEGER4,npchcbf,6*mbloc*ngrid,1,ierr)
        call PVMFpack (INTEGER4,nswpchb,mbloc*ngrid,1,ierr)
c
        call PVMFpack (INTEGER4,ipatchc,mxtpchc*ngrid,1,ierr)
        call PVMFpack (INTEGER4,jpatchc,mxtpchc*ngrid,1,ierr)
        call PVMFpack (INTEGER4,kpatchc,mxtpchc*ngrid,1,ierr)
c
        call PVMFpack (INTEGER4,npchitm,mxpchs*ngrid,1,ierr)
        call PVMFpack (INTEGER4,lspchb1,mxpchs*ngrid,1,ierr)
        call PVMFpack (INTEGER4,lspchf1,mxpchs*ngrid,1,ierr)
        call PVMFpack (INTEGER4,lspchb2,mxpchs*ngrid,1,ierr)
        call PVMFpack (INTEGER4,lspchf2,mxpchs*ngrid,1,ierr)
c
        call PVMFpack (INTEGER4,ipitmb1,mpchitm*ngrid,1,ierr)
        call PVMFpack (INTEGER4,ipitmb2,mpchitm*ngrid,1,ierr)
        call PVMFpack (INTEGER4,jpitmb1,mpchitm*ngrid,1,ierr)
        call PVMFpack (INTEGER4,jpitmb2,mpchitm*ngrid,1,ierr)
        call PVMFpack (INTEGER4,kpitmb1,mpchitm*ngrid,1,ierr)
        call PVMFpack (INTEGER4,kpitmb2,mpchitm*ngrid,1,ierr)
c
        call PVMFpack (RTYPE,frc,mpchitm*ngrid,1,ierr)
        call PVMFmcast (nnodes,nodes(ndlist+1),TAG_PATCH,ierr)
#  else
        do inode=1,nnodes
          call MPI_Send (ntpchcb,mbloc*ngrid,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
c
          call MPI_Send (npchcbf,6*mbloc*ngrid,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (nswpchb,mbloc*ngrid,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
c
cnewpch   call MPI_Send (ipatchc,mxtpchc*ngrid,MPI_INTEGER,
          call MPI_Send (ipatchc,mxtpchc,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
cnewpch   call MPI_Send (jpatchc,mxtpchc*ngrid,MPI_INTEGER,
          call MPI_Send (jpatchc,mxtpchc,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
cnewpch   call MPI_Send (kpatchc,mxtpchc*ngrid,MPI_INTEGER,
          call MPI_Send (kpatchc,mxtpchc,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
c
          call MPI_Send (npchitm,mxpchs*ngrid,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (lspchb1,mxpchs*ngrid,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (lspchf1,mxpchs*ngrid,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (lspchb2,mxpchs*ngrid,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (lspchf2,mxpchs*ngrid,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
c
cnewpch   call MPI_Send (ipitmb1,mpchitm*ngrid,MPI_INTEGER,
          call MPI_Send (ipitmb1,mpchitm,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
cnewpch   call MPI_Send (ipitmb2,mpchitm*ngrid,MPI_INTEGER,
          call MPI_Send (ipitmb2,mpchitm,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
cnewpch   call MPI_Send (jpitmb1,mpchitm*ngrid,MPI_INTEGER,
          call MPI_Send (jpitmb1,mpchitm,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
cnewpch   call MPI_Send (jpitmb2,mpchitm*ngrid,MPI_INTEGER,
          call MPI_Send (jpitmb2,mpchitm,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
cnewpch   call MPI_Send (kpitmb1,mpchitm*ngrid,MPI_INTEGER,
          call MPI_Send (kpitmb1,mpchitm,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
cnewpch   call MPI_Send (kpitmb2,mpchitm*ngrid,MPI_INTEGER,
          call MPI_Send (kpitmb2,mpchitm,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
c
cnewpch   call MPI_Send (frc,mpchitm*ngrid,RTYPE,
          call MPI_Send (frc,mpchitm,RTYPE,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
c
cnewpch (vatsa: Feb. 97) additional pointers needed for new patcher
c
          call MPI_Send (m1pch1,mbloc*ngrid,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
c
          call MPI_Send (lswpchb,mbloc*ngrid,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
          call MPI_Send (lswpche,mbloc*ngrid,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
c
          call MPI_Send (lpchitm,mxpchs*ngrid,MPI_INTEGER,
     .                   nodes(ndlist+inode),TAG_PATCH,
     .                   mycomm,ierr)
        enddo
#  endif
c
      else
c
# ifdef BUILD_PVM
        call PVMFprecv (master,-1,ntpchcb,mbloc*ngrid,INTEGER4,
     .                  itid,itag,ilen,ierr)
# else
        call MPI_Recv (ntpchcb,mbloc*ngrid,MPI_INTEGER,master,
     .                 MPI_ANY_TAG,mycomm,istat,ierr)
        itag = istat(MPI_TAG)
# endif
        if (itag.ne.TAG_PATCH)
     .    call ERREXIT (nodes)
c
# ifdef BUILD_PVM
        call PVMFrecv (master,TAG_PATCH,ierr)
        call PVMFunpack (INTEGER4,npchcbf,6*mbloc*ngrid,1,ierr)
        call PVMFunpack (INTEGER4,nswpchb,mbloc*ngrid,1,ierr)
c
        call PVMFunpack (INTEGER4,ipatchc,mxtpchc*ngrid,1,ierr)
        call PVMFunpack (INTEGER4,jpatchc,mxtpchc*ngrid,1,ierr)
        call PVMFunpack (INTEGER4,kpatchc,mxtpchc*ngrid,1,ierr)
c
        call PVMFunpack (INTEGER4,npchitm,mxpchs*ngrid,1,ierr)
        call PVMFunpack (INTEGER4,lspchb1,mxpchs*ngrid,1,ierr)
        call PVMFunpack (INTEGER4,lspchf1,mxpchs*ngrid,1,ierr)
        call PVMFunpack (INTEGER4,lspchb2,mxpchs*ngrid,1,ierr)
        call PVMFunpack (INTEGER4,lspchf2,mxpchs*ngrid,1,ierr)
c
        call PVMFunpack (INTEGER4,ipitmb1,mpchitm*ngrid,1,ierr)
        call PVMFunpack (INTEGER4,ipitmb2,mpchitm*ngrid,1,ierr)
        call PVMFunpack (INTEGER4,jpitmb1,mpchitm*ngrid,1,ierr)
        call PVMFunpack (INTEGER4,jpitmb2,mpchitm*ngrid,1,ierr)
        call PVMFunpack (INTEGER4,kpitmb1,mpchitm*ngrid,1,ierr)
        call PVMFunpack (INTEGER4,kpitmb2,mpchitm*ngrid,1,ierr)
c
        call PVMFunpack (RTYPE,frc,mpchitm*ngrid,1,ierr)
# else
        call MPI_Recv (npchcbf,6*mbloc*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
        call MPI_Recv (nswpchb,mbloc*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
c
cnewpch call MPI_Recv (ipatchc,mxtpchc*ngrid,MPI_INTEGER,
        call MPI_Recv (ipatchc,mxtpchc,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
cnewpch call MPI_Recv (jpatchc,mxtpchc*ngrid,MPI_INTEGER,
        call MPI_Recv (jpatchc,mxtpchc,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
cnewpch call MPI_Recv (kpatchc,mxtpchc*ngrid,MPI_INTEGER,
        call MPI_Recv (kpatchc,mxtpchc,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
c
        call MPI_Recv (npchitm,mxpchs*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
        call MPI_Recv (lspchb1,mxpchs*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
        call MPI_Recv (lspchf1,mxpchs*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
        call MPI_Recv (lspchb2,mxpchs*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
        call MPI_Recv (lspchf2,mxpchs*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
c
cnewpch call MPI_Recv (ipitmb1,mpchitm*ngrid,MPI_INTEGER,
        call MPI_Recv (ipitmb1,mpchitm,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
cnewpch call MPI_Recv (ipitmb2,mpchitm*ngrid,MPI_INTEGER,
        call MPI_Recv (ipitmb2,mpchitm,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
cnewpch call MPI_Recv (jpitmb1,mpchitm*ngrid,MPI_INTEGER,
        call MPI_Recv (jpitmb1,mpchitm,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
cnewpch call MPI_Recv (jpitmb2,mpchitm*ngrid,MPI_INTEGER,
        call MPI_Recv (jpitmb2,mpchitm,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
cnewpch call MPI_Recv (kpitmb1,mpchitm*ngrid,MPI_INTEGER,
        call MPI_Recv (kpitmb1,mpchitm,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
cnewpch call MPI_Recv (kpitmb2,mpchitm*ngrid,MPI_INTEGER,
        call MPI_Recv (kpitmb2,mpchitm,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
c
cnewpch call MPI_Recv (frc,mpchitm*ngrid,RTYPE,
        call MPI_Recv (frc,mpchitm,RTYPE,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
c
cnewpch (vatsa: Feb. 97) additional pointers needed for new patcher
c
        call MPI_Recv (m1pch1,mbloc*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
c
        call MPI_Recv (lswpchb,mbloc*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
c
        call MPI_Recv (lswpche,mbloc*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
c
        call MPI_Recv (lpchitm,mxpchs*ngrid,MPI_INTEGER,
     .                 master,TAG_PATCH,mycomm,istat,ierr)
c
c
# endif
c
      endif
c
#endif
      return
      end
