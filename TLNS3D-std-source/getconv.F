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
      subroutine getconv (mbloc,nbloc,wcount,rtrmsg,hrmsg,
     .                    nsupg,cxg,cyg,czg,cxvg,cyvg,czvg,
     .                    cmxg,cmyg,cmzg,
     .                    cx,cy,cz,cxv,cyv,czv,cmx,cmy,cmz,totw,wrk,
     .                    nodes )
c=======================================================================
c
c     get convergence parameters from all blocks in parallel code
c     Initial coding: by Bruce Wedan (1995)
c     Latest revision: by v.n.vatsa (Jan. 1998) 
c     (To account for effect of yaw angle on force computation)
c

#if defined(BUILD_PVM) || defined(BUILD_MPI)
#ifdef BUILD_PVM
#     include "fpvm3.h"
#  ifdef TLN_DBL
#    define RTYPE REAL8
#  else
#    define RTYPE REAL4
#  endif
#else
      include "mpif.h"
#  ifdef TLN_DBL
#    define RTYPE MPI_DOUBLE_PRECISION
#  else
#    define RTYPE MPI_REAL
#  endif
      dimension istat(MPI_STATUS_SIZE)
#endif
#     include "tags.h"
c
      character*128  errmsg
      common/dstcom/ errmsg,master,myrank,mycomm,ndlist,nnodes
      dimension      nodes(1)
c
      common/conv/  rtmaxg,hmaxg,rtrms2, hrms2, totpts,
     .              irtmxg,jrtmxg,krtmxg,mrtbloc,
     .               ihmxg, jhmxg, khmxg, mhbloc
c
      common/ma/    amachg
c
      common/rkerr/ erkmaxl(5),erkmaxg(5),erkrms2(5),erkrmsg(5),
     .              totptse
c
      dimension cx(mbloc),cy(mbloc),cmx(mbloc),cmy(mbloc),cmz(mbloc),
     .          cz(mbloc),cxv(mbloc),cyv(mbloc),czv(mbloc)
      dimension wrk(*)
c
      nconv = 25 + 9*nbloc
      nconv = nconv + 1
c
      if (myrank.eq.master) then
c
c     save stuff off from the master node before
c     destroying it all
c
        wrk(1)  = wcount - totw
        wrk(2)  = totpts
        wrk(3)  = rtrms2
        wrk(4)  = hrms2
        wrk(5)  = nsupg
        wrk(6)  = rtmaxg
        wrk(7)  = irtmxg
        wrk(8)  = jrtmxg
        wrk(9)  = krtmxg
        wrk(10) = mrtbloc
        wrk(11) = hmaxg
        wrk(12) = ihmxg
        wrk(13) = jhmxg
        wrk(14) = khmxg
        wrk(15) = mhbloc
        wrk(16) = amachg
c
        ncount  = 15
        ncount  = ncount + 1
        do nv   = 1,5
          ncount      = ncount + 1
          wrk(ncount) = erkrms2(nv)
        enddo
        do nv   = 1,5
          ncount      = ncount + 1
          wrk(ncount) = erkmaxg(nv)
        enddo
c
        n       = ncount
        do ibloc=1,nbloc
          wrk(n+1) = cx(ibloc)
          wrk(n+2) = cy(ibloc)
          wrk(n+3) = cz(ibloc)
          wrk(n+4) = cxv(ibloc)
          wrk(n+5) = cyv(ibloc)
          wrk(n+6) = czv(ibloc)
          wrk(n+7) = cmx(ibloc)
          wrk(n+8) = cmy(ibloc)
          wrk(n+9) = cmz(ibloc)
          n        = n +9
        enddo
c
c       initialize
c
        wcount = 0.
        totpts = 0.
        rtrms2 = 0.
        hrms2  = 0.
        nsupg  = 0
        rtmaxg = 0.
        hmaxg  = 0.
        amachg = 0.
        do nv = 1,5
           erkrms2(nv) = 0.
        enddo
        do nv = 1,5
           erkmaxg(nv) = 0.
        enddo
        do ibloc=1,nbloc
          cx(ibloc)  = 0.
          cy(ibloc)  = 0.
          cz(ibloc)  = 0.
          cxv(ibloc) = 0.
          cyv(ibloc) = 0.
          czv(ibloc) = 0.
          cmx(ibloc) = 0.
          cmy(ibloc) = 0.
          cmz(ibloc) = 0.
        enddo
c
c       collect results from the nodes
c
        do nn=1,nnodes
          inode  = nodes(ndlist+nn)
#ifdef BUILD_PVM
          call PVMFprecv (inode,TAG_CONV,
     .                    wrk,nconv,RTYPE,
     .                    itid,itag,ilen,ierr)
#else
cjm
      if (inode-1.ne.master) then
          call MPI_Recv (wrk,nconv,RTYPE,
     .                   inode-1,TAG_CONV,
     .                   mycomm,istat,ierr)
      endif
#endif
c
          wcount = wcount +wrk(1)
          totpts = totpts +wrk(2)
          rtrms2 = rtrms2 +wrk(3)
          hrms2  = hrms2  +wrk(4)
          nsupg  = nsupg  +wrk(5)
          do nv = 1,5
             erkrms2(nv) = erkrms2(nv) + wrk(nv+15)
          enddo
c
          do nv = 1,5
            if (erkmaxg(nv).lt.wrk(nv+20)) then
                erkmaxg(nv) =  wrk(nv+20)
            endif
          enddo
         
c
          if (rtmaxg.lt.wrk(6)) then
            rtmaxg  = wrk(6)
            irtmxg  = wrk(7)
            jrtmxg  = wrk(8)
            krtmxg  = wrk(9)
            mrtbloc = wrk(10)
            jbloc   = 0
            jcheck  = 0
            do ibloc=1,nbloc
              if (nodes(ibloc).eq.inode) then
                jbloc = jbloc + 1
                if (mrtbloc.eq.jbloc.and.jcheck.eq.0) then
                   mrtbloc = ibloc
                   jcheck  = 1
                endif
              endif
            enddo
          endif
c
          if (hmaxg.lt.wrk(11)) then
            hmaxg  = wrk(11)
            ihmxg  = wrk(12)
            jhmxg  = wrk(13)
            khmxg  = wrk(14)
            mhbloc = wrk(15)
            jbloc   = 0
            jcheck  = 0
            do ibloc=1,nbloc
              if (nodes(ibloc).eq.inode) then
                jbloc = jbloc + 1
                if (mhbloc.eq.jbloc.and.jcheck.eq.0) then
                   mhbloc = ibloc
                   jcheck = 1
                endif
              endif
            enddo
          endif
c
          if (amachg.lt.wrk(16)) then
            amachg  = wrk(16)
          endif
c
          n        = 25
          n        = n + 1
          do ibloc=1,nbloc
            cx(ibloc)  = cx(ibloc)  +wrk(n+1)
            cy(ibloc)  = cy(ibloc)  +wrk(n+2)
            cz(ibloc)  = cz(ibloc)  +wrk(n+3)
            cxv(ibloc) = cxv(ibloc) +wrk(n+4)
            cyv(ibloc) = cyv(ibloc) +wrk(n+5)
            czv(ibloc) = czv(ibloc) +wrk(n+6)
            cmx(ibloc) = cmx(ibloc) +wrk(n+7)
            cmy(ibloc) = cmy(ibloc) +wrk(n+8)
            cmz(ibloc) = cmz(ibloc) +wrk(n+9)
            n          = n +9
          enddo
        enddo
c
c       compute totals
c
        wcount = totw + wcount / real(nnodes)
        rtrmsg = sqrt (rtrms2 / totpts)
        hrmsg  = sqrt (hrms2  / totpts)
        do nv = 1,5
           erkrmsg (nv) = sqrt(erkrms2(nv) / totpts)
        enddo

        cxg    = 0.
        cyg    = 0.
        czg    = 0.
        cxvg   = 0.
        cyvg   = 0.
        czvg   = 0.
        cmxg   = 0.
        cmyg   = 0.
        cmzg   = 0.
        do ibloc=1,nbloc
          cxg    = cxg  +cx(ibloc)
          cyg    = cyg  +cy(ibloc)
          czg    = czg  +cz(ibloc)
          cxvg   = cxvg +cxv(ibloc)
          cyvg   = cyvg +cyv(ibloc)
          czvg   = czvg +czv(ibloc)
          cmxg   = cmxg +cmx(ibloc)
          cmyg   = cmyg +cmy(ibloc)
          cmzg   = cmzg +cmz(ibloc)
        enddo
c
      else
c
        wrk(1)  = wcount - totw
        wrk(2)  = totpts
        wrk(3)  = rtrms2
        wrk(4)  = hrms2
        wrk(5)  = nsupg
        wrk(6)  = rtmaxg
        wrk(7)  = irtmxg
        wrk(8)  = jrtmxg
        wrk(9)  = krtmxg
        wrk(10) = mrtbloc
        wrk(11) = hmaxg
        wrk(12) = ihmxg
        wrk(13) = jhmxg
        wrk(14) = khmxg
        wrk(15) = mhbloc
        wrk(16) = amachg
c
        ncount  = 15
        ncount  = ncount + 1
        do nv   = 1,5
          ncount      = ncount + 1
          wrk(ncount) = erkrms2(nv)
        enddo
        do nv   = 1,5
          ncount      = ncount + 1
          wrk(ncount) = erkmaxg(nv)
        enddo
c
        n       = ncount
c
        do ibloc=1,nbloc
          wrk(n+1) = cx(ibloc)
          wrk(n+2) = cy(ibloc)
          wrk(n+3) = cz(ibloc)
          wrk(n+4) = cxv(ibloc)
          wrk(n+5) = cyv(ibloc)
          wrk(n+6) = czv(ibloc)
          wrk(n+7) = cmx(ibloc)
          wrk(n+8) = cmy(ibloc)
          wrk(n+9) = cmz(ibloc)
          n        = n +9
        enddo

#ifdef BUILD_PVM
        call PVMFpsend (master,TAG_CONV,
     .                  wrk,nconv,RTYPE,ierr)
#else
        call MPI_Send (wrk,nconv,RTYPE,master,
     .                 TAG_CONV,mycomm,ierr)
#endif
c
      endif
c
#endif
c
# ifdef BUILD_MPI
      call MPI_Bcast(rtmaxg,1,RTYPE,master,mycomm,ierr)
# endif
      return
      end
