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
      subroutine  savdata(ibeg,iend,iskip,jbeg,jend,jskip,kbeg,kend,
     .                    kskip,
     .                    iblock,ngrp,
     .                    imp1,jmp1,kmp1,
     .                    imp2,jmp2,kmp2,
     .                    mxsizc,mxsiz5c,mxsiz3n,
     .                    mgrlev,mbloc,msegt,msegn,
     .                    m5cc,m1cc,m3np,
     .                    imap,nbloc,nseg,
# if defined BUILD_MPI
     .                   mxszfg,buffw,buffp,buff,
# endif
     .                    x,w,p,wn,alpha,totptim,ifirst,nodes       )
c=======================================================================
c
c     object: To save solutions at selected points at all time-steps
c     original coding : M.D. Sanetrik
c     MPI (distributed) version : By V.N.Vatsa (May, 2000)
# if defined BUILD_MPI
      include "mpif.h"
#     include "tags.h"
#  ifdef TLN_DBL
#    define RTYPE MPI_DOUBLE_PRECISION
#  else
#    define RTYPE MPI_REAL
#  endif
      dimension istat(MPI_STATUS_SIZE)
#endif
      character*128  errmsg
      common/dstcom/ errmsg,master,myrank,mycomm,ndlist,nnodes
      dimension      nodes(1)
c
      dimension   w(mxsiz5c),p(mxsizc),x(mxsiz3n)
c
      dimension  wn(1)
c
      dimension                      imp1(mgrlev,mbloc),
     .           jmp1(mgrlev,mbloc), kmp1(mgrlev,mbloc),
     .           imp2(mgrlev,mbloc), jmp2(mgrlev,mbloc),
     .           kmp2(mgrlev,mbloc)
c
      dimension  m1cc(mgrlev,mbloc), m5cc(mgrlev,mbloc),
     .           m3np(mgrlev,mbloc)
c
      dimension  imap(msegt,msegn,mbloc), nseg(mbloc)
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
c
# if defined BUILD_MPI
      dimension buffw(5*mxszfg),buffp(mxszfg),buff(mxsiz3n)
# endif
c
c     common block containing limits for groups of points to be output
c     at every time step for time-accurate runs
c
      dimension     ibeg(ngrp),iend(ngrp),iskip(ngrp),
     .              jbeg(ngrp),jend(ngrp),jskip(ngrp),
     .              kbeg(ngrp),kend(ngrp),kskip(ngrp),
     .              iblock(ngrp)
c
      do 120 igrp = 1,ngrp
      ibloc = iblock(igrp)
c
c     write out flow solutions/grid at selected blocks/points
c
c     For distributed version, we first need to collect the
c     flow variables from different nodes
c
      npts  = imp2(1,ibloc)*jmp2(1,ibloc)*kmp2(1,ibloc)
      nptsg = imp1(1,ibloc)*jmp1(1,ibloc)*kmp1(1,ibloc)*3
c
# if defined BUILD_MPI
c     receive the flow data at the master node
      if (myrank.eq.master) then
        if (nodes(ibloc)-1.ne.myrank) then
c       data is nonlocal
          call MPI_Recv (buffw,5*npts,RTYPE,
     .                   nodes(ibloc)-1,TAG_SAVE,mycomm,istat,ierr)
          call MPI_Recv (buffp,npts,RTYPE,
     .                   nodes(ibloc)-1,TAG_SAVE,mycomm,istat,ierr)
          call MPI_Recv (buff,nptsg,RTYPE,
     .                   nodes(ibloc)-1,TAG_GRID,mycomm,istat,ierr)
c
          call datamb(ibeg(igrp),iend(igrp),iskip(igrp),
     .            jbeg(igrp),jend(igrp),jskip(igrp),
     .            kbeg(igrp),kend(igrp),kskip(igrp),
     .            imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
     .            imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
     .            mbloc,msegt,msegn,nseg,imap,ibloc,
     .            buffw,buffp,buff,
     .            wn,alpha,totptim,ifirst                 )
        else
c
c       data is local
        call datamb(ibeg(igrp),iend(igrp),iskip(igrp),
     .            jbeg(igrp),jend(igrp),jskip(igrp),
     .            kbeg(igrp),kend(igrp),kskip(igrp),
     .            imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
     .            imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
     .            mbloc,msegt,msegn,nseg,imap,ibloc,
     .            w(m5cc(1,ibloc)),p(m1cc(1,ibloc)),x(m3np(1,ibloc)),
     .            wn,alpha,totptim,ifirst                 )
c
        endif
c
c     send the flow data from other nodes to master
      else if (nodes(ibloc)-1.eq.myrank) then
c
        call MPI_Send (w(m5cc(1,ibloc)),5*npts,RTYPE,
     .                 master,TAG_SAVE,mycomm,ierr)
        call MPI_Send (p(m1cc(1,ibloc)),npts,RTYPE,
     .                 master,TAG_SAVE,mycomm,ierr)
        call MPI_Send (x(m3np(1,ibloc)),nptsg,RTYPE,
     .                 master,TAG_GRID,mycomm,ierr)
c
      endif
c
#else   /* IN-CORE version */
c
      call datamb(ibeg(igrp),iend(igrp),iskip(igrp),
     .            jbeg(igrp),jend(igrp),jskip(igrp),
     .            kbeg(igrp),kend(igrp),kskip(igrp),
     .            imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
     .            imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
     .            mbloc,msegt,msegn,nseg,imap,ibloc,
     .            w(m5cc(1,ibloc)),p(m1cc(1,ibloc)),x(m3np(1,ibloc)),
     .            wn,alpha,totptim,ifirst                 )
c
#endif
c
  120 continue
c
      write(idatfil,'(" ")')
c
      return
      end 
