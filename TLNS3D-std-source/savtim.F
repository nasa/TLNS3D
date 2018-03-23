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
      subroutine  savtim(imn, jmn, kmn, im,  jm,  km,  
     .                   imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .                   mxsizc,mxsiz5c,mxsiz3n,
     .                   mgrlev,mbloc,msegt,msegn,
     .                   m5cc,m1cc,m3np,
     .                   imap,nbloc,nseg,
# if defined BUILD_MPI
     .                   mxszfg,buffw,buffp,buff,
# endif
     .                   x,w,p,wn,alpha,totptim,ifirst,nodes    )
c=======================================================================
c 
c     object: To save the solution at prespecified time intervals
c     for time-dependent solutions
c     grid is saved only at the initial time step (ifirst=1)
c     original coding : M.D. Sanetrik
c     MPI (distributed) version : By V.N.Vatsa (May, 2000)
c
c
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
      dimension  im  (mgrlev,mbloc), jm  (mgrlev,mbloc),
     .           km  (mgrlev,mbloc), imp1(mgrlev,mbloc),
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
      if (myrank.eq.master) then
        if(ifirst.eq.1) rewind igpfil
        if(ifirst.eq.1) rewind iqpfil
c
        write (iqpfil) nbloc
        write (iqpfil) (im(1,ibloc)-1,jm(1,ibloc)-1,km(1,ibloc)-1,
     .                      ibloc=1,nbloc)
c
        if(ifirst .eq. 1) then
          rewind igpfil
          write (igpfil) nbloc
          write (igpfil) (im(1,ibloc)-1,jm(1,ibloc)-1,km(1,ibloc)-1,
     .                        ibloc=1,nbloc)
        end if
c
      endif
c
      do 120 ibloc = 1,nbloc
c
c     write out flow solution at a given time level
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
          call flowmb(imn,jmn,kmn,
     .            im  (1,ibloc),jm  (1,ibloc),km  (1,ibloc),
     .            imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
     .            imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
     .            mbloc,msegt,msegn,nseg,imap,ibloc,
     .            buffw,buffp,buff,
     .            wn,alpha,totptim,ifirst                 )
        else
c
c       data is local
        call flowmb(imn,jmn,kmn,
     .            im  (1,ibloc),jm  (1,ibloc),km  (1,ibloc),
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
      call flowmb(imn,jmn,kmn,
     .            im  (1,ibloc),jm  (1,ibloc),km  (1,ibloc),
     .            imp1(1,ibloc),jmp1(1,ibloc),kmp1(1,ibloc),
     .            imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),
     .            mbloc,msegt,msegn,nseg,imap,ibloc,
     .            w(m5cc(1,ibloc)),p(m1cc(1,ibloc)),x(m3np(1,ibloc)),
     .            wn,alpha,totptim,ifirst                 )
c
#endif
  120 continue
c
c
      return
      end 
