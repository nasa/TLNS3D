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
      subroutine savsol2(imn,jmn,kmn,imp2,jmp2,kmp2,m1cc,m5cc,mxsizc,
     .                   mxsiz5c,mgrlev,mbloc,nbloc,mres,nres,
     .                   w,p,eomu,turv1,turv2,smin,wcount,rtrmsg,
     .                   hrmsg,clg,cmzg,cdtg,cdvg,nsupg,
     .                   mxszfg,wold,tv1old,tv2old,m5fgt,m1fgt,mtorder,
     .                   ntorder,iturv1,iturv2,iunsteady,totptim,
# if defined BUILD_MPI
     .                   scal,xref,yref,zref,sref,
     .                   buffw,buffp,buffe,buff1,buff2,buffs,
     .                   mx3dwk,buffwo,bufftv1,bufftv2,nodes )
# else
     .                   scal,xref,yref,zref,sref,nodes )
# endif
c=======================================================================
c 
c            saves solution for restart
c
c     modified for distributed computing : April, 1995
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
#endif
      character*128  errmsg
      common/dstcom/ errmsg,master,myrank,mycomm,ndlist,nnodes
      dimension      nodes(1)
c 
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension imp2(mgrlev,mbloc), jmp2(mgrlev,mbloc),
     .          kmp2(mgrlev,mbloc)
c
      dimension m1cc(mgrlev,mbloc), m5cc(mgrlev,mbloc)
c
      dimension w(mxsiz5c), p(mxsizc), eomu(mxsizc ), turv1(mxsizc),
     .                                 smin(mxsizc),  turv2(mxsizc)
c
      dimension wcount(mres),rtrmsg(mres),hrmsg(mres),nsupg(mres),
     .             clg(mres),  cmzg(mres), cdtg(mres), cdvg(mres)

      dimension     m1fgt(mtorder,mbloc), m5fgt(mtorder,mbloc)
      dimension     wold (mxszfg*5*ntorder+1),
     .              tv1old(mxszfg*iturv1*ntorder+1),
     .              tv2old(mxszfg*iturv2*ntorder+1)
c
# if defined BUILD_MPI
       dimension buffw(5*mxszfg),buffp(mxszfg),buffe(mxszfg),
     . buff1(mxszfg),buff2(mxszfg),buffs(mxszfg)

       dimension buffwo(mx3dwk*5,ntorder),
     .           bufftv1(mx3dwk*iturv1,ntorder),
     .           bufftv2(mx3dwk*iturv2,ntorder)
# endif
c
      if (myrank.eq.master) then
c
        rewind isavfil
        write(isavfil) nres,iturb,i2dfl,ntorder,totptim
        write(isavfil) gamma,rm,acos(ca),pr,prt,
     .                 rey/(1.e+06*scal),suthc,tref
        write(isavfil) scal,xref,yref,zref,sref
        write(isavfil)
     .    (wcount(ires),rtrmsg(ires),hrmsg(ires),  clg(ires),
     .     cmzg(ires),  cdtg(ires), cdvg(ires),nsupg(ires),ires=1,nres)
        write(isavfil) nbloc
        write(isavfil)
     .    (imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),ibloc=1,nbloc)
c
      endif
c
      do 100 ibloc=1,nbloc
      npts  = imp2(1,ibloc)*jmp2(1,ibloc)*kmp2(1,ibloc)
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
      if (myrank.eq.master) then
#ifdef BUILD_PVM
        call PVMFrecv (nodes(ibloc),TAG_SAVE,ierr)
c
        call PVMFunpack (RTYPE,w,5*npts,1,ierr)
        call PVMFunpack (RTYPE,p,npts,1,ierr)
        call PVMFunpack (RTYPE,eomu,npts,1,ierr)
        write (isavfil)
     .    (w(n),n=1,5*npts),(p(n),n=1,npts),(eomu(n),n=1,npts)
#else
      if (nodes(ibloc)-1.ne.myrank) then
        call MPI_Recv (buffw,5*npts,RTYPE,nodes(ibloc)-1,
     .                 TAG_SAVE,mycomm,istat,ierr)
        call MPI_Recv (buffp,npts,RTYPE,nodes(ibloc)-1,
     .                 TAG_SAVE,mycomm,istat,ierr)
        call MPI_Recv (buffe,npts,RTYPE,nodes(ibloc)-1,
     .                 TAG_SAVE,mycomm,istat,ierr)
        write (isavfil)
     .    (buffw(n),n=1,5*npts),(buffp(n),n=1,npts),(buffe(n),n=1,npts)
      else
        write (isavfil)
     .    (w(m5cc(1,ibloc)+n-1),n=1,5*npts),(p(m1cc(1,ibloc)+n-1),
     .    n=1,npts),(eomu(m1cc(1,ibloc)+n-1),n=1,npts)
      end if
#endif
c
        if (iturb.eq.2) then
#ifdef BUILD_PVM
          call PVMFunpack (RTYPE,smin,npts,1,ierr)
          call PVMFunpack (RTYPE,turv1,npts,1,ierr)
          write (isavfil)
     .      (smin(n),n=1,npts),(turv1(n),n=1,npts)
#else
      if (nodes(ibloc)-1.ne.myrank) then
          call MPI_Recv (buffs,npts,RTYPE,nodes(ibloc)-1,
     .                   TAG_SAVE,mycomm,istat,ierr)
          call MPI_Recv (buff1,npts,RTYPE,nodes(ibloc)-1,
     .                   TAG_SAVE,mycomm,istat,ierr)
          write (isavfil)
     .      (buffs(n),n=1,npts),(buff1(n),n=1,npts)
      else
          write (isavfil)
     .      (smin(m1cc(1,ibloc)+n-1),n=1,npts),
     .      (turv1(m1cc(1,ibloc)+n-1),n=1,npts)
      end if
#endif
        endif
c
        if (iturb.eq.3) then
#ifdef BUILD_PVM
          call PVMFunpack (RTYPE,smin,npts,1,ierr)
          call PVMFunpack (RTYPE,turv1,npts,1,ierr)
          call PVMFunpack (RTYPE,turv2,npts,1,ierr)
          write (isavfil)
     .      (smin(n),n=1,npts),(turv1(n),n=1,npts),(turv2(n),n=1,npts)
#else
      if(nodes(ibloc)-1.ne.myrank) then
          call MPI_Recv (buffs,npts,RTYPE,nodes(ibloc)-1,
     .                   TAG_SAVE,mycomm,istat,ierr)
          call MPI_Recv (buff1,npts,RTYPE,nodes(ibloc)-1,
     .                   TAG_SAVE,mycomm,istat,ierr)
          call MPI_Recv (buff2,npts,RTYPE,nodes(ibloc)-1,
     .                   TAG_SAVE,mycomm,istat,ierr)
          write (isavfil)
     .      (buffs(n),n=1,npts),(buff1(n),n=1,npts),(buff2(n),n=1,npts)
      else
          write (isavfil)
     .      (smin(m1cc(1,ibloc)+n-1),n=1,npts),
     .      (turv1(m1cc(1,ibloc)+n-1),n=1,npts),
     .      (turv2(m1cc(1,ibloc)+n-1),n=1,npts)
      end if
#endif
        endif
c
c       end of logic for master node, begining of logic
c       for node on which the data resides
c
      else if (nodes(ibloc)-1.eq.myrank) then
c
#ifdef BUILD_PVM
        call PVMFinitsend (PvmDataInplace,ierr)
        call PVMFpack (RTYPE,w(m5cc(1,ibloc)),5*npts,1,ierr)
        call PVMFpack (RTYPE,p(m1cc(1,ibloc)),npts,1,ierr)
        call PVMFpack (RTYPE,eomu(m1cc(1,ibloc)),npts,1,ierr)
c
        if (iturb.eq.2 .or. iturb.eq.3) then
          call PVMFpack (RTYPE,smin(m1cc(1,ibloc)),npts,1,ierr)
          call PVMFpack (RTYPE,turv1(m1cc(1,ibloc)),npts,1,ierr)
          if (iturb.eq.3)
     .      call PVMFpack (RTYPE,turv2(m1cc(1,ibloc)),npts,1,ierr)
        endif
c
        call PVMFsend (master,TAG_SAVE,ierr)
#else
        call MPI_Send (w(m5cc(1,ibloc)),5*npts,RTYPE,
     .                 master,TAG_SAVE,mycomm,ierr)
        call MPI_Send (p(m1cc(1,ibloc)),npts,RTYPE,
     .                 master,TAG_SAVE,mycomm,ierr)
        call MPI_Send (eomu(m1cc(1,ibloc)),npts,RTYPE,
     .                 master,TAG_SAVE,mycomm,ierr)
c
        if (iturb.eq.2 .or. iturb.eq.3) then
          call MPI_Send (smin(m1cc(1,ibloc)),npts,RTYPE,
     .                   master,TAG_SAVE,mycomm,ierr)
          call MPI_Send (turv1(m1cc(1,ibloc)),npts,RTYPE,
     .                   master,TAG_SAVE,mycomm,ierr)
          if (iturb.eq.3)
     .      call MPI_Send (turv2(m1cc(1,ibloc)),npts,
     .                     RTYPE,master,
     .                     TAG_SAVE,mycomm,ierr)
        endif
#endif
c
      endif
c
#else	/* IN-CORE version */
c
      nwbeg = m5cc(1,ibloc)
      nwend = nwbeg + 5*npts - 1
      npbeg = m1cc(1,ibloc)
      npend = npbeg + npts - 1
c
      write (isavfil) (   w(n),n=nwbeg,nwend),
     .                (   p(n),n=npbeg,npend),
     .                (eomu(n),n=npbeg,npend)
c
      if (iturb.eq.2)
     .  write (isavfil) ( smin(n),n=npbeg,npend),
     .                  (turv1(n),n=npbeg,npend)
c
      if (iturb.eq.3)
     .  write (isavfil) ( smin(n),n=npbeg,npend),
     .                  (turv1(n),n=npbeg,npend),
     .                  (turv2(n),n=npbeg,npend)
c
#endif
c
c     logic for time-dependent terms
c
      if(iunsteady.eq.0.or.ntorder.lt.1) go to 100
      do 200 iorder=1,ntorder
c
# if defined BUILD_MPI
c begin coding for distributed (MPI) version
      if (myrank.eq.master) then
c     begin coding for master node
c
      if (nodes(ibloc)-1.ne.myrank) then
        call MPI_Recv (buffwo(1,iorder),5*npts,RTYPE,
     .                 nodes(ibloc)-1,TAG_SAVE,mycomm,istat,ierr)
        write (isavfil)
     .  (buffwo(n,iorder),n=1,5*npts)
        
        if(iturb.ge.2) then
          call MPI_Recv (bufftv1(1,iorder),npts,RTYPE,
     .                   nodes(ibloc)-1,TAG_SAVE,mycomm,istat,ierr)
          write (isavfil)
     .    (bufftv1(n,iorder),n=1,npts)
        endif
        if(iturb.eq.3) then 
          call MPI_Recv (bufftv2(1,iorder),npts,RTYPE,
     .                 nodes(ibloc)-1,TAG_SAVE,mycomm,istat,ierr)
          write (isavfil)
     .    (bufftv2(n,iorder),n=1,npts)
        endif
c
      else
        write (isavfil)
     .  (wold(m5fgt(iorder,ibloc)+n-1),n=1,5*npts)

        if(iturb.ge.2)
     .    write (isavfil)
     .    (tv1old(m1fgt(iorder,ibloc)+n-1),n=1,npts)

        if(iturb.eq.3)
     .    write (isavfil)
     .    (tv2old(m1fgt(iorder,ibloc)+n-1),n=1,npts)

      end if
c
c       end of logic for master node, begining of logic
c       for node on which the data resides
c
      else if (nodes(ibloc)-1.eq.myrank) then

        call MPI_Send (wold(m5fgt(iorder,ibloc)),5*npts,
     .                   RTYPE,
     .                   master,TAG_SAVE,mycomm,ierr)
        if(iturb.ge.2)
     .    call MPI_Send (tv1old(m1fgt(iorder,ibloc)),npts,
     .                   RTYPE,
     .                   master,TAG_SAVE,mycomm,ierr)
      
        if(iturb.eq.3)
     .    call MPI_Send (tv2old(m1fgt(iorder,ibloc)),npts,
     .                   RTYPE,
     .                   master,TAG_SAVE,mycomm,ierr)
      
      endif
c
c
#else	/* IN-CORE version */

        write (isavfil)
     .  (wold(m5fgt(iorder,ibloc)+n-1),n=1,5*npts)

        if(iturb.ge.2)
     .    write (isavfil)
     .    (tv1old(m1fgt(iorder,ibloc)+n-1),n=1,5*npts)

        if(iturb.eq.3)
     .    write (isavfil)
     .    (tv2old(m1fgt(iorder,ibloc)+n-1),n=1,5*npts)
      
#endif
  200 continue
  100 continue
c
c
      return
      end 
