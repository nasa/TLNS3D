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
      subroutine savsolavg(imn,jmn,kmn,imp2,jmp2,kmp2,m1fg,m5fg,
     .                   mxszfg,mxsz5fg,mgrlev,mbloc,nbloc,
     .                   w,p,eomu,wcount,rtrmsg,
     .                   hrmsg,clg,cmzg,cdtg,cdvg,nsupg,
     .                   iunsteady,totptim,
# if defined BUILD_MPI
     .                   scal,xref,yref,zref,sref,
     .                   buffw,buffp,buffe,nodes )
# else
     .                   scal,xref,yref,zref,sref,nodes )
# endif
c=======================================================================
c 
c            saves average flow solution for unsteady flows
c
c     Initial coding: v.n.vatsa (June, 2001)
c

#if defined(BUILD_MPI)
      include "mpif.h"
#  ifdef TLN_DBL
#    define RTYPE MPI_DOUBLE_PRECISION
#  else
#    define RTYPE MPI_REAL
#  endif
      dimension istat(MPI_STATUS_SIZE)
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
      dimension m1fg(mbloc), m5fg(mbloc)
c
      dimension w(mxsz5fg), p(mxszfg), eomu(mxszfg )
c
c
# if defined BUILD_MPI
       dimension buffw(5*mxszfg),buffp(mxszfg),buffe(mxszfg)

# endif
c
      if (myrank.eq.master) then
c
        nres    = 1
        ntorder = 0
        iturbl  = 1
        rewind iavgfil
        write(iavgfil) nres,iturbl,i2dfl,ntorder,totptim
        write(iavgfil) gamma,rm,acos(ca),pr,prt,
     .                 rey/(1.e+06*scal),suthc,tref
        write(iavgfil) scal,xref,yref,zref,sref
        write(iavgfil)
     .     wcount,rtrmsg,hrmsg,clg,cmzg,cdtg,cdvg,nsupg
        write(iavgfil) nbloc
        write(iavgfil)
     .    (imp2(1,ibloc),jmp2(1,ibloc),kmp2(1,ibloc),ibloc=1,nbloc)
c
      endif
c
      do 100 ibloc=1,nbloc
      npts  = imp2(1,ibloc)*jmp2(1,ibloc)*kmp2(1,ibloc)
c
#if defined(BUILD_MPI)
c
      if (myrank.eq.master) then
      if (nodes(ibloc)-1.ne.myrank) then
        call MPI_Recv (buffw,5*npts,RTYPE,nodes(ibloc)-1,
     .                 TAG_SAVE,mycomm,istat,ierr)
        call MPI_Recv (buffp,npts,RTYPE,nodes(ibloc)-1,
     .                 TAG_SAVE,mycomm,istat,ierr)
        call MPI_Recv (buffe,npts,RTYPE,nodes(ibloc)-1,
     .                 TAG_SAVE,mycomm,istat,ierr)
        write (iavgfil)
     .    (buffw(n),n=1,5*npts),(buffp(n),n=1,npts),(buffe(n),n=1,npts)
      else
        write (iavgfil)
     .    (w(m5fg(ibloc)+n-1),n=1,5*npts),(p(m1fg(ibloc)+n-1),
     .    n=1,npts),(eomu(m1fg(ibloc)+n-1),n=1,npts)
      end if
c
c
c       end of logic for master node, begining of logic
c       for node on which the data resides
c
      else if (nodes(ibloc)-1.eq.myrank) then
c
        call MPI_Send (w(m5fg(ibloc)),5*npts,RTYPE,
     .                 master,TAG_SAVE,mycomm,ierr)
        call MPI_Send (p(m1fg(ibloc)),npts,RTYPE,
     .                 master,TAG_SAVE,mycomm,ierr)
        call MPI_Send (eomu(m1fg(ibloc)),npts,RTYPE,
     .                 master,TAG_SAVE,mycomm,ierr)
c
      endif
c
#else	/* IN-CORE version */
c
      nwbeg = m5fg(ibloc)
      nwend = nwbeg + 5*npts - 1
      npbeg = m1fg(ibloc)
      npend = npbeg + npts - 1
c
      write (iavgfil) (   w(n),n=nwbeg,nwend),
     .                (   p(n),n=npbeg,npend),
     .                (eomu(n),n=npbeg,npend)
c
#endif
c
  100 continue
c
c
      return
      end 
