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
      subroutine initpvm (mbloc,nnodes,nodeid,iwrit)
c====================================================================
c     start up PVM node processes
c
#     include "fpvm3.h"
c
      dimension nodeid(mbloc)
      character*80 fname,host,proc
c
      input = 13
      iread = 14
c
c---- get command line argument
c     format: rundist -pg[...] procfile < input
c         or: rundist -pr[...] procfile < input
c
      narg = iargc ()
      if (narg.gt.1) then
        do n=1,narg-1
          call getarg (n,fname)
          if (fname(1:3).eq.'-pg' .or.
     .        fname(1:3).eq.'-p4') then
            iarg = n+1
            go to 10
          endif
        enddo
      endif
c
      write (iwrit,'("need to specify [-pg procfile] argument")')
      go to 300
c
c---- open procfile and preprocess
c
   10 call getarg (iarg,fname)
      len = istlen (fname,80)
      open (input,file=fname(1:len),form='formatted')
      call pinput (input,iread,'#',.true.,.true.,.true.)
      close (input)
c
c---- read nodes from procfile
c
      nnodes = 0
   20 continue
      read (iread,'(a)',end=100,err=200) host
      read (iread,'(i10)',end=200,err=200) nproc
      read (iread,'(a)',end=200,err=200) proc
      if (nproc.le.0) go to 20
c
      if (nnodes+nproc.gt.mbloc)
     .  nproc = mbloc - nnodes
c
c     add host if neccessary
c
      call PVMFmstat (host,mstat)
      if (mstat.ne.PvmOk) then
        call PVMFaddhost (host,ierr)
        if (ierr.lt.1) then
          len = istlen (host,80)
          write (iwrit,'("PVM error: failed to add host ",a)')
     .      host(1:len)
          go to 210
        endif
      endif
c
c     spawn the task(s)
c
      call PVMFspawn (proc,PvmHost,host,nproc,nodeid(nnodes+1),numt)
      if (numt.lt.nproc) then
        lh = istlen (host,80)
        lp = istlen (proc,80)
        write (iwrit,'("PVM error: failed to start ",i3,
     .    " processes of ",a," on ",a)') nproc,proc(1:lp),host(1:lh)
        if (numt.gt.0) nnodes = nnodes + numt
        go to 210
      endif
c
      nnodes = nnodes + nproc
      if (nnodes.lt.mbloc) go to 20
c
  100 close (iread)
      if (nnodes.eq.0) then
        write (iwrit,'("no node processes specified")')
        go to 300
      endif
      return
c
  200 close (iread)
      write (iwrit,'("invalid procfile or missing data")')
c
  210 if (nnodes.gt.0) then
        do inode=1,nnodes
          call PVMFkill (nodeid(inode),ierr)
        enddo
      endif
c
  300 call PVMFexit (ierr)
      stop
      end
