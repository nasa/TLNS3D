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
      program mapper
c
c     prints block to process mapping for different number
c     of processes
c
      parameter (mbloc  = 100)
      parameter (mnode  = 32)
      parameter (nwords = 56)
c

      dimension ic(mbloc),jc(mbloc),kc(mbloc),nc(mbloc)
      dimension nodes(mbloc),iwrk(mbloc,3)
      dimension nlist(3,mnode,mnode)
c
      character*80 title,gridin
c
c---- check for output range command line arguments
c
      ipbeg = 1
      ipend = 0
      narg  = iargc ()
      if (narg.gt.0) then
        call getarg (1,title)
        read (title,'(i10)') ipend
        if (ipend.gt.0 .and. narg.gt.1) then
          call getarg (2,title)
          read (title,'(i10)') ipbeg
          if (ipbeg.le.0) ipbeg = 1
          if (ipbeg.gt.ipend) then
            narg  = ipbeg
            ipbeg = ipend
            ipend = narg
          endif
        endif
      endif
c
c---- read TLNS3D-MB input file
c
      igrid = 2
      iread = 3
      input = 5
      iwrit = 6
c
      call pinput (input,iread,'#',.true.,.false.,.false.)
c
      read (iread,'(a80)') title
      len = istlen (title,80)
      write (iwrit,'(a)') title(1:len)
      read (iread,'(1x)')
      read (iread,'(8f10.6)') grdfmt,fstrt,fsave,fplot,fcplpr
      read (iread,'(1x)')
      read (iread,'(8f10.6)') flev,fgrids,fmgtyp
      igrdfmt = grdfmt
      nlev    = flev
      if (nlev.le.0) nlev = 1
      ngrlev  = fgrids
      read (iread,'(1x)')
      do n=1,nlev
        read (iread,'(1x)')
      enddo
      do n=1,24
        read (iread,'(1x)')
      enddo
c
c     read grid file name
c
      read (iread,'(a80)') gridin
      call fixstr (gridin,80)
      len = istlen (gridin,80)
      if (igrdfmt.eq.1) then
        write (iwrit,'("formatted grid from ",a)') gridin(1:len)
        open (igrid,file=gridin(1:len),form='formatted')
        rewind (igrid)
        read (igrid,*,end=2) nbloc
        if (nbloc.gt.mbloc) then
          write (iwrit,'("block dimension exceeded")')
          stop
        endif
        read (igrid,*,end=2) (ic(n),jc(n),kc(n),n=1,nbloc)
      else
        write (iwrit,'("unformatted grid from ",a)') gridin(1:len)
        open (igrid,file=gridin(1:len),form='unformatted')
        rewind (igrid)
        read (igrid,end=2) nbloc
        if (nbloc.gt.mbloc) then
          write (iwrit,'("block dimension exceeded")')
          stop
        endif
        read (igrid,end=2) (ic(n),jc(n),kc(n),n=1,nbloc)
      endif
      close (igrid)
      goto 3
c
    2 write (iwrit,'("grid file ",a," not found or EOF")')
     .  gridin(1:len)
      stop
c
c---- print grid block sizes
c
    3 write (iwrit,'(/,"Grid Block Sizes")')
      write (iwrit,'("block   imax   jmax   kmax     total")')
      write (iwrit,'("-----   ----   ----   ----     -----")')
c
      ncmax = 0
      do n=1,nbloc
        nc(n) = ic(n)*jc(n)*kc(n)
        ncmax = ncmax + nc(n)
        write (iwrit,'(i5,3i7,i10)') n,ic(n),jc(n),kc(n),nc(n)
      enddo
c
c---- do block to node mapping for various number of nodes
c
      write (iwrit,'(/,
     .  "nodes   maxpts   minpts   avgpts   %avgdev   ",
     .  "megawords  exetime",/,
     .  "-----   ------   ------   ------   -------   ",
     .  "---------  -------")')

      do nnodes=1,nbloc
        call mapblks (1,mbloc,nbloc,ic,jc,kc,
     .                nnodes,nodes,iwrk,0,iwrit,ierr)
        if (ierr.ne.0) stop
c
        do n=1,nnodes
          iwrk(n,1) = 0
          iwrk(n,2) = 0
        enddo
c
        do n=1,nbloc
          m = nodes(n)
          iwrk(m,1) = iwrk(m,1) + 1
          iwrk(m,2) = iwrk(m,2) + nc(n)
        enddo
c
        navg = ncmax / nnodes
        ndev = 0
        nmin = iwrk(1,2)
        nmax = iwrk(1,2)
        do n=1,nnodes
          iwrk(n,3) = iwrk(n,2) - navg
          ndev      = ndev + iabs(iwrk(n,3))
          if (nmin.gt.iwrk(n,2)) nmin = iwrk(n,2)
          if (nmax.lt.iwrk(n,2)) nmax = iwrk(n,2)
        enddo
        dev = 100.0 * real(ndev) / real(ncmax)
c
        write (iwrit,'(i5,3i9,3f10.3)') nnodes,nmax,nmin,navg,dev,
     .    real(nmax*nwords)/1.e+06,real(nmax)/real(ncmax)
c
        if (nnodes.le.mnode) then
          do n=1,nnodes
            do m=1,3
              nlist(m,n,nnodes) = iwrk(n,m)
            enddo
          enddo
        endif
c
      enddo
c
      if (ipend.gt.0) then
        if (ipend.gt.nbloc) ipend = nbloc
        if (ipbeg.gt.ipend) ipbeg = ipend
        do nnodes=ipbeg,ipend
          write (iwrit,'(/,"number nodes = ",i2)') nnodes
          write (iwrit,'("node blocks  points  deviation",/,
     .               "---- ------  ------  ---------")')
          do n=1,nnodes
            write (iwrit,'(i3,i6,2i10)')
     .        n,(nlist(m,n,nnodes),m=1,3)
          enddo
        enddo
      endif
c
      stop
      end
