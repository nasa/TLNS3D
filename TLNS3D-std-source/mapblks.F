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
c======================================================================
      subroutine mapblks (mgrlev,mbloc,nbloc,im,jm,km,
     .                    nnodes,nodes,iwrk,input,iwrit,ierr)
c======================================================================
c     map grid blocks to node processes
c

      dimension im(mgrlev,mbloc),jm(mgrlev,mbloc),km(mgrlev,mbloc)
      dimension nodes(mbloc),iwrk(mbloc,3)
c
      character*132 line
      character*32 buf
c
      ierr = 0
c
c---- simple case of 1 node
c
      if (nnodes.eq.1) then
        do ibloc=1,nbloc
          nodes(ibloc) = 1
        enddo
        return
      endif
c
c---- zero out block assignments
c
      do ibloc=1,nbloc
        nodes(ibloc)  = 0
        iwrk(ibloc,3) = 0
      enddo
c
c---- check for -b2n or -btn option
c
      if (input.gt.0) then
c
      narg = iargc ()
      if (narg.gt.1) then
        do iarg=1,narg-1
          call getarg (iarg,line)
          if (line(1:4).eq.'-b2n' .or.
     .        line(1:4).eq.'-btn') go to 10
        enddo
      endif
      go to 100
c
c     read block-to-node mapfile
c
   10 call getarg (iarg+1,line)
      open (input+1,file=line,form='formatted')
      rewind (input+1)
      call pinput (input+1,input,'#',.true.,.true.,.false.)
      close (input+1)
c
   20 continue
      read (input,'(a132)',end=90) line
      i1 = 1
      i2 = nextwrd (line,132,i1,ilen)
      if (i2.eq.0) go to 20
      if (ilen.gt.32) go to 80
      do i=1,32
        buf(i:i) = ' '
      enddo
      do i=1,ilen
        buf(i:i) = line(i1+i-1:i1+i-1)
      enddo
      read (buf,'(i32)',err=80) inode
      if (inode.le.0 .or. inode.gt.nnodes) then
        if (iwrit.gt.0)
     .    write (iwrit,'(" node ",i4,
     .      " from block-to-node file out of range")') inode
        ierr = ierr + 1
        go to 20
      endif
c
   30 continue
      i1 = i2
      i2 = nextwrd (line,132,i1,ilen)
      if (ilen.gt.32) go to 80
      do i=1,32
        buf(i:i) = ' '
      enddo
      do i=1,ilen
        buf(i:i) = line(i1+i-1:i1+i-1)
      enddo
      if (i2.eq.0) then
        read (buf,'(i32)',err=80) nb1
        nb2 = nb1
      else
        read (buf,'(i32)',err=80) nb1
        if (line(i2:i2).eq.'-') then
          i1 = nextwrd (line,132,i2,ilen)
          if (i1.eq.0) then
            if (iwrit.gt.0)
     .        write (iwrit,'("missing block range value",
     .          " from block-to-node file")')
            ierr = ierr + 1
            go to 20
          endif
          i2 = nextwrd (line,132,i1,ilen)
          if (ilen.gt.32) go to 80
          do i=1,32
            buf(i:i) = ' '
          enddo
          do i=1,ilen
            buf(i:i) = line(i1+i-1:i1+i-1)
          enddo
          read (buf,'(i32)',err=80) nb2
        else
          nb2 = nb1
        endif
      endif
c
      if (nb1.lt.1   .or. nb1.gt.nbloc .or.
     .    nb2.lt.nb1 .or. nb2.gt.nbloc) then
        if (iwrit.gt.0)
     .    write (iwrit,'("block range ",i4," to ",i4,
     .      " from block-to-node file out of range")') nb1,nb2
        ierr = ierr + 1
        go to 20
      endif
c
c     assign blocks
c
      do nb=nb1,nb2
        nodes(nb) = inode
      enddo
      if (i2.eq.0) go to 20
      go to 30
c
   80 if (iwrit.gt.0)
     .  write (iwrit,'("invalid block-to-node map file")')
      ierr = ierr + 1
c
   90 close (input)
c
      endif
c
c---- find block sizes and number of un-mapped blocks
c
  100 icnt = 0
      do ibloc=1,nbloc
        iwrk(ibloc,1) = ibloc
        iwrk(ibloc,2) = im(1,ibloc)*jm(1,ibloc)*km(1,ibloc)
        if (nodes(ibloc).eq.0) then
          icnt = icnt + 1
        else
          inode         = nodes(ibloc)
          iwrk(inode,3) = iwrk(inode,3) + iwrk(ibloc,2)
        endif
      enddo
c
c     if all blocks are assigned, return
c
      if (icnt.eq.0) return
c
c---- sort blocks largest to smallest
c
      do ibloc=1,nbloc-1
        do jbloc=ibloc+1,nbloc
          if (iwrk(ibloc,2).lt.iwrk(jbloc,2)) then
            itmp          = iwrk(ibloc,1)
            iwrk(ibloc,1) = iwrk(jbloc,1)
            iwrk(jbloc,1) = itmp
            itmp          = iwrk(ibloc,2)
            iwrk(ibloc,2) = iwrk(jbloc,2)
            iwrk(jbloc,2) = itmp
          endif
        enddo
      enddo
c
c---- assign unassigned blocks in order (largest to smallest)
c     to the node with the fewest points
c
      do jbloc=1,nbloc
        ibloc = iwrk(jbloc,1)
        if (nodes(ibloc).eq.0) then
          imin  = 1
          do inode=2,nnodes
            if (iwrk(inode,3).lt.iwrk(imin,3)) imin = inode
          enddo
          nodes(ibloc) = imin
          iwrk(imin,3) = iwrk(imin,3) + iwrk(jbloc,2)
        endif
      enddo
c
      return
      end
c======================================================================
      integer function nextwrd (str,len,is,ilen)
c======================================================================
c     get next word from string and length of current word
c
      character*(*) str
      character*11  digits
c
      data digits /'0123456789-'/
c
      nextwrd = 0
      ilen    = 1
      if (is.lt.1 .or. is.ge.len) return
c
c     skip current word
c
      if (str(is:is).eq.'-') then
        i = is
        go to 20
      endif
c
      do 10 i=is+1,len
      if (str(i:i).eq.'-') then
        nextwrd = i
        return
      endif
      do j=1,10
        if (str(i:i).eq.digits(j:j)) go to 10
      enddo
      go to 20
   10 ilen = ilen + 1
      i  = len
c
   20 ii = i + 1
      if (ii.ge.len) return
c
c     find next word
c
      do i=ii,len
        do j=1,11
          if (str(i:i).eq.digits(j:j)) then
            nextwrd = i
            return
          endif
        enddo
      enddo
c
      return
      end
