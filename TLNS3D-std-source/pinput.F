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
      subroutine pinput (iread,iwrit,cmt,bline,lspace,brksp)
c=======================================================================
c     strips comments and trailing space and optionally
c     blank lines and leading space from input
c     arguments:
c         iread  - unit to read input from
c         iwrit  - unit to write output to. Opens a scratch
c                  file to write to if not unit 6
c         cmt    - comment character
c         bline  - boolean to remove blank lines (.true. removes)
c         lspace - boolean to remove leading space (.true. removes)
c         brksp  - boolean to break lines at space (.true. breaks)
c
c
      character*1 sp,ht,nl,cmt
      character*1 buf(132)
      character*7 format
      logical bline,lspace,brksp
c
      data sp /' '/,ht /'\t'/, nl /'\n'/
c
      if (iwrit.ne.6)
     .    open (unit=iwrit,status='scratch')
c
   10 continue
      read (iread,'(132a1)',end=110) buf
      ii  = 1
      len = 132
c
c     strip anything following comment character
c
      do i=1,len
        if (buf(i).eq.cmt) then
          if (i.eq.1) go to 10
          len = i - 1
          go to 30
        endif
      enddo
c
c     strip all trailing white space
c
   30 do ie=len,1,-1
        if (buf(ie).ne.sp .and.
     .      buf(ie).ne.ht .and.
     .      buf(ie).ne.nl) go to 50
      enddo
      if (.not.bline)
     .    write (iwrit,'(/)')
      go to 10
c
c     skip leading white space
c
   50 do is=ii,ie
        if (buf(is).ne.sp .and.
     .      buf(is).ne.ht) go to 70
      enddo
      if (.not.bline)
     .   write (iwrit,'(/)')
      go to 10
c
c     break lines at interior space
c
   70 if (brksp) then
        do ii=is,ie
          if (buf(ii).eq.sp .or.
     .        buf(ii).eq.ht) go to 80
        enddo
        go to 90
   80   len = ii - is
        if (len.gt.99) then
          write (format,'("(",i3,"a1)")') len
        else if (len.gt.9) then
          write (format,'("(",i2,"a1)")') len
        else
          write (format,'("(",i1,"a1)")') len
        endif
        write (iwrit,format) (buf(i),i=is,ii-1)
        go to 50
      endif
c
c     ouput remaining string
c
   90 if (.not.lspace) is = 1
      len = ie - is + 1
      if (len.gt.99) then
        write (format,'("(",i3,"a1)")') len
      else if (len.gt.9) then
        write (format,'("(",i2,"a1)")') len
      else
        write (format,'("(",i1,"a1)")') len
      endif
      write (iwrit,format) (buf(i),i=is,ie)
      go to 10
c
  110 continue
      if (iwrit.ne.6) rewind (iwrit)
      return
      end
