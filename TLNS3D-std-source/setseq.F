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
      subroutine setseq( mstrt,mcyc,mgrlev,ngrid,nlev,mgtype,mgrleg,
     .                   nrleg,isoln,icoars,ifine,ibegr )
c
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      dimension icoars(mgrleg,mgrlev),ifine(mgrleg,mgrlev),
     .          ibegr (mgrleg,mgrlev),nrleg(mgrlev),
     .          mcyc(mgrlev),isoln(mgrlev)
c
c     mgtype is used for identifying the type of multigrid employed
c     currently 3 options are included
c    
c     mgtype = 0  use only single grid, no multigrid used
c              1  multigrid with v-cycle
c              2  multigrid with w-cycle
c
c     full multigrid (fmg) is activated if nlev.gt.1
c
c     making sure that nlev doesnot exceed dimensions
cBWW  this was already checked in main
c
c     if( nlev .gt. mgrlev ) then
c      write (iwrit,*)
c      write (iwrit,'(2x,"The no. of fmg levels input for this run")')
c      write (iwrit,'(2x,"excedes the number of grids the program has",
c    . 1x, "been set up to run")')
c      write (iwrit,'(2x,"the no. of fmg levels in the input file = ",
c    .   i2)')nlev
c      write (iwrit,'(2x,"the max no. of grids the program can run = ",
c    .   i2)')mgrlev
c      stop
c     endif
c
c     check to see that ngrid is within dimensions
c
c     mxgrdav = mgrlev
c
cBWW  also checked in main
c     if( ngrid .gt. mxgrdav ) then
c      write (iwrit,*)
c      write (iwrit,'(2x,"The number of grids input ")')
c      write (iwrit,'(2x,"excedes the number of grids the program has",
c    . 1x, "been set up to run")')
c      write (iwrit,'(2x,"the no. of grids in the input file = ",
c    .   i2)')ngrid
c      write (iwrit,'(2x,"the max no. of grids the program can run = ",
c    .   1x, "at this fmg level is ",i2)') mxgrdav
c      stop
c     endif
c
cBWW  also checked now in main
c     if( mstrt .gt. 1 .and. nlev .ge. 1 ) then
c      write (iwrit,'(2x,"Can Not Restart and use Full Multigrid")')
c      write (iwrit,'(2x,"Program Terminated")')
c      stop
c     endif
c
c
c     identify the solution grids for all the fmg levels
c
      do 20 ilev  = 1,nlev
        isoln(ilev) = nlev - ilev + 1
  20  continue
c
c
c     check to see that soln. grid on all fmg levels is within dimensions
cBWW  this should never happen since nlev was already checked
cBWW  against mgrlev in the main
c
      mxgrdav    = mgrlev
c
      do 30 ilev = 1,nlev
c     write(iwrit,'(2i5)') isoln(ilev),mxgrdav
      if( isoln(ilev) .gt. mxgrdav ) then
       write (iwrit,*)
       write (iwrit,'(2x,"The solution grid indicated for this level")')
       write (iwrit,'(2x,"excedes the number of grids the program has",
     . 1x, "been set up to run")')
       write (iwrit,'(2x,"the solution grid = ",
     .   i2)')isoln(ilev)
       write (iwrit,'(2x,"the max no. of grids the program can run = ",
     .   1x, "at this fmg level is ",i2)') mxgrdav
       stop
      endif
c
      mxgrdav  = mxgrdav - 1
  30  continue
c
c
      if( mgtype .eq. 0 ) then
c  No Multigrid
        do 200 ilev    = 1,nlev
        nrleg (ilev)   = 1
        ibegr (1,ilev) = isoln(ilev)
        icoars(1,ilev) = isoln(ilev)
        ifine (1,ilev) = isoln(ilev)
  200   continue
c
      else if( mgtype .eq. 1 ) then
c  V-cycle
       do 400 ilev     = 1,nlev
        nrleg(ilev)    = 1
        ibegr (1,ilev) = isoln(ilev)
        icoars(1,ilev) = ngrid
        ifine (1,ilev) = isoln(ilev)
  400 continue
c
      else if( mgtype .ge. 2 ) then
c  W-cycle
        do 415 ilev  = 1,nlev
         n           = ngrid - isoln(ilev) + 1
         if(n.eq.1) then
            nrleg (ilev)   = 1
            ibegr (1,ilev) = isoln(ilev)
            icoars(1,ilev) = isoln(ilev)
            ifine (1,ilev) = isoln(ilev)
         else
            nrleg(ilev) = 2**( n-2)
            do 410 n  = isoln(ilev),ngrid -1
            do 410 il = 1,nrleg(ilev),2**( n-isoln(ilev) )
             icoars(il,ilev) = ngrid
             ibegr (il,ilev) = ngrid - n + ( isoln(ilev) - 1 )
             ifine(nrleg(ilev)-il+1,ilev) = ibegr(il,ilev)
  410       continue
         end if
  415   continue
c
      endif
c
c
cBWW  this output is now done in main at the start of
cBWW  the iterations on each grid level
c
c     write (iwrit,*)
c     if (nlev.gt.1) go to 2000
c      if( mgtype .eq. 1 ) then
c        write (iwrit,'(2x,"V-Cycles")')
c        write (iwrit,'(2x,"No. of Cycles =",i4)')mcyc(nlev)
c        write (iwrit,'(2x,"Solution grid =",i4)')isoln(nlev)
c      else if( mgtype .eq. 2 ) then
c        write (iwrit,'(2x,"W-Cycles")')
c        write (iwrit,'(2x,"No. of Cycles = ",i2)')mcyc(nlev)
c        write (iwrit,'(2x,"Solution grid = ",i2)')isoln(nlev)
c      endif
c      goto 2100
c2000  continue
c      if( mgtype .eq. 1 ) then
c        write (iwrit,'(2x,"Full Multigrid with V-Cycles")')
c      else
c        write (iwrit,'(2x,"Full Multigrid with W-Cycles")')
c      endif
c      write (iwrit,'("Number of Full Multigrid Levels =",i2)')nlev
c      write (iwrit,'("level #   Cycles    Soln-Grid   ")')
c       do 420 ilev=1,nlev
c       write (iwrit,'(6x,i4,2(i4,5x))')ilev,mcyc(ilev),isoln(ilev)
c       write (iwrit,'(i5,2(5x,i5))')ilev,mcyc(ilev),isoln(ilev)
c 420   continue
c2100 continue
c     write (iwrit,*)
c       do 435 ilev=1,nlev
c     write (iwrit,'(2x,"ibegr",2x,"icoars",2x,"ifine")')
c       do 430 irleg=1,nrleg(ilev)
c     write (iwrit,'(3x,i2,6x,i2,5x,i2)')
c    . ibegr(irleg,ilev),icoars(irleg,ilev),ifine(irleg,ilev)
c 430 continue
c     write (iwrit,*)
c 435 continue
c
c
      return
      end
