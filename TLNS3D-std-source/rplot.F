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
      subroutine  rplot(iwrit,nres,nsup,res,count,resout,supout)
c=======================================================================
c 
c     plots convergence rate and number supersonic points 
c
      parameter (nlines=12)
c

      dimension   res(1),nsup(1),count(1),resout(1),supout(1) 
      dimension   label(6) 
      character*1 line(101),kode(5)
      data        kode/' ','-','|','*','+'/
c 
      rmin      = 0.
      rmax      = 0.
      res1      = res(1)
      fsupf     = real(nsup(nres))
      if (fsupf.gt.1.e-06) fsupf = 1./fsupf 
      do 10 n=1,nres
      resout(n) = log10(res(n)/res1) 
      rmin      = min(rmin,resout(n)) 
      rmax      = max(rmax,resout(n)) 
   10 supout(n) = fsupf*real(nsup(n))
c
      imin      = rmin
      imax      = rmax
      if ((rmin-real(imin)).lt.0.) imin = imin - 1
      if ((rmax-real(imax)).gt.0.) imax = imax + 1
      if (mod(imax - imin,2).eq.0)  imax = imax + 1
      rmin      = imin
      rmax      = imax
c
      mlines    = (imax - imin)*nlines
      mpos      = mlines/2
      dr        = (rmax - rmin)/real(mlines)
      ru        = rmax + .5*dr
      lr        = imax
      fr        = 1.2
      if (imax-imin.gt.1)
     .  fr      = 1. + 1./(rmax - rmin - 1.)
      df        = fr/real(mlines)
      fu        = fr + .5*df 
      dfr       = fr/(rmax -rmin) 
c
      nc        = count(nres)
      lc        = nc/50 + 1
      if (mod(nc,50).eq.0) lc = lc -1
      cmax      = lc*50
      lc        = lc*10
      do k=1,6
        label(k) = lc*(k -1)
      enddo
c
      write (iwrit,'("1",//," convergence plot",/)')
      do k=1,101 
        line(k) = kode(2)
      enddo
      do k=1,101,20
        line(k) = kode(3)
      enddo
      write (iwrit,'(1x,i9,2x,101a1,f6.3)') lr,line,fr
c
      kk        = 1
      do 100 m=1,mlines
      if (mod(m,nlines).eq.0) kk = 2
      do k=1,101
        line(k) = kode(kk)
      enddo
      do k=1,101,20
        line(k) = kode(3)
      enddo
c
      ru        = ru  -dr 
      rl        = ru  -dr 
      fu        = fu  -df 
      fl        = fu  -df 
      do 70 n=1,nres
      if (supout(n).ge.fl .and. supout(n).le.fu) then
        k       = int(1.  +100.*count(n)/cmax)
        k       = min(k,101) 
        line(k) = kode(5)
      endif
c
      if (resout(n).ge.rl .and. resout(n).le.ru) then
        k       = int(1.  +100.*count(n)/cmax)
        k       = min(k,101)
        line(k) = kode(4)
      endif
   70 continue
c
      if (m.eq.mpos) then
        write (iwrit,'(5x,a1,6x,101a1,5x,a1)') kode(4),line,kode(5)
      else if (m.eq.mpos+1) then
        write (iwrit,'(" log(rn/r1) ",101a1,"  sup/supf")') line
      else if (kk.eq.1) then
        write (iwrit,'(12x,101a1)') line
      else
        kk      = 1
        lr      = lr  -1
        fr      = fr  -dfr
        write (iwrit,'(1x,i9,2x,101a1,f6.3)') lr,line,fr
      endif
  100 continue
c
      write (iwrit,'(3x,i10,5i20,/,60x,"work")') label 
      return
      end 
