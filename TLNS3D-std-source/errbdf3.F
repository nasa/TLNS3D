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
      subroutine  errbdf3 (imn,jmn,kmn,imx,jmx,kmx,imp2,jmp2,kmp2,
     .                    delw,dtl,vol,dtphy,mx1dwk,errw,ibloc,nodes,
     .                    w,wold,ntorder)
c=======================================================================
c 
c     objective : computes temporal error in BDF3
c                 time-stepping using solutions 
c                 with different order accuracies
c
c     iniitial coding : v.n.vatsa (Nov. 2003)
c 
c 

      dimension delw(imp2,jmp2,kmp2,5),dtl(imp2,jmp2,kmp2)
c
      dimension errw(mx1dwk),vol(imp2,jmp2,kmp2),
     .          w(imp2,jmp2,kmp2,5),
     .          wold (imp2,jmp2,kmp2,5,ntorder)
c
      common/rkerr/ erkmaxl(5),erkmaxg(5),erkrms2(5),erkrmsg(5),
     .              totptse
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      common/rk/    cfl,c(6),qfil(6),beta(6),vt,hm,mstage
c
c     initialize global values of errors for ibloc = 1
c
      if (ibloc .eq. 1) then
         do nv = 1,5
           erkrms2(nv)      = 0.
           erkrmsg(nv)      = 0.
           erkmaxg(nv)      = 0.
           totptse          = 0
         enddo
      endif
c
      do 100 nv=1,5
c
      erkmaxl(nv)      = 0.
      c1 = 4./3. - 18./11.
      c2 = 9./11. - 1./3.
      c3 =-2./11.
      c4 = 2./3. - 6./11.
      do 90 k=kmn+1,kmx
      do 90 j=jmn+1,jmx
      iloc   = 0
      do 50 i=imn+1,imx
      iloc   = iloc + 1
      dt        = cfl*((1.-vt)*dtl(i,j,k) + vt*dtmin)
         errw(iloc) =  c1*wold(i,j,k,nv,1) + c2*wold(i,j,k,nv,2)
     .              +  c3*wold(i,j,k,nv,3) + c4*delw(i,j,k,nv)*dtphy
         erkrms2(nv)     = erkrms2(nv) + errw(iloc)**2
   50 continue
      i         = isamax(iloc,errw,1)
      if (abs(errw(i)).gt.abs(erkmaxl(nv))) then 
        erkmaxl(nv)     = abs(errw(i))
      endif 
   90 continue
c
c     keep track of total no. of interior points on all blocks and
c     evaluate global rms values including points at current block
c
      if(nv.eq.1)
     .totptse   = totptse+ (imx - imn) * (jmx - jmn) * (kmx - kmn)
      erkrmsg(nv)    = sqrt(erkrms2(nv)/totptse)
c
      if (ibloc .eq. 1) then
         erkmaxg(nv)   = erkmaxl(nv)
cdbg     return
      endif
c
c     check if erkmaxl on the current bloc is greater than erkmaxg
c
      if(abs(erkmaxl(nv)) .gt. abs(erkmaxg(nv))) then
         erkmaxg(nv)   = erkmaxl(nv)
      endif
c
c     check for Nan's
c
      erkmx = erkmaxl(nv)
      if ( (erkmx+erkmx.eq.erkmx.and.abs(erkmx).gt.1.d-100).or.
     .      .not.(erkmx.lt.abs(erkmx).or.erkmx.ge.0.d0)) then
       write (iwrit,901) ibloc
 901   format(2x,'Nan detected in erkmax at ibloc = ',i5)
c      stop 1001
       call ERREXIT(nodes) 
      endif
c
      erkrms = erkrmsg(nv)
      if ( (erkrms+erkrms.eq.erkrms.and.abs(erkrms).gt.1.d-100).or.
     .      .not.(erkrms.lt.abs(erkrms).or.erkrms.ge.0.d0)) then
       write (iwrit,902) ibloc
 902   format(2x,'Nan detected in erkrms at ibloc = ',i5)
c      stop 1001
       call ERREXIT(nodes) 
      endif
 100  continue
c
      return
      end
