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
      subroutine  ckhconv(imn,jmn,kmn,imx,jmx,kmx,imp2,jmp2,kmp2,
     .                    w,p,nsupg,hrmsg,mx1dwk,h,ibloc,nodes )
c=======================================================================
c 
c     objective : checks convergence of enthalpy
c                 quantities with "g" at the end refer to the
c                 global values over all the blocks
c     iniitial coding : v.n.vatsa (March 1990)
c 
c 

      dimension   w(imp2,jmp2,kmp2,5), p(imp2,jmp2,kmp2)
c
      dimension h(mx1dwk)
c
      common/conv/  rtmaxg,hmaxg,rtrms2, hrms2, totpts,
     .              irtmxg,jrtmxg,krtmxg,mrtbloc,
     .               ihmxg, jhmxg, khmxg, mhbloc
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      common/ma/    amachg
      common/rk/    cfl,c(6),qfil(6),beta(6),vt,hm,mstage
c
      sgrmi     = 1./(rm*sqrt(gamma))
c
c     initialize global values of residuals for ibloc = 1
c
      if (ibloc .eq. 1) then
         hrmsg      = 0.
         hmaxg      = 0.
         amachg     = 0.
         nsupg      = 0
         ihmxg      = 0
         jhmxg      = 0
         khmxg      = 0
         mhbloc     = 1
         totpts     = 0
         hrms2      = 0.
      endif
c
      hmax       = 0.
      amach      = 0.
c            loop over all blocks
      do 90 k=kmn+1,kmx
      do 90 j=jmn+1,jmx
      iloc   = 0
      do i=imn+1,imx
        iloc     = iloc + 1
        qq       = (w(i,j,k,2)**2 + w(i,j,k,3)**2 + w(i,j,k,4)**2)
     .           /  w(i,j,k,1) 
        am       = qq/(gamma*p(i,j,k))
        amach    = max(am,amach)
c       nsupg    = nsupg + min(1,int(qq/(gamma*p(i,j,k))))
        nsupg    = nsupg + min(1,int(am))
        h(iloc)  = (w(i,j,k,5) + p(i,j,k))/w(i,j,k,1) - h0
        hrms2    = hrms2 + h(iloc)**2
      enddo
      i          = isamax(iloc,h,1)
      if (abs(h(i)).gt.abs(hmax)) then 
        hmax     = abs(h(i))
        ih       = i + imn
        jh       = j
        kh       = k
      endif 
   90 continue
c
c     keep track of total no. of interior points on all blocks and
c     evaluate global rms values including points at current block
c
      totpts    = totpts + (imx - imn) * (jmx - jmn) * (kmx - kmn)
      hrmsg     = sqrt(hrms2 /totpts)
c
      if (ibloc .eq. 1) then
         ihmxg    = ih
         jhmxg    = jh
         khmxg    = kh
         hmaxg    = hmax
         amachg   = amach
         return
      endif
c
c     check if hmax on the current bloc is greater than hmaxg
c
      if(abs(hmax) .gt. abs(hmaxg)) then
         ihmxg    = ih
         jhmxg    = jh
         khmxg    = kh
         hmaxg    = hmax
         mhbloc   = ibloc
      endif

      if(abs(amach) .gt. abs(amachg)) then
         amachg   = amach
      endif
c
c     check for Nan's
c
      if ( (hmax+hmax.eq.hmax.and.abs(hmax).gt.1.d-100).or.
     .      .not.(hmax.lt.abs(hmax).or.hmax.ge.0.d0)) then
       write (iwrit,901) ibloc
 901   format(2x,'Nan detected in hmax at ibloc = ',i5)
c      stop 1001
       call ERREXIT(nodes)
      endif
c
      if ( (hrmsg+hrmsg.eq.hrmsg.and.abs(hrmsg).gt.1.d-100).or.
     .      .not.(hrmsg.lt.abs(hrmsg).or.hrmsg.ge.0.d0)) then
       write (iwrit,902) ibloc
 902   format(2x,'Nan detected in hrmsg at ibloc = ',i5)
c      stop 1001
       call ERREXIT(nodes)
      endif
c
      return
      end
