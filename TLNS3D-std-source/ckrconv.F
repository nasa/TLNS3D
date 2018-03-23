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
      subroutine  ckrconv(imn,jmn,kmn,imx,jmx,kmx,imp2,jmp2,kmp2,
     .                    delw,dtl,dtmin,rtrmsg,mx1dwk,rt,ibloc,nodes)
c=======================================================================
c 
c     objective : checks convergence of residuals
c                 quantities with "g" at the end refer to the
c                 global values over all the blocks
c     iniitial coding : v.n.vatsa (March 1990)
c 
c 

      dimension delw(imp2,jmp2,kmp2,5),dtl(imp2,jmp2,kmp2)
c
      dimension rt(mx1dwk)
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
      common/rk/    cfl,c(6),qfil(6),beta(6),vt,hm,mstage
c
c     initialize global values of residuals for ibloc = 1
c
      if (ibloc .eq. 1) then
         rtrmsg     = 0.
         rtmaxg     = 0.
         irtmxg     = 0
         jrtmxg     = 0
         krtmxg     = 0
         mrtbloc    = 1
         totpts     = 0
         rtrms2     = 0.
      endif
c
      rtmax     = 0.
      do 90 k=kmn+1,kmx
      do 90 j=jmn+1,jmx
      iloc   = 0
      do 50 i=imn+1,imx
      iloc   = iloc + 1
      dt        = cfl*((1.-vt)*dtl(i,j,k) + vt*dtmin)
      rt(iloc)  = delw(i,j,k,1)/dt
      rtrms2    = rtrms2 + rt(iloc)**2
   50 continue
      i         = isamax(iloc,rt,1)
      if (abs(rt(i)).gt.abs(rtmax)) then 
        rtmax     = abs(rt(i))
        irt       = i + imn
        jrt       = j
        krt       = k
      endif 
   90 continue
c
c     keep track of total no. of interior points on all blocks and
c     evaluate global rms values including points at current block
c
      totpts    = totpts + (imx - imn) * (jmx - jmn) * (kmx - kmn)
      rtrmsg    = sqrt(rtrms2/totpts)
c
      if (ibloc .eq. 1) then
         irtmxg   = irt
         jrtmxg   = jrt
         krtmxg   = krt
         rtmaxg   = rtmax
         return
      endif
c
c     check if rtmax on the current bloc is greater than rtmaxg
c
      if(abs(rtmax) .gt. abs(rtmaxg)) then
         irtmxg   = irt
         jrtmxg   = jrt
         krtmxg   = krt
         rtmaxg   = rtmax
         mrtbloc  = ibloc
      endif
c
c     check for Nan's
c
      if ( (rtmax+rtmax.eq.rtmax.and.abs(rtmax).gt.1.d-100).or.
     .      .not.(rtmax.lt.abs(rtmax).or.rtmax.ge.0.d0)) then
       write (iwrit,901) ibloc
 901   format(2x,'Nan detected in rtmax at ibloc = ',i5)
c      stop 1001
       call ERREXIT(nodes) 
      endif
c
      if ( (rtrmsg+rtrmsg.eq.rtrmsg.and.abs(rtrmsg).gt.1.d-100).or.
     .      .not.(rtrmsg.lt.abs(rtrmsg).or.rtrmsg.ge.0.d0)) then
       write (iwrit,902) ibloc
 902   format(2x,'Nan detected in rtrmsg at ibloc = ',i5)
c      stop 1001
       call ERREXIT(nodes) 
      endif
c
      return
      end
