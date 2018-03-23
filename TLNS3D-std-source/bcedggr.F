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
      subroutine  bcedggr(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                    imp2,jmp2,kmp2,vol,igrid,isoln   )
c=======================================================================
c 
c
c     object : to fill-in the flow variables at block corners
c              using extrapolation
c
c     initial coding : by v.n.vatsa, nasa langley ( Oct., 1994)
c
c

      dimension      vol(imp2,jmp2,kmp2)
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      imnp1   = imn+1
      jmnp1   = jmn+1
      kmnp1   = kmn+1
c
      do 100 k=kmn,kmp1,kmp1-kmn
c
c$dir novector
      do 110 j=jmn,jmp1,jmp1-jmn
c
      vol(imn  ,j,k) = vol(imnp1,j,k)
c
      vol(imp1,j,k)  = vol(imx  ,j,k)
 110  continue
c
 100  continue
c
      if(igrid.ne.isoln) return
c
c
      imnm1   = imn-1
      jmnm1   = jmn-1
      kmnm1   = kmn-1
c
      do 200 k=kmnm1,kmp2,kmp2-kmnm1
c
c$dir novector
      do 210 j=jmnm1,jmp2,jmp2-jmnm1
c
      vol(imnm1 ,j,k) = vol(imnp1,j,k)
c
      vol(imp2,j,k)   = vol(imp1,j,k)
 210  continue
c
 200  continue
c
      return
      end
