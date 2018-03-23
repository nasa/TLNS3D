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
      subroutine gridin( igrdfmt,imn,jmn,kmn,imx,jmx,kmx,
     .                   imp1,jmp1,kmp1,x)
c
c   reading from unit "igrdfil" the x,y and z coordinates
c   in the array x(i,j,k,3)
c
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension x(imp1,jmp1,kmp1,3)
c
c
      if(i2dfl.eq.1) go to 100
c
      if (igrdfmt .eq. 1) then
          read (igrdfil,*)
     .    (((x(i,j,k,1),i=imn,imx),j=jmn,jmx),k=kmn,kmx),
     .    (((x(i,j,k,2),i=imn,imx),j=jmn,jmx),k=kmn,kmx),
     .    (((x(i,j,k,3),i=imn,imx),j=jmn,jmx),k=kmn,kmx)
      else
          read (igrdfil)
     .    (((x(i,j,k,1),i=imn,imx),j=jmn,jmx),k=kmn,kmx),
     .    (((x(i,j,k,2),i=imn,imx),j=jmn,jmx),k=kmn,kmx),
     .    (((x(i,j,k,3),i=imn,imx),j=jmn,jmx),k=kmn,kmx)
      endif
      return
c
  100 continue
c     read 2-d grid
c
      if (igrdfmt .eq. 1) then
          read (igrdfil,*)
     .    ((x(i,j,kmn,1),i=imn,imx),j=jmn,jmx),
     .    ((x(i,j,kmn,2),i=imn,imx),j=jmn,jmx)
      else
          read (igrdfil)
     .    ((x(i,j,kmn,1),i=imn,imx),j=jmn,jmx),
     .    ((x(i,j,kmn,2),i=imn,imx),j=jmn,jmx)
      endif
c
      do 120 k=kmn,kmx
      zvar   = k-1.
      do 120 j=jmn,jmx
      do 120 i=imn,imx
      x(i,j,k,1) = x(i,j,kmn,1)
      x(i,j,k,2) = x(i,j,kmn,2)
      x(i,j,k,3) = zvar
  120 continue
c
      return
      end
