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
      subroutine  bcedgev(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                    imp2,jmp2,kmp2,eomu,turv1,turv2,igrid,isoln)
c=======================================================================
c 
c
c     object : to fill-in the flow variables at block corners
c              using extrapolation
c
c     initial coding : by v.n.vatsa, nasa langley ( Nov., 1994)
c
c

      dimension      turv1(imp2,jmp2,kmp2), turv2(imp2,jmp2,kmp2),
     .               eomu (imp2,jmp2,kmp2)
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
c*******************************************************************
c     mods for filling corner points at i-faces
c     By : vatsa (Jan. 1997)
c
c*******************************************************************
c     fill in jmn,kmn corner

      k  = kmn
      j  = jmn

      do i=imnp1,imx
        eomu(i,j,k)  = .5*(eomu(i,j+1,k)+eomu(i,j,k+1))
      enddo
c     
      if(iturb.ge.2) then
        do i=imnp1,imx
        turv1(i,j,k)  = .5*(turv1(i,j+1,k)+turv1(i,j,k+1))
        enddo
      endif
c     
      if(iturb.eq.3) then
        do i=imnp1,imx
        turv2(i,j,k)  = .5*(turv2(i,j+1,k)+turv2(i,j,k+1))
        enddo
      endif
c
c
c     fill in jmp1,kmn corner

      k  = kmn
      j  = jmp1

      do i=imnp1,imx
        eomu(i,j,k)  = .5*(eomu(i,j-1,k)+eomu(i,j,k+1))
      enddo
c
      if(iturb.ge.2) then
        do i=imnp1,imx
        turv1(i,j,k)  = .5*(turv1(i,j-1,k)+turv1(i,j,k+1))
        enddo
      endif
c
      if(iturb.eq.3) then
        do i=imnp1,imx
        turv2(i,j,k)  = .5*(turv2(i,j-1,k)+turv2(i,j,k+1))
        enddo
      endif
c
c
c     fill in jmn,kmp1 corner

      k  = kmp1
      j  = jmn

      do i=imnp1,imx
        eomu(i,j,k)  = .5*(eomu(i,j+1,k)+eomu(i,j,k-1))
      enddo
c
      if(iturb.ge.2) then
c
        do i=imnp1,imx
        turv1(i,j,k)  = .5*(turv1(i,j+1,k)+turv1(i,j,k-1))
        enddo
c
      endif
c
      if(iturb.eq.3) then
c
        do i=imnp1,imx
        turv2(i,j,k)  = .5*(turv2(i,j+1,k)+turv2(i,j,k-1))
        enddo
c
      endif
c
c     fill in jmp1,kmp1 corner

      k  = kmp1
      j  = jmp1

      do i=imnp1,imx
        eomu(i,j,k)  = .5*(eomu(i,j-1,k)+eomu(i,j,k-1))
      enddo
c
      if(iturb.ge.2) then
c
        do i=imnp1,imx
        turv1(i,j,k)  = .5*(turv1(i,j-1,k)+turv1(i,j,k-1))
        enddo
      endif
c
      if(iturb.eq.3) then
c
        do i=imnp1,imx
        turv2(i,j,k)  = .5*(turv2(i,j-1,k)+turv2(i,j,k-1))
        enddo
      endif
c
c
c             end of logic for corner points
c*******************************************************************
c
      do 100 k=kmn,kmp1,kmp1-kmn
c
cdir$ novector
      do 110 j=jmn,jmp1,jmp1-jmn
c
      eomu(imn ,j,k)  = eomu(imnp1,j,k)
      eomu(imp1,j,k)  = eomu(imx  ,j,k)
 110  continue
c
      if(iturb.eq.2) then
cdir$ novector
        do 111 j=jmn,jmp1,jmp1-jmn
c
        turv1(imn ,j,k)  = turv1(imnp1,j,k)
        turv1(imp1,j,k)  = turv1(imx  ,j,k)
 111    continue
      endif
c
      if(iturb.eq.3) then
cdir$ novector
        do 112 j=jmn,jmp1,jmp1-jmn
c
        turv1(imn ,j,k)  = turv1(imnp1,j,k)
        turv1(imp1,j,k)  = turv1(imx  ,j,k)
        turv2(imn ,j,k)  = turv2(imnp1,j,k)
        turv2(imp1,j,k)  = turv2(imx  ,j,k)
 112    continue
      endif
c
 100  continue
c
      return
      end
