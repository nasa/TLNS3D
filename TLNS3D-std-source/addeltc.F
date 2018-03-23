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
      subroutine addeltc (imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,ntorder,delw,vol,w,dtphy,
     .  rkphysa,mrkstg,nrkstg,irkstg )
c#ifdef  RKPHY
c# else
c    .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,ntorder,delw,vol,w,dtphy)
c# endif
c
c     code for accelerated time-accurate calculations
c
c
c     Time derivatives are added to the convective fluxes (delw)
c     and viscous fluxes for all 5 equations
c     in this subroutine
c
c     time difference is made between the old value from the previous 
c     time step and the value at the zeroth stage of the RK (previous
c     iteration (or cycle))
c
c     Note: on the coarse grid only the current stage of the RK is accounted for
c
c     initial coding :  Aug. 1992
c
c     latest mods    : by v.n.vatsa (May, 2000) for distributed version
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/timek/ ct(9)
c
      dimension delw(imp2,jmp2,kmp2,5), vol(imp2,jmp2,kmp2),
     .             w(imp2,jmp2,kmp2,5)
c
      dimension rkphysa(mrkstg,mrkstg)
c
#ifdef  RKPHY
c
c     contribution of fluxes due to physical R-K
c
      if (nrkstg .gt. 2) then

      do 70 nv=1,5
      do 70 k = kmn+1,kmx
      do 70 j = jmn+1,jmx
      do 70 i = imn+1,imx
         delw(i,j,k,nv) = rkphysa(irkstg,irkstg) * delw(i,j,k,nv)
   70 continue

      endif
#endif
c
c
c     backward difference scheme of order "ntorder"
c
      do k = kmn+1,kmx
      do j = jmn+1,jmx
      do i = imn+1,imx
        vdt           = ct(ntorder+1)*vol(i,j,k)/dtphy
        delw(i,j,k,1) = delw(i,j,k,1) + w(i,j,k,1)*vdt
        delw(i,j,k,2) = delw(i,j,k,2) + w(i,j,k,2)*vdt
        delw(i,j,k,3) = delw(i,j,k,3) + w(i,j,k,3)*vdt
        delw(i,j,k,4) = delw(i,j,k,4) + w(i,j,k,4)*vdt
        delw(i,j,k,5) = delw(i,j,k,5) + w(i,j,k,5)*vdt
      enddo
      enddo
      enddo
        if (ntorder .gt. 6) then
c
        write(iwrit,'(
     .         " this order of time differencing is not supported")')
      stop
      endif
c
      return
      end
