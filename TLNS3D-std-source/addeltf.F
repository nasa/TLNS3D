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
      subroutine addeltf (imn,jmn,kmn,imx,jmx,kmx,
     .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,ntorder,delw,vol,w,wold,dtphy,
     .   mebdf,ibdf2opt,icyc,rkphysa,mrkstg,nrkstg,irkstg,delwork,wbar,
     .   itime )
c#ifdef  RKPHY
c# else
c    .  imp1,jmp1,kmp1,imp2,jmp2,kmp2,ntorder,delw,vol,w,wold,dtphy)
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
c     initial coding :  Aug. 1992
c     latest mods    : by v.n.vatsa (May, 2000) for distributed version
c
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      dimension delw(imp2,jmp2,kmp2,5), vol(imp2,jmp2,kmp2),
     .             w(imp2,jmp2,kmp2,5),
     .          wold(imp2,jmp2,kmp2,5,ntorder)
c
      dimension delwork(imp2,jmp2,kmp2,5,nrkstg),
     .          rkphysa(mrkstg,mrkstg)
c
      dimension wbar (imp2,jmp2,kmp2,5,ntorder-1)

#ifdef  RKPHY
c
ca    if a physical Runge Kutta time integration is used, add fluxes at
c     previous stages and multiply present flux by diagonal coefficient
c
c     fill in the rhs residual for initial solution in delwork array
      if(irkstg.eq.2.and.icyc.eq.1) then
        do nv=1,5
          do k = kmn+1,kmx
            do j = jmn+1,jmx
              do i = imn+1,imx
              delwork(i,j,k,nv,1) = delw(i,j,k,nv)
              enddo
            enddo
          enddo
        enddo
      endif
c
c     save the rhs residual for current stage for later use
      if(irkstg.lt.nrkstg) then
        do nv=1,5
          do k = kmn+1,kmx
            do j = jmn+1,jmx
              do i = imn+1,imx
              delwork(i,j,k,nv,irkstg) = delw(i,j,k,nv)
              enddo
            enddo
          enddo
        enddo
      endif
c
      if (nrkstg .gt. 2) then

      do 70 nv=1,5
      do 70 k = kmn+1,kmx
      do 70 j = jmn+1,jmx
      do 70 i = imn+1,imx
         tdir = 0.
         do 71 l = 1, irkstg-1
             tdir = tdir + rkphysa(irkstg,l) * delwork(i,j,k,nv,l)
   71    continue
         delw(i,j,k,nv) = rkphysa(irkstg,irkstg) * delw(i,j,k,nv) +
     .                    tdir
   70 continue

      endif
#endif
c
      if (mebdf.gt.0) go to 1001
c
      if (ntorder .eq. 1) then
c
c     first order backward difference
c
      do 10 nv=1,5
      do 10 k = kmn+1,kmx
      do 10 j = jmn+1,jmx
      do 10 i = imn+1,imx
      tdir    =  (w(i,j,k,nv) - wold(i,j,k,nv,1)) * 
     .          vol(i,j,k) / dtphy
      delw(i,j,k,nv) = delw(i,j,k,nv) + tdir
  10  continue
c
      elseif (ntorder .eq. 2) then
c
c     coefficients for constant dt
c     a4 = 1.5/dt, a3 = -2./dt, a2=.5/dt
c
c     second order backward difference : BDF2
c     variable dt
c
c     H0=dtv(1)
c     H1=dtv(2)
      H0  = dtphy
      H1  = dtphy
c 
      A4 =  (H1+2*H0)/(H0*(H1+H0))
      A3 = -(H1+H0)/(H0*H1)
      A2 =   H0/(H1*(H1+H0))
c
      do 20 nv=1,5
      do 20 k = kmn+1,kmx
      do 20 j = jmn+1,jmx
      do 20 i = imn+1,imx
      tdir    =  vol(i,j,k) *
     .      (A4*w(i,j,k,nv) + A3*wold(i,j,k,nv,1) + A2*wold(i,j,k,nv,2))
      delw(i,j,k,nv) = delw(i,j,k,nv) + tdir
  20  continue
c
      elseif (ntorder .eq. 3) then
c
c     logic for ibdf2opt
c
      if(ibdf2opt.eq.1.and.itime.ge.3) then
c       a4  =  1.69/dtphy
c       a3  = -2.58/dtphy
c       a2  =  1.08/dtphy
c       a1  = -0.19/dtphy

        bdfoptc = -0.58d0/3.0d0
        a4  = ( 1.5d0 - bdfoptc)/dtphy
        a3  = (-2.0d0 + bdfoptc*3.0d0)/dtphy
        a2  = ( 0.5d0 - bdfoptc*3.0d0)/dtphy
        a1  = bdfoptc/dtphy
      else
c
c     coefficients for constant dt
c     a4 = 11./(6.*dt), a3 = -3./dt, a2=1.5/dt, a1=-1./(3.*dt)
c
c     third order backward difference : BDF3
c     variable dt
c
c     H0=dtv(1)
c     H1=dtv(2)
c     H2=dtv(3)
c
      H0 = dtphy
      H1 = dtphy
      H2 = dtphy
      if(itime.le.2) then
         dt0 = dtphy
         dt1 = dtphy
         dt2 = dtphy
c        dt1 = 0.015625d0
c        dt2 = 0.015625d0
         if(itime.eq.2) dt1 = dtphy
         h0  = dt0
         h1  = dt1
         h2  = dt2
       endif
c
      A4 = (H1*H2+2*H0*H2+H1**2+4*H0*H1+3*H0**2)/(H0*(H1+H0)*(H2+H1+H0))
      A3 = -(H1*H2+H0*H2+H1**2+2*H0*H1+H0**2)/(H0*H1*(H2+H1))
      A2 = (H0*H2+H0*H1+H0**2)/(H1*(H1+H0)*H2)
      A1 = -(H0*H1+H0**2)/(H2*(H2+H1)*(H2+H1+H0))
c
      endif

      do 30 nv=1,5
      do 30 k = kmn+1,kmx
      do 30 j = jmn+1,jmx
      do 30 i = imn+1,imx
          tdir = vol(i,j,k) * 
     .         (  A4    *    w(i,j,k,nv) 
     .          + A3    * wold(i,j,k,nv,1) 
     .          + A2    * wold(i,j,k,nv,2)
     .          + A1    * wold(i,j,k,nv,3) )
          delw(i,j,k,nv) = delw(i,j,k,nv) + tdir
   30 continue
c
      elseif (ntorder .eq. 4) then
c
c         fourth order backward difference : BDF4
c         constant dt
c
      do 40 nv=1,5
      do 40 k = kmn+1,kmx
      do 40 j = jmn+1,jmx
      do 40 i = imn+1,imx
          tdir = (
     .            25./12.  *    w(i,j,k,nv) 
     .          - 4.       * wold(i,j,k,nv,1) 
     .          + 3.       * wold(i,j,k,nv,2)
     .          - 4./3.    * wold(i,j,k,nv,3)
     .          + 1./4.    * wold(i,j,k,nv,4)
     .                   ) *  vol(i,j,k) / dtphy
          delw(i,j,k,nv) = delw(i,j,k,nv) + tdir
   40 continue
c
      elseif (ntorder .eq. 5) then
c
c         fifth order backward difference
c
      do 50 nv=1,5
      do 50 k = kmn+1,kmx
      do 50 j = jmn+1,jmx
      do 50 i = imn+1,imx
          tdir = (
     .            137./60. *    w(i,j,k,nv) 
     .          - 5.       * wold(i,j,k,nv,1) 
     .          + 5.       * wold(i,j,k,nv,2)
     .          - 10./3.   * wold(i,j,k,nv,3)
     .          + 5./4.    * wold(i,j,k,nv,4)
     .          - 1./5.    * wold(i,j,k,nv,5)
     .                   ) *  vol(i,j,k) / dtphy
          delw(i,j,k,nv) = delw(i,j,k,nv) + tdir
   50 continue
c
      elseif (ntorder .eq. 6) then
c
c         sixth order backward difference
c
      do 60 nv=1,5
      do 60 k = kmn+1,kmx
      do 60 j = jmn+1,jmx
      do 60 i = imn+1,imx
          tdir = (
     .            49./20.  *    w(i,j,k,nv) 
     .          - 6.       * wold(i,j,k,nv,1) 
     .          + 15./2.   * wold(i,j,k,nv,2)
     .          - 20./3.   * wold(i,j,k,nv,3)
     .          + 15./4.   * wold(i,j,k,nv,4)
     .          - 6./5.    * wold(i,j,k,nv,5)
     .          + 1./6.    * wold(i,j,k,nv,6)
     .                   ) *  vol(i,j,k) / dtphy
          delw(i,j,k,nv) = delw(i,j,k,nv) + tdir
   60     continue
c
        elseif (ntorder .gt. 6) then
c
        write(iwrit,'(
     .         " this order of time differencing is not supported")')
      stop
      endif
      go to 1002
c
 1001 continue
c
c     MEBDF4 scheme
c
c     H0=dtv(1)
c     H1=dtv(2)
c     H2=dtv(3)
c
      H0 = dtphy
      H1 = dtphy
      H2 = dtphy
      if(itime.le.2) then
         H0 = dtphy
         H1 = dtphy
         H2 = dtphy
c        H1 = 0.015625d0
c        H2 = 0.015625d0
         if(itime.eq.2) H1 = dtphy
       endif
c
      if(irkstg.eq.1) then
c
c     MEBDF4 : stage 1
c
c     Save the fluxes for later stages before adding time terms
c
      do nv=1,5
        do k = kmn+1,kmx
          do j = jmn+1,jmx
            do i = imn+1,imx
c
            delwork(i,j,k,nv,irkstg) = delw(i,j,k,nv)
c
            enddo
          enddo
        enddo
      enddo
c
      A4 = (H1*H2+2*H0*H2+H1**2+4*H0*H1+3*H0**2)/(H0*(H1+H0)*(H2+H1+H0))
      A3 = -(H1*H2+H0*H2+H1**2+2*H0*H1+H0**2)/(H0*H1*(H2+H1))
      A2 = (H0*H2+H0*H1+H0**2)/(H1*(H1+H0)*H2)
      A1 = -(H0*H1+H0**2)/(H2*(H2+H1)*(H2+H1+H0))

      do nv=1,5
        do k = kmn+1,kmx
          do j = jmn+1,jmx
            do i = imn+1,imx
            tdir = vol(i,j,k) *
     .           ( A4    *    w(i,j,k,nv) 
     .           + A3    * wold(i,j,k,nv,1) 
     .           + A2    * wold(i,j,k,nv,2)
     .           + A1    * wold(i,j,k,nv,3) )
            delw(i,j,k,nv) = delw(i,j,k,nv) + tdir
c
            enddo
          enddo
        enddo
      enddo
c
      else if(irkstg.eq.2) then
c     MEBDF4 : stage 2
c
c     Save the fluxes for later stages before adding time terms
c
      do nv=1,5
        do k = kmn+1,kmx
          do j = jmn+1,jmx
            do i = imn+1,imx
c
            delwork(i,j,k,nv,irkstg) = delw(i,j,k,nv)
c
            enddo
          enddo
        enddo
      enddo
c
      HF = dtphy
      A4 = (H0*H1+2*HF*H1+H0**2+4*HF*H0+3*HF**2)/(HF*(H0+HF)*(H1+H0+HF))
      A3 = -(H0*H1+HF*H1+H0**2+2*HF*H0+HF**2)/(HF*H0*(H1+H0))
      A2 = (HF*H1+HF*H0+HF**2)/(H0*(H0+HF)*H1)
      A1 = -(HF*H0+HF**2)/(H1*(H1+H0)*(H1+H0+HF))

      do nv=1,5
        do k = kmn+1,kmx
          do j = jmn+1,jmx
            do i = imn+1,imx
            tdir = vol(i,j,k) *
     .           ( A4    *    w(i,j,k,nv) 
     .           + A3    * wbar(i,j,k,nv,1) 
     .           + A2    * wold(i,j,k,nv,1)
     .           + A1    * wold(i,j,k,nv,2) )
            delw(i,j,k,nv) = delw(i,j,k,nv) + tdir
c
            enddo
          enddo
        enddo
      enddo
c
      else if(irkstg.eq.3) then
c     MEBDF4 : stage 3
c
        RbetaK = (H1*H2+2*H0*H2+H1**2+4*H0*H1+3*H0**2)/(H0*(H1+H0)*
     &           (H2+H1+H0))
      HF  = dtphy
      t1  = 0.
      t2  = t1 + h2
      t3  = t2 + h1
      t4  = t3 + h0
      t5  = t4 + hf
      DENOMIN = (t5-t4)*(12*t4**2*t5**2-8*t3*t4*t5**2-8*t2*t4*t5**2+4*t2
     1   *t3*t5**2-8*t3*t4**2*t5-8*t2*t4**2*t5+6*t3**2*t4*t5+10*t2*t3*t4
     2   *t5+6*t2**2*t4*t5-3*t2*t3**2*t5-3*t2**2*t3*t5+4*t2*t3*t4**2-3*t
     3   2*t3**2*t4-3*t2**2*t3*t4+2*t2**2*t3**2)
      betaKH=(t4*(t4-t2)*(t4-t3)*(4*t5**3-3*t4*t5**2-3*t3*t5**2-3*t2*t5*
     1   *2+2*t3*t4*t5+2*t2*t4*t5+2*t2*t3*t5-t2*t3*t4) )/ denomin
      betaKP1H=(-t4**2*(t4-t2)**2*(t4-t3)**2 )/ denomin
      AH4 = 1
      AH3 =(-t4**2*(t4-t2)**2*(t5-t4)*(4*t5**2-2*t4*t5-3*t2*t5+t2*t4)/(t
     1   3*(t3-t2)) )/ denomin
      AH2 =(t4**2*(t4-t3)**2*(t5-t4)*(4*t5**2-2*t4*t5-3*t3*t5+t3*t4)/(t2
     1   *(t3-t2)) )/ denomin
      AH1 =(-(t4-t2)**2*(t4-t3)**2*(t5-t4)*(4*t5**2-2*t4*t5-3*t3*t5-3*t2
     1   *t5+t3*t4+t2*t4+2*t2*t3)/(t2*t3) )/ denomin

      A4 = aH4*Rbetak
      A3 = aH3*Rbetak
      A2 = aH2*Rbetak
      A1 = aH1*Rbetak
      B1  = betaKH*RbetaK - 1
      B2  = betaKp1H*RbetaK
c
      do nv=1,5
        do k = kmn+1,kmx
          do j = jmn+1,jmx
            do i = imn+1,imx
            tdir = vol(i,j,k) *
     .           ( A4 *    w(i,j,k,nv) 
     .           + A3 * wold(i,j,k,nv,1) 
     .           + A2 * wold(i,j,k,nv,2)
     .           + A1 * wold(i,j,k,nv,3) )
     .           + B1 * delwork(i,j,k,nv,1) + B2 * delwork(i,j,k,nv,2)
            delw(i,j,k,nv) = delw(i,j,k,nv) + tdir
c
            enddo
          enddo
        enddo
      enddo
c
      endif
c     end of logic for MEBDF4
c
 1002 continue
c
      return
      end
