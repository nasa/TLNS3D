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
      subroutine  navier(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,imp2,
     .                   jmp2,kmp2,imap,msegt,msegn,mbloc,nseg,w,p,vol,
     .                   x,six,siy,siz,sjx,sjy,sjz,skx,sky,skz,eomu,
     .                  tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .                   gs1,gs2,gs3,gs4,gs5,hs1,hs2,hs3,hs4,hs5,
     .                   flux2_mn,flux3_mn,flux4_mn,flux5_mn,
     .                   flux2_mx,flux3_mx,flux4_mx,flux5_mx,mx1dwk,
     .                   up,vp,wp,ts,amu,fv,ivisb                    )
c=======================================================================
c 
c 
c       object : compute viscous stress terms based on thin layer 
c                navier-stokes equations
c       viscous terms associated with both eta(normal to wing) and
c       zeta(spanwise) directions are included
c       initial coding : by v.n.vatsa of nasa langley (may-1984)
c       multiblock version : Dec. 1989, by v.n.vatsa, nasa langley
c 
c     w(i,j,k,1)  = density 
c     w(i,j,k,2)  = momentum in x direction 
c     w(i,j,k,3)  = momentum in y direction 
c     w(i,j,k,4)  = momentum in z direction 
c     w(i,j,k,5)  = total energy
c     p(i,j,k)    = pressure
c 
c 

      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c 
c
      dimension  ivisb(3)
c 
      dimension 
     .    w(imp2,jmp2,kmp2,5), p(imp2,jmp2,kmp2),vol(imp2,jmp2,kmp2),
     .   fv(imp2,jmp2,kmp2,4),x(imp1,jmp1,kmp1,3),
     .   up(imp1,jmp1,kmp1),vp(imp1,jmp1,kmp1),wp(imp1,jmp1,kmp1),
     .   ts(imp1,jmp1,kmp1),amu(imp1,jmp1,kmp1),eomu(imp2,jmp2,kmp2),
     .  six(imp1,jmp1,kmp1),siy(imp1,jmp1,kmp1),siz(imp1,jmp1,kmp1),
     .  sjx(imp1,jmp1,kmp1),sjy(imp1,jmp1,kmp1),sjz(imp1,jmp1,kmp1),
     .  skx(imp1,jmp1,kmp1),sky(imp1,jmp1,kmp1),skz(imp1,jmp1,kmp1),
     .  gs1(imp1),gs2(imp1),gs3(imp1),gs4(imp1),gs5(imp1), 
     .  hs1(imp1,jmp1),hs2(imp1,jmp1),hs3(imp1,jmp1),
     .  hs4(imp1,jmp1),hs5(imp1,jmp1)
c
      dimension    flux2_mn(mx1dwk,mx1dwk), flux3_mn(mx1dwk,mx1dwk),
     .             flux4_mn(mx1dwk,mx1dwk), flux5_mn(mx1dwk,mx1dwk),
     .             flux2_mx(mx1dwk,mx1dwk), flux3_mx(mx1dwk,mx1dwk),
     .             flux4_mx(mx1dwk,mx1dwk), flux5_mx(mx1dwk,mx1dwk)
c
      dimension    imap(msegt,msegn),
     .             tauwfi1(jmp1,kmp1)  , tauwfi2(jmp1,kmp1),
     .             tauwfj1(imp1,kmp1)  , tauwfj2(imp1,kmp1),
     .             tauwfk1(imp1,jmp1)  , tauwfk2(imp1,jmp1)
c
c
      amsgore   = sqrt(gamma)*rm/rey
      gogm1pr   = gamma/(gm1*pr)
      threei    = 1./3. 
c
      imnp  = imn + 1
      jmnp  = jmn + 1
      kmnp  = kmn + 1
      imxp  = imx + 1
      jmxp  = jmx + 1
      kmxp  = kmx + 1
      imxm1 = imx - 1
      jmxm1 = jmx - 1
      kmxm1 = kmx - 1
c 
c--------- compute velocities and molecular viscosity
c
      do 10 k=kmn,kmp1
      do 10 j=jmn,jmp1
      do  5 i=imn,imp1
      up(i,j,k) =  w(i,j,k,2)/w(i,j,k,1)
      vp(i,j,k) =  w(i,j,k,3)/w(i,j,k,1)
      wp(i,j,k) =  w(i,j,k,4)/w(i,j,k,1)
      ts(i,j,k) =  p(i,j,k)/w(i,j,k,1) 
      amu(i,j,k)= ts(i,j,k) 
    5 continue
      if(suthc.gt.0.) then
c
        visc1     = suthc/tref
        do 6 i=imn,imp1
        amu(i,j,k)=((1.+visc1)/(ts(i,j,k)+visc1))*ts(i,j,k)**1.5
    6   continue
      endif
c
   10 continue
c
c
c     initialize viscous fluxes
c
      do 11 nv=1,4
      do 11 k=1,kmp1
      do 11 j=1,jmp1
      do 11 i=1,imp1
      fv(i,j,k,nv) = 0.
 11   continue
c
c     initialize temporary flux arrays 
c
      do j=1,mx1dwk
        do i=1,mx1dwk
        flux2_mn(i,j) = 0.
        flux3_mn(i,j) = 0.
        flux4_mn(i,j) = 0.
        flux5_mn(i,j) = 0.
        flux2_mx(i,j) = 0.
        flux3_mx(i,j) = 0.
        flux4_mx(i,j) = 0.
        flux5_mx(i,j) = 0.
        enddo
      enddo
c
      if (ivisb(1).eq.0.and.ivisb(2).eq.0.and.ivisb(3).eq.0) return
c 
c********* begin evaluation of viscous terms in k (zeta) direction
c          k=kmn plane
c 
      if (ivisb(3).eq.0) go to 121
c
      do 110 iseg = 1,nseg
c
      nface    =  imap(2 ,iseg)
      if(nface .ne. 5) go to 110
      nbctype =  imap(1 ,iseg)
      ibeg    =  imap(3 ,iseg) +1
      iend    =  imap(4 ,iseg)
      jbeg    =  imap(5 ,iseg) +1
      jend    =  imap(6 ,iseg)
c
      
      if(nbctype.eq.102) then
c
c       use tau from wall function to construct fluxes
c
        do j=jbeg,jend
          do i=ibeg,iend
          flux2_mn(i,j) =  0. 
          flux3_mn(i,j) =  0.
          flux4_mn(i,j) =  0.
          flux5_mn(i,j) =  0.
          enddo
        enddo
c
      else
c
c       standard viscous fluxes
c
        do j=jbeg,jend
          do i=ibeg,iend
          b13       = skx(i,j,kmn)
          b23       = sky(i,j,kmn)
          b33       = skz(i,j,kmn)
c 
          ua        = .5*(up(i,j,kmnp) + up(i,j,kmn))
          va        = .5*(vp(i,j,kmnp) + vp(i,j,kmn))
          wa        = .5*(wp(i,j,kmnp) + wp(i,j,kmn))
          uzeta     =     up(i,j,kmnp) - up(i,j,kmn)
          vzeta     =     vp(i,j,kmnp) - vp(i,j,kmn)
          wzeta     =     wp(i,j,kmnp) - wp(i,j,kmn)
          tzeta     =     ts(i,j,kmnp) - ts(i,j,kmn)
c 
          amuoj     = (amu(i,j,kmnp) + amu(i,j,kmn))/(vol(i,j,kmnp)+ 
     .                                              vol(i,j,kmn ))
          eps       = .5*(eomu(i,j,kmnp) + eomu(i,j,kmn))
          amuojb    = amuoj*(1.+eps)
          amuojh    = amuoj*(1.+eps*pr/prt) 
c 
          tauki     = amuojb*((4.*threei*b13*b13+b23*b23+b33*b33)*uzeta+
     .                         threei*b13*(b23*vzeta+b33*wzeta))
          taukj     = amuojb*((b13*b13+4.*threei*b23*b23+b33*b33)*vzeta+
     .                         threei*b23*(b13*uzeta+b33*wzeta))
          taukk     = amuojb*((b13*uzeta+b23*vzeta)*threei*b33+wzeta* 
     .                        (b13*b13+b23*b23+4.*threei*b33*b33)) 
          hs2(i,j)  = tauki 
          hs3(i,j)  = taukj 
          hs4(i,j)  = taukk 
          hs5(i,j)  = ua*tauki+va*taukj+wa*taukk+gogm1pr*amuojh*
     .                (b13*b13+b23*b23+b33*b33)*tzeta 
          flux2_mn(i,j) = hs2(i,j)
          flux3_mn(i,j) = hs3(i,j)
          flux4_mn(i,j) = hs4(i,j)
          flux5_mn(i,j) = hs5(i,j)
          enddo
        enddo
c       end of viscous fluxes at kmn
      endif
  110 continue
c
c     begin flux evaluation at k=kmx plane
      do 120 iseg = 1,nseg
c
      nface    =  imap(2 ,iseg)
      if(nface .ne. 6) go to 120
      nbctype =  imap(1 ,iseg)
      ibeg    =  imap(3 ,iseg) +1
      iend    =  imap(4 ,iseg)
      jbeg    =  imap(5 ,iseg) +1
      jend    =  imap(6 ,iseg)
c
      
      if(nbctype.eq.102) then
c
c       use tau from wall function to construct fluxes
c
        do j=jbeg,jend
          do i=ibeg,iend
          flux2_mx(i,j) =  0. 
          flux3_mx(i,j) =  0.
          flux4_mx(i,j) =  0.
          flux5_mx(i,j) =  0.
          enddo
        enddo
c
      else
c
c       standard viscous fluxes
c
        do j=jbeg,jend
          do i=ibeg,iend
          b13       = skx(i,j,kmx)
          b23       = sky(i,j,kmx)
          b33       = skz(i,j,kmx)
c 
          ua        = .5*(up(i,j,kmxp) + up(i,j,kmx))
          va        = .5*(vp(i,j,kmxp) + vp(i,j,kmx))
          wa        = .5*(wp(i,j,kmxp) + wp(i,j,kmx))
          uzeta     =     up(i,j,kmxp) - up(i,j,kmx)
          vzeta     =     vp(i,j,kmxp) - vp(i,j,kmx)
          wzeta     =     wp(i,j,kmxp) - wp(i,j,kmx)
          tzeta     =     ts(i,j,kmxp) - ts(i,j,kmx)
c 
          amuoj     = (amu(i,j,kmxp) + amu(i,j,kmx))/(vol(i,j,kmxp)+ 
     .                                              vol(i,j,kmx ))
          eps       = .5*(eomu(i,j,kmxp) + eomu(i,j,kmx))
          amuojb    = amuoj*(1.+eps)
          amuojh    = amuoj*(1.+eps*pr/prt) 
c 
          tauki     = amuojb*((4.*threei*b13*b13+b23*b23+b33*b33)*uzeta+
     .                         threei*b13*(b23*vzeta+b33*wzeta))
          taukj     = amuojb*((b13*b13+4.*threei*b23*b23+b33*b33)*vzeta+
     .                         threei*b23*(b13*uzeta+b33*wzeta))
          taukk     = amuojb*((b13*uzeta+b23*vzeta)*threei*b33+wzeta* 
     .                        (b13*b13+b23*b23+4.*threei*b33*b33)) 
          hs2(i,j)  = tauki 
          hs3(i,j)  = taukj 
          hs4(i,j)  = taukk 
          hs5(i,j)  = ua*tauki+va*taukj+wa*taukk+gogm1pr*amuojh*
     .                (b13*b13+b23*b23+b33*b33)*tzeta 
          flux2_mx(i,j) = hs2(i,j)
          flux3_mx(i,j) = hs3(i,j)
          flux4_mx(i,j) = hs4(i,j)
          flux5_mx(i,j) = hs5(i,j)
          enddo
        enddo
c       end of viscous fluxes at kmx
      endif
  120 continue
c
c     flux computation at kmin and kmax boundary faces is now complete
c
c     compute viscous fluxes in k-direction in the interior
c
c     First store the fluxes at kmn in hs arrays
c
      do j=jmn+1,jmx
        do i=imn+1,imx
        hs2(i,j) = flux2_mn(i,j)
        hs3(i,j) = flux3_mn(i,j)
        hs4(i,j) = flux4_mn(i,j)
        hs5(i,j) = flux5_mn(i,j)
        enddo
      enddo
c
c     start the k-loop
c
      do k=kmn+1,kmx-1
c       
        do j=jmn+1,jmx
          do i=imn+1,imx
          hs2km1    = hs2(i,j)
          hs3km1    = hs3(i,j)
          hs4km1    = hs4(i,j)
          hs5km1    = hs5(i,j)
c
          b13       = skx(i,j,k)
          b23       = sky(i,j,k)
          b33       = skz(i,j,k)
c 
          ua        = .5*(up(i,j,k+1) + up(i,j,k))
          va        = .5*(vp(i,j,k+1) + vp(i,j,k))
          wa        = .5*(wp(i,j,k+1) + wp(i,j,k))
          uzeta     =     up(i,j,k+1) - up(i,j,k)
          vzeta     =     vp(i,j,k+1) - vp(i,j,k)
          wzeta     =     wp(i,j,k+1) - wp(i,j,k)
          tzeta     =     ts(i,j,k+1) - ts(i,j,k)
c 
          amuoj     = (amu(i,j,k+1) + amu(i,j,k))/(vol(i,j,k+1)+ 
     .                                              vol(i,j,k ))
          eps       = .5*(eomu(i,j,k+1) + eomu(i,j,k))
          amuojb    = amuoj*(1.+eps)
          amuojh    = amuoj*(1.+eps*pr/prt) 
c 
          tauki     = amuojb*((4.*threei*b13*b13+b23*b23+b33*b33)*uzeta+
     .                         threei*b13*(b23*vzeta+b33*wzeta))
          taukj     = amuojb*((b13*b13+4.*threei*b23*b23+b33*b33)*vzeta+
     .                         threei*b23*(b13*uzeta+b33*wzeta))
          taukk     = amuojb*((b13*uzeta+b23*vzeta)*threei*b33+wzeta* 
     .                        (b13*b13+b23*b23+4.*threei*b33*b33)) 
          hs2(i,j)  = tauki 
          hs3(i,j)  = taukj 
          hs4(i,j)  = taukk 
          hs5(i,j)  = ua*tauki+va*taukj+wa*taukk+gogm1pr*amuojh*
     .                (b13*b13+b23*b23+b33*b33)*tzeta 
c
c--------- accumulate viscous-fluxes in zeta-direction
c
c
          fv(i,j,k,1) =  fv(i,j,k,1) - amsgore*(hs2(i,j)-hs2km1)
          fv(i,j,k,2) =  fv(i,j,k,2) - amsgore*(hs3(i,j)-hs3km1)
          fv(i,j,k,3) =  fv(i,j,k,3) - amsgore*(hs4(i,j)-hs4km1)
          fv(i,j,k,4) =  fv(i,j,k,4) - amsgore*(hs5(i,j)-hs5km1)
c
          enddo
        enddo
c       
      enddo
c
c       add the contribution from kmx face
c
        k = kmx
c
        do j=jmn+1,jmx
          do i=imn+1,imx
          fv(i,j,k,1) =  fv(i,j,k,1) - amsgore*(flux2_mx(i,j) -hs2(i,j))
          fv(i,j,k,2) =  fv(i,j,k,2) - amsgore*(flux3_mx(i,j) -hs3(i,j))
          fv(i,j,k,3) =  fv(i,j,k,3) - amsgore*(flux4_mx(i,j) -hs4(i,j))
          fv(i,j,k,4) =  fv(i,j,k,4) - amsgore*(flux5_mx(i,j) -hs5(i,j))
          enddo
        enddo
c
  121 continue
c 
c********* begin evaluation of viscous terms in j (eta) direction
c          j=jmn plane 
c 
      if (ivisb(2).eq.0) go to 221
c
c     initialize flux arrays for complete i-k range
c     before starting the computations on segments
c
        do k=kmnp,kmx
          do i=imnp,imx
          flux2_mn(i,k) =  0. 
          flux3_mn(i,k) =  0.
          flux4_mn(i,k) =  0.
          flux5_mn(i,k) =  0.
          enddo
        enddo
c
      do 210 iseg = 1,nseg
c
      nface    =  imap(2 ,iseg)
      if(nface .ne. 3) go to 210
      nbctype =  imap(1 ,iseg)
      kbeg    =  imap(3 ,iseg) +1
      kend    =  imap(4 ,iseg)
      ibeg    =  imap(5 ,iseg) +1
      iend    =  imap(6 ,iseg)
c
      
      if(nbctype.eq.102) then
c
c       use tau from wall function to construct fluxes
c
        do k=kbeg,kend
          do i=ibeg,iend
c         uslip   = w(i,jmn,k,2)/w(i,jmn,k,1)
c         vslip   = w(i,jmn,k,3)/w(i,jmn,k,1)
c         wslip   = w(i,jmn,k,4)/w(i,jmn,k,1)

cvn       uslip   = w(i,jmnp,k,2)/w(i,jmnp,k,1)
cvn       vslip   = w(i,jmnp,k,3)/w(i,jmnp,k,1)
cvn       wslip   = w(i,jmnp,k,4)/w(i,jmnp,k,1)

          uslip   = (w(i,jmnp,k,2) + w(i,jmn,k,2))/(w(i,jmnp,k,1) +
     .                                              w(i,jmn ,k,1) )
          vslip   = (w(i,jmnp,k,3) + w(i,jmn,k,3))/(w(i,jmnp,k,1) +
     .                                              w(i,jmn ,k,1) )
          wslip   = (w(i,jmnp,k,4) + w(i,jmn,k,4))/(w(i,jmnp,k,1) +
     .                                              w(i,jmn ,k,1) )
          vtslip  = sqrt(uslip**2 + vslip**2 + wslip**2)
c         vtslip  = max (1.e-10,vtslip)
c
          b12       = sjx(i,jmn,k)
          b22       = sjy(i,jmn,k)
          b32       = sjz(i,jmn,k)
c
c         surf      = (b12**2+b22**2+b32**2)/vol(i,jmnp,k)
          surf      = sqrt(b12**2+b22**2+b32**2)
          etauw   = tauwfj1(i,k)*surf/vtslip
          flux2_mn(i,k) =  uslip*etauw
          flux3_mn(i,k) =  vslip*etauw
          flux4_mn(i,k) =  wslip*etauw
c
c         heat flux
c
          teta      =  ts(i,jmnp,k) - ts(i,jmn,k)
c
c         set heat transfer = 0 for adiabatic wall
          teta      =  0.
c
          amuoj     = (amu(i,jmnp,k)+amu(i,jmn,k))/(vol(i,jmnp,k)+
     .                                              vol(i,jmn ,k))
          eps       = .5*(eomu(i,jmnp,k)  +eomu(i,jmn,k))
          amuojb    = amuoj*(1.+eps)
          amuojh    = amuoj*(1.+eps*pr/prt)
          flux5_mn(i,k) =  uslip*flux2_mn(i,k) + vslip*flux3_mn(i,k)
     .                  +  wslip*flux4_mn(i,k)
     .                  +  gogm1pr*amuojh*
     .                    (b12*b12+ b22*b22+b32*b32)*teta
c
c         set flux for energy to zero
cvn       flux5_mn(i,k) = 0.
          enddo
        enddo
c
      else
c
c       standard viscous fluxes
c
        do k=kbeg,kend
          do i=ibeg,iend
          b12       = sjx(i,jmn,k)
          b22       = sjy(i,jmn,k)
          b32       = sjz(i,jmn,k)
c
          ua        = .5*(up(i,jmnp,k) + up(i,jmn,k))
          va        = .5*(vp(i,jmnp,k) + vp(i,jmn,k))
          wa        = .5*(wp(i,jmnp,k) + wp(i,jmn,k))
          ueta      =     up(i,jmnp,k) - up(i,jmn,k)
          veta      =     vp(i,jmnp,k) - vp(i,jmn,k)
          weta      =     wp(i,jmnp,k) - wp(i,jmn,k)
          teta      =     ts(i,jmnp,k) - ts(i,jmn,k)
c
          amuoj     = (amu(i,jmnp,k)+amu(i,jmn,k))/(vol(i,jmnp,k)+
     .                                              vol(i,jmn ,k))
          eps       = .5*(eomu(i,jmnp,k)  +eomu(i,jmn,k))
          amuojb    = amuoj*(1.+eps)
          amuojh    = amuoj*(1.+eps*pr/prt)
c
          tauji     = amuojb*((4.*threei*b12*b12+b22*b22+b32*b32)*ueta+
     .                         threei*b12*(b22*veta+b32*weta))
          taujj     = amuojb*(threei*b22*(b12*ueta+b32*weta)+(b12*b12+
     .                        4.*threei*b22*b22+b32*b32)*veta)
          taujk     = amuojb*(threei*b32*(b12*ueta+b22*veta)+
     .                        (b12*b12+b22*b22+4.*threei*b32*b32)*weta)
          gs2(i)    = tauji
          gs3(i)    = taujj
          gs4(i)    = taujk
          gs5(i)    = ua*tauji+va*taujj+wa*taujk+gogm1pr*amuojh*(b12*b12
     .              + b22*b22+b32*b32)*teta
          flux2_mn(i,k) = gs2(i)
          flux3_mn(i,k) = gs3(i)
          flux4_mn(i,k) = gs4(i)
          flux5_mn(i,k) = gs5(i)
          enddo
        enddo
c       end of viscous fluxes at jmn
      endif
  210 continue
c
c     begin flux evaluation at jmx
      do 220 iseg = 1,nseg
c
      nface    =  imap(2 ,iseg)
      if(nface .ne. 4) go to 220
      nbctype =  imap(1 ,iseg)
      kbeg    =  imap(3 ,iseg) +1
      kend    =  imap(4 ,iseg)
      ibeg    =  imap(5 ,iseg) +1
      iend    =  imap(6 ,iseg)
c
      
      if(nbctype.eq.102) then
c
c       use tau from wall function to construct fluxes
c
c
        do k=kbeg,kend
          do i=ibeg,iend
          uslip   = w(i,jmxp,k,2)/w(i,jmxp,k,1)
          vslip   = w(i,jmxp,k,3)/w(i,jmxp,k,1)
          wslip   = w(i,jmxp,k,4)/w(i,jmxp,k,1)
          vtslip  = sqrt(uslip**2 + vslip**2 + wslip**2)
c         vtslip  = max (1.e-10,vtslip)
c
          b12       = sjx(i,jmx,k)
          b22       = sjy(i,jmx,k)
          b32       = sjz(i,jmx,k)
c
          surf      = sqrt(b12**2+b22**2+b32**2)
          etauw   = tauwfj2(i,k)*surf/vtslip
          flux2_mx(i,k) =  uslip*etauw
          flux3_mx(i,k) =  vslip*etauw
          flux4_mx(i,k) =  wslip*etauw
c
c         heat flux
c
          teta      =  ts(i,jmxp,k) - ts(i,jmx,k)
c
c         set heat transfer = 0 for adiabatic wall
          teta      =  0.
c
          amuoj     = (amu(i,jmxp,k)+amu(i,jmx,k))/(vol(i,jmxp,k)+
     .                                              vol(i,jmx ,k))
          eps       = .5*(eomu(i,jmxp,k)  +eomu(i,jmx,k))
          amuojb    = amuoj*(1.+eps)
          amuojh    = amuoj*(1.+eps*pr/prt)
          flux5_mx(i,k) =  uslip*flux2_mx(i,k) + vslip*flux3_mx(i,k)
     .                  +  wslip*flux4_mx(i,k)
     .                  +  gogm1pr*amuojh*
     .                    (b12*b12+ b22*b22+b32*b32)*teta
c
c         set flux for energy to zero
          flux5_mx(i,k) = 0.
          enddo
        enddo
c
c
      else
c
c       standard viscous fluxes
c
        do k=kbeg,kend
          do i=ibeg,iend
          b12       = sjx(i,jmx,k)
          b22       = sjy(i,jmx,k)
          b32       = sjz(i,jmx,k)
c
          ua        = .5*(up(i,jmxp,k) + up(i,jmx,k))
          va        = .5*(vp(i,jmxp,k) + vp(i,jmx,k))
          wa        = .5*(wp(i,jmxp,k) + wp(i,jmx,k))
          ueta      =     up(i,jmxp,k) - up(i,jmx,k)
          veta      =     vp(i,jmxp,k) - vp(i,jmx,k)
          weta      =     wp(i,jmxp,k) - wp(i,jmx,k)
          teta      =     ts(i,jmxp,k) - ts(i,jmx,k)
c
          amuoj     = (amu(i,jmxp,k)+amu(i,jmx,k))/(vol(i,jmxp,k)+
     .                                              vol(i,jmx ,k))
          eps       = .5*(eomu(i,jmxp,k)  +eomu(i,jmx,k))
          amuojb    = amuoj*(1.+eps)
          amuojh    = amuoj*(1.+eps*pr/prt)
c
          tauji     = amuojb*((4.*threei*b12*b12+b22*b22+b32*b32)*ueta+
     .                         threei*b12*(b22*veta+b32*weta))
          taujj     = amuojb*(threei*b22*(b12*ueta+b32*weta)+(b12*b12+
     .                        4.*threei*b22*b22+b32*b32)*veta)
          taujk     = amuojb*(threei*b32*(b12*ueta+b22*veta)+
     .                        (b12*b12+b22*b22+4.*threei*b32*b32)*weta)
          gs2(i)    = tauji
          gs3(i)    = taujj
          gs4(i)    = taujk
          gs5(i)    = ua*tauji+va*taujj+wa*taujk+gogm1pr*amuojh*(b12*b12
     .              + b22*b22+b32*b32)*teta
          flux2_mx(i,k) = gs2(i)
          flux3_mx(i,k) = gs3(i)
          flux4_mx(i,k) = gs4(i)
          flux5_mx(i,k) = gs5(i)

          enddo
        enddo
c       end of viscous fluxes at jmx
      endif
  220 continue
c
c     flux computation at jmin and jmax boundary faces is now complete
c
c     compute viscous fluxes in j-direction in the interior
c
c     First store the fluxes at jmn in hs arrays
c
      do k=kmn+1,kmx
        do i=imn+1,imx
        hs2(i,k) = flux2_mn(i,k)
        hs3(i,k) = flux3_mn(i,k)
        hs4(i,k) = flux4_mn(i,k)
        hs5(i,k) = flux5_mn(i,k)
        enddo
      enddo
c
c     start the j-loop
c
      do j=jmn+1,jmx-1
c       
        do k=kmn+1,kmx
          do i=imn+1,imx
          hs2jm1    = hs2(i,k)
          hs3jm1    = hs3(i,k)
          hs4jm1    = hs4(i,k)
          hs5jm1    = hs5(i,k)
c
          b12       = sjx(i,j,k)
          b22       = sjy(i,j,k)
          b32       = sjz(i,j,k)
c
          ua        = .5*(up(i,j+1,k) + up(i,j,k))
          va        = .5*(vp(i,j+1,k) + vp(i,j,k))
          wa        = .5*(wp(i,j+1,k) + wp(i,j,k))
          ueta      =     up(i,j+1,k) - up(i,j,k)
          veta      =     vp(i,j+1,k) - vp(i,j,k)
          weta      =     wp(i,j+1,k) - wp(i,j,k)
          teta      =     ts(i,j+1,k) - ts(i,j,k)
c
          amuoj     = (amu(i,j+1,k)+amu(i,j,k))/(vol(i,j+1,k)+
     .                                              vol(i,j ,k))
          eps       = .5*(eomu(i,j+1,k)  +eomu(i,j,k))
          amuojb    = amuoj*(1.+eps)
          amuojh    = amuoj*(1.+eps*pr/prt)
c
          tauji     = amuojb*((4.*threei*b12*b12+b22*b22+b32*b32)*ueta+
     .                         threei*b12*(b22*veta+b32*weta))
          taujj     = amuojb*(threei*b22*(b12*ueta+b32*weta)+(b12*b12+
     .                        4.*threei*b22*b22+b32*b32)*veta)
          taujk     = amuojb*(threei*b32*(b12*ueta+b22*veta)+
     .                        (b12*b12+b22*b22+4.*threei*b32*b32)*weta)
          hs2(i,k)  = tauji
          hs3(i,k)  = taujj
          hs4(i,k)  = taujk
          hs5(i,k)  = ua*tauji+va*taujj+wa*taujk+gogm1pr*amuojh*(b12*b12
     .              + b22*b22+b32*b32)*teta
c
c--------- accumulate viscous-fluxes in eta-direction
c
c
          fv(i,j,k,1) =  fv(i,j,k,1) - amsgore*(hs2(i,k)-hs2jm1)
          fv(i,j,k,2) =  fv(i,j,k,2) - amsgore*(hs3(i,k)-hs3jm1)
          fv(i,j,k,3) =  fv(i,j,k,3) - amsgore*(hs4(i,k)-hs4jm1)
          fv(i,j,k,4) =  fv(i,j,k,4) - amsgore*(hs5(i,k)-hs5jm1)
c
          enddo
        enddo
c       
      enddo
c
c       add the contribution from jmx face
c
        j = jmx
c
        do k=kmn+1,kmx
          do i=imn+1,imx
          fv(i,j,k,1) =  fv(i,j,k,1) - amsgore*(flux2_mx(i,k) -hs2(i,k))
          fv(i,j,k,2) =  fv(i,j,k,2) - amsgore*(flux3_mx(i,k) -hs3(i,k))
          fv(i,j,k,3) =  fv(i,j,k,3) - amsgore*(flux4_mx(i,k) -hs4(i,k))
          fv(i,j,k,4) =  fv(i,j,k,4) - amsgore*(flux5_mx(i,k) -hs5(i,k))
          enddo
        enddo
c
  221 continue
c 
c********* begin evaluation of viscous terms for i (xi) direction
c          i=imn plane 
c 
      if (ivisb(1).eq.0) go to 321
c
c     initialize flux arrays for complete j-k range
c     before starting the computations on segments
c
        do k=kmnp,kmx
          do j=imnp,imx
          flux2_mn(j,k) =  0. 
          flux3_mn(j,k) =  0.
          flux4_mn(j,k) =  0.
          flux5_mn(j,k) =  0.
          enddo
        enddo
c
c
      do 310 iseg = 1,nseg
c
      nface    =  imap(2 ,iseg)
      if(nface .ne. 1) go to 310
      nbctype =  imap(1 ,iseg)
      jbeg    =  imap(3 ,iseg) +1
      jend    =  imap(4 ,iseg)
      kbeg    =  imap(5 ,iseg) +1
      kend    =  imap(6 ,iseg)
c
      
      if(nbctype.eq.102) then
c
c       use tau from wall function to construct fluxes
c
        do k=kbeg,kend
          do j=jbeg,jend
          flux2_mn(j,k) =  0. 
          flux3_mn(j,k) =  0.
          flux4_mn(j,k) =  0.
          flux5_mn(j,k) =  0.
          enddo
        enddo
c
      else
c
c       standard viscous fluxes
c
        i = imn
        do k=kbeg,kend
          do j=jbeg,jend
          b11     = six(i,j,k)
          b21     = siy(i,j,k)
          b31     = siz(i,j,k)
c
          ua      = .5*(up(i+1,j,k) + up(i,j,k))
          va      = .5*(vp(i+1,j,k) + vp(i,j,k))
          wa      = .5*(wp(i+1,j,k) + wp(i,j,k))
          uxi     = up(i+1,j,k) - up(i,j,k)
          vxi     = vp(i+1,j,k) - vp(i,j,k)
          wxi     = wp(i+1,j,k) - wp(i,j,k)
          txi     = ts(i+1,j,k) - ts(i,j,k)
c
          amuoj   = (amu(i,j,k)+amu(i+1,j,k))/(vol(i,j,k)+vol(i+1,j,k))
          eps     = .5*(eomu(i,j,k)+eomu(i+1,j,k))
          amuojb  = amuoj*(1.+eps)
          amuojh  = amuoj*(1.+eps*pr/prt)
c
          tauii     = amuojb*((4.*threei*b11*b11+b21*b21+b31*b31)*uxi+
     .                         threei*b11*(b21*vxi+b31*wxi))
          tauij     = amuojb*(threei*b21*(b11*uxi+b31*wxi)+(b11*b11+
     .                        4.*threei*b21*b21+b31*b31)*vxi)
          tauik     = amuojb*(threei*b31*(b11*uxi+b21*vxi)+(b11*b11+
     .                        b21*b21+4.*threei*b31*b31)*wxi)
c
          flux2_mn(j,k) = tauii
          flux3_mn(j,k) = tauij
          flux4_mn(j,k) = tauik
          flux5_mn(j,k) = ua*tauii+va*tauij+wa*tauik +gogm1pr*amuojh*
     .                   (b11*b11+b21*b21+b31*b31)*txi
          enddo
        enddo
c       end of viscous fluxes at imn face
      endif
  310 continue
c
c     begin flux evaluation at jmx face
      do 320 iseg = 1,nseg
c
      nface    =  imap(2 ,iseg)
      if(nface .ne. 2) go to 320
      nbctype =  imap(1 ,iseg)
      jbeg    =  imap(3 ,iseg) +1
      jend    =  imap(4 ,iseg)
      kbeg    =  imap(5 ,iseg) +1
      kend    =  imap(6 ,iseg)
c
      
      if(nbctype.eq.102) then
c
c       use tau from wall function to construct fluxes
c
        do k=kbeg,kend
          do j=jbeg,jend
          flux2_mx(j,k) =  0. 
          flux3_mx(j,k) =  0.
          flux4_mx(j,k) =  0.
          flux5_mx(j,k) =  0.
          enddo
        enddo
c
      else
c
c       standard viscous fluxes
c
        i = imx
        do k=kbeg,kend
          do j=jbeg,jend
          b11     = six(i,j,k)
          b21     = siy(i,j,k)
          b31     = siz(i,j,k)
c
          ua      = .5*(up(i+1,j,k) + up(i,j,k))
          va      = .5*(vp(i+1,j,k) + vp(i,j,k))
          wa      = .5*(wp(i+1,j,k) + wp(i,j,k))
          uxi     = up(i+1,j,k) - up(i,j,k)
          vxi     = vp(i+1,j,k) - vp(i,j,k)
          wxi     = wp(i+1,j,k) - wp(i,j,k)
          txi     = ts(i+1,j,k) - ts(i,j,k)
c
          amuoj   = (amu(i,j,k)+amu(i+1,j,k))/(vol(i,j,k)+vol(i+1,j,k))
          eps     = .5*(eomu(i,j,k)+eomu(i+1,j,k))
          amuojb  = amuoj*(1.+eps)
          amuojh  = amuoj*(1.+eps*pr/prt)
c
          tauii     = amuojb*((4.*threei*b11*b11+b21*b21+b31*b31)*uxi+
     .                         threei*b11*(b21*vxi+b31*wxi))
          tauij     = amuojb*(threei*b21*(b11*uxi+b31*wxi)+(b11*b11+
     .                        4.*threei*b21*b21+b31*b31)*vxi)
          tauik     = amuojb*(threei*b31*(b11*uxi+b21*vxi)+(b11*b11+
     .                        b21*b21+4.*threei*b31*b31)*wxi)
c
          flux2_mx(j,k) = tauii
          flux3_mx(j,k) = tauij
          flux4_mx(j,k) = tauik
          flux5_mx(j,k) = ua*tauii+va*tauij+wa*tauik +gogm1pr*amuojh*
     .                   (b11*b11+b21*b21+b31*b31)*txi
          enddo
        enddo
c       end of viscous fluxes at imx
      endif
  320 continue
c
c     flux computation at imn and imx boundary faces is now complete
c
c     compute viscous fluxes in i-direction in the interior
c
c     First store the fluxes at imn in hs arrays
c
      do k=kmn+1,kmx
        do j=jmn+1,jmx
        hs2(j,k) = flux2_mn(j,k)
        hs3(j,k) = flux3_mn(j,k)
        hs4(j,k) = flux4_mn(j,k)
        hs5(j,k) = flux5_mn(j,k)
        enddo
      enddo
c
c     start the i-loop
c
      do i=imn+1,imx-1
c       
        do k=kmn+1,kmx
          do j=jmn+1,jmx
          hs2im1    = hs2(j,k)
          hs3im1    = hs3(j,k)
          hs4im1    = hs4(j,k)
          hs5im1    = hs5(j,k)
c
          b11     = six(i,j,k)
          b21     = siy(i,j,k)
          b31     = siz(i,j,k)
c
          ua      = .5*(up(i+1,j,k) + up(i,j,k))
          va      = .5*(vp(i+1,j,k) + vp(i,j,k))
          wa      = .5*(wp(i+1,j,k) + wp(i,j,k))
          uxi     = up(i+1,j,k) - up(i,j,k)
          vxi     = vp(i+1,j,k) - vp(i,j,k)
          wxi     = wp(i+1,j,k) - wp(i,j,k)
          txi     = ts(i+1,j,k) - ts(i,j,k)
c
          amuoj   = (amu(i,j,k)+amu(i+1,j,k))/(vol(i,j,k)+vol(i+1,j,k))
          eps     = .5*(eomu(i,j,k)+eomu(i+1,j,k))
          amuojb  = amuoj*(1.+eps)
          amuojh  = amuoj*(1.+eps*pr/prt)
c
          tauii     = amuojb*((4.*threei*b11*b11+b21*b21+b31*b31)*uxi+
     .                         threei*b11*(b21*vxi+b31*wxi))
          tauij     = amuojb*(threei*b21*(b11*uxi+b31*wxi)+(b11*b11+
     .                        4.*threei*b21*b21+b31*b31)*vxi)
          tauik     = amuojb*(threei*b31*(b11*uxi+b21*vxi)+(b11*b11+
     .                        b21*b21+4.*threei*b31*b31)*wxi)
c
          hs2(j,k)  = tauii
          hs3(j,k)  = tauij
          hs4(j,k)  = tauik
          hs5(j,k)  = ua*tauii+va*tauij+wa*tauik +gogm1pr*amuojh*
     .               (b11*b11+b21*b21+b31*b31)*txi
c
c--------- accumulate viscous-fluxes in xi-direction
c
c
          fv(i,j,k,1) =  fv(i,j,k,1) - amsgore*(hs2(j,k)-hs2im1)
          fv(i,j,k,2) =  fv(i,j,k,2) - amsgore*(hs3(j,k)-hs3im1)
          fv(i,j,k,3) =  fv(i,j,k,3) - amsgore*(hs4(j,k)-hs4im1)
          fv(i,j,k,4) =  fv(i,j,k,4) - amsgore*(hs5(j,k)-hs5im1)
c
          enddo
        enddo
c       
      enddo
c
c       add the contribution from imx face
c
        i = imx
c
        do k=kmn+1,kmx
          do j=jmn+1,jmx
          fv(i,j,k,1) =  fv(i,j,k,1) - amsgore*(flux2_mx(j,k) -hs2(j,k))
          fv(i,j,k,2) =  fv(i,j,k,2) - amsgore*(flux3_mx(j,k) -hs3(j,k))
          fv(i,j,k,3) =  fv(i,j,k,3) - amsgore*(flux4_mx(j,k) -hs4(j,k))
          fv(i,j,k,4) =  fv(i,j,k,4) - amsgore*(flux5_mx(j,k) -hs5(j,k))
          enddo
        enddo
c
  321 continue
c
c
      return
      end 
