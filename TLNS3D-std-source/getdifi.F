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
      subroutine getdifi(imn,jmn,kmn,im,jm,km,
     .                   imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .                   m1cc,mgrlev,igrid,mbloc,ibloc,
     .                   mxsizc,mx1dwk,
     .                   imap,msegt,msegn,nseg,p,
     .                   imx,jmx,kmx,imxp1,jmxp1,kmxp1,
     .                   imxp2,jmxp2,kmxp2,wn,pn,fv,
     .                   wk2da,wk2db,wk2dc,wk2dd,wk3d    )
c=======================================================================
c
c     object : To get the dp & dw terms needed for the dissipation
c              in the i-direction.
c
c     initial coding : by m.d.sanetrik (Sept, 1990)
c     latest mods    : by v.n.vatsa to accommodate patching (Nov. 94)
c                      by v.n.vatsa for preconditioning (Feb. 1998)
c
 
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      common/rkdis/ vis0,vis2,vis4,zeta,vepsn,vepsl,enteps,icau
c
c
c-------   preconditioning related quantities      -----------------
c
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubcon2,iprec
c
c
      dimension     im  (mgrlev,mbloc), jm  (mgrlev,mbloc),
     .              km  (mgrlev,mbloc), imp1(mgrlev,mbloc),
     .              jmp1(mgrlev,mbloc), kmp1(mgrlev,mbloc),
     .              imp2(mgrlev,mbloc), jmp2(mgrlev,mbloc),
     .              kmp2(mgrlev,mbloc),
     .              imap(msegt,msegn)
c
      dimension     m1cc(mgrlev,mbloc)
c
      dimension     p(mxsizc )
c
      dimension    wn(imxp2,jmxp2,kmxp2,5),pn(imxp2,jmxp2,kmxp2),
     .             fv(imxp2,jmxp2,kmxp2)
c
c
c     array declaration for temporary work-space
c
c     2-d & 3-d arrays
c
      dimension   wk2da(1),wk2db(1),wk2dc(1),wk2dd(1)
      dimension   wk3d(imxp2,jmxp2,kmxp2,5)
c
c
c   calculate dw's & dp in the i-dir for interior points
c   (note: dw's are stored in wk3d, and dp's are stored in fv)
c
      do 20 k=kmn+1,kmx
      do 20 j=jmn+1,jmx
c
      if (iprec.le.1) then
c
        do i=imn-1,imxp1
          wk3d(i,j,k,1) = wn(i+1,j,k,1) -wn(i,j,k,1)
          wk3d(i,j,k,2) = wn(i+1,j,k,2) -wn(i,j,k,2)
          wk3d(i,j,k,3) = wn(i+1,j,k,3) -wn(i,j,k,3)
          wk3d(i,j,k,4) = wn(i+1,j,k,4) -wn(i,j,k,4)
          wk3d(i,j,k,5) = wn(i+1,j,k,5) -wn(i,j,k,5)
     .                  + (1-icau)*(pn(i+1,j,k) -pn(i,j,k))
        enddo
c
      else
c
c         use first derivative of "p,u,v,w,T" variables
c
        do i=imn-1,imxp1
          rhoip         = 1./wn(i+1,j,k,1)
          rhoi          = 1./wn(i  ,j,k,1)
          wk3d(i,j,k,1) = pn(i+1,j,k  )       - pn(i,j,k  )
          wk3d(i,j,k,2) = wn(i+1,j,k,2)*rhoip - wn(i,j,k,2)*rhoi
          wk3d(i,j,k,3) = wn(i+1,j,k,3)*rhoip - wn(i,j,k,3)*rhoi
          wk3d(i,j,k,4) = wn(i+1,j,k,4)*rhoip - wn(i,j,k,4)*rhoi
          wk3d(i,j,k,5) = pn(i+1,j,k  )*rhoip - pn(i,j,k  )*rhoi
        enddo
c
      endif
c
      do i=imn,imx+1
        dplus     = pn(i+1,j,k)-pn(i,j,k)
        dminus    = pn(i,j,k)-pn(i-1,j,k)
        top       = abs(dplus - dminus)
        b1        = abs(pn(i+1,j,k) + 2.*pn(i,j,k) + pn(i-1,j,k))
        b2        = abs(dplus) + abs(dminus)
        bot       = enteps*b1 + (1.-enteps)*b2
        fv(i,j,k) = top/bot
      enddo
   20 continue
c
c
c----------  begin outer loop on the segments  ----------------------------
c
      do 100 iseg = 1,nseg
c
      nface    =  imap(2 ,iseg)
c
      if(nface .eq. 1) then
        nbctype =  imap(1 ,iseg)
        jbeg    =  imap(3 ,iseg) +1
        jend    =  imap(4 ,iseg)
        kbeg    =  imap(5 ,iseg) +1
        kend    =  imap(6 ,iseg)
c
       if(nbctype .le. 1) then
c
        do k=kbeg,kend
        do j=jbeg,jend
          fv(imn-1,j,k)   =   2.*fv(imn,j,k) -fv(imn+1,j,k)
        enddo
        enddo
c
c
        else if(nbctype.gt.1.and.nbctype.le.4) then
          do k=kbeg,kend
          do j=jbeg,jend
            fv(  imn+1,j,k)   =       fv(imn+2,j,k)
            fv(  imn  ,j,k)   =       0.0
            fv(  imn-1,j,k)   =       0.0
          enddo
          enddo
c
        else if(nbctype.gt.4) then
          do k=kbeg,kend
          do j=jbeg,jend
            fv(  imn  ,j,k)   =   fv(imn+1,j,k)
            fv(  imn-1,j,k)   =   fv(imn  ,j,k)
          enddo
          enddo
        end if
c
      else if(nface .eq. 2) then
        nbctype  =  imap(1 ,iseg)
        jbeg     =  imap(3 ,iseg) +1
        jend     =  imap(4 ,iseg)
        kbeg     =  imap(5 ,iseg) +1
        kend     =  imap(6 ,iseg)
c
       if(nbctype .le. 1) then
c
        do k=kbeg,kend
        do j=jbeg,jend
          fv(  imxp2,j,k)   =   2.*fv(imxp1,j,k) - fv(imx,j,k)
        enddo
        enddo
c
        else if(nbctype.gt.1.and.nbctype.le.4) then
          do k=kbeg,kend
          do j=jbeg,jend
            fv(  imx  ,j,k)   =       fv(imx-1,j,k)
            fv(  imxp1,j,k)   =       0.0
            fv(  imxp2,j,k)   =       0.0
          enddo
          enddo
c
        else if(nbctype.gt.4) then
          do k=kbeg,kend
          do j=jbeg,jend
            fv(  imxp1,j,k)   =   fv(imx  ,j,k)
            fv(  imxp2,j,k)   =   fv(imxp1,j,k)
          enddo
          enddo
        end if
      end if
  100 continue
c
c
       return
       end
