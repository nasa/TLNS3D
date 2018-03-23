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
      subroutine  errphrk(imn,jmn,kmn,imx,jmx,kmx,imp2,jmp2,kmp2,
     .                    delw,dtl,vol,dtphy,mx1dwk,errw,ibloc,nodes,
     .                    delwork,rkphysa,rkphyse,mrkstg,nrkstg,irkstg)
c=======================================================================
c 
c     objective : computes temporal error in physical R-K
c                 time-stepping using solutions with different
c                 number of stages
c
c     iniitial coding : v.n.vatsa (Feb. 2001)
c 
c 

      dimension delw(imp2,jmp2,kmp2,5),dtl(imp2,jmp2,kmp2)
c
      dimension errw(mx1dwk),vol(imp2,jmp2,kmp2),
     .          delwork(imp2,jmp2,kmp2,5,nrkstg),
     .          rkphysa(mrkstg,mrkstg),rkphyse(mrkstg)
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
      do 90 k=kmn+1,kmx
      do 90 j=jmn+1,jmx
      iloc   = 0
      do 50 i=imn+1,imx
      iloc   = iloc + 1
         sdir   = 0.
         fac    = dtphy/vol(i,j,k)
         do lstg = 1,irkstg-1
         sdir    = sdir + (rkphysa(irkstg,lstg)-rkphyse(lstg))*
     .                     delwork(i,j,k,nv,lstg)*fac
         enddo
         errw(iloc) = (rkphysa(irkstg,irkstg)-rkphyse(irkstg))*
     .                 delw(i,j,k,nv)*fac + sdir
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
c
      subroutine  errmebdf(imn,jmn,kmn,imx,jmx,kmx,imp2,jmp2,kmp2,
     .                    delw,dtl,vol,dtphy,mx1dwk,errw,ibloc,nodes,
     .                    w,wbar,ntorder)
c=======================================================================
c 
c     objective : computes temporal error in MEBDF
c                 time-stepping using solutions at
c                 different stages
c
c     iniitial coding : v.n.vatsa (Nov. 2003)
c 
c 

      dimension delw(imp2,jmp2,kmp2,5),dtl(imp2,jmp2,kmp2)
c
      dimension errw(mx1dwk),vol(imp2,jmp2,kmp2),
     .          w(imp2,jmp2,kmp2,5),
     .          wbar (imp2,jmp2,kmp2,5,ntorder-1)
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
      do 90 k=kmn+1,kmx
      do 90 j=jmn+1,jmx
      iloc   = 0
      do 50 i=imn+1,imx
      iloc   = iloc + 1
cvn   dt        = cfl*((1.-vt)*dtl(i,j,k) + vt*dtmin)
      dt        = cfl*dtl(i,j,k)
         errw(iloc) = (w(i,j,k,nv) - wbar(i,j,k,nv,1))/dt
cvn      errw(iloc) =  w(i,j,k,nv) - wbar(i,j,k,nv,1)
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
c
      subroutine  errbdf3 (imn,jmn,kmn,imx,jmx,kmx,imp2,jmp2,kmp2,
     .                    delw,dtl,vol,dtphy,mx1dwk,errw,ibloc,nodes,
     .                    w,wold,delwork,ntorder)
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
     .          wold (imp2,jmp2,kmp2,5,ntorder),
     .          delwork(imp2,jmp2,kmp2,5,1)
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
c     c1 = 11./6. - 1.5
c     c2 = -(3.-2.)
c     c3 = 1.5 - .5
c     c4 = -1./3.
      c1 = 11./6. - 1.6933333
      c2 = -(3.-2.58)
      c3 = 1.5 - 1.08
      c4 = -(1./3. -0.1933333)
      do 90 k=kmn+1,kmx
      do 90 j=jmn+1,jmx
      iloc   = 0
      do 50 i=imn+1,imx
      iloc   = iloc + 1
cvn   dt        = cfl*((1.-vt)*dtl(i,j,k) + vt*dtmin)
      dt        = cfl*dtl(i,j,k) 
         errw(iloc) =  c1*w   (i,j,k,nv  ) + c2*wold(i,j,k,nv,1)
     .              +  c3*wold(i,j,k,nv,2) + c4*wold(i,j,k,nv,3)
         errw(iloc) = errw(iloc)/dt
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
c
      subroutine  errbdf2opt (imn,jmn,kmn,imx,jmx,kmx,imp2,jmp2,kmp2,
     .                    delw,dtl,vol,dtphy,mx1dwk,errw,ibloc,nodes,
     .                    w,wold,delwork,ntorder)
c=======================================================================
c 
c     objective : computes temporal error in BDF2OPT
c                 time-stepping using solutions 
c                 with different order accuracies
c
c     iniitial coding : v.n.vatsa (Dec. 2003)
c 
c 

      dimension delw(imp2,jmp2,kmp2,5),dtl(imp2,jmp2,kmp2)
c
      dimension errw(mx1dwk),vol(imp2,jmp2,kmp2),
     .          w(imp2,jmp2,kmp2,5),
     .          wold (imp2,jmp2,kmp2,5,ntorder),
     .          delwork(imp2,jmp2,kmp2,5,1)
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
c     c1 = 2.58/1.6933333 - 18./11.
c     c2 = 9./11. - 1.08/1.6933333
c     c3 = 0.1933333/1.6933333 - 2./11.
c     c4 = 1./1.6933333 - 6./11.
cc    c1 = 11./6. - 1.6933333
cc    c2 = -(3.-2.58)
cc    c3 = 1.5 - 1.08
cc    c4 = -(1./3. -0.1933333)
      c1 = 1.6933333 -1.5
      c2 = -(2.58 -2.)
      c3 = 1.08 - 0.5
      c4 = -0.1933333
      do 90 k=kmn+1,kmx
      do 90 j=jmn+1,jmx
      iloc   = 0
      do 50 i=imn+1,imx
      iloc   = iloc + 1
cvn   dt        = cfl*((1.-vt)*dtl(i,j,k) + vt*dtmin)
      dt        = cfl*dtl(i,j,k) 
         errw(iloc) =  c1*w   (i,j,k,nv  ) + c2*wold(i,j,k,nv,1)
     .              +  c3*wold(i,j,k,nv,2) + c4*wold(i,j,k,nv,3)
         errw(iloc) = errw(iloc)/dt
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
c
      subroutine  errbdf2 (imn,jmn,kmn,imx,jmx,kmx,imp2,jmp2,kmp2,
     .                    delw,dtl,vol,dtphy,mx1dwk,errw,ibloc,nodes,
     .                    w,wold,delwork,ntorder)
c=======================================================================
c 
c     objective : computes temporal error in BDF2
c                 time-stepping using solutions 
c                 with different order accuracies
c
c     iniitial coding : v.n.vatsa (Feb. 2004)
c 
c 

      dimension delw(imp2,jmp2,kmp2,5),dtl(imp2,jmp2,kmp2)
c
      dimension errw(mx1dwk),vol(imp2,jmp2,kmp2),
     .          w(imp2,jmp2,kmp2,5),
     .          wold (imp2,jmp2,kmp2,5,ntorder),
     .          delwork(imp2,jmp2,kmp2,5,1)
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
      c1 = 1.5 -1.
      c2 = -(2.-1.0)
      c3 = 0.5 
      do 90 k=kmn+1,kmx
      do 90 j=jmn+1,jmx
      iloc   = 0
      do 50 i=imn+1,imx
      iloc   = iloc + 1
cvn   dt        = cfl*((1.-vt)*dtl(i,j,k) + vt*dtmin)
      dt        = cfl*dtl(i,j,k) 
         errw(iloc) =  c1*w   (i,j,k,nv  ) + c2*wold(i,j,k,nv,1)
     .              +  c3*wold(i,j,k,nv,2) 
         errw(iloc) = errw(iloc)/dt
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
