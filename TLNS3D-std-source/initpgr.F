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
      subroutine  initpgr(imn,jmn,kmn,imx,jmx,kmx,imp2,jmp2,kmp2,
     .                    vol,npatchf,n1dir,n2dir,n3dir,
     .                    igrid,isoln )
c=======================================================================
c 
c
c     object : to initialize the ghost cell values for flow variables
c              at patched block interfaces
c
c     initial coding : by v.n.vatsa, nasa langley ( Nov., 1994)
c
c

      dimension      vol(imp2,jmp2,kmp2)
c
      dimension      npatchf(6), n1dir(1), n2dir(1),
     .                           n3dir(1)
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      common/facetr/ ifacetr(6),ipatchg
c
c
c      loop over all the block faces
c
       lpatchb = 1
       do 1000 nnface=1,6
c
       lpatche  = lpatchb + npatchf(nnface) - 1
       if(npatchf(nnface).le.0) go to 1000
c
c      convert face numbers to tlns3d's convention
c
       nface = ifacetr(nnface)
c
c      branch-cut for different type of boundary faces
c
      go to (100,200,300,400,500,600) nface
c
c
  100 continue
c--------------------- i = imn  face ----------------------------
c
c
      i   =   imn
c
      do 110 lpatch= lpatchb,lpatche
c
      j   =   n2dir(lpatch)
      k   =   n3dir(lpatch)
c
      vol(i,j,k)  =    0.
  110 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         i   =  imn - 1
c
         do 120 lpatch= lpatchb,lpatche
c
         j   =   n2dir(lpatch)
         k   =   n3dir(lpatch)
c
         vol(i,j,k)  =    0.
  120    continue
      endif
      go to 1000
c
c
  200 continue
c--------------------- i = imx  face ----------------------------
c
c
      i   =   imx + 1
c
      do 210 lpatch= lpatchb,lpatche
c
      j   =   n2dir(lpatch)
      k   =   n3dir(lpatch)
c
      vol(i,j,k)  =    0.
  210 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         i   =  imx + 2
c
         do 220 lpatch= lpatchb,lpatche
c
         j   =   n2dir(lpatch)
         k   =   n3dir(lpatch)
c
         vol(i,j,k)  =    0.
  220    continue
      endif
      go to 1000
c
c
  300 continue
c--------------------- j = jmn  face ----------------------------
c
c
      j   =   jmn
c
      do 310 lpatch= lpatchb,lpatche
c
      k   =   n3dir(lpatch)
      i   =   n1dir(lpatch)
c
      vol(i,j,k)  =    0.
  310 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         j   =  jmn - 1
c
         do 320 lpatch= lpatchb,lpatche
c
         k   =   n3dir(lpatch)
         i   =   n1dir(lpatch)
c
         vol(i,j,k)  =    0.
  320    continue
      endif
      go to 1000
c
c
  400 continue
c--------------------- j = jmx  face ----------------------------
c
c
      j   =   jmx + 1
c
      do 410 lpatch= lpatchb,lpatche
c
      k   =   n3dir(lpatch)
      i   =   n1dir(lpatch)
c
      vol(i,j,k)  =    0.
  410 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         j   =  jmx + 2
c
         do 420 lpatch= lpatchb,lpatche
c
         k   =   n3dir(lpatch)
         i   =   n1dir(lpatch)
c
         vol(i,j,k)  =    0.
  420    continue
      endif
      go to 1000
c
c
  500 continue
c--------------------- k = kmn  face ----------------------------
c
c
      k   =   kmn
c
      do 510 lpatch= lpatchb,lpatche
c
      i   =   n1dir(lpatch)
      j   =   n2dir(lpatch)
c
      vol(i,j,k)  =    0.
  510 continue
c
c
      if( igrid .eq. isoln ) then
c
c        second layer of ghost cells for fine grids
c
         k   =  kmn - 1
c
         do 520 lpatch= lpatchb,lpatche
c
         i   =   n1dir(lpatch)
         j   =   n2dir(lpatch)
c
         vol(i,j,k)  =    0.
  520    continue
      endif
      go to 1000
c
  600 continue
c--------------------- k = kmx  face ----------------------------
c
c
      k   =   kmx + 1
c
      do 610 lpatch= lpatchb,lpatche
c
      i   =   n1dir(lpatch)
      j   =   n2dir(lpatch)
c
      vol(i,j,k)  =    0.
  610 continue
c
      if( igrid .eq. isoln ) then
c
c        second layer of ghost cells for fine grids
c
c
         k   =  kmx + 2
c
         do 620 lpatch= lpatchb,lpatche
c
         i   =   n1dir(lpatch)
         j   =   n2dir(lpatch)
c
         vol(i,j,k)  =    0.
  620    continue
      endif
c
 1000    lpatchb = lpatchb + npatchf(nnface)
c
      return
      end
