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
      subroutine  initpev(imn,jmn,kmn,imx,jmx,kmx,imp2,jmp2,kmp2,
     .                  eomu,turv1,turv2,npatchf,n1dir,n2dir,n3dir,
     .                  igrid,isoln )
c=======================================================================
c 
c
c     object : to initialize the ghost cell values for turbulent
c              quantities at patched block interfaces
c
c     initial coding : by v.n.vatsa, nasa langley ( Nov., 1994)
c
c

      dimension      turv1(imp2,jmp2,kmp2), turv2(imp2,jmp2,kmp2),
     .               eomu (imp2,jmp2,kmp2)
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
      eomu(i,j,k)  =    0.
  110 continue
c
      if(iturb.eq.2) then
        do 111 lpatch= lpatchb,lpatche
c
        j   =   n2dir(lpatch)
        k   =   n3dir(lpatch)
c
        turv1(i,j,k)  =    0.
  111   continue
      endif
c
      if(iturb.eq.3) then
        do 112 lpatch= lpatchb,lpatche
c
        j   =   n2dir(lpatch)
        k   =   n3dir(lpatch)
c
        turv1(i,j,k)  =    0.
        turv2(i,j,k)  =    0.
  112   continue
      endif
c
      go to 1000
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
      eomu(i,j,k)  =    0.
  210 continue
c
      if(iturb.eq.2) then
        do 211 lpatch= lpatchb,lpatche
c
        j   =   n2dir(lpatch)
        k   =   n3dir(lpatch)
c
        turv1(i,j,k)  =    0.
  211   continue
      endif
c
      if(iturb.eq.3) then
        do 212 lpatch= lpatchb,lpatche
c
        j   =   n2dir(lpatch)
        k   =   n3dir(lpatch)
c
        turv1(i,j,k)  =    0.
        turv2(i,j,k)  =    0.
  212   continue
      endif
c
      go to 1000
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
      eomu(i,j,k)  =    0.
  310 continue
c
      if(iturb.eq.2) then
        do 311 lpatch= lpatchb,lpatche
c
        k   =   n3dir(lpatch)
        i   =   n1dir(lpatch)
c
        turv1(i,j,k)  =    0.
  311   continue
      endif
c
      if(iturb.eq.3) then
        do 312 lpatch= lpatchb,lpatche
c
        k   =   n3dir(lpatch)
        i   =   n1dir(lpatch)
c
        turv1(i,j,k)  =    0.
        turv2(i,j,k)  =    0.
  312   continue
      endif
c
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
      eomu(i,j,k)  =    0.
  410 continue
c
      if(iturb.eq.2) then
        do 411 lpatch= lpatchb,lpatche
c
        k   =   n3dir(lpatch)
        i   =   n1dir(lpatch)
c
        turv1(i,j,k)  =    0.
  411   continue
      endif
c
      if(iturb.eq.3) then
        do 412 lpatch= lpatchb,lpatche
c
        k   =   n3dir(lpatch)
        i   =   n1dir(lpatch)
c
        turv1(i,j,k)  =    0.
        turv2(i,j,k)  =    0.
  412   continue
      endif
c
      go to 1000
c
 500  continue
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
      eomu(i,j,k)  =    0.
  510 continue
c
      if(iturb.eq.2) then
        do 511 lpatch= lpatchb,lpatche
c
        i   =   n1dir(lpatch)
        j   =   n2dir(lpatch)
c
        turv1(i,j,k)  =    0.
  511   continue
      endif
c
      if(iturb.eq.3) then
        do 512 lpatch= lpatchb,lpatche
c
        i   =   n1dir(lpatch)
        j   =   n2dir(lpatch)
c
        turv1(i,j,k)  =    0.
        turv2(i,j,k)  =    0.
  512   continue
      endif
c
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
      eomu(i,j,k)  =    0.
  610 continue
c
      if(iturb.eq.2) then
        do 611 lpatch= lpatchb,lpatche
c
        i   =   n1dir(lpatch)
        j   =   n2dir(lpatch)
c
        turv1(i,j,k)  =    0.
  611   continue
      endif
c
      if(iturb.eq.3) then
        do 612 lpatch= lpatchb,lpatche
c
        i   =   n1dir(lpatch)
        j   =   n2dir(lpatch)
c
        turv1(i,j,k)  =    0.
        turv2(i,j,k)  =    0.
  612   continue
      endif
c
 1000    lpatchb = lpatchb + npatchf(nnface)
c
      return
      end
