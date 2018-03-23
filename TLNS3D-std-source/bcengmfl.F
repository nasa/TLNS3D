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
      subroutine  bcengmfl(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                    imp2,jmp2,kmp2,six,siy,siz,sjx,sjy,sjz,
     .                    skx,sky,skz,w,p,x,eomu,turv1,nface,n1beg,
     .                    n1end,n2beg,n2end,igrid,isoln,extyp,
     .                    epseng,relax,nodes                      )
c=======================================================================
c
c     modified by T. Reyhner using N.J.Yu's code as a model to set
c     a specified mass flow at duct exit - March 1997
c     switched from varing momentum to modifying pressure - seems
c     to work well (tested  by T. Reyhner)
c
c     Note: changes only along imin and imax faces
c 
c     object : to compute the ghost cell values for flow variables
c              at downstream planes by extrapolation
c
c     latest mods : by v.n.vatsa, nasa langley ( Jan. 1998)
c    ( To make it compatible with current map file. Also added
c      error exit code if this b.c. is applied to j and k faces)
c
c     w(i,j,k,1)  = density 
c     w(i,j,k,2)  = momentum in x direction 
c     w(i,j,k,3)  = momentum in y direction 
c     w(i,j,k,4)  = momentum in z direction 
c     w(i,j,k,5)  = total energy
c     p(i,j,k)    = pressure
c     eomu(i,j,k) = eddy-viscosity
c
c
c

      dimension      w(imp2,jmp2,kmp2,5),
     .               p(imp2,jmp2,kmp2)  ,  eomu(imp2,jmp2,kmp2),
     .                                    turv1(imp2,jmp2,kmp2)
      dimension      x(imp1,jmp1,kmp1,3)
      dimension      six(imp1,jmp1,kmp1), siy(imp1,jmp1,kmp1),
     .               siz(imp1,jmp1,kmp1), sjx(imp1,jmp1,kmp1),
     .               sjy(imp1,jmp1,kmp1), sjz(imp1,jmp1,kmp1),
     .               skx(imp1,jmp1,kmp1), sky(imp1,jmp1,kmp1),
     .               skz(imp1,jmp1,kmp1)
c
      dimension nodes(1)
c
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c
c------------------------------------------------------------------
c
c     extyp  =   0  zero order extrapolation 
c                1  linear     extrapolation  (not yet coded)
c
c------------------------------------------------------------------
c
c     extyp = 0.0
c
c      check the permissibility of the 'extyp' 
c
       coef=sqrt(gamma)*rm
c      the code has worked with relax=0.1, but there are overshoots
c      and a oscillatory pattern of convergence.
c      relax=0.025
c
cBWW moved this test to imapin.f
cBWW  if(extyp .gt.1.01 .or. extyp .lt. -0.01 )   then
cBWW     write (iwrit,'(2x," extyp = ",f8.4, " not permissible")')
cBWW     stop
cBWW  endif
c
c      branch-cut for different type of boundary faces
c
      go to (100,200,300,400,500,600) nface
c
c
  100 continue
c--------------------- i = imn  face ----------------------------
c
      if( n1beg .ne. jmn) n1beg = n1beg +1
      if( n1end .eq. jmx) n1end = n1end +1
      if( n2beg .ne. kmn) n2beg = n2beg +1
      if( n2end .eq. kmx) n2end = n2end +1
c
      i   =   imn

      if( igrid .eq. isoln) then

        dxtest=x(i,n1beg,n2beg,1)-x(i+1,n1beg,n2beg,1)
        signx=1.0
        if(dxtest .lt. 0.0) signx=-1.0
c
c       compute mass flux at fan face 
c
        rqa=0.0
        afanc=0.0
c
        its=imn
        do 10 k=n2beg+1,n2end-1
        do 10 j=n1beg+1,n1end-1
        ru=w(its,j,k,2)+w(i+1,j,k,2)
        rv=w(its,j,k,3)+w(i+1,j,k,3)
        rw=w(its,j,k,4)+w(i+1,j,k,4)
        rqn=-(ru*six(its,j,k)+rv*siy(its,j,k)+rw*siz(its,j,k))/coef
        rqa=rqa+rqn*signx
        afanc=afanc-six(i,j,k)*signx
  10    continue
        fact=2.0*epseng*afanc/rqa
        fact2=min(1.02,fact)
        fact2=max(0.98,fact2)
        fact2=1.0+relax*(1.0-fact2)
        rqb=rqa
c
c       debug output
c       kkk=(n2beg+n2end)/2
c       jjj=(n1beg+n1end)/2
c       write (81,14)
c 14    format(1x,'    mass flow      ainf/afan           afan  ',
c    1            '         fact          fact2      i    j    k',
c    2            '     p(i,j,k)  ')
c       write (81,15) rqb,epseng,afanc,fact,fact2,i,jjj,kkk,p(i,jjj,kkk)
c 15    format(1x,e15.7,f15.7,e15.7,2f15.7,3i5,e15.7)
c
      else
c
        fact2=1.0
c
      endif
c
      do 110 k= n2beg, n2end
      do 110 j= n1beg, n1end
c
      w(i,j,k,1)  = w(i+1,j,k,1)
      w(i,j,k,2)  = w(i+1,j,k,2)
      w(i,j,k,3)  = w(i+1,j,k,3)
      w(i,j,k,4)  = w(i+1,j,k,4)
      w(i,j,k,5)  = w(i+1,j,k,5)
      p(i,j,k)    = fact2*p(i,j,k)
c
      eomu(i,j,k) = eomu(i+1,j,k)
ccc   turv1(i,j,k)=turv1(i+1,j,k)
  110 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         i   =  imn - 1
c
         do 120 k= n2beg, n2end
         do 120 j= n1beg, n1end
         w(i,j,k,1)   =    w(i+1,j,k,1)
         w(i,j,k,2)   =    w(i+1,j,k,2)
         w(i,j,k,3)   =    w(i+1,j,k,3)
         w(i,j,k,4)   =    w(i+1,j,k,4)
         w(i,j,k,5)   =    w(i+1,j,k,5)
         p(i,j,k)     =    p(i+1,j,k)
         eomu(i,j,k)  =    eomu(i+1,j,k)
ccc      turv1(i,j,k) =   turv1(i+1,j,k)
  120    continue
      endif
      return
c
c
  200 continue
c--------------------- i = imx  face ----------------------------
c
      if( n1beg .ne. jmn) n1beg = n1beg +1
      if( n1end .eq. jmx) n1end = n1end +1
      if( n2beg .ne. kmn) n2beg = n2beg +1
      if( n2end .eq. kmx) n2end = n2end +1
c
      i   =   imx + 1

      if( igrid .eq. isoln) then

        dxtest=x(i-2,n1beg,n2beg,1)-x(i-1,n1beg,n2beg,1)
        signx=1.0
        if(dxtest .lt. 0.0) signx=-1.0
c
c       compute mass flux at fan face 
c
        its=imx
        rqa=0.0
        afanc=0.0
        do 201 k=n2beg+1,n2end-1
        do 201 j=n1beg+1,n1end-1
        ru=w(its,j,k,2)+w(i-1,j,k,2)
        rv=w(its,j,k,3)+w(i-1,j,k,3)
        rw=w(its,j,k,4)+w(i-1,j,k,4)
        rqn=ru*six(its,j,k)+rv*siy(its,j,k)+rw*siz(its,j,k)
        rqa=rqa-signx*rqn/coef
        afanc=afanc-signx*six(its,j,k)
  201   continue
        fact=2.0*epseng*afanc/rqa
        fact2=min(1.02,fact)
        fact2=max(0.98,fact2)
        fact2=1.0+relax*(1.0-fact2)
        rqb=rqa
c       kkk=(n2beg+n2end)/2
c       jjj=(n1beg+n1end)/2
c       write (81,14)
c       write (81,15) rqb,epseng,afanc,fact,fact2,i,jjj,kkk,p(i,jjj,kkk)
c
      else
c
        fact2=1.0
c
      endif
c
      do 206 k=n2beg,n2end
      do 202 j=n1beg,n1end
      w(i,j,k,1)  = w(i-1,j,k,1)
      w(i,j,k,2)  = w(i-1,j,k,2)
      w(i,j,k,3)  = w(i-1,j,k,3)
      w(i,j,k,4)  = w(i-1,j,k,4)
      w(i,j,k,5)  = w(i-1,j,k,5)
      p(i,j,k)    = fact2*p(i,j,k)
      eomu(i,j,k) = eomu(i-1,j,k)
ccc   turv1(i,j,k)=  turv1(i-1,j,k)
  202 continue
  206 continue
c     k=n2beg+1
c     write(82,207) k,i,imx,fact
c 207 format(' output from bcfan1c, k,i,imx,fac=',3i5,e12.5,/,
c    1  '  j,r,u,v,w,e,u,v,w at i=imx,imx+1')
c     do 208 j=n1beg,n1end
c     write(82,209) j,w(i-1,j,k,1),w(i-1,j,k,2),w(i-1,j,k,3),
c    1  w(i-1,j,k,4),w(i-1,j,k,5),w(i,j,k,2),w(i,j,k,3),w(i,j,k,4)
c 208 continue
c 209 format(2x,i5,5f10.4,2x,3f10.4)
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         i   =  imx + 2
c
         do 220 k= n2beg, n2end
         do 220 j= n1beg, n1end
         w(i,j,k,1)   =    w(i-1,j,k,1)
         w(i,j,k,2)   =    w(i-1,j,k,2)
         w(i,j,k,3)   =    w(i-1,j,k,3)
         w(i,j,k,4)   =    w(i-1,j,k,4)
         w(i,j,k,5)   =    w(i-1,j,k,5)
         p(i,j,k)     =    p(i-1,j,k)
         eomu(i,j,k)  =    eomu(i-1,j,k)
ccc      turv1(i,j,k) =   turv1(i-1,j,k)
  220    continue
      endif
      return
c
c
  300 continue
c--------------------- j = jmn  face ----------------------------
c
      if( n1beg .ne. kmn) n1beg = n1beg +1
      if( n1end .eq. kmx) n1end = n1end +1
      if( n2beg .ne. imn) n2beg = n2beg +1
      if( n2end .eq. imx) n2end = n2end +1
c
c***********************       warning    ************************
        write (iwrit,'(" engine b.c. not coded for jmn face")')
        call ERREXIT(nodes)
c*****************************************************************
c
      j   =   jmn
c
      do 310 k= n1beg, n1end
      do 310 i= n2beg, n2end
      w(i,j,k,1)   =    w(i,j+1,k,1)
      w(i,j,k,2)   =    w(i,j+1,k,2)
      w(i,j,k,3)   =    w(i,j+1,k,3)
      w(i,j,k,4)   =    w(i,j+1,k,4)
      w(i,j,k,5)   =    w(i,j+1,k,5)
c     p(i,j,k)     =    p(i,j+1,k)
      p(i,j,k)     =      pfu
      eomu(i,j,k)  =    eomu(i,j+1,k)
ccc   turv1(i,j,k) =   turv1(i,j+1,k)
  310 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         j   =  jmn - 1
c
         do 320 k= n1beg, n1end
         do 320 i= n2beg, n2end
         w(i,j,k,1)   =    w(i,j+1,k,1)
         w(i,j,k,2)   =    w(i,j+1,k,2)
         w(i,j,k,3)   =    w(i,j+1,k,3)
         w(i,j,k,4)   =    w(i,j+1,k,4)
         w(i,j,k,5)   =    w(i,j+1,k,5)
         p(i,j,k)     =    p(i,j+1,k)
         eomu(i,j,k)  =    eomu(i,j+1,k)
ccc      turv1(i,j,k) =   turv1(i,j+1,k)
  320    continue
      endif
      return
c
c
  400 continue
c--------------------- j = jmx  face ----------------------------
c
      if( n1beg .ne. kmn) n1beg = n1beg +1
      if( n1end .eq. kmx) n1end = n1end +1
      if( n2beg .ne. imn) n2beg = n2beg +1
      if( n2end .eq. imx) n2end = n2end +1
c
c***********************       warning    ************************
        write (iwrit,'(" engine b.c. not coded for jmx face")')
        call ERREXIT(nodes)
c*****************************************************************
c
      j   =   jmx + 1
c
      do 410 k= n1beg, n1end
      do 410 i= n2beg, n2end
      w(i,j,k,1)   =    w(i,j-1,k,1)
      w(i,j,k,2)   =    w(i,j-1,k,2)
      w(i,j,k,3)   =    w(i,j-1,k,3)
      w(i,j,k,4)   =    w(i,j-1,k,4)
      w(i,j,k,5)   =    w(i,j-1,k,5)
c     p(i,j,k)     =    p(i,j-1,k)
      p(i,j,k)     =      pfu
      eomu(i,j,k)  =    eomu(i,j-1,k)
ccc   turv1(i,j,k) =   turv1(i,j-1,k)
  410 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         j   =  jmx + 2
c
         do 420 k= n1beg, n1end
         do 420 i= n2beg, n2end
         w(i,j,k,1)   =    w(i,j-1,k,1)
         w(i,j,k,2)   =    w(i,j-1,k,2)
         w(i,j,k,3)   =    w(i,j-1,k,3)
         w(i,j,k,4)   =    w(i,j-1,k,4)
         w(i,j,k,5)   =    w(i,j-1,k,5)
         p(i,j,k)     =    p(i,j-1,k)
         eomu(i,j,k)  =    eomu(i,j-1,k)
ccc      turv1(i,j,k) =   turv1(i,j-1,k)
  420    continue
      endif
      return
c
c
  500 continue
c--------------------- k = kmn  face ----------------------------
c
      if( n1beg .ne. imn) n1beg = n1beg +1
      if( n1end .eq. imx) n1end = n1end +1
      if( n2beg .ne. jmn) n2beg = n2beg +1
      if( n2end .eq. jmx) n2end = n2end +1
c
c***********************       warning    ************************
        write (iwrit,'(" engine b.c. not coded for kmn face")')
        call ERREXIT(nodes)
c*****************************************************************
c
      k   =   kmn
c
      do 510 j= n2beg, n2end
      do 510 i= n1beg, n1end
      w(i,j,k,1)   =    w(i,j,k+1,1)
      w(i,j,k,2)   =    w(i,j,k+1,2)
      w(i,j,k,3)   =    w(i,j,k+1,3)
      w(i,j,k,4)   =    w(i,j,k+1,4)
      w(i,j,k,5)   =    w(i,j,k+1,5)
c     p(i,j,k)     =    p(i,j,k+1)
      p(i,j,k)     =      pfu
      eomu(i,j,k)  =    eomu(i,j,k+1)
ccc   turv1(i,j,k) =   turv1(i,j,k+1)
  510 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         k   =  kmn - 1
c
         do 520 j= n2beg, n2end
         do 520 i= n1beg, n1end
         w(i,j,k,1)   =    w(i,j,k+1,1)
         w(i,j,k,2)   =    w(i,j,k+1,2)
         w(i,j,k,3)   =    w(i,j,k+1,3)
         w(i,j,k,4)   =    w(i,j,k+1,4)
         w(i,j,k,5)   =    w(i,j,k+1,5)
         p(i,j,k)     =    p(i,j,k+1)
         eomu(i,j,k)  =    eomu(i,j,k+1)
ccc      turv1(i,j,k) =   turv1(i,j,k+1)
  520    continue
      endif
      return
c
c
  600 continue
c--------------------- k = kmx  face ----------------------------
c
      if( n1beg .ne. jmn) n1beg = n1beg +1
      if( n1end .eq. jmx) n1end = n1end +1
      if( n2beg .ne. kmn) n2beg = n2beg +1
      if( n2end .eq. kmx) n2end = n2end +1
c
c***********************       warning    ************************
        write (iwrit,'(" engine b.c. not coded for kmx face")')
        call ERREXIT(nodes)
c*****************************************************************
c
      k   =   kmx + 1
c
      do 610 j= n2beg, n2end
      do 610 i= n1beg, n1end
      w(i,j,k,1)   =    w(i,j,k-1,1)
      w(i,j,k,2)   =    w(i,j,k-1,2)
      w(i,j,k,3)   =    w(i,j,k-1,3)
      w(i,j,k,4)   =    w(i,j,k-1,4)
      w(i,j,k,5)   =    w(i,j,k-1,5)
c     p(i,j,k)     =    p(i,j,k-1)
      p(i,j,k)     =      pfu
      eomu(i,j,k)  =    eomu(i,j,k-1)
ccc   turv1(i,j,k) =   turv1(i,j,k-1)
  610 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
c
         k   =  kmx + 2
c
         do 620 j= n2beg, n2end
         do 620 i= n1beg, n1end
         w(i,j,k,1)   =    w(i,j,k-1,1)
         w(i,j,k,2)   =    w(i,j,k-1,2)
         w(i,j,k,3)   =    w(i,j,k-1,3)
         w(i,j,k,4)   =    w(i,j,k-1,4)
         w(i,j,k,5)   =    w(i,j,k-1,5)
         p(i,j,k)     =    p(i,j,k-1)
         eomu(i,j,k)  =    eomu(i,j,k-1)
ccc      turv1(i,j,k) =   turv1(i,j,k-1)
  620    continue
      endif
c
c
      return
      end
