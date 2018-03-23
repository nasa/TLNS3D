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
      subroutine  bcvwev (imn,jmn,kmn,
     .                    imx ,jmx ,kmx ,imp2 ,jmp2 ,kmp2 , w, p,
     .                                  eomu,smin, turv1 , turv2,
     .                    nface ,n1beg ,n1end ,n2beg ,n2end ,isoln ,
     .                    mx1dwk, amu                               )
c=======================================================================
c 
c     object : to compute the ghost cell values for eddy-viscosity
c              on viscous wall
c
c     initial coding : by v.n.vatsa, nasa langley ( April, 1993)
c
c
c

      dimension      eomu(imp2,jmp2,kmp2), smin(imp2,jmp2,kmp2)
      dimension     turv1(imp2,jmp2,kmp2), turv2(imp2,jmp2,kmp2)
      dimension     w(imp2,jmp2,kmp2,5),   p(imp2,jmp2,kmp2)
c
      dimension     amu(mx1dwk)
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c
       sfac      = sqrt(gamma)*rm/rey
       beta1     = 0.075
       visc1     = suthc/tref
c
c
c      branch-cut for different type of boundary faces
c
c
c     write(6,901) nface,imx,jmx,kmx
c901  format(2x,'nface,i,j,k',4i5)
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
c
      do 110 k= n2beg, n2end
c
      if(iturb.eq.1) then
        do 112 j= n1beg, n1end
        eomu(i,j,k)   =    eomu(i+1,j,k)
  112   continue
      endif
c
      if(iturb.eq.2) then
        do 114 j= n1beg, n1end
        eomu(i,j,k)   =    eomu(i+1,j,k)
        turv1(i,j,k)  =    turv1(i+1,j,k)
  114   continue
      endif
c
      if(iturb.eq.3) then
        if(suthc.le.0.) then
           do 117 j=n1beg, n1end
           ts   = p(i+1,j,k)/w(i+1,j,k,1)
           amu(j) = ts
  117      continue
        else
           do 118 j=n1beg, n1end
           ts     = p(i+1,j,k)/w(i+1,j,k,1)
           amu(j) = ((1.+visc1)/(ts+visc1))*ts**1.5
  118      continue
        endif
c
        do 116 j= n1beg, n1end
c       eomu(i,j,k)   =    eomu(i+1,j,k)
c       turv1(i,j,k)  =   -turv1(i+1,j,k)
        turv1(i,j,k)  =    0.
        turv2(i,j,k)  =    60.*amu(j)*sfac*sfac/(w(i+1,j,k,1)*
     .                       beta1*smin(i+1,j,k)*smin(i+1,j,k) )
        eomu(i,j,k)   =    w(i,j,k,1)*turv1(i,j,k)/(turv2(i,j,k)*amu(j))
  116   continue
      endif
c
  110 continue
c
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
c
      do 210 k= n2beg, n2end
c
      if(iturb.eq.1) then
        do 212 j= n1beg, n1end
        eomu(i,j,k)   =    eomu(i-1,j,k)
  212   continue
      endif
c
      if(iturb.eq.2) then
        do 214 j= n1beg, n1end
        eomu(i,j,k)   =    eomu(i-1,j,k)
        turv1(i,j,k)  =    turv1(i-1,j,k)
  214   continue
      endif
c
      if(iturb.eq.3) then
        if(suthc.le.0.) then
           do 217 j=n1beg, n1end
           ts   = p(i-1,j,k)/w(i-1,j,k,1)
           amu(j) = ts
  217      continue
        else
           do 218 j=n1beg, n1end
           ts     = p(i-1,j,k)/w(i-1,j,k,1)
           amu(j) = ((1.+visc1)/(ts+visc1))*ts**1.5
  218      continue
        endif
c
        do 216 j= n1beg, n1end
c       eomu(i,j,k)   =    eomu(i-1,j,k)
c       turv1(i,j,k)  =   -turv1(i-1,j,k)
        turv1(i,j,k)  =    0.
        turv2(i,j,k)  =    60.*amu(j)*sfac*sfac/(w(i-1,j,k,1)*
     .                       beta1*smin(i-1,j,k)*smin(i-1,j,k) )
        eomu(i,j,k)   =    w(i,j,k,1)*turv1(i,j,k)/(turv2(i,j,k)*amu(j))
  216   continue
      endif
c
  210 continue
c
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
      j   =   jmn
c
      do 310 k= n1beg, n1end
c
      if(iturb.eq.1) then
        do 312 i= n2beg, n2end
        eomu(i,j,k)   =    eomu(i,j+1,k)
  312   continue
      endif
c
      if(iturb.eq.2) then
        do 314 i= n2beg, n2end
        eomu(i,j,k)   =    eomu(i,j+1,k)
        turv1(i,j,k)  =    turv1(i,j+1,k)
  314   continue
      endif
c
      if(iturb.eq.3) then
        if(suthc.le.0.) then
           do 317 i=n2beg, n2end
           ts   = p(i,j+1,k)/w(i,j+1,k,1)
           amu(i) = ts
  317      continue
        else
           do 318 i=n2beg, n2end
           ts     = p(i,j+1,k)/w(i,j+1,k,1)
           amu(i) = ((1.+visc1)/(ts+visc1))*ts**1.5
  318      continue
        endif
c
        do 316 i= n2beg, n2end
c       eomu(i,j,k)   =    eomu(i,j+1,k)
c       turv1(i,j,k)  =   -turv1(i,j+1,k)
        turv1(i,j,k)  =    0.
        turv2(i,j,k)  =    60.*amu(i)*sfac*sfac/(w(i,j+1,k,1)*
     .                       beta1*smin(i,j+1,k)*smin(i,j+1,k) )
        eomu(i,j,k)   =    w(i,j,k,1)*turv1(i,j,k)/(turv2(i,j,k)*amu(j))
  316   continue
      endif
c
  310   continue
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
      j   =   jmx + 1
c
      do 410 k= n1beg, n1end
c
      if(iturb.eq.1) then
        do 412 i= n2beg, n2end
        eomu(i,j,k)   =    eomu(i,j-1,k)
  412   continue
      endif
c
      if(iturb.eq.2) then
        do 414 i= n2beg, n2end
        eomu(i,j,k)   =    eomu(i,j-1,k)
        turv1(i,j,k)  =    turv1(i,j-1,k)
  414   continue
      endif
c
      if(iturb.eq.3) then
        if(suthc.le.0.) then
           do 417 i=n2beg, n2end
           ts   = p(i,j-1,k)/w(i,j-1,k,1)
           amu(i) = ts
  417      continue
        else
           do 418 i=n2beg, n2end
           ts     = p(i,j-1,k)/w(i,j-1,k,1)
           amu(i) = ((1.+visc1)/(ts+visc1))*ts**1.5
  418      continue
        endif
c
        do 416 i= n2beg, n2end
c       eomu(i,j,k)   =    eomu(i,j-1,k)
c       turv1(i,j,k)  =   -turv1(i,j-1,k)
        turv1(i,j,k)  =    0.
        turv2(i,j,k)  =    60.*amu(i)*sfac*sfac/(w(i,j-1,k,1)*
     .                       beta1*smin(i,j-1,k)*smin(i,j-1,k) )
        eomu(i,j,k)   =    w(i,j,k,1)*turv1(i,j,k)/(turv2(i,j,k)*amu(j))
  416   continue
      endif
c
  410 continue
c
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
      k   =   kmn
c
      do 510 j= n2beg, n2end
c
      if(iturb.eq.1) then
        do 512 i= n1beg, n1end
        eomu(i,j,k)   =    eomu(i,j,k+1)
  512   continue
      endif
c
      if(iturb.eq.2) then
        do 514 i= n1beg, n1end
        eomu(i,j,k)   =    eomu(i,j,k+1)
        turv1(i,j,k)  =    turv1(i,j,k+1)
  514   continue
      endif
c
      if(iturb.eq.3) then
        if(suthc.le.0.) then
           do 517 i=n1beg, n1end
           ts   = p(i,j,k+1)/w(i,j,k+1,1)
           amu(i) = ts
  517      continue
        else
           do 518 i=n1beg, n1end
           ts     = p(i,j,k+1)/w(i,j,k+1,1)
           amu(i) = ((1.+visc1)/(ts+visc1))*ts**1.5
  518      continue
        endif
c
        do 516 i= n1beg, n1end
c       eomu(i,j,k)   =    eomu(i,j,k+1)
c       turv1(i,j,k)  =   -turv1(i,j,k+1)
        turv1(i,j,k)  =    0.
        turv2(i,j,k)  =    60.*amu(i)*sfac*sfac/(w(i,j,k+1,1)*
     .                       beta1*smin(i,j,k+1)*smin(i,j,k+1) )
        eomu(i,j,k)   =    w(i,j,k,1)*turv1(i,j,k)/(turv2(i,j,k)*amu(j))
  516   continue
      endif
c
  510 continue
c
      return
c
c
  600 continue
c--------------------- k = kmx  face ----------------------------
c
      if( n1beg .ne. imn) n1beg = n1beg +1
      if( n1end .eq. imx) n1end = n1end +1
      if( n2beg .ne. jmn) n2beg = n2beg +1
      if( n2end .eq. jmx) n2end = n2end +1
c
      k   =   kmx + 1
c
      do 610 j= n2beg, n2end
c
      if(iturb.eq.1) then
        do 612 i= n1beg, n1end
        eomu(i,j,k)   =    eomu(i,j,k-1)
  612   continue
      endif
c
      if(iturb.eq.2) then
        do 614 i= n1beg, n1end
        eomu(i,j,k)   =    eomu(i,j,k-1)
        turv1(i,j,k)  =    turv1(i,j,k-1)
  614   continue
      endif
c
      if(iturb.eq.3) then
        if(suthc.le.0.) then
           do 617 i=n1beg, n1end
           ts   = p(i,j,k-1)/w(i,j,k-1,1)
           amu(i) = ts
  617      continue
        else
           do 618 i=n1beg, n1end
           ts     = p(i,j,k-1)/w(i,j,k-1,1)
           amu(i) = ((1.+visc1)/(ts+visc1))*ts**1.5
  618      continue
        endif
c
        do 616 i= n1beg, n1end
c       eomu(i,j,k)   =    eomu(i,j,k-1)
c       turv1(i,j,k)  =   -turv1(i,j,k-1)
        turv1(i,j,k)  =    0.
        turv2(i,j,k)  =    60.*amu(i)*sfac*sfac/(w(i,j,k-1,1)*
     .                       beta1*smin(i,j,k-1)*smin(i,j,k-1) )
        eomu(i,j,k)   =    w(i,j,k,1)*turv1(i,j,k)/(turv2(i,j,k)*amu(j))
  616   continue
      endif
c
  610   continue
c
c
      return
      end
