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
      subroutine  bcfarev(imn,jmn,kmn,
     .                    imx ,jmx ,kmx ,imp2 ,jmp2 ,kmp2 , eomu ,
     .                                              turv1 , turv2,
     .                    nface ,n1beg ,n1end ,n2beg ,n2end ,isoln )
c=======================================================================
c 
c     object : to compute the ghost cell values for eddy-viscosity
c              by specifying far-field values
c
c     initial coding : by v.n.vatsa, nasa langley ( April, 1993)
c
c
c
c

      common/trbfs/ anutinf,zkinf,ominf
c
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension      eomu(imp2,jmp2,kmp2)
      dimension     turv1(imp2,jmp2,kmp2),turv2(imp2,jmp2,kmp2)
c
      evinf     =  0.0
c
      if(iturb.eq.2) then
         chi    =  anutinf
         fv1    =  chi**3/(chi**3 + 7.1**3)
         evinf  =  fv1*anutinf
      endif
c
      if(iturb.eq.3) then
         evinf  =  zkinf/ominf
      endif
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
c
      do 110 k= n2beg, n2end
c
      if(iturb.eq.1) then
        do 112 j= n1beg, n1end
        eomu(i,j,k)   =  evinf
  112   continue
      endif
c
      if(iturb.eq.2) then
        do 113 j= n1beg, n1end
        eomu(i,j,k)   =  evinf
        turv1(i,j,k)  =  anutinf
  113   continue
      endif
c
      if(iturb.eq.3) then
        do 114 j= n1beg, n1end
        eomu(i,j,k)   =  evinf
        turv1(i,j,k)  =  zkinf
        turv2(i,j,k)  =  ominf
  114   continue
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
        eomu(i,j,k)   =  evinf
  212   continue
      endif
c
      if(iturb.eq.2) then
        do 213 j= n1beg, n1end
        eomu(i,j,k)   =  evinf
        turv1(i,j,k)  =  anutinf
  213   continue
      endif
c
      if(iturb.eq.3) then
        do 214 j= n1beg, n1end
        eomu(i,j,k)   =  evinf
        turv1(i,j,k)  =  zkinf
        turv2(i,j,k)  =  ominf
  214   continue
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
        eomu(i,j,k)   =  evinf
  312   continue
      endif
c
      if(iturb.eq.2) then
        do 313 i= n2beg, n2end
        eomu(i,j,k)   =  evinf
        turv1(i,j,k)  =  anutinf
  313   continue
      endif
c
      if(iturb.eq.3) then
        do 314 i= n2beg, n2end
        eomu(i,j,k)   =  evinf
        turv1(i,j,k)  =  zkinf
        turv2(i,j,k)  =  ominf
  314   continue
      endif
c
  310 continue
c
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
        eomu(i,j,k)   =  evinf
  412   continue
      endif
c
      if(iturb.eq.2) then
        do 413 i= n2beg, n2end
        eomu(i,j,k)   =  evinf
        turv1(i,j,k)  =  anutinf
  413   continue
      endif
c
      if(iturb.eq.3) then
        do 414 i= n2beg, n2end
        eomu(i,j,k)   =  evinf
        turv1(i,j,k)  =  zkinf
        turv2(i,j,k)  =  ominf
  414   continue
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
        eomu(i,j,k)   =  evinf
  512   continue
      endif
c
      if(iturb.eq.2) then
        do 513 i= n1beg, n1end
        eomu(i,j,k)   =  evinf
        turv1(i,j,k)  =  anutinf
  513   continue
      endif
c
      if(iturb.eq.3) then
        do 514 i= n1beg, n1end
        eomu(i,j,k)   =  evinf
        turv1(i,j,k)  =  zkinf
        turv2(i,j,k)  =  ominf
  514   continue
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
        eomu(i,j,k)   =  evinf
  612   continue
      endif
c
      if(iturb.eq.2) then
        do 613 i= n1beg, n1end
        eomu(i,j,k)   =  evinf
        turv1(i,j,k)  =  anutinf
  613   continue
      endif
c
      if(iturb.eq.3) then
        do 614 i= n1beg, n1end
        eomu(i,j,k)   =  evinf
        turv1(i,j,k)  =  zkinf
        turv2(i,j,k)  =  ominf
  614   continue
      endif
c
  610 continue
c
c
      return
      end
