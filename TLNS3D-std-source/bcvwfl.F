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
      subroutine  bcvwfl(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                   imp2,jmp2,kmp2,w,p,eomu,turv1,turv2,
     .                   six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                   ri1,ri2,rj1,rj2,rk1,rk2,mx1dwk,dp,
     .                   nface,n1beg,n1end,n2beg,n2end,igrid,
     .                   isoln,extyp,twotref,cmassmx,cporous,pchambr,
     .                   tperiod,totptim,iunsteady,nodes,ibloc,myrank)
c=======================================================================
c 
c     object : to compute the ghost cell values for flow variables
c              at viscous walls
c
c     initial coding : by v.n.vatsa, nasa langley ( March, 1990)
c     Latest Mods: Oscillatory mass flow at imin face (Sept. 2001) 
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
c
      data init/0/

      dimension      w(imp2,jmp2,kmp2,5),
     .               p(imp2,jmp2,kmp2)   ,  eomu(imp2,jmp2,kmp2),
     .               turv1(imp2,jmp2,kmp2),turv2(imp2,jmp2,kmp2)
c
      dimension      six(imp1,jmp1,kmp1), siy(imp1,jmp1,kmp1),
     .               siz(imp1,jmp1,kmp1), sjx(imp1,jmp1,kmp1),
     .               sjy(imp1,jmp1,kmp1), sjz(imp1,jmp1,kmp1),
     .               skx(imp1,jmp1,kmp1), sky(imp1,jmp1,kmp1),
     .               skz(imp1,jmp1,kmp1)
c
      dimension      ri1(jmp1,kmp1,3)   , ri2(jmp1,kmp1,3),
     .               rj1(imp1,kmp1,3)   , rj2(imp1,kmp1,3),
     .               rk1(imp1,jmp1,3)   , rk2(imp1,jmp1,3)
c
      dimension      dp(mx1dwk,mx1dwk)
c
c     temporary storage for variable at boundary-faces
c
      dimension      wbound(mx1dwk,mx1dwk,5),pbound(mx1dwk,mx1dwk)
c
      dimension      nodes(1)
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      common/trbfs/ anutinf,zkinf,ominf
c
c
c------------------------------------------------------------------
c     extyp  =   0  zero order extrapolation for p at wall boundary
c                1  linear     extrapolation for p at wall boundary
c               -1  normal momentum for pressure at wall boundary
c
c------------------------------------------------------------------
c
c
c     extyp = 0.0
c
c      check the permissibility of the 'extyp' 
c
cBWW moved this test to imapin.f
cBWW  if(extyp .gt.1.01 .or. extyp .lt. -1.01 )   then
cBWW     write (iwrit,'(2x," extyp = ",f8.4, " not permissible")')
cBWW     stop
cBWW  endif
c
c      define a switch for adiabatic or specified wall temperature
c      condition
c
c----- mods by vatsa (10-3-94) to compute only interior cells
c
cvnv   n1beg = n1beg + 1
cvnv   n2beg = n2beg + 1
c-----
c
c      define free-stream turbulence quantities
c
c
      sqrtgam   = sqrt(gamma)
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
         pi     = 3.1415927
         pic    = 3.1415927
         cmassw = cmassmx
         iadb = 0
         if (twotref.le.0.0) iadb = 1
c
c      branch-cut for different type of boundary faces
c
      go to (100,200,300,400,500,600) nface
c
c
  100 continue
c--------------------- i = imn  face ----------------------------
c
c----- mods by vatsa (10-3-94) to compute only interior cells
      if( n1beg .ne. jmn) n1beg = n1beg +1
      if( n1end .eq. jmx) n1end = n1end +1
      if( n2beg .ne. kmn) n2beg = n2beg +1
      if( n2end .eq. kmx) n2end = n2end +1
c
cvn   if(cmassw.gt.0.) then
cvn     write (iwrit,'("suction/injection  not coded for imn face",
cvn  .                 2x,"cmassw = ",e15.6)')
cvn     call ERREXIT (nodes)
cvn   endif
c
      if(abs(cporous).gt.0.) then
        write (iwrit,'("porosity  not coded for imn face",
     .                 2x,"cporous = ",e15.6)')
        call ERREXIT (nodes)
      endif
c
      i    =    imn
c
c     save old values at boundary
c
      do k=n2beg, n2end
        do j=n1beg, n1end
          wbound(j,k,1) = w(i,j,k,1) 
          wbound(j,k,2) = w(i,j,k,2) 
          wbound(j,k,3) = w(i,j,k,3) 
          wbound(j,k,4) = w(i,j,k,4) 
          wbound(j,k,5) = w(i,j,k,5) 
          pbound(j,k)   = p(i,j,k)
        enddo
      enddo
c
      if (extyp .lt. 0.) go to 150
      if (extyp .gt. 0.) go to 140
c
c--------     zero order extrapolation for surface pressure
      do 130 k= n2beg, n2end
      do 130 j= n1beg, n1end
      dp(j,k)      = 0.
 130  continue
      go to 169
c
c--------     linear extrapolation for surface pressure
 140  continue
      do 145 k= n2beg, n2end
      do 145 j= n1beg, n1end
      dp(j,k)      = p(i+2,j,k) - p(i+1,j,k)
 145  continue
      go to 169
 150  continue
c
c--------- extrapolation using normal momentum equation at surface
c
      do 155 k= n2beg, n2end
      kp        = k + 1
      km        = k - 1
      do 155 j= n1beg, n1end
      jp        = j + 1
      jm        = j - 1
      skxa      = (skx(i+1,j,k) + skx(i+1,j,k-1))
      skya      = (sky(i+1,j,k) + sky(i+1,j,k-1))
      skza      = (skz(i+1,j,k) + skz(i+1,j,k-1))
      sjxa      = (sjx(i+1,j,k) + sjx(i+1,j-1,k))
      sjya      = (sjy(i+1,j,k) + sjy(i+1,j-1,k))
      sjza      = (sjz(i+1,j,k) + sjz(i+1,j-1,k))
      sixa      = 2.*six(i,j,k)
      siya      = 2.*siy(i,j,k)
      siza      = 2.*siz(i,j,k)
      rxk       = .5*(ri1(j,kp,1) - ri1(j,km,1))
      ryk       = .5*(ri1(j,kp,2) - ri1(j,km,2))
      rzk       = .5*(ri1(j,kp,3) - ri1(j,km,3))
      pk        = .5*(p(i+1,j,kp)   - p(i+1,j,km))
      rxj       = .5*(ri1(jp,k,1) - ri1(jm,k,1))
      ryj       = .5*(ri1(jp,k,2) - ri1(jm,k,2))
      rzj       = .5*(ri1(jp,k,3) - ri1(jm,k,3))
      pj        = .5*(p(i+1,jp,k)   - p(i+1,jm,k))
      si        = ri1(j,k,1)*sixa + ri1(j,k,2)*siya + ri1(j,k,3)*siza 
      sj        = ri1(j,k,1)*sjxa + ri1(j,k,2)*sjya + ri1(j,k,3)*sjza 
      sk        = ri1(j,k,1)*skxa + ri1(j,k,2)*skya + ri1(j,k,3)*skza 
      qk        = w(i+1,j,k,2)*skxa + w(i+1,j,k,3)*skya 
     .          + w(i+1,j,k,4)*skza
      qj        = w(i+1,j,k,2)*sjxa + w(i+1,j,k,3)*sjya
     .          + w(i+1,j,k,4)*sjza
      dp(j,k)   = (qk*(w(i+1,j,k,2)*rxk + w(i+1,j,k,3)*ryk
     .                +w(i+1,j,k,4)*rzk)/ w(i+1,j,k,1) 
     .            +qj*(w(i+1,j,k,2)*rxj + w(i+1,j,k,3)*ryj
     .                +w(i+1,j,k,4)*rzj)/ w(i+1,j,k,1) 
     .            -sk*pk - sj*pj)/si
  155 continue
      if (extyp .eq. -1.) go to 169
c
c--------- improved extrapolation using quadratic fit 
c          and normal pressure gradient at the surface
c
      do 160 k= n2beg,n2end
      do 160 j= n1beg,n1end
      dp(j,k)   = .75*dp(j,k) + .25*(p(i+2,j,k) - p(i+1,j,k))
  160 continue
c
c
  169  continue
c
c     modified for normal velocity sepcification (July, 2001)
c
      do j= n1beg, n1end
        do k= n2beg, n2end
          p(i,j,k)     =    dim(p(i+1,j,k),dp(i,k))
          p(i,j,k)     =    max(p(i,j,k),0.001)
c
          w(i,j,k,1)   =    w(i+1,j,k,1)
        enddo
      enddo
c
c     sinusoidal variation (used in CFDVAL2004 workshop)
cvn   if(iunsteady.gt.0.and.abs(tperiod).gt.0.0000001)
cvn  .                   cmassw = cmassmx*sin(2.*pi*totptim/tperiod)
c     write(81,16) iunsteady,nodes(ibloc),myrank,
c    .             cmassmx,pi,totptim,tperiod,cmassw
c16   format(2x,'iuns,cmax,pi,tot,tp,cmw',3i3,5e16.5)
c
c     curve fit from George Huang (UKY)
c     PIV
c     if(iunsteady.gt.0.and.abs(tperiod).gt.0.0000001) then
c       a1     = 0.3747050
c       a2     = 1.7783791e-02
c       a3     = 2.7846638e-02
c       a4     = 5.5155195e-03
c       a5     = 1.7597616e-02
c       a6     = 9.7805797e-04
c
c       theta3 = 52.
c       theta4 = 144.4
c       theta5 = 108.8
c       theta6 = 75.2
c
c       tratio = totptim/tperiod
c       sum1   = -    a1*sin(2.*pi*tratio) + a2*cos(2.*pi*tratio)
c    .           +    a3*sin(2.*pi*(tratio - theta3/360.))
c    .           - 2.*a4*sin(4.*pi*(tratio - theta4/360.))
c    .           - 3.*a5*sin(6.*pi*(tratio - theta5/360.))
c    .           - 4.*a6*sin(8.*pi*(tratio - theta6/360.))
c
c     FFT curve fit from PIV velocity data (with phase shift) at jet exit
c     (July 2004)
c
      if(iunsteady.gt.0.and.abs(tperiod).gt.0.0000001) then
c
        tratio = totptim/tperiod

        a1     = 0.7715591E+01
        a2     = 0.3625974E+00
        a3     = -0.2882919E+01
        a4     = -0.5151815E+00
        a5     = -0.8232565E+00
        b1     = 0.2149226E+02
        b2     = 0.5741831E+01
        b3     = 0.1585448E+01
        b4     = 0.6679176E+00
        b5     = 0.6231386E+00

        sum1   =      a1*cos(2.*pic*tratio) + b1*sin(2.*pic*tratio)
     .         +      a2*cos(4.*pic*tratio) + b2*sin(4.*pic*tratio)
     .         +      a3*cos(6.*pic*tratio) + b3*sin(6.*pic*tratio)
     .         +      a4*cos(8.*pic*tratio) + b4*sin(8.*pic*tratio)
c       convert to m/sec
        cmassw = 2.*pic*sum1*cmassmx*447.5/1000.
      endif
c

      if(abs(cmassw).ge.0.0000001) go to 171
c
      do 170 j= n1beg, n1end
      do 170 k= n2beg, n2end
ccc   w(i,j,k,1)   =    w(i+1,j,k,1)
ccc   w(i,j,k,2)   =   -w(i+1,j,k,2)
ccc   w(i,j,k,3)   =   -w(i+1,j,k,3)
ccc   w(i,j,k,4)   =   -w(i+1,j,k,4)
      eomu(i,j,k)  =    eomu(i+1,j,k)
ccc   turv1(i,j,k) =   turv1(i+1,j,k)
      p(i,j,k)     =    dim(p(i+1,j,k),dp(j,k))
      rhoi         =    1./w(i+1,j,k,1)
      ugc          =   -w(i+1,j,k,2)*rhoi
      vgc          =   -w(i+1,j,k,3)*rhoi
      wgc          =   -w(i+1,j,k,4)*rhoi
      temp         =    p(i+1,j,k)*rhoi
      tsgc         =    max((2.*twotref - temp),0.8*twotref)*(1-iadb)
     .             +    iadb*temp
      p(i,j,k)     =    max(p(i,j,k),0.001)
      w(i,j,k,1)   =    p(i,j,k)/tsgc
      w(i,j,k,2)   =    w(i,j,k,1)*ugc
      w(i,j,k,3)   =    w(i,j,k,1)*vgc
      w(i,j,k,4)   =    w(i,j,k,1)*wgc
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  170 continue
c
      go to 181
c
c     mass transfer (blowing if cmassw>0, suction if cmassw<0) normal to wall
c
  171 continue
      vnorm        = sqrt(gamma)*rm*cmassw
c     specify momentum
c     vnorm        = sqrt(gamma)*rm*vnorm
      do 172 j= n1beg, n1end
      do 172 k= n2beg, n2end
cvn   if(cmassw.lt.0.) p(i,j,k) = pchambr
c(fix temp. of cavity)
      w(i,j,k,1)   = p(i,j,k)/pchambr
      w(i,j,k,2)   = 2.*vnorm*ri1(j,k,1) - w(i+1,j,k,2)
      w(i,j,k,3)   = 2.*vnorm*ri1(j,k,2) - w(i+1,j,k,3)
      w(i,j,k,4)   = 2.*vnorm*ri1(j,k,3) - w(i+1,j,k,4)
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
c     if(cmassw.gt.0) go to 172
c     w(i,j,k,1)  = w(i+1,j,k,1)
c     w(i,j,k,2)  = w(i+1,j,k,2)
c     w(i,j,k,3)  = w(i+1,j,k,3)
c     w(i,j,k,4)  = w(i+1,j,k,4)
c     w(i,j,k,5)  = w(i+1,j,k,5)
c     p(i,j,k)    = p(i+1,j,k)
  172 continue
c
      go to 181

cif   if( igrid .eq. isoln) then

        signx=1.0
c
c       compute momentum flux at fan face
c
        rqa=0.0
        rqb=0.0
        afanc=0.0
        coef=sqrt(gamma)*rm
c
        its=imn
        do 10 k=n2beg+1,n2end-1
        do 10 j=n1beg+1,n1end-1
c       ru=w(its,j,k,2)+w(i+1,j,k,2)
c       rv=w(its,j,k,3)+w(i+1,j,k,3)
c       rw=w(its,j,k,4)+w(i+1,j,k,4)
        ru=wbound(j,k,2)+w(i+1,j,k,2)
        rv=wbound(j,k,3)+w(i+1,j,k,3)
        rw=wbound(j,k,4)+w(i+1,j,k,4)
cvn     rqn=-(ru*six(its,j,k)+rv*siy(its,j,k)+rw*siz(its,j,k))/coef
        rqn=-(ri1(j,k,1)*six(its,j,k)+ri1(j,k,2)*siy(its,j,k)
     .                               +ri1(j,k,3)*siz(its,j,k))/coef
cvn     rqb=rqb-rqn*signx*.5
        rqb=rqb-rqn*signx
        rhol = 0.5*(wbound(j,k,1)+w(i+1,j,k,1))
        vtots = 0.25*(ru**2+rv**2+rw**2)/(coef*coef*rhol)
        rqa= rqa + vtots*six(i,j,k)/rhol
        afanc=afanc-six(i,j,k)*signx
  10    continue
c       write (81,15) n2beg,n2end,rqb,rqa,afanc,fact
  15    format(2x,'k,j,rqb,rqa,afanc,six=',2i5,4e18.6)
c
c       adjust flow at boundary to match momentum for suction case
c
c     if(cmassw.gt.0) go to 173
        relax = 0.1
cmom    fact =  sqrt(abs(cmassw/rqa ))
cvn     fact =  abs(cmassw/rqb )
        fact =  cmassw/rqb
c       write (81,15) n2beg,n2end,rqb,rqa,afanc,fact
        do j= n1beg, n1end
          do k= n2beg, n2end
c         w(i,j,k,2)  = (2.*fact-1.)*w(i+1,j,k,2)
c         w(i,j,k,3)  = (2.*fact-1.)*w(i+1,j,k,3)
c         w(i,j,k,4)  = (2.*fact-1.)*w(i+1,j,k,4)
c
          w(i,j,k,1)  = w(i+1,j,k,1)
          p(i,j,k)    = p(i+1,j,k)
c
c overspecify the field similar to CFL3D
c
c         w(i,j,k,1)   = 2.0*0.77
c         w(i,j,k,2)   = -2.0*0.1902*sqrtgam*.77
c         w(i,j,k,3)   = -2.0*0.1102*sqrtgam*.77 
c
          w(i,j,k,2)   = 2.0*fact*ri1(j,k,1) - w(i+1,j,k,2)
          w(i,j,k,3)   = 2.0*fact*ri1(j,k,2) - w(i+1,j,k,3)
          w(i,j,k,4)   = 2.0*fact*ri1(j,k,3) - w(i+1,j,k,4)
c         w(i,j,k,4)   = 0.
c 
c specify temp. for inflow
c
c via pressure specification
c         p(i,j,k)     = 0.763 * 2.0 - p(i+1,j,k)
c(9-4-02)     p(i,j,k)     = w(i,j,k,1)*pchambr
c
c         if(cmassw.gt.0.0000001)
c    .    w(i,j,k,1)   = p(i,j,k)/pchambr
          w(i,j,k,1)   = p(i,j,k)/pchambr

          vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                      w(i,j,k,4)**2
          w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
c
cvn       p(i,j,k)     = ( w(i,j,k,5) - .5*vtots / w(i,j,k,1) )/dgm1
          enddo
c
c         fix the turb. quantities to free-stream
c
          do k= n2beg, n2end
            eomu(i,j,k)   =  evinf
            eomu(i,j,k)   =  eomu(i+1,j,k)
          enddo
c
          if(iturb.eq.2) then
            do k= n2beg, n2end
            turv1(i,j,k)  =  anutinf
            turv1(i,j,k)  =  turv1(i+1,j,k)
            enddo
          endif
          if(iturb.eq.3) then
            do k= n2beg, n2end
            turv1(i,j,k)  =  zkinf
            turv2(i,j,k)  =  ominf
c           turv1(i,j,k)  =  turv1(i+1,j,k)
c           turv2(i,j,k)  =  turv2(i+1,j,k)
            enddo
          endif
        enddo
 173  continue 
cif   endif
c
  181 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         i   =  imn - 1
c
         do 180 j= n1beg, n1end
         do 180 k= n2beg, n2end
         w(i,j,k,1)  = 3.*w(i+1,j,k,1) -3.*w(i+2,j,k,1) +w(i+3,j,k,1)
         w(i,j,k,2)  = 3.*w(i+1,j,k,2) -3.*w(i+2,j,k,2) +w(i+3,j,k,2)
         w(i,j,k,3)  = 3.*w(i+1,j,k,3) -3.*w(i+2,j,k,3) +w(i+3,j,k,3)
         w(i,j,k,4)  = 3.*w(i+1,j,k,4) -3.*w(i+2,j,k,4) +w(i+3,j,k,4)
         w(i,j,k,5)  = 3.*w(i+1,j,k,5) -3.*w(i+2,j,k,5) +w(i+3,j,k,5)
         p(i,j,k)    = 3.*p(i+1,j,k)   -3.*p(i+2,j,k)   +p(i+3,j,k)
         p(i,j,k)    =    max(p(i,j,k),0.001)
         eomu(i,j,k) = eomu(i+1,j,k)
ccc      turv1(i,j,k)=turv1(i+1,j,k)
  180    continue
      endif
      return
c
c
  200 continue
c--------------------- i = imx  face ----------------------------
c
c----- mods by vatsa (10-3-94) to compute only interior cells
      if( n1beg .ne. jmn) n1beg = n1beg +1
      if( n1end .eq. jmx) n1end = n1end +1
      if( n2beg .ne. kmn) n2beg = n2beg +1
      if( n2end .eq. kmx) n2end = n2end +1
c
c
      if(cmassw.gt.0.) then
        write (iwrit,'("suction/injection  not coded for imx face",
     .                 2x,"cmassw = ",e15.6)')
        call ERREXIT (nodes)
      endif
c
      if(abs(cporous).gt.0.) then
        write (iwrit,'("porosity  not coded for imx face",
     .                 2x,"cporous = ",e15.6)')
        call ERREXIT (nodes)
      endif
c
      i   =   imx + 1
c
      if (extyp .lt. 0.) go to 250
      if (extyp .gt. 0.) go to 240
c
c--------     zero order extrapolation for surface pressure
      do 230 k= n2beg, n2end
      do 230 j= n1beg, n1end
      dp(j,k)      = 0.
 230  continue
      go to 269
c
c--------     linear extrapolation for surface pressure
 240  continue
      do 245 k= n2beg, n2end
      do 245 j= n1beg, n1end
      dp(j,k)      = p(i-2,j,k) - p(i-1,j,k)
 245  continue
      go to 269
 250  continue
c
c--------- extrapolation using normal momentum equation at surface
c
      do 255 k= n2beg, n2end
      kp        = k + 1
      km        = k - 1
      do 255 j= n1beg, n1end
      jp        = j + 1
      jm        = j - 1
      skxa      = (skx(i-1,j,k) + skx(i-1,j,k-1))
      skya      = (sky(i-1,j,k) + sky(i-1,j,k-1))
      skza      = (skz(i-1,j,k) + skz(i-1,j,k-1))
      sjxa      = (sjx(i-1,j,k) + sjx(i-1,j-1,k))
      sjya      = (sjy(i-1,j,k) + sjy(i-1,j-1,k))
      sjza      = (sjz(i-1,j,k) + sjz(i-1,j-1,k))
      sixa      = 2.*six(i,j,k)
      siya      = 2.*siy(i,j,k)
      siza      = 2.*siz(i,j,k)
      rxk       = .5*(ri2(j,kp,1) - ri2(j,km,1))
      ryk       = .5*(ri2(j,kp,2) - ri2(j,km,2))
      rzk       = .5*(ri2(j,kp,3) - ri2(j,km,3))
      pk        = .5*(p(i-1,j,kp)   - p(i-1,j,km))
      rxj       = .5*(ri2(jp,k,1) - ri2(jm,k,1))
      ryj       = .5*(ri2(jp,k,2) - ri2(jm,k,2))
      rzj       = .5*(ri2(jp,k,3) - ri2(jm,k,3))
      pj        = .5*(p(i-1,jp,k)   - p(i-1,jm,k))
      si        = ri2(j,k,1)*sixa + ri2(j,k,2)*siya + ri2(j,k,3)*siza 
      sj        = ri2(j,k,1)*sjxa + ri2(j,k,2)*sjya + ri2(j,k,3)*sjza 
      sk        = ri2(j,k,1)*skxa + ri2(j,k,2)*skya + ri2(j,k,3)*skza 
      qk        = w(i-1,j,k,2)*skxa + w(i-1,j,k,3)*skya 
     .          + w(i-1,j,k,4)*skza
      qj        = w(i-1,j,k,2)*sjxa + w(i-1,j,k,3)*sjya
     .          + w(i-1,j,k,4)*sjza
      dp(j,k)   = (qk*(w(i-1,j,k,2)*rxk + w(i-1,j,k,3)*ryk
     .                +w(i-1,j,k,4)*rzk)/ w(i-1,j,k,1) 
     .            +qj*(w(i-1,j,k,2)*rxj + w(i-1,j,k,3)*ryj
     .                +w(i-1,j,k,4)*rzj)/ w(i-1,j,k,1) 
     .            -sk*pk - sj*pj)/si
  255 continue
      if (extyp .eq. -1.) go to 269
c
c--------- improved extrapolation using quadratic fit 
c          and normal pressure gradient at the surface
c
      do 260 k= n2beg,n2end
      do 260 j= n1beg,n1end
      dp(j,k)   = .75*dp(j,k) + .25*(p(i-2,j,k) - p(i-1,j,k))
  260 continue
c
  269 continue
c
c
c
      do 270 j= n1beg, n1end
      do 270 k= n2beg, n2end
ccc   w(i,j,k,1)   =    w(i-1,j,k,1)
ccc   w(i,j,k,2)   =   -w(i-1,j,k,2)
ccc   w(i,j,k,3)   =   -w(i-1,j,k,3)
ccc   w(i,j,k,4)   =   -w(i-1,j,k,4)
      eomu(i,j,k)  =    eomu(i-1,j,k)
ccc   turv1(i,j,k) =   turv1(i-1,j,k)
      p(i,j,k)     =    dim(p(i-1,j,k),dp(j,k))
      rhoi         =    1./w(i-1,j,k,1)
      ugc          =   -w(i-1,j,k,2)*rhoi
      vgc          =   -w(i-1,j,k,3)*rhoi
      wgc          =   -w(i-1,j,k,4)*rhoi
      temp         =    p(i-1,j,k)*rhoi
      tsgc         =    max((2.*twotref - temp),0.8*twotref)*(1-iadb)
     .             +    iadb*temp
      p(i,j,k)     =    max(p(i,j,k),0.001)
      w(i,j,k,1)   =    p(i,j,k)/tsgc
      w(i,j,k,2)   =    w(i,j,k,1)*ugc
      w(i,j,k,3)   =    w(i,j,k,1)*vgc
      w(i,j,k,4)   =    w(i,j,k,1)*wgc
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  270 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         i   =  imx + 2
c
         do 280 j= n1beg, n1end
         do 280 k= n2beg, n2end
         w(i,j,k,1)  = 3.*w(i-1,j,k,1) -3.*w(i-2,j,k,1) +w(i-3,j,k,1)
         w(i,j,k,2)  = 3.*w(i-1,j,k,2) -3.*w(i-2,j,k,2) +w(i-3,j,k,2)
         w(i,j,k,3)  = 3.*w(i-1,j,k,3) -3.*w(i-2,j,k,3) +w(i-3,j,k,3)
         w(i,j,k,4)  = 3.*w(i-1,j,k,4) -3.*w(i-2,j,k,4) +w(i-3,j,k,4)
         w(i,j,k,5)  = 3.*w(i-1,j,k,5) -3.*w(i-2,j,k,5) +w(i-3,j,k,5)
         p(i,j,k)    = 3.*p(i-1,j,k)   -3.*p(i-2,j,k)   +p(i-3,j,k)
         eomu(i,j,k) = eomu(i-1,j,k)
ccc      turv1(i,j,k)=turv1(i-1,j,k)
  280    continue
      endif
      return
c
c
  300 continue
c--------------------- j = jmn  face ----------------------------
c
c----- mods by vatsa (10-3-94) to compute only interior cells
      if( n1beg .ne. kmn) n1beg = n1beg +1
      if( n1end .eq. kmx) n1end = n1end +1
      if( n2beg .ne. imn) n2beg = n2beg +1
      if( n2end .eq. imx) n2end = n2end +1
c     n1beg = n1beg+1
c     n2beg = n2beg+1
c
      j   =   jmn
c
      if (extyp .lt. 0.) go to 350
      if (extyp .gt. 0.) go to 340
c
c--------     zero order extrapolation for surface pressure
      do 330 k= n1beg, n1end
      do 330 i= n2beg, n2end
      dp(i,k)      = 0.
 330  continue
      go to 369
c
c--------     linear extrapolation for surface pressure
 340  continue
      do 345 k= n1beg, n1end
      do 345 i= n2beg, n2end
      dp(i,k)      = p(i,j+2,k) - p(i,j+1,k)
 345  continue
      go to 369
 350  continue
c
c--------     normal momentum equation for surface pressure
c
      do 355 k= n1beg, n1end
      kp      = k + 1
      km      = k - 1
      do 355 i= n2beg, n2end
c
      sixa      = (six(i,j+1,k)  +six(i-1,j+1,k))
      siya      = (siy(i,j+1,k)  +siy(i-1,j+1,k))
      siza      = (siz(i,j+1,k)  +siz(i-1,j+1,k))
      sjxa      = 2.*sjx(i,j,k)
      sjya      = 2.*sjy(i,j,k)
      sjza      = 2.*sjz(i,j,k)
      skxa      = (skx(i,j+1,k)  +skx(i,j+1,k-1))
      skya      = (sky(i,j+1,k)  +sky(i,j+1,k-1))
      skza      = (skz(i,j+1,k)  +skz(i,j+1,k-1))
      rxi       = .5*(rj1(i+1,k,1)   -rj1(i-1,k,1)) 
      ryi       = .5*(rj1(i+1,k,2)   -rj1(i-1,k,2)) 
      rzi       = .5*(rj1(i+1,k,3)   -rj1(i-1,k,3)) 
      pi        = .5*(p(i+1,j+1,k)     -p(i-1,j+1,k))
      rxk       = .5*(rj1(i,kp,1)    -rj1(i,km,1))
      ryk       = .5*(rj1(i,kp,2)    -rj1(i,km,2))
      rzk       = .5*(rj1(i,kp,3)    -rj1(i,km,3))
      pk        = .5*(p(i,j+1,kp)      -p(i,j+1,km))
      si        = rj1(i,k,1)*sixa  +rj1(i,k,2)*siya  +rj1(i,k,3)*siza 
      sj        = rj1(i,k,1)*sjxa  +rj1(i,k,2)*sjya  +rj1(i,k,3)*sjza 
      sk        = rj1(i,k,1)*skxa  +rj1(i,k,2)*skya  +rj1(i,k,3)*skza 
      qi        = w(i,j+1,k,2)*sixa  +w(i,j+1,k,3)*siya
     .          + w(i,j+1,k,4)*siza
      qk        = w(i,j+1,k,2)*skxa  +w(i,j+1,k,3)*skya 
     .          + w(i,j+1,k,4)*skza
      dp(i,k)   = (qi*(w(i,j+1,k,2)*rxi  +w(i,j+1,k,3)*ryi
     .                +w(i,j+1,k,4)*rzi) /w(i,j+1,k,1) 
     .            +qk*(w(i,j+1,k,2)*rxk  +w(i,j+1,k,3)*ryk
     .                +w(i,j+1,k,4)*rzk) /w(i,j+1,k,1) 
     .            -si*pi  -sk*pk)/sj
  355 continue
      if (extyp.eq.-1.) go to 369
c
c--------- improved extrapolation using quadratic fit 
c          and normal pressure gradient at the surface
c
      do 360 k= n1beg, n1end
      do 360 i= n2beg,n2end
      dp(i,k)   = .75*dp(i,k) + .25*(p(i,j+2,k) - p(i,j+1,k))
  360 continue
c
c
  369 continue
c
c     modified for normal velocity sepcification (oct., 1999)
c
      do k= n1beg, n1end
        do i= n2beg, n2end
          p(i,j,k)     =    dim(p(i,j+1,k),dp(i,k))
          p(i,j,k)     =    max(p(i,j,k),0.001)
c
          w(i,j,k,1)   =    w(i,j+1,k,1)
        enddo
      enddo
c
c     FFT curve fit from PIV velocity data (with phase shift) at jet exit
c     (July 2004)
c
      if(iunsteady.gt.0.and.abs(tperiod).gt.0.0000001) then
c
        tratio = totptim/tperiod

        a1     = 0.7715591E+01
        a2     = 0.3625974E+00
        a3     = -0.2882919E+01
        a4     = -0.5151815E+00
        a5     = -0.8232565E+00
        b1     = 0.2149226E+02
        b2     = 0.5741831E+01
        b3     = 0.1585448E+01
        b4     = 0.6679176E+00
        b5     = 0.6231386E+00

        sum1   =      a1*cos(2.*pic*tratio) + b1*sin(2.*pic*tratio)
     .         +      a2*cos(4.*pic*tratio) + b2*sin(4.*pic*tratio)
     .         +      a3*cos(6.*pic*tratio) + b3*sin(6.*pic*tratio)
     .         +      a4*cos(8.*pic*tratio) + b4*sin(8.*pic*tratio)
c       convert to m/sec
        cmassw = 2.*pic*sum1*cmassmx*447.5/1000.
      endif
c
      if(abs(cmassw).ge.0.0000001) go to 371
      if(cporous.ge.0.0000001)     go to 376
      if(cporous.lt.0.)            go to 385
c
      do 370 k= n1beg, n1end
      do 370 i= n2beg, n2end
ccc   w(i,j,k,1)   =    w(i,j+1,k,1)
ccc   w(i,j,k,2)   =   -w(i,j+1,k,2)
ccc   w(i,j,k,3)   =   -w(i,j+1,k,3)
ccc   w(i,j,k,4)   =   -w(i,j+1,k,4)
      eomu(i,j,k)  =    eomu(i,j+1,k)
ccc   turv1(i,j,k) =   turv1(i,j+1,k)
      p(i,j,k)     =    dim(p(i,j+1,k),dp(i,k))
      rhoi         =    1./w(i,j+1,k,1)
      ugc          =   -w(i,j+1,k,2)*rhoi
      vgc          =   -w(i,j+1,k,3)*rhoi
      wgc          =   -w(i,j+1,k,4)*rhoi
      temp         =    p(i,j+1,k)*rhoi
      tsgc         =    max((2.*twotref - temp),0.8*twotref)*(1-iadb)
     .             +    iadb*temp
      p(i,j,k)     =    max(p(i,j,k),0.001)
      w(i,j,k,1)   =    p(i,j,k)/tsgc
      w(i,j,k,2)   =    w(i,j,k,1)*ugc
      w(i,j,k,3)   =    w(i,j,k,1)*vgc
      w(i,j,k,4)   =    w(i,j,k,1)*wgc
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  370 continue
c
      go to 381
c
c     mass transfer (blowing if cmassw>0, suction if cmassw<0) normal to wall
c
  371 continue
      vnorm        = sqrt(gamma)*rm*cmassw
      do 372 k= n1beg, n1end
      do 372 i= n2beg, n2end
      w(i,j,k,1)   = p(i,j,k)/pchambr
      w(i,j,k,2)   = 2.*vnorm*rj1(i,k,1) - w(i,j+1,k,2)
      w(i,j,k,3)   = 2.*vnorm*rj1(i,k,2) - w(i,j+1,k,3)
      w(i,j,k,4)   = 2.*vnorm*rj1(i,k,3) - w(i,j+1,k,4)
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  372 continue
      go to 381
c
c     wall porosity condition (Darcy's law)
c
  376 continue
      if(igrid.eq.isoln) init = init + 1
      cvnorm        = sqrt(gamma)*rm*cporous
c
c     compute the chamber pressure in cavity as the average
c     of surface pressure
c
      psurf   = 0.
      stot    = 0.
c
      do k= n1beg, n1end
         do i = n2beg, n2end
         sjt   = sqrt(sjx(i,j,k)**2 + sjy(i,j,k)**2 + sjz(i,j,k)**2)
         psurf = psurf + sjt*p(i,j+1,k)
         stot  = stot  + sjt
         enddo
      enddo
c
      if (stot.le.0.) then
c       write (iwrit,'("surface area on the porous wall is collapsed",
c    .                 2x,"stot = ",e15.6)')
        call ERREXIT (nodes)
      endif
      p3      =  psurf/stot
      iterp  = 0
  383 continue
      iterp  = iterp + 1
      pchambr =  p3
c
      do 378 k= n1beg, n1end
      do 378 i= n2beg, n2end
      vnorm        = cvnorm*(pchambr - .5*(p(i,j,k)+p(i,j+1,k)))
c     if( (pchambr - .5*(p(i,j,k)+p(i,j+1,k))).gt.0.) vnorm=3.*vnorm
      w(i,j,k,2)   = 2.*vnorm*rj1(i,k,1) - w(i,j+1,k,2)
      w(i,j,k,3)   = 2.*vnorm*rj1(i,k,2) - w(i,j+1,k,3)
      w(i,j,k,4)   = 2.*vnorm*rj1(i,k,3) - w(i,j+1,k,4)
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  378 continue
      stot  = 0.
      psurf = 0.
      amass = 0.
      do k= n1beg, n1end
         do i = n2beg, n2end
         sjt   = sqrt(sjx(i,j,k)**2 + sjy(i,j,k)**2 + sjz(i,j,k)**2)
c        psurf = psurf + sjt*p(i,j+1,k)
         psurf = psurf + sjt*(p(i,j+1,k)+p(i,j,k))*0.5
         ravg    = .5*(w(i,j,k,1)+w(i,j+1,k,1))
         uavg    = .5*(w(i,j,k,2)+w(i,j+1,k,2))/ravg
         vavg    = .5*(w(i,j,k,3)+w(i,j+1,k,3))/ravg
         wavg    = .5*(w(i,j,k,4)+w(i,j+1,k,4))/ravg
         vn1     =  uavg*rj1(i,k,1) + vavg*rj1(i,k,2) + wavg*rj1(i,k,3)
c        if(vn1.gt.0.) psurf = psurf + ravg*vn1*vn1*sjt
c        if(vn1.lt.0.) psurf = psurf - ravg*vn1*vn1*sjt
         amass = amass + ravg*vn1*sjt
         stot  = stot  + sjt
         enddo
      enddo
      p3      =  psurf/stot
      amass   =  amass/stot
ccc   p3      =  p3 - 0.667*amass/cvnorm
      p3      =  p3 - 0.52*amass/cvnorm
        if(igrid.eq.isoln)
     .  write (11,'("init,igrid,p3,p3o,amass",2x,2i5,4f12.5)')
     .               init,igrid,p3,pchambr,amass
c       if(iterp.le.1) go to 383
      go to 381
  385 continue
      if( igrid .eq. isoln) init = init + 1
      relax = 1.0
c     relax = 0.5
      relax = 0.25
c
c     wall porosity condition (a.k.a. R.H.Bush, AIAA Paper 97-2076)
c
      solidity     = 1. -abs(cporous)
      s            = solidity
c
c     compute the chamber pressure in cavity as the average 
c     of surface pressure
c
      psurf   = 0.
      stot    = 0.
c
      do k= n1beg, n1end
         do i = n2beg, n2end
         sjt   = sqrt(sjx(i,j,k)**2 + sjy(i,j,k)**2 + sjz(i,j,k)**2)
         psurf = psurf + sjt*p(i,j+1,k)
c     ravg    = .5*(w(i,j,k,1)+w(i,j+1,k,1))
c     uavg    = .5*(w(i,j,k,2)+w(i,j+1,k,2))/ravg
c     vavg    = .5*(w(i,j,k,3)+w(i,j+1,k,3))/ravg
c     wavg    = .5*(w(i,j,k,4)+w(i,j+1,k,4))/ravg
      ravg    =  w(i,j,k,1)
      uavg    =  w(i,j,k,2)/ravg
      vavg    =  w(i,j,k,3)/ravg
      wavg    =  w(i,j,k,4)/ravg
      vn1     =  uavg*rj1(i,k,1) + vavg*rj1(i,k,2) + wavg*rj1(i,k,3)
cc    if(vn1.gt.0.) psurf = psurf + ravg*vn1*vn1*sjt
cc    if(vn1.lt.0.) psurf = psurf - ravg*vn1*vn1*sjt
c     if(vn1.gt.0.) psurf = psurf - ravg*vn1*vn1*sjt
c     if(vn1.lt.0.) psurf = psurf + ravg*vn1*vn1*sjt
         stot  = stot  + sjt
         enddo
      enddo
c
      if (stot.le.0.) then
c       write (iwrit,'("surface area on the porous wall is collapsed",
c    .                 2x,"stot = ",e15.6)')
        call ERREXIT (nodes)
      endif
      p3      =  psurf/stot
ctemp por2-9.body
c     p3      = 0.994
ccraig
c     p3      = 1.0276
      p3      = 0.995
      p3      = 0.9898
      p3      = 0.993
c     p3      = 1.0
c     p3      = 0.992
c
c     temp (fix p3)
c     p3  = 0.9954
ccc
c       if(igrid.eq.isoln)
c    .  write (11,'("init,p3,relax",2x,i5,2f12.5)') init,p3,relax
c
      iterp  = 0
 384  continue
      iterp  = iterp + 1
      do 386 k = n1beg, n1end
      do 386 i = n2beg, n2end
c
cc    if(init.gt.20) go to 387
ccc   if(init.gt.00) go to 387
      if(init.lt.00) go to 387
      rvtots  = (w(i,j+1,k,2)**2 + w(i,j+1,k,3)**2 + w(i,j+1,k,4)**2)/
     .           w(i,j+1,k,1)
      am1sq   =  rvtots/(gamma*p(i,j+1,k))
      fm1     =  1. + 0.5*gm1*am1sq
      ts1     =  p(i,j+1,k)/w(i,j+1,k,1)
      tt1     =  ts1*fm1
      pt1     =  p(i,j+1,k)*fm1**gogm1
      rvn1    =  w(i,j+1,k,2)*rj1(i,k,1) + w(i,j+1,k,3)*rj1(i,k,2) +
     .           w(i,j+1,k,4)*rj1(i,k,3)
      vn1     =  rvn1/w(i,j+1,k,1)
      am1n    =  vn1/sqrt(ts1*gamma)
      ravg    =  w(i,j+1,k,1)
      pavg    =  p(i,j+1,k)
      tavg    =  ts1
      go to 388
c
 387  continue
      ravg    = .5*(w(i,j,k,1)+w(i,j+1,k,1))
      uavg    = .5*(w(i,j,k,2)+w(i,j+1,k,2))/ravg
      vavg    = .5*(w(i,j,k,3)+w(i,j+1,k,3))/ravg
      wavg    = .5*(w(i,j,k,4)+w(i,j+1,k,4))/ravg
      vtsq    =  uavg**2 + vavg**2 + wavg**2
      pavg    = .5*(p(i,j,k)+p(i,j+1,k))
      tavg    =  pavg/ravg
      am1sq   =  vtsq/(tavg*gamma)
      fm1     =  1. + 0.5*gm1*am1sq
      tt1     =  tavg*fm1
      pt1     =  pavg*fm1**gogm1
      vn1     =  uavg*rj1(i,k,1) + vavg*rj1(i,k,2) + wavg*rj1(i,k,3)
      rvn1    =  ravg*vn1
      am1n    =  vn1/sqrt(tavg*gamma)
      rm1s    =  am1n**2
      phi0    =  0.04137/(0.0982+s) + 0.57323 + 0.005786*(1-s)
 388  continue
c
ccc   if(abs(am1n).le.1.e-07.or.init.le.20) then
      if(abs(am1n).le.1.e-07) then
      vn2  = 0.
      p2   = p(i,j+1,k)
      rho2 = w(i,j+1,k,1)
      p2   = pavg
      rho2 = ravg
c     relax2 = 0.25
c     write(11,'("velocity too low for porous b.c. at init,i,k",
c    .            2x,3i5)')init,i,k
      go to 1012
      endif
      
c  PULLED FROM CFL3D . . .
c              the magnitude of the original normal velocity
c              is kept. The direction is given by Steve's assumption
cvn            if(p3.ge.p(i,j+1,k)) then
cvn               if(vn1.ge.0.) vn1= vn1
cvn               if(vn1.lt.0.) vn1=-vn1
cvn            else if(p3.lt.p(i,j+1,k)) then
cvn               if(vn1.ge.0.) vn1=-vn1
cvn               if(vn1.lt.0.) vn1= vn1
cvn            endif

               rvn1 = ravg*vn1       ! normal mass flow cons.
               rvn2 = rvn1
               am1n    =  vn1/sqrt(tavg*gamma)
c
c     test for inflow/outflow
c
      if(rvn1.le.0.) then
c        flow going into cavity
         pstat  = p(i,j+1,k)
      else
c        flow coming out of cavity
         pstat  = p3
      endif
c
      pt1a   = pstat*(1. + 0.5*gm1*am1n*am1n)**gogm1
      rvn2   = rvn1
c
c     start using Darryl Bonhaus logic
c
      run2    = rvn1**2
      unorm1  = vn1
      rm1s    =  am1n**2
      phi0    =  0.04137/(0.0982+s) + 0.57323 + 0.005786*(1-s)
      gp1     =  gamma + 1.
      rm2s    =  0.0
cvn        if(unorm1.gt.0.0)then
           if(unorm1.le.0.0)then
               rm3s = rm1s
cvn            phi = phi0 + 0.185*s**0.25*(0.5*gm1)**(gamma/gm1)
               fm = 1. + .5*gm1*rm3s
c
c    compute chokes mass flow
               phi = phi0 + 0.185*s**0.25*((1+.5*gm1)**gogm1-1.)
               do 3000 iter=1,30
                  fm = 1 + 0.5*gm1*rm3s
                  ff = 0.5*gp1*phi*phi*(1-s)*(1-s)
     .               * ((1+gamma*rm3s)/(1+gamma*phi*(1-s)))**2
     .               - rm3s*fm
                  fp = gamma*gp1*phi*phi*(1-s)*(1-s)
     .               * ((1+gamma*rm3s)/(1+gamma*phi*(1-s))**2)
     .               - 2*fm + 1
                  rm3s = rm3s - ff/fp
 3000          continue
               if(abs(ff/fp).gt.1.e-07) 
     .write(11,'("failed to conv in 3000,init,i,k,rm1s,rm3s,ff,fp",
     .            2x,3i5,5e15.6)')init,i,k,rm1s,rm3s,ff,fp,ff/fp
c              endif
             
               fm = 1 + 0.5*gm1*rm3s
cvn            run2m = gamma*gamma/tt1*p3*p3*rm3s*fm
               run2m =       gamma/tt1*p3*p3*rm3s*fm
               if (run2.gt.run2m) then
                  run2 = run2m
                  rm2s = 1.0
               else
cvn               rm3s = (sqrt(1+2*gm1/gamma/gamma*tt1*run2/p3/p3)-1)/gm1
                  rm3s = (sqrt(1.+2.*gm1/gamma *tt1*run2/p3/p3)-1.)/gm1
cvn               cf = run2*tt1/(gamma*gamma*p3*p3)
                  cf = run2*tt1/(      gamma*p3*p3)
      rm2s    =  0.0
                  do 3010 iter=1,20
               fm = 1. + .5*gm1*rm2s
               phi = phi0 + 0.185*s**0.25*(fm**gogm1-1.)
cvn                  phi = phi0 + 0.185*s**0.25
cvn  .                   * (0.5*gm1*rm2s)**(gamma/gm1)
                     dphi = gamma*0.0925*s**0.25*(fm)**(1/gm1)
                     fm = 1 + 0.5*gm1*rm2s
                     ff = rm2s*fm*phi*phi*(1-s)*(1-s)
     .                  * ((1+gamma*rm3s)/(1+gamma*rm2s*phi*(1-s)))**2
     .                  - cf
                     fp = ((1+gamma*rm3s)/(1+gamma*rm2s*phi*(1-s)))**2
     .                  * phi*phi*(1-s)*(1-s) * (0.5*gm1*rm2s + fm
     .                  * (1-2*gamma*(1-s)*rm2s*(phi+rm2s*dphi)
     .                  /  (1+gamma*rm2s*phi*(1-s)) + 2*rm2s/phi*dphi))
                     rm2s = rm2s - ff/fp
 3010             continue
               if(abs(ff/fp).gt.1.e-12)
     .write(11,'("failed to conv in 3010, init,i,k,cf,rm2s,ff,fp",
     .            2x,3i5,5e15.6)')init,i,k,cf,rm2s,ff,fp,ff/fp
               end if
cvn            phi = phi0 + 0.185*s**0.25*(0.5*gm1*rm2s)**(gamma/gm1)
               fm = 1 + 0.5*gm1*rm2s
               phi = phi0 + 0.185*s**0.25*(fm**gogm1-1.)
               pt1 = p3*((1+gamma*rm3s)/(1+gamma*rm2s*phi*(1-s)))
     .             * fm**(gamma/gm1)
               rm1s = rm3s
cvn            cf = run2*tt1/(gamma*gamma*pt1*pt1)
               cf = run2*tt1/(gamma*      pt1*pt1)
               do 3020 iter=1,20
                  ff = rm1s*(1+0.5*gm1*rm1s)**(-gp1/gm1) - cf
                  fp = (1+0.5*gm1*rm1s)**(-gp1/gm1)
     .               * (1-0.5*gp1*rm1s/(1+0.5*gm1*rm1s))
                  rm1s = rm1s - ff/fp
 3020          continue
               if(abs(ff/fp).gt.1.e-12) then
      write(11,'("failed to conv in 3020, init,i,k,cf,rm1s,ff,fp",
     .            2x,2i5,5e15.6)')init,i,k,cf,rm1s,ff,fp,ff/fp
               endif
            else
               rm3s = rm1s
cvn            phi = phi0 + 0.185*s**0.25*(0.5*gm1)**(gamma/gm1)
               fm = 1 + 0.5*gm1*rm3s
c
c           compute choked flow mass
               phi = phi0 + 0.185*s**0.25*((1.+.5*gm1)**gogm1-1.)
               do 3100 iter=1,20
                  fm = 1 + 0.5*gm1*rm3s
                  ff = (2/gp1*fm)**(gp1/gm1)*phi*phi*(1-s)*(1-s)-rm3s
                  fp = (2/gp1*fm)**(2/gm1)*phi*phi*(1-s)*(1-s) - 1
                  rm3s = rm3s - ff/fp
 3100          continue
               if(abs(ff/fp).gt.1.e-12)
     .write(11,'("failed to conv in 3100, init,i,k,cf,rm3s,ff,fp",
     .            2x,3i5,5e15.6)')init,i,k,cf,rm3s,ff,fp,ff/fp
               fm = 1 + 0.5*gm1*rm3s
cvn            run2m = gamma*gamma/tt1*p3*p3*rm3s*fm
               run2m =       gamma/tt1*p3*p3*rm3s*fm
               if (run2.gt.run2m) then
                  run2 = run2m
                  rm2s = 1.0
               else
cvn               rm3s = (sqrt(1+2*gm1/gamma/gamma*tt1*run2/p3/p3)-1)/gm1
                  rm3s = (sqrt(1.+2.*gm1/gamma *tt1*run2/p3/p3)-1.)/gm1
cvn               cf = run2*tt1/(gamma*gamma*p3*p3)
                  cf = run2*tt1/(gamma*      p3*p3)
      rm2s    =  0.0
                  do 3110 iter=1,20
                     fm = 1 + 0.5*gm1*rm2s
                     phi = phi0 + 0.185*s**0.25*(fm**gogm1-1.)
cvn                  phi = phi0 + 0.185*s**0.25
cvn  .                   * (0.5*gm1*rm2s)**(gamma/gm1)
                     dphi = gamma*0.0925*s**0.25*(fm**(1/gm1))
                     ff = rm2s*(1+0.5*gm1*rm3s)**(2*gamma/gm1)
     .                  * fm**(-gp1/gm1)*phi*phi*(1-s)*(1-s) - cf
                     fp = (1+0.5*gm1*rm3s)**(2*gamma/gm1)
     .                  * fm**(-gp1/gm1)*phi*phi*(1-s)*(1-s)
     .                  * (1-0.5*gp1*rm2s/fm+2*rm2s/phi*dphi)
                     rm2s = rm2s - ff/fp
 3110             continue
               if(abs(ff/fp).gt.1.e-12)
     .write(11,'("failed to conv in 3110, init,i,k,cf,rm2s,ff,fp",
     .            2x,3i5,5e15.6)')init,i,k,cf,rm2s,ff,fp,ff/fp
               end if
cvn            phi = phi0 + 0.185*s**0.25*(0.5*gm1*rm2s)**(gamma/gm1)
               fm = 1 + 0.5*gm1*rm2s
               phi = phi0 + 0.185*s**0.25*(fm**gogm1-1.)
               p2 = p3*((1+0.5*gm1*rm3s)/(1+0.5*gm1*rm2s))**(gamma/gm1)
               rm1s = rm3s
cvn            cf = run2*tt1/(gamma*gamma*p2*p2)
               cf = run2*tt1/(      gamma*p2*p2)
               do 3120 iter=1,20
                  fm = 1 + 0.5*gm1*rm1s
                  ff = ((1+gamma*rm2s*phi*(1-s))/(1+gamma*rm1s))**2
     .               * rm1s*fm - cf
                  fp = ((1+gamma*rm2s*phi*(1-s))/(1+gamma*rm1s))**2
     .               * (fm*(1-2*gamma*rm1s/(1+gamma*rm1s))+0.5*gm1*rm1s)
                  rm1s = rm1s - ff/fp
 3120          continue
               if(abs(ff/fp).gt.1.e-12) then
      write(11,'("failed to conv in 3120, init,i,k,cf,rm1s,ff,fp",
     .            2x,3i5,5e15.6)')init,i,k,cf,rm1s,ff,fp,ff/fp
               endif
               p1 = p2*(1+gamma*rm2s*phi*(1-s))/(1+gamma*rm1s)
               pt1 = p1*(1 + 0.5*gm1*rm1s)**(gamma/gm1)
            end if
c
c now compute fluxes from run2, rm1s and pt1
c
            am2s = rm1s
            fm   = 1 + 0.5*gm1*am2s
            p2   = pt1*fm**(-gamma/gm1)
            ts2  = tt1/fm
            rho2 = p2/ts2
cvn         vn2  = sqrt(am2s*ts2*gamma)
            vn2  = sqrt(run2)/rho2
            if(vn1.lt.0.) vn2 = -vn2
c           vn2   = vn2*rho2
c
c end Darryl's coding
c
c    values of normal velocity & pressure from last iteration
c

      pstato  = .5*(p(i,j,k)  +p(i,j+1,k))
      ravg    = .5*(w(i,j,k,1)+w(i,j+1,k,1))
      uavg    = .5*(w(i,j,k,2)+w(i,j+1,k,2))/ravg
      vavg    = .5*(w(i,j,k,3)+w(i,j+1,k,3))/ravg
      wavg    = .5*(w(i,j,k,4)+w(i,j+1,k,4))/ravg
      vno     =  uavg*rj1(i,k,1) + vavg*rj1(i,k,2) + wavg*rj1(i,k,3)
c     vno     =  ravg*vno
c
c     underrelax the porous b.c.
c
c     relax   = 1.0
c     relax   = 0.5
c     relax   = 0.25
cvn   vn2     = relax*vn2 + (1.-relax)*vno
      p2      = relax*p2  + (1.-relax)*pstato
      rho2    = relax*rho2+ (1.-relax)*ravg
 1012 continue
c
c COMPUTE NORMAL VELOCITY CONTRIBUTION

      w(i,j,k,1)   = 2.*rho2           - w(i,j+1,k,1)
c     w(i,j,k,1)   = rho2
cvn   w(i,j,k,1)   = w(i,j+1,k,1)
c
c     new logic to negate tangential vel & impose normal vel at surface
c
      ughost   = 2.*(vn2*rj1(i,k,1)+0.)-w(i,j+1,k,2)/w(i,j+1,k,1)
      vghost   = 2.*(vn2*rj1(i,k,2)+0.)-w(i,j+1,k,3)/w(i,j+1,k,1)
      wghost   = 2.*(vn2*rj1(i,k,3)+0.)-w(i,j+1,k,4)/w(i,j+1,k,1)

      rvnjp        =(w(i,j+1,k,2)*rj1(i,k,1) +
     .               w(i,j+1,k,3)*rj1(i,k,2) +
     .               w(i,j+1,k,4)*rj1(i,k,3) ) /w(i,j+1,k,1)
      rvnxjp       = rvnjp*rj1(i,k,1)
      rvnyjp       = rvnjp*rj1(i,k,2)
      rvnzjp       = rvnjp*rj1(i,k,3)
      utan1        = w(i,j+1,k,2)/w(i,j+1,k,1) - rvnxjp
      vtan1        = w(i,j+1,k,3)/w(i,j+1,k,1) - rvnyjp
      wtan1        = w(i,j+1,k,4)/w(i,j+1,k,1) - rvnzjp

c     w(i,j,k,2)   = ughost * w(i,j,k,1) - utan1 * w(i,j,k,1)
c     w(i,j,k,3)   = vghost * w(i,j,k,1) - vtan1 * w(i,j,k,1)
c     w(i,j,k,4)   = wghost * w(i,j,k,1) - wtan1 * w(i,j,k,1)
      w(i,j,k,2)   = ughost * w(i,j,k,1)
      w(i,j,k,3)   = vghost * w(i,j,k,1) 
      w(i,j,k,4)   = wghost * w(i,j,k,1)
c

      p(i,j,k)     = 2.*p2             - p(i,j+1,k)
cvn   p(i,j,k)     =    p2
      p(i,j,k)     = max(p(i,j,k),0.001)
c     w(i,j,k,1)   = p(i,j,k)*w(i,j+1,k,1)/p(i,j+1,k)
      eomu(i,j,k)  = eomu(i,j+1,k)
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
c
  386 continue
      

      stot  = 0.
      psurf = 0.
      amass = 0.
      rsurf = 0.
      do k= n1beg, n1end
         do i = n2beg, n2end
         sjt   = sqrt(sjx(i,j,k)**2 + sjy(i,j,k)**2 + sjz(i,j,k)**2)
c        psurf = psurf + sjt*p(i,j+1,k)
         psurf = psurf + sjt*(p(i,j+1,k)+p(i,j,k))*0.5
         ravg    = .5*(w(i,j,k,1)+w(i,j+1,k,1))
         rsurf   = rsurf + sjt*ravg
         uavg    = .5*(w(i,j,k,2)+w(i,j+1,k,2))/ravg
         vavg    = .5*(w(i,j,k,3)+w(i,j+1,k,3))/ravg
         wavg    = .5*(w(i,j,k,4)+w(i,j+1,k,4))/ravg
         vn1     =  uavg*rj1(i,k,1) + vavg*rj1(i,k,2) + wavg*rj1(i,k,3)
c        if(vn1.gt.0.) psurf = psurf + ravg*vn1*vn1*sjt
c        if(vn1.lt.0.) psurf = psurf - ravg*vn1*vn1*sjt
         amass = amass + ravg*vn1*sjt
         stot  = stot  + sjt
         enddo
      enddo
      p3o     =  p3
      p3      =  relax*p3o + (1.-relax)*psurf/stot
      amass   =  amass/stot
      rsurf   =  rsurf/stot
      p3      =  p3 - amass/(sqrt(gamma)*rm)
        if(igrid.eq.isoln)
     .  write (11,'("init,igrid,p3,p3o,amass,relax,r",2x,2i5,5f12.5)') 
     .               init,igrid,p3,p3o,amass,relax,rsurf
        call flush (11)
c     if (iterp.le.1.and.init.gt.20) go to 384

  381 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         j   =  jmn - 1
c
         do 380 k= n1beg, n1end
         do 380 i= n2beg, n2end
         w(i,j,k,1)  = 3.*w(i,j+1,k,1) -3.*w(i,j+2,k,1) +w(i,j+3,k,1)
         w(i,j,k,2)  = 3.*w(i,j+1,k,2) -3.*w(i,j+2,k,2) +w(i,j+3,k,2)
         w(i,j,k,3)  = 3.*w(i,j+1,k,3) -3.*w(i,j+2,k,3) +w(i,j+3,k,3)
         w(i,j,k,4)  = 3.*w(i,j+1,k,4) -3.*w(i,j+2,k,4) +w(i,j+3,k,4)
         w(i,j,k,5)  = 3.*w(i,j+1,k,5) -3.*w(i,j+2,k,5) +w(i,j+3,k,5)
         p(i,j,k)    = 3.*p(i,j+1,k)   -3.*p(i,j+2,k)   +p(i,j+3,k)
         p(i,j,k)    =    max(p(i,j,k),0.001)
         eomu(i,j,k) = eomu(i,j+1,k)
         if(abs(cporous).le.0.) go to 380
         w(i,j,k,1)   = 1.5*w(i,j+1,k,1) - 0.5*w(i,j+2,k,1)
         w(i,j,k,2)   = 1.5*w(i,j+1,k,2) - 0.5*w(i,j+2,k,2)
         w(i,j,k,3)   = 1.5*w(i,j+1,k,3) - 0.5*w(i,j+2,k,3)
         w(i,j,k,4)   = 1.5*w(i,j+1,k,4) - 0.5*w(i,j+2,k,4)
         w(i,j,k,5)   = 1.5*w(i,j+1,k,5) - 0.5*w(i,j+2,k,5)
         p(i,j,k  )   = 1.5*p(i,j+1,k  ) - 0.5*p(i,j+2,k)
  380    continue
c
c        do 380 k= n1beg, n1end
c        do 380 i= n2beg, n2end
c        w(i,j,k,1)   =    w(i,j+3,k,1)
c        w(i,j,k,2)   =   -w(i,j+3,k,2)
c        w(i,j,k,3)   =   -w(i,j+3,k,3)
c        w(i,j,k,4)   =   -w(i,j+3,k,4)
c        w(i,j,k,5)   =    w(i,j+3,k,5)
c        eomu(i,j,k)  =    eomu(i,j+3,k)
c        turv1(i,j,k) =   turv1(i,j+3,k)
c        p(i,j,k)     =    dim(p(i,j+1,k),dp(i,k))
c        vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
c    .                     w(i,j,k,4)**2
c        w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
c 380    continue
      endif
      return
c
c
  400 continue
c--------------------- j = jmx  face ----------------------------
c
c----- mods by vatsa (10-3-94) to compute only interior cells
      if( n1beg .ne. kmn) n1beg = n1beg +1
      if( n1end .eq. kmx) n1end = n1end +1
      if( n2beg .ne. imn) n2beg = n2beg +1
      if( n2end .eq. imx) n2end = n2end +1
c
      if(cmassw.gt.0.) then
        write (iwrit,'("suction/injection  not coded for jmx face",
     .                 2x,"cmassw = ",e15.6)')
        call ERREXIT (nodes)
      endif
c
      if(abs(cporous).gt.0.) then
        write (iwrit,'("porosity  not coded for jmx face")')
        write (iwrit,'("suction/injection  not coded for jmx face",
     .                 2x,"cporous = ",e15.6)')
        call ERREXIT (nodes)
      endif
c
      j   =   jmx + 1
c
      if (extyp .lt. 0.) go to 450
      if (extyp .gt. 0.) go to 440
c
c--------     zero order extrapolation for surface pressure
      do 430 k= n1beg, n1end
      do 430 i= n2beg, n2end
      dp(i,k)      = 0.
 430  continue
      go to 469
c
c--------     linear extrapolation for surface pressure
 440  continue
      do 445 k= n1beg, n1end
      do 445 i= n2beg, n2end
      dp(i,k)      = p(i,j-2,k) - p(i,j-1,k)
 445  continue
      go to 469
 450  continue
c
c--------     normal momentum equation for surface pressure
c
      do 455 k= n1beg, n1end
      kp      = k + 1
      km      = k - 1
      do 455 i= n2beg, n2end
c
      sixa      = (six(i,j-1,k)  +six(i-1,j-1,k))
      siya      = (siy(i,j-1,k)  +siy(i-1,j-1,k))
      siza      = (siz(i,j-1,k)  +siz(i-1,j-1,k))
      sjxa      = 2.*sjx(i,j,k)
      sjya      = 2.*sjy(i,j,k)
      sjza      = 2.*sjz(i,j,k)
      skxa      = (skx(i,j-1,k)  +skx(i,j-1,k-1))
      skya      = (sky(i,j-1,k)  +sky(i,j-1,k-1))
      skza      = (skz(i,j-1,k)  +skz(i,j-1,k-1))
      rxi       = .5*(rj2(i+1,k,1)   -rj2(i-1,k,1)) 
      ryi       = .5*(rj2(i+1,k,2)   -rj2(i-1,k,2)) 
      rzi       = .5*(rj2(i+1,k,3)   -rj2(i-1,k,3)) 
      pi        = .5*(p(i+1,j-1,k)     -p(i-1,j-1,k))
      rxk       = .5*(rj2(i,kp,1)    -rj2(i,km,1))
      ryk       = .5*(rj2(i,kp,2)    -rj2(i,km,2))
      rzk       = .5*(rj2(i,kp,3)    -rj2(i,km,3))
      pk        = .5*(p(i,j-1,kp)      -p(i,j-1,km))
      si        = rj2(i,k,1)*sixa  +rj2(i,k,2)*siya  +rj2(i,k,3)*siza 
      sj        = rj2(i,k,1)*sjxa  +rj2(i,k,2)*sjya  +rj2(i,k,3)*sjza 
      sk        = rj2(i,k,1)*skxa  +rj2(i,k,2)*skya  +rj2(i,k,3)*skza 
      qi        = w(i,j-1,k,2)*sixa  +w(i,j-1,k,3)*siya
     .          + w(i,j-1,k,4)*siza
      qk        = w(i,j-1,k,2)*skxa  +w(i,j-1,k,3)*skya
     .          + w(i,j-1,k,4)*skza
      dp(i,k)   = (qi*(w(i,j-1,k,2)*rxi  +w(i,j-1,k,3)*ryi
     .                +w(i,j-1,k,4)*rzi) /w(i,j-1,k,1) 
     .            +qk*(w(i,j-1,k,2)*rxk  +w(i,j-1,k,3)*ryk
     .                +w(i,j-1,k,4)*rzk) /w(i,j-1,k,1) 
     .            -si*pi  -sk*pk)/sj
  455 continue
      if (extyp.eq.-1.) go to 469
c
c--------- improved extrapolation using quadratic fit 
c          and normal pressure gradient at the surface
c
      do 460 k= n1beg, n1end
      do 460 i= n2beg,n2end
      dp(i,k)   = .75*dp(i,k) + .25*(p(i,j-2,k) - p(i,j-1,k))
  460 continue
c
c
  469 continue
      do 470 k= n1beg, n1end
      do 470 i= n2beg, n2end
ccc   w(i,j,k,1)   =    w(i,j-1,k,1)
ccc   w(i,j,k,2)   =   -w(i,j-1,k,2)
ccc   w(i,j,k,3)   =   -w(i,j-1,k,3)
ccc   w(i,j,k,4)   =   -w(i,j-1,k,4)
      eomu(i,j,k)  =    eomu(i,j-1,k)
ccc   turv1(i,j,k) =   turv1(i,j-1,k)
      p(i,j,k)     =    dim(p(i,j-1,k),dp(i,k))
      rhoi         =    1./w(i,j-1,k,1)
      ugc          =   -w(i,j-1,k,2)*rhoi
      vgc          =   -w(i,j-1,k,3)*rhoi
      wgc          =   -w(i,j-1,k,4)*rhoi
      temp         =    p(i,j-1,k)*rhoi
      tsgc         =    max((2.*twotref - temp),0.8*twotref)*(1-iadb)
     .             +    iadb*temp
      p(i,j,k)     =    max(p(i,j,k),0.001)
      w(i,j,k,1)   =    p(i,j,k)/tsgc
      w(i,j,k,2)   =    w(i,j,k,1)*ugc
      w(i,j,k,3)   =    w(i,j,k,1)*vgc
      w(i,j,k,4)   =    w(i,j,k,1)*wgc
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  470 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         j   =  jmx + 2
c
         do 480 k= n1beg, n1end
         do 480 i= n2beg, n2end
         w(i,j,k,1)  = 3.*w(i,j-1,k,1) -3.*w(i,j-2,k,1) +w(i,j-3,k,1)
         w(i,j,k,2)  = 3.*w(i,j-1,k,2) -3.*w(i,j-2,k,2) +w(i,j-3,k,2)
         w(i,j,k,3)  = 3.*w(i,j-1,k,3) -3.*w(i,j-2,k,3) +w(i,j-3,k,3)
         w(i,j,k,4)  = 3.*w(i,j-1,k,4) -3.*w(i,j-2,k,4) +w(i,j-3,k,4)
         w(i,j,k,5)  = 3.*w(i,j-1,k,5) -3.*w(i,j-2,k,5) +w(i,j-3,k,5)
         p(i,j,k)    = 3.*p(i,j-1,k)   -3.*p(i,j-2,k)   +p(i,j-3,k)
         eomu(i,j,k) = eomu(i,j-1,k)
ccc      turv1(i,j,k)=turv1(i,j-1,k)
  480    continue
      endif
      return
c
c
  500 continue
c--------------------- k = kmn  face ----------------------------
c
c----- mods by vatsa (10-3-94) to compute only interior cells
      if( n1beg .ne. imn) n1beg = n1beg +1
      if( n1end .eq. imx) n1end = n1end +1
      if( n2beg .ne. jmn) n2beg = n2beg +1
      if( n2end .eq. jmx) n2end = n2end +1
c
c     if(cmassw.gt.0.) then
c       write (iwrit,'("suction/injection  not coded for kmn face",
c    .                 2x,"cmassw = ",e15.6)')
c       call ERREXIT (nodes)
c     endif
c
      if(abs(cporous).gt.0.) then
        write (iwrit,'("porosity  not coded for kmn face",
     .                 2x,"cporous = ",e15.6)')
        call ERREXIT (nodes)
      endif
c
      k   =   kmn
c
      if (extyp .lt. 0.) go to 550
      if (extyp .gt. 0.) go to 540
c
c--------     zero order extrapolation for surface pressure
      do 530 j= n2beg, n2end
      do 530 i= n1beg, n1end
      dp(i,j)      = 0.
 530  continue
      go to 569
c
c--------     linear extrapolation for surface pressure
 540  continue
      do 545 j= n2beg, n2end
      do 545 i= n1beg, n1end
      dp(i,j)      = p(i,j,k+2) - p(i,j,k+1)
 545  continue
      go to 569
c
 550  continue
c
c--------- extrapolation using normal momentum equation at surface
c
      do 555 j= n2beg, n2end
      jp        = j + 1
      jm        = j - 1
      do 555 i= n1beg, n1end
      ip        = i + 1
      im        = i - 1
      sixa      = (six(i,j,k+1) + six(i-1,j,k+1))
      siya      = (siy(i,j,k+1) + siy(i-1,j,k+1))
      siza      = (siz(i,j,k+1) + siz(i-1,j,k+1))
      sjxa      = (sjx(i,j,k+1) + sjx(i,j-1,k+1))
      sjya      = (sjy(i,j,k+1) + sjy(i,j-1,k+1))
      sjza      = (sjz(i,j,k+1) + sjz(i,j-1,k+1))
      skxa      = 2.*skx(i,j,k)
      skya      = 2.*sky(i,j,k)
      skza      = 2.*skz(i,j,k)
      rxi       = .5*(rk1(ip,j,1) - rk1(im,j,1))
      ryi       = .5*(rk1(ip,j,2) - rk1(im,j,2))
      rzi       = .5*(rk1(ip,j,3) - rk1(im,j,3))
      pi        = .5*(p(ip,j,k+1)   - p(im,j,k+1))
      rxj       = .5*(rk1(i,jp,1) - rk1(i,jm,1))
      ryj       = .5*(rk1(i,jp,2) - rk1(i,jm,2))
      rzj       = .5*(rk1(i,jp,3) - rk1(i,jm,3))
      pj        = .5*(p(i,jp,k+1)   - p(i,jm,k+1))
      si        = rk1(i,j,1)*sixa + rk1(i,j,2)*siya + rk1(i,j,3)*siza 
      sj        = rk1(i,j,1)*sjxa + rk1(i,j,2)*sjya + rk1(i,j,3)*sjza 
      sk        = rk1(i,j,1)*skxa + rk1(i,j,2)*skya + rk1(i,j,3)*skza 
      qi        = w(i,j,k+1,2)*sixa + w(i,j,k+1,3)*siya 
     .          + w(i,j,k+1,4)*siza
      qj        = w(i,j,k+1,2)*sjxa + w(i,j,k+1,3)*sjya
     .          + w(i,j,k+1,4)*sjza
      dp(i,j)   = (qi*(w(i,j,k+1,2)*rxi + w(i,j,k+1,3)*ryi
     .                +w(i,j,k+1,4)*rzi)/ w(i,j,k+1,1) 
     .            +qj*(w(i,j,k+1,2)*rxj + w(i,j,k+1,3)*ryj
     .                +w(i,j,k+1,4)*rzj)/ w(i,j,k+1,1) 
     .            -si*pi - sj*pj)/sk
  555 continue
      if (extyp .eq. -1.) go to 569
c
c--------- improved extrapolation using quadratic fit 
c          and normal pressure gradient at the surface
c
      do 560 j= n2beg,n2end
      do 560 i= n1beg,n1end
      dp(i,j)   = .75*dp(i,j) + .25*(p(i,j,k+2) - p(i,j,k+1))
  560 continue
c
  569 continue
c                                                             
c     modified for normal velocity sepcification (May, 2003) 
c                                                             
      do j= n2beg, n2end                                      
        do i= n1beg, n1end                                    
          p(i,j,k)     =    dim(p(i,j,k+1),dp(i,k))           
          p(i,j,k)     =    max(p(i,j,k),0.001)               
c                                                             
          w(i,j,k,1)   =    w(i,j,k+1,1)                      
        enddo                                                 
      enddo                                                   
c                                                             
      if(iunsteady.gt.0.and.abs(tperiod).gt.0.0000001)        
     .                   cmassw = cmassmx*sin(2.*pic*totptim/tperiod)        
c     write(81,16) iunsteady,nodes(ibloc),myrank,             
c    .             cmassmx,pi,totptim,tperiod,cmassw          
c16   format(2x,'iuns,cmax,pi,tot,tp,cmw',3i3,5e16.5)         
      if(abs(cmassw).ge.0.0000001) go to 571                  
c 
c
      do 570 i= n1beg, n1end
      do 570 j= n2beg, n2end
ccc   w(i,j,k,1)   =    w(i,j,k+1,1)
ccc   w(i,j,k,2)   =   -w(i,j,k+1,2)
ccc   w(i,j,k,3)   =   -w(i,j,k+1,3)
ccc   w(i,j,k,4)   =   -w(i,j,k+1,4)
      eomu(i,j,k)  =    eomu(i,j,k+1)
ccc   turv1(i,j,k) =   turv1(i,j,k+1)
      p(i,j,k)     =    dim(p(i,j,k+1),dp(i,j))
      rhoi         =    1./w(i,j,k+1,1)
      ugc          =   -w(i,j,k+1,2)*rhoi
      vgc          =   -w(i,j,k+1,3)*rhoi
      wgc          =   -w(i,j,k+1,4)*rhoi
      temp         =    p(i,j,k+1)*rhoi
      tsgc         =    max((2.*twotref - temp),0.8*twotref)*(1-iadb)
     .             +    iadb*temp
      p(i,j,k)     =    max(p(i,j,k),0.001)
      w(i,j,k,1)   =    p(i,j,k)/tsgc
      w(i,j,k,2)   =    w(i,j,k,1)*ugc
      w(i,j,k,3)   =    w(i,j,k,1)*vgc
      w(i,j,k,4)   =    w(i,j,k,1)*wgc
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  570 continue
c
      go to 581
c
c     mass transfer (blowing if cmassw>0, suction if cmassw<0) normal to wall
c
  571 continue
      vnorm        = sqrt(gamma)*rm*cmassw
c     specify momentum
c     vnorm        = sqrt(gamma)*rm*vnorm
      do 572 i= n1beg, n1end
      do 572 j= n2beg, n2end
cvn   if(cmassw.lt.0.) p(i,j,k) = pchambr
c(fix temp. of cavity)
      w(i,j,k,1)   = p(i,j,k)/pchambr
      w(i,j,k,2)   = 2.*vnorm*rk1(i,j,1) - w(i,j,k+1,2)
      w(i,j,k,3)   = 2.*vnorm*rk1(i,j,2) - w(i,j,k+1,3)
      w(i,j,k,4)   = 2.*vnorm*rk1(i,j,3) - w(i,j,k+1,4)
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
c     if(cmassw.gt.0) go to 172
c     w(i,j,k,1)  = w(i,j,k+1,1)
c     w(i,j,k,2)  = w(i,j,k+1,2)
c     w(i,j,k,3)  = w(i,j,k+1,3)
c     w(i,j,k,4)  = w(i,j,k+1,4)
c     w(i,j,k,5)  = w(i,j,k+1,5)
c     p(i,j,k)    = p(i,j,k+1)
  572 continue
c
      go to 581
cif   if( igrid .eq. isoln) then                              
c                                                             
        signx=1.0                                             
c                                                             
c       compute momentum flux at fan face                     
c                                                             
        rqa=0.0                                               
        rqb=0.0                                               
        afanc=0.0                                             
        coef=sqrt(gamma)*rm                                   
c                                                             
        do 573 i= n1beg, n1end
        do 573 j= n2beg, n2end
        rqn=-(rk1(i,j,1)*skx(i,j,kmn)+rk1(i,j,2)*sky(i,j,kmn) 
     .                               +rk1(i,j,3)*skz(i,j,kmn))/coef         
cvn     rqb=rqb-rqn*signx*.5                                  
        rqb=rqb-rqn*signx                                     
        afanc=afanc-six(i,j,k)*signx
  573   continue
        fact =  cmassw/rqb
c       write (81,15) n2beg,n2end,rqb,rqa,afanc,fact
c
        do i= n1beg, n1end
          do j= n2beg, n2end
          w(i,j,k,1)  = w(i,j,k+1,1)
          p(i,j,k)    = p(i,j,k+1)                            
c                                                             
          w(i,j,k,2)   = 2.0*fact*rk1(i,j,1) - w(i,j,k+1,2)   
          w(i,j,k,3)   = 2.0*fact*rk1(i,j,2) - w(i,j,k+1,3)   
          w(i,j,k,4)   = 2.0*fact*rk1(i,j,3) - w(i,j,k+1,4)   
c 
c specify temp. for inflow
c
c via pressure specification
c                                                             
c         if(cmassw.gt.0.0000001)                             
c    .    w(i,j,k,1)   = p(i,j,k)/pchambr                     
          w(i,j,k,1)   = p(i,j,k)/pchambr                     
                                                              
          vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +   
     .                      w(i,j,k,4)**2                     
          w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
c                                                             
cvn       p(i,j,k)     = ( w(i,j,k,5) - .5*vtots / w(i,j,k,1) )/dgm1        
          enddo
c                                                             
c         fix the turb. quantities to free-stream             
c                                                             
          do j=n2beg, n2end
            eomu(i,j,k)   =  evinf                            
            eomu(i,j,k)   =  eomu(i,j,k+1)                    
          enddo
c                                                             
          if(iturb.eq.2) then                                 
            do j= n2beg, n2end                                
c           turv1(i,j,k)  =  anutinf                          
            turv1(i,j,k)  =  turv1(i,j,k+1)                   
            enddo                                             
          endif                                               
          if(iturb.eq.3) then                                 
            do j= n2beg, n2end                                
            turv1(i,j,k)  =  zkinf                            
            turv2(i,j,k)  =  ominf                            
c           turv1(i,j,k)  =  turv1(i,j,k+1)                   
c           turv2(i,j,k)  =  turv2(i,j,k+1)                   
            enddo                                             
          endif                                               
        enddo
  581   continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         k   =  kmn - 1
c
         do 580 i= n1beg, n1end
         do 580 j= n2beg, n2end
         w(i,j,k,1)  = 3.*w(i,j,k+1,1) -3.*w(i,j,k+2,1) +w(i,j,k+3,1)
         w(i,j,k,2)  = 3.*w(i,j,k+1,2) -3.*w(i,j,k+2,2) +w(i,j,k+3,2)
         w(i,j,k,3)  = 3.*w(i,j,k+1,3) -3.*w(i,j,k+2,3) +w(i,j,k+3,3)
         w(i,j,k,4)  = 3.*w(i,j,k+1,4) -3.*w(i,j,k+2,4) +w(i,j,k+3,4)
         w(i,j,k,5)  = 3.*w(i,j,k+1,5) -3.*w(i,j,k+2,5) +w(i,j,k+3,5)
         p(i,j,k)    = 3.*p(i,j,k+1)   -3.*p(i,j,k+2)   +p(i,j,k+3)
         eomu(i,j,k) = eomu(i,j,k+1)
ccc      turv1(i,j,k)=turv1(i,j,k+1)
  580    continue
      endif
      return
c
c
  600 continue
c--------------------- k = kmx  face ----------------------------
c
c----- mods by vatsa (10-3-94) to compute only interior cells
      if( n1beg .ne. imn) n1beg = n1beg +1
      if( n1end .eq. imx) n1end = n1end +1
      if( n2beg .ne. jmn) n2beg = n2beg +1
      if( n2end .eq. jmx) n2end = n2end +1
c
      if(cmassw.gt.0.) then
        write (iwrit,'("suction/injection  not coded for kmx face",
     .                 2x,"cmassw = ",e15.6)')
        call ERREXIT (nodes)
      endif
c
      if(abs(cporous).gt.0.) then
        write (iwrit,'("porosity  not coded for kmx face")')
        write (iwrit,'("suction/injection  not coded for kmx face",
     .                 2x,"cporous = ",e15.6)')
        call ERREXIT (nodes)
      endif
c
      k   =   kmx + 1
c
      if (extyp .lt. 0.) go to 650
      if (extyp .gt. 0.) go to 640
c
c--------     zero order extrapolation for surface pressure
      do 630 j= n2beg, n2end
      do 630 i= n1beg, n1end
      dp(i,j)      = 0.
 630  continue
      go to 669
c
c--------     linear extrapolation for surface pressure
 640  continue
      do 645 j= n2beg, n2end
      do 645 i= n1beg, n1end
      dp(i,j)      = p(i,j,k-2) - p(i,j,k-1)
 645  continue
      go to 669
c
 650  continue
c
c--------- extrapolation using normal momentum equation at surface
c
      do 655 j= n2beg, n2end
      jp        = j + 1
      jm        = j - 1
      do 655 i= n1beg, n1end
      ip        = i + 1
      im        = i - 1
      sixa      = (six(i,j,k-1) + six(i-1,j,k-1))
      siya      = (siy(i,j,k-1) + siy(i-1,j,k-1))
      siza      = (siz(i,j,k-1) + siz(i-1,j,k-1))
      sjxa      = (sjx(i,j,k-1) + sjx(i,j-1,k-1))
      sjya      = (sjy(i,j,k-1) + sjy(i,j-1,k-1))
      sjza      = (sjz(i,j,k-1) + sjz(i,j-1,k-1))
      skxa      = 2.*skx(i,j,k)
      skya      = 2.*sky(i,j,k)
      skza      = 2.*skz(i,j,k)
      rxi       = .5*(rk2(ip,j,1) - rk2(im,j,1))
      ryi       = .5*(rk2(ip,j,2) - rk2(im,j,2))
      rzi       = .5*(rk2(ip,j,3) - rk2(im,j,3))
      pi        = .5*(p(ip,j,k-1)   - p(im,j,k-1))
      rxj       = .5*(rk2(i,jp,1) - rk2(i,jm,1))
      ryj       = .5*(rk2(i,jp,2) - rk2(i,jm,2))
      rzj       = .5*(rk2(i,jp,3) - rk2(i,jm,3))
      pj        = .5*(p(i,jp,k-1)   - p(i,jm,k-1))
      si        = rk2(i,j,1)*sixa + rk2(i,j,2)*siya + rk2(i,j,3)*siza 
      sj        = rk2(i,j,1)*sjxa + rk2(i,j,2)*sjya + rk2(i,j,3)*sjza 
      sk        = rk2(i,j,1)*skxa + rk2(i,j,2)*skya + rk2(i,j,3)*skza 
      qi        = w(i,j,k-1,2)*sixa + w(i,j,k-1,3)*siya 
     .          + w(i,j,k-1,4)*siza
      qj        = w(i,j,k-1,2)*sjxa + w(i,j,k-1,3)*sjya
     .          + w(i,j,k-1,4)*sjza
      dp(i,j)   = (qi*(w(i,j,k-1,2)*rxi + w(i,j,k-1,3)*ryi
     .                +w(i,j,k-1,4)*rzi)/ w(i,j,k-1,1) 
     .            +qj*(w(i,j,k-1,2)*rxj + w(i,j,k-1,3)*ryj
     .                +w(i,j,k-1,4)*rzj)/ w(i,j,k-1,1) 
     .            -si*pi - sj*pj)/sk
  655 continue
      if (extyp .eq. -1.) go to 669
c
c--------- improved extrapolation using quadratic fit 
c          and normal pressure gradient at the surface
c
      do 660 j= n2beg,n2end
      do 660 i= n1beg,n1end
      dp(i,j)   = .75*dp(i,j) + .25*(p(i,j,k-2) - p(i,j,k-1))
  660 continue
c
c
  669 continue
c
      do 670 i= n1beg, n1end
      do 670 j= n2beg, n2end
      w(i,j,k,1)   =    w(i,j,k-1,1)
      w(i,j,k,2)   =   -w(i,j,k-1,2)
      w(i,j,k,3)   =   -w(i,j,k-1,3)
      w(i,j,k,4)   =   -w(i,j,k-1,4)
c     w(i,j,k,5)   =    w(i,j,k-1,5)
      eomu(i,j,k)  =    eomu(i,j,k-1)
ccc   turv1(i,j,k) =   turv1(i,j,k-1)
      p(i,j,k)     =    dim(p(i,j,k-1),dp(i,j))
      rhoi         =    1./w(i,j,k-1,1)
      ugc          =   -w(i,j,k-1,2)*rhoi
      vgc          =   -w(i,j,k-1,3)*rhoi
      wgc          =   -w(i,j,k-1,4)*rhoi
      temp         =    p(i,j,k-1)*rhoi
      tsgc         =    max((2.*twotref - temp),0.8*twotref)*(1-iadb)
     .             +    iadb*temp
      p(i,j,k)     =    max(p(i,j,k),0.001)
      w(i,j,k,1)   =    p(i,j,k)/tsgc
      w(i,j,k,2)   =    w(i,j,k,1)*ugc
      w(i,j,k,3)   =    w(i,j,k,1)*vgc
      w(i,j,k,4)   =    w(i,j,k,1)*wgc
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
  670 continue
c
      if( igrid .eq. isoln) then
c
c        second layer of ghost cells for fine grids
c
         k   =  kmx + 2
c
         do 680 i= n1beg, n1end
         do 680 j= n2beg, n2end
         w(i,j,k,1)  = 3.*w(i,j,k-1,1) -3.*w(i,j,k-2,1) +w(i,j,k-3,1)
         w(i,j,k,2)  = 3.*w(i,j,k-1,2) -3.*w(i,j,k-2,2) +w(i,j,k-3,2)
         w(i,j,k,3)  = 3.*w(i,j,k-1,3) -3.*w(i,j,k-2,3) +w(i,j,k-3,3)
         w(i,j,k,4)  = 3.*w(i,j,k-1,4) -3.*w(i,j,k-2,4) +w(i,j,k-3,4)
         w(i,j,k,5)  = 3.*w(i,j,k-1,5) -3.*w(i,j,k-2,5) +w(i,j,k-3,5)
         p(i,j,k)    = 3.*p(i,j,k-1)   -3.*p(i,j,k-2)   +p(i,j,k-3)
         eomu(i,j,k) = eomu(i,j,k-1)
ccc      turv1(i,j,k)=turv1(i,j,k-1)
  680    continue
      endif
c
c
      return
      end
