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
      subroutine us_iter(utot,ynorm,rhow,pw,tw,visl,ustar,tauw,uplusc,
     .                   yplusc,i,j,k)
c
c     Object: To solve for ustar using Newton Raphson Scheme
c             and Eq. 6-62 of White
c
c     Initial Coding: By V.N. vatsa of NASA Langley (Jan. 2003)
c     Patterned after USM3D
c
      implicit none
c
      integer i2dfl,iturb
      integer i,j,k
c
      real gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .     ca,sa,pr,prt,rey,suthc,tref
c
      real utot,ynorm,rhow,pw,tw,visl,ustar,tauw,uplusc,yplusc
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c     local variables
c
      integer icnt,iter
c
      real akarman,akumax,tbeta,expbk,sfac,sfaci,uwf,duwf,usold,usoldi,
     .     usisq,akupl,akuplp,expkupl,fus,fusp 
c
c     iniitialization
c
      icnt    =  0
      iter    =  0
      ustar   =  0.05*rm
      akarman =  0.41
      akumax  =  4000.
      tbeta   =  5.
      expbk   =  exp(-tbeta*akarman)
      sfac    =  sqrt(gamma)*rm/rey
      sfaci   =  1./sfac
c
c     underrelaxation parameters
      uwf     =  0.1
      duwf    =  0.1 
c
 110  continue
      uwf     =  min(0.5,uwf + duwf)
c
c     limit the value of ustar
c
      usold   =  min(ustar,rm)
      usoldi  =  1./usold
      akupl   =  min(akarman*utot*usoldi, akumax)
      usoldi  =  akupl/(akarman*utot)
      usold   =  1./usoldi
      usisq   =  usoldi*usoldi
      akuplp  = -akarman*utot*usisq
      expkupl =  exp(akupl)
c
      fus     =  utot*usoldi
     .        +  expbk*(expkupl - 1. - akupl - 0.5*akupl*akupl
     .                               - akupl*akupl*akupl/6.  )
     .        -  rhow*ynorm*usold/(visl*sfac) 

      fusp    = -utot*usisq
     .        +  expbk*(expkupl*akuplp - (1.+akupl)*akuplp 
     .                                 - 0.5*akuplp*akupl**2 )
     .        -  rhow*ynorm/(visl*sfac)

c 
      ustar   = usold - uwf*fus/fusp
c
      iter    = iter + 1
c
      if(iter.gt.100) go to 112
      if(abs(ustar-usold).gt.1.e-9) go to 110
      go to 111
c
 112  continue
c
      write (6,'("us_iter failed to converge",3i4,7e15.6)')i,j,k,
     .            ustar,usold,fus,fusp,utot,ynorm,visl
 111  continue

      if(abs(ustar).le.1.e-12 ) ustar = 1.e-12
      tauw    = rhow*ustar*ustar*sfaci
      uplusc  = utot/ustar
      yplusc  = rhow*ynorm*ustar/(visl*sfac)

      return
      end
      subroutine us2_iter(utot,ynorm,rhow,pw,tw,visl,ustar,tauw,uplusc,
     .                   yplusc,dpdx,i,j,k)
c
c     Object: To solve for ustar using Newton Raphson Scheme
c             and Eq. 24 of Wilcox AIAA Paper89-0611
c
c     Initial Coding: By V.N. vatsa of NASA Langley (Feb. 2003)
c     Patterned after USM3D
c
      implicit none
c
      integer i2dfl,iturb
      integer i,j,k
      integer init
      data init/0/
c
      real gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .     ca,sa,pr,prt,rey,suthc,tref
c
      real utot,ynorm,rhow,pw,tw,visl,ustar,tauw,uplusc,yplusc,
     .     dpdx,phi,term1,termdp,phiyp
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c     local variables
c
      integer icnt,iter
c
      real akarman,akumax,tbeta,sfac,sfaci,uwf,duwf,usold,
     .     usisq,fus,fusp 
c
c     iniitialization
c
      init    = init+1
      icnt    =  0
      iter    =  0
      ustar   =  0.05*rm
      akarman =  0.41
      akumax  =  4000.
      tbeta   =  5.
      sfac    =  sqrt(gamma)*rm/rey
      sfaci   =  1./sfac
c
c     underrelaxation parameters
      uwf     =  0.1
      duwf    =  0.1 
c
      usold   =  min(ustar,rm)
 110  continue
      uwf     =  min(0.4,uwf + duwf)
c     uwf     =  min(0.5,uwf + duwf)
c     uwf     =  min(1.0,uwf + duwf)
c
c     limit the value of ustar
c
      usold   =  abs(ustar)
      uplusc  =  utot/usold
c     uplusc  =  min(utot/usold,17.5)
c     uplusc  =  max(uplusc,3.)
c     usold   =  utot/uplusc
cvn   usold   =  min(usold,rm)
      yplusc  =  rhow*ynorm*usold/(visl*sfac)
c     limit phi value
      phi     =  sfac*visl*dpdx/(rhow*rhow*usold**3)
      phiyp   =  phi*yplusc
      phiyp   =  sign(min (0.50,abs(phiyp)),phiyp)
c     phiyp   =  sign(min (0.22,abs(phiyp)),phiyp)
      phi     =  phiyp/yplusc
c
c     phi     = 0.
c
      term1   =  log(rhow*ynorm*usold/(visl*sfac))/akarman+tbeta
      termdp  =  1.13*phi*rhow*ynorm*usold*sfaci/(visl*akarman)
c     termdp  =  0.48*phi*rhow*ynorm*usold*sfaci/(visl*akarman)
      if(termdp.gt.0.) termdp = min(0.99*term1,termdp)
      fus     =  usold*(term1 - termdp) - utot
c     fus     =  usold*(log(rhow*ynorm*usold/(visl*sfac))/akarman+tbeta
c    .                  -termdp) 
c    .                  -0.48*phi*rhow*ynorm*usold*sfaci/(visl*akarman))
c    .                  - utot

      fusp    = term1
     .        - termdp + 1.
c     fusp    = log(rhow*ynorm*usold/(visl*sfac))/akarman + tbeta 
c    .        + usold*rhow*ynorm*sfaci/
c    .         (visl*akarman*rhow*ynorm*usold/(visl*sfac) )
c    .        + 0.48*phi*rhow*ynorm*usold*sfaci/(visl*akarman)
c    .         (visl*akarman*log(rhow*ynorm*usold/(visl*sfac)) )
c    .        - 2.0*0.48*phi*rhow*ynorm*usold*sfaci/(visl*akarman)

c 
      ustar   = usold - uwf*fus/fusp
c
      iter    = iter + 1
c
      if(iter.gt.100) go to 112
      if(abs(ustar-usold).gt.1.e-9) go to 110
      go to 111
c
 112  continue
c
      if(abs(ustar-usold).gt.1.e-4)
     .write (6,'("us_iter failed to converge",3i4,8e14.5)')i,j,k,
     .            ustar,usold,fus,fusp,utot,phiyp,phi,termdp
 111  continue
      if(init.eq.300) 
     .write (6,'("us_iter dpdx,phi",3i4,8e14.5)')i,j,k,
     .            ustar,usold,fus,fusp,sfac,rey,dpdx,phi
      if(abs(ustar).le.1.e-12 ) ustar = 1.e-12
      tauw    = rhow*ustar*ustar*sfaci
c     uplusc  = utot/ustar
      uplusc  = min(utot/ustar,17.5)
      uplusc  = max(uplusc,3.)
      yplusc  = rhow*ynorm*ustar/(visl*sfac)

      return
      end
