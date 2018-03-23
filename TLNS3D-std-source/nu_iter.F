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
      subroutine nu_iter(rhow,pw,tw,visl,uplusc,tnuslip,vist )
c
c     Object: To solve for nutidel of S-A turb. model using
c             Newton Raphson Scheme
c
c     Initial Coding: By V.N. vatsa of NASA Langley (Jan. 2003)
c     Patterned after USM3D
c
      implicit none
c
      integer i2dfl,iturb
c
      real gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .     ca,sa,pr,prt,rey,suthc,tref
c
      real rhow,pw,tw,visl,uplusc,tnuslip
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c     local variables
c
      integer icnt,iter
c
      real akarman,cv1,tbeta,expbk,sfac,sfaci,akupl,expkupl,vist,
     .     anut,anul,tnuold,chi,chi3,fv1,ftnu,ftnup
c
c     iniitialization
c
      icnt    =  0
      iter    =  0
      tnuslip =  3.0
      akarman =  0.41
      cv1     =  7.1
      tbeta   =  5.
      expbk   =  exp(-tbeta*akarman)
      sfac    =  sqrt(gamma)*rm/rey
      sfaci   =  1./sfac
c
      akupl   =  akarman*uplusc
      expkupl =  exp(akupl)
      vist    =  visl*akarman*expbk*(expkupl - 1. - akupl
     .                                       - 0.50*akupl*akupl)
      if(vist.le.0..or.vist.gt.300)
     .write (6,'("vist in nu_iter unbounded",5e15.6)') uplusc,vist,
     .                                                 visl
      vist    =  max(visl,vist)
      anut    =  vist/rhow 
      anul    =  visl/rhow 
      
c
c
 110  continue
c
c     limit the value of ustar
c
c     tnuold  =  min(tnuslip,100.)
      tnuold  =  min(tnuslip,150.)
c
      chi     =  tnuold/anul
      chi3    =  chi**3
      fv1     =  chi3/(chi3 + cv1**3)
c
      ftnu    =  tnuold*fv1 - anut
      ftnup   =  fv1 + 3.0*chi3*(cv1**3)/((chi3+cv1**3)**2)

c 
      tnuslip = tnuold - ftnu/ftnup
c
      iter    = iter + 1
c
      if(iter.gt.100) go to 112
      if(abs(tnuslip-tnuold).gt.1.e-9) go to 110
      go to 111
c
 112  continue
c
      write (6,'("nu_iter failed to converge",5e15.6)') tnuslip,tnuold,
     .            ftnu,ftnup,uplusc
 111  continue

      return
      end
