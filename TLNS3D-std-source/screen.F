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
      subroutine screen (ru,tt,pt,m,gamma,r,s,pt2,ru2,j,nodes)
c*********************************************************************
c     purpose: calculates total pressure losses for flow through 
c              screens. 
c     Provided by S. Bauer (from R. Bush)
c
c     Latest mods: By V.N.Vatsa for TLNS3D-MB (Oct. 1999)
c*********************************************************************
 
      dimension nodes(1)

      real m,ms1,ms2,ms3,ms4
      iflg=0
      idebug = 0

c     check for no flow through screen
c     write(6,*)ru
      if ( abs(ru) .lt. 1.e-12 .or. s.eq.0. ) then
        pt2 = pt
        ru2 = ru
        return
      end if

c     constants needed later
      gami = gamma - 1.
      fg1 = gamma/gami

c     incompressible contraction coefficient for this solidity
c       limits: s=0;phi0=1.    s=1.;phi0=.611
      phi0 = .041366/(1.098216 + s - 1.) + 0.573265 + .005786*(1.-s)

c     recompute mach to ensure consistency, (use incoming m to define branch)
      ms1 = m*m
      if(idebug.gt.0)
     .write(46,*)j,ms1,ru

c----------------------------------------------------------------------------
c     corrfl commented out because the resultant ms1 is always slightly
c     larger than the starting value. Goes nuts when things start to degrade
c------------------------------------------------------AFT012199------------

c     call corrfl(ru,tt,pt,gamma,r,0.,ms1,ru2,j)
      if(idebug.gt.0)
     .write(46,*)j,ms1,ru2
c     reset mass flow in case flow would choke with zero solidity
c       (if unchoked, this does nothing...)
      ru = ru2

c     determine direction
      if ( ru .le. 0 ) then
c       flow leaving zone, know upstream conditons, find downstream

c       find conditions at screen (2)

c       ensure subsonic flow entering screen
c         take normal shock loss if m>1.
        if ( ms1 .gt. 1. ) then
c         supersonic flow
          pt1 = pt*( (gamma+1.)*ms1/(gami*ms1+2) )**(gamma/gami)
     .            *( (gamma+1.)/( 2.*gamma*ms1 - gami ) )**(1./gami)
          ms2 = ( gami*ms1 + 2. ) / ( 2.*gamma*ms1 - gami )
        else
          pt1 = pt
          ms2 = ms1
        end if

c       compute Mach at screen(2), and reset mass flow if choked
c         (limits: s=0;unchanged     s=1;ms2=choke,ru2=0)
        call corrfl(ru,tt,pt1,gamma,r,s,ms2,ru2,j)

c       Compute pressure at screen(2), using isentropic relations
c         (limit: s=0;p2=p1     s=1;p2=choke)
        fm = 1. + .5*gami*ms2
        p2 = pt1*fm**(-gamma/gami)

c       iterate for exit Mach given screen conditions

        phi = phi0 + .185*(s**.25)*(-1. + fm**fg1 )
c         (limit: s=0;phi=1.     s=1;phi=choke~.8)
        phi = min(phi,1.)
cvn     cf  = ru2*ru2*r*tt/( gamma*p2*p2*( 1.+gamma*ms2*phi*(1.-s) )**2)
        cf  = ru2*ru2*tt/( gamma*p2*p2*( 1.+gamma*ms2*phi*(1.-s) )**2)
c         (limit: s=0;cf=        s=1;cf=x*1.e-10)

c       begin on subsonic branch

        ms3 = .25
        do n = 1,50
          fm = 1. + .5*gami*ms3
          ff = -cf + ms3*fm/((1.+gamma*ms3)**2)
          fp = .5*gami*ms3/((1.+gamma*ms3)**2)
     .        - 2.*gamma*fm*ms3 /((1.+gamma*ms3)**3)
     .        + fm/((1.+gamma*ms3)**2)
cvn       co = max(ff/fp,-.5*ms3)
cvn  .        + fm/((1.+gamma*ms3)**2)
          co = max(ff/fp,-.5*ms3)
          co = min(co,.5*ms3)
          ms4 = ms3 - co
          if ( abs((ms4-ms3)/ms4) .lt. 1.e-6 ) go to 100
          if ( abs((ms4-ms3)    ) .lt. 1.e-10) go to 100
          ms3 = ms4
        end do
        iflg = 1
        if (iflg.eq.1) then
           write(11,*)'Stopping in bcporous (screen:up):'
           write(11,*)'could not converge to Mach number',m,s,
     .                                       ms3,ms4
           call ERREXIT (nodes)
        end if

100     ms3 = ms4
c         (limit; s=0;ms3=ms2=ms1=ms       s=1;ms3=x.e-10)
        fm = 1. + .5*gami*ms3
        p3 = abs(ru2)*sqrt(r*tt/(gamma*ms3*fm))
        if(idebug.gt.0)
     .  write(47,*)ru,p3,ms3
        pt2 = p3*fm**(gamma/gami)
c         (limit: s=0;pt2=pt1              s=1;pt2=p3~p2~choke)

c-------------------------------------------------------------
      else

c       flow entering zone, conditions downstream known
        ru2 = ru
        ms3 = ms1

c       compute static pressure downstream
        fm = 1. + .5*gami*ms3
        p3 = pt*fm**(-gamma/gami)

        if(idebug.gt.0)
     .  write(47,*)ru,p3,ms3

c       iterate for conditions at screen (2)
        fg1 = gamma/gami
cvn     cf  = ru*ru*r*tt/(gamma*p3*p3*(1.+gamma*ms3)**2)
        cf  = ru*ru*tt/(gamma*p3*p3*(1.+gamma*ms3)**2)
        ms2 = .25
        do n = 1,50
          fm = 1. + .5*gami*ms2
          phi = phi0 + .185*(s**.25)*(-1. + fm**fg1 )
          phi = min(phi,1.)
          dphidx = .185*(s**.25)*gamma*.5*fm**(1./gami)
          ff = -cf + phi*phi*(1.-s)*(1.-s)*ms2*fm
     .                  /( (1.+gamma*ms2*phi*(1.-s))**2 )
          fp = 2.*phi*dphidx*(1.-s)*(1.-s)*ms2*fm
     .                  /( (1.+gamma*ms2*phi*(1.-s))**2 )
     .        +phi*phi*(1.-s)*(1.-s)*fm
     .                  /( (1.+gamma*ms2*phi*(1.-s))**2 )
     .        +.5*gami*phi*phi*(1.-s)*(1.-s)*ms2
     .                  /( (1.+gamma*ms2*phi*(1.-s))**2 )
     .        -2.*phi*phi*(1.-s)*(1.-s)*ms2*fm
     .                  *(gamma*phi*(1.-s)+gamma*ms2*dphidx*(1.-s))
     .                  /( (1.+gamma*ms2*phi*(1.-s))**3 )
          co = max(ff/fp,-.5*ms2)
          co = min(co,.5*ms2)
          ms1 = ms2 - co

          if(idebug.gt.0)
     .    write(48,301)n,ff,fp,ff/fp,co,ms2,abs((ms1-ms2)/ms1)
301       format(i3,6(1x,e11.4))

          if ( abs((ms1-ms2)/ms1) .lt. 1.e-5 ) go to 207
          if ( abs((ms1-ms2)    ) .lt. 1.e-10) go to 207
          ms2 = ms1
        end do

        write(11,*) 'screen:down:could not converge to mach number',m,s
        write(11,*) 'ru,tt,pt,m,m3,ms1,2)',ru,tt,pt,m,sqrt(ms3),ms1,ms2

c--------------------------------------------------------------------------
c       coment out next section. Makes no sense to set ms1 to 1.4 when
c       choked. Arbitrarily reset to 1 
c----------------------------------------------------AFT012199-------------
c       set mach consistent with choked flow
c       if ( s .lt. .14 ) then
c         ms1 = (1. + 2.85*s)**2
c       else
c         ms1 = 1.4*1.4
c       end if
        ms1=1.

207     ms2 = ms1
        fm = 1. + .5*gami*ms2
        phi = phi0 + .185*(s**.25)*(-1. + fm**fg1 )
        phi = min(phi,1.)
        p2 = p3*(1.+gamma*ms3)/(1.+gamma*ms2*phi*(1.-s) )
        pt2 = p2*fm**(gamma/gami)

c       reset mass flow in case flow choked
        ru2 = sqrt( gamma/(r*tt) * ms2 * fm )
     .       * p2 * phi*(1.-s)

      end if

      return
      end


      subroutine corrfl (ru,tt,pt,gamma,r,s,ms,ru2,j)
c************************************************************************
c     Purpose:  use newton iteration to compute mach from mass flow using
c     		the corrected flow function.
c     		returns negative if flow chokes (specified too much 
c		corrected flow)
c     Provided by S. Bauer (from R. Bush)
c************************************************************************
c     ms is mach squared

      real ms,ms1,ms2

      idebug = 0
      gami = gamma - 1.
      rgami = 1./gami
      fg1 = - (1. + gamma) * rgami
      fg2 = - 2.*gamma *rgami
      fg3 = -(gamma) * rgami
      phi0 = .041366/(1.098216 + s - 1.) + 0.573265 + .005786*(1.-s)

c     compute choked flow for this pt,tt,s
c       note: after s=.15 phi,s make corr flow function unbounded,
c             can solve for any mass flow
c             curve fit is from looking at peak of flow function vs x
c             for various s's

      if ( s .eq. 0. ) then
c       bypass check for choked flow
        ru2 = ru
      else
c       check for choked flow and limit mass flow if needed
c       if ( s .lt. .14 ) then
c         ms1 = (1. + 2.85*s)**2
c       else
c         ms1 = 1.4*1.4
c       end if
        ms1=1.
        fm = 1. + .5*gami*ms1
        phi = phi0 + .185*(s**.25)*(-1. + fm**(-.5*fg2) )
        phi = min(phi,1.)
cvn     ru2 = sqrt(ms1*gamma/(r*tt) ) * pt * phi*(1.-s)* fm**(.5*fg1)
        ru2 = sqrt(ms1*gamma/tt ) * pt * phi*(1.-s)* fm**(.5*fg1)
c       if no mass flow, set ru2 to some small number to avoid singularity
        ru2 = max(ru2,1.e-8)
        if ( ru .lt. 0. ) ru2 = -ru2
        if ( abs(ru2) .le. abs(ru) ) then
c         limit mass flow
          ms = ms1
          return
        else
c         reset same mass flow
          ru2 = ru
        end if
      end if

      corrs = ru*ru*tt/(pt*pt)
      fs = 1.
      ms1 = ms

      do n = 1,100
        fm = 1. + .5*gami*ms1
        if ( s .gt. 0. )  then
          phi = phi0 + .185*(s**.25)*(-1. + fm**(-.5*fg2) )
          phi = min(phi,1.)
          fs = phi*phi*(1.-s)*(1.-s)
        end if
cvn     ff = corrs - fs * gamma/r * ms1 * fm**fg1
        ff = corrs - fs * gamma   * ms1 * fm**fg1
        fp = -gamma/r *fs*( fm**fg1 - .5*(gamma + 1.)*ms1 *fm**fg2 )
        if ( s .gt. 0. ) then
cvn       fp = fp - gamma/r *(1.-s)*(1.-s)*.185*(s**.25) *2.*phi
          fp = fp - gamma   *(1.-s)*(1.-s)*.185*(s**.25) *2.*phi
     .                      *.5*gamma*ms1*fm**fg3
        end if
        co = ff/fp
        co = max(co,-.5*ms1)
        co = min(co,.5*ms1)
        ms2 = ms1 - co

        if(idebug.gt.0)
     .  write(49,301)j,n,ff,fp,ff/fp,co,ms2,abs((ms2-ms1)/ms2)
301     format(2i3,6(1x,e11.4))

        if ( abs((ms2-ms1)/ms2) .lt. 1.e-5 ) go to 100
        ms1 = ms2
      end do

c     set mach consistent with choked flow
c     if ( s .lt. .14 ) then
c       ms2 = (1. + 2.85*s)**2
c     else
c       ms2 = 1.4*1.4
c     end if
      ms2=1. 

100   ms = ms2

      return
      end
