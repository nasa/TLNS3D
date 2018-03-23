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
      subroutine  blkstp(imn, jmn, kmn, imx, jmx, kmx,
     .                   imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .                   vol,dti,dtj,dtk,dtvi,dtvj,dtvk,
     .                   dtl,ratioij,ratiojk,eprs,
#ifdef  RKPHY
     .       ntorder, dtphy ,rkphyd )
# else
     .       ntorder, dtphy )
# endif

c=======================================================================
c
c     object :  compute local time step and eigen-values for
c               dissipation and resudual smoothing routines
c
c     initial coding : by M.D.Sanetrik of AS&M (Mar. 1992)
c     modifications  : by V.N.Vatsa of NASA Langley (June, 1993)
c                     (to improve efficiency for zeta=0.5
c                      following Eli Turkel's suggestions)
c     latest mods    : by v.n.vatsa for preconditioning (July 95)
c
c
      dimension dtl(imp2,jmp2,kmp2),    dti(imp2,jmp2,kmp2),
     .          dtj(imp2,jmp2,kmp2),    dtk(imp2,jmp2,kmp2),
     .         dtvi(imp2,jmp2,kmp2),   dtvj(imp2,jmp2,kmp2),
     .         dtvk(imp2,jmp2,kmp2),
     .      ratioij(imp2,jmp2,kmp2),ratiojk(imp2,jmp2,kmp2),
     .         eprs(imp2,jmp2,kmp2,3),  vol(imp2,jmp2,kmp2)
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c-------   preconditioning related quantities      -----------------
c
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubfil2,iprec
c
      common/rk/    cfl,c(6),qfil(6),beta(6),vt,hm,mstage
c
      common/rkdis/ vis0,vis2,vis4,zeta,vepsn,vepsl,enteps,icau
c
      common/rkrsm/ smoopi,smoopj,smoopk,smoopic,smoopjc,smoopkc
c
cvn   common/smooie/ omega,ier
      common/timek/  ct(9)
      common/vis   / diffac,rsvis,omega
      common/unscom1/ anin,aninb,bninb,betvis,betvisv
      common/unsprec/ npcons,nexp,iprecg,ijac,icord,inav
c
c     cfls   = (real(mstage)-1.)*.937
      cflrat = cfl/cfls
c     cflrat = max(cflrat,1.6)
      difac  = diffac
c
c          for 2d implementation i2dfl=1
      ai2  = i2dfl
      ai2m = 1. - ai2
c
c         for 2d a12 = 1  ai2m = 0
c         for sd a12 = 0  ai2m = 1
c
c         omega = 1     original Martinelli
c         omega = 2     not presenetly used
c         omega = 3     Swanson sum in 3D
      if (omega .eq. 1) then
          om1 = 1.
c         om2 = 0.
          om3 = 0.
      endif
c     if (omega .eq. 2) then
c         om1 = 0.
c         om2 = 1.
c         om3 = 0.
c     endif
      if (omega .eq. 3) then
          om1 = 0.
c         om2 = 0.
          om3 = 1.
      endif
   
c     omega2 = .03125*(1. + ai2)
c     omega2 = .03125
      omega3 = .02*(12.5*ai2 + ai2m)
c
        if(abs(zeta-0.5).le.0.0001) then
c
c     utilize sqrt functions when zeta=0.5
c
      do 20 k=kmn,kmx+1
      do 20 j=jmn,jmx+1
      do 20 i=imn,imx+1
c
c--------- step limit for complete difference equation
c
c     dtli     =  1./dti(i,j,k)
c     dtlj     =  1./dtj(i,j,k)
c     dtlk     =  1./dtk(i,j,k)
c
c     dtl(i,j,k)  = vol(i,j,k)*dtli*dtlj*dtlk/
c    .             (dtli*dtlj  +dtlj*dtlk  +dtlk*dtli)
c
c     dtl(i,j,k)  = ai2m*dtl(i,j,k)
c    .            + ai2*vol(i,j,k)*dtli*dtlj/(dtli+dtlj)
c
      dtli     =  dti(i,j,k) + difac*dtvi(i,j,k)
      dtlj     =  dtj(i,j,k) + difac*dtvj(i,j,k)
      dtlk     =  dtk(i,j,k) + difac*dtvk(i,j,k)
c               anin is factor for influence of
c         physical time derivative on artificial time step
      eigt     =  anin*ct(ntorder+1)*vol(i,j,k) / dtphy
      dtl(i,j,k)  = vol(i,j,k)/(eigt + dtli + dtlj + ai2m*dtlk)
#ifdef  RKPHY
      dtl(i,j,k)  = dtl(i,j,k) /rkphyd
# endif
c     dtld     =  dti(i,j,k) + dtj(i,j,k) + ai2m*dtk(i,j,k)
c
      reigij   = sqrt(dtli/dtlj)
      reigji   =  1./reigij
      reigjk   = sqrt(dtlj/dtlk)
      reigik   = reigij*reigjk
c
c     2-d implementation
      reigjk   = ai2m*reigjk + ai2
      reigik   = ai2m*reigik + ai2
      dtlk     = dtlk*ai2m
c
      reigkj   =  1./reigjk
      reigki   =  1./reigik
      ratioij(i,j,k) = reigij
      ratiojk(i,j,k) = reigjk
c
c         constant coefficients for point Jacobi
c
      scali    =  1. + reigji + reigki*ai2m
      scalj    =  1. + reigij + reigkj*ai2m
      scalk    =  1. + reigik + reigjk
c     scali    =  1. + max(reigji , reigki*ai2m)
c     scalj    =  1. + max(reigij , reigkj*ai2m)
c     scalk    =  1. + ai2m*max(reigik , reigjk)
c
c         dtratx1 = Martinelli scaling
c         dtratx3 = Swanson scaling   sum
c
      dtratx1  = dtli*scali/(dtli+dtlj+dtlk)
c     dtm      = ai2m*min(dtlj,dtlk) + ai2*dtlj
c     radval   = 1./(dtli+omega2*dtm)
c     dtratx2  = dtli*radval
      dtratx3  = dtli/(dtli+omega3*(dtlj+dtlk))
c     dtratx3  = 1.
      dtratx   = om1*dtratx1 +               om3*dtratx3
      ex       = .25*((cflrat*dtratx)**2 -1.)
      ex       = min(smoopj,ex)
     
      exv      = rsvis*dtvi(i,j,k)/dti(i,j,k)
c     exv      = rsvis*dtvi(i,j,k)/dtld
      ex       = max(ex,exv)
      eprs(i,j,k,1) = max(smoopi,ex)
c
      dtraty1  = dtlj*scalj/(dtli+dtlj+dtlk)
c     dtm      = ai2m*min(dtli,dtlk) + ai2*dtli
c     radval   = 1./(dtlj+omega2*dtm)
c     dtraty2  = dtlj*radval
      dtraty3  = dtlj/(dtlj+omega3*(dtli+dtlk))
      dtraty   = om1*dtraty1 +               om3*dtraty3
      ey       = .25*((cflrat*dtraty)**2 -1.)
      ey       = min(smoopj,ey)
      eyv      = rsvis*dtvj(i,j,k)/dtj(i,j,k)
c     eyv      = rsvis*dtvj(i,j,k)/dtld
      ey       = max(ey,eyv)
      eprs(i,j,k,2) = max(smoopi,ey)
c
      dtratz1  = dtlk*scalk/(dtli+dtlj+dtlk)
c     radval   = 1./(dtlk+omega2*min(dtli,dtlj))
c     dtratz2  = dtlk*radval
      dtratz3  = dtlk/(dtlk+omega3*(dtli+dtlj))
      dtratz   = om1*dtratz1 +               om3*dtratz3
      ez       = .25*((cflrat*dtratz)**2 -1.)
      ez       = min(smoopj,ez)
      ezv      = rsvis*dtvk(i,j,k)/dtk(i,j,k)
c     ezv      = rsvis*dtvk(i,j,k)/dtld
      ez       = max(ez,ezv)
      eprs(i,j,k,3) = ai2m*max(smoopi,ez)
c
      sum      = eprs(i,j,k,1)+eprs(i,j,k,2)+eprs(i,j,k,3)
      sum1     = min(smoopk,sum)
      scalr    = sum1/sum
      eprs(i,j,k,1) = scalr*eprs(i,j,k,1)
      eprs(i,j,k,2) = scalr*eprs(i,j,k,2)
      eprs(i,j,k,3) = scalr*eprs(i,j,k,3)
c
c     dtld     = dti(i,j,k) + dtj(i,j,k) + ai2m*dtk(i,j,k)
c     exv      =  rsvis*dtvi(i,j,k)/dti(i,j,k)
c     exv      =  rsvis*dtvi(i,j,k)/dtld
c     ex       =  max(ex,exv)
c     eprs(i,j,k,1) = max(smoopi,ex)
c
c     eyv      =  rsvis*dtvj(i,j,k)/dtj(i,j,k)
c     eyv      =  rsvis*dtvj(i,j,k)/dtld
c     ey       =  max(ey,eyv)
c     eprs(i,j,k,2) = max(smoopi,ey)
c     ezv      =  rsvis*dtvk(i,j,k)/dtk(i,j,k)
c     ezv      =  rsvis*dtvk(i,j,k)/dtld
c     ez       =  max(ez,ezv)
c     eprs(i,j,k,3) = ai2m*max(smoopi,ez)
c
   20 continue
c 
         else
c
c     zeta .ne. 0.5 case
c
      do 30 k=kmn,kmx+1
      do 30 j=jmn,jmx+1
      do 30 i=imn,imx+1
c
c--------- step limit for complete difference equation
c
c     dtli     =  1./dti(i,j,k)
c     dtlj     =  1./dtj(i,j,k)
c     dtlk     =  1./dtk(i,j,k)
c
c     dtl(i,j,k)  = vol(i,j,k)*dtli*dtlj*dtlk/
c    .             (dtli*dtlj  +dtlj*dtlk  +dtlk*dtli)
c
c     dtl(i,j,k)  = ai2m*dtl(i,j,k)
c    .            + ai2*vol(i,j,k)*dtli*dtlj/(dtli+dtlj)
c
      dtli     =  dti(i,j,k) + difac*dtvi(i,j,k)
      dtlj     =  dtj(i,j,k) + difac*dtvj(i,j,k)
      dtlk     =  dtk(i,j,k) + difac*dtvk(i,j,k)
c               anin is factor for influence of
c         physical time derivative on artificial time step
      eigt     =  anin*ct(ntorder+1)*vol(i,j,k) / dtphy
      dtl(i,j,k)  = vol(i,j,k)/(eigt + dtli + dtlj + ai2m*dtlk)
#ifdef  RKPHY
      dtl(i,j,k)  = dtl(i,j,k) /rkphyd
# endif
c     dtld     =  dti(i,j,k) + dtj(i,j,k) + ai2m*dtk(i,j,k)
c
      reigij   = (dtli/dtlj)**zeta
      reigji   =  1./reigij
      reigjk   = (dtlj/dtlk)**zeta
      reigik   = reigij*reigjk
c
c     2-d implementation
      reigjk   = ai2m*reigjk + ai2
      reigik   = ai2m*reigik + ai2
      dtlk     = dtlk*ai2m
c
      reigkj   =  1./reigjk
      reigki   =  1./reigik
      ratioij(i,j,k) = reigij
      ratiojk(i,j,k) = reigjk
c
c         constant coefficients for point Jacobi
c     if (ijac.gt.0) go to 30
c
      scali    =  1. + reigji + reigki*ai2m
      scalj    =  1. + reigij + reigkj*ai2m
      scalk    =  1. + reigik + reigjk
c     scali    =  1. + max(reigji , reigki*ai2m)
c     scalj    =  1. + max(reigij , reigkj*ai2m)
c     scalk    =  1. + ai2m*max(reigik , reigjk)
c
c         dtratx1 = Martinelli scaling
c         dtratx3 = Swanson scaling   sum
c
      dtratx1  = dtli*scali/(dtli+dtlj+dtlk)
c     dtm      = ai2m*min(dtlj,dtlk) + ai2*dtlj
c     radval   = 1./(dtli+omega2*dtm)
c     dtratx2  = dtli*radval
      dtratx3  = dtli/(dtli+omega3*(dtlj+dtlk))
c     dtratx3  = 1.
      dtratx   = om1*dtratx1 +               om3*dtratx3
      ex       = .25*((cflrat*dtratx)**2 -1.)
      ex       = min(smoopj,ex)
     
      exv      = rsvis*dtvi(i,j,k)/dti(i,j,k)
c     exv      = rsvis*dtvi(i,j,k)/dtld
      ex       = max(ex,exv)
      eprs(i,j,k,1) = max(smoopi,ex)
c
      dtraty1  = dtlj*scalj/(dtli+dtlj+dtlk)
c     dtm      = ai2m*min(dtli,dtlk) + ai2*dtli
c     radval   = 1./(dtlj+omega2*dtm)
c     dtraty2  = dtlj*radval
      dtraty3  = dtlj/(dtlj+omega3*(dtli+dtlk))
      dtraty   = om1*dtraty1 +               om3*dtraty3
      ey       = .25*((cflrat*dtraty)**2 -1.)
      ey       = min(smoopj,ey)
      eyv      = rsvis*dtvj(i,j,k)/dtj(i,j,k)
c     eyv      = rsvis*dtvj(i,j,k)/dtld
      ey       = max(ey,eyv)
      eprs(i,j,k,2) = max(smoopi,ey)
c
      dtratz1  = dtlk*scalk/(dtli+dtlj+dtlk)
c     radval   = 1./(dtlk+omega2*min(dtli,dtlj))
c     dtratz2  = dtlk*radval
      dtratz3  = dtlk/(dtlk+omega3*(dtli+dtlj))
      dtratz   = om1*dtratz1 +               om3*dtratz3
      ez       = .25*((cflrat*dtratz)**2 -1.)
      ez       = min(smoopj,ez)
      ezv      = rsvis*dtvk(i,j,k)/dtk(i,j,k)
c     ezv      = rsvis*dtvk(i,j,k)/dtld
      ez       = max(ez,ezv)
      eprs(i,j,k,3) = ai2m*max(smoopi,ez)
c
      sum      = eprs(i,j,k,1)+eprs(i,j,k,2)+eprs(i,j,k,3)
      sum1     = min(smoopk,sum)
      scalr    = sum1/sum
      eprs(i,j,k,1) = scalr*eprs(i,j,k,1)
      eprs(i,j,k,2) = scalr*eprs(i,j,k,2)
      eprs(i,j,k,3) = scalr*eprs(i,j,k,3)
c
c     dtld     = dti(i,j,k) + dtj(i,j,k) + ai2m*dtk(i,j,k)
c     exv      =  rsvis*dtvi(i,j,k)/dti(i,j,k)
c     exv      =  rsvis*dtvi(i,j,k)/dtld
c     ex       =  max(ex,exv)
c     eprs(i,j,k,1) = max(smoopi,ex)
c
c     eyv      =  rsvis*dtvj(i,j,k)/dtj(i,j,k)
c     eyv      =  rsvis*dtvj(i,j,k)/dtld
c     ey       =  max(ey,eyv)
c     eprs(i,j,k,2) = max(smoopi,ey)
c     ezv      =  rsvis*dtvk(i,j,k)/dtk(i,j,k)
c     ezv      =  rsvis*dtvk(i,j,k)/dtld
c     ez       =  max(ez,ezv)
c     eprs(i,j,k,3) = ai2m*max(smoopi,ez)
   30 continue
c
         endif

      return
      end
