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
      subroutine  init(imn,jmn,kmn,imp2,jmp2,kmp2,w,p,eomu,
     .                     smin,x,turv1,turv2,ibloc ,wold,
     .                     tv1old,tv2old,ntorder,iunsteady,igrid ,
     .                     dtphy )
c=======================================================================
c 
c              initial flow field
c
c

      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c 
c
      common/trbfs/ anutinf,zkinf,ominf
c
      dimension   w(imp2,jmp2,kmp2,5),p(imp2,jmp2,kmp2),
     .            eomu(imp2,jmp2,kmp2), smin(imp2,jmp2,kmp2),
     .            turv1(imp2,jmp2,kmp2),turv2(imp2,jmp2,kmp2),
     .            x(imp2-1,jmp2-1,kmp2-1,3)
c
      dimension   wold(imp2,jmp2,kmp2,5,ntorder),
     .            tv1old(imp2,jmp2,kmp2,ntorder),
     .            tv2old(imp2,jmp2,kmp2,ntorder)
c
c   ***data of "typical" k,omega profiles for initial guess***
      parameter (mxtrbin = 32)
      dimension zks(mxtrbin),oms(mxtrbin)
      data zks /.182e-8,.751e-6,.998e-5,.158e-3,.116e-2,
     +          .232e-2,.419e-2,.694e-2,.108e-1,.151e-1,
     +          .139e-1,.157e-2,.125e-4,.420e-7,.900e-8,
     +          .900e-8,.900e-8,.900e-8,.900e-8,.900e-8,
     +          .900e-8,.900e-8,.900e-8,.900e-8,.900e-8,
     +          .900e-8,.900e-8,.900e-8,.900e-8,.900e-8,
     +          .900e-8,.900e-8 /
      data oms /.104e+1,.387e-1,.797e-2,.152e-2,.291e-3,
     +          .648e-4,.236e-4,.125e-4,.902e-5,.823e-5,
     +          .711e-5,.144e-5,.100e-5,.100e-5,.100e-5,
     +          .100e-5,.100e-5,.100e-5,.100e-5,.100e-5,
     +          .100e-5,.100e-5,.100e-5,.100e-5,.100e-5,
     +          .100e-5,.100e-5,.100e-5,.100e-5,.100e-5,
     +          .100e-5,.100e-5 /
c   ********
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
c
      do 20 k=1,kmp2
      do 20 j=1,jmp2
      do 20 i=1,imp2
      w(i,j,k,1)  = rho0
      w(i,j,k,2)  = rho0*u0 
      w(i,j,k,3)  = rho0*v0 
      w(i,j,k,4)  = rho0*w0 
      w(i,j,k,5)  = -p0 + h0*rho0
      p(i,j,k)    = p0
      vtots        =    w(i,j,k,2)**2 + w(i,j,k,3)**2 +
     .                  w(i,j,k,4)**2
      w(i,j,k,5)   =    p(i,j,k) *dgm1 + .5*vtots / w(i,j,k,1)
      eomu(i,j,k) = evinf
   20 continue
c
      if(iturb.eq.2) then
c     spalart's model
c
        do 30 k=1,kmp2
        do 30 j=1,jmp2
        do 30 i=1,imp2
        turv1(i,j,k)  =  anutinf
 30     continue
      endif
c
      if(iturb.eq.3) then
c     sst (k-omega) model
c
cnew (feb. 99)
        const1 = 45.8
        const2 = 1.68
        smax   = const2/(2.*const1)
        tmax   =-const1*smax**2 + const2*smax
        cmu    = 1.
        do 40 k=1,kmp2
        do 40 i=1,imp2
c
        turv1(i,1,k)  =  zks(1)
        turv1(i,2,k)  =  zks(1)
        turv2(i,1,k)  =  oms(1)
        turv2(i,2,k)  =  oms(1)
        turv1(i,jmp2-1,k)  =  zks(mxtrbin)
        turv1(i,jmp2  ,k)  =  zks(mxtrbin)
        turv2(i,jmp2-1,k)  =  oms(mxtrbin)
        turv2(i,jmp2  ,k)  =  oms(mxtrbin)
c
c       do 40 j=3,jmp2-2
c       dist   = real(j-2)/real(jmp2-3)
c       yjj    = 1. + dist*real(mxtrbin-1)
c       jj     = int(yjj)
c       fraction = yjj - real(jj)
c       turv1(i,j,k) = fraction*zks(jj+1) + (1.-fraction)*zks(jj)
c       turv2(i,j,k) = fraction*oms(jj+1) + (1.-fraction)*oms(jj)
c40     continue
        do 40 j=3,jmp2-2
        ak1    = 0.9e-08
        ak2    = -45.8*smin(i,j,k)*smin(i,j,k) + 1.68*smin(i,j,k)
        turv1(i,j,k) = max(ak1,ak2)
        aw1    = -12444.0*smin(i,j,k) + 0.54
c       aw2    = -0.00007*smin(i,j,k) + 0.00001
        v3d    = turv1(i,j,k)*100/tmax
        aw2    =  cmu*turv1(i,j,k)/v3d
        aw     =  max(aw1,aw2)
c       aw3    =  0.1e-05
c       turv2(i,j,k) = max(aw,aw3)
        turv2(i,j,k) = max(aw1,aw2)
 40     continue
      endif
c
      if(igrid.gt.1.or.iunsteady.eq.0) return
c
      do 50 iorder = 1, ntorder
      do 50 k=1,kmp2
      do 50 j=1,jmp2
      do 50 i=1,imp2
      wold(i,j,k,1,iorder)  =     w(i,j,k,1)
      wold(i,j,k,2,iorder)  =     w(i,j,k,2)
      wold(i,j,k,3,iorder)  =     w(i,j,k,3)
      wold(i,j,k,4,iorder)  =     w(i,j,k,4)
      wold(i,j,k,5,iorder)  =     w(i,j,k,5)
      if(iturb.ge.2)
     .tv1old(i,j,k,iorder)  = turv1(i,j,k)
      if(iturb.eq.3)
     .tv2old(i,j,k,iorder)  = turv2(i,j,k)
   50 continue
c
      return
      end
