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
      subroutine start( mstrt,imn,jmn,kmn,im,jm,km,
     .                  imp2,jmp2,kmp2,w,p,eomu,turv1,turv2,smin,
     .                  wold,tv1old,tv2old,ntorder,iunsteady,jturb,
     .                  nnorder )
c
c   reading from unit "irstfil" the flow quantities
c   from previous run needed for restarting the solution
c

      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension  w(imp2,jmp2,kmp2,5), p(imp2,jmp2,kmp2),
     .        eomu(imp2,jmp2,kmp2),turv1(imp2,jmp2,kmp2),
     .        smin(imp2,jmp2,kmp2),turv2(imp2,jmp2,kmp2)
c
      dimension  wold(imp2,jmp2,kmp2,5,ntorder),
     .           tv1old(imp2,jmp2,kmp2,ntorder),
     .           tv2old(imp2,jmp2,kmp2,ntorder)
c
          read (irstfil)
     .    (((w(i,j,k,1), i=1,imp2),j=1,jmp2),k=1,kmp2),
     .    (((w(i,j,k,2), i=1,imp2),j=1,jmp2),k=1,kmp2),
     .    (((w(i,j,k,3), i=1,imp2),j=1,jmp2),k=1,kmp2),
     .    (((w(i,j,k,4), i=1,imp2),j=1,jmp2),k=1,kmp2),
     .    (((w(i,j,k,5), i=1,imp2),j=1,jmp2),k=1,kmp2),
     .    (((p(i,j,k  ), i=1,imp2),j=1,jmp2),k=1,kmp2),
     .    (((eomu(i,j,k),i=1,imp2),j=1,jmp2),k=1,kmp2)
c
cvn       if(iturb.eq.2) read(irstfil)
          if(jturb.eq.2) read(irstfil)
     .   (((smin (i,j,k),i=1,imp2),j=1,jmp2),k=1,kmp2),
     .   (((turv1(i,j,k),i=1,imp2),j=1,jmp2),k=1,kmp2)
c
cvn       if(iturb.eq.3) read(irstfil)
          if(jturb.eq.3) read(irstfil)
     .   (((smin (i,j,k),i=1,imp2),j=1,jmp2),k=1,kmp2),
     .   (((turv1(i,j,k),i=1,imp2),j=1,jmp2),k=1,kmp2),
     .   (((turv2(i,j,k),i=1,imp2),j=1,jmp2),k=1,kmp2)
c
cvn
        if(iturb.ne.jturb) 
     .  call init_turb(imn,jmn,kmn,imp2,jmp2,kmp2,w,p,eomu,
     .                 smin,turv1,turv2,ibloc )
cvn
      if (iunsteady.eq.0.or.ntorder.lt.1) return
c
c     code for accelerated time-accurate calculations
c
      if (mstrt .eq. 2) then
c
c     set wold values equal to values from `steady' restart
c
      do iorder=1,ntorder
        do k=1,kmp2
          do j=1,jmp2
            do i=1,imp2
              do nv=1,5
                wold(i,j,k,nv,iorder) =     w(i,j,k,nv)
              enddo
            enddo
c
            if (iturb.ge.2) then
              do i=1,imp2
                tv1old(i,j,k,iorder) = turv1(i,j,k)
              enddo
            endif

            if (iturb.eq.3) then
              do i=1,imp2
                tv2old(i,j,k,iorder) = turv2(i,j,k)
              enddo
            endif

          enddo
        enddo
      enddo
c
c
      else
c
c    read in wold values from `unsteady' restart file
c
        maxorder = ntorder
        if (nnorder.lt.ntorder) maxorder = nnorder
        do iorder=1,maxorder
          read (irstfil)
     .    (((wold(i,j,k,1,iorder), i=1,imp2),j=1,jmp2),k=1,kmp2),
     .    (((wold(i,j,k,2,iorder), i=1,imp2),j=1,jmp2),k=1,kmp2),
     .    (((wold(i,j,k,3,iorder), i=1,imp2),j=1,jmp2),k=1,kmp2),
     .    (((wold(i,j,k,4,iorder), i=1,imp2),j=1,jmp2),k=1,kmp2),
     .    (((wold(i,j,k,5,iorder), i=1,imp2),j=1,jmp2),k=1,kmp2)
          if(iturb.ge.2) 
     .      read (irstfil)
     .      (((tv1old(i,j,k,iorder), i=1,imp2),j=1,jmp2),k=1,kmp2)
          if(iturb.eq.3)
     .      read (irstfil)
     .      (((tv2old(i,j,k,iorder), i=1,imp2),j=1,jmp2),k=1,kmp2)
        enddo
       
        if (nnorder.lt.ntorder) then
c
          do iorder=nnorder+1,ntorder
            do k=1,kmp2
              do j=1,jmp2
                do i=1,imp2
                  do nv=1,5
                    wold(i,j,k,nv,iorder) =  wold(i,j,k,nv,iorder-1)
                  enddo
                enddo
c
                if (iturb.ge.2) then
                  do i=1,imp2
                    tv1old(i,j,k,iorder) = tv1old(i,j,k,iorder-1)
                  enddo
                endif

                if (iturb.eq.3) then
                  do i=1,imp2
                    tv2old(i,j,k,iorder) = tv2old(i,j,k,iorder-1)
                  enddo
                endif

              enddo
            enddo
          enddo
        endif
c

      endif
     
     
      return
      end
cvn
      subroutine  init_turb(imn,jmn,kmn,imp2,jmp2,kmp2,w,p,eomu,
     .                      smin,turv1,turv2,ibloc )
c=======================================================================
c 
c             initial flow field
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
     .            turv1(imp2,jmp2,kmp2),turv2(imp2,jmp2,kmp2)
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
      flow_dir = 1.0
      if(ibloc.ge.16) flow_dir = -1.0
      do 20 k=1,kmp2
      do 20 j=1,jmp2
      do 20 i=1,imp2
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
c
      return
      end
