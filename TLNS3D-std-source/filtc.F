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
      subroutine  filtc (nstage,imn,jmn,kmn,imx,jmx,kmx,
     .                   imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .                   w,p,vol,dti,dtj,dtk,dtl,ratioij,ratiojk,
     .                   mx1dwk,fs1,fs2,fs3,fs4,fs5,phi,h,fw,
     .                   dtvi,dtvj,dtvk,fbeta2 )
c=======================================================================
c 
c           artificial dissipation on coarse grids
c           second differences with fixed coefficient
c
c     initial coding : by v.n.vatsa (1990-1991)
c     latest mods    : by v.n.vatsa for preconditioning (Feb. 1998)
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c
c-------   preconditioning related quantities      -----------------
c
      common/precon/ cfls,ulprec,ubmin2,upk,alppre,deltpre,ubcon2,iprec
c
      common/rk/    cfl,c(6),qfil(6),beta(6),vt,hm,mstage
c
      common/rkdis/ vis0,vis2,vis4,zeta,vepsn,vepsl,enteps,icau
c
      common/rkrsm/ smoopi,smoopj,smoopk,smoopic,smoopjc,smoopkc
      common/unsprec/ npcons,nexp,iprecg,ijac,icord,inav
      common/vis   / diffac,rsvis,omega
c
c
      dimension    fw(imp2,jmp2,kmp2,5),
     .              w(imp2,jmp2,kmp2,5),
     .              p(imp2,jmp2,kmp2),vol(imp2,jmp2,kmp2),
     .            dtl(imp2,jmp2,kmp2),dti(imp2,jmp2,kmp2),
     .            dtj(imp2,jmp2,kmp2),dtk(imp2,jmp2,kmp2),
     .            ratioij(imp2,jmp2,kmp2),ratiojk(imp2,jmp2,kmp2),
     .       dtvi(imp2,jmp2,kmp2),dtvj(imp2,jmp2,kmp2),
     .       dtvk(imp2,jmp2,kmp2),
     .     fbeta2(imp2,jmp2,kmp2)
c
c     space allocation for working arrays
c
      dimension   phi(mx1dwk),  h(mx1dwk),fs1(mx1dwk),fs2(mx1dwk),
     .            fs3(mx1dwk),fs4(mx1dwk),fs5(mx1dwk)
c 
c=======================================================================
c
c 
      fis0      = vis0/32.
      fil1      = 1.
      if (qfil(1).lt.0.) fil1 = 0.
      bet       = beta(nstage)
      cp        = gamma/gm1
      difac     = diffac
      if (inav.eq.1) difac = 0.
c       keeps viscous contribution in scalar viscosity
      delta     = deltpre
c
c     initialize fw during first call on a given mesh
c
      if (nstage .eq. 1) then
        do nv = 1,5
        do k=1,kmp2
        do j=1,jmp2
        do i=1,imp2
          fw(i,j,k,nv) = 0.
        enddo
        enddo
        enddo
        enddo
      endif
c
c--------- dissipation in i direction 
c
      do 99 k=kmn+1,kmx
      do 99 j=jmn+1,jmx
      do i=imn,imp1
        reigji   =  1./ratioij(i,j,k)
        reigki   =  1./(ratioij(i,j,k)*ratiojk(i,j,k))
        phi(i)   =  1. + reigji + reigki*(1.-i2dfl)
      enddo
c  
      if (iprec.le.1) then
c
c            form first differences of conservative variables
c
        do i=imn,imx
          fila    = 0.5*((dti(i+1,j,k)+difac*dtvi(i+1,j,k))*phi(i+1) 
     .                 + (dti(i  ,j,k)+difac*dtvi(i  ,j,k))*phi(i  ))
          dis     = fis0*fila
          fs1(i)  = dis*(w(i+1,j,k,1) - w(i,j,k,1)) 
          fs2(i)  = dis*(w(i+1,j,k,2) - w(i,j,k,2)) 
          fs3(i)  = dis*(w(i+1,j,k,3) - w(i,j,k,3)) 
          fs4(i)  = dis*(w(i+1,j,k,4) - w(i,j,k,4)) 
          fs5(i)  = dis*(w(i+1,j,k,5) - w(i,j,k,5)  
     .       + (1-icau)*(p(i+1,j,k  ) - p(i,j,k) ))
        enddo
c
      else
c
c            form first differences of "u,v,w,p,T" variables
c
        rhoi      = 1./w(imn,j,k,1)
        do i=imn,imx
          fila    = 0.5*((dti(i+1,j,k)+difac*dtvi(i+1,j,k))*phi(i+1) 
     .                 + (dti(i  ,j,k)+difac*dtvi(i  ,j,k))*phi(i  ))
          dis     = fis0*fila
          rhoip   = 1./w(i+1,j,k,1)
          fs1(i)  = dis*(p(i+1,j,k  )       - p(i,j,k  )     ) 
          fs2(i)  = dis*(w(i+1,j,k,2)*rhoip - w(i,j,k,2)*rhoi) 
          fs3(i)  = dis*(w(i+1,j,k,3)*rhoip - w(i,j,k,3)*rhoi) 
          fs4(i)  = dis*(w(i+1,j,k,4)*rhoip - w(i,j,k,4)*rhoi) 
          fs5(i)  = dis*(p(i+1,j,k  )*rhoip - p(i,j,k  )*rhoi)  
          rhoi    = rhoip
        enddo
      endif
c
      if (iprec.eq.1) then
c            for conservative formulation multiply by P^{-1}
c            assumes alpha = 0
      do 70 i=imn,imx
        rhoi      = 1./w(i,j,k,1)
        uc        = w(i,j,k,2)*rhoi
        vc        = w(i,j,k,3)*rhoi
        wc        = w(i,j,k,4)*rhoi
        q2        = 0.5*(uc**2 + vc**2 + wc**2)
        ts        = p(i,j,k)*rhoi
        c2        = gamma*ts
        hc        = cp*ts + q2
        bc        = fbeta2(i,j,k)
        delta     = deltpre
        z         = gm1*(q2*fs1(i) - uc*fs2(i)
     .            - vc*fs3(i) - wc*fs4(i)
     .            + fs5(i) )/c2
        dw1       = ((1. - bc)*z - delta*(z - fs1(i)))/bc
        fs1(i)    = fs1(i) + dw1
        fs2(i)    = fs2(i) + uc*dw1
        fs3(i)    = fs3(i) + vc*dw1
        fs4(i)    = fs4(i) + wc*dw1
        fs5(i)    = fs5(i) + hc*dw1
   70 continue
      endif
      if (iprec.eq.2) then
c            for conservative formulation multiply by Gamma^{-1}
      if (delta.eq.0.) then
      do 80 i=imn,imx
        rhoa    = 0.5*(w(i,j,k,1)+w(i+1,j,k,1))
        rhoi    = 1./rhoa
        rhou    = 0.5*(w(i,j,k,2) + w(i+1,j,k,2))
        rhov    = 0.5*(w(i,j,k,3) + w(i+1,j,k,3))
        rhow    = 0.5*(w(i,j,k,4) + w(i+1,j,k,4))
        ts      = 0.5*(p(i,j,k)   + p(i+1,j,k))*rhoi
        hqs     = 0.5*(rhou*rhou+rhov*rhov+rhow*rhow)*rhoi*rhoi
        c2      = gamma*ts
        dc2     = 1./c2
        hh      = c2/gm1 + hqs
c
c               evaluate alpha and beta for preconditioning
c
c       bc      = 0.5*(fbeta2(i,j,k) + fbeta2(i+1,j,k))
        bc      = max (fbeta2(i,j,k) , fbeta2(i+1,j,k))
        alp     = max(ulprec*(1.-bc),0.)
        dbeta2  = 1./bc
        theta   = (gm1 + dbeta2)*dc2
        rtq5    = rhoa*fs5(i)/ts
        aq1     = alp*fs1(i)*dbeta2*dc2
c
        dw1     = theta*fs1(i) - rtq5
        dwpt    = rhoi*(dw1 + aq1)
        dw2     = rhoa*fs2(i) + rhou*dwpt
        dw3     = rhoa*fs3(i) + rhov*dwpt
        dw4     = rhoa*fs4(i) + rhow*dwpt
        dw5     = (theta*hh-1.)*fs1(i)
     .          + rhou*fs2(i) + rhov*fs3(i) + rhow*fs4(i)
     .          + hqs*(2.*aq1 - rtq5)
c
        fs1(i) = dw1
        fs2(i) = dw2
        fs3(i) = dw3
        fs4(i) = dw4
        fs5(i) = dw5
   80 continue
      else
c             delta=1
      do 85 i=imn,imx
        rhoa    = 0.5*(w(i,j,k,1)+w(i+1,j,k,1))
        rhoi    = 1./rhoa
        uu      = 0.5*(w(i,j,k,2) + w(i+1,j,k,2))*rhoi
        vv      = 0.5*(w(i,j,k,3) + w(i+1,j,k,3))*rhoi
        ww      = 0.5*(w(i,j,k,4) + w(i+1,j,k,4))*rhoi
        ww      = 0.5*(w(i,j,k,4) + w(i+1,j,k,4))*rhoi
        ts      = 0.5*(p(i,j,k)   + p(i+1,j,k))*rhoi
        hh      = 0.5*(w(i,j,k,5) + w(i+1,j,k,5))*rhoi + ts
        c2      = gamma*ts
        us      = uu*uu + vv*vv + ww*ww
c
c               evaluate alpha and beta for preconditioning
c
c       bc      = 0.5*(fbeta2(i,j,k) + fbeta2(i+1,j,k))
        bc      = max (fbeta2(i,j,k) , fbeta2(i+1,j,k))
        beta2   = bc*c2
        alp     = max(ulprec*(1.-bc),0.)
        dbeta2  = 1./beta2
        alp1    = (1.+alp)*dbeta2*fs1(i)
c
        dw1     = dbeta2*fs1(i)
        dw2     = rhoa*fs2(i) + alp1*uu
        dw3     = rhoa*fs3(i) + alp1*vv
        dw4     = rhoa*fs4(i) + alp1*ww
        g41     = dbeta2*(hh + alp*us) -1.
        dw5     = g41*fs1(i) +
     .            rhoa*(uu*fs2(i) + vv*fs3(i) + ww*fs4(i)
     .          + cp*fs5(i))
        fs1(i) = dw1
        fs2(i) = dw2
        fs3(i) = dw3
        fs4(i) = dw4
        fs5(i) = dw5
   85 continue
      endif
      endif
c
      do 99 i=imn+1,imx
      fw(i,j,k,1)  = (1.-bet)*fw(i,j,k,1) + bet*(fs1(i-1) - fs1(i))
      fw(i,j,k,2)  = (1.-bet)*fw(i,j,k,2) + bet*(fs2(i-1) - fs2(i))
      fw(i,j,k,3)  = (1.-bet)*fw(i,j,k,3) + bet*(fs3(i-1) - fs3(i))
      fw(i,j,k,4)  = (1.-bet)*fw(i,j,k,4) + bet*(fs4(i-1) - fs4(i))
      fw(i,j,k,5)  = (1.-bet)*fw(i,j,k,5) + bet*(fs5(i-1) - fs5(i))
   99 continue
c
c--------- dissipation in j direction 
c
      do 199 k=kmn+1,kmx
      do 199 i=imn+1,imx
      do j=jmn,jmp1
        reigij   = ratioij(i,j,k)               
        reigkj   =  1./ratiojk(i,j,k)
        phi(j)   =  1. + reigij + reigkj*(1.-i2dfl)
      enddo
c
      if (iprec.le.1) then
c
c            form first differences of conservative variables
c
        do j=jmn,jmx
          fila    = 0.5*((dtj(i,j+1,k)+difac*dtvj(i,j+1,k))*phi(j+1) 
     .                 + (dtj(i  ,j,k)+difac*dtvj(i,j  ,k))*phi(j  ))
          dis     = fis0*fila
          fs1(j)  = dis*(w(i,j+1,k,1) - w(i,j,k,1)) 
          fs2(j)  = dis*(w(i,j+1,k,2) - w(i,j,k,2)) 
          fs3(j)  = dis*(w(i,j+1,k,3) - w(i,j,k,3)) 
          fs4(j)  = dis*(w(i,j+1,k,4) - w(i,j,k,4)) 
          fs5(j)  = dis*(w(i,j+1,k,5) - w(i,j,k,5)
     .       + (1-icau)*(p(i,j+1,k)   - p(i,j,k) ))
        enddo
c
      else
c
c            form first differences of "u,v,w,p,T" variables
c
        rhoi      = 1./w(i,jmn,k,1)
        do j=jmn,jmx
          fila    = 0.5*((dtj(i,j+1,k)+difac*dtvj(i,j+1,k))*phi(j+1) 
     .                 + (dtj(i  ,j,k)+difac*dtvj(i,j  ,k))*phi(j  ))
          dis     = fis0*fila
          rhoip   = 1./w(i,j+1,k,1)
c         rhoi    = 1./w(i,j  ,k,1)
          fs1(j)  = dis*(p(i,j+1,k  )       - p(i,j,k  )     ) 
          fs2(j)  = dis*(w(i,j+1,k,2)*rhoip - w(i,j,k,2)*rhoi) 
          fs3(j)  = dis*(w(i,j+1,k,3)*rhoip - w(i,j,k,3)*rhoi) 
          fs4(j)  = dis*(w(i,j+1,k,4)*rhoip - w(i,j,k,4)*rhoi) 
          fs5(j)  = dis*(p(i,j+1,k  )*rhoip - p(i,j,k  )*rhoi)  
          rhoi    = rhoip
        enddo
c
      endif
      if (iprec.eq.1) then
c            for conservative formulation multiply by P^{-1}
      do 170 j=jmn,jmx
        rhoi      = 1./w(i,j,k,1)
        uc        = w(i,j,k,2)*rhoi
        vc        = w(i,j,k,3)*rhoi
        wc        = w(i,j,k,4)*rhoi
        q2        = 0.5*(uc**2 + vc**2 + wc**2)
        ts        = p(i,j,k)*rhoi
        c2        = gamma*ts
        hc        = cp*ts + q2
        bc        = fbeta2(i,j,k)
        delta     = deltpre
        z         = gm1*(q2*fs1(j) - uc*fs2(j)
     .            - vc*fs3(j) - wc*fs4(j)
     .            + fs5(j) )/c2
        dw1       = ((1. - bc)*z - delta*(z - fs1(j)))/bc
        fs1(j)    = fs1(j) + dw1
        fs2(j)    = fs2(j) + uc*dw1
        fs3(j)    = fs3(j) + vc*dw1
        fs4(j)    = fs4(j) + wc*dw1
        fs5(j)    = fs5(j) + hc*dw1
  170 continue
      endif
      if (iprec.eq.2) then
c         for conservative formulation multiply by Gamma^{-1}
      if (delta.eq.0.) then
      do 180 j=jmn,jmx
        rhoa    = 0.5*(w(i,j,k,1)+w(i,j+1,k,1))
        rhoi    = 1./rhoa
        rhou    = 0.5*(w(i,j,k,2) + w(i,j+1,k,2))
        rhov    = 0.5*(w(i,j,k,3) + w(i,j+1,k,3))
        rhow    = 0.5*(w(i,j,k,4) + w(i,j+1,k,4))
        ts      = 0.5*(p(i,j,k)   + p(i,j+1,k))*rhoi
        hqs     = 0.5*(rhou*rhou+rhov*rhov+rhow*rhow)*rhoi*rhoi
        c2      = gamma*ts
        dc2     = 1./c2
        hh      = c2/gm1 + hqs
c
c               evaluate alpha and beta for preconditioning
c
c       bc      = 0.5*(fbeta2(i,j,k) + fbeta2(i,j+1,k))
        bc      = max (fbeta2(i,j,k) , fbeta2(i,j+1,k))
        alp     = max(ulprec*(1.-bc),0.)
        dbeta2  = 1./bc
        theta   = (gm1 + dbeta2)*dc2
        rtq5    = rhoa*fs5(j)/ts
        aq1     = alp*fs1(j)*dbeta2*dc2
c
        dw1     = theta*fs1(j) - rtq5
        dwpt    = rhoi*(dw1 + aq1)
        dw2     = rhoa*fs2(j) + rhou*dwpt
        dw3     = rhoa*fs3(j) + rhov*dwpt
        dw4     = rhoa*fs4(j) + rhow*dwpt
        dw5     = (theta*hh-1.)*fs1(j)
     .          + rhou*fs2(j) + rhov*fs3(j) + rhow*fs4(j)
     .          + hqs*(2.*aq1 - rtq5)
c
        fs1(j) = dw1
        fs2(j) = dw2
        fs3(j) = dw3
        fs4(j) = dw4
        fs5(j) = dw5
  180 continue
      else
      do 185 j=jmn,jmx
        rhoa    = 0.5*(w(i,j,k,1)+w(i,j+1,k,1))
        rhoi    = 1./rhoa
        uu      = 0.5*(w(i,j,k,2) + w(i,j+1,k,2))*rhoi
        vv      = 0.5*(w(i,j,k,3) + w(i,j+1,k,3))*rhoi
        ww      = 0.5*(w(i,j,k,4) + w(i,j+1,k,4))*rhoi
        ts      = 0.5*(p(i,j,k)   + p(i,j+1,k))*rhoi
        hh      = 0.5*(w(i,j,k,5) + w(i,j+1,k,5))*rhoi + ts
        c2      = gamma*ts
        us      = uu*uu + vv*vv + ww*ww
c
c               evaluate alpha and beta for preconditioning
c
c       bc      = 0.5*(fbeta2(i,j,k) + fbeta2(i,j+1,k))
        bc      = max (fbeta2(i,j,k) , fbeta2(i,j+1,k))
        beta2   = bc*c2
        alp     = max(ulprec*(1.-bc),0.)
        dbeta2  = 1./beta2
        alp1    = (1.+alp)*dbeta2*fs1(j)
c
        dw1     = dbeta2*fs1(j)
        dw2     = rhoa*fs2(j) + alp1*uu
        dw3     = rhoa*fs3(j) + alp1*vv
        dw4     = rhoa*fs4(j) + alp1*ww
        g41     = dbeta2*(hh + alp*us) -1.
        dw5     = g41*fs1(j)
     .          + rhoa*(uu*fs2(j) + vv*fs3(j) + ww*fs4(j)
     .          + cp*fs5(j))
        fs1(j) = dw1
        fs2(j) = dw2
        fs3(j) = dw3
        fs4(j) = dw4
        fs5(j) = dw5
  185 continue
      endif
      endif
c
      do 199 j=jmn+1,jmx
      fw(i,j,k,1)  = fw(i,j,k,1)  + bet*(fs1(j-1) - fs1(j))
      fw(i,j,k,2)  = fw(i,j,k,2)  + bet*(fs2(j-1) - fs2(j))
      fw(i,j,k,3)  = fw(i,j,k,3)  + bet*(fs3(j-1) - fs3(j))
      fw(i,j,k,4)  = fw(i,j,k,4)  + bet*(fs4(j-1) - fs4(j))
      fw(i,j,k,5)  = fw(i,j,k,5)  + bet*(fs5(j-1) - fs5(j))
  199 continue
c
c--------- dissipation in k direction 
c
      if(i2dfl.eq.1) return
c
      do 299 j=jmn+1,jmx
      do 299 i=imn+1,imx
      do k=kmn,kmp1
        reigjk   = ratiojk(i,j,k)                  
        reigik   = ratioij(i,j,k)*ratiojk(i,j,k)
        phi(k)   =  1. + reigik + reigjk
      enddo
c
      if (iprec.le.1) then
c
c            form first differences of conservative variables
c
      do k=kmn,kmx
        fila      = 0.5*((dtk(i,j,k+1)+difac*dtvk(i,j,k+1))*phi(k+1)       
     .                 + (dtk(i,j,k  )+difac*dtvk(i,j,k  ))*phi(k  ))
        dis       = fis0*fila
        fs1(k)    = dis*(w(i,j,k+1,1) - w(i,j,k,1)) 
        fs2(k)    = dis*(w(i,j,k+1,2) - w(i,j,k,2)) 
        fs3(k)    = dis*(w(i,j,k+1,3) - w(i,j,k,3)) 
        fs4(k)    = dis*(w(i,j,k+1,4) - w(i,j,k,4)) 
        fs5(k)    = dis*(w(i,j,k+1,5) - w(i,j,k,5)  
     .       + (1-icau)*(p(i,j,k+1)   - p(i,j,k)  ))
      enddo
c
        else
c
c     form first differences of "u,v,w,p,T" variables
c
c
        rhoi    = 1./w(i,j,kmn,1)
        do k=kmn,kmx
          fila    = 0.5*((dtk(i,j,k+1)+difac*dtvk(i,j,k+1))*phi(k+1)       
     .                 + (dtk(i,j,k  )+difac*dtvk(i,j,k  ))*phi(k  ))
          dis     = fis0*fila
          rhoip   = 1./w(i,j,k+1,1)
          fs1(k)  = dis*(p(i,j,k+1  )       - p(i,j,k  )     ) 
          fs2(k)  = dis*(w(i,j,k+1,2)*rhoip - w(i,j,k,2)*rhoi) 
          fs3(k)  = dis*(w(i,j,k+1,3)*rhoip - w(i,j,k,3)*rhoi) 
          fs4(k)  = dis*(w(i,j,k+1,4)*rhoip - w(i,j,k,4)*rhoi) 
          fs5(k)  = dis*(p(i,j,k+1  )*rhoip - p(i,j,k  )*rhoi)  
          rhoi    = rhoip
        enddo
c
        endif
c
      if (iprec.eq.1) then
c            for conservative formulation multiply by P^{-1}
      do 270 k=kmn,kmx
        rhoi      = 1./w(i,j,k,1)
        uc        = w(i,j,k,2)*rhoi
        vc        = w(i,j,k,3)*rhoi
        wc        = w(i,j,k,4)*rhoi
        q2        = 0.5*(uc**2 + vc**2 + wc**2)
        ts        = p(i,j,k)*rhoi
        c2        = gamma*ts
        hc        = cp*ts + q2
        bc        = fbeta2(i,j,k)
        delta     = deltpre
        z         = gm1*(q2*fs1(k) - uc*fs2(k)
     .            - vc*fs3(k) - wc*fs4(k)
     .            + fs5(k) )/c2
        dw1       = ((1. - bc)*z - delta*(z - fs1(k)))/bc
        fs1(k)    = fs1(k) + dw1
        fs2(k)    = fs2(k) + uc*dw1
        fs3(k)    = fs3(k) + vc*dw1
        fs4(k)    = fs4(k) + wc*dw1
        fs5(k)    = fs5(k) + hc*dw1
  270 continue
      endif
      if (iprec.eq.2) then
c         for conservative formulation multiply by Gamma^{-1}
      if (delta.eq.0.) then
      do 280 k=kmn,kmx
        rhoa    = 0.5*(w(i,j,k,1)+w(i,j,k+1,1))
        rhoi    = 1./rhoa
        rhou    = 0.5*(w(i,j,k,2) + w(i,j,k+1,2))
        rhov    = 0.5*(w(i,j,k,3) + w(i,j,k+1,3))
        rhow    = 0.5*(w(i,j,k,4) + w(i,j,k+1,4))
        ts      = 0.5*(p(i,j,k) + p(i,j,k+1))*rhoi
        hqs     = 0.5*(rhou*rhou+rhov*rhov+rhow*rhow)*rhoi*rhoi
        c2      = gamma*ts
        dc2     = 1./c2
        hh      = c2/gm1 + hqs
c
c               evaluate alpha and beta for preconditioning
c
c       bc      = 0.5*(fbeta2(i,j,k) + fbeta2(i,j,k+1))
        bc      = max (fbeta2(i,j,k) , fbeta2(i,j,k+1))
        alp     = max(ulprec*(1.-bc),0.)
        dbeta2  = 1./bc
        theta   = (gm1 + dbeta2)*dc2
        rtq5    = rhoa*fs5(k)/ts
        aq1     = alp*fs1(k)*dbeta2*dc2
c
        dw1     = theta*fs1(k) - rtq5
        dwpt    = rhoi*(dw1 + aq1)
        dw2     = rhoa*fs2(k) + rhou*dwpt
        dw3     = rhoa*fs3(k) + rhov*dwpt
        dw4     = rhoa*fs4(k) + rhow*dwpt
        dw5     = (theta*hh-1.)*fs1(k)
     .          + rhou*fs2(k) + rhov*fs3(k) + rhow*fs4(k)
     .          + hqs*(2.*aq1 - rtq5)
c
        fs1(k) = dw1
        fs2(k) = dw2
        fs3(k) = dw3
        fs4(k) = dw4
        fs5(k) = dw5
  280 continue
      else
      do 285 k=kmn,kmx
        rhoa    = 0.5*(w(i,j,k,1)+w(i,j,k+1,1))
        rhoi    = 1./rhoa
        uu      = 0.5*(w(i,j,k,2) + w(i,j,k+1,2))*rhoi
        vv      = 0.5*(w(i,j,k,3) + w(i,j,k+1,3))*rhoi
        ww      = 0.5*(w(i,j,k,4) + w(i,j,k+1,4))*rhoi
        ts      = 0.5*(p(i,j,k) + p(i,j,k+1))*rhoi
        hh      = 0.5*(w(i,j,k,5) + w(i+1,j,k,5))*rhoi + ts
        c2      = gamma*ts
        us      = uu*uu + vv*vv + ww*ww
c
c               evaluate alpha and beta for preconditioning
c
c       bc      = 0.5*(fbeta2(i,j,k) + fbeta2(i,j,k+1))
        bc      = max (fbeta2(i,j,k) , fbeta2(i,j,k+1))
        beta2   = bc*c2
        alp     = max(ulprec*(1.-bc),0.)
        dbeta2  = 1./beta2
        alp1    = (1.+alp)*dbeta2*fs1(k)
c
        dw1     = dbeta2*fs1(k)
        dw2     = rhoa*fs2(k) + alp1*uu
        dw3     = rhoa*fs3(k) + alp1*vv
        dw4     = rhoa*fs4(k) + alp1*ww
        g41     = dbeta2*(hh + alp*us) -1.
        dw5     = g41*fs1(k)
     .          + rhoa*(uu*fs2(k) + vv*fs3(k) + ww*fs4(k)
     .          + cp*fs5(k))
        fs1(k) = dw1
        fs2(k) = dw2
        fs3(k) = dw3
        fs4(k) = dw4
        fs5(k) = dw5
  285 continue
      endif
      endif
c
      do 299 k=kmn+1,kmx
      fw(i,j,k,1)  = fw(i,j,k,1)  + bet*(fs1(k-1) - fs1(k))
      fw(i,j,k,2)  = fw(i,j,k,2)  + bet*(fs2(k-1) - fs2(k))
      fw(i,j,k,3)  = fw(i,j,k,3)  + bet*(fs3(k-1) - fs3(k))
      fw(i,j,k,4)  = fw(i,j,k,4)  + bet*(fs4(k-1) - fs4(k))
      fw(i,j,k,5)  = fw(i,j,k,5)  + bet*(fs5(k-1) - fs5(k))
  299 continue
c
      return
      end
