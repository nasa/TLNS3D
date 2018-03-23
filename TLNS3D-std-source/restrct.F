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
      subroutine restrct(iftf,igrid,imn,jmn,kmn,
     .                   imxf, jmxf, kmxf, imp1f,jmp1f,kmp1f,
     .                   imp2f,jmp2f,kmp2f,imxc, jmxc, kmxc,
     .                   imp1c,jmp1c,kmp1c,imp2c,jmp2c,kmp2c,
     .                   wf,wc,pf,pc,delwf,wrc,volf,ibloc              )
c
c

      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension wf(imp2f,jmp2f,kmp2f,5),  wc(imp2c,jmp2c,kmp2c,5),
     .          pf(imp2f,jmp2f,kmp2f  ),  pc(imp2c,jmp2c,kmp2c  ),
     .       delwf(imp2f,jmp2f,kmp2f,5), wrc(imp2c,jmp2c,kmp2c,5),
     .        volf(imp2f,jmp2f,kmp2f  )
c
c
      if(ibloc.eq.1) iftf = iftf +1
c     write(*,'("in restrct, igrid from",i4," to",i4," iftf=",i4)') 
c    &           igrid-1,igrid,iftf
c
      kmnp1 = kmn + 1
      jmnp1 = jmn + 1
      imnp1 = imn + 1
c
c     restrict the conserved variables from the fine to the coarse
c     grid using volume weighting, at the interior points
c
        do 5 n=1,5
        do 5 kc=kmnp1,kmxc
c
c     because of the memory shift necessary to allow each block
c     to have two ghost cells
c     use        kf=  kc*2 -kmnp1
c     instead of kf=  kc*2 -2, and so on for jc and ic
c
      kf= kc*2 - kmnp1
c
        do 5 jc=jmnp1,jmxc
      jf= jc*2 - jmnp1
c
        do 5 ic=imnp1,imxc
      if= ic*2 - imnp1
c
      v1 = volf(if  ,jf  ,kf  )
      v2 = volf(if+1,jf  ,kf  )
      v3 = volf(if  ,jf+1,kf  )
      v4 = volf(if+1,jf+1,kf  )
      v5 = volf(if  ,jf  ,kf+1)
      v6 = volf(if+1,jf  ,kf+1)
      v7 = volf(if  ,jf+1,kf+1)
      v8 = volf(if+1,jf+1,kf+1)
      vc = v1 +v2 +v3 +v4 +v5 +v6 +v7 +v8
c
      wc(ic,jc,kc,n)=(wf(if,jf  ,kf  ,n)*v1+wf(if+1,jf  ,kf  ,n)*v2+
     .                wf(if,jf+1,kf  ,n)*v3+wf(if+1,jf+1,kf  ,n)*v4+
     .                wf(if,jf  ,kf+1,n)*v5+wf(if+1,jf  ,kf+1,n)*v6+
     .                wf(if,jf+1,kf+1,n)*v7+wf(if+1,jf+1,kf+1,n)*v8)/vc
c
    5 continue
c
c     evaluate the coarse-grid flow-variables at all face
c     boundaries based on fine grid solutions
c     note that :
c     except for far-field type of boundaries, these variables
c     will be recomputed with the appropriate b.c
c     through the call to b.c. routines later in the main program
c
      do 10 n = 1,5
      do 10 kc=kmn,kmp1c
      do 10 jc=jmn,jmp1c
c
      wc(imn  ,jc,kc,n) = wc(imnp1,jc,kc,n)
      wc(imp1c,jc,kc,n) = wc(imxc ,jc,kc,n)
c
  10  continue
c
      do 15 n = 1,5
      do 15 kc=kmn,kmp1c
      do 15 ic=imn,imp1c
c
      wc(ic,jmn  ,kc,n) = wc(ic,jmnp1,kc,n)
      wc(ic,jmp1c,kc,n) = wc(ic,jmxc ,kc,n)
c
  15  continue
c
      do 20 n = 1,5
      do 20 jc=jmn,jmp1c
      do 20 ic=imn ,imp1c
c
      wc(ic,jc,kmn  ,n) = wc(ic,jc,kmnp1,n)
      wc(ic,jc,kmp1c,n) = wc(ic,jc,kmxc ,n)
c
  20  continue
c
c     calculate the pressure on the coarse grid
c
      do 30 kc=kmn,kmp1c
      do 30 jc=jmn,jmp1c
      do 30 ic=imn,imp1c
      pc(ic,jc,kc) = gm1*(wc(ic,jc,kc,5)
     .               -.5*(wc(ic,jc,kc,2)**2 + wc(ic,jc,kc,3)**2 +
     .                    wc(ic,jc,kc,4)**2)/ wc(ic,jc,kc,1)     )
      pc(ic,jc,kc) = max(pc(ic,jc,kc),0.001)
c
   30 continue
c
c     initialize the coarse grid forcing function
c
      do 40 n=1,5
      do 40 kc=1,kmp2c
      do 40 jc=1,jmp2c
      do 40 ic=1,imp2c
c
      wrc(ic,jc,kc,n) = 0.
c
  40  continue
c
c     restrict the fine grid residual to the coarse grid and
c     store in the coarse grid forcing function
c
        do 45 n=1,5
        do 45 kc=kmnp1,kmxc
      kf= kc*2 -kmnp1
      kfp1 = kf + 1
c
      if(i2dfl.eq.1) then
        kf  = kmnp1
        kfp1= kf
      endif
c
        do 45 jc=jmnp1,jmxc
      jf= jc*2 -jmnp1
        do 45 ic=imnp1,imxc
      if= ic*2 -imnp1
c
      wrc(ic,jc,kc,n) = 
     .            delwf(if,jf  ,kf  ,n) + delwf(if+1,jf  ,kf  ,n) +
     .            delwf(if,jf+1,kf  ,n) + delwf(if+1,jf+1,kf  ,n) +
     . (1-i2dfl)*(delwf(if,jf  ,kfp1,n) + delwf(if+1,jf  ,kfp1,n) +
     .            delwf(if,jf+1,kfp1,n) + delwf(if+1,jf+1,kfp1,n) )
c
   45 continue
c
c
      return
      end
