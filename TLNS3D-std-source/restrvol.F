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
      subroutine restrvol(igrid,imn,jmn,kmn,
     .                   imxf, jmxf, kmxf, imp1f,jmp1f,kmp1f,
     .                   imp2f,jmp2f,kmp2f,imxc, jmxc, kmxc,
     .                   imp1c,jmp1c,kmp1c,imp2c,jmp2c,kmp2c,
     .                   volf,volc)
c
c

      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension volf(imp2f,jmp2f,kmp2f),  volc(imp2c,jmp2c,kmp2c)
c
c
c
      kmnp1 = kmn + 1
      jmnp1 = jmn + 1
      imnp1 = imn + 1
c
c     restrict the volume from the fine to the coarse grid
c
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
      volc(ic,jc,kc)= vc
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
      do 10 kc=kmn,kmp1c
      do 10 jc=jmn,jmp1c
c
      volc(imn  ,jc,kc) =  volc(imnp1,jc,kc)
      volc(imp1c,jc,kc) =  volc(imxc ,jc,kc)
c
  10  continue
c
      do 15 kc=kmn,kmp1c
      do 15 ic=imn,imp1c
c
      volc(ic,jmn  ,kc) =  volc(ic,jmnp1,kc)
      volc(ic,jmp1c,kc) =  volc(ic,jmxc ,kc)
c
  15  continue
c
      do 20 jc=jmn,jmp1c
      do 20 ic=imn ,imp1c
c
      volc(ic,jc,kmn  ) =  volc(ic,jc,kmnp1)
      volc(ic,jc,kmp1c) =  volc(ic,jc,kmxc )
c
  20  continue
c
c
      return
      end
