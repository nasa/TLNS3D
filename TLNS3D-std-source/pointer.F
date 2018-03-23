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
      subroutine pointer( ngrid,m1np,m2np,m3np,m4np,m5np,
     .                          m1cc,m2cc,m3cc,m4cc,m5cc,m5cg,
     .                          m1is,m1js,m1ks,
     .                          m3is,m3js,m3ks,im,jm,km,
     .                          imp1,jmp1,kmp1,imp2,jmp2,kmp2 )
c

      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension      im(ngrid),  jm(ngrid),  km(ngrid),
     .             imp1(ngrid),jmp1(ngrid),kmp1(ngrid),
     .             imp2(ngrid),jmp2(ngrid),kmp2(ngrid),
     . m1np(ngrid),m2np(ngrid),m3np(ngrid),m4np(ngrid),m5np(ngrid),
     . m1cc(ngrid),m2cc(ngrid),m3cc(ngrid),m4cc(ngrid),m5cc(ngrid),
     . m1is(ngrid),m1js(ngrid),m1ks(ngrid),
     . m3is(ngrid),m3js(ngrid),m3ks(ngrid),m5cg(ngrid)
c
c     generate the node point limits (imp1 etc.)
c     and the cell-center limits     (imp2 etc.)
c
      imp1(1) =   im(1) + 1
      jmp1(1) =   jm(1) + 1
      kmp1(1) =   km(1) + 1
      imp2(1) = imp1(1) + 1
      jmp2(1) = jmp1(1) + 1
      kmp2(1) = kmp1(1) + 1
c
c     i2dfl = 0
c
        do 5 igrid = 2,ngrid
c
c
c     because the memory had to be shifted to account for two
c     memory cells at each block boundary
c   # of interior cell-centers in physical domain are: im-2
c
c
        im(igrid) =  (im(igrid-1) - 2)/2 + 2
      imp1(igrid) =   im(igrid)   + 1
      imp2(igrid) = imp1(igrid)   + 1
c
        jm(igrid) =  (jm(igrid-1) - 2)/2 + 2
      jmp1(igrid) =   jm(igrid)   + 1
      jmp2(igrid) = jmp1(igrid)   + 1
c
        km(igrid) =  (km(igrid-1) - 2)/2 + 2
c
      if (i2dfl .eq.1) km(igrid) = km(igrid-1)
      kmp1(igrid) =   km(igrid)   + 1
      kmp2(igrid) = kmp1(igrid)   + 1
c
    5 continue
c
c   generate the pointers for node-point & cell-center 3-d arrays
c
c
      do 10 igrid = 2,ngrid
c
      ltempnp = imp1(igrid-1)*jmp1(igrid-1)*kmp1(igrid-1)
      ltempcc = imp2(igrid-1)*jmp2(igrid-1)*kmp2(igrid-1)
c
        m1np(igrid) = m1np(igrid-1) + ltempnp
        m2np(igrid) = m2np(igrid-1) + ltempnp*2
        m3np(igrid) = m3np(igrid-1) + ltempnp*3
        m4np(igrid) = m4np(igrid-1) + ltempnp*4
        m5np(igrid) = m5np(igrid-1) + ltempnp*5
c
        m1cc(igrid) = m1cc(igrid-1) + ltempcc
        m2cc(igrid) = m2cc(igrid-1) + ltempcc*2
        m3cc(igrid) = m3cc(igrid-1) + ltempcc*3
        m4cc(igrid) = m4cc(igrid-1) + ltempcc*4
        m5cc(igrid) = m5cc(igrid-1) + ltempcc*5
c
   10 continue
c
c
      do 20 igrid = 2,ngrid
c
        ltempcc = imp2(igrid-1)*jmp2(igrid-1)*kmp2(igrid-1)
c
        m5cg(igrid) = m5cg(igrid-1) + ltempcc*5
c
c       reset this pointer at igrid=2, since it
c       is not used for the finest grid (igrid=1)
c
        if(igrid .eq. 2) m5cg(igrid) = m5cg(igrid-1)
c
   20 continue
c
c   generate pointers for 2-d or surface arrays
c
      do 30 igrid = 2,ngrid
c
      itemp = jmp1(igrid-1)*kmp1(igrid-1)
      jtemp = imp1(igrid-1)*kmp1(igrid-1)
      ktemp = imp1(igrid-1)*jmp1(igrid-1)
c
      m1is(igrid) = m1is(igrid-1) + itemp
      m1js(igrid) = m1js(igrid-1) + jtemp
      m1ks(igrid) = m1ks(igrid-1) + ktemp
c
      m3is(igrid) = m3is(igrid-1) + itemp*3
      m3js(igrid) = m3js(igrid-1) + jtemp*3
      m3ks(igrid) = m3ks(igrid-1) + ktemp*3
c
  30  continue
c
c
      return
      end
