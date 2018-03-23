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
      subroutine bbdist(ng,grid,nsurf,surf,
     .           nbb,bbdef,ipv,vlist,dist,idist,ncalc,wk3d5)
      dimension grid(ng,3),surf(nsurf,3),dist(ng),idist(ng)
      dimension bbdef(6,nbb),ipv(2,nbb)
      dimension wk3d5(*)
      integer   vlist(nsurf)
      integer wrk,tsurf,test
c      pointer (p_wrk,wrk(nbb)),(p_tsurf,tsurf(nsurf,3))
c      pointer (p_test,test(nsurf))
c      p_wrk=my_malloc(nbb)
c      p_tsurf=my_malloc(3*nsurf)
c      p_test=my_malloc(nsurf)
      wrk = ifalloc(nbb)
      tsurf = ifalloc(3*nsurf)
      test = ifalloc(nsurf)
      call bbdst1(ng,grid,nsurf,surf,
     .           nbb,bbdef,ipv,vlist,dist,idist,ncalc,
     .           wk3d5(wrk),wk3d5(tsurf),wk3d5(test))
      call ffree(nsurf)
      call ffree(3*nsurf)
      call ffree(nbb)
c
      return
      end
c
      subroutine bbdst1(ng,grid,nsurf,surf,
     .           nbb,bbdef,ipv,vlist,dist,idist,ncalc,
     .           wrk,tsurf,test)
      dimension grid(ng,3),surf(nsurf,3),dist(ng),idist(ng)
      dimension bbdef(6,nbb),ipv(2,nbb)
      dimension wrk(nbb), tsurf(nsurf,3), test(nsurf)
      integer   vlist(nsurf)
      do i=1,nsurf
      kk=vlist(i)
      tsurf(i,1)=surf(kk,1)
      tsurf(i,2)=surf(kk,2)
      tsurf(i,3)=surf(kk,3)
      end do
      ncalc=0
      do 100 i=1,ng
      x=grid(i,1)
      y=grid(i,2)
      z=grid(i,3)
      smin=1.0d34
      jp=0
c   build table of distances from grid to bounding boxes
      do 120 j=1,nbb
      px=x
      if (px .le. bbdef(1,j)) px=bbdef(1,j)
      if (px .ge. bbdef(2,j)) px=bbdef(2,j)
      py=y
      if (py .le. bbdef(3,j)) py=bbdef(3,j)
      if (py .ge. bbdef(4,j)) py=bbdef(4,j)
      pz=z
      if (pz .le. bbdef(5,j)) pz=bbdef(5,j)
      if (pz .ge. bbdef(6,j)) pz=bbdef(6,j)
      wrk(j)=(x-px)**2+(y-py)**2+(z-pz)**2
 120  continue
      do 200 j=1,nbb
c   find nearest bounding box that has not been searched
      bbmin=wrk(1)
      do ii=2,nbb
      bbmin=min(bbmin,wrk(ii))
      end do
      jj=isrcheq(nbb,wrk,1,bbmin)
      wrk(jj)=2.0d34
c   stop searching when nearest bounding box is too far away
      if (bbmin .gt. smin) goto 201
      n=ipv(1,jj)
      l=ipv(2,jj)
      testmin=1.0e34
      do 160 k=1,n
      xs=tsurf(l,1)
      ys=tsurf(l,2)     
      zs=tsurf(l,3)
      l=l+1
      test(k)=(x-xs)**2+(y-ys)**2+(z-zs)**2
      testmin=min(testmin,test(k))
 160  continue
      ncalc=ncalc+n
      if (testmin .lt. smin) then
        kmin=isrcheq(n,test,1,testmin)
        smin=testmin
        jp=vlist(ipv(2,jj)+kmin-1)
      end if
 200  continue
 201  continue
      dist(i)=sqrt(smin)
      idist(i)=jp
 100  continue
c      call my_free(p_wrk,nbb)
c      call my_free(p_tsurf,3*nsurf)
c      call my_free(p_test,nsurf)
      return
      end
