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
      subroutine calc_dist (imn,jmn,kmn,
     .    imx  , jmx , kmx  ,
     .    imp1, jmp1, kmp1,
     .    imp2, jmp2, kmp2,
     .       x, smin,
     .    nsurf, surf, nbb, bbdef, ipv, vlist,
     .    ntri, iptri ,wk3d5,iwrk)
      dimension x(imp1,jmp1,kmp1,3),smin(imp2,jmp2,kmp2)
      dimension wk3d5(*), iwrk(*)
      integer vlist,grid,dist
      dimension surf(nsurf,3),bbdef(6,nbb),ipv(2,nbb),vlist(nsurf)
      dimension ntri(nsurf),iptri(nsurf,8)
c      pointer (p_grid,grid(1)),(p_dist,dist(1))
c      pointer (p_idist,idist(1))
      ng=(imx-imn+1)*(jmx-jmn+1)*(kmx-kmn+1)
c      print *,"ng=",ng
c      p_grid=my_malloc(3*ng)
c      p_dist=my_malloc(ng)
c      p_idist=my_malloc(ng)
      grid = ifalloc(3*ng)
      dist = ifalloc(ng)
      idist = iialloc(ng)
      ig=0
      do k=kmn,kmx
      do j=jmn,jmx
CDIR$ IVDEP
      do i=imn,imx
        ig=ig+1
        wk3d5(grid+ig-1)=x(i,j,k,1)
        wk3d5(grid+ig+ng-1)=x(i,j,k,2)
        wk3d5(grid+ig+2*ng-1)=x(i,j,k,3)
      end do
      end do
      end do
c      time1=second()
      time1=0.0
      call bbdist(ng,wk3d5(grid),nsurf,surf,
     .            nbb,bbdef,ipv,vlist,wk3d5(dist),iwrk(idist),
     .            ncalc,wk3d5)
c      time2=second()
      time2=1.0
c      print *,"time for finding nearest points=",time2-time1
c
c  put in calculation to triangles
c
c      time1=second()
      time1=0.0
      do i=1,ng
      isurf=iwrk(idist+i-1)
      xp=wk3d5(grid+i-1)
      yp=wk3d5(grid+i+ng-1)
      zp=wk3d5(grid+i+2*ng-1)
      numtri=ntri(isurf)
      if (numtri .eq. 4) then
        do itri=1,numtri
        it1=iptri(isurf,2*itri-1)
        it2=iptri(isurf,2*itri)
        call triang(wk3d5(dist+i-1),xp,yp,zp,
     .              surf(isurf,1),surf(isurf,2),surf(isurf,3),
     .              surf(it1,1),surf(it1,2),surf(it1,3),
     .              surf(it2,1),surf(it2,2),surf(it2,3))
        end do
      else
        xps=surf(isurf,1)
        yps=surf(isurf,2)
        zps=surf(isurf,3)
        do ixmin=isurf,1,-1
        if (surf(ixmin,1) .lt. xps) go to 10
        end do
        ixmin=0
 10     continue
        ixmin=ixmin+1
        do ixmax=isurf,nsurf
        if (surf(ixmax,1) .gt. xps) go to 20
        end do
        ixmax=nsurf+1
 20     continue
        ixmax=ixmax-1
        do isurf=ixmin,ixmax
          if (surf(isurf,2) .eq. yps .and. surf(isurf,3) .eq. zps
     .        .and. ntri(isurf) .ne. 4) then
            numtri=ntri(isurf)
            do itri=1,numtri
            it1=iptri(isurf,2*itri-1)
            it2=iptri(isurf,2*itri)
            call triang(wk3d5(dist+i-1),xp,yp,zp,
     .              surf(isurf,1),surf(isurf,2),surf(isurf,3),
     .              surf(it1,1),surf(it1,2),surf(it1,3),
     .              surf(it2,1),surf(it2,2),surf(it2,3))
            end do
          end if
        end do
      end if
      end do
c      time2=second()
      time2=1.0
c      print *,"time for triangles=",time2-time1
      ig=0
      do k=kmn,kmx
      do j=jmn,jmx
      do i=imn,imx
      ig=ig+1
      smin(i,j,k)=min(smin(i,j,k),wk3d5(dist+ig-1))
      end do
      end do
      end do
c      write(ismfil) (((smin(i,j,k),i=imn,imx),j=jmn,jmx),k=kmn,kmx)
c      call my_free(p_grid,3*ng)
c      call my_free(p_dist,ng)
c      call my_free(p_idist,ng)
      call ifree(ng)
      call ffree(ng)
      call ffree(3*ng)
      return
      end
