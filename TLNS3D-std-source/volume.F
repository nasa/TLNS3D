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
      subroutine  volume(imn, jmn, kmn, imx,jmx,kmx,imp1,jmp1,kmp1,
     .                   imp2,jmp2,kmp2,vol,x,igrid,ibloc)
c=======================================================================
c 
c             calculates cell volumes
c 
c 

      dimension   vol(imp2,jmp2,kmp2),x(imp1,jmp1,kmp1,3)
c
c--------- volume of a pyramid
c
      volpym(xp,yp,zp,xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd) =
     .      (xp  -.25*(xa  +xb  +xc  +xd))
     .      *((ya  -yc)*(zb  -zd)  -(za  -zc)*(yb  -yd))
     .      +(yp  -.25*(ya  +yb  +yc  +yd))
     .      *((za  -zc)*(xb  -xd)  -(xa  -xc)*(zb  -zd))
     .      +(zp  -.25*(za  +zb  +zc  +zd))
     .      *((xa  -xc)*(yb  -yd)  -(ya  -yc)*(xb  -xd))
c
c
c     do  5 k=kmn,kmp1
c     do  5 j=jmn,jmp1
c     do  5 i=imn,imp1
      do  5 k=1,kmp2
      do  5 j=1,jmp2
      do  5 i=1,imp2
      vol(i,j,k)  =  1.0
    5 continue
      vol(1,1,1) = -3.5
c
c
      negv = 0
      do 10 k=kmn+1,kmx
      do 10 j=jmn+1,jmx
      do 10 i=imn+1,imx
      xp        = .125*(x(i,j  ,k  ,1) + x(i-1,j  ,k  ,1) 
     .                 +x(i,j-1,k  ,1) + x(i-1,j-1,k  ,1) 
     .                 +x(i,j  ,k-1,1) + x(i-1,j  ,k-1,1) 
     .                 +x(i,j-1,k-1,1) + x(i-1,j-1,k-1,1))
      yp        = .125*(x(i,j  ,k  ,2) + x(i-1,j  ,k  ,2) 
     .                 +x(i,j-1,k  ,2) + x(i-1,j-1,k  ,2) 
     .                 +x(i,j  ,k-1,2) + x(i-1,j  ,k-1,2) 
     .                 +x(i,j-1,k-1,2) + x(i-1,j-1,k-1,2))
      zp        = .125*(x(i,j  ,k  ,3) + x(i-1,j  ,k  ,3) 
     .                 +x(i,j-1,k  ,3) + x(i-1,j-1,k  ,3) 
     .                 +x(i,j  ,k-1,3) + x(i-1,j  ,k-1,3) 
     .                 +x(i,j-1,k-1,3) + x(i-1,j-1,k-1,3))
c 
      vp1       = volpym(xp,yp,zp,
     .            x(i  ,j  ,k  ,1),x(i  ,j  ,k  ,2),x(i  ,j  ,k  ,3), 
     .            x(i  ,j-1,k  ,1),x(i  ,j-1,k  ,2),x(i  ,j-1,k  ,3), 
     .            x(i  ,j-1,k-1,1),x(i  ,j-1,k-1,2),x(i  ,j-1,k-1,3), 
     .            x(i  ,j  ,k-1,1),x(i  ,j  ,k-1,2),x(i  ,j  ,k-1,3)) 
      vp2       = volpym(xp,yp,zp,
     .            x(i-1,j  ,k  ,1),x(i-1,j  ,k  ,2),x(i-1,j  ,k  ,3), 
     .            x(i-1,j-1,k  ,1),x(i-1,j-1,k  ,2),x(i-1,j-1,k  ,3), 
     .            x(i-1,j-1,k-1,1),x(i-1,j-1,k-1,2),x(i-1,j-1,k-1,3), 
     .            x(i-1,j  ,k-1,1),x(i-1,j  ,k-1,2),x(i-1,j  ,k-1,3)) 
      vp3       = volpym(xp,yp,zp,
     .            x(i  ,j  ,k  ,1),x(i  ,j  ,k  ,2),x(i  ,j  ,k  ,3), 
     .            x(i  ,j  ,k-1,1),x(i  ,j  ,k-1,2),x(i  ,j  ,k-1,3), 
     .            x(i-1,j  ,k-1,1),x(i-1,j  ,k-1,2),x(i-1,j  ,k-1,3), 
     .            x(i-1,j  ,k  ,1),x(i-1,j  ,k  ,2),x(i-1,j  ,k  ,3)) 
      vp4       = volpym(xp,yp,zp,
     .            x(i  ,j-1,k  ,1),x(i  ,j-1,k  ,2),x(i  ,j-1,k  ,3), 
     .            x(i  ,j-1,k-1,1),x(i  ,j-1,k-1,2),x(i  ,j-1,k-1,3), 
     .            x(i-1,j-1,k-1,1),x(i-1,j-1,k-1,2),x(i-1,j-1,k-1,3), 
     .            x(i-1,j-1,k  ,1),x(i-1,j-1,k  ,2),x(i-1,j-1,k  ,3)) 
      vp5       = volpym(xp,yp,zp,
     .            x(i  ,j  ,k  ,1),x(i  ,j  ,k  ,2),x(i  ,j  ,k  ,3), 
     .            x(i-1,j  ,k  ,1),x(i-1,j  ,k  ,2),x(i-1,j  ,k  ,3), 
     .            x(i-1,j-1,k  ,1),x(i-1,j-1,k  ,2),x(i-1,j-1,k  ,3), 
     .            x(i  ,j-1,k  ,1),x(i  ,j-1,k  ,2),x(i  ,j-1,k  ,3)) 
      vp6       = volpym(xp,yp,zp,
     .            x(i  ,j  ,k-1,1),x(i  ,j  ,k-1,2),x(i  ,j  ,k-1,3), 
     .            x(i-1,j  ,k-1,1),x(i-1,j  ,k-1,2),x(i-1,j  ,k-1,3), 
     .            x(i-1,j-1,k-1,1),x(i-1,j-1,k-1,2),x(i-1,j-1,k-1,3), 
     .            x(i  ,j-1,k-1,1),x(i  ,j-1,k-1,2),x(i  ,j-1,k-1,3)) 
      vol(i,j,k)    = (-vp1  +vp2  -vp3  +vp4  -vp5  +vp6)/6.

      if(vol(i,j,k).le.0.0.and.negv.le.20) then
        negv  = negv + 1
        write (6,901) i,j,k,igrid,ibloc,vol(i,j,k)
      endif
 901  format(2x,'neg vol at i j k igrid ibloc', 5i5, e15.6)
c       vol(i,j,k) = abs(vol(i,j,k))
   10 continue
c 
c
      return
      end
