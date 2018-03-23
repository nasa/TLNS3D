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
      subroutine distcc(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                  imp2,jmp2,kmp2, smin,wk3d              )
c
c     object : evaluate distances to cell-centers from walls 
c              by averaging node-center values
c
c     initial coding : by v.n.vatsa of nasa langley (March, 1993)
c                      modified to reduce sqrt operations (Nov., 1993)
c

      dimension smin(imp2,jmp2,kmp2),wk3d(imp1,jmp1,kmp1)
c
c
c     accomodates two cell faces at each boundary
c
c
        do 10 k= kmn,kmx
        do 10 j= jmn,jmx
        do 10 i= imn,imx
c
c      (vatsa sqrt is already taken in  computation of dismin)
        wk3d(i,j,k) = smin(i,j,k)
ccc     wk3d(i,j,k) = sqrt(smin(i,j,k))
c
   10 continue
c
        do 20 k= kmn+1,kmx
        do 20 j= jmn+1,jmx
        do 20 i= imn+1,imx
      smin(i,j,k)  = .125*(wk3d(i  ,j,k  ) + wk3d(i  ,j-1,k  )
     .                    +wk3d(i  ,j,k-1) + wk3d(i  ,j-1,k-1)
     .                    +wk3d(i-1,j,k  ) + wk3d(i-1,j-1,k  )
     .                    +wk3d(i-1,j,k-1) + wk3d(i-1,j-1,k-1))
   20 continue
c
c     fill-in ghost-cell values for b..c routines
      do 30 k=kmn+1,kmx
      do 30 j=jmn+1,jmx
      smin(imn,j,k)   = smin(imn+1,j,k)
      smin(imx+1,j,k) = smin(imx  ,j,k)
 30   continue
      do 40 k=kmn+1,kmx
      do 40 i=imn,imx+1
      smin(i,jmn,k)   = smin(i,jmn+1,k)
      smin(i,jmx+1,k) = smin(i,jmx  ,k)
 40   continue
      do 50 j=jmn,jmx+1
      do 50 i=imn,imx+1
      smin(i,j,kmn)   = smin(i,j,kmn+1)
      smin(i,j,kmx+1) = smin(i,j,kmx  )
 50   continue
c
c
      return
      end
