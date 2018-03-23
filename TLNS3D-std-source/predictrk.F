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
      subroutine predictrk(imn,jmn,kmn,imx,jmx,kmx,
     .                     imp1,jmp1,kmp1,imp2,jmp2,kmp2,irkcase,
     .                     nrk,dto,dtn,rkphysa,vol,wold,delwold)
#ifdef  RKPHY
      dimension bd(10,3),c(10),rkphysa(nrk,nrk),bb(10,10)
      dimension wEX(10)
      dimension wold(imp2,jmp2,kmp2,5),delwold(imp2,jmp2,kmp2,5,nrk)
      dimension vol (imp2,jmp2,kmp2)
c
        if(nrk.gt.10) return
        if(irkcase.ne.5) return

c     3rd-order intrp coef valid for method ``5''

        bd(1,1) = 6943876665148.d0/7220017795957.d0
        bd(2,1) = 0.d0
        bd(3,1) = 7640104374378.d0/9702883013639.d0
        bd(4,1) = -20649996744609.d0/7521556579894.d0
        bd(5,1) = 8854892464581.d0/2390941311638.d0
        bd(6,1) =-11397109935349.d0/6675773540249.d0
        bd(1,2) =-54480133.d0/30881146.d0
        bd(2,2) = 0.d0
        bd(3,2) =-11436875.d0/14766696.d0
        bd(4,2) = 174696575.d0/18121608.d0
        bd(5,2) =-12120380.d0/966161.d0
        bd(6,2) = 3843.d0/706.d0
        bd(1,3) = 6818779379841.d0/7100303317025.d0
        bd(2,3) = 0.d0
        bd(3,3) = 2173542590792.d0/12501825683035.d0
        bd(4,3) =-31592104683404.d0/5083833661969.d0
        bd(5,3) =+61146701046299.d0/7138195549469.d0
        bd(6,3) =-17219254887155.d0/4939391667607.d0

      do i = 1,nrk
        c(i) = 0.0
        do j = 1,nrk
          c(i) = c(i) + rkphysa(i,j)
          bb(i,j) = 0.0
        enddo
      enddo

      ratio = dtn/dto
      do L = 2,nrk
        theta = 1.0 + ratio*c(L)
        do LL = 1,nrk
        bb(L,LL) = bD(LL,1)*theta 
     &           + bD(LL,2)*theta*theta 
     &           + bD(LL,3)*theta*theta*theta
        enddo
      enddo


c   M  : variable, (rho, rho*u . . . . )
c   id : zone
c   i  : 3-D array of data

      do 200 nv = 1,5
        do 100 k = kmn+1,kmx
        do 100 j = jmn+1,jmx
        do 100 i = imn+1,imx
        voli           = 1./vol(i,j,k)
        voli           = 1.

            do irk = 2,nrk
              wex(irk) = wold(i,j,k,nv)
              do lrk = 1,nrk
                wex(irk) = wex(irk) 
     .                   - dto*bb(irk,lrk)*delwold(i,j,k,nv,lrk)*voli
cvn  .                   + dto*bb(irk,lrk)*delwold(i,j,k,nv,lrk)*voli
              enddo
            enddo

c
c      next line could be used for saving flux computation
c      on top of a new stage
c           delwold(i,j,k,nv,1) = delwold(i,j,k,nv,nrk)
c
            do irk = 2,nrk
c      save the extrapolated value of conserved variables in delwold
              delwold(i,j,k,nv,irk) = wex(irk) 
            enddo

 100    continue
 200  continue
      return
#else
      return
#endif
      end
