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
      subroutine  bcedgfl(imn,jmn,kmn,imx,jmx,kmx,imp1,jmp1,kmp1,
     .                    imp2,jmp2,kmp2,w,p,eomu,igrid,isoln   )
c=======================================================================
c
c
c     object : to fill-in the flow variables at block corners
c              using extrapolation
c
c     initial coding : by v.n.vatsa, nasa langley ( Oct., 1994)
c
c

      dimension      w(imp2,jmp2,kmp2,5),
     .               p(imp2,jmp2,kmp2)  ,  eomu(imp2,jmp2,kmp2)
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      imnp1   = imn+1
      jmnp1   = jmn+1
      kmnp1   = kmn+1
c
c     mods for filling corner points at i-faces
c     By : vatsa (Jan. 1997)
c
c*******************************************************************
c     fill in jmn,kmn corner

      k  = kmn
      j  = jmn

      do i=imnp1,imx
        w(i,j,k,1)   = .5*(w(i,j+1,k,1)+w(i,j,k+1,1))
        w(i,j,k,2)   = .5*(w(i,j+1,k,2)+w(i,j,k+1,2))
        w(i,j,k,3)   = .5*(w(i,j+1,k,3)+w(i,j,k+1,3))
        w(i,j,k,4)   = .5*(w(i,j+1,k,4)+w(i,j,k+1,4))
        w(i,j,k,5)   = .5*(w(i,j+1,k,5)+w(i,j,k+1,5))
        p(i,j,k  )   = .5*(p(i,j+1,k  )+p(i,j,k+1  ))
        eomu(i,j,k)  = .5*(eomu(i,j+1,k)+eomu(i,j,k+1))
      enddo

c
      if(igrid.eq.isoln) then
c
      do i=imnp1,imx
      w(i,j-1,k-1,1)   = w(i,j,k,1)
      w(i,j-1,k-1,2)   = w(i,j,k,2)
      w(i,j-1,k-1,3)   = w(i,j,k,3)
      w(i,j-1,k-1,4)   = w(i,j,k,4)
      w(i,j-1,k-1,5)   = w(i,j,k,5)
      p(i,j-1,k-1  )   = p(i,j,k)
      eomu(i,j-1,k-1)  = eomu(i,j,k)
c
      w(i,j  ,k-1,1)   = w(i,j,k,1)
      w(i,j  ,k-1,2)   = w(i,j,k,2)
      w(i,j  ,k-1,3)   = w(i,j,k,3)
      w(i,j  ,k-1,4)   = w(i,j,k,4)
      w(i,j  ,k-1,5)   = w(i,j,k,5)
      p(i,j  ,k-1  )   = p(i,j,k)
      eomu(i,j  ,k-1)  = eomu(i,j,k)
c
      w(i,j-1,k  ,1)   = w(i,j,k,1)
      w(i,j-1,k  ,2)   = w(i,j,k,2)
      w(i,j-1,k  ,3)   = w(i,j,k,3)
      w(i,j-1,k  ,4)   = w(i,j,k,4)
      w(i,j-1,k  ,5)   = w(i,j,k,5)
      p(i,j-1,k    )   = p(i,j,k)
      eomu(i,j-1,k  )  = eomu(i,j,k)
c
      enddo
     
      endif
c
c
c     fill in jmp1,kmn corner

      k  = kmn
      j  = jmp1

      do i=imnp1,imx
        w(i,j,k,1)   = .5*(w(i,j-1,k,1)+w(i,j,k+1,1))
        w(i,j,k,2)   = .5*(w(i,j-1,k,2)+w(i,j,k+1,2))
        w(i,j,k,3)   = .5*(w(i,j-1,k,3)+w(i,j,k+1,3))
        w(i,j,k,4)   = .5*(w(i,j-1,k,4)+w(i,j,k+1,4))
        w(i,j,k,5)   = .5*(w(i,j-1,k,5)+w(i,j,k+1,5))
        p(i,j,k  )   = .5*(p(i,j-1,k  )+p(i,j,k+1  ))
        eomu(i,j,k)  = .5*(eomu(i,j-1,k)+eomu(i,j,k+1))
      enddo
c
c
      if(igrid.eq.isoln) then
      do i=imnp1,imx
c
      w(i,j+1,k-1,1)   = w(i,j,k,1)
      w(i,j+1,k-1,2)   = w(i,j,k,2)
      w(i,j+1,k-1,3)   = w(i,j,k,3)
      w(i,j+1,k-1,4)   = w(i,j,k,4)
      w(i,j+1,k-1,5)   = w(i,j,k,5)
      p(i,j+1,k-1  )   = p(i,j,k)
      eomu(i,j+1,k-1)  = eomu(i,j,k)
c
      w(i,j  ,k-1,1)   = w(i,j,k,1)
      w(i,j  ,k-1,2)   = w(i,j,k,2)
      w(i,j  ,k-1,3)   = w(i,j,k,3)
      w(i,j  ,k-1,4)   = w(i,j,k,4)
      w(i,j  ,k-1,5)   = w(i,j,k,5)
      p(i,j  ,k-1  )   = p(i,j,k)
      eomu(i,j  ,k-1)  = eomu(i,j,k)
c
      w(i,j+1,k  ,1)   = w(i,j,k,1)
      w(i,j+1,k  ,2)   = w(i,j,k,2)
      w(i,j+1,k  ,3)   = w(i,j,k,3)
      w(i,j+1,k  ,4)   = w(i,j,k,4)
      w(i,j+1,k  ,5)   = w(i,j,k,5)
      p(i,j+1,k    )   = p(i,j,k)
      eomu(i,j+1,k  )  = eomu(i,j,k)
c
      enddo
     
      endif
c
c
c     fill in jmn,kmp1 corner

      k  = kmp1
      j  = jmn

      do i=imnp1,imx
        w(i,j,k,1)   = .5*(w(i,j+1,k,1)+w(i,j,k-1,1))
        w(i,j,k,2)   = .5*(w(i,j+1,k,2)+w(i,j,k-1,2))
        w(i,j,k,3)   = .5*(w(i,j+1,k,3)+w(i,j,k-1,3))
        w(i,j,k,4)   = .5*(w(i,j+1,k,4)+w(i,j,k-1,4))
        w(i,j,k,5)   = .5*(w(i,j+1,k,5)+w(i,j,k-1,5))
        p(i,j,k  )   = .5*(p(i,j+1,k  )+p(i,j,k-1  ))
        eomu(i,j,k)  = .5*(eomu(i,j+1,k)+eomu(i,j,k-1))
      enddo
c
c
      if(igrid.eq.isoln) then
c
      do i=imnp1,imx
      w(i,j-1,k+1,1)   = w(i,j,k,1)
      w(i,j-1,k+1,2)   = w(i,j,k,2)
      w(i,j-1,k+1,3)   = w(i,j,k,3)
      w(i,j-1,k+1,4)   = w(i,j,k,4)
      w(i,j-1,k+1,5)   = w(i,j,k,5)
      p(i,j-1,k+1  )   = p(i,j,k)
      eomu(i,j-1,k+1)  = eomu(i,j,k)
c
      w(i,j  ,k+1,1)   = w(i,j,k,1)
      w(i,j  ,k+1,2)   = w(i,j,k,2)
      w(i,j  ,k+1,3)   = w(i,j,k,3)
      w(i,j  ,k+1,4)   = w(i,j,k,4)
      w(i,j  ,k+1,5)   = w(i,j,k,5)
      p(i,j  ,k+1  )   = p(i,j,k)
      eomu(i,j  ,k+1)  = eomu(i,j,k)
c
      w(i,j-1,k  ,1)   = w(i,j,k,1)
      w(i,j-1,k  ,2)   = w(i,j,k,2)
      w(i,j-1,k  ,3)   = w(i,j,k,3)
      w(i,j-1,k  ,4)   = w(i,j,k,4)
      w(i,j-1,k  ,5)   = w(i,j,k,5)
      p(i,j-1,k    )   = p(i,j,k)
      eomu(i,j-1,k  )  = eomu(i,j,k)
      enddo

      endif
c
c
c     fill in jmp1,kmp1 corner

      k  = kmp1
      j  = jmp1

      do i=imnp1,imx
        w(i,j,k,1)   = .5*(w(i,j-1,k,1)+w(i,j,k-1,1))
        w(i,j,k,2)   = .5*(w(i,j-1,k,2)+w(i,j,k-1,2))
        w(i,j,k,3)   = .5*(w(i,j-1,k,3)+w(i,j,k-1,3))
        w(i,j,k,4)   = .5*(w(i,j-1,k,4)+w(i,j,k-1,4))
        w(i,j,k,5)   = .5*(w(i,j-1,k,5)+w(i,j,k-1,5))
        p(i,j,k  )   = .5*(p(i,j-1,k  )+p(i,j,k-1  ))
        eomu(i,j,k)  = .5*(eomu(i,j-1,k)+eomu(i,j,k-1))
      enddo
c
      if(igrid.eq.isoln) then
c
      do i=imnp1,imx
      w(i,j+1,k+1,1)   = w(i,j,k,1)
      w(i,j+1,k+1,2)   = w(i,j,k,2)
      w(i,j+1,k+1,3)   = w(i,j,k,3)
      w(i,j+1,k+1,4)   = w(i,j,k,4)
      w(i,j+1,k+1,5)   = w(i,j,k,5)
      p(i,j+1,k+1  )   = p(i,j,k)
      eomu(i,j+1,k+1)  = eomu(i,j,k)
c
      w(i,j  ,k+1,1)   = w(i,j,k,1)
      w(i,j  ,k+1,2)   = w(i,j,k,2)
      w(i,j  ,k+1,3)   = w(i,j,k,3)
      w(i,j  ,k+1,4)   = w(i,j,k,4)
      w(i,j  ,k+1,5)   = w(i,j,k,5)
      p(i,j  ,k+1  )   = p(i,j,k)
      eomu(i,j  ,k+1)  = eomu(i,j,k)
c
      w(i,j+1,k  ,1)   = w(i,j,k,1)
      w(i,j+1,k  ,2)   = w(i,j,k,2)
      w(i,j+1,k  ,3)   = w(i,j,k,3)
      w(i,j+1,k  ,4)   = w(i,j,k,4)
      w(i,j+1,k  ,5)   = w(i,j,k,5)
      p(i,j+1,k    )   = p(i,j,k)
      eomu(i,j+1,k  )  = eomu(i,j,k)
      enddo
c
      endif
c
c*******************************************************************

      do 100 k=kmn,kmp1,kmp1-kmn
c
c$dir novector
      do 110 j=jmn,jmp1,jmp1-jmn
c
      w(imn ,j,k,1)   = w(imnp1,j,k,1)
      w(imn ,j,k,2)   = w(imnp1,j,k,2)
      w(imn ,j,k,3)   = w(imnp1,j,k,3)
      w(imn ,j,k,4)   = w(imnp1,j,k,4)
      w(imn ,j,k,5)   = w(imnp1,j,k,5)
      p(imn ,j,k)     = p(imnp1,j,k)
      eomu(imn  ,j,k) = eomu(imnp1,j,k)
c
      w(imp1,j,k,1)   = w(imx  ,j,k,1)
      w(imp1,j,k,2)   = w(imx  ,j,k,2)
      w(imp1,j,k,3)   = w(imx  ,j,k,3)
      w(imp1,j,k,4)   = w(imx  ,j,k,4)
      w(imp1,j,k,5)   = w(imx  ,j,k,5)
      p(imp1,j,k)     = p(imx  ,j,k)
      eomu(imp1,j,k)  = eomu(imx  ,j,k)
 110  continue
c
 100  continue
c
      if(igrid.ne.isoln) return
c
c
      imnm1   = imn-1
      jmnm1   = jmn-1
      kmnm1   = kmn-1
c
      do 120 k=kmnm1,kmp2,kmp2-kmnm1
c
c$dir novector
      do 130 j=jmnm1,jmp2,jmp2-jmnm1
c
      w(imn ,j,k,1)   = w(imnp1,j,k,1)
      w(imn ,j,k,2)   = w(imnp1,j,k,2)
      w(imn ,j,k,3)   = w(imnp1,j,k,3)
      w(imn ,j,k,4)   = w(imnp1,j,k,4)
      w(imn ,j,k,5)   = w(imnp1,j,k,5)
      p(imn ,j,k)     = p(imnp1,j,k)
      eomu(imn  ,j,k) = eomu(imnp1,j,k)
c
      w(imp1,j,k,1)   = w(imx  ,j,k,1)
      w(imp1,j,k,2)   = w(imx  ,j,k,2)
      w(imp1,j,k,3)   = w(imx  ,j,k,3)
      w(imp1,j,k,4)   = w(imx  ,j,k,4)
      w(imp1,j,k,5)   = w(imx  ,j,k,5)
      p(imp1,j,k)     = p(imx  ,j,k)
      eomu(imp1,j,k)  = eomu(imx  ,j,k)
 130  continue
c
 120  continue
c
      do 200 k=kmnm1,kmp2,kmp2-kmnm1
c
c$dir novector
      do 210 j=jmnm1,jmp2,jmp2-jmnm1
c
      w(imnm1,j,k,1)   = w(imn,j,k,1)
      w(imnm1,j,k,2)   = w(imn,j,k,2)
      w(imnm1,j,k,3)   = w(imn,j,k,3)
      w(imnm1,j,k,4)   = w(imn,j,k,4)
      w(imnm1,j,k,5)   = w(imn,j,k,5)
      p(imnm1,j,k)     = p(imn,j,k)
      eomu(imnm1 ,j,k) = eomu(imnp1,j,k)
c
      w(imp2,j,k,1)    = w(imp1,j,k,1)
      w(imp2,j,k,2)    = w(imp1,j,k,2)
      w(imp2,j,k,3)    = w(imp1,j,k,3)
      w(imp2,j,k,4)    = w(imp1,j,k,4)
      w(imp2,j,k,5)    = w(imp1,j,k,5)
      p(imp2,j,k)      = p(imp1,j,k)
      eomu(imp2,j,k)   = eomu(imp1,j,k)
 210  continue
c
 200  continue
c
c     fill the remaining corner points in second layer
c
c
      j = jmn
      k = kmn
c
      w(imn,j  ,k-1,1)   = w(imn,j,k,1)
      w(imn,j  ,k-1,2)   = w(imn,j,k,2)
      w(imn,j  ,k-1,3)   = w(imn,j,k,3)
      w(imn,j  ,k-1,4)   = w(imn,j,k,4)
      w(imn,j  ,k-1,5)   = w(imn,j,k,5)
      p(imn,j  ,k-1  )   = p(imn,j,k)
      eomu(imn,j  ,k-1)  = eomu(imn,j,k)
c
      w(imn,j-1,k  ,1)   = w(imn,j,k,1)
      w(imn,j-1,k  ,2)   = w(imn,j,k,2)
      w(imn,j-1,k  ,3)   = w(imn,j,k,3)
      w(imn,j-1,k  ,4)   = w(imn,j,k,4)
      w(imn,j-1,k  ,5)   = w(imn,j,k,5)
      p(imn,j-1,k    )   = p(imn,j,k)
      eomu(imn,j-1,k  )  = eomu(imn,j,k)
c
      w(imp1  ,j  ,k-1,1)   = w(imp1  ,j,k,1)
      w(imp1  ,j  ,k-1,2)   = w(imp1  ,j,k,2)
      w(imp1  ,j  ,k-1,3)   = w(imp1  ,j,k,3)
      w(imp1  ,j  ,k-1,4)   = w(imp1  ,j,k,4)
      w(imp1  ,j  ,k-1,5)   = w(imp1  ,j,k,5)
      p(imp1  ,j  ,k-1  )   = p(imp1  ,j,k  )
      eomu(imp1  ,j  ,k-1)  = eomu(imp1  ,j,k)
c
      w(imp1  ,j-1,k  ,1)   = w(imp1  ,j,k,1)
      w(imp1  ,j-1,k  ,2)   = w(imp1  ,j,k,2)
      w(imp1  ,j-1,k  ,3)   = w(imp1  ,j,k,3)
      w(imp1  ,j-1,k  ,4)   = w(imp1  ,j,k,4)
      w(imp1  ,j-1,k  ,5)   = w(imp1  ,j,k,5)
      p(imp1  ,j-1,k    )   = p(imp1  ,j,k  )
      eomu(imp1  ,j-1,k  )  = eomu(imp1  ,j,k)
c
      w(imnm1,j  ,k-1,1)   = w(imnm1,j,k,1)
      w(imnm1,j  ,k-1,2)   = w(imnm1,j,k,2)
      w(imnm1,j  ,k-1,3)   = w(imnm1,j,k,3)
      w(imnm1,j  ,k-1,4)   = w(imnm1,j,k,4)
      w(imnm1,j  ,k-1,5)   = w(imnm1,j,k,5)
      p(imnm1,j  ,k-1  )   = p(imnm1,j,k)
      eomu(imnm1,j  ,k-1)  = eomu(imnm1,j,k)
c
      w(imnm1,j-1,k  ,1)   = w(imnm1,j,k,1)
      w(imnm1,j-1,k  ,2)   = w(imnm1,j,k,2)
      w(imnm1,j-1,k  ,3)   = w(imnm1,j,k,3)
      w(imnm1,j-1,k  ,4)   = w(imnm1,j,k,4)
      w(imnm1,j-1,k  ,5)   = w(imnm1,j,k,5)
      p(imnm1,j-1,k    )   = p(imnm1,j,k)
      eomu(imnm1,j-1,k  )  = eomu(imnm1,j,k)
c
      w(imp2  ,j  ,k-1,1)   = w(imp2  ,j,k,1)
      w(imp2  ,j  ,k-1,2)   = w(imp2  ,j,k,2)
      w(imp2  ,j  ,k-1,3)   = w(imp2  ,j,k,3)
      w(imp2  ,j  ,k-1,4)   = w(imp2  ,j,k,4)
      w(imp2  ,j  ,k-1,5)   = w(imp2  ,j,k,5)
      p(imp2  ,j  ,k-1  )   = p(imp2  ,j,k  )
      eomu(imp2  ,j  ,k-1)  = eomu(imp2  ,j,k)
c
      w(imp2  ,j-1,k  ,1)   = w(imp2  ,j,k,1)
      w(imp2  ,j-1,k  ,2)   = w(imp2  ,j,k,2)
      w(imp2  ,j-1,k  ,3)   = w(imp2  ,j,k,3)
      w(imp2  ,j-1,k  ,4)   = w(imp2  ,j,k,4)
      w(imp2  ,j-1,k  ,5)   = w(imp2  ,j,k,5)
      p(imp2  ,j-1,k    )   = p(imp2  ,j,k  )
      eomu(imp2  ,j-1,k  )  = eomu(imp2  ,j,k)
c
      j  =  jmp1
      k  =  kmn
c
      w(imn,j  ,k-1,1)   = w(imn,j,k,1)
      w(imn,j  ,k-1,2)   = w(imn,j,k,2)
      w(imn,j  ,k-1,3)   = w(imn,j,k,3)
      w(imn,j  ,k-1,4)   = w(imn,j,k,4)
      w(imn,j  ,k-1,5)   = w(imn,j,k,5)
      p(imn,j  ,k-1  )   = p(imn,j,k)
      eomu(imn,j  ,k-1)  = eomu(imn,j,k)
c
      w(imn,j+1,k  ,1)   = w(imn,j,k,1)
      w(imn,j+1,k  ,2)   = w(imn,j,k,2)
      w(imn,j+1,k  ,3)   = w(imn,j,k,3)
      w(imn,j+1,k  ,4)   = w(imn,j,k,4)
      w(imn,j+1,k  ,5)   = w(imn,j,k,5)
      p(imn,j+1,k    )   = p(imn,j,k)
      eomu(imn,j+1,k  )  = eomu(imn,j,k)
c
      w(imnm1,j  ,k-1,1)   = w(imnm1,j,k,1)
      w(imnm1,j  ,k-1,2)   = w(imnm1,j,k,2)
      w(imnm1,j  ,k-1,3)   = w(imnm1,j,k,3)
      w(imnm1,j  ,k-1,4)   = w(imnm1,j,k,4)
      w(imnm1,j  ,k-1,5)   = w(imnm1,j,k,5)
      p(imnm1,j  ,k-1  )   = p(imnm1,j,k)
      eomu(imnm1,j  ,k-1)  = eomu(imnm1,j,k)
c
      w(imnm1,j+1,k  ,1)   = w(imnm1,j,k,1)
      w(imnm1,j+1,k  ,2)   = w(imnm1,j,k,2)
      w(imnm1,j+1,k  ,3)   = w(imnm1,j,k,3)
      w(imnm1,j+1,k  ,4)   = w(imnm1,j,k,4)
      w(imnm1,j+1,k  ,5)   = w(imnm1,j,k,5)
      p(imnm1,j+1,k    )   = p(imnm1,j,k)
      eomu(imnm1,j+1,k  )  = eomu(imnm1,j,k)
c
      w(imp1  ,j  ,k-1,1)   = w(imp1  ,j,k,1)
      w(imp1  ,j  ,k-1,2)   = w(imp1  ,j,k,2)
      w(imp1  ,j  ,k-1,3)   = w(imp1  ,j,k,3)
      w(imp1  ,j  ,k-1,4)   = w(imp1  ,j,k,4)
      w(imp1  ,j  ,k-1,5)   = w(imp1  ,j,k,5)
      p(imp1  ,j  ,k-1  )   = p(imp1  ,j,k)
      eomu(imp1  ,j  ,k-1)  = eomu(imp1  ,j,k)
c
      w(imp1  ,j+1,k  ,1)   = w(imp1  ,j,k,1)
      w(imp1  ,j+1,k  ,2)   = w(imp1  ,j,k,2)
      w(imp1  ,j+1,k  ,3)   = w(imp1  ,j,k,3)
      w(imp1  ,j+1,k  ,4)   = w(imp1  ,j,k,4)
      w(imp1  ,j+1,k  ,5)   = w(imp1  ,j,k,5)
      p(imp1  ,j+1,k    )   = p(imp1  ,j,k)
      eomu(imp1  ,j+1,k  )  = eomu(imp1  ,j,k)
c
      w(imp2  ,j  ,k-1,1)   = w(imp2  ,j,k,1)
      w(imp2  ,j  ,k-1,2)   = w(imp2  ,j,k,2)
      w(imp2  ,j  ,k-1,3)   = w(imp2  ,j,k,3)
      w(imp2  ,j  ,k-1,4)   = w(imp2  ,j,k,4)
      w(imp2  ,j  ,k-1,5)   = w(imp2  ,j,k,5)
      p(imp2  ,j  ,k-1  )   = p(imp2  ,j,k)
      eomu(imp2  ,j  ,k-1)  = eomu(imp2  ,j,k)
c
      w(imp2  ,j+1,k  ,1)   = w(imp2  ,j,k,1)
      w(imp2  ,j+1,k  ,2)   = w(imp2  ,j,k,2)
      w(imp2  ,j+1,k  ,3)   = w(imp2  ,j,k,3)
      w(imp2  ,j+1,k  ,4)   = w(imp2  ,j,k,4)
      w(imp2  ,j+1,k  ,5)   = w(imp2  ,j,k,5)
      p(imp2  ,j+1,k    )   = p(imp2  ,j,k)
      eomu(imp2  ,j+1,k  )  = eomu(imp2  ,j,k)
c
      j  = jmn
      k  = kmp1
c
      w(imn,j  ,k+1,1)   = w(imn,j,k,1)
      w(imn,j  ,k+1,2)   = w(imn,j,k,2)
      w(imn,j  ,k+1,3)   = w(imn,j,k,3)
      w(imn,j  ,k+1,4)   = w(imn,j,k,4)
      w(imn,j  ,k+1,5)   = w(imn,j,k,5)
      p(imn,j  ,k+1  )   = p(imn,j,k)
      eomu(imn,j  ,k+1)  = eomu(imn,j,k)
c
      w(imn,j-1,k  ,1)   = w(imn,j,k,1)
      w(imn,j-1,k  ,2)   = w(imn,j,k,2)
      w(imn,j-1,k  ,3)   = w(imn,j,k,3)
      w(imn,j-1,k  ,4)   = w(imn,j,k,4)
      w(imn,j-1,k  ,5)   = w(imn,j,k,5)
      p(imn,j-1,k    )   = p(imn,j,k)
      eomu(imn,j-1,k  )  = eomu(imn,j,k)
c
      w(imnm1,j  ,k+1,1)   = w(imnm1,j,k,1)
      w(imnm1,j  ,k+1,2)   = w(imnm1,j,k,2)
      w(imnm1,j  ,k+1,3)   = w(imnm1,j,k,3)
      w(imnm1,j  ,k+1,4)   = w(imnm1,j,k,4)
      w(imnm1,j  ,k+1,5)   = w(imnm1,j,k,5)
      p(imnm1,j  ,k+1  )   = p(imnm1,j,k)
      eomu(imnm1,j  ,k+1)  = eomu(imnm1,j,k)
c
      w(imnm1,j-1,k  ,1)   = w(imnm1,j,k,1)
      w(imnm1,j-1,k  ,2)   = w(imnm1,j,k,2)
      w(imnm1,j-1,k  ,3)   = w(imnm1,j,k,3)
      w(imnm1,j-1,k  ,4)   = w(imnm1,j,k,4)
      w(imnm1,j-1,k  ,5)   = w(imnm1,j,k,5)
      p(imnm1,j-1,k    )   = p(imnm1,j,k)
      eomu(imnm1,j-1,k  )  = eomu(imnm1,j,k)
c
      w(imp1  ,j  ,k+1,1)   = w(imp1  ,j,k,1)
      w(imp1  ,j  ,k+1,2)   = w(imp1  ,j,k,2)
      w(imp1  ,j  ,k+1,3)   = w(imp1  ,j,k,3)
      w(imp1  ,j  ,k+1,4)   = w(imp1  ,j,k,4)
      w(imp1  ,j  ,k+1,5)   = w(imp1  ,j,k,5)
      p(imp1  ,j  ,k+1  )   = p(imp1  ,j,k  )
      eomu(imp1  ,j  ,k+1)  = eomu(imp1  ,j,k)
c
      w(imp1  ,j-1,k  ,1)   = w(imp1  ,j,k,1)
      w(imp1  ,j-1,k  ,2)   = w(imp1  ,j,k,2)
      w(imp1  ,j-1,k  ,3)   = w(imp1  ,j,k,3)
      w(imp1  ,j-1,k  ,4)   = w(imp1  ,j,k,4)
      w(imp1  ,j-1,k  ,5)   = w(imp1  ,j,k,5)
      p(imp1  ,j-1,k    )   = p(imp1  ,j,k  )
      eomu(imp1  ,j-1,k  )  = eomu(imp1  ,j,k)
c
      w(imp2  ,j  ,k+1,1)   = w(imp2  ,j,k,1)
      w(imp2  ,j  ,k+1,2)   = w(imp2  ,j,k,2)
      w(imp2  ,j  ,k+1,3)   = w(imp2  ,j,k,3)
      w(imp2  ,j  ,k+1,4)   = w(imp2  ,j,k,4)
      w(imp2  ,j  ,k+1,5)   = w(imp2  ,j,k,5)
      p(imp2  ,j  ,k+1  )   = p(imp2  ,j,k  )
      eomu(imp2  ,j  ,k+1)  = eomu(imp2  ,j,k)
c
      w(imp2  ,j-1,k  ,1)   = w(imp2  ,j,k,1)
      w(imp2  ,j-1,k  ,2)   = w(imp2  ,j,k,2)
      w(imp2  ,j-1,k  ,3)   = w(imp2  ,j,k,3)
      w(imp2  ,j-1,k  ,4)   = w(imp2  ,j,k,4)
      w(imp2  ,j-1,k  ,5)   = w(imp2  ,j,k,5)
      p(imp2  ,j-1,k    )   = p(imp2  ,j,k  )
      eomu(imp2  ,j-1,k  )  = eomu(imp2  ,j,k)
c
      j  = jmp1
      k  = kmp1
c
      w(imn,j  ,k+1,1)   = w(imn,j,k,1)
      w(imn,j  ,k+1,2)   = w(imn,j,k,2)
      w(imn,j  ,k+1,3)   = w(imn,j,k,3)
      w(imn,j  ,k+1,4)   = w(imn,j,k,4)
      w(imn,j  ,k+1,5)   = w(imn,j,k,5)
      p(imn,j  ,k+1  )   = p(imn,j,k)
      eomu(imn,j  ,k+1)  = eomu(imn,j,k)
c
      w(imn,j+1,k  ,1)   = w(imn,j,k,1)
      w(imn,j+1,k  ,2)   = w(imn,j,k,2)
      w(imn,j+1,k  ,3)   = w(imn,j,k,3)
      w(imn,j+1,k  ,4)   = w(imn,j,k,4)
      w(imn,j+1,k  ,5)   = w(imn,j,k,5)
      p(imn,j+1,k    )   = p(imn,j,k)
      eomu(imn,j+1,k  )  = eomu(imn,j,k)
c
      w(imnm1,j  ,k+1,1)   = w(imnm1,j,k,1)
      w(imnm1,j  ,k+1,2)   = w(imnm1,j,k,2)
      w(imnm1,j  ,k+1,3)   = w(imnm1,j,k,3)
      w(imnm1,j  ,k+1,4)   = w(imnm1,j,k,4)
      w(imnm1,j  ,k+1,5)   = w(imnm1,j,k,5)
      p(imnm1,j  ,k+1  )   = p(imnm1,j,k)
      eomu(imnm1,j  ,k+1)  = eomu(imnm1,j,k)
c
      w(imnm1,j+1,k  ,1)   = w(imnm1,j,k,1)
      w(imnm1,j+1,k  ,2)   = w(imnm1,j,k,2)
      w(imnm1,j+1,k  ,3)   = w(imnm1,j,k,3)
      w(imnm1,j+1,k  ,4)   = w(imnm1,j,k,4)
      w(imnm1,j+1,k  ,5)   = w(imnm1,j,k,5)
      p(imnm1,j+1,k    )   = p(imnm1,j,k)
      eomu(imnm1,j+1,k  )  = eomu(imnm1,j,k)
c
      w(imp1  ,j  ,k+1,1)   = w(imp1  ,j,k,1)
      w(imp1  ,j  ,k+1,2)   = w(imp1  ,j,k,2)
      w(imp1  ,j  ,k+1,3)   = w(imp1  ,j,k,3)
      w(imp1  ,j  ,k+1,4)   = w(imp1  ,j,k,4)
      w(imp1  ,j  ,k+1,5)   = w(imp1  ,j,k,5)
      p(imp1  ,j  ,k+1  )   = p(imp1  ,j,k  )
      eomu(imp1  ,j  ,k+1)  = eomu(imp1  ,j,k)
c
      w(imp1  ,j+1,k  ,1)   = w(imp1  ,j,k,1)
      w(imp1  ,j+1,k  ,2)   = w(imp1  ,j,k,2)
      w(imp1  ,j+1,k  ,3)   = w(imp1  ,j,k,3)
      w(imp1  ,j+1,k  ,4)   = w(imp1  ,j,k,4)
      w(imp1  ,j+1,k  ,5)   = w(imp1  ,j,k,5)
      p(imp1  ,j+1,k    )   = p(imp1  ,j,k  )
      eomu(imp1  ,j+1,k  )  = eomu(imp1  ,j,k)
c
      w(imp2  ,j  ,k+1,1)   = w(imp2  ,j,k,1)
      w(imp2  ,j  ,k+1,2)   = w(imp2  ,j,k,2)
      w(imp2  ,j  ,k+1,3)   = w(imp2  ,j,k,3)
      w(imp2  ,j  ,k+1,4)   = w(imp2  ,j,k,4)
      w(imp2  ,j  ,k+1,5)   = w(imp2  ,j,k,5)
      p(imp2  ,j  ,k+1  )   = p(imp2  ,j,k  )
      eomu(imp2  ,j  ,k+1)  = eomu(imp2  ,j,k)
c
      w(imp2  ,j+1,k  ,1)   = w(imp2  ,j,k,1)
      w(imp2  ,j+1,k  ,2)   = w(imp2  ,j,k,2)
      w(imp2  ,j+1,k  ,3)   = w(imp2  ,j,k,3)
      w(imp2  ,j+1,k  ,4)   = w(imp2  ,j,k,4)
      w(imp2  ,j+1,k  ,5)   = w(imp2  ,j,k,5)
      p(imp2  ,j+1,k    )   = p(imp2  ,j,k  )
      eomu(imp2  ,j+1,k  )  = eomu(imp2  ,j,k)
c
      return
      end
