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
      subroutine flowmb(imn,jmn,kmn,
     .                  imx,jmx,kmx,imp1,jmp1,kmp1,imp2,jmp2,kmp2,
     .                  mbloc,msegt,msegn,nseg,imap,ibloc,
     .                  w,p,x,wn,alpha,totptim,ifirst               )
c***********************************************************************
c
c     object         : To move flow quantities from cell-centers to nodes
c     initial coding : v.n.vatsa (may, 1991)
c     latest mods    : remove iblank array (May, 2000)
c
      dimension w (imp2,jmp2,kmp2,5),  p(imp2,jmp2,kmp2),
     .          wn(imp2,jmp2,kmp2,5),  x(imp1,jmp1,kmp1,3)
c
      dimension  imap(msegt,msegn,mbloc), nseg(mbloc)
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
c
c     evaluate convenient groupings
c
      eighti=  1.000/8.000
      scalv = eighti/sqrt(gamma)
      scalp = eighti/gamma
c
      do 20 k=kmn,kmx
      kp = k+1
      do 20 j=jmn,jmx
      do 20 i=imn,imx
      wn(i,j,k,1) = eighti*(    w(i+1,j+1,k ,1)   +w(i,j+1,k ,1)
     .                         +w(i+1,j  ,k ,1)   +w(i,j  ,k ,1)
     .                         +w(i+1,j+1,kp,1)   +w(i,j+1,kp,1)
     .                         +w(i+1,j,  kp,1)   +w(i,j  ,kp,1) )
      wn(i,j,k,2) = scalv *(    w(i+1,j+1,k ,2)   +w(i,j+1,k ,2)
     .                         +w(i+1,j  ,k ,2)   +w(i,j  ,k ,2)
     .                         +w(i+1,j+1,kp,2)   +w(i,j+1,kp,2)
     .                         +w(i+1,j  ,kp,2)   +w(i,j  ,kp,2) )
      wn(i,j,k,3) = scalv *(    w(i+1,j+1,k ,3)   +w(i,j+1,k ,3)
     .                         +w(i+1,j  ,k ,3)   +w(i,j  ,k ,3)
     .                         +w(i+1,j+1,kp,3)   +w(i,j+1,kp,3)
     .                         +w(i+1,j  ,kp,3)   +w(i,j  ,kp,3) )
      wn(i,j,k,4) = scalv *(    w(i+1,j+1,k ,4)   +w(i,j+1,k ,4)
     .                         +w(i+1,j  ,k ,4)   +w(i,j  ,k ,4)
     .                         +w(i+1,j+1,kp,4)   +w(i,j+1,kp,4)
     .                         +w(i+1,j  ,kp,4)   +w(i,j  ,kp,4) )
      pn          = scalp *(    p(i+1,j+1,k )     +p(i,j+1,k )
     .                         +p(i+1,j  ,k )     +p(i,j  ,k )
     .                         +p(i+1,j+1,kp)     +p(i,j+1,kp)
     .                         +p(i+1,j  ,kp)     +p(i,j  ,kp) )
      wn(i,j,k,5) = dgm1*pn +.500/wn(i,j,k,1)*
     .                          ( wn(i,j,k,2)*wn(i,j,k,2) +
     .                            wn(i,j,k,3)*wn(i,j,k,3) +
     .                            wn(i,j,k,4)*wn(i,j,k,4) )
  20  continue
c
c     write the solution vector on ipltq
c
      write(iqpfil) rm,alpha,rey,totptim
      write(iqpfil)((((wn(i,j,k,nv),i=imn,imx),j=jmn,jmx),k=kmn,kmx),
     .                                                      nv=1,5 )
c
c     return if not first time step
c
      if(ifirst .ne. 1) return
c
c
c       write the coordinates of reduced mesh size on file 'ipltg'
c		
        write (igpfil)
     .  (((x(i,j,k,1),i=imn,imx),j=jmn,jmx),k=kmn,kmx),
     .  (((x(i,j,k,2),i=imn,imx),j=jmn,jmx),k=kmn,kmx),
     .  (((x(i,j,k,3),i=imn,imx),j=jmn,jmx),k=kmn,kmx)
c
c
      return
      end
