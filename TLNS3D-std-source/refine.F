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
      subroutine refine( imn,jmn,kmn,imp2f,jmp2f,kmp2f,
     .                   imp2c,jmp2c,kmp2c,
     .                   wf,wc,pf,pc,eomuf,eomuc,turv1f,turv1c,
     .                                           turv2f,turv2c)
c
c

      common/trbfs/ anutinf,zkinf,ominf
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension wf(imp2f,jmp2f,kmp2f,5),  wc(imp2c,jmp2c,kmp2c,5),
     .          pf(imp2f,jmp2f,kmp2f),    pc(imp2c,jmp2c,kmp2c),
     .       eomuf(imp2f,jmp2f,kmp2f), eomuc(imp2c,jmp2c,kmp2c),
     .      turv1f(imp2f,jmp2f,kmp2f),turv1c(imp2c,jmp2c,kmp2c),
     .      turv2f(imp2f,jmp2f,kmp2f),turv2c(imp2c,jmp2c,kmp2c)
c
c
      kmnp1 = kmn   + 1
      jmnp1 = jmn   + 1
      imnp1 = imn   + 1
c
      kmp1f = kmp2f - 1
      jmp1f = jmp2f - 1
      imp1f = imp2f - 1
      kmp1c = kmp2c - 1
      jmp1c = jmp2c - 1
      imp1c = imp2c - 1
c
      kmf   = kmp2f - 2
      jmf   = jmp2f - 2
      imf   = imp2f - 2
      kmc   = kmp2c - 2
      jmc   = jmp2c - 2
      imc   = imp2c - 2
c
      sfi=1./64.
c
      do 10 n=1,5
      kc        = kmn -1
      do 10 k=kmn,kmp2f,2
      kc        = kc  +1
      jc        = jmn -1
      do 10 j=jmn,jmp2f,2
      jc        = jc  +1
      ic        = imn -1
      do 10 i=imn,imp2f,2
      ic        = ic  +1
      wf(i,j,k,n) = wc(ic,jc,kc,n)
  10  continue
c
c  for the interior points
c 
c--------- fill in fine grid in i direction 
c
      do 50 n=1,5
      do 40 k=kmn,kmp2f,2
c 
      do 30 j=jmn  ,jmp2f,2
      do 30 i=imnp1,imp1f,2
      wf(i  ,j,k,n) = .25*wf(i-1,j,k,n) + .75*wf(i+1,j,k,n) 
      wf(i-1,j,k,n) = .75*wf(i-1,j,k,n) + .25*wf(i+1,j,k,n) 
   30 continue
c 
c--------- fill in fine grid in j direction 
c 
      do 40 j=jmnp1,jmp1f,2
      do 40 i=imn  ,imp1f
      wf(i,j  ,k,n) = .25*wf(i,j-1,k,n) + .75*wf(i,j+1,k,n) 
      wf(i,j-1,k,n) = .75*wf(i,j-1,k,n) + .25*wf(i,j+1,k,n) 
   40 continue
c 
c--------- fill in fine grid in k direction 
c 
      do 50 j=jmn,jmp1f
      do 50 k=kmnp1,kmp1f,2
      do 50 i=imn  ,imp1f
      wf(i,j,k  ,n) = .25*wf(i,j,k-1,n) + .75*wf(i,j,k+1,n) 
      wf(i,j,k-1,n) = .75*wf(i,j,k-1,n) + .25*wf(i,j,k+1,n) 
   50 continue
c
c       do the interpolation for the eddy-viscosiy
c
        do 45 kf=kmn,kmf,2
c
      kc= ( kf + 2 ) / 2
        do 45 jf=jmn,jmf,2
      jc= ( jf + 2 ) / 2
c
        do 46 iff=imn,imf,2
      ic= ( iff + 2 ) / 2
c
      eomuf(iff  ,jf  ,kf  )=sfi*(
     .        27.*eomuc(ic  ,jc  ,kc  ) + 3.*eomuc(ic+1,jc+1,kc  )   +
     .         9.*eomuc(ic  ,jc  ,kc+1) +    eomuc(ic+1,jc+1,kc+1)   +
     .       9.*( eomuc(ic+1,jc  ,kc  ) +    eomuc(ic  ,jc+1,kc  ) ) +
     .       3.*( eomuc(ic+1,jc  ,kc+1) +    eomuc(ic  ,jc+1,kc+1) ) )
      eomuf(iff+1,jf  ,kf  )=sfi*(
     .        27.*eomuc(ic+1,jc  ,kc  ) + 3.*eomuc(ic  ,jc+1,kc  )   +
     .         9.*eomuc(ic+1,jc  ,kc+1) +    eomuc(ic  ,jc+1,kc+1)   +
     .       9.*( eomuc(ic  ,jc  ,kc  ) +    eomuc(ic+1,jc+1,kc  ) ) +
     .       3.*( eomuc(ic  ,jc  ,kc+1) +    eomuc(ic+1,jc+1,kc+1) ) )
      eomuf(iff  ,jf+1,kf  )=sfi*(
     .        27.*eomuc(ic  ,jc+1,kc  ) + 3.*eomuc(ic+1,jc  ,kc  )   +
     .         9.*eomuc(ic  ,jc+1,kc+1) +    eomuc(ic+1,jc  ,kc+1)   +
     .       9.*( eomuc(ic  ,jc  ,kc  ) +    eomuc(ic+1,jc+1,kc  ) ) +
     .       3.*( eomuc(ic  ,jc  ,kc+1) +    eomuc(ic+1,jc+1,kc+1) ) )
      eomuf(iff+1,jf+1,kf  )=sfi*(
     .        27.*eomuc(ic+1,jc+1,kc  ) + 3.*eomuc(ic  ,jc  ,kc  )   +
     .         9.*eomuc(ic+1,jc+1,kc+1) +    eomuc(ic  ,jc  ,kc+1)   +
     .       9.*( eomuc(ic  ,jc+1,kc  ) +    eomuc(ic+1,jc  ,kc  ) ) +
     .       3.*( eomuc(ic  ,jc+1,kc+1) +    eomuc(ic+1,jc  ,kc+1) ) )
c
c  next higher k location
c
      eomuf(iff  ,jf  ,kf+1)=sfi*(
     .        27.*eomuc(ic  ,jc  ,kc+1) + 3.*eomuc(ic+1,jc+1,kc+1)   +
     .         9.*eomuc(ic  ,jc  ,kc  ) +    eomuc(ic+1,jc+1,kc  )   +
     .       9.*( eomuc(ic+1,jc  ,kc+1) +    eomuc(ic  ,jc+1,kc+1) ) +
     .       3.*( eomuc(ic+1,jc  ,kc  ) +    eomuc(ic  ,jc+1,kc  ) ) )
      eomuf(iff+1,jf  ,kf+1)=sfi*(
     .        27.*eomuc(ic+1,jc  ,kc+1) + 3.*eomuc(ic  ,jc+1,kc+1)   +
     .         9.*eomuc(ic+1,jc  ,kc  ) +    eomuc(ic  ,jc+1,kc  )   +
     .       9.*( eomuc(ic  ,jc  ,kc+1) +    eomuc(ic+1,jc+1,kc+1) ) +
     .       3.*( eomuc(ic  ,jc  ,kc  ) +    eomuc(ic+1,jc+1,kc  ) ) )
      eomuf(iff  ,jf+1,kf+1)=sfi*(
     .        27.*eomuc(ic  ,jc+1,kc+1) + 3.*eomuc(ic+1,jc  ,kc+1)   +
     .         9.*eomuc(ic  ,jc+1,kc  ) +    eomuc(ic+1,jc  ,kc  )   +
     .       9.*( eomuc(ic  ,jc  ,kc+1) +    eomuc(ic+1,jc+1,kc+1) ) +
     .       3.*( eomuc(ic  ,jc  ,kc  ) +    eomuc(ic+1,jc+1,kc  ) ) )
      eomuf(iff+1,jf+1,kf+1)=sfi*(
     .        27.*eomuc(ic+1,jc+1,kc+1) + 3.*eomuc(ic  ,jc  ,kc+1)   +
     .         9.*eomuc(ic+1,jc+1,kc  ) +    eomuc(ic  ,jc  ,kc  )   +
     .       9.*( eomuc(ic  ,jc+1,kc+1) +    eomuc(ic+1,jc  ,kc+1) ) +
     .       3.*( eomuc(ic  ,jc+1,kc  ) +    eomuc(ic+1,jc  ,kc  ) ) )
c
 46   continue
c
        if(iturb.eq.2) then
        do 47 iff=imn,imf,2
      ic= ( iff + 2 ) / 2
c
      turv1f(iff  ,jf  ,kf  )=sfi*(
     .        27.*turv1c(ic  ,jc  ,kc  ) + 3.*turv1c(ic+1,jc+1,kc  )   +
     .         9.*turv1c(ic  ,jc  ,kc+1) +    turv1c(ic+1,jc+1,kc+1)   +
     .       9.*( turv1c(ic+1,jc  ,kc  ) +    turv1c(ic  ,jc+1,kc  ) ) +
     .       3.*( turv1c(ic+1,jc  ,kc+1) +    turv1c(ic  ,jc+1,kc+1) ) )
      turv1f(iff+1,jf  ,kf  )=sfi*(
     .        27.*turv1c(ic+1,jc  ,kc  ) + 3.*turv1c(ic  ,jc+1,kc  )   +
     .         9.*turv1c(ic+1,jc  ,kc+1) +    turv1c(ic  ,jc+1,kc+1)   +
     .       9.*( turv1c(ic  ,jc  ,kc  ) +    turv1c(ic+1,jc+1,kc  ) ) +
     .       3.*( turv1c(ic  ,jc  ,kc+1) +    turv1c(ic+1,jc+1,kc+1) ) )
      turv1f(iff  ,jf+1,kf  )=sfi*(
     .        27.*turv1c(ic  ,jc+1,kc  ) + 3.*turv1c(ic+1,jc  ,kc  )   +
     .         9.*turv1c(ic  ,jc+1,kc+1) +    turv1c(ic+1,jc  ,kc+1)   +
     .       9.*( turv1c(ic  ,jc  ,kc  ) +    turv1c(ic+1,jc+1,kc  ) ) +
     .       3.*( turv1c(ic  ,jc  ,kc+1) +    turv1c(ic+1,jc+1,kc+1) ) )
      turv1f(iff+1,jf+1,kf  )=sfi*(
     .        27.*turv1c(ic+1,jc+1,kc  ) + 3.*turv1c(ic  ,jc  ,kc  )   +
     .         9.*turv1c(ic+1,jc+1,kc+1) +    turv1c(ic  ,jc  ,kc+1)   +
     .       9.*( turv1c(ic  ,jc+1,kc  ) +    turv1c(ic+1,jc  ,kc  ) ) +
     .       3.*( turv1c(ic  ,jc+1,kc+1) +    turv1c(ic+1,jc  ,kc+1) ) )
c
c  next higher k location
c
c
      turv1f(iff  ,jf  ,kf+1)=sfi*(
     .        27.*turv1c(ic  ,jc  ,kc+1) + 3.*turv1c(ic+1,jc+1,kc+1)   +
     .         9.*turv1c(ic  ,jc  ,kc  ) +    turv1c(ic+1,jc+1,kc  )   +
     .       9.*( turv1c(ic+1,jc  ,kc+1) +    turv1c(ic  ,jc+1,kc+1) ) +
     .       3.*( turv1c(ic+1,jc  ,kc  ) +    turv1c(ic  ,jc+1,kc  ) ) )
      turv1f(iff+1,jf  ,kf+1)=sfi*(
     .        27.*turv1c(ic+1,jc  ,kc+1) + 3.*turv1c(ic  ,jc+1,kc+1)   +
     .         9.*turv1c(ic+1,jc  ,kc  ) +    turv1c(ic  ,jc+1,kc  )   +
     .       9.*( turv1c(ic  ,jc  ,kc+1) +    turv1c(ic+1,jc+1,kc+1) ) +
     .       3.*( turv1c(ic  ,jc  ,kc  ) +    turv1c(ic+1,jc+1,kc  ) ) )
      turv1f(iff  ,jf+1,kf+1)=sfi*(
     .        27.*turv1c(ic  ,jc+1,kc+1) + 3.*turv1c(ic+1,jc  ,kc+1)   +
     .         9.*turv1c(ic  ,jc+1,kc  ) +    turv1c(ic+1,jc  ,kc  )   +
     .       9.*( turv1c(ic  ,jc  ,kc+1) +    turv1c(ic+1,jc+1,kc+1) ) +
     .       3.*( turv1c(ic  ,jc  ,kc  ) +    turv1c(ic+1,jc+1,kc  ) ) )
      turv1f(iff+1,jf+1,kf+1)=sfi*(
     .        27.*turv1c(ic+1,jc+1,kc+1) + 3.*turv1c(ic  ,jc  ,kc+1)   +
     .         9.*turv1c(ic+1,jc+1,kc  ) +    turv1c(ic  ,jc  ,kc  )   +
     .       9.*( turv1c(ic  ,jc+1,kc+1) +    turv1c(ic+1,jc  ,kc+1) ) +
     .       3.*( turv1c(ic  ,jc+1,kc  ) +    turv1c(ic+1,jc  ,kc  ) ) )
 47   continue
        endif
c
        if(iturb.eq.3) then
        do 48 iff=imn,imf,2
      ic= ( iff + 2 ) / 2
c
      turv2f(iff  ,jf  ,kf  )=sfi*(
     .        27.*turv2c(ic  ,jc  ,kc  ) + 3.*turv2c(ic+1,jc+1,kc  )   +
     .         9.*turv2c(ic  ,jc  ,kc+1) +    turv2c(ic+1,jc+1,kc+1)   +
     .       9.*( turv2c(ic+1,jc  ,kc  ) +    turv2c(ic  ,jc+1,kc  ) ) +
     .       3.*( turv2c(ic+1,jc  ,kc+1) +    turv2c(ic  ,jc+1,kc+1) ) )
      turv2f(iff+1,jf  ,kf  )=sfi*(
     .        27.*turv2c(ic+1,jc  ,kc  ) + 3.*turv2c(ic  ,jc+1,kc  )   +
     .         9.*turv2c(ic+1,jc  ,kc+1) +    turv2c(ic  ,jc+1,kc+1)   +
     .       9.*( turv2c(ic  ,jc  ,kc  ) +    turv2c(ic+1,jc+1,kc  ) ) +
     .       3.*( turv2c(ic  ,jc  ,kc+1) +    turv2c(ic+1,jc+1,kc+1) ) )
      turv2f(iff  ,jf+1,kf  )=sfi*(
     .        27.*turv2c(ic  ,jc+1,kc  ) + 3.*turv2c(ic+1,jc  ,kc  )   +
     .         9.*turv2c(ic  ,jc+1,kc+1) +    turv2c(ic+1,jc  ,kc+1)   +
     .       9.*( turv2c(ic  ,jc  ,kc  ) +    turv2c(ic+1,jc+1,kc  ) ) +
     .       3.*( turv2c(ic  ,jc  ,kc+1) +    turv2c(ic+1,jc+1,kc+1) ) )
      turv2f(iff+1,jf+1,kf  )=sfi*(
     .        27.*turv2c(ic+1,jc+1,kc  ) + 3.*turv2c(ic  ,jc  ,kc  )   +
     .         9.*turv2c(ic+1,jc+1,kc+1) +    turv2c(ic  ,jc  ,kc+1)   +
     .       9.*( turv2c(ic  ,jc+1,kc  ) +    turv2c(ic+1,jc  ,kc  ) ) +
     .       3.*( turv2c(ic  ,jc+1,kc+1) +    turv2c(ic+1,jc  ,kc+1) ) )
c
c  next higher k location
c
c
      turv2f(iff  ,jf  ,kf+1)=sfi*(
     .        27.*turv2c(ic  ,jc  ,kc+1) + 3.*turv2c(ic+1,jc+1,kc+1)   +
     .         9.*turv2c(ic  ,jc  ,kc  ) +    turv2c(ic+1,jc+1,kc  )   +
     .       9.*( turv2c(ic+1,jc  ,kc+1) +    turv2c(ic  ,jc+1,kc+1) ) +
     .       3.*( turv2c(ic+1,jc  ,kc  ) +    turv2c(ic  ,jc+1,kc  ) ) )
      turv2f(iff+1,jf  ,kf+1)=sfi*(
     .        27.*turv2c(ic+1,jc  ,kc+1) + 3.*turv2c(ic  ,jc+1,kc+1)   +
     .         9.*turv2c(ic+1,jc  ,kc  ) +    turv2c(ic  ,jc+1,kc  )   +
     .       9.*( turv2c(ic  ,jc  ,kc+1) +    turv2c(ic+1,jc+1,kc+1) ) +
     .       3.*( turv2c(ic  ,jc  ,kc  ) +    turv2c(ic+1,jc+1,kc  ) ) )
      turv2f(iff  ,jf+1,kf+1)=sfi*(
     .        27.*turv2c(ic  ,jc+1,kc+1) + 3.*turv2c(ic+1,jc  ,kc+1)   +
     .         9.*turv2c(ic  ,jc+1,kc  ) +    turv2c(ic+1,jc  ,kc  )   +
     .       9.*( turv2c(ic  ,jc  ,kc+1) +    turv2c(ic+1,jc+1,kc+1) ) +
     .       3.*( turv2c(ic  ,jc  ,kc  ) +    turv2c(ic+1,jc+1,kc  ) ) )
      turv2f(iff+1,jf+1,kf+1)=sfi*(
     .        27.*turv2c(ic+1,jc+1,kc+1) + 3.*turv2c(ic  ,jc  ,kc+1)   +
     .         9.*turv2c(ic+1,jc+1,kc  ) +    turv2c(ic  ,jc  ,kc  )   +
     .       9.*( turv2c(ic  ,jc+1,kc+1) +    turv2c(ic+1,jc  ,kc+1) ) +
     .       3.*( turv2c(ic  ,jc+1,kc  ) +    turv2c(ic+1,jc  ,kc  ) ) )
 48   continue
        endif
c
   45 continue
c
c
c     compute pressure on the fine grid
c
      do 60 k = kmn,kmf+1
      do 60 j = jmn,jmf+1
      do 60 i = imn,imf+1
c
      pf(i,j,k) = gm1*(wf(i,j,k,5)
     .          - 0.5*(wf(i,j,k,2)**2 + wf(i,j,k,3)**2 + wf(i,j,k,4)**2)
     .          /      wf(i,j,k,1) )
      pf(i,j,k) = max(pf(i,j,k),0.001)
  60  continue
c
c
      return
      end
