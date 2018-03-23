      subroutine prolong(iftf,igrid,imn,jmn,kmn,
     .               imxf,jmxf,kmxf,imp1f,jmp1f,kmp1f,imp2f,jmp2f,kmp2f,
     .               imxc,jmxc,kmxc,imp1c,jmp1c,kmp1c,imp2c,jmp2c,kmp2c,
     .               wf,wc,dwf,wk2da,wk2db,pf,w1c,ibloc )
c
c

      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      common/rkrsm/ smoopi,smoopj,smoopk,smoopic,smoopjc,smoopkc
c
c
      dimension wf(imp2f,jmp2f,kmp2f,5),  wc(imp2c,jmp2c,kmp2c,5),
     .          pf(imp2f,jmp2f,kmp2f  ), w1c(imp2c,jmp2c,kmp2c,5),
     .         dwf(imp2f,jmp2f,kmp2f,5),
     .          wk2da(1),wk2db(1)
c
c
      if(ibloc.eq.1) iftf = iftf +1
c
      kmnp1 = kmn + 1
      jmnp1 = jmn + 1
      imnp1 = imn + 1
c
c         2h    2h       2h
c     delw   = w(new) - w(old)
c
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
      dwf(i,j,k,n) = wc(ic,jc,kc,n) - w1c(ic,jc,kc,n)
  10  continue
c 
c--------- fill in fine grid in i direction 
c
      do 50 n=1,5
      do 40 k=kmn,kmp2f,2
c     do 24 j=jmn,jmp2f,2
c     dwf(imn  ,j,k,n)  = 0.
c     dwf(imp2f,j,k,n)  = 0.
c  24 continue
c 
      do 30 j=jmn  ,jmp2f,2
      do 30 i=imnp1,imp1f,2
      dwf(i  ,j,k,n) = .25*dwf(i-1,j,k,n) + .75*dwf(i+1,j,k,n) 
      dwf(i-1,j,k,n) = .75*dwf(i-1,j,k,n) + .25*dwf(i+1,j,k,n) 
   30 continue
c 
c--------- fill in fine grid in j direction 
c 
c     do 34 i=imn,imp1f
c     dwf(i,jmp2f,k,n)  = 0.
c  34 continue
c 
      do 40 j=jmnp1,jmp1f,2
      do 40 i=imn  ,imp1f
      dwf(i,j  ,k,n) = .25*dwf(i,j-1,k,n) + .75*dwf(i,j+1,k,n) 
      dwf(i,j-1,k,n) = .75*dwf(i,j-1,k,n) + .25*dwf(i,j+1,k,n) 
   40 continue
c 
c--------- fill in fine grid in k direction 
c 
      do 50 j=jmn,jmp1f
c     do 44 i=imn,imp1f
c     dwf(i,j,kmp2f,n)  = 0.
c  44 continue
c 
      if(i2dfl.eq.1) then
        k = kmnp1
        do 55 i=imn  ,imp1f
  55    dwf(i,j,k,n) = dwf(i,j,kmn,n)
      else
c
        do 54 k=kmnp1,kmp1f,2
        do 54 i=imn  ,imp1f
        dwf(i,j,k  ,n) = .25*dwf(i,j,k-1,n) + .75*dwf(i,j,k+1,n) 
        dwf(i,j,k-1,n) = .75*dwf(i,j,k-1,n) + .25*dwf(i,j,k+1,n) 
   54   continue
      endif
   50 continue
c
c
      call psmooc (imn,  jmn,  kmn,  imxf,  jmxf,kmxf,
     .             imp1f,jmp1f,kmp1f,imp2f,jmp2f,kmp2f,
     .             wk2da,wk2db,dwf,smoopic,smoopjc,smoopkc )
c
c
c   h        h        h    2h
c  w(new) = w(old) + I  { delw }
c                     2h
c
c  apply correction at only interior points
c
c     first check to make sure rho doesnot change by more than
c     half of its previous value anywhere in the field
c
      delmx  = 0.25
      tol    = 1.e-12
c
      do 60 k=kmnp1,kmxf
      do 60 j=jmnp1,jmxf
      do 60 i=imnp1,imxf
      rho         = wf(i,j,k,1)
      rhoi        = 1./rho
      ut          = wf(i,j,k,2)*rhoi
      vt          = wf(i,j,k,3)*rhoi
      wt          = wf(i,j,k,4)*rhoi
      drho        = dwf(i,j,k,1)
c
      dratio      = abs(drho*rhoi) + tol 
      check       = dratio - delmx 
      sfac        = delmx/dratio
      sfac1       = cvmgm(1.,sfac,check)
c
      dratio      = abs(dwf(i,j,k,5)/wf(i,j,k,5)) + tol 
      check       = dratio - delmx 
      sfac        = delmx/dratio
      sfac2       = cvmgm(1.,sfac,check)
      sfac1       = min(sfac1,sfac2)
c
      qq          = 0.5*(ut*ut +vt*vt +wt*wt)
      ps          = gm1*(wf(i,j,k,5) - rho*qq)
      dp          = gm1*(qq*dwf(i,j,k,1) -ut*dwf(i,j,k,2)
     .                  -vt*dwf(i,j,k,3) -wt*dwf(i,j,k,4)
     .                     +dwf(i,j,k,5))
c
      dratio      = abs(dp/ps) + tol
      check       = dratio - delmx
      sfac        = delmx/dratio
      sfac2       = cvmgm(1.,sfac,check)
c
      tmp1         = min(sfac1,sfac2)
      dwf(i,j,k,1) = dwf(i,j,k,1)*tmp1
      dwf(i,j,k,2) = dwf(i,j,k,2)*tmp1
      dwf(i,j,k,3) = dwf(i,j,k,3)*tmp1
      dwf(i,j,k,4) = dwf(i,j,k,4)*tmp1
      dwf(i,j,k,5) = dwf(i,j,k,5)*tmp1
   60 continue
c
      do 70 n=1,5
      do 70 k=kmnp1,kmxf
      do 70 j=jmnp1,jmxf
      do 70 i=imnp1,imxf
      wf(i,j,k,n)= wf(i,j,k,n) + dwf(i,j,k,n)
   70 continue
c
c  computing the pressures at the interior points
c
      do 75 k=kmnp1,kmxf
      do 75 j=jmnp1,jmxf
      do 75 i=imnp1,imxf
      pf(i,j,k) = gm1*(wf(i,j,k,5)
     .               -.5*(wf(i,j,k,2)**2 + wf(i,j,k,3)**2 +
     .                    wf(i,j,k,4)**2)/ wf(i,j,k,1)     )
      pf(i,j,k) = max(pf(i,j,k),0.001)
   75 continue
c
c
      return
      end
