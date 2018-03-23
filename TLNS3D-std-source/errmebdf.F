      subroutine  errmebdf(imn,jmn,kmn,imx,jmx,kmx,imp2,jmp2,kmp2,
     .                    delw,dtl,vol,dtphy,mx1dwk,errw,ibloc,nodes,
     .                    w,wbar,ntorder)
c=======================================================================
c 
c     objective : computes temporal error in MEBDF
c                 time-stepping using solutions at
c                 different stages
c
c     iniitial coding : v.n.vatsa (Nov. 2003)
c 
c 

      dimension delw(imp2,jmp2,kmp2,5),dtl(imp2,jmp2,kmp2)
c
      dimension errw(mx1dwk),vol(imp2,jmp2,kmp2),
     .          w(imp2,jmp2,kmp2,5),
     .          wbar (imp2,jmp2,kmp2,5,ntorder-1)
c
      common/rkerr/ erkmaxl(5),erkmaxg(5),erkrms2(5),erkrmsg(5),
     .              totptse
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      common/rk/    cfl,c(6),qfil(6),beta(6),vt,hm,mstage
c
c     initialize global values of errors for ibloc = 1
c
      if (ibloc .eq. 1) then
         do nv = 1,5
           erkrms2(nv)      = 0.
           erkrmsg(nv)      = 0.
           erkmaxg(nv)      = 0.
           totptse          = 0
         enddo
      endif
c
      do 100 nv=1,5
c
      erkmaxl(nv)      = 0.
      do 90 k=kmn+1,kmx
      do 90 j=jmn+1,jmx
      iloc   = 0
      do 50 i=imn+1,imx
      iloc   = iloc + 1
      dt        = cfl*((1.-vt)*dtl(i,j,k) + vt*dtmin)
cvn      errw(iloc) = (w(i,j,k,nv) - wbar(i,j,k,nv,1))/dt
         errw(iloc) =  w(i,j,k,nv) - wbar(i,j,k,nv,1)
         erkrms2(nv)     = erkrms2(nv) + errw(iloc)**2
   50 continue
      i         = isamax(iloc,errw,1)
      if (abs(errw(i)).gt.abs(erkmaxl(nv))) then 
        erkmaxl(nv)     = abs(errw(i))
      endif 
   90 continue
c
c     keep track of total no. of interior points on all blocks and
c     evaluate global rms values including points at current block
c
      if(nv.eq.1)
     .totptse   = totptse+ (imx - imn) * (jmx - jmn) * (kmx - kmn)
      erkrmsg(nv)    = sqrt(erkrms2(nv)/totptse)
c
      if (ibloc .eq. 1) then
         erkmaxg(nv)   = erkmaxl(nv)
cdbg     return
      endif
c
c     check if erkmaxl on the current bloc is greater than erkmaxg
c
      if(abs(erkmaxl(nv)) .gt. abs(erkmaxg(nv))) then
         erkmaxg(nv)   = erkmaxl(nv)
      endif
c
c     check for Nan's
c
      erkmx = erkmaxl(nv)
      if ( (erkmx+erkmx.eq.erkmx.and.abs(erkmx).gt.1.d-100).or.
     .      .not.(erkmx.lt.abs(erkmx).or.erkmx.ge.0.d0)) then
       write (iwrit,901) ibloc
 901   format(2x,'Nan detected in erkmax at ibloc = ',i5)
c      stop 1001
       call ERREXIT(nodes) 
      endif
c
      erkrms = erkrmsg(nv)
      if ( (erkrms+erkrms.eq.erkrms.and.abs(erkrms).gt.1.d-100).or.
     .      .not.(erkrms.lt.abs(erkrms).or.erkrms.ge.0.d0)) then
       write (iwrit,902) ibloc
 902   format(2x,'Nan detected in erkrms at ibloc = ',i5)
c      stop 1001
       call ERREXIT(nodes) 
      endif
 100  continue
c
      return
      end
