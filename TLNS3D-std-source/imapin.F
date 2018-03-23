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
      subroutine imapin(imn,jmn,kmn,im,jm,km,
     .                  mgrlev,ngrid,igrid,mbloc,nbloc,
     .                  imap,msegt,msegn,nseg,mx1dwk,const1,
     .                  const2,const3,const4,const5,
     .                  wk2da,wk2db,wk2dc,wk2dd,wk2de,wk2df,ivisb,
     .                  iturbb,itrb1,itrb2,jtrb1,jtrb2,ktrb1,ktrb2,
     .                  nodes )
c-----------------------------------------------------------------------
c
c     object : to read in and set up the topological mapping 
c              funcions, 'imap'
c
c     initial coding : by m.d.sanetrik (August, 1990)

c     latest mods for  new engine and multiple jet-stream
c     conditions : by v.n.vatsa, Jan. 1998
c     engine b.c.'s: based on Ted Reyhner's  work
c     
c     Note: The third line entries in map file for segment b.c.'s
c           have different meaning. Read the comments later on
c           in this subroutine for new varaibales const1... const5
c
c
      dimension     nodes(1)
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      common/fld/   gamma,gm1,dgm1,gogm1,rm,rho0,p0,ei0,h0,c0,u0,v0,w0,
     .              ca,sa,pr,prt,rey,suthc,tref,i2dfl,iturb
c
      dimension     im(mgrlev,mbloc),  jm(mgrlev,mbloc),
     .              km(mgrlev,mbloc),
     .              imap(msegt,msegn,mbloc,mgrlev), nseg(mbloc)
c
      dimension     itrb1(mgrlev,mbloc),  itrb2(mgrlev,mbloc),
     .              jtrb1(mgrlev,mbloc),  jtrb2(mgrlev,mbloc),
     .              ktrb1(mgrlev,mbloc),  ktrb2(mgrlev,mbloc),
     .              iturbb(mbloc)
c
c***********************************************************************
c
c
c     Instead of using the variable names twotref,cmassw,
c     cporous, pchambr and epsrad, which are meaningful only
c     for viscous walls, we are renaming these 5 varibales
c     that are read from map file as const1.....const5
c     with the following meanings (Jan. 1998)
c
c     For viscous walls (bctyp=3) :
c         const1   = twotref
c         const2   = cmassw
c         const3   = cporous
c         const4   = pchambr
c         const5   = epsrad
c
c     For extrapolation condition (bctyp=5) :
c
c         const1   = pback (back pressure)
c        (if const1 is zero, pressure is also extrapolated)
c
c     For fan-upstream engine boundary condition (bctyp=8):
c
c         const1   = epsfu (stream-tube contraction ratio)
c
c     For powered engine (massflow) boundary condition (bctyp=9):
c
c         const1   = epseng (engine/fan massflow to free-stream ratio)
c         const2   = relax (relaxation factor >= 0.05 to 1.0)
c
c     For jet/nozzle boundary condition (bctyp=10) :
c
c         const1   = ptjet (jet total pressure/ref. static pressure)
c         const2   = ttjet (jet total temp./ref. static temp.)
c         const3   = rmjet (nominal jet Mach number)
c         const4   = alfaj (jet inclination angle in x-y plane)
c         const5   = betaj (jet inclination angle in x-z plane)
c        (alfaj, betaj are usually set to zero)
c
c     For flow-streams with different total pressure & temp. (bctyp=11)
c         and specified inflow angles
c
c         const1   = ptstr (stream total pressure/ref. static pressure)
c         const2   = ttstr (stream total temp./ref. static temp.)
c         const3   = relax (relaxation factor : default = 1.0)
c         const4   = alfstr (stream inclination angle in x-y plane)
c         const5   = betstr (stream inclination angle in x-z plane)
c        (alfstr, betstr are usually set to zero)
c
c     For flow-streams with different total pressure & temp. (bctyp=12)
c         and inflow angles aligned with grid lines
c
c         const1   = ptstr (stream total pressure/ref. static pressure)
c         const2   = ttstr (stream total temp./ref. static temp.)
c         const3   = relax (relaxation factor : default = 1.0)
c
      dimension     ivisb(3,mbloc),      const1(msegn,mbloc),
     .              const2(msegn,mbloc),const3(msegn,mbloc),
     .              const4(msegn,mbloc),const5(msegn,mbloc)
c****************************************************************************
c
      dimension nmn1(6),nmx1(6),nmn2(6),nmx2(6)
c
      dimension wk2da(mx1dwk,mx1dwk),wk2db(mx1dwk,mx1dwk),
     .          wk2dc(mx1dwk,mx1dwk),wk2dd(mx1dwk,mx1dwk),
     .          wk2de(mx1dwk,mx1dwk),wk2df(mx1dwk,mx1dwk)
c
c
c*************************************************************************
c--------------------  mapping function description ---------------------
c
c      imap  : mapping function containing topological information
c      msegt : maximum no. of types of operations/boundary conditons
c              (currently set to 20)
c      msegn : maximum no. of segments permitted on all faces of a block
c              (currently set to 20)
c
c      imap(1 , ) : specifies boundary/operation type
c      imap(2 , ) : specifies face number (1-6)
c      imap(3 , ) : specifies beginning of direction 1
c      imap(4 , ) : specifies end       of direction 1
c      imap(5 , ) : specifies beginning of direction 2
c      imap(6 , ) : specifies end       of direction 2
c
c      imap(7,  ) : for bctype 4 & 5, specifies order of extrapolation
c                   0  for zero  order extrapolation
c                   1  for first order extrapolation
c
c      imap(8,  ) : if the boundary type is symmetry
c                   takes the value 1,2 or 3 to indicate symmetry 
c                   about x, y or z = constant planes, respectively
c 
c  if the boundary type is a cut/interface
c  additional information described below is required 
c
c      imap(7 , ) : specifies block no. of source segment
c     |imap(8 , )|: specifies face no. of source segment
c      imap(8 , ) > 0 corresponds to direction 1 of source segment
c                      matching with direction 1 of target segment
c                      and same  for direction 2
c      imap(8 , ) < 0 corresponds to direction 2 of source segment
c                      matching with direction 1 of target segment
c                      and vice-versa
c
c      imap(9 , ) : specifies beginning of direction 1 of source segment
c      imap(10, ) : specifies end       of direction 1 of source segment
c      imap(11, ) : specifies beginning of direction 2 of source segment
c      imap(12, ) : specifies end       of direction 2 of source segment
c
c      turbulence/transitional related information
c    
c      imap(13, ) : specifies if there is turbulent flow on this segment
c                   1  corresponds to turbulent flow
c                   0  corresponds to laminar   flow
c      imap(14, ) : begining index in direction 1 for turbulent flow
c      imap(15, ) : ending   index in direction 1 for turbulent flow
c      imap(16, ) : begining index in direction 2 for turbulent flow
c      imap(17, ) : ending   index in direction 2 for turbulent flow
c      imap(18, ) : begining index in normal direction for fmax 
c      imap(19, ) : ending   index in normal direction for fmax 
c      imap(20, ) : ending   index in normal direction for turb. flow 
c
c
c--------------------  boundary/operation type description ---------------------
c
c      nbctype    = imap(1 , )
c                 =-1  corresponds to a cut that patches to another block
c                      in a non-C0 fashion
c                    (no. of points on donor and receiving segments differ)
c                 = 0  corresponds to an inner cut that maps to same block
c                 = 1  corresponds to a cut that maps to another block
c                 = 2  corresponds to a slip (inviscid) wall
c                 = 3  corresponds to a noslip (viscous) wall
c                 = 4  symmetry condition 
c                 = 5  extrapolation condition (with specified back pressure)
c                 = 6  far-field condition (Riemann invariants)
c                 = 7  singular line
c                 = 8  fan-upstream engine boundary(streamtube contraction)
c                 = 9  engine/fan boundary (specified mass flow)
c                 =10  jet/nozzle boundary (specify M, total p,T & flow angle )
c                 =11  Total pressure/temperature & flow angle b.c.
c                 =12  Total pressure/temperature & flow aligned with grid b.c.
c
c***************************************************************************
c
c     ivisb(1)    = 1  viscous fluxes evaluated along i-direction
c     ivisb(2)    = 1  viscous fluxes evaluated along j-direction
c     ivisb(3)    = 1  viscous fluxes evaluated along k-direction
c
c     iturbb      = 1  turbulent flow in this bloc, laminar otherwise
c     itrb1,itrb2      index range for turbulent flow in i-direction
c     jtrb1,jtrb2      index range for turbulent flow in j-direction
c     ktrb1,ktrb2      index range for turbulent flow in k-direction
c                      if itrb1,itrb2 = 0, complete block is turbulent
c                      if itrb1,itrb2 =negative, complete block is laminar
c                      same is true for jtrb1,jtrb2,ktrb1,ktrb2
c***************************************************************************
c
      rad       = 45./atan(1.)

      read(imapfil,'(1x)')
      read(imapfil,*) idum
      if(idum.ne.nbloc) then
         write(iwrit,'(2x,"no. of blocks on map file ",i3,
     .                "donot match the no. of blocks on grid file",
     .                 i3)') idum,nbloc
         call ERREXIT (nodes)
      endif
c
c     initialize imap array (vatsa 3-17-97)
c
      do m1=1,msegt
      do m2=1,msegn
      do m3=1,mbloc
      do m4=1,mgrlev
      imap(m1,m2,m3,m4) = 0
      enddo
      enddo
      enddo
      enddo
c
c----------  read in imap array and shift due to ghost cells  -------------
c
      do 100 ibloc = 1,nbloc
      read(imapfil,'(1x)')
      read(imapfil,'(1x)')
      read(imapfil,*)  nseg(ibloc),ivisb(1,ibloc),ivisb(2,ibloc),
     .                             ivisb(3,ibloc) ,
     .                 iturbb(ibloc),itrb1(1,ibloc),itrb2(1,ibloc),
     .                               jtrb1(1,ibloc),jtrb2(1,ibloc),
     .                               ktrb1(1,ibloc),ktrb2(1,ibloc)
      if(i2dfl.eq.1) ivisb(3,ibloc) = 0
      read(imapfil,'(1x)')
      read(imapfil,'(1x)')
c
c     check to see if nseg is within allowed dimension
c
      if(nseg(ibloc) .gt. msegn) then
       write(iwrit,'(1x)')
       write(iwrit,'("**Error: the number of segments for block #",i3,
     .               " (",i3,")"/8x," is larger than the allowed ",
     .               "dimension  (",i3,")")') ibloc,nseg(ibloc),msegn
         write(iwrit,'(1x)')
         call ERREXIT (nodes)
      end if
c
      write (iwrit,'(//10x,"***  Boundary conditions for block  #",
     .               5x,i5,5x,"***")') ibloc
c
      do 100 iseg  = 1,nseg(ibloc)
c
      read(imapfil,'(1x)')
      read(imapfil,*) idum,idum,(imap(n,iseg,ibloc,1),n= 1,12)
      read(imapfil,*)           (imap(n,iseg,ibloc,1),n=13,20)
      read(imapfil,*)            const1(iseg,ibloc) 
     .                          ,const2 (iseg,ibloc) 
     .                          ,const3 (iseg,ibloc) 
     .                          ,const4 (iseg,ibloc) 
     .                          ,const5 (iseg,ibloc) 
c
c
      if(imap(2,iseg,ibloc,1).eq.1) 
     .  write (iwrit,'(/" On imin face: j= ",i5," to ",i5," and k= ",
     .                                     i5," to ",i5)')
     .  imap(3,iseg,ibloc,1),imap(4,iseg,ibloc,1),
     .  imap(5,iseg,ibloc,1),imap(6,iseg,ibloc,1)
c
      if(imap(2,iseg,ibloc,1).eq.2) 
     .  write (iwrit,'(/" On imax face: j= ",i5," to ",i5," and k= ",
     .                                     i5," to ",i5)')
     .  imap(3,iseg,ibloc,1),imap(4,iseg,ibloc,1),
     .  imap(5,iseg,ibloc,1),imap(6,iseg,ibloc,1)
c
      if(imap(2,iseg,ibloc,1).eq.3) 
     .  write (iwrit,'(/" On jmin face: k= ",i5," to ",i5," and i= ",
     .                                     i5," to ",i5)')
     .  imap(3,iseg,ibloc,1),imap(4,iseg,ibloc,1),
     .  imap(5,iseg,ibloc,1),imap(6,iseg,ibloc,1)
c
      if(imap(2,iseg,ibloc,1).eq.4) 
     .  write (iwrit,'(/" On jmax face: k= ",i5," to ",i5," and i= ",
     .                                     i5," to ",i5)')
     .  imap(3,iseg,ibloc,1),imap(4,iseg,ibloc,1),
     .  imap(5,iseg,ibloc,1),imap(6,iseg,ibloc,1)
c
      if(imap(2,iseg,ibloc,1).eq.5) 
     .  write (iwrit,'(/" On kmin face: i= ",i5," to ",i5," and j= ",
     .                                     i5," to ",i5)')
     .  imap(3,iseg,ibloc,1),imap(4,iseg,ibloc,1),
     .  imap(5,iseg,ibloc,1),imap(6,iseg,ibloc,1)
c
      if(imap(2,iseg,ibloc,1).eq.6) 
     .  write (iwrit,'(/" On kmax face: i= ",i5," to ",i5," and j= ",
     .                                     i5," to ",i5)')
     .  imap(3,iseg,ibloc,1),imap(4,iseg,ibloc,1),
     .  imap(5,iseg,ibloc,1),imap(6,iseg,ibloc,1)
c
      if(imap(1,iseg,ibloc,1).le.1) then
          write (iwrit,'("   interface condition",
     .           5x,"bctype  ",i5)') imap(1,iseg,ibloc,1)
          if(imap(1,iseg,ibloc,1).ge.0)
     .    write (iwrit,'(3x,"one-to-one matched  interface with")')
          if(imap(1,iseg,ibloc,1).lt.0)
     .    write (iwrit,'(3x,"patched  interface with")')
c
          if(imap(2,iseg,ibloc,1).eq.1) 
     .    write (iwrit,'(3x,"imin  face: j= ",i5," to ",
     .                  i5," and k= ",i5," to ",i5," of block",i5)')
     .    imap( 9,iseg,ibloc,1),imap(10,iseg,ibloc,1),
     .    imap(11,iseg,ibloc,1),imap(12,iseg,ibloc,1),
     .    imap( 8,iseg,ibloc,1)
c
          if(imap(2,iseg,ibloc,1).eq.2) 
     .    write (iwrit,'(3x,"imax  face: j= ",i5," to ",
     .                  i5," and k= ",i5," to ",i5," of block",i5)')
     .    imap( 9,iseg,ibloc,1),imap(10,iseg,ibloc,1),
     .    imap(11,iseg,ibloc,1),imap(12,iseg,ibloc,1),
     .    imap( 8,iseg,ibloc,1)
c
          if(imap(2,iseg,ibloc,1).eq.3) 
     .    write (iwrit,'(3x,"jmin  face: k= ",i5," to ",
     .                  i5," and i= ",i5," to ",i5," of block",i5)')
     .    imap( 9,iseg,ibloc,1),imap(10,iseg,ibloc,1),
     .    imap(11,iseg,ibloc,1),imap(12,iseg,ibloc,1),
     .    imap( 8,iseg,ibloc,1)
c
          if(imap(2,iseg,ibloc,1).eq.4) 
     .    write (iwrit,'(3x,"jmax  face: j= ",i5," to ",
     .                  i5," and k= ",i5," to ",i5," of block",i5)')
     .    imap( 9,iseg,ibloc,1),imap(10,iseg,ibloc,1),
     .    imap(11,iseg,ibloc,1),imap(12,iseg,ibloc,1),
     .    imap( 8,iseg,ibloc,1)
c
          if(imap(2,iseg,ibloc,1).eq.5) 
     .    write (iwrit,'(3x,"kmin  face: i= ",i5," to ",
     .                  i5," and j= ",i5," to ",i5," of block",i5)')
     .    imap( 9,iseg,ibloc,1),imap(10,iseg,ibloc,1),
     .    imap(11,iseg,ibloc,1),imap(12,iseg,ibloc,1),
     .    imap( 8,iseg,ibloc,1)
c
          if(imap(2,iseg,ibloc,1).eq.6) 
     .    write (iwrit,'(3x,"kmax  face: i= ",i5," to ",
     .                  i5," and j= ",i5," to ",i5," of block",i5)')
     .    imap( 9,iseg,ibloc,1),imap(10,iseg,ibloc,1),
     .    imap(11,iseg,ibloc,1),imap(12,iseg,ibloc,1),
     .    imap( 8,iseg,ibloc,1)
c
       endif
c
      if(imap(1,iseg,ibloc,1).eq.2) 
     .    write (iwrit,'(2x," inviscid or slip wall condition",
     .           5x,"bctype  ",i5)') imap(1,iseg,ibloc,1)
c
      if(imap(1,iseg,ibloc,1).eq.3) then
          write (iwrit,'(2x," viscous wall condition",
     .           5x,"bctype  ",i5)') imap(1,iseg,ibloc,1)
          write (iwrit,'(2x," twall/tref = ",e12.5,2x,
     .                   " cmassw  = ",e12.5,2x,
     .                   " cporous = ",e12.5,4x,
     .                   " pchamber   = ",e12.5,2x,
     .                   " epsrad  = ",e12.5,2x)')
     .    const1(iseg,ibloc), const2(iseg,ibloc),
     .    const3(iseg,ibloc), const4(iseg,ibloc),
     .    const5(iseg,ibloc)
      endif
c
      if(imap(1,iseg,ibloc,1).eq.4) 
     .    write (iwrit,'(2x," symmetry condition",
     .           5x,"bctype  ",i5)') imap(1,iseg,ibloc,1)
c
      if(imap(1,iseg,ibloc,1).eq.5) then
          write (iwrit,'(2x," extrapolation condition",
     .           5x,"bctype  ",i5)') imap(1,iseg,ibloc,1)
          write (iwrit,'(2x," pback/pref = ",e12.5,2x)')
     .    const1(iseg,ibloc)
      endif
c
      if(imap(1,iseg,ibloc,1).eq.6) 
     .    write (iwrit,'(2x," far-field condition",
     .           5x,"bctype  ",i5)') imap(1,iseg,ibloc,1)
c
      if(imap(1,iseg,ibloc,1).eq.7) 
     .    write (iwrit,'(2x," singular line condition",
     .           5x,"bctype  ",i5)') imap(1,iseg,ibloc,1)
c
      if(imap(1,iseg,ibloc,1).eq.8) then
          write (iwrit,'(2x," fan-upstream stream tube contr. ratio",
     .           2x,"condition",
     .           5x,"bctype  ",i5)') imap(1,iseg,ibloc,1)
          write (iwrit,'(2x,"epsfu=",e12.5/)')
     .    const1(iseg,ibloc)
c
          if(const1(iseg,ibloc).le.0.) then
          write (iwrit,'(2x," non-physical fan-upstream b.c.",/)')
          call ERREXIT (nodes)
          endif
      endif
c
      if(imap(1,iseg,ibloc,1).eq.9) then
          write (iwrit,'(2x,"engine/fan (mass-flow) condition",
     .           5x,"bctype  ",i5)') imap(1,iseg,ibloc,1)
          write (iwrit,'(2x,"epseng =",
     .            e12.5,2x,"relax. factor=",e12.5/)')
     .    const1(iseg,ibloc), const2(iseg,ibloc)
c
          if(const1(iseg,ibloc).le.0..or.const2(iseg,ibloc).le.0.)then
          write (iwrit,'(2x," non-physical engine-massflow b.c.",/)')
          call ERREXIT (nodes)
          endif
      endif
c
      if(imap(1,iseg,ibloc,1).eq.10) then
          write (iwrit,'(2x," jet/nozzle inflow/outflow condition",
     .           5x,"bctype  ",i5)') imap(1,iseg,ibloc,1)
          write (iwrit,'(2x," ptjopref= ",e12.5,2x,
     .                   " ttjotref= ",e12.5,2x,
     .                   " rmjet  = ",e12.5,9x,
     .                   "alfaj   = ",e12.5,2x,
     .                   " betaj   = ",e12.5,2x,/)')
     .    const1(iseg,ibloc), const2(iseg,ibloc),
     .    const3(iseg,ibloc), const4(iseg,ibloc),
     .    const5(iseg,ibloc)
c
c     convert alfaj, betaj to radians
c
      const4(iseg,ibloc) = const4(iseg,ibloc)/rad
      const5(iseg,ibloc) = const5(iseg,ibloc)/rad
c
          if(const1(iseg,ibloc).le.0..or.const2(iseg,ibloc).le.0.or.
     .    const3(iseg,ibloc).le.0.) then
          write (iwrit,'(2x," non-physical jet-nozzle b.c.",/)')
          call ERREXIT (nodes)
          endif
      endif
c
      if(imap(1,iseg,ibloc,1).eq.11) then
          write (iwrit,'(2x," total head b.c. at inflow streams",
     .           5x,"bctype  ",i5)') imap(1,iseg,ibloc,1)
          write (iwrit,'(2x," ptstr/pref = ",e12.5,2x,
     .                   " ttstr/tref = ",e12.5,2x,
     .                   " relax  = ",e12.5,3x,
     .                   "alfastr    = ",e12.5,2x,
     .                   " betastr    = ",e12.5,2x,/)')
     .    const1(iseg,ibloc), const2(iseg,ibloc),
     .    const3(iseg,ibloc), const4(iseg,ibloc),
     .    const5(iseg,ibloc)
c
c     convert alfastr, betastr to radians
c
      const4(iseg,ibloc) = const4(iseg,ibloc)/rad
      const5(iseg,ibloc) = const5(iseg,ibloc)/rad
c
          if(const1(iseg,ibloc).le.0..or.const2(iseg,ibloc).le.0)
     .    then
          write (iwrit,'(2x," non-physical total-head b.c.",/)')
          call ERREXIT (nodes)
          endif
      endif
c
      if(imap(1,iseg,ibloc,1).eq.12) then
          write (iwrit,'(2x," total head b.c. at inflow stream",
     .           2x,"which is aligned with grid lines",
     .           5x,"bctype  ",i5)') imap(1,iseg,ibloc,1)
          write (iwrit,'(2x," ptstr/pref = ",e12.5,2x,
     .                   " ttstr/tref = ",e12.5,2x,
     .                   " relax  = ",e12.5,9x,/)')
     .    const1(iseg,ibloc), const2(iseg,ibloc),
     .    const3(iseg,ibloc)
c    
          if(const1(iseg,ibloc).le.0..or.const2(iseg,ibloc).le.0.or.
     .    const3(iseg,ibloc).le.0.) then
          write (iwrit,'(2x," non-physical total-head b.c.",/)')
          call ERREXIT (nodes)
          endif
      endif
c
c     setup for 2-d flows ( zero order extrpolation on k-faces)
c
      if(i2dfl.eq.1) then
        if(imap( 2,iseg,ibloc,1).ge.5) then
           imap( 1,iseg,ibloc,1) = 5
           imap( 7,iseg,ibloc,1) = 0
           imap(13,iseg,ibloc,1) = 0
        endif
      endif
c
      imap( 3,iseg,ibloc,1) = imap( 3,iseg,ibloc,1) +imn -1
      imap( 4,iseg,ibloc,1) = imap( 4,iseg,ibloc,1) +imn -1
      imap( 5,iseg,ibloc,1) = imap( 5,iseg,ibloc,1) +imn -1
      imap( 6,iseg,ibloc,1) = imap( 6,iseg,ibloc,1) +imn -1
      imap( 9,iseg,ibloc,1) = imap( 9,iseg,ibloc,1) +imn -1
      imap(10,iseg,ibloc,1) = imap(10,iseg,ibloc,1) +imn -1
      imap(11,iseg,ibloc,1) = imap(11,iseg,ibloc,1) +imn -1
      imap(12,iseg,ibloc,1) = imap(12,iseg,ibloc,1) +imn -1
c
      imap(14,iseg,ibloc,1) = imap(14,iseg,ibloc,1) +imn -1
      imap(15,iseg,ibloc,1) = imap(15,iseg,ibloc,1) +imn -1
      imap(16,iseg,ibloc,1) = imap(16,iseg,ibloc,1) +imn -1
      imap(17,iseg,ibloc,1) = imap(17,iseg,ibloc,1) +imn -1
      imap(18,iseg,ibloc,1) = imap(18,iseg,ibloc,1) +imn -1
      imap(19,iseg,ibloc,1) = imap(19,iseg,ibloc,1) +imn -1
      imap(20,iseg,ibloc,1) = imap(20,iseg,ibloc,1) +imn -1
c
  100 continue
      write (iwrit,'(//10x,"***  Boundary conditions printout  ",
     .               "is now completed",5x,"***")') 
c
c
c     consistency checks
c
c      check to see if the dimensions of any segment exceed the
c            dimensions for that face
c
      do 200 ibloc=1,nbloc
      do 200 iseg =1,nseg(ibloc)
      if(imap(2,iseg,ibloc,1) .le. 2) then
         nmaxd1 = jm(igrid,ibloc)
         nmaxd2 = km(igrid,ibloc)
      else if(imap(2,iseg,ibloc,1) .ge. 5) then
         nmaxd1 = im(igrid,ibloc)
         nmaxd2 = jm(igrid,ibloc)
      else
         nmaxd1 = km(igrid,ibloc)
         nmaxd2 = im(igrid,ibloc)
      end if
c
      if(imap(4,iseg,ibloc,1) .gt. nmaxd1) then
         write(iwrit,'(2x,"size of seg ",i3," of block ",i3,
     .                 " too large in the 1st direction")') iseg,ibloc
cBWW  added stop statement
         call ERREXIT (nodes)
      end if
      if(imap(6,iseg,ibloc,1) .gt. nmaxd2) then
         write(iwrit,'(2x,"size of seg ",i3," of block ",i3,
     .                 " too large in the 2nd direction")') iseg,ibloc
         call ERREXIT (nodes)
      end if
  200 continue
c
c      check to see if indices for turbulent region exceed the
c      dimensions of the segment
c
      do 250 ibloc=1,nbloc
      do 250 iseg=1,nseg(ibloc)
      if(imap(13,iseg,ibloc,1) .eq. 1) then
        if(imap(15,iseg,ibloc,1) .gt. imap(4,iseg,ibloc,1). or .
     .     imap(14,iseg,ibloc,1) .lt. imap(3,iseg,ibloc,1) ) then
           write(iwrit,'(2x,"turbulent region ",i3," of block ",i3,
     .                 " too large in the 1st direction")') iseg,ibloc
           call ERREXIT (nodes)
        end if
c
        if(imap(17,iseg,ibloc,1) .gt. imap(6,iseg,ibloc,1). or .
     .     imap(16,iseg,ibloc,1) .lt. imap(5,iseg,ibloc,1) ) then
           write(iwrit,'(2x,"turbulent region ",i3," of block ",i3,
     .                 " too large in the 2nd direction")') iseg,ibloc
           call ERREXIT (nodes)
        end if
c
      endif
c
  250 continue
c
c      check to see if each face is completely filled
c
      do 300 ibloc=1,nbloc
c
c      find the min and max values for each face
c
      nmn1(1) =  jmn
      nmx1(1) = jm(igrid,ibloc) +1
      nmn2(1) =  kmn
      nmx2(1) = km(igrid,ibloc) +1
      nmn1(2) =  jmn
      nmx1(2) = jm(igrid,ibloc) +1
      nmn2(2) =  kmn
      nmx2(2) = km(igrid,ibloc) +1
      nmn1(3) =  kmn
      nmx1(3) = km(igrid,ibloc) +1
      nmn2(3) =  imn
      nmx2(3) = im(igrid,ibloc) +1
      nmn1(4) =  kmn
      nmx1(4) = km(igrid,ibloc) +1
      nmn2(4) =  imn
      nmx2(4) = im(igrid,ibloc) +1
      nmn1(5) =  imn
      nmx1(5) = im(igrid,ibloc) +1
      nmn2(5) =  jmn
      nmx2(5) = jm(igrid,ibloc) +1
      nmn1(6) =  imn
      nmx1(6) = im(igrid,ibloc) +1
      nmn2(6) =  jmn
      nmx2(6) = jm(igrid,ibloc) +1
c
c        initialize all face points to -1
c
      do 211 nd1=nmn1(1),nmx1(1)
      do 211 nd2=nmn2(1),nmx2(1)
      wk2da(nd1,nd2) = -1.
  211 continue
c
      do 212 nd1=nmn1(2),nmx1(2)
      do 212 nd2=nmn2(2),nmx2(2)
      wk2db(nd1,nd2) = -1.
  212 continue
c
      do 213 nd1=nmn1(3),nmx1(3)
      do 213 nd2=nmn2(3),nmx2(3)
      wk2dc(nd1,nd2) = -1.
  213 continue
c
      do 214 nd1=nmn1(4),nmx1(4)
      do 214 nd2=nmn2(4),nmx2(4)
      wk2dd(nd1,nd2) = -1.
  214 continue
c
      do 215 nd1=nmn1(5),nmx1(5)
      do 215 nd2=nmn2(5),nmx2(5)
      wk2de(nd1,nd2) = -1.
  215 continue
c
      do 216 nd1=nmn1(6),nmx1(6)
      do 216 nd2=nmn2(6),nmx2(6)
      wk2df(nd1,nd2) = -1.
  216 continue
c
c        fill in face points for each segment
c
      do 220 iseg=1,nseg(ibloc)
      nf     = imap(2,iseg,ibloc,1)
      nstrt1 = imap(3,iseg,ibloc,1)
      nstop1 = imap(4,iseg,ibloc,1) +1
      nstrt2 = imap(5,iseg,ibloc,1)
      nstop2 = imap(6,iseg,ibloc,1) +1
c
c     if(nstrt1.ne.nmn1(nf)) nstrt1 = nstrt1 +1
c     if(nstop1.ne.nmx1(nf)) nstop1 = nstop1 -1
c     if(nstrt2.ne.nmn2(nf)) nstrt2 = nstrt2 +1
c     if(nstop2.ne.nmx2(nf)) nstop2 = nstop2 -1
c---
c     adjust segment dimensions to affect interior portion (vatsa 10-3-94)
       nstrt1 = nstrt1 +1
       nstop1 = nstop1 -1
       nstrt2 = nstrt2 +1
       nstop2 = nstop2 -1
c
      do 220 nd1=nstrt1,nstop1
      do 220 nd2=nstrt2,nstop2
      if(nf .eq. 1) wk2da(nd1,nd2) = 1.
      if(nf .eq. 2) wk2db(nd1,nd2) = 1.
      if(nf .eq. 3) wk2dc(nd1,nd2) = 1.
      if(nf .eq. 4) wk2dd(nd1,nd2) = 1.
      if(nf .eq. 5) wk2de(nd1,nd2) = 1.
      if(nf .eq. 6) wk2df(nd1,nd2) = 1.
  220 continue
c
c        check to see if any face points have not been set
c
c     adjust segment dimensions to affect interior portion (vatsa 10-3-94)
c     do 231 nd1=nmn1(1),nmx1(1)
c     do 231 nd2=nmn2(1),nmx2(1)
      do 231 nd1=nmn1(1)+1,nmx1(1)-1
      do 231 nd2=nmn2(1)+1,nmx2(1)-1
      if(wk2da(nd1,nd2) .lt. 0.) then
         write(iwrit,'(2x,"boundary info for point (",i3,",",i3,
     .                    ") of imin face not defined for block",
     .                    i3)') nd1,nd2,ibloc
         call ERREXIT (nodes)
      end if
  231 continue
c
c
c     adjust segment dimensions to affect interior portion (vatsa 10-3-94)
c     do 232 nd1=nmn1(2),nmx1(2)
c     do 232 nd2=nmn2(2),nmx2(2)
      do 232 nd1=nmn1(2)+1,nmx1(2)-1
      do 232 nd2=nmn2(2)+1,nmx2(2)-1
      if(wk2db(nd1,nd2) .lt. 0.) then
         write(iwrit,'(2x,"boundary info for point (",i3,",",i3,
     .                    ") of imax face not defined for block",
     .                    i3)') nd1,nd2,ibloc
         call ERREXIT (nodes)
      end if
  232 continue
c
c
c     adjust segment dimensions to affect interior portion (vatsa 10-3-94)
c     do 233 nd1=nmn1(3),nmx1(3)
c     do 233 nd2=nmn2(3),nmx2(3)
      do 233 nd1=nmn1(3)+1,nmx1(3)-1
      do 233 nd2=nmn2(3)+1,nmx2(3)-1
      if(wk2dc(nd1,nd2) .lt. 0.) then
         write(iwrit,'(2x,"boundary info for point (",i3,",",i3,
     .                    ") of jmin face not defined for block",
     .                    i3)') nd1,nd2,ibloc
         call ERREXIT (nodes)
      end if
  233 continue
c
c
c     adjust segment dimensions to affect interior portion (vatsa 10-3-94)
c     do 234 nd1=nmn1(4),nmx1(4)
c     do 234 nd2=nmn2(4),nmx2(4)
      do 234 nd1=nmn1(4)+1,nmx1(4)-1
      do 234 nd2=nmn2(4)+1,nmx2(4)-1
      if(wk2dd(nd1,nd2) .lt. 0.) then
         write(iwrit,'(2x,"boundary info for point (",i3,",",i3,
     .                    ") of jmax face not defined for block",
     .                    i3)') nd1,nd2,ibloc
         call ERREXIT (nodes)
      end if
  234 continue
c
c
c     adjust segment dimensions to affect interior portion (vatsa 10-3-94)
c     do 235 nd1=nmn1(5),nmx1(5)
c     do 235 nd2=nmn2(5),nmx2(5)
      do 235 nd1=nmn1(5)+1,nmx1(5)-1
      do 235 nd2=nmn2(5)+1,nmx2(5)-1
      if(wk2de(nd1,nd2) .lt. 0.) then
         write(iwrit,'(2x,"boundary info for point (",i3,",",i3,
     .                    ") of kmin face not defined for block",
     .                    i3)') nd1,nd2,ibloc
         call ERREXIT (nodes)
      end if
  235 continue
c
c
c     adjust segment dimensions to affect interior portion (vatsa 10-3-94)
c     do 236 nd1=nmn1(6),nmx1(6)
c     do 236 nd2=nmn2(6),nmx2(6)
      do 236 nd1=nmn1(6)+1,nmx1(6)-1
      do 236 nd2=nmn2(6)+1,nmx2(6)-1
      if(wk2df(nd1,nd2) .lt. 0.) then
         write(iwrit,'(2x,"boundary info for point (",i3,",",i3,
     .                    ") of kmax face not defined for block",
     .                    i3)') nd1,nd2,ibloc
         call ERREXIT (nodes)
      end if
  236 continue
c
  300 continue
c
c      check to see if each interface matches correctly with its target
c
      do 400 ibloc=1,nbloc
      do 400 iseg=1,nseg(ibloc)
cpch  if(imap(1,iseg,ibloc,1) .le. 1) then
c
cBWW  cut condition
c
      if (imap(1,iseg,ibloc,1).eq.0 .or.
     .    imap(1,iseg,ibloc,1).eq.1) then
         nbs  = imap(7,iseg,ibloc,1)
         nfs  = imap(8,iseg,ibloc,1)
         mns1 = min(imap( 9,iseg,ibloc,1),imap(10,iseg,ibloc,1))
         mxs1 = max(imap( 9,iseg,ibloc,1),imap(10,iseg,ibloc,1))
         mns2 = min(imap(11,iseg,ibloc,1),imap(12,iseg,ibloc,1))
         mxs2 = max(imap(11,iseg,ibloc,1),imap(12,iseg,ibloc,1))
         nfs  = abs(nfs)
c
c        check to see if the source info is correct
c
         match = 0
         do 310 ns=1,nseg(nbs)
         if(imap(1,ns,nbs,1).le.1 .and. imap(2,ns,nbs,1).eq.nfs) then
          if((imap(3,ns,nbs,1).eq.mns1 .and. 
     .        imap(4,ns,nbs,1).eq.mxs1) .and.
     .       (imap(5,ns,nbs,1).eq.mns2 .and.
     .        imap(6,ns,nbs,1).eq.mxs2))     then
             match = 1
          end if
         end if
  310    continue
         if(match .eq. 0) then
            write(iwrit,'(2x,"segment ",i3," of block ",i3,
     .                    " does not match target info")') iseg,ibloc
            call ERREXIT (nodes)
         end if
c
cBWW  wall condition (2 - inviscid, 3 - viscous)
cBWW  check on extyp = imap(7,...) moved here from
cBWW  bcinwfl, bcvwfl and bcvwfl2
c
      else if (imap(1,iseg,ibloc,1).eq.2 .or.
     .         imap(1,iseg,ibloc,1).eq.3) then
        if (imap(7,iseg,ibloc,1).lt.-1 .or.
     .      imap(7,iseg,ibloc,1).gt. 1) then
          write (iwrit,'(2x," extyp =",i5,
     .      " not valid for bctype =",i5)')
     .      imap(7,iseg,ibloc,1),imap(1,iseg,ibloc,1)
          call ERREXIT (nodes)
        endif
c
cBWW  symmetry condition
cBWW  checks on extyp = imap(7,...) and nsym = imap(8,...)
cBWW  moved here from bcsymfl
c
      else if (imap(1,iseg,ibloc,1).eq.4) then
        if (imap(7,iseg,ibloc,1).lt.0 .or.
     .      imap(7,iseg,ibloc,1).gt.1) then
          write (iwrit,'(2x," extyp =",i5,
     .      " not valid for bctype =",i5)')
     .      imap(7,iseg,ibloc,1),imap(1,iseg,ibloc,1)
          call ERREXIT (nodes)
        endif
c
c
cBWW  all other conditions (5 and above)
cBWW  checks on extyp = imap(7,...) moved here from
cBWW  bcextfl,bcfndfl,bcfnufl and bcsinfl
c
      else if (imap(1,iseg,ibloc,1).ge.5) then
        if (imap(7,iseg,ibloc,1).lt.0 .or.
     .      imap(7,iseg,ibloc,1).gt.1) then
          write (iwrit,'(2x," extyp =",i5,
     .      " not valid for bctype = ",i5)')
     .      imap(7,iseg,ibloc,1),imap(1,iseg,ibloc,1)
          call ERREXIT (nodes)
        endif
      end if
c
  400 continue
c
c     set limits for turbulent zone
c
cBWW  changed output format for turbulent zones
cBWW  output now done in main
cBWW  numout = 0
      do 450 ibloc=1,nbloc
c
      imd = im(1,ibloc)
      if( itrb1(1,ibloc).gt.imd.or.itrb2(1,ibloc).gt.imd) then
         write(iwrit,'(2x,"turbulent zone in i-direction",2i5,
     .                    "larger than block size",i5)')
     .                 itrb1(1,ibloc),itrb2(1,ibloc),imd
         call ERREXIT (nodes)
      endif
c
      if(itrb1(1,ibloc).le.0) then
         itrb1(1,ibloc) = imn
      else
         itrb1(1,ibloc) = itrb1(1,ibloc) + imn -1
      endif
c
      if(itrb2(1,ibloc).le.0) then
         if(itrb2(1,ibloc).lt.0) itrb2(1,ibloc) = imn
         if(itrb2(1,ibloc).eq.0) itrb2(1,ibloc) = im(1,ibloc)
      else
         itrb2(1,ibloc) = itrb2(1,ibloc) + imn -1
      endif
c
      jmd = jm(1,ibloc) 
      if( jtrb1(1,ibloc).gt.jmd.or.jtrb2(1,ibloc).gt.jmd) then
         write(iwrit,'(2x,"turbulent zone in j-direction",2i5,
     .                    "larger than block size",i5)')
     .                 jtrb1(1,ibloc),jtrb2(1,ibloc),jmd
         call ERREXIT (nodes)
      endif
c
      if(jtrb1(1,ibloc).le.0) then
         jtrb1(1,ibloc) = jmn
      else
         jtrb1(1,ibloc) = jtrb1(1,ibloc) + jmn -1
      endif
c
      if(jtrb2(1,ibloc).le.0) then
         if(jtrb2(1,ibloc).lt.0) jtrb2(1,ibloc) = jmn
         if(jtrb2(1,ibloc).eq.0) jtrb2(1,ibloc) = jm(1,ibloc)
      else
         jtrb2(1,ibloc) = jtrb2(1,ibloc) + jmn -1
      endif
c
      kmd = km(1,ibloc)
      if(i2dfl.eq.1) then
         ktrb1(1,ibloc) = kmn
         ktrb2(1,ibloc) = km(1,ibloc)
         go to 440
      endif
c
      if( ktrb1(1,ibloc).gt.kmd.or.ktrb2(1,ibloc).gt.kmd) then
         write(iwrit,'(2x,"turbulent zone in k-direction",2i5,
     .                    "larger than block size",i5)')
     .                 ktrb1(1,ibloc),ktrb2(1,ibloc),kmd
         call ERREXIT (nodes)
      endif
c
      if(ktrb1(1,ibloc).le.0) then
         ktrb1(1,ibloc) = kmn
      else
         ktrb1(1,ibloc) = ktrb1(1,ibloc) + kmn -1
      endif
c
      if(ktrb2(1,ibloc).le.0) then
         if(ktrb2(1,ibloc).lt.0) ktrb2(1,ibloc) = kmn
         if(ktrb2(1,ibloc).eq.0) ktrb2(1,ibloc) = km(1,ibloc)
      else
         ktrb2(1,ibloc) = ktrb2(1,ibloc) + kmn -1
      endif
c
 440  continue
c
 450  continue
c
c     set limits for segments for all multigrid levels
c
      do 500 ng=2,ngrid
      do 500 ibloc=1,nbloc
      do 510 iseg =1,nseg(ibloc)
      imap(1,iseg,ibloc,ng) =  imap(1,iseg,ibloc,ng-1)
      imap(2,iseg,ibloc,ng) =  imap(2,iseg,ibloc,ng-1)
      imap(3,iseg,ibloc,ng) = (imap(3,iseg,ibloc,ng-1)-imn)/2 +imn
      imap(4,iseg,ibloc,ng) = (imap(4,iseg,ibloc,ng-1)-imn)/2 +imn
      imap(5,iseg,ibloc,ng) = (imap(5,iseg,ibloc,ng-1)-imn)/2 +imn
      imap(6,iseg,ibloc,ng) = (imap(6,iseg,ibloc,ng-1)-imn)/2 +imn
      imap(7,iseg,ibloc,ng) =  imap(7,iseg,ibloc,ng-1)
      imap(8,iseg,ibloc,ng) =  imap(8,iseg,ibloc,ng-1)
c
      imap(13,iseg,ibloc,ng) =  imap(13,iseg,ibloc,ng-1)
      if(imap(13,iseg,ibloc,ng) .eq. 1) then
         imap(14,iseg,ibloc,ng) = (imap(14,iseg,ibloc,ng-1)-imn)/2 +imn
         imap(15,iseg,ibloc,ng) = (imap(15,iseg,ibloc,ng-1)-imn)/2 +imn
         imap(16,iseg,ibloc,ng) = (imap(16,iseg,ibloc,ng-1)-imn)/2 +imn
         imap(17,iseg,ibloc,ng) = (imap(17,iseg,ibloc,ng-1)-imn)/2 +imn
         imap(18,iseg,ibloc,ng) = (imap(18,iseg,ibloc,ng-1)-imn)/2 +imn
         imap(19,iseg,ibloc,ng) = (imap(19,iseg,ibloc,ng-1)-imn)/2 +imn
         imap(20,iseg,ibloc,ng) = (imap(20,iseg,ibloc,ng-1)-imn)/2 +imn
      end if
c
      if(imap( 1,iseg,ibloc,ng) .le. 1) then
         imap( 9,iseg,ibloc,ng) = (imap( 9,iseg,ibloc,ng-1)-imn)/2 +imn
         imap(10,iseg,ibloc,ng) = (imap(10,iseg,ibloc,ng-1)-imn)/2 +imn
         imap(11,iseg,ibloc,ng) = (imap(11,iseg,ibloc,ng-1)-imn)/2 +imn
         imap(12,iseg,ibloc,ng) = (imap(12,iseg,ibloc,ng-1)-imn)/2 +imn
      end if
c
c     setup for 2-d flows
c
      if(i2dfl.eq.1) then
         if(imap( 2,iseg,ibloc,ng).eq.1.or.imap( 2,iseg,ibloc,ng).eq.2)
     .   then
            imap( 5,iseg,ibloc,ng) = imap( 5,iseg,ibloc,ng-1)
            imap( 6,iseg,ibloc,ng) = imap( 6,iseg,ibloc,ng-1)
c
            if(imap( 1,iseg,ibloc,ng).le.1) then
               if(imap( 8,iseg,ibloc,ng).gt.0) then
                  imap(11,iseg,ibloc,ng) = imap(11,iseg,ibloc,ng-1)
                  imap(12,iseg,ibloc,ng) = imap(12,iseg,ibloc,ng-1)
               else
                  imap( 9,iseg,ibloc,ng) = imap( 9,iseg,ibloc,ng-1)
                  imap(10,iseg,ibloc,ng) = imap(10,iseg,ibloc,ng-1)
               end if
            endif
c
            if(imap(13,iseg,ibloc,ng).eq.1) then
               imap(16,iseg,ibloc,ng) = imap(16,iseg,ibloc,ng-1)
               imap(17,iseg,ibloc,ng) = imap(17,iseg,ibloc,ng-1)
            endif
c
         endif
c
         if(imap( 2,iseg,ibloc,ng).eq.3.or.imap( 2,iseg,ibloc,ng).eq.4)
     .   then
            imap( 3,iseg,ibloc,ng) = imap( 3,iseg,ibloc,ng-1)
            imap( 4,iseg,ibloc,ng) = imap( 4,iseg,ibloc,ng-1)
c
            if(imap( 1,iseg,ibloc,ng).le.1) then
               if(imap( 8,iseg,ibloc,ng).gt.0) then
                  imap( 9,iseg,ibloc,ng) = imap( 9,iseg,ibloc,ng-1)
                  imap(10,iseg,ibloc,ng) = imap(10,iseg,ibloc,ng-1)
               else
                  imap(11,iseg,ibloc,ng) = imap(11,iseg,ibloc,ng-1)
                  imap(12,iseg,ibloc,ng) = imap(12,iseg,ibloc,ng-1)
               end if
            endif
c
            if(imap(13,iseg,ibloc,ng).eq.1) then
               imap(14,iseg,ibloc,ng) = imap(14,iseg,ibloc,ng-1)
               imap(15,iseg,ibloc,ng) = imap(15,iseg,ibloc,ng-1)
            endif
c
         endif
c
      endif
c
  510 continue
c
      itrb1(ng,ibloc) = (itrb1(ng-1,ibloc)-imn)/2 + imn
      itrb2(ng,ibloc) = (itrb2(ng-1,ibloc)-imn)/2 + imn
      jtrb1(ng,ibloc) = (jtrb1(ng-1,ibloc)-jmn)/2 + jmn
      jtrb2(ng,ibloc) = (jtrb2(ng-1,ibloc)-jmn)/2 + jmn
      if(i2dfl.eq.1) then
         ktrb1(ng,ibloc) = ktrb1(ng-1,ibloc)
         ktrb2(ng,ibloc) = ktrb2(ng-1,ibloc)
      else
         ktrb1(ng,ibloc) = (ktrb1(ng-1,ibloc)-kmn)/2 + kmn
         ktrb2(ng,ibloc) = (ktrb2(ng-1,ibloc)-kmn)/2 + kmn
      endif
c
c     check to see if each segment is multigrid-able to ngrid levels
c
      if(i2dfl .ne. 1) then
        ifirst = 0
        do 520 iseg =1,nseg(ibloc)
        n1dimf = abs(imap(4,iseg,ibloc,ng-1) -imap(3,iseg,ibloc,ng-1))
        n2dimf = abs(imap(6,iseg,ibloc,ng-1) -imap(5,iseg,ibloc,ng-1))
        n1dimc = abs(imap(4,iseg,ibloc,ng  ) -imap(3,iseg,ibloc,ng  ))
        n2dimc = abs(imap(6,iseg,ibloc,ng  ) -imap(5,iseg,ibloc,ng  ))
        if(n1dimf .ne. 2*n1dimc .or. n2dimf .ne. 2*n2dimc) then
           if (ifirst .eq. 0) write(iwrit,'(1x//)')
           ifirst = 1
           write(iwrit,'("****WARNING: seg ",i3," of block ",i3," is "
     .                   "not dimensioned for ",i3," levels of "
     .                   "multigrid***")') iseg,ibloc,ng
        end if
  520   continue
      else
        ifirst = 0
        do 530 iseg =1,nseg(ibloc)
        nface  =     imap(2,iseg,ibloc,ng)
        if(nface .le. 2) then
          n1dimf = abs(imap(4,iseg,ibloc,ng-1) -imap(3,iseg,ibloc,ng-1))
          n1dimc = abs(imap(4,iseg,ibloc,ng  ) -imap(3,iseg,ibloc,ng  ))
          if(n1dimf .ne. 2*n1dimc) then
             if (ifirst .eq. 0) write(iwrit,'(1x//)')
             ifirst = 1
             write(iwrit,'("****WARNING: seg ",i3," of block ",i3," is "
     .                     "not dimensioned for ",i3," levels of "
     .                     "multigrid***")') iseg,ibloc,ng
          end if
        else if(nface .gt. 2 .and. nface .lt. 5) then
          n2dimf = abs(imap(6,iseg,ibloc,ng-1) -imap(5,iseg,ibloc,ng-1))
          n2dimc = abs(imap(6,iseg,ibloc,ng  ) -imap(5,iseg,ibloc,ng  ))
          if(n2dimf .ne. 2*n2dimc) then
             if (ifirst .eq. 0) write(iwrit,'(1x//)')
             ifirst = 1
             write(iwrit,'("****WARNING: seg ",i3," of block ",i3," is "
     .                     "not dimensioned for ",i3," levels of "
     .                     "multigrid***")') iseg,ibloc,ng
          end if
        end if
  530   continue
      end if
c
  500 continue
c
cBWW     write(iwrit,'(1x/)')
c
c
      return
      end
