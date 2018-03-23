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
c
      subroutine sizer_g0(ndb,nbloc,ndmg,ngrlev,ic,jc,kc,id1,id2,
     .               idum,nnodes,nodes,iwrk,npchc,npchitm,npchits,
     .               mxid,mxjd,mxkd,mx1d,mxsizc,mxsizn,mxszcg,mxszfg,
     .               mxs2dc,mxsurf,mxpchs,mxtpchc,mpchitm,mxpchc,
     .               mxpitm,ipatchg,mxijkp3,mx2dwk,nbufct,nbufsz,i2dfl)
c

c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      dimension ic(ndb),jc(ndb),kc(ndb),id1(ndb),id2(ndb)
      dimension idum(21),nodes(ndb),iwrk(3,ndb)
      dimension npchc(ndb,ndmg),npchitm(ndb,ndmg),npchits(ndb,ndmg)
c
c---- set up block to process mapping
c
c
c     find the largest single block dimension sizes for 3-d temp arrays
c     (mods by : v.n.vatsa : 5-29-97 )
c
      ijkmax  = 0
      do n=1,nbloc
        ijkp3loc = (ic(n)+3)*(jc(n)+3)*(kc(n)+3)
        if (ijkp3loc.gt.ijkmax) then
          ijkmax = ijkp3loc
        endif
      enddo
      nmxijkp3 = ijkmax
c
c     find the largest single block dimension sizes for 2-d temp arrays
c     (mods by : v.n.vatsa : 5-29-97 )
c
      ijmax  = 0
      do n=1,nbloc
         ijmx1 = (ic(n)+3)*(jc(n)+3)
         ijmx2 = (jc(n)+3)*(kc(n)+3)
         ijmx  = max(ijmx1,ijmx2)
         ijmx3 = (kc(n)+3)*(ic(n)+3)
         ijmxloc = max(ijmx,ijmx3) 
         if (ijmxloc.gt.ijmax) then
          ijmax = ijmxloc
         endif
      enddo
      nmx2dwk = ijmax
c
c
c     read dimension information from map file 
c
      npsurf = 0
      rewind (imapfil)
      npsurf = impsrf(imapfil,maxnseg,msegmch,msegpch)
      close(imapfil)
c
      if (nnodes.le.0 .or. nnodes.gt.nbloc)
     .  nnodes = nbloc
      call mapblks (1,ndb,nbloc,ic,jc,kc,
     .              nnodes,nodes,iwrk,65,iwrit,ierr)
      if (ierr.ne.0) stop
c
      do inode=1,nnodes
        iwrk(1,inode) = 0
        iwrk(2,inode) = 0
      enddo
      do nb=1,nbloc
        inode         = nodes(nb)
        iwrk(1,inode) = iwrk(1,inode) + 1
        iwrk(2,inode) = iwrk(2,inode) + ic(nb)*jc(nb)*kc(nb)
      enddo
      nodemax = 0
      nodemxa = 0
      do inode=1,nnodes
        if (nodemax.lt.iwrk(2,inode))  nodemax  = iwrk(2,inode)
        if (nodemxa.lt.iwrk(1,inode))  nodemxa  = iwrk(1,inode)
      enddo
c
c---- work with patched grid files and get max dimensions
c
      nxtpchc = 1
      nxpchs  = 1
      nxtpitm = 1
      nxtpits = 1
c
c       initialize npchc array
c
        do igrid=1,ndmg
          do ib   = 1,ndb
             npchc(ib,igrid) = 0
          enddo
        enddo
c
c       initialize npchitm array
c
        do igrid=1,ndmg
          do ib   = 1,ndb
             npchitm(ib,igrid) = 0
             npchits(ib,igrid) = 0
          enddo
        enddo
c
      if (ipatchg.gt.0) then
c
c       read first patch file
c
        rewind ipfil1
        read (ipfil1,*,end=4) ngridp
c
c
        do igrid=1,ngridp
          read (ipfil1,'(1x)')
          lpchcb = 0
          do ibloc=1,nbloc
            do n=1,7
              read (ipfil1,'(1x)')
            enddo
            read (ipfil1,*) (idum(n),n=1,21)
            do n=10,15
              if (idum(n).gt.0) then
                do nn=1,idum(n)
                  read (ipfil1,*) iptch,jptch,kptch
                enddo
                lpchcb = lpchcb + idum(n)
                npchc(ibloc,igrid) = npchc(ibloc,igrid)+idum(n)
              endif
            enddo
          enddo
cvn       if (nxtpchc.lt.lpchcb) nxtpchc = lpchcb
        enddo
        close (ipfil1)
c
c       read second patch file
c
        rewind ipfil2
        read (ipfil2,*,end=4) ngridp
c
c
cfix for patched asynchronous array sizes (vatsa: nov. 98)
      mxitms  = 1
c
        do igrid=1,ngridp
          read (ipfil2,'(1x)')
          read (ipfil2,*) nswpcht
          if (nxpchs.lt.nswpcht) nxpchs = nswpcht
cvn       litmnum = 0
          do nn=1,nswpcht
            read (ipfil2,'(1x)')
            read (ipfil2,'(1x)')
            read (ipfil2,'(1x)')
            read (ipfil2,*) (idum(n),n=1,4),nitems
cfix for patched asynchronous array sizes (vatsa: nov. 98)
            mxitms  = max(mxitms,nitems)
c
            iblock1 = idum(1)
            iblock2 = idum(3)
            npchitm(iblock1,igrid) = npchitm(iblock1,igrid) + nitems
            npchits(iblock2,igrid) = npchits(iblock2,igrid) + nitems
cvn         litmnum = litmnum + nitems
            do item=1,nitems
              read (ipfil2,*) (idum(n),n=1,6),dum
            enddo
          enddo
cvn       if (npchitm.lt.litmnum) npchitm = litmnum

        enddo
        close (ipfil2)

        go to 5
    4   write (iwrit,'("patch file not found or empty")')
        stop
      endif
c
c---- print node assignments and grid block dimensions
c
    5 continue
      write (iwrit,'(/,"Node Assignments")')
      write (iwrit,'("node  blocks  points  pts/max")')
      write (iwrit,'("----  ------  ------  -------")')
      do inode=1,nnodes
        write (iwrit,'(i4,i6,i10,f9.4)') inode,iwrk(1,inode),
     .    iwrk(2,inode),real(iwrk(2,inode))/real(nodemax)
      enddo
c
      write (iwrit,'(/,"Grid Block Sizes")')
      write (iwrit,'("block  node   imax   jmax   kmax     total")')
      write (iwrit,'("-----  ----   ----   ----   ----     -----")')
c
      do 45 nb=1,nbloc
      write (iwrit,'(i5,i6,3i7,i10)') nb,nodes(nb),
     .  ic(nb),jc(nb),kc(nb),ic(nb)*jc(nb)*kc(nb)
      ic(nb) = ic(nb) -1
      jc(nb) = jc(nb) -1
      kc(nb) = kc(nb) -1
c
      id1(nb) = max(ic(nb),jc(nb),kc(nb))
      if(id1(nb).eq.ic(nb)) then
        id2(nb) = max(jc(nb),kc(nb))
      else if(id1(nb).eq.jc(nb)) then
        id2(nb) = max(ic(nb),kc(nb))
      else
        id2(nb) = max(ic(nb),jc(nb))
      end if
   45 continue
c
c     find maximum dimensions in i,j,k across all blocks
c
      nxid   = ic(1)
      nxjd   = jc(1)
      nxkd   = kc(1)
      do 65 n=2,nbloc
      if(ic(n).gt.nxid) nxid = ic(n)
      if(jc(n).gt.nxjd) nxjd = jc(n)
      if(kc(n).gt.nxkd) nxkd = kc(n)
 65   continue
      nxid   = nxid + 1
      nxjd   = nxjd + 1
      nxkd   = nxkd + 1
      nx1d   = max(nxid,nxjd,nxkd)
c
cold  i2dfl = 0
cold  if(kc(1).le.2) i2dfl = 1
c
c---- find max dimensions across all nodes
c
      nxsizc = 1
      nxsizn = 1
      nxs2dc = 1
      nxszfg = 1
cvnv  --- initialize nxszcg -- (by vatsa 7-18-95)
      nxszcg = 1
cvnew dimensions for distributed patched arrays (3-3-97)
      nxtpchc = 1
      nxtpitm = 1
      nxtpits = 1
c
      do 120 inode=1,nnodes
      lxsizc = 1
      lxsizn = 1
      lxs2dc = 1
      lxszfg = 1
cvnew (5-10-96)
      lxszcg = 1
cvnew dimensions for distributed patched arrays (3-3-97)
      lxszpcc = 1
      lxszitm = 1
      lxszits = 1
c
c     mods to identify max. values of patch cells/items for any block
c     (mods by vatsa: Jan. 99)
c
      lmpchc  = 1
      lmpitms = 1
c
      do 110 ng=1,ngrlev
      do 100 nb=1,nbloc
c
c
        lmpchc  = max(lmpchc,npchc(nb,ng) )
        lmpitms = max(lmpitms,npchitm(nb,ng) )
c
c     only blocks assigned to this node
c
      if (nodes(nb).eq.inode) then
c
        lxsizc = lxsizc +(ic(nb) +4)*(jc(nb) +4)*(kc(nb) +4)
        lxsizn = lxsizn +(ic(nb) +3)*(jc(nb) +3)*(kc(nb) +3)
        lxs2dc = lxs2dc +(id1(nb)+3)*(id2(nb)+3)
cvnew (5-10-96)
        if(ng.gt.1)
     .  lxszcg = lxszcg +(ic(nb) +4)*(jc(nb) +4)*(kc(nb) +4)
c
        lxszpcc = lxszpcc + npchc(nb,ng)
        lxszitm = lxszitm + npchitm(nb,ng)
        lxszits = lxszits + npchits(nb,ng)
        
c
        ic(nb)  = ic(nb)/2
        jc(nb)  = jc(nb)/2
        if (i2dfl.ne.1) kc(nb) = kc(nb)/2
        id1(nb) = id1(nb)/2
        id2(nb) = id2(nb)/2
      endif
  100 continue
c
      if (ng.eq.1) lxszfg = lxsizc
  110 continue
c
      if (nxsizc.lt.lxsizc) nxsizc = lxsizc
      if (nxsizn.lt.lxsizn) nxsizn = lxsizn
      if (nxs2dc.lt.lxs2dc) nxs2dc = lxs2dc
      if (nxszfg.lt.lxszfg) nxszfg = lxszfg
cvnv  --- (mod by vatsa 7-18-95)
cvnv  nxszcg = nxsizc - nxszfg +1
cvnew (5-10-96)
cvnew nxszcg = max (nxszcg,(nxsizc - lxszfg +1))
      if (nxszcg.lt.lxszcg) nxszcg = lxszcg
c
cvnew dimensions for distributed patched arrays (3-3-97)
c
      if(nxtpchc.lt.lxszpcc) nxtpchc = lxszpcc
      if(nxtpitm.lt.lxszitm) nxtpitm = lxszitm
      if(nxtpits.lt.lxszits) nxtpits = lxszits
c
  120 continue
c
      if(nxtpits.gt.nxtpitm) nxtpitm = nxtpits
c
c---- gather the sizing information computed so far
c
      icnt    = 0

      mxid    =   nxid
      mxjd    =   nxjd
      mxkd    =   nxkd
      mx1d    =   nx1d
      mxsizc  =   nxsizc
      mxsizn  =   nxsizn
      mxszcg  =   nxszcg
      mxszfg  =   nxszfg
      mxs2dc  =   nxs2dc
      mxsurf  =   npsurf
      mxpchs  =   nxpchs
      mxtpchc =   nxtpchc
      mpchitm =   nxtpitm
      mxpchc  =   lmpchc
      mxpitm  =   lmpitms
      mxijkp3 =   nmxijkp3
c
c     fix the 2-d arry size to avoid problems for patched cases
c     (4-25-00)
c     mx2dwk  =  (mx1d+3)**2
      mx2dwk1 =  (mx1d+3)**2
      mx2dwk2 =   12*mpchitm/15+1
      mx2dwk  =   max(mx2dwk1,mx2dwk2)

c
c 011497 - Added ASYNC communication for bcflow/bcturb


      nbufct = 1
      nbufsz = 1
      if (nnodes.ge.4) then
cvnv     nbufct=(maxnseg*nbloc)/nnodes
cvnv     nbufct= maxnseg*int((real(nbloc)/real(nnodes)+0.999))
cvnv     nbufct= maxnseg*nodemxa
c98oct nbufsz=max((mxid+3)*(mxjd+3),(mxid+3)*(mxkd+3),(mxjd+3)*(mxkd+3))
cOct00
cvnv  nbufsz=max((nxid+3)*(nxjd+3),(nxid+3)*(nxkd+3),(nxjd+3)*(nxkd+3))
      nbufct  = max(mxsegmch,mxsegpch)*nodemxa
      nbufsz1 = mxsegmch*nmx2dwk
      nbufsz2 = mxsegpch*mxitms
         write(iwrit,'("maxnseg,nbloc,nnodes,nbufct,nodemxa",8i9)')
     .                  maxnseg,nbloc,nnodes,nbufct,nodemxa,
     .                  mxitms,nbufsz1,nbufsz2
         nbufsz = max(nbufsz1,nbufsz2)
         nbufsz=(nodemxa*nbufsz*10)+1
      endif
c
c
      return
      end
c--------------------------------------------------------------------
      subroutine prtprm (iwrit,name,mval,nval,icnt)
c--------------------------------------------------------------------
      character*7 name
c
      if (mval.eq.nval) then
        write (iwrit,'(a7,2x,i10,i10,3x,"SAME")') name,mval,nval
      else if (mval.gt.nval) then
        write (iwrit,'(a7,2x,i10,i10,3x,"OK")') name,mval,nval
      else
        write (iwrit,'(a7,2x,i10,i10,3x,"TOO SMALL")') name,mval,nval
        icnt = icnt + 1
      endif
      return
      end
c-----------------------------------------------------------------------
      integer function impsrf(imapfil,maxnseg,msegmch,msegpch)
c-----------------------------------------------------------------------
c
c     object : to read in and set up the topological mapping 
c              funcions, 'imap'
c
c     initial coding : by m.d.sanetrik (August, 1990)
c     latest  mods   : by vatsa (Oct, 2000) for finding
c                      max. no. of segments with point matched
c                      and patched b.c.'s (msegmch, msegpch)
c
c
c
      dimension     imap(55)
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
c                 = 0  corresponds to an inner cut that maps to same block
c                 = 1  corresponds to a cut that maps to another block
c                 = 2  corresponds to a slip (inviscid) wall
c                 = 3  corresponds to a noslip (viscous) wall
c                 = 4  symmetry condition (imap(8) tells about which plane)
c                 = 5  downstream/extrapolation condition
c                 = 6  far-field condition (Riemann invariants)
c                 = 7  singular line
c                 = 8  fan-upstream engine boundary
c                 = 9  fan-downstream engine boundary
c
c***************************************************************************
c
c     ivisb(1)    = 1  viscous fluxes evaluated along i-direction
c     ivisb(2)    = 1  viscous fluxes evaluated along j-direction
c     ivisb(3)    = 1  viscous fluxes evaluated along k-direction
c
c***************************************************************************
c
      maxnseg = 0
      npsurf  = 0
      msegmch = 1
      msegpch = 1
c
      read(imapfil,'(1x)')
      read(imapfil,*) nbloc
c
c
c----------  read in imap array and shift due to ghost cells  -------------
c
      do 100 ibloc = 1,nbloc
      read(imapfil,'(1x)')
      read(imapfil,'(1x)')
      read(imapfil,*) nseg
      if (nseg .gt. maxnseg) maxnseg = nseg
      read(imapfil,'(1x)')
      read(imapfil,'(1x)')
      lsegmch  = 0
      lsegpch  = 0
      do  99 iseg  = 1,nseg
c
      read(imapfil,'(1x)')
      read(imapfil,*) idum,idum,(imap(n),n= 1,12)
      read(imapfil,*)           (imap(n),n=13,20)
      read(imapfil,'(1x)')
c
      if (imap(1).eq.2 .or. imap(1).eq.3. or. imap(1).eq.102 ) then
         ni = 1 + (imap(4)-imap(3))
         nj = 1 + (imap(6)-imap(5))
         npsurf = npsurf + ni*nj
      endif
      if (imap(1).eq.0 .or. imap(1).eq.1 ) lsegmch = lsegmch + 1
      if (imap(1).eq.-1 ) lsegpch = lsegpch + 1
c
   99 continue
      msegmch = max (msegmch,lsegmch)
      msegpch = max (msegpch,lsegpch)
c
  100 continue
c
      impsrf = npsurf
c
      return
      end
