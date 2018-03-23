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
      program tlns3d_driver
c
c     calculates size of workspace needed for a given
c     run of the distributed version of TLNS3D-MB
c

#if defined(BUILD_PVM) || defined(BUILD_MPI)
c
c---- distributed computing related ----------------------------
c

#  ifdef BUILD_PVM
#     include "fpvm3.h"
#     ifdef TLN_DBL
#       define RTYPE REAL8
#     else
#       define RTYPE REAL4
#     endif
#  else
      include "mpif.h"
#     ifdef TLN_DBL
#       define RTYPE MPI_DOUBLE_PRECISION
#     else
#       define RTYPE MPI_REAL
#     endif
      dimension istat(MPI_STATUS_SIZE)
#  endif
#     include "tags.h"
c
#else
c
c     need these defined for in-core version
c
      integer TAG_ERROR,TAG_OK
      parameter (TAG_ERROR=0,TAG_OK=1)
#endif

c
cNOTE:The variable nodes() must be the last variable in the common
c     block dstcom, since in the subroutines in which it is used
c     the size is unknown and declared as nodes(1). Only in the main
c     routine is the size given as nodes(mbloc*2). This means if you
c     want to add something to the common block, add it BEFORE nodes.
c
cSep. 99 : In the dynamic version of code nodes will be passed through
c          argument list of subroutines instead of common dstcom
c
      character*128  errmsg
      common/dstcom/ errmsg,master,myrank,mycomm,ndlist,nnodes
cdyn .               errmsg,nodes(mbloc*2)
cdyn  nodes will be passed through argument list
c
c***************************************************************
c     data statements to initialize selected variables
c
      data          imn/2/, jmn/2/, kmn/2/
      data mres/1350000/, msegt/20/, msegn/50/
cvn   data mres/350000/, msegt/20/, msegn/50/
c     data outnam/'parms.dc'/,outasync/'async.dc'/
c****************************************************************
c
c     data statements for R-K schemes
c
      dimension    idumrk(12)
      data         idumrk/4,4,4,6,6,6,5,8,7,7,6,6/
      data         mrkstg/8/
c*****************************************************************
c 
      character*80 title,grdfil,rstfil,pltfil,savfil,mapfil,
     .             pchfil1,pchfil2,rhtall,rhtphy,
     .             gpltfil,qpltfil,datfil,pntfil,avgfil
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      dimension idum(21)
      dimension iparamg(35)

      dimension nodes(1), iwrk(1)
      dimension npchc(1),npchitm(1),npchits(1)
      dimension ic(1),jc(1),kc(1),id1(1),id2(1)
      pointer (i_ic,ic),(i_jc,jc),(i_kc,kc),(i_id1,id1),(i_id2,id2)
      pointer (i_nodes,nodes),(i_iwrk,iwrk)
      pointer (i_npchc,npchc),(i_npchitm,npchitm),(i_npchits,npchits)
c
c---- unit assignments
c
c     igrdfil = unit for grid input
c     irstfil = unit for restart file input
c     isavfil = unit for solution output
c     iread   = unit for input data
c     iwrit   = unit for print output
c     ipltfil = unit for plot output (used for R-K error estimates)
c     imapfil = unit for map file input
c
c     additional files for time-dependent version
c     irhtall = unit fot complete residual history at all time/pseudo times
c     irhtphy = unit for residual history at physical time steps
c     igpfil  = unit for storing plot3d grid for time-dependent solution
c               (all node points at selected physical time-steps)
c     iqpfil  = unit for storing plot3d solution for time-dependent solution
c               (all node points at selected physical time-steps)
c     idatfil = unit for storing x,y,z,w(1-5),p
c               (selected node points at all physical time-steps)
c     ipntfil = unit for specifying at which ibloc,i,j,k the info. is
c               to be written on idatfil (ngroup set of values)
c     iavgfil = unit for writing average flow quantities
c
      igrdfil    =  2 
      irstfil    =  3 
      isavfil    =  4 
      input      =  5
      iread      =  9
      iwrit      =  6
      ipltfil    =  7 
      imapfil    =  8 
      iparm      = 20
      ipfil1     = 21
      ipfil2     = 22
c
c     time-dependent related files
c
      irhtall    = 23
      irhtphy    = 24
      igpfil     = 25
      iqpfil     = 26
      idatfil    = 27
      ipntfil    = 28
      iavgfil    = 29
      
c
      nnodes     = 0
c
#if defined(BUILD_PVM) || defined(BUILD_MPI)
#  ifdef BUILD_PVM
c
c---- PVM hook
c
        write (iwrit,'("PVM not supported")')
        stop
#  else
c
c---- MPI initialization
c
      master = 0
      call MPI_Init (ierr)
      if (ierr.ne.MPI_SUCCESS) then
        write (iwrit,'("failed to initialize MPI")')
        stop
      endif
      mycomm = MPI_COMM_WORLD
      call MPI_Comm_rank (mycomm,myrank,ierr)
      call MPI_Comm_size (mycomm,numprocs,ierr)
cvn   if (numprocs.lt.2) then
      if (numprocs.lt.1) then
        write (iwrit,'("no node processes started")')
        call MPI_Finalize (ierr)
        stop
      endif
cjm      nnodes = numprocs - 1
      nnodes = numprocs

c     call xgetpid (mypid)
c     write (iwrit, 6105) myrank, mypid
c6105  format ("rank ", i2, ": mypid = ", i8)
c     call flush (iwrit)
c     call sleep (300)
c
#  endif
#else
c
c---- IN-CORE initialization
c
      master = 0
      myrank = 0
      nnodes = 1

#endif


c####################################################################
c master process (rank 0)
c     read input and send to nodes
c####################################################################
c
      if (myrank.eq.master) then ! begin if 1001
cif1001 begin
c
c
      call pinput (input,iread,'#',.true.,.false.,.false.)
c
      read (iread,'(a80)') title
      len = istlen (title,80)
      write (iwrit,'(a)') title(1:len)
      read (iread,'(1x)')
      read (iread,'(8f10.6)') grdfmt,fstrt,fsave,fplot,fcplpr
c     (for unsteday flows, fstrt=2 starts from steady-state solution)
      read (iread,'(1x)')
c
c     ftorder : order for temporal accuracy (0 for steady-state)
c     dtphy   : physical time-step 
c     frkcase : identifier no. for the R-K scheme picked in rkphil
c             : identifier for MEBDF scheme if frkcase < 0
c
      read (iread,'(8f10.6)') flev,fgrids,fmgtyp,ftorder,dtphy,frkcase
      igrdfmt = grdfmt
      nlev    = flev
      if (nlev.le.0) nlev = 1
      ntorder = abs(ftorder)
      ibdf2opt= 0
      if(ftorder.lt.0) ibdf2opt = 1
c
c     ibdf2opt will run as a BDF3 with different coefficients
c
c     set the mebdf flag if frkcase < 0
c     
      if(frkcase.lt.0.0) then
        mebdf   = 1
        irkcase = 0
      else
        irkcase = frkcase
        mebdf   = 0
      endif
c
c     set the no. of R-K stages for time-accurate scheme
c
      nrkstg  = 1
      if (irkcase.ge.1.and.irkcase.le.12) nrkstg = idumrk(irkcase)
      if (irkcase.gt.12) then
        write (iwrit,'(2x,"irkcase.gt.12 not supported",/)')
        write (iwrit,'(2x,"irkcase = ",i5)') irkcase
        stop 1002
      endif
c
c     set the no. of stgaes for MEBDF4 scheme - use infrastructure of R-K
c
      if (mebdf.eq.1) nrkstg  = 3
c

      if (nrkstg.gt.mrkstg) then
        write (iwrit,'(2x,"nrkstg.gt.mrkstg not supported",/)')
        write (iwrit,'(2x,"nrkstg, mrkstg = ",2i5)') nrkstg,mrkstg
        stop 1002
      endif
      if(irkcase.ge.1) ntorder = 1

      mtorder = max(ntorder,1)
      iunsteady = 0
      if(ntorder.gt.0) iunsteady = 1
      if(ntorder.ge.1.and.nlev.gt.1) then
        write (iwrit,'(2x,"FMG not supported for unsteady flow",/)')
        write (iwrit,'(2x,"nlev = ",i5,2x,"ntorder = ",i5)')
     .                     nlev,ntorder
        stop 1001
      endif
      if(ntorder.gt.6) then
        write (iwrit,'(2x,"ntorder.gt.6 not supported",/)')
        write (iwrit,'(2x,"ntorder = ",i5)') ntorder
        stop 1001
      endif

      ngrlev  = fgrids
c
c     set max no. of grid levels based on the input value of fgrids
c
      mgrlev  = ngrlev
c
c     skip to patch grid flag
c
      do n=1,nlev+18
        read (iread,'(1x)')
      enddo
c
      read (iread,'(8f10.6)') dum1,dum1,dum1,dum1,dum1,fi2dfl
      i2dfl = fi2dfl
      read (iread,'(1x)')
c
      read (iread,'(8f10.6)') gamma,pr,prt,suthc,tref,fiturb,fpatchg
      iturv1  = 0
      iturv2  = 0
      iturb   = fiturb
      if(iturb.ge.2) iturv1 = 1
      if(iturb.eq.3) iturv2 = 1
      ipatchg = fpatchg
      if(ntorder.gt.3.and.iturv1.gt.0) then
        write (iwrit,'(2x,"ntorder.gt.3 unsupported for turb flows",/)')
        write (iwrit,'(2x,"ntorder = ",i5, "iturb = ",i5 )') 
     .  ntorder,iturb
        stop 1001
      endif
c
      do n=1,4
        read (iread,'(1x)')
      enddo
c
c     accomodate additional input related to implicit smoothing
      do n=1,4
        read (iread,'(1x)')
      enddo
c
c     read grid file name
c
      read (iread,'(a80)') grdfil
      call fixstr (grdfil,80)
      len = istlen (grdfil,80)
c
c     get block grid dimensions
c
      if (igrdfmt.eq.1) then
        write (iwrit,'("formatted grid from ",a)') grdfil(1:len)
        open (igrdfil,file=grdfil(1:len),form='formatted')
        rewind (igrdfil)
        read (igrdfil,*,end=2) nbloc
      else
        write (iwrit,'("unformatted grid from ",a)') grdfil(1:len)
        open (igrdfil,file=grdfil(1:len),form='unformatted')
        rewind (igrdfil)
        read (igrdfil,end=2) nbloc
      endif
c
c       set the mbloc parameter based on grid file
c
        mbloc   = nbloc
      if (nnodes.le.0 .or. nnodes.gt.nbloc)
     .  nnodes = nbloc
c
c     allocate memory to grid dimensions (ic,jc,kc)
c
c     iname will be incremented on every call to umalloc
c     icumsiz is a cumulative count of storage allocated in bytes
c
      iname    = 10
      icumsiz = 0
c
      call umalloc (i_ic,mbloc,1,iname,icumsiz)
      call umalloc (i_jc,mbloc,1,iname,icumsiz)
      call umalloc (i_kc,mbloc,1,iname,icumsiz)
      call umalloc (i_id1,mbloc,1,iname,icumsiz)
      call umalloc (i_id2,mbloc,1,iname,icumsiz)
      call umalloc (i_nodes,2*mbloc,1,iname,icumsiz)
c
c     compute parameters for intval and nodeid
c
ctemp  mxint1  = 10*mbloc+2
ctemp  mxint2  = 6*mbloc*mgrlev
ctemp  mxintln = 
ctemp .  ((mxint1/mxint2)*mxint1 + (mxint2/mxint1)*mxint2) /
ctemp .  ((mxint1/mxint2)        + (mxint2/mxint1)       ))
ctemp  call umalloc (i_nodeid,mbloc,1,iname,icumsiz)
ctemp  call umalloc (i_intval,mxintln,1,iname,icumsiz)
      call umalloc (i_iwrk,3*mbloc,1,iname,icumsiz)
      call umalloc (i_npchc,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (i_npchitm,mgrlev*mbloc,1,iname,icumsiz)
      call umalloc (i_npchits,mgrlev*mbloc,1,iname,icumsiz)
c
c
      if (igrdfmt.eq.1) then
        if (i2dfl.eq.1) then
          read (igrdfil,*,end=2) (ic(n),jc(n),n=1,nbloc)
          do n=1,nbloc
cvn         kc(n) = 2
            kc(n) = 3
          enddo
        else
          read (igrdfil,*,end=2) (ic(n),jc(n),kc(n),n=1,nbloc)
        endif
      else
        if (i2dfl.eq.1) then
          read (igrdfil,end=2) (ic(n),jc(n),n=1,nbloc)
          do n=1,nbloc
cvn         kc(n) = 2
            kc(n) = 3
          enddo
        else
          read (igrdfil,end=2) (ic(n),jc(n),kc(n),n=1,nbloc)
        endif
      endif
      close (igrdfil)
      goto 3
c
    2 write (iwrit,'("grid file ",a," not found or EOF")')
     .  grdfil(1:len)
      stop
c
    3 continue
c
c------- read the remaining file names
c
      read  (iread,'(a80)') mapfil
      call fixstr (mapfil,80)
      read  (iread,'(a80)') rstfil
      call fixstr (rstfil,80)
      read  (iread,'(a80)') pltfil
      call fixstr (pltfil,80)
      read  (iread,'(a80)') savfil
      call fixstr (savfil,80)
      read  (iread,'(a80)') pchfil1
      call fixstr (pchfil1,80)
      read  (iread,'(a80)') pchfil2
      call fixstr (pchfil2,80)
c
c     time-dependent related files
c
      read  (iread,'(a80)') rhtall
      call fixstr (rhtall,80)
      read  (iread,'(a80)') rhtphy
      call fixstr (rhtphy,80)
      read  (iread,'(a80)') gpltfil
      call fixstr (gpltfil,80)
      read  (iread,'(a80)') qpltfil
      call fixstr (qpltfil,80)
      read  (iread,'(a80)') datfil
      call fixstr (datfil,80)
      read  (iread,'(a80)') pntfil
      call fixstr (pntfil,80)
      read  (iread,'(a80)') avgfil
      call fixstr (avgfil,80)
c---- map file
c
      len = istlen (mapfil,80)
      open (imapfil,file=mapfil(1:len),form='formatted')
      write (iwrit,'(" topological map file  : ",a)')
     .               mapfil(1:len)
      rewind imapfil
c
c---- patched grid files
c
      if (ipatchg.gt.0) then
        len = istlen (pchfil1,80)
        open (ipfil1,file=pchfil1(1:len),form='formatted')
        write (iwrit,'(" grid patching file 1  : ",a)')
     .                 pchfil1(1:len)
        rewind ipfil1
c
        len = istlen (pchfil2,80)
        open (ipfil2,file=pchfil2(1:len),form='formatted')
        write (iwrit,'(" grid patching file 2  : ",a)')
     .                 pchfil2(1:len)
        rewind ipfil2
      endif
c
c     time-dependent flow related files
c
        len = istlen (rhtall,80)
        open (irhtall,file=rhtall(1:len),form='formatted')
        write (iwrit,'(" complete residual  file : ",a)')
     .                   rhtall(1:len)
        rewind irhtall
c
        len = istlen (rhtphy,80)
        open (irhtphy,file=rhtphy(1:len),form='formatted')
        write (iwrit,'(" physical time residual  file : ",a)')
     .                   rhtphy(1:len)
        rewind irhtphy
c
        len = istlen (gpltfil,80)
        open (igpfil,file=gpltfil(1:len),form='unformatted')
        write (iwrit,'(" plot3d - time dependent grid  file : ",a)')
     .                   gpltfil(1:len)
        rewind igpfil
c
        len = istlen (qpltfil,80)
        open (iqpfil,file=qpltfil(1:len),form='unformatted')
        write (iwrit,'(" plot3d - time dependent solution file: ",a)')
     .                   qpltfil(1:len)
        rewind iqpfil
c
        len = istlen (datfil,80)
        open (idatfil,file=datfil(1:len),form='formatted')
        write (iwrit,'("time dependent data - selected nodes): ",a)')
     .                   datfil(1:len)
        rewind idatfil
c
        len = istlen (pntfil,80)
        open (ipntfil,file=pntfil(1:len),form='formatted')
        write (iwrit,'("node specification file for datfil : ",a)')
     .                   pntfil(1:len)
        rewind ipntfil
c
        len = istlen (avgfil,80)
        open (iavgfil,file=avgfil(1:len),form='unformatted')
        write (iwrit,'("average flow-field file for UNS : ",a)')
     .                   avgfil(1:len)
        rewind iavgfil
c
c     read and store the pointer data
c
      ngroup = 1
      if(ntorder.eq.0) go to 98
      do 95 igrp=1,100000
      read(ipntfil,*,err=96,end=97)
     .    iblk,ibeg,iend,iskip,jbeg,jend,jskip,
     .         kbeg,kend,kskip
   95 continue
   96 write(6,'(2x,"Error reading pointer data for special ",/,
     . "time-dependant output for selected points at each time step",/,
     . 2x,"Check value of parameter ngroup or data on pointer file:",/,
     . a80)' ) ipntfil
   97 continue
      ngroup = igrp - 1
      if(ngroup.le.0) ngroup = 1
      rewind ipntfil
   98 continue
c
c
      call sizer_g0 (mbloc,nbloc,mgrlev,ngrlev,ic,jc,kc,id1,id2,
     .               idum,nnodes,nodes,iwrk,npchc,npchitm,npchits,
     .               mxid,mxjd,mxkd,mx1d,mxsizc,mxsizn,mxszcg,mxszfg,
     .               mxs2dc,mxsurf,mxpchs,mxtpchc,mpchitm,mxpchc,
     .               mxpitm,ipatchg,mxijkp3,mx2dwk,nbufct,nbufsz,i2dfl)
c
c     set the no. of multigrid legs
c
      mmgrlev = max(2,mgrlev)
      mgrleg  = (2**(mmgrlev-2))
c
c     write the parameters
c
      write (iwrit,'("parameters for sizing for nnodes =",//i6/)')nnodes
      write (iwrit,'("mbloc = ",i9,2x,"mgrlev =",i9,2x,"icumsiz =",i9)')
     .                mbloc,mgrlev,icumsiz
      write (iwrit,'("mxid  = ",i9,2x,"mxjd  = ",i9)')
     .                mxid,mxjd
      write (iwrit,'("mxkd  = ",i9,2x,"mx1d  = ",i9)')
     .                mxkd,mx1d
      write (iwrit,'("mxsizc = ",i9,2x,"mxsizn = ",i9)')
     .                mxsizc,mxsizn
      write (iwrit,'("mxszcg = ",i9,2x,"mxszfg = ",i9)')
     .                mxszcg,mxszfg
      write (iwrit,'("mxs2dc = ",i9,2x,"mxsurf = ",i9)')
     .                mxs2dc,mxsurf
      write (iwrit,'("mxpchs = ",i9,2x,"mxtpchc =",i9)')
     .                mxpchs,mxtpchc
      write (iwrit,'("mpchitm =",i9,2x,"mxpchc =",i9,2x,"mxpitm =",i9)')
     .                mpchitm,mxpchc,mxpitm
      write (iwrit,'("mxijkp3 =",i9,2x,"mx2dwk = ",i9)')
     .                mxijkp3,mx2dwk
      write (iwrit,'("nbufct = ",i9,2x,"nbufsz = ",i9)')
     .                nbufct,nbufsz
      write (iwrit,'("mgrleg = ",i9,2x,"mres   = ",i9)')
     .                mgrleg,mres
      write (iwrit,'("msegt  = ",i9,2x,"msegn  = ",i9)')
     .                msegt,msegn
      write (iwrit,'("mtorder =",i9,2x,"ntorder =",i9,2x,"ngroup=",i9)')
     .                mtorder,ntorder,ngroup
      write (iwrit,'("mebdf = "i9,2x,"ibdf2opt =",i9)')  mebdf,ibdf2opt
# ifdef BUILD_MPI
c
c     put all the dimension related parameters in a single array 
c
      iparamg( 1)  = mbloc
      iparamg( 2)  = mgrlev
      iparamg( 3)  = mxid
      iparamg( 4)  = mxjd
      iparamg( 5)  = mxkd
      iparamg( 6)  = mx1d
      iparamg( 7)  = mxsizc
      iparamg( 8)  = mxsizn
      iparamg( 9)  = mxszcg
      iparamg(10)  = mxszfg
      iparamg(11)  = mxs2dc
      iparamg(12)  = mxsurf
      iparamg(13)  = mxpchs
      iparamg(14)  = mxtpchc
      iparamg(15)  = mpchitm
      iparamg(16)  = mxpchc
      iparamg(17)  = mxpitm
      iparamg(18)  = mxijkp3
      iparamg(19)  = mx2dwk
      iparamg(20)  = nbufct
      iparamg(21)  = nbufsz
      iparamg(22)  = mres
      iparamg(23)  = msegt
      iparamg(24)  = msegn
      iparamg(25)  = mgrleg

      iparamg(26)  = mtorder
      iparamg(27)  = ntorder
      iparamg(28)  = ngroup
      iparamg(29)  = iturv1
      iparamg(30)  = iturv2
      iparamg(31)  = iunsteady

      iparamg(32)  = mrkstg
      iparamg(33)  = nrkstg
      iparamg(34)  = mebdf
      iparamg(35)  = ibdf2opt
c
      do inode = 2,nnodes
c
        itag     =   inode -1
c       print *,'before send: inode,TAG,myrank', inode,itag,myrank
        call MPI_Send (iparamg, 35, MPI_INTEGER,
     .                 inode-1,itag,
     .                 mycomm,ierr)
      enddo
c
      else ! correspond to if(myrank.eq.myhost)  statement no.if 1001
c     print *,'start else before mpirecv with myrank',myrank
c
      itag  = myrank
c       print *,'before recv: myrank,TAG', myrank,itag
      call MPI_Recv (iparamg, 35, MPI_INTEGER, master,
     .               itag, mycomm, istat, ierr)

       mbloc   = iparamg( 1)
       mgrlev  = iparamg( 2)
       mxid    = iparamg( 3)
       mxjd    = iparamg( 4)
       mxkd    = iparamg( 5)
       mx1d    = iparamg( 6)
       mxsizc  = iparamg( 7)
       mxsizn  = iparamg( 8)
       mxszcg  = iparamg( 9)
       mxszfg  = iparamg(10)
       mxs2dc  = iparamg(11)
       mxsurf  = iparamg(12)
       mxpchs  = iparamg(13)
       mxtpchc = iparamg(14)
       mpchitm = iparamg(15)
       mxpchc  = iparamg(16)
       mxpitm  = iparamg(17)
       mxijkp3 = iparamg(18)
       mx2dwk  = iparamg(19)
       nbufct  = iparamg(20)
       nbufsz  = iparamg(21)
       mres    = iparamg(22)
       msegt   = iparamg(23)
       msegn   = iparamg(24)
       mgrleg  = iparamg(25)
c
       mtorder = iparamg(26)
       ntorder = iparamg(27)
       ngroup  = iparamg(28)
       iturv1  = iparamg(29)
       iturv2  = iparamg(30)
       iunsteady=iparamg(31)

       mrkstg  = iparamg(32)
       nrkstg  = iparamg(33)
       mebdf   = iparamg(34)
       ibdf2opt= iparamg(35)
# endif
c
       endif ! correspond to if(myrank.eq.myhost)  statement no. 1001
c
c*********************** start logic for all nodes  ************************
c
c     compute some of the secondary parameters needed for dimensioning
c     arrays in the tlns3d program
c
      mxidp3     = mxid+3
      mxjdp3     = mxjd+3
      mxkdp3     = mxkd+3
      mx1dwk     = mx1d+3
      mx3dwk     = mxijkp3
      mxijp3     = mxidp3*mxjdp3
      mx1dwk5    = mx1dwk*5
      mx3dwk3    = mx3dwk*3
      mx3dwk4    = mx3dwk*4
      mx3dwk5    = mx3dwk*5
      mxsiz2c    = mxsizc*2
      mxsiz3c    = mxsizc*3
      mxsiz4c    = mxsizc*4
      mxsiz5c    = mxsizc*5
      mxsz5cg    = mxszcg*5
      mxsz4fg    = mxszfg*4
      mxsz5fg    = mxszfg*5
      mxsiz2n    = mxsizn*2
      mxsiz3n    = mxsizn*3
      mxsiz4n    = mxsizn*4
      mxsiz5n    = mxsizn*5
      mxs2d3c    = mxs2dc*3
c     mxdstf     = mxsz4fg + 14*mxsurf
      mxdsti     = mxszfg  + 12*mxsurf
cvnv (april 01)
cvnv  mxdstf     = mxsz4fg + 15*mxsurf
cvnv (aug. 03)
      mxdstf     = mxsz4fg + 16*mxsurf
c     mxdsti     = mxdstf
cvn ( next 4 lines suggested by Biedron)
c     minbox     = sqrt(float(mxsurf))
c     minbox     = max(minbox,50)
c     nbb        = 3*mxsurf/minbox
c     mxdstf     = mxsz4fg + 9*mxsurf + 7*nbb
c
      write (iwrit,'("mxdsti  = ",i9,2x,"mxdstf  = ",i9)')
     .                mxdsti,mxdstf
cvnv
cvnv  mxinpln    = 15+mgrlev
      mxinpln    = 17+mgrlev
      mxint1     = 10*mbloc+2
      mxint2     = 6*mbloc*mgrlev
      mxintln    =
     .            ((mxint1/mxint2)*mxint1 + (mxint2/mxint1)*mxint2) /
     .            ((mxint1/mxint2)        + (mxint2/mxint1)       )
c
c     rewind and close all the files
c
      rewind igrdfil
      rewind imapfil
      rewind ipfil1
      rewind ipfil2
c
      close (igrdfil)
      close (imapfil)
      close (ipfil1)
      close (ipfil2)
     
     
      call tlns3d (mbloc,mgrlev,mxid,mxjd,mxkd,mx1d,mxsizc,mxsizn,
     .             mxszcg,mxszfg,mxs2dc,mxsurf,mxpchs,mxtpchc,
     .             mpchitm,mxpchc,mxpitm,mx3dwk,mx2dwk,
     .             nbufct,nbufsz,mres,msegt,msegn,mgrleg,
     .             mxidp3,mxjdp3,mxkdp3,mx1dwk,mxijp3,
     .             mx1dwk5,mx3dwk3,mx3dwk4,mx3dwk5,mxsiz2c,
     .             mxsiz3c,mxsiz4c,mxsiz5c,mxsz5cg,mxsz4fg,
     .             mxsz5fg,mxsiz2n,mxsiz3n,mxsiz4n,mxsiz5n,
     .             mxs2d3c,mxdstf,mxdsti,mxinpln,mxint1,mxint2,mxintln,
     .             mtorder,ntorder,iunsteady,ngroup,iturv1,iturv2,
     .             mrkstg,nrkstg,mebdf,ibdf2opt )
c
# ifdef BUILD_MPI
c---- terminate MPI session
c
      call MPI_Finalize (ierr)
#  endif
c
       stop
       end
