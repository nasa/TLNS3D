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
      subroutine force(imn,jmn,kmn,im,jm,km,imp1,jmp1,kmp1,
     .                 imp2,jmp2,kmp2,m1cc,m5cc,m1np,m3np,
     .                 m1is,m1js,m1ks,
     .                 m3is,m3js,m3ks,mgrlev,igrid,mbloc,nbloc,
     .                 mxsizc,mxsiz5c,mxsizn,mxsiz3n,mx1dwk,
     .                 imap,msegt,msegn,nseg,w,p,vol,x,
     .                 six,siy,siz,sjx,sjy,sjz,skx,sky,skz,
     .                 ri1,ri2,rj1,rj2,rk1,rk2,
     .                 tauwfi1,tauwfi2,tauwfj1,tauwfj2,tauwfk1,tauwfk2,
     .                 wk2da,wk2db,wk2dc,wk2dd,
     .                 scal,xref,yref,zref,sref,
     .                 cxg,cyg,czg,cxvg,cyvg,czvg,cmxg,cmyg,cmzg,
     .                 cx,cy,cz,cxv,cyv,czv,cmx,cmy,cmz,ivisb, iforce,
     .                 nodes )
c-----------------------------------------------------------------------
c
c     object : To compute force/moment coefficients on
c              selected blocks/faces
c
c     initial coding : by v.n.vatsa (April, 1990)
c
c     modified for distributed computing : April, 1995
c

      character*128  errmsg
      common/dstcom/ errmsg,master,myrank,mycomm,ndlist,nnodes
      dimension      nodes(1)
c
      dimension     im  (mgrlev,mbloc), jm  (mgrlev,mbloc),
     .              km  (mgrlev,mbloc), imp1(mgrlev,mbloc),
     .              jmp1(mgrlev,mbloc), kmp1(mgrlev,mbloc),
     .              imp2(mgrlev,mbloc), jmp2(mgrlev,mbloc),
     .              kmp2(mgrlev,mbloc),
     .              imap(msegt,msegn,mbloc), nseg(mbloc)
c
      dimension     m1np(mgrlev,mbloc), m3np(mgrlev,mbloc),
     .              m1cc(mgrlev,mbloc), m5cc(mgrlev,mbloc),
     .              m1is(mgrlev,mbloc), m1js(mgrlev,mbloc),
     .              m1ks(mgrlev,mbloc),
     .              m3is(mgrlev,mbloc), m3js(mgrlev,mbloc),
     .              m3ks(mgrlev,mbloc)
c
      dimension     w (mxsiz5c), p(mxsizc ),
     .              vol(mxsizc), x(mxsiz3n),
     .              six(mxsizn),siy(mxsizn),siz(mxsizn),
     .              sjx(mxsizn),sjy(mxsizn),sjz(mxsizn),
     .              skx(mxsizn),sky(mxsizn),skz(mxsizn)
c
      dimension     cx(mbloc),cy(mbloc),cz(mbloc),
     .             cxv(mbloc),cyv(mbloc),czv(mbloc),
     .             cmx(mbloc),cmy(mbloc),cmz(mbloc),ivisb(3,mbloc)
c
      dimension     ri1(1),ri2(1),rj1(1),rj2(1),rk1(1),rk2(1)
      dimension     tauwfi1(1),tauwfi2(1),tauwfj1(1),tauwfj2(1),
     .              tauwfk1(1),tauwfk2(1)
c
c
c     array declaration for temporary work-space
c
c     2-d arrays
c
      dimension   wk2da(1), wk2db(1),
     .            wk2dc(1), wk2dd(1) 
c
c
***********************************************************************
c------------------  mapping function description ---------------------
c
c      imap  : mapping function containing topological information
c      msegt : maximum no. of types of operations/boundary conditons
c              (currently set to 13)
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
c      imap(13, ) : > 0  compute force/moment coefficients
c                   = 0  skip force/moment computation
c
c
c---------------  boundary/operation type description ------------------
c
c      nbctype    = imap(1 , )
c                 = 2  corresponds to a slip (inviscid) wall
c                 = 3  corresponds to a noslip (viscous) wall
c
c***********************************************************************
c
       ref      = scal*scal/sref
c
c      initialize the global(total) force/moment coefficients
c
       cxg      =  0.
       cyg      =  0.
       czg      =  0.
       cxvg     =  0.
       cyvg     =  0.
       czvg     =  0.
       cmxg     =  0.
       cmyg     =  0.
       cmzg     =  0.
c
c      initialize the force/moment coefficients for all blocks
c
       do 10 ibloc = 1,nbloc
       cx (ibloc)   =  0.
       cy (ibloc)   =  0.
       cz (ibloc)   =  0.
       cxv(ibloc)   =  0.
       cyv(ibloc)   =  0.
       czv(ibloc)   =  0.
       cmx(ibloc)   =  0.
       cmy(ibloc)   =  0.
       cmz(ibloc)   =  0.
  10   continue
c
       if (iforce .eq. 0) return
c
c----------  begin outer loop on the blocks  ----------------------------
c
       do 1000 ibloc = 1,nbloc
# if defined BUILD_MPI
       if (nodes(ibloc)-1.eq.myrank) then
# else
       if (nodes(ibloc).eq.myrank) then
# endif
c
       ns       =  nseg(ibloc)
c
c----------  begin outer loop on the segments  ----------------------------
c
       do 100 iseg = 1,ns
c
       nbctype  =  imap(1 ,iseg ,ibloc)
       nface    =  imap(2 ,iseg ,ibloc)
       n1beg    =  imap(3 ,iseg ,ibloc)
       n1end    =  imap(4 ,iseg ,ibloc)
       n2beg    =  imap(5 ,iseg ,ibloc)
       n2end    =  imap(6 ,iseg ,ibloc)
c
c      check to see if forces are to be computed on this segment
c
c      if(imap(13,iseg,ibloc) .eq. 0) go to 100
c
c
c
c TBD
c---   temporary change to exclude contribution of sting to forces
c      for refh-10-block wing/body (By vatsa 8-2-95)
c
c      if(ibloc.eq.10.and.iseg.eq.4) go to 100
c      if(ibloc.eq.06.and.iseg.eq.5) go to 100
c
c
       if(nbctype.eq.2.or.nbctype.eq. 3.or.nbctype.eq.102) then
         if(imap(8 ,iseg ,ibloc).ge.0) then
c
          call forcblk (imn,jmn,kmn,
     .    im  (igrid,ibloc), jm  (igrid,ibloc), km  (igrid,ibloc),
     .    imp1(igrid,ibloc), jmp1(igrid,ibloc), kmp1(igrid,ibloc),
     .    imp2(igrid,ibloc), jmp2(igrid,ibloc), kmp2(igrid,ibloc),
     .       w(m5cc(igrid,ibloc)),   p(m1cc(igrid,ibloc)),
     .     vol(m1cc(igrid,ibloc )),  x(m3np(igrid,ibloc)),
     .     six(m1np(igrid,ibloc)), siy(m1np(igrid,ibloc)),
     .     siz(m1np(igrid,ibloc)), sjx(m1np(igrid,ibloc)),
     .     sjy(m1np(igrid,ibloc)), sjz(m1np(igrid,ibloc)),
     .     skx(m1np(igrid,ibloc)), sky(m1np(igrid,ibloc)),
     .     skz(m1np(igrid,ibloc)),
     .     tauwfi1(m1is(igrid,ibloc)), tauwfi2(m1is(igrid,ibloc)),
     .     tauwfj1(m1js(igrid,ibloc)), tauwfj2(m1js(igrid,ibloc)),
     .     tauwfk1(m1ks(igrid,ibloc)), tauwfk2(m1ks(igrid,ibloc)),
     .     ri1(m3is(igrid,ibloc)), ri2(m3is(igrid,ibloc)),
     .     rj1(m3js(igrid,ibloc)), rj2(m3js(igrid,ibloc)),
     .     rk1(m3ks(igrid,ibloc)), rk2(m3ks(igrid,ibloc)),
     .     mx1dwk ,wk2da, wk2db, wk2dc, wk2dd,
     .     nface ,n1beg ,n1end ,n2beg ,n2end , 
     .     scal,  xref,  yref,  zref,
     .     cx (ibloc), cy (ibloc), cz (ibloc),
     .     cxv(ibloc), cyv(ibloc), czv(ibloc),
     .     cmx(ibloc), cmy(ibloc), cmz(ibloc), nbctype, ivisb(1,ibloc))
c
          endif
        endif
c
c----      end loop on segments
  100  continue
c
c      normalize the force/moment coefficients
c
       cx (ibloc)   =  cx (ibloc)*ref
       cy (ibloc)   =  cy (ibloc)*ref
       cz (ibloc)   =  cz (ibloc)*ref
       cxv(ibloc)   =  cxv(ibloc)*ref
       cyv(ibloc)   =  cyv(ibloc)*ref
       czv(ibloc)   =  czv(ibloc)*ref
       cmx(ibloc)   =  cmx(ibloc)*ref/scal
       cmy(ibloc)   =  cmy(ibloc)*ref/scal
       cmz(ibloc)   =  cmz(ibloc)*ref/scal
c
c      add the force/moment coefficients to obtain total (global)
c      values for the complete configuration 
       cxg      =  cxg  + cx (ibloc)
       cyg      =  cyg  + cy (ibloc)
       czg      =  czg  + cz (ibloc)
       cxvg     =  cxvg + cxv(ibloc)
       cyvg     =  cyvg + cyv(ibloc)
       czvg     =  czvg + czv(ibloc)
       cmxg     =  cmxg + cmx(ibloc)
       cmyg     =  cmyg + cmy(ibloc)
       cmzg     =  cmzg + cmz(ibloc)
c
      endif
 1000 continue
c
c
      return
      end
