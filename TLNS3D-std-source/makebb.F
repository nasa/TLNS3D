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
	subroutine makebb(minpb0,maxlv0,nwv0,wx,wy,wz,
     .                 nbb0,bblen,bbdef,iv,vlist,iwrk)

	dimension wx(*), wy(*), wz(*)
	dimension iwrk(*)
        dimension bbdef(6,*)
	integer bblen,iv(2,*), vlist(*)

	integer totbbv
	common /bbcom/ nbbv,minpbb,maxlev,nwv,nbb,totbbv,
     .                 bbmin(3), bbmax(3)
	integer bbi

	minpbb = minpb0
	maxlev = maxlv0
	nwv = nwv0
	nbb = 0
	totbbv = 0
c
c       The array bbi is a list of the surface points in wx,wy,wz.
c
        bbi = iialloc(nwv)
	do 10 i=1,nwv
	   iwrk(bbi+i-1) = i
 10	continue
c
c       Create a bounding box containing all the surface points.
c
	call calcbb(bbmin,bbmax,nwv,iwrk(bbi),wx,wy,wz)

c
c       Now subdivide the bounding box. Subdivision occurs maxlev
c       times, or until there are less than minpbb points in the
c       bounding box, whichever comes first.
c
	level = 0
	call spltbb(level,bbdef,iv,vlist,iwrk(bbi),wx,wy,wz)
	call ifree(nwv)
        nbb0 = nbb
	bblen = totbbv
c	print *,'nbb0,bblen: ',nbb0,bblen

	return
	end

	subroutine calcbb(bbmin,bbmax,nvi,vi,wx,wy,wz)
	parameter(bval=1.0e20, sval=-1.0e20)
	
	dimension bbmin(3), bbmax(3)
	integer vi(*)
	dimension wx(*), wy(*), wz(*)

c
c       vi is an array of nvi indices into the wx,wy,wz arrays.
c       These indices represent nvi points. This subroutine
c       creates a bounding box that contains these nvi points.
c
	do 10 i=1,3
	   bbmin(i) = bval
	   bbmax(i) = sval
 10	continue

	do 20 j=1,nvi
	   i = vi(j)
	   if(wx(i).lt.bbmin(1)) bbmin(1) = wx(i)
	   if(wy(i).lt.bbmin(2)) bbmin(2) = wy(i)
	   if(wz(i).lt.bbmin(3)) bbmin(3) = wz(i)
	   if(wx(i).gt.bbmax(1)) bbmax(1) = wx(i)
	   if(wy(i).gt.bbmax(2)) bbmax(2) = wy(i)
	   if(wz(i).gt.bbmax(3)) bbmax(3) = wz(i)
 20	continue

	return
	end

	subroutine getvrt(bbmin,bbmax,nwv,nbbv,bbi,wx,wy,wz)
	dimension bbmin(3), bbmax(3)
	integer bbi(*)
	dimension wx(*), wy(*), wz(*)

c
c       This subroutine searches the wx,wy,wz arrays to find
c       all points that fall within a bounding box.
c
	nbbv = 0
	do 10 i=1,nwv
	   if( (wx(i).ge.bbmin(1)).and.(wx(i).le.bbmax(1)).and.
     .         (wy(i).ge.bbmin(2)).and.(wy(i).le.bbmax(2)).and.
     .         (wz(i).ge.bbmin(3)).and.(wz(i).le.bbmax(3)))then
		nbbv = nbbv + 1
		bbi(nbbv) = i
	   endif
 10	continue

	return
	end

	subroutine shells(n,is,s)
	dimension is(*), s(*)

c
c       Shell sort (see Knuth)
c
	integer gap,temp

	gap = n / 2
 100	if(gap.le.0) goto 999
	   do 10 i=gap+1,n
		j = i - gap
 101		if (j.le.0) goto 10
		if (s(is(j)).le.s(is(j+gap))) goto 10
		   temp = is(j)
		   is(j) = is(j+gap)
		   is(j+gap) = temp
		   j = j - gap
		   goto 101
 10	   continue
	   gap = gap / 2
	   goto 100
 999	continue

	return
	end

	subroutine spltbb(level0,bbdef,iv,vlist,bbi,wx,wy,wz)
	parameter(MAXBB=200)
	dimension bbdef(6,*)
	integer bbi(*), iv(2,*), vlist(*)
	dimension wx(*), wy(*), wz(*)
	
	integer totbbv
	common /bbcom/ nbbv,minpbb,maxlev,nwv,nbb,totbbv,
     .                 bbmin(3), bbmax(3)

	dimension stack(7,MAXBB)
	dimension parmin(3), parmax(3)
	dimension plmin(3), plmax(3)
	dimension prmin(3), prmax(3)

	isptr = 1
	call push(stack,isptr,bbmin,bbmax,level0)
	k = 1

 100	if (isptr.le.1) goto 999
	   call pop(stack,isptr,parmin,parmax,level)
	   call getvrt(parmin,parmax,nwv,nbbv,bbi,wx,wy,wz)
	   if (nbbv.ge.minpbb.and.level.lt.maxlev) then
		bx = parmax(1) - parmin(1)
		by = parmax(2) - parmin(2)
		bz = parmax(3) - parmin(3)
		if (bx.gt.by.and.bx.gt.bz) then
		   call shells(nbbv,bbi,wx)
		   iright = nbbv/2
		   i1 = bbi(iright)
		   i2 = bbi(iright+1)
 101		   if(wx(i1).ne.wx(i2).or.iright.ge.nbbv)goto 102
			iright = iright + 1
		        i1 = bbi(iright)
		        i2 = bbi(iright+1)
			goto 101
 102		   continue
		elseif (by.gt.bx.and.by.gt.bz) then
		   call shells(nbbv,bbi,wy)
		   iright = nbbv/2
		   i1 = bbi(iright)
		   i2 = bbi(iright+1)
 103		   if(wy(i1).ne.wy(i2).or.iright.ge.nbbv)goto 104
			iright = iright + 1
		        i1 = bbi(iright)
		        i2 = bbi(iright+1)
			goto 103
 104		   continue
		else
		   call shells(nbbv,bbi,wz)
		   iright = nbbv/2
		   i1 = bbi(iright)
		   i2 = bbi(iright+1)
 105		   if(wz(i1).ne.wz(i2).or.iright.ge.nbbv)goto 106
			iright = iright + 1
		        i1 = bbi(iright)
		        i2 = bbi(iright+1)
			goto 105
 106		   continue
		endif
		call calcbb(plmin,plmax,iright,bbi,wx,wy,wz)
		call calcbb(prmin,prmax,nbbv-iright,bbi(iright+1),
     .			    wx,wy,wz)
		call push(stack,isptr,prmin,prmax,level+1)
		call push(stack,isptr,plmin,plmax,level+1)
	   else
		if (nbbv.gt.0) then
		   nbb = nbb + 1
		   bbdef(1,nbb) = parmin(1)
		   bbdef(2,nbb) = parmax(1)
		   bbdef(3,nbb) = parmin(2)
		   bbdef(4,nbb) = parmax(2)
		   bbdef(5,nbb) = parmin(3)
		   bbdef(6,nbb) = parmax(3)
		   iv(1,nbb) = nbbv
		   iv(2,nbb) = k
		   do 10 i=1,nbbv
			vlist(k) = bbi(i)
			k = k + 1
			totbbv = totbbv + 1
 10		   continue
		endif
	   endif
	   goto 100
 999	continue

	return
	end

	subroutine push(st,is,bbmin,bbmax,lev)
	dimension st(7,*), bbmin(3),bbmax(3)

	do 10 i=1,3
	   ii = i+3
	   st(i,is) = bbmin(i)
	   st(ii,is) = bbmax(i)
 10	continue
	st(7,is) = lev + 0.1
	is = is + 1

	return
	end

	subroutine pop(st,is,bbmin,bbmax,lev)
	dimension st(7,*), bbmin(3),bbmax(3)

	is = is - 1
	do 10 i=1,3
	   ii = i+3
	   bbmin(i) = st(i,is)
	   bbmax(i) = st(ii,is)
 10	continue
	lev = st(7,is)

	return
	end

