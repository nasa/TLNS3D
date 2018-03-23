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
c======================================================================
      subroutine savplt(imn,jmn,kmn,im,jm,km,
     .                  mgrlev,mbloc,nbloc,mres,nres,
     .                  wcount,rtrmsg,hrmsg,nsupg,
     .                  clg,cmxg,cmyg,cmzg,cdtg,cdvg,czg,
     .                  cl ,cmx ,cmy ,cmz ,cdt ,cdv, cz )
c=======================================================================
c
c     save convergence history as ViGraph plot file
c 
c
      common/files/ iread,iwrit,igrdfil,irstfil,isavfil,ipltfil
     .                   ,imapfil,ipfil1,ipfil2,irhtall,irhtphy
     .                   ,igpfil,iqpfil,idatfil,ipntfil,iavgfil
c
      dimension im(mgrlev,mbloc), jm(mgrlev,mbloc),
     .          km(mgrlev,mbloc)
c
      dimension wcount(mres),rtrmsg(mres),hrmsg(mres),nsupg(mres),
     .             clg(mres),  cmxg(mres), cmyg(mres), cmzg(mres),
     .            cdtg(mres),  cdvg(mres),  czg(mres)
c
      dimension    cl(mbloc),  cmx(mbloc), cmy(mbloc), cmz(mbloc),
     .            cdt(mbloc),  cdv(mbloc),  cz(mbloc)
c
      write (ipltfil,'("12")')
      write (ipltfil,'("iter")')
      write (ipltfil,'("work")')
      write (ipltfil,'("rrms")')
      write (ipltfil,'("hrms")')
      write (ipltfil,'("nsup")')
      write (ipltfil,'("clg")')
      write (ipltfil,'("cmxg")')
      write (ipltfil,'("cmyg")')
      write (ipltfil,'("cmzg")')
      write (ipltfil,'("cdtg")')
      write (ipltfil,'("cdvg")')
      write (ipltfil,'("czg")')
      write (ipltfil,'(2i5)') nres,1
c
      do ires=1,nres
        sup = nsupg(ires)
        write (ipltfil,'(i5,11e12.4)') ires,wcount(ires),
     .    rtrmsg(ires),hrmsg(ires),sup,clg(ires),cmxg(ires),
     .    cmyg(ires),cmzg(ires),cdtg(ires),cdvg(ires),czg(ires)
      enddo
c
c     write (ipltfil,'("9")')
c     write (ipltfil,'("block")')
c     write (ipltfil,'("points")')
c     write (ipltfil,'("cl")')
c     write (ipltfil,'("cmx")')
c     write (ipltfil,'("cmy")')
c     write (ipltfil,'("cmz")')
c     write (ipltfil,'("cdt")')
c     write (ipltfil,'("cdv")')
c     write (ipltfil,'("cz")')
c     write (ipltfil,'(2i5)') nbloc,1
c
c     do ibloc=1,nbloc
c       npts = (im(1,ibloc) - imn + 1) *
c    .         (jm(1,ibloc) - jmn + 1) *
c    .         (km(1,ibloc) - kmn + 1)
c       write (ipltfil,'(i5,i8,7e12.4)') ibloc,npts,
c    .    cl(ibloc),cmx(ibloc),cmy(ibloc),cmz(ibloc),
c    .    cdt(ibloc),cdv(ibloc),cz(ibloc)
c     enddo
c
      return
      end 
