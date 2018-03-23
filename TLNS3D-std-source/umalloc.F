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
      subroutine umalloc(ipoint,n_words,intflag,icall,memuse)
c***********************************************************************
c     Purpose:  universal malloc...allocates memory via fortran 77
c     pointer extension.
c
c     NOTE: in this special-purpose version, allocates memory
c     for REAL (not complex real) arrays, as well as integer arrays.
c
c     icall   - integer counter to aid in debugging
c     intflag - flag to set variable type:
c               0...real variable
c               2...integer variable
c               3...character variable
c     memuse  - running total of memory that has been allocated
c
c     original version provied courtesy Steven J. Massey, AS&M Inc.
c     modified by r. t. biedron, January, 1999.
c***********************************************************************
c
      parameter(len_i=4,len_c=1)
c
      integer errcode,abort
c
      dimension dum(1)
c
      pointer (ipoint, dum)
c
#ifdef CRAY_TIME
c     cray_time implies cray
      len_r = 8
#else
#  ifdef TLN_DBL
      len_r = 8
#  else 
      len_r = 4
#  endif
#endif
c
      icall = icall+1
c
      if (intflag.eq.1) then
        len = len_i
      else if (intflag.eq.0) then
        len = len_r
      else
        len = len_c
      end if
c
      n_bytes = n_words * len
c
#ifdef CRAY_TIME
c     cray_time implies cray
      abort = 0
      errcode = 0
      nwords = n_words
      call hpalloc( ipoint, nwords, errcode, abort)
#else
      ipoint = malloc(n_bytes)
#endif
c
      if (ipoint.eq.0) then
        write(6,'(a,i4)') 
     .  'stopping...malloc failed to allocate the memory on call ',
     .  icall
        if (n_bytes.lt.1e10) then 
          write(6,12) n_words, n_bytes
        else if (n_bytes.lt.1e15) then 
          write(6,13) n_words, n_bytes
        else 
          write(6,*) 'requested more than 10^(15) bytes'
        end if
        write(6,14) memuse
        stop
      else
c        write(6,*) 'call',icall,' umalloc: allocated ',
c    .   n_bytes,' bytes'
      end if
c
  12  format('requested ', i10, ' words, (',i10, ' bytes)')
  13  format('requested ', i15, ' words, (',i15, ' bytes)')
  14  format('prior to this call, total allocation was',
     .        i15, ' bytes')
c
      memuse  = memuse + n_bytes
c
      return
      end
