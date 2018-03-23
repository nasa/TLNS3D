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
      SUBROUTINE HEAP_SORT(N,ARRIN,INDX)
      REAL ARRIN(N),Q
      INTEGER INDX(N)
      
      DO 11 J=1,N
         INDX(J) = J
 11   CONTINUE

      if(n .eq. 1)return
      
      L = N/2 + 1
      IR= N
 10   CONTINUE
      IF(L.GT.1)THEN
         L=L-1
         INDXT = INDX(L)
         Q     = ARRIN(INDXT)
      ELSE
         INDXT = INDX(IR)
         Q     = ARRIN(INDXT)
         INDX(IR) = INDX(1)
         IR = IR - 1
         IF(IR.EQ.1)THEN
            INDX(1) = INDXT
            RETURN
         ENDIF
      ENDIF
      I=L
      J=L+L
 20   IF(J.LE.IR)THEN
         IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
         ENDIF
         IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I) = INDX(J)
            I = J
            J = J+J
         ELSE
            J = IR+1
         ENDIF
         GO TO 20
      ENDIF
      INDX(I) = INDXT
      GO TO 10
      END
