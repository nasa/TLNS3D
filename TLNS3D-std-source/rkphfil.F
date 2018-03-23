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
      subroutine rkphfil(irkcase,mrkmax,rkphysa,rkphyse,
     .                   sorder)
c
c     object : according to runge kutta case fill
c              parameters for implicit runge kutta
c              multi stage method for physical time
c
      dimension rkphysa(mrkmax,mrkmax),rkphyse(mrkmax)
c
c
      if (irkcase.eq.1) then
  
       rkphysa(2,1) = 1767732205903.d0/4055673282236.d0
       rkphysa(2,2) = 1767732205903.d0/4055673282236.d0
       rkphysa(3,1) = 2746238789719.d0/10658868560708.d0
       rkphysa(3,2) =-640167445237.d0/6845629431997.d0
       rkphysa(3,3) = 1767732205903.d0/4055673282236.d0
       rkphysa(4,1) = 1471266399579.d0/7840856788654.d0
       rkphysa(4,2) =-4482444167858.d0/7529755066697.d0
       rkphysa(4,3) = 11266239266428.d0/11593286722821.d0
       rkphysa(4,4) = 1767732205903.d0/4055673282236.d0
       c1 = 0.
       c2 = rkphysa(2,1) + rkphysa(2,2)
       c3 = rkphysa(3,1) + rkphysa(3,2) + rkphysa(3,3)
       c4 = rkphysa(4,1) + rkphysa(4,2) + rkphysa(4,3) + rkphysa(4,4)
c      write(*,*)c1,c2,c3,c4

       rkphyse(1)   = 2756255671327.d0/12835298489170.d0 
       rkphyse(2)   =-10771552573575.d0/22201958757719.d0
       rkphyse(3)   = 9247589265047.d0/10645013368117.d0
       rkphyse(4)   = 2193209047091.d0/5459859503100.d0

       sorder = 2.0

       elseif (irkcase.eq.2) then

        rkphysa(2,1) = 2.d0/3.d0
        rkphysa(3,1) = 1.d0/3.d0
        rkphysa(3,2) = 1.d0/3.d0
        rkphysa(4,1) = 1.d0/4.d0
        rkphysa(4,2) = 0.d0
        rkphysa(4,3) = 3.d0/4.d0

        rkphyse(1) = 1.d0/4.d0
        rkphyse(2) = 0.d0
        rkphyse(3) = 3.d0/4.d0

        sorder = 2.0 

       elseif (irkcase.eq.3) then

       gamma = 0.435866521508458999
       b2    = -(3.d0/2.d0)*gamma**2 + 4.d0*gamma - 1.d0/4.d0
       b3    =  (3.d0/2.d0)*gamma**2 - 5.d0*gamma + 5.d0/4.d0
c
       rkphyse(1) = 1./6.
       rkphyse(2) = 2./6.
       rkphyse(3) = 2./6.
       rkphyse(4) = 1./6.
c
       sorder = 2.0
c
       rkphysa(2,1) = 0.d0
       rkphysa(2,2) = gamma
       rkphysa(3,1) = 0.d0
       rkphysa(3,2) = (1.d0 - gamma)/2.d0
       rkphysa(3,3) = gamma
       rkphysa(4,1) = 0.
       rkphysa(4,2) = b2
       rkphysa(4,3) = b3
       rkphysa(4,4) = gamma

       elseif (irkcase.eq.4) then
c      print*,
c   & 'M.P. Calvo et al., LIRK4, Rep. 1999/3, Valladolid, Spain (1999)'
c       write(*,*)'getting padded SDIRK'

        sorder = 3.0

        rkphyse(1) = 0.d0
        rkphyse(2) = 59.d0/48.d0
        rkphyse(3) = -17.d0/96.d0
        rkphyse(4) = 225.d0/32.d0
        rkphyse(5) =-85.d0/12.d0
        rkphyse(6) = 0.d0
c
        rkphysa(1,1) = 0.d0
        rkphysa(2,1) = 0.d0
        rkphysa(2,2) = 1.d0/4.d0
        rkphysa(3,1) = 0.d0
        rkphysa(3,2) = 1.d0/2.d0
        rkphysa(3,3) = 1.d0/4.d0
        rkphysa(4,1) = 0.d0
        rkphysa(4,2) = 17.d0/50.d0
        rkphysa(4,3) =-1.d0/25.d0
        rkphysa(4,4) = 1.d0/4.d0
        rkphysa(5,1) = 0.d0
        rkphysa(5,2) = 371.d0/1360.d0
        rkphysa(5,3) =-137.d0/2720.d0
        rkphysa(5,4) = 15.d0/544.d0
        rkphysa(5,5) = 1.d0/4.d0
        rkphysa(6,1) = 0.d0
        rkphysa(6,2) = 25.d0/24.d0
        rkphysa(6,3) = -49.d0/48.d0
        rkphysa(6,4) = 125.d0/16.d0
        rkphysa(6,5) =-85.d0/12.d0
        rkphysa(6,6) = 1.d0/4.d0
c

       elseif (irkcase.eq.5) then
c       write(*,*)'getting ESDIRK'

c      print*,'ARK4(3)6L[2]SAc'

       sorder = 3.0

       rkphyse(1) = +4586570599.d0/29645900160.d0
       rkphyse(2) = 0.d0
       rkphyse(3) = +178811875.d0/945068544.d0
       rkphyse(4) = +814220225.d0/1159782912.d0
       rkphyse(5) = -3700637.d0/11593932.d0
       rkphyse(6) = +61727.d0/225920.d0

       rkphysa(2,1) = 1.d0/4.d0
       rkphysa(2,2) = 1.d0/4.d0
       rkphysa(3,1) = 8611.d0/62500.d0
       rkphysa(3,2) =-1743.d0/31250.d0
       rkphysa(3,3) = 1.d0/4.d0
       rkphysa(4,1) = 5012029.d0/34652500.d0
       rkphysa(4,2) =-654441.d0/2922500.d0
       rkphysa(4,3) = 174375.d0/388108.d0
       rkphysa(4,4) = 1.d0/4.d0
       rkphysa(5,1) = 15267082809.d0/155376265600.d0
       rkphysa(5,2) =-71443401.d0/120774400.d0
       rkphysa(5,3) = 730878875.d0/902184768.d0
       rkphysa(5,4) = 2285395.d0/8070912.d0
       rkphysa(5,5) = 1.d0/4.d0
       rkphysa(6,1) = 82889.d0/524892.d0
       rkphysa(6,2) = 0.d0
       rkphysa(6,3) = 15625.d0/83664.d0
       rkphysa(6,4) = 69875.d0/102672.d0
       rkphysa(6,5) =-2260.d0/8211.d0
       rkphysa(6,6) = 1.d0/4.d0

       elseif (irkcase.eq.6) then

       rkphyse(1) = 5974964756349.d0/12345176340353.d0
       rkphyse(2) = 1445407564659.d0/ 8949752891087.d0
       rkphyse(3) =-2098652516771.d0/ 4895426490835.d0
       rkphyse(4) = 6245407545074.d0/ 9990007552309.d0
       rkphyse(5) =   30734188585.d0/10315892450192.d0
       rkphyse(6) = 1128743561416.d0/ 7279518462557.d0
c
       rkphysa(2,1) =  3282482714977.d0/11805205429139.d0
       rkphysa(2,2) =  3282482714977.d0/11805205429139.d0
       rkphysa(3,1) = -1748659331805.d0/ 9772184669482.d0
       rkphysa(3,2) = - 275228881984.d0/ 9116607863767.d0
       rkphysa(3,3) =  3282482714977.d0/11805205429139.d0
       rkphysa(4,1) = - 613979224052.d0/  905787279025.d0
       rkphysa(4,2) = -3072163075577.d0/18154595953427.d0
       rkphysa(4,3) =  8639835010822.d0/ 8478664872845.d0
       rkphysa(4,4) =  3282482714977.d0/11805205429139.d0
       rkphysa(5,1) = 23868300505303.d0/13380250758717.d0
       rkphysa(5,2) = -3225335589184.d0/ 5272415062559.d0
       rkphysa(5,3) =-42661694614627.d0/20157616768973.d0
       rkphysa(5,4) = 24229026669092.d0/13526367824531.d0
       rkphysa(5,5) =  3282482714977.d0/11805205429139.d0
       rkphysa(6,1) = -1288628682006.d0/17478471813649.d0
       rkphysa(6,2) =  3215146213920.d0/ 9235068495727.d0
       rkphysa(6,3) =  2831045640924.d0/ 8795426369095.d0
       rkphysa(6,4) =   564478078147.d0/ 2818192429028.d0
       rkphysa(6,5) = - 552280568451.d0/ 7398512356043.d0
       rkphysa(6,6) =  3282482714977.d0/11805205429139.d0

       elseif (irkcase.eq.7) then

       rkphyse(1) = 1./5.
       rkphyse(2) = 1./5.
       rkphyse(3) = 1./5.
       rkphyse(4) = 1./5.
       rkphyse(5) = 1./5.

       rkphysa(2,1) =  11./40.
       rkphysa(2,2) =  11./40.
       rkphysa(3,1) =  8751008911./25635583160.
       rkphysa(3,2) =  100865621./640889579.
       rkphysa(3,3) =  11./40.
       rkphysa(4,1) =  78553077359./161791326400.
       rkphysa(4,2) = -8302193139./46785526040.
       rkphysa(4,3) =  1760428191262./4222346697671.
       rkphysa(4,4) =  11./40.
       rkphysa(5,1) =  3658861./15591840.
       rkphysa(5,2) =  1297565./3381543.
       rkphysa(5,3) =  3113033735959./7601233198560.
       rkphysa(5,4) = -34243./113040.
       rkphysa(5,5) =  11./40.

       elseif (irkcase.eq.8) then

       sorder = 4.0

       rkphyse(1) = 354702883506.d0/3789842051159.d0
       rkphyse(2) = 0.d0
       rkphyse(3) =-471861028701.d0/9261158205962.d0
       rkphyse(4) = 3832214029186.d0/7971181559071.d0
       rkphyse(5) =-19522742290.d0/5106068051629.d0
       rkphyse(6) = 16391434965105.d0/17160319691969.d0
       rkphyse(7) =-14302054201217.d0/23091863285126.d0
       rkphyse(8) = 5188720167746.d0/35887481870821.d0

       rkphysa(2,1) = 29.d0/200.d0
       rkphysa(2,2) = 29.d0/200.d0
       rkphysa(3,1) = 29.d0/400.d0
       rkphysa(3,2) =-350623571328.d0/10420535143571.d0
       rkphysa(3,3) = 29.d0/200.d0
       rkphysa(4,1) = 807226490794.d0/11294700695621.d0
       rkphysa(4,2) =-874668641845.d0/7098011249891.d0
       rkphysa(4,3) = 2182184544357.d0/11684570159401.d0
       rkphysa(4,4) = 29.d0/200.d0
       rkphysa(5,1) = 145253674553.d0/8663875643876.d0
       rkphysa(5,2) = 15128786102297.d0/11906071161648.d0
       rkphysa(5,3) = 6152117303669.d0/10771869990990.d0
       rkphysa(5,4) =-15616452376054.d0/9444072785965.d0
       rkphysa(5,5) = 29.d0/200.d0
       rkphysa(6,1) = 6091186570339.d0/15702166824390.d0
       rkphysa(6,2) = 3455924488897.d0/7171035336684.d0
       rkphysa(6,3) =-10979430586178.d0/7250509791871.d0
       rkphysa(6,4) = 22232248077499.d0/13875476324503.d0
       rkphysa(6,5) =-3121948903746.d0/10309592483189.d0
       rkphysa(6,6) = 29.d0/200.d0
       rkphysa(7,1) = 5104696209616.d0/9675258774483.d0
       rkphysa(7,2) = 10291357491335.d0/15459444993913.d0
       rkphysa(7,3) =-28707109686536.d0/13058909764101.d0
       rkphysa(7,4) = 34789977277050.d0/15818991491101.d0
       rkphysa(7,5) =-26433723305593.d0/53141752197286.d0
       rkphysa(7,6) = 1359046691.d0/74921600000.d0
       rkphysa(7,7) = 29.d0/200.d0
       rkphysa(8,1) = 69549.d0/842800.d0
       rkphysa(8,2) = 0.d0
       rkphysa(8,3) = 0.d0
       rkphysa(8,4) = 63475.d0/147784.d0
       rkphysa(8,5) = 1408.d0/337365.d0
       rkphysa(8,6) = 81553.d0/84240.d0
       rkphysa(8,7) =-360200.d0/572373.d0
       rkphysa(8,8) = 29.d0/200.d0

       elseif (irkcase.eq.9) then

       sorder = 4.0

       rkphyse(1) = 674948148687.d0/7791690761692.d0
       rkphyse(2) = 0.d0
       rkphyse(3) = 1977890653369.d0/10850705122729.d0
       rkphyse(4) = 3307959803092.d0/13406079769963.d0
       rkphyse(5) = 10225964142456.d0/17021390495047.d0
       rkphyse(6) =-26981889059942.d0/90301476281381.d0
       rkphyse(7) = 612291744556.d0/3357424677617.d0

       rkphysa(2,1) = 23.d0/125.d0
       rkphysa(2,2) = 23.d0/125.d0
       rkphysa(3,1) = 23.d0/250.d0
       rkphysa(3,2) =-404036019743.d0/9462789448531.d0
       rkphysa(3,3) = 23.d0/125.d0
       rkphysa(4,1) = 1039574512792.d0/10823206779923.d0
       rkphysa(4,2) =-915977319726.d0/6883742396117.d0
       rkphysa(4,3) = 494332862300.d0/2701075432643.d0
       rkphysa(4,4) = 23.d0/125.d0
       rkphysa(5,1) = 1350755389980.d0/12393594817193.d0
       rkphysa(5,2) =-12603447208977.d0/8749348280120.d0
       rkphysa(5,3) =-6628209800229.d0/12303071059028.d0
       rkphysa(5,4) = 2782135010663.d0/1132672415741.d0
       rkphysa(5,5) = 23.d0/125.d0
       rkphysa(6,1) = 2456820185761.d0/13309240171219.d0
       rkphysa(6,2) =-44178692907383.d0/14974654786062.d0
       rkphysa(6,3) =-24266808557751.d0/15857662065622.d0
       rkphysa(6,4) = 70602821193752.d0/14269244229287.d0
       rkphysa(6,5) = 625101859940.d0/9763661771393.d0
       rkphysa(6,6) = 23.d0/125.d0
       rkphysa(7,1) = 1332710774222.d0/19650076401203.d0
       rkphysa(7,2) = 0.d0
       rkphysa(7,3) = 2767517002471.d0/9006457202997.d0
       rkphysa(7,4) = 4311620031746.d0/35864909658309.d0
       rkphysa(7,5) = 4352693668673.d0/6738567275899.d0
       rkphysa(7,6) =-3441835818133.d0/10581815184094.d0
       rkphysa(7,7) = 23.d0/125.d0

       elseif (irkcase.eq.10) then

       rkphyse(1) = 1332710774222.d0/19650076401203.d0
       rkphyse(2) = 0.d0
       rkphyse(3) = 2767517002471.d0/9006457202997.d0
       rkphyse(4) = 4311620031746.d0/35864909658309.d0
       rkphyse(5) = 4352693668673.d0/6738567275899.d0
       rkphyse(6) =-3441835818133.d0/10581815184094.d0
       rkphyse(7) = 23.d0/125.d0

       rkphysa(2,1) =  29.d0/200.d0
       rkphysa(2,2) =  29.d0/200.d0
       rkphysa(3,1) =  29.d0/400.d0
       rkphysa(3,2) = -350623571328.d0/10420535143571.d0
       rkphysa(3,3) =  29.d0/200.d0
       rkphysa(4,1) =  322855878611.d0/13195927472715.d0
       rkphysa(4,2) = -1017779885597.d0/5299915360689.d0
       rkphysa(4,3) =  11833025699139.d0/28002476865679.d0
       rkphysa(4,4) =  29.d0/200.d0
       rkphysa(5,1) =  903925235668.d0/15674108167895.d0
       rkphysa(5,2) =  1963760237863.d0/4344126218749.d0
       rkphysa(5,3) =  590306497559.d0/7872932197458.d0
       rkphysa(5,4) = -959221702753.d0/7395767254095.d0
       rkphysa(5,5) =  29.d0/200.d0
       rkphysa(6,1) = -1703616817423.d0/5678278513847.d0
       rkphysa(6,2) =  36860600553986.d0/9468284938995.d0
       rkphysa(6,3) = -4994439221027.d0/15135774984203.d0
       rkphysa(6,4) = -44361687063259.d0/12665432616375.d0
       rkphysa(6,5) =  29726446095475.d0/33231764886716.d0
       rkphysa(6,6) =  29.d0/200.d0
       rkphysa(7,1) = -83353732129.d0/3036615448297.d0
       rkphysa(7,2) =  0.d0
       rkphysa(7,3) =  7897536467593.d0/12786817195313.d0
       rkphysa(7,4) = -4113767841355.d0/7573121155078.d0
       rkphysa(7,5) =  8274481613094.d0/8817160230589.d0
       rkphysa(7,6) = -950106813193.d0/7284576383586.d0
       rkphysa(7,7) =  29.d0/200.d0
c
      elseif(irkcase.eq.11) then

c      print*,'ARK4(3)6L[2]SA Imp -> C^{[E]}(2)+C^{[I]}(3) '

       rkphysa(2,1) =  1/4.
       rkphysa(2,2) =  1/4.
       rkphysa(3,1) =  1/16.
       rkphysa(3,2) = -1/16.
       rkphysa(3,3) =  1/4.
       rkphysa(4,1) =  -53/100.
       rkphysa(4,2) =  -7/10.
       rkphysa(4,3) =  42/25.
       rkphysa(4,4) =  1/4.
       rkphysa(5,1) = -5/24.
       rkphysa(5,2) = -13/48.
       rkphysa(5,3) =  19/18.
       rkphysa(5,4) =  25/144.
       rkphysa(5,5) =  1/4.
       rkphysa(6,1) =  20/147.
       rkphysa(6,2) =  19/42.
       rkphysa(6,3) =  8/63.
       rkphysa(6,4) =  125/882.
       rkphysa(6,5) = -3/28.
       rkphysa(6,6) =  1/4.

       rkphyse(1) = rkphysa(5,1)
       rkphyse(2) = rkphysa(5,2)
       rkphyse(3) = rkphysa(5,3)
       rkphyse(4) = rkphysa(5,4)
       rkphyse(5) = rkphysa(5,5)
       rkphyse(6) = 0

      elseif(irkcase.eq.12) then

c      print*,'ARK4(3)6L[2]SAe -> C^{[E]}(2)+C^{[I]}(3) '

       rkphysa(2,1) = 1.d0/4.d0
       rkphysa(2,2) = 1.d0/4.d0
       rkphysa(3,1) = 1.d0/8.d0
       rkphysa(3,2) =-685497852586.d0/11816340736199.d0
       rkphysa(3,3) = 1.d0/4.d0
       rkphysa(4,1) = 52920034407.d0/6115188259195.d0
       rkphysa(4,2) =-2581862900803.d0/8162419824178.d0
       rkphysa(4,3) = 15401743608433.d0/18581139755562.d0
       rkphysa(4,4) = 1.d0/4.d0
       rkphysa(5,1) =-898051810951.d0/6555027953449.d0
       rkphysa(5,2) =-2300599048126.d0/3548526458713.d0
       rkphysa(5,3) = 11299807919202.d0/8447519174845.d0
       rkphysa(5,4) = 500561058752.d0/3389519422739.d0
       rkphysa(5,5) = 1.d0/4.d0
       rkphysa(6,1) = 2834623777653.d0/26906817174524.d0
       rkphysa(6,2) = 0.d0
       rkphysa(6,3) = 4079485946257.d0/9358623902187.d0
       rkphysa(6,4) = 1949886071195.d0/4030519482993.d0
       rkphysa(6,5) =-419653551822.d0/1525809906703.d0
       rkphysa(6,6) = 1.d0/4.d0

       rkphyse(1) = +850233530267.d0/7617317741865.d0
       rkphyse(2) = 0.d0
       rkphyse(3) = +2353973301757.d0/5618343433312.d0
       rkphyse(4) = +4600875836615.d0/8897330496886.d0
       rkphyse(5) = -2420296691065.d0/7524082998448.d0
       rkphyse(6) = +3133574984401.d0/11437765310009.d0



       end if
      
c
c
c
      return
      end
