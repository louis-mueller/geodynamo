x!real(dp)::rb,rtc,Tb,Tc,km,rhom,rhoc,Hm                     ! dimensionalized variables
      !real(dp)::cm0,cc0,cm1,cc1                                  ! integration constants
      !integer,parameter::n_la=5
      !integer::ipiv(n_la),info
      !double precision, dimension(n_la,n_la)::A                  ! linear system variables
      !double precision, dimension(n_la)::b


if (double_TBL > 0) then
        ! Secular cooling of the core
core%QS = -midpoint_integration(core%r,core%r**2.0_dp*core%rho*core%Cp*core%T)/core%Tc1*4.0_dp*pi
        
core%dT_old = core%Tc1_old - core%Tc1
core%Tc1_old = core%Tc1
! print *, "dT_old: ", core%dT_old, "K"
! latent heat release when core cools 
if (core%dT_old > 0.0_dp) then
  core%QL = core%ri**2.0_dp*core%LH*core%Ti/((core%dTm_dp-core%dT_dp)*core%gi*core%Tc1)*4.0_dp*pi
else
  core%QL = 0.0_dp ! Todo: add latent heat release when inner core partially melts
end if

! radiogenic heat release
core%QR = midpoint_integration(core%r,core%rho*core%r**2.0_dp)*core%Hc*4.0_dp*pi

! pressure heating
core%QP = midpoint_integration(core%r,core%alpha*core%T*core%r**2.0_dp)*core%PT*4.0_dp*pi

! total heat flux at the CMB redimentionalized
core%QCMB = core%kc*(core%Tc-core%Tc1)/core%delta*4.0_dp*pi*core%rc**2.0_dp

! total heat flux at the core-mantle boundary after eq. 78 
! in Nimmo et al. 2015 an update of Gubbins et al. 2003: dTc/dt = (Qcmb-QR)/~QT
core%Tc1 = core%Tc1 - (core%QCMB-core%QR)/(core%QS+core%QP)*dt !core%QL

! ToDo: find the correct definitions of rb, rc, Tb, Ttc, km, kc, rhom, Hm, rhoc, Hc.
rhom = 5536.0_dp ! ToDo: which rho (mantle) should be used here, if rho not constant in mantle? reference value, average value or value directly above CMB?
rhoc = core%rho(core%n) ! Same question as above
Hm = 2.0e-12_dp !radiogenic_heating(params,params%pI%H0,params%pI%lambda,t)*params%pF%k*params%pF%DeltaT/params%pF%D**2.0_dp

core%rc1 = core%rc - core%delta
core%rm1 = core%rc + (mesh%rc(1)-mesh%rmin)*params%pF%D !thickness of first cell above the cmb + core radius

! radial heat conduction in a sphere or sphereical shell 
! dq/dr + 2q/r = pho*H --> T = -rho*H/(6*k) + c1/r + c0

!solve linear system of integration constants and core-mantle-boundary temperature
!  1   2   3   4   5
!(cm1,cm0,cc1,cc0,Tc) 

A = 0.0_dp

! row cm1
A(1,1) = 1.0_dp
A(1,2) = 1.0_dp/core%rm1
b(1) = core%Tm1 + rhom*Hm/(6.0_dp*core%km)*core%rm1**2

! row cm0
A(2,2) = 1.0_dp
A(2,5) = -core%rc*core%rm1/(core%rm1-core%rc)
b(2) = core%Tm1*core%rm1*core%rc/(core%rm1-core%rc) - rhom*Hm/(6.0_dp*core%km)*(core%rc**2*core%rm1 + core%rc*core%rm1**2)

! row cc1
A(3,3) = 1.0_dp
A(3,4) = 1.0_dp/core%rc
A(3,5) = -1.0_dp
b(3) = rhoc*core%Hc/(6.0_dp*core%kc)*core%rc**2

! row cc0
A(4,4) = 1.0_dp
A(4,5) = core%rc1*core%rc/(core%rc-core%rc1)
b(4) = core%Tc1*core%rc1*core%rc/(core%rc-core%rc1) - rhoc*core%Hc/(6.0_dp*core%kc)*(core%rc**2*core%rc1 + core%rc*core%rc1**2) 

! row Tc
A(5,5) = 1.0_dp
b(5) = (core%Tc1 + core%Tm1*core%km/core%kc*(core%rc-core%rc1)/(core%rm1-core%rc))/ &
        (1 + core%km/core%kc*(core%rc-core%rc1)/(core%rm1-core%rc))

! Call gsev from LAPACK to solve the system A * x = b
call dgesv(n_la, 1, A, n_la, ipiv, b, n_la, info)

! Check if the solution was successful
if (info /= 0) then
  print *, "Error: The solution could not be computed."
  stop
end if

!B is overwritten by the result in dgesv() 
cm1     = b(1)
cm0     = b(2)
cc1     = b(3)
cc0     = b(4)
core%Tc = b(5) !new temperature at the core-mantle-boundary

! pass cmb temp. to CHIC 
params%pT%Tbottom = (core%Tc-params%pF%Ts)/params%pF%DeltaT

core%dTc1 = abs(core%T(core%n-1) - core%Tc1) / dt !change to last time step dTc/dt = (T(n)-Tc)/dt 

core%T(core%n-1) = core%Tc1

do i = core%n-1, 2, -1
  core%T(i-1) = core%T(i) + core%alpha(i) * core%g(i) / core%Cp(i) * core%T(i) * core%dr
end do

! heat output units [TW]
core%QR = core%QR*1.0e-12_dp
core%QS = core%QS*core%dTc1*1.0e-12_dp
core%QL = core%QL*core%dTc1*1.0e-12_dp
core%QP = core%QP*core%dTc1*1.0e-12_dp
else