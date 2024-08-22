module core
  use precision
  use mesher, only:mesh_cp
  use variables
  use parameters_total, only:params_tot 
  use heating
  
  implicit none

  type core_properties
    ! note: the subscript old denotes the previous time step
    integer::n                                  ! core resolution
    real(dp)::dr  !m                            ! core radius step size (n dependant)
    real(dp)::ri  !m                            ! radius of the inner core
    real(dp)::rc  !m                            ! core radius
    real(dp),allocatable::g(:)  !m/s^2
    real(dp),allocatable::p(:)  !Pa
    real(dp),allocatable::rho(:)  !kg/m^3
    real(dp),allocatable::r(:)  !m
    real(dp),allocatable::T(:)  !K              ! core temperature profile (adiabatic)
    real(dp),allocatable::Tm(:) !K              ! pure iron melting temperature profile  
    real(dp),allocatable::Cp(:) !J/(kg*K)       ! specific heat capacity
    real(dp),allocatable::alpha(:)  !K^-1       ! thermal expansion coefficient
    real(dp),allocatable::mat(:)                ! material phase coeff. (7.0: lower mantle, 8.0: core)
    real(dp)::Hc  !W/kg                         ! volumetric heating rate (constant because equally distributed?)
    real(dp)::LH  !J/kg                         ! latent heat of iron crystallization
    real(dp)::PT                                ! numerical coeff. to relate pressure change at ICB to core cooling (=1.0_dp?)
    real(dp)::kc  !W/(m*K)                      ! thermal conductivity of the core
    real(dp)::rhoi  !kg/m^3                     ! density of the inner core
    real(dp)::Ttc  !K                           ! temperature at the bottom of the core TBL
    real(dp)::Tc
    real(dp),allocatable::T_old(:) !K           
    real(dp)::Ttc_old  !K                       
    real(dp)::ri_old !m    
    real(dp)::QR,QL,QS,QP,QCMB  ! W             ! heat fluxes                     
  end type core_properties

  contains
    function linspace(start,end,num,endpoint,step) result(samples)
      ! Author: Ivan Pribec
      ! https://gist.github.com/ivan-pi/f4b4741d7ed54ceff787c85d6ba22a5a (20.08.2024)
      real(dp), intent(in) :: start
      real(dp), intent(in) :: end
      integer, intent(in), optional :: num !default 50
      logical, intent(in), optional :: endpoint ! default is true (include endpoint)
      real(dp), intent(out), optional :: step
      real(dp), allocatable :: samples(:)
      integer :: num_, i
      logical :: endpoint_ 
      real(dp) :: step_

      num_ = 50
      if (present(num)) num_ = num

      endpoint_ = .true.
      if (present(endpoint)) endpoint_ = endpoint

      ! find step size
      if (endpoint_) then
          step_ = (end - start)/real(num_-1,dp)
      else
          step_ = (end - start)/real(num_,dp)
      end if

      if (present(step)) step = step_

      allocate(samples(num_))
      do i = 1, num_
          samples(i) = start + (i-1)*step_
      end do
    end function linspace

    function reverse_array(array) result(reversed)
      real(dp), intent(in) :: array(:)
      real(dp), allocatable :: reversed(:)
      integer :: i, n
  
      n = size(array)
      allocate(reversed(n))
      
      do i = 1, n
          reversed(i) = array(n - i + 1)
      end do
    end function reverse_array

    function midpoint_integration(x, y) result(integrand)
      real(dp), intent(in)::x(:),y(:) 
      integer::i,n
      real(dp)::integrand

      n = size(x)
      integrand = 0.0_dp

      if (size(x) /= size(y)) then 
        print *, "Error: Arrays x and y must have the same length."
        stop
      end if
      
      do i = 1, n-1
        integrand = integrand + (y(i) + y(i+1)) / 2.0_dp * (x(i+1) - x(i))
      end do
    end function midpoint_integration

    function cubic_spline_interpolation(x,y,x_interp) result(y_interp)
      real(dp),allocatable,intent(in)::x(:),y(:),x_interp(:)
      real(dp),allocatable::y_interp(:)
      real(dp),allocatable::y2(:),u(:)
      real(dp)::p,sig,qn,un,h,a,b
      integer::i,j,k,klo,khi,n,m

      if (size(x) /= size(y)) then ! Check if the lengths of x and y are equal
        print *, "Error: Arrays x and y must have the same length."
        stop
      end if

      n = size(x); m = size(x_interp)

      allocate(y2(n),u(n-1),y_interp(m))

      ! Step 1: Compute the second derivatives (y2)
      y2(1) = 0.0_8; u(1) = 0.0_8; qn = 0.0_dp; un = 0.0_dp
      do i = 2, n-1
        sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
        p = sig * y2(i-1) + 2.0_dp
        y2(i) = (sig - 1.0_dp) / p
        u(i) = (y(i+1) - y(i)) / (x(i+1) - x(i)) - (y(i) - y(i-1)) / (x(i) - x(i-1))
        u(i) = (6.0_dp * u(i) / (x(i+1) - x(i-1)) - sig * u(i-1)) / p
      end do
      
      y2(n) = (un - qn * u(n-1)) / (qn * y2(n-1) + 1.0_dp) ! Boundary condition

      do k = n-1, 1, -1 ! Back-substitution loop
        y2(k) = y2(k) * y2(k+1) + u(k)
      end do

      ! Step 2: Interpolation using the second derivatives, with binary search
      do j = 1, m
        klo = 1
        khi = n
        do while (khi - klo > 1)
            k = (khi + klo) / 2
            if (x(k) > x_interp(j)) then
                khi = k
            else
                klo = k
            end if
        end do

        h = x(khi) - x(klo)
        if (h == 0.0_dp) stop 'Bad input in cubic_spline_interpolation'

        a = (x(khi) - x_interp(j)) / h
        b = (x_interp(j) - x(klo)) / h
        y_interp(j) = a * y(klo) + b * y(khi) + ((a**3 - a) * y2(klo) + (b**3 - b) * y2(khi)) * (h**2) / 6.0_dp
      end do  
    end function cubic_spline_interpolation

    function QR_func(r,rho,H,dt) result(Q_R)
      real(dp),allocatable,intent(in)::r(:),rho(:)
      real(dp),intent(in)::H,dt
      real(dp),allocatable::qr(:)
      real(dp)::Q_R, pi_

      pi_ = acos(-1.0_dp)

      qr = rho*H*r**2

      Q_R = midpoint_integration(r,qr)*4*pi_/dt
    end function QR_func

    function QL_func(ri,ri_before,L_H,rhoi,dt) result(Q_L)
      real(dp),intent(in)::L_H,ri_before,ri,rhoi,dt
      real(dp)::Q_L, pi_, dri

      pi_ = acos(-1.0_dp)

      dri = ri - ri_before

      Q_L = 4*pi_*ri**2.0_dp*L_H*rhoi*dri/dt
    end function QL_func

    function QS_func(r,rho,Cp,T,T_before,dt) result(Q_S)
      real(dp),allocatable,intent(in)::r(:),rho(:),Cp(:),T(:),T_before(:)
      real(dp),intent(in)::dt
      real(dp)::Q_S,pi_
      real(dp),allocatable::qs(:)

      pi_ = acos(-1.0_dp)

      qs = r**2.0_dp*rho*Cp*(T-T_before)

      Q_S = midpoint_integration(r,qs)*4*pi_/dt
    end function QS_func

    function QP_func(r,alpha,T,PT,Ttc,Ttc_old,dt) result(Q_P)
      real(dp),allocatable,intent(in)::r(:),alpha(:),T(:)
      real(dp),intent(in)::PT,Ttc,Ttc_old,dt
      real(dp)::Q_P,pi_
      real(dp),allocatable::qp(:)

      pi_ = acos(-1.0_dp)

      qp = r**2.0_dp*alpha*T

      Q_P = midpoint_integration(r,qp)*4*PT*(Ttc-Ttc_old)*pi_/dt
    end function QP_func

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Initialize core parameters to be used in core cooling !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine init_core(params, core)
      type(core_properties),intent(inout)::core
      type(params_tot),intent(inout)::params
      real(dp),allocatable::vals_t(:,:) 
      real(dp),allocatable::g_init(:),p_init(:),rho_init(:),r_init(:),T_init(:),Cp_init(:),alpha_init(:),mat_init(:)
      integer::n_cmb,i,m
      
      ! length of arrays will depend on the core size
      n_cmb = params%pI%prof_nm+1 !cmb index (initialize.f90 line 585, 07.08.2024) !+1 and we start in mat:8 without we start in mat:7
      m = 1000 - n_cmb + 1 ! size of init arrays
      core%n = 15000 ! change core resolution (e.g. r_init(2)-r_init(-1))

      ! check that core%n is gt m
      if (core%n < m) then
        print *, "Error: The core resolution must be greater than the number of core data points."
        stop
      end if

      allocate(vals_t(1000,14),g_init(m),p_init(m),rho_init(m),r_init(m),T_init(m),Cp_init(m),alpha_init(m),mat_init(m),core%T_old(core%n))
      
      open(70,file="profs.res",form='formatted',status='unknown')
        do i=1,1000
          read(70, *) vals_t(i,:) ! surface values at n_t, read here from cmb to core
        enddo
      close(70)
      
      ! store and reverse interior structure arrays for interpolation
      g_init = reverse_array(vals_t(n_cmb:1000,1))
      p_init = reverse_array(vals_t(n_cmb:1000,2))*1.0e9_dp ! GPa to Pa
      rho_init = reverse_array(vals_t(n_cmb:1000,3))
      r_init = reverse_array(vals_t(n_cmb:1000,4)) 
      T_init = reverse_array(vals_t(n_cmb:1000,5))
      Cp_init = reverse_array(vals_t(n_cmb:1000,7))
      alpha_init = reverse_array(vals_t(n_cmb:1000,8))
      mat_init = reverse_array(vals_t(n_cmb:1000,14))
      
      ! initialize radii
      !---------------------- r_init(m+1) - mantle,  mat = 7.0
      !---------------------- r_init(m) - cmb,     mat = 8.0
      !---------------------- r_init(m-1) - core,   mat = 8.0 
      !          ...
      !---------------------- r_init(1) - core = 0.0,   mat = 8.0  
      core%ri = 0.0_dp
      core%r = linspace(0.0_dp,r_init(m),core%n) 
      core%rc = core%r(core%n)-core%r(1)
      core%dr = (core%rc)/real(core%n-1,dp) ! core radius step size

      ! initiallize const. parameters
      core%Hc = 1.0e-3_dp   ! (W/kg)
      core%LH = 750000.0_dp ! (J/kg) (Nimmo et al. 2015)
      core%PT = 1.0_dp !                                                      ToDo: find correct value
      core%kc = 125.07_dp ! W/m/K       !125.07 to 216.18 W/m/K (Li et al. 2021) from Nimmo et al. 2015 130 W/m/K is given

      ! build core arrays with cubic spline interpolation and core radius
      core%T_old = 0.0_dp
      core%g = cubic_spline_interpolation(r_init,g_init,core%r)
      core%p = cubic_spline_interpolation(r_init,p_init,core%r)
      core%rho = cubic_spline_interpolation(r_init,rho_init,core%r)
      core%T = cubic_spline_interpolation(r_init,T_init,core%r) 
      core%Cp = cubic_spline_interpolation(r_init,Cp_init,core%r)
      core%alpha = cubic_spline_interpolation(r_init,alpha_init,core%r)

      core%Ttc_old = core%T(core%n)+0.0004_dp ! initial guess
      core%Ttc = core%T(core%n-2) ! find correct value
      ! ToDo: define a more soffisticated way to calculate Tad of the previous step
      core%T_old = core%T+0.0004_dp ! initial guess     ! T = T_t exp( alpha_t g / Cp_t * dr )
      
      ! pure iron melting temperature (Sixtrude et al. 2014)
      ! core%Tm = 6500.0_dp * (core%p*1.0e-9_dp / 340.0_dp)**0.515_dp 
      ! pure iron melting temperature (Gonzalez-Cataldo & Militzer 2023)
      ! core%Tm = 6469.0_dp * (1.0_dp + (core%p*1.0e-9_dp - 300.0_dp)/434.82_dp)**0.54369_dp 

      ! pure iron melting temperature (Stevenson et al. 1983)
      core%Tm = 2060.0_dp * (1 + 6.14_dp*1.0e-12_dp * core%p - 4.5_dp*1.0e-24_dp * core%p**2.0_dp) 
      
      deallocate(vals_t,g_init,p_init,rho_init,r_init,T_init,Cp_init,alpha_init,mat_init)

      !ToDo: Do I need to set up and consider a core thermal boundary layer initially? 
      !ToDo: how and where do I save the mantle thermal boundary data?

      !reveal the two core radii Lena - Louis
      print *, "core%rc: ", core%rc
      print *, "dimentionalized --> params%pI%Rc: ", params%pI%Rc*params%pF%D
      print *, "what about: ", (params%pI%Rp-params%pI%Rc)*params%pF%D
    end subroutine init_core
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Core cooling depending on heat flux out of core into mantle !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine core_cooling(params,core,field,mesh,km_,dt_,t)
      type(core_properties),intent(inout)::core
      type(params_tot),intent(inout)::params
      type(variables_unknowns),intent(in)::field
      type(mesh_cp),intent(in)::mesh
      real(dp),intent(in)::km_,dt_,t
      real(dp)::rho_c,Cp_c,Vol,Tav                          ! local variables stored in and passed from params
      real(dp)::rb,rtc,Tb,Tc,km,rhom,rhoc,Hm,dt                ! dimensionalized variables
      real(dp)::cm0,cc0,cm1,cc1                             ! integration constants
      integer,parameter::n_la=5
      integer::ipiv(n_la),info,i_lo,i_hi,i,j                ! iteration vars
      double precision, dimension(n_la,n_la)::A             ! linear system variables
      double precision, dimension(n_la)::b

      ! redimentionalize time step variable
      dt = dt_/params%pF%time_yr*365.0_dp*24.0_dp*3600.0_dp ! s 

      ! Block 1 - Finding Inner Core Radius
      !--------------------------------------------------------------------------------------
      ! 1st check if the core is completely solid 
      if (core%T(core%n) < core%Tm(core%n)) then
        core%ri = core%r(core%n)
        core%rhoi = core%rho(core%n)
        print *, "Core is completely solid"
      else
        ! 2nd find where the melting temperature is reached
        i_lo = 1; i_hi = core%n
        do while (i_hi - i_lo > 1) 
          i = (i_hi + i_lo)/2
          if (core%T(i) > core%Tm(i)) then
            i_hi = i
          else
            i_lo = i
          end if
        end do
        ! 3rd then either the melting temperature is not reached
        if ((i < 3).and.(core%T(i) > core%Tm(i))) then
          core%ri = 0.0_dp
          core%rhoi = core%rho(1)
          print *, "Iron at the center of the core has not cooled enough to crystallize,"
          print *, "because core%T(n) = ", core%T(core%n)
          print *, "is greater than core%Tm(n) = ", core%Tm(core%n)
        else
          ! 4th or it is reached
          core%ri = core%r(i)
          core%rhoi = core%rho(i)
        ! then find the exact melting point by linear interpolation:
        !elseif ((core%T(i) < core%Tm(i)).and.(core%T(i+1) > core%Tm(i+1))) then
        !  core%ri = (core%r(i)+core%r(i+1))/2.0_dp ! linear interpolation
        !else
        !  core%ri = (core%r(i)+core%r(i-1))/2.0_dp ! s.a.
        end if
      end if
      
      print *, "core%T(1): ", core%T(1), "K"
      print *, "core%Tm(1): ", core%Tm(1), "K"
      print *, "core%T(n): ", core%T(core%n), "K"
      print *, "core%Tm(n): ", core%Tm(core%n), "K"
      print *, "Inner core radius lies at a depth of: ", core%ri, "m"
      print *, "Density at the inner core: ", core%rhoi, "kg/m^3"
      !--------------------------------------------------------------------------------------


      ! Block 2 - Calculate Energy Balance at the Core-Mantle Boundary
      !--------------------------------------------------------------------------------------
      core%ri_old = core%ri-0.5_dp

      ! radiogenic heat release
      core%QR = QR_func(core%r,core%rho,core%Hc,dt)
      ! latent heat release
      core%QL = QL_func(core%ri,core%ri_old,core%LH,core%rhoi,dt)
      ! Secular cooling of the core
      core%QS = QS_func(core%r,core%rho,core%Cp,core%T,core%T_old,dt)
      ! pressure heating
      core%QP = QP_func(core%r,core%alpha,core%T,core%PT,core%Ttc,core%Ttc_old,dt)
      ! total heat flux at the core-mantle boundary
      core%QCMB = core%QR + core%QL + core%QS + core%QP

      core%Tc = core%QCMB/(-core%kc*4*acos(-1.0_dp)*core%rc**2.0_dp)*core%dr + core%Ttc

      print *, "QR: ", core%QR*1.0e-12_dp, "TW"
      print *, "QL: ", core%QL*1.0e-12_dp, "TW"
      print *, "QS: ", core%QS*1.0e-12_dp, "TW"
      print *, "QP: ", core%QP*1.0e-12_dp, "TW"
      print *, "QCMB: ", core%QCMB*1.0e-12_dp, "TW"
      print *, "Tc: ", core%Tc, "K"
      !--------------------------------------------------------------------------------------

      ! Block 3 - save and update core properties
      !--------------------------------------------------------------------------------------
      core%T_old = core%T
      !...ToDo: update further core properties
      !--------------------------------------------------------------------------------------

      ! Block 3 - Updated Core Mantle Boundary Temperature Tc
      !--------------------------------------------------------------------------------------
      ! To find the correct sollution for Tc we have rb, rc, Tb, Tc, km, kc, rhom, Hm given.
      rb = core%r(core%n) ! find correct value
      rtc = core%r(core%n-2) ! find correct value
      Tb = core%T(core%n) ! find correct value
      km = params%pF%k
      rhom = 5536.0_dp ! ToDo: which rho (mantle) should be used here, if rho not constant in mantle? reference value, average value or value directly above CMB?
      rhoc = core%rho(core%n) ! Same question as above
      Hm = radiogenic_heating(params,params%pI%H0,params%pI%lambda,t)*params%pF%k*params%pF%DeltaT/params%pF%D**2.0_dp
      core%Hc = Hm*0.4_dp ! ToDo: find correct value for Hc, Hm* 0.4 is a placeholder value
      
      print *, "param%pF%k: ", params%pF%k
      print *,  "rb  	  rc"
      print *, rb,core%rc
      print *,  "Tb  	  Ttc"
      print *, Tb,core%Ttc
      print *, km,core%kc,rhom
      print *, "rhoc  	  Hm  	  Hc"
      print *, rhoc,Hm,core%Hc

      ! Placeholder values for now
      cm0=0.0_dp; cc0=0.0_dp; cm1=0.0_dp; cc1=0.0_dp; Tc=1.0_dp

      ! radial heat conduction in a sphere or sphereical shell 
      ! dq/dr + 2q/r = pho*H --> T = -rho*H/(6*k) + c1/r + c0

      !solve linear system of integration constants and core-mantle-boundary temperature
      !  1   2   3   4   5
      !(cm0,cc0,cm1,cc1,Tc) 
      A = reshape([1.0_dp, 0.0_dp, 1.0_dp/core%rc, 0.0_dp, -1.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, -core%rc*rb/(rb-core%rc), &
        0.0_dp, 1.0_dp, 1.0_dp/core%rc, 0.0_dp, -1.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, rtc*core%rc/(core%rc-rtc), &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [n_la, n_la])
  
      b(1) = rhom*Hm/(6.0_dp*km)*core%rc**2
      b(2) = rhom*Hm/(6.0_dp*km)*(core%rc**2*rb + core%rc*rb**2) - Tb*rb*core%rc/(rb-core%rc)
      b(3) = rhoc*core%Hc/(6.0_dp*core%kc)*core%rc**2
      b(4) = rhoc*core%Hc/(6.0_dp*core%kc)*(core%rc**2*rtc + core%rc*rtc**2) + core%Ttc*rtc*core%rc/(core%rc-rtc)
      b(5) = (Tb + core%Ttc*core%kc/km*(core%rc-rb)/(rtc-core%rc))/(1+core%Ttc*core%kc/km*(core%rc-rb)/(rtc-core%rc))
  
      ! Call gsev from LAPACK to solve the system A * x = b
      call dgesv(n_la, 1, A, n_la, ipiv, b, n_la, info)
  
      ! Check if the solution was successful
      if (info /= 0) then
        print *, "Error: The solution could not be computed."
        stop
      end if
      
      !B is overwritten by the result in dgesv() 
      cm0 = b(1)
      cc0 = b(2)
      cm1 = b(3)
      cc1 = b(4)
      Tc = b(5) !new temperature at the core-mantle-boundary
      !-------------------------------------------------------------------------------------
      
      rho_c=params%pI%rho_c
      Cp_c=params%pI%Cp_c
  
      Tav = 0.0_dp
      Vol = 0.0_dp
      do i=1,mesh%nl
        do j=1,mesh%ny
          Tav = Tav + mesh%dVi(i,j,1)*( field%T(i,j,1) - params%pT%Tbottom )/(mesh%rc(1)-mesh%rmin)
          Vol = Vol + mesh%dVi(i,j,1)
        enddo
      enddo
      !write(46,*)Tav,Vol
      Tav = Tav/Vol
  
      !Tbottom = Tbottom + 3.0_dp*(rho*Cp/(rho_c*Cp_c))*dt*Tav / mesh%rmin
      params%pT%Tbottom = params%pT%Tbottom + 3.0_dp*dt_*km_*Tav / (mesh%rmin*rho_c*Cp_c) ! this is how it should be, -km*dT/dr is heat flux, rho_c*Cp_c is multiplied to V_c -> version above simplification for kappa=const
      !write(46,*)"New Tbottom:",Tbottom
      !write(75,*)

      !print *, "mesh%rc(1): ", mesh%rc(1)*params%pF%D
      !print *, "mesh%rmin: ", mesh%rmin*params%pF%D
      print *, "Tav: ", Tav*params%pF%DeltaT + params%pF%Ts, "K"
      print *, "calculated Tc: ", Tc, "K"
      print *, "new Tbottom: ", params%pT%Tbottom*params%pF%DeltaT + params%pF%Ts, "K"

    end subroutine core_cooling
    
    !ToDo: implement core_write subroutine in main.f90 #!16.8.24#
    !-------------------------------------------------------------------------------------
    subroutine core_write(core,params)
      type(core_properties),intent(in)::core
      type(params_tot),intent(in)::params
      integer::n_cmb, k

      n_cmb = size(core%g)

      if (n_cmb /= 1000-params%pI%prof_nm+1) then 
        print *, "Error: The integer n_cmb does not match the length of the loaded core data."
        stop
      end if

      ! ToDo: write resulting core energy data here #!7.8.24#

      write(75, *) 'Gravity [m/s] ','Pressure [GPa] ', 'Density [kg/m^3] ', 'Radius [m] ', &
                  'Temperature [K] ', 'Melt Temperature [K] ', 'Specific Heat Capacity [J/(K*kg)] ', &
                  'Thermal expansion [K^-1] ', 'Material Phase Number '

      do k = 1, n_cmb
        write(75, '(9(f35.15))') core%g(k),core%p(k),core%rho(k),core%r(k),core%T(k),core%Tm(k),&
                              core%Cp(k),core%alpha(k),core%mat(k)
      end do
    end subroutine core_write
    !-------------------------------------------------------------------------------------
  
end module core