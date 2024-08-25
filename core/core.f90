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
    real(dp)::Tb  !K                            ! temperature above core-mantle boundary
    real(dp)::Tc    !K                          ! temperature at the core-mantle boundary
    real(dp)::Ttc  !K                           ! temperature at the bottom of the core TBL
    real(dp),allocatable::T_old(:) !K  
    real(dp)::Tb_old !K
    real(dp)::Tc_old  !K                      
    real(dp)::ri_old !m    
    real(dp)::QR,QL,QS,QP,QCMB  ! W             ! heat fluxes (out of the core)                
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

    function QS_func(r,rho,Cp,T,T_before,dt) result(Q_S)
      real(dp),allocatable,intent(in)::r(:),rho(:),Cp(:),T(:),T_before(:)
      real(dp),intent(in)::dt
      real(dp)::Q_S,pi_
      real(dp),allocatable::qs(:)

      pi_ = acos(-1.0_dp)

      qs = r**2.0_dp*rho*Cp*(T-T_before)

      Q_S = midpoint_integration(r,qs)*4*pi_/dt
    end function QS_func

    function QL_func(ri,ri_before,L_H,rhoi,dt) result(Q_L)
      real(dp),intent(in)::L_H,ri_before,ri,rhoi,dt
      real(dp)::Q_L, pi_, dri

      pi_ = acos(-1.0_dp)

      dri = ri - ri_before

      Q_L = 4*pi_*ri**2.0_dp*L_H*rhoi*dri/dt
    end function QL_func

    function QR_func(r,rho,H,dt) result(Q_R)
      real(dp),allocatable,intent(in)::r(:),rho(:)
      real(dp),intent(in)::H,dt
      real(dp),allocatable::qr(:)
      real(dp)::Q_R, pi_

      pi_ = acos(-1.0_dp)

      qr = rho*H*r**2

      Q_R = midpoint_integration(r,qr)*4*pi_/dt
    end function QR_func

    function QP_func(r,alpha,T,PT,Ttc,Ttc_old,dr,dt) result(Q_P)
      real(dp),allocatable,intent(in)::r(:),alpha(:),T(:)
      real(dp),intent(in)::PT,Ttc,Ttc_old,dr,dt
      real(dp)::Q_P,pi_
      real(dp),allocatable::qp(:)

      pi_ = acos(-1.0_dp)

      qp = r**2.0_dp*alpha*T

      Q_P = midpoint_integration(r,qp)*4*PT*pi_*(Ttc-Ttc_old)/(dt*dr)
    end function QP_func

    function T_melt(P) result(T_m)
      real(dp),intent(in)::P
      real(dp)::T_m
      ! pure iron melting temperature (Sixtrude et al. 2014)
      ! core%Tm = 6500.0_dp * (core%p*1.0e-9_dp / 340.0_dp)**0.515_dp 
      ! pure iron melting temperature (Gonzalez-Cataldo & Militzer 2023)
      ! core%Tm = 6469.0_dp * (1.0_dp + (core%p*1.0e-9_dp - 300.0_dp)/434.82_dp)**0.54369_dp
      ! pure iron melting temperature (Stevenson et al. 1983)
      T_m = 2060.0_dp * (1 + 6.14_dp*1.0e-12_dp * P - 4.5_dp*1.0e-24_dp * P**2.0_dp)
    end function T_melt

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
      core%n = 20000 ! change core resolution (e.g. r_init(2)-r_init(-1))

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
      core%Hc = 1.0e-12_dp  ! (W/kg)
      core%LH = 750000.0_dp ! (J/kg) (Nimmo et al. 2015)
      core%PT = 1.0_dp !                                                      ToDo: find correct value
      core%kc = 125.07_dp ! W/m/K       !125.07 to 216.18 W/m/K (Li et al. 2021) from Nimmo et al. 2015 130 W/m/K is given

      ! build core arrays with cubic spline interpolation and core radius
      core%g = cubic_spline_interpolation(r_init,g_init,core%r)
      core%p = cubic_spline_interpolation(r_init,p_init,core%r)
      core%rho = cubic_spline_interpolation(r_init,rho_init,core%r)
      core%T = cubic_spline_interpolation(r_init,T_init,core%r) 
      core%Cp = cubic_spline_interpolation(r_init,Cp_init,core%r)
      core%alpha = cubic_spline_interpolation(r_init,alpha_init,core%r)

      core%Tb = params%pT%Tbottom*params%pF%DeltaT + params%pF%Ts !ToDo check when it is defined?
      core%ri_old = 0.0_dp ! initial guess
      core%Ttc = core%T(core%n-1) ! find correct value
      core%Tc_old = core%T(core%n)+0.00004_dp ! initial guess
      core%Tc = core%T(core%n) ! find correct value
      core%T_old = core%T+0.00004_dp ! initial guess     ! T = T_t exp( alpha_t g / Cp_t * dr )
      core%Tm = 2060.0_dp * (1 + 6.14_dp*1.0e-12_dp * core%p - 4.5_dp*1.0e-24_dp * core%p**2.0_dp)
      
      !headers = ['dt  ','ri  ','rhoi','Tc  ','Tb  ','QCMB','QS  ','QL  ','QR  ','QP  ']
      ! initiate header of file 75 for core properties
      !write(75, '(9(A4))') headers
      
      deallocate(vals_t,g_init,p_init,rho_init,r_init,T_init,Cp_init,alpha_init,mat_init)
    end subroutine init_core
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Core cooling depending on heat flux out of core into mantle !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine core_cooling(params,core,field,mesh,km_,dt_,dt_old_,t_)
      type(core_properties),intent(inout)::core
      type(params_tot),intent(inout)::params
      type(variables_unknowns),intent(in)::field
      type(mesh_cp),intent(in)::mesh
      real(dp),intent(in)::km_,dt_,dt_old_,t_
      real(dp)::rho_c,Cp_c,Vol,Tav,dt,dt_old,t                    ! local variables stored in and passed from params
      character(len=27) :: headers(10)                           ! headers for file 75
      !real(dp)::rb,rtc,Tb,Tc,km,rhom,rhoc,Hm               ! dimensionalized variables
      !real(dp)::cm0,cc0,cm1,cc1                            ! integration constants
      !integer,parameter::n_la=5
      !integer::ipiv(n_la),info
      integer::i_lo,i_hi,i,j,check                              ! iteration vars
      !double precision, dimension(n_la,n_la)::A            ! linear system variables
      !double precision, dimension(n_la)::b

      !view results 1: 
      check = 0

      ! redimentionalize time step variable
      dt = dt_/params%pF%time_yr  !yr
      dt_old = dt_old_/params%pF%time_yr !yr
      t = t_/params%pF%time_yr !yr

      if (t < 1050) then              ! ToDo: find better statement for first time step
        dt = dt*365.0_dp*24.0_dp*3600.0_dp ! time step in s

        ! Define the headers, padded with spaces to fit the format
        headers(1) = 't [yr]                    '
        headers(2) = 'ri [m]                    '
        headers(3) = 'rhoi [kg/m^3]             '
        headers(4) = 'Tc [K]                    '
        headers(5) = 'Tb [K]                    '
        headers(6) = 'QCMB [TW]                 '
        headers(7) = 'QS  [TW]                  '
        headers(8) = 'QL  [TW]                  '
        headers(9) = 'QR  [TW]                  '
        headers(10)= 'QP  [TW]                  '
        write(75, '(10(A27))') headers
      else
        dt = (dt-dt_old)*365.0_dp*24.0_dp*3600.0_dp ! time step in s
      end if
      
      ! Block 1 - Finding Inner Core Radius
      !--------------------------------------------------------------------------------------
      ! 1st check if the core is completely solid 
      if (core%T(core%n) < core%Tm(core%n)) then          ! remember center at r(1) and cmb at r(n)
        core%ri = core%r(core%n)
        core%rhoi = core%rho(core%n)
        print *, "Core is completely solid"
        print *, "core%T(n): ", core%T(core%n), "K"
        print *, "core%Tm(n): ", core%Tm(core%n), "K"
      ! 2nd check if the solidus is not reached anywhere in the core
      elseif ((core%T(1) > core%Tm(1))) then                     
          core%ri = 0.0_dp
          core%rhoi = core%rho(1)
          print *, "Iron at the center of the core has not cooled enough to crystallize,"
          print *, "because core%T(n) = ", core%T(1)
          print *, "is greater than core%Tm(n) = ", core%Tm(1)
      ! 3rd then find the interval where the melting temperature is reached (binary search)
      else       
        i_lo = 1; i_hi = core%n
        do while (i_hi - i_lo > 1) 
          i = (i_hi + i_lo)/2
          if ((core%T(i_lo) - core%Tm(i_lo) < 0.0_dp).and.(core%T(i)-core%Tm(i) > 0.0_dp)) then !confines of the intersection point
            i_hi = i
          else
            i_lo = i
          end if
        end do 
        ! 4th use linear interpolation to find the exact radius and density                               !ToDo: check maximum error here for decreasing n
        core%ri = core%r(i_lo) + (core%r(i_hi) - core%r(i_lo)) * (core%T(i_lo) - core%Tm(i_lo)) / ((core%Tm(i_hi) - core%Tm(i_lo)) - (core%T(i_hi) - core%T(i_lo)))
        core%rhoi = core%rho(i_lo) + (core%rho(i_hi) - core%rho(i_lo)) * (core%T(i_lo) - core%Tm(i_lo)) / ((core%Tm(i_hi) - core%Tm(i_lo)) - (core%T(i_hi) - core%T(i_lo)))
      end if
      !--------------------------------------------------------------------------------------

      ! Block 2 - Update Core Mantle Boundary Temperature Tc (passed from  core)
      !--------------------------------------------------------------------------------------
       !here the very first timestep calculation is done with initial values set in core_init
      core%ri_old = core%ri - 0.0001_dp !m !initial guess
      
      ! Secular cooling of the core
      core%QS = QS_func(core%r,core%rho,core%Cp,core%T,core%T_old,dt)/1.0e12_dp
      ! latent heat release
      core%QL = QL_func(core%ri,core%ri_old,core%LH,core%rhoi,dt)/1.0e12_dp
      ! radiogenic heat release
      core%QR = QR_func(core%r,core%rho,core%Hc,dt)/1.0e12_dp
      ! pressure heating
      core%QP = QP_func(core%r,core%alpha,core%T,core%PT,core%Tc,core%Tc_old,core%dr,dt)/1.0e12_dp ! W to TW	
      ! total heat flux at the core-mantle boundary
      core%QCMB = core%QR + core%QL + core%QS + core%QP

      ! => QCMB = -kc*(Ttc-Tc)/dr*4*pi*rc^2  and  QCMB = QR+QL+QS+QP

      ! => Tc = dr/(4*pi*rc^2)*(QR+QL+QS+QP) + Ttc

      core%Tc = 1.0e12*core%QCMB/(core%kc*4*acos(-1.0_dp)*core%rc**2.0_dp)*core%dr + core%Ttc
      !--------------------------------------------------------------------------------------

      !Block 3 - Update Core Mantle Boundary Temperature Tb (passed to mantle)
      !--------------------------------------------------------------------------------------
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
      Tav = Tav/Vol
  
      ! Tbottom = Tbottom + 3.0_dp*(rho*Cp/(rho_c*Cp_c))*dt*Tav / mesh%rmin
      ! this is how it should be, -km*dT/dr is heat flux, rho_c*Cp_c is multiplied to V_c -> version above simplification for kappa=const
      params%pT%Tbottom = params%pT%Tbottom + 3.0_dp*dt_*km_*Tav / (mesh%rmin*rho_c*Cp_c)
      !--------------------------------------------------------------------------------------

      ! Block 4 - Update core properties
      !--------------------------------------------------------------------------------------
      ! save old values
      core%ri_old = core%ri
      core%T_old = core%T
      core%Tc_old = core%Tc
      core%Tb_old = core%Tb

      ! store new values
      do i = 1, core%n
        core%T(i) = core%T(i) + params%pT%Tbottom*params%pF%DeltaT + params%pF%Ts - core%Tc   ! adiabaic core cooling dependent on heat flux out of core
      end do
      core%Tb = params%pT%Tbottom*params%pF%DeltaT + params%pF%Ts
      core%Tc = core%T(core%n)
      core%Ttc = core%T(core%n-1)
      !...ToDo: update further core properties?
      !--------------------------------------------------------------------------------------

      ! BLock 5 - Convergence Check
      !--------------------------------------------------------------------------------------
      if (core%Tc > 1.0e+24) then
        print *, "Error: CMB temperature is too high and is likely diverging."
        stop
      elseif (core%Tb > 1.0e+24) then
        print *, "Error: Temperature at the bottom of the mantle is too high and is likely diverging."
        stop
      elseif (core%Tc < 0.0_dp) then
        print *, "Error: CMB temperature is too cold."
        stop
      elseif (core%Tb < 0.0_dp) then
        print *, "Error: Temperature at the bottom of the mantle is too cold."
        stop
      end if
      !--------------------------------------------------------------------------------------

      ! Block 6 - Write core properties to file 75 and print results
      !-------------------------------------------------------------------------------------
      if (check > 0) then
        print *, "times: ", dt_/params%pF%time_yr, "yr"
        print *, "QR: ", core%QR, "TW", " QL: ", core%QL, "TW", " QS: ", core%QS, "TW", " QP: ", core%QP, "TW"
        print *, "QCMB: ", core%QCMB, "TW"
        print *, "After the calculation of the core cooling subroutine:"
        print *, "old CMB Temp.: ", core%Tc_old, "K"
        print *, "new CMB Temp.: ", core%Tc, "K"
        print *, "old Tbottom: ", core%Tb_old, "K"
        print *, "new Tbottom: ", params%pT%Tbottom*params%pF%DeltaT + params%pF%Ts, "K"
        print *, "new bottom core TBL Temp.: ", core%Ttc, "K"
      end if

      write(75, '(10(x e23.15,x) )') t,core%ri,core%rhoi,core%Tc,core%Tb,core%QCMB,core%QS,core%QL,core%QR,core%QP
      !-------------------------------------------------------------------------------------
    end subroutine core_cooling
end module core