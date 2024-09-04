module core
  use precision
  use mesher, only:mesh_cp
  use variables
  use parameters_total, only:params_tot 
  use heating
  
  implicit none

  type core_properties
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
    real(dp)::kappa  !m^2/s                     ! thermal diffusivity
    real(dp)::beta                              ! saturation constant for fast rotating dynamos
    real(dp)::mu  !H/m                          ! magnetic permeability
    real(dp)::gi  !m/s^2                        ! gravity at the inner core
    real(dp)::pi  !Pa                           ! pressure at the inner core
    real(dp)::Ra                                ! Rayleigh number
    real(dp)::Ra_crit                           ! critical Rayleigh number
    real(dp)::delta                             ! thickness of the thermal boundary layer
    real(dp)::eta                               ! viscosity of the thermal boundary layer
    real(dp)::Tb  !K                            ! temperature above core-mantle boundary
    real(dp)::Tc  !K                            ! temperature at the core-mantle boundary
    real(dp)::dTc !K                            ! temperature at the core-mantle boundary at the previous time step
    real(dp)::Ti  !K                            ! pure iron melting temperature
    real(dp)::dTm_dp !K/Pa                      ! derivative of the pure iron melting temperature with respect to pressure
    real(dp)::dT_dp !K/Pa                       ! derivative of the core temperature with respect to pressure
    real(dp)::FT  !W/m^2                        ! thermal buoyancy flux
    real(dp)::m  !A*m^2                         ! magnetic moment
    real(dp)::Bc  !T                            ! magnetic field strength at the core-mantle boundary
    real(dp)::Bs  !T                            ! magnetic field strength at the surface
    real(dp)::TD,TR,TS,TP  !K                             
    real(dp)::QR,QL,QS,QP,QC,QT,QD  ! W         ! heat fluxes (out of the core)                
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
      real(dp)::pi
      integer::n_cmb,i,m
      
      ! length of arrays will depend on the core size
      n_cmb = params%pI%prof_nm+1 !cmb index (initialize.f90 line 585, 07.08.2024) !+1 and we start in mat:8 without we start in mat:7
      m = 1000 - n_cmb + 1 ! size of init arrays
      core%n = 5000 ! change core resolution (e.g. r_init(2)-r_init(-1))
      pi = acos(-1.0_dp)

      ! check that core%n is gt m
      if (core%n < m) then
        print *, "Error: The core resolution must be greater than the number of core data points."
        stop
      end if

      allocate(vals_t(1000,14),g_init(m),p_init(m),rho_init(m),r_init(m),T_init(m),Cp_init(m),alpha_init(m),mat_init(m))
      
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
      alpha_init = reverse_array(vals_t(n_cmb:1000,8))*1.0e-5_dp ! 1/K
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
      core%Hc = 0.0_dp  ! (W/kg)
      core%LH = 750000.0_dp ! (J/kg) (Nimmo et al. 2015)
      core%PT = 0.0_dp !                                                      ToDo: find correct value
      core%kc = 125.07_dp ! W/m/K       !125.07 to 216.18 W/m/K (Li et al. 2021) from Nimmo et al. 2015 130 W/m/K is given
      core%gi = 0.0_dp ! (m/s^2)     
      core%eta = 1.0e+5_dp ! Pas         de Wijs et al. 1998 is 0.006 Pas but this value changes with 14 magnitudes in the literature so what to choose? 
      core%beta = 0.2_dp ! saturation constant for fast rotating dynamos (Bonati et al. 2021)
      core%mu = 4.0_dp*pi*1.0e-7_dp ! (H/m)           magnetic permeability 

      ! build core arrays with cubic spline interpolation and core radius
      core%g = cubic_spline_interpolation(r_init,g_init,core%r)
      core%p = cubic_spline_interpolation(r_init,p_init,core%r)
      core%rho = cubic_spline_interpolation(r_init,rho_init,core%r)
      core%T = cubic_spline_interpolation(r_init,T_init,core%r) 
      core%Cp = cubic_spline_interpolation(r_init,Cp_init,core%r)
      core%alpha = cubic_spline_interpolation(r_init,alpha_init,core%r)

      core%Tb = params%pT%Tbottom*params%pF%DeltaT + params%pF%Ts !ToDo check when it is defined?
      core%Tc = core%T(core%n-1)
      core%Tm = 2060.0_dp * (1 + 6.14_dp*1.0e-12_dp * core%p - 4.5_dp*1.0e-24_dp * core%p**2.0_dp)
      core%kappa = core%kc/(core%rho(core%n)*core%Cp(core%n)) ! thermal diffusivity
      
      deallocate(vals_t,g_init,p_init,rho_init,r_init,T_init,Cp_init,alpha_init,mat_init)
    end subroutine init_core
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Core cooling depending on heat flux out of core into mantle !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine core_cooling(params,core,field,mesh,km_,dt_,t_)
      type(core_properties),intent(inout)::core
      type(params_tot),intent(inout)::params
      type(variables_unknowns),intent(in)::field
      type(mesh_cp),intent(in)::mesh
      real(dp),intent(in)::km_,dt_,t_
      real(dp)::rho_c,Cp_c,Vol,Tav,dt,t,pi                        ! local variables stored in and passed from params
      character(len=27) :: headers(13)                            ! headers for file 75
      integer::i_lo,i_hi,i,j                                      ! iteration vars

      ! redimentionalize time step variable
      dt = dt_/params%pF%time_yr  !yr
      t = t_/params%pF%time_yr !yr
      dt = dt*365.0_dp*24.0_dp*3600.0_dp ! time step in s
      pi = acos(-1.0_dp)

      if (t < 1050) then              ! ToDo: find better statement for first time step
        ! Define the headers, padded with spaces to fit the format
        headers(1) = 't[Ma]                     '
        headers(2) = 'ri[m]                     '
        headers(3) = 'Ra[]                      '
        headers(4) = 'Tb[K]                     '
        headers(5) = 'Tc[K]                     '
        headers(6) = 'delta[m]                  '
        headers(7) = 'Ti[K]                     '
        headers(8) = 'QC[TW]                    '
        headers(9) = 'QS[TW]                    '
        headers(10)= 'QL[TW]                    '
        headers(11)= 'QD[TW]                    '
        headers(12)= 'Bc[T]                     '
        headers(13)= 'Bs[T]                     '
        
        write(75, '(13(A27))') headers
      end if
      
      ! Block 1 - Finding Inner Core Radius
      !--------------------------------------------------------------------------------------
      ! 1st check if the core is completely solid 
      if (core%T(core%n) < core%Tm(core%n)) then          ! remember center at r(1) and cmb at r(n)
        core%ri = core%r(core%n)
        core%gi = core%rho(core%n)
        print *, "Core is completely solid"
        print *, "core%T(n): ", core%T(core%n), "K"
        print *, "core%Tm(n): ", core%Tm(core%n), "K"
        stop

      ! 2nd check if the solidus is not reached anywhere in the core
      elseif ((core%T(1) > core%Tm(1))) then                     
          core%ri = 0.0_dp
          core%gi = core%g(1)
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
        core%gi = core%g(i_lo) + (core%g(i_hi) - core%g(i_lo)) * (core%T(i_lo) - core%Tm(i_lo)) / ((core%Tm(i_hi) - core%Tm(i_lo)) - (core%T(i_hi) - core%T(i_lo)))
        core%pi = core%p(i_lo) + (core%p(i_hi) - core%p(i_lo)) * (core%T(i_lo) - core%Tm(i_lo)) / ((core%Tm(i_hi) - core%Tm(i_lo)) - (core%T(i_hi) - core%T(i_lo)))
        core%Ti = core%T(i_lo) + (core%T(i_hi) - core%T(i_lo)) * (core%ri - core%r(i_lo)) / (core%r(i_hi) - core%r(i_lo))

        ! 5th calculate the derivative of melting and adiabatic temperature at the ICB
        core%dTm_dp = ((core%Tm(i_hi) - core%Tm(i_lo)) / (core%p(i_hi) - core%p(i_lo))) ! add weighted difference for non-uniform grid P if n low?
        core%dT_dp = ((core%T(i_hi) - core%T(i_lo)) / (core%p(i_hi) - core%p(i_lo)))  
      end if
      !--------------------------------------------------------------------------------------

      !Block 2 - Update Core Mantle Boundary Temperature Tb (passed to mantle)
      !--------------------------------------------------------------------------------------
      rho_c=params%pI%rho_c ! ToDo: define better values as in init_core for rho and Cp
      Cp_c=params%pI%Cp_c
  
      Tav = 0.0_dp
      Vol = 0.0_dp
      do i=1,mesh%nl
        do j=1,mesh%ny
          Tav = Tav + mesh%dVi(i,j,1)*( field%T(i,j,1) - params%pT%Tbottom )/(mesh%rc(1)-mesh%rmin)
          Vol = Vol + mesh%dVi(i,j,1)
        end do
      end do
      Tav = Tav/Vol
  
      ! Tbottom = Tbottom + 3.0_dp*(rho*Cp/(rho_c*Cp_c))*dt*Tav / mesh%rmin
      ! this is how it should be, -km*dT/dr is heat flux, rho_c*Cp_c is multiplied to V_c -> version above simplification for kappa=const
      params%pT%Tbottom = params%pT%Tbottom + 3.0_dp*dt_*km_*Tav / (mesh%rmin*rho_c*Cp_c)
      core%Tb = params%pT%Tbottom*params%pF%DeltaT + params%pF%Ts
      !--------------------------------------------------------------------------------------

      ! Block 3 - CMB Temp. Tc derived from Energy Budget
      !--------------------------------------------------------------------------------------
      !Define the rayleigh number
      core%Ra = core%alpha(core%n) * core%g(core%n) * (core%Ti-core%T(core%n)) * (core%rc-core%ri)**3.0_dp / (core%eta*core%kappa)
      core%Ra_crit = 0.28_dp * core%Ra**0.21_dp ! Breuer et al. 2010
      core%delta = (core%rc-core%ri) * (core%Ra_crit/core%Ra)**(1.0_dp/3.0_dp)
      
      ! Secular cooling of the core
      core%QS = -midpoint_integration(core%r,core%r**2.0_dp*core%rho*core%Cp*core%T)/core%Tc
      
      ! latent heat release
      core%QL = core%ri**2.0_dp*core%LH*core%Ti/((core%dTm_dp-core%dT_dp)*core%gi*core%Tc)

      ! radiogenic heat release
      core%QR = midpoint_integration(core%r,core%rho*core%r**2.0_dp)*core%Hc

      ! pressure heating
      core%QP = midpoint_integration(core%r,core%alpha*core%T*core%r**2.0_dp)*core%PT

      ! total heat flux at the core-mantle boundary
      core%QC = -params%pF%k*(core%Tc-core%Tb)/core%delta*core%rc**2.0_dp

      ! total heat flux at the core-mantle boundary after eq. 78 
      !in Nimmo et al. 2015 an update of Gubbins et al. 2003: dTc/dt = (Qcmb-QR)/~QT
      core%Tc = core%Tc + (core%QC-core%QR)/(core%QS+core%QL+core%QP)*dt
      !--------------------------------------------------------------------------------------

      ! Block 4 - Update core temperature profile
      !--------------------------------------------------------------------------------------
      ! check convergence
      if ((core%Tc > 1.0e+10_dp).or.(core%Tb > 1.0e+10_dp)) then
        print *, "Error: CMB temperature is too high and is likely diverging."
        stop
      elseif ((core%Tc < 0.0_dp).or.(core%Tb < 0.0_dp)) then
        print *, "Error: CMB temperature is too cold."
        stop
      end if

      core%dTc = abs(core%T(core%n) - core%Tc) / dt !change to last time step dTc/dt
      
      core%T(core%n) = core%Tc ! updated cmb boundary cond.

      do i = core%n, 2, -1
        core%T(i-1) = core%T(i) + core%alpha(i) * core%g(i) / core%Cp(i) * core%T(i) * core%dr
      end do

      ! ToDo: update viscosity based on temperature?
      ! eta = eta0 * e^(2.6*R*Tm(p)/(kB * T))
      
      !...ToDo: update further core properties?
      !--------------------------------------------------------------------------------------

      ! Block 5 - Find the Dispersion Energy
      !--------------------------------------------------------------------------------------
      ! heat output units [W]
      core%QC = core%QC*4.0_dp*pi
      core%QS = core%QS*4.0_dp*pi*core%dTc
      core%QL = core%QL*4.0_dp*pi*core%dTc
      core%QR = core%QR*4.0_dp*pi
      core%QP = core%QP*4.0_dp*pi*core%dTc

      ! speceific tempeatures
      core%TD = core%T(core%n-2)
      core%TS = core%T(core%n-1)
      core%TR = core%T(core%n-1)
      core%TP = core%T(core%n-1)

      core%QD = core%TD/core%Tc * ((1-core%Tc/core%TS)*core%QS + (1-core%Tc/core%Ti)*core%QL &
                 + (1-core%Tc/core%TR)*core%QR + (1-core%Tc/core%TP)*core%QP)
      !--------------------------------------------------------------------------------------

      ! Block 6 - Define the magnetic field strength
      !--------------------------------------------------------------------------------------
      ! after Bonati et al. 2021
      ! thermal buoyancy flux =>    FT = alpha*g*Qc/(rho*Cp)
      core%FT = core%alpha(core%n)*core%g(core%n)/(core%rho(core%n)*core%Cp(core%n))*(-core%QC/(4.0_dp*pi*core%rc**2.0_dp))

      ! magnetic moment =>    m = 4*pi*rc^3*beta*(rho/mu)*FT^(1/3) 
      core%m = 4.0_dp*pi*core%rc**3.0_dp*core%beta*(core%rho(core%n)/core%mu)*(core%FT*(core%rc-core%ri))**(1.0_dp/3.0_dp)

      ! magnetic field strength at cmb =>    Bc = beta*(rho*mu)^(1/2)*FT^(1/3)*(rc-ri)^(1/3)
      core%Bc = core%beta*(core%rho(core%n)*core%mu)**0.5_dp*(core%FT*(core%rc-core%ri))**(1.0_dp/3.0_dp)

      ! magnetic field strength at the surface =>    Bs = Bc*(rc/rp)^(3)
      core%Bs = core%Bc*(core%rc/(params%pF%D*params%pI%Rp))**3.0_dp
      !--------------------------------------------------------------------------------------

      ! Block 7 - Write Core Properties to File 75
      !-------------------------------------------------------------------------------------
      ! bring the first flux values in order for visibility when plotting:
      if (t < 10000) then
        core%QS = 0.0_dp
        core%QL = 0.0_dp
        core%QR = 0.0_dp
        core%QP = 0.0_dp
        core%QD = 0.0_dp
      end if

      ! heat output units [TW]
      write(75, '(13(x e23.15,x) )') t*1.0e-6_dp,core%ri,core%Ra,core%Tb,core%Tc,core%delta, &
                                    core%Ti,core%QC*1.0e-12_dp,core%QS*1.0e-12_dp,core%QL*1.0e-12_dp, &
                                    core%QD*1.0e-12_dp,core%Bc,core%Bs
      !-------------------------------------------------------------------------------------
    end subroutine core_cooling
end module core