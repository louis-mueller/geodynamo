! Block 5 - Updated Core Mantle Boundary Temperature Tc
      !--------------------------------------------------------------------------------------
      ! ToDo: find the correct definitions of rb, rc, Tb, Ttc, km, kc, rhom, Hm, rhoc, Hc.
      rb = core%r(core%n)+core%dr ! find correct value
      rtc = core%r(core%n-1) ! find correct value
      km = params%pF%k
      rhom = 5536.0_dp ! ToDo: which rho (mantle) should be used here, if rho not constant in mantle? reference value, average value or value directly above CMB?
      rhoc = core%rho(core%n) ! Same question as above
      Hm = radiogenic_heating(params,params%pI%H0,params%pI%lambda,t)*params%pF%k*params%pF%DeltaT/params%pF%D**2.0_dp
      core%Hc = Hm*0.4_dp ! ToDo: find correct value for Hc, Hm* 0.4 is a placeholder value

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
      print *, "found Tc: ", Tc, "K"
      !-------------------------------------------------------------------------------------

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