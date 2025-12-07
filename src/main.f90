program stevia 
    use constants, only: dp
    use stellar_evolution
    ! use random_uniform
    implicit none
    integer, parameter :: nsteps = 100
    integer, parameter :: nstars = 5
    integer :: i, step
    
    real(dp) :: masses(nstars)
    real(dp) :: lifetimes(nstars)
    real(dp) :: time, dt
    real(dp) :: tend
    real(dp) :: current_mass = 0.0, current_teff = 0.0, current_lum = 0.0, current_mdot = 0.0, current_gamma_ed = 0.0, current_v_term = 0.0
    character(len=60) :: output_directory = "stellar_evolution_output/"
    character(len=32) :: mass_str
    
    ! Define stellar masses in solar masses
    masses = [1.0_dp, 5.0_dp, 10.0_dp, 20.0_dp, 50.0_dp]
    lifetimes = 0.0_dp
    time = 0.0_dp

    ! Calculate lifetimes for each star
    call load_stellar_data()
    do i = 1, nstars
        call SN_time(masses(i), lifetimes(i))
        write(*,'(A,F6.2,A,F8.2,A)') 'Star mass: ', masses(i), ' M_sun, Lifetime: ', lifetimes(i), ' Myr'
    enddo  
    ! Create output directory if it doesn't exist
    call execute_command_line('mkdir -p '//trim(output_directory), wait=.true.)

    open(unit=10, file=trim(output_directory)//'stellar_evolution_output.txt', status='replace')

    write(10,*) '# Stellar Evolution Output for masses: ', (masses(i), i=1,nstars)
    write(10,*) '# Mass (M_sun), Time (Myr), Current Mass (M_sun), Teff (K), Luminosity (L_sun), Mdot (M_sun/yr), Vterm (km/s)'
    
    ! Evolve stars over time accorse their rescpective lifetimes and print information at each time step
    do i = 1, nstars
        write(*, '(F6.2)') masses(i)
        write(mass_str,'(F6.2)') masses(i)
        open(unit=20+i, file=trim(output_directory)//'stellar_evolution_mass_'//trim(adjustl(mass_str))//'.txt', status='replace')
        write(20+i,*) '# Stellar Evolution Output for mass [M_sun] and lifetime [Myr]: '
        ! write(20+i,*) masses(i), ' M_sun'
        write(20+i,*) masses(i), lifetimes(i)
        ! write(20+i,*) '# Stellar Evolution Data for mass: ', masses(i), ' M_sun'
        write(20+i,*) '# steps: '
        write(20+i,*) nsteps
        write(20+i,*) '# Time (Myr), Current Mass (M_sun), Teff (K), Luminosity (L_sun), Mdot (M_sun/yr), Vterm (km/s)'
    enddo 
    do step = 0, nsteps
        write(*,*) 'Step: ', step
        write(10,*) 'Step: ', step
        do i = 1, nstars
            tend = lifetimes(i)
            dt = tend / nsteps
            time = step * dt
            call evolve_star(masses(i), time, current_mass, current_teff, current_lum, current_mdot, current_gamma_ed, current_v_term)
            write(*,'(A,F6.2,A,F8.2,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,F6.3)') &
             '  Star mass: ', masses(i), ' M_sun, Time: ', time, ' Myr, Mass: ', current_mass, ' M_sun, Teff: ', &
             current_teff, ' K, Luminosity: ', current_lum, ' L_sun, Mdot: ', current_mdot, ' M_sun/yr, Vterm: ', &
             current_v_term/kms, ' km/s, Gamma: ', current_gamma_ed

            write(10,'(A,F6.2,A,F8.2,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,F6.3)') &
             '  Star mass: ', masses(i), ' M_sun, Time: ', time, ' Myr, Mass: ', current_mass, ' M_sun, Teff: ', &
             current_teff, ' K, Luminosity: ', current_lum, ' L_sun, Mdot: ', current_mdot, ' M_sun/yr, Vterm: ', &
             current_v_term/kms, ' km/s, Gamma: ', current_gamma_ed
            write(20+i,*) time, current_mass, current_teff, current_lum, current_mdot, current_v_term/kms, current_gamma_ed
        enddo
        write(*,*) "-------------------------------------------------------"
        write(10,*) "-------------------------------------------------------"
    enddo


end program stevia