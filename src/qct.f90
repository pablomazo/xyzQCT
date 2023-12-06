program QCT
    use xyzqct_constants, only: dp, autofs, autouma, autocm_1, autoA, &
        sal_unit, xyz_unit, end_unit, as_unit
    use xyzqct_settings, only: initial_settings, ndim, nA, mass, XP, XPini, potential_mode, &
        initcond_mode, temperature, propagator_mode, nB, nat, rfin, Ts, sysA, sysB
    use xyzqct_hamiltonian, only: derivs, get_potential, total_ener
    use xyzqct_initial_conditions, only: set_init_cond, get_init_cond, write_end_cond
    use xyzqct_physics, only: get_COM, get_LMOM_AMOM, get_angular_velocity, add_angular_velocity, &
        get_inertia_moments, matrix_rotation
    use xyzqct_propagator, only: set_propagator, propagate
    use ddeabm_module, wp => ddeabm_rk
    use xyzqct_utils, only: code_starter
    implicit none

    character(len=80) :: traj_file
    integer :: ntrajs, itraj, nini, seed_size, iat
    real(dp) :: tottime, kener, &
                potener, print_time, tprev, Eini, Eend, E_s,&
                init_cond_print, final_t, &
                QCOM(3), PCOM(3), LMOM(3), AMOM(3), elapsed, &
                inertia(3), inertia_vec(3,3), omega(3), dum
    logical :: open_unit, prop
    integer, allocatable :: seed(:)

    namelist /input/ &
        nini, &
        ntrajs, &
        nA, &
        nB, &
        print_time, &
        tottime,&
        rfin, &
        potential_mode, &
        initcond_mode, &
        propagator_mode, &
        temperature, &
        init_cond_print

    ! Defaults
    nini = 1
    nB = 0
    potential_mode = 0
    initcond_mode = 0
    propagator_mode = 0
    init_cond_print = 0._dp
    print_time = 0._dp
    rfin=200._dp

    open(sal_unit, file="sal", status="replace")
    open(end_unit, file="end_conditions", status="replace")
    open(10,file="input.dat", status="old")
    read(10, nml=input)
    call code_starter()
    write(sal_unit, nml=input)
    ! Convert times
    tottime=tottime/autofs
    print_time=print_time/autofs
    init_cond_print=init_cond_print/autofs

    call initial_settings()
    call set_init_cond(initcond_mode) ! set initial conditions.
    call get_potential(potential_mode)
    call set_propagator(propagator_mode)
    close(10)

    call RANDOM_SEED(size=seed_size)
    allocate(seed(seed_size))
    do itraj=nini, ntrajs
        final_t = 0._dp
        elapsed = 0._dp
        seed = itraj
        call RANDOM_SEED(PUT=seed) ! To guarantee that trajectories are reproducible.
                                   ! There are probably better ways...
        call RANDOM_NUMBER(dum) ! Improves the first few random numbers
        flush(sal_unit)
        write(sal_unit,"(/A)") '---------------'
        write(sal_unit,*) "Starting traj =", itraj
        tprev = 0._dp
        if (print_time > 0._dp) then
            write(traj_file,'("traj_",i6.6,".xyz")')itraj
            open(xyz_unit, file=trim(traj_file), status="replace")
        end if
        call get_init_cond(XPini, prop)
        XP = XPini

        call total_ener(0._dp, XP, kener, potener)
        call get_COM(ndim, XP, 1, nat, mass, QCOM, PCOM)
        call get_LMOM_AMOM(ndim, XP, 1, nat, mass, QCOM, PCOM, LMOM, AMOM)
        Eini = (kener + potener) * autocm_1
        write(sal_unit,*) "Ener(initial) / cm-1:", Eini
        write(sal_unit,*) " COM / au:", QCOM
        write(sal_unit,*) " |P| / au:", sqrt(sum(PCOM)**2)
        write(sal_unit,*) " |L| / au:", sqrt(sum(AMOM**2))
        call flush(sal_unit)
        if (prop) then
            if (potential_mode .eq. 2) then
                ! Adiabatic switching propagation
                call propagate(XP, 0.0_dp, Ts, print_time, init_cond_print, rfin, final_t, elapsed)
                call total_ener(final_t, XP, kener, potener)
                E_s = (kener + potener) * autocm_1
                write(sal_unit,*) "Ener after switch / cm-1  :", E_s

                ! Remove spurious angular momenta
                call get_COM(ndim, XP,1, nat, mass, QCOM, PCOM)
                do iat=1,nat
                    XP(3*(iat-1)+1:3*iat) = XP(3*(iat-1)+1:3*iat) - QCOM
                end do
                call get_LMOM_AMOM(ndim, XP,1, nat, mass, QCOM, PCOM, LMOM, AMOM)
                write(sal_unit,*) "Total angular momentum after switch: ", sqrt(sum(AMOM**2))
                write(sal_unit,*) "Filtering angular momentum..."
                call get_inertia_moments(nat, XP(:3*nat), mass, inertia, inertia_vec)
                call matrix_rotation(3, nat, XP(:3*nat), inertia_vec)
                call matrix_rotation(3, nat, XP(3*nat+1:), inertia_vec)
                call get_LMOM_AMOM(ndim, XP,1, nat, mass, QCOM, PCOM, LMOM, AMOM)
                call get_angular_velocity(inertia, AMOM, omega)
                call add_angular_velocity(nat, XP, mass, omega, -1._dp)
                call total_ener(final_t, XP, kener, potener)
                Eend = (kener + potener) * autocm_1
                write(sal_unit,*) "Scaling momenta to match energy after switch..."
                do while (abs(E_s - Eend) / E_s > 1e-3)
                    XP(3*nat+1:) = XP(3*nat+1:) * (E_s / Eend)
                    call total_ener(final_t, XP, kener, potener)
                    Eend = (kener + potener) * autocm_1
                end do
                write(sal_unit,*) "Ener after angular momentum removal / cm-1  :", Eend
            end if
            call propagate(XP, Ts, tottime, print_time, init_cond_print, rfin, final_t, elapsed)
        end if
        call total_ener(final_t, XP, kener, potener)
        call get_COM(ndim, XP, 1, nat, mass, QCOM, PCOM)
        call get_LMOM_AMOM(ndim, XP, 1, nat, mass, QCOM, PCOM, LMOM, AMOM)
        Eend = (kener + potener) * autocm_1
        write(sal_unit,*) "Ener(final) / cm-1  :", Eend
        write(sal_unit,*) " COM / au:", QCOM
        write(sal_unit,*) " |P| / au:", sqrt(sum(PCOM)**2)
        write(sal_unit,*) " |L| / au:", sqrt(sum(AMOM**2))
        write(sal_unit,*) "Ener(delta) / cm-1  :", Eend - Eini
        write(sal_unit,*) "Final time / fs:", final_t * autofs
        write(sal_unit,*) "End of traj =", itraj
        write(sal_unit,*) "Time elapsed / s:", elapsed
        call write_end_cond(itraj, final_t * autofs, XPini, XP)
        inquire(unit=xyz_unit, opened=open_unit)
        if (open_unit) close(xyz_unit)
        inquire(unit=as_unit, opened=open_unit)
        if (open_unit) close(as_unit)
        write(sal_unit,"(A/)") '---------------'
    end do
    close(sal_unit)
    close(end_unit)
end program
