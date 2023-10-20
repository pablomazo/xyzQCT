program QCT
    use constants, only: dp, autofs, autouma, autocm_1, autoA, &
        sal_unit, xyz_unit, end_unit, as_unit
    use settings, only: initial_settings, ndim, nA, mass, XP, XPini, potential_mode, &
        initcond_mode, Ts, temperature, propagator_mode, nB, nat
    use hamiltonian, only: derivs, get_potential, total_ener
    use initial_conditions, only: set_init_cond, get_init_cond
    use physics, only: get_COM, get_LMOM_AMOM
    use propagator, only: set_propagator, propagate, reset_propagator
    use ddeabm_module, wp => ddeabm_rk
    use utils, only: code_starter
    implicit none

    character(len=80) :: traj_file
    integer :: ntrajs, itraj, nini
    real(dp) :: tottime, timein, timeout, kener, &
                potener, print_time, tprev, rfin, Eini, Eend, &
                init_cond_print, final_t, &
                QCOM(3), PCOM(3), LMOM(3), AMOM(3), elapsed
    logical :: open_unit

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
        Ts, &
        temperature, &
        init_cond_print

    ! Defaults
    nini = 1
    nB = 0
    potential_mode = 0
    initcond_mode = 0
    propagator_mode = 0
    Ts = 0._dp
    init_cond_print = 0._dp
    rfin=200._dp

    open(sal_unit, file="sal", status="replace")
    open(end_unit, file="end_conditions", status="replace")
    open(10,file="input.dat", status="old")
    read(10, nml=input)
    call code_starter()
    write(sal_unit, nml=input)
    ! Convert times
    Ts = Ts/autofs
    tottime=tottime/autofs
    print_time=print_time/autofs
    call initial_settings()
    call set_init_cond(initcond_mode) ! set initial conditions.
    call get_potential(potential_mode)
    call set_propagator(propagator_mode, tottime, print_time, rfin)
    close(10)


    do itraj=nini, ntrajs
        flush(sal_unit)
        write(sal_unit,"(/A)") '---------------'
        write(sal_unit,*) "Starting traj =", itraj
        tprev = 0._dp
        if (print_time > 0._dp) then
            write(traj_file,'("traj_",i6.6,".xyz")')itraj
            open(xyz_unit, file=trim(traj_file), status="replace")
        end if
        call reset_propagator()
        call get_init_cond(XPini)
        XP = XPini

        call total_ener(timein, XP, kener, potener)
        call get_COM(ndim, XP, 1, nat, mass, QCOM, PCOM)
        call get_LMOM_AMOM(ndim, XP, 1, nat, mass, QCOM, PCOM, LMOM, AMOM)
        Eini = (kener + potener) * autocm_1
        write(sal_unit,*) "Ener(initial) / cm-1:", Eini
        write(sal_unit,*) "              COM / au:", QCOM
        write(sal_unit,*) "  Linear momentum / au:", PCOM
        write(sal_unit,*) " Angular momentum / au:", AMOM
        call flush(sal_unit)
        call propagate(XP, final_t, elapsed)
        call total_ener(timeout, XP, kener, potener)
        call get_COM(ndim, XP, 1, nat, mass, QCOM, PCOM)
        call get_LMOM_AMOM(ndim, XP, 1, nat, mass, QCOM, PCOM, LMOM, AMOM)
        Eend = (kener + potener) * autocm_1
        write(sal_unit,*) "Ener(final) / cm-1  :", Eend
        write(sal_unit,*) "              COM / au:", QCOM
        write(sal_unit,*) "  Linear momentum / au:", LMOM
        write(sal_unit,*) " Angular momentum / au:", AMOM
        write(sal_unit,*) "Ener(delta) / cm-1  :", Eend - Eini
        write(sal_unit,*) "Final time / fs:", final_t * autofs
        write(sal_unit,*) "End of traj =", itraj
        write(sal_unit,*) "Time elapsed / s:", elapsed
        write(end_unit,*) itraj, XPini, XP
        inquire(unit=xyz_unit, opened=open_unit)
        if (open_unit) close(xyz_unit)
        inquire(unit=as_unit, opened=open_unit)
        if (open_unit) close(as_unit)
        write(sal_unit,"(A/)") '---------------'
    end do
    close(sal_unit)
    close(end_unit)
end program
