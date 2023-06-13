program QCT
    use constants, only: dp, autofs, autouma, autocm_1, autoA, sal_unit, xyz_unit, end_unit
    use settings, only: initial_settings, ndim, nA, massA, atnameA, XP, XPini, propagation_mode, &
        initcond_mode, Ts, temperature
    use hamiltonian, only: derivs, get_potential, total_ener
    use initial_conditions, only: set_init_cond, get_init_cond, initcond_file
    use physics, only: get_COM, get_LMOM_AMOM
    use ddeabm_module, wp => ddeabm_rk
    implicit none

    type(ddeabm_with_event_class) :: s
    character(len=80) :: traj_file
    integer :: ntrajs, itraj, maxcond, totalsteps, idid, time_init, time_end, time_rate
    real(dp) :: tottime, t, timein, timeout, kener, max_step_factor, &
                potener, print_time, tprev, rfin, gval, Eini, Eend, &
                relerr,abserr, &
                QCOM(3), PCOM(3), LMOM(3), AMOM(3)
    logical :: open_unit

    namelist /input/ &
        ntrajs, &
        nA, &
        tottime, &
        print_time, &
        max_step_factor, &
        relerr, &
        abserr, &
        rfin, &
        initcond_file, &
        propagation_mode, &
        initcond_mode, &
        Ts, &
        temperature

    ! Defaults
    gval = 0._dp
    relerr = 1.e-8_dp
    abserr = 1.e-8_dp
    print_time = 0._dp
    max_step_factor = 10._dp
    propagation_mode = 0
    initcond_mode = 0
    Ts = 0._dp
    initcond_file = ""

    open(sal_unit, file="sal", status="replace")
    open(end_unit, file="end_conditions", status="replace")
    open(10,file="input.dat", status="old")
    read(10, nml=input)
    write(sal_unit, nml=input)
    call initial_settings()
    call set_init_cond(initcond_mode) ! set initial conditions.
    call get_potential(propagation_mode)
    close(10)

    ! Convert times
    totalsteps = int(max_step_factor * tottime)
    tottime=tottime/autofs
    print_time=print_time/autofs
    call s%initialize_event(ndim, maxnum=totalsteps, df=derivs, rtol=[relerr], atol=[abserr], &
        report=xyz_report, g=checkend, root_tol=1e-10_dp)

    do itraj=1, ntrajs
        flush(sal_unit)
        write(sal_unit,"(/A)") '---------------'
        tprev = 0._dp
        if (print_time > 0._dp) then
            write(traj_file,'("traj_",i6.6,".xyz")')itraj
            open(xyz_unit, file=trim(traj_file), status="replace")
        end if
        call s%first_call()
        t = 0._dp

        write(sal_unit,*) "Starting traj =", itraj
        call get_init_cond(XPini)
        XP = XPini

        timein = 0._dp
        timeout = tottime
        call total_ener(timein, XP, kener, potener)
        call get_COM(ndim, XP, 1, nA, massA, QCOM, PCOM)
        call get_LMOM_AMOM(ndim, XP, 1, nA, massA, QCOM, PCOM, LMOM, AMOM)
        Eini = (kener + potener) * autocm_1
        write(sal_unit,*) "Ener(initial) / cm-1:", Eini
        write(sal_unit,*) "              COM / au:", QCOM
        write(sal_unit,*) "  Linear momentum / au:", PCOM
        write(sal_unit,*) " Angular momentum / au:", AMOM
        call system_clock(time_init, time_rate)
        call s%integrate_to_event(timein, XP, timeout, idid=idid, integration_mode=2, gval=gval)
        call system_clock(time_end)
        if (idid .eq. -1) then
            write(sal_unit,"(/A)") "*******************************************"
            write(sal_unit,"(A)") "ERROR: Maximum number of integration steps reached."
            write(sal_unit,"(A)") "Maximum number of iterations computed as int(max_step_factor * tottime)"
            write(sal_unit,"(A)") "Increase either parameter in the input to increase the maximum number of steps"
            write(sal_unit,"(A/)") "*******************************************"
        end if
        call total_ener(timeout, XP, kener, potener)
        call get_COM(ndim, XP, 1, nA, massA, QCOM, PCOM)
        call get_LMOM_AMOM(ndim, XP, 1, nA, massA, QCOM, PCOM, LMOM, AMOM)
        Eend = (kener + potener) * autocm_1
        write(sal_unit,*) "Ener(final) / cm-1  :", Eend
        write(sal_unit,*) "              COM / au:", QCOM
        write(sal_unit,*) "  Linear momentum / au:", LMOM
        write(sal_unit,*) " Angular momentum / au:", AMOM
        write(sal_unit,*) "Ener(delta) / cm-1  :", Eend - Eini
        write(sal_unit,*) "Final time / fs:", timein * autofs
        write(sal_unit,*) "End of traj =", itraj
        write(sal_unit,*) "Time elapsed / s:", real(time_end - time_init) / real(time_rate)
        write(end_unit,*) itraj, XPini, XP
        inquire(unit=xyz_unit, opened=open_unit)
        if (open_unit) close(xyz_unit)
        write(sal_unit,"(A/)") '---------------'
    end do
    close(sal_unit)
    close(end_unit)

    contains 

    subroutine  xyz_report(me, t, XP)
        implicit none
        class(ddeabm_class), intent(inout) :: me
        integer :: iat
        real(dp), intent(in) :: t, XP(:)
        real(dp) :: kener, potener

        if (t - tprev > print_time .and. print_time > 0._dp) then
            call total_ener(t, XP, kener, potener)
            tprev = t
            write(xyz_unit,*) nA
            write(xyz_unit,*) t * autofs, kener, potener, " = t/fs, kinetic (au), pot (au)"
            do iat=1,nA
                write(xyz_unit,*) atnameA(iat), XP(3*(iat-1)+1:3*iat) * autoA
            end do
        end if

    end subroutine

    subroutine checkend(me, t, XP, is_end)
        implicit none
        class(ddeabm_with_event_class), intent(inout) :: me
        integer :: iat1, iat2, ix
        real(dp), intent(in) :: t, XP(:)
        real(dp), intent(out) :: is_end
        real(dp) :: d

        is_end = 1._dp
        outer: do iat1=1,nA
            do iat2=iat1+1,nA
                d = 0._dp
                do ix=1,3
                    d = d + (XP(3*(iat1-1)+ix) - XP(3*(iat2-1)+ix))**2
                end do
                d = sqrt(d)
                if (d > rfin) then
                    is_end = 0._dp
                    exit outer
                end if
            end do
        end do outer
    end subroutine
end program
