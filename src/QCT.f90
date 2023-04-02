program QCT
    use constants, only: dp, autofs, autouma, autocm_1, autoA, sal_unit, xyz_unit, end_unit
    use settings, only: initial_settings, ndim, nA, massA, atnameA, XP, XPini, propagation_mode
    use hamiltonian, only: derivs, get_derivs
    use ddeabm_module, wp => ddeabm_rk
    implicit none

    type(ddeabm_with_event_class) :: s
    character(len=80) :: initcond_file, traj_file
    integer :: ntrajs, itraj, maxcond, totalsteps, idid
    real(dp) :: tottime, t, timein, timeout, kener, max_step_factor, &
                potener, print_time, tprev, rfin, gval, Eini, Eend
    integer, parameter :: cond_unit = 11
    logical :: open_unit

    !variables for ddeabm
    real(8) :: relerr,abserr !absolute and relative error in propagation

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
        propagation_mode


    open(sal_unit, file="sal", status="replace")
    open(end_unit, file="end_conditions", status="replace")
    gval = 0._dp
    relerr = 1.e-8_dp
    abserr = 1.e-8_dp
    print_time = 0._dp
    max_step_factor = 10._dp
    propagation_mode = 0

    open(10,file="input.dat", status="old")
    read(10, nml=input)
    write(sal_unit, nml=input)
    call initial_settings()
    close(10)

    ! Set X and P derivatives routine
    call get_derivs(propagation_mode)

    ! Convert times
    totalsteps = int(max_step_factor * tottime)
    tottime=tottime/autofs
    print_time=print_time/autofs
    call s%initialize_event(ndim, maxnum=totalsteps, df=derivs, rtol=[relerr], atol=[abserr], &
        report=xyz_report, g=checkend, root_tol=1e-10_dp)


    open(cond_unit, file=trim(initcond_file), status="old")
    read(cond_unit,*) maxcond
    rewind(cond_unit)
    write(sal_unit,*) "Total number of initial conditions =", maxcond

    do itraj=1, ntrajs
        write(sal_unit,"(/A)") '---------------'
        tprev = 0._dp
        if (print_time > 0._dp) then
            write(traj_file,'("traj_",i6.6,".xyz")')itraj
            open(xyz_unit, file=trim(traj_file), status="replace")
        end if
        call s%first_call()
        t = 0._dp

        write(sal_unit,*) "Starting traj =", itraj
        call get_init_cond(ndim, XPini, maxcond, cond_unit)
        XP = XPini

        timein = 0._dp
        timeout = tottime
        call total_ener(XP, kener, potener)
        Eini = (kener + potener) * autocm_1
        call s%integrate_to_event(timein, XP, timeout, idid=idid, integration_mode=2, gval=gval)
        if (idid .eq. -1) then
            write(sal_unit,"(/A)") "*******************************************"
            write(sal_unit,"(A)") "ERROR: Maximum number of integration steps reached."
            write(sal_unit,"(A)") "Maximum number of iterations computed as int(max_step_factor * tottime)"
            write(sal_unit,"(A)") "Increase either parameter in the input to increase the maximum number of steps"
            write(sal_unit,"(A/)") "*******************************************"
        end if
        call total_ener(XP, kener, potener)
        Eend = (kener + potener) * autocm_1
        write(sal_unit,*) "Ener(initial) / cm-1:", Eini
        write(sal_unit,*) "Ener(final) / cm-1  :", Eend
        write(sal_unit,*) "Ener(delta) / cm-1  :", Eend - Eini
        write(sal_unit,*) "Final time / fs:", timein * autofs
        write(sal_unit,*) "End of traj =", itraj
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
            call total_ener(XP, kener, potener)
            tprev = t
            write(xyz_unit,*) nA
            write(xyz_unit,*) "t/fs, kinetic (au), pot (au)=", t * autofs, kener, potener
            do iat=1,nA
                write(xyz_unit,*) atnameA(iat), XP(3*(iat-1)+1:3*iat) * autoA, XP(ndim/2 + 3*(iat-1)+1:ndim/2 + 3*iat)
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

subroutine get_init_cond(ndim, XP, maxcond, cond_unit)
    use constants, only: dp, sal_unit
    implicit none
    integer, intent(in) :: ndim, maxcond, cond_unit
    real(dp), intent(inout) :: XP(ndim)
    real(dp) :: r
    integer :: icond, i 

    XP(:) = 0._dp

    call RANDOM_NUMBER(r)
    icond = floor(maxcond * r + 1)
    write(sal_unit,*) "Using icond =", icond

    read(cond_unit, *)
    do i=1, icond
        read(cond_unit, *) XP
    end do
    rewind(cond_unit)
end subroutine


subroutine total_ener(XP, k, pot)
    use constants, only: dp
    use settings, only: ndim
    implicit none
    real(dp), intent(in) :: XP(ndim)
    real(dp), intent(out) :: k, pot
    real(dp) :: der(ndim/2)

    call kinetic_ener(XP(ndim/2+1:), k)
    call potxyz(XP(:ndim/2), pot, der)
end subroutine

subroutine kinetic_ener(P, E)
    use constants, only: dp, sal_unit
    use settings, only: ndim, nA, massA
    implicit none
    integer :: iat, ix
    real(dp), intent(in) :: P(ndim/2)
    real(dp), intent(out) :: E

    E = 0._dp

    do iat=1,nA
        do ix=1,3
            E = E + P(3*(iat-1)+ix)**2 / massA(iat)
        end do
    end do
    E = E / 2._dp
end subroutine


