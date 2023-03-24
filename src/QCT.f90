program QCT
    use constants, only: dp, autofs, autouma, autocm_1, autoA, sal_unit, xyz_unit, end_unit
    use settings, only: ndim, nA, massA, atnameA
    use ddeabm_module, wp => ddeabm_rk
    implicit none

    type(ddeabm_class) :: s
    character(len=80) :: initcond_file, traj_file
    integer :: ntrajs, itraj, maxcond, totalsteps
    real(dp) :: tottime, tstep, t, timein, timeout, ener, print_time, tprev
    real(dp), allocatable :: XP(:)
    integer, parameter :: cond_unit = 11

    !variables for ddeabm
    integer :: idid,lrw,liw
    real(8), allocatable :: rwork(:)
    real(8) :: relerr,abserr   !absolute and relative error in propagation

    namelist /input/ &
        ntrajs, &
        nA, &
        tottime, &
        tstep, &
        print_time, &
        relerr, &
        abserr, &
        initcond_file

    namelist /mass/ &
        massA, &
        atnameA


    open(sal_unit, file="sal", status="replace")
    open(end_unit, file="end_conditions", status="replace")
    relerr = 1.e-8_dp
    abserr = 1.e-8_dp
    print_time = 0._dp

    open(10,file="input.dat", status="old")
    read(10, nml=input)
    ndim = 2 * 3 * nA
    allocate(XP(ndim), rwork(lrw), massA(nA), atnameA(nA))
    read(10, nml=mass)
    close(10)

    !variables for ddeab subroutine
    lrw=130+21*ndim
    liw=51

    ! Convert times
    tottime=tottime/autofs
    tstep=tstep/autofs
    print_time=print_time/autofs
    totalsteps = int(tottime / tstep)
    call s%initialize(ndim, maxnum=totalsteps, df=derivs, rtol=[relerr], atol=[abserr], &
        report=xyz_report)
    !call s%initialize(ndim, maxnum=50, df=derivs, rtol=[relerr], atol=[abserr])

    ! Convert mass
    massA = massA/autouma

    open(cond_unit, file=trim(initcond_file), status="old")
    read(cond_unit,*) maxcond
    rewind(cond_unit)
    write(sal_unit,*) "Total number of initial conditions =", maxcond

    do itraj=1, ntrajs
        tprev = 0._dp
        write(traj_file,'("traj_",i6.6,".xyz")')itraj
        open(xyz_unit, file=trim(traj_file), status="replace")
        call s%first_call()
        t = 0._dp

        write(sal_unit,*) "Starting traj =", itraj
        call get_init_cond(ndim, XP, maxcond, cond_unit)

        timein = 0._dp
        timeout = tottime
        idid = 0
        call total_ener(XP, ener)
        write(sal_unit,*) "Ener init/cm-1 =", ener * autocm_1

        call s%integrate(timein, XP, timeout, idid=idid, integration_mode=2)

        write(sal_unit,*) "Final time / fs:", timein * autofs
        call total_ener(XP, ener)
        write(sal_unit,*) "Ener end/cm-1 =", ener * autocm_1
        write(sal_unit,*) "End of traj =", itraj
        write(end_unit,*) itraj, XP
        close(xyz_unit)
    end do
    close(sal_unit)
    close(end_unit)

    contains 

    subroutine derivs(me, t, XP, XPder)
        implicit none

        class(ddeabm_class), intent(inout) :: me
        integer ::iat, ix
        real(dp), intent(in) :: t, XP(:)
        real(dp), intent(out) :: XPder(:)
        real(dp) :: pot, posxyz(ndim/2), P(ndim/2), derxyz(ndim/2)

        XPder = 0._dp
        pot = 0._dp
        derxyz = 0._dp
        posxyz = 0._dp

        posxyz = XP(:ndim/2)
        P = XP(ndim/2+1:)
        call potxyz(posxyz, pot, derxyz)
        do iat=1,nA
            do ix=1,3
                XPder(3 * (iat-1) + ix) = P(3 * (iat-1) + ix) / massA(iat)
                XPder(ndim/2 + 3 * (iat-1) + ix) = -derxyz(3 * (iat-1) + ix)
            end do
        end do
    end subroutine derivs

    subroutine  xyz_report(me, t, XP)
        implicit none
        class(ddeabm_class), intent(inout) :: me
        integer :: iat
        real(dp), intent(in) :: t, XP(:)

        if (t - tprev > print_time) then
            tprev = t
            write(xyz_unit,*) nA
            write(xyz_unit,*) "t=", t * autofs
            do iat=1,nA
                write(xyz_unit,*) atnameA(iat), XP(3*(iat-1)+1:3*iat) * autoA
            end do
        end if

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
    icond = 1
    write(sal_unit,*) "Using icond =", icond
    write(sal_unit,*) "Initial condition in unit =", cond_unit

    read(cond_unit, *)
    do i=1, icond
        read(cond_unit, *) XP
    end do
    rewind(cond_unit)
end subroutine


subroutine total_ener(XP, ener)
    use constants, only: dp, sal_unit
    use settings, only: ndim
    implicit none
    real(dp), intent(in) :: XP(ndim)
    real(dp), intent(out) :: ener
    real(dp) :: der(ndim/2), pot, k

    call kinetic_ener(XP(ndim/2+1:), k)
    call potxyz(XP(:ndim/2), pot, der)
    ener = k + pot
    write(sal_unit,*) "Kinetic =", k
    write(sal_unit,*) "Potential =", pot
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
    write(sal_unit,*) '---------------'
    E = E / 2._dp
end subroutine


