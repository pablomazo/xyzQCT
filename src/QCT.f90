program QCT
    use constants, only: dp, autofs, autouma, autocm_1
    use settings, only: ndim, nA, massA
    use ddeabm_module, wp => ddeabm_rk
    implicit none

    type(ddeabm_class) :: s
    character(len=80) :: initcond_file
    integer :: ntrajs, itraj, maxcond, totalsteps
    real(dp) :: tottime, tstep, t, timein, timeout, ener
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
        relerr, &
        abserr, &
        initcond_file

    namelist /mass/ &
        massA


    relerr = 1.e-8_dp
    abserr = 1.e-8_dp

    open(10,file="input.dat", status="old")
    read(10, nml=input)
    ndim = 2 * 3 * nA
    allocate(XP(ndim), rwork(lrw), massA(nA))
    read(10, nml=mass)
    close(10)

    !variables for ddeab subroutine
    lrw=130+21*ndim
    liw=51

    ! Convert times
    tottime=tottime/autofs
    tstep=tstep/autofs
    totalsteps = int(tottime / tstep)
    call s%initialize(ndim, maxnum=totalsteps, df=derivs, rtol=[relerr], atol=[abserr])
    !call s%initialize(ndim, maxnum=50, df=derivs, rtol=[relerr], atol=[abserr])

    ! Convert mass
    massA = massA/autouma

    open(cond_unit, file=trim(initcond_file), status="old")
    read(cond_unit,*) maxcond
    rewind(cond_unit)
    write(*,*) "Total number of initial conditions =", maxcond

    do itraj=1, ntrajs
        call s%first_call()
        t = 0._dp

        write(*,*) "Starting traj =", itraj
        call get_init_cond(ndim, XP, maxcond, cond_unit)

        timein = 0._dp
        timeout = tottime
        idid = 0
        call total_ener(XP, ener)
        write(*,*) "Ener/cm-1 =", ener * autocm_1

        call s%integrate(timein, XP, timeout, idid=idid)
        write(*,*) "Final time / fs:", timein * autofs
        call total_ener(XP, ener)
        write(*,*) XP(:ndim/2)
        write(*,*) "Ener/cm-1 =", ener * autocm_1
        write(*,*) "End of traj =", itraj
    end do

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
end program

subroutine get_init_cond(ndim, XP, maxcond, cond_unit)
    use constants, only: dp
    implicit none
    integer, intent(in) :: ndim, maxcond, cond_unit
    real(dp), intent(inout) :: XP(ndim)
    real(dp) :: r
    integer :: icond, i 

    XP(:) = 0._dp

    call RANDOM_NUMBER(r)
    icond = floor(maxcond * r + 1)
    icond = 1
    write(*,*) "Using icond =", icond
    write(*,*) "Initial condition in unit =", cond_unit

    read(cond_unit, *)
    do i=1, icond
        read(cond_unit, *) XP
    end do
    rewind(cond_unit)
end subroutine


subroutine total_ener(XP, ener)
    use constants, only: dp
    use settings, only: ndim
    implicit none
    real(dp), intent(in) :: XP(ndim)
    real(dp), intent(out) :: ener
    real(dp) :: der(ndim/2), pot, k

    call kinetic_ener(XP(ndim/2+1:), k)
    call potxyz(XP(:ndim/2), pot, der)
    ener = k + pot
    write(*,*) "Kinetic =", k
    write(*,*) "Potential =", pot
end subroutine

subroutine kinetic_ener(P, E)
    use constants, only: dp
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
    write(*,*) '---------------'
    E = E / 2._dp
end subroutine


