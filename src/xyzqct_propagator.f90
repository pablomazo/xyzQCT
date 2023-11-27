module xyzqct_propagator
    use xyzqct_constants, only: dp, sal_unit, autofs
    use ddeabm_module, wp => ddeabm_rk
    use xyzqct_settings, only: ndim
    type(ddeabm_with_event_class) :: s
    real(dp) :: gval, relerr, abserr, tottime, tprev, rfin, print_time, deltat, deltat2, tprev_as, &
        print_time_as
    integer :: max_step_factor, iprop

    public :: xyz_report, checkend
    private:: iprop
    contains
        subroutine set_propagator(ipropagator)
            ! Select propagator:
            implicit none
            integer, intent(in) :: ipropagator
            iprop = ipropagator
            select case(iprop)
                case(0)
                    write(sal_unit,"(/A/)") "Using DDEABM propagator."
                    call set_DDEABM_propagator()
                case(1)
                    write(sal_unit,"(/A/)") "Using Verlet propagator."
                    call set_verlet_propagator()
                case default
                    write(sal_unit,"(/A/)") "Unknown propagator."
                    stop
            end select
        end subroutine

        subroutine propagate(XP, ti, tf, ptime, ptime_as, r, final_t, elapsed)
            implicit none
            real(dp), intent(in) :: ti, tf, ptime, ptime_as, r
            real(dp), intent(inout) :: XP(ndim)
            real(dp), intent(out) :: final_t, elapsed
            print_time = ptime
            print_time_as = ptime_as
            rfin = r
            tprev = 0._dp
            select case(iprop)
                case(0)
                    call propagate_DDEABM(XP, ti, tf, final_t, elapsed)
                case(1)
                    call propagate_verlet(XP, ti, tf, final_t, elapsed)
            end select
        end subroutine

        subroutine set_DDEABM_propagator()
            use xyzqct_hamiltonian, only: derivs
            implicit none
            integer :: ios, totalsteps
            namelist /propagator/ relerr, abserr, max_step_factor
            gval = 0._dp
            relerr = 1.e-8_dp
            abserr = 1.e-8_dp
            max_step_factor = 10._dp
            rewind(10)
            read(10, nml=propagator, iostat=ios)
            write(sal_unit, nml=propagator)
        end subroutine

        subroutine propagate_DDEABM(XP, ti, tf, final_t, elapsed)
            use xyzqct_hamiltonian, only: derivs
            implicit none
            real(dp), intent(in) :: ti, tf
            real(dp), intent(inout) :: XP(ndim)
            real(dp), intent(out) :: elapsed, final_t
            real(dp) :: timein, timeout
            integer :: time_init, time_rate, time_end, idid, totalsteps

            totalsteps = int(max_step_factor * tf)
            timein = ti
            timeout = tf
            call s%first_call()
            call s%initialize_event(ndim, maxnum=totalsteps, df=derivs, rtol=[relerr], atol=[abserr], &
                report=xyz_report, g=checkend, root_tol=1e-10_dp)
            call system_clock(time_init, time_rate)
            call s%integrate_to_event(timein, XP, timeout, idid=idid, integration_mode=2, gval=gval)
            call system_clock(time_end)
            if (idid < 0) then
                write(sal_unit,"(/A)") "*******************************************"
                write(sal_unit,"(A)") "ERROR during DDEABM propagation. Routine exited with idid:"
                write(sal_unit,*) idid
                if (idid .eq. -1) then
                    write(sal_unit,"(/A)") "ERROR: Maximum number of integration steps reached."
                    write(sal_unit,"(A)") "Maximum number of iterations computed as int(max_step_factor * tottime)"
                    write(sal_unit,"(A)") "Increase either parameter in the input to increase the maximum number of steps"
                end if
            write(sal_unit,"(A/)") "*******************************************"
            end if
            elapsed = real(time_end - time_init) / real(time_rate)
            final_t = timein
        end subroutine

        subroutine set_verlet_propagator()
            use xyzqct_constants, only: autofs
            implicit none
            integer :: ios
            namelist/propagator/ deltat
            deltat = 1._dp
            rewind(10)
            read(10, nml=propagator, iostat=ios)
            write(sal_unit, nml=propagator)
            deltat = deltat / autofs
            deltat2 = deltat / 2._dp
        end subroutine

        subroutine propagate_verlet(XP, ti, tf, final_t, elapsed)
            use xyzqct_hamiltonian, only: derivs
            use xyzqct_settings, only: nat, mass
            implicit none
            real(dp), intent(in) :: ti, tf
            real(dp), intent(inout) :: XP(ndim)
            real(dp), intent(out) :: elapsed, final_t
            real(dp) :: time, XPder(ndim), is_end, mass_(ndim/2)
            integer :: numsteps, istep, iat
            integer :: time_init, time_rate, time_end

            numsteps = int((tf - ti) / deltat) + 1
            do iat=1,nat
                mass_(3*(iat-1)+1:3*iat) = mass(iat)
            end do

            call system_clock(time_init, time_rate)
            time = ti
            call derivs(s, time, XP, XPder)
            do istep=1, numsteps
                time = (istep - 1) * deltat
                XP(ndim/2+1:) = XP(ndim/2+1:) + deltat2 * XPder(ndim/2+1:)
                XP(:ndim/2) = XP(:ndim/2) + deltat * XP(ndim/2+1:) / mass_
                call derivs(s, time, XP, XPder)
                XP(ndim/2+1:) = XP(ndim/2+1:) + deltat2 * XPder(ndim/2+1:)
                call xyz_report(s, time, XP)
                call checkend(s, time, XP, is_end)
                if (is_end .eq. 0._dp) exit
            end do
            call system_clock(time_end)
            elapsed = real(time_end - time_init) / real(time_rate)
            final_t = time
        end subroutine

        subroutine xyz_report(me, t, XP)
            use xyzqct_constants, only: xyz_unit, autoA, as_unit
            use xyzqct_settings, only: nat, atname, Ts
            use xyzqct_hamiltonian, only: total_ener
            implicit none
            class(ddeabm_class), intent(inout) :: me
            integer :: iat
            real(dp), intent(in) :: t, XP(:)
            real(dp) :: kener, potener

            if (t - tprev > print_time .and. print_time > 0._dp) then
                call total_ener(t, XP, kener, potener)
                tprev = t
                write(xyz_unit,*) nat
                write(xyz_unit,*) t * autofs, kener, potener, " = t/fs, kinetic (au), pot (au)"
                do iat=1,nat
                    write(xyz_unit,*) atname(iat), XP(3*(iat-1)+1:3*iat) * autoA
                end do
            end if

            if (Ts > 0 .and. t > Ts .and. t - tprev_as > print_time_as) then
                tprev_as = t
                write(as_unit,*) XP
            end if
        end subroutine

        subroutine checkend(me, t, XP, is_end)
            use xyzqct_settings, only: nat
            implicit none
            class(ddeabm_with_event_class), intent(inout) :: me
            integer :: iat1, iat2, ix
            real(dp), intent(in) :: t, XP(:)
            real(dp), intent(out) :: is_end
            real(dp) :: d

            is_end = 1._dp
            outer: do iat1=1,nat
                do iat2=iat1+1,nat
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
end module
