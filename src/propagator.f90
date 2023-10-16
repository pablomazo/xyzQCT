module propagator
    use constants, only: dp, sal_unit, autofs
    use ddeabm_module, wp => ddeabm_rk
    use settings, only: ndim
    type(ddeabm_with_event_class) :: s
    real(dp) :: gval, relerr, abserr, tottime, tprev, rfin, print_time
    integer :: max_step_factor, iprop

    public :: xyz_report, checkend
    private:: iprop
    contains
        subroutine set_propagator(ipropagator, ttime, ptime, r)
            ! Select propagator:
            implicit none
            integer, intent(in) :: ipropagator
            real(dp), intent(in) :: ttime, ptime, r
            iprop = ipropagator
            tottime = ttime
            print_time = ptime
            rfin = r
            select case(iprop)
                case(0)
                    write(sal_unit,"(/A/)") "Using DDEABM propagator."
                    call set_DDEABM_propagator()
            end select
        end subroutine

        subroutine propagate(XP, final_t, elapsed)
            implicit none
            real(dp), intent(inout) :: XP(ndim)
            real(dp), intent(out) :: final_t, elapsed
            select case(iprop)
                case(0)
                    call propagate_DDEABM(XP, final_t, elapsed)
            end select
        end subroutine

        subroutine reset_propagator()
            implicit none
            tprev = 0._dp
            select case(iprop)
                case(0)
                    call s%first_call()
            end select
        end subroutine

        subroutine set_DDEABM_propagator()
            use hamiltonian, only: derivs
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
            totalsteps = int(max_step_factor * tottime)
            call s%initialize_event(ndim, maxnum=totalsteps, df=derivs, rtol=[relerr], atol=[abserr], &
                report=xyz_report, g=checkend, root_tol=1e-10_dp)
        end subroutine

        subroutine propagate_DDEABM(XP, final_t, elapsed)
            implicit none
            real(dp), intent(inout) :: XP(ndim)
            real(dp), intent(out) :: elapsed, final_t
            real(dp) :: timein, timeout
            integer :: time_init, time_rate, time_end, idid

            timein = 0._dp
            timeout = tottime
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
            elapsed = real(time_end - time_init) / real(time_rate)
            final_t = timein
        end subroutine

        subroutine xyz_report(me, t, XP)
            use constants, only: xyz_unit, autoA
            use settings, only: potential_mode, nA, atnameA
            use hamiltonian, only: total_ener
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
            use settings, only: nA
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
end module