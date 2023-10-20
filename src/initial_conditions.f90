module initial_conditions
    use constants, only: dp, sal_unit, cond_unit
    use settings, only : ndim, Qnum, amp, nfreqs, freqs, CXQ, massA, Xeq, atnameA, nA, Qmax, temperature
    implicit none
    integer :: init_cond_mode, max_cond

    private :: init_cond_mode

    namelist /Qvib/ &
        Qnum

    ! init_cond_mode:
    ! 0 => Read from file
    ! 1 => NM initial condition
    ! 2 => NM initial condition (sampled at temperature T)

    contains
        subroutine set_init_cond(mode)
            use settings, only: initcond_fileA, initcond_fileB
            implicit none
            integer :: mode
            init_cond_mode = mode
            select case(mode)
                case(0)
                    write(sal_unit, "(/A)") "Reading initial conditions from file:", trim(initcond_fileA)
                    open(cond_unit, file=trim(initcond_fileA), status="old")
                    read(cond_unit,*) max_cond
                    rewind(cond_unit)
                    write(sal_unit,*) "Total number of initial conditions =", max_cond
                case(1)
                    write(sal_unit, "(/A)") "Using NM initial conditions"
                    open(11, file=trim(initcond_fileA), status="old")
                    call read_Data4NM(11)
                    close(11)
                case(2)
                    write(sal_unit, "(/A)") "Using NM initial conditions (sampled at temperature T)"
                    open(11, file=initcond_fileA, status="old")
                    call read_Data4NM(11)
                    close(11)
                    call compute_Qmax(temperature, Qmax)
            end select
        end subroutine

        subroutine get_init_cond(XP)
            implicit none
            real(dp), intent(out) :: XP(ndim)
            select case(init_cond_mode)
                case(0)
                    call from_file_init_cond(XP)
                case(1)
                    call NM_init_cond(XP)
                case(2)
                    call NM_init_cond_T(XP)
            end select 
        end subroutine

        subroutine from_file_init_cond(XP)
            implicit none
            real(dp), intent(out) :: XP(ndim)
            real(dp) :: r
            integer :: icond, i 

            XP(:) = 0._dp

            call RANDOM_NUMBER(r)
            icond = floor(max_cond * r + 1)
            write(sal_unit,*) "Using icond =", icond

            read(cond_unit, *)
            do i=1, icond
                read(cond_unit, *) XP
            end do
            rewind(cond_unit)
        end subroutine

        subroutine NM_init_cond(XP)
            use constants, only: pi
            implicit none
            integer :: ix, ifreq
            real(dp), intent(out) :: XP(ndim)
            real(dp) :: phase(nfreqs), Q(nfreqs), P(nfreqs), mass

            amp = sqrt((2._dp * Qnum +1._dp) / freqs) ! Maximum NM amplitudes

            call RANDOM_NUMBER(phase)
            XP = 0._dp
            phase = 2._dp * pi * phase
            Q = amp * sin(freqs * phase)
            P = freqs * amp * cos(freqs * phase)

            do ix=1, ndim/2
                mass = sqrt(massA((ix-1)/ 3 + 1))
                do ifreq=1, nfreqs
                    XP(ix) = XP(ix) + Q(ifreq) * CXQ(ix, ifreq)
                    XP(ix+ndim/2) = XP(ix+ndim/2) + P(ifreq) * CXQ(ix, ifreq)
                end do
                XP(ix) = XP(ix) / mass
                XP(ndim/2 + ix) = XP(ndim/2 + ix) * mass
            end do
            XP(:ndim/2) = XP(:ndim/2) + Xeq
        end subroutine

        subroutine NM_init_cond_T(XP)
            use constants, only: h, c, kb, autocm_1
            use settings, only: Qmax, temperature
            implicit none
            integer :: imode, ii
            real(dp), intent(out) :: XP(ndim)
            real(dp) :: rand1, rand2, factor, prob

            ! Sample the values of Qnum for temperature T
            do imode = 1, nfreqs
                factor = h *  freqs(imode) * autocm_1 * 100 * c / (kb * temperature)
                do ii=1,1000
                call random_number(rand1)
                call random_number(rand2)
                Qnum(imode) = floor((Qmax(imode) + 1) * rand1)
                prob = exp(-Qnum(imode) * factor) * (1. - exp(-factor))
                if (rand2 < prob) exit
                end do
            end do
            call NM_init_cond(XP)
        end subroutine

      subroutine read_Data4NM(read_unit)
          use constants, only: autoA, autocm_1
          implicit none
          integer, intent(in) :: read_unit
          integer :: iat, ix, ios

          read(read_unit, *)
          read(read_unit, *)
          do iat=1,nA
              read(read_unit, *) atnameA(iat), Xeq(3*(iat-1)+1:3*iat)
          end do
          Xeq = Xeq / autoA
          read(read_unit, *) nfreqs
          allocate(freqs(nfreqs), CXQ(3*nA, nfreqs), Qnum(nfreqs), amp(nfreqs), Qmax(nfreqs))
          read(read_unit, *) freqs
          do ix=1, 3*nA
              read(read_unit,*) CXQ(ix,:)
          end do
          freqs = freqs / autocm_1
          Qnum = 0 ! Vibrational quantum numbers
          read(10, nml=Qvib, iostat=ios)
          if (ios .ne. 0) write(sal_unit, *) "Namelist Qvib not found"
          write(sal_unit, nml=Qvib)
      end subroutine read_Data4NM

      subroutine compute_Qmax(T, Qmax)
          use constants, only: h, c, kb, autocm_1
          implicit none
          real(dp), intent(in) :: T
          integer, intent(out) :: Qmax(nfreqs)
          integer :: ifreq, iv
          real(dp) :: prob, factor
          integer, parameter :: maxv = 100
          real(dp), parameter :: threshold = 1.e-3_dp

          Qmax = 0
          do ifreq=1, nfreqs
              factor = h *  freqs(ifreq) * autocm_1 * 100 * c / (kb * T)
              do iv=1,maxv
                  prob = exp(-iv * factor) * (1. - exp(-factor))
                  if (prob .lt. threshold) exit
              end do
              Qmax(ifreq) = iv
          end do
          write(sal_unit, *) "Qmax =", Qmax
      end subroutine
end module initial_conditions
