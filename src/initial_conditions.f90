module initial_conditions
    use constants, only: dp, sal_unit, cond_unit, pi
    use settings, only : ndim, Qnum, amp, nfreqs, freqs, CXQ, massA, Xeq, atnameA, nA
    implicit none
    character(len=80) :: initcond_file
    integer :: init_cond_mode, max_cond

    namelist /Qvib/ &
        Qnum


    ! init_cond_mode:
    ! 0 => Read from file
    ! 1 => NM initial condition

    abstract interface
        subroutine get_init_cond_base(XP)
            use constants, only : dp 
            use settings, only : ndim
            real(dp), intent(out) :: XP(ndim)
        end subroutine get_init_cond_base
    end interface
    procedure(get_init_cond_base), pointer :: get_init_cond => null()

    contains
        subroutine set_init_cond(mode)
            implicit none
            integer :: mode

            select case(mode)
                case(0)
                    write(sal_unit, "(/A)") "Reading initial conditions from file:", trim(initcond_file)
                    open(cond_unit, file=trim(initcond_file), status="old")
                    read(cond_unit,*) max_cond
                    rewind(cond_unit)
                    write(sal_unit,*) "Total number of initial conditions =", max_cond
                    get_init_cond => from_file_init_cond
                case(1)
                    write(sal_unit, "(/A)") "Using NM initial conditions"
                    open(11, file="Data4NM.dat", status="old")
                    call read_Data4NM(11)
                    close(11)
                    get_init_cond => NM_init_cond
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
            implicit none
            integer :: ix, ifreq
            real(dp), intent(out) :: XP(ndim)
            real(dp) :: phase(nfreqs), Q(nfreqs), P(nfreqs), mass

            call RANDOM_NUMBER(phase)
            XP = 0._dp
            phase = 2._dp * pi * phase
            phase = 0._dp
            write(sal_unit, *) "Initial phase:", phase
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

      subroutine read_Data4NM(read_unit)
          use constants, only: autoA, autocm_1
          implicit none
          integer, intent(in) :: read_unit
          integer :: iat, ix, ifreq, ios

          read(read_unit, *)
          read(read_unit, *)
          do iat=1,nA
              read(read_unit, *) atnameA(iat), Xeq(3*(iat-1)+1:3*iat)
          end do
          Xeq = Xeq / autoA
          read(read_unit, *) nfreqs
          allocate(freqs(nfreqs), CXQ(3*nA, nfreqs), Qnum(nfreqs), amp(nfreqs))
          read(read_unit, *) freqs
          do ix=1, 3*nA
              read(read_unit,*) CXQ(ix,:)
          end do
          freqs = freqs / autocm_1
          Qnum = 0 ! Vibrational quantum numbers
          read(10, nml=Qvib, iostat=ios)
          if (ios .ne. 0) write(sal_unit, *) "Namelist Qvib not found"
          write(sal_unit, nml=Qvib)
          amp = sqrt((2._dp * Qnum +1._dp) / freqs) ! Maximum NM amplitudes
          write(sal_unit, *) "Qnum, freqs, amp "
          do ifreq=1,nfreqs
              write(sal_unit, *) Qnum(ifreq), freqs(ifreq), amp(ifreq)
          end do
      end subroutine read_Data4NM
end module initial_conditions
