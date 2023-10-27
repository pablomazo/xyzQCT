module xyzqct_initial_conditions
    use xyzqct_constants, only: dp, sal_unit, end_unit, cond_unitA, cond_unitB, autoeV, pi
    use xyzqct_settings, only : ndim, Qnum, amp, nfreqs, freqs, CXQ, massA, Xeq, atnameA, nA, Qmax, temperature, rfin
    implicit none
    integer :: init_cond_mode, max_condA, max_condB
    real(dp) :: Ecoll, Trot, rini, bmax, bmin, bparam

    private :: init_cond_mode, bparam

    namelist /Qvib/ &
        Qnum
    namelist /collision/ &
        Ecoll, &
        Trot, &
        rini, &
        bmax, &
        bmin

    ! init_cond_mode:
    ! 0 => Read from file
    ! 1 => NM initial condition
    ! 2 => NM initial condition (sampled at temperature T)

    contains
        subroutine set_init_cond(mode)
            use xyzqct_settings, only: initcond_fileA, initcond_fileB
            implicit none
            integer :: mode
            init_cond_mode = mode
            select case(mode)
                case(0)
                    write(sal_unit, "(/A)") "Reading initial conditions from file:", trim(initcond_fileA)
                    open(cond_unitA, file=trim(initcond_fileA), status="old")
                    read(cond_unitA,*) max_condA
                    rewind(cond_unitA)
                    write(sal_unit,*) "Total number of initial conditions =", max_condA
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
                case(3)
                    write(sal_unit, "(/A)") "Using initial conditions for A+B collision"
                    Ecoll = 0._dp
                    Trot = 0._dp
                    bmin = 0._dp
                    read(10, nml=collision)
                    write(sal_unit, nml=collision)
                    Ecoll = Ecoll/autoeV
                    if (rini > rfin) then
                        write(sal_unit,*) "ERROR: rini must be smaller than rfin"
                        stop
                    end if

                    if (bmax > rini) then
                        write(sal_unit,*) "ERROR: bmax must be smaller or equal than rfin"
                        stop
                    end if

                    open(cond_unitA, file=trim(initcond_fileA), status="old")
                    read(cond_unitA,*) max_condA
                    rewind(cond_unitA)
                    write(sal_unit,*) "Total number of initial conditions (A) =", max_condA
                    open(cond_unitB, file=trim(initcond_fileB), status="old")
                    read(cond_unitB,*) max_condB
                    rewind(cond_unitB)
                    write(sal_unit,*) "Total number of initial conditions (B) =", max_condB
                case default
                    write(sal_unit,"(/A/)") "Unknown initial condition mode."
                    stop
            end select
            ! write info to end_cond file
            select case(init_cond_mode)
                case (:2)
                    write(end_unit,*) "itraj, XP0, XP, time"
                case (3)
                    write(end_unit,*) "itraj, XP0, XP, time, bmax"
            end select
        end subroutine

        subroutine write_end_cond(itraj, time, XP0, XP)
            implicit none
            integer :: itraj
            real(dp), intent(in) :: time, XP0(ndim), XP(ndim)
            select case(init_cond_mode)
                case (:2)
                    write(end_unit,*) itraj, XP0, XP, time
                case (3)
                    write(end_unit,*) itraj, XP0, XP, time, bparam
            end select
        end subroutine

        subroutine get_init_cond(XP)
            implicit none
            real(dp), intent(out) :: XP(ndim)
            select case(init_cond_mode)
                case(0)
                    call from_file_init_cond(max_condA, cond_unitA, nA, XP)
                case(1)
                    call NM_init_cond(XP)
                case(2)
                    call NM_init_cond_T(XP)
                case(3)
                    call AplusB_init_cond(XP)
            end select 
        end subroutine

        subroutine from_file_init_cond(max_cond, cond_unit, n, XP)
            implicit none
            integer :: max_cond, cond_unit, n
            real(dp), intent(out) :: XP(3*2*n)
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
            use xyzqct_constants, only: pi
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
            use xyzqct_constants, only: h, c, kb, autocm_1
            use xyzqct_settings, only: Qmax, temperature
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
          use xyzqct_constants, only: autoA, autocm_1
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
          use xyzqct_constants, only: h, c, kb, autocm_1
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

      subroutine AplusB_init_cond(XP)
          use xyzqct_physics, only: get_COM, get_inertia_moments, matrix_rotation, get_angular_velocity, &
              add_angular_velocity, get_LMOM_AMOM, rotate_euler
          use xyzqct_settings, only: nA, nB, massA, massB, nat, mass
          implicit none
          real(dp), intent(out) :: XP(ndim)
          real(dp) :: XPA(3*2*nA), XPB(3*2*nB), QCOM(3), PCOM(3), mtot, r, ang, &
              inertia(3), inertia_vec(3,3), LMOM(3), AMOM(3), omega(3), phi, theta, chi, mred, mA, mB
          integer :: iat

          XP = 0.0_dp
          call from_file_init_cond(max_condA, cond_unitA, nA, XPA)
          call from_file_init_cond(max_condB, cond_unitB, nB, XPB)

          mA = sum(massA)
          mB = sum(massB)
          mtot = mA + mB
          !-------------------------------
          ! Remove COM and momentum (A and B)
          call get_COM(3*2*nA, XPA, 1, nA, massA, QCOM, PCOM)
          do iat=1,nA
              XPA(3*(iat-1)+1:3*iat) = XPA(3*(iat-1)+1:3*iat) - QCOM
              XPA(3*nA+3*(iat-1)+1:3*nA+3*iat) = &
                  XPA(3*nA+3*(iat-1)+1:3*nA+3*iat) - PCOM * massA(iat) / mA
          end do
          call get_COM(3*2*nA, XPA, 1, nA, massA, QCOM, PCOM)

          call get_COM(3*2*nB, XPB, 1, nB, massB, QCOM, PCOM)
          do iat=1,nB
              XPB(3*(iat-1)+1:3*iat) = XPB(3*(iat-1)+1:3*iat) - QCOM
              XPB(3*nB+3*(iat-1)+1:3*nB+3*iat) = &
                  XPB(3*nB+3*(iat-1)+1:3*nB+3*iat) - PCOM * massB(iat) / mB
          end do
          call get_COM(3*2*nB, XPB, 1, nB, massB, QCOM, PCOM)
          !-------------------------------

          !-------------------------------
          ! Align with inertia axis and remove angular momentum + euler rotation
          ! A
          call get_inertia_moments(nA, XPA(:3*nA), massA, inertia, inertia_vec)
          call matrix_rotation(3, nA, XPA(:3*nA), inertia_vec)
          call matrix_rotation(3, nA, XPA(3*nA+1:), inertia_vec)
          call get_inertia_moments(nA, XPA(:3*nA), massA, inertia, inertia_vec)
          call get_LMOM_AMOM(3*2*nA, XPA,1, nA, massA, QCOM, PCOM, LMOM, AMOM)
          call get_angular_velocity(inertia, AMOM, omega)
          call add_angular_velocity(nA, XPA, massA, omega, -1._dp)

          call RANDOM_NUMBER(r)
          phi = 2 * pi * r
          call RANDOM_NUMBER(r)
          theta = pi * r
          call RANDOM_NUMBER(r)
          chi = 2 * pi * r
          call rotate_euler(nA, XPA, phi, theta, chi)
          write(sal_unit,*) "Setting euler phi, theta, chi (A) = ", phi, theta, chi

          ! B
          call get_inertia_moments(nB, XPB(:3*nB), massB, inertia, inertia_vec)
          call matrix_rotation(3, nB, XPB(:3*nB), inertia_vec)
          call matrix_rotation(3, nB, XPB(3*nB+1:), inertia_vec)
          call get_inertia_moments(nB, XPB(:3*nB), massB, inertia, inertia_vec)
          call get_LMOM_AMOM(3*2*nB, XPB,1, nB, massB, QCOM, PCOM, LMOM, AMOM)
          call get_angular_velocity(inertia, AMOM, omega)
          call add_angular_velocity(nB, XPB, massB, omega, -1._dp)

          call RANDOM_NUMBER(r)
          phi = 2 * pi * r
          call RANDOM_NUMBER(r)
          theta = pi * r
          call RANDOM_NUMBER(r)
          chi = 2 * pi * r
          call rotate_euler(nB, XPB, phi, theta, chi)
          write(sal_unit,*) "Setting euler phi, theta, chi (B) = ", phi, theta, chi
          !-------------------------------

          !-------------------------------
          ! Move system B to rini with bparam
          call RANDOM_NUMBER(r)
          bparam = bmin + (bmax - bmin) * r
          write(sal_unit,*) "Setting bmax / au = ", bparam
          call RANDOM_NUMBER(r)
          ang = 2 * pi * r
          QCOM = 0.0_dp
          QCOM(1) = bparam * cos(ang)
          QCOM(2) = bparam * sin(ang)
          QCOM(3) = sqrt(abs(rini**2-bparam**2))
          do iat=1,nB
              XPB(3*(iat-1)+1:3*iat) = XPB(3*(iat-1)+1:3*iat) + QCOM
          end do
          !-------------------------------

          !-------------------------------
          ! Add pZ to system b
          write(sal_unit,*) "Setting Ecoll / au = ", Ecoll
          mred = mA * mB / (mA + mB)
          PCOM = 0.0_dp
          PCOM(3)=sqrt(2._dp * mred * Ecoll)
          do iat=1,nA
              XPA(3*nA+3*(iat-1)+1:3*nA+3*iat) = &
                  XPA(3*nA+3*(iat-1)+1:3*nA+3*iat) + PCOM * massA(iat) / mA
          end do
          do iat=1,nB
              XPB(3*nB+3*(iat-1)+1:3*nB+3*iat) = &
                  XPB(3*nB+3*(iat-1)+1:3*nB+3*iat) - PCOM * massB(iat) / mB
          end do
          !-------------------------------

          XP(1:3*nA) = XPA(1:3*nA)
          XP(3*nA+1:3*nA+3*nB) = XPB(1:3*nB)
          XP(3*nA+3*nB+1:3*2*nA+3*nB) = XPA(3*nA+1:)
          XP(3*2*nA+3*nB+1:) = XPB(3*nB+1:)

          ! Set COM at [0,0,0]
          call get_COM(ndim, XP, 1, nat, mass, QCOM, PCOM)
          mtot = sum(mass)
          do iat=1,nat
              XP(3*(iat-1)+1:3*iat) = XP(3*(iat-1)+1:3*iat) - QCOM
          end do
          flush(sal_unit)
      end subroutine
end module xyzqct_initial_conditions
