module xyzqct_initial_conditions
    use xyzqct_constants, only: dp, sal_unit, end_unit, cond_unitA, cond_unitB, autoeV, pi, kb, autoJ
    use xyzqct_settings, only : ndim, Qnum, amp, nfreqs, freqs, CXQ, massA, Xeq, atnameA, nA, Qmax, temperature, rfin
    implicit none
    integer :: init_cond_mode, max_condA, max_condB
    real(dp) :: Ecoll, Trot, rini, bmax, bmin, bparam, Ttrans, A_capture, n_capture

    private :: init_cond_mode, bparam

    namelist /Qvib/ &
        Qnum
    namelist /collision/ &
        Ecoll, &
        Trot, &
        Ttrans, &
        rini, &
        bmax, &
        bmin, &
        A_capture, &
        n_capture

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
                    Ttrans = 0._dp
                    bmin = 0._dp
                    A_capture = 0._dp
                    n_capture = 0._dp
                    rewind(10)
                    read(10, nml=collision)
                    write(sal_unit, nml=collision)
                    Ecoll = Ecoll/autoeV
                    Trot = Trot * kb / autoJ
                    Ttrans = Ttrans * kb / autoJ
                    if (Ecoll > 0._dp .and. Ttrans > 0._dp) then
                        write(sal_unit,*) "ERROR: Either Ecoll or Ttrans must be zero."
                        stop
                    end if
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

        subroutine get_init_cond(XP, propagate)
            implicit none
            real(dp), intent(out) :: XP(ndim)
            logical, intent(out) :: propagate
            propagate = .true.
            select case(init_cond_mode)
                case(0)
                    call from_file_init_cond(max_condA, cond_unitA, nA, XP)
                case(1)
                    call NM_init_cond(XP)
                case(2)
                    call NM_init_cond_T(XP)
                case(3)
                    call AplusB_init_cond(XP, propagate)
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
            Q = amp * sin(phase)
            P = freqs * amp * cos(phase)

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

      subroutine AplusB_init_cond(XP, propagate)
          use xyzqct_physics, only: get_COM, get_inertia_moments, matrix_rotation, get_angular_velocity, &
              add_angular_velocity, get_LMOM_AMOM, rotate_euler
          use xyzqct_settings, only: nA, nB, massA, massB, nat, mass
          implicit none
          real(dp), intent(out) :: XP(ndim)
          logical, intent(out) :: propagate
          real(dp) :: XPA(3*2*nA), XPB(3*2*nB), QCOM(3), PCOM(3), mtot, r, ang, &
              inertia(3), inertia_vec(3,3), LMOM(3), AMOM(3), omega(3), phi, theta, chi, mred, mA, mB, &
              J(4), erel, erelmax, prob, Erot
          integer :: iat

          XP = 0.0_dp
          propagate = .true.
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
          if (nA > 1) then
              call get_inertia_moments(nA, XPA(:3*nA), massA, inertia, inertia_vec)
              call matrix_rotation(3, nA, XPA(:3*nA), inertia_vec)
              call matrix_rotation(3, nA, XPA(3*nA+1:), inertia_vec)
              call get_inertia_moments(nA, XPA(:3*nA), massA, inertia, inertia_vec)
              call get_LMOM_AMOM(3*2*nA, XPA,1, nA, massA, QCOM, PCOM, LMOM, AMOM)
              call get_angular_velocity(inertia, AMOM, omega)
              call add_angular_velocity(nA, XPA, massA, omega, -1._dp)
              if (Trot .gt. 0._dp) then
                  call sample_J(inertia, Trot, J, Erot)
                  write(sal_unit,*) "Setting Jx, Jy, Jz, J (A):", J
                  write(sal_unit,*) "Rotational energy /au (A):", Erot
                  call get_angular_velocity(inertia, J(1:3), omega)
                  call add_angular_velocity(nA, XPA, massA, omega, 1._dp)
              end if

              call RANDOM_NUMBER(r)
              phi = 2 * pi * r
              call RANDOM_NUMBER(r)
              theta = pi * r
              call RANDOM_NUMBER(r)
              chi = 2 * pi * r
              call rotate_euler(nA, XPA, phi, theta, chi)
              write(sal_unit,*) "Setting euler phi, theta, chi (A) = ", phi, theta, chi
          end if

          ! B
          if (nB > 1) then
              call get_inertia_moments(nB, XPB(:3*nB), massB, inertia, inertia_vec)
              call matrix_rotation(3, nB, XPB(:3*nB), inertia_vec)
              call matrix_rotation(3, nB, XPB(3*nB+1:), inertia_vec)
              call get_inertia_moments(nB, XPB(:3*nB), massB, inertia, inertia_vec)
              call get_LMOM_AMOM(3*2*nB, XPB,1, nB, massB, QCOM, PCOM, LMOM, AMOM)
              call get_angular_velocity(inertia, AMOM, omega)
              call add_angular_velocity(nB, XPB, massB, omega, -1._dp)
              if (Trot .gt. 0._dp) then
                  call sample_J(inertia, Trot, J, Erot)
                  write(sal_unit,*) "Setting Jx, Jy, Jz, J (B):", J
                  write(sal_unit,*) "Rotational energy /au (B):", Erot
                  call get_angular_velocity(inertia, J(1:3), omega)
                  call add_angular_velocity(nB, XPB, massB, omega, 1._dp)
              end if

              call RANDOM_NUMBER(r)
              phi = 2 * pi * r
              call RANDOM_NUMBER(r)
              theta = pi * r
              call RANDOM_NUMBER(r)
              chi = 2 * pi * r
              call rotate_euler(nB, XPB, phi, theta, chi)
              write(sal_unit,*) "Setting euler phi, theta, chi (B) = ", phi, theta, chi
          end if
          !-------------------------------

          !-------------------------------
          ! Add pZ to system b
          if (Ecoll .ne. 0._dp) then
              erel = Ecoll
          else
              erelmax = 15._dp * Ttrans
              r = 1.0_dp
              prob = 0.0_dp
              do while (r > prob)
                  call RANDOM_NUMBER(r)
                  erel = erelmax * r
                  prob = erel / Ttrans * exp(-erel / Ttrans) * exp(1.0_dp)
                  call RANDOM_NUMBER(r)
              end do
          end if
          write(sal_unit,*) "Setting relative translational energy / au = ", erel
          mred = mA * mB / (mA + mB)
          PCOM = 0.0_dp
          PCOM(3)=sqrt(2._dp * mred * erel)
          do iat=1,nA
              XPA(3*nA+3*(iat-1)+1:3*nA+3*iat) = &
                  XPA(3*nA+3*(iat-1)+1:3*nA+3*iat) + PCOM * massA(iat) / mA
          end do
          do iat=1,nB
              XPB(3*nB+3*(iat-1)+1:3*nB+3*iat) = &
                  XPB(3*nB+3*(iat-1)+1:3*nB+3*iat) - PCOM * massB(iat) / mB
          end do
          !-------------------------------

          !-------------------------------
          ! Move system B to rini with bparam
          call RANDOM_NUMBER(r)
          bparam = bmin + (bmax - bmin) * r
          if (A_capture .ne. 0._dp) then
              propagate = bparam <= A_capture / erel**n_capture
              if (.not. propagate) then
                  write(sal_unit,*) &
                  "Trajectory will not propagate because bmax is too high according to capture model."
              end if
          end if
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

      subroutine sample_J(inertia, T, J, Erot)
          ! Faraday Discuss. Chem. Soc., 1973,55, 93-99
          ! VENUS manual
          implicit none
          real(dp), intent(in) :: inertia(3), T
          real(dp), intent(out) :: J(4)
          integer :: nlin, ix
          real(dp) :: prob, r, diff1, diff2, Jmax, inertia_mean, Erot

          J = 0.0_dp
          nlin = 0
          diff1 = abs(inertia(1) - inertia(2))
          diff2 = abs(inertia(3) - inertia(2))

          if (diff1 > diff2) then
              ix = 1 ! sample Jx
          else
              ix = 3 ! sample Jz
          end if

          if (abs(inertia(1)) < 1e-10 .and. diff2 < 1e-10) nlin = 1 ! linear molecule

          if (nlin == 0) then
              Jmax = 25._dp * inertia(ix) * T
              r = 1.0_dp
              prob = 0.0_dp
              do while (r > prob)
                  call RANDOM_NUMBER(r)
                  J(ix) = Jmax * r
                  call RANDOM_NUMBER(r)
                  if (r > 0.5_dp) J(ix) = -J(ix)
                  prob = exp(-J(ix)**2 / (2._dp * inertia(ix) * T))
                  call RANDOM_NUMBER(r)
              end do
              Erot = J(ix)**2 / inertia(ix)
          end if

          if (nlin == 1 .or. ix == 1) then
              inertia_mean = sqrt(inertia(2) * inertia(3))
              call RANDOM_NUMBER(r)
              J(4) = sqrt(J(1)**2 - 2 * inertia_mean * T * log(1._dp - r))
              call RANDOM_NUMBER(r)
              J(2) = sqrt(J(4)**2 - J(1)**2) * sin(2._dp * pi * r)
              J(3) = sqrt(J(4)**2 - J(1)**2) * cos(2._dp * pi * r)
              Erot = (J(2)**2 / inertia(2) + J(3)**2 / inertia(3) + Erot) / 2._dp
          else
              inertia_mean = sqrt(inertia(1) * inertia(2))
              call RANDOM_NUMBER(r)
              J(4) = sqrt(J(3)**2 - 2 * inertia_mean * T * log(1._dp - r))
              call RANDOM_NUMBER(r)
              J(1) = sqrt(J(4)**2 - J(3)**2) * sin(2._dp * pi * r)
              J(2) = sqrt(J(4)**2 - J(3)**2) * cos(2._dp * pi * r)
              Erot = (J(1)**2 / inertia(1) + J(2)**2 / inertia(2) + Erot) / 2._dp
          end if
      end subroutine
end module xyzqct_initial_conditions
