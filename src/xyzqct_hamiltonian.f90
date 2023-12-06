module xyzqct_hamiltonian
    use xyzqct_constants, only: dp, sal_unit, pi, as_unit, autofs
    use xyzqct_settings, only: Ts
    use ddeabm_module, wp => ddeabm_rk
    use xyzqct_settings, only: ndim, nat, mass
    implicit none
    namelist /adiabatic_switch/ &
        Ts

    abstract interface
        subroutine potential_base(t, posxyz, pot, derxyz)
            use xyzqct_constants, only: dp
            use xyzqct_settings, only: ndim
            real(dp), intent(in) :: t, posxyz(ndim/2)
            real(dp), intent(out) :: pot, derxyz(ndim/2)
        end subroutine potential_base
    end interface
    procedure(potential_base), pointer :: potential => null()

    contains

    subroutine get_potential(mode)
        integer, intent(in) :: mode
        integer :: ios
        write(sal_unit, "(/A)") "************************************"
        write(sal_unit, *) "SETTING POTENTIAL..."
        call setpotxyz
        select case(mode)
            case(0)
                write(sal_unit,"(/A/)") "Using user defined potential."
                potential => userpot
            case(1)
                write(sal_unit,"(/A/)") "Using NM potential."
                potential => NMpotential
            case(2)
                write(sal_unit,"(/A/)") "Adiabatic switching"
                potential => adiabatic_switching
                Ts=0._dp
                rewind(10)
                read(10, nml=adiabatic_switch, iostat=ios)
                if (ios .ne. 0) then
                    write(sal_unit,"(/A/)") "Namelist adiabatic_switch not found"
                    write(sal_unit,"(/A/)") "You must define the switching time in it."
                    stop
                end if
                write(sal_unit, nml=adiabatic_switch)
                Ts = Ts/autofs
                open(as_unit, file="initial_conditions", status="replace")
            case default
                write(sal_unit,"(/A/)") "Unknown potential mode."
                stop
        end select
        write(sal_unit, "(/A)") "************************************"
    end subroutine get_potential

    subroutine derivs(me, t, XP, XPder)
        implicit none
        class(ddeabm_class), intent(inout) :: me
        integer ::iat, ix
        real(dp), intent(in) :: t, XP(:)
        real(dp), intent(out) :: XPder(:)
        real(dp) :: pot, posxyz(ndim/2), P(ndim/2), derxyz(ndim/2)

        XPder = 0._dp
        posxyz = 0._dp

        posxyz = XP(:ndim/2)
        P = XP(ndim/2+1:)
        call potential(t, posxyz, pot, derxyz)
        do iat=1,nat
            do ix=1,3
                XPder(3 * (iat-1) + ix) = P(3 * (iat-1) + ix) / mass(iat)
                XPder(ndim/2 + 3 * (iat-1) + ix) = -derxyz(3 * (iat-1) + ix)
            end do
        end do
    end subroutine derivs

    subroutine userpot(t, posxyz, pot, derxyz)
        implicit none
        real(dp), intent(in) :: t, posxyz(ndim/2)
        real(dp), intent(out) :: pot, derxyz(ndim/2)

        pot = 0._dp
        derxyz = 0._dp
        call potxyz(posxyz, pot, derxyz)
    end subroutine userpot

    subroutine NM_pot(posxyz, pot)
        use  xyzqct_settings, only: sysA
        use xyzqct_physics, only: matrix_rotation, rotate_to_eckart
        implicit none
        integer :: ifreq, i, j
        real(dp), intent(in) :: posxyz(ndim/2)
        real(dp), intent(out) :: pot
        real(dp) :: Q(sysA % nfreqs), q_(ndim/2), mass_(ndim/2), Teck(3,3), rot_pos(ndim/2)

        Q = 0._dp
        q_ = 0._dp
        pot = 0._dp

        call rotate_to_eckart(sysA % nat, posxyz, sysA % Xeq, sysA % mass, Teck)
        rot_pos = posxyz
        call matrix_rotation(3, sysA % nat, rot_pos, Teck)
        do i=1,ndim/2
            mass_(i) = sqrt(sysA % mass((i-1)/ 3 + 1))
            q_(i) = (rot_pos(i) - sysA % Xeq(i)) * mass_(i)
            do ifreq=1, sysA % nfreqs
                    Q(ifreq) = Q(ifreq) + sysA % CXQ(i, ifreq) * q_(i)
            end do
        end do
        pot = 5e-1 * sum(sysA % freqs**2 * Q**2)
    end subroutine

    subroutine NMpotential(t, posxyz, pot, derxyz)
        implicit none
        integer :: i
        real(dp), intent(in) :: t, posxyz(ndim/2)
        real(dp), intent(out) :: pot, derxyz(ndim/2)
        real(dp) :: pos_(ndim/2), pot1, pot2
        real(dp), parameter :: delta = 1e-6

        pot = 0._dp
        derxyz = 0._dp
        call NM_pot(posxyz, pot)
        do i=1,ndim/2
            pos_ = posxyz
            pos_(i) = posxyz(i) + delta
            call NM_pot(pos_, pot1)
            pos_(i) = posxyz(i) - delta
            call NM_pot(pos_, pot2)
            derxyz(i) = (pot1 - pot2) / (2 * delta)
        end do

    end subroutine NMpotential

    !subroutine NMpotential_ana(t, posxyz, pot, derxyz)
    !    use  xyzqct_settings, only: nfreqs, freqs, CXQ, Xeq
    !    implicit none
    !    integer :: ifreq, i, j
    !    real(dp), intent(in) :: t, posxyz(ndim/2)
    !    real(dp), intent(out) :: pot, derxyz(ndim/2)
    !    real(dp) :: Q(nfreqs), aux, q_(ndim/2), mass_(ndim/2)

    !    Q = 0._dp
    !    q_ = 0._dp
    !    pot = 0._dp
    !    derxyz = 0._dp

    !    do i=1,ndim/2
    !        mass_(i) = sqrt(mass((i-1)/ 3 + 1))
    !        q_(i) = (posxyz(i) - Xeq(i)) * mass_(i)
    !        do ifreq=1, nfreqs
    !                Q(ifreq) = Q(ifreq) + CXQ(i, ifreq) * q_(i)
    !        end do
    !    end do

    !    do i =1, ndim/2
    !        do ifreq=1, nfreqs
    !            derxyz(i) = derxyz(i) + CXQ(i, ifreq) * freqs(ifreq)**2 * Q(ifreq)
    !        end do
    !        derxyz(i) = derxyz(i) * mass_(i)
    !    end do
    !    pot = 5e-1 * sum(freqs**2 * Q**2)
    !end subroutine NMpotential_ana

    subroutine adiabatic_switching(t, posxyz, pot, derxyz)
        use xyzqct_settings, only: Ts
        implicit none
        real(dp), intent(in) :: t, posxyz(ndim/2)
        real(dp), intent(out) :: pot, derxyz(ndim/2)
        real(dp) :: der1(ndim/2), der2(ndim/2), pot1, pot2, s

        call userpot(t, posxyz, pot1, der1)
        if (t <= Ts) then
            call NMpotential(t, posxyz, pot2, der2)
            s = t / Ts - 1._dp / (2._dp * pi) * sin(2._dp * pi * t / Ts)
        else
            pot2 = 0._dp
            der2 = 0._dp
            s = 1._dp
        end if

        pot = (1._dp - s) * pot2 + s * pot1
        derxyz = (1._dp - s) * der2 + s * der1
    end subroutine adiabatic_switching

    subroutine kinetic_ener(P, E)
        use xyzqct_constants, only: dp, sal_unit
        implicit none
        integer :: iat, ix
        real(dp), intent(in) :: P(ndim/2)
        real(dp), intent(out) :: E

        E = 0._dp

        do iat=1,nat
            do ix=1,3
                E = E + P(3*(iat-1)+ix)**2 / mass(iat)
            end do
        end do
        E = E / 2._dp
    end subroutine

    subroutine total_ener(t, XP, k, pot)
        use xyzqct_constants, only: dp
        use xyzqct_settings, only: ndim
        implicit none
        real(dp), intent(in) :: t, XP(ndim)
        real(dp), intent(out) :: k, pot
        real(dp) :: der(ndim/2)

        call kinetic_ener(XP(ndim/2+1:), k)
        call potential(t, XP(:ndim/2), pot, der)
    end subroutine
end module xyzqct_hamiltonian
