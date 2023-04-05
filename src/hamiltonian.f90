module hamiltonian
    use constants, only: dp, sal_unit
    use ddeabm_module, wp => ddeabm_rk
    use settings, only: ndim, nA, massA
    implicit none

    abstract interface
        subroutine potential_base(t, posxyz, pot, derxyz)
            use constants, only: dp
            use settings, only: ndim
            real(dp), intent(in) :: t, posxyz(ndim/2)
            real(dp), intent(out) :: pot, derxyz(ndim/2)
        end subroutine potential_base
    end interface
    procedure(potential_base), pointer :: potential => null()

    contains

    subroutine get_potential(mode)
        integer, intent(in) :: mode
        select case(mode)
            case(0)
                write(sal_unit,"(/A/)") "Using user defined potential."
                potential => userpot
            case(1)
                write(sal_unit,"(/A/)") "Using NM potential."
                potential => NMpotential
        end select
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
        do iat=1,nA
            do ix=1,3
                XPder(3 * (iat-1) + ix) = P(3 * (iat-1) + ix) / massA(iat)
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

    subroutine NMpotential(t, posxyz, pot, derxyz)
        use  settings, only: nfreqs, freqs, CXQ, Xeq
        implicit none
        integer :: ifreq, i, j
        real(dp), intent(in) :: t, posxyz(ndim/2)
        real(dp), intent(out) :: pot, derxyz(ndim/2)
        real(dp) :: Q(nfreqs), P(nfreqs), aux, q_(ndim/2), mass(ndim/2)

        Q = 0._dp
        q_ = 0._dp
        P = 0._dp
        pot = 0._dp
        derxyz = 0._dp

        do i=1,ndim/2
            mass(i) = sqrt(massA((i-1)/ 3 + 1))
            q_(i) = (posxyz(i) - Xeq(i)) * mass(i)
        end do
        do i =1, ndim/2
            do ifreq=1, nfreqs
                do j=1, ndim/2
                    derxyz(i) = derxyz(i) + CXQ(j, ifreq) * q_(j) * freqs(ifreq)**2 * CXQ(i, ifreq)
                end do
            end do
            derxyz(i) = derxyz(i) * mass(i)
        end do
        do ifreq=1, nfreqs
            aux = 0._dp
            do i=1,ndim/2
                aux = aux + CXQ(i, ifreq) * q_(i)
            end do
            pot = pot + (freqs(ifreq) * aux)**2
        end do
        pot = pot / 2._dp
    end subroutine NMpotential

    subroutine kinetic_ener(P, E)
        use constants, only: dp, sal_unit
        use settings, only: ndim, nA
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
end module hamiltonian
