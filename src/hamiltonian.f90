module hamiltonian
    use constants, only: dp, sal_unit
    use ddeabm_module, wp => ddeabm_rk
    use settings, only: ndim, nA, massA
    implicit none

    abstract interface 
        subroutine derivs_base(me, t, XP, XPder)
            use constants, only: dp
            import :: ddeabm_class
            class(ddeabm_class), intent(inout) :: me
            real(dp), intent(in) :: t, XP(:)
            real(dp), intent(out) :: XPder(:)
        end subroutine
    end interface
    procedure(derivs_base), pointer :: derivs => null()

    contains

    subroutine get_derivs(mode)
        integer, intent(in) :: mode
        select case(mode)
            case(0)
                write(sal_unit,"(/A/)") "Using regular xyz hamiltonian in propagation."
                derivs => Hxyz
        end select
    end subroutine get_derivs

    subroutine Hxyz(me, t, XP, XPder)
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
    end subroutine Hxyz
end module hamiltonian
