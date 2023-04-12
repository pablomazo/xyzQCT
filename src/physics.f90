module physics
    use constants, only: dp
    contains

    subroutine get_COM(ndim, XP, nat1, nat2, mass, QCOM, PCOM)
        implicit none
        integer, intent(in) :: ndim, nat1, nat2
        real(dp), intent(in) :: XP(ndim), mass(nat2 - nat1 + 1)
        real(dp), intent(out) :: QCOM(3), PCOM(3)
        real(dp) :: mtot
        integer :: ix, iat

        QCOM = 0._dp
        PCOM = 0._dp
        mtot = 0._dp
        do iat=nat1, nat2
            mtot = mtot + mass(iat)
            do ix=1, 3
                QCOM(ix) = QCOM(ix) + XP(3 * (iat-1) + ix) * mass(iat)
                PCOM(ix) = PCOM(ix) + XP(3 * (iat-1) + ix + ndim / 2)
            end do
        end do
        QCOM = QCOM / mtot
    end subroutine get_COM

    subroutine get_LMOM_AMOM(ndim, XP, nat1, nat2, mass, QCOM, PCOM, LMOM, AMOM)
        implicit none
        integer, intent(in) :: ndim, nat1, nat2
        real(dp), intent(in) :: XP(ndim), mass(nat2 - nat1 + 1), QCOM(3), PCOM(3)
        real(dp), intent(out) :: LMOM(3), AMOM(3)
        real(dp) :: mtot, Q(3), P(3), prod(3)
        integer :: iat, ix

        mtot = 0._dp
        LMOM = 0._dp 
        AMOM = 0._dp
        do iat=nat1, nat2
            mtot = mtot + mass(iat)
        end do

        do iat=nat1, nat2
            do ix=1,3
                Q(ix) = XP(3 * (iat-1) + ix) - QCOM(ix)
                P(ix) = &
                    XP(ndim/2 + 3 * (iat-1) + ix) - PCOM(ix) * mass(iat) / mtot
                LMOM(ix) = LMOM(ix) + XP(ndim/2 + 3 * (iat-1) + ix)
            end do
            call cross_prod(Q, P, prod)
            AMOM = AMOM + prod
        end do
    end subroutine get_LMOM_AMOM

    subroutine cross_prod(v, w, p)
        implicit none
        real(dp), intent(in) :: v(3), w(3)
        real(dp), intent(out) :: p(3)

        p(1) = v(2) * w(3) - w(2) * v(3)
        p(2) = v(3) * w(1) - w(3) * v(1)
        p(3) = v(1) * w(2) - w(1) * v(2)

    end subroutine cross_prod
end module physics
