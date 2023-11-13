module xyzqct_physics
    use xyzqct_constants, only: dp
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

    subroutine get_inertia_moments(nat, pos, mass, inertia, inertia_mat)
        use xyzqct_lapack, only: eigh
        implicit none
        integer, intent(in) :: nat
        real(dp), intent(in) :: pos(3*nat), mass(nat)
        real(dp), intent(out) :: inertia(3), inertia_mat(3,3) ! inertia_mat holds the eigenvector of the inertia matrix on out
        integer :: iat

        inertia_mat = 0._dp
        do iat=1,nat
            inertia_mat(1,1) = inertia_mat(1,1) + mass(iat) * (pos(3*(iat-1)+2)**2 + pos(3*(iat-1)+3)**2)
            inertia_mat(2,2) = inertia_mat(2,2) + mass(iat) * (pos(3*(iat-1)+1)**2 + pos(3*(iat-1)+3)**2)
            inertia_mat(3,3) = inertia_mat(3,3) + mass(iat) * (pos(3*(iat-1)+1)**2 + pos(3*(iat-1)+2)**2)
            inertia_mat(1,2) = inertia_mat(1,2) - mass(iat) * pos(3*(iat-1)+1) * pos(3*(iat-1)+2)
            inertia_mat(1,3) = inertia_mat(1,3) - mass(iat) * pos(3*(iat-1)+1) * pos(3*(iat-1)+3)
            inertia_mat(2,3) = inertia_mat(2,3) - mass(iat) * pos(3*(iat-1)+2) * pos(3*(iat-1)+3)
        end do
        inertia_mat(2,1) = inertia_mat(1,2)
        inertia_mat(3,1) = inertia_mat(1,3)
        inertia_mat(3,2) = inertia_mat(2,3)
        call eigh(3, inertia_mat, inertia)
    end subroutine

    subroutine matrix_rotation(n,m, A, R)
        implicit none
        integer, intent(in) :: n, m
        real(dp), intent(inout) :: A(n*m)
        real(dp), intent(in) :: R(n,n)
        real(dp) :: aux(m,n), tmp(m,n)
        integer :: i, j

        do i=1,m ! nat
            do j=1,n ! ncoor
                aux(i,j) = A(3 * (i-1)+j)
            end do
        end do
        tmp = matmul(aux, R)
        do i=1,m
            do j=1,n
                A(3*(i-1)+j) = tmp(i,j)
            end do
        end do
    end subroutine

    subroutine get_angular_velocity(inertia, AMOM, omega)
        implicit none
        real(dp), intent(in) :: inertia(3), AMOM(3)
        real(dp), intent(out) :: omega(3)
        real(dp) :: inertia_inv(3)

        inertia_inv = 0._dp
        where (inertia > 1e-10_dp)
            inertia_inv = 1._dp / inertia
        end where
        omega = AMOM * inertia_inv
    end subroutine

    subroutine add_angular_velocity(n, XP, mass, omega, sg)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(inout) :: XP(3*2*n)
        real(dp), intent(in) :: mass(n), omega(3), sg
        real(dp) :: x(3), rotvel(3)
        integer :: iat

        do iat=1,n
            x = XP(3*(iat-1)+1:3*iat)
            call cross_prod(omega, x, rotvel)
            XP(3*n+3*(iat-1)+1:3*n+3*iat) = XP(3*n+3*(iat-1)+1:3*n+3*iat) + rotvel * mass(iat) * sg
        end do
    end subroutine

    subroutine rotate_euler(nat, XP, phi, theta, chi)
        implicit none
        integer, intent(in) :: nat
        real(dp), intent(inout) :: XP(3*2*nat)
        real(dp), intent(in) :: phi, theta, chi
        real(dp) :: EulerR(3,3), c(3), s(3)
        c(1) = cos(phi)
        c(2) = cos(theta)
        c(3) = cos(chi)
        s(1) = sin(phi)
        s(2) = sin(theta)
        s(3) = sin(chi)

        EulerR(1,1) = c(1) * c(2) * c(3) - s(1) * s(3)
        EulerR(1,2) = s(1) * c(2) * c(3) + c(1) * s(3)
        EulerR(1,3) = -s(2) * c(3)
        EulerR(2,1) = -c(1) * c(2) * s(3) - s(1) * c(3)
        EulerR(2,2) = c(1) * c(3) - c(2) * s(1) * s(3)
        EulerR(2,3) = s(2) * s(3)
        EulerR(3,1) = c(1) * s(2)
        EulerR(3,2) = s(1) * s(2)
        EulerR(3,3) = c(2)
        call matrix_rotation(3, nat, XP(:3*nat), EulerR)
        call matrix_rotation(3, nat, XP(3*nat+1:), EulerR)
    end subroutine
end module xyzqct_physics
