subroutine setpotxyz
    write(*,*) "quadratic potential for 3 particles"
    return
end subroutine

subroutine potxyz(pos, ener, der)
    implicit none
    integer, parameter :: nat = 3
    integer :: i
    real(8), intent(in) :: pos(nat * 3)
    real(8), intent(out) :: ener, der(nat * 3)
    real(8) :: ener1, ener2, pos_(nat*3)
    real(8), parameter :: delta = 1e-5

    call quad_pot(pos, ener)
    do i=1,3*nat
        pos_ = pos 
        pos_(i) = pos(i) + delta
        call quad_pot(pos_, ener1)
        pos_(i) = pos(i) - delta
        call quad_pot(pos_, ener2)
        der(i) = (ener1 - ener2) / ( 2.d0 * delta)
    end do
end subroutine

subroutine quad_pot(pos, ener)
    implicit none
    integer, parameter :: nat = 3
    real(8), intent(in) :: pos(nat * 3)
    real(8), intent(out) :: ener
    real(8) :: r(3)
    real(8), parameter :: alpha = 5.0d-3
    r(1) = sqrt(sum((pos(1:3) - pos(4:6))**2))
    r(2) = sqrt(sum((pos(1:3) - pos(7:9))**2))
    r(3) = sqrt(sum((pos(4:6) - pos(7:9))**2))

    ener = alpha * ( (r(1) - 1.0d0)**2 + (r(2) - 2.0d0)**2 + (r(3) - 3.0d0)**2)
end subroutine
