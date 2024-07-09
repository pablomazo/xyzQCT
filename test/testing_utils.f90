module testutils
    implicit none
    private
    public :: all_close

    contains

    logical function all_close(x, y, tol) result(res)
        implicit none
        real(8), intent(in) :: x(:), y(:), tol
        res = all(abs(x-y) < tol)
    end function
end module
