module xyzqct_lapack
   use xyzqct_constants, only: dp
   implicit none

contains

   subroutine eigh(s, A, w)
      implicit none
      integer, intent(in) :: s
      real(dp), intent(inout) :: A(s, s)
      real(dp), intent(out) :: w(s)
      real(dp) :: n_work_arr(1)
      real(dp), allocatable :: work_arr(:)
      integer :: error_flag

      call DSYEV("V", 'U', s, A, s, w, n_work_arr, -1, error_flag)
      allocate (work_arr(nint(n_work_arr(1))))
      call DSYEV("V", 'U', s, A, s, w, work_arr(1), size(work_arr), error_flag)
   end subroutine
end module
