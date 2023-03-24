module settings
      use constants, only: dp
      implicit none
      integer, public :: nA, ndim
      real(dp), allocatable :: massA(:)
      character(len=2), allocatable :: atnameA(:)
end module
