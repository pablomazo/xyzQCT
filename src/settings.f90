module settings
      use constants, only: dp
      implicit none
      integer, public :: nA, ndim
      real(dp), allocatable :: XP(:), XPini(:), massA(:)
      character(len=2), allocatable :: atnameA(:)

      contains
      subroutine allocate_variables()
          implicit none
          allocate(XP(ndim), XPini(ndim), massA(nA), atnameA(nA))
      end subroutine allocate_variables
end module
