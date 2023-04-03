module settings
      use constants, only: dp, autouma, sal_unit
      implicit none
      integer, public :: nA, ndim, nfreqs, propagation_mode, initcond_mode
      integer :: ios
      integer, allocatable :: Qnum(:)
      real(dp), allocatable :: XP(:), XPini(:), massA(:), Xeq(:), CXQ(:,:), freqs(:), amp(:)
      character(len=2), allocatable :: atnameA(:)

      namelist /mass/ &
            massA, &
            atnameA
      contains
      subroutine initial_settings()
          implicit none
          ndim = 2 * 3 * nA
          allocate(XP(ndim), XPini(ndim), massA(nA), atnameA(nA), Xeq(ndim/2))

          read(10, nml=mass, iostat=ios)
          if (ios .ne. 0) write(sal_unit, *) "Namelist mass not found"
          write(sal_unit, nml=mass)
          massA = massA/autouma
      end subroutine initial_settings

end module
