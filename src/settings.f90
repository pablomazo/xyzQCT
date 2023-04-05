module settings
      use constants, only: dp, autouma, sal_unit
      implicit none
      integer, public :: nA, ndim, nfreqs, propagation_mode, initcond_mode
      integer :: ios
      integer, allocatable :: Qnum(:)
      real(dp), allocatable :: XP(:), XPini(:), massA(:), Xeq(:), CXQ(:,:), freqs(:), amp(:)
      real(dp) :: Ts
      character(len=2), allocatable :: atnameA(:)

      namelist /systemA/ &
            massA, &
            atnameA
      contains
      subroutine initial_settings()
          implicit none
          ndim = 2 * 3 * nA
          allocate(XP(ndim), XPini(ndim), massA(nA), atnameA(nA), Xeq(ndim/2))

          read(10, nml=systemA, iostat=ios)
          if (ios .ne. 0) write(sal_unit, *) "Namelist systemA not found"
          write(sal_unit, nml=systemA)
          massA = massA/autouma
      end subroutine initial_settings

end module
