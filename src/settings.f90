module settings
      use constants, only: dp, autouma, sal_unit, get_mass
      implicit none
      integer, public :: nA, ndim, nfreqs, potential_mode, initcond_mode, propagator_mode
      integer :: ios
      integer, allocatable :: Qnum(:), Qmax(:)
      real(dp), allocatable :: XP(:), XPini(:), massA(:), Xeq(:), CXQ(:,:), freqs(:), amp(:)
      real(dp) :: Ts, temperature
      character(len=2), allocatable :: atnameA(:)

      namelist /systemA/ &
            massA, &
            atnameA
      contains
      subroutine initial_settings()
          implicit none
          integer :: iat
          ndim = 2 * 3 * nA
          allocate(XP(ndim), XPini(ndim), massA(nA), atnameA(nA), Xeq(ndim/2))
          massA = 0._dp

          call setpotxyz
          read(10, nml=systemA, iostat=ios)
          if (ios .ne. 0) write(sal_unit, *) "Namelist systemA not found"
          do iat=1, nA
             if (massA(iat) == 0._dp) then
                 call get_mass(atnameA(iat), massA(iat))
                 if (massA(iat) == -1.0_dp) stop
             end if
          end do
          write(sal_unit, nml=systemA)
          massA = massA/autouma
      end subroutine initial_settings

end module
