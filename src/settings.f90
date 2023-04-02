module settings
      use constants, only: dp, autouma, sal_unit
      implicit none
      integer, public :: nA, ndim, nfreqs, propagation_mode, initcond_mode
      real(dp), allocatable :: XP(:), XPini(:), massA(:), Xeq(:), CXQ(:,:), freqs(:)
      character(len=2), allocatable :: atnameA(:)

      namelist /mass/ &
            massA, &
            atnameA

      contains
      subroutine initial_settings()
          implicit none
          ndim = 2 * 3 * nA
          allocate(XP(ndim), XPini(ndim), massA(nA), atnameA(nA), Xeq(ndim/2))

          read(10, nml=mass)
          write(sal_unit, nml=mass)
          massA = massA/autouma

          if (propagation_mode == 1) then
              ! Adiabatic Switching
              open(11, file="Data4NM.dat", status="old")
              call read_Data2NM(11)
              close(11)
          end if
      end subroutine initial_settings

      subroutine read_Data2NM(read_unit)
          use constants, only: autoA
          implicit none
          integer, intent(in) :: read_unit
          integer :: iat, ix

          read(read_unit, *)
          read(read_unit, *)
          do iat=1,nA
              read(read_unit, *) atnameA(iat), Xeq(3*(iat-1)+1:3*iat)
          end do
          Xeq = Xeq / autoA
          read(read_unit, *) nfreqs
          allocate(freqs(nfreqs), CXQ(3*nA, nfreqs))
          read(read_unit, *) freqs
          do ix=1, 3*nA
              read(read_unit,*) CXQ(ix,:)
          end do
      end subroutine read_Data2NM
end module
