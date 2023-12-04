module xyzqct_settings
      use xyzqct_constants, only: dp, autouma, sal_unit, get_mass
      implicit none
      integer, public :: nA, nB, nat, ndim, nfreqs, potential_mode, initcond_mode, propagator_mode
      integer :: ios
      integer, allocatable :: Qnum(:), Qmax(:)
      real(dp), allocatable :: XP(:), XPini(:), massA(:), Xeq(:), CXQ(:,:), freqs(:), amp(:), &
          massB(:), mass(:), XeqA(:), XeqB(:)
      real(dp) :: Ts, temperature, rfin
      character(len=2), allocatable :: atnameA(:), atnameB(:), atname(:)
      character(len=80) :: initcond_fileA, initcond_fileB

      namelist /systemA/ &
            massA, &
            atnameA, &
            initcond_fileA
      namelist /systemB/ &
            massB, &
            atnameB, &
            initcond_fileB
      contains
      subroutine initial_settings()
          implicit none
          integer :: iat
          nat = nA + nB
          ndim = 2 * 3 * (nA + nB)
          allocate(XP(ndim), XPini(ndim), massA(nA), atnameA(nA), Xeq(ndim/2), &
              massB(nB), atnameB(nB), mass(nat), atname(nat))

          ! Defaults
          initcond_fileA = ""
          initcond_fileB = ""
          massA = 0._dp
          massB = 0._dp

          call setpotxyz
          read(10, nml=systemA, iostat=ios)
          if (ios .ne. 0) write(sal_unit, *) "Namelist systemA not found"
          read(10, nml=systemB, iostat=ios)
          if (ios .ne. 0) write(sal_unit, *) "Namelist systemB not found"
          do iat=1, nA
             if (massA(iat) == 0._dp) then
                 call get_mass(atnameA(iat), massA(iat))
                 if (massA(iat) == -1.0_dp) stop
             end if
          end do
          do iat=1, nB
             if (massB(iat) == 0._dp) then
                 call get_mass(atnameB(iat), massB(iat))
                 if (massB(iat) == -1.0_dp) stop
             end if
          end do
          write(sal_unit, nml=systemA)
          if (ios == 0) write(sal_unit, nml=systemB)
          massA = massA/autouma
          massB = massB/autouma
          mass(:nA) = massA
          mass(nA+1:) = massB
          atname(:nA) = atnameA
          atname(nA+1:) = atnameB
      end subroutine initial_settings

end module
