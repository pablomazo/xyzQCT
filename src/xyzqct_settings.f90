module xyzqct_settings
   use xyzqct_constants, only: dp, autouma, sal_unit, get_mass, autoA
   implicit none
   type :: System
      integer :: nat, nfreqs
      real(dp), allocatable :: Xeq(:), mass(:)
      ! For NM
      real(dp), allocatable :: CXQ(:, :), freqs(:)
      integer, allocatable :: Qnum(:), Qmax(:)
      real(dp) :: Tvib
      character(len=2), allocatable :: atname(:)
      character(len=80) :: initcond_file
   end type
   integer, public :: nA, nB, nat, ndim, nfreqs, potential_mode, initcond_mode, propagator_mode
   type(System), public :: sysA, sysB
   real(dp) :: Ts, temperature, rfin
   integer :: ios
   integer, allocatable :: Qnum(:), Qmax(:)
   real(dp), allocatable :: XP(:), XPini(:), mass(:)
   character(len=2), allocatable :: atname(:)

contains
   subroutine initial_settings()
      implicit none
      nat = nA + nB
      ndim = 2*3*nat
      allocate (XP(ndim), XPini(ndim), mass(nat), atname(nat))

      call setup_system('A', nA, sysA)
      call setup_system('B', nB, sysB)
      mass(:nA) = sysA%mass
      mass(nA + 1:) = sysB%mass
      atname(:nA) = sysA%atname
      atname(nA + 1:) = sysB%atname
   end subroutine initial_settings

   subroutine setup_system(s, nat, sys)
      use xyzqct_constants, only: iunit, autoJ, kb
      implicit none
      character(len=1), intent(in) :: s
      integer, intent(in) :: nat
      type(System), intent(inout) :: sys
      integer :: ios, iat
      real(dp), allocatable :: mass(:), Xeq(:)
      real(dp) :: Tvib
      character(len=2), allocatable :: atname(:)
      character(len=80) :: initcond_file
      namelist /systemA/ &
         mass, atname, Xeq, initcond_file, Tvib
      namelist /systemB/ &
         mass, atname, Xeq, initcond_file, Tvib

      sys%nat = nat
      allocate (mass(nat), atname(nat), Xeq(3*nat), &
                sys%mass(sys%nat), sys%atname(sys%nat), sys%Xeq(3*sys%nat))
      initcond_file = ""
      mass = 0._dp
      Xeq = 0._dp
      Tvib = 0.0_dp
      rewind (10)
      if (s == "A") then
         read (iunit, nml=systemA, iostat=ios)
      else
         read (iunit, nml=systemB, iostat=ios)
      end if
      if (ios .ne. 0) write (sal_unit, "('Namelist system',a1,' not found')") s

      do iat = 1, nat
         if (mass(iat) == 0._dp) then
            call get_mass(atname(iat), mass(iat))
            if (mass(iat) == -1.0_dp) stop
         end if
      end do

      if (s == "A") then
         if (ios == 0) write (sal_unit, nml=systemA)
      else
         if (ios == 0) write (sal_unit, nml=systemB)
      end if

      sys%nat = nat
      sys%mass = mass/autouma
      sys%atname = atname
      sys%Xeq = Xeq/autoA
      sys%initcond_file = initcond_file
      sys%Tvib = Tvib * kb / autoJ
      deallocate (mass, atname, Xeq)
   end subroutine
end module
