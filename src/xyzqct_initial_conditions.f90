module xyzqct_initial_conditions
   use xyzqct_constants, only: dp, sal_unit, end_unit, cond_unitA, cond_unitB, autoeV, pi, kb, autoJ, autoA
   use xyzqct_settings, only: ndim, temperature, rfin, &
                              sysA, sysB, System
   use xyzqct_rand, only: ran2
   implicit none
   integer :: init_cond_mode, max_condA, max_condB
   real(dp) :: Ecoll, Trot, rini, bmax, bmin, bparam, Ttrans, A_capture, n_capture, inertiaA(3), inertiaB(3)

   private :: init_cond_mode, bparam

   namelist /collision/ &
      Ecoll, &
      Trot, &
      Ttrans, &
      rini, &
      bmax, &
      bmin, &
      A_capture, &
      n_capture
   ! init_cond_mode:
   ! 0 => Read from file
   ! 1 => NM initial condition
   ! 2 => NM initial condition (sampled at temperature T)

contains
   subroutine set_init_cond(mode)
      use xyzqct_physics, only: get_inertia_moments
      implicit none
      integer :: mode
      real(dp) :: inertia_vec(3, 3)
      init_cond_mode = mode
      write (sal_unit, "(/A)") "************************************"
      write (sal_unit, *) "SETTING INITIAL CONDITIONS..."
      select case (mode)
      case (0)
         write (sal_unit, "(/A)") "Reading initial conditions from file:", trim(sysA%initcond_file)
         open (cond_unitA, file=trim(sysA%initcond_file), status="old")
         read (cond_unitA, *) max_condA
         rewind (cond_unitA)
         write (sal_unit, *) "Total number of initial conditions =", max_condA
      case (1)
         write (sal_unit, "(/A)") "Using NM initial conditions"
         call setup_NM(sysA)
      case (2)
         write (sal_unit, "(/A)") "Using NM initial conditions (sampled at temperature T)"
         call setup_NM(sysA)
         call compute_Qmax(temperature, sysA)
      case (3)
         write (sal_unit, "(/A)") "Using initial conditions for A+B collision"
         Ecoll = 0._dp
         Trot = 0._dp
         Ttrans = 0._dp
         bmin = 0._dp
         A_capture = 0._dp
         n_capture = 0._dp
         rewind (10)
         read (10, nml=collision)
         write (sal_unit, nml=collision)
         Ecoll = Ecoll/autoeV
         Trot = Trot*kb/autoJ
         Ttrans = Ttrans*kb/autoJ
         if (Ecoll > 0._dp .and. Ttrans > 0._dp) then
            write (sal_unit, *) "ERROR: Either Ecoll or Ttrans must be zero."
            stop
         end if
         if (rini > rfin) then
            write (sal_unit, *) "ERROR: rini must be smaller than rfin"
            stop
         end if

         if (bmax > rini) then
            write (sal_unit, *) "ERROR: bmax must be smaller or equal than rfin"
            stop
         end if

         open (cond_unitA, file=trim(sysA%initcond_file), status="old")
         read (cond_unitA, *) max_condA
         rewind (cond_unitA)
         write (sal_unit, *) "Total number of initial conditions (A) =", max_condA
         open (cond_unitB, file=trim(sysB%initcond_file), status="old")
         read (cond_unitB, *) max_condB
         rewind (cond_unitB)
         write (sal_unit, *) "Total number of initial conditions (B) =", max_condB
         if (sysA%nat > 1 .and. Trot > 0._dp) then
            if (all(sysA%Xeq == 0._dp)) then
               write (sal_unit, *) "Unknown equilibrium position of system A. Cannot perform rotational analysis"
               stop
            end if
            write (sal_unit, "('Rotational analysis of system A:')")
            call get_inertia_moments(sysA%nat, sysA%Xeq, sysA%mass, inertiaA, inertia_vec)
            write (sal_unit, *) "B (cm-1):", 1._dp/(4*pi*inertiaA*autoA*1e-8_dp*137)
         end if
         if (sysB%nat > 1 .and. Trot > 0._dp) then
            if (all(sysB%Xeq == 0._dp)) then
               write (sal_unit, *) "Unknown equilibrium position of system B. Cannot perform rotational analysis"
               stop
            end if
            write (sal_unit, "('Rotational analysis of system B:')")
            call get_inertia_moments(sysB%nat, sysB%Xeq, sysB%mass, inertiaB, inertia_vec)
            write (sal_unit, *) "B (cm-1):", 1._dp/(4*pi*inertiaB*autoA*1e-8_dp*137)
         end if
      case default
         write (sal_unit, "(/A/)") "Unknown initial condition mode."
         stop
      end select
      ! write info to end_cond file
      select case (init_cond_mode)
      case (:2)
         write (end_unit, *) "itraj, XP0, XP, time"
      case (3)
         write (end_unit, *) "itraj, XP0, XP, time, bmax"
      end select
      write (sal_unit, "(/A/)") "************************************"
   end subroutine

   subroutine write_end_cond(itraj, time, XP0, XP)
      implicit none
      integer :: itraj
      real(dp), intent(in) :: time, XP0(ndim), XP(ndim)
      select case (init_cond_mode)
      case (:2)
         write (end_unit, *) itraj, XP0, XP, time
      case (3)
         write (end_unit, *) itraj, XP0, XP, time, bparam
      end select
   end subroutine

   subroutine get_init_cond(XP, propagate)
      implicit none
      real(dp), intent(out) :: XP(ndim)
      logical, intent(out) :: propagate
      propagate = .true.
      select case (init_cond_mode)
      case (0)
         call from_file_init_cond(max_condA, cond_unitA, sysA%nat, XP)
      case (1)
         call NM_init_cond(XP)
      case (2)
         call NM_init_cond_T(XP)
      case (3)
         call AplusB_init_cond(XP, propagate)
      end select
   end subroutine

   subroutine from_file_init_cond(max_cond, cond_unit, n, XP)
      implicit none
      integer :: max_cond, cond_unit, n
      real(dp), intent(out) :: XP(3*2*n)
      real(dp) :: r
      integer :: icond, i

      XP(:) = 0._dp

      call ran2(r)
      icond = floor(max_cond*r + 1)
      write (sal_unit, *) "Using icond =", icond

      read (cond_unit, *)
      do i = 1, icond
         read (cond_unit, *) XP
      end do
      rewind (cond_unit)
   end subroutine

   subroutine NM_init_cond(XP)
      use xyzqct_constants, only: pi
      use xyzqct_physics, only: get_COM, get_LMOM_AMOM, matrix_rotation, get_angular_velocity, &
                                add_angular_velocity, get_inertia_moments
      use xyzqct_hamiltonian, only: total_ener
      implicit none
      integer :: ix, ifreq
      real(dp), intent(out) :: XP(ndim)
      real(dp) :: phase(sysA%nfreqs), Q(sysA%nfreqs), P(sysA%nfreqs), mass, &
                  QCOM(3), PCOM(3), LMOM(3), AMOM(3), inertia(3), inertia_vec(3, 3), omega(3), &
                  kener, potener, E0, E

      sysA%amp = sqrt((2._dp*sysA%Qnum + 1._dp)/sysA%freqs) ! Maximum NM amplitudes

      call ran2(phase)
      XP = 0._dp
      phase = 2._dp*pi*phase
      Q = sysA%amp*sin(phase)
      P = sysA%freqs*sysA%amp*cos(phase)

      do ix = 1, ndim/2
         mass = sqrt(sysA%mass((ix - 1)/3 + 1))
         do ifreq = 1, sysA%nfreqs
            XP(ix) = XP(ix) + Q(ifreq)*sysA%CXQ(ix, ifreq)
            XP(ix + ndim/2) = XP(ix + ndim/2) + P(ifreq)*sysA%CXQ(ix, ifreq)
         end do
         XP(ix) = XP(ix)/mass
         XP(ndim/2 + ix) = XP(ndim/2 + ix)*mass
      end do
      XP(:ndim/2) = XP(:ndim/2) + sysA%Xeq

      call total_ener(0._dp, XP, kener, potener)
      E0 = kener + potener
      call get_inertia_moments(sysA%nat, XP(:3*sysA%nat), sysA%mass, inertia, inertia_vec)
      call matrix_rotation(3, sysA%nat, XP(:3*sysA%nat), inertia_vec)
      call matrix_rotation(3, sysA%nat, XP(3*sysA%nat + 1:), inertia_vec)
      call get_COM(ndim, XP, 1, sysA%nat, sysA%mass, QCOM, PCOM)
      call get_LMOM_AMOM(ndim, XP, 1, sysA%nat, sysA%mass, QCOM, PCOM, LMOM, AMOM)
      write (sal_unit, *) "Total angular of initial condition: ", sqrt(sum(AMOM**2))
      call get_angular_velocity(inertia, AMOM, omega)
      call add_angular_velocity(sysA%nat, XP, sysA%mass, omega, -1._dp)
      E = 1e10_dp
      do while (abs(E - E0)/E0 > 1e-4)
         XP(3*sysA%nat + 1:) = XP(3*sysA%nat + 1:)*E0/E
         call total_ener(0._dp, XP, kener, potener)
         E = (kener + potener)
      end do
      write (sal_unit, *) "Filtered angular momentum from initial condition ..."
      write (sal_unit, *) "Added vibrational energy to compensate rotational energy removal."
   end subroutine

   subroutine NM_init_cond_T(XP)
      use xyzqct_constants, only: h, c, kb, autocm_1
      use xyzqct_settings, only: temperature
      implicit none
      integer :: imode, ii
      real(dp), intent(out) :: XP(ndim)
      real(dp) :: rand1, rand2, factor, prob

      ! Sample the values of Qnum for temperature T
      do imode = 1, sysA%nfreqs
         factor = h*sysA%freqs(imode)*autocm_1*100*c/(kb*temperature)
         do ii = 1, 1000
            call random_number(rand1)
            call random_number(rand2)
            sysA%Qnum(imode) = floor((sysA%Qmax(imode) + 1)*rand1)
            prob = exp(-sysA%Qnum(imode)*factor)*(1.-exp(-factor))
            if (rand2 < prob) exit
         end do
      end do
      call NM_init_cond(XP)
   end subroutine

   subroutine read_Data4NM(read_unit, sys)
      use xyzqct_constants, only: autoA, autocm_1
      implicit none
      integer, intent(in) :: read_unit
      type(System), intent(inout) :: sys
      integer :: ix

      !read(read_unit, *)
      !read(read_unit, *)
      !do iat=1,nA
      !    read(read_unit, *) atnameA(iat), Xeq(3*(iat-1)+1:3*iat)
      !end do
      !Xeq = Xeq / autoA
      read (read_unit, *) sys%nfreqs
      allocate (sys%freqs(sys%nfreqs), &
                sys%CXQ(3*sys%nat, sys%nfreqs))
      read (read_unit, *) sys%freqs
      do ix = 1, 3*sys%nat
         read (read_unit, *) sys%CXQ(ix, :)
      end do
      sys%freqs = sys%freqs/autocm_1
   end subroutine read_Data4NM

   subroutine compute_Qmax(T, sys)
      use xyzqct_constants, only: h, c, kb, autocm_1
      implicit none
      real(dp), intent(in) :: T
      type(System), intent(inout) :: sys
      integer :: ifreq, iv
      real(dp) :: prob, factor
      integer, parameter :: maxv = 100
      real(dp), parameter :: threshold = 1.e-3_dp

      sys%Qmax = 0
      do ifreq = 1, sys%nfreqs
         factor = h*sys%freqs(ifreq)*autocm_1*100*c/(kb*T)
         do iv = 1, maxv
            prob = exp(-iv*factor)*(1.-exp(-factor))
            if (prob .lt. threshold) exit
         end do
         sys%Qmax(ifreq) = iv
      end do
      write (sal_unit, *) "Qmax =", sys%Qmax
   end subroutine

   subroutine AplusB_init_cond(XP, propagate)
      use xyzqct_physics, only: get_COM, get_inertia_moments, matrix_rotation, get_angular_velocity, &
                                add_angular_velocity, get_LMOM_AMOM, rotate_euler
      use xyzqct_settings, only: nat, mass
      implicit none
      real(dp), intent(out) :: XP(ndim)
      logical, intent(out) :: propagate
      real(dp) :: XPA(3*2*sysA%nat), XPB(3*2*sysB%nat), QCOM(3), PCOM(3), mtot, r, ang, &
                  inertia(3), inertia_vec(3, 3), LMOM(3), AMOM(3), omega(3), phi, theta, chi, mred, mA, mB, &
                  J(4), erel, erelmax, prob, Erot
      integer :: iat

      XP = 0.0_dp
      propagate = .true.
      call from_file_init_cond(max_condA, cond_unitA, sysA%nat, XPA)
      call from_file_init_cond(max_condB, cond_unitB, sysB%nat, XPB)

      mA = sum(sysA%mass)
      mB = sum(sysB%mass)
      mtot = mA + mB
      !-------------------------------
      ! Remove COM and momentum (A and B)
      call get_COM(3*2*sysA%nat, XPA, 1, sysA%nat, sysA%mass, QCOM, PCOM)
      do iat = 1, sysA%nat
         XPA(3*(iat - 1) + 1:3*iat) = XPA(3*(iat - 1) + 1:3*iat) - QCOM
         XPA(3*sysA%nat + 3*(iat - 1) + 1:3*sysA%nat + 3*iat) = &
            XPA(3*sysA%nat + 3*(iat - 1) + 1:3*sysA%nat + 3*iat) - PCOM*sysA%mass(iat)/mA
      end do
      call get_COM(3*2*sysA%nat, XPA, 1, sysA%nat, sysA%mass, QCOM, PCOM)

      call get_COM(3*2*sysB%nat, XPB, 1, sysB%nat, sysB%mass, QCOM, PCOM)
      do iat = 1, sysB%nat
         XPB(3*(iat - 1) + 1:3*iat) = XPB(3*(iat - 1) + 1:3*iat) - QCOM
         XPB(3*sysB%nat + 3*(iat - 1) + 1:3*sysB%nat + 3*iat) = &
            XPB(3*sysB%nat + 3*(iat - 1) + 1:3*sysB%nat + 3*iat) - PCOM*sysB%mass(iat)/mB
      end do
      call get_COM(3*2*sysB%nat, XPB, 1, sysB%nat, sysB%mass, QCOM, PCOM)
      !-------------------------------

      !-------------------------------
      ! Align with inertia axis and remove angular momentum + euler rotation
      ! A
      if (sysA%nat > 1) then
         call get_inertia_moments(sysA%nat, XPA(:3*sysA%nat), sysA%mass, inertia, inertia_vec)
         call matrix_rotation(3, sysA%nat, XPA(:3*sysA%nat), inertia_vec)
         call matrix_rotation(3, sysA%nat, XPA(3*sysA%nat + 1:), inertia_vec)
         call get_LMOM_AMOM(3*2*sysA%nat, XPA, 1, sysA%nat, sysA%mass, QCOM, PCOM, LMOM, AMOM)
         call get_angular_velocity(inertia, AMOM, omega)
         call add_angular_velocity(sysA%nat, XPA, sysA%mass, omega, -1._dp)
         if (Trot .gt. 0._dp) then
            call sample_J(inertiaA, Trot, J, Erot)
            write (sal_unit, *) "Setting Jx, Jy, Jz, J (A):", J
            write (sal_unit, *) "Rotational energy /au (A):", Erot
            call get_angular_velocity(inertia, J(1:3), omega)
            call add_angular_velocity(sysA%nat, XPA, sysA%mass, omega, 1._dp)
         end if

         call ran2(r)
         phi = 2*pi*r
         call ran2(r)
         theta = pi*r
         call ran2(r)
         chi = 2*pi*r
         call rotate_euler(sysA%nat, XPA, phi, theta, chi)
         write (sal_unit, *) "Setting euler phi, theta, chi (A) = ", phi, theta, chi
      end if

      ! B
      if (sysB%nat > 1) then
         call get_inertia_moments(sysB%nat, XPB(:3*sysB%nat), sysB%mass, inertia, inertia_vec)
         call matrix_rotation(3, sysB%nat, XPB(:3*sysB%nat), inertia_vec)
         call matrix_rotation(3, sysB%nat, XPB(3*sysB%nat + 1:), inertia_vec)
         call get_LMOM_AMOM(3*2*sysB%nat, XPB, 1, sysB%nat, sysB%mass, QCOM, PCOM, LMOM, AMOM)
         call get_angular_velocity(inertia, AMOM, omega)
         call add_angular_velocity(sysB%nat, XPB, sysB%mass, omega, -1._dp)
         if (Trot .gt. 0._dp) then
            call sample_J(inertiaB, Trot, J, Erot)
            write (sal_unit, *) "Setting Jx, Jy, Jz, J (B):", J
            write (sal_unit, *) "Rotational energy /au (B):", Erot
            call get_angular_velocity(inertia, J(1:3), omega)
            call add_angular_velocity(sysB%nat, XPB, sysB%mass, omega, 1._dp)
         end if

         call ran2(r)
         phi = 2*pi*r
         call ran2(r)
         theta = pi*r
         call ran2(r)
         chi = 2*pi*r
         call rotate_euler(sysB%nat, XPB, phi, theta, chi)
         write (sal_unit, *) "Setting euler phi, theta, chi (B) = ", phi, theta, chi
      end if
      !-------------------------------

      !-------------------------------
      ! Add pZ to system b
      if (Ecoll .ne. 0._dp) then
         erel = Ecoll
      else
         erelmax = 15._dp*Ttrans
         r = 1.0_dp
         prob = 0.0_dp
         do while (r > prob)
            call ran2(r)
            erel = erelmax*r
            prob = erel/Ttrans*exp(-erel/Ttrans)*exp(1.0_dp)
            call ran2(r)
         end do
      end if
      write (sal_unit, *) "Setting relative translational energy / au = ", erel
      mred = mA*mB/(mA + mB)
      PCOM = 0.0_dp
      PCOM(3) = sqrt(2._dp*mred*erel)
      do iat = 1, sysA%nat
         XPA(3*sysA%nat + 3*(iat - 1) + 1:3*sysA%nat + 3*iat) = &
            XPA(3*sysA%nat + 3*(iat - 1) + 1:3*sysA%nat + 3*iat) + PCOM*sysA%mass(iat)/mA
      end do
      do iat = 1, sysB%nat
         XPB(3*sysB%nat + 3*(iat - 1) + 1:3*sysB%nat + 3*iat) = &
            XPB(3*sysB%nat + 3*(iat - 1) + 1:3*sysB%nat + 3*iat) - PCOM*sysB%mass(iat)/mB
      end do
      !-------------------------------

      !-------------------------------
      ! Move system B to rini with bparam
      call ran2(r)
      bparam = bmin + (bmax - bmin)*sqrt(r)
      if (A_capture .ne. 0._dp) then
         propagate = bparam <= A_capture/erel**n_capture
         if (.not. propagate) then
            write (sal_unit, *) &
               "Trajectory will not propagate because bmax is too high according to capture model."
         end if
      end if
      write (sal_unit, *) "Setting bparam / au = ", bparam
      call ran2(r)
      ang = 2*pi*r
      QCOM = 0.0_dp
      QCOM(1) = bparam*cos(ang)
      QCOM(2) = bparam*sin(ang)
      QCOM(3) = sqrt(abs(rini**2 - bparam**2))
      do iat = 1, sysB%nat
         XPB(3*(iat - 1) + 1:3*iat) = XPB(3*(iat - 1) + 1:3*iat) + QCOM
      end do
      !-------------------------------

      XP(1:3*sysA%nat) = XPA(1:3*sysA%nat)
      XP(3*sysA%nat + 1:3*sysA%nat + 3*sysB%nat) = XPB(1:3*sysB%nat)
      XP(3*sysA%nat + 3*sysB%nat + 1:3*2*sysA%nat + 3*sysB%nat) = XPA(3*sysA%nat + 1:)
      XP(3*2*sysA%nat + 3*sysB%nat + 1:) = XPB(3*sysB%nat + 1:)

      ! Set COM at [0,0,0]
      call get_COM(ndim, XP, 1, nat, mass, QCOM, PCOM)
      mtot = sum(mass)
      do iat = 1, nat
         XP(3*(iat - 1) + 1:3*iat) = XP(3*(iat - 1) + 1:3*iat) - QCOM
      end do
      flush (sal_unit)
   end subroutine

   subroutine sample_J(inertia, T, J, Erot)
      ! Faraday Discuss. Chem. Soc., 1973,55, 93-99
      ! VENUS manual
      implicit none
      real(dp), intent(in) :: inertia(3), T
      real(dp), intent(out) :: J(4)
      integer :: nlin, ix
      real(dp) :: prob, r, diff1, diff2, Jmax, inertia_mean, Erot

      J = 0.0_dp
      nlin = 0
      diff1 = abs(inertia(1) - inertia(2))
      diff2 = abs(inertia(3) - inertia(2))

      if (diff1 > diff2) then
         ix = 1 ! sample Jx
      else
         ix = 3 ! sample Jz
      end if

      if (abs(inertia(1)) < 1e-10 .and. diff2 < 1e-10) nlin = 1 ! linear molecule

      if (nlin == 0) then
         Jmax = 25._dp*inertia(ix)*T
         r = 1.0_dp
         prob = 0.0_dp
         do while (r > prob)
            call ran2(r)
            J(ix) = real(floor(Jmax*r), dp)
            call ran2(r)
            if (r > 0.5_dp) J(ix) = -J(ix)
            prob = exp(-J(ix)**2/(2._dp*inertia(ix)*T))
            call ran2(r)
         end do
         Erot = J(ix)**2/inertia(ix)
      end if

      if (nlin == 1 .or. ix == 1) then
         inertia_mean = sqrt(inertia(2)*inertia(3))
         call ran2(r)
         J(4) = real(floor(sqrt(J(1)**2 - 2*inertia_mean*T*log(1._dp - r))), dp)
         call ran2(r)
         J(2) = sqrt(J(4)**2 - J(1)**2)*sin(2._dp*pi*r)
         J(3) = sqrt(J(4)**2 - J(1)**2)*cos(2._dp*pi*r)
         Erot = (J(2)**2/inertia(2) + J(3)**2/inertia(3) + Erot)/2._dp
      else
         inertia_mean = sqrt(inertia(1)*inertia(2))
         call ran2(r)
         J(4) = real(floor(sqrt(J(3)**2 - 2*inertia_mean*T*log(1._dp - r))), dp)
         call ran2(r)
         J(1) = sqrt(J(4)**2 - J(3)**2)*sin(2._dp*pi*r)
         J(2) = sqrt(J(4)**2 - J(3)**2)*cos(2._dp*pi*r)
         Erot = (J(1)**2/inertia(1) + J(2)**2/inertia(2) + Erot)/2._dp
      end if
   end subroutine

   subroutine setup_NM(sys)
      use xyzqct_physics, only: NM_analysis, get_COM, get_inertia_moments, matrix_rotation
      use xyzqct_utils, only: write_freq_NM, write_xyz
      implicit none
      type(system), intent(inout) :: sys
      integer :: nfreqs, ios, iat
      real(dp), allocatable :: freqs(:), CXQ(:, :)
      integer, allocatable :: Qnum(:)
      real(dp) :: com(3), dum(3), Iaxis(3, 3)
      namelist /Qvib/ Qnum

      if (all(sys%Xeq == 0._dp)) then
         write (sal_unit, *) "Equilibrium geometry of system not found. Please provide it in the system namelist"
         stop
      end if
      if (sys%initcond_file == "") then
         write (sal_unit, "(/A)") "Perfoming vibrational analysis of system:"
         write (sal_unit, *) "Moving system to COM and orient according to inertia moments"
         call get_COM(3*sys%nat, sys%Xeq, 1, sys%nat, sys%mass, com, dum)
         do iat = 1, sys%nat
            sys%Xeq(3*(iat - 1) + 1:3*iat) = sys%Xeq(3*(iat - 1) + 1:3*iat) - com
         end do
         call get_inertia_moments(sys%nat, sys%Xeq, sys%mass, dum, Iaxis)
         call matrix_rotation(3, sys%nat, sys%Xeq, Iaxis)
         call write_xyz(sal_unit, sys%nat, sys%Xeq, sys%atname, "Oriented system")
         allocate (freqs(3*sys%nat), CXQ(3*sys%nat, 3*sys%nat))
         call NM_analysis(sys%nat, sys%Xeq, sys%mass, nfreqs, freqs, CXQ)
         sys%nfreqs = nfreqs
         allocate (sys%freqs(nfreqs), sys%CXQ(3*sys%nat, nfreqs))
         sys%freqs = freqs(3*sys%nat - nfreqs + 1:)
         sys%CXQ = CXQ(:, 3*sys%nat - nfreqs + 1:)
         deallocate (freqs, CXQ)
         call write_freq_NM(sys%nat, sys%nfreqs, sys%freqs, sys%CXQ)
      else
         open (11, file=trim(sys%initcond_file), status="old")
         call read_Data4NM(11, sys)
         close (11)
      end if
      allocate (sys%Qnum(sys%nfreqs), &
                sys%amp(sys%nfreqs), &
                sys%Qmax(sys%nfreqs), &
                Qnum(sys%nfreqs))
      rewind (10)
      Qnum = 0
      read (10, nml=Qvib, iostat=ios)
      if (ios .ne. 0) write (sal_unit, *) "Namelist Qvib not found"
      write (sal_unit, nml=Qvib)
      sys%Qnum = Qnum
      deallocate (Qnum)
   end subroutine
end module xyzqct_initial_conditions
