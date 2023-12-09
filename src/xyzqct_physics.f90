module xyzqct_physics
   use xyzqct_constants, only: dp
contains

   subroutine get_COM(ndim, XP, nat1, nat2, mass, QCOM, PCOM)
      implicit none
      integer, intent(in) :: ndim, nat1, nat2
      real(dp), intent(in) :: XP(ndim), mass(nat2 - nat1 + 1)
      real(dp), intent(out) :: QCOM(3), PCOM(3)
      real(dp) :: mtot
      integer :: ix, iat

      QCOM = 0._dp
      PCOM = 0._dp
      mtot = 0._dp
      do iat = nat1, nat2
         mtot = mtot + mass(iat)
         do ix = 1, 3
            QCOM(ix) = QCOM(ix) + XP(3*(iat - 1) + ix)*mass(iat)
            PCOM(ix) = PCOM(ix) + XP(3*(iat - 1) + ix + ndim/2)
         end do
      end do
      QCOM = QCOM/mtot
   end subroutine get_COM

   subroutine get_LMOM_AMOM(ndim, XP, nat1, nat2, mass, QCOM, PCOM, LMOM, AMOM)
      implicit none
      integer, intent(in) :: ndim, nat1, nat2
      real(dp), intent(in) :: XP(ndim), mass(nat2 - nat1 + 1), QCOM(3), PCOM(3)
      real(dp), intent(out) :: LMOM(3), AMOM(3)
      real(dp) :: mtot, Q(3), P(3), prod(3)
      integer :: iat, ix

      mtot = 0._dp
      LMOM = 0._dp
      AMOM = 0._dp
      do iat = nat1, nat2
         mtot = mtot + mass(iat)
      end do

      do iat = nat1, nat2
         do ix = 1, 3
            Q(ix) = XP(3*(iat - 1) + ix) - QCOM(ix)
            P(ix) = &
               XP(ndim/2 + 3*(iat - 1) + ix) - PCOM(ix)*mass(iat)/mtot
            LMOM(ix) = LMOM(ix) + XP(ndim/2 + 3*(iat - 1) + ix)
         end do
         call cross_prod(Q, P, prod)
         AMOM = AMOM + prod
      end do
   end subroutine get_LMOM_AMOM

   subroutine cross_prod(v, w, p)
      implicit none
      real(dp), intent(in) :: v(3), w(3)
      real(dp), intent(out) :: p(3)

      p(1) = v(2)*w(3) - w(2)*v(3)
      p(2) = v(3)*w(1) - w(3)*v(1)
      p(3) = v(1)*w(2) - w(1)*v(2)

   end subroutine cross_prod

   subroutine get_inertia_moments(nat, pos, mass, inertia, inertia_mat)
      use xyzqct_lapack, only: eigh
      implicit none
      integer, intent(in) :: nat
      real(dp), intent(in) :: pos(3*nat), mass(nat)
      real(dp), intent(out) :: inertia(3), inertia_mat(3, 3) ! inertia_mat holds the eigenvector of the inertia matrix on out
      integer :: iat

      inertia_mat = 0._dp
      do iat = 1, nat
         inertia_mat(1, 1) = inertia_mat(1, 1) + mass(iat)*(pos(3*(iat - 1) + 2)**2 + pos(3*(iat - 1) + 3)**2)
         inertia_mat(2, 2) = inertia_mat(2, 2) + mass(iat)*(pos(3*(iat - 1) + 1)**2 + pos(3*(iat - 1) + 3)**2)
         inertia_mat(3, 3) = inertia_mat(3, 3) + mass(iat)*(pos(3*(iat - 1) + 1)**2 + pos(3*(iat - 1) + 2)**2)
         inertia_mat(1, 2) = inertia_mat(1, 2) - mass(iat)*pos(3*(iat - 1) + 1)*pos(3*(iat - 1) + 2)
         inertia_mat(1, 3) = inertia_mat(1, 3) - mass(iat)*pos(3*(iat - 1) + 1)*pos(3*(iat - 1) + 3)
         inertia_mat(2, 3) = inertia_mat(2, 3) - mass(iat)*pos(3*(iat - 1) + 2)*pos(3*(iat - 1) + 3)
      end do
      inertia_mat(2, 1) = inertia_mat(1, 2)
      inertia_mat(3, 1) = inertia_mat(1, 3)
      inertia_mat(3, 2) = inertia_mat(2, 3)
      call eigh(3, inertia_mat, inertia)
      inertia_mat = transpose(inertia_mat)
   end subroutine

   subroutine matrix_rotation(n, m, A, R)
      implicit none
      integer, intent(in) :: n, m
      real(dp), intent(inout) :: A(n*m)
      real(dp), intent(in) :: R(n, n)
      real(dp) :: aux(n, m), tmp(n, m)
      integer :: i, j

      do i = 1, m ! nat
         do j = 1, n ! ncoor
            aux(j, i) = A(3*(i - 1) + j)
         end do
      end do
      tmp = matmul(R, aux)
      do i = 1, m
         do j = 1, n
            A(3*(i - 1) + j) = tmp(j, i)
         end do
      end do
   end subroutine

   subroutine get_angular_velocity(inertia, AMOM, omega)
      implicit none
      real(dp), intent(in) :: inertia(3), AMOM(3)
      real(dp), intent(out) :: omega(3)
      real(dp) :: inertia_inv(3)

      inertia_inv = 0._dp
      where (inertia > 1e-10_dp)
         inertia_inv = 1._dp/inertia
      end where
      omega = AMOM*inertia_inv
   end subroutine

   subroutine add_angular_velocity(n, XP, mass, omega, sg)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(inout) :: XP(3*2*n)
      real(dp), intent(in) :: mass(n), omega(3), sg
      real(dp) :: x(3), rotvel(3)
      integer :: iat

      do iat = 1, n
         x = XP(3*(iat - 1) + 1:3*iat)
         call cross_prod(omega, x, rotvel)
         XP(3*n + 3*(iat - 1) + 1:3*n + 3*iat) = XP(3*n + 3*(iat - 1) + 1:3*n + 3*iat) + rotvel*mass(iat)*sg
      end do
   end subroutine

   subroutine rotate_euler(nat, XP, phi, theta, chi)
      implicit none
      integer, intent(in) :: nat
      real(dp), intent(inout) :: XP(3*2*nat)
      real(dp), intent(in) :: phi, theta, chi
      real(dp) :: EulerR(3, 3), c(3), s(3)
      c(1) = cos(phi)
      c(2) = cos(theta)
      c(3) = cos(chi)
      s(1) = sin(phi)
      s(2) = sin(theta)
      s(3) = sin(chi)

      EulerR(1, 1) = c(1)*c(2)*c(3) - s(1)*s(3)
      EulerR(1, 2) = s(1)*c(2)*c(3) + c(1)*s(3)
      EulerR(1, 3) = -s(2)*c(3)
      EulerR(2, 1) = -c(1)*c(2)*s(3) - s(1)*c(3)
      EulerR(2, 2) = c(1)*c(3) - c(2)*s(1)*s(3)
      EulerR(2, 3) = s(2)*s(3)
      EulerR(3, 1) = c(1)*s(2)
      EulerR(3, 2) = s(1)*s(2)
      EulerR(3, 3) = c(2)
      call matrix_rotation(3, nat, XP(:3*nat), EulerR)
      call matrix_rotation(3, nat, XP(3*nat + 1:), EulerR)
   end subroutine

   subroutine rotate_to_eckart(nat, x, xeq, m, Tbest)
      use xyzqct_lapack, only: eigh
      ! A. Y. Dymarsky and K. N. Kudin, J. Chem. Phys. 122, 124103 (2005)
      implicit none
      integer, intent(in) :: nat
      real(dp), intent(in) :: x(3*nat), xeq(3*nat), m(nat)
      real(dp), intent(out) :: Tbest(3, 3)
      real(dp) :: error, error_best, A(3, 3), A1(3, 3), A2(3, 3), v1(3), v2(3), w1(3, 3), xrot(3*nat), T(3, 3)
      integer :: iat, icoor, jcoor, iv, ii
      logical :: invert(8, 3)

      invert(1, :) = [.False., .False., .False.]
      invert(2, :) = [.True., .False., .False.]
      invert(3, :) = [.False., .True., .False.]
      invert(4, :) = [.False., .False., .True.]
      invert(5, :) = [.True., .True., .False.]
      invert(6, :) = [.True., .False., .True.]
      invert(7, :) = [.False., .True., .True.]
      invert(8, :) = [.True., .True., .True.]

      Tbest = 0._dp
      error_best = 1e10_dp

      A = 0._dp
      do iat = 1, nat
         do icoor = 1, 3
            do jcoor = 1, 3
               A(icoor, jcoor) = A(icoor, jcoor) + &
                                 m(iat)*x(3*(iat - 1) + icoor)*xeq(3*(iat - 1) + jcoor)
            end do
         end do
      end do
      A1 = matmul(A, transpose(A))
      A2 = matmul(transpose(A), A)

      call eigh(3, A1, v1)
      call eigh(3, A2, v2)
      call cross_prod(A1(:, 1), A1(:, 2), A1(:, 3))
      call cross_prod(A2(:, 1), A2(:, 2), A2(:, 3))

      do iv = 1, 8
         T = 0._dp
         error = 0._dp
         w1 = A1
         xrot = x
         do ii = 1, 3
            if (invert(iv, ii)) w1(:, ii) = -w1(:, ii)
         end do
         T = matmul(A2, transpose(w1))

         call matrix_rotation(3, nat, xrot, T)
         do ii = 1, nat
            error = error + &
                    sqrt(sum((xrot(3*(ii - 1) + 1:3*ii) - xeq(3*(ii - 1) + 1:3*ii))**2))
         end do
         if (error < error_best) then
            error_best = error
            Tbest = T
         end if
      end do
   end subroutine

   subroutine NM_analysis(n, x, mass, nfreqs, freqs, CXQ)
      use xyzqct_lapack, only: eigh
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: x(3*n), mass(n)
      real(dp), intent(out) :: freqs(3*n), CXQ(3*n, 3*n)
      integer, intent(out) :: nfreqs
      real(dp) :: sm(n)
      integer :: i, j

      freqs = 0._dp
      CXQ = 0._dp
      nfreqs = 0
      sm = sqrt(mass)

      call eval_hessian(3*n, x, CXQ)
      call filter_hessian(n, mass, x, CXQ)
      do i = 1, 3*n
         do j = 1, 3*n
            CXQ(i, j) = CXQ(i, j)/(sm((i - 1)/3 + 1)*sm((j - 1)/3 + 1))
         end do
      end do
      flush (6)
      call eigh(3*n, CXQ, freqs)

      do i = 1, 3*n
         if (freqs(i) > 1e-10_dp) nfreqs = nfreqs + 1
      end do
      freqs = sqrt(abs(freqs))
   end subroutine

   subroutine eval_hessian(n, x, hessian)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: hessian(n, n)

      integer :: i, j
      real(dp) :: pot, &
                  x1(n), &
                  grad1(n), grad2(n), &
                  deltag(n, n)

      real(dp), parameter :: delta = 1e-4_dp, fourdelta = 4*delta

      ! delta calculation taken from:
      !https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
      ! delta = (1e-r)**(1/3) with r the precision.

      hessian = 0._dp
      deltag = 0._dp
      do i = 1, n
         x1 = x
         x1(i) = x(i) + delta
         call potxyz(x1, pot, grad1)
         x1(i) = x(i) - delta
         call potxyz(x1, pot, grad2)
         deltag(i, :) = grad1 - grad2
         do j = 1, i
            hessian(i, j) = (deltag(j, i) + deltag(i, j))/fourdelta
            hessian(j, i) = hessian(i, j)
         end do
      end do
   end subroutine

   subroutine filter_hessian(nat, mass, xin, h)
      ! Tries to project translation and rotations out of the hessian.
      integer, intent(in) :: nat
      real(dp), intent(in) :: xin(3*nat), mass(nat)
      real(dp), intent(inout) :: h(3*nat, 3*nat)

      integer :: iat, i, j
      real(dp) :: x(3*nat), com(3), Iaxis(3, 3), &
                  tmp(3*nat, 3*nat), P(3*nat, 3*nat), dum(3), &
                  sm(nat)

      sm = sqrt(mass)
      ! Evaluate COM of the system.
      call get_COM(3*nat, xin, 1, nat, mass, com, dum)

      ! Set origin in COM
      do iat = 1, nat
         x(3*(iat - 1) + 1:3*iat) = xin(3*(iat - 1) + 1:3*iat) - com
      end do

      ! Evaluate inertia axis.
      call get_inertia_moments(nat, x, mass, dum, Iaxis)

      ! Generate projector out of translations and rotations.
      call trans_rot_projector(nat, x, mass, Iaxis, P)

      ! Mass scale hessian:
      do i = 1, 3*nat
         do j = i, 3*nat
            h(i, j) = h(i, j)/(sm((i - 1)/3 + 1)*sm((j - 1)/3 + 1))
            h(j, i) = h(i, j)
         end do
      end do

      ! Project translation and rotation out of hessian
      tmp = matmul(transpose(P), h)
      h = matmul(tmp, P)

      ! Mass scale hessian back:
      do i = 1, 3*nat
         do j = i, 3*nat
            h(i, j) = h(i, j)*(sm((i - 1)/3 + 1)*sm((j - 1)/3 + 1))
            h(j, i) = h(i, j)
         end do
      end do
   end subroutine

   subroutine trans_rot_projector(nat, x, mass, Iaxis, proj)
      ! Generates a basis of internal coordinates.
      ! Input:
      !   - x: Point in which the hessian was calculated.
      !   - Iaxis: Principal axis of inertia.

      ! Output:
      !  - proj: Projector to remove rotations and translations.
      integer, intent(in) :: nat
      real(dp), intent(in) :: x(3*nat), Iaxis(3, 3), mass(nat)
      real(dp), intent(out) :: proj(3*nat, 3*nat)

      integer :: iat, icoor, imax, i, j
      real(dp) :: norm, P(nat, 3), xyz(nat, 3), D(6, 3*nat), &
                  tmp(3*nat, 3*nat), sm(nat)
      real(dp), allocatable :: redD(:, :)

      real(dp), parameter :: thres = 1e-5_dp

      D = 0e0_dp
      xyz = 0e0_dp
      sm = sqrt(mass)

      do iat = 1, nat
         xyz(iat, :) = x(3*(iat - 1) + 1:3*iat)
      end do

      P = matmul(xyz, Iaxis)

      do iat = 1, nat
         ! Translational modes:
         D(1, 3*(iat - 1) + 1) = sm(iat)
         D(2, 3*(iat - 1) + 2) = sm(iat)
         D(3, 3*(iat - 1) + 3) = sm(iat)
         do icoor = 1, 3
            ! Rotational modes:
            D(4, 3*(iat - 1) + icoor) = (P(iat, 2)*Iaxis(icoor, 3) - P(iat, 3)*Iaxis(icoor, 2))*sm(iat)
            D(5, 3*(iat - 1) + icoor) = (P(iat, 3)*Iaxis(icoor, 1) - P(iat, 1)*Iaxis(icoor, 3))*sm(iat)
            D(6, 3*(iat - 1) + icoor) = (P(iat, 1)*Iaxis(icoor, 2) - P(iat, 2)*Iaxis(icoor, 1))*sm(iat)
         end do
      end do

      ! Normalize translational and rotational modes:
      imax = 0
      do icoor = 1, 6
         norm = sqrt(sum(D(icoor, :)**2))
         ! If norm of mode is too small discard it.
         if (norm > thres) then
            imax = imax + 1
            D(imax, :) = D(icoor, :)/norm
         end if
      end do

      ! Allocate reduced D matrix
      allocate (redD(imax, 3*nat))
      do i = 1, imax
         redD(i, :) = D(i, :)
      end do

      ! Generate projector to remove rotation and translation.
      tmp = matmul(transpose(redD), redD)
      proj = 0e0_dp
      do i = 1, 3*nat
         proj(i, i) = 1e0_dp
         do j = 1, 3*nat
            proj(i, j) = proj(i, j) - tmp(i, j)
         end do
      end do
   end subroutine
end module xyzqct_physics
