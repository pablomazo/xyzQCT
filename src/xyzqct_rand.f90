MODULE xyzqct_rand
   use xyzqct_constants, only: dp, I4B, sal_unit
   IMPLICIT NONE
   INTEGER, PARAMETER :: K4B = selected_int_kind(9)
   INTEGER(K4B), PARAMETER :: hg = huge(1_K4B), hgm = -hg, hgng = hgm - 1
   INTEGER(K4B), SAVE :: lenran = 0, seq = 0
   INTEGER(K4B), SAVE :: iran0, jran0, kran0, nran0, mran0, rans
   INTEGER(K4B), DIMENSION(:, :), POINTER, SAVE :: ranseeds
   INTEGER(K4B), DIMENSION(:), POINTER, SAVE :: iran, jran, kran, &
                                                nran, mran, ranv
   INTEGER(I4B), PARAMETER :: NPAR_ARTH = 16, NPAR2_ARTH = 8
   REAL(dp), SAVE :: amm
   INTERFACE ran_hash
      MODULE PROCEDURE ran_hash_s, ran_hash_v
   END INTERFACE

   !----------------------------------------
   ! Add by P. Mazo to avoid use of nrutils
   INTERFACE reallocate
      MODULE PROCEDURE reallocate_iv, reallocate_im
   END INTERFACE

   INTERFACE arth
      MODULE PROCEDURE arth_i
   END INTERFACE

   INTERFACE ran2
      SUBROUTINE ran2_s(harvest)
         USE xyzqct_constants, only: dp
         REAL(dp), INTENT(OUT) :: harvest
      END SUBROUTINE ran2_s
      SUBROUTINE ran2_v(harvest)
         USE xyzqct_constants, only: dp
         REAL(dp), DIMENSION(:), INTENT(OUT) :: harvest
      END SUBROUTINE ran2_v
   END INTERFACE
   !----------------------------------------

CONTAINS
   !----------------------------------------
   ! Add by P. Mazo to avoid use of nrutils
   FUNCTION reallocate_iv(p, n)
      INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
      INTEGER(I4B), INTENT(IN) :: n
      INTEGER(I4B) :: nold, ierr
      allocate (reallocate_iv(n), stat=ierr)
      if (ierr /= 0) &
         write (sal_unit, *) 'reallocate_iv: problem in attempt to allocate memory'
      if (.not. associated(p)) RETURN
      nold = size(p)
      reallocate_iv(1:min(nold, n)) = p(1:min(nold, n))
      deallocate (p)
   END FUNCTION reallocate_iv
   FUNCTION reallocate_im(p, n, m)
      INTEGER(I4B), DIMENSION(:, :), POINTER :: p, reallocate_im
      INTEGER(I4B), INTENT(IN) :: n, m
      INTEGER(I4B) :: nold, mold, ierr
      allocate (reallocate_im(n, m), stat=ierr)
      if (ierr /= 0) &
         write (sal_unit, *) 'reallocate_im: problem in attempt to allocate memory'
      if (.not. associated(p)) RETURN
      nold = size(p, 1)
      mold = size(p, 2)
      reallocate_im(1:min(nold, n), 1:min(mold, m)) = &
         p(1:min(nold, n), 1:min(mold, m))
      deallocate (p)
   END FUNCTION reallocate_im

   FUNCTION arth_i(first, increment, n)
      INTEGER(I4B), INTENT(IN) :: first, increment, n
      INTEGER(I4B), DIMENSION(n) :: arth_i
      INTEGER(I4B) :: k, k2, temp
      if (n > 0) arth_i(1) = first
      if (n <= NPAR_ARTH) then
      do k = 2, n
         arth_i(k) = arth_i(k - 1) + increment
      end do
      else
      do k = 2, NPAR2_ARTH
         arth_i(k) = arth_i(k - 1) + increment
      end do
      temp = increment*NPAR2_ARTH
      k = NPAR2_ARTH
      do
         if (k >= n) exit
         k2 = k + k
         arth_i(k + 1:min(k2, n)) = temp + arth_i(1:min(k, n - k))
         temp = temp + temp
         k = k2
      end do
      end if
   END FUNCTION arth_i
   !----------------------------------------
   SUBROUTINE ran_init(length)
      IMPLICIT NONE
      INTEGER(K4B), INTENT(IN) :: length
      INTEGER(K4B) :: new, j, hgt
      if (length < lenran) RETURN
      hgt = hg
      if (hg /= 2147483647) write (sal_unit, *) 'ran_init: arith assump 1 fails'
      if (hgng >= 0) write (sal_unit, *) 'ran_init: arith assump 2 fails'
      if (hgt + 1 /= hgng) write (sal_unit, *) 'ran_init: arith assump 3 fails'
      if (not(hg) >= 0) write (sal_unit, *) 'ran_init: arith assump 4 fails'
      if (not(hgng) < 0) write (sal_unit, *) 'ran_init: arith assump 5 fails'
      if (hg + hgng >= 0) write (sal_unit, *) 'ran_init: arith assump 6 fails'
      if (not(-1_k4b) < 0) write (sal_unit, *) 'ran_init: arith assump 7 fails'
      if (not(0_k4b) >= 0) write (sal_unit, *) 'ran_init: arith assump 8 fails'
      if (not(1_k4b) >= 0) write (sal_unit, *) 'ran_init: arith assump 9 fails'
      if (lenran > 0) then
         ranseeds => reallocate(ranseeds, length, 5)
         ranv => reallocate(ranv, length - 1)
         new = lenran + 1
      else
         allocate (ranseeds(length, 5))
         allocate (ranv(length - 1))
         new = 1
         amm = nearest(1.0_dp, -1.0_dp)/hgng
         if (amm*hgng >= 1.0 .or. amm*hgng <= 0.0) &
            write (sal_unit, *) 'ran_init: arth assump 10 fails'
      end if
      ranseeds(new:, 1) = seq
      ranseeds(new:, 2:5) = spread(arth(new, 1, size(ranseeds(new:, 1))), 2, 4)
      do j = 1, 4
         call ran_hash(ranseeds(new:, j), ranseeds(new:, j + 1))
      end do
      where (ranseeds(new:, 1:3) < 0) &
         ranseeds(new:, 1:3) = not(ranseeds(new:, 1:3))
      where (ranseeds(new:, 4:5) == 0) ranseeds(new:, 4:5) = 1
      if (new == 1) then
         iran0 = ranseeds(1, 1)
         jran0 = ranseeds(1, 2)
         kran0 = ranseeds(1, 3)
         mran0 = ranseeds(1, 4)
         nran0 = ranseeds(1, 5)
         rans = nran0
      end if
      if (length > 1) then
         iran => ranseeds(2:, 1)
         jran => ranseeds(2:, 2)
         kran => ranseeds(2:, 3)
         mran => ranseeds(2:, 4)
         nran => ranseeds(2:, 5)
         ranv = nran
      end if
      lenran = length
   END SUBROUTINE ran_init
   SUBROUTINE ran_deallocate
   if (lenran > 0) then
      deallocate (ranseeds, ranv)
      nullify (ranseeds, ranv, iran, jran, kran, mran, nran)
      lenran = 0
   end if
   END SUBROUTINE ran_deallocate

   SUBROUTINE ran_seed(sequence, size, put, get)
      IMPLICIT NONE
      INTEGER, OPTIONAL, INTENT(IN) :: sequence
      INTEGER, OPTIONAL, INTENT(OUT) :: size
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: put
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: get
      if (present(size)) then
         size = 5*lenran
      else if (present(put)) then
         if (lenran == 0) RETURN
         ranseeds = reshape(put, shape(ranseeds))
         where (ranseeds(:, 1:3) < 0) ranseeds(:, 1:3) = not(ranseeds(:, 1:3))
         where (ranseeds(:, 4:5) == 0) ranseeds(:, 4:5) = 1
         iran0 = ranseeds(1, 1)
         jran0 = ranseeds(1, 2)
         kran0 = ranseeds(1, 3)
         mran0 = ranseeds(1, 4)
         nran0 = ranseeds(1, 5)
      else if (present(get)) then
         if (lenran == 0) RETURN
         ranseeds(1, 1:5) = (/iran0, jran0, kran0, mran0, nran0/)
         get = reshape(ranseeds, shape(get))
      else if (present(sequence)) then
         call ran_deallocate
         seq = sequence
      end if
   END SUBROUTINE ran_seed

   SUBROUTINE ran_hash_s(il, ir)
      IMPLICIT NONE
      INTEGER(K4B), INTENT(INOUT) :: il, ir
      INTEGER(K4B) :: is, j
      do j = 1, 4
         is = ir
         ir = ieor(ir, ishft(ir, 5)) + 1422217823
         ir = ieor(ir, ishft(ir, -16)) + 1842055030
         ir = ieor(ir, ishft(ir, 9)) + 80567781
         ir = ieor(il, ir)
         il = is
      end do
   END SUBROUTINE ran_hash_s

   SUBROUTINE ran_hash_v(il, ir)
      IMPLICIT NONE
      INTEGER(K4B), DIMENSION(:), INTENT(INOUT) :: il, ir
      INTEGER(K4B), DIMENSION(size(il)) :: is
      INTEGER(K4B) :: j
      do j = 1, 4
         is = ir
         ir = ieor(ir, ishft(ir, 5)) + 1422217823
         ir = ieor(ir, ishft(ir, -16)) + 1842055030
         ir = ieor(ir, ishft(ir, 9)) + 80567781
         ir = ieor(il, ir)
         il = is
      end do
   END SUBROUTINE ran_hash_v
END MODULE xyzqct_rand

SUBROUTINE ran2_s(harvest)
   USE xyzqct_constants, only: dp
   USE xyzqct_rand, ONLY: K4B, amm, lenran, ran_init, &
                          iran0, jran0, kran0, nran0, mran0, rans
   IMPLICIT NONE
   REAL(dp), INTENT(OUT) :: harvest
   if (lenran < 1) call ran_init(1)
   rans = iran0 - kran0
   if (rans < 0) rans = rans + 2147483579_k4b
   iran0 = jran0
   jran0 = kran0
   kran0 = rans
   nran0 = ieor(nran0, ishft(nran0, 13))
   nran0 = ieor(nran0, ishft(nran0, -17))
   nran0 = ieor(nran0, ishft(nran0, 5))
   rans = iand(mran0, 65535)
   mran0 = ishft(3533*ishft(mran0, -16) + rans, 16) + &
           3533*rans + 820265819_k4b
   rans = ieor(nran0, kran0) + mran0
   harvest = amm*merge(rans, not(rans), rans < 0)
END SUBROUTINE ran2_s

SUBROUTINE ran2_v(harvest)
   USE xyzqct_constants, only: dp
   USE xyzqct_rand, ONLY: K4B, amm, lenran, ran_init, &
                          iran, jran, kran, nran, mran, ranv
   IMPLICIT NONE
   REAL(dp), DIMENSION(:), INTENT(OUT) :: harvest
   INTEGER(K4B) :: n
   n = size(harvest)
   if (lenran < n + 1) call ran_init(n + 1)
   ranv(1:n) = iran(1:n) - kran(1:n)
   where (ranv(1:n) < 0) ranv(1:n) = ranv(1:n) + 2147483579_k4b
   iran(1:n) = jran(1:n)
   jran(1:n) = kran(1:n)
   kran(1:n) = ranv(1:n)
   nran(1:n) = ieor(nran(1:n), ishft(nran(1:n), 13))
   nran(1:n) = ieor(nran(1:n), ishft(nran(1:n), -17))
   nran(1:n) = ieor(nran(1:n), ishft(nran(1:n), 5))
   ranv(1:n) = iand(mran(1:n), 65535)
   mran(1:n) = ishft(3533*ishft(mran(1:n), -16) + ranv(1:n), 16) + &
               3533*ranv(1:n) + 820265819_k4b
   ranv(1:n) = ieor(nran(1:n), kran(1:n)) + mran(1:n)
   harvest = amm*merge(ranv(1:n), not(ranv(1:n)), ranv(1:n) < 0)
END SUBROUTINE ran2_v
