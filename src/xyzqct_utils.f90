module xyzqct_utils
   use xyzqct_constants, only: sal_unit, dp
   implicit none

contains

   subroutine code_starter()
      implicit none
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      character(len=30) :: hostname
      integer, dimension(8) :: values

      call date_and_time(date, time, zone, values)
      call hostnm(hostname)

      call write_header()
      write (sal_unit, *)
      write (sal_unit, *) "Hostname: ", trim(hostname)
      write (sal_unit, "(X,A5,2X,2(I0.2,A1),I4,A1,3(I0.2,A1)/)") &
         "Date:", values(3), "/", values(2), "/", values(1), "-", values(5), ":", values(6), ":", values(7)
      write (sal_unit, *)

   end subroutine

   subroutine write_header()
      implicit none
      write (sal_unit, *) " ___    ___ ___    ___ ________                ________  ________ _________   "
      write (sal_unit, *) "|\  \  /  /|\  \  /  /|\_____  \              |\   __  \|\   ____\\___   ___\ "
      write (sal_unit, *) "\ \  \/  / | \  \/  / /\|___/  /| ____________\ \  \|\  \ \  \___\|___ \  \_| "
      write (sal_unit, *) " \ \    / / \ \    / /     /  / /|\____________\ \  \\\  \ \  \       \ \  \  "
      write (sal_unit, *) "  /     \/   \/  /  /     /  /_/_\|____________|\ \  \\\  \ \  \____   \ \  \ "
      write (sal_unit, *) " /  /\   \ __/  / /      |\________\             \ \_____  \ \_______\  \ \__\"
      write (sal_unit, *) "/__/ /\ __\\___/ /        \|_______|              \|___| \__\|_______|   \|__|"
      write (sal_unit, *) "|__|/ \|__\|___|/                                       \|__|                 "
   end subroutine

   subroutine write_freq_NM(nat, nfreq, freq, NM)
      use xyzqct_constants, only: autocm_1
      implicit none
      integer, intent(in) :: nfreq, nat
      real(dp), dimension(nfreq), intent(in) :: freq
      real(dp), dimension(3*nat, nfreq), intent(in) :: NM

      integer, parameter :: max_in_line = 5
      integer :: i, j, k, index, nlines, nprint, modul
      character(len=100) :: FORMAT1, COOR_NAME

      nlines = nfreq/max_in_line
      modul = mod(nfreq, max_in_line)
      if (modul .ne. 0) nlines = nlines + 1

      nprint = max_in_line
      write (sal_unit, "('ZPE:', F15.2, ' cm-1')") 0.5_dp*sum(freq)*autocm_1
      do i = 1, nlines
         ! This only affects the last line.
         if (i == nlines .and. modul .ne. 0) nprint = modul

         ! Write line with frequencies:
         write (FORMAT1, '(A5,I1,A10)') "(A17,", nprint, "(F15.2,X))"
         write (sal_unit, FORMAT1) 'Wavenumber/cm-1:', &
            freq(max_in_line*(i - 1) + 1:max_in_line*(i - 1) + nprint)*autocm_1

         ! Write normal mode coordinates:
         write (FORMAT1, '(A5,I1,A10)') "(A17,", nprint, "(F15.9,X))"
         do j = 1, nat
            do k = 1, 3
               index = 3*(j - 1) + k
               if (k == 3) then
                  write (COOR_NAME, '(A2,I3)') 'QZ', j
               elseif (k == 2) then
                  write (COOR_NAME, '(A2,I3)') 'QY', j
               else
                  write (COOR_NAME, '(A2,I3)') 'QX', j
               end if

               write (sal_unit, FORMAT1) &
                  trim(COOR_NAME), NM(index, max_in_line*(i - 1) + 1:max_in_line*(i - 1) + nprint)
            end do
         end do
         write (sal_unit, *)
      end do
   end subroutine

   subroutine write_xyz(u, nat, x, atname, message)
      use xyzqct_constants, only: dp, autoA
      implicit none
      integer, intent(in) :: nat, u
      real(dp), intent(in) :: x(3*nat)
      character(len=2), intent(in) :: atname(nat)
      character(len=*), intent(in) :: message
      integer :: i

      write (u, *) nat
      write (u, *) message
      do i = 1, nat
         write (u, '(A2,X,3(F17.12,X))') atname(i), x(3*(i - 1) + 1:3*i)*autoA
      end do
   end subroutine
end module
