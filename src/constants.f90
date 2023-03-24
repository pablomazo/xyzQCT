module constants
      implicit none
      integer, parameter :: dp = kind(1d0), &
                            sal_unit = 90, &
                            xyz_unit = 91, &
                            end_unit = 92
      real(dp),parameter,public :: autofs=2.41888e-2_dp
      real(dp),parameter,public :: autoA=0.529177
      real(dp),parameter,public :: autokg=9.10939e-31_dp
      real(dp),parameter,public :: autouma=5.4858e-4_dp
      real(dp),parameter,public :: autocm_1=219474.625_dp
end module constants
