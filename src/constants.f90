module constants
      implicit none
      integer, parameter :: dp = kind(1d0), &
                            sal_unit = 90, &
                            xyz_unit = 91, &
                            end_unit = 92
      real(dp),parameter,public :: autofs=2.41888e-2_dp, &
                                   autoA=0.529177, &
                                   autokg=9.10939e-31_dp, &
                                   autouma=5.4858e-4_dp, &
                                   autocm_1=219474.625_dp, &
                                   pi=acos(-1._dp)
end module constants
