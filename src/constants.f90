module constants
      implicit none
      integer, parameter :: dp = kind(1d0), &
                            sal_unit = 90, &
                            xyz_unit = 91, &
                            end_unit = 92, &
                            cond_unit = 100
      real(dp),parameter,public :: autofs=2.41888e-2_dp, &
                                   autoA=0.529177_dp, &
                                   autokg=9.10939e-31_dp, &
                                   autouma=5.4858e-4_dp, &
                                   autocm_1=219474.625_dp, &
                                   pi=acos(-1._dp)
      contains
          subroutine get_mass(atname, mass)
              implicit none
              character(len=2), intent(in) :: atname
              real(dp), intent(out) :: mass

              select case (atname)
                case('H')
                    mass = 1.0080_dp
                case('He')
                    mass = 4.00260_dp
                case('Li')
                    mass = 7.0_dp
                case('Be')
                    mass = 9.012183_dp
                case('B')
                    mass = 10.81_dp
                case('C')
                    mass = 12.011_dp
                case('N')
                    mass = 14.007_dp
                case('O')
                    mass = 15.999_dp
                case('F')
                    mass = 18.99840316_dp
                case('Ne')
                    mass = 20.180_dp
                case('Na')
                    mass = 22.9897693_dp
                case('Mg')
                    mass = 24.305_dp
                case('Al')
                    mass = 26.981538_dp
                case('Si')
                    mass = 28.085_dp
                case('P')
                    mass = 30.97376200_dp
                case('S')
                    mass = 32.07_dp
                case('Cl')
                    mass = 35.45_dp
                case('Ar')
                    mass = 39.9_dp
                case('K')
                    mass = 39.0983_dp
                case('Ca')
                    mass = 40.08_dp
                case('Sc')
                    mass = 44.95591_dp
                case('Ti')
                    mass = 47.867_dp
                case('V')
                    mass = 50.9415_dp
                case('Cr')
                    mass = 51.996_dp
                case('Mn')
                    mass = 54.93804_dp
                case('Fe')
                    mass = 55.84_dp
                case('Co')
                    mass = 58.93319_dp
                case('Ni')
                    mass = 58.693_dp
                case('Cu')
                    mass = 63.55_dp
                case('Zn')
                    mass = 65.4_dp
                case('Ga')
                    mass = 69.723_dp
                case('Ge')
                    mass = 72.63_dp
                case('As')
                    mass = 74.92159_dp
                case('Se')
                    mass = 78.97_dp
                case('Br')
                    mass = 79.90_dp
                case('Kr')
                    mass = 83.80_dp
                case('Rb')
                    mass = 85.468_dp
                case('Sr')
                    mass = 87.62_dp
                case('Y')
                    mass = 88.90584_dp
                case('Zr')
                    mass = 91.22_dp
                case('Nb')
                    mass = 92.90637_dp
                case('Mo')
                    mass = 95.95_dp
                case('Tc')
                    mass = 96.90636_dp
                case('Ru')
                    mass = 101.1_dp
                case('Rh')
                    mass = 102.9055_dp
                case('Pd')
                    mass = 106.42_dp
                case('Ag')
                    mass = 107.868_dp
                case('Cd')
                    mass = 112.41_dp
                case('In')
                    mass = 114.818_dp
                case('Sn')
                    mass = 118.71_dp
                case('Sb')
                    mass = 121.760_dp
                case('Te')
                    mass = 127.6_dp
                case('I')
                    mass = 126.9045_dp
                case('Xe')
                    mass = 131.29_dp
                case('Cs')
                    mass = 132.9054520_dp
                case('Ba')
                    mass = 137.33_dp
                case('La')
                    mass = 138.9055_dp
                case('Ce')
                    mass = 140.116_dp
                case('Pr')
                    mass = 140.90766_dp
                case('Nd')
                    mass = 144.24_dp
                case('Pm')
                    mass = 144.91276_dp
                case('Sm')
                    mass = 150.4_dp
                case('Eu')
                    mass = 151.964_dp
                case('Gd')
                    mass = 157.2_dp
                case('Tb')
                    mass = 158.92535_dp
                case('Dy')
                    mass = 162.500_dp
                case('Ho')
                    mass = 164.93033_dp
                case('Er')
                    mass = 167.26_dp
                case('Tm')
                    mass = 168.93422_dp
                case('Yb')
                    mass = 173.05_dp
                case('Lu')
                    mass = 174.9668_dp
                case('Hf')
                    mass = 178.49_dp
                case('Ta')
                    mass = 180.9479_dp
                case('W')
                    mass = 183.84_dp
                case('Re')
                    mass = 186.207_dp
                case('Os')
                    mass = 190.2_dp
                case('Ir')
                    mass = 192.22_dp
                case('Pt')
                    mass = 195.08_dp
                case('Au')
                    mass = 196.96657_dp
                case('Hg')
                    mass = 200.59_dp
                case('Tl')
                    mass = 204.383_dp
                case('Pb')
                    mass = 207_dp
                case('Bi')
                    mass = 208.98040_dp
                case('Po')
                    mass = 208.98243_dp
                case('At')
                    mass = 209.98715_dp
                case('Rn')
                    mass = 222.01758_dp
                case('Fr')
                    mass = 223.01973_dp
                case('Ra')
                    mass = 226.02541_dp
                case('Ac')
                    mass = 227.02775_dp
                case('Th')
                    mass = 232.038_dp
                case('Pa')
                    mass = 231.03588_dp
                case('U')
                    mass = 238.0289_dp
                case('Np')
                    mass = 237.048172_dp
                case('Pu')
                    mass = 244.06420_dp
                case('Am')
                    mass = 243.061380_dp
                case('Cm')
                    mass = 247.07035_dp
                case('Bk')
                    mass = 247.07031_dp
                case('Cf')
                    mass = 251.07959_dp
                case('Es')
                    mass = 252.0830_dp
                case('Fm')
                    mass = 257.09511_dp
                case('Md')
                    mass = 258.09843_dp
                case('No')
                    mass = 259.10100_dp
                case('Lr')
                    mass = 266.120_dp
                case('Rf')
                    mass = 267.122_dp
                case('Db')
                    mass = 268.126_dp
                case('Sg')
                    mass = 269.128_dp
                case('Bh')
                    mass = 270.133_dp
                case('Hs')
                    mass = 269.1336_dp
                case('Mt')
                    mass = 277.154_dp
                case('Ds')
                    mass = 282.166_dp
                case('Rg')
                    mass = 282.169_dp
                case('Cn')
                    mass = 286.179_dp
                case('Nh')
                    mass = 286.182_dp
                case('Fl')
                    mass = 290.192_dp
                case('Mc')
                    mass = 290.196_dp
                case('Lv')
                    mass = 293.205_dp
                case('Ts')
                    mass = 294.211_dp
                case('Og')
                    mass = 295.216_dp
                case default
                    write(sal_unit, ("(/A)")) '*******************************'
                    write(sal_unit, *) 'Error: Unknown mass for symbol ', atname
                    write(sal_unit, ("(A/)")) '*******************************'
                    mass = -1.0_dp
              end select
          end subroutine
end module constants
