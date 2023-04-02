module initial_conditions
    use constants, only: dp, sal_unit
    use settings, only : ndim
    implicit none
    character(len=80) :: initcond_file
    integer, parameter :: cond_unit = 100
    integer :: init_cond_mode, max_cond

    ! init_cond_mode:
    ! 0 => Read from file
    ! 1 => NM initial condition

    abstract interface
        subroutine get_init_cond_base(XP)
            use constants, only : dp 
            use settings, only : ndim
            real(dp), intent(out) :: XP(ndim)
        end subroutine get_init_cond_base
    end interface
    procedure(get_init_cond_base), pointer :: get_init_cond => null()

    contains
        subroutine set_init_cond(mode)
            implicit none
            integer :: mode
            open(cond_unit, file=trim(initcond_file), status="old")
            read(cond_unit,*) max_cond
            rewind(cond_unit)
            write(sal_unit,*) "Total number of initial conditions =", max_cond

            select case(mode)
                case(0)
                    write(sal_unit, "(/A)") "Reading initial conditions from file:", trim(initcond_file)
                    get_init_cond => from_file_init_cond
            end select 
        end subroutine

        subroutine from_file_init_cond(XP)
            use constants, only: dp, sal_unit
            implicit none
            real(dp), intent(out) :: XP(ndim)
            real(dp) :: r
            integer :: icond, i 

            XP(:) = 0._dp

            call RANDOM_NUMBER(r)
            icond = floor(max_cond * r + 1)
            write(sal_unit,*) "Using icond =", icond

            read(cond_unit, *)
            do i=1, icond
                read(cond_unit, *) XP
            end do
            rewind(cond_unit)
        end subroutine
end module initial_conditions
