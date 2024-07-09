program test
    implicit none
    character(len=100) :: testname
    logical :: ok, passed
    ok = .true.
    testname = "NM_DDEABM_userdefined"
    call run_traj_sysA(testname, passed)
    ok = ok .and. passed

    testname = "NM_verlet_userdefined"
    call run_traj_sysA(testname, passed)
    ok = ok .and. passed

    testname = "NM_DDEABM_NM"
    call run_traj_sysA(testname, passed)
    ok = ok .and. passed

    testname = "NM_verlet_NM"
    call run_traj_sysA(testname, passed)
    ok = ok .and. passed

    testname = "NM_DDEABM_AS"
    call run_traj_sysA(testname, passed)
    ok = ok .and. passed

    testname = "userdefined_DDEABM_userdefined"
    call run_traj_sysA(testname, passed)
    ok = ok .and. passed

    if (ok) then
        write(*,*) " :) All test passed successfully :)"
    else
        write(*,*) ":( One or more tests failed :("
    end if
end program

subroutine run_traj_sysA(testname, passed)
    use testutils, only: all_close
    implicit none
    character(len=100), intent(in) :: testname
    character(len=100) :: cmd, outf
    logical, intent(out) :: passed
    real(8), allocatable :: XP0ref(:), XPref(:), XP0(:), XP(:)
    real(8) :: tref, t
    integer :: u, v, nat, d
    logical :: ok0, ok

    write(cmd, '(a,a,a)') "./bin/QCT.x < ./inp/", trim(testname), ".inp"
    write(outf, '(a,a,a)') "./refs/end_conditions_", trim(testname)
    call system("rm -rf sal end_conditions")
    call system(trim(cmd))

    open(newunit=u, file=trim(outf), status="old")
    open(newunit=v, file="end_conditions", status="old")
    read(u,*) nat
    read(v,*)
    allocate(XP0ref(2*3*nat), XPref(2*3*nat), XP0(2*3*nat), XP(2*3*nat))
    read(u,*) d, XP0ref, XPref, tref
    read(v,*) d, XP0, XP, t
    close(u)
    close(v)

    ok0 = all_close(XP0, XP0ref, 1d-6) 
    ok = all_close(XP, XPref, 1d-6) 
    passed = ok0 .and. ok
    write(*,'(a)') "----------------------------"
    write(*,'(a,a)') "test: ", trim(testname)
    write(*,'(a,l)') " pass: ", passed
    if (.not. passed) then
        write(*,'(a,l)') " initial conditions: ", ok0
        write(*,'(a,l)') " final conditions: ", ok
    end if
    write(*,'(a)') "----------------------------"
end subroutine
