module utils
    use constants, only: sal_unit
    implicit none

    contains

    subroutine code_starter()
        implicit none
        character(8)  :: date
        character(10) :: time
        character(5)  :: zone
        character(len=30) :: hostname
        integer,dimension(8) :: values

        call date_and_time(date,time,zone,values)
        call hostnm(hostname)

        call write_header()
        write(sal_unit, *) 
        write(sal_unit, *) "Hostname:", trim(hostname)
        write(sal_unit, "(X,A5,2X,2(I0.2,A1),I4,A1,3(I0.2,A1)/)") &
            "Date:", values(3),"/",values(2),"/",values(1),"-",values(5),":",values(6),":",values(7)
        write(sal_unit, *) 

    end subroutine

    subroutine write_header()
        implicit none
write(sal_unit,*) " ___    ___ ___    ___ ________                ________  ________ _________   "
write(sal_unit,*) "|\  \  /  /|\  \  /  /|\_____  \              |\   __  \|\   ____\\___   ___\ "
write(sal_unit,*) "\ \  \/  / | \  \/  / /\|___/  /| ____________\ \  \|\  \ \  \___\|___ \  \_| "
write(sal_unit,*) " \ \    / / \ \    / /     /  / /|\____________\ \  \\\  \ \  \       \ \  \  "
write(sal_unit,*) "  /     \/   \/  /  /     /  /_/_\|____________|\ \  \\\  \ \  \____   \ \  \ "
write(sal_unit,*) " /  /\   \ __/  / /      |\________\             \ \_____  \ \_______\  \ \__\"
write(sal_unit,*) "/__/ /\ __\\___/ /        \|_______|              \|___| \__\|_______|   \|__|"
write(sal_unit,*) "|__|/ \|__\|___|/                                       \|__|                 "
    end subroutine
end module
