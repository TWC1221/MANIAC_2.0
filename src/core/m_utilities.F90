module m_utilities
    implicit none

    public :: print_splash, derive_case_nametag, int_to_str
    contains


    pure function int_to_str(ii) result(res)
        integer, intent(in) :: ii
        character(len=12) :: res
        write(res, '(I0)') ii
        res = adjustl(res)
    end function int_to_str

    ! Replace the file extension with .vtk for output naming.
    pure function derive_case_nametag(filename) result(tag)
        character(len=*), intent(in) :: filename
        character(len=:), allocatable :: tag
        integer :: pos, dot_pos

        pos = index(filename, '/', back=.true.)
        if (pos == 0) pos = index(filename, '\', back=.true.)
        tag = filename(pos+1:)
        dot_pos = index(tag, '.', back=.true.)
        if (dot_pos > 0) then
            tag = tag(1:dot_pos-1) // ".vtk"
        else
            tag = tag // ".vtk"
        end if
    end function derive_case_nametag

    subroutine print_splash()
        implicit none
        integer :: i
        character(len=148), dimension(9) :: logo
        character(len=8)  :: date_str
        character(len=10) :: time_str

        logo(1) = "__/\\\\____________/\\\\___________/\\\\\\\\\___________/\\\\\_____/\\\________/\\\\\\\\\\\___________/\\\\\\\\\_________________/\\\\\\\\\_        "
        logo(2) = " _\/\\\\\\________/\\\\\\_________/\\\\\\\\\\\\\________\/\\\\\\___\/\\\_______\/////\\\///__________/\\\\\\\\\\\\\____________/\\\////////__       "
        logo(3) = "  _\/\\\//\\\____/\\\//\\\________/\\\/////////\\\_______\/\\\/\\\__\/\\\___________\/\\\____________/\\\/////////\\\_________/\\\/___________      "
        logo(4) = "   _\/\\\\///\\\/\\\/_\/\\\_______\/\\\_______\/\\\_______\/\\\//\\\_\/\\\___________\/\\\___________\/\\\_______\/\\\________/\\\_____________     "
        logo(5) = "    _\/\\\__\///\\\/___\/\\\_______\/\\\\\\\\\\\\\\\_______\/\\\\//\\\\/\\\___________\/\\\___________\/\\\\\\\\\\\\\\\_______\/\\\_____________    "
        logo(6) = "     _\/\\\____\///_____\/\\\_______\/\\\/////////\\\_______\/\\\_\//\\\/\\\___________\/\\\___________\/\\\/////////\\\_______\//\\\____________   "
        logo(7) = "      _\/\\\_____________\/\\\_______\/\\\_______\/\\\_______\/\\\__\//\\\\\\___________\/\\\___________\/\\\_______\/\\\________\///\\\__________  "
        logo(8) = "       _\/\\\_____________\/\\\__/\\\_\/\\\_______\/\\\__/\\\_\/\\\___\//\\\\\__/\\\__/\\\\\\\\\\\__/\\\_\/\\\_______\/\\\__/\\\____\////\\\\\\\\\_ "
        logo(9) = "        _\///______________\///__\///__\///________\///__\///__\///_____\/////__\///__\///////////__\///__\///________\///__\///________\/////////__"

        call date_and_time(date_str, time_str)

        write(*, *) ""
        do i = 1, 9
            write(*, '(A)') logo(i)
        end do

        write(*, *) ""
        write(*, '(A)') "                M U L T I G R O U P   A N G U L A R   N E U T R O N I C S : I S O G E O M E T R IC  A N A L Y T I C   C O M P U T A N T"
        write(*, '(A)') " ----------------------------------------------------------------------"
        
        ! --- SYSTEM METADATA ---
        write(*, '(A, A4, "-", A2, "-", A2, A, A2, ":", A2, ":", A2)') &
            " [ SYSTEM TIMESTAMP ] : ", date_str(1:4), date_str(5:6), date_str(7:8), &
            " @ ", time_str(1:2), time_str(3:4), time_str(5:6)
        
        write(*, '(A)') " [ OPERATOR TYPE    ] : DETERMINISTIC BOLTZMANN TRANSPORT"
        write(*, '(A)') " [ SPATIAL BASIS    ] : "
        write(*, '(A)') " [ ANGULAR BASIS    ] : DISCRETE ORDINATES (SN) QUADRATURE"
        write(*, '(A)') " [ ENERGY MODEL     ] : MULTIGROUP DISCRETIZATION"
        write(*, '(A)') " ----------------------------------------------------------------------"
        write(*, *) ""

    end subroutine print_splash
end module m_utilities
