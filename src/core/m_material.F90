! Material database reader.  Decoupled from any mesh type: callers extract
! material IDs from their mesh and pass them in as a plain integer array.
module m_material
    use m_constants
    use m_types
    implicit none
    private
    public :: InitialiseMaterials

    ! Number of scattering-group columns in the MATS.txt file format.
    ! The C5G7 benchmark writes 7 columns regardless of the group count
    ! used in the simulation; we read all and store the first n_groups.
    integer, parameter :: N_SCAT_IN_FILE = 7

contains

    ! ------------------------------------------------------------------
    ! Top-level entry point.
    !
    !   mat_ids  – every material ID present in the mesh (duplicates OK;
    !              only unique IDs matching MATS.txt entries are loaded)
    !   n_groups – energy-group count used by this run
    ! ------------------------------------------------------------------
    subroutine InitialiseMaterials(mats, mat_ids, n_groups, filename, printout)
        type(t_material), allocatable, intent(out) :: mats(:)
        integer,          intent(in)  :: mat_ids(:)
        integer,          intent(in)  :: n_groups
        character(len=*), intent(in)  :: filename
        logical,          intent(in)  :: printout

        integer :: i, max_id

        if (size(mat_ids) == 0) return
        max_id = maxval(mat_ids)
        if (max_id < 1) return

        allocate(mats(max_id))
        call ParseMaterialDeck(mats, mat_ids, n_groups, filename)

        do i = 1, max_id
            if (allocated(mats(i)%SigA)) call UpdateComputables(mats(i), n_groups)
        end do

        if (printout) call PrintMaterialSummary(mats, n_groups)
    end subroutine InitialiseMaterials

    ! ------------------------------------------------------------------
    ! Parse the MATS.txt file and populate mats(:).
    ! Only IDs listed in requested_ids are loaded.
    ! ------------------------------------------------------------------
    subroutine ParseMaterialDeck(mats, requested_ids, n_groups, filename)
        type(t_material), intent(inout) :: mats(:)
        integer,          intent(in)    :: requested_ids(:), n_groups
        character(len=*), intent(in)    :: filename

        integer            :: u, id, g, ios, i, max_id, n_scat
        character(len=1024):: line
        real(dp)           :: temp_scat(N_SCAT_IN_FILE)
        logical, allocatable :: processed(:)

        open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(A)') " WARNING: Could not open material file: " // trim(filename)
            return
        end if

        max_id = size(mats)
        n_scat = min(n_groups, N_SCAT_IN_FILE)
        allocate(processed(max_id)); processed = .false.

        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)

            ! Skip headers and separators
            if (line == '' .or. line(1:1) == '-' .or. index(line, 'matID') > 0) cycle

            read(line, *, iostat=ios) id
            if (ios /= 0 .or. id < 1 .or. id > max_id) cycle
            if (.not. any(requested_ids == id) .or. processed(id)) cycle

            ! Allocate storage on first encounter of this ID
            if (.not. allocated(mats(id)%SigA)) then
                allocate(mats(id)%SigA(n_groups),   mats(id)%Nu(n_groups),     &
                         mats(id)%SigF(n_groups),   mats(id)%NuSigF(n_groups), &
                         mats(id)%Chi(n_groups),    mats(id)%Src(n_groups),    &
                         mats(id)%SigmaS(n_groups, n_groups),                  &
                         mats(id)%D(n_groups),      mats(id)%SigmaT(n_groups), &
                         mats(id)%SigmaR(n_groups))
                mats(id)%SigmaS = 0.0_dp
                mats(id)%name   = 'Unknown'
            end if

            do g = 1, n_groups
                if (g > 1) then
                    read(u, '(A)', iostat=ios) line
                    if (ios /= 0) then
                        write(*,'(A,I3,A,I3)') " ERROR: unexpected EOF reading material ", id, " group ", g
                        stop "Material parsing failed."
                    end if
                end if

                temp_scat = 0.0_dp
                if (g == 1) then
                    read(line, *, iostat=ios) id, mats(id)%SigA(g), mats(id)%Nu(g),   &
                                  mats(id)%SigF(g), mats(id)%Chi(g), mats(id)%Src(g), &
                                  (temp_scat(i), i=1, N_SCAT_IN_FILE), mats(id)%name
                else
                    read(line, *, iostat=ios) id, mats(id)%SigA(g), mats(id)%Nu(g),   &
                                  mats(id)%SigF(g), mats(id)%Chi(g), mats(id)%Src(g), &
                                  (temp_scat(i), i=1, N_SCAT_IN_FILE)
                end if

                if (ios /= 0) then
                    write(*,'(A,I3,A,I3)') " ERROR: parse failure for material ", id, " group ", g
                    write(*,'(A,A)') "   Line: ", trim(line)
                    stop "Material parsing failed."
                end if
                mats(id)%SigmaS(g, 1:n_scat) = temp_scat(1:n_scat)
            end do

            processed(id) = .true.
        end do
        close(u)
    end subroutine ParseMaterialDeck

    ! ------------------------------------------------------------------
    ! Derive computed cross sections from the raw parsed values.
    ! ------------------------------------------------------------------
    subroutine UpdateComputables(mat, n_groups)
        type(t_material), intent(inout) :: mat
        integer,          intent(in)    :: n_groups
        integer :: g

        if (sum(mat%Chi) > 1.0e-12_dp) mat%Chi = mat%Chi / sum(mat%Chi)

        do g = 1, n_groups
            mat%NuSigF(g) = mat%Nu(g) * mat%SigF(g)
            mat%SigmaT(g) = mat%SigA(g) + sum(mat%SigmaS(g, :))
            mat%SigmaR(g) = mat%SigmaT(g) - mat%SigmaS(g, g)

            if (mat%SigmaT(g) > 1.0e-10_dp) then
                mat%D(g) = 1.0_dp / (3.0_dp * mat%SigmaT(g))
            else
                write(*,'(A,I3)') " WARNING: SigmaT ~ 0 for group ", g
                mat%D(g) = 0.0_dp
            end if
        end do
    end subroutine UpdateComputables

    ! ------------------------------------------------------------------
    ! Pretty-print material summary table to stdout.
    ! ------------------------------------------------------------------
    subroutine PrintMaterialSummary(mats, n_groups)
        type(t_material), intent(in) :: mats(:)
        integer,          intent(in) :: n_groups
        integer :: i, g, gp, total_width
        character(len=256) :: fmt_data
        character(len=170) :: line, sep

        total_width = 61 + n_groups * 13

        line = " |" // repeat('=', total_width - 3) // "|"
        sep  = " |-----+------------+------------+------------+------------+" // &
               repeat('-', total_width - 61) // "|"

        write(fmt_data, '(A,I2,A)') '(A,I3,X,4(" | ",ES10.3)," | ",', n_groups, '(ES11.4," |"),A)'

        write(*, '(A)') line
        write(*, '(A)') " |" // repeat(' ', (total_width - 28) / 2) // &
                        "MATERIAL DATABASE SUMMARY" // &
                        repeat(' ', (total_width - 28) / 2) // "|"
        write(*, '(A)') line

        do i = 1, size(mats)
            if (.not. allocated(mats(i)%SigA)) cycle
            write(*, '(A,I3,A,A,A,T' // int2str(total_width) // ',A)') &
                " | [M] Material ID: ", i, " (", trim(mats(i)%name), ")", "|"
            write(*, '(A)') " |" // repeat('-', total_width - 3) // "|"
            write(*, '(A,T63,A,T' // int2str(total_width) // ',A)') &
                " | Grp |    SigT    |    SigA    |   NuSigF   |     Chi    |", &
                merge("Scattering Matrix (From Row \\ To Column)", &
                      "Scattering                               ", n_groups >= 5), "|"
            write(*, '(A)') sep
            do g = 1, n_groups
                write(*, fmt_data) " |", g, mats(i)%SigmaT(g), mats(i)%SigA(g), &
                    mats(i)%NuSigF(g), mats(i)%Chi(g), &
                    (mats(i)%SigmaS(g, gp), gp=1, n_groups), ""
            end do
            write(*, '(A)') line
        end do
    end subroutine PrintMaterialSummary

    pure function int2str(n) result(s)
        integer, intent(in) :: n
        character(len=4) :: s
        write(s, '(I4)') n
        s = adjustl(s)
    end function int2str

end module m_material
