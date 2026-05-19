#include <petsc/finclude/petscsys.h>

! MANIAC 2.0 unified driver.
! All three physics modes read ASMG mesh files.
! FEM meshes are ASMG files without knot vectors; the reader auto-generates
! trivial open knot vectors, making FEM a degenerate IGA case.
!
!   physics_mode = "transport_iga"   -- IGA DG-SN transport
!   physics_mode = "diffusion_iga"   -- IGA CG diffusion
!   physics_mode = "diffusion_fem"   -- FEM CG diffusion (ASMG with trivial knots)
program maniac
    use m_constants
    use m_types
    use m_types_iga
    use m_utilities
    use m_quadrature, only: LinearQuadrature, QuadrilateralQuadrature, &
                             HexahedralQuadrature, AngularQuadrature
    use m_material
    use m_asmg,                 only: read_asmg_mesh, write_mesh_to_files
    use m_basis_iga,            only: InitialiseBasis
    use m_sweep_order,          only: InitialiseGeometry
    use m_transport_precompute, only: InitialiseTransport
    use m_transport_iga,        only: SolveTransport
    use m_diffusion_iga,        only: SolveDiffusion
    use m_output_iga,           only: export_transport_vtk, export_diffusion_vtk
    use petscsys
    implicit none

    ! ---- Configuration (namelist defaults) ----------------------------
    character(len=64)  :: physics_mode   = "transport_iga"
    character(len=256) :: mesh_file      = "/home/tom/Documents/MANIAC/MANIAC_2.0/input/rod_test.asmg"
    character(len=256) :: mat_file       = "/home/tom/Documents/MANIAC/MANIAC_2.0/input/MATS.txt"
    character(len=256) :: output_dir     = "/home/tom/Documents/MANIAC/MANIAC_2.0/output"
    integer            :: sn_order       = 16
    integer            :: max_outer      = 900
    integer            :: solver_type    = SOLVER_KSP_GMRES
    integer            :: preconditioner = PRECON_GAMG
    integer            :: vtk_refine     = 4
    real(dp)           :: tol            = 1.0e-7_dp
    logical            :: is_eigenvalue  = .true.
    logical            :: is_adjoint     = .false.
    logical            :: ang_out        = .false.
    integer            :: ref_ids(8)     = 0
    integer            :: n_ref_ids      = 0

    namelist /maniac_config/ physics_mode, mesh_file, mat_file, output_dir, &
        sn_order, max_outer, solver_type, preconditioner, vtk_refine, &
        tol, is_eigenvalue, is_adjoint, ang_out, ref_ids, n_ref_ids

    ! ---- Common objects -----------------------------------------------
    type(t_material), allocatable :: mats(:)
    integer,          allocatable :: ref_id_list(:)

    ! ---- IGA/FEM objects (unified via ASMG reader) --------------------
    type(t_mesh_iga),  allocatable :: mesh
    type(t_finite_iga)             :: FE
    type(t_quadrature)             :: QuadVol, QuadFace
    type(t_sn_quadrature)          :: QuadSn

    ! ---- Solver outputs -----------------------------------------------
    type(t_transport_iga) :: TD
    real(dp), allocatable :: ang_flux(:,:,:)   ! transport: (n_dof, n_angles, n_groups)
    real(dp), allocatable :: scalar_flux(:,:)  ! transport: (n_dof, n_groups)
    real(dp), allocatable :: phi(:,:)          ! diffusion: (n_nodes, n_groups)
    real(dp)              :: k_eff
    integer,  allocatable :: sweep_order(:,:)

    ! ---- PETSc --------------------------------------------------------
    PetscErrorCode :: ierr

    ! ---- Misc ---------------------------------------------------------
    integer :: u_cfg, ios
    character(len=256) :: cfg_file, nametag, run_dir

    ! ------------------------------------------------------------------
    ! Read config
    ! ------------------------------------------------------------------
    cfg_file = "/home/tom/Documents/MANIAC/MANIAC_2.0/config.nml"
    if (command_argument_count() >= 1) call get_command_argument(1, cfg_file)

    open(newunit=u_cfg, file=trim(cfg_file), status='old', action='read', iostat=ios)
    if (ios == 0) then
        read(u_cfg, nml=maniac_config, iostat=ios)
        close(u_cfg)
        if (ios /= 0) write(*,'(A)') "Warning: error reading namelist; using defaults."
    else
        write(*,'(A)') "Config file not found; using defaults."
    end if

    call print_splash()
    write(*,'(A)') "  Physics: " // trim(physics_mode)
    write(*,'(A)') "  Mesh:    " // trim(mesh_file)
    write(*,'(A)') "============================================================"

    allocate(ref_id_list(max(1, n_ref_ids)))
    if (n_ref_ids > 0) then
        ref_id_list(1:n_ref_ids) = ref_ids(1:n_ref_ids)
    else
        ref_id_list = [0]
    end if

    call execute_command_line("mkdir -p " // trim(output_dir))

    nametag = derive_case_nametag(mesh_file)
    nametag = nametag(1:index(nametag,'.vtk')-1)

    ! ------------------------------------------------------------------
    ! Dispatch
    ! ------------------------------------------------------------------
    select case (trim(physics_mode))

    ! ------------------------------------------------------------------
    case ("transport_iga")

        allocate(mesh)
        call read_asmg_mesh(trim(mesh_file), mesh)

        run_dir = trim(output_dir) // "/transport/IGA/" // trim(nametag)
        call execute_command_line("mkdir -p " // trim(run_dir) // "/mesh_cache")
        call write_mesh_to_files(mesh, trim(run_dir) // "/mesh_cache")

        call InitialiseBasis(FE, mesh)

        QuadSn%order = sn_order
        if (mesh%dim == 2) then
            call LinearQuadrature      (QuadFace, 2*FE%order**mesh%dim + 1)
            call QuadrilateralQuadrature(QuadVol,  2*FE%order**mesh%dim + 1)
        else
            call QuadrilateralQuadrature(QuadFace, 2*FE%order**mesh%dim + 1)
            call HexahedralQuadrature  (QuadVol,  2*FE%order**mesh%dim + 1)
        end if
        call AngularQuadrature(mesh%dim, sn_order, QuadSn, is_adjoint)

        call InitialiseMaterials(mats, mesh%material_ids, mesh%n_groups, trim(mat_file), .true.)
        call InitialiseGeometry (mesh, FE, QuadSn, QuadFace, TD, sweep_order)
        call InitialiseTransport(mesh, FE, QuadSn, QuadVol, QuadFace, mats, TD)

        call SolveTransport(mesh, mats, FE, QuadSn, TD, &
                             scalar_flux, ang_flux, k_eff, &
                             sweep_order, ref_id_list, max_outer, tol, &
                             is_adjoint, is_eigenvalue)

        call export_transport_vtk(trim(run_dir), trim(nametag), &
                                   mesh, FE, QuadSn, scalar_flux, &
                                   mesh%n_groups, vtk_refine, ang_flux, ang_out)

    ! ------------------------------------------------------------------
    case ("diffusion_iga")

        PetscCall(PetscInitialize(PETSC_NULL_CHARACTER, ierr))

        allocate(mesh)
        call read_asmg_mesh(trim(mesh_file), mesh)

        run_dir = trim(output_dir) // "/diffusion/IGA/" // trim(nametag)
        call execute_command_line("mkdir -p " // trim(run_dir) // "/mesh_cache")
        call write_mesh_to_files(mesh, trim(run_dir) // "/mesh_cache")

        call InitialiseBasis(FE, mesh)

        if (mesh%dim == 2) then
            call LinearQuadrature      (QuadFace, 2*FE%order + 1)
            call QuadrilateralQuadrature(QuadVol,  2*FE%order + 1)
        else
            call QuadrilateralQuadrature(QuadFace, 2*FE%order + 1)
            call HexahedralQuadrature  (QuadVol,  2*FE%order + 1)
        end if

        call InitialiseMaterials(mats, mesh%material_ids, mesh%n_groups, trim(mat_file), .true.)

        call SolveDiffusion(mesh, FE, QuadVol, QuadFace, mats,        &
                             solver_type, preconditioner, ref_id_list, &
                             max_outer, tol, is_eigenvalue, is_adjoint, &
                             phi, k_eff)

        call export_diffusion_vtk(trim(run_dir), trim(nametag), &
                                   mesh, FE, phi, mesh%n_groups, vtk_refine)

        PetscCall(PetscFinalize(ierr))

    ! ------------------------------------------------------------------
    case default
        write(*,'(A)') "Unknown physics_mode: " // trim(physics_mode)
        stop 1

    end select

    write(*,'(A,F12.8)') "  k_eff = ", k_eff
    write(*,'(A)') ">>> Simulation Complete."
end program maniac
