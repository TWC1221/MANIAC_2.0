#include <petsc/finclude/petscsys.h>
! MANIAC 2.0 unified driver.
!
! 'solver'  is set in config.nml: "transport" or "diffusion"
! 'method'  is detected automatically from the Problem_type field in the
!           ASMG mesh file: "iga" (contains "IGA") or "fem" (contains "FEM")
!
! Dispatch structure:
!   1. Detect method, read mesh, check quality
!   2. Build quadrature (n_pts depends on solver)
!   3. Build basis (IGA: NURBS, FEM: Lagrange)
!   4. PETSc init (diffusion only)
!   5. Load materials
!   6. Solve + output  (nested select: solver x method)
!   7. PETSc finalize (diffusion only)

program maniac
    use m_constants
    use m_types
    use m_types_iga
    use m_utilities
    use m_quadrature, only: LinearQuadrature, QuadrilateralQuadrature, &
                             HexahedralQuadrature, AngularQuadrature
    use m_material
    use m_asmg,                     only: read_asmg_mesh, write_mesh_to_files, detect_mesh_type
    use m_basis_iga,                only: InitialiseNurbsBasis
    use m_types_fem
    use m_basis_fem,                only: InitialiseLagrangeBasis
    use m_sweep_order,              only: InitialiseGeometry
    use m_transport_precompute_iga, only: InitialiseTransport
    use m_transport_precompute_fem, only: InitialiseTransport_FEM, InitialiseGeometry_FEM
    use m_transport_iga,            only: SolveTransport
    use m_transport_fem,            only: SolveTransport_FEM
    use m_diffusion_iga,            only: SolveDiffusion
    use m_diffusion_fem,            only: SolveDiffusion_FEM
    use m_output_iga,               only: export_transport_vtk_iga, export_diffusion_vtk_iga
    use m_output_fem,               only: export_transport_vtk_fem, export_diffusion_vtk_fem
    use petscsys
    implicit none

    ! ---- Configuration (namelist defaults) ----------------------------
    character(len=32)  :: solver         = "transport"   ! "transport" | "diffusion"
    character(len=256) :: mesh_file      = "input/rod_test.asmg"
    character(len=256) :: mat_file       = "input/MATS.txt"
    character(len=256) :: output_dir     = "output"
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

    namelist /maniac_config/ solver, mesh_file, mat_file, output_dir, &
        sn_order, max_outer, solver_type, preconditioner, vtk_refine, &
        tol, is_eigenvalue, is_adjoint, ang_out, ref_ids, n_ref_ids

    ! ---- Common objects -----------------------------------------------
    type(t_material), allocatable :: mats(:)
    integer,          allocatable :: ref_id_list(:)

    ! ---- Mesh objects — target so a class(t_mesh) pointer can alias them
    type(t_mesh_iga), allocatable, target :: mesh_iga
    type(t_mesh_fem), allocatable, target :: mesh_fem
    class(t_mesh),    pointer             :: mesh_ptr
    type(t_finite_iga)                    :: FE
    type(t_finite_fem)                    :: FE_fem
    type(t_quadrature)                    :: QuadVol, QuadFace
    type(t_sn_quadrature)                 :: QuadSn

    ! ---- Solver outputs -----------------------------------------------
    type(t_transport_data) :: TD
    real(dp), allocatable :: ang_flux(:,:,:)   ! (n_dof, n_angles, n_groups)
    real(dp), allocatable :: scalar_flux(:,:)  ! (n_dof, n_groups)
    real(dp), allocatable :: phi(:,:)          ! diffusion: (n_nodes, n_groups)
    real(dp)              :: k_eff
    integer,  allocatable :: sweep_order(:,:)

    ! ---- PETSc --------------------------------------------------------
    PetscErrorCode :: ierr = 0

    ! ---- Misc ---------------------------------------------------------
    integer            :: u_cfg, ios, dim_mesh, order_mesh, n_quad
    character(len=256) :: cfg_file, nametag, run_dir
    character(len=8)   :: method         ! auto-detected: "iga" or "fem"
    character(len=3)   :: method_upper   ! "IGA" or "FEM"

    nullify(mesh_ptr)

    ! ------------------------------------------------------------------
    ! Read config namelist
    ! ------------------------------------------------------------------
    cfg_file = "config.nml"
    if (command_argument_count() >= 1) call get_command_argument(1, cfg_file)

    open(newunit=u_cfg, file=trim(cfg_file), status='old', action='read', iostat=ios)
    if (ios == 0) then
        read(u_cfg, nml=maniac_config, iostat=ios)
        close(u_cfg)
        if (ios /= 0) write(*,'(A)') "Warning: error reading namelist; using defaults."
    else
        write(*,'(A)') "Config file not found; using defaults."
    end if

    ! Detect IGA vs FEM from the mesh file Problem_type header
    method = detect_mesh_type(trim(mesh_file))
    method_upper = merge("IGA", "FEM", trim(method) == "iga")

    call print_splash(solver, method, mesh_file, mat_file, output_dir, &
                      max_outer, tol, is_eigenvalue, is_adjoint, &
                      sn_order, solver_type, preconditioner)

    allocate(ref_id_list(max(1, n_ref_ids)))
    if (n_ref_ids > 0) then
        ref_id_list(1:n_ref_ids) = ref_ids(1:n_ref_ids)
    else
        ref_id_list = [0]
    end if

    call execute_command_line("mkdir -p " // trim(output_dir))

    ! ------------------------------------------------------------------
    ! Read mesh; alias via class pointer for shared downstream steps
    ! ------------------------------------------------------------------
    select case (trim(method))
    case ("iga")
        allocate(mesh_iga)
        call read_asmg_mesh(trim(mesh_file), mesh_iga)
        mesh_ptr => mesh_iga
    case ("fem")
        allocate(mesh_fem)
        call read_asmg_mesh(trim(mesh_file), mesh_fem)
        mesh_ptr => mesh_fem
    case default
        write(*,'(A)') "Unknown Problem_type in mesh file; expected IGA or FEM."
        stop 1
    end select

    call check_mesh_quality(mesh_ptr)
    dim_mesh   = mesh_ptr%dim
    order_mesh = mesh_ptr%order

    ! ------------------------------------------------------------------
    ! Quadrature  (2*p+1 Gauss points integrates degree-4p exactly)
    ! ------------------------------------------------------------------
    n_quad = 2*order_mesh + 1
    if (dim_mesh == 2) then
        call LinearQuadrature      (QuadFace, n_quad)
        call QuadrilateralQuadrature(QuadVol,  n_quad)
    else
        call QuadrilateralQuadrature(QuadFace, n_quad)
        call HexahedralQuadrature  (QuadVol,  n_quad)
    end if
    if (trim(solver) == "transport") then
        QuadSn%order = sn_order
        call AngularQuadrature(dim_mesh, sn_order, QuadSn, is_adjoint)
    end if

    ! ------------------------------------------------------------------
    ! Basis
    ! ------------------------------------------------------------------
    select case (trim(method))
    case ("iga")
        call InitialiseNurbsBasis(FE, mesh_iga)
    case ("fem")
        call InitialiseLagrangeBasis(FE_fem, dim_mesh, order_mesh, QuadVol, QuadFace)
    end select

    ! ------------------------------------------------------------------
    ! PETSc (diffusion only)
    ! ------------------------------------------------------------------
    if (trim(solver) == "diffusion") then
        PetscCall(PetscInitialize(PETSC_NULL_CHARACTER, ierr))
    end if

    ! ------------------------------------------------------------------
    ! Materials and output directories
    ! ------------------------------------------------------------------
    call InitialiseMaterials(mats, mesh_ptr%material_ids, mesh_ptr%n_groups, trim(mat_file), .true.)

    nametag = derive_case_nametag(mesh_file)
    nametag = nametag(1:index(nametag,'.vtk')-1)
    run_dir = trim(output_dir) // "/" // trim(solver) // "/" // trim(method_upper) // "/" // trim(nametag)
    call execute_command_line("mkdir -p " // trim(run_dir))

    call execute_command_line("mkdir -p " // trim(run_dir) // "/mesh_cache")
    call write_mesh_to_files(mesh_ptr, trim(run_dir) // "/mesh_cache")

    ! ------------------------------------------------------------------
    ! Solve + output  (solver x method)
    ! ------------------------------------------------------------------

    select case (trim(solver))
    case ("transport")
        select case (trim(method))

        case ("iga")
            call InitialiseGeometry (mesh_iga, FE, QuadSn, QuadFace, TD, sweep_order)
            call InitialiseTransport(mesh_iga, FE, QuadSn, QuadVol, QuadFace, mats, TD)
            call SolveTransport(mesh_iga, mats, FE, QuadSn, TD, &
                                 scalar_flux, ang_flux, k_eff, &
                                 sweep_order, ref_id_list, max_outer, tol, &
                                 is_adjoint, is_eigenvalue)
            call export_transport_vtk_iga(trim(run_dir), trim(nametag), &
                                          mesh_iga, FE, QuadSn, scalar_flux, &
                                          mesh_iga%n_groups, vtk_refine, ang_flux, ang_out)

        case ("fem")
            call InitialiseGeometry_FEM(mesh_fem, FE_fem, QuadSn, QuadFace, TD, sweep_order)
            call InitialiseTransport_FEM(mesh_fem, FE_fem, QuadSn, QuadVol, QuadFace, mats, TD)
            call SolveTransport_FEM(mesh_fem, mats, FE_fem, QuadSn, TD, &
                                     scalar_flux, ang_flux, k_eff, &
                                     sweep_order, ref_id_list, max_outer, tol, &
                                     is_adjoint, is_eigenvalue)
            call export_transport_vtk_fem(trim(run_dir), trim(nametag), &
                                           mesh_fem, FE_fem, QuadSn, scalar_flux, &
                                           mesh_fem%n_groups, vtk_refine, ang_flux, ang_out)

        end select

    case ("diffusion")
        select case (trim(method))

        case ("iga")
            call SolveDiffusion(mesh_iga, FE, QuadVol, QuadFace, mats, &
                                 solver_type, preconditioner, ref_id_list, &
                                 max_outer, tol, is_eigenvalue, is_adjoint, &
                                 phi, k_eff)
            call export_diffusion_vtk_iga(trim(run_dir), trim(nametag), &
                                          mesh_iga, FE, phi, mesh_iga%n_groups, vtk_refine)

        case ("fem")
            call SolveDiffusion_FEM(mesh_fem, FE_fem, QuadVol, QuadFace, mats, &
                                     solver_type, preconditioner, ref_id_list, &
                                     max_outer, tol, is_eigenvalue, is_adjoint, &
                                     phi, k_eff)
            call export_diffusion_vtk_fem(trim(run_dir), trim(nametag), &
                                           mesh_fem, FE_fem, phi, mesh_fem%n_groups, vtk_refine)

        end select

    case default
        write(*,'(A)') "Unknown solver: " // trim(solver) // " (expected 'transport' or 'diffusion')"
        stop 1

    end select

    if (trim(solver) == "diffusion") then
        PetscCall(PetscFinalize(ierr))
    end if

end program maniac
