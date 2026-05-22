! Shared derived types used across all solvers (diffusion/transport, FEM/IGA).
!
! Mesh hierarchy:    t_mesh  ←  t_mesh_iga  /  t_mesh_fem
! Transport data:    t_fem_dg  (span-level DG, IGA and FEM)
!                    t_patch_dg (patch-level DG, IGA only)
module m_types
    use m_constants
    implicit none

    ! ------------------------------------------------------------------
    ! Boundary surface: node IDs and BC tag only.
    ! IGA surfaces carry additional NURBS data in t_surface_iga.
    ! ------------------------------------------------------------------
    type :: t_surface_mesh
        integer, allocatable :: cp_ids(:)
        integer              :: bc_id = -1
    end type t_surface_mesh

    ! ------------------------------------------------------------------
    ! Common geometry and element topology, shared by all physics modes.
    ! IGA extends this with NURBS fields; FEM uses this directly.
    ! ------------------------------------------------------------------
    type :: t_mesh
        integer :: dim               = 0
        integer :: n_groups          = 1
        integer :: n_nodes           = 0
        integer :: order             = 1
        integer :: n_faces_per_elem  = 4
        integer :: n_elems           = 0
        logical :: DG                = .false.
        real(dp), allocatable :: nodes(:,:)         ! (n_nodes, 3)
        integer,  allocatable :: elems(:,:)         ! (n_elems, n_basis)
        integer,  allocatable :: material_ids(:)
        type(t_surface_mesh), allocatable :: surfaces(:)
    end type t_mesh

    ! ------------------------------------------------------------------
    ! Nuclear cross-section data for one material, all energy groups.
    ! ------------------------------------------------------------------
    type t_material
        character(len=32)     :: name
        real(dp), allocatable :: D(:)        ! diffusion coefficient [cm]
        real(dp), allocatable :: SigmaT(:)   ! total XS [cm^-1]  (transport)
        real(dp), allocatable :: SigmaR(:)   ! removal XS [cm^-1] (diffusion)
        real(dp), allocatable :: SigA(:)     ! absorption XS [cm^-1]
        real(dp), allocatable :: SigF(:)     ! fission XS [cm^-1]
        real(dp), allocatable :: Nu(:)       ! average fission production [-]
        real(dp), allocatable :: NuSigF(:)   ! nu*fission XS [cm^-1]
        real(dp), allocatable :: Chi(:)      ! fission spectrum [-]
        real(dp), allocatable :: Src(:)      ! fixed source [cm^-3 s^-1]
        real(dp), allocatable :: SigmaS(:,:) ! scatter matrix SigmaS(from_g, to_g)
    end type t_material

    ! ------------------------------------------------------------------
    ! Spatial quadrature rule.
    ! ------------------------------------------------------------------
    type t_quadrature
        integer :: n_points
        real(dp), allocatable :: xi(:), eta(:), zeta(:), weights(:)
    end type t_quadrature

    ! ------------------------------------------------------------------
    ! Discrete-ordinates (Sn) angular quadrature set.
    ! ------------------------------------------------------------------
    type t_sn_quadrature
        integer :: order
        integer :: n_angles
        real(dp), allocatable :: dirs(:,:)    ! (n_angles, 3)
        real(dp), allocatable :: weights(:)   ! (n_angles)
    end type t_sn_quadrature

    ! ------------------------------------------------------------------
    ! Boundary condition for one physical surface ID.
    ! ------------------------------------------------------------------
    type t_bc_config
        integer  :: mat_id
        integer  :: bc_type   ! BC_VACUUM, BC_REFLECTIVE, BC_DIRICHLET, BC_ALBEDO
        real(dp) :: value
    end type t_bc_config

    ! ------------------------------------------------------------------
    ! Span-level DG transport precomputed data.
    ! Populated by InitialiseGeometry (connectivity/normals) and
    ! InitialiseTransport (integrals, reflective map, LU factors).
    ! Shared by IGA (span = knot-span element) and FEM.
    ! ------------------------------------------------------------------
    type t_fem_dg
        ! Face topology and outward normals
        integer,  allocatable :: face_connectivity(:,:,:) ! (4, n_faces, n_elems)
        real(dp), allocatable :: face_normals(:,:,:)      ! (3, n_faces, n_elems)
        integer,  allocatable :: upwind_idx(:,:,:)        ! (n_nodes_per_face, n_faces, n_elems)
        integer,  allocatable :: reflect_map(:,:,:)       ! (n_angles, n_faces, n_elems)

        ! Volume and face integrals
        real(dp), allocatable :: elem_mass_matrix(:,:,:)  ! (n_basis, n_basis, n_elems)
        real(dp), allocatable :: elem_stiffness_x(:,:,:)
        real(dp), allocatable :: elem_stiffness_y(:,:,:)
        real(dp), allocatable :: elem_stiffness_z(:,:,:)
        real(dp), allocatable :: face_mass_x(:,:,:,:)     ! (n_basis, n_basis, n_faces, n_elems)
        real(dp), allocatable :: face_mass_y(:,:,:,:)
        real(dp), allocatable :: face_mass_z(:,:,:,:)
        real(dp), allocatable :: basis_integrals_vol(:,:) ! (n_basis, n_elems)

        ! Per-element-per-angle LU factors
        real(dp), allocatable :: local_lu(:,:,:,:,:)      ! (n_basis, n_basis, n_elems, n_angles, n_groups)
        integer,  allocatable :: local_pivots(:,:,:,:)    ! (n_basis, n_elems, n_angles, n_groups)
    end type t_fem_dg

    ! ------------------------------------------------------------------
    ! Patch-level DG-IGA transport precomputed data.
    ! Each NURBS patch is one DG "element"; C^(p-1) smoothness is
    ! retained within the patch, DG upwind coupling at patch interfaces.
    !
    ! DOF layout: patch pp → global DOFs (pp-1)*n_basis_patch+1 .. pp*n_basis_patch
    ! ------------------------------------------------------------------
    type t_patch_dg
        integer :: n_basis_patch     ! DOFs per patch (n_cp_xi * n_cp_eta [* n_cp_zeta])
        integer :: n_face_basis_max  ! max DOFs per patch face across all face directions

        integer, allocatable :: n_face_basis_f(:)          ! (n_faces_per_patch)
        integer, allocatable :: face_node_map_patch(:,:)   ! (n_face_basis_max, n_faces)
        integer, allocatable :: elem_to_patch_dof(:,:)     ! (n_basis, n_elems)

        ! Patch element lists (COO-style, sorted by patch)
        integer, allocatable :: patch_elem_start(:)        ! (n_patches+1)
        integer, allocatable :: patch_elem_list(:)         ! (n_elems)

        ! Volume matrices (n_basis_patch, n_basis_patch, n_patches)
        real(dp), allocatable :: patch_mass(:,:,:)
        real(dp), allocatable :: patch_stiff_x(:,:,:)
        real(dp), allocatable :: patch_stiff_y(:,:,:)
        real(dp), allocatable :: patch_stiff_z(:,:,:)
        real(dp), allocatable :: basis_integrals_vol(:,:)  ! (n_basis_patch, n_patches)

        ! Patch-boundary face matrices (n_basis_patch, n_basis_patch, n_faces, n_patches)
        real(dp), allocatable :: face_mass_x(:,:,:,:)
        real(dp), allocatable :: face_mass_y(:,:,:,:)
        real(dp), allocatable :: face_mass_z(:,:,:,:)

        ! Patch-level topology
        integer,  allocatable :: face_connectivity(:,:,:)  ! (4, n_faces, n_patches)
        real(dp), allocatable :: face_normals(:,:,:)       ! (3, n_faces, n_patches)
        integer,  allocatable :: upwind_idx(:,:,:)         ! (n_face_basis_max, n_faces, n_patches)
        integer,  allocatable :: reflect_map(:,:,:)        ! (n_angles, n_faces, n_patches)

        ! Per-patch-per-angle LU factors
        real(dp), allocatable :: local_lu(:,:,:,:,:)       ! (nb, nb, n_patches, n_angles, n_groups)
        integer,  allocatable :: local_pivots(:,:,:,:)     ! (nb, n_patches, n_angles, n_groups)
    end type t_patch_dg

end module m_types
