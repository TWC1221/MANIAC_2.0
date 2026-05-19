! IGA-specific derived types.
!
! t_patch_iga    -- NURBS volume patch (holds knot vectors + CP connectivity)
! t_surface_iga  -- NURBS boundary entity (surface in 3D, edge in 2D)
! t_mesh_iga     -- full IGA/FEM mesh: geometry and element topology only
! t_finite_iga   -- basis descriptor: orders, n_basis, face node map
! t_transport_iga-- transport-specific precomputed data (separate from mesh)
!
! Keeping mesh and transport data separate means the same mesh object
! works for both diffusion (no transport data needed) and transport.
module m_types_iga
    use m_constants
    implicit none
    public

    ! ------------------------------------------------------------------
    ! NURBS volume patch.  knots_zeta allocated only for 3D meshes.
    ! ------------------------------------------------------------------
    type t_patch_iga
        integer,  allocatable :: cp_ids(:)
        integer               :: material_id = -1
        real(dp), allocatable :: knots_xi(:)
        real(dp), allocatable :: knots_eta(:)
        real(dp), allocatable :: knots_zeta(:)
        integer               :: face_to_surface(6) = 0
    end type t_patch_iga

    ! ------------------------------------------------------------------
    ! NURBS boundary entity.
    ! 3D: 2D face patch (knots_xi and knots_eta allocated).
    ! 2D: 1D edge patch (only knots_xi allocated).
    ! ------------------------------------------------------------------
    type t_surface_iga
        integer,  allocatable :: cp_ids(:)
        integer               :: bc_id = -1
        real(dp), allocatable :: knots_xi(:)
        real(dp), allocatable :: knots_eta(:)
        ! Surface elements — precomputed by m_asmg for 3D surfaces.
        ! Each non-zero-measure span pair becomes one surface element.
        integer               :: n_elements = 0
        integer,  allocatable :: elem_span_indices(:,:)  ! (2, n_elements): xi/eta span indices
        integer,  allocatable :: elems(:,:)              ! (n_elements, n_basis): global CP IDs
    end type t_surface_iga

    ! ------------------------------------------------------------------
    ! Full IGA (or FEM-via-ASMG) mesh.
    ! Contains geometry and element topology only — no solver-specific
    ! precomputed arrays.  Transport precomputed data lives in t_transport_iga.
    ! ------------------------------------------------------------------
    type t_mesh_iga
        integer :: dim               ! 2 or 3
        integer :: n_groups          ! energy groups (from ASMG header)
        integer :: n_nodes           ! total control points / nodes
        integer :: order             ! global polynomial order
        integer :: n_faces_per_elem  ! 4 (2D) or 6 (3D)
        integer :: n_elems

        real(dp), allocatable :: nodes(:,:)    ! (n_nodes, 3)
        real(dp), allocatable :: weights(:)    ! (n_nodes) NURBS weights

        type(t_patch_iga),   allocatable :: patches(:)
        type(t_surface_iga), allocatable :: surfaces(:)

        ! Element topology — one entry per non-zero-measure knot span
        integer, allocatable :: elems(:,:)             ! (n_elems, n_basis)
        integer, allocatable :: material_ids(:)
        integer, allocatable :: elem_patch_id(:)
        integer, allocatable :: elem_span_indices(:,:) ! (dim, n_elems)

        ! Patch-span → element lookup (2D or 3D, one allocated)
        integer, allocatable :: elem_map_2d(:,:,:)     ! (n_patches, max_k, max_k)
        integer, allocatable :: elem_map_3d(:,:,:,:)   ! (n_patches, max_k, max_k, max_k)

        ! Parametric knot-span ranges per element
        real(dp), allocatable :: elem_u_min(:), elem_u_max(:)
        real(dp), allocatable :: elem_v_min(:), elem_v_max(:)
        real(dp), allocatable :: elem_w_min(:), elem_w_max(:)
    end type t_mesh_iga

    ! ------------------------------------------------------------------
    ! NURBS/B-spline basis descriptor.
    ! Basis is evaluated on-the-fly per element (not precomputed globally).
    ! ------------------------------------------------------------------
    type t_finite_iga
        integer :: dim
        integer :: order
        integer :: p_order, q_order, r_order  ! per-direction orders
        integer :: n_basis
        integer :: n_nodes_per_face
        integer, allocatable :: face_node_map(:,:)  ! (n_nodes_per_face, n_faces_per_elem)
    end type t_finite_iga

    ! ------------------------------------------------------------------
    ! Transport-specific precomputed data for an IGA mesh.
    ! Populated by InitialiseGeometry (connectivity/normals) and
    ! InitialiseTransport (integrals, reflective map, LU factors).
    ! ------------------------------------------------------------------
    type t_transport_iga
        ! Face topology and outward normals (from InitialiseGeometry)
        integer,  allocatable :: face_connectivity(:,:,:) ! (4, n_faces, n_elems)
        real(dp), allocatable :: face_normals(:,:,:)      ! (3, n_faces, n_elems)
        integer,  allocatable :: upwind_idx(:,:,:)        ! (n_nodes_per_face, n_faces, n_elems)
        integer,  allocatable :: reflect_map(:,:,:)       ! (n_angles, n_faces, n_elems)

        ! Volume and face integrals (from InitialiseTransport)
        real(dp), allocatable :: elem_mass_matrix(:,:,:)  ! (n_basis, n_basis, n_elems)
        real(dp), allocatable :: elem_stiffness_x(:,:,:)
        real(dp), allocatable :: elem_stiffness_y(:,:,:)
        real(dp), allocatable :: elem_stiffness_z(:,:,:)
        real(dp), allocatable :: face_mass_x(:,:,:,:)     ! (n_basis, n_basis, n_faces, n_elems)
        real(dp), allocatable :: face_mass_y(:,:,:,:)
        real(dp), allocatable :: face_mass_z(:,:,:,:)
        real(dp), allocatable :: basis_integrals_vol(:,:) ! (n_basis, n_elems)

        ! Per-element-per-angle LU factors (from InitialiseTransport)
        real(dp), allocatable :: local_lu(:,:,:,:,:)      ! (n_basis, n_basis, n_elems, n_angles, n_groups)
        integer,  allocatable :: local_pivots(:,:,:,:)    ! (n_basis, n_elems, n_angles, n_groups)
    end type t_transport_iga

end module m_types_iga
