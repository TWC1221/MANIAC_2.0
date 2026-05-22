! IGA-specific derived types.
!
! t_patch_iga   -- NURBS volume patch (holds knot vectors + CP connectivity)
! t_surface_iga -- NURBS boundary entity (surface in 3D, edge in 2D)
! t_mesh_iga    -- IGA mesh extending t_mesh with NURBS-specific fields
! t_basis_iga   -- NURBS/B-spline basis descriptor: orders, n_basis, face node map
!
! Transport precomputed data (t_fem_dg, t_patch_dg) lives in m_types so it
! is shared across IGA and FEM without duplication.
module m_types_iga
    use m_constants
    use m_types
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
        integer               :: n_elements = 0
        integer,  allocatable :: elem_span_indices(:,:)  ! (2, n_elements)
        integer,  allocatable :: elems(:,:)              ! (n_elements, n_basis)
    end type t_surface_iga

    ! ------------------------------------------------------------------
    ! IGA mesh, extending the common t_mesh base with NURBS-specific
    ! fields.  Common fields (dim, n_nodes, nodes, elems, etc.) are
    ! inherited from t_mesh.
    ! ------------------------------------------------------------------
    type, extends(t_mesh) :: t_mesh_iga
        real(dp), allocatable :: weights(:)    ! (n_nodes) NURBS weights

        type(t_patch_iga),   allocatable :: patches(:)
        type(t_surface_iga), allocatable :: iga_surfaces(:)

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
    type t_basis_iga
        integer :: dim
        integer :: order
        integer :: p_order, q_order, r_order  ! per-direction polynomial degrees
        integer :: n_basis
        integer :: n_nodes_per_face
        integer, allocatable :: face_node_map(:,:)  ! (n_nodes_per_face, n_faces_per_elem)
    end type t_basis_iga

end module m_types_iga
