! FEM-specific derived types.
!
! t_mesh_fem   -- Lagrange/FEM mesh, extends t_mesh (no extra fields)
! t_finite_fem -- Lagrange basis descriptor with precomputed reference-element arrays
!
! Transport precomputed data lives in t_transport_data (m_types) so it is
! shared across IGA and FEM without duplication.
module m_types_fem
    use m_constants
    use m_types
    implicit none
    public

    ! ------------------------------------------------------------------
    ! FEM mesh: inherits all common fields from t_mesh.
    ! No IGA-specific data (knot vectors, NURBS weights, span maps).
    ! mesh%surfaces holds boundary entities with cp_ids and bc_id.
    ! ------------------------------------------------------------------
    type, extends(t_mesh) :: t_mesh_fem
    end type t_mesh_fem

    ! ------------------------------------------------------------------
    ! Lagrange basis descriptor.
    ! Unlike IGA, basis functions are uniform across all elements and are
    ! precomputed once on the reference element [-1,1]^d, stored here.
    ! GetMapping2D/3D_FEM index directly into these arrays by quad point,
    ! avoiding per-element knot-span evaluation.
    ! ------------------------------------------------------------------
    type t_finite_fem
        integer :: dim
        integer :: order
        integer :: p_order, q_order, r_order   ! per-direction (all equal for FEM)
        integer :: n_basis                      ! (p+1)^dim
        integer :: n_nodes_per_face             ! (p+1)^(dim-1)

        real(dp), allocatable :: node_roots(:)       ! equispaced Lagrange nodes in [-1,1]
        integer,  allocatable :: face_node_map(:,:)  ! (n_nodes_per_face, n_faces_per_elem)

        ! Precomputed at volume quadrature points -- (n_quad_vol, n_basis)
        real(dp), allocatable :: basis_at_quad(:,:)
        real(dp), allocatable :: dbasis_dxi(:,:)
        real(dp), allocatable :: dbasis_deta(:,:)
        real(dp), allocatable :: dbasis_dzeta(:,:)  ! 3D only

        ! Precomputed at face quadrature points -- (n_quad_face, n_nodes_per_face)
        real(dp), allocatable :: basis_at_face_quad(:,:)
        real(dp), allocatable :: dbasis_face_dxi(:,:)
        real(dp), allocatable :: dbasis_face_deta(:,:)  ! 3D only
    end type t_finite_fem

end module m_types_fem
