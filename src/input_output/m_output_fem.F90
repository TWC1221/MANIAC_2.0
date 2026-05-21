! VTK output for FEM transport and diffusion.
! Supports 2D (quad cells, VTK type 9) and 3D (hex cells, VTK type 12).
!
! export_transport_vtk_fem  -- element-DG scalar flux on a FEM mesh
! export_diffusion_vtk_fem  -- CG nodal flux on a FEM mesh
module m_output_fem
    use m_constants
    use m_types_fem
    use m_utilities
    use m_basis_fem, only: EvalLagrange2D, EvalLagrange3D
    use m_types,      only: t_sn_quadrature
    implicit none
    private
    public :: export_transport_vtk_fem
    public :: export_diffusion_vtk_fem

contains

    subroutine export_transport_vtk_fem(outdir, tag, mesh, FE, QuadSn, scalar_flux, &
                                         n_groups, refine_level, ang_flux, ang_out)
        character(len=*),      intent(in)        :: outdir, tag
        type(t_mesh_fem),      intent(in)        :: mesh
        type(t_finite_fem),    intent(in)        :: FE
        type(t_sn_quadrature), intent(in)        :: QuadSn
        real(dp),              intent(in)        :: scalar_flux(:,:)
        integer,               intent(in)        :: n_groups, refine_level
        real(dp), optional,    intent(in)        :: ang_flux(:,:,:)
        logical,  optional,    intent(in)        :: ang_out

        character(len=512) :: fpath
        integer :: unit_v, gid

        fpath = trim(outdir) // "/" // trim(tag) // &
                "_fem_n=" // trim(int_to_str(FE%order)) // &
                "_sn=" // trim(int_to_str(QuadSn%order)) // ".vtk"
        unit_v = 103
        open(unit_v, file=trim(fpath), status='replace', action='write')

        if (mesh%dim == 2) then
            call write_2d_vtk_body_fem(unit_v, mesh, FE, refine_level, &
                                        "FEM Transport Scalar Flux", n_groups, gid, &
                                        flux_dg=scalar_flux)
        else
            call write_3d_vtk_body_fem(unit_v, mesh, FE, refine_level, &
                                        "FEM Transport Scalar Flux", n_groups, gid, &
                                        flux_dg=scalar_flux)
        end if

        ! Angular flux export (first few angles)
        if (present(ang_flux) .and. ang_out) then
            call write_ang_flux_fem(unit_v, mesh, FE, QuadSn, ang_flux, n_groups, refine_level)
        end if

        close(unit_v)
        write(*,'(A)') "  Written: " // trim(fpath)
    end subroutine export_transport_vtk_fem

    ! ------------------------------------------------------------------
    ! 2D VTK body: Lagrange quad mesh, refine_level^2 sub-nodes per element.
    ! ------------------------------------------------------------------
    subroutine write_2d_vtk_body_fem(unit_v, mesh, FE, refine_level, title, &
                                      n_groups, gid_out, flux_dg)
        integer,            intent(in)           :: unit_v, refine_level, n_groups
        type(t_mesh_fem),   intent(in)           :: mesh
        type(t_finite_fem), intent(in)           :: FE
        character(len=*),   intent(in)           :: title
        integer,            intent(out)          :: gid_out
        real(dp), optional, intent(in)           :: flux_dg(:,:)

        integer  :: ee, g, ii, jj, gid, cid, n_sub_nodes, n_sub_cells, basep
        integer  :: n00, n10, n11, n01, s, k
        real(dp) :: xi_val, eta_val
        real(dp), allocatable :: xi_g(:), N_e(:), dNdxi(:), dNdeta(:)
        real(dp), allocatable :: Xp(:,:), Up(:,:)
        integer,  allocatable :: Cells(:,:), node_bc(:), elem_bc(:)
        real(dp) :: local_nodes(FE%n_basis, 3)

        n_sub_nodes = mesh%n_elems * refine_level**2
        n_sub_cells = mesh%n_elems * (refine_level-1)**2

        allocate(xi_g(refine_level))
        allocate(N_e(FE%n_basis), dNdxi(FE%n_basis), dNdeta(FE%n_basis))
        allocate(Xp(n_sub_nodes, 3), Up(n_sub_nodes, n_groups), Cells(n_sub_cells, 4))
        Xp(:,3) = 0.0_dp

        ! BC data from base mesh%surfaces
        allocate(node_bc(mesh%n_nodes), elem_bc(mesh%n_elems))
        node_bc = 0
        do s = 1, size(mesh%surfaces)
            do k = 1, size(mesh%surfaces(s)%cp_ids)
                node_bc(mesh%surfaces(s)%cp_ids(k)) = mesh%surfaces(s)%bc_id
            end do
        end do
        do ee = 1, mesh%n_elems
            elem_bc(ee) = maxval(node_bc(mesh%elems(ee, 1:FE%n_basis)))
        end do
        deallocate(node_bc)

        do ii = 1, refine_level
            xi_g(ii) = -1.0_dp + 2.0_dp*real(ii-1,dp)/real(refine_level-1,dp)
        end do

        gid = 0; cid = 0
        do ee = 1, mesh%n_elems
            basep = (ee-1)*FE%n_basis
            local_nodes = mesh%nodes(mesh%elems(ee, 1:FE%n_basis), :)
            do jj = 1, refine_level
                do ii = 1, refine_level
                    gid    = gid + 1
                    xi_val  = xi_g(ii)
                    eta_val = xi_g(jj)
                    call EvalLagrange2D(FE, xi_val, eta_val, N_e, dNdxi, dNdeta)
                    Xp(gid,1) = dot_product(N_e, local_nodes(:,1))
                    Xp(gid,2) = dot_product(N_e, local_nodes(:,2))
                    do g = 1, n_groups
                        if (present(flux_dg)) &
                            Up(gid,g) = dot_product(N_e, flux_dg(basep+1:basep+FE%n_basis, g))
                    end do
                end do
            end do

            basep = (ee-1)*refine_level**2
            do jj = 1, refine_level-1
                do ii = 1, refine_level-1
                    n00 = basep + (jj-1)*refine_level + ii
                    n10 = n00 + 1
                    n11 = n00 + refine_level + 1
                    n01 = n00 + refine_level
                    cid = cid + 1
                    Cells(cid,:) = [n00, n10, n11, n01]
                end do
            end do
        end do
        gid_out = gid

        write(unit_v,'(A)') "# vtk DataFile Version 3.0"
        write(unit_v,'(A)') trim(title)
        write(unit_v,'(A)') "ASCII"
        write(unit_v,'(A)') "DATASET UNSTRUCTURED_GRID"
        write(unit_v,'(A,I0,A)') "POINTS ", n_sub_nodes, " double"
        do gid = 1, n_sub_nodes; write(unit_v,'(3F18.10)') Xp(gid,:); end do
        write(unit_v,'(A,2I10)') "CELLS ", n_sub_cells, n_sub_cells*5
        do cid = 1, n_sub_cells
            write(unit_v,'(5I10)') 4, Cells(cid,1)-1, Cells(cid,2)-1, Cells(cid,3)-1, Cells(cid,4)-1
        end do
        write(unit_v,'(A,I10)') "CELL_TYPES ", n_sub_cells
        do cid = 1, n_sub_cells; write(unit_v,'(I2)') 9; end do  ! VTK_QUAD

        write(unit_v,'(A,I10)') "CELL_DATA ", n_sub_cells
        write(unit_v,'(A)') "SCALARS Material_ID int 1"
        write(unit_v,'(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            do ii = 1, (refine_level-1)**2; write(unit_v,'(I10)') mesh%material_ids(ee); end do
        end do
        write(unit_v,'(A)') "SCALARS BC_ID int 1"
        write(unit_v,'(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            do ii = 1, (refine_level-1)**2; write(unit_v,'(I10)') elem_bc(ee); end do
        end do
        deallocate(elem_bc)

        write(unit_v,'(A,I10)') "POINT_DATA ", n_sub_nodes
        do g = 1, n_groups
            write(unit_v,'(A,I0)') "SCALARS Flux_Group_", g
            write(unit_v,'(A)') "double 1"; write(unit_v,'(A)') "LOOKUP_TABLE default"
            do gid = 1, n_sub_nodes; write(unit_v,'(F18.10)') Up(gid,g); end do
        end do

        deallocate(xi_g, N_e, dNdxi, dNdeta, Xp, Up, Cells)
    end subroutine write_2d_vtk_body_fem

    ! ------------------------------------------------------------------
    ! 3D VTK body: Lagrange hex mesh, refine_level^3 sub-nodes per element.
    ! ------------------------------------------------------------------
    subroutine write_3d_vtk_body_fem(unit_v, mesh, FE, refine_level, title, &
                                      n_groups, gid_out, flux_dg)
        integer,            intent(in)           :: unit_v, refine_level, n_groups
        type(t_mesh_fem),   intent(in)           :: mesh
        type(t_finite_fem), intent(in)           :: FE
        character(len=*),   intent(in)           :: title
        integer,            intent(out)          :: gid_out
        real(dp), optional, intent(in)           :: flux_dg(:,:)

        integer  :: ee, g, ii, jj, kk, gid, cid, n_sub_nodes, n_sub_cells, basep, s, k
        integer  :: n000, n100, n110, n010, n001, n101, n111, n011
        real(dp) :: xi_val, eta_val, zeta_val
        real(dp), allocatable :: xi_g(:), N_e(:), dNdxi(:), dNdeta(:), dNdzeta(:)
        real(dp), allocatable :: Xp(:,:), Up(:,:)
        integer,  allocatable :: Cells(:,:), node_bc(:), elem_bc(:)
        real(dp) :: local_nodes(FE%n_basis, 3)

        n_sub_nodes = mesh%n_elems * refine_level**3
        n_sub_cells = mesh%n_elems * (refine_level-1)**3

        allocate(xi_g(refine_level))
        allocate(N_e(FE%n_basis), dNdxi(FE%n_basis), dNdeta(FE%n_basis), dNdzeta(FE%n_basis))
        allocate(Xp(n_sub_nodes, 3), Up(n_sub_nodes, n_groups), Cells(n_sub_cells, 8))

        allocate(node_bc(mesh%n_nodes), elem_bc(mesh%n_elems))
        node_bc = 0
        do s = 1, size(mesh%surfaces)
            do k = 1, size(mesh%surfaces(s)%cp_ids)
                node_bc(mesh%surfaces(s)%cp_ids(k)) = mesh%surfaces(s)%bc_id
            end do
        end do
        do ee = 1, mesh%n_elems
            elem_bc(ee) = maxval(node_bc(mesh%elems(ee, 1:FE%n_basis)))
        end do
        deallocate(node_bc)

        do ii = 1, refine_level
            xi_g(ii) = -1.0_dp + 2.0_dp*real(ii-1,dp)/real(refine_level-1,dp)
        end do

        gid = 0; cid = 0
        do ee = 1, mesh%n_elems
            basep = (ee-1)*FE%n_basis
            local_nodes = mesh%nodes(mesh%elems(ee, 1:FE%n_basis), :)
            do kk = 1, refine_level; do jj = 1, refine_level; do ii = 1, refine_level
                gid      = gid + 1
                xi_val   = xi_g(ii)
                eta_val  = xi_g(jj)
                zeta_val = xi_g(kk)
                call EvalLagrange3D(FE, xi_val, eta_val, zeta_val, N_e, dNdxi, dNdeta, dNdzeta)
                Xp(gid,1) = dot_product(N_e, local_nodes(:,1))
                Xp(gid,2) = dot_product(N_e, local_nodes(:,2))
                Xp(gid,3) = dot_product(N_e, local_nodes(:,3))
                do g = 1, n_groups
                    if (present(flux_dg)) &
                        Up(gid,g) = dot_product(N_e, flux_dg(basep+1:basep+FE%n_basis, g))
                end do
            end do; end do; end do

            basep = (ee-1)*refine_level**3
            do kk = 1, refine_level-1; do jj = 1, refine_level-1; do ii = 1, refine_level-1
                n000 = basep + (kk-1)*refine_level**2 + (jj-1)*refine_level + ii
                n100=n000+1; n110=n000+refine_level+1; n010=n000+refine_level
                n001=n000+refine_level**2; n101=n001+1; n111=n001+refine_level+1; n011=n001+refine_level
                cid = cid + 1
                Cells(cid,:) = [n000,n100,n110,n010,n001,n101,n111,n011]
            end do; end do; end do
        end do
        gid_out = gid

        write(unit_v,'(A)') "# vtk DataFile Version 3.0"
        write(unit_v,'(A)') trim(title)
        write(unit_v,'(A)') "ASCII"
        write(unit_v,'(A)') "DATASET UNSTRUCTURED_GRID"
        write(unit_v,'(A,I0,A)') "POINTS ", n_sub_nodes, " double"
        do gid = 1, n_sub_nodes; write(unit_v,'(3F18.10)') Xp(gid,:); end do
        write(unit_v,'(A,2I10)') "CELLS ", n_sub_cells, n_sub_cells*9
        do cid = 1, n_sub_cells
            write(unit_v,'(9I10)') 8, Cells(cid,1)-1, Cells(cid,2)-1, Cells(cid,3)-1, Cells(cid,4)-1, &
                                      Cells(cid,5)-1, Cells(cid,6)-1, Cells(cid,7)-1, Cells(cid,8)-1
        end do
        write(unit_v,'(A,I10)') "CELL_TYPES ", n_sub_cells
        do cid = 1, n_sub_cells; write(unit_v,'(I2)') 12; end do  ! VTK_HEXAHEDRON

        write(unit_v,'(A,I10)') "CELL_DATA ", n_sub_cells
        write(unit_v,'(A)') "SCALARS Material_ID int 1"
        write(unit_v,'(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            do ii = 1, (refine_level-1)**3; write(unit_v,'(I10)') mesh%material_ids(ee); end do
        end do
        write(unit_v,'(A)') "SCALARS BC_ID int 1"
        write(unit_v,'(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            do ii = 1, (refine_level-1)**3; write(unit_v,'(I10)') elem_bc(ee); end do
        end do
        deallocate(elem_bc)

        write(unit_v,'(A,I10)') "POINT_DATA ", n_sub_nodes
        do g = 1, n_groups
            write(unit_v,'(A,I0)') "SCALARS Flux_Group_", g
            write(unit_v,'(A)') "double 1"; write(unit_v,'(A)') "LOOKUP_TABLE default"
            do gid = 1, n_sub_nodes; write(unit_v,'(F18.10)') Up(gid,g); end do
        end do

        deallocate(xi_g, N_e, dNdxi, dNdeta, dNdzeta, Xp, Up, Cells)
    end subroutine write_3d_vtk_body_fem

    ! ------------------------------------------------------------------
    ! Append angular flux scalars to an already-open VTK file.
    ! Writes the first min(5, n_angles) directions.
    ! ------------------------------------------------------------------
    subroutine write_ang_flux_fem(unit_v, mesh, FE, QuadSn, ang_flux, n_groups, refine_level)
        integer,               intent(in) :: unit_v, n_groups, refine_level
        type(t_mesh_fem),      intent(in) :: mesh
        type(t_finite_fem),    intent(in) :: FE
        type(t_sn_quadrature), intent(in) :: QuadSn
        real(dp),              intent(in) :: ang_flux(:,:,:)

        integer  :: mm, g, ee, ii, jj, kk, gid, basep, n_angles_export
        real(dp) :: xi_val, eta_val, zeta_val
        real(dp), allocatable :: xi_g(:), N_e(:), dNdxi(:), dNdeta(:), dNdzeta(:)
        real(dp) :: local_nodes(FE%n_basis, 3)

        n_angles_export = min(5, QuadSn%n_angles)
        allocate(xi_g(refine_level), N_e(FE%n_basis))
        allocate(dNdxi(FE%n_basis), dNdeta(FE%n_basis), dNdzeta(FE%n_basis))
        do ii = 1, refine_level
            xi_g(ii) = -1.0_dp + 2.0_dp*real(ii-1,dp)/real(refine_level-1,dp)
        end do

        do mm = 1, n_angles_export
            do g = 1, n_groups
                write(unit_v,'(A,I0,A,I0)') "SCALARS Ang_Flux_G",g,"_A",mm
                write(unit_v,'(A)') "double 1"; write(unit_v,'(A)') "LOOKUP_TABLE default"
                gid = 0
                do ee = 1, mesh%n_elems
                    basep = (ee-1)*FE%n_basis
                    local_nodes = mesh%nodes(mesh%elems(ee, 1:FE%n_basis), :)
                    if (mesh%dim == 2) then
                        do jj = 1, refine_level; do ii = 1, refine_level
                            gid     = gid + 1
                            xi_val  = xi_g(ii); eta_val = xi_g(jj)
                            call EvalLagrange2D(FE, xi_val, eta_val, N_e, dNdxi, dNdeta)
                            write(unit_v,'(F18.10)') dot_product(N_e, ang_flux(basep+1:basep+FE%n_basis,mm,g))
                        end do; end do
                    else
                        do kk = 1, refine_level; do jj = 1, refine_level; do ii = 1, refine_level
                            gid      = gid + 1
                            xi_val   = xi_g(ii); eta_val = xi_g(jj); zeta_val = xi_g(kk)
                            call EvalLagrange3D(FE, xi_val, eta_val, zeta_val, N_e, dNdxi, dNdeta, dNdzeta)
                            write(unit_v,'(F18.10)') dot_product(N_e, ang_flux(basep+1:basep+FE%n_basis,mm,g))
                        end do; end do; end do
                    end if
                end do
            end do
        end do
        deallocate(xi_g, N_e, dNdxi, dNdeta, dNdzeta)
    end subroutine write_ang_flux_fem

    ! ------------------------------------------------------------------
    ! CG diffusion VTK output.
    ! phi(n_nodes, n_groups) is gathered into DG element layout so the
    ! existing sub-grid VTK body writers (which already produce correct
    ! VTK node ordering) can be reused unchanged.
    ! ------------------------------------------------------------------
    subroutine export_diffusion_vtk_fem(outdir, tag, mesh, FE, phi, n_groups, refine_level)
        character(len=*),   intent(in) :: outdir, tag
        type(t_mesh_fem),   intent(in) :: mesh
        type(t_finite_fem), intent(in) :: FE
        real(dp),           intent(in) :: phi(:,:)
        integer,            intent(in) :: n_groups, refine_level

        character(len=512) :: fpath
        integer :: unit_v, ee, k, gid_dummy
        real(dp), allocatable :: phi_dg(:,:)

        ! Gather CG nodal flux into DG element layout:
        !   phi_dg((ee-1)*n_basis + k, g) = phi(global_node_k_of_ee, g)
        allocate(phi_dg(mesh%n_elems * FE%n_basis, n_groups))
        do ee = 1, mesh%n_elems
            do k = 1, FE%n_basis
                phi_dg((ee-1)*FE%n_basis + k, :) = phi(mesh%elems(ee,k), :)
            end do
        end do

        fpath = trim(outdir) // "/" // trim(tag) // &
                "_fem_diff_n=" // trim(int_to_str(FE%order)) // ".vtk"
        unit_v = 104
        open(unit_v, file=trim(fpath), status='replace', action='write')

        if (mesh%dim == 2) then
            call write_2d_vtk_body_fem(unit_v, mesh, FE, refine_level, &
                                        "FEM Diffusion Scalar Flux", n_groups, gid_dummy, &
                                        flux_dg=phi_dg)
        else
            call write_3d_vtk_body_fem(unit_v, mesh, FE, refine_level, &
                                        "FEM Diffusion Scalar Flux", n_groups, gid_dummy, &
                                        flux_dg=phi_dg)
        end if

        close(unit_v)
        deallocate(phi_dg)
        write(*,'(A)') "  Written: " // trim(fpath)
    end subroutine export_diffusion_vtk_fem

end module m_output_fem
