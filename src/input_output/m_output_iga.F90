! VTK output for IGA solvers (CG diffusion and DG transport).
! Supports 2D (quad cells, VTK type 9) and 3D (hex cells, VTK type 12).
!
! export_diffusion_vtk_iga      -- nodal (CG) flux on IGA mesh
! export_transport_vtk_iga      -- element-DG flux on IGA mesh
! export_transport_vtk_patchiga -- patchwise-DG flux (remapped to element layout)
module m_output_iga
    use m_constants
    use m_types
    use m_types_iga
    use m_utilities
    use m_basis_iga, only: EvalNURBS2D, EvalNURBS3D
    use m_quadrature
    implicit none
    private
    public :: export_diffusion_vtk_iga, export_transport_vtk_iga, &
              export_transport_vtk_patchiga

contains

    subroutine export_diffusion_vtk_iga(outdir, tag, mesh, FE, X_cg, n_groups, refine_level)
        character(len=*),   intent(in) :: outdir, tag
        type(t_mesh_iga),   intent(in) :: mesh
        type(t_basis_iga), intent(in) :: FE
        real(dp),           intent(in) :: X_cg(:,:)
        integer,            intent(in) :: n_groups, refine_level

        character(len=512) :: fpath
        integer :: unit_v, gid

        fpath = trim(outdir) // "/" // trim(tag) // &
                "_iga_p=" // trim(int_to_str(FE%order)) // ".vtk"
        unit_v = 101
        open(unit_v, file=trim(fpath), status='replace', action='write')

        if (mesh%dim == 2) then
            call write_2d_vtk_body(unit_v, mesh, FE, refine_level, "IGA Diffusion Flux", &
                                    n_groups, gid, &
                                    flux_cg=X_cg, write_mat=.true., write_bc=.true.)
        else
            call write_3d_vtk_body(unit_v, mesh, FE, refine_level, "IGA Diffusion Flux", &
                                    n_groups, gid, &
                                    flux_cg=X_cg, write_mat=.true., write_bc=.true.)
        end if

        close(unit_v)
        write(*,'(A)') "  Written: " // trim(fpath)
    end subroutine export_diffusion_vtk_iga

    subroutine export_transport_vtk_iga(outdir, tag, mesh, FE, QuadSn, scalar_flux, &
                                     n_groups, refine_level, label, ang_flux, ang_out)
        character(len=*),      intent(in)           :: outdir, tag
        type(t_mesh_iga),      intent(in)           :: mesh
        type(t_basis_iga),    intent(in)           :: FE
        type(t_sn_quadrature), intent(in)           :: QuadSn
        real(dp),              intent(in)           :: scalar_flux(:,:)
        integer,               intent(in)           :: n_groups, refine_level
        character(len=*), optional, intent(in)      :: label
        real(dp), optional,    intent(in)           :: ang_flux(:,:,:)
        logical,  optional,    intent(in)           :: ang_out

        character(len=512) :: fpath
        character(len=32)  :: lbl
        integer :: unit_v, g, mm, gid, ee, ii, jj, kk, n_sub_nodes, nbasis
        real(dp) :: xi_val, eta_val, zeta_val
        real(dp), allocatable :: xi_g(:), N_e(:), dRu(:), dRv(:), dRw(:)
        integer :: n_angles_export, basep
        logical :: write_ang

        lbl = "iga_fdg"
        if (present(label)) lbl = trim(label)

        ! Evaluate optional ang_out safely: Fortran .and. has no short-circuit,
        ! so accessing an absent optional directly causes SIGSEGV.
        write_ang = .false.
        if (present(ang_out)) write_ang = ang_out

        fpath = trim(outdir) // "/" // trim(tag) // &
                "_" // trim(lbl) // "_p=" // trim(int_to_str(FE%order)) // &
                "_sn=" // trim(int_to_str(QuadSn%order)) // ".vtk"
        unit_v = 102
        open(unit_v, file=trim(fpath), status='replace', action='write')

        if (mesh%dim == 2) then
            call write_2d_vtk_body(unit_v, mesh, FE, refine_level, "IGA Transport Scalar Flux", &
                                    n_groups, gid, &
                                    flux_dg=scalar_flux, write_mat=.true., write_bc=.true.)
        else
            call write_3d_vtk_body(unit_v, mesh, FE, refine_level, "IGA Transport Scalar Flux", &
                                    n_groups, gid, &
                                    flux_dg=scalar_flux, write_mat=.true., write_bc=.true.)
        end if

        ! Angular flux export (shared 2D/3D path using scalar_flux evaluation pattern)
        if (present(ang_flux) .and. write_ang) then
            nbasis = FE%n_basis
            n_angles_export = min(5, QuadSn%n_angles)
            if (mesh%dim == 2) then
                n_sub_nodes = mesh%n_elems * refine_level**2
            else
                n_sub_nodes = mesh%n_elems * refine_level**3
            end if
            allocate(xi_g(refine_level), N_e(nbasis), dRu(nbasis), dRv(nbasis), dRw(nbasis))
            do ii = 1, refine_level
                xi_g(ii) = -1.0_dp + 2.0_dp*real(ii-1,dp)/real(refine_level-1,dp)
            end do

            do mm = 1, n_angles_export
                do g = 1, n_groups
                    write(unit_v,'(A,I0,A,I0)') "SCALARS Ang_Flux_G",g,"_A",mm
                    write(unit_v,'(A)') "double 1"; write(unit_v,'(A)') "LOOKUP_TABLE default"
                    gid = 0
                    do ee = 1, mesh%n_elems
                        basep = (ee-1)*nbasis
                        if (mesh%dim == 2) then
                            do jj = 1, refine_level; do ii = 1, refine_level
                                gid = gid + 1
                                xi_val  = 0.5_dp*((mesh%elem_u_max(ee)-mesh%elem_u_min(ee))*xi_g(ii) + &
                                                   (mesh%elem_u_max(ee)+mesh%elem_u_min(ee)))
                                eta_val = 0.5_dp*((mesh%elem_v_max(ee)-mesh%elem_v_min(ee))*xi_g(jj) + &
                                                   (mesh%elem_v_max(ee)+mesh%elem_v_min(ee)))
                                call EvalNURBS2D(FE, ee, mesh, xi_val, eta_val, N_e, dRu, dRv)
                                write(unit_v,'(F18.10)') dot_product(N_e, ang_flux(basep+1:basep+nbasis,mm,g))
                            end do; end do
                        else
                            do kk = 1, refine_level; do jj = 1, refine_level; do ii = 1, refine_level
                                gid = gid + 1
                                xi_val   = 0.5_dp*((mesh%elem_u_max(ee)-mesh%elem_u_min(ee))*xi_g(ii) + &
                                                    (mesh%elem_u_max(ee)+mesh%elem_u_min(ee)))
                                eta_val  = 0.5_dp*((mesh%elem_v_max(ee)-mesh%elem_v_min(ee))*xi_g(jj) + &
                                                    (mesh%elem_v_max(ee)+mesh%elem_v_min(ee)))
                                zeta_val = 0.5_dp*((mesh%elem_w_max(ee)-mesh%elem_w_min(ee))*xi_g(kk) + &
                                                    (mesh%elem_w_max(ee)+mesh%elem_w_min(ee)))
                                call EvalNURBS3D(FE, ee, mesh, xi_val, eta_val, zeta_val, N_e, dRu, dRv, dRw)
                                write(unit_v,'(F18.10)') dot_product(N_e, ang_flux(basep+1:basep+nbasis,mm,g))
                            end do; end do; end do
                        end if
                    end do
                end do
            end do
            deallocate(xi_g, N_e, dRu, dRv, dRw)
        end if

        close(unit_v)
        write(*,'(A)') "  Written: " // trim(fpath)
    end subroutine export_transport_vtk_iga

    ! ------------------------------------------------------------------
    ! Patchwise-DG VTK output.  Remaps from patch DOF layout to element
    ! layout then delegates to export_transport_vtk_iga.
    ! ------------------------------------------------------------------
    subroutine export_transport_vtk_patchiga(outdir, tag, mesh, FE, QuadSn, &
                                              scalar_flux_patch, PD, n_groups, &
                                              refine_level, ang_flux, ang_out)
        character(len=*),             intent(in)           :: outdir, tag
        type(t_mesh_iga),             intent(in)           :: mesh
        type(t_basis_iga),           intent(in)           :: FE
        type(t_sn_quadrature),        intent(in)           :: QuadSn
        real(dp),                     intent(in)           :: scalar_flux_patch(:,:)
        type(t_patch_dg), intent(in)           :: PD
        integer,                      intent(in)           :: n_groups, refine_level
        real(dp), optional,           intent(in)           :: ang_flux(:,:,:)
        logical,  optional,           intent(in)           :: ang_out

        real(dp), allocatable :: flux_elem(:,:)
        integer :: ee, pp, a, pdof, nb, g

        nb = PD%n_basis_patch
        allocate(flux_elem(mesh%n_elems * FE%n_basis, n_groups))

        do ee = 1, mesh%n_elems
            pp = mesh%elem_patch_id(ee)
            do a = 1, FE%n_basis
                pdof = PD%elem_to_patch_dof(a, ee)
                do g = 1, n_groups
                    flux_elem((ee-1)*FE%n_basis + a, g) = scalar_flux_patch((pp-1)*nb + pdof, g)
                end do
            end do
        end do

        call export_transport_vtk_iga(outdir, tag, mesh, FE, QuadSn, flux_elem, &
                                       n_groups, refine_level, label="iga_pdg")

        deallocate(flux_elem)
    end subroutine export_transport_vtk_patchiga

    ! ------------------------------------------------------------------
    ! 2D VTK body: quad mesh, refine_level^2 nodes per element.
    ! Writes header, POINTS, CELLS, CELL_TYPES, then flux data.
    ! ------------------------------------------------------------------
    subroutine write_2d_vtk_body(unit_v, mesh, FE, refine_level, title, &
                                   n_groups, gid_out, flux_cg, flux_dg, write_mat, write_bc)
        integer,            intent(in)           :: unit_v, refine_level, n_groups
        type(t_mesh_iga),   intent(in)           :: mesh
        type(t_basis_iga), intent(in)           :: FE
        character(len=*),   intent(in)           :: title
        integer,            intent(out)          :: gid_out
        real(dp), optional, intent(in)           :: flux_cg(:,:), flux_dg(:,:)
        logical,  optional, intent(in)           :: write_mat, write_bc

        integer  :: ee, g, ii, jj, gid, cid, n_sub_nodes, n_sub_cells, basep, s, k
        integer  :: n00, n10, n11, n01
        logical  :: do_mat, do_bc
        real(dp) :: u1, u2, v1, v2, xi_val, eta_val
        real(dp), allocatable :: xi_g(:), N_e(:), dRu(:), dRv(:)
        real(dp), allocatable :: Xp(:,:), Up(:,:)
        integer,  allocatable :: Cells(:,:), node_bc(:), elem_bc(:)

        do_mat = present(write_mat) .and. write_mat
        do_bc  = present(write_bc)  .and. write_bc

        n_sub_nodes = mesh%n_elems * refine_level**2
        n_sub_cells = mesh%n_elems * (refine_level-1)**2

        allocate(xi_g(refine_level), N_e(FE%n_basis), dRu(FE%n_basis), dRv(FE%n_basis))
        allocate(Xp(n_sub_nodes, 3), Up(n_sub_nodes, n_groups), Cells(n_sub_cells, 4))
        Xp(:,3) = 0.0_dp

        if (do_bc) then
            allocate(node_bc(mesh%n_nodes), elem_bc(mesh%n_elems))
            node_bc = 0
            do s = 1, size(mesh%iga_surfaces)
                do k = 1, size(mesh%iga_surfaces(s)%cp_ids)
                    node_bc(mesh%iga_surfaces(s)%cp_ids(k)) = mesh%iga_surfaces(s)%bc_id
                end do
            end do
            do ee = 1, mesh%n_elems
                elem_bc(ee) = maxval(node_bc(mesh%elems(ee, 1:FE%n_basis)))
            end do
            deallocate(node_bc)
        end if

        do ii = 1, refine_level
            xi_g(ii) = -1.0_dp + 2.0_dp*real(ii-1,dp)/real(refine_level-1,dp)
        end do

        gid = 0; cid = 0
        do ee = 1, mesh%n_elems
            u1=mesh%elem_u_min(ee); u2=mesh%elem_u_max(ee)
            v1=mesh%elem_v_min(ee); v2=mesh%elem_v_max(ee)
            basep = (ee-1)*FE%n_basis
            do jj = 1, refine_level
                do ii = 1, refine_level
                    gid = gid + 1
                    xi_val  = 0.5_dp*((u2-u1)*xi_g(ii) + (u2+u1))
                    eta_val = 0.5_dp*((v2-v1)*xi_g(jj) + (v2+v1))
                    call EvalNURBS2D(FE, ee, mesh, xi_val, eta_val, N_e, dRu, dRv)
                    Xp(gid,1) = dot_product(N_e, mesh%nodes(mesh%elems(ee,1:FE%n_basis), 1))
                    Xp(gid,2) = dot_product(N_e, mesh%nodes(mesh%elems(ee,1:FE%n_basis), 2))
                    do g = 1, n_groups
                        if (present(flux_cg)) then
                            Up(gid,g) = dot_product(N_e, flux_cg(mesh%elems(ee,1:FE%n_basis), g))
                        else if (present(flux_dg)) then
                            Up(gid,g) = dot_product(N_e, flux_dg(basep+1:basep+FE%n_basis, g))
                        end if
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

        if (do_mat .or. do_bc) then
            write(unit_v,'(A,I10)') "CELL_DATA ", n_sub_cells
            if (do_mat) then
                write(unit_v,'(A)') "SCALARS Material_ID int 1"
                write(unit_v,'(A)') "LOOKUP_TABLE default"
                do ee = 1, mesh%n_elems
                    do ii = 1, (refine_level-1)**2; write(unit_v,'(I10)') mesh%material_ids(ee); end do
                end do
            end if
            if (do_bc) then
                write(unit_v,'(A)') "SCALARS BC_ID int 1"
                write(unit_v,'(A)') "LOOKUP_TABLE default"
                do ee = 1, mesh%n_elems
                    do ii = 1, (refine_level-1)**2; write(unit_v,'(I10)') elem_bc(ee); end do
                end do
                deallocate(elem_bc)
            end if
        end if

        write(unit_v,'(A)') "SCALARS PATCH_ID int 1"
        write(unit_v,'(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            do ii = 1, (refine_level-1)**2; write(unit_v,'(I10)') mesh%elem_patch_id(ee); end do
        end do

        write(unit_v,'(A,I10)') "POINT_DATA ", n_sub_nodes
        do g = 1, n_groups
            write(unit_v,'(A,I0)') "SCALARS Flux_Group_", g
            write(unit_v,'(A)') "double 1"; write(unit_v,'(A)') "LOOKUP_TABLE default"
            do gid = 1, n_sub_nodes; write(unit_v,'(F18.10)') Up(gid,g); end do
        end do

        deallocate(xi_g, N_e, dRu, dRv, Xp, Up, Cells)
    end subroutine write_2d_vtk_body

    ! ------------------------------------------------------------------
    ! 3D VTK body: hex mesh, refine_level^3 nodes per element.
    ! ------------------------------------------------------------------
    subroutine write_3d_vtk_body(unit_v, mesh, FE, refine_level, title, &
                                   n_groups, gid_out, flux_cg, flux_dg, write_mat, write_bc)
        integer,            intent(in)           :: unit_v, refine_level, n_groups
        type(t_mesh_iga),   intent(in)           :: mesh
        type(t_basis_iga), intent(in)           :: FE
        character(len=*),   intent(in)           :: title
        integer,            intent(out)          :: gid_out
        real(dp), optional, intent(in)           :: flux_cg(:,:), flux_dg(:,:)
        logical,  optional, intent(in)           :: write_mat, write_bc

        integer  :: ee, g, ii, jj, kk, gid, cid, n_sub_nodes, n_sub_cells, basep
        integer  :: n000, n100, n110, n010, n001, n101, n111, n011
        real(dp) :: u1, u2, v1, v2, w1, w2, xi_val, eta_val, zeta_val
        real(dp), allocatable :: xi_g(:), N_e(:), dRu(:), dRv(:), dRw(:)
        real(dp), allocatable :: Xp(:,:), Up(:,:)
        integer,  allocatable :: Cells(:,:), node_bc(:), elem_bc(:)
        real(dp) :: local_cps(FE%n_basis, 3)
        logical  :: do_mat, do_bc

        do_mat = present(write_mat) .and. write_mat
        do_bc  = present(write_bc)  .and. write_bc

        if (do_bc) then
            allocate(node_bc(mesh%n_nodes), elem_bc(mesh%n_elems))
            node_bc = 0
            do ii = 1, size(mesh%iga_surfaces)
                do jj = 1, size(mesh%iga_surfaces(ii)%cp_ids)
                    node_bc(mesh%iga_surfaces(ii)%cp_ids(jj)) = mesh%iga_surfaces(ii)%bc_id
                end do
            end do
            do ee = 1, mesh%n_elems
                elem_bc(ee) = maxval(node_bc(mesh%elems(ee, 1:FE%n_basis)))
            end do
            deallocate(node_bc)
        end if


        n_sub_nodes = mesh%n_elems * refine_level**3
        n_sub_cells = mesh%n_elems * (refine_level-1)**3

        allocate(xi_g(refine_level), N_e(FE%n_basis), dRu(FE%n_basis), dRv(FE%n_basis), dRw(FE%n_basis))
        allocate(Xp(n_sub_nodes, 3), Up(n_sub_nodes, n_groups), Cells(n_sub_cells, 8))

        do ii = 1, refine_level
            xi_g(ii) = -1.0_dp + 2.0_dp*real(ii-1,dp)/real(refine_level-1,dp)
        end do

        gid = 0; cid = 0
        do ee = 1, mesh%n_elems
            u1=mesh%elem_u_min(ee); u2=mesh%elem_u_max(ee)
            v1=mesh%elem_v_min(ee); v2=mesh%elem_v_max(ee)
            w1=mesh%elem_w_min(ee); w2=mesh%elem_w_max(ee)
            basep = (ee-1)*FE%n_basis
            do ii = 1, FE%n_basis; local_cps(ii,:) = mesh%nodes(mesh%elems(ee,ii),:); end do
            do kk = 1, refine_level; do jj = 1, refine_level; do ii = 1, refine_level
                gid = gid + 1
                xi_val   = 0.5_dp*((u2-u1)*xi_g(ii) + (u2+u1))
                eta_val  = 0.5_dp*((v2-v1)*xi_g(jj) + (v2+v1))
                zeta_val = 0.5_dp*((w2-w1)*xi_g(kk) + (w2+w1))
                call EvalNURBS3D(FE, ee, mesh, xi_val, eta_val, zeta_val, N_e, dRu, dRv, dRw)
                Xp(gid,1) = dot_product(N_e, local_cps(:,1))
                Xp(gid,2) = dot_product(N_e, local_cps(:,2))
                Xp(gid,3) = dot_product(N_e, local_cps(:,3))
                do g = 1, n_groups
                    if (present(flux_cg)) then
                        Up(gid,g) = dot_product(N_e, flux_cg(mesh%elems(ee,1:FE%n_basis), g))
                    else if (present(flux_dg)) then
                        Up(gid,g) = dot_product(N_e, flux_dg(basep+1:basep+FE%n_basis, g))
                    end if
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

        if (do_mat .or. do_bc) then
            write(unit_v,'(A,I10)') "CELL_DATA ", n_sub_cells
            if (do_mat) then
                write(unit_v,'(A)') "SCALARS Material_ID int 1"
                write(unit_v,'(A)') "LOOKUP_TABLE default"
                do ee = 1, mesh%n_elems
                    do ii = 1, (refine_level-1)**3; write(unit_v,'(I10)') mesh%material_ids(ee); end do
                end do
            end if
            if (do_bc) then
                write(unit_v,'(A)') "SCALARS BC_ID int 1"
                write(unit_v,'(A)') "LOOKUP_TABLE default"
                do ee = 1, mesh%n_elems
                    do ii = 1, (refine_level-1)**3; write(unit_v,'(I10)') elem_bc(ee); end do
                end do
                deallocate(elem_bc)
            end if
        end if

        write(unit_v,'(A)') "SCALARS PATCH_ID int 1"
        write(unit_v,'(A)') "LOOKUP_TABLE default"
        do ee = 1, mesh%n_elems
            do ii = 1, (refine_level-1)**3; write(unit_v,'(I10)') mesh%elem_patch_id(ee); end do
        end do

        write(unit_v,'(A,I10)') "POINT_DATA ", n_sub_nodes
        do g = 1, n_groups
            write(unit_v,'(A,I0)') "SCALARS Flux_Group_", g
            write(unit_v,'(A)') "double 1"; write(unit_v,'(A)') "LOOKUP_TABLE default"
            do gid = 1, n_sub_nodes; write(unit_v,'(F18.10)') Up(gid,g); end do
        end do

        deallocate(xi_g, N_e, dRu, dRv, dRw, Xp, Up, Cells)
    end subroutine write_3d_vtk_body

end module m_output_iga
