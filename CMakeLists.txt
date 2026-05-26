cmake_minimum_required(VERSION 3.18)
project(MANIAC LANGUAGES Fortran C)

# ---------------------------------------------------------------------------
# Build type
# ---------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# MPI (let CMake find the right libraries automatically)
# ---------------------------------------------------------------------------
find_package(MPI REQUIRED COMPONENTS Fortran C)

# ---------------------------------------------------------------------------
# PETSc
# ---------------------------------------------------------------------------
set(PETSC_DIR "/opt/homebrew/opt/petsc" CACHE PATH "PETSc installation prefix")
set(PETSC_INC "${PETSC_DIR}/include")
set(PETSC_LIB_DIR "${PETSC_DIR}/lib")

# ---------------------------------------------------------------------------
# Compiler flags
# ---------------------------------------------------------------------------
set(CMAKE_Fortran_FLAGS
    "-Wall -ffree-line-length-none -I${PETSC_INC}")
set(CMAKE_Fortran_FLAGS_DEBUG
    "-fcheck=all -Og -fopenmp -g")
set(CMAKE_Fortran_FLAGS_RELEASE
    "-O3 -fopenmp -march=native")
set(CMAKE_EXE_LINKER_FLAGS
    "-L${PETSC_LIB_DIR} -Wl,-rpath,${PETSC_LIB_DIR}")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mods)

# ---------------------------------------------------------------------------
# Core library
# ---------------------------------------------------------------------------
add_library(maniac_core STATIC
    src/core/m_constants.F90
    src/core/m_types.F90
    src/core/m_material.F90
    src/core/m_quadrature.F90
    src/core/m_utilities.F90
)
set_target_properties(maniac_core PROPERTIES Fortran_PREPROCESS ON)
target_include_directories(maniac_core PUBLIC ${PETSC_INC})
target_link_libraries(maniac_core PUBLIC MPI::MPI_Fortran)

# ---------------------------------------------------------------------------
# IGA layer
# ---------------------------------------------------------------------------
add_library(maniac_iga STATIC
    src/iga/m_basis_iga.F90
)
set_target_properties(maniac_iga PROPERTIES Fortran_PREPROCESS ON)
target_link_libraries(maniac_iga PUBLIC maniac_core)

# ---------------------------------------------------------------------------
# FEM layer
# ---------------------------------------------------------------------------
add_library(maniac_fem STATIC
    src/fem/m_basis_fem.F90
)
set_target_properties(maniac_fem PROPERTIES Fortran_PREPROCESS ON)
target_link_libraries(maniac_fem PUBLIC maniac_core)

# ---------------------------------------------------------------------------
# Solver utilities
# ---------------------------------------------------------------------------
add_library(maniac_solver STATIC
    src/solver/m_power_iteration.F90
)
set_target_properties(maniac_solver PROPERTIES Fortran_PREPROCESS ON)
target_include_directories(maniac_solver PUBLIC ${PETSC_INC})
target_link_libraries(maniac_solver PUBLIC maniac_core)

# ---------------------------------------------------------------------------
# Diffusion solver
# ---------------------------------------------------------------------------
add_library(maniac_diffusion STATIC
    src/diffusion/m_petsc.F90
    src/diffusion/m_diffusion_iga.F90
    src/diffusion/m_diffusion_fem.F90
)
set_target_properties(maniac_diffusion PROPERTIES Fortran_PREPROCESS ON)
target_include_directories(maniac_diffusion PUBLIC ${PETSC_INC})
target_link_libraries(maniac_diffusion PUBLIC
    maniac_solver maniac_iga maniac_fem maniac_core)

# ---------------------------------------------------------------------------
# Transport solver
# ---------------------------------------------------------------------------
add_library(maniac_transport STATIC
    src/transport/m_sweep_order.F90
    src/transport/m_precompute_iga.F90
    src/transport/m_transport_iga_pdg.F90
    src/transport/m_precompute_fem.F90
    src/transport/m_transport_fem.F90
)
set_target_properties(maniac_transport PROPERTIES Fortran_PREPROCESS ON)
target_link_libraries(maniac_transport PUBLIC maniac_solver maniac_iga maniac_fem maniac_core)

# ---------------------------------------------------------------------------
# Output utilities
# ---------------------------------------------------------------------------
add_library(maniac_input_output STATIC
    src/input_output/m_asmg.F90
    src/input_output/m_output_iga.F90
    src/input_output/m_output_fem.F90
)
set_target_properties(maniac_input_output PROPERTIES Fortran_PREPROCESS ON)
target_link_libraries(maniac_input_output PUBLIC maniac_iga maniac_fem maniac_core)

# ---------------------------------------------------------------------------
# Main executable
# ---------------------------------------------------------------------------
add_executable(maniac src/main.F90)
set_target_properties(maniac PROPERTIES Fortran_PREPROCESS ON)
target_include_directories(maniac PRIVATE ${PETSC_INC})
target_link_libraries(maniac
    maniac_diffusion
    maniac_transport
    maniac_input_output
    maniac_solver
    maniac_iga
    maniac_fem
    maniac_core
    petsc lapack blas m dl
)