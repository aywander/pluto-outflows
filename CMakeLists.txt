cmake_minimum_required(VERSION 3.12)

project(PLUTO)

# ______________________________________________________________________________________________________________________
# Set these

set(PROBLEM                     Outflows)
set(SETUP                       mwagn)

set(PROBLEM_TYPE                Real)

# ______________________________________________________________________________________________________________________
# PLUTO settings

# -- Library paths --

set(HDF5_PATH_PAR Lib/hdf5-1.10.4-parallel)
set(HDF5_PATH_SER Lib/hdf5-1.10.4-serial)

# -- Compilers --

# Only if C compiler is not found
#set(WITH_C_COMPILER /usr/bin/gcc)

# Only if not using find_package(MPI REQUIRED)
#set(WITH_MPI_COMPILER /usr/local/bin/mpicc)
#set(WITH_MPIXX_COMPILER /usr/local/bin/mpicxx)

# -- For Chombo --

set(CHOMBO_CONFIG 3d.Darwin.64.g++.gfortran.OPT)

# ______________________________________________________________________________________________________________________
# PLUTO source code files


# NOTE: The dependences of which files to build depending on which module can be found in the individual makefiles
# in the module folders, and in AppendAdditionalFiles() in Tools/Python/define_problem.py .


# Directoiries

set(PROBLEM_DIR  ${PROBLEM_TYPE}_Problems/${PROBLEM})
set(SRC_DIR      Src)


# Included in all cases

set(INCLUDE_DIRECTORIES
        ${PROBLEM_DIR}
        ${SRC_DIR}
        ${SRC_DIR}/States
        ${SRC_DIR}/Math_Tools
        ${SRC_DIR}/Templates
        ${SRC_DIR}/Time_Stepping
        )

set(LINK_DIRECTORIES)

#set(LINK_LIBRARIES m z)
set(LINK_LIBRARIES m)


# Setup-dependent / Problem-type dependent files
# and include problem-dependent CMakeLists.txt

if(SETUP)
    set(SOURCE_FILES
            ${PROBLEM_DIR}/setup/${SETUP}/definitions.h
            ${PROBLEM_DIR}/setup/${SETUP}/init.c
            )

    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES}
            ${PROBLEM_DIR}/setup/${SETUP}
            )

    include(${PROBLEM_DIR}/setup/${SETUP}/CMakeLists.txt)

else()

    set(SOURCE_FILES
            ${PROBLEM_DIR}/definitions.h
            ${PROBLEM_DIR}/init.c
            )

    include(${PROBLEM_DIR}/CMakeLists.txt)

endif()



# Include all new files
# TODO: instead of foreach, can also do string replacement (see https://cmake.org/cmake/help/v3.0/command/string.html)

foreach(FILE ${NEW_FILES})

    set(SOURCE_FILES ${SOURCE_FILES} ${PROBLEM_DIR}/${FILE})

endforeach()


# Include all source files and overridden files

if (WITH_CHOMBO)
    set(SRC_FILES    Chombo/amrPluto.cpp)
else()
    set(SRC_FILES    main.c)
endif()


set(SRC_FILES ${SRC_FILES}
        adv_flux.c
        arrays.c
        bin_io.c
        boundary.c
        #cartcoord.c                        # Never used
        check_states.c
        cmd_line_opt.c
        colortable.c
        doxy-comment.c
        entropy_switch.c
        #failsafe.c                         # Never used
        flag_shock.c
        get_nghost.c
        globals.h                           # Modified
        initialize.c
        input_data.c
        int_bound_reset.c
        jet_domain.c
        macros.h
        mappers3D.c
        mean_mol_weight.c
        output_log.c
        parse_file.c
        parabolic_update.c
        pluto.h                             # Modified
        prototypes.h
        rbox.c
        restart.c
#        rkc.c                              # Buggy (Never used?)
        runtime_setup.c
        set_geometry.c
        set_grid.c
        set_image.c
        set_indexes.c
        set_output.c
        show_config.c
        split_source.c
        startup.c
        #startup_old.c
        structs.h                           # Modified
        tools.c
        userdef_output.c
        #userdef_output.dt.c                # Never used?
        #userdef_output.shock_tracking.c    # Never used?
        var_names.c
        vec_pot_diff.c
        #vec_pot_diff_old.c
        vec_pot_update.c
        write_data.c
        write_img.c
        write_tab.c
        write_vtk.c
        write_vtk_proc.c
        )



# Module/Physics dependent part

# Make flags options

if (WITH_PARALLEL)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPARALLEL")

    # ----

#    set(CMAKE_C_COMPILER ${WITH_MPI_COMPILER})
#    set(CMAKE_CXX_COMPILER ${WITH_MPIXX_COMPILER})

    # OR (https://cmake.org/pipermail/cmake/2015-April/060252.html)

    find_package(MPI REQUIRED)
    set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS})
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
    set(LINK_LIBRARIES ${LINK_LIBRARIES} ${MPI_LIBRARIES})

    # ----

    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Parallel)

    set(SRC_FILES  ${SRC_FILES}
            Parallel/al.h
            Parallel/al_alloc.c
            Parallel/al_boundary.c
            Parallel/al_codes.h
            Parallel/al_decomp_.c
            Parallel/al_decompose.c
            Parallel/al_defs.h
            Parallel/al_exchange.c
            Parallel/al_exchange_dim.c
            Parallel/al_finalize.c
            Parallel/al_hidden.h
            Parallel/al_init.c
            Parallel/al_io.c
            Parallel/al_proto.h
            Parallel/al_sort_.c
            Parallel/al_subarray_.c
            Parallel/al_sz_free.c
            Parallel/al_sz_get.c
            Parallel/al_sz_init.c
            Parallel/al_sz_set.c
            Parallel/al_szptr_.c
            Parallel/al_write_array_async.c
#            mpi_sub_decomp.c     # Probably shouldn't be there
            )
else()

    set(CMAKE_C_COMPILER ${WITH_C_COMPILER})
endif ()


if (WITH_PNG)
    find_package(PNG REQUIRED)

    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${PNG_INCLUDE_DIRS})
    set(LINK_LIBRARIES ${LINK_LIBRARIES} ${PNG_LIBRARIES})

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_PNG")

endif ()


if (WITH_HDF5)

    if (WITH_PARALLEL)
        set(HDF5_PREFER_PARALLEL TRUE)
        if (DEFINED HDF5_PATH_PAR)
            set(HDF5_ROOT ${HDF5_PATH_PAR})
        endif()
    else()
        if (DEFINED HDF5_PATH_SER)
            set(HDF5_ROOT ${HDF5_PATH_SER})
        endif ()
    endif ()

    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    find_package(HDF5 REQUIRED)


    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${HDF5_INCLUDE_DIRS})
    set(LINK_LIBRARIES ${LINK_LIBRARIES} ${HDF5_LIBRARIES})

    #message(STATUS "HDF_LIBRARIES = ${HDF5_LIBRARIES}")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_HDF5")
    set(SRC_FILES  ${SRC_FILES}
            hdf5_io.c
            )
endif ()


if (WITH_ASYNC_IO)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_ASYNC_IO")

endif ()


if (WITH_CHOMBO)

    # Parse the CHOMBO_CONFIG variable to set some parameters

    string(SUBSTRING ${CHOMBO_CONFIG} 0 1 CHOMBO_DIM)

    string(REGEX MATCH MPI CHOMBO_MPI ${CHOMBO_CONFIG})
    if (CHOMBO_MPI STREQUAL MPI)
        set(CHOMBO_PARALLEL ON)
    else()
        set(CHOMBO_PARALLEL OFF)
    endif()

    string(REGEX MATCH DEBUG CHOMBO_DEBUG ${CHOMBO_CONFIG})
    if (CHOMBO_DEBUG STREQUAL DEBUG)
        set(CHOMBO_DEBUG ON)
    else()
        set(CHOMBO_DEBUG OFF)
    endif()


    # Add definitions

    add_definitions(-DCH_LANG_CC)
    add_definitions(-DUSE_HDF5)
    add_definitions(-DCH_SPACEDIM=${CHOMBO_DIM})
    if (CHOMBO_PARALLEL)
        add_definitions(-DCH_MPI)
    endif()


    # Libraries include

    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Chombo)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} Lib/Chombo-3.2/lib/include)

    # The following is in Makefile.chombo, but not needed. Just keeping for reference.
    # set(CHOMBO_LIBNAMES AMRTimeDependent AMRTools BoxTools BaseTools)
    #foreach(LIBNAME ${CHOMBO_LIBNAMES})
    #    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} Lib/Chombo-3.2/lib/src/${LIBNAME})
    #endforeach()


    # Libraries link directory

    set(LINK_DIRECTORIES ${LINK_DIRECTORIES} Lib/Chombo-3.2/lib)


    # Link libraries

    # TODO: just use string replacement rather than loop (see https://cmake.org/cmake/help/v3.0/command/string.html)
    set(CHOMBO_LIBNAMES AMRTimeDependent AMRTools BoxTools BaseTools)
    foreach(LIBNAME ${CHOMBO_LIBNAMES})
        string(TOLOWER ${LIBNAME} LIBNAME)
        set(LINK_LIBRARIES ${LINK_LIBRARIES} ${LIBNAME}${CHOMBO_CONFIG})
    endforeach()


    # PLUTO CHOMBO source files

    set(SRC_FILES  ${SRC_FILES}
            Chombo/AMRLevelPluto.cpp
            Chombo/AMRLevelPluto.H
            Chombo/AMRLevelPlutoFactory.cpp
            Chombo/AMRLevelPlutoFactory.H
            Chombo/LevelPluto.cpp
            Chombo/LevelPluto.H
#            Chombo/PatchCTU.cpp   # Not used
            Chombo/PatchEuler.cpp
            Chombo/PatchGrid.H
            Chombo/PatchPluto.cpp
            Chombo/PatchPluto.H
            Chombo/PatchStartup.cpp
            Chombo/PatchTools.cpp
            Chombo/TagCells.cpp
            )

endif ()


# Physics

if (WITH_PHYSICS STREQUAL HD)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/HD ${SRC_DIR}/MHD)
    set(SRC_FILES  ${SRC_FILES}
            HD/ausm.c
#            HD/ausm_up.c           # Buggy (don't need)
            HD/eigenv.c
            HD/fluxes.c
#            MHD/get_current.c      # Don't need?
            HD/hll.c
            HD/hll_speed.c
            HD/hllc.c
            HD/mappers.c
            HD/mod_defs.h
            HD/prim_eqn.c            # Only used for ChTr / Hancock (?)
            MHD/rhs.c
            MHD/rhs_source.c
#            HD/riemann_full.c       # Buggy (not used?)
            HD/roe.c
#            HD/rusanov-dw.c         # Don't need
            HD/set_solver.c
            HD/tvdlf.c
            HD/two_shock.c
            )


elseif (WITH_PHYSICS STREQUAL RHD)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/RHD ${SRC_DIR}/RMHD)
    set(SRC_FILES  ${SRC_FILES}
            RHD/eigenv.c
            RHD/energy_solve.c
            RHD/entropy_solve.c
            RHD/fluxes.c
            RMHD/four_vel.c
            RHD/hll.c
            RHD/hll_speed.c
            RHD/hllc.c
            RHD/mappers.c
            RHD/mod_defs.h
            RHD/pressure_fix.c
            RHD/prim_eqn.c
            RMHD/rhs.c
            RHD/set_solver.c
            RHD/tvdlf.c
            RHD/two_shock.c
            RMHD/vel_limiter.c
            )


elseif (WITH_PHYSICS STREQUAL MHD)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/MHD)
    set(SRC_FILES  ${SRC_FILES}
            MHD/bckgrnd_field.c
            MHD/eigenv.c
            MHD/fluxes.c
            MHD/get_current.c
            MHD/hll.c
            MHD/hll_speed.c
            MHD/hllc.c
            MHD/hlld.c
            MHD/mappers.c
            MHD/mod_defs.h
            MHD/prim_eqn.c
            MHD/rhs.c
            MHD/rhs_source.c
#            MHD/roe.bckfld.c          # Don't need
            MHD/roe.c
            MHD/set_solver.c
            MHD/source.c
            MHD/tvdlf.c
            )

elseif (WITH_PHYSICS STREQUAL RMHD)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/RMHD)
    set(SRC_FILES  ${SRC_FILES}
            RMHD/eigenv.c
            RMHD/energy_solve.c        # Modified: some bugs for TAUB
            RMHD/entropy_solve.c
            RMHD/fluxes.c
            RMHD/gmusta1.c
            RMHD/four_vel.c
            RMHD/hll.c
            RMHD/hll_speed.c
            RMHD/hllc.c
            RMHD/hlld.c
            RMHD/mappers.c
            RMHD/mod_defs.h
            RMHD/pressure_fix.c         # Modified: some bugs for TAUB
            RMHD/prim_eqn.c
            RMHD/quartic.c
            RMHD/rhs.c
#            RMHD/rusanov-dw.c           # Not used
            RMHD/set_solver.c
            RMHD/source.c
            RMHD/tvdlf.c
            RMHD/vel_limiter.c
            )

endif ()


if (WITH_THERMAL_CONDUTION)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Thermal_Conduction)
    set(SRC_FILES  ${SRC_FILES}
            sts.c
            rkl.c
            Thermal_Conduction/tc.h
            Thermal_Conduction/tc_flux.c
            Thermal_Conduction/tc_functions.c
            Thermal_Conduction/tc_kappa.c
            Thermal_Conduction/tc_rhs.c
            )
endif ()

if (WITH_VISCOSITY)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Viscosity)
    set(SRC_FILES  ${SRC_FILES}
            sts.c
            rkl.c
            Viscosity/visc_nu.c
            Viscosity/viscosity.h
            Viscosity/viscous_flux.c
            Viscosity/viscous_rhs.c
            )
endif ()

if (WITH_SHEARING_BOX)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/MHD/ShearingBox)
    set(SRC_FILES  ${SRC_FILES}
            MHD/ShearingBox/sb_boundary.c
            MHD/ShearingBox/sb_flux.c
            MHD/ShearingBox/sb_tools.c
            MHD/ShearingBox/shearingbox.h
            )
endif ()

if (WITH_RESISTIVITY)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/MHD/Resistive)
    set(SRC_FILES  ${SRC_FILES}
            sts.c
            rkl.c
            MHD/Resistive/res.h
            MHD/Resistive/res_eta.c
            MHD/Resistive/res_flux.c
            MHD/Resistive/res_functions.c
            MHD/Resistive/res_rhs.c
            )
endif ()

if (WITH_CT)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/MHD/CT)
    set(SRC_FILES  ${SRC_FILES}
            MHD/CT/ct.h
            MHD/CT/ct_emf.c
            MHD/CT/ct_emf_average.c
            MHD/CT/ct_field_average.c
            MHD/CT/ct_fill_mag_field.c
            MHD/CT/ct_stag_slopes.c
            MHD/CT/ct_update.c
            MHD/CT/emf_boundary.c
            MHD/CT/emf_userdef.c
            )
endif ()

if (WITH_GLM)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/MHD/GLM)
    set(SRC_FILES  ${SRC_FILES}
            MHD/GLM/analysis.c
            MHD/GLM/glm.c
            MHD/GLM/glm.h
            )
endif ()

if (WITH_HALL_MHD)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/MHD/Hall_MHD)
    set(SRC_FILES  ${SRC_FILES}
            MHD/Hall_MHD/hall_flux.c
            MHD/Hall_MHD/hall_mhd.h
            MHD/Hall_MHD/hall_ne.c
            MHD/Hall_MHD/hall_whistler_speed.c
            MHD/Hall_MHD/hall_whistler_speed.c
            )
endif ()


# Cooling options

if (WITH_COOLING STREQUAL H2_COOL)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Cooling/H2_COOL)
    set(SRC_FILES  ${SRC_FILES}
            Cooling/H2_COOL/comp_equil.c
            Cooling/H2_COOL/cooling.h
            Cooling/H2_COOL/jacobian.c
            Cooling/H2_COOL/maxrate.c
            Cooling/H2_COOL/radiat.c
            )

elseif (WITH_COOLING STREQUAL TABULATED)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Cooling/TABULATED)
    set(SRC_FILES ${SRC_FILES}
            cooling_ode_solver.c
            cooling_source.c
            Cooling/TABULATED/cooling.h
#            Cooling/TABULATED/extract.c ???
            Cooling/TABULATED/jacobian.c
            Cooling/TABULATED/maxrate.c
            Cooling/TABULATED/radiat.c
            )

elseif (WITH_COOLING STREQUAL MINEQ)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Cooling/MINEq)
    set(SRC_FILES  ${SRC_FILES}
            Cooling/MINEq/comp_equil.c
            Cooling/MINEq/cooling.h
            Cooling/MINEq/cooling_defs.h
            Cooling/MINEq/ion_init.c
            Cooling/MINEq/jacobian.c
            Cooling/MINEq/make_tables.c
            Cooling/MINEq/maxrate.c
            Cooling/MINEq/radiat.c
            Cooling/MINEq/test.c
            )


elseif (WITH_COOLING STREQUAL POWERLAW)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Cooling/Power_Law)
    set(SRC_FILES  ${SRC_FILES}
            Cooling/Power_Law/cooling.c
            Cooling/Power_Law/cooling.h
            )

elseif (WITH_COOLING STREQUAL SNEQ)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Cooling/SNEq)
    set(SRC_FILES  ${SRC_FILES}
            Cooling/SNEq/cooling.h
            Cooling/SNEq/jacobian.c
            Cooling/SNEq/maxrate.c
            Cooling/SNEq/radiat.c
            Cooling/SNEq/radiat.modified_by_Bhargav.c
            )

endif ()


# Particles

if (WITH_PARTICLES)

    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Particles)
    set(SRC_FILES  ${SRC_FILES}
            Particles/particles.h
            Particles/particles_boundary.c
            Particles/particles_deposit.c
            Particles/particles_distrib_regular.c
            Particles/particles_init.c
            Particles/particles_load.c
            Particles/particles_mpi_datatype.c
            Particles/particles_restart.c
            Particles/particles_set.c
            Particles/particles_set_output.c
            Particles/particles_tools.c
            Particles/particles_write_bin.c
            Particles/particles_write_data.c
            Particles/particles_write_trajectory.c
            Particles/particles_write_vtk.c
            Particles/plist_tools.c
            )

if (WITH_PARTICLES STREQUAL LAGRANGIAN)

    set(SRC_FILES  ${SRC_FILES}
            Particles/particles_lp_tools.c
            Particles/particles_lp_update.c
            )

elseif (WITH_PARTICLES STREQUAL COSMIC_RAYS)

    set(SRC_FILES  ${SRC_FILES}
            Particles/particles_cr_feedback.c
            Particles/particles_cr_force.c
            Particles/particles_cr_update.c
            Particles/particles_cr_update_test.c
            )

elseif (WITH_PARTICLES STREQUAL DUST)

    set(SRC_FILES  ${SRC_FILES}
            )

endif ()

endif ()


# Dust Fluid

if (WITH_DUST_FLUID)

    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Dust)
    set(SRC_FILES  ${SRC_FILES}
            Dust/dust.h
            Dust/dust.c
            )
endif ()


# Forced turbulence

if (WITH_FORCED_TURB)

    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Forced_Turb)
    set(SRC_FILES  ${SRC_FILES}
            Forced_Turb/forced_turb.h
            Forced_Turb/forced_turb.c
            )

endif ()


# EOS

if (WITH_EOS STREQUAL IDEAL)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/EOS/Ideal)
    set(SRC_FILES  ${SRC_FILES}
            EOS/Ideal/eos.c
            EOS/Ideal/eos.h
            )


elseif (WITH_EOS STREQUAL ISOTHERMAL)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/EOS/Isothermal)
    set(SRC_FILES  ${SRC_FILES}
            EOS/Isothermal/eos.c
            EOS/Isothermal/eos.h
            )


elseif (WITH_EOS STREQUAL PVTE)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/EOS/PVTE)
    set(SRC_FILES  ${SRC_FILES}
            EOS/PVTE/eos.h
            EOS/PVTE/fundamental_derivative.c
            EOS/PVTE/internal_energy.c
            EOS/PVTE/pvte_law.c
            EOS/PVTE/pvte_law_dAngelo.c
            EOS/PVTE/pvte_law_H+.c
            EOS/PVTE/pvte_law_template.c
            EOS/PVTE/scvh.c
            EOS/PVTE/thermal_eos.c
            EOS/PVTE/zeta_tables.c
            )

elseif (WITH_EOS STREQUAL TAUB)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/EOS/Taub)
    set(SRC_FILES  ${SRC_FILES}
            EOS/Taub/eos.c
            EOS/Taub/eos.h
            )

endif ()


# Time stepping

set(SRC_FILES  ${SRC_FILES}
        Time_Stepping/update_stage.c
        )

if ((WITH_TIME_STEPPING STREQUAL RK2) OR (WITH_TIME_STEPPING STREQUAL RK3))
    set(SRC_FILES  ${SRC_FILES}
            Time_Stepping/rk_step.c
            )

elseif (WITH_TIME_STEPPING STREQUAL CHARACTERISTIC_TRACING)
    set(SRC_FILES  ${SRC_FILES}
            States/char_tracing.c
            Time_Stepping/ctu_step.c
            )

elseif (WITH_TIME_STEPPING STREQUAL HANCOCK)
    set(SRC_FILES  ${SRC_FILES}
            States/hancock.c
            Time_Stepping/ctu_step.c
            )
endif()


# Interpolation

set(SRC_FILES  ${SRC_FILES}
        States/flatten.c
        )

if (WITH_INTERPOLATION STREQUAL FLAT)
    set(SRC_FILES  ${SRC_FILES}
            States/flat_states.c
            )

elseif (WITH_INTERPOLATION STREQUAL LINEAR)
    set(SRC_FILES  ${SRC_FILES}
            States/plm_coeffs.c
            States/plm_coeffs.h
            States/plm_states.c
            )

elseif (WITH_INTERPOLATION STREQUAL PARABOLIC)
    set(SRC_FILES  ${SRC_FILES}
            States/plm_coeffs.c
            States/plm_coeffs.h
            States/ppm_coeffs.c
            States/ppm_coeffs.h
            States/ppm_states.c
            )

elseif (WITH_INTERPOLATION STREQUAL WENO3)
    set(SRC_FILES  ${SRC_FILES}
            States/weno3_states.c
            )

elseif (WITH_INTERPOLATION STREQUAL LIMO3)
    set(SRC_FILES  ${SRC_FILES}
            States/limo3_states.c
            )

endif ()


# Other modules

if (WITH_FARGO)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SRC_DIR}/Fargo)
    set(SRC_FILES  ${SRC_FILES}
            Fargo/fargo.c
            Fargo/fargo.h
            Fargo/fargo_source.c
            Fargo/fargo_velocity.c
            )
endif ()

if (WITH_FD)
    set(SRC_FILES  ${SRC_FILES}
            fd_flux.c
            fd_reconstruct.c
            States/fd_states.c
            )
endif ()


# Math tools

set(SRC_FILES  ${SRC_FILES}
        Math_Tools/math_interp.c
        Math_Tools/math_lu_decomp.c
        Math_Tools/math_misc.c
        Math_Tools/math_ode.c
        Math_Tools/math_qr_decomp.c
        Math_Tools/math_quadrature.c
        Math_Tools/math_rand_mt19937-64.c
        Math_Tools/math_random.c
        Math_Tools/math_root_finders.c
        Math_Tools/math_table2D.c
        Math_Tools/math_tools.h
        )


# Override SRC files with files in problem directory

foreach (FILE ${SRC_FILES})

    if (${FILE} IN_LIST OVERRIDE_FILES)

        set(SOURCE_FILES  ${SOURCE_FILES}  ${PROBLEM_DIR}/${FILE})

    elseif (setup/${SETUP}/${FILE} IN_LIST OVERRIDE_FILES)

        set(SOURCE_FILES  ${SOURCE_FILES}  ${PROBLEM_DIR}/setup/${SETUP}/${FILE})

    else()

        set(SOURCE_FILES  ${SOURCE_FILES}  ${SRC_DIR}/${FILE})

    endif()

endforeach()


# Other flags


# Include, link, and executable

include_directories(${INCLUDE_DIRECTORIES})

link_directories(${LINK_DIRECTORIES})

link_libraries(${LINK_LIBRARIES})

add_executable(pluto ${SOURCE_FILES})


