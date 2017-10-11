cmake_minimum_required(VERSION 3.4)

project(PLUTO)

# ______________________________________________________________________________________________________________________
# Set these

set(PROBLEM                     DynamicalFriction)
set(PROBLEM_TYPE                Real)


# ______________________________________________________________________________________________________________________
# PLUTO settings

set(PLUTO_SOURCE_DIR            /Users/ayw/dev/PLUTO)

# -- Library paths --

set(HDF5_PATH_PAR ${PLUTO_SOURCE_DIR}/Lib/hdf5-parallel)
set(HDF5_PATH_SER ${PLUTO_SOURCE_DIR}/Lib/hdf5-serial)

# -- Compilers --

set(WITH_C_COMPILER /usr/bin/gcc)
set(WITH_MPI_COMPILER /usr/local/bin/mpicc)
set(WITH_MPIXX_COMPILER /usr/local/bin/mpicxx)


# ______________________________________________________________________________________________________________________
# PLUTO source code files


# NOTE: The dependences of which files to build depending on which module can be found in the individual makefiles
# in the module folders, and in AppendAdditionalFiles() in Tools/Python/define_problem.py .


# Directoiries

set(PROBLEM_DIR  ${PLUTO_SOURCE_DIR}/${PROBLEM_TYPE}_Problems/${PROBLEM})
set(SOURCE_DIR   ${PLUTO_SOURCE_DIR}/Src)


# Included in all cases

set(INCLUDE_DIRECTORIES
        ${PROBLEM_DIR}
        ${SOURCE_DIR}
        ${SOURCE_DIR}/States
        ${SOURCE_DIR}/Math_Tools
        ${SOURCE_DIR}/Templates
        ${SOURCE_DIR}/Time_Stepping
        )

set(LINK_DIRECTORIES)

set(LINK_LIBRARIES m z)


# Setup-dependent / Problem-type dependent files

set(SOURCE_FILES
        ${PROBLEM_DIR}/definitions.h
        ${PROBLEM_DIR}/init.c
        )


# Include problem-dependent CMakeLists.txt

include(${PROBLEM_DIR}/CMakeLists.txt)


# Include all new files

foreach(FILE ${NEW_FILES})

    set(SOURCE_FILES ${SOURCE_FILES} ${PROBLEM_DIR}/${FILE})

endforeach()


# Include all source files and overrideen files

set(SRC_FILES
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
        main.c
        mappers3D.c
        mean_mol_weight.c
        parabolic_flux.c
        parabolic_rhs.c
        parse_file.c
        pluto.h                             # Modified
        prototypes.h
        rbox.c
        restart.c
        #rkc.c                              # Buggy (Never used?)
        runtime_setup.c
        set_geometry.c
        set_grid.c
        set_image.c
        set_indexes.c
        set_output.c
        show_config.c
        split_source.c
        startup.c
        structs.h                           # Modified
        sts.c
        tools.c
        userdef_output.c
        #userdef_output.dt.c                # Never used?
        #userdef_output.shock_tracking.c    # Never used?
        var_names.c
        vec_pot_diff.c
        vec_pot_update.c
        write_data.c
        write_img.c
        write_tab.c
        write_vtk.c
        )

foreach(FILE ${SRC_FILES})

    if(FILE IN_LIST OVERRIDE_FILES)

        set(SOURCE_FILES  ${SOURCE_FILES}  ${PROBLEM_DIR}/${FILE})

    else()

        set(SOURCE_FILES  ${SOURCE_FILES}  ${SOURCE_DIR}/${FILE})

    endif()

endforeach()


# Module/Physics dependent part

# Make flags options

if (WITH_PARALLEL)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPARALLEL")

    # ----

    set(CMAKE_C_COMPILER ${WITH_MPI_COMPILER})
    set(CMAKE_CXX_COMPILER ${WITH_MPIXX_COMPILER})

    # OR (https://cmake.org/pipermail/cmake/2015-April/060252.html)

    #find_package(MPI REQUIRED)
    #set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    #set(CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS})
    #set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    #set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
    #set(LINK_LIBRARIES ${LINK_LIBRARIES} ${MPI_LIBRARIES})

    # ----

    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Parallel)

    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/Parallel/al.h
            ${SOURCE_DIR}/Parallel/al_alloc.c
            ${SOURCE_DIR}/Parallel/al_boundary.c
            ${SOURCE_DIR}/Parallel/al_codes.h
            ${SOURCE_DIR}/Parallel/al_decomp_.c
            ${SOURCE_DIR}/Parallel/al_decompose.c
            ${SOURCE_DIR}/Parallel/al_defs.h
            ${SOURCE_DIR}/Parallel/al_exchange.c
            ${SOURCE_DIR}/Parallel/al_exchange_dim.c
            ${SOURCE_DIR}/Parallel/al_finalize.c
            ${SOURCE_DIR}/Parallel/al_hidden.h
            ${SOURCE_DIR}/Parallel/al_init.c
            ${SOURCE_DIR}/Parallel/al_io.c
            ${SOURCE_DIR}/Parallel/al_proto.h
            ${SOURCE_DIR}/Parallel/al_sort_.c
            ${SOURCE_DIR}/Parallel/al_subarray_.c
            ${SOURCE_DIR}/Parallel/al_sz_free.c
            ${SOURCE_DIR}/Parallel/al_sz_get.c
            ${SOURCE_DIR}/Parallel/al_sz_init.c
            ${SOURCE_DIR}/Parallel/al_sz_set.c
            ${SOURCE_DIR}/Parallel/al_szptr_.c
            ${SOURCE_DIR}/Parallel/al_write_array_async.c
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
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/hdf5_io.c
            )
endif ()


if (WITH_ASYNC_IO)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_ASYNC_IO")

endif ()


if (WITH_CHOMBO)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_HDF5")
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Chombo)
    set(LINK_DIRECTORIES ${LINK_DIRECTORIES} ${PLUTO_SOURCE_DIR}/Lib/Chombo-3.2/lib)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/hdf5_io.vChombo3.1.c
            ${SOURCE_DIR}/hdf5_io.c
            ${SOURCE_DIR}/Chombo/AMRLevelPluto.cpp
            ${SOURCE_DIR}/Chombo/AMRLevelPluto.H
            ${SOURCE_DIR}/Chombo/AMRLevelPlutoFactory.cpp
            ${SOURCE_DIR}/Chombo/AMRLevelPlutoFactory.H
            ${SOURCE_DIR}/Chombo/amrPluto.cpp
            ${SOURCE_DIR}/Chombo/LevelPluto.cpp
            ${SOURCE_DIR}/Chombo/LevelPluto.H
            ${SOURCE_DIR}/Chombo/PatchCTU.cpp
            ${SOURCE_DIR}/Chombo/PatchEuler.cpp
            ${SOURCE_DIR}/Chombo/PatchGrid.H
            ${SOURCE_DIR}/Chombo/PatchPluto.cpp
            ${SOURCE_DIR}/Chombo/PatchPluto.H
            ${SOURCE_DIR}/Chombo/PatchStartup.cpp
            ${SOURCE_DIR}/Chombo/PatchTools.cpp
            ${SOURCE_DIR}/Chombo/TagCells.cpp
            )
endif ()


# Physics

if (WITH_PHYSICS STREQUAL HD)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/HD ${SOURCE_DIR}/MHD)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/HD/ausm.c
            #${SOURCE_DIR}/HD/ausm_up.c           # Buggy (don't need)
            ${SOURCE_DIR}/HD/eigenv.c
            ${SOURCE_DIR}/HD/fluxes.c
            ${SOURCE_DIR}/MHD/get_area_flux.c
            ${SOURCE_DIR}/HD/hll.c
            ${SOURCE_DIR}/HD/hll_speed.c
            ${SOURCE_DIR}/HD/hllc.c
            ${SOURCE_DIR}/HD/mappers.c
            ${SOURCE_DIR}/HD/mod_defs.h
            ${SOURCE_DIR}/HD/prim_eqn.c            # Only used for ChTr / Hancock (?)
            ${SOURCE_DIR}/MHD/rhs.c
            ${SOURCE_DIR}/MHD/rhs_source.c
            #${SOURCE_DIR}/HD/riemann_full.c       # Buggy (not used?)
            ${SOURCE_DIR}/HD/roe.c
            #${SOURCE_DIR}/HD/rusanov-dw.c         # Don't need
            ${SOURCE_DIR}/HD/set_solver.c
            ${SOURCE_DIR}/HD/tvdlf.c
            ${SOURCE_DIR}/HD/two_shock.c
            )


elseif (WITH_PHYSICS STREQUAL RHD)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/RHD ${SOURCE_DIR}/RMHD)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/RHD/eigenv.c
            ${SOURCE_DIR}/RHD/energy_solve.c
            ${SOURCE_DIR}/RHD/entropy_solve.c
            ${SOURCE_DIR}/RHD/fluxes.c
            ${SOURCE_DIR}/RMHD/four_vel.c
            ${SOURCE_DIR}/RHD/hll.c
            ${SOURCE_DIR}/RHD/hll_speed.c
            ${SOURCE_DIR}/RHD/hllc.c
            ${SOURCE_DIR}/RHD/mappers.c
            ${SOURCE_DIR}/RHD/mod_defs.h
            ${SOURCE_DIR}/RHD/pressure_fix.c
            ${SOURCE_DIR}/RHD/prim_eqn.c
            ${SOURCE_DIR}/RMHD/rhs.c
            ${SOURCE_DIR}/RHD/set_solver.c
            ${SOURCE_DIR}/RHD/tvdlf.c
            ${SOURCE_DIR}/RHD/two_shock.c
            ${SOURCE_DIR}/RMHD/vel_limiter.c
            )


elseif (WITH_PHYSICS STREQUAL MHD)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/MHD)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/MHD/bckgrnd_field.c
            ${SOURCE_DIR}/MHD/eigenv.c
            ${SOURCE_DIR}/MHD/fluxes.c
            ${SOURCE_DIR}/MHD/get_area_flux.c
            ${SOURCE_DIR}/MHD/hll.c
            ${SOURCE_DIR}/MHD/hll_speed.c
            ${SOURCE_DIR}/MHD/hllc.c
            ${SOURCE_DIR}/MHD/hlld.c
            ${SOURCE_DIR}/MHD/mappers.c
            ${SOURCE_DIR}/MHD/mod_defs.h
            ${SOURCE_DIR}/MHD/prim_eqn.c
            ${SOURCE_DIR}/MHD/rhs.c
            ${SOURCE_DIR}/MHD/rhs_source.c
            #${SOURCE_DIR}/MHD/roe.bckfld.c          # Don't need
            ${SOURCE_DIR}/MHD/roe.c
            ${SOURCE_DIR}/MHD/set_solver.c
            ${SOURCE_DIR}/MHD/source.c
            ${SOURCE_DIR}/MHD/tvdlf.c
            )

elseif (WITH_PHYSICS STREQUAL RMHD)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/RMHD)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/RMHD/eigenv.c
            ${SOURCE_DIR}/RMHD/energy_solve.c
            ${SOURCE_DIR}/RMHD/entropy_solve.c
            ${SOURCE_DIR}/RMHD/fluxes.c
            ${SOURCE_DIR}/RMHD/four_vel.c
            ${SOURCE_DIR}/RMHD/hll.c
            ${SOURCE_DIR}/RMHD/hll_speed.c
            ${SOURCE_DIR}/RMHD/hllc.c
            ${SOURCE_DIR}/RMHD/hlld.c
            ${SOURCE_DIR}/RMHD/mappers.c
            ${SOURCE_DIR}/RMHD/mod_defs.h
            ${SOURCE_DIR}/RMHD/pressure_fix.c
            ${SOURCE_DIR}/RMHD/prim_eqn.c
            ${SOURCE_DIR}/RMHD/quartic.c
            #${SOURCE_DIR}/RMHD/quartic.devel1.c       # Not used
            #${SOURCE_DIR}/RMHD/quartic.new.c          # Not used
            #${SOURCE_DIR}/RMHD/quartic.old.c          # Not used
            ${SOURCE_DIR}/RMHD/rhs.c
            #${SOURCE_DIR}/RMHD/rusanov-dw.c           # Not used
            ${SOURCE_DIR}/RMHD/set_solver.c
            ${SOURCE_DIR}/RMHD/source.c
            ${SOURCE_DIR}/RMHD/tvdlf.c
            ${SOURCE_DIR}/RMHD/vel_limiter.c
            )

endif ()


if (WITH_THERMAL_CONDUTION)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Thermal_Conduction)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/Thermal_Conduction/tc.h
            ${SOURCE_DIR}/Thermal_Conduction/tc_flux.c
            ${SOURCE_DIR}/Thermal_Conduction/tc_functions.c
            ${SOURCE_DIR}/Thermal_Conduction/tc_kappa.c
            )
endif ()

if (WITH_VISCOSITY)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Viscosity)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/Viscosity/visc_nu.c
            ${SOURCE_DIR}/Viscosity/viscosity.h
            ${SOURCE_DIR}/Viscosity/viscous_flux.c
            )
endif ()

if (WITH_SHEARING_BOX)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/MHD/ShearingBox)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/MHD/ShearingBox/sb_boundary.c
            ${SOURCE_DIR}/MHD/ShearingBox/sb_flux.c
            ${SOURCE_DIR}/MHD/ShearingBox/sb_tools.c
            ${SOURCE_DIR}/MHD/ShearingBox/shearingbox.h
            )
endif ()

if (WITH_RESISTIVITY)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/MHD/Resistive)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/MHD/Resistive/res.h
            ${SOURCE_DIR}/MHD/Resistive/res_eta.c
            ${SOURCE_DIR}/MHD/Resistive/res_flux.c
            ${SOURCE_DIR}/MHD/Resistive/res_functions.c
            )
endif ()

if (WITH_CT)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/MHD/CT)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/MHD/CT/ct.c
            ${SOURCE_DIR}/MHD/CT/ct.h
            ${SOURCE_DIR}/MHD/CT/ct_emf.c
            ${SOURCE_DIR}/MHD/CT/ct_emf_average.c
            ${SOURCE_DIR}/MHD/CT/ct_field_average.c
            ${SOURCE_DIR}/MHD/CT/ct_fill_mag_field.c
            ${SOURCE_DIR}/MHD/CT/ct_stag_slopes.c
            ${SOURCE_DIR}/MHD/CT/emf_boundary.c
            ${SOURCE_DIR}/MHD/CT/emf_userdef.c
            )
endif ()

if (WITH_GLM)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/MHD/GLM)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/MHD/GLM/analysis.c
            ${SOURCE_DIR}/MHD/GLM/glm.c
            ${SOURCE_DIR}/MHD/GLM/glm.h
            )
endif ()


# Cooling options

if (WITH_COOLING STREQUAL H2_COOL)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Cooling/H2_COOL)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/Cooling/H2_COOL/comp_equil.c
            ${SOURCE_DIR}/Cooling/H2_COOL/cooling.h
            ${SOURCE_DIR}/Cooling/H2_COOL/jacobian.c
            ${SOURCE_DIR}/Cooling/H2_COOL/maxrate.c
            ${SOURCE_DIR}/Cooling/H2_COOL/pluto.ini
            ${SOURCE_DIR}/Cooling/H2_COOL/radiat.c
            )

elseif (WITH_COOLING STREQUAL TABULATED)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Cooling/TABULATED)
    set(SOURCE_FILES ${SOURCE_FILES}
            ${PROBLEM_DIR}/cooling_source.c
            ${PROBLEM_DIR}/radiat.c
            ${SOURCE_DIR}/cooling_ode_solver.c
            ${SOURCE_DIR}/Cooling/TABULATED/cooling.h
            #${SOURCE_DIR}/Cooling/TABULATED/extract.c
            ${SOURCE_DIR}/Cooling/TABULATED/jacobian.c
            ${SOURCE_DIR}/Cooling/TABULATED/maxrate.c
            #${SOURCE_DIR}/Cooling/TABULATED/radiat.c           # Override
            )

elseif (WITH_COOLING STREQUAL MINEQ)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Cooling/MINEq)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/Cooling/MINEq/comp_equil.c
            ${SOURCE_DIR}/Cooling/MINEq/cooling.h
            ${SOURCE_DIR}/Cooling/MINEq/cooling_defs.h
            ${SOURCE_DIR}/Cooling/MINEq/ion_init.c
            ${SOURCE_DIR}/Cooling/MINEq/jacobian.c
            ${SOURCE_DIR}/Cooling/MINEq/make_tables.c
            ${SOURCE_DIR}/Cooling/MINEq/maxrate.c
            ${SOURCE_DIR}/Cooling/MINEq/radiat.c
            ${SOURCE_DIR}/Cooling/MINEq/test.c
            )


elseif (WITH_COOLING STREQUAL POWERLAW)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Cooling/Power_Law)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/Cooling/Power_Law/cooling.c
            ${SOURCE_DIR}/Cooling/Power_Law/cooling.h
            )

elseif (WITH_COOLING STREQUAL SNEQ)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Cooling/SNEq)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/Cooling/SNEq/cooling.h
            ${SOURCE_DIR}/Cooling/SNEq/jacobian.c
            ${SOURCE_DIR}/Cooling/SNEq/maxrate.c
            ${SOURCE_DIR}/Cooling/SNEq/radiat.c
            ${SOURCE_DIR}/Cooling/SNEq/radiat.modified_by_Bhargav.c
            )

else ()
    if (PROBLEM_TYPE STREQUAL Real)
        set(SOURCE_FILES ${SOURCE_FILES}
                )
    endif()

endif ()


# EOS

if (WITH_EOS STREQUAL IDEAL)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/EOS/Ideal)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/EOS/Ideal/eos.c
            ${SOURCE_DIR}/EOS/Ideal/eos.h
            )


elseif (WITH_EOS STREQUAL ISOTHERMAL)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/EOS/Isothermal)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/EOS/Isothermal/eos.c
            ${SOURCE_DIR}/EOS/Isothermal/eos.h
            )


elseif (WITH_EOS STREQUAL PVTE)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/EOS/PVTE)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/EOS/PVTE/eos.h
            ${SOURCE_DIR}/EOS/PVTE/fundamental_derivative.c
            ${SOURCE_DIR}/EOS/PVTE/internal_energy.c
            ${SOURCE_DIR}/EOS/PVTE/pvte_law.c
            ${SOURCE_DIR}/EOS/PVTE/pvte_law_dAngelo.c
            ${SOURCE_DIR}/EOS/PVTE/pvte_law_H+.c
            ${SOURCE_DIR}/EOS/PVTE/pvte_law_template.c
            ${SOURCE_DIR}/EOS/PVTE/scvh.c
            ${SOURCE_DIR}/EOS/PVTE/thermal_eos.c
            ${SOURCE_DIR}/EOS/PVTE/zeta_tables.c
            )

elseif (WITH_EOS STREQUAL TAUB)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/EOS/Taub)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/EOS/Taub/eos.c
            ${SOURCE_DIR}/EOS/Taub/eos.h
            )

endif ()


# Time stepping

set(SOURCE_FILES  ${SOURCE_FILES}
        ${SOURCE_DIR}/Time_Stepping/update_stage.c
        )

if ((WITH_TIME_STEPPING STREQUAL RK2) OR (WITH_TIME_STEPPING STREQUAL RK3))
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/Time_Stepping/rk_step.c
            )

elseif (WITH_TIME_STEPPING STREQUAL CHARACTERISTIC_TRACING)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/States/char_tracing.c
            ${SOURCE_DIR}/Time_Stepping/ctu_step.c
            )

elseif (WITH_TIME_STEPPING STREQUAL HANCOCK)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/States/hancock.c
            ${SOURCE_DIR}/Time_Stepping/ctu_step.c
            )
endif()


# Interpolation

set(SOURCE_FILES  ${SOURCE_FILES}
        ${SOURCE_DIR}/States/flatten.c
        )

if (WITH_INTERPOLATION STREQUAL FLAT)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/States/flat_states.c
            )

elseif (WITH_INTERPOLATION STREQUAL LINEAR)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/States/plm_coeffs.c
            ${SOURCE_DIR}/States/plm_coeffs.h
            ${SOURCE_DIR}/States/plm_states.c
            )

elseif (WITH_INTERPOLATION STREQUAL PARABOLIC)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/States/plm_coeffs.c
            ${SOURCE_DIR}/States/plm_coeffs.h
            ${SOURCE_DIR}/States/ppm_coeffs.c
            ${SOURCE_DIR}/States/ppm_coeffs.h
            ${SOURCE_DIR}/States/ppm_states.c
            )

elseif (WITH_INTERPOLATION STREQUAL WENO3)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/States/weno3_states.c
            )

elseif (WITH_INTERPOLATION STREQUAL LIMO3)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/States/limo3_states.c
            )

endif ()


# Other modules

if (WITH_FARGO)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Fargo)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/Fargo/fargo.c
            ${SOURCE_DIR}/Fargo/fargo.h
            ${SOURCE_DIR}/Fargo/fargo_source.c
            ${SOURCE_DIR}/Fargo/fargo_velocity.c
            )
endif ()

if (WITH_PARTICLES)
    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Particles)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/Particles/particles.h
            ${SOURCE_DIR}/Particles/setparticles.c
            )
endif ()

if (WITH_FD)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/fd_flux.c
            ${SOURCE_DIR}/fd_reconstruct.c
            ${SOURCE_DIR}/States/fd_states.c
            )
endif ()


# Math tools

set(SOURCE_FILES  ${SOURCE_FILES}
        ${SOURCE_DIR}/Math_Tools/math_interp.c
        ${SOURCE_DIR}/Math_Tools/math_lu_decomp.c
        ${SOURCE_DIR}/Math_Tools/math_misc.c
        ${SOURCE_DIR}/Math_Tools/math_ode.c
        ${SOURCE_DIR}/Math_Tools/math_qr_decomp.c
        ${SOURCE_DIR}/Math_Tools/math_quadrature.c
        ${SOURCE_DIR}/Math_Tools/math_root_finders.c
        ${SOURCE_DIR}/Math_Tools/math_table2D.c
        ${SOURCE_DIR}/Math_Tools/math_tools.h
        )



# Other flags




# Include, link, and executable

include_directories(${INCLUDE_DIRECTORIES})

link_directories(${LINK_DIRECTORIES})

link_libraries(${LINK_LIBRARIES})

add_executable(pluto ${SOURCE_FILES})


