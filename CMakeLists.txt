cmake_minimum_required(VERSION 3.3)

project(PLUTO)

#_______________________________________________________________________________________________________________________
# Set these

set(PLUTO_SOURCE_DIR            /Users/ayw/dev/PLUTO)
#set(PLUTO_SOURCE_DIR            C:/Users/Utilisateur/Desktop/PLUTO)
set(PROBLEM_TYPE                Real)
set(PROBLEM                     Outflows)

set(WITH_PARALLEL               OFF)
set(WITH_PNG                    OFF)
set(WITH_HDF5                   ON)
set(WITH_ASYNC_IO               OFF)
set(WITH_CHOMBO                 OFF)

set(WITH_PHYSICS                HD)
set(WITH_COOLING                TABULATED)
set(WITH_INTERPOLATION          PARABOLIC)
set(WITH_TIME_STEPPING          RK3)

set(WITH_PARTICLES              OFF)
set(WITH_FARGO                  OFF)
set(WITH_SHEARING_BOX           OFF)
set(WITH_RESISTIVITY            OFF)
set(WITH_FD                     OFF)

# -- physics dependent declarations --

set(WITH_EOS                    IDEAL)
set(WITH_THERMAL_CONDUCTION     OFF)
set(WITH_VISCOSITY              OFF)
set(WITH_ROTATING_FRAME         OFF)
set(WITH_CT                     OFF)
set(WITH_GLM                    OFF)


# -- Library paths --

set(HDF5_PATH ${PLUTO_SOURCE_DIR}/Lib)

#_______________________________________________________________________________________________________________________
# Main bit

# The dependences of which files to build depending on which module can be found in the individual makefiles in the
# module folders, and in AppendAdditionalFiles() in Tools/Python/define_problem.py .



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

set(LINK_LIBRARIES)


# Setup-dependent / Problem-type dependent files

set(SOURCE_FILES
        ${PROBLEM_DIR}/definitions.h
        ${PROBLEM_DIR}/init.c
        )

if (PROBLEM_TYPE STREQUAL Real)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${PROBLEM_DIR}/definitions_usr.h
            ${PROBLEM_DIR}/macros_usr.h
            ${PROBLEM_DIR}/pluto_usr.h
            )
endif()


# Problem source files in problem directory
# NOTE: Module dependent overwritten data can be found
#       in the Module/Physics dependent part below.

if (PROBLEM STREQUAL Outflows)

#     TODO: Automate necessary file link to run simulations
    # Executable output location
#    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PLUTO_SOURCE_DIR}/Lab/Accretion)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY /Users/ayw/lab/AGN_Feedback/Accretion/wind)
#    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY /Users/ayw/lab/AGN_Feedback/IC5063)
#    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY /Users/ayw/lab/AGN_Feedback/Outflow)

    # New files
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${PROBLEM_DIR}/accretion.c
            ${PROBLEM_DIR}/accretion.h
            ${PROBLEM_DIR}/clouds.c
            ${PROBLEM_DIR}/clouds.h
            ${PROBLEM_DIR}/grid_geometry.c
            ${PROBLEM_DIR}/grid_geometry.h
            ${PROBLEM_DIR}/hot_halo.c
            ${PROBLEM_DIR}/hot_halo.h
            ${PROBLEM_DIR}/idealEOS.c
            ${PROBLEM_DIR}/idealEOS.h
            ${PROBLEM_DIR}/init_tools.c
            ${PROBLEM_DIR}/init_tools.h
            ${PROBLEM_DIR}/interpolation.c
            ${PROBLEM_DIR}/interpolation.h
            ${PROBLEM_DIR}/mean_mol_weight.c
            ${PROBLEM_DIR}/multicloud_init.c
            ${PROBLEM_DIR}/multicloud_init.h
            ${PROBLEM_DIR}/outflow.c
            ${PROBLEM_DIR}/outflow.h
            ${PROBLEM_DIR}/random.c
            ${PROBLEM_DIR}/random.h
            ${PROBLEM_DIR}/read_grav_table.c
            ${PROBLEM_DIR}/read_grav_table.h
            ${PROBLEM_DIR}/read_hot_table.c
            ${PROBLEM_DIR}/read_hot_table.h
            ${PROBLEM_DIR}/read_mu_table.c
            ${PROBLEM_DIR}/read_mu_table.h
#            ${PROBLEM_DIR}/supernovae.c
#            ${PROBLEM_DIR}/supernovae.h
            )

    # Overridden files
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${PROBLEM_DIR}/cmd_line_opt.c
            ${PROBLEM_DIR}/input_data.c
            ${PROBLEM_DIR}/jet_domain.c
            ${PROBLEM_DIR}/main.c
            ${PROBLEM_DIR}/mappers3D.c       # Not confident about the modifications inside
            ${PROBLEM_DIR}/set_grid.c
            ${PROBLEM_DIR}/startup.c
            ${PROBLEM_DIR}/userdef_output.c
            )


elseif (PROBLEM STREQUAL SupernovaFeedback)

    # Executable output location
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY /Users/ayw/lab/SupernovaFeedback/Supernova)
#    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY /Users/ayw/lab/SupernovaFeedback/Starburst)

    # New files
    set(SOURCE_FILES  ${SOURCE_FILES}
            )

    # Overridden files
    set(SOURCE_FILES  ${SOURCE_FILES}
            )


elseif (PROBLEM STREQUAL RadiativeShock)

    # Executable output location
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY /Users/ayw/lab/Radiative_Shocks/PLUTO_simulations/runs/ovi)

    # New files
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${PROBLEM_DIR}/init_tools.c
            ${PROBLEM_DIR}/init_tools.h
            ${PROBLEM_DIR}/interpolation.c
            ${PROBLEM_DIR}/interpolation.h
            ${PROBLEM_DIR}/read_mu_table.c
            ${PROBLEM_DIR}/read_mu_table.h
            ${PROBLEM_DIR}/shock.c
            ${PROBLEM_DIR}/shock.h
            )

    # Overridden files
    # TODO: This was wrong (PROBLEM_DIR was SOURCE_DIR)! Might not compile now
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${PROBLEM_DIR}/cmd_line_opt.c
            ${PROBLEM_DIR}/input_data.c
            ${PROBLEM_DIR}/jet_domain.c
            ${PROBLEM_DIR}/main.c
            #${PROBLEM_DIR}/mappers3D.c
            ${PROBLEM_DIR}/set_grid.c
            ${PROBLEM_DIR}/startup.c
            ${PROBLEM_DIR}/userdef_output.c
            )

elseif (PROBLEM STREQUAL HD/Jet)

    # Executable output location
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PLUTO_SOURCE_DIR}/Test_Problems/${PROBLEM})

    # TODO: Change the way problem-dependent overrides are handled. Not trivial as files in modules may be overridden.
    # Non-overridden Src files
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/cmd_line_opt.c
            ${SOURCE_DIR}/input_data.c
            ${SOURCE_DIR}/jet_domain.c
            ${SOURCE_DIR}/main.c
            ${SOURCE_DIR}/set_grid.c
            ${SOURCE_DIR}/startup.c
            )

    # Overridden files
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${PROBLEM_DIR}/userdef_output.c
            )

elseif (PROBLEM STREQUAL HD/Stellar_Wind)

    # Executable output location
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PLUTO_SOURCE_DIR}/Test_Problems/${PROBLEM})

    # TODO: Change the way problem-dependent overrides are handled. Not trivial as files in modules may be overridden.
    # Non-overridden Src files
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/cmd_line_opt.c
            ${SOURCE_DIR}/input_data.c
            ${SOURCE_DIR}/jet_domain.c
            ${SOURCE_DIR}/main.c
            ${SOURCE_DIR}/set_grid.c
            ${SOURCE_DIR}/startup.c
            ${SOURCE_DIR}/userdef_output.c
            )

endif()



# Files under SRC. Used for all problems

# TODO: Change the way problem-dependent overrides are handled. Not trivial as files in modules may be overridden.
set(SOURCE_FILES  ${SOURCE_FILES}
        ${SOURCE_DIR}/adv_flux.c
        ${SOURCE_DIR}/arrays.c
        ${SOURCE_DIR}/bin_io.c
        ${SOURCE_DIR}/boundary.c
        #${SOURCE_DIR}/cartcoord.c                        # Never used
        ${SOURCE_DIR}/check_states.c
        #${SOURCE_DIR}/cmd_line_opt.c                     # Override
        ${SOURCE_DIR}/colortable.c
        #${SOURCE_DIR}/cooling_source.c                   # Override
        ${SOURCE_DIR}/doxy-comment.c
        ${SOURCE_DIR}/entropy_switch.c
        #${SOURCE_DIR}/failsafe.c                         # Never used
        ${SOURCE_DIR}/flag_shock.c
        ${SOURCE_DIR}/get_nghost.c
        ${SOURCE_DIR}/globals.h                           # Modified
        ${SOURCE_DIR}/initialize.c
        #${SOURCE_DIR}/input_data.c                       # Override
        ${SOURCE_DIR}/int_bound_reset.c
        #${SOURCE_DIR}/jet_domain.c                       # Override
        ${SOURCE_DIR}/macros.h
        #${SOURCE_DIR}/main.c                             # Override
#        ${SOURCE_DIR}/mappers3D.c                         # Override
#        ${SOURCE_DIR}/mean_mol_weight.c                  # Override
        ${SOURCE_DIR}/parabolic_flux.c
        ${SOURCE_DIR}/parabolic_rhs.c
        ${SOURCE_DIR}/parse_file.c
        ${SOURCE_DIR}/pluto.h                             # Modified
        ${SOURCE_DIR}/prototypes.h
        ${SOURCE_DIR}/rbox.c
        ${SOURCE_DIR}/restart.c
#        ${SOURCE_DIR}/rkc.c                              # Buggy (not used?)
        ${SOURCE_DIR}/runtime_setup.c
        ${SOURCE_DIR}/set_geometry.c
        #${SOURCE_DIR}/set_grid.c                         # Overwritten
        ${SOURCE_DIR}/set_image.c
        ${SOURCE_DIR}/set_indexes.c
        ${SOURCE_DIR}/set_output.c
        ${SOURCE_DIR}/show_config.c
        ${SOURCE_DIR}/split_source.c
        #${SOURCE_DIR}/startup.c                          # Override
        ${SOURCE_DIR}/structs.h                           # Modified
        ${SOURCE_DIR}/sts.c
        ${SOURCE_DIR}/tools.c
        #${SOURCE_DIR}/userdef_output.c                   # Override
        #${SOURCE_DIR}/userdef_output.dt.c                # Never used?
        #${SOURCE_DIR}/userdef_output.shock_tracking.c    # Never used?
        ${SOURCE_DIR}/var_names.c
        ${SOURCE_DIR}/vec_pot_diff.c
        ${SOURCE_DIR}/vec_pot_update.c
        ${SOURCE_DIR}/write_data.c
        ${SOURCE_DIR}/write_img.c
        ${SOURCE_DIR}/write_tab.c
        ${SOURCE_DIR}/write_vtk.c

        ${SOURCE_DIR}/Math_Tools/math_interp.c
        ${SOURCE_DIR}/Math_Tools/math_lu_decomp.c
        ${SOURCE_DIR}/Math_Tools/math_misc.c
        ${SOURCE_DIR}/Math_Tools/math_ode.c
        ${SOURCE_DIR}/Math_Tools/math_qr_decomp.c
        ${SOURCE_DIR}/Math_Tools/math_quadrature.c
        ${SOURCE_DIR}/Math_Tools/math_root_finders.c
        ${SOURCE_DIR}/Math_Tools/math_table2D.c
        ${SOURCE_DIR}/Math_Tools/math_tools.h

        ${SOURCE_DIR}/Time_Stepping/update_stage.c

        ${SOURCE_DIR}/States/flatten.c
        )





# Module/Physics dependent part

# Make flags options

if (WITH_PARALLEL)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPARALLEL")

    # ----

    set(CMAKE_C_COMPILER /usr/local/bin/mpicc)
    set(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)

    # OR (https://cmake.org/pipermail/cmake/2015-April/060252.html)

#    find_package(MPI REQUIRED)
#    set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
#    set(CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_LINK_FLAGS})
#    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
#    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
#    set(LINK_LIBRARIES ${LINK_LIBRARIES} ${MPI_LIBRARIES})

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
        set(HDF5_ROOT ${HDF5_PATH}/hdf5-parallel)
    else()
        set(HDF5_ROOT ${HDF5_PATH}/hdf5-serial)
    endif ()

    find_package(HDF5 REQUIRED)


    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${HDF5_INCLUDE_DIRS})
    set(LINK_LIBRARIES ${LINK_LIBRARIES} ${HDF5_LIBRARIES})

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
#            ${SOURCE_DIR}/HD/riemann_full.c       # Buggy (not used?)
            ${SOURCE_DIR}/HD/roe.c
#            ${SOURCE_DIR}/HD/rusanov-dw.c         # Don't need
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
#            ${SOURCE_DIR}/MHD/roe.bckfld.c          # Don't need
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
#            ${SOURCE_DIR}/RMHD/quartic.devel1.c       # Not used
#            ${SOURCE_DIR}/RMHD/quartic.new.c          # Not used
#            ${SOURCE_DIR}/RMHD/quartic.old.c          # Not used
            ${SOURCE_DIR}/RMHD/rhs.c
#            ${SOURCE_DIR}/RMHD/rusanov-dw.c           # Not used
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



# Other flags




# Include, link, and executable

include_directories(${INCLUDE_DIRECTORIES})

link_directories(${LINK_DIRECTORIES})

link_libraries(${LINK_LIBRARIES})

add_executable(pluto ${SOURCE_FILES})


