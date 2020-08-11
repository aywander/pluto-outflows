    set(SOURCE_FILES  ${SOURCE_FILES}
            )

endif ()

endif ()


# Dust Fluid

if (WITH_DUST_FLUID)

    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Dust)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/Dust/dust.h
            ${SOURCE_DIR}/Dust/dust.c
            )
endif ()


# Forced turbulence

if (WITH_FORCED_TURB)

    set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${SOURCE_DIR}/Forced_Turb)
    set(SOURCE_FILES  ${SOURCE_FILES}
            ${SOURCE_DIR}/Forced_Turb/forced_turb.h
            ${SOURCE_DIR}/Forced_Turb/forced_turb.c
            )

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
        ${SOURCE_DIR}/Math_Tools/math_rand_mt19937-64.c
        ${SOURCE_DIR}/Math_Tools/math_random.c
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

