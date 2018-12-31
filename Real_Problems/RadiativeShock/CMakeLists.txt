# -- Simulation dependent settings --

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
set(WITH_FORCED_TURB            OFF)
set(WITH_FARGO                  OFF)
set(WITH_SHEARING_BOX           OFF)
set(WITH_RESISTIVITY            OFF)
set(WITH_FD                     OFF)


# -- physics dependent declarations --

set(WITH_EOS                    IDEAL)
set(WITH_THERMAL_CONDUCTION     OFF)
set(WITH_VISCOSITY              OFF)
set(WITH_AMBIPOLAR_DIFFUSION    OFF)
set(WITH_ROTATING_FRAME         OFF)
set(WITH_CT                     OFF)
set(WITH_GLM                    OFF)
set(WITH_HALL_MHD               OFF)
set(WITH_DUST_FLUID             OFF)


# Executable output location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY /Users/ayw/lab/Radiative_Shocks/PLUTO_simulations/runs/ovi)

# New files
set(NEW_FILES
        init_tools.c
        init_tools.h
        interpolation.c
        interpolation.h
        read_mu_table.c
        read_mu_table.h
        shock.c
        shock.h
        )

# Overridden files
set(OVERRIDE_FILES
        cmd_line_opt.c
        input_data.c
        mappers3D.c
        mean_mol_weight.c
        runtime_setup.c
        userdef_output.c
        )

