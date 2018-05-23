# -- Simulation dependent settings --

set(WITH_PARALLEL               OFF)
set(WITH_PNG                    OFF)
set(WITH_HDF5                   ON)
set(WITH_ASYNC_IO               OFF)
set(WITH_CHOMBO                 OFF)

set(WITH_PHYSICS                HD)
set(WITH_COOLING                NONE)
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


# Executable output location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY /Users/ayw/lab/DynamicalFriction/PLUTO_simulations/runs)

# New files
set(NEW_FILES
        pluto_usr.h
        definitions_usr.h
        friction.h
        friction.c
        macros_usr.h
        init_tools.c
        init_tools.h
        interpolation.c
        interpolation.h
        io_tools.c
        io_tools.h
        )

# Overridden files
set(OVERRIDE_FILES
        cmd_line_opt.c
        input_data.c
        jet_domain.c
        main.c
        #mappers3D.c
        startup.c
        set_grid.c
        userdef_output.c
        )

