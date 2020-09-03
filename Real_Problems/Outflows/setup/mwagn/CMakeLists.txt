# -- Simulation dependent settings --

set(WITH_HDF5                   OFF)
set(WITH_PNG                    ON)
set(WITH_PARALLEL               ON)
set(WITH_ASYNC_IO               OFF)
set(WITH_CHOMBO                 OFF)

set(WITH_PHYSICS                RHD)
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
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY /scratch/n72/ayw100/milky_way_agn/inner/)

# New files
set(NEW_FILES
        accretion.c
        accretion.h
        clouds.c
        clouds.h
        setup/${SETUP}/definitions_usr.h
        grid_geometry.c
        grid_geometry.h
        hot_halo.c
        hot_halo.h
        idealEOS.c
        idealEOS.h
        init_tools.c
        init_tools.h
        interpolation.c
        interpolation.h
        io_tools.c
        io_tools.h
        macros_usr.h
        multicloud_init.c
        multicloud_init.h
        nozzle.c
        nozzle.h
        outflow.c
        outflow.h
        pluto_usr.h
        random.c
        random.h
        read_grav_table.c
        read_grav_table.h
        read_hot_table.c
        read_hot_table.h
        read_mu_table.c
        read_mu_table.h
        sfr.c
        sfr.h
        #supernovae.c
        #supernovae.h
        )

# Overridden files
set(OVERRIDE_FILES
        Cooling/TABULATED/radiat.c
        cmd_line_opt.c
        colortable.c
        cooling_source.c
        input_data.c
        main.c
        mappers3D.c       # Not confident about the modifications inside
        mean_mol_weight.c
        runtime_setup.c
        setup/${SETUP}/userdef_output.c
        )

