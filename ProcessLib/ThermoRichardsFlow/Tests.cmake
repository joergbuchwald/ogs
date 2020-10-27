if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE ThermoRichardsFlow/RichardsFlow2D/RichardsFlow_2d_small.prj)
endif()

AddTest(
    NAME ThermoRichardsFlow_OgataBanks
    PATH ThermoRichardsFlow/OgataBanks
    EXECUTABLE ogs
    EXECUTABLE_ARGS OgataBanks.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    HT_OgataBanks_ts_50_t_50000.000000.vtu OgataBanks_ts_50_t_50000.000000.vtu temperature  temperature_interpolated 5e-3 1e-10
    HT_OgataBanks_ts_50_t_50000.000000.vtu OgataBanks_ts_50_t_50000.000000.vtu pressure  pressure_interpolated 1e-10 1e-10
)
