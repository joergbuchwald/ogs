# FEBEX1D
AddTest(
    NAME NonIsoTHM_fabex1D
    PATH NonIsothermalRichardsFlowMechanics/FEBEX1DModel
    EXECUTABLE ogs
    EXECUTABLE_ARGS febex1Dmodel.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    febex1D_ts_t_3153600.000000.vtu  febex1D_ts_t_3153600.000000.vtu  temperature temperature 1e-14 1e-12
    febex1D_ts_t_3153600.000000.vtu  febex1D_ts_t_3153600.000000.vtu  saturation saturation 1e-14 1e-12
    febex1D_ts_t_3153600.000000.vtu  febex1D_ts_t_3153600.000000.vtu  displacement displacement 1e-14 1e-12
    febex1D_ts_t_3153600.000000.vtu  febex1D_ts_t_3153600.000000.vtu   sigma sigma 1e-14 1e-12
    febex1D_ts_t_3153600.000000.vtu  febex1D_ts_t_3153600.000000.vtu   epsilon epsilon 1e-14 1e-12
    febex1D_ts_t_31536000000000.000000.vtu febex1D_ts_t_31536000000000.000000.vtu temperature temperature 1e-14 1e-12
    febex1D_ts_t_31536000000000.000000.vtu febex1D_ts_t_31536000000000.000000.vtu  saturation saturation 1e-14 1e-12
    febex1D_ts_t_31536000000000.000000.vtu febex1D_ts_t_31536000000000.000000.vtu  displacement displacement 1e-14 1e-12
    febex1D_ts_t_31536000000000.000000.vtu febex1D_ts_t_31536000000000.000000.vtu   sigma sigma 1e-14 1e-12
    febex1D_ts_t_31536000000000.000000.vtu febex1D_ts_t_31536000000000.000000.vtu   epsilon epsilon 1e-14 1e-12
)
