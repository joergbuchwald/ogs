AddTest(
    NAME ThermoRichardsFlow_PressureDiffusionTemperatureDiffusion
    PATH ThermoRichardsFlow/HT/SimpleSynthetics
    EXECUTABLE ogs
    EXECUTABLE_ARGS PressureDiffusionTemperatureDiffusion.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    PressureDiffusionTemperatureDiffusion_expected.vtu PressureDiffusionTemperatureDiffusion_ts_1_t_1.000000.vtu T T 1e-5 1e-8
    PressureDiffusionTemperatureDiffusion_expected.vtu PressureDiffusionTemperatureDiffusion_ts_1_t_1.000000.vtu p p 1e-5 1e-8
    PressureDiffusionTemperatureDiffusion_expected.vtu PressureDiffusionTemperatureDiffusion_ts_1_t_1.000000.vtu darcy_velocity darcy_velocity 1e-5 1e-8
)
AddTest(
    NAME ThermoRichardsFlow_HeatTransportInStationaryFlow
    PATH ThermoRichardsFlow/HT/HeatTransportInStationaryFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS HeatTransportInStationaryFlow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu temperature  temperature 5e-3 1e-8
    HT_HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu HeatTransportInStationaryFlow_ts_50_t_50000.000000.vtu pressure  pressure 5e-3 1e-8
)
AddTest(
    NAME ThermoRichardsFlow_RichardsFlow2DSmall
    PATH ThermoRichardsFlow/RichardsFlow2D
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_small.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    expected_Richards_2D_small_pcs_ts_1100_t_1600.000000.vtu Richards_2D_small_pcs_ts_1100_t_1600.000000.vtu pressure pressure 5e-3 1e-8
    expected_Richards_2D_small_pcs_ts_1100_t_1600.000000.vtu Richards_2D_small_pcs_ts_1100_t_1600.000000.vtu saturation saturation 5e-3 1e-8
)
AddTest(
    NAME ThermoRichardsFlow_RichardsFlow2DSmall_ogs5
    PATH ThermoRichardsFlow/RichardsFlow2D
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_compare_ogs5.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 17
    DIFF_DATA
    h_us_quad_1000.vtu richards_ogs5_pcs_ts_100_t_100.000000.vtu PRESSURE1 pressure 1e-1 1e-1
)
AddTest(
    NAME ThermoRichardsFlow_comp_TRMuni_saturated-TRuni_saturated
    PATH ThermoRichardsFlow/SimplifiedMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS TRuni_saturated.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    expected_TRMuni_sat_ts_10_t_1.000000.vtu TRuni_sat_ts_10_t_1.000000.vtu temperature temperature 5e-5 1e-10
    expected_TRMuni_sat_ts_10_t_1.000000.vtu TRuni_sat_ts_10_t_1.000000.vtu pressure pressure 5e-5 1e-10
    expected_TRMuni_sat_ts_10_t_1.000000.vtu TRuni_sat_ts_10_t_1.000000.vtu saturation saturation 5e-5 1e-10
)
AddTest(
    NAME ThermoRichardsFlow_comp_TRMuni_unsaturated-TRuni_unsaturated
    PATH ThermoRichardsFlow/SimplifiedMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS TRuni_unsaturated.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    expected_TRMuni_unsat_ts_10_t_1.000000.vtu TRuni_unsat_ts_10_t_1.000000.vtu temperature temperature 5e-5 1e-10
    expected_TRMuni_unsat_ts_10_t_1.000000.vtu TRuni_unsat_ts_10_t_1.000000.vtu pressure pressure 5e-3 1e-6
    expected_TRMuni_unsat_ts_10_t_1.000000.vtu TRuni_unsat_ts_10_t_1.000000.vtu saturation saturation 5e-5 1e-10
)
AddTest(
    NAME ThermoRichardsFlow_comp_TRMiso_saturated-TRiso_saturated
    PATH ThermoRichardsFlow/SimplifiedMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS TRhyd_saturated.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    expected_TRMiso_sat_ts_10_t_1.000000.vtu TRiso_sat_ts_10_t_1.000000.vtu temperature temperature 5e-5 1e-10
    expected_TRMiso_sat_ts_10_t_1.000000.vtu TRiso_sat_ts_10_t_1.000000.vtu pressure pressure 5e-5 1e-10
    expected_TRMiso_sat_ts_10_t_1.000000.vtu TRiso_sat_ts_10_t_1.000000.vtu saturation saturation 5e-5 1e-10
)
AddTest(
    NAME ThermoRichardsFlow_comp_TRMiso_unsaturated-TRiso_unsaturated
    PATH ThermoRichardsFlow/SimplifiedMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS TRhyd_unsaturated.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    expected_TRMiso_unsat_ts_10_t_1.000000.vtu TRiso_unsat_ts_10_t_1.000000.vtu temperature temperature 5e-5 1e-10
    expected_TRMiso_unsat_ts_10_t_1.000000.vtu TRiso_unsat_ts_10_t_1.000000.vtu pressure pressure 5e-3 1e-6
    expected_TRMiso_unsat_ts_10_t_1.000000.vtu TRiso_unsat_ts_10_t_1.000000.vtu saturation saturation 5e-5 1e-10
)
AddTest(
    NAME ThermoRichardsFlow_comp_TRMuni_bishopstest-TRuni_bishopstest
    PATH ThermoRichardsFlow/SimplifiedMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS TRuni_unsaturated.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    expected_TRMuni_unsat_bishopstest_ts_10_t_1.000000.vtu TRuni_unsat_ts_10_t_1.000000.vtu temperature temperature 5e-5 1e-10
    expected_TRMuni_unsat_bishopstest_ts_10_t_1.000000.vtu TRuni_unsat_ts_10_t_1.000000.vtu pressure pressure 5e-3 1e-6
    expected_TRMuni_unsat_bishopstest_ts_10_t_1.000000.vtu TRuni_unsat_ts_10_t_1.000000.vtu saturation saturation 5e-5 1e-10
)
AddTest(
    NAME ThermoRichardsFlow_comp_TRMiso_bishopstest-TRiso_bishopstest
    PATH ThermoRichardsFlow/SimplifiedMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS TRhyd_unsaturated.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    RUNTIME 1
    DIFF_DATA
    expected_TRMiso_unsat_bishopstest_ts_10_t_1.000000.vtu TRiso_unsat_ts_10_t_1.000000.vtu temperature temperature 5e-5 1e-10
    expected_TRMiso_unsat_bishopstest_ts_10_t_1.000000.vtu TRiso_unsat_ts_10_t_1.000000.vtu pressure pressure 5e-3 1e-6
    expected_TRMiso_unsat_bishopstest_ts_10_t_1.000000.vtu TRiso_unsat_ts_10_t_1.000000.vtu saturation saturation 5e-5 1e-10
)
