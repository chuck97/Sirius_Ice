#pragma once
#include "external.h"
#include "string_trim.h"
#include "enum_classes.h"

// mesh params
class MeshParams
{
public:
    MeshParams();
    MeshParams(const std::string& json_path);
    std::string GetMeshPath() const;
    CoordsType GetCoordsType() const;
    SurfaceType GetSurfaceType() const;
    void Log() const;

private:
    std::string mesh_path;
    CoordsType coords_type;
    SurfaceType surface_type;
};

// model params
class ModelParams
{
public:
    ModelParams();
    ModelParams(const std::string& json_path);
    double GetTimeStepHours() const;
    double GetTotalTimeHours() const;
    double GetGravity() const;
    double GetWaterDensity() const;
    double GetAirDensity() const;
    double GetIceDensity() const;
    double GetWaterDragCoeff() const;
    double GetAirDragCoeff() const;
    double GetPressureCoeff() const;
    double GetPressureStar() const;
    double GetEccentricity() const;
    double GetDeltaMin() const;
    double GetCoriolisParam() const;
    double GetEarthAngularVel() const;
    double GetEarthRadius() const;
    double GetMinConcentration() const;
    void Log() const;

private:
    double time_step_hours;
    double total_time_hours;
    double g_gravity;
    double w_density;
    double a_density;
    double i_density;
    double w_drag_coeff;
    double a_drag_coeff;
    double pressure_coeff;
    double pressure_star;
    double eccentricity;
    double delta_min;
    double Coriolis_parameter;
    double earth_angular_vel;
    double earth_radius;
    double min_conc;
};

// advection params
class AdvectionParams
{
public:
    AdvectionParams();
    AdvectionParams(const std::string& json_path);
    AdvectionParams(AdvectionSolverType ast_, bool is_fct_, double cd_);
    AdvectionSolverType GetAdvectionSolverType() const;
    bool GetIsFct() const;
    double GetFctCd() const;
    void Log() const;

private:
    AdvectionSolverType advection_solver_type;
    bool is_fct;
    double fct_cd;
};

// output params
class OutputParams
{
public:
    OutputParams();
    OutputParams(const std::string& json_path);
    std::string GetOutputPvtuDirectory() const;
    std::string GetOutputNetcdfDirectory() const;
    std::string GetOutputErrorsFile() const;
    int GetEveryNthScreenshot() const;
    std::map<ModelVariableNotation, bool> GetDisplayedVariables() const;
    bool GetIsVerposeOutput() const;
    bool GetIsPvtuOutput() const;
    bool GetIsNetcdfOutput() const;
    bool GetIsErrorsOutput() const;
    void Log() const;

private:
    std::string output_pvtu_directory;
    std::string output_netcdf_directory;
    std::string output_errors_file;
    bool is_pvtu_output = false;
    bool is_netcdf_output = false;
    bool is_errors_output = false;
    int every_nth_screenshot = 1;
    std::map<ModelVariableNotation, bool> displayed_variables;
    bool is_verbose_output = true;
};

// momentum params
struct MomentumParam
{
    int Int;
    double Double;
    std::string String;
    bool Bool;    
};

class MomentumParams
{
public:
    MomentumParams();
    MomentumParams(const std::string& json_path);
    MomentumSolverType GetMomentumSolverType() const;
    bool GetIsCoriolis() const;
    bool GetIsWaterDrag() const;
    bool GetIsAirDrag() const;
    std::map<std::string, MomentumParam> GetSolverParams() const;
    MomentumBC GetMomentumBC() const;
    void Log() const;

private:
    MomentumSolverType momentum_solver_type;
    std::map<std::string, MomentumParam> solver_params;
    bool is_Coriolis = true;
    bool is_water_drag = true;
    bool is_air_drag = true;
    MomentumBC boundary_conditions = MomentumBC::no_slip;
};

// dynamics test params
class DynamicsTestParams
{
public:
    DynamicsTestParams();
    DynamicsTestParams(const std::string& json_path);
    double GetDomainSize() const;
    double GetMaxWaterSpeed() const;
    double GetMaxAirSpeed() const;
    double GetAirReductionFactor() const;
    double GetAirConvergenceAngle() const;
    double GetInitialIceConcentration() const;
    double GetInitialHeightBackground() const;
    double GetInitialHeightScaleFactor() const;
    double GetInitialHeightXfactor() const;
    double GetInitialHeightYfactor() const;
    void Log() const;

private:
    double domain_size = 512000.0;
    double v_max_water = 0.01;
    double v_max_air = 15.0;
    double reduction_factor_air = 0.02;
    double convergence_angle_air = 72.0*(M_PI/180.0);
    double initial_a_ice = 1.0;
    double background_h_ice = 0.3;
    double scale_h_ice = 0.005;
    double x_factor_h_ice = 0.00006;
    double y_factor_h_ice = 0.00003;
};

// all config params
class ConfigParams
{
public:
    ConfigParams(const std::string& json_path);
    MeshParams& GetMeshParams();
    ModelParams& GetModelParams();
    AdvectionParams& GetAdvectionParams();
    OutputParams& GetOutputParams();
    MomentumParams& GetMomentumParams();

    // specific dynamics test params
    DynamicsTestParams& GetDynamicsTestParams();
    
private:
    MeshParams mesh_params;
    ModelParams model_params;
    AdvectionParams advection_params;
    OutputParams output_params;
    MomentumParams momentum_params;
    DynamicsTestParams dynamics_test_params;
};

