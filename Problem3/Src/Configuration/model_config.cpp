#include "config.h"

using namespace std;

ModelParams::ModelParams()
{};

ModelParams::ModelParams(const string& json_path)
{
    string no_spaces_json = json_path;
    rtrim(no_spaces_json);
    if(no_spaces_json.substr(no_spaces_json.size()-5, 5) != ".json")
    {
        INMOST_ICE_ERR("input file shoud be ended by .json!");
    }
    ifstream ifs(no_spaces_json);
    nlohmann::json j_input = nlohmann::json::parse(ifs);
    BARRIER
    // Parse time step in hours
    if (j_input["time step (h)"].empty())
    {
        INMOST_ICE_ERR("time step in hours should be given!");
    }
    else
    {
        time_step_hours = j_input["time step (h)"];
    }
    BARRIER
    // Parse total time in hours
    if (j_input["total time (h)"].empty())
    {
        INMOST_ICE_ERR("total time in hours should be given!");
    }
    else
    {
        total_time_hours = j_input["total time (h)"];
    }
    BARRIER
    // Parse gravity acceleration coefficient
    if (j_input["gravity acceleration (m/s^2)"].empty())
    {
        INMOST_ICE_ERR("gravity acceleration (m/s^2) should be given!");
    }
    else
    {
        g_gravity = j_input["gravity acceleration (m/s^2)"];
    }
    BARRIER
    // Parse water density
    if (j_input["water density (kg/m^3)"].empty())
    {
        INMOST_ICE_ERR("water density (kg/m^3) should be given!");
    }
    else
    {
        w_density = j_input["water density (kg/m^3)"];
    }
    BARRIER
    // Parse air density
    if (j_input["air density (kg/m^3)"].empty())
    {
        INMOST_ICE_ERR("air density (kg/m^3) should be given!");
    }
    else
    {
        a_density = j_input["air density (kg/m^3)"];
    }
    BARRIER
    // Parse ice density
    if (j_input["ice density (kg/m^3)"].empty())
    {
        INMOST_ICE_ERR("ice density (kg/m^3) should be given!");
    }
    else
    {
        i_density = j_input["ice density (kg/m^3)"];
    }
    BARRIER
    // Parse water-ice drag coefficient
    if (j_input["water-ice drag coefficient"].empty())
    {
        INMOST_ICE_ERR("water-ice drag coefficient should be given!");
    }
    else
    {
        w_drag_coeff = j_input["water-ice drag coefficient"];
    }
    BARRIER
    // Parse air-ice drag coefficient
    if (j_input["air-ice drag coefficient"].empty())
    {
        INMOST_ICE_ERR("air-ice drag coefficient should be given!");
    }
    else
    {
        a_drag_coeff = j_input["air-ice drag coefficient"];
    }
    BARRIER
    // Parse pressure coefficient
    if (j_input["pressure coefficient"].empty())
    {
        INMOST_ICE_ERR("pressure coefficient should be given!");
    }
    else
    {
        pressure_coeff = j_input["pressure coefficient"];
    }
    BARRIER
    // Parse pressure star coefficient
    if (j_input["pressure star coefficient (Pa)"].empty())
    {
        INMOST_ICE_ERR("pressure star coefficient (Pa) should be given!");
    }
    else
    {
        pressure_star = j_input["pressure star coefficient (Pa)"];
    }
    BARRIER
    // Parse ellipse eccentricity parameter
    if (j_input["stress ellipse eccentricity"].empty())
    {
        INMOST_ICE_ERR("stress ellipse eccentricity should be given!");
    }
    else
    {
        eccentricity = j_input["stress ellipse eccentricity"];
    }
    BARRIER
    // Parse minimal delta  
    if (j_input["minimal delta (1/s)"].empty())
    {
        INMOST_ICE_ERR("minimal delta (1/s) should be given!");
    }
    else
    {
        delta_min = j_input["minimal delta (1/s)"];
    }
    BARRIER
    // Parse Coriolis parameter  
    if (j_input["Coriolis parameter (1/s)"].empty())
    {
        INMOST_ICE_ERR("Coriolis parameter (1/s) should be given!");
    }
    else
    {
        Coriolis_parameter = j_input["Coriolis parameter (1/s)"];
    }
    BARRIER

    // Parse min concentration   
    if (j_input["minimal concentration"].empty())
    {
        INMOST_ICE_ERR("minimal concentration should be given!");
    }
    else
    {
        min_conc = j_input["minimal concentration"];
    }
    BARRIER

    // Parse earth angular velocity  
    //if (j_input["earth angular velocity (1/s)"].empty())
    //{
    //    INMOST_ICE_ERR("earth angular velocity (1/s) should be given!");
    //}
    //else
    //{
    //    earth_angular_vel = j_input["earth angular velocity (1/s)"];
    //}
    //BARRIER
    // Parse earth radius  
    //if (j_input["earth radius (m)"].empty())
    //{
    //    INMOST_ICE_ERR("earth radius (m) should be given!");
    //}
    //else
    //{
    //    earth_radius = j_input["earth radius (m)"];
    //}
    //BARRIER
};

double ModelParams::GetTimeStepHours() const
{
    return time_step_hours;
};

double ModelParams::GetTotalTimeHours() const
{
    return total_time_hours;
};

double ModelParams::GetGravity() const
{
    return g_gravity;
};

double ModelParams::GetWaterDensity() const
{
    return w_density;
};

double ModelParams::GetAirDensity() const
{
    return a_density;
};

double ModelParams::GetIceDensity() const
{
    return i_density;
};

double ModelParams::GetWaterDragCoeff() const
{
    return w_drag_coeff;
};

double ModelParams::GetAirDragCoeff() const
{
    return a_drag_coeff;
};

double ModelParams::GetPressureCoeff() const
{
    return pressure_coeff;
};

double ModelParams::GetPressureStar() const
{
    return pressure_star;
};

double ModelParams::GetEccentricity() const
{
    return eccentricity;
};

double ModelParams::GetDeltaMin() const
{
    return delta_min;
};

double ModelParams::GetCoriolisParam() const
{
    return Coriolis_parameter;
};

double ModelParams::GetEarthAngularVel() const
{
    return earth_angular_vel;
};

double ModelParams::GetEarthRadius() const
{
    return earth_radius;
};

double ModelParams::GetMinConcentration() const
{
    return min_conc;
};

void ModelParams::Log() const
{
    cout << "Time step (h): " << time_step_hours << endl;
    cout << "Total time (h): " << total_time_hours << endl;
    cout << "Gravity acceleration (m/s^2): " << g_gravity << endl;
    cout << "Water density (kg/m^3): " << w_density << endl;
    cout << "Air density (kg/m^3): " << a_density << endl;
    cout << "Ice density (kg/m^3): " << i_density << endl;
    cout << "Water-Ice drag coefficient: " << w_drag_coeff << endl;
    cout << "Air-Ice drag coefficient: " << a_drag_coeff << endl;
    cout << "Pressure coefficient: " << pressure_coeff << endl;
    cout << "Pressure star coefficient (Pa): " << pressure_star << endl;
    cout << "Stress ellipse eccentricity: " << eccentricity << endl;
    cout << "Minimal delta value (1/s): " << delta_min << endl;
    cout << "Coriolis parameter (1/s): " << Coriolis_parameter << endl;
    cout << "Minimal concentration: " << min_conc << endl;
    //cout << "Earth angular velocity (1/s): " << earth_angular_vel << endl;
    //cout << "Earth radius (m): " << earth_radius << endl;
};