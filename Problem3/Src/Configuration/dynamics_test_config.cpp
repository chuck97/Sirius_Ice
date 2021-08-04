#include "config.h"

using namespace std;

DynamicsTestParams::DynamicsTestParams()
{};

DynamicsTestParams::DynamicsTestParams(const std::string& json_path)
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

    // Parse domain size
    if (!j_input["domain size (m)"].empty())
    {
        domain_size = j_input["domain size (m)"];
    }
    BARRIER

    // Parse max water speed
    if (!j_input["v max water (m/s)"].empty())
    {
        v_max_water = j_input["v max water (m/s)"];
    }
    BARRIER   

    // Parse max air speed
    if (!j_input["v max air (m/s)"].empty())
    {
        v_max_air = j_input["v max air (m/s)"];
    }
    BARRIER

    // Parse air reduction factor 
    if (!j_input["reduction factor air"].empty())
    {
        reduction_factor_air = j_input["reduction factor air"];
    }
    BARRIER

    // Parse air convergence angle 
    if (!j_input["convergence angle air (deg)"].empty())
    {
        convergence_angle_air = j_input["convergence angle air (deg)"]*(M_PI/180.0);
    }
    BARRIER

    // Parse initial ice concentration
    if (!j_input["initial ice concentration"].empty())
    {
        initial_a_ice = j_input["initial ice concentration"];
    }
    BARRIER

    // Parse initial ice concentration
    if (!j_input["initial ice height background (m)"].empty())
    {
        background_h_ice = j_input["initial ice height background (m)"];
    }
    BARRIER

    // Parse initial ice hight scale factor
    if (!j_input["initial ice height scale factor (m)"].empty())
    {
        scale_h_ice = j_input["initial ice height scale factor (m)"];
    }
    BARRIER

    // Parse x ice hight scale factor
    if (!j_input["initial ice height x factor"].empty())
    {
        x_factor_h_ice = j_input["initial ice height x factor"];
    }
    BARRIER

    // Parse y ice hight scale factor
    if (!j_input["initial ice height y factor"].empty())
    {
        y_factor_h_ice = j_input["initial ice height y factor"];
    }
    BARRIER
};

double DynamicsTestParams::GetDomainSize() const
{
    return domain_size;
};

double DynamicsTestParams::GetMaxWaterSpeed() const
{
    return v_max_water;
};

double DynamicsTestParams::GetMaxAirSpeed() const
{
    return v_max_air;
};

double DynamicsTestParams::GetAirReductionFactor() const
{
    return reduction_factor_air;
};

double DynamicsTestParams::GetAirConvergenceAngle() const
{
    return convergence_angle_air;
};

double DynamicsTestParams::GetInitialIceConcentration() const
{
    return initial_a_ice;
};

double DynamicsTestParams::GetInitialHeightBackground() const
{
    return background_h_ice;
};

double DynamicsTestParams::GetInitialHeightScaleFactor() const
{
    return scale_h_ice;
};

double DynamicsTestParams::GetInitialHeightXfactor() const
{
    return x_factor_h_ice;
};

double DynamicsTestParams::GetInitialHeightYfactor() const
{
    return y_factor_h_ice;
};

void DynamicsTestParams::Log() const
{
    cout << "domain size (m): " << domain_size << endl;
    cout << "v max water (m/s): " << v_max_water << endl;
    cout << "v max air (m/s): " << v_max_air << endl;
    cout << "reduction factor air: " << reduction_factor_air << endl;
    cout << "air convergence angle (deg): " << convergence_angle_air*(180.0/M_PI) << endl;
    cout << "initial ice concentration: " << initial_a_ice << endl;
    cout << "initial ice height background (m): " << background_h_ice << endl;
    cout << "initial ice height scale factor (m): " << scale_h_ice << endl;
    cout << "initial ice height x factor: " << x_factor_h_ice << endl;
    cout << "initial ice height y factor: " << y_factor_h_ice << endl;
};