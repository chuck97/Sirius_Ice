#include "config.h"

using namespace std;

ConfigParams::ConfigParams(const string& json_path)
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

    // Parse mesh params .json path
    if (j_input["mesh params"].empty())
    {
        INMOST_ICE_ERR("Mesh params .json path should be given!");
    }
    else
    {
        mesh_params = MeshParams(j_input["mesh params"]);
    }
    BARRIER

    // Parse model params .json path
    if (j_input["model params"].empty())
    {
        INMOST_ICE_ERR("Model params .json path should be given!");
    }
    else
    {
        model_params = ModelParams(j_input["model params"]);
    }
    BARRIER

    // Parse advection params .json path
    if (j_input["advection params"].empty())
    {
        cout << "No advection params! Dafault advection solver will be used (TG2+FCT, cd = 0.5)" << endl;
        advection_params = {AdvectionSolverType::TG2, true, 0.5};
    }
    else
    {
        advection_params = AdvectionParams(j_input["advection params"]);
    }
    BARRIER

    // Parse output params .json path
    if (j_input["output params"].empty())
    {
        INMOST_ICE_ERR("Output params .json path should be given!");
    }
    else
    {
        output_params = OutputParams(j_input["output params"]);
    }
    BARRIER

    // Parse momentum params .json path
    if (j_input["momentum params"].empty())
    {
        INMOST_ICE_ERR("Momentum params .json path should be given!");
    }
    else
    {
        momentum_params = MomentumParams(j_input["momentum params"]);
    }
    BARRIER

    // Parse dynamics test params .json path
    if (!j_input["dynamics test params"].empty())
    {
        dynamics_test_params = DynamicsTestParams(j_input["dynamics test params"]);
    }
    BARRIER
};

MeshParams& ConfigParams::GetMeshParams()
{
    return mesh_params;
};

ModelParams& ConfigParams::GetModelParams()
{
    return model_params;
};

AdvectionParams& ConfigParams::GetAdvectionParams()
{
    return advection_params;
};

OutputParams& ConfigParams::GetOutputParams()
{
    return output_params;
};

MomentumParams& ConfigParams::GetMomentumParams()
{
    return momentum_params;
};

DynamicsTestParams& ConfigParams::GetDynamicsTestParams()
{
    return dynamics_test_params;
};