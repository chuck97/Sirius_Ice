#include "config.h"

using namespace std;

MomentumParams::MomentumParams()
{};

MomentumParams::MomentumParams(const string& json_path)
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

    // Parse momentum solver 

    // mEVP solver
    if (!j_input["mEVP"].empty())
    {
        momentum_solver_type = MomentumSolverType::mEVP;
        double alpha, beta, rel_resid;
        int nit;
        bool is_iterations = false;
        bool is_residual = false;

        if (!j_input["mEVP"]["alpha"].empty())
        {
            alpha = j_input["mEVP"]["alpha"];
            MomentumParam param;
            param.Double = alpha;
            solver_params["alpha"] = param;
        }
        else
        {
            INMOST_ICE_ERR("alpha value for mEVP should be given!");
        }

        if (!j_input["mEVP"]["beta"].empty())
        {
            beta = j_input["mEVP"]["beta"];
            MomentumParam param;
            param.Double = beta;
            solver_params["beta"] = param;
        }
        else
        {
            INMOST_ICE_ERR("beta value for mEVP should be given!");
        }

        if (!j_input["mEVP"]["number of iterations"].empty())
        {
            is_iterations = true;
            nit = j_input["mEVP"]["number of iterations"];
            MomentumParam param1;
            param1.Bool = is_iterations;
            MomentumParam param2;
            param2.Int = nit;
            solver_params["number of iterations"] = param2;
            solver_params["is iteration number limit"] = param1;
        }
        else
        {
            is_iterations = false;
            MomentumParam param;
            param.Bool = is_iterations;
            solver_params["is iteration number limit"] = param;
        }

        if (!j_input["mEVP"]["relative residual"].empty())
        {
            is_residual = true;
            rel_resid = j_input["mEVP"]["relative residual"];
            MomentumParam param1;
            param1.Bool = is_residual;
            MomentumParam param2;
            param2.Double = rel_resid;
            solver_params["relative residual"] = param2;
            solver_params["is relative residual limit"] = param1;
        }
        else
        {
            is_residual = false;
            MomentumParam param;
            param.Bool = is_residual;
            solver_params["is relative residual limit"] = param;
        }

        if (!is_iterations and !is_residual)
        {
            INMOST_ICE_ERR("relative resiaual or maximal iteration number should be given for mEVP!");
        }
    }
    else // default is mEVP-500
    {
        momentum_solver_type = MomentumSolverType::mEVP;
        MomentumParam param1, param2, param3, param4, param5;
        param1.Double = 500.0;
        param2.Double = 500.0;
        param3.Int = 500;
        param4.Bool = false;
        param5.Bool = true;
        solver_params["alpha"] = param1;
        solver_params["beta"] = param2;
        solver_params["number of iterations"] = param3;
        solver_params["is relative residual limit"] = param4;
        solver_params["is iteration number limit"] = param5;
    }
    BARRIER

    // Parse is Coriolis
    if (!j_input["is Coriolis"].empty())
    {
        is_Coriolis = j_input["is Coriolis"];
    }
    BARRIER

    // Parse is water drag
    if (!j_input["is water drag"].empty())
    {
        is_water_drag = j_input["is water drag"];
    }
    BARRIER

    // Parse is air drag
    if (!j_input["is air drag"].empty())
    {
        is_air_drag = j_input["is air drag"];
    }
    BARRIER

    // Parse momentum BC type
    if (!j_input["boundary conditions"].empty())
    {
        string b_c = j_input["boundary conditions"];
        if (MomentumNameToBC.count(b_c) == 0)
        {
            INMOST_ICE_ERR(b_c + " no such BC");
        }
        else
        {
            boundary_conditions = MomentumNameToBC[b_c];
        }
    }
    BARRIER
};

MomentumSolverType MomentumParams::GetMomentumSolverType() const
{
    return momentum_solver_type;
};

map<string, MomentumParam> MomentumParams::GetSolverParams() const
{
    return solver_params;
};

bool MomentumParams::GetIsCoriolis() const
{
    return is_Coriolis;
};

bool MomentumParams::GetIsWaterDrag() const
{
    return is_water_drag;
};

bool MomentumParams::GetIsAirDrag() const
{
    return is_air_drag;
};

MomentumBC MomentumParams::GetMomentumBC() const
{
    return boundary_conditions;
};

void MomentumParams::Log() const
{
    cout << "Momentum solver: " << MomentumSolverTypeToName[momentum_solver_type] << endl; 
    cout << "Is Coriolis: " << is_Coriolis << endl;
    cout << "Is water drag: " << is_water_drag << endl; 
    cout << "Is air drag: " << is_air_drag << endl;
    cout << "Momentum boundary condition: " << MomentumBCToName[boundary_conditions] << endl;
};