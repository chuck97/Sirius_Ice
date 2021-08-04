#include "config.h"

using namespace std;

AdvectionParams::AdvectionParams()
{};

AdvectionParams::AdvectionParams(AdvectionSolverType ast_, bool is_fct_, double cd_)
    : advection_solver_type(ast_),
      is_fct(is_fct_),
      fct_cd(cd_)
{};

AdvectionParams::AdvectionParams(const std::string& json_path)
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

    // Parse advection solver type
    if (j_input["type"].empty())
    {
        INMOST_ICE_ERR("type should be given!");
    }
    else
    {
        string type_name = j_input["type"];
        if (AdvectionSolverNameToType.count(type_name) == 0)
        {
            INMOST_ICE_ERR("Possible advection solver types: TG2, CG2, TTG2, TTG3, TTG4!");
        }
        else
        {
            advection_solver_type = AdvectionSolverNameToType[type_name];
        }
    }
    BARRIER

    // Parse is_fct
    if (j_input["is fct"].empty())
    {
        INMOST_ICE_ERR("\"is fct\" should be given!");
    }
    else
    {
        is_fct = j_input["is fct"];
    }
    BARRIER

    // Parse cd_fct
    if (j_input["fct cd"].empty())
    {
        INMOST_ICE_ERR("\"fct cd\" should be given!");
    }
    else
    {
        fct_cd = j_input["fct cd"];
    }
    BARRIER
};

AdvectionSolverType AdvectionParams::GetAdvectionSolverType() const
{
    return advection_solver_type;
};

bool AdvectionParams::GetIsFct() const
{
    return is_fct;
};

double AdvectionParams::GetFctCd() const
{
    return fct_cd;
};

void AdvectionParams::Log() const
{
    cout << "Advection solver type: " << AdvectionSolverTypeToName[advection_solver_type] << endl;
    cout << "Is FCT: " << is_fct << endl;
    cout << "FCT cd: " << fct_cd << endl;
}