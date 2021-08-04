#include "config.h"

using namespace std;

MeshParams::MeshParams()
{};

MeshParams::MeshParams(const string& json_path)
{
    std::string no_spaces_json = json_path;
    rtrim(no_spaces_json);
    if(no_spaces_json.substr(no_spaces_json.size()-5, 5) != ".json")
    {
        INMOST_ICE_ERR("input file shoud be ended by .json!");
    }
    std::ifstream ifs(no_spaces_json);
    nlohmann::json j_input = nlohmann::json::parse(ifs);
    BARRIER
    // Parse mesh path
    if (j_input["mesh"].empty())
    {
        INMOST_ICE_ERR("mesh path should be given (\"mesh\": \"path_to_mesh.pmf\")!");
    }
    else
    {
        std::string m_path = j_input["mesh"];
        rtrim(m_path);
        mesh_path = m_path;
    }
    BARRIER

    // Parse surface type
    if (j_input["surface type"].empty())
    {
        INMOST_ICE_ERR("surface type should be given!");
    }
    else
    {
        std::string s_t = j_input["surface type"];
        if (SurfaceNameToType.count(s_t) == 0)
        {
            INMOST_ICE_ERR("no such surface type:" + s_t);
        }
        else
        {
            surface_type = SurfaceNameToType[s_t];
        }   
    }
    BARRIER

    // Parse coords type
    if (j_input["coords type"].empty())
    {
        INMOST_ICE_ERR("coords type should be given!");
    }
    else
    {
        std::string c_t = j_input["coords type"];
        if (CoordsNameToType.count(c_t) == 0)
        {
            INMOST_ICE_ERR("no such coords type:" + c_t);
        }
        else
        {
            coords_type = CoordsNameToType[c_t];
        }   
    }
    BARRIER
};

string MeshParams::GetMeshPath() const
{
    return mesh_path; 
};

CoordsType MeshParams::GetCoordsType() const
{
    return coords_type; 
};

SurfaceType MeshParams::GetSurfaceType() const
{
    return surface_type; 
};

void MeshParams::Log() const
{
    cout << "Mesh: " << mesh_path << endl;
    cout << "surface type: " << SurfaceTypeToName[surface_type] << endl;
    cout << "coords type: " << CoordsTypeToName[coords_type] << endl;
}