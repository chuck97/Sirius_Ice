#include "config.h"

using namespace std;

OutputParams::OutputParams()
{};

OutputParams::OutputParams(const std::string& json_path)
{
    string no_spaces_json = json_path;
    rtrim(no_spaces_json);
    if(no_spaces_json.substr(no_spaces_json.size()-5, 5) != ".json")
    {
        INMOST_ICE_ERR("input file shoud be ended by .json!");
    }
    std::ifstream ifs(no_spaces_json);
    nlohmann::json j_input = nlohmann::json::parse(ifs);
    BARRIER

    // Parse pvtu directory
    if (j_input.count("pvtu output directory") != 0)
    {
        is_pvtu_output = true;
        string output_pvtu_dir = j_input["pvtu output directory"];
        rtrim(output_pvtu_dir);
        output_pvtu_directory = output_pvtu_dir;
    }

    // Parse netcdf directory
    if (j_input.count("netcdf output directory") != 0)
    {
        is_netcdf_output = true;
        string output_netcdf_dir = j_input["netcdf output directory"];
        rtrim(output_netcdf_dir);
        output_netcdf_directory = output_netcdf_dir;
    }
    BARRIER

    // Parse errors file
    if (j_input.count("errors output file") != 0)
    {
        is_errors_output = true;
        string output_errors_fil = j_input["errors output file"];
        rtrim(output_errors_fil);
        output_errors_file = output_errors_fil;
    }
    BARRIER

    // Parse output gap
    if (j_input.count("every nth screenshot") != 0)
    {
        every_nth_screenshot = j_input["every nth screenshot"];
    }
    BARRIER

    // Parse verbosity of output
    if (j_input["verbose output"] != 0)
    {
        is_verbose_output = j_input["verbose output"];
    }
    BARRIER

    // fill all displayed variables to false
    for (auto item : ModelVariableNotationList)
    {
        displayed_variables[item] = false;
    }
    
    // Parse displayed model variables
    if (is_pvtu_output or is_netcdf_output)
    {
        if (j_input["displayed variables"].empty())
        {
            INMOST_ICE_ERR("should be at least one displayed variable!");
        }
        else
        {
            vector<string> variable_names; 
            for (auto& item: j_input["displayed variables"])
            {
                variable_names.push_back(string(item));
            }
            
            for (auto& item: variable_names)
            {
                if (ModelVariableNameToNotation.count(item) == 0)
                {
                    INMOST_ICE_ERR("variable " + item + " doesn't exist");
                }
                else
                {
                    displayed_variables[ModelVariableNameToNotation[item]] = true;
                }
            }
        }
    }
};

string OutputParams::GetOutputPvtuDirectory() const
{
    return output_pvtu_directory;
};

string OutputParams::GetOutputNetcdfDirectory() const
{
    return output_netcdf_directory;
};

string OutputParams::GetOutputErrorsFile() const
{
    return output_errors_file;
};

int OutputParams::GetEveryNthScreenshot() const
{
    return every_nth_screenshot;
};

std::map<ModelVariableNotation, bool> OutputParams::GetDisplayedVariables() const
{
    return displayed_variables;
};

bool OutputParams::GetIsVerposeOutput() const
{
    return is_verbose_output;
};

bool OutputParams::GetIsPvtuOutput() const
{
    return is_pvtu_output;
};

bool OutputParams::GetIsErrorsOutput() const
{
    return is_errors_output;
};

bool OutputParams::GetIsNetcdfOutput() const
{
    return is_netcdf_output;
};

void OutputParams::Log() const
{
    if (is_pvtu_output)
    {
        cout << "output pvtu dir: " << output_pvtu_directory << endl;
    }

    if (is_netcdf_output)
    {
        cout << "output netcdf dir: " << output_netcdf_directory << endl;
    }

    if (is_errors_output)
    {
        cout << "output errors file: " << output_errors_file << endl;
    }
    
    cout << "every nth screenshot will be printed: " << every_nth_screenshot << endl;
    
    if (is_verbose_output)
    {
        cout << "Verbose output" << endl;
    }

    cout << "list of displayed variables: ";
    for (auto [key, value]: displayed_variables)
    {
        if (value)
        {
            cout << ModelVariableNotationToName[key] << ", ";
        }
    }
    cout << endl;
};