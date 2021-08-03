#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <map>

#define M_Assert(Expr, Msg) \
    if (!Expr)  \
    { \
    std::cerr << "Assert failed:\t" << Msg << "\n" \
        << "Expected:\t" << #Expr << "\n" \
        << "Source:\t\t" << __FILE__ << ", line " << __LINE__ << "\n"; \
        abort(); \
    }

void WriteTrInfoToFile(std::vector<std::vector<double>>& verts,
					   std::vector<std::vector<int>>& trians,
					   const std::string& vrtoutput,
					   const std::string& trianoutput)
{
	std::ofstream output_vrt(vrtoutput);
	std::ofstream output_tri(trianoutput);
	
	size_t nv = verts.size();
	size_t nt = trians.size();

	int tri[3], bnd[2];

	// vrt coords
	output_vrt << "number of verts: " << nv << std::endl;
	output_vrt.precision(10);
	
	for (size_t i = 0; i < nv; ++i)
	{
		output_vrt << std::fixed;
		output_vrt << verts[i][0] << ' '
				   << verts[i][1] << ' '
				   << verts[i][2] << std::endl;
	}
	
	// tri connenctivity list
	output_tri << "number of triangles: " << nt << std::endl;
	for (size_t i = 0; i < nt; ++i)
	{
		output_tri << trians[i][0] << ' '
				   << trians[i][1] << ' '
				   << trians[i][2] << std::endl;
	}
}


void ReadFromOutWriteTxt(const std::string& filepath,
						 const std::string& vrtoutput,
						 const std::string& trianoutput)
{
	std::ifstream file;
    file.open(filepath);
	M_Assert((file.is_open()),"cannto open file to read");
	
    //read all verts
    int nv;
    file >> nv;
    std::cout << "nv = " << nv << std::endl;
    std::vector<std::vector<double>> all_verts(nv);
    std::vector<int> new_pos(nv);
    std::vector<bool> is_vert(nv);
    for (int i = 0; i < nv; ++i)
    {
		std::vector<double> tmp_v(3);
		file >> tmp_v[0];
		file >> tmp_v[1];
		file >> tmp_v[2];
		
		all_verts[i] = tmp_v;
		is_vert[i] = false;
	}
	
	// skip all tetrahedrons
	int n_th;
	file >> n_th;
	std::cout << "n_th = " << n_th << std::endl;
	int foo_int;
	for(int i = 0; i < n_th; ++i)
	{
		file >> foo_int;
		file >> foo_int;
		file >> foo_int;
		file >> foo_int;
		file >> foo_int;
	}
	
	// read all boundary triangles
	int nessv = 0;
	int nt;
	file >> nt;
	std::cout << "nt = " << nt << std::endl;
	
	std::vector<std::vector<double>> verts;
	
	std::vector<std::vector<int>> trians;
	
	for (int i = 0; i < nt; ++i)
    {
		std::vector<int> current_trian(3);
		std::vector<int> tmp_t(3);
		file >> tmp_t[0];
		file >> tmp_t[1];
		file >> tmp_t[2];
		file >> foo_int;
		
		// 1
		if (!is_vert[tmp_t[0] - 1])
		{
			verts.push_back(all_verts[tmp_t[0] - 1]);
			new_pos[tmp_t[0] - 1] = nessv + 1;
			is_vert[tmp_t[0] - 1] = true;
			current_trian[0] = nessv + 1;
			++nessv;
		}
		else
		{
			current_trian[0] = new_pos[tmp_t[0] - 1];
		}
		
		// 2
		if (!is_vert[tmp_t[1] - 1])
		{
			verts.push_back(all_verts[tmp_t[1] - 1]);
			new_pos[tmp_t[1] - 1] = nessv + 1;
			is_vert[tmp_t[1] - 1] = true;
			current_trian[1] = nessv + 1;
			++nessv;
		}
		else
		{
			current_trian[1] = new_pos[tmp_t[1] - 1];
		}
		
		// 3
		if (!is_vert[tmp_t[2] - 1])
		{
			verts.push_back(all_verts[tmp_t[2] - 1]);
			new_pos[tmp_t[2] - 1] = nessv + 1;
			is_vert[tmp_t[2] - 1] = true;
			current_trian[2] = nessv + 1;
			++nessv;
		}
		else
		{
			current_trian[2] = new_pos[tmp_t[2] - 1];
		}
		
		trians.push_back(current_trian);
	}
	
	file.close();
	
	std::cout << "there are " << nessv << " verts and " << nt << " trians" << std::endl;
	
	 WriteTrInfoToFile(verts, trians, vrtoutput, trianoutput);
}

int main(int argc, char* argv[]) 
{	
	std::string current_exec_name = argv[0]; 
    std::vector<std::string> all_args;

    std::string config_file;

    if (argc > 1) 
    {
        all_args.assign(argv + 1, argv + argc);
    }

    if (all_args.size() > 1)
    {
        std::cout << "should be only one *.txt file on input" << std::endl;
		std::terminate();
    }
    else
    {
        config_file = all_args.back();
    }
    
    std::ifstream conf(config_file);
    
    std::string line;
	
	// parse .out file
	std::string out_file; 
	getline(conf, line);
	std::stringstream ss1(line);
	ss1 >> line;
	ss1 >> out_file;
	
	// parse postfix
	std::string postfix; 
	getline(conf, line);
	std::stringstream ss2(line);
	ss2 >> line;
	ss2 >> postfix;
	
	// parse output folder
	std::string output_folder; 
	getline(conf, line);
	std::stringstream ss3(line);
	ss3 >> line;
	ss3 >> output_folder;
	
    ReadFromOutWriteTxt(out_file, output_folder+"vrt"+postfix+".txt",
								  output_folder+"tri"+postfix+".txt");
}
