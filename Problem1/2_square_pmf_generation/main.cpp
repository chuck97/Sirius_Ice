#include <utility>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "inmost.h"

using namespace INMOST;

#define M_Assert(Expr, Msg) \
    if (!Expr)  \
    { \
    std::cerr << "Assert failed:\t" << Msg << "\n" \
        << "Expected:\t" << #Expr << "\n" \
        << "Source:\t\t" << __FILE__ << ", line " << __LINE__ << "\n"; \
        abort(); \
    }



int main(int argc, char* argv[]) 
{	
	std::string current_exec_name = argv[0]; 
    std::vector<std::string> all_args;

    std::string input_external_path;

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
        input_external_path = all_args.back();
    }
    
    std::ifstream input_external(input_external_path); 
	
	std::string line;
	
	// parse vrt
	getline(input_external, line);
	std::stringstream ss_vrt(line);
	std::string name_vrt;
	ss_vrt >> line;
	ss_vrt >> name_vrt;
    
    // parse  tri
    getline(input_external, line);
	std::stringstream ss_tri(line);
	std::string name_tri;
	ss_tri >> line;
	ss_tri >> name_tri;
	
	// parse  output folder
	getline(input_external, line);
	std::stringstream ss_output(line);
	std::string name_output_folder;
	ss_output >> line;
	ss_output >> name_output_folder;
	
	// parse output postfix
	getline(input_external, line);
	std::stringstream ss_postfix(line);
	std::string name_postfix;
	ss_postfix >> line;
	ss_postfix >> name_postfix;
    
    std::cout << name_vrt << std::endl <<
				 name_bnd << std::endl <<
				 name_tri << std::endl <<
				 name_output_folder << std::endl<<
				 name_postfix << std::endl;
    
    Mesh::Initialize(&argc, &argv);

    int nv,nt,nf;

    Mesh *m;
    m = new Mesh();
    
    
    std::ifstream file;

    std::cout << "Start to read mesh!" << std::endl;

    // Read verticies
    file.open(name_vrt);
    M_Assert((file.is_open()), "can not open vrt file");

    std::string tmp_STR;

    file >> tmp_STR;
    file >> tmp_STR;
    file >> tmp_STR;
    file >> nv;
    ElementArray<Node> newverts(m);
    newverts.reserve(nv);
    std::cout << "num of verts: " << nv << std::endl;
    for (int i = 0; i < nv; i++) 
    {
        Storage::real xyz[3];

        file >> xyz[0];
        file >> xyz[1];
        xyz[2] = 0.0;

        Node nod = m -> CreateNode(xyz);
        newverts.push_back(nod);
    }
    file.close();

    // Read triangles
    file.open(name_tri);
    M_Assert((file.is_open()), "can not open tri file");
    file >> tmp_STR;
    file >> tmp_STR;
    file >> tmp_STR;
    file >> nt;
    std::cout << "num of triangles: " << nt << std::endl;
    int num_node_tmp;
    for (int i = 0; i < nt; i++) 
    {
        ElementArray<Node> verts(m);
        for (int j = 0; j < 3; j++) 
        {
            file >> num_node_tmp;
            verts.push_back(newverts[num_node_tmp - 1]);
        }

        const INMOST_DATA_INTEGER_TYPE ne_face_nodes[6] = {0, 1, 1, 2, 2, 0};
        const INMOST_DATA_INTEGER_TYPE ne_num_nodes[3] = {2, 2, 2};

        std::pair<Cell, bool> pair = m -> CreateCell(verts, ne_face_nodes, ne_num_nodes, 3);

    }
    file.close();
    
    ////////////////////////////////////// load and repartition mesh //////////////////////////
    m->Save(name_output_folder + "/mesh" + name_postfix + ".pmf");
    m->Save(name_output_folder + "/mesh" + name_postfix + ".vtk");
    //////////////////////////////////create discretization/////////////////////////////

    return 0;
}
