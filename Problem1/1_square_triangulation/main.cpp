#include <iostream>
#include <exception>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <cmath>
#include <cstring>

extern "C" 
{
	#include "ani2d_routines.h"
}


static const int    nvmax        = 10'000'000;
static const int    ntmax        = 2*nvmax;
static const int    nbmax        = 1'000'000;

template<typename ValueType>
ValueType* allocate1d(size_t size)
{
	ValueType* array = new ValueType[size];
	return array;
}

template<typename ValueType>
void deallocate1d(ValueType* array)
{
	delete[] array;
}

template<typename ValueType>
ValueType** allocate2d(size_t nrows, size_t ncols)
{
	ValueType** array = new ValueType*[nrows];
	for (size_t i = 0; i < nrows; ++i)
	{
		array[i] = new ValueType[ncols];
	}
	return array;
}

template<typename ValueType>
void deallocate2d(ValueType** array, size_t nrows)
{
	for (size_t i = 0; i < nrows; ++i)
	{
		delete[] array[i];
	} 
}

template<typename ValueType>
ValueType* reshape_2d_2_1d(ValueType** array_2d, size_t nrows, size_t ncols)
{
	ValueType* array_1d = new ValueType[nrows*ncols];
	for(size_t i = 0; i < nrows; ++i)
	{
		for (size_t j = 0; j < ncols; ++j)
		{
			array_1d[i*ncols + j] = array_2d[i][j];
		}
	}
	return array_1d;
}

template<typename ValueType>
ValueType** reshape_1d_2_2d(ValueType* array_1d, size_t nrows, size_t ncols)
{
	ValueType** array_2d = new ValueType*[nrows];
	for(size_t i = 0; i < nrows; ++i)
	{
		array_2d[i] = new ValueType[ncols];
	}
	
	for (size_t i = 0; i < nrows; ++i)
	{
		for (size_t j = 0; j < ncols; ++j)
		{
			array_2d[i][j] = array_1d[i*ncols + j];
		}
	}
	return array_2d;
}

void WriteTrInfoToFile(int nv, double* vrt, int nt, int* tri, int nb, int* bnd,
					   const std::string& vrt_output_path,
					   const std::string& tri_output_path,
					   const std::string& bnd_output_path)
{
	double** vrt2d = reshape_1d_2_2d<double>(vrt, nv, 2);
	int** tri2d    = reshape_1d_2_2d<int>   (tri, nt, 3);
	int** bnd2d    = reshape_1d_2_2d<int>	(bnd, nb, 2);
	
	std::ofstream output_vrt(vrt_output_path);
	std::ofstream output_tri(tri_output_path);
	std::ofstream output_bnd(bnd_output_path);
	
	// vrt coords
	output_vrt << "number of verts: " << nv << std::endl;
	output_vrt.precision(10);
	
	for (size_t i = 0; i < (size_t)nv; ++i)
	{
		output_vrt << std::fixed;
		output_vrt << vrt2d[i][0] << ' ' << vrt2d[i][1] << std::endl;
	}
	
	// tri connenctivity list
	output_tri << "number of triangles: " << nt << std::endl;
	
	for (size_t i = 0; i < (size_t)nt; ++i)
	{
		output_tri << tri2d[i][0] << ' ' << tri2d[i][1] << ' ' << tri2d[i][2] << std::endl;
	}
	
	// bnd edges connectivity list
	output_bnd << "number of boundary edges: " << nb << std::endl;
	
	for (size_t i = 0; i < (size_t)nb; ++i)
	{
		output_bnd << bnd2d[i][0] << ' ' << bnd2d[i][1] << std::endl;
	}
    
    deallocate2d<double>(vrt2d, nv);
    deallocate2d<int>(tri2d, nt);
    deallocate2d<int>(bnd2d, nb);
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
	
	// skip first line
	getline(input_external, line);
	
	// parse output directory
	getline(input_external, line);
	std::stringstream ss_outputdir(line);
	std::string output_folder_path;
	ss_outputdir >> line;
	ss_outputdir >> output_folder_path; 
	
	// parse triangle size
	getline(input_external, line);
	std::stringstream ss_tr_size(line);
	std::string tr_size_str;
	ss_tr_size >> line;
	double tr_size;
	ss_tr_size >> tr_size;
	
	std::vector<std::pair<double, double>> points_external;
	
	// parse domain geometry
	while(getline(input_external, line))
	{
		std::stringstream ss(line);
		double first, second;
		ss >> first >> second;
		points_external.push_back({first, second});
	}
	
	//remove last point (it is the same as first)
	points_external.pop_back();
	
	// Make input for AFT routine	
	int Nbv = points_external.size();
	int Nbl = Nbv;
	
	std::cout << "boundary verticies: " << Nbv << std::endl;
	std::cout << "boundary edges: " << Nbl << std::endl;
	
	
	double**     bv2d = allocate2d<double>(Nbv, 2);
	int**        bl2d = allocate2d<int>(Nbl, 7);
	double** bltail2d = allocate2d<double>(Nbl, 2);
	
	// fill bv and bl
	size_t k = 0;  // global index
	
	// external points
	
	for (size_t i = 0; i < points_external.size()-1; ++i)
	{
		std::pair<double, double> p = points_external[i];
		bv2d[k][0]  = p.first;
		bv2d[k][1]  = p.second;

 		bl2d[k][0]  = k+1;
 		bl2d[k][1]  = k+2;
 		bl2d[k][2]  = 0;
 		bl2d[k][3]  = -1;
 		bl2d[k][4]  =  k+1;
 		bl2d[k][5]  =  1;
 		bl2d[k][6]  =  0;
 
		bltail2d[k][0] = 0.0;
		bltail2d[k][1] = 0.0;
		++k;
	}
	
	//last point is special
	bv2d[k][0] = points_external.back().first;
	bv2d[k][1] = points_external.back().second;
	
	
	bl2d[k][0] = k+1;
	bl2d[k][1] = 1;
	bl2d[k][2] = 0;
	bl2d[k][3] = -1;
	bl2d[k][4] = k+1;
	bl2d[k][5] = 1;
	bl2d[k][6] = 0;
	
	bltail2d[k][0] = 0.0;
	bltail2d[k][1] = 0.0;
	
	++k;
	
	// reshape input
	double*     bv = reshape_2d_2_1d<double>(bv2d, Nbv, 2);
	int*        bl = reshape_2d_2_1d<int>(bl2d, Nbl, 7);
	double* bltail = reshape_2d_2_1d<double>(bltail2d, Nbl, 2);
		
	 
	//output data
	int nv, nt, nb, nc;
	
	double** vrt2d = allocate2d<double>(nvmax, 2);
	double*  vrt   = reshape_2d_2_1d<double>(vrt2d, nvmax, 2);
	
	int**    tri2d = allocate2d<int>(ntmax, 3);
	int*     tri   = reshape_2d_2_1d<int>(tri2d, ntmax, 3);
	
	int*     labelT = allocate1d<int>(ntmax);
	
	int**    bnd2d = allocate2d<int>(nbmax, 2);
	int*     bnd   = reshape_2d_2_1d<int>(bnd2d, nbmax, 2);
	  
	int*     labelB = allocate1d<int>(nbmax);
	
	double** crv2d  = allocate2d<double>(nbmax, 2);
	double*  crv    = reshape_2d_2_1d<double>(crv2d, nbmax, 2);
	
	int*     iFNC   = allocate1d<int>(nbmax);
	
	double h = tr_size;
	int ierr;
	
	//triangulation
	ierr = ani::aft2dboundary_(&Nbv, bv, &Nbl, bl, bltail, &h,   // geometry
																 // mesh data on output
                          &nv, vrt,        						 // coordinates of nodes
                          &nt, tri, labelT,   					 // triangles and their material
                          &nb, bnd, labelB,    					 // boundary edges and their labels
                          &nc, crv, iFNC);        				 // curved edges and parametrization (dummy)
  
    //possible errors
    if (ierr != 0)
    {
		std::cout << "error in function aft2dboundary, error code: " << ierr << std::endl;
		std::terminate();
	}
	
	std::cout << "mesh: number of triangles/vertices: " << nt << " " << nv << std::endl;
	
	if (nv > nvmax)
	{
		std::cout << "to many nodes" << std::endl;
		std::terminate();
	}
	
	if (nt > ntmax)
	{
		std::cout << "to many triangles" << std::endl;
		std::terminate();
	}
	
	if (nb > nbmax)
	{
		std::cout << "to many boundary edges" << std::endl;
		std::terminate();
	}
	
	std::ostringstream strs;
	strs << tr_size;
	std::string str_tr_size = strs.str();
	std::string ps_name_string = "mesh" + str_tr_size + ".ps";
	
	int n = ps_name_string.size();
    char ps_name_char[n + 1];
    strcpy(ps_name_char, ps_name_string.c_str());
	
	
	//visualization
	ani::graph_(&nv, vrt, &nt, tri, ps_name_char);
	
	WriteTrInfoToFile(nv, vrt, nt, tri, nb, bnd,
					  output_folder_path + "vrt.txt",
					  output_folder_path + "tri.txt",
					  output_folder_path + "bnd.txt");
	
	deallocate2d<double>(bv2d, Nbv);
    deallocate1d<double>(bv);
    deallocate2d<int>(bl2d, Nbl);
    deallocate1d<int>(bl);
    deallocate2d<double>(bltail2d, Nbl);
    deallocate1d<double>(bltail);
    deallocate2d<double>(vrt2d, nvmax);
    deallocate1d<double>(vrt);
	deallocate2d<int>(tri2d, ntmax);
	deallocate1d<int>(tri);
	deallocate1d<int>(labelT);
	deallocate2d<int>(bnd2d, nbmax);
	deallocate1d<int>(bnd);
	deallocate1d<int>(labelB);
	deallocate2d<double>(crv2d, nbmax);
	deallocate1d<double>(crv);
	deallocate1d<int>(iFNC);
	
    return 0;
}
