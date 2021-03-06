Here is code that creates the **.pmf (internal INMOST format)** and **.vtk (universal format for unstructured meshes) for generated mesh**. 

Before launching make sure that you have built INMOST library.

1) You need to specify location of INMOST SOURSE (line 10) and INMOST BUILD (line 11) directories in CMakeLists.txt. 

3) To build the code go to the build folder (cd build), configure cmake (cmake ..) and build project (cmake --build .):
``` 
cd build
cmake ..
cmake --build .
``` 
5) After building the code you need to specify configuration file (config.txt)

*first line*: full path to vrt.txt file containing the nodes coordinates (after whitespace); 

*second line*: full path to tri.txt file containing the triangles adjacency list (after whitespace);

*third lin*e: output folder for generated .vtk and .pmf file. I recomend to use MESHES folder in root/Problem1/MESHES;

*fourth line*: specify the postfix for .vtk and .pmf files (e.g. for 8km grid good prefix is "8km" so your files would be named like "square8km.vtk" and "square8km.pmf"). Also whitespace after ":" should be made. 
 
4) run code and provide .exe with configuration file 
``` 
./GRID_TO_PMF ../config.txt)
``` 

**The output of the program**: 

square + specified postfix +.vtk in specified folder (you can visualize the mesh using ParaView); 

square + specified postfix +.pmf in specified folder (internal INMOST mesh format - will be used further).
