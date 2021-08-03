Here is code that creates the .pmf (internal INMOST format) and .vtk (universal format for unstructured meshes) for generated mesh. \

Before launching make sure that you have built INMOST library. \
1) You need to specify location of INMOST SOURSE (line 10) and INMOST BUILD (line 11) directories in CMakeLists.txt. \
2) To build the code go to the build folder (cd build), configure cmake (cmake ..) and build project (cmake --build .). \
3) After building the code you need to specify configuration file (config.txt) \

first line: full path to vrt.txt file containing the nodes voordinates (after whitespace); \
second line: full path to tri.txt file containing the triangles adjacency list (after whitespace). \
third line: output folder for generated .vtk and .pmf file. It is better to create folder for each resolution (e.g. mkdir square8km). After that copy the full path to this folder (with "/") and add it after whitespace. \
fourth line: specify the postfix for .vtk and .pmf files (e.g. for 8km grid good prefix is "8km" so your files would be named like "square8km.vtk" and "square8km.pmf"). Also whitespace after ":" should be made. \
 
4) run code and provide .exe with configuration file (./GRID_TO_PMF ../config.txt) \

The output of the program: \
square + specified postfix +.vtk on specified folder (you can visualize the mesh using ParaView); \
square + specified postfix +.pmf on specified folder (internal INMOST mesh format - will be used further).
