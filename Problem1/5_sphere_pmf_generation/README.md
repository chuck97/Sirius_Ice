Here is code that creates the **.pmf (internal INMOST format) and .vtk (universal format for unstructured meshes) for generated sphere triangulation**. 

Before launching make sure that you have built INMOST library. 
1) You need to specify location of INMOST SOURSE (line 10) and INMOST BUILD (line 11) directories in CMakeLists.txt. 
2) Build the code in "build" folder
``` 
cd build 
cmake .. 
cmake --build .
``` 
3) After building the code you need to specify configuration file (config.txt) 

first line: specify full path to nodes (smth+vrt.txt)

second line: specify full path to triangles (smth + vrt.txt)

third line: specify output folder for .pmf and .vtk file containing final sphere triangulation (I recomend use MESH folder in root/Problem1)

fourth line: specify postfix (I recomend use mesh resolution)
 
4) run the code and provide it with configuration file, e.g.
``` 
./GRID_TO_PMF ../config.txt
``` 

The output of the program:

sphere + specified postfix +.vtk in specified folder (you can visualize the mesh using ParaView); 

sphere + specified postfix +.pmf in specified folder (internal INMOST mesh format - will be used further).
