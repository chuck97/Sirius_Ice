Here is code that performs **square triangulation using Ani2D** lib. 

Before launching make sure that you have built Ani2D library. 

1) You need to specify location of ani2d library folder in Makefile (line 3); 

2) For building the code, just type: 
```
make  
```
3) Specify configuration file (square.txt): 

*second line*: output directory for triangular mesh data. It is better to make a directory like "mesh8km" (mkdir mesh8km). After that copy the full path to this directory 
```
cd mesh8km && pwd
```
and paste it after whitespace. Every folder path should be ended by "/"; 

*third line*: the size of uniform triangular grid in meters (smt like 8000.0 for 8km grid); 

*fourth - eighth line*: coordinates in meters of square verticies (we will use 512 km. square domain, the last point should be the same as first, the internal of domain should be on the right side).

3) run code and provide .exe with configuration file:
```
./triangulation.exe ./square.txt 
```
**The output of the program**: 

picture of domain (.ps file); 

vrt.txt - the file containing the coordinates of verts in meters; 

tri.txt - the file containing the adjacency list for triangles; 

bnd.txt - the file containing the adjacency list for boundary edges.
