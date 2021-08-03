Here is code that performs converts sphere tetration in Ani3D internal format to .txt files.

1) You need to specify configuration file (./config.txt):

first line: path to previously generated .out file

second line: postfix to the file

third line: output folder for generated .txt files (I recommend use full path to build folder in current directory)

2) Build code in "build" directory
``` 
cd build
cmake ..
cmake --build .
```
3) To run the program you need to pass the path to config.txt file, e.g.:
```
./OUT_TO_TXT ../config.txt
``` 

The output of the program:

"vrt" + prefix + ".txt"  - file containing node coordinates;

"tri" + prefix + ".txt"  - file containing the adjacency list for triangles;

