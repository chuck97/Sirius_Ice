Here is code that performs **sphere tetration** using Ani3D lib. 

Before launching make sure that you have built Ani3D library.

1) You need to specify location of ani3d library folder in CMakeLists.txt (line 7);

3) Build code in "build" directory
``` 
cd build
cmake ..
cmake --build .
```
3) To run the program you need to pass the triangle size in km and postfix (e.g.
```
./SPHERE_TETRATION 200.0 200km
``` 
would make Earth sphere tetration with 200km edges size and produce file mesh200km.out in build directory). 

**The output of the program**:

"mesh" + prefix + ".out"  - the internal Ani3D format with tetration;

