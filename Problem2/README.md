The second problem is to test different types of Finite Element Taylor-Galerkin class advection schemes on sphere.

The folder consists of:

1) *./Src* - the folder with source code for advection test;
2) *Problem2.pdf* - the description of problem (will be discussed in class);
3) *config.json* - the example of configuration file;
4) *launcher.qs* - the example of launching file for slurm queue;
5) *PlanePics.ipynb* - the example of python code (should be opened with jupyter notebook) that makes plane pictures (you can use it, or write your own).

Before building the code make sure, that you have built INMOST with Petsc/Parmetis support!

To build the code, edit your *./Src/CMakeLists.txt*. The instructions of running the code are given in *Problem2.pdf*.


The description of the problem is presented in the end of *Problem2.pdf*.
