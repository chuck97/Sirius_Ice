#ifndef H_MEMORY3_MESH3D  
#define H_MEMORY3_MESH3D

#include <stdlib.h>

/* exported  function  function */
void *myAlloc(size_t n);
int myReallocv(int needed, int *pmax, int default_inc, int n, ...);
void initMemory(surface_mesh *pm);
void freeMemory(surface_mesh *pm);

#endif

