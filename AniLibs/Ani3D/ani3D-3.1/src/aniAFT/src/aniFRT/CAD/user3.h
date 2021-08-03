#ifndef H_USER3_MESH3D
#define H_USER3_MESH3D

/* exported  function  function */
int bounSurf(surface_mesh *pm,  int i, double u, double v, double *x, double *y, double *z );
void V_U(surface_mesh *pm,  int i, double u, double *v );
int bounLine(surface_mesh *pm, int i, double t, double *u, double *v);
double periodic(surface_mesh *pm, int dir);

#endif
