/* surface meshing */
int set_surface_size_function(double (*f) (double, double, double));
int set_param_functions(int (*bounsurf) (int, double, double, double*, double*, double*), void (*v_u) (int, double, double*));

//int surface_boundary (int   nVVert, double *VVertxyz, int   nLine, int *LineD, int *LineP, double *LineT, int   nSurface, int *SurfL, int *SurfI, double *SurfT, int *pnV, double *vertex, int *pnF, int *face, int *facecolor, int   maxnV, int   maxnF);
int surface_boundary_(int *pnVVert, double *VVertxyz, int *pnLine, int *LineD, int *LineP, double *LineT, int *pnSurface, int *SurfL, int *SurfI, double *SurfT, int *pnV, double *vertex, int *pnF, int *face, int *facecolor, int *pmaxnV, int *pmaxnF);

/* polyhedrons */
int polyhedron_add_faceloop(int *pnE, int *edge, int nP, ...);
int polyhedron_make_front (int *pnV, double *vertex, int nS, int *nE, int **edge, int *color1, int *color2, int *pnF, int *face, int *facematerial, int nnV, int nnF);
