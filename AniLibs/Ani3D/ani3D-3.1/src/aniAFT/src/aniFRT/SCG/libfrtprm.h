#define ani3d_surface_boundary ani3d_surface_boundary_
#define ani3d_surface_edges_boundary ani3d_surface_edges_boundary_

/* surface meshing */

int ani3d_surface_boundary_0 (
	int nVVert, double *VVertxyz,  /* V vertices */
	int nLine, int *LineD, int *LineP, double *LineT,  /* CRV edges */
	int nSurface, int *SurfL, int *SurfI, double *SurfT,  /* CRV faces */
	void (*v_u) (int, double, double*),  /* CRV parametric function for edges */
	int (*bounsurf) (int, double, double, double*, double*, double*),  /* CRV parametric function for faces */
	int (*bounline) (int, double, double*, double*),  /* CRV parametric function for edges */
	double (*periodicfunction) (int, int),  /* Periodic parametization */
	double (*fsize) (double, double, double),  /* Size function for surface */
	int *pnV, double *vertex,  /* Mesh vertices */
	int *pnF, int *face, int *facecolor,  /* Mesh triangles */
	int maxnV, int maxnF  /* Size of output arrays */
);
int ani3d_surface_boundary (
	int *pnVVert, double *VVertxyz,  /* V vertices */
	int *pnLine, int *LineD, int *LineP, double *LineT,  /* CRV edges */
	int *pnSurface, int *SurfL, int *SurfI, double *SurfT,  /* CRV faces */
	void (*v_u) (int, double, double*),  /* CRV parametric function for edges */
	int (*bounsurf) (int, double, double, double*, double*, double*),  /* CRV parametric function for faces */
	int (*bounline) (int, double, double*, double*),  /* CRV parametric function for edges */
	double (*periodicfunction) (int, int),  /* Periodic parametization */
	double (*fsize) (double, double, double),  /* Size function for surface */
	int *pnV, double *vertex,  /* Mesh vertices */
	int *pnF, int *face, int *facecolor,  /* Mesh triangles */
	int *pmaxnV, int *pmaxnF  /* Size of output arrays */
);

int ani3d_surface_edges_boundary (
	int *pnVVert, double *VVertxyz,  /* V vertices */
	int *pnLine, int *LineD, int *LineP, double *LineT,  /* CRV edges */
	int *pnSurface, int *SurfL, int *SurfI, double *SurfT,  /* CRV faces */
	void (*v_u) (int, double, double*),  /* CRV parametric function for edges */
	int (*bounsurf) (int, double, double, double*, double*, double*),  /* CRV parametric function for faces */
	int (*bounline) (int, double, double*, double*),  /* CRV parametric function for edges */
	double (*periodicfunction) (int, int),  /* Periodic parametization */
	double (*fsize) (double, double, double),  /* Size function for surface */
	int *exportCurves,  /* Curve colors */
	int *pnV, double *vertex,  /* Mesh vertices */
	int *pnF, int *face, int *facecolor,  /* Mesh triangles */
	int *pnE, int *edge, int *edgecolor,  /* Mesh edges */
	int *pmaxnV, int *pmaxnF, int *pmaxnE  /* Size of output arrays */
);


/* polyhedrons */

int ani3d_polyhedron_add_faceloop(int *pnE, int *edge, int nP, ...);
int ani3d_polyhedron_make_front (
	int *pnV, double *vertex,
	int nS, int *nE, int **edge, int *color1, int *color2,
	double (*fsize) (double, double, double),
	int *pnF, int *face, int *facematerial,
	int nnV, int nnF
);


/* for compatibility with Ani3D-2.3 */

/* surface meshing */
int set_surface_size_function(double (*f) (double, double, double));
int set_param_functions(int (*bounsurf) (int, double, double, double*, double*, double*), void (*v_u) (int, double, double*));

//int surface_boundary (int   nVVert, double *VVertxyz, int   nLine, int *LineD, int *LineP, double *LineT, int   nSurface, int *SurfL, int *SurfI, double *SurfT, int *pnV, double *vertex, int *pnF, int *face, int *facecolor, int   maxnV, int   maxnF);
int surface_boundary_(int *pnVVert, double *VVertxyz, int *pnLine, int *LineD, int *LineP, double *LineT, int *pnSurface, int *SurfL, int *SurfI, double *SurfT, int *pnV, double *vertex, int *pnF, int *face, int *facecolor, int *pmaxnV, int *pmaxnF);

/* polyhedrons */
#define polyhedron_add_faceloop ani3d_polyhedron_add_faceloop
int polyhedron_make_front (int *pnV, double *vertex, int nS, int *nE, int **edge, int *color1, int *color2, int *pnF, int *face, int *facematerial, int nnV, int nnF);
