#define ani3d_surface_cgm_model ani3d_surface_cgm_model_

int ani3d_surface_CGM_model_0 (
	CGMmodel model,  /* CGM model */
	double (*fsize) (double, double, double),  /* Size function for surface */
	int *pnV, double *vertex,  /* Mesh vertices */
	int *pnF, int *face, int *facecolor,  /* Mesh triangles */
	int maxnV, int maxnF  /* Size of output arrays */
);
int ani3d_surface_cgm_model (
	CGMmodel model,  /* CGM model */
	double (*fsize) (double, double, double),  /* Size function for surface */
	int *pnV, double *vertex,  /* Mesh vertices */
	int *pnF, int *face, int *facecolor,  /* Mesh triangles */
	int *pmaxnV, int *pmaxnF  /* Size of output arrays */
);

/* for compatibility with Ani3D-2.3 */

int surface_CGM_model(CGMmodel model, int *pnV, double *vertex, int *pnF, int *face, int *facecolor, int maxnV, int maxnF);
int surface_CGM_model_(CGMmodel model, int *pnV, double *vertex, int *pnF, int *face, int *facecolor, int *pmaxnV, int *pmaxnF);
