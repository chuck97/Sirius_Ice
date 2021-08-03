#include "cgmwrap.h"
#include "supportc.h"
#include "libfrtcad.h"

static int ani3d_surface_CGM_model_uni (
	CGMmodel model,  /* CGM model */
	double (*fsize) (double, double, double),  /* Size function for surface */
	int *pnV, double *vertex,  /* Mesh vertices */
	int *pnF, int *face, int *facecolor,  /* Mesh triangles */
	int maxnV, int maxnF,  /* Size of output arrays */
	int shift  /* Shift indices in arrays */
) {
    return aft3dmodel(
	    model,
	    fsize,
	    pnV, vertex,
	    pnF, face, facecolor,
	    maxnV, maxnF,
	    shift
    );
}

int ani3d_surface_CGM_model_0 (
	CGMmodel model,  /* CGM model */
	double (*fsize) (double, double, double),  /* Size function for surface */
	int *pnV, double *vertex,  /* Mesh vertices */
	int *pnF, int *face, int *facecolor,  /* Mesh triangles */
	int maxnV, int maxnF  /* Size of output arrays */
) {
    return ani3d_surface_CGM_model_uni(
	    model,
	    fsize,
	    pnV, vertex,
	    pnF, face, facecolor,
	    maxnV, maxnF,
	    0
    );
}

int ani3d_surface_cgm_model (
	CGMmodel model,  /* CGM model */
	double (*fsize) (double, double, double),  /* Size function for surface */
	int *pnV, double *vertex,  /* Mesh vertices */
	int *pnF, int *face, int *facecolor,  /* Mesh triangles */
	int *pmaxnV, int *pmaxnF  /* Size of output arrays */
) {
    return ani3d_surface_CGM_model_uni(
	    model,
	    fsize,
	    pnV, vertex,
	    pnF, face, facecolor,
	    *pmaxnV, *pmaxnF,
	    1
    );
}

extern double (*ani3d_2_X_comp_internal_fsize) (double, double, double);

int surface_CGM_model (
		CGMmodel model,
		int *pnV, double *vertex,
		int *pnF, int *face, int *facecolor,
		int maxnV, int maxnF) {
    return ani3d_surface_CGM_model_uni(
	    model,
	    ani3d_2_X_comp_internal_fsize,
	    pnV, vertex,
	    pnF, face, facecolor,
	    maxnV, maxnF,
	    0
    );
}
int surface_CGM_model_(CGMmodel model, int *pnV, double *vertex, int *pnF, int *face, int *facecolor, int *pmaxnV, int *pmaxnF) {
    return ani3d_surface_CGM_model_uni(
	    model,
	    ani3d_2_X_comp_internal_fsize,
	    pnV, vertex,
	    pnF, face, facecolor,
	    *pmaxnV, *pmaxnF,
	    1
    );
}
