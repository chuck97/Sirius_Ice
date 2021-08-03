#include <math.h>
#include "common.h"
#include "tree3.h"
#include "memory3.h"

double distanceS(double x,double y,double z, double xc,double yc,double zc) {
	return (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc);
} /*distanceS*/

PStrucFace3 addFace(surface_mesh *pm, int v1, int v2, int v3, int twin, int color) {
    PStrucFace3 face;

    myReallocv(pm->nFace+1, &pm->maxFace, 1024, 1, (void**)&pm->face, sizeof(StrucFace3));

    face = &pm->face[pm->nFace];
    face->color = color;
    face->v1 = v1;
    face->v2 = v2;
    face->v3 = v3;
    pm->nFace++;

    return face;
    (void) twin;
} /*addFace*/

