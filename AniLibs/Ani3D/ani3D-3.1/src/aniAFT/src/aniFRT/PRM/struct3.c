#include "common.h"
#include "memory3.h"

void addPoint(surface_mesh *pm, double x, double y, double z) {
    myReallocv(pm->nPoint+1, &pm->maxPoint, 1024, 1, (void**)&pm->vert, sizeof(StrucVert3));
    pm->vert[pm->nPoint].x = x;
    pm->vert[pm->nPoint].y = y;
    pm->vert[pm->nPoint].z = z;
    pm->nPoint++;
    return;
} /*addPoint*/

