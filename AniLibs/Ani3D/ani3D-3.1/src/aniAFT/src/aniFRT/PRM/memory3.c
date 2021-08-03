#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "common.h"
#include "error3.h"
#include "memory3.h"
#include "tree32.h"

void *myAlloc (size_t n) {
    void  *p;

    p = (void *)malloc(n);
    if (p == NULL)
	errorExit3(2,"myAlloc");

    return  p;
} /*myAlloc*/

int myReallocv(int needed, int *pmax, int default_inc, int n, ...) {
    int nmax = *pmax;
    void *np;
    void **p;
    va_list ap;
    int i;
    size_t ps;
    
    if (needed < nmax)  return 0;
    nmax = needed + default_inc;
    va_start(ap, n);
    for (i=0; i<n; i++) {
	p = va_arg(ap, void**);
	ps = va_arg(ap, size_t);
	np = realloc(*p, ps*nmax);
	if (!np) {
	    errorExit3(2,"myRealloc");
	    return -1;
	}
	*p = np;
    }
    *pmax = nmax;
    return nmax;
}

void initMemory(surface_mesh *pm) {
    pm->nPoint = 0;
    pm->nFace = 0;
    pm->gmesh = myAlloc(sizeof(mesh32));
    memset(pm->gmesh, 0, sizeof(mesh32));
    pm->cf = 0.0;
    pm->sizelim = -1.0;
    pm->badVert = NULL,  pm->nBadVert = 0,  pm->maxBadVert = 0;
    pm->sphereVert = NULL,  pm->nSphereVert = 0,  pm->maxSphereVert = 0;
    pm->tvicinityFace = NULL,  pm->tvicinityFace = 0,  pm->tmaxVicinityFace = 0;
    pm->face = NULL,  pm->nFace = 0,  pm->maxFace = 0;
    pm->vert = NULL,  pm->nPoint = 0,  pm->maxPoint = 0;
    pm->edge = NULL,  pm->nEdge = 0,  pm->nnEdge = 0;
    return;
}

void freeMemory(surface_mesh *pm) {
    int i;
    for (i=pm->nEdge-1; i>=0; i--)  remFace32(pm, pm->edge[i]);
    free(pm->vert),  pm->vert = NULL,  pm->maxPoint = 0;
    free(pm->badVert),  pm->badVert = NULL,  pm->maxBadVert = 0;
    free(pm->sphereVert),  pm->sphereVert = NULL,  pm->maxSphereVert = 0;
    free(pm->tvicinityFace),  pm->tvicinityFace = NULL,  pm->tmaxVicinityFace = 0;
    free(pm->face),  pm->face = NULL,  pm->maxFace = 0;
    free(pm->edge),  pm->edge = NULL,  pm->nnEdge = 0;
    if (pm->gmesh) {
	free(pm->gmesh->e);
	free(pm->gmesh->heap);
	free(pm->gmesh->pack);
	free(pm->gmesh);
	pm->gmesh = NULL;
    }
    return;
} /*freeMemory*/


