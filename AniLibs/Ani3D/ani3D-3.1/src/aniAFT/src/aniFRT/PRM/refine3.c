#include <math.h>
#include "common.h"
#include "memory3.h"
#include "user3.h"
#include "error3.h"
#include "region3.h"
#include "refine3.h"

static double distance(double x,double y,double z, double xc,double yc,double zc) {
	return sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) );
} /*distance*/

typedef struct {
    int p;
    int n;
} plist;

static int add_glist(int a, int b, int *png, plist *glist, int *s) {
    int c = s[a];
    while (c>=0) {
	if (glist[c].p == b) return 0;
	c = glist[c].n;
    }
    glist[*png].p = b;
    glist[*png].n = s[a];
    s[a] = *png;
    (*png)++;
    return 1;
}

static int fill_eadj(surface_mesh *pm) {
    int *s;
    int ng;
    plist *glist;
    int i, j, n, l;
    int a, b, c;
    
    glist = (plist*)myAlloc(sizeof(plist)*6*pm->nTria); /* Safe allocation, same function */
    ng = 0;
    s = (int*)myAlloc(sizeof(int)*pm->nPoint); /* Safe allocation, same function */
    for (j=0; j<pm->nPoint; j++)  s[j] = -1;
    
    for (i=0; i<pm->nTria; i++) {
	a = pm->tria[3*i+0],  b = pm->tria[3*i+1],  c = pm->tria[3*i+2];
	add_glist(a, b, &ng, glist, s);
	add_glist(a, c, &ng, glist, s);
	add_glist(b, c, &ng, glist, s);
	add_glist(b, a, &ng, glist, s);
	add_glist(c, a, &ng, glist, s);
	add_glist(c, b, &ng, glist, s);
    }
    n = 0;
    for (j=0; j<pm->nPoint; j++) {
	pm->eadj.ia[j] = n;
	l = s[j];
	while (l>=0) {
	    pm->eadj.ja[n++] = glist[l].p;
	    l = glist[l].n;
	}
    }
    pm->eadj.ia[pm->nPoint] = n;
    
    free(s),  free(glist);
    return n;
}


/* Simple gaussian smoothing */
void smoothingSurf(surface_mesh *pm) {
    int    i, j, k, n;
    double x0, y0, z0, xx, yy, zz, x, y, z, u, v;
    double f1, du, dv, du1, dv1, du2, dv2, maxdu, maxdv, d1u, d1v, size;

    pm->eadj.ia   = (int*)myAlloc(sizeof(int )*(pm->nPoint + 1)); /* Safe allocation, same function */
    pm->eadj.ja   = (int*)myAlloc(sizeof(int )*(6*pm->nTria)); /* Safe allocation, same function */

    fill_eadj(pm);

    for (i=pm->nLinePoint; i<pm->nPoint; i++) {
	n = pm->eadj.ia[i+1] - pm->eadj.ia[i];
	if (n <= 0)  continue;
	xx = 0.0,  yy = 0.0,  zz = 0.0;
	for (j=pm->eadj.ia[i]; j<pm->eadj.ia[i+1]; j++) {
	    k = pm->eadj.ja[j];
            x = pm->vert[k].x;
            y = pm->vert[k].y;
            z = pm->vert[k].z;
            xx += x;
            yy += y;
            zz += z;
	}
        xx /= n;
        yy /= n;
        zz /= n;

        u = pm->vert[i].u;
        v = pm->vert[i].v;
        bounSurf(pm, pm->iSurf, u, v, &x, &y, &z);
        size = sizeFace(pm, x, y, z);
        f1 = distance(xx, yy, zz, x, y, z);
        for (j=0; j<15; j++) { /* NEW METH */
            du = 0.2*(pm->uMax - pm->uMin);
            do {
                du /= 2.0;
                bounSurf(pm, pm->iSurf, u+du, v, &x0, &y0, &z0);
            } while (distance(x, y, z, x0, y0, z0) > 0.01*size);
            d1u = (distance(xx, yy, zz, x0, y0, z0) - f1)/du;

            dv = 0.2*(pm->vMax - pm->vMin);
            do {
                dv /= 2.0;
                bounSurf(pm, pm->iSurf, u, v+dv, &x0, &y0, &z0);
            } while (distance(x, y, z, x0, y0, z0) > 0.01*size);
            d1v = (distance(xx, yy, zz, x0, y0, z0) - f1)/dv;

            du = 0.6*(pm->uMax - pm->uMin);
            do {
                du /= 2.0;
                bounSurf(pm, pm->iSurf, u+du, v, &x0, &y0, &z0);
            } while (distance(x, y, z, x0, y0, z0) > 0.1*f1);

            dv = 0.6*(pm->vMax - pm->vMin);
            do {
                dv /= 2.0;
                bounSurf(pm, pm->iSurf, u, v+dv, &x0, &y0, &z0);
            } while (distance(x, y, z, x0, y0, z0) > 0.1*f1);

            maxdu = du;
            maxdv = dv;
            if (d1u == 0. && d1v == 0)
                errorExit3(3, " der == 0.   in  smooth ");
            if (fabs(d1u) > fabs(d1v)) {
                dv1 = maxdv;
                du1 = (-f1 - d1v*dv)/d1u;
                dv2 = -maxdv;
                du2 = (-f1 - d1v*dv)/d1u;
                if (fabs(du1) < fabs(du2)) {
                    du = du1;
                    dv = dv1;
                } else {
                    du = du2;
                    dv = dv2;
                }
                if (du > maxdu) {
                    dv *= (maxdu/du);
                    du = maxdu;
                }
                if (du < -maxdu) {
                    dv *= (-maxdu/du);
                    du = -maxdu;
                }
            } else {
                du1 = maxdu;
                dv1 = (-f1 - d1u*du)/d1v;
                du2 = maxdu;
                dv2 = (-f1 - d1u*du)/d1v;
                if (fabs(dv1) < fabs(dv2)) {
                    du = du1;
                    dv = dv1;
                } else {
                    du = du2;
                    dv = dv2;
                }
                if (dv > maxdv) {
                    du *= (maxdv/dv);
                    dv = maxdv;
                }
                if (dv < -maxdv) {
                    du *= (-maxdv/dv);
                    dv = -maxdv;
                }
            }
            u += du;
            v += dv;
            bounSurf(pm, pm->iSurf, u, v, &x, &y, &z);
            f1 = distance(xx, yy, zz, x, y, z);
        }/*NEW METH*/
        pm->vert[i].u = u;
        pm->vert[i].v = v;
        pm->vert[i].x = x;
        pm->vert[i].y = y;
        pm->vert[i].z = z;
    }/* for i */

    free(pm->eadj.ia),  free(pm->eadj.ja);

    return;
} /*smoothingSurf*/


