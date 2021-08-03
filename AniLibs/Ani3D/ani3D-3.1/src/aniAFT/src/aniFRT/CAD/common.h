#ifndef AFT25D_COMMON
#define AFT25D_COMMON included

#ifndef CGMWRAP_H
#define CGMface void*
#endif

typedef struct{
	double    x,y,z;           /*   */
	double    u,v;             /*   */
}  StrucVert3;
typedef  StrucVert3  *PStrucVert3;

typedef struct {
    int *ia;
    int *ja;
} neigh_edges;

typedef struct{
	int     v1,v2;   /**/
	int     f;          /* number  of  face       */
	double  x,y,z,s;    /**/
	int     color;
	int     fail;
}  StrucEdge3;
typedef  StrucEdge3  *PStrucEdge3;

typedef struct{
	int     v1,v2,v3;   /**/
	int     color;
}  StrucFace3;
typedef  StrucFace3  *PStrucFace3;

typedef struct _toctreenode{ /* Octree node */
    int ne, nn;
    int *e;
    struct _toctreenode *child[8];
} octreenode;

typedef struct { /* edge3d */
    int v[2];
    double x[3];
    double r;
    double s;
    int h;
    octreenode *node;
} edge3d;

typedef struct { /* mesh32 */
    int  nE, nnE;
    edge3d *e;
    int  *heap;
    int  npack,  *pack;
    double bb[6];
    octreenode *root;
} mesh32;

typedef struct { /* curve */
    int nV, nnV;
    double *vertex;
    int nE, nnE, *edge;
    double (*fsize) (double *x, void *pdsize);
    int (*fparam) (double t, double *x, void *pdparam);
    void *pdsize, *pdparam;
    int v0, v1;
    double t0, t1;
    unsigned char state; /* 0 – void, 15 – params defined, 128 – mesh constructed */
} curve_mesh;

typedef struct {
    void (*v_u) (int, double, double*);  /* CRV parametric function for edges */
    int (*bounsurf) (int, double, double, double*, double*, double*);  /* CRV parametric function for faces */
    int (*bounline) (int, double, double*, double*);  /* CRV parametric function for edges */
    double (*fsize) (double, double, double);  /* Size function for surface */
    double (*periodicfunction) (int, int);  /* Periodic parametization */
    int (*cgmfunction) (int, double, double, double*, double*, double*);  /* CRV parametric function from CGM */
    int (*cgmsurfNormal) (double,  double, double, double *);  /* Get normal function from CGM */

    double x0, y0, z0, x1, y1, z1, x2, y2, z2;  /* Coeff. for plane parametrization */

    double S0;  /* Default surface mesh size */

    /* From region3.h : reg3 */
    double tBegin,tEnd;  /* Parametric space */
    double uMin,uMax,vMin,vMax,uSave,vSave;
    double suMin,suMax,svMin,svMax;
    int *tria, nTria, maxTria;
    int vv0,vv1,vv2;
    int iSurf,iNorm;

    /* From struct3.h : mesh3 */

    int     maxPoint;             /**/
    int     nPoint;    
    PStrucVert3   vert;

    int     nLinePoint;
    neigh_edges eadj;
   
    /* From tree3.h : tree3 */
    
    PStrucFace3 face; /*  array of  the  faces  in  advanced  front */
    int         nFace;
    int         maxFace; /*  number  &  max  ...  of  faces  in  ... */

    double	xc,yc,zc,side;	/*  global  for  recursive  remove  &  insert  function  */
    double	x,y,z;	/*  global  for  recursive  remove  &  insert  function  */
    double      boxcx, boxcy, boxcz;
    double      boxsize;


    /* From section3.h : sect3 */

    StrucVert3   sv;
    int          sworkVert;        /*   */
    int          sneigBool[3];     /*   */
    PStrucEdge3  sneigFace[3];     /*   */
    PStrucEdge3  swork;            /*   */
    PStrucEdge3  ssect;            /*   */


    /* tree32 section */

    PStrucEdge3	*edge;	/*  array of  the  faces  in  advanced  front */
    int		nEdge;
    int		nnEdge;	/*  number  &  max  ...  of  faces  in  ... */

    PStrucEdge3	*tvicinityFace;	/*  array of  faces  the  vicinity  */
    int		tnVicinityFace,tmaxVicinityFace;	/*  number  &  max  ...  of  faces  in  ... */



    /* From tetra3.c */

    int nBadVert,  nSphereVert,  *badVert,  *sphereVert,  maxBadVert,  maxSphereVert;

    /* From tree32.c */
    mesh32 *gmesh;

    double cf, sizelim;

    /* From cgmwrap.h */
    CGMface cgmface;

    /* New section */
    int fresh;
    StrucVert3 gp;

    int *exportCurves, nexpEdge, *expEdge, *expEdgeColor, nnexpEdge;

} surface_mesh;

#endif
