#ifndef H_REGION3_MESH3D
#define H_REGION3_MESH3D

#ifndef CGMWRAP_H
#define CGMface void*
#endif

typedef struct{
	int     iSurf,iV_U;
	double  tBegin,tEnd,*u,*v;
}  StrucSub;

typedef struct{
	int      nVert,*vert;
	int      vBegin,vEnd;
	int      nSub;
	StrucSub *sub;
}  StrucLine3;

typedef struct{
	int     nLine,*line,*inverse;
	int     iSurf,iNorm,iLabel,bCut;
	double  uMin,uMax,vMin,vMax;
}  StrucSurface;

/* exported  functions */
double sizeFace(surface_mesh *pm,  double x, double y, double z );

void initAFS_ (
	surface_mesh *pm,
	int *pnVVert, double *VVertxyz,
	int *pnLine, int *LineD, int *LineP, double *LineT,
	int *pnSurface, int *SurfL, int *SurfI, double *SurfT
);
#endif

