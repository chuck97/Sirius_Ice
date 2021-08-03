C ==========================================================
      Program main
C ==========================================================
C Program reads a mesh form file data/ramp.ani and defines 
C mesh function SOL. Then it generates a metric based on this
C discrete function. Then it builds a mesh quasi-optimal in
C this metric.
C ==========================================================
c Maximum number of nodes, elements and boundary edges
      Integer   MaxP, MaxE, MaxF
      Parameter(MaxP = 30000, MaxE = 4 * MaxP, MaxF = 14000)

c Basket capacity and the maximum number of local modifications 
      Integer   MaxSkipE, MaxQItr 
      Parameter(MaxSkipE = 500, MaxQItr = 100 000) 
 
c Available memory 
      Integer   MaxWr, MaxWi 
      Parameter(MaxWr = 4 500 000, MaxWi = 6 000 000) 
 
c Desired number of elements in the final mesh 
      Integer   nEStar 
      Parameter(nEStar = 20 000) 
 
c Desired final mesh quality 
      Real*8    Quality 
      Parameter(Quality = 5D-1) 
 
C ==========================================================
C group (M)
      Real*8  XYP(3, MaxP)
      Integer IPE(4, MaxE), IPF(3, MaxF)
      Integer IPV(MaxP), lbF(MaxF), lbE(MaxE)

C Group (Dev)
      Integer IFV(MaxF), IEV(MaxE) 
      Logical flagAuto
      Integer status

C group (Q)
      Real*8  rQuality, Metric(6, MaxP)

C group (W)
      Real*8  rW(MaxWr)
      Integer iW(MaxWi)

C local variables
      Real*8  x, y, z, det

C ==========================================================
c ... load the initial mesh
      Call loadMani(
     &     MaxP, MaxF, MaxE,
     &     nP, nF, nE,
     &     XYP, IPF, IPE, lbF, lbE,
     &     nPv, nFv, nEv, IPV, IFV, IEV,
     &     iW, iW, "../data/ramp.ani")


C ... draw the initial mesh
      Call saveMgmv(
     &     nP, nF, nE, XYP, IPF, IPE, lbF, lbE, 
     &     'ramp.gmv', iW)


c ... define a piecewise-linear metric as Hessian of function
c     (x**8)*(1-x) * y*(1-y) * z*(1-z)
c
      Do n = 1, nP
         x = XYP(1, n)
         y = XYP(2, n)
         z = XYP(3, n)

         Metric(1, n) = (56*(x**6) - 72*(x**7)) * y*(1-y) * z*(1-z)
         Metric(2, n) = (x**8 - x**9) * (-2) * z*(1-z)
         Metric(3, n) = (x**8 - x**9) * y*(1-y) * (-2)

         Metric(4, n) = (8*(x**7) - 9*(x**8)) * (1-2*y) * z*(1-z)
         Metric(5, n) = (8*(x**7) - 9*(x**8)) * y*(1-y) * (1-2*z)
         Metric(6, n) = (x**8 - x**9) * (1-2*y) * (1-2*z)
 
         Call SpectralModule(Metric(1, n), det)
      End do


c ... generate adaptive mesh
      flagAuto = .TRUE.  ! default mesh generation options
      status = 1         ! forbid boundary tetrahedra
      iPrint = 1         ! average level of output information


      Call mbaNodal(
c group (M)
     &     nP, MaxP, nF, MaxF, nE, MaxE,
     &     XYP, IPF, IPE, lbF, lbE,
     &     nEStar, 
c group (Dev)
     &     nPv, nFv, nEv, IPV, IFV, IEV, 
     &     flagAuto, status,
c group (Q)
     &     MaxSkipE, MaxQItr,
     &     Metric, Quality, rQuality, 
c group (W)
     &     MaxWr, MaxWi, rW, iW,
     &     iPrint, iERR)


c ... draw the final mesh
      Call saveMgmv(
     &     nP, nF, nE, XYP, IPF, IPE, lbF, lbE, 
     &     "save.gmv", iW)


c ... save the final mesh
      Call saveMani(
     &     nP, nF, nE,  
     &     XYP, IPF, IPE, lbF, lbE,
     &     nPv, nFv, nEv, IPV, IFV, IEV,
     &     .FALSE., iW, iW, "save.ani")

      Stop
      End


