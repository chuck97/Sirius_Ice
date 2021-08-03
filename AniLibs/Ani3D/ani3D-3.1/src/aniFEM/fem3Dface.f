C ======================================================================
      Subroutine FEM3Dface(
C ======================================================================
     &      XY1, XY2, XY3, XY4, idface, dot_with_normal,
     &      operatorA, FEMtypeA, operatorB, FEMtypeB, 
     &      label, D, DDATA, IDATA, iSYS, order, 
     &      LDA, A, nRow, nCol)
C ======================================================================
      implicit none
      include 'fem3Dtri.fd'
C ======================================================================
C  The routine computes the elemental matrix for the surface integral 
C
C              \int_f [D OpA(u) . OpB(v)] dx
C
C  if dot_with_normal=0, or
C  
C              \int_f [(D OpA(u) . N) . OpB(v)] dx
C
C  if dot_with_normal is not zero. Here where D is a tensor function,
C  OpA and OpB are linear operators, u and v are finite element functions, 
C  u in space "A", and v in space "B", and N is the exterior normal to 
C  face f. The operations marked with "."  must be well defined.
C
C  Parameter idface is consistent with the predefined enumeration of faces
C  in the tetrahedron, see the documentation.
C
C  idface=1 denotes face {XY1, XY2, XY3}
C  idface=2 denotes face {XY2, XY3, XY4}
C  idface=3 denotes face {XY3, XY4, XY1}
C  idface=4 denotes face {XY4, XY1, XY2}
C
C  See fem3Dtet.f for detailed description of the remaining input 
C  parameters. 
C   
C ======================================================================
      Real*8   XY1(*), XY2(*), XY3(*), XY4(*)

      Real*8   DDATA(*)
      Integer  IDATA(*)
      Integer  label, iSYS(*), D

      EXTERNAL D

      Integer  FEMtypeA, FEMtypeB, operatorA, operatorB
      Integer  order, LDA, nRow, nCol
      Integer  idface, dot_with_normal

      Real*8   A(LDA, *)
C ======================================================================
c Local variables
      Integer i, i1, i2, i3, i4, iGauss, iL
   
      Real*8  U(9,MaxSize2D,MaxPnt2DGauss), V(9,MaxSize2D,MaxPnt2DGauss)
      Real*8  Diff(9, 9, MaxPnt2DGauss),   DU(9,MaxSize2D,MaxPnt2DGauss)

      Real*8  wG(MaxPnt2DGauss), XYG(3, MaxPnt2DGauss)
      Real*8  XYL2D(3, AllPnt2DGauss), XYL(4, AllPnt2DGauss)

      Real*8  XYP(3, 4), normal(3), area
      Integer iref(5)

      EXTERNAL tri_area, calNorm
      Real*8   tri_area, calNorm

C ======================================================================
      DATA XYL2D/3 * U1A, 
c ... 3 points (order 2)
     &         U2A,U2B,U2B,  U2B,U2A,U2B,  U2B,U2B,U2A,
c ... 7 points (order 5)
     &         U5A,U5A,U5A,  U5B,U5C,U5C,  U5C,U5B,U5C,  U5C,U5C,U5B,
     &                       U5D,U5E,U5E,  U5E,U5D,U5E,  U5E,U5E,U5D,
c ... 12 points (order 6)
     &         U6A,U6B,U6B,  U6B,U6A,U6B,  U6B,U6B,U6A,  
     &         U6C,U6D,U6D,  U6D,U6C,U6D,  U6D,U6D,U6C,  
     &         U6E,U6F,U6G,  U6G,U6E,U6F,  U6F,U6G,U6E,  
     &         U6E,U6G,U6F,  U6G,U6F,U6E,  U6F,U6E,U6G/

      DATA iref /1,2,3,4,1/
C ================================================================
      If(order.LE.0 .OR. order.EQ.4 .OR. order.GT.6) 
     &   Call errMesFEM(2001, 'fem3Dface', 
     &        'input data are incorrect: order = 1,2,3,5, or 6')

c ... weights and points
      i1 = idface
      i2 = iref(i1 + 1)
      i3 = iref(i2 + 1)
      i4 = iref(i3 + 1)

      Do i = 1, 3
        XYP(i, 1) = XY1(i)
        XYP(i, 2) = XY2(i)
        XYP(i, 3) = XY3(i)
        XYP(i, 4) = XY4(i)
      End do
  
      area = tri_area(XYP(1, i1), XYP(1, i2), XYP(1, i3)) 
      Call WeightsPnt3D(XYP(1, i1), XYP(1, i2), XYP(1, i3), area, order,
     &                  XYG, wG, iGauss, iL)

      Do i = 1, iGauss
         XYL(i1, i) = XYL2D(1, iL + i - 1)
         XYL(i2, i) = XYL2D(2, iL + i - 1)
         XYL(i3, i) = XYL2D(3, iL + i - 1)
         XYL(i4, i) = 0D0
      End do

c ... calculate optional exterior normal
      If (dot_with_normal.NE.0) Then
         Call tri_normal_ext(XYP(1, i1), XYP(1, i2),
     &                       XYP(1, i3), XYP(1, i4), normal)

         area = calNorm(normal)
         Do i = 1, 3
            normal(i) = normal(i) / area
         End do
      End if

c ... integration over volume
      Call FEM3Dinternal(
     &     XY1, XY2, XY3, XY4,
     &     operatorA, FEMtypeA, operatorB, FEMtypeB, 
     &     label, D, DDATA, IDATA, iSYS,
     &     iGauss, XYG, wG, XYL, dot_with_normal, normal,
     &     LDA, A, nRow, nCol)

      Return
      End



c ======================================================================
      Subroutine WeightsPnt3D(XY1, XY2, XY3, area, order,
     &                        XYG, w, iGauss, iLref)
c ======================================================================
      implicit none
      include 'fem3Dtri.fd'
c ======================================================================
C The procedure is used for effective computing points for numerical
C integration, b/c of a symmetry.
c ======================================================================
      Real*8  XY1(3), XY2(3), XY3(3), area
      Integer i1, i2, i3, i4, order

      Real*8  XYG(3, *), w(*)
      Integer iGauss, iLref, i

c ======================================================================
      If(order.EQ.1) Then
         iGauss = KDG1
         iLref = 1
         w(1) = S1A * area

         Do i = 1, 3
            XYG(i, 1) = U1A * (XY1(i) + XY2(i) + XY3(i))
         End do

      Else If(order.EQ.2) Then
         iGauss = KDG2
         iLref = KDG1 + 1
         Do i = 1, iGauss
            w(i) = S2A * area
         End do

         Do i = 1, 3
            XYG(i, 1) = U2B * (XY2(i) + XY3(i)) 
            XYG(i, 2) = U2B * (XY1(i) + XY3(i)) 
            XYG(i, 3) = U2B * (XY1(i) + XY2(i)) 
         End do

      Else If(order.LE.5) Then
         iGauss = KDG5
         iLref = KDG1 + KDG2 + 1

         w(1) = S5A * area
         Do i = 1, 3
            XYG(i, 1) = U5A * (XY1(i) + XY2(i) + XY3(i))
         End do

         Do i = 2, 4
            w(i) = S5B * area
         End do
         Do i = 1, 3
            XYG(i, 2) = U5B * XY1(i) + U5C * (XY2(i) + XY3(i)) 
            XYG(i, 3) = U5B * XY2(i) + U5C * (XY1(i) + XY3(i)) 
            XYG(i, 4) = U5B * XY3(i) + U5C * (XY1(i) + XY2(i)) 
         End do

         Do i = 5, 7
            w(i) = S5C * area
         End do
         Do i = 1, 3
            XYG(i, 5) = U5D * XY1(i) + U5E * (XY2(i) + XY3(i)) 
            XYG(i, 6) = U5D * XY2(i) + U5E * (XY1(i) + XY3(i)) 
            XYG(i, 7) = U5D * XY3(i) + U5E * (XY1(i) + XY2(i)) 
         End do

      Else If(order.LE.6) Then
         iGauss = KDG6
         iLref = KDG1 + KDG2 + KDG5 + 1

         Do i = 1, 3
            w(i) = S6A * area
         End do
         Do i = 1, 3
            XYG(i, 1) = U6A * XY1(i) + U6B * (XY2(i) + XY3(i)) 
            XYG(i, 2) = U6A * XY2(i) + U6B * (XY1(i) + XY3(i)) 
            XYG(i, 3) = U6A * XY3(i) + U6B * (XY1(i) + XY2(i)) 
         End do

         Do i = 4, 6
            w(i) = S6B * area
         End do
         Do i = 1, 3
            XYG(i, 4) = U6C * XY1(i) + U6D * (XY2(i) + XY3(i)) 
            XYG(i, 5) = U6C * XY2(i) + U6D * (XY1(i) + XY3(i)) 
            XYG(i, 6) = U6C * XY3(i) + U6D * (XY1(i) + XY2(i)) 
         End do

         Do i = 7, 12
            w(i) = S6C * area
         End do
         Do i = 1, 3
            XYG(i, 7) = U6E * XY1(i) + U6F * XY2(i) + U6G * XY3(i)
            XYG(i, 8) = U6G * XY1(i) + U6E * XY2(i) + U6F * XY3(i)
            XYG(i, 9) = U6F * XY1(i) + U6G * XY2(i) + U6E * XY3(i)

            XYG(i,10) = U6E * XY1(i) + U6G * XY2(i) + U6F * XY3(i)
            XYG(i,11) = U6G * XY1(i) + U6F * XY2(i) + U6E * XY3(i)
            XYG(i,12) = U6F * XY1(i) + U6E * XY2(i) + U6G * XY3(i)
         End do
      End if

      Return
      End

