C ======================================================================
      Subroutine fem3Dsub(XY1, XY2, XY3, XY4,
     &                    P1, P2, P3, P4,  Q1, Q2, Q3, Q4,
     &                    operatorA, FEMtypeA, operatorB, FEMtypeB, 
     &                    label, D, DDATA, IDATA, iSYS, order,  
     &                    LDA, A, nRow, nCol)
C ======================================================================
      implicit none
      include 'fem3Dtet.fd'
C ======================================================================
C  *** ITRODUCTION ***
C  The routine computes elemental matrix for the bilinear form 
C
C  (1)            <OpA(u), OpB(v)>               
C
C  where OpA and OpB are linear operators, and u and v are finete
C  element functions, u in "A", and v in "B". The integration is made
C  over a sub-tetrahedron XY1,XY2,XY3,XY4 which is a part 
C  of tetrahedra P1,P2,P3,P4 and Q1,Q2,Q3,Q4.
C
C  Remarks.
C    1. Update this routine together with the family
C 
C ======================================================================
      Real*8   XY1(*), XY2(*), XY3(*), XY4(*)
      Real*8   P1(*), P2(*), P3(*), P4(*),  Q1(*), Q2(*), Q3(*), Q4(*)

      Integer  FEMtypeA, FEMtypeB, operatorA, operatorB
      Integer  label, order, LDA, D, nRow, nCol

      Real*8   DDATA(*)
      Integer  IDATA(*), iSYS(*)
      EXTERNAL D

      Real*8   A(LDA, *)

C ======================================================================
c Local variables
      Real*8  XYP(3, 4), XYQ(3, 4)
      Real*8  detp, detq, cdetp, cdetq, vol, s
   
      Real*8  PSI(3, 3), QSI(3, 3)
      Real*8  U(9, MaxSize, MaxPointGauss), V(9, MaxSize, MaxPointGauss)
      Real*8  Diff(9, 9, MaxPointGauss), DU(9, MaxSize, MaxPointGauss)

      Real*8  w(MaxPointGauss), XYG(3, MaxPointGauss)
      Real*8  XYL(4, AllPointGauss)
      Real*8  XPL(4, MaxPointGauss), XQL(4, MaxPointGauss)

      Integer i,j,k,n, nfa,nfb, idim,jdim, iD,jD, tensor, iGauss, iL
      Logical ifXtensor
C ======================================================================
      DATA XYL/4 * T1A, 
c ... 4 points (order 2)
     &         T2B,T2A,T2A,T2A,  T2A,T2B,T2A,T2A,
     &         T2A,T2A,T2B,T2A,  T2A,T2A,T2A,T2B,
c ... 8 points (order 3)
     &     T3B,T3A,T3A,T3A,  T3A,T3B,T3A,T3A,  T3A,T3A,T3B,T3A,  
     &     T3A,T3A,T3A,T3B,  T3D,T3C,T3C,T3C,  T3C,T3D,T3C,T3C,
     &         T3C,T3C,T3D,T3C,  T3C,T3C,T3C,T3D,
c ... 15 points (order 5)
     &         4 * T5A,
     &     T5C,T5B,T5B,T5B,  T5B,T5C,T5B,T5B,  T5B,T5B,T5C,T5B,  
     &     T5B,T5B,T5B,T5C,  T5E,T5D,T5D,T5D,  T5D,T5E,T5D,T5D,
     &     T5D,T5D,T5E,T5D,  T5D,T5D,T5D,T5E,  T5F,T5F,T5G,T5G,  
     &     T5G,T5G,T5F,T5F,  T5F,T5G,T5G,T5F,  T5G,T5F,T5F,T5G,
     &     T5F,T5G,T5F,T5G,  T5G,T5F,T5G,T5F,
c ... 24 point (order 6)
     &     T6A,T6B,T6B,T6B,  T6B,T6A,T6B,T6B,  T6B,T6B,T6A,T6B,
     &     T6B,T6B,T6B,T6A,  T6C,T6D,T6D,T6D,  T6D,T6C,T6D,T6D,
     &     T6D,T6D,T6C,T6D,  T6D,T6D,T6D,T6C,  T6E,T6F,T6F,T6F,
     &     T6F,T6E,T6F,T6F,  T6F,T6F,T6E,T6F,  T6F,T6F,T6F,T6E,
     &     T6G,T6H,T6I,T6I,  T6G,T6I,T6H,T6I,  T6G,T6I,T6I,T6H,
     &     T6H,T6G,T6I,T6I,  T6I,T6G,T6H,T6I,  T6I,T6G,T6I,T6H,
     &     T6H,T6I,T6G,T6I,  T6I,T6H,T6G,T6I,  T6I,T6I,T6G,T6H,
     &     T6H,T6I,T6I,T6G,  T6I,T6H,T6I,T6G,  T6I,T6I,T6H,T6G/
C ================================================================

      If(order.LE.0 .OR. order.EQ.4 .OR. order.GT.6) 
     &   Call errMesFEM(2001, 'fem3Dtet', 
     &        'input data are incorrect: order = 1,2,3,5, or 6')

c ... transformation of variables y = PSI * (x - x_0)
      Do i = 1, 3
         XYP(i, 1) = 0D0
         XYP(i, 2) = XY2(i) - XY1(i)
         XYP(i, 3) = XY3(i) - XY1(i)
         XYP(i, 4) = XY4(i) - XY1(i)
      End do


      Call solve3x3(XYP(1, 2), XYP(2, 2), XYP(3, 2), PSI(1, 1),
     &              XYP(1, 3), XYP(2, 3), XYP(3, 3), PSI(1, 2),
     &              XYP(1, 4), XYP(2, 4), XYP(3, 4), PSI(1, 3), cdetp)

      Call solve3x3(XYP(1, 3), XYP(2, 3), XYP(3, 3), PSI(2, 1),
     &              XYP(1, 2), XYP(2, 2), XYP(3, 2), PSI(2, 2),
     &              XYP(1, 4), XYP(2, 4), XYP(3, 4), PSI(2, 3), detp)

      Call solve3x3(XYP(1, 4), XYP(2, 4), XYP(3, 4), PSI(3, 1), 
     &              XYP(1, 3), XYP(2, 3), XYP(3, 3), PSI(3, 2),
     &              XYP(1, 2), XYP(2, 2), XYP(3, 2), PSI(3, 3), detp)

c ... transformation of variables y = QSI * (x - x_0)
      Do i = 1, 3
         XYQ(i, 1) = 0D0
         XYQ(i, 2) = XY2(i) - XY1(i)
         XYQ(i, 3) = XY3(i) - XY1(i)
         XYQ(i, 4) = XY4(i) - XY1(i)
      End do

      Call solve3x3(XYQ(1, 2), XYQ(2, 2), XYQ(3, 2), QSI(1, 1),
     &              XYQ(1, 3), XYQ(2, 3), XYQ(3, 3), QSI(1, 2),
     &              XYQ(1, 4), XYQ(2, 4), XYQ(3, 4), QSI(1, 3), cdetq)

      Call solve3x3(XYQ(1, 3), XYQ(2, 3), XYQ(3, 3), QSI(2, 1),
     &              XYQ(1, 2), XYQ(2, 2), XYQ(3, 2), QSI(2, 2),
     &              XYQ(1, 4), XYQ(2, 4), XYQ(3, 4), QSI(2, 3), detq)

      Call solve3x3(XYQ(1, 4), XYQ(2, 4), XYQ(3, 4), QSI(3, 1), 
     &              XYQ(1, 3), XYQ(2, 3), XYQ(3, 3), QSI(3, 2),
     &              XYQ(1, 2), XYQ(2, 2), XYQ(3, 2), QSI(3, 3), detq)


c ... weights and points
      vol = dabs(detp) / 6

      Call WeightsPoints(XY1, XY2, XY3, XY4, vol, order, XYL,
     &                   XYG, w, iGauss, iL)

c ... calculate barycentric coordinates for parents
      Call barycentric3D (XY1, XY2, XY3, XY4, 
     &                    P1, P2, P3, P4,
     &                    iGauss, XYL(1,iL), XPL)
      Call barycentric3D (XY1, XY2, XY3, XY4, 
     &                    Q1, Q2, Q3, Q4,
     &                    iGauss, XYL(1,iL), XQL)

c ... compute operatorA * FEMtypeA
      If(operatorA.EQ.GRAD) Then
         Call applyGRAD(iGauss, XPL, PSI, FEMtypeA, nfa, idim, U)

      Else If(operatorA.EQ.DIV) Then
         Call applyDIV( iGauss, XPL, PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, detp)

      Else If(operatorA.EQ.IDEN) Then
         Call applyIDEN(iGauss, XPL, PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, detp)

      Else If(operatorA.EQ.CURL) Then
         Call applyCURL(iGauss, XPL, PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, detp)

      Else If(operatorA.EQ.DUDX) Then
         Call applyDUDX(iGauss, XPL, PSI, FEMtypeA, nfa, idim, U)

      Else
         Call errMesFEM(2001, 'fem3Dsub', 'operatorA is not supported')
      End if

      If(nfa.GT.LDA) Call errMesFEM(2001, 'fem3Dsub',
     &     'the local matrix leading dimension, LDA, is too small')

c ... compute operatorB * FEMtypeB
      If(operatorB.EQ.GRAD) Then
         Call applyGRAD(iGauss, XQL, QSI, FEMtypeB, nfb, jdim, V)

      Else If(operatorB.EQ.DIV) Then
         Call applyDIV( iGauss, XQL, QSI, FEMtypeB, 
     &                  nfb, jdim, V, XYQ, detq)

      Else If(operatorB.EQ.IDEN) Then
         Call applyIDEN(iGauss, XQL, QSI, FEMtypeB, 
     &                  nfb, jdim, V, XYQ,  detq)

      Else If(operatorB.EQ.CURL) Then
         Call applyCURL(iGauss, XQL, QSI, FEMtypeB, 
     &                  nfb, jdim, V, XYQ, detq)

      Else If(operatorB.EQ.DUDX) Then
         Call applyDUDX(iGauss, XQL, QSI, FEMtypeB, nfb, jdim, V)

      Else
         Call errMesFEM(2001, 'fem3Dsub', 'operatorB is not supported')
      End if

      If(nfb.GT.LDA) Call errMesFEM(2001, 'fem3Dsub',
     &     'the local matrix second dimension, LDA, is too small')
 

c ... compute D * U
      iD = idim
      jD = jdim

      Do n = 1, iGauss
         tensor = D(XYG(1,n), XYG(2,n), XYG(3,n), label, 
     &              DDATA, IDATA, iSYS, Diff(1, 1, n))
         If(ifXtensor(tensor, TENSOR_NULL)) Goto 200
      End do
      tensor = TENSOR_GENERAL

 200  Continue
      If(ifXtensor(tensor, TENSOR_NULL)) Then
         If(idim.NE.jdim) Call errMesFEM(2001, 
     &        'fem3Dsub', 'Operators A and B are not compatible')

         Do n = 1, iGauss
            Do i = 1, idim
               Do k = 1, nfa
                  DU(i, k, n) = U(i, k, n) * w(n)
               End do
            End do
         End do
      Else If(ifXtensor(tensor, TENSOR_SCALAR)) Then
         If(idim.NE.jdim) Call errMesFEM(2001, 'fem3Dsub', 
     &        'Operators A and B are not compatible')

         Do n = 1, iGauss
            s = Diff(1, 1, n) * w(n) 
            Do i = 1, idim
               Do k = 1, nfa
                  DU(i, k, n) = U(i, k, n) * s
               End do
            End do
         End do
      Else If(ifXtensor(tensor, TENSOR_SYMMETRIC) .OR. 
     &        ifXtensor(tensor, TENSOR_GENERAL)) Then
         iD = iSYS(1)
         jD = iSYS(2)

         If(jD.NE.idim .OR. iD.NE.jdim) Call errMesFEM(2001, 
     &        'fem3Dsub', 'the operators A and B are not compatible')

         Do n = 1, iGauss
            Do i = 1, iD
               Do k = 1, nfa
                  s = 0D0
                  Do j = 1, jD
                     s = s + Diff(i, j, n) * U(j, k, n)
                  End do
                  DU(i, k, n) = s * w(n)
               End do
            End do
         End do
      Else
         Call errMesFEM(2001, 'fem3Dsub', 
     &        'the given tensor is not supported') 
      End if


c ... compute <D U, V>
      Do i = 1, nfa
         Do j = 1, nfb
            s = 0D0
            Do k = 1, iD
               Do n = 1, iGauss
                  s = s + DU(k, i, n) * V(k, j, n)
               End do
            End do
            A(i, j) = s
         End do
      End do

      nRow = nfa
      nCol = nfb

      Return
      End



C ======================================================================
      Subroutine barycentric3D (XY1, XY2, XY3, XY4, 
     &                          P1, P2, P3, P4,
     &                          iGauss, XYL, XPL)
C ======================================================================
C Routine transforms local (subtetrahedron XY1,XY2,XY3,XY4) barycentric 
C coordinates XYL to global (tetrahedron P1,P2,P3,P4) barycentric 
C coordinates XPL
C ======================================================================
      implicit none

      External calVol
      Real*8   calVol

      Integer  iGauss
      Real*8   XY1(3), XY2(3), XY3(3), XY4(3)
      Real*8   P1(3), P2(3), P3(3), P4(3)
      Real*8   XYL(4, *), XPL(4, *)
      Real*8   ptet(3,4), vol


      Integer  i, j, n
      Real*8   XYsub(3, 4), L(4, 4), v, s
C ======================================================================

      Do i = 1, 3
         ptet(i, 1) = XY1(i)
         ptet(i, 2) = XY2(i)
         ptet(i, 3) = XY3(i)
         ptet(i, 4) = XY4(i)
      End do

      vol = calVol(P1, P2, P3, P4)

      Do i = 1, 4

         v = calVol(ptet(1, i), P2, P3, P4) 
         L(1, i) = dabs(v / vol)

         v = calVol(P1, ptet(1, i), P3, P4) 
         L(2, i) = dabs(v / vol)

         v = calVol(P1, P2, ptet(1, i), P4) 
         L(3, i) = dabs(v / vol)

         v = calVol(P1, P2, P3, ptet(1, i)) 
         L(4, i) = dabs(v / vol)

      End do

      Do n = 1, iGauss
         Do i = 1, 4
            s = 0D0
            Do j = 1, 4
               s = s + L(i, j) * XYL(j, n)
            End do
            XPL(i, n) = s
         End do
      End do

      Return
      End

C ======================================================================
      Logical Function ifXtensor(clr, iXtensor)
C ======================================================================
      Integer clr, iXtensor

      ifXtensor = IAND(clr, iXtensor) .EQ. iXtensor

      Return
      End



