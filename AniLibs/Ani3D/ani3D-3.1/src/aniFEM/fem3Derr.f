C ======================================================================
      Subroutine FEM3Derr(
C ======================================================================
     &      XY1, XY2, XY3, XY4, Lp,
     &      operatorA, FEMtypeA, Uh, Fu, DDATAFU, IDATAFU,
     &      label, D, DDATAFEM, IDATAFEM, iSYS, order, ERR)

C ======================================================================
      implicit none
      include 'fem3Dtet.fd'
C ======================================================================
C  *** ITRODUCTION ***
C  The routine computes the weighted norm of error:
C
C        Integral (D (OpA(u) - Fu) (OpA(u) - Fu))^{Lp/2}
C
C  where D is a tensor, OpA is a linear operator, u is the discrete
C  solution and Fu is the target (analytic or exact) function.
C
C  Details of input the parameters are described in file fem3Dtet.f.
C  Here we focus on additional parameters, Lp, Uh and Fu, only.
C
C  Lp    - norm for which the metric is to be adjusted:
C             Lp > 0  means  L_p     norm
C             Lp = 0  means  maximum norm (L_infinity)

C  Uh - Real*8 vector of discrete unknowns. It must have the same
C       size as the matrix of bilinear form <OpA(u), OpA(u)>
C
C  Fu - Integer function Futh the standard format:
C       Fu(x, y, label, DDATAFU, IDATAFU, iSYS, Diff)
C
C       The function returns type of the tendor Diff which is the
C       value of this function. See fem3Dtet.f for more details.
C
C ======================================================================
      Real*8   XY1(*), XY2(*), XY3(*), XY4(*), Lp
      Integer  operatorA, FEMtypeA, label, order, D, Fu

      Real*8   Uh(*), DDATAFEM(*), DDATAFU(*), ERR
      Integer  IDATAFEM(*), IDATAFU(*), iSYS(*)

      EXTERNAL D, Fu

C ======================================================================
c Local variables
      Real*8  XYP(3, 4)
      Real*8  det, cdet, vol, s, p
   
      Real*8  PSI(3, 3)
      Real*8  U(9, MaxSize, MaxPointGauss), V(9, MaxSize, MaxPointGauss)
      Real*8  Diff(9, 9, MaxPointGauss)
      Real*8  UUh(9, MaxPointGauss),     DUUh(9,   MaxPointGauss)

      Real*8  w(MaxPointGauss), XYG(3, MaxPointGauss)
      Real*8  XYL(4, AllPointGauss)


      Integer i,j,k,n, iD,jD, nfa, idim, iL, iGauss, tensor
C ======================================================================
      DATA XYL/4 * T1A, 
c ... 4 points (order 2)
     &     T2B,T2A,T2A,T2A,  T2A,T2B,T2A,T2A,
     &     T2A,T2A,T2B,T2A,  T2A,T2A,T2A,T2B,
c ... 8 points (order 3)
     &     T3B,T3A,T3A,T3A,  T3A,T3B,T3A,T3A,  T3A,T3A,T3B,T3A,  
     &     T3A,T3A,T3A,T3B,  T3D,T3C,T3C,T3C,  T3C,T3D,T3C,T3C,
     &     T3C,T3C,T3D,T3C,  T3C,T3C,T3C,T3D,
c ... 15 points (order 5)
     &     4 * T5A,
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
     &   Call errMesFEM(2001, 'fem3Derr', 
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
     &              XYP(1, 4), XYP(2, 4), XYP(3, 4), PSI(1, 3), cdet)

      Call solve3x3(XYP(1, 3), XYP(2, 3), XYP(3, 3), PSI(2, 1),
     &              XYP(1, 2), XYP(2, 2), XYP(3, 2), PSI(2, 2),
     &              XYP(1, 4), XYP(2, 4), XYP(3, 4), PSI(2, 3), det)

      Call solve3x3(XYP(1, 4), XYP(2, 4), XYP(3, 4), PSI(3, 1), 
     &              XYP(1, 3), XYP(2, 3), XYP(3, 3), PSI(3, 2),
     &              XYP(1, 2), XYP(2, 2), XYP(3, 2), PSI(3, 3), det)

c ... weights and points
      vol = dabs(det) / 6
      Call WeightsPoints(XY1, XY2, XY3, XY4, vol, order, XYL,
     &                   XYG, w, iGauss, iL)

c ... compute operatorA * FEMtypeA
      If(operatorA.EQ.GRAD) Then
         Call applyGRAD(iGauss, XYL(1, iL), PSI, FEMtypeA, nfa, idim, U)
         If(nfa.EQ.0) 
     &   Call applyGRA2(iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, cdet)

      Else If(operatorA.EQ.DUDX) Then
         Call applyDUDX(iGauss, XYL(1, iL), PSI, FEMtypeA, nfa, idim, U)

      Else If(operatorA.EQ.DIV) Then
         Call applyDIV( iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, cdet)

      Else If(operatorA.EQ.IDEN) Then
         Call applyIDEN(iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, cdet)

      Else If(operatorA.EQ.CURL) Then
         Call applyCURL(iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, cdet)

      Else
         Call errMesFEM(2001, 'fem3Derr', 'operatorA is not supported')
      End if

c ... compute U * Uh and save in UUh
      Do n = 1, iGauss
         Do i = 1, idim
            UUh(i, n) = 0D0
            Do k = 1, nfa
               UUh(i, n) = UUh(i, n) + U(i, k, n) * Uh(k)
            End do
         End do
      End do


c ... compute values of Fu at quadrature points
      Do n = 1, iGauss
         k = Fu(XYG(1,n),XYG(2,n),XYG(3,n), label,
     &          DDATAFU,IDATAFU,iSYS, V(1,1,n))

         If(iSYS(1).ne.idim) Call errMesFEM(2001, 'fem3Derr',
     &      'Analytic function must return a column vector, iSYS(2)=1')

         Do i = 1, idim
            UUh(i, n) = UUh(i, n) - V(i, 1, n)
         End do
      End do

c ... compute values of D * UUh and save it in DUUh
      Do n = 1, iGauss
         tensor = D(XYG(1,n), XYG(2,n), XYG(3,n), label,
     &              DDATAFEM, IDATAFEM, iSYS, Diff(1, 1, n))
         If(tensor.EQ.TENSOR_NULL) Goto 200
      End do

 200  Continue
      If(tensor.EQ.TENSOR_NULL) Then
         Do n = 1, iGauss
            Do i = 1, idim
               DUUh(i, n) = UUh(i, n)
            End do
         End do
      Else If(tensor.EQ.TENSOR_SCALAR) Then
         Do n = 1, iGauss
            Do i = 1, idim
               DUUh(i, n) = UUh(i, n) * Diff(1,1,n)
            End do
         End do
      Else If(tensor.EQ.TENSOR_SYMMETRIC .OR.
     &        tensor.EQ.TENSOR_GENERAL) Then
         iD = iSYS(1)
         jD = iSYS(2)

         Do n = 1, iGauss
            Do i = 1, iD
               s = 0D0
               Do j = 1, jD
                  s = s + Diff(i, j, n) * UUh(j, n)
               End do
               DUUh(i, n) = s
            End do
         End do
      Else
         Call errMesFEM(2001, 'femr3Derr',
     &        'the given tensor is not supported')
      End if

c ... compute |D (OpA(u) - Fu) (OpA(u) - Fu)|^p
      p = Lp / 2

      ERR = 0D0
      Do n = 1, iGauss
         s = 0D0
         Do k = 1, idim
            s = s + DUUh(k, n) * UUh(k, n)
         End do

         If(p.EQ.0D0) Then
            ERR = max(ERR, s)
         Else
            ERR = ERR + (s**p) * w(n)
         End if
      End do

      Return
      End



