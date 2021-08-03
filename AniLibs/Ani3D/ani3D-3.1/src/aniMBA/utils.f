C ================================================================
      Subroutine listP2P(nP, nE, MaxList, IPE, nPP, IPP, iW)
C ================================================================
C  The routine creates connectivity lists P->P for mesh points.
C
C  *** Remarks:
C         1. iW(*) - working memory of size 4 * nE + nP
C ================================================================
      Integer IPE(4, *), nPP(*), IPP(*)
      Integer iW(*)

C ================================================================
      inEP = 0
      iIEP = inEP + nP

      Call backReferences(nP, nE, 4,4, IPE, iW(inEP + 1), iW(iIEP + 1))
 
c ... main algorithm: array nEP is overloaded inside
      nL = 0

      i2 = 0
      Do n = 1, nP
         nLo = nL

         i1 = i2 + 1
         i2 = iW(inEP + n)

         Do m = i1, i2
            iE = iW(iIEP + m)

            Do j = 1, 4
               iPt = IPE(j, iE)

               If(iW(inEP + iPt).GT.0) Then
                  nL = nL + 1
                  If(nL.GT.MaxList) Call errMes(2011, 'listP2P',
     &                             'user parameter MaxList is small')

                  IPP(nL) = iPt

                  iW(inEP + iPt) = -iW(inEP + iPt)
               End if
            End do
         End do

         nPP(n) = nL

c  ...  recovering values of array nEP
         Do m = nLo + 1, nL
            iPt = IPP(m)
            iW(inEP + iPt) = -iW(inEP + iPt)
         End do
      End do

      Return
      End


 
C ================================================================
      Subroutine listR2P(nP, nR, nE, MaxR, IPE, IPR, iW)
C ================================================================
C  The routine reates connectivity list R -> P for mesh edges. 
C  The algorithm has a linear arithmetical complexity. 
C
C  Parameters:
C      nP, nR, nE - the number of points, edges and elements
C      MaxR       - the maximal number of edges
C
C      IPE(4, nE) - connectivity list of tetrahedra: E -> P
C      IPR(2, nR) - connectivity list of edges:      R -> P
C
C      iW(*) - working memory of size 4 * nE + nP
C ================================================================
      Integer IPE(4, *), IPR(2, *)
      Integer iW(*)

C ================================================================
      inEP = 0
      iIEP = inEP + nP

      Call backReferences(nP, nE, 4, 4, IPE, iW(inEP + 1), iW(iIEP + 1))


c ... main algorithm: array nEP is overloaded inside
      nR = 0
   
      i2 = 0
      Do n = 1, nP
         nRo = nR

         i1 = i2 + 1
         i2 = iW(inEP + n)

         Do m = i1, i2
            iE = iW(iIEP + m)

            Do j = 1, 4
               iPt = IPE(j, iE)

               If(iPt.GT.n .AND. iW(inEP + iPt).GT.0) Then
                  nR = nR + 1
                  If(nR.GT.MaxR) Then
                     Call errMes(2011, 'listR2P', 
     &                          'user parameter MaxR is small')
                  End if

                  IPR(1, nR) = n
                  IPR(2, nR) = iPt

                  iW(inEP + iPt) = -iW(inEP + iPt)
               End if
            End do
         End do

c  ...  recovering values of array nEP
         Do m = nRo + 1, nR
            iPt = IPR(2, m)
            iW(inEP + iPt) = -iW(inEP + iPt)
         End do
      End do

      Return
      End



C ================================================================
      Subroutine listR2R(nP, nR, nE, MaxL, IPE, nRR, IRR, iW, iERR)
C ================================================================
C  The routine creates connectivity lists R->R for mesh edges.
C  Routine returns 0 upon successful completion.
C
C  *** Remarks:
C         1. iW(*) - working memory of size 12 * nE + nR
C ================================================================
      Integer IPE(4, *), nRR(*), IRR(*)
      Integer iW(*)

C ================================================================
      iERR = 0

      iIRE = 1
      inEP = iIRE + 6 * nE
      iIEP = inEP + nP
      iEnd = iIEP + 4 * nE

      Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))

      inER = inEP
      iIER = inER + nR
      iEnd = iIER + 6 * nE
      Call backReferences(nR, nE, 6,6, iW(iIRE), iW(inER), iW(iIER))

      nL = 0

      i2 = 0
      Do n = 1, nR
         nLo = nL

         i1 = i2 + 1
         i2 = iW(inER + n - 1)

         Do m = i1, i2
            iE = iW(iIER + m - 1)

            Do 100 j = 1, 6
               iRt = iW(iIRE + 6 * (iE - 1) + j - 1)

               Do k = nLo + 1, nL
                  If(iRt.EQ.IRR(k)) Goto 100
               End do

               nL = nL + 1
               If(nL.GT.MaxL) Then
                  iERR = n
                  Goto 9000
               End if

               IRR(nL) = iRt
 100        Continue
         End do

         nRR(n) = nL
      End do

 9000 Return
      End



C ================================================================
      Subroutine listE2R(nP, nR, nE, IPE, IRE, nEP, IEP)
C ================================================================
C  The routine computes connectivity lists for mesh edges
C
C  nEP(*) - working array of size nP
C  IEP(*) - working array of size 4*nE
C ================================================================
      Integer IPE(4, *), IRE(6, *)
      Integer IEP(*), nEP(*)
C ================================================================
C group (Local variables)
      Integer ipr(2, 6)
      Logical check22

      DATA ipr /1,2, 1,3, 1,4, 2,3, 2,4, 3,4/

C ================================================================
      Call backReferences(nP, nE, 4, 4, IPE, nEP, IEP)

      Do n = 1, nE
         Do i = 1, 6
            IRE(i, n) = 0
         End do
      End do

      nR = 0
      Do n = 1, nE
         Do 20 i = 1, 6
            If(IRE(i, n).EQ.0) Then
               nR = nR + 1
               IRE(i, n) = nR

               ip1 = IPE(ipr(1, i), n)
               ip2 = IPE(ipr(2, i), n)
               ip3 = max(ip1, ip2)

               m1 = 1
               If(ip3.GT.1) m1 = nEP(ip3 - 1) + 1

               Do m = m1, nEP(ip3)
                  iE = IEP(m)
                  Do j = 1, 6
                     jp1 = IPE(ipr(1, j), iE)
                     jp2 = IPE(ipr(2, j), iE)

                     If(check22(ip1, ip2, jp1, jp2)) Then
                        If(IRE(j, iE).EQ.0) IRE(j, iE) = nR
                     End if
                  End do
               End do
            End if
 20      Continue
      End do

      Return
      End



C ======================================================================
      Subroutine listE2F(nP, nF, nE, IPE, IFE, nEP, IEP)
C ======================================================================
C  The routine creates connectivity list E->F for mesh elements
C
C  *** Remarks:
C         1. Working memory is nEP(nP), IEP(4 * nE)
C ======================================================================
      Integer IPE(4, *), IFE(4, *)
      Integer IEP(*), nEP(*)

C ======================================================================
C group (Local variables)
      Integer iref(5)
      Logical cmpE, check33

      DATA    iref /1,2,3,4,1/

C ======================================================================
      Call backReferences(nP, nE, 4, 4, IPE, nEP, IEP)

      Do n = 1, nE
         Do i = 1, 4
            IFE(i, n) = 0
         End do
      End do

      nF = 0
      Do n = 1, nE
         Do 10 i1 = 1, 4
            If(IFE(i1, n).EQ.0) Then
               nF = nF + 1
               IFE(i1, n) = nF

               i2 = iref(i1 + 1)
               i3 = iref(i2 + 1)

               iP1 = IPE(i1, n)
               iP2 = IPE(i2, n)
               iP3 = IPE(i3, n)

               If(cmpE(iP1, iP2, iP3, IEP, nEP, n, iE2)) Then
                  Do j1 = 1, 4
                     j2 = iref(j1 + 1)
                     j3 = iref(j2 + 1)

                     jP1 = IPE(j1, iE2)
                     jP2 = IPE(j2, iE2)
                     jP3 = IPE(j3, iE2)
                     If(check33(iP1, iP2, iP3, jP1, jP2, jP3)) Then
                        IFE(j1, iE2) = nF
                        Goto 10
                     End if
                  End do
               End if
            End if
 10      Continue
      End do

      Return
      End



C ================================================================
      Subroutine listE2Fb(nP, nFb, nE, IPF, IPE, IFE, nEP, IEP)
C ================================================================
C  The routine computes connectivity list E->F for BOUNDARY faces
C
C  *** Remarks:
C         1. Working memory is nEP(nP), IEP(3 * nFb)
C ================================================================
      Integer IPF(3, *), IPE(4, *), IFE(4, *)
      Integer IEP(*), nEP(*)

C ================================================================
C group (Local variables)
      Integer iref(5)
      Logical cmpE

      DATA    iref /1,2,3,4,1/

C ================================================================
      Call backReferences(nP, nFb, 3, 3, IPF, nEP, IEP)

      Do n = 1, nE
         Do i1 = 1, 4
            IFE(i1, n) = 0

            i2 = iref(i1 + 1)
            i3 = iref(i2 + 1)

            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)
            iP3 = IPE(i3, n)

            If(cmpE(iP1, iP2, iP3, IEP, nEP, 0, iF)) Then
               IFE(i1, n) = iF
            End if
         End do
      End do

      Return
      End

C ======================================================================
      Subroutine listE2E(nP, nE, IPE, IEE, nEP, IEP)
C ======================================================================
C  The routine computes connectivity lists E->E for neighboring
C  teterahedra.
C ======================================================================
      Implicit none
      Integer nP, nE
      Integer IPE(4, *), IEE(4, *), IEP(*), nEP(*)
C group (Local variables)
      Integer j(3), ib(3), ie(3)
      Integer n, i, i1, i2, i3, j1, j2, j3, k, m
C ======================================================================
      Call backReferences(nP, nE, 4, 4, IPE, nEP, IEP)

      Do n = 1, nE
        Do i = 1, 4

          IEE(i, n) = 0

          m = 0
          Do k = 1, 4
            If (k .NE. i) Then
              m = m + 1
              j(m) = IPE(k, n)
            End if
          End do

          Do k = 1, 3
            If(j(k) .EQ. 1) Then
              ib(k) = 1
             Else
              ib(k) = nEP(j(k) - 1) + 1
            End if
            ie(k) = nEP(j(k))
          End do

          Do i1 = ib(1), ie(1)
            j1 = IEP(i1)
            If (j1 .EQ. n) Go to 30

            Do i2 = ib(2), ie(2)
              j2 = IEP(i2)
              If (j2 .EQ. n) Go to 20

              Do i3 = ib(3), ie(3)
                j3 = IEP(i3)
                If (j3 .EQ. n) Go to 10

                If (j1 .EQ. j2 .and. j1 .EQ. j3) Then
                   IEE(i, n) = j1
                   Goto 40
                End if
10              Continue
              End do

20            Continue
            End do

30          Continue
          End do

40        Continue
        End do
      End do


      Return
      End




C ================================================================
      Subroutine listConv(
     &           nP, nR, nE, nEP, IEP, L, IRE, 
     &           nX, MaxX, nRP, IRP, iW, iERR)
C ================================================================
C  The routine convolutes unstructured map X->Y and structured map 
C  Y->Z to get the map X->Z. For examples, if X means points (P), 
C  Y means elements (E), and Z means edges (R), we get the map from 
C  a point to all edges in the elements having this point.  
C
C  Routine returns 0 upon successful completion.
C
C  *** Remarks:
C         1. Working memory is iW(nR)
C ================================================================
      Integer nP, nR, nE, L, MaxX
      Integer nEP(*), IEP(*), IRE(L, *), nRP(*), IRP(*), iW(*)
C ================================================================
      iERR = 0

      Do n = 1, nR
         iW(n) = 1
      End do

      nX = 0
      i2 = 0

      Do n = 1, nP
         nX0 = nX

         i1 = i2 + 1
         i2 = nEP(n)

         Do i = i1, i2
            iE = IEP(i)

            Do j = 1, L
               iR = IRE(j, iE)  
               If(iW(iR).GT.0) Then
                  nX = nX + 1
                  If(nX.GT.MaxX) Then
                     iERR = n
                     Goto 9000
                  End if

                  IRP(nX) = iR
                  iW(iR) = -iW(iR)
               End if
            End do
         End do 

         nRP(n) = nX
        
c ...    restore array iW
         Do i = nX0 + 1, nX
            iR = IRP(i)
            iW(iR) = -iW(iR)
         End do 
      End do      

 9000 Return
      End



C ================================================================
      Subroutine reverseMap(nP, nR, nRP, IRP, nPR, IPR)
C ================================================================
C Routine creates map R->P reverse to the map P->R.
C ================================================================
      Integer nRP(*), IRP(*), nPR(*), IPR(*)

      Do n = 1, nR
         nPR(n) = 0
      End do

      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = nRP(n)

         Do i = i1, i2
            iR = IRP(i)
            nPR(iR) = nPR(iR) + 1
         End do
      End do

      Do n = 2, nR
         nPR(n) = nPR(n) + nPR(n - 1)
      End do

      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = nRP(n)

         Do i = i1, i2
            iR = IRP(i)

            iref = nPR(iR)
            IPR(iref) = n

            nPR(iR) = iref - 1
         End do
      End do

      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = nRP(n)

         Do i = i1, i2
            iR = IRP(i)
            nPR(iR) = nPR(iR) + 1
         End do
      End do

      Return
      End



C ======================================================================
      Subroutine backReferences(nP, nE, L, M, IPE, nEP, IEP)
C ======================================================================
C Routine creates map P->E reverse to the map E->P.
C     nEP(P) - nEP(P-1) = number of elements having common
C                         point P.
C     IPE([nEP(P-1) + 1 : nEP(P)]) = list of elements having
C                         common point P.
C ======================================================================
      Integer IPE(M, *), nEP(*), IEP(*), iLast

      Do n = 1, nP 
         nEP(n) = 0
      End do

      Do n = 1, nE
         Do i = 1, L
            i1 = IPE(i, n)
            nEP(i1) = nEP(i1) + 1
         End do
      End do

      Do n = 2, nP
         nEP(n) = nEP(n) + nEP(n - 1)
      End do
      iLast = nEP(nP)

      Do n = 1, nE
         Do i = 1, L
            i1 = IPE(i, n)

            iref = nEP(i1)
            IEP(iref) = n

            nEP(i1) = iref - 1
         End do
      End do

      Do n = 1, nP - 1
         nEP(n) = nEP(n + 1)
      End do
      nEP(nP) = iLast

      Return
      End



C ======================================================================
      Subroutine backReferencesOld(nP, nE, L, M, IPE, nEP, IEP)
C ======================================================================
C Routine creates map P->E reverse to the map E->P.
C     nEP(P) - nEP(P-1) = number of elements having common
C                         point P.
C     IPE([nEP(P-1) + 1 : nEP(P)]) = list of elements having
C                         common point P.
C ======================================================================
      Integer IPE(M, *), nEP(*), IEP(*)

      Do n = 1, nP 
         nEP(n) = 0
      End do

      Do n = 1, nE
         Do i = 1, L
            i1 = IPE(i, n)
            nEP(i1) = nEP(i1) + 1
         End do
      End do

      Do n = 2, nP
         nEP(n) = nEP(n) + nEP(n - 1)
      End do

      Do n = 1, nE
         Do i = 1, L
            i1 = IPE(i, n)

            iref = nEP(i1)
            IEP(iref) = n

            nEP(i1) = iref - 1
         End do
      End do

      Do n = 1, nE
         Do i = 1, L
            i1 = IPE(i, n)
            nEP(i1) = nEP(i1) + 1
         End do
      End do

      Return
      End



C ================================================================
      Subroutine dualNormals(nP, nR, nE, IPE, IPR, XYP, NRM, iW)
C ================================================================
C  The routine computes normals to surfaces of a dual mesh
C  which separate points of a primary mesh.
C
C  Parameters:
C      nP, nR, nE - the number of points, edges and elements in the
C                   primary mesh
C
C      IPE(4, nE) - connectivity list of tetrahedra: E -> P
C      IPR(2, nR) - connectivity list of edges:      R -> P
C
C      XYP(3, nP) - Cartesian coordinates of mesh points
C      NRM(3, nR) - array of unit vectors normal to surfaces of a
C                   dual mesh separating point of the primary mesh
C
C      iW(*) - working memory of size 4 * nE + nP
C ================================================================
      Integer IPE(4, *), IPR(2, *)
      Integer iW(*)

      Real*8  XYP(3, *), NRM(3, *)

C ================================================================
C (Local variables)
      Real*8  XYA(3), XYB(3), XYE(3), XYC(3), XYD(3)
      Real*8  XYM(3), XYN(3), XYR(3)
      Real*8  s, calNorm, DotMul
 
C ================================================================
      inEP = 0
      iIEP = inEP + nP

      Call backReferences(nP, nE, 4, 4, IPE, iW(inEP + 1), iW(iIEP + 1))

      Do n = 1, nR
         iP1 = IPR(1, n)
         iP2 = IPR(2, n)

         Do i = 1, 3
            NRM(i, n) = 0D0
            XYE(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2

            XYR(i) = XYP(i, iP2) - XYP(i, iP1)
         End do

         i1 = 1
         If(iP1.GT.1) i1 = iW(inEP + iP1 - 1) + 1

         i2 = iW(inEP + iP1)

         s = 0D0
         Do 10 m = i1, i2
            iE = iW(iIEP + m)

            icnt = 0
            Do j = 1, 4
               iPt = IPE(j, iE)
               If(iPt.NE.iP1 .AND. iPt.NE.iP2) Then
                  icnt = icnt + 1
                  If(icnt.EQ.1) iP3 = iPt
                  If(icnt.EQ.2) iP4 = iPt
               End if
            End do

            If(icnt.NE.2) Goto 10

c   ...   computing vertices of tet differ from iP1 & iP2  
            Do i = 1, 3
               XYC(i) = (XYP(i, iP1) + XYP(i, iP2) 
     &                 + XYP(i, iP3) + XYP(i, iP4)) / 4

               XYA(i) = (XYP(i, iP1) + XYP(i, iP2) 
     &                               + XYP(i, iP3)) / 3 - XYC(i)

               XYB(i) = (XYP(i, iP1) + XYP(i, iP2) 
     &                               + XYP(i, iP4)) / 3 - XYC(i)

               XYD(i) = XYE(i) - XYC(i)
            End do

            Call VecMul(XYD, XYA, XYM)
            Call VecMul(XYD, XYB, XYN)

            s = s + calNorm(XYM) + calNorm(XYN)
            
c    ...    orienting the normal vector from iP1 to iP2
            If(DotMul(XYR, XYM).LT.0D0) Then
               Do i = 1, 3
                  NRM(i, n) = NRM(i, n) - XYM(i)
               End do
            Else
               Do i = 1, 3
                  NRM(i, n) = NRM(i, n) + XYM(i)
               End do
            End if

            If(DotMul(XYR, XYN).LT.0D0) Then
               Do i = 1, 3
                  NRM(i, n) = NRM(i, n) - XYN(i)
               End do
            Else
               Do i = 1, 3
                  NRM(i, n) = NRM(i, n) + XYN(i)
               End do
            End if
 10      Continue

c  ...  rescaling the area of the interface separating iP1 & iP2
c        s = s / calNorm(NRM(1, n))
c        Do i = 1, 3
c           NRM(i, n) = NRM(i, n) * s
c        End do
      End do

      Return
      End


