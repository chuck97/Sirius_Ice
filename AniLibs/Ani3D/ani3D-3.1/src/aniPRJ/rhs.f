C ======================================================================
      Subroutine assemble_rhs(nv, vrt, nt, tet, nv2, vrt2, nt2, tet2,
     &                        nv12, vrt12, nt12, tet12, parents, 
     &                        operatorA, FEMtypeA, operatorB, FEMtypeB,
     &                        U1, U2, MaxWi, iW)
C ======================================================================
C Routine calculates the right-hand side vector U1 = M12 U2, where
C Ui is the function on the i-th mesh and M12 is the mass or stiffness
C (for interpolation in general normed spaces) matrix corresponing
C to cross-integration of basis functions v1_i and v2_j on different 
C meshes:
C
C    \Int_\Omega v1_i v2_j dx
C
C The matrix M12 is not assembled and multiplication is done using the
C meta-mesh.
C
C At the moment, the integration weigth is a constant, but it may be 
C generalized to spatially-varying weigths.
C ======================================================================
      implicit none
      include 'fem3Dtet.fd'
      include 'assemble.fd'

      Integer  nv, nt, tet(4, *), nv2, nt2, tet2(4, *)
      Real*8   vrt(3, *), vrt2(3, *)

      Integer  nv12, nt12, tet12(4, *), parents(2, *)
      Real*8   vrt12(3, *), U1(*), U2(*)

      Integer  operatorA, FEMtypeA, operatorB, FEMtypeB
      Integer  MaxWi, iW(*)


C LOCAL VARIABLES
      Integer  IDATA(1), iSYS(MAXiSYS), label, order, nRow, nCol
      Real*8   DDATA(1), M12(MaxSize, MaxSize), s, tvol

      Integer  dof(MaxSize), dof2(MaxSize), ip(MaxSize), iq(MaxSize)
      Integer  ANI_Dnull, orderDOF
      EXTERNAL ANI_Dnull, orderDOF

      Integer  i, j, k1, k2, n, n1, n2
      Integer  iv1, iv2, iv3, iv4
      Integer  ip1, ip2, ip3, ip4, iq1, iq2, iq3, iq4
      Integer  nr, nr2, nf, nf2 
      Integer  iIRE, iIRE2, iIFE, iIFE2, inEP, iIEP, iEnd
      
C ======================================================================
      label = 1   !dummy: is not used by the ANI_Dnull routine
      order = orderDOF(FEMtypeA) + orderDOF(FEMtypeB)  !maximal quadrature

c memory allocation for mapping
      iIRE  = 1
      iIRE2 = iIRE  + 6 * nt
      iIFE  = iIRE2 + 6 * nt2
      iIFE2 = iIFE  + 4 * nt
      inEP  = iIFE2 + 4 * nt2
      iIEP  = inEP  + max(nv, nv2)
      iEnd  = iIEP  + 4 * max(nt, nt2)

      If(iEnd.GT.MaxWi) Call  errMesFEM(1001, 
     &   'assemble_rhs', 'Not enough integer memory')

c maps E->R for both meshes
      Call listE2R(nv,  nr,  nt,  tet,  iW(iIRE),  iW(inEP), iW(iIEP))
      Call listE2R(nv2, nr2, nt2, tet2, iW(iIRE2), iW(inEP), iW(iIEP))

c maps E->F for both meshes
      Call listE2F(nv,  nf,  nt,  tet,  iW(iIFE),  iW(inEP), iW(iIEP))
      Call listE2F(nv2, nf2, nt2, tet2, iW(iIFE2), iW(inEP), iW(iIEP))

c create local -> global maps
      Call listDOF(FEMtypeA, k1, dof)  
      Call listDOF(FEMtypeB, k2, dof2)  


c assemble the righ-hand side
      Call dof2nRow(k1, dof, nv, nr, nf, nt, nRow)
c
      tvol = 0.0d0
      Do n = 1, nRow
         U1(n) = 0
      End do

      Do n = 1, nt12
        iv1 = tet12(1, n)
        iv2 = tet12(2, n)
        iv3 = tet12(3, n)
        iv4 = tet12(4, n)

        n1 = parents(1, n)
        n2 = parents(2, n)

c vertices of a tetrahedron on mesh #1
        ip1 = tet(1, n1)
        ip2 = tet(2, n1)
        ip3 = tet(3, n1)
        ip4 = tet(4, n1)

c vertices of a tetrahedron on mesh #2
        iq1 = tet2(1, n2)
        iq2 = tet2(2, n2)
        iq3 = tet2(3, n2)
        iq4 = tet2(4, n2)

c calling integration routine from library aniFEM
         Call fem3Dsub(vrt12(1,iv1), vrt12(1,iv2),
     &                 vrt12(1,iv3), vrt12(1,iv4),  
     &                 vrt(  1,ip1), vrt(  1,ip2),
     &                 vrt(  1,ip3), vrt(  1,ip4), 
     &                 vrt2( 1,iq1), vrt2( 1,iq2),
     &                 vrt2( 1,iq3), vrt2( 1,iq4),  
     &                 operatorA, FEMtypeA, operatorB, FEMtypeB,
     &                 label, ANI_Dnull, DDATA, IDATA, iSYS, order,
     &                 MaxSize, M12, nRow, nCol)

c      Do i = 1, nRow
c         Write(*, '(100F8.3)') (M12(i, j), j = 1, nCol)
c      End do

c integrating over tetrahedron iv1-iv2-iv3-iv4
        Call dof2global(k1, dof,  nv,  nr,  nf,  
     &                 n1, tet, iW(iIRE), iW(iIFE),  ip)
        Call dof2global(k2, dof2, nv2, nr2, nf2,
     &                  n2, tet2, iW(iIRE2), iW(iIFE2), iq)

        Do i = 1, nRow
          s = 0D0
          Do j = 1, nCol
            s = s + M12(i, j) * U2(iq(j))
            tvol = tvol + M12(i, j)
          End do
          U1(ip(i)) = U1(ip(i)) + s
c      Write(*,'(A,E14.8)') 'PRJ: DEBUG of tvol =', tvol
        End do

      End do
      

      Return
      End



C ======================================================================
      Subroutine dof2global(N, dof, nP, nR, nF, kE, IPE, IRE, IFE, ind)
C ======================================================================
C Default rules are used for local->global mapping of DOF for element iE
C ======================================================================
      implicit none
      include 'fem3Dtet.fd'

      Integer N, dof(*), nP, nR, nF
      Integer IPE(4, *), IRE(6, *), IFE(4, *), ind(*), kE

      Integer i, j, ip, ifs, ir, ie, isR, isF, isE

C ======================================================================
c cound gropus of unknowns
 
      ip  = 0
      ir  = 0
      ifs = 0
      ie  = 0

      Do i = 1, N
        If (dof(i).EQ.Vdof) Then
          ip = ip + 1
        Else If(dof(i).EQ.Rdof .OR. dof(i).EQ.RdofOrient) Then
          ir = ir + 1
        Else If(dof(i).EQ.Fdof .OR. dof(i).EQ.FdofOrient) Then
          ifs = ifs + 1
        Else If(dof(i).EQ.Edof) Then
          ie = ie + 1
        End if
      End do
      ip  = ip  / 4 
      ir  = ir  / 6 
      ifs = ifs / 4 
 
c define global shifts
      isR = ip * nP
      isF = isR + ir  * nR
      isE = isF + ifs * nF

      ip  = 0
      ir  = 0
      ifs = 0
c
      i = 1
      Do While(i.LE.N)
         If(dof(i).EQ.Vdof) Then
            Do j = 1, 4
              ind(i + j - 1)  = ip * nP + IPE(j, kE) 
            End do 

            i  = i + 4
            ip = ip + 1
c
         Else If(dof(i).EQ.Rdof .OR. dof(i).EQ.RdofOrient) Then
            Do j = 1, 6
              ind(i + j - 1)  = isR + ir * nR + IRE(j, kE)
            End do
c
            i  = i + 6
            ir = ir + 1

         Else If(dof(i).EQ.Fdof .OR. dof(i).EQ.FdofOrient) Then
            Do j = 1, 4
              ind(i + j - 1)  = isF + ifs * nF + IFE(j, kE)
            End do
c
            i   = i   + 4
            ifs = ifs + 1

         Else If(dof(i).EQ.Edof) Then
            Do j = 1, ie
               ind(i + j - 1) = isE + ie * (kE - 1) + j 
            End do

            i = i + ie
         End if
      End do

      Return
      End



C ======================================================================
      Subroutine dof2nRow(N, dof, nP, nR, nF, nE, nRow)
C ======================================================================
C Calculates the size of the array U1.
C ======================================================================
      implicit none
      include 'fem3Dtet.fd'

      Integer N, dof(*), nP, nF, nR, nE, nRow
      Integer i, ip, ifs, ir, ie
C ======================================================================
c count groups of unknowns

      ip  = 0
      ir  = 0
      ifs = 0
      ie  = 0
      Do i = 1, N
        If (dof(i).EQ.Vdof) Then
          ip = ip + 1
        Else If(dof(i).EQ.Rdof .OR. dof(i).EQ.RdofOrient) Then
          ir = ir + 1
        Else If(dof(i).EQ.Fdof .OR. dof(i).EQ.FdofOrient) Then
          ifs = ifs + 1
        Else If(dof(i).EQ.Edof) Then
          ie = ie + 1
        End if
      End do
      ip  = ip  / 4 
      ir  = ir  / 6 
      ifs = ifs / 4 

      nRow = ip * nP + ir * nR + ifs * nF + ie * nE

      Return
      End


