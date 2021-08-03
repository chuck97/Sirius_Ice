c =====================================================================
      Program UnsteadyConvectionDiffusion
c =====================================================================
c The finite element solution of the boundary value problem:
c
c  du/dt - div (D grad u) + v * grad u = 0  in  Omega
c                                    u = g  on  dOmega
c
c where Omega = [0;1]x[0;1]x[0;1]
c           D = diag{0.0001, 0.0001, 0.0001} 
c           v = (1, 0, 0)^T
c
c SUPG stabilization is used in this example.
c =====================================================================
      implicit none

c nvmax - maximum number of mesh nodes
c ntmax - maximum number of mesh tetrahedra
c nbmax - maximum number of boundary faces
c namax - maximum number of non-zero matrix entries
      Integer   nvmax,ntmax,nbmax,namax
      parameter(nvmax = 200 000, ntmax = 6*nvmax, nbmax = 50 000)
      parameter(namax = 5 000 000)

c work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 5 000 000, MaxWi = 20 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)


c ======================================================================
c Mesh definition
c ======================================================================
      Integer  nv, nvfix, labelV(nvmax), fixedV(1)
      Real*8   vrt(3,nvmax)

      Integer  nb, nbfix, bnd(3,nbmax), labelB(nbmax), fixedB(1)
      Integer  nt, ntfix, tet(4,ntmax), labelT(ntmax), fixedT(1)

      DATA     nvfix/0/,  nbfix/0/, ntfix/0/


c =====================================================================
c for library aniRCB
c =====================================================================
c array  keeps mappings of each coarse tetrahedron to make it equilateral
      Real*8 MapMtr(3,3,ntmax)

c array keeps references of each fine cell to MapMtr
      Integer Ref2MapMtr(ntmax)


c =====================================================================
c for library aniFEM
c =====================================================================
      include 'fem3Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nvmax), JA(namax)
      Real*8    A(namax), F(nvmax), G(nvmax)

      Integer  IB(nvmax), JB(namax)
      Real*8    B(namax), U(nvmax), Uprev(nvmax)

      Integer  iSYS(MAXiSYS)
      Real*8   DDATA(5+ntmax), DATAMASS(1)
      Integer  IDATA(1), IDATANULL(1), iPrint

      Integer  Dbc 
      Real*8   calEdge
      EXTERNAL Dbc, FEM3DextStif, FEM3DextMass, calEdge


c =====================================================================
c for library aniILU
c =====================================================================
      Integer  imatvec(1), iprevec(1), luinfo, iter, nunit
      Real*8   RESID

      EXTERNAL matvec, prevec0


c LOCAL VARIABLEs
      Integer   i, i1,i2,i3,i4, n, ibc, iTime, status, nRow,nCol, iERR
      Integer   nz, NewWi, NewWr, iiEnd, irEnd, iiLU
      Real*8    x, y, z, eBC(1), Time, DeltaT

c =====================================================================
c === load a mesh
      Call loadMani(
     &     nvmax, nbmax, ntmax,
     &     nv, nb, nt,
     &     vrt, bnd, tet, labelB, labelT,
     &     nvfix, nbfix, ntfix, fixedV, fixedB, fixedT,
     &     iW, iW, "../data/cube.ani")


c === initialize the refinement (filling MapMtr, Ref2MapMtr)
      Call initializeRefinement(nv, nt, vrt, tet, MapMtr, Ref2MapMtr)

c refine the mesh uniformly 
      Do i = 1, 4
         Call uniformRefinement(
     &        nv, nvmax, nb, nbmax, nt, ntmax,
     &        vrt, bnd, tet, labelB, labelT,
     &        MapMtr, Ref2MapMtr,
     &        iW, MaxWi)
      End do
      Write(*,'(A,I6,A)') 'The refined mesh has ',nt,' tetrahedra'


c time integration step
      DeltaT = 0.015D0

c data provided for the user subroutines Ddiff Dconv used in FEM3DextStif
      DDATA(1) = 1d-4 ! D  
      DDATA(2) = 1    ! v_x
      DDATA(3) = 0    ! v_y
      DDATA(4) = 0    ! v_y
      DDATA(5) = DeltaT  

      Do i = 1, nt
         DDATA(5+i) = 0D0
      End do

      IDATA(1) = 5 + nt ! length of DDATA


c === assemble the stifness matrix
c mark the Dirichlet points on a side of a unit cube
      Do n = 1, nv
         labelV(n) = 0

         x = vrt(1, n)
         y = vrt(2, n)
         z = vrt(3, n)

         If(x.EQ.0D0) labelV(n) = 1

         If(x.EQ.0D0 .AND. 
     &      y.GE.0.25D0 .AND. y.LE.0.75D0 .AND.
     &      z.GE.0.25D0 .AND. z.LE.0.75D0) Then
            labelV(n) = 2
         End if
      End do

c general sparse matrix in the AMG format
      status = IOR(MATRIX_GENERAL, FORMAT_AMG)
      iPrint = 1

      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelV, bnd, labelB, tet, labelT,
     &     FEM3DextStif, DDATA, IDATA, status,
     &     nvmax, namax, IA, JA, A, F, nRow, nCol,
     &     MaxWi, iW, iPrint)


c === call the driver for ILU solution
      nz = IA(nRow + 1) - 1

      irEnd = nz + 1
      NewWr = MaxWr - irEnd

      If(NewWr.LT.0) Then
         Write(*,*) 'Real*8 memory is not sufficient'
         Stop
      End if

      iiLU  = nRow + 2
      iiEnd = iiLU + nz
      NewWi = MaxWi - iiEnd

      If(NewWi.LT.0) Then
         Write(*,*) 'Integer memory is not sufficient'
         Stop
      End if

      Call ilu0(nRow, A,JA,IA, rW, iW(iiLU), iW, iW(iiEnd), iERR)

      If(iERR.NE.0) Then
         Write(*,*) 'Initialization of ilu0 has failed, iERR =', iERR
         Stop
      End if
c Keep data in rW(1:irEnd-1) and iW(1:iiEnd-1) 


c === make the mass matrix scaled by 1.5/DeltaT
c general sparce matrix in the AMG format
      status = IOR(MATRIX_GENERAL, FORMAT_AMG)
      iPrint = 1

      DATAMASS(1) = 1.5d0 / DeltaT

c generate bilinear form using build-in function ANI_D1x1_scalar
      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelV, bnd, labelB, tet, labelT,
     &     FEM3DextMass, DATAMASS, IDATANULL, status,
     &     nvmax, namax, IB, JB, B, U, nRow, nCol,
     &     NewWi, iW(iiEnd), iPrint)


c === set the initial state U = 0
      Call dcopy(nv, 0D0, 0, U, 1)

      Do i = 1, nv
         if (labelV(i).GT.0) Then
            x = vrt(1, i)
            y = vrt(2, i)
            z = vrt(3, i)

            ibc = Dbc(x, y, z, labelV(i), DDATA, IDATA, iSYS, eBC)

            U(i)= eBC(1) 
         End if
      End do
         
c copy U to Uprev
      Call dcopy(nv, U, 1, Uprev, 1)


c ======================================================================
c Time iteration loop
c ======================================================================
      Time  = 0
      iTime = 0

 2    Continue
         iTime = iTime + 1
         Time  = Time  + DeltaT

         write(*,'(A,I2,A,F5.2)') 'Time step ', iTime, ' time=', Time

c === regenerate RHS due to SUPG
         Do i = 1, nt
            i1 = tet(1, i)
            i2 = tet(2, i)
            i3 = tet(3, i)
            i4 = tet(4, i)
            DDATA(5+i) = 2*(U(i1) + U(i2) + U(i3) + U(i4)) / 4 - 
     &       (Uprev(i1) + Uprev(i2) + Uprev(i3) + Uprev(i4)) / 8
         End do

         status = IOR(MATRIX_GENERAL, FORMAT_AMG)
         iPrint = 1

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tet, labelT,
     &        FEM3DextStif, DDATA, IDATA, status,
     &        nvmax, namax, IA, JA, A, F, nRow, nCol,
     &        NewWi, iW(iiEnd), iPrint)


c === BFD time step: (1.5 u_i - 2 u_{i-1} + 0.5 u_{i-2}) / DeltaT
c note: mass matrix is scaled with 1.5/DeltaT
         Call daxpy(nRow, -4d0, U,1, Uprev,1)
         Call dscal(nRow, -1d0/3d0, Uprev,1)

         Call mulAgen(nRow, IB, JB, B, Uprev, G)

         Call daxpy(nRow, 1d0, F,1, G,1)
         Call dcopy(nRow, U,1, Uprev,1)


c solve Ax=b, without iterative refinement
         iter   = 1000      ! max number of iterations
         RESID  = 1D-12     ! threshold for \|RESID\|
         luinfo = 0         ! no troubles on input
         nunit  = 6         ! output channel
         iprevec(1) = nRow  ! single entry required: system size 
         imatvec(1) = nRow  ! single entry required: system size 
     
         Call slpbcgs(prevec0, iprevec, iW, rW,
     &                matvec, imatvec, IA,JA,A,
     &                rW(irEnd), nRow, 8,
     &                nRow, G, U,
     &                iter, RESID,
     &                luinfo, nunit)

c assign the boundary condition
         Do i = 1, nv
            If(labelV(i).GT.0) then
               x =  vrt(1, i)
               y =  vrt(2, i)
               z =  vrt(3, i)

               ibc  = Dbc(x, y, z, labelV(i), DDATA, IDATA, iSYS, eBC)
               U(i) = eBC(1) 
            End if
         End do


c === save the initial solution and the mesh
         If(iTime.EQ.1) 
     &      Call GMVscalarVrt(U, "solution_initial.gmv", 10,
     &                        nv, vrt, nt, tet, nb, bnd, labelB)


      If(Time.LT.0.5) Goto 2
c END OF TIME INTEGRATION


c === save the final solution and the mesh
      Call GMVscalarVrt(U, "solution_final.gmv", 10,
     &                  nv, vrt, nt, tet, nb, bnd, labelB)

      Stop 
      End





