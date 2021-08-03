C =====================================================================
      Program Main
C =====================================================================
      implicit none
 
c nvmax - maximum number of mesh nodes
c ntmax - maximum number of mesh tetrahedron
c nbmax - maximum number of boundary edges
c namax - maximum number of non-zero matrix entries
      Integer   nvmax,ntmax,nbmax,namax, nvMetaMax,ntMetaMax
      parameter(nvmax = 10 000, ntmax = 6*nvmax, nbmax = 5 000)
      parameter(namax = 2 000 000)
      parameter(nvMetaMax = 50 * nvmax, ntMetaMax = 50 * ntmax)

c work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 3 000 000, MaxWi = 5 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)

C ======================================================================
C Mesh definition (first mesh)
C ======================================================================
      Integer  nv, labelV(nvmax)
      Real*8   vrt(3,nvmax)
      Integer  nvfix, fixedV(nvmax)

      Integer  nb, bnd(3,nbmax), labelB(nbmax)
      Integer  nbfix, fixedB(nbmax)

      Integer  nt, tet(4,ntmax), labelT(ntmax)
      Integer  ntfix, fixedT(ntmax)

C ======================================================================
C Mesh definition (second mesh)
C ======================================================================
      Integer  nv2, labelV2(nvmax)
      Real*8   vrt2(3,nvmax)
      Integer  nvfix2, fixedV2(nvmax)

      Integer  nb2, bnd2(3,nbmax), labelB2(nbmax)
      Integer  nbfix2, fixedB2(nbmax)

      Integer  nt2, tet2(4,ntmax), labelT2(ntmax)
      Integer  ntfix2, fixedT2(ntmax)

C ======================================================================
C Meta mesh (intermidiate data structure)
C ======================================================================
      Integer  nv12, nt12, tet12(4,ntMetaMax), parents(2,ntMetaMax)
      Real*8   vrt12(3,nvMetaMax)

C ======================================================================
C For library aniFEM
C ======================================================================
      include 'fem3Dtet.fd'
      include 'assemble.fd'

      Integer  IA(nvmax), JA(namax), controlFEM(3), nRow, nCol
      Real*8   A(namax)
      Real*8   DDATA(1)
      Integer  IDATA(1), iPrint

      EXTERNAL FEM3Dext

C ======================================================================
C For library aniILU
C ======================================================================
      EXTERNAL matvec, prevec2
      Integer  imatvec(1), iprevec(1)

      Integer  iter, info, nunit, verb, UsedWr, UsedWi
      Real*8   resid, tau1,tau2, partlur,partlurout 

C ======================================================================
C For uniformRefinement from library aniMBA
C ======================================================================
c ... array  keeps mappings of each coarse tetrahedron to make it equilateral
      Real*8 MapMtr(3,3,ntmax)
c ... array keeps references of each fine cell to MapMtr
      Integer Ref2MapMtr(ntmax)

C ======================================================================
C For library aniPRJ
C ======================================================================
      Real*8  U1(nvmax), U2(nvmax)

C LOCAL VARIABLEs
      Integer  i, j, n, iERR, ipBCG, status
      Real*8   x, y,  z, errMax, s

C ======================================================================
c ... load the first mesh. 
      Call loadMaft(
     &      nvmax, nbmax, ntmax, nv, nb, nt,   
     &      vrt, bnd, tet, labelB, labelT,  
     &      nvfix, nbfix, ntfix, fixedV, fixedB, fixedT,
     &      iW, iW, "../data/cub6.out")
      Write(*,'(A,I6,A)')
     &      'The loaded mesh 1 has ',nt,' tetrahedra'
c ... Initialize the refinement (filling MapMtr, Ref2MapMtr)
      Call initializeRefinement(
     &      nv, nt, vrt, tet,
     &      MapMtr, Ref2MapMtr)
c ... Refine the first mesh uniformly several times
      Do i = 1, 2
         Call uniformRefinement(
     &      nv, nvmax, nb, nbmax, nt, ntmax,
     &      vrt, bnd, tet, labelB, labelT,
     &      MapMtr, Ref2MapMtr,
     &      iW, MaxWi)
      End do
      Write(*,'(A,I6,A)')
     &      'The uniformly refined mesh 1 has ',nt,' tetrahedra'


c ... load the second mesh.
      Call loadMani(
     &      nvmax, nbmax, ntmax, nv2, nb2, nt2,   
     &      vrt2, bnd2, tet2, labelB2, labelT2,  
     &      nvfix2, nbfix2, ntfix2, fixedV2, fixedB2, fixedT2,
     &      iW, iW, "../data/cube.ani")
      Write(*,'(A,I6,A)')
     &      'The loaded mesh 2 has ',nt2,' tetrahedra'
c ... Initialize the refinement (filling MapMtr, Ref2MapMtr)
      Call initializeRefinement(
     &      nv2, nt2, vrt2, tet2,
     &      MapMtr, Ref2MapMtr)
c ... Refine the first mesh uniformly several times
      Do i = 1, 1
         Call uniformRefinement(
     &      nv2, nvmax, nb2, nbmax, nt2, ntmax,
     &      vrt2, bnd2, tet2, labelB2, labelT2,
     &      MapMtr, Ref2MapMtr,
     &      iW, MaxWi)
      End do
      Write(*,'(A,I6,A)')
     &      'The uniformly refined mesh 2 has ',nt2,' tetrahedra'

c ... rescale meshes for bigger code coverage
      Do n = 1, nv
         Do i = 1, 3
            vrt(i, n) = vrt(i, n) * 2
         End do
      End do
      Do n = 1, nv2
         Do i = 1, 3
            vrt2(i, n) = vrt2(i, n) * 2
         End do
      End do

c === draw the meshes (no labels of mesh points is available)
c The name must have extension with .gmv
      Call GMVmesh("../bin/mesh1.gmv", 10, 
     &             nv, vrt, nt, tet, 
     &             nb, bnd, labelB)

      Call GMVmesh("../bin/mesh2.gmv", 10, 
     &             nv2, vrt2, nt2, tet2, 
     &             nb2, bnd2, labelB2)

c === no extra data is provided for the user subroutine Ddiff
      DDATA(1) = 0
      IDATA(1) = 0

c general sparse matrix in the AMG format 
      status = IOR(MATRIX_GENERAL, FORMAT_AMG)
      iPrint = 1

      Do i = 1, nv
         labelV(i) = 0
      End do  
     
      Do i = 1, nb
         Do j = 1, 3
            labelV(bnd(j,i)) = 1
         End do   
      End do  
      
      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelV, bnd, labelB, tet, labelT,
     &     fem3Dext, DDATA, IDATA, status,
     &     nvmax, namax, IA, JA, A, U2, nRow, nCol,
     &     MaxWi, iW, iPrint)


c === library aniPRJ: rescale meshes and create a meta-mesh
      Call MetaMesh(nv, vrt, nt, tet, 
     &              nv2, vrt2, nt2, tet2,
     &              nv12, nvMetaMax, vrt12, 
     &              nt12, ntMetaMax, tet12, parents,
     &              MaxWi, MaxWr, iW, rW, iERR)
c
c === draw the meta-mesh  (no labels of mesh points is available)
      Call GMVmesh("../bin/metamesh.gmv", 10, 
     &             nv12, vrt12, nt12, tet12, 
     &             0, iW, iW)

c assemble the right-hand side
      Do i = 1, nv2
         x = vrt2(1, i)
         y = vrt2(2, i)
         z = vrt2(3, i)
         U2(i) = x
      End do
 
      Call assemble_rhs(nv, vrt, nt, tet, nv2, vrt2, nt2, tet2,
     &                  nv12, vrt12, nt12, tet12, parents, 
     &                  IDEN, FEM_P2, IDEN, FEM_P1,
     &                  U1, U2, MaxWi, iW)
 
c === ILU SOLVER 
c initialization of the preconditioner
      verb    = 0     ! verbose no
      tau1    = 1d-2  ! absolute threshold for L,U
      tau2    = 1d-3  ! absolute threshold for T,R
      partlur = 0.5   ! even partition of memory between LU and R
      iERR    = 0     ! error code
 
      Call iluoo(nRow, IA, JA, A, tau1, tau2, verb,
     &           rW, iW, MaxWr, MaxWi, partlur, partlurout,
     &           UsedWr, UsedWi, iERR)
 
      If(iERR.NE.0) Then
         Write(*,*) 'Initialization(1) of iluoo failed, iERR=', iERR
         Stop
      End if

c iterative solution
      If(UsedWr + 8*nRow.GT.MaxWr) Then
         Write(*,'(A,I7)') 'Increase MaxWr to ', UsedWr + 8*nRow
         Stop
      End if
 
      ipBCG = UsedWr + 1
 
      iter  = 10000  ! max number of iterations
      resid = 1d-12  ! final residual
      info  = 0      ! no troubles on input
      nunit = 6      ! output to display 
      Do i = 1, nRow ! initial guess
         U2(i) = 1d0
      End do
 
      iprevec(1) = nRow
      imatvec(1) = nRow
      Call slpbcgs(prevec2, iprevec, iW, rW,
     &             matvec,  imatvec, IA, JA, A,
     &             rW(ipBCG), nRow, 8,
     &             nRow, U1, U2,
     &             iter, resid, info, nunit)

      If(info.NE.0) Stop 'BiCGStab had failed'
 
 
c === draw remapped function
      errMax = 0D0
      Do i = 1, nv
         x = vrt(1, i)
         y = vrt(2, i)
         z = vrt(3, i)
c        write (*,*) 'i=', i, 'U2(i)=',U2(i),x
         errMax = max(errMax, dabs(x - U2(i)))
      End do
      Write(*,'(A,E14.7)') 'Maximal remap error =', errMax 
 
      Call GMVscalarVrt(U2, "remap.gmv", 10,
     &                  nv, vrt, nt, tet, nb, bnd, labelB)

      Stop
      End


C ======================================================================
C  The user defined routines required above
C ======================================================================
c Templated routine for the elemental matrix. It calls standard bilinear
c forms and imposes the boundary conditions using the provided labels.
C ======================================================================
      Subroutine FEM3Dext(XY1, XY2, XY3, XY4,
     &                    lbE, lbF, lbR, lbP, DDATA, IDATA, iSYS,
     &                    LDA, A, F, nRow, nCol,
     &                    templateR, templateC)
C ======================================================================
      implicit none
      Include 'fem3Dtet.fd'
      Include 'assemble.fd'

      Real*8  XY1(*), XY2(*), XY3(*), XY4(*)
      
      Integer lbE, lbF(4), lbR(6), lbP(4)
      Real*8  DDATA(*)
      Integer IDATA(*)
      Integer iSYS(*), LDA, nRow, nCol

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C LOCAL VARIABLEs
      Integer  Ddiff, Dreact, Drhs, Dbc, DbcRobCoef
      External Ddiff, Dreact, Drhs, Dbc, DbcRobCoef

      Real*8   B(4, 4), C(3, 3), G(3), XYP(3, 4)
      Real*8   x, y, z, eBC(1)

      Integer  i,j,k,l,m, ir, ic, label, ibc
  
      Integer  iref(5), ip(4)
      DATA     iref/1,2,3,4,1/

      Integer  ANI_Dnull
      EXTERNAL ANI_Dnull
 
C ======================================================================
      nRow = 10
      nCol = 10

c ... set up templated degrees of freedom for rows and columns. 
c     used convention 'V'=vertex d.o.f. and 'R'=edge d.o.f.
      Do i = 1, 4
         templateR(i) = Vdof
      End do
      Do i = 1, 6
        templateR(i+4) = Rdof
      End do
      Do i = 1, 10     
        templateC(i) = templateR(i)
      End do

c ... compute the stiffness matrix (M)
      label = lbE

c     A(1:10,1:10) is elemental vector elliptic operator;
c     in other words, for the bilinear form <iden(P2), iden(P2)>
      Call fem3Dtet(XY1, XY2, XY3, XY4,
     &              IDEN, FEM_P2, IDEN, FEM_P2,
     &              label, ANI_Dnull, DDATA, IDATA, iSYS, 5,
     &              LDA, A, ir, ic)

      Return
      End

