      Program  main
c ==========================================================
c The program demonstrates the use of local refinement/coarsening library librcb3D.a based
c on the Marked Edge Bisection.  It starts from a simple coarse tetrahedral meshes,
c (the mesh may be arbitrary conformal)  
c refines it nlevel times according to the rule given in RefineRule,
c and coarse the refined mesh 4 times according to the rule given in CoarseRule.
c Example of RefineRule, CoarseRule are at the file end.
c Prior the refinement/coarsening the routine  InitializeRCB has to be called (once).
c The rules are designed so that the mesh is refined towards a point source moving along
c the edge of the computational polyhedral domain.
c
c Caution! Curve-linear boundary faces are processed as straight:
c          no crv data are on input/output.
c The program uses  routines [GMVmesh] from libview3D.a and [error] from libmba3D.a.
c
c The output  is  three gmv-pictures of the initial mesh (ini.gmv) and the mesh refined towards
c the initial position of the source (fst.gmv) and the final position of the source (lst.gmv)
c ==========================================================
      implicit none

c ... user defined procedures 
      external  RefineRule,  CoarseRule 

c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary faces
      Integer   nvmax,ntmax,nbmax
      Parameter (nvmax = 200 000, ntmax = 6*nvmax, nbmax = 100 000)

c ... standard mesh arrays 
c ... number of points, tetrahedra and boundary faces
      Integer   nv, nt, nb

c ... coordinates of mesh points 
      Double precision   vrt(3,nvmax)

c ... connectivity table for tetrahedra, and their label
      Integer   tet(4,ntmax),material(ntmax)

c ... connectivity table for boundary faces, and their labels
      Integer   bnd(3,nbmax),labelF(nbmax)

c ... additional dummy arrays  needed to call loadMaft() from libmba3D.a 
      Integer  nvfix, ivfix(1), nbfix, ibfix(1), ntfix, itfix(1)

c ... maxlevel - maximum number of levels for refinement and coarsening
      Integer   maxlevel
      Parameter (maxlevel = 20)

c ... work memory (at least 15*ntmax+41)
      Integer   MaxWi
      Parameter (MaxWi =  20 000 000)
      Integer   iW(MaxWi), WiW

      Integer   MaxWr
      Parameter (MaxWr =   7 000 000) 
      Integer   rW(MaxWr)

c ... history of bisection
      Logical history(4*maxlevel*ntmax)

c ... list of tetrahedra
      Integer   tetpmax
      Parameter (tetpmax = 290)
      Integer   listtet((tetpmax+1)*ntmax)

c ... number of levels of refinement/coarsening
      integer  nlevel


c ... user data to be passed to RefineRule 
      Double precision RefineRuleData(1)
c ... user data to be passed to CoarseRule 
      Double precision CoarseRuleData(1)

c ... local variables
      Integer   ilevel, iERR, iPrint
      character*20 fname

c ==========================================================

c Step 1: load the initial mesh

      Call loadMaft( nvmax, nbmax, ntmax,  nv, nb, nt,
     &     vrt, bnd, tet, labelF, material,
     &     nvfix, nbfix, ntfix, ivfix, ibfix, itfix,
     &     iW, iW, "../src/Tutorials/PackageRCB/mesh_ini.out")


      Write(*,'(A,3I7)')
     &     '  Initial mesh: numbers of nodes, tets and boundary faces:',
     &      nv, nt,nb

      fname = 'ini.gmv'
      Call GMVmesh(fname, 99,nv,vrt, nt,tet, nb,bnd,labelF)

c Step 2: initialize data structure (work memory is at least 15*ntmax+41)

      iERR     = 0
      Call InitializeRCB (nt, ntmax, nv, nvmax, vrt, tet, 
     &      MaxWi, iW, MaxWr, rW, listtet, tetpmax, iERR)
      If(iERR.GT.0) Call errMes(iERR, 'main', 'error in InitializeRCB')
      iPrint   = 1       ! low level of output information

      nlevel = 6


c Step 3: refine initial mesh nlevel times by local bisection. 
c         The rule is defined in RefineRule
c         RefineRuleData is dummy array in this version of routine RefineRule
      
      Do ilevel = 1, nlevel
        Call LocalRefine (
     &          nv, nvmax, nb, nbmax, nt, ntmax,
     &          vrt, tet, bnd, material,labelF,
     &          RefineRule, ilevel,
     &          maxlevel, history, 
     &          listtet, tetpmax, RefineRuleData,
     &          MaxWi, iW, 
     &          iPrint, iERR)
        If(iERR.GT.0) Call errMes(iERR, 'main',
     &                 'error in LocalRefine')
      End do

      Write(*,'(A,3I7)')
     &     '  Refined mesh: numbers of nodes, tets and boundary faces:',
     &      nv, nt, nb
c
      fname = 'rfn.gmv'
      call GMVmesh(fname, 99,nv,vrt, nt,tet, nb,bnd,labelF)
c
c Step 4: coarse the refined mesh 3 times by local coarsening. 
c         The rule is defined in CoarseRule
c         CoarseRuleData is dummy array in this version of routine CoarseRule

      Do ilevel = nlevel, nlevel-2, -1
         Call LocalCoarse (
     &        nv, nvmax, nb, nbmax, nt, ntmax,
     &        vrt, tet, bnd, material,labelF,
     &        CoarseRule, ilevel,
     &        maxlevel, history,
     &        listtet, tetpmax, CoarseRuleData,
     &        MaxWi, iW,
     &        iPrint, iERR)
         If(iERR.GT.0) Call errMes(iERR, 'main',
     &                     'error in LocalCoarse')
      End do

c      
      Write(*,'(A,3I7)')
     &    'Coarsened mesh: numbers of nodes, tets and boundary faces:',
     &      nv, nt,nb
c
      fname = 'crs.gmv'
      call GMVmesh(fname, 99,nv,vrt, nt,tet, nb,bnd,labelF)

      Stop 
      End

C ================================================================
c User routine  defines the rule for local refinement depending on 
c current level.
c On input: nt  current number of elements
c           tet current connectivity table
c           vrt current coordinates
c           ilevel current level of refinement
c           RefineRuleData  data array for refining
c On output: verf marker for refinement of each element
c            (0 - no need to refine, 1 - refine by single bisection,
c             2 - refine by two   levels of bisection,
c             3 - refine by three levels of bisection preserving the shape)
C ================================================================
      Subroutine RefineRule (nt, tet, vrt, verf, ilevel, RefineRuleData)
C ================================================================
      implicit none
      
      Integer                nt
      Integer                tet(4,*)
      Double precision       vrt(3,*)
      Integer                verf(*)
      Integer                ilevel
      Double precision       RefineRuleData(*) ! dummy in this case

      Integer                i, k
C ================================================================
c ... refine all cells
      Do i = 1, nt
          verf(i) =  1! one level of bisection
      End do 


      Return
      End

C ================================================================
c User routine  defines the rule for local coarsening depending on 
c current level.
c On input: nt  current number of elements
c           tet current connectivity table
c           vrt current coordinates
c           ilevel current level of refinement
c           CoarseRuleData  data array for coarsening
c On output: verf marker for coarsening of each element
c            (0 - no need to coarse, 1 - coarse by single merging,
c             2 - coarse by two   levels of merging,
c             3 - coarse by three levels of merging   preserving the shape)
C ================================================================
      Subroutine CoarseRule (nt, tet, vrt, verf, ilevel, CoarseRuleData)
C ================================================================
      implicit none
      
      Integer                nt
      Integer                tet(4,*)
      Double precision       vrt(3,*)
      Integer                verf(*)
      Integer                ilevel
      Double precision       CoarseRuleData(*) ! dummy in this case

      Integer                i
      Double precision       xy, xy1, xy2, xy3
C ================================================================
c Uniform coarsening with moderate rate of coarsening
c
c Here we exploit the feature of coarsening procedure:
c coarsening of a coarse tetrahedron causes multiple coarsenings of
c finer tetrahedra due to mesh conformity
      Do i = 1, nt
          verf(i) =  1 ! one level  of merging
      End do


      Return
      End


