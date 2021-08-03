c ==========================================================
      Program  main
c ==========================================================
c This program loads a  mesh (data/cube_hole.out) 
c split it into two local submeshes with frozen interface, updates
c each submesh, and glues the updated submeshes back to global mesh
c ==========================================================
      implicit none

      integer nvmax,ntmax,nbmax
c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh tetrahedra
c ... nbmax - maximum number of boundary faces
c ... MaxP = nvmax, MaxE = 6 * MaxP, MaxF = nbmax, MaxP, MaxE, MaxF
      parameter(nvmax = 10000, ntmax = 6*nvmax, nbmax = 10000)


c Available memory. NOTE: rW(maxWr) is needed ONLY for call mbaFixShape 
c                         rW(maxWr) is not needed for splitting/uniting submeshes
      Integer    MaxWi, MaxWr
      Parameter( MaxWi =  25 000 000, MaxWr = 1 000 000) 
      Integer iW(maxWi)
      Real*8  rW(maxWr)

c Global mesh
c ... standard mesh arrays (see doc/aft_guide.pdf for more detail)
c ... number of points, tets, and boundary faces
      Integer  nv, nt, nb

c ... coordinates of mesh points 
      Real*8   vrt(3,nvmax)
c ... connectivity table for tets and their labels
      Integer  tet(4,ntmax),mat(ntmax)

c ... connectivity table for boundary faces, and their labels
      Integer  bnd(3,nbmax),lab(nbmax)

c     ivfix(nvfix) is array of fixed (never touched) points
c     ibfix(nbfix) is array of fixed (never touched) faces
c     itfix(ntfix) is array of fixed (never touched) elements
      Integer   nvfix, ivfix(1), nbfix, ibfix(1), ntfix, itfix(1) 

c ... colors of tets used for definition of submeshes
      Integer TetColors(ntmax)

c Two local submeshes
      Integer  nvL(2), ntL(2), nbL(2)      
      Real*8   vrtL(3,nvmax,2)
      Integer  tetL(4,ntmax,2),matL(ntmax,2) 
      Integer  bndL(3,nbmax,2),labL(nbmax,2)
c     Integer  nvfixL(2), ivfixL(1,2), nbfixL(2), ibfixL(1,2), 
      Integer  nvfixL(2), ivfixL(1,2), nbfixL(2), ibfixL(nbmax,2), 
     &         ntfixL(2), itfixL(1,2) 

c     IPPl(nPl) is references to a global enumeration necessary to glue submeshes in one global mesh
c     IFFl(nFl) is references to global enumeration of surface triangles
c     nFvi      is the number of interface triangles
      Integer  nFviL(2), IPPL(nvmax,2), IFFL(nbmax,2)

c Work arrays keeping concatinated data of local meshes
      Real*8   vrtall(3,nvmax*2)
      Integer  tetall(4,ntmax*2),matall(ntmax*2) 
      Integer  bndall(3,nbmax*2),laball(nbmax*2)
c     Integer  ivfixall(1*2), ibfixall(1*2), itfixall(1*2) 
      Integer  ivfixall(1*2), ibfixall(nbmax*2), itfixall(1*2) 
      Integer  IPPall(nvmax*2), IFFall(nbmax*2)

c ... library MBA, in order to call mbaFixShape
c Basket capacity and the maximum number of local modifications
      Integer   MaxSkipE, MaxQItr
      Parameter(MaxSkipE = 300, MaxQItr = 500 000)
c Desired final mesh quality
      Real*8    Quality
      Parameter(Quality = 4D-1)
      Real*8    rQuality
      Logical   flagAuto
      Integer   iPrint, status, iERR
c The routine which defines the nodal metric to be Euclidean
      Integer   ANI_Metric_Eucl
      External  ANI_Metric_Eucl


c Local variables
      Integer   i,j,k,m,iClr,ip1,ip2,ip3,ip4,iEnd
      Real*8    xc
c ==========================================================

c ... Load the mesh in aft-format
      Call loadMaft(
     &     nvmax, nbmax, ntmax,
     &     nv, nb, nt,
     &     vrt, bnd, tet, lab, mat,
     &     nvfix, nbfix, ntfix, ivfix, ibfix, itfix,
     &     iW, iW, "../data/cube_hole.out")


      Write(*,*)
      Write(*,'(A,3(I8,A))')'The loaded mesh has ',
     &                      nt,' tetrahedra, ',
     &                      nv,' nodes, ',
     &                      nb,' boundary faces '

c ... Check the input mesh
      ip1 = 1
      ip2 = ip1 + nv
      ip3 = ip2 + nt*4
      ip4 = ip3 + nt*4
      iEnd = ip4+ nt*4
      if (iEnd.gt.MaxWi) stop 'increase MaxWi'

      call check_mesh(
     &            nv, nb, nt,
     &            vrt, bnd, tet, lab, mat,
     &            iW(ip1),iW(ip2),iW(ip3),iW(ip4) )

      Write(*,'(A)')'The loaded mesh is checked'

c ... Split tetrahedra into 2 subsets according to position with respect to plane x=0.5
      Do i = 1, nt
         xc = 0.25*(vrt(1,tet(1,i))+vrt(1,tet(2,i))
     &             +vrt(1,tet(3,i))+vrt(1,tet(4,i)))
         If (xc.le.0.5d0) Then
             TetColors(i) = 1
         Else
             TetColors(i) = 2
         End If
      End do

c ... extract submeshes 
      Do iClr = 1, 2

        ip1 = 1
        ip2 = ip1 + 2*nv
        ip3 = ip2 + nb
        iEnd= ip3 + 2*nv + 3*nb + 4*nt
        if (iEnd.gt.MaxWi) stop 'increase MaxWi'

        call global2local(
     &           iClr, TetColors,      ! extract submesh with TetColors=iClr
c group (Mg)
     &           nv, nb, nbmax, nt,
     &           vrt, bnd, tet, lab, mat,
     &           nvfix, nbfix, ntfix, ivfix, ibfix, itfix, 
c group (Ml)
     &           nvL(iClr), nbL(iClr), ntL(iClr),
     &           vrtL(1,1,iClr), bndL(1,1,iClr), tetL(1,1,iClr), 
     &           labL(1,iClr), matL(1,iClr),
     &           nvfixL(iClr),   nbfixL(iClr), ntfixL(iClr), 
     &           ivfixL(1,iClr), ibfixL(1,iClr), itfixL(1,iClr), 
c group (I)
     &           nFviL(iClr), IPPL(1,iClr), IFFL(1,iClr), 
c group (W)
     &           iW(ip1),iW(ip2),iW(ip3) ) 

        Write(*,*)
        Write(*,'(A,4(I8,A))')'The extracted mesh has ',
     &                      ntL(iClr),' tetrahedra, ',
     &                      nvL(iClr),' nodes, ',
     &                      nbL(iClr),' boundary faces, ',
     &                      nFviL(iClr),' interfaces'

c ... Check  submesh 
        ip1 = 1
        ip2 = ip1 + nvL(iClr)
        ip3 = ip2 + ntL(iClr)*4
        ip4 = ip3 + ntL(iClr)*4
        iEnd = ip4+ ntL(iClr)*4
        if (iEnd.gt.MaxWi) stop 'increase MaxWi'

        call check_mesh(
     &            nvL(iClr), nbL(iClr), ntL(iClr),
     &            vrtL(1,1,iClr), bndL(1,1,iClr), tetL(1,1,iClr),
     &            labL(1,iClr), matL(1,iClr),
     &            iW(ip1),iW(ip2),iW(ip3),iW(ip4) )
        Write(*,'(A,I2,A)')'Submesh ',iClr,' is checked'
        Write(*,*)

c ... Modify local meshes

c Important restriction of mbaFixShape: labL and matL MUST be POSITIVE !!!
        Do i = 1, ntL(iClr)
           If (matL(i,iClr).le.0) stop 'non-positive tet material'
        End do
        Do i = 1, nbL(iClr)
           If (labL(i,iClr).le.0) stop 'non-positive bnd face label'
        End do
c freeze all boundary faces in order to preserve volume of the domain
        nbfixL(iClr) = nbL(iClr)
        Do i = 1, nbfixL(iClr)
         ibfixL(i,iClr) = i
        End do
c mesh fixed points (not to be touched in adaptation) may be detected automatically since
c they separate boundary edges with different colors
        nvfixL(iClr) = 0
        ntfixL(iClr) = 0
c make up the local mesh: generate a shape regular mesh with ntL tetrahedra.
c Shape regularity is understood in the metric defined in subroutine
c ANI_Metric_Eucl (see below).
        flagAuto = .TRUE.  ! default mesh generation options
        status = 1         ! forbid boundary elements
        iPrint = 1         ! average level of output information

        Call mbaFixShape(
c group (M)
     &     nvL(iClr), nvmax, nbL(iClr), nbmax, ntL(iClr), ntmax,
     &     vrtL(1,1,iClr), bndL(1,1,iClr), tetL(1,1,iClr), 
     &     labL(1,iClr), matL(1,iClr),
c group (Dev)
     &     nvfixL(iClr), nbfixL(iClr), ntfixL(iClr), 
     &     ivfixL(1,iClr), ibfixL(1,iClr), itfixL(1,iClr),
     &     flagAuto, status,
c group (Q)
     &     MaxSkipE, MaxQItr,
     &     ANI_Metric_Eucl, Quality, rQuality,
c group (W)
     &     MaxWr, MaxWi, rW, iW,
     &     iPrint, iERR)
        Write(*,'(A,I2,A)')'Submesh ',iClr,' is modified'

      End Do


 
c Create a GMV-file with the mesh. The name must have extension .gmv
c     Call saveMgmv(nvL(2), nbL(2), ntL(2), 
c    &         vrtL(1,1,2),bndL(1,1,2),tetL(1,1,2),labL(1,2),matL(1,2), 
c    &         "save.gmv", iW)


c ... Unite submeshes back to global mesh

c prepare concatinated arrays of local meshes
      j = 0
      k = 0
      m = 0
      Do iClr = 1, 2
       Do i = 1, nbL(iClr)
          j = j + 1
          bndall(1,j) = bndL(1,i,iClr) 
          bndall(2,j) = bndL(2,i,iClr) 
          bndall(3,j) = bndL(3,i,iClr) 
          laball(  j) = labL(  i,iClr) 
          IFFall(  j) = IFFL(  i,iClr)
       End Do
       Do i = 1, nvL(iClr)
          k = k + 1
          vrtall(1,k) = vrtL(1,i,iClr) 
          vrtall(2,k) = vrtL(2,i,iClr) 
          vrtall(3,k) = vrtL(3,i,iClr) 
          IPPall(  k) = IPPL(  i,iClr)
       End Do
       Do i = 1, ntL(iClr)
          m = m + 1
          tetall(1,m) = tetL(1,i,iClr) 
          tetall(2,m) = tetL(2,i,iClr) 
          tetall(3,m) = tetL(3,i,iClr) 
          tetall(4,m) = tetL(4,i,iClr) 
          matall(  m) = matL(  i,iClr) 
       End Do
      End Do

      j = 0
      k = 0
      m = 0
      Do iClr = 1, 2
       Do i = 1, nbfixL(iClr)
          j = j + 1
          ibfixall(j) = ibfixL(i,iClr) 
       End Do
       Do i = 1, nvfixL(iClr)
          k = k + 1
          ivfixall(k) = ivfixL(i,iClr) 
       End Do
       Do i = 1, ntfixL(iClr)
          m = m + 1
          itfixall(m) = itfixL(i,iClr) 
       End Do
      End Do

c ... unite local meshes
      ip1 = 1
      ip2 = ip1 + nvmax
      ip3 = ip2 + nvmax
      ip4 = ip3 + nbmax
      iEnd = ip4+ ntmax*4
      if (iEnd.gt.MaxWi) stop 'increase MaxWi'

      Call local2global(
     &           2,    ! number of meshes to be united
c group (Mg)           ! global mesh on output
     &           nvmax, nbmax, ntmax,
     &           nv, nb, nt,
     &           vrt, bnd, tet, lab, mat,
     &           nvfix, nbfix, ntfix, ivfix, ibfix, itfix, 
c group (Ml)           ! local meshes (concatinated arrays)
     &           nvL, nbL, ntL,
     &           vrtall, bndall, tetall, laball, matall,
     &           nvfixL, nbfixL, ntfixL, 
     &           ivfixall, ibfixall, itfixall, 
c group (I)
     &           IPPall, IFFall,
c group (W)
     &           iW(ip1),iW(ip2),iW(ip3),iW(ip4) )

      Write(*,*)
      Write(*,'(A,3(I8,A))')'The united mesh has ',
     &                      nt,' tetrahedra, ',
     &                      nv,' nodes, ',
     &                      nb,' boundary faces '

c ... Check the united mesh
      ip1 = 1
      ip2 = ip1 + nv
      ip3 = ip2 + nt*4
      ip4 = ip3 + nt*4
      iEnd = ip4+ nt*4
      if (iEnd.gt.MaxWi) stop 'increase MaxWi'

      call check_mesh(
     &            nv, nb, nt,
     &            vrt, bnd, tet, lab, mat,
     &            iW(ip1),iW(ip2),iW(ip3),iW(ip4) )

      Write(*,'(A)')'The united mesh is checked'


c Create a GMV-file with the mesh. The name must have extension .gmv
      Call saveMgmv(nv, nb, nt, 
     &         vrt,bnd,tet,lab,mat, 
     &         "save.gmv", iW)


      Stop
      End

