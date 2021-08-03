C ================================================================
      Real*8 Function calVol(xy1, xy2, xy3, xy4)
C ================================================================
C The oriented volume of the tetrahedron given by 4 vertices
C ================================================================
      Real*8 xy1(3), xy2(3), xy3(3), xy4(3)

C group (Local variables)
      Real*8 v1(3), v2(3), v3(3)

      Do i = 1, 3
         v1(i) = xy1(i) - xy4(i)
         v2(i) = xy2(i) - xy4(i)
         v3(i) = xy3(i) - xy4(i)
      End do

      calVol = (v1(1) * (v2(2) * v3(3) - v2(3) * v3(2)) +
     &          v1(2) * (v2(3) * v3(1) - v2(1) * v3(3)) +
     &          v1(3) * (v2(1) * v3(2) - v2(2) * v3(1))) /
     &          6D0
      Return
      End



c ======================================================================
      Real*8 Function tet_diameter(xy1, xy2, xy3, xy4)
c ======================================================================
c  Routines computes diameter of the tetrahedron 
c ======================================================================
      implicit none
      Real*8   xy1(3), xy2(3), xy3(3), xy4(3)
      Real*8   sqrEdge, d1, d2, d3, d4, d5, d6
c ======================================================================
      d1 = sqrEdge(xy1, xy2)
      d2 = sqrEdge(xy1, xy3)
      d1 = max(d1, d2)

      d3 = sqrEdge(xy1, xy4)
      d4 = sqrEdge(xy2, xy3)
      d3 = max(d3, d4)

      d5 = sqrEdge(xy2, xy4)
      d6 = sqrEdge(xy3, xy4)
      d5 = max(d5, d6)

      d1 = max(d1, d3)
      d1 = max(d1, d5)

      tet_diameter = dsqrt(d1)

      Return
      End

 

c ======================================================================
      Subroutine tet_center(xy1, xy2, xy3, xy4, x, y, z)
c ======================================================================
c  Routines computes center of the tetrahedron
c ======================================================================
      implicit none
      Real*8   xy1(3), xy2(3), xy3(3), xy4(3), x, y, z
c ======================================================================
      x = (xy1(1) + xy2(1) + xy3(1) + xy4(1)) / 4
      y = (xy1(2) + xy2(2) + xy3(2) + xy4(2)) / 4
      z = (xy1(3) + xy2(3) + xy3(3) + xy4(3)) / 4

      Return
      End



c ======================================================================
      Subroutine tet_heights(xy1, xy2, xy3, xy4, H)
c ======================================================================
c  Routines computes four heights of the tetrahedron
c ======================================================================
      implicit none
      Real*8   xy1(3), xy2(3), xy3(3), xy4(3), H(4)
      Real*8   calVol, tri_area, v, s

c ======================================================================
      v = calVol(xy1, xy2, xy3, xy4)
      v = dabs(v)

      s = tri_area(xy1, xy2, xy3)
      H(1) = 3 * v / s

      s = tri_area(xy2, xy3, xy4)
      H(2) = 3 * v / s

      s = tri_area(xy3, xy4, xy1)
      H(3) = 3 * v / s

      s = tri_area(xy4, xy1, xy2)
      H(4) = 3 * v / s

      Return
      End



C ======================================================================
      Subroutine tet_normals(xy1, xy2, xy3, xy4, normal)
C ======================================================================
C Routines calculates unnormalized external normals to tetrahedron faces
C ======================================================================
      implicit none
      Real*8   xy1(3), xy2(3), xy3(3), xy4(3), normal(3, 4)

      Integer  i, n
      Real*8   DotMul, edge(3), vol

C ======================================================================
      Call tri_normal(xy1, xy2, xy3, normal(1,1))
      Call tri_normal(xy2, xy4, xy3, normal(1,2))
      Call tri_normal(xy3, xy4, xy1, normal(1,3))
      Call tri_normal(xy1, xy4, xy2, normal(1,4))

c define orientation
      Do i = 1, 3
         edge(i) = xy4(i) - xy1(i)
      End do

      vol = DotMul(normal(1,1), edge)

      If(vol.GT.0D0) Then
         Do n = 1, 4
            Do i = 1, 3
               normal(i, n) = -normal(i, n)
            End do
         End do
      End if

      Return 
      End



C ======================================================================
      Subroutine tet_dihedrals(xy1, xy2, xy3, xy4, diH)
C ======================================================================
C Routines calculates dihedral angles between tetrahedron faces
C ======================================================================
      implicit none
      Real*8   xy1(3), xy2(3), xy3(3), xy4(3), diH(6), normal(3, 4)

      Integer  i, j, k
      Real*8   calNorm, DotMul, area(4), c, pi

C ======================================================================
      pi = 3.14159265358979D0

c calculate normals and areas
      Call tet_normals(xy1, xy2, xy3, xy4, normal)

      Do i = 1, 4
         area(i) = calNorm(normal(1,i))
      End do

c calculate dihedral angles
      k = 0
      Do i = 1, 4
         Do j = i+1, 4
            k = k + 1
            c = DotMul(normal(1,i), normal(1, j))
            c = c / (area(i) * area(j))
            diH(k) = dacos(-c) / pi * 180  
         End do
      End do 

      Return
      End

