C ======================================================================
      Real*8 Function tri_area(xy1, xy2, xy3)
C ======================================================================
C  Routine calculates area of a triangle formed by three points
C ======================================================================
      Real*8 xy1(3), xy2(3), xy3(3)

C Local variables
      Real*8 ax, ay, az, bx, by, bz

      ax = xy1(1) - xy3(1)
      ay = xy1(2) - xy3(2)
      az = xy1(3) - xy3(3)

      bx = xy2(1) - xy3(1)
      by = xy2(2) - xy3(2)
      bz = xy2(3) - xy3(3)

      tri_area = dsqrt((ay * bz - az * by) ** 2 +
     &                 (ax * bz - az * bx) ** 2 +
     &                 (ax * by - ay * bx) ** 2) / 2
      Return
      End



C ======================================================================
      Subroutine tri_normal(xy1, xy2, xy3, normal)
C ======================================================================
C Routine calculates unnormalized normal to the triangle
C ======================================================================
      Real*8 xy1(3), xy2(3), xy3(3), normal(3), a(3), b(3)

      Do i = 1, 3
         a(i) = xy2(i) - xy1(i)
         b(i) = xy3(i) - xy1(i)
      End do

      Call VecMul(a, b, normal)

      Return
      End



C ======================================================================
      Subroutine tri_normal_ext(xy1, xy2, xy3, xy4, normal)
C ======================================================================
C Routine calculates unnormalized exterior normal to the triangle 123
C ======================================================================
      Real*8 xy1(3), xy2(3), xy3(3), xy4(3), normal(3)
      Real*8 DotMul, a(3), b(3), vol

      Do i = 1, 3
         a(i) = xy2(i) - xy1(i)
         b(i) = xy3(i) - xy1(i)
      End do

      Call VecMul(a, b, normal)

c ... fix orientation
      Do i = 1, 3
         a(i) = xy4(i) - xy1(i)
      End do

      vol = DotMul(normal, a)

      If(vol.GT.0D0) Then
         Do i = 1, 3
            normal(i) = -normal(i)
         End do
      End if

      Return 
      End


