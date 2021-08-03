C ================================================================
      Logical Function check1j(i, j)
C ================================================================
C check1j = TRUE if i belongs to the set j(4)}.
C ================================================================
      Integer i, j(4)

      check1j = .TRUE.
      If(i.EQ.j(1)) Goto 1000
      If(i.EQ.j(2)) Goto 1000
      If(i.EQ.j(3)) Goto 1000
      If(i.EQ.j(4)) Goto 1000

      check1j = .FALSE.
 1000 Return
      End



C ================================================================
      Logical Function check2j(i1, i2, j)
C ================================================================
C check2j = TRUE if {i1, i2} is a subset of j(3).
C ================================================================
      Integer j(3)

      check2j = .FALSE.
      If(i1.NE.j(1) .AND. i1.NE.j(2) .AND. i1.NE.j(3)) Goto 1000
      If(i2.NE.j(1) .AND. i2.NE.j(2) .AND. i2.NE.j(3)) Goto 1000

      check2j = .TRUE.
 1000 Return
      End



C ================================================================
      Logical Function check22(i1, i2, j1, j2)
C ================================================================
C check22 = TRUE if pair {i1, i2} coinsides with {j1, j2}.
C ================================================================
      check22 = .FALSE.
      If(i1.NE.j1 .AND. i1.NE.j2) Goto 1000
      If(i2.NE.j1 .AND. i2.NE.j2) Goto 1000

      check22 = .TRUE.
 1000 Return
      End



C ================================================================
      Logical Function check3j(i1, i2, i3, j)
C ================================================================
C check3j = TRUE if {i1, i2, i3} is a subset of j(4).
C ================================================================
      Integer j(4)

      If(i1.EQ.j(1) .OR. i1.EQ.j(2) .OR.
     &   i1.EQ.j(3) .OR. i1.EQ.j(4)) Then
         If(i2.EQ.j(1) .OR. i2.EQ.j(2) .OR.
     &      i2.EQ.j(3) .OR. i2.EQ.j(4)) Then
            If(i3.EQ.j(1) .OR. i3.EQ.j(2) .OR.
     &         i3.EQ.j(3) .OR. i3.EQ.j(4)) Then
               check3j = .TRUE.
               Goto 1000
            End if
         End if
      End if

      check3j = .FALSE.
 1000 Return
      End



C ================================================================
      Logical Function check33(i1, i2, i3, j1, j2, j3)
C ================================================================
C check33 = TRUE if triple {i1, i2, i3} coinsides with {j1, j2, j3}
C ================================================================
      If(i1.EQ.j1 .OR. i1.EQ.j2 .OR. i1.EQ.j3) Then
         If(i2.EQ.j1 .OR. i2.EQ.j2 .OR. i2.EQ.j3) Then
            If(i3.EQ.j1 .OR. i3.EQ.j2 .OR. i3.EQ.j3) Then
               check33 = .TRUE.
               Goto 1000
            End if
         End if
      End if

      check33 = .FALSE.
 1000 Return
      End



C ================================================================
      Logical Function check13(i1, j1, j2, j3)
C ================================================================
C check13 = TRUE if i1 belongs to the set {j1, j2, j3}.
C ================================================================
      check13 = .TRUE.
      If(i1.EQ.j1) Goto 1000
      If(i1.EQ.j2) Goto 1000
      If(i1.EQ.j3) Goto 1000

      check13 = .FALSE.
 1000 Return
      End



C ================================================================
      Logical Function check14(i1, j1, j2, j3, j4)
C ================================================================
C check14 = TRUE if i1 belongs to the set {j1, j2, j3, j4}.
C ================================================================
      check14 = .TRUE.
      If(i1.EQ.j1) Goto 1000
      If(i1.EQ.j2) Goto 1000
      If(i1.EQ.j3) Goto 1000
      If(i1.EQ.j4) Goto 1000

      check14 = .FALSE.
 1000 Return
      End


