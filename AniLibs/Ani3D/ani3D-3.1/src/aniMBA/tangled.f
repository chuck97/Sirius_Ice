C ======================================================================
      Subroutine updQa(n, XYP, IPE, IEE, qE)
C ======================================================================
C The initial quality modification for tangled elements and 
C their closest (common face) neighboors.
C ======================================================================
      Real*8  XYP(3, *), qE(*)
      Integer IPE(5, *), IEE(4, *)

C (Local variables)
      Integer ip(5)
      Real*8  vv, mutualOrientation
      Logical check33

C ======================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      Do 100 i1 = 1, 4
         iE = IEE(i1, n)
         If(iE.LE.0) Goto 100

         i2 = ip(i1 + 1)
         i3 = ip(i2 + 1)
         i4 = ip(i3 + 1)

         iP1 = IPE(i1, n)
         iP2 = IPE(i2, n)
         iP3 = IPE(i3, n)
         iP4 = IPE(i4, n)

         Do j1 = 1, 4
            j2 = ip(j1 + 1)
            j3 = ip(j2 + 1)

            jP1 = IPE(j1, iE)
            jP2 = IPE(j2, iE)
            jP3 = IPE(j3, iE)

            If(check33(iP1, iP2, iP3, jP1, jP2, jP3)) Then
               j4  = ip(j3 + 1)
               jP4 = IPE(j4, iE)

               vv = mutualOrientation(
     &              XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), 
     &                           XYP(1, iP4), XYP(1, jP4))

               If(vv.GE.0D0) Then
                  qE(n)  = -dabs(qE(n))
                  qE(iE) = -dabs(qE(iE))
               End if
               Goto 100
            End if
         End do
 100  Continue

      Return
      End


C ======================================================================
      Subroutine updQb(nEs, lE, iEs, XYP, IPEs, qEs)
C ======================================================================
C A dynamic quality modification for tangled elements inside
C a super-element.
C
C Remark: inefficient, time-consuming, but robust algorithm.
C ======================================================================
      Real*8  XYP(3, *), qEs(*)
      Integer iEs(*), IPEs(5, *)

C group (Local variables)
      Integer ip(5)
      Real*8  vv, mutualOrientation
      Logical check33

C ======================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      Do 100 i1 = 1, 4
         i2 = ip(i1 + 1)
         i3 = ip(i2 + 1)
         i4 = ip(i3 + 1)

         iP1 = IPEs(i1, nEs)
         iP2 = IPEs(i2, nEs)
         iP3 = IPEs(i3, nEs)
         iP4 = IPEs(i4, nEs)

         Do 20 k = 1, lE
            If(iEs(k).LT.0) Goto 20
            If(k.EQ.nEs)    Goto 20

            Do j1 = 1, 4
               j2 = ip(j1 + 1)
               j3 = ip(j2 + 1)

               jP1 = IPEs(j1, k)
               jP2 = IPEs(j2, k)
               jP3 = IPEs(j3, k)

               If(check33(iP1, iP2, iP3, jP1, jP2, jP3)) Then
                  j4  = ip(j3 + 1)
                  jP4 = IPEs(j4, k)

                  vv = mutualOrientation(
     &                 XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), 
     &                              XYP(1, iP4), XYP(1, jP4))

                  If(vv.GE.0D0) Then
                     qEs(nEs) = -dabs(qEs(nEs))
                     qEs(k)   = -dabs(qEs(k))
                  End if

                  Goto 100
               End if
            End do
 20      Continue
 100  Continue

      Return
      End




