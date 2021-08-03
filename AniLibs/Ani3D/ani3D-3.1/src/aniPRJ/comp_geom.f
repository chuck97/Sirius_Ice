c ----------------------------------------------------------
      Logical function inter_otr_pl(
     &               ptpl1, ptpl2, ptpl3,
     &               pt1, pt2, eps, pt3)
c ==========================================================
c Function detects the intersection of a segment and a plane
c ==========================================================
      Implicit none
c
c Input
c
      Double precision  ptpl1(3), ptpl2(3), ptpl3(3) 
      Double precision  pt1(3), pt2(3), eps 
c
c Output
c
      Double precision  pt3(3)
c
c External
c      
      External          sgn_pt_plk, pt_pl
      Integer           sgn_pt_plk
      Double precision  pt_pl
c
c Local 
c
      Double precision  kp(4), pt(3), s, s1
      Integer           ipt1, ipt2
c ==========================================================
c
      Call koef_pl (ptpl1,ptpl2,ptpl3,kp)
      
      ipt1=sgn_pt_plk(kp,pt1,eps)
      ipt2=sgn_pt_plk(kp,pt2,eps)
 
      If (ipt1*ipt2.gt.0) then
            inter_otr_pl=.false.
            Return
      End if
 
      If ((ipt1.eq.0d0).and.(ipt2.eq.0d0)) then
            inter_otr_pl=.false.
            Return
      End if
      
      If (ipt1.eq.0d0) then
            inter_otr_pl=.true.
            pt3(1)=pt1(1)
            pt3(2)=pt1(2)
            pt3(3)=pt1(3)
            Return
      End if

      If (ipt2.eq.0d0) then
            inter_otr_pl=.true.
            pt3(1)=pt2(1)
            pt3(2)=pt2(2)
            pt3(3)=pt2(3)
            Return
      End if

      inter_otr_pl=.true.
      s= pt_pl(pt1,kp)  
      pt(1)=pt2(1)-pt1(1)
      pt(2)=pt2(2)-pt1(2)
      pt(3)=pt2(3)-pt1(3)
      s1= pt_pl(pt2,kp)-pt_pl(pt1,kp)

      pt3(1)=pt1(1)-pt(1)*s/s1
      pt3(2)=pt1(2)-pt(2)*s/s1
      pt3(3)=pt1(3)-pt(3)*s/s1
c
      Return
      End
c
c -----------------------------------------------------------
      Integer function sgn_pt_pl(ptpl1, ptpl2, ptpl3, pt, eps)
c ==========================================================
c The function detects in which halfspace lies the point with
c respect to plane given by the three points.
c ==========================================================
      Implicit none
c
c Input
c
      Double precision  ptpl1(3), ptpl2(3),ptpl3(3), pt(3), eps 
c
c External
c     
      External          sgn_pt_plk
      Integer           sgn_pt_plk
c
c Local
c
      Double precision  kp(4)
c =========================================================
c
      Call koef_pl (ptpl1,ptpl2,ptpl3,kp)
      sgn_pt_pl=sgn_pt_plk(kp,pt,eps)
c
      Return
      End
c
c-----------------------------------------------------------
      Integer function sgn_pt_plk (kp, pt, eps)
c =========================================================
c The function detects in which halfspace lies the point with
c respect to plane given by the coefficients
c =========================================================
      Implicit none
c
c Input
c
      Double precision  kp(4), pt(3), eps 
c
c External
c      
      External          pt_pl
      Double precision  pt_pl
c
c Local
c
      Double precision  s
      Integer           i
c =========================================================
c
      s=pt_pl(pt,kp)
 
      If (s.gt.eps) then
        sgn_pt_plk=1
        Return
      End if  
 
      If (s.lt.-eps) then
        sgn_pt_plk=-1
        Return
      End if  
      
      sgn_pt_plk=0
c
      Return
      End
c
c-----------------------------------------------------------
c
      Double precision function pt_pl(pt, kp)
c =========================================================
c The function computes  the right hand side of 
c equation of the  plane given by  coefficients kp
c =========================================================
      Implicit none
c
c Input
c
      Double precision  pt(3),kp(4)
c
c Local
c
      Double precision  s
      Integer           i
c =========================================================
c
      s=kp(4)  
      Do i=1,3
       s = s + kp(i)*pt(i)
      End do
      pt_pl=s/sqrt(kp(1)**2+kp(2)**2+kp(3)**2)
c       
      Return
      End
c
c-----------------------------------------------------------
      Subroutine koef_pl (pt1, pt2, pt3, kp)
c =========================================================
c The routine computes the equation of plane passing through tree points
c =========================================================
      Implicit none
c
c Input
c
      Double precision  pt1(3), pt2(3), pt3(3)
c
c Output
c
      Double precision  kp(4), s
c =========================================================
c
      kp(1)=(pt2(2)-pt1(2))*(pt3(3)-pt1(3))-
     &      (pt2(3)-pt1(3))*(pt3(2)-pt1(2))
 
      kp(2)=-(pt2(1)-pt1(1))*(pt3(3)-pt1(3))+
     &       (pt2(3)-pt1(3))*(pt3(1)-pt1(1))
 
      kp(3)=(pt2(1)-pt1(1))*(pt3(2)-pt1(2))-
     &      (pt2(2)-pt1(2))*(pt3(1)-pt1(1))
 
      kp(4)=-pt1(1)*kp(1)-pt1(2)*kp(2)-pt1(3)*kp(3)
c
      Return
      End
c
c-----------------------------------------------------------
      Subroutine vest_prod (a, b, c)
c =========================================================
c The routine computes the vector product of two vectors          
c =========================================================
      Implicit none
c
c Input
c   
      Double precision  a(3), b(3), c(3)
c =========================================================
c
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
c
      Return
      End
c
c-----------------------------------------------------------
      Subroutine sgn_triple_prod (a, b, c, sgn, eps)
c =========================================================
c The routine computes the sign of triple product of three vectors          
c =========================================================
      Implicit none
c
c Input
c   
      Double precision  a(3), b(3), c(3), eps
c
c Ouput
c
      Integer           sgn
c
c Local
c
      Double precision  s
c =========================================================
c
      s =     a(1)*(b(2)*c(3) - b(3)*c(2))
      s = s + a(2)*(b(3)*c(1) - b(1)*c(3))
      s = s + a(3)*(b(1)*c(2) - b(2)*c(1))

      If (s .gt. eps) then
        sgn = 1
        Return
      End if  
 
      If (s .lt. -eps) then
        sgn = -1
        Return
      End if  
      
      sgn = 0
c
      Return
      End
c
c-----------------------------------------------------------
      Logical function point_in_triang(pt1, pt2, pt3, pt, eps)
c =========================================================
c The function detects: does points belong to the triangle
c =========================================================
      Implicit none
c
c Input
c
      Double precision  pt1(3), pt2(3), pt3(3), pt(3), eps 
c
c Local
c
      Double precision  a(3), b(3), c(3)
      Integer           sgn(3), absssgn, sabssgn, i
c =========================================================
c
      Do i = 1, 3
        a(i) = pt2(i) - pt1(i)
        b(i) = pt3(i) - pt1(i)
      End do
      Call vest_prod (a, b, c)

      Do i = 1, 3
        a(i) = pt1(i) - pt(i)
        b(i) = pt2(i) - pt(i)
      End do
      Call sgn_triple_prod (a, b, c, sgn(1), eps)

      Do i = 1, 3
        a(i) = pt2(i) - pt(i)
        b(i) = pt3(i) - pt(i)
      End do
      Call sgn_triple_prod (a, b, c, sgn(2), eps)

      Do i = 1, 3
        a(i) = pt3(i) - pt(i)
        b(i) = pt1(i) - pt(i)
      End do
      Call sgn_triple_prod (a, b, c, sgn(3), eps)

      absssgn = abs(sgn(1) + sgn(2) + sgn(3))
      sabssgn = abs(sgn(1)) + abs(sgn(2)) + abs(sgn(3))

      point_in_triang = absssgn .eq. sabssgn
c
      Return
      End
c
c-----------------------------------------------------------
      Double precision function dist3d (pt1,pt2)
c =========================================================
c The function detects lenght of the segment in 3d
c =========================================================
      Implicit none
c
c Input
c
      Double precision  pt1(3),pt2(3)
c =========================================================
c
      dist3d =dsqrt((pt1(1)-pt2(1))**2+(pt1(2)-pt2(2))**2
     &                               +(pt1(3)-pt2(3))**2)
c
      Return
      End
c
c-----------------------------------------------------------------
      Subroutine oder_point (nv, np, list1, list2, x, ptet, eps)
c =========================================================
c The routine order points          
c =========================================================
      Implicit none
c
c Input
c   
      Integer           nv, np, list1(*)
      Double precision  x(3,*), ptet(3,4), eps
c
c Externl
c
      External          sgn_pt_pl
      Integer           sgn_pt_pl
c
c Output
c   
      Integer           list2(*)
c
c Local
c      
      Integer           i,j,k,sgn
      Double precision  a(3)
c =========================================================
c
      If (nv .eq. np) Then
        Do j = 1, 4
          sgn = sgn_pt_pl(x(1,list1(1)),x(1,list1(2)),
     &                    x(1,list1(3)),ptet(1,j),eps)
          If (sgn .ne. 0) Then
            Do k = 1, 3
              a(k) = ptet(k, j)
            End do  
            Go to 20
          End if  
        End do
      End if

      Do i = 1, nv
        Do j = 1, np
          If (list1(j) .eq. i) Go to 10
        End do
        Do k = 1, 3
          a(k) = x(k, i)
        End do  
        Go to 20
10      Continue
      End do

20    Continue

      list2(1) = list1(1)
      list1(1) = - list1(1)

      Do i = 1, np - 1

        Do j = 1, np
          If (list1(j) .lt. 0) Go to 40

          Do k = 1, np
            If (k .eq. j) Go to 30
            If (abs(list1(k)) .eq. abs(list2(i))) Go to 30

            sgn = sgn_pt_pl(x(1,list2(i)),x(1,abs(list1(j))),
     &                      a,x(1,abs(list1(k))),eps)
            If (sgn .lt. 0) Go to 40
30          Continue
          End do  
          list2(i+1) = list1(j)
          list1(j) = - list1(j)
          Go to 50

40        Continue
        End do
 
50      Continue
      End do  
c
      Return
      End
c
c-----------------------------------------------------------
      
C ======================================================================
      Subroutine dcentre (nP, xP, x)
C ======================================================================
C  The routine computes coordinates of center mass of cloud of points
C ======================================================================
      Implicit none
      Integer nP
      Real*8  xP(3,*), x(3)
C group (Local variables)
      Integer i, j
C ======================================================================

      Do j = 1, 3
        x(j) = 0d0
      End do

      Do i = 1, nP
        Do j = 1, 3
          x(j) = x(j) + xP(j,i)
        End do
      End do

      Do j = 1, 3
        x(j) = x(j)/nP
      End do

      Return
      End

C ======================================================================
      Subroutine fcentre (nPF, PF, xP, x)
C ======================================================================
C  The routine computes coordinates of faces center mass
C ======================================================================
      Implicit none
      Integer nPF, PF(*)
      Real*8  xP(3,*), x(3)

C group (Local variables)
      Integer i, j
C ======================================================================

      Do j = 1, 3
        x(j) = 0d0
      End do

      Do i = 1, nPF
        Do j = 1, 3
          x(j) = x(j) + xP(j,PF(i))
        End do
      End do

      Do j = 1, 3
        x(j) = x(j)/nPF
      End do

      Return
      End

