c ----------------------------------------------------------
      Subroutine   tet2tet  (
c      
     &                       ptet1, ptet2,
c
     &                       typeintrsec,
c                     
     &                       nv, x, nf, npf, pf)
c ==========================================================
c The subroutine finds intersection of two tetrahedra
c ==========================================================
      Implicit none
c                          
c Input                    
c                          
      Double precision       ptet1(3,4),ptet2(3,4)
c                          
c Output                    
c                          
      Integer                typeintrsec
      Integer                nv
      Double precision       x(3,*)
      Integer                nf,npf(*),pf(*)
c
c External
c      
      External               sgn_pt_pl
      External               dist3d
      External               point_in_triang
      External               inter_otr_pl
      Integer                sgn_pt_pl
      Double precision       dist3d
      Logical                point_in_triang
      Logical                inter_otr_pl
c                          
c Local                    
c                        
      Double precision       eps0, eps
      Parameter              (eps0 = 1d-12)

      Double precision       ptet(3,4,2), xloc(3)
      Double precision       s, s1
      Integer                k, m
      Integer                i1,j1,k1,m1
      Integer                i2,j2
      Integer                nt1, nt2
      Integer                sgn1, sgn2
      Logical                flag
      Integer                ipf

      Integer                nsa, nsp(10)
      Integer                sp(10,10), orsp(10,10)
      Logical                flagsp
c =========================================================
c
c .. inicialization
c
      typeintrsec = 0
      nv = 0
      nf = 0
      ipf = 0 
      nsa    = 0

      Do k = 1, 4
        Do m = 1, 3
          ptet(m,k,1) = ptet1(m,k)
          ptet(m,k,2) = ptet2(m,k)
        End do
      End do
c
c .. define eps
c
      s = dist3d (ptet(1,1,1),ptet(1,2,1))
      Do nt1 = 1, 2
        Do i1 = 1, 3
          Do j1 = i1 + 1, 4
            s1 = dist3d (ptet(1,i1,nt1),ptet(1,j1,nt1)) 
            If (s1 .gt. s) s = s1
          End do  
        End do
      End do
      eps = max(eps0 * s,1d-12)
c
c .. one tetrahedron in other tetrahedron
c
      Do nt1 = 1, 2
        nt2 = 3 - nt1

        Do i1 = 1, 2
          Do j1 = i1 + 1, 3
            Do k1 = j1 + 1, 4

              Do m = 1, 4
                If (m .eq. i1) Go to 10
                If (m .eq. j1) Go to 10
                If (m .eq. k1) Go to 10
                m1 = m
10              Continue
              End do     

              sgn1 = sgn_pt_pl(ptet(1,i1,nt1),
     &                         ptet(1,j1,nt1),  
     &                         ptet(1,k1,nt1),  
     &                         ptet(1,m1,nt1), eps)

              Do m = 1, 4
                sgn2 = sgn_pt_pl(ptet(1,i1,nt1),
     &                           ptet(1,j1,nt1),  
     &                           ptet(1,k1,nt1),  
     &                           ptet(1, m,nt2), eps)
                If (sgn1*sgn2 .eq. -1) Go to 20
              End do
            End do
          End do
        End do

        nv = 4
        Do i1 = 1, 4
          Do j1 = 1, 3
             x(j1, i1) = ptet(j1,i1,nt2)
          End do
        End do

        Do i1 = 1, 2
          Do j1 = i1 + 1, 3
            Do k1 = j1 + 1, 4
              nf = nf + 1
              npf(nf) = 3
              pf(ipf + 1) = i1
              pf(ipf + 2) = j1
              pf(ipf + 3) = k1
              ipf = ipf + 3
            End do
          End do
        End do


        typeintrsec = 4
        Go to 100

20      Continue
      End do       

c
c .. form set of points
c

c
c .. point one tetrahedron inside other tetrahedron
c     
      Do nt1 = 1, 2
        nt2 = 3 - nt1

        Do j2 = 1, 4

          Do i1 = 1, 2
            Do j1 = i1 + 1, 3
              Do k1 = j1 + 1, 4
       
                Do m = 1, 4
                  If (m .eq. i1) Go to 30
                  If (m .eq. j1) Go to 30
                  If (m .eq. k1) Go to 30
                  m1 = m
30                Continue
                End do     
       
                sgn1 = sgn_pt_pl(ptet(1,i1,nt1),
     &                         ptet(1,j1,nt1),  
     &                         ptet(1,k1,nt1),  
     &                         ptet(1,m1,nt1), eps)

                sgn2 = sgn_pt_pl(ptet(1,i1,nt1),
     &                           ptet(1,j1,nt1),  
     &                           ptet(1,k1,nt1),  
     &                           ptet(1,j2,nt2), eps)
                If (sgn1*sgn2 .eq. -1) Go to 40

              End do
            End do
          End do

          Do m = 1, nv
            s = dist3d (ptet(1,j2,nt2),x(1,m))
            If (s .lt. eps) go to 40
          End do  

          nv = nv + 1
          Do m = 1, 3
             x(m, nv) = ptet(m,j2,nt2)
          End do

40        Continue

        End do
      End do

c
c .. point from intersections one tetrahedron with other tetrahedron
c     

      Do nt1 = 1, 2
        nt2 = 3 - nt1

        Do i1 = 1, 2
          Do j1 = i1 + 1, 3
            Do k1 = j1 + 1, 4

              Do i2 = 1, 3
                Do j2 = i2 + 1, 4

                  If (.not. inter_otr_pl(
     &                 ptet(1,i1,nt1), ptet(1,j1,nt1),ptet(1,k1,nt1),
     &                 ptet(1,i2,nt2), ptet(1,j2,nt2), eps, xloc))
     &                 Go to 50

                  If (.not. point_in_triang (
     &                 ptet(1,i1,nt1), ptet(1,j1,nt1),ptet(1,k1,nt1),
     &                 xloc, eps)) Go to 50

                  Do m = 1, nv
                    s = dist3d (xloc,x(1,m))
                    If (s .lt. eps) go to 50
                  End do  

                  nv = nv + 1
                  Do m = 1, 3
                    x(m, nv) = xloc(m)
                  End do  

50                Continue

                End do     
              End do

            End do
          End do
        End do

      End do
c
c .. form set of triangles and polygons
c
      Do i1 = 1, nv - 2
        Do j1 = i1 + 1, nv - 1
          Do k1 = j1 + 1, nv

            Do i2 = 1, nsa
              k = 0
              Do j2 = 1, nsp(i2)
                If (sp(j2,i2) .eq. i1) k = k + 1
                If (sp(j2,i2) .eq. j1) k = k + 1
                If (sp(j2,i2) .eq. k1) k = k + 1
              End do  
              If (k .eq. 3) Go to 80
            End do  

            flagsp = .false.    
            flag = .false.
            Do m = 1, nv
              If (m .eq. i1) Go to 70
              If (m .eq. j1) Go to 70
              If (m .eq. k1) Go to 70

              sgn2 = sgn_pt_pl(
     &                x(1, i1), x(1, j1), x(1, k1), x(1, m), eps)

              If (sgn2 .ne. 0) Then
                If (.not. flag) Then
                  sgn1 = sgn2
                  flag = .true.
                 Else
                  If (sgn2 * sgn1 .eq. -1) Then
                    If (flagsp) nsa = nsa -1
                    Go to 80
                  End if  
                End if 
               Else
                If (.not. flagsp) Then
                  flagsp = .true.
                  nsa = nsa + 1
                  nsp(nsa) = 0
                End if 

                k = nsp(nsa) + 1
                nsp(nsa) = k
                sp(k,nsa) = m
              End if

70            Continue
            End do                            

            If (flagsp) Then
               k = nsp(nsa) 
               sp(k+1, nsa) = i1
               sp(k+2, nsa) = j1
               sp(k+3, nsa) = k1
               nsp(nsa) = k + 3
             Else
               nf = nf + 1
               npf(nf) = 3
               pf(ipf + 1) = i1
               pf(ipf + 2) = j1
               pf(ipf + 3) = k1
               ipf = ipf + 3
            End if

80          Continue
          End do
        End do
      End do

100   Continue
c
c .. orderring of points on faces
c
      If (nv .le. 3) Goto 200
c
c .. triangles
c
      ipf  = 0
      Do i1 = 1, nf
        Call oder_point (nv, 3, pf(ipf+1), 
     &                   orsp(1,1), x, ptet1, eps)
        pf(ipf+1)=orsp(1,1)
        pf(ipf+2)=orsp(2,1)
        pf(ipf+3)=orsp(3,1)
        ipf = ipf + 3
      End do  
c
c .. polygons
c
      If (nsa .eq. 0) Goto 200
      Do i1 = 1, nsa
        Call oder_point (nv, nsp(i1), sp(1,i1), 
     &                   orsp(1,i1), x, ptet1, eps)
        nf = nf + 1
        npf(nf) = nsp(i1)
        Do j1 = 1, nsp(i1)
          ipf = ipf + 1 
          pf(ipf) = orsp(j1,i1)
        End do  
      End do  

200   Continue

      If (nv .eq. 1) typeintrsec = 1
      If (nv .eq. 2) typeintrsec = 2
      If (nv .ge. 3) Then 
        If (nf .eq. 1) Then
          typeintrsec = 3
         Else 
          typeintrsec = 4
        End if
      End if  

      Return
      End
c
c-----------------------------------------------------------
                                        
