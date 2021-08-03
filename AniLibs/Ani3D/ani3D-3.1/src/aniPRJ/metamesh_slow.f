C =====================================================================
      Subroutine MetaMesh_Slow(nv, vrt, nt, tet, 
     &                         nv2, vrt2, nt2, tet2,
     &                         nv12, nv12max, vrt12, 
     &                         nt12, nt12max, tet12, parents,
     &                         MaxWi, MaxWr, iW, rW, iERR)
C =====================================================================
      implicit none

      External calVol 
      Real*8   calVol 
 
      Integer  nv, nt, tet(4, *)
      Integer  nv2, nt2, tet2(4, *)
      Real*8   vrt(3, *), vrt2(3, *)

      Integer  nv12, nv12max, nt12, nt12max
      Integer  tet12(4, *), parents(2, *)
      Real*8   vrt12(3, *)

      Integer  iW(*), MaxWi, MaxWr, iERR
      Real*8   rW(*)


c LOCAL VARIABLEs
      Real*8   v,v1,dV
      Integer  i, j, k, k1, l, n, n2, n2t, m,  pc
      Integer  nOverlaps
c intersection
      Integer  typeintrsec
      Real*8   ptet1(3,4), ptet2(3,4)
      Integer  nvint, nfint, npfint(10), pfint(50)
      Real*8   xint(3,50), xdc(3), xfc(3)
C =====================================================================
      iERR = 0  ! reserved for future

c ... loop over elements of the first mesh and create a meta-mesh
      nv12 = 0
      nt12 = 0
      nOverlaps = 0
      dV = 0
 
      Do n = 1, nt
         v1 = dabs(calVol(vrt(1,tet(1,n)),
     &                    vrt(1,tet(2,n)),
     &                    vrt(1,tet(3,n)),
     &                    vrt(1,tet(4,n))))
       
         Do i = 1, 4
            Do j = 1, 3
               ptet1(j,i) = vrt(j, tet(i,n))
            End do  
         End do  

c check all tetrahedra in mesh 2 for intersection with tetrahedron #n in mesh 1
         Do n2 = 1, nt2

            Do i = 1, 4
               Do j = 1, 3
                  ptet2(j,i) = vrt2(j, tet2(i,n2))
               End do  
            End do  

            Call tet2tet(ptet1, ptet2, typeintrsec,        
     &                   nvint, xint,
     &                   nfint, npfint, pfint)

            If (typeintrsec .le. 3) Goto 200

            nOverlaps = nOverlaps + 1

c check of memory
            m = 0
            Do i = 1, nfint
              m = m + npfint(i)
            End do

            If(nv12 + nvint + m + 1 .GT. nv12max) Then  
                Call errMesFEM(1001, 'metamesh', 
     &                           'Not enough memory for meta-mesh')
            End if

            If(nt12 + m.GT.nt12max) Then  
               Call errMesFEM(1001, 'metamesh', 
     &                        'Not enough memory for meta-mesh')
            End if
c
c working with intersection
c
c add the points
            Call dcentre (nvint, xint, xdc)

            nv12 = nv12 + 1
            pc = nv12
            Do j = 1, 3
               vrt12(j, nv12) = xdc(j) 
            End do

            Do i = 1, nvint
               nv12 = nv12 + 1
               Do j = 1, 3 
                  vrt12(j, nv12) = xint(j, i) 
               End do  
            End do 
 
c add the tetrahedra
            m = 1
            Do i = 1, nfint
              Call fcentre(npfint(i), pfint(m), xint, xfc)
              nv12 = nv12 + 1
              Do j = 1, 3
                 vrt12(j, nv12) = xfc(j) 
              End do
              
              Do k = 1, npfint(i)
                If (k .lt. npfint(i)) Then
                   k1 = k + 1
                Else
                   k1 = 1
                End if

                nt12 = nt12 + 1
                tet12(1, nt12) = pc
                tet12(2, nt12) = nv12 
                tet12(3, nt12) = pc + pfint(m + k - 1)
                tet12(4, nt12) = pc + pfint(m + k1 - 1)

                parents(1, nt12) = n
                parents(2, nt12) = n2

                v = dabs(calVol(vrt12(1,tet12(1,nt12)),
     &                          vrt12(1,tet12(2,nt12)),
     &                          vrt12(1,tet12(3,nt12)),
     &                          vrt12(1,tet12(4,nt12))))
                v1 = v1 - v
              End do
              m = m + npfint(i)
            End do
c 
200         Continue

         End do
         dV = dV + v1
      End do
c 
      m = 100 * nOverlaps / (nt * nt2)
 
      Write(*,'(A,2I7)') 'PRJ: 1st mesh: nv/nt:', nv,  nt
      Write(*,'(A,2I7)') '     2nd mesh: nv/nt:', nv2, nt2
      Write(*,'(A,2I7,A,I8,A,I2,A)') '     metamesh: nv/nt:',
     &                         nv12,nt12,
     &         '   after', nOverlaps, ' overlaps (', m, '%)'
      Write(*,'(A,E12.4)') '   volume mismatch is', dV

      Return
      End




