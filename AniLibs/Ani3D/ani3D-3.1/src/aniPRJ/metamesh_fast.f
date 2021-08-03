C =====================================================================
      Subroutine MetaMesh(nv, vrt, nt, tet, 
     &                    nv2, vrt2, nt2, tet2,
     &                    nv12, nv12max, vrt12, 
     &                    nt12, nt12max, tet12, parents,
     &                    MaxWi, MaxWr, iW, rW, iERR)
C =====================================================================
C    nv,   vrt(2,nv1), nt,   tet(3,nt1) - the first mesh
C    nv2, vrt2(2,nv2), nt2, tet2(3,nt2) - the second mesh
C    nv12, vrt12(2,nv12), nt12, tet12(3,nt12) - intersection of two meshes
C    parents(2,nt12) - two parents of new tetrahedra
C
C    rW(MaxWr) - Real*8 working memory of size greater than 5*nv2 ??
C    iW(MaxWi) - Integer working memory of size greater than 11*nv2 + 16*nt2
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

      Integer  HostTetForXYZ
      External HostTetForXYZ

c LOCAL VARIABLEs
      Real*8   XYZ(3), XYP(2, 6), v,v1,dV
      Integer  i, j, k, k1, l, n, n2, n2t, m, l1, l2, l2t, pc
      Integer  iLST, iMRK, iIEE, inEP, iIEP, iLIN, iEnd, iControl
      Integer  rLIN, nOverlaps
      Logical  flagOVERLAP
      Logical  ENCLOSE
      External ENCLOSE
c intersection
      Integer  typeintrsec
      Real*8   ptet1(3,4), ptet2(3,4)
      Integer  nvint, nfint, npfint(10), pfint(50)
      Real*8   xint(3,50), xdc(3), xfc(3)
C =====================================================================
      
      iERR = 0  ! reserved for future

c ... memory management
      iLST = 1
      iMRK = iLST + nt2
      iIEE = iMRK + nt2
      inEP = iIEE + 4 * nt2
      iIEP = inEP + nv2
      iEnd = iIEP + 6 * nt2

      iLIN = inEP
      rLIN = 1

      If(iEnd + 10*nv2+4*nt2.GT.MaxWi) Then
         iERR = 1001 
         Stop 9901
      End if
      If(5*nv2.GT.MaxWr) Then
         iERR = 1002
         Stop 9902
      End if

c ... create a map E2 -> E2
      Call listE2E(nv2, nt2, tet2, iW(iIEE), iW(inEP), iW(iIEP))

c ... populate the octtree sructure
      Call scale2Meshes(nv, vrt, nv2, vrt2, .TRUE.)

      iControl = 100 
      XYZ(1) = vrt2(1,1)
      XYZ(2) = vrt2(2,1)
      XYZ(3) = vrt2(3,1)
      n2 = HostTetForXYZ(nt2, tet2, nv2, vrt2, XYZ,  
     &     iW(iLIN), MaxWi-iLIN, rW(rLIN), MaxWr-rLIN, iControl)

c ... clean the list
      Do n2 = 1, nt2
         iW(iMRK + n2-1) = 0
      End do
 
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
 
c find tetrahedron #n2 in mesh 2 that intersects tetrahedron #n in mesh 1
         Do i = 1, 4
            Do j = 1, 3
               ptet1(j,i) = vrt(j, tet(i,n))
            End do  
         End do  

         Do j = 1, 3
            XYZ(j) = (ptet1(j,1) + ptet1(j,2) +
     &                ptet1(j,3) + ptet1(j,4)) / 4
         End do  

         iControl = 200 
         n2 = HostTetForXYZ (nt2, tet2, nv2, vrt2, XYZ,  
     &        iW(iLIN), MaxWi-iLIN, rW(rLIN), MaxWr-rLIN, iControl)
         if (.not.ENCLOSE( XYZ, vrt2, tet2(1, n2), 1d-6)) then
            write(*,*) (XYZ(j), j=1,3),' is not in tet ',n2
            stop
         end if

         l1 = 1
         l2 = 1
         iW(iLST) = n2 
         iW(iMRK + n2-1) = 1

 100     Continue

c loop over neighboors of the previous tetrahedron
         flagOVERLAP = .FALSE.
         Do l = l1, l2
            n2 = iW(iLST + l-1)

            Do i = 1, 4
              Do j = 1, 3
                ptet2(j,i) = vrt2(j, tet2(i,n2))
              End do  
            End do  

            Call tet2tet(ptet1, ptet2, typeintrsec,        
     &                   nvint, xint,
     &                   nfint, npfint, pfint)

            If (typeintrsec .le. 3) Goto 200

            flagOVERLAP = .TRUE. 
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
c working with intrsection
c
c add the points
            Call dcentre(nvint, xint, xdc)

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
 
200         Continue
         End do
 
c add more neighbors if the intersection is not complete
c otherwise, clean the marked tetrahedra
         If(flagOVERLAP) Then
            l2t = l2
            Do l = l1, l2t
               n2 = iW(iLST + l-1)
 
               Do i = 1, 4
                  n2t = iW(iIEE + 4*(n2-1) + i-1)
                  If(n2t.GT.0 .AND. iW(iMRK + n2t-1).EQ.0) Then
                     l2 = l2 + 1
                     iW(iLST + l2-1)  = n2t
                     iW(iMRK + n2t-1) = 1
                  End if
               End do
            End do
 
            l1 = l2t + 1
            Goto 100
         Else
            dV = dV + v1
            Do l = 1, l2
               n2 = iW(iLST + l-1)
               iW(iMRK + n2-1) = 0
            End do
         End if
      End do

c scale the meshes back (the flag should be false)
      Call scale2Meshes(nv, vrt, nv2, vrt2, .FALSE.)
      Call scale2Meshes(nv12, vrt12, 0, vrt2, .FALSE.)

      m = 100 * nOverlaps / (nt * nt2)
 
      Write(*,'(A,2I7)') 'PRJ: 1st mesh: nv/nt:', nv,  nt
      Write(*,'(A,2I7)') '     2nd mesh: nv/nt:', nv2, nt2
      Write(*,'(A,2I7,A,I8,A,I2,A)') '     metamesh: nv/nt:',
     &                         nv12,nt12,
     &         '   after', nOverlaps, ' overlaps (', m, '%)'
      Write(*,'(A,E12.4)') '   volume mismatch is', dV

      Return
      End



C ================================================================
      Integer Function HostTetForXYZ(
C ================================================================
     &                 NT, tet, NV, vrt,  xyz, 
     &                 imem,Nimem, dmem,Ndmem, iControl )
C ================================================================
C  Function HostTetForXYZ returns the index of tetrahedron which contains point XYZ.
C
C  Parameters of subroutine in order to locate (without typing):
C       NT         - number of tetrahedra
C       tet(4, NT) - list of tetrahedra ( 4 vertices )
C       NV         - number of nodes (vetices) of  triangulation
C       vrt(3, NV) - coords of vertices
C
C       xyz(3)     - coords of the point
C
C       imem(Nimem) - auxiliary array of integers of length Nimem
C       dmem(Ndmem) - auxiliary array of double precision numbers of 
C                     length Ndmem
C
C       iContol - 3-digit integer representing 2 control parameters
C            1st digit = 1 means "Initialisation is needed"
C                      = 2 means "Initialisation is not needed"
C                      otherwise the input is erroneous
C            2,3 digits = 00 means "No output information is provided"
C                  ab   > 00 means "Important warnings are written to 
C                                   the channel (NUNIT=ab )"
C ================================================================
C  REMARK: 
C    It is assumed that the given point belongs to the domain covered by the mesh.
C    If for the given point there is no an 
C    element containing it within tolerance PREC=10^{-6} , then the ERROR 
C    is assumed. The code could not find an element within PREC tolerance 
C    for MaxTrials.
C ================================================================
      implicit none
      include 'lintrp.fd'
C ================================================================
      integer NT, NV, Nimem, Ndmem
      integer tet(4, NT)
      double precision vrt(3, NV)
      double precision xyz(3)
      integer imem(*)
      double precision dmem(*)
      integer iControl
      
      double precision h(MaxH)
      integer i, NQT, QT, REF, ITET
      logical flag
      integer NUNIT, iHostTet

C ================================================================
c ... decoding control parameters
      i = iControl/100
      if ( i .eq. 1 ) then
         flag = .true.
      else if ( i .eq. 2 ) then
         flag = .false.
      else
         Call errMes(6101, 'HostTetForXYZ',
     .              'Wrong 1st digit in iControl')
      end if 

      NUNIT = mod(iControl,100)


c ... distributing working memory
      if(flag) then
         QT = 5
         call INITQT( NQT, imem(QT), h, dmem(MaxH + 1) )
         do i = 1, NV
            call DROPS( NQT, imem(QT), Nimem-4, h,dmem(MaxH + 1),vrt,i)
         enddo
         REF  = QT + 8 * NQT
         ITET = REF + NV + 1
      else
         QT  = imem(1)
         NQT = imem(2)
         REF = imem(3)
         ITET =imem(4)
         
         do i = 1, MaxH
            h(i) = dmem(i)
         end do
      end if
      
      if (ITET+4*NT.gt.Nimem) then
         Call errMes(1001, 'HostTetForXYZ',
     .              'not enough memory for Integer arrays')
      end if
      if (MAXH+3*NQT.gt.Ndmem) then
         Call  errMes(1002, 'HostTetForXYZ',
     .        'not enough memory for Real*8 arrays')
      end if


c ... calling the main module 
      call RESTORE2(NT, tet, NV, vrt, xyz, iHostTet, 
     .              imem(QT), dmem(MaxH + 1),
     .              imem(REF), imem(ITET),
     .              h, flag, imem(ITET+4*NT), Nimem-ITET-4*NT,
     .              NUNIT ) 

      HostTetForXYZ = iHostTet

c ... saving memory distribution
      if(flag) then
         imem(1) = QT
         imem(2) = NQT
         imem(3) = REF
         imem(4) = ITET
         
         do i = 1, MaxH
            dmem(i) = h(i)
         end do
      end if
      return
      end


C ================================================================
      Subroutine RESTORE2(NT, tet, NV, vrt, xyz, iHostTet, 
     .                    QT, XYZc, ref, itet,
     .                    h, flag, buf, Nbuf, NUNIT )
C ================================================================
C  Code builds the transposed list of tetrahedra grid and restores
C  the values of vector-function F in the points of interest with the aid
C  of the octtree structure
C ================================================================
      implicit none
      include 'lintrp.fd'
C ================================================================
c Predefined tolerance for checking whether the point belongs to an element
      double precision PREC
      parameter( PREC = 1D-6 )

c Predefined number of possible relaxations of the above tolerance (6 orders of 10, here)
      integer    MaxRelax
      parameter( MaxRelax = 7 )

C ================================================================
      integer NV, NT, iHostTet
      integer QT(2, 2, 2, *), tet(4, NT), ref(*), itet(*)
      double precision XYZc(3, *)
      double precision vrt(3, NV), xyz(3)
      double precision h(*)
      logical          flag
      integer          NUNIT
      integer          buf(*)
      integer          Nbuf

C ================================================================
C (Local variables)
      integer          i, j, idx, ip(4), BASETET

C ================================================================
      if(flag) then
         do i = 1, NV+1
            ref(i) = 0
         enddo
         do i = 1, NT
            do j = 1, 4
               idx = tet(j, i)
               ref(idx) = ref(idx) + 1
            enddo
         enddo
         ref(1) = ref(1) + 1
         do i = 2, NV+1
            ref(i) = ref(i-1) + ref(i)
         enddo
         do i = 1, NT
            do j = 1, 4
               idx = tet(j, i)
               ref(idx) = ref(idx) - 1
               itet(ref(idx)) = i
            enddo
         enddo
      end if
      
 11      idx = BASETET( QT, XYZc, vrt, NT, tet, ref, itet, xyz,
     &                  h,  buf, Nbuf, PREC )
         if(idx.le.0) then
            if (NUNIT.GT.0) then
               write(NUNIT,22) xyz(1),xyz(2),xyz(3)
               write(NUNIT,23) PREC
               write(NUNIT,24)
            end if
            Call errMes(6104, 'HostTetForXYZ',
     .           'Failed to find element within PREC tolerance')
         end if
         
         iHostTet = idx

 22   Format('the point ',3F9.5,' is probably out of the mesh')
 23   Format('After MaxTrials failed to find a host within ',F8.6,
     &       ' tolerance') 
 24   Format('NEED TO REINSTALL constants in lintrp.fd')

      return
      end




