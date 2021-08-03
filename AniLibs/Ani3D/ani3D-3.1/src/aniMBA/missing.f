C ================================================================
      Subroutine addBoundaryFaces(
C ================================================================
     &      nP, nF, MaxF, nE,  
     &      XYP, IPF, IPE, lbF, lbE, 
     &      iW)
C ================================================================
C     iW(*) - working memory of size 2 * nP + 3 * nF + 4 * nE
C ================================================================
      include 'makS.fd'
C ================================================================
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*)

      Integer iW(*)

C ================================================================
C group (Functions)
      Logical cmpE

C group (Local variables)
      Integer ip(5)
      Logical flagE, flagF

C ================================================================
      iERR = 0
 
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      iIFP = 1
      iIEP = iIFP + 3 * nF 
      inFP = iIEP + 4 * nE
      inEP = inFP + nP 

c ... checking that [iVface, MaxS] is clear
      Do n = 1, nF
         If(lbF(n).GE.iVface) Call errMes(1011, 'addBoundaryFaces',
     &                  'reserved boundary identificator is used')
      End do


C ... creating an auxiliary structure
      Call backReferences(nP, nF, 3, 3, IPF, iW(inFP), iW(iIFP))
      Call backReferences(nP, nE, 4, 4, IPE, iW(inEP), iW(iIEP))


C ... creating material and missing boundaries
      k = 0
      kMax = MaxS - iMface
      Do n = 1, nE
         Do i1 = 1, 4
            i2 = ip(i1 + 1)
            i3 = ip(i2 + 1)

            ip1 = IPE(i1, n)
            ip2 = IPE(i2, n)
            ip3 = IPE(i3, n)

            flagF = cmpE(ip1, ip2, ip3, iW(iIFP), iW(inFP), 0, iFt)
            flagE = cmpE(ip1, ip2, ip3, iW(iIEP), iW(inEP), n, iEt)

            If(.NOT.flagF .AND. .NOT.flagE) Then
               nF = nF + 1
               If(nF.GT.MaxF) Call errMes(1004, 'addBoundaryFaces',
     &                            'local parameter MaxF is small')

               IPF(1, nF) = IPE(i1, n)
               IPF(2, nF) = IPE(i2, n)
               IPF(3, nF) = IPE(i3, n)
               lbF(nF) = iVface
            End if
         End do
      End do

      Return
      End



C ================================================================
      Subroutine addMaterialFaces(
C ================================================================
     &      nP, nF, MaxF, nE,  
     &      XYP, IPF, IPE, lbF, lbE, 
     &      iW)
C ================================================================
C     iW(*) - working memory of size 2 * nP + 3 * nF + 4 * nE
C ================================================================
      include 'makS.fd'
C ================================================================
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*)

      Integer iW(*)

C ================================================================
C group (Functions)
      Logical cmpE

C group (Local variables)
      Integer ip(5)
      Logical flagE, flagF

      Integer Mlist(2, MaxS-iMface+1)

C ================================================================
      iERR = 0
 
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 4
      ip(5) = 1

      iIFP = 1
      iIEP = iIFP + 3 * nF 
      inFP = iIEP + 4 * nE
      inEP = inFP + nP 

c ... checking that [iVface, MaxS] is clear
      Do n = 1, nF
         If(lbF(n).GE.iVface) Call errMes(1011, 'addMaterialFaces',
     &                  'reserved boundary identificator is used')
      End do


C ... creating an auxiliary structure
      Call backReferences(nP, nF, 3, 3, IPF, iW(inFP), iW(iIFP))
      Call backReferences(nP, nE, 4, 4, IPE, iW(inEP), iW(iIEP))


C ... creating material and missing boundaries
      k = 0
      kMax = MaxS - iMface
      Do n = 1, nE
         Do i1 = 1, 4
            i2 = ip(i1 + 1)
            i3 = ip(i2 + 1)

            ip1 = IPE(i1, n)
            ip2 = IPE(i2, n)
            ip3 = IPE(i3, n)

            flagF = cmpE(ip1, ip2, ip3, iW(iIFP), iW(inFP), 0, iFt)
            flagE = cmpE(ip1, ip2, ip3, iW(iIEP), iW(inEP), n, iEt)

            If(.NOT.flagF .AND. flagE .AND. iEt.GT.n) Then
               mat1 = lbE(n)
               mat2 = lbE(iEt)
               If(mat1.NE.mat2) Then
c  ...  searching for this pair in the list of material interfaces
                  Do i = 1, k
                     If((Mlist(1, i).EQ.mat1 .AND.
     &                   Mlist(2, i).EQ.mat2) .OR.
     &                  (Mlist(1, i).EQ.mat2 .AND.
     &                   Mlist(2, i).EQ.mat1)) Then
                        iCface = iMface + i - 1
                        Goto 1
                     End if
                  End do

c  ...  making the new material interface
                  k = k + 1
                  If(k.GT.kMax) Call errMes(1010, 'addMaterialFaces',
     &                        'not enough memory for material faces')
                  Mlist(1, k) = mat1
                  Mlist(2, k) = mat2
                  iCface = iMface + k - 1

 1                nF = nF + 1
                  If(nF.GT.MaxF) Call errMes(1004, 'addMaterialFaces',
     &                               'local parameter MaxF is small')

                  IPF(1, nF) = IPE(i1, n)
                  IPF(2, nF) = IPE(i2, n)
                  IPF(3, nF) = IPE(i3, n)
                  lbF(nF) = iCface
               End if
            End if
         End do
      End do

      Return
      End



C ======================================================================
      Subroutine global2local(
C ======================================================================
     &           myID, ICE,
c group (Mg)
     &           nP, nF, MaxF, nE, 
     &           XYP, IPF, IPE, lbF, lbE,
     &           nPv, nFv, nEv, IPV, IFV, IEV,
c group (Ml)
     &           nPl, nFl, nEl, 
     &           XYPl, IPFl, IPEl, lbFl, lbEl,
     &           nPvl, nFvl, nEvl, IPVl, IFVl, IEVl, 
c group (I)
     &           nFvi, IPPl, IFFl, 
c group (W)
     &           iW)
C ======================================================================
C  The subroutine extracts a submesh from the global mesh using
C  tetrahedra colored by myID color in array ICE. 
C
C  The fixed triangles are placed at the beginning of corresponding 
C  list. The interfaces triangles are created with utility 
C  addBoundaryFaces and are placed right after the fixed surface 
C  triangles. The number of these triangles equals to nFvi. Note 
C  that these triangles are not added to the list of fixed triangles. 
C
C  The intereface points are also placed at the beginning of list.
C  The user don't need to include them in the list of fixed trianges. 
C  It will be done automatically when the interface triangles will 
C  be added to the list of fixed triangles.
C
C  IPFl(MaxF) - since we need to add interface triangles, nFl may be
C               bigger than nF.
C
C  IPPl(nPl)  - references to a global enumeration necessary to
C               glue submeshes in one global mesh:
C                  0 - point is located outside the DD interfaces;
C                 >0 - splitted points from different subdomains
C                      have the same value of IPPl.
C
C IFFl(nFl)   - references to global enumeration of surface triangles:
C                  0 - the auxiliary interface triangle
C                 >0 - splitted triangles from different meshes have
C                      the same value of IFFl.
C
C  In order to use the information inside IPPl, the interface 
C  triangles should not be modified and be always at the beginning
C  of the corresponding list. The user has to add the interface 
C  triangles (nFvi triangles) to the list of fixed triangles.
C
C  Working memory: iW(4*nP + 5*nF + 4*nE)
C
C  Remarks: The number of interface triangles is correct if and 
C  only if the global mesh is complete, i.e. there are no missing 
C  surface triangles. If necessary, the user has to use utility 
C  addBoundaryFaces before splitting the global mesh.
C
C ======================================================================
      include 'makS.fd'
      include 'color.fd'
      include 'magic.fd'
C ======================================================================
C group (Mg)
      Integer myID, ICE(*)

      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *), lbF(*), lbE(*)
      Integer IPV(*), IFV(*), IEV(*)

C group (Ml)
      Real*8  XYPl(3, *)
      Integer IPEl(4, *), IPFl(3, *), lbFl(*), lbEl(*)
      Integer IPVl(*), IFVl(*), IEVl(*)

C group (I)
      Integer IPPl(*), IFFl(*)

C group (W)
      Integer iW(*)

c group (Local variables)
      Logical cmpE, flag 

      Integer iref(5)
      DATA    iref/1, 2, 3, 4, 1/

C ======================================================================
      iP1  = 0
      iP2  = iP1  + nP
      iF1  = iP2  + nP
      iF2  = iF1  + nF
      iE1  = iF2  + nF
      inEP = iE1  + nE
      iIEP = inEP + nP

      Call makTnode(nP, nP, nE, iW(iP2 + 1), IPE, ICE)


c ... marking points of the subdomain with myID
      Do n = 1, nP
         iW(n) = 0
      End do

      nEl = 0
      Do n = 1, nE
         If(ICE(n).EQ.myID) Then
            nEl = nEl + 1
            iW(iE1 + n) = nEl

            Do i = 1, 4
               iP = IPE(i, n)
               iW(iP) = iP

               IPEl(i, nEl) = iP
            End do

            lbEl(nEl) = lbE(n)
         End if
      End do

      nEvl = 0
      Do n = 1, nEv
         iE = IEV(n)
         If(ICE(iE).EQ.myID) Then
            nEvl = nEvl + 1
            IEVl(nEvl) = iW(iE1 + iE)
         End if
      End do


c ... mapping from boundary faces to elements
      Do n = 1, nF
         iW(iF2 + n) = 0
      End do

      Call backReferences(nP, nF, 3, 3, IPF, iW(inEP+1), iW(iIEP+1))

      Do n = 1, nE
         If(ICE(n).EQ.myID) Then
            Do i1 = 1, 4
               i2 = iref(i1 + 1)
               i3 = iref(i2 + 1)

               jp1 = IPE(i1, n)
               jp2 = IPE(i2, n)
               jp3 = IPE(i3, n)

               If(cmpE(jp1, jp2, jp3, iW(iIEP+1), 
     &                                iW(inEP+1), 0, iF)) Then
                  iW(iF2 + iF) = 1
               End if
            End do
         End if
      End do


c ... coping points
      nPl = 0
      Do m = 1, 2
         Do n = 1, nP
            iPt = iP2 + n
            flag = (m.EQ.1 .AND. iW(iPt).GT.0) .OR. 
     &             (m.EQ.2 .AND. iW(iPt).EQ.0)  

            If(iW(n).NE.0 .AND. flag) Then
               nPl = nPl + 1
               iW(n) = nPl

               Do i = 1, 3
                  XYPl(i, nPl) = XYP(i, n)
               End do
               IPPl(nPl) = iW(iP2 + n)
            End if
         End do
      End do


      nPvl = 0
      Do n = 1, nPv
         If(iW(IPV(n)).NE.0) Then
            nPvl = nPvl + 1
            IPVl(nPvl) = iW(IPV(n))
         End if
      End do


C ... coping fixed faces to on-processor array
      Do n = 1, MaxF
         IFFl(n) = 0
      End do

      Do n = 1, nF
         iW(iF1 + n) = 0
      End do

      nFl = 0
      Do 100 n = 1, nFv
         iF = IFV(n)
         If(iW(iF2 + iF).EQ.0) Goto 100

         nFl = nFl + 1
         iW(iF1 + iF) = nFl
         Do i = 1, 3
            IPFl(i, nFl) = iW(IPF(i, iF))
         End do

         lbFl(nFl) = lbF(iF)
         IFFl(nFl) = iF
 100  Continue


c ... coping surface triangles to the list of interface triangles 
      Do n = 1, nP
         iW(iP2 + n) = iW(n)
      End do

      Do n = 1, nE
         If(ICE(n).NE.myID) Then 
            Do i = 1, 4
               iPt = iP2 + IPE(i, n)
               iW(iPt) = -iW(iPt)
            End do     
         End if
      End do


c ... coping the rest of faces
      Do 300 n = 1, nF
         If(iW(iF2 + n).EQ.0) Goto 300

         nFl = nFl + 1
         iW(iF1 + n) = nFl
         Do i = 1, 3
            IPFl(i, nFl) = iW(IPF(i, n))
         End do

         lbFl(nFl) = lbF(n)
         IFFl(nFl) = n
 300  Continue

      nFvl = 0
      Do n = 1, nFv
         iF = IFV(n)
         iFt = iF1 + iF
         If(iW(iFt).NE.0) Then
            nFvl = nFvl + 1
            IFVl(nFvl) = iW(iFt)
         End if
      End do


c ... computing local coordinates for elements
      Do n = 1, nEl
         Do i = 1, 4
            IPEl(i, n) = iW(IPEl(i, n))
         End do
      End do


c ... computing the interface triangles
      nFlo = nFl
      Call addBoundaryFaces(
     &      nPl, nFl, MaxF, nEl,  
     &      XYPl, IPFl, IPEl, lbFl, lbEl, 
     &      iW(iE1 + 1))

      nFvi = nFl - nFlo


c ... moving the fixed and interface triangles
      m = nFvl + 1  
      Do n = nFl, nFlo + 1, -1
         Do i = 1, 3
            Call swapii(IPFl(i, n), IPFl(i, m))
         End do

         Call swapii(lbFl(n), lbFl(m))
         Call swapii(IFFl(n), IFFl(m))
         m = m + 1
         If(m.GT.nFlo) Goto 500
      End do


c ... changing color of interface triangles (if possible)
 500  nFvi = nFvi

      icfree = iVface
      Do 600 ic = 1, iVface - 1
         Do n = 1, nFl
            If(lbFl(n).EQ.ic) Goto 600
         End do
         icfree = ic
         Goto 700
 600  Continue

 700  Continue
      Do n = 1, nFl
         If(lbFl(n).EQ.iVface) lbFl(n) = icfree
      End do

      Return
      End



C ================================================================
      Subroutine local2global(
C ================================================================
     &           nMeshes, 
c group (Mg)
     &           MaxP, MaxF, MaxE, 
     &           nP, nF, nE, 
     &           XYP, IPF, IPE, lbF, lbE,
     &           nPv, nFv, nEv, IPV, IFV, IEV,
c group (Ml)
     &           nPl, nFl, nEl, 
     &           XYPl, IPFl, IPEl, lbFl, lbEl,
     &           nPvl, nFvl, nEvl, IPVl, IFVl, IEVl,
c group (I)
     &           IPPl, IFFl, 
c group (W)
     &           IPs, IPw, IFw, IEw)
C ================================================================
C  The subroutine gathers submeshes into a global mesh using
C  references to global enumeration of mesh points kept in array 
C  IPPl(*) and the global enumeration of faces in array IFFl(*).
C
C  The meshes are kept as continuous lists of data. For instance,
C  array XYP(1:3, 1:nP1) keeps coordinates of mesh points in the 
C  first grid. Continuation of the array, XYP(1:3, nP1 : nP1 + nP2),
C  keeps coordinates of the mesh points of the second grid. And
C  so on. The array having such a structure are:
C
C    XYPl(3, *) - mesh coordinates
C    IPFl(3, *) - connectivity list for triangles
C    IPEl(4, *) - connectivity list for tetrahedra
C
C    IPVl(*) - list of fixed points
C    IFVl(*) - list of fixed triangles
C    IEVl(*) - list of fixed tetrahedra  
C
C    IPPl(*) - references to global enumeration of mesh point
C    IFFl(*) - references to global enumeration of mesh triangles 
C
C
C  The necessary information to extract submesh data are kept in
C  the following arrays:
C
C    nPl(nMeshes)  - numbers of points
C    nFl(nMeshes)  - numbers of faces
C    nEl(nMeshes)  - numbers of tetrahedra 
C
C    nPvl(nMeshes) - numbers of fixed points
C    nFvl(nMeshes) - numbers of fixed triangles
C    nEvl(nMeshes) - numbers of fixed tetrahedra 
C
C  Working memory: IPs(MaxP), IPw(MaxP), IFw(MaxF), IEw(4*MaxE)
C 
C  Remark: Array IPFl is destroyed in the algorithm.
C ================================================================
C group (Mg)
      Real*8  XYP(3, *)
      Integer IPF(3, *), IPE(4, *) 
      Integer lbF(*), lbE(*)

      Integer IPV(*), IFV(*), IEV(*)

C group (Ml)
      Integer nPl(*), nFl(*), nEl(*)

      Real*8  XYPl(3, *)
      Integer IPFl(3, *), IPEl(4, *)
      Integer lbFl(*), lbEl(*)

      Integer nPvl(*), nFvl(*), nEvl(*), IPVl(*), IFVl(*), IEVl(*)

C group (I)
      Integer IPPl(*), IFFl(*)

C group (W)
      Integer IPs(*), IPw(*), IFw(*), IEw(*)

c group (Local variables)
      Logical cmpE

C ================================================================
      mP = 0
      mF = 0
      mE = 0     
      Do n = 1, nMeshes
         mP = mP + nPl(n)
         mF = mF + nFl(n)
         mE = mE + nEl(n)
      End do
      If(mP.GT.MaxP) Call errMes(1003, 'local2global',
     &                   'local parameter MaxP is small')
      If(mF.GT.MaxF) Call errMes(1004, 'local2global',
     &                   'local parameter MaxF is small')
      If(mE.GT.MaxE) Call errMes(1006, 'local2global',
     &                   'local parameter MaxE is small')

      Do n = 1, mP
         IPs(n) = 0
         IPw(n) = 0
      End do

c ... counting points on interfeices 
      kP = 0
      Do n = 1, mP
         If(IPPl(n).NE.0) Then
            If(IPs(IPPl(n)).EQ.0) Then
               kP = kP + 1
               IPs(IPPl(n)) = kP
            End if
         End if
      End do

c ... counting the rest of points
      Do n = 1, mP
         If(IPPl(n).EQ.0) Then
            kP = kP + 1
            IPw(n) = kP
         Else
            IPw(n) = IPs(IPPl(n))
         End if
      End do


c ... making points
      nP = kP
      Do n = 1, mP
         Do i = 1, 3
            XYP(i, IPw(n)) = XYPl(i, n)
         End do
      End do


c ... making vertices
      Do n = 1, nP
         IPs(n) = 0
      End do

      mP = 0
      mV = 0
      Do k = 1, nMeshes
         Do n = 1, nPvl(k)
            mV = mV + 1
            IPs(IPw(mP + IPVl(mV))) = 1
         End do
         mP = mP + nPl(k)
      End do

      nPv = 0
      Do n = 1, nP
         If(IPs(n).NE.0) Then
            nPv = nPv + 1
            IPV(nPv) = n
         End if
      End do


c ... making faces
      mP = 0
      mF = 0
      Do k = 1, nMeshes
         Do n = 1, nFl(k)
            mF = mF + 1

            IFw(mF) = 0

            Do i = 1, 3
               IPFl(i, mF) = IPw(mP + IPFl(i, mF))
            End do
         End do
         mP = mP + nPl(k)
      End do

      Call backReferences(nP, mF, 3, 3, IPFl, IPs, IEw)

      Do n = 1, mF
         If(IFw(n).EQ.0) Then
            ip1 = IPFl(1, n)
            ip2 = IPFl(2, n)
            ip3 = IPFl(3, n)

            If(cmpE(ip1, ip2, ip3, IEw, IPs, n, iF)) Then
               If(iF.GT.n) IFw(n) = iF
            End if
         End if
      End do

      nF = 0
      Do n = 1, mF
         If(IFw(n).EQ.0 .AND. IFFl(n).GT.0) Then
            nF = nF + 1
            IFw(n) = nF

            Do i = 1, 3
               IPF(i, nF) = IPFl(i, n)
            End do
            lbF(nF) = lbFl(n)
         Else
            IFw(n) = 0
         End if
      End do


c ... making fixed faces
      Do n = 1, nF
         IEw(n) = 0
      End do

      mF = 0
      mV = 0
      Do k = 1, nMeshes
         Do n = 1, nFvl(k)
            mV = mV + 1
            IEw(IFw(mF + IFVl(mV))) = 1
         End do
         mF = mF + nFl(k)
      End do

      nFv = 0
      Do n = 1, nF
         If(IEw(n).NE.0) Then
            nFv = nFv + 1
            IFV(nFv) = n
         End if
      End do


c ... making elements
      mP = 0
      nE = 0
      Do k = 1, nMeshes
         Do n = 1, nEl(k)
            nE = nE + 1
            
            Do i = 1, 4
               IPE(i, nE) = IPw(mP + IPEl(i, nE))
            End do
            lbE(nE) = lbEl(nE)
         End do
         mP = mP + nPl(k)
      End do


c ... making fixed elements
      nEv = 0

      mE = 0
      mV = 0
      Do k = 1, nMeshes
         Do n = 1, nEvl(k)
            nEv = nEv + 1
            IEV(nEv) = mE + IEVl(mV + n)
         End do
         mE = mE + nEl(k)
         mV = mV + nEvl(k)
      End do

      Return
      End



C ================================================================
      Subroutine makTnode(nP, MaxP, nE, IPP, IPE, ICE)
C ================================================================
C   The T-nodes are marked. The original connectivity
C   list IPE is used. 
C ================================================================
      Integer IPP(*), IPE(4, *), ICE(*)

      Do n = 1, MaxP
         IPP(n) = 0
      End do

      Do n = 1, nE
         ict = ICE(n) + 1
         Do i = 1, 4
            iP = IPE(i, n)
            icn = IPP(iP)
            If(icn.EQ.0) Then
               IPP(iP) = ict
            Else If(icn.NE.ict) Then
               IPP(iP) = -1
            End if
         End do
      End do


      Do n = 1, nP
         If(IPP(n).EQ.-1) Then
            IPP(n) = n
         Else
            IPP(n) = 0
         End if
      End do

      Return
      End

