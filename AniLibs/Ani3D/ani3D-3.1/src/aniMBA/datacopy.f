C ======================================================================
      Subroutine copyEdge2Tetrahedron(nP, nE, IPE, Rdata, Edata, iW)
C ======================================================================
C Routines copies data from edges to triangles 
C The working memory is nP + 4*nE + 6*nE  =  nP + 10*nE
C ======================================================================
      implicit none

      Integer  nP, nE, IPE(4, *), iW(*)
      Real*8   Rdata(*), Edata(6, *)

C LOCAL VARIABLES
      Integer  i, n, inEP, iIEP, iIRE, iR, nR

C ======================================================================
      inEP = 1
      iIEP = inEP + nP
      iIRE = iIEP + 4 * nE

      Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))

      Do n = 1, nE
         Do i = 1, 6
            iR = iW(iIRE + 6 * (n-1) + i-1)
            Edata(i, n) = Rdata(iR)
         End do
      End do

      Return
      End

