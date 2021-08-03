C ======================================================================
      Subroutine countBadElements(nE, L1E, L2E, qE, Quality, 
     &                            nLines, output, flagFILE)
C ======================================================================
      include 'output.fd'
C ======================================================================
      Integer L1E(2, *), L2E(*)
      Real*8  qE(*), Quality

      Character*(LineLenght) output(*)
      Logical flagFILE

c group (Local variables)
      Integer nBad(10), mE
      Real*8  q
      Character*(LineLenght) message

C ======================================================================
      Do n = 1, 10
         nBad(n) = 0
      End do

      mE = 0
      icnt = 0
      iE = L2E(1)
      Do n = 1, nE
         If(qE(iE).LT.Quality) icnt = icnt + 1

         Do k = 1, 9
            q  = dlog10(qE(iE))

            If(-k.LE.q .AND. q.LT.1-k) Then
               nBad(k) = nBad(k) + 1
               mE = mE + 1
               Goto 500
            End if
         End do
 500     iE = L1E(2, iE)
      End do
      nBad(10) = nE - mE


      If(flagFILE) Then
         Write(message,6005)
         Call addOut(nLines, message, output)
         Write(message,6006) (1D-1 * i, i = 1, 10)
         Call addOut(nLines, message, output)
         Write(message,6008) (nBad(i), i = 1, 10)
         Call addOut(nLines, message, output)
         Write(message,6007) icnt
         Call addOut(nLines, message, output)
      Else
         Write(*, 5005) (i, i = 1, 9),
     &                  (nBad(i), i = 1, 10), icnt
      End if


c ... screen output
 5005 Format(
     & ' *** Distribution of tetrahedra by their quality',/,
     &   9('    1e-',I1),'   1e-10',/, 10I8,/, 
     & ' *** Number of bad tetrahedrons =', I8/)

c ... parallel output
 6005 Format(' *** Distribution of tetrahedra by their quality')
 6006 Format(9('    1e-',I1),'   1e-10')
 6008 Format(10I8)
 6007 Format(' *** Number of bad tetrahedrons =', I8)

      Return
      End



C ======================================================================
      Subroutine addOut(nLines, message, output)
C ======================================================================
      include 'output.fd'
C ======================================================================
      Character*(*) message, output(*)

      nLines = min(nLines + 1, MaxLines)
      output(nLines) = message

      Return
      End



