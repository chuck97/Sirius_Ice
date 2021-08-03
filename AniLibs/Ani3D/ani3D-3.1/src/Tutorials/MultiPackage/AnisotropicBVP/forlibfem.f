C ================================================================
C  Diffusion tensor             
C ================================================================
      Integer Function Ddiff(x,y,z, label, DDATA,IDATA, iSYS,Coef)
      include 'fem3Dtet.fd'
C ================================================================
      Real*8  x, y, z, DDATA(*), Coef(9, 9)
      Integer label, IDATA(*), iSYS(*)
C ================================================================
      iSYS(1) = 3
      iSYS(2) = 3

      Do i = 1, 3
         Do j = 1, 3
            Coef(i, j) = 0D0
         End do
         Coef(i, i) = DDATA(1)
      End do
      Coef(3, 3) = DDATA(2)


      Ddiff = TENSOR_SYMMETRIC

      Return
      End



C ================================================================
C Boundary condition
C ================================================================
      Integer Function Dbc(x, y, z, label, DDATA,IDATA, iSYS, eBC)
      Include 'assemble.fd'
C ================================================================
      Real*8  x, y, z, DDATA(*), eBC(*)
      Integer label, IDATA(*), iSYS(*)
C ================================================================
      iSYS(1) = 1
      iSYS(2) = 1

      If (label.eq.1.or.label.eq.2) then
          eBC(1) = 0D0
          Dbc = BC_NEUMANN
      Else
      eBC(1) = 0D0
      Dbc = BC_DIRICHLET
      End if

      Return
      End



C ================================================================
C Right hand side
C ================================================================
      Integer Function Drhs(x, y, z, label, DDATA,IDATA, iSYS, F)
      Include 'fem3Dtet.fd'
C ================================================================
      Real*8  x, y, z, DDATA(*), F(*)
      Integer label, IDATA(*), iSYS(*)
C ================================================================
      iSYS(1) = 1
      iSYS(2) = 1

      F(1) = 1D0
      Drhs = TENSOR_SCALAR

      Return
      End


