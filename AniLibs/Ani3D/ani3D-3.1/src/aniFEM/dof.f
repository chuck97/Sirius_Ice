C ======================================================================
      Subroutine listDOF(FEMtype, N, dof)
C ======================================================================
      implicit none
      include 'fem3Dtet.fd'
 
      Integer FEMtype, N, dof(*), i

      If(FEMtype.EQ.FEM_P0) Then
         N = 1
         dof(1) = Edof

      Else If(FEMtype.EQ.FEM_P1) Then
         N = 4
         Do i = 1, 4
            dof(i) = Vdof
         End do

      Else If(FEMtype.EQ.FEM_P2) Then
         N = 10
         Do i = 1, 4
            dof(i)   = Vdof
         End do
         Do i = 1, 6
            dof(i+4) = Rdof
         End do

      Else If(FEMtype.EQ.FEM_P3) Then
         N = 20
         Do i = 1, 4
            dof(i)     = Vdof
         End do
         Do i = 1, 6
            dof(i+4)   = Rdof
            dof(i+10)  = Rdof
         End do
         Do i = 1, 4
            dof(i+16)  = Fdof
         End do
       
      Else If(FEMtype.EQ.FEM_P1vector) Then
         N = 12
         Do i = 1, 4
            dof(i)      = Vdof
            dof(i + 4)  = Vdof
            dof(i + 8)  = Vdof
         End do

      Else If(FEMtype.EQ.FEM_P2vector) Then
         N = 30
         Do i = 1, 4
            dof(i)    = Vdof
            dof(i+10) = Vdof
            dof(i+20) = Vdof
         End do
         Do i = 1, 6
            dof(i+4)  = Rdof
            dof(i+14) = Rdof
            dof(i+24) = Rdof
         End do

      Else If(FEMtype.EQ.FEM_RT0) Then
         N = 4
         Do i = 1, 4
           dof(i) = FdofOrient
         End do

      Else If(FEMtype.EQ.FEM_ND0) Then
         N = 6
         Do i = 1, 6
           dof(i) = RdofOrient
         End do

      Else If(FEMtype.EQ.FEM_CR1) Then
         N = 3
         Do i = 1, 3
            dof(i) = Rdof
         End do

      Else
         N = 0
      End if

      Return 
      End


C ======================================================================
      Integer Function orderDOF(FEMtype)
C ======================================================================
C Calculate maximal order for integration
C ======================================================================
      implicit none
      include 'fem3Dtet.fd'
 
      Integer FEMtype

      If(FEMtype.EQ.FEM_P0) Then
         orderDOF = 1

      Else If(FEMtype.EQ.FEM_P1) Then
         orderDOF = 1

      Else If(FEMtype.EQ.FEM_P2) Then
         orderDOF = 2

      Else If(FEMtype.EQ.FEM_P3) Then
         orderDOF = 3
       
      Else If(FEMtype.EQ.FEM_P1vector) Then
         orderDOF = 1

      Else If(FEMtype.EQ.FEM_RT0) Then
         orderDOF = 1

      Else If(FEMtype.EQ.FEM_ND0) Then
         orderDOF = 2

      Else If(FEMtype.EQ.FEM_CR1) Then
         orderDOF = 1

      Else
         orderDOF = 1
      End if

      Return 
      End


