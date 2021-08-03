C =====================================================================
      Subroutine GMVmesh(fname, Chan,
     &                   Nvrt,vrt, Ntet,tet, Nbnd,bnd,lbbnd)
C =====================================================================
C Routine saves mesh in a GMV file.
C =====================================================================
      Implicit None

      Real*8        vrt(3,*)
      Integer       tet(4,*), bnd(3,*), lbbnd(*)
      Integer       Nvrt, Ntet, Nbnd, Chan
      Character*(*) fName

      Real*8        DETofTET, det
      Integer       i, j, i1, i2, i3, i4
      Character*70  fNameExt

C =====================================================================
      fNameExt = fName
      Open(Chan, file=fNameExt, status='UNKNOWN')

      Write(Chan,'(A)') 'gmvinput ascii'
      Write(Chan,'(A)') 'codename aniVIEW'

      Write(Chan,*) 'nodev ',nvrt
      Do i = 1, nvrt
         Write(Chan,'(3f14.6)')(vrt(j,i),j=1,3)
      End do

      Write(Chan,*) 'cells ',ntet

      Do i = 1, ntet
         i1 = tet(1,i)
         i2 = tet(2,i)
         i3 = tet(3,i)
         i4 = tet(4,i)

         det = DETofTET(vrt(1,i1), vrt(1,i2), vrt(2,i3), vrt(2,i4))

         Write(Chan,*) ' tet 4'
         If(det.GT.0) Then
            Write(Chan,*) i1, i2, i4, i3
         Else
            Write(Chan,*) i1, i2, i3, i4
         End if
      End do

      Write(Chan,*)
      Write(Chan,'(A)') 'polygons '

      Do i = 1, nbnd
         Write(Chan,'(2I6,3F10.6)') lbbnd(i),3,(vrt(1,bnd(j,i)),j=1,3)
         Write(Chan,'(3F10.6)') (vrt(2,bnd(j,i)),j=1,3)
         Write(Chan,'(3F10.6)') (vrt(3,bnd(j,i)),j=1,3)
      End do

      Write(Chan,'(A)') 'endpoly '
      Write(Chan,*)
      Write(Chan,'(A)') 'endgmv'

      Close(Chan)

      Return
      End



C =====================================================================
      Subroutine GMVscalarVrt(U, fname, Chan,
     &                        Nvrt,vrt, Ntet,tet, Nbnd,bnd,lbbnd)
C =====================================================================
C Routine saves mesh and nodal solution U in a GMV file.
C =====================================================================
      Implicit None

      Real*8        U(*), vrt(3,*)
      Integer       tet(4,*), bnd(3,*), lbbnd(*)
      Integer       Nvrt, Ntet, Nbnd, Chan
      Character*(*) fName

      Real*8        DETofTET, det
      Integer       i, j, i1, i2, i3, i4
      Character*70  fNameExt

C =====================================================================
      fNameExt = fName
      Open(Chan, file=fNameExt, status='UNKNOWN')

      Write(Chan,'(A)') 'gmvinput ascii'
      Write(Chan,'(A)') 'codename aniVIEW'

      Write(Chan,*) 'nodev ', Nvrt
      Do i = 1, nvrt
         Write(Chan,'(3f14.6)') (vrt(j,i),j=1,3)
      End do

      Write(Chan,*) 'cells ', Ntet

      Do i = 1, Ntet
         i1 = tet(1,i)
         i2 = tet(2,i)
         i3 = tet(3,i)
         i4 = tet(4,i)

         det = DETofTET(vrt(1,i1), vrt(1,i2), vrt(2,i3), vrt(2,i4))

         Write(Chan,*) ' tet 4'
         If(det.GT.0) Then
            Write(Chan,*) i1, i2, i4, i3
         Else
            Write(Chan,*) i1, i2, i3, i4
         End if
      End do

      Write(Chan,*)
      Write(Chan,'(A)') 'polygons '

      Do i = 1, nbnd
         Write(Chan,'(2I6,3F10.6)') lbbnd(i),3,(vrt(1,bnd(j,i)),j=1,3)
         Write(Chan,'(3F10.6)') (vrt(2,bnd(j,i)),j=1,3)
         Write(Chan,'(3F10.6)') (vrt(3,bnd(j,i)),j=1,3)
      End do

      Write(Chan,'(A)') 'endpoly '
      Write(Chan,*)
      Write(Chan,'(A)') 'variable'
      Write(Chan,'(A)') 'solution 1'
      Do i = 1, Nvrt
         Write(Chan,'(F14.5)') U(i)
      End do

      Write(Chan,'(A)') 'endvars'
      Write(Chan,'(A)') 'endgmv'

      Close(Chan)

      Return
      End



C =====================================================================
      Subroutine GMVscalarTet(U, fname, Chan,
     &                        Nvrt,vrt, Ntet,tet, Nbnd,bnd,lbbnd)
C =====================================================================
C Routine saves mesh and cellwise solution U in a GMV file.
C =====================================================================
      Implicit None

      Real*8        U(*), vrt(3,*)
      Integer       tet(4,*), bnd(3,*), lbbnd(*)
      Integer       Nvrt, Ntet, Nbnd, Chan
      Character*(*) fName

      Real*8        DETofTET, det
      Integer       i, j, i1, i2, i3, i4
      Character*70  fNameExt

C =====================================================================
      fNameExt = fName
      Open(Chan, file=fNameExt, status='UNKNOWN')

      Write(Chan,'(A)') 'gmvinput ascii'
      Write(Chan,'(A)') 'codename aniVIEW'

      Write(Chan,*) 'nodev ', Nvrt
      Do i = 1, nvrt
         Write(Chan,'(3f14.6)') (vrt(j,i),j=1,3)
      End do

      Write(Chan,*) 'cells ', Ntet

      Do i = 1, ntet
         i1 = tet(1,i)
         i2 = tet(2,i)
         i3 = tet(3,i)
         i4 = tet(4,i)

         det = DETofTET(vrt(1,i1), vrt(1,i2), vrt(2,i3), vrt(2,i4))

         Write(Chan,*) ' tet 4'
         If(det.GT.0) Then
            Write(Chan,*) i1, i2, i4, i3
         Else
            Write(Chan,*) i1, i2, i3, i4
         End if
      End do

      Write(Chan,*)
      Write(Chan,'(A)') 'polygons '

      Do i = 1, Nbnd
         Write(Chan,'(2I6,3F10.6)') lbbnd(i),3,(vrt(1,bnd(j,i)),j=1,3)
         Write(Chan,'(3F10.6)') (vrt(2,bnd(j,i)),j=1,3)
         Write(Chan,'(3F10.6)') (vrt(3,bnd(j,i)),j=1,3)
      End do

      Write(Chan,'(A)') 'endpoly '
      Write(Chan,*)
      Write(Chan,'(A)') 'variable'
      Write(Chan,'(A)') 'solution 0'
      Do i = 1, Ntet
         Write(Chan,'(F14.5)') U(i)
      End do

      Write(Chan,'(A)') 'endvars'
      Write(Chan,'(A)') 'endgmv'

      Close(Chan)

      Return
      End



C =====================================================================
      Subroutine GMVvectorTet(U, fname, Chan,
     &                        Nvrt,vrt, Ntet,tet, Nbnd,bnd,lbbnd)
C =====================================================================
C Routine saves mesh and cellwise vector solution U in a GMV file.
C =====================================================================
      Implicit None

      Real*8        U(3,*), vrt(3,*)
      Integer       tet(4,*), bnd(3,*), lbbnd(*)
      Integer       Nvrt, Ntet, Nbnd, Chan
      Character*(*) fName

      Real*8        DETofTET, det
      Integer       i, j, i1, i2, i3, i4
      Character*70  fNameExt

C =====================================================================
      fNameExt = fName
      Open(Chan, file=fNameExt, status='UNKNOWN')

      Write(Chan,'(A)') 'gmvinput ascii'
      Write(Chan,'(A)') 'codename aniVIEW'

      Write(Chan,*) 'nodev ', Nvrt
      Do i = 1, nvrt
         Write(Chan,'(3f14.6)') (vrt(j,i),j=1,3)
      End do

      Write(Chan,*) 'cells ', Ntet

      Do i = 1,Ntet
         i1 = tet(1,i)
         i2 = tet(2,i)
         i3 = tet(3,i)
         i4 = tet(4,i)

         det = DETofTET(vrt(1,i1), vrt(1,i2), vrt(2,i3), vrt(2,i4))

         Write(Chan,*) ' tet 4'
         If(det.GT.0) Then
            Write(Chan,*) i1, i2, i4, i3
         Else
            Write(Chan,*) i1, i2, i3, i4
         End if
      End do

      Write(Chan,*)
      Write(Chan,'(A)') 'polygons '

      Do i = 1, nbnd
         Write(Chan,'(2I6,3F10.6)') lbbnd(i),3,(vrt(1,bnd(j,i)),j=1,3)
         Write(Chan,'(3F10.6)') (vrt(2,bnd(j,i)),j=1,3)
         Write(Chan,'(3F10.6)') (vrt(3,bnd(j,i)),j=1,3)
      End do

      Write(Chan,'(A)') 'endpoly '
      Write(Chan,*)
      Write(Chan,'(A)') 'velocity 0'
      Write(Chan,'(10000000F14.5)') (U(1,i),i=1,Ntet)
      Write(Chan,'(10000000F14.5)') (U(2,i),i=1,Ntet)
      Write(Chan,'(10000000F14.5)') (U(3,i),i=1,Ntet)

      Write(Chan,'(A)') 'endgmv'

      Close(Chan)

      Return
      End


      
C =====================================================================
      Real*8 Function DETofTET(v1, v2, v3, v4)
C =====================================================================
      Real*8  v1(3), v2(3), v3(3), v4(3)
    
      Real*8  dx12,dy12,dz12, dx13,dy13,dz13, dx14,dy14,dz14
C =====================================================================

      dx12 = v2(1) - v1(1)
      dy12 = v2(2) - v1(2)
      dz12 = v2(3) - v1(3)
      dx13 = v3(1) - v1(1)
      dy13 = v3(2) - v1(2)
      dz13 = v3(3) - v1(3)
      dx14 = v4(1) - v1(1)
      dy14 = v4(2) - v1(2)
      dz14 = v4(3) - v1(3)

      DETofTET = dx12*dy13*dz14 + dy12*dz13*dx14 + dz12*dx13*dy14
     &         - dx14*dy13*dz12 - dy14*dz13*dx12 - dz14*dx13*dy12

      Return
      End

