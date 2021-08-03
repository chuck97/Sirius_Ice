c ==============================================================
c Initialization routine
c ==============================================================
      Subroutine InitializeRCB (
c ==============================================================
     &           nE, nEmax, nP, nPmax, XYP, IPE,
     &           MaxWi, iW, MaxWr, rW, listtet, tetpmax, iERR)
      implicit none
c ==============================================================
      Integer             nE, nEmax, nP, nPmax
      Double precision    XYP(3,*)
      Integer             IPE(4,*)
      Integer             tetpmax
      Integer             listtet(tetpmax+1,*)
      Integer             MaxWi, MaxWr, iERR
      Integer             iW(*) 
      Double precision    rW(*)

      Integer             iP_pbisE, iP_flag
      Integer             iP_level 
      Integer             iP_fpe, ip_edp

      Integer             iP_ire,  iP_mr, iP_mrn
      Integer             rP_mr, rP_led
      Integer             nRmax
c ==============================================================

      nRmax   = 2*nE

1     Continue

      If (MaxWi.lt.11*nEmax+2*nRmax+41) Then
          iERR=1 ! LocalRefine needs at least 11*nEmax+2*nRmax+41
          write (*,*) 'MaxWi too small'
          return
      end if

      iP_pbisE   =              1
      iP_flag    = iP_pbisE   + 3*nEmax
      iP_level   = iP_flag    + nEmax
      iP_fpe     = iP_level   + nEmax
      ip_edp     = iP_fpe     + 24
      iP_ire     = iP_edp     + 16
      iP_mr      = iP_ire     + 6*nE
      iP_mrn     = iP_mr      + nRmax

      If (MaxWr.lt.iP_mrn + nRmax) Then
          iERR=1 ! LocalRefine needs more
          write (*,*) 'MaxWr too small, needs', iP_mrn + nRmax
          return
      end if

      rP_led    =               1
  
      call InitializeMeshData (nE, nP,  XYP, IPE, iW(iP_pbisE),
     &        iW(iP_ire), iW(iP_mr), iW(iP_mrn), rW(rP_led),
     &        iW(iP_flag), iW(iP_level), iW(iP_fpe), iW(iP_edp),
     &        listtet, tetpmax, nRmax, iERR)
      If (iERR .eq. 2) Then
        nRmax = nRmax  + nE
        Goto 1
      End if  

      Return
      End
c ==============================================================
c ==============================================================
      Subroutine InitializeMeshData (
c ==============================================================
     &           nE, nP, XYP, IPE, pbisE, IRE, mr, mrn,
     &           led, flag, level, 
     &           fpe, edp, listtet, tetpmax, nRmax, iERR)
      implicit none
c ==============================================================
      Integer             nE, nP
      Double precision    XYP(3,*)
      Integer             IPE(4,*)
      Integer             pbisE(3,*)
      Integer             IRE(6,*) 
      Integer             mr(*)
      Integer             mrn(*)
      Double precision    led(*)
      Integer             flag(*)
      Integer             level(*)
      Integer             fpe(4,6)
      Integer             edp(4,4)
      Integer             tetpmax
      Integer             listtet(tetpmax+1,*)
      Integer             iERR
      
      Integer             nR, nRmax 
      Integer             i, i1, j, j1, k, m
      Integer             pt1, pt2
      Integer             tt1, tt2
      Integer             k1, k2, kk, kk1, kk2
c ==============================================================
      iERR = 0
c
c ... define fpe
c
      fpe(1,1) = 1
      fpe(2,1) = 2
      fpe(3,1) = 1
      fpe(4,1) = 2 

      fpe(1,2) = 1
      fpe(2,2) = 3
      fpe(3,2) = 1
      fpe(4,2) = 3 

      fpe(1,3) = 1
      fpe(2,3) = 4
      fpe(3,3) = 2
      fpe(4,3) = 3 

      fpe(1,4) = 2
      fpe(2,4) = 3
      fpe(3,4) = 1
      fpe(4,4) = 4 

      fpe(1,5) = 2
      fpe(2,5) = 4
      fpe(3,5) = 2
      fpe(4,5) = 4 

      fpe(1,6) = 3
      fpe(2,6) = 4
      fpe(3,6) = 3
      fpe(4,6) = 4 
c
c ... define edp 
c
      edp(1,1) = 0
      edp(2,1) = 1
      edp(3,1) = 2
      edp(4,1) = 3

      edp(1,2) = 1
      edp(2,2) = 0
      edp(3,2) = 4
      edp(4,2) = 5

      edp(1,3) = 2
      edp(2,3) = 4 
      edp(3,3) = 0
      edp(4,3) = 6

      edp(1,4) = 3
      edp(2,4) = 5
      edp(3,4) = 6
      edp(4,4) = 0
c
c ... define level
c
      Do i = 1, nE
        level(i)=0
      End do
c
c ... generate list of tet
c
      Call Initializelisttet (nE, nP, IPE, listtet, tetpmax, iERR)
c
c ... define IRE
c
      Do i = 1, nE
        Do j = 1, 6
          IRE (j,i) = 0
        End do  
      End do  

      nR = 0

      Do i = 1, nE
        Do  j = 1, 3
          Do k = j + 1, 4
            If (IRE (edp(j,k),i) .eq. 0) Then

              nR = nR + 1
              If (nR .ge. nRmax) Then
                iERR = 2
                Return
              End if  
              IRE (edp(j,k),i) = nR
              pt1 = IPE (j, i)
              pt2 = IPE (k, i)
              Call Dist3D (XYP(1, pt1), XYP(1, pt2), led(nR))       

              Do k1 = 1, listtet(1, pt1)
                tt1 = listtet(k1 + 1, pt1)
                Do k2 = 1, listtet(1, pt2)
                  tt2 = listtet(k2 + 1, pt2)
                  If (tt1 .eq. i) Go to 10
                  If (tt1 .ne. tt2) Go to 10
                  Do kk = 1, 4
                    If (IPE(kk, tt1) .eq. pt1) kk1 = kk 
                    If (IPE(kk, tt1) .eq. pt2) kk2 = kk 
                  End do

                  IRE (edp(kk1,kk2),tt1) = nR
10                Continue
                End do
              End do

            End if
          End do
        End do
      End do

      Do i = 1, nR
        mr(i)=i
      End do 

      Call DSORT_RCB(led, mr, nR, 2) 

      Do i = 1, nR
        mrn(mr(i))=i
      End do  
   
c
c ... define pbisE
c
      Do i = 1, nE
        Do j = 1, 3
          pbisE(j, i) = 0
        End do
        flag(i) = 0
      End do

      Do i = 1, nE

        m = mrn(IRE(edp(1,2),i))
        j = 1
        kk1 = 1
        kk2 = 2
        Do k1 = 1, 3
          Do k2 = k1 + 1, 4
            j1 = edp(k1,k2)
            kk = mrn(IRE(j1,i))
            If (kk .gt. m) Then
              m = kk
              j = j1
              kk1 = k1
              kk2 = k2
            End if  
          End do  
        End do
        pbisE(1, i) = j

        m = 0 
        j = 0
        Do k1 = 1, 3
          If (k1 .eq. kk2) Go to 30
          Do k2 = k1 + 1, 4
            If (k2 .eq. kk2) Go to 20
            j1 = edp(k1,k2)
            kk = mrn(IRE(j1,i))
            If (kk .gt. m) Then
              m = kk
              j = j1
            End if 
20          Continue            
          End do  
30        Continue            
        End do
        pbisE(2, i) = j

        m = 0
        j = 0
        Do k1 = 1, 3
          If (k1 .eq. kk1) Go to 50
          Do k2 = k1 + 1, 4
            If (k2 .eq. kk1) Go to 40
            j1 = edp(k1,k2)
            kk = mrn(IRE(j1,i))
            If (kk .gt. m) Then
              m = kk
              j = j1
            End if 
40          Continue            
          End do  
50        Continue            
        End do
        pbisE(3, i) = j
      End do 


      Return
      End
c ==============================================================
      subroutine Ordertriple (a,b,iERR)
c ==============================================================
c ==============================================================
      implicit none
      
      Integer      a(3), b(3), iERR
      Integer      c(3), d(2), i, j, k
c ==============================================================

       If (a(1) .eq. a(2)) Then
           iERR = 1050   !?
           call errMes(iERR, 'Ordertriple','Bag faces')
       End if    
       If (a(1) .eq. a(3)) Then
           iERR = 1050   !?
           call errMes(iERR, 'Ordertriple','Bag faces')
       End if    
       If (a(2) .eq. a(3)) Then
           iERR = 1050   !?
           call errMes(iERR, 'Ordertriple','Bag faces')
       End if 

       Do i = 1, 3 
         c(i) = a(i)
       End do

       k = 1
       Do i = 2, 3
         If (c(i) .lt. c(k)) k = i
       End do

       j = 0
       Do i = 1, 3
         If (i .ne. k) Then
             j = j + 1
             d(j) = i
         End if
       End do

       b(1) = c(k)
       If (c(d(1)) .lt. c(d(2))) Then
          b(2) = c(d(1))  
          b(3) = c(d(2))
        Else
          b(2) = c(d(2))  
          b(3) = c(d(1))
       End if

      Return
      End
c ==============================================================
c Addition auxiliary procedure.
c ==============================================================
      subroutine OrientandMakeNeighbors (
c ==============================================================
     &           nP, nPmax, nF, nFmax, nE, nEmax, 
     &           IPE, IEF, IPF, labelF, listtet, tetpmax,
     &           iERR)
      implicit none
c ==============================================================
      Integer             nP, nPmax, nF, nFmax, nE, nEmax
      Integer             IPE(4,*),IEF(4,*),IPF(3,*)
      Integer             labelF(*)
      Integer             tetpmax
      Integer             listtet((tetpmax + 1), *)
      Integer             iERR

      Integer             pf(3)
      Integer             i,j,i1,i2,i3
      Integer             k, m, maxlabelF
c ==============================================================
c
c  Check of the orientation of boundary faces 
c

       iERR = 0
       Do i = 1, nF
         call Ordertriple (IPF(1,i),IPF(1,i),iERR)
       end do 
   
c
c ... generate list of tetetrahedra
c
       Call Initializelisttet (nE, nP, IPE, listtet, tetpmax, iERR)

c        
c ... form IEF
c
       maxlabelF=0
       Do i = 1, nF
         j = labelF(i)
         If (j .gt. maxlabelF) maxlabelF = j
       End do
       maxlabelF = maxlabelF +1

       Do i = 1, nE
         Do j = 1, 4
           IEF(j,i)=-maxlabelF
         End do
       End do

       Do i = 1, nE
         Do j = 1, 4
           If (IEF(j,i) .eq. -maxlabelF) Then

             m = 0
             Do k = 1, 4
               If (k .ne. j) Then
                  m = m + 1
                  pf(m) =  IPE(k,i)
               End if
             End do

             Do i1 = 1, listtet(1, pf(1))
               If (listtet(i1 + 1, pf(1)) .ne. i) Then
                 Do i2 = 1, listtet(1, pf(2))
                   If (listtet(i2 + 1, pf(2)) .ne. i) Then
                     Do i3 = 1, listtet(1, pf(3))
                       If ((listtet(i1 + 1, pf(1)) .eq.
     &                      listtet(i2 + 1, pf(2))) .and.
     &                     (listtet(i1 + 1, pf(1)) .eq.
     &                      listtet(i3 + 1, pf(3)))) Then
                             IEF(j,i) = listtet(i1 + 1, pf(1))
                             Go to 10
                       End if                      
                     End do
                   End if    
                 End do
               End if
             End do
           
             call Ordertriple (pf,pf,iERR)
             
             Do k = 1, nF
               If ((pf(1) .eq. IPF(1,k)) .and.
     &             (pf(2) .eq. IPF(2,k)) .and.
     &             (pf(3) .eq. IPF(3,k))) Then 
                    IEF(j,i) = - labelF(k)
                    Go to 10
               End if
             End do

             write (*,*) 'i=',i, ' versh=',j
             write (*,*) ' IPE=',IPE(1,i),IPE(2,i),IPE(3,i),IPE(4,i)
             write (*,*) ' IEF=',IEF(1,i),IEF(2,i),IEF(3,i),IEF(4,i)
             call errMes(iERR, 'Init_mesh_gen',
     &                'Mistake in the bound date')
           End if
10         Continue

         End do
       End do 
c


      Return
      End
      
c ==============================================================
c Recover boundary edges.
c ==============================================================
      Subroutine RestoreBndEdges (
     &           nF, nFmax, nE, 
     &           IPE, IEF, IPF, labelF, iERR)
      implicit none
c ==============================================================
      Integer             nF, nFmax, nE
      Integer             IPE(4,*), IEF(4,*)
      Integer             IPF(3,*), labelF(*),  iERR

      Integer             i,j,k,m,m1
c ==============================================================

      nF = 0

      Do i = 1, nE
        Do j =1, 4
        
          k = IEF(j,i)
          If (k .le. 0) Then
            nF = nF + 1
            If (nF .gt. nFmax) Then 
             iERR = 1004
             call errMes(iERR, 'Post_procesing','nFmax is too small')
            End if
            
            m1 = 0
            Do m = 1, 4
              If (m .ne. j) Then
                 m1 = m1 + 1
                 IPF (m1, nF) = IPE(m, i)
              End if
            End do
            labelF(nF) = -k            
            call Ordertriple (IPF(1,nF),IPF(1,nF),iERR)
          End if        
          
        End do
      End do
        
      Return
      End
c ==============================================================
c Codering the number of the point tet
c ==============================================================
      Subroutine Code (k, fl1, fl2)
c ==============================================================
      implicit none
      
      Integer        k
      Logical        fl1, fl2
c ==============================================================
      
      If (k .ge. 3) Then
        fl1 = .true.
       Else
        fl1 = .false.
      End if

      If ((k .eq. 2) .or. (k .eq. 4)) Then
        fl2 = .true.
       Else
        fl2 = .false.
      End if

      Return
      End
c      
c ==============================================================
c Uncodering the number of the point tet
c ==============================================================
      Subroutine Uncode (fl1, fl2, k)
c ==============================================================
      implicit none
      
      Logical        fl1, fl2
      Integer        k
c ==============================================================
     
      If ((.not. fl1) .and. (.not. fl2)) k = 1
      If ((.not. fl1) .and. (      fl2)) k = 2
      If ((      fl1) .and. (.not. fl2)) k = 3
      If ((      fl1) .and. (      fl2)) k = 4

      Return
      End
c      
c ==============================================================
c ==============================================================
      Subroutine Initializelisttet (
c ==============================================================
     &           nE, nP, IPE, listtet, tetpmax, iERR)
      implicit none
c ==============================================================
      Integer             nE, nP
      Integer             IPE(4,*)
      Integer             tetpmax
      Integer             listtet(tetpmax+1,*)
      
      Integer             i, j, k, m, iERR
c ==============================================================

       Do i = 1, nP
         listtet(1, i) = 0
       End do

       Do i = 1, nE
         Do j = 1, 4
           k = IPE(j, i)
           listtet(1, k) = listtet(1, k) + 1
           If (listtet(1, k) .gt. tetpmax) Then
             iERR = 1007 
             call errMes(iERR, 'Initializelisttet',
     &                         'tetpmax is too small')
           End if  
           m = listtet(1, k)
           listtet(m + 1, k) = i
         End do
       End do

      Return
      End
c ==============================================================
      subroutine Dist3D (pt1,pt2,res)
c ==============================================================
c ==============================================================
      implicit none
      
      Double precision  pt1(3), pt2(3), res
c ==============================================================

       res = dsqrt((pt1(1) - pt2(1)) **2 +
     &             (pt1(2) - pt2(2)) **2 +
     &             (pt1(3) - pt2(3)) **2)
                         
      Return
      End
      

*DECK DSORT
      SUBROUTINE DSORT_RCB(DX, DY, N, KFLAG)
C***BEGIN PROLOGUE  DSORT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2B
C***TYPE      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   DSORT sorts array DX and optionally makes the same interchanges in
C   array DY.  The array DX may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      DX - array of values to be sorted   (usually abscissas)
C      DY - array to be (optionally) carried along
C      N  - number of values in array DX to be sorted
C      KFLAG - control parameter
C            =  2  means sort DX in increasing order and carry DY along.
C            =  1  means sort DX in increasing order (ignoring DY)
C            = -1  means sort DX in decreasing order (ignoring DY)
C            = -2  means sort DX in decreasing order and carry DY along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***REVISION HISTORY  (YYMMDD)
C   761101  DATE WRITTEN
C   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced statement labels.  (WRB)
C   891024  Changed category.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901012  Declared all variables; changed X,Y to DX,DY; changed
C           code to parallel SSORT. (M. McClain)
C   920501  Reformatted the REFERENCES section.  (DWL, WRB)
C   920519  Clarified error messages.  (DWL)
C   920801  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
C***END PROLOGUE  DSORT
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
      INTEGER          DY(*)
C     .. Local Scalars ..
      DOUBLE PRECISION R, T, TT
      INTEGER          TTY, TY
      INTEGER I, IJ, J, K, KK, L, M, NN
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. External Subroutines ..
      EXTERNAL XERMSG
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***FIRST EXECUTABLE STATEMENT  DSORT
      NN = N
      IF (NN .LT. 1) THEN
         WRITE(*,*)'DSORT_RCB: number of values to be sorted < 1'
         RETURN
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         WRITE(*,*)'DSORT_RCB: control parameter K is not 2,1,-1,-2'
         RETURN
      ENDIF
C
C     Alter array DX to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            DX(I) = -DX(I)
   10    CONTINUE
      ENDIF
C
      IF (KK .EQ. 2) GO TO 100
C
C     Sort DX only
C
      M = 1
      I = 1
      J = NN
      R = 0.375D0
C
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
C
   30 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = DX(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (DX(I) .GT. T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than than T, interchange with T
C
      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (DX(I) .GT. T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
   40 L = L-1
      IF (DX(L) .GT. T) GO TO 40
C
C     Find an element in the first half of the array which is greater
C     than T
C
   50 K = K+1
      IF (DX(K) .LT. T) GO TO 50
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         GO TO 40
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
C
C     Begin again on another portion of the unsorted array
C
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
C
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = DX(I+1)
      IF (DX(I) .LE. T) GO TO 80
      K = I
C
   90 DX(K+1) = DX(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 90
      DX(K+1) = T
      GO TO 80
C
C     Sort DX and carry DY along
C
  100 M = 1
      I = 1
      J = NN
      R = 0.375D0
C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
C
  120 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = DX(IJ)
      TY = DY(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (DX(I) .GT. T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
         DY(IJ) = DY(I)
         DY(I) = TY
         TY = DY(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than T, interchange with T
C
      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)
         DY(IJ) = DY(J)
         DY(J) = TY
         TY = DY(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (DX(I) .GT. T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
            DY(IJ) = DY(I)
            DY(I) = TY
            TY = DY(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 L = L-1
      IF (DX(L) .GT. T) GO TO 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 K = K+1
      IF (DX(K) .LT. T) GO TO 140
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         TTY = DY(L)
         DY(L) = DY(K)
         DY(K) = TTY
         GO TO 130
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
C
C     Begin again on another portion of the unsorted array
C
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
C
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = DX(I+1)
      TY = DY(I+1)
      IF (DX(I) .LE. T) GO TO 170
      K = I
C
  180 DX(K+1) = DX(K)
      DY(K+1) = DY(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 180
      DX(K+1) = T
      DY(K+1) = TY
      GO TO 170
C
C     Clean up
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            DX(I) = -DX(I)
  200    CONTINUE
      ENDIF
      RETURN
      END
