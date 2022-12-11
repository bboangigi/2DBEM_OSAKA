!! ---------------------------------------------------------------------
Module FORCE_M
!! ---------------------------------------------------------------------
      Use glb2DBEM
!! ---------------------------------------------------------------------
      IMPLICIT None
!! ---------------------------------------------------------------------
Contains
!! ---------------------------------------------------------------------
SUBROUTINE FORCE(IDIF, IRAD, NB, NEQ, IPRINT, AKB, XQ, YQ, VN, ZFI, ZAB, ZEXF)

      Implicit None
!! -----------Input variables from app2DBEM_OSAKA-----------------------
      Integer(4), Intent(in)          :: IDIF, IRAD, NB, NEQ, IPRINT     ! from app2DBEM_OSAKA
      Double Precision, Intent(in) :: AKB                                 ! from app2DBEM_OSAKA
!! ---------------------------------------------------------------------
!! -----------Input variables&arrays from OFFSET_M-----------------------------
      Double precision, Allocatable, Intent(in) :: XQ(:), YQ(:)
      Double Precision, Allocatable, Intent(in) :: VN(:,:)
!! -----------Input variables&arrays from SOLVE_M-----------------------------
      Complex(8), Allocatable, Intent(in) :: ZFI(:,:)
!! ----------Arrays declaration in FORCE subroutine------------------
      Integer(4) :: I, J, K
      Double Precision :: D, DX, DY, C1, C2, CHK
      Double Precision, Allocatable :: A(:,:), B(:,:), BE(:,:)
      Double Precision, Allocatable :: EAMP(:), EPHA(:)
      Complex(8), Allocatable, Intent(out) :: ZAB(:,:)
      Complex(8), Allocatable, Intent(out) :: ZEXF(:)

      Allocate(A(3,3), B(3,3), BE(3,3))
      Allocate(EAMP(3), EPHA(3))
      Allocate(ZAB(3,3))
      Allocate(ZEXF(3))

      Do I=1,3
            Do J=1,3
                  ZAB(I,J)=Z0
            End Do

            ZEXF( I)=Z0
      End Do

      DO K=1,NB
            DX=XQ(K+1)-XQ(K)
            DY=YQ(K+1)-YQ(K)
            D =DSQRT(DX*DX+DY*DY)
      
            Do I=1,3
                  Do J=1,3
                        ZAB(I,J)=ZAB(I,J)-ZFI(J,K)*VN(I,K)*D
                  End Do
                        
                  ZEXF(I )=ZEXF(I )+ZFI(4,K)*VN(I,K)*D

            End Do
      End Do
      
      Do I=1,3
            Do J=1,3
                  A (I,J)= DREAL(ZAB(I,J))
                  B (I,J)=-DIMAG(ZAB(I,J))
            End Do

            EAMP(I)=CDABS(ZEXF(I))
            EPHA(I)=DATAN2(DIMAG(ZEXF(I)),DREAL(ZEXF(I)))*180.0D0/PI
      End Do

      Write(IDIF, 650) AKB, (EAMP(I), I=1,3), (EPHA(I), I=1,3)
      Write(IRAD, 650) AKB

      Do I = 1, 3
            Write(IRAD, 650) ( A(I, J), J = 1, 3)
      End Do

      Do I = 1, 3
            Write(IRAD, 650) ( B(I, J), J = 1, 3)
      End Do

      IF(IPRINT.EQ.0) RETURN
      WRITE(6,600) NB,AKB

!!------------------------------------------------------------
      DO I=1,3
            DO J=1,3
                  IF(I.EQ.J) then
                        C1=B (I,I)
                        C2=BE(I,I)
                        CHK=DABS(C1-C2)/DABS(C1+C2)*200.0D0
                        Write(6,610) I,I,A(I,I),B(I,I)
                  ENDIF
            End Do
      End Do
      WRITE(6,*)

      DO I=1,3
            DO J=1,3
                  IF(I.NE.J) then
                        Write(6,610) I,J,A(I,J),B(I,J)             
                  ENDIF
            End Do
      End Do
      WRITE(6,630)
!------------------------------------------------------------
      Do I=1,3
            Write(6,640) I, ZEXF(I), EAMP(I), EPHA(I)
      End Do

      600 FORMAT(//5X,'++++++++ ADDED-MASS & DAMPING COEFF. ( ', &
            'NB=',I3,', K*B/2=',F8.4,' )  +++++++',//10X, &
            'I  J',8X,'ADDED-MASS',6X,'DAMPING')
      610 FORMAT(8X,'(',I2,',',I2,')',3X,E13.4,3(2X,E13.4))
      615 FORMAT(' ')
      620 FORMAT(8X,'(',I2,',',I2,')',3X,E13.4,2(2X,E13.4))
      630 FORMAT(//5X,'+++++ WAVE EXCITING FORCE +++++', &
       //17X,'PRESSURE INTEGRAL',12X,'AMP',5X,'PHASE(DEG)')
      640 FORMAT(8X,I2,2E13.4,2X,2E13.4,3X,E11.4,2X,F9.3)
      650 FORMAT(99(1pe15.6))

      END Subroutine
!! ---------------------------------------------------------------
End Module