!! ---------------------------------------------------------------------
Module MOTION_M
!! ---------------------------------------------------------------------
    Use glb2DBEM
!! ---------------------------------------------------------------------
      Implicit None
!! ---------------------------------------------------------------------
Contains
!! ---------------------------------------------------------------------
      SUBROUTINE MOTION(IRAO, IPRINT, OGD, KZZB, SML, AKB, C22, CMAS, GM, ZAB, ZEXF)

      Implicit None
      !! -----------Input variables from app2DBEM_OSAKA-----------------------
      Integer(4), Intent(in)          :: IRAO, IPRINT     ! from app2DBEM_OSAKA
      Double Precision, Intent(in) :: OGD, KZZB, SML, AKB ! from app2DBEM_OSAKA
      !! -----------Input variables from OFFSET_M-----------------------------
      Double Precision, Intent(in) :: C22, CMAS, GM     ! from OFFSET_M
      !! -----------Input variables from FORCE_M-----------------------------
      Complex(8), Allocatable, Intent(in) :: ZAB(:,:), ZEXF(:)
      !! ----------Arrays declaration in MOTION subroutine------------------     
      Integer(4), Parameter :: NDIM = 3, N = 3, NEP = 1
      Integer(4) :: I, J, K
      Double Precision, Allocatable :: AMPG(:), PHAG(:)
      Complex(8), Allocatable :: ZAA_local(:,:), ZBB_local(:,:)
      Complex(8), Allocatable :: ZMTNG(:), ZMTNO(:)

      Allocate(AMPG(3), PHAG(3)) 
      Allocate(ZAA_local(3,3), ZBB_local(3,1))
      Allocate(ZMTNG(3), ZMTNO(3))

      ZAA_local(1,1)=-AKB*(CMAS+ZAB(1,1))
      ZAA_local(1,2)=-AKB* ZAB(1,2)
      ZAA_local(1,3)=-AKB*(ZAB(1,3)+OGD*ZAB(1,1))
      ZBB_local(1,1)= ZEXF(1)

      ZAA_local(2,1)=-AKB*ZAB(2,1)
      ZAA_local(2,2)=-AKB*(CMAS+ZAB(2,2))+C22
      ZAA_local(2,3)=-AKB*(ZAB(2,3)+OGD*ZAB(2,1))
      ZBB_local(2,1)= ZEXF(2)

      ZAA_local(3,1)=-AKB*(ZAB(3,1)+OGD*ZAB(1,1))
      ZAA_local(3,2)=-AKB*(ZAB(3,2)+OGD*ZAB(1,2))
      ZAA_local(3,3)=-AKB*(CMAS*KZZB**2+ZAB(3,3)+OGD*ZAB(1,3)+OGD*(ZAB(3,1)+OGD*ZAB(1,1)))+CMAS*GM
      ZBB_local(3,1)= ZEXF(3)+OGD*ZEXF(1)

      CALL ZSWEEP(NDIM,N,ZAA_local,ZBB_local,NEP,SML)
      ! NDIM = 3, N=3, NEP = 1

      IF(CDABS(ZAA_local(1,1)).LT.SML) then
        Write(6,600)
      ENDIF
      
  600 FORMAT(///10X,'+++ ERROR: ZSWEEP IN (MOTION) +++'///)

      DO I=1,3
          ZMTNG(I)=ZBB_local(I,1)
      End Do

      ZMTNO(1)=ZMTNG(1)+OGD*ZMTNG(3)
      ZMTNO(2)=ZMTNG(2)
      ZMTNO(3)=ZMTNG(3)

      DO I=1,3
          AMPG(I)=CDABS(ZMTNG(I))
          IF(I.EQ.3) AMPG(I)=AMPG(I)/AKB
          PHAG(I)=DATAN2(DIMAG(ZMTNG(I)),DREAL(ZMTNG(I)))*180.0D0/PI
      End Do

      Write(IRAO, 620) AKB, AMPG(1), AMPG(2), AMPG(3), PHAG(1), PHAG(2), PHAG(3)   ! Textfile에 결과값 기술
      Write(   6, 610) AKB, AMPG(1), PHAG(1), AMPG(2), PHAG(2), AMPG(3), PHAG(3)   ! Command창에 결과값 기술

  610 FORMAT(//5X,'+++++ MOTIONS ABOUT ''G'' FOR K*B/2=',F7.3, '+++++', &
       //21X,'AMP.',7X,'PHASE',                                         &
          /9X, 'SWAY  ',E11.4, 2X, F9.3,' (DEG)',                        &
          /9X, 'HEAVE ',E11.4, 2X, F9.3,' (OEG)',                        &
          /9X, 'ROLL  ',E11.4, 2X, F9.3,' (DEG)')

  620 FORMAT(99(1pe15.6))

      END Subroutine
!! ---------------------------------------------------------------
!! !! ---------------------------------------------------------------------
SUBROUTINE ZSWEEP(NDIM,N,ZAA_local,ZBB_local,NEP,EPS)

      Implicit None
      Integer(4) :: I, J, K, N, NDIM, IP, NEP 
      Double Precision :: EPS, P 

      Complex(8) :: ZW
      Complex(8), Allocatable :: ZAA_local(:,:), ZBB_local(:,:)

      DO 5 K=1,N
      P=0.0D0

            DO 1 I=K,N
            IF(P.GE.CDABS(ZAA_local(I,K))) GOTO 1
            P=CDABS(ZAA_local(I,K))
            IP=I
          1 CONTINUE
            IF(P.LE.EPS) GOTO 6
            IF(IP.EQ.K)  GOTO 7
              DO 2 J=K,N
              ZW=ZAA_local(K,J)
              ZAA_local(K,J)=ZAA_local(IP,J)
          2   ZAA_local(IP,J)=ZW
                DO 20 J=1,NEP
                ZW=ZBB_local(K,J)
                ZBB_local(K,J)=ZBB_local(IP,J)
        20     ZBB_local(IP,J)=ZW
          7 CONTINUE

      IF(K.EQ.N) GOTO 70
      DO 3 J=K+1,N
    3 ZAA_local(K,J)=ZAA_local(K,J)/ZAA_local(K,K)
    70   DO 30 J=1,NEP
    30   ZBB_local(K,J)=ZBB_local(K,J)/ZAA_local(K,K)
      DO 5 I=1,N
      IF(I.EQ.K) GOTO 5
      IF(K.EQ.N) GOTO 40
        DO 4 J=K+1,N
    4   ZAA_local(I,J)=ZAA_local(I,J)-ZAA_local(I,K)*ZAA_local(K,J)
    40   CONTINUE
        DO 45 J=1,NEP
    45   ZBB_local(I,J)=ZBB_local(I,J)-ZAA_local(I,K)*ZBB_local(K,J)
    5 CONTINUE
      ZAA_local(1,1)=(1.0D0,0.0D0)
      RETURN
    6 ZAA_local(1,1)=DCMPLX(DABS(P),0.0D0)
      RETURN
      END Subroutine
!!---------------------------------------------------------------------
End Module
