!! ---------------------------------------------------------------------
Module OFFSET_M
!! ---------------------------------------------------------------------
      Use glb2DBEM
!! ---------------------------------------------------------------------
      Implicit None
!! ---------------------------------------------------------------------
Contains
!! ---------------------------------------------------------------------
!! ---------------------------------------------------------------------
SUBROUTINE OFFSET(NB, NT, H0, SIG1, SIG2, OGD, KZZB, IPRINT, XP, YP, XQ, YQ, VN, C22, CMAS, GM)
      Implicit None
      Integer(4), Intent(in)          :: NB, NT, IPRINT
      Double Precision, Intent(in) :: H0, SIG1, SIG2, OGD, KZZB
!! ---------------------------------------------------------------------
      !! NB: 단면 패널 수 (80)                                     (1)
      !! NT: Irregular freq.현상 방지를 위한 외부 소스 3개 추가(83) (1)
      !! IPRINT:                                                  (0)
      !! H0: Half-Breadth(0.5B) to Draft(D): B/(2D) (1.0)         (1)
      !! Sigma: Sig1(0.8) / Sig2(0.8)                             (0)
      !! OGD: 무게 중심 / 흘수, OG / D: 0.05                       (1)
      !! KZZB: Gyration of radius /Half-breadth (KZZ/(0.5B))      (0)

      Integer(4) :: IOFF1, IOFF2, IAD, I, J, II, K
      !! IOFF 텍스트 생성 관련 호출 번호                           (1)
      !! IAD NT - NB
      Double Precision :: DTH, Sigma, RSUB, AMD, A1, A3, AMB, TH, D, DS, DX, DY, C22, CMAS, KZZ, OG, SUM, S1, S2, S3, OBM, GM

      !! DTH: 180도 단면을 Panel 개수로 나눈 값                    (1)
      !! RSUB:                                                    (0)
      !! AMD:                                                     (0)
      !! A1:                                                      (0)
      !! A3:                                                      (0)
      !! AMB:                                                     (0)
      !! TH:                                                      (1)
      !! D: 현 Panel 다음 Panel 중점 간 거리 sqrt(DX^2+DY2)        (1)
      !! DS:
      !! DX: 현 Panel 다음 Panel 중점 간 X 축 거리                 (1)
      !! DY: 현 Panel 다음 Panel 중점 간 Y 축 거리                 (1)
      !! C22:
      !! CMAS:
      !! KZZ:
      !! OG:
      !! KZZ:
      !! SUM:
      !! S1:
      !! S2:
      !! S3:
      !! OBM:
      !! GM:
      Double precision, Allocatable :: LIST(:), XP(:), YP(:), XQ(:), YQ(:)
      Double Precision, Allocatable :: VN(:,:)
      !! 공간 할당을 최소 3개 이상은 해야, Solver가 정상으로 작동함!
      !! NB+1이나 NB가 아니라 왜 NB+3을 해야 하는지 모르겠음.
      !! 혹시 Irregular frequency 때문에 위치하는 3개의 Source 때문에?
      Allocate(LIST(NB+3), XP(NB+3), YP(NB+3), XQ(NB+3), YQ(NB+3))
      Allocate(VN(3,NB+3))
      !! XP: 단면 곡선 각 Panel의 중점 X 방향 위치                (1)
      !! YP: 단면 곡선 각 Panel의 중점 Y 방향 위치                (1)
      !! XQ: 단면 곡선 각 Panel의 Node 점 X 방향 위치             (1)
      !! YQ: 단면 곡선 각 Panel의 Node 점 Y 방향 위치             (1)
      !! VN: Generalized Normal Vector 1:Sway 2:Panel 2: Cos 값  (0)

      Open(newunit = IOFF1, file='Offset_default.dat',status='old', Action='Read')

      Do K = 1, NB+1
            Read(IOFF1,*) LIST(K), XQ(K),YQ(K)
      End Do
      close(IOFF1)

      IAD=NT-NB
      DTH=PI/DFLOAT(NB)

      Do I=1, NB
            XP(I)=(XQ(I+1)+XQ(I))/2.0D0
            YP(I)=(YQ(I+1)+YQ(I))/2.0D0
            DX=XQ(I+1)-XQ(I)
            DY=YQ(I+1)-YQ(I)
            D =DSQRT(DX*DX+DY*DY)
            VN(1,I)= DY/D
            VN(2,I)=-DX/D
            VN(3,I)=XP(I)*VN(2,I)-YP(I)*VN(1,I)
      End Do

      IF(IAD.EQ.0) GOTO 130
      DS=(XQ(1)-XQ(NB+1))/DFLOAT(IAD+1)

      !! 
      DO I=1,IAD
            II=NB+I
            XP(II)=XQ(NB+1)+DS*DFLOAT(I)
            YP(II)=0.0D0
      End Do

      130 CMAS=(SIG1+SIG2)/H0
      C22 =(XQ(1)-XQ(NB+1))/XQ(1)
      OG  =OGD/H0
      KZZ =KZZB
      SUM=0.0D0

      !! 면적 관련 계산값(?) 사다리꼴 정리(?)
      DO J=1,NB
            S1 =YQ(J+1)-YQ(J)
            S2 =XQ(J  )*(2.0D0*YQ(J  )+YQ(J+1))
            S3 =XQ(J+1)*(2.0D0*YQ(J+1)+YQ(J  ))
      SUM=SUM+S1*(S2+S3)
      End Do

      !! OBM과 GM 계산 때문에 SUM, S1, S2, S3 계산했음!       (0)
      OBM=SUM/6.0D0
      GM =(2.0D0/3.0D0-OBM)/CMAS+OG

      Open(newunit = IOFF2, file='Offset.dat',status='replace', Action='Write')

      Do J = 1,NB+1
            write(IOFF2, "(i5,99(1pe15.6))") J, XQ(J), YQ(J)
      End do
      close(IOFF2)

      WRITE(6,600) CMAS,C22,OGD,KZZ,GM
      IF(IPRINT.EQ.0) RETURN
      WRITE(6,610)
      Do 300 J=1,NB+1
      300 WRITE(6,620) J,XQ(J),YQ(J),XP(J),YP(J)
      600 FORMAT( &
            15X,'NONDIMENSIONAL MASS------- S/(B/2)**2=',F8.5,/ &
            /15X,'HEAVE RESTORING FORCE COEFF--AW/(B/2)=',F8.5,/ &
            /15X,'CENTER OF GRAVITY----------------OG/D=',F8.5,/ &
            /15X,'GYRATIONAL RADIUS-----------KZZ/(B/2)=',F8.5,/ &
            /15X,'METACENTRIC HEIGHT-----------GM/(B/2)=',F8.5/)
      610 FORMAT(/15X,'***** CHECK OF ORDINATES *****' &
      /8X,'J',6X,'XQ',8X,'YQ',10X,'XP',8X,'YP')
      620 FORMAT(7X,I2,1X,2F10.5,2X,2F10.5)

      END Subroutine
!! ---------------------------------------------------------------------
End Module
!! ---------------------------------------------------------------------

