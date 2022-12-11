!! -----------------------------------------------------------------------
Module SOLVE_M
!! -----------------------------------------------------------------------
      Use glb2DBEM
!! -----------------------------------------------------------------------
      Implicit None
!! -----------------------------------------------------------------------
Contains
!! -----------------------------------------------------------------------
      SUBROUTINE SOLVE(NB, NT, NEQ, AKB, SML, XP, YP, XQ, YQ, VN, ZFI)
      Implicit None
!! -----------Input variables from app2DBEM_OSAKA--------------------------
      Integer(4), Intent(in)          :: NB, NT, NEQ                       ! from app2DBEM_OSAKA
      Double Precision, Intent(in)    :: AKB, SML                          ! from app2DBEM_OSAKA
      !! NB: 단면 패널 수(80)                                            (1)
      !! NT: 단면 패널 수(80) + 3 (Irreg. freq. 방지용 외부 source 3개)   (1)
      !! NEQ: 포텐셜 개수 (Radiation + Diffraction Potential)   (4)      (1)
      !! AKB: Wave Number             (0 ~ 5)                            (1)
      !! SML: Round-off Error 기준값   (1e-14)                           (1)
      !! IPRINT: Round-off Error 확인을 위한 기준 변수                    (1)
!! -----------Variables declaration in SOLVE subroutine---------------------
      Integer(4) :: I, J, K, M
      !! Variables declaration in SDSUB subroutine -------------------------
      Double Precision :: XPI, YPI
      !! XPI: 단면 곡선 각 Panel의 중점 X 방향 위치                       (1)
      !! YPI: 단면 곡선 각 Panel의 중점 Y 방향 위치                       (1)
!! -----------Input variables&arrays from OFFSET_M---------------------------
      Double precision, Allocatable, Intent(in) :: XP(:), YP(:), XQ(:), YQ(:)
      Double Precision, Allocatable, Intent(in) :: VN(:,:)
      !! XP: 단면 곡선 각 Panel의 중점 X 방향 위치                        (1)
      !! YP: 단면 곡선 각 Panel의 중점 Y 방향 위치                        (1)
      !! XQ: 단면 곡선 각 Panel의 Node 점 X 방향 위치                     (1)
      !! YQ: 단면 곡선 각 Panel의 Node 점 Y 방향 위치                     (1)
      !! VN: [d0/dn]                              Matrix                 (1)
!! ----------Arrays declaration in SOLVE subroutine---------------------------
      Complex(8), Allocatable, Intent(out) :: ZFI(:,:)
      Complex(8), Allocatable :: ZAA(:,:), ZBB(:,:), ZSA(:,:), ZSB(:,:)
      Complex(8), Allocatable :: ZS(:), ZD(:)
      Double Precision, Allocatable :: SS(:), DD(:)
      Allocate(ZSA(NT,NB), ZSB(NT,NEQ), ZAA(NB,NB), ZBB(NB,NEQ), ZFI(NEQ,NB))
      Allocate(ZS(NB), ZD(NB), SS(NB), DD(NB))
      !! ZFI: FI 1 ~ FI 4 (Radiation + Diffraction potential)           (1)
      !! ZAA: Modified ZSA: [C(delta) + D]       Matrix                 (1)
      !! ZBB: Modified ZSB: [S][d0/dn]           Vector                 (1)
      !! ZSA: [C(delta) + D]                     Matrix                 (1)
      !! ZSB: [S][d0/dn]                         Matrix                 (1) 
      !! ZS:  [Sij] 또는 [Fij]   Green Fn.: Free-surface term           (1)
      !! ZD:  [Dij] 또는 [Eij]   Green Fn.: Free-surface term           (1)
      !! SS:  [Lij]              Green Fn.: Source-sink term            (1)
      !! DD:  [Tij]              Green Fn.: Source-sink term            (1)
      !! VN: [d0/dn]                                                    (1)
! -------------------------------------------------------------------------
      !! ZAA, ZBB, ZSA, ZSB 행렬 Size 및 복소수 Complex(8) 0행렬 정의    (1)

      DO I=1,NB
            DO J=1,NB
                  ZAA(I,J)=Z0
            End Do

            DO M=1,NEQ
                  ZBB(I,M)=Z0
            End Do
      End Do

      DO I=1,NT
            Do J=1,NB
                  ZSA(I,J)=Z0
            End Do

            DO M=1,NEQ
                  ZSB(I,M)=Z0
            End Do

            IF(I.LE.NB) then
                  ZSA(I,I)=DCMPLX(PI,0.0D0)
            ENDIF
      End Do
!! -------------------------------------------------------------------------
      !!
      Do I=1,NT

            !! 각 패널의 중점 위치를 기준으로 ......
            !! SDSUB에서 Source-sink term(SS / DD) Output 계산
            !! SDCAL에서 Free surface wave term(ZS/ZD) Output 계산
            XPI = XP(I)
            YPI = YP(I)

            !! Source-sink term 1(SS)/2(DD)를 구하는 Subroutine을 Calling.
            Call SDSUB(NEQ, NB, XP, YP, XQ, YQ, VN, XPI, YPI, SS, DD)
            !! Local + Prog. wave term 1(ZS)/2(ZD)를 구하는 Subroutine을 Calling.
            Call SDCAL(NEQ, NB, AKB, XPI, YPI, XQ, YQ, ZS, ZD)

            !! ZSA는 DD(Local wave 1)와 ZD (Prog. wave 1)의 합으로 구성됨.
            DO J=1,NB
                  ZSA(I,J)=ZSA(I,J)+DD(J)+ZD(J)
            End Do

            !! ZSB는 SS(Local wave 2)와 ZS (Prog. wave 2)의 합으로 구성됨.
            !! 여기서, cos(x), sin(x)와 같은 연산이 추가로 이루어짐.
            DO M=1,3
                  DO J=1,NB
                        ZSB(I,M)=ZSB(I,M)+(SS(J)+ZS(J))*VN(M,J)
                  End Do
            End Do

            !! ZSB는 3개의 열을 가지고 있었는데, 여기서 복소지수함수를 반영한 값을 4번째 열을 만들어 추가함.
            ZSB(I,4)=PI2*ZEXP(-AKB*(YP(I)-ZI*XP(I)))

      End Do

      !! 앞서 구한 ZSA의 특정 off diagonal 값들을 곱하고 더해서 ZAA를 만든다!
      !! ZAA: Modified ZSA: [C(delta) + D]       Matrix                 (1)
      Do I=1,NB
            Do J=1,NB
                  Do K=1,NT
                        ZAA(I,J)=ZAA(I,J)+ZSA(K,I)*ZSA(K,J)
                  End Do
            End Do

      !! 앞서 구한 ZSA와 ZSB의 특정 off-diagonal 값들을 곱하고 더해서 ZBB를 만든다!
      !! ZBB: Modified ZSB: [S][d0/dn]           Vector                 (1)

            Do M=1,NEQ
                  Do K=1,NT
                        ZBB(I,M)=ZBB(I,M)+ZSA(K,I)*ZSB(K,M)
                  End Do
            End Do
      End Do

      !! Inverse matrix를 게산하여 실제 해를 구하는 과정임.
      Call ZSWEEP(NB,ZAA,ZBB,NEQ,SML)
      ! Call ZGESV(N, NRHS, ZAA, LDA, IPIV, ZBB, LDB, INFO)
      ! Call ZGESV(100, 100, ZAA, 100, 100, ZBB, 100, 100)
      ! NP = 100, NB = 80

      IF(CDABS(ZAA(1,1)).LT.SML) then
            Write(6,600)
      ENDIF

      600 FORMAT(//10X,'*** ERROR: ZSWEEP IN SUBROUTINE (SOLVE)', &
            ' WAS ABNORMALLY DONE.',/23X,'PLEASE CHECK!'///)

      !! ZBB의 역행렬이 곧 ZFI이며, SOLVE_M의 최종 output 변수임.            (1)
      Do M=1,NEQ
            Do I=1,NB
                  ZFI(M,I)=ZBB(I,M)
            End Do
      End Do

      END Subroutine
! !! ---------------------------------------------------------------------
!!! ---------------------------------------------------------------------
      SUBROUTINE SDSUB(NEQ, NB, XP, YP, XQ, YQ, VN, XPI, YPI, SS, DD)
      !     Kernel Function: Rankine Source
!! -----------Input variables from app2DBEM_OSAKA------------------------------
      Integer :: NEQ, NB
      Double Precision :: AKB
      !! NB: 단면 패널 수(80)                                             (1)
      !! NEQ: 포텐셜 개수 (Radiation + Diffraction Potential)   (4)       (1)
      !! AKB: Wave Number                                      (0 ~ 5)   (1)
!! -----------Input variables in SDSUB-----------------------------------------
      Integer :: J, L
      Double Precision :: XPI, YPI
      Double Precision :: DX, DY, D, SL, CDEL, SDEL, COEF, XA, XB, YA, YB, SUBA, SUBB, WA1, WA2, WA3, ABSC, SWA, DWA ! SDSUB 내 주요 변수들!
      !! XPI: 단면 곡선 각 Panel의 중점 X 방향 위치                       (1)
      !! YPI: 단면 곡선 각 Panel의 중점 Y 방향 위치                       (1)
      !! D: 현 Panel 다음 Panel 중점 간 거리 sqrt(DX^2+DY2)              (1)
      !! DX: 현 Panel 다음 Panel 중점 간 X 축 거리                       (1)
      !! DY: 현 Panel 다음 Panel 중점 간 Y 축 거리                       (1)
      !! SL: once -1, once 1                                            (1)
      !! CDEL: D와 DX에 대한 cos(x) 값                                   (1)
      !! SDEL: D와 DY에 대한 sin(x) 값                                   (1)

      !! XA: 해당 Panel 중점과 왼쪽 노드 사이의 X 방향 차이값             (1)
      !! XB: 해당 Panel 중점과 오른쪽 노드 사이의 X 방향 차이값           (1)
      !! SUBA: cos0*XA + sin0*YA                                        (1)
      !! SUBB: cos0*XB + sin0*YB                                        (1)
      !! COEF: sin0*XA - cos0*YA                                        (1)
      !! ABSC: abs(COEF)                                                (1)
      !! WA1: SS[Lij]'s 2nd term                                        (1)
      !! WA2: SS[Lij]'s 3rd term                                        (1)
      !! WA3: SS[Lij]'s 3rd term / COEF                                 (1)
      !! SWA: [Lij]'s value at each panel                               (1)
      !! DWA: [Tij]'s value at each panel                               (1)
!! -----------Input variables&arrays from OFFSET_M-------------------------
      Double Precision, Allocatable :: XP(:), YP(:), XQ(:), YQ(:)
      Double Precision, Allocatable :: VN(:,:)
      !! XP: 단면 곡선 각 Panel의 중점 X 방향 위치                        (1)
      !! YP: 단면 곡선 각 Panel의 중점 Y 방향 위치                        (1)
      !! XQ: 단면 곡선 각 Panel의 Node 점 X 방향 위치                     (1)
      !! YQ: 단면 곡선 각 Panel의 Node 점 Y 방향 위치                     (1)
      !! VN: [d0/dn]                                                     (1)
!! ----------Arrays declaration in SDSUB subroutine-------------------------
      Double Precision, Allocatable :: SS(:), DD(:)
      !! SS:  [Lij]              Green Fn.: Source-sink term            (1)
      !! DD:  [Tij]              Green Fn.: Source-sink term            (1)

      DO J=1,NB

            SWA=0.0D0
            DWA=0.0D0

            IF(DABS(YPI).LT.1.0D-8) then
                  SS(J)=0.0D0
                  DD(J)=0.0D0
            ENDIF

            DX=XQ(J+1)-XQ(J)
            DY=YQ(J+1)-YQ(J)
            D=DSQRT(DX*DX+DY*DY)
            CDEL=DX/D
            SDEL=DY/D
            XA=XPI-XQ(J  )
            XB=XPI-XQ(J+1)

            SL=-1.0D0

            DO L=1,2
                  SL=-SL
                  YA=SL*YPI-YQ(J  )
                  YB=SL*YPI-YQ(J+1)
                  SUBA=XA*CDEL+YA*SDEL
                  SUBB=XB*CDEL+YB*SDEL
                  COEF=XA*SDEL-YA*CDEL
                  ABSC=DABS(COEF)

                  !! Kernal Function - Source term의 두번째 항!
                  WA1=0.5D0*(SUBB*DLOG(XB*XB+YB*YB)-SUBA*DLOG(XA*XA+YA*YA))

                  IF(ABSC.LT.1.0D-10) THEN
                  WA2=0.0D0
                  WA3=0.0D0
                  ELSE
                  !! Kernal Function - Source term의 세번째 항!
                  WA2=ABSC*(DATAN(SUBB/ABSC)-DATAN(SUBA/ABSC))
                  WA3=WA2/COEF
                  ENDIF

                  SWA=SWA-(WA1+WA2)*SL
                  DWA=DWA+ WA3*SL

            End Do

            SS(J)=SWA
            DD(J)=DWA

      End Do
      End Subroutine
!!! --------------------------------------------------------------------
!! ---------------------------------------------------------------------
      SUBROUTINE SDCAL(NEQ, NB, AKB, XPI, YPI, XQ, YQ, ZS, ZD)
      !     Kernel Function: Free surface Term
!! -----------Input variables from app2DBEM_OSAKA------------------------------
      Integer :: NEQ, NB
      Double Precision :: AKB
      !! NB: 단면 패널 수(80)                                             (1)
      !! NEQ: 포텐셜 개수 (Radiation + Diffraction Potential)   (4)       (1)
      !! AKB: Wave Number                                      (0 ~ 5)   (1)
!! -----------Input variables in SDSUB-----------------------------------------
      Integer :: J, L
      Double Precision :: XPI, YPI
      Double Precision, Allocatable :: XQ(:), YQ(:)
      !! XPI: 단면 곡선 각 Panel의 중점 X 방향 위치                       (1)
      !! YPI: 단면 곡선 각 Panel의 중점 Y 방향 위치                       (1)
      !! XQ: 단면 곡선 각 Panel의 Node 점 X 방향 위치                     (1)
      !! YQ: 단면 곡선 각 Panel의 Node 점 Y 방향 위치                     (1)
!! -----------Input variables in SDCAL-----------------------------------------
      Double Precision :: XX, YY, XE, YE, RFL1, RFT1, RFL2, RFT2, EC, ES, SUB             ! SDCAL 내 주요 변수들!
      Double Precision :: DX, DY, D, CDEL, SDEL, SGNX                                     ! SDCAL 내 주요 변수들!
      Complex(8) :: ZETA, ZFC1, ZFC2, ZFS1, ZFS2, ZSUB                                    ! SDCAL 내 주요 변수들!
      !! XX: |X - zeta|                                                  (1)
      !! YY: (Y + ita)                                                   (1)
      !! XE: -K*(Y + ita)                                                (1)
      !! YE: -K*|X - zeta|                                               (1)
      !! RFL1:                                                           (0)
      !! RFT1:                                                           (0)
      !! RFL2: 0.5*log({K(Y+ita)}^2 - {K*|X-zeta|}^2)                    (1)
      !! RFT2: atan((-i*K*|X - zeta|)/(K*|Y + ita|))                     (1)
      !! EC: EZE1Z function output variable                              (0)
      !! ES: EZE1Z function output variable                              (0)
      !! SUB: SDEL*(RFL2-RFL1)+CDEL*(RFT2-RFT1)                          (1)
      !! D: 현 Panel 다음 Panel 중점 간 거리 sqrt(DX^2+DY2)               (1)
      !! DX: 현 Panel 다음 Panel 중점 간 X 축 거리                        (1)
      !! DY: 현 Panel 다음 Panel 중점 간 Y 축 거리                        (1)
      !! CDEL: D와 DX에 대한 cos(x) 값                                    (1)
      !! SDEL: D와 DY에 대한 sin(x) 값                                    (1)
      !! SGNX: SIGN(X-zeta)                                              (1)
      !! ZETA: -Z: -K*(Y+ita)-i*K*|X-zeta|                               (1)
      !! ZFC1:                                                           (0)
      !! ZFC2: Re[e^(-Z)*E1(-Z)] - pi*i*e^(-Z)              Gf           (1)
      !! ZFS1:                                                           (0)
      !! ZFS2: {Im[e^(-Z)*E1(-Z)] - pi*e^(-Z)}*sgn(x-zeta)  Gfs          (1)
      !! ZSUB: ZSUB=SDEL*(ZFC2-ZFC1)-CDEL*(ZFS2-ZFS1)                    (1)
!! -----------Input variables&arrays from OFFSET_M----------------------------
      Double Precision, Allocatable :: XP(:), YP(:)
      Double Precision, Allocatable :: VN(:,:)
      !! XP: 단면 곡선 각 Panel의 중점 X 방향 위치                         (1)
      !! YP: 단면 곡선 각 Panel의 중점 Y 방향 위치                         (1)
      !! VN: [d0/dn]                                                     (1)
!! -----------Input variables in SDSUB-----------------------------------------
      Double Precision, Allocatable :: SS(:), DD(:)
      !! SS:  [Lij]              Green Fn.: Source-sink term             (1)
      !! DD:  [Tij]              Green Fn.: Source-sink term             (1)
!! ----------Arrays declaration in SDCAL subroutine----------------------------
      Complex(8), Allocatable :: ZS(:), ZD(:)
      !! ZS:  [Sij] 또는 [Fij]   Green Fn.: Free-surface term            (1)
      !! ZD:  [Dij] 또는 [Eij]   Green Fn.: Free-surface term            (1)

      !! ZD, ZS로 복소 영벡터를 정의

      DO J=1,NB
            ZS(J) = DCMPLX(0.0D0, 0.d0)
            ZD(J) = DCMPLX(0.0D0, 0.d0)
      End Do

      XX=XPI-XQ(1)
      YY=YPI+YQ(1)

      SGNX=DSIGN(1.0D0,XX)
      IF(DABS(XX).LT.1.0D-10) then
            SGNX=0.0D0
      ENDIF

      XE=-AKB*YY
      YE=-AKB*DABS(XX)
      ZETA=DCMPLX(XE,YE)

      CALL EZE1Z(XE,YE,EC,ES)
      
      RFL1=0.5D0*DLOG(XX**2+YY**2)
      RFT1=DATAN2(YY,XX)
      ZFC1= EC-(PI*ZEXP(ZETA)*ZI)
      ZFS1=(ES-PI*ZEXP(ZETA))*SGNX

      DO J=1,NB
            XX=XPI-XQ(J+1)
            YY=YPI+YQ(J+1)
            SGNX=DSIGN(1.0D0,XX)

            IF(DABS(XX).LT.1.0D-10) then
                  SGNX=0.0D0
            ENDIF

            !! -K*(Y + ita), -K*|X - zeta|,                                                                  (1)  
            !! -K*(Y + ita)  + i * (-K*|X - zeta|)                                                           (1)
            XE=-AKB*YY
            YE=-AKB*DABS(XX)
            ZETA=DCMPLX(XE,YE)

            CALL EZE1Z(XE,YE,EC,ES)

            !! 0.5*log({K(Y+ita)}^2 - {K*|X-zeta|}^2)                                                        (1)  
            !! atan((-i*K*|X - zeta|)/(K*|Y + ita|))                                                         (1)
            RFL2=0.5D0*DLOG(XX**2+YY**2)
            RFT2=DATAN2(YY,XX)
            !! Re[e^(-Z)*E1(-Z)] - pi*i*e^(-Z)              Gf                                               (1)
            !! {Im[e^(-Z)*E1(-Z)] - pi*e^(-Z)}*sgn(x-zeta)  Gfs                                              (1)
            ZFC2= EC-(PI*ZEXP(ZETA)*ZI)
            ZFS2=(ES-PI*ZEXP(ZETA))*SGNX

            DX=XQ(J+1)-XQ(J)
            DY=YQ(J+1)-YQ(J)
            D =DSQRT(DX*DX+DY*DY)
            CDEL=DX/D
            SDEL=DY/D
            !! ith that - (i-1) th that [sin(0)*ln(sqrt({X-zeta}^2+{Y-ita}^2))+cos0*atan({Y+ita}/{X-zeta})] (1)
            !! ith that - (i-1) th that [sin(0)*Gf -cos0*Gfs]                                               (1)
            SUB =SDEL*(RFL2-RFL1)+CDEL*(RFT2-RFT1)
            ZSUB=SDEL*(ZFC2-ZFC1)-CDEL*(ZFS2-ZFS1)

            !! [Sij] 또는 [Fij]     Green Fn.: Free-surface term                                            (1)
            !! SUB: 1/K*[sin(0)*ln(sqrt({X-zeta}^2+{Y-ita}^2))+cos0*atan({Y+ita}/{X-zeta})]                 (1)
            !! ZSUB: 1/K*[sin(0)*Gf -cos0*Gfs]
            ZS(J)=ZS(J)+2.0D0/AKB*(SUB+ZSUB)
            !! [Dij] 또는 [Eij]   Green Fn.: Free-surface term                                              (1)
            !! {ith Gfs} - {i-1 th Gfs}                                                                     (1)
            ZD(J)=ZD(J)-2.0D0*(ZFS2-ZFS1)

            RFL1=RFL2
            RFT1=RFT2
            ZFC1=ZFC2
            ZFS1=ZFS2
      End Do
      END Subroutine
!! !!-------------------------------------------------------------------
!! !! ------------------------------------------------------------------
SUBROUTINE ZSWEEP(N,ZAA,ZBB,NEQ,EPS)

      Implicit None
      Integer(4), intent(in) :: NEQ
      Integer(4) :: I, J, K, N, IP
      Double Precision :: EPS, P

      Complex(8) :: ZW
      Complex(8), Allocatable :: ZAA(:,:), ZBB(:,:)

      DO 5 K=1,N
            P=0.0D0
            DO 1 I=K,N
            IF(P.GE.CDABS(ZAA(I,K))) GOTO 1
            P=CDABS(ZAA(I,K))
            IP=I
            1 CONTINUE

            IF(P.LE.EPS) GOTO 6
            IF(IP.EQ.K)  GOTO 7

            DO 2 J=K,N
            ZW=ZAA(K,J)
            ZAA(K,J)=ZAA(IP,J)
      2   ZAA(IP,J)=ZW
            DO 20 J=1,NEQ
            ZW=ZBB(K,J)
            ZBB(K,J)=ZBB(IP,J)
      20    ZBB(IP,J)=ZW
      7 CONTINUE

            IF(K.EQ.N) GOTO 70

            DO 3 J=K+1,N
            3 ZAA(K,J)=ZAA(K,J)/ZAA(K,K)

      70   DO 30 J=1,NEQ
      30   ZBB(K,J)=ZBB(K,J)/ZAA(K,K)

            DO 5 I=1,N
            IF(I.EQ.K) GOTO 5
            IF(K.EQ.N) GOTO 40

            DO 4 J=K+1,N
      4   ZAA(I,J)=ZAA(I,J)-ZAA(I,K)*ZAA(K,J)
      40   CONTINUE

            DO 45 J=1,NEQ
      45   ZBB(I,J)=ZBB(I,J)-ZAA(I,K)*ZBB(K,J)
      5 CONTINUE

      ZAA(1,1)=(1.0D0,0.0D0)
      RETURN
      6 ZAA(1,1)=DCMPLX(DABS(P),0.0D0)

      RETURN
      END Subroutine
!!---------------------------------------------------------------------
!! --------------------------------------------------------------------
SUBROUTINE EZE1Z(XX,YY,EC,ES)

      Implicit None

      Integer(4) :: I, J, N
      Double Precision :: XX, YY, X, Y, R, C, ER, EI, SB, FN, CN, CC, SSS, OLD, EXC, EXS, NEW, EC, ES
      Complex(8), allocatable :: Z, Z1, ZSUB, ZS
      !! XX: |X - zeta|                                                  (1)
      !! YY: (Y + ita)                                                   (1)
      !! X:  |X - zeta|                                                  (1) 
      !! Y:  |Y + ita|                                                   (1)
      !! R: sqrt({X - zeta}^2 + {Y + ita}^2)                             (1)
      !! C: atan(|Y + ita|/|X - zeta|)                                   (1)

      X =XX
      Y =DABS(YY)
      R =DSQRT(X*X+Y*Y)
      C =DATAN2(Y,X)

      IF(R.GT.25.0D0) then
            OLD=-1.0D0/R
            EXC=OLD*DCOS(C)
            EXS=OLD*DSIN(C)
      ENDIF

      IF(X.GT.0.0D0.AND.R.GT.8.0D0) then
            Z =DCMPLX(X,Y)
            Z1=(1.0D0,0.0D0)
            ZSUB=(10.0D0,0.0D0)
            ZS  =Z+ZSUB/(Z1+ZSUB/Z)
      ENDIF

      IF(X.LE.0.0D0.AND.Y.GT.10.0D0) then
            Z =DCMPLX(X,Y)
            Z1=(1.0D0,0.0D0)
            ZSUB=(10.0D0,0.0D0)
            ZS  =Z+ZSUB/(Z1+ZSUB/Z)
      ENDIF

      !! R: sqrt({X - zeta}^2 + {Y + ita}^2)                                  (1)
      !! C: atan(|Y + ita|/|X - zeta|)                                        (1)
      !! ER: gamma - log(R) - R * cos(atan(|Y + ita|/|X - zeta|) )            (1)
      !! EI: -atan(|Y + ita|/|X - zeta|) + R*COS(atan(|Y + ita|/|X - zeta|))  (1)
      !! SB: sqrt({X - zeta}^2 + {Y + ita}^2)                                 (1)
      ER=-GAMMA-DLOG(R)+R*DCOS(C)
      EI=-C+R*DSIN(C)
      SB=-R

      Do 100 N=2,100
            FN=DFLOAT(N)
            CN=C*FN
            SB=-SB*R*(FN-1.0D0)/FN/FN
            ER=ER-SB*DCOS(CN)
            EI=EI-SB*DSIN(CN)

            IF(N.EQ.100)  GO TO 1
            IF(EI.EQ.0.0D0)  GO TO 10
            IF(DABS(SB/EI).LE.1.0D-8) GO TO 10
            GO TO 100
            10   IF(DABS(SB/ER).LE.1.0D-8) GO TO 1
      100   CONTINUE

      1     CC=DEXP(X)*DCOS(Y)
            SSS=DEXP(X)*DSIN(Y)
            EC=CC*ER-SSS*EI
            ES=CC*EI+SSS*ER

      IF(YY.LT.0.0D0) ES=-ES
      RETURN

      Do J=1,9
	      ZSUB=DCMPLX(DFLOAT(10-J),0.0D0)
            ZS  =Z+ZSUB/(Z1+ZSUB/ZS)
      End Do

      ZSUB=Z1/ZS
      EC=DREAL(ZSUB)
      ES=DIMAG(ZSUB)

      IF(YY.LT.0.0D0) ES=-ES
      RETURN

      DO 300 N=2,100
            NEW=-OLD/R*DFLOAT(N-1)
            IF(EXS.EQ.0.0D0) GO TO 31
            IF(DABS(NEW/EXS).LE.1.0D-8) GO TO 31
            GO TO 32
            31    IF(EXC.EQ.0.0D0) GO TO 32
                  IF(DABS(NEW/EXC).LE.1.0D-8) GO TO 33
            32    IF(DABS(OLD).LT.DABS(NEW))  GO TO 33
            OLD=NEW
            EXC=EXC+OLD*DCOS(C*DFLOAT(N))
            EXS=EXS+OLD*DSIN(C*DFLOAT(N))
      300   CONTINUE

      33    EC=-EXC
            ES=EXS

      IF(DABS(PI-DABS(C)).LT.1.0D-10) ES=-PI*DEXP(X)
      IF(YY.LT.0.0D0) ES=-ES
      RETURN
      
      END Subroutine

!! --------------------------------------------------------------------
! *  =====================================================================
!       SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
! *
! *  -- LAPACK driver routine --
! *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
! *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! *
! *     .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, N, NRHS
! *     ..
! *     .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX(8)        A( LDA, * ), B( LDB, * )
! *     ..
! *
! *  =====================================================================
! *
! *     .. External Subroutines ..
!       EXTERNAL           XERBLA, ZGETRF, ZGETRS
! *     ..
! *     .. Intrinsic Functions ..
!       INTRINSIC          MAX
! *     ..
! *     .. Executable Statements ..
! *
! *     Test the input parameters.
! *
!       INFO = 0
!       IF( N.LT.0 ) THEN
!          INFO = -1
!       ELSE IF( NRHS.LT.0 ) THEN
!          INFO = -2
!       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
!          INFO = -4
!       ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
!          INFO = -7
!       END IF
!       IF( INFO.NE.0 ) THEN
!          CALL XERBLA( 'ZGESV ', -INFO )
!          RETURN
!       END IF
! *
! *     Compute the LU factorization of A.
! *
!       CALL ZGETRF( N, N, A, LDA, IPIV, INFO )
!       IF( INFO.EQ.0 ) THEN
! *
! *        Solve the system A*X = B, overwriting B with X.
! *
!          CALL ZGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
!      $                INFO )
!       END IF
!       RETURN
! *
! *     End of ZGESV
! *
!       END
! !! --------------------------------------------------------------------



!! --------------------------------------------------------------------
End Module