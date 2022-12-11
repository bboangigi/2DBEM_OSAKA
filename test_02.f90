    Program computes

!     .. Parameters ..
      INTEGER          N    !! 행렬크기 N by N
      PARAMETER        ( N = 5 )
      INTEGER          LDA, LDVL, LDVR  !! 별다른 것 없으면 행렬크기랑 같게
      PARAMETER        ( LDA = N, LDVL = N, LDVR = N )
      INTEGER          LWMAX  !! 임의로 큰 숫자
      PARAMETER        ( LWMAX = 1000 )
!
!     .. Local Scalars ..
      INTEGER          INFO, LWORK
!
!     .. Local Arrays ..
      DOUBLE PRECISION A( LDA, N ), VL( LDVL, N ), VR( LDVR, N ),&   !! array 지정
                      WR( N ), WI( N ), WORK( LWMAX )
      DATA             A/&  !! input 행렬
      -1.01, 3.98, 3.30, 4.43, 7.31,&
       0.86, 0.53, 8.26, 4.96,-6.43,&
      -4.60,-7.04,-3.89,-7.66,-6.16,&
       3.31, 5.29, 8.20,-7.33, 2.47,&
      -4.81, 3.55,-1.51, 6.18, 5.58 &
                       /
!
!     .. External Subroutines ..
      EXTERNAL         DGEEV
!
!     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
!
!     Query the optimal workspace.  최적의 LWORK값을 구하는 과정
!

      LWORK = -1     !! -1로 주면 다른 계산 안하고 최적의 LWORK값만 구함
      CALL DGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) ) !! 임시로 줬던 값과 결과로 나온 WORK( 1 ) 값중 작은 것으로
!
!     Solve eigenproblem. 이제 LWORK값을 구했으니 계산 시작
!
      CALL DGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
!
!     Check for convergence.
!
    !   IF( INFO.GT.0 ) THEN  !! 계산이 실패했는지를 판단

    !     WRITE(*,*)'The algorithm failed to compute eigenvalues.'
    !     STOP
    !   ENDIF

    ! Print *, 'LDVL'
    ! Print *, LDVL
    ! Print *, 'LDVR'
    ! Print *, LDVR
    ! Print *, 'WORK'
    ! Print *, WORK
    ! Print *, 'LWORK'
    ! Print *, LWORK 
    ! Print *, 'INFO'    
    ! Print *, INFO

    STOP
    End Program computes