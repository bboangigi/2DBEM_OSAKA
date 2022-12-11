!! -------------------------------------------------------------------------- !!
Program app2DBEM_OSAKA
!! -------------------------------------------------------------------------- !!
      Use pkg2DBEM_OSAKA
!! -------------------------------------------------------------------------- !!
      Implicit None
      Integer :: IDIF, IRAD, IRAO
!! ----------Initial Variables Argument--------------------------------------- !!
      Integer(4), Parameter :: IPRINT=1, NPRINT=0, NB = 80, NTPlus = 3, NT = 83, NEQ = 4
      Double Precision, Parameter :: H0 = 1.0d0, SIG1 = 0.8d0, SIG2 = 0.8d0, OGD = 0.05d0, KZZB = 0.35d0, SML = 1.0D-14
!! ----------Input / Output Variables Argument related to OFFSET_M-------------!!
      ! Integer(4), Parameter :: MX=105, NP=100, NQ=101   !--Input
      Double Precision :: C22, CMAS, GM
      Double Precision, Allocatable :: XP(:), YP(:), XQ(:), YQ(:)     !--Output
      Double Precision, Allocatable :: VN(:,:)                        !--Output
      !! -----Input / Output Variables Argument related to SDSUB----------!!
      Double Precision, Allocatable :: SS(:), DD(:)                   !--Output
      !! -----Input / Output Variables Argument related to SDCAD for SOLVE----------!!
      Complex(8), Allocatable :: ZS(:), ZD(:)                   !--Output
      !! -----Input / Output Variables Argument related to ZSWEEP for SOLVE----------!!
      Complex(8), Allocatable :: ZAA(:), ZBB(:)                   !--Output
      !! -----Input / Output Variables Argument related to SOLVE for FORCE----------!!
      Complex(8), Allocatable :: ZFI(:,:)                         !--Output

!! ----------Input Variables Argument related to SOLVE_M, Output text files----!!
      Double Precision :: dk
      Double Precision :: AKB
      ! Variables related to SDSUB and SDCAL

      Integer(4) :: iK
      Integer(4), Parameter :: nK  = 300
      Double Precision, Parameter :: kStart = 0.01D0, kEnd = 5.0D0
!! -----------Output Variables Argument related to FORCE_M------------------------
      Complex(8), Allocatable :: ZAB(:,:)
      Complex(8), Allocatable :: ZEXF(:)
!! -------------Input / Output Variables Argument related to MOTION_M-------------
      !! -----Input / Output Variables Argument related to ZSWEEP for MOTION---!!
      Complex(8), Allocatable :: ZAA_local(:,:), ZBB_local(:,:)
!! -------------------------------------------------------------------------------
!! ----------OFFSET Module Calling (Input/output variables corrected!)----------!!
      Call OFFSET(NB, NT, H0, SIG1, SIG2, OGD, KZZB, IPRINT, XP, YP, XQ, YQ, VN, C22, CMAS, GM)
!! -------------------------------------------------------------------------------
!! ----------Output Text file writing command-----------------------------------!!

      open(newunit=IRAO,file="MotionRAO.dat",status="replace")
       write(IRAO,'(a)') "# kB/2 |X2|/A |X3|/A |X4|" // "(kA) Phs(X2) Phs(X3) Phs(X4)"

      open(newunit=IDIF,file="WaveExtForce.dat",status="replace")
       write(IRAO,'(a)') "# kB/2 |F2|/A |F3|/A |F4|" // "(kA) Phs(F2) Phs(F3) Phs(F4)"

      open(newunit=IRAD,file="RadiationForce.dat",status="replace")
      write(IRAD,'(a)') "# nFreq "
      write(IRAD,'(a)') "# kB/2 "
      write(IRAD,'(a)') "# AddedMass(3, 3) "
      write(IRAD,'(a)') "# Damping(3, 3) "
      write(IRAD,'(i5)') nK

!! -------------------------------------------------------------------------------
      dk = (kEnd - kStart) / ( nK - 1.D0 )
      AKB = kStart

      Do iK = 1, nK
            Call SOLVE(NB, NT, NEQ, AKB, SML, XP, YP, XQ, YQ, VN, ZFI)
            Call FORCE (IDIF, IRAD, NB, NEQ, IPRINT, AKB, XQ, YQ, VN, ZFI, ZAB, ZEXF)
            Call MOTION(IRAO, IPRINT, OGD, KZZB, SML, AKB, C22, CMAS, GM, ZAB, ZEXF)
            AKB = AKB + dk
      End Do
! -------------------------------------------------------------------------------
Contains
!! -------------------------------------------------------------------------- !!
End Program app2DBEM_OSAKA
