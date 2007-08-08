
module ModElectrodynamics

  use ModSizeGitm

  ! This is the divergence of the Neutral wind current
  real, dimension(-1:nLons+2,-1:nLats+2,1:nAlts) :: DivJu

  ! This is the height integral of the divergence
  real, dimension(-1:nLons+2,-1:nLats+2) :: DivJuAlt

  ! This is the field-line integral of the conductance and divergence
  real, dimension(-1:nLons+2,-1:nLats+2) :: &
       HallFieldLine, PedersenFieldLine, DivJuFieldLine, LengthFieldLine

  ! This is the field-line integral of the conductance and divergence
  real, dimension(-1:nLons+2,-1:nLats+2) :: &
       SigmaPP, SigmaLL, SigmaHH, SigmaCC, SigmaPL, SigmaLP, &
       KDpm, KDlm

  ! This is the field aligned integral in magnetic coordinates
  real, dimension(:,:), allocatable :: DivJuAltMC

  ! These are the conductances in magnetic coordinates
  real, dimension(:,:), allocatable :: SigmaHallMC
  real, dimension(:,:), allocatable :: SigmaPedersenMC

  real, dimension(:,:), allocatable :: SigmaPPMC
  real, dimension(:,:), allocatable :: SigmaLLMC
  real, dimension(:,:), allocatable :: SigmaHHMC
  real, dimension(:,:), allocatable :: SigmaCCMC
  real, dimension(:,:), allocatable :: SigmaPLMC
  real, dimension(:,:), allocatable :: SigmaLPMC

  real, dimension(:,:), allocatable :: KDpmMC
  real, dimension(:,:), allocatable :: KDlmMC

  ! These are the magnetic coordinates
  real, dimension(:,:), allocatable :: MagLatMC
  real, dimension(:,:), allocatable :: MagLocTimeMC
  real, dimension(:,:), allocatable :: MagLonMC
  real, dimension(:,:), allocatable :: GeoLonMC
  real, dimension(:,:), allocatable :: GeoLatMC

  real, dimension(:,:), allocatable :: MagBufferMC
  real, dimension(:,:), allocatable :: LengthMC

  real, dimension(:,:), allocatable :: &
       solver_a_mc, solver_b_mc, solver_c_mc, solver_d_mc, solver_e_mc, &
       solver_s_mc, deltalmc, deltapmc, &
       dSigmaLLdlMC, dSigmaLPdlMC, dSigmaPLdpMC, dSigmaPPdpMC, &
       dKDpmdpMC, dKDlmdlMC, DynamoPotentialMC

  real, dimension(:,:), allocatable :: oldpotmc
  
  integer :: nMagLats = 90  ! 2 degrees
  integer :: nMagLons = 72  ! 5 degrees

  !----------------------------------------------------------------------
  ! These are in geographic coordinates : 

  real, dimension(-1:nLons+2,-1:nLats+2, nBlocksMax) :: &
       HallConductance, PedersenConductance

  real, dimension(-1:nLons+2,-1:nLats+2, -1: nAlts+2) :: &
       Sigma_0, Sigma_Pedersen, Sigma_Hall

  real, dimension(-1:nLons+2,-1:nLats+2, -1: nAlts+2, 3) :: &
       UxB, Ju

  real :: SigmaR(-1:nLons+2,-1:nLats+2, -1: nAlts+2, 3, 3)

end module ModElectrodynamics
