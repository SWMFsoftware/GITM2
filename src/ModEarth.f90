
module ModPlanet

  use ModConstants

  implicit none

  integer, parameter :: nSpecies = 4
  integer, parameter :: iO_  = 1
  integer, parameter :: iO2_ = 2
  integer, parameter :: iN2_ = 3
  integer, parameter :: iN_4S_ =  4
  integer, parameter :: iNO_   =  5

  integer, parameter :: nSpeciesTotal = 11
  integer, parameter :: iN_2D_ =  6
  integer, parameter :: iN_2P_ =  7
  integer, parameter :: iH_    =  8
  integer, parameter :: iHe_   =  9
  integer, parameter :: iAr_   = 10
  integer, parameter :: iO_1D_ = 11

  integer, parameter  :: iO_4SP_ = 1
  integer, parameter  :: iO2P_   = 2
  integer, parameter  :: iN2P_   = 3
  integer, parameter  :: iNP_    = 4
  integer, parameter  :: iNOP_   = 5
  integer, parameter  :: iO_2DP_ = 6
  integer, parameter  :: iO_2PP_ = 7
  integer, parameter  :: iHP_    = 8
  integer, parameter  :: iHeP_   = 9
  integer, parameter  :: ie_     = 10
  integer, parameter  :: nIons   = ie_
  integer, parameter  :: nIonsAdvect = 1

  real :: Mass(nSpeciesTotal), MassI(nIons)

  real :: Vibration(nSpeciesTotal)

  integer, parameter :: iE2470_ = 1
  integer, parameter :: iE7320_ = 2
  integer, parameter :: iE3726_ = 3
  integer, parameter :: iE5200_ = 4
  integer, parameter :: iE10400_ = 5
  integer, parameter :: nEmissions = 10

  real, parameter :: GC_Earth               = 9.8                    ! m/s^2
  real, parameter :: RP_Earth               = 24.0*3600.0            ! seconds
  real, parameter :: R_Earth                = 6372.0*1000.0          ! meters
  real, parameter :: DP_Earth               = -31100.0e-9            ! nT

  real, parameter :: Gravitational_Constant = GC_Earth
  real, parameter :: Rotation_Period        = RP_Earth
  real, parameter :: RBody                  = R_Earth
  real, parameter :: DipoleStrength         = DP_Earth

  real, parameter :: OMEGABody              = 2.00*pi/Rotation_Period  ! rad/s

  real, parameter :: HoursPerDay = Rotation_Period / 3600.0
  real, parameter :: Tilt = 23.5

  real, parameter :: DaysPerYear = 365.25
  real, parameter :: SecondsPerYear = DaysPerYear * Rotation_Period

  integer, parameter :: iVernalYear   = 1999
  integer, parameter :: iVernalMonth  =    3
  integer, parameter :: iVernalDay    =   21
  integer, parameter :: iVernalHour   =    0
  integer, parameter :: iVernalMinute =    0
  integer, parameter :: iVernalSecond =    0

  real, parameter :: SunOrbit_A = 1.000110
  real, parameter :: SunOrbit_B = 0.034221
  real, parameter :: SunOrbit_C = 0.001280
  real, parameter :: SunOrbit_D = 0.000719
  real, parameter :: SunOrbit_E = 0.000077

  logical :: IsEarth = .true.
  character (len=10) :: cPlanet = "Earth"

  ! These are for the neutral friction routine...

  ! These are the numerical coefficients in Table 1 in m^2 instead of cm^2
  real, parameter, dimension(4, 4) :: Diff0 = 1.0e4 * reshape( (/ &
       ! 0      02     N2      N     NO
       !---------------------------------+
       0.00,  0.260, 0.260, 0.300, &            ! O
       0.26,  0.000, 0.181, 0.220, &            ! O2
       0.26,  0.181, 0.000, 0.220, &            ! N2
       0.30,  0.220, 0.220, 0.000 /), (/4,4/) )  ! N

  ! These are the exponents
  real, parameter, dimension(4, 4) :: DiffExp = reshape( (/ &
       ! 0      02     N2
       !---------------------------------+
       0.00,  0.75,  0.75, 0.75, &             ! O
       0.75,  0.00,  0.75, 0.75, &             ! O2
       0.75,  0.75,  0.00, 0.75, &             ! N2
       0.75,  0.75,  0.75, 0.00 /), (/4,4/) )  ! N

contains

  subroutine init_planet

    use ModTime

    integer :: itime(7)

    Mass(iH_)    = 1.0 * AMU
    Mass(iHe_)   = 4.0 * AMU
    Mass(iN_4S_) = 14.0 * AMU
    Mass(iO_)    = 16.0 * AMU
    Mass(iN_2D_) = Mass(iN_4S_)
    Mass(iN_2P_) = Mass(iN_4S_)
    Mass(iN2_)   = 2*Mass(iN_4S_)
    Mass(iO2_)   = 2*Mass(iO_)

    Vibration(iO_)    = 5.0
    Vibration(iO2_)   = 7.0
    Vibration(iN2_)   = 7.0
    if (nSpecies > 3) Vibration(iN_4S_) = 5.0
    if (nSpecies > 4) Vibration(iNO_)   = 7.0

    MassI(iO_4SP_) = Mass(iO_)
    MassI(iO_2DP_) = Mass(iO_)
    MassI(iO_2PP_) = Mass(iO_)
    MassI(iO2P_) = Mass(iO2_)
    MassI(iNP_) = Mass(iN_2D_)
    MassI(iN2P_) = Mass(iN2_)
    MassI(iHP_) = Mass(iH_)
    MassI(iHeP_) = Mass(iHe_)
    MassI(iNOP_) = Mass(iN_4S_) + Mass(iO_)
    MassI(ie_) = Mass_Electron

    itime = 0
    itime(1) = iVernalYear
    itime(2) = iVernalMonth
    itime(3) = iVernalDay
    itime(4) = iVernalHour
    itime(5) = iVernalMinute
    itime(6) = iVernalSecond
    call time_int_to_real(itime, VernalTime)

  end subroutine init_planet

end module ModPlanet
