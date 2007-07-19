
subroutine euv_ionization_heat(iBlock)

  use ModGITM
  use ModEUV
  use ModPlanet
  use ModConstants
  use ModInputs
  use ModSources
  use ModTime, only : tSimulation, CurrentTime
  use ModIndicesInterfaces

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iWave, iSpecies, iIon, iError
  real, dimension(nLons,nLats) :: Tau, Intensity

  logical :: IsFirstTime(nBlocksMax) = .true.

  real :: photoion(Num_WaveLengths_High, nIons-1)
  real :: photoabs(Num_WaveLengths_High, nSpecies)
  real :: photodis(Num_WaveLengths_High, nSpecies)
  real :: NeutralDensity(nLons, nLats, nSpecies)
  real :: ChapmanLittle(nLons, nLats, nSpecies)
  real :: EHeat(nLons, nLats)

  if (IsFirstTime(iBlock)) then
     IsFirstTime(iBlock) = .false.
  else
     if (floor((tSimulation - dT)/dTAurora) == &
          floor(tSimulation/dTAurora)) return
  endif

  call report("euv_ionization_heat",2)
  call start_timing("euv_ionization_heat")

  iError = 0
  call get_f107(CurrentTime, f107, iError)
  if (iError /= 0) then
     write(*,*) "Error in getting F107 value.  Is this set?"
     write(*,*) "Code : ",iError
     call stop_gitm("Stopping in euv_ionization_heat")
  endif

  call get_f107a(CurrentTime, f107a, iError)
  if (iError /= 0) then
     write(*,*) "Error in getting F107a value.  Is this set?"
     write(*,*) "Code : ",iError
     call stop_gitm("Stopping in euv_ionization_heat")
  endif

  call chapman_integrals(iBlock)

  EuvIonRate = 0.0
  EuvHeating(:,:,:,iBlock)= 0.0
  eEuvHeating(:,:,:,iBlock) = 0.0
  EuvIonRateS(:,:,:,:,iBlock) = 0.0 
  EuvDissRateS(:,:,:,:,iBlock) = 0.0

  photoion = 0.0
  photoabs = 0.0
  photodis = 0.0

  ! This transfers the specific photo absorption and ionization cross
  ! sections into general variables, so we can use loops...

  call fill_photo(photoion, photoabs, photodis)

  do iAlt = 1, nAlts

     NeutralDensity = NDensityS(1:nLons,1:nLats,iAlt,1:nSpecies,iBlock)
     ChapmanLittle  = Chapman(:,:,iAlt,1:nSpecies,iBlock)
     EHeat = 0.0

     do iWave = 1, Num_WaveLengths_High

        Tau = 0.0
        do iSpecies = 1, nSpecies
           Tau = Tau + &
                photoabs(iWave, iSpecies) * ChapmanLittle(:,:,iSpecies)
        enddo

        Intensity = Flux_of_EUV(iWave) * exp(-1.0*Tau)

        do iIon = 1, nIons-1
           EuvIonRateS(:,:,iAlt,iIon,iBlock) = &
                EuvIonRateS(:,:,iAlt,iIon,iBlock) + &
                Intensity*PhotoIon(iWave,iIon)
        enddo

        do iSpecies = 1, nSpecies
           EuvDissRateS(:,:,iAlt,iSpecies,iBlock) = &
                EuvDissRateS(:,:,iAlt,iSpecies,iBlock) + &
                Intensity*PhotoDis(iWave,iSpecies)
        enddo

        do iSpecies = 1, nSpecies
           EHeat = EHeat + &
                Intensity*PhotonEnergy(iWave)* &
                photoabs(iWave, iSpecies) * NeutralDensity(:,:,iSpecies)
        enddo

     enddo

     EuvHeating(:,:,iAlt,iBlock)  = EHeat*HeatingEfficiency_CB(:,:,iAlt,iBlock)
     eEuvHeating(:,:,iAlt,iBlock) = EHeat*eHeatingEfficiency_CB(:,:,iAlt,iBlock)

  enddo

  !\
  ! Zero out EuvHeating if specified not to use it.
  !/

  if (UseSolarHeating) then
     do iAlt = 1, nAlts
        EuvHeating(:,:,iAlt,iBlock) = EuvHeating(:,:,iAlt,iBlock) / &
           Rho(1:nLons,1:nLats,iAlt,iBlock) / &
           cp(1:nLons,1:nLats,iAlt,iBlock) / &
           TempUnit(1:nLons,1:nLats,iAlt)
     enddo
  else
     EuvHeating = 0.0
  endif

  call end_timing("euv_ionization_heat")

end subroutine euv_ionization_heat


!-------------------------------------------------------------------
! Subroutine for calculating the EUV flux in a vacuum.
!-------------------------------------------------------------------

subroutine calc_euv

  use ModEUV
  use ModInputs

  implicit none

  integer :: i
  real    :: flxfac, wavelength_ave

  !:::::::::::::::::::::::::::::::: EUVAC :::::::::::::::::::::::
  !------ This EUV flux model uses the F74113 solar reference spectrum and
  !------ ratios determined from Hinteregger's SERF1 model. It uses the daily
  !------ F10.7 flux (F107) and the 81 day mean (F107A) as a proxy for 
  !------ scaling
  !------ The fluxes are returned in EUVFLX and correspond to the 37
  !------ wavelength
  !------ bins of Torr et al. [1979] Geophys. Res. Lett. p771.
  !------ See Richards et al. [1994] J. Geophys. Res. p8981 for details.
  !
  !...... F107   = input daily 10.7 cm flux index. (e.g. 74)
  !...... F107A  = input 81 day ave. of daily F10.7 centered on current day
  !...... EUVFLX = output array for EUV flux in units of photons/cm2/sec.
  !
  !
  !----- loop through the wavelengths calculating the scaling factors and
  !----- the resulting solar flux.
  !----- The scaling factors are restricted to be greater than 0.8
  !

  do i = 1, Num_waveLengths_Low

     FLXFAC=(1.0 + AFAC(I) * (0.5*(F107+F107A) - 80.0))
     IF(FLXFAC.LT.0.8) FLXFAC=0.8
     EUV_Flux(i) = F74113(I) * FLXFAC * 1.0E9 * 10000.

 enddo

end subroutine calc_euv

!-------------------------------------------------------------------
! Subroutine for calculating scaled solar flux.
!-------------------------------------------------------------------

subroutine calc_scaled_euv

  use ModEUV
  use ModInputs
  use ModTime

  implicit none

  integer, parameter :: Hinteregger_Contrast_Ratio  = 0
  integer, parameter :: Hinteregger_Linear_Interp   = 1
  integer, parameter :: Tobiska_EUV91               = 2
  integer, parameter :: WoodsAndRottman_10Nov88     = 3
  integer, parameter :: WoodsAndRottman_20Jun89     = 4

  integer :: N, NN
  real    :: f107_Ratio, r1, r2, hlybr, fexvir, hlya, heiew
  real    :: xuvfac, hlymod, heimod, xuvf, wavelength_ave
  real (Real8_) :: rtime
  integer, dimension(7) :: Time_Array

 !DAVES:
  real :: wvavg(Num_WaveLengths_High)
  character (len=2) :: dday, dhour, dminute 
  character (len=7) :: dtime

  ! regression coefficients which reduce to solar min. spectrum:
  ! for Hinteregger_Contrast_Ratio model:

  real, dimension(1:3) :: B1, B2
  real, dimension(Num_WaveLengths_High) :: Timed_Flux
  data B1/1.0, 0.0138, 0.005/
  data B2/1.0, 0.59425, 0.3811/

  ! 'best fit' regression coefficients, commented out, for reference:
  !     DATA B1/1.31, 0.01106, 0.00492/, B2/-6.618, 0.66159, 0.38319/

  call calc_euv

  hlybr = 0.
  fexvir = 0.
  hlya = 3.E+11 + 0.4E+10 * (f107-70.)
  heiew = 0.

  f107_ratio = (f107-68.0) / (243.0-68.0)

  xuvfac = 4.0 - f107_ratio
  if (xuvfac < 1.0) xuvfac = 1.0

  do N = 1, Num_WaveLengths_High
     Solar_Flux(N) = RFLUX(N) + (XFLUX(N)-RFLUX(N)) * f107_Ratio
  enddo

  iModelSolar = Tobiska_EUV91
!  iModelSolar = Hinteregger_Contrast_Ratio

  select case(iModelSolar)

  case (Hinteregger_Contrast_Ratio)

     if (hlybr > 0.001) then
        r1 = hlybr
     else
        r1 =  B1(1) + B1(2)*(f107A-71.5) + B1(3)*(f107-f107A+3.9)
     endif
     if (fexvir > 0.001) THEN
        r2 = fexvir
     else
        r2 =  B2(1) + B2(2)*(f107A-71.5) + B2(3)*(f107-f107A+3.9)
     endif
     do N = 13, Num_WaveLengths_High
        Solar_Flux(N) = (RFLUX(N) + ((r1-1.)*SCALE1(N)              &
             +  (r2-1.)*SCALE2(N)) / 1000.)
     enddo

  case (Hinteregger_Linear_Interp)

     ! do nothing

  case (Tobiska_EUV91)

     if (HLYA > 0.001) then

        hlymod = hlya

     else

        if (heiew > 0.001) then
           hlymod = heiew * 3.77847e9 + 8.40317e10
        else
           hlymod = 8.70e8 * F107 + 1.90e11
        endif

     endif

     if (heiew > 0.001) then
        heimod = heiew * 3.77847e9 + 8.40317e10
     else
        heimod = hlymod
     endif

     do N=16,55
        Solar_Flux(N) = TCHR0(N)        + &
             TCHR1(N)*hlymod + &
             TCHR2(N)*heimod + &
             TCOR0(N)        + &
             TCOR1(N)*f107   + &
             TCOR2(N)*f107A
     enddo

  case (WoodsAndRottman_10Nov88)

     DO N=15,55
        Solar_Flux(N) = WAR1(N)
     enddo

  case (WoodsAndRottman_20Jun89)

     DO N=15,55
        Solar_Flux(N) = WAR2(N)
     enddo

  end select

  !
  ! Substitute in H Lyman-alpha and XUVFAC if provided:
  !

  if (hlya > 0.001) Solar_Flux(12) = hlya / 1.E9
  if (xuvfac > 0.001) THEN
     xuvf = xuvfac
  else
     xuvf = 1.0
  endif

  !
  ! Convert from gigaphotons to photons, cm^-2 to m^-2, etc.:
  !

  do N=1,Num_WaveLengths_High

     IF (Solar_Flux(N) < 0.0) Solar_Flux(N) = 0.0

     IF ((WAVEL(N) < 251.0) .AND. (WAVES(N) > 15.0)) then
        Solar_Flux(N) = Solar_Flux(N)*xuvf
     endif

     !
     ! Convert to photons/m^2/s
     !

     Solar_Flux(N) = Solar_Flux(N) * 1.E9 * 10000.0

     !
     ! Convert to eV/m^2/s
     !


     wavelength_ave = (WAVEL(N) + WAVES(N))/2.0
     PhotonEnergy(N)= 6.626e-34*2.998e8/(wavelength_ave*1.0e-10)

     !
     !       Solar_Flux(N) = Solar_Flux(N) * &
     !            (Planck_Constant * Speed_Light) / &
     !            (wavelength_ave * 1.0e-10 * Element_Charge)

  enddo

  do N = 16,Num_WaveLengths_High
     Flux_of_EUV(N) = Solar_Flux(N)
  enddo

  do N = 1,15
     Flux_of_EUV(N) = Solar_Flux(N)
  enddo

  do N=1,Num_WaveLengths_Low
     NN = N+15
     Flux_of_EUV(NN) = 0.5*(EUV_Flux(N)+Solar_Flux(NN))
  enddo

  Flux_of_EUV = Flux_of_EUV/(SunOrbitEccentricity**2)

  do N=1,Num_WaveLengths_High
     wvavg(N)=(WAVEL(N)+WAVES(N))/2.
  enddo

  if (UseEUVData) then

     call start_timing("new_euv")
     open(unit = iInputUnit_, file=cEUVFile)
     rtime=0
     do while(rtime < CurrentTime)
        read(iInputUnit_,*) Time_Array(1:6),Timed_Flux
        Time_Array(7)=0
        call time_int_to_real(Time_Array,rtime)
     enddo
     close(iInputUnit_)
     !!need to convert from W/m^2 to photons/m^2/s
     do N=1,Num_WaveLengths_High 
        Flux_of_EUV(N) = Timed_Flux(N)*wvavg(N)*1.0e-10/(6.626e-34*2.998e8)
     enddo
     call end_timing("new_euv")

  endif

  ! Second Spectra, provided by Steve Bougher....

!  do n = 1, nS2WaveLengths
!     S2PhotonEnergy(N) = 6.626e-34*2.998e8/(S2WaveLengths(n)*1.0e-10)
!  enddo

end subroutine calc_scaled_euv

!-------------------------------------------------------------------
! Subroutine for initializing photoabsorption, photoionization,
! and branching ratio quantities.
!-------------------------------------------------------------------

subroutine init_euv

  use ModEUV

  implicit none

  integer :: N, NN

  call report("init_euv",2)

  EUVEFF = 0.05

  do N = 1,15

     RLMSRC(N) = RLMEUV(N)

  enddo

  PhotoAbs_CO  = 0.0
  PhotoAbs_CO2 = 0.0
  PhotoIon_CO  = 0.0
  PhotoIon_CO2 = 0.0

  DO N = 1, Num_WaveLengths_Low

     NN = N+15    ! 16:52

     BranchingRatio_OPlus2P(N) = 0.
     IF (N.GT.14) then
        BranchingRatio_OPlus2P(N) =                                   &
             1. - BranchingRatio_OPlus2D(N) - BranchingRatio_OPlus4S(N)
     endif

     PhotoAbs_O2(NN)      = Photoabsorption_O2(N)
     PhotoAbs_O(NN)       = Photoabsorption_O(N) 
     PhotoAbs_N2(NN)      = Photoabsorption_N2(N)

     PhotoIon_O2(NN)      = Photoionization_O2(N)
     PhotoIon_OPlus4S(NN) = Photoionization_O(N)*BranchingRatio_OPlus4S(N)
     PhotoIon_N2(NN)      = Photoionization_N2(N)
     PhotoIon_N(NN)       = Photoionization_N(N)
     PhotoIon_OPlus2D(NN) = Photoionization_O(N)*BranchingRatio_OPlus2D(N)
     PhotoIon_OPlus2P(NN) = Photoionization_O(N)*BranchingRatio_OPlus2P(N)

     BranchingRatio_N2(NN) = BranchingRatio_N2_to_NPlus(N)
     BranchingRatio_O2(NN) = BranchingRatio_O2_to_OPlus(N)

     ! Since these are the 37 Torr bands, we should put the cross
     ! sections provided by Steve Bougher here...
     !  ** I use 80 wavelength centered bins (old scheme). 
     !    -- last two (79 and 80) : soft xrays (10-50A)
     !    -- 42-78 bins : Torr 37 wavelength bands or lines (50-1050 A)
     !    -- 1-41 : some lines (including Lyman-alpha) plus SR bands
     !              and continuum
     !---------------------------------------------------------------
     ! There are differences in the Spectrum!!!  We need to investigate
     ! this further!!!
     !---------------------------------------------------------------

     PhotoAbs_CO(NN)  = S2PhotoAbsCO(41+N) * 1.0e-18 / 10000.0
     PhotoAbs_CO2(NN) = S2PhotoAbsCO2(41+N) * 1.0e-18 / 10000.0
     PhotoIon_CO(NN)  = S2PhotoIonCO(41+N) * 1.0e-18 / 10000.0
     PhotoIon_CO2(NN) = S2PhotoIonCO2(41+N) * 1.0e-18 / 10000.0

  enddo

  ! Convert everything from /cm2 to /m2

  DO N = 1, Num_WaveLengths_High

     PhotoAbs_O2(N)      = PhotoAbs_O2(N)  / 10000.0
     PhotoAbs_O(N)       = PhotoAbs_O(N)   / 10000.0
     PhotoAbs_N2(N)      = PhotoAbs_N2(N)  / 10000.0

     PhotoIon_O2(N)      = PhotoIon_O2(N)        / 10000.0
     PhotoIon_OPlus4S(N) = PhotoIon_OPlus4S(N)   / 10000.0
     PhotoIon_N2(N)      = PhotoIon_N2(N)        / 10000.0
     PhotoIon_N(N)       = PhotoIon_N(N)         / 10000.0
     PhotoIon_OPlus2D(N) = PhotoIon_OPlus2D(N)   / 10000.0
     PhotoIon_OPlus2P(N) = PhotoIon_OPlus2P(N)   / 10000.0

  enddo

!  do n = 1, nS2WaveLengths
!
!     S2PhotoAbsCO  = S2PhotoAbsCO * 1.0e-18 / 10000.0
!     S2PhotoAbsCO2 = S2PhotoAbsCO2 * 1.0E-18 / 10000.0
!     S2PhotoAbsN2  = S2PhotoAbsN2 * 1.0E-18 / 10000.0
!     S2PhotoAbsO2  = S2PhotoAbsO2 * 1.0E-18 / 10000.0
!     S2PhotoIonN2  = S2PhotoIonN2 * 1.0E-18 / 10000.0
!     S2PhotoIonCO  = S2PhotoIonCO * 1.0E-18 / 10000.0
!     S2PhotoIonCO2 = S2PhotoIonCO2 * 1.0E-18 / 10000.0
!     S2PhotoIonO2  = S2PhotoIonO2 * 1.0E-18 / 10000.0
!
!  enddo

end subroutine init_euv

