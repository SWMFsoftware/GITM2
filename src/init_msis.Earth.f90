
subroutine init_msis

  use ModGITM
  use ModInputs
  use ModConstants
  use ModPlanet
  use ModTime

  implicit none

  ! msis variables

  real, dimension(1:2) :: msis_temp
  real, dimension(1:8) :: msis_dens

  integer, dimension(25) :: sw

  integer :: iBlock, iAlt, iLat, iLon, iSpecies
  real :: geo_lat, geo_lon, geo_alt, geo_lst
  real :: ap = 10.0

  call report("init_msis",1)

  !--------------------------------------------------------------------------
  !
  !  From the msis90 library:
  !
  !     OUTPUT:
  !        D(1) - HE NUMBER DENSITY(CM-3)
  !        D(2) - O NUMBER DENSITY(CM-3)
  !        D(3) - N2 NUMBER DENSITY(CM-3)
  !        D(4) - O2 NUMBER DENSITY(CM-3)
  !        D(5) - AR NUMBER DENSITY(CM-3)                       
  !        D(6) - TOTAL MASS DENSITY(GM/CM3)
  !        D(7) - H NUMBER DENSITY(CM-3)
  !        D(8) - N NUMBER DENSITY(CM-3)
  !        T(1) - EXOSPHERIC TEMPERATURE
  !        T(2) - TEMPERATURE AT ALT
  !
  !      TO GET OUTPUT IN M-3 and KG/M3:   CALL METER6(.TRUE.) 
  !
  !      O, H, and N set to zero below 72.5 km
  !      Exospheric temperature set to average for altitudes below 120 km.
  !
  !--------------------------------------------------------------------------

  ! We want units of /m3 and not /cm3

  call meter6(.true.)

  sw = 1

  !  call tselec(sw)

  !           The following is for test and special purposes:
  !            TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW)
  !               WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
  !               FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
  !               FOR THE FOLLOWING VARIATIONS
  !               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
  !               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
  !               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
  !               7 - DIURNAL               8 - SEMIDIURNAL
  !               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
  !              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
  !              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
  !              15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
  !              16 - ALL TINF VAR         17 - ALL TLB VAR
  !              18 - ALL TN1 VAR           19 - ALL S VAR
  !              20 - ALL TN2 VAR           21 - ALL NLB VAR
  !              22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR

  ! Initialize data

  do iBlock = 1, nBlocks
     do iAlt = -1, nAlts+2
        do iLon=1,nLons
           do iLat=1,nLats

              geo_lat = Latitude(iLat,iBlock)*180.0/pi
              geo_lon = Longitude(iLon,iBlock)*180.0/pi

              geo_alt = altitude(iAlt)/1000.0
              geo_lst = mod(utime/3600.0+geo_lon/15.0,24.0)

              !
              ! Call MSIS (results will be im mks units)
              !

              CALL GTD6(iJulianDay,utime,geo_alt,geo_lat,geo_lon,geo_lst, &
                   F107A,F107,AP,48,msis_dens,msis_temp)

              NDensityS(iLon,iLat,iAlt,iH_,iBlock)          = msis_dens(1)
              NDensityS(iLon,iLat,iAlt,iO_,iBlock)          = msis_dens(2)
              NDensityS(iLon,iLat,iAlt,iN2_,iBlock)         = msis_dens(3)
              NDensityS(iLon,iLat,iAlt,iO2_,iBlock)         = msis_dens(4)
              NDensityS(iLon,iLat,iAlt,iAr_,iBlock)         = msis_dens(5)
              NDensityS(iLon,iLat,iAlt,iHe_,iBlock)         = msis_dens(7)
              NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)       = msis_dens(8)
              NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)       = &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)/100.0

              Temperature(iLon,iLat,iAlt,iBlock) = msis_temp(2)/TempUnit
              Rho(iLon,iLat,iAlt,iBlock) = msis_dens(6)

              ! The initial profile of [NO] is refered to:
              !  [Charles A. Barth, AGU, 1995]

              if (geo_alt < 120.) then
                 NDensityS(iLon,iLat,iAlt,iNO_,iBlock)=  &
                      10**(-0.003*(geo_alt-105.)**2 +14+LOG10(3.))
              else 
                 NDensityS(iLon,iLat,iAlt,iNO_,iBlock)=  &
                      MAX(10**(13.-LOG10(3.)*(geo_alt-165.)/35.),1.0)
              endif

              NDensity(iLon,iLat,iAlt,iBlock) = &
                   sum(NDensityS(iLon,iLat,iAlt,:,iBlock))

           enddo
        enddo
     enddo

     Rho(:,:,:,iBlock) = 0.0
     NDensity(:,:,:,iBlock) = 0.0

     do iSpecies=1,nSpecies

        NDensity(:,:,:,iBlock) = NDensity(:,:,:,iBlock) + &
             NDensityS(:,:,:,iSpecies,iBlock)

        Rho(:,:,:,iBlock) = Rho(:,:,:,iBlock) + &
             Mass(iSpecies)*NDensityS(:,:,:,iSpecies,iBlock)
        
     enddo

  enddo
 
end subroutine init_msis

subroutine msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
             F107A,F107,AP,LogNS, Temp, LogRho)

  use ModPlanet
  use ModGITM, only: TempUnit

  implicit none

  real, dimension(1:2) :: msis_temp
  real, dimension(1:8) :: msis_dens

  integer, dimension(25) :: sw

  integer, intent(in) :: iJulianDay
  real, intent(in) :: uTime, Alt, Lat, Lon, LST, f107a, f107, ap
  real, intent(out) :: LogNS(nSpecies), Temp, LogRho
  real :: fac
  integer :: iSpecies

  CALL GTD6(iJulianDay,utime,Alt,Lat,Lon,LST, &
       F107A,F107,AP,48,msis_dens,msis_temp)

  fac = 50.0/sqrt(f107)
  fac = 4.75
!  LogNS(iO_)  = alog(3.0*msis_dens(2))
!  LogNS(iO2_) = alog(msis_dens(4)/2.0)
!  LogNS(iN2_) = alog(msis_dens(3)/2.0)
  LogNS(iO_)  = alog(msis_dens(2))
  LogNS(iO2_) = alog(msis_dens(4))
  if (nSpecies > 2) then
     ! This tricks the compiler...
     iSpecies = iN2_
     LogNS(iSpecies) = alog(msis_dens(3))
  endif
  if (nSpecies > 3) then
     ! This tricks the compiler...
     iSpecies = iN_4S_
     LogNS(iSpecies) = alog(msis_dens(8))
  endif
  if (nSpecies > 4) then
     ! This tricks the compiler...
     iSpecies = iNO_
     LogNS(iSpecies) = alog(8.0e12)
  endif

  Temp        = msis_temp(2)/TempUnit
  LogRho      = alog(msis_dens(6))

end subroutine msis_bcs
