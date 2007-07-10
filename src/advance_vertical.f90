subroutine advance_vertical(iLon,iLat,iBlock)

  use ModGITM
  use ModPlanet, only: nSpecies, OmegaBody, nIonsAdvect
  use ModConstants, only: pi
  use ModSources, only: EUVHeating, VerticalTempSource, KappaEddyDiffusion
  use ModInputs, only: UseIonAdvection, iDebugLevel
  use ModVertical, ONLY: &
       LogRho, &
       cMax1      => cMax,&
       KappaTemp1 => KappaTemp, &
       LogNS1     => LogNS, &
       Heating, &
       Vel_GD,  &
       Temp, &
       Centrifugal, Coriolis, &
       LogINS, &
       IVel, Lat, Lon, &
       VertVel, &
       MeanMajorMass_1d, &
       gamma_1d, &
       EddyCoef_1d

  implicit none

  integer, intent(in) :: iLon, iLat, iBlock

  integer :: iIon, iSpecies, iAlt, iDim

!  KappaTemp1 = KappaTemp(iLon,iLat,:,iBlock)

  EddyCoef_1d(1:nAlts) = KappaEddyDiffusion(iLon,iLat,1:nAlts,iBlock)

  
  if (minval(NDensityS(iLon,iLat,:,1:nSpecies,iBlock)) <= 0.0) then
     write(*,*) "negative density found!"
     call stop_gitm("Can't Continue")
  endif

  Heating     = EuvHeating(iLon,iLat,:,iBlock)
  Centrifugal = (cos(Latitude(iLat,iBlock)) * OmegaBody)**2
  Coriolis    = 2 * cos(Latitude(iLat,iBlock)) * OmegaBody
  LogRho  = log(Rho(iLon,iLat,:,iBlock))
  do iDim = 1, 3 
     Vel_GD(:,iDim)  = Velocity(iLon,iLat,:,iDim,iBlock)
  enddo

  !!!! CHANGE !!!!

  VerticalTempSource(iLon,iLat,1:nAlts) = &
       Temp(1:nAlts)/TempUnit(iLon,iLat,1:nAlts)-Temperature(iLon,iLat,1:nAlts,iBlock)

  Temp    = Temperature(iLon,iLat,:,iBlock)*TempUnit(iLon,iLat,:)
  do iSpecies = 1, nSpecies 
     LogNS1(:,iSpecies)  = log(NDensityS(iLon,iLat,:,iSpecies,iBlock))
     VertVel(:,iSpecies) = VerticalVelocity(iLon,iLat,:,iSpecies,iBlock)
  enddo

  cMax1   = cMax_GDB(iLon,iLat,:,iUp_,iBlock)

  do iDim = 1, 3 
     IVel(:,iDim) = IVelocity(iLon,iLat,:,iDim,iBlock)
  enddo

  do iSpecies = 1, nIonsAdvect
     LogINS(:,iSpecies)  = log(IDensityS(iLon,iLat,:,iSpecies,iBlock))
  enddo

  MeanMajorMass_1d = MeanMajorMass(iLon,iLat,:)
  gamma_1d=gamma(ilon,ilat,:,iBlock)     

!!!!  LogINS  = IDensityS(iLon,iLat,:,1:nIonsAdvect,iBlock)

  Lat = Latitude(iLat, iBlock) * 180.0/pi
  Lon = Longitude(iLon, iBlock) * 180.0/pi

  call advance_vertical_1d

  Rho(iLon,iLat,:,iBlock)                  = exp(LogRho)

  do iDim = 1, 3 
     Velocity(iLon,iLat,:,iDim,iBlock)           = Vel_GD(:,iDim)
  enddo

  !!!! CHANGE !!!!
  Temperature(iLon,iLat,:,iBlock)          = Temp/TempUnit(iLon,iLat,:)

  do iSpecies = 1, nSpecies 
     LogNS(iLon,iLat,:,iSpecies,iBlock)              = LogNS1(:,iSpecies)
     VerticalVelocity(iLon,iLat,:,iSpecies,iBlock) = VertVel(:,iSpecies)
  enddo

  if (minval(Temp) < 0.0) then
     write(*,*) "Temperature is negative!!!"
     do iAlt = -1,nAlts+2
        if (Temp(iAlt) < 0.0) then
           write(*,*) "iAlt : ", iAlt, Temp(iAlt)
        endif
     enddo
     call stop_gitm("Can't continue")
  endif

  if (Maxval(LogNS1) > 75.0) then
     write(*,*) "Maxval of LogNS too high!!!"
     do iAlt = -1,nAlts+2
        do iSpecies = 1, nSpecies
           if (LogNS1(iAlt,iSpecies) > 75.0) then
              write(*,*) "iSpecies, iBlock, Alt,Lon, Lat, maxval : ",&
                   iSpecies,iBlock,&
                   altitude(iAlt)/1000.0, &
                   longitude(iLon,iBlock)*180/pi, &
                   latitude(iLat,iBlock)*180/pi, LogNS1(iAlt,ispecies)
           endif
        enddo
     enddo
     call stop_gitm("Can't continue")
  endif

  do iSpecies = 1, nSpecies 
     nDensityS(iLon,iLat,:,iSpecies,iBlock) = exp(LogNS1(:,iSpecies))
  enddo

  if (UseIonAdvection) then

     do iIon = 1, nIonsAdvect
        IDensityS(iLon,iLat,:,iIon,iBlock) = exp(LogINS(:,iIon))
     enddo

!     if (Maxval(LogINS) > 75.0) then
!        write(*,*) "Maxval of Ion LogINS too high!!!"
!        do iAlt = -1,nAlts+2
!           do iSpecies = 1, nIonsAdvect
!              if (LogINS(iAlt,iSpecies) > 75.0) then
!                 write(*,*) "iSpecies, iBlock, Alt,Lon, Lat, maxval : ",&
!                      iSpecies,iBlock,&
!                      altitude(iAlt)/1000.0, &
!                      longitude(iLon,iBlock)*180/pi, &
!                      latitude(iLat,iBlock)*180/pi, LogINS(iAlt,ispecies)
!              endif
!           enddo
!        enddo
!        call stop_gitm("Can't continue")
!     endif
!

!!!!!     IDensityS(iLon,iLat,:,1:nIonsAdvect,iBlock) = LogINS
     !\
     ! New Electron Density
     !/
     IDensityS(iLon,iLat,:,ie_,iBlock) = 0.0
     do iIon = 1, nIons-1
        IDensityS(iLon,iLat,:,ie_,iBlock) = &
             IDensityS(iLon,iLat,:,ie_,iBlock) + &
             IDensityS(iLon,iLat,:,iIon,iBlock)
     enddo
  endif

end subroutine advance_vertical
