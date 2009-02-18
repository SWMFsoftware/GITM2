subroutine init_b0

  use ModSizeGitm
  use ModGITM
  use ModInputs
  use ModTime

  implicit none

  real :: date, GeoLat, GeoLon, GeoAlt, ALat, ALon, LShell
  real :: xmag, ymag, zmag, r3, MagPot, bmag, MagneticPoleStrength, cD
  real, dimension(3) :: d1,d2,d3,e1,e2,e3
  integer :: iLat, iLon, iBlock, iAlt

  ! We need to fill in the Ghost Cells for the NonChanging Variables, such
  ! as the Magnetic Field:

  call report("init_b0",1)
  call start_timing("init_b0")

  do iBlock=1,nBlocks
     if (nBlocks > 1 .and. iDebugLevel > 1) write(*,*) "==> Block : ", iBlock
     do iAlt=-1,nAlts+2
        do iLat=-1,nLats+2
           do iLon=-1,nLons+2

              GeoLat = Latitude(iLat,iBlock)*180.0/pi
              GeoLon = Longitude(iLon,iBlock)*180.0/pi

              GeoAlt = Altitude_GB(iLon, iLat, iAlt, iBlock)/1000.0

              if (GeoLat > 90.0) then
                 GeoLat = 180.0 - GeoLat
                 GeoLon = mod(GeoLon + 180.0,360.0)
              endif

              if (GeoLat < -90.0) then
                 GeoLat = -180.0 - GeoLat
                 GeoLon = mod(GeoLon + 180.0,360.0)
              endif

              alat = 0.0
              alon = 0.0
              xmag = 0.0
              ymag = 0.0
              zmag = 0.0

              call get_magfield_all(GeoLat,GeoLon,GeoAlt,alat,alon, &
                   xmag,ymag,zmag,d1,d2,d3,e1,e2,e3,cD)

              mLatitude(iLon,iLat,iAlt,iBlock)  = alat
              mLongitude(iLon,iLat,iAlt,iBlock) = alon
              B0(iLon,iLat,iAlt,iEast_,iBlock)  = ymag
              B0(iLon,iLat,iAlt,iNorth_,iBlock) = xmag
              B0(iLon,iLat,iAlt,iUp_,iBlock)    = zmag
              B0(iLon,iLat,iAlt,iMag_,iBlock)   = &
                   sqrt(xmag*xmag + ymag*ymag + zmag*zmag)

              if(UseDynamo)then
                 b0_d1(iLon,iLat,iAlt,:,iBlock) = d1
                 b0_d2(iLon,iLat,iAlt,:,iBlock) = d2
                 b0_d3(iLon,iLat,iAlt,:,iBlock) = d3
                 b0_e1(iLon,iLat,iAlt,:,iBlock) = e1
                 b0_e2(iLon,iLat,iAlt,:,iBlock) = e2
                 b0_e3(iLon,iLat,iAlt,:,iBlock) = e3
                 b0_cD(iLon,iLat,iAlt,iBlock)   = cD
                 b0_Be3(iLon,iLat,iAlt,iBlock)  = &
                      B0(iLon,iLat,iAlt,iMag_,iBlock)/cD
              end if

              if (B0(iLon,iLat,iAlt,iMag_,iBlock) == 0.0) then
                 B0(iLon,iLat,iAlt,iMag_,iBlock) = 1.0e-10
                 DipAngle(iLon,iLat,iAlt, iBlock) = 0.0
                 DecAngle(iLon,iLat,iAlt, iBlock) = 0.0
              else

                 ! Calculate the magnetic dip angle, the magnetic 
                 ! declination angle, and sines and cosines
                 ! For now, only have positive sin of the dip angle 
                 ! (and the magnetic dip angle)
              
                 DipAngle(iLon,iLat,iAlt, iBlock) = &
                      atan (abs(zmag)/sqrt(xmag**2+ymag**2) )
                 DecAngle(iLon,iLat,iAlt, iBlock) = atan2(ymag,xmag) * 180.0/pi

              endif

           enddo
        enddo
     enddo
  enddo

  if (UseApex .and. IsEarth) &
       call dypol( &
       MagneticPoleColat, MagneticPoleLon, MagneticPoleStrength)


  !    GyroFrequency_Electron(:,:,:,iBlock) = &
  !         Element_Charge * B0(:,:,:,Magnitude,iBlock) /  Mass_Electron

  call end_timing("init_b0")
  
end subroutine init_b0

subroutine get_magfield(GeoLat,GeoLon,GeoAlt,xmag,ymag,zmag)

  use ModGITM
  use ModPlanet
  use ModInputs, only: UseApex
  use ModConstants, only: Pi

  implicit none

  real, intent(in) :: GeoLat, GeoLon, GeoAlt
  real, intent(out) :: xmag, ymag, zmag
  real :: r3, bmag

  if (UseApex .and. IsEarth) then

     CALL FELDG(1,GeoLat,GeoLon,GeoALT,XMAG,YMAG,ZMAG,bMag)

     xmag =  xmag * 1.0e-9
     ymag =  ymag * 1.0e-9
     zmag = -zmag * 1.0e-9

  else

     r3 = (RBody / (RBody + GeoAlt*1000.0)) ** 3

     ymag = 0.0
     xmag =     - DipoleStrength * cos(GeoLat*pi/180.0) * r3
     zmag = 2.0 * DipoleStrength * sin(GeoLat*pi/180.0) * r3

  endif

end subroutine get_magfield

!==============================================================================

subroutine get_magfield_all(GeoLat,GeoLon,GeoAlt,alat,alon,xmag,ymag,zmag, &
     d1,d2,d3,e1,e2,e3,CapDMag)

  use ModGITM
  use ModTime
  use ModPlanet
  use ModInputs
  use ModConstants, only: Pi

  implicit none

  real, intent(in) :: GeoLat, GeoLon, GeoAlt
  real, intent(out) :: alat, alon, xmag, ymag, zmag

  real, intent(out) :: d1(3), d2(3), d3(3), CapDMag
  real, intent(out) :: e1(3), e2(3), e3(3)

  real :: CapD(3)
  real :: date, bmag, LShell, r3, MagPot, rBelow, LShell0, mag
  real :: bx, by, bz, twodegrees
  real :: alatp, alatm, alonp, alonm, sinIm, londiff
  integer, external :: jday
  !--------------------------------------------------------------------------
  date = iStartTime(1) + float(iJulianDay)/float(jday(iStartTime(1),12,31))

  twodegrees = 2.0 * pi / 180.0

  rBelow = (2*RadialDistance_GB(1,1,-1,1) - RadialDistance_GB(1,1,1,1)) / RBody

  if (NonMagnetic) then
      xmag = 0.0
      ymag = 0.0
      zmag = 0.0
  endif

  if (UseApex .and. IsEarth) then

     call APEX(DATE,GeoLat,GeoLon,GeoAlt,LShell, &
          alat,alon,bmag,xmag,ymag,zmag,MagPot)

     LShell0 = LShell

     alat = acos(sqrt(rBelow/LShell))*180.0/pi * sign(1.0,aLat)
     xmag =  xmag * 1.0e-9
     ymag =  ymag * 1.0e-9
     zmag = -zmag * 1.0e-9
     mag = sqrt(xmag*xmag + ymag*ymag + zmag*zmag)
     bx = xmag/mag
     by = ymag/mag
     bz = zmag/mag

     ! need to compute the gradient of different apex coordinates

     ! Latitudinal component
     call APEX(DATE,GeoLat+1.0,GeoLon,GeoAlt,LShell, &
          alatp,alonp,bmag,xmag,ymag,zmag,MagPot)
     alatp = acos(sqrt(rBelow/LShell))*180.0/pi * sign(1.0,aLatp)
     
     call APEX(DATE,GeoLat-1.0,GeoLon,GeoAlt,LShell, &
          alatm,alonm,bmag,xmag,ymag,zmag,MagPot)
     alatm = acos(sqrt(rBelow/LShell))*180.0/pi * sign(1.0,aLatm)

     londiff = alonp - alonm
     do while (londiff > 70.0)
        londiff = londiff - 90.0
     enddo
     do while (londiff < -70.0)
        londiff = londiff + 90.0
     enddo
     
     d1(iNorth_) = londiff/twodegrees
     d2(iNorth_) = (alatp - alatm)/twodegrees

     ! Longitudinal component
     call APEX(DATE,GeoLat,mod(GeoLon+1,360.0),GeoAlt,LShell, &
          alatp,alonp,bmag,xmag,ymag,zmag,MagPot)
     alatp = acos(sqrt(rBelow/LShell))*180.0/pi * sign(1.0,aLatp)
     
     call APEX(DATE,GeoLat,mod(GeoLon-1+360,360.0),GeoAlt,LShell, &
          alatm,alonm,bmag,xmag,ymag,zmag,MagPot)
     alatm = acos(sqrt(rBelow/LShell))*180.0/pi * sign(1.0,aLatm)

     londiff = alonp - alonm
     do while (londiff > 70.0)
        londiff = londiff - 90.0
     enddo
     do while (londiff < -70.0)
        londiff = londiff + 90.0
     enddo
     
     d1(iEast_) = (londiff)/(twodegrees * cos(GeoLat*pi/180.0))
     d2(iEast_) = (alatp - alatm)/(twodegrees * cos(GeoLat*pi/180.0))

     ! Altitude component
     call APEX(DATE,GeoLat,GeoLon,GeoAlt+1,LShell, &
          alatp,alonp,bmag,xmag,ymag,zmag,MagPot)
     alatp = acos(sqrt(rBelow/LShell))*180.0/pi * sign(1.0,aLatp)
     
     call APEX(DATE,GeoLat,GeoLon,GeoAlt-1,LShell, &
          alatm,alonm,bmag,xmag,ymag,zmag,MagPot)
     alatm = acos(sqrt(rBelow/LShell))*180.0/pi * sign(1.0,aLatm)

     if (rBelow/LShell > 1.0) then

        write(*,*) 'Reference Altitude in init_b0 : ',rBelow*RBody/1000.0,' km'
        write(*,*) 'This seems to be too high.  Please change the first few'
        write(*,*) 'lines in get_magfield_all.'
        call stop_gitm("Must Stop!!")

     endif

     londiff = alonp - alonm
     do while (londiff > 70.0)
        londiff = londiff - 90.0
     enddo
     do while (londiff < -70.0)
        londiff = londiff + 90.0
     enddo
     
     d1(iUp_) = (RBody+GeoAlt*1000.0)*(londiff)/(2000.0)
     d2(iUp_) = (RBody+GeoAlt*1000.0)*(alatp - alatm)/(2000.0)

     ! Redo the initial calculation to get everything right before
     ! leaving the subroutine....
     call APEX(DATE,GeoLat,GeoLon,GeoAlt,LShell, &
          alat,alon,bmag,xmag,ymag,zmag,MagPot)
     alat = acos(sqrt(rBelow/LShell))*180.0/pi * sign(1.0,aLat)
     xmag =  xmag * 1.0e-9
     ymag =  ymag * 1.0e-9
     zmag = -zmag * 1.0e-9

  else

     if (DipoleStrength /= 0) then

        r3 = (RBody / (RBody + GeoAlt*1000.0)) ** 3
        LShell =  (RBody + GeoAlt*1000.0) / RBody / (sin(pi/2 - GeoLat*pi/180.0))**2.0
        alat = acos(1.0/sqrt(LShell))*180.0/pi * sign(1.0,GeoLat)
        alon = GeoLon
        ymag = 0.0
        xmag =     - DipoleStrength * cos(GeoLat*pi/180.0) * r3
        zmag = 2.0 * DipoleStrength * sin(GeoLat*pi/180.0) * r3

        mag = sqrt(xmag*xmag + ymag*ymag + zmag*zmag)
        bx = xmag/mag
        by = ymag/mag
        bz = zmag/mag

        LShell0 = LShell

        LShell =  (RBody + GeoAlt*1000.0) / RBody / (sin(pi/2-(GeoLat+1)*pi/180.0))**2.0
        alatp = acos(1.0/sqrt(LShell))*180.0/pi * sign(1.0,(GeoLat+1))
        LShell =  (RBody + GeoAlt*1000.0) / RBody / (sin(pi/2-(GeoLat-1)*pi/180.0))**2.0
        alatm = acos(1.0/sqrt(LShell))*180.0/pi * sign(1.0,(GeoLat-1))

        d1(iNorth_) = 0.0
        d2(iNorth_) = (alatp - alatm)/twodegrees

        d1(iEast_) = (2.0)/(twodegrees * cos(GeoLat*pi/180.0))
        d2(iEast_) = 0.0

        LShell =  (RBody + GeoAlt*1000.0+1000.0) / RBody / (sin(pi/2-GeoLat*pi/180.0))**2.0
        alatp = acos(1.0/sqrt(LShell))*180.0/pi * sign(1.0,GeoLat)
        LShell =  (RBody + GeoAlt*1000.0-1000.0) / RBody / (sin(pi/2-GeoLat*pi/180.0))**2.0
        alatm = acos(1.0/sqrt(LShell))*180.0/pi * sign(1.0,GeoLat)

        d1(iUp_) = 0.0
        d2(iUp_) = (RBody+GeoAlt*1000.0)*(alatp - alatm)/(2000.0)

     else

        alat = GeoLat
        alon = GeoLon

        d1 = 0.0
        d2 = 0.0
        LShell0 = rBelow

     endif

  end if

  ! Finish Eqn. 3.8-3.10 from Richmond 1995

  sinIm = 2 * sin(alat*pi/180.0) / sqrt(4.0 - 3.0 * rBelow/LShell0)

  d1 = (d1 * pi / 180.0) * sqrt(rBelow/LShell0) !* sign(1.0,aLat)
  d2 = - (d2 * pi / 180.0) * sinIm

  CapD(iEast_)  =    d1(iNorth_)*d2(iUp_   ) - d2(iNorth_)*d1(iUp_   )
  CapD(iNorth_) = - (d1(iEast_ )*d2(iUp_   ) - d2(iEast_ )*d1(iUp_   ))
  CapD(iUp_)    =    d1(iEast_ )*d2(iNorth_) - d2(iEast_ )*d1(iNorth_)
  CapDMag = sqrt(sum(CapD*CapD))

  if (CapDMag > 0.0) then
     d3(iEast_ ) = by / CapDMag
     d3(iNorth_) = bx / CapDMag
     d3(iUp_   ) = bz / CapDMag
  else
     d3 = 0.0
  endif

  e1(iEast_)  =    d2(iNorth_)*d3(iUp_   ) - d3(iNorth_)*d2(iUp_   )
  e1(iNorth_) = - (d2(iEast_ )*d3(iUp_   ) - d3(iEast_ )*d2(iUp_   ))
  e1(iUp_)    =    d2(iEast_ )*d3(iNorth_) - d3(iEast_ )*d2(iNorth_)

  e2(iEast_)  =    d3(iNorth_)*d1(iUp_   ) - d1(iNorth_)*d3(iUp_   )
  e2(iNorth_) = - (d3(iEast_ )*d1(iUp_   ) - d1(iEast_ )*d3(iUp_   ))
  e2(iUp_)    =    d3(iEast_ )*d1(iNorth_) - d1(iEast_ )*d3(iNorth_)
  
  e3(iEast_)  =    d1(iNorth_)*d2(iUp_   ) - d2(iNorth_)*d1(iUp_   )
  e3(iNorth_) = - (d1(iEast_ )*d2(iUp_   ) - d2(iEast_ )*d1(iUp_   ))
  e3(iUp_)    =    d1(iEast_ )*d2(iNorth_) - d2(iEast_ )*d1(iNorth_)

!!  if (GeoLon > 69 .and. GeoLon < 71) then
!!     write(*,*) GeoLon, GeoLat, GeoAlt
!!     write(*,*) "b0 : ", by,bx,bz
!!     write(*,*) "d1 : ", d1
!!     write(*,*) "d2 : ", d2
!!     write(*,*) "d3 : ", d3
!!     write(*,*) "e1 : ", e1
!!     write(*,*) "e2 : ", e2
!!     write(*,*) "e3 : ", e3
!!  endif

end subroutine get_magfield_all
