
subroutine calc_GITM_sources(iBlock)

  use ModInputs
  use ModSources
  use ModGITM
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon, iSpecies
  integer :: iiAlt, iiLat, iiLon
  real :: tmp(nLons, nLats, nAlts+1)
  real :: tmp2(nLons, nLats, nAlts), diff(nLons, nLats, nAlts)
  real :: tmp3(nLons, nLats, -1:nAlts+2)
  real :: tmp4(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: tmp5(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: mmm(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: mmr(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: inside(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: Omega(nLons, nLats, nAlts)
  real :: RhoI(nLons, nLats, nAlts)
  real :: ReducedMass(nLons, nLats, nAlts)
  real :: ScaleHeight(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: TmpGradient(nLons, nLats, nAlts,3)

  real :: NSoverN(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: NSoverNGradient(1:nLons, 1:nLats, 1:nAlts, 3)
  real:: LogNum(-1:nLons+2, -1:nLats+2, -1:nAlts+2)

  real :: oldtemp(-1:nLons+2, -1:nLats+2, -1:nAlts+2)

  real :: dtsub, dttotal, dtsubmax, dtr

  real :: diffusion_velocity(nLons, nLats, 0:nAlts+1,nspecies)

  real :: ConS(nLons,nLats,-1:nAlts+2,1:nSpecies)
  real :: LogConS(nLons,nLats,-1:nAlts+2,1:nSpecies)
  real :: GradLogConS(nLons,nLats,1:nAlts,1:nSpecies)

  real :: nVel(1:nAlts, nSpecies)
  real :: NF_Eddy(1:nAlts), NF_NDen(1:nAlts), NF_Temp(1:nAlts)
  real :: NF_NDenS(1:nAlts,1:nSpecies), NF_EddyRatio(1:nAlts,1:nSpecies)
  real :: NF_Gravity(1:nAlts)
  real :: NF_GradLogCon(1:nAlts,1:nSpecies)
  real :: Prandtl(nLons,nLats,0:nalts+1)

! Temporary
  real :: EddyCoefRatio(nLons, nLats, 1:nAlts,nSpecies)

! Potential Temperature
! Used in New Eddy Conduction Calculations:  Bell 1-15-2009
  real :: Theta(nLons, nLats, -1:nAlts+2)
  real :: GammaScale(nLons, nLats, -1:nAlts+2)
  real :: P0(nLons, nLats, -1:nAlts+2)
  real :: Ones(nLons, nLats, -1:nAlts+2)


  call report("calc_GITM_sources",1)

  ! calc_rate is used to determine reaction rates, heating coefficients,
  ! ion-neutral collision frequency, lambdas for ion drag, etc.

  if (iDebugLevel > 4) write(*,*) "=====> going into calc_rates", iproc

  ChemicalHeatingRate = 0.0
  ChemicalHeatingSpecies = 0.0
  call calc_eddy_diffusion_coefficient(iBlock)  
  call calc_rates(iBlock)
  call calc_collisions(iBlock)

  RhoI = IDensityS(1:nLons,1:nLats,1:nAlts,ie_,iBlock) * &
       MeanIonMass(1:nLons,1:nLats,1:nAlts)

  !\
  ! Gravity is a source term which is calculated in initialize.f90
  !/

  !\
  ! ---------------------------------------------------------------
  ! These terms are for Neutral Temperature
  ! ---------------------------------------------------------------
  !/

  !\
  ! Solar Heating -------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> solar heating", iproc

  if (UseSolarHeating .or. UseIonChemistry) then

     ! So far, calc_physics only has stuff that is needed for solar
     ! euv, such as solar zenith angles, and local time.

     call calc_physics(iBlock)

     call euv_ionization_heat(iBlock)

  endif

  if (.not. UseSolarHeating) EuvHeating = 0.0

  ! The auroral heating is specified below, after the aurora is described
  ! in get_potential

  !\
  ! Joule Heating --------------------------------------------------
  !/

  ! If you have turned off IonDrag, then Joule Heating should NOT be
  ! included either.  This is because the Neutrals can become HIGHLY
  ! seperated from the ions, causing MASSIVE Joule heating!!!

  if (UseJouleHeating .and. UseIonDrag) then

     tmp2 = Collisions(1:nLons,1:nLats,1:nAlts,iVIN_) * &
          RhoI(1:nLons,1:nLats,1:nAlts)/ &
          Rho(1:nLons,1:nLats,1:nAlts,iBlock) 

     ! After reading the Thayer paper, this term needs to be taken
     ! out, since it is approximately 1/2, which is the contribution
     ! from ions heating the neutrals through conduction.
     ! * &
     !     (MeanIonMass(1:nLons,1:nLats,1:nAlts)/AMU) /  &
     !     (MeanIonMass(1:nLons,1:nLats,1:nAlts)/AMU + &
     !     MeanMajorMass(1:nLons,1:nLats,1:nAlts)/AMU)

     JouleHeating = 0.0

     do iDir = 1, 3

        JouleHeating(:,:,:) = JouleHeating(:,:,:) + tmp2 * &
             (IVelocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock) - &
             Velocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock))**2

     enddo

     JouleHeating = JouleHeating / &
          TempUnit(1:nLons,1:nLats,1:nAlts) / &
          cp(:,:,1:nAlts,iBlock)

  else

     JouleHeating = 0.0

  endif

  !\
  ! Conduction ----------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> conduction", iproc

  if(UseConduction)then

     tmp2 = Rho(1:nLons, 1:nLats,1:nAlts, iBlock) * &
          cp(1:nLons, 1:nLats,1:nAlts, iBlock)
     

     call calc_conduction(iBlock, &
          Temperature(1:nLons, 1:nLats,-1:nAlts+2, iBlock) * &
          TempUnit(1:nLons, 1:nLats,-1:nAlts+2), &
          KappaTemp(1:nLons, 1:nLats,0:nAlts+1, iBlock), &
          tmp2, &
          MoleConduction)
     
      Conduction = MoleConduction/TempUnit(1:nLons, 1:nLats,1:nAlts) 


!!!
! This term was added by Jared to replace Yue's eddy conduction term in vertical solver.  
!!!

if (UseTurbulentCond) then

    if (UseUpdatedTurbulentCond) then

         Ones(1:nLons,1:nLats,1:nAlts) = 1.0

         do iAlt = -1,nAlts+2
          P0(1:nLons,1:nLats,iAlt) = Pressure(1:nLons,1:nLats,0,iBlock)
         enddo


!! Set the Exponent for the Potential Temperature as (Gamma - 1/Gamma)
!! Must span the range from -1 to nAlts + 2  (same as Theta)

         do iAlt = 1,nAlts
          GammaScale(1:nLons,1:nLats,iAlt) = &
                 ( Gamma(1:nLons,1:nLats,iAlt,iBlock) - Ones(1:nLons,1:nLats,iAlt) )/Gamma(1:nLons,1:nLats,iAlt,iBlock)
         enddo

         do iAlt = -1,0
            GammaScale(1:nLons,1:nLats,iAlt) = GammaScale(1:nLons,1:nLats,1) 
         enddo

         do iAlt = nAlts,nAlts+2
            GammaScale(1:nLons,1:nLats,iAlt) = GammaScale(1:nLons,1:nLats,nAlts) 
         enddo
        

         Theta(1:nLons,1:nLats,-1:nAlts+2) = &
            Temperature(1:nLons,1:nLats,-1:nAlts+2,iBlock)*TempUnit(1:nLons,1:nLats,-1:nAlts+2)*&
            (P0(1:nLons,1:nLats,-1:nAlts+2)/ &
       Pressure(1:nLons,1:nLats,-1:nAlts+2,iBlock))**GammaScale(1:nLons,1:nLats,-1:nAlts+2)
 
!! Prandtl is the Eddy Heat Conduction Coefficient After Hickey et al [2000] 

           Prandtl(1:nLons,1:nLats,0:nAlts+1) = &
             KappaEddyDiffusion(1:nLons,1:nLats,0:nAlts+1,iBlock)*Rho(1:nLons,1:nLats,0:nAlts+1,iBlock)

           tmp2(1:nLons,1:nLats,1:nAlts) = &
              Rho(1:nLons,1:nLats,1:nAlts,iBlock)/Gamma(1:nLons,1:nLats,1:nAlts,iBlock)

          call calc_conduction(iBlock, &
               Theta(1:nLons, 1:nLats,-1:nAlts+2), &
               Prandtl(1:nLons, 1:nLats,0:nAlts+1), &
               tmp2,EddyCond)

!! Eddy Scaling is Set in UAM.in.  Defaults to 1.0

           EddyCondAdia(1:nLons,1:nLats,1:nAlts) = &
               (1.0/EddyScaling)*(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)/Theta(1:nLons,1:nLats,1:nAlts))*&
                EddyCond(1:nLons,1:nLats,1:nAlts)

           UserData3D(1:nLons,1:nLats,1:nAlts,1,iBlock) = EddyCondAdia(1:nLons,1:nLats,1:nAlts)

           Conduction = Conduction + EddyCondAdia
         

    else  !! Use the Old Turbulent Conduction
          Prandtl = 0.0
  
          Prandtl = &
               KappaEddyDiffusion(1:nLons, 1:nLats,0:nAlts+1, iBlock) * &
               Rho(1:nLons, 1:nLats,0:nAlts+1, iBlock) * &
               Cp(1:nLons, 1:nLats,0:nAlts+1, iBlock)

          call calc_conduction(iBlock, &
               Temperature(1:nLons, 1:nLats,-1:nAlts+2, iBlock) * &
               TempUnit(1:nLons, 1:nLats,-1:nAlts+2), &
               Prandtl(1:nLons, 1:nLats,0:nAlts+1), &
               tmp2, &
               EddyCond)

          Conduction = Conduction + EddyCond/TempUnit(1:nLons, 1:nLats,1:nAlts) 

        
        call calc_conduction(iBlock, &
             Pressure(1:nLons, 1:nLats,-1:nAlts+2,iBlock), &
             Prandtl(1:nLons, 1:nLats,0:nAlts+1)/&
             Gamma(1:nLons,1:nLats, 0:nAlts+1,iBlock)/Rho(1:nLons, 1:nLats,0:nAlts+1, iBlock) / &
               Cp(1:nLons, 1:nLats,0:nAlts+1, iBlock),  &
             tmp2, EddyCondAdia)
   
        Conduction = Conduction - &
             EddyCondAdia/TempUnit(1:nLons, 1:nLats,1:nAlts)

     endif ! New Updated Conduction Check
   
 endif


 else
     Conduction = 0.0
 end if





  !\
  ! ---------------------------------------------------------------
  ! These terms are for Neutral Winds
  ! ---------------------------------------------------------------
  !/

  if (UseIonDrag) then

     tmp2 = Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)*&
          RhoI/Rho(1:nLons,1:nLats,1:nAlts,iBlock)

     do iDir = 1, 3
        IonDrag(:,:,:,iDir) = tmp2 * &
             (IVelocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock) - &
             Velocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock))
     enddo

     ! F_iondrag = rho_i/rho * Vis * (Ui-Un)
     ! where Vis = Vin *(Ns/N)  

     do iSpecies = 1, nSpecies
        tmp2 = Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)*&
             RhoI / &
             (Mass(iSpecies) * &
             NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock)) * &
             (NDensityS(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) / &
             NDensity(1:nLons,1:nLats,1:nAlts,iBlock))

        VerticalIonDrag(:,:,:,iSpecies) = tmp2 * &
             (IVelocity(1:nLons,1:nLats,1:nAlts,iUp_,iBlock) - &
             VerticalVelocity(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock))

     enddo

  else

     IonDrag = 0.0

  endif

  !\
  ! ---------------------------------------------------------------
  ! These terms are for Vertical Neutral Wind drag
  ! ---------------------------------------------------------------
  !/

  if (UseNeutralFriction .and. .not.UseNeutralFrictionInSolver) then

     if (UseBoquehoAndBlelly) then

         do iLat = 1, nLats
           do iLon = 1, nLons
             do iAlt = 1, nAlts
               do iSpecies = 1, nSpecies

               GradLogConS(iLon,iLat,iAlt,iSpecies) = &
                  -1.0*Gravity_GB(iLon,iLat,iAlt,iBlock)*&
                   (1.0 - MeanMajorMass(iLon,iLat,iAlt)/Mass(iSpecies))

               enddo 
             enddo 
           enddo 
         enddo 

     else

         do iLat = 1, nLats
           do iLon = 1, nLons
             do iAlt = -1, nAlts+2
               do iSpecies = 1, nSpecies

               ConS(iLon,iLat,iAlt,iSpecies) = &
                  NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)/&
                   NDensity(iLon,iLat,iAlt,iBlock)

               enddo 
             enddo 
           enddo 
         enddo 

         LogConS(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies) = &
             alog(ConS(1:nLons,1:nLats,-1:nAlts+2,1:nSpecies) )

         do iSpecies = 1, nSpecies
           do iAlt = 1, nAlts
                 GradLogConS(1:nLons,1:nLats,iAlt,  iSpecies) = &
               (-1.0*LogConS(1:nLons,1:nLats,iAlt+2,iSpecies) + & 
                 8.0*LogConS(1:nLons,1:nLats,iAlt+1,iSpecies) - &
                 8.0*LogConS(1:nLons,1:nLats,iAlt-1,iSpecies) + &
                 1.0*LogConS(1:nLons,1:nLats,iAlt-2,iSpecies) )/&
                  (12.0*dAlt_GB(1:nLons,1:nLats,iAlt,iBlock))
           enddo 
         enddo 

     endif



!     write(*,*) '==========> Now Entering Neutral Friction Calculation!!'
     do iLat = 1, nLats
        do iLon = 1, nLons

           do iAlt = 1, nAlts
                  NF_NDen(iAlt) = NDensity(iLon,iLat,iAlt,iBlock)
                  NF_Temp(iAlt) = Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt)
                  NF_Eddy(iAlt) = KappaEddyDiffusion(iLon,iLat,iAlt,iBlock)
                  NF_Gravity(iAlt) = Gravity_GB(iLon,iLat,iAlt,iBlock)

             do iSpecies = 1, nSpecies
                  nVel(iAlt,iSpecies) = VerticalVelocity(iLon,iLat,iAlt,iSpecies,iBlock)
                  NF_NDenS(iAlt,iSpecies) = NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)
                  NF_EddyRatio(iAlt,iSpecies) = 0.0
                  NF_GradLogCon(iAlt,iSpecies) = GradLogConS(iLon,iLat,iAlt,iSpecies)

             enddo !iSpecies = 1, nSpecies

           enddo !iAlt = 1, nAlts

           call calc_neutral_friction(nVel(1:nAlts,1:nSpecies), &
                                      NF_Eddy(1:nAlts), &
                                      NF_NDen(1:nAlts), &
                                      NF_NDenS(1:nAlts,1:nSpecies), &
                                      NF_GradLogCon(1:nAlts,1:nSpecies), &
                                      NF_EddyRatio(1:nAlts,1:nSpecies), &
                                      NF_Temp(1:nAlts), NF_Gravity(1:nAlts) )

           do iAlt = 1, nAlts
              NeutralFriction(iLon, iLat, iAlt, 1:nSpecies) = &
                   nVel(iAlt,1:nSpecies) - VerticalVelocity(iLon,iLat,iAlt,1:nSpecies,iBlock) 
!              
!              EddyCoefRatio(iLon, iLat, iAlt, 1:nSpecies,iBlock) = &
!                    NF_EddyRatio(iAlt,1:nSpecies)
!
!              EddyCoefRatio(iLon, iLat, iAlt, 1:nSpecies) = &
!                    NF_EddyRatio(iAlt,1:nSpecies)
           enddo

        enddo
     enddo

  else

     NeutralFriction = 0.0

  endif
!     write(*,*) '==========> Now Exiting Neutral Friction Calculation!!'

  !\
  ! Viscosity ----------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> conduction", iproc

  if(UseViscosity)then

     call calc_viscosity(iBlock)

     call calc_conduction(iBlock, &
          Velocity(1:nLons, 1:nLats,-1:nAlts+2, iNorth_, iBlock), &
          ViscCoef(1:nLons, 1:nLats,0:nAlts+1), &
          Rho(1:nLons, 1:nLats,1:nAlts, iBlock), &
          Viscosity(1:nLons, 1:nLats,1:nAlts, iNorth_))

     call calc_conduction(iBlock, &
          Velocity(1:nLons, 1:nLats,-1:nAlts+2, iNorth_, iBlock), &
          ViscCoef(1:nLons, 1:nLats,0:nAlts+1), &
          Rho(1:nLons, 1:nLats,1:nAlts, iBlock), &
          Viscosity(1:nLons, 1:nLats,1:nAlts, iEast_))


        if (UseEddyViscosity) then

!! Formulation taken from Hickey et al [2000],
!! Waves in Jupiter's Thermosphere

            EddyViscosity = 0.0

           call calc_conduction(iBlock, &
                Velocity(1:nLons, 1:nLats,-1:nAlts+2, iEast_, iBlock), &
                Rho(1:nLons, 1:nLats,0:nAlts+1,iBlock)*ViscCoef(1:nLons, 1:nLats,0:nAlts+1), &
                Rho(1:nLons, 1:nLats,1:nAlts, iBlock), &
                EddyViscosity(1:nLons, 1:nLats,1:nAlts, iEast_))

           call calc_conduction(iBlock, &
                Velocity(1:nLons, 1:nLats,-1:nAlts+2, iNorth_, iBlock), &
                Rho(1:nLons, 1:nLats,0:nAlts+1,iBlock)*ViscCoef(1:nLons, 1:nLats,0:nAlts+1), &
                Rho(1:nLons, 1:nLats,1:nAlts, iBlock), &
                EddyViscosity(1:nLons, 1:nLats,1:nAlts, iNorth_))

            Viscosity(1:nLons,1:nLats,1:nAlts,iEast_) =  &
            Viscosity(1:nLons,1:nLats,1:nAlts,iEast_) - EddyViscosity(1:nLons,1:nLats,1:nAlts,iEast_)

            Viscosity(1:nLons,1:nLats,1:nAlts,iNorth_) =  &
            Viscosity(1:nLons,1:nLats,1:nAlts,iNorth_) - EddyViscosity(1:nLons,1:nLats,1:nAlts,iNorth_)

        endif  ! Use Eddy Viscosity

       Viscosity(:,:,:,iUp_) = 0.0

!! Please Note that the Vertical Viscosity Terms should have Species-Dependent
!! Viscosity Coefficients, the Below Formulation is meant as an approximation

        if (UseVerticalViscosity) then 
            VerticalViscosity = 0.0

            do iSpecies = 1,nSpecies
             call calc_conduction(iBlock, &
                VerticalVelocity(1:nLons, 1:nLats,-1:nAlts+2, iSpecies,iBlock), &
                (4.0/3.0)*ViscCoef(1:nLons, 1:nLats,0:nAlts+1), &
                NDensityS(1:nLons, 1:nLats,1:nAlts, iSpecies,iBlock)*Mass(iSpecies), &
                VerticalViscosity(1:nLons, 1:nLats,1:nAlts, iSpecies))
            enddo

            if (UseEddyViscosity) then
            VerticalEddyViscosity = 0.0

            do iSpecies = 1,nSpecies

             call calc_conduction(iBlock, &
                VerticalVelocity(1:nLons, 1:nLats,-1:nAlts+2, iSpecies,iBlock), &
                NDensityS(1:nLons, 1:nLats,0:nAlts+1, iSpecies,iBlock)*Mass(iSpecies)*&
                KappaEddyDiffusion(1:nLons, 1:nLats,0:nAlts+1,iBlock), &
                NDensityS(1:nLons, 1:nLats,1:nAlts, iSpecies,iBlock)*Mass(iSpecies), &
                VerticalEddyViscosity(1:nLons, 1:nLats,1:nAlts, iSpecies))

             VerticalViscosity(1:nLons,1:nLats,1:nAlts,iSpecies) =  & 
             VerticalViscosity(1:nLons,1:nLats,1:nAlts,iSpecies) -  &
             VerticalEddyViscosity(1:nLons,1:nLats,1:nAlts,iSpecies) 

            enddo

            endif

        endif  ! Test Vertical Viscosity


  else
     Viscosity = 0.0
  end if

  !\
  ! ---------------------------------------------------------------
  ! These terms are for Ion Densities
  ! ---------------------------------------------------------------
  !/

  if (iDebugLevel > 4) write(*,*) "=====> IonChemistry", iproc, UseIonChemistry

  if (iDebugLevel > 4) write(*,*) "=====> calc_ion_v", iproc

  call get_potential(iBlock)
  call calc_efield(iBlock)
  call aurora(iBlock)
  if (UseAuroralHeating) then
     AuroralHeating = AuroralHeatingRate(:,:,:,iBlock) / &
          TempUnit(1:nLons,1:nLats,1:nAlts) / cp(:,:,1:nAlts,iBlock) / &
          rho(1:nLons,1:nLats,1:nAlts, iBlock)
  else
     AuroralHeating = 0.0
  endif

  call calc_ion_v(iBlock)

  ! This includes Radiative Cooling....
  RadCooling = 0.0
  LowAtmosRadRate(:,:,:,iBlock) = 0.0
  call calc_planet_sources(iBlock)

  ! The Emissions array was never set. Should this be here or earlier ????
  Emissions(:,:,:,:,iBlock) = 0.0
  call calc_chemistry(iBlock)

  ChemicalHeatingRate(:,:,:) = &
       ChemicalHeatingRate(:,:,:) * Element_Charge / &
       TempUnit(1:nLons,1:nLats,1:nAlts) / cp(1:nLons,1:nLats,1:nAlts,iBlock)/&
       rho(1:nLons,1:nLats,1:nAlts,iBlock)

  ChemicalHeatingSpecies = ChemicalHeatingSpecies * Element_Charge
end subroutine calc_GITM_sources
