!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module UA_wrapper

  ! Wrapper for GITM Upper Atmosphere (UA) component

  implicit none

  private ! except

  public:: UA_set_param
  public:: UA_init_session
  public:: UA_run
  public:: UA_save_restart
  public:: UA_finalize

  ! CON_coupler_points
  public:: UA_find_points
  public:: UA_get_grid_info

  ! GM coupler
  public:: UA_get_for_gm
  public:: UA_put_from_gm_dt
  public:: UA_put_from_gm_init
  public:: UA_put_from_gm

contains

  !============================================================================
  subroutine UA_set_param(CompInfo, TypeAction)

    use ModInputs, only: cInputText
    use ModReadParam, only: read_text, n_line_read

    use ModTime, ONLY: StartTime, tSimulation, CurrentTime
    use ModInputs, only: iStartTime, IsFramework, iOutputUnit_, set_defaults, &
         nInputLines
    use ModTimeConvert, ONLY: time_real_to_int
    use CON_physics,    ONLY: get_time
    use ModIoUnit
    use ModProcUA
    use ModGITM, only: iCommGITM, nProcs, iProcGITM => iProc
    use ModPlanet, only: init_planet
    use CON_comp_info
    use ModUtilities, ONLY: check_dir

    character (len=*), parameter :: NameSub='UA_set_param'

    ! Arguments
    type(CompInfoType), intent(inout):: CompInfo   ! Information for this comp.
    character (len=*), intent(in)    :: TypeAction ! What to do

    integer :: iError

    iError = 0

    !-------------------------------------------------------------------------
    write(*,*) "-->Starting UA_set_param..."
    select case(TypeAction)
    case('VERSION')

       call put(CompInfo,&
            Use=.true.,                                      &
            NameVersion='Global Iono-Thermo Model (Ridley)', &
            Version=2.0)

    case('MPI')

       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)

       iCommGITM = iComm
       iProcGITM = iProc
       nProcs    = nProc

       if(iProc==0)then
          call check_dir("UA/DataIn")
          call check_dir("UA/data")
          call check_dir("UA/RestartOUT")
       end if

       IsFramework = .true.

       call init_planet
       call set_defaults

    case('READ')

       call read_text(cInputText)
       cInputText(n_line_read()+1) = "#END"
       nInputLines=n_line_read()+1

       call set_inputs

    case('CHECK')

       iError = 0

    case('STDOUT')

       !     iUnitStdOut=STDOUT_
       !     if(nProc==1)then
       !        StringPrefix='UA:'
       !     else
       !        write(StringPrefix,'(a,i3.3,a)')'UA',iProc,':'
       !     end if

    case('FILEOUT')

       call get(CompInfo,iUnitOut=iOutputUnit_)
       !     StringPrefix=''

    case('GRID')

       call UA_set_grid

    case default

       call CON_stop(NameSub//' UA_ERROR: invalid TypeAction='//TypeAction)

    end select

    if (iError /= 0) &
         call CON_stop(NameSub//' UA_ERROR in TypeAction='//TypeAction)

  end subroutine UA_set_param

  !============================================================================

  subroutine UA_set_grid

    ! Set the grid descriptor for UA
    ! Since UA has a static grid the descriptor has to be set once.
    ! There can be many couplers that attempt to set the descriptor,
    ! so we must check IsInitialized.
    use ModProcUA
    use CON_Coupler
    use CON_comp_info
    use ModNumConst
    use ModSizeGitm
    use ModSphereInterface, only: iStartBLK
    use ModInputs, only: nBlocksLat, nBlocksLon
    use ModUtilities, ONLY: check_allocate

    character (len=*), parameter :: NameSub='UA_set_grid'
    logical :: IsInitialized=.false.
    integer, allocatable :: iProc_A(:), iProcPE_A(:)

    integer :: iError, iBlock, iBlockPE

    real, allocatable :: CoLat_I(:), Lon_I(:), Alt_I(:)
    real, allocatable :: LatPE_I(:), LonPE_I(:)

    logical :: DoTest, DoTestMe, Done

    !------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest, DoTestMe)
    if(DoTest) write(*,*)NameSub,' IsInitialized=',IsInitialized

    write(*,*) "-->Starting UA_set_grid..."

    if(IsInitialized) return

    IsInitialized=.true.

    if(iProc>=0)then
       allocate(CoLat_I(nBlocksLat*nLats), &
            Lon_I(nBlocksLon*nLons), &
            iProc_A(nBlocksLat*nBlocksLon), &
            LatPE_I(nBlocksLat*nBlocksLon*nLats), &
            LonPE_I(nBlocksLat*nBlocksLon*nLons), &
            iProcPE_A(nBlocksLat*nBlocksLon),&
            Alt_I(-1:nAlts+2),&
            stat = iError)

       if (iError /= 0) then
          write(*,*) NameSub, " Error in allocating variables"
          write(*,*) " Lat_I, Lon_I, iProc_A, LatPE_I, LonPE_I, iProcPE_A"
          call CON_stop(NameSub//' UA_ERROR')
       endif

       LatPE_I   = -1.0e32
       LonPE_I   = -1.0e32
       iProcPE_A = -1

       do iBlockPE = 1, nBlocks

          iBlock = iStartBLK + iBlockPE

                  !LatPE_I((iBlock-1)*nLats+1:iBlock*nLats) = &
                  !     Latitude(1:nLats,iBlockPE)
                  !LonPE_I((iBlock-1)*nLons+1:iBlock*nLons) = &
                  !     Longitude(1:nLons,iBlockPE)
                  write(*,*) 'iProc: ',iProc
          iProcPE_A(iBlock) = iProc

       enddo

       ! call MPI_allreduce( LatPE_I, CoLat_I, nBlocksLon*nBlocksLat*nLats, &
       !          MPI_REAL, MPI_MAX, iComm, iError)
       !     ! Save into colatitudes instead of latitude
       !     CoLat_I = cHalfPi - CoLat_I
       !
       !     call MPI_allreduce( LonPE_I, Lon_I, nBlocksLon*nBlocksLat*nLons, &
       !          MPI_REAL, MPI_MAX, iComm, iError)

       call MPI_allreduce( iProcPE_A, iProc_A, nBlocksLon*nBlocksLat, &
            MPI_INTEGER, MPI_MAX, iComm, iError)
       !     Alt_I=Altitude(:)
    else
       allocate( CoLat_I(1), Lon_I(1),iProc_A(1),Alt_I(1),stat=iError)
       call check_allocate(iError,NameSub)
    end if

    write(*,*) 'nCell: ',nLats,nLons,nAlts

    call set_grid_descriptor(                        &
         UA_,                                        &! component index
         nDim=3,                                     &! dimensionality
         nRootBlock_D=(/nBlocksLat,nBlocksLon,1/),     &! blocks
         nCell_D =(/nLats,nLons,nAlts/),             &! size of node based grid
         XyzMin_D=(/cHalf,cHalf,cHalf/),                   &! generalize coord
         XyzMax_D=(/nLats-cHalf,nLons-cHalf,nAlts-cHalf/), &! generalize coord
         TypeCoord='GEO',                            &! magnetic coordinates
         !Coord1_I= CoLat_I,                          &! colatitudes
         !Coord2_I= Lon_I,                            &! longitudes
         !Coord3_I= Alt_I,                            &! radial size in meters
         iProc_A = iProc_A,                          &! processor assigment
         IsPeriodic_D=(/.false.,.true.,.false./))     ! periodic in longitude

    !call set_grid_descriptor(                        &
    !     UA_,                                        &! component index
    !     nDim=3,                                     &! dimensionality
    !     nRootBlock_D=(/nBlocksLat,nBlocksLon,1/),     &! blocks
    !     nCell_D =(/nLats,nLons,nAlts/),             &! size of node based grid
    !     XyzMin_D=(/cHalf,cHalf,cHalf/),                   &! generalize coord
    !     XyzMax_D=(/nLats-cHalf,nLons-cHalf,nAlts-cHalf/), &! generalize coord
    !     TypeCoord='GEO',                            &! magnetic coordinates
    !     NameVar = 'Bx BY Bz CO2n', &
         !Coord1_I= CoLat_I,                          &! colatitudes
         !Coord2_I= Lon_I,                            &! longitudes
         !Coord3_I= Alt_I,                            &! radial size in meters
    !     iProc_A = iProc_A,                          &! processor assigment
    !     IsPeriodic_D=(/.false.,.true.,.false./))     ! periodic in longitude

    !call set_coord_system(UA_, &
    !     TypeCoord = 'GEO', &
         !UnitX     = No2Si_V(UnitX_), &
    !     nVar      = 4, &
    !     NameVar   = 'Bx BY Bz CO2n', &
    !     TypeGeometry = 'spherical')!, &
         !Coord1_I     = CoLat_I, &!(/ RadiusMin, RadiusMax /), &
         !Coord2_I     = Lon_I, &!(/ CoordMin_D(2), CoordMax_D(2) /), &
         !Coord3_I     = Alt_I )!(/ CoordMin_D(3), CoordMax_D(3) /)  )

  end subroutine UA_set_grid

  !============================================================================

  subroutine UA_init_session(iSession, SWMFTime)

    use CON_physics,    ONLY: get_time
    use ModTime, only : StartTime, iTimeArray, CurrentTime

    real, intent(in)    :: SWMFTime
    integer, intent(in) :: iSession

    logical :: IsFirstTime = .true.

    write(*,*) "-->Starting UA_init_session..."
    if (IsFirstTime) then

       ! Set time related variables for UA
       call get_time(tStartOut = StartTime)

       CurrentTime = StartTime + SWMFTime
       call time_real_to_int(StartTime, iTimeArray)

       call fix_vernal_time

       call initialize_gitm(CurrentTime)
       call write_output

       IsFirstTime = .false.

    endif

  end subroutine UA_init_session

  !==========================================================================

  subroutine UA_run(SWMFTime, SWMFTimeLimit)

    use ModGITM
    use ModTime
    use ModInputs, only: iDebugLevel, Is1D
    use ModTimeConvert, ONLY: time_real_to_int, n_day_of_year

    save

    real, intent(in)    :: SWMFTimeLimit
    real, intent(inout) :: SWMFTime

    integer :: index,lat,long,s,i,j,a, k, jj
    logical :: done, status_ok
    integer :: ierror, CLAWiter, n
    real    :: maxi,tt
    integer :: time_array(7)
    !real :: nDenNuSpecies_srcInterface_CBI

    logical :: exist, IsDone

    CurrentTime = StartTime + SWMFTime
    EndTime     = StartTime + SWMFTimeLimit

    write(*,*) "-->Starting UA_run..."
    if (iDebugLevel > 1) then
       call time_real_to_int(CurrentTime, time_array)
       write(*,"(a,i5,5i3,i3)") "> Running UA from time : ",time_array(1:7)
    endif

    call calc_pressure

    Dt = 1.e32

    call calc_timestep_vertical
    if (.not. Is1D) call calc_timestep_horizontal

    if (iDebugLevel > 1) write(*,"(a,f13.5)") "> UA_run Dt : ",Dt

    call advance

    iStep = iStep + 1

    call write_output

    SWMFTime = CurrentTime - StartTime

  end subroutine UA_run

  !==========================================================================

  subroutine UA_save_restart(TimeSimulation)

    use ModInputs

    real, intent(in) :: TimeSimulation

    character(len=*), parameter :: NameSub='UA_save_restart'
    !--------------------------------------------------------------------------
    call write_restart("UA/restartOUT/")

  end subroutine UA_save_restart

  !==========================================================================

  subroutine UA_finalize(TimeSimulation)

    real, intent(in) :: TimeSimulation

    character(len=*), parameter :: NameSub='UA_finalize'
    !--------------------------------------------------------------------------
    call finalize_gitm

  end subroutine UA_finalize

  !============================================================================
  subroutine UA_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

    !use BATL_lib,   ONLY: MaxDim, find_grid_block
    use ModPhysics, ONLY: Si2No_V, UnitX_
    use ModSizeGitm
    use ModSphereInterface
    use ModCoordTransform, ONLY: xyz_to_rlonlat

    integer, intent(in) :: nDimIn                ! dimension of position vector
    integer, intent(in) :: nPoint                ! number of positions
    real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions // Xyz_DI is the same as Pos_DI!!!
    integer, intent(out):: iProc_I(nPoint)       ! processor owning position

    ! Find array of points and return processor indexes owning them
    ! Could be generalized to return multiple processors...

    ! COPIED FROM EE !!!!
    !real:: Xyz_D(MaxDim) = 0.0
    integer:: iPoint, iBlock
    real :: LatFind, LonFind, AltFind
    integer :: iiLat, iiLon, iiAlt, iiBlock, iiProc,iAlt
    real :: rLon, rLat, rAlt
    integer :: iLons,iLats,iAlts
    real :: r,Lon,Lat,Alt
    real :: rLonLat_D(3),Xyz_D(3)

    character(len=*), parameter:: NameSub = 'UA_find_points'
    !--------------------------------------------------------------------------
    write(*,*) "-->Starting UA_find_points..."
    write(*,*) 'nPoint in UA_find_points: ', nPoint, nBlocks
    write(*,*) 'Positions: ',Xyz_DI(:,10),maxval(Xyz_DI),minval(Xyz_DI) ! GM positions in Mars radii
    !LatFind = 0.1
    !LonFind = 0.1
    !call BlockLocationIndex(LonFind,LatFind,1,iiLon,iiLat,rLon,rLat)
    !write(*,*) '-->-->--> iiLon: ',iiLon,iiLat,rLon,rLat
    !AltFind = 150e3
    !call BlockAltIndex(AltFind,1,iiLon,iiLat,iiAlt,rAlt)
    !write(*,*) '-->-->--> iiAlt: ',iiLon,iiLat,iiAlt,rLon,rLat

    write(*,*) 'nLons,nLats,nAlts: ',nLons,nLats,nAlts

    do iPoint = 1, nPoint
       Xyz_D(1:nDimIn) = Xyz_DI(:,iPoint)*R_Mars!Si2No_V(UnitX_)
       ! Convert to spherical... xyz_to_sph ModCoordTransform
       write(*,*) '-.-.->',Xyz_D
       call xyz_to_rlonlat(Xyz_D, rLonLat_D)
       r = rLonLat_D(1)
       Alt = r - R_Mars
       Lon = rLonLat_D(2)
       Lat = rLonLat_D(3)
       write(*,*) 'rLonLat: ',Alt,Lon,Lat
       !write(*,*) 'MinAlt, MaxAlt: ',minval(Altitude_GB),maxval(Altitude_GB)
       call LocationProcIndex(Lon, Lat, Alt, iiBlock, iiLon, iiLat, iAlt, rLon, rLat, rAlt, iiProc)
       !call LocationIndex(Lon, Lat, iiBlock, iiLon, iiLat, rLon, rLat)
       write(*,*) 'Indices: ',iiProc
       iProc_I(iPoint) = iiProc

       ! Inputs: LonFind = Desired longitude
       !         LatFind = Desired latitude
       !
       ! Outputs: iiBlock = Block index containing the desired location
       !          iiLon   = Longitude index for LonFind
       !          iiLat   = Latitude index for LatFind
       !          rLon    = Longitude interpolation scaling factor
       !          rLat    = Latitude interpolation scaling factor


       !iProc_I(iPoint) =

       !write(*,*) "Xyz_D: ",Xyz_D
       !call find_grid_block(Xyz_D, iProc_I(iPoint), iBlock)  Some day...
    end do
    write(*,*) 'Finished iProc: ',minval(iProc_I),maxval(iProc_I)
    ! COPIED FROM EE !!!!
  !--------------------------------------------------------------------------
  end subroutine UA_find_points
  !============================================================================
  subroutine UA_get_grid_info(nDimOut, iGridOut, iDecompOut)

    use ModInputs, ONLY: Is1D,IsFullSphere
    use ModPlanet,  ONLY: iNewGridGITM, iNewDecompositionGITM

    integer, intent(out):: nDimOut    ! grid dimensionality
    integer, intent(out):: iGridOut   ! grid index
    integer, intent(out):: iDecompOut ! decomposition index

    logical:: IsFirstTime = .true.
    integer:: nDimGITM = -1

    character(len=*), parameter :: NameSub = 'UA_get_grid_info'

    ! Return basic grid information useful for model coupling.
    ! The decomposition index increases with load balance and AMR.
    !--------------------------------------------------------------------------
    write(*,*) "-->Starting UA_get_grid_info...",Is1D,IsFullSphere!, nDim, iNewGrid, iNewDecomposition
    if (Is1D) nDimGITM = 1
    if (IsFullSphere) nDimGITM = 3
    write(*,*) 'nDimGITM: ',nDimGITM
    nDimOut    = nDimGITM
    iGridOut   = iNewGridGITM
    iDecompOut = iNewDecompositionGITM

  end subroutine UA_get_grid_info
  !============================================================================
  subroutine UA_get_for_gm(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)

    use ModGITM,  ONLY: iProc,nDimGITM
    use ModPhysics, ONLY: Si2No_V, UnitX_, No2Si_V, iUnitCons_V
    use ModAdvance, ONLY: State_VGB, Bx_, Bz_, nVar
    use ModVarIndexes, ONLY: nVar
    use ModB0,      ONLY: UseB0, get_b0
    use BATL_lib,   ONLY: nDim, MaxDim, MinIJK_D, MaxIJK_D, find_grid_block
    use ModIO, ONLY: iUnitOut
    use CON_coupler, ONLY: nVarBuffer, iVarSource_V
    !use BATL_lib,   ONLY: MaxDim, MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    !use ModInterpolate, ONLY: interpolate_vector
    use ModInterpolate, ONLY: trilinear
    use ModInputs, ONLY: AltMin, AltMax
    use ModPlanet, ONLY: RBody, nSpeciesTotal
    use ModGeometry, ONLY: R_BLK
    use ModCoordTransform, ONLY: xyz_to_rlonlat
    use ModSizeGitm
    use ModSphereInterface
    use ModUACoupling

    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI / Xyz_DI is the same as Pos_DI!!!

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors / Xyz_DI is the same as Pos_DI!!!
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    real:: Xyz_D(MaxDim), B0_D(MaxDim)
    real:: Dist_D(MaxDim), State_V(1)!State_V(nVar)
    integer:: iCell_D(MaxDim)

    integer :: i

    integer, allocatable, save:: iBlockCell_DI(:,:)
    real,    allocatable, save:: Dist_DI(:,:)
    real,    allocatable:: NDensityS_reshaped(:,:,:,:,:)

    integer:: iPoint, iBlock, iProcFound, iVarBuffer, iVar
    real :: LatFind, LonFind, AltFind
    integer :: iiLat, iiLon, iiAlt, iiBlock, iiProc,iAlt
    real :: rLon, rLat, rAlt
    integer :: iLons,iLats,iAlts,iSpecies
    real :: r,Lon,Lat,Alt
    real :: rLonLat_D(3)

    logical:: DoTest, DoTestMe

    character(len=*), parameter :: NameSub='UA_get_for_gm'
    !--------------------------------------------------------------------------
    !write(*,*) "-->Starting UA_get_for_gm..."
    write(*,*) "Calling UA_get_for_gm...", nPoint, nVarIn
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! AltMin = 0.0 and AltMax = 300.0 (from set_inputs), user RBody + ...
    write(*,*) 'AltMin: ',AltMin
    write(*,*) 'AltMax: ',AltMax
    write(*,*) 'RBody: ',RBody
    write(*,*) 'NDim: ',nDim
    write(*,*) 'NPoint: ',nPoint,nLons,nLats,nAlts    ! why is nPoint = 0???
    ! If nDim < MaxDim, make sure that all elements are initialized
    Dist_D = -1.0
    Xyz_D  =  0.0
    if(.not.allocated(NDensityS_reshaped)) then
       allocate(NDensityS_reshaped(nSpeciesTotal,nLons,nLats,nAlts,nBlocks))
       NDensityS_reshaped = 0.0
       do iSpecies = 1, nSpeciesTotal
         do iLons = 1, nLons
           do iLats = 1, nLats
             do iAlts = 1, nAlts
               do iBlock = 1, nBlocks
                 NDensityS_reshaped(iSpecies,iLons,iLats,iAlts,iBlock) = &
                    NDensityS(iLons,iLats,iAlts,iSpecies,iBlock)
               end do!iBlock
             end do!iAlts
           end do!iLats
         end do!iLons
       end do!iSpecies
    end if

    ! Need to go through all points in nPoint and write data from GITM into an
    ! array that can be used later on by BATS-R-US (GM uses Data_VI). Region in
    ! space needs to be defined as well in cartesian!!!

    do iPoint = 1, nPoint
       Xyz_D(1:nDimIn) = Xyz_DI(:,iPoint)*R_Mars!Si2No_V(UnitX_)
       ! Convert to spherical... xyz_to_sph ModCoordTransform
       !write(*,*) '-.-.->',Xyz_D
       call xyz_to_rlonlat(Xyz_D, rLonLat_D)
       r = rLonLat_D(1)
       Alt = r - R_Mars
       Lon = rLonLat_D(2)
       Lat = rLonLat_D(3)
       !write(*,*) 'rLonLat: ',Alt,Lon,Lat
       !write(*,*) 'MinAlt, MaxAlt: ',minval(Altitude_GB),maxval(Altitude_GB)
       call LocationProcIndex(Lon, Lat, Alt, iiBlock, iiLon, iiLat, iAlt, rLon, rLat, rAlt, iiProc)
       Dist_D(1) = rLon
       Dist_D(2) = rLat
       Dist_D(3) = rAlt
       iCell_D(1) = iiLon
       iCell_D(2) = iiLat
       iCell_D(3) = iAlt
       !call LocationIndex(Lon, Lat, iiBlock, iiLon, iiLat, rLon, rLat)
       write(*,*) 'Indices: ',iiProc,iBlock,iiBlock
       write(*,*) 'iCell_D: ',iCell_D(1),iCell_D(2),iCell_D(3)
       write(*,*) 'Dist_D: ',Dist_D(1),Dist_D(2),Dist_D(3)
       write(*,*) 'Point1: ',rLonLat_D,Alt
       write(*,*) 'Point2: ',Longitude(iiLon,iiBlock),Latitude(iiLat,iiBlock),Altitude_GB(iiLon,iiLat,iAlt,iiBlock)

       !iProc_I(iPoint) = iiProc
       if(iiProc /= iProc)then
          write(*,*)NameSub,' ERROR: Xyz_D, iProcFound=', Xyz_D, iProcFound
          call stop_mpi(NameSub//' could not find position on this proc')
       end if

       !State_V = trilinear(State_VGB(:,:,:,:,iBlock), &
      !      nVar, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, Xyz_D, &
      !      iCell_D = iCell_D, Dist_D = Dist_D)
      write(*,*) 'Pre-trilinear...',iPoint,size(Xyz_DI),nVar,iiBlock,maxval(NDensityS_reshaped)
      write(*,*) '111: ',NDensityS_reshaped(:,2,6,42,iiBlock)
      !write(*,*) 'Pre-trilinear...',State_VGB(:,:,:,:,iBlock)!NDensityS(:,:,:,iCO2_,iiBlock)
       !State_V = trilinear(State_VGB(:,:,:,:,iBlock), &
       ! Array must be (nvar,i,j,k,iblock)
       State_V = trilinear(NDensityS_reshaped(:,:,:,:,iiBlock), &         ! Why is nVar 12?!!!
            !nVar, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, Xyz_D, &
            1, 1, nLons, 1, nLats, 1, nAlts, Xyz_D, &
            iCell_D = iCell_D, Dist_D = Dist_D)
       write(*,*) 'Post-trilinear...',iPoint,size(State_V),size(Data_VI(:,iPoint)),Size(Xyz_D),Size(Data_VI)
       Data_VI(:,iPoint) = State_V!*No2Si_V(iUnitCons_V)
       !NdensityS(1:nLons,1:nLats,1:nAlts,iCO2_,iBlock)

       ! Inputs: LonFind = Desired longitude
       !         LatFind = Desired latitude
       !
       ! Outputs: iiBlock = Block index containing the desired location
       !          iiLon   = Longitude index for LonFind
       !          iiLat   = Latitude index for LatFind
       !          rLon    = Longitude interpolation scaling factor
       !          rLat    = Latitude interpolation scaling factor


       !iProc_I(iPoint) =

       !write(*,*) "Xyz_D: ",Xyz_D
       !call find_grid_block(Xyz_D, iProc_I(iPoint), iBlock)  Some day...
    end do
    if(dummyvar2 .eq. 0) then
      open(unit=12,file='testoutput2.txt',action="write",status="replace")
      do iSpecies = 1, nSpeciesTotal
        do iLons = 1, nLons
          do iLats = 1, nLats
            do iAlts = 1, nAlts
              do iBlock = 1, nBlocks
                write(12,*) Longitude(iLons,iBlock),Latitude(iLats,iBlock),&
                Altitude_GB(iLons,iLats,iAlts,iBlock),NDensityS_reshaped(iSpecies,iLons,iLats,iAlts,iBlock)
              end do!iBlock
            end do!iAlts
          end do!iLats
        end do!iLons
      end do!iSpecies
      !do i=1,nPoint
        !write(12,*) NDensityS_reshaped
      !end do
      close(unit=12)
      dummyvar2 = 1
    end if
    write(*,*) 'Maxval Data: ',maxval(Data_VI)

  !   do iPoint = 1, nPoint
  !
  !
  !   Xyz_D = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)
  !   call find_grid_block(Xyz_D, iProcFound, iBlock, iCell_D, Dist_D, &
  !        UseGhostCell = .true.)
  !        !write(*,*) "Proc and Block: ", iProcFound, iBlock
  !
  !   if(iProcFound /= iProc)then
  !      write(*,*)NameSub,' ERROR: Xyz_D, iProcFound=', Xyz_D, iProcFound
  !      call stop_mpi(NameSub//' could not find position on this proc')
  !   end if
  !   !State_V = trilinear(State_VGB(:,:,:,:,iBlock), &
  !   !     nVar, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, Xyz_D, &
  !   !     iCell_D = iCell_D, Dist_D = Dist_D)
  !
  !        !write(*,*) 'State_V: ',State_V
  !
  !   !if(UseB0)then
  !  !    call get_b0(Xyz_D, b_D)
  !  !    State_V(Bx_:Bz_) = State_V(Bx_:Bz_) + b_D
  !   !else
  !  !    b_D = 0.0
  !   !end if
  !   Data_VI(:,iPoint) = State_V*No2Si_V(iUnitCons_V)
  !
  ! end do


    ! write(*,*) 'Inside UA_get_for_gm: ',IsNew ! IsNew is T
    ! if(IsNew)then
    !    !if(DoTest)write(iUnitOut,*) NameSub,': iProc, nPoint=', iProc, nPoint
    !
    !    !if(allocated(iBlockCell_DI)) deallocate(iBlockCell_DI, Dist_DI)
    !    !allocate(iBlockCell_DI(0:nDim,nPoint), Dist_DI(nDim,nPoint))
    !
    !    do iPoint = 1, nPoint
    !
    !       Xyz_D(1:nDim) = Xyz_DI(:,iPoint)!*Si2No_V(UnitX_)
    !       write(*,*) "iProcFound: ",iProcFound,Xyz_D
    !       call find_grid_block(Xyz_D, iProcFound, iBlock, iCell_D, Dist_D, &
    !            UseGhostCell = .true.)
    !       write(*,*) "iProcFound: ",iProcFound
    !       if(iProcFound /= iProc)then
    !          write(*,*) NameSub,' ERROR: Xyz_D, iProcFound=', Xyz_D, iProcFound
    !          call stop_mpi(NameSub//' could not find position on this proc')
    !       end if
    !
    !       ! Store block and cell indexes and distances for interpolation
    !       !iBlockCell_DI(0,iPoint)      = iBlock
    !       !iBlockCell_DI(1:nDim,iPoint) = iCell_D(1:nDim)
    !       !Dist_DI(:,iPoint)            = Dist_D(1:nDim)
    !
    !    end do
    ! end if

    ! do iPoint = 1, nPoint
    !
    !    Xyz_D(1:nDim) = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)
    !
    !    ! Use stored block and cell indexes and distances
    !    iBlock          = iBlockCell_DI(0,iPoint)
    !    iCell_D(1:nDim) = iBlockCell_DI(1:nDim,iPoint)
    !    Dist_D(1:nDim)  = Dist_DI(:,iPoint)
    !
    !    ! Write STATE_VGB on STATE_V to later write it on DATA_VI!!!
    !    State_V = interpolate_vector(State_VGB(:,:,:,:,iBlock), nVar, nDim, &
    !         MinIJK_D, MaxIJK_D, iCell_D=iCell_D, Dist_D=Dist_D)
    !
    !    if(UseB0)then
    !       call get_b0(Xyz_D, B0_D)
    !       State_V(Bx_:Bz_) = State_V(Bx_:Bz_) + B0_D
    !    end if
    !
    !    do iVarBuffer = 1, nVarBuffer
    !       iVar = iVarSource_V(iVarBuffer)
    !       Data_VI(iVarBuffer,iPoint) = State_V(iVar)*No2Si_V(iUnitCons_V(iVar))
    !    end do
    !
    ! end do

  end subroutine UA_get_for_gm

  !============================================================================
    subroutine UA_get_gm_region(numcall,NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

      ! This function will be called 3 times :
      !
      ! 1) Count outer boundary ghost cells to be obtained from GM
      !
      ! 2) Return the Xyz_DGB coordinates of these cells
      !
      ! 3) Recieve Data_VI from GM and put them into State_VGB.
      !    The indexing array iPoint_I is needed to maintain the same order as
      !    the original position array Pos_DI was given in 2)

      !use BATL_lib,     ONLY: Xyz_DGB, nBlock, Unused_B, &
      !     MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
      !     CoordMin_D, CoordMax_D, CoordMin_DB, CellSize_DB
      use ModGeometry,  ONLY: far_field_bcs_blk, r_BLK
      use ModPhysics,   ONLY: No2Si_V, UnitX_, Si2No_V, iUnitCons_V
      use ModMain,      ONLY: UseB0
      use ModB0,        ONLY: B0_DGB
      use ModAdvance,   ONLY: State_VGB, Bx_, Bz_
      use ModMultiFluid,ONLY: nIonFluid
      use ModEnergy,    ONLY: calc_energy
      use CON_coupler,  ONLY: Grid_C, GM_, nVarBuffer, iVarTarget_V
      use ModGITM,      ONLY: State_VGB_GM
      !use ModInputs,    ONLY: AltMin, AltMax
      !use ModPlanet,    ONLY: RBody
      use ModSizeGitm
      use ModSphereInterface
      !use AB_SPH_module
      !use AB_module
      !use AB_COMM_module
      !use AB_ERROR_module

      character(len=*), intent(inout):: NameVar ! List of variables
      integer,          intent(inout):: nVar    ! Number of variables in Data_VI
      integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
      integer,          intent(in)::numcall
      real, intent(inout), allocatable, optional :: Pos_DI(:,:)  ! Positions

      real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
      integer, intent(in), optional:: iPoint_I(nPoint)! Order of data

      logical :: DoCountOnly
      integer :: i, j, k, iBlock, iPoint, iVarBuffer, iVar, iLons, iLats, iAlts
      real    :: Coord_D(3), Xyz_D(3)
      real    :: rMinUA, rMaxUA, temp1, temp2

      character(len=*), parameter :: NameSub='UA_get_gm_region'
      !--------------------------------------------------------------------------
      write(*,*) "-->Starting UA_get_gm_region: ",numcall,nPoint
       DoCountOnly = nPoint < 1
       write(*,*) "DoCountOnly: ",DoCountOnly
      ! Altitude_GB contains altitude (in meters, from 0 to 300k)
      ! Latitude and Longitude contain values in radians
      ! Altitude_GB(iLons,iLats,iAlts,iBlock)
      ! Longitude(iLons,iBlock)
      ! Latitude(iLats,iBlock)
       !rMinUA = RBody + AltMin ! R min
       !rMaxUA = RBody + AltMax ! R max
       iPoint = 0
       do iBlock = 1, nBlocks
         do iLons = 1, nLons
           do iLats = 1, nLats
             do iAlts = 1, nAlts
               !Coord_D(1) = Altitude_GB(iLons,iLats,iAlts,iBlock)/1.e3*sin(Longitude(iLons,iBlock))*cos(Latitude(iLats,iBlock))
               !Coord_D(2) = Altitude_GB(iLons,iLats,iAlts,iBlock)/1.e3*sin(Longitude(iLons,iBlock))*sin(Latitude(iLats,iBlock))
               !Coord_D(3) = Altitude_GB(iLons,iLats,iAlts,iBlock)/1.e3*cos(Longitude(iLons,iBlock))
               !temp1 = sin(Longitude(iLon,iBlock))*cos(Latitude(iLat,iBlock))
               !temp2 = sin(Longitude(iLon,iBlock))*sin(Latitude(iLat,iBlock))
               !Xyz_gitm(1,iLon,iLat,iAlt,iBlock) = Altitude_GB(iLon,iLat,iAlt,iBlock)/1.e3*temp1!sin(Longitude(iLon,iBlock))*cos(Latitude(iLat,iBlock))
               !Xyz_gitm(2,iLon,iLat,iAlt,iBlock) = Altitude_GB(iLon,iLat,iAlt,iBlock)/1.e3*temp2!sin(Longitude(iLon,iBlock))*sin(Latitude(iLat,iBlock))
               !Xyz_gitm(3,iLon,iLat,iAlt,iBlock) = Altitude_GB(iLon,iLat,iAlt,iBlock)/1.e3*cos(Longitude(iLon,iBlock))
               !write(*,*) 'Coord_D: ',Coord_D

               ! Exclude points below the UA-GM bottom boundary (100 km)
               if(Altitude_GB(iLons,iLats,iAlts,iBlock)/1.e3 < 100) CYCLE

               ! Exclude points above the UA-GM bottom boundary (300 km)
               !if(Altitude_GB(iLons,iLats,iAlts,iBlock)/1.e3 > 300) CYCLE

               ! Found a point to be set by UA
               iPoint = iPoint + 1
               if(DoCountOnly) CYCLE
               !write(*,*) 'Writing POS_DI for: ',iPoint,nPoint
               if(present(Data_VI))then
                 State_VGB_GM(:,iLons,iLats,iAlts,iBlock) = Data_VI(:,iPoint)

                  !allocate(State_VGB_GM(nVarGM,nLons,nLats,nAlts,nBlocks))
                  !State_VGB_GM = Data_VI;
                         ! Put Data_VI obtained from GM into State_VGB
                         !write(*,*) "nVarBuffer: ", nVarBuffer
                         !write(*,*) 'Data_VI: ',Data_VI
                         !do iVarBuffer = 1, nVarBuffer
                        !    iVar = iVarTarget_V(iVarBuffer)
                      !      State_VGB(iVar,i,j,k,iBlock) = &
                    !             Data_VI(iVarBuffer,iPoint_I(iPoint))! &
                  !               !*Si2No_V(iUnitCons_V(iVar))
                !            !write(*,*) "Data: ", State_VGB(iVar,i,j,k,iBlock)
              !           end do
                         ! Set variables not defined by coupling
                         !if(UseElectronPressure .and. &
                         !     .not. DoCoupleVar_V(ElectronPressure_)then
                         !if(UseB0) State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                        !      State_VGB(Bx_:Bz_,i,j,k,iBlock) - B0_DGB(:,i,j,k,iBlock)
                         !call calc_energy(i,i,j,j,k,k,iBlock,1,nIonFluid)
                      else
                         !write(*,*) 'Pos_DI in get_gm_region: ', iPoint,iLons,iLats,iAlts,iBlock
                         !write(*,*) 'Altitude in get_gm_region: ', Altitude_GB(iLons,iLats,iAlts,iBlock)/1.e3
                         !write(*,*) 'Pos_DI in get_gm_region: ', Xyz_gitm(:,iLons,iLats,iAlts,iBlock)
                         Pos_DI(:,iPoint) = Xyz_gitm(:,iLons,iLats,iAlts,iBlock)*1e3!Xyz_DGB(:,i,j,k,iBlock)*No2Si_V(UnitX_)
                         !write(*,*) 'Xyz_gitm: ',minval(Xyz_gitm),maxval(Xyz_gitm)
                end if
             end do
           end do
         end do
       end do

       write(*,*) "iPoint: ",iPoint
       if(DoCountOnly) nPoint = iPoint
      !
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      ! ! Find ghost cells in the UA domain
      ! iPoint = 0
      ! do iBlock = 1, nBlock
      !    do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
      !       ! Set generalized coordinate
      !       Coord_D = &
      !       CoordMin_DB(:,iBlock) + ((/i,j,k/)-0.5)*CellSize_DB(:,iBlock)
      !
      !       ! Exclude points that are inside the domain
      !       !write(*,*) "Exclude points inside domain..."
      !       !if(all(Coord_D > CoordMin_D) .and. all(Coord_D < CoordMax_D)) CYCLE
      !
      !       ! Exclude points below the UA bottom boundary
      !       if(r_BLK(i,j,k,iBlock) < rMinGm) CYCLE
      !
      !       ! Found a point to be set by UA
      !       iPoint = iPoint + 1
      !       if(DoCountOnly) CYCLE
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      !       if(present(Data_VI))then
      !          ! Put Data_VI obtained from GM into State_VGB
      !          !write(*,*) "nVarBuffer: ", nVarBuffer
      !          do iVarBuffer = 1, nVarBuffer
      !             iVar = iVarTarget_V(iVarBuffer)
      !             State_VGB(iVar,i,j,k,iBlock) = &
      !                  Data_VI(iVarBuffer,iPoint_I(iPoint))! &
      !                  !*Si2No_V(iUnitCons_V(iVar))
      !             !write(*,*) "Data: ", State_VGB(iVar,i,j,k,iBlock)
      !          end do
      !          ! Set variables not defined by coupling
      !          !if(UseElectronPressure .and. &
      !          !     .not. DoCoupleVar_V(ElectronPressure_)then
      !          !if(UseB0) State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
      !         !      State_VGB(Bx_:Bz_,i,j,k,iBlock) - B0_DGB(:,i,j,k,iBlock)
      !          call calc_energy(i,i,j,j,k,k,iBlock,1,nIonFluid)
      !       else
      !          ! Provide position to GM
                !Pos_DI(:,iPoint) = 1.1!Xyz_DGB(:,i,j,k,iBlock)*No2Si_V(UnitX_)
      !       end if
      !
      !    end do; end do; end do
      ! end do
      !
      ! if(DoCountOnly) nPoint = iPoint


    end subroutine UA_get_gm_region
  !============================================================================
  subroutine UA_put_from_gm( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    use BATL_lib, ONLY: nDim
    use ModGITM,      ONLY: State_VGB_GM

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real,    intent(in), optional  :: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional  :: iPoint_I(nPoint)! Order of data
    real, intent(out), allocatable, optional:: Pos_DI(:,:)  ! Position vectors

    ! Finish the initialization after the first coupling
    logical:: IsFirstTime = .true.

    logical:: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='UA_put_from_gm'
    !--------------------------------------------------------------------------
    write(*,*) "-->Starting UA_put_from_gm..."
    !call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(.not. present(Data_VI))then
       nPoint=0;
       ! get nPoint
       write(*,*) 'Before get_gm_region: ',NameVar
       write(*,*) 'Before get_gm_region: ',nVar,nPoint
       !write(*,*) 'Before get_gm_region: ',Pos_DI
       call UA_get_gm_region(1,NameVar, nVar, nPoint, Pos_DI)
       write(*,*) "After get_gm_region: ",NameVar
       write(*,*) 'After get_gm_region: ',nVar,nPoint,nDim
       !write(*,*) 'After get_gm_region: ',Pos_DI

       if(allocated(Pos_DI)) deallocate(Pos_DI)
       allocate(Pos_DI(nDim,nPoint))
       Pos_DI = 0.0
       ! get Pos_DI
       !write(*,*) 'Before get_gm_region 2: ',Pos_DI
       call UA_get_gm_region(2,NameVar, nVar, nPoint, Pos_DI)
       !write(*,*) 'After get_gm_region 2: ',Pos_DI
       write(*,*) 'Max,Min in Pos_DI: ',maxval(Pos_DI),minval(Pos_DI)
       !write(*,*) "Pos11: ", Pos_DI(1,1)
       !write(*,*) "Pos12: ", Pos_DI(1,2)
       !write(*,*) "Pos21: ", Pos_DI(2,1)

       RETURN
    end if

    ! set State variables
    write(*,*) "Data before: ", maxval(State_VGB_GM) ! It is getting what is sent from GM
    call UA_get_gm_region(3,NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)
    write(*,*) "Data after: ", maxval(State_VGB_GM)
    !write(*,*) "Data after: ", State_VGB_GM(:,1,1,1,1)
    write(*,*) 'nPoint at the end of UA_put_from_ua: ',nPoint

  end subroutine UA_put_from_gm
  !============================================================================
  subroutine UA_put_from_gm_dt(DtSiIn)

    real, intent(in) :: DtSiIn

    character(len=*), parameter :: NameSub = 'UA_put_from_gm_dt'
    !--------------------------------------------------------------------------
    write(*,*) "-->Starting UA_put_from_gm_dt..."
  end subroutine UA_put_from_gm_dt
  !============================================================================
  subroutine UA_put_from_gm_init(nParamInt, nParamReal, iParam_I, Param_I, &
       NameVar)

    integer, intent(in)         :: nParamInt, nParamReal! number of parameters
    integer, intent(in)         :: iParam_I(nParamInt)  ! integer parameters
    real,    intent(in)         :: Param_I(nParamReal)  ! real parameters
    character(len=*), intent(in):: NameVar              ! names of variables

    character(len=*), parameter :: NameSub = 'UA_put_from_gm_init'
    !--------------------------------------------------------------------------
    ! store GM's nDim, so it is reported as PC's nDim for the point coupler
    write(*,*) "-->Starting UA_put_from_gm_init..."
  end subroutine UA_put_from_gm_init
  !============================================================================

end module UA_wrapper
