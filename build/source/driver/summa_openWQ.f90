module summa_openWQ
  USE nrtype
  USE openWQ, only:ClassWQ_OpenWQ
  USE data_types, only:gru_hru_doubleVec
  implicit none
  private
  ! Subroutines
  public :: init_openwq
  public :: run_time_start
  public :: run_time_start_go
  public :: run_space_step
  public :: run_time_end

  ! Global Data for prognostic Variables of HRUs
  type(gru_hru_doubleVec),save,public   :: progStruct_timestep_start ! copy of progStruct at the start of timestep for passing fluxes

  contains

  ! Subroutine to initalize the openWQ object
  ! putting it here to keep the SUMMA_Driver clean
subroutine init_openwq(err, message)

  USE globalData,only:openWQ_obj
  USE globalData,only:gru_struc                               ! gru-hru mapping structures
  USE globalData,only:prog_meta
  USE allocspace_module,only:allocGlobal                      ! module to allocate space for global data structures

  implicit none

  integer(i4b),intent(inout)                      :: err
  character(*),intent(inout)                      :: message         ! error messgage
  integer(i4b)                                    :: hruCount
  integer(i4b)                                    :: nCanopy_2openwq
  integer(i4b)                                    :: nSnow_2openwq
  integer(i4b)                                    :: nSoil_2openwq
  integer(i4b)                                    :: nAquifer_2openwq
  integer(i4b)                                    :: nYdirec_2openwq     ! number of layers in the y-dir (not used in summa)
  integer(i4b)                                    :: iGRU, iHRU          ! indices of GRUs and HRUs

  openwq_obj = ClassWQ_OpenWQ() ! initalize openWQ object

  ! nx -> num of HRUs)
  hruCount = sum( gru_struc(:)%hruCount )

  ! ny -> this seems to be fixes because SUMMA is based on the HRU concept, so grids are serialized)
  nYdirec_2openwq = 1

  ! Openwq nz (number of layers)
  nCanopy_2openwq = 1       ! Cannopy has only 1 layer
  nAquifer_2openwq = 1      ! GW has only 1 layer
  nSoil_2openwq = 0   ! Soil may have multiple layers, and gru-hrus may have different values
  nSnow_2openwq = 0   ! Snow has multiple layers, and gru-hrus may have different values (up to 5 layers)
  do iGRU = 1, size(gru_struc(:))
    do iHRU = 1, gru_struc(iGRU)%hruCount
      nSoil_2openwq = max( gru_struc(iGRU)%hruInfo(iHRU)%nSoil, nSoil_2openwq )
      !num_layers_volFracWat = max( gru_struc(iGRU)%hruInfo(iHRU)%nSoil, num_layers_volFracWat )
    enddo
  enddo
  nSnow_2openwq = 5 ! maximum number of snow layers

  ! intialize openWQ
  err=openwq_obj%decl(    &
    hruCount,             & ! num HRU
    nCanopy_2openwq,      & ! num layers of canopy (fixed to 1)
    nSnow_2openwq,        & ! num layers of snow (fixed to max of 5 because it varies)
    nSoil_2openwq,        & ! num layers of snoil (variable)
    nAquifer_2openwq,     & ! num layers of aquifer (fixed to 1)
    nYdirec_2openwq)             ! num of layers in y-dir (set to 1 because not used in summa)
  
  ! Create copy of state information, needed for passing to openWQ with fluxes that require
  ! the previous time_steps volume
  call allocGlobal(prog_meta, progStruct_timestep_start, err, message) 

end subroutine init_openwq
  
! Subroutine that SUMMA calls to pass varialbes that need to go to
! openWQ - the copy of progStruct is done in here
subroutine run_time_start(  &
    openWQ_obj,             & ! passing openwq object
    summa1_struc)

  USE summa_type, only: summa1_type_dec            ! master summa data type
  USE globalData, only: gru_struc

  implicit none

  ! Dummy Varialbes
  class(ClassWQ_OpenWQ), intent(in)  :: openWQ_obj
  type(summa1_type_dec), intent(in)  :: summa1_struc
  ! local variables
  integer(i4b)                       :: iGRU
  integer(i4b)                       :: iHRU
  integer(i4b)                       :: nSoil_2openwq ! maximum number of layers for soil
  integer(i4b)                       :: nSnow_2openwq ! maximum number of layers for snow)
  integer(i4b)                       :: err

  ! Get number of soil and snow layers
  ! Needs to be isolated because explicit-shaped arrays can only be defined with parameters 
  ! or int passed as an argument
  nSoil_2openwq = 0
  nSnow_2openwq = 0
  do iGRU = 1, size(gru_struc(:))
    do iHRU = 1, gru_struc(iGRU)%hruCount
      ! nSnow_2openwq = max( gru_struc(iGRU)%hruInfo(iHRU)%nSnow, nSnow_2openwq )
      nSoil_2openwq = max( gru_struc(iGRU)%hruInfo(iHRU)%nSoil, nSoil_2openwq )
    enddo
  enddo

  call run_time_start_go( &
    openwq_obj,           &
    summa1_struc,         &
    nSnow_2openwq,        &
    nSoil_2openwq)

end subroutine

subroutine run_time_start_go( &
    openWQ_obj,               &
    summa1_struc,             &
    nSnow_2openwq,            &
    nSoil_2openwq)

  USE summa_type, only: summa1_type_dec            ! master summa data type
  USE globalData, only: gru_struc
  USE var_lookup, only: iLookPROG  ! named variables for state variables
  USE var_lookup, only: iLookTIME  ! named variables for time data structure
  USE var_lookup, only: iLookATTR  ! named variables for real valued attribute data structure
  USE multiconst,only:&
                        iden_ice,       & ! intrinsic density of ice             (kg m-3)
                        iden_water        ! intrinsic density of liquid water    (kg m-3)

  implicit none

  ! Dummy Varialbes
  class(ClassWQ_OpenWQ), intent(in)   :: openWQ_obj
  type(summa1_type_dec), intent(in)   :: summa1_struc
  ! local variables
  integer(i4b), intent(in)            :: nSnow_2openwq
  integer(i4b), intent(in)            :: nSoil_2openwq
  integer(i4b)                        :: iGRU
  integer(i4b)                        :: iHRU
  integer(i4b)                        :: ilay
  integer(i4b)                        :: iVar
  integer(i4b)                        :: iDat
  integer(i4b)                        :: openWQArrayIndex
  integer(i4b)                        :: simtime(5) ! 5 time values yy-mm-dd-hh-min
  real(rkind)                         :: airTemp_K_depVar(sum(gru_struc(:)%hruCount))
  real(rkind)                         :: canopyWatVol_stateVar(sum(gru_struc(:)%hruCount))
  real(rkind)                         :: aquiferWatVol_stateVar(sum(gru_struc(:)%hruCount))
  real(rkind)                         :: sweWatVol_stateVar(sum(gru_struc(:)%hruCount), nSnow_2openwq)
  real(rkind)                         :: soilWatVol_stateVar(sum(gru_struc(:)%hruCount), nSoil_2openwq)
  real(rkind)                         :: soilTemp_K_depVar(sum(gru_struc(:)%hruCount), nSoil_2openwq)
  real(rkind)                         :: soilMoist_depVar(sum(gru_struc(:)%hruCount), nSoil_2openwq)
  integer(i4b)                        :: err
  real(rkind),parameter              :: valueMissing=-9999._rkind   ! seems to be SUMMA's default value for missing data

  summaVars: associate(&
      progStruct     => summa1_struc%progStruct             , &
      timeStruct     => summa1_struc%timeStruct             , &
      attrStruct     => summa1_struc%attrStruct               &
  )

  ! Update dependencies and storage volumes
  ! Assemble the data to send to openWQ

  openWQArrayIndex = 0 ! index into the arrays that are being passed to openWQ

  do iGRU = 1, size(gru_struc(:))
      do iHRU = 1, gru_struc(iGRU)%hruCount

        GeneralVars: associate(&
          hru_area_m2                 => attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)                      ,&
          Tair_summa_K                => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanairTemp)%dat(1)      ,&
          CanopyStorWat_summa_kg_m2   => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyWat)%dat(1)        ,&
          AquiferStorWat_summa_m      => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarAquiferStorage)%dat(1)   &
        )

        openWQArrayIndex = openWQArrayIndex + 1 

        ! ############################
        ! Update unlayered variables and dependencies 
        ! (1 layer only)
        ! ############################

        ! Tair 
        ! (Summa in K)
        if(Tair_summa_K /= valueMissing) then
          airTemp_K_depVar(openWQArrayIndex) =  Tair_summa_K 
        endif
          
        ! Vegetation
        ! unit for volume = m3 (summa-to-openwq unit conversions needed)
        ! scalarCanopyWat [kg m-2], so needs to  to multiply by hru area [m2] and divide by water density
        if(CanopyStorWat_summa_kg_m2 /= valueMissing) then
          canopyWatVol_stateVar(openWQArrayIndex) = CanopyStorWat_summa_kg_m2 * hru_area_m2 / iden_water
        endif

        ! Aquifer
        ! unit for volume = m3 (summa-to-openwq unit conversions needed)
        ! scalarAquiferStorage [m], so needs to  to multiply by hru area [m2] only
        if(AquiferStorWat_summa_m /= valueMissing) then
          aquiferWatVol_stateVar(openWQArrayIndex) = AquiferStorWat_summa_m * hru_area_m2
        endif 
        
        ! ############################
        ! Update layered variables and dependenecies
        ! ############################

        ! Soil
        do ilay = 1, nSoil_2openwq
          
          SoilVars: associate(&
            Tsoil_summa_K => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerTemp)%dat(ilay)        ,&
            Wsoil_summa_m => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerMatricHead)%dat(ilay)   &
          )
          ! Tsoil
          ! (Summa in K)
          if(Tsoil_summa_K /= valueMissing) then
            soilTemp_K_depVar(openWQArrayIndex, ilay) = Tsoil_summa_K
          endif

          soilMoist_depVar(openWQArrayIndex, ilay) = 0     ! TODO: Find the value for this varaibles

          ! Soil
          ! unit for volume = m3 (summa-to-openwq unit conversions needed)
          ! mLayerMatricHead [m], so needs to  to multiply by hru area [m2]
          if(Wsoil_summa_m /= valueMissing) then
            soilWatVol_stateVar(openWQArrayIndex, ilay) = Wsoil_summa_m * hru_area_m2
          endif

          end associate SoilVars

        enddo

        ! Snow
        do ilay = 1, nSnow_2openwq

          SnowVars: associate(&
            mLayerDepth      => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat(ilay)         , &    ! depth of each layer (m)
            mLayerVolFracIce => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracIce)%dat(ilay)    , &    ! volumetric fraction of ice in each layer  (-)
            mLayerVolFracLiq => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracLiq)%dat(ilay)      &    ! volumetric fraction of liquid water in each layer (-)
          )
          
          ! Snow
          ! unit for volume = m3 (summa-to-openwq unit conversions needed)
          ! mLayerVolFracIce and mLayerVolFracLiq [-], so needs to  to multiply by hru area [m2] and divide by water density
          ! But needs to account for both ice and liquid, and convert to liquid volumes
          if(mLayerVolFracIce /= valueMissing .or. &
             mLayerVolFracLiq /= valueMissing) then

            sweWatVol_stateVar(openWQArrayIndex, ilay) =                                              &
              (max(mLayerVolFracIce, 0._rkind) * iden_ice + &
              max(mLayerVolFracLiq, 0._rkind) * iden_water) / iden_water  &
              * mLayerDepth * hru_area_m2
            else
              sweWatVol_stateVar(openWQArrayIndex, ilay) = 0._rkind

          endif

          end associate SnowVars

        enddo

        ! **************  
        ! Fluxes
        !************** 

        ! TODO

        ! Copy the prog structure
        do iVar = 1, size(progStruct%gru(iGRU)%hru(iHRU)%var)
          do iDat = 1, size(progStruct%gru(iGRU)%hru(iHRU)%var(iVar)%dat)
            progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iVar)%dat(iDat) = progStruct%gru(iGRU)%hru(iHRU)%var(iVar)%dat(iDat)
          end do
        end do

        end associate GeneralVars

      end do
  end do

  ! add the time values to the array
  simtime(1) = timeStruct%var(iLookTIME%iyyy)  ! Year
  simtime(2) = timeStruct%var(iLookTIME%im)    ! month
  simtime(3) = timeStruct%var(iLookTIME%id)    ! hour
  simtime(4) = timeStruct%var(iLookTIME%ih)    ! day
  simtime(5) = timeStruct%var(iLookTIME%imin)  ! minute
  
  err=openWQ_obj%run_time_start(&
        sum(gru_struc(:)%hruCount),             & ! total HRUs
        nSnow_2openwq,                          &
        nSoil_2openwq,                          &
        simtime,                                &
        soilMoist_depVar,                       &                    
        soilTemp_K_depVar,                      &
        airTemp_K_depVar,                       &
        sweWatVol_stateVar,                     &
        canopyWatVol_stateVar,                  &
        soilWatVol_stateVar,                    &
        aquiferWatVol_stateVar)

  ! copy progStruct values to progStruct_timestep_start

  end associate summaVars

end subroutine


subroutine run_space_step(  &
    timeStruct,             &
    summa1_struc,           &
    fluxStruct,             &
    nGRU)

  USE var_lookup,   only: iLookPROG  ! named variables for state variables
  USE var_lookup,   only: iLookTIME  ! named variables for time data structure
  USE var_lookup,   only: iLookFLUX  ! named varaibles for flux data
  USE var_lookup, only: iLookATTR  ! named variables for real valued attribute data structure
  USE summa_type,   only: summa1_type_dec            ! master summa data type
  USE globalData,   only: openWQ_obj
  USE data_types,   only: var_dlength,var_i
  USE globalData,   only: gru_struc
  USE globalData,   only: data_step   ! time step of forcing data (s)
  USE multiconst,   only:&
                        iden_ice,       & ! intrinsic density of ice             (kg m-3)
                        iden_water        ! intrinsic density of liquid water    (kg m-3)

  implicit none

  type(var_i),             intent(in)    :: timeStruct 
  type(gru_hru_doubleVec), intent(in)    :: fluxStruct
  type(summa1_type_dec),   intent(in)    :: summa1_struc
  
  integer(i4b),            intent(in)    :: nGRU

  integer(i4b)                           :: hru_index ! needed because openWQ saves hrus as a single array
  integer(i4b)                           :: iHRU      ! variable needed for looping
  integer(i4b)                           :: iGRU      ! variable needed for looping
  integer(i4b)                           :: iLayer    ! varaible needed for looping

  integer(i4b)                           :: simtime(5) ! 5 time values yy-mm-dd-hh-min
  integer(i4b)                           :: err
  real(rkind),parameter                  :: valueMissing=-9999._rkind   ! seems to be SUMMA's default value for missing data

  ! compartment indexes
  integer(i4b)                           :: scalarCanopyWat_index=0 ! SUMMA Side units: kg m-2
  integer(i4b)                           :: mLayerMatricHead_index=1 ! SUMMA Side units: m
  integer(i4b)                           :: scalarAquifer_index=2 ! SUMMA Side units: m
  integer(i4b)                           :: mLayerVolFracWat_index=3 ! SUMMA Side units: ????
  integer(i4b)                           :: iy_r
  integer(i4b)                           :: iz_r
  integer(i4b)                           :: iy_s
  integer(i4b)                           :: iz_s

  real(rkind)                            :: hru_area_m2
  real(rkind)                            :: canopyStorWat_kg_m3
  real(rkind)                            :: scalarThroughfallRain_summa_m3
  real(rkind)                            :: scalarThroughfallSnow_summa_m3
  real(rkind)                            :: fluxOUT_scalarCanopySnowUnloading_summa_m3
  real(rkind)                            :: fluxOUT_scalarCanopyLiqDrainage_summa_m3
  real(rkind)                            :: fluxOUT_scalarCanopyTranspiration_summa_m3
  real(rkind)                            :: fluxOUT_scalarCanopyEvaporation_summa_m3
  real(rkind)                            :: fluxOUT_scalarCanopySublimation_summa_m3      

  simtime(1) = timeStruct%var(iLookTIME%iyyy)  ! Year
  simtime(2) = timeStruct%var(iLookTIME%im)    ! month
  simtime(3) = timeStruct%var(iLookTIME%id)    ! hour
  simtime(4) = timeStruct%var(iLookTIME%ih)    ! day
  simtime(5) = timeStruct%var(iLookTIME%imin)  ! minute

  hru_index = 0

  do iGRU=1,nGRU
    do iHRU=1,gru_struc(iGRU)%hruCount

      hru_index = hru_index + 1

      ! ####################################################################
      ! CANOPY Fluxes
      ! ####################################################################

      ! Associate relevant variables
      CanopyVars: associate( &     
        hru_area_m2                               => summa1_struc%attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)                   ,&
        ! Canopy           
        canopyStorWat_summa_kg_m2                 => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyWat)%dat(1)  ,&
        scalarThroughfallRain_summa_kg_m2_s1      => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarThroughfallRain)%dat(1)           ,&
        scalarThroughfallSnow_summa_kg_m2_s1      => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarThroughfallSnow)%dat(1)           ,&
        scalarCanopySnowUnloading_summa_kg_m2_s1  => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1)       ,&
        scalarCanopyLiqDrainage_summa_kg_m2_s1    => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)         ,&
        scalarCanopyTranspiration_summa_kg_m2_s1  => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyTranspiration)%dat(1)       ,& 
        scalarCanopyEvaporation_summa_kg_m2_s1    => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)         ,&
        scalarCanopySublimation_summa_kg_m2_s1    => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopySublimation)%dat(1)         ,& 
        ! Snow + Soil - Control Volume
        mLayerVolFracWat_summa_kg_m2              => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracWat)%dat(:) ,&
        ! Snow Fluxes
        nSnow                                     => gru_struc(iGRU)%hruInfo(iHRU)%nSnow                                                  ,&
        mLayerLiqFluxSnow_s1                      => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%mLayerLiqFluxSnow)%dat(:)               ,&
        ! Soil Fluxes
        nSoil                                     => gru_struc(iGRU)%hruInfo(iHRU)%nSoil                                                  ,&


        ! Aquifer
        scalarAquiferStorage_summa_kg_m2          => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarAquiferStorage)%dat(1) &
      )

        ! If no canopy, then skip
        !if(canopyStorWat_summa_kg_m2 =/ valueMissing) then

        ! Convert Canopy water volume to m3
        canopyStorWat_kg_m3 = CanopyStorWat_summa_kg_m2 * hru_area_m2 / iden_water

        ! Unit conversions for the fluxes (summa data_step is in seconds)
        ! Dividing both by iden_water because we want in vol in m3 of liquid water
        scalarThroughfallRain_summa_m3 = scalarThroughfallRain_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water
        scalarThroughfallSnow_summa_m3 = scalarThroughfallSnow_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water
        fluxOUT_scalarCanopySnowUnloading_summa_m3 = scalarCanopySnowUnloading_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water
        fluxOUT_scalarCanopyLiqDrainage_summa_m3 = scalarCanopyLiqDrainage_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water
        fluxOUT_scalarCanopyTranspiration_summa_m3 = scalarCanopyTranspiration_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water
        fluxOUT_scalarCanopyEvaporation_summa_m3 = scalarCanopyEvaporation_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water
        fluxOUT_scalarCanopySublimation_summa_m3 = scalarCanopySublimation_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water

        ! .......................................................
        ! Call run_space ........................................
        ! .......................................................
        
        ! 1) Input flux (Precipitation)
        ! scalarThroughfallRain + scalarThroughfallSnow
        iy_r = 1; iz_r = 1
        err=openwq_obj%run_space_in(                      &
          simtime,                                        &
          'PRECIP',                                       &
          mLayerVolFracWat_index, hru_index, iy_r, iz_r,  &
          scalarThroughfallRain_summa_m3                  &
            + scalarThroughfallSnow_summa_m3)

        ! 2) fluxOUT to Upper soil/snow layer 
        ! scalarCanopySnowUnloading + flux_scalarCanopyLiqDrainage_m3
        ! They can be together because the source and sink compartments and ix,iy,iz are the same
        iy_s = 1; iz_s = 1; 
        iy_r = 1; iz_r = 1
        err=openwq_obj%run_space(                         &
          simtime,                                        &
          scalarCanopyWat_index, hru_index, iy_s, iz_s,   &
          mLayerVolFracWat_index, hru_index, iy_r, iz_r,  &
          fluxOUT_scalarCanopySnowUnloading_summa_m3      & 
            + fluxOUT_scalarCanopyLiqDrainage_summa_m3,   &
          canopyStorWat_kg_m3)

        ! 3) fluxOUT (lost from system) -> leads to incresed concentration
        ! fluxOUT_scalarCanopyTranspiration_summa_m3 + fluxOUT_scalarCanopyEvaporation_summa_m3 + fluxOUT_scalarCanopySublimation_summa_m3
        iy_s = 1; iz_s = 1; 
        iy_r = -1; iz_r = -1 ! -1 is the flag for no recipient inside the system (so lost from model)
        err=openwq_obj%run_space(                         &
          simtime,                                        &
          scalarCanopyWat_index, hru_index, iy_s, iz_s,   &
          mLayerVolFracWat_index, hru_index, iy_r, iz_r,  &
          fluxOUT_scalarCanopyTranspiration_summa_m3      &  
            + fluxOUT_scalarCanopyEvaporation_summa_m3    &
            + fluxOUT_scalarCanopySublimation_summa_m3,   &
          canopyStorWat_kg_m3)

      !endif
        

      ! ####################################################################
      ! Snow Fluxes
      ! ####################################################################
      do iLayer = 1, nSnow-1 ! last layer of snow becomes different fluxes
        iz_s = iLayer
        iz_r = iLayer + 1
        err=openwq_obj%run_space(                         &
          simtime,                                        &
          mLayerVolFracWat_index, hru_index, iy_s, iz_s,  &
          mLayerVolFracWat_index, hru_index, iy_r, iz_r,  &
          mLayerLiqFluxSnow_s1(iLayer),                   &
          mLayerVolFracWat_summa_kg_m2(iLayer))
      end do
      

      ! ####################################################################
      ! Soil Fluxes
      ! ####################################################################

      ! Kyle...

      ! ####################################################################
      ! Aquifer Fluxes
      ! ####################################################################

      ! Kyle...

      end associate CanopyVars

    end do
  end do

end subroutine run_space_step


subroutine run_time_end( &
  openWQ_obj, &
  summa1_struc)

  USE summa_type, only:summa1_type_dec            ! master summa data type
  
  USE var_lookup, only: iLookTIME  ! named variables for time data structure

  implicit none

  ! Dummy Varialbes
  class(ClassWQ_OpenWQ), intent(in)  :: openWQ_obj
  type(summa1_type_dec), intent(in)  :: summa1_struc

  ! Local Variables
  integer(i4b)                       :: simtime(5) ! 5 time values yy-mm-dd-hh-min
  integer(i4b)                       :: err ! error control

  summaVars: associate(&
      timeStruct     => summa1_struc%timeStruct       &       
  )

  simtime(1) = timeStruct%var(iLookTIME%iyyy)  ! Year
  simtime(2) = timeStruct%var(iLookTIME%im)    ! month
  simtime(3) = timeStruct%var(iLookTIME%id)    ! hour
  simtime(4) = timeStruct%var(iLookTIME%ih)    ! day
  simtime(5) = timeStruct%var(iLookTIME%imin)  ! minute

  err=openwq_obj%run_time_end(simtime)           ! minute

  end associate summaVars
end subroutine



end module summa_openWQ