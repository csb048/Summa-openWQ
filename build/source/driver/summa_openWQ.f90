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
      nSnow_2openwq = max( gru_struc(iGRU)%hruInfo(iHRU)%nSnow, nSnow_2openwq )
    enddo
  enddo

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
          scalarCanopyWat_summa_kg_m2   => progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyWat)%dat(1)        ,&
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
        if(scalarCanopyWat_summa_kg_m2 /= valueMissing) then
          canopyWatVol_stateVar(openWQArrayIndex) = scalarCanopyWat_summa_kg_m2 * hru_area_m2 / iden_water
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
  USE var_lookup,   only: iLookATTR  ! named variables for real valued attribute data structure
  USE var_lookup,   only: iLookINDEX 
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

  ! compartment indexes in OpenWQ (defined in the hydrolink)
  integer(i4b)                           :: canopy_index_openwq    = 0
  integer(i4b)                           :: snowSoil_index_openwq  = 1
  integer(i4b)                           :: aquifer_index_openwq   = 2
  integer(i4b)                           :: OpenWQindex_s
  integer(i4b)                           :: OpenWQindex_r
  integer(i4b)                           :: iy_r
  integer(i4b)                           :: iz_r
  integer(i4b)                           :: iy_s
  integer(i4b)                           :: iz_s
  real(rkind)                            :: wflux_s2r
  real(rkind)                            :: wmass_source

  ! Summa to OpenWQ units
  ! DomainVars
  real(rkind)                            :: hru_area_m2
  ! PrecipVars
  real(rkind)                            :: scalarRainfall_summa_m3
  real(rkind)                            :: scalarSnowfall_summa_m3
  real(rkind)                            :: scalarThroughfallRain_summa_m3
  real(rkind)                            :: scalarThroughfallSnow_summa_m3
  ! CanopyVars
  real(rkind)                            :: canopyStorWat_kg_m3
  real(rkind)                            :: scalarCanopySnowUnloading_summa_m3
  real(rkind)                            :: scalarCanopyLiqDrainage_summa_m3
  real(rkind)                            :: scalarCanopyTranspiration_summa_m3
  real(rkind)                            :: scalarCanopyEvaporation_summa_m3
  real(rkind)                            :: scalarCanopySublimation_summa_m3
  ! Snow_SoilVars
  real(rkind)                            :: mLayerLiqFluxSnow_summa_m3
  real(rkind)                            :: mLayerLiqFluxSoil_summa_m3
  real(rkind)                            :: mLayerVolFracWat_summa_m3
  real(rkind)                            :: scalarSnowSublimation_summa_m3
  real(rkind)                            :: scalarGroundEvaporation_summa_m3
  ! AquiferVars
  real(rkind)                            :: scalarAquiferBaseflow_summa_m3
  real(rkind)                            :: scalarAquiferRecharge_summa_m3
  real(rkind)                            :: scalarAquiferStorage_summa_m3
  real(rkind)                            :: scalarAquiferTranspire_summa_m3

  simtime(1) = timeStruct%var(iLookTIME%iyyy)  ! Year
  simtime(2) = timeStruct%var(iLookTIME%im)    ! month
  simtime(3) = timeStruct%var(iLookTIME%id)    ! hour
  simtime(4) = timeStruct%var(iLookTIME%ih)    ! day
  simtime(5) = timeStruct%var(iLookTIME%imin)  ! minute

  hru_index = 0

  ! Summa does not have a y-direction, 
  ! so the dimension will always be 1
  iy_r = 1 
  iy_s = 1

  do iGRU=1,nGRU
    do iHRU=1,gru_struc(iGRU)%hruCount
      print*, hru_index
      hru_index = hru_index + 1

      ! ####################################################################
      ! Associate relevant variables
      ! ####################################################################

      DomainVars: associate( &
        ! General Summa info
        hru_area_m2 => summa1_struc%attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)  &
      )

      PrecipVars: associate( &
        ! Precipitation 
        scalarRainfall_summa_kg_m2_s1             => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarRainfall)%dat(1)           ,&
        scalarSnowfall_summa_kg_m2_s1             => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSnowfall)%dat(1)           ,&
        scalarThroughfallRain_summa_kg_m2_s1      => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarThroughfallRain)%dat(1)           ,&
        scalarThroughfallSnow_summa_kg_m2_s1      => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarThroughfallSnow)%dat(1)           &
      )

      CanopyVars: associate( &     
        ! Canopy           
        scalarCanopyWat_summa_kg_m2               => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyWat)%dat(1)  ,&
        scalarCanopySnowUnloading_summa_kg_m2_s1  => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1)       ,&
        scalarCanopyLiqDrainage_summa_kg_m2_s1    => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)         ,&
        scalarCanopyTranspiration_summa_kg_m2_s1  => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyTranspiration)%dat(1)       ,& 
        scalarCanopyEvaporation_summa_kg_m2_s1    => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)         ,&
        scalarCanopySublimation_summa_kg_m2_s1    => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopySublimation)%dat(1)         &
      )

      Snow_SoilVars: associate(&
        ! Snow + Soil - Control Volume
        current_nSnow                             => summa1_struc%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)             ,&
        current_nSoil                             => summa1_struc%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)             ,&
        nSnow                                     => gru_struc(iGRU)%hruInfo(iHRU)%nSnow                                                  ,&
        nSoil                                     => gru_struc(iGRU)%hruInfo(iHRU)%nSoil                                                  ,& 
        ! Layer depth and water frac
        mLayerDepth_summa_m                       => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat(:)      ,&
        mLayerVolFracWat_summa_frac               => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracWat)%dat(:) ,&
        ! Snow Fluxes
        scalarSnowSublimation_summa_kg_m2_s1      => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSnowSublimation)%dat(1)         ,& 
        iLayerLiqFluxSnow_summa_m_s               => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%iLayerLiqFluxSnow)%dat(:)               ,&
        ! Soil Fluxes
        scalarGroundEvaporation_summa_kg_m2_s1    => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarGroundEvaporation)%dat(1)         ,&
        iLayerLiqFluxSoil_summa_m_s               => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%iLayerLiqFluxSoil)%dat(:)               &
      )

      AquiferVars: associate(&
        ! Aquifer
        scalarAquiferStorage_summa_m              => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarAquiferStorage)%dat(1), &
        scalarAquiferRecharge_summa_m_s           => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferRecharge)%dat(1)              , &        
        scalarAquiferBaseflow_summa_m_s           => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferBaseflow)%dat(1)              , &        
        scalarAquiferTranspire_summa_m_s          => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferTranspire)%dat(1)               &
      )

      ! ####################################################################
      ! Converte associate variable units: from SUMMA to OpenWQ units
      ! Here only scalar/unlayered variables
      ! OpenWQ: Volume (m3), Time (sec)
      ! Where: Vol in kg/m2, then convert to m3 by multipling by (hru_area_m2 / iden_water)
      ! Where: Flux in kg/m2/s, then convert to m3/time_step by multiplying by (hru_area_m2 * data_step / iden_water)
      ! ####################################################################

      ! PrecipVars
      scalarRainfall_summa_m3 = scalarRainfall_summa_kg_m2_s1               * hru_area_m2 * data_step / iden_water
      scalarSnowfall_summa_m3 = scalarSnowfall_summa_kg_m2_s1               * hru_area_m2 * data_step / iden_water
      scalarThroughfallRain_summa_m3 = scalarThroughfallRain_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water ! flux
      scalarThroughfallSnow_summa_m3 = scalarThroughfallSnow_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water ! flux

      ! CanopyVars
      canopyStorWat_kg_m3 = scalarCanopyWat_summa_kg_m2 * hru_area_m2 / iden_water ! vol
      scalarCanopySnowUnloading_summa_m3 = scalarCanopySnowUnloading_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water ! flux
      scalarCanopyLiqDrainage_summa_m3 = scalarCanopyLiqDrainage_summa_kg_m2_s1     * hru_area_m2 * data_step / iden_water ! flux
      scalarCanopyTranspiration_summa_m3 = scalarCanopyTranspiration_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water ! flux
      scalarCanopyEvaporation_summa_m3 = scalarCanopyEvaporation_summa_kg_m2_s1     * hru_area_m2 * data_step / iden_water ! flux
      scalarCanopySublimation_summa_m3 = scalarCanopySublimation_summa_kg_m2_s1     * hru_area_m2 * data_step / iden_water ! flux

      ! Snow_SoilVars (unlayered variables)
      ! Other variables are layered and added below as needed
      scalarSnowSublimation_summa_m3 = scalarSnowSublimation_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water
      scalarGroundEvaporation_summa_m3 = scalarGroundEvaporation_summa_kg_m2_s1 * hru_area_m2 * data_step / iden_water

      ! AquiferVars
      scalarAquiferStorage_summa_m3 = scalarAquiferStorage_summa_m        * hru_area_m2
      scalarAquiferRecharge_summa_m3 = scalarAquiferRecharge_summa_m_s    * hru_area_m2 * data_step
      scalarAquiferBaseflow_summa_m3 = scalarAquiferBaseflow_summa_m_s    * hru_area_m2 * data_step
      scalarAquiferTranspire_summa_m3 = scalarAquiferTranspire_summa_m_s  * hru_area_m2 * data_step
      
      ! ####################################################################
      ! Apply Fluxes
      ! Call RunSpaceStep
      ! ####################################################################

      ! --------------------------------------------------------------------
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! 1 Fluxes involving the canopy
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! --------------------------------------------------------------------

      if(scalarCanopyWat_summa_kg_m2 /= valueMissing) then
        
        ! ====================================================
        ! 1.1 precipitation -> canopy 
        ! ====================================================
        ! *Source*:
        ! PRECIP (external flux, so need call run_space_in) 
        ! *Recipient*: canopy (only 1 z layer)
        iz_r = 1 
        ! *Flux*: the portion of rainfall and snowfall not throughfall
        wflux_s2r = (scalarRainfall_summa_m3 - scalarThroughfallRain_summa_m3) &
                    + (scalarSnowfall_summa_m3 - scalarThroughfallSnow_summa_m3)
        ! *Call run_space_in*
        err=openwq_obj%run_space_in(                                            &
          simtime,                                                              &
          'PRECIP',                                                             &
          canopy_index_openwq, hru_index, iy_r, iz_r,                           &
          wflux_s2r)

        ! ====================================================
        ! 1.2 canopy -> upper snow/soil upper layer
        ! scalarCanopySnowUnloading + scalarCanopyLiqDrainage
        ! ====================================================
        ! *Source*
        ! canopy (only 1 z layer)
        OpenWQindex_s = canopy_index_openwq
        iz_s          = 1
        wmass_source = canopyStorWat_kg_m3
        ! *Recipient*
        ! snow+soil (upper layer: z = 1)
        OpenWQindex_r = snowSoil_index_openwq
        iz_r = 1
        ! *Flux*
        ! snow uloading + liq drainage
        wflux_s2r = scalarCanopySnowUnloading_summa_m3 &
                      + scalarCanopyLiqDrainage_summa_m3
        ! *Call run_space*
        err=openwq_obj%run_space(                                                 &
          simtime,                                                                &
          OpenWQindex_s, hru_index, iy_s, iz_s,                     &
          OpenWQindex_r, hru_index, iy_r, iz_r,                    &
          wflux_s2r,  &
          wmass_source)
        
        ! ====================================================
        ! 1.3 canopy -> OUT (lost from model) (Transp + Evap + Subl)
        ! ====================================================
        ! *Source*:
        ! canopy (only 1 z layer)
        OpenWQindex_s = canopy_index_openwq
        iz_s = 1
        wmass_source = canopyStorWat_kg_m3
        ! *Recipient*: 
        ! lost from system
        OpenWQindex_r = -1
        iz_r = -1
        ! *Flux*
        ! transpiration + evaporation + sublimation
        wflux_s2r = scalarCanopyTranspiration_summa_m3    &
                      + scalarCanopyEvaporation_summa_m3  &
                      + scalarCanopySublimation_summa_m3
        ! *Call run_space*
        err=openwq_obj%run_space(                                                   &
          simtime,                                                                  &
          OpenWQindex_s, hru_index, iy_s, iz_s,                       &
          OpenWQindex_r, hru_index, iy_r, iz_r,                                   &
          wflux_s2r,                                                              &
          wmass_source)

      endif
      
      ! --------------------------------------------------------------------
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! 2. Snow + Soil Fluxes
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! --------------------------------------------------------------------
      
      ! ====================================================
      ! 2.1 precicipitation -> upper snow/soil layer
      ! scalarThroughfallRain + scalarThroughfallSnow
      ! ====================================================

      ! *Source*:
      ! PRECIP (external flux, so need call run_space_in)
      ! *Recipient*: 
      ! snow+soil (upper layer)
      OpenWQindex_r = snowSoil_index_openwq
      iz_r = 1
      ! *Flux*
      ! throughfall rain and snow
      wflux_s2r = scalarThroughfallRain_summa_m3 &
                    + scalarThroughfallSnow_summa_m3
      ! *Call run_space_in*
      err=openwq_obj%run_space_in(                                            &
        simtime,                                                              &
        'PRECIP',                                                             &
        OpenWQindex_r, hru_index, iy_r, iz_r,                  &
        wflux_s2r       &
        )
      
      ! Snow fluxes
      if (current_nSnow /= 0) then

        ! ====================================================
        ! 2.2 snow -> OUT (lost from model) (sublimation)
        ! only top snow layer! ->>>> NEED TO CONFIRM THIS
        ! ====================================================
        ! *Source*:
        ! snow (upper layer)
        OpenWQindex_s = snowSoil_index_openwq
        iz_s = 1
        mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(1) * hru_area_m2 * mLayerDepth_summa_m(1)
        wmass_source = mLayerVolFracWat_summa_m3
        ! *Recipient*: 
        ! lost from system
        OpenWQindex_r = -1
        iz_r  = -1
        ! *Flux*
        ! snow sublimation
        wflux_s2r = scalarSnowSublimation_summa_m3
        ! *Call run_space*
        err=openwq_obj%run_space(                                     &
          simtime,                                                    &
          OpenWQindex_s, hru_index, iy_s, iz_s,        &
          -1, -1, -1, -1,                                             & ! lost
          wflux_s2r,                                      &
          wmass_source)

        ! ====================================================
        ! 2.3 snow internal fluxes
        ! ====================================================
        do iLayer = 1, nSnow-1 ! last layer of snow becomes different fluxes 

          ! *Source*: 
          OpenWQindex_s = snowSoil_index_openwq
          iz_s = iLayer
          mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(iLayer) * hru_area_m2 * mLayerDepth_summa_m(iLayer)
          wmass_source = mLayerVolFracWat_summa_m3
          ! *Recipient*: 
          OpenWQindex_r = snowSoil_index_openwq
          iz_r = iLayer + 1
          ! *Flux*
          mLayerLiqFluxSnow_summa_m3 = iLayerLiqFluxSnow_summa_m_s(iLayer) * hru_area_m2 * data_step
          wflux_s2r = mLayerLiqFluxSnow_summa_m3 
          ! *Call run_space*
          err=openwq_obj%run_space(                               &
            simtime,                                              &
            OpenWQindex_s, hru_index, iy_s, iz_s,  &
            OpenWQindex_r, hru_index, iy_r, iz_r,  &
            wflux_s2r,                 & 
            wmass_source)

        end do

      end if
      
      ! ====================================================
      ! 2.4 soil fluxes
      ! upper soil -> OUT (lost from system) (evaporatoon)
      ! ====================================================

      ! *Source*: 
      OpenWQindex_s = snowSoil_index_openwq
      iz_s = 1
      mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(nSnow+1) * hru_area_m2 * mLayerDepth_summa_m(nSnow+1)
      wmass_source = mLayerVolFracWat_summa_m3
      ! *Recipient*: 
      ! lost from system
      OpenWQindex_r = -1
      iz_r = -1
      ! *Flux*
      wflux_s2r = scalarGroundEvaporation_summa_m3
      ! *Call run_space*
      err=openwq_obj%run_space(                                     &
        simtime,                                                    &
        OpenWQindex_s, hru_index, iy_s, iz_s,        &
        OpenWQindex_r, hru_index, iy_r, iz_r,  &
        wflux_s2r,                           &
        wmass_source)

      ! ====================================================
      ! 2.5 Infiltrations
      ! which is: scalarRainPlusMelt = iLayerLiqFluxSnow(nSnow)
      ! ====================================================
      ! *Source*:
      ! snow lower layer
      OpenWQindex_s = snowSoil_index_openwq
      iz_s = nSnow; 
      mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(nSnow) * hru_area_m2 * mLayerDepth_summa_m(nSnow)
      wmass_source = mLayerVolFracWat_summa_m3
      ! *Recipient*:
      ! soil top layer
      iz_r = nSnow+1;
      OpenWQindex_r = snowSoil_index_openwq
      iz_r = 1           
      ! *Flux*
      mLayerLiqFluxSnow_summa_m3 = iLayerLiqFluxSnow_summa_m_s(nSnow) * hru_area_m2 * data_step
      wflux_s2r = mLayerLiqFluxSnow_summa_m3
      ! *Call run_space*
      err=openwq_obj%run_space(                               &
        simtime,                                              &
        OpenWQindex_s, hru_index, iy_s, iz_s,  &
        OpenWQindex_r, hru_index, iy_r, iz_r,  &
        wflux_s2r,                 & 
        wmass_source)

      ! ====================================================
      ! 2.6 soil internal fluxes
      ! ====================================================
      do iLayer = 1, nSoil - 1 ! last layer of soil becomes different fluxes

        ! *Source*:
        ! soil layer iLayer
        OpenWQindex_s = snowSoil_index_openwq
        iz_s = nSnow + iLayer
        mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(iLayer+nSnow) * hru_area_m2 * mLayerDepth_summa_m(iLayer+nSnow)
        wmass_source = mLayerVolFracWat_summa_m3
        ! *Recipient*: 
        ! soi layer iLayer+1
        OpenWQindex_r = snowSoil_index_openwq
        iz_r = nSnow + iLayer + 1
        ! *Flux*
        ! flux between soil layer (it's -1 because the first layer gets)
        mLayerLiqFluxSoil_summa_m3 = iLayerLiqFluxSoil_summa_m_s(iLayer-1) * hru_area_m2 * data_step
        wflux_s2r = mLayerLiqFluxSoil_summa_m3 
        ! *Call run_space*
        err=openwq_obj%run_space(                               &
          simtime,                                              &
          OpenWQindex_s, hru_index, iy_s, iz_s,  &
          OpenWQindex_r, hru_index, iy_r, iz_r,  & 
          wflux_s2r,                 & 
          wmass_source)
      end do

      ! --------------------------------------------------------------------
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! 3 Aquifer Fluxes
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! --------------------------------------------------------------------
      
      ! ====================================================
      ! 3.1 Lower soil layer -> Aquifer
      ! ====================================================

      ! *Source*: 
      ! soil lower layer
      OpenWQindex_s = snowSoil_index_openwq
      iz_s = nSoil + nSoil 
      mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(nSoil + nSoil) * hru_area_m2 * mLayerDepth_summa_m(nSoil + nSoil)
      wmass_source = mLayerVolFracWat_summa_m3
      ! *Recipient*: 
      ! aquifer (only 1 layer)
      OpenWQindex_r = aquifer_index_openwq
      iz_r = 1
      ! *Flux*
      wflux_s2r = scalarAquiferRecharge_summa_m3
      ! *Call run_space*
      err=openwq_obj%run_space(                           &
        simtime,                                          &
        OpenWQindex_s, hru_index, iy_s, iz_s,    &
        OpenWQindex_r, hru_index, iy_r, iz_r,       &
        wflux_s2r,                 & 
        wmass_source)

      ! ====================================================
      ! 3.2 Aquifer -> OUT (lost from model) (baseflow) 
      ! ->>>> NEED TO CONFIRM THAT THIS IS REALLY LOST
      ! ====================================================

      ! *Source*: 
      ! aquifer (only 1 z layer)
      OpenWQindex_s = aquifer_index_openwq
      iz_s = 1
      wmass_source = scalarAquiferStorage_summa_m3
      ! *Recipient*: 
      ! lost from system
      OpenWQindex_r = -1
      iz_r = -1
      ! *Flux*
      wflux_s2r = scalarAquiferBaseflow_summa_m3
      ! *Call run_space*
      err=openwq_obj%run_space(                           &
        simtime,                                          &
        OpenWQindex_s, hru_index, iy_s, iz_s,    &
        OpenWQindex_r, hru_index, iy_r, iz_r,       &
        wflux_s2r,                 & 
        wmass_source)

      ! ====================================================
      ! 3.3 Aquifer -> OUT (lost from model) (transpiration) 
      ! ====================================================
      
      ! *Source*: 
      ! aquifer (only 1 z layer)
      OpenWQindex_s = aquifer_index_openwq
      iz_s = 1
      wmass_source = scalarAquiferStorage_summa_m3
      ! *Recipient*: 
      ! lost from system
      OpenWQindex_r = -1
      iz_r = -1
      ! *Flux*
      wflux_s2r = scalarAquiferTranspire_summa_m3
      ! *Call run_space*
      err=openwq_obj%run_space(                         &
        simtime,                                        &
        OpenWQindex_s, hru_index, iy_s, iz_s,    &
        OpenWQindex_r, hru_index, iy_r, iz_r,       &
        wflux_s2r,                 & 
        wmass_source)

      
      end associate AquiferVars
      end associate Snow_SoilVars
      end associate CanopyVars
      end associate PrecipVars
      end associate DomainVars
      
      
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